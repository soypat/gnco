package nozzle

import (
	"errors"
	"math"

	"github.com/soypat/geometry/md1"
)

const (
	// Universal gas constant [m3⋅Pa⋅K−1⋅mol−1] also [kg⋅m^2⋅s^−2⋅K^−1⋅mol^−1]
	_Ru = 8.31446261815324
	// Standard gravity [m/s^2]
	_g0 = 9.80665
)

// Laval defines a Laval (convergent-divergent) nozzle.
// All methods are for isoentropic processes unless specified otherwise in method name.
type Laval struct {
	areaThroat float64 // [m^2]
	areaExit   float64 // [m^2]
}

// NewLaval creates a Laval nozzle with throat area at and exit area ae.
func NewLaval(at, ae float64) Laval {
	if !(at < ae) {
		panic("throat area must be less than exit area")
	}
	return Laval{areaThroat: at, areaExit: ae}
}

// Dims returns the throat and exit areas [m^2].
func (laval Laval) Dims() (at, ae float64) { return laval.areaThroat, laval.areaExit }

// EstimateExpansionRatio estimates the expansion ratio (Ae/At) from exhaust
// specific heat ratio (gamma=cp/cv) and chamber/exit pressure.
func (laval Laval) EstimateExpansionRatio(pc, pe, gamma float64) float64 {
	AtdAe := math.Pow((gamma+1)/2, 1/(gamma-1)) *
		math.Pow(pe/pc, 1/gamma) *
		math.Sqrt((gamma+1)/(gamma-1)*(1-math.Pow(pe/pc, (gamma-1)/gamma)))
	return 1 / AtdAe
}

// Cf is the thrust coefficient from chamber pressure pc,
// exit pressure pe, atmospheric pressure pa, and exhaust gamma.
// Returns 0 when pc is not greater than pe and pa.
//
//	Cf = sqrt( f(γ, pc, pe) ) + ε*(pe-pa)/pc
func (laval Laval) Cf(pc, pe, pa, gamma float64) float64 {
	if pc <= pe || pc <= pa {
		return 0
	}
	ga := gamma
	res := 2 * ga * ga / (ga - 1) * math.Pow(2/(ga+1), (ga+1)/(ga-1)) *
		(1 - math.Pow(pe/pc, (ga-1)/ga))
	er := laval.EstimateExpansionRatio(pc, pe, gamma)
	return math.Sqrt(res) + (pe-pa)*er/pc
}

// Thrust calculates ideal thrust [N] from chamber pressure pc, exit pressure pe,
// and atmospheric pressure pa.
//
//	T = At * pc * Cf
func (laval Laval) Thrust(pc, pe, pa, gamma float64) float64 {
	return laval.areaThroat * pc * laval.Cf(pc, pe, pa, gamma)
}

// Cstar is the characteristic exhaust velocity [m/s].
//
// Reference: Rocket Propulsion Elements, 8th Edition, Equation 3-32.
func (laval Laval) Cstar(ex ExhaustIsoentropic) float64 {
	ga := ex.Gamma
	return ex.SpeedOfSound(ex.MaxBurnTemperature) / ga /
		math.Pow(2/(ga+1), (ga+1)/(2*ga-2))
}

// SpecificImpulse calculates the ideal specific impulse [s].
func (laval Laval) SpecificImpulse(pc, pe, pa float64, ex ExhaustIsoentropic) float64 {
	return laval.Cf(pc, pe, pa, ex.Gamma) * laval.Cstar(ex) / _g0
}

// MassFlow in [kg/s] from chamber pressure pc and stagnation conditions.
//
// Reference: Rocket Propulsion Elements, 8th Edition, Equation 3-24.
func (laval Laval) MassFlow(pc float64, ex ExhaustIsoentropic) float64 {
	return (laval.areaThroat * pc * ex.Gamma / ex.SpeedOfSound(ex.MaxBurnTemperature)) *
		math.Pow(2/(ex.Gamma+1), (ex.Gamma+1)/(2*ex.Gamma-2))
}

// IsChoked returns true if the nozzle flow is choked (Mach 1 at throat).
func (laval Laval) IsChoked(pc, pe float64, ex ExhaustIsoentropic) bool {
	return pe/pc < math.Pow(2/(ex.Gamma+1), ex.Gamma/(ex.Gamma-1))
}

// SolvePressureRatioFromExpansion returns the exit-to-chamber pressure ratio pe/pc
// that produces the given area expansion ratio, using Newton-Raphson iteration.
func (laval Laval) SolvePressureRatioFromExpansion(expansion, gamma float64) (float64, bool) {
	x0 := 1e-3 / expansion
	residual := func(x float64) float64 {
		return expansion - laval.EstimateExpansionRatio(1, x, gamma)
	}
	solver := md1.DefaultNewtonRaphsonSolver()
	prat, converged := solver.Root(x0, residual)
	return prat, converged > 0
}

// ExitMach returns the supersonic exit Mach number using the Majdalani-Maickie
// series approximation (order 6).
//
// Reference: J. Majdalani and B. A. Maickie, http://maji.utsi.edu/publications/pdf/HT02_11.pdf
func (laval Laval) ExitMach(gamma float64) (float64, error) {
	const N = 6
	e := laval.areaThroat / laval.areaExit
	ga := gamma
	B := (ga + 1) / (ga - 1)
	k := math.Sqrt(.5 * (ga - 1))
	u := math.Pow(e, 1/B) / math.Sqrt(1+k*k)
	M := math.Pow(u*k, B/(1-B))

	k2, u2, B2 := k*k, u*u, B*B
	for i := 1; i < N; i++ {
		M2 := M * M
		lamb := 1 / (math.Pow(2*M, 2/B)*(B-2) + M2*B2*k2*u2)
		rt := math.Pow(M, 2+2/B)*k2*u2*(B2-4*B+4) - M2*B2*k2*u2*u2 +
			math.Pow(M, 4/B)*(2*B-3) + 2*math.Pow(M, 2/B)*u2*(2-B)
		if rt < 0 {
			return M, errors.New("imaginary part in exit Mach iteration")
		}
		M += lamb * M * B * (math.Pow(M, 2/B) - M2*B*k2*u2 + math.Sqrt(rt))
	}
	return M, nil
}

// ExhaustIsoentropic describes the thermodynamic properties of combustion exhaust.
type ExhaustIsoentropic struct {
	// Molar mass of exhaust gas [kg/mol]
	MolarMass float64
	// Specific heat ratio cp/cv
	Gamma float64
	// Stagnation (chamber) temperature [K]
	MaxBurnTemperature float64
}

// SpeedOfSound returns the isentropic speed of sound [m/s] at temperature T.
//
//	a = sqrt(γ·R·T),  R = Ru/M
func (ex ExhaustIsoentropic) SpeedOfSound(T float64) float64 {
	return math.Sqrt(ex.Gamma * _Ru / ex.MolarMass * T)
}
