package nozzle

import "math"

const (
	//Boltzmann's constant [W.s/K]
	boltzmann = 1.38e-23
	// Avogadro's constant [molecules/mol]
	avogadro = 6.0221409e+23
	// Universal gas constant [m3⋅Pa⋅K−1⋅mol−1] also [kg⋅m^2⋅s^−2⋅K^−1⋅mol^−1]
	_Ru = 8.31446261815324
)

// Laval defines a Laval nozzle.
// Reference: Rocket Propulsion Elements, 8th Edition, Equation
type Laval struct {
	areaThroat float64 // [m^2]
	areaExit   float64 // [m^2]
}

// EstimateExpansionRatio estimates the expansion ratio from exhaust
// specific heat ratio (gamma=cp/cv ratio) and chamber/exit pressure.
func (laval Laval) EstimateExpansionRatio(pc, pe, exhaustGamma float64) float64 {
	AtdAe := math.Pow((exhaustGamma+1)/2, 1/(exhaustGamma-1)) *
		math.Pow(pe/pc, 1/exhaustGamma) *
		math.Sqrt((exhaustGamma+1)/(exhaustGamma-1)*(1-math.Pow(pe/pc, (exhaustGamma-1)/exhaustGamma)))
	return 1 / AtdAe
}

func (laval Laval) Cf(pc, pe, pa, exhaustGamma float64) float64 {
	if pc <= pe || pc <= pa {
		// avoid NaN values. Thrust is zero for negative relative pressures.
		return 0
	}
	ga := exhaustGamma // short hand gamma declaration
	res := 2 * ga * ga / (ga - 1) * math.Pow(2/(ga+1), (ga+1)/(ga-1)) *
		(1 - math.Pow(pe/pc, (ga-1)/ga))
	er := laval.EstimateExpansionRatio(pc, pe, exhaustGamma)
	return math.Sqrt(res) + (pe-pa)*er/pc
}

type ExhaustIsoentropic struct {
	// Molar mass of exhaust gas [kg/mol]
	MolarMass float64
	// Specific heat ratio cp/cv
	Gamma float64
	// Maximum burn temperature [K]
	MaxBurnTemperature float64
}

// MassFlow in [kg/s] from chamber pressure pc and given
// stagnation conditions.
//
// Reference: Rocket Propulsion Elements, 8th Edition, Equation 3-24.
func (n Laval) MassFlow(pc float64, ex ExhaustIsoentropic) float64 {
	return (n.areaThroat * pc * ex.Gamma / math.Sqrt(ex.Gamma*_Ru/ex.MolarMass*ex.MaxBurnTemperature)) * math.Pow(2/(ex.Gamma+1), (ex.Gamma+1)/(2*ex.Gamma-2))
}
