package nozzles_test

import (
	"math"
	"testing"

	"github.com/soypat/gnco/rocket/nozzles"
)

// keroseneOx is a representative LOX/kerosene exhaust for testing.
var keroseneOx = nozzles.DefaultLOXKeroseneExhaust()

// TestThrustIspMassFlowIdentity checks T = ṁ · Isp · g₀ across a range of
// chamber pressures. The identity follows from the definitions of Isp and Cf
// and must hold exactly regardless of operating point.
func TestThrustIspMassFlowIdentity(t *testing.T) {
	const g0 = 9.80665
	noz := nozzles.NewLaval(0.01, 0.16)
	ex := keroseneOx
	pe := 10e3  // [Pa] near-vacuum exit
	pa := 101e3 // [Pa] sea level

	for _, pc := range []float64{2e6, 5e6, 10e6, 20e6} {
		thrust := noz.Thrust(pc, pe, pa, ex.Gamma)
		mdot := noz.MassFlow(pc, ex)
		isp := noz.SpecificImpulse(pc, pe, pa, ex)

		got := mdot * isp * g0
		if !relClose(got, thrust, 1e-9) {
			t.Errorf("pc=%.0f Pa: T=%.4f N, ṁ·Isp·g₀=%.4f N (rel err %.2e)",
				pc, thrust, got, math.Abs(got-thrust)/thrust)
		}
	}
}

// TestExpansionRatioRoundTrip checks that SolvePressureRatioFromExpansion is
// the numerical inverse of EstimateExpansionRatio: feeding the output of one
// into the other must recover the original pressure ratio.
func TestExpansionRatioRoundTrip(t *testing.T) {
	noz := nozzles.NewLaval(0.005, 0.08)
	gamma := keroseneOx.Gamma
	pc := 8e6

	for _, peFrac := range []float64{0.005, 0.01, 0.05, 0.1, 0.2} {
		pe := pc * peFrac
		er := noz.EstimateExpansionRatio(pc, pe, gamma)

		recovered, ok := noz.SolvePressureRatioFromExpansion(er, gamma)
		if !ok {
			t.Fatalf("pe/pc=%.3f: solver did not converge", peFrac)
		}
		if !relClose(recovered, peFrac, 1e-6) {
			t.Errorf("pe/pc=%.4f: recovered=%.8f (rel err %.2e)",
				peFrac, recovered, math.Abs(recovered-peFrac)/peFrac)
		}
	}
}

// TestExitMachAreaMachConsistency checks that ExitMach satisfies the isentropic
// area-Mach relation independently of how it was computed:
//
//	Ae/At = (1/Me) · [(2/(γ+1))·(1 + (γ-1)/2·Me²)]^((γ+1)/(2(γ-1)))
func TestExitMachAreaMachConsistency(t *testing.T) {
	gamma := keroseneOx.Gamma

	for _, ratio := range []float64{4, 8, 16, 32} {
		at := 0.01
		noz := nozzles.NewLaval(at, at*ratio)

		Me, err := noz.ExitMach(gamma)
		if err != nil {
			t.Fatalf("Ae/At=%.0f: ExitMach error: %v", ratio, err)
		}

		erFromMach := areaMachRatio(Me, gamma)
		if !relClose(erFromMach, ratio, 1e-4) {
			t.Errorf("Ae/At=%.0f: area-Mach gives %.6f (rel err %.2e)",
				ratio, erFromMach, math.Abs(erFromMach-ratio)/ratio)
		}
	}
}

// areaMachRatio computes Ae/At from exit Mach number via the isentropic relation.
func areaMachRatio(Me, gamma float64) float64 {
	exp := (gamma + 1) / (2 * (gamma - 1))
	inner := 2 / (gamma + 1) * (1 + (gamma-1)/2*Me*Me)
	return math.Pow(inner, exp) / Me
}

func relClose(a, b, tol float64) bool {
	if b == 0 {
		return math.Abs(a) < tol
	}
	return math.Abs(a-b)/math.Abs(b) <= tol
}
