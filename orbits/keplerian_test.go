package orbits

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md3"
)

const muEarthVallado = 398600.4418e9 // [m³/s²]

// Vallado, Example 2-6 (COE2RV): p=11067.790 km, e=0.83285, i=87.87°,
// raan=227.89°, aop=53.38°, ta=92.335°.
func TestKeplerianRVVallado(t *testing.T) {
	const d = math.Pi / 180
	const p, e = 11067.790e3, 0.83285
	sma := p / (1 - e*e)
	k, err := NewKeplerian(sma, e, 87.87*d, 227.89*d, 53.38*d)
	if err != nil {
		t.Fatal(err)
	}
	r, v := k.RV(muEarthVallado, 92.335*d)
	// Published result (Vallado 4th ed.): differs from the Ex 2-5 input
	// vector by ~25 m because the elements are rounded between examples.
	wantR := md3.Vec{X: 6525.368e3, Y: 6861.532e3, Z: 6449.119e3}
	wantV := md3.Vec{X: 4.902279e3, Y: 5.533140e3, Z: -1.975710e3}
	if dr := md3.Norm(md3.Sub(r, wantR)); dr > 1.0 {
		t.Errorf("position error %.3f m: got %v want %v", dr, r, wantR)
	}
	if dv := md3.Norm(md3.Sub(v, wantV)); dv > 1e-3 {
		t.Errorf("velocity error %.6f m/s: got %v want %v", dv, v, wantV)
	}
}

// Vallado, Example 2-5 (RV2COE): r=(6524.834, 6862.875, 6448.296) km,
// v=(4.901327, 5.533756, -1.976341) km/s.
func TestKeplerianFromRVVallado(t *testing.T) {
	const d = math.Pi / 180
	r := md3.Vec{X: 6524.834e3, Y: 6862.875e3, Z: 6448.296e3}
	v := md3.Vec{X: 4.901327e3, Y: 5.533756e3, Z: -1.976341e3}
	k, ta, err := KeplerianFromRV(muEarthVallado, r, v)
	if err != nil {
		t.Fatal(err)
	}
	checks := []struct {
		name, unit string
		got, want  float64
		tol        float64
	}{
		{"sma", "km", k.SemiMajorAxis() / 1e3, 36127.343, 0.2},
		{"ecc", "", k.Eccentricity(), 0.832853, 1e-5},
		{"inc", "deg", k.Inclination() / d, 87.870, 1e-3},
		{"raan", "deg", k.RAAN() / d, 227.898, 1e-3},
		{"aop", "deg", k.ArgumentOfPeriapsis() / d, 53.38, 1e-2},
		{"ta", "deg", ta / d, 92.335, 1e-2},
	}
	for _, c := range checks {
		if math.Abs(c.got-c.want) > c.tol {
			t.Errorf("%s = %.6f %s, want %.6f", c.name, c.got, c.unit, c.want)
		}
	}
}

// Round trip with the SolarCalc orbit elements.
func TestKeplerianRoundTrip(t *testing.T) {
	const d = math.Pi / 180
	const mu = 3.986004415e14
	k, err := NewKeplerian(6928.5e3, 0.0011013, 97.794*d, 215.84*d, 191.1*d)
	if err != nil {
		t.Fatal(err)
	}
	for _, taDeg := range []float64{0, 45, 118.33, 179.99, 251, 359} {
		r, v := k.RV(mu, taDeg*d)
		k2, ta2, err := KeplerianFromRV(mu, r, v)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(k2.SemiMajorAxis()-k.SemiMajorAxis()) > 1e-3 ||
			math.Abs(k2.Eccentricity()-k.Eccentricity()) > 1e-9 ||
			math.Abs(k2.Inclination()-k.Inclination()) > 1e-9 ||
			math.Abs(k2.RAAN()-k.RAAN()) > 1e-9 ||
			math.Abs(k2.ArgumentOfPeriapsis()-k.ArgumentOfPeriapsis()) > 1e-6 ||
			math.Abs(mod2pi(ta2-taDeg*d)) > 1e-6 && math.Abs(mod2pi(ta2-taDeg*d)-2*math.Pi) > 1e-6 {
			t.Errorf("round trip mismatch at ta=%g°: %+v ta=%g", taDeg, k2, ta2/d)
		}
	}
	// Energy/period cross-check against Elliptical bridge.
	ell, err := k.Elliptical()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ell.Period(mu)-k.Period(mu)) > 1e-6 {
		t.Error("Period mismatch between Keplerian and Elliptical")
	}
}
