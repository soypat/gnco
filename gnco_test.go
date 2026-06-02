package gnco_test

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md1"
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/orbits"
)

// tvgRow extracts row n of a TVG matrix via its transpose action on a basis vector.
// TVG^T * e_n  ==  n-th row of TVG.
func tvgRow(tvg md3.Mat3, n int) md3.Vec {
	e := [3]md3.Vec{{X: 1}, {Y: 1}, {Z: 1}}
	return md3.MulMatVecTrans(tvg, e[n])
}

// checkTVGOrthonormal asserts that tvg is a proper rotation matrix (orthonormal rows,
// right-hand system). tol is the absolute tolerance for each check.
func checkTVGOrthonormal(t *testing.T, tvg md3.Mat3, tol float64) {
	t.Helper()
	r0, r1, r2 := tvgRow(tvg, 0), tvgRow(tvg, 1), tvgRow(tvg, 2)

	for i, r := range []md3.Vec{r0, r1, r2} {
		if !md1.EqualWithinAbs(md3.Norm(r), 1, tol) {
			t.Errorf("row%d not unit length: |r|=%.6g", i, md3.Norm(r))
		}
	}
	for _, pair := range [][2]md3.Vec{{r0, r1}, {r0, r2}, {r1, r2}} {
		if !md1.EqualWithinAbs(md3.Dot(pair[0], pair[1]), 0, tol) {
			t.Errorf("rows not orthogonal: dot=%.6g", md3.Dot(pair[0], pair[1]))
		}
	}
	// Right-hand system: r1 × r2 == r0.
	diff := md3.Sub(md3.Cross(r1, r2), r0)
	if !md1.EqualWithinAbs(md3.Norm(diff), 0, tol) {
		t.Errorf("not right-hand system: r1×r2 - r0 = %v (norm %.6g)", diff, md3.Norm(diff))
	}
}

func TestTVGFromGeographicVelocity(t *testing.T) {
	const tol = 1e-12

	t.Run("zero velocity returns identity", func(t *testing.T) {
		tvg := gnco.TVGFromGeographicVelocity(md3.Vec{})
		for _, e := range []md3.Vec{{X: 1}, {Y: 1}, {Z: 1}} {
			got := md3.MulMatVec(tvg, e)
			if d := md3.Norm(md3.Sub(got, e)); d > tol {
				t.Errorf("TVG(0)*%v = %v, want identity", e, got)
			}
		}
	})

	// Cases defined by elevation and bearing, consistent with
	// GeographicVectorFromElevationAndBearing so the two functions are inverses.
	cases := []struct {
		name          string
		elev, bearing float64
	}{
		{"north level", 0, 0},
		{"east level", 0, math.Pi / 2},
		{"south level", 0, math.Pi},
		{"west level", 0, -math.Pi / 2},
		{"45deg up north", math.Pi / 4, 0},
		{"45deg up east", math.Pi / 4, math.Pi / 2},
		{"30deg up southwest", math.Pi / 6, 5 * math.Pi / 4},
		{"straight up", math.Pi / 2, 0},    // vertical flight, special case
		{"straight down", -math.Pi / 2, 0}, // vertical flight, special case
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			vbg := gnco.GeographicVectorFromElevationAndBearing(tc.elev, tc.bearing, 500)
			tvg := gnco.TVGFromGeographicVelocity(vbg)

			checkTVGOrthonormal(t, tvg, tol)

			// Round-trip: TVG^T * e_x (V-frame forward) must map back to the
			// normalised geographic velocity direction.
			want := md3.Unit(vbg)
			got := md3.MulMatVecTrans(tvg, md3.Vec{X: 1})
			if d := md3.Norm(md3.Sub(got, want)); d > tol {
				t.Errorf("TVG^T*e_x = %v, want %v (diff %.3g)", got, want, d)
			}
		})
	}

	t.Run("scale invariance", func(t *testing.T) {
		// Multiplying velocity magnitude should not change TVG.
		vbg := gnco.GeographicVectorFromElevationAndBearing(0.3, 1.1, 1)
		tvg1 := gnco.TVGFromGeographicVelocity(vbg)
		tvg2 := gnco.TVGFromGeographicVelocity(md3.Scale(1000, vbg))
		for n := 0; n < 3; n++ {
			diff := md3.Sub(tvgRow(tvg1, n), tvgRow(tvg2, n))
			if !md1.EqualWithinAbs(md3.Norm(diff), 0, tol) {
				t.Errorf("row%d differs under scaling: %v vs %v", n, tvgRow(tvg1, n), tvgRow(tvg2, n))
			}
		}
	})
}

func TestWorldHASL(t *testing.T) {
	earth := gnco.NewEarth()

	// Round-trip: HASL(HASLToElevation(h)) == h for representative altitudes.
	for _, h := range []float64{0, 100, 10_000, 400_000, 35_786_000} {
		elev := earth.HASLToElevation(h)
		got := earth.HASL(elev)
		if !md1.EqualWithinAbs(got, h, 1e-6) {
			t.Errorf("HASL(HASLToElevation(%g)) = %g, want %g", h, got, h)
		}
	}

	// Negative altitudes (underground / impact scenarios) round-trip correctly.
	for _, h := range []float64{-50, -1000} {
		elev := earth.HASLToElevation(h)
		got := earth.HASL(elev)
		if !md1.EqualWithinAbs(got, h, 1e-6) {
			t.Errorf("HASL(HASLToElevation(%g)) = %g, want %g", h, got, h)
		}
	}
}

func TestTEI(t *testing.T) {
	const tol = 1e-12
	earth := gnco.NewEarth()

	t.Run("identity at t=0", func(t *testing.T) {
		tei := earth.TEI(0)
		for _, e := range []md3.Vec{{X: 1}, {Y: 1}, {Z: 1}} {
			got := md3.MulMatVec(tei, e)
			if d := md3.Norm(md3.Sub(got, e)); d > tol {
				t.Errorf("TEI(0)*%v = %v, want identity", e, got)
			}
		}
	})

	t.Run("identity after one sidereal day", func(t *testing.T) {
		tei := earth.TEI(earth.Day())
		for _, e := range []md3.Vec{{X: 1}, {Y: 1}, {Z: 1}} {
			got := md3.MulMatVec(tei, e)
			if d := md3.Norm(md3.Sub(got, e)); d > tol {
				t.Errorf("TEI(day)*%v = %v, want identity", e, got)
			}
		}
	})

	// After a quarter rotation, the ECEF frame has rotated π/2 CCW relative to ECI.
	// An ECI +X point is therefore at ECEF -Y (Earth has moved past it).
	t.Run("quarter rotation ECI+X maps to ECEF -Y", func(t *testing.T) {
		quarterDay := math.Pi / 2 / earth.Rotation()
		tei := earth.TEI(quarterDay)
		got := md3.MulMatVec(tei, md3.Vec{X: 1})
		want := md3.Vec{Y: -1}
		if d := md3.Norm(md3.Sub(got, want)); d > tol {
			t.Errorf("TEI(quarterDay)*{X:1} = %v, want %v", got, want)
		}
	})

	// Inverse: an ECEF +X surface point appears at ECI +Y after a quarter rotation.
	t.Run("quarter rotation ECEF+X maps to ECI +Y", func(t *testing.T) {
		quarterDay := math.Pi / 2 / earth.Rotation()
		tei := earth.TEI(quarterDay)
		got := md3.MulMatVecTrans(tei, md3.Vec{X: 1})
		want := md3.Vec{Y: 1}
		if d := md3.Norm(md3.Sub(got, want)); d > tol {
			t.Errorf("TEI(quarterDay)^T*{X:1} = %v, want %v", got, want)
		}
	})

	t.Run("orthonormal at various times", func(t *testing.T) {
		for _, tt := range []float64{0, 100, 3600, 86400} {
			tei := earth.TEI(tt)
			// TEI * TEI^T = I: check each basis vector round-trips.
			for _, e := range []md3.Vec{{X: 1}, {Y: 1}, {Z: 1}} {
				v := md3.MulMatVecTrans(tei, e)
				got := md3.MulMatVec(tei, v)
				if d := md3.Norm(md3.Sub(got, e)); d > tol {
					t.Errorf("TEI(%g)*TEI^T*%v = %v, want %v", tt, e, got, e)
				}
			}
		}
	})
}

func TestGeocentricECIRoundTrip(t *testing.T) {
	const angTol = 1e-10 // rad (~1 mm on Earth's surface)
	const elevTol = 1e-3 // m
	earth := gnco.NewEarth()

	cases := []struct {
		name      string
		longDeg   float64
		latDeg    float64
		hasl      float64
		epochTime float64
	}{
		{"equator t=0", 0, 0, 400e3, 0},
		{"equator t=1h", 0, 0, 400e3, 3600},
		{"NE quadrant t=1h", 45, 30, 200e3, 3600},
		{"SW quadrant t=30min", -60, -45, 100e3, 1800},
		{"prime meridian t=1day", 0, 60, 0, 86400},
		{"high latitude t=6h", 15, 75, 500e3, 6 * 3600},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			g := earth.GeocentricFromDegrees(tc.longDeg, tc.latDeg, earth.HASLToElevation(tc.hasl))
			// EarthFixedCoords returns the ECI position of this surface-fixed point at time t.
			sBII := g.EarthFixedCoords(tc.epochTime)
			// GeocentricFromEarthFixedCoords inverts: given ECI position + epoch, recover ECEF coords.
			g2 := earth.GeocentricFromEarthFixedCoords(sBII, tc.epochTime)
			if d := math.Abs(g2.Long - g.Long); d > angTol {
				t.Errorf("Long round-trip: got %g, want %g (diff %g)", g2.Long, g.Long, d)
			}
			if d := math.Abs(g2.Lat - g.Lat); d > angTol {
				t.Errorf("Lat round-trip: got %g, want %g (diff %g)", g2.Lat, g.Lat, d)
			}
			if d := math.Abs(g2.Elev - g.Elev); d > elevTol {
				t.Errorf("Elev round-trip: got %g, want %g (diff %g m)", g2.Elev, g.Elev, d)
			}
		})
	}
}

func TestPhysicsKeplerianEnergy(t *testing.T) {
	const (
		perigeeHASL = 400e3 // m
		apogeeHASL  = 500e3 // m
		dt          = 300.0 // s
		nOrbits     = 5
		energyTol   = 1.5e-14
	)
	earth := gnco.NewEarth()
	mu := earth.G()
	rP := earth.Radius() + earth.HASLToElevation(perigeeHASL)
	rA := earth.Radius() + earth.HASLToElevation(apogeeHASL)
	orbit, err := orbits.NewElliptical(rA, rP)
	if err != nil {
		t.Fatal(err)
	}

	_, vT := orbit.Velocity(mu, 0) // vRadial=0 at periapsis; vT is tangential velocity
	SBI0 := md3.Vec{X: rP}
	VBI0 := md3.Vec{Y: vT}
	E0 := orbit.SpecificEnergy(mu)
	T := orbit.Period(mu)

	coords := earth.GeocentricFromEarthFixedCoords(SBI0, 0)
	integrator := gnco.NewPhysicsPointIntegrator(&coords, 0, SBI0, VBI0)

	tt, SBI, VBI := integrator.State()
	var maxErrE float64
	for tt < float64(nOrbits)*T {
		r := md3.Norm(SBI)
		v := md3.Norm(VBI)
		E := 0.5*v*v - mu/r
		if errE := math.Abs((E - E0) / E0); errE > maxErrE {
			maxErrE = errE
		}
		tt, SBI, VBI = integrator.Step(dt, md3.Vec{})
	}
	if maxErrE > energyTol {
		t.Errorf("max |ΔE/E₀| = %.2e over %d orbits, want < %.2e", maxErrE, nOrbits, energyTol)
	}
}

func TestCoordsHASL(t *testing.T) {
	earth := gnco.NewEarth()

	for _, hasl := range []float64{0, 100, 10_000, 400_000} {
		c := earth.GeocentricFromDegrees(45, 90, earth.HASLToElevation(hasl))

		// GeocentricCoords.HASL
		if got := c.HASL(); !md1.EqualWithinAbs(got, hasl, 1e-6) {
			t.Errorf("GeocentricCoords.HASL() = %g, want %g", got, hasl)
		}

		// GeodesicCoords.HASL delegates to the geocentric value.
		if got := c.Geodesic().HASL(); !md1.EqualWithinAbs(got, hasl, 1e-6) {
			t.Errorf("GeodesicCoords.HASL() = %g, want %g", got, hasl)
		}
	}
}

func BenchmarkPhysicsPointIntegrator(b *testing.B) {
	const (
		perigeeHASL = 500e3
		apogeeHASL  = 1000e3
		dt          = 100.0
	)
	earth := gnco.NewEarth()
	mu := earth.G()
	rP := earth.Radius() + earth.HASLToElevation(perigeeHASL)
	rA := earth.Radius() + earth.HASLToElevation(apogeeHASL)
	orbit, err := orbits.NewElliptical(rA, rP)
	if err != nil {
		b.Fatal(err)
	}
	_, vT := orbit.Velocity(mu, 0)
	SBI0 := md3.Vec{X: rP}
	VBI0 := md3.Vec{Y: vT}
	E0 := orbit.SpecificEnergy(mu)

	b.Run("Step", func(b *testing.B) {
		coords := earth.GeocentricFromEarthFixedCoords(SBI0, 0)
		integrator := gnco.NewPhysicsPointIntegrator(&coords, 0, SBI0, VBI0)
		var maxErrE float64
		for b.Loop() {
			_, SBI, VBI := integrator.Step(dt, md3.Vec{})
			r, v := md3.Norm(SBI), md3.Norm(VBI)
			E := 0.5*v*v - mu/r
			if errE := math.Abs((E - E0) / E0); errE > maxErrE {
				maxErrE = errE
			}
		}
		b.ReportMetric(maxErrE, "|ΔE/E₀|")
	})
	b.Run("StepFast", func(b *testing.B) {
		coords := earth.GeocentricFromEarthFixedCoords(SBI0, 0)
		integrator := gnco.NewPhysicsPointIntegrator(&coords, 0, SBI0, VBI0)
		var maxErrE float64
		for b.Loop() {
			_, SBI, VBI := integrator.StepFast(dt, md3.Vec{})
			r, v := md3.Norm(SBI), md3.Norm(VBI)
			E := 0.5*v*v - mu/r
			if errE := math.Abs((E - E0) / E0); errE > maxErrE {
				maxErrE = errE
			}
		}
		b.ReportMetric(maxErrE, "|ΔE/E₀|")
	})
}
