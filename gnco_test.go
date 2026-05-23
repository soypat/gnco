package gnco_test

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md1"
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
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
		{"straight up", math.Pi / 2, 0},   // vertical flight, special case
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
