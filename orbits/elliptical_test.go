package orbits

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md1"
)

var (
	// universal gravitational constant - [N.m^2.kg^-2]
	bigG = 6.673e-11

	earthMass           = 5.973332e24
	earthc20            = -4.8416685e-4
	earthsemiMajorAxis  = 6378137                  // WGS84 [m]
	earthweii           = 7.292114999999999893e-05 // earth angular velocity in inertial coordinates [rad/s]
	earthradius         = 6370987.
	earthseaLevelRadius = 6371146
	earthflattening     = 3.33528106e-3
	earthcelestialLong  = 0
	earthGravParam      = bigG * earthMass
	earthDay            = 2 * math.Pi / earthweii
)

func TestElliptical_geocentric(t *testing.T) {
	el, err := NewElliptical(21_000e3, 9600e3)
	if err != nil {
		t.Fatal(err)
	}

	e := el.Eccentricity()
	if !md1.EqualWithinAbs(e, 0.37255, 0.0001) {
		t.Errorf("wanted %f, got %f", 0.37255, e)
	}
	const trueAnomaly = 120. * math.Pi / 180
	const wantEccAnomaly = 1.7281
	const wantMeanAnomaly = 1.3601
	const wantElapsed = 1.132 * 60 * 60
	ea := el.EccentricAnomaly(trueAnomaly)
	if !md1.EqualWithinAbs(ea, wantEccAnomaly, 0.001) {
		t.Errorf("wanted %f, got %f", wantEccAnomaly, ea)
	}

	Me := el.MeanAnomaly(trueAnomaly)
	if !md1.EqualWithinAbs(Me, wantMeanAnomaly, 0.001) {
		t.Errorf("wanted %f, got %f", wantMeanAnomaly, Me)
	}
	elapsed := el.ElapsedSincePeriapsis(earthGravParam, trueAnomaly)
	if !md1.EqualWithinAbs(elapsed, wantElapsed, 0.001*60*60) {
		t.Errorf("wanted %f, got %f", wantElapsed, elapsed)
	}
	const elapsed10800 = 10800.0
	gotTrueAnomaly := el.TrueAnomalyFromElapsedSincePeriapsis(earthGravParam, elapsed10800, 0.001)
	if roundTrip := el.ElapsedSincePeriapsis(earthGravParam, gotTrueAnomaly); !md1.EqualWithinAbs(roundTrip, elapsed10800, 0.001*60*60) {
		t.Errorf("round-trip at t=10800s: got ν=%.6f rad, ElapsedSincePeriapsis=%.1fs", gotTrueAnomaly, roundTrip)
	}
}

// TestTrueAnomalyRoundtrip verifies TrueAnomalyFromElapsedSincePeriapsis is the
// inverse of ElapsedSincePeriapsis across the full [0, 2π) range, including the
// descending half (ν > π) where math.Abs-based reconstruction produces wrong results.
func TestTrueAnomalyRoundtrip(t *testing.T) {
	el, err := NewElliptical(21_000e3, 9600e3)
	if err != nil {
		t.Fatal(err)
	}
	const tol = 1e-10
	cases := []struct {
		name string
		nu   float64
	}{
		{"periapsis", 0},
		{"ascending 60deg", 60 * math.Pi / 180},
		{"ascending 120deg", 120 * math.Pi / 180},
		{"apoapsis", math.Pi},
		{"descending 240deg", 240 * math.Pi / 180},
		{"descending 300deg", 300 * math.Pi / 180},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			elapsed := el.ElapsedSincePeriapsis(earthGravParam, tc.nu)
			got := el.TrueAnomalyFromElapsedSincePeriapsis(earthGravParam, elapsed, tol)
			if !md1.EqualWithinAbs(got, tc.nu, 1e-6) {
				t.Errorf("nu=%.4f rad: round-trip got %.6f", tc.nu, got)
			}
		})
	}
}

// TestTrueAnomalyMultiPeriod verifies elapsed times beyond one orbital period
// are handled correctly via modulo reduction.
func TestTrueAnomalyMultiPeriod(t *testing.T) {
	el, err := NewElliptical(21_000e3, 9600e3)
	if err != nil {
		t.Fatal(err)
	}
	const tol = 1e-10
	const nu = 120 * math.Pi / 180
	T := el.Period(earthGravParam)
	elapsed := el.ElapsedSincePeriapsis(earthGravParam, nu)
	for _, periods := range []float64{1, 2, 5} {
		got := el.TrueAnomalyFromElapsedSincePeriapsis(earthGravParam, elapsed+periods*T, tol)
		if !md1.EqualWithinAbs(got, nu, 1e-6) {
			t.Errorf("elapsed+%gT: want %.4f rad, got %.6f", periods, nu, got)
		}
	}
}
