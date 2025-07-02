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
	const wantTA = 193.2 * math.Pi / 180
	gotTrueAnomaly := el.TrueAnomalyFromElapsedSincePeriapsis(earthGravParam, 10800, 0.001)
	if !md1.EqualWithinAbs(gotTrueAnomaly, wantTA, 0.01) {
		t.Errorf("wanted %f, got %f", wantTA, gotTrueAnomaly)
	}
}
