package cosmos

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// IAU-76/FK5 Earth orientation: precession + 1980 nutation + apparent
// sidereal rotation, matching GMAT's reduction for EarthFixed axes (polar
// motion excluded: sub-arcsecond, ~10 m at the surface).
// References: Vallado, Fundamentals of Astrodynamics and Applications,
// sec. 3.7 (Eqns 3-68, 3-79, 3-82 to 3-88).

const arcsecToRad = math.Pi / (180 * 3600)

// precessionMOD returns the [MOD ← J2000] rotation: R3(-z)·R2(theta)·R3(-zeta)
// with the IAU-76 precession angles (Vallado Eqn 3-88).
func precessionMOD(tTT float64) md3.Mat3 {
	zeta := (2306.2181 + (0.30188+0.017998*tTT)*tTT) * tTT * arcsecToRad
	theta := (2004.3109 - (0.42665+0.041833*tTT)*tTT) * tTT * arcsecToRad
	z := (2306.2181 + (1.09468+0.018203*tTT)*tTT) * tTT * arcsecToRad
	return md3.MulMat3(rot3(-z), md3.MulMat3(rot2(theta), rot3(-zeta)))
}

// nutation1980Angles returns the IAU-1980 nutation in longitude and obliquity
// and the mean obliquity of the ecliptic, all in radians (Vallado Eqns 3-68,
// 3-82, 3-83; full 106-term series from GMAT's NUTATION.DAT).
func nutation1980Angles(tTT float64) (dPsi, dEps, epsBar float64) {
	const d2r = math.Pi / 180
	t := tTT
	// Mean obliquity, Vallado Eqn (3-68) [arcsec → rad].
	epsBar = (84381.448 - (46.8150+(0.00059-0.001813*t)*t)*t) * arcsecToRad
	// Delaunay fundamental arguments, Vallado Eqn (3-82) [deg → rad].
	r := 360.0
	l := 134.96298139 + (1325*r+198.8673981)*t + 0.0086972*t*t + 1.78e-5*t*t*t
	lp := 357.52772333 + (99*r+359.0503400)*t - 0.0001603*t*t - 3.3e-6*t*t*t
	f := 93.27191028 + (1342*r+82.0175381)*t - 0.0036825*t*t + 3.1e-6*t*t*t
	dd := 297.85036306 + (1236*r+307.1114800)*t - 0.0019142*t*t + 5.3e-6*t*t*t
	om := 125.04452222 - (5*r+134.1362608)*t + 0.0020708*t*t + 2.2e-6*t*t*t
	l, lp, f, dd, om = l*d2r, lp*d2r, f*d2r, dd*d2r, om*d2r

	for i := range nutation1980 {
		term := &nutation1980[i]
		arg := float64(term.l)*l + float64(term.lp)*lp + float64(term.f)*f +
			float64(term.d)*dd + float64(term.om)*om
		sin, cos := math.Sincos(arg)
		dPsi += (term.a + term.b*t) * sin
		dEps += (term.c + term.dD*t) * cos
	}
	// Series coefficients are in 0.0001 arcsec.
	dPsi *= 1e-4 * arcsecToRad
	dEps *= 1e-4 * arcsecToRad
	return dPsi, dEps, epsBar
}

// GAST returns the Greenwich apparent sidereal time [rad] in [0, 2pi):
// GMST plus the equation of the equinoxes (Vallado Eqn 3-79, including the
// post-1997 kinematic terms).
func (e Epoch) GAST() float64 {
	tTT := e.secsTT / (36525 * secsPerDay)
	dPsi, _, epsBar := nutation1980Angles(tTT)
	const d2r = math.Pi / 180
	om := (125.04452222 - (5*360+134.1362608)*tTT + 0.0020708*tTT*tTT + 2.2e-6*tTT*tTT*tTT) * d2r
	eqEquinox := dPsi*math.Cos(epsBar) + (0.00264*math.Sin(om)+0.000063*math.Sin(2*om))*arcsecToRad
	gast := math.Mod(e.GMST()+eqEquinox, 2*math.Pi)
	if gast < 0 {
		gast += 2 * math.Pi
	}
	return gast
}

// rot1, rot2, rot3 are frame rotation matrices about X, Y, Z (Vallado Eqn 3-15).
func rot1(a float64) md3.Mat3 {
	s, c := math.Sincos(a)
	return md3.NewMat3([]float64{
		1, 0, 0,
		0, c, s,
		0, -s, c,
	})
}

func rot2(a float64) md3.Mat3 {
	s, c := math.Sincos(a)
	return md3.NewMat3([]float64{
		c, 0, -s,
		0, 1, 0,
		s, 0, c,
	})
}

func rot3(a float64) md3.Mat3 {
	s, c := math.Sincos(a)
	return md3.NewMat3([]float64{
		c, s, 0,
		-s, c, 0,
		0, 0, 1,
	})
}
