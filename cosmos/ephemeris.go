package cosmos

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// AU is the astronomical unit [m] (GMAT/DE405 value).
const AU = 1.49597870691e11

// Ephemeris yields a body's position in the central body's MJ2000Eq inertial
// frame [m] at a given epoch.
type Ephemeris interface {
	Position(e Epoch) md3.Vec
}

// NewSun returns the Sun with GMAT-compatible constants (point-mass
// perturber and eclipse occultation geometry; most Body fields do not apply).
func NewSun() *Body {
	return &Body{
		name:          "Sun",
		mu:            1.32712440018e20, // GMAT Sun.Mu [m³/s²]
		mass:          1.9891e30,
		radius:        695990e3, // GMAT Sun.EquatorialRadius [m]
		semiMajorAxis: 695990e3,
	}
}

// AnalyticSun is the low-precision analytic solar ephemeris of Vallado,
// Algorithm 29 (sec. 5.1): apparent accuracy ~0.01° (36″), ample for beta
// angle and eclipse timing (36″ ≈ 0.16 s of LEO eclipse boundary timing).
// The series yields a mean-of-date vector; it is rotated by the transpose of
// the IAU-76 precession into MJ2000Eq. Light-time and stellar aberration
// (~20″, GMAT EclipseLocator applies them) are not modeled.
type AnalyticSun struct{}

// NewAnalyticSun returns the Vallado Alg. 29 solar ephemeris.
func NewAnalyticSun() *AnalyticSun { return &AnalyticSun{} }

// Position returns the geocentric Sun position in MJ2000Eq [m].
func (*AnalyticSun) Position(e Epoch) md3.Vec {
	const d2r = math.Pi / 180
	t := (e.JulianDateUT1() - jdJ2000) / 36525
	// Mean longitude and mean anomaly of the Sun [deg].
	lamM := 280.460 + 36000.771*t
	M := (357.5291092 + 35999.05034*t) * d2r
	sinM, cosM := math.Sincos(M)
	sin2M, cos2M := math.Sincos(2 * M)
	// Ecliptic longitude (apparent, low precision) and distance.
	lam := (lamM+1.914666471*sinM+0.019994643*sin2M)*d2r + 0 // [rad]
	rmag := (1.000140612 - 0.016708617*cosM - 0.000139589*cos2M) * AU
	eps := (23.439291 - 0.0130042*t) * d2r
	sinL, cosL := math.Sincos(lam)
	sinE, cosE := math.Sincos(eps)
	rMOD := md3.Vec{
		X: rmag * cosL,
		Y: rmag * cosE * sinL,
		Z: rmag * sinE * sinL,
	}
	// Mean-of-date → J2000 (inverse precession). Nutation (≤17″) is below
	// the series accuracy and is not applied.
	tTT := e.secsTT / (36525 * secsPerDay)
	return md3.MulMatVecTrans(precessionMOD(tTT), rMOD)
}

// Shadow returns the angular eclipse margins [rad] of a point at inertial
// position posBI [m] occulted by the central body (sphere of occRadius [m]
// at the origin) with the Sun at sunBI [m], per the apparent-radius cone
// geometry (Vallado sec. 5.3 / GMAT EclipseLocator):
//
//	a = asin(sunRadius/|sun-pos|)   apparent Sun radius
//	b = asin(occRadius/|pos|)       apparent body radius
//	c = angle between (-pos) and (sun-pos)
//
// penumbra = c-(a+b): negative while any sunlight is blocked.
// umbra = c-(b-a): negative while the Sun is fully occulted.
func Shadow(sunBI, posBI md3.Vec, sunRadius, occRadius float64) (penumbra, umbra float64) {
	sunRel := md3.Sub(sunBI, posBI)
	dSun := md3.Norm(sunRel)
	dOcc := md3.Norm(posBI)
	a := math.Asin(math.Min(1, sunRadius/dSun))
	b := math.Asin(math.Min(1, occRadius/dOcc))
	cosc := -md3.Dot(posBI, sunRel) / (dOcc * dSun)
	c := math.Acos(math.Max(-1, math.Min(1, cosc)))
	return c - (a + b), c - (b - a)
}
