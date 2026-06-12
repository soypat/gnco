// Package cosmos provides celestial body constants, astrodynamic time scales
// and gravity field models for orbit propagation.
package cosmos

import (
	"github.com/soypat/geometry/md3"
)

// Body describes a celestial body's physical constants for orbit propagation.
// Unlike gnco.World, Body carries the gravitational parameter mu directly since
// high-precision gravity models (JGM-2, GMAT's Earth.Mu) do not factor cleanly
// into G*mass.
type Body struct {
	name           string
	mu             float64 // gravitational parameter GM [m³/s²]
	mass           float64 // total mass of body [kg]
	c20            float64 // second degree zonal gravitational coefficient [Adim]
	semiMajorAxis  float64 // semi-major axis of body [m], also known as "equatorial radius"
	rotation       float64 // Angular rotation of body w.r.t inertial frame [rad/s]
	radius         float64 // Radius of body [m]
	seaLevelRadius float64 // If body stopped rotating the sea level would take this distance from center [m] https://www.esri.com/news/arcuser/0703/geoid3of3.html
	flattening     float64 // Flattening of body, (WGS84 for earth) [Adim]
	celestialLong  float64 // Celestial longitude. For earth is Greenwich meridian. Will indicate start of epoch [rad]

	// SGP4 parameters:
	ke float64 // kₑ: square root of gravitational parameter in body radii³ min⁻²
	j2 float64 // J₂: un-normalised second zonal harmonic
	j3 float64 // J₃: un-normalised third zonal harmonic
	j4 float64 // J₄: un-normalised fourth zonal harmonic
}

// NewEarth returns Earth with GMAT-compatible gravitational parameter
// (mu = 3.986004415e14 m³/s², GMAT Earth.Mu, identical to the JGM-2 potential
// file value) and WGS84 constants matching gnco.NewEarth for remaining fields.
func NewEarth() *Body {
	return &Body{
		name:           "Earth",
		mu:             3.986004415e14,
		mass:           5.973332e24,
		c20:            -4.8416685e-4,
		semiMajorAxis:  6378137, // WGS84 [m]
		rotation:       7.292114999999999893e-05,
		radius:         6370987.,
		seaLevelRadius: 6371146,
		flattening:     3.33528106e-3,
		celestialLong:  0,

		// SGP4 according to WGS84. Recommended by IAU to propagate orbits
		ke: 0.07436685316871385,
		j2: 0.00108262998905,
		j3: -0.00000253215306,
		j4: -0.00000161098761,
	}
}

// Name returns the body's name, i.e. "Earth".
func (b *Body) Name() string { return b.name }

// Mu returns the gravitational parameter GM [m³/s²].
func (b *Body) Mu() float64 { return b.mu }

// Mass returns the total mass of the body [kg].
func (b *Body) Mass() float64 { return b.mass }

// C20 returns the second degree zonal gravitational coefficient (dimensionless).
func (b *Body) C20() float64 { return b.c20 }

// SemiMajorAxis returns the equatorial radius [m] (WGS84 semi-major axis for Earth).
func (b *Body) SemiMajorAxis() float64 { return b.semiMajorAxis }

// Rotation returns the angular rotation rate of the body w.r.t the inertial frame [rad/s].
func (b *Body) Rotation() float64 { return b.rotation }

// Radius returns the mean radius of the body [m].
func (b *Body) Radius() float64 { return b.radius }

// Flattening returns the body's ellipsoidal flattening (dimensionless, WGS84 for Earth).
func (b *Body) Flattening() float64 { return b.flattening }

// Ke returns the SGP4 square root of the body's gravitational parameter [body radii³ min⁻²].
func (b *Body) Ke() float64 { return b.ke }

// J2 returns the SGP4 un-normalised second zonal harmonic.
func (b *Body) J2() float64 { return b.j2 }

// J3 returns the SGP4 un-normalised third zonal harmonic.
func (b *Body) J3() float64 { return b.j3 }

// J4 returns the SGP4 un-normalised fourth zonal harmonic.
func (b *Body) J4() float64 { return b.j4 }

// TEI returns the [T]^{EI} body-fixed ← inertial (MJ2000Eq) rotation tensor
// at epoch e using the IAU-76/FK5 reduction (Vallado sec. 3.7), as GMAT does
// for its EarthFixed axes:
//
//	TEI = R3(GAST) · [nutation 1980] · [precession IAU-76]
//
// Polar motion is excluded (sub-arcsecond, ~10 m at the surface; GMAT reads
// it from measured EOP data). The remaining gap to GMAT is the package-wide
// ΔUT1 = 0 assumption inside GMST (≤0.9 s of rotation ≈ ≤420 m equatorial
// displacement of the body-fixed frame, |2026 values| ≈ 0.1 s).
func (b *Body) TEI(e Epoch) md3.Mat3 {
	tTT := e.secsTT / (36525 * secsPerDay)
	prec := precessionMOD(tTT)
	dPsi, dEps, epsBar := nutation1980Angles(tTT)
	nut := md3.MulMat3(rot1(-(epsBar + dEps)), md3.MulMat3(rot3(-dPsi), rot1(epsBar)))
	return md3.MulMat3(rot3(e.GAST()+b.celestialLong), md3.MulMat3(nut, prec))
}
