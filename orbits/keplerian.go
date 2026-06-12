package orbits

import (
	"fmt"
	"math"

	"github.com/soypat/geometry/md3"
)

// Keplerian is the classical osculating orbit geometry about a central body:
// the in-plane shape (semi-major axis, eccentricity) plus the 3D orientation
// (inclination, RAAN, argument of periapsis) w.r.t an equatorial inertial
// frame (GMAT's EarthMJ2000Eq). Like Elliptical, the position along the orbit
// is not part of the type — methods take a true anomaly argument.
// Conversions per Vallado, Fundamentals of Astrodynamics and Applications,
// Algorithms 9 (RV2COE) and 10 (COE2RV). Parabolic and hyperbolic orbits
// (eccentricity >= 1) are not supported, as with Elliptical.
type Keplerian struct {
	sma  float64 // semi-major axis [m]
	ecc  float64 // eccentricity [adim]
	inc  float64 // inclination [rad]
	raan float64 // right ascension of the ascending node [rad]
	aop  float64 // argument of periapsis [rad]
}

// NewKeplerian creates an orbit geometry. semiMajorAxis in meters, angles in
// radians. inclination must lie in [0, pi]; raan and argPeriapsis are
// normalized into [0, 2pi).
func NewKeplerian(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis float64) (Keplerian, error) {
	switch {
	case math.IsNaN(semiMajorAxis + eccentricity + inclination + raan + argPeriapsis):
		return Keplerian{}, fmt.Errorf("NaN argument to NewKeplerian")
	case semiMajorAxis <= 0:
		return Keplerian{}, fmt.Errorf("bad semi-major axis %.5g km, must be positive", semiMajorAxis/1e3)
	case eccentricity < 0 || eccentricity >= 1:
		return Keplerian{}, fmt.Errorf("bad eccentricity %.3g, parabolic/hyperbolic orbits not supported", eccentricity)
	case inclination < 0 || inclination > math.Pi:
		return Keplerian{}, fmt.Errorf("bad inclination %.3g rad, must be in [0, pi]", inclination)
	}
	return Keplerian{
		sma:  semiMajorAxis,
		ecc:  eccentricity,
		inc:  inclination,
		raan: mod2pi(raan),
		aop:  mod2pi(argPeriapsis),
	}, nil
}

// KeplerianFromRV computes the osculating orbit geometry and the true anomaly
// at which the body sits on it from an inertial equatorial position r [m] and
// velocity v [m/s], given the gravitational parameter mu [m³/s²].
// Vallado Alg. 9 (RV2COE).
//
// Degenerate geometries follow the Vallado/GMAT conventions: for equatorial
// orbits raan=0 and the node is taken along +X; for circular orbits aop=0 and
// the true anomaly is measured from the node (argument of latitude).
func KeplerianFromRV(mu float64, r, v md3.Vec) (k Keplerian, trueAnomaly float64, err error) {
	const tol = 1e-11 // dimensionless threshold for degenerate geometry
	rn, vn := md3.Norm(r), md3.Norm(v)
	if mu <= 0 || rn == 0 || vn == 0 {
		return Keplerian{}, 0, fmt.Errorf("bad arguments to KeplerianFromRV: mu=%.3g r=%.3g v=%.3g", mu, rn, vn)
	}
	h := md3.Cross(r, v)
	hn := md3.Norm(h)
	if hn <= tol*rn*vn {
		return Keplerian{}, 0, fmt.Errorf("degenerate rectilinear orbit: r and v are parallel")
	}
	n := md3.Vec{X: -h.Y, Y: h.X} // node vector k×h
	nn := md3.Norm(n)
	// Eccentricity vector and energy.
	rdotv := md3.Dot(r, v)
	evec := md3.Sub(md3.Scale(vn*vn-mu/rn, r), md3.Scale(rdotv, v))
	evec = md3.Scale(1/mu, evec)
	ecc := md3.Norm(evec)
	if ecc >= 1 {
		return Keplerian{}, 0, fmt.Errorf("bad orbit eccentricity, got %.3g. Orbit does not expect parabolic/hyperbolic orbits", ecc)
	}
	xi := vn*vn/2 - mu/rn
	sma := -mu / (2 * xi)
	inc := math.Acos(clamp1(h.Z / hn))

	circular := ecc < tol
	equatorial := nn < tol*hn

	var raan, aop, ta float64
	switch {
	case equatorial && circular:
		// True longitude measured from +X.
		ta = math.Acos(clamp1(r.X / rn))
		if r.Y < 0 {
			ta = 2*math.Pi - ta
		}
	case equatorial: // elliptical equatorial: aop is the true longitude of periapsis.
		aop = math.Acos(clamp1(evec.X / ecc))
		if evec.Y < 0 {
			aop = 2*math.Pi - aop
		}
		ta = anomalyFromPeriapsis(evec, ecc, r, rn, rdotv)
	case circular: // circular inclined: ta is the argument of latitude from the node.
		raan = raanFromNode(n, nn)
		ta = math.Acos(clamp1(md3.Dot(n, r) / (nn * rn)))
		if r.Z < 0 {
			ta = 2*math.Pi - ta
		}
	default:
		raan = raanFromNode(n, nn)
		aop = math.Acos(clamp1(md3.Dot(n, evec) / (nn * ecc)))
		if evec.Z < 0 {
			aop = 2*math.Pi - aop
		}
		ta = anomalyFromPeriapsis(evec, ecc, r, rn, rdotv)
	}
	k, err = NewKeplerian(sma, ecc, inc, raan, aop)
	return k, ta, err
}

func raanFromNode(n md3.Vec, nn float64) float64 {
	raan := math.Acos(clamp1(n.X / nn))
	if n.Y < 0 {
		raan = 2*math.Pi - raan
	}
	return raan
}

func anomalyFromPeriapsis(evec md3.Vec, ecc float64, r md3.Vec, rn, rdotv float64) float64 {
	ta := math.Acos(clamp1(md3.Dot(evec, r) / (ecc * rn)))
	if rdotv < 0 {
		ta = 2*math.Pi - ta
	}
	return ta
}

// RV returns the inertial equatorial position [m] and velocity [m/s] at a
// true anomaly position, given the gravitational parameter mu [m³/s²].
// Vallado Alg. 10 (COE2RV): perifocal state rotated by R3(-raan)R1(-inc)R3(-aop).
func (k Keplerian) RV(mu, trueAnomaly float64) (r, v md3.Vec) {
	p := k.sma * (1 - k.ecc*k.ecc) // semi-latus rectum [m]
	sinta, costa := math.Sincos(trueAnomaly)
	rmag := p / (1 + k.ecc*costa)
	rpqw := md3.Vec{X: rmag * costa, Y: rmag * sinta}
	sqmup := math.Sqrt(mu / p)
	vpqw := md3.Vec{X: -sqmup * sinta, Y: sqmup * (k.ecc + costa)}

	sinO, cosO := math.Sincos(k.raan)
	sini, cosi := math.Sincos(k.inc)
	sinw, cosw := math.Sincos(k.aop)
	// Rows of the PQW→IJK rotation, Vallado Eqn (2-99).
	px := md3.Vec{X: cosO*cosw - sinO*sinw*cosi, Y: -cosO*sinw - sinO*cosw*cosi, Z: sinO * sini}
	py := md3.Vec{X: sinO*cosw + cosO*sinw*cosi, Y: -sinO*sinw + cosO*cosw*cosi, Z: -cosO * sini}
	pz := md3.Vec{X: sinw * sini, Y: cosw * sini, Z: cosi}
	r = md3.Vec{X: md3.Dot(px, rpqw), Y: md3.Dot(py, rpqw), Z: md3.Dot(pz, rpqw)}
	v = md3.Vec{X: md3.Dot(px, vpqw), Y: md3.Dot(py, vpqw), Z: md3.Dot(pz, vpqw)}
	return r, v
}

// SemiMajorAxis returns the semi-major axis [m].
func (k Keplerian) SemiMajorAxis() float64 { return k.sma }

// Eccentricity returns the orbit eccentricity [adim].
func (k Keplerian) Eccentricity() float64 { return k.ecc }

// Inclination returns the orbit inclination [rad] in range [0, pi].
func (k Keplerian) Inclination() float64 { return k.inc }

// RAAN returns the right ascension of the ascending node [rad] in range [0, 2pi).
func (k Keplerian) RAAN() float64 { return k.raan }

// ArgumentOfPeriapsis returns the argument of periapsis [rad] in range [0, 2pi).
func (k Keplerian) ArgumentOfPeriapsis() float64 { return k.aop }

// Elliptical returns the in-plane orbit shape, bridging to the
// anomaly/period/Kepler-equation machinery of Elliptical.
func (k Keplerian) Elliptical() (Elliptical, error) {
	return NewElliptical(k.sma*(1+k.ecc), k.sma*(1-k.ecc))
}

// Period returns the orbital period [s] given the gravitational parameter mu [m³/s²].
func (k Keplerian) Period(mu float64) float64 {
	return 2 * math.Pi * math.Sqrt(k.sma*k.sma*k.sma/mu) // Curtis Eqn (2.83)
}

// mod2pi normalizes an angle into [0, 2pi).
func mod2pi(rad float64) float64 {
	rad = math.Mod(rad, 2*math.Pi)
	if rad < 0 {
		rad += 2 * math.Pi
	}
	return rad
}

// clamp1 limits acos arguments against roundoff outside [-1, 1].
func clamp1(x float64) float64 { return math.Max(-1, math.Min(1, x)) }
