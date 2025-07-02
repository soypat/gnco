package orbits

import (
	"fmt"
	"math"

	"github.com/soypat/geometry/md1"
)

// Orbit defines a typical earthbound circular or elliptic orbit at
// an inclination plane. Most, if not all logic, implemented with
// Curtis, Howard's Orbital Mechanics for Mechanical Engineering Students - Third edition.
type Elliptical struct {
	// ra and rp [m] are the apoapsis and periapsis radius of an elliptical orbit.
	// Setting periapsis radius rp to 0 implies a circular orbit.
	ra, rp float64
	// plane [rad] defines the orbital inclination w.r.t the equator.
	// A plane of 0 is an equatorial orbit. A plane of pi or 90 degrees is a polar orbit.
	// plane float64
}

func NewCircular(r float64) (Elliptical, error) {
	return NewElliptical(r, r)
}

func NewElliptical(ra, rp float64) (Elliptical, error) {
	if ra < rp || ra < 0 {
		return Elliptical{}, fmt.Errorf("got bad argument to NewOrbit: ra=%.5gkm, rp=%.5gkm", ra/1e3, rp/1e3)
	}
	if ra == rp {
		rp = 0 // Circular orbit case, ra is non-zero in backend.
	}
	o := Elliptical{ra: ra, rp: rp}
	e := o.Eccentricity()
	if e >= 1 || e < 0 {
		return Elliptical{}, fmt.Errorf("bad orbit eccentricity, got %.3g. Orbit does not expect parabolic orbits", e)
	}
	return o, nil
}

// Apoapsis returns the maximum distance from center of earth o reaches. [m]
// To calculate maximum orbit height:
//
//	maxHeight := o.Apoapsis() - world.Radius
func (o Elliptical) Apoapsis() float64 { return o.ra }

// Apoapsis returns the minimum distance from earth o reaches. [m]
// To calculate minimum orbit height:
//
//	minHeight := o.Periapsis() - world.Radius
func (o Elliptical) Periapsis() float64 {
	if o.rp == 0 {
		return o.ra // Is circular orbit.
	}
	return o.rp
}

// IsCircular returns true if the orbit is circular to within a given tolerance.
func (o Elliptical) IsCircular(tol float64) bool {
	return o.rp == 0 || o.ra-o.rp < tol
}

// Eccentricity returns parameter of eccentricity of orbit. More elliptical orbits have higher eccentricity.
// A circular orbit has an eccentricity of 1. [adim]
func (o Elliptical) Eccentricity() float64 {
	a, c := o.a(), o.c() // On a separate line for debugging purposes.
	return c / a
}

// Ellipse semimajor axis length [m]. a=(Ra + Rp)/2
func (o Elliptical) a() float64 { return 0.5 * (o.Apoapsis() + o.Periapsis()) }

// Ellipse semiminor axis length [m]. b=a*sqrt(1-e^2)
func (o Elliptical) b() float64 {
	e := o.c() / o.a() // Direct eccentricity calculation.
	return o.a() * math.Sqrt(1-e*e)
}

// Ellipse distance between focus and center [m].
// WARNING: Some bibliogaphies define this as distance between focii, so double this c.
func (o Elliptical) c() float64 { return o.a() - o.Periapsis() }

// Cartesian returns x and y coordinate [m] for the orbit at a
// true anomaly position with the coordinates centered on the ellipse center
// and the x axis aligned with the semimajor pointing towards earth (periapsis).
func (o Elliptical) CartesianCoordinates(trueAnomaly float64) (x, y float64) {
	e := o.Eccentricity()
	sint, cost := math.Sincos(trueAnomaly)
	x = o.a() * (e + cost) / (1 + e*cost)              // Eqn (2.77)
	y = o.b() * sint * math.Sqrt(1-e*e) / (1 + e*cost) // Eqn (2.78)
	return x, y
}

// FlightPathAngle is defined as the angle the velocity vector makes with
// the normal to the position vector. The normal to position vector points
// in the direction of the tangential velocity.
//
//	tan(gamma) = v_radial / v_tangential  // where gamma is the flight path angle.
func (o Elliptical) FlightPathAngle(trueAnomaly float64) float64 {
	e := o.Eccentricity()
	sint, cost := math.Sincos(trueAnomaly)
	return math.Atan2(e*sint, 1+e*cost) // Eqn. (2.52)
}

// MeanAnomaly returns the fraction of an elliptical orbit's period that has elapsed since the orbiting body passed periapsis
// as an angle [rad]. For reference: the mean anomaly increases uniformly during an elliptic orbit. The true anomaly does not.
// For circular orbits the mean anomaly matches the true anomaly.
// MeanAnomaly shall return a value in range [0, 2pi).
func (o Elliptical) MeanAnomaly(trueAnomaly float64) float64 {
	if o.IsCircular(0) {
		return trueAnomaly
	}
	e := o.Eccentricity()
	sint, cost := math.Sincos(trueAnomaly)
	sqrt1esin := math.Sqrt(1-e*e) * sint
	return math.Atan2(-sqrt1esin, -e-cost) + math.Pi - e*sqrt1esin/(1+e*cost) // Eqn (3.3)
}

// EccentricAnomaly returns the eccentric anomaly angular parameter that defines an orbit.
// Usually stylized as upper case E in literature.
func (o Elliptical) EccentricAnomaly(trueAnomaly float64) float64 {
	if o.IsCircular(0) {
		return trueAnomaly
	}
	const iter = 30 // TODO(soypat): How many iterations is enough? Depends on eccentricity I think.
	e := o.Eccentricity()
	Me := o.MeanAnomaly(trueAnomaly)
	sum := 0.0
	n := 1.0
	for i := 1; i < iter; i++ {
		// Solve with fourier series.
		// TODO(soypat): Use modified newton raphson https://academic.oup.com/mnras/article/467/2/1702/2929272 An efficient code to solve the Kepler equation. Elliptic case  V. Raposo-Pulido, J. Peláez
		sum += math.Jn(i, n*e) / n * math.Sin(n*Me)
		n++
	}
	return Me + 2*sum
}

/*
 * Elliptical Physics.
 */

// AngularMomentum returns the angular momentum of the orbit given the gravitational parameter of the world [m/s^2].
func (o Elliptical) AngularMomentum(gravParam float64) float64 {
	// We evaluate the orbit equation at perigee where trueAnomaly==0
	// and solve for h.l
	return math.Sqrt(gravParam * o.rp * (1 + o.Eccentricity()))
}

func (o Elliptical) SpecificEnergy(gravParam float64) float64 {
	return -gravParam / (2 * o.a()) // Eqn (2.80). See also Eqn (2.60)
}

// Period returns the period of the orbit, or the amount of time it takes to complete a single orbit around world. [s]
func (o Elliptical) Period(gravParam float64) float64 {
	if o.IsCircular(0) {
		return 2 * math.Pi * o.ra / math.Hypot(o.Velocity(gravParam, 0))
	}
	// Not validated!
	a := o.a()
	return 2 * math.Pi * math.Sqrt(a*a*a/gravParam) // Eqn (2.83)
}

// ElapsedSincePeriapsis returns the seconds elapsed since periapsis.
// This function will return a value in range [0, T] where T is the orbit period.
func (o Elliptical) ElapsedSincePeriapsis(gravParam, trueAnomaly float64) float64 {
	M := o.MeanAnomaly(trueAnomaly)
	T := o.Period(gravParam)
	return M * T / (2 * math.Pi)
}

func (o Elliptical) TrueAnomalyFromElapsedSincePeriapsis(gravParam, elapsedSincePeriapsis, tol float64) float64 {
	T := o.Period(gravParam)
	Me := 2 * math.Pi * elapsedSincePeriapsis / T
	e := o.Eccentricity()
	solver := md1.DefaultNewtonRaphsonSolver()
	solver.Tolerance = tol
	E, convergedIn := solver.Root(Me-e/2, func(xGuess float64) float64 {
		return xGuess - e*math.Sin(xGuess) - Me
	})
	if convergedIn < 0 {
		return math.NaN()
	}
	rhs := math.Sqrt((1+e)/(1-e)) * math.Tan(E/2)
	trueAnomaly := 2 * math.Atan(rhs) // En (3.10a)
	return math.Abs(trueAnomaly)
}

// DistanceToCenter solves the orbit equation as given by Curtis, Howard in
// Orbital Mechanics for Mechanical Engineering Students, Eqn. (2.71).
func (o Elliptical) DistanceToCenter(gravParam, trueAnomaly float64) float64 {
	h := o.AngularMomentum(gravParam)
	return h * h / (gravParam * (1 + o.Eccentricity()*math.Cos(trueAnomaly)))
}

// Velocity returns the velocity components of the orbit given a trueAnomaly position. [m/s]
func (o Elliptical) Velocity(gravParam, trueAnomaly float64) (vRadial, vTangential float64) {
	sint, cost := math.Sincos(trueAnomaly)
	gdivh := gravParam / o.AngularMomentum(gravParam)
	e := o.Eccentricity()
	vRadial = gdivh * e * sint         // Eqn (2.49)
	vTangential = gdivh * (1 + e*cost) // Eqn (2.48)
	return vRadial, vTangential
}
