package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// NewEarth returns the world as we know it according to WGS84.
func NewEarth() *World {
	return &World{
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

// Important constants.
const (
	// universal gravitational constant - [N.m^2.kg^-2]
	bigG = 6.673e-11
)

type World struct {
	mass           float64 // total mass of planet. [kg]
	c20            float64 // second degree zonal gravitational coefficient [Adim]
	semiMajorAxis  float64 // semi-major axis of planet [m], also known as "equatorial radius"
	rotation       float64 // Angular rotation of planet w.r.t inertial frame [rad/s]
	radius         float64 // Radius of planet [m]
	seaLevelRadius float64 // If earth stopped rotating the sea level would take this distance from center of earth [m] https://www.esri.com/news/arcuser/0703/geoid3of3.html
	flattening     float64 // Flattening of planet, (WGS84) [Adim]
	celestialLong  float64 // Celestial longitude, for earth is Greenwich meridian. Will indicate start of epoch [rad]

	// SGP4 parameters:
	ke float64 // kₑ: square root of gravitational parameter in Earth radii³ min⁻²
	j2 float64 // J₂: un-normalised second zonal harmonic
	j3 float64 // J₃: un-normalised third zonal harmonic
	j4 float64 // J₄: un-normalised fourth zonal harmonic
}

// Mass returns the total mass of the planet [kg].
func (w *World) Mass() float64 { return w.mass }

// C20 returns the second degree zonal gravitational coefficient (dimensionless).
func (w *World) C20() float64 { return w.c20 }

// SemiMajorAxis returns the equatorial radius [m] (WGS84 semi-major axis).
func (w *World) SemiMajorAxis() float64 { return w.semiMajorAxis }

// Rotation returns the angular rotation rate of the planet w.r.t the inertial frame [rad/s].
func (w *World) Rotation() float64 { return w.rotation }

// Radius returns the mean radius of the planet [m].
func (w *World) Radius() float64 { return w.radius }

// Ke returns the SGP4 square root of World's gravitational parameter [world radii³ min⁻²].
func (w *World) Ke() float64 { return w.ke }

// J2 returns the SGP4 un-normalised second zonal harmonic.
func (w *World) J2() float64 { return w.j2 }

// J3 returns the SGP4 un-normalised third zonal harmonic.
func (w *World) J3() float64 { return w.j3 }

// J4 returns the SGP4 un-normalised fourth zonal harmonic.
func (w *World) J4() float64 { return w.j4 }

func (w *World) GeocentricFromEarthFixedCoords(sBIE md3.Vec, epochTime float64) GeocentricCoords {
	dbi := md3.Norm(sBIE)
	lat := math.Asin(sBIE.Z / dbi)
	elev := dbi - w.radius
	// longitude calculation using specialized quadrant algorithm and total earth rotation.
	long := asinlong(sBIE.Y, sBIE.X) - w.rotation*epochTime + w.celestialLong
	long = clampLongLat(long)
	return GeocentricCoords{
		w:    w,
		Long: long,
		Lat:  lat,
		Elev: elev,
	}
}

func (w *World) GeocentricFromDegrees(longDeg, latDeg, elevationAboveRefSphere float64) (longlat GeocentricCoords) {
	if elevationAboveRefSphere < -w.radius {
		panic("bad elevatiojn")
	}
	return GeocentricCoords{
		Long: clampLongLat(math.Pi / 180. * longDeg),
		Lat:  clampLongLat(math.Pi / 180. * latDeg),
		Elev: elevationAboveRefSphere,
		w:    w,
	}
}

// G Gravitational parameter=G*mass of world [m^3.s^-2]
func (w *World) G() float64 {
	return w.mass * bigG
}

// seaLevelHeight is the height of sea level above earth reference sphere.
func (w *World) seaLevelHeight() float64 {
	return w.seaLevelRadius - w.radius
}

func (w *World) HASLToElevation(hasl float64) float64 {
	if w.seaLevelRadius == 0 {
		return hasl
	}
	return w.seaLevelHeight() + hasl
}

// HASL returns height above sea level [m] for an elevation above the reference sphere.
// It is the inverse of HASLToElevation.
func (w *World) HASL(elev float64) float64 {
	if w.seaLevelRadius == 0 {
		return elev
	}
	return elev - w.seaLevelHeight()
}

// TEI returns the [T]^{EI} transformation tensor given the epochTime in seconds.
func (w *World) TEI(epochTime float64) md3.Mat3 {
	// argument to sincos: [s]*[rad/s]=[rad]
	sin, cos := math.Sincos(epochTime * w.rotation)
	return mat3(
		cos, sin, 0,
		-sin, cos, 0,
		0, 0, 1,
	)
}

// Day returns amount of seconds in a day.
func (w *World) Day() float64 {
	return 2 * math.Pi / w.rotation
}

// GeographicFromElevBearing returns a vector in geographic coordinates
// pointing in direction given by an elevation and bearing in radians.
// When obtaining Geographic coordinates bearing can be thought of as North/West/East/South
// parameter, while the elevation describes whether direction is Up or Down, with
//
// Geographic coordinates:
//
//	X: North
//	Y: East // <- TODO this looks wrong...
//	Z: Center of earth
//
// Elevation:
//
//	Pi/2: Pointing up.
//	0: Pointing towards horizon.
//	-Pi/2: Pointing down.
//
// Bearing:
//
//	0: Pointing North.
//	Pi/2: Pointing East.
//	Pi: Pointing South.
func GeographicVectorFromElevationAndBearing(elevation, bearing, NormOfVector float64) (dirG md3.Vec) {
	// See CADAC matcar routine.
	sine, cose := math.Sincos(elevation)
	sinb, cosb := math.Sincos(bearing)
	dirG = md3.Vec{
		X: cosb * cose,
		Y: -sinb * cose,
		Z: -sine,
	}
	dirG = md3.Scale(NormOfVector, md3.Unit(dirG)) // TODO: does this need to be normalized before scaling?
	return dirG
}
