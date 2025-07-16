package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// NewEarth returns the world as we know it according to WGS84.
func NewEarth() *World {
	return &World{
		Mass:           5.973332e24,
		C20:            -4.8416685e-4,
		SemiMajorAxis:  6378137, // WGS84 [m]
		Rotation:       7.292114999999999893e-05,
		Radius:         6370987.,
		seaLevelRadius: 6371146,
		flattening:     3.33528106e-3,
		celestialLong:  0,

		// SGP4 according to WGS84. Recommended by IAU to propagate orbits
		Ke: 0.07436685316871385,
		J2: 0.00108262998905,
		J3: -0.00000253215306,
		J4: -0.00000161098761,
	}
}

// Important constants.
const (
	// universal gravitational constant - [N.m^2.kg^-2]
	bigG = 6.673e-11
)

type World struct {
	Mass           float64 // total mass of planet. [kg]
	C20            float64 // second degree zonal gravitational coefficient [Adim]
	SemiMajorAxis  float64 // semi-major axis of planet [m], also known as "equatorial radius"
	Rotation       float64 // Angular rotation of planet w.r.t inertial frame [rad/s]
	Radius         float64 // Radius of planet [m]
	seaLevelRadius float64 // If earth stopped rotating the sea level would take this distance from center of earth [m] https://www.esri.com/news/arcuser/0703/geoid3of3.html
	flattening     float64 // Flattening of planet, (WGS84) [Adim]
	celestialLong  float64 // Celestial longitude, for earth is Greenwich meridian. Will indicate start of epoch [rad]

	// SGP4 parameters:

	// kₑ Square root of earth's gravitational parameter in earth [radii³ min⁻²].
	Ke float64
	// J₂ Un-normalised second zonal harmonic.
	J2 float64
	// J₃ Un-normalised third zonal harmonic.
	J3 float64
	// J₄ Un-normalised fourth zonal harmonic.
	J4 float64
}

func (w *World) GeocentricFromEarthFixedCoords(sBIE md3.Vec, epochTime float64) GeocentricCoords {
	dbi := md3.Norm(sBIE)
	lat := math.Asin(sBIE.Z / dbi)
	elev := dbi - w.Radius
	// longitude calculation using specialized quadrant algorithm and total earth rotation.
	long := asinlong(sBIE.Y, sBIE.X) - w.Rotation*epochTime + w.celestialLong
	long = clampLongLat(long)
	return GeocentricCoords{
		w:    w,
		Long: long,
		Lat:  lat,
		Elev: elev,
	}
}

func (w *World) GeocentricFromDegrees(longDeg, latDeg, elevationAboveRefSphere float64) (longlat GeocentricCoords) {
	if elevationAboveRefSphere < -w.Radius {
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
	return w.Mass * bigG
}

// seaLevelHeight is the height of sea level above earth reference sphere.
func (w *World) seaLevelHeight() float64 {
	return w.seaLevelRadius - w.Radius
}

func (w *World) HASLToElevation(hasl float64) float64 {
	if w.seaLevelRadius == 0 {
		return hasl
	}
	return w.seaLevelHeight() + hasl
}

// TEI returns the [T]^{EI} transformation tensor given the epochTime in seconds.
func (w *World) TEI(epochTime float64) md3.Mat3 {
	daysPerSec := 1. / w.Day()
	sin, cos := math.Sincos(epochTime * daysPerSec)
	return mat3(
		cos, sin, 0,
		-sin, cos, 0,
		0, 0, 1,
	)
}

// Day returns amount of seconds in a day.
func (w *World) Day() float64 {
	return 2 * math.Pi / w.Rotation
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
