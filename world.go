package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// GeographicVectorFromElevationAndBearing returns a vector in geographic coordinates
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
