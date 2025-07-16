package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

type Coordinates interface {
	AGravG() md3.Vec
	TGE() md3.Mat3
	SetFromEarthFixedCoords(SBIE md3.Vec, epochTime float64)
	World() *World
}

var (
	_ Coordinates = (*GeocentricCoords)(nil)
	_ Coordinates = (*GeodesicCoords)(nil)
)

// Geocentric latitude, longitude, and elevation (height above earth reference sphere). See https://en.wikipedia.org/wiki/Geographic_coordinate_system
//
// See NewGeocentricFromDD for examples of usage.
type GeocentricCoords struct {
	// Longitude [rad]
	// The longitudinal lines on a map go from South to North. East is positive direction.
	Long float64
	// Latitude [rad]
	// The latitude lines on a map go from West to East. North is positive direction.
	Lat float64
	// Elev is the height above the reference sphere [m].
	//
	// Points at a distance of Rearth from center
	// of earth have elevation=0. HASL=Elev-seaLevelHeight
	Elev float64
	w    *World
}

// Geodesic as ellipsoidal coordinates using geodetic latitude, longitude, and elevation using WGS84 ellipsoidal model.
type GeodesicCoords struct {
	c GeocentricCoords
}

func (g GeocentricCoords) Geodesic() GeodesicCoords {
	return GeodesicCoords{c: g}
}

func (g GeocentricCoords) Degrees() (longitude float64, latitude float64) {
	return g.Long * 180 / math.Pi, g.Lat * 180 / math.Pi
}

// InertialCoords returns the planet-centered absolute inertial frame (ECI) of reference coordinates. See Earth-centered inertial.
func (g GeocentricCoords) InertialCoords(epochTime float64) (sBII md3.Vec, TGI md3.Mat3) {
	TEI := g.w.TEI(epochTime)
	TGE := g.TGE()
	TGI = md3.MulMat3(TGE, TEI)
	sBIE := g.EarthFixedCoords(epochTime)
	sBII = md3.MulMatVecTrans(TEI, sBIE)
	return sBII, TGI
}

// EarthFixedCoords returns the planet-centerd, planet-fixed (ECEF) frame of reference coordinates. These rotate with the planet. See Earth-centered, earth fixed.
func (g GeocentricCoords) EarthFixedCoords(epochTime float64) (sBIE md3.Vec) {
	slon, clon := math.Sincos(g.Long + g.w.Rotation*epochTime)
	slat, clat := math.Sincos(g.Lat)
	sBIE.X = clat * clon
	sBIE.Y = clat * slon
	sBIE.Z = slat
	radius := g.Elev + g.w.Radius
	return md3.Scale(radius, sBIE)
}

func (g *GeocentricCoords) SetFromEarthFixedCoords(sBIE md3.Vec, epochTime float64) {
	if g.w == nil {
		panic("nil world")
	}
	*g = g.w.GeocentricFromEarthFixedCoords(sBIE, epochTime)
}

func (g GeocentricCoords) World() *World { return g.w }

func (g GeocentricCoords) TGI(epochTime float64) md3.Mat3 {
	TEI := g.w.TEI(epochTime)
	TGE := g.TGE()
	TGI := md3.MulMat3(TGE, TEI)
	return TGI
}

func (g GeocentricCoords) TGE() md3.Mat3 {
	slo, clo := math.Sincos(g.Long)
	sla, cla := math.Sincos(g.Lat)
	return mat3(
		-sla*clo, -sla*slo, cla,
		-slo, clo, 0,
		-cla*clo, -cla*slo, -sla,
	)
}

func (g GeocentricCoords) Radius() float64 {
	return g.w.Radius + g.Elev
}

// AGravG returns gravity acceleration in geographic coordinates. [m.s^-2]
func (g GeocentricCoords) AGravG() (gravityVec md3.Vec) {
	dbi := g.Radius()
	gravityVec.Z = g.w.G() / (dbi * dbi)
	return gravityVec
}

func (g GeodesicCoords) AGravG() (gravityVec md3.Vec) {
	// Sqrt(0.5)
	const sqrtHalf = 0.7071067811865475244008443621048490392848359376884740365883398689
	const dum2 = 3 * sqrtHalf
	w := g.c.w
	dbi := g.c.Radius()
	dum1 := w.G() / (dbi * dbi)
	dum3 := w.SemiMajorAxis / dbi
	dum3 *= dum3 // square it, much faster than Pow
	sinlat, coslat := math.Sincos(g.c.Lat)
	gravityVec.X = -dum1 * dum2 * w.C20 * dum3 * sinlat * coslat
	gravityVec.Z = dum1 * (1 + dum2/2*w.C20*dum3*(3*sinlat*sinlat-1))
	return gravityVec
}

func (g *GeodesicCoords) SetFromEarthFixedCoords(sBIE md3.Vec, epochTime float64) {
	g.c.SetFromEarthFixedCoords(sBIE, epochTime)
}

func (g GeodesicCoords) World() *World { return g.c.w }

func (g GeodesicCoords) Geocentric() GeocentricCoords { return g.c }

func (g GeodesicCoords) TGE() md3.Mat3 { return g.c.TGE() }

// clampLongLat limits the value of rad to within range [-pi,pi] such that
//
//	sin(rad) == sin(clampLongLat(rad))
//	cos(rad) == cos(clampLongLat(rad))
func clampLongLat(rad float64) float64 {
	// TODO(soypat): maybe enclosing both checks in a single if statement for better branch prediction- write benchmarks before trying this: `if math.Abs(rad) > math.Pi ...`
	// TODO(soypat): Maybe better to replace this function with a more precise clamper? See Go's standard library math.satan(yes, that's the actual name) used in the Atan,Asin,Acos family .
	// we'd like to work in range [-pi,pi] for greatest precision of geometric math functions
	if rad < -math.Pi {
		rad += 2 * math.Pi
	} else if rad > math.Pi {
		rad -= 2 * math.Pi
	}
	return rad
}

// asinlong returns the longitude given x and y SBII coordinates
// without taking into account epoch time.
// TODO(pato) This function can be optimized.
func asinlong(y, x float64) float64 {
	long := math.Asin(y / math.Hypot(x, y))
	switch {
	// case x >= 0 && y >= 0:
	// Quadrant I.
	// Do nothing.
	case x < 0 && y >= 0:
		// Quadrant II.
		long = math.Pi - long

	case x < 0 && y < 0: // TODO merge this with above case clause to optimize after writing tests.
		// Quadrant III.
		long = math.Pi - long

	case x >= 0 && y < 0:
		// Quadrant IV.
		long = 2*math.Pi + long
	}
	return long
}

func mat3(a, b, c, d, e, f, g, h, i float64) md3.Mat3 {
	return md3.NewMat3([]float64{
		a, b, c,
		d, e, f,
		g, h, i,
	})
}
