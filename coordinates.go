package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/cosmos"
)

// Coordinates represents location-related operations in geographic and
// body-fixed reference frames.
type Coordinates interface {
	// AGravG returns gravity acceleration in geographic coordinates.
	AGravG() md3.Vec
	// TGE returns the transform from geographic to body-fixed frame.
	TGE() md3.Mat3
	// SetFromEarthFixedCoords updates coordinates from a planet-centered
	// inertial (ECI) vector at the given epoch.
	SetFromEarthFixedCoords(SBI md3.Vec, epoch cosmos.Epoch)
	// Body returns the associated celestial body model.
	Body() *cosmos.Body
}

// Compile-time guarantee of interface implementation.
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
	w    *cosmos.Body
}

// Geodesic as ellipsoidal coordinates using geodetic latitude, longitude, and elevation using WGS84 ellipsoidal model.
type GeodesicCoords struct {
	c GeocentricCoords
}

// Geodesic returns the equivalent [GeodesicCoords] representation.
func (g GeocentricCoords) Geodesic() GeodesicCoords {
	return GeodesicCoords{c: g}
}

// Degrees returns the longitude and latitude in degrees.
func (g GeocentricCoords) Degrees() (longitude, latitude float64) {
	return g.Long * 180 / math.Pi, g.Lat * 180 / math.Pi
}

// InertialCoords returns the planet-centered absolute inertial frame (ECI) of reference coordinates. See Earth-centered inertial.
func (g GeocentricCoords) InertialCoords(epochTime cosmos.Epoch) (sBII md3.Vec, TGI md3.Mat3) {
	TEI := g.w.TEI(epochTime)
	TGE := g.TGE()
	TGI = md3.MulMat3(TGE, TEI)
	sBIE := g.EarthFixedCoords()
	// sBIE is body-fixed (ECEF); TEI maps inertial→body-fixed, so its transpose
	// rotates the body-fixed vector back into the inertial frame.
	sBII = md3.MulMatVecTrans(TEI, sBIE)
	return sBII, TGI
}

// EarthFixedCoords returns the planet-centered, planet-fixed (ECEF) frame of reference coordinates. These rotate with the planet. See Earth-centered, earth fixed.
//
// The body-fixed position depends only on the stored Long/Lat/Elev and is therefore
// time-independent; the epoch dependence lives entirely in [GeocentricCoords.InertialCoords].
func (g GeocentricCoords) EarthFixedCoords() (sBIE md3.Vec) {
	slon, clon := math.Sincos(g.Long)
	slat, clat := math.Sincos(g.Lat)
	sBIE.X = clat * clon
	sBIE.Y = clat * slon
	sBIE.Z = slat
	radius := g.Elev + g.w.Radius()
	return md3.Scale(radius, sBIE)
}

// SetFromEarthFixedCoords updates coordinates from a planet-centered inertial (ECI)
// vector at the given epoch. Implements [Coordinates].
func (g *GeocentricCoords) SetFromEarthFixedCoords(sBII md3.Vec, epochTime cosmos.Epoch) {
	if g.w == nil {
		panic("nil body")
	}
	*g = NewGeocentricFromEarthFixed(g.w, sBII, epochTime)
}

// Body returns the reference [cosmos.Body] model used by these coordinates.
func (g GeocentricCoords) Body() *cosmos.Body { return g.w }

// TGI returns the transform matrix from geographic coordinates to the
// planet-centered inertial frame at the given epoch time.
//
//	TGI = TGE*TEI
func (g GeocentricCoords) TGI(epochTime cosmos.Epoch) md3.Mat3 {
	TEI := g.w.TEI(epochTime)
	TGE := g.TGE()
	TGI := md3.MulMat3(TGE, TEI)
	return TGI
}

// TGE returns the transform matrix from geographic coordinates to the
// planet-fixed Earth-centered frame for the current geocentric coordinates.
func (g GeocentricCoords) TGE() md3.Mat3 {
	slo, clo := math.Sincos(g.Long)
	sla, cla := math.Sincos(g.Lat)
	return mat3(
		-sla*clo, -sla*slo, cla,
		-slo, clo, 0,
		-cla*clo, -cla*slo, -sla,
	)
}

// Radius returns the distance from body center to geocentric coordinate [m].
func (g GeocentricCoords) Radius() float64 {
	return g.w.Radius() + g.Elev
}

// HASL returns height above sea level [m].
func (g GeocentricCoords) HASL() float64 {
	return g.w.HASL(g.Elev)
}

// AGravG returns gravity acceleration in geographic coordinates. [m.s^-2]
func (g GeocentricCoords) AGravG() (gravityVec md3.Vec) {
	dbi := g.Radius()
	gravityVec.Z = g.w.Mu() / (dbi * dbi)
	return gravityVec
}

// AGravG returns gravity acceleration in geographic coordinates for a geodesic point.
// The vector is expressed in the geographic frame and accounts for ellipsoidal gravity
// variations using the current world model.
func (g GeodesicCoords) AGravG() (gravityVec md3.Vec) {
	// Sqrt(0.5)
	const sqrtHalf = 0.7071067811865475244008443621048490392848359376884740365883398689
	const dum2 = 3 * sqrtHalf
	w := g.c.w
	dbi := g.c.Radius()
	dum1 := w.Mu() / (dbi * dbi)
	dum3 := w.SemiMajorAxis() / dbi
	dum3 *= dum3 // square it, much faster than Pow
	c20 := w.C20()
	sinlat, coslat := math.Sincos(g.c.Lat)
	gravityVec.X = -dum1 * dum2 * c20 * dum3 * sinlat * coslat
	gravityVec.Z = dum1 * (1 + dum2/2*c20*dum3*(3*sinlat*sinlat-1))
	return gravityVec
}

// SetFromEarthFixedCoords updates coordinates from a planet-centered inertial (ECI)
// vector at the given epoch. Implements [Coordinates].
func (g *GeodesicCoords) SetFromEarthFixedCoords(sBII md3.Vec, epochTime cosmos.Epoch) {
	g.c.SetFromEarthFixedCoords(sBII, epochTime)
}

// Body returns the reference [cosmos.Body] model used by these coordinates.
func (g GeodesicCoords) Body() *cosmos.Body { return g.c.w }

// HASL returns height above sea level [m].
func (g GeodesicCoords) HASL() float64 { return g.c.HASL() }

// Geocentric returns the equivalent [GeocentricCoords] representation.
func (g GeodesicCoords) Geocentric() GeocentricCoords { return g.c }

// TGE returns the geographic-to-earth-fixed transform matrix.
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

// NewGeocentricFromDegrees builds geocentric coordinates for body b from a
// body-fixed longitude and latitude in degrees and an elevation above the body
// reference sphere [m].
func NewGeocentricFromDegrees(b *cosmos.Body, longDeg, latDeg, elevationAboveRefSphere float64) GeocentricCoords {
	if elevationAboveRefSphere < -b.Radius() {
		panic("bad elevation")
	}
	return GeocentricCoords{
		Long: clampLongLat(math.Pi / 180. * longDeg),
		Lat:  clampLongLat(math.Pi / 180. * latDeg),
		Elev: elevationAboveRefSphere,
		w:    b,
	}
}

// NewGeocentricFromEarthFixed builds geocentric coordinates for body b from a
// planet-centered inertial (ECI) vector at epoch e. The inertial vector is
// rotated into the body-fixed frame via b.TEI before extracting longitude,
// latitude and elevation, so the conversion is consistent with the full
// IAU-76/FK5 reduction used everywhere else.
func NewGeocentricFromEarthFixed(b *cosmos.Body, sBII md3.Vec, e cosmos.Epoch) GeocentricCoords {
	sBF := md3.MulMatVec(b.TEI(e), sBII)
	dbf := md3.Norm(sBF)
	return GeocentricCoords{
		Long: clampLongLat(math.Atan2(sBF.Y, sBF.X)),
		Lat:  math.Asin(sBF.Z / dbf),
		Elev: dbf - b.Radius(),
		w:    b,
	}
}

func mat3(a, b, c, d, e, f, g, h, i float64) md3.Mat3 {
	return md3.NewMat3([]float64{
		a, b, c,
		d, e, f,
		g, h, i,
	})
}
