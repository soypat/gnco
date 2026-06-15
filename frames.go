package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

// Frame identifies a coordinate frame used for orientation conversions.
type Frame rune

const (
	FrameInertial   Frame = 'I' // inertial
	FrameGeographic Frame = 'G' // geographic
	FrameVelocity   Frame = 'V' // velocity
	FrameBody       Frame = 'B' // body
)

// Orientation holds the DCMs that chain body → velocity → geographic → inertial.
// Each matrix T_XY satisfies T_XY * v_Y = v_X (subscript convention: destination
// first, source second), so MulMatVec converts source→destination and MulMatVecTrans
// converts destination→source. Full inertial-to-body chain: TBV * TVG * TGI.
type Orientation struct {
	TBV md3.Mat3 // velocity → body (v_B = TBV * v_V).
	TVG md3.Mat3 // geographic → velocity (v_V = TVG * v_G).
	TGI md3.Mat3 // inertial → geographic (v_G = TGI * v_I).
}

// ToInertial converts frameVec from frame F into inertial frame of reference.
func (F Frame) ToInertial(dir Orientation, frameVec md3.Vec) md3.Vec {
	switch F {
	case FrameBody:
		frameVec = md3.MulMatVecTrans(dir.TBV, frameVec)
		fallthrough
	case FrameVelocity:
		frameVec = md3.MulMatVecTrans(dir.TVG, frameVec)
		fallthrough
	case FrameGeographic:
		frameVec = md3.MulMatVecTrans(dir.TGI, frameVec)
	case FrameInertial:
		// no conversion needed
	default:
		panic("unknown frame")
	}
	return frameVec
}

// ToGeographic converts frameVec from frame F into geographic frame of reference.
func (F Frame) ToGeographic(v Orientation, frameVec md3.Vec) md3.Vec {
	switch F {
	case FrameBody:
		frameVec = md3.MulMatVecTrans(v.TBV, frameVec)
		fallthrough
	case FrameVelocity:
		frameVec = md3.MulMatVecTrans(v.TVG, frameVec)
	case FrameGeographic:
		// no conversion needed
	case FrameInertial:
		frameVec = md3.MulMatVec(v.TGI, frameVec)
	default:
		panic("unknown frame")
	}
	return frameVec
}

// ToVelocity converts frameVec from frame F into velocity frame of reference.
func (F Frame) ToVelocity(v Orientation, frameVec md3.Vec) md3.Vec {
	switch F {
	case FrameBody:
		frameVec = md3.MulMatVecTrans(v.TBV, frameVec)
	case FrameVelocity:
		// no conversion needed
	case FrameInertial:
		frameVec = md3.MulMatVec(v.TGI, frameVec)
		fallthrough
	case FrameGeographic:
		frameVec = md3.MulMatVec(v.TVG, frameVec)
	default:
		panic("unknown frame")
	}
	return frameVec
}

// ToBody converts frameVec from frame F to body frame of reference.
func (F Frame) ToBody(v Orientation, frameVec md3.Vec) md3.Vec {
	switch F {
	case FrameInertial:
		frameVec = md3.MulMatVec(v.TGI, frameVec)
		fallthrough
	case FrameGeographic:
		frameVec = md3.MulMatVec(v.TVG, frameVec)
		fallthrough
	case FrameVelocity:
		frameVec = md3.MulMatVec(v.TBV, frameVec)
	case FrameBody:
		// No conversion needed.
	default:
		panic("unknown frame")
	}
	return frameVec
}

// TVGFromGeographicVelocity returns the geographic-to-velocity transform matrix.
// vbg is the vehicle velocity in geographic coordinates.
//
//	TVG * v_G = v_V; transpose(TVG) converts V to G.
//
// The V-frame x-axis aligns with velocity. Near-vertical motion uses north as azimuth.
func TVGFromGeographicVelocity(vbg md3.Vec) (TVG md3.Mat3) {
	vnorm := md3.Norm(vbg)
	if vnorm == 0 {
		return md3.IdentityMat3()
	}
	vx := vbg.X / vnorm
	vy := vbg.Y / vnorm
	vz := vbg.Z / vnorm
	hspeed := math.Hypot(vx, vy)
	if hspeed < 1e-9 {
		// Near-vertical: azimuth undefined; default to North (ψ=0).
		// V_y = East = [0,1,0]; V_z = [-vz,0,0] (South when going up, North when going down).
		return mat3(
			0, 0, vz,
			0, 1, 0,
			-vz, 0, 0,
		)
	}
	// Rows are the V-frame basis vectors expressed in geographic (NED) coordinates:
	//   Row 0 (V_x): velocity direction
	//   Row 1 (V_y): right wing — normalize(ẑ_down × V_x), East when flying North
	//   Row 2 (V_z): down — normalize(V_x × V_y)
	return mat3(
		vx, vy, vz,
		-vy/hspeed, vx/hspeed, 0,
		-vz*vx/hspeed, -vz*vy/hspeed, hspeed,
	)
}

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
