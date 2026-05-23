package gnco

import (
	"math"

	"github.com/soypat/geometry/md3"
)

type Frame rune

const (
	FrameInertial   Frame = 'I'
	FrameGeographic Frame = 'G'
	FrameVelocity   Frame = 'V'
	FrameBody       Frame = 'B'
)

type Orientation struct {
	TBV md3.Mat3 // Rotation tensor: body to velocity coordinates.
	TVG md3.Mat3 // Rotation tensor: velocity to geographical coordinates.
	TGI md3.Mat3 // Rotation tensor: geographical to inertial coordinates.
}

// ToInertial converts frameVec in the given F frame to inertial frame of reference.
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

// ToGeographic converts frameVec in the given F frame to geographic frame of reference.
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

// ToVelocity converts frameVec in the given F frame to velocity frame of reference.
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

// TVGFromGeographicVelocity returns the transformation matrix from geographic (G)
// to velocity (V) frame given vbg, the vehicle velocity in geographic coordinates.
// TVG satisfies TVG * v_G = v_V; its transpose converts the opposite direction (V→G).
//
// The V-frame X-axis is aligned with the velocity direction. For near-vertical flight
// (horizontal component < 1e-9 of total speed) azimuth defaults to zero (North).
func TVGFromGeographicVelocity(vbg md3.Vec) md3.Mat3 {
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

// ToBody converts frameVec in the given F frame to body frame of reference.
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
