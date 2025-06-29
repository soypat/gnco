package gnco

import "github.com/soypat/geometry/md3"

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
