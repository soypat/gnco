package physics

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/physics/ode"
)

// RigidBodyIntegrator integrates a rigid body's full 6-DOF state —
// inertial position and velocity plus attitude (quaternion) and body-frame
// angular velocity — under gravity and externally supplied body-frame force and
// torque.
type RigidBodyIntegrator struct {
	// Unlike the point-mass integrator, attitude dynamics are velocity-dependent
	// (the gyroscopic term ω×Iω and the quaternion kinematics q̇=½q⊗ω), so the
	// velocity-independent RKN1210 cannot be used. The whole state is integrated as
	// a single first-order system with RKF78.
	//
	// The 13-element state vector is laid out as:
	//
	//	[0:3]   r   inertial position
	//	[3:6]   v   inertial velocity
	//	[6:10]  q   attitude quaternion body→inertial, stored {I,J,K,W}
	//	[10:13] ω   angular velocity, body frame
	integrator     ode.RKF78
	coord          gnco.Coordinates
	mass           float64
	inertia        md3.Mat3
	inertiaInv     md3.Mat3
	lastForceBody  md3.Vec // external body-frame force, held constant across a Step
	lastTorqueBody md3.Vec // external body-frame torque, held constant across a Step
}

// Configure initializes the rigid-body integrator.
// SBI0/VBI0 are the initial inertial position and velocity, att0 is the initial
// body→inertial attitude (normalized internally), omega0 is the initial
// body-frame angular velocity, mass is the body mass [kg] and inertia is the
// body-frame inertia tensor [kg·m²].
func (rbi *RigidBodyIntegrator) Configure(coord gnco.Coordinates, t0 float64, SBI0, VBI0 md3.Vec, att0 md3.Quat, omega0 md3.Vec, mass float64, inertia md3.Mat3) error {
	*rbi = RigidBodyIntegrator{
		coord:      coord,
		mass:       mass,
		inertia:    inertia,
		inertiaInv: inertia.Inverse(),
	}
	err := rbi.integrator.Configure(ode.Parameters{
		AbsTolerance: 0,
		MinStep:      0,
		MaxStep:      0,
	})
	if err != nil {
		return err
	}
	q := att0.Unit()
	rbi.integrator.Init(ode.IVP1{
		T0: t0,
		Y0: []float64{
			SBI0.X, SBI0.Y, SBI0.Z,
			VBI0.X, VBI0.Y, VBI0.Z,
			q.I, q.J, q.K, q.W,
			omega0.X, omega0.Y, omega0.Z,
		},
		Func: rbi.rates,
	})
	return nil
}

// State returns the current time, inertial position and velocity, body→inertial
// attitude and body-frame angular velocity.
func (rbi *RigidBodyIntegrator) State() (t float64, SBI, VBI md3.Vec, att md3.Quat, omega md3.Vec) {
	t, y := rbi.integrator.State()
	return t,
		md3.Vec{X: y[0], Y: y[1], Z: y[2]},
		md3.Vec{X: y[3], Y: y[4], Z: y[5]},
		md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]},
		md3.Vec{X: y[10], Y: y[11], Z: y[12]}
}

// Step advances the state by dt using external body-frame force and torque,
// which are held constant across the step. Gravity is computed internally per
// stage. The attitude quaternion is renormalized after the step.
func (rbi *RigidBodyIntegrator) Step(dt float64, forceBody, torqueBody md3.Vec) (t float64, SBI, VBI md3.Vec, att md3.Quat, omega md3.Vec) {
	rbi.lastForceBody = forceBody
	rbi.lastTorqueBody = torqueBody
	rbi.integrator.Step(dt)
	// Post-step projection: renormalize the attitude quaternion to undo the
	// small norm drift accumulated over the step.
	tNow, y := rbi.integrator.State()
	q := md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]}.Unit()
	y[6], y[7], y[8], y[9] = q.I, q.J, q.K, q.W
	rbi.integrator.SetState(tNow, y)
	return rbi.State()
}

// rates is the first-order right-hand side of the 13-state rigid-body system.
func (rbi *RigidBodyIntegrator) rates(dst, y []float64, t float64) {
	r := md3.Vec{X: y[0], Y: y[1], Z: y[2]}
	v := md3.Vec{X: y[3], Y: y[4], Z: y[5]}
	q := md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]}
	omega := md3.Vec{X: y[10], Y: y[11], Z: y[12]}

	// Translation: gravity (inertial) + body-frame force rotated to inertial.
	aGrav := gravInertial(rbi.coord, r, t)
	aExt := md3.Scale(1/rbi.mass, q.Rotate(rbi.lastForceBody))
	a := md3.Add(aGrav, aExt)

	// Attitude kinematics: q̇ = ½ q⊗[0,ω] with ω in the body frame.
	qd := q.Mul(md3.Quat{}.WithIJK(omega)).Scale(0.5)

	// Rotational dynamics (Euler's equation): ω̇ = I⁻¹(τ − ω×Iω).
	Iw := md3.MulMatVec(rbi.inertia, omega)
	wd := md3.MulMatVec(rbi.inertiaInv, md3.Sub(rbi.lastTorqueBody, md3.Cross(omega, Iw)))

	dst[0], dst[1], dst[2] = v.X, v.Y, v.Z
	dst[3], dst[4], dst[5] = a.X, a.Y, a.Z
	dst[6], dst[7], dst[8], dst[9] = qd.I, qd.J, qd.K, qd.W
	dst[10], dst[11], dst[12] = wd.X, wd.Y, wd.Z
}
