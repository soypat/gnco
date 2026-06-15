package physics

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/physics/ode"
)

// RigidState bundles the kinematic state handed to a WrenchSource each
// integrator stage.
type RigidState struct {
	SBI   md3.Vec  // inertial position [m]
	VBI   md3.Vec  // inertial velocity [m/s]
	Att   md3.Quat // body→inertial attitude (unit)
	Omega md3.Vec  // body-frame angular velocity [rad/s]
}

// WrenchSource provides the total force and torque acting on the integrated
// rigid body. It is evaluated once per integrator stage at state st and absolute
// epoch, and must return the total force [N] expressed in the inertial frame
// (gravity plus any external forcing) and the total torque [N·m] expressed in
// the body frame. Because it sees the full state it may depend on attitude and
// velocity (aerodynamics, thrust, gravity-gradient, …). The epoch is the
// integrator's epoch0 anchor advanced by the elapsed integration time.
//
// This is the rigid-body analogue of AccelSource.
type WrenchSource func(st RigidState, epoch cosmos.Epoch) (forceInertial, torqueBody md3.Vec)

// RigidBodyIntegrator integrates a rigid body's full 6-DOF state —
// inertial position and velocity plus attitude (quaternion) and body-frame
// angular velocity — under a WrenchSource supplying inertial force and
// body-frame torque.
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
	integrator ode.RKF78
	src        WrenchSource
	coordSrc   *coordWrenchSource // non-nil when ConfigureCoord is used
	epoch0     cosmos.Epoch       // absolute anchor; integration time is relative to it
	mass       float64
	inertia    md3.Mat3
	inertiaInv md3.Mat3
	hNext      float64 // suggested next internal step [s]; 0 → take the full Step interval
}

// coordWrenchSource is the default WrenchSource backing ConfigureCoord: gravity
// is taken from a gnco.Coordinates model and combined with a per-step external
// body-frame force and torque. extForceBody/extTorqueBody are updated by Step and
// held constant across a step (mirroring coordAccelSource).
type coordWrenchSource struct {
	coord         gnco.Coordinates
	mass          float64
	extForceBody  md3.Vec // external body-frame force, held constant across a Step
	extTorqueBody md3.Vec // external body-frame torque, held constant across a Step
}

func (s *coordWrenchSource) Wrench(st RigidState, epoch cosmos.Epoch) (md3.Vec, md3.Vec) {
	// Gravity (inertial) plus body-frame external force rotated to inertial.
	fGrav := md3.Scale(s.mass, gravInertial(s.coord, st.SBI, epoch))
	fExt := st.Att.Rotate(s.extForceBody)
	return md3.Add(fGrav, fExt), s.extTorqueBody
}

// Configure initializes the integrator with an arbitrary WrenchSource and
// integrator parameters cfg, anchored at absolute epoch epoch0. SBI0/VBI0 are the
// initial inertial position and velocity, att0 is the initial body→inertial
// attitude (normalized internally), omega0 is the initial body-frame angular
// velocity, mass is the body mass [kg] and inertia is the body-frame inertia
// tensor [kg·m²]. The source returns the full inertial force and body torque per
// stage, so external forcing is folded into the source rather than passed to
// Step. Non-zero cfg.RelTolerance/AbsTolerance with cfg.MaxStep>0 enable adaptive
// substepping within Step.
func (rbi *RigidBodyIntegrator) Configure(src WrenchSource, cfg ode.Parameters, epoch0 cosmos.Epoch, SBI0, VBI0 md3.Vec, att0 md3.Quat, omega0 md3.Vec, mass float64, inertia md3.Mat3) error {
	*rbi = RigidBodyIntegrator{
		src:        src,
		epoch0:     epoch0,
		mass:       mass,
		inertia:    inertia,
		inertiaInv: inertia.Inverse(),
		hNext:      cfg.MaxStep, // 0 when non-adaptive → single full-dt step
	}
	err := rbi.integrator.Configure(cfg)
	if err != nil {
		return err
	}
	q := att0.Unit()
	rbi.integrator.Init(ode.IVP1{
		T0: 0, // ODE works relative to epoch0; absolute epoch rebuilt in rates.
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

// ConfigureCoord initializes the integrator with the default gravity model backed
// by coord, anchored at absolute epoch epoch0. Gravity is computed internally;
// Step supplies external body-frame force and torque. The remaining arguments
// match Configure. Stepping is non-adaptive (fixed full-dt steps).
func (rbi *RigidBodyIntegrator) ConfigureCoord(coord gnco.Coordinates, epoch0 cosmos.Epoch, SBI0, VBI0 md3.Vec, att0 md3.Quat, omega0 md3.Vec, mass float64, inertia md3.Mat3) error {
	cs := &coordWrenchSource{coord: coord, mass: mass}
	if err := rbi.Configure(cs.Wrench, ode.Parameters{}, epoch0, SBI0, VBI0, att0, omega0, mass, inertia); err != nil {
		return err
	}
	rbi.coordSrc = cs
	return nil
}

// State returns the current absolute epoch, inertial position and velocity,
// body→inertial attitude and body-frame angular velocity.
func (rbi *RigidBodyIntegrator) State() (epoch cosmos.Epoch, SBI, VBI md3.Vec, att md3.Quat, omega md3.Vec) {
	t, y := rbi.integrator.State()
	return rbi.epoch0.Add(t),
		md3.Vec{X: y[0], Y: y[1], Z: y[2]},
		md3.Vec{X: y[3], Y: y[4], Z: y[5]},
		md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]},
		md3.Vec{X: y[10], Y: y[11], Z: y[12]}
}

// Step advances the state by dt seconds, internally substepping under adaptive
// step control when configured. forceBody and torqueBody are the external
// body-frame force and torque; they are held constant across the step and are
// only meaningful when the integrator was set up with ConfigureCoord (the
// gravity model). With a custom WrenchSource the source provides all forcing and
// these arguments are ignored. The attitude quaternion is renormalized after each
// internal substep.
func (rbi *RigidBodyIntegrator) Step(dt float64, forceBody, torqueBody md3.Vec) (epoch cosmos.Epoch, SBI, VBI md3.Vec, att md3.Quat, omega md3.Vec) {
	if rbi.coordSrc != nil {
		rbi.coordSrc.extForceBody = forceBody
		rbi.coordSrc.extTorqueBody = torqueBody
	}
	t0, _ := rbi.integrator.State()
	target := t0 + dt
	// Guard against float stagnation near the target.
	const eps = 1e-9
	for {
		tNow, _ := rbi.integrator.State()
		remaining := target - tNow
		if remaining <= eps {
			break
		}
		h := remaining
		if rbi.hNext > 0 && rbi.hNext < h {
			h = rbi.hNext
		}
		hSuggest := rbi.integrator.Step(h)
		if hSuggest > 0 {
			rbi.hNext = hSuggest
		}
		// Post-substep projection: renormalize the attitude quaternion to undo the
		// small norm drift accumulated over the substep. The Baumgarte term in
		// rates keeps |q|≈1 during the substep; this projection is the hard
		// guarantee.
		tCur, y := rbi.integrator.State()
		q := md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]}.Unit()
		y[6], y[7], y[8], y[9] = q.I, q.J, q.K, q.W
		rbi.integrator.SetState(tCur, y)
	}
	return rbi.State()
}

// rates is the first-order right-hand side of the 13-state rigid-body system.
func (rbi *RigidBodyIntegrator) rates(dst, y []float64, t float64) {
	st := RigidState{
		SBI:   md3.Vec{X: y[0], Y: y[1], Z: y[2]},
		VBI:   md3.Vec{X: y[3], Y: y[4], Z: y[5]},
		Att:   md3.Quat{I: y[6], J: y[7], K: y[8], W: y[9]},
		Omega: md3.Vec{X: y[10], Y: y[11], Z: y[12]},
	}

	force, torque := rbi.src(st, rbi.epoch0.Add(t))

	// Translation: inertial acceleration a = F/m.
	a := md3.Scale(1/rbi.mass, force)

	// Attitude kinematics: q̇ = ½ q⊗[0,ω] with ω in the body frame, plus a
	// Baumgarte normalization feedback term k(1−|q|²)q that continuously drives
	// the unit-norm constraint during integration.
	q := st.Att
	qd := q.Mul(md3.Quat{}.WithIJK(st.Omega)).Scale(0.5)
	const kQuat = 1.0
	qd = qd.Add(q.Scale(kQuat * (1 - q.Dot(q))))

	// Rotational dynamics (Euler's equation): ω̇ = I⁻¹(τ − ω×Iω).
	Iw := md3.MulMatVec(rbi.inertia, st.Omega)
	wd := md3.MulMatVec(rbi.inertiaInv, md3.Sub(torque, md3.Cross(st.Omega, Iw)))

	dst[0], dst[1], dst[2] = st.VBI.X, st.VBI.Y, st.VBI.Z
	dst[3], dst[4], dst[5] = a.X, a.Y, a.Z
	dst[6], dst[7], dst[8], dst[9] = qd.I, qd.J, qd.K, qd.W
	dst[10], dst[11], dst[12] = wd.X, wd.Y, wd.Z
}
