package physics

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/physics/ode"
)

// AccelSource provides the total inertial acceleration acting on the integrated
// point mass. It is evaluated once per integrator stage at inertial position sbi
// [m] and integration time t [s], and must return the full acceleration [m/s²]
// (gravity plus any external forcing) expressed in the inertial frame.
type AccelSource func(t float64, sbi md3.Vec) (abi md3.Vec)

// PointIntegrator integrates a point mass with inertial kinematics under an
// AccelSource. The acceleration is evaluated per integrator stage, so a source
// may depend on both position and time (e.g. a spherical-harmonic gravity field).
type PointIntegrator struct {
	integrator      ode.RKN1210
	integratorFast  ode.RK45
	src             AccelSource
	coordSrc        *coordAccelSource // non-nil when Configure(coord,...) is used
	lastStepWasFast bool
}

// coordAccelSource is the default AccelSource backing Configure: gravity is taken
// from a gnco.Coordinates model in the geographic frame and combined with a
// per-step external acceleration (also geographic, excluding gravity) before
// being rotated into the inertial frame. extAccelGeo is updated by Step/StepFast.
type coordAccelSource struct {
	coord       gnco.Coordinates
	extAccelGeo md3.Vec
}

func (s *coordAccelSource) AccelInertial(t float64, sbi md3.Vec) md3.Vec {
	coord := s.coord
	w := coord.World()
	TEI := w.TEI(t)
	coord.SetFromEarthFixedCoords(sbi, t)
	// Calculate TM geographic wrt earth coordinates.
	TGE := coord.TGE()
	// Calculate TM of geographic wrt inertial coordinates.
	TGI := md3.MulMat3(TGE, TEI)
	// Combine external acceleration and gravity, then transform to inertial.
	abi := md3.Add(s.extAccelGeo, coord.AGravG())
	return md3.MulMatVecTrans(TGI, abi)
}

// Configure initializes the integrator with the default geographic-gravity model
// backed by coord. t0 is the initial time, SBI0 is the initial position in the
// inertial frame, and VBI0 is the initial velocity in the inertial frame.
// Gravity is computed internally; Step/StepFast supply external acceleration in
// the geographic frame excluding gravity.
func (pi *PointIntegrator) Configure(coord gnco.Coordinates, t0 float64, SBI0, VBI0 md3.Vec) error {
	cs := &coordAccelSource{coord: coord}
	if err := pi.ConfigureSource(cs.AccelInertial, ode.Parameters{}, t0, SBI0, VBI0); err != nil {
		return err
	}
	pi.coordSrc = cs
	return nil
}

// ConfigureSource initializes the integrator with an arbitrary AccelSource and
// integrator parameters cfg. The source returns the full inertial acceleration
// per stage, so external forcing is folded into the source rather than passed to
// Step. Use StepRKN to advance and drive an adaptive substep loop.
func (pi *PointIntegrator) ConfigureSource(src AccelSource, cfg ode.Parameters, t0 float64, SBI0, VBI0 md3.Vec) error {
	*pi = PointIntegrator{src: src}
	err := pi.integrator.Configure(ode.DefaultRelaxFactor, ode.DefaultPreconditioner, cfg)
	if err != nil {
		return err
	}
	err = pi.integratorFast.Configure(cfg)
	if err != nil {
		return err
	}
	pi.integrator.Init(ode.IVP2{
		T0:   t0,
		Y0:   []float64{SBI0.X, SBI0.Y, SBI0.Z},
		DY0:  []float64{VBI0.X, VBI0.Y, VBI0.Z},
		Func: pi.accelRK12,
	})
	return nil
}

// State returns the current time, inertial position and velocity.
func (phys *PointIntegrator) State() (t float64, SBI, VBI md3.Vec) {
	t, y, dy := phys.integrator.State()
	return t, md3.Vec{X: y[0], Y: y[1], Z: y[2]}, md3.Vec{X: dy[0], Y: dy[1], Z: dy[2]}
}

// Step advances the integrator by dt using external geographic-frame acceleration.
// The supplied acceleration must exclude gravity; gravity is computed internally.
// It is only meaningful when the integrator was set up with Configure (the
// geographic-gravity model); with ConfigureSource the source provides all
// acceleration and the argument is ignored.
func (phys *PointIntegrator) Step(dt float64, externalAccelGeographicFrameNoGravity md3.Vec) (t float64, SBI, VBI md3.Vec) {
	if phys.coordSrc != nil {
		phys.coordSrc.extAccelGeo = externalAccelGeographicFrameNoGravity
	}
	phys.syncFromFast()
	phys.integrator.Step(dt)
	return phys.State()
}

// StepRKN advances the RKN12(10) integrator by dt seconds and returns the
// integrator's suggested next step [s] and any error. Acceleration comes
// entirely from the configured AccelSource; unlike Step it takes no external
// geographic acceleration and surfaces the adaptive step controller's feedback,
// which a caller driving its own substep loop needs.
func (phys *PointIntegrator) StepRKN(dt float64) (suggested float64, err error) {
	phys.syncFromFast()
	return phys.integrator.Step(dt)
}

// syncFromFast copies the RK45 fast-path state back into the RKN1210 integrator
// when the previous step was taken with StepFast, so a subsequent high-order
// step resumes from the latest state.
func (phys *PointIntegrator) syncFromFast() {
	if phys.lastStepWasFast {
		tNow, y := phys.integratorFast.State()
		phys.integrator.SetState(tNow,
			[]float64{y[0], y[1], y[2]},
			[]float64{y[3], y[4], y[5]},
		)
		phys.lastStepWasFast = false
	}
}

// StepFast advances the integrator by dt using RK45 instead of RKN1210.
// It avoids the overhead of the high-order method for situations where
// lower accuracy is acceptable. Gravity is still computed per stage.
func (phys *PointIntegrator) StepFast(dt float64, externalAccelGeographicFrameNoGravity md3.Vec) (t float64, SBI, VBI md3.Vec) {
	if phys.coordSrc != nil {
		phys.coordSrc.extAccelGeo = externalAccelGeographicFrameNoGravity
	}
	if !phys.lastStepWasFast {
		tNow, sbi, vbi := phys.integrator.State()
		phys.integratorFast.Init(ode.IVP1{
			T0:   tNow,
			Y0:   []float64{sbi[0], sbi[1], sbi[2], vbi[0], vbi[1], vbi[2]},
			Func: phys.accelFast,
		})
		phys.lastStepWasFast = true
	}
	phys.integratorFast.Step(dt)
	tFinal, y := phys.integratorFast.State()
	return tFinal,
		md3.Vec{X: y[0], Y: y[1], Z: y[2]},
		md3.Vec{X: y[3], Y: y[4], Z: y[5]}
}

func (phys *PointIntegrator) accelRK12(yppDst, y []float64, t float64) {
	ABII := phys.src(t, md3.Vec{X: y[0], Y: y[1], Z: y[2]})
	yppDst[0], yppDst[1], yppDst[2] = ABII.X, ABII.Y, ABII.Z
}

// accelFast is the rates function for integratorFast.
// State y = [x, y, z, vx, vy, vz]; dst = [vx, vy, vz, ax, ay, az].
func (phys *PointIntegrator) accelFast(dst, y []float64, t float64) {
	ABII := phys.src(t, md3.Vec{X: y[0], Y: y[1], Z: y[2]})
	dst[0], dst[1], dst[2] = y[3], y[4], y[5]
	dst[3], dst[4], dst[5] = ABII.X, ABII.Y, ABII.Z
}

// gravInertial returns the gravitational acceleration in the inertial frame at
// inertial position sbi and time t, using coord's gravity model. coord is
// stateful: SetFromEarthFixedCoords mutates it, so callers must invoke this
// sequentially (one ODE stage at a time).
func gravInertial(coord gnco.Coordinates, sbi md3.Vec, t float64) md3.Vec {
	w := coord.World()
	TEI := w.TEI(t)
	coord.SetFromEarthFixedCoords(sbi, t)
	// TM of geographic wrt inertial coordinates.
	TGI := md3.MulMat3(coord.TGE(), TEI)
	return md3.MulMatVecTrans(TGI, coord.AGravG())
}
