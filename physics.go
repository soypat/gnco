package gnco

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/internal/ode"
)

// PhysicsPointIntegrator integrates a point mass with inertial kinematics.
// External acceleration is provided in the geographic frame without gravity.
type PhysicsPointIntegrator struct {
	integrator        ode.RKN1210
	integratorFast    ode.RK45
	coord             Coordinates
	lastInternalAccel md3.Vec
	lastStepWasFast   bool
}

// NewPhysicsPointIntegrator creates and initializes a physics integrator.
// coord provides the coordinate system, t0 is the initial time, SBI0 is the initial
// position in inertial frame, and VBI0 is the initial velocity in inertial frame.
func NewPhysicsPointIntegrator(coord Coordinates, t0 float64, SBI0, VBI0 md3.Vec) *PhysicsPointIntegrator {
	p := &PhysicsPointIntegrator{
		coord: coord,
	}
	err := p.integrator.Configure(ode.DefaultRelaxFactor, ode.DefaultPreconditioner, ode.Parameters{
		AbsTolerance: 0,
		MinStep:      0,
		MaxStep:      0,
	})
	if err != nil {
		panic(err)
	}
	err = p.integratorFast.Configure(ode.Parameters{
		AbsTolerance: 0,
		MinStep:      0,
		MaxStep:      0,
	})
	if err != nil {
		panic(err)
	}
	p.integrator.Init(ode.IVP2{
		T0:   t0,
		Y0:   SBI0,
		DY0:  VBI0,
		Func: p.accelRK12,
	})
	return p
}

// State returns the current time, inertial position and velocity.
func (phys *PhysicsPointIntegrator) State() (t float64, SBI, VBI md3.Vec) {
	return phys.integrator.State()
}

// Step advances the integrator by dt using external geographic-frame acceleration.
// The supplied acceleration must exclude gravity; gravity is computed internally.
func (phys *PhysicsPointIntegrator) Step(dt float64, externalAccelGeographicFrameNoGravity md3.Vec) (t float64, SBI, VBI md3.Vec) {
	phys.lastInternalAccel = externalAccelGeographicFrameNoGravity
	if phys.lastStepWasFast {
		tNow, y := phys.integratorFast.State()
		phys.integrator.SetState(tNow,
			md3.Vec{X: y[0], Y: y[1], Z: y[2]},
			md3.Vec{X: y[3], Y: y[4], Z: y[5]},
		)
		phys.lastStepWasFast = false
	}
	phys.integrator.Step(dt)
	return phys.integrator.State()
}

// StepFast advances the integrator by dt using RK45 instead of RKN1210.
// It avoids the overhead of the high-order method for situations where
// lower accuracy is acceptable. Gravity is still computed per stage.
func (phys *PhysicsPointIntegrator) StepFast(dt float64, externalAccelGeographicFrameNoGravity md3.Vec) (t float64, SBI, VBI md3.Vec) {
	phys.lastInternalAccel = externalAccelGeographicFrameNoGravity
	if !phys.lastStepWasFast {
		tNow, sbi, vbi := phys.integrator.State()
		phys.integratorFast.Init(ode.IVP1{
			T0:   tNow,
			Y0:   []float64{sbi.X, sbi.Y, sbi.Z, vbi.X, vbi.Y, vbi.Z},
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

func (phys *PhysicsPointIntegrator) accelRK12(yppDst []md3.Vec, tv []float64, yv []md3.Vec) {
	for i := range yppDst {
		t, SBII := tv[i], yv[i]
		ABII := phys.accelMain(SBII, t)
		yppDst[i] = ABII
	}
}

// accelFast is the rates function for integratorFast.
// State y = [x, y, z, vx, vy, vz]; dst = [vx, vy, vz, ax, ay, az].
func (phys *PhysicsPointIntegrator) accelFast(dst, y []float64, t float64) {
	SBII := md3.Vec{X: y[0], Y: y[1], Z: y[2]}
	vel := md3.Vec{X: y[3], Y: y[4], Z: y[5]}
	ABII := phys.accelMain(SBII, t)
	dst[0], dst[1], dst[2] = vel.X, vel.Y, vel.Z
	dst[3], dst[4], dst[5] = ABII.X, ABII.Y, ABII.Z
}

func (phys *PhysicsPointIntegrator) accelMain(sbi md3.Vec, t float64) (abi md3.Vec) {
	coord := phys.coord
	w := coord.World()
	TEI := w.TEI(t)
	coord.SetFromEarthFixedCoords(sbi, t)
	// Calculate TM geographic wrt earth coordinates.
	TGE := coord.TGE()
	// Calculate TM of geographic wrt inertial coordinates.
	TGI := md3.MulMat3(TGE, TEI)

	accelGravity := coord.AGravG()
	// Combine internal acceleration and gravity, then transform to inertial.
	abi = md3.Add(phys.lastInternalAccel, accelGravity)
	abi = md3.MulMatVecTrans(TGI, abi)
	return abi
}
