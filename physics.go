package gnco

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/internal/ode"
)

// PhysicsPointIntegrator integrates a point mass with inertial kinematics.
// External acceleration is provided in the geographic frame without gravity.
type PhysicsPointIntegrator struct {
	integrator        ode.RKN1210
	coord             Coordinates
	lastInternalAccel md3.Vec
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
	p.integrator.Init(ode.IVP2{
		T0:   t0,
		Y0:   SBI0,
		DY0:  VBI0,
		Func: p.accel,
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
	phys.integrator.Step(dt)
	return phys.integrator.State()
}

func (phys *PhysicsPointIntegrator) accel(yppDst []md3.Vec, tv []float64, yv []md3.Vec) {
	coord := phys.coord
	w := coord.World()
	accelInternalG := phys.lastInternalAccel
	for i := range yppDst {
		t, SBII := tv[i], yv[i]
		TEI := w.TEI(t)
		coord.SetFromEarthFixedCoords(SBII, t)
		// Calculate TM geographic wrt earth coordinates.
		TGE := coord.TGE()
		// Calculate TM of geographic wrt inertial coordinates.
		TGI := md3.MulMat3(TGE, TEI)

		accelGravity := coord.AGravG()
		// Combine internal acceleration and gravity, then transform to inertial.
		ABII := md3.Add(accelInternalG, accelGravity)
		ABII = md3.MulMatVecTrans(TGI, ABII)
		yppDst[i] = ABII
	}
}
