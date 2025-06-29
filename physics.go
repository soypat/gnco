package gnco

import (
	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/internal/ode"
)

type PhysicsPointIntegrator struct {
	integrator        ode.RKN1210
	coord             Coordinates
	lastInternalAccel md3.Vec
}

func NewPhysicsPointIntegrator(coord Coordinates, t0 float64, SBI0, VBI0 md3.Vec) *PhysicsPointIntegrator {
	p := &PhysicsPointIntegrator{
		coord: coord,
		integrator: *ode.NewRKN1210(ode.DefaultRelaxFactor, ode.DefaultPreconditioner, ode.Parameters{
			AbsTolerance: 0,
			MinStep:      0,
			MaxStep:      0,
		}),
	}
	p.integrator.Init(ode.IVP2{
		T0:   t0,
		Y0:   SBI0,
		DY0:  VBI0,
		Func: p.accel,
	})
	return p
}

// Step steps the physics engine with the external acceleration in geographical frame which is obtained by TVG*ABV.
// Gravity should not be included in the external acceleration as it is obtained from the coordinate system [Coordinates] AGravG method.
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
		SBIE := md3.MulMatVec(TEI, SBII)
		coord.SetFromEarthFixedCoords(SBIE, t)
		// Calculate TM geographic wrt earth coordinates.
		TGE := coord.TGE()
		// Calculate TM of geographic wrt inertial coordinates.
		TGI := md3.MulMat3(TGE, TEI)

		accelGravity := coord.AGravG()
		// TM to Geographic coordinates. ABIIn = [0 0 -Az]
		ABII := md3.Add(accelInternalG, accelGravity)
		ABII = md3.MulMatVecTrans(TGI, ABII)
		yppDst[i] = ABII
	}
}
