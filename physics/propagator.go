package physics

import (
	"fmt"
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/physics/ode"
)

// ForceModel accumulates the accelerations acting on an orbiting point mass,
// mirroring GMAT's ForceModel resource. Central-body gravity is always
// present: a point mass with the central body's mu until SetHarmonics
// enables a spherical-harmonic field.
type ForceModel struct {
	central   *cosmos.Body
	harmonics *cosmos.Harmonics        // nil → point-mass central body
	oriCache  *cosmos.OrientationCache // nil → recompute orientation every call
}

// NewForceModel creates a force model with point-mass gravity of the central body.
func NewForceModel(central *cosmos.Body) *ForceModel {
	if central == nil {
		panic("nil central body")
	}
	return &ForceModel{central: central}
}

// SetHarmonics enables spherical-harmonic central-body gravity, evaluated in
// the body-fixed frame rotated by the central body's orientation at epoch.
// The two-body term is then included by the harmonic evaluation using the
// potential file's mu (GMAT GravityField behavior).
func (fm *ForceModel) SetHarmonics(h *cosmos.Harmonics) { fm.harmonics = h }

// SetNutationInterval enables GMAT-style caching of the Earth orientation
// reduction: the nutation/precession matrix is re-evaluated only every
// intervalSec seconds of propagation time (GMAT's Nutation Update Interval,
// default 60 s), while the fast sidereal rotation is still applied every call.
// intervalSec <= 0 disables caching (recompute every call). When enabled, use one
// ForceModel per goroutine. Recommended: 60.
func (fm *ForceModel) SetNutationInterval(intervalSec float64) {
	if intervalSec <= 0 {
		fm.oriCache = nil
		return
	}
	fm.oriCache = cosmos.NewOrientationCache(intervalSec)
}

// Accel returns the total acceleration [m/s²] on an orbiting point mass at
// inertial (MJ2000Eq) position sBI [m] at absolute epoch e.
func (fm *ForceModel) Accel(sBI md3.Vec, e cosmos.Epoch) md3.Vec {
	if fm.harmonics != nil {
		TEI := fm.central.TEICached(fm.oriCache, e)
		sBF := md3.MulMatVec(TEI, sBI)
		aBF := fm.harmonics.AccelBodyFixed(sBF)
		return md3.MulMatVecTrans(TEI, aBF)
	}
	r := md3.Norm(sBI)
	return md3.Scale(-fm.central.Mu()/(r*r*r), sBI)
}

// PropagatorConfig mirrors GMAT's Propagator settings (SolarCalc template:
// Accuracy=1e-12, MinStep=0.001, MaxStep=2700).
type PropagatorConfig struct {
	// Accuracy is the relative per-step error tolerance, analogous to GMAT's
	// Accuracy with RSSStep error control. Zero disables adaptive stepping:
	// the propagator then advances in fixed steps of MaxStep (or the full
	// Step argument when MaxStep is zero).
	Accuracy float64
	// MinStep and MaxStep bound the internal integration step [s].
	MinStep, MaxStep float64
	// InitialStep is the first attempted internal step [s]; 0 → MaxStep.
	InitialStep float64
}

// OrbitPropagator integrates an orbiting point-mass state under a ForceModel
// in the central body's MJ2000Eq inertial frame, anchored to an absolute
// epoch. Wraps a PointIntegrator (RKN12(10) plus an RK45 fast path) with
// GMAT-like adaptive stepping.
type OrbitPropagator struct {
	fm     *ForceModel
	epoch0 cosmos.Epoch
	integ  PointIntegrator
	hNext  float64 // suggested next internal step [s]
}

// forceModelSource adapts a ForceModel anchored at epoch0 into an AccelSource:
// integration time t maps to absolute epoch epoch0+t, evaluated per stage.
type forceModelSource struct {
	fm     *ForceModel
	epoch0 cosmos.Epoch
}

func (s forceModelSource) AccelInertial(sbi md3.Vec, epoch cosmos.Epoch) md3.Vec {
	return s.fm.Accel(sbi, s.epoch0.Add(epoch.SecondsTT()))
}

// NewOrbitPropagator creates a propagator with initial inertial position
// rBI [m] and velocity vBI [m/s] at absolute epoch epoch0.
func NewOrbitPropagator(fm *ForceModel, epoch0 cosmos.Epoch, rBI, vBI md3.Vec, cfg PropagatorConfig) (*OrbitPropagator, error) {
	switch {
	case fm == nil:
		return nil, fmt.Errorf("nil force model")
	case cfg.Accuracy < 0 || cfg.MinStep < 0 || cfg.MaxStep < cfg.MinStep:
		return nil, fmt.Errorf("bad propagator config %+v", cfg)
	case cfg.Accuracy > 0 && cfg.MaxStep <= 0:
		return nil, fmt.Errorf("adaptive stepping requires MaxStep > 0")
	}
	p := &OrbitPropagator{
		fm:     fm,
		epoch0: epoch0,
	}
	err := p.integ.ConfigureSource(forceModelSource{fm: fm, epoch0: epoch0}.AccelInertial, ode.Parameters{
		RelTolerance: cfg.Accuracy,
		MinStep:      cfg.MinStep,
		MaxStep:      cfg.MaxStep,
	}, 0, rBI, vBI)
	if err != nil {
		return nil, err
	}
	p.hNext = cfg.InitialStep
	if p.hNext <= 0 {
		p.hNext = cfg.MaxStep
	}
	return p, nil
}

// State returns the current absolute epoch and inertial position [m] and
// velocity [m/s].
func (p *OrbitPropagator) State() (e cosmos.Epoch, r, v md3.Vec) {
	t, r, v := p.integ.State()
	return p.epoch0.Add(t), r, v
}

// Elapsed returns seconds integrated since the initial epoch.
func (p *OrbitPropagator) Elapsed() float64 {
	t, _, _ := p.integ.State()
	return t
}

// Step advances the state by exactly dt seconds, internally substepping with
// adaptive step control. It returns the new state; on error the state
// remains at the last successful internal step.
func (p *OrbitPropagator) Step(dt float64) (e cosmos.Epoch, r, v md3.Vec, err error) {
	if dt <= 0 || math.IsNaN(dt) {
		e, r, v = p.State()
		return e, r, v, fmt.Errorf("bad step %g", dt)
	}
	t0 := p.Elapsed()
	target := t0 + dt
	// Guard against float stagnation near the target.
	const eps = 1e-9
	for {
		remaining := target - p.Elapsed()
		if remaining <= eps {
			break
		}
		h := math.Min(p.hNext, remaining)
		if h <= 0 {
			h = remaining
		}
		hSuggest, err := p.integ.StepRKN(h)
		if err != nil {
			e, r, v = p.State()
			return e, r, v, err
		}
		if hSuggest > 0 {
			p.hNext = hSuggest
		}
	}
	e, r, v = p.State()
	return e, r, v, nil
}
