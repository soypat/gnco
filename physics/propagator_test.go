package physics

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

// TestOrbitPropagatorClosesOrbit propagates a point-mass Keplerian orbit for
// exactly one period and checks that the state returns to its start and that
// specific orbital energy is conserved. Point-mass gravity yields a closed
// ellipse, so position and velocity must come back to the initial conditions.
func TestOrbitPropagatorClosesOrbit(t *testing.T) {
	earth := cosmos.NewEarth()
	mu := earth.Mu()

	const (
		perigee = 6378e3 + 500e3  // [m]
		apogee  = 6378e3 + 1000e3 // [m]
	)
	orbit, err := orbits.NewElliptical(apogee, perigee)
	if err != nil {
		t.Fatal(err)
	}
	rp := orbit.Periapsis()
	_, vT := orbit.Velocity(mu, 0) // tangential speed at periapsis
	r0 := md3.Vec{X: rp}
	v0 := md3.Vec{Y: vT}
	T := orbit.Period(mu)
	E0 := specificEnergy(mu, r0, v0)

	fm := NewForceModel(earth) // point-mass central body
	prop, err := NewOrbitPropagator(fm, cosmos.Epoch{}, r0, v0, PropagatorConfig{
		Accuracy: 1e-12, MinStep: 0.001, MaxStep: 2700, InitialStep: 60,
	})
	if err != nil {
		t.Fatal(err)
	}

	_, r, v, err := prop.Step(T)
	if err != nil {
		t.Fatal(err)
	}

	if got := math.Abs(prop.Elapsed() - T); got > 1e-6 {
		t.Errorf("elapsed off by %.2e s, want exact period", got)
	}
	// Position/velocity should return to start to a small fraction of orbit size.
	if rel := md3.Norm(md3.Sub(r, r0)) / rp; rel > 1e-7 {
		t.Errorf("position closure error |Δr|/rp = %.2e, want < 1e-7", rel)
	}
	if rel := md3.Norm(md3.Sub(v, v0)) / vT; rel > 1e-7 {
		t.Errorf("velocity closure error |Δv|/v = %.2e, want < 1e-7", rel)
	}
	if rel := math.Abs((specificEnergy(mu, r, v) - E0) / E0); rel > 1e-9 {
		t.Errorf("energy drift |ΔE/E₀| = %.2e, want < 1e-9", rel)
	}
}

// specificEnergy returns the two-body specific orbital energy v²/2 − μ/r.
func specificEnergy(mu float64, r, v md3.Vec) float64 {
	return 0.5*md3.Dot(v, v) - mu/md3.Norm(r)
}
