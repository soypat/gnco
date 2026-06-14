package gnco

import (
	"cmp"
	"fmt"
	"math"
	"slices"

	"github.com/soypat/geometry/md1"
	"github.com/soypat/geometry/md3"
	"github.com/soypat/geometry/ms3"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

// State is one propagated sample: the complete simulation state at one instant.
// Quantities a plot or power calculation needs — Sun direction, orbital
// elements, beta angle, face cosines, altitude — are DERIVED from these fields
// by Trajectory methods and are never stored, keeping the trajectory
// non-redundant. Attitude is carried explicitly because the orbit propagator
// integrates only position and velocity.
type State struct {
	T   cosmos.Epoch // absolute timestamp of this sample
	R   md3.Vec      // inertial position in the central body's MJ2000Eq frame [m]
	V   md3.Vec      // inertial velocity [m/s]
	Att md3.Quat     // body→inertial attitude (rotates a body-frame vector to inertial)
}

// QuatToEulerAngles returns the 3-2-1 (yaw e1 about Z, pitch e2 about Y, roll e3
// about X) intrinsic Euler angles [rad] of a attitude quaternion, the sequence
// GMAT reports by default as EulerAngle1..3. Validate the sequence against the
// GMAT attitude configuration before relying on it for a strict comparison.
func QuatToEulerAngles(q md3.Quat) (e1, e2, e3 float64) {
	// Tait-Bryan ZYX from a unit quaternion (q rotates body→inertial).
	e3 = math.Atan2(2*(q.W*q.I+q.J*q.K), 1-2*(q.I*q.I+q.J*q.J)) // roll about X
	e2 = math.Asin(clampUnit(2 * (q.W*q.J - q.K*q.I)))          // pitch about Y
	e1 = math.Atan2(2*(q.W*q.K+q.I*q.J), 1-2*(q.J*q.J+q.K*q.K)) // yaw about Z
	return e1, e2, e3
}

// AttNadirPointing returns the body→inertial attitude of a nadir-pointing
// spacecraft in the GMAT NadirPointing convention SolarCalc uses: body +Z is
// aligned with the nadir (toward the central body, −R̂), body +X with the
// velocity direction as closely as the +Z lock allows (the Velocity attitude
// constraint), and +Y completes the right-handed frame. It is derived from
// position and velocity alone and is undefined for a purely radial velocity.
func (s State) AttNadirPointing() md3.Quat {
	zb := md3.Unit(md3.Scale(-1, s.R))                            // nadir: toward the central body centre
	xb := md3.Unit(md3.Sub(s.V, md3.Scale(md3.Dot(s.V, zb), zb))) // velocity ⟂ nadir → body +X

	// Align body +Z with the nadir, then roll about it to bring body +X onto xb.
	align := md3.RotationBetweenVecs(md3.Vec{Z: 1}, zb)
	x1 := align.Rotate(md3.Vec{X: 1}) // where body +X lands after the alignment
	roll := math.Atan2(md3.Dot(zb, md3.Cross(x1, xb)), clampUnit(md3.Dot(x1, xb)))
	// q1.Mul(q2) applies q2 then q1, so this rolls after aligning.
	return md3.Rotation(roll, zb).Mul(align)
}

// OrbitElements returns the osculating Keplerian elements and true anomaly [rad] at
// sample i, derived from position and velocity.
func (s State) OrbitElements(mu float64) (orbits.Keplerian, float64, error) {
	return orbits.KeplerianFromRV(mu, s.R, s.V)
}

// OrbitSemiMajorAxis returns the osculating semi-major axis [m] at sample from the
// vis-viva relation, avoiding the full element solve.
func (s State) OrbitSemiMajorAxis(mu float64) float64 {
	r := md3.Norm(s.R)
	v2 := md3.Dot(s.V, s.V)
	return 1 / (2/r - v2/mu)
}

// OrbitPeriod returns the osculating orbital period [s] at sample i.
func (s State) OrbitPeriod(mu float64) float64 {
	a := s.OrbitSemiMajorAxis(mu)
	return 2 * math.Pi * math.Sqrt(a*a*a/mu)
}

// OrbitBetaAngle returns the orbital beta angle [rad]: the elevation of
// the Sun above the instantaneous orbit plane, sin(β) = ĥ·ŝ with ĥ the orbit
// normal and ŝ the geocentric Sun direction.
func (s State) OrbitBetaAngle(sun cosmos.Ephemeris) float64 {
	h := md3.Unit(md3.Cross(s.R, s.V))
	su := md3.Unit(sun.Position(s.T))
	return math.Asin(clampUnit(md3.Dot(h, su)))
}

// Trajectory is the propagation output: an array of self-contained samples in
// time order. Each sample carries its own absolute timestamp, so sampling may
// be non-uniform and no separate time slice is needed. The trajectory stores
// only the raw integrated state; orbital elements, Sun geometry and panel
// illumination are recovered through the methods below.
type Trajectory struct {
	Samples []State
}

// Len reports the number of samples.
func (t *Trajectory) Len() int { return len(t.Samples) }

// Span returns the elapsed seconds from the first to the last sample.
func (t *Trajectory) Span() float64 {
	if len(t.Samples) < 2 {
		return 0
	}
	return t.Samples[len(t.Samples)-1].T.Sub(t.Samples[0].T)
}

// catmullRom is the interpolating cubic used by PositionAt.
var catmullRom = ms3.SplineCatmullRom()

// PositionAt returns the inertial position interpolated to absolute epoch e
// with a Catmull-Rom spline through the surrounding samples, and reports
// whether e lies within the trajectory span.
//
// Interpolation is physical only when the sampling step resolves the orbit
// (samples spaced well under the period); on a coarse survey trajectory whose
// step exceeds the period the spline is not a meaningful orbit and the result
// must not be trusted. Positions are interpolated as offsets from the left
// bracketing sample so the float32 spline keeps sub-metre precision at the
// large absolute coordinates of an inertial frame.
func (t *Trajectory) PositionAt(e cosmos.Epoch) (md3.Vec, bool) {
	n := len(t.Samples)
	switch n {
	case 0:
		return md3.Vec{}, false
	case 1:
		ok := e.Sub(t.Samples[0].T) == 0
		return t.Samples[0].R, ok
	}
	te := e.Sub(t.Samples[0].T)
	if te < 0 || te > t.Span() {
		return md3.Vec{}, false
	}
	// Bracket: largest i with Samples[i].T <= e.
	i, found := slices.BinarySearchFunc(t.Samples, te, func(s State, target float64) int {
		return cmp.Compare(s.T.Sub(t.Samples[0].T), target)
	})
	if !found {
		i-- // BinarySearchFunc returns the first sample after e; step back to it
	}
	if i < 0 {
		i = 0
	}
	if i >= n-1 {
		return t.Samples[n-1].R, true
	}
	dt := t.Samples[i+1].T.Sub(t.Samples[i].T)
	u := float32(0)
	if dt > 0 {
		u = float32(e.Sub(t.Samples[i].T) / dt)
	}
	origin := t.Samples[i].R
	// Phantom control points outside the array are linearly extrapolated rather
	// than duplicated, which gives the Catmull-Rom end segments usable tangents.
	sampleR := func(k int) md3.Vec {
		switch {
		case k < 0:
			return md3.Sub(md3.Scale(2, t.Samples[0].R), t.Samples[1].R)
		case k >= n:
			return md3.Sub(md3.Scale(2, t.Samples[n-1].R), t.Samples[n-2].R)
		default:
			return t.Samples[k].R
		}
	}
	off := func(k int) ms3.Vec {
		d := md3.Sub(sampleR(k), origin)
		return ms3.Vec{X: float32(d.X), Y: float32(d.Y), Z: float32(d.Z)}
	}
	r := catmullRom.Evaluate(u, off(i-1), off(i), off(i+1), off(i+2))
	return md3.Add(origin, md3.Vec{X: float64(r.X), Y: float64(r.Y), Z: float64(r.Z)}), true
}

// AttitudeFunc supplies the body→inertial attitude during trajectory
// construction. The orbit propagator integrates only position and velocity, so
// attitude is provided here — e.g. nadir-pointing derived from r and v, or an
// externally propagated kinematic attitude.
type AttitudeFunc func(e cosmos.Epoch, r, v md3.Vec) md3.Quat

// Propagate runs p forward over duration seconds, sampling every step seconds
// (the initial state is the first sample), and returns the resulting
// Trajectory. Attitude for each sample comes from att; pass nil for identity.
// On a propagation error the partial trajectory and the error are returned.
func Propagate(p *OrbitPropagator, step, duration float64, att AttitudeFunc) (*Trajectory, error) {
	if step <= 0 || duration <= 0 || math.IsNaN(step) || math.IsNaN(duration) {
		return nil, fmt.Errorf("bad step %g or duration %g", step, duration)
	}
	if att == nil {
		att = func(cosmos.Epoch, md3.Vec, md3.Vec) md3.Quat { return md3.QuatIdent() }
	}
	nstep := int(math.Floor(duration/step + 1e-9))
	tr := &Trajectory{Samples: make([]State, 0, nstep+1)}
	e, r, v := p.State()
	tr.Samples = append(tr.Samples, State{T: e, R: r, V: v, Att: att(e, r, v)})
	for i := 0; i < nstep; i++ {
		e, r, v, err := p.Step(step)
		if err != nil {
			return tr, err
		}
		tr.Samples = append(tr.Samples, State{T: e, R: r, V: v, Att: att(e, r, v)})
	}
	return tr, nil
}

// Eclipses returns the penumbra and umbra phases along the trajectory in time
// order by repeatedly driving cosmos.FindNextEclipse, which interpolates
// position with PositionAt. The trajectory must resolve the orbit (see
// PositionAt). The scan resolution is the mean sample spacing.
func (t *Trajectory) Eclipses(sun cosmos.Ephemeris, occRadius, occFlattening, sunRadius float64) []cosmos.Eclipse {
	if t.Len() < 2 {
		return nil
	}
	maxStep := t.Span() / float64(t.Len()-1)
	var out []cosmos.Eclipse
	start := t.Samples[0].T
	for {
		phases, n := cosmos.FindNextEclipse(start, maxStep, t.PositionAt, sun, occRadius, occFlattening, sunRadius)
		if n == 0 {
			break
		}
		out = append(out, phases[:n]...)
		if next := phases[n-1].Exit; next.Sub(start) > 0 {
			start = next
		} else {
			break // guard against a non-advancing degenerate event
		}
	}
	return out
}

// clampUnit clamps x into [-1, 1] for safe asin/acos.
func clampUnit(x float64) float64 { return md1.Clamp(x, -1, 1) }
