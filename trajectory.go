package gnco

import (
	"fmt"
	"math"
	"sort"

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

// EulerAngles returns the 3-2-1 (yaw e1 about Z, pitch e2 about Y, roll e3
// about X) intrinsic Euler angles [rad] of the stored attitude, the sequence
// GMAT reports by default as EulerAngle1..3. Validate the sequence against the
// GMAT attitude configuration before relying on it for a strict comparison.
func (s State) EulerAngles() (e1, e2, e3 float64) {
	q := s.Att
	// Tait-Bryan ZYX from a unit quaternion (q rotates body→inertial).
	e3 = math.Atan2(2*(q.W*q.I+q.J*q.K), 1-2*(q.I*q.I+q.J*q.J)) // roll about X
	e2 = math.Asin(clampUnit(2 * (q.W*q.J - q.K*q.I)))          // pitch about Y
	e1 = math.Atan2(2*(q.W*q.K+q.I*q.J), 1-2*(q.J*q.J+q.K*q.K)) // yaw about Z
	return e1, e2, e3
}

// Trajectory is the propagation output: an array of self-contained samples in
// time order. Each sample carries its own absolute timestamp, so sampling may
// be non-uniform and no separate time slice is needed. The trajectory stores
// only the raw integrated state; orbital elements, Sun geometry and panel
// illumination are recovered through the methods below.
type Trajectory struct {
	Mu      float64 // central-body gravitational parameter [m³/s²], for element & period derivations
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

// Elements returns the osculating Keplerian elements and true anomaly [rad] at
// sample i, derived from position and velocity.
func (t *Trajectory) Elements(i int) (k orbits.Keplerian, trueAnomaly float64, _ error) {
	s := t.Samples[i]
	return orbits.KeplerianFromRV(t.Mu, s.R, s.V)
}

// SemiMajorAxis returns the osculating semi-major axis [m] at sample i from the
// vis-viva relation, avoiding the full element solve.
func (t *Trajectory) SemiMajorAxis(i int) float64 {
	s := t.Samples[i]
	r := md3.Norm(s.R)
	v2 := md3.Dot(s.V, s.V)
	return 1 / (2/r - v2/t.Mu)
}

// Period returns the osculating orbital period [s] at sample i.
func (t *Trajectory) Period(i int) float64 {
	a := t.SemiMajorAxis(i)
	return 2 * math.Pi * math.Sqrt(a*a*a/t.Mu)
}

// SunDirection returns the unit vector from the satellite toward the Sun and
// the satellite–Sun distance [m] at sample i. The Sun position is computed from
// sun on demand rather than stored.
func (t *Trajectory) SunDirection(i int, sun cosmos.Ephemeris) (dir md3.Vec, dist float64) {
	s := t.Samples[i]
	rel := md3.Sub(sun.Position(s.T), s.R)
	dist = md3.Norm(rel)
	return md3.Scale(1/dist, rel), dist
}

// BetaAngle returns the orbital beta angle [rad] at sample i: the elevation of
// the Sun above the instantaneous orbit plane, sin(β) = ĥ·ŝ with ĥ the orbit
// normal and ŝ the geocentric Sun direction.
func (t *Trajectory) BetaAngle(i int, sun cosmos.Ephemeris) float64 {
	s := t.Samples[i]
	h := md3.Unit(md3.Cross(s.R, s.V))
	su := md3.Unit(sun.Position(s.T))
	return math.Asin(clampUnit(md3.Dot(h, su)))
}

// FaceCosines returns cos(angle) between each body-frame face normal and the
// satellite→Sun direction at sample i. faces are unit normals expressed in the
// body frame; a negative result means the face points away from the Sun.
// Eclipse masking (forcing shadowed values to zero) is left to the caller, so
// that the geometric cosine and the illumination state stay independent.
func (t *Trajectory) FaceCosines(i int, faces []md3.Vec, sun cosmos.Ephemeris) []float64 {
	dir, _ := t.SunDirection(i, sun)
	att := t.Samples[i].Att
	out := make([]float64, len(faces))
	for j, f := range faces {
		out[j] = md3.Dot(att.Rotate(f), dir) // f and dir are unit vectors
	}
	return out
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
	if n == 0 {
		return md3.Vec{}, false
	}
	if n == 1 {
		ok := e.Sub(t.Samples[0].T) == 0
		return t.Samples[0].R, ok
	}
	te := e.Sub(t.Samples[0].T)
	if te < 0 || te > t.Span() {
		return md3.Vec{}, false
	}
	// Bracket: largest i with Samples[i].T <= e.
	i := sort.Search(n, func(k int) bool { return t.Samples[k].T.Sub(t.Samples[0].T) > te }) - 1
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
func Propagate(p *OrbitPropagator, mu, step, duration float64, att AttitudeFunc) (*Trajectory, error) {
	if step <= 0 || duration <= 0 || math.IsNaN(step) || math.IsNaN(duration) {
		return nil, fmt.Errorf("bad step %g or duration %g", step, duration)
	}
	if att == nil {
		att = func(cosmos.Epoch, md3.Vec, md3.Vec) md3.Quat { return md3.QuatIdent() }
	}
	nstep := int(math.Floor(duration/step + 1e-9))
	tr := &Trajectory{Mu: mu, Samples: make([]State, 0, nstep+1)}
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
// order, with boundaries refined between samples by bisection. Each partial-
// shadow pass yields up to three phases: an entry Penumbra, the Umbra (when the
// Sun is fully occulted), and an exit Penumbra; a grazing pass that never
// reaches umbra is a single Penumbra. occRadius [m] and occFlattening describe
// the occulting central body (pass flattening 0 for a sphere); sunRadius [m] is
// the Sun's radius.
//
// Both the margin sampling and the boundary refinement rely on
// Trajectory.PositionAt, so the trajectory must resolve the orbit (see
// PositionAt): on an undersampled survey pass the sign changes are still
// detected but the refined boundaries are not physical.
func (t *Trajectory) Eclipses(sun cosmos.Ephemeris, occRadius, occFlattening, sunRadius float64) []cosmos.Eclipse {
	n := len(t.Samples)
	if n < 2 {
		return nil
	}
	// margin selects the penumbra (umbra=false) or umbra (umbra=true) shadow
	// margin from cosmos.Shadow; negative means inside that shadow region.
	margin := func(pos md3.Vec, e cosmos.Epoch, umbra bool) float64 {
		pen, umb := cosmos.Shadow(sun.Position(e), pos, sunRadius, occRadius, occFlattening)
		if umbra {
			return umb
		}
		return pen
	}
	// negativeIntervals brackets every span where the selected margin is
	// negative, refining each crossing on the interpolated position.
	negativeIntervals := func(umbra bool) [][2]cosmos.Epoch {
		var iv [][2]cosmos.Epoch
		prev := margin(t.Samples[0].R, t.Samples[0].T, umbra)
		in := prev < 0
		enter := t.Samples[0].T
		for i := 0; i+1 < n; i++ {
			cur := prev
			next := margin(t.Samples[i+1].R, t.Samples[i+1].T, umbra)
			prev = next
			if (cur < 0) == (next < 0) {
				continue
			}
			cross := bisectMargin(t.Samples[i].T, t.Samples[i+1].T, func(e cosmos.Epoch) float64 {
				pos, _ := t.PositionAt(e)
				return margin(pos, e, umbra)
			})
			if next < 0 {
				enter, in = cross, true
			} else {
				iv = append(iv, [2]cosmos.Epoch{enter, cross})
				in = false
			}
		}
		if in {
			iv = append(iv, [2]cosmos.Epoch{enter, t.Samples[n-1].T})
		}
		return iv
	}

	shadow := negativeIntervals(false) // penumbra-or-deeper (any occultation)
	umbra := negativeIntervals(true)   // fully occulted

	// Split each shadow pass into Penumbra / Umbra / Penumbra around its
	// contained umbra (LEO: at most one umbra per pass).
	var out []cosmos.Eclipse
	ui := 0
	for _, sh := range shadow {
		var inner *[2]cosmos.Epoch
		for ui < len(umbra) && umbra[ui][0].Sub(sh[1]) < 0 {
			if umbra[ui][0].Sub(sh[0]) >= 0 {
				inner = &umbra[ui]
			}
			ui++
		}
		if inner == nil {
			out = append(out, cosmos.Eclipse{Enter: sh[0], Exit: sh[1], Kind: cosmos.Penumbra})
			continue
		}
		out = append(out,
			cosmos.Eclipse{Enter: sh[0], Exit: inner[0], Kind: cosmos.Penumbra},
			cosmos.Eclipse{Enter: inner[0], Exit: inner[1], Kind: cosmos.Umbra},
			cosmos.Eclipse{Enter: inner[1], Exit: sh[1], Kind: cosmos.Penumbra})
	}
	return out
}

// bisectMargin finds the epoch in [lo, hi] where f changes sign, to ~0.1 ms.
// f(lo) and f(hi) are assumed to straddle zero.
func bisectMargin(lo, hi cosmos.Epoch, f func(cosmos.Epoch) float64) cosmos.Epoch {
	flo := f(lo)
	for k := 0; k < 60; k++ {
		mid := lo.Add(hi.Sub(lo) * 0.5)
		if hi.Sub(lo) < 1e-4 {
			return mid
		}
		fm := f(mid)
		if (fm < 0) == (flo < 0) {
			lo, flo = mid, fm
		} else {
			hi = mid
		}
	}
	return lo.Add(hi.Sub(lo) * 0.5)
}

// clampUnit clamps x into [-1, 1] for safe asin/acos.
func clampUnit(x float64) float64 { return math.Max(-1, math.Min(1, x)) }
