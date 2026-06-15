package gnco

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

// earthMu is Earth's gravitational parameter [m³/s²], used as a generic central
// body parameter for the orbit-math tests below.
const earthMu = 3.986004415e14

// fixedSun is a stationary Sun ephemeris for geometry tests.
type fixedSun struct{ p md3.Vec }

func (f fixedSun) Position(cosmos.Epoch) md3.Vec { return f.p }

// circularTrajectory samples a circular equatorial orbit of radius r0 over one
// period into n+1 samples, attitude identity.
func circularTrajectory(r0 float64, n int) (*Trajectory, float64) {
	vc := math.Sqrt(earthMu / r0)
	period := 2 * math.Pi * math.Sqrt(r0*r0*r0/earthMu)
	tr := &Trajectory{Samples: make([]State, 0, n+1)}
	for i := 0; i <= n; i++ {
		tt := period * float64(i) / float64(n)
		th := 2 * math.Pi * float64(i) / float64(n)
		sin, cos := math.Sincos(th)
		tr.Samples = append(tr.Samples, State{
			T:   cosmos.EpochFromTT(tt),
			R:   md3.Vec{X: r0 * cos, Y: r0 * sin},
			V:   md3.Vec{X: -vc * sin, Y: vc * cos},
			Att: md3.QuatIdent(),
		})
	}
	return tr, period
}

func TestTrajectoryElements(t *testing.T) {
	k, err := orbits.NewKeplerian(7000e3, 0.01, 0.9, 0.5, 1.0)
	if err != nil {
		t.Fatal(err)
	}
	tr := &Trajectory{}
	tas := []float64{0.1, 1.2, 2.5, 4.0, 5.5}
	for i, ta := range tas {
		r, v := k.RV(earthMu, ta)
		tr.Samples = append(tr.Samples, State{T: cosmos.EpochFromTT(float64(i) * 100), R: r, V: v, Att: md3.QuatIdent()})
	}
	for i, ta := range tas {
		sample := tr.Samples[i]
		got, gotTA, err := sample.OrbitElements(earthMu)
		if err != nil {
			t.Fatalf("sample %d: %v", i, err)
		}
		if math.Abs(got.SemiMajorAxis()-k.SemiMajorAxis()) > 1 {
			t.Errorf("sample %d SMA = %g, want %g", i, got.SemiMajorAxis(), k.SemiMajorAxis())
		}
		if math.Abs(got.Eccentricity()-k.Eccentricity()) > 1e-9 {
			t.Errorf("sample %d ECC = %g, want %g", i, got.Eccentricity(), k.Eccentricity())
		}
		if math.Abs(got.Inclination()-k.Inclination()) > 1e-9 {
			t.Errorf("sample %d INC = %g, want %g", i, got.Inclination(), k.Inclination())
		}
		if d := math.Abs(gotTA - ta); d > 1e-6 {
			t.Errorf("sample %d TA = %g, want %g", i, gotTA, ta)
		}
	}
}

func TestTrajectoryPeriod(t *testing.T) {
	tr, period := circularTrajectory(7000e3, 360)
	sample := tr.Samples[0]
	if got := sample.OrbitPeriod(earthMu); math.Abs(got-period) > 1e-3 {
		t.Errorf("Period = %g, want %g", got, period)
	}
	if got := sample.OrbitSemiMajorAxis(earthMu); math.Abs(got-7000e3) > 1 {
		t.Errorf("SemiMajorAxis = %g, want 7000e3", got)
	}
}

func TestPositionAtInterpolation(t *testing.T) {
	r0 := 7000e3
	tr, period := circularTrajectory(r0, 720)
	// Interpolate at the midpoint of interior segments (Catmull-Rom end
	// segments are inherently less accurate; eclipse refinement is interior).
	for i := 100; i < 105; i++ {
		frac := (float64(i) + 0.5) / 720
		e := cosmos.EpochFromTT(period * frac)
		got, ok := tr.PositionAt(e)
		if !ok {
			t.Fatalf("segment %d: PositionAt out of range", i)
		}
		th := 2 * math.Pi * frac
		sin, cos := math.Sincos(th)
		want := md3.Vec{X: r0 * cos, Y: r0 * sin}
		if d := md3.Norm(md3.Sub(got, want)); d > 1 {
			t.Errorf("segment %d: interpolation off by %g m", i, d)
		}
	}
	if _, ok := tr.PositionAt(cosmos.EpochFromTT(-1)); ok {
		t.Error("PositionAt accepted epoch before span")
	}
}

func TestEclipseCircular(t *testing.T) {
	r0 := 7000e3
	tr, period := circularTrajectory(r0, 1440)
	sun := fixedSun{p: md3.Vec{X: cosmos.AU}} // Sun toward +X
	const earthRadius = 6371e3
	ecl := tr.Eclipses(sun, earthRadius, 0, cosmos.NewSun().Radius())
	// One shadow pass split into entry Penumbra, Umbra, exit Penumbra.
	if len(ecl) != 3 {
		t.Fatalf("got %d phases, want 3 (penumbra/umbra/penumbra)", len(ecl))
	}
	if ecl[0].Kind != cosmos.Penumbra || ecl[1].Kind != cosmos.Umbra || ecl[2].Kind != cosmos.Penumbra {
		t.Errorf("phase kinds = %v/%v/%v, want Penumbra/Umbra/Penumbra", ecl[0].Kind, ecl[1].Kind, ecl[2].Kind)
	}
	umbra := ecl[1]
	if d := umbra.Duration(); d <= 0 || d > period/2 {
		t.Errorf("umbra duration = %g s, want in (0, %g)", d, period/2)
	}
	// Total shadow should be contiguous and centred on the anti-Sun point θ = π.
	shadowDur := ecl[2].Exit.Sub(ecl[0].Enter)
	if shadowDur <= umbra.Duration() || shadowDur > period/2 {
		t.Errorf("shadow duration = %g s, want in (umbra, %g)", shadowDur, period/2)
	}
	mid := ecl[0].Enter.Add(shadowDur / 2)
	thMid := 2 * math.Pi * mid.Sub(tr.Samples[0].T) / period
	if math.Abs(thMid-math.Pi) > 0.2 {
		t.Errorf("eclipse centre θ = %g rad, want ~π", thMid)
	}
}
