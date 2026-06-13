package ode

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md3"
)

// Harmonic oscillator: y” = −ω²y.
// Exact solution with y(0)=1, y'(0)=0:  y(t)=cos(ωt),  y'(t)=−ω·sin(ωt).
const oscOmega = 1.0

// oscRates1 implements [y, v]' = [v, −ω²y] for RK45.
func oscRates1(dst, y []float64, _ float64) {
	dst[0] = y[1]
	dst[1] = -oscOmega * oscOmega * y[0]
}

// oscRates2 implements y” = −ω²y for RKN1210 (uses X component).
func oscRates2(ypp []md3.Vec, _ []float64, yv []md3.Vec) {
	for i := range yv {
		ypp[i] = md3.Vec{X: -oscOmega * oscOmega * yv[i].X}
	}
}

// oscExact returns the exact position and velocity at time t.
func oscExact(t float64) (y, v float64) {
	return math.Cos(oscOmega * t), -oscOmega * math.Sin(oscOmega*t)
}

// stepRK45 advances rk to tf using adaptive steps.
func stepRK45(rk *RK45, tf float64) {
	h := rk.hMax
	for {
		t, _ := rk.State()
		if t >= tf {
			break
		}
		if t+h > tf {
			h = tf - t
		}
		h = rk.Step(h)
	}
}

// stepRKF78 advances rk to tf using adaptive steps.
func stepRKF78(rk *RKF78, tf float64) {
	h := rk.hMax
	for {
		t, _ := rk.State()
		if t >= tf {
			break
		}
		if t+h > tf {
			h = tf - t
		}
		h = rk.Step(h)
	}
}

// stepVerner9 advances rk to tf using adaptive steps.
func stepVerner9(rk *Verner9, tf float64) {
	h := rk.hMax
	for {
		t, _ := rk.State()
		if t >= tf {
			break
		}
		if t+h > tf {
			h = tf - t
		}
		h = rk.Step(h)
	}
}

// stepRKN advances rk to tf using adaptive steps.
func stepRKN(t testing.TB, rk *RKN1210, tf float64) {
	t.Helper()
	h := rk.maxStep
	for {
		tNow, _, _ := rk.State()
		if tNow >= tf {
			break
		}
		if tNow+h > tf {
			h = tf - tNow
		}
		var err error
		if h, err = rk.Step(h); err != nil {
			t.Fatal("Step:", err)
		}
	}
}

// BenchmarkRK45Step exercises the adaptive Step path (the one containing the
// error-norm sqrt) by integrating the harmonic oscillator one period per Loop.
func BenchmarkRK45Step(b *testing.B) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 2 * math.Pi
	)
	var rk RK45
	if err := rk.Configure(Parameters{
		AbsTolerance: atol, RelTolerance: rtol, MinStep: 1e-8, MaxStep: 0.5,
	}); err != nil {
		b.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})
	b.ReportAllocs()
	for b.Loop() {
		rk.SetState(0, []float64{1, 0})
		stepRK45(&rk, tf)
	}
}

func TestRK45HarmonicOscillator(t *testing.T) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 10 * math.Pi // 5 full periods
	)
	var rk RK45
	if err := rk.Configure(Parameters{
		AbsTolerance: atol,
		RelTolerance: rtol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})

	stepRK45(&rk, tf)

	_, y := rk.State()
	wantY, wantV := oscExact(tf)
	errY := math.Abs(y[0] - wantY)
	errV := math.Abs(y[1] - wantV)

	const limit = 1e4 * atol
	if errY > limit {
		t.Errorf("position error %g exceeds limit %g", errY, limit)
	}
	if errV > limit {
		t.Errorf("velocity error %g exceeds limit %g", errV, limit)
	}
	t.Logf("steps=%d  errY=%.2e  errV=%.2e", rk.StepCount, errY, errV)
}

func TestRKF78HarmonicOscillator(t *testing.T) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 10 * math.Pi // 5 full periods
	)
	var rk RKF78
	if err := rk.Configure(Parameters{
		AbsTolerance: atol,
		RelTolerance: rtol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})

	stepRKF78(&rk, tf)

	_, y := rk.State()
	wantY, wantV := oscExact(tf)
	errY := math.Abs(y[0] - wantY)
	errV := math.Abs(y[1] - wantV)

	const limit = 1e4 * atol
	if errY > limit {
		t.Errorf("position error %g exceeds limit %g", errY, limit)
	}
	if errV > limit {
		t.Errorf("velocity error %g exceeds limit %g", errV, limit)
	}
	t.Logf("steps=%d  errY=%.2e  errV=%.2e", rk.StepCount, errY, errV)
}

// TestRKF78SelectInitialStep checks that SelectInitialStep returns a finite positive value.
func TestRKF78SelectInitialStep(t *testing.T) {
	var rk RKF78
	if err := rk.Configure(Parameters{
		AbsTolerance: 1e-9,
		RelTolerance: 1e-9,
		MaxStep:      1.0,
	}); err != nil {
		t.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})

	h := rk.SelectInitialStep()
	if h <= 0 || math.IsInf(h, 0) || math.IsNaN(h) {
		t.Errorf("SelectInitialStep returned %g, want finite positive", h)
	}
	t.Logf("initial step h=%.6g", h)
}

func TestVerner9HarmonicOscillator(t *testing.T) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 10 * math.Pi // 5 full periods
	)
	var rk Verner9
	if err := rk.Configure(Parameters{
		AbsTolerance: atol,
		RelTolerance: rtol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})

	stepVerner9(&rk, tf)

	_, y := rk.State()
	wantY, wantV := oscExact(tf)
	errY := math.Abs(y[0] - wantY)
	errV := math.Abs(y[1] - wantV)

	const limit = 1e4 * atol
	if errY > limit {
		t.Errorf("position error %g exceeds limit %g", errY, limit)
	}
	if errV > limit {
		t.Errorf("velocity error %g exceeds limit %g", errV, limit)
	}
	t.Logf("steps=%d  errY=%.2e  errV=%.2e", rk.StepCount, errY, errV)
}

func TestRKN1210HarmonicOscillator(t *testing.T) {
	const (
		atol = 1e-9
		tf   = 10 * math.Pi // 5 full periods
	)
	var rk RKN1210
	if err := rk.Configure(DefaultRelaxFactor, DefaultPreconditioner, Parameters{
		AbsTolerance: atol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal("Configure:", err)
	}
	rk.Init(IVP2{Y0: md3.Vec{X: 1}, DY0: md3.Vec{}, T0: 0, Func: oscRates2})

	stepRKN(t, &rk, tf)

	tFinal, yFinal, dyFinal := rk.State()
	wantY, wantV := oscExact(tFinal)
	errY := math.Abs(yFinal.X - wantY)
	errV := math.Abs(dyFinal.X - wantV)

	const limit = 1e4 * atol
	if errY > limit {
		t.Errorf("position error %g exceeds limit %g", errY, limit)
	}
	if errV > limit {
		t.Errorf("velocity error %g exceeds limit %g", errV, limit)
	}
	t.Logf("errY=%.2e  errV=%.2e", errY, errV)
}

// TestOscillatorComparison runs both integrators on identical tolerance settings
// and checks that both achieve comparable accuracy.
func TestOscillatorComparison(t *testing.T) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 2 * math.Pi // one period keeps global error predictable
	)

	// --- RK45 ---
	var rk45 RK45
	if err := rk45.Configure(Parameters{
		AbsTolerance: atol,
		RelTolerance: rtol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal(err)
	}
	rk45.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})
	stepRK45(&rk45, tf)
	_, y45 := rk45.State()
	wantY, wantV := oscExact(tf)
	errY45 := math.Abs(y45[0] - wantY)
	errV45 := math.Abs(y45[1] - wantV)

	// --- RKN1210 ---
	var rkn RKN1210
	if err := rkn.Configure(DefaultRelaxFactor, DefaultPreconditioner, Parameters{
		AbsTolerance: atol,
		MinStep:      1e-8,
		MaxStep:      0.5,
	}); err != nil {
		t.Fatal(err)
	}
	rkn.Init(IVP2{Y0: md3.Vec{X: 1}, DY0: md3.Vec{}, T0: 0, Func: oscRates2})
	stepRKN(t, &rkn, tf)
	tFinal, yFinal, dyFinal := rkn.State()
	errYN := math.Abs(yFinal.X - math.Cos(oscOmega*tFinal))
	errVN := math.Abs(dyFinal.X - (-oscOmega * math.Sin(oscOmega*tFinal)))

	const limit = 1e4 * atol
	if errY45 > limit || errV45 > limit {
		t.Errorf("RK45 errors (y=%.2e, v=%.2e) exceed %g", errY45, errV45, limit)
	}
	if errYN > limit || errVN > limit {
		t.Errorf("RKN1210 errors (y=%.2e, v=%.2e) exceed %g", errYN, errVN, limit)
	}
	t.Logf("RK45   steps=%d  errY=%.2e  errV=%.2e", rk45.StepCount, errY45, errV45)
	t.Logf("RKN1210       errY=%.2e  errV=%.2e", errYN, errVN)
}

// TestRK45SelectInitialStep checks that SelectInitialStep returns a finite positive value.
func TestRK45SelectInitialStep(t *testing.T) {
	var rk RK45
	if err := rk.Configure(Parameters{
		AbsTolerance: 1e-9,
		RelTolerance: 1e-9,
		MaxStep:      1.0,
	}); err != nil {
		t.Fatal(err)
	}
	rk.Init(IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1})

	h := rk.SelectInitialStep()
	if h <= 0 || math.IsInf(h, 0) || math.IsNaN(h) {
		t.Errorf("SelectInitialStep returned %g, want finite positive", h)
	}
	t.Logf("initial step h=%.6g", h)
}

func BenchmarkIVP_noadaptivestep(b *testing.B) {
	const step = 0.1
	params := Parameters{}
	ivp := IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1}
	ivp2 := IVP2{Y0: md3.Vec{X: 1}, DY0: md3.Vec{}, T0: 0, Func: oscRates2}
	b.Run("RK4(5)", func(b *testing.B) {
		var integ RK45
		err := integ.Configure(params)
		if err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		for b.Loop() {
			integ.Step(step)
		}
	})
	b.Run("RKF7(8)", func(b *testing.B) {
		var integ RKF78
		err := integ.Configure(params)
		if err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		for b.Loop() {
			integ.Step(step)
		}
	})
	b.Run("Verner9", func(b *testing.B) {
		var integ Verner9
		err := integ.Configure(params)
		if err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		for b.Loop() {
			integ.Step(step)
		}
	})
	b.Run("RKN12(10)", func(b *testing.B) {
		var integ RKN1210
		err := integ.Configure(DefaultRelaxFactor, DefaultPreconditioner, params)
		if err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp2)
		for b.Loop() {
			integ.Step(step)
		}
	})
}

// BenchmarkIVP_adaptiveconvergence integrates the harmonic oscillator over one period
// with adaptive step control at a fixed tolerance. Unlike the non-adaptive
// benchmark, this exercises the error-norm / step-controller path (where the
// sqrt optimization lives) and reports cost-per-accuracy: ns/op is the time to
// reach tf, while the steps/op and err metrics show how many steps and how much
// final position error each method needed to get there. Compare ns/op together
// with err — RK7(8) costs more per step but takes far fewer of them.
func BenchmarkIVP_adaptiveconvergence(b *testing.B) {
	const (
		atol = 1e-9
		rtol = 1e-9
		tf   = 2 * math.Pi // one period
	)
	params := Parameters{AbsTolerance: atol, RelTolerance: rtol, MinStep: 1e-8, MaxStep: 0.5}
	ivp := IVP1{Y0: []float64{1, 0}, T0: 0, Func: oscRates1}
	ivp2 := IVP2{Y0: md3.Vec{X: 1}, DY0: md3.Vec{}, T0: 0, Func: oscRates2}
	wantY, _ := oscExact(tf)

	b.Run("RK4(5)", func(b *testing.B) {
		var integ RK45
		if err := integ.Configure(params); err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		b.ReportAllocs()
		var errY float64
		for b.Loop() {
			integ.SetState(0, []float64{1, 0})
			integ.StepCount = 0
			stepRK45(&integ, tf)
			_, y := integ.State()
			errY = math.Abs(y[0] - wantY)
		}
		b.ReportMetric(float64(integ.StepCount), "steps/op")
		b.ReportMetric(errY*1e9, "errY×1e-9")
	})
	b.Run("RKF7(8)", func(b *testing.B) {
		var integ RKF78
		if err := integ.Configure(params); err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		b.ReportAllocs()
		var errY float64
		for b.Loop() {
			integ.SetState(0, []float64{1, 0})
			integ.StepCount = 0
			stepRKF78(&integ, tf)
			_, y := integ.State()
			errY = math.Abs(y[0] - wantY)
		}
		b.ReportMetric(float64(integ.StepCount), "steps/op")
		b.ReportMetric(errY*1e9, "errY×1e-9")
	})
	b.Run("Verner9", func(b *testing.B) {
		var integ Verner9
		if err := integ.Configure(params); err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp)
		b.ReportAllocs()
		var errY float64
		for b.Loop() {
			integ.SetState(0, []float64{1, 0})
			integ.StepCount = 0
			stepVerner9(&integ, tf)
			_, y := integ.State()
			errY = math.Abs(y[0] - wantY)
		}
		b.ReportMetric(float64(integ.StepCount), "steps/op")
		b.ReportMetric(errY*1e9, "errY×1e-9")
	})
	b.Run("RKN12(10)", func(b *testing.B) {
		var integ RKN1210
		p2 := params
		p2.RelTolerance = 0
		if err := integ.Configure(DefaultRelaxFactor, DefaultPreconditioner, p2); err != nil {
			b.Fatal(err)
		}
		integ.Init(ivp2)
		b.ReportAllocs()
		var errY float64
		for b.Loop() {
			integ.SetState(0, ivp2.Y0, ivp2.DY0)
			integ.StepCount = 0
			stepRKN(b, &integ, tf)
			_, y, _ := integ.State()
			errY = math.Abs(y.X - wantY)
		}
		b.ReportMetric(float64(integ.StepCount), "steps/op")
		b.ReportMetric(errY*1e9, "errY×1e-9")
	})
}
