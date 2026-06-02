package ode

import "math"

// IVP1 defines a first-order initial value problem with an n-dimensional slice state.
type IVP1 struct {
	Y0 []float64
	T0 float64
	// Func evaluates the derivative: dst[i] = f_i(t, y).
	Func func(dst, y []float64, t float64)
}

// RK45 is a Runge-Kutta-Dormand-Prince 4(5) integrator for first-order ODE systems.
// Its step controller matches scipy's RK45 exactly.
type RK45 struct {
	// Step control — set before calling Init or Step.
	hMin, hMax float64
	atol, rtol float64

	// Diagnostics.
	LastErrNorm float64
	StepCount   int

	t  float64
	y  []float64
	fn func(dst, y []float64, t float64)

	// Workspace (allocated once by Init).
	k      [7][]float64
	tmp    []float64
	yNew   []float64
	errBuf []float64
}

func (rk45 *RK45) Configure(cfg Parameters) error {
	if err := cfg.Validate(); err != nil {
		return err
	}
	rk45.atol = cfg.AbsTolerance
	rk45.rtol = cfg.RelTolerance
	rk45.hMax = cfg.MaxStep
	rk45.hMin = cfg.MinStep
	return nil
}

func (rk *RK45) Init(ivp IVP1) {
	n := len(ivp.Y0)
	rk.t = ivp.T0
	rk.y = append(make([]float64, 0, n), ivp.Y0...)
	rk.fn = ivp.Func
	for i := range rk.k {
		rk.k[i] = make([]float64, n)
	}
	rk.tmp = make([]float64, n)
	rk.yNew = make([]float64, n)
	rk.errBuf = make([]float64, n)
	rk.LastErrNorm = 0
	rk.StepCount = 0
}

// State returns the current time and a view of the current state vector.
func (rk *RK45) State() (t float64, y []float64) {
	return rk.t, rk.y
}

// SetState overwrites the current integration state.
func (rk *RK45) SetState(t float64, y []float64) {
	rk.t = t
	copy(rk.y, y)
}

// SelectInitialStep computes a safe first step size following Hairer, Norsett &
// Wanner §II.4, the same algorithm used by scipy's select_initial_step.
func (rk *RK45) SelectInitialStep() float64 {
	const errOrder = 4.0
	n := len(rk.y)
	nf := float64(n)
	fn := rk.fn
	y := rk.y
	t := rk.t
	atol, rtol := rk.atol, rk.rtol

	f0 := make([]float64, n)
	fn(f0, y, t)

	d0, d1 := 0.0, 0.0
	for i, yi := range y {
		sc := atol + math.Abs(yi)*rtol
		d0 += (yi / sc) * (yi / sc)
		d1 += (f0[i] / sc) * (f0[i] / sc)
	}
	d0 = math.Sqrt(d0 / nf)
	d1 = math.Sqrt(d1 / nf)

	var h0 float64
	if d0 < 1e-5 || d1 < 1e-5 {
		h0 = 1e-6
	} else {
		h0 = 0.01 * d0 / d1
	}

	// Euler probe step to estimate second derivative.
	y1 := make([]float64, n)
	for i, yi := range y {
		y1[i] = yi + h0*f0[i]
	}
	f1 := make([]float64, n)
	fn(f1, y1, t+h0)

	d2 := 0.0
	for i := range f1 {
		sc := atol + math.Abs(y[i])*rtol
		dd := (f1[i] - f0[i]) / sc
		d2 += dd * dd
	}
	d2 = math.Sqrt(d2/nf) / h0

	var h1 float64
	if maxD := math.Max(d1, d2); maxD <= 1e-5 {
		h1 = math.Max(1e-6, h0*1e-3)
	} else {
		h1 = math.Pow(0.01/maxD, 1.0/(errOrder+1))
	}
	return math.Min(100*h0, h1)
}

// Step performs a single RK45 step. Returns the suggested step size for the next call.
// In non-adaptive mode (atol=rtol=0) the step is always accepted and h is returned unchanged.
func (rk *RK45) Step(h float64) float64 {
	const (
		safety = 0.9
		minFac = 0.2
		maxFac = 10.0
		order  = 5.0
	)
	adaptive := rk.atol > 0 || rk.rtol > 0
	t := rk.t
	y := rk.y
	n := len(y)
	fn := rk.fn
	sqrtN := math.Sqrt(float64(n))

	stepRejected := false
	for {
		tNext := t + h
		hEff := tNext - t

		// Stage 0.
		fn(rk.k[0], y, t)

		// Stages 1..6.
		for i := 1; i < 7; i++ {
			for c := range n {
				sum := 0.0
				for j := range i {
					sum += dpA[i][j] * rk.k[j][c]
				}
				rk.tmp[c] = y[c] + hEff*sum
			}
			fn(rk.k[i], rk.tmp, t+dpC[i]*hEff)
		}

		// 5th-order solution and error estimate.
		for c := range n {
			sumB, sumE := 0.0, 0.0
			for i := range 7 {
				sumB += dpB[i] * rk.k[i][c]
				sumE += dpE[i] * rk.k[i][c]
			}
			rk.yNew[c] = y[c] + hEff*sumB
			rk.errBuf[c] = hEff * sumE
		}

		if !adaptive {
			rk.t = tNext
			copy(rk.y, rk.yNew)
			rk.StepCount++
			return hEff
		}

		// RMS error norm with per-component scaling — matches scipy's rms_norm.
		sumSq := 0.0
		for c := range n {
			sc := rk.atol + rk.rtol*math.Max(math.Abs(y[c]), math.Abs(rk.yNew[c]))
			e := rk.errBuf[c] / sc
			sumSq += e * e
		}
		errNorm := math.Sqrt(sumSq) / sqrtN

		if errNorm <= 1 {
			rk.t = tNext
			copy(rk.y, rk.yNew)
			rk.LastErrNorm = errNorm
			rk.StepCount++
			if errNorm == 0 {
				return math.Min(hEff*maxFac, rk.hMax)
			}
			factorClamped := math.Max(minFac, math.Min(maxFac, safety*math.Pow(errNorm, -1.0/order)))
			if stepRejected {
				factorClamped = math.Min(1.0, factorClamped)
			}
			return math.Max(rk.hMin, math.Min(rk.hMax, hEff*factorClamped))
		}

		factor := math.Max(minFac, safety*math.Pow(errNorm, -1.0/order))
		h = math.Max(rk.hMin, hEff*factor)
		stepRejected = true

		if h <= rk.hMin {
			rk.t = tNext
			copy(rk.y, rk.yNew)
			return rk.hMin
		}
	}
}

// Dormand-Prince RK45 coefficients — match scipy's RK45 exactly.
var (
	dpC = [7]float64{0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1, 1}

	dpA = [7][6]float64{
		{},
		{1.0 / 5},
		{3.0 / 40, 9.0 / 40},
		{44.0 / 45, -56.0 / 15, 32.0 / 9},
		{19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
		{9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
		{35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84},
	}

	dpB = [7]float64{35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0}

	// Error coefficients in direct form (avoids cancellation from b-b* form).
	dpE = [7]float64{
		71.0 / 57600,
		0,
		-71.0 / 16695,
		71.0 / 1920,
		-17253.0 / 339200,
		22.0 / 525,
		-1.0 / 40,
	}
)
