package ode

import "math"

// RKF78 is a Runge-Kutta-Fehlberg 7(8) integrator for first-order ODE systems.
// It uses the 13-stage Fehlberg tableau and advances with the 8th-order solution,
// estimating the local error from Fehlberg's embedded 7th-order term.
type RKF78 struct {
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
	k      [13][]float64
	tmp    []float64
	yNew   []float64
	errBuf []float64
}

// Configure validates and applies step-control parameters. Call before Init.
func (rk *RKF78) Configure(cfg Parameters) error {
	if err := cfg.Validate(); err != nil {
		return err
	}
	rk.atol = cfg.AbsTolerance
	rk.rtol = cfg.RelTolerance
	rk.hMax = cfg.MaxStep
	rk.hMin = cfg.MinStep
	return nil
}

// Init binds the initial value problem and allocates workspace.
func (rk *RKF78) Init(ivp IVP1) {
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
func (rk *RKF78) State() (t float64, y []float64) {
	return rk.t, rk.y
}

// SetState overwrites the current integration state.
func (rk *RKF78) SetState(t float64, y []float64) {
	rk.t = t
	copy(rk.y, y)
}

// SelectInitialStep computes a safe first step size following Hairer, Norsett &
// Wanner §II.4, the same algorithm used by scipy's select_initial_step.
func (rk *RKF78) SelectInitialStep() float64 {
	const errOrder = 7.0
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

// Step performs a single RKF78 step. Returns the suggested step size for the next call.
// In non-adaptive mode (atol=rtol=0) the step is always accepted and h is returned unchanged.
func (rk *RKF78) Step(h float64) float64 {
	const (
		safety = 0.9
		minFac = 0.2
		maxFac = 10.0
		order  = 8.0
	)
	adaptive := rk.atol > 0 || rk.rtol > 0
	t := rk.t
	y := rk.y
	n := len(y)
	fn := rk.fn

	stepRejected := false
	for {
		tNext := t + h
		hEff := tNext - t

		// Stage 0.
		fn(rk.k[0], y, t)

		// Stages 1..12.
		for i := 1; i < 13; i++ {
			for c := range n {
				sum := 0.0
				for j := range i {
					sum += fbA[i][j] * rk.k[j][c]
				}
				rk.tmp[c] = y[c] + hEff*sum
			}
			fn(rk.k[i], rk.tmp, t+fbC[i]*hEff)
		}

		// 8th-order solution and embedded 7th-order error estimate.
		for c := range n {
			sumB, sumE := 0.0, 0.0
			for i := range 13 {
				sumB += fbB[i] * rk.k[i][c]
				sumE += fbE[i] * rk.k[i][c]
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

		// Squared RMS error norm with per-component scaling — matches scipy's
		// rms_norm. We keep it squared and fold the missing sqrt into the Pow
		// exponent below, since errNorm^(-1/order) == errNorm2^(-0.5/order).
		sumSq := 0.0
		for c := range n {
			sc := rk.atol + rk.rtol*math.Max(math.Abs(y[c]), math.Abs(rk.yNew[c]))
			e := rk.errBuf[c] / sc
			sumSq += e * e
		}
		errNorm2 := sumSq / float64(n)

		if errNorm2 <= 1 {
			rk.t = tNext
			copy(rk.y, rk.yNew)
			rk.LastErrNorm = math.Sqrt(errNorm2)
			rk.StepCount++
			if errNorm2 == 0 {
				return math.Min(hEff*maxFac, rk.hMax)
			}
			factorClamped := math.Max(minFac, math.Min(maxFac, safety*math.Pow(errNorm2, -0.5/order)))
			if stepRejected {
				factorClamped = math.Min(1.0, factorClamped)
			}
			return math.Max(rk.hMin, math.Min(rk.hMax, hEff*factorClamped))
		}

		factor := math.Max(minFac, safety*math.Pow(errNorm2, -0.5/order))
		h = math.Max(rk.hMin, hEff*factor)
		stepRejected = true

		if h <= rk.hMin {
			rk.t = tNext
			copy(rk.y, rk.yNew)
			return rk.hMin
		}
	}
}

// Runge-Kutta-Fehlberg 7(8) coefficients (Fehlberg 1968, 13 stages).
// The method advances with the 8th-order weights fbB; the local error is
// estimated from the embedded 7th-order pair via fbE = b8 − b7, which reduces
// to ±41/840 on stages 0, 10, 11, 12.
var (
	fbC = [13]float64{
		0, 2.0 / 27, 1.0 / 9, 1.0 / 6, 5.0 / 12, 1.0 / 2,
		5.0 / 6, 1.0 / 6, 2.0 / 3, 1.0 / 3, 1, 0, 1,
	}

	fbA = [13][12]float64{
		{},
		{2.0 / 27},
		{1.0 / 36, 1.0 / 12},
		{1.0 / 24, 0, 1.0 / 8},
		{5.0 / 12, 0, -25.0 / 16, 25.0 / 16},
		{1.0 / 20, 0, 0, 1.0 / 4, 1.0 / 5},
		{-25.0 / 108, 0, 0, 125.0 / 108, -65.0 / 27, 125.0 / 54},
		{31.0 / 300, 0, 0, 0, 61.0 / 225, -2.0 / 9, 13.0 / 900},
		{2, 0, 0, -53.0 / 6, 704.0 / 45, -107.0 / 9, 67.0 / 90, 3},
		{-91.0 / 108, 0, 0, 23.0 / 108, -976.0 / 135, 311.0 / 54, -19.0 / 60, 17.0 / 6, -1.0 / 12},
		{2383.0 / 4100, 0, 0, -341.0 / 164, 4496.0 / 1025, -301.0 / 82, 2133.0 / 4100, 45.0 / 82, 45.0 / 164, 18.0 / 41},
		{3.0 / 205, 0, 0, 0, 0, -6.0 / 41, -3.0 / 205, -3.0 / 41, 3.0 / 41, 6.0 / 41},
		{-1777.0 / 4100, 0, 0, -341.0 / 164, 4496.0 / 1025, -289.0 / 82, 2193.0 / 4100, 51.0 / 82, 33.0 / 164, 12.0 / 41, 0, 1},
	}

	// 8th-order solution weights.
	fbB = [13]float64{
		0, 0, 0, 0, 0, 34.0 / 105, 9.0 / 35, 9.0 / 35,
		9.0 / 280, 9.0 / 280, 0, 41.0 / 840, 41.0 / 840,
	}

	// Error weights (b8 − b7); only stages 0, 10, 11, 12 are nonzero.
	fbE = [13]float64{
		-41.0 / 840, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		-41.0 / 840, 41.0 / 840, 41.0 / 840,
	}
)
