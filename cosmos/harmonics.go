package cosmos

import (
	"fmt"
	"math"

	"github.com/soypat/geometry/md3"
)

// Harmonics is a spherical-harmonic gravitational potential model such as
// JGM-2, truncated on load to the requested degree and order.
//
// The Stokes coefficients C̄nm/S̄nm are not calculable: they are empirical
// geodesy products specific to each body, fit from satellite tracking
// (JGM-2: Joint Gravity Model 2, NASA GSFC/UT Austin 1994). They arrive as
// data: the built-in JGM-2 table (jgm2.go, to 21×21) or caller-supplied
// coefficients via NewHarmonics. What gnco implements is the evaluation:
// body-fixed acceleration from the coefficients via the V,W recursion of
// Montenbruck & Gill, Satellite Orbits, sec. 3.2.4.
type Harmonics struct {
	name          string
	mu            float64 // gravitational parameter from the potential file [m³/s²]
	radius        float64 // reference equatorial radius from the file [m]
	degree, order int     // truncation set at load / used in evaluation
	// Normalized Stokes coefficients in triangular-packed slices,
	// index(n,m) = n(n+1)/2 + m, n = 0..degree, m = 0..min(n, order).
	// c[index(2,0)] is the normalized C̄20 (J2 = -√5·C̄20).
	c, s []float64
	// Unnormalized coefficients used by the V,W evaluation recursion
	// (numerically safe in float64 to JGM-2's full degree 70).
	cu, su []float64
	// Scratch space for the V,W recursion, sized for degree+1.
	v, w []float64
}

// triIndex indexes triangular-packed coefficient storage.
func triIndex(n, m int) int { return n*(n+1)/2 + m }

// JGM2 returns the JGM-2 Earth gravity field truncated to degree×order from
// the built-in coefficient table in jgm2.go (available up to 21×21).
func JGM2(degree, order int) (*Harmonics, error) {
	if degree > jgm2Degree {
		return nil, fmt.Errorf("requested degree %d exceeds built-in JGM-2 table (max %d)", degree, jgm2Degree)
	}
	want := triIndex(degree, degree) + 1
	return NewHarmonics("JGM-2", jgm2Mu, jgm2Radius, degree, order, jgm2C[:want], jgm2S[:want])
}

// NewHarmonics builds a field from caller-supplied normalized coefficients in
// triangular-packed order, index(n,m) = n(n+1)/2 + m, both slices of length
// (degree+1)(degree+2)/2. The degree 0 entry must be the central term (1).
func NewHarmonics(name string, mu, refRadius float64, degree, order int, c, s []float64) (*Harmonics, error) {
	want := triIndex(degree, degree) + 1
	if degree < 0 || order < 0 || order > degree || mu <= 0 || refRadius <= 0 ||
		len(c) != want || len(s) != want {
		return nil, fmt.Errorf("bad arguments to NewHarmonics: degree=%d order=%d len(c)=%d len(s)=%d want %d", degree, order, len(c), len(s), want)
	}
	nvw := triIndex(degree+1, degree+1) + 1
	h := &Harmonics{
		name: name, mu: mu, radius: refRadius, degree: degree, order: order,
		c: append([]float64(nil), c...), s: append([]float64(nil), s...),
		cu: make([]float64, want), su: make([]float64, want),
		v: make([]float64, nvw), w: make([]float64, nvw),
	}
	h.unnormalize()
	return h, nil
}

// unnormalize fills cu, su with conventional (unnormalized) coefficients:
// Cnm = N(n,m)·C̄nm with N = √((2-δ0m)(2n+1)(n-m)!/(n+m)!).
// Montenbruck & Gill Eqn (3.13).
func (h *Harmonics) unnormalize() {
	for n := 0; n <= h.degree; n++ {
		for m := 0; m <= n; m++ {
			// prod = (n+m)!/(n-m)! computed without factorial overflow.
			prod := 1.0
			for j := n - m + 1; j <= n+m; j++ {
				prod *= float64(j)
			}
			k := 2.0
			if m == 0 {
				k = 1
			}
			N := math.Sqrt(k * float64(2*n+1) / prod)
			i := triIndex(n, m)
			h.cu[i] = N * h.c[i]
			h.su[i] = N * h.s[i]
		}
	}
}

// Name returns the gravity model name, i.e. "JGM-2".
func (h *Harmonics) Name() string { return h.name }

// Mu returns the gravitational parameter from the potential file [m³/s²].
func (h *Harmonics) Mu() float64 { return h.mu }

// ReferenceRadius returns the reference equatorial radius from the potential file [m].
func (h *Harmonics) ReferenceRadius() float64 { return h.radius }

// Degree returns the truncation degree used in evaluation.
func (h *Harmonics) Degree() int { return h.degree }

// Order returns the truncation order used in evaluation.
func (h *Harmonics) Order() int { return h.order }

// Coefficient returns the normalized Stokes coefficients C̄nm, S̄nm.
func (h *Harmonics) Coefficient(n, m int) (c, s float64) {
	if n < 0 || n > h.degree || m < 0 || m > n {
		return math.NaN(), math.NaN()
	}
	i := triIndex(n, m)
	return h.c[i], h.s[i]
}

// AccelBodyFixed returns the total gravitational acceleration [m/s²] at
// body-fixed position sBF [m], INCLUDING the central mu/r² term (the n=0
// coefficient), using the potential file's mu — GMAT's GravityField behavior.
//
// V,W recursion per Montenbruck & Gill, Satellite Orbits, Eqns (3.29)-(3.33).
// Not safe for concurrent use (shared scratch space).
func (h *Harmonics) AccelBodyFixed(sBF md3.Vec) md3.Vec {
	R := h.radius
	x, y, z := sBF.X, sBF.Y, sBF.Z
	r2 := x*x + y*y + z*z
	// Auxiliary quantities scaled by R/r².
	x0, y0, z0 := R*x/r2, R*y/r2, R*z/r2
	rho := R * R / r2
	nmax := h.degree + 1 // V,W needed one degree above the field.

	V, W := h.v, h.w
	V[0] = R / math.Sqrt(r2)
	W[0] = 0
	// Diagonal recursion, M&G Eqn (3.29).
	for m := 1; m <= nmax; m++ {
		i, ip := triIndex(m, m), triIndex(m-1, m-1)
		V[i] = float64(2*m-1) * (x0*V[ip] - y0*W[ip])
		W[i] = float64(2*m-1) * (x0*W[ip] + y0*V[ip])
	}
	// Vertical recursion, M&G Eqn (3.30)-(3.31).
	for m := 0; m <= nmax; m++ {
		for n := m + 1; n <= nmax; n++ {
			i, i1 := triIndex(n, m), triIndex(n-1, m)
			V[i] = float64(2*n-1) * z0 * V[i1]
			W[i] = float64(2*n-1) * z0 * W[i1]
			if n-m >= 2 {
				i2 := triIndex(n-2, m)
				V[i] -= float64(n+m-1) * rho * V[i2]
				W[i] -= float64(n+m-1) * rho * W[i2]
			}
			inv := 1 / float64(n-m)
			V[i] *= inv
			W[i] *= inv
		}
	}

	// Acceleration sums, M&G Eqn (3.33).
	var ax, ay, az float64
	for n := 0; n <= h.degree; n++ {
		mtop := min(n, h.order)
		for m := 0; m <= mtop; m++ {
			i := triIndex(n, m)
			c, s := h.cu[i], h.su[i]
			if c == 0 && s == 0 {
				continue
			}
			if m == 0 {
				ax += -c * V[triIndex(n+1, 1)]
				ay += -c * W[triIndex(n+1, 1)]
			} else {
				f := float64((n - m + 1) * (n - m + 2))
				ip, im := triIndex(n+1, m+1), triIndex(n+1, m-1)
				ax += 0.5 * (-c*V[ip] - s*W[ip] + f*(c*V[im]+s*W[im]))
				ay += 0.5 * (-c*W[ip] + s*V[ip] + f*(-c*W[im]+s*V[im]))
			}
			iz := triIndex(n+1, m)
			az += float64(n-m+1) * (-c*V[iz] - s*W[iz])
		}
	}
	scale := h.mu / (R * R)
	return md3.Vec{X: scale * ax, Y: scale * ay, Z: scale * az}
}
