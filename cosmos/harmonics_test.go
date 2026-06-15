package cosmos

import (
	"math"
	"testing"

	"github.com/soypat/geometry/md3"
)

func TestJGM2Header(t *testing.T) {
	h, err := JGM2(4, 4)
	if err != nil {
		t.Fatal(err)
	}
	if h.Mu() != 3.986004415e14 {
		t.Errorf("mu = %v", h.Mu())
	}
	if h.ReferenceRadius() != 6378136.3 {
		t.Errorf("radius = %v", h.ReferenceRadius())
	}
	c20, _ := h.Coefficient(2, 0)
	if math.Abs(c20 - -4.8416539e-4) > 1e-12 {
		t.Errorf("C̄20 = %v", c20)
	}
	// J2 = -√5·C̄20 ≈ 1.0826e-3.
	if j2 := -math.Sqrt(5) * c20; math.Abs(j2-1.0826e-3) > 1e-6 {
		t.Errorf("J2 = %v", j2)
	}
	if _, err = JGM2(22, 22); err == nil {
		t.Error("expected error above embedded degree 21")
	}
}

// Degree 0 truncation must reproduce the central point mass exactly.
func TestHarmonicsCentralTerm(t *testing.T) {
	h, err := JGM2(0, 0)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range []md3.Vec{
		{X: 7000e3}, {Y: -6500e3, Z: 100e3}, {X: 1500e3, Y: 2500e3, Z: -6300e3},
	} {
		got := h.AccelBodyFixed(p)
		r := md3.Norm(p)
		want := md3.Scale(-h.Mu()/(r*r*r), p)
		if d := md3.Norm(md3.Sub(got, want)); d > 1e-12*md3.Norm(want) {
			t.Errorf("central accel at %v: got %v want %v", p, got, want)
		}
	}
}

// Degree 2 order 0 must match the closed-form J2 acceleration
// (Vallado Eqn 8-30 / Curtis Eqn 12.30).
func TestHarmonicsJ2ClosedForm(t *testing.T) {
	h, err := JGM2(2, 0)
	if err != nil {
		t.Fatal(err)
	}
	c20, _ := h.Coefficient(2, 0)
	j2 := -math.Sqrt(5) * c20 // unnormalized J2
	mu, R := h.Mu(), h.ReferenceRadius()
	for _, p := range []md3.Vec{
		{X: 6928.5e3},
		{X: 4000e3, Y: 3000e3, Z: 4500e3},
		{Z: 7000e3},
		{X: -2000e3, Y: 6300e3, Z: -2500e3},
	} {
		x, y, z := p.X, p.Y, p.Z
		r := md3.Norm(p)
		z2r2 := z * z / (r * r)
		k := 1.5 * j2 * mu * R * R / (r * r * r * r * r) // 3/2·J2·μR²/r⁵
		want := md3.Vec{
			X: -mu*x/(r*r*r) + k*x*(5*z2r2-1),
			Y: -mu*y/(r*r*r) + k*y*(5*z2r2-1),
			Z: -mu*z/(r*r*r) + k*z*(5*z2r2-3),
		}
		got := h.AccelBodyFixed(p)
		if d := md3.Norm(md3.Sub(got, want)); d > 1e-12*md3.Norm(want) {
			t.Errorf("J2 accel at %v:\n got %v\nwant %v", p, got, want)
		}
	}
}

func TestHarmonicsRoundTripNewHarmonics(t *testing.T) {
	src, _ := JGM2(4, 4)
	h, err := NewHarmonics("copy", src.Mu(), src.ReferenceRadius(), 4, 4, src.c, src.s)
	if err != nil {
		t.Fatal(err)
	}
	p := md3.Vec{X: 4000e3, Y: -3000e3, Z: 4500e3}
	if d := md3.Norm(md3.Sub(h.AccelBodyFixed(p), src.AccelBodyFixed(p))); d != 0 {
		t.Errorf("NewHarmonics copy differs by %g", d)
	}
}
