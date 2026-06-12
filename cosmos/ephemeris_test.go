package cosmos

import (
	"math"
	"testing"
	"time"

	"github.com/soypat/geometry/md3"
)

// Solar invariants: equinox/solstice geometry and orbital distance bounds.
func TestAnalyticSunInvariants(t *testing.T) {
	sun := NewAnalyticSun()
	const d = 180 / math.Pi

	// March 2026 equinox: 20 Mar 2026 ~14:46 UTC. Declination crosses zero on
	// the equator OF DATE; in J2000 the 2026 equinox sits at ~-0.146° from
	// 26 years of precession. Check both: MOD-frame dec ~0, J2000 dec ~-0.146°.
	e, _ := NewEpochUTC(2026, time.March, 20, 14, 46, 0)
	p := sun.Position(e)
	dec := math.Asin(p.Z/md3.Norm(p)) * d
	if math.Abs(dec - -0.146) > 0.05 {
		t.Errorf("equinox J2000 declination = %.4f°, want ~-0.146 (precession)", dec)
	}
	tTT := e.SecondsTT() / (36525 * 86400)
	pMOD := md3.MulMatVec(precessionMOD(tTT), p)
	if decMOD := math.Asin(pMOD.Z/md3.Norm(pMOD)) * d; math.Abs(decMOD) > 0.05 {
		t.Errorf("equinox of-date declination = %.4f°, want ~0", decMOD)
	}

	// June 2026 solstice: 21 Jun 2026 ~08:25 UTC. Max declination ≈ +23.44°.
	e, _ = NewEpochUTC(2026, time.June, 21, 8, 25, 0)
	p = sun.Position(e)
	dec = math.Asin(p.Z/md3.Norm(p)) * d
	if math.Abs(dec-23.44) > 0.05 {
		t.Errorf("solstice declination = %.4f°, want ~23.44", dec)
	}

	// Distance: perihelion (early Jan) ~0.9833 AU, aphelion (early Jul) ~1.0167 AU.
	e, _ = NewEpochUTC(2026, time.January, 3, 12, 0, 0)
	if r := md3.Norm(sun.Position(e)) / AU; math.Abs(r-0.9833) > 0.002 {
		t.Errorf("perihelion distance = %.5f AU, want ~0.9833", r)
	}
	e, _ = NewEpochUTC(2026, time.July, 5, 12, 0, 0)
	if r := md3.Norm(sun.Position(e)) / AU; math.Abs(r-1.0167) > 0.002 {
		t.Errorf("aphelion distance = %.5f AU, want ~1.0167", r)
	}
}

func TestShadowGeometry(t *testing.T) {
	sunR := NewSun().Radius()
	const earthR = 6378136.3
	sunPos := md3.Vec{X: AU}
	// Directly behind Earth at LEO altitude: deep umbra.
	pen, umb := Shadow(sunPos, md3.Vec{X: -6928.5e3}, sunR, earthR)
	if pen >= 0 || umb >= 0 {
		t.Errorf("antisolar LEO point: penumbra=%.4f umbra=%.4f, want both negative", pen, umb)
	}
	// Sunlit side: fully illuminated.
	pen, umb = Shadow(sunPos, md3.Vec{X: 6928.5e3}, sunR, earthR)
	if pen <= 0 || umb <= 0 {
		t.Errorf("subsolar LEO point: penumbra=%.4f umbra=%.4f, want both positive", pen, umb)
	}
	// Quadrature: illuminated.
	pen, _ = Shadow(sunPos, md3.Vec{Y: 6928.5e3}, sunR, earthR)
	if pen <= 0 {
		t.Errorf("quadrature point: penumbra=%.4f, want positive", pen)
	}
	// Umbra cone length ~1.385e6 km: beyond it the eclipse is annular, the
	// umbra margin turns positive while penumbra stays negative.
	pen, umb = Shadow(sunPos, md3.Vec{X: -2e9}, sunR, earthR)
	if umb <= 0 || pen >= 0 {
		t.Errorf("beyond umbra cone: penumbra=%.6f umbra=%.6f, want negative/positive (annular)", pen, umb)
	}
}
