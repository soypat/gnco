package cosmos

import "github.com/soypat/geometry/md3"

// OrientationCache memoizes the slowly-varying part of the IAU-76/FK5 Earth
// reduction (the nutation·precession matrix and the nutation angles GAST needs),
// mirroring GMAT's Nutation Update Interval. The 106-term nutation series is
// re-evaluated only when the epoch advances past intervalSec seconds since the
// last evaluation; the fast sidereal rotation (GAST) is still applied every call.
// Not safe for concurrent use: keep one cache per propagation.
type OrientationCache struct {
	intervalSec  float64
	valid        bool
	lastSecsTT   float64
	dPsi, epsBar float64  // nutation angles GAST's equation of equinoxes needs
	m            md3.Mat3 // nut · prec at lastSecsTT
}

// NewOrientationCache returns a cache that re-evaluates the nutation/precession
// reduction only every intervalSec seconds of propagation time.
func NewOrientationCache(intervalSec float64) *OrientationCache {
	return &OrientationCache{intervalSec: intervalSec}
}

// refresh recomputes M and the cached angles for the epoch at secsTT.
func (c *OrientationCache) refresh(secsTT float64) {
	tTT := secsTT / (36525 * secsPerDay)
	prec := precessionMOD(tTT)
	dPsi, dEps, epsBar := nutation1980Angles(tTT)
	nut := md3.MulMat3(rot1(-(epsBar + dEps)), md3.MulMat3(rot3(-dPsi), rot1(epsBar)))
	c.m = md3.MulMat3(nut, prec)
	c.dPsi, c.epsBar = dPsi, epsBar
	c.lastSecsTT, c.valid = secsTT, true
}
