package cosmos

import "github.com/soypat/geometry/md3"

// ShadowKind classifies a satellite's illumination during an eclipse.
type ShadowKind uint8

const (
	_        ShadowKind = iota // unknown
	Sunlit                     // fully illuminated
	Penumbra                   // Sun partially occulted
	Umbra                      // Sun fully occulted
)

func (k ShadowKind) String() string {
	switch k {
	case Sunlit:
		return "Sunlit"
	case Penumbra:
		return "Penumbra"
	case Umbra:
		return "Umbra"
	default:
		return "ShadowKind(?)"
	}
}

// Eclipse is a single shadow interval. Enter and Exit are the two boundary
// crossings; Duration is derived from them rather than stored, mirroring the
// trajectory's non-redundancy rule.
type Eclipse struct {
	Enter, Exit Epoch
	Kind        ShadowKind
}

// Duration returns the eclipse length [s].
func (e Eclipse) Duration() float64 { return e.Exit.Sub(e.Enter) }

// maxEclipsePhases is the most illumination phases one contiguous shadow event
// can hold: an entry Penumbra, the Umbra, and an exit Penumbra.
const maxEclipsePhases = 3

// FindNextEclipse scans forward from start and returns the next contiguous
// shadow event as its ordered phases, with n the number of phases (0 once the
// data ends with no further event). A full pass is entry Penumbra, Umbra, exit
// Penumbra; a grazing pass that never reaches umbra is a single Penumbra. The
// fixed-size array is returned by value, so the caller drives the search (resume
// at phases[n-1].Exit) and owns any accumulation — no allocation per event.
//
// The satellite position at any epoch is posAt(e); ok=false marks the end of the
// available data and bounds the scan. posAt is sampled in increments of at most
// maxStep [s] to bracket shadow boundaries, which are then refined by bisection;
// maxStep is thus the bracketing resolution, independent of how posAt is sourced
// (it may interpolate, so sampling need not be uniform) and must resolve the
// orbit. A shadow flank thinner than maxStep may be missed. occRadius [m] and
// occFlattening describe the occulter (flattening 0 = sphere), sunRadius [m] is
// the Sun's radius, and sun gives the Sun position at any epoch.
func FindNextEclipse(start Epoch, maxStep float64, posAt func(Epoch) (md3.Vec, bool),
	sun Ephemeris, occRadius, occFlattening, sunRadius float64) (phases [maxEclipsePhases]Eclipse, n uint8) {
	if maxStep <= 0 {
		return phases, 0
	}
	// classify reports the illumination at e from the two shadow margins.
	classify := func(e Epoch) (ShadowKind, bool) {
		pos, ok := posAt(e)
		if !ok {
			return Sunlit, false
		}
		pen, umb := Shadow(sun.Position(e), pos, sunRadius, occRadius, occFlattening)
		switch {
		case umb < 0:
			return Umbra, true
		case pen < 0:
			return Penumbra, true
		default:
			return Sunlit, true
		}
	}
	// margin selects the shadow margin whose zero separates kinds a and b: the
	// umbra margin whenever Umbra is involved, else the penumbra margin.
	margin := func(a, b ShadowKind) func(Epoch) float64 {
		umbra := a == Umbra || b == Umbra
		return func(e Epoch) float64 {
			pos, _ := posAt(e)
			pen, umb := Shadow(sun.Position(e), pos, sunRadius, occRadius, occFlattening)
			if umbra {
				return umb
			}
			return pen
		}
	}

	e := start
	k, ok := classify(e)
	if !ok {
		return phases, 0
	}
	// curKind/curEnter track the phase currently open. A shadow phase already
	// active at start is clipped to start (contiguous with a prior event's Exit).
	curKind, curEnter := k, start
	for {
		prevE := e
		e = e.Add(maxStep)
		k, ok = classify(e)
		if !ok { // data ended: close any open shadow phase at the range end
			if curKind != Sunlit {
				phases[n] = Eclipse{Enter: curEnter, Exit: rangeEnd(prevE, e, posAt), Kind: curKind}
				n++
			}
			return phases, n
		}
		if k == curKind {
			continue
		}
		// The illumination may jump two levels in one step (Sunlit↔Umbra) when a
		// penumbra flank is thinner than maxStep. Walk one level at a time so the
		// hidden flank's boundaries are both refined within [prevE, e].
		for curKind != k {
			next := curKind + 1
			if k < curKind {
				next = curKind - 1
			}
			b := refineCrossing(prevE, e, margin(curKind, next))
			if curKind != Sunlit {
				if int(n) >= maxEclipsePhases {
					return phases, n // degenerate multi-umbra pass; array full
				}
				phases[n] = Eclipse{Enter: curEnter, Exit: b, Kind: curKind}
				n++
			}
			curKind, curEnter = next, b
		}
		if curKind == Sunlit && n > 0 {
			return phases, n // returned to sunlight: event complete
		}
	}
}

// refineCrossing finds the epoch in [lo, hi] where f changes sign, to ~0.1 ms,
// returning the endpoint on the non-negative (brighter) side. Classifying at
// that epoch therefore never re-enters the shadow, keeping event resume clean.
// f(lo) and f(hi) are assumed to straddle zero.
func refineCrossing(lo, hi Epoch, f func(Epoch) float64) Epoch {
	flo := f(lo)
	for k := 0; k < 60 && hi.Sub(lo) >= 1e-4; k++ {
		mid := lo.Add(hi.Sub(lo) * 0.5)
		fm := f(mid)
		if (fm < 0) == (flo < 0) {
			lo, flo = mid, fm
		} else {
			hi = mid
		}
	}
	if flo >= 0 { // lo carries sign flo; return whichever side is non-negative
		return lo
	}
	return hi
}

// rangeEnd locates the last epoch in (lo, hi] for which posAt is in range, where
// posAt is valid at lo and out of range at hi.
func rangeEnd(lo, hi Epoch, posAt func(Epoch) (md3.Vec, bool)) Epoch {
	for k := 0; k < 60 && hi.Sub(lo) >= 1e-4; k++ {
		mid := lo.Add(hi.Sub(lo) * 0.5)
		if _, ok := posAt(mid); ok {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo
}
