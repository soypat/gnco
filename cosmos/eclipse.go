package cosmos

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
