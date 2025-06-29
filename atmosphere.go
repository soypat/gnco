package gnco

import (
	"math"
	"slices"

	"github.com/soypat/geometry/md1"
)

func InternationalStandardAtmosphere(zAltitude float64, T0seaLevel float64) (T, P, Rho float64) {
	const (
		g                = 9.79
		R                = 8.314472
		M                = 28.97e-3
		hydrogenAtomMass = 1.6735575e-27
	)
	var lambda, P0, rho0 float64
	switch {
	case zAltitude > 80_000:
		// Outer space medium
		P = 1.322e-11                                // https://en.wikipedia.org/wiki/Orders_of_magnitude_(pressure)
		Rho = 4 * hydrogenAtomMass * 100 * 100 * 100 // supposedly 4 atom per cm^2 in solar medium https://en.wikipedia.org/wiki/Outer_space
		T = 178                                      // Ballpark figure

	case zAltitude < 11_000:
		lambda = -6.5e-3
		P0 = 101325
		rho0 = 1.225

		T = T0seaLevel + lambda*zAltitude // linear
		P = P0 * math.Pow(T/T0seaLevel, -g*M/(R*lambda))
		Rho = rho0 * math.Pow(T/T0seaLevel, -g*M/(R*lambda)-1)

	case zAltitude < 25_000:
		P0 = 22552
		rho0 = 0.3629

		T = 216.65
		exp := math.Exp(-g * M * (zAltitude - 11000) / (R * T))
		P = P0 * exp
		Rho = rho0 * exp

	case zAltitude < 47_000:
		lambda = 3e-3
		T0 := 216.65
		P0 = 2481
		rho0 = 0.0399

		T = T0 + lambda*(zAltitude-25000)
		P = P0 * math.Pow(T/T0, -g*M/(R*lambda))
		Rho = rho0 * math.Pow(T/T0, -g*M/(R*lambda)-1)

	case zAltitude <= 80_000:
		T0 := 270.0
		lambda = (200 - T0) / (80_000 - 47_000)
		T = T0 + lambda*(zAltitude-47_000)
		// interp data from https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
		idx, _ := slices.BinarySearch(_tblAlt, zAltitude)
		if idx >= len(_tblAlt)-1 {
			P = _tblPressure[len(_tblPressure)-1]
			Rho = _tblRho[len(_tblRho)-1]
		} else {
			interp := (zAltitude - _tblAlt[idx]) / (_tblAlt[idx+1] - _tblAlt[idx])
			P = md1.Interp(_tblPressure[idx], _tblPressure[idx+1], interp)
			Rho = md1.Interp(_tblRho[idx], _tblRho[idx+1], interp)
		}
	}
	return T, P, Rho
}

var (
	_tblAlt      = []float64{4e4, 5e4, 6e4, 7e4, 8e4}
	_tblPressure = []float64{2.87e2, 7.978e1, 2.196e1, 5.2, 1.1}
	_tblRho      = []float64{3.996e-3, 1.027e-3, 3.996e-4, 8.283e-5, 1.846e-5}
)
