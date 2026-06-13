package cosmos

import (
	"math"
	"time"
)

// Epoch is an absolute instant in time for astrodynamics computations, stored
// as seconds elapsed since the J2000.0 epoch (2000-01-01 12:00:00.000
// Terrestrial Time). float64 resolution at year-2030 distances from J2000 is
// ~0.1 µs, ample for GMAT comparison (GMAT reports ~ms).
//
// Time scales (Vallado sec. 3.5): UTC + ΔAT = TAI, TAI + 32.184 s = TT.
// ΔAT is the published leap second count (see leapSeconds). UT1 is assumed
// equal to UTC (ΔUT1 = 0): GMAT reads measured ΔUT1 from an EOP file,
// |ΔUT1| < 0.9 s, equivalent to under 420 m of Earth-fixed longitude.
type Epoch struct {
	secsTT float64 // seconds since J2000.0 TT
}

const (
	jdJ2000    = 2451545.0 // Julian date of J2000.0 epoch (2000-01-01 12:00:00 TT)
	secsPerDay = 86400.0
	ttMinusTAI = 32.184 // TT − TAI offset, fixed by definition [s]
)

// leapSeconds is the ΔAT = TAI−UTC step table since 1972, mirroring GMAT's
// tai-utc.dat. mjd is the UTC Modified Julian Date (JD − 2400000.5) at which
// the value takes effect. Dates before 1972 are not supported and clamp to
// the first entry.
var leapSeconds = []struct{ mjd, dat float64 }{
	{41317, 10}, {41499, 11}, {41683, 12}, {42048, 13}, {42413, 14},
	{42778, 15}, {43144, 16}, {43509, 17}, {43874, 18}, {44239, 19},
	{44786, 20}, {45151, 21}, {45516, 22}, {46247, 23}, {47161, 24},
	{47892, 25}, {48257, 26}, {48804, 27}, {49169, 28}, {49534, 29},
	{50083, 30}, {50630, 31}, {51179, 32}, {53736, 33}, {54832, 34},
	{56109, 35}, {57204, 36}, {57754, 37},
}

// deltaAT returns TAI−UTC [s] for a UTC Modified Julian Date.
func deltaAT(mjdUTC float64) float64 {
	for i := len(leapSeconds) - 1; i > 0; i-- {
		if mjdUTC >= leapSeconds[i].mjd {
			return leapSeconds[i].dat
		}
	}
	return leapSeconds[0].dat
}

var j2000UTCNoon = time.Date(2000, time.January, 1, 12, 0, 0, 0, time.UTC)

// J2000UTCNoon is the J2000 reference instant expressed on the UTC civil clock
// (a uniform 86400 s/day count, no leap seconds). Counting Go time.Duration from
// here and then adding deltaAT + ttMinusTAI yields seconds since J2000.0 TT.
func J2000UTCNoon() time.Time {
	return j2000UTCNoon
}

// EpochFromTime builds an Epoch from a time.Time, interpreting it as a UTC
// civil instant. Any location is normalized to UTC and the monotonic clock
// reading, if present, is irrelevant (the civil/wall-clock value is used).
//
// time.Time follows the Unix convention of 86400 s/day and cannot represent a
// leap second (no 23:59:60); the leap-second bridge to TAI/TT lives here via the
// deltaAT table. The result is therefore exact except at a leap second instant,
// where the time.Time input is itself ill-defined.
func EpochFromTime(t time.Time) Epoch {
	u := t.UTC()
	// Uniform (leap-second-free) UTC seconds from the J2000 civil reference.
	secsUTC := u.Sub(j2000UTCNoon).Seconds()
	mjdUTC := jdJ2000 + secsUTC/secsPerDay - 2400000.5 // table resolution insensitive to time of day
	return Epoch{secsTT: secsUTC + deltaAT(mjdUTC) + ttMinusTAI}
}

// EpochFromTT builds an Epoch directly from seconds elapsed since J2000.0 TT.
func EpochFromTT(secondsSinceJ2000TT float64) Epoch {
	return Epoch{secsTT: secondsSinceJ2000TT}
}

// SecondsTT returns seconds elapsed since J2000.0 (2000-01-01 12:00:00.000 TT).
func (e Epoch) SecondsTT() float64 { return e.secsTT }

// JulianDateTT returns the julian date in the TT time scale [days].
func (e Epoch) JulianDateTT() float64 { return jdJ2000 + e.secsTT/secsPerDay }

// secsUTC returns seconds since J2000 in the UTC scale. ΔAT is looked up with
// the TT-based MJD then refined with the resulting UTC-based MJD, which
// corrects instants falling in the ~69 s before a leap second boundary.
func (e Epoch) secsUTC() float64 {
	mjdTT := e.JulianDateTT() - 2400000.5
	dat := deltaAT(mjdTT)
	mjdUTC := mjdTT - (ttMinusTAI+dat)/secsPerDay
	dat = deltaAT(mjdUTC)
	return e.secsTT - ttMinusTAI - dat
}

// JulianDateUT1 returns the julian date in the UT1 time scale [days],
// under the package assumption ΔUT1 = 0 (UT1 = UTC).
func (e Epoch) JulianDateUT1() float64 { return jdJ2000 + e.secsUTC()/secsPerDay }

// GMST returns the Greenwich mean sidereal time [rad] in range [0, 2π).
// IAU-82 model, Vallado Alg. 15 / Eqn (3-47).
func (e Epoch) GMST() float64 {
	// Julian centuries of UT1 since J2000.
	t := (e.JulianDateUT1() - jdJ2000) / 36525
	// GMST in seconds, where 86400 s correspond to one full turn.
	gmst := 67310.54841 + ((876600*3600+8640184.812866)+(0.093104-6.2e-6*t)*t)*t
	gmst = math.Mod(gmst, secsPerDay)
	if gmst < 0 {
		gmst += secsPerDay
	}
	return gmst * (2 * math.Pi / secsPerDay)
}

// Add returns the epoch advanced by the given seconds (negative to go back).
func (e Epoch) Add(seconds float64) Epoch { return Epoch{secsTT: e.secsTT + seconds} }

// Sub returns e − o in seconds.
func (e Epoch) Sub(o Epoch) float64 { return e.secsTT - o.secsTT }

// Time returns the epoch as a UTC time.Time, the inverse of EpochFromTime.
// The conversion goes through the UTC scale (TT − ttMinusTAI − ΔAT) using the
// same leap-second handling as EpochFromTime, so it is exact away from a leap second.
// Two caveats are inherent to time.Time: an instant on a leap second maps to the
// following second (no 23:59:60), and sub-microsecond detail below float64's
// ~0.1 µs resolution near year 2030 is not recoverable.
func (e Epoch) Time() time.Time {
	secsUTC := e.secsUTC()
	// Split whole/fractional seconds so the large integer part stays exact and
	// only the sub-second remainder is rounded to time.Time's ns granularity.
	whole, frac := math.Modf(secsUTC)
	d := time.Duration(int64(whole))*time.Second + time.Duration(math.Round(frac*1e9))*time.Nanosecond
	return j2000UTCNoon.Add(d)
}
