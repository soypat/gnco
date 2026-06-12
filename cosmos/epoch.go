package cosmos

import (
	"fmt"
	"math"
	"strings"
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

// NewEpochUTC builds an Epoch from a UTC Gregorian civil date (GMAT
// DateFormat=UTCGregorian). Valid for years 1972-2099 (Vallado Alg. 14 date
// range intersected with the leap second table).
func NewEpochUTC(year int, month time.Month, day, hour, min int, sec float64) (Epoch, error) {
	if year < 1972 || year > 2099 || month < 1 || month > 12 ||
		day < 1 || day > 31 || hour < 0 || hour > 23 || min < 0 || min > 59 ||
		sec < 0 || sec >= 61 || math.IsNaN(sec) {
		return Epoch{}, fmt.Errorf("invalid UTC date %d-%d-%d %d:%d:%g (years 1972-2099 supported)", year, month, day, hour, min, sec)
	}
	// Whole julian day number per Vallado Alg. 14 in exact integer arithmetic;
	// jdInt + 0.5 is the julian date of the civil date's midnight.
	y, m := year, int(month)
	jdInt := 367*y - 7*(y+(m+9)/12)/4 + 275*m/9 + day + 1721013
	// Seconds from J2000 noon to the civil date's midnight (jdInt+0.5−2451545 days), exactly.
	secsUTC := float64(jdInt-2451545)*secsPerDay + secsPerDay/2
	secsUTC += float64(hour*3600+min*60) + sec
	mjdUTC := float64(jdInt) + 0.5 - 2400000.5 // time of day irrelevant at table resolution
	return Epoch{secsTT: secsUTC + deltaAT(mjdUTC) + ttMinusTAI}, nil
}

// ParseEpochUTC parses GMAT's UTCGregorian format, e.g. "12 Nov 2026 21:36:00.000".
func ParseEpochUTC(s string) (Epoch, error) {
	var (
		day, year, hour, min int
		monStr               string
		sec                  float64
	)
	n, err := fmt.Sscanf(strings.TrimSpace(s), "%d %3s %d %d:%d:%f", &day, &monStr, &year, &hour, &min, &sec)
	if err != nil || n != 6 {
		return Epoch{}, fmt.Errorf("cannot parse UTCGregorian epoch %q: want \"DD Mon YYYY HH:MM:SS.sss\"", s)
	}
	mon := monthFromName(monStr)
	if mon == 0 {
		return Epoch{}, fmt.Errorf("cannot parse UTCGregorian epoch %q: unknown month %q", s, monStr)
	}
	return NewEpochUTC(year, mon, day, hour, min, sec)
}

var monthNames = [12]string{"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"}

func monthFromName(s string) time.Month {
	for i, name := range monthNames {
		if strings.EqualFold(s, name) {
			return time.Month(i + 1)
		}
	}
	return 0
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

// String formats the epoch as GMAT UTCGregorian, e.g. "12 Nov 2026 21:36:00.000",
// with millisecond resolution.
func (e Epoch) String() string {
	// Round to milliseconds first so 59.9996 s carries into the next minute.
	ms := math.Round(e.secsUTC() * 1000)
	jdUTC := jdJ2000 + ms/1000/secsPerDay
	// Civil date from julian day number (Fliegel & Van Flandern; Vallado Alg. 22 equivalent).
	jdn := int64(math.Floor(jdUTC + 0.5))
	tod := jdUTC + 0.5 - float64(jdn) // [days] since midnight
	l := jdn + 68569
	n := 4 * l / 146097
	l -= (146097*n + 3) / 4
	yy := 4000 * (l + 1) / 1461001
	l -= 1461*yy/4 - 31
	mm := 80 * l / 2447
	day := l - 2447*mm/80
	l = mm / 11
	month := mm + 2 - 12*l
	year := 100*(n-49) + yy + l

	todMS := int64(math.Round(tod * secsPerDay * 1000))
	if todMS >= secsPerDay*1000 {
		todMS = secsPerDay*1000 - 1 // guard float roundoff at midnight
	}
	hour := todMS / 3600000
	min := todMS / 60000 % 60
	sec := float64(todMS%60000) / 1000
	return fmt.Sprintf("%02d %s %d %02d:%02d:%06.3f", day, monthNames[month-1], year, hour, min, sec)
}
