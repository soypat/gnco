package cosmos

import (
	"math"
	"testing"
	"time"
)

// utcGregorianLayout mirrors GMAT's UTCGregorian format for test fixtures.
const utcGregorianLayout = "02 Jan 2006 15:04:05.000"

// parseUTC parses a GMAT UTCGregorian string into an Epoch, failing the test on
// a malformed fixture.
func parseUTC(t *testing.T, s string) Epoch {
	t.Helper()
	tm, err := time.Parse(utcGregorianLayout, s)
	if err != nil {
		t.Fatalf("parse %q: %v", s, err)
	}
	return EpochFromTime(tm)
}

func TestEpochJ2000(t *testing.T) {
	// J2000.0 in UTC is 2000-01-01 11:58:55.816 (TT − 64.184 s with ΔAT=32).
	e := EpochFromTime(time.Date(2000, time.January, 1, 11, 58, 55, 816e6, time.UTC))
	if math.Abs(e.SecondsTT()) > 1e-9 {
		t.Errorf("J2000 secsTT = %g, want 0", e.SecondsTT())
	}
	if jd := e.JulianDateTT(); jd != 2451545.0 {
		t.Errorf("J2000 JD_TT = %v, want 2451545", jd)
	}
}

func TestGMSTVallado(t *testing.T) {
	// Vallado, Example 3-5: Aug 20 1992 12:14:00 UT1 → GMST 152.578788°.
	// ΔUT1=0 in this package so the UT1 instant is entered as UTC.
	e := EpochFromTime(time.Date(1992, time.August, 20, 12, 14, 0, 0, time.UTC))
	gotDeg := e.GMST() * 180 / math.Pi
	const wantDeg = 152.578788
	if math.Abs(gotDeg-wantDeg) > 1e-4 {
		t.Errorf("GMST = %.6f°, want %.6f°", gotDeg, wantDeg)
	}
}

func TestDeltaATBoundary(t *testing.T) {
	before := EpochFromTime(time.Date(2016, time.December, 31, 12, 0, 0, 0, time.UTC))
	after := EpochFromTime(time.Date(2017, time.January, 2, 12, 0, 0, 0, time.UTC))
	// Two civil days apart but one extra leap second of TT elapsed.
	if d := after.Sub(before); math.Abs(d-(2*86400+1)) > 1e-9 {
		t.Errorf("TT elapsed across 2017 leap second = %v, want %v", d, 2*86400+1)
	}
}

func TestEpochTimeRoundTrip(t *testing.T) {
	for _, s := range []string{
		"12 Nov 2026 21:36:00.000",
		"01 Jan 2000 11:59:28.000",
		"29 Feb 2024 23:59:59.999",
		"31 Dec 2016 23:59:30.500",
		"20 Aug 1992 12:14:00.000",
	} {
		e := parseUTC(t, s)
		// Time() must invert EpochFromTime to within float64 resolution.
		back := EpochFromTime(e.Time())
		if d := back.Sub(e); math.Abs(d) > 1e-6 {
			t.Errorf("%s: Epoch→Time→Epoch drifted %g s", s, d)
		}
		// Time() must land on the same UTC civil instant the string parsed to,
		// within float64's sub-µs resolution at these distances from J2000.
		want, _ := time.Parse(utcGregorianLayout, s)
		if dt := e.Time().Sub(want); math.Abs(dt.Seconds()) > 1e-6 {
			t.Errorf("%s: Time() off by %v (got %s)", s, dt, e.Time().Format(utcGregorianLayout))
		}
	}
}

func TestEpochAddSub(t *testing.T) {
	e := parseUTC(t, "12 Nov 2026 21:36:00.000")
	e2 := e.Add(3600)
	if e2.Sub(e) != 3600 {
		t.Error("Add/Sub mismatch")
	}
	if got := e2.Time().Format(utcGregorianLayout); got != "12 Nov 2026 22:36:00.000" {
		t.Errorf("Add 3600s → %q", got)
	}
}
