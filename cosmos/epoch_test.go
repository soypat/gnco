package cosmos

import (
	"math"
	"testing"
	"time"
)

func TestEpochJ2000(t *testing.T) {
	// J2000.0 in UTC is 2000-01-01 11:58:55.816 (TT − 64.184 s with ΔAT=32).
	e, err := NewEpochUTC(2000, time.January, 1, 11, 58, 55.816)
	if err != nil {
		t.Fatal(err)
	}
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
	e, err := NewEpochUTC(1992, time.August, 20, 12, 14, 0)
	if err != nil {
		t.Fatal(err)
	}
	gotDeg := e.GMST() * 180 / math.Pi
	const wantDeg = 152.578788
	if math.Abs(gotDeg-wantDeg) > 1e-4 {
		t.Errorf("GMST = %.6f°, want %.6f°", gotDeg, wantDeg)
	}
}

func TestEpochStringRoundTrip(t *testing.T) {
	for _, s := range []string{
		"12 Nov 2026 21:36:00.000",
		"01 Jan 2000 11:59:28.000",
		"29 Feb 2024 23:59:59.999",
		"31 Dec 2016 23:59:30.500",
	} {
		e, err := ParseEpochUTC(s)
		if err != nil {
			t.Fatal(err)
		}
		if got := e.String(); got != s {
			t.Errorf("round trip %q → %q", s, got)
		}
	}
}

func TestDeltaATBoundary(t *testing.T) {
	before, _ := NewEpochUTC(2016, time.December, 31, 12, 0, 0)
	after, _ := NewEpochUTC(2017, time.January, 2, 12, 0, 0)
	// Two civil days apart but one extra leap second of TT elapsed.
	if d := after.Sub(before); math.Abs(d-(2*86400+1)) > 1e-9 {
		t.Errorf("TT elapsed across 2017 leap second = %v, want %v", d, 2*86400+1)
	}
}

func TestEpochAddSub(t *testing.T) {
	e, _ := ParseEpochUTC("12 Nov 2026 21:36:00.000")
	e2 := e.Add(3600)
	if e2.Sub(e) != 3600 {
		t.Error("Add/Sub mismatch")
	}
	if got := e2.String(); got != "12 Nov 2026 22:36:00.000" {
		t.Errorf("Add 3600s → %q", got)
	}
}
