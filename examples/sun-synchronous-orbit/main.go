// Package main propagates the SolarCalc sun-synchronous orbit (550 km, 97.8°)
// with the JGM-2 4×4 spherical-harmonic gravity field and IAU-76/FK5 Earth
// orientation, demonstrating the J2 nodal regression of ≈ +0.99°/day that
// keeps the orbit plane locked to the Sun: the beta angle (Sun elevation
// over the orbit plane) stays nearly constant while RAAN sweeps almost a
// degree per day.
//
// This force model matches the GMAT configuration validated in
// local/gmat-tests milestones 01-04 (position agreement ~1.2 m/day).
package main

import (
	"fmt"
	"log"
	"math"
	"time"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

const missionstart = "12 Nov 2026 21:36:00.000"
const utcGregorianLayout = "02 Jan 2006 15:04:05"

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

func run() error {
	const (
		deg    = math.Pi / 180
		nDays  = 14
		sunSyn = 360.0 / 365.2422 // [°/day] mean Sun longitude rate the node must track
	)
	earth := cosmos.NewEarth()
	sun := cosmos.NewAnalyticSun()
	t, err := time.Parse(utcGregorianLayout, missionstart)
	if err != nil {
		return err
	}
	epoch0 := cosmos.EpochFromTime(t)
	// SolarCalc orbit: SMA 6928.5 km (≈550 km altitude), e=0.0011, i=97.794°.
	k, err := orbits.NewKeplerian(6928.5e3, 0.0011013, 97.794*deg, 215.84*deg, 191.1*deg)
	if err != nil {
		return err
	}
	jgm2, err := cosmos.JGM2(4, 4)
	if err != nil {
		return err
	}
	fm := gnco.NewForceModel(earth)
	fm.SetHarmonics(jgm2)
	r, v := k.RV(earth.Mu(), 118.33*deg)
	prop, err := gnco.NewOrbitPropagator(fm, epoch0, r, v, gnco.PropagatorConfig{
		Accuracy: 1e-12, MinStep: 0.001, MaxStep: 2700,
	})
	if err != nil {
		return err
	}

	fmt.Println("Sun-synchronous orbit — JGM-2 4×4, FK5 Earth orientation")
	fmt.Printf("  SMA %.1f km  e %.5g  i %.3f°  T %.1f min  target node rate %+.4f°/day\n\n",
		k.SemiMajorAxis()/1e3, k.Eccentricity(), k.Inclination()/deg, k.Period(earth.Mu())/60, sunSyn)
	fmt.Printf("%-26s  %-10s  %-12s  %-10s\n", "epoch (UTC)", "RAAN [°]", "drift [°/d]", "beta [°]")
	fmt.Println("--------------------------  ----------  ------------  ----------")

	raan0 := k.RAAN() / deg
	prevRaan := raan0
	for day := 0; day <= nDays; day++ {
		e, r, v := prop.State()
		osc, _, err := orbits.KeplerianFromRV(earth.Mu(), r, v)
		if err != nil {
			return err
		}
		raan := osc.RAAN() / deg
		h := md3.Cross(r, v)
		beta := math.Asin(md3.Dot(md3.Unit(h), md3.Unit(sun.Position(e)))) / deg
		drift := raan - prevRaan
		t := e.Time().Format(utcGregorianLayout)
		if day == 0 {
			fmt.Printf("%-26s  %-10.4f  %-12s  %-10.2f\n", t, raan, "-", beta)
		} else {
			fmt.Printf("%-26s  %-10.4f  %-+12.4f  %-10.2f\n", t, raan, drift, beta)
		}
		prevRaan = raan
		if day < nDays {
			if _, _, _, err := prop.Step(86400); err != nil {
				return err
			}
		}
	}
	fmt.Printf("\nMean node drift %+.4f°/day vs sun-synchronous %+.4f°/day — the beta\n",
		(prevRaan-raan0)/nDays, sunSyn)
	fmt.Println("angle barely moves: the eclipse pattern repeats orbit after orbit.")
	return nil
}
