// Package main demonstrates RKN12(10) energy conservation for a Keplerian
// elliptical orbit using large integration steps.
//
// A satellite is placed in a elliptical orbit and integrated with a fixed step of 300s,
// about five minutes per step (~19 steps per orbit).  The specific orbital energy
//
//	E = v²/2 − μ/r
//
// is a conserved quantity in Keplerian motion.  The table below tracks |ΔE/E₀|
// to quantify how well the integrator preserves this invariant over time.
package main

import (
	"errors"
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
	"github.com/soypat/gnco/physics"
)

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

type Flags struct {
	UseRK45 bool
}

func run() error {
	var flags Flags
	flag.BoolVar(&flags.UseRK45, "rk45", false, "Use faster but more error prone RK45 integrator.")
	flag.Parse()
	const (
		perigeeHASL = 400e3 // [m] 400 km: ISS-like altitude
		apogeeHASL  = 500e3 // [m] 500 km: slightly elliptical
		dt          = 300.0 // [s] integration step: 5 minutes
		nOrbits     = 5
	)
	earth := gnco.NewEarth()
	rA := earth.Radius() + earth.HASLToElevation(apogeeHASL)
	rP := earth.Radius() + earth.HASLToElevation(perigeeHASL)
	orbit, err := orbits.NewElliptical(rA, rP)
	if err != nil {
		return err
	}
	mu := earth.G() // gravitational parameter [m³/s²]
	T := orbit.Period(mu)
	E0 := orbit.SpecificEnergy(mu)

	anoP := orbit.TrueAnomalyFromElapsedSincePeriapsis(mu, 0, 1e-8)
	vR, vT := orbit.Velocity(mu, anoP)
	if vR != 0 {
		return fmt.Errorf("expected no radial velocity at periapsis, got %f", vR)
	}
	SBI0 := md3.Vec{X: rP}
	VBI0 := md3.Vec{Y: vT}
	coords := earth.GeocentricFromEarthFixedCoords(SBI0, cosmos.EpochFromTT(0))
	a := 0.5 * (orbit.Apoapsis() + orbit.Periapsis())

	fmt.Println("Keplerian orbit — RKN12(10) energy conservation")
	fmt.Printf("  Perigee: %g km   Apogee: %g km   e = %.3f   a = %.0f km\n",
		perigeeHASL/1e3, apogeeHASL/1e3, orbit.Eccentricity(), a/1e3)
	fmt.Printf("  Period T = %.0f s (%.2f h)   dt = %.0f s   %.1f steps/orbit\n\n",
		T, T/3600, dt, T/dt)

	fmt.Printf("%-8s  %-12s  %-12s  %-12s  %-12s\n", "t [h]", "radius [km]", "speed [m/s]", "rk45 |ΔE/E₀|", "rk1210 |ΔE/E₀|")
	fmt.Println("--------  ------------  ------------  ------------  -----------")
	var integrator, integratorFast physics.PointIntegrator
	err1 := integrator.Configure(&coords, 0, SBI0, VBI0)
	err2 := integratorFast.Configure(&coords, 0, SBI0, VBI0)
	if err1 != nil || err2 != nil {
		return errors.Join(err1, err2)
	}
	t, SBI, VBI := integrator.State()
	SBIfast, VBIfast := SBI, VBI // copy for fast integration comparison.
	tfast := t
	nextPrint := 0.0
	totalSteps := 0
	var maxErrE float64
	for t < float64(nOrbits)*T {
		r, v := md3.Norm(SBI), md3.Norm(VBI)
		rfast, vfast := md3.Norm(SBIfast), md3.Norm(VBIfast)
		Efast := 0.5*vfast*vfast - mu/rfast
		E := 0.5*v*v - mu/r

		errE := math.Abs((E - E0) / E0)
		errEfast := math.Abs((Efast - E0) / E0)
		if errE > maxErrE {
			maxErrE = errE
		}

		if t >= nextPrint {
			fmt.Printf("%-8.3f  %-12.1f  %-12.1f  %-13.2e %-12.2e\n",
				t/3600, r/1e3, v, errEfast, errE)
			nextPrint += T / 4 // four samples per orbit
		}
		// Zero external forces other than gravity.
		// Gravity is calculated within Step from the gnco.Coordinates system provided.
		t, SBI, VBI = integrator.Step(dt, md3.Vec{})
		tfast, SBIfast, VBIfast = integratorFast.StepFast(dt, md3.Vec{})
		totalSteps++
	}

	fmt.Println()
	fmt.Printf("Completed %.0f orbits in %d steps (dt = %.0f s). t=%.3fh tfast=%.3fh\n",
		float64(nOrbits), totalSteps, dt, t/3600, tfast/3600)
	fmt.Printf("Maximum specific energy error |ΔE/E₀|: %.2e\n", maxErrE)
	return nil
}
