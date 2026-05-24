// Package main demonstrates RKN12(10) energy conservation for a Keplerian
// elliptical orbit using large integration steps.
//
// A satellite is placed in a highly elliptical orbit (e ≈ 0.59, perigee 400 km,
// apogee 20 000 km) and integrated for five full orbits with a fixed step of 300 s
// — about five minutes per step (~71 steps per orbit).  The specific orbital energy
//
//	E = v²/2 − μ/r
//
// is a conserved quantity in Keplerian motion.  The table below tracks |ΔE/E₀|
// to quantify how well the integrator preserves this invariant over time.
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/orbits"
)

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

func run() error {
	const (
		perigeeHASL = 400e3 // [m] 400 km: ISS-like altitude
		apogeeHASL  = 500e3 // [m] 500 km: slightly elliptical
		dt          = 300.0 // [s] integration step: 5 minutes
		nOrbits     = 5
	)
	earth := gnco.NewEarth()
	rA := earth.Radius + earth.HASLToElevation(apogeeHASL)
	rP := earth.Radius + earth.HASLToElevation(perigeeHASL)
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
	coords := earth.GeocentricFromEarthFixedCoords(SBI0, 0)
	a := 0.5 * (orbit.Apoapsis() + orbit.Periapsis())

	fmt.Println("Keplerian orbit — RKN12(10) energy conservation")
	fmt.Printf("  Perigee: %g km   Apogee: %g km   e = %.3f   a = %.0f km\n",
		perigeeHASL/1e3, apogeeHASL/1e3, orbit.Eccentricity(), a/1e3)
	fmt.Printf("  Period T = %.0f s (%.2f h)   dt = %.0f s   %.1f steps/orbit\n\n",
		T, T/3600, dt, T/dt)

	fmt.Printf("%-8s  %-12s  %-12s  %-12s\n", "t [h]", "radius [km]", "speed [m/s]", "|ΔE/E₀|")
	fmt.Println("--------  ------------  ------------  ------------")
	SBI0 = md3.Vec{X: 6.771146e+06}
	VBI0 = md3.Vec{Y: 7700.584943721379}
	integrator := gnco.NewPhysicsPointIntegrator(&coords, 0, SBI0, VBI0)
	t, SBI, VBI := integrator.State()
	nextPrint := 0.0
	totalSteps := 0
	var maxErrE float64
	for t < float64(nOrbits)*T {
		r := md3.Norm(SBI)
		v := md3.Norm(VBI)
		E := 0.5*v*v - mu/r
		errE := math.Abs((E - E0) / E0)
		if errE > maxErrE {
			maxErrE = errE
		}

		if t >= nextPrint {
			fmt.Printf("%-8.3f  %-12.1f  %-12.1f  %-12.2e\n",
				t/3600, r/1e3, v, errE)
			nextPrint += T / 4 // four samples per orbit
		}

		t, SBI, VBI = integrator.Step(dt, md3.Vec{})
		totalSteps++
		if totalSteps == 1 {
			fmt.Printf("gnco Step %d: %+v\n\n", totalSteps, integrator.RK())
		}
	}

	fmt.Println()
	fmt.Printf("Completed %.0f orbits in %d steps (dt = %.0f s).\n",
		float64(nOrbits), totalSteps, dt)
	fmt.Printf("Maximum specific energy error |ΔE/E₀|: %.2e\n", maxErrE)
	return nil
}
