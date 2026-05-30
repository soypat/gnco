package main

import (
	"fmt"
	"log"
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
)

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

func run() error {
	const (
		// Vehicle — single-stage sounding rocket.
		massWet     = 1800.0   // [kg] initial mass (structure + propellant)
		massDry     = 1000.0   // [kg] structural mass after burnout
		thrustForce = 40_000.0 // [N]
		burnTime    = 40.0     // [s]  →  mass flow = 20 kg/s
		refArea     = 0.07     // [m²] reference cross-section (~30 cm Ø)
		dragCoeff   = 0.35     // [-]
		liftCoeff   = 0.0      // [-]
		// Launch geometry — 85° elevation, heading North.
		launchElevDeg = 85.0
		launchBearing = 0.0 // [rad] bearing from North

		// rampLength matches trajectory-sim's dof5 convention: TVG is held fixed at
		// the launch orientation until the rocket has traveled this distance from the
		// launch pad. This prevents the severe gravity turn that occurs when the
		// velocity-aligned thrust model is applied at near-zero ground-launch speed.

		// ramp needed to maintain initial bearing, else gravity will pull nose down
		// while rocket still slow and vulnerable to gravity's impact.
		// To simulate the ramp we keep TVG (direction tensor) constant until ramp is cleared.
		rampLength = 20.0 // [m] launch rail length

		// Integration.
		dt   = 0.5    // [s] fixed step
		maxT = 1200.0 // [s] abort cutoff

		T0sea = 288.15 // [K] ISA sea-level temperature
	)

	earth := gnco.NewEarth()
	launchSite := earth.GeocentricFromDegrees(-34.6, -58.4, earth.HASLToElevation(25))

	SBI0, TGI0 := launchSite.InertialCoords(0)

	// Tiny seed speed so TVGFromGeographicVelocity returns the correct launch-
	// direction orientation. Bearing is North, elevation 85°.
	launchElevRad := launchElevDeg * math.Pi / 180
	VBG0 := gnco.GeographicVectorFromElevationAndBearing(launchElevRad, launchBearing, 0.01)

	// Earth rotation vector (ECI z-axis).
	weii := md3.Vec{Z: earth.Rotation()}

	// Initial inertial velocity: geographic velocity rotated to inertial + Earth surface velocity.
	// VBI = TGI^T * VBG + (ω × SBI)
	VBI0 := md3.Add(md3.Cross(weii, SBI0), md3.MulMatVecTrans(TGI0, VBG0))

	// coords tracks the rocket position; the integrator updates it at every RKN stage.
	coords := launchSite
	integrator := gnco.NewPhysicsPointIntegrator(&coords, 0, SBI0, VBI0)

	mass := massWet
	massFlow := (massWet - massDry) / burnTime

	t := 0.0
	SBI := SBI0
	VBI := VBI0

	// TVG is held at the initial launch orientation until the rocket clears the ramp.
	initialTVG := gnco.TVGFromGeographicVelocity(VBG0)
	TVG := initialTVG
	var sbeg md3.Vec // geographic displacement from launch pad [m]

	var apogeeHASL, apogeeTime float64
	nextPrint := 0.0

	fmt.Printf("%-8s  %-11s  %-11s  %-10s  %-6s  %-6s\n",
		"t [s]", "HASL [km]", "speed [m/s]", "mass [kg]", "Mach", "Q [Pa]")
	fmt.Println("--------  -----------  -----------  ----------  ------  ------")

	for t < maxT {
		hasl := coords.HASL()
		TGI := coords.TGI(t)

		// Earth-relative velocity in geographic frame: VBEG = TGI * (VBI - ω×SBI).
		// Aero forces and TVG updates use Earth-relative speed, not inertial speed.
		VBEG := md3.MulMatVec(TGI, md3.Sub(VBI, md3.Cross(weii, SBI)))
		speed := md3.Norm(VBEG)

		if hasl > apogeeHASL {
			apogeeHASL = hasl
			apogeeTime = t
		}

		Tatm, _, rho := gnco.InternationalStandardAtmosphere(hasl, T0sea)
		Fdrag, Flift, Q, mach := AeroForces(speed, rho, Tatm, dragCoeff, liftCoeff, refArea)
		if t >= nextPrint {
			fmt.Printf("%-8.1f  %-11.3f  %-11.1f  %-10.1f  %-6.2f  %-6.0f\n",
				t, hasl/1000, speed, mass, mach, Q)
			nextPrint += 20.0
		}

		// Stop once descending back to launch altitude (after a grace period post-burnout).
		if t > burnTime+10 && hasl < launchSite.HASL() {
			break
		}

		// Accelerations in the velocity frame.
		// Gravity is handled internally by PhysicsPointIntegrator; do not add it here.
		// See [gnco.Orientation] and [gnco.Frame] to understand velocity frame orientation.
		var FpsV = md3.Vec{X: -Fdrag, Z: Flift}
		AFpsV := md3.Scale(1/mass, FpsV) // F=ma
		if t < burnTime {
			AFpsV.X += thrustForce / mass
			mass = max(massDry, mass-massFlow*dt)
		}

		accelGeog := md3.MulMatVecTrans(TVG, AFpsV)
		t, SBI, VBI = integrator.Step(dt, accelGeog)
		// Post-step: integrate geographic displacement and update TVG.
		newTGI := coords.TGI(t)
		newVBEG := md3.MulMatVec(newTGI, md3.Sub(VBI, md3.Cross(weii, SBI)))
		sbeg = md3.Add(sbeg, md3.Scale(dt/2, md3.Add(newVBEG, VBEG)))
		if t > 10 || md3.Norm(sbeg) > rampLength {
			// If ramp cleared we update direction.
			// The time comparison avoids expensive Norm calculation by which ramp should definetely be cleared.
			TVG = gnco.TVGFromGeographicVelocity(newVBEG)
		}
	}

	fmt.Println()
	fmt.Printf("Apogee: %.1f km  at t = %.0f s\n", apogeeHASL/1000, apogeeTime)
	return nil
}

// AeroForces computes aerodynamic drag and lift forces, dynamic pressure, and
// Mach number for a body moving at vel [m/s] through air of density rho [kg/m³]
// and temperature T [K], given drag coefficient Cd, lift coefficient Cl, and
// reference area Aref [m²].
//
//	Q     = ½·ρ·v²
//	Mach  = v / sqrt(γ·R/M·T)
//	Fdrag = Cd · Q · Aref
//	Flift = Cl · Q · Aref
func AeroForces(vel, rho, T, Cd, Cl, Aref float64) (Fdrag, Flift, Q, Mach float64) {
	const (
		gamma  = 1.4      // specific heat ratio for air
		airMol = 28.97e-3 // molar mass of air [kg/mol]
		Ru     = 8.314472 // universal gas constant [J/(mol·K)]
	)
	Q = 0.5 * rho * vel * vel
	Mach = vel / math.Sqrt(gamma*Ru/airMol*T)
	Fdrag = Cd * Q * Aref
	Flift = Cl * Q * Aref
	return
}
