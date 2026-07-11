package main

import (
	"cmp"
	"fmt"
	"log"
	"math"
	"slices"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/physics"
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

	earth := cosmos.NewEarth()
	launchSite := gnco.NewGeocentricFromDegrees(earth, -34.6, -58.4, earth.HASLToElevation(25))

	SBI0, TGI0 := launchSite.InertialCoords(cosmos.EpochFromTT(0))

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
	epoch0 := cosmos.Epoch{} // J2000 anchor; only elapsed time matters here.
	var integrator physics.PointIntegrator
	err := integrator.ConfigureCoord(&coords, epoch0, SBI0, VBI0)
	if err != nil {
		return err
	}
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
		TGI := coords.TGI(cosmos.EpochFromTT(t))

		if hasl > apogeeHASL {
			apogeeHASL = hasl
			apogeeTime = t
		}

		Tatm, _, rho := gnco.InternationalStandardAtmosphere(hasl, T0sea)

		// Earth-relative velocity in geographic frame: VBEG = TGI * (VBI - ω×SBI).
		// Aero forces and TVG updates use Earth-relative speed, not inertial speed.
		VBEG := md3.MulMatVec(TGI, md3.Sub(VBI, md3.Cross(weii, SBI)))
		// Zero angle-of-attack assumption: body axis aligned with Earth-relative
		// velocity, so body frame == velocity frame. For an axisymmetric rocket at
		// zero AoA aero reduces to pure drag along velocity (no lift by symmetry).
		// Only approximate on the ramp, where TVG is held at launch orientation.
		speed := md3.Norm(VBEG)
		Fdrag, Q, mach := AeroForces(speed, rho, Tatm, hasl, refArea)
		Flift := 0.0
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
		var e cosmos.Epoch
		e, SBI, VBI = integrator.Step(dt, accelGeog)
		t = e.Sub(epoch0)
		// Post-step: integrate geographic displacement and update TVG.
		newTGI := coords.TGI(e)
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
	if dragExtrapolatedHeight > 0 {
		// If AeroForces extrapolates outside drag coefficient table then this will print.
		fmt.Printf("Drag extrapolated @ %.1fkm\n", dragExtrapolatedHeight/1e3)
	}
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
func AeroForces(vel, rho, T, hasl, Aref float64) (Fdrag, Q, Mach float64) {
	const (
		gamma  = 1.4      // specific heat ratio for air
		airMol = 28.97e-3 // molar mass of air [kg/mol]
		Ru     = 8.314472 // universal gas constant [J/(mol·K)]
	)
	Q = 0.5 * rho * vel * vel
	Mach = vel / math.Sqrt(gamma*Ru/airMol*T)
	var Cd float64
	if hasl < 60e3 { // If in space drag=0.
		Cd = dragCoeff.Interp1(hasl, Mach)
		if math.IsNaN(Cd) {
			if dragExtrapolatedHeight < 0 {
				dragExtrapolatedHeight = hasl
			}
			// height or Mach out of table bounds, do a extrapolation.
			Cd = dragCoeff.Extrap1(hasl, Mach)
		}
	}
	Fdrag = Cd * Q * Aref
	return
}

var (
	dragExtrapolatedHeight float64 = -100
	dragCoeff                      = dualEntryTable{
		X: []float64{0, 10e3, 20e3, 30e3, 40e3, 50e3}, // altitude [m]
		// Records Y is mach number.
		Records: []dualEntryRecord{
			{Y: 0., F: []float64{0.500, 0.500, 0.500, 0.600, 0.800, 1.300}},
			{Y: .5, F: []float64{0.404, 0.449, 0.542, 0.687, 0.910, 1.231}},
			{Y: .8, F: []float64{0.383, 0.424, 0.508, 0.636, 0.833, 1.111}},
			{Y: .9, F: []float64{0.396, 0.437, 0.521, 0.649, 0.843, 1.118}},
			{Y: 1.1, F: []float64{0.505, 0.541, 0.614, 0.726, 0.897, 1.137}},
			{Y: 1.2, F: []float64{0.497, 0.530, 0.597, 0.701, 0.858, 1.080}},
			{Y: 1.5, F: []float64{0.480, 0.510, 0.572, 0.668, 0.813, 1.016}},
			{Y: 2, F: []float64{0.414, 0.441, 0.496, 0.581, 0.690, 0.893}},
			{Y: 3, F: []float64{0.303, 0.325, 0.370, 0.440, 0.548, 0.703}},
			{Y: 4, F: []float64{0.235, 0.253, 0.290, 0.350, 0.440, 0.581}},
		},
	}
)

// TODO: consider creating a package for dual entry tables?

type dualEntryRecord struct {
	Y float64
	F []float64
}

type dualEntryTable struct {
	X       []float64
	Records []dualEntryRecord
}

// Interp1 provides a bilinear interpolation for a dual entry table.
func (d *dualEntryTable) Interp1(x, y float64) (f float64) {
	i, oobx := d.lowerXIdx(x)
	j, ooby := d.lowerYIdx(y)
	if oobx || ooby { // We do not extrapolate in this routine.
		return math.NaN()
	}
	return d.interpRaw1(i, j, x, y)
}

// Extrap1 provides a bilinear extrapolation for a dual entry table.
// Extrap1 expects to extrapolate and will panic if both x and y are inside
// table domain.
func (d *dualEntryTable) Extrap1(x, y float64) (f float64) {
	i, oobx := d.lowerXIdx(x)
	j, ooby := d.lowerYIdx(y)
	if !oobx && !ooby {
		panic("Extrap1 expected values outside domain")
	}
	return d.interpRaw1(i, j, x, y)
}

// interpRaw1 does raw bilinear interpolation or extrapolation depending on x,y values.
func (d *dualEntryTable) interpRaw1(i, j int, x, y float64) (f float64) {
	f11 := d.Records[j].F[i]
	f12 := d.Records[j+1].F[i]
	f21 := d.Records[j].F[i+1]
	f22 := d.Records[j+1].F[i+1]
	x1, x2 := d.X[i], d.X[i+1]
	y1, y2 := d.Records[j].Y, d.Records[j+1].Y
	x2mx := (x2 - x) / (x2 - x1)
	xmx1 := (x - x1) / (x2 - x1)
	return (y2-y)/(y2-y1)*(x2mx*f11+xmx1*f21) + (y-y1)/(y2-y1)*(x2mx*f12+xmx1*f22)
}

// lowerXIdx returns the index of the lower x value which contains x0.
//
//	x[idx] <= x0 < x[idx+1]
//
// idx will always be less than len(x)-1. If x0 is outside table x bounds
// oob is true and idx is the nearest edge interval, ready for extrapolation.
func (d *dualEntryTable) lowerXIdx(x0 float64) (xidx int, oob bool) {
	switch {
	case x0 < d.X[0]:
		return 0, true
	case x0 > d.X[len(d.X)-1]:
		return len(d.X) - 2, true
	}
	i, _ := slices.BinarySearch(d.X, x0)
	return max(i-1, 0), false
}

// lowerYIdx returns the index of the lower y value which contains y0.
//
//	y[idx] <= y0 < y[idx+1]
//
// idx will always be less than len(y)-1. If y0 is outside table y bounds
// oob is true and idx is the nearest edge interval, ready for extrapolation.
func (d *dualEntryTable) lowerYIdx(y0 float64) (yidx int, oob bool) {
	switch {
	case y0 < d.Records[0].Y:
		return 0, true
	case y0 > d.Records[len(d.Records)-1].Y:
		return len(d.Records) - 2, true
	}
	j, _ := slices.BinarySearchFunc(d.Records, y0, func(r dualEntryRecord, y float64) int {
		return cmp.Compare(r.Y, y)
	})
	return max(j-1, 0), false
}
