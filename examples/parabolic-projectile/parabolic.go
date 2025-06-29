package main

import (
	"fmt"
	"log"
	"math"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
)

func main() {
	err := run()
	if err != nil {
		log.Fatal(err)
	}
}

func run() error {
	// The "World" type provides fixed-frame facilities as well as
	// simple or geodesic gravity calculation.
	earth := gnco.NewEarth()
	buenosAires := earth.GeocentricFromDegrees(34.6, 58.4, earth.HASLToElevation(25))
	// We declare our initial conditions for the integrator.
	// Note we integrate in inertial coordinates to avoid ficticious forces.
	const (
		t0              = 0.0 // [s] epoch time
		initialVelocity = 60  // [m/s]
		launchAngleElev = 45  // [degrees]

		projectileAngleRad = launchAngleElev * math.Pi / 180 // [rad]
	)
	SBI0, TGI := buenosAires.InertialCoords(t0)
	// Calculate the velocity in inertial frame of reference.
	// We start out with velocity in geographical frame since declaring it
	// with an elevation angle relative to our horizon makes it easier to reason about.
	VBG0 := gnco.GeographicVectorFromElevationAndBearing(projectileAngleRad, 0, initialVelocity)
	orient0 := gnco.Orientation{
		// Disregard body and velocity transformations. We are interested in a single point.
		TBV: md3.IdentityMat3(),
		TVG: md3.IdentityMat3(),
		TGI: TGI,
	}
	VBI0 := gnco.FrameGeographic.ToInertial(orient0, VBG0)

	// Instantiate physics engine. The physics engine integrates for a point.
	// For the simplicity of the example there is no external/internal force other than gravity
	// so we can omit force/mass calculations.
	projectileCoords := buenosAires // projectileCoords will store coordinates of our projectile over course of simulation.
	integrator := gnco.NewPhysicsPointIntegrator(&projectileCoords, t0, SBI0, VBI0)
	dt := 0.0001
	t := t0
	wantTime := parabolicTimeOfFlight(initialVelocity, projectileAngleRad, md3.Norm(buenosAires.AGravG()))
	var SBI, VBI md3.Vec
	for projectileCoords.Elev >= buenosAires.Elev && t-t0 < wantTime*2 {
		// No internal acceleration other than coordinate system gravity.
		// We should get parabolic trajectory.
		accelGeographical := md3.Vec{X: 0, Y: 0, Z: 0}
		t, SBI, VBI = integrator.Step(dt, accelGeographical)
	}
	simulationDuration := t - t0
	fmt.Println("total flight duration", simulationDuration, "with final velocity", md3.Norm(VBI))
	fmt.Println("expected flight duration", wantTime, "error +/-", dt)
	// we can also compare distance with typical parabolic trajectory distance.
	distance := md3.Norm(md3.Sub(SBI, SBI0))
	wantDistance := parabolicDistanceOfFlight(initialVelocity, projectileAngleRad, md3.Norm(buenosAires.AGravG()))
	fmt.Println("distance from launch:", distance)
	fmt.Println("expected distance:", wantDistance)
	return nil
}

func parabolicTimeOfFlight(v0, elevation, g float64) float64 {
	return 2 * v0 * math.Sin(elevation) / g
}

func parabolicDistanceOfFlight(v0, elevation, g float64) float64 {
	return v0 * v0 * math.Sin(2*elevation) / g
}
