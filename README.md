# gnco
[![go.dev reference](https://pkg.go.dev/badge/github.com/soypat/gnco)](https://pkg.go.dev/github.com/soypat/gnco)
[![Go Report Card](https://goreportcard.com/badge/github.com/soypat/gnco)](https://goreportcard.com/report/github.com/soypat/gnco)
[![codecov](https://codecov.io/gh/soypat/gnco/branch/main/graph/badge.svg)](https://codecov.io/gh/soypat/gnco)
[![Go](https://github.com/soypat/gnco/actions/workflows/go.yml/badge.svg)](https://github.com/soypat/gnco/actions/workflows/go.yml)
[![sourcegraph](https://sourcegraph.com/github.com/soypat/gnco/-/badge.svg)](https://sourcegraph.com/github.com/soypat/gnco?badge)

gnco provides logic for projectile trajectory calculation on a rotating Earth model.

## Examples

- [`parabolic-projectile`](./examples/parabolic-projectile/parabolic.go) — minimal example: parabolic trajectory of a point mass with no external forces, demonstrating the integrator and coordinate system basics.

- [`5dof-rocket`](./examples/5dof-rocket/main.go) — 5-DOF sounding rocket simulation. A single-stage solid-fuel rocket launches from a geographic site at 85° elevation with Earth rotation, ISA atmosphere, drag, and a launch-tower hold to prevent a gravity turn at low speed.

- [`keplerian-orbit`](./examples/keplerian-orbit/main.go) — energy-conservation demo: an elliptical orbit propagated with the RKN12(10) integrator at large 300 s steps (~19 steps/orbit), tracking the specific orbital energy `|ΔE/E₀|` to show how well the integrator preserves this invariant.

- [`sun-synchronous-orbit`](./examples/sun-synchronous-orbit/main.go) — sun-synchronous orbit (550 km, 97.8°) propagated with the JGM-2 4×4 spherical-harmonic gravity field and IAU-76/FK5 Earth orientation. Force model validated against GMAT (~1.2 m/day).

## Install
How to install package with newer versions of Go (+1.16):
```sh
go mod download github.com/soypat/gnco@latest
```

### About the integrator
The physics integrator used is a state of the art Runge-Kutta-Nyström 12(10) second order integrator and presents very well behaved energy conservation for
elliptical orbits for very large integration steps in the order of the hundreds of seconds, given no external forces other than gravity are acting.

Below are some oscillator benchmarks for integrators in this project. Notably RKN12(10) has near ULP precision and is far more precise for adaptive stepping compared with RKF7(8) and RK4(5) and converges in about the same amount of wall-clock time as RK7(8) and much faster than RK4(5).
```
go test ./physics/ode -bench=. -benchmem
goos: linux
goarch: amd64
pkg: github.com/soypat/gnco/physics/ode
cpu: 12th Gen Intel(R) Core(TM) i5-12400F
BenchmarkIVP_noadaptivestep/RK4(5)-12   16229380                70.98 ns/op            0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/RKF7(8)-12   7191062               173.0 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/Verner9-12   5540967               226.4 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/Feagin12-12                  2364345               498.0 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/RKN12(10)-12                 6000030               192.3 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RK4(5)-12                 83816             13752 ns/op                 2.581 errY×1e-9     93.00 steps/op            0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RKF7(8)-12               267368              4120 ns/op                 0.6264 errY×1e-9     16.00 steps/op            0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/Verner9-12               293768              3795 ns/op                 0.09069 errY×1e-9     13.00 steps/op            0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/Feagin12-12              155316              7209 ns/op                 0.05304 errY×1e-9     13.00 steps/op            0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RKN12(10)-12             323271              3561 ns/op                 0.0000001 errY×1e-9     13.00 steps/op            0 B/op          0 allocs/op
PASS
ok      github.com/soypat/gnco/physics/ode      11.628s
```