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
go test ./internal/ode -bench=. -benchmem
goos: linux
goarch: amd64
pkg: github.com/soypat/gnco/internal/ode
cpu: 12th Gen Intel(R) Core(TM) i5-12400F
BenchmarkRK45Step-12                       82034             14125 ns/op               0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/RK4(5)-12   15809316                73.90 ns/op            0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/RKF7(8)-12    6944835               171.5 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_noadaptivestep/RKN12(10)-12                 4857708               242.5 ns/op             0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RK4(5)-12                 83965             14250 ns/op                 2.581 errY×1e-9                93.00 steps/op         0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RKF7(8)-12                270642              4290 ns/op                 0.6264 errY×1e-9               16.00 steps/op         0 B/op          0 allocs/op
BenchmarkIVP_adaptiveconvergence/RKN12(10)-12             254018              4570 ns/op                 0.0000001 errY×1e-9            13.00 steps/op         0 B/op          0 allocs/op
PASS
ok      github.com/soypat/gnco/internal/ode     8.220s
```