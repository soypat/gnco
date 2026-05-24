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