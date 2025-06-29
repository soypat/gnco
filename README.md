# gnco
[![go.dev reference](https://pkg.go.dev/badge/github.com/soypat/gnco)](https://pkg.go.dev/github.com/soypat/gnco)
[![Go Report Card](https://goreportcard.com/badge/github.com/soypat/gnco)](https://goreportcard.com/report/github.com/soypat/gnco)
[![codecov](https://codecov.io/gh/soypat/gnco/branch/main/graph/badge.svg)](https://codecov.io/gh/soypat/gnco)
[![Go](https://github.com/soypat/gnco/actions/workflows/go.yml/badge.svg)](https://github.com/soypat/gnco/actions/workflows/go.yml)
[![sourcegraph](https://sourcegraph.com/github.com/soypat/gnco/-/badge.svg)](https://sourcegraph.com/github.com/soypat/gnco?badge)

gnco provides logic for projectile trajectory calculation. See [`parabolic-projectile`](./examples/parabolic-projectile/parabolic.go)
for a basic example of use to recreate a parabolic trajectory of a point mass with no external forces.



## Install
How to install package with newer versions of Go (+1.16):
```sh
go mod download github.com/soypat/gnco@latest
```

### About the integrator
The physics integrator used is a state of the art Runge-Kutta-Nyström 12(10) second order integrator and presents very well behaved energy conservation for
elliptical orbits for very large integration steps in the order of the hundreds of seconds, given no external forces other than gravity are acting.