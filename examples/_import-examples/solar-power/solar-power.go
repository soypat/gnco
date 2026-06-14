package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"time"

	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

var (
	// Simulation start date: 12 Nov 2026 21:36:00.000
	startDate = time.Date(2026, time.November, 12, 21, 36, 0, 0, time.UTC)
)

const (
	day      = 24 * 60 * 60
	propDT   = 60.0        // [s] coarse propagation sampling.
	propSpan = 366.0 * day // propagation span
	// Define Orbital elements.
	deg          = math.Pi / 180
	kepOrbitSMA  = 6928.5e3 // [m]
	kepOrbitECC  = 0.0011013
	kepOrbitINC  = 97.794 * deg
	kepOrbitRAAN = 215.84 * deg
	kepOrbitAOP  = 191.1 * deg
	kepOrbitTA   = 118.33 * deg
)

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

func run() error {
	earth := cosmos.NewEarth()
	sun := cosmos.NewAnalyticSun()
	epoch0 := cosmos.EpochFromTime(startDate)
	k, err := orbits.NewKeplerian(kepOrbitSMA, kepOrbitECC, kepOrbitINC, kepOrbitRAAN, kepOrbitAOP)
	if err != nil {
		return err
	}
	jgm2, err := cosmos.JGM2(4, 4)
	if err != nil {
		return err
	}
	fm := gnco.NewForceModel(earth)
	fm.SetHarmonics(jgm2)
	fm.SetNutationInterval(120) // GMAT Nutation Update Interval default. Main optimization.
	r0, v0 := k.RV(earth.Mu(), kepOrbitTA)
	prop, err := gnco.NewOrbitPropagator(fm, epoch0, r0, v0, gnco.PropagatorConfig{
		Accuracy: 1e-12, MinStep: 0.001, MaxStep: 2700, InitialStep: 60,
	})
	if err != nil {
		return err
	}
	// Propagate the whole span once; attitude is irrelevant to eclipse geometry
	// so the identity attitude (nil) is used.
	traj, err := gnco.Propagate(prop, propDT, propSpan, nil)
	if err != nil {
		return err
	}

	// Beta angle range, sampled hourly.
	i := 0
	betaMin, betaMax := math.Inf(1), math.Inf(-1)
	for t := 0.0; t < traj.Span(); t += 60 * 60 {
		for traj.Samples[i].T.Sub(traj.Samples[0].T) < t && i < traj.Len()-1 {
			i++
		}
		beta := traj.Samples[i].OrbitBetaAngle(sun) / deg
		betaMin, betaMax = math.Min(betaMin, beta), math.Max(betaMax, beta)
	}
	fmt.Printf("%.1f days, beta angle range [%.2f°, %.2f°]\n", traj.Span()/3600/24, betaMin, betaMax)
	// Eclipse phases (penumbra/umbra), boundaries refined between samples
	// and write the phases in GMAT EclipseLocator report format.
	var buf bytes.Buffer
	ecls := traj.Eclipses(sun, earth.Radius(), earth.Flattening(), sun.Radius())
	writeEclipseReport(&buf, "Sat", "Earth", ecls)
	os.WriteFile("eclipses.txt", buf.Bytes(), 0777)
	fmt.Println("wrote eclipses to file")

	return nil
}

// writeEclipseReport writes the eclipse phases in GMAT EclipseLocator report
// format. Contiguous phases of one shadow pass (each Enter equal to the previous
// phase's Exit) share an Event Number and a Total Duration, matching GMAT.
func writeEclipseReport(w io.Writer, spacecraft, occBody string, ecl []cosmos.Eclipse) {
	const layout = "02 Jan 2006 15:04:05.000"
	fmt.Fprintf(w, "Spacecraft: %s\n\n", spacecraft)
	fmt.Fprintf(w, "%-28s%-30s%-16s%-16s%-12s%-14s%s\n",
		"Start Time (UTC)", "Stop Time (UTC)", "Duration (s)",
		"Occ Body", "Type", "Event Number", "Total Duration (s)")
	for i, ev := 0, 1; i < len(ecl); ev++ {
		// Gather the contiguous phases of this pass.
		j := i + 1
		for j < len(ecl) && ecl[j].Enter.Sub(ecl[j-1].Exit) == 0 {
			j++
		}
		total := ecl[j-1].Exit.Sub(ecl[i].Enter)
		for k := i; k < j; k++ {
			e := ecl[k]
			fmt.Fprintf(w, "%-28s%-30s%-16.11g%-16s%-12s%-14d%.11g\n",
				e.Enter.Time().Format(layout), e.Exit.Time().Format(layout),
				e.Duration(), occBody, e.Kind.String(), ev, total)
		}
		i = j
	}
}
