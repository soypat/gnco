// Solar-power example: propagate a nadir-pointing CubeSat, locate the worst-case
// orbit (the one period of least generated power over the year), and reproduce
// GMAT SolarCalc's per-face power/energy for that period (Power_Calc reference).
//
// Power model (from SolarCalc's potencia.py), per sample and per face f:
//
//	cos_f   = max(0, n̂_f · ŝ)     n̂_f body-frame face normal rotated to inertial
//	                              ŝ   satellite→Sun unit direction
//	cos_f   = 0 while the satellite is in Earth shadow (penumbra or umbra)
//	P_f     = Pmp · cos_f         [W per cell]   (irradiance baked into Pmp)
//	P_total = Σ_f  cells_f · P_f
//
// Pmp is one scalar per case (Lab/BOL), derived from the cell electricals at a
// given temperature with degradation factors. Nothing is stored on the
// trajectory: power is derived from the raw state (R, V, attitude) on demand.
package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"time"

	"github.com/soypat/geometry/md3"
	"github.com/soypat/gnco"
	"github.com/soypat/gnco/cosmos"
	"github.com/soypat/gnco/orbits"
)

var (
	// Simulation start date: 12 Nov 2026 21:36:00.000 (SolarCalc worst-case set).
	startDate = time.Date(2026, time.November, 12, 21, 36, 0, 0, time.UTC)
	earth     = cosmos.NewEarth()
	sun       = cosmos.NewAnalyticSun()
)

const (
	day      = 24 * 60 * 60
	propDT   = 60.0         // [s] coarse propagation/scan step.
	propSpan = 366.0 * day  // propagation span.
	fineDT   = 0.5e-5 * day // [s] fine step over the worst-case period (GMAT stepWC = 0.432 s).
	// Orbital elements (worst-case set).
	deg          = math.Pi / 180
	kepOrbitSMA  = 6928.5e3 // [m]
	kepOrbitECC  = 0.0011013
	kepOrbitINC  = 97.794 * deg
	kepOrbitRAAN = 215.84 * deg
	kepOrbitAOP  = 191.1 * deg
	kepOrbitTA   = 118.33 * deg

	// Shadow geometry matching GMAT / milestone 05 (JGM2 equatorial radius, Sun radius).
	earthShadowR = 6378136.3 // [m]
	sunShadowR   = 695990e3  // [m]

	// Solar cell electricals (settings.ini [celda solar]; densities × area → A, V).
	cellAreaCm2 = 30.18   // [cm²]
	cellVmp     = 2.409   // [V] max-power-point voltage at reference temperature
	cellJmp     = 0.01666 // [A/cm²] max-power-point current density
	cellTcVmp   = -0.0067 // [V/°C] voltage temperature coefficient
	cellTcImp   = 8e-6    // [A/°C] current temperature coefficient
	cellTref    = 28.0    // [°C] reference (laboratory) temperature

	// Case temperatures. BOL uses the on-orbit cell temperature: 55 °C reproduces
	// the committed Power_Calc output (the run folder's stale settings T=60 °C
	// would not — it yields Pmp_BOL≈0.997 vs the file's ≈1.12 per cell).
	labTempC = 28.0
	bolTempC = 55.0
)

// Cube faces in SolarCalc order: +X, +Y, +Z, -X, -Y, -Z.
var (
	faceNames    = [6]string{"+X", "+Y", "+Z", "-X", "-Y", "-Z"}
	faceNormals  = [6]md3.Vec{{X: 1}, {Y: 1}, {Z: 1}, {X: -1}, {Y: -1}, {Z: -1}}
	cellsPerFace = [6]int{2, 2, 0, 2, 1, 2}
)

// Reference per-cell face energies [J] and mean power [W] from GMAT Power_Calc
// (periodoWC report), used to validate the near-epoch worst case.
var (
	refLab = [6]float64{1538.69, 673.06, 189.95, 1544.15, 0, 2190.22}
	refBOL = [6]float64{1423.28, 622.58, 175.71, 1428.34, 0, 2025.95}
)

const (
	refMeanLab = 2.07
	refMeanBOL = 1.92
)

// attNadir is the GMAT NadirPointing attitude: body +Z to nadir, +X along velocity.
func attNadir(_ cosmos.Epoch, r, v md3.Vec) md3.Quat {
	return gnco.State{R: r, V: v}.AttNadirPointing()
}

// SolarCell holds a cell's max-power-point electricals and temperature
// coefficients. Pmp (per-cell power) is derived, never stored.
type SolarCell struct {
	Imp, Vmp     float64 // max-power-point current [A] / voltage [V] at Tref
	TcImp, TcVmp float64 // temperature coefficients [A/°C], [V/°C]
	Tref         float64 // reference temperature [°C]
}

// Pmp returns the per-cell maximum power [W] at cell temperature tempC, after
// applying the multiplicative current/voltage degradation factors degI, degV
// (1,1 = none). Mirrors SolarCalc add_deg followed by add_degT.
func (c SolarCell) Pmp(tempC, degI, degV float64) float64 {
	dT := tempC - c.Tref
	imp := c.Imp*degI + c.TcImp*dT
	vmp := c.Vmp*degV + c.TcVmp*dT
	if p := imp * vmp; p > 0 {
		return p
	}
	return 0
}

// Satellite is the panel layout: a shared cell type, the body-frame face normals,
// and the cell count on each face (parallel slices, no per-face redundancy).
type Satellite struct {
	Cell    SolarCell
	Normals []md3.Vec
	Cells   []int
}

// facePowerPerCell fills dst[f] with the per-cell power [W] of each face at state
// s for per-cell max power pmp: Pmp·max(0, n̂·ŝ), or 0 for every face in shadow.
func (sat Satellite) facePowerPerCell(s gnco.State, sun cosmos.Ephemeris, pmp float64, dst []float64) {
	if inShadow(s, sun) {
		for f := range dst {
			dst[f] = 0
		}
		return
	}
	sd := sunDir(s, sun)
	for f, n := range sat.Normals {
		c := md3.Dot(s.Att.Rotate(n), sd) // n̂ rotated body→inertial, dotted with Sun direction
		if c < 0 {
			c = 0 // back-facing
		}
		dst[f] = pmp * c
	}
}

// totalPower sums per-cell face powers weighted by the cell count on each face.
func (sat Satellite) totalPower(perCell []float64) (total float64) {
	for f, p := range perCell {
		total += float64(sat.Cells[f]) * p
	}
	return total
}

// sunDir returns the unit satellite→Sun direction in the inertial frame.
func sunDir(s gnco.State, sun cosmos.Ephemeris) md3.Vec {
	return md3.Unit(md3.Sub(sun.Position(s.T), s.R))
}

// inShadow reports whether the satellite is within Earth's penumbra or umbra.
func inShadow(s gnco.State, sun cosmos.Ephemeris) bool {
	pen, _ := cosmos.Shadow(sun.Position(s.T), s.R, sunShadowR, earthShadowR, earth.Flattening())
	return pen < 0
}

func main() {
	if err := run(); err != nil {
		log.Fatal(err)
	}
}

func run() error {
	epoch0 := cosmos.EpochFromTime(startDate)
	mu := earth.Mu()
	k, err := orbits.NewKeplerian(kepOrbitSMA, kepOrbitECC, kepOrbitINC, kepOrbitRAAN, kepOrbitAOP)
	if err != nil {
		return err
	}
	jgm2, err := cosmos.JGM2(4, 4)
	if err != nil {
		return err
	}
	cfg := gnco.PropagatorConfig{Accuracy: 1e-12, MinStep: 0.001, MaxStep: 2700, InitialStep: 60}
	r0, v0 := k.RV(mu, kepOrbitTA)

	// Coarse propagation over the whole span (used to scan for the worst case).
	prop, err := newProp(earth, jgm2, epoch0, r0, v0, cfg)
	if err != nil {
		return err
	}
	traj, err := gnco.Propagate(prop, propDT, propSpan, attNadir)
	if err != nil {
		return err
	}
	fmt.Printf("%.1f days propagated, %d coarse samples\n", traj.Span()/day, traj.Len())

	sat := Satellite{
		Cell:    SolarCell{Imp: cellJmp * cellAreaCm2, Vmp: cellVmp, TcImp: cellTcImp, TcVmp: cellTcVmp, Tref: cellTref},
		Normals: faceNormals[:],
		Cells:   cellsPerFace[:],
	}
	pmpLab := sat.Cell.Pmp(labTempC, 1, 1)
	pmpBOL := sat.Cell.Pmp(bolTempC, 1, 1)
	fmt.Printf("Pmp per cell: Lab(%.0f°C)=%.4f W   BOL(%.0f°C)=%.4f W\n", labTempC, pmpLab, bolTempC, pmpBOL)

	// Two worst-case definitions over one orbital period:
	//   - near-epoch: least power within GMAT's WC sub-simulation (~8.5 orbits
	//     after epoch); reproduces GMAT's periodoWC report.
	//   - global: least power anywhere in the year — the true mission worst case
	//     (beta drifts toward 0 over the year, so this is months out and lower).
	period := traj.Samples[0].OrbitPeriod(mu)
	const wcSubSimSpan = 0.5664573915369207 * day // GMAT WC sub-sim duration

	nearIdx := worstCasePeriod(sat, traj, sun, pmpLab, period, epoch0, wcSubSimSpan)
	near, err := analyzePeriod(traj.Samples[nearIdx], period, earth, jgm2, cfg, sat, sun, pmpLab, pmpBOL)
	if err != nil {
		return err
	}
	printCase("near-epoch worst (matches GMAT periodoWC)", near, sat, refLab[:], refBOL[:])

	globalIdx := worstCasePeriod(sat, traj, sun, pmpLab, period, epoch0, 0)
	global, err := analyzePeriod(traj.Samples[globalIdx], period, earth, jgm2, cfg, sat, sun, pmpLab, pmpBOL)
	if err != nil {
		return err
	}
	printCase("global worst over the year (true mission worst case)", global, sat, nil, nil)

	// Per-sample power series (Power_Calc-like table) for both cases.
	if err := writeSeriesFile("power-near-epoch.txt", near.traj, sun, sat, pmpLab, pmpBOL); err != nil {
		return err
	}
	if err := writeSeriesFile("power-global.txt", global.traj, sun, sat, pmpLab, pmpBOL); err != nil {
		return err
	}
	fmt.Println("\nwrote per-sample series to power-near-epoch.txt and power-global.txt")
	return nil
}

// caseResult holds the analysis of one worst-case orbital period.
type caseResult struct {
	span             float64
	eLab, eBOL       []float64 // per-cell face energies [J]
	meanLab, meanBOL float64   // orbital mean power [W]
	traj             *gnco.Trajectory
}

// analyzePeriod re-propagates one orbital period from start state ws at the fine
// step and computes per-face Lab/BOL energies and mean power over it.
func analyzePeriod(ws gnco.State, period float64, earth *cosmos.Body, jgm2 *cosmos.Harmonics, cfg gnco.PropagatorConfig, sat Satellite, sun cosmos.Ephemeris, pmpLab, pmpBOL float64) (caseResult, error) {
	prop, err := newProp(earth, jgm2, ws.T, ws.R, ws.V, cfg)
	if err != nil {
		return caseResult{}, err
	}
	wc, err := gnco.Propagate(prop, fineDT, period, attNadir)
	if err != nil {
		return caseResult{}, err
	}
	eLab, totLab := faceEnergies(sat, wc, sun, pmpLab)
	eBOL, totBOL := faceEnergies(sat, wc, sun, pmpBOL)
	span := wc.Span()
	return caseResult{span, eLab, eBOL, totLab / span, totBOL / span, wc}, nil
}

// printCase prints a case's per-face energies and mean power; when ref slices are
// non-nil the GMAT reference is shown alongside for comparison.
func printCase(label string, r caseResult, sat Satellite, refLab, refBOL []float64) {
	start := r.traj.Samples[0].T.Time().Format("02 Jan 2006 15:04:05")
	fmt.Printf("\n=== %s ===\nstart %s, period %.2f s, %d samples\n", label, start, r.span, r.traj.Len())
	if refLab != nil {
		fmt.Printf("per-cell energy over the period [J]   (ours vs GMAT Power_Calc)\n")
		fmt.Printf("  face  cells     Lab ours    Lab ref     BOL ours    BOL ref\n")
		for f := range faceNames {
			fmt.Printf("  %-3s   %2d    %10.2f %10.2f   %10.2f %10.2f\n",
				faceNames[f], sat.Cells[f], r.eLab[f], refLab[f], r.eBOL[f], refBOL[f])
		}
		fmt.Printf("mean power [W]: Lab %.2f (ref %.2f)   BOL %.2f (ref %.2f)\n",
			r.meanLab, refMeanLab, r.meanBOL, refMeanBOL)
		return
	}
	fmt.Printf("per-cell energy over the period [J]\n  face  cells     Lab         BOL\n")
	for f := range faceNames {
		fmt.Printf("  %-3s   %2d    %10.2f %10.2f\n", faceNames[f], sat.Cells[f], r.eLab[f], r.eBOL[f])
	}
	fmt.Printf("mean power [W]: Lab %.2f   BOL %.2f\n", r.meanLab, r.meanBOL)
}

// writeSeriesFile renders the per-sample power table to path.
func writeSeriesFile(path string, traj *gnco.Trajectory, sun cosmos.Ephemeris, sat Satellite, pmpLab, pmpBOL float64) error {
	var buf bytes.Buffer
	writePowerSeries(&buf, traj, sun, sat, pmpLab, pmpBOL)
	return os.WriteFile(path, buf.Bytes(), 0777)
}

// newProp builds a fresh force model (JGM2 4×4 with the nutation cache) and
// propagator anchored at epoch with state r, v.
func newProp(earth *cosmos.Body, jgm2 *cosmos.Harmonics, epoch cosmos.Epoch, r, v md3.Vec, cfg gnco.PropagatorConfig) (*gnco.OrbitPropagator, error) {
	fm := gnco.NewForceModel(earth)
	fm.SetHarmonics(jgm2)
	fm.SetNutationInterval(120)
	return gnco.NewOrbitPropagator(fm, epoch, r, v, cfg)
}

// worstCasePeriod returns the index of the coarse sample that starts the
// one-orbital-period window of least generated energy (Lab). Window starts are
// restricted to the first limitSec seconds after epoch0 (0 = the whole span).
func worstCasePeriod(sat Satellite, traj *gnco.Trajectory, sun cosmos.Ephemeris, pmp, period float64, epoch0 cosmos.Epoch, limitSec float64) int {
	n := traj.Len()
	t := make([]float64, n)   // elapsed seconds
	cum := make([]float64, n) // cumulative generated energy [J]
	perCell := make([]float64, len(sat.Normals))
	prevP := 0.0
	for i, s := range traj.Samples {
		sat.facePowerPerCell(s, sun, pmp, perCell)
		p := sat.totalPower(perCell)
		t[i] = s.T.Sub(epoch0)
		if i > 0 {
			cum[i] = cum[i-1] + 0.5*(p+prevP)*(t[i]-t[i-1])
		}
		prevP = p
	}
	best, bestI := math.Inf(1), 0
	for i := 0; i < n; i++ {
		if limitSec > 0 && t[i] > limitSec {
			break
		}
		tEnd := t[i] + period
		if tEnd > t[n-1] {
			break
		}
		if e := interp(t, cum, tEnd) - cum[i]; e < best {
			best, bestI = e, i
		}
	}
	return bestI
}

// interp linearly interpolates y(x) for sorted xs; x is assumed within [xs[0], xs[len-1]].
func interp(xs, ys []float64, x float64) float64 {
	hi := len(xs) - 1
	for lo := 0; lo < hi; {
		mid := (lo + hi) / 2
		if xs[mid] < x {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	if hi == 0 {
		return ys[0]
	}
	x0, x1 := xs[hi-1], xs[hi]
	if x1 == x0 {
		return ys[hi]
	}
	u := (x - x0) / (x1 - x0)
	return ys[hi-1] + u*(ys[hi]-ys[hi-1])
}

// faceEnergies integrates each face's per-cell power [W] over the trajectory by
// trapezoid, returning per-cell energies [J] and the cell-weighted total [J].
func faceEnergies(sat Satellite, traj *gnco.Trajectory, sun cosmos.Ephemeris, pmp float64) (perCellE []float64, total float64) {
	nf := len(sat.Normals)
	perCellE = make([]float64, nf)
	prev := make([]float64, nf)
	cur := make([]float64, nf)
	sat.facePowerPerCell(traj.Samples[0], sun, pmp, prev)
	for i := 1; i < traj.Len(); i++ {
		dt := traj.Samples[i].T.Sub(traj.Samples[i-1].T)
		sat.facePowerPerCell(traj.Samples[i], sun, pmp, cur)
		for f := 0; f < nf; f++ {
			perCellE[f] += 0.5 * (cur[f] + prev[f]) * dt
		}
		prev, cur = cur, prev
	}
	for f := 0; f < nf; f++ {
		total += float64(sat.Cells[f]) * perCellE[f]
	}
	return perCellE, total
}

// writePowerSeries emits a Power_Calc-style table: per-cell face power [W] and
// the cell-weighted total, for the Laboratory and BOL cases.
func writePowerSeries(w io.Writer, traj *gnco.Trajectory, sun cosmos.Ephemeris, sat Satellite, pmpLab, pmpBOL float64) {
	const layout = "2006-01-02 15:04:05.000"
	fmt.Fprint(w, "fecha\ttiempo")
	for _, c := range []string{"Lab", "BOL"} {
		for _, f := range faceNames {
			fmt.Fprintf(w, "\t%s-%s-W", c, f)
		}
		fmt.Fprintf(w, "\t%s-Total-W", c)
	}
	fmt.Fprintln(w)

	t0 := traj.Samples[0].T
	lab := make([]float64, len(faceNames))
	bol := make([]float64, len(faceNames))
	for _, s := range traj.Samples {
		fmt.Fprintf(w, "%s\t%.3f", s.T.Time().Format(layout), s.T.Sub(t0))
		sat.facePowerPerCell(s, sun, pmpLab, lab)
		sat.facePowerPerCell(s, sun, pmpBOL, bol)
		for f := range faceNames {
			fmt.Fprintf(w, "\t%.3f", lab[f])
		}
		fmt.Fprintf(w, "\t%.3f", sat.totalPower(lab))
		for f := range faceNames {
			fmt.Fprintf(w, "\t%.3f", bol[f])
		}
		fmt.Fprintf(w, "\t%.3f\n", sat.totalPower(bol))
	}
}
