// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.20
//

package main

import (
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"time"

	m "github.com/mkhts/gortk"
)

func main() {

	// Parse command line arguments
	args, err := parseArgs()
	if err != nil {
		flag.Usage()
		os.Exit(1)
	}

	// Run the main application
	if err := runApplication(args); err != nil {
		m.PrintE(err)
		os.Exit(1)
	}
}

// Main application processing
func runApplication(args cmdOpt) error {

	// Load input files
	obs, nav, baseObs, err := loadInputFiles(args)
	if err != nil {
		return fmt.Errorf("failed to load input files: %w", err)
	}

	if m.DBG_ >= 1 {
		m.PrintA("--- obs data (%s)---\n", filepath.Base(args.obsFn))
		fmt.Println(obs)
		if baseObs != nil {
			m.PrintA("--- obs data (%s)---\n", filepath.Base(args.baseObsFn))
			fmt.Println(baseObs)
		}
	}
	if m.DBG_ >= 2 {
		m.PrintA("--- nav data (%s)---\n", filepath.Base(args.navFn))
		fmt.Println(nav)
	}

	// Prepare output file
	pos, err := prepareOutput(args)
	if err != nil {
		return fmt.Errorf("failed to prepare output: %w", err)
	}
	defer closeOutput(pos)

	// Print header
	if !args.noPosHeader {
		printPosHeader(pos, os.Args[0], args.mode, args.obsFn, args.navFn, args.baseObsFn, args.ts, args.te, args.basePos, obs)
	}

	// Process epochs
	return processEpochs(args, obs, nav, baseObs, pos)
}

// Load input files
func loadInputFiles(args cmdOpt) (*m.Obs, *m.Nav, *m.Obs, error) {

	obs, err := readObs(args.obsFn)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("failed to read observation file: %w", err)
	}

	nav, err := readNav(args.navFn)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("failed to read navigation file: %w", err)
	}

	var baseObs *m.Obs
	if args.mode == m.DGPS || args.mode == m.RTK {
		baseObs, err = readObs(args.baseObsFn)
		if err != nil {
			return nil, nil, nil, fmt.Errorf("failed to read base observation file: %w", err)
		}
	}

	return obs, nav, baseObs, nil
}

// Prepare output file
func prepareOutput(args cmdOpt) (io.WriteCloser, error) {

	// Use stdout if no output file is specified
	if len(args.posFn) == 0 {
		return &nopCloser{os.Stdout}, nil
	}

	// Create output file
	posf, err := os.Create(args.posFn)
	if err != nil {
		return nil, fmt.Errorf("failed to create output file: %w", err)
	}
	return posf, nil
}

// Close output file
func closeOutput(pos io.WriteCloser) {
	if pos != nil {
		pos.Close()
	}
}

// Process epochs
func processEpochs(args cmdOpt, obs *m.Obs, nav *m.Nav, baseObs *m.Obs, pos io.Writer) error {

	// Variables to maintain previous information
	state := &epochState{
		prevFloatSol: nil,
		prevRovObsE:  nil,
		prevBaseObsE: nil,
		prevAmbSol:   nil,
	}

	// Process each epoch
	for _, obse := range obs.DatE {
		if err := processSingleEpoch(args, obse, nav, baseObs, state, pos); err != nil {
			m.PrintB(obse.Time, "Error processing epoch: %s\n", err.Error())
			continue
		}
	}

	return nil
}

// Process single epoch
func processSingleEpoch(args cmdOpt, obse *m.ObsE, nav *m.Nav, baseObs *m.Obs, state *epochState, pos io.Writer) error {

	// Filter epochs
	if !shouldProcessEpoch(obse, args) {
		return nil
	}

	m.PrintD(2, "\n>>> %s\n", obse.Time.ToTime().UTC())

	// Process base station
	baseSpp, baseObsE, err := processBaseStation(args, obse, nav, baseObs)
	if err != nil {
		return fmt.Errorf("base station processing failed: %w", err)
	}

	// Process rover station
	rovSpp, err := processRoverStation(args, obse, nav, baseSpp)
	if err != nil {
		return fmt.Errorf("rover station processing failed: %w", err)
	}

	// RTK processing
	floatSol, ambSol, err := processRTK(args, obse, rovSpp, baseSpp, state)
	if err != nil {
		return fmt.Errorf("RTK processing failed: %w", err)
	}

	// Output results
	printPos(args.mode, rovSpp.Time, rovSpp, baseSpp, floatSol, ambSol, args.ratioThres, args.nddToNs, pos)

	// Update state
	state.prevRovObsE = obse
	state.prevBaseObsE = baseObsE

	return nil
}

// State management for epoch processing
type epochState struct {
	prevFloatSol *m.FloatSol
	prevRovObsE  *m.ObsE
	prevBaseObsE *m.ObsE
	prevAmbSol   *m.AmbSol
}

// Filter epochs
func shouldProcessEpoch(obse *m.ObsE, args cmdOpt) bool {

	// Skip epochs before processing start time
	if obse.Time.Before(args.ts, true) {
		return false
	}

	// Stop after processing end time
	if obse.Time.After(args.te, true) {
		return false
	}

	// Skip epochs that are not divisible by the specified time interval
	if args.ti > 0 && !obse.Time.Divisible(args.ti) {
		return false
	}

	return true
}

// Process base station
func processBaseStation(args cmdOpt, obse *m.ObsE, nav *m.Nav, baseObs *m.Obs) (*m.SppSol, *m.ObsE, error) {

	// Skip if not in a mode that requires base station processing
	if args.mode != m.DGPS && args.mode != m.RTK {
		return nil, nil, nil
	}

	// Get the temporally nearest base station data for the current epoch
	baseObsE, err := baseObs.GetNearest(obse.Time)
	if err != nil {
		return nil, nil, fmt.Errorf("no base data found")
	}

	// Calculate single point positioning for base station
	m.PrintD(2, "\n\t--- spp for base ---\n")
	sppOpt := setSppOpt(&args)
	sppOpt.EpheSelT = &obse.Time
	sppOpt.IsBase = true

	baseSpp, err := m.CalcSpp(baseObsE, nav, sppOpt)
	if err != nil {
		return nil, nil, fmt.Errorf("spp for base station failed: %w", err)
	}

	return baseSpp, baseObsE, nil
}

// Process rover station
func processRoverStation(args cmdOpt, obse *m.ObsE, nav *m.Nav, baseSpp *m.SppSol) (*m.SppSol, error) {

	m.PrintD(2, "\n\t--- spp for rover ---\n")
	sppOpt := setSppOpt(&args)
	sppOpt.IsBase = false

	// Set DGPS correction values obtained from base station SPP calculation to rover SPP calculation options
	if baseSpp != nil {
		sppOpt.DgpsCorr = baseSpp.DgpsCorr
	}

	rovSpp, err := m.CalcSpp(obse, nav, sppOpt)
	if err != nil {
		return nil, fmt.Errorf("spp for rover station failed: %w", err)
	}

	return rovSpp, nil
}

// RTK processing
func processRTK(args cmdOpt, obse *m.ObsE, rovSpp, baseSpp *m.SppSol, state *epochState) (*m.FloatSol, *m.AmbSol, error) {

	if args.mode != m.RTK {
		return nil, nil, nil
	}

	// Calculate float solution
	floatSol, err := calculateFloatSolution(args, obse, rovSpp, baseSpp, state)
	if err != nil {
		return nil, nil, err
	}

	// Solve ambiguity
	ambSol, err := solveAmbiguity(args, obse, rovSpp, baseSpp, floatSol, state)
	if err != nil {
		return nil, nil, err
	}

	// If not fixed, try removing Glonass
	if args.sys.Contains('R') && args.ratioThres > 0 && ambSol.Ratio < args.ratioThres {
		m.PrintD(2, "\n\t--- RETRY W/O GLONASS ---\n")
		removeGlo(rovSpp)
		removeGlo(baseSpp)
		// Calculate float solution (2nd attempt)
		floatSol, err = calculateFloatSolution(args, obse, rovSpp, baseSpp, state)
		if err != nil {
			return nil, nil, err
		}
		// Solve ambiguity (2nd attempt)
		ambSol, err = solveAmbiguity(args, obse, rovSpp, baseSpp, floatSol, state)
		if err != nil {
			return nil, nil, err
		}
	}

	// If still not fixed, try single frequency
	if args.numFreq >= 2 && args.ratioThres > 0 && ambSol.Ratio < args.ratioThres {
		m.PrintD(2, "\n\t--- RETRY WITH SINGLE FREQ ---\n")
		nf := args.numFreq
		args.numFreq = 1
		defer func() {
			args.numFreq = nf
		}()
		// Calculate float solution (3rd attempt)
		floatSol, err = calculateFloatSolution(args, obse, rovSpp, baseSpp, state)
		if err != nil {
			return nil, nil, err
		}
		// Solve ambiguity (3rd attempt)
		ambSol, err = solveAmbiguity(args, obse, rovSpp, baseSpp, floatSol, state)
		if err != nil {
			return nil, nil, err
		}
	}

	return floatSol, ambSol, nil
}

// Calculate float solution
func calculateFloatSolution(args cmdOpt, obse *m.ObsE, rovSpp, baseSpp *m.SppSol, state *epochState) (*m.FloatSol, error) {

	m.PrintD(2, "\n\t--- calc float sol. ---\n")
	m.PrintD(2, "\trover time            : %s\n", obse.Time.ToTime().UTC().Format("2006-01-02T15:04:05.000000"))
	m.PrintD(2, "\trover time (corrected): %s\n", rovSpp.Time.ToTime().UTC().Format("2006-01-02T15:04:05.000000"))
	m.PrintD(2, "\tbase  time (corrected): %s\n", baseSpp.Time.ToTime().UTC().Format("2006-01-02T15:04:05.000000"))

	floatOpt := setFloatOpt(&args)
	floatOpt.PrevFloatSol = state.prevFloatSol
	floatOpt.PrevRovObsE = state.prevRovObsE
	floatOpt.PrevBaseObsE = state.prevBaseObsE

	floatSol, err := m.CalcFloat(rovSpp, baseSpp, &args.basePos, floatOpt)
	if err != nil {
		return nil, fmt.Errorf("failed to calculate float solution: %w", err)
	}

	// Update state
	state.prevFloatSol = floatSol

	return floatSol, nil
}

// Solve ambiguity
func solveAmbiguity(args cmdOpt, obse *m.ObsE, rovSpp, baseSpp *m.SppSol, floatSol *m.FloatSol, state *epochState) (*m.AmbSol, error) {

	if args.ratioThres <= 0 {
		return nil, nil
	}

	m.PrintD(2, "\n\t--- solve ambiguity ---\n")
	m.PrintD(2, "\trover time            : %s\n", obse.Time.ToTime().UTC().Format("2006-01-02T15:04:05.000000"))
	m.PrintD(2, "\trover time (corrected): %s\n", rovSpp.Time.ToTime().UTC().Format("2006-01-02T15:04:05.000000"))

	ambOpt := setAmbOpt(&args)
	ambOpt.PrevAmbSol = state.prevAmbSol

	ambSol, err := m.SolveAmb(rovSpp, baseSpp, floatSol, &args.basePos, ambOpt)
	if err != nil {
		return nil, fmt.Errorf("failed to solve ambiguity: %w", err)
	}

	// Update state
	state.prevAmbSol = ambSol

	return ambSol, nil
}

// nopCloser - WriteCloser that ignores close operations
type nopCloser struct {
	io.Writer
}

func (nopCloser) Close() error { return nil }

// Structure to hold command line argument information
type cmdOpt struct {
	obsFn           string
	navFn           string
	baseObsFn       string
	posFn           string
	mode            m.Mode
	ts, te          time.Time
	ti              int
	noPosHeader     bool
	sys             m.SysVar
	cnMask          float64
	elMask          float64
	basePos         m.PosXYZ
	exSats          m.SatVar
	wghMode         int
	noChiTest       bool
	maxDop          float64
	maxRes          float64
	numFreq         int
	gByXyz          bool
	wghModeF        int
	instan          bool
	gjMix           bool
	stdCp           float64
	stdPr           float64
	ratioThres      float64
	nddToNs         bool
	noTrop          bool
	retNum          int
	maxElevToRemove float64
}

// Parse command line arguments
func parseArgs() (a cmdOpt, err error) {
	flag.Usage = func() {
		m.PrintA(`
[Usage]
	%s [Options] [-p 0]                                 rover.obs          nav_file.nav (for SPP)
	%s [Options]  -p 1  -l "base_lat base_lon base_hei" rover.obs base.obs nav_file.nav (for DGPS)
	%s [Options] [-p 2] -l "base_lat base_lon base_hei" rover.obs base.obs nav_file.nav (for RTK)

[Options]
`, filepath.Base(os.Args[0]), filepath.Base(os.Args[0]), filepath.Base(os.Args[0]))
		flag.PrintDefaults()
	}
	sOpt := m.NewSppOpt()
	flag.Var(&a.sys, "sys", "Satellite systems to use for calculation. G(GPS), J(QZSS), E(Galileo), R(Glonass), C(Beidou). Comma-separated without spaces. Default: G,J,E,R,C")
	flag.Var(&a.mode, "p", "Calculation mode. 0(SPP), 1(DGPS), 2(RTK)")
	var ts_, te_ m.TimeStr
	flag.TextVar(&ts_, "ts", m.NewTimeStr(time.Time{}), "Start epoch specification. Enclose in quotes like -ts \"2023/01/01 00:00:00.000\"")
	flag.TextVar(&te_, "te", m.NewTimeStr(time.Now().UTC()), "End epoch specification. Enclose in quotes like -te \"2023/01/02 00:00:00.000\". This epoch is also included.")
	flag.IntVar(&a.ti, "ti", 0, "Calculation interval. Calculation is executed when the epoch's second value is divisible by the specified value. Integer only. Omit or set to 0 to calculate all epochs.")
	flag.StringVar(&a.posFn, "o", "", "Output pos file path. If not specified, output to stdout.")
	flag.BoolVar(&a.noPosHeader, "nh", false, "Do not output header section of pos file.")
	flag.Var(&a.exSats, "ex", "List of satellites to exclude. Comma-separated satellite names without spaces like C02,E14.")
	flag.Float64Var(&a.cnMask, "cn", sOpt.CnMask, "Signal strength mask [dB]. Set to 0 for no mask.")
	flag.Float64Var(&a.elMask, "m", sOpt.ElMask, "Elevation mask [deg]. Set to 0 for no mask.")
	flag.IntVar(&a.wghMode, "w", sOpt.WghMode, "Weighting method for SPP calculation. 0(no weighting),1(RTKLIB method),2(RTK core method),3(GPS practical programming book method)")
	flag.BoolVar(&a.gByXyz, "gx", sOpt.GByXyz, "Construct design matrix in XYZ coordinate system. Usually ENU coordinate system. Results are the same. For development to compare with RTKLIB XYZ system.")
	flag.Float64Var(&a.maxDop, "d", sOpt.MaxDop, "Skip calculation and output no results when GDOP exceeds this value. Set to 0 to always calculate regardless of GDOP.")
	flag.BoolVar(&a.noChiTest, "nx2", sOpt.NoChiTest, "Specify to not perform solution evaluation (exclusion) by chi-square test. Default is to perform.")
	flag.Float64Var(&a.maxRes, "mr", sOpt.MaxRes, "Threshold residual for excluding satellite with maximum residual in SPP calculation. Set to 0 to not exclude. Default is no exclusion.")
	var basePosLLH m.PosLLH
	flag.Var(&basePosLLH, "l", "Base station latitude/longitude/ellipsoidal height. Enclose in quotes like -l \"35.73101206 139.7396917 80.33\"")
	fOpt := m.NewFloatOpt()
	flag.IntVar(&a.numFreq, "f", fOpt.NumFreq, "Number of carrier frequencies to use for RTK calculation")
	flag.IntVar(&a.wghModeF, "r", fOpt.WghMode, "Weighting method (observation error covariance matrix) for RTK calculation. 0(no weighting),1(RTKLIB method),2(RTK core method)")
	flag.BoolVar(&a.instan, "i", fOpt.Instan, "Calculate float solution independently for each epoch. Do not use previous epoch results.")
	flag.BoolVar(&a.gjMix, "gj", fOpt.GJMix, "Treat G and J as the same satellite system")
	flag.Float64Var(&a.stdCp, "stdCp", fOpt.StdCp, "Carrier phase noise (standard deviation) specification [m] when weighting method is RTK core method")
	flag.Float64Var(&a.stdPr, "stdPr", fOpt.StdPr, "Pseudorange noise (standard deviation) specification [m] when weighting method is RTK core method")
	aOpt := m.NewAmbOpt()
	flag.Float64Var(&a.ratioThres, "v", aOpt.RatioThres, "Ratio test threshold for FIX determination. Set to 0 to output float solution without AR.")
	var dbg int
	flag.IntVar(&dbg, "x", 0, "Debug information display. Specify level value. 0(OFF), 1(display), 2(detailed display), 3(more detailed), 4(most detailed)")
	flag.BoolVar(&a.nddToNs, "ndd", false, "Output number of double-difference pairs when solving ambiguity instead of number of satellites (for development)")
	flag.IntVar(&a.retNum, "rn", aOpt.MaxRetNum, "Maximum retry count when solving integer ambiguity")
	flag.BoolVar(&a.noTrop, "ntr", false, "Do not perform tropospheric correction")
	flag.Float64Var(&a.maxElevToRemove, "me", aOpt.MaxElevToRemove, "Maximum elevation when retrying by removing low elevation satellites (do not remove satellites with elevation higher than this)")
	flag.Parse()
	switch flag.NArg() {
	case 2:
		a.obsFn = flag.Arg(0)
		a.navFn = flag.Arg(1)
		a.mode = m.SPP
	case 3:
		a.obsFn = flag.Arg(0)
		a.baseObsFn = flag.Arg(1)
		a.navFn = flag.Arg(2)
		if basePosLLH.Lat == 0 {
			return a, fmt.Errorf("the base station position must be specified! (-l option)")
		}
	default:
		return a, fmt.Errorf("too less or many arguments")
	}
	a.ts = time.Time(ts_)
	a.te = time.Time(te_)
	a.basePos = basePosLLH.ToXYZ()
	a.numFreq = min(a.numFreq, m.NFREQ)
	m.DBG_ = dbg
	if m.DBG_ >= 1 && a.mode > 0 {
		m.PrintA("rpos(llh, xyz): %14.9f %14.9f %10.4f, %10.4f %10.4f %10.4f\n", basePosLLH.Lat, basePosLLH.Lon, basePosLLH.Hei, a.basePos.X, a.basePos.Y, a.basePos.Z)
	}
	return
}

// Read observation file
func readObs(fn string) (*m.Obs, error) {
	f, err := os.Open(fn)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	obs, err := m.ReadObs(f)
	if err != nil {
		return nil, err
	}
	return obs, nil
}

// Read navigation file
func readNav(fn string) (*m.Nav, error) {
	f, err := os.Open(fn)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	nav, err := m.ReadNav(f)
	if err != nil {
		return nil, err
	}
	return nav, nil
}

// Print pos file header
func printPosHeader(pos io.Writer, cmd string, mode m.Mode, obsFn, navFn, baseObsFn string, ts, te time.Time, basePos m.PosXYZ, obs *m.Obs) {
	fmt.Fprintf(pos, "%% program   : %s\n", filepath.Base(cmd))
	fmt.Fprintf(pos, "%% inp file  : %s\n", obsFn)
	fmt.Fprintf(pos, "%% inp file  : %s\n", navFn)
	switch mode {
	case m.SPP:
		fmt.Fprintf(pos, "%% obs start : %s\n", getObsStart(obs, ts))
		fmt.Fprintf(pos, "%% obs end   : %s\n", getObsEnd(obs, te))
		fmt.Fprintf(pos, "%%  GPST                 latitude(deg) longitude(deg)  height(m)   Q  ns      clk_bias(s)      isb(E)(s)      isb(R)(s)      isb(C)(s)       gdop       pdop       hdop       vdop\n")
	case m.DGPS:
		fmt.Fprintf(pos, "%% inp file  : %s\n", baseObsFn)
		fmt.Fprintf(pos, "%% obs start : %s\n", getObsStart(obs, ts))
		fmt.Fprintf(pos, "%% obs end   : %s\n", getObsEnd(obs, te))
		llh := basePos.ToLLH()
		fmt.Fprintf(pos, "%% ref pos   : %.8f %.8f %.3f\n", m.ToDeg(llh.Lat), m.ToDeg(llh.Lon), llh.Hei)
		fmt.Fprintf(pos, "%%  GPST                 latitude(deg) longitude(deg)  height(m)   Q  ns      clk_bias(s)      isb(E)(s)      isb(R)(s)      isb(C)(s)       gdop       pdop     age(s)       vdop\n")
	case m.RTK:
		fmt.Fprintf(pos, "%% inp file  : %s\n", baseObsFn)
		fmt.Fprintf(pos, "%% obs start : %s\n", getObsStart(obs, ts))
		fmt.Fprintf(pos, "%% obs end   : %s\n", getObsEnd(obs, te))
		llh := basePos.ToLLH()
		fmt.Fprintf(pos, "%% ref pos   : %.8f %.8f %.3f\n", m.ToDeg(llh.Lat), m.ToDeg(llh.Lon), llh.Hei)
		fmt.Fprintf(pos, "%%  GPST                 latitude(deg) longitude(deg)  height(m)   Q  ns      clk_bias(s)      isb(E)(s)      isb(R)(s)      isb(C)(s)       gdop       pdop     age(s)      ratio\n")
	}
}

// Return processing start epoch date and time as string
func getObsStart(obs *m.Obs, ts time.Time) string {
	t := obs.DatE[0].Time
	for _, obse := range obs.DatE {
		if obse.Time.Before(ts, true) {
			continue
		}
		t = obse.Time
		break
	}
	return fmt.Sprintf("%s(UTC) (week%d %7.1fs)(GPST)", t.ToTime().UTC().Format("2006/01/02 15:04:05.000"), t.Week, t.Sec)
}

// Return processing end epoch date and time as string
func getObsEnd(obs *m.Obs, te time.Time) string {
	t := obs.DatE[len(obs.DatE)-1].Time
	for _, obse := range obs.DatE {
		if obse.Time.After(te, true) {
			break
		}
		t = obse.Time
	}
	return fmt.Sprintf("%s(UTC) (week%d %7.1fs)(GPST)", t.ToTime().UTC().Format("2006/01/02 15:04:05.000"), t.Week, t.Sec)
}

// Output POS file
func printPos(mode m.Mode, rcvt m.GTime, uspp, rspp *m.SppSol, fsol *m.FloatSol, asol *m.AmbSol, rath float64, nddToNs bool, pos io.Writer) {
	gdop := uspp.Dop["gdop"]
	pdop := uspp.Dop["pdop"]
	hdop := uspp.Dop["hdop"]
	vdop := uspp.Dop["vdop"]
	ns := len(uspp.Sats)
	llh := uspp.Pos.ToLLH()
	rcvtStr := rcvt.ToTime().UTC().Format("2006/01/02 15:04:05.000000")
	roundSec := true // Whether to round time to milliseconds
	if roundSec {
		rcvt2 := m.GTime{
			Week: rcvt.Week,
			Sec:  math.Round(rcvt.Sec*1000) / 1000,
		}
		rcvtStr = rcvt2.ToTime().UTC().Format("2006/01/02 15:04:05.000")
	}
	switch mode {
	case m.SPP:
		Q := 5
		fmt.Fprintf(pos, "%s %13.9f %14.9f %10.4f %3d %3d %16.4f %14.4f %14.4f %14.4f %10.3f %10.3f %10.3f %10.3f\n", rcvtStr, m.ToDeg(llh.Lat), m.ToDeg(llh.Lon), llh.Hei, Q, ns, uspp.Clk[0], uspp.Clk[1], uspp.Clk[2], uspp.Clk[3], gdop, pdop, hdop, vdop)
	case m.DGPS:
		Q := 4
		age := uspp.Time.ToTime().Sub(rspp.Time.ToTime()).Seconds()
		fmt.Fprintf(pos, "%s %13.9f %14.9f %10.4f %3d %3d %16.4f %14.4f %14.4f %14.4f %10.3f %10.3f %10.7f %10.3f\n", rcvtStr, m.ToDeg(llh.Lat), m.ToDeg(llh.Lon), llh.Hei, Q, ns, uspp.Clk[0], uspp.Clk[1], uspp.Clk[2], uspp.Clk[3], gdop, pdop, age, vdop)
	case m.RTK:
		Q := 2
		ratio := 0.0
		ns = fsol.NumSats
		if nddToNs && asol != nil {
			ns = len(asol.DDPairs) // Output number of double-difference pairs used when solving integer ambiguity instead of number of satellites
		}
		if rath > 0 {
			llh = asol.Pos.ToLLH()
			ratio = asol.Ratio
			if ratio >= rath {
				Q = 1
			}
		} else {
			llh = fsol.Pos.ToLLH()
		}
		fmt.Fprintf(pos, "%s %13.9f %14.9f %10.4f %3d %3d %16.4f %14.4f %14.4f %14.4f %10.3f %10.3f %10.7f %10.3f\n", rcvtStr, m.ToDeg(llh.Lat), m.ToDeg(llh.Lon), llh.Hei, Q, ns, uspp.Clk[0], uspp.Clk[1], uspp.Clk[2], uspp.Clk[3], gdop, pdop, fsol.Age, ratio)
	}

}

func setSppOpt(args *cmdOpt) *m.SppOpt {
	opt := m.NewSppOpt()
	opt.Sys = args.sys
	opt.ExSats = args.exSats
	opt.CnMask = args.cnMask
	opt.ElMask = args.elMask
	opt.WghMode = args.wghMode
	opt.NoChiTest = args.noChiTest
	opt.MaxDop = args.maxDop
	opt.MaxRes = args.maxRes
	opt.BasePos = &args.basePos
	opt.DgpsCorr = nil
	opt.GByXyz = args.gByXyz
	return opt
}

func setFloatOpt(args *cmdOpt) *m.FloatOpt {
	opt := m.NewFloatOpt()
	opt.NumFreq = args.numFreq
	opt.WghMode = args.wghModeF
	opt.CnMask = args.cnMask
	opt.ElMask = args.elMask
	opt.Instan = args.instan
	opt.GJMix = args.gjMix
	opt.StdCp = args.stdCp
	opt.StdPr = args.stdPr
	opt.SkipTrop = args.noTrop
	return opt
}

func setAmbOpt(args *cmdOpt) *m.AmbOpt {
	opt := m.NewAmbOpt()
	opt.RatioThres = args.ratioThres
	opt.SkipTrop = args.noTrop
	opt.MaxRetNum = args.retNum
	opt.MaxElevToRemove = args.maxElevToRemove
	return opt
}

func removeGlo(sol *m.SppSol) {
	sats := []m.SatType{}
	for _, s := range sol.Sats {
		if s.Sys() != 'R' {
			sats = append(sats, s)
		}
	}
	sol.Sats = sats
}
