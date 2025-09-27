// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.14
//

// Implements single point positioning (SPP) calculation for GNSS positioning.

package gortk

import (
	"fmt"
	"math"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/mat"
)

// CalcSpp performs single point positioning calculation using pseudorange observations
// It computes receiver position, clock bias, and system biases through iterative least squares
//
// Parameters:
//   - obse: Single epoch observation data
//   - nav: Navigation data (all epochs)
//   - opt: SPP calculation options and parameters
//
// Returns:
//   - SppSol: SPP solution containing position, clock bias, and quality metrics
//   - error: Any error encountered during processing
func CalcSpp(
	obse *ObsE, // Single epoch observation data
	nav *Nav, // Navigation data (all epochs)
	opt *SppOpt, // Calculation options
) (*SppSol, error) {

	// Initialize result structure
	rslt := NewSppSol()

	// Select valid satellites for calculation
	err := selectValidSatellites(obse, nav, opt, rslt)
	if err != nil {
		return nil, fmt.Errorf("selectValidSatellites() failed, err=%v", err)
	}

	// Solve observation equations iteratively
	err = solveSppEquations(obse, opt, rslt)
	if err != nil {
		return nil, fmt.Errorf("solveSppEquations() failed, err=%v", err)
	}

	// Validate solution quality
	err = validateSppSol(rslt, opt)
	if err != nil {
		return nil, fmt.Errorf("validateSppSol() failed, err=%v", err)
	}

	return rslt, nil
}

// SPP uses only single frequency (L1) data
const FREQ = 0

// Calculation constants for SPP processing
const (
	MAX_LOOP_COUNT           = 10    // Maximum number of iteration loops
	EARLY_LOOP_SKIP          = 3     // Number of initial loops to skip elevation calculation
	CONVERGENCE_THRESHOLD    = 0.001 // Convergence threshold [m]
	MIN_WEIGHT               = 0.001 // Minimum weight value
	MIN_ELEVATION_FOR_WEIGHT = 5.0   // Minimum elevation angle for weight calculation [deg]
)

// Satellite system indices for clock bias parameter counting
const (
	GPS_GLONASS_INDEX = 0 // GPS/QZSS systems
	GALILEO_INDEX     = 1 // Galileo system
	GLONASS_INDEX     = 2 // GLONASS system
	BEIDOU_INDEX      = 3 // BeiDou system
)

// SppOpt contains options and parameters for single point positioning calculation
// These parameters control satellite selection, weighting schemes, and quality control
type SppOpt struct {
	Sys       []SysType           // List of satellite systems to use in calculation
	ExSats    []SatType           // List of satellites to exclude from calculation
	CnMask    float64             // Signal strength mask [dB]
	ElMask    float64             // Elevation mask [deg]
	WghMode   int                 // Weighting scheme: 0(OFF), 1(RTKLIB), 2(RTK Core), 3(GPS Programming)
	NoChiTest bool                // If true, skip chi-square test for solution validation (rover only)
	MaxDop    float64             // Maximum allowed GDOP value. Skip calculation if exceeded (rover only)
	MaxRes    float64             // Maximum residual threshold for satellite exclusion. 0 means no exclusion
	BasePos   *PosXYZ             // Base station position for DGPS correction calculation (optional)
	DgpsCorr  map[SatType]float64 // DGPS correction values for each satellite (optional)
	GByXyz    bool                // If true, build design matrix in XYZ coordinates, otherwise ENU (debug use)
	EpheSelT  *GTime              // Time for ephemeris selection. If nil, use observation time
	IsBase    bool                // Whether this is base station SPP calculation
}

// NewSppOpt creates a new SppOpt with default values
// Default values are tuned for typical GNSS positioning applications
func NewSppOpt() *SppOpt {
	return &SppOpt{
		Sys:       []SysType{}, // Use all available systems
		ExSats:    []SatType{}, // No excluded satellites
		CnMask:    35,          // Signal strength threshold [dB]
		ElMask:    15,          // Elevation mask [deg]
		WghMode:   1,           // RTKLIB weighting scheme
		NoChiTest: false,       // Perform chi-square test
		MaxDop:    35,          // GDOP threshold
		MaxRes:    0,           // No residual-based exclusion
		BasePos:   nil,         // No base position
		DgpsCorr:  nil,         // No DGPS corrections
		GByXyz:    false,       // Use ENU coordinates
		EpheSelT:  nil,         // Use observation time
		IsBase:    false,       // Rover station
	}
}

// SppSol contains the results of single point positioning calculation
// It holds position estimates, clock biases, quality metrics, and observation data
type SppSol struct {
	Pos      PosXYZ                      // Receiver position (positioning result)
	Time     GTime                       // Time after receiver clock error correction
	Clk      []float64                   // Receiver clock bias and inter-system biases (ISB): 0:receiver, 1:E, 2:R, 3:C
	Sats     []SatType                   // List of satellites used in calculation
	Cov      [4][4]float64               // Estimation error covariance matrix ((G^T W G)^-1)
	Dop      map[string]float64          // Dilution of precision values: 'gdop', 'pdop', 'hdop', 'vdop'
	DgpsCorr map[SatType]float64         // DGPS correction values for pseudorange (base station only)
	SatPos   map[SatType]PosXYZ          // 3D satellite positions
	Elev     map[SatType]float64         // Satellite elevation angles
	Res      map[SatType]float64         // Final residuals at convergence
	Ephe     map[SatType]*Ephe           // Satellite ephemeris data
	Pr       map[SatType][NFREQ]float64  // Pseudorange observations
	Cp       map[SatType][NFREQ]float64  // Carrier phase observations
	Sn       map[SatType][NFREQ]float64  // Signal strength values
	Freq     map[SatType][NFREQ]float64  // Carrier frequencies
	LLI      map[SatType][NFREQ]byte     // Loss of lock indicators
	Code     map[SatType][NFREQ]CodeType // Observation codes
	SatClk   map[SatType]float64         // Satellite clock errors
	DesMat   mat.Matrix                  // Design matrix
	ResVec   mat.Vector                  // Residual vector
	WghMat   mat.Matrix                  // Weight matrix
}

// NewSppSol creates a new empty SppSol structure with initialized maps
func NewSppSol() *SppSol {
	return &SppSol{
		Pos:  PosXYZ{},
		Time: GTime{},
		Clk:  []float64{},
		Sats: []SatType{},
		Cov:  [4][4]float64{},
		Dop: map[string]float64{
			"gdop": 0,
			"pdop": 0,
			"hdop": 0,
			"vdop": 0,
		},
		DgpsCorr: map[SatType]float64{},
		SatPos:   map[SatType]PosXYZ{},
		Elev:     map[SatType]float64{},
		Res:      map[SatType]float64{},
		Ephe:     map[SatType]*Ephe{},
		Pr:       map[SatType][NFREQ]float64{},
		Cp:       map[SatType][NFREQ]float64{},
		Sn:       map[SatType][NFREQ]float64{},
		Freq:     map[SatType][NFREQ]float64{},
		LLI:      map[SatType][NFREQ]byte{},
		Code:     map[SatType][NFREQ]CodeType{},
		SatClk:   map[SatType]float64{},
		DesMat:   nil,
		ResVec:   nil,
		WghMat:   nil,
	}
}

// selectValidSatellites selects valid satellites for SPP calculation
// It applies various filtering criteria including system selection, signal strength, and health status
func selectValidSatellites(obse *ObsE, nav *Nav, opt *SppOpt, rslt *SppSol) (err error) {

	// Loop through satellites in observation data
	for _, sat := range Sorted(obse.Sats()) {

		// Get observation data for this satellite
		obss := obse.DatS[sat]

		// Skip if satellite system not in specified list
		if opt.Sys != nil && !slices.Contains(opt.Sys, sat.Sys()) {
			continue
		}

		// Skip if satellite in exclusion list
		if opt.ExSats != nil && slices.Contains(opt.ExSats, sat) {
			PrintD(3, "\t%s: Exclude satellite\n", sat)
			continue
		}

		// Get ephemeris data
		var eph *Ephe
		if opt.EpheSelT == nil {
			eph, err = nav.GetEphe(sat, obse.Time) // Use observation time
		} else {
			eph, err = nav.GetEphe(sat, *opt.EpheSelT) // Use specified time
		}

		// Skip if no ephemeris available
		if err != nil {
			PrintD(3, "\t%s: No ephemeris\n", sat)
			continue
		}

		// Skip if signal strength below threshold
		if opt.CnMask > 0 && obss.Sn[FREQ] < opt.CnMask {
			PrintD(3, "\t%s: C/N Mask (c/n=%f < %f)\n", sat, obss.Sn[FREQ], opt.CnMask)
			continue
		}

		// Skip if satellite not healthy
		svh := eph.Svh
		if sat.Sys() == 'J' { // Ignore health flag for QZSS
			svh &= 0xfffffffe
			PrintD(3, "\t%s: Not healthy, but ignored.\n", sat)
		}
		if svh != 0 {
			PrintD(3, "\t%s: Not healthy\n", sat)
			continue
		}

		// Skip if no valid observation data
		if obss.Pr[FREQ] == 0 || obss.Cp[FREQ] == 0 { // SPP doesn't use carrier phase, but exclude for RTK compatibility
			if DBG_ >= 3 {
				PrintA("\t%s: No pseudorange or carrier phase data.\n", sat)
			}
			continue
		}

		// Skip if DGPS correction not available when required
		if opt.DgpsCorr != nil {
			if _, ok := opt.DgpsCorr[sat]; !ok {
				PrintD(3, "\t%s: No DGPS corr data.\n", sat)
				continue
			}
		}

		if DBG_ >= 3 {
			toe := eph.Toe.ToTime().UTC().Format("2006-01-02 15:04:05.000")
			toc := eph.Toc.ToTime().UTC().Format("2006-01-02 15:04:05.000")
			tot := eph.Tot.ToTime().UTC().Format("2006-01-02 15:04:05.000")
			PrintA("\t%s: Toe: %s, Toc: %s, Tot: %s, Iode: %d\n", sat, toe, toc, tot, eph.Iode)
		}

		// Select satellite that passed all criteria
		rslt.Pr[sat] = obss.Pr
		rslt.Cp[sat] = obss.Cp
		rslt.Sn[sat] = obss.Sn
		rslt.LLI[sat] = obss.LLI
		rslt.Code[sat] = obss.Code
		rslt.Ephe[sat] = eph
		rslt.Sats = append(rslt.Sats, sat)
		rslt.Freq[sat] = obss.Freq
		// GLONASS frequency correction
		if sat.Sys() == 'R' {
			rslt.Freq[sat] = getGloFreq(eph, obss.Freq)
		}
	}

	PrintD(2, "\tsat: %d / %d\n", len(rslt.Sats), len(obse.DatS))

	// Error if no satellites selected
	if len(rslt.Sats) == 0 {
		return fmt.Errorf("no valid satellites")
	}

	return nil
}

// solveSppEquations sets up observation equations and solves them iteratively
// It performs least squares adjustment with convergence checking and quality control
func solveSppEquations(obse *ObsE, opt *SppOpt, rslt *SppSol) (err error) {

	// Receiver position (initial value: 0, 0, 0)
	var upos PosXYZ

	if DBG_ >= 2 {
		llh := upos.ToLLH()
		llh.Lat = ToDeg(llh.Lat)
		llh.Lon = ToDeg(llh.Lon)
		PrintA("\tupos(init): LLH= %.8f %.8f %.3f, XYZ= %.3f %.3f %.3f\n", llh.Lat, llh.Lon, llh.Hei, upos.X, upos.Y, upos.Z)
	}

	// Receiver clock bias and inter-system biases (ISB)
	clkb := []float64{0, 0, 0, 0, 0} // receiver clock bias, E ISB, R ISB, C ISB, S ISB
	clkbs := 0.0                     // sum of above

	// List of satellites to exclude from calculation
	// Satellites with low elevation or large residuals are added during iteration
	xs := []SatType{}

	// Covariance matrix
	var cov mat.Matrix

	// Design matrix
	var G2 mat.Matrix

	// Residual vector
	var dr2 mat.Vector

	// Weight matrix
	var W mat.Matrix

	// Flag to exit iteration loop
	exitLoop := false

	// Solve observation equations iteratively
	for loop := 0; loop < MAX_LOOP_COUNT; loop++ {

		PrintAIf(DBG_ >= 3 && !exitLoop, "\t--- LOOP: %d ---\n", loop+1)

		// ---------------------------------
		// Setup equations
		// ---------------------------------

		// Count satellites by system (needed to determine number of unknowns)
		nSys := countSatellitesBySystem(rslt.Sats, xs)

		// Number of clock bias parameters (receiver + ISB)
		// Total unknowns = 3 (position) + nClk
		nClk := countClockBiasParameters(nSys)
		PrintAIf(DBG_ >= 3 && !exitLoop, "\tnSys: %v, nClk: %d\n", nSys, nClk)

		// Number of satellites used in calculation (= number of equations)
		n := len(rslt.Sats) - len(xs)

		// Number of unknowns
		nx := 3 + nClk

		// Error if not enough satellites
		if n < nx {
			return fmt.Errorf("not enough satellites :%d < %d", n, nx)
		}

		// Design matrix (recreated each loop as equation count may change)
		G := mat.NewDense(n, nx, nil) // n x nx

		// Residual vector (recreated each loop)
		dr := mat.NewVecDense(n, nil) // n x 1

		// Receiver time corrected for current clock bias estimate
		rcvt2 := GTime{Week: obse.Time.Week, Sec: obse.Time.Sec - clkb[0]/C}
		PrintD(3, "\trcvt2: %v\n", rcvt2.ToTime().UTC().Format("2006/01/02 15:04:05.000000"))

		// Maximum residual in observation equations (calculated later)
		var dmax float64
		var dsat SatType

		// Weight for each satellite
		w := make([]float64, 0, n)

		// Process each satellite (set design matrix and residual vector)
		i := 0 // Counter (len(rslt.Sats) ≠ number of satellites used, so don't use range i)
		for _, sat := range rslt.Sats {

			// Skip if satellite in exclusion list
			if slices.Contains(xs, sat) {
				continue
			}

			// Get ephemeris for this satellite
			eph := rslt.Ephe[sat]

			// Pseudorange observation
			psr := rslt.Pr[sat][FREQ]

			// Apply DGPS correction if available
			if opt.DgpsCorr != nil {
				if prc, ok := opt.DgpsCorr[sat]; ok {
					psr += prc
				}
			}

			// Correct pseudorange for receiver clock bias (keep original for later use)
			psr2 := psr - clkb[0]

			// Correct pseudorange for satellite clock error
			psr2 += satClk(eph, rcvt2, psr2) * C

			// Calculate satellite position
			spos := SatPos(eph, rcvt2, psr2)

			// Recalculate satellite clock error with updated pseudorange
			sclk := satClk(eph, rcvt2, psr2)
			rslt.SatClk[sat] = sclk

			// Calculate DGPS correction if base position available
			if opt.BasePos != nil {
				rib := EucDist(&spos, opt.BasePos) // Distance from base to satellite
				rslt.DgpsCorr[sat] = -(psr2 - rib) // Correction = -(pseudorange - distance), negative for rover
			}

			// Euclidean distance from receiver to satellite
			ri := EucDist(&spos, &upos)

			// Convert satellite position to ENU coordinates
			sposEnu := spos.ToENU(upos)

			// Elevation mask
			if loop >= EARLY_LOOP_SKIP && opt.ElMask > 0 { // Elevation calculation needs receiver position, so skip early loops
				elv := ToDeg(sposEnu.Elevation())
				if elv < opt.ElMask {
					if !slices.Contains(xs, sat) {
						xs = append(xs, sat)
						PrintD(3, "\t%s: elev=%f < %f\n", sat, elv, opt.ElMask)
					}
					continue
				}
			}

			// Set values in the design matrix
			if opt.GByXyz { // If constructing the design matrix in the XYZ system
				G.Set(i, 0, DistDx(&spos, &upos))
				G.Set(i, 1, DistDy(&spos, &upos))
				G.Set(i, 2, DistDz(&spos, &upos))
			} else { // If constructing the design matrix in the ENU system (the positioning result does not change regardless of the system used. Switch as needed for debugging, etc.)
				G.Set(i, 0, -sposEnu.E/ri)
				G.Set(i, 1, -sposEnu.N/ri)
				G.Set(i, 2, -sposEnu.U/ri)
			}
			G.Set(i, 3, 1)

			// Below, set 1 in the design matrix. Branching is necessary depending on the combination of satellite systems.
			for j := 1; j < nClk; j++ {
				G.Set(i, 3+j, 0)
			}
			clkbs = clkb[0] // The value of the receiver clock error to be set in the residual vector also requires branching if there are multiple satellite systems
			if nClk > 1 {   // Below, branching process assuming the design matrix is constructed in the order of G/J (these two satellite systems are not distinguished), E, R, C
				switch sat.Sys() {
				case 'E':
					if nSys[GPS_GLONASS_INDEX] > 0 { // G,E
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					}
				case 'R':
					if nSys[GPS_GLONASS_INDEX] == 0 && nSys[GALILEO_INDEX] > 0 { // E,R
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					} else if nSys[GPS_GLONASS_INDEX] > 0 && nSys[GALILEO_INDEX] == 0 { // G,R
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					} else if nSys[GPS_GLONASS_INDEX] > 0 && nSys[GALILEO_INDEX] > 0 { // G,E,R
						G.Set(i, 5, 1)
						clkbs += clkb[2]
					}
				case 'C':
					if nSys[GPS_GLONASS_INDEX] > 0 && nSys[GALILEO_INDEX] == 0 && nSys[GLONASS_INDEX] == 0 && nSys[BEIDOU_INDEX] > 0 { // G,C
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					} else if nSys[GPS_GLONASS_INDEX] == 0 && nSys[GALILEO_INDEX] > 0 && nSys[GLONASS_INDEX] == 0 && nSys[BEIDOU_INDEX] > 0 { // E,C
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					} else if nSys[GPS_GLONASS_INDEX] == 0 && nSys[GALILEO_INDEX] == 0 && nSys[GLONASS_INDEX] > 0 && nSys[BEIDOU_INDEX] > 0 { // R,C
						G.Set(i, 4, 1)
						clkbs += clkb[1]
					} else if nSys[GPS_GLONASS_INDEX] > 0 && nSys[GALILEO_INDEX] > 0 && nSys[GLONASS_INDEX] == 0 && nSys[BEIDOU_INDEX] > 0 { // G,E,C
						G.Set(i, 5, 1)
						clkbs += clkb[2]
					} else if nSys[GPS_GLONASS_INDEX] > 0 && nSys[GALILEO_INDEX] == 0 && nSys[GLONASS_INDEX] > 0 && nSys[BEIDOU_INDEX] > 0 { // G,R,C
						G.Set(i, 5, 1)
						clkbs += clkb[2]
					} else if nSys[GPS_GLONASS_INDEX] == 0 && nSys[GALILEO_INDEX] > 0 && nSys[GLONASS_INDEX] > 0 && nSys[BEIDOU_INDEX] > 0 { // E,R,C
						G.Set(i, 5, 1)
						clkbs += clkb[2]
					} else { // G,E,R,C
						G.Set(i, 6, 1)
						clkbs += clkb[3]
					}
				}
			}

			// Set residual vector value
			dr.SetVec(i, psr+sclk*C-(ri+clkbs)) // Use original pseudorange observation (psr)

			// Calculate weight for this satellite
			elv := sposEnu.Elevation()
			wg := getWeight(opt.WghMode, sat, elv, eph)
			w = append(w, wg)

			if DBG_ >= 3 && !exitLoop {
				tk_ := rcvt2.ToTime().Sub(eph.Toc.ToTime()).Seconds() - psr/C
				dt_ := eph.Af0 + eph.Af1*tk_ + eph.Af2*tk_*tk_
				rcvt_ := (GTime{Week: obse.Time.Week, Sec: obse.Time.Sec - clkb[0]/C - psr/C - dt_})
				ts_ := rcvt_.ToTime().Format("2006/01/02 15:04:05.000000")
				PrintA("\t%s: trsmt=%s, elev=%8.3f, weight=%8.3f, azim=%8.3f, x=%16.3f, y=%16.3f, z=%16.3f, psr=%14.3f, psr2=%14.3f, sclk*C=%12.3f, dr=%12.3f\n", sat, ts_, ToDeg(elv), wg, ToDeg(sposEnu.Azimuth()), spos.X, spos.Y, spos.Z, psr, psr2, sclk*C, dr.AtVec(i))
			}

			// Store calculation results
			rslt.Elev[sat] = elv
			rslt.Res[sat] = dr.AtVec(i)
			rslt.SatPos[sat] = spos

			// Increment index
			i += 1

		} // for sat

		// Exit loop if already converged
		if exitLoop {
			break
		}

		// Error if not enough equations
		if i < nx {
			return fmt.Errorf("not enough number of equations: %d < %d", i, nx)
		}

		// ---------------------------------
		// Solve equations (least squares)
		// ---------------------------------

		var dx mat.Vector               // Unknown vector
		G2 = G.Slice(0, i, 0, nx)       // Remove excluded satellite rows
		dr2 = dr.SliceVec(0, i)         // Remove excluded satellite rows
		W = mat.NewDiagDense(len(w), w) // Weight matrix
		if DBG_ >= 4 {
			PrintA("G=\n")
			PrintMat(G2)
			PrintA("dr=\n")
			PrintMat(dr2)
			wd := mat.NewVecDense(len(w), w)
			PrintA("W(diag)=\n")
			PrintMat(wd)
		}
		dx, cov, err = SolveLS(G2, dr2, W)
		if err != nil {
			PrintD(2, "\tSolve() failed., err= %s\n", err.Error())
			return err
		}

		if DBG_ >= 4 {
			PrintA("dx=\n")
			PrintMat(dx)
		}

		// Check maximum residual
		dmax, dsat = getMaxResidualAndSatellite(dr2, rslt.Sats, xs)

		// Update receiver position
		updateUserPosition(&upos, dx, opt)

		// Update clock bias and inter-system biases
		updateClockBias(clkb, dx, nClk)

		if DBG_ >= 2 {
			llh := upos.ToLLH()
			llh.Lat = ToDeg(llh.Lat)
			llh.Lon = ToDeg(llh.Lon)
			PrintAIf(nClk == 1, "\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f, s[0]=%.3f\n", loop+1, llh.Lat, llh.Lon, llh.Hei, upos.X, upos.Y, upos.Z, clkb[0])
			PrintAIf(nClk == 2, "\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f, s[0]=%.3f, s[1]=%.3f\n", loop+1, llh.Lat, llh.Lon, llh.Hei, upos.X, upos.Y, upos.Z, clkb[0], clkb[1])
			PrintAIf(nClk == 3, "\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f, s[0]=%.3f, s[1]=%.3f, s[2]=%.3f\n", loop+1, llh.Lat, llh.Lon, llh.Hei, upos.X, upos.Y, upos.Z, clkb[0], clkb[1], clkb[2])
			PrintAIf(nClk == 4, "\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f, s[0]=%.3f, s[1]=%.3f, s[2]=%.3f, s[3]=%.3f\n", loop+1, llh.Lat, llh.Lon, llh.Hei, upos.X, upos.Y, upos.Z, clkb[0], clkb[1], clkb[2], clkb[3])
		}

		// Check convergence (position update < 1mm)
		if isSppConverged(dx, CONVERGENCE_THRESHOLD) {
			PrintD(2, "\tdiff: %f, %f, %f < %f\n", math.Abs(dx.AtVec(0)), math.Abs(dx.AtVec(1)), math.Abs(dx.AtVec(2)), CONVERGENCE_THRESHOLD)
			exitLoop = true // Don't break, run one more loop to update rslt values
		}

		// Exclude satellite with maximum residual if threshold exceeded
		// This has side effects, so disabled by default
		if loop >= EARLY_LOOP_SKIP && opt.MaxRes > 0 && dmax > opt.MaxRes {
			xs = append(xs, dsat)
			PrintD(3, "\t%s: residual=%f > %f\n", dsat, dmax, opt.MaxRes)
		}

		// Error if maximum loop count reached
		if loop+1 == MAX_LOOP_COUNT {
			return fmt.Errorf("number of loop reached max")
		}

	} // for LOOP

	// Error if no calculations performed
	if clkb[0] == 0 {
		return fmt.Errorf("no valid calculations")
	}

	// Remove excluded satellites from data
	for _, xsat := range xs {
		delete(rslt.Pr, xsat)
		delete(rslt.Cp, xsat)
		delete(rslt.Sn, xsat)
		delete(rslt.Freq, xsat)
		delete(rslt.DgpsCorr, xsat)
		delete(rslt.SatPos, xsat)
		delete(rslt.Elev, xsat)
		delete(rslt.Res, xsat)
		delete(rslt.DgpsCorr, xsat)
		delete(rslt.LLI, xsat)
		delete(rslt.Code, xsat)
		delete(rslt.SatClk, xsat)
	}
	sats := make([]SatType, 0, len(rslt.Sats))
	for _, sat := range rslt.Sats {
		if slices.Contains(xs, sat) {
			continue
		} else {
			sats = append(sats, sat)
		}
	}
	rslt.Sats = sats

	// Calculate DOP values
	var GtG mat.Dense
	GtG.Mul(G2.T(), G2)
	var cov2 mat.Dense
	err = cov2.Inverse((&GtG))
	if err != nil {
		return fmt.Errorf("\tfailed to calculate inverse of matrix, G^T G")
	}
	rslt.Dop["gdop"] = math.Sqrt(cov2.At(0, 0) + cov2.At(1, 1) + cov2.At(2, 2) + cov2.At(3, 3))
	rslt.Dop["pdop"] = math.Sqrt(cov2.At(0, 0) + cov2.At(1, 1) + cov2.At(2, 2))
	rslt.Dop["hdop"] = math.Sqrt(cov2.At(0, 0) + cov2.At(1, 1))
	rslt.Dop["vdop"] = math.Sqrt(cov2.At(2, 2))

	// Set values in result structure
	rslt.Pos = upos
	rslt.Time = GTime{Week: obse.Time.Week, Sec: obse.Time.Sec - clkb[0]/C} // Time corrected for receiver clock bias
	rslt.Clk = make([]float64, len(clkb))
	copy(rslt.Clk, clkb)
	for j := 0; j < 4; j++ {
		for k := 0; k < 4; k++ {
			rslt.Cov[j][k] = cov.At(j, k)
		}
	}
	rslt.DesMat = G2
	rslt.ResVec = dr2
	rslt.WghMat = W

	return nil
}

// validateSppSol validates the SPP solution quality using chi-square test and GDOP check
func validateSppSol(rslt *SppSol, opt *SppOpt) (err error) {

	// Chi-square test
	if !opt.NoChiTest && !opt.IsBase {
		nM, nX := rslt.DesMat.Dims() // nM: num of measurements, nX: num of parameters
		if nM > nX {
			v := mat.NewVecDense(rslt.ResVec.Len(), nil)
			for k := range v.Len() {
				v.SetVec(k, rslt.ResVec.AtVec(k)*math.Sqrt(rslt.WghMat.At(k, k)))
			}
			vv := mat.Dot(v, v)
			v.MulVec(rslt.WghMat, rslt.ResVec)
			if DBG_ >= 4 {
				PrintA("v=\n")
				PrintMat(v)
			}
			if vv > ChiSqr(nM-nX-1) {
				return fmt.Errorf("chi-Square test failed: nM=%d, nX=%d, |w dr|=%f > %f", nM, nX, vv, ChiSqr(nM-nX-1))
			}
			PrintD(3, "\tChi-Square Test: |W dr|^2=%.3f <= %f\n", vv, ChiSqr(nM-nX-1))
		} else {
			PrintD(3, "\tChi-Square Test: undone. nM=%d <= nX=%d\n", nM, nX)
		}
	}

	// Check if GDOP exceeds threshold
	if opt.MaxDop > 0 && !opt.IsBase {
		if rslt.Dop["gdop"] > opt.MaxDop {
			return fmt.Errorf("GDOP exceeded threshold, GDOP=%.3f > %f", rslt.Dop["gdop"], opt.MaxDop)
		} else {
			PrintD(3, "\tGDOP Test: %.3f <= %f\n", rslt.Dop["gdop"], opt.MaxDop)
		}
	}

	return nil
}

// countSatellitesBySystem counts satellites by system type
func countSatellitesBySystem(sats []SatType, excluded []SatType) [4]int {
	nSys := [4]int{0, 0, 0, 0} // G/J, E, R, C
	for _, sat := range sats {
		if slices.Contains(excluded, sat) {
			continue
		}
		sys := sat.Sys()
		switch {
		case sys == 'G' || sys == 'J':
			nSys[GPS_GLONASS_INDEX] += 1
		case sys == 'E':
			nSys[GALILEO_INDEX] += 1
		case sys == 'R':
			nSys[GLONASS_INDEX] += 1
		case sys == 'C':
			nSys[BEIDOU_INDEX] += 1
		}
	}
	return nSys
}

// countClockBiasParameters calculates the number of clock bias parameters needed
func countClockBiasParameters(nSys [4]int) int {
	nClk := 0
	for _, v := range nSys {
		if v > 0 {
			nClk += 1
		}
	}
	return nClk
}

// isSppConverged checks if the SPP solution has converged
func isSppConverged(dx mat.Vector, threshold float64) bool {
	return math.Abs(dx.AtVec(0)) < threshold &&
		math.Abs(dx.AtVec(1)) < threshold &&
		math.Abs(dx.AtVec(2)) < threshold
}

// getMaxResidualAndSatellite finds the satellite with maximum residual
func getMaxResidualAndSatellite(dr mat.Vector, sats []SatType, excluded []SatType) (float64, SatType) {
	dmax := 0.0
	var dsat SatType

	i := 0
	for _, sat := range sats {
		if slices.Contains(excluded, sat) {
			continue
		}
		if math.Abs(dr.AtVec(i)) > dmax {
			dmax = math.Abs(dr.AtVec(i))
			dsat = sat
		}
		i++
	}
	return dmax, dsat
}

// updateUserPosition updates receiver position based on coordinate system
func updateUserPosition(upos *PosXYZ, dx mat.Vector, opt *SppOpt) {
	if opt.GByXyz {
		upos.X += dx.AtVec(0)
		upos.Y += dx.AtVec(1)
		upos.Z += dx.AtVec(2)
	} else {
		dxyz := NewPosENU(dx.AtVec(0), dx.AtVec(1), dx.AtVec(2)).ToXYZ(*upos)
		upos.X = dxyz.X
		upos.Y = dxyz.Y
		upos.Z = dxyz.Z
	}
}

// updateClockBias updates clock bias and inter-system bias parameters
func updateClockBias(clkb []float64, dx mat.Vector, nClk int) {
	for j := 0; j < nClk; j++ {
		clkb[j] += dx.AtVec(3 + j)
	}
}

// satClk calculates satellite clock correction
func satClk(e *Ephe, rcvt GTime, psr float64) (dts float64) {
	Mue := 3.986005e14 // Earth gravitational constant [m^3/s^2]
	sys := e.Sat.Sys()
	if sys == 'E' || sys == 'C' {
		Mue = 3.986004418e14
	}
	switch {
	case sys == 'G' || sys == 'J' || sys == 'E' || sys == 'C':
		// Relativistic correction
		tk := rcvt.ToTime().Sub(e.Toe.ToTime()).Seconds() - psr/C
		n := math.Sqrt(Mue)/e.SqrtA/e.SqrtA/e.SqrtA + e.DeltaN
		mk := e.M0 + n*tk
		ek := mk
		for i := 0; i < 10; i++ {
			ek = mk + e.Ecc*math.Sin(ek)
		}
		tr := -2 * math.Sqrt(Mue) / C / C * e.Ecc * e.SqrtA * math.Sin(ek)
		// Clock correction coefficients
		tk = rcvt.ToTime().Sub(e.Toc.ToTime()).Seconds() - psr/C
		dt := e.Af0 + e.Af1*tk + e.Af2*tk*tk
		// Group delay
		tg := e.Tgd
		if sys == 'E' {
			tg = e.Tgd2 // E1/E5b
		}
		dts = tr + dt - tg
	case sys == 'R':
		tk := rcvt.ToTime().Sub(e.Toe.ToTime()).Seconds() - psr/C // GLONASS uses Toe
		dt := -e.TauN + e.GammaN*tk
		dts = dt
	}
	return
}

// getWeight calculates observation weight based on elevation angle and system type
func getWeight(mode int, sat SatType, elv float64, eph *Ephe) (wg float64) {
	if elv > 0 {
		switch mode {
		case 0: // No weighting (equal weights)
			return 1.0
		case 1: // RTKLIB weighting
			fact := 1.0
			freq := L1
			switch sat.Sys() {
			case 'R':
				fact = 1.5
				freq = G1 + G1d*float64(eph.FreqN)
			case 'E':
				freq = E1
			case 'C':
				freq = B1
			}
			el := elv
			if el < ToRad(MIN_ELEVATION_FOR_WEIGHT) {
				el = ToRad(MIN_ELEVATION_FOR_WEIGHT)
			}
			varr := SQ(100) * (SQ(0.003) + SQ(0.003)/math.Sin(el))
			wg = SQ(fact) * varr
			wg += uraEphe(sat, eph.Sva) // SISA
			wg += SQ(0.3)               // code bias error std
			wg += SQ(5.0) * SQ(L1/freq) // ionospheric delay (L1) variance
			wg += SQ(3.0)               // tropospheric delay variance
			wg = 1.0 / wg
		case 2: // RTK Core weighting
			wg = ToDeg(elv) / 90.0
		case 3: // GPS Programming book weighting
			const VER_ZNH = 0.8 * 0.8
			wg = math.Sin(elv) * math.Sin(elv) / VER_ZNH
		}
		if wg < MIN_WEIGHT {
			wg = MIN_WEIGHT
		}
		return
	} else {
		return 1.0
	}
}

// uraEphe calculates User Range Accuracy (URA) variance for ephemeris
func uraEphe(sat SatType, ura int) float64 {
	uraVal := [...]float64{2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0}
	switch sat.Sys() {
	case 'E': // Galileo SIS ICD v2.1
		if ura <= 49 {
			return SQ(float64(ura) * 0.01)
		} else if ura <= 74 {
			return SQ(0.5 + (float64(ura)-50)*0.02)
		} else if ura <= 99 {
			return SQ(1.0 + (float64(ura)-75)*0.04)
		} else if ura <= 125 {
			return SQ(2.0 + (float64(ura)-100)*0.16)
		} else {
			return SQ(500.0)
		}
	case 'R':
		return SQ(5.0)
	default:
		if ura < 0 || ura > 14 {
			return SQ(6144.0)
		} else {
			return SQ(uraVal[ura])
		}
	}
}

// getGloFreq calculates GLONASS carrier frequencies based on frequency number
func getGloFreq(eph *Ephe, freq [NFREQ]float64) [NFREQ]float64 {
	freq2 := [NFREQ]float64{}
	for f := range NFREQ {
		switch freq[f] {
		case G1:
			freq2[f] = G1 + G1d*float64(eph.FreqN)
		case G2:
			freq2[f] = G2 + G2d*float64(eph.FreqN)
		default:
			freq2[f] = freq[f]
		}
	}
	return freq2
}
