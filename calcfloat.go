// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.15
//

// Implements float solution calculation for RTK (Real-Time Kinematic) positioning.

package gortk

import (
	"fmt"
	"math"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/mat"
)

// CalcFloat computes float solution from single point positioning results of rover and base stations
// It performs double-difference processing and Kalman filtering to estimate receiver position
// and float ambiguity values with high precision
//
// Parameters:
//   - rovSpp: Single point positioning solution for rover station
//   - baseSpp: Single point positioning solution for base station
//   - basePos: Known position of the base station
//   - opt: Float solution calculation options and parameters
//
// Returns:
//   - FloatSol: Float solution containing position and ambiguity estimates
//   - error: Any error encountered during processing
func CalcFloat(
	rovSpp *SppSol, // Single point positioning result for rover station
	baseSpp *SppSol, // Single point positioning result for base station
	basePos *PosXYZ, // Base station position
	opt *FloatOpt, // Calculation options
) (*FloatSol, error) {

	// Create satellite pairs for double-difference processing
	satPairF, satF, numSats, err := createSatPairs(rovSpp, baseSpp, opt)
	if err != nil {
		return nil, fmt.Errorf("createSatPairs() failed, err= %s", err.Error())
	}

	// Setup observation equations and initial state vector
	x, dy, ex, err := setupEquations(rovSpp, baseSpp, basePos, opt, satPairF, satF)
	if err != nil {
		return nil, fmt.Errorf("setupEquations() failed, err= %s", err.Error())
	}

	// Remove outlier observations to improve solution quality
	satPairF, satF, numSats, x, dy, err = removeOutliers(satPairF, satF, numSats, ex, x, dy, rovSpp, baseSpp, basePos, opt)
	if err != nil {
		return nil, fmt.Errorf("removeOutliers() failed, err= %s", err.Error())
	}

	// Solve using Kalman filter for optimal estimation
	finalX, finalP, age, err := solveWithKalmanFilter(x, dy, rovSpp, baseSpp, basePos, satPairF, satF, opt)
	if err != nil {
		return nil, fmt.Errorf("solveWithKalmanFilter() failed, err= %s", err.Error())
	}

	return setFloatSol(satPairF, satF, numSats, finalX, finalP, age, rovSpp), nil
}

// FloatOpt contains options and parameters for float solution calculation
// These parameters control the behavior of double-difference processing and Kalman filtering
type FloatOpt struct {
	NumFreq      int       // Number of frequencies to use in calculation
	WghMode      int       // Weighting mode for observation error covariance matrix: 0(OFF), 1(RTKLIB), 2(RTK Core)
	CnMask       float64   // Signal strength mask [dB]
	ElMask       float64   // Elevation mask [deg] (for debug info, not used in float calculation)
	StdCp        float64   // Carrier phase noise standard deviation
	StdPr        float64   // Pseudorange noise standard deviation
	VarP         float64   // Initial diagonal value for error covariance matrix in Kalman filter
	Instan       bool      // If true, calculate float solution independently for each epoch
	GJMix        bool      // If true, treat GPS and QZSS as the same satellite system
	MaxInnovPr   float64   // Threshold for pseudorange residual to detect outliers
	MaxInnovCp   float64   // Threshold for carrier phase residual to detect outliers
	ProcNoiseStd float64   // Process noise standard deviation for state estimation
	PrevFloatSol *FloatSol // Float solution from previous epoch for continuity
	PrevRovObsE  *ObsE     // Previous epoch rover observation data (for cycle slip detection)
	PrevBaseObsE *ObsE     // Previous epoch base observation data (for cycle slip detection)
	SkipTrop     bool      // If true, skip tropospheric correction
}

// FloatSol contains the results of float solution calculation
// It holds position estimates, float ambiguity values, and associated covariance information
type FloatSol struct {
	Pos      PosXYZ         // Float solution for receiver position
	Bias     []float64      // Float solution for single-difference integer ambiguities (rover - base)
	SatPairF []SatPairFType // List of satellite pairs used for double-difference processing
	SatF     []SatFType     // List of satellites and frequencies corresponding to ambiguity values
	P        mat.Matrix     // Error covariance matrix for position and ambiguity estimates
	NumSats  int            // Number of satellites used in calculation
	Age      float64        // Time difference between rover and base observations [sec] (rover - base)
	Time     GTime          // Rover observation time (after clock error correction)
}

// SatPairFType represents a satellite pair and frequency for double-difference processing
// It holds the reference satellite (highest elevation) and paired satellite information
type SatPairFType struct {
	S1 SatType // Reference satellite (highest elevation satellite)
	S2 SatType // Paired satellite for double-difference
	F  int     // Frequency index (0, 1, 2, 3)
}

// SatFType represents a satellite and frequency for single-difference ambiguity processing
// It holds satellite identification and frequency information
type SatFType struct {
	Sat SatType // Satellite identifier
	F   int     // Frequency index (0, ..., NFREQ-1)
}

// NewFloatOpt creates a new FloatOpt with default values
// Default values are tuned for typical RTK applications with good performance
func NewFloatOpt() *FloatOpt {
	return &FloatOpt{
		NumFreq:      1,       // Use single frequency (L1)
		WghMode:      1,       // RTKLIB weighting mode
		CnMask:       35,      // Signal strength threshold [dB]
		ElMask:       15,      // Elevation mask [deg]
		StdCp:        0.005,   // Carrier phase noise [m]
		StdPr:        0.5,     // Pseudorange noise [m]
		Instan:       false,   // Use epoch-to-epoch continuity
		GJMix:        false,   // Separate GPS and QZSS systems
		VarP:         30 * 30, // Initial position variance [m^2]
		MaxInnovPr:   30,      // Pseudorange outlier threshold [m]
		MaxInnovCp:   5,       // Carrier phase outlier threshold [m]
		ProcNoiseStd: 0.0001,  // Process noise [m/sqrt(s)]
		PrevRovObsE:  nil,     // No previous rover data
		PrevBaseObsE: nil,     // No previous base data
		PrevFloatSol: nil,     // No previous float solution
		SkipTrop:     false,   // Include tropospheric correction
	}
}

// NewFloatSol creates a new empty FloatSol structure
func NewFloatSol() *FloatSol {
	return &FloatSol{}
}

// createSatPairs creates satellite pairs for double-difference processing
// It selects reference satellites (highest elevation) for each system and forms pairs
func createSatPairs(
	rovSpp, baseSpp *SppSol,
	opt *FloatOpt,
) ([]SatPairFType, []SatFType, int, error) {

	// Extract highest elevation satellite for each satellite system
	maxElevSat, maxElev := getMaxElev(rovSpp, baseSpp, opt)

	// Create double-difference satellite pairs
	satPairF, satF, numSats := getSatPair(rovSpp, baseSpp, opt, maxElevSat)

	// Check minimum number of double-difference pairs required for solution
	if len(satPairF) < 3 {
		return nil, nil, 0, fmt.Errorf("number of dds: %d < %d", len(satPairF), 4)
	}

	if DBG_ >= 2 {
		PrintA("\trover sats:\n")
		for _, sat := range rovSpp.Sats {
			PrintA("\t%s: sx=%16.6f, sy=%16.6f, sz=%16.6f", sat, rovSpp.SatPos[sat].X, rovSpp.SatPos[sat].Y, rovSpp.SatPos[sat].Z)
			PrintA(", Elev:")
			if ToDeg(rovSpp.Elev[sat]) < opt.ElMask {
				PrintA("(%4.1f)", ToDeg(rovSpp.Elev[sat]))
			} else {
				PrintA(" %4.1f ", ToDeg(rovSpp.Elev[sat]))
			}
			PrintA(", C/N:")
			for f := range opt.NumFreq {
				if rovSpp.Sn[sat][f] < opt.CnMask {
					PrintA("(%4.1f)", rovSpp.Sn[sat][f])
				} else {
					PrintA(" %4.1f ", rovSpp.Sn[sat][f])
				}
			}
			PrintA(", LLI:")
			for f := range opt.NumFreq {
				if rovSpp.LLI[sat][f] != 0 {
					PrintA("(%d)", rovSpp.LLI[sat][f])
				} else {
					PrintA(" %d ", rovSpp.LLI[sat][f])
				}
			}
			PrintA(", Code:")
			for f := range opt.NumFreq {
				PrintA(" %s", rovSpp.Code[sat][f])
			}
			PrintA("\n")
		}
		PrintA("\tbase sats:\n")
		for _, sat := range baseSpp.Sats {
			PrintA("\t%s: sx=%16.6f, sy=%16.6f, sz=%16.6f", sat, baseSpp.SatPos[sat].X, baseSpp.SatPos[sat].Y, baseSpp.SatPos[sat].Z)
			PrintA(", Elev:")
			if ToDeg(baseSpp.Elev[sat]) < opt.ElMask {
				PrintA("(%4.1f)", ToDeg(baseSpp.Elev[sat]))
			} else {
				PrintA(" %4.1f ", ToDeg(baseSpp.Elev[sat]))
			}
			PrintA(", C/N:")
			for f := range opt.NumFreq {
				if baseSpp.Sn[sat][f] < opt.CnMask {
					PrintA("(%4.1f)", baseSpp.Sn[sat][f])
				} else {
					PrintA(" %4.1f ", baseSpp.Sn[sat][f])
				}
			}
			PrintA(", LLI:")
			for f := range opt.NumFreq {
				if baseSpp.LLI[sat][f] != 0 {
					PrintA("(%d)", baseSpp.LLI[sat][f])
				} else {
					PrintA(" %d ", baseSpp.LLI[sat][f])
				}
			}
			PrintA(", Code:")
			for f := range opt.NumFreq {
				PrintA(" %s", baseSpp.Code[sat][f])
			}
			PrintA("\n")
		}
		PrintA("\tmax-elev sats:\n")
		msats := []SatType{}
		for _, sat := range maxElevSat {
			if len(sat) > 0 {
				msats = append(msats, sat)
			}
		}
		for _, sat := range Sorted(msats) {
			PrintA("\t%c: %s (%.1f deg)\n", sat.Sys(), sat, ToDeg(maxElev[sat.Sys()]))
		}
		PrintA("\tdd pairs: (%d)", len(satPairF))
		var s1 SatType
		f := 0
		for _, sp := range satPairF {
			if s1 != sp.S1 || f != sp.F {
				PrintA("\n\t%s(%d):", sp.S1, sp.F+1)
				s1 = sp.S1
				f = sp.F
			}
			PrintA(" %s(%d)", sp.S2, sp.F+1)
		}
		PrintA("\n")
		PrintA("\tsd sats: (%d)", len(satF))
		var sys SysType
		f = -1
		for _, sf := range satF {
			if sys != sf.Sat.Sys() || f != sf.F {
				PrintA("\n\t%s(%d)", sf.Sat, sf.F+1)
			} else {
				PrintA(" %s(%d)", sf.Sat, sf.F+1)
			}
			sys = sf.Sat.Sys()
			f = sf.F
		}
		PrintA("\n")
	}

	return satPairF, satF, numSats, nil
}

// setupEquations sets up unknown vector and residual vector, detects outliers
// It initializes the state vector and computes initial residuals for Kalman filtering
func setupEquations(
	rovSpp, baseSpp *SppSol,
	basePos *PosXYZ,
	opt *FloatOpt,
	satPairF []SatPairFType,
	satF []SatFType,
) (*mat.VecDense, *mat.VecDense, []SatFType, error) {

	// Number of unknowns (3 for position + ambiguities)
	nx := 3 + len(satF)

	// Number of equations (2 per satellite pair: carrier phase + pseudorange)
	nv := len(satPairF) * 2

	// Setup unknown vector (position + ambiguities)
	x := makeX(nx, rovSpp, baseSpp, opt, satF)

	// Setup residual vector and detect outliers
	dy, ex := makeY(nv, x, basePos, rovSpp, baseSpp, opt, satPairF, satF)

	if DBG_ >= 2 {
		xyz := NewPosXYZ(x.AtVec(0), x.AtVec(1), x.AtVec(2))
		llh := xyz.ToLLH()
		llh.Lat = ToDeg(llh.Lat)
		llh.Lon = ToDeg(llh.Lon)
		PrintA("\tinit pos: LLH= %.8f %.8f %.3f, XYZ= %.3f %.3f %.3f\n", llh.Lat, llh.Lon, llh.Hei, xyz.X, xyz.Y, xyz.Z)
	}

	if DBG_ >= 4 {
		fmt.Printf("--- x(%d) ---\n", x.Len())
		for i := range x.Len() {
			fmt.Printf("%16.6f\n", x.AtVec(i))
		}
		fmt.Printf("--- dy(%d) ---\n", dy.Len())
		for i := range dy.Len() {
			fmt.Printf("%16.6f\n", dy.AtVec(i))
		}
	}

	return x, dy, ex, nil
}

// removeOutliers removes outlier observations and recalculates equations
// It filters out satellites with large residuals and updates the system accordingly
func removeOutliers(
	satPairF []SatPairFType,
	satF []SatFType,
	numSats int,
	ex []SatFType,
	x *mat.VecDense,
	dy *mat.VecDense,
	rovSpp, baseSpp *SppSol,
	basePos *PosXYZ,
	opt *FloatOpt,
) ([]SatPairFType, []SatFType, int, *mat.VecDense, *mat.VecDense, error) {

	if len(ex) > 0 {

		// Update satellite pairs by removing outliers
		newSatPairF, newSatF, newNumSats := renewSatPair(satPairF, satF, ex)
		if len(newSatPairF) < 4 {
			return nil, nil, 0, nil, nil, fmt.Errorf("number of dds: %d < %d", len(newSatPairF), 4)
		}

		// Recalculate unknown vector and residual vector with filtered satellites
		newNX := 3 + len(newSatF)
		newX := makeX(newNX, rovSpp, baseSpp, opt, newSatF)
		newNV := len(newSatPairF) * 2
		newDy, _ := makeY(newNV, newX, basePos, rovSpp, baseSpp, opt, newSatPairF, newSatF)

		if DBG_ >= 2 {
			PrintA("\t--- recalculated ---\n")
			PrintA("\tdd pairs: (%d)", len(newSatPairF))
			var s1 SatType
			f := 0
			for _, sp := range newSatPairF {
				if s1 != sp.S1 || f != sp.F {
					PrintA("\n\t%s(%d):", sp.S1, sp.F+1)
					s1 = sp.S1
					f = sp.F
				}
				PrintA(" %s(%d)", sp.S2, sp.F+1)
			}
			PrintA("\n")
			PrintA("\tsd sats: (%d)", len(newSatF))
			var sys SysType
			f = -1
			for _, sf := range newSatF {
				if sys != sf.Sat.Sys() || f != sf.F {
					PrintA("\n\t%s(%d)", sf.Sat, sf.F+1)
				} else {
					PrintA(" %s(%d)", sf.Sat, sf.F+1)
				}
				sys = sf.Sat.Sys()
				f = sf.F
			}
			PrintA("\n")
		}

		if DBG_ >= 4 {
			fmt.Printf("--- x(%d) ---\n", newX.Len())
			for i := range x.Len() {
				fmt.Printf("%16.6f\n", x.AtVec(i))
			}
			fmt.Printf("--- dy(%d) ---\n", newDy.Len())
			for i := range dy.Len() {
				fmt.Printf("%16.6f\n", dy.AtVec(i))
			}
		}

		return newSatPairF, newSatF, newNumSats, newX, newDy, nil
	}

	return satPairF, satF, numSats, x, dy, nil
}

// solveWithKalmanFilter solves float solution using Kalman filtering
// It performs optimal estimation of position and ambiguities using double-difference observations
func solveWithKalmanFilter(
	x *mat.VecDense,
	dy *mat.VecDense,
	rovSpp *SppSol,
	baseSpp *SppSol,
	basePos *PosXYZ,
	satPairF []SatPairFType,
	satF []SatFType,
	opt *FloatOpt,
) (*mat.VecDense, *mat.Dense, float64, error) {

	nx := x.Len()
	nv := dy.Len()

	// Create Jacobian matrix (design matrix)
	H := makeH(nv, nx, x, rovSpp, satPairF, satF)

	// Calculate time difference between rover and base observations
	age := rovSpp.Time.ToTime().Sub(baseSpp.Time.ToTime()).Seconds()
	PrintD(2, "\tage: %f\n", age)

	// Check if time difference exceeds maximum allowed value
	const MAX_AGE = 30.0
	if math.Round(age) > MAX_AGE {
		PrintD(2, "\tage exceeds max: %f(%f) > %f\n", math.Round(age), age, MAX_AGE)
		return nil, nil, age, fmt.Errorf("\tage exceeds max: %f(%f) > %f", math.Round(age), age, MAX_AGE)
	}

	// Create observation error covariance matrix
	R := makeR(rovSpp, opt, nv, satF, satPairF, age)

	// Create estimation error covariance matrix
	P := makeP(nx, opt.VarP, opt, satF, rovSpp, baseSpp)

	if DBG_ >= 4 {
		r, c := P.Dims()
		fmt.Printf("--- P(%d,%d) ---\n", r, c)
		for i := range r {
			for j := range c {
				fmt.Printf("%12.8f ", P.At(i, j))
			}
			fmt.Printf("\n")
		}
		r, c = H.Dims()
		fmt.Printf("--- H(%d,%d) ---\n", r, c)
		for i := range r {
			for j := range c {
				fmt.Printf("%12.8f ", H.At(i, j))
			}
			fmt.Printf("\n")
		}
		r, c = R.Dims()
		fmt.Printf("--- R(%d,%d) ---\n", r, c)
		for i := range r {
			for j := range c {
				fmt.Printf("%12.8f ", R.At(i, j))
			}
			fmt.Printf("\n")
		}
	}

	// Kalman filter computation
	currentX := x
	currentP := P
	currentDy := dy
	currentH := H
	nLoop := 1
	if opt.Instan {
		// In instantaneous mode, repeat observation updates to effectively nullify P matrix
		// This gives the same result as solving observation equations iteratively
		nLoop = 3
	}
	for loop := range nLoop {
		// Calculate Kalman gain
		K := makeK(currentP, currentH, R)

		if DBG_ >= 4 {
			if opt.Instan {
				PrintA("--- loop: %d ---\n", loop+1)
			}
			r, c := K.Dims()
			fmt.Printf("--- K(%d,%d) ---\n", r, c)
			for i := range r {
				for j := range c {
					fmt.Printf("%12.8f ", K.At(i, j))
				}
				fmt.Printf("\n")
			}
		}

		// Observation update: calculate float solution
		currentX, _ = updateX(currentX, K, currentDy)

		// Update covariance matrix
		currentP = updateP(K, currentH, currentP)

		if DBG_ >= 4 {
			fmt.Printf("--- x (after) ---\n")
			for i := range currentX.Len() {
				fmt.Printf("%16.6f\n", currentX.AtVec(i))
			}
			fmt.Printf("--- P (after) ---\n")
			r, c := currentP.Dims()
			for i := range r {
				for j := range c {
					fmt.Printf("%12.8f ", currentP.At(i, j))
				}
				fmt.Printf("\n")
			}
		}

		if DBG_ >= 2 {
			rovPos := NewPosXYZ(currentX.AtVec(0), currentX.AtVec(1), currentX.AtVec(2))
			llh := rovPos.ToLLH()
			llh.Lat = ToDeg(llh.Lat)
			llh.Lon = ToDeg(llh.Lon)
			msg := fmt.Sprintf("\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f, N=", loop+1, llh.Lat, llh.Lon, llh.Hei, rovPos.X, rovPos.Y, rovPos.Z)
			for j := 3; j < currentX.Len(); j++ {
				msg += fmt.Sprintf(" %.3f", currentX.AtVec(j))
			}
			msg += "\n"
			PrintA(msg)
		}

		// In instantaneous mode, update residuals and Jacobian
		if opt.Instan {
			currentDy, _ = makeY(nv, currentX, basePos, rovSpp, baseSpp, opt, satPairF, satF)
			currentH = makeH(nv, nx, currentX, rovSpp, satPairF, satF)
		}
	}

	return currentX, currentP, age, nil

}

// setFloatSol extracts position and integer ambiguities from unknown vector x
// It creates the final FloatSol structure with computed results
func setFloatSol(satPairF []SatPairFType, satF []SatFType, numSats int, x *mat.VecDense, P *mat.Dense, age float64, rovSpp *SppSol) *FloatSol {
	rslt := NewFloatSol()
	rslt.Time = rovSpp.Time
	rslt.SatPairF = satPairF
	rslt.SatF = satF
	rslt.NumSats = numSats
	rslt.P = P
	rslt.Age = age
	// Extract receiver position solution (XYZ coordinates)
	rslt.Pos.X = x.AtVec(0)
	rslt.Pos.Y = x.AtVec(1)
	rslt.Pos.Z = x.AtVec(2)
	// Extract integer ambiguity solutions
	// First 3 elements of x vector are position (X, Y, Z)
	// Elements 4 and beyond are single-difference integer ambiguities
	rslt.Bias = make([]float64, x.Len()-3)
	for k, j := 0, 3; j < x.Len(); j++ {
		rslt.Bias[k] = x.AtVec(j)
		k += 1
	}
	return rslt
}

// getMaxElev finds the highest elevation satellite for each satellite system
// It selects reference satellites based on elevation angle and data availability
func getMaxElev(rovSpp, baseSpp *SppSol, opt *FloatOpt) (map[SysType]SatType, map[SysType]float64) {
	// Count available observation data for each satellite and find maximum per system
	numF := map[SatType]int{}
	maxF := map[SysType]int{}
	if true {
		for _, sat := range rovSpp.Sats {
			if slices.Contains(baseSpp.Sats, sat) {
				c := 0
				for f := range opt.NumFreq {
					if isDataOK(rovSpp, baseSpp, sat, f, opt.CnMask) {
						c += 1
					}
				}
				numF[sat] = c
				sys := sat.Sys()
				if opt.GJMix && sys == 'J' {
					sys = 'G'
				}
				if c > maxF[sys] {
					maxF[sys] = c
				}
			}
		}
		// if DBG_ >= 4 {
		// 	PrintA("\tnumF:")
		// 	for _, sat := range rovSpp.Sats {
		// 		if slices.Contains(baseSpp.Sats, sat) {
		// 			PrintA(" %s:%d", sat, numF[sat])
		// 		}
		// 	}
		// 	PrintA("\n")
		// 	PrintA("\tmaxF:")
		// 	for _, sys := range []SysType{'G', 'E', 'R', 'C', 'S'} {
		// 		if maxF[sys] > 0 {
		// 			PrintA(" %c:%d", sys, maxF[sys])
		// 		}
		// 	}
		// 	PrintA("\n")
		// }
	}
	// Find highest elevation satellite for each system
	maxElevSat := map[SysType]SatType{}
	maxElev := map[SysType]float64{}
	for _, sat := range rovSpp.Sats {
		if slices.Contains(baseSpp.Sats, sat) {
			sys := sat.Sys()
			if opt.GJMix && sys == 'J' {
				sys = 'G'
			}
			el := rovSpp.Elev[sat]
			mx := maxElev[sys]
			if el > mx {
				if true {
					// Only select satellites with maximum number of available frequencies
					if numF[sat] == maxF[sys] {
						maxElevSat[sys] = sat
						maxElev[sys] = el
					}
				} else {
					maxElevSat[sys] = sat
					maxElev[sys] = el
				}
			}
		}
	}
	return maxElevSat, maxElev
}

// getSatPair generates list of double-difference satellite pairs
// It creates pairs with reference satellites and other satellites in the same system
func getSatPair(rovSpp, baseSpp *SppSol, opt *FloatOpt, maxElevSat map[SysType]SatType) ([]SatPairFType, []SatFType, int) {
	sp := []SatPairFType{}
	sf := []SatFType{}
	ns := 0 // Count number of satellites
	for f := range opt.NumFreq {
		for _, sat := range rovSpp.Sats {
			if slices.Contains(baseSpp.Sats, sat) { // Satellite must exist in both rover and base data
				sys := sat.Sys()
				if opt.GJMix && sys == 'J' {
					sys = 'G'
				}
				// Extract satellite pairs with both pseudorange and carrier phase data meeting signal strength criteria
				msat := maxElevSat[sys]
				if sat != msat {
					if isDataOK(rovSpp, baseSpp, sat, f, opt.CnMask) && isDataOK(rovSpp, baseSpp, msat, f, opt.CnMask) {
						sp = append(sp, SatPairFType{S1: msat, S2: sat, F: f})
						sf = append(sf, SatFType{Sat: sat, F: f})
						if f == 0 {
							ns += 1
						}
					}
				} else {
					if isDataOK(rovSpp, baseSpp, msat, f, opt.CnMask) {
						sf = append(sf, SatFType{Sat: msat, F: f})
						if f == 0 {
							ns += 1
						}
					}
				}
			}
		}
	}
	return sp, sf, ns
}

// isDataOK checks if observation data is valid for both rover and base stations
// It verifies pseudorange, carrier phase, and signal strength requirements
func isDataOK(rovSpp, baseSpp *SppSol, sat SatType, f int, cnMask float64) bool {
	return rovSpp.Pr[sat][f] != 0 && rovSpp.Cp[sat][f] != 0 && rovSpp.Sn[sat][f] >= cnMask && baseSpp.Pr[sat][f] != 0 && baseSpp.Cp[sat][f] != 0 && baseSpp.Sn[sat][f] >= cnMask
}

// renewSatPair removes satellites from exclusion list from satellite pair and frequency lists
// It filters out outlier satellites and updates the satellite configuration
func renewSatPair(satPairF []SatPairFType, satF []SatFType, exSats []SatFType) (sp []SatPairFType, sf []SatFType, ns int) {
	// Remove excluded satellites from frequency list
	for _, s := range satF {
		if !slices.Contains(exSats, s) {
			sf = append(sf, s)
			if s.F == 0 {
				ns += 1
			}
		}
	}
	// Remove satellite pairs containing excluded satellites
	for _, s := range satPairF {
		if !slices.Contains(exSats, SatFType{Sat: s.S2, F: s.F}) {
			sp = append(sp, s)
		}
	}
	return
}

// makeX sets up the unknown vector for Kalman filtering
// It initializes position and ambiguity estimates with appropriate starting values
func makeX(nx int, rovSpp, baseSpp *SppSol, opt *FloatOpt, sats []SatFType) *mat.VecDense {
	x := mat.NewVecDense(nx, nil)
	// Set first 3 elements to rover single point positioning results
	x.SetVec(0, rovSpp.Pos.X)
	x.SetVec(1, rovSpp.Pos.Y)
	x.SetVec(2, rovSpp.Pos.Z)
	// Set remaining elements to initial estimates of single-difference integer ambiguities
	for i, sf := range sats {
		sat := sf.Sat
		f := sf.F
		lmu := C / rovSpp.Freq[sat][f]
		lmr := C / baseSpp.Freq[sat][f]
		// Initial estimate of single-difference integer ambiguity
		b := (rovSpp.Cp[sat][f] - baseSpp.Cp[sat][f]) - (rovSpp.Pr[sat][f]/lmu - baseSpp.Pr[sat][f]/lmr)
		// Use previous epoch float solution ambiguity values if available
		if opt.PrevFloatSol != nil {
			// Only if no cycle slip detected
			if !isSlipDetected(sat, f, rovSpp, opt.PrevRovObsE) && !isSlipDetected(sat, f, baseSpp, opt.PrevBaseObsE) {
				l := slices.Index(opt.PrevFloatSol.SatF, sf)
				if l >= 0 {
					b = opt.PrevFloatSol.Bias[l]
				}
			}
		}
		x.SetVec(3+i, b)
	}
	return x
}

// makeP sets up the estimation error covariance matrix for Kalman filtering
// It initializes the covariance matrix with appropriate values and handles continuity
func makeP(nx int, varP float64, opt *FloatOpt, sats []SatFType, rovSpp, baseSpp *SppSol) (P *mat.Dense) {
	P = mat.NewDense(nx, nx, nil)
	// Set diagonal elements
	for j := 0; j < nx; j++ {
		P.Set(j, j, varP)
	}
	// If instantaneous mode, return here
	if opt.Instan {
		return
	}
	// Use previous epoch float solution error covariance matrix values if available
	if opt.PrevFloatSol != nil {
		pSol := opt.PrevFloatSol
		for j := 3; j < nx; j++ {
			sat := sats[j-3].Sat
			f := sats[j-3].F
			// Only if no cycle slip detected
			if !isSlipDetected(sat, f, rovSpp, opt.PrevRovObsE) && !isSlipDetected(sat, f, baseSpp, opt.PrevBaseObsE) {
				i := slices.Index(pSol.SatF, sats[j-3])
				if i >= 0 {
					dt := rovSpp.Time.ToTime().Sub(opt.PrevFloatSol.Time.ToTime()).Seconds() // Time elapsed since previous epoch
					vs := math.Abs(dt) * opt.ProcNoiseStd * opt.ProcNoiseStd                 // Variance increase due to time elapsed
					P.Set(j, j, pSol.P.At(3+i, 3+i)+vs)                                      // Add to diagonal only
					for k := 3; k < j; k++ {
						l := slices.Index(pSol.SatF, sats[k-3])
						if l >= 0 {
							P.Set(k, j, pSol.P.At(3+l, 3+i))
							P.Set(j, k, pSol.P.At(3+l, 3+i))
						}
					}
				}
			}
		}
		// Reset covariance for satellites with cycle slips
		for j := 3; j < nx; j++ {
			sat := sats[j-3].Sat
			f := sats[j-3].F
			if isSlipDetected(sat, f, rovSpp, opt.PrevRovObsE) || isSlipDetected(sat, f, baseSpp, opt.PrevBaseObsE) {
				P.Set(j, j, varP)
				for k := 3; k < nx; k++ {
					if k != j {
						P.Set(k, j, 0)
						P.Set(j, k, 0)
					}
				}
			}
		}
	}
	return
}

// makeY calculates double-differences and creates residual vector
// It computes observation residuals and detects outliers for quality control
func makeY(nv int, x *mat.VecDense, basePos *PosXYZ, rovSpp, baseSpp *SppSol, opt *FloatOpt, satPair []SatPairFType, sf []SatFType) (*mat.VecDense, []SatFType) {
	dy := mat.NewVecDense(nv, nil)
	rovPos := NewPosXYZ(x.AtVec(0), x.AtVec(1), x.AtVec(2))
	// Calculate zenith hydrostatic delay for tropospheric correction
	zhdu := TropModel(rovPos)
	zhdr := TropModel(basePos)
	exSats := []SatFType{}
	for i, sp := range satPair {
		k := sp.S1                    // Reference satellite
		j := sp.S2                    // Paired satellite for double-difference
		f := sp.F                     // Carrier frequency index
		Xkr := baseSpp.SatPos[k]      // Reference satellite position (from base SPP)
		Xjr := baseSpp.SatPos[j]      // Paired satellite position (from base SPP)
		Xku := rovSpp.SatPos[k]       // Reference satellite position (from rover SPP)
		Xju := rovSpp.SatPos[j]       // Paired satellite position (from rover SPP)
		Rkr := EucDist(&Xkr, basePos) // Distance from base to reference satellite [m]
		Rjr := EucDist(&Xjr, basePos) // Distance from base to paired satellite [m]
		Rku := EucDist(&Xku, rovPos)  // Distance from rover to reference satellite [m]
		Rju := EucDist(&Xju, rovPos)  // Distance from rover to paired satellite [m]
		Pkr := baseSpp.Pr[k][f]       // Pseudorange from base to reference satellite [m]
		Pjr := baseSpp.Pr[j][f]       // Pseudorange from base to paired satellite [m]
		Pku := rovSpp.Pr[k][f]        // Pseudorange from rover to reference satellite [m]
		Pju := rovSpp.Pr[j][f]        // Pseudorange from rover to paired satellite [m]
		Ckr := baseSpp.Cp[k][f]       // Carrier phase from base to reference satellite [cycle]
		Cjr := baseSpp.Cp[j][f]       // Carrier phase from base to paired satellite [cycle]
		Cku := rovSpp.Cp[k][f]        // Carrier phase from rover to reference satellite [cycle]
		Cju := rovSpp.Cp[j][f]        // Carrier phase from rover to paired satellite [cycle]
		// Add satellite clock error correction (needed when observation times differ)
		Pkr += C * baseSpp.SatClk[k]
		Pjr += C * baseSpp.SatClk[j]
		Pku += C * rovSpp.SatClk[k]
		Pju += C * rovSpp.SatClk[j]
		lmkr := C / baseSpp.Freq[k][f] // Reference satellite wavelength at base [m]
		lmjr := C / baseSpp.Freq[j][f] // Paired satellite wavelength at base [m]
		lmku := C / rovSpp.Freq[k][f]  // Reference satellite wavelength at rover [m]
		lmju := C / rovSpp.Freq[j][f]  // Paired satellite wavelength at rover [m]
		// Add satellite clock error correction for carrier phase (needed when observation times differ)
		Ckr += C * baseSpp.SatClk[k] / lmkr
		Cjr += C * baseSpp.SatClk[j] / lmjr
		Cku += C * rovSpp.SatClk[k] / lmku
		Cju += C * rovSpp.SatClk[j] / lmju
		// GLONASS specific correction: receiver clock error doesn't cancel in double-difference
		if k.Sys() == 'R' {
			Ckr -= baseSpp.Clk[0] / lmkr
			Cjr -= baseSpp.Clk[0] / lmjr
			Cku -= rovSpp.Clk[0] / lmku
			Cju -= rovSpp.Clk[0] / lmju
		}
		Cdd := (Cku*lmku - Ckr*lmkr) - (Cju*lmju - Cjr*lmjr) // Carrier phase double-difference [m]
		Pdd := (Pku - Pkr) - (Pju - Pjr)                     // Pseudorange double-difference [m]
		Rdd := (Rku - Rkr) - (Rju - Rjr)                     // Geometric distance double-difference [m]
		ki := slices.Index(sf, SatFType{Sat: k, F: f})
		ji := slices.Index(sf, SatFType{Sat: j, F: f})
		Bkd := x.AtVec(3 + ki)     // Reference satellite single-difference integer ambiguity
		Bjd := x.AtVec(3 + ji)     // Paired satellite single-difference integer ambiguity
		Bdd := Bkd*lmku - Bjd*lmju // Integer ambiguity double-difference [m] (wavelength separation)
		ddCp := Cdd - Rdd - Bdd    // Carrier phase double-difference residual
		ddPr := Pdd - Rdd          // Pseudorange double-difference residual
		// Apply tropospheric correction if enabled
		if !opt.SkipTrop {
			Tkr := TropMapf(&baseSpp.Time, basePos, baseSpp.Elev[k]) * zhdr
			Tjr := TropMapf(&baseSpp.Time, basePos, baseSpp.Elev[j]) * zhdr
			Tku := TropMapf(&rovSpp.Time, rovPos, rovSpp.Elev[k]) * zhdu
			Tju := TropMapf(&rovSpp.Time, rovPos, rovSpp.Elev[j]) * zhdu
			Tdd := (Tku - Tkr) - (Tju - Tjr)
			ddCp -= Tdd
			ddPr -= Tdd
		}
		// Check for outliers and remove if necessary
		if math.Abs(ddCp) > opt.MaxInnovCp || math.Abs(ddPr) > opt.MaxInnovPr {
			exSats = append(exSats, SatFType{Sat: sp.S2, F: sp.F})
			PrintD(2, "\tOutlier: sat=%s, f=%d, ddCp=%f, ddPr=%f (thres: %f, %f)\n", sp.S2, sp.F, ddCp, ddPr, opt.MaxInnovCp, opt.MaxInnovPr)
		} else {
			dy.SetVec(i, ddCp)      // Set carrier phase residual
			dy.SetVec(nv/2+i, ddPr) // Set pseudorange residual
		}
	}
	return dy, exSats
}

// makeH creates the Jacobian matrix (H = dh(x)/dx) for Kalman filtering
// It computes partial derivatives of observations with respect to unknowns
func makeH(nv, nx int, x *mat.VecDense, rovSpp *SppSol, satPair []SatPairFType, sf []SatFType) *mat.Dense {
	H := mat.NewDense(nv, nx, nil)
	rovPos := NewPosXYZ(x.AtVec(0), x.AtVec(1), x.AtVec(2))
	for i, sp := range satPair {
		k := sp.S1              // Reference satellite
		j := sp.S2              // Paired satellite for double-difference
		f := sp.F               // Carrier frequency index
		Xku := rovSpp.SatPos[k] // Reference satellite position (from rover SPP)
		Xju := rovSpp.SatPos[j] // Paired satellite position (from rover SPP)
		// Calculate partial derivatives with respect to position
		Ex := -(DistDx(&Xju, rovPos) - DistDx(&Xku, rovPos))
		Ey := -(DistDy(&Xju, rovPos) - DistDy(&Xku, rovPos))
		Ez := -(DistDz(&Xju, rovPos) - DistDz(&Xku, rovPos))
		// Set position partial derivatives for carrier phase and pseudorange
		H.Set(i, 0, Ex)
		H.Set(i, 1, Ey)
		H.Set(i, 2, Ez)
		H.Set(nv/2+i, 0, Ex)
		H.Set(nv/2+i, 1, Ey)
		H.Set(nv/2+i, 2, Ez)
		lmku := C / rovSpp.Freq[k][f] // Reference satellite carrier wavelength
		lmju := C / rovSpp.Freq[j][f] // Paired satellite carrier wavelength (differs for GLONASS)
		ki := slices.Index(sf, SatFType{Sat: k, F: f})
		ji := slices.Index(sf, SatFType{Sat: j, F: f})
		// Set ambiguity partial derivatives (only for carrier phase)
		H.Set(i, 3+ki, lmku)
		H.Set(i, 3+ji, -lmju)
	}
	return H
}

// makeR creates the observation error covariance matrix
// It applies different weighting schemes based on elevation angle and system type
func makeR(rovSpp *SppSol, opt *FloatOpt, nv int, sats []SatFType, satPair []SatPairFType, age float64) *mat.Dense {
	R := mat.NewDense(nv, nv, nil)
	w_cp := 4.0  // Carrier phase observation error (diagonal)
	w_pr := 4.0  // Pseudorange observation error (diagonal)
	w_cp2 := 2.0 // Carrier phase observation error (off-diagonal)
	w_pr2 := 2.0 // Pseudorange observation error (off-diagonal)
	// Create column indices for setting off-diagonal terms
	ki := []int{0}
	sys0 := sats[0].Sat.Sys()
	if opt.GJMix && sys0 == 'J' {
		sys0 = 'G'
	}
	f0 := sats[0].F
	for i, sf := range sats {
		sys := sf.Sat.Sys()
		if opt.GJMix && sys == 'J' {
			sys = 'G'
		}
		if sys != sys0 || sf.F != f0 {
			ki = append(ki, i)
		}
		sys0 = sys
		f0 = sf.F
	}
	ki = append(ki, len(sats))
	// Calculate start and end indices for off-diagonal terms
	cs, ce := []int{}, []int{}
	for i := range len(ki) - 1 {
		n := ki[i+1] - ki[i] - 1
		l1 := ki[i] - i
		l2 := l1 + n
		for range n {
			cs = append(cs, l1)
			ce = append(ce, l2)
		}
	}
	PrintD(5, "\rcs=%v, len=%d\n", cs, len(cs))
	PrintD(5, "\rce=%v, len=%d\n", ce, len(ce))
	for i, sp := range satPair {
		k := sp.S1                          // Reference satellite
		j := sp.S2                          // Paired satellite for double-difference
		sinel_k := math.Sin(rovSpp.Elev[k]) // Elevation-dependent weight for reference satellite
		sinel_j := math.Sin(rovSpp.Elev[j]) // Elevation-dependent weight for paired satellite
		switch opt.WghMode {
		case 1:
			// RTKLIB weighting scheme
			a_ := 0.003           // Base term
			b_ := 0.003           // Elevation term
			c_ := 0.0             // Baseline length term (not used)
			d_ := C * 5e-12 * age // Clock error term
			if j[0] == 'R' {
				a_ = 0.0045
				b_ = 0.0045
			}
			w_cp_k := 2.0*(a_*a_+b_*b_/sinel_k/sinel_k+c_*c_) + d_*d_ // Single-difference observation error
			w_cp_j := 2.0*(a_*a_+b_*b_/sinel_j/sinel_j+c_*c_) + d_*d_
			a_ = 0.3
			b_ = 0.3
			if j[0] == 'R' {
				a_ = 0.45
				b_ = 0.45
			}
			w_pr_k := 2.0*(a_*a_+b_*b_/sinel_k/sinel_k+c_*c_) + d_*d_
			w_pr_j := 2.0*(a_*a_+b_*b_/sinel_j/sinel_j+c_*c_) + d_*d_
			w_cp = w_cp_k + w_cp_j // Double-difference observation error
			w_pr = w_pr_k + w_pr_j
			w_cp2 = w_cp_k // Off-diagonal terms use reference satellite values
			w_pr2 = w_pr_k
		case 2:
			// RTK Core weighting scheme
			w_cp = 4 * (opt.StdCp + opt.StdCp/sinel_j) * (opt.StdCp + opt.StdCp/sinel_j)
			w_pr = 4 * (opt.StdPr + opt.StdPr/sinel_j) * (opt.StdPr + opt.StdPr/sinel_j)
			w_cp2 = 2 * opt.StdCp * opt.StdCp
			w_pr2 = 2 * opt.StdPr * opt.StdPr
		}
		// Set off-diagonal terms
		for m := cs[i]; m < ce[i]; m++ {
			R.Set(i, m, w_cp2)
			R.Set(nv/2+i, nv/2+m, w_pr2)
		}
		// Set diagonal terms
		R.Set(i, i, w_cp)
		R.Set(nv/2+i, nv/2+i, w_pr)
	}
	return R
}

// makeK calculates Kalman gain K = P H^T (H P H^T + R)^-1
// This is the core computation for optimal state estimation
func makeK(P, H, R *mat.Dense) *mat.Dense {
	var A, B, C, D, K mat.Dense
	A.Mul(H, P)
	B.Mul(&A, H.T())
	C.Add(&B, R)
	C.Inverse(&C)
	D.Mul(P, H.T())
	K.Mul(&D, &C)
	return &K
}

// updateX calculates x' = x + dx, where dx = K dy
// This performs the state update in Kalman filtering
func updateX(x *mat.VecDense, K *mat.Dense, dy *mat.VecDense) (*mat.VecDense, mat.Vector) {
	var A mat.Dense
	A.Mul(K, dy)
	nx, _ := A.Dims()
	dx := mat.NewVecDense(nx, nil)
	for j := range nx {
		dx.SetVec(j, A.At(j, 0))
	}
	var x2 mat.VecDense
	x2.AddVec(x, dx)
	return &x2, dx
}

// updateP calculates P' = (I - K H) P
// This updates the error covariance matrix after observation update
func updateP(K, H, P *mat.Dense) *mat.Dense {
	nx, _ := K.Dims()
	I := mat.NewDense(nx, nx, nil)
	for j := range nx {
		for k := range nx {
			if j == k {
				I.Set(j, k, 1)
			}
		}
	}
	var A, B, C mat.Dense
	A.Mul(K, H)
	B.Sub(I, &A)
	C.Mul(&B, P)
	return &C
}

// isSlipDetected checks for cycle slip detection using LLI (Loss of Lock Indicator)
// It detects both half-cycle slips and full cycle slips using different LLI bits
func isSlipDetected(sat SatType, f int, currSol *SppSol, prevObsE *ObsE) bool {
	// Check LLI Bit 0 for half-cycle slip detection
	if currSol != nil && currSol.LLI[sat][f]&1 == 1 {
		PrintD(4, "slip detected by LLI, sat=%s, f=%d, LLI=%d\n", sat, f, currSol.LLI[sat][f])
		return true
	}
	// Check LLI Bit 1 change for full cycle slip detection
	if prevObsE != nil {
		if obss, ok := prevObsE.DatS[sat]; ok {
			prev := (obss.LLI[f]&2 == 2)
			curr := (currSol.LLI[sat][f]&2 == 2)
			if (prev && !curr) || (!prev && curr) {
				PrintD(4, "slip detected by LLI change, sat=%s, f=%d, prevLLI=%d, currLLI=%d\n", sat, f, obss.LLI[f], currSol.LLI[sat][f])
				return true
			}
		}
	}
	return false
}
