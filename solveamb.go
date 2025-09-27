// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.21
//

// Implements integer ambiguity resolution for RTK (Real-Time Kinematic) positioning.

package gortk

import (
	"fmt"
	"math"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/mat"
)

// SolveAmb performs integer ambiguity resolution for RTK positioning.
// It takes float solution results and attempts to fix carrier phase ambiguities to integer values
// using the LAMBDA method with various retry strategies.
//
// Parameters:
//   - rovSpp: Single Point Positioning solution for rover station
//   - baseSpp: Single Point Positioning solution for base station
//   - floatSol: Float solution results containing ambiguity estimates
//   - basePos: Known position of the base station
//   - opt: Ambiguity resolution options and parameters
//
// Returns:
//   - AmbSol: Ambiguity solution with fixed integer values and position
//   - error: Any error encountered during processing
func SolveAmb(
	rovSpp *SppSol, // Single point positioning solution for rover station
	baseSpp *SppSol, // Single point positioning solution for base station
	floatSol *FloatSol, // Float solution results
	basePos *PosXYZ, // Base station position
	opt *AmbOpt, // Calculation options
) (*AmbSol, error) {

	if DBG_ >= 3 && opt.PrevAmbSol != nil {
		PrintA("\tdd pairs (prev epoch): (%d), ratio=%8.2f", len(opt.PrevAmbSol.DDPairs), opt.PrevAmbSol.Ratio)
		printDDpairs(opt.PrevAmbSol.DDPairs, rovSpp)
	}

	// Select double-difference pairs for ambiguity resolution
	ddPairs := selectDDpairs(rovSpp, baseSpp, floatSol, opt)

	// Check if we have enough double-difference pairs (minimum 3 required)
	if len(ddPairs) < 3 {
		return nil, fmt.Errorf("not enough dd pairs, num of dd pairs=%d", len(ddPairs))
	}

	// Solve integer ambiguities using LAMBDA method
	ratio, b, bf, res, err := solveAmbByLambda(ddPairs, floatSol)
	if err != nil {
		return nil, fmt.Errorf("solveAmbByLambda() failed, err= %s", err.Error())
	}
	PrintD(2, "\tratio: %8.2f (%.3f/%.3f)\n", ratio, res[1], res[0])

	// Validate solution and perform retry strategies if needed
	finalRatio, finalB, finalBf, finalRes, finalDDPairs, err := validateAndRetry(ratio, b, bf, res, ddPairs, rovSpp, floatSol, opt)
	if err != nil {
		return nil, fmt.Errorf("validateAndRetry() failed, err= %s", err.Error())
	}

	// If retry strategies fail to achieve fix, return float solution
	if finalRatio < opt.RatioThres {
		PrintD(2, "\tratio: %8.2f < %3.2f\n", finalRatio, opt.RatioThres)
		return setAmbSol(finalRatio, finalB, finalBf, finalRes, finalDDPairs, floatSol.Pos), nil
	}

	// Recalculate receiver position using fixed integer ambiguities
	finalPos, err := recalcPos(finalDDPairs, finalB, rovSpp, baseSpp, basePos, floatSol.Pos, opt)
	if err != nil {
		return nil, fmt.Errorf("recalcPos() failed, err= %s", err.Error())
	}

	return setAmbSol(finalRatio, finalB, finalBf, finalRes, finalDDPairs, finalPos), nil
}

// Constants for ambiguity resolution and position calculation
const (
	MAX_LOOP          = 10    // Maximum number of iterations for position convergence
	CONVERGENCE_THRES = 0.001 // Convergence threshold for position updates (meters)
)

// AmbOpt contains options and parameters for ambiguity resolution
// These parameters control the behavior of the LAMBDA method and retry strategies
type AmbOpt struct {
	RatioThres      float64 // Threshold for ratio test to determine successful fix
	PrevAmbSol      *AmbSol // Ambiguity solution from previous epoch for continuity
	MaxRetNum       int     // Maximum number of retry attempts
	SkipTrop        bool    // Skip tropospheric correction if true
	MaxElevToRemove float64 // Maximum elevation angle for satellite removal during retry (degrees)
}

// AmbSol contains the results of integer ambiguity resolution
// This structure holds both the fixed solution and diagnostic information
type AmbSol struct {
	Pos     PosXYZ         // Final position solution with fixed ambiguities
	Ratio   float64        // Ratio test value (second best / best solution)
	B       []float64      // Fixed integer ambiguity values
	Bf      []float64      // Float ambiguity values before fixing
	Res     []float64      // Residuals used for ratio test calculation
	DDPairs []SatPairFType // Double-difference pairs used in ambiguity resolution
}

// NewAmbOpt creates a new AmbOpt with default values
// Default values are tuned for typical RTK applications
func NewAmbOpt() *AmbOpt {
	return &AmbOpt{
		RatioThres:      3.0,   // Standard ratio threshold for fix acceptance
		MaxRetNum:       10,    // Allow up to 10 retry attempts
		SkipTrop:        false, // Include tropospheric correction by default
		MaxElevToRemove: 45,    // Remove satellites below 45 degrees elevation
	}
}

// NewAmbSol creates a new empty AmbSol structure
func NewAmbSol() *AmbSol {
	return &AmbSol{}
}

// selectDDpairs selects double-difference pairs for ambiguity resolution
// It filters out pairs with cycle slips and applies continuity constraints
func selectDDpairs(rovSpp, baseSpp *SppSol, floatSol *FloatSol, opt *AmbOpt) []SatPairFType {
	ddPairs := []SatPairFType{}
	for _, sp := range floatSol.SatPairF {
		// For first epoch, select all available pairs unconditionally
		if opt.PrevAmbSol == nil {
			ddPairs = append(ddPairs, sp)
			continue
		}
		// For subsequent epochs, only select pairs without cycle slips
		// Check both rover and base station for cycle slip indicators
		// LLI bit 2: cycle slip detected, LLI bit 1: half-cycle slip
		j := sp.S2
		f := sp.F
		if rovSpp.LLI[j][f]&2 != 2 && baseSpp.LLI[j][f]&2 != 2 && rovSpp.LLI[j][f]&1 != 1 && baseSpp.LLI[j][f]&1 != 1 {
			ddPairs = append(ddPairs, sp)
		}
	}
	if DBG_ >= 2 {
		PrintA("\tdd pairs: (%d)", len(ddPairs))
		printDDpairs(ddPairs, rovSpp)
	}
	return ddPairs
}

// solveAmbByLambda solves integer ambiguities using the LAMBDA method
// It performs the core ambiguity resolution algorithm and returns the ratio test result
func solveAmbByLambda(ddPairs []SatPairFType, floatSol *FloatSol) (ratio float64, B []float64, Bf []float64, res []float64, err error) {

	// Calculate double-difference ambiguities and covariance matrix
	Bdd, Pdd := calcDD(ddPairs, floatSol)

	// Number of double-difference pairs
	ndd := len(ddPairs)

	if DBG_ >= 4 {
		fmt.Printf("--- Bdd(%d) ---\n", ndd)
		for i := range ndd {
			fmt.Printf("%12.6f\n", Bdd[i])
		}
		fmt.Printf("--- Pdd(%d,%d) ---\n", ndd, ndd)
		for i := range ndd {
			for j := range ndd {
				fmt.Printf("%12.6f ", Pdd[i*ndd+j])
			}
			PrintA("\n")
		}
		fmt.Printf("--- Pdd(diag) ---\n")
		for i := range ndd {
			fmt.Printf("%12.6f\n", Pdd[i*ndd+i])
		}
	}

	// Solve ambiguities using LAMBDA method
	B = make([]float64, ndd*2) // Array to store fixed integer ambiguities
	res = make([]float64, 2)   // Array for ratio test values
	err = LAMBDA(ndd, 2, Bdd, Pdd, B, res)
	if err != nil {
		return 0, nil, nil, nil, fmt.Errorf("LAMBDA() failed, err= %s", err.Error())
	}
	Bf = Bdd

	// Calculate ratio test value (second best / best solution)
	if res[0] > 0 {
		ratio = res[1] / res[0]
	} else {
		return 0, nil, nil, nil, fmt.Errorf("invalid results from LAMBDA(), res[0]= %f", res[0])
	}

	return
}

// calcDD calculates double-difference ambiguities and covariance matrix
// It extracts the ambiguity part from the float solution covariance matrix
func calcDD(ddPairs []SatPairFType, floatSol *FloatSol) (Bdd, Pdd []float64) {

	// Extract ambiguity covariance matrix from float solution (remove first 3x3 position block)
	r, c := floatSol.P.Dims()
	cov := make([]float64, (r-3)*(c-3))
	for i, j := 0, 3; j < r; j++ {
		for k := 3; k < c; k++ {
			cov[i] = floatSol.P.At(j, k)
			i += 1
		}
	}
	// Calculate double-difference for each row of the covariance matrix
	tmp := []float64{}
	n := len(floatSol.SatF)
	for i := range n {
		for _, sp := range floatSol.SatPairF {
			if slices.Contains(ddPairs, sp) {
				// Find indices for satellites in the pair
				ki := slices.Index(floatSol.SatF, SatFType{Sat: sp.S1, F: sp.F})
				ji := slices.Index(floatSol.SatF, SatFType{Sat: sp.S2, F: sp.F})
				// Calculate single-difference covariance
				tmp = append(tmp, cov[i*n+ki]-cov[i*n+ji])
			}
		}
		// Alternative calculation order (commented out)
		// for j := range len(floatSol.SatPairF) {
		// 	tmp = append(tmp, cov[i*n]-cov[i*n+j+1]) // This order gives same result
		// }
	}
	// Calculate double-difference for each column to form final covariance matrix
	ndd := len(ddPairs)
	Pdd = make([]float64, ndd*ndd)
	for j := range ndd {
		for k, sp := range ddPairs {
			// Find indices for satellites in the pair
			ki := slices.Index(floatSol.SatF, SatFType{Sat: sp.S1, F: sp.F})
			ji := slices.Index(floatSol.SatF, SatFType{Sat: sp.S2, F: sp.F})
			// Calculate double-difference covariance
			Pdd[k*ndd+j] = tmp[ki*ndd+j] - tmp[ji*ndd+j]
		}
		// Alternative calculation order (commented out)
		// for k := range ndd {
		// 	Pdd[k*ndd+j] = tmp[j] - tmp[(k+1)*ndd+j] // This order gives same result
		// }
	}

	// Create double-difference float ambiguities
	Bdd = []float64{}
	for _, sp := range ddPairs {
		// Find indices for satellites in the pair
		ki := slices.Index(floatSol.SatF, SatFType{Sat: sp.S1, F: sp.F})
		ji := slices.Index(floatSol.SatF, SatFType{Sat: sp.S2, F: sp.F})
		// Calculate double-difference float ambiguity
		Bdd = append(Bdd, floatSol.Bias[ki]-floatSol.Bias[ji])
	}
	// Alternative calculation order (commented out)
	// for i := range len(floatSol.SatPairF) {
	// 	Bdd = append(Bdd, floatSol.Bias[0]-floatSol.Bias[i+1]) // This order gives same result
	// }
	return
}

// printDDpairs prints double-difference pairs in a formatted way
// Shows satellite pairs with their frequency and elevation angles
func printDDpairs(ddPairs []SatPairFType, rovSpp *SppSol) {
	var s1 SatType
	f := 0
	for _, sp := range ddPairs {
		if s1 != sp.S1 || f != sp.F {
			PrintA("\n\t%s(%d,%.1f):", sp.S1, sp.F+1, ToDeg(rovSpp.Elev[sp.S1]))
			s1 = sp.S1
			f = sp.F
		}
		PrintA(" %s(%d,%.1f)", sp.S2, sp.F+1, ToDeg(rovSpp.Elev[sp.S2]))
	}
	PrintA("\n")
}

// selectDDpairsWithLowestElevSatRemoved removes the satellite with lowest elevation
// from double-difference pairs for retry strategy
func selectDDpairsWithLowestElevSatRemoved(ddPairs []SatPairFType, rovSpp *SppSol, maxElevToRemove float64) ([]SatPairFType, SatType) {
	minElev := 90.0
	var minElevSat SatType
	// Find satellite with lowest elevation angle
	for _, spf := range ddPairs {
		elev := ToDeg(rovSpp.Elev[spf.S2])
		if elev < minElev && elev < maxElevToRemove {
			minElev = elev
			minElevSat = spf.S2
		}
	}
	// If no satellite found below threshold, return original pairs
	if minElev == 90.0 {
		return ddPairs, minElevSat
	}
	// Remove pairs containing the lowest elevation satellite
	ddp := []SatPairFType{}
	for _, spf := range ddPairs {
		if spf.S2 != minElevSat {
			ddp = append(ddp, spf)
		}
	}
	return ddp, minElevSat
}

// foundIn checks if a satellite pair exists in the given list
// It matches by frequency and either satellite in the pair
func foundIn(ddp []SatPairFType, spf SatPairFType) bool {
	for _, spf2 := range ddp {
		// Check if satellites match (either direction) and frequency matches
		// This allows for flexible matching regardless of which satellite is S1 or S2
		if (spf.S2 == spf2.S1 || spf.S2 == spf2.S2) && spf.F == spf2.F {
			return true
		}
	}
	return false
}

// validateAndRetry validates the initial ambiguity solution and applies retry strategies
// if the ratio test fails to meet the threshold
func validateAndRetry(
	initialRatio float64,
	initialB, initialBf, initialRes []float64,
	initialDDPairs []SatPairFType,
	rovSpp *SppSol,
	floatSol *FloatSol,
	opt *AmbOpt,
) (float64, []float64, []float64, []float64, []SatPairFType, error) {

	currentRatio := initialRatio
	currentB := initialB
	currentBf := initialBf
	currentRes := initialRes
	currentDDPairs := initialDDPairs

	// Retry strategy 1: Remove low elevation satellites sequentially
	if currentRatio < opt.RatioThres {
		PrintD(2, "\tratio: %8.2f < %3.2f\n", currentRatio, opt.RatioThres)
		newRatio, newB, newBf, newRes, newDDPairs, err := retryWithLowElevSatRemoval(
			currentDDPairs, rovSpp, floatSol, opt)
		if err == nil && newRatio >= opt.RatioThres {
			currentRatio = newRatio
			currentB = newB
			currentBf = newBf
			currentRes = newRes
			currentDDPairs = newDDPairs
		}
	}

	// Retry strategy 2: Remove newly appeared satellites (commented out - low effectiveness)
	if currentRatio < opt.RatioThres && opt.PrevAmbSol != nil && false {
		PrintD(2, "\tratio: %8.2f < %3.2f\n", currentRatio, opt.RatioThres)
		newRatio, newB, newBf, newRes, newDDPairs, err := retryWithNewSatRemoval(
			currentDDPairs, rovSpp, floatSol, opt)
		if err == nil && newRatio >= opt.RatioThres {
			currentRatio = newRatio
			currentB = newB
			currentBf = newBf
			currentRes = newRes
			currentDDPairs = newDDPairs
		}
	}

	// Retry strategy 3: Remove GLONASS satellites if present (commented out)
	if currentRatio < opt.RatioThres && false {
		PrintD(2, "\tratio: %8.2f < %3.2f\n", currentRatio, opt.RatioThres)
		newRatio, newB, newBf, newRes, newDDPairs, err := retryWithGlonassRemoval(
			currentDDPairs, rovSpp, floatSol, opt)
		if err == nil && newRatio >= opt.RatioThres {
			currentRatio = newRatio
			currentB = newB
			currentBf = newBf
			currentRes = newRes
			currentDDPairs = newDDPairs
		}
	}

	if DBG_ >= 2 {
		PrintA("\t--- final ---\n")
		PrintA("\tdd pairs: (%d)", len(currentDDPairs))
		printDDpairs(currentDDPairs, rovSpp)
		PrintA("\tbiases: (%d)\n", len(currentB)/2)
		for i := range len(currentDDPairs) {
			PrintA("\t%10.3f ---> %10.1f\n", currentBf[i], currentB[i])
		}
	}

	// Check if number of DD pairs is too small compared to previous epoch
	// Even if fix is achieved, reject if too many pairs were removed (quality control)
	if currentRatio >= opt.RatioThres && opt.PrevAmbSol != nil {
		const NUM_DDPAIR_MIN_RATIO = 1 / 3.0 // Minimum 1/3 of previous epoch pairs
		minNumDDPairs := float64(len(opt.PrevAmbSol.DDPairs)) * NUM_DDPAIR_MIN_RATIO
		if float64(len(currentDDPairs)) < minNumDDPairs {
			return 0, currentB, currentBf, currentRes, currentDDPairs, fmt.Errorf("number of dd pairs is too small, ndd(prev)=%d, ndd(curr)=%d < %.1f", len(opt.PrevAmbSol.DDPairs), len(currentDDPairs), minNumDDPairs)
		}
	}

	return currentRatio, currentB, currentBf, currentRes, currentDDPairs, nil
}

// retryWithLowElevSatRemoval implements retry strategy by removing low elevation satellites
// It iteratively removes the satellite with lowest elevation and retries ambiguity resolution
func retryWithLowElevSatRemoval(
	ddPairs []SatPairFType,
	rovSpp *SppSol,
	floatSol *FloatSol,
	opt *AmbOpt,
) (float64, []float64, []float64, []float64, []SatPairFType, error) {

	currentPairs := make([]SatPairFType, len(ddPairs))
	copy(currentPairs, ddPairs)

	// Iteratively remove low elevation satellites and retry
	for i := 0; i < opt.MaxRetNum; i++ {
		newPairs, removedSat := selectDDpairsWithLowestElevSatRemoved(currentPairs, rovSpp, opt.MaxElevToRemove)

		// Check if any satellite was removed
		if len(newPairs) == len(currentPairs) {
			PrintD(2, "\tno sat pair to remove\n")
			break
		}

		// Ensure minimum number of DD pairs (at least 3 required for solution)
		if len(newPairs) < 3 {
			PrintD(2, "\tnot enough dd pairs, num of dd pairs=%d", len(newPairs))
			break
		}

		currentPairs = newPairs

		if DBG_ >= 2 {
			//PrintA("\t--- retry(1-%d) ---\n", i+1)
			PrintA("\t--- retry(%d/%d) ---\n", i+1, opt.MaxRetNum)
			PrintA("\tsat removed: %s(%.1f)\n", removedSat, ToDeg(rovSpp.Elev[removedSat]))
			PrintA("\tdd pairs: (%d)", len(currentPairs))
			printDDpairs(currentPairs, rovSpp)
		}

		// Solve ambiguities with reduced satellite set
		ratio, b, bf, res, err := solveAmbByLambda(currentPairs, floatSol)
		if err != nil {
			return 0, nil, nil, nil, currentPairs, fmt.Errorf("solveAmbByLambda() failed, err= %s", err.Error())
		}

		PrintD(2, "\tratio(%d/%d): %8.2f (%.3f/%.3f)\n", i+1, opt.MaxRetNum, ratio, res[1], res[0])

		// Check if ratio test passes (successful fix achieved)
		if ratio >= opt.RatioThres {
			return ratio, b, bf, res, currentPairs, nil
		}
	}

	if opt.MaxRetNum > 0 {
		return 0, nil, nil, nil, currentPairs, fmt.Errorf("retry failed")
	} else {
		return 0, nil, nil, nil, currentPairs, nil
	}
}

// retryWithNewSatRemoval implements retry strategy by removing newly appeared satellites
// This strategy filters out satellite pairs that were not present in the previous epoch
func retryWithNewSatRemoval(
	ddPairs []SatPairFType,
	rovSpp *SppSol,
	floatSol *FloatSol,
	opt *AmbOpt,
) (float64, []float64, []float64, []float64, []SatPairFType, error) {

	PrintD(2, "\t--- retry(2) ---\n")

	filteredPairs := filterPairsByPreviousEpoch(ddPairs, opt.PrevAmbSol.DDPairs)

	if len(filteredPairs) == len(ddPairs) {
		PrintD(2, "\tno sat pair to remove\n")
		return 0, nil, nil, nil, ddPairs, fmt.Errorf("no pairs to remove")
	}

	if len(filteredPairs) < 3 {
		PrintD(2, "\tnot enough dd pairs, num of dd pairs=%d", len(filteredPairs))
		return 0, nil, nil, nil, nil, fmt.Errorf("not enough pairs")
	}

	if DBG_ >= 2 {
		PrintA("\tdd pairs: (%d)", len(filteredPairs))
		printDDpairs(filteredPairs, rovSpp)
	}

	ratio, b, bf, res, err := solveAmbByLambda(filteredPairs, floatSol)
	if err != nil {
		return 0, nil, nil, nil, ddPairs, err
	}

	PrintD(2, "\tratio(2): %8.2f (%.3f/%.3f)\n", ratio, res[1], res[0])

	if ratio >= opt.RatioThres {
		return ratio, b, bf, res, filteredPairs, nil
	}

	return 0, nil, nil, nil, ddPairs, fmt.Errorf("retry(1) failed")
}

// retryWithGlonassRemoval implements retry strategy by removing GLONASS satellites
// This strategy filters out satellite pairs containing GLONASS satellites
func retryWithGlonassRemoval(
	ddPairs []SatPairFType,
	rovSpp *SppSol,
	floatSol *FloatSol,
	opt *AmbOpt,
) (float64, []float64, []float64, []float64, []SatPairFType, error) {

	PrintD(2, "\t--- retry(3) ---\n")

	filteredPairs := filterGlonassPairs(ddPairs)

	if len(filteredPairs) == len(ddPairs) {
		PrintD(2, "\tno sat pair to remove\n")
		return 0, nil, nil, nil, nil, fmt.Errorf("no pairs to remove")
	}

	if len(filteredPairs) < 3 {
		PrintD(2, "\tnot enough dd pairs, num of dd pairs=%d", len(filteredPairs))
		return 0, nil, nil, nil, nil, fmt.Errorf("not enough pairs")
	}

	if DBG_ >= 2 {
		PrintA("\tdd pairs: (%d)", len(filteredPairs))
		printDDpairs(filteredPairs, rovSpp)
	}

	ratio, b, bf, res, err := solveAmbByLambda(filteredPairs, floatSol)
	if err != nil {
		return 0, nil, nil, nil, nil, err
	}

	PrintD(2, "\tratio(3): %8.2f (%.3f/%.3f)\n", ratio, res[1], res[0])

	if ratio >= opt.RatioThres {
		return ratio, b, bf, res, filteredPairs, nil
	}

	return 0, nil, nil, nil, nil, fmt.Errorf("retry failed")
}

// filterPairsByPreviousEpoch filters satellite pairs to keep only those present in previous epoch
// This ensures continuity in ambiguity resolution across epochs
func filterPairsByPreviousEpoch(currentPairs, prevPairs []SatPairFType) []SatPairFType {
	filtered := []SatPairFType{}
	for _, spf := range currentPairs {
		if foundIn(prevPairs, spf) {
			filtered = append(filtered, spf)
		} else {
			PrintD(3, "\tsat pair removed: %s-%s(%d)\n", spf.S1, spf.S2, spf.F+1)
		}
	}
	return filtered
}

// filterGlonassPairs filters out satellite pairs containing GLONASS satellites
// GLONASS satellites are identified by system code 'R'
func filterGlonassPairs(ddPairs []SatPairFType) []SatPairFType {
	filtered := []SatPairFType{}
	for _, sp := range ddPairs {
		if sp.S1.Sys() != 'R' && sp.S2.Sys() != 'R' {
			filtered = append(filtered, sp)
		} else {
			PrintD(3, "\tsat pair removed: %s-%s(%d)\n", sp.S1, sp.S2, sp.F+1)
		}
	}
	return filtered
}

// recalcPos recalculates receiver position using fixed integer ambiguities
// It performs iterative least squares adjustment with carrier phase observations
func recalcPos(
	ddPairs []SatPairFType,
	B []float64,
	rovSpp, baseSpp *SppSol,
	basePos *PosXYZ,
	initialPos PosXYZ,
	opt *AmbOpt,
) (PosXYZ, error) {

	rovPos := NewPosXYZ(initialPos.X, initialPos.Y, initialPos.Z)

	if DBG_ >= 2 {
		llh := rovPos.ToLLH()
		llh.Lat = ToDeg(llh.Lat)
		llh.Lon = ToDeg(llh.Lon)
		PrintA("\tupos(init): LLH= %.8f %.8f %.3f, XYZ= %.3f %.3f %.3f\n",
			llh.Lat, llh.Lon, llh.Hei, rovPos.X, rovPos.Y, rovPos.Z)
	}

	// Iterative least squares adjustment
	for loop := 0; loop < MAX_LOOP; loop++ {
		// Build observation equations
		dy, H, R, err := buildObsEqn(ddPairs, B, rovSpp, baseSpp, basePos, rovPos, opt)
		if err != nil {
			return PosXYZ{}, err
		}

		if DBG_ >= 4 {
			PrintA("H=\n")
			PrintMat(H)
			PrintA("dy=\n")
			PrintMat(dy)
			PrintA("R=\n")
			PrintMat(R)
		}

		// Solve least squares equations
		dx, _, err := SolveLS(H, dy, R)
		if err != nil {
			return PosXYZ{}, fmt.Errorf("SolveLS() failed., err= %s", err.Error())
		}

		if DBG_ >= 4 {
			PrintA("dx=\n")
			PrintMat(dx)
		}

		// Update receiver position
		rovPos.X += dx.AtVec(0)
		rovPos.Y += dx.AtVec(1)
		rovPos.Z += dx.AtVec(2)

		if DBG_ >= 2 {
			llh := rovPos.ToLLH()
			llh.Lat = ToDeg(llh.Lat)
			llh.Lon = ToDeg(llh.Lon)
			msg := fmt.Sprintf("\tLOOP %d: LLH= %.9f %.9f %.4f, XYZ= %.3f %.3f %.3f\n",
				loop+1, llh.Lat, llh.Lon, llh.Hei, rovPos.X, rovPos.Y, rovPos.Z)
			PrintA(msg)
		}

		// Check convergence
		if isConverged(dx) {
			PrintD(2, "\tdiff: %f, %f, %f < %f\n",
				math.Abs(dx.AtVec(0)), math.Abs(dx.AtVec(1)), math.Abs(dx.AtVec(2)), CONVERGENCE_THRES)
			break
		}
	}

	return *rovPos, nil
}

// buildObsEqn builds observation equations for position calculation
// It constructs the design matrix, observation vector, and weight matrix
func buildObsEqn(
	ddPairs []SatPairFType,
	B []float64,
	rovSpp, baseSpp *SppSol,
	basePos *PosXYZ,
	rovPos *PosXYZ,
	opt *AmbOpt,
) (*mat.VecDense, *mat.Dense, *mat.Dense, error) {

	nv := len(ddPairs)
	dy := mat.NewVecDense(nv, nil)
	H := mat.NewDense(nv, 3, nil)
	R := mat.NewDense(nv, nv, nil)

	// Build equations using only carrier phase observations
	for i, sp := range ddPairs {

		residual, err := calcDDRes(sp, rovSpp, baseSpp, basePos, rovPos, B[i], opt)
		if err != nil {
			return nil, nil, nil, err
		}
		dy.SetVec(i, residual)

		designVector, err := calcDesignVec(sp, rovSpp, rovPos)
		if err != nil {
			return nil, nil, nil, err
		}
		H.Set(i, 0, designVector.X)
		H.Set(i, 1, designVector.Y)
		H.Set(i, 2, designVector.Z)

		// Set weights based on elevation angle (simple model)
		// Higher elevation satellites get higher weights
		elev := ToDeg(rovSpp.Elev[sp.S2])
		R.Set(i, i, math.Sqrt(90.0/elev))
	}

	return dy, H, R, nil
}

// calcDDRes calculates double-difference carrier phase residual
// It computes the difference between observed and computed carrier phase values
func calcDDRes(
	sp SatPairFType,
	rovSpp, baseSpp *SppSol,
	basePos *PosXYZ,
	rovPos *PosXYZ,
	B float64,
	opt *AmbOpt,
) (float64, error) {

	k, j, f := sp.S1, sp.S2, sp.F

	Xkr, Xjr := baseSpp.SatPos[k], baseSpp.SatPos[j]
	Xku, Xju := rovSpp.SatPos[k], rovSpp.SatPos[j]
	Rkr := EucDist(&Xkr, basePos)
	Rjr := EucDist(&Xjr, basePos)
	Rku := EucDist(&Xku, rovPos)
	Rju := EucDist(&Xju, rovPos)
	Ckr, Cjr := baseSpp.Cp[k][f], baseSpp.Cp[j][f]
	Cku, Cju := rovSpp.Cp[k][f], rovSpp.Cp[j][f]
	lmkr := C / baseSpp.Freq[k][f]
	lmjr := C / baseSpp.Freq[j][f]
	lmku := C / rovSpp.Freq[k][f]
	lmju := C / rovSpp.Freq[j][f]
	Ckr += C * baseSpp.SatClk[k] / lmkr
	Cjr += C * baseSpp.SatClk[j] / lmjr
	Cku += C * rovSpp.SatClk[k] / lmku
	Cju += C * rovSpp.SatClk[j] / lmju
	if k.Sys() == 'R' {
		Ckr -= baseSpp.Clk[0] / lmkr
		Cjr -= baseSpp.Clk[0] / lmjr
		Cku -= rovSpp.Clk[0] / lmku
		Cju -= rovSpp.Clk[0] / lmju
	}
	// Calculate double-difference in cycle units
	// (For GLONASS compatibility, equations after ambiguity fixing use cycle units, not meters)
	Cdd := (Cku - Ckr) - (Cju - Cjr)
	Rdd := (Rku/lmku - Rkr/lmkr) - (Rju/lmju - Rjr/lmjr)
	Bdd := B
	// Apply tropospheric correction if enabled
	if !opt.SkipTrop {
		// Calculate zenith hydrostatic delay for rover and base
		zhdu := TropModel(rovPos)
		zhdr := TropModel(basePos)
		// Calculate mapping function and tropospheric delay for each satellite
		Tkr := TropMapf(&baseSpp.Time, basePos, baseSpp.Elev[k]) * zhdr
		Tjr := TropMapf(&baseSpp.Time, basePos, baseSpp.Elev[j]) * zhdr
		Tku := TropMapf(&rovSpp.Time, rovPos, rovSpp.Elev[k]) * zhdu
		Tju := TropMapf(&rovSpp.Time, rovPos, rovSpp.Elev[j]) * zhdu
		// Calculate double-difference tropospheric correction
		Tdd := (Tku/lmku - Tkr/lmkr) - (Tju/lmju - Tjr/lmjr)
		Rdd += Tdd
	}
	return Cdd - Rdd - Bdd, nil
}

// calcDesignVec calculates the design vector for double-difference observations
// It computes the partial derivatives of the observation with respect to position
func calcDesignVec(
	sp SatPairFType,
	rovSpp *SppSol,
	rovPos *PosXYZ,
) (PosXYZ, error) {

	k, j, f := sp.S1, sp.S2, sp.F
	Xku, Xju := rovSpp.SatPos[k], rovSpp.SatPos[j]

	lmku := C / rovSpp.Freq[k][f]
	lmju := C / rovSpp.Freq[j][f]

	Ex := -(DistDx(&Xju, rovPos)/lmju - DistDx(&Xku, rovPos)/lmku)
	Ey := -(DistDy(&Xju, rovPos)/lmju - DistDy(&Xku, rovPos)/lmku)
	Ez := -(DistDz(&Xju, rovPos)/lmju - DistDz(&Xku, rovPos)/lmku)

	return PosXYZ{X: Ex, Y: Ey, Z: Ez}, nil
}

// isConverged checks if the position update has converged
// Returns true if all position components are below the convergence threshold
func isConverged(dx mat.Vector) bool {
	return math.Abs(dx.AtVec(0)) < CONVERGENCE_THRES &&
		math.Abs(dx.AtVec(1)) < CONVERGENCE_THRES &&
		math.Abs(dx.AtVec(2)) < CONVERGENCE_THRES
}

// setAmbSol creates an AmbSol structure with the given parameters
// This is a helper function to construct the final ambiguity solution
func setAmbSol(ratio float64, B []float64, Bf []float64, res []float64, ddp []SatPairFType, pos PosXYZ) *AmbSol {
	rslt := NewAmbSol()
	rslt.Pos = pos
	rslt.Ratio = ratio
	rslt.B = B
	rslt.Bf = Bf
	rslt.Res = res
	rslt.DDPairs = ddp
	return rslt
}
