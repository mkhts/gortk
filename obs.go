// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.27
//

package gortk

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
	"time"

	"golang.org/x/exp/slices"
)

// Type representing satellite name like "G10"
type SatType string

// Type representing satellite system like 'G'
type SysType byte

// Extract satellite system from satellite name
func (p *SatType) Sys() SysType {
	return SysType((*p)[0])
}

// Check validity of satellite system
func (p *SysType) IsValid() bool {
	return *p == 'G' || *p == 'J' || *p == 'E' || *p == 'R' || *p == 'C' || *p == 'S'
}

// Extract satellite number from satellite name
func (p *SatType) Num() int {
	i, err := strconv.Atoi(string((*p)[1:3]))
	if err != nil {
		return 0
	}
	return i
}

// Number of carrier frequencies
const NFREQ = 4

// Structure to store observation data for one satellite for one epoch
type ObsS struct {
	Pr   [NFREQ]float64  // Pseudorange
	Cp   [NFREQ]float64  // Carrier phase
	Dp   [NFREQ]float64  // Doppler frequency
	Sn   [NFREQ]float64  // Signal strength
	LLI  [NFREQ]byte     // LLI (Loss-of-Lock Indicator) (0: OK, 1: cycle slip, 2: possible half-cycle slip, 3: other problems)
	Freq [NFREQ]float64  // Carrier frequency
	Code [NFREQ]CodeType // Observation code (1C,2X,5I etc.)
}

// Constructor for the above structure
func NewObsS() *ObsS {
	return &ObsS{
		Pr:   [NFREQ]float64{},
		Cp:   [NFREQ]float64{},
		Dp:   [NFREQ]float64{},
		Sn:   [NFREQ]float64{},
		LLI:  [NFREQ]uint8{},
		Freq: [NFREQ]float64{},
		Code: [NFREQ]CodeType{},
	}
}

// Structure to store observation data for all satellites in one epoch
type ObsE struct {
	Time GTime             // Epoch time
	DatS map[SatType]*ObsS // Observation data for each satellite
}

// Return map keys as slice
func (p *ObsE) Sats() []SatType {
	s := make([]SatType, 0)
	for k := range p.DatS {
		s = append(s, k)
	}
	return s
}

// Structure to store observation data for all epochs
type Obs struct {
	DatE  []*ObsE                // Satellite data for each time (sorted by time in ascending order)
	Codes map[SysType][]CodeType // List of observation codes contained in file
}

// Display observation data overview
func (p *Obs) String() string {
	if len(p.DatE) == 0 {
		return "NO DATA"
	}
	// Satellite list
	sl := map[SysType][]SatType{}
	for _, obse := range p.DatE {
		for sat, _ := range obse.DatS {
			if a, ok := sl[sat.Sys()]; ok {
				if !slices.Contains(a, sat) {
					sl[sat.Sys()] = append(sl[sat.Sys()], sat)
				}
			} else {
				sl[sat.Sys()] = []SatType{sat}
			}
		}
	}
	// Observation codes (first epoch)
	// cl := map[SysType][NFREQ]CodeType{}
	// for sat, obss := range p.DatE[0].DatS {
	// 	if _, ok := cl[sat.Sys()]; !ok {
	// 		cl[sat.Sys()] = [NFREQ]CodeType{}
	// 	}
	// 	t := cl[sat.Sys()]
	// 	for i := range NFREQ {
	// 		if obss.Freq[i] > 0 {
	// 			t[i] = obss.Code[i]
	// 		}
	// 	}
	// 	cl[sat.Sys()] = t
	// }
	// String conversion
	var sb, sb2 strings.Builder
	for _, sys := range []SysType{'G', 'J', 'E', 'R', 'C', 'S'} {
		if a, ok := sl[sys]; ok {
			if len(a) > 0 {
				sb.WriteString(fmt.Sprintf("\t%c (%2d):", sys, len(a)))
				sort.Slice(a, func(i, j int) bool {
					return a[i] < a[j]
				})
				for _, b := range a {
					sb.WriteString(fmt.Sprintf(" %s", b[1:]))
				}
				sb.WriteString("\n")
			}
		}
		if a, ok := p.Codes[sys]; ok {
			if len(p.Codes[sys]) > 0 {
				sb2.WriteString(fmt.Sprintf("\t%c (%2d):", sys, len(p.Codes[sys])))
				for _, b := range a {
					sb2.WriteString(fmt.Sprintf(" %s", b))
				}
				sb2.WriteString("\n")
			}
		}
		// if a, ok := cl[sys]; ok {
		// 	sb3.WriteString(fmt.Sprintf("\t%c:", sys))
		// 	for i := range NFREQ {
		// 		var b CodeType
		// 		if a[i] != b {
		// 			sb3.WriteString(fmt.Sprintf(" %s", a[i]))
		// 		} else {
		// 			sb3.WriteString("   ")
		// 		}
		// 	}
		// 	sb3.WriteString("\n")
		// }
	}
	a := `
datetime:
	%s - %s (%d)

sats:
%s
codes:
%s`
	//return fmt.Sprintf(a, p.DatE[0].Time.ToTime().UTC().Format("2006/01/02 15:04:05.000"), p.DatE[len(p.DatE)-1].Time.ToTime().UTC().Format("2006/01/02 15:04:05.000"), len(p.DatE), sb.String(), sb2.String(), sb3.String())
	return fmt.Sprintf(a, p.DatE[0].Time.ToTime().UTC().Format("2006/01/02 15:04:05.000"), p.DatE[len(p.DatE)-1].Time.ToTime().UTC().Format("2006/01/02 15:04:05.000"), len(p.DatE), sb.String(), sb2.String())
}

// Return data for the epoch closest in time to the specified time
func (p *Obs) GetNearest(t GTime) (obse *ObsE, err error) {
	if len(p.DatE) > 0 {
		const MAX_AGE = float64(time.Second * 30)
		m := MAX_AGE + 1
		for _, s := range p.DatE {
			if false { // Set to false if allowing base station data time to be in the future relative to rover data time
				if t.Less(s.Time, true) {
					continue
				}
			}
			d := math.Abs(t.ToTime().Sub(s.Time.ToTime()).Seconds())
			if d <= m {
				obse = s
				m = d
			}
		}
		if m > MAX_AGE {
			return nil, fmt.Errorf("no nearest data is found within %d seconds. t=%s, m=%f", int(MAX_AGE), t.ToTime(), m)
		} else {
			return obse, nil
		}
	} else {
		return nil, fmt.Errorf("the container is empty")
	}
}
