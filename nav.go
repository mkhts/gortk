// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.5.30
//

package gortk

import (
	"fmt"
	"strings"
	"time"
)

// Structure to store ephemeris (navigation data for one satellite, one issue)
type Ephe struct {

	// Common for G,J,E,C,R
	Sat  SatType
	Toc  GTime // Reference time for satellite clock error correction
	Toe  GTime // Reference time for satellite orbit calculation
	Tot  GTime // Transmission time
	Iode int

	// for GPS, QZSS, GALILEO, BEIDOU
	Af0    float64
	Af1    float64
	Af2    float64
	Crs    float64
	DeltaN float64
	M0     float64
	Cuc    float64
	Ecc    float64
	Cus    float64
	SqrtA  float64
	Cic    float64
	Omega0 float64
	Cis    float64
	I0     float64
	Crc    float64
	Omega  float64
	OmegaD float64
	Idot   float64
	Code   int
	Week   int
	Flag   int
	Sva    int
	Svh    int
	Tgd    float64 // GPS, QZS, GAL(E5a/E1), BDS(B1/B3)
	Tgd2   float64 // GAL(E5b/E1), BDS(B2/B3) // Galileo and Beidou have two group delay parameters
	Iodc   int     // GPS, QZS, BDS
	Fit    float64 // GPS, QZS

	// for GLONASS
	TauN   float64
	GammaN float64
	PosX   float64
	VecX   float64
	AccX   float64
	PosY   float64
	VecY   float64
	AccY   float64
	FreqN  int
	PosZ   float64
	VecZ   float64
	AccZ   float64
	Age    int

	// for SBAS
	Gf0  float64
	Gf1  float64
	Iodn int
}

func (e *Ephe) String() string {
	str := ""
	str += fmt.Sprintf("### Nav. for %s (%c, %d)\n", e.Sat, e.Sat.Sys(), e.Sat.Num())
	str += fmt.Sprintf("    Toc: %v (%v)\n", e.Toc.ToTime().UTC(), e.Toc)
	str += fmt.Sprintf("    Toe: %v (%v)\n", e.Toe.ToTime().UTC(), e.Toe)
	str += fmt.Sprintf("    Tot: %v (%v)\n", e.Tot.ToTime().UTC(), e.Tot)
	str += fmt.Sprintf("   Iode: %v\n", e.Iode)
	str += fmt.Sprintf("    Af0: %v\n", e.Af0)
	str += fmt.Sprintf("    Af1: %v\n", e.Af1)
	str += fmt.Sprintf("    Af2: %v\n", e.Af2)
	str += fmt.Sprintf("    Crs: %v\n", e.Crs)
	str += fmt.Sprintf(" DeltaN: %v\n", e.DeltaN)
	str += fmt.Sprintf("     M0: %v\n", e.M0)
	str += fmt.Sprintf("    Cuc: %v\n", e.Cuc)
	str += fmt.Sprintf("    Ecc: %v\n", e.Ecc)
	str += fmt.Sprintf("    Cus: %v\n", e.Cus)
	str += fmt.Sprintf("  SqrtA: %v\n", e.SqrtA)
	str += fmt.Sprintf("    Cic: %v\n", e.Cic)
	str += fmt.Sprintf(" Omega0: %v\n", e.Omega0)
	str += fmt.Sprintf("    Cis: %v\n", e.Cis)
	str += fmt.Sprintf("     I0: %v\n", e.I0)
	str += fmt.Sprintf("    Crc: %v\n", e.Crc)
	str += fmt.Sprintf("  Omega: %v\n", e.Omega)
	str += fmt.Sprintf(" OmegaD: %v\n", e.OmegaD)
	str += fmt.Sprintf("   Idot: %v\n", e.Idot)
	str += fmt.Sprintf("   Code: %v\n", e.Code)
	str += fmt.Sprintf("   Week: %v\n", e.Week)
	str += fmt.Sprintf("   Flag: %v\n", e.Flag)
	str += fmt.Sprintf("    Sva: %v\n", e.Sva)
	str += fmt.Sprintf("    Svh: %v\n", e.Svh)
	str += fmt.Sprintf("    Tgd: %v\n", e.Tgd)
	str += fmt.Sprintf("   Tgd2: %v\n", e.Tgd2)
	str += fmt.Sprintf("   Iodc: %v\n", e.Iodc)
	str += fmt.Sprintf("    Fit: %v\n", e.Fit)
	str += fmt.Sprintf("   TauN: %v\n", e.TauN)
	str += fmt.Sprintf(" GammaN: %v\n", e.GammaN)
	str += fmt.Sprintf("   PosX: %v\n", e.PosX)
	str += fmt.Sprintf("   VecX: %v\n", e.VecX)
	str += fmt.Sprintf("   AccX: %v\n", e.AccX)
	str += fmt.Sprintf("   PosY: %v\n", e.PosY)
	str += fmt.Sprintf("   VecY: %v\n", e.VecY)
	str += fmt.Sprintf("   AccY: %v\n", e.AccY)
	str += fmt.Sprintf("  FreqN: %v\n", e.FreqN)
	str += fmt.Sprintf("   PosZ: %v\n", e.PosZ)
	str += fmt.Sprintf("   VecZ: %v\n", e.VecZ)
	str += fmt.Sprintf("   AccZ: %v\n", e.AccZ)
	str += fmt.Sprintf("    Age: %v\n", e.Age)
	str += fmt.Sprintf("    Gf0: %v\n", e.Gf0)
	str += fmt.Sprintf("    Gf1: %v\n", e.Gf1)
	str += fmt.Sprintf("   Iodn: %v\n", e.Iodn)
	return str
}

// Structure to store navigation data for each satellite at each time
// - Map with satellite name as Key and slice sorted by transmission time (Tot) in ascending order as Value
type Nav map[SatType][]*Ephe

// Determine if ephemeris is valid for the specified time
func (navs *Ephe) isValid(dt time.Time) bool {
	// Time difference between specified time and TOC time must be within 2 hours
	td := dt.Sub(navs.Toc.ToTime())
	//fmt.Println(dt.UTC(), navs.Toc.ToTime().UTC(), navs.Tot.ToTime().UTC(), td, (td >= time.Duration(time.Second*-7201)), (td <= time.Duration(time.Second*7201)), (navs.Tot.ToTime().Sub(dt) < 1))
	if td >= time.Duration(time.Second*-7201) && td <= time.Duration(time.Second*7201) { // (1 second buffer is set to handle rounding errors)
		// And transmission time must be in the past relative to specified time
		//if navs.Tot.ToTime().Sub(dt) < 1 {
		if navs.Tot.Before(dt, true) {
			return true
		}
	}
	return false
}

// Function to select ephemeris by specifying satellite and time: for GPS, QZSS, GLONASS, BEIDOU
func (nav *Nav) GetEphe(sat SatType, gt GTime) (*Ephe, error) {
	if false {
		// 1) Select ephemeris where the difference between specified time and ToC is within specified time and time is after transmission time (ToT) (Practical programming method)
		dt := gt.ToTime()
		if navs, ok := (*nav)[sat]; ok {
			// Search from newest transmission time and return the first valid ephemeris found
			for i := len(navs) - 1; i >= 0; i-- {
				if navs[i].isValid(dt) {
					return (*nav)[sat][i], nil
				}
			}
		} else {
			return nil, fmt.Errorf("can't find %s", sat)
		}
		return nil, fmt.Errorf("can't find a valid ephemeris for %s", sat)
	} else {
		// 2) Select ephemeris closest to ToE including future, where the difference between specified time and ToE is within specified time (RTKLIB, RTK core method)
		var diffMax float64
		switch sat.Sys() {
		case 'E':
			diffMax = 14400 // Following RTKLIB's MAXDTOE_GAL
		case 'C':
			diffMax = 21601 // Following RTKLIB's MAXDTOE_CMP
		default:
			diffMax = 7201
		}
		dt := gt.ToTime()
		j := -1
		if navs, ok := (*nav)[sat]; ok {
			for i, eph := range navs {
				// For GALILEO, future ToE is not allowed (RTKLIB does this)
				if sat.Sys() == 'E' && eph.Toe.ToTime().Sub(dt).Seconds() >= 0 {
					//if sat.Sys() == 'E' && eph.Toe.After(dt, true) {
					continue
				}
				diff := eph.Toe.ToTime().Sub(dt).Abs().Seconds()
				if diff < diffMax {
					diffMax = diff
					j = i
				}
			}
			if j >= 0 {
				return (*nav)[sat][j], nil
			} else {
				return nil, fmt.Errorf("can't find a valid ephemeris for %s", sat)
			}
		}
		return nil, fmt.Errorf("can't find %s", sat)
	}
}

// Display navigation data overview
func (p *Nav) String() string {
	keys := []SatType{}
	for k := range *p {
		keys = append(keys, k)
	}
	keys = Sorted(keys)
	var sb strings.Builder
	sb.WriteString("toc:\n")
	for _, sat := range keys {
		sb.WriteString(fmt.Sprintf("\t%s: ", sat))
		if len((*p)[sat]) > 0 {
			st := (*p)[sat][0].Toc
			et := (*p)[sat][len((*p)[sat])-1].Toc
			sb.WriteString(fmt.Sprintf("%s - %s (%d)\n",
				st.ToTime().UTC().Format("2006/01/02 15:04:05.000"), et.ToTime().UTC().Format("2006/01/02 15:04:05.000"), len((*p)[sat])))
		} else {
			sb.WriteString("\n")
		}
	}
	return sb.String()
}
