// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.20
//

package gortk

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"
)

// RINEX 3.04 specification
// https://files.igs.org/pub/data/format/rinex304.pdf
//

// Type representing observation codes like C1C (3 or 2 characters)
type CodeType string

// Returns observation type (C,L,D,S)
func (p *CodeType) T() byte {
	return (*p)[0]
}

// Returns frequency band and attributes of observation (1C,2P,5I etc.)
func (p *CodeType) NA() CodeType {
	return CodeType(*p)[1:]
}

// Priority and corresponding frequency settings for observation codes used in calculation
var CODE_ASSIGNS = map[SysType]map[CodeType]struct {
	priority int
	freqIdx  int
	freq     float64
}{
	'G': {
		"1C": {0, 0, 1.57542e9}, // L1
		"1P": {1, 0, 1.57542e9},
		"1Y": {2, 0, 1.57542e9},
		"1W": {3, 0, 1.57542e9},
		"1M": {4, 0, 1.57542e9},
		"1N": {5, 0, 1.57542e9},
		"1S": {6, 0, 1.57542e9},
		"1L": {7, 0, 1.57542e9},
		"1X": {8, 0, 1.57542e9},
		"2C": {0, 1, 1.22760e9}, // L2
		"2P": {1, 1, 1.22760e9},
		"2Y": {2, 1, 1.22760e9},
		"2W": {3, 1, 1.22760e9},
		"2M": {4, 1, 1.22760e9},
		"2N": {5, 1, 1.22760e9},
		"2D": {6, 1, 1.22760e9},
		"2L": {7, 1, 1.22760e9},
		"2S": {8, 1, 1.22760e9},
		"2X": {9, 1, 1.22760e9},
		"5I": {0, 2, 1.17645e9}, // L5
		"5Q": {1, 2, 1.17645e9},
		"5X": {2, 2, 1.17645e9},
	},
	'J': {
		"1C": {0, 0, 1.57542e9}, // L1
		"1L": {1, 0, 1.57542e9},
		"1S": {2, 0, 1.57542e9},
		"1X": {3, 0, 1.57542e9},
		"1Z": {4, 0, 1.57542e9},
		"2L": {5, 1, 1.22760e9}, // L2
		"2S": {6, 1, 1.22760e9},
		"2X": {7, 1, 1.22760e9},
		"5I": {8, 2, 1.17645e9}, // L5
		"5Q": {9, 2, 1.17645e9},
		"5X": {10, 2, 1.17645e9},
		"5D": {11, 2, 1.17645e9},
		"5P": {12, 2, 1.17645e9},
		"5Z": {13, 2, 1.17645e9},
	},
	'E': {
		"1C": {0, 0, 1.57542e9}, // E1
		"1A": {1, 0, 1.57542e9},
		"1B": {2, 0, 1.57542e9},
		"1X": {3, 0, 1.57542e9},
		"1Z": {4, 0, 1.57542e9},
		"7X": {5, 1, 1.20714e9}, // E5b
		"7I": {6, 1, 1.20714e9},
		"7Q": {7, 1, 1.20714e9},
		"5X": {8, 2, 1.17645e9}, // E5a
		"5I": {9, 2, 1.17645e9},
		"5Q": {10, 2, 1.17645e9},
		"8I": {11, 3, 1.191795e9}, // E5a+E5b
		"8Q": {12, 3, 1.191795e9},
		"8X": {13, 3, 1.191795e9},
		"6A": {14, 4, 1.27875e9}, // E6
		"6B": {15, 4, 1.27875e9},
		"6C": {16, 4, 1.27875e9},
		"6X": {17, 4, 1.27875e9},
		"6Z": {18, 4, 1.27875e9},
	},
	'R': {
		"1C": {0, 0, 1.60200e9}, // G1 FDMA
		"1P": {1, 0, 1.60200e9},
		"4A": {2, 0, 1.600995e9}, // G1a
		"4B": {3, 0, 1.600995e9},
		"4X": {4, 0, 1.600995e9},
		"2C": {5, 1, 1.24600e9}, // G2 FDMA
		"2P": {6, 1, 1.24600e9},
		"6A": {7, 1, 1.248060e9}, // G2a
		"6B": {8, 1, 1.248060e9},
		"6X": {9, 1, 1.248060e9},
		"3I": {10, 2, 1.202025e9}, // G3
		"3Q": {11, 2, 1.202025e9},
		"3X": {12, 2, 1.202025e9},
	},
	'C': {
		"2I": {0, 0, 1.561098e9}, // B1-2
		"2Q": {1, 0, 1.561098e9},
		"2X": {2, 0, 1.561098e9},
		"1D": {3, 0, 1.57542e9}, // B1
		"1P": {4, 0, 1.57542e9},
		"1X": {5, 0, 1.57542e9},
		"1A": {6, 0, 1.57542e9},
		"1N": {7, 0, 1.57542e9},
		"7I": {8, 1, 1.20714e9}, // B2b
		"7Q": {9, 1, 1.20714e9},
		"7X": {10, 1, 1.20714e9},
		"7D": {11, 1, 1.20714e9},
		"7P": {12, 1, 1.20714e9},
		"7Z": {13, 1, 1.20714e9},
		"6I": {14, 2, 1.26852e9}, // B3
		"6Q": {15, 2, 1.26852e9},
		"6X": {16, 2, 1.26852e9},
		"6A": {17, 2, 1.26852e9},
		"5D": {18, 3, 1.17645e9}, // B2a
		"5P": {19, 3, 1.17645e9},
		"5X": {20, 3, 1.17645e9},
		"8D": {21, 4, 1.191795e9}, // B2a+B2b
		"8P": {22, 4, 1.191795e9},
		"8X": {23, 4, 1.191795e9},
	},
	'S': {
		"1C": {0, 0, 1.57542e9}, // L1
		"5I": {1, 1, 1.17645e9}, // L5
		"5Q": {2, 1, 1.17645e9},
		"5X": {3, 1, 1.17645e9},
	},
}

// Extract HEADER LABEL string from observation data file header line
func getHeaderLabel(l string) string {
	if len(l) < 60 {
		return ""
	}
	return strings.TrimSpace(l[60:])
}

// Fix Beidou B1 observation codes in RINEX 3.02
func fixRnx302BeidouCode(la []string) []string {
	la2 := []string{}
	for _, a := range la {
		if a[1:3] == "1I" || a[1:3] == "1Q" || a[1:3] == "1X" {
			// In RINEX 3.04, B1(1561.098 MHz) observation codes {C|L|D|S}1{I|Q|X} have been changed to {C|L|D|S}2{I|Q|X}. Match 3.04.
			la2 = append(la2, a[:1]+"2"+a[2:3])
		} else {
			la2 = append(la2, a)
		}
	}
	return la2
}

// Read date and time from observation data file epoch line
func getObsTime(l string) (gt GTime, ns int, err error) {
	la := strings.Fields(l)
	if len(la) > 8 {
		year, err := strconv.ParseInt(la[1], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		month, err := strconv.ParseInt(la[2], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		day, err := strconv.ParseInt(la[3], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		hour, err := strconv.ParseInt(la[4], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		minute, err := strconv.ParseInt(la[5], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		la2 := strings.Split(la[6], ".")
		var sec, nsec int64
		if len(la2) == 2 {
			sec, err = strconv.ParseInt(la2[0], 10, 0)
			if err != nil {
				return gt, 0, err
			}
			nsec, err = strconv.ParseInt(la2[1], 10, 0)
			if err != nil {
				return gt, 0, err
			}
			nsec *= 100
		} else {
			return gt, 0, fmt.Errorf("invalid format in epoch line: %s (%s)", l, la[6])
		}
		ns, err := strconv.ParseInt(la[8], 10, 0)
		if err != nil {
			return gt, 0, err
		}
		return *NewGTime(time.Date(int(year), time.Month(month), int(day), int(hour), int(minute), int(sec), int(nsec), time.UTC)), int(ns), nil
	} else {
		return gt, 0, fmt.Errorf("not enought fields in epoch line: %s (%d)", l, len(la))
	}
}

// Set values according to observation code
func setValObsS(val float64, lli byte, sys SysType, code CodeType, out *ObsS) {
	if a, ok := CODE_ASSIGNS[sys][code.NA()]; ok {
		if a.freqIdx < NFREQ {
			if out.Freq[a.freqIdx] != 0 { // If value is already written
				b := CODE_ASSIGNS[sys][out.Code[a.freqIdx]] // Observation code of already written value
				if a.priority > b.priority {                // If priority is lower than that, don't write and exit
					return
				}
			}
			// Write value
			out.Freq[a.freqIdx] = a.freq
			if lli > 0 {
				out.LLI[a.freqIdx] = lli
			}
			out.Code[a.freqIdx] = code.NA()
			switch code.T() {
			case 'C':
				out.Pr[a.freqIdx] = val
			case 'L':
				out.Cp[a.freqIdx] = val
			case 'D':
				out.Dp[a.freqIdx] = val
			case 'S':
				out.Sn[a.freqIdx] = val
			}
		}
	}
}

// Read each observation value from observation data line
func getObsData(l string, oc map[SysType][]CodeType) (SatType, *ObsS, error) {
	if len(l) > 0 {
		sys := l[:1]
		num, _ := strconv.Atoi(strings.TrimSpace(l[1:3]))
		sat := SatType(fmt.Sprintf("%s%02d", sys, num))
		b := sat.Sys()
		if !b.IsValid() {
			return sat, nil, fmt.Errorf("unknown satellite system, '%c'", b)
		}
		n := len(oc[sat.Sys()])
		if len(l) < n*16+3 { // Fill in blanks if omitted to end of line
			l = l + strings.Repeat(" ", n*16+3-len(l))
		}
		// Set values according to observation code
		obsS := NewObsS()
		for i, code := range oc[sat.Sys()] {
			j := 3 + 16*i
			v, err := strconv.ParseFloat(strings.TrimSpace(l[j:j+14]), 64)
			if err != nil {
				continue
			}
			lli, err := strconv.ParseUint(strings.TrimSpace(l[j+14:j+15]), 10, 8)
			if err != nil {
				lli = 0
			}
			setValObsS(v, byte(lli), sat.Sys(), code, obsS)
		}
		return sat, obsS, nil
	} else {
		return "", nil, fmt.Errorf("can't read data. the given line is empty")
	}
}

// Read observation data
func ReadObs(r io.Reader) (*Obs, error) {

	// Flag indicating header reading is complete
	headerDone := false

	// RINEX version
	var ver string

	// List of observation codes in header
	oc := map[SysType][]CodeType{
		'G': make([]CodeType, 0),
		'J': make([]CodeType, 0),
		'E': make([]CodeType, 0),
		'R': make([]CodeType, 0),
		'C': make([]CodeType, 0),
		'S': make([]CodeType, 0),
	}

	// Variable to hold satellite data for one epoch during reading
	obsE := &ObsE{}

	// Temporarily store all epoch data in map to eliminate duplicates
	me := map[GTime]*ObsE{}

	// Reader to read line by line with newline as delimiter
	s := bufio.NewScanner(r)

	// Read line by line
	for s.Scan() {

		// Read line
		line := s.Text()

		// Process header lines
		if !headerDone {

			// Check version and file type
			if getHeaderLabel(line) == "RINEX VERSION / TYPE" {
				ver = line[5:9]
				if ver != "3.02" && ver != "3.04" {
					return nil, fmt.Errorf("unsupported RINEX version. RINEX version must be ether 3.02 or 3.04 (ver=%s)", ver)
				}
				typ := line[20:21]
				if typ != "O" {
					return nil, fmt.Errorf("not an observation data file (typ=%s)", typ)
				}
			}

			// Read Observation Code
			if getHeaderLabel(line) == "SYS / # / OBS TYPES" {
				sys := SysType(line[:1][0])
				if _, ok := oc[sys]; ok {
					la := strings.Fields(line[6:60])
					if ver == "3.02" && sys == 'C' {
						la = fixRnx302BeidouCode(la)
					}
					for _, code := range la {
						oc[sys] = append(oc[sys], CodeType(code))
					}
					nc, err := strconv.ParseInt(strings.TrimSpace(line[1:6]), 10, 0)
					if err == nil && nc >= 14 && s.Scan() { // When Code spans 2 lines
						line = s.Text()
						la = strings.Fields(line[6:60])
						if ver == "3.02" && sys == 'C' {
							la = fixRnx302BeidouCode(la)
						}
						for _, code := range la {
							oc[sys] = append(oc[sys], CodeType(code))
						}
					}
				}
			}

			// Check end of header lines
			if getHeaderLabel(line) == "END OF HEADER" {
				headerDone = true
			}

		} else { // Process observation data lines

			switch line[:1] {
			case ">":
				// Store satellite data read so far if any
				if len(obsE.DatS) > 0 {
					me[obsE.Time] = obsE
				}
				// Read epoch date and time
				t, ns, err := getObsTime(line)
				if err != nil {
					PrintD(2, "getObsTime() failed. err=%s", err.Error())
					continue
				}
				obsE = &ObsE{
					Time: t,
					DatS: make(map[SatType]*ObsS, ns),
				}
			default: // Read satellite data
				sat, obsS, err := getObsData(line, oc)
				if err != nil {
					PrintD(2, "getObsData() failed. err=%s", err.Error())
					continue
				}
				if obsE.DatS != nil {
					obsE.DatS[sat] = obsS
				}
			}
		}
	}

	// Check if reading completed without error
	if err := s.Err(); err != nil {
		return nil, err
	}

	// Set last epoch data
	if len(obsE.DatS) > 0 {
		me[obsE.Time] = obsE
	}

	// Sort data by date and time
	t := make([]GTime, 0, len(me))
	for k := range me {
		t = append(t, k)
	}
	sort.Slice(t, func(i, j int) bool {
		return t[i].Before(t[j].ToTime(), false)
	})

	// Set return value
	obs := &Obs{
		DatE:  make([]*ObsE, 0, len(me)),
		Codes: oc,
	}
	for _, a := range t {
		obs.DatE = append(obs.DatE, me[a])
	}

	return obs, nil
}

// Read satellite name and ToC from navigation data epoch line
func getNavTime(l string) (gt GTime, sat SatType, err error) {
	m := regexp.MustCompile(`^([GJERCS])([0-9 ][0-9]) (\d{4}) ([ \d]{2}) ([ \d]{2}) ([ \d]{2}) ([ \d]{2}) ([ \d]{2})`)
	ms := m.FindAllStringSubmatch(l, -1)
	if ms == nil {
		return gt, sat, fmt.Errorf("regexp match failed. l=%s", l)
	}
	sys := SysType(ms[0][1][0])
	num, err := strconv.ParseInt(strings.TrimSpace(ms[0][2]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	sat = SatType(fmt.Sprintf("%c%02d", sys, num))
	year, err := strconv.ParseInt(strings.TrimSpace(ms[0][3]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	month, err := strconv.ParseInt(strings.TrimSpace(ms[0][4]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	day, err := strconv.ParseInt(strings.TrimSpace(ms[0][5]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	hour, err := strconv.ParseInt(strings.TrimSpace(ms[0][6]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	min, err := strconv.ParseInt(strings.TrimSpace(ms[0][7]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	sec, err := strconv.ParseInt(strings.TrimSpace(ms[0][8]), 10, 0)
	if err != nil {
		return gt, sat, err
	}
	if sys == 'C' {
		sec += 14
	}
	gt = *NewGTime(time.Date(int(year), time.Month(month), int(day), int(hour), int(min), int(sec), 0, time.UTC))
	return
}

// Read navigation data
func ReadNav(r io.Reader) (*Nav, error) {

	// Flag indicating header reading is complete
	headerDone := false

	// RINEX version
	var ver string

	// Variable to store navigation data
	nav := Nav{}

	// Variable to hold ephemeris information during reading
	var eph *Ephe

	// Variable to hold which satellite system the data being read belongs to
	var sys SysType

	// Current line number being read, counted from satellite name and ToC line
	var lineCount int = 0

	// Reader to read line by line with newline as delimiter
	s := bufio.NewScanner(r)

	// Read line by line
	for s.Scan() {

		// Read line
		line := s.Text()

		// Process header lines
		if !headerDone {
			// Check version and file type
			if getHeaderLabel(line) == "RINEX VERSION / TYPE" {
				ver = line[5:9]
				if ver != "3.02" && ver != "3.04" {
					return nil, fmt.Errorf("unsupported RINEX version. RINEX version must be ether 3.02 or 3.04 (ver=%s)", ver)
				}
				typ := line[20:21]
				if typ != "N" {
					return nil, fmt.Errorf("not a navigation message file (typ=%s)", typ)
				}
			}

			// TODO: Read ionospheric correction information
			if getHeaderLabel(line) == "IONOSPHERIC CORR" {
			}

			// TODO: Read time correction information
			if getHeaderLabel(line) == "TIME SYSTEM CORR" {
			}

			// Check end of header lines
			if getHeaderLabel(line) == "END OF HEADER" {
				headerDone = true
			}

		} else { // Process navigation message lines

			m := regexp.MustCompile(`[- +\d]{2}\.\d{12}[DE][-+]\d{2}`)
			//fmt.Println(m.MatchString(line), line)
			if !m.MatchString(line) {
				continue
			}

			var err error
			sys_ := SysType(line[:1][0])
			switch sys_ {
			case 'G', 'J', 'E', 'C', 'R', 'S':
				if len(line) >= 80 {
					sys = sys_
					eph = &Ephe{}
					eph.Toc, eph.Sat, err = getNavTime(line)
					if err != nil {
						return nil, fmt.Errorf("failed to read time of clock in navigation message. err=%s", err.Error())
					}
					switch sys {
					case 'G', 'J', 'E', 'C':
						eph.Af0 = parseFloat(line[23:42])
						eph.Af1 = parseFloat(line[42:61])
						eph.Af2 = parseFloat(line[61:80])
					case 'R':
						eph.TauN = -parseFloat(line[23:42])
						eph.GammaN = parseFloat(line[42:61])
						toc15 := GTime{Week: eph.Toc.Week, Sec: math.Floor((eph.Toc.Sec+450)/900) * 900}
						dow := math.Floor(eph.Toc.Sec / 86400.0)
						tod := math.Mod(parseFloat(line[61:80]), 86400)
						eph.Tot = GTime{Week: eph.Toc.Week, Sec: tod + dow*86400}
						if eph.Tot.Sec-toc15.Sec < -43200 {
							eph.Tot.Sec += 86400
						} else if eph.Tot.Sec-toc15.Sec > 43200 {
							eph.Tot.Sec -= 86400
						}
						//eph.Toe = GTime{Week: eph.Toc.Week, Sec: eph.Toc.Sec + LS} // Glonass is in UTC time, so convert to GPS time
						eph.Toe = GTime{Week: toc15.Week, Sec: toc15.Sec + LS} // RTKLIB uses the time rounded to every 15 minutes.
						eph.Tot = GTime{Week: eph.Tot.Week, Sec: eph.Tot.Sec + LS}
						eph.Iode = int(math.Mod(eph.Toc.Sec+10800.0, 86400.0)/900.0 + 0.5)
					case 'S':
						eph.Gf0 = parseFloat(line[23:42])
						eph.Gf1 = parseFloat(line[42:61])
						eph.Tot = GTime{Week: eph.Toc.Week, Sec: parseFloat(line[61:80])}
						if eph.Tot.Sec-eph.Toc.Sec < -302400 {
							eph.Tot.Sec += 604800
						} else if eph.Tot.Sec-eph.Toc.Sec > 302400 {
							eph.Tot.Sec -= 604800
						}
						// The specification of SBAS's ToE is unclear, so it is processed the same as R.
						toc15 := GTime{Week: eph.Toc.Week, Sec: math.Floor((eph.Toc.Sec+450)/900) * 900}
						eph.Toe = GTime{Week: toc15.Week, Sec: toc15.Sec + LS} // RTKLIB uses the time rounded to every 15 minutes.
					}
					lineCount = 0
				}
			case ' ':
				var v0, v1, v2, v3 float64
				if len(line) < 80 {
					line = line + strings.Repeat(" ", 80-len(line))
				}
				v0 = parseFloat(line[4:23])
				v1 = parseFloat(line[23:42])
				v2 = parseFloat(line[42:61])
				v3 = parseFloat(line[61:80])
				lineCount += 1
				switch sys {
				case 'G', 'J', 'E', 'C':
					switch lineCount {
					case 1:
						eph.Iode = int(v0)
						eph.Crs = v1
						eph.DeltaN = v2
						eph.M0 = v3
					case 2:
						eph.Cuc = v0
						eph.Ecc = v1
						eph.Cus = v2
						eph.SqrtA = v3
					case 3:
						if sys == 'C' {
							eph.Toe = GTime{Week: eph.Toc.Week, Sec: v0 + 14} // The value of Week has not been read yet, so it is temporarily filled.
						} else {
							eph.Toe = GTime{Week: eph.Toc.Week, Sec: v0} // The value of Week has not been read yet, so it is temporarily filled.
						}
						eph.Cic = v1
						eph.Omega0 = v2
						eph.Cis = v3
					case 4:
						eph.I0 = v0
						eph.Crc = v1
						eph.Omega = v2
						eph.OmegaD = v3
					case 5:
						eph.Idot = v0
						eph.Code = int(v1)
						eph.Week = int(v2)
						if sys == 'C' {
							eph.Week += 1356 // BDT Week -> GPS Week
						}
						eph.Toe.Week = eph.Week
						if eph.Toe.Sec-eph.Toc.Sec < -302400 {
							eph.Toe.Sec += 604800
						} else if eph.Toe.Sec-eph.Toc.Sec > 302400 {
							eph.Toe.Sec -= 604800
						}
						eph.Flag = int(v3)
					case 6:
						if sys != 'E' {
							eph.Sva = getURAIndex(v0)
						} else {
							eph.Sva = getSISAIndex(v0)
						}
						eph.Svh = int(v1)
						eph.Tgd = v2
						eph.Iodc = int(v3)
						eph.Tgd2 = v3
					case 7:
						if sys == 'C' {
							eph.Tot = GTime{Week: eph.Week, Sec: v0 + 14}
						} else {
							eph.Tot = GTime{Week: eph.Week, Sec: v0}
						}
						if eph.Tot.Sec-eph.Toc.Sec < -302400 {
							eph.Tot.Sec += 604800
						} else if eph.Tot.Sec-eph.Toc.Sec > 302400 {
							eph.Tot.Sec -= 604800
						}
						switch sys {
						case 'G':
							eph.Fit = v1
						case 'J':
							if v1 == 0.0 {
								eph.Fit = 1
							} else {
								eph.Fit = 2
							}
						case 'C':
							eph.Iodc = int(v1)
						}
						nav[eph.Sat] = append(nav[eph.Sat], eph)
						// if _, ok := nav[eph.Sat]; ok {
						// 	nav[eph.Sat] = append(nav[eph.Sat], eph)
						// } else {
						// 	nav[eph.Sat] = []*Ephe{eph}
						// }
					}
				case 'R':
					switch lineCount {
					case 1:
						eph.PosX = v0 * 1000
						eph.VecX = v1 * 1000
						eph.AccX = v2 * 1000
						eph.Svh = int(v3)
					case 2:
						eph.PosY = v0 * 1000
						eph.VecY = v1 * 1000
						eph.AccY = v2 * 1000
						eph.FreqN = int(v3)
						if eph.FreqN > 128 {
							eph.FreqN -= 256
						}
					case 3:
						eph.PosZ = v0 * 1000
						eph.VecZ = v1 * 1000
						eph.AccZ = v2 * 1000
						eph.Age = int(v3)
						nav[eph.Sat] = append(nav[eph.Sat], eph)
						// if _, ok := nav[eph.Sat]; ok {
						// 	nav[eph.Sat] = append(nav[eph.Sat], eph)
						// } else {
						// 	nav[eph.Sat] = []*Ephe{eph}
						// }
					}
				case 'S':
					switch lineCount {
					case 1:
						eph.PosX = v0 * 1000
						eph.VecX = v1 * 1000
						eph.AccX = v2 * 1000
						eph.Svh = int(v3)
					case 2:
						eph.PosY = v0 * 1000
						eph.VecY = v1 * 1000
						eph.AccY = v2 * 1000
						eph.Sva = getURAIndex(v3)
					case 3:
						eph.PosZ = v0 * 1000
						eph.VecZ = v1 * 1000
						eph.AccZ = v2 * 1000
						eph.Iodn = int(v3)
						nav[eph.Sat] = append(nav[eph.Sat], eph)
						// if _, ok := nav[eph.Sat]; ok {
						// 	nav[eph.Sat] = append(nav[eph.Sat], eph)
						// } else {
						// 	nav[eph.Sat] = []*Ephe{eph}
						// }
					}
				}
			}

		}
	}

	// Check if reading completed without error
	if err := s.Err(); err != nil {
		return nil, err
	}

	// Sort by transmission time
	for k, v := range nav {
		if len(v) > 1 {
			sort.Slice(nav[k], func(i, j int) bool { return nav[k][i].Tot.Less(nav[k][j].Tot, false) })
		}
	}

	return &nav, nil
}

// Read real values by absorbing variations in exponential notation within RINEX files
func parseFloat(str string) float64 {
	s := strings.TrimSpace(str)
	if strings.ContainsAny(s, "Dd") {
		s = strings.Replace(s, "D", "E", 1)
		s = strings.Replace(s, "d", "e", 1)
	}
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

// Return URA index for specified value
func getURAIndex(x float64) int {
	if x > 0 && x <= 2.4 {
		return 0
	} else if x > 2.4 && x <= 3.4 {
		return 1
	} else if x > 3.4 && x <= 4.85 {
		return 2
	} else if x > 4.85 && x <= 6.85 {
		return 3
	} else if x > 6.85 && x <= 9.65 {
		return 4
	} else if x > 9.65 && x <= 13.65 {
		return 5
	} else if x > 13.65 && x <= 24.0 {
		return 6
	} else if x > 24.0 && x <= 48.0 {
		return 7
	} else if x > 48.0 && x <= 96.0 {
		return 8
	} else if x > 96.0 && x <= 192.0 {
		return 9
	} else if x > 192.0 && x <= 384.0 {
		return 10
	} else if x > 384.0 && x <= 768.0 {
		return 11
	} else if x > 768.0 && x <= 1536.0 {
		return 12
	} else if x > 1536.0 && x <= 3072.0 {
		return 13
	} else if x > 3072.0 && x <= 6144.0 {
		return 14
	} else {
		return 15
	}
}

// Return Galileo SISA index for specified value
func getSISAIndex(x float64) int {
	if x >= 0 && x <= 0.5 {
		return int(x / 0.01)
	} else if x > 0.5 && x <= 1.0 {
		return int((x-0.5)/0.02) + 50
	} else if x > 1.0 && x <= 2.0 {
		return int((x-1.0)/0.04) + 75
	} else if x > 2.0 && x <= 6.0 {
		return int((x-2.0)/0.16) + 100
	} else {
		return 255
	}
}
