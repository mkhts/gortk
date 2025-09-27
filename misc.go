// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.21
//

package gortk

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"

	"gonum.org/v1/gonum/mat"
)

// ------------------------------------
// Mini functions
// ------------------------------------

func SQ(x float64) float64 {
	return x * x
}

func EucDist(a, b *PosXYZ) float64 {
	return math.Sqrt(SQ(a.X-b.X) + SQ(a.Y-b.Y) + SQ(a.Z-b.Z))
}

func DistDx(a, b *PosXYZ) float64 {
	return (b.X - a.X) / EucDist(a, b)
}

func DistDy(a, b *PosXYZ) float64 {
	return (b.Y - a.Y) / EucDist(a, b)
}

func DistDz(a, b *PosXYZ) float64 {
	return (b.Z - a.Z) / EucDist(a, b)
}

func ToDeg(rad float64) float64 {
	return rad / PI * 180.0
}

func ToRad(deg float64) float64 {
	return deg / 180.0 * PI
}

// ------------------------------------
// Debug print function
// ------------------------------------

func PrintMat(X mat.Matrix) {
	r, c := X.Dims()
	fmt.Fprintf(os.Stderr, "(%d x %d)\n", r, c)
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Fprintf(os.Stderr, "%v\n", fa)
}

func PrintA(format string, a ...any) {
	fmt.Fprintf(os.Stderr, format, a...)
}

func PrintAIf(cond bool, format string, a ...any) {
	if cond {
		PrintA(format, a...)
	}
}

func PrintB(t GTime, format string, a ...any) {
	fmt.Fprintf(os.Stderr, t.ToTime().UTC().Format("2006-01-02T15:04:05.000000")+"\t"+format, a...)
}

// Debug display level
var DBG_ int

// Debug display
func PrintD(v int, format string, a ...any) {
	PrintAIf(DBG_ >= v, format, a...)
}

func PrintE(err error) {
	fmt.Fprintf(os.Stderr, "err=%s\n", err.Error())
}

// ------------------------------------
// For command argument parsing
// ------------------------------------

type SysVar []SysType

func (p *SysVar) Set(s string) error {
	*p = []SysType{}
	for _, a := range strings.Split(s, ",") {
		*p = append(*p, SysType(a[:][0]))
	}
	return nil
}

func (p *SysVar) String() string {
	*p = []SysType{'G', 'J', 'E', 'R', 'C'}
	return ""
}

func (p *SysVar) Contains(s SysType) bool {
	for _, v := range *p {
		if s == v {
			return true
		}
	}
	return false
}

type SatVar []SatType

func (p *SatVar) Set(s string) error {
	*p = []SatType{}
	for _, a := range strings.Split(s, ",") {
		*p = append(*p, SatType(a))
	}
	return nil
}

func (p *SatVar) String() string {
	return ""
}

// Date and Time Parser (for command arguments)
type TimeStr time.Time

func (p *TimeStr) MarshalText() (text []byte, err error) {
	text, err = time.Time(*p).MarshalText()
	if err != nil {
		return nil, err
	}
	return text, nil
}

func (p *TimeStr) UnmarshalText(text []byte) error {
	s := string(text)
	t, err := time.Parse("2006/01/02 15:04:05", s)
	if err != nil {
		return err
	}
	*p = TimeStr(t)
	return nil
}
func NewTimeStr(t time.Time) *TimeStr {
	m := new(TimeStr)
	*m = TimeStr(t)
	return m
}

// Processing mode (0: SPP, 1: DGPS, 2: RTK)
type Mode int

const (
	SPP = iota
	DGPS
	RTK
)

func (p *Mode) Set(s string) error {
	i, err := strconv.ParseInt(s, 10, 0)
	if err != nil {
		return err
	}
	*p = Mode(i)
	return nil
}

func (p *Mode) String() string {
	switch *p {
	case SPP:
		return "SPP"
	case DGPS:
		return "DGPS"
	case RTK:
		return "RTK"
	default:
		return "UNKNOWN!"
	}
}

// ------------------------------------
// Others
// ------------------------------------

// Sort the list of satellite names
func Sorted(s []SatType) []SatType {
	s2 := make([]SatType, len(s))
	copy(s2, s)
	sort.Slice(s2, func(i, j int) bool {
		m := map[byte]int{'G': 0, 'J': 1, 'E': 2, 'R': 3, 'C': 4, 'S': 5}
		if m[s2[i][0]] == m[s2[j][0]] {
			return s2[i] < s2[j]
		} else {
			return m[s2[i][0]] < m[s2[j][0]]
		}
	})
	return s2
}

func SortedSatF(s []SatFType) []SatFType {
	s2 := make([]SatFType, len(s))
	copy(s2, s)
	sort.Slice(s2, func(i, j int) bool {
		m := map[byte]int{'G': 0, 'J': 1, 'E': 2, 'R': 3, 'C': 4, 'S': 5}
		if s2[i].F == s2[j].F {
			if s2[i].Sat.Sys() == s2[j].Sat.Sys() {
				return s2[i].Sat.Num() < s2[j].Sat.Num()
			} else {
				return m[byte(s2[i].Sat.Sys())] < m[byte(s2[j].Sat.Sys())]
			}
		} else {
			return s2[i].F < s2[j].F
		}
	})
	return s2
}

// Chi-squared test (α=0.001)
func ChiSqr(i int) float64 {
	v := [...]float64{
		10.8, 13.8, 16.3, 18.5, 20.5, 22.5, 24.3, 26.1, 27.9, 29.6,
		31.3, 32.9, 34.5, 36.1, 37.7, 39.3, 40.8, 42.3, 43.8, 45.3,
		46.8, 48.3, 49.7, 51.2, 52.6, 54.1, 55.5, 56.9, 58.3, 59.7,
		61.1, 62.5, 63.9, 65.2, 66.6, 68.0, 69.3, 70.7, 72.1, 73.4,
		74.7, 76.0, 77.3, 78.6, 80.0, 81.3, 82.6, 84.0, 85.4, 86.7,
		88.0, 89.3, 90.6, 91.9, 93.3, 94.7, 96.0, 97.4, 98.7, 100,
		101, 102, 103, 104, 105, 107, 108, 109, 110, 112,
		113, 114, 115, 116, 118, 119, 120, 122, 123, 125,
		126, 127, 128, 129, 131, 132, 133, 134, 135, 137,
		138, 139, 140, 142, 143, 144, 145, 147, 148, 149}
	if i < len(v) {
		return v[i]
	} else {
		return 0
	}
}
