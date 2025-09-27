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
	"strconv"
	"strings"
)

//-------------------------------------------------------------------
// PosLLH
//-------------------------------------------------------------------

type PosLLH struct {
	Lat float64
	Lon float64
	Hei float64
}

func NewPosLLH(lat, lon, hei float64) *PosLLH {
	return &PosLLH{
		Lat: lat,
		Lon: lon,
		Hei: hei,
	}
}

func (llh *PosLLH) ToXYZ() PosXYZ {
	// Ellipsoid parameters
	f := Fe                     // Flattening
	a := Re                     // Semi-major axis
	e := math.Sqrt(f * (2 - f)) // Eccentricity

	// Conversion to Cartesian coordinates
	n := a / math.Sqrt(1-e*e*math.Sin(llh.Lat)*math.Sin(llh.Lat))
	return PosXYZ{
		X: (n + llh.Hei) * math.Cos(llh.Lat) * math.Cos(llh.Lon),
		Y: (n + llh.Hei) * math.Cos(llh.Lat) * math.Sin(llh.Lon),
		Z: (n*(1-e*e) + llh.Hei) * math.Sin(llh.Lat),
	}
}

func (llh *PosLLH) ToENU(base PosXYZ) PosENU {
	xyz := llh.ToXYZ()
	return xyz.ToENU(base)
}

func (usr *PosLLH) Elevation(sat PosXYZ) float64 {
	// Convert to ENU coordinates to calculate elevation angle
	xyz := usr.ToXYZ()
	enu := sat.ToENU(xyz)
	return math.Atan2(enu.U, math.Sqrt(enu.E*enu.E+enu.N*enu.N))
}

func (usr *PosLLH) Azimuth(sat PosXYZ) float64 {
	// Convert to ENU coordinates to calculate azimuth
	xyz := usr.ToXYZ()
	enu := sat.ToENU(xyz)
	return math.Atan2(enu.E, enu.N)
}

// Read from string
func (llh *PosLLH) Set(s string) error {
	var err error
	f := strings.Fields(s)
	llh.Lat, err = strconv.ParseFloat(f[0], 64)
	if err != nil {
		return err
	}
	llh.Lon, err = strconv.ParseFloat(f[1], 64)
	if err != nil {
		return err
	}
	llh.Hei, err = strconv.ParseFloat(f[2], 64)
	if err != nil {
		return err
	}
	llh.Lat *= math.Pi / 180
	llh.Lon *= math.Pi / 180
	return nil
}

// Convert to string
func (llh *PosLLH) String() string {
	return fmt.Sprintf("%.8f %.8f %.4f", llh.Lat, llh.Lon, llh.Hei)
}

//-------------------------------------------------------------------
// PosXYZ
//-------------------------------------------------------------------

type PosXYZ struct {
	X float64
	Y float64
	Z float64
}

func NewPosXYZ(x, y, z float64) *PosXYZ {
	return &PosXYZ{
		X: x,
		Y: y,
		Z: z,
	}
}

func (pos *PosXYZ) ToLLH() PosLLH {
	// In case of origin
	if pos.X == 0 && pos.Y == 0 && pos.Z == 0 {
		return PosLLH{Lat: 0, Lon: 0, Hei: -Re}
	}

	// Ellipsoid parameters
	f := Fe                     // Flattening
	a := Re                     // Semi-major axis
	b := a * (1 - f)            // Semi-minor axis
	e := math.Sqrt(f * (2 - f)) // Eccentricity

	// Parameters for coordinate transformation
	h := a*a - b*b
	p := math.Sqrt(pos.X*pos.X + pos.Y*pos.Y)
	t := math.Atan2(pos.Z*a, p*b)
	sint := math.Sin(t)
	cost := math.Cos(t)

	// Conversion to latitude and longitude
	lat := math.Atan2(pos.Z+h/b*sint*sint*sint, p-h/a*cost*cost*cost)
	lon := math.Atan2(pos.Y, pos.X)
	n := a / math.Sqrt(1-e*e*math.Sin(lat)*math.Sin(lat)) // Radius of curvature in the prime vertical
	hei := p/math.Cos(lat) - n
	return PosLLH{Lat: lat, Lon: lon, Hei: hei}
}

func (pos *PosXYZ) ToENU(base PosXYZ) PosENU {
	// Relative position from the reference location
	x := pos.X - base.X
	y := pos.Y - base.Y
	z := pos.Z - base.Z

	// Latitude and longitude of the reference location
	llh := base.ToLLH()
	s1 := math.Sin(llh.Lon)
	c1 := math.Cos(llh.Lon)
	s2 := math.Sin(llh.Lat)
	c2 := math.Cos(llh.Lat)

	// Rotate the relative position to convert to ENU coordinates
	return PosENU{
		E: -x*s1 + y*c1,
		N: -x*c1*s2 - y*s1*s2 + z*c2,
		U: x*c1*c2 + y*s1*c2 + z*s2,
	}
}

func (usr *PosXYZ) Elevation(sat PosXYZ) float64 {
	// Convert to ENU coordinates to calculate elevation angle
	enu := sat.ToENU(*usr)
	return math.Atan2(enu.U, math.Sqrt(enu.E*enu.E+enu.N*enu.N))
}

func (usr *PosXYZ) Azimuth(sat PosXYZ) float64 {
	// Convert to ENU coordinates to calculate azimuth
	enu := sat.ToENU(*usr)
	return math.Atan2(enu.E, enu.N)
}

//-------------------------------------------------------------------
// PosENU
//-------------------------------------------------------------------

type PosENU struct {
	E float64
	N float64
	U float64
}

func NewPosENU(e, n, u float64) *PosENU {
	return &PosENU{
		E: e,
		N: n,
		U: u,
	}
}

func (enu *PosENU) ToXYZ(base PosXYZ) PosXYZ {
	// Latitude and longitude of the reference location
	llh := base.ToLLH()
	s1 := math.Sin(llh.Lon)
	c1 := math.Cos(llh.Lon)
	s2 := math.Sin(llh.Lat)
	c2 := math.Cos(llh.Lat)

	// Rotate the ENU coordinates to convert to relative position
	x := -enu.E*s1 - enu.N*c1*s2 + enu.U*c1*c2
	y := enu.E*c1 - enu.N*s1*s2 + enu.U*s1*c2
	z := enu.N*c2 + enu.U*s2

	// Add to the reference location
	x += base.X
	y += base.Y
	z += base.Z
	return PosXYZ{
		X: x,
		Y: y,
		Z: z,
	}
}

func (enu *PosENU) Elevation() float64 {
	return math.Atan2(enu.U, math.Sqrt(enu.E*enu.E+enu.N*enu.N))
}

func (enu *PosENU) Azimuth() float64 {
	return math.Atan2(enu.E, enu.N)
}
