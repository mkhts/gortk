// This code is adapted from RTKLIB.
// The author gratefully acknowledges T.Takasu for his outstanding contribution in developing RTKLIB.
//
// Last modified: 2025.9.14
//

package gortk

import (
	"math"
)

func TropModel(pos *PosXYZ) float64 {
	const TEMP0 = 15.0
	const ZELEV = 90.0 * math.Pi / 180.0
	const HUMI = 0.0
	llh := pos.ToLLH()
	if llh.Hei < -100.0 || 1e4 < llh.Hei {
		return 0.0
	}
	if llh.Hei < 0.0 {
		llh.Hei = 0.0
	}
	hgt := llh.Hei
	pres := 1013.25 * math.Pow(1.0-2.2557e-5*hgt, 5.2568)
	temp := TEMP0 - 6.5e-3*hgt + 273.16
	e := 6.108 * HUMI * math.Exp((17.15*temp-4684.0)/(temp-38.45))
	// saastamoinen model
	z := math.Pi/2.0 - ZELEV
	trph := 0.0022768 * pres / (1.0 - 0.00266*math.Cos(2.0*llh.Lat) - 0.00028*hgt/1e3) / math.Cos(z)
	trpw := 0.002277 * (1255.0/temp + 0.05) * e / math.Cos(z)
	return trph + trpw
}

func TropMapf(time *GTime, pos *PosXYZ, elev float64) float64 {
	llh := pos.ToLLH()
	if llh.Hei < -1000.0 || llh.Hei > 20000.0 {
		return 0.0
	}
	return nmf(time, &llh, elev)
}

func nmf(time *GTime, pos *PosLLH, elev float64) float64 {
	coef := [9][5]float64{
		{1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3},
		{2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3},
		{62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3},

		{0.0000000e-0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5},
		{0.0000000e-0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5},
		{0.0000000e-0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5},

		{5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4},
		{1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3},
		{4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2},
	}
	aht := [3]float64{2.53e-5, 5.49e-3, 1.14e-3}
	lat := ToDeg(pos.Lat)
	hgt := pos.Hei
	if elev <= 0.0 {
		return 0.0
	}
	doy := time.ToTime().YearDay()
	y := (float64(doy) - 28.0) / 365.25
	if lat < 0.0 {
		y += 0.5
	}
	cosy := math.Cos(2 * math.Pi * y)
	lat = math.Abs(lat)
	ah := [3]float64{}
	for i := range 3 {
		ah[i] = interpc(coef[i], lat) - interpc(coef[i+3], lat)*cosy
	}
	dm := (1.0/math.Sin(elev) - mapf(elev, aht[0], aht[1], aht[2])) * hgt / 1e3
	return mapf(elev, ah[0], ah[1], ah[2]) + dm
}

func interpc(coef [5]float64, lat float64) float64 {
	i := int(lat / 15.0)
	if i < 1 {
		return coef[0]
	} else if i > 4 {
		return coef[4]
	}
	return coef[i-1]*(1.0-lat/15.0+float64(i)) + coef[i]*(lat/15.0-float64(i))
}

func mapf(el, a, b, c float64) float64 {
	sinel := math.Sin(el)
	return (1.0 + a/(1.0+b/(1.0+c))) / (sinel + (a / (sinel + b/(sinel+c))))
}
