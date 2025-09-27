// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.20
//

package gortk

import (
	"math"
)

// Calculate satellite position from receiver data reception time and pseudorange
func SatPos(e *Ephe, rcvt GTime, psr float64) (xyz PosXYZ) {
	dOMGe := 7.2921151467e-5 // Earth rotation angular velocity [rad/s]
	Mue := 3.986005e14       // Earth gravitational constant [m^3/s^2]
	if e.Sat.Sys() == 'E' {
		dOMGe = 7.2921151467e-5
		Mue = 3.986004418e14
	} else if e.Sat.Sys() == 'C' {
		dOMGe = 7.292115e-5
		Mue = 3.986004418e14
	}
	switch {
	case e.Sat.Sys() == 'G' || e.Sat.Sys() == 'J' || e.Sat.Sys() == 'E' || e.Sat.Sys() == 'C':
		tk0 := rcvt.ToTime().Sub(e.Toe.ToTime()).Seconds()
		tk := tk0 - psr/C
		n := math.Sqrt(Mue)/e.SqrtA/e.SqrtA/e.SqrtA + e.DeltaN
		mk := e.M0 + n*tk
		ek := mk
		for i := 0; i < 10; i++ {
			ek = mk + e.Ecc*math.Sin(ek)
		}
		rk := e.SqrtA * e.SqrtA * (1 - e.Ecc*math.Cos(ek))
		vk := math.Atan2(math.Sqrt(1-e.Ecc*e.Ecc)*math.Sin(ek), math.Cos(ek)-e.Ecc)
		pk := vk + e.Omega
		d_uk := e.Cus*math.Sin(2*pk) + e.Cuc*math.Cos(2*pk)
		d_rk := e.Crs*math.Sin(2*pk) + e.Crc*math.Cos(2*pk)
		d_ik := e.Cis*math.Sin(2*pk) + e.Cic*math.Cos(2*pk)
		uk := pk + d_uk
		rk = rk + d_rk
		ik := e.I0 + d_ik + e.Idot*tk
		xk := rk * math.Cos(uk)
		yk := rk * math.Sin(uk)
		omk := e.Omega0 + (e.OmegaD-dOMGe)*tk0 - dOMGe*e.Toe.Sec // Including Sagnac effect
		//omk := e.Omega0 + (e.OmegaD-dOMGe)*tk - dOMGe*e.Toe.Sec // Not including Sagnac effect
		if e.Sat.Sys() == 'C' {
			omk = e.Omega0 + (e.OmegaD-dOMGe)*tk0 - dOMGe*(e.Toe.Sec-14) // Including Sagnac effect
			//omk = e.Omega0 + (e.OmegaD-dOMGe)*tk - dOMGe*(e.Toe.Sec-14) // Not including Sagnac effect
		}
		xyz.X = xk*math.Cos(omk) - yk*math.Sin(omk)*math.Cos(ik)
		xyz.Y = xk*math.Sin(omk) + yk*math.Cos(omk)*math.Cos(ik)
		xyz.Z = yk * math.Sin(ik)
		if e.Sat.Sys() == 'C' && (e.Sat.Num() <= 5 || e.Sat.Num() >= 59) { // Beidou geostationary
			omk = e.Omega0 + e.OmegaD*tk0 - dOMGe*(e.Toe.Sec-14)
			xg := xk*math.Cos(omk) - yk*math.Sin(omk)*math.Cos(ik)
			yg := xk*math.Sin(omk) + yk*math.Cos(omk)*math.Cos(ik)
			zg := yk * math.Sin(ik)
			sino := math.Sin(dOMGe * tk0)
			coso := math.Cos(dOMGe * tk0)
			cos5 := math.Cos(-5 * math.Pi / 180.0)
			sin5 := math.Sin(-5 * math.Pi / 180.0)
			xyz.X = xg*coso + yg*sino*cos5 + zg*sino*sin5
			xyz.Y = -xg*sino + yg*coso*cos5 + zg*coso*sin5
			xyz.Z = -yg*sin5 + zg*cos5
		}
	case e.Sat.Sys() == 'R':
		fallthrough
	case e.Sat.Sys() == 'S':
		tk0 := rcvt.ToTime().Sub(e.Toe.ToTime()).Seconds()
		tk := tk0 - psr/C
		var x [6]float64
		x[0], x[1], x[2] = e.PosX, e.PosY, e.PosZ
		x[3], x[4], x[5] = e.VecX, e.VecY, e.VecZ
		var acc [3]float64
		acc[0], acc[1], acc[2] = e.AccX, e.AccY, e.AccZ
		const TSTEP = 60.0
		tt := TSTEP
		if tk < 0 {
			tt = -TSTEP
		}
		for math.Abs(tk) > 1e-9 {
			if math.Abs(tk) < TSTEP {
				tt = tk
			}
			glorbit(tt, &x, acc)
			tk -= tt
		}
		// With Sagnac effect
		omk := dOMGe * psr / C
		xyz.X = x[0]*math.Cos(omk) + x[1]*math.Sin(omk)
		xyz.Y = -x[0]*math.Sin(omk) + x[1]*math.Cos(omk)
		xyz.Z = x[2]
		// Without Sagnac effect
		//xyz.X = x[0]
		//xyz.Y = x[1]
		//xyz.Z = x[2]
	}
	return
}

// Function used when calculating GLONASS satellite position (1)
func deq(x [6]float64, xdot *[6]float64, acc [3]float64) {
	const dOMGeR = 7.292115e-5 // Earth rotation angular velocity [rad/s] for GLONASS
	const OMG2 = dOMGeR * dOMGeR
	const J2_GLO = 1.0826257e-3
	const MU_GLO = 3.9860044e14
	const RE_GLO = 6378136.0

	r2 := x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
	r3 := r2 * math.Sqrt(r2)
	if r2 <= 0 {
		xdot[0], xdot[1], xdot[2], xdot[3], xdot[4], xdot[5] = 0, 0, 0, 0, 0, 0
		return
	}
	a := 1.5 * J2_GLO * MU_GLO * (RE_GLO * RE_GLO) / r2 / r3
	b := 5.0 * x[2] * x[2] / r2
	c := -MU_GLO/r3 - a*(1.0-b)
	xdot[0] = x[3]
	xdot[1] = x[4]
	xdot[2] = x[5]
	xdot[3] = (c+OMG2)*x[0] + 2.0*dOMGeR*x[4] + acc[0]
	xdot[4] = (c+OMG2)*x[1] - 2.0*dOMGeR*x[3] + acc[1]
	xdot[5] = (c-2.0*a)*x[2] + acc[2]
}

// Function used when calculating GLONASS satellite position (2)
func glorbit(t float64, x *[6]float64, acc [3]float64) {
	var k1, k2, k3, k4, w [6]float64
	deq(*x, &k1, acc)
	for i := 0; i < 6; i++ {
		w[i] = x[i] + k1[i]*t/2.0
	}
	deq(w, &k2, acc)
	for i := 0; i < 6; i++ {
		w[i] = x[i] + k2[i]*t/2.0
	}
	deq(w, &k3, acc)
	for i := 0; i < 6; i++ {
		w[i] = x[i] + k3[i]*t
	}
	deq(w, &k4, acc)
	for i := 0; i < 6; i++ {
		x[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) * t / 6.0
	}
}
