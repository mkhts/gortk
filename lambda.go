// This code is adapted from RTKLIB.
// The author gratefully acknowledges T.Takasu for his outstanding contribution in developing RTKLIB.
//
// Last modified: 2023.12.9
//

package gortk

import (
	"fmt"
	"math"
)

func sgn(x float64) float64 {
	if x <= 0.0 {
		return -1.0
	} else {
		return 1.0
	}
}

func ld(n int, Q []float64, L []float64, D []float64, info *int) error {
	A := make([]float64, n*n)
	for i := range A {
		A[i] = Q[i]
	}
	*info = 0
	for i := n - 1; i > -1; i-- {
		D[i] = A[i+i*n]
		if D[i] <= 0 {
			*info = -1
			break
		}
		a := math.Sqrt(D[i])
		for j := 0; j < i+1; j++ {
			L[i+j*n] = A[i+j*n] / a
		}
		for j := 0; j < i; j++ {
			for k := 0; k < j+1; k++ {
				A[j+k*n] -= L[i+k*n] * L[i+j*n]
			}
		}
		for j := 0; j < i+1; j++ {
			L[i+j*n] /= L[i+i*n]
		}
	}
	if *info != 0 {
		return fmt.Errorf("LD factorization error")
	}
	return nil
}

func gauss(n int, L, Z []float64, i, j int) {
	mu := int(math.Round(L[i+j*n]))
	if mu != 0 {
		for k := i; k < n; k++ {
			L[k+n*j] -= float64(mu) * L[k+i*n]
		}
		for k := 0; k < n; k++ {
			Z[k+n*j] -= float64(mu) * Z[k+i*n]
		}
	}
}

func perm(n int, L, D []float64, j int, Del float64, Z []float64) {
	eta := D[j] / Del
	lam := D[j+1] * L[j+1+j*n] / Del
	D[j] = eta * D[j+1]
	D[j+1] = Del
	for k := 0; k < j; k++ {
		a0 := L[j+k*n]
		a1 := L[j+1+k*n]
		L[j+k*n] = -L[j+1+j*n]*a0 + a1
		L[j+1+k*n] = eta*a0 + lam*a1
	}
	L[j+1+j*n] = lam
	for k := j + 2; k < n; k++ {
		L[k+j*n], L[k+(j+1)*n] = L[k+(j+1)*n], L[k+j*n]
	}
	for k := 0; k < n; k++ {
		Z[k+j*n], Z[k+(j+1)*n] = Z[k+(j+1)*n], Z[k+j*n]
	}
}

func reduction(n int, L, D, Z []float64) {
	j := n - 2
	k := n - 2
	for j >= 0 {
		if j <= k {
			for i := j + 1; i < n; i++ {
				gauss(n, L, Z, i, j)
			}
		}
		Del := D[j] + L[j+1+j*n]*L[j+1+j*n]*D[j+1]
		if (Del + 1e-6) < D[j+1] {
			perm(n, L, D, j, Del, Z)
			k = j
			j = n - 2
		} else {
			j -= 1
		}
	}
}

func matmul(tr string, n, k, m int, alpha float64, A, B []float64, beta float64, C []float64) {
	var f int
	if tr[0] == 'N' {
		if tr[1] == 'N' {
			f = 1
		} else {
			f = 2
		}
	} else {
		if tr[1] == 'N' {
			f = 3
		} else {
			f = 4
		}
	}
	for i := 0; i < n; i++ {
		for j := 0; j < k; j++ {
			d := 0.0
			if f == 1 {
				for x := 0; x < m; x++ {
					d += A[i+x*n] * B[x+j*m]
				}
			} else if f == 2 {
				for x := 0; x < m; x++ {
					d += A[i+x*n] * B[j+x*k]
				}
			} else if f == 3 {
				for x := 0; x < m; x++ {
					d += A[x+i*m] * B[x+j*m]
				}
			} else if f == 4 {
				for x := 0; x < m; x++ {
					d += A[x+i*m] * B[j+x*k]
				}
			}
			if beta == 0.0 {
				C[i+j*n] = alpha * d
			} else {
				C[i+j*n] = alpha*d + beta*C[i+j*n]
			}
		}
	}
}

func search(n, m int, L, D []float64, zs, zn, s []float64, info *int) error {
	const LOOPMAX = 2000000
	nn := 0
	imax := 0
	maxdist := 1e99
	S := make([]float64, n*n)
	dist := make([]float64, n)
	zb := make([]float64, n)
	z := make([]float64, n)
	step := make([]float64, n)
	k := n - 1
	dist[k] = 0.0
	zb[k] = zs[k]
	z[k] = math.Round(zb[k])
	y := zb[k] - z[k]
	step[k] = sgn(y)
	c := 0
	for c = 0; c < LOOPMAX; c++ {
		newdist := dist[k] + y*y/D[k]
		if newdist < maxdist {
			if k != 0 {
				k -= 1
				dist[k] = newdist
				for i := 0; i < k+1; i++ {
					S[k+i*n] = S[k+1+i*n] + (z[k+1]-zb[k+1])*L[k+1+i*n]
				}
				zb[k] = zs[k] + S[k+k*n]
				z[k] = math.Round(zb[k])
				y = zb[k] - z[k]
				step[k] = sgn(y)
			} else {
				if nn < m {
					if nn == 0 || newdist > s[imax] {
						imax = nn
					}
					for i := 0; i < n; i++ {
						zn[i+nn*n] = z[i]
					}
					s[nn] = newdist
					nn += 1
				} else {
					if newdist < s[imax] {
						for i := 0; i < n; i++ {
							zn[i+imax*n] = z[i]
						}
						s[imax] = newdist
						imax = 0
						for i := 0; i < m; i++ {
							if s[imax] < s[i] {
								imax = i
							}
						}
					}
					maxdist = s[imax]
				}
				z[0] += step[0]
				y = zb[0] - z[0]
				step[0] = -step[0] - sgn(step[0])
			}
		} else {
			if k == n-1 {
				break
			} else {
				k += 1
				z[k] += step[k]
				y = zb[k] - z[k]
				step[k] = -step[k] - sgn(step[k])
			}
		}
	}
	for i := 0; i < m-1; i++ {
		for j := i + 1; j < m; j++ {
			if s[i] < s[j] {
				continue
			}
			s[i], s[j] = s[j], s[i]
			for k := 0; k < n; k++ {
				zn[k+i*n], zn[k+j*n] = zn[k+j*n], zn[k+i*n]
			}
		}
	}
	if c >= LOOPMAX {
		*info = -1
		return fmt.Errorf("search loop count overflow")
	}
	*info = 0
	return nil
}

func ludcmp(A []float64, n int, indx []int) int {
	vv := make([]float64, n)
	imax := 0
	d := 1.0
	for i := 0; i < n; i++ {
		big := 0.0
		for j := 0; j < n; j++ {
			tmp := math.Abs(A[i+j*n])
			if tmp > big {
				big = tmp
			}
		}
		if big > 0.0 {
			vv[i] = 1.0 / big
		} else {
			return -1
		}
	}
	for j := 0; j < n; j++ {
		for i := 0; i < j; i++ {
			s := A[i+j*n]
			for k := 0; k < i; k++ {
				s -= A[i+k*n] * A[k+j*n]
			}
			A[i+j*n] = s
		}
		big := 0.0
		for i := j; i < n; i++ {
			s := A[i+j*n]
			for k := 0; k < j; k++ {
				s -= A[i+k*n] * A[k+j*n]
			}
			A[i+j*n] = s
			tmp := vv[i] * math.Abs(s)
			if tmp >= big {
				big = tmp
				imax = i
			}
		}
		if j != imax {
			for k := 0; k < n; k++ {
				tmp := A[imax+k*n]
				A[imax+k*n] = A[j+k*n]
				A[j+k*n] = tmp
			}
			d = -d
			vv[imax] = vv[j]
		}
		indx[j] = imax
		if A[j+j*n] == 0.0 {
			return -1
		}
		if j != n-1 {
			tmp := 1.0 / A[j+j*n]
			for i := j + 1; i < n; i++ {
				A[i+j*n] *= tmp
			}
		}
	}
	return 0
}

func lubksb2(A []float64, n int, indx []int, b []float64, k int) {
	ii := -1
	for i := 0; i < n; i++ {
		ip := indx[i]
		s := b[k+ip]
		b[k+ip] = b[k+i]
		if ii >= 0 {
			for j := ii; j < i; j++ {
				s -= A[i+j*n] * b[k+j]
			}
		} else if s != 0 {
			ii = i
		}
		b[k+i] = s
	}
	for i := n - 1; i > -1; i-- {
		s := b[k+i]
		for j := i + 1; j < n; j++ {
			s -= A[i+j*n] * b[k+j]
		}
		b[k+i] = s / A[i+i*n]
	}
}

func matinv(A []float64, n int) int {
	indx := make([]int, n)
	B := make([]float64, n*n)
	for i := 0; i < n*n; i++ {
		B[i] = A[i]
	}
	if ludcmp(B, n, indx) != 0 {
		return -1
	}
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			A[i+j*n] = 0.0
		}
		A[j+j*n] = 1.0
		lubksb2(B, n, indx, A, j*n)
	}
	return 0
}

func solve(tr string, A, Y []float64, n, m int, X []float64, info *int) error {
	B := make([]float64, n*n)
	for i := 0; i < n*n; i++ {
		B[i] = A[i]
	}
	*info = matinv(B, n)
	if *info == 0 {
		if tr[0] == 'N' {
			matmul("NN", n, m, n, 1.0, B, Y, 0.0, X)
		} else {
			matmul("TN", n, m, n, 1.0, B, Y, 0.0, X)
		}
	}
	return nil
}

func LAMBDA(n, m int, a, Q, F, s []float64) error {
	if n <= 0 || m <= 0 {
		return fmt.Errorf("n <= 0 || m <= 0")
	}
	L := make([]float64, n*n)
	D := make([]float64, n)
	Z := make([]float64, n*n)
	for i := 0; i < n; i++ {
		Z[i+i*n] = 1.0
	}
	z := make([]float64, n)
	E := make([]float64, n*m)
	var info int
	err := ld(n, Q, L, D, &info)
	if err != nil {
		return err
	}
	if info == 0 {
		reduction(n, L, D, Z)
		matmul("TN", n, 1, n, 1.0, Z, a, 0.0, z)
		err = search(n, m, L, D, z, E, s, &info)
		if err != nil {
			return err
		}
		if info == 0 {
			err = solve("T", Z, E, n, m, F, &info)
			if err != nil {
				return err
			}
		}
	}
	return nil
}
