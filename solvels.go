// Copyright (c) 2025 hitoshi.mukai.b@gmail.com. All rights reserved.
// You are free to use this source code for any purpose. The copyright remains with the author.
// The author accepts no liability for any damages arising from the use of this source code.
//
// Last modified: 2025.9.21
//

package gortk

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// Solve the observation equation using weighted least squares
// - dx = (G^t W G)^-1 G^t W dr
// - Return the error covariance matrix (G^t W G)^-1 as cov
func SolveLS(G mat.Matrix, dr mat.Vector, W mat.Matrix) (dx mat.Vector, cov mat.Matrix, err error) {

	n1, m1 := G.Dims()
	n2, m2 := W.Dims()
	if n1 != n2 {
		return nil, nil, fmt.Errorf("invalid matrix size. G^T(%d x %d), W(%d x %d)", m1, n1, n2, m2)
	}
	l1 := dr.Len()
	if l1 != m2 {
		return nil, nil, fmt.Errorf("invalid matrix size. W(%d x %d), dr(%d x 1)", n2, m2, l1)
	}

	// A（G^t W G)
	var WG mat.Dense
	WG.Mul(W, G)
	var A mat.Dense
	A.Mul(G.T(), &WG)
	// tol := 1e-10
	// rank := matrixRank(&A, tol)
	// fmt.Printf("Rank: %d\n", rank)

	// b（G^t W dr）
	var GtW mat.Dense
	GtW.Mul(G.T(), W)
	var b mat.VecDense
	b.MulVec(&GtW, dr)

	// Solve for x (x = A^-1 b)
	var x mat.VecDense
	err = x.SolveVec(&A, &b)
	if err != nil {
		return nil, nil, err
	}
	dx = &x

	// Set (G^T W G)^-1 as the covariance matrix
	var c mat.Dense
	err = c.Inverse(&A)
	if err != nil {
		return nil, nil, err
	}
	cov = &c

	return
}

// func matrixRank(A *mat.Dense, tol float64) int {
// 	var svd mat.SVD
// 	svd.Factorize(A, mat.SVDThin)

// 	// Retrieve singular values
// 	s := svd.Values(nil)

// 	// Count singular values that are greater than or equal to tol
// 	rank := 0
// 	for _, v := range s {
// 		if v > tol {
// 			rank++
// 		}
// 	}
// 	return rank
// }
