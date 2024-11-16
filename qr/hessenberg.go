package main

import (
	"math/cmplx"
)

func hessenberg(A [][]complex128) [][]complex128 {
	n := len(A)
	H := make([][]complex128, n)
	for i := range H {
		H[i] = make([]complex128, n)
		copy(H[i], A[i])
	}

	for k := 0; k < n-2; k++ {
		u := make([]complex128, n-k-1)
		for i := 0; i < n-k-1; i++ {
			u[i] = H[k+i+1][k]
		}

		normU := cmplx.Sqrt(complex(norm(u), 0))

		if cmplx.Abs(normU) < 1e-15 {
			continue
		}

		if cmplx.Abs(u[0]) < 1e-15 {
			u[0] += 0
		} else {
			u[0] += (cmplx.Conj(u[0]) / complex(cmplx.Abs(u[0]), 0)) * normU
		}

        normalizingFactor := cmplx.Sqrt(complex(norm(u), 0))
		for i := 0; i < len(u); i++ {
			u[i] /= normalizingFactor
		}

		for j := k; j < n; j++ {
			sum := complex(0, 0)
			for i := 0; i < len(u); i++ {
				sum += cmplx.Conj(u[i]) * H[k+i+1][j]
			}

			for i := 0; i < len(u); i++ {
				H[k+i+1][j] -= 2 * u[i] * sum
			}
		}

		for i := 0; i < n; i++ {
			sum := complex(0, 0)
			for j := 0; j < len(u); j++ {
				sum += H[i][k+1+j] * u[j]
			}

			for j := 0; j < len(u); j++ {
				H[i][k+j+1] -= 2 * cmplx.Conj(u[j]) * sum
			}
		}
	}

	return H
}

