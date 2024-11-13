package main

import (
	"fmt"
	"math/cmplx"
)

// Compute the squared norm of a complex vector (sum of squared absolute values)
func norm(v []complex128) float64 {
	var sum float64
	for _, val := range v {
		sum += cmplx.Abs(val) * cmplx.Abs(val)
	}
	return sum
}

// Hessenberg transformation using Householder reflections
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
			u[0] += normU
		} else {
			u[0] += cmplx.Conj(u[0]) / complex(cmplx.Abs(u[0]), 0) * normU
		}

		//beta := complex(cmplx.Abs(u[0]), 0)
        normalizingFactor := cmplx.Sqrt(complex(norm(u), 0))
		for i := 0; i < len(u); i++ {
			u[i] /= normalizingFactor
		}

		for j := k; j < n; j++ {
			sum := complex(0, 0)
			for i := range u {
				sum += cmplx.Conj(u[i]) * H[k+1+i][j]
			}

			for i := range u {
				H[k+1+i][j] -= 2 * u[i] * sum
			}
		}

		for i := 0; i < n; i++ {
			// Compute the scalar product (H[i][k+1:n] * u)
			sum := complex(0, 0)
			for j := range u {
				sum += H[i][k+1+j] * u[j]
			}

			for j := range u {
				H[i][k+1+j] -= 2 * cmplx.Conj(u[j]) * sum
			}
		}
	}

	return H
}

// Print a complex matrix in a readable format
func printMatrix(matrix [][]complex128) {
	for _, row := range matrix {
		for _, value := range row {
			fmt.Printf("%6.2f +%6.2fi ", real(value), imag(value))
		}
		fmt.Println()
	}
}

func main() {
	matrix := [][]complex128{
		{complex(4, 0), complex(1, 0), complex(3, 0), complex(4, 0)},
		{complex(1, 0), complex(2, 0), complex(4, 0), complex(0, 0)},
		{complex(3, 0), complex(4, 0), complex(5, 0), complex(1, 0)},
		{complex(4, 0), complex(0, 0), complex(1, 0), complex(4, 0)},
	}

	fmt.Println("Original Matrix:")
	printMatrix(matrix)

	hessenbergMatrix := hessenberg(matrix)

	fmt.Println("Upper Hessenberg form:")
	printMatrix(hessenbergMatrix)
}

