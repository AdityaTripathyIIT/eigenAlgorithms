package main

import (
    "fmt"
    "math/cmplx"
    "math"
)
func givens(A [][]complex128, k int) (complex128, complex128){
    c_k := complex(0, 0)
    s_k := complex(0, 0)
    f := cmplx.Abs(A[k][k])
    g := cmplx.Abs(A[k+1][k])
    c_k = -cmplx.Conj(A[k][k])/complex((math.Sqrt(f * f + g * g)), 0)
    s_k = -cmplx.Conj(A[k + 1][k])/complex((math.Sqrt(f * f + g * g)), 0)
    return c_k, s_k
}

func applyGivens(A [][]complex128) [][]complex128 {
    n := len(A)
    for k := 0; k < n - 1 ; k++{
        c_k , s_k := givens(A, k)
        for j := k; j < n ; j++{
            temp1 := c_k * A[k][j] + s_k * A[k + 1][j]
            temp2 := -cmplx.Conj(s_k) * A[k][j] + cmplx.Conj(c_k) * A[k + 1][j]
            A[k][j] = temp1
            A[k + 1][j] = temp2
        }
        fmt.Printf("\nK = %d\n", k)
        printMatrix(A)
    }
    for k := 0; k < n - 1; k++{
        c_k , s_k := givens(A, k)
        for j := 0; j < k + 1; j++{
            temp1 := A[j][k] * cmplx.Conj(c_k) + A[j][k + 1] * cmplx.Conj(s_k)
            temp2 :=  -A[j][k] * s_k + A[j][k + 1] * c_k
            A[j][k] = temp1
            A[j][k + 1] = temp2
        }
    }
    return A
}
