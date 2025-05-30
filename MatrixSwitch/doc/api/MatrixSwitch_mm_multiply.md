# Overview

The module `MatrixSwitch_mm_multiply` implements various subroutines for performing matrix-matrix multiplication, `C := alpha*op(A)*op(B) + beta*C`. It caters to different combinations of matrix storage formats (serial dense, serial sparse CSC/CSR, parallel dense block-cyclic, parallel sparse CSC), data types (real and complex), and operations on A and B (transpose, conjugate transpose). These specific implementations are invoked by the generic `mm_multiply` interface in the main `MatrixSwitch` module.

# Key Components

*   **Module `MatrixSwitch_mm_multiply`**: Houses the matrix-matrix multiplication implementations.
    *   **Serial Dense x Dense -> Dense**:
        *   `mm_multiply_sddenref(A, trA, B, trB, C, alpha, beta)`: Real.
        *   `mm_multiply_szdenref(A, tcA, B, tcB, C, alpha, beta)`: Complex.
    *   **Serial Sparse x Dense -> Dense**:
        *   `mm_multiply_sdcscsddensddenref(A, trA, B, trB, C, alpha, beta)`: A is CSC, B,C are Dense (Real).
        *   `mm_multiply_sddensdcscsddenref(A, trA, B, trB, C, alpha, beta)`: B is CSC, A,C are Dense (Real).
        *   `mm_multiply_sdcsrsddensddenref(A, trA, B, trB, C, alpha, beta)`: A is CSR, B,C are Dense (Real).
        *   `mm_multiply_sddensdcsrsddenref(A, trA, B, trB, C, alpha, beta)`: B is CSR, A,C are Dense (Real).
        *   `mm_multiply_szcscszdenszdenref(A, tcA, B, tcB, C, alpha, beta)`: A is CSC, B,C are Dense (Complex).
        *   `mm_multiply_szdenszcscszdenref(A, tcA, B, tcB, C, alpha, beta)`: B is CSC, A,C are Dense (Complex).
        *   `mm_multiply_szcsrszdenszdenref(A, tcA, B, tcB, C, alpha, beta)`: A is CSR, B,C are Dense (Complex).
        *   `mm_multiply_szdenszcsrszdenref(A, tcA, B, tcB, C, alpha, beta)`: B is CSR, A,C are Dense (Complex).
    *   **Serial Dense x Dense -> Sparse**:
        *   `mm_multiply_sddensddensdcscref(A, trA, B, trB, C, alpha, beta)`: C is CSC, A,B are Dense (Real).
        *   `mm_multiply_sddensddensdcsrref(A, trA, B, trB, C, alpha, beta)`: C is CSR, A,B are Dense (Real).
        *   `mm_multiply_szdenszdenszcscref(A, tcA, B, tcB, C, alpha, beta)`: C is CSC, A,B are Dense (Complex).
        *   `mm_multiply_szdenszdenszcsrref(A, tcA, B, tcB, C, alpha, beta)`: C is CSR, A,B are Dense (Complex).
    *   **Parallel (pspBLAS based, 1D distribution - "t1D")**: (Requires `HAVE_PSPBLAS` and `HAVE_MPI`)
        *   `mm_multiply_pddbcpdcscpddbct1D(A, B, trB, C, alpha, beta)`: A (pddbc), B (pdcsc), C (pddbc) - Real. op(A)=A.
        *   `mm_multiply_pdcscpddbcpddbct1D(A, trA, B, C, alpha, beta)`: A (pdcsc), B (pddbc), C (pddbc) - Real. op(B)=B.
        *   `mm_multiply_pddbcpddbcpdcsct1D(A, trA, B, trB, C, alpha, beta)`: A (pddbc), B (pddbc), C (pdcsc) - Real. op(A)!=op(B).
        *   `mm_multiply_pzdbcpzcscpzdbct1D(A, B, tcB, C, alpha, beta)`: A (pzdbc), B (pzcsc), C (pzdbc) - Complex. op(A)=A.
        *   `mm_multiply_pzcscpzdbcpzdbct1D(A, tcA, B, C, alpha, beta)`: A (pzcsc), B (pzdbc), C (pzdbc) - Complex. op(B)=B.
        *   `mm_multiply_pzdbcpzdbcpzcsct1D(A, tcA, B, tcB, C, alpha, beta)`: A (pzdbc), B (pzdbc), C (pzcsc) - Complex. op(A)!=op(B) (transpose sense).

# Important Variables/Constants

*   `A`, `B`: Input matrices.
*   `C`: Output matrix.
*   `alpha`, `beta`: Scalar multipliers.
*   `trA`, `trB` (logical): Transposition flags for real matrices A and B.
*   `tcA`, `tcB` (integer): Transposition/conjugation flags for complex matrices A and B (0: N, 1: C/H, 2: T).
*   `dp`: Double precision kind (from `MatrixSwitch_ops`).
*   `cmplx_0`: Complex zero (from `MatrixSwitch_ops`).
*   Matrix structure fields: `dval`, `zval`, `iaux3`, `iaux4` (for sparse CSR/CSC indices/pointers), `spm` (pspBLAS structure), `dim1`, `dim2`.
*   `HAVE_CONFIG_H`, `HAVE_PSPBLAS`, `HAVE_MPI`: Preprocessor macros.
*   `ms_mpi_size`, `ms_mpi_rank`, `ms_mpi_comm`: MPI related variables used in parallel routines.
*   `indxl2g`: External function (likely ScaLAPACK tool) to convert local index to global.

# Usage Examples

These subroutines are the workhorses for the `mm_multiply` interface in `MatrixSwitch.F90`. Users would typically call the interface:
```fortran
! C = 0.5 * A * B + 0.0 * C (effectively C = 0.5*A*B)
! Assuming A, B, C are TYPE(matrix) with compatible formats and dimensions.
CALL mm_multiply(A, 'N', B, 'N', C, 0.5_dp, 0.0_dp)
```
The `MatrixSwitch` module would then dispatch to one of the specific implementations in `MatrixSwitch_mm_multiply` based on the properties of A, B, and C.

# Dependencies and Interactions

*   **Uses `MatrixSwitch_ops`**: For `TYPE(matrix)` definition, `dp`, `cmplx_0`, and potentially other constants.
*   **MPI Library (`mpif.h`)**: Required for parallel implementations (`HAVE_MPI` and `HAVE_PSPBLAS` routines). Used for `mpi_bcast`, `mpi_reduce`.
*   **pspBLAS**: Parallel sparse implementations rely on the pspBLAS matrix structure (`spm`) within the `TYPE(matrix)`. The actual computation in these t1D routines is implemented directly using MPI communication, not via direct calls to pspBLAS multiplication routines.
*   **ScaLAPACK tools**: `indxl2g` is used in parallel "t1D" routines.
*   These routines form the core computational backend for the generic matrix multiplication interface. They handle the loops and logic for different storage schemes and transpose operations.
