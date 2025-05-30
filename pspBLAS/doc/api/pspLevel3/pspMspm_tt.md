# Overview

The `pspMspm_tt.F90` file provides the `pspMspm_tt` module, a specialized component of the pspBLAS Level 3 library. This module handles the multiplication of a transposed dense matrix `A` (`op(A)=A^T` or `A^H`) by a transposed sparse matrix `B` (`op(B)=B^T` or `B^H`). The operation is formulated as `C := alpha*op(A)*op(B) + beta*C`, where `C` is a dense matrix. This module contains the "_tt" (Transpose-Transpose) specific implementations for the `psp_gemspm` (General Matrix - Sparse Matrix product) functionality.

The current implementation strategy for this double-transpose case involves:
1.  Explicitly transposing the dense matrix `A` (i.e., computing `A^T` or `A^H`) into a temporary distributed dense matrix.
2.  Calling the `psp_gemspm_nt` routines (which handle Normal-Transpose products) using the explicitly transposed `A` (now effectively in 'Normal' form) and the original sparse matrix `B` with its transpose operator `opB`.

A "TODO" comment in the source code indicates that this approach (explicitly transposing `A`) could be optimized, for example, by implementing a method that directly handles `op(A)*op(B)` or by providing an on-the-fly transpose capability for sparse matrices.

# Key Components

*   **Module `pspMspm_tt`**:
    *   **Public Interface `psp_gemspm_tt`**: This interface exposes the routines for multiplying a transposed dense matrix by a transposed sparse matrix.
        *   `psp_dgemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Dense matrix `A` is globally K-by-M (so `op(A)` is M-by-K). Sparse matrix `B` (of `TYPE(psp_matrix_spm)`) is globally N-by-K (so `op(B)` is K-by-N). The resulting dense matrix `C` is M-by-N.
            *   `A`: A 2D Fortran array for the local part of the distributed K-by-M dense input matrix `A`.
            *   `opA`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `B`: A variable of `TYPE(psp_matrix_spm)`, representing the local part of the distributed N-by-K sparse input matrix `B`.
            *   `opB`: Character flag, must be 'T' or 'C' for these routines.
            *   `C`: A 2D Fortran array (inout), the local part of the distributed M-by-N dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm**:
            1.  **Grid Information**: Fetches BLACS process grid details using `blacs_gridinfo` and `psp_icontxt`.
            2.  **Explicit Transpose of Dense Matrix A**:
                *   Determines local dimensions for the input matrix `A` (KxM) and the temporary transposed matrix `tmp` (MxK) using `numroc`.
                *   Initializes BLACS array descriptors (`desc_before` for `A`, `desc_after` for `tmp`).
                *   Allocates memory for the temporary distributed matrix `tmp`.
                *   Performs a parallel transpose of `A` into `tmp` using ScaLAPACK routines:
                    *   `pdtran` for real matrix `A`.
                    *   `pztranc` (for `opA='C'`) or `pztranu` (for `opA='T'`) for complex matrix `A`.
                The matrix `tmp` now holds `op(A)`.
            3.  **Call `psp_gemspm_nt`**: Invokes the "Normal-Transpose" routine `psp_gemspm_nt` (from the `pspMspm_nt` module).
                *   `tmp` (which is `op(A)`) is passed as the non-transposed dense matrix argument.
                *   The original sparse matrix `B` and its transpose operator `opB` are passed as is.
                *   The call effectively becomes: `C = alpha * tmp * op(B) + beta * C`, which is equivalent to `C = alpha * op(A) * op(B) + beta * C`.
            4.  **Cleanup**: Deallocates the temporary matrix `tmp`.
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Provides essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes).
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `tmp`: An allocatable 2D array used to store the explicit transpose of the dense matrix `A`.
*   `desc_before`, `desc_after`: BLACS array descriptors used for the `pdtran`/`pztranc`/`pztranu` calls.

# Usage Examples

The `psp_dgemspm_tt` and `psp_zgemspm_tt` subroutines are typically not called directly by end-users. They are invoked by the higher-level dispatcher routines (`psp_dgemspm`, `psp_zgemspm`) within the `pspMspm` module when both `opA` and `opB` arguments indicate a transpose operation ('T' or 'C'). For end-user examples, please refer to the documentation for the `pspMspm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `B`.
*   **`pspUtility`**: Uses `psp_process_opM` (for complex case to distinguish 'T' vs 'C' for `pztranc`/`pztranu`).
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/`.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements, but their routines are not directly called in this file.
*   **`pspMspm_nt`**: This is a key dependency, as the `_tt` routines delegate the core multiplication to `psp_gemspm_nt` after transposing the dense matrix `A`.
*   **MPI Library**: `mpif.h` is included, and MPI calls are made implicitly through ScaLAPACK routines and the called `psp_gemspm_nt` routine.
*   **BLACS Library**: `blacs_gridinfo` and `descinit` are called to set up for ScaLAPACK's parallel transpose.
*   **ScaLAPACK Library**: Uses `pdtran` (for real matrices) or `pztranc`/`pztranu` (for complex matrices) to perform the out-of-place transpose of the dense matrix `A`.
*   The efficiency of this `_tt` implementation is partly tied to the performance of the ScaLAPACK parallel transpose operation and the subsequent `_nt` multiplication. The "TODO" comment in the code suggests that a more direct approach without explicitly forming `op(A)` might be more optimal.
