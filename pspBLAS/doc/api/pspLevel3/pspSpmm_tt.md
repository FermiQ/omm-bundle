# Overview

The `pspSpmm_tt.F90` file provides the `pspSpmm_tt` module, a specialized component of the pspBLAS Level 3 library. This module is responsible for computing the product of a transposed sparse matrix `A` (`op(A)=A^T` or `A^H`) and a transposed dense matrix `B` (`op(B)=B^T` or `B^H`). The operation is formulated as `C := alpha*op(A)*op(B) + beta*C`, where `C` is a dense matrix. This module contains the "_tt" (Transpose-Transpose) specific implementations for the `psp_gespmm` (General Sparse Matrix - Dense Matrix product) functionality.

The current implementation strategy for this double-transpose case involves:
1.  Explicitly transposing the dense matrix `B` (i.e., computing `B^T` or `B^H`) into a temporary distributed dense matrix `tmp`.
2.  Calling the `psp_gespmm_tn` routines (which handle Transpose-Normal products: `op(A)_sparse * B_dense_normal`) using the original `op(A)` for the sparse matrix and the newly created `tmp` (which is `op(B)`) as the "normal" dense matrix.

A "TODO" comment in the source code suggests that this approach (explicitly transposing dense `B`) could be optimized, perhaps by a more direct method or by introducing an on-the-fly transpose capability for the sparse matrix `A` if that were more beneficial.

# Key Components

*   **Module `pspSpmm_tt`**:
    *   **Public Interface `psp_gespmm_tt`**: This interface exposes the routines for transposed-sparse-matrix times transposed-dense-matrix multiplication.
        *   `psp_dgespmm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Sparse matrix `A` (of `TYPE(psp_matrix_spm)`) is globally K-by-M (so `op(A)` is M-by-K). Dense matrix `B` is globally N-by-K (so `op(B)` is K-by-N). The resulting dense matrix `C` is M-by-N.
            *   `A`: `TYPE(psp_matrix_spm)`, the input K-by-M sparse matrix.
            *   `opA`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `B`: A 2D Fortran array for the local part of the distributed N-by-K dense input matrix `B`.
            *   `opB`: Character flag, must be 'T' or 'C' for these routines.
            *   `C`: A 2D Fortran array (inout), the local part of the distributed M-by-N dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm**:
            1.  **Grid Information**: Fetches BLACS 2D process grid details using `blacs_gridinfo` and `psp_icontxt`.
            2.  **Explicit Transpose of Dense Matrix B**:
                *   Determines local dimensions for the input dense matrix `B` (N-by-K) and the temporary transposed dense matrix `tmp` (K-by-N) using `numroc`.
                *   Initializes BLACS array descriptors (`desc_before` for `B`, `desc_after` for `tmp`).
                *   Allocates memory for the temporary distributed matrix `tmp`.
                *   Performs a parallel transpose of `B` into `tmp` using ScaLAPACK routines:
                    *   `pdtran` for real matrix `B`.
                    *   `pztranc` (if `opB='C'`) or `pztranu` (if `opB='T'`) for complex matrix `B`.
                The matrix `tmp` now holds `op(B)`.
            3.  **Call `psp_gespmm_tn`**: Invokes the "Transpose-Normal" routine `psp_gespmm_tn` (from the `pspSpmm_tn` module).
                *   The original sparse matrix `A` and its operator `opA` are passed.
                *   `tmp` (which is `op(B)`) is passed as the non-transposed dense matrix argument (its effective `op` becomes 'N').
                *   The output matrix `C` and scalars `alpha`, `beta` are passed through.
                *   The call effectively becomes: `C = alpha * op(A) * tmp + beta * C`, which is mathematically equivalent to `C = alpha * op(A) * op(B) + beta * C`.
            4.  **Cleanup**: Deallocates the temporary matrix `tmp`.
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes).
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `tmp`: An allocatable 2D array used to store the explicit transpose of the dense matrix `B`.
*   `desc_before`, `desc_after`: BLACS array descriptors used for the parallel transpose operation on `B`.

# Usage Examples

The `psp_dgespmm_tt` and `psp_zgespmm_tt` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmm`, `psp_zgespmm`) located within the `pspSpmm` module. This occurs when both `opA` and `opB` arguments indicate a transpose operation ('T' or 'C'). For end-user examples, please refer to the documentation for the `pspSpmm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `A`.
*   **`pspUtility`**: Uses `psp_process_opM` (for the complex case to determine if `opB` implies conjugate transpose for `pztranc`).
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/`.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but their routines are not directly called in this file.
*   **`pspSpmm_tn`**: This is a key dependency, as the `_tt` routines delegate the core multiplication to `psp_gespmm_tn` after transposing the dense matrix `B`.
*   **MPI Library**: `mpif.h` is included. MPI calls are made implicitly through the ScaLAPACK routines and the called `psp_gespmm_tn` routine.
*   **BLACS Library**: `blacs_gridinfo` and `descinit` are called to set up for ScaLAPACK's parallel transpose operation.
*   **ScaLAPACK Library**: Uses `pdtran` (for real matrices) or `pztranc`/`pztranu` (for complex matrices) to perform the out-of-place transpose of the dense matrix `B`.
*   The efficiency of this `_tt` implementation is partly dependent on the performance of the ScaLAPACK parallel transpose operation on the dense matrix `B` and the subsequent `_tn` multiplication. The "TODO" comment in the code suggests that alternative approaches might be considered for optimization.
