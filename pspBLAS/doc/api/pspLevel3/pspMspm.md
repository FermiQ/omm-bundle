# Overview

The `pspMspm.F90` file defines the `pspMspm` module, which is a component of the pspBLAS Level 3 library. This module specializes in "Matrix-Sparse Matrix Product" operations, specifically computing `C := alpha*op(A)*op(B) + beta*C`, where `A` is a dense matrix, `B` is a sparse matrix (of `TYPE(psp_matrix_spm)`), and the resulting matrix `C` is dense.

The module provides a public generic interface `psp_gemspm` that handles both real and complex data types. This interface, in turn, dispatches the call to specific internal routines (`psp_dgemspm` for double precision real and `psp_zgemspm` for double precision complex). These type-specific routines further delegate the computation to specialized worker subroutines based on the transpose operations applied to `A` and `B`. These worker routines (`psp_gemspm_nn`, `psp_gemspm_nt`, `psp_gemspm_tn`, `psp_gemspm_tt`) are imported from other dedicated modules (e.g., `pspMspm_nn.F90`).

# Key Components

*   **Module `pspMspm`**:
    *   **Public Interface `psp_gemspm`**: This is the user-callable interface for dense matrix - sparse matrix multiplication.
        *   `psp_dgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation. `op(A)` is an M-by-K dense matrix, `op(B)` is a K-by-N sparse matrix, and `C` is an M-by-N dense matrix.
            *   `A`: A 2D Fortran array representing the local part of the distributed dense input matrix `A`.
            *   `opA`: A character ('N', 'T', 'C') specifying the operation on matrix `A` (No transpose, Transpose, or Conjugate transpose).
            *   `B`: A variable of `TYPE(psp_matrix_spm)` representing the local part of the distributed sparse input matrix `B`.
            *   `opB`: A character specifying the operation on matrix `B`.
            *   `C`: A 2D Fortran array (inout) representing the local part of the distributed dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Functionality**: These subroutines first process the `opA` and `opB` flags using `psp_process_opM` (from `pspUtility`). Based on the combination of transpose operations, they then select and call the appropriate specialized computation routine (e.g., `psp_gemspm_nn` if no transposes, `psp_gemspm_nt` if `A` is normal and `B` is transposed, etc.). These specialized routines are imported from other modules. If `alpha` is zero, the multiplication is skipped, and only `C = beta*C` is performed.
        *   *Note*: Commented-out code within these routines suggests an earlier consideration to internally convert the sparse matrix `B` to CSC format if it's initially in COO format. This conversion is not active in the current code.

    *   **Subroutine `die(message)`** (private): An error handling routine that writes a message to "MatrixSwitch.log" and terminates the program.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros for build configuration and MPI support.
*   `numroc`: Declared as an `EXTERNAL` function, typically from ScaLAPACK, used for calculating local dimensions of distributed arrays.
*   Input `A` is a dense array, input `B` is `TYPE(psp_matrix_spm)`, and output `C` is a dense array.

# Usage Examples

To perform a dense matrix - sparse matrix multiplication, the user would call `psp_gemspm` after initializing the pspBLAS environment (including MPI and BLACS if applicable) and the necessary matrices.

```fortran
USE pspBLAS  ! Provides access to psp_gemspm via pspLevel3 & pspMspm
IMPLICIT NONE

REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: local_dense_A, local_dense_C
TYPE(psp_matrix_spm)                 :: local_sparse_B
INTEGER                              :: global_M, global_N, global_K
REAL(dp)                             :: alpha_val, beta_val

! 1. Initialize pspBLAS (MPI, BLACS grid via psp_gridinit_2D, etc.).
! 2. Define global matrix dimensions (M, N, K).
! 3. Determine local dimensions for A and C (dense) and B (sparse)
!    using numroc and pspBLAS distribution schemes.
! 4. Allocate and initialize local_dense_A, local_sparse_B, and local_dense_C.
!    (local_sparse_B needs its psp_matrix_spm structure filled appropriately).

alpha_val = 2.0_dp
beta_val = 1.0_dp

! Example: Compute C_dense = alpha * A_dense^T * B_sparse + beta * C_dense
CALL psp_gemspm(global_M, global_N, global_K, &
                local_dense_A, 'T', &    ! op(A) = A^T
                local_sparse_B, 'N', &   ! op(B) = B
                local_dense_C, &
                alpha_val, beta_val)

! local_dense_C now contains the updated result for the current process.

! 5. Finalize pspBLAS environment.
```

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition used for the sparse matrix `B`.
*   **`pspUtility`**: Relies on `psp_process_opM` for parsing the transpose operator flags.
*   **`pspMPI`**: Included via `USE` statement, suggesting that the underlying specialized routines (`psp_gemspm_nn`, etc.) are MPI-aware and likely use information from the `/psp_grid2D/` common block (though not directly accessed in `pspMspm.F90` itself).
*   **`pspLevel1`, `pspLevel2`**: These modules are `USE`d, but their specific functionalities are not directly called within this file. They might be dependencies for the specialized implementation modules.
*   **Specialized Implementation Modules**:
    *   `pspMspm_nn`: Provides `psp_gemspm_nn` (for no transposes).
    *   `pspMspm_nt`: Provides `psp_gemspm_nt` (for `A*op(B)` where `op(B)` is transposed/conjugated).
    *   `pspMspm_tn`: Provides `psp_gemspm_tn` (for `op(A)*B` where `op(A)` is transposed/conjugated).
    *   `pspMspm_tt`: Provides `psp_gemspm_tt` (for `op(A)*op(B)` where both are transposed/conjugated).
*   **`mpif.h`**: Included if `HAVE_MPI` is defined, providing MPI constants and interfaces for the underlying parallel routines.
*   **`numroc` (External Function)**: Its declaration indicates that calculations of local dimensions for distributed arrays (handled by the specialized routines) are based on ScaLAPACK conventions.

The `pspMspm` module acts as a high-level dispatcher for dense-matrix times sparse-matrix multiplications, delegating the complex work of handling different transpose operations and parallel execution (if any) to a set of specific backend modules.
