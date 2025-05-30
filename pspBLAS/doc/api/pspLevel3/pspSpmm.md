# Overview

The `pspSpmm.F90` file defines the `pspSpmm` module, which is a component of the pspBLAS Level 3 library. This module is dedicated to operations involving the multiplication of a sparse matrix (`A`) by a dense matrix (`B`), resulting in a dense matrix (`C`). The fundamental computation performed is `C := alpha*op(A)*op(B) + beta*C`, where `op(A)` and `op(B)` represent the matrices `A` and `B` or their transposes (or conjugate transposes for complex data).

The module provides a public generic interface named `psp_gespmm` (General Sparse Matrix - Dense Matrix product). This interface handles both real and complex data types by dispatching calls to type-specific internal routines (`psp_dgespmm` for real and `psp_zgespmm` for complex). These routines, in turn, further delegate the computational work to specialized subroutines based on the transpose attributes specified for `op(A)` and `op(B)`. These specialized routines (e.g., `psp_gespmm_nn`, `psp_gespmm_nt`) are imported from other dedicated modules within the pspBLAS structure (e.g., `pspSpmm_nn.F90`).

# Key Components

*   **Module `pspSpmm`**:
    *   **Public Interface `psp_gespmm`**: This is the primary user-callable interface for sparse matrix - dense matrix multiplication.
        *   `psp_dgespmm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation. `op(A)` (sparse) is an M-by-K matrix, `op(B)` (dense) is a K-by-N matrix, and the resulting dense matrix `C` is M-by-N.
            *   `A`: `TYPE(psp_matrix_spm)`, the input sparse matrix. Although declared `INTENT(IN)`, commented-out code suggests it might have been intended as `INTENT(INOUT)` for potential internal format conversions.
            *   `opA`: A character ('N', 'T', 'C') specifying the operation on sparse matrix `A`.
            *   `B`: A 2D Fortran array representing the local part of the distributed dense input matrix `B`.
            *   `opB`: A character specifying the operation on dense matrix `B`.
            *   `C`: A 2D Fortran array (inout), representing the local part of the distributed dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Functionality**:
            1.  If `alpha` is zero, the multiplication part `alpha*op(A)*op(B)` is skipped, and only the scaling `C = beta*C` is performed on the output matrix `C`.
            2.  The `opA` and `opB` character flags are processed by `psp_process_opM` (from `pspUtility`) to determine integer transpose codes.
            3.  Based on these transpose codes, one of four specialized routines (`psp_gespmm_nn`, `psp_gespmm_nt`, `psp_gespmm_tn`, or `psp_gespmm_tt`) is called to perform the actual computation. These routines are imported from other `pspSpmm_xx` modules.
            4.  *Note*: Commented-out sections in the code (`!if (A%str_type=='coo') then ...`) suggest an earlier design might have included automatic conversion of the input sparse matrix `A` from COO (Coordinate) format to CSC (Compressed Sparse Column) format before computation, and potentially back afterwards. This functionality is not active in the provided code.

    *   **Subroutine `die(message)`** (private): A standard error handling routine that writes a message to "MatrixSwitch.log" and terminates the program.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros for build configuration and MPI support.
*   `numroc`: Declared as an `EXTERNAL` function, typically from ScaLAPACK, used for calculating local dimensions of distributed arrays.
*   The input sparse matrix `A` is of `TYPE(psp_matrix_spm)`, while dense matrices `B` and `C` are standard Fortran arrays.

# Usage Examples

To perform a sparse matrix - dense matrix multiplication, the user would call the `psp_gespmm` interface. This requires that the pspBLAS environment (including MPI and BLACS if applicable) and the matrices themselves have been properly initialized and distributed.

```fortran
USE pspBLAS  ! Provides access to psp_gespmm via pspLevel3 & pspSpmm
IMPLICIT NONE

TYPE(psp_matrix_spm)                 :: local_sparse_A
REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: local_dense_B, local_dense_C
INTEGER                              :: global_M, global_N, global_K
REAL(dp)                             :: alpha_val, beta_val

! 1. Initialize pspBLAS environment (MPI, BLACS grid via psp_gridinit_2D, etc.).
! 2. Define global matrix dimensions (M, N, K).
! 3. Initialize local_sparse_A (TYPE(psp_matrix_spm) structure and data).
! 4. Determine local dimensions for dense matrices B and C using numroc and
!    pspBLAS distribution schemes. Allocate and initialize local_dense_B and local_dense_C.

alpha_val = 1.0_dp
beta_val = 0.0_dp  ! For C = A*B

! Example: Compute C_dense = alpha * op(A)_sparse * B_dense + beta * C_dense
! Assuming op(A) is A (no transpose) and op(B) is B (no transpose)
CALL psp_gespmm(global_M, global_N, global_K, &
                local_sparse_A, 'N', &  ! op(A) = A
                local_dense_B,  'N', &  ! op(B) = B
                local_dense_C, &
                alpha_val, beta_val)

! local_dense_C now holds the result for the current process's portion of C.

! 5. Finalize pspBLAS environment.
```

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition used for the sparse matrix `A`.
*   **`pspUtility`**: Relies on `psp_process_opM` for parsing transpose flags. The commented-out code also suggests potential use of `psp_coo2csc` and `psp_csc2coo` for format conversions.
*   **`pspMPI`**: Included via `USE`, indicating that the underlying specialized routines (`psp_gespmm_nn`, etc.) are likely MPI-aware and may use information from the `/psp_grid2D/` common block.
*   **`pspLevel1`, `pspLevel2`**: These modules are `USE`d, but their specific functionalities are not directly called within this file. They might be dependencies for the specialized implementation modules.
*   **Specialized Implementation Modules**:
    *   `pspSpmm_nn`: Provides `psp_gespmm_nn` (for A normal, B normal).
    *   `pspSpmm_nt`: Provides `psp_gespmm_nt` (for A normal, B transposed).
    *   `pspSpmm_tn`: Provides `psp_gespmm_tn` (for A transposed, B normal).
    *   `pspSpmm_tt`: Provides `psp_gespmm_tt` (for A transposed, B transposed).
*   **`mpif.h`**: Included if `HAVE_MPI` is defined, providing MPI constants and interfaces for the underlying parallel routines.
*   **`numroc` (External Function)**: Its declaration suggests that calculations of local dimensions for distributed arrays (handled by the specialized routines) are based on ScaLAPACK conventions.

The `pspSpmm` module functions as a high-level dispatcher for sparse-matrix times dense-matrix multiplications. It abstracts the complexities of handling different transpose combinations and data types, delegating these to more specialized modules. The commented-out format conversion logic for the sparse matrix `A` hints that the core computational kernels might have a preference for the CSC format.
