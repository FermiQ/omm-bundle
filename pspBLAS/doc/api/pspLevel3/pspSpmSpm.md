# Overview

The `pspSpmSpm.F90` file defines the `pspSpmSpm` module, which is a component of the pspBLAS Level 3 library. This module is dedicated to performing Sparse Matrix-Sparse Matrix Multiplications (SpMSpM). The core operation is `C := alpha*op(A)*op(B) + beta*C`, where `A`, `B`, and `C` are all sparse matrices represented by the `TYPE(psp_matrix_spm)` derived type from the `pspVariable` module.

The module provides a public generic interface `psp_gespmspm` (General Sparse Matrix-Sparse Matrix product) that handles both real and complex data types. This interface then dispatches the call to specific internal routines (`psp_dgespmspm` for double precision real and `psp_zgespmspm` for double precision complex). These type-specific routines further delegate the computation to specialized worker subroutines based on the transpose attributes of `op(A)` and `op(B)` (e.g., `_nn` for Normal-Normal, `_nt` for Normal-Transpose, etc.). These specialized workers are imported from other modules (e.g., `pspSpmSpm_nn.F90`).

A significant feature of this module is its handling of sparse matrix formats. It checks if the input sparse matrices `A` and `B` (and potentially `C` if it's being overwritten with `beta=0`) are in the Coordinate (COO) format. If so, they are temporarily converted to Compressed Sparse Column (CSC) format using `psp_coo2csc` before the computation. After the computation, if a matrix was converted, it is converted back to its original COO format using `psp_csc2coo`. This implies that the underlying computational kernels for SpMSpM prefer or are optimized for the CSC format.

# Key Components

*   **Module `pspSpmSpm`**:
    *   **Public Interface `psp_gespmspm`**: This is the primary user-callable interface for sparse matrix-sparse matrix multiplication.
        *   `psp_dgespmspm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmspm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation. `op(A)` is an M-by-K sparse matrix, `op(B)` is a K-by-N sparse matrix, and `C` is an M-by-N sparse matrix.
            *   `A, B`: `TYPE(psp_matrix_spm)` (inout), the input sparse matrices. They are `inout` because they might be temporarily converted between COO and CSC formats.
            *   `opA, opB`: Character flags ('N', 'T', 'C') specifying the operation on sparse matrices `A` and `B`.
            *   `C`: `TYPE(psp_matrix_spm)` (inout), the output sparse matrix. It's also `inout` due to potential format conversions and accumulation.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Functionality**:
            1.  **Alpha Check**: If `alpha` is zero, the multiplication `alpha*op(A)*op(B)` is skipped. Only the scaling `C = beta*C` is performed (by scaling the `C%dval` or `C%zval` array directly).
            2.  **Format Conversion (COO to CSC)**:
                *   For matrix `A`: If `A%str_type` is 'coo', `psp_coo2csc(A)` is called. A flag `changeFmtA` is set.
                *   For matrix `B`: If `B%str_type` is 'coo', `psp_coo2csc(B)` is called. A flag `changeFmtB` is set.
                *   For matrix `C`: If `C%str_type` is 'coo' (and `alpha` is non-zero, indicating `C` will be modified by the product), `psp_coo2csc(C)` is called. A flag `changeFmtC` is set.
            3.  **Transpose Processing**: The `opA` and `opB` flags are processed using `psp_process_opM` (from `pspUtility`) to determine integer transpose codes (`trA`, `trB`).
            4.  **Dispatch**: Based on `trA` and `trB`, one of the specialized routines (`psp_gespmspm_nn`, `psp_gespmspm_nt`, `psp_gespmspm_tn`, or `psp_gespmspm_tt`) is called to perform the actual multiplication. These routines are imported from other `pspSpmSpm_xx` modules.
            5.  **Format Back-Conversion (CSC to COO)**:
                *   If `changeFmtA` is true, `psp_csc2coo(A)` is called.
                *   If `changeFmtB` is true, `psp_csc2coo(B)` is called.
                *   If `changeFmtC` is true, `psp_csc2coo(C)` is called.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   `numroc`: Declared as an `EXTERNAL` function (typically from ScaLAPACK).
*   `changeFmtA, changeFmtB, changeFmtC`: Logical flags used internally to track whether a sparse matrix format conversion (COO to CSC and back) was performed.
*   `A%str_type`, `B%str_type`, `C%str_type`: Fields within the `psp_matrix_spm` type that store the sparse matrix format identifier ('coo', 'csc', etc.).
*   `psp_coo2csc`, `psp_csc2coo`: Utility functions (from `pspUtility`) for converting between COO and CSC sparse matrix formats.

# Usage Examples

To perform a sparse matrix-sparse matrix multiplication, the user would call the `psp_gespmspm` interface. This requires prior initialization of the pspBLAS environment and the sparse matrices.

```fortran
USE pspBLAS  ! This provides access to psp_gespmspm via pspLevel3 & pspSpmSpm
IMPLICIT NONE

TYPE(psp_matrix_spm) :: sparse_A, sparse_B, sparse_C
INTEGER              :: M_global, N_global, K_global
REAL(dp)             :: alpha_val, beta_val

! 1. Initialize pspBLAS environment (MPI, BLACS grid if applicable).
! 2. Define global matrix dimensions (M_global, N_global, K_global).
! 3. Initialize sparse_A, sparse_B, and sparse_C. This includes setting their
!    dimensions, format type (e.g., 'coo' or 'csc'), and populating their
!    index and value arrays for their respective local portions.
!    Ensure sparse_C is initialized correctly, especially its str_type,
!    as it might be converted if it's 'coo'.

alpha_val = 1.0_dp
beta_val = 0.0_dp ! For C = A*B

! Example: Compute C_sparse = alpha * A_sparse * B_sparse^T + beta * C_sparse
! Assuming A, B, C are already appropriately distributed if in a parallel environment.
CALL psp_gespmspm(M_global, N_global, K_global, &
                  sparse_A, 'N', &      ! op(A) = A
                  sparse_B, 'T', &      ! op(B) = B^T
                  sparse_C, &
                  alpha_val, beta_val)

! sparse_C now holds the result of the multiplication. If any input/output
! matrices were 'coo', they would have been converted to 'csc' for computation
! and then converted back.

! 4. Finalize pspBLAS environment.
```

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition used for all sparse matrices `A, B, C`.
*   **`pspUtility`**: Relies on `psp_process_opM` for parsing transpose flags, and crucially on `psp_coo2csc` and `psp_csc2coo` for sparse matrix format conversions.
*   **`pspMPI`**: Included via `USE`, suggesting that the underlying specialized SpMSpM routines (`_nn`, `_nt`, etc.) are likely MPI-aware and may use information from the `/psp_grid2D/` common block for distributed computations.
*   **`pspLevel1`, `pspLevel2`, `pspMatSum`**: These modules are `USE`d, but their specific functionalities are not directly invoked within this `pspSpmSpm.F90` file. They might be dependencies for the specialized `pspSpmSpm_xx` modules or provide a common operational context.
*   **Specialized Implementation Modules**:
    *   `pspSpmSpm_nn`: Provides `psp_gespmspm_nn`.
    *   `pspSpmSpm_nt`: Provides `psp_gespmspm_nt`.
    *   `pspSpmSpm_tn`: Provides `psp_gespmspm_tn`.
    *   `pspSpmSpm_tt`: Provides `psp_gespmspm_tt`.
*   **`mpif.h`**: Included if `HAVE_MPI` is defined, providing MPI constants and interfaces for the underlying parallel routines.
*   **`numroc` (External Function)**: Its declaration suggests that dimension calculations for distributed data (handled by the specialized routines) might be based on ScaLAPACK conventions.

The `pspSpmSpm` module acts as a high-level dispatcher for sparse matrix-sparse matrix multiplications. It standardizes the input matrix formats to CSC for the computational kernels and then delegates the actual work to specific modules based on the transpose requirements of the operation.
