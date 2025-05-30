# Overview

The `pspMatSum.F90` file defines the `pspMatSum` module, which is a component of the pspBLAS Level 3 library. This module is dedicated to providing routines for performing matrix addition, specifically focusing on operations that involve at least one sparse matrix. The primary operations provided are of the form `C = alpha*A + beta*B`.

The module offers interfaces for two main scenarios:
1.  Summing a sparse matrix (`A`) and a dense matrix (`B`) to produce a dense matrix (`C`).
2.  Summing two sparse matrices (`A` and `B`) to produce a sparse matrix (`C`).

These operations are available for both real and complex data types. The actual computations are delegated to lower-level "Subroutine Subprogram Template" (sst) routines, which are not defined within this file but are expected to be provided elsewhere (likely within the `pspUtility` or a more specialized `psp_spBLAS_Level3` module).

# Key Components

*   **Module `pspMatSum`**:
    *   **Public Interface `psp_sum_spmm`**:
        *   Performs `C = alpha*A + beta*B`, where `A` is sparse, `B` is dense, and `C` is dense.
        *   `psp_dsum_spmm(A, B, C, alpha, beta)`: Real double precision version.
        *   `psp_zsum_spmm(A, B, C, alpha, beta)`: Complex double precision version.
        *   **Arguments**:
            *   `A`: `TYPE(psp_matrix_spm)`, the input sparse matrix.
            *   `B`: `REAL(dp)` or `COMPLEX(dp)` array, the input dense matrix.
            *   `C`: `REAL(dp)` or `COMPLEX(dp)` array (inout), the output dense matrix.
            *   `alpha, beta`: Scalar multipliers.
        *   **Functionality**: These routines check the `str_type` of sparse matrix `A` (either 'coo' for coordinate or 'csc' for compressed sparse column) and then call a corresponding `psp_sst_sum_spmm` routine to perform the computation.

    *   **Public Interface `psp_sum_spmspm`**:
        *   Performs `C = alpha*A + beta*B`, where `A`, `B`, and `C` are all sparse matrices.
        *   `psp_dsum_spmspm(A, B, C, alpha, beta)`: Real double precision version.
        *   `psp_zsum_spmspm(A, B, C, alpha, beta)`: Complex double precision version.
        *   **Arguments**:
            *   `A, B`: `TYPE(psp_matrix_spm)`, the input sparse matrices.
            *   `C`: `TYPE(psp_matrix_spm)` (inout), the output sparse matrix. `C` must be properly initialized (e.g., its `str_type` defined) before calling.
            *   `alpha, beta`: Scalar multipliers.
        *   **Functionality**: These routines use nested `SELECT CASE` statements based on the `str_type` ('coo' or 'csc') of `A`, `B`, and `C` to dispatch the call to the appropriate `psp_sst_sum_spmspm` routine. After the call, `C%nnz` (number of non-zeros in C) is updated.

    *   **Subroutine `die(message)`** (private): An error handling routine that writes to "MatrixSwitch.log" and stops execution.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   Matrix arguments:
    *   For `psp_sum_spmm`: `A` is `TYPE(psp_matrix_spm)`; `B` and `C` are standard Fortran arrays.
    *   For `psp_sum_spmspm`: `A`, `B`, and `C` are all of `TYPE(psp_matrix_spm)`.
*   `A%str_type`, `A%loc_dim1`, `A%loc_dim2`, `A%row_ind`, `A%col_ind`, `A%col_ptr`, `A%dval`, `A%zval`, `A%nnz`: These are fields of the `psp_matrix_spm` derived type, crucial for describing the sparse matrix structure and accessing its data.

# Usage Examples

To use the matrix summation routines, one would typically call the public interfaces `psp_sum_spmm` or `psp_sum_spmspm` after initializing the pspBLAS environment and the matrices.

```fortran
USE pspBLAS  ! This provides access to pspMatSum routines via pspLevel3
IMPLICIT NONE

TYPE(psp_matrix_spm) :: sparse_matrix_A, sparse_matrix_B, sparse_matrix_C
REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: dense_matrix_B, dense_matrix_C
REAL(dp) :: alpha_val, beta_val

! 1. Initialize pspBLAS environment.
! 2. Initialize sparse_matrix_A, sparse_matrix_B (if used), dense_matrix_B (if used).
!    For output matrices (dense_matrix_C, sparse_matrix_C), ensure they are allocated
!    to the correct dimensions. For sparse_matrix_C, its type ('coo' or 'csc')
!    should also be set before calling psp_sum_spmspm.

alpha_val = 1.0_dp
beta_val = 1.0_dp

! Example 1: Sum a sparse matrix and a dense matrix
! dense_C = alpha * sparse_A + beta * dense_B
! CALL psp_sum_spmm(sparse_matrix_A, dense_matrix_B, dense_matrix_C, alpha_val, beta_val)

! Example 2: Sum two sparse matrices
! sparse_C = alpha * sparse_A + beta * sparse_B
! sparse_matrix_C%str_type = 'csc' ! or 'coo', depending on desired output format
! CALL psp_sum_spmspm(sparse_matrix_A, sparse_matrix_B, sparse_matrix_C, alpha_val, beta_val)

! ... use dense_matrix_C or sparse_matrix_C ...

! 3. Finalize pspBLAS environment.
```

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition, which describes sparse matrices (their dimensions, format type, index arrays, and value arrays).
*   **`pspUtility`**: This module is expected to provide the actual implementations of the summation logic through the `psp_sst_sum_spmm` and `psp_sst_sum_spmspm` subroutines. The `pspMatSum` module acts as a dispatcher to these routines.
*   **`pspMPI`**: The `USE pspMPI` statement suggests that the underlying summation operations might be parallelized using MPI, or that they require information from an MPI environment (e.g., process rank for data distribution or output control).
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements, but their functionalities are not directly invoked within the `pspMatSum.F90` code. They might be dependencies for the `_sst_` routines or provide a common context.
*   **`mpif.h`**: Included if `HAVE_MPI` is defined, providing access to MPI constants and interfaces.
*   **`numroc` (External Function)**: Declared as `EXTERNAL`, but not directly used in the routines shown in this file. It might be used by the underlying `_sst_` implementations if they deal with distributed sparse matrices.

The `pspMatSum` module provides a structured interface for matrix addition involving sparse matrices, handling different combinations of sparse and dense inputs and outputs, and delegating the core computational work to specialized template subroutines.
