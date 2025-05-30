# Overview

The `psp_spBLAS.F90` file defines the `psp_spBLAS` module. This module functions as a high-level organizational unit within the pspBLAS library, specifically designed to aggregate and provide access to different levels of core Sparse Basic Linear Algebra Subprograms (spBLAS) routines. These core routines, often referred to as "Subroutine Subprogram Templates" (SST) in other parts of pspBLAS documentation (e.g., `psp_sst_gespmm`), are the fundamental computational kernels that operate directly on the array representations of sparse matrices (e.g., values, indices, pointers).

The `psp_spBLAS` module itself does not implement these kernels but rather uses other modules (`psp_spBLAS_Level1`, `psp_spBLAS_Level2`, `psp_spBLAS_Level3`) where these kernels are expected to be defined.

# Key Components

*   **Module `psp_spBLAS`**:
    *   This module acts as an aggregator. Its primary role is to `USE` and thereby re-export the public entities from the following specialized sparse BLAS level modules:
        *   `psp_spBLAS_Level1`
        *   `psp_spBLAS_Level2`
        *   `psp_spBLAS_Level3`
    *   By doing so, it provides a single point of access to the various levels of core sparse BLAS computational routines.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined during compilation, it includes a `config.h` file. This file is typically used for managing project-wide build configurations and preprocessor definitions.

# Usage Examples

The `psp_spBLAS` module is primarily intended for internal use by other higher-level pspBLAS modules (such as those in the `pspLevel3` directory like `pspSpmm_nn`, `pspSpmSpm_nt`, etc.). These higher-level modules would `USE psp_spBLAS` (often via `USE pspUtility`) to gain access to the core computational kernels needed for their algorithms.

End-users of the pspBLAS library would typically not interact with `psp_spBLAS` directly but would use the more abstract and user-friendly interfaces provided by modules like `pspSpmm`, `pspMspm`, `pspSpmSpm`, or the main `pspBLAS` module.

```fortran
! Conceptual usage within another pspBLAS module, for example,
! in a routine that implements a parallel sparse matrix multiplication.

MODULE pspSomeParallelSparseMultiply
  USE pspUtility ! This would typically make psp_spBLAS contents available too
  ! OR
  ! USE psp_spBLAS ! If direct access to kernels is preferred

  IMPLICIT NONE
  CONTAINS
  SUBROUTINE do_local_sparse_multiplication(...)
    ! ... arguments for local sparse panels A_panel and B_panel ...
    ! ... and output C_panel_accumulator ...

    ! Conceptual call to a Level 3 sparse kernel from psp_spBLAS_Level3
    ! CALL psp_sst_gespmspm( M_panel, N_panel, K_panel, opA, opB, alpha, &
    !                        A_panel_row_ind, A_panel_col_ptr, A_panel_val, &
    !                        B_panel_row_ind, B_panel_col_ptr, B_panel_val, &
    !                        C_panel_accumulator_row_ind, C_panel_accumulator_col_ptr, &
    !                        C_panel_accumulator_val, beta_accumulate)
  END SUBROUTINE do_local_sparse_multiplication
END MODULE pspSomeParallelSparseMultiply
```

# Dependencies and Interactions

The `psp_spBLAS` module has direct dependencies on the following modules, which it aggregates:

*   **`psp_spBLAS_Level1`**: This module is `USE`d and is expected to contain the core implementations or "Subroutine Subprogram Templates" (SST) for Level 1 sparse BLAS operations. These typically involve operations on individual sparse vectors or components of sparse matrices (e.g., manipulating arrays of non-zero values or their indices).
*   **`psp_spBLAS_Level2`**: This module is `USE`d and is expected to house the SST routines for Level 2 sparse BLAS operations. These generally involve operations between a sparse matrix and a dense vector (e.g., sparse matrix-vector multiplication, solving a sparse triangular system for a dense vector).
*   **`psp_spBLAS_Level3`**: This module is `USE`d and is expected to provide the SST routines for Level 3 sparse BLAS operations. These are typically matrix-matrix operations where at least one of the matrices is sparse (e.g., sparse matrix - dense matrix product, or sparse matrix - sparse matrix product kernels like `psp_sst_gespmspm`).

The `psp_spBLAS` module itself is part of the `pspUtility` module's hierarchy (since `pspUtility` uses `psp_spBLAS`). This layered structure helps in organizing the library by separating high-level interfaces from the core computational kernels. The routines aggregated by `psp_spBLAS` are the workhorses that perform the actual numerical computations on sparse data arrays.
