# Overview

The `psp_spBLAS_Level1.F90` file defines the `psp_spBLAS_Level1` module. This module is part of the pspBLAS utility suite and is intended to provide fundamental "Subroutine Subprogram Templates" (SST) for Level 1 sparse BLAS operations. Level 1 operations typically involve vector-vector computations. The routines in this specific file focus on computing dot products of sparse vectors.

A sparse vector is represented here by two arrays: one for the indices of non-zero elements and one for their corresponding values. The indices are assumed to be sorted in ascending order for efficient computation.

# Key Components

*   **Module `psp_spBLAS_Level1`**:
    *   **Public Interface `psp_sst_DOT`**:
        *   Computes the standard dot product of two sparse vectors: `sol = a^T * b`.
        *   `psp_sst_dDOT(N,idxa,vala,idxb,valb,sol,numa,numb)`: Real double precision version.
        *   `psp_sst_zDOT(N,idxa,vala,idxb,valb,sol,numa,numb)`: Complex double precision version.
        *   **Arguments**:
            *   `N`: `INTEGER, INTENT(IN)`, the conceptual total length of the dense vectors from which sparse `a` and `b` are derived.
            *   `idxa(:)`: `INTEGER, INTENT(IN)`, array of sorted 1-based indices of non-zero elements in vector `a`.
            *   `vala(:)`: `REAL(dp)` or `COMPLEX(dp), INTENT(IN)`, array of non-zero values of vector `a`.
            *   `idxb(:)`: `INTEGER, INTENT(IN)`, array of sorted 1-based indices of non-zero elements in vector `b`.
            *   `valb(:)`: `REAL(dp)` or `COMPLEX(dp), INTENT(IN)`, array of non-zero values of vector `b`.
            *   `sol`: `REAL(dp)` or `COMPLEX(dp), INTENT(INOUT)`, the result of the dot product. It is initialized to zero within the routine.
            *   `numa` (Optional): `INTEGER, INTENT(IN)`, number of non-zero elements in vector `a`. Defaults to `SIZE(vala)`.
            *   `numb` (Optional): `INTEGER, INTENT(IN)`, number of non-zero elements in vector `b`. Defaults to `SIZE(valb)`.
        *   **Method**: The routines iterate through both index lists (`idxa`, `idxb`) simultaneously (similar to a merge step in merge-sort). When an index matches in both vectors, the corresponding values from `vala` and `valb` are multiplied and added to `sol`.

    *   **Public Interface `psp_sst_DOTC`**:
        *   Computes the dot product involving the conjugate transpose of the first sparse vector: `sol = a^H * b`. For real data, this is equivalent to `psp_sst_DOT`.
        *   `psp_sst_dDOTC(...)`: Real double precision version (behaves identically to `psp_sst_dDOT`).
        *   `psp_sst_zDOTC(...)`: Complex double precision version. When an index match occurs, `CONJG(vala_element) * valb_element` is added to `sol`.
        *   The arguments are identical to those for `psp_sst_DOT`.

# Important Variables/Constants

*   `dp`: Integer parameter defining the kind for double precision real numbers (`selected_real_kind(15,300)`).
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros. `HAVE_MPI` leads to the inclusion of `mpif.h`, though these specific Level 1 routines do not appear to have MPI-parallel logic themselves.

# Usage Examples

These routines are low-level "Subroutine Subprogram Templates" and are typically called by other, higher-level pspBLAS routines or within more complex algorithms that require sparse vector dot products.

```fortran
MODULE example_dot_product
  USE psp_spBLAS_Level1
  IMPLICIT NONE

  CONTAINS
  SUBROUTINE test_sparse_dot
    REAL(dp) :: result_dot
    COMPLEX(dp) :: result_dotc
    INTEGER, PARAMETER :: N_GLOBAL = 10

    ! Sparse vector X: (1.0 at index 2, 3.0 at index 5, 5.0 at index 8)
    INTEGER, DIMENSION(3) :: x_indices = [2, 5, 8]
    REAL(dp), DIMENSION(3) :: x_values_real  = [1.0_dp, 3.0_dp, 5.0_dp]
    COMPLEX(dp), DIMENSION(3) :: x_values_complex = [(1.0_dp, 0.5_dp), (3.0_dp, -1.0_dp), (5.0_dp, 0.0_dp)]

    ! Sparse vector Y: (2.0 at index 5, 4.0 at index 8, 6.0 at index 9)
    INTEGER, DIMENSION(3) :: y_indices = [5, 8, 9]
    REAL(dp), DIMENSION(3) :: y_values_real  = [2.0_dp, 4.0_dp, 6.0_dp]
    COMPLEX(dp), DIMENSION(3) :: y_values_complex = [(2.0_dp, 1.0_dp), (4.0_dp, -0.5_dp), (6.0_dp, 2.0_dp)]

    ! Real dot product: X^T * Y
    CALL psp_sst_dDOT(N_GLOBAL, x_indices, x_values_real, y_indices, y_values_real, result_dot)
    ! Expected: (3.0 * 2.0) + (5.0 * 4.0) = 6.0 + 20.0 = 26.0
    PRINT *, "Real DOT product:", result_dot

    ! Complex dot product (conjugate X): X^H * Y
    CALL psp_sst_zDOTC(N_GLOBAL, x_indices, x_values_complex, y_indices, y_values_complex, result_dotc)
    ! Expected: CONJG(3.0-1.0i)*(2.0+1.0i) + CONJG(5.0+0.0i)*(4.0-0.5i)
    !         = (3.0+1.0i)*(2.0+1.0i)   + (5.0)*(4.0-0.5i)
    !         = (6+3i+2i-1)             + (20-2.5i)
    !         = (5+5i)                  + (20-2.5i) = 25 + 2.5i
    PRINT *, "Complex DOTC product:", result_dotc

  END SUBROUTINE test_sparse_dot
END MODULE example_dot_product
```

# Dependencies and Interactions

*   **`pspVariable`**: The module `USE`s `pspVariable`. While these specific dot product routines operate on raw arrays, `pspVariable` might define `dp` or other constants if not locally defined, or be a common dependency for other potential Level 1 routines.
*   **`pspBasicTool`**: `USE pspBasicTool`. This could be for utility functions or type definitions used if this module were more extensive.
*   **`pspListTool`**: `USE pspListTool`. Unlikely to be used by these specific dot product routines but included perhaps for other planned Level 1 operations.
*   The efficiency of these dot product routines hinges on the assumption that the input index arrays (`idxa` and `idxb`) are sorted. This allows for a linear time comparison (similar to merging sorted lists).
*   These routines are foundational and would be building blocks for more complex sparse linear algebra algorithms, such as iterative solvers (e.g., Conjugate Gradient) or other matrix operations.
