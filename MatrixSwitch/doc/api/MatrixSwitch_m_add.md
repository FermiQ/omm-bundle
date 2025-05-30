# Overview

The module `MatrixSwitch_m_add` provides specific implementations for the matrix addition operation `C := alpha*op(A) + beta*C`. It includes subroutines for various combinations of matrix storage formats (dense, compressed sparse column (CSC), compressed sparse row (CSR)), data types (real, complex), and distributions (serial, parallel). These subroutines are called by the generic `m_add` interface in the main `MatrixSwitch` module.

# Key Components

*   **Module `MatrixSwitch_m_add`**: Contains subroutines for matrix addition.
    *   `m_add_pdcscpddbcref(A, C, alpha, beta)`: Real matrix addition where A is parallel CSC (pdcsc), C is parallel dense block cyclic (pddbc). Reference implementation. (Requires `HAVE_PSPBLAS`)
    *   `m_add_pzcscpzdbcref(A, C, alpha, beta)`: Complex matrix addition where A is parallel CSC (pzcsc), C is parallel dense block cyclic (pzdbc). Reference implementation. (Requires `HAVE_PSPBLAS`)
    *   `m_add_sddenref(A, trA, C, alpha, beta)`: Real matrix addition where A and C are serial dense (sdden). Reference implementation. `trA` indicates if A should be transposed.
    *   `m_add_szdenref(A, tcA, C, alpha, beta)`: Complex matrix addition where A and C are serial dense (szden). Reference implementation. `tcA` indicates if A should be transposed or conjugate transposed.
    *   `m_add_sdcscsddenref(A, trA, C, alpha, beta)`: Real matrix addition where A is serial CSC (sdcsc) and C is serial dense (sdden). Reference implementation.
    *   `m_add_sdcsrsddenref(A, trA, C, alpha, beta)`: Real matrix addition where A is serial CSR (sdcsr) and C is serial dense (sdden). Reference implementation.
    *   `m_add_szcscszdenref(A, tcA, C, alpha, beta)`: Complex matrix addition where A is serial CSC (szcsc) and C is serial dense (szden). Reference implementation.
    *   `m_add_szcsrszdenref(A, tcA, C, alpha, beta)`: Complex matrix addition where A is serial CSR (szcsr) and C is serial dense (szden). Reference implementation.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: Preprocessor macro for configuration.
*   `HAVE_PSPBLAS`: Preprocessor macro, enables pspBLAS-specific subroutines.
*   Input parameters in subroutines:
    *   `A`: Input matrix A.
    *   `C`: Input/Output matrix C.
    *   `alpha`: Scalar multiplier for matrix A.
    *   `beta`: Scalar multiplier for matrix C.
    *   `trA` (logical): If true, transpose matrix A (for real versions).
    *   `tcA` (integer): Transposition flag for complex matrix A (0: No transpose, 1: Conjugate transpose, 2: Transpose).
*   Matrix type fields used: `dval` (real data), `zval` (complex data), `spm` (pspBLAS sparse matrix structure), `col_ptr`, `row_ind`, `dim1`, `dim2`, `iaux3`, `iaux4`.

# Usage Examples

These subroutines are typically not called directly but through the `m_add` interface in the `MatrixSwitch` module. An example of how the interface would be used (which in turn calls one of these implementations):
```fortran
! Assuming A_sparse_csc, C_dense are TYPE(matrix) of appropriate kinds
! C_dense = 2.0 * A_sparse_csc + 0.5 * C_dense
CALL m_add(A_sparse_csc, 'N', C_dense, 2.0_dp, 0.5_dp)
```

# Dependencies and Interactions

*   **Uses `MatrixSwitch_ops`**: Likely for the `matrix` type definition, `dp` kind parameter, and potentially other constants or utility functions.
*   **Includes `mpif.h`**: For MPI definitions, used in parallel versions (e.g., `m_add_pdcscpddbcref`).
*   **Preprocessor Directives**:
    *   `#if defined HAVE_CONFIG_H`: Includes `config.h`.
    *   `#ifdef HAVE_PSPBLAS`: Conditionally compiles pspBLAS-dependent routines.
*   Interacts with the `m_add` interface in the `MatrixSwitch` module, which selects the appropriate subroutine based on matrix properties.
