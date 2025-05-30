# Overview

The module `MatrixSwitch_m_copy` provides implementations for the `m_copy` operation, which is responsible for copying one matrix to another. This can involve changing the storage format (e.g., dense to sparse, sparse to dense, or between different sparse formats like COO, CSC, CSR), data type (real/complex, though most direct copies maintain the type), and distribution (serial/parallel). It includes routines for both direct copies and copies with thresholding (zeroing out elements below a certain magnitude). It also provides interfaces and subroutines for copying data from external, pre-existing matrix data (e.g., from a pspBLAS dense block cyclic format to a sparse format managed by `MatrixSwitch`).

# Key Components

*   **Module `MatrixSwitch_m_copy`**: Contains various subroutines for copying matrices.
    *   **Interfaces (for external data copying, require `HAVE_PSPBLAS`)**:
        *   `m_copy_external_pdbcpcoo`: Copies from external parallel dense block cyclic (pdbc) to parallel coordinate (pcoo) format.
            *   `m_copy_external_pddbcpdcoo` (real)
            *   `m_copy_external_pzdbcpzcoo` (complex)
        *   `m_copy_external_pdbcpcsc`: Copies from external parallel dense block cyclic (pdbc) to parallel compressed sparse column (pcsc) format.
            *   `m_copy_external_pddbcpdcsc` (real)
            *   `m_copy_external_pzdbcpzcsc` (complex)
    *   **Internal Copy Subroutines (selection based on input/output formats)**:
        *   `m_copy_sdensdenref(m_name, A)`: Serial dense to serial dense.
        *   `m_copy_sdensdenref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial dense to serial dense with thresholding.
        *   `m_copy_pdbcpdbcref(m_name, A)`: Parallel dense block cyclic to parallel dense block cyclic.
        *   `m_copy_pdbcpdbcref_thre(m_name, A, abs_threshold, soft_threshold)`: Parallel dense block cyclic to parallel dense block cyclic with thresholding.
        *   `m_copy_scooscooref(m_name, A)`: Serial COO to serial COO.
        *   `m_copy_sdenscooref(m_name, A)`: Serial dense to serial COO.
        *   `m_copy_sdenscooref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial dense to serial COO with thresholding.
        *   `m_copy_scscscscref(m_name, A)`: Serial CSC to serial CSC.
        *   `m_copy_scsrscsrref(m_name, A)`: Serial CSR to serial CSR.
        *   `m_copy_sdenscscref(m_name, A)`: Serial dense to serial CSC.
        *   `m_copy_sdenscscref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial dense to serial CSC with thresholding.
        *   `m_copy_sdenscsrref(m_name, A)`: Serial dense to serial CSR.
        *   `m_copy_sdenscsrref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial dense to serial CSR with thresholding.
        *   `m_copy_scoosdenref(m_name, A)`: Serial COO to serial dense.
        *   `m_copy_scoosdenref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial COO to serial dense with thresholding.
        *   `m_copy_scscsdenref(m_name, A)`: Serial CSC to serial dense.
        *   `m_copy_scscsdenref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial CSC to serial dense with thresholding.
        *   `m_copy_scsrsdenref(m_name, A)`: Serial CSR to serial dense.
        *   `m_copy_scsrsdenref_thre(m_name, A, abs_threshold, soft_threshold)`: Serial CSR to serial dense with thresholding.
        *   `m_copy_scscscooref(m_name, A)`: Serial CSC to serial COO.
        *   `m_copy_scsrscooref(m_name, A)`: Serial CSR to serial COO.
        *   `m_copy_scooscscref(m_name, A)`: Serial COO to serial CSC (involves sorting).
        *   `m_copy_scooscsrref(m_name, A)`: Serial COO to serial CSR (involves sorting).

# Important Variables/Constants

*   `HAVE_CONFIG_H`: Preprocessor macro for configuration.
*   `HAVE_PSPBLAS`: Preprocessor macro, enables pspBLAS-specific interfaces and calls to `psp_den2sp_m`.
*   `m_name`: Output matrix (being created/copied into).
*   `A` (in internal copies): Input matrix (source of copy).
*   `A` (in external copies): Input array of matrix values.
*   `desc` (in external copies): BLACS array descriptor for the input parallel matrix.
*   `abs_threshold`: Absolute threshold for zeroing elements during conversions or copies.
*   `soft_threshold`: Value used for soft thresholding (shifting values instead of just zeroing).
*   Matrix type fields used extensively: `is_real`, `dval`, `zval`, `dim1`, `dim2`, `iaux1` (stores BLACS descriptor for pdbc), `iaux2` (stores nnz for sparse or local dims for pdbc), `iaux3`, `iaux4` (sparse matrix index arrays), `spm` (pspBLAS structure).

# Usage Examples

These subroutines are primarily called by the `m_copy` or `m_convert` routines in the main `MatrixSwitch` module.
Direct usage of an external copy routine might look like:
```fortran
! Assuming m_new_coo_matrix is TYPE(matrix) and A_pdbc_values, A_pdbc_desc are existing data
! and HAVE_PSPBLAS is defined.
CALL m_copy_external_pddbcpdcoo(m_new_coo_matrix, A_pdbc_values, A_pdbc_desc, threshold=0.001_dp)
```
An internal copy (usually through `MatrixSwitch.F90`'s `m_copy`):
```fortran
! Convert dense_matrix to sparse_csc_matrix
! CALL m_copy(sparse_csc_matrix, dense_matrix, label='scsc', threshold=1.0E-5_dp)
```

# Dependencies and Interactions

*   **Uses `MatrixSwitch_ops`**: For the `matrix` type definition, `dp` kind parameter, `cmplx_0`, and potentially other constants or error handling routines.
*   **pspBLAS library** (if `HAVE_PSPBLAS` is defined):
    *   `psp_den2sp_m`: Called by the `m_copy_external_pdbcpcsc` and `m_copy_external_pdbcpcoo` routines to convert dense parallel matrices to pspBLAS sparse formats.
*   These routines are the core workers for the `m_copy` and `m_convert` functionalities offered by the main `MatrixSwitch` module. The main module selects the appropriate routine from this module based on the source and target matrix specifications.
*   Memory allocation (`allocate`) for the target matrix (`m_name`) arrays is performed within these subroutines.
