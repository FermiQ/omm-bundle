# Overview

The module `MatrixSwitch_m_register` provides functionality to "register" pre-existing, externally allocated matrix data into a `TYPE(MATRIX)` variable used by the MatrixSwitch library. This allows MatrixSwitch to operate on matrices that were not allocated by its own `m_allocate` routines. It defines interfaces and corresponding subroutines for different matrix storage formats (serial dense, parallel dense block cyclic, parallel sparse COO, parallel sparse CSC) and data types (real and complex). The key idea is to point the internal data pointers of the `TYPE(MATRIX)` variable (e.g., `dval`, `zval`, `iaux1`) to the existing external arrays, rather than allocating new memory and copying data for some fields.

# Key Components

*   **Module `MatrixSwitch_m_register`**: Contains interfaces and subroutines for registering external matrices.
    *   **Interface `m_register_sden`**: For serial dense matrices.
        *   `m_register_sdden(m_name, A)`: Registers a real serial dense matrix. `A` is the 2D real array. Data pointer `m_name%dval` points to `A`.
        *   `m_register_szden(m_name, A)`: Registers a complex serial dense matrix. `A` is the 2D complex array. Data pointer `m_name%zval` points to `A`.
    *   **Interface `m_register_pdbc`** (Requires `HAVE_MPI`): For parallel dense block cyclic matrices.
        *   `m_register_pddbc(m_name, A, desc)`: Registers a real parallel dense block cyclic matrix. `A` is the 2D real local array, `desc` is the BLACS descriptor. `m_name%dval` points to `A`, and `m_name%iaux1` points to `desc`.
        *   `m_register_pzdbc(m_name, A, desc)`: Registers a complex parallel dense block cyclic matrix. `A` is the 2D complex local array, `desc` is the BLACS descriptor. `m_name%zval` points to `A`. `m_name%iaux1` is allocated and `desc` is copied into it.
    *   **Interface `m_register_pcoo`** (Requires `HAVE_PSPBLAS`): For parallel sparse coordinate (COO) matrices.
        *   `m_register_pdcoo(m_name, idx1, idx2, val, desc)`: Registers a real parallel COO matrix. `idx1` (row indices), `idx2` (column indices), `val` (values) are 1D arrays. `desc` is BLACS descriptor. `m_name%iaux1` points to `desc`. Calls `psp_register_spm`.
        *   `m_register_pzcoo(m_name, idx1, idx2, val, desc)`: Registers a complex parallel COO matrix. `m_name%iaux1` points to `desc`. Calls `psp_register_spm`.
    *   **Interface `m_register_pcsc`** (Requires `HAVE_PSPBLAS`): For parallel compressed sparse column (CSC) matrices.
        *   `m_register_pdcsc(m_name, idx1, idx2, val, desc)`: Registers a real parallel CSC matrix. `idx1` (row indices), `idx2` (column pointers), `val` (values) are 1D arrays. `desc` is BLACS descriptor. `m_name%iaux1` points to `desc`. Calls `psp_register_spm`.
        *   `m_register_pzcsc(m_name, idx1, idx2, val, desc)`: Registers a complex parallel CSC matrix. `m_name%iaux1` points to `desc`. Calls `psp_register_spm`.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: Preprocessor macro for configuration.
*   `HAVE_MPI`: Preprocessor macro, enables registration of MPI-based parallel matrices (pdbc).
*   `HAVE_PSPBLAS`: Preprocessor macro, enables registration of pspBLAS-compatible parallel sparse matrices (pcoo, pcsc).
*   `m_name`: The `TYPE(MATRIX)` variable to be initialized by registering the external data.
*   `A` (for dense types): The pre-existing 2D array holding matrix data. Associated using pointer assignment (`=>`).
*   `desc` (for parallel types): The BLACS array descriptor for the distributed matrix. For `pddbc`, `pdcoo`, `pzcoo`, `pdcsc`, `pzcsc` the `m_name%iaux1` points to `desc`. For `pzdbc`, `desc` is copied into an allocated `m_name%iaux1`.
*   `idx1`, `idx2`, `val` (for sparse types): Pre-existing 1D arrays for sparse data (indices, pointers, values), passed to `psp_register_spm`.
*   Matrix type fields set: `dim1`, `dim2`, `is_square`, `str_type`, `is_serial`, `is_real`, `is_sparse`, `dval`, `zval`, `iaux1`, `iaux2` (allocated, stores local dimensions or other info), `spm` (for pspBLAS types), `is_initialized`.
*   `ms_lap_icontxt`, `ms_lap_nprow`, `ms_lap_npcol`: BLACS context variables (used by `psp_register_spm`).

**Usage Examples:**

These subroutines are intended to be called via the public interfaces in the main `MatrixSwitch` module.
```fortran
! Example: Registering an existing serial dense real matrix 'my_array'
TYPE(matrix) :: my_matrix_wrapper
REAL(dp), TARGET, DIMENSION(100,100) :: my_serial_dense_data
! ... populate my_serial_dense_data ...
CALL m_register_sden(my_matrix_wrapper, my_serial_dense_data)

! Example: Registering an existing parallel dense block cyclic real matrix
! Requires my_pdbc_data and my_blacs_desc to be set up appropriately.
#ifdef HAVE_MPI
TYPE(matrix) :: my_pdbc_wrapper
REAL(dp), TARGET, DIMENSION(:,:) :: my_pdbc_data_local_part
INTEGER, TARGET :: my_blacs_descriptor(9)
! ... initialize my_pdbc_data_local_part and my_blacs_descriptor ...
CALL m_register_pdbc(my_pdbc_wrapper, my_pdbc_data_local_part, my_blacs_descriptor)
#endif
```

# Dependencies and Interactions

*   **Uses `MatrixSwitch_ops`**: For `TYPE(matrix)` definition, `dp` kind parameter.
*   **MPI Library** (if `HAVE_MPI` is defined): BLACS (and thus pdbc format) depends on MPI.
*   **pspBLAS Library** (if `HAVE_PSPBLAS` is defined): For registering pcoo and pcsc matrices, `psp_register_spm` is called. This function likely populates the `m_name%spm` structure.
*   **BLACS/ScaLAPACK**: The `desc` array is the BLACS descriptor. Functions `blacs_gridinfo` and `numroc` are used internally when registering pspBLAS types to determine local matrix dimensions.
*   Pointer assignment (`=>`) is used for `dval`, `zval`, and often `iaux1` (the BLACS descriptor) to avoid deep copies. `iaux2` is typically allocated to store local dimensions or other metadata.
*   Sets `is_initialized = .true.` for the `m_name` matrix. The `*_is_allocated` flags for pointer components (`dval`, `zval`, `iaux1` if pointed) should be implicitly false or handled by `MatrixSwitch_ops` type definition to ensure `m_deallocate` only nullifies them. Components allocated within these registration routines (like `iaux2`) will have their corresponding `*_is_allocated` flag set to true.
