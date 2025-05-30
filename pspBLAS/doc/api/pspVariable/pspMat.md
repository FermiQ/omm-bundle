# Overview

The `pspMat.F90` file defines the `pspMat` module, which is a crucial part of the `pspVariable` collection within the pspBLAS library. This module's primary contribution is the definition of the `psp_matrix_spm` derived data type. This type is the fundamental structure used throughout pspBLAS to represent sparse matrices, including their data, format, and distribution information for parallel processing.

In addition to the type definition, the `pspMat` module provides essential routines for:
1.  **Registering Data**: `psp_register_spm` allows existing raw sparse matrix data (indices, values, and descriptor information) to be populated into a `psp_matrix_spm` variable.
2.  **Deallocation**: `psp_deallocate_spm` provides a way to free the allocatable components of a `psp_matrix_spm` variable.

The `psp_matrix_spm` type is designed to support distributed sparse matrices and can accommodate different storage formats, primarily COO (Coordinate) and CSC (Compressed Sparse Column).

# Key Components

*   **Module `pspMat`**:
    *   **`TYPE psp_matrix_spm`**: This is the central derived type for representing sparse matrices in pspBLAS. Its components are:
        *   `str_type`: `CHARACTER(3)`, a label identifying the storage format (e.g., 'csc', 'coo'). The comments also mention 'dym' (dynamic COO using linked lists), but the direct components for this (like list pointers) are commented out in the type definition itself.
        *   `is_initialized`: `LOGICAL`, flag indicating if the matrix structure has been properly set up (defaults to `.FALSE.`).
        *   `is_serial`: `LOGICAL`, flag indicating if the matrix is serial or distributed. The `psp_register_spm` routines in this module set it to `.FALSE.`, implying a focus on distributed matrices.
        *   `is_real`: `LOGICAL`, flag indicating if the matrix data is real (`.TRUE.`) or complex (`.FALSE.`).
        *   `is_square`: `LOGICAL`, flag indicating if the global dimensions of the matrix are square (`glb_dim1 == glb_dim2`).
        *   `nnz`: `INTEGER`, stores the number of non-zero entries held by the local process.
        *   `glb_dim1`, `glb_dim2`: `INTEGER`, global row and column dimensions of the sparse matrix.
        *   `loc_dim1`, `loc_dim2`: `INTEGER`, local row and column dimensions of the matrix portion held by the current process.
        *   `row_ind(:)`: `INTEGER, ALLOCATABLE`, array for row indices of non-zero elements.
        *   `col_ind(:)`: `INTEGER, ALLOCATABLE`, array used for column indices if `str_type` is 'coo'.
        *   `col_ptr(:)`: `INTEGER, ALLOCATABLE`, array used for column start/end pointers if `str_type` is 'csc'.
        *   `desc(9)`: `INTEGER` array, the BLACS (Basic Linear Algebra Communication Subprograms) array descriptor for the distributed sparse matrix.
        *   `dval(:)`: `REAL(dp), ALLOCATABLE`, array for storing real double precision non-zero values.
        *   `zval(:)`: `COMPLEX(dp), ALLOCATABLE`, array for storing complex double precision non-zero values.

    *   **Public Interface `psp_register_spm`**: This interface provides routines to populate a `psp_matrix_spm` variable from existing raw sparse data arrays.
        *   `psp_register_dspm(m_name, idx1, idx2, val, desc, fmt, loc_dims, nprow, npcol)`: For real sparse matrices.
        *   `psp_register_zspm(m_name, idx1, idx2, val, desc, fmt, loc_dims, nprow, npcol)`: For complex sparse matrices.
        *   **Arguments**:
            *   `m_name`: The `TYPE(psp_matrix_spm)` variable to be populated.
            *   `idx1(:)`: Array of row indices.
            *   `idx2(:)`: Array of column indices (for 'coo') or column pointers (for 'csc').
            *   `val(:)`: Array of non-zero values.
            *   `desc(9)`: BLACS array descriptor.
            *   `fmt`: Character string ('coo' or 'csc') specifying the format of input `idx1`, `idx2`, `val`.
            *   `loc_dims(2)`: Local dimensions `[local_rows, local_cols]` on the current process.
            *   `nprow, npcol`: Number of processor rows and columns in the grid (used for context, though `desc` should contain this).
        *   **Functionality**: These routines copy the provided descriptor and dimensions into `m_name`. They allocate and copy the `val`, `row_ind`, and either `col_ind` (for 'coo') or `col_ptr` (for 'csc') arrays.

    *   **Public Interface `psp_deallocate_spm`**:
        *   `psp_deallocate_spm(m_name)`: This subroutine deallocates all allocatable array components within the `m_name` (`psp_matrix_spm`) variable (i.e., `row_ind`, `col_ind`, `col_ptr`, `dval`, `zval`) and sets `is_initialized` to `.FALSE.`.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`: Integer parameter defining the kind for double precision real numbers.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   `numroc`: Declared as an `EXTERNAL` function (typically from ScaLAPACK), though not directly used in `psp_register_spm` which takes `loc_dims` as input. It's used by `psp_spm_zeros` (in `pspBasicTool`), which might be called before registering some types of data.
*   The `psp_matrix_spm` type itself is the most crucial definition in this module, serving as the standard representation for sparse matrices across the pspBLAS library.

# Usage Examples

```fortran
MODULE example_psp_matrix_usage
  USE pspMat         ! For psp_matrix_spm, psp_register_spm, psp_deallocate_spm
  USE pspMPI         ! For psp_gridinit_2D and common block variables like psp_icontxt
                     ! (assuming psp_bs_def_row, psp_bs_def_col are also set there)
  IMPLICIT NONE

  SUBROUTINE test_matrix_registration
    TYPE(psp_matrix_spm) :: A_sparse
    INTEGER :: global_rows = 100, global_cols = 100
    INTEGER :: local_rows, local_cols, nnz_local
    INTEGER, ALLOCATABLE :: r_indices(:), c_pointers(:)
    REAL(dp), ALLOCATABLE :: values(:)
    INTEGER :: blacs_descriptor(9), blacs_context_handle
    INTEGER :: num_procs_rows, num_procs_cols, my_proc_row, my_proc_col
    INTEGER :: ierr

    ! 1. Initialize MPI and BLACS grid (e.g., using psp_gridinit_2D)
    !    This sets psp_icontxt, psp_nprow, psp_npcol from the /psp_grid2D/ common block.
    !    For this example, assume these are already set and retrieve them:
    blacs_context_handle = psp_icontxt
    CALL BLACS_GRIDINFO(blacs_context_handle, num_procs_rows, num_procs_cols, my_proc_row, my_proc_col)

    ! 2. Determine local dimensions and create sample local sparse data (CSC format)
    local_rows = NUMROC(global_rows, psp_bs_def_row, my_proc_row, 0, num_procs_rows)
    local_cols = NUMROC(global_cols, psp_bs_def_col, my_proc_col, 0, num_procs_cols)

    ! (Simplified example: create a diagonal local matrix if my_proc_row == my_proc_col)
    IF (my_proc_row == my_proc_col) THEN
      nnz_local = MIN(local_rows, local_cols)
      ALLOCATE(r_indices(nnz_local), c_pointers(local_cols + 1), values(nnz_local))
      c_pointers = 1
      DO i = 1, nnz_local
        r_indices(i) = i ! Assuming 1-based local indexing for this example part
        values(i) = REAL(i, dp)
        c_pointers(i+1) = c_pointers(i) + 1
      END DO
      DO i = nnz_local + 1, local_cols
         c_pointers(i+1) = c_pointers(i)
      ENDDO
    ELSE
      nnz_local = 0
      ALLOCATE(r_indices(0), c_pointers(local_cols + 1), values(0))
      c_pointers = 1
    END IF

    ! 3. Create BLACS descriptor (details depend on actual distribution)
    CALL DESCINIT(blacs_descriptor, global_rows, global_cols, &
                  psp_bs_def_row, psp_bs_def_col, 0, 0, &
                  blacs_context_handle, local_rows, ierr)

    ! 4. Register the local sparse data into the psp_matrix_spm structure
    CALL psp_register_spm(A_sparse, r_indices, c_pointers, values, blacs_descriptor, &
                          'csc', [local_rows, local_cols], num_procs_rows, num_procs_cols)

    PRINT *, "Sparse matrix registered. Local NNZ:", A_sparse%nnz

    ! 5. Use A_sparse with other pspBLAS routines...

    ! 6. Deallocate when done
    CALL psp_deallocate_spm(A_sparse)
    DEALLOCATE(r_indices, c_pointers, values)

  END SUBROUTINE test_matrix_registration
END MODULE example_psp_matrix_usage
```

# Dependencies and Interactions

*   **`pspListType`**: The `USE pspListType` statement is present, likely because the `pspVariable` module (which `pspMat` is conceptually part of) aims to provide all variable and type definitions, including list types, even if `pspMat.F90` itself doesn't directly use the list interfaces. The commented-out `dymdcoo` and `dymzcoo` components in `psp_matrix_spm` might have intended to use these list types.
*   **BLACS/ScaLAPACK**: The `psp_matrix_spm` type includes a BLACS descriptor (`desc`), and the `psp_register_spm` routines expect this to be correctly set up for distributed matrices. The `numroc` function (external, from ScaLAPACK) is relevant for calculating local dimensions based on global dimensions and distribution parameters, although `psp_register_spm` directly takes `loc_dims`.
*   This module is foundational for sparse matrix operations in pspBLAS. Nearly all other modules that deal with sparse matrices will use the `psp_matrix_spm` type defined here. Modules like `pspBasicTool` and those in `pspLevel1`, `pspLevel2`, and `pspLevel3` will operate on or produce instances of `psp_matrix_spm`.
*   The `psp_register_spm` routines copy the data from the input arrays into the allocated components of the `psp_matrix_spm` structure. It does not take ownership of the input arrays via pointer assignment.
