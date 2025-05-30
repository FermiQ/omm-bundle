# Overview

The `pspListTool.F90` file defines the `pspListTool` module, which provides a collection of utility subroutines for working with linked lists that represent sparse matrix data. In sparse matrix computations, especially operations like matrix multiplication where the resulting non-zero pattern can be complex and not known in advance, linked lists offer a flexible way to store and accumulate non-zero elements (typically as `(row_index, column_index, value)` triplets) before converting them into a final, more structured sparse format like Compressed Sparse Column (CSC) or Coordinate (COO).

This module facilitates:
*   Conversion of sparse matrix data (from raw array representations) into linked lists, with options to apply index offsets (useful for assembling a larger matrix from sub-blocks).
*   Conversion of these linked lists back into standard sparse matrix array formats.
*   Creation of a linked list directly from sparse matrix arrays.
*   Merging or combining two linked lists, or a linked list with sparse matrix data from arrays, effectively performing an element-wise addition where indices overlap.
*   Basic list utilities such as finding the last element of a list or printing its contents for debugging purposes.

The routines are generally provided for both real (`dList`) and complex (`zList`) data types, leveraging underlying linked list type definitions and primitive operations from the `pspVariable` module.

# Key Components

*   **Module `pspListTool`**:
    *   **Interfaces & Public Subroutines**:
        *   `psp_spm2list(m,n,idx1,idx2,val,nnz,fmt,list)`: Converts sparse matrix data given in raw arrays (`idx1`, `idx2`, `val` in format `fmt`) into a single linked list.
            *   `psp_dspm2list` (real), `psp_zspm2list` (complex).
        *   `psp_spm2list_shift(m,n,idx1,idx2,val,nnz,fmt,list,last,disp)`: Similar to `psp_spm2list`, but applies row and column index shifts based on the optional `disp(:)` argument before inserting into the list. Updates `last` to point to the last element added.
            *   `psp_dspm2list_shift` (real), `psp_zspm2list_shift` (complex).
        *   `psp_spm2lists_shift(m,n,idx1,idx2,val,nnz,fmt,listArraySt,listArrayEd,listCreated,disp)`: Converts sparse matrix data into an array of linked lists (`listArraySt`, `listArrayEd`). Each list in the array typically corresponds to a column (or block of columns) of the sparse matrix. Applies shifts via `disp`. `listCreated` tracks which lists in the array have been initiated.
            *   `psp_dspm2lists_shift` (real), `psp_zspm2lists_shift` (complex).
        *   `psp_list2spm(m,n,idx1,idx2,val,fmt,list,nnz,ifCount)`: Converts a linked list `list` back into raw sparse matrix arrays (`idx1`, `idx2`, `val`). The output format is specified by `fmt`. If `fmt` is 'csc', it internally converts from COO (the list natural format) to CSC using `psp_sst_coo2csc`. `nnz` can be pre-calculated or counted using `list_count`.
            *   `psp_dlist2spm` (real), `psp_zlist2spm` (complex).
        *   `psp_list_create_mat(listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)`: Creates a new linked list `listA` by adding elements from sparse matrix arrays (`idxi`, `idxj`, `val`), scaling values by `beta`. `numAddInB` specifies how many elements to add.
            *   `psp_dlist_create_mat` (real), `psp_zlist_create_mat` (complex).
        *   `psp_list_combine_listMat(m,n,alpha,listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)`: Combines (merges/adds) elements from sparse matrix arrays (`idxi`, `idxj`, `beta*val`) into an existing linked list `listA` (which is first scaled by `alpha`). Assumes elements are sorted by position for merging.
            *   `psp_dlist_combine_listMat` (real), `psp_zlist_combine_listMat` (complex).
        *   `psp_list_combine_listList(m,n,alpha,listA,lenA,beta,listB,lenB,lastElem)`: Combines two linked lists, `listA = alpha*listA + beta*listB`. Elements are merged based on their row/column indices.
            *   `psp_dlist_combine_listList` (real), `psp_zlist_combine_listList` (complex).
        *   `psp_list_getLast(list,lastElem)`: Traverses `list` to find and return a pointer to its last element in `lastElem`.
            *   `psp_dlist_getLast` (real), `psp_zlist_getLast` (complex).
        *   `psp_list_print(text, list, len)`: Prints the elements (row, col, value) of a linked list for debugging, prefixed by `text`. Optionally prints list length `len`.
            *   `psp_dlist_print` (real), `psp_zlist_print` (complex).
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   **Linked List Types**: The module extensively uses pointer types for linked lists like `dList`, `zList` and node data types `dNodeData`, `zNodeData`. It also uses `dListPtrArray` and `zListPtrArray` (arrays of pointers to list elements), which are expected to be defined in the `pspVariable` module.
*   **Basic List Primitives**: The routines rely on fundamental linked list operations such as `list_create`, `list_insert`, `list_insert_head`, `list_next`, `list_get_data`, `list_count`, and `list_destroy`. These are assumed to be provided by the `pspVariable` module.

# Usage Examples

These list manipulation tools are primarily intended for internal use within pspBLAS, particularly in Level 3 sparse matrix operations like SpMSpM (Sparse Matrix-Sparse Matrix Multiplication), where the resulting sparsity pattern is generated dynamically.

```fortran
! Conceptual example of building a sparse matrix C using list tools,
! perhaps as part of a sparse matrix multiplication C = A*B.
USE pspListTool
USE pspVariable
IMPLICIT NONE

TYPE(dList), POINTER :: C_list_head, C_list_tail, panel_list_head, panel_list_tail
TYPE(psp_matrix_spm) :: final_C_sparse
INTEGER :: m_dim, n_dim, nnz_panel, total_nnz_C
INTEGER, DIMENSION(2) :: panel_offset
REAL(dp), ALLOCATABLE :: panel_val(:)
INTEGER, ALLOCATABLE :: panel_row_idx(:), panel_col_idx(:)

! Initialize list for C
C_list_head => NULL(); C_list_tail => NULL(); total_nnz_C = 0

! Loop over panels or blocks contributing to C
  ! ... compute a partial product panel, resulting in coo arrays:
  !     panel_row_idx, panel_col_idx, panel_val, with nnz_panel elements ...
  ! ... set panel_offset for global positioning ...

  CALL psp_dspm2list_shift(m_panel_dim, n_panel_dim, panel_row_idx, panel_col_idx, &
                           panel_val, nnz_panel, 'coo', panel_list_head, panel_list_tail, panel_offset)

  ! Combine panel_list into the main C_list
  IF (.NOT. ASSOCIATED(C_list_head)) THEN
    C_list_head => panel_list_head
    C_list_tail => panel_list_tail
    total_nnz_C = nnz_panel
  ELSE
    CALL psp_dlist_combine_listList(m_dim, n_dim, 1.0_dp, C_list_head, total_nnz_C, &
                                    1.0_dp, panel_list_head, nnz_panel, C_list_tail)
    CALL list_destroy(panel_list_head) ! Assuming panel_list_head is a distinct list
  END IF
! End loop

! Convert the final assembled list to a CSC sparse matrix
CALL psp_dlist2spm(m_dim, n_dim, final_C_sparse%row_ind, final_C_sparse%col_ptr, &
                   final_C_sparse%dval, 'csc', C_list_head, final_C_sparse%nnz, .TRUE.)
CALL list_destroy(C_list_head)
```

# Dependencies and Interactions

*   **`pspVariable`**: This is a critical dependency. `pspListTool` relies on `pspVariable` for the definitions of linked list types (`dList`, `zList`, `dNodeData`, `zNodeData`, `dListPtrArray`, `zListPtrArray`) and the fundamental procedures to manipulate these lists (e.g., `list_create`, `list_insert`, `list_get_data`, `list_next`, `list_count`, `list_destroy`).
*   **`pspBasicTool`**: Uses `psp_sst_coo2csc` from `pspBasicTool` when converting a linked list (which inherently stores COO-like data) to the CSC sparse matrix format in the `psp_list2spm` routines.
*   **MPI**: While `mpif.h` is included if `HAVE_MPI` is defined, the routines within `pspListTool.F90` itself do not appear to make direct MPI calls. Parallelism involving these list tools would typically be managed by the calling routines, for instance, by having each MPI process build its own local list of sparse matrix elements, which might then be gathered or processed in a distributed manner.
*   The routines in this module are essential for algorithms that construct sparse matrices dynamically, especially when the number and location of non-zero elements are not known beforehand. They are likely heavily utilized by the sparse matrix multiplication (SpMSpM) routines in `pspLevel3`.
