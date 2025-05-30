# Overview

The `psp_spBLAS_Level3.F90` file defines the `psp_spBLAS_Level3` module. This module is a key part of the pspBLAS utility suite, providing fundamental "Subroutine Subprogram Templates" (SST) or core sequential computational kernels for Level 3 sparse BLAS operations. These operations involve matrix-matrix computations where at least one matrix is sparse. The routines in this module operate directly on raw array representations of sparse matrices (e.g., CSC format: column pointers, row indices, and values) and dense matrices (standard Fortran arrays).

The module includes implementations for:
*   Summation of a sparse matrix and a dense matrix.
*   Summation of two sparse matrices (often involving linked lists for accumulation).
*   Sparse matrix - dense matrix multiplication (`C_dense = op(A)_sparse * op(B)_dense`).
*   Dense matrix - sparse matrix multiplication (`C_dense = op(A)_dense * op(B)_sparse`).
*   Sparse matrix - sparse matrix multiplication (`C_sparse = op(A)_sparse * op(B)_sparse`), which itself dispatches to specific routines based on transpose operations and often uses linked lists extensively to construct the result.

# Key Components

*   **Module `psp_spBLAS_Level3`**:
    *   **Public Interfaces and Subroutines**:
        *   `psp_sst_sum_spmm(m,n,alpha,idx1A,idx2A,valA,fmtA,beta,B_dense,C_dense)`:
            *   Computes `C_dense = alpha*A_sparse + beta*B_dense`.
            *   `psp_sst_dsum_spmm` (real), `psp_sst_zsum_spmm` (complex).
        *   `psp_sst_sum_spmspm(m,n,alpha,idx1A,idx2A,valA,fmtA,beta,idx1B,idx2B,valB,fmtB,idx1C,idx2C,valC,fmtC,nnzA,nnzB)`:
            *   Computes `C_sparse = alpha*A_sparse + beta*B_sparse`. The result `valC, idx1C, idx2C` is built using linked lists (via `pspListTool`) to handle dynamic allocation and merging of non-zero elements.
            *   `psp_sst_dsum_spmspm` (real), `psp_sst_zsum_spmspm` (complex).
        *   `psp_sst_gespmm(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A,B_dense,IB,JB,C_dense,IC,JC,beta)`:
            *   Computes `C_dense(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A_sparse)*opB(B_dense(IB:IB+K-1,JB:JB+N-1)) + beta*C_dense(...)`. `A` is sparse CSC.
            *   `psp_sst_dgespmm` (real), `psp_sst_zgespmm` (complex). Implements logic for all four `opA`/`opB` transpose combinations.
        *   `psp_sst_gemspm(M,N,K,opA,opB,alpha,A_dense,IA,JA,row_ind_B,col_ptr_B,val_B,C_dense,IC,JC,beta)`:
            *   Computes `C_dense(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A_dense(IA:IA+M-1,JA:JA+K-1))*opB(B_sparse) + beta*C_dense(...)`. `B` is sparse CSC.
            *   `psp_sst_dgemspm` (real), `psp_sst_zgemspm` (complex). Implements logic for all four `opA`/`opB` transpose combinations. This is the counterpart to `psp_sst_gespmm`.
        *   `psp_sst_gespmspm(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A,row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)`:
            *   General dispatcher for `C_sparse = alpha*opA(A_sparse)*opB(B_sparse) + beta*C_sparse`.
            *   `psp_sst_dgespmspm` (real), `psp_sst_zgespmspm` (complex). These select specific `_nn`, `_nt`, `_tn`, `_tt` routines.
        *   **Specialized `_gespmspm_xx` Kernels (CSC format assumed for inputs)**:
            *   `psp_sst_dgespmspm_nn`, `psp_sst_zgespmspm_nn`: Computes `A*B`. Iterates through columns of B. For each non-zero `B(k,j)`, it scales column `k` of `A` and adds it to a temporary representation of column `j` of `C` (managed via linked lists).
            *   `psp_sst_dgespmspm_nt`, `psp_sst_zgespmspm_nt`: Computes `A*op(B)`. Iterates through columns of `op(B)` (rows of B). For each non-zero in `op(B)(k,j)`, it scales column `k` of `A` and adds to list for column `j` of `C`. Uses an array of lists for `C`.
            *   `psp_sst_dgespmspm_tn`, `psp_sst_zgespmspm_tn`: Computes `op(A)*B`. Iterates columns of B. For each column of `B` and each column of `op(A)` (rows of A), computes dot product and adds to list for `C`.
            *   `psp_sst_dgespmspm_tt`, `psp_sst_zgespmspm_tt`: Computes `op(A)*op(B)`. Similar to `_tn` but `op(B)` is first converted to an array of lists representing its columns.
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   Sparse matrices are primarily handled using their raw components in CSC format (e.g., `row_ind_A`, `col_ptr_A`, `val_A`).
*   `opA`, `opB`: Character flags ('N', 'T', 'C') determining transpose/conjugate operations.
*   Linked list types (`dList`, `zList`, `dNodeData`, `zNodeData`, `dListPtrArray`, `zListPtrArray`) are extensively used in sparse-sparse multiplication routines to dynamically build the result matrix `C`.

# Usage Examples

These are low-level SST routines, designed to be called by higher-level, potentially parallel, pspBLAS routines in the `pspLevel3` directory (e.g., `pspSpmm`, `pspMspm`, `pspSpmSpm` modules and their `_nn`, `_nt`, `_tn`, `_tt` sub-modules). End-users would typically not call these SST routines directly.

```fortran
! Conceptual call from a higher-level pspBLAS routine (e.g., within pspSpmm_nn)
! MODULE pspSpmm_nn_implementation
!   USE psp_spBLAS_Level3 ! To access psp_sst_dgespmm
!   IMPLICIT NONE
! CONTAINS
!   SUBROUTINE some_local_computation_kernel(...)
!     ! ... A_row_idx, A_col_ptr, A_vals represent a local sparse A panel ...
!     ! ... B_dense_panel is a local dense B panel ...
!     ! ... C_dense_panel is the local output dense C panel ...
!     CALL psp_sst_dgespmm(M_panel, N_panel, K_panel, 'N', 'N', alpha, &
!                          A_row_idx, A_col_ptr, A_vals, &
!                          B_dense_panel, 1, 1, & ! Assuming B starts at its (1,1) for this call
!                          C_dense_panel, 1, 1, beta) ! Assuming C starts at its (1,1)
!   END SUBROUTINE
! END MODULE
```

# Dependencies and Interactions

*   **`pspVariable`**: Crucial for linked list type definitions (`dList`, `zList`, `dNodeData`, `zNodeData`, `dListPtrArray`, `zListPtrArray`) and their associated primitive operations (`list_create`, `list_insert`, `list_get_data`, `list_destroy`, etc.) which are heavily used in sparse-sparse products.
*   **`pspListTool`**: Provides higher-level list manipulation routines like `psp_spm2list`, `psp_list2spm`, `psp_list_create_mat`, and `psp_list_combine_listMat`, which are used by `psp_sst_dsum_spmspm` and the various `psp_sst_dgespmspm_xx` routines.
*   **`pspBasicTool`**: For `psp_process_opM` to parse transpose flags, and `psp_sst_fmtCnvt` used in `psp_sst_dsum_spmspm`.
*   **`psp_spBLAS_Level1`**: The `psp_sst_DOT` and `psp_sst_DOTC` routines from Level 1 are used in `psp_sst_dgespmspm_tn` and `psp_sst_dgespmspm_tt` for computing individual elements of the output sparse matrix through dot products of sparse vectors.
*   **`psp_spBLAS_Level2`**: This module is `USE`d, but direct calls to its routines are not immediately apparent in the provided Level 3 SST code. It might be a dependency for other planned routines or for completeness.
*   These SST routines form the computational backbone for sequential operations on sparse data within pspBLAS. The extensive use of linked lists for sparse-sparse products highlights a strategy to manage dynamic fill-in before finalizing the output matrix structure. These kernels are then wrapped by parallel algorithms in other modules.
