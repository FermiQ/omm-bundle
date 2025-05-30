# Overview

The `pspSpmSpm_nt.F90` file provides the `pspSpmSpm_nt` module, a component of the pspBLAS Level 3 library. This module is specifically designed to compute the product of a sparse matrix `A` (in its normal, non-transposed form) and the transpose (or conjugate transpose) of another sparse matrix `B`. The operation performed is `C := alpha*A*op(B) + beta*C`, where `A`, `B`, and `C` are all sparse matrices of `TYPE(psp_matrix_spm)`. This module contains the "_nt" (Normal-Transpose) implementations for the `psp_gespmspm` (General Sparse Matrix - Sparse Matrix product) functionality.

The algorithm employed is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. This involves broadcasting panels of the sparse matrix `B` (which correspond to columns of `op(B)`), performing local sparse-matrix times sparse-matrix-panel multiplications (`A * op(B_panel)`), and then assembling the resulting sparse matrix `C`. A key feature of this implementation is the use of linked lists to gather the distributed segments of the product matrix `C` before converting it to the final sparse matrix format.

# Key Components

*   **Module `pspSpmSpm_nt`**:
    *   **Public Interface `psp_gespmspm_nt`**: This interface exposes the routines for normal-sparse times transposed-sparse matrix multiplication.
        *   `psp_dgespmspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Sparse matrix `A` is M-by-K. Sparse matrix `B` is N-by-K (so `op(B)` is K-by-N). The resulting sparse matrix `C` is M-by-N. All are of `TYPE(psp_matrix_spm)`.
            *   `A, B`: `TYPE(psp_matrix_spm)` (inout). Input sparse matrices. `A` is used as is. `B` is used to form `op(B)`. They are `inout` because the parent `pspSpmSpm` module might have converted them from COO to CSC.
            *   `opA`: Character flag, must be 'N' (or 'n') for these routines.
            *   `opB`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `C`: `TYPE(psp_matrix_spm)` (inout). The output sparse matrix.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid details.
            2.  **Initialization**: Initializes linked list pointers (`list`, `lastElem`) for collecting segments of `C`. `C_loc_total_nnz` is initialized to 0.
            3.  **Outer Loop (`kloop`)**: Iterates over blocks of dimension `N` (columns of `op(B)` and `C`) in chunks of size `psp_update_rank`.
                *   **Panel B Extraction and Broadcast**: The process row (`idx_prow`) owning the current panel of `B` (rows of `B` which correspond to columns of `op(B)`) extracts this panel into temporary CSC-like arrays (`B_loc_idx1`, `B_loc_idx2`, `B_loc_val`) using `psp_copy_spm2st`. These arrays are then broadcast along the process column using `MPI_Bcast` on `psp_mpi_comm_col`.
                *   **Local Sparse-Sparse Multiplication (`psp_sst_gespmspm`)**: Each process computes its local contribution `C_loc_panel = A_local_part * op(B_panel)`. This is done by calling `psp_sst_gespmspm` (an external kernel). The result is stored in a local sparse structure `C_loc` (of `TYPE(psp_MPI_dspm)` or `psp_MPI_zspm)` which holds data in `idx1, idx2, val` arrays.
                *   **Reduce Sparse Panels (`psp_MPI_REDUCE_spm_packed`)**: The computed sparse panels `C_loc` are reduced (summed) across each process row. The result is collected by the process column `idx_pcol` designated to own this panel of the global `C`.
                *   **Assemble C via Linked List**: The process in column `idx_pcol` that received the reduced `C_loc` panel converts this sparse panel data into a linked list segment using `psp_spm2list_shift` (which applies global column offsets `disp`) and appends it to the main linked list (`list`). `C_loc_total_nnz` is updated.
            4.  **Convert List to Sparse Matrix `C_loc_full`**: After the loop, the entire linked list `list` (which now contains all non-zero elements of `A*op(B)` correctly shifted to their global positions) is converted into a single sparse matrix structure (CSC format) using `psp_list2spm`. This result is stored in `C_loc%idx1, C_loc%idx2, C_loc%val` (reusing the `C_loc` variable name, but it now effectively represents the full product `A*op(B)` for the current process's portion of C).
            5.  **Final Scaling and Addition**: The final output matrix `C` is computed by `C = alpha*C_loc_full + beta*C_input` using `psp_sst_sum_spmspm`. This routine sums the newly computed product (now in `C_loc` arrays) with the initial `C` matrix data.
            6.  **Cleanup**: Deallocates temporary arrays for `B_loc_panel` and the linked list structure.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `B_loc_idx1, B_loc_idx2, B_loc_val`: Temporary arrays for holding panels of sparse matrix `B` in CSC format for broadcasting.
*   `C_loc`: A variable of `TYPE(psp_MPI_dspm)` (for real) or `TYPE(psp_MPI_zspm)` (for complex). It's used to store the result of local `A * op(B_panel)` computations and later, after `psp_list2spm`, to hold the full product `A*op(B)` for the current process's part of C.
*   `list, lastElem, listLoc, lastElemLoc`: Pointers of `TYPE(dList)` or `TYPE(zList)` used for managing the linked list that accumulates segments of the product matrix `C`.
*   Key utility/kernel routines:
    *   `psp_copy_spm2st`: Extracts a panel from a `psp_matrix_spm`.
    *   `psp_sst_gespmspm`: Performs local multiplication of `A_local_part * op(B_panel)`.
    *   `psp_MPI_REDUCE_spm_packed`: Custom MPI reduction for packed sparse matrix data.
    *   `psp_spm2list_shift`: Converts a sparse panel to a linked list segment with global index shifting.
    *   `psp_list2spm`: Converts the aggregated linked list back to a CSC sparse matrix representation.
    *   `psp_sst_sum_spmspm`: Performs the final sparse matrix sum `alpha*RESULT + beta*C_initial`.
    *   `list_destroy`: Deallocates a linked list.

# Usage Examples

The `psp_dgespmspm_nt` and `psp_zgespmspm_nt` subroutines are typically called by the higher-level dispatcher routines (`psp_dgespmspm`, `psp_zgespmspm`) in the `pspSpmSpm` module when `opA` is 'N' and `opB` is 'T' or 'C'. For end-user examples, refer to the documentation for the `pspSpmSpm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for `TYPE(psp_matrix_spm)`, `TYPE(psp_MPI_dspm/zspm)`, and `TYPE(dList/zList)` definitions, as well as list manipulation utilities like `list_destroy`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st`, `psp_idx_glb2loc`, `psp_spm2list_shift`, `psp_list2spm`. The core computational kernels `psp_sst_gespmspm` and `psp_sst_sum_spmspm` are also expected from here or `psp_spBLAS_Level3`.
*   **`pspMPI`**: Depends on the initialized `/psp_grid2D/` common block and provides `psp_MPI_REDUCE_spm_packed`.
*   **`pspLevel1`, `pspLevel2`**: Included via `USE`, but not directly invoked.
*   **`pspMatSum`**: While `psp_sum_spmspm` is logically from `pspMatSum`, the call here is to `psp_sst_sum_spmspm`, suggesting it's a lower-level kernel.
*   **MPI Library**: `MPI_Bcast` is used directly. Custom reduction logic is encapsulated in `psp_MPI_REDUCE_spm_packed`.
*   **BLACS Library**: `blacs_gridinfo` is called.
*   This module implements `C_sparse = A_sparse * op(B)_sparse + C_sparse` where `op(B)` is transposed. It uses a complex distributed algorithm involving panel broadcasts of sparse data, local sparse-sparse multiplications, a custom packed sparse matrix reduction, and a linked-list based assembly of the final sparse result. This approach is designed to handle the fill-in that can occur with sparse matrix products.
