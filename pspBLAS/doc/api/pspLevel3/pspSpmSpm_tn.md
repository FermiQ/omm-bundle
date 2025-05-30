# Overview

The `pspSpmSpm_tn.F90` file provides the `pspSpmSpm_tn` module, a specialized component of the pspBLAS Level 3 library. This module is responsible for computing the product of a transposed (or conjugate-transposed) sparse matrix `A` and a sparse matrix `B` in its normal (non-transposed) form. The operation is `C := alpha*op(A)*B + beta*C`, where `A, B, C` are all sparse matrices of `TYPE(psp_matrix_spm)`. This module contains the "_tn" (Transpose-Normal) specific implementations for the `psp_gespmspm` (General Sparse Matrix - Sparse Matrix product) functionality.

The algorithm is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. It involves broadcasting panels of the sparse matrix `A` (which correspond to rows of `op(A)`), performing local multiplications of these transposed panels with the appropriate parts of sparse matrix `B` (`op(A_panel) * B_local_part`), and then assembling the resulting sparse matrix `C`. A key feature of this implementation is the use of an array of linked lists, where each list corresponds to a block of columns in the output matrix `C`, to manage the accumulation of results before final conversion to the `psp_matrix_spm` format.

# Key Components

*   **Module `pspSpmSpm_tn`**:
    *   **Public Interface `psp_gespmspm_tn`**: This interface exposes the routines for transposed-sparse times normal-sparse matrix multiplication.
        *   `psp_dgespmspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Sparse matrix `A` is K-by-M (so `op(A)` is M-by-K). Sparse matrix `B` is K-by-N. The resulting sparse matrix `C` is M-by-N. All are of `TYPE(psp_matrix_spm)`.
            *   `A, B`: `TYPE(psp_matrix_spm)` (inout). Input sparse matrices. `A` is used to form `op(A)`. They are `inout` because the parent `pspSpmSpm` module might have converted them from COO to CSC.
            *   `opA`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `opB`: Character flag, must be 'N' (or 'n') for these routines.
            *   `C`: `TYPE(psp_matrix_spm)` (inout). The output sparse matrix.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid details.
            2.  **Initialization**: Initializes an array of linked list structures (`listArraySt`, `listArrayEd`) and corresponding logical flags (`listCreated`) for each block of columns of `C` that the current process will handle. `C_loc_total_nnz` is initialized to 0.
            3.  **Outer Loop (`kloop`)**: Iterates over blocks of dimension `M` (rows of `op(A)` and `C`) in chunks of size `psp_update_rank`.
                *   **Panel A Extraction and Broadcast**: The process column (`idx_pcol`) that owns the current panel of `A` (columns of `A` which correspond to rows of `op(A)`) extracts this panel into temporary CSC-like arrays (`A_loc_idx1`, `A_loc_idx2`, `A_loc_val`) using `psp_copy_spm2st`. These arrays are then broadcast row-wise across the process grid using `MPI_Bcast` on `psp_mpi_comm_row`.
                *   **Local Sparse-Sparse Multiplication (`psp_sst_gespmspm`)**: Each process computes its local contribution `C_loc_panel = op(A_panel) * B_local_part`. This is done by calling `psp_sst_gespmspm` (an external kernel). The result is stored in a local sparse structure `C_loc` (of `TYPE(psp_MPI_dspm)` or `psp_MPI_zspm)` which holds data in `idx1, idx2, val` arrays. The `B` matrix components (`B%row_ind`, `B%col_ptr`, `B%dval`/`%zval`) are passed directly to `psp_sst_gespmspm`, which implies it handles accessing the relevant local portion of `B`.
                *   **Reduce Sparse Panels (`psp_MPI_REDUCE_spm_packed`)**: The computed sparse panels `C_loc` are reduced (summed) across each process column. The result is collected by the process row `idx_prow` designated to own this panel of the global `C`.
                *   **Assemble C via Array of Linked Lists**: The process in row `idx_prow` that received the reduced `C_loc` panel (which is a horizontal strip of `C`) further processes it. It uses `psp_spm2lists_shift` to convert columns of this `C_loc` panel into segments and appends them to the appropriate linked list in the `listArraySt`/`listArrayEd` arrays (each list corresponding to a global column or block of columns of `C`). The `disp` argument in `psp_spm2lists_shift` handles global row offsets.
            4.  **Combine Columnar Linked Lists & Convert to Global `C_loc`**: After the loop, the array of linked lists (`listArraySt`) is consolidated into a single master linked list (`list`). This master list, containing all computed non-zero elements of `op(A)*B` for the current process's portion of `C`, is then converted into a CSC sparse matrix representation using `psp_list2spm`. This result is stored in `C_loc%idx1, C_loc%idx2, C_loc%val` (reusing `C_loc` variable components).
            5.  **Final Scaling and Addition**: The final output matrix `C` is computed by `C = alpha*C_loc_result + beta*C_input` using `psp_sst_sum_spmspm`.
            6.  **Cleanup**: Deallocates temporary arrays and the linked list structures.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `A_loc_idx1, A_loc_idx2, A_loc_val`: Temporary arrays for holding panels of sparse matrix `A` in CSC format for broadcasting.
*   `C_loc`: A variable of `TYPE(psp_MPI_dspm)` (for real) or `TYPE(psp_MPI_zspm)` (for complex). It's used to store the result of local `op(A_panel) * B_local_part` computations and later, after `psp_list2spm`, to hold the full product `op(A)*B` for the current process's part of `C`.
*   `listArraySt, listArrayEd, listArrayStLoc, listArrayEdLoc`: Arrays of `TYPE(dListPtrArray)` or `TYPE(zListPtrArray)` used to manage multiple linked lists, each corresponding to a local column or block of columns of the result matrix `C`.
*   `listCreated, listCreatedLoc`: Logical arrays to track if a linked list for a particular column block has been initiated.
*   Key utility/kernel routines:
    *   `psp_copy_spm2st`: Extracts a panel from a `psp_matrix_spm`.
    *   `psp_sst_gespmspm`: Performs local multiplication `op(A_panel) * B_local_part`.
    *   `psp_MPI_REDUCE_spm_packed`: Custom MPI reduction for packed sparse matrix data.
    *   `psp_spm2lists_shift`: Converts columns of a sparse panel to linked list segments with global index shifting.
    *   `psp_list2spm`: Converts the aggregated linked list(s) back to a CSC sparse matrix representation.
    *   `psp_sst_sum_spmspm`: Performs the final sparse matrix sum `alpha*RESULT + beta*C_initial`.
    *   `list_destroy`: Deallocates a linked list.

# Usage Examples

The `psp_dgespmspm_tn` and `psp_zgespmspm_tn` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmspm`, `psp_zgespmspm`) located within the `pspSpmSpm` module. This occurs when the `opA` argument is 'T' (transpose) or 'C' (conjugate transpose) and `opB` is 'N' (normal). For end-user examples, please refer to the documentation for the `pspSpmSpm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for `TYPE(psp_matrix_spm)`, `TYPE(psp_MPI_dspm/zspm)`, and `TYPE(dList/zList/dListPtrArray/zListPtrArray)` definitions, as well as list manipulation utilities like `list_destroy`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st`, `psp_idx_glb2loc`, `psp_spm2lists_shift`, `psp_list2spm`. The core computational kernels `psp_sst_gespmspm` and `psp_sst_sum_spmspm` are also expected from here or `psp_spBLAS_Level3`.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and provides `psp_MPI_REDUCE_spm_packed`.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but not directly invoked.
*   **`pspMatSum`**: While `psp_sum_spmspm` is logically from `pspMatSum`, the call here is to `psp_sst_sum_spmspm`, suggesting it's a lower-level kernel from `pspUtility` or similar.
*   **MPI Library**: `MPI_Bcast` is used directly. Custom reduction logic is encapsulated in `psp_MPI_REDUCE_spm_packed`.
*   **BLACS Library**: `blacs_gridinfo` is called.
*   This module implements `C_sparse = op(A)_sparse * B_sparse + C_sparse` where `op(A)` involves a transpose. It uses a complex distributed algorithm with panel broadcasts of sparse matrix `A`, local sparse-sparse multiplications, a custom sparse panel reduction, and an array of linked lists for assembling column blocks of the result matrix `C`.
