# Overview

The `pspMspm_nn.F90` file contains the `pspMspm_nn` module, which provides specialized routines for a specific type of Level 3 pspBLAS operation: the multiplication of a dense matrix `A` by a sparse matrix `B`, where both matrices are used in their normal (non-transposed) form. The operation performed is `C := alpha*A*B + beta*C`, where `C` is a dense matrix. This module implements the "nn" (Normal-Normal) case for the `psp_gemspm` (General Matrix - Sparse Matrix product) functionality.

The implementation uses a SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach, adapted for the dense-sparse case. This involves distributing panels of the dense matrix `A` and the sparse matrix `B` across a 2D process grid and performing local computations, followed by necessary communications.

# Key Components

*   **Module `pspMspm_nn`**:
    *   **Public Interface `psp_gemspm_nn`**: This interface exposes the dense-sparse multiplication routines for the non-transposed case.
        *   `psp_dgemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation: dense matrix `A` is M-by-K, sparse matrix `B` (of `TYPE(psp_matrix_spm)`) is K-by-N, and dense matrix `C` is M-by-N.
            *   `A`: A 2D Fortran array, representing the local segment of the distributed dense input matrix `A`.
            *   `opA, opB`: Character flags. For this `_nn` module, they are expected to be 'N' (or 'n'), indicating no transpose.
            *   `B`: A variable of `TYPE(psp_matrix_spm)`, representing the local segment of the distributed sparse input matrix `B`.
            *   `C`: A 2D Fortran array (inout), representing the local segment of the distributed dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Algorithm (SUMMA-like for A_dense * B_sparse)**:
            1.  **Grid Information**: Retrieves the current process's coordinates (`iprow`, `ipcol`) and the total grid dimensions (`nprow`, `npcol`) from the BLACS context (`psp_icontxt`) stored in the `/psp_grid2D/` common block.
            2.  **Local Allocations**: Allocates local arrays: `A_loc` for holding panels of `A`, `C_loc` for accumulating local results for `C`, and `B_loc_val`, `B_loc_idx1`, `B_loc_idx2` to temporarily store panels of the sparse matrix `B` in a CSC-like (Compressed Sparse Column) format for broadcasting.
            3.  **Outer Loop (`kloop`)**: Iterates over the common dimension `K` in blocks of size `psp_update_rank` (a parameter from `/psp_grid2D/`).
                *   **Panel Preparation & Broadcast of A**: The process column that owns the current panel of `A` copies its local portion into `A_loc`. This `A_loc` panel is then broadcast row-wise across the process grid using `MPI_Bcast` on the `psp_mpi_comm_row` communicator.
                *   **Panel Preparation & Broadcast of B**: The process row that owns the corresponding panel of sparse matrix `B` (rows `glb_st` to `glb_ed`) extracts this panel into the temporary CSC-like arrays (`B_loc_idx1`, `B_loc_idx2`, `B_loc_val`) using `psp_copy_spm2st`. These arrays (and the panel width `width`) are then broadcast column-wise using `MPI_Bcast` on the `psp_mpi_comm_col` communicator.
                *   **Local Computation**: Each process performs a local multiplication of its received `A_loc` panel with the received sparse `B` panel. This is not a direct call to a dense BLAS `dgemm` for the `A_loc * B_sparse_panel` part but is implemented as loops over the columns and non-zero elements of the sparse `B` panel:
                    For each column `idx_col` in the local `C_loc` block:
                    For each non-zero element `B_loc_val(i)` in that column of the sparse `B` panel (with row index `B_loc_idx1(i)` relative to the panel's start):
                    `C_loc(:,idx_col) = C_loc(:,idx_col) + B_loc_val(i) * A_loc(:, B_loc_idx1(i))`
                    This effectively performs `C_loc = C_loc + A_loc * B_sparse_panel`.
            4.  **Final Update**: After all `kloop` iterations, the accumulated local result `C_loc` is scaled by `alpha`, and added to the input `C` (scaled by `beta`): `C = beta*C + alpha*C_loc`.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm, such as BLACS context (`psp_icontxt`), process grid dimensions (`psp_nprow`, `psp_npcol`), MPI communicators (`psp_mpi_comm_row`, `psp_mpi_comm_col`), and block sizes (`psp_bs_def_row`, `psp_bs_def_col`, `psp_update_rank`). This common block is expected to be initialized by `pspMPI` routines.
*   `numroc`: External ScaLAPACK function for calculating local dimensions of distributed arrays.
*   `psp_copy_spm2st`: A utility function (from `pspUtility`) used to extract a panel from a `psp_matrix_spm` object into separate arrays representing a standard sparse format (likely CSC for the panel).
*   `B_loc_idx1`, `B_loc_idx2`, `B_loc_val`: Temporary arrays used to hold and broadcast panels of the sparse matrix `B` in a CSC-like format.

# Usage Examples

The `psp_dgemspm_nn` and `psp_zgemspm_nn` subroutines are typically not called directly by end-users. They are invoked by the `psp_dgemspm` or `psp_zgemspm` dispatcher routines within the `pspMspm` module when both `opA` and `opB` arguments are 'N' (indicating no transpose). For usage, refer to the example provided for the `pspMspm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Required for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `B`.
*   **`pspUtility`**: Uses `psp_idx_glb2loc` (implicitly via `psp_copy_spm2st` or directly, for mapping global indices to local) and `psp_copy_spm2st` (to extract sparse matrix panels).
*   **`pspMPI`**: Relies on the variables (BLACS context, communicators, grid dimensions) defined and initialized within this module, accessed via the `/psp_grid2D/` common block.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements, but their routines are not directly called in this specific file. They might be used by utility functions or provide a broader context.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for data distribution in the SUMMA-like algorithm. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to get information about the process grid.
*   This module implements the core computational logic for the non-transposed dense-matrix times sparse-matrix multiplication using a distributed, panel-based algorithm. It forms a key part of the pspBLAS Level 3 functionality.
