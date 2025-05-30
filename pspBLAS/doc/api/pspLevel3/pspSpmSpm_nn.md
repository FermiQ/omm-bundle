# Overview

The `pspSpmSpm_nn.F90` file provides the `pspSpmSpm_nn` module, a specialized component of the pspBLAS Level 3 library. This module is dedicated to the multiplication of two sparse matrices, `A` and `B`, where both matrices are used in their normal (non-transposed) forms. The operation is `C := alpha*A*B + beta*C`, where `A`, `B`, and `C` are all sparse matrices of `TYPE(psp_matrix_spm)`. This module contains the "_nn" (Normal-Normal) specific implementations for the `psp_gespmspm` (General Sparse Matrix - Sparse Matrix product) functionality.

The algorithm is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. It involves breaking down the matrices into panels, broadcasting these panels across a 2D process grid, performing local sparse matrix-sparse matrix multiplications on these panels, and then accumulating the results into the output sparse matrix `C`.

# Key Components

*   **Module `pspSpmSpm_nn`**:
    *   **Public Interface `psp_gespmspm_nn`**: This interface exposes the routines for non-transposed sparse-sparse matrix multiplication.
        *   `psp_dgespmspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation. Sparse matrix `A` is M-by-K, sparse matrix `B` is K-by-N, and sparse matrix `C` is M-by-N. All are of `TYPE(psp_matrix_spm)`.
            *   `A, B`: `TYPE(psp_matrix_spm)` (inout). Input sparse matrices. Marked `inout` because the parent `pspSpmSpm` module might convert them from COO to CSC format before calling these routines.
            *   `opA, opB`: Character flags. For this `_nn` module, they are expected to be 'N' (or 'n'), indicating no transpose for `A` and `B`.
            *   `C`: `TYPE(psp_matrix_spm)` (inout). The output sparse matrix.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Algorithm (SUMMA-like for A_sparse * B_sparse)**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid information (`iprow`, `ipcol`, `nprow`, `npcol`) using `psp_icontxt`.
            2.  **Local Initialization**: A local sparse matrix `C_loc` of `TYPE(psp_matrix_spm)` is initialized to zero with the global dimensions M-by-N and 'csc' format using `psp_spm_zeros`. This matrix will accumulate the local products `A_panel * B_panel`.
            3.  **Outer Loop (`kloop`)**: Iterates over the common dimension `K` in blocks of size `psp_update_rank` (a parameter from the `/psp_grid2D/` common block).
                *   **Panel A Extraction and Broadcast**:
                    *   The process column (`idx_pcol`) that owns the current panel of `A` (columns `glb_st` to `glb_ed` of `A`) extracts this panel into temporary CSC-like arrays (`A_loc_idx1`, `A_loc_idx2`, `A_loc_val`) using `psp_copy_spm2st`.
                    *   These arrays representing the sparse panel of `A` are then broadcast row-wise across the process grid using `MPI_Bcast` on the `psp_mpi_comm_row` communicator.
                *   **Panel B Extraction and Broadcast**:
                    *   Similarly, the process row (`idx_prow`) owning the current panel of `B` (rows `glb_st` to `glb_ed` of `B`) extracts it into CSC-like arrays (`B_loc_idx1`, `B_loc_idx2`, `B_loc_val`) using `psp_copy_spm2st`.
                    *   These arrays are broadcast column-wise using `MPI_Bcast` on the `psp_mpi_comm_col` communicator.
                *   **Local Sparse-Sparse Multiplication (`psp_sst_gespmspm`)**: Each process performs a local multiplication of the received sparse panel of `A` (`A_loc_...`) and the sparse panel of `B` (`B_loc_...`). The result is accumulated into the local sparse matrix `C_loc`. This is done by calling `psp_sst_gespmspm` (an external "Subroutine Subprogram Template") with `alpha_local=1.0_dp` and `beta_local=1.0_dp` (for accumulation).
            4.  **Final Combination**: After all `kloop` iterations, `C_loc` holds the result of `A*B`. The final output `C` is computed by `C = alpha*C_loc + beta*C` using the `psp_sum_spmspm` routine from the `pspMatSum` module.
            5.  **Cleanup**: Deallocates `C_loc` and the temporary arrays used for broadcasting panels.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes like `psp_update_rank`).
*   `numroc`: External ScaLAPACK function for calculating local dimensions.
*   `A_loc_idx1, A_loc_idx2, A_loc_val`: Temporary arrays for holding panels of sparse matrix `A` in CSC format for broadcasting.
*   `B_loc_idx1, B_loc_idx2, B_loc_val`: Temporary arrays for holding panels of sparse matrix `B` in CSC format for broadcasting.
*   `C_loc`: A local `TYPE(psp_matrix_spm)` variable used to accumulate the sum of local panel products `A_panel * B_panel`.
*   `psp_copy_spm2st`: A utility function (from `pspUtility`) to extract a panel from a `psp_matrix_spm` object into standard sparse arrays (CSC format).
*   `psp_spm_zeros`: A utility function (from `pspVariable`) to initialize a `psp_matrix_spm` object to represent a zero matrix of specified dimensions and format.
*   `psp_sst_gespmspm`: An external subroutine (expected to be in `pspUtility` or `psp_spBLAS_Level3`) that performs the core local multiplication of two sparse matrices, both provided in CSC format.
*   `psp_sum_spmspm`: A routine from the `pspMatSum` module used for the final scaling and addition: `C = alpha*C_loc + beta*C`.
*   `psp_deallocate_spm`: A utility function (from `pspVariable`) to deallocate the components of a `psp_matrix_spm` object.

# Usage Examples

The `psp_dgespmspm_nn` and `psp_zgespmspm_nn` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmspm`, `psp_zgespmspm`) located within the `pspSpmSpm` module. This occurs when both `opA` and `opB` arguments are 'N' (indicating no transpose). For end-user examples, please refer to the documentation for the `pspSpmSpm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition and utilities like `psp_spm_zeros` and `psp_deallocate_spm`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st` for extracting sparse matrix panels. The core computational kernel `psp_sst_gespmspm` is also expected to be provided by or through `pspUtility` or a related module like `psp_spBLAS_Level3`.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and the MPI communicators it defines.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements, but their routines are not directly called in this file.
*   **`pspMatSum`**: Uses `psp_sum_spmspm` for the final combination of the computed product `C_loc` with the input/output matrix `C`.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for distributing panels of the sparse matrices. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to query information about the process grid.
*   This module implements the `C_sparse = A_sparse * B_sparse + C_sparse` operation (all non-transposed) using a distributed SUMMA-like algorithm. It involves extracting and broadcasting sparse matrix panels, performing local sparse-sparse multiplications via `psp_sst_gespmspm`, and then accumulating the results using `psp_sum_spmspm`.
