# Overview

The `pspSpmm_nn.F90` file provides the `pspSpmm_nn` module, a specialized component of the pspBLAS Level 3 library. This module focuses on the multiplication of a sparse matrix `A` by a dense matrix `B`, where both matrices are used in their normal (non-transposed) forms. The operation is `C := alpha*A*B + beta*C`, where `C` is a dense matrix. This module contains the "_nn" (Normal-Normal) specific implementations for the `psp_gespmm` (General Sparse Matrix - Dense Matrix product) functionality.

The algorithm employed is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. This involves distributing panels of the sparse matrix `A` (in CSC format) and the dense matrix `B` across a 2D process grid. Local multiplications of these panels are then performed on each process, and the results are accumulated into the corresponding local portion of the output matrix `C`. The module includes an optimization to use Intel MKL's sparse-dense matrix multiplication routines (`mkl_dcscmm` or `mkl_zcscmm`) for the local computations if the `HAVE_MKL` preprocessor macro is defined.

# Key Components

*   **Module `pspSpmm_nn`**:
    *   **Public Interface `psp_gespmm_nn`**: This interface exposes the routines for non-transposed sparse-matrix times dense-matrix multiplication.
        *   `psp_dgespmm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing the global dimensions of the operation. Sparse matrix `A` (of `TYPE(psp_matrix_spm)`) is M-by-K. Dense matrix `B` is K-by-N. The resulting dense matrix `C` is M-by-N.
            *   `A`: `TYPE(psp_matrix_spm)`, the input sparse matrix. Assumed to be in CSC format by the time it's used in the core loop (parent module `pspSpmm` handles COO to CSC conversion).
            *   `opA, opB`: Character flags. For this `_nn` module, they are expected to be 'N' (or 'n'), indicating no transpose for `A` and `B`.
            *   `B`: A 2D Fortran array, representing the local part of the distributed dense input matrix `B`.
            *   `C`: A 2D Fortran array (inout), representing the local part of the distributed dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers of the appropriate type (real or complex).
        *   **Algorithm (SUMMA-like for A_sparse * B_dense)**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid information (`iprow`, `ipcol`, `nprow`, `npcol`) using `psp_icontxt`.
            2.  **Local Allocations**: Allocates local arrays: `B_loc` for holding panels of dense matrix `B`, `C_loc` for accumulating the local result `A_panel*B_panel`, and `A_loc_val`, `A_loc_idx1`, `A_loc_idx2` for storing panels of sparse matrix `A` in CSC format. If `HAVE_MKL` is defined, `A_loc_idx3` (for MKL's CSC end pointers) is also allocated.
            3.  **Outer Loop (`kloop`)**: Iterates over the common dimension `K` in blocks of size `psp_update_rank`.
                *   **Panel A (Sparse) Extraction and Broadcast**:
                    *   The process column (`idx_pcol`) that owns the current panel of sparse `A` (columns `glb_st` to `glb_ed` of `A`) extracts this panel into the CSC arrays (`A_loc_val`, `A_loc_idx1`, `A_loc_idx2`) using `psp_copy_spm2st`.
                    *   These arrays representing the sparse panel of `A` are then broadcast row-wise across the process grid using `MPI_Bcast` on `psp_mpi_comm_row`.
                *   **Panel B (Dense) Extraction and Broadcast**:
                    *   The process row (`idx_prow`) owning the current panel of dense `B` (rows `glb_st` to `glb_ed` of `B`) copies its local portion into `B_loc`.
                    *   This `B_loc` panel is then broadcast column-wise using `MPI_Bcast` on `psp_mpi_comm_col`.
                *   **Local Computation**: Each process computes its contribution: `C_loc_panel = C_loc_panel + A_loc_sparse_panel * B_loc_dense_panel`.
                    *   **If `HAVE_MKL` is defined**: The MKL sparse BLAS routine (`mkl_dcscmm` for real, `mkl_zcscmm` for complex) is called to perform the multiplication of the received sparse `A` panel and dense `B` panel, accumulating the result into `C_loc`. The `matdescr` array is set up for MKL, and `A_loc_idx3` (column end pointers) is derived from `A_loc_idx2` (column start pointers).
                    *   **Else (Generic Implementation)**: A nested loop structure is used. The outer loop iterates from `cnt1 = 1` to `width` (columns in the `A_loc` panel / rows in `B_loc`). The inner loop iterates over the non-zero elements `i` in column `cnt1` of `A_loc_sparse_panel`. The update is `C_loc(A_loc_idx1(i), :) = C_loc(A_loc_idx1(i), :) + A_loc_val(i) * B_loc(cnt1, :)`.
            4.  **Final Update**: After all `kloop` iterations, the accumulated local result `C_loc` is scaled by `alpha`, and added to the input `C` (scaled by `beta`): `C = beta*C + alpha*C_loc`.
            5.  **Cleanup**: Deallocates all temporary local arrays.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes).
*   `numroc`: External ScaLAPACK function for calculating local dimensions of distributed arrays.
*   `A_loc_val, A_loc_idx1, A_loc_idx2`: Temporary arrays for holding panels of sparse matrix `A` in CSC format for broadcasting. `A_loc_idx1` are row indices, `A_loc_idx2` are start pointers for columns.
*   `A_loc_idx3` (Conditional, if `HAVE_MKL`): Temporary array for MKL CSC format (column end pointers).
*   `B_loc`: Temporary array for holding panels of dense matrix `B`.
*   `C_loc`: Temporary array for accumulating the local product `A_panel * B_panel`.
*   `psp_copy_spm2st`: A utility function (from `pspUtility`) to extract a panel from a `psp_matrix_spm` object into standard sparse arrays (CSC format).
*   `mkl_dcscmm`, `mkl_zcscmm`: Intel MKL library routines for sparse matrix-dense matrix multiplication (used if `HAVE_MKL` is defined).

# Usage Examples

The `psp_dgespmm_nn` and `psp_zgespmm_nn` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmm`, `psp_zgespmm`) located within the `pspSpmm` module. This occurs when both `opA` and `opB` arguments are 'N' (indicating no transpose). For end-user examples, please refer to the documentation for the `pspSpmm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `A`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st` for extracting sparse matrix panels and `psp_idx_glb2loc` for global-to-local index conversions.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and the MPI communicators it defines.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but their routines are not directly called in this file.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for distributing panels of matrices. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to query information about the process grid.
*   **Intel MKL (Optional)**: If `HAVE_MKL` is defined at compile time, this module will link against and use MKL's `mkl_dcscmm` or `mkl_zcscmm` for the local sparse-dense matrix products, which can offer significant performance benefits. Otherwise, a generic loop-based implementation is used.
*   This module implements the `C_dense = A_sparse * B_dense + C_dense` operation (where both A and B are non-transposed) using a distributed SUMMA-like algorithm. It handles the panel-wise broadcasting of both sparse and dense matrices and performs local computations, with an optional path for MKL acceleration.
