# Overview

The `pspSpmm_tn.F90` file provides the `pspSpmm_tn` module, a specialized component of the pspBLAS Level 3 library. This module focuses on the multiplication of a transposed (or conjugate-transposed) sparse matrix `A` by a dense matrix `B` in its normal (non-transposed) form. The operation is `C := alpha*op(A)*B + beta*C`, where `C` is a dense matrix. `op(A)` represents `A^T` or `A^H`. This module contains the "_tn" (Transpose-Normal) specific implementations for the `psp_gespmm` (General Sparse Matrix - Dense Matrix product) functionality.

The algorithm is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. Panels of the sparse matrix `A` (which correspond to rows of `op(A)`) are extracted and broadcast across a 2D process grid. Each process then performs a local multiplication of its received transposed sparse panel `op(A_panel)` with its local portion of the dense matrix `B`. These partial results are subsequently reduced across process columns to form the final result matrix `C`. The module includes an option to use Intel MKL routines for the local sparse-dense matrix multiplication if `HAVE_MKL` is defined.

# Key Components

*   **Module `pspSpmm_tn`**:
    *   **Public Interface `psp_gespmm_tn`**: This interface exposes the routines for transposed-sparse-matrix times normal-dense-matrix multiplication.
        *   `psp_dgespmm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Sparse matrix `A` (of `TYPE(psp_matrix_spm)`) is K-by-M (so `op(A)` is M-by-K). Dense matrix `B` is K-by-N. The resulting dense matrix `C` is M-by-N.
            *   `A`: `TYPE(psp_matrix_spm)`, the input K-by-M sparse matrix. Assumed to be in CSC format (or converted by the calling `pspSpmm` module).
            *   `opA`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `B`: A 2D Fortran array, representing the local part of the distributed K-by-N dense input matrix `B`.
            *   `opB`: Character flag, must be 'N' (or 'n') for these routines.
            *   `C`: A 2D Fortran array (inout), representing the local part of the distributed M-by-N dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm (SUMMA-like for op(A_sparse) * B_dense)**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid details (`iprow`, `ipcol`, `nprow`, `npcol`).
            2.  **Local Allocations**: Allocates local work arrays: `C_loc` (for a panel of `C` computed locally before reduction), `CC_loc` (to store results from `MPI_REDUCE`), and `A_loc_val`, `A_loc_idx1`, `A_loc_idx2` for storing panels of sparse matrix `A` in CSC format. If `HAVE_MKL`, `A_loc_idx3` is also allocated. Note that `B_loc` is declared but not used for broadcasting B; the local part of B is used directly.
            3.  **Outer Loop (`kloop`)**: Iterates over the dimension `M` (rows of `op(A)` and `C`) in blocks of size `psp_update_rank`.
                *   **Panel A (Sparse, Transposed) Extraction and Broadcast**:
                    *   The process column (`idx_pcol`) that owns the current panel of `A` (columns of `A`, which correspond to rows of `op(A)`) extracts this panel into the CSC arrays (`A_loc_val`, `A_loc_idx1`, `A_loc_idx2`) using `psp_copy_spm2st`.
                    *   These arrays representing the sparse panel of `A` are then broadcast row-wise across the process grid using `MPI_Bcast` on `psp_mpi_comm_row`.
                *   **Local Computation**: Each process computes its partial contribution to a panel of `C`, storing it in `C_loc`. This involves multiplying the received (and to be operated on as transposed) sparse `A` panel by the relevant local block of dense matrix `B`.
                    *   **If `HAVE_MKL`**: Calls Intel MKL's `mkl_dcscmm` (real) or `mkl_zcscmm` (complex), passing `opA`. The `matdescr` array is set for MKL. `A_loc_idx3` (column end pointers for MKL CSC) is derived from `A_loc_idx2`.
                    *   **Else (Generic Implementation)**: Calls `psp_sst_gespmm` (an external "Subroutine Subprogram Template"). This routine performs the local `op(A_sparse_panel) * B_dense_local_part` multiplication.
                *   **Reduction**: The partial results `C_loc` (computed by processes in the same process column, corresponding to different blocks of `B` contributing to the same panel of `C`) are summed using `MPI_REDUCE` along the `psp_mpi_comm_col`. The summed result for the panel is stored in `CC_loc` on the specific process row (`idx_prow`) that owns this panel of the global matrix `C`.
                *   **Global Update**: The processes in the designated row `idx_prow` update their local portion of the global matrix `C` using the reduced result from `CC_loc`: `C_panel = beta*C_panel + alpha*CC_loc_panel`.
            4.  **Cleanup**: Deallocates all temporary local arrays.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `A_loc_val, A_loc_idx1, A_loc_idx2`: Temporary arrays for holding panels of sparse matrix `A` in CSC format for broadcasting.
*   `A_loc_idx3` (Conditional, if `HAVE_MKL`): Temporary array for MKL CSC format (column end pointers for matrix A panel).
*   `C_loc`, `CC_loc`: Temporary arrays for local accumulation and then for receiving reduced results of panels of `C`.
*   `psp_copy_spm2st`: A utility function (from `pspUtility`) to extract a panel from a `psp_matrix_spm` object into standard sparse arrays (CSC format).
*   `mkl_dcscmm`, `mkl_zcscmm`: Intel MKL library routines for sparse matrix-dense matrix multiplication.
*   `psp_sst_gespmm`: An external subroutine (expected from `pspUtility` or `psp_spBLAS_Level3`) that performs the core local multiplication of a (transposed) sparse matrix panel by a dense matrix block if MKL is not used.

# Usage Examples

The `psp_dgespmm_tn` and `psp_zgespmm_tn` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmm`, `psp_zgespmm`) located within the `pspSpmm` module. This occurs when the `opA` argument is 'T' (transpose) or 'C' (conjugate transpose) and `opB` is 'N' (normal). For end-user examples, please refer to the documentation for the `pspSpmm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `A`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st` for extracting sparse matrix panels and `psp_idx_glb2loc` for index conversions. The kernel `psp_sst_gespmm` (if MKL not used) is also expected from here or a related module.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and its MPI communicators.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but not directly invoked in this file.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for distributing panels of matrix `A` and `MPI_Reduce` for accumulating partial results of `C`. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to query process grid information.
*   **Intel MKL (Optional)**: If `HAVE_MKL` is defined at compile time, this module will use MKL's `mkl_dcscmm` or `mkl_zcscmm`.
*   This module implements the `C_dense = op(A)_sparse * B_dense + C_dense` operation (where A is transposed, B is normal) using a distributed SUMMA-like algorithm. It handles panel-wise broadcasting of the sparse matrix `A` and uses either MKL or a generic kernel for local sparse-dense computations, followed by result reduction.
