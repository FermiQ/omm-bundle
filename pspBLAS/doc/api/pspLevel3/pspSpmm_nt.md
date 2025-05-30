# Overview

The `pspSpmm_nt.F90` file provides the `pspSpmm_nt` module, a specialized component of the pspBLAS Level 3 library. This module is designed for multiplying a sparse matrix `A` (in its normal, non-transposed form) by a dense matrix `B` which is transposed or conjugate-transposed. The operation is `C := alpha*A*op(B) + beta*C`, where `C` is a dense matrix. This module contains the "_nt" (Normal-Transpose) specific implementations for the `psp_gespmm` (General Sparse Matrix - Dense Matrix product) functionality.

The core algorithm is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach. This involves broadcasting panels of the (transposed) dense matrix `B` across a 2D process grid. Each process then performs a local multiplication of its portion of sparse matrix `A` with the received panel of `op(B)`. These partial results are subsequently reduced across process rows to form the final result panels of matrix `C`. The module includes an optimization path to use Intel MKL routines for the local sparse-dense matrix products if `HAVE_MKL` is defined.

# Key Components

*   **Module `pspSpmm_nt`**:
    *   **Public Interface `psp_gespmm_nt`**: This interface exposes the routines for normal-sparse times transposed-dense matrix multiplication.
        *   `psp_dgespmm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgespmm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Sparse matrix `A` (of `TYPE(psp_matrix_spm)`) is M-by-K. Dense matrix `B` is N-by-K (so `op(B)` is K-by-N). The resulting dense matrix `C` is M-by-N.
            *   `A`: `TYPE(psp_matrix_spm)`, the input sparse matrix. Assumed to be in CSC format (or converted by the calling `pspSpmm` module).
            *   `opA`: Character flag, must be 'N' (or 'n') for these routines.
            *   `B`: A 2D Fortran array, representing the local part of the distributed N-by-K dense input matrix `B`.
            *   `opB`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `C`: A 2D Fortran array (inout), representing the local part of the distributed M-by-N dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm (SUMMA-like for A_sparse * op(B_dense))**:
            1.  **Grid Information**: Retrieves BLACS 2D process grid details (`iprow`, `ipcol`, `nprow`, `npcol`).
            2.  **Local Allocations**: Allocates local work arrays: `B_loc` (for a panel of `op(B)`), `C_loc` (for a panel of `C` computed locally before reduction), and `CC_loc` (to store results from `MPI_REDUCE`). Temporary arrays for sparse `A` panels (`A_loc_val`, etc.) are declared but not used for broadcasting `A` in this `_nt` variant (as `A` is used locally by each process row).
            3.  **Outer Loop (`kloop`)**: Iterates over the dimension `N` (columns of `op(B)` and `C`) in blocks of size `psp_update_rank`.
                *   **Panel B (Dense, Transposed) Preparation and Broadcast**:
                    *   The process row (`idx_prow`) that owns the current panel of `B` (rows of `B`, which form columns of `op(B)`) copies its local data into `B_loc`.
                    *   **If `HAVE_MKL` is defined**: `B_loc` is explicitly transposed (and conjugated if `opB='C'`) from the corresponding segment of `B` to match MKL's input requirements for `op(B)`.
                    *   This `B_loc` panel (now representing `op(B)_panel`) is then broadcast column-wise using `MPI_Bcast` on `psp_mpi_comm_col`.
                *   **Local Computation**: Each process computes its partial contribution to a panel of `C`, storing it in `C_loc`. This involves multiplying the local part of sparse matrix `A` by the received `B_loc` (panel of `op(B)`).
                    *   **If `HAVE_MKL`**: Calls Intel MKL's `mkl_dcscmm` (real) or `mkl_zcscmm` (complex). The sparse matrix `A` (using `A%dval`/`A%zval`, `A%row_ind`, `A%col_ptr`, and derived `A_loc_idx3` for end pointers) is multiplied by the dense `B_loc` panel.
                    *   **Else (Generic Implementation)**: Calls `psp_sst_gespmm` (an external "Subroutine Subprogram Template", likely from `pspUtility` or `psp_spBLAS_Level3`). This routine performs the local `A_sparse_local_part * op(B_panel)` multiplication, where `op(B)` is handled by `psp_sst_gespmm` using the non-MKL `B_loc` (which is not pre-transposed in this path).
                *   **Reduction**: The partial results `C_loc` (computed by processes in the same process row, corresponding to different blocks of `A`) are summed using `MPI_REDUCE` along the `psp_mpi_comm_row`. The summed result for the panel is stored in `CC_loc` on the specific process column (`idx_pcol`) that owns this panel of the global matrix `C`.
                *   **Global Update**: The processes in the designated column `idx_pcol` update their local portion of the global matrix `C` using the reduced result from `CC_loc`: `C_panel = beta*C_panel + alpha*CC_loc_panel`.
            4.  **Cleanup**: Deallocates all temporary local arrays.

    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `B_loc`: Temporary array for holding panels of the (potentially pre-transposed for MKL) dense matrix `B`.
*   `C_loc`, `CC_loc`: Temporary arrays for local accumulation and then for receiving reduced results of panels of `C`.
*   `A%dval`/`A%zval`, `A%row_ind`, `A%col_ptr`: Components of the input sparse matrix `A`, assumed to be in CSC format for local operations.
*   `A_loc_idx3` (Conditional, if `HAVE_MKL`): Temporary array for MKL CSC format (column end pointers for matrix A).
*   `mkl_dcscmm`, `mkl_zcscmm`: Intel MKL library routines for sparse matrix-dense matrix multiplication.
*   `psp_sst_gespmm`: An external subroutine (expected from `pspUtility` or `psp_spBLAS_Level3`) that performs the core local multiplication of a sparse matrix by a (transposed) dense matrix panel if MKL is not used.

# Usage Examples

The `psp_dgespmm_nt` and `psp_zgespmm_nt` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgespmm`, `psp_zgespmm`) located within the `pspSpmm` module. This occurs when the `opA` argument is 'N' (normal) and `opB` is 'T' (transpose) or 'C' (conjugate transpose). For end-user examples, please refer to the documentation for the `pspSpmm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `A`.
*   **`pspUtility`**: Relies on `psp_idx_glb2loc` for index conversions. The kernel `psp_sst_gespmm` (if MKL not used) is also expected from here or a related module.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and its MPI communicators.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but not directly invoked in this file.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for distributing panels of matrix `B` and `MPI_Reduce` for accumulating partial results of `C`. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to query process grid information.
*   **Intel MKL (Optional)**: If `HAVE_MKL` is defined at compile time, this module will use MKL's `mkl_dcscmm` or `mkl_zcscmm`. Note that if MKL is used, this routine explicitly transposes the `B_loc` panel before calling MKL, as MKL's CSCMM typically expects the dense matrix in normal form.
*   This module implements the `C_dense = A_sparse * op(B)_dense + C_dense` operation (where A is normal, B is transposed) using a distributed SUMMA-like algorithm. It handles panel-wise broadcasting of the dense matrix `B` and uses either MKL or a generic kernel for local sparse-dense computations, followed by result reduction.
