# Overview

The `pspMspm_nt.F90` file provides the `pspMspm_nt` module, a specialized component within the pspBLAS Level 3 library. This module implements the dense-matrix times sparse-matrix multiplication where the dense matrix `A` is in its normal (non-transposed) form and the sparse matrix `B` is transposed or conjugate-transposed. The operation is `C := alpha*A*op(B) + beta*C`, with `op(B)` being `B^T` or `B^H`, and `C` being a dense matrix. This module provides the "_nt" (Normal-Transpose) specific implementations called by the generic `psp_gemspm` interface.

The algorithm employed is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach, adapted for the specific requirements of multiplying a dense matrix by a transposed sparse matrix on distributed memory architectures.

# Key Components

*   **Module `pspMspm_nt`**:
    *   **Public Interface `psp_gemspm_nt`**: This interface makes the "normal-transpose" dense-sparse multiplication routines available.
        *   `psp_dgemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. `A` is an M-by-K dense matrix. `B` is an N-by-K sparse matrix (so `op(B)` is K-by-N). `C` is an M-by-N dense matrix.
            *   `A`: A 2D Fortran array for the local part of the distributed dense input matrix `A`.
            *   `opA`: Character flag, must be 'N' (or 'n') for these routines.
            *   `B`: `TYPE(psp_matrix_spm)`, the local part of the distributed N-by-K sparse input matrix `B`.
            *   `opB`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `C`: A 2D Fortran array (inout), the local part of the distributed dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm (SUMMA-like for A_dense * op(B_sparse))**:
            1.  **Grid Setup**: Retrieves BLACS process grid information (`iprow`, `ipcol`, `nprow`, `npcol`) from `psp_icontxt`.
            2.  **Local Allocations**: Allocates local work arrays: `C_loc` (for a panel of C being computed by the current process row/column before reduction), `CC_loc` (to receive reduced results for `C_loc`), and `B_loc_idx1`, `B_loc_idx2`, `B_loc_val` for storing and broadcasting panels of the sparse matrix `B` (likely in CSC format for the panel).
            3.  **Outer Loop (`kloop`)**: Iterates over blocks of the dimension `N` (columns of `op(B)` and `C`) in chunks of size `psp_update_rank`.
                *   **Sparse Panel Broadcast**:
                    *   The process row owning the current panel of `B` (rows of `B` that become columns of `op(B)`) uses `psp_copy_spm2st` to extract this panel into the `B_loc_*` arrays.
                    *   These arrays representing the sparse panel are broadcast along the process column using `MPI_Bcast` via `psp_mpi_comm_col`.
                *   **Local Computation (`psp_sst_gemspm`)**: Each process computes its partial contribution to the output panel `C_loc`. It calls `psp_sst_gemspm` (an external "Subroutine Subprogram Template") which performs the local dense matrix (`A`) times transposed sparse panel (`B_loc_panel`) multiplication: `C_loc = A * op(B_loc_panel)`.
                *   **Reduction**: The partial results in `C_loc` (from processes in the same process row, corresponding to different blocks of `A`) are summed up using `MPI_REDUCE` into `CC_loc` on the process column `idx_pcol` that owns the current panel of `C`. The reduction happens along `psp_mpi_comm_row`.
                *   **Global Update**: The process column `idx_pcol` updates its portion of the global matrix `C` using the reduced result `CC_loc`: `C_panel = beta*C_panel + alpha*CC_loc_panel`.
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Provides essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes), expected to be initialized by `pspMPI`.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `psp_copy_spm2st`: Utility function (from `pspUtility`) to extract a panel from a `psp_matrix_spm` object into a standard sparse array format for broadcasting.
*   `psp_sst_gemspm`: An external subroutine (not defined in this module) that performs the core local multiplication of a dense matrix by a transposed sparse matrix panel. This is the computational kernel for the local block operations.
*   `B_loc_idx1`, `B_loc_idx2`, `B_loc_val`: Temporary arrays for holding and broadcasting panels of the sparse matrix `B`.

# Usage Examples

The `psp_dgemspm_nt` and `psp_zgemspm_nt` subroutines are primarily intended to be called by the dispatcher routines (`psp_dgemspm`, `psp_zgemspm`) within the `pspMspm` module when `opA` is 'N' and `opB` is 'T' or 'C'. For end-user examples, refer to the documentation for the `pspMspm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Necessary for the `TYPE(psp_matrix_spm)` definition of the sparse matrix `B`.
*   **`pspUtility`**: Relies on `psp_copy_spm2st` for extracting sparse matrix panels and potentially `psp_idx_glb2loc` for index conversions. The crucial `psp_sst_gemspm` routine is also expected to be provided by or through `pspUtility` or a related module like `psp_spBLAS_Level3`.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and the MPI communicators defined therein.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements but their routines are not directly called in this file. They might be used by `psp_sst_gemspm`.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for data distribution and `MPI_Reduce` for accumulating partial results. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to obtain process grid information.
*   This module implements the `C = A*op(B)_sparse + C` operation for `op(B)` being transposed/conjugated. It uses a distributed algorithm where panels of the sparse matrix `B` are broadcast, local dense-sparse products are computed (via `psp_sst_gemspm`), and results are then reduced.
