# Overview

The `pspMspm_tn.F90` file provides the `pspMspm_tn` module, which is a specialized part of the pspBLAS Level 3 library. It focuses on the multiplication of a transposed (or conjugate-transposed) dense matrix `A` with a sparse matrix `B` that is in its normal (non-transposed) form. The operation is `C := alpha*op(A)*B + beta*C`, where `op(A)` represents `A^T` or `A^H`, and `C` is a dense matrix. This module contains the "_tn" (Transpose-Normal) implementations for the `psp_gemspm` (General Matrix - Sparse Matrix product) functionality.

The underlying algorithm is a parallel SUMMA-like (Scalable Universal Matrix Multiplication Algorithm) approach, adapted for distributed memory systems to handle the dense-transpose times sparse multiplication.

# Key Components

*   **Module `pspMspm_tn`**:
    *   **Public Interface `psp_gemspm_tn`**: This interface provides the routines for multiplying a transposed dense matrix by a non-transposed sparse matrix.
        *   `psp_dgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Real double precision version.
        *   `psp_zgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)`: Complex double precision version.
        *   **Arguments**:
            *   `M, N, K`: Integers representing global dimensions. Dense matrix `A` is globally K-by-M (so `op(A)` is M-by-K). Sparse matrix `B` (of `TYPE(psp_matrix_spm)`) is K-by-N. The resulting dense matrix `C` is M-by-N.
            *   `A`: A 2D Fortran array for the local part of the distributed K-by-M dense input matrix `A`.
            *   `opA`: Character flag, must be 'T' (transpose) or 'C' (conjugate transpose) for these routines.
            *   `B`: A variable of `TYPE(psp_matrix_spm)`, representing the local part of the distributed K-by-N sparse input matrix `B`.
            *   `opB`: Character flag, must be 'N' (or 'n') for these routines.
            *   `C`: A 2D Fortran array (inout), the local part of the distributed M-by-N dense output matrix `C`.
            *   `alpha, beta`: Scalar multipliers (real or complex).
        *   **Algorithm (SUMMA-like for op(A_dense) * B_sparse)**:
            1.  **Grid Information**: Fetches BLACS process grid details (`iprow`, `ipcol`, `nprow`, `npcol`) via `psp_icontxt`.
            2.  **Local Allocations**: Allocates local work arrays: `A_loc` (for a panel of `A`), `C_loc` (for a panel of `C` computed locally before reduction), and `CC_loc` (to store results after `MPI_REDUCE`). Note that unlike `pspMspm_nn` or `pspMspm_nt`, temporary arrays for broadcasting parts of sparse matrix `B` (`B_loc_val`, etc.) are allocated but not explicitly populated by panel extraction and broadcast in this top-level routine; `B`'s components are passed directly to `psp_sst_gemspm`.
            3.  **Outer Loop (`kloop`)**: Iterates over blocks of the dimension `M` (which corresponds to the rows of `op(A)` and `C`) in chunks of size `psp_update_rank`.
                *   **Dense Panel Broadcast**:
                    *   The process column (`idx_pcol`) that owns the current panel of `A` (columns of `A` which form rows of `op(A)`) copies its local data into `A_loc`.
                    *   This `A_loc` panel is then broadcast row-wise across the process grid using `MPI_Bcast` on the `psp_mpi_comm_row` communicator.
                *   **Local Computation (`psp_sst_gemspm`)**: Each process computes its partial contribution to a panel of the output matrix `C`, storing it in `C_loc`. This is done by calling `psp_sst_gemspm`. This external subroutine is expected to perform the local multiplication of the transposed `A_loc` panel with the relevant local part of the sparse matrix `B` (i.e., `op(A_loc) * B_local_part`). The full `B%row_ind`, `B%col_ptr`, `B%dval` (or `B%zval`) are passed to `psp_sst_gemspm`, suggesting it handles accessing the correct portion of the sparse matrix `B` that interacts with `op(A_loc)`.
                *   **Reduction**: The partial results `C_loc` (computed by processes in the same process column, corresponding to different blocks of `B` contributing to the same panel of `C`) are summed using `MPI_REDUCE` along the `psp_mpi_comm_col`. The summed result for the panel is stored in `CC_loc` on the specific process row (`idx_prow`) that owns this panel of the global matrix `C`.
                *   **Global Update**: The processes in the designated row `idx_prow` update their local portion of the global matrix `C` using the reduced result from `CC_loc`: `C_panel = beta*C_panel + alpha*CC_loc_panel`.
    *   **Subroutine `die(message)`** (private): Standard error handling routine.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   **COMMON Block `/psp_grid2D/`**: Contains essential parameters for the distributed algorithm (BLACS context, grid dimensions, MPI communicators, block sizes), expected to be initialized by `pspMPI`.
*   `numroc`: External ScaLAPACK function for local dimension calculations.
*   `psp_sst_gemspm`: An external subroutine (not defined in this module) that performs the core local multiplication of a transposed dense matrix panel by a sparse matrix (or its relevant local block). This is the computational kernel.

# Usage Examples

The `psp_dgemspm_tn` and `psp_zgemspm_tn` subroutines are primarily intended to be called by the higher-level dispatcher routines (`psp_dgemspm`, `psp_zgemspm`) located within the `pspMspm` module. This occurs when the `opA` argument is 'T' (transpose) or 'C' (conjugate transpose) and `opB` is 'N' (normal). For end-user examples, please refer to the documentation for the `pspMspm` module.

# Dependencies and Interactions

*   **`pspVariable`**: Essential for the definition of `TYPE(psp_matrix_spm)` used for the sparse matrix `B`.
*   **`pspUtility`**: Likely for `psp_idx_glb2loc` (for index conversions) and the definition/interface of the crucial `psp_sst_gemspm` routine.
*   **`pspMPI`**: Depends on the initialized common block `/psp_grid2D/` and the MPI communicators it defines.
*   **`pspLevel1`, `pspLevel2`**: These modules are included via `USE` statements, but their specific routines are not directly called within this file. They might be dependencies for `psp_sst_gemspm`.
*   **MPI Library**: Makes direct calls to `MPI_Bcast` for distributing panels of matrix `A` and `MPI_Reduce` for accumulating partial results of `C`. `mpif.h` is included if `HAVE_MPI` is defined.
*   **BLACS Library**: `blacs_gridinfo` is called to query information about the process grid.
*   This module implements the `C = op(A)_dense * B_sparse + C` operation where `op(A)` involves a transpose. It uses a distributed algorithm where panels of the dense matrix `A` are broadcast. The local multiplication `op(A_panel) * B_local_part` is handled by `psp_sst_gemspm`, and these partial results are then reduced across process columns.
