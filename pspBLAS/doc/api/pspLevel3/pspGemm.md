# Overview

The `pspGemm.F90` file is a core part of the pspBLAS library, focusing on Level 3 BLAS functionality. Specifically, it provides implementations for the general matrix-matrix multiplication (GEMM) operation: `C := alpha*op(A)*op(B) + beta*C`. Here, `op(X)` can be `X`, `X^T` (transpose), or `X^H` (conjugate transpose for complex matrices).

The module is designed for distributed memory parallel systems and employs the SUMMA (Scalable Universal Matrix Multiplication Algorithm) for its parallel GEMM routines. It offers interfaces for both real (`psp_dgemm`) and complex (`psp_zgemm`) arithmetic. These top-level routines then dispatch to specialized internal routines based on the transposition specified for matrices `A` and `B` (e.g., `_nn` for no transposes, `_nt` for `A` no transpose and `B` transpose, etc.).

# Key Components

*   **Module `pspGemm`**:
    *   **Public Interface `psp_gemm`**: This is the primary user-facing interface for GEMM operations.
        *   `psp_dgemm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Performs real double-precision GEMM.
        *   `psp_zgemm(M,N,K,A,opA,B,opB,C,alpha,beta)`: Performs complex double-precision GEMM.
        *   **Arguments**:
            *   `M, N, K`: Global dimensions of the matrices after applying `opA` and `opB`. `op(A)` is M-by-K, `op(B)` is K-by-N, `C` is M-by-N.
            *   `A, B`: Local portions of the distributed input matrices.
            *   `opA, opB`: Character flags ('N' for normal, 'T' for transpose, 'C' for conjugate transpose) indicating the operation on A and B.
            *   `C`: Local portion of the distributed output matrix.
            *   `alpha, beta`: Scalar multipliers.
        *   **Functionality**: These routines parse `opA` and `opB` (using `psp_process_opM` from `pspUtility`) and then call one of the specialized `_nn, _nt, _tn, _tt` routines.

    *   **Private Interfaces & SUMMA Implementations**: For each combination of `op(A)` and `op(B)`:
        *   `psp_gemm_nn` (`A*B`): `psp_dgemm_nn`, `psp_zgemm_nn`
        *   `psp_gemm_nt` (`A*B^T`): `psp_dgemm_nt`, `psp_zgemm_nt`
        *   `psp_gemm_tn` (`A^T*B`): `psp_dgemm_tn`, `psp_zgemm_tn`
        *   `psp_gemm_tt` (`A^T*B^T`): `psp_dgemm_tt`, `psp_zgemm_tt`
        *   **SUMMA Algorithm**: These routines implement the SUMMA algorithm. Key steps include:
            1.  Retrieving 2D process grid information via `blacs_gridinfo`.
            2.  Iterating (`kloop`) over blocks of the common dimension `K` in chunks of size `psp_update_rank`.
            3.  In each iteration:
                *   The process column (or row, depending on transpose) that owns the current block of `A` broadcasts it along its process row.
                *   The process row (or column) that owns the current block of `B` broadcasts it along its process column.
                *   Each process performs a local DGEMM/ZGEMM on the received blocks `A_loc` and `B_loc`, accumulating the result into a local portion of `C` (`C_loc`).
            4.  The `psp_dgemm_tt` and `psp_zgemm_tt` routines may choose an optimal strategy, potentially transposing one of the matrices first (using ScaLAPACK's `pdtran` or `pztranc`/`pztranu`) if it's deemed more efficient, and then calling another `psp_gemm_xx` variant.

    *   **Subroutine `die(message)`** (private): An error handling routine that writes to "MatrixSwitch.log" and stops execution.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   **COMMON Block `/psp_grid2D/`**: This is critical for the SUMMA implementations. It stores global variables describing the BLACS 2D process grid and communication setup, initialized by `pspMPI` routines (e.g., `psp_gridinit_2D`):
    *   `psp_mpi_comm_world`, `psp_mpi_size`, `psp_nprow` (number of process rows), `psp_npcol` (number of process columns).
    *   `psp_bs_def_row`, `psp_bs_def_col`: Default block sizes for matrix distribution.
    *   `psp_update_rank`: The block size used in the SUMMA algorithm's main loop for panel broadcasts.
    *   `psp_icontxt`: The BLACS context handle for the 2D grid.
    *   `psp_mpi_comm_cart`, `psp_mpi_comm_row`, `psp_mpi_comm_col`: MPI communicators for the Cartesian grid, and for row-wise and column-wise broadcasts, respectively.
*   `numroc`: External ScaLAPACK function used to calculate the local number of rows/columns a process owns for a distributed matrix.

# Usage Examples

The end-user of pspBLAS would typically initialize the pspBLAS environment (including MPI and BLACS setup via `psp_gridinit_2D` or similar from `pspMPI`), distribute their matrices according to the defined block sizes, and then call `psp_gemm`:

```fortran
USE pspBLAS ! This makes psp_gemm available (as it's public in pspLevel3, which is used by pspBLAS)
IMPLICIT NONE

REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: A_local_data, B_local_data, C_local_data
INTEGER :: global_M, global_N, global_K
INTEGER :: local_rows_A, local_cols_A, local_rows_B, local_cols_B, local_rows_C, local_cols_C
REAL(dp) :: alpha_val, beta_val

! 1. Initialize pspBLAS (MPI, BLACS grid using psp_gridinit_2D from pspMPI)
!    This will populate the /psp_grid2D/ common block.

! 2. Define global matrix dimensions (global_M, global_N, global_K)
!    and determine local dimensions for A, B, C using numroc and psp_bs_def_row/col.
!    Allocate and initialize A_local_data, B_local_data, C_local_data.

alpha_val = 2.0_dp
beta_val = -1.0_dp

! Example: Compute C = alpha * A^T * B + beta * C
CALL psp_gemm(global_M, global_N, global_K, &
              A_local_data, 'T', &
              B_local_data, 'N', &
              C_local_data, &
              alpha_val, beta_val)

! ... C_local_data now contains the updated result ...

! 3. Finalize pspBLAS (e.g., psp_gridexit_2D)
```

# Dependencies and Interactions

*   **`pspVariable`**: Likely used for internal type definitions related to distributed matrices, although not directly visible in the subroutine signatures which take Fortran arrays.
*   **`pspUtility`**: For utility functions like `psp_process_opM` (to parse 'N', 'T', 'C' flags) and `psp_idx_glb2loc` (to convert global indices/dimensions to local ones for panel copying).
*   **`pspMPI`**: Essential for providing the initialized `/psp_grid2D/` common block variables and MPI communicators (`psp_mpi_comm_row`, `psp_mpi_comm_col`).
*   **`pspLevel1`, `pspLevel2`**: These modules are `USE`d but their functionalities are not directly invoked within the `pspGemm.F90` code itself. They might be used by other Level 3 routines or utilities.
*   **MPI Library**: Direct calls to `MPI_Bcast` and `MPI_Reduce` are made for inter-process communication during the SUMMA algorithm. The `mpif.h` header is included.
*   **BLACS Library**: `blacs_gridinfo` is called to get process coordinates and grid dimensions. `descinit` is used in `psp_dgemm_tt`/`psp_zgemm_tt` for ScaLAPACK calls.
*   **Standard BLAS (Level 3)**: Local matrix multiplications within each process in SUMMA are performed by calling standard `dgemm` or `zgemm`.
*   **ScaLAPACK**: The `numroc` function is used for dimension calculations. The `_tt` variants use ScaLAPACK's parallel transpose routines (`pdtran`, `pztranc`, `pztranu`).
*   The module heavily relies on a correctly initialized BLACS 2D process grid and the corresponding MPI communicators, which are expected to be set up by routines in `pspMPI` and stored in the `/psp_grid2D/` common block.
