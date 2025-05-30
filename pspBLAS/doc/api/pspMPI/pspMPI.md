# Overview

The `pspMPI.F90` file defines the `pspMPI` module, which is a cornerstone of the pspBLAS library's capabilities for parallel computing. This module encapsulates functionalities related to the Message Passing Interface (MPI) and its interplay with BLACS (Basic Linear Algebra Communication Subprograms) for grid-based parallel linear algebra operations.

Key responsibilities of the `pspMPI` module include:
1.  **Processor Grid Initialization**: Provides routines (`psp_gridinit_2D`, `psp_gridinit_3D`) to set up 2D or 3D process grids using MPI and initialize BLACS contexts. These routines populate common blocks with grid parameters and communicators that are then used by other pspBLAS routines.
2.  **Custom Data Types for MPI**: Defines derived types (`psp_MPI_dspm`, `psp_MPI_zspm`) for packaging sparse matrix components (non-zero values and their indices) for efficient MPI communication.
3.  **Custom MPI Reduction Operations**: Implements specialized MPI reduction (summation) operations for both sparse matrices (using packed data or MPI derived types) and dense matrices. These are crucial for parallel algorithms like SUMMA where partial results from different processes need to be combined.

# Key Components

*   **Module `pspMPI`**:
    *   **Derived Types**:
        *   `psp_MPI_dspm`: Used to bundle data for real sparse matrices for MPI operations. It contains:
            *   `nnz`: Integer, local number of non-zero entries.
            *   `loc_dim1`, `loc_dim2`: Integers, local row and column dimensions.
            *   `idx1(:)`, `idx2(:)`: Allocatable integer arrays for row and column indices (or column pointers for CSC).
            *   `val(:)`: Allocatable `REAL(dp)` array for non-zero matrix values.
        *   `psp_MPI_zspm`: Similar to `psp_MPI_dspm`, but `val(:)` is `COMPLEX(dp)` for complex sparse matrices.

    *   **Public Interfaces and Subroutines**:
        *   `psp_gridinit_2D(mpi_comm_in, mpi_size, nprow, order, bs_def_row, bs_def_col, icontxt)`:
            *   Initializes a 2D process grid.
            *   Takes an existing MPI communicator (`mpi_comm_in`), total processes (`mpi_size`), desired number of process rows (`nprow`), grid ordering ('R' for row-major, 'C' for column-major), default block sizes, and an optional existing BLACS context.
            *   Populates the `/psp_grid2D/` common block with details like BLACS context (`psp_icontxt`), grid dimensions (`psp_nprow`, `psp_npcol`), and MPI communicators for the Cartesian grid (`psp_mpi_comm_cart`), rows (`psp_mpi_comm_row`), and columns (`psp_mpi_comm_col`).
        *   `psp_gridinit_3D(mpi_comm_in, mpi_size, nprow, npcol, nplay, order, bs_def_row, bs_def_col, icontxt)`:
            *   Initializes a 3D process grid.
            *   Similar to `psp_gridinit_2D` but adds a layer dimension (`nplay`). Populates `/psp_grid2D/` and `/psp_grid3D/` (which includes `psp_nplay` and `psp_mpi_comm_lay`).
        *   `psp_MPI_REDUCE_spm_packed`: Interface for reducing sparse matrices by packing data into a buffer.
            *   `psp_MPI_REDUCE_dspm_packed(A, idx_proc, isRow)`: Real version. `A` is `TYPE(psp_MPI_dspm)`.
            *   `psp_MPI_REDUCE_zspm_packed(A, idx_proc, isRow)`: Complex version. `A` is `TYPE(psp_MPI_zspm)`.
            *   These routines perform a sum reduction of sparse matrix data (packed using `MPI_PACK`/`MPI_UNPACK`) along either process rows (if `isRow` is true) or columns. The result is accumulated on the process with coordinate `idx_proc` within that row/column. The reduction uses a butterfly-like send/receive pattern.
        *   `psp_MPI_REDUCE_dspm_struct(A, idx_proc, isRow)` and `psp_MPI_REDUCE_zspm_struct(A, idx_proc, isRow)`:
            *   Interfaces intended for reducing sparse matrices of `TYPE(psp_matrix_spm)` using MPI derived data types (`MPI_Type_struct`).
            *   **Note**: These routines are marked with a "TODO: bug inside" comment in the source, indicating they may not be fully functional or correct.
        *   `psp_MPI_REDUCE_den`: Interface for reducing dense matrices.
            *   `psp_MPI_REDUCE_dden(H, m, n)`: Real version.
            *   `psp_MPI_REDUCE_zden(H, m, n)`: Complex version.
            *   These routines perform a sum reduction of a dense matrix `H` (local portion) along process columns (hardcoded, `isRow` parameter is not used for direction). The result is accumulated on process row 0 of each column.

# Important Variables/Constants

*   `dp`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard pspBLAS precision and complex number parameters.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros.
*   **COMMON Blocks**: These are crucial for storing the pspBLAS parallel environment configuration.
    *   `/psp_grid2D/`: Contains `psp_mpi_comm_world`, `psp_mpi_size`, `psp_nprow`, `psp_npcol`, `psp_bs_def_row`, `psp_bs_def_col`, `psp_update_rank` (block size for SUMMA updates), `psp_bs_num`, `psp_icontxt` (BLACS context), `psp_mpi_comm_cart` (Cartesian communicator), `psp_mpi_comm_row` (row communicator), `psp_mpi_comm_col` (column communicator).
    *   `/psp_grid3D/`: Contains `psp_nplay` (number of layers in 3D grid) and `psp_mpi_comm_lay` (layer communicator). These are extensions to the `/psp_grid2D/` common block declarations if a 3D grid is used.

# Usage Examples

1.  **Grid Initialization**: This is a fundamental first step in any pspBLAS parallel application.
    ```fortran
    USE pspMPI
    IMPLICIT NONE
    INTEGER :: world_comm, num_processes, num_p_rows, blacs_ctxt_handle

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, world_comm, ierr) ! Or use existing app communicator
    CALL MPI_COMM_SIZE(world_comm, num_processes, ierr)

    num_p_rows = 2 ! Example: try to create a PxQ grid with P=2 rows
    ! Default block sizes for matrix distribution
    INTEGER, PARAMETER :: default_block_size_r = 64
    INTEGER, PARAMETER :: default_block_size_c = 64

    CALL psp_gridinit_2D(world_comm, num_processes, num_p_rows, 'R', &
                         default_block_size_r, default_block_size_c, blacs_ctxt_handle)
    ! The BLACS context and grid info are now stored in /psp_grid2D/ common block
    ! and blacs_ctxt_handle will also hold the created context.
    ! Other pspBLAS routines can now use this initialized grid.
    ! ...
    CALL BLACS_GRIDEXIT(blacs_ctxt_handle)
    CALL MPI_FINALIZE(ierr)
    ```

2.  **Custom Reductions**: These are typically used internally by parallel pspBLAS algorithms (like SUMMA in Level 3 operations) to combine partial results. For example, `psp_MPI_REDUCE_spm_packed` might be called by a parallel sparse matrix multiplication routine. End-users would not usually call these reduction functions directly.

# Dependencies and Interactions

*   **`pspUtility`**: The `psp_MPI_REDUCE_dspm_packed` and `psp_MPI_REDUCE_zspm_packed` routines call `psp_sst_sum_spmspm` (from `pspUtility` or `psp_spBLAS_Level3`) to sum sparse matrix data after receiving it.
*   **MPI Library**: This module is heavily dependent on MPI. It uses the `mpi` module (if `HAVE_MPI` is defined) and makes numerous direct MPI calls (e.g., `MPI_Cart_create`, `MPI_Cart_Sub`, `MPI_ISEND`, `MPI_IRECV`, `MPI_WAIT`, `MPI_PACK_SIZE`, `MPI_PACK`, `MPI_UNPACK`, `MPI_ADDRESS`, `MPI_Type_struct`, `MPI_Type_commit`).
*   **BLACS Library**: Routines `blacs_gridinfo` and `blacs_gridinit` are used for setting up and querying the process grid. The BLACS context handle (`psp_icontxt`) is a key piece of information stored and used by pspBLAS.
*   This `pspMPI` module is foundational for all parallel operations within the pspBLAS library. Higher-level pspBLAS routines (especially in Level 3) rely on the grid and communicators established here and utilize the custom reduction operations for their parallel algorithms.
*   The "TODO: bug inside" comment for the `_struct` versions of sparse matrix reduction suggests that the `_packed` versions are the currently stable methods for this task.
