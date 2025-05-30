# Overview

The `pspBasicTool.F90` file defines the `pspBasicTool` module, which serves as a foundational component of the pspBLAS library. It aggregates a collection of essential utility subroutines and functions that support various low-level operations on sparse and dense matrices, as well as other general-purpose tasks needed by the library.

These tools include routines for:
*   Initializing sparse matrix structures (setting them to represent zero matrices).
*   Converting matrices between dense and sparse representations.
*   Converting sparse matrices between different storage formats, primarily Coordinate (COO) and Compressed Sparse Column (CSC).
*   Copying segments of sparse matrices.
*   Calculating global-to-local index mappings for distributed data.
*   Processing character flags for matrix operations (e.g., transpose, no transpose).
*   Initializing the random number generator.

Many of these utilities are provided for both real and complex data types through generic interfaces.

# Key Components

*   **Module `pspBasicTool`**:
    *   **Public Interfaces and Subroutines**:
        *   `psp_spm_zeros(m_name, M, N, fmt, isReal)`: Initializes a `TYPE(psp_matrix_spm)` variable `m_name` to represent a zero sparse matrix. It sets global dimensions (M, N), the specified format `fmt` ('coo' or 'csc'), data type (real/complex via `isReal`), local dimensions (using `numroc` based on `/psp_grid2D/` info), and initializes internal arrays (e.g., `col_ptr` for CSC to represent no non-zeros).
        *   `psp_den2sp_m`: Interface for converting a dense matrix (Fortran array) to a `TYPE(psp_matrix_spm)` sparse matrix.
            *   `psp_den2sp_dm(denMat, desc, spMat, fmt, thre)`: Real version.
            *   `psp_den2sp_zm(denMat, desc, spMat, fmt, thre)`: Complex version.
            *   These routines take the dense matrix `denMat`, its BLACS descriptor `desc`, an output `spMat` of `TYPE(psp_matrix_spm)`, the target sparse format `fmt` ('coo' or 'csc'), and an optional threshold `thre` to prune small values. They first generate COO data from `denMat` based on `thre`, then register this data into `spMat` (using `psp_register_spm`), and finally convert `spMat` to the desired `fmt` if it's 'csc'.
        *   `psp_sp2den_m`: Interface for converting a `TYPE(psp_matrix_spm)` sparse matrix to a dense Fortran array.
            *   `psp_sp2den_dm(spMat, denMat, desc)`: Real version.
            *   `psp_sp2den_zm(spMat, denMat, desc)`: Complex version.
            *   These routines allocate the `denMat` based on `spMat%loc_dim1` and `spMat%loc_dim2`, copy the BLACS descriptor `spMat%desc` to `desc`, and then populate `denMat` from the sparse data.
        *   `psp_sst_den2sp_m`: ("Subroutine Subprogram Template") Interface for converting a dense Fortran array directly into raw sparse format arrays (values, row indices, column indices/pointers).
            *   `psp_sst_den2sp_dm(denMat, idx1, idx2, val, fmt, thre)`: Real version.
            *   `psp_sst_den2sp_zm(denMat, idx1, idx2, val, fmt, thre)`: Complex version.
            *   `idx1` gets row indices, `val` gets values. `idx2` gets column indices if `fmt` is 'coo', or becomes column pointers if `fmt` is 'csc' (via internal call to `psp_sst_coo2csc`).
        *   `psp_sst_sp2den_m`: Interface for converting raw sparse format arrays (values, row indices, column indices/pointers) to a dense Fortran array.
            *   `psp_sst_sp2den_dm(m,n,idx1,idx2,val,fmt,denMat)`: Real version.
            *   `psp_sst_sp2den_zm(m,n,idx1,idx2,val,fmt,denMat)`: Complex version.
        *   `psp_sst_coo2csc(m, n, nnz, idx)`: Converts COO column indices stored in `idx` (representing `col_ind`) into CSC column pointers (CSR `row_ptr` equivalent) in-place within `idx`.
        *   `psp_sst_csc2coo(m, n, nnz, col_ptr)`: Converts CSC column pointers stored in `col_ptr` into COO column indices in-place within `col_ptr`.
        *   `psp_sst_fmtCnvt(m, n, nnz, idx, fmt1, fmt2)`: A generic SST routine to convert `idx` between COO column indices and CSC column pointers based on `fmt1` and `fmt2`.
        *   `psp_coo2csc(spMat)`: Converts a `TYPE(psp_matrix_spm)` variable `spMat` from COO format to CSC format in-place. Modifies `spMat%str_type` and converts `spMat%col_ind` to `spMat%col_ptr`.
        *   `psp_csc2coo(spMat)`: Converts a `TYPE(psp_matrix_spm)` variable `spMat` from CSC format to COO format in-place. Modifies `spMat%str_type` and converts `spMat%col_ptr` to `spMat%col_ind`.
        *   `psp_copy_spm2st`: Interface for copying a submatrix from a `TYPE(psp_matrix_spm)` (assumed to be in CSC format) to separate CSC arrays.
            *   `psp_copy_dspm2st(M,N,A,IA,JA,B_idx1,B_idx2,B_val,B_dim1,B_dim2,IB,JB,beta)`: Real version. Copies `A(IA:IA+M-1, JA:JA+N-1)` to `B` arrays.
            *   `psp_copy_zspm2st(...)`: Complex version.
        *   `psp_idx_glb2loc(glb, bs, npproc, loc)`: Calculates the local index `loc` for a 1D global index `glb`, given a block size `bs` and the number of processes `npproc` in that dimension of a block-cyclically distributed array.
        *   `psp_process_opM`: Interface for converting character operation flags ('N', 'T', 'C') into internal representations.
            *   `psp_process_lopM(opM, trM)`: Output `trM` is `LOGICAL` (for real matrices, 'T'/'C' map to .TRUE.).
            *   `psp_process_iopM(opM, tcM)`: Output `tcM` is `INTEGER` (for complex matrices: 0 for 'N', 1 for 'C', 2 for 'T').
        *   `init_random_seed()`: Initializes the seed for Fortran's intrinsic `RANDOM_SEED` subroutine. It attempts to read from `/dev/urandom`. If unsuccessful, it uses a combination of system clock time and MPI process rank (if available) to generate a seed, employing a simple Linear Congruential Generator (LCG) as part of this fallback.
    *   **Subroutine `die(message)`** (private): Standard error handling routine that writes to "MatrixSwitch.log" and stops.

# Important Variables/Constants

*   `dp`: Integer parameter for double precision kind.
*   `cmplx_1, cmplx_i, cmplx_0`: Standard complex constants of kind `dp`.
*   `HAVE_CONFIG_H`, `HAVE_MPI`: Preprocessor macros affecting conditional compilation.
*   **COMMON Block `/psp_grid2D/`**: Used by `psp_spm_zeros` and `psp_den2sp_dm`/`zm` to get BLACS grid information (`psp_icontxt`, `psp_nprow`, `psp_npcol`, `psp_bs_def_row`, `psp_bs_def_col`) for setting up `psp_matrix_spm` descriptors and local dimensions.
*   `numroc`: Declared as an `EXTERNAL` function (from ScaLAPACK), used by `psp_spm_zeros` for calculating local dimensions.

# Usage Examples

These routines are primarily utilities called by other higher-level pspBLAS routines.
```fortran
USE pspBasicTool
USE pspVariable ! For TYPE(psp_matrix_spm)
IMPLICIT NONE

TYPE(psp_matrix_spm) :: sparse_A
REAL(dp), DIMENSION(100,50) :: dense_A_data
INTEGER :: blacs_desc(9), M_global, N_global

! Example: Initialize a sparse matrix structure to represent a 1000x500 zero matrix in CSC format
! (Assumes BLACS grid has been initialized and info is in /psp_grid2D/)
CALL psp_spm_zeros(sparse_A, 1000, 500, 'csc', .TRUE.)

! Example: Convert a dense matrix to a sparse matrix structure
! (Requires dense_A_data, blacs_desc to be set up)
! CALL psp_den2sp_m(dense_A_data, blacs_desc, sparse_A, 'csc', 1.0E-8_dp)

! Example: Convert a psp_matrix_spm from COO to CSC format
! IF (sparse_A%str_type == 'coo') THEN
!   CALL psp_coo2csc(sparse_A)
! END IF

! Example: Initialize random seed
CALL init_random_seed()
```

# Dependencies and Interactions

*   **`pspVariable`**: This module is `USE`d, primarily for the `TYPE(psp_matrix_spm)` definition which is manipulated by many routines in `pspBasicTool`.
*   **MPI Library**: Included if `HAVE_MPI` is defined. `init_random_seed` uses `MPI_COMM_RANK` (via `psp_mpi_comm_world` from `/psp_grid2D/`). Several routines rely on grid parameters from `/psp_grid2D/` which is MPI/BLACS dependent.
*   **BLACS Library**: `blacs_gridinfo` is called by `psp_spm_zeros`. `descinit` is also called by `psp_spm_zeros` to initialize the BLACS descriptor within the `psp_matrix_spm` type. The `psp_den2sp_m` routines take a BLACS descriptor as input.
*   **ISO_FORTRAN_ENV**: Used in `init_random_seed` for `int64` type.
*   This module provides many of the fundamental building blocks for handling sparse matrix data structures, format conversions, and other basic tasks required by the computational routines in pspBLAS Level 1, 2, and 3 modules. The `_sst_` routines suggest a layer of abstraction for raw data manipulation, while other routines work directly with the `psp_matrix_spm` type.
