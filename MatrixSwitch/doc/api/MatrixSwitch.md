# Overview

The file `MatrixSwitch.F90` defines the main module `MatrixSwitch`. This module provides a high-level interface for various matrix operations, including allocation, deallocation, copy, conversion, multiplication, addition, trace, scaling, and setting/getting elements. It acts as a central hub, utilizing other specialized `MatrixSwitch_*` modules for specific functionalities. It also includes conditional compilation for MPI, ScaLAPACK, DBCSR, and pspBLAS libraries to enable parallel and specialized operations.

# Key Components

*   **Module `MatrixSwitch`**: The main module.
    *   **Interfaces**:
        *   `mm_multiply`: Matrix-matrix multiplication (C := alpha\*op(A)\*op(B) + beta\*C). Procedures: `mm_dmultiply`, `mm_zmultiply`.
        *   `m_add`: Matrix addition (C := alpha\*op(A) + beta\*C). Procedures: `m_dadd`, `m_zadd`.
        *   `m_trace`: Matrix trace (alpha := tr(A)). Procedures: `m_dtrace`, `m_ztrace`.
        *   `mm_trace`: Matrix product trace (alpha := tr(A^H\*B)). Procedures: `mm_dtrace`, `mm_ztrace`.
        *   `m_scale`: Scale matrix (C := beta\*C). Procedures: `m_dscale`, `m_zscale`.
        *   `m_set`: Set matrix elements (C_ij := {alpha (i /= j), beta (i == j)}). Procedures: `m_dset`, `m_zset`.
        *   `m_set_element`: Set a single matrix element (C_ij := alpha + beta\*C_ij). Procedures: `m_dset_element`, `m_zset_element`, `m_dset_block`.
        *   `m_get_element`: Get a single matrix element (alpha := C_ij). Procedures: `m_dget_element`, `m_zget_element`, `m_dget_block`.
        *   `m_allocate`: Allocate matrix. Procedures: `m_allocate_elements`, `m_allocate_blocks`.
    *   **Public Subroutines/Functions (from module itself and other modules)**:
        *   `matrix`: (Type definition, likely from `MatrixSwitch_ops` or a wrapper) - Represents the matrix data structure.
        *   `m_allocate_elements(m_name, i, j, label)`: Allocates a matrix by elements, initializing basic information and arrays.
        *   `m_allocate_blocks(m_name, row_sizes, col_sizes, label)`: Allocates a blocked matrix using block-cycling distribution.
        *   `m_deallocate(m_name)`: Deallocates arrays within a `TYPE(MATRIX)` variable.
        *   `m_copy(m_name, A, label, threshold, threshold_is_soft)`: Copies matrix A to `m_name`, allowing format conversion and thresholding.
        *   `m_convert(m_name, label, threshold, threshold_is_soft)`: Converts a matrix to a new storage format in-place.
        *   `mm_dmultiply(A, opA, B, opB, C, alpha, beta, label)`: Real matrix-matrix multiplication.
        *   `mm_zmultiply(A, opA, B, opB, C, alpha, beta, label)`: Complex matrix-matrix multiplication.
        *   `m_dadd(A, opA, C, alpha, beta, label)`: Real matrix addition.
        *   `m_zadd(A, opA, C, alpha, beta, label)`: Complex matrix addition.
        *   `m_dtrace(A, alpha, label)`: Real matrix trace.
        *   `m_ztrace(A, alpha, label)`: Complex matrix trace.
        *   `mm_dtrace(A, B, alpha, label)`: Real matrix product trace.
        *   `mm_ztrace(A, B, alpha, label)`: Complex matrix product trace.
        *   `m_dscale(C, beta, label)`: Real matrix scaling.
        *   `m_zscale(C, beta, label)`: Complex matrix scaling.
        *   `m_dset(C, seC, alpha, beta, label)`: Set real matrix elements based on diagonal/off-diagonal.
        *   `m_zset(C, seC, alpha, beta, label)`: Set complex matrix elements based on diagonal/off-diagonal.
        *   `m_dset_element(C, i, j, alpha, beta, label)`: Set a specific real matrix element.
        *   `m_zset_element(C, i, j, alpha, beta, label)`: Set a specific complex matrix element.
        *   `m_dset_block(C, i, j, block_data, beta)`: Set a block of a real matrix (DBCSR specific).
        *   `m_dget_element(C, i, j, alpha, label)`: Get a specific real matrix element.
        *   `m_zget_element(C, i, j, alpha, label)`: Get a specific complex matrix element.
        *   `m_dget_block(C, i, j, block_data, found)`: Get a block of a real matrix (DBCSR specific).
        *   `m_register_sden`: (Likely from `MatrixSwitch_m_register`) Registers a serial dense matrix.
        *   `m_register_pdbc` (MPI): (Likely from `MatrixSwitch_m_register`) Registers a parallel dense block cyclic matrix.
        *   `ms_lap_icontxt` (MPI): BLACS context handle for ScaLAPACK.
        *   `ms_scalapack_setup(...)` (MPI & ScaLAPACK): Sets up ScaLAPACK environment.
        *   `ms_dbcsr_setup(...)` (MPI & DBCSR): Sets up DBCSR environment.
        *   `ms_dbcsr_finalize()` (MPI & DBCSR): Finalizes DBCSR environment.
        *   `m_copy_external_pdbcpcoo` (pspBLAS): Copies external data to parallel COO format.
        *   `m_copy_external_pdbcpcsc` (pspBLAS): Copies external data to parallel CSC format.
        *   `m_register_pcoo` (pspBLAS): (Likely from `MatrixSwitch_m_register`) Registers a parallel COO matrix.
        *   `m_register_pcsc` (pspBLAS): (Likely from `MatrixSwitch_m_register`) Registers a parallel CSC matrix.

# Important Variables/Constants

*   Preprocessor Macros: `HAVE_CONFIG_H`, `HAVE_MPI`, `HAVE_DBCSR`, `HAVE_SCALAPACK`, `HAVE_PSPBLAS`, `CONV`. These control conditional compilation for various features and libraries.
*   Module-level variables for ScaLAPACK setup (within `ms_scalapack_setup` and used by `ms_scalapack_allocate`):
    *   `ms_mpi_comm`: MPI communicator.
    *   `ms_mpi_size`: Size of the MPI communicator.
    *   `ms_mpi_rank`: Rank in the MPI communicator.
    *   `ms_lap_nprow`: Number of process rows in BLACS grid.
    *   `ms_lap_npcol`: Number of process columns in BLACS grid.
    *   `ms_lap_order`: Process grid ordering ('C' or 'R').
    *   `ms_lap_bs_def`: Default block size for ScaLAPACK.
    *   `ms_lap_bs_list`: List of specific block sizes for matrix dimensions.
    *   `ms_lap_icontxt`: BLACS context handle.
*   Module-level variables for DBCSR setup:
    *   `ms_dbcsr_init`: Logical flag, true if DBCSR is initialized.
    *   `ms_dbcsr_group`: MPI communicator for the DBCSR 2D grid.

# Usage Examples

Usage examples are not directly available in the comments of this specific file. However, the detailed comments for each interface and subroutine describe the parameters and their purpose, which guides their usage. For example, to perform a matrix multiplication `C = A * B`:
```fortran
! Assuming A, B, C are of TYPE(matrix) and already allocated and initialized
CALL mm_multiply(A, 'N', B, 'N', C, 1.0_dp, 0.0_dp)
```

# Dependencies and Interactions

*   **Uses other `MatrixSwitch_*` modules**:
    *   `MatrixSwitch_ops`: Likely contains the `matrix` type definition and basic operations or constants like `dp`, `cmplx_0`, `die`.
    *   `MatrixSwitch_mm_multiply`: Provides specific implementations for matrix multiplication.
    *   `MatrixSwitch_m_add`: Provides specific implementations for matrix addition.
    *   `MatrixSwitch_m_set`: Provides specific implementations for setting matrix values.
    *   `MatrixSwitch_m_copy`: Provides specific implementations for copying matrices.
    *   `MatrixSwitch_m_register`: Provides procedures for registering externally allocated matrices.
*   **External Libraries (conditional)**:
    *   `pspBLAS` (if `HAVE_PSPBLAS` is defined): Used for sparse matrix operations.
    *   `dbcsr_api` (if `HAVE_MPI` and `HAVE_DBCSR` are defined): Used for distributed block compressed sparse row (DBCSR) matrix operations.
    *   MPI (`mpif.h`): Used for parallel processing if `HAVE_MPI` is defined. ScaLAPACK and DBCSR typically depend on MPI.
    *   LAPACK (implicitly via `dgemm`, `zgemm`, `ddot` etc.): Standard library for linear algebra operations, often used as a backend for serial operations.
    *   ScaLAPACK (if `HAVE_MPI` and `HAVE_SCALAPACK` are defined): Used for parallel dense linear algebra operations (e.g., `pdgemm`, `pzgeadd`, `pdlaset`). BLACS (`blacs_gridinit`, `blacs_gridinfo`) is its communication layer.
*   **Includes**:
    *   `config.h` (if `HAVE_CONFIG_H` is defined): Project-specific configuration.
    *   `mpif.h` (in some subroutines if `HAVE_MPI` is defined): MPI Fortran interface.

This module orchestrates calls to various specialized subroutines based on matrix types (real/complex, dense/sparse, serial/parallel) and requested operations. The `die` subroutine (likely from `MatrixSwitch_ops`) is used for error handling.
