# Overview

The module `MatrixSwitch_ops` serves as a foundational part of the MatrixSwitch library. It defines the core `matrix` derived type, which is used throughout the library to represent matrices in various formats. It also provides fundamental parameters (like double precision kind, complex constants), global variables for MPI/BLACS and DBCSR setup (if enabled), and utility subroutines for parameter processing and error handling.

# Key Components

*   **Module `MatrixSwitch_ops`**:
    *   **Type `matrix`**: The central derived data type for all MatrixSwitch operations. It includes:
        *   `str_type`: Character(3) label for storage format (e.g., "den", "csc", "dbc").
        *   `is_initialized`, `is_serial`, `is_real`, `is_square`, `is_sparse`: Logical flags describing matrix properties.
        *   `*_is_allocated`: Logical flags indicating if data/index arrays are allocated by MatrixSwitch or are pointers to external data. Crucial for memory management, especially with registered matrices.
        *   `dim1`, `dim2`: Integer dimensions of the matrix.
        *   `iaux1`, `iaux2`, `iaux3`, `iaux4`: Pointers to integer arrays for auxiliary data (e.g., sparse indices, BLACS descriptors for parallel dense, local dimensions for parallel sparse).
        *   `dval`, `zval`: Pointers to real(dp) or complex(dp) 2D arrays for matrix element values.
        *   `spm` (if `HAVE_PSPBLAS`): `psp_matrix_spm` type for pspBLAS internal storage for parallel sparse matrices.
        *   `dbcsr_dist`, `dbcsr_mat` (if `HAVE_MPI` and `HAVE_DBCSR`): DBCSR specific types for distribution and matrix objects for block sparse parallel matrices.
    *   **Parameters**:
        *   `dp`: Integer parameter for double precision kind (`selected_real_kind(15,300)`).
        *   `cmplx_1`, `cmplx_i`, `cmplx_0`: Complex(dp) constants for 1, i (imaginary unit), and 0.
    *   **Global Variables (MPI/BLACS/DBCSR related, conditionally compiled if `HAVE_MPI` is defined)**:
        *   `ms_lap_order`: Character(1), ordering of BLACS grid ('R' or 'C').
        *   `ms_mpi_comm`: Integer, MPI communicator used by MatrixSwitch.
        *   `ms_mpi_size`: Integer, number of processes in `ms_mpi_comm`.
        *   `ms_mpi_rank`: Integer, rank of the current process in `ms_mpi_comm`.
        *   `ms_lap_nprow`, `ms_lap_npcol`: Integers, number of process rows/columns in BLACS grid.
        *   `ms_lap_bs_def`: Integer, default block size for ScaLAPACK.
        *   `ms_lap_bs_num`: Integer, number of custom block size exceptions.
        *   `ms_lap_icontxt`: Integer, BLACS context handle.
        *   `ms_lap_bs_list(:,:)`: Allocatable integer array for (dimension, block_size) pairs.
        *   `ms_dbcsr_init` (if `HAVE_DBCSR`): Logical, true if DBCSR is initialized.
        *   `ms_dbcsr_group` (if `HAVE_DBCSR`): Integer, MPI communicator for DBCSR grid.
    *   **Interface `process_opM`**: Converts character operation codes ('N' for None, 'T' for Transpose, 'C' for Conjugate transpose) to internal logical or integer flags.
        *   `process_lopM(opM, trM)`: For real matrices (output `trM` is logical: .true. for T/C, .false. for N).
        *   `process_iopM(opM, tcM)`: For complex matrices (output `tcM` is integer: 0 for N, 1 for C, 2 for T).
    *   **Subroutine `process_seM(seM, luM)`**: Converts character codes for matrix part selection ('L' for Lower, 'U' for Upper, other for Full matrix) into an integer flag `luM` (0 for full, 1 for upper, 2 for lower). Used by `m_set` operations.
    *   **Subroutine `die(message)`**: Standardized error handling routine. Prints an optional message and MPI rank (if applicable) to "MatrixSwitch.err" and stops execution.

# Important Variables/Constants

*   `dp`: Defines the double precision kind used for all real and complex matrix data throughout MatrixSwitch.
*   `matrix` type: The core data structure. Its flexible pointer-based design and descriptive flags allow it to handle a wide variety of matrix types and storage schemes, including those managed externally.
*   MPI/BLACS/DBCSR global variables: These store the necessary context for parallel operations, initialized by setup routines (e.g., `ms_scalapack_setup`) and used by various computational kernels.

# Usage Examples

This module primarily provides definitions and utilities that are used internally by other `MatrixSwitch_*` modules and the main `MatrixSwitch` module. Direct user interaction is generally not expected.
The `die` subroutine is a crucial part of error handling:
```fortran
IF (size(matrix_A%dval, 1) /= expected_rows) THEN
  CALL die("Matrix A has incorrect number of rows.")
END IF
```

# Dependencies and Interactions

*   **pspBLAS** (if `HAVE_PSPBLAS` is defined): The `matrix` type includes `psp_matrix_spm` from the `pspBLAS` module.
*   **MPI** (if `HAVE_MPI` is defined): The `matrix` type and global variables include MPI related handles and communicators (e.g. `ms_dbcsr_group` uses `mpi_comm_null` from `mpi` module).
*   **DBCSR** (if `HAVE_MPI` and `HAVE_DBCSR` are defined): The `matrix` type includes `dbcsr_distribution_type` and `dbcsr_type` from the `dbcsr_api` module.
*   This module is foundational and is `USE`d by virtually all other modules within the MatrixSwitch library due to its definition of the `matrix` type, precision kinds, and common utility routines like `process_opM`, `process_seM`, and `die`.
*   The global variables for parallel environments are intended to be set by specific setup routines (e.g., `ms_scalapack_setup`, `ms_dbcsr_setup` found in `MatrixSwitch.F90`) and then utilized by computational routines across different modules when dealing with parallel matrix operations.
