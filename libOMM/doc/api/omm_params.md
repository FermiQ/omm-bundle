# Overview

The `omm_params` module in the libOMM library serves as a centralized repository for defining fundamental parameters, mathematical constants, and global (module-level) state variables that are utilized across various parts of the library. This includes specifying the precision for floating-point calculations, defining common constants like Pi and complex numbers, and holding global variables related to the execution environment, such as MPI context and logging file units.

# Key Components

*   **Module `omm_params`**:
    *   **Parameters (Compile-time constants)**:
        *   `dp`: An integer parameter derived from `selected_real_kind(15,300)`. It defines the kind for double-precision floating-point numbers, ensuring consistent precision throughout the libOMM library.
        *   `Pi`: A `REAL(dp)` parameter holding the value of Pi to a high degree of precision.
        *   `cmplx_1`: A `COMPLEX(dp)` parameter representing the complex number (1.0, 0.0).
        *   `cmplx_i`: A `COMPLEX(dp)` parameter representing the imaginary unit (0.0, 1.0).
        *   `cmplx_0`: A `COMPLEX(dp)` parameter representing the complex number (0.0, 0.0).
        *   `i64`: An integer parameter derived from `selected_int_kind(18)`, defining a kind for 64-bit integers suitable for large integer counts or indices.
    *   **Module-Level Variables (Global State with `SAVE` attribute)**:
        *   `ms_scalapack_running`: A logical variable, initialized to `.FALSE.`. This flag is likely used by the MatrixSwitch interface or libOMM itself to track whether ScaLAPACK has been initialized and is currently operational.
        *   `log_unit`: An integer variable. This is intended to store the Fortran I/O unit number that is assigned to the libOMM log file (e.g., "libOMM.log"). Its value is not initialized within this module but is expected to be set during an initial setup phase of the library or application.
        *   `mpi_size`: An integer variable. It is designed to hold the total number of processes in the MPI communicator being used.
        *   `mpi_rank`: An integer variable. It is designed to hold the rank of the current MPI process within its communicator.

# Important Variables/Constants

*   `dp`: This is arguably the most critical parameter, as it standardizes double-precision calculations across the entire libOMM and associated MatrixSwitch operations.
*   Mathematical Constants (`Pi`, `cmplx_1`, `cmplx_i`, `cmplx_0`): Provide readily available high-precision values for common mathematical needs.
*   `i64`: Offers a portable way to declare 64-bit integers when necessary.
*   Global State Variables:
    *   `ms_scalapack_running`: Acts as a global flag for ScaLAPACK's status.
    *   `log_unit`: Essential for consistent logging output. Its proper initialization elsewhere is critical.
    *   `mpi_size` and `mpi_rank`: Fundamental for any MPI-parallelized code sections to understand their execution context. Like `log_unit`, these are expected to be initialized externally (e.g., during MPI setup and potentially passed into a libOMM initialization routine).

# Usage Examples

The `omm_params` module itself does not contain executable procedures beyond declarations. Its parameters and variables are intended for use in other modules within the libOMM library.

```fortran
! Example of using omm_params in another Fortran module or subroutine
MODULE another_libOMM_routine
  USE omm_params ! Gains access to dp, Pi, log_unit, mpi_rank, etc.

  IMPLICIT NONE

  SUBROUTINE do_calculation(input_value)
    REAL(dp), INTENT(IN) :: input_value
    REAL(dp) :: result

    result = input_value * Pi / 2.0_dp

    IF (mpi_rank == 0) THEN
      WRITE(log_unit, *) "Calculation result on rank 0: ", result
    END IF
  END SUBROUTINE do_calculation

END MODULE another_libOMM_routine
```

# Dependencies and Interactions

*   This module is self-contained in its definitions and does not `USE` other libOMM or MatrixSwitch modules to define its parameters. The choice of `dp` is, however, designed to be consistent with the precision used in `MatrixSwitch_ops`.
*   It serves as a foundational module for the rest of libOMM. Many other modules (e.g., `omm_ops`, `omm`, `omm_callback`) will typically include a `USE omm_params` statement to access these shared definitions.
*   The `SAVE` attribute on variables like `log_unit`, `mpi_size`, `mpi_rank`, and `ms_scalapack_running` ensures that their values persist across different procedure calls within the library, effectively establishing them as global state variables. The responsibility for initializing these global state variables correctly lies outside this module, typically in higher-level setup or initialization routines within the application or the libOMM library itself.
