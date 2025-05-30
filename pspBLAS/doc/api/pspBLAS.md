# Overview

The `pspBLAS.F90` file defines the main `pspBLAS` module. This module serves as the primary entry point or umbrella module for the pspBLAS library. It consolidates various components of the library by using other pspBLAS submodules, making the library's functionalities accessible to the end-user through a single `USE pspBLAS` statement.

# Key Components

*   **Module `pspBLAS`**: The top-level module for the pspBLAS library. It does not define any procedures or variables itself but rather aggregates other core modules of the library.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined, it includes a `config.h` file, which is typically used for managing build configurations and preprocessor definitions across the project.

# Usage Examples

Usage examples are not explicitly provided within this file. A user intending to use the pspBLAS library would typically include `USE pspBLAS` in their Fortran code. This statement would then provide access to all public entities (subroutines, functions, types, parameters) defined within the modules that `pspBLAS` itself uses.

```fortran
USE pspBLAS

! After this statement, one can call routines from the sub-modules like:
! - pspLevel1 (e.g., for vector operations)
! - pspLevel2 (e.g., for matrix-vector operations)
! - pspLevel3 (e.g., for matrix-matrix operations)
! - pspUtility (e.g., for helper functions)
! - pspVariable (e.g., for matrix type definitions)
! - pspMPI (e.g., for MPI-related parallel functionalities)

! Example (conceptual):
! TYPE(psp_matrix_type) :: A, B, C
! CALL psp_some_level3_operation(A, B, C)
```

# Dependencies and Interactions

The `pspBLAS` module has dependencies on several other modules within the pspBLAS library:

*   **`pspVariable`**: This module is `USE`d, indicating it likely provides definitions for fundamental data structures, such as the specific `TYPE` for matrices (e.g., sparse matrix representations) and other variables used throughout pspBLAS.
*   **`pspUtility`**: This module is `USE`d, suggesting it contains various utility functions, helper routines, or common tools that support the core operations of the library.
*   **`pspLevel1`**: This module is `USE`d and is expected to provide Level 1 BLAS-like operations, which typically include vector operations like dot products, norms, scaling, and copy.
*   **`pspLevel2`**: This module is `USE`d and is expected to provide Level 2 BLAS-like operations, which generally involve matrix-vector operations such as matrix-vector multiplication or solving triangular systems for a single vector.
*   **`pspLevel3`**: This module is `USE`d and is expected to provide Level 3 BLAS-like operations, primarily focusing on matrix-matrix operations like matrix multiplication.
*   **`pspMPI`**: This module is `USE`d, indicating that pspBLAS has capabilities for parallel execution using the Message Passing Interface (MPI). This module would contain MPI-related initializations, communication routines, or parallel data distributions.

The `pspBLAS` module acts as a high-level aggregator, simplifying the inclusion of library components for the end-user. By using `pspBLAS`, the user gains access to the functionalities of all the aforementioned modules.
