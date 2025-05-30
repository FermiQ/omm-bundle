# Overview

The `pspUtility.F90` file defines the `pspUtility` module. This module serves as a central umbrella or aggregator for a variety of utility routines, data structures, and parameters that are used throughout the pspBLAS library. Instead of implementing functionalities directly, it leverages other, more specialized modules by using them and thus re-exporting their public entities. This provides a convenient way for other parts of the pspBLAS library to access a wide range of tools through a single `USE pspUtility` statement.

# Key Components

*   **Module `pspUtility`**:
    *   The primary purpose of this module is to consolidate various utility-focused modules within the pspBLAS ecosystem. It achieves this by `USE`ing other modules such as `pspVariable`, `pspBasicTool`, `pspListTool`, and `psp_spBLAS`.
    *   It does not define any new procedures, variables, or types itself.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined at compile time, it includes a `config.h` file. This file is typically used to manage project-specific build configurations, preprocessor definitions, and other settings that might influence how the library is compiled or operates.

# Usage Examples

A developer working on or using the pspBLAS library would include `USE pspUtility` in their Fortran code to gain access to the collective public entities (subroutines, functions, derived types, parameters, etc.) defined within the modules that `pspUtility` itself uses.

```fortran
USE pspUtility

! After this statement, all public entities from the following modules
! become accessible:
! - pspVariable (e.g., TYPE(psp_matrix_spm), linked list types like dList, zList)
! - pspBasicTool (e.g., psp_coo2csc, psp_den2sp_m, psp_idx_glb2loc)
! - pspListTool (e.g., psp_spm2list, psp_list2spm, psp_list_combine_listList)
! - psp_spBLAS (e.g., core sparse BLAS operations or further utilities)

! Conceptual Example (assuming psp_coo2csc is made public by pspBasicTool):
! TYPE(psp_matrix_spm) :: my_matrix
! ! ... initialize my_matrix, potentially in COO format ...
! IF (my_matrix%str_type == 'coo') THEN
!   CALL psp_coo2csc(my_matrix) ! Convert to CSC format
! END IF
```

# Dependencies and Interactions

The `pspUtility` module acts as a high-level aggregator and therefore has dependencies on the following modules:

*   **`pspVariable`**: This module is `USE`d, indicating that `pspUtility` provides access to fundamental data structure definitions (such as `TYPE(psp_matrix_spm)` for sparse matrices, various linked list types like `dList`, `zList`, `dNodeData`, `zNodeData`, and `dListPtrArray`), as well as associated constants and basic manipulation routines for these types.
*   **`pspBasicTool`**: This module is `USE`d, meaning `pspUtility` offers access to a range of basic utility functions. These typically include routines for converting between dense and sparse matrix representations, changing sparse matrix formats (e.g., COO to CSC), copying matrix portions, calculating global-to-local indices for distributed arrays, and processing operator flags (like transpose options).
*   **`pspListTool`**: This module is `USE`d, providing a collection of tools for creating, managing, and converting linked lists that are often employed in the dynamic construction of sparse matrices (e.g., during sparse matrix multiplication).
*   **`psp_spBLAS`**: This module is `USE`d. It likely contains the interfaces or implementations of core sparse BLAS operations (Level 1, 2, 3 kernels) or another layer of specialized utility routines that build upon the more basic tools.

By consolidating these modules, `pspUtility` simplifies the development within the pspBLAS library by reducing the number of `USE` statements needed in higher-level computational modules. These higher-level modules can then readily access a comprehensive suite of tools for data structure manipulation, format conversion, and other foundational operations.
