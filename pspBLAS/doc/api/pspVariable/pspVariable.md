# Overview

The `pspVariable.F90` file defines the `pspVariable` module. This module acts as a central hub or an aggregator within the pspBLAS library, primarily focused on providing convenient access to various custom data types and variable-related definitions. Instead of defining these types itself, it `USE`s other specialized modules (`pspMat` and `pspListType`) and thereby re-exports their public entities. This allows other parts of the pspBLAS library to access key data structures, such as the sparse matrix type and linked list definitions, through a single, consolidated `USE pspVariable` statement.

# Key Components

*   **Module `pspVariable`**:
    *   The core function of this module is to serve as an aggregator for modules that define important data types used in pspBLAS.
    *   It achieves this by containing `USE` statements for `pspMat` and `pspListType`.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined during compilation, it includes a `config.h` file. This file is typically used for managing project-wide build configurations and preprocessor definitions that might affect how the library is compiled or behaves.

# Usage Examples

The `pspVariable` module is intended to be `USE`d by other modules within the pspBLAS library that need to declare or operate on the custom data types it provides access to.

```fortran
USE pspVariable  ! This statement grants access to entities from pspMat and pspListType

IMPLICIT NONE

! Declare a sparse matrix using the type from pspMat (via pspVariable)
TYPE(psp_matrix_spm) :: my_sparse_matrix_A, my_sparse_matrix_B

! Declare a real-valued linked list and a node data item (from pspListType via pspVariable)
TYPE(dList), POINTER :: head_of_my_real_list => NULL()
TYPE(dNodeData)      :: real_list_item

! Initialize the sparse matrix (routine from pspMat via pspVariable)
! CALL psp_register_spm(my_sparse_matrix_A, ...)

! Create a linked list (routine from pspListType via pspVariable)
! real_list_item%row_ind = 1
! real_list_item%col_ind = 1
! real_list_item%val     = 10.0_dp
! CALL list_create(head_of_my_real_list, real_list_item)
```
By using `pspVariable`, a developer avoids needing to explicitly `USE pspMat` and `USE pspListType` separately when both sets of definitions are required.

# Dependencies and Interactions

The `pspVariable` module has direct dependencies on the following modules, whose public entities it effectively re-exports:

*   **`pspMat`**: This module is `USE`d by `pspVariable`. `pspMat` is responsible for defining the core sparse matrix data structure, `TYPE(psp_matrix_spm)`, and associated routines like `psp_register_spm` (for populating the structure) and `psp_deallocate_spm` (for freeing its components).
*   **`pspListType`**: This module is also `USE`d by `pspVariable`. `pspListType` provides generic interfaces for linked list operations (e.g., `list_create`, `list_insert`, `list_destroy`) and makes available the specialized linked list types (`dList`, `zList`), node data types (`dNodeData`, `zNodeData`), and array-of-list-pointer types (`dListPtrArray`, `zListPtrArray`).

The `pspVariable` module itself is likely to be a foundational module `USE`d by many other higher-level modules within the pspBLAS library, including utility modules (like `pspUtility`) and the modules implementing Level 1, 2, and 3 BLAS operations, as they all will need to work with the sparse matrix and linked list data structures. This aggregation simplifies module dependencies throughout the library.
