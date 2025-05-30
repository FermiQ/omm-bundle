# Overview

The `pspListType.F90` file defines the `pspListType` module. This module serves as a crucial abstraction layer within the pspBLAS library for handling linked list operations in a type-generic manner. It achieves this by providing a set of generic interfaces for common linked list manipulations (such as creation, insertion, deletion, counting, etc.). These generic interfaces then dispatch to specialized implementations for lists containing real (`dList`) or complex (`zList`) data, which are imported from the `psp_dList` and `psp_zList` modules, respectively.

This approach allows higher-level pspBLAS routines to use a common set of function names for list operations, regardless of the underlying data type of the list elements, thereby simplifying code and improving maintainability.

Additionally, the `pspListType` module defines simple derived types, `dListPtrArray` and `zListPtrArray`. These types are designed to hold pointers to `dList` and `zList` respectively, and are likely used to create and manage arrays of linked lists. Such structures are particularly useful in sparse matrix algorithms, for example, when accumulating contributions to different columns (or rows) of a sparse matrix product into separate lists before final assembly.

# Key Components

*   **Module `pspListType`**:
    *   **Imports and Renaming**:
        *   From `psp_dList`: It imports the specialized `LINKED_LIST` type (renaming it to `dList`), the specialized `LIST_DATA` type (renaming it to `dNodeData`), and all the generic list operation procedures (e.g., `list_create` is renamed to `dlist_create`, `list_insert` to `dlist_insert`, and so on).
        *   From `psp_zList`: Similarly, it imports and renames the complex counterparts: `LINKED_LIST` becomes `zList`, `LIST_DATA` becomes `zNodeData`, and procedures like `list_create` become `zlist_create`.

    *   **Public Derived Types**:
        *   `dListPtrArray`: A simple type definition containing a single component: `ptr`, which is a `POINTER` to `dList`. This allows for creating arrays where each element can point to the head of a real-valued linked list.
        *   `zListPtrArray`: Similar to `dListPtrArray`, but its `ptr` component is a `POINTER` to `zList`, for arrays of complex-valued linked lists.

    *   **Public Generic Interfaces**: For each fundamental linked list operation, a generic interface is defined. This interface maps a common name (e.g., `list_create`) to the specific module procedures imported from `psp_dList` and `psp_zList` (e.g., `dlist_create` and `zlist_create`). The Fortran compiler then resolves calls to these generic interfaces to the correct type-specific version based on the arguments provided. The exposed generic interfaces include:
        *   `list_create(list, data)`
        *   `list_count(list)`
        *   `list_destroy(list)`
        *   `list_insert(elem, data)`
        *   `list_insert_head(list, data)`
        *   `list_next(elem)`
        *   `list_get_data(elem)`
        *   `list_put_data(elem, data)`
        *   `list_delete_element(list, elem)`

    *   **Public Type Aliases**:
        *   `dList`, `zList`: These are made public, providing convenient aliases for the specialized `LINKED_LIST` types for real and complex data.
        *   `dNodeData`, `zNodeData`: These are made public, providing aliases for the specialized `LIST_DATA` types that define the structure of data held in real and complex list nodes (typically including row index, column index, and value).

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined, it includes a `config.h` file, generally used for build-time configurations.
*   The module itself does not define local parameters like `dp` (double precision kind). It relies on these being consistently defined and used within the imported `psp_dList` and `psp_zList` modules (which, in turn, likely inherit them from a common source like `pspVariable` or `pspBasicTool`).

# Usage Examples

The `pspListType` module allows other parts of pspBLAS to work with linked lists in a more abstract, type-generic way.

```fortran
USE pspListType  ! Provides generic interfaces like list_create,
                 ! and types like dList, dNodeData, zList, zNodeData.
IMPLICIT NONE

TYPE(dList), POINTER :: real_elements_list => NULL()
TYPE(dNodeData)      :: real_data_item
TYPE(zList), POINTER :: complex_elements_list => NULL()
TYPE(zNodeData)      :: complex_data_item
INTEGER              :: count_real, count_complex

! Working with a list of real sparse matrix elements
real_data_item%row_ind = 10
real_data_item%col_ind = 20
real_data_item%val     = 3.14_dp
CALL list_create(real_elements_list, real_data_item) ! Dispatches to dlist_create

real_data_item%row_ind = 15
real_data_item%col_ind = 25
real_data_item%val     = 6.28_dp
CALL list_insert(real_elements_list, real_data_item) ! Inserts after the first element

count_real = list_count(real_elements_list)
PRINT *, "Number of real elements in list:", count_real

! Working with a list of complex sparse matrix elements
complex_data_item%row_ind = 5
complex_data_item%col_ind = 7
complex_data_item%val     = (1.0_dp, -2.5_dp)
CALL list_create(complex_elements_list, complex_data_item) ! Dispatches to zlist_create

count_complex = list_count(complex_elements_list)
PRINT *, "Number of complex elements in list:", count_complex

! Clean up
CALL list_destroy(real_elements_list)    ! Dispatches to dlist_destroy
CALL list_destroy(complex_elements_list) ! Dispatches to zlist_destroy
```

# Dependencies and Interactions

*   **`psp_dList`**: This module is a primary dependency. It provides the concrete definitions for `dList` (as `LINKED_LIST` specialized for `dNodeData`) and `dNodeData` (as `LIST_DATA`), along with the implementations of all list operations for real data (e.g., `dlist_create`).
*   **`psp_zList`**: Similarly, this module provides the concrete definitions for `zList` and `zNodeData`, and the implementations of list operations for complex data (e.g., `zlist_create`).
*   The `psp_dList` and `psp_zList` modules are themselves expected to `INCLUDE 'pspLinkedlist.F90'`, which contains the generic, data-type-agnostic linked list logic. `psp_dList` defines `LIST_DATA` as `dNodeData` before inclusion, and `psp_zList` defines `LIST_DATA` as `zNodeData`.
*   `pspListType` serves as an abstraction layer that is likely used by modules in `pspUtility` (like `pspListTool`) and various sparse matrix computation routines in `pspLevel3` that need to dynamically build lists of sparse matrix elements.
*   The `dListPtrArray` and `zListPtrArray` types are used in more complex scenarios, such as the sparse matrix-sparse matrix multiplication routines (`pspSpmSpm_tn`, `pspSpmSpm_nt`), where results for different parts of a matrix might be accumulated in separate lists before being combined.
