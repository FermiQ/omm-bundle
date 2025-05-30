# Overview

The `psp_zList.F90` file defines the `psp_zList` module. In parallel with `psp_dList.F90` (which handles real-valued lists), the `psp_zList` module specializes the generic linked list implementation from `pspLinkedlist.F90` for handling sparse matrix elements that have complex (double precision) values.

This specialization is achieved by defining the `LIST_DATA` type—a placeholder in the generic linked list code—to be equivalent to `zNode`. The `zNode` type is imported from the `pspNode` module and is structured to store a row index, a column index, and a complex value. By including the `pspLinkedlist.F90` source file after this type definition, all the generic list operations (like `list_create`, `list_insert`, etc.) become available and are specifically tailored for linked lists of these complex-valued sparse elements.

# Key Components

*   **Module `psp_zList`**:
    *   **Type Specialization via `USE` Statement**:
        *   `use pspNode, ONLY: LIST_DATA => zNode`: This line is central to the module's function. It imports the `zNode` derived type from the `pspNode` module. The `zNode` type (defined in `pspNode.F90`) includes fields for `row_ind` (integer), `col_ind` (integer), and `val` (complex double precision). By renaming/aliasing `zNode` to `LIST_DATA` within the scope of the `psp_zList` module, it provides the concrete data type definition that the generic linked list code included next will operate upon.

    *   **Inclusion of Generic Linked List Implementation**:
        *   `include "pspLinkedlist.F90"`: This Fortran directive textually includes the source code from the `pspLinkedlist.F90` file. This included file contains the definition of `TYPE LINKED_LIST` and a suite of generic subroutines and functions for list manipulation (e.g., `list_create`, `list_destroy`, `list_insert`, `list_count`). Because `LIST_DATA` has been defined as `zNode` immediately prior to this `INCLUDE` statement, all these generic routines are effectively instantiated to work specifically with linked lists where each node's data payload is of type `zNode`.

    *   **Exposed Entities**: Although the `psp_zList.F90` file itself does not list explicit `PUBLIC` attributes for the list operations (as these are defined in the included `pspLinkedlist.F90`), the effective public interface of the `psp_zList` module consists of:
        *   `TYPE(LINKED_LIST)`: Now specialized as a list of `zNode` elements.
        *   All the list manipulation procedures from `pspLinkedlist.F90` (e.g., `list_create`, `list_insert`), now operating on lists of `zNode` data.
        These entities are typically imported, renamed (e.g., `LINKED_LIST` to `zList`, `list_create` to `zlist_create`), and re-exported by the `pspListType` module to provide a unified and generic API to the rest of the pspBLAS library.

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined at compile time, it includes a `config.h` file, commonly used for managing build-specific configurations and preprocessor definitions.
*   The most critical definition within this module is `LIST_DATA => zNode`. This directive is what tailors the generic linked list code from `pspLinkedlist.F90` to handle complex sparse matrix elements. The `dp` kind parameter for the complex values is implicitly part of the `zNode` type definition (imported from `pspNode`).

# Usage Examples

The `psp_zList` module is primarily an internal building block for the pspBLAS library. Higher-level modules, most notably `pspListType`, will `USE psp_zList` to make its specialized types and list operation routines available. End-users of the pspBLAS library would typically interact with these complex-valued linked lists through the generic interfaces provided by `pspListType`.

If `psp_zList` were to be used directly (which is less common than using `pspListType`):

```fortran
USE psp_zList  ! This makes LINKED_LIST (as a list of zNode)
               ! and list_create, etc. (specialized for zNode) available.
               ! LIST_DATA here is effectively zNode.
USE pspNode    ! To explicitly use zNode if needed, or ensure dp is in scope
IMPLICIT NONE

TYPE(LINKED_LIST), POINTER :: head_complex_list => NULL()
TYPE(LIST_DATA)    :: complex_element_data
! LIST_DATA is an alias for zNode in this context.
! Alternatively, explicitly: TYPE(zNode) :: complex_element_data

complex_element_data%row_ind = 7
complex_element_data%col_ind = 12
complex_element_data%val     = (2.5_dp, -3.1_dp) ! Assuming dp is accessible

CALL list_create(head_complex_list, complex_element_data)

! ... perform other list operations like list_insert, list_count, etc. ...

CALL list_destroy(head_complex_list)
```
More typically, the `pspListType` module is used, which renames `LINKED_LIST` (from `psp_zList`) to `zList` and `list_create` (from `psp_zList`) to `zlist_create`, etc., providing a clearer, generic interface.

# Dependencies and Interactions

*   **`pspNode`**: This module is a direct and crucial dependency. `psp_zList` uses `pspNode` to import the `zNode` type and assign it as the alias for `LIST_DATA`. The `zNode` type defines the structure (row index, column index, complex value) of the data stored in each node of the linked lists created using this module.
*   **`pspLinkedlist.F90` (via `INCLUDE` statement)**: This uncompiled source file is textually included and provides the generic, data-type-agnostic Fortran code for the linked list structure (`TYPE LINKED_LIST`) and all associated list manipulation subroutines and functions. The behavior of these included routines is specialized by the `LIST_DATA => zNode` definition established within `psp_zList`.
*   **`pspListType`**: This higher-level module is a primary consumer of `psp_zList`. `pspListType` `USE`s `psp_zList` to obtain the complex-specific list type (which it typically renames to `zList`) and complex-specific list procedures (which it renames to `zlist_create`, `zlist_insert`, etc.). It then provides generic interfaces that map to these (and their real counterparts from `psp_dList`).
*   Other pspBLAS modules that need to work with linked lists of complex-valued sparse matrix elements will typically do so through the abstractions provided by `pspListType`, which in turn relies on `psp_zList`.
