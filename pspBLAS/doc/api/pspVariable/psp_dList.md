# Overview

The `psp_dList.F90` file defines the `psp_dList` module. The primary purpose of this module is to instantiate and specialize a generic linked list implementation for handling sparse matrix elements that have real (double precision) values. It does this by leveraging a common linked list template file (`pspLinkedlist.F90`) and providing a specific data type for the elements to be stored in the list.

The module effectively creates a concrete linked list implementation for "double precision nodes" (`dNode`), where each node typically stores a row index, a column index, and a real-valued matrix element.

# Key Components

*   **Module `psp_dList`**:
    *   **Type Specialization via `USE` Statement**:
        *   `use pspNode, ONLY: LIST_DATA => dNode`: This crucial line imports the `dNode` derived type from the `pspNode` module. The `dNode` type (defined in `pspNode.F90`) contains fields for `row_ind` (integer), `col_ind` (integer), and `val` (real double precision). By renaming `dNode` to `LIST_DATA` within the scope of the `psp_dList` module, it provides the necessary concrete data type definition that the generic linked list code expects.

    *   **Inclusion of Generic Linked List Implementation**:
        *   `include "pspLinkedlist.F90"`: This directive textually includes the source code from `pspLinkedlist.F90`. The code in `pspLinkedlist.F90` defines `TYPE LINKED_LIST` and a suite of list manipulation subroutines and functions (e.g., `list_create`, `list_insert`, `list_destroy`). Because `LIST_DATA` has been defined as `dNode` immediately before this inclusion, all these generic routines become specialized operations for linked lists of `dNode` elements.

    *   **Exposed Entities**: Although this module file itself does not contain explicit `PUBLIC` statements for the list operations, the effect of including `pspLinkedlist.F90` is that `TYPE(LINKED_LIST)` (now effectively a list of `dNode`s) and all the list manipulation procedures like `list_create`, `list_insert`, etc. (now specialized for `dNode` data) become part of the `psp_dList` module. These are then typically renamed and re-exported by the `pspListType` module to provide a more generic interface to the rest of the library (e.g., `dList` for the type, `dlist_create` for the creation subroutine).

# Important Variables/Constants

*   `HAVE_CONFIG_H`: A preprocessor macro. If defined, it includes a `config.h` file, used for build configurations.
*   The module implicitly uses the `dp` parameter (for double precision kind) defined within `pspNode` (as `dNode` uses it) and also within `pspLinkedlist.F90` if it were to define its own, though `pspLinkedlist.F90` is a template expecting `LIST_DATA` to be fully defined.

# Usage Examples

The `psp_dList` module is primarily an internal building block for the pspBLAS library. Higher-level modules, particularly `pspListType`, will `USE psp_dList` to make its specialized types and routines available. Users of pspBLAS would typically interact with these lists via the generic interfaces provided by `pspListType`.

If one were to use `psp_dList` directly (which is less common than using `pspListType`):

```fortran
USE psp_dList  ! This makes LINKED_LIST (as a list of dNode)
               ! and list_create, etc. (specialized for dNode) available.
               ! LIST_DATA here is effectively dNode.
IMPLICIT NONE

TYPE(LINKED_LIST), POINTER :: head_of_list => NULL()
TYPE(LIST_DATA)    :: list_item_data ! LIST_DATA is dNode from pspNode via this module

list_item_data%row_ind = 10
list_item_data%col_ind = 20
list_item_data%val     = 12.345_dp ! Assuming dp is available

CALL list_create(head_of_list, list_item_data)

! ... perform other list operations using list_insert, list_count, etc. ...

CALL list_destroy(head_of_list)
```
More typically, one would use `pspListType` which renames `LINKED_LIST` to `dList` and `list_create` to `dlist_create` etc., for clarity.

# Dependencies and Interactions

*   **`pspNode`**: This module is a direct dependency. `psp_dList` uses `pspNode` to import the `dNode` type and assign it as the alias for `LIST_DATA`. The `dNode` type defines the structure (row index, column index, real value) of the data stored in each node of the linked lists created using this module.
*   **`pspLinkedlist.F90` (via `INCLUDE` statement)**: This uncompiled source file is textually included. It provides the generic, data-type-agnostic Fortran code for the linked list structure (`TYPE LINKED_LIST`) and all associated list manipulation subroutines and functions. The behavior of these included routines is specialized by the `LIST_DATA => dNode` definition within `psp_dList`.
*   **`pspListType`**: This module is a primary consumer of `psp_dList`. `pspListType` uses `psp_dList` to obtain the real-specific list type (which it calls `dList`) and real-specific list procedures (which it calls `dlist_create`, `dlist_insert`, etc.). It then provides generic interfaces that map to these (and their complex counterparts from `psp_zList`).
*   Other pspBLAS modules that need to work with linked lists of real-valued sparse matrix elements will typically do so through the abstractions provided by `pspListType`, which in turn relies on `psp_dList`.
