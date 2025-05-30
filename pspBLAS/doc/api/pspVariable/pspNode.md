# Overview

The `pspNode.F90` file defines the `pspNode` module, which is a part of the `pspVariable` collection within the pspBLAS library. The primary and specific purpose of this module is to define the data structures (derived types) that represent individual non-zero elements of sparse matrices when they are stored in linked lists. These structures serve as the `LIST_DATA` payload within the generic linked list implementation provided by `pspLinkedlist.F90`.

The module provides two distinct types for this purpose:
*   `dNode`: For sparse matrix elements with real (double precision) values.
*   `zNode`: For sparse matrix elements with complex (double precision) values.

These types simply group the row index, column index, and the numerical value of a sparse matrix element.

# Key Components

*   **Module `pspNode`**:
    *   **Public Derived Type `dNode`**:
        *   **Purpose**: To define the structure for storing a single real-valued non-zero element of a sparse matrix.
        *   **Components**:
            *   `row_ind`: `INTEGER`, stores the 1-based row index of the non-zero element.
            *   `col_ind`: `INTEGER`, stores the 1-based column index of the non-zero element.
            *   `val`: `REAL(dp)`, stores the double precision real value of the non-zero element. The `dp` kind parameter is defined within this module.

    *   **Public Derived Type `zNode`**:
        *   **Purpose**: To define the structure for storing a single complex-valued non-zero element of a sparse matrix.
        *   **Components**:
            *   `row_ind`: `INTEGER`, stores the 1-based row index of the non-zero element.
            *   `col_ind`: `INTEGER`, stores the 1-based column index of the non-zero element.
            *   `val`: `COMPLEX(dp)`, stores the double precision complex value of the non-zero element.

# Important Variables/Constants

*   `dp`: An integer parameter defined as `selected_real_kind(15,300)`, specifying the kind for double precision real numbers. This ensures consistent precision for the `val` component in `dNode` and `zNode`.
*   `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard complex constants (`(1.0_dp,0.0_dp)`, `(0.0_dp,1.0_dp)`, `(0.0_dp,0.0_dp)` respectively) of kind `dp`. These are defined but not directly used within the type definitions themselves.
*   `HAVE_CONFIG_H`: A preprocessor macro. If defined, it includes a `config.h` file, typically used for managing build configurations.

# Usage Examples

The `dNode` and `zNode` types defined in this module are primarily intended to be used as the specific structure for `LIST_DATA` when the generic linked list code (`pspLinkedlist.F90`) is instantiated. Modules like `psp_dList.F90` and `psp_zList.F90` make use of these types (or define identical structures often named `dNodeData` and `zNodeData`).

```fortran
! Conceptual use within a routine that populates a list of real sparse elements:
! (Assuming list operations are available from an instantiated linked list module)
USE pspNode ! To get the definition of dNode (and zNode)
IMPLICIT NONE

TYPE(dNode) :: new_sparse_element
TYPE(dList), POINTER :: my_list_head => NULL() ! Assuming dList is the specialized LINKED_LIST type

! Populate a dNode variable
new_sparse_element%row_ind = 10
new_sparse_element%col_ind = 25
new_sparse_element%val     = 1.2345_dp

! Add it to a linked list (conceptual call, actual call depends on specialized list module)
! CALL list_add_element(my_list_head, new_sparse_element)
```

In practice, `psp_dList` and `psp_zList` define their own `dNodeData` and `zNodeData` types that mirror `dNode` and `zNode` structure, and then include `pspLinkedlist.F90`.

# Dependencies and Interactions

*   This module is self-contained in terms of its type definitions but defines the `dp` parameter locally.
*   **`psp_dList` and `psp_zList` Modules**: These modules are the primary consumers of the type definitions from `pspNode`. They typically rename `dNode` to `dNodeData` (and `zNode` to `zNodeData`) and use this as the `LIST_DATA` type when they `INCLUDE 'pspLinkedlist.F90'`. This mechanism allows the generic linked list operations in `pspLinkedlist.F90` to work with nodes containing sparse matrix element data (row, column, value).
*   **`pspListTool`**: Higher-level routines in `pspListTool` that operate on `dList` or `zList` types will indirectly be working with data structured according to `dNode` or `zNode`.

The `pspNode` module plays a simple but vital role by providing standardized structures for the data content of linked list nodes used for sparse matrix computations in pspBLAS.
