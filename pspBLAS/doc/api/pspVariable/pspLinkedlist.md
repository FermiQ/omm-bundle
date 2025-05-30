# Overview

The `pspLinkedlist.F90` file provides a generic template for implementing singly linked lists in Fortran. It is not a self-contained, directly usable module in the typical Fortran sense (i.e., via a `USE module_name` statement for `pspLinkedlist` itself). Instead, it's designed as a piece of source code to be `INCLUDE`d within other Fortran modules. The including module must define a specific data type named `LIST_DATA`, which specifies the actual data to be stored in each node of the linked list. Once `LIST_DATA` is defined, including `pspLinkedlist.F90` effectively instantiates a suite of linked list operations for that specific data type.

The comments in the file indicate that this code is a modified version of `linkedlist.f90` by Arjen Markus and emphasize that variables holding these lists should typically be pointers.

# Key Components

*   **`TYPE LINKED_LIST`**:
    *   This is the core derived type defining the structure of a node in the linked list. Each node comprises:
        *   `next`: A `POINTER` to another object of `TYPE(LINKED_LIST)`. This pointer links the current node to the subsequent node in the list, or to `NULL()` if it's the last node.
        *   `data`: A variable of `TYPE(LIST_DATA)`. This component holds the actual data payload for the node. The concrete definition of `LIST_DATA` must be provided by the context (e.g., the module) that includes `pspLinkedlist.F90`.

*   **Subroutines and Functions for List Operations**:
    The file defines a standard set of procedures to manage these linked lists:
    *   `list_create(list, data)`: Allocates a new list (a single node), initializes its `data` field with the provided `data` argument, and sets its `next` pointer to `NULL()`. The `list` argument becomes a pointer to this new node.
    *   `list_destroy(list)`: Deallocates all nodes in an entire linked list, starting from the node pointed to by `list`. It traverses the list and deallocates each node encountered. This routine assumes that the `LIST_DATA` components do not contain pointers or allocatable arrays that require separate, deep deallocation.
    *   `list_count(list)`: An integer function that traverses the list and returns the total number of elements (nodes) it contains. Returns 0 if the list is not associated (i.e., `NULL()`).
    *   `list_next(elem)`: A function that takes a pointer to a list element (`elem`) and returns a pointer to the next element in the list (`elem%next`).
    *   `list_insert(elem, data)`: Inserts a new element, initialized with `data`, immediately after the specified existing element `elem` in the list.
    *   `list_insert_head(list, data)`: Inserts a new element, initialized with `data`, at the beginning of the `list`. The `list` pointer is updated to point to this new head element.
    *   `list_delete_element(list, elem)`: Removes the specified element `elem` from the `list`. It correctly handles cases where `elem` is the head of the list or an intermediate node.
    *   `list_get_data(elem)`: A function that takes a pointer to a list element `elem` and returns a copy of the `LIST_DATA` stored within that element.
    *   `list_put_data(elem, data)`: Updates the `data` field of an existing list element `elem` with the new `data` provided.

*   **Commented-out `interface assignment(=)`**:
    *   The file contains a commented-out attempt to define a private interface for the assignment operator (`=`) for `LINKED_LIST` types. The purpose was likely to prevent direct pointer assignment (`list_left = list_right`) which would lead to aliasing rather than creating a new copy of the list, or to flag accidental use of `=` where `=>` (pointer assignment) was intended for list manipulation. The comment indicates this approach had a "private/public conflict," suggesting issues with visibility or type resolution in its original context.

# Important Variables/Constants

*   This file itself does not define module-level constants like `dp` (double precision kind) or specific complex number constants. It relies on these being available from the scope where `pspLinkedlist.F90` is included, or through `USE` association if `LIST_DATA` requires them (e.g., if `LIST_DATA` contains `REAL(dp)` fields).
*   The most critical external "variable" is the **`TYPE(LIST_DATA)` definition**. The entire set of routines is generic with respect to this type, which must be defined before `pspLinkedlist.F90` is included.

# Usage Examples

`pspLinkedlist.F90` provides a generic template. Its actual use is instantiated within other modules that define a specific `LIST_DATA` type. For example, `psp_dList.F90` (for lists of real sparse matrix elements) and `psp_zList.F90` (for complex sparse matrix elements) in the pspBLAS library define `dNodeData` and `zNodeData` as their respective `LIST_DATA` types and then include `pspLinkedlist.F90` to get these list operations tailored for their specific data.

A conceptual illustration (assuming `LIST_DATA` is defined as a type holding a single integer):

```fortran
! MODULE my_custom_integer_list
!   ! Assume dp or other kinds are available if LIST_DATA needs them
!   IMPLICIT NONE
!   PRIVATE
!
!   TYPE :: LIST_DATA  ! Definition of the data stored in each list node
!     INTEGER :: item_value
!   END TYPE LIST_DATA
!
!   INCLUDE 'pspLinkedlist.F90' ! This makes list_create, list_insert, etc.,
!                               ! available, specialized for the above LIST_DATA.
!
!   PUBLIC :: list_create       ! Expose the desired list operations
!   PUBLIC :: list_insert
!   PUBLIC :: list_count
!   PUBLIC :: list_destroy
!   ! ... other public procedures derived from the INCLUDEd file ...
!
! CONTAINS
!   ! No additional subroutines needed here as they come from the INCLUDE
! END MODULE my_custom_integer_list
!
! PROGRAM test_integer_list
!   USE my_custom_integer_list
!   IMPLICIT NONE
!
!   TYPE(LINKED_LIST), POINTER :: head_ptr => NULL()
!   TYPE(LIST_DATA) :: current_data
!   INTEGER :: count
!
!   current_data%item_value = 100
!   CALL list_create(head_ptr, current_data)
!
!   current_data%item_value = 200
!   CALL list_insert(head_ptr, current_data) ! Inserts after the first element
!
!   count = list_count(head_ptr)
!   PRINT *, "Number of items in list:", count ! Output: 2
!
!   CALL list_destroy(head_ptr)
! END PROGRAM test_integer_list
```

# Dependencies and Interactions

*   **External `TYPE(LIST_DATA)` Definition**: This is the primary dependency. The code in `pspLinkedlist.F90` is a template that becomes concrete only when `LIST_DATA` is defined by the including module or scope.
*   **Inclusion Mechanism**: This file is designed to be incorporated into other modules using the Fortran `INCLUDE` statement. This allows the generic list logic to be "stamped out" for different data types without code duplication. Examples in pspBLAS are `psp_dList.F90` and `psp_zList.F90`.
*   **Memory Management**: The routines handle allocation and deallocation of `LINKED_LIST` nodes. They assume shallow copies for `LIST_DATA`; if `LIST_DATA` contains pointers or allocatable components requiring deep copies or specialized memory management, the user of this template must handle that externally or by modifying these routines.
*   The `use pspVariable` statement at the top of the file in the provided context suggests that even this generic template might be intended to be included in an environment where `pspVariable` (providing things like `dp`) is already available, or it's placed there for structural consistency within the pspBLAS library. The core linked-list logic itself doesn't directly use entities from `pspVariable` unless `LIST_DATA` does.
