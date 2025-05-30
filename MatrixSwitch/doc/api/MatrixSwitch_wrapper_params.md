# Overview

The module `MatrixSwitch_wrapper_params` serves as a dedicated parameter and storage manager for the `MatrixSwitch_wrapper` module. Its main purpose is to handle the association between human-readable string keys (typically used by C code calling the wrapper) and the internal Fortran `TYPE(matrix)` objects that MatrixSwitch operates on. It defines the arrays to store these keys and matrices and provides the crucial lookup function to find a matrix by its key.

# Key Components

*   **Module `MatrixSwitch_wrapper_params`**:
    *   **Imports**:
        *   `MatrixSwitch, only: matrix`: Imports the core `TYPE(matrix)` definition from the main MatrixSwitch module.
        *   `MatrixSwitch_ops, only: dp, die`: Imports the `dp` kind parameter for double precision and the `die` subroutine for standardized error termination.
    *   **Parameters**:
        *   `max_key_length`: An integer parameter (set to 10) that defines the maximum number of characters allowed in a string key. This imposes a limit on the length of names that can be used to identify matrices in the wrapper system.
    *   **Module Variables (Public, Allocatable)**:
        *   `ms_keys(:)`: A 1D allocatable array of `CHARACTER(LEN=max_key_length)`. This array stores the unique string keys that identify the matrices managed by the wrapper.
        *   `ms_num_matrices`: An integer variable that holds the total number of matrices that the wrapper has been initialized to manage. This determines the allocated size of `ms_keys` and `ms_matrices`.
        *   `ms_matrices(:)`: A 1D allocatable array of `TYPE(matrix)`. This array is the actual container for the MatrixSwitch matrix objects being managed by the wrapper.
    *   **Public Function**:
        *   `ms_lookup(key)`: An integer function that takes a `CHARACTER(*)` `key` as input. It searches for this `key` within the `ms_keys` array. If the key is found, `ms_lookup` returns its 1-based index. This index directly corresponds to the position of the associated `TYPE(matrix)` object in the `ms_matrices` array. If the `key` is not found in `ms_keys` after checking all entries, the function calls `die` with the error message "ms_lookup: key not found", terminating the program.

# Important Variables/Constants

*   `max_key_length`: This constant dictates the maximum length of string keys. C code interacting with the wrapper must ensure its keys conform to this length.
*   `ms_keys`: Stores the string names for each matrix.
*   `ms_num_matrices`: Defines the capacity of the wrapper.
*   `ms_matrices`: Holds the actual Fortran matrix objects. The integrity of the wrapper system relies on `ms_keys(i)` always corresponding to `ms_matrices(i)`.
*   `ms_lookup` function: This is the central mechanism enabling the string-keyed access provided by `MatrixSwitch_wrapper`. Its correctness and efficiency are important for the wrapper's performance. The error handling (calling `die`) ensures that attempts to use invalid keys are caught.

# Usage Examples

This module is primarily for internal use by `MatrixSwitch_wrapper`. Its components are manipulated during the lifecycle of the wrapper:
1.  **Initialization (`ms_wrapper_open` in `MatrixSwitch_wrapper`)**:
    *   `ms_num_matrices` is assigned the number of matrices to be managed.
    *   `ms_keys` and `ms_matrices` are allocated to this size.
    *   The string keys provided by the user (e.g., from C) are stored in `ms_keys`.
2.  **Matrix Operations (e.g., `m_allocate` in `MatrixSwitch_wrapper`)**:
    *   A string key (e.g., "myMatrix1") is passed to a wrapper function.
    *   The wrapper function calls `idx = ms_lookup("myMatrix1")`.
    *   The actual operation is then performed on `ms_matrices(idx)`.
3.  **Finalization (`ms_wrapper_close` in `MatrixSwitch_wrapper`)**:
    *   `ms_matrices` elements are deallocated.
    *   `ms_keys` and `ms_matrices` arrays are deallocated.

# Dependencies and Interactions

*   **`MatrixSwitch`**: Dependency for the `TYPE(matrix)` definition.
*   **`MatrixSwitch_ops`**: Dependency for the `dp` precision kind and the `die` error handling routine.
*   **`MatrixSwitch_wrapper`**: This is the sole intended client of `MatrixSwitch_wrapper_params`. `MatrixSwitch_wrapper` uses the variables and the `ms_lookup` function provided by this module to implement its key-based API. The public accessibility of `ms_keys`, `ms_num_matrices`, and `ms_matrices` allows `MatrixSwitch_wrapper` to initialize and manage these central storage arrays.
