# Overview

The `MatrixSwitch_wrapper` module provides a Fortran-based API that is designed to be easily callable from C, effectively acting as C bindings for the MatrixSwitch library. It achieves this by managing an internal array of MatrixSwitch `TYPE(matrix)` objects. These objects are not manipulated directly by the C caller but are referenced through character string keys. The wrapper then translates operations specified by these keys and parameters into calls to the core MatrixSwitch Fortran routines.

# Key Components

*   **Module `MatrixSwitch_wrapper`**:
    *   **Imports**:
        *   `MatrixSwitch_ops`: For the `dp` (double precision) kind parameter.
        *   `MatrixSwitch_wrapper_params`: This crucial module provides the storage for matrices (`ms_matrices`), their corresponding keys (`ms_keys`), the total number of manageable matrices (`ms_num_matrices`), and the `ms_lookup` function which maps a string key to an integer index for the `ms_matrices` array.
        *   `MatrixSwitch`: The core library. Most public procedures from `MatrixSwitch` are imported with an `_orig` suffix (e.g., `mm_multiply_orig`) to distinguish them from the wrapper's own interface procedures. Conditional compilation includes MPI, ScaLAPACK, and pspBLAS related functionalities.
    *   **Public Wrapper Interfaces and Subroutines**:
        *   `ms_wrapper_open(num_matrices, keys)`: Initializes the wrapper. It allocates internal storage for `num_matrices` `TYPE(matrix)` objects and stores their associated string `keys`.
        *   `ms_wrapper_close()`: Finalizes the wrapper. It deallocates all `TYPE(matrix)` objects stored within `ms_matrices` (by calling `m_deallocate_orig` on each) and then deallocates the `ms_matrices` and `ms_keys` arrays themselves.
        *   **Property Query Functions**: `ms_is_initialized(m_name)`, `ms_is_serial(m_name)`, `ms_is_real(m_name)`, `ms_is_square(m_name)`, `ms_is_sparse(m_name)`, `ms_dim1(m_name)`, `ms_dim2(m_name)`. These functions take a matrix string key (`m_name`), look up the corresponding matrix object, and return the requested property.
        *   **Core Operation Wrappers**:
            *   `m_allocate(m_name, i, j, label)`
            *   `m_deallocate(m_name)`
            *   `m_copy(m_name, A_key, label, threshold, threshold_is_soft)`
            *   `m_convert(m_name, label, threshold, threshold_is_soft)`
            *   These subroutines take string key(s) for the matrix/matrices involved, look up the actual `TYPE(matrix)` object(s) using `ms_lookup`, and then call the corresponding `_orig` procedure from the `MatrixSwitch` module.
        *   **Computational Operation Wrappers**: Interfaces like `mm_multiply`, `m_add`, `m_trace`, `mm_trace`, `m_scale`, `m_set`, `m_set_element`, `m_get_element` are provided. Each interface typically has real (`_d`) and complex (`_z`) versions (e.g., `mm_dmultiply`, `mm_zmultiply`). These also use `ms_lookup` to get the target matrix objects before calling the `_orig` routines.
        *   **Registration and External Copy Wrappers**: Similar wrappers are provided for `m_register_sden`, `m_register_pdbc` (MPI), `m_copy_external_pdbcpcoo` (pspBLAS), `m_copy_external_pdbcpcsc` (pspBLAS), `m_register_pcoo` (pspBLAS), `m_register_pcsc` (pspBLAS), forwarding calls to their `_orig` counterparts.
        *   **Parallel Setup Re-exports**: `ms_scalapack_setup` and `ms_lap_icontxt` are re-exported if MPI and relevant libraries are enabled.

# Important Variables/Constants

*   `dp`: Double precision kind parameter (from `MatrixSwitch_ops`).
*   Implicitly from `MatrixSwitch_wrapper_params`:
    *   `ms_matrices()`: An array of `TYPE(matrix)` from `MatrixSwitch_ops`, holding the actual matrix data structures managed by the wrapper.
    *   `ms_keys()`: A character array storing the string keys corresponding to matrices in `ms_matrices`.
    *   `ms_num_matrices`: An integer specifying the capacity of the wrapper.
    *   `ms_lookup(character_key)`: A function (critical for the wrapper) that takes a string key and returns the integer index of the corresponding matrix in `ms_matrices`.
*   Preprocessor Flags: `HAVE_CONFIG_H`, `HAVE_MPI`, `HAVE_SCALAPACK`, `HAVE_PSPBLAS` determine which wrapper functionalities are compiled, especially for parallel operations.

# Usage Examples

This module is intended to be the Fortran side of a C-Fortran interoperability layer. A C program would call functions (defined with `BIND(C)`) which in turn call the procedures in this module.
1.  **Initialization (from C, via a BIND(C) Fortran function)**:
    `call ms_wrapper_open(5, ["matrixA", "matrixB", "matrixC", "temp1", "temp2"])`
2.  **Operations (from C, via BIND(C) functions)**:
    `call m_allocate("matrixA", 100, 100, "sdden")`
    `call mm_dmultiply("matrixA", 'N', "matrixB", 'N', "matrixC", 1.0_dp, 0.0_dp)`
3.  **Cleanup (from C, via a BIND(C) Fortran function)**:
    `call ms_wrapper_close()`

The Fortran procedures in this module use `ms_lookup(m_name)` to get the actual `TYPE(matrix)` object from the `ms_matrices` array based on the provided string key.

# Dependencies and Interactions

*   **`MatrixSwitch_ops`**: Required for the `dp` kind parameter. The `TYPE(matrix)` itself is used via `ms_matrices` from `MatrixSwitch_wrapper_params`.
*   **`MatrixSwitch_wrapper_params`**: This is a core dependency, providing the matrix storage array (`ms_matrices`), key array (`ms_keys`), and the essential `ms_lookup` function.
*   **`MatrixSwitch`**: The primary functional dependency. All matrix operations are ultimately delegated to the `_orig` procedures imported from this module.
*   **C Standard Library / Fortran `ISO_C_BINDING`**: Although not directly visible in this file's code, for C interoperability, the public procedures would typically be exposed to C using `BIND(C, name="c_function_name")` and by using C-compatible types from `ISO_C_BINDING`. This file defines the Fortran target for such C calls.
*   **MPI, ScaLAPACK, pspBLAS**: Conditional dependencies based on build configurations, affecting which parallel operations are wrapped and exposed.
