# Overview

The `omm_wrapper.F90` file provides a single public subroutine, also named `omm_wrapper`. This subroutine serves as an interface layer, primarily intended to be called from a C environment (or other languages via C bindings). Its main function is to bridge the gap between an API where matrices are identified by string keys and the core Fortran `omm` routine (from `omm.F90`), which operates on native Fortran `TYPE(matrix)` objects provided by the `MatrixSwitch` library.

The `omm_wrapper` subroutine accepts matrix arguments as character strings. It then uses the lookup mechanism provided by `MatrixSwitch_wrapper_params` (specifically, `ms_lookup` and `ms_matrices`) to retrieve the actual `TYPE(matrix)` objects associated with these keys. Once these matrix objects are resolved, `omm_wrapper` calls the main `omm` calculation routine, passing along the resolved matrices and all other computational parameters. It also handles the potential type conversion for logical arguments when C interoperability (via `iso_c_binding`) is active.

# Key Components

*   **Subroutine `omm_wrapper(...)`**:
    *   **Purpose**: To act as an intermediary, allowing the core `omm` functionality to be invoked using string keys for matrix identification.
    *   **Input Arguments**:
        *   `m`: Integer, size of the basis set.
        *   `n`: Integer, number of occupied states (wavefunctions).
        *   `H`: `CHARACTER(*)`, string key for the Hamiltonian matrix.
        *   `S`: `CHARACTER(*)`, string key for the Overlap matrix (or its Cholesky factor).
        *   `D_min`: `CHARACTER(*)` (inout via underlying `omm`), string key for the Density matrix (or energy-weighted density matrix).
        *   `C_min`: `CHARACTER(*)` (inout via underlying `omm`), string key for the Wavefunction coefficients matrix.
        *   `T`: `CHARACTER(*)`, string key for the Kinetic energy matrix (used for preconditioning).
        *   `new_S`: `LOGICAL` (or `LOGICAL(c_bool)` if `CBIND` is defined), flag indicating if `S` is new for the current spin/k-point.
        *   `calc_ED`: `LOGICAL` (or `LOGICAL(c_bool)`), flag to calculate energy-weighted density matrix.
        *   `init_C`: `LOGICAL` (or `LOGICAL(c_bool)`), flag indicating if `C_min` is externally initialized.
        *   `long_out`: `LOGICAL` (or `LOGICAL(c_bool)`), flag for enabling detailed logging.
        *   `dealloc`: `LOGICAL` (or `LOGICAL(c_bool)`), flag to deallocate internal persistent matrices.
        *   `flavour`: Integer, specifies the OMM variant.
        *   `np`: Integer, total number of spin/k-points.
        *   `ip`: Integer, current spin/k-point index.
        *   `eta`: `REAL(dp)`, OMM energy shift parameter.
        *   `cg_tol`: `REAL(dp)`, convergence tolerance for CG minimization.
        *   `scale_T`: `REAL(dp)`, scale factor for kinetic energy in preconditioner.
        *   `m_storage`: `CHARACTER(5)`, MatrixSwitch label for storage format.
        *   `m_operation`: `CHARACTER(3)`, MatrixSwitch label for operation implementation.
    *   **Output Arguments**:
        *   `e_min`: `REAL(dp)`, the calculated OMM functional energy.
    *   **Method**:
        1.  **Logical Argument Conversion**: The input logical flags (which might be `LOGICAL(c_bool)` if compiled with `CBIND` for C compatibility) are copied to local Fortran `LOGICAL` variables (e.g., `new_S_conv = new_S`). This ensures type safety when these flags are passed to the main `omm` Fortran subroutine.
        2.  **Matrix Resolution and Call to `omm`**: For each character string key representing a matrix (`H`, `S`, `D_min`, `C_min`, `T`), the wrapper uses `ms_lookup(key)` (from `MatrixSwitch_wrapper_params`) to get an integer index. This index is then used to access the corresponding `TYPE(matrix)` object from the `ms_matrices` array (also from `MatrixSwitch_wrapper_params`). These resolved `TYPE(matrix)` objects, along with all other scalar and converted logical parameters, are then passed to the main `omm` subroutine (from `omm.F90`) to perform the actual OMM calculation.

# Important Variables/Constants

*   **`CBIND`**: A preprocessor macro. If defined, `iso_c_binding` is used, and logical arguments are declared as `LOGICAL(c_bool)` to match C's `bool` type, facilitating direct calls from C.
*   **Matrix String Keys**: The input arguments `H, S, D_min, C_min, T` are character strings that serve as identifiers for matrices managed by the `MatrixSwitch_wrapper` system.
*   **Converted Logical Flags**: Local variables like `new_S_conv`, `calc_ED_conv`, `init_C_conv`, `long_out_conv`, `dealloc_conv` are used to hold standard Fortran `LOGICAL` versions of the input logical flags.

# Usage Examples

The `omm_wrapper` subroutine is typically not called directly from other Fortran code within the libOMM project. Instead, it is designed to be the Fortran procedure called by a C function that acts as the C Application Programming Interface (API) for libOMM's `omm` functionality.

A C program would first use other C-wrapper functions (interfacing with `MatrixSwitch_wrapper`) to initialize the matrix management system, allocate and register matrices with specific keys. Then, it would call a C function, say `omm_c_interface(...)`, which internally invokes this Fortran `omm_wrapper` subroutine using the C interoperability features.

```c
/* Conceptual C code calling a C function that wraps the Fortran omm_wrapper */
// Assume matrices "H_matrix_key", "S_matrix_key", etc., have been previously
// set up and populated using other C wrapper calls.

double energy_min_out;
// ... other C variables for parameters ...

int result = omm_c_interface(
    m_val, n_val,
    "H_matrix_key", "S_matrix_key",
    is_S_new_flag,
    &energy_min_out,
    "D_matrix_key",
    calc_ED_flag,
    eta_val,
    "C_matrix_key",
    is_C_initialized_flag,
    "T_matrix_key",
    scale_T_val,
    flavour_val,
    np_val, ip_val,
    cg_tol_val,
    long_out_flag,
    dealloc_flag,
    "sdden", /* m_storage */
    "ref"    /* m_operation */
);
// Now energy_min_out and the matrix associated with "D_matrix_key" are updated.
```

# Dependencies and Interactions

*   **`omm_params`**: Used for the `dp` (double precision) kind parameter.
*   **`MatrixSwitch_wrapper`**: This module is `USE`d to bring the main `omm` subroutine (presumably the one from `omm.F90`, possibly aliased or re-exported) into the current scope.
*   **`MatrixSwitch_wrapper_params`**: This is a critical dependency. `omm_wrapper` relies on `ms_matrices` (the array holding the actual `TYPE(matrix)` objects) and `ms_lookup` (the function that converts a string key into an index for `ms_matrices`) from this module.
*   **`omm` subroutine (from `omm.F90`)**: This is the core computational routine that `omm_wrapper` ultimately calls after resolving the matrix string keys into their corresponding `TYPE(matrix)` objects.
*   **`iso_c_binding`** (Module): This Fortran standard module is used if the `CBIND` preprocessor macro is defined. It provides `c_bool` and other tools for ensuring that Fortran procedures can be correctly called from C.
*   The `omm_wrapper` acts as a vital translation layer, enabling external code (typically C) that identifies matrices by string names to interface seamlessly with the underlying Fortran `omm` routine which requires direct `TYPE(matrix)` objects. It abstracts the matrix object management handled by the `MatrixSwitch_wrapper` system.
