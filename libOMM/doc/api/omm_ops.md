# Overview

The `omm_ops` module serves as a library of specialized matrix and algebraic operations tailored for the Orbital Minimization Method (OMM) implemented in libOMM. These operations are the computational workhorses used within the main `omm` and `omm_callback` subroutines. The module includes functions for calculating the OMM gradient, transforming matrices between different bases (Atomic Orbital to Wavefunction basis and vice-versa), performing Cholesky factorization, matrix inversion, applying Cholesky-related transformations, computing coefficients for the quartic line search used in OMM's CG minimizer, and constructing preconditioners. Many routines are generic interfaces that delegate to specific real or complex arithmetic versions, which in turn may call standard libraries like LAPACK and ScaLAPACK if available.

# Key Components

*   **Module `omm_ops`**:
    *   `calc_G(HW, SW, G, HC, SC, m_operation)`: Computes the gradient of the OMM functional: `G = 2*(2*H*C - S*C*HW - H*C*SW)`. `HC` and `SC` are precomputed `H*C` and `S*C`.
    *   `calc_AW(A, C, AW, AC, m_operation)`: Calculates an operator `A` in the Wavefunction (WF) basis: `AW = C^T*A*C`. `AC` is a workspace for `A*C`.
    *   `calc_HW_callback(H_callback, C, HW, HC, m_operation)`: Similar to `calc_AW`, but specifically for the Hamiltonian. It uses a callback function `H_callback` (passed as an argument) to compute `HC = H*C`.
        *   **Interface `H_callback(C, HC)`**: Defines the signature for the Hamiltonian application callback: `SUBROUTINE H(C, HC)`.
    *   `calc_A(AW, C1, A, CAW, m_operation)`: Transforms an operator `AW` from WF basis back to Atomic Orbital (AO) basis: `A = C1^T*AW*C1`. `CAW` is workspace.
    *   `calc_A2(AW, C1, C2, A, CAW, m_operation)`: Computes `A = C1^T*AW*C2`, a more general AO basis transformation.
    *   `calc_coeff(HW, SW, HWd, SWd, HWdd, SWdd, SWdH, coeff, m_operation)`: Calculates the five coefficients `coeff(0:4)` of the quartic polynomial `E(lambda) = coeff(4)*lambda^4 + ... + coeff(0)` for the OMM line search. Inputs are various trace terms involving `H` and `S` matrices in different bases/projections.
    *   `m_factorize(C, label)`: Interface for Cholesky factorization `C = U^T*U` (stores `U` in `C`).
        *   `m_dfactorize(C, label)`: Real version (LAPACK `dpotrf` / ScaLAPACK `pdpotrf`).
        *   `m_zfactorize(C, label)`: Complex version (LAPACK `zpotrf` / ScaLAPACK `pzpotrf`).
    *   `m_reduce(A, C, label)`: Interface to transform `C` to `A^-T*C*A^-1` where `A` is an upper triangular Cholesky factor. Used to reduce a generalized eigenvalue problem.
        *   `m_dreduce(A, C, label)`: Real version (LAPACK `dsygst` / ScaLAPACK `pdsygst`).
        *   `m_zreduce(A, C, label)`: Complex version (LAPACK `zhegst` / ScaLAPACK `pzhegst`).
    *   `m_transform(A, C, label)`: Interface for `C := C*A^-1` (where `A` is typically `U^T` from Cholesky `S=U^T*U`).
        *   `m_dtransform(A, C, label)`: Real version (LAPACK `dtrmm` / ScaLAPACK `pdtrmm`).
        *   `m_ztransform(A, C, label)`: Complex version (LAPACK `ztrmm` / ScaLAPACK `pztrmm`).
    *   `m_back_transform(A, C, label)`: Interface for `C := C*A` (where `A` is typically `U^T`).
        *   `m_dback_transform(A, C, label)`: Real version (LAPACK `dtrsm` / ScaLAPACK `pdtrsm`). Note: `dtrsm` solves a triangular system; this effectively multiplies by the inverse if `A` is the matrix in `dtrsm`. The comment `C:=C*A` might be an oversimplification if `dtrsm(A_inv, C)` is intended. However, `dtrsm` with `alpha=1.0` and `SIDE='R'` computes `X B = C -> X = C B^-1`. If `A` is `U^T`, then `C_new = C_old * (U^T)^-1`. This seems to be consistent with `m_transform` being `C*U^-1`. The typical sequence is `C_ortho = C_ao * U^-1` and `C_ao = C_ortho * U`. `dtrsm` for `C_new = C_old * U` would be `SIDE='R', TRANSA='N'`. The current `dtrsm` usage (`SIDE='R', TRANSA='T'`) computes `C_new = C_old * (U^T)^-1`. This needs careful check against typical usage.
        *   `m_zback_transform(A, C, label)`: Complex version (LAPACK `ztrsm` / ScaLAPACK `pztrsm`). Same consideration for `ztrsm`.
    *   `m_inverse(C, label)`: Interface for in-place matrix inversion `C := C^-1`.
        *   `m_dinverse(C, label)`: Real symmetric/Hermitian (LAPACK `dsytrf`/`dsytri` or ScaLAPACK `pdgetrf`/`pdgetri` for general square).
        *   `m_zinverse(C, label)`: Complex Hermitian (LAPACK `zhetrf`/`zhetri` or ScaLAPACK `pzgetrf`/`pzgetri` for general square).
    *   `calc_PW_precon(T, scale_T, P)`: Interface for plane-wave preconditioner construction.
        *   `dcalc_PW_precon(T, scale_T, P)`: Real. For sparse `T`, `P_ii = (27 + 18s + 12s^2 + 8s^3) / (27 + 18s + 12s^2 + 8s^3 + 16s^4)` where `s = 2*T_ii/scale_T`.
        *   `zcalc_PW_precon(T, scale_T, P)`: Complex. For sparse `T`, `P_ii = 1 / (1 + abs(T_ii)/scale_T)`.
    *   `die(message)`: Error termination subroutine, writes to "libOMM.err".

# Important Variables/Constants

*   `HAVE_CONFIG_H`, `HAVE_LAPACK`, `HAVE_SCALAPACK`, `HAVE_PSPBLAS`, `HAVE_MPI`: Preprocessor flags that enable use of corresponding libraries and code paths.
*   Input/output matrices are all of `TYPE(matrix)` from `MatrixSwitch`.
*   `m_operation`: Character(3) argument in many routines, passed to `MatrixSwitch` calls to select specific implementations.
*   `dp`, `cmplx_1`: Precision and complex unit constant (from `omm_params`).

# Usage Examples

These subroutines are designed for internal use by the `omm` and `omm_callback` procedures. For example, within the OMM CG loop:
1.  `calc_AW` (or `calc_HW_callback`) is used to form `HW = C^T*H*C` and `SW = C^T*S*C`.
2.  `calc_G` is called to get the gradient.
3.  `calc_coeff` is called to get quartic coefficients for the line search.
4.  `m_factorize`, `m_reduce`, etc., are used during setup phases depending on the OMM `flavour`.

# Dependencies and Interactions

*   **`MatrixSwitch`**: This is the primary dependency. All actual matrix manipulations (multiplications, additions, traces, allocations, LAPACK/ScaLAPACK calls for factorizations/inversions/solves) are performed through `MatrixSwitch` library calls.
*   **`omm_params`**: Used for fundamental constants like `dp`, `cmplx_1`, and for global parameters like `log_unit` and `mpi_rank` needed by the local `die` routine.
*   **LAPACK/ScaLAPACK**: These are runtime dependencies if `HAVE_LAPACK` or `HAVE_SCALAPACK` are defined. `omm_ops` contains the logic to call the appropriate routines (e.g., `dpotrf`, `pdpotrf`) via `MatrixSwitch`.
*   **pspBLAS**: If `HAVE_PSPBLAS` is defined, `calc_PW_precon` uses pspBLAS specific matrix components (`T%spm%dval`, `T%spm%zval`).
*   The module provides a clear separation of OMM-specific high-level operations from the general-purpose matrix operations in `MatrixSwitch`.
*   The `calc_HW_callback` with its `INTERFACE` for `H` allows for a flexible way to define the action of the Hamiltonian, potentially enabling matrix-free methods if the callback can compute `H*C` without explicitly forming `H`.
