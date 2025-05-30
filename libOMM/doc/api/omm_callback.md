# Overview

The file `libOMM/src/omm_callback.F90` contains the `omm_callback` subroutine, which is a variant of the main Orbital Minimization Method (OMM) routine found in `omm.F90`. The primary conceptual difference is its design to potentially accommodate a callback mechanism for applying the Hamiltonian operator (i.e., calculating `H*C` products), rather than solely relying on an explicit Hamiltonian matrix `H`. This is indicated by an `INTERFACE` block defining a procedure `H(C, HC)`.

However, in the current implementation within this specific file:
1.  The `omm_callback` subroutine itself still accepts an explicit `TYPE(matrix)` named `H` as an argument, just like `omm.F90`.
2.  The defined callback procedure `H(C,HC)` is not actually an argument to `omm_callback`.
3.  Operations involving the Hamiltonian, such as `H*C` or `H*G`, are performed by calling `calc_HW_callback` (from `omm_ops.F90`), which is presumably designed to use such a callback if it were passed down to it. The top-level `omm_callback` does not pass a specific Hamiltonian procedure to `calc_HW_callback`.

Apart from this structural difference related to Hamiltonian application, the `omm_callback` subroutine follows the same OMM algorithm as `omm.F90`, including CG minimization, quartic line search, state management for multiple spin/k-points, and support for different OMM "flavours" (though flavours 1 and 2 for Cholesky are commented out in this version's logic).

# Key Components

*   **Subroutine `omm_callback(...)`**: The OMM algorithm, largely mirroring `omm.F90`.
    *   **Input Arguments**: Mostly identical to those in `omm.F90`. This includes `m`, `n`, an explicit `TYPE(matrix) :: H`, `S`, `new_S`, `eta`, `C_min`, `init_C`, `T`, `scale_T`, `flavour` (note: logic for `flavour=1` and `flavour=2` is commented out), `np`, `ip`, `cg_tol`, `long_out`, `dealloc`, `m_storage`, `m_operation`.
    *   **Output Arguments**: `e_min` (OMM energy), `D_min` (density matrix or energy-weighted density matrix).
    *   **Defined Callback Interface (not a subroutine argument)**:
        ```fortran
        INTERFACE
          SUBROUTINE H(C, HC)
            USE MatrixSwitch
            IMPLICIT NONE
            TYPE(matrix), INTENT(IN) :: C    ! WF coeffs. matrix
            TYPE(matrix), INTENT(INOUT) :: HC ! work matrix (to store H*C)
          END SUBROUTINE H
        END INTERFACE
        ```
        This interface describes a procedure that would take a coefficient matrix `C` and compute `H*C`, returning it in `HC`.
    *   **Algorithm**: The sequence of operations (initialization, state management, CG loop, finalization) is the same as in `omm.F90`. The main difference in internal calls is the use of `calc_HW_callback` for operations like `H*C` and `H*G`, instead of `calc_AW` which is used for `S*C` and `S*G` (and for H-products in `omm.F90`).

# Important Variables/Constants

*   The set of important local and `SAVE`d module-level variables is identical to that in `omm.F90` (e.g., `HW(np)`, `SW(np)`, `G`, `D`, `coeff(0:4)` etc.).
*   The logic for `flavour` handling is modified:
    *   `flavour=0` (basic) and `flavour=3` (preconditioning with S, T optional) are handled.
    *   The code paths for `flavour=1` (Cholesky, S provided) and `flavour=2` (Cholesky, U provided) are commented out. This means these specific Cholesky modes are not functional in this version of the OMM routine.

# Usage Examples

The intended usage is similar to `omm.F90`, as part of an iterative electronic structure calculation. If the callback mechanism were fully plumbed through the argument list, a user would pass their own Hamiltonian application routine.
```fortran
! Conceptual call - NOTE: The actual omm_callback does not take a procedure for H.
! It still takes an explicit H matrix.
!
! EXTERNAL my_hamiltonian_action_subroutine ! Matches the defined H interface
!
! CALL omm_callback( &
!   m=size_m, n=size_n, H=H_matrix_explicit, S=S_matrix, ..., &
!   ! H_action_callback=my_hamiltonian_action_subroutine, ! This part is missing
!   ... other arguments ...
! )
```
In its current form, `H_matrix_explicit` (an actual `TYPE(matrix)`) is passed and used by `calc_HW_callback` (which itself might be designed to take a procedure, but isn't given one by `omm_callback`).

# Dependencies and Interactions

*   **`omm_ops`**: Key dependency. Uses:
    *   `calc_A`, `calc_A2` (for density matrix construction).
    *   `calc_AW` (for S-related products like `S*C`, `S*G`).
    *   `calc_HW_callback` (for H-related products like `H*C`, `H*G`). This is the main operational difference from `omm.F90`'s direct use of `calc_AW` for H.
    *   `calc_G` (gradient calculation).
    *   `calc_coeff` (quartic coefficients).
    *   `dp` (precision) and `die` (error handling).
*   **`MatrixSwitch`**: All underlying matrix operations are performed using this library.
*   **`omm_rand`**: For initializing wavefunction coefficients randomly.
*   **`MatrixSwitch_ops`**: For MPI related variables if `HAVE_MPI` is defined.
*   **`omm_quartic`** (implicitly): Dependency for `omm_solve_quartic` to solve the line search polynomial.
*   The primary functional difference from `omm.F90` is the call to `calc_HW_callback` instead of `calc_AW` for Hamiltonian products and the disabled Cholesky flavours (1 and 2). The callback interface `H` is defined but not used as a parameter to `omm_callback`.
