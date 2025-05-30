# Overview

The file `omm.F90` contains the primary subroutine `omm`, which is the core of the libOMM (Orbital Minimization Method) library. This routine implements a direct energy minimization scheme with respect to the wavefunction (WF) coefficients, represented as a matrix `C_min`. The minimization typically employs a conjugate gradient (CG) algorithm coupled with a quartic line search to find the optimal step size along search directions.

The `omm` subroutine is designed to be versatile, supporting several "flavours" of the OMM functional. These include:
*   Basic OMM.
*   OMM with Cholesky factorization of the overlap matrix `S` (where `S` itself or its Cholesky factor `U` can be provided).
*   Preconditioned OMM, where a preconditioner based on the overlap matrix `S` and optionally a kinetic energy matrix `T` can be used.

The routine is stateful, capable of managing calculations for multiple spin and k-points (`np`, `ip`) by saving intermediate results. It can be instructed to calculate either the standard density matrix `D_min` or an energy-weighted density matrix. It's intended to be called iteratively, such as within a Self-Consistent Field (SCF) cycle in electronic structure calculations.

# Key Components

*   **Subroutine `omm(...)`**: The main computational routine of the libOMM library.
    *   **Input Arguments**:
        *   `m, n`: Integers, size of the atomic orbital basis and the number of occupied states (or WFs), respectively.
        *   `H`: `TYPE(matrix)`, the Hamiltonian matrix (m x m).
        *   `S`: `TYPE(matrix)`, the overlap matrix (m x m). Can also be its Cholesky factor `U` if `flavour=2`. Can be uninitialized if `no_S` is true.
        *   `new_S`: Logical. If `.TRUE.`, the `S` matrix (or `T` for preconditioner) is new for the current `ip` and needs reprocessing (e.g., Cholesky factorization, preconditioner construction).
        *   `calc_ED`: Logical. If `.TRUE.`, the routine computes the energy-weighted density matrix using existing WF coefficients and returns.
        *   `eta`: Real(dp), an energy shift parameter used in the OMM functional and energy-weighted density calculation.
        *   `C_min`: `TYPE(matrix)` (inout), the wavefunction coefficient matrix (n x m). Input if `init_C` is true or from previous calls; output with optimized coefficients.
        *   `init_C`: Logical. If `.TRUE.`, `C_min` is assumed to be initialized externally and will be used as the starting point.
        *   `T`: `TYPE(matrix)`, the kinetic energy matrix (m x m), used if `flavour=3` or `flavour=4` for preconditioning.
        *   `scale_T`: Real(dp), a scale factor for the kinetic energy matrix in the preconditioner.
        *   `flavour`: Integer, specifies the OMM variant:
            *   `0`: Basic OMM (S is used as is, or identity if `no_S`).
            *   `1`: Cholesky factorization, `S` provided (S is factorized to U).
            *   `2`: Cholesky factorization, `U` provided (S is assumed to be U).
            *   `3`: Preconditioning, `S` provided (T optional). Preconditioner `P = (S + T/scale_T)^-1`.
            *   `4`: Preconditioning, `S` not used directly for P. `P` calculated via `calc_PW_precon` (likely for plane-wave specific preconditioner using T).
        *   `np`: Integer, total number of independent calculations (e.g., spin points * k-points).
        *   `ip`: Integer, the current calculation index (1 to `np`).
        *   `cg_tol`: Real(dp), convergence tolerance for the CG minimization based on relative energy change.
        *   `long_out`: Logical. If `.TRUE.`, detailed iteration information is printed to "libOMM.log".
        *   `dealloc`: Logical. If `.TRUE.`, all internal `SAVE`d matrices associated with `ip` (or all `np` points if called after the loop) are deallocated.
        *   `m_storage`: Character(5), MatrixSwitch label for storage format of internal matrices.
        *   `m_operation`: Character(3), MatrixSwitch label for implementation of operations.
    *   **Output Arguments**:
        *   `e_min`: Real(dp), the minimized OMM functional energy (spin degeneracy not included).
        *   `D_min`: `TYPE(matrix)` (inout), the resulting density matrix `C_min * (2*I - S*C_min^T*C_min) * C_min^T` or energy-weighted density matrix if `calc_ED` is true.

    *   **Core Algorithm**:
        1.  **Setup**: Handles MPI information, logging, OMM flavour logic, and CG tolerance.
        2.  **Energy-Weighted Density (Early Exit)**: If `calc_ED` is true, computes `D_min = C*[(2*I-SW)*(HW+eta*SW)]*C^T` using data from the previous successful minimization for `ip`, then returns.
        3.  **State Management**: Allocates/accesses `SAVE`d arrays for storing matrices like `HW(ip)`, `SW(ip)`, `C_Chl(ip)`, `P(ip)` etc., which persist between calls for a given `ip`.
        4.  **Preprocessing**:
            *   If Cholesky flavour (`1,2`), performs Cholesky factorization `S = U^T*U` (if `S` is new and `flavour=1`) and transforms `H` to `U^-T*H*U^-1`.
            *   If preconditioning flavour (`3,4`) and `S` (or `T`) is new, constructs the preconditioner matrix `P(ip)`.
        5.  **WF Initialization**: If it's the first call for `ip` and `C_min` is not externally initialized, `C_min` is filled with random values and scaled. If Cholesky flavour, `C_min` is transformed to `C_Chl(ip) = C_min * U^-1`.
        6.  **Initial Calculations**: Computes initial `HW = C^T*H*C`, `SW = C^T*S*C`, gradient `G`, preconditioned gradient `PG` (if used), and terms for the quartic line search (`HWd`, `SWd`, `HWdd`, `SWdd`). The initial energy `e_min` is the 0th order coefficient from `calc_coeff`.
        7.  **CG Loop**: Iteratively refines `C_min` (or `C_Chl(ip)`).
            *   Determines search direction `D` (from `G` or `PG`, and previous `D`).
            *   Calculates coefficients for a quartic polynomial describing energy along `D`.
            *   Solves the quartic equation for optimal step `x_min(ip)` using `omm_solve_quartic`.
            *   Handles line search failures by rescaling `C_min` and restarting the CG direction.
            *   Updates `C_min = C_min + x_min(ip)*D`.
            *   Updates `HW`, `SW`, `G`, `PG`, and `lambda` (Polak-Ribiere factor).
            *   Checks for convergence based on `e_diff` (relative energy change).
        8.  **Post-CG**:
            *   Calculates `QW(ip) = 2*I - SW(ip)`.
            *   If Cholesky, back-transforms `C_min` from `C_Chl(ip)`.
            *   Computes the final density matrix `D_min = C_min * QW(ip) * C_min^T`.
            *   Calculates `Tr[QW*SW]` for an occupancy check.
            *   Adjusts `e_min` with the `eta` term: `e_min = e_min + Tr[QW*SW]*eta`.
        9.  **Cleanup**: Deallocates temporary work matrices. If `dealloc` is true, deallocates persistent `SAVE`d arrays.

# Important Variables/Constants

*   **Stateful `SAVE`d Arrays**: `HW(np)`, `SW(np)`, `SC(np)`, `SG(np)`, `C_Chl(np)`, `P(np)`, `QW(np)`, `CD(np)`, `x_min(np)`, `first_call(np)`. These arrays store matrix data and state for each of the `np` spin/k-points across multiple calls to the `omm` subroutine. This is critical for iterative procedures like SCF.
*   `log_unit`: File unit for "libOMM.log".
*   `mpi_size`, `mpi_rank`: MPI variables (from `MatrixSwitch_ops` if `HAVE_MPI`).
*   Local `TYPE(matrix)` variables: `SWd`, `SWdd`, `G`, `G_p`, `PG`, `PG_p`, `D`, `HC`, `HG`, `work1`, `work2(ip)` (though `work2` is also in `SAVE`d arrays). These are used for intermediate steps in the CG algorithm.
*   `coeff(0:4)`: Stores coefficients of the quartic polynomial for the line search.

# Usage Examples

The `omm` subroutine is the main computational engine of libOMM. It's intended to be called repeatedly within an iterative framework, like an SCF cycle in quantum chemistry/physics codes.

```fortran
! Simplified conceptual usage in an SCF loop
DO iter = 1, max_scf_iterations
  DO ispin_kpoint = 1, num_total_spin_kpoints
    ! Prepare H_current, S_current, C_coeffs_guess for ispin_kpoint
    ! Set flags like is_S_new, is_C_init based on SCF state

    CALL omm( &
      m=basis_size, n=num_occupied_states, &
      H=H_current, S=S_current, new_S=is_S_new, &
      e_min=energy_omm, D_min=density_matrix_out, calc_ED=.FALSE., &
      eta=energy_shift_eta, C_min=C_coeffs_guess, init_C=is_C_init, &
      T=T_kinetic_matrix, scale_T=precon_scale_factor, &
      flavour=omm_calculation_flavour, &
      np=num_total_spin_kpoints, ip=ispin_kpoint, &
      cg_tol=my_cg_tolerance, long_out=.TRUE., dealloc=.FALSE., &
      m_storage="sdden", m_operation="ref" &
    )
    ! Use density_matrix_out to update H_current for next SCF iteration
    ! Store C_coeffs_guess for next iteration or for calc_ED
  END DO
  ! Check for SCF convergence
  IF (scf_converged) EXIT
END DO

! Optional: Final calculation of energy-weighted density matrix
! CALL omm(..., calc_ED=.TRUE., ...)

! Optional: Cleanup OMM's internal saved matrices
! CALL omm(..., dealloc=.TRUE., ...)
! (Need to call for one ip, or ensure a final call that triggers deallocation for all ip if intended)
```

# Dependencies and Interactions

*   **`omm_ops`**: This is a major dependency. `omm` calls many subroutines from `omm_ops` for specific mathematical operations:
    *   `calc_A`, `calc_A2`: For constructing density-like matrices.
    *   `calc_AW`: For `C^T*M*C` type operations.
    *   `calc_G`: For calculating the OMM gradient.
    *   `calc_coeff`: For calculating coefficients of the quartic line search.
    *   `calc_PW_precon`: For specific preconditioner calculation (flavour 4).
    *   It also uses `dp` (double precision kind) and `die` (error handling) from `omm_ops` (or `MatrixSwitch_ops`).
*   **`MatrixSwitch`**: All low-level matrix manipulations (allocation, deallocation, arithmetic operations like `m_add`, `mm_multiply`, `m_scale`, `m_set_element`, and more complex operations like `m_factorize`, `m_reduce`, `m_inverse`, `m_transform`, `m_back_transform`, `mm_trace`) are delegated to the `MatrixSwitch` library.
*   **`omm_rand`**: Uses `omm_rand_seed()` to initialize a seed and `omm_bsd_lcg()` to generate random numbers for the initial guess of `C_min` when `init_C` is false and it's the `first_call` for that `ip`.
*   **`omm_quartic`** (external/not shown in `USE`): The line `call omm_solve_quartic(coeff(0:4),x_min(ip),ls_fail)` indicates a dependency on a subroutine `omm_solve_quartic` which is responsible for finding the roots of the quartic polynomial and determining the line search minimum. This routine is likely located in `omm_quartic.F90`.
*   **MPI**: If `HAVE_MPI` is defined, it uses `ms_mpi_size` and `ms_mpi_rank` (presumably from `MatrixSwitch_ops` or a similar MPI utility module accessed via `MatrixSwitch`).
*   **File I/O**: Writes logs to "libOMM.log". The `log_unit` is likely defined in `omm_ops` or `omm_params`.
*   **Stateful Nature**: The use of `SAVE` for several key allocatable arrays (`HW`, `SW`, `SC`, `C_Chl`, `P`, `QW`, `CD`, `x_min`, `first_call`) indexed by `ip` makes the `omm` subroutine stateful between calls for the same `ip`. This is essential for the iterative SCF process. The `dealloc` flag provides a way to reset this state.
