# Overview

The `omm_quartic.F90` file provides subroutines related to solving a quartic polynomial equation, which is a critical step in the line search minimization phase of the Orbital Minimization Method (OMM) algorithm. When minimizing the OMM energy functional along a given search direction, the energy can be expressed or approximated as a quartic function of the step length lambda (`E(lambda) = c4*lambda^4 + c3*lambda^3 + c2*lambda^2 + c1*lambda + c0`). Finding the minimum of this quartic involves finding the roots of its derivative, which is a cubic polynomial.

This file contains:
1.  `omm_fit_quartic`: A routine to determine the coefficients of a quartic polynomial by fitting to three points and two gradients at those points. (Note: The main OMM routines in `omm.F90` and `omm_callback.F90` use an analytical method `calc_coeff` based on matrix traces, not this fitting routine, to get the quartic coefficients.)
2.  `omm_solve_quartic`: A routine that takes the coefficients of a quartic polynomial and finds the step `x_min` (lambda) that minimizes it. This involves solving the derivative (a cubic equation) analytically.

# Key Components

*   **Module `omm_quartic`** (implicitly, as the file contains these procedures):
    *   **Subroutine `omm_fit_quartic(x, y, g, c)`**:
        *   **Purpose**: Calculates the five coefficients `c(0:4)` of a quartic polynomial `y(x) = c(4)x^4 + c(3)x^3 + c(2)x^2 + c(1)x + c(0)`.
        *   **Inputs**:
            *   `x(1:3)`: A `REAL(dp)` array containing the x-coordinates of three distinct points.
            *   `y(1:3)`: A `REAL(dp)` array containing the y-values (function values) at the corresponding `x` points.
            *   `g(1:2)`: A `REAL(dp)` array containing the gradients (dy/dx) at the first two points, `x(1)` and `x(2)`.
        *   **Output**:
            *   `c(0:4)`: A `REAL(dp)` array storing the calculated coefficients of the quartic polynomial.
        *   **Method**: The coefficients are determined by solving a system of five linear equations derived from the five input conditions (3 function values and 2 gradient values). The explicit formulas in the code are stated to be "produced automatically using Maple 12."

    *   **Subroutine `omm_solve_quartic(c, x_min, fail)`**:
        *   **Purpose**: Finds the value `x_min` that minimizes the quartic polynomial defined by coefficients `c`.
        *   **Inputs**:
            *   `c(0:4)`: A `REAL(dp)` array of the quartic coefficients.
        *   **Outputs**:
            *   `x_min`: A `REAL(dp)` scalar, the step length that minimizes the quartic.
            *   `fail`: A logical flag. It is set to `.TRUE.` if a suitable minimum cannot be found (e.g., if the quartic is unbounded below in the search direction relevant to OMM, or if other conditions indicate failure).
        *   **Method**:
            1.  **Derivative**: The minimum of the quartic `E(x)` occurs where its derivative `dE/dx = 4*c(4)*x^3 + 3*c(3)*x^2 + 2*c(2)*x + c(1)` is zero. This is a cubic equation.
            2.  **Special Case (Quadratic Derivative)**: If `c(4)` is very small (or `c(2)` very large, making the quartic effectively quadratic in its derivative's dominant terms for finding minima), it approximates the minimum as `x_min = -0.5*c(1)/c(2)` (minimum of a quadratic `c(2)x^2 + c(1)x + c(0)`).
            3.  **Cubic Solver**: Otherwise, it solves the cubic equation for its real roots `t(1:3)` using the general analytical solution for cubic equations (Cardano's method or similar, as described by Numerical Recipes). This involves:
                *   Transforming to a depressed cubic form.
                *   Calculating discriminants `Q` and `R`.
                *   If `R^2 < Q^3` (three distinct real roots): Calculate roots using trigonometric form (`cos(theta/3)`).
                *   If `R^2 >= Q^3` (one real root, possibly multiple if `R^2=Q^3`): Calculate root(s) using `S` and `U` terms.
            4.  **Minimum Selection**:
                *   Evaluates the quartic function `c(4)*t^4 + ... + c(0)` at the real roots `t`.
                *   If `c(4) > 0` (quartic opens upwards, desired for minimization): Implements specific logic to choose among the roots, often preferring a positive root that leads to a valid energy decrease. The selection logic (`x_order`) aims to find the most relevant minimum for the OMM line search.
                *   If `c(4) <= 0` (quartic opens downwards or is of lower order): Selects the root that gives the lowest function value. If `c(4) < 0` and only one real root for the derivative, it may set `fail = .TRUE.`.

# Important Variables/Constants

*   `dp`: The double precision kind parameter, imported from `omm_params`.
*   `Pi`: The mathematical constant Pi, imported from `omm_params`, used in the trigonometric solution for cubic roots.
*   In `omm_fit_quartic`:
    *   `x, y, g`: Input arrays defining the fitting conditions.
    *   `c`: Output array of quartic coefficients.
*   In `omm_solve_quartic`:
    *   `c`: Input array of quartic coefficients.
    *   `x_min`, `fail`: Output variables for the minimum position and success status.
    *   `a, b, d`: Coefficients of the normalized cubic derivative equation.
    *   `Q, R, theta, S, U`: Intermediate values in the analytical cubic root-finding formulas.
    *   `t(1:3)`: Array to store the real roots of the cubic derivative.
    *   `z(1:3)`: Array to store the quartic function's values at the roots `t`.

# Usage Examples

These subroutines are primarily internal tools for the main OMM algorithms (`omm.F90`, `omm_callback.F90`).
*   `omm_fit_quartic`: This routine is *not* directly used by the main `omm.F90` or `omm_callback.F90` as they use an analytical method (`calc_coeff`) to derive the quartic coefficients. `omm_fit_quartic` would be useful in alternative line search schemes that rely on fitting to numerically evaluated points.
*   `omm_solve_quartic`: This is a core part of the line search in `omm.F90` and `omm_callback.F90`. After `calc_coeff` determines the quartic coefficients `c`, `omm_solve_quartic` is called:
    ```fortran
    ! Inside omm.F90 or omm_callback.F90
    REAL(dp) :: quartic_coeffs(0:4)
    REAL(dp) :: step_length
    LOGICAL  :: linesearch_has_failed

    ! ... (quartic_coeffs are computed by calc_coeff) ...

    CALL omm_solve_quartic(quartic_coeffs, step_length, linesearch_has_failed)

    IF (linesearch_has_failed) THEN
      ! Handle failure: e.g., try a different step or rescale matrices
    ELSE
      ! Update wavefunction coefficients using step_length
    END IF
    ```

# Dependencies and Interactions

*   **`omm_params`**: Relies on this module for the `dp` kind parameter (double precision) and the constant `Pi`.
*   The `omm_fit_quartic` subroutine is mathematically self-contained once its input data (points and gradients) are provided.
*   The `omm_solve_quartic` subroutine is a critical component of the line minimization step within the conjugate gradient algorithm employed by the main OMM routines. The success and accuracy of `omm_solve_quartic` directly impact the convergence and stability of the OMM energy minimization.
*   The selection logic for `x_min` when multiple roots exist for the derivative is tailored to the context of an energy minimization, typically seeking a step that leads to the lowest possible energy within the quartic approximation. The `fail` flag is important for the calling routine to handle problematic line search steps.
