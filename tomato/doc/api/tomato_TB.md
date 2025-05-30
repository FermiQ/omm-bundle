# Overview

The `tomato_TB.F90` file contains a single, comprehensive subroutine named `tomato_TB`. This subroutine is designed to construct Tight-Binding (TB) Hamiltonian (`H`) and, optionally, Overlap (`S`) matrices for a given physical system. The construction process is template-based, meaning it relies on pre-calculated data for a base system, which is then scaled, repeated, and modified according to user-specified parameters to model different system sizes, basis sets, boundary conditions (k-points), and the presence of defects.

The routine offers significant flexibility through various input switches and parameters, allowing users to:
*   Define the system based on a `system_label` and a `template_basedir` where template files are stored.
*   Control the basis set size by either specifying a desired fractional occupancy of states or by targeting a specific number of orbitals per atom. The routine selects the closest available basis from the template.
*   Define the supercell size by specifying the number of unit cell repetitions in three directions, or by targeting a total number of orbitals.
*   Adjust matrix sparsity by selecting an orbital interaction cutoff radius, or by targeting a specific sparsity percentage.
*   Perform Gamma-point calculations (real matrices) or k-point specific calculations (complex matrices) for periodic systems.
*   Optionally introduce a point defect by applying a perturbation to the on-site Hamiltonian elements of a specific atom.
*   Choose the storage format for the output matrices `H` and `S` via the `MatrixSwitch` library.

The subroutine reads system parameters and matrix element templates from external files, processes this information, and then constructs the `H` and `S` matrices using the `MatrixSwitch` library for matrix operations. It also includes an internal random number generator for estimating the sparsity of the generated matrices.

# Key Components

*   **Subroutine `tomato_TB(...)`**:
    *   **Input Arguments**:
        *   `template_basedir`: `CHARACTER(*)`, path to the directory containing template files.
        *   `system_label`: `CHARACTER(*)`, label identifying the system (used in template filenames).
        *   `switch1, switch2, switch3`: `LOGICAL` flags that control how `num_orbs_per_atom`, `num_cells_dir`, and `orb_r_cut` are determined. If true, the routine tries to match `frac_occ`, `num_orbs`, or `sparsity` respectively. If false, it uses the provided `num_orbs_per_atom`, `num_cells_dir`, or `orb_r_cut`.
        *   `frac_occ`: `REAL(dp)` (inout), desired fractional occupancy (num_occ_states / num_orbs). Used if `switch1` is true. Updated by the routine.
        *   `num_orbs_per_atom`: `INTEGER` (inout), desired/actual number of orbitals per atom.
        *   `num_orbs`: `INTEGER` (inout), desired/actual total number of orbitals in the supercell.
        *   `num_cells_dir(3)`: `INTEGER` (inout), number of unit cell repetitions in x, y, z directions.
        *   `sparsity`: `REAL(dp)` (inout), desired/actual matrix sparsity.
        *   `orb_r_cut`: `REAL(dp)` (inout), orbital interaction cutoff radius.
        *   `k_point(3)`: `REAL(dp)` (inout), k-point coordinates for reciprocal space sampling (fractional coordinates). Adjusted for non-periodic directions.
        *   `gamma_point`: `LOGICAL`, if `.TRUE.`, forces a Gamma-point calculation (k_point = [0,0,0]), resulting in real matrices.
        *   `defect`: `LOGICAL`, if `.TRUE.`, applies a defect perturbation.
        *   `defect_perturbation`: `REAL(dp)`, the scaling factor for the defect perturbation on Hamiltonian elements.
        *   `m_storage`: `CHARACTER(5)`, specifies the `MatrixSwitch` storage format (e.g., "sdden", "sdcsc").
        *   `build_matrix`: `LOGICAL`, if `.TRUE.`, the H and S matrices are constructed. If `.FALSE.`, the routine only calculates and outputs parameters like `num_orbs`, `num_occ_states`, and `sparsity`.
    *   **Output Arguments**:
        *   `num_occ_states`: `INTEGER`, total number of occupied states calculated based on `num_occ_states_per_atom` and supercell size.
    *   **Inout Arguments**:
        *   `H, S`: `TYPE(matrix)` (from `MatrixSwitch` library), the constructed Hamiltonian and Overlap matrices.
    *   **Internal Random Number Generator**:
        *   `omm_rand_seed()`: A contained function to generate a seed. It uses system time unless the `NORAND` preprocessor macro is defined (in which case it's fixed).
        *   `omm_bsd_lcg(x, r)`: A contained subroutine implementing a BSD Linear Congruential Generator. Used for the final sparsity estimation step.
    *   **Local Derived Type `t_basis_subset`**:
        *   `ref`: `INTEGER`, reference number of orbitals for a full basis set.
        *   `index(:)`: `INTEGER, ALLOCATABLE`, array mapping subset orbital indices to full basis indices.
    *   **Core Functionality**:
        1.  **MPI Setup**: Initializes local MPI rank and size variables if `HAVE_MPI` is defined (using values from `MatrixSwitch_ops`).
        2.  **Read Template Information**: Reads a `*_template.info` file. This file contains metadata about the system: number of periodic dimensions, atoms per cell, occupied states per atom, and lists of available orbital cutoff radii and basis sizes (including information about whether a basis is a subset of a larger one). This data is broadcast to all MPI processes.
        3.  **Determine Basis Size**: Based on `switch1` and input `frac_occ` or `num_orbs_per_atom`, the routine selects the most appropriate `num_orbs_per_atom` from the available template basis sizes. The actual `frac_occ` is then updated. It also determines if the selected basis is a subset of a reference basis.
        4.  **Determine Supercell Dimensions**: Based on `switch2` and input `num_orbs` or `num_cells_dir`, the routine determines the number of unit cell repetitions in each direction (`num_cells_dir`). It then calculates the total number of cells, atoms, orbitals (`num_orbs`), and occupied states (`num_occ_states`).
        5.  **Determine Orbital Cutoff and Sparsity**: Based on `switch3` and input `sparsity` or `orb_r_cut`, it selects the most appropriate `orb_r_cut_i` from the template's cutoff list. If `build_matrix` is false, it estimates the sparsity based on the chosen parameters by reading the number of non-zero elements from the corresponding template data file.
        6.  **K-Point Adjustment**: If the system is not periodic in certain directions, the corresponding `k_point` components are set to zero.
        7.  **Matrix Construction (if `build_matrix` is `.TRUE.`)**:
            *   Reads the main template data file (e.g., `system_label_cutoff_basis.template`), which contains non-zero Hamiltonian (and optionally overlap) matrix elements along with their orbital indices and relative cell offsets `(i, j, k)` for interactions.
            *   Allocates `H` and (if `overlap_present` in template) `S` matrices using `MatrixSwitch%m_allocate` with the specified `m_storage` format.
            *   Creates a mapping `cell_index` from 3D cell coordinates to a base orbital index for that cell.
            *   Iterates through each non-zero element entry in the template file. For each template entry:
                *   It iterates through all cells `(i,j,k)` in the defined supercell.
                *   It calculates the target cell `(i2,j2,k2)` interacting with the current cell `(i,j,k)` based on the relative cell offsets from the template entry. Periodic boundary conditions are applied using `modulo`.
                *   The global orbital indices `a` (in cell `(i,j,k)`) and `b` (in cell `(i2,j2,k2)`) are computed. If a basis subset is used, a `subset_convert` mapping is applied to the template orbital indices.
                *   If `gamma_point` is true, real matrix elements are set in `H` (and `S`) using `MatrixSwitch%m_set_element`.
                *   If not `gamma_point` (k-point calculation), a phase factor `exp(I * k.T)` is calculated based on `k_point` and the cell offset vector `(i3,j3,k3)` (derived from how many times the interaction crossed a supercell boundary). Complex matrix elements (template value * phase factor) are set in `H` (and `S`).
            *   **Defect Application**: If `defect` is true, it iterates through the template entries again. If an entry corresponds to an on-site interaction within the first atom of the primary cell (cell 0,0,0), its Hamiltonian value is perturbed by `defect_perturbation`.
            *   **Sparsity Estimation**: After construction, it estimates the sparsity of matrix `H` by randomly sampling a number of its elements (equal to `num_orbs`) using the internal LCG and `MatrixSwitch%m_get_element`.
        8.  **Cleanup**: Deallocates internal temporary arrays.

# Important Variables/Constants

*   `dp`, `i64`: Integer parameters for double precision real and 64-bit integer kinds.
*   `Pi`, `twoPi`, `cmplx_1`, `cmplx_i`, `cmplx_0`: Standard mathematical constants.
*   `HAVE_MPI`, `NORAND`: Preprocessor macros. `NORAND` makes the random seed fixed.
*   File I/O: The routine reads two types of template files:
    *   `[template_basedir]/[system_label]_template.info`: Contains metadata like periodicity, atom counts, and available basis/cutoff options.
    *   `[template_basedir]/[system_label]_[orb_r_cut_i]_[basis_ref].template`: Contains the actual tight-binding parameters (orbital indices, cell offsets, H and S values).
*   `basis_subset(:)`: Array of `TYPE(t_basis_subset)` storing information about how different basis sets might be subsets of larger ones.
*   `cell_index(:,:, :)`: A 3D array mapping supercell cell coordinates `(i,j,k)` to a base global orbital index for that cell.
*   `subset_convert(:)`: A temporary mapping array used when a selected basis is a subset of the reference basis from the template.

# Usage Examples

The `tomato_TB` subroutine is a high-level tool for generating tight-binding Hamiltonians and overlap matrices.
```fortran
! Assuming tomato_TB is accessible, perhaps via a module
USE tomato_module ! Hypothetical module containing tomato_TB
USE MatrixSwitch  ! For TYPE(matrix)
IMPLICIT NONE

CHARACTER(LEN=256) :: template_dir = "/path/to/my/templates"
CHARACTER(LEN=50)  :: system_name  = "graphene_primitive"
TYPE(matrix)       :: hamiltonian_matrix, overlap_matrix
INTEGER            :: num_orbitals, num_occupied, num_cells_xyz(3)
REAL(dp)           :: k_vec(3), target_sparsity

! Initialize parameters
num_cells_xyz = [10, 10, 1] ! Example: 10x10x1 supercell
k_vec         = [0.0_dp, 0.0_dp, 0.0_dp] ! Gamma point calculation
target_sparsity = 0.95_dp

CALL tomato_TB(template_basedir=TRIM(template_dir), &
               system_label=TRIM(system_name), &
               switch1=.FALSE., num_orbs_per_atom=1, & ! Use 1 orb/atom (e.g., pz for graphene)
               switch2=.FALSE., num_cells_dir=num_cells_xyz, &
               switch3=.TRUE., sparsity=target_sparsity, & ! Target a specific sparsity
               num_occ_states=num_occupied, &
               gamma_point=(K_VEC(1)==0.0_dp .AND. K_VEC(2)==0.0_dp .AND. K_VEC(3)==0.0_dp), &
               k_point=k_vec, &
               defect=.FALSE., defect_perturbation=0.0_dp, &
               H=hamiltonian_matrix, S=overlap_matrix, m_storage="sdcsc", & ! Store as serial dense CSC
               build_matrix=.TRUE., &
               num_orbs=num_orbitals, frac_occ=target_frac_occ, orb_r_cut=target_orb_cut)

PRINT *, "Total orbitals:", num_orbitals
PRINT *, "Total occupied states:", num_occupied
PRINT *, "Final sparsity:", target_sparsity
PRINT *, "Final frac_occ:", target_frac_occ
PRINT *, "Final orb_r_cut:", target_orb_cut

! ... use H_matrix and S_matrix ...

CALL m_deallocate(hamiltonian_matrix)
CALL m_deallocate(overlap_matrix)
```

# Dependencies and Interactions

*   **`MatrixSwitch`**: This is a primary dependency. The subroutine uses `MatrixSwitch` for allocating (`m_allocate`) the `H` and `S` matrices and for setting their elements (`m_set_element`). It also uses `m_get_element` for the final sparsity estimation.
*   **`MatrixSwitch_ops`** (Conditional): If `HAVE_MPI` is defined, this module (part of `MatrixSwitch`) is used to access MPI communicator details like `ms_mpi_comm`, `ms_mpi_size`, and `ms_mpi_rank`.
*   **MPI Library (`mpif.h`)**: If `HAVE_MPI` is defined, MPI routines like `MPI_BCAST` are used to distribute data read by the root process (rank 0) from template files to all other processes.
*   **File System**: The routine depends on the presence and correct formatting of template files (`*_template.info` and `_cutoff_basis.template`) in the specified `template_basedir`.
*   The internal random number generator (`omm_rand_seed`, `omm_bsd_lcg`) is defined locally within the subroutine and does not depend on external RNG modules.
*   The logic for constructing matrices involves careful handling of 1D and 3D indexing, periodic boundary conditions (phase factors for k-points), and potentially mapping orbital indices if a subset of a larger template basis is used.
