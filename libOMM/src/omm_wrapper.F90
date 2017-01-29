#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine omm_wrapper(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,&
                       m_storage,m_operation)
  use omm_params, only : dp
  use MatrixSwitch_wrapper
  use MatrixSwitch_wrapper_params, only : ms_matrices, ms_lookup
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use
  character(*), intent(in) :: H ! Hamiltonian matrix
  character(*), intent(in) :: S ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  character(*), intent(in) :: D_min ! density (or energy-weighted density) matrix
  character(*), intent(in) :: C_min ! WF coeffs.
  character(*), intent(in) :: T ! kinetic energy matrix

#ifdef CBIND
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** LOCAL ***********************************!

  logical :: new_S_conv
  logical :: calc_ED_conv
  logical :: init_C_conv
  logical :: long_out_conv
  logical :: dealloc_conv

  !**********************************************!

  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc

  call omm(m, &
           n, &
           ms_matrices(ms_lookup(H)), &
           ms_matrices(ms_lookup(S)), &
           new_S_conv, &
           e_min, &
           ms_matrices(ms_lookup(D_min)), &
           calc_ED_conv, &
           eta, &
           ms_matrices(ms_lookup(C_min)), &
           init_C_conv, &
           ms_matrices(ms_lookup(T)), &
           scale_T, &
           flavour, &
           np, &
           ip, &
           cg_tol, &
           long_out_conv, &
           dealloc_conv, &
           m_storage, &
           m_operation)

end subroutine omm_wrapper
