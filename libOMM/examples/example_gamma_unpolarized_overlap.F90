#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==================================================================================================!
! example1 : Gamma-point only, spin unpolarized                                                    !
!                                                                                                  !
! This example demonstrates a typical calculation with an outer loop of two MD iterations, and an  !
! inner loop of five SCF iterations per MD step. For convenience, the Hamiltonian and overlap      !
! matrices at each point have been pre-generated and are read in from file; in a real code,        !
! however, they should be calculated at each step based on the output density matrix given by      !
! libOMM.                                                                                          !
!                                                                                                  !
! This example is for a system with 832 basis orbitals and 128 occupied states (generated from a   !
! 64-atom supercell of bulk Si). At the end of each SCF iteration, the Kohn-Sham energy 2*e_min is !
! printed out, together with the first element of the density matrix as an extra check. At the end !
! of each MD iteration, the energy-weighted density matrix is also calculated, and its first       !
! element is printed out.                                                                          !
!                                                                                                  !
! Things to note:                                                                                  !
!   1. The eigenspectrum shift parameter eta has to be set larger than 0 for convergence, since    !
!      the occupied spectrum extends beyond 0 (this is typically a sign of a poor basis). Try      !
!      changing eta to see how the convergence speed is affected. However, also take into account: !
!   2. The Hamiltonian has to be provided to libOMM *already shifted by eta* (H -> H-eta*S)        !
!   3. There are no optional arguments in the call to libOMM. Therefore, matrices which are not    !
!      needed should simpy be passed without having been allocated by MatrixSwitch (m_allocate     !
!      routine). Other variables which are not needed will be ignored by libOMM.                   !
!   4. Try enabling Cholesky factorization (precon=1) or preconditioning (precon=3) in the call to !
!      libOMM to see how the convergence speed is affected. Preconditioning is even more effective !
!      if a T matrix is provided (scale_T should be set around 10 Ry). Take care that, for         !
!      Cholesky factorization, S and H will be overwritten by U and U^(-T)*H*U^(-1).               !
!   5. The dealloc variable should only be .true. for the very last call. This is because libOMM   !
!      stores and reuses internal information from one call to the next.                           !
!                                                                                                  !
! Sample output:                                                                                   !
!--------------------------------------------------------------------------------------------------!
! e_min :   -63.812664276045751                                                                    !
! D_11  :     0.500000000325070                                                                    !
! e_min :   -65.685566689065666                                                                    !
! D_11  :     0.500000000000009                                                                    !
! ED_11 :    -1.642139167231525                                                                    !
!--------------------------------------------------------------------------------------------------!
!==================================================================================================!
program example_gamma_unpolarized_overlap
  use MatrixSwitch

  implicit none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  real(dp), parameter :: e_min_check(2)=(/-63.812664_dp,-65.685567_dp/)
  real(dp), parameter :: D_el_check(2)=(/0.5_dp,0.5_dp/)
  real(dp), parameter :: ED_el_check=-1.642139_dp

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: new_S

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: num_pairs, m, n, i, j, k, f, flavour

  real(dp) :: he1, he2, se, e_min, el

  type(matrix) :: H, S, D_min, ED_min, C_min, T

  !**********************************************!

#ifdef HAVE_MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  call ms_scalapack_setup(mpi_comm_world,1,'c',3)

  m_storage='pddbc'
  m_operation='lap'
#else
  mpi_size=1
  mpi_rank=0

  m_storage='sdden'
  m_operation='lap'
#endif

  num_pairs=20

  m=num_pairs*2
  n=num_pairs

  do f=1,3

    select case(f)
    case (1)
      flavour=0
    case (2)
      flavour=1
    case(3)
      flavour=3
    end select

    call m_allocate(H,m,m,m_storage)
    call m_allocate(S,m,m,m_storage)
    call m_allocate(D_min,m,m,m_storage)
    call m_allocate(ED_min,m,m,m_storage)
    call m_allocate(C_min,n,m,m_storage)

    se=0.1_dp

    do i=1,2

      select case (i)
      case (1)
        he1=-1.5_dp
        he2=-3.0_dp
        new_S=.true.
      case (2)
        he1=-1.5_dp
        he2=-3.1_dp
        new_S=.false.
      end select

      call m_set(S,'a',0.0_dp,1.0_dp)

      call m_set_element(H,1,m,he1,0.0_dp)
      call m_set_element(H,1,2,he2,0.0_dp)
      call m_set_element(S,1,m,se,0.0_dp)
      call m_set_element(S,1,2,se,0.0_dp)
      do j=2,m-1
        if (mod(j,2)==0) then
          call m_set_element(H,j,j-1,he2,0.0_dp)
          call m_set_element(H,j,j+1,he1,0.0_dp)
        else
          call m_set_element(H,j,j-1,he1,0.0_dp)
          call m_set_element(H,j,j+1,he2,0.0_dp)
        end if
        call m_set_element(S,j,j-1,se,0.0_dp)
        call m_set_element(S,j,j+1,se,0.0_dp)
      end do
      call m_set_element(H,m,m-1,he2,0.0_dp)
      call m_set_element(H,m,1,he1,0.0_dp)
      call m_set_element(S,m,m-1,se,0.0_dp)
      call m_set_element(S,m,1,se,0.0_dp)

      call omm(m, &         ! m
               n, &         ! n
               H, &         ! H 
               S, &         ! S
               new_S, &     ! new_S
               e_min, &     ! e_min
               D_min, &     ! D_min
               .false., &   ! calc_ED
               0.0_dp, &    ! eta
               C_min, &     ! C_min
               .false., &   ! init_C
               T, &         ! T
               0.0_dp, &    ! scale_T
               flavour, &   ! flavour
               1, &         ! np
               1, &         ! ip
               -1.0_dp, &   ! cg_tol
               .true., &    ! long_out
               .false., &   ! dealloc
               m_storage, & ! m_storage
               m_operation) ! m_operation

      if (mpi_rank==0) print('(a,f21.15)'), 'e_min : ', e_min
      !call assert_equal_dp(e_min, e_min_check(i))

      call m_get_element(D_min,1,1,el)
      if (mpi_rank==0) print('(a,f21.15)'), 'D_11  : ', el
      !call assert_equal_dp(el, D_el_check(i))

    end do

    call omm(m, &         ! m
             n, &         ! n
             H, &         ! H 
             S, &         ! S
             .false., &   ! new_S
             e_min, &     ! e_min
             ED_min, &    ! D_min
             .true., &    ! calc_ED
             0.0_dp, &    ! eta
             C_min, &     ! C_min
             .false., &   ! init_C
             T, &         ! T
             0.0_dp, &    ! scale_T
             0, &         ! flavour
             1, &         ! np
             1, &         ! ip
             -1.0_dp, &   ! cg_tol
             .true., &    ! long_out
             .true., &    ! dealloc
             m_storage, & ! m_storage
             m_operation) ! m_operation

    call m_get_element(ED_min,1,1,el)
    if (mpi_rank==0) print('(a,f21.15)'), 'ED_11 : ', el
    !call assert_equal_dp(el, ED_el_check)

    call m_deallocate(C_min)
    call m_deallocate(ED_min)
    call m_deallocate(D_min)
    call m_deallocate(S)
    call m_deallocate(H)

  end do

#ifdef HAVE_MPI
  call mpi_finalize(mpi_err)
#endif

  contains

  subroutine assert_equal_dp(value1, value2)
    implicit none

    !**** PARAMS **********************************!

    real(dp), parameter :: tolerance=1.0d-5

    !**** INPUT ***********************************!

    real(dp), intent(in) :: value1
    real(dp), intent(in) :: value2

    !**********************************************!

    if (abs(value1-value2)>tolerance) stop 1

  end subroutine assert_equal_dp

end program example_gamma_unpolarized_overlap
