#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==================================================================================================!
! example : Gamma-point only program with real matrices in simple dense (serial)/dense block       !
!           cyclic (parallel) format                                                               !
!                                                                                                  !
! This example demonstrates how to calculate the first two iterations in an SCF loop using a toy   !
! system. At the end of the loop, the energy-weighted density matrix is build from the current     !
! coefficients matrix.                                                                             !
!                                                                                                  !
! The calculation is repeat three times with different flavours of the OMM functional:             !
!   1. basic                                                                                       !
!   2. Cholesky factorization with S matrix provided                                               !
!   3. preconditioning                                                                             !
! All three give the same output, but require different numbers of line searches to reach          !
! convergence.                                                                                     !
!                                                                                                  !
! Sample output:                                                                                   !
!--------------------------------------------------------------------------------------------------!
! e_min [1] :   -24.438247165207933                                                                !
! D_11 [1]  :     0.375373030916752                                                                !
! e_min [2] :   -24.458492988091574                                                                !
! D_11 [2]  :     0.376073904787868                                                                !
! ED_11     :    -0.432150986741063                                                                !
! e_min [1] :   -24.438247214944731                                                                !
! D_11 [1]  :     0.375337656047806                                                                !
! e_min [2] :   -24.458493047887327                                                                !
! D_11 [2]  :     0.376091689034213                                                                !
! ED_11     :    -0.432152431346589                                                                !
! e_min [1] :   -24.438247217082530                                                                !
! D_11 [1]  :     0.375338986730707                                                                !
! e_min [2] :   -24.458493048136315                                                                !
! D_11 [2]  :     0.376091631495715                                                                !
! ED_11     :    -0.432151505859273                                                                !
!--------------------------------------------------------------------------------------------------!
!==================================================================================================!
program example_gamma
  use MatrixSwitch

  implicit none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  real(dp), parameter :: e_min_check(2)=(/-24.43825_dp,-24.45849_dp/)
  real(dp), parameter :: D_el_check(2)=(/0.37534_dp,0.37608_dp/)
  real(dp), parameter :: ED_el_check=-0.43215_dp

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: new_S, init_C

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: m, n, num_elements, i, j, k, f, flavour
  integer, allocatable :: elements_indices(:,:)

  real(dp) :: e_min, el
  real(dp), allocatable :: elements_values(:,:)

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

  ! Read input data
  if (mpi_rank==0) then
    open(10,file=SRCDIR//'/example_gamma.dat')
    read(10,'(i2,1x,i2,1x,i4)') m, n, num_elements
    allocate(elements_indices(2,num_elements))
    allocate(elements_values(2,num_elements))
    do i=1,num_elements
      read(10,'(i2,1x,i2,2(1x,es10.3e2))') elements_indices(1:2,i), elements_values(1:2,i)
    end do
    close(10)
  end if
#ifdef HAVE_MPI
  call mpi_bcast(m,1,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(n,1,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(num_elements,1,mpi_int,0,mpi_comm_world,mpi_err)
  if (mpi_rank/=0) then
    allocate(elements_indices(2,num_elements))
    allocate(elements_values(2,num_elements))
  end if
  call mpi_bcast(elements_indices,2*num_elements,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(elements_values,2*num_elements,mpi_double_precision,0,mpi_comm_world,mpi_err)
#endif

  ! Repeat test with different OMM flavours
  do f=1,3

    select case(f)
    case (1)
      ! Basic
      flavour=0
    case (2)
      ! Cholesky factorization with S matrix provided
      flavour=1
    case(3)
      ! Preconditioning
      flavour=3
    end select

    call m_allocate(H,m,m,m_storage)
    call m_allocate(S,m,m,m_storage)
    call m_allocate(D_min,m,m,m_storage)
    call m_allocate(ED_min,m,m,m_storage)
    call m_allocate(C_min,n,m,m_storage)

    do i=1,2

      select case(i)
      case (1)
        ! First step: build H and S from input data
        do j=1,num_elements
          call m_set_element(H, &
                             elements_indices(1,j), &
                             elements_indices(2,j), &
                             elements_values(1,j), &
                             0.0_dp, &
                             m_operation)
          call m_set_element(S, &
                             elements_indices(1,j), &
                             elements_indices(2,j), &
                             elements_values(2,j), &
                             0.0_dp, &
                             m_operation)
        end do
        new_S=.true.
        init_C=.false.
      case (2)
        ! Second step: rebuild H and introduce a small perturbation
        call m_set(H, &
                   'a', &
                   0.0_dp, &
                   0.0_dp, &
                   m_operation)
        do j=1,num_elements
          call m_set_element(H, &
                             elements_indices(1,j), &
                             elements_indices(2,j), &
                             elements_values(1,j), &
                             0.0_dp, &
                             m_operation)
        end do
        do j=1,m
          call m_set_element(H, &
                             j, &
                             j, &
                             -0.001_dp, &
                             1.0_dp, &
                             m_operation)
        end do
        new_S=.false.
        init_C=.true.
      end select

      ! Solve with OMM
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
               init_C, &    ! init_C
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

      if (mpi_rank==0) print('(a,i1,a,f21.15)'), 'e_min [', i, '] : ', e_min
      call assert_equal_dp(e_min, e_min_check(i))

      call m_get_element(D_min,1,1,el)
      if (mpi_rank==0) print('(a,i1,a,f21.15)'), 'D_11 [', i, ']  : ', el
      call assert_equal_dp(el, D_el_check(i))

    end do

    ! Build the energy-weighted density matrix
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
             .true., &    ! init_C
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
    if (mpi_rank==0) print('(a,f21.15)'), 'ED_11     : ', el
    call assert_equal_dp(el, ED_el_check)

    call m_deallocate(C_min)
    call m_deallocate(ED_min)
    call m_deallocate(D_min)
    call m_deallocate(S)
    call m_deallocate(H)

  end do

  deallocate(elements_values)
  deallocate(elements_indices)

#ifdef HAVE_MPI
  call mpi_finalize(mpi_err)
#endif

  contains

  subroutine assert_equal_dp(value1, value2)
    implicit none

    !**** PARAMS **********************************!

    real(dp), parameter :: tolerance=1.0d-4

    !**** INPUT ***********************************!

    real(dp), intent(in) :: value1
    real(dp), intent(in) :: value2

    !**********************************************!

    if (abs(value1-value2)>tolerance) stop 1

  end subroutine assert_equal_dp

end program example_gamma
