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
program example_kpoints
  use MatrixSwitch_wrapper

  implicit none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  real(dp), parameter :: e_min_check(4)=(/-24.43825_dp,-24.45849_dp, &
                                          -24.43941_dp,-24.45965_dp/)
  real(dp), parameter :: D_el_check(4)=(/0.37534_dp,0.37608_dp, &
                                         0.37573_dp,0.37648_dp/)
  real(dp), parameter :: ED_el_check(2)=(/-0.43215_dp,-0.43221_dp/)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation
  character(2), allocatable :: keys(:)

  logical :: new_S, init_C, dealloc

  integer :: num_matrices
  integer :: mpi_err, mpi_size, mpi_rank
  integer :: m, n, num_elements, num_kpoints, i, j, k, f, flavour
  integer, allocatable :: elements_indices(:,:)

  real(dp) :: e_min
  real(dp), allocatable :: elements_values(:,:)

  complex(dp) :: el

  !**********************************************!

#ifdef HAVE_MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  call ms_scalapack_setup(mpi_comm_world,1,'c',3)

  m_storage='pzdbc'
  m_operation='lap'
#else
  mpi_size=1
  mpi_rank=0

  m_storage='szden'
  m_operation='lap'
#endif

  ! Read input data
  if (mpi_rank==0) then
    open(10,file=SRCDIR//'/example_gamma.dat')
    read(10,'(i2,1x,i2,1x,i4)') m, n, num_elements
    allocate(elements_indices(4,num_elements))
    allocate(elements_values(8,num_elements))
    elements_values=0.0_dp
    do i=1,num_elements
      read(10,'(i2,1x,i2,4(1x,es10.3e2))') elements_indices(1:2,i), elements_values(1,i), elements_values(5,i)
    end do
    close(10)
    open(10,file=SRCDIR//'/example_kpoints.dat')
    read(10,'()')
    do i=1,num_elements
      read(10,'(i2,1x,i2,4(1x,es10.3e2))') elements_indices(3:4,i), elements_values(3:4,i), elements_values(7:8,i)
    end do
    close(10)
  end if
#ifdef HAVE_MPI
  call mpi_bcast(m,1,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(n,1,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(num_elements,1,mpi_int,0,mpi_comm_world,mpi_err)
  if (mpi_rank/=0) then
    allocate(elements_indices(4,num_elements))
    allocate(elements_values(8,num_elements))
  end if
  call mpi_bcast(elements_indices,4*num_elements,mpi_int,0,mpi_comm_world,mpi_err)
  call mpi_bcast(elements_values,8*num_elements,mpi_double_precision,0,mpi_comm_world,mpi_err)
#endif

  num_kpoints=2
  num_matrices=6*num_kpoints
  allocate(keys(num_matrices))
  do k=1,num_kpoints
    keys((k-1)*6+1:(k-1)*6+6)=(/key('H',k),key('S',k),key('D',k),key('E',k),key('C',k),key('T',k)/)
  end do
  call ms_wrapper_open(num_matrices,keys)

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

    do k=1,num_kpoints
      call m_allocate(key('H',k),m,m,m_storage)
      call m_allocate(key('S',k),m,m,m_storage)
      call m_allocate(key('D',k),m,m,m_storage)
      call m_allocate(key('E',k),m,m,m_storage)
      call m_allocate(key('C',k),n,m,m_storage)
    end do

    do i=1,2

      do k=1,num_kpoints

        select case(i)
        case (1)
          ! First step: build H and S from input data
          do j=1,num_elements
            call m_set_element(key('H',k), &
                               elements_indices((k-1)*2+1,j), &
                               elements_indices((k-1)*2+2,j), &
                               cmplx(elements_values((k-1)*2+1,j), &
                                     elements_values((k-1)*2+2,j), &
                                     dp), &
                               cmplx_0, &
                               m_operation)
            call m_set_element(key('S',k), &
                               elements_indices((k-1)*2+1,j), &
                               elements_indices((k-1)*2+2,j), &
                               cmplx(elements_values((k-1)*2+5,j), &
                                     elements_values((k-1)*2+6,j), &
                                     dp), &
                               cmplx_0, &
                               m_operation)
          end do
          new_S=.true.
          init_C=.false.
        case (2)
          ! Second step: rebuild H and introduce a small perturbation
          call m_set(key('H',k), &
                     'a', &
                     0.0_dp, &
                     0.0_dp, &
                     m_operation)
          do j=1,num_elements
            call m_set_element(key('H',k), &
                               elements_indices((k-1)*2+1,j), &
                               elements_indices((k-1)*2+2,j), &
                               cmplx(elements_values((k-1)*2+1,j), &
                                     elements_values((k-1)*2+2,j), &
                                     dp), &
                               cmplx_0, &
                               m_operation)
          end do
          do j=1,m
            call m_set_element(key('H',k), &
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
        call omm_wrapper(m, &         ! m
                         n, &         ! n
                         key('H',k), &      ! H 
                         key('S',k), &      ! S
                         new_S, &     ! new_S
                         e_min, &     ! e_min
                         key('D',k), &  ! D_min
                         .false., &   ! calc_ED
                         0.0_dp, &    ! eta
                         key('C',k), &  ! C_min
                         init_C, &    ! init_C
                         key('T',k), &      ! T
                         0.0_dp, &    ! scale_T
                         flavour, &   ! flavour
                         2, &         ! np
                         k, &         ! ip
                         -1.0_dp, &   ! cg_tol
                         .true., &    ! long_out
                         .false., &   ! dealloc
                         m_storage, & ! m_storage
                         m_operation) ! m_operation

        if (mpi_rank==0) print('(2(a,i1),a,f21.15)'), 'e_min [', i, ',', k, '] : ', e_min
        call assert_equal_dp(e_min, e_min_check((k-1)*2+i))

        call m_get_element(key('D',k),1,1,el)
        if (mpi_rank==0) print('(2(a,i1),2(a,f21.15))'), 'D_11 [', i, ',', k, '] : ', real(el,dp), ' , ', aimag(el)
        call assert_equal_dp(real(el,dp), D_el_check((k-1)*2+i))
        call assert_equal_dp(aimag(el), 0.0_dp)

      end do

    end do

    do k=1,num_kpoints

      select case(k)
      case (1)
        dealloc=.false.
      case (2)
        dealloc=.true.
      end select

      ! Build the energy-weighted density matrix
      call omm_wrapper(m, &         ! m
                       n, &         ! n
                       key('H',k), &      ! H 
                       key('S',k), &      ! S
                       .false., &   ! new_S
                       e_min, &     ! e_min
                       key('E',k), & ! D_min
                       .true., &    ! calc_ED
                       0.0_dp, &    ! eta
                       key('C',k), &  ! C_min
                       .true., &    ! init_C
                       key('T',k), &      ! T
                       0.0_dp, &    ! scale_T
                       0, &         ! flavour
                       2, &         ! np
                       k, &         ! ip
                       -1.0_dp, &   ! cg_tol
                       .true., &    ! long_out
                       dealloc, &   ! dealloc
                       m_storage, & ! m_storage
                       m_operation) ! m_operation

      call m_get_element(key('E',k),1,1,el)
      if (mpi_rank==0) print('(a,i1,2(a,f21.15))'), 'ED_11 [', k, '] : ', real(el,dp), ' , ', aimag(el)
      call assert_equal_dp(real(el,dp), ED_el_check(k))
      call assert_equal_dp(aimag(el), 0.0_dp)

      call m_deallocate(key('C',k))
      call m_deallocate(key('E',k))
      call m_deallocate(key('D',k))
      call m_deallocate(key('S',k))
      call m_deallocate(key('H',k))

    end do

  end do

  call ms_wrapper_close()
  deallocate(keys)
  deallocate(elements_values)
  deallocate(elements_indices)

#ifdef HAVE_MPI
  call mpi_finalize(mpi_err)
#endif

  contains

  character(2) function key(name, k)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: name

    integer, intent(in) :: k

    !**** INTERNAL ********************************!

    character(1) :: k_char

    !**********************************************!

    write(k_char,'(i1)') k
    key=name//k_char

  end function key

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

end program example_kpoints
