#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==================================================================================================!
! example : parallel program with real matrices in CSR and dense block cyclic formats              !
!                                                                                                  !
! This example demonstrates how to use DBCSR for the CSR format and it compares the results to     !
! SCALAPACK. Note that in DBCSR case we use 10x10 blocks, where each block has 2x2 elements. The   !
! same block size is used in SCALAPACK. The occupancy of the matrices is 20%.                      !
!                                                                                                  !
! This example evaluates:                                                                          !
!   res1: tr(A^T+B)                                                                                !
!   res2: tr(A*B)                                                                                  !
!   res3: tr(A*B^T)                                                                                !
!   res4: tr(C*alpha)                                                                              !
! Finally, as an extra check, the values on the result matrix are printed out.                     !
!                                                                                                  !
! Sample results on 4 MPI ranks:                                                                   !
!--------------------------------------------------------------------------------------------------!
! res1 :    15.030509283314089                                                                     !
! res2 :    17.674348547066977                                                                     !
! res3 :    12.760482513544938                                                                     !
! res4 :    40.833544043343800                                                                     !
!--------------------------------------------------------------------------------------------------!
!==================================================================================================!
program example_pdcsr_pddbc
  use MatrixSwitch

  implicit none
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK) && defined(HAVE_DBCSR)
  include 'mpif.h'

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  integer, parameter :: block_size = 2
  integer, parameter :: nblocks = 10
  real(dp), parameter :: occupancy = 0.2_dp

  !**** VARIABLES *******************************!

  character(5) :: m_storage

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: full_row_size, full_col_size, ii, jj, ll, kk
  integer, dimension(2) :: dims
  integer, pointer, dimension(:) :: row_blk_sizes, col_blk_sizes
  type(matrix) :: A, B, C, CC, D, E, F, FF
  real(dp) :: val
  real(dp), dimension(5) :: res_sparse, res_dense
  real(dp), dimension(:, :), pointer :: myblock
  logical :: found_block

  !**********************************************!

  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  ! Set the 2D grid dimensions for DBCSR and SCALAPACK
  dims(:) = 0
  call mpi_dims_create(mpi_size, 2, dims, mpi_err)

  ! Setup DBCSR
  call ms_dbcsr_setup(mpi_comm_world)

  ! Set SCALAPACK
  ! Se use row-order because it is used in DBCSR
  call ms_scalapack_setup(mpi_comm_world, dims(1), 'r', block_size)

  ! Set the block sizes
  allocate(row_blk_sizes(nblocks), col_blk_sizes(nblocks))
  row_blk_sizes(:) = block_size
  col_blk_sizes(:) = block_size

  ! Allocate the DBCSR matrices
  m_storage='pdcsr'
  call m_allocate(A,row_blk_sizes,col_blk_sizes,m_storage)
  call m_allocate(B,row_blk_sizes,col_blk_sizes,m_storage)
  call m_allocate(C,row_blk_sizes,col_blk_sizes,m_storage)
  
  ! Allocate the Scalapack matrices
  full_row_size = SUM(row_blk_sizes)
  full_col_size = SUM(col_blk_sizes)
  m_storage='pddbc'
  call m_allocate(D,full_row_size,full_col_size,m_storage)
  call m_allocate(E,full_row_size,full_col_size,m_storage)
  call m_allocate(F,full_row_size,full_col_size,m_storage)

  ! Fill the matrices
  do ii=1,SIZE(row_blk_sizes)
     do jj=1,SIZE(col_blk_sizes)
        call fill_blocks(A, D, ii, jj)
        call fill_blocks(B, E, ii, jj)
        call fill_blocks(C, F, ii, jj)
     end do
  end do
  
  ! DBCSR
  call run_operations(A, B, C, res_sparse)
  ! dbc
  call run_operations(D, E, F, res_dense)

  ! Test Copies
  ! CSR -> CSR
  call m_copy(CC, C)
  call m_trace(CC, res_sparse(5))
  ! dbc -> dbc
  call m_copy(FF, F)
  call m_trace(FF, res_dense(5))

  ! Print results
  if (mpi_rank==0) then
     print('(a,f21.15)'), 'res1 : ', res_sparse(1)
     print('(a,f21.15)'), 'res2 : ', res_sparse(2)
     print('(a,f21.15)'), 'res3 : ', res_sparse(3)
     print('(a,f21.15)'), 'res4 : ', res_sparse(4)
  endif

  do ii=1, SIZE(res_sparse)
     call assert_equal_dp(res_sparse(ii), res_dense(ii))
  enddo

  ! Print values on a given MPI rank
  if (mpi_rank==0) then
     do ii=1,SIZE(row_blk_sizes)
        do jj=1,SIZE(col_blk_sizes)
           call m_get_element(C,ii,jj,myblock,found_block)
           if (found_block) then
              print *, "block found at (",ii,",",jj,") on rank 0. Values= ", myblock
              do ll=1, block_size
                 do kk=1, block_size
                    call m_get_element(F, (ii-1)*block_size+kk, (jj-1)*block_size+ll, val)
                    call assert_equal_dp(myblock(kk, ll), val)
                 enddo
              enddo
           endif
        enddo
     enddo
  endif

  call m_deallocate(FF)
  call m_deallocate(F)
  call m_deallocate(E)
  call m_deallocate(D)
  call m_deallocate(CC)
  call m_deallocate(C)
  call m_deallocate(B)
  call m_deallocate(A)

  deallocate(row_blk_sizes, col_blk_sizes)

  ! Finalize DBCSR
  call ms_dbcsr_finalize()

  call mpi_finalize(mpi_err)

  contains

    subroutine fill_blocks(mat_sparse, mat_dense, irow, icol)
      type(matrix), intent(inout) :: mat_sparse, mat_dense
      integer, intent(in)         :: irow, icol

      real(dp), dimension(block_size, block_size) :: block_data
      integer :: kk, ll

      if (RAND() .lt. occupancy) then
         call RANDOM_NUMBER(block_data)
         ! Fill sparse matrix with the entire block
         call m_set_element(mat_sparse, irow, icol, block_data, 0.0_dp)
      else
         ! Fill the dense matrix with zeros
         block_data = 0
      endif
      ! Fill dense matrix with each element of the block
      do ll=1, block_size
         do kk=1, block_size
            call m_set_element(mat_dense, (irow-1)*block_size+kk, (icol-1)*block_size+ll, block_data(kk,ll), 0.0_dp)
         enddo
      enddo
    end subroutine fill_blocks

    subroutine assert_equal_dp(value1, value2)
      implicit none

      !**** PARAMS **********************************!

      real(dp), parameter :: tolerance=1.0d-10

      !**** INPUT ***********************************!

      real(dp), intent(in) :: value1
      real(dp), intent(in) :: value2

      !**********************************************!

      if (abs(value1-value2)>tolerance) stop 1
    end subroutine assert_equal_dp

    subroutine run_operations(A, B, C, res)
      implicit none

      type(matrix), intent(in) :: A, B
      type(matrix), intent(inout) :: C
      real(dp), dimension(:), intent(out) :: res

      !**********************************************!

      call m_add(A,'t',C,2.0_dp,3.0_dp)
      call m_trace(C,res(1))
    
      call mm_trace(A,C,res(2))
      
      call mm_multiply(A,'n',B,'t',C, 2.0_dp, 0.5_dp)
      call m_trace(C,res(3))

      call m_scale(C,3.2_dp)
      call m_trace(C,res(4))
    end subroutine run_operations

#else
  print('(a,f21.15)'), 'To run this example, compile with MPI, SCALAPACK and DBCSR.'
#endif

end program example_pdcsr_pddbc
