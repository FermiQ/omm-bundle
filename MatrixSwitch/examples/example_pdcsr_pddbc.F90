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


!!$! This example demonstrates how to calculate:                                                      !
!!$!   alpha = tr(A^T*B)                                                                              !
!!$! for two NxM matrices A and B, in five different ways. Each way should return the same result.    !
!!$! They are:                                                                                        !
!!$!   1. performing the matrix product trace res1 := tr(A^T*B) directly as a single operation        !
!!$!   2. performing the multiplication D := A^T*B, and then the trace res2 := tr(D)                  !
!!$!   3. performing E := B*A^T, and then res3 := tr(E)                                               !
!!$!   4. performing the transpose C := A^T, then D := C*B, and then res4 := tr(D)                    !
!!$!   5. performing C := A^T, then E := B*C, and then res5 := tr(E)                                  !
!!$! Finally, as an extra check, the first element of E is printed out.                               !
!!$!                                                                                                  !
!!$! Note the difference in the code in how matrix B is handled compared with the other matrices.     !
!!$! While the others are allocated directly by MatrixSwitch (i.e., all the data is contained within  !
!!$! the type(matrix) variable), B is a wrapper for MyMatrix (i.e., MyMatrix has been registered as B !
!!$! for use with MatrixSwitch, and B contains a pointer to MyMatrix).                                !
!!$!                                                                                                  !
!!$! Sample output:                                                                                   !
!!$!--------------------------------------------------------------------------------------------------!
!!$! res1 :    31.194861937321843                                                                     !
!!$! res2 :    31.194861937321846                                                                     !
!!$! res3 :    31.194861937321846                                                                     !
!!$! res4 :    31.194861937321846                                                                     !
!!$! res5 :    31.194861937321846                                                                     !
!!$! E_11 :     2.779421970931733                                                                     !
!!$!--------------------------------------------------------------------------------------------------!
!!$!==================================================================================================!
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
  type(matrix) :: A, B, C, D, E, F
  real(dp) :: val
  real(dp), dimension(4) :: res_sparse, res_dense
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
!!$
!!$  ! DBCSR
!!$  call m_add(A,'t',C,2.0_dp,3.0_dp)
!!$  call m_trace(C,res1)
!!$
!!$  call mm_trace(A,C,res2)
!!$
!!$  call mm_multiply(A,'n',B,'t',C,1.0_dp,1.0_dp)
!!$  call m_trace(C,res3)
!!$
!!$  call m_scale(C,3.0_dp)
!!$  call m_trace(C,res4)
!!$

  !   
  ! dbc
!  call m_add(D,'t',F,2.0_dp,3.0_dp)
!  call m_trace(F,res_dense(1))
  
!  call mm_trace(D,F,res_dense(2))
  
!!$  call mm_multiply(D,'n',E,'t',F,1.0_dp,1.0_dp)
!!$  call m_trace(F,res7)
!!$  
!!$  call m_scale(F,3.0_dp)
!!$  call m_trace(F,res8)
!!$
  if (mpi_rank==0) then
!!$     print('(a,f21.15)'), 'res1 : ', res1
!!$     print('(a,f21.15)'), 'res2 : ', res2
!!$     print('(a,f21.15)'), 'res3 : ', res3
!!$     print('(a,f21.15)'), 'res4 : ', res4
  endif

  if (mpi_rank==0) then
     print('(a)'), 'SCALAPACK'
     print('(a,f21.15)'), 'res1 : ', res_dense(1)
     print('(a,f21.15)'), 'res2 : ', res_dense(2)
!!$     print('(a,f21.15)'), 'res7 : ', res7
!!$     print('(a,f21.15)'), 'res8 : ', res8
  endif

!!$  call assert_equal_dp(res1, res5)
!!$  call assert_equal_dp(res2, res6)
!!$  call assert_equal_dp(res3, res7)
!!$  call assert_equal_dp(res4, res8)
!!$

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

  call m_deallocate(F)
  call m_deallocate(E)
  call m_deallocate(D)
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

#else
  print('(a,f21.15)'), 'To run this example, compile with MPI, SCALAPACK and DBCSR.'
#endif

end program example_pdcsr_pddbc
