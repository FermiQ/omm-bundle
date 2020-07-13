#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief Main MatrixSwitch module.
!==============================================================================!
module MatrixSwitch
  use MatrixSwitch_ops
  use MatrixSwitch_mm_multiply
  use MatrixSwitch_m_add
  use MatrixSwitch_m_set
  use MatrixSwitch_m_copy
  use MatrixSwitch_m_register
#ifdef HAVE_PSPBLAS
  use pspBLAS
#endif

#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
  use dbcsr_api
#endif

  implicit none

  private

  !**** INTERFACES ********************************!

  !============================================================================!
  !> @brief Matrix-matrix multiplication.
  !!
  !! Performs the operation:
  !! C := alpha*op(A)*op(B) + beta*C, where op(M) = {M, M^T, M^H}
  !!
  !! @param[in]    A     Matrix A. Note that the definition of the matrix
  !!                     (real/complex) needs to be the same as for the other
  !!                     matrices.
  !! @param[in]    opA   Form of op(A):
  !!                     \arg \c n / \c N A
  !!                     \arg \c t / \c T A^T
  !!                     \arg \c c / \c C A^H (equivalent to A^T for a real
  !!                     matrix)
  !! @param[in]    B     Matrix B. Note that the definition of the matrix
  !!                     (real/complex) needs to be the same as for the other
  !!                     matrices.
  !! @param[in]    opB   Form of op(B):
  !!                     \arg \c n / \c N B
  !!                     \arg \c t / \c T B^T
  !!                     \arg \c c / \c C B^H (equivalent to B^T for a real
  !!                     matrix)
  !! @param[inout] C     Matrix C. Note that the definition of the matrix
  !!                     (real/complex) needs to be the same as for the other
  !!                     matrices.
  !! @param[in]    alpha Scalar alpha. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrices (real/complex); otherwise, it only has
  !!                     to match the type of \p beta, and will be
  !!                     automatically converted to match the matrices.
  !! @param[in]    beta  Scalar beta. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrices (real/complex); otherwise, it only has
  !!                     to match the type of \p alpha, and will be
  !!                     automatically converted to match the matrices.
  !! @param[in]    label Implementation of the operation to use. See online
  !!                     documentation for the list of available
  !!                     implementations.
  !============================================================================!
  interface mm_multiply
     module procedure mm_dmultiply
     module procedure mm_zmultiply
  end interface mm_multiply

  !============================================================================!
  !> @brief Matrix addition.
  !!
  !! Performs the operation:
  !! C := alpha*op(A) + beta*C, where op(M) = {M, M^T, M^H}
  !!
  !! @param[in]    A     Matrix A. Note that the definition of the matrix
  !!                     (real/complex) needs to be the same as for the other
  !!                     matrices.
  !! @param[in]    opA   Form of op(A):
  !!                     \arg \c n / \c N A
  !!                     \arg \c t / \c T A^T
  !!                     \arg \c c / \c C A^H (equivalent to A^T for a real
  !!                     matrix)
  !! @param[inout] C     Matrix C. Note that the definition of the matrix
  !!                     (real/complex) needs to be the same as for the other
  !!                     matrices.
  !! @param[in]    alpha Scalar alpha. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrices (real/complex); otherwise, it only has
  !!                     to match the type of \p beta, and will be
  !!                     automatically converted to match the matrices.
  !! @param[in]    beta  Scalar beta. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrices (real/complex); otherwise, it only has
  !!                     to match the type of \p alpha, and will be
  !!                     automatically converted to match the matrices.
  !! @param[in]    label Implementation of the operation to use. See online
  !!                     documentation for the list of available
  !!                     implementations.
  !============================================================================!
  interface m_add
     module procedure m_dadd
     module procedure m_zadd
  end interface m_add

  !============================================================================!
  !> @brief Matrix trace.
  !!
  !! Performs the operation:
  !! alpha := tr(A)
  !!
  !! @param[in]  A     Matrix A.
  !! @param[out] alpha Scalar alpha. If the library is compiler without the
  !!                   \c -DCONV flag, the type has to match the definition of
  !!                   the matrix (real/complex); otherwise, it will be
  !!                   automatically converted to match it.
  !! @param[in]  label Implementation of the operation to use. See online
  !!                   documentation for the list of available implementations.
  !============================================================================!
  interface m_trace
     module procedure m_dtrace
     module procedure m_ztrace
  end interface m_trace

  !============================================================================!
  !> @brief Matrix product trace.
  !!
  !! Performs the operation:
  !! alpha := tr(A^H*B) = tr(B*A^H)
  !!
  !! @param[in]  A     Matrix A. Note that the definition of the matrix
  !!                   (real/complex) needs to be the same as for the other
  !!                   matrix.
  !! @param[in]  B     Matrix B. Note that the definition of the matrix
  !!                   (real/complex) needs to be the same as for the other
  !!                   matrix.
  !! @param[out] alpha Scalar alpha. If the library is compiler without the
  !!                   \c -DCONV flag, the type has to match the definition of
  !!                   the matrix (real/complex); otherwise, it will be
  !!                   automatically converted to match it.
  !! @param[in]  label Implementation of the operation to use. See online
  !!                   documentation for the list of available implementations.
  !============================================================================!
  interface mm_trace
     module procedure mm_dtrace
     module procedure mm_ztrace
  end interface mm_trace

  !============================================================================!
  !> @brief Scale matrix.
  !!
  !! Performs the operation:
  !! C := beta*C
  !!
  !! @param[in]  C     Matrix C.
  !! @param[out] beta  Scalar beta. If the library is compiler without the
  !!                   \c -DCONV flag, the type has to match the definition of
  !!                   the matrix (real/complex); otherwise, it will be
  !!                   automatically converted to match it.
  !! @param[in]  label Implementation of the operation to use. See online
  !!                   documentation for the list of available implementations.
  !============================================================================!
  interface m_scale
     module procedure m_dscale
     module procedure m_zscale
  end interface m_scale

  !============================================================================!
  !> @brief Set matrix.
  !!
  !! Performs the operation:
  !! C_ij := {alpha (i /= j), beta (i == j)}
  !!
  !! @param[inout] C     Matrix C.
  !! @param[in]    seC   Form of the operation:
  !!                     \arg \c l / \c L lower triangle
  !!                     \arg \c u / \c U upper triangle
  !!                     \arg other: complete matrix
  !! @param[in]    alpha Scalar alpha. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrix (real/complex); otherwise, it only has
  !!                     to match the type of \p beta, and will be
  !!                     automatically converted to match the matrix.
  !! @param[in]    beta  Scalar beta. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrix (real/complex); otherwise, it only has
  !!                     to match the type of \p alpha, and will be
  !!                     automatically converted to match the matrix.
  !! @param[in]    label Implementation of the operation to use. See online
  !!                     documentation for the list of available
  !!                     implementations.
  !============================================================================!
  interface m_set
     module procedure m_dset
     module procedure m_zset
  end interface m_set

  !============================================================================!
  !> @brief Set matrix element.
  !!
  !! Performs the operation:
  !! C_ij := alpha + beta*C_ij
  !!
  !! @param[inout] C     Matrix C.
  !! @param[in]    i     Row index of the element.
  !! @param[in]    j     Column index of the element.
  !! @param[in]    alpha Scalar alpha. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrix (real/complex); otherwise, it only has
  !!                     to match the type of \p beta, and will be
  !!                     automatically converted to match the matrix.
  !! @param[in]    beta  Scalar beta. If the library is compiler without the
  !!                     \c -DCONV flag, the type has to match the definition
  !!                     of the matrix (real/complex); otherwise, it only has
  !!                     to match the type of \p alpha, and will be
  !!                     automatically converted to match the matrix.
  !! @param[in]    label Implementation of the operation to use. See online
  !!                     documentation for the list of available
  !!                     implementations.
  !============================================================================!
  interface m_set_element
     module procedure m_dset_element
     module procedure m_zset_element
     module procedure m_dset_block
  end interface m_set_element

  !============================================================================!
  !> @brief Get matrix element.
  !!
  !! Performs the operation:
  !! alpha := C_ij
  !!
  !! @param[in]  C     Matrix C.
  !! @param[in]  i     Row index of the element.
  !! @param[in]  j     Column index of the element.
  !! @param[out] alpha Scalar alpha. If the library is compiler without the
  !!                   \c -DCONV flag, the type has to match the definition of
  !!                   the matrix (real/complex); otherwise, it will be
  !!                   automatically converted to match it.
  !! @param[in]  label Implementation of the operation to use. See online
  !!                   documentation for the list of available implementations.
  !============================================================================!
  interface m_get_element
     module procedure m_dget_element
     module procedure m_zget_element
     module procedure m_dget_block
  end interface m_get_element

  !============================================================================!
  !> @brief Allocate matrix.
  !!
  !! Two procedures are provided for elements and blocked matrices.
  !!
  !============================================================================!
  interface m_allocate
     module procedure m_allocate_elements
     module procedure m_allocate_blocks
     module procedure m_allocate_eq_blocks
  end interface m_allocate

  !************************************************!

  public :: matrix
  public :: m_allocate
  public :: m_deallocate
  public :: m_copy
  public :: m_read
  public :: m_write
  public :: m_convert
  public :: mm_multiply
  public :: m_add
  public :: m_trace
  public :: mm_trace
  public :: m_scale
  public :: m_set
  public :: m_setzero
  public :: m_set_element
  public :: m_get_element
  public :: m_register_sden
#ifdef HAVE_MPI
  public :: m_register_pdbc
  public :: ms_lap_icontxt
#ifdef HAVE_SCALAPACK
  public :: ms_scalapack_setup
#endif
#ifdef HAVE_DBCSR
  public :: ms_dbcsr_setup
  public :: ms_dbcsr_finalize
  public :: m_convert_dbcsrcsr
  public :: m_convert_csrdbcsr
  public :: m_register_pdcsr
  public :: m_set_diag_dbcsr
#endif
#endif
#ifdef HAVE_PSPBLAS
  public :: m_copy_external_pdbcpcoo
  public :: m_copy_external_pdbcpcsc
  public :: m_register_pcoo
  public :: m_register_pcsc
#endif

contains

  !============================================================================!
  !> @brief Allocate matrix.
  !!
  !! Initializes a TYPE(MATRIX) variable by saving some basic information about
  !! the matrix, and allocating the necessary arrays for the requested storage
  !! format. Matrix elements are set to zero.
  !!
  !! @param[inout] m_name The matrix to be allocated.
  !! @param[in]    i      Row dimension size of the matrix.
  !! @param[in]    j      Column dimension size of the matrix.
  !! @param[in]    label  Storage format to use. See online documentation for
  !!                      the list of available formats. Default is \c sdden.
  !============================================================================!
  subroutine m_allocate_elements(m_name,i,j,label)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    character(1) :: c1, c2

    integer :: st

    !**********************************************!

    m_name%dim1=i
    m_name%dim2=j
    if (i==j) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if

    if (present(label)) then
       read(label,'(a1,a1,a3)') c1, c2, m_name%str_type
    else
       c1='s'
       c2='d'
       m_name%str_type='den'
    end if

    if (c1 .eq. 's') then
       m_name%is_serial=.true.
    else if (c1 .eq. 'p') then
       m_name%is_serial=.false.
#ifndef HAVE_MPI
       call die('m_allocate_elements: compile with MPI')
#endif
    else
       call die('m_allocate_elements: invalid label')
    end if
    if (c2 .eq. 'd') then
       m_name%is_real=.true.
    else if (c2 .eq. 'z') then
       m_name%is_real=.false.
    else
       call die('m_allocate_elements: invalid label')
    end if

    ! storage type
    if (m_name%is_serial) then
       if (m_name%str_type .eq. 'den') then
          m_name%is_sparse=.false.
          st=1
       else if (m_name%str_type .eq. 'coo') then
          m_name%is_sparse=.true.
          st=3
       else if (m_name%str_type .eq. 'csc') then
          m_name%is_sparse=.true.
          st=3
       else if (m_name%str_type .eq. 'csr') then
          m_name%is_sparse=.true.
          st=3
       else
          call die('m_allocate_elements: invalid label')
       end if
    else
       if (m_name%str_type .eq. 'dbc') then
          m_name%is_sparse=.false.
          st=2
       else
          call die('m_allocate_elements: invalid label')
       end if
    end if

    select case (st)
    case (1)
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
       end if
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call ms_scalapack_allocate(m_name)
#else
       call die('m_allocate_elements: compile with ScaLAPACK')
#endif
    case (3)
       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
    end select

    m_name%is_initialized=.true.

  end subroutine m_allocate_elements

  !============================================================================!
  !> @brief Allocate blocked matrix.
  !!
  !! Initializes a TYPE(MATRIX) variable by saving some basic information about
  !! the matrix, and allocating the necessary arrays for the requested storage
  !! format. It uses a block-cycling distribution for the blocks.
  !! The matrix is created empty.
  !!
  !! @param[inout] m_name    The matrix to be allocated.
  !! @param[in]    row_sizes Row block dimensions of the matrix.
  !! @param[in]    col_sizes Column block dimensions of the matrix.
  !! @param[in]    label  Storage format to use. See online documentation for
  !!                      the list of available formats. Default is \c pdcsr.
  !============================================================================!

  subroutine m_allocate_blocks(m_name,row_sizes,col_sizes,label,use2D)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label
    logical, intent(in), optional :: use2D

    integer, intent(in), dimension(:), pointer :: row_sizes
    integer, intent(in), dimension(:), pointer :: col_sizes

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    character(1) :: c1, c2

    integer :: st

    !**********************************************!

    m_name%dim1=SUM(row_sizes)
    m_name%dim2=SUM(col_sizes)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if

    if (present(label)) then
       read(label,'(a1,a1,a3)') c1, c2, m_name%str_type
    else
       c1='p'
       c2='d'
       m_name%str_type='csr'
    end if

    if (c1 .eq. 's') then
       call die('m_allocate_blocks: serial algorithm with blocks is not implemented yet')
    else if (c1 .eq. 'p') then
       m_name%is_serial=.false.
#ifndef HAVE_MPI
       call die('m_allocate_blocks: compile with MPI')
#endif
#ifndef HAVE_DBCSR
       call die('m_allocate_blocks: DBCSR is not linked')
#endif
    else
       call die('m_allocate_blocks: invalid label')
    end if
    !
    if (c2 .eq. 'd') then
       m_name%is_real=.true.
    else if (c2 .eq. 'z') then
       m_name%is_real=.false.
    else
       call die('m_allocate: invalid label')
    end if
    !
    ! storage type
    if (.not. m_name%is_serial) then
       if (m_name%str_type .eq. 'csr') then
          m_name%is_sparse=.true.
          st=1
       else
          call die('m_allocate_blocks: invalid label (only csr available)')
       end if
    end if

    select case (st)
    case (1)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       if(present(use2D)) then
         call ms_dbcsr_allocate(m_name, row_sizes, col_sizes, use2D)
         m_name%Use2D = use2D
       else
         call ms_dbcsr_allocate(m_name, row_sizes, col_sizes)
       end if
#else
       call die('m_allocate_blocks: compile with MPI and DBCSR')
#endif
    end select

    m_name%is_initialized=.true.

  end subroutine m_allocate_blocks

  subroutine m_allocate_eq_blocks(m_name,dim1,dim2,blocksize1,blocksize2,label,use2D)

    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label
    logical, intent(in), optional :: use2D

    integer, intent(in) :: dim1, dim2
    integer, intent(in) :: blocksize1, blocksize2

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer:: nblocks1, nblocks2
    integer, dimension(:), pointer :: row_sizes, col_sizes
    integer :: i

    !**********************************************!

    m_name%blk_size1=blocksize1
    m_name%blk_size2=blocksize2
    nblocks1 = ceiling(real(dim1,dp)/blocksize1)
    nblocks2 = ceiling(real(dim2,dp)/blocksize2)
    allocate(row_sizes(1:nblocks1))
    row_sizes(:) = blocksize1
    allocate(col_sizes(1:nblocks2))
    col_sizes(:) = blocksize2

    if(present(use2D)) then
      if(present(label)) then
        call m_allocate_blocks(m_name,row_sizes,col_sizes,label=label,use2D=use2D)
      else
        call m_allocate_blocks(m_name,row_sizes,col_sizes,use2D=use2D)
      end if
    else
       if(present(label)) then
        call m_allocate_blocks(m_name,row_sizes,col_sizes,label=label)
      else
        call m_allocate_blocks(m_name,row_sizes,col_sizes)
      end if
    end if

  end subroutine m_allocate_eq_blocks

  !============================================================================!
  !> @brief Deallocate matrix.
  !!
  !! Deallocates any allocated arrays in a TYPE(MATRIX) variable. For a
  !! registered matrix, the pointers are nullified.
  !!
  !! @param[inout] m_name The matrix to be deallocated.
  !============================================================================!
  subroutine m_deallocate(m_name)
    implicit none

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

    if (associated(m_name%iaux1)) then
       if (m_name%iaux1_is_allocated) then
          deallocate(m_name%iaux1)
          m_name%iaux1_is_allocated=.false.
       else
          nullify(m_name%iaux1)
       end if
    end if
    if (associated(m_name%iaux2)) then
       if (m_name%iaux2_is_allocated) then
          deallocate(m_name%iaux2)
          m_name%iaux2_is_allocated=.false.
       else
          nullify(m_name%iaux2)
       end if
    end if
    if (associated(m_name%iaux3)) then
       if (m_name%iaux3_is_allocated) then
          deallocate(m_name%iaux3)
          m_name%iaux3_is_allocated=.false.
       else
          nullify(m_name%iaux3)
       end if
    end if
    if (associated(m_name%iaux4)) then
       if (m_name%iaux4_is_allocated) then
          deallocate(m_name%iaux4)
          m_name%iaux4_is_allocated=.false.
       else
          nullify(m_name%iaux4)
       end if
    end if
    if (associated(m_name%dval)) then
       if (m_name%dval_is_allocated) then
          deallocate(m_name%dval)
          m_name%dval_is_allocated=.false.
       else
          nullify(m_name%dval)
       end if
    end if
    if (associated(m_name%zval)) then
       if (m_name%zval_is_allocated) then
          deallocate(m_name%zval)
          m_name%zval_is_allocated=.false.
       else
          nullify(m_name%zval)
       end if
    end if

#ifdef HAVE_PSPBLAS
    if (((m_name%str_type .eq. 'coo') .or. &
         (m_name%str_type .eq. 'csc')) .and. &
        (.not. m_name%is_serial)) call psp_deallocate_spm(m_name%spm)
#endif

#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
    if ((m_name%str_type .eq. 'csr') .and. &
         (.not. m_name%is_serial)) then
       call dbcsr_release(m_name%dbcsr_mat)
       call dbcsr_distribution_release(m_name%dbcsr_dist)
       if(associated(m_name%csr_val)) nullify(m_name%csr_val)
    endif
#endif

    m_name%is_initialized=.false.

  end subroutine m_deallocate

  !============================================================================!
  !> @brief Copy matrix.
  !!
  !! Initializes a TYPE(MATRIX) variable by copying the information and element
  !! values from an existing one. The storage format of the new matrix can be
  !! different from that of the old matrix, thus allowing to convert between
  !! formats. This can be useful, e.g., to convert from a dense to a sparse
  !! format or vice-versa, or between different sparse formats. The optional
  !! thresholding variables allow to shave an existing matrix, increasing the
  !! level of sparsity.
  !!
  !! @param[inout] m_name            The matrix to be created.
  !! @param[in]    A                 The matrix to be copied.
  !! @param[in]    label             Storage format to use for \p m_name.
  !! @param[in]    threshold         Tolerance for zeroing elements. Elements
  !!                                 with an absolute value below this
  !!                                 threshold will be omitted for sparse
  !!                                 storage formats, and set to zero for dense
  !!                                 storage formats. For blocks, the threshold
  !!                                 applies to the Frobenius norm of the blocks.
  !! @param[in]    threshold_is_soft Is the thresholding soft?
  !!                                 \arg \c .true. Values above \p threshold
  !!                                 are shifted down to remove the jump
  !!                                 discontinuity.
  !!                                 \arg \c .false. Values are not shifted
  !!                                 (default).
  !============================================================================!
  subroutine m_copy(m_name,A,label,threshold,threshold_is_soft,keep_sparsity)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    logical, intent(in), optional :: keep_sparsity

    real(dp), intent(in), optional :: threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    character(1) :: c1, c2

    integer :: st, i, j, k

    real(dp) :: abs_threshold, soft_threshold

    logical :: my_threshold_is_soft

    !**********************************************!

    m_name%dim1=A%dim1
    m_name%dim2=A%dim2
    m_name%is_square=A%is_square

    if (present(label)) then
       read(label,'(a1,a1,a3)') c1, c2, m_name%str_type
       if (c1 .eq. 's') then
          m_name%is_serial=.true.
       else if (c1 .eq. 'p') then
          m_name%is_serial=.false.
#ifndef HAVE_MPI
          call die('m_copy: compile with MPI')
#endif
       else
          call die('m_copy: invalid label')
       end if
       if (c2 .eq. 'd') then
          m_name%is_real=.true.
       else if (c2 .eq. 'z') then
          m_name%is_real=.false.
       else
          call die('m_copy: invalid label')
       end if
    else
       m_name%is_serial=A%is_serial
       m_name%is_real=A%is_real
       m_name%str_type=A%str_type
    end if

    if (m_name%is_real .neqv. A%is_real) call die('m_copy: invalid label')

    ! storage type
    if ((m_name%str_type .eq. 'den') .and. &
         (m_name%is_serial)) then
       m_name%is_sparse=.false.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=9
       else if ((A%str_type .eq. 'csc') .and. &
                (A%is_serial)) then
          st=10
       else if ((A%str_type .eq. 'csr') .and. &
                (A%is_serial)) then
          st=11
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=1
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'dbc') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.false.
       if ((A%str_type .eq. 'dbc') .and. &
             (.not. A%is_serial)) then
          st=2
       endif
    else if ((m_name%str_type .eq. 'coo') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=3
       else if ((A%str_type .eq. 'csc') .and. &
                (A%is_serial)) then
          st=12
       else if ((A%str_type .eq. 'csr') .and. &
                (A%is_serial)) then
          st=13
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=4
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'coo') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
          st=16
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csc') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=14
       else if ((A%str_type .eq. 'csc') .and. &
           (A%is_serial)) then
          st=5
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=6
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csc') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
          st=17
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csr') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=15
       else if ((A%str_type .eq. 'csr') .and. &
           (A%is_serial)) then
          st=7
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=8
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csr') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'csr') .and. &
            (.not. A%is_serial)) then
          st=18
       else
          call die('m_copy: invalid label')
       endif
    else
       call die('m_copy: invalid label')
    end if

    my_threshold_is_soft = .false.

    if (present(threshold)) then
       abs_threshold=abs(threshold)
       if (present(threshold_is_soft)) then
          if (threshold_is_soft) then
             my_threshold_is_soft = .true.
             soft_threshold=abs_threshold
          else
             soft_threshold=0.0_dp
          end if
       else
          soft_threshold=0.0_dp
       end if
    else
       abs_threshold=0.0_dp
       soft_threshold=0.0_dp
    end if

    select case (st)
    case (1)
       if (present(threshold)) then
          call m_copy_sdensdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdensdenref(m_name,A)
       end if
    case (2)
       if (present(threshold)) then
          call m_copy_pdbcpdbcref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_pdbcpdbcref(m_name,A)
       end if
    case (3)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscooref(m_name,A)
       end if
    case (4)
       if (present(threshold)) then
          call m_copy_sdenscooref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscooref(m_name,A)
       end if
    case (5)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scscscscref(m_name,A)
       end if
    case (6)
       if (present(threshold)) then
          call m_copy_sdenscscref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscscref(m_name,A)
       end if
    case (7)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scsrscsrref(m_name,A)
       end if
    case (8)
       if (present(threshold)) then
          call m_copy_sdenscsrref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscsrref(m_name,A)
       end if
    case (9)
       if (present(threshold)) then
          call m_copy_scoosdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scoosdenref(m_name,A)
       end if
    case (10)
       if (present(threshold)) then
          call m_copy_scscsdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scscsdenref(m_name,A)
       end if
    case (11)
       if (present(threshold)) then
          call m_copy_scsrsdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scsrsdenref(m_name,A)
       end if
    case (12)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scscscooref(m_name,A)
       end if
    case (13)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scsrscooref(m_name,A)
       end if
    case (14)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscscref(m_name,A)
       end if
    case (15)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscsrref(m_name,A)
       end if
    case (16)
       if (.not. present(threshold)) then
          call die('m_copy: threshold must be specified')
       else
          if (my_threshold_is_soft) then
             call die('m_copy: soft thresholding not yet implemented')
          else
#ifdef HAVE_PSPBLAS
             if (A%is_real) then
                call m_copy_external_pdbcpcoo(m_name,A%dval,A%iaux1,threshold)
             else
                call m_copy_external_pdbcpcoo(m_name,A%zval,A%iaux1,threshold)
             end if
#else
             call die('m_copy: compile with pspBLAS')
#endif
          end if
       end if
    case (17)
       if (.not. present(threshold)) then
          call die('m_copy: threshold must be specified')
       else
          if (my_threshold_is_soft) then
             call die('m_copy: soft thresholding not yet implemented')
          else
#ifdef HAVE_PSPBLAS
             if (A%is_real) then
                call m_copy_external_pdbcpcsc(m_name,A%dval,A%iaux1,threshold)
             else
                call m_copy_external_pdbcpcsc(m_name,A%zval,A%iaux1,threshold)
             end if
#else
             call die('m_copy: compile with pspBLAS')
#endif
          end if
       end if
       case (18)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
          if (my_threshold_is_soft) then
             call die('m_copy: soft thresholding not yet implemented')
          else
             ! First check if A is valid
             if (.not. dbcsr_valid_index(A%dbcsr_mat)) then
                call die('m_copy: Invalid DBCSR matrix to copy from')
             endif
             ! Copy DBCSR in DBCSR (deep-copy)
             if(present(keep_sparsity)) then
               call dbcsr_copy(m_name%dbcsr_mat, A%dbcsr_mat, keep_sparsity=keep_sparsity)
             else
               call dbcsr_copy(m_name%dbcsr_mat, A%dbcsr_mat)
             end if
             if (present(threshold)) then
                call dbcsr_filter(m_name%dbcsr_mat, abs_threshold)
             endif
          endif
#else
          call die('m_copy: compile with MPI and DBCSR')
#endif
    end select

    m_name%is_initialized=.true.

  end subroutine m_copy

  !============================================================================!
  !> @brief Convert matrix format.
  !!
  !! This routine facilitates an in-place conversion between storage formats.
  !! Internally it uses \a m_copy to produce a temporary matrix with the new
  !! format, then overwrites the original matrix with this information and
  !! finally deletes the temporary matrix.
  !!
  !! @param[inout] m_name            The matrix to be converted.
  !! @param[in]    label             The new storage format to use.
  !! @param[in]    threshold         Tolerance for zeroing elements. Elements
  !!                                 with an absolute value below this
  !!                                 threshold will be omitted for sparse
  !!                                 storage formats, and set to zero for dense
  !!                                 storage formats.
  !! @param[in]    threshold_is_soft Is the thresholding soft?
  !!                                 \arg \c .true. Values above \p threshold
  !!                                 are shifted down to remove the jump
  !!                                 discontinuity.
  !!                                 \arg \c .false. Values are not shifted
  !!                                 (default).
  !============================================================================!
  subroutine m_convert(m_name,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    type(matrix) :: temp_matrix

    !**********************************************!

    if (present(label)) then
       if (present(threshold)) then
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,label,threshold,threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,label,threshold)
          end if
       else
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,label,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,label)
          end if
       end if
    else
       if (present(threshold)) then
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,threshold=threshold,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,threshold=threshold)
          end if
       else
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name)
          end if
       end if
    end if
    call m_deallocate(m_name)
    call m_copy(m_name,temp_matrix)
    call m_deallocate(temp_matrix)

  end subroutine m_convert

  subroutine m_read(m_name,filepath,file_exist,keep_sparsity)

    implicit none
    !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath
    logical, intent(in), optional :: keep_sparsity

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** OUT ************************************!

    logical, intent(out) :: file_exist

    !**** INTERNAL ********************************!

    integer :: st

    !**********************************************!

    if(m_name%is_serial) then
      call die('m_read: not implemented for serial')
    else
#ifndef HAVE_MPI
      call die('m_read: compile with MPI')
#endif
      if(m_name%is_real) then
        if(m_name%str_type .eq. 'csr') then
          st = 1
        else
          if(m_name%str_type .eq. 'dbc') then
            st = 2
          else
            call die('m_read: not implemented for this matrix type')
          end if
        end if
      else
        call die('m_read: not implemented for complex')
      end if
    end if

    select case (st)
    case (1)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
      if(present(keep_sparsity)) then
        call ms_dbcsr_read(m_name,filepath,file_exist,keep_sparsity=keep_sparsity)
      else
        call ms_dbcsr_read(m_name,filepath,file_exist)
      end if
#else
      call die('m_read: compile with MPI and DBCSR')
#endif
   case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
      call ms_pddbc_read(m_name,filepath,file_exist)
#else
      call die('m_read: compile with MPI and SCALAPACK')
#endif
    end select
  end subroutine m_read


  subroutine m_write(m_name,filepath)
    implicit none

    !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: st

    !**********************************************!

     if(m_name%is_serial) then
        call die('m_write: not implemented for serial')
     else
#ifndef HAVE_MPI
     call die('m_write: compile with MPI')
#endif
     if(m_name%is_real) then
        if(m_name%str_type .eq. 'csr') then
          st = 1
        else
          if(m_name%str_type .eq. 'dbc') then
            st = 2
          else
            call die('m_write: not implemented for this matrix type')
          end if
        end if
      else
        call die('m_write: not implemented for complex')
      end if
    end if

    select case (st)
    case (1)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
      call ms_dbcsr_write(m_name,filepath)
#else
      call die('m_write: compile with MPI and DBCSR')
#endif
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
      call ms_pddbc_write(m_name,filepath)
#else
      call die('m_write: compile with MPI and SCALAPACK')
#endif
    end select
  end subroutine m_write


  !============================================================================!
  !> @brief Matrix-matrix multiplication (real version).
  !============================================================================!
  subroutine mm_dmultiply(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    logical :: trA, trB

    integer :: ot, i

#ifdef CONV
    complex(dp) :: cmplx_alpha, cmplx_beta
#endif

    !**********************************************!

#ifdef CONV
    if ((.not. A%is_real) .and. (.not. B%is_real) .and. (.not. C%is_real)) then
       cmplx_alpha=cmplx(alpha,0.0_dp,dp)
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call mm_zmultiply(A,opA,B,opB,C,cmplx_alpha,cmplx_beta,label)
       return
    end if
#endif
    if (.not. A%is_real) call die('mm_dmultiply: matrix A is complex')
    if (.not. B%is_real) call die('mm_dmultiply: matrix B is complex')
    if (.not. C%is_real) call die('mm_dmultiply: matrix C is complex')
    call process_opM(opA,trA)
    call process_opM(opB,trB)
    if ((.not. trA) .and. (.not. trB)) then
       if ((A%dim1/=C%dim1) .or. &
            (B%dim2/=C%dim2) .or. &
            (A%dim2/=B%dim1)) call die('mm_dmultiply: matrices A, B and C are not compatible')
    else if (trA .and. (.not. trB)) then
       if ((A%dim2/=C%dim1) .or. &
            (B%dim2/=C%dim2) .or. &
            (A%dim1/=B%dim1)) call die('mm_dmultiply: matrices A, B and C are not compatible')
    else if ((.not. trA) .and. trB) then
       if ((A%dim1/=C%dim1) .or. &
            (B%dim1/=C%dim2) .or. &
            (A%dim2/=B%dim2)) call die('mm_dmultiply: matrices A, B and C are not compatible')
    else if (trA .and. trB) then
       if ((A%dim2/=C%dim1) .or. &
            (B%dim1/=C%dim2) .or. &
            (A%dim1/=B%dim2)) call die('mm_dmultiply: matrices A, B and C are not compatible')
    end if

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=2
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'lap') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'psp') then
             ot=4
          else if (label .eq. 't1D') then
             ot=6
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'psp') then
             ot=5
          else if (label .eq. 't1D') then
             ot=7
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=14
       else
          if (label .eq. 't1D') then
             ot=14
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=8
       else
          if (label .eq. 'ref') then
             ot=8
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=9
       else
          if (label .eq. 'ref') then
             ot=9
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=10
       else
          if (label .eq. 'ref') then
             ot=10
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=11
       else
          if (label .eq. 'ref') then
             ot=11
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=12
       else
          if (label .eq. 'ref') then
             ot=12
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csr') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=13
       else
          if (label .eq. 'ref') then
             ot=13
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (.not. B%is_serial)) then
       ot=15
    else
       call die('mm_dmultiply: invalid implementation')
    end if

    select case (ot)
    case (1)
       call mm_multiply_sddenref(A,trA,B,trB,C,alpha,beta)
    case (2)
#ifdef HAVE_LAPACK
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call dgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%dval,A%dim1,B%dval,B%dim1,beta,C%dval,C%dim1)
#else
       call die('mm_dmultiply: compile with LAPACK')
#endif
    case (3)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call pdgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%dval,1,1,A%iaux1,B%dval,1,1,B%iaux1,beta,C%dval,1,1,C%iaux1)
#else
       call die('mm_dmultiply: compile with ScaLAPACK')
#endif
    case (4)
#ifdef HAVE_PSPBLAS
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gespmm(C%dim1,C%dim2,i,A%spm,opA,B%dval,opB,C%dval,alpha,beta)
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (5)
#ifdef HAVE_PSPBLAS
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gemspm(C%dim1,C%dim2,i,A%dval,opA,B%spm,opB,C%dval,alpha,beta)
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (6)
#ifdef HAVE_PSPBLAS
       if (trB) then
          call die('mm_dmultiply: not implemented for transposed B')
       else
          call mm_multiply_pdcscpddbcpddbct1D(A,trA,B,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef HAVE_PSPBLAS
       if (trA) then
          call die('mm_dmultiply: not implemented for transposed A')
       else
          call mm_multiply_pddbcpdcscpddbct1D(A,B,trB,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (8)
       call mm_multiply_sdcscsddensddenref(A,trA,B,trB,C,alpha,beta)
    case (9)
       call mm_multiply_sddensdcscsddenref(A,trA,B,trB,C,alpha,beta)
    case (10)
       call mm_multiply_sddensddensdcscref(A,trA,B,trB,C,alpha,beta)
    case (11)
       call mm_multiply_sdcsrsddensddenref(A,trA,B,trB,C,alpha,beta)
    case (12)
       call mm_multiply_sddensdcsrsddenref(A,trA,B,trB,C,alpha,beta)
    case (13)
       call mm_multiply_sddensddensdcsrref(A,trA,B,trB,C,alpha,beta)
    case (14)
#ifdef HAVE_PSPBLAS
       if (trA .eqv. trB) then
          call die('mm_dmultiply: not implemented for combination of op(A) and op(B)')
       else
          call mm_multiply_pddbcpddbcpdcsct1D(A,trA,B,trB,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (15)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR) 
       call dbcsr_multiply(opA,opB,alpha,A%dbcsr_mat,B%dbcsr_mat,beta,C%dbcsr_mat)
#else
       call die('mm_dmultiply: compile with MPI and DBCSR')
#endif
    end select

  end subroutine mm_dmultiply

  !============================================================================!
  !> @brief Matrix-matrix multiplication (complex version).
  !============================================================================!
  subroutine mm_zmultiply(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: tcA, tcB, ot, i

#ifdef CONV
    real(dp) :: real_alpha, real_beta
#endif

    !**********************************************!

#ifdef CONV
    if (A%is_real .and. B%is_real .and. C%is_real) then
       real_alpha=real(alpha,dp)
       real_beta=real(beta,dp)
       call mm_dmultiply(A,opA,B,opB,C,real_alpha,real_beta,label)
       return
    end if
#endif
    if (A%is_real) call die('mm_zmultiply: matrix A is real')
    if (B%is_real) call die('mm_zmultiply: matrix B is real')
    if (C%is_real) call die('mm_zmultiply: matrix C is real')
    call process_opM(opA,tcA)
    call process_opM(opB,tcB)
    if ((tcA==0) .and. (tcB==0)) then
       if ((A%dim1/=C%dim1) .or. &
            (B%dim2/=C%dim2) .or. &
            (A%dim2/=B%dim1)) call die('mm_zmultiply: matrices A, B and C are not compatible')
    else if ((tcA>0) .and. (tcB==0)) then
       if ((A%dim2/=C%dim1) .or. &
            (B%dim2/=C%dim2) .or. &
            (A%dim1/=B%dim1)) call die('mm_zmultiply: matrices A, B and C are not compatible')
    else if ((tcA==0) .and. (tcB>0)) then
       if ((A%dim1/=C%dim1) .or. &
            (B%dim1/=C%dim2) .or. &
            (A%dim2/=B%dim2)) call die('mm_zmultiply: matrices A, B and C are not compatible')
    else if ((tcA>0) .and. (tcB>0)) then
       if ((A%dim2/=C%dim1) .or. &
            (B%dim1/=C%dim2) .or. &
            (A%dim1/=B%dim2)) call die('mm_zmultiply: matrices A, B and C are not compatible')
    end if

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=2
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'lap') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'psp') then
             ot=4
          else if (label .eq. 't1D') then
             ot=6
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'psp') then
             ot=5
          else if (label .eq. 't1D') then
             ot=7
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=14
       else
          if (label .eq. 't1D') then
             ot=14
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=8
       else
          if (label .eq. 'ref') then
             ot=8
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=9
       else
          if (label .eq. 'ref') then
             ot=9
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=10
       else
          if (label .eq. 'ref') then
             ot=10
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=11
       else
          if (label .eq. 'ref') then
             ot=11
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=12
       else
          if (label .eq. 'ref') then
             ot=12
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csr') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=13
       else
          if (label .eq. 'ref') then
             ot=13
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else
       call die('mm_zmultiply: invalid implementation')
    end if

    select case (ot)
    case (1)
       call mm_multiply_szdenref(A,tcA,B,tcB,C,alpha,beta)
    case (2)
#ifdef HAVE_LAPACK
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call zgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%zval,A%dim1,B%zval,B%dim1,beta,C%zval,C%dim1)
#else
       call die('mm_zmultiply: compile with LAPACK')
#endif
    case (3)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call pzgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%zval,1,1,A%iaux1,B%zval,1,1,B%iaux1,beta,C%zval,1,1,C%iaux1)
#else
       call die('mm_zmultiply: compile with ScaLAPACK')
#endif
    case (4)
#ifdef HAVE_PSPBLAS
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gespmm(C%dim1,C%dim2,i,A%spm,opA,B%zval,opB,C%zval,alpha,beta)
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (5)
#ifdef HAVE_PSPBLAS
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gemspm(C%dim1,C%dim2,i,A%zval,opA,B%spm,opB,C%zval,alpha,beta)
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (6)
#ifdef HAVE_PSPBLAS
       if (tcB>0) then
          call die('mm_zmultiply: not implemented for transposed B')
       else
          call mm_multiply_pzcscpzdbcpzdbct1D(A,tcA,B,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef HAVE_PSPBLAS
       if (tcA>0) then
          call die('mm_zmultiply: not implemented for transposed A')
       else
          call mm_multiply_pzdbcpzcscpzdbct1D(A,B,tcB,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (8)
       call mm_multiply_szcscszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (9)
       call mm_multiply_szdenszcscszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (10)
       call mm_multiply_szdenszdenszcscref(A,tcA,B,tcB,C,alpha,beta)
    case (11)
       call mm_multiply_szcsrszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (12)
       call mm_multiply_szdenszcsrszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (13)
       call mm_multiply_szdenszdenszcsrref(A,tcA,B,tcB,C,alpha,beta)
    case (14)
#ifdef HAVE_PSPBLAS
       if ((tcA==tcB) .or. (tcA+tcB>2)) then
          call die('mm_zmultiply: not implemented for combination of op(A) and op(B)')
       else
          call mm_multiply_pzdbcpzdbcpzcsct1D(A,tcA,B,tcB,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    end select

  end subroutine mm_zmultiply

  !============================================================================!
  !> @brief Matrix addition (real version).
  !============================================================================!
  subroutine m_dadd(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    logical :: trA

    integer :: ot

#ifdef CONV
    complex(dp) :: cmplx_alpha, cmplx_beta
#endif

#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
    type(dbcsr_type) :: dbcsr_trA
#endif

    !**********************************************!

#ifdef CONV
    if ((.not. A%is_real) .and. (.not. C%is_real)) then
       cmplx_alpha=cmplx(alpha,0.0_dp,dp)
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call m_zadd(A,opA,C,cmplx_alpha,cmplx_beta,label)
       return
    end if
#endif
    if (.not. A%is_real) call die('m_dadd: matrix A is complex')
    if (.not. C%is_real) call die('m_dadd: matrix C is complex')
    call process_opM(opA,trA)
    if (.not. (trA)) then
       if ((A%dim1/=C%dim1) .or. &
            (A%dim2/=C%dim2)) call die('m_dadd: matrices A and C are not compatible')
    else
       if ((A%dim1/=C%dim2) .or. &
            (A%dim2/=C%dim1)) call die('m_dadd: matrices A and C are not compatible')
    end if

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'ref') then
             ot=4
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'ref') then
             ot=5
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'csr') .and. &
         (.not. C%is_serial)) then
       ot=6
    else
       call die('m_dadd: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_add_sddenref(A,trA,C,alpha,beta)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pdgeadd(opA,C%dim1,C%dim2,alpha,A%dval,1,1,A%iaux1,beta,C%dval,1,1,C%iaux1)
#else
       call die('m_dadd: compile with ScaLAPACK')
#endif
    case (3)
#ifdef HAVE_PSPBLAS
       if (trA) call die('m_dadd: implementation only valid for opA=''n''')
       if ((A%spm%loc_dim1/=C%iaux2(1)) .or. &
           (A%spm%loc_dim2/=C%iaux2(2))) call die('m_dadd: matrices A and C must have identical parallel distributions')
       call m_add_pdcscpddbcref(A,C,alpha,beta)
#else
       call die('m_dadd: compile with pspBLAS')
#endif
    case (4)
       call m_add_sdcscsddenref(A,trA,C,alpha,beta)
    case (5)
       call m_add_sdcsrsddenref(A,trA,C,alpha,beta)
    case (6)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       if (trA) then
          call dbcsr_transposed(dbcsr_trA, A%dbcsr_mat)
          call dbcsr_add(C%dbcsr_mat, dbcsr_trA, beta, alpha)
          call dbcsr_release(dbcsr_trA)
       else
          call dbcsr_add(C%dbcsr_mat, A%dbcsr_mat, beta, alpha)
       endif
#else
       call die('m_dadd: compile with MPI and DBCSR')
#endif
    end select

  end subroutine m_dadd

  !============================================================================!
  !> @brief Matrix addition (complex version).
  !============================================================================!
  subroutine m_zadd(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: tcA, ot

#ifdef CONV
    real(dp) :: real_alpha, real_beta
#endif

    !**********************************************!

#ifdef CONV
    if (A%is_real .and. C%is_real) then
       real_alpha=real(alpha,dp)
       real_beta=real(beta,dp)
       call m_dadd(A,opA,C,real_alpha,real_beta,label)
       return
    end if
#endif
    if (A%is_real) call die('m_zadd: matrix A is real')
    if (C%is_real) call die('m_zadd: matrix C is real')
    call process_opM(opA,tcA)
    if (tcA==0) then
       if ((A%dim1/=C%dim1) .or. &
            (A%dim2/=C%dim2)) call die('m_zadd: matrices A and C are not compatible')
    else
       if ((A%dim1/=C%dim2) .or. &
            (A%dim2/=C%dim1)) call die('m_zadd: matrices A and C are not compatible')
    end if

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'ref') then
             ot=4
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'ref') then
             ot=5
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else
       call die('m_zadd: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_add_szdenref(A,tcA,C,alpha,beta)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pzgeadd(opA,C%dim1,C%dim2,alpha,A%zval,1,1,A%iaux1,beta,C%zval,1,1,C%iaux1)
#else
       call die('m_zadd: compile with ScaLAPACK')
#endif
    case (3)
#ifdef HAVE_PSPBLAS
       if (tcA>0) call die('m_zadd: implementation only valid for opA=''n''')
       if ((A%spm%loc_dim1/=C%iaux2(1)) .or. &
           (A%spm%loc_dim2/=C%iaux2(2))) call die('m_zadd: matrices A and C must have identical parallel distributions')
       call m_add_pzcscpzdbcref(A,C,alpha,beta)
#else
       call die('m_zadd: compile with pspBLAS')
#endif
    case (4)
       call m_add_szcscszdenref(A,tcA,C,alpha,beta)
    case (5)
       call m_add_szcsrszdenref(A,tcA,C,alpha,beta)
    end select

  end subroutine m_zadd

  !============================================================================!
  !> @brief Matrix trace (real version).
  !============================================================================!
  subroutine m_dtrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    type(matrix), intent(in) :: A

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, i

#ifdef CONV
    complex(dp) :: cmplx_alpha
#endif

    !**** EXTERNAL ********************************!

#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
    real(dp), external :: pdlatra
#endif

    !**********************************************!

#ifdef CONV
    if (.not. A%is_real) then
       call m_ztrace(A,cmplx_alpha,label)
       alpha=real(cmplx_alpha,dp)
       return
    end if
#else
    if (.not. A%is_real) call die('m_dtrace: matrix A is complex')
#endif
    if (.not. A%is_square) call die('m_dtrace: matrix A is not square')

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_dtrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_dtrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (.not. A%is_serial)) then
       ot=3
    else
       call die('m_dtrace: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=0.0_dp
       do i=1,A%dim1
          alpha=alpha+A%dval(i,i)
       end do
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       alpha=pdlatra(A%dim1,A%dval,1,1,A%iaux1)
#else
       call die('m_dtrace: compile with ScaLAPACK')
#endif
    case (3)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       CALL dbcsr_trace(A%dbcsr_mat, alpha)
#else
       call die('m_dtrace: compile with MPI and DBCSR')
#endif
    end select

  end subroutine m_dtrace

  !============================================================================!
  !> @brief Matrix trace (complex version).
  !============================================================================!
  subroutine m_ztrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    type(matrix), intent(in) :: A

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, i

#ifdef CONV
    real(dp) :: real_alpha
#endif

    !**** EXTERNAL ********************************!

#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
    complex(dp), external :: pzlatra
#endif

    !**********************************************!

#ifdef CONV
    if (A%is_real) then
       call m_dtrace(A,real_alpha,label)
       alpha=cmplx(real_alpha,0.0_dp,dp)
       return
    end if
#else
    if (A%is_real) call die('m_ztrace: matrix A is real')
#endif
    if (.not. A%is_square) call die('m_ztrace: matrix A is not square')

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_ztrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_ztrace: invalid implementation')
          end if
       end if
    else
       call die('m_ztrace: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=cmplx_0
       do i=1,A%dim1
          alpha=alpha+A%zval(i,i)
       end do
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       alpha=pzlatra(A%dim1,A%zval,1,1,A%iaux1)
#else
       call die('m_ztrace: compile with ScaLAPACK')
#endif
    end select

  end subroutine m_ztrace

  !============================================================================!
  !> @brief Matrix product trace (real version).
  !============================================================================!
  subroutine mm_dtrace(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, i, j

#ifdef CONV
    complex(dp) :: cmplx_alpha
#endif

    !**** EXTERNAL ********************************!

#ifdef HAVE_MPI
    integer :: info

    real(dp) :: alpha_loc
#endif
#ifdef HAVE_LAPACK
    real(dp), external :: ddot
#endif

    !**********************************************!

#ifdef CONV
    if ((.not. A%is_real) .and. (.not. B%is_real)) then
       call mm_ztrace(A,B,cmplx_alpha,label)
       alpha=real(cmplx_alpha,dp)
       return
    end if
#endif
    if (.not. A%is_real) call die('mm_dtrace: matrix A is complex')
    if (.not. B%is_real) call die('mm_dtrace: matrix B is complex')
    if ((A%dim1/=B%dim1) .or. &
         (A%dim2/=B%dim2)) call die('mm_dtrace: matrices A and B are not compatible')

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=2
          else
             call die('mm_dtrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'lap') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_dtrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (.not. B%is_serial)) then
       ot=4
    else
       call die('mm_dtrace: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=0.0_dp
       do i=1,A%dim1
          do j=1,A%dim2
             alpha=alpha+A%dval(i,j)*B%dval(i,j)
          end do
       end do
    case (2)
#ifdef HAVE_LAPACK
       alpha=ddot(A%dim1*A%dim2,A%dval,1,B%dval,1)
#else
       call die('mm_dtrace: compile with LAPACK')
#endif
    case (3)
#if defined(HAVE_MPI) && defined(HAVE_LAPACK)
       if ((A%iaux2(1)/=B%iaux2(1)) .or. &
           (A%iaux2(2)/=B%iaux2(2))) call die('mm_dtrace: matrices A and B must have identical parallel distributions')
       alpha_loc=ddot(A%iaux2(1)*A%iaux2(2),A%dval,1,B%dval,1)
       call mpi_allreduce(alpha_loc,alpha,1,mpi_double_precision,mpi_sum,ms_mpi_comm,info)
       if (info/=0) call die('mm_dtrace: error in mpi_allreduce')
#else
       call die('mm_dtrace: compile with MPI + LAPACK')
#endif
    case (4)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       call dbcsr_trace(A%dbcsr_mat, B%dbcsr_mat, alpha)
#else
       call die("mm_dtrace: compile with MPI and DBCSR")
#endif
    end select

  end subroutine mm_dtrace

  !============================================================================!
  !> @brief Matrix product trace (complex version).
  !============================================================================!
  subroutine mm_ztrace(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, i, j

#ifdef CONV
    real(dp) :: real_alpha
#endif

    !**** EXTERNAL ********************************!

#ifdef HAVE_MPI
    integer :: info

    complex(dp) :: alpha_loc
#endif

    !**********************************************!

#ifdef CONV
    if (A%is_real .and. B%is_real) then
       call mm_dtrace(A,B,real_alpha,label)
       alpha=cmplx(real_alpha,0.0_dp,dp)
       return
    end if
#endif
    if (A%is_real) call die('mm_ztrace: matrix A is real')
    if (B%is_real) call die('mm_ztrace: matrix B is real')
    if ((A%dim1/=B%dim1) .or. &
         (A%dim2/=B%dim2)) call die('mm_ztrace: matrices A and B are not compatible')

    ! operation table
    if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('mm_ztrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'lap') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_ztrace: invalid implementation')
          end if
       end if
    else
       call die('mm_ztrace: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=cmplx_0
       do i=1,A%dim1
          do j=1,A%dim2
             alpha=alpha+conjg(A%zval(i,j))*B%zval(i,j)
          end do
       end do
    case (3)
#ifdef HAVE_MPI
       if ((A%iaux2(1)/=B%iaux2(1)) .or. &
           (A%iaux2(2)/=B%iaux2(2))) call die('mm_ztrace: matrices A and B must have identical parallel distributions')
       alpha_loc=cmplx_0
       do i=1,A%iaux2(1)
          do j=1,A%iaux2(2)
             alpha_loc=alpha_loc+conjg(A%zval(i,j))*B%zval(i,j)
          end do
       end do
       call mpi_allreduce(alpha_loc,alpha,1,mpi_double_complex,mpi_sum,ms_mpi_comm,info)
       if (info/=0) call die('mm_ztrace: error in mpi_allreduce')
#else
       call die('mm_ztrace: compile with MPI')
#endif
    end select

  end subroutine mm_ztrace

  !============================================================================!
  !> @brief Scale matrix (real version).
  !============================================================================!
  subroutine m_dscale(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    real(dp), intent(in) :: beta

    !**** INTOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: ot

#ifdef CONV
    complex(dp) :: cmplx_beta
#endif

    !**********************************************!

#ifdef CONV
    if (.not. C%is_real) then
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call m_zscale(C,cmplx_beta,label)
       return
    end if
#else
    if (.not. C%is_real) call die('m_dscale: matrix C is complex')
#endif

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
             (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csc') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csr') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csr') .and. &
         (.not. C%is_serial)) then
       ot=2
    else
       call die('m_dscale: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%dval=beta*C%dval
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       call dbcsr_scale(C%dbcsr_mat, beta)
#else
       call die('m_dscale: compile with DBCSR')
#endif
    end select

  end subroutine m_dscale

  !============================================================================!
  !> @brief Scale matrix (complex version).
  !============================================================================!
  subroutine m_zscale(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: beta

    !**** INTOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: ot

#ifdef CONV
    real(dp) :: real_beta
#endif

    !**********************************************!

#ifdef CONV
    if (C%is_real) then
       real_beta=real(beta,dp)
       call m_dscale(C,real_beta,label)
       return
    end if
#else
    if (C%is_real) call die('m_zscale: matrix C is complex')
#endif

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
             (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csc') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csr') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else
       call die('m_zscale: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%zval=beta*C%zval
    end select

  end subroutine m_zscale

  subroutine m_setzero(m_name, label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

    if((.not. m_name%is_serial) .and. (m_name%is_real) .and.&
      (m_name%str_type .eq. 'csr')) then
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
      call ms_setzero_dbcsr(m_name)
#else
      call die('m_write: compile with MPI and DBCSR')
#endif
    else
      call m_set(m_name,'a',0.0_dp,0.0_dp,label)
    endif

  end subroutine m_setzero

  !============================================================================!
  !> @brief Set matrix (real version).
  !============================================================================!
  subroutine m_dset(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: ot, i, j

#ifdef CONV
    complex(dp) :: cmplx_alpha, cmplx_beta
#endif

    !**********************************************!

#ifdef CONV
    if (.not. C%is_real) then
       cmplx_alpha=cmplx(alpha,0.0_dp,dp)
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call m_zset(C,seC,cmplx_alpha,cmplx_beta,label)
       return
    end if
#else
    if (.not. C%is_real) call die('m_dset: matrix C is complex')
#endif

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_dset: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_dset: invalid implementation')
          end if
       end if
    else
       call die('m_dset: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_set_sddenref(C,seC,alpha,beta)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pdlaset(seC,C%dim1,C%dim2,alpha,beta,C%dval,1,1,C%iaux1)
#else
       call die('m_dset: compile with ScaLAPACK')
#endif
    end select

  end subroutine m_dset

  !============================================================================!
  !> @brief Set matrix (complex version).
  !============================================================================!
  subroutine m_zset(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: ot, i, j

#ifdef CONV
    real(dp) :: real_alpha, real_beta
#endif

    !**********************************************!

#ifdef CONV
    if (C%is_real) then
       real_alpha=real(alpha,dp)
       real_beta=real(beta,dp)
       call m_dset(C,seC,real_alpha,real_beta,label)
       return
    end if
#else
    if (C%is_real) call die('m_zset: matrix C is real')
#endif

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_zset: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_zset: invalid implementation')
          end if
       end if
    else
       call die('m_zset: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_set_szdenref(C,seC,alpha,beta)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pzlaset(seC,C%dim1,C%dim2,alpha,beta,C%zval,1,1,C%iaux1)
#else
       call die('m_zset: compile with ScaLAPACK')
#endif
    end select

  end subroutine m_zset

  !============================================================================!
  !> @brief Set matrix element (real version).
  !============================================================================!
  subroutine m_dset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    logical :: el_present

    integer :: ot, k, buffer
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

    real(dp) :: el
    real(dp), allocatable :: dval_temp(:,:)

#ifdef CONV
    complex(dp) :: cmplx_alpha, cmplx_beta
#endif

    !**********************************************!

#ifdef CONV
    if (.not. C%is_real) then
       cmplx_alpha=cmplx(alpha,0.0_dp,dp)
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call m_zset_element(C,i,j,cmplx_alpha,cmplx_beta,label)
       return
    end if
#else
    if (.not. C%is_real) call die('m_dset_element: matrix C is complex')
#endif
    if ((i<1) .or. &
         (i>C%dim1) .or. &
         (j<1) .or. &
         (j>C%dim2)) call die('m_dset_element: element out of range')

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
            call die('m_dset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
            call die('m_dset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
            call die('m_dset_element: invalid implementation')
          end if
       end if
    else
       call die('m_dset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%dval(i,j)=alpha+beta*C%dval(i,j)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       if (beta==0.0_dp) then
          call pdelset(C%dval,i,j,C%iaux1,alpha)
       else
          call pdelget('p',' ',el,C%dval,i,j,C%iaux1)
          call pdelset(C%dval,i,j,C%iaux1,alpha+beta*el)
       end if
#else
       call die('m_dset_element: compile with ScaLAPACK')
#endif
    case (3)
       if (C%iaux2(1)==0) then
          C%iaux2(1)=1
          buffer=min(C%dim1,C%dim2)
          allocate(C%iaux3(buffer))
          C%iaux3_is_allocated=.true.
          allocate(C%iaux4(buffer))
          C%iaux4_is_allocated=.true.
          allocate(C%dval(buffer,1))
          C%dval_is_allocated=.true.
          C%iaux3(1)=i
          C%iaux4(1)=j
          C%dval(1,1)=alpha
       else
          el_present=.false.
          do k=1,C%iaux2(1)
             if ((C%iaux3(k)==i) .and. &
                 (C%iaux4(k)==j)) then
                C%dval(k,1)=alpha+beta*C%dval(k,1)
                el_present=.true.
                exit
             end if
          end do
          if (.not. el_present) then
             if (C%iaux2(1)==size(C%iaux3)) then
                allocate(iaux3_temp(C%iaux2(1)))
                allocate(iaux4_temp(C%iaux2(1)))
                allocate(dval_temp(C%iaux2(1),1))
                iaux3_temp=C%iaux3
                iaux4_temp=C%iaux4
                dval_temp=C%dval
                deallocate(C%dval)
                C%dval_is_allocated=.false.
                deallocate(C%iaux4)
                C%iaux4_is_allocated=.false.
                deallocate(C%iaux3)
                C%iaux3_is_allocated=.false.
                buffer=C%iaux2(1)+min(C%dim1,C%dim2)
                allocate(C%iaux3(buffer))
                C%iaux3_is_allocated=.true.
                allocate(C%iaux4(buffer))
                C%iaux4_is_allocated=.true.
                allocate(C%dval(buffer,1))
                C%dval_is_allocated=.true.
                C%iaux3(1:C%iaux2(1))=iaux3_temp(1:C%iaux2(1))
                C%iaux4(1:C%iaux2(1))=iaux4_temp(1:C%iaux2(1))
                C%dval(1:C%iaux2(1),1)=dval_temp(1:C%iaux2(1),1)
                deallocate(dval_temp)
                deallocate(iaux4_temp)
                deallocate(iaux3_temp)
             end if
             C%iaux2(1)=C%iaux2(1)+1
             C%iaux3(C%iaux2(1))=i
             C%iaux4(C%iaux2(1))=j
             C%dval(C%iaux2(1),1)=alpha
          end if
       end if
    end select

  end subroutine m_dset_element

  !============================================================================!
  !> @brief Set matrix element (complex version).
  !============================================================================!
  subroutine m_zset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    logical :: el_present

    integer :: ot, k, buffer
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

#ifdef CONV
    real(dp) :: real_alpha, real_beta
#endif

    complex(dp) :: el
    complex(dp), allocatable :: zval_temp(:,:)

    !**********************************************!

#ifdef CONV
    if (C%is_real) then
       real_alpha=real(alpha,dp)
       real_beta=real(beta,dp)
       call m_dset_element(C,i,j,real_alpha,real_beta,label)
       return
    end if
#else
    if (C%is_real) call die('m_zset_element: matrix C is real')
#endif
    if ((i<1) .or. &
         (i>C%dim1) .or. &
         (j<1) .or. &
         (j>C%dim2)) call die('m_zset_element: element out of range')

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_zset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_zset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_zset_element: invalid implementation')
          end if
       end if
    else
       call die('m_zset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%zval(i,j)=alpha+beta*C%zval(i,j)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       if (beta==cmplx_0) then
          call pzelset(C%zval,i,j,C%iaux1,alpha)
       else
          call pzelget('p',' ',el,C%zval,i,j,C%iaux1)
          call pzelset(C%zval,i,j,C%iaux1,alpha+beta*el)
       end if
#else
       call die('m_zset_element: compile with ScaLAPACK')
#endif
    case (3)
       if (C%iaux2(1)==0) then
          C%iaux2(1)=1
          buffer=min(C%dim1,C%dim2)
          allocate(C%iaux3(buffer))
          C%iaux3_is_allocated=.true.
          allocate(C%iaux4(buffer))
          C%iaux4_is_allocated=.true.
          allocate(C%zval(buffer,1))
          C%zval_is_allocated=.true.
          C%iaux3(1)=i
          C%iaux4(1)=j
          C%zval(1,1)=alpha
       else
          el_present=.false.
          do k=1,C%iaux2(1)
             if ((C%iaux3(k)==i) .and. &
                 (C%iaux4(k)==j)) then
                C%zval(k,1)=alpha+beta*C%zval(k,1)
                el_present=.true.
                exit
             end if
          end do
          if (.not. el_present) then
             if (C%iaux2(1)==size(C%iaux3)) then
                allocate(iaux3_temp(C%iaux2(1)))
                allocate(iaux4_temp(C%iaux2(1)))
                allocate(zval_temp(C%iaux2(1),1))
                iaux3_temp=C%iaux3
                iaux4_temp=C%iaux4
                zval_temp=C%zval
                deallocate(C%zval)
                C%zval_is_allocated=.false.
                deallocate(C%iaux4)
                C%iaux4_is_allocated=.false.
                deallocate(C%iaux3)
                C%iaux3_is_allocated=.false.
                buffer=C%iaux2(1)+min(C%dim1,C%dim2)
                allocate(C%iaux3(buffer))
                C%iaux3_is_allocated=.true.
                allocate(C%iaux4(buffer))
                C%iaux4_is_allocated=.true.
                allocate(C%zval(buffer,1))
                C%zval_is_allocated=.true.
                C%iaux3(1:C%iaux2(1))=iaux3_temp(1:C%iaux2(1))
                C%iaux4(1:C%iaux2(1))=iaux4_temp(1:C%iaux2(1))
                C%zval(1:C%iaux2(1),1)=zval_temp(1:C%iaux2(1),1)
                deallocate(zval_temp)
                deallocate(iaux4_temp)
                deallocate(iaux3_temp)
             end if
             C%iaux2(1)=C%iaux2(1)+1
             C%iaux3(C%iaux2(1))=i
             C%iaux4(C%iaux2(1))=j
             C%zval(C%iaux2(1),1)=alpha
          end if
       end if
    end select

  end subroutine m_zset_element

  !============================================================================!
  !> @brief Set matrix block (real version).
  !============================================================================!
  subroutine m_dset_block(C,i,j,block_data,beta)
    implicit none

    type(matrix), intent(inout) :: C
    integer, intent(in) :: i, j
    real(dp), dimension(:, :), intent(in) :: block_data
    real(dp), intent(in) :: beta

    integer :: st, myproc, proc_holds_blk
    real(dp), dimension(:, :), pointer :: myblock
    logical :: found

    if ((C%str_type .eq. 'csr') .and. &
        (.not. C%is_serial) .and. &
        (C%is_sparse)) then
       st = 1
    else
       call die('m_dset_block: invalid operation')
    endif

    select case (st)
    case (1)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       myproc = ms_mpi_rank
       call dbcsr_get_stored_coordinates(C%dbcsr_mat, i, j, proc_holds_blk)
       if (myproc .eq. proc_holds_blk) then
          if (beta==0.0_dp) then
             call dbcsr_put_block(C%dbcsr_mat, i, j, block_data)
          else
             call dbcsr_get_block_p(C%dbcsr_mat, i, j, myblock, found)
             call dbcsr_put_block(C%dbcsr_mat, i, j, block_data)
             if (found) then
                call dbcsr_put_block(C%dbcsr_mat, i, j, myblock, .true., beta)
             endif
          endif
          call dbcsr_finalize(C%dbcsr_mat)
       endif
#endif
    end select
  end subroutine m_dset_block

  !============================================================================!
  !> @brief Get matrix element (real version).
  !============================================================================!
  subroutine m_dget_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    type(matrix), intent(in) :: C

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, k

#ifdef CONV
    complex(dp) :: cmplx_alpha
#endif

    !**********************************************!

#ifdef CONV
    if (.not. C%is_real) then
       call m_zget_element(C,i,j,cmplx_alpha,label)
       alpha=real(cmplx_alpha,dp)
       return
    end if
#else
    if (.not. C%is_real) call die('m_dget_element: matrix C is complex')
#endif
    if ((i<1) .or. &
         (i>C%dim1) .or. &
         (j<1) .or. &
         (j>C%dim2)) call die('m_dget_element: element out of range')

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_dget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_dget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_dget_element: invalid implementation')
          end if
       end if
    else
       call die('m_dget_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=C%dval(i,j)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pdelget('a',' ',alpha,C%dval,i,j,C%iaux1)
#else
       call die('m_dget_element: compile with ScaLAPACK')
#endif
    case (3)
       alpha=0.0_dp
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
              (C%iaux4(k)==j)) then
             alpha=C%dval(k,1)
             exit
          end if
       end do
    end select

  end subroutine m_dget_element

  !============================================================================!
  !> @brief Get matrix element (complex version).
  !============================================================================!
  subroutine m_zget_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    type(matrix), intent(in) :: C

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**** INTERNAL ********************************!

    integer :: ot, k

#ifdef CONV
    real(dp) :: real_alpha
#endif

    !**********************************************!

#ifdef CONV
    if (C%is_real) then
       call m_dget_element(C,i,j,real_alpha,label)
       alpha=cmplx(real_alpha,0.0_dp,dp)
       return
    end if
#else
    if (C%is_real) call die('m_zget_element: matrix C is real')
#endif
    if ((i<1) .or. &
         (i>C%dim1) .or. &
         (j<1) .or. &
         (j>C%dim2)) call die('m_zget_element: element out of range')

    ! operation table
    if ((C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'ref') then
             ot=1
          else if (label .eq. 'lap') then
             ot=1
          else
             call die('m_zget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=2
       else
          if (label .eq. 'lap') then
             ot=2
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_zget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_zget_element: invalid implementation')
          end if
       end if
    else
       call die('m_zget_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=C%zval(i,j)
    case (2)
#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
       call pzelget('a',' ',alpha,C%zval,i,j,C%iaux1)
#else
       call die('m_zget_element: compile with ScaLAPACK')
#endif
    case (3)
       alpha=cmplx_0
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
              (C%iaux4(k)==j)) then
             alpha=C%zval(k,1)
             exit
          end if
       end do
    end select

  end subroutine m_zget_element

  !============================================================================!
  !> @brief Get matrix block (real version).
  !============================================================================!
  subroutine m_dget_block(C,i,j,block_data,found)
    implicit none

    type(matrix), intent(inout) :: C
    integer, intent(in) :: i, j
    real(dp), dimension(:, :), pointer, intent(out) :: block_data
    logical, intent(out) :: found

    integer :: st

    if ((C%str_type .eq. 'csr') .and. &
        (.not. C%is_serial) .and. &
        (C%is_sparse)) then
       st = 1
    else
       call die('m_dget_block: invalid operation')
    endif

    select case (st)
    case (1)
#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
       call dbcsr_get_block_p(C%dbcsr_mat, i, j, block_data, found)
#endif
    end select
  end subroutine m_dget_block

#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
  !============================================================================!
  !> @brief ScaLAPACK MPI setup.
  !!
  !! Sets up everything needed to use \c p?dbc matrices with ScaLAPACK. Has to
  !! be called once at the start of the code.
  !!
  !! @param[in] mpi_comm MPI communicator to use.
  !! @param[in] nprow    Row dimension of the process grid (has to be a divisor
  !!                     of the size of the group defined by \p mpi_comm).
  !! @param[in] order    Ordering of the process grid:
  !!                     \arg \c c / \c C column-major ordering
  !!                     \arg \c r / \c R row-major ordering
  !! @param[in] bs_def   Default block size to use when allocating \c p?dbc
  !!                     matrices.
  !! @param[in] bs_list  List of exceptions to \p bs_def to use for specific
  !!                     matrix dimension sizes. Has to be formatted as
  !!                     (\c dim_1, \c bs_1, \c dim_2, \c bs_2, etc.), where
  !!                     \c dim_x is the matrix dimension size, and \c bs_x is
  !!                     the corresponding block size to use for it.
  !! @param[in] icontxt  BLACS context handle, if already initialized.
  !============================================================================!
  subroutine ms_scalapack_setup(mpi_comm,nprow,order,bs_def,bs_list,icontxt)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    character(1), intent(in) :: order

    integer, intent(in) :: mpi_comm
    integer, intent(in) :: nprow
    integer, intent(in) :: bs_def
    integer, intent(in), optional :: bs_list(:)
    integer, intent(in), optional :: icontxt

    !**** INTERNAL ********************************!

    integer :: i, mpi_err

    !**********************************************!

    ms_mpi_comm=mpi_comm
    call mpi_comm_size(ms_mpi_comm,ms_mpi_size,mpi_err)
    call mpi_comm_rank(ms_mpi_comm,ms_mpi_rank,mpi_err)
    ms_lap_nprow=nprow
    ms_lap_npcol=ms_mpi_size/nprow
    if (ms_lap_nprow*ms_lap_npcol/=ms_mpi_size) call die('ms_scalapack_setup: invalid nprow')
    ms_lap_order=order
    ms_lap_bs_def=bs_def
    if (present(bs_list)) then
       if ((size(bs_list)/2)*2/=size(bs_list)) call die('ms_scalapack_setup: invalid bs_list')
       ms_lap_bs_num=size(bs_list)/2
       allocate(ms_lap_bs_list(2,ms_lap_bs_num))
       do i=1,ms_lap_bs_num
          ms_lap_bs_list(1,i)=bs_list((i-1)*2+1)
          ms_lap_bs_list(2,i)=bs_list((i-1)*2+2)
       end do
    end if

    if (present(icontxt)) then
       ms_lap_icontxt=icontxt
    else
       ms_lap_icontxt=ms_mpi_comm
       call blacs_gridinit(ms_lap_icontxt,ms_lap_order,ms_lap_nprow,ms_lap_npcol)
    end if

#ifdef HAVE_PSPBLAS
    ! initialize grid information in pspBLAS
    call psp_gridinit_2D(ms_mpi_comm,ms_mpi_size,ms_lap_nprow,ms_lap_order,ms_lap_bs_def,ms_lap_bs_def,ms_lap_icontxt)
#endif

  end subroutine ms_scalapack_setup
#endif

#if defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
  !============================================================================!
  !> @brief Allocate matrix (dense block cyclic, parallel distribution).
  !!
  !! Initializes a TYPE(MATRIX) variable for the \c p?dbc storage format. Note
  !! that all the common metadata for the matrix must already have been set,
  !! and \a ms_scalapack_setup must have been called.
  !!
  !! @param[inout] A The matrix to be allocated.
  !============================================================================!
  subroutine ms_scalapack_allocate(A)
    implicit none

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, bs1, bs2, info

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    allocate(A%iaux1(9))
    A%iaux1_is_allocated=.true.
    allocate(A%iaux2(2))
    A%iaux2_is_allocated=.true.
    call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
    bs1=ms_lap_bs_def
    bs2=ms_lap_bs_def
    do i=1,ms_lap_bs_num
       if (ms_lap_bs_list(1,i)==A%dim1) then
          bs1=ms_lap_bs_list(2,i)
          exit
       end if
    end do
    do i=1,ms_lap_bs_num
       if (ms_lap_bs_list(1,i)==A%dim2) then
          bs2=ms_lap_bs_list(2,i)
          exit
       end if
    end do
    A%iaux2(1)=numroc(A%dim1,bs1,k,0,ms_lap_nprow)
    A%iaux2(2)=numroc(A%dim2,bs2,l,0,ms_lap_npcol)
    call descinit(A%iaux1,A%dim1,A%dim2,bs1,bs2,0,0,ms_lap_icontxt,max(1,A%iaux2(1)),info)
    if (info/=0) call die('ms_scalapack_allocate: error in descinit')
    if (A%is_real) then
       allocate(A%dval(A%iaux2(1),A%iaux2(2)))
       A%dval_is_allocated=.true.
       A%dval=0.0_dp
    else
       allocate(A%zval(A%iaux2(1),A%iaux2(2)))
       A%zval_is_allocated=.true.
       A%zval=cmplx_0
    end if

  end subroutine ms_scalapack_allocate

  subroutine ms_pddbc_write(m_name, filepath)

    implicit none
    !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: fout
    integer :: blk_r, blk_c, nblks(2)
    integer :: i, j, k, i_node
    integer :: stat(mpi_status_size), mpi_err
    integer :: i0, j0, size, size_i, size_j
    real(dp), allocatable :: block(:)

    !**********************************************!

    if(ms_mpi_rank==0) then
      open(newunit=fout,file=filepath,form="unformatted",status="replace",action="write")
      write(fout) m_name%iaux1
    end if

    nblks(1)=(m_name%iaux1(3)-1)/m_name%iaux1(5)
    nblks(2)=(m_name%iaux1(4)-1)/m_name%iaux1(6)
    size=m_name%iaux1(5)*m_name%iaux1(6)
    allocate(block(1:size))

    do blk_r=0, nblks(1)
      do blk_c=0, nblks(2)
        i_node=mod(blk_r,ms_lap_nprow)+mod(blk_c,ms_lap_npcol)*ms_lap_nprow
        if(ms_lap_order .eq. 'r') i_node=mod(blk_r,ms_lap_nprow)*ms_lap_npcol+mod(blk_c,ms_lap_npcol)
        if((i_node==ms_mpi_rank) .or. (ms_mpi_rank==0)) then
          block(:)=0.0_dp
          i0=blk_r/ms_lap_nprow*m_name%iaux1(5)
          j0=blk_c/ms_lap_npcol*m_name%iaux1(6)
          size_i=m_name%iaux1(5)
          size_j=m_name%iaux1(6)
          if((blk_r+1)*size_i>m_name%iaux1(3)) size_i=m_name%iaux1(3)-blk_r*size_i
          if((blk_c+1)*size_j>m_name%iaux1(4)) size_j=m_name%iaux1(4)-blk_c*size_j
          !size=size_i*size_j
          k=0
          if(i_node==ms_mpi_rank) then
            do i=1, size_i
              do j=1, size_j
                k=k+1
                block(k)=m_name%dval(i0+i,j0+j)
              end do
            end do
            if(i_node==0) then
              write(fout) block(1:size)
            else
              call mpi_send(block(1:size),size,mpi_double_precision,0,1,ms_mpi_comm,mpi_err)
              if (mpi_err/=0) call die('ms_pddbc_write: error in mpi_send')
            end if
          else if(ms_mpi_rank==0) then
            call mpi_recv(block(1:size),size,mpi_double_precision,i_node,1,ms_mpi_comm,stat,mpi_err)
            if (mpi_err/=0) call die('ms_pddbc_write: error in mpi_recv')
            write(fout) block(1:size)
          end if
        end if
      end do
    end do

    deallocate(block)
    if(ms_mpi_rank==0) close(fout)
  end subroutine ms_pddbc_write

  subroutine ms_pddbc_read(m_name, filepath, file_exist)

    implicit none
    !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** OUT ***********************************!

    logical, intent(out) :: file_exist


    !**** INTERNAL ********************************!

    integer :: fout
    integer :: blk_r, blk_c, nblks(2), dim_loc(2)
    integer :: i, j, k, l, i_node, desc(9), desc_old(9)
    integer :: stat(mpi_status_size), mpi_err
    integer :: i0, j0, size, size_i, size_j
    real(dp), allocatable :: block(:), c_read(:,:)
    integer :: info

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    inquire(file=filepath,exist=file_exist)
    if(file_exist) then

      if(ms_mpi_rank==0) then
        open(newunit=fout,file=filepath,form="unformatted",status="old",action="read")
        read(fout) desc
      end if

      call mpi_bcast(desc,9,mpi_integer,0,ms_mpi_comm,mpi_err)

      if((m_name%iaux1(3) .ne. desc(3)) .or. (m_name%iaux1(4) .ne. desc(4))) then
        call die('ms_pddbc_read: Wrong matrix size')
      end if

      call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
      dim_loc(1)=numroc(desc(3),desc(5),k,0,ms_lap_nprow)
      dim_loc(2)=numroc(desc(4),desc(6),l,0,ms_lap_npcol)

      allocate(c_read(1:dim_loc(1),1:dim_loc(2)))
      c_read(:,:) = 0.0_dp

      nblks(1)=(desc(3)-1)/desc(5)
      nblks(2)=(desc(4)-1)/desc(6)
      size=desc(5)*desc(6)
      allocate(block(1:size))

      do blk_r=0, nblks(1)
        do blk_c=0, nblks(2)
          i_node=mod(blk_r,ms_lap_nprow)+mod(blk_c,ms_lap_npcol)*ms_lap_nprow
          if(ms_lap_order .eq. 'r') i_node=mod(blk_r,ms_lap_nprow)*ms_lap_npcol+mod(blk_c,ms_lap_npcol)
          if((i_node==ms_mpi_rank) .or. (ms_mpi_rank==0)) then
            block(:)=0.0_dp
            i0=blk_r/ms_lap_nprow*desc(5)
            j0=blk_c/ms_lap_npcol*desc(6)
            size_i=desc(5)
            size_j=desc(6)
            if((blk_r+1)*size_i>desc(3)) size_i=desc(3)-blk_r*size_i
            if((blk_c+1)*size_j>desc(4)) size_j=desc(4)-blk_c*size_j
            size=size_i*size_j
            k=0
            if(i_node==ms_mpi_rank) then
              if(i_node==0) then
                read(fout) block(1:size)
              else
                call mpi_recv(block(1:size),size,mpi_double_precision,0,1,ms_mpi_comm,stat,mpi_err)
                if (mpi_err/=0) call die('ms_pddbc_read: error in mpi_recv')
              end if
              do i=1, size_i
                do j=1, size_j
                  k=k+1
                  c_read(i0+i,j0+j)=block(k)
                end do
              end do
            else if(ms_mpi_rank==0) then
              read(fout) block(1:size)
              call mpi_send(block(1:size),size,mpi_double_precision,i_node,1,ms_mpi_comm,mpi_err)
              if (mpi_err/=0) call die('ms_pddbc_read: error in mpi_send')
            end if
          end if
       end do
      end do

      if((m_name%iaux1(5)/=desc(5)) .or. (m_name%iaux1(6)/=desc(6))) then
        call descinit(desc_old,desc(3),desc(4),desc(5),desc(6),0,0,ms_lap_icontxt,max(1,dim_loc(1)),info)
        if (info/=0) call die('ms_pddbc_read: error in descinit')

        call pdgemr2d(desc(3),desc(4),c_read,1,1,desc_old,m_name%dval,1,1,m_name%iaux1,ms_lap_icontxt)
      else
        m_name%dval(:,:)=c_read(:,:)
      end if

      deallocate(block)
      if(allocated(c_read)) deallocate(c_read)
      if(ms_mpi_rank==0) close(fout)
    end if
  end subroutine ms_pddbc_read
#endif

#if defined(HAVE_MPI) && defined(HAVE_DBCSR)
  !============================================================================!
  !> @brief DBCSR setup.
  !!
  !! Initialize the DBCSR library. It declares a MPI 2D grid.
  !!
  !! @param[in] MPI communicator.
  !============================================================================!
  subroutine ms_dbcsr_setup(mpi_comm)
    integer, intent(in) :: mpi_comm

    integer                 :: mpi_error
    integer, dimension(2)   :: dims
    logical, dimension(2)   :: period  = .true.
    logical                 :: reorder = .false.

    if (ms_dbcsr_init) then
       call die("ms_dbcsr_setup: DBCSR already initialized")
    endif
    ms_dbcsr_init = .true.

    ms_mpi_comm=mpi_comm
    call mpi_comm_size(ms_mpi_comm,ms_mpi_size,mpi_error)
    call mpi_comm_rank(ms_mpi_comm,ms_mpi_rank,mpi_error)

    ! Set the 2D grid
    dims(:) = 0
    call mpi_dims_create(ms_mpi_size, 2, dims, mpi_error)
    call mpi_cart_create(ms_mpi_comm, 2, dims, period, reorder, ms_dbcsr_group, mpi_error)

    ! Set the 1D grid
    dims(1) = ms_mpi_size
    dims(2) = 1
    call mpi_cart_create(ms_mpi_comm, 2, dims, period, reorder, ms_dbcsr_brd_group, mpi_error)

    ! initialize libdbcsr
    call dbcsr_init_lib()
  end subroutine ms_dbcsr_setup

  !============================================================================!
  !> @brief DBCSR finalization.
  !!
  !! Finalize DBCSR library. It removes the 2D grid communicator.
  !!
  !============================================================================!
  subroutine ms_dbcsr_finalize()

    integer                 :: mpi_error

    if (.not. ms_dbcsr_init) then
       call die("ms_dbcsr_finalize: DBCSR not initialized")
    endif

    ! finalize libdbcsr
    call dbcsr_finalize_lib(ms_dbcsr_group)

    call mpi_comm_free(ms_dbcsr_group, mpi_error)
    ms_dbcsr_group = mpi_comm_null

    call mpi_comm_free(ms_dbcsr_brd_group, mpi_error)
    ms_dbcsr_brd_group = mpi_comm_null

    ms_dbcsr_init = .false.

  end subroutine ms_dbcsr_finalize

  !============================================================================!
  !> @brief Allocate an empty DBCSR matrix.
  !!
  !! Allocate an empty DBCSR matrix. It creates a 2D block-cycling
  !! distribution over processors.
  !!
  !! @param[inout] Matrix.
  !! @param[in] Sizes of the row-blocks
  !! @param[in] Sizes of the column-blocks
  !============================================================================!

  subroutine ms_dbcsr_allocate(A, row_sizes, col_sizes, use2D)
    implicit none

    !**** INOUT ***********************************!
    type(matrix), intent(inout) :: A

    !**** IN **************************************!
    integer, dimension(:), intent(in), pointer :: row_sizes, col_sizes
    logical, intent(in), optional :: use2D

    integer                           :: mpi_error, my_group
    integer, dimension(2)             :: dims, coords
    logical, dimension(2)             :: period
    integer, dimension(:), pointer    :: cd_data, rd_data
    logical                           :: my_use2D

    if (.not. ms_dbcsr_init) then
       call die("ms_dbcsr_allocate: DBCSR not inizialized")
    endif
    
    my_use2D = .true.
    if(present(use2D)) my_use2D = use2D

    if(my_use2D) then
      my_group = ms_dbcsr_group
    else
      my_group = ms_dbcsr_brd_group
    end if
    call mpi_cart_get(my_group, 2, dims, period, coords, mpi_error)

    ! Create distribution arrays (rows and columns)
    ! Block-cycling distribution
    call dbcsr_make_dist(rd_data, SIZE(row_sizes), dims(1))
    call dbcsr_make_dist(cd_data, SIZE(col_sizes), dims(2))

    call dbcsr_distribution_new (A%dbcsr_dist,&
                                 group=my_group,&
                                 row_dist=rd_data,&
                                 col_dist=cd_data,&
                                 reuse_arrays=.true.)

    if (A%is_real) then
       call dbcsr_create(A%dbcsr_mat, "Matrix", A%dbcsr_dist, &
            dbcsr_type_no_symmetry, row_sizes, col_sizes, &
            reuse_arrays=.false., data_type=dbcsr_type_real_8)
    else
       call die("ms_dbcsr_allocate: DBCSR Complex type not implemented!")
    endif

    contains

      subroutine dbcsr_make_dist(dist_array, dist_size, nbins)
        integer, dimension(:), intent(out), pointer        :: dist_array
        integer, intent(in)                                :: dist_size, nbins
    
        integer                                            :: i

        allocate (dist_array(dist_size))
        do i = 1, dist_size
           dist_array(i) = MODULO(i-1, nbins)
        end do
      end subroutine dbcsr_make_dist

  end subroutine ms_dbcsr_allocate

  subroutine m_convert_csrdbcsr(m_name,threshold,bl_size)

    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), optional :: threshold
    integer, intent(in), optional :: bl_size

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    real(dp), dimension(:,:), allocatable :: block_data
    integer :: io, j, jo, ind, nrow, nrows_loc, n_bl, i
    integer :: bsize1, bsize2, dim1, dim2, iblk, nblk
    logical :: use2D = .false., found
    character(5) :: m_storage
    type(matrix) :: brd_m
    integer :: nblk_zero, nblk_tot

    !**********************************************!

    bsize1 = m_name%csr_blk
    bsize2 = bsize1
    if(present(bl_size)) bsize2 = bl_size

    dim1 = m_name%dim1
    dim2 = m_name%dim2

    m_storage ='pdcsr'
    call m_allocate(brd_m,dim1,dim2,bsize1,bsize2,label=m_storage,use2D=use2D)

    allocate(block_data(1:bsize1,1:bsize2))
    block_data(:,:) = 0.0_dp
    nblk = ceiling(real(brd_m%dim2,dp)/bsize2)

    nrows_loc = m_name%csr_nrows
    found = .false.
    nblk_zero = 0
    nblk_tot= 0
    do iblk = 1, nblk
      i = 0
      do io = 1, nrows_loc
        i = i + 1
        do j = 1, m_name%iaux3(io)
           ind = m_name%iaux1(io) + j
           jo = m_name%iaux2(ind)
           jo = jo - (iblk-1) * bsize2
           if((jo .gt. 0) .and. (jo .le. bsize2)) then
             block_data(i,jo) = m_name%csr_val(ind)
             found = .true.
           end if
        end do
        if((mod(io, bsize1) .eq. 0) .or. (io == nrows_loc)) then

          i = 0
          if(found) then
            n_bl = (io - 1)/bsize1
            nrow = ms_mpi_size * n_bl + ms_mpi_rank + 1
            call m_set_element(brd_m, nrow, iblk, block_data,0.0_dp)
            block_data(:,:) = 0.0_dp
            found = .false.
          else
            nblk_zero = nblk_zero + 1
          end if
          nblk_tot = nblk_tot + 1
        end if
      end do
    end do

    deallocate(block_data)
    call dbcsr_complete_redistribute(brd_m%dbcsr_mat, m_name%dbcsr_mat)
    if (present(threshold)) then
      call dbcsr_filter(m_name%dbcsr_mat,abs(threshold))
    endif
    call m_deallocate(brd_m)

  end subroutine m_convert_csrdbcsr

  subroutine m_convert_dbcsrcsr(m_name,threshold,bl_size)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), optional :: threshold
    integer, intent(in), optional :: bl_size

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer   :: io, j, jo, ind, nrow, nrows_loc, n_bl, i
    real(dp),pointer :: myblock(:,:)
    logical :: found
    integer :: bsize1, bsize2, dim1, dim2, iblk, nblk
    logical :: use2D = .false.
    character(5) :: m_storage
    type(matrix) :: brd_m

    !**********************************************!

    bsize1 = m_name%csr_blk
    bsize2 = bsize1
    if(present(bl_size)) bsize2 = bl_size

    dim1 = m_name%dim1
    dim2 = m_name%dim2

    m_storage ='pdcsr'
    call m_allocate(brd_m,dim1,dim2,bsize1,bsize2,label=m_storage,use2D=use2D)

    if (present(threshold)) then
      call dbcsr_filter(m_name%dbcsr_mat,abs(threshold))
    endif
    call dbcsr_complete_redistribute(m_name%dbcsr_mat,brd_m%dbcsr_mat)

    nrows_loc = m_name%csr_nrows

    nblk = ceiling(real(brd_m%dim2,dp)/bsize2)
    found = .false.
    do iblk = 1, nblk
      do io = 1, nrows_loc
        if(mod(io, bsize1) .eq. 1) then
          found = .false.
          n_bl = (io - 1)/bsize1
          nrow = ms_mpi_size * n_bl + ms_mpi_rank + 1
          call m_get_element(brd_m,nrow,iblk,myblock,found)
          i = 0
        end if
        i = i + 1
        do j = 1, m_name%iaux3(io)
          ind = m_name%iaux1(io) + j
          jo = m_name%iaux2(ind)
          jo = jo - bsize2 * (iblk - 1)
          if((jo .gt. 0) .and. (jo .le. bsize2)) then
            if(found) then
              m_name%csr_val(ind) = myblock(i,jo)
            else
              m_name%csr_val(ind) = 0.0_dp
            end if
          end if
        end do
      end do
    end do
    call m_deallocate(brd_m)

  end subroutine m_convert_dbcsrcsr

  subroutine m_set_diag_dbcsr(m_name, alpha)

    !**** INPUT ***********************************!

    real(dp), intent(in) :: alpha

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

    call ms_setzero_dbcsr(m_name)
    call dbcsr_add_on_diag(m_name%dbcsr_mat,alpha)

  end subroutine m_set_diag_dbcsr

  subroutine ms_setzero_dbcsr(m_name)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

    call dbcsr_set(m_name%dbcsr_mat,0.0_dp)
    call dbcsr_filter(m_name%dbcsr_mat,1.0d-15)

  end subroutine ms_setzero_dbcsr

  subroutine ms_dbcsr_write(m_name,filepath)

    !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: fout
    character(len=120) :: filepath_bl

    !**********************************************!

    filepath_bl=trim(filepath)//'_DATA'
    if(ms_mpi_rank==0) then
      open(newunit=fout,file=filepath_bl,form="unformatted",status="replace",action="write")
      write(fout) m_name%dim1, m_name%dim2, m_name%blk_size1, m_name%blk_size2
      close(fout)
    end if
    call dbcsr_binary_write(m_name%dbcsr_mat,filepath)

  end subroutine ms_dbcsr_write

  subroutine ms_dbcsr_read(m_name,filepath,file_exist,keep_sparsity)

   !**** INPUT ***********************************!

    character(len=*), intent(in) :: filepath
    logical, intent(in), optional :: keep_sparsity

   !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

   !**** OUT *************************************!

    logical, intent(out) :: file_exist

   !**** INTERNAL ********************************!

    logical :: file_exist_bl, my_keep_sparsity
    type(matrix) :: C_restart
    integer :: C_data(4), fout, mpi_err
    character(len=120) :: filepath_bl

   !**********************************************!

    filepath_bl=trim(filepath)//'_DATA'
    inquire(file=filepath_bl,exist=file_exist_bl)

    if(file_exist_bl) then
      if(ms_mpi_rank==0) then
        open(newunit=fout,file=filepath_bl,form="unformatted",status="old",action="read")
        read(fout) C_data
        print'(a,i5,a,i5,a,i5,a,i5)',' C_data ',C_data(1),' ',C_data(2),' ',C_data(3),' ',C_data(4)
        close(fout)
      end if

      call mpi_bcast(C_data,4,mpi_integer,0,ms_mpi_comm,mpi_err)
    else
      C_data(1)=m_name%dim1
      C_data(2)=m_name%dim2
      C_data(3)=m_name%blk_size1
      C_data(4)=m_name%blk_size2
    end if

    my_keep_sparsity=.false.
    if(present(keep_sparsity)) my_keep_sparsity=keep_sparsity
    if((m_name%dim1 .ne. C_data(1)) .or. (m_name%dim2 .ne. C_data(2)) &
        .or. (m_name%blk_size1 .ne. C_data(3)) .or. (m_name%blk_size2 .ne. C_data(4)) &
        .or. my_keep_sparsity) then
      inquire(file=filepath,exist=file_exist)
      if(file_exist) then
        if (.not. C_restart%is_initialized) call m_allocate(C_restart,C_data(1),C_data(2),&
          C_data(3),C_data(4),label='pdcsr',use2D=m_name%Use2D)
        call dbcsr_binary_read(filepath,C_restart%dbcsr_dist,ms_dbcsr_group,C_restart%dbcsr_mat)
        if(my_keep_sparsity) then
          call dbcsr_complete_redistribute(C_restart%dbcsr_mat,m_name%dbcsr_mat,&
            keep_sparsity=my_keep_sparsity)
        else
          call dbcsr_complete_redistribute(C_restart%dbcsr_mat,m_name%dbcsr_mat)
        end if
        if(C_restart%is_initialized) call m_deallocate(C_restart)
      end if
    else
      inquire(file=filepath,exist=file_exist)
      if(file_exist) then
        call dbcsr_binary_read(filepath,m_name%dbcsr_dist,ms_dbcsr_group,m_name%dbcsr_mat)
      end if
    end if

  end subroutine ms_dbcsr_read
#endif

end module MatrixSwitch
