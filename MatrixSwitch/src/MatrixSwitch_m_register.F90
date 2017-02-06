#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief \a m_register functions for different storage types.
!==============================================================================!
module MatrixSwitch_m_register
  use MatrixSwitch_ops

  implicit none

  !**** INTERFACES ********************************!

  !============================================================================!
  !> @brief Register matrix (simple dense, serial distribution).
  !!
  !! Registers pre-existing matrix data into a TYPE(MATRIX) variable with
  !! \c s?den format.
  !!
  !! @param[inout] m_name The matrix to be allocated.
  !! @param[in]    A      The values of the matrix elements, stored as a
  !!                      two-dimensional array.
  !============================================================================!
  interface m_register_sden
     module procedure m_register_sdden
     module procedure m_register_szden
  end interface m_register_sden

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Register matrix (dense block cyclic, parallel distribution).
  !!
  !! Registers pre-existing matrix data into a TYPE(MATRIX) variable with
  !! \c p?dbc format.
  !!
  !! @param[inout] m_name The matrix to be allocated.
  !! @param[in]    A      The values of the local matrix elements, stored as a
  !!                      two-dimensional array.
  !! @param[in]    desc   BLACS array descriptor.
  !============================================================================!
  interface m_register_pdbc
     module procedure m_register_pddbc
     module procedure m_register_pzdbc
  end interface m_register_pdbc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list/compressed sparse column
  !!        from dense block cyclic, parallel distribution).
  !!
  !! Registers pre-existing matrix data into a TYPE(MATRIX) variable with
  !! \c p?coo or \c p?csc format.
  !!
  !! @param[inout] m_name      The matrix to be allocated.
  !! @param[in]    A           The values of the local matrix elements, stored
  !!                           as a two-dimensional array.
  !! @param[in]    desc        BLACS array descriptor.
  !! @param[in]    spm_storage Storage format to use:
  !!                           \arg \c coo Sparse coordinate list.
  !!                           \arg \c csc Compressed sparse column.
  !! @param[in]    thre        Tolerance for zeroing elements. Elements with an
  !!                           absolute value below this threshold will be
  !!                           omitted.
  !============================================================================!
  interface m_register_psp_thre
     module procedure m_register_pdsp_thre
     module procedure m_register_pzsp_thre
  end interface m_register_psp_thre
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list, parallel distribution).
  !!
  !! Registers pre-existing matrix data into a TYPE(MATRIX) variable with
  !! \c p?coo format.
  !!
  !! @param[inout] m_name      The matrix to be allocated.
  !! @param[in]    idx1        The local row indices, stored as a
  !!                           one-dimensional array.
  !! @param[in]    idx2        The local column indices, stored as a
  !!                           one-dimensional array.
  !! @param[in]    val         The values of the local matrix elements, stored
  !!                           as a one-dimensional array.
  !! @param[in]    desc        BLACS array descriptor.
  !============================================================================!
  interface m_register_pcoo
     module procedure m_register_pdcoo
     module procedure m_register_pzcoo
  end interface m_register_pcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (compressed sparse column, parallel distribution).
  !!
  !! Registers pre-existing matrix data into a TYPE(MATRIX) variable with
  !! \c p?csc format.
  !!
  !! @param[inout] m_name      The matrix to be allocated.
  !! @param[in]    idx1        The local row indices, stored as a
  !!                           one-dimensional array.
  !! @param[in]    idx2        The local column pointers, stored as a
  !!                           one-dimensional array.
  !! @param[in]    val         The values of the local matrix elements, stored
  !!                           as a one-dimensional array.
  !! @param[in]    desc        BLACS array descriptor.
  !============================================================================!
  interface m_register_pcsc
     module procedure m_register_pdcsc
     module procedure m_register_pzcsc
  end interface m_register_pcsc
#endif

  !************************************************!

contains

  !============================================================================!
  !> @brief Register matrix (simple dense, serial distribution, real version).
  !============================================================================!
  subroutine m_register_sdden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    dim=shape(A)
    m_name%dim1=dim(1)
    m_name%dim2=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='den'
    m_name%is_serial=.true.
    m_name%is_real=.true.
    m_name%is_sparse=.false.

    m_name%dval => A

    m_name%is_initialized=.true.

  end subroutine m_register_sdden

  !============================================================================!
  !> @brief Register matrix (simple dense, serial distribution, complex
  !!        version).
  !============================================================================!
  subroutine m_register_szden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    dim=shape(A)
    m_name%dim1=dim(1)
    m_name%dim2=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='den'
    m_name%is_serial=.true.
    m_name%is_real=.false.
    m_name%is_sparse=.false.

    m_name%zval => A

    m_name%is_initialized=.true.

  end subroutine m_register_szden

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Register matrix (dense block cyclic, parallel distribution, real
  !!        version).
  !============================================================================!
  subroutine m_register_pddbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='dbc'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.false.

    m_name%dval => A

    m_name%is_initialized=.true.

  end subroutine m_register_pddbc
#endif

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Register matrix (dense block cyclic, parallel distribution, complex
  !!        version).
  !============================================================================!
  subroutine m_register_pzdbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in), target :: A(:,:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(9))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='dbc'
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.false.

    m_name%zval => A

    m_name%is_initialized=.true.

  end subroutine m_register_pzdbc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list/compressed sparse column
  !!        from dense block cyclic, parallel distribution, real version).
  !============================================================================!
  subroutine m_register_pdsp_thre(m_name,A,desc,spm_storage,thre)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.
    m_name%is_initialized=.true.

    call psp_den2sp_m(A,desc,m_name%spm,spm_storage,thre)

  end subroutine m_register_pdsp_thre
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list/compressed sparse column
  !!        from dense block cyclic, parallel distribution, complex version).
  !============================================================================!
  subroutine m_register_pzsp_thre(m_name,A,desc,spm_storage,thre)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.true.

    m_name%is_initialized=.true.

    call psp_den2sp_m(A,desc,m_name%spm,spm_storage,thre)

  end subroutine m_register_pzsp_thre
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list, parallel distribution,
  !!        real version).
  !============================================================================!
  subroutine m_register_pdcoo(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2), i, j, k, l

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
    dim(1)=numroc(m_name%dim1,desc(5),k,0,ms_lap_nprow)
    dim(2)=numroc(m_name%dim2,desc(6),l,0,ms_lap_npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='coo'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,m_name%str_type,m_name%iaux2,ms_lap_nprow,ms_lap_npcol)

    m_name%is_initialized=.true.

  end subroutine m_register_pdcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (sparse coordinate list, parallel distribution,
  !!        complex version).
  !============================================================================!
  subroutine m_register_pzcoo(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2), i, j, k, l

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
    dim(1)=numroc(m_name%dim1,desc(5),k,0,ms_lap_nprow)
    dim(2)=numroc(m_name%dim2,desc(6),l,0,ms_lap_npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='coo'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,m_name%str_type,m_name%iaux2,ms_lap_nprow,ms_lap_npcol)

    m_name%is_initialized=.true.

  end subroutine m_register_pzcoo
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (compressed sparse column, parallel distribution,
  !!        real version).
  !============================================================================!
  subroutine m_register_pdcsc(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2), i, j, k, l

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
    dim(1)=numroc(m_name%dim1,desc(5),k,0,ms_lap_nprow)
    dim(2)=numroc(m_name%dim2,desc(6),l,0,ms_lap_npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='csc'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,m_name%str_type,m_name%iaux2,ms_lap_nprow,ms_lap_npcol)

    m_name%is_initialized=.true.

  end subroutine m_register_pdcsc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Register matrix (compressed sparse column, parallel distribution,
  !!        complex version).
  !============================================================================!
  subroutine m_register_pzcsc(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: dim(2), i, j, k, l

    !**** EXTERNAL ********************************!

    integer, external :: numroc

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
    dim(1)=numroc(m_name%dim1,desc(5),k,0,ms_lap_nprow)
    dim(2)=numroc(m_name%dim2,desc(6),l,0,ms_lap_npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='csc'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,m_name%str_type,m_name%iaux2,ms_lap_nprow,ms_lap_npcol)

    m_name%is_initialized=.true.

  end subroutine m_register_pzcsc
#endif

end module MatrixSwitch_m_register
