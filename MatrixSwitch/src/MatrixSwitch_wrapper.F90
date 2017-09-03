#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!==============================================================================!
!> @brief MatrixSwitch C bindings.
!==============================================================================!
module MatrixSwitch_wrapper
  use MatrixSwitch_ops, only: dp
  use MatrixSwitch_wrapper_params
  use MatrixSwitch, only: &
    mm_multiply_orig => mm_multiply, &
    m_add_orig => m_add, &
    m_trace_orig => m_trace, &
    mm_trace_orig => mm_trace, &
    m_scale_orig => m_scale, &
    m_set_orig => m_set, &
    m_set_element_orig => m_set_element, &
    m_get_element_orig => m_get_element, &
    m_register_sden_orig => m_register_sden, &
    m_allocate_orig => m_allocate, &
    m_deallocate_orig => m_deallocate, &
    m_copy_orig => m_copy, &
    m_convert_orig => m_convert
#ifdef HAVE_MPI
  use MatrixSwitch, only: &
    m_register_pdbc_orig => m_register_pdbc, &
#ifdef HAVE_SCALAPACK
    ms_scalapack_setup, &
#endif
    ms_lap_icontxt
#endif
#ifdef HAVE_PSPBLAS
  use MatrixSwitch, only: &
    m_copy_external_pdbcpcoo_orig => m_copy_external_pdbcpcoo, &
    m_copy_external_pdbcpcsc_orig => m_copy_external_pdbcpcsc, &
    m_register_pcoo_orig => m_register_pcoo, &
    m_register_pcsc_orig => m_register_pcsc
#endif

  implicit none

  private

  !**** INTERFACES ********************************!

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_initialized).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_initialized
  !============================================================================!
  interface ms_is_initialized
     module procedure ms_is_initialized_key
     module procedure ms_is_initialized_index
  end interface ms_is_initialized

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_serial).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_serial
  !============================================================================!
  interface ms_is_serial
     module procedure ms_is_serial_key
     module procedure ms_is_serial_index
  end interface ms_is_serial

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_real).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_real
  !============================================================================!
  interface ms_is_real
     module procedure ms_is_real_key
     module procedure ms_is_real_index
  end interface ms_is_real

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_square).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_square
  !============================================================================!
  interface ms_is_square
     module procedure ms_is_square_key
     module procedure ms_is_square_index
  end interface ms_is_square

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_is_sparse).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_is_sparse
  !============================================================================!
  interface ms_is_sparse
     module procedure ms_is_sparse_key
     module procedure ms_is_sparse_index
  end interface ms_is_sparse

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_dim1).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_dim1
  !============================================================================!
  interface ms_dim1
     module procedure ms_dim1_key
     module procedure ms_dim1_index
  end interface ms_dim1

  !============================================================================!
  !> @brief Unpack matrix data type (\p ms_dim2).
  !!
  !! @param[in]  m_name The matrix to unpack.
  !! @return            \p ms_dim2
  !============================================================================!
  interface ms_dim2
     module procedure ms_dim2_key
     module procedure ms_dim2_index
  end interface ms_dim2

  !============================================================================!
  !> @brief Wrapper to allocate matrix.
  !============================================================================!
  interface m_allocate
     module procedure m_allocate_key
     module procedure m_allocate_index
  end interface m_allocate

  !============================================================================!
  !> @brief Wrapper to deallocate matrix.
  !============================================================================!
  interface m_deallocate
     module procedure m_deallocate_key
     module procedure m_deallocate_index
  end interface m_deallocate

  !============================================================================!
  !> @brief Wrapper to copy matrix.
  !============================================================================!
  interface m_copy
     module procedure m_copy_key
     module procedure m_copy_index
  end interface m_copy

  !============================================================================!
  !> @brief Wrapper to convert matrix format.
  !============================================================================!
  interface m_convert
     module procedure m_convert_key
     module procedure m_convert_index
  end interface m_convert

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication.
  !============================================================================!
  interface mm_multiply
     module procedure mm_dmultiply_key
     module procedure mm_dmultiply_index
     module procedure mm_zmultiply_key
     module procedure mm_zmultiply_index
  end interface mm_multiply

  !============================================================================!
  !> @brief Wrapper to matrix addition.
  !============================================================================!
  interface m_add
     module procedure m_dadd_key
     module procedure m_dadd_index
     module procedure m_zadd_key
     module procedure m_zadd_index
  end interface m_add

  !============================================================================!
  !> @brief Wrapper to matrix trace.
  !============================================================================!
  interface m_trace
     module procedure m_dtrace_key
     module procedure m_dtrace_index
     module procedure m_ztrace_key
     module procedure m_ztrace_index
  end interface m_trace

  !============================================================================!
  !> @brief Wrapper to matrix product trace.
  !============================================================================!
  interface mm_trace
     module procedure mm_dtrace_key
     module procedure mm_dtrace_index
     module procedure mm_ztrace_key
     module procedure mm_ztrace_index
  end interface mm_trace

  !============================================================================!
  !> @brief Wrapper to scale matrix.
  !============================================================================!
  interface m_scale
     module procedure m_dscale_key
     module procedure m_dscale_index
     module procedure m_zscale_key
     module procedure m_zscale_index
  end interface m_scale

  !============================================================================!
  !> @brief Wrapper to set matrix.
  !============================================================================!
  interface m_set
     module procedure m_dset_key
     module procedure m_dset_index
     module procedure m_zset_key
     module procedure m_zset_index
  end interface m_set

  !============================================================================!
  !> @brief Wrapper to set matrix element.
  !============================================================================!
  interface m_set_element
     module procedure m_dset_element_key
     module procedure m_dset_element_index
     module procedure m_zset_element_key
     module procedure m_zset_element_index
  end interface m_set_element

  !============================================================================!
  !> @brief Wrapper to get matrix element.
  !============================================================================!
  interface m_get_element
     module procedure m_dget_element_key
     module procedure m_dget_element_index
     module procedure m_zget_element_key
     module procedure m_zget_element_index
  end interface m_get_element

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution).
  !============================================================================!
  interface m_register_sden
     module procedure m_register_sdden_key
     module procedure m_register_sdden_index
     module procedure m_register_szden_key
     module procedure m_register_szden_index
  end interface m_register_sden

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pdbc
     module procedure m_register_pddbc_key
     module procedure m_register_pddbc_index
     module procedure m_register_pzdbc_key
     module procedure m_register_pzdbc_index
  end interface m_register_pdbc
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution).
  !============================================================================!
  interface m_copy_external_pdbcpcoo
     module procedure m_copy_external_pddbcpdcoo_key
     module procedure m_copy_external_pddbcpdcoo_index
     module procedure m_copy_external_pzdbcpzcoo_key
     module procedure m_copy_external_pzdbcpzcoo_index
  end interface m_copy_external_pdbcpcoo

  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution).
  !============================================================================!
  interface m_copy_external_pdbcpcsc
     module procedure m_copy_external_pddbcpdcsc_key
     module procedure m_copy_external_pddbcpdcsc_index
     module procedure m_copy_external_pzdbcpzcsc_key
     module procedure m_copy_external_pzdbcpzcsc_index
  end interface m_copy_external_pdbcpcsc

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pcoo
     module procedure m_register_pdcoo_key
     module procedure m_register_pdcoo_index
     module procedure m_register_pzcoo_key
     module procedure m_register_pzcoo_index
  end interface m_register_pcoo

  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution).
  !============================================================================!
  interface m_register_pcsc
     module procedure m_register_pdcsc_key
     module procedure m_register_pdcsc_index
     module procedure m_register_pzcsc_key
     module procedure m_register_pzcsc_index
  end interface m_register_pcsc
#endif

  !************************************************!

  public :: ms_wrapper_open
  public :: ms_wrapper_close
  public :: ms_is_initialized
  public :: ms_is_serial
  public :: ms_is_real
  public :: ms_is_square
  public :: ms_is_sparse
  public :: ms_dim1
  public :: ms_dim2
  public :: mm_multiply
  public :: m_add
  public :: m_trace
  public :: mm_trace
  public :: m_scale
  public :: m_set
  public :: m_set_element
  public :: m_get_element
  public :: m_register_sden
  public :: m_allocate
  public :: m_deallocate
  public :: m_copy
  public :: m_convert
#ifdef HAVE_MPI
  public :: m_register_pdbc
#ifdef HAVE_SCALAPACK
  public :: ms_scalapack_setup
#endif
  public :: ms_lap_icontxt
#endif
#ifdef HAVE_PSPBLAS
  public :: m_copy_external_pdbcpcoo
  public :: m_copy_external_pdbcpcsc
  public :: m_register_pcoo
  public :: m_register_pcsc
#endif

contains

  !============================================================================!
  !> @brief MatrixSwitch wrapper setup.
  !!
  !! Sets up the MatrixSwitch wrapper with a fixed number of matrices, each of
  !! which is optionally identified by a unique key.
  !!
  !! @param[in] num_matrices The number of matrices to store in the wrapper.
  !! @param[in] keys         The optional array of keys for accessing the
  !!                         stored matrices. Each key cannot exceed
  !!                         \p max_key_length characters.
  !============================================================================!
  subroutine ms_wrapper_open(num_matrices,keys)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: keys(num_matrices)

    integer, intent(in) :: num_matrices

    !**********************************************!

    allocate(ms_matrices(num_matrices))
    ms_num_matrices=num_matrices
    
    if (present(keys)) then
      allocate(ms_keys(num_matrices))
      ms_keys=keys
    end if

  end subroutine ms_wrapper_open

  !============================================================================!
  !> @brief Close the MatrixSwitch wrapper.
  !!
  !! All MatrixSwitch matrices stored in the wrapper are deallocated, and the
  !! map is deallocated.
  !============================================================================!
  subroutine ms_wrapper_close()
    implicit none

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

    if (allocated(ms_keys)) deallocate(ms_keys)
    if (allocated(ms_matrices)) then
      do i=1,ms_num_matrices
        call m_deallocate_orig(ms_matrices(i))
      end do
      deallocate(ms_matrices)
    end if

  end subroutine ms_wrapper_close

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  logical function ms_is_initialized_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_initialized_key=ms_is_initialized_index(ms_lookup(m_name))

  end function ms_is_initialized_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  logical function ms_is_initialized_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_is_initialized_index=ms_matrices(m_name)%is_initialized

  end function ms_is_initialized_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  logical function ms_is_serial_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_serial_key=ms_is_serial_index(ms_lookup(m_name))

  end function ms_is_serial_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  logical function ms_is_serial_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_is_serial_index=ms_matrices(m_name)%is_serial

  end function ms_is_serial_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  logical function ms_is_real_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_real_key=ms_is_real_index(ms_lookup(m_name))

  end function ms_is_real_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  logical function ms_is_real_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_is_real_index=ms_matrices(m_name)%is_real

  end function ms_is_real_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  logical function ms_is_square_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_square_key=ms_is_square_index(ms_lookup(m_name))

  end function ms_is_square_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  logical function ms_is_square_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_is_square_index=ms_matrices(m_name)%is_square

  end function ms_is_square_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  logical function ms_is_sparse_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_is_sparse_key=ms_is_sparse_index(ms_lookup(m_name))

  end function ms_is_sparse_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  logical function ms_is_sparse_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_is_sparse_index=ms_matrices(m_name)%is_sparse

  end function ms_is_sparse_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  integer function ms_dim1_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_dim1_key=ms_dim1_index(ms_lookup(m_name))

  end function ms_dim1_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  integer function ms_dim1_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_dim1_index=ms_matrices(m_name)%dim1

  end function ms_dim1_index

  !============================================================================!
  !> @brief Unpack matrix data type (key version).
  !============================================================================!
  integer function ms_dim2_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    ms_dim2_key=ms_dim2_index(ms_lookup(m_name))

  end function ms_dim2_key

  !============================================================================!
  !> @brief Unpack matrix data type (index version).
  !============================================================================!
  integer function ms_dim2_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    ms_dim2_index=ms_matrices(m_name)%dim2

  end function ms_dim2_index

  !============================================================================!
  !> @brief Wrapper to allocate matrix (key version).
  !============================================================================!
  subroutine m_allocate_key(m_name,i,j,label)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label
    character(*), intent(in) :: m_name

    integer, intent(in) :: i
    integer, intent(in) :: j

    !**********************************************!

    call m_allocate_index(ms_lookup(m_name),i,j,label)

  end subroutine m_allocate_key

  !============================================================================!
  !> @brief Wrapper to allocate matrix (index version).
  !============================================================================!
  subroutine m_allocate_index(m_name,i,j,label)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label
    integer, intent(in) :: m_name

    integer, intent(in) :: i
    integer, intent(in) :: j

    !**********************************************!

    call m_allocate_orig(ms_matrices(m_name),i,j,label)

  end subroutine m_allocate_index

  !============================================================================!
  !> @brief Wrapper to deallocate matrix (key version).
  !============================================================================!
  subroutine m_deallocate_key(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_deallocate_index(ms_lookup(m_name))

  end subroutine m_deallocate_key

  !============================================================================!
  !> @brief Wrapper to deallocate matrix (index version).
  !============================================================================!
  subroutine m_deallocate_index(m_name)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: m_name

    !**********************************************!

    call m_deallocate_orig(ms_matrices(m_name))

  end subroutine m_deallocate_index

  !============================================================================!
  !> @brief Wrapper to copy matrix (key version).
  !============================================================================!
  subroutine m_copy_key(m_name,A,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: A
    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_index(ms_lookup(m_name),ms_lookup(A),label,threshold,threshold_is_soft)

  end subroutine m_copy_key

  !============================================================================!
  !> @brief Wrapper to copy matrix (index version).
  !============================================================================!
  subroutine m_copy_index(m_name,A,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    integer, intent(in) :: A
    integer, intent(in) :: m_name

    !**********************************************!

    call m_copy_orig(ms_matrices(m_name),ms_matrices(A),label,threshold,threshold_is_soft)

  end subroutine m_copy_index

  !============================================================================!
  !> @brief Wrapper to convert matrix format (key version).
  !============================================================================!
  subroutine m_convert_key(m_name,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_convert_index(ms_lookup(m_name),label,threshold,threshold_is_soft)

  end subroutine m_convert_key

  !============================================================================!
  !> @brief Wrapper to convert matrix format (index version).
  !============================================================================!
  subroutine m_convert_index(m_name,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label

    logical, intent(in), optional :: threshold_is_soft

    real(dp), intent(in), optional :: threshold

    integer, intent(in) :: m_name

    !**********************************************!

    call m_convert_orig(ms_matrices(m_name),label,threshold,threshold_is_soft)

  end subroutine m_convert_index

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (key real version).
  !============================================================================!
  subroutine mm_dmultiply_key(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: B
    character(*), intent(in) :: C

    !**********************************************!

    call mm_dmultiply_index(ms_lookup(A),opA,ms_lookup(B),opB,ms_lookup(C),alpha,beta,label)

  end subroutine mm_dmultiply_key

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (index real version).
  !============================================================================!
  subroutine mm_dmultiply_index(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    integer, intent(in) :: A
    integer, intent(in) :: B
    integer, intent(in) :: C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(A),opA,ms_matrices(B),opB,ms_matrices(C),alpha,beta,label)

  end subroutine mm_dmultiply_index

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (key complex version).
  !============================================================================!
  subroutine mm_zmultiply_key(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: B
    character(*), intent(in) :: C

    !**********************************************!

    call mm_zmultiply_index(ms_lookup(A),opA,ms_lookup(B),opB,ms_lookup(C),alpha,beta,label)

  end subroutine mm_zmultiply_key

  !============================================================================!
  !> @brief Wrapper to matrix-matrix multiplication (index complex version).
  !============================================================================!
  subroutine mm_zmultiply_index(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(1), intent(in) :: opB
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    integer, intent(in) :: A
    integer, intent(in) :: B
    integer, intent(in) :: C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(A),opA,ms_matrices(B),opB,ms_matrices(C),alpha,beta,label)

  end subroutine mm_zmultiply_index

  !============================================================================!
  !> @brief Wrapper to matrix addition (key real version).
  !============================================================================!
  subroutine m_dadd_key(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: C

    !**********************************************!

    call m_dadd_index(ms_lookup(A),opA,ms_lookup(C),alpha,beta,label)

  end subroutine m_dadd_key

  !============================================================================!
  !> @brief Wrapper to matrix addition (index real version).
  !============================================================================!
  subroutine m_dadd_index(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    integer, intent(in) :: A
    integer, intent(in) :: C

    !**********************************************!

    call m_add_orig(ms_matrices(A),opA,ms_matrices(C),alpha,beta,label)

  end subroutine m_dadd_index

  !============================================================================!
  !> @brief Wrapper to matrix addition (key complex version).
  !============================================================================!
  subroutine m_zadd_key(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: A
    character(*), intent(in) :: C

    !**********************************************!

    call m_zadd_index(ms_lookup(A),opA,ms_lookup(C),alpha,beta,label)

  end subroutine m_zadd_key

  !============================================================================!
  !> @brief Wrapper to matrix addition (index complex version).
  !============================================================================!
  subroutine m_zadd_index(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    integer, intent(in) :: A
    integer, intent(in) :: C

    !**********************************************!

    call m_add_orig(ms_matrices(A),opA,ms_matrices(C),alpha,beta,label)

  end subroutine m_zadd_index

  !============================================================================!
  !> @brief Wrapper to matrix trace (key real version).
  !============================================================================!
  subroutine m_dtrace_key(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_dtrace_index(ms_lookup(A),alpha,label)

  end subroutine m_dtrace_key

  !============================================================================!
  !> @brief Wrapper to matrix trace (index real version).
  !============================================================================!
  subroutine m_dtrace_index(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: A

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(A),alpha,label)

  end subroutine m_dtrace_index

  !============================================================================!
  !> @brief Wrapper to matrix trace (key complex version).
  !============================================================================!
  subroutine m_ztrace_key(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_ztrace_index(ms_lookup(A),alpha,label)

  end subroutine m_ztrace_key

  !============================================================================!
  !> @brief Wrapper to matrix trace (index complex version).
  !============================================================================!
  subroutine m_ztrace_index(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: A

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(A),alpha,label)

  end subroutine m_ztrace_index

  !============================================================================!
  !> @brief Wrapper to matrix product trace (key real version).
  !============================================================================!
  subroutine mm_dtrace_key(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A
    character(*), intent(in) :: B

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call mm_dtrace_index(ms_lookup(A),ms_lookup(B),alpha,label)

  end subroutine mm_dtrace_key

  !============================================================================!
  !> @brief Wrapper to matrix product trace (index real version).
  !============================================================================!
  subroutine mm_dtrace_index(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: A
    integer, intent(in) :: B

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(A),ms_matrices(B),alpha,label)

  end subroutine mm_dtrace_index

  !============================================================================!
  !> @brief Wrapper to matrix product trace (key complex version).
  !============================================================================!
  subroutine mm_ztrace_key(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    character(*), intent(in) :: A
    character(*), intent(in) :: B

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call mm_ztrace_index(ms_lookup(A),ms_lookup(B),alpha,label)

  end subroutine mm_ztrace_key

  !============================================================================!
  !> @brief Wrapper to matrix product trace (index complex version).
  !============================================================================!
  subroutine mm_ztrace_index(A,B,alpha,label)
    implicit none
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: A
    integer, intent(in) :: B

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(A),ms_matrices(B),alpha,label)

  end subroutine mm_ztrace_index

  !============================================================================!
  !> @brief Wrapper to scale matrix (key real version).
  !============================================================================!
  subroutine m_dscale_key(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_dscale_index(ms_lookup(C),beta,label)

  end subroutine m_dscale_key

  !============================================================================!
  !> @brief Wrapper to scale matrix (index real version).
  !============================================================================!
  subroutine m_dscale_index(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    real(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_scale_orig(ms_matrices(C),beta,label)

  end subroutine m_dscale_index

  !============================================================================!
  !> @brief Wrapper to scale matrix (key complex version).
  !============================================================================!
  subroutine m_zscale_key(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_zscale_index(ms_lookup(C),beta,label)

  end subroutine m_zscale_key

  !============================================================================!
  !> @brief Wrapper to scale matrix (index complex version).
  !============================================================================!
  subroutine m_zscale_index(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_scale_orig(ms_matrices(C),beta,label)

  end subroutine m_zscale_index

  !============================================================================!
  !> @brief Wrapper to set matrix (key real version).
  !============================================================================!
  subroutine m_dset_key(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_dset_index(ms_lookup(C),seC,alpha,beta,label)

  end subroutine m_dset_key

  !============================================================================!
  !> @brief Wrapper to set matrix (index real version).
  !============================================================================!
  subroutine m_dset_index(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_set_orig(ms_matrices(C),seC,alpha,beta,label)

  end subroutine m_dset_index

  !============================================================================!
  !> @brief Wrapper to set matrix (key complex version).
  !============================================================================!
  subroutine m_zset_key(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_zset_index(ms_lookup(C),seC,alpha,beta,label)

  end subroutine m_zset_key

  !============================================================================!
  !> @brief Wrapper to set matrix (index complex version).
  !============================================================================!
  subroutine m_zset_index(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC
    character(3), intent(in), optional :: label

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_set_orig(ms_matrices(C),seC,alpha,beta,label)

  end subroutine m_zset_index

  !============================================================================!
  !> @brief Wrapper to set matrix element (key real version).
  !============================================================================!
  subroutine m_dset_element_key(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_dset_element_index(ms_lookup(C),i,j,alpha,beta,label)

  end subroutine m_dset_element_key

  !============================================================================!
  !> @brief Wrapper to set matrix element (index real version).
  !============================================================================!
  subroutine m_dset_element_index(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_set_element_orig(ms_matrices(C),i,j,alpha,beta,label)

  end subroutine m_dset_element_index

  !============================================================================!
  !> @brief Wrapper to set matrix element (key complex version).
  !============================================================================!
  subroutine m_zset_element_key(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    character(*), intent(in) :: C

    !**********************************************!

    call m_zset_element_index(ms_lookup(C),i,j,alpha,beta,label)

  end subroutine m_zset_element_key

  !============================================================================!
  !> @brief Wrapper to set matrix element (index complex version).
  !============================================================================!
  subroutine m_zset_element_index(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    integer, intent(in) :: C

    !**********************************************!

    call m_set_element_orig(ms_matrices(C),i,j,alpha,beta,label)

  end subroutine m_zset_element_index

  !============================================================================!
  !> @brief Wrapper to get matrix element (key real version).
  !============================================================================!
  subroutine m_dget_element_key(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    character(*), intent(in) :: C

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_dget_element_index(ms_lookup(C),i,j,alpha,label)

  end subroutine m_dget_element_key

  !============================================================================!
  !> @brief Wrapper to get matrix element (index real version).
  !============================================================================!
  subroutine m_dget_element_index(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    integer, intent(in) :: C

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(C),i,j,alpha,label)

  end subroutine m_dget_element_index

  !============================================================================!
  !> @brief Wrapper to get matrix element (key complex version).
  !============================================================================!
  subroutine m_zget_element_key(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    character(*), intent(in) :: C

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_zget_element_index(ms_lookup(C),i,j,alpha,label)

  end subroutine m_zget_element_key

  !============================================================================!
  !> @brief Wrapper to get matrix element (index complex version).
  !============================================================================!
  subroutine m_zget_element_index(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label

    integer, intent(in) :: i
    integer, intent(in) :: j

    integer, intent(in) :: C

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(C),i,j,alpha,label)

  end subroutine m_zget_element_index

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution, key
  !!        real version).
  !============================================================================!
  subroutine m_register_sdden_key(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_sdden_index(ms_lookup(m_name),A)

  end subroutine m_register_sdden_key

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution,
  !!        index real version).
  !============================================================================!
  subroutine m_register_sdden_index(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_sden_orig(ms_matrices(m_name),A)

  end subroutine m_register_sdden_index

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution,
  !!        key complex version).
  !============================================================================!
  subroutine m_register_szden_key(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_szden_index(ms_lookup(m_name),A)

  end subroutine m_register_szden_key

  !============================================================================!
  !> @brief Wrapper to register matrix (simple dense, serial distribution,
  !!        index complex version).
  !============================================================================!
  subroutine m_register_szden_index(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_sden_orig(ms_matrices(m_name),A)

  end subroutine m_register_szden_index

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, key real version).
  !============================================================================!
  subroutine m_register_pddbc_key(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pddbc_index(ms_lookup(m_name),A,desc)

  end subroutine m_register_pddbc_key
#endif

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, index real version).
  !============================================================================!
  subroutine m_register_pddbc_index(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(m_name),A,desc)

  end subroutine m_register_pddbc_index
#endif

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, key complex version).
  !============================================================================!
  subroutine m_register_pzdbc_key(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pzdbc_index(ms_lookup(m_name),A,desc)

  end subroutine m_register_pzdbc_key
#endif

#ifdef HAVE_MPI
  !============================================================================!
  !> @brief Wrapper to register matrix (dense block cyclic, parallel
  !!        distribution, index complex version).
  !============================================================================!
  subroutine m_register_pzdbc_index(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9)

    complex(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(m_name),A,desc)

  end subroutine m_register_pzdbc_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, key real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcoo_key(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pddbcpdcoo_index(ms_lookup(m_name),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcoo_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, index real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcoo_index(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    integer, intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcoo_orig(ms_matrices(m_name),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcoo_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, key complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcoo_key(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pzdbcpzcoo_index(ms_lookup(m_name),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcoo_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (sparse coordinate list from dense
  !!        block cyclic, parallel distribution, index complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcoo_index(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcoo_orig(ms_matrices(m_name),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcoo_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, key real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcsc_key(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pddbcpdcsc_index(ms_lookup(m_name),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcsc_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, index real version).
  !============================================================================!
  subroutine m_copy_external_pddbcpdcsc_index(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), target :: A(:,:)
    real(dp), intent(in), optional :: threshold

    integer, intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcsc_orig(ms_matrices(m_name),A,desc,threshold)

  end subroutine m_copy_external_pddbcpdcsc_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, key complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcsc_key(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pzdbcpzcsc_index(ms_lookup(m_name),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcsc_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to copy external matrix (compressed sparse column from
  !!        dense block cyclic, parallel distribution, index complex version).
  !============================================================================!
  subroutine m_copy_external_pzdbcpzcsc_index(m_name,A,desc,threshold)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)

    real(dp), intent(in), optional :: threshold

    complex(dp), intent(in), target :: A(:,:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_copy_external_pdbcpcsc_orig(ms_matrices(m_name),A,desc,threshold)

  end subroutine m_copy_external_pzdbcpzcsc_index
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, key real version).
  !============================================================================!
  subroutine m_register_pdcoo_key(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pdcoo_index(ms_lookup(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pdcoo_key
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, index real version).
  !============================================================================!
  subroutine m_register_pdcoo_index(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pcoo_orig(ms_matrices(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pdcoo_index
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, key complex version).
  !============================================================================!
  subroutine m_register_pzcoo_key(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pzcoo_index(ms_lookup(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pzcoo_key
#endif

#ifdef HAVE_PSPBLAS

  !============================================================================!
  !> @brief Wrapper to register matrix (sparse coordinate list, parallel
  !!        distribution, index complex version).
  !============================================================================!
  subroutine m_register_pzcoo_index(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pcoo_orig(ms_matrices(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pzcoo_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, key real version).
  !============================================================================!
  subroutine m_register_pdcsc_key(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pdcsc_index(ms_lookup(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pdcsc_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, index real version).
  !============================================================================!
  subroutine m_register_pdcsc_index(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    real(dp), intent(in), target :: val(:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pcsc_orig(ms_matrices(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pdcsc_index
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, key complex version).
  !============================================================================!
  subroutine m_register_pzcsc_key(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    character(*), intent(in) :: m_name

    !**********************************************!

    call m_register_pzcsc_index(ms_lookup(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pzcsc_key
#endif

#ifdef HAVE_PSPBLAS
  !============================================================================!
  !> @brief Wrapper to register matrix (compressed sparse column, parallel
  !!        distribution, index complex version).
  !============================================================================!
  subroutine m_register_pzcsc_index(m_name,idx1,idx2,val,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9)
    integer, intent(in), target :: idx1(:)
    integer, intent(in), target :: idx2(:)

    complex(dp), intent(in), target :: val(:)

    integer, intent(in) :: m_name

    !**********************************************!

    call m_register_pcsc_orig(ms_matrices(m_name),idx1,idx2,val,desc)

  end subroutine m_register_pzcsc_index
#endif

end module MatrixSwitch_wrapper
