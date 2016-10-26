module MatrixSwitch_m_add
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

#ifdef PSP
  subroutine m_add_pdcscpddbcref(A,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%dval=beta*C%dval

    do i=1,A%spm%loc_dim2
       do j=0,A%spm%col_ptr(i+1)-A%spm%col_ptr(i)-1
          l=A%spm%col_ptr(i)+j
          C%dval(A%spm%row_ind(l),i)=C%dval(A%spm%row_ind(l),i)+alpha*A%spm%dval(l)
       end do
    end do

  end subroutine m_add_pdcscpddbcref
#endif

  subroutine m_add_sddenref(A,trA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

    C%dval=beta*C%dval

    if (.not. trA) then
       C%dval=C%dval+alpha*A%dval
    else if (trA) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(j,i)
          end do
       end do
    end if

  end subroutine m_add_sddenref

  subroutine m_add_szdenref(A,tcA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

    C%zval=beta*C%zval

    if (tcA==0) then
       C%zval=C%zval+alpha*A%zval
    else if (tcA==1) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(j,i))
          end do
       end do
    else if (tcA==2) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(j,i)
          end do
       end do
    end if

  end subroutine m_add_szdenref

end module MatrixSwitch_m_add
