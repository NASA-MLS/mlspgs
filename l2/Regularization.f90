! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Regularization
!=============================================================================

! Apply a regularization condition to the Jacobian of the least-squares
! problem for retrieval.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MatrixModule_0, only: M_Absent, M_Banded, MatrixElement_T, CreateBlock
  use MatrixModule_1, only: Matrix_T, CreateBlock
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  Regularize  -----
  subroutine Regularize ( A, Order )
  ! Apply regularization conditions of the form "Differences of degree
  ! Order in dX are approximately zero" to the blocks of A.
  ! Order specifies the order of polynomial.  If size(Order) == 1, it
  ! applies to all of the blocks of A.  Otherwise, size(Order) shall
  ! equal A%col%nb (except the extra block, if present, is set to zero).
  ! It is necessary that the total number of rows in the first row block
  ! of A be large enough to accomodate the regularization -- roughly at
  ! least (number of columns of a) - min(order).  The number of columns
  ! of each block shall be at least one more than the regularization
  ! order for that block.

    type(Matrix_T), intent(inout) :: A
    integer, intent(in), dimension(:) :: Order

    integer :: I, J                ! Subscripts, Loop inductors
    integer :: IB                  ! Which block is being regularized
    integer :: NB                  ! Number of column blocks of A
    integer :: NCOL                ! Number of columns in a block of A
    integer :: NROW                ! Number of rows in first row block of A
    integer :: NV                  ! Next element in VALUES component
    integer :: Ord                 ! Order for the current block
    real(r8), pointer, dimension(:) :: Pascal ! Regularization coefficients
    integer :: Rows                ! Next row to use
    real :: S ! NOT R8!            ! Sign of regularization coefficient.

    ! The regularization order is never going to be so large that
    ! SBINOM cannot calculate coefficients.
!     interface
!       real function SBINOM ( N, K )
!         integer, intent(in) :: N, K
!       end function SBINOM
!     end interface

    nb = a%col%nb
    if ( a%col%extra ) nb = nb - 1
    if ( size(order) /= 1 .and. size(order) /= nb ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, "Wrong size for Order" )

    nrow = a%row%nelts(1)
    rows = 1
    do ib = 1, nb
      ord = order(1)
      if ( size(order) > 1 ) ord = order(ib)
      ncol = a%col%nelts(ib)
      if ( rows + ncol - ord - 1 > a%row%nelts(1) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Not enough rows in the matrix to do regularization" )
      if ( ord == 0 ) then
        call createBlock ( a%block(1,ib), nrow, ncol,  m_absent )
      else
        if ( ncol < ord+1 ) call MLSMessage ( MLSMSG_Error, &
          & moduleName, "Not enough columns for specified regularization order" )
        ! Construct regularization coefficients using Pascal's triangle
        call allocate_test ( Pascal, ord, "Pascal's triangle", moduleName, &
          & lowBound=0 )
        s = -1.0
        do i = 0, ord
          pascal(i) = s * sbinom(ord,i)
          s = - s
        end do
        call createBlock ( a%block(1,ib), nrow, ncol, m_banded, &
          & (ncol-ord)*(ord+1) )

        ! Each row has the binomial coefficients.  Therefore, each column
        ! has the binomial coefficients in referse order (which wouldn't
        ! matter if the signs didn't alternate).  Except the first and last
        ! ord columns have 1, 2, ..., ord and ord, ord-1, ..., 1 elements.

        nv = 1
        do i = 1, ord
          a%block(1,ib)%r1(i) = rows
          a%block(1,ib)%r2(i) = nv
          a%block(1,ib)%values(nv:nv+i-1,1) = pascal(ord-i+1:ord)
          nv = nv + i
        end do
        do i = ord+1, ncol-ord
          a%block(1,ib)%r1(i) = rows
          a%block(1,ib)%r2(i) = nv
          a%block(1,ib)%values(nv:nv+ord,1) = pascal
          nv = nv + ord + 1
          rows = rows + 1
        end do
        j = ord-1
        do i = ncol-ord+1, ncol
          a%block(1,ib)%r1(i) = a%block(1,ib)%r1(i-1)+1
          a%block(1,ib)%r2(i) = nv
          a%block(1,ib)%values(nv:nv+j,1) = pascal(0:j)
          nv = nv + j + 1
        end do
        call deallocate_test ( Pascal, "Pascal's triangle", moduleName )
      end if
      rows = rows + ncol - ord
    end do
    if ( a%col%extra ) &
      & call createBlock ( a%block(1,a%col%nb), nrow, 1,  m_absent )
  end subroutine Regularize

  real function SBINOM ( N, K )
    integer, intent(in) :: N, K
  end function SBINOM

end module Regularization



! $Log$
! Revision 2.2  2001/06/02 16:58:46  livesey
! Temporary fix to let it compile before Paul has a go.
! (commented out interface to sbinom and wrote empty routine instead).
!
! Revision 2.1  2001/06/02 01:40:29  vsnyder
! Initial commit
!
