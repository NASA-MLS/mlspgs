! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Regularization
!=============================================================================

! Apply Tikhonov regularization to the matrix of the least-squares
! problem for retrieval.

  use MatrixModule_0, only: M_Absent, M_Banded, MatrixElement_T, CreateBlock
  use MatrixModule_1, only: Matrix_T, CreateBlock
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, parameter, public :: MaxRegOrd = 33 ! Maximum regularization
  ! order.  33!/(16!)**2 < HUGE(0) < 34!/(17!)**2 for 31-bit integers.
  !         66!/(33!)**2 < HUGE(0) < 67!/(34!**2) for 63-bit integers.

contains

  ! -------------------------------------------------  Regularize  -----
  subroutine Regularize ( A, Order )
  !{Apply Tikhonov regularization conditions of the form ``Differences of
  ! degree Order in dX are approximately zero'' to the blocks of A. Order
  ! specifies the order of polynomial.  The regularization condition is
  ! imposed by an Order'th degree difference operator, which has binomial
  ! coefficients with alternating sign.  If size(Order) == 1, it applies to
  ! all of the blocks of A.  Otherwise, size(Order) shall equal A\%col\%nb
  ! (except the extra block, if present, is set to zero). It is necessary
  ! that the total number of rows in the first row block of A be large
  ! enough to accomodate the regularization~-- roughly at least (number of
  ! columns of a)~- min(order).  The number of columns of each block shall
  ! be at least one more than the regularization order for that block.

    type(Matrix_T), intent(inout) :: A
    integer, intent(in), dimension(:) :: Order

    integer, dimension(0:maxRegOrd) :: C ! Binomial regularization coefficients
    integer :: I, J                ! Subscripts, Loop inductors
    integer :: IB                  ! Which block is being regularized
    character(len=2) :: MSG        ! In case of an error message
    integer :: NB                  ! Number of column blocks of A
    integer :: NCOL                ! Number of columns in a block of A
    integer :: NROW                ! Number of rows in first row block of A
    integer :: NV                  ! Next element in VALUES component
    integer :: Ord                 ! Order for the current block
    integer :: Rows                ! Next row to use
    integer :: S                   ! Sign of regularization coefficient.

    ! The regularization order is never going to be so large that
    ! the coefficients cannot be represented by integers.

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
        if ( ord > maxRegOrd ) then
          write ( msg, '(i2)' ) maxRegOrd
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Regularization order exceeds " // trim(adjustl(msg)) )
        end if
        !{ Calculate binomial coefficients $C_i^n = \frac{n!}{i! (n-i)!}$
        !  by the recursion
        !  $C_0^n = 1\text{, } C_i^n = (n-i+1) C_{i-1}^n / i$.
        ! Notice that $C_i^n = C_{n-i}^n$, so we only need
        ! to go halfway through the array.
        c(0) = 1
        c(ord) = 1
        do i = 1, ord / 2
          c(i) = ( (ord-i+1) * c(i-1) ) / i
          c(ord-i) = c(i)
        end do
        s = -1
        do i = 0, ord ! Now alternate the signs
          c(i) = s * c(i)
          s = -s
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
          a%block(1,ib)%values(nv:nv+i-1,1) = c(ord-i+1:ord)
          nv = nv + i
        end do
        do i = ord+1, ncol-ord
          a%block(1,ib)%r1(i) = rows
          a%block(1,ib)%r2(i) = nv
          a%block(1,ib)%values(nv:nv+ord,1) = c
          nv = nv + ord + 1
          rows = rows + 1
        end do
        j = ord-1
        do i = ncol-ord+1, ncol
          a%block(1,ib)%r1(i) = a%block(1,ib)%r1(i-1)+1
          a%block(1,ib)%r2(i) = nv
          a%block(1,ib)%values(nv:nv+j,1) = c(0:j)
          nv = nv + j + 1
        end do
      end if
      rows = rows + ncol - ord
    end do
    if ( a%col%extra ) &
      & call createBlock ( a%block(1,a%col%nb), nrow, 1,  m_absent )
  end subroutine Regularize

end module Regularization

! $Log$
! Revision 2.3  2001/06/22 00:41:54  vsnyder
! Replace use of SBINOM by in-line calculation of binomial coefficients
!
! Revision 2.2  2001/06/02 16:58:46  livesey
! Temporary fix to let it compile before Paul has a go.
! (commented out interface to sbinom and wrote empty routine instead).
!
! Revision 2.1  2001/06/02 01:40:29  vsnyder
! Initial commit
!
