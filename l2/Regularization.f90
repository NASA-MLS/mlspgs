! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Regularization
!=============================================================================

! Apply Tikhonov regularization to the matrix of the least-squares
! problem for retrieval.

  use Expr_M, only: EXPR
  use Intrinsic, only: PHYQ_DIMENSIONLESS
  use Lexer_Core, only: Print_Source
  use MatrixModule_0, only: M_Absent, M_Banded, MatrixElement_T, CreateBlock
  use MatrixModule_1, only: Matrix_T, CreateBlock
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Output_M, only: Output
  use Tree, only: DECORATION, NSONS, SOURCE_REF, SUBTREE

  private

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

  public :: Regularize

  integer, parameter :: FieldSizes = 1   ! size(regOrders) /= size(regQuants)
  integer, parameter :: OrderTooBig = FieldSizes + 1
  integer, parameter :: RegQuantsReq = OrderTooBig + 1 ! RegQuants required
  integer, parameter :: TooFewCols = RegQuantsReq + 1  ! Won't fit
  integer, parameter :: TooFewRows = TooFewCols + 1    ! Won't fit
  integer, parameter :: Unitless = TooFewRows + 1      ! Orders must be unitless

  integer :: Error       ! non-zero if an error occurs

contains

  ! -------------------------------------------------  Regularize  -----
  subroutine Regularize ( A, Orders, Quants, Weight )

  !{Apply Tikhonov regularization conditions of the form ``Differences of
  ! degree Order in dX are approximately zero'' to the blocks of A, where
  ! Order is the order of a polynomial.  The regularization condition is
  ! imposed by an Order'th degree difference operator, which has binomial
  ! coefficients with alternating sign.
  !
  ! The Orders argument is the index in the tree of the orders field.
  ! The Quants argument is the index in the tree of the quants field.
  ! The values of the quants field are quantity type names, to which the
  ! corresponding order applies.   The number of orders shall be one (in
  ! which case it applies to all quants), or the same as the number of
  ! quants.  If a column block of the A matrix is of a quantity type that
  ! is not represented in the quants field, no regularization is applied
  ! to that block/
  !
  ! It is necessary that the total number of rows in the first row block of
  ! A be large enough to accomodate the regularization~-- roughly at least
  ! (number of columns of a)~- min(order).  The number of columns of each
  ! block shall be at least one more than the regularization order for that
  ! block.
  !
  ! The Weight argument is a scalar that multiplies the regularization
  ! matrix.

    type(Matrix_T), intent(inout) :: A
    integer, intent(in) :: Orders
    integer, intent(in) :: Quants
    real(r8), intent(in) :: Weight

    integer, dimension(0:maxRegOrd) :: C ! Binomial regularization coefficients
    ! The regularization order is never going to be so large that
    ! the coefficients cannot be represented by integers.
    integer :: I, J                ! Subscripts, Loop inductors
    integer :: IB                  ! Which block is being regularized
    integer :: NB                  ! Number of column blocks of A
    integer :: NCOL                ! Number of columns in a block of A
    integer :: NROW                ! Number of rows in first row block of A
    integer :: NV                  ! Next element in VALUES component
    integer :: Ord                 ! Order for the current block
    integer :: Rows                ! Next row to use
    integer :: S                   ! Sign of regularization coefficient.
    integer :: Type                ! Type of value returned by EXPR
    integer :: Units(2)            ! Units of value returned by EXPR
    double precision :: Value(2)   ! Value returned by EXPR

    error = 0
    nb = a%col%nb
    if ( a%col%extra ) nb = nb - 1
    if ( nsons(orders) /= 2 ) then
      if ( quants == 0 ) call announceError ( regQuantsReq, orders )
      if ( nsons(orders) /= nsons(quants) ) &
        & call announceError ( fieldSizes, orders )
    end if

    if ( error == 0 ) then
      nrow = a%row%nelts(1)
      rows = 1
o:    do ib = 1, nb
        ord = 0
        if ( nsons(orders) == 2 ) then
          call expr ( subtree(2,orders), units, value, type )
          ord = value(1)
        else
          do i = 2, nsons(quants)
            if ( decoration(subtree(i,quants)) == &
              & a%col%vec%quantities(a%col%quant(i))%template%quantityType ) then
              call expr ( subtree(i,orders), units, value, type )
              ord = value(1)
              exit
            end if
          end do
        end if
        ncol = a%col%nelts(ib)
        if ( rows + ncol - ord - 1 > a%row%nelts(1) ) &
          & call announceError ( tooFewRows, orders )
        if ( ord == 0 ) then
          call createBlock ( a%block(1,ib), nrow, ncol,  m_absent )
        else
          if ( ncol < ord+1 ) call announceError ( tooFewCols, orders )
          if ( ord > maxRegOrd ) call announceError ( orderTooBig, orders )
          if ( error > 0 ) exit o
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
          ! has the binomial coefficients in reverse order (which wouldn't
          ! matter if the signs didn't alternate).  Except the first and last
          ! ord columns have 1, 2, ..., ord and ord, ord-1, ..., 1 elements.

          nv = 1
          do i = 1, ord
            a%block(1,ib)%r1(i) = rows
            a%block(1,ib)%r2(i) = nv
            a%block(1,ib)%values(nv:nv+i-1,1) = weight * c(ord-i+1:ord)
            nv = nv + i
          end do
          do i = ord+1, ncol-ord
            a%block(1,ib)%r1(i) = rows
            a%block(1,ib)%r2(i) = nv
            a%block(1,ib)%values(nv:nv+ord,1) = weight * c
            nv = nv + ord + 1
            rows = rows + 1
          end do
          j = ord-1
          do i = ncol-ord+1, ncol
            a%block(1,ib)%r1(i) = a%block(1,ib)%r1(i-1)+1
            a%block(1,ib)%r2(i) = nv
            a%block(1,ib)%values(nv:nv+j,1) = weight * c(0:j)
            nv = nv + j + 1
          end do
        end if
        rows = rows + ncol - ord
      end do o
    end if ! error == 0
    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Regularization failed." )
    if ( a%col%extra ) &
      & call createBlock ( a%block(1,a%col%nb), nrow, 1,  m_absent )
  end subroutine Regularize

  subroutine AnnounceError ( code, where )
    integer, intent(in) :: Code    ! The message number
    integer, intent(in) :: Where   ! Where in the tree

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ' RetrievalModule complained: ' )
    select case ( code )
    case ( fieldSizes )       ! size(regOrders) /= size(regQuants)
      call output ( "Number of values of regOrders shall be 1 or the same as " )
      call output ( "for regQuants.", advance="yes" )
    case ( orderTooBig )
      call output ( "Regularization order exceeds " )
      call output ( maxRegOrd, advance="yes" )
    case ( regQuantsReq )     ! RegQuants required if size(regOrders) /= 1
      call output ( "The regQuants field is required if more than one order " )
      call output ( "is specified.", advance="yes" )
    case ( tooFewCols )
      call output ( "Not enough columns for specified regularization order.", &
        & advance="yes" )
    case ( tooFewRows )
      call output ( "Not enough rows in the matrix to do regularization.", &
        & advance="yes" )
    case ( unitless )         ! regOrders must be unitless
      call output ( "The orders shall be unitless.", advance="yes" )
    end select

  end subroutine AnnounceError

end module Regularization

! $Log$
! Revision 2.5  2001/06/26 20:13:24  vsnyder
! Fixed a blunder -- call call announceError
!
! Revision 2.4  2001/06/26 19:01:00  vsnyder
! Specify regularization orders according to quantities
!
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
