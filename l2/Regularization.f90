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
  subroutine Regularize ( A, Orders, Quants, Weights, Rows )

  !{Apply Tikhonov regularization conditions of the form $\Delta^k \delta
  ! \bf{x} \simeq 0$ to the blocks of A, where $\Delta^k$ is the central
  ! difference operator of order $k$.  $\Delta^k$ has binomial
  ! coefficients with alternating sign.
  !
  ! $k$ is given by the {\tt Orders} argument, which is the index in the tree
  ! of the {\tt regOrders} field of the {\tt retrieve} specification. The
  ! {\tt Weights} argument is the index in the tree of the {\tt regWeights}
  ! field.  It gives a weight for the regularization. The {\tt Quants}
  ! argument is the index in the tree of the {\tt regQuants} field. The
  ! values of the {\tt regQuants} field are quantity type names.   The number
  ! of values of {\tt regOrders} and {\tt regWeights} shall be one (in which
  ! case they apply to all quants), or the same as the number of {\tt
  ! regQuants}, in which case the corresponding order or weight applies to
  ! the corresponding quantity.  If a column block of the A matrix is of a
  ! quantity type that is not represented in the {\tt regQuants} field, no
  ! regularization is applied to that block.
  !
  ! It is necessary that the total number of rows in the first row block of
  ! A be large enough to accomodate the regularization~-- roughly at least
  ! (number of columns of a)~- min(order)).  The number of columns of each
  ! block shall be at least one more than the regularization order for that
  ! block.

    type(Matrix_T), intent(inout) :: A
    integer, intent(in) :: Orders
    integer, intent(in) :: Quants
    integer, intent(in) :: Weights
    integer, intent(out) :: Rows   ! Last row used; ultimately, number of rows

    integer, dimension(0:maxRegOrd) :: C ! Binomial regularization coefficients
    ! The regularization order is never going to be so large that
    ! the coefficients cannot be represented by integers.
    integer :: I, J                ! Subscripts, Loop inductors
    integer :: IB                  ! Which block is being regularized
    integer :: MaxRow              ! Maximum row to be filled = ncol - ord
    integer :: NB                  ! Number of column blocks of A
    integer :: NCOL                ! Number of columns in a block of A
    integer :: NROW                ! Number of rows in first row block of A
    integer :: NV                  ! Next element in VALUES component
    integer :: Ord                 ! Order for the current block
    integer :: S                   ! Sign of regularization coefficient.
    integer :: Type                ! Type of value returned by EXPR
    integer :: Units(2)            ! Units of value returned by EXPR
    double precision :: Value(2)   ! Value returned by EXPR
    real(r8) :: Wt                 ! The weight for the current block

    error = 0
    nb = a%col%nb
    if ( nsons(orders) /= 2 ) then
      if ( quants == 0 ) then
        call announceError ( regQuantsReq, orders )
      else if ( nsons(orders) /= nsons(quants) ) then
        call announceError ( fieldSizes, orders )
      end if
    end if
    if ( nsons(weights) /= 2 ) then
      if ( quants == 0 ) then
        call announceError ( regQuantsReq, weights )
      else if ( nsons(weights) /= nsons(quants) ) then
        call announceError ( fieldSizes, weights )
      end if
    end if
    rows = 0

    if ( error == 0 ) then
      nrow = a%row%nelts(1)
      do ib = 1, nb             ! Loop over matrix blocks = quantities
        ord = 0
        if ( quants == 0 ) then ! only one order and weight allowed
          call expr ( subtree(2,orders), units, value, type )
          ord = value(1)
          call expr ( subtree(2,weights), units, value, type )
          wt = value(1)
        else
          do i = 2, nsons(quants)
            if ( decoration(subtree(i,quants)) == &
              & a%col%vec%quantities(a%col%quant(ib))%template%quantityType ) then
              j = min(i,nsons(orders))
              call expr ( subtree(j,orders), units, value, type )
              ord = value(1)
              j = min(i,nsons(weights))
              call expr ( subtree(j,weights), units, value, type )
              wt = value(1)
              exit
            end if
          end do
        end if
        ncol = a%col%nelts(ib)
        if ( rows + ncol - ord - 1 > a%row%nelts(1) ) &
          & call announceError ( tooFewRows, orders )
        if ( ncol < ord+1 ) call announceError ( tooFewCols, orders )
        if ( ord > maxRegOrd ) call announceError ( orderTooBig, orders )
        if ( ord == 0 .or. wt <= 0.0_r8 .or. error /= 0 ) then
          call createBlock ( a%block(1,ib), nrow, ncol,  m_absent )
        else
          call createBlock ( a%block(1,ib), nrow, ncol, m_banded, &
            & (ncol-ord)*(ord+1) )
          rows = rows + 1

          !{ Calculate binomial coefficients with alternating sign,
          ! $(-1)^i C_i^n = (-1)^i \frac{n!}{i! (n-i)!}$
          !  by the recursion
          !  $C_0^n = 1\text{, } C_i^n = -(n-i+1) C_{i-1}^n / i$.
          ! Notice that $C_i^n = C_{n-i}^n$, so we only need
          ! to go halfway through the array.
          s = 1 - 2*mod(ord,2) ! +1 for even order, -1 for odd order
          c(0) = s
          c(ord) = 1
          do i = 1, ord / 2
            c(i) = ( -(ord-i+1) * c(i-1) ) / i
            c(ord-i) = s * c(i)
          end do

          ! Each row has the binomial coefficients.  Therefore, each column
          ! has the binomial coefficients in reverse order (which wouldn't
          ! matter if the signs didn't alternate).  Except the first and
          ! last ord columns have 1, 2, ..., ord + 1 (or less if ncol-ord <
          ! ord) and ord, ord-1, ..., 1 elements.

          maxRow = ncol - ord
          ! Fill in coefficients from the end of Ord (but no more than
          ! maxRow of them)
          do i = 1, ord+1
            nv = a%block(1,ib)%r2(i-1) + 1
            j = min(i,maxRow) ! Number of coefficients
            a%block(1,ib)%r1(i) = rows
            a%block(1,ib)%r2(i) = nv+j-1
            a%block(1,ib)%values(nv:nv+j-1,1) = c(ord-i+1:ord-i+j)
          end do
          ! Fill in coefficients from all of Ord
          do i = ord+2, maxRow
            nv = a%block(1,ib)%r2(i-1) + 1
            rows = rows + 1
            a%block(1,ib)%r1(i) = rows
            a%block(1,ib)%r2(i) = nv+ord
            a%block(1,ib)%values(nv:nv+ord,1) = c(0:ord)
          end do
          ! Fill in coefficients from the beginning of Ord
          j = min(maxrow-1,ord) - 1 ! Index of last coefficient
          do i = max(ord+2,maxRow+1), ncol
            nv = a%block(1,ib)%r2(i-1) + 1
            rows = rows + 1
            a%block(1,ib)%r1(i) = rows
            a%block(1,ib)%r2(i) = nv+j
            a%block(1,ib)%values(nv:nv+j,1) = c(0:j)
            j = j - 1
          end do
        end if
      end do  
    end if ! error == 0
    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Regularization failed." )
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
      call output ( "Number of values of regOrders or regWeights shall be 1 " )
      call output ( "or the same as for regQuants.", advance="yes" )
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
! Revision 2.12  2002/05/22 19:15:06  vsnyder
! Output the number of rows used for regularization conditions
!
! Revision 2.11  2002/05/11 01:10:16  vsnyder
! Big revision... Maybe it's right this time.
!
! Revision 2.10  2002/05/07 01:03:36  vsnyder
! Allow different weight for each quantity, or one weight for the whole
! shebang.  Don't regularize a block if its weight is <= zero.
!
! Revision 2.9  2002/02/01 00:48:26  vsnyder
! Get rid of 'extra' field of RC_Info
!
! Revision 2.8  2001/10/09 20:36:04  vsnyder
! Repair calculation of banded representation
!
! Revision 2.7  2001/06/28 20:42:42  vsnyder
! Update comments
!
! Revision 2.6  2001/06/26 20:16:41  vsnyder
! Don't look at nsons(0)
!
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
