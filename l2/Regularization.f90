! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Regularization
!=============================================================================

! Apply Tikhonov regularization to the matrix of the least-squares
! problem for retrieval.

  private

  public :: Regularize

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  Regularize  -----
  subroutine Regularize ( A, Orders, Quants, Weights, WeightVec, Rows )

  !{Apply Tikhonov regularization conditions of the form $\Delta^k \delta
  ! \bf{x} \simeq 0$ to the blocks of A, where $\Delta^k$ is the central
  ! difference operator of order $k$.  $\Delta^k$ has binomial
  ! coefficients with alternating sign.  The operator is normalized so that
  ! the sum of the coefficients is one.
  !
  ! The order $k$ is given by the {\tt Orders} argument, which is the index in
  ! the tree of the {\tt regOrders} field of the {\tt retrieve} specification.
  ! If one value is supplied for {\tt regOrders} it applies to all quantities.
  ! If several values are supplied, the {\tt regWeights} (if present) and {\tt
  ! regQuants} fields must have the same number of values, and corresponding
  ! elements of the {\tt regOrders} and {\tt regWeights} fields apply to the
  ! quantity specified by the corresponding element of the {\tt regQuants}
  ! field.
  !
  ! The {\tt Quants} argument is the index in the tree of the {\tt regQuants}
  ! field. The values of the {\tt regQuants} field are quantity type names.
  ! If no {\tt regQuants} field is given, all quantities are regularized, and
  ! only one value may be given for each of the {\tt regOrders} and {\tt
  ! regWeights} fields.
  !
  ! The {\tt Weights} argument is the index in the tree of the {\tt regWeights}
  ! field.  It gives a weight for the regularization.  If one value is given,
  ! it is used for every quantity.  If several are given, the number given must
  ! be the same as the number of quantities, and corresponding weights are
  ! applied to corresponding quantities.  The weights can also be given by {\tt
  ! WeightVec}, which if present must have the same template as the column
  ! template for A.  For each quantity, the weight vector is averaged to the
  ! number of rows of regularization (i.e., (number of columns)~- (order of
  ! regularization operator)) using the absolute value of the regularization
  ! operator.  If both {\tt regWeights} and {\tt regWeightVec} are provided,
  ! their product is used.
  !
  ! It is necessary that the number of rows in the row block of A that has
  ! the most rows be large enough to accomodate the regularization~--
  ! roughly at least (number of columns of a)~- min(order)).  If the
  ! regularization order for a block is less than one less than the number of
  ! columns of that block, the regularization order is set to one less than
  ! the number of columns for that block, and a warning message is emitted.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_M, only: EXPR
    use MatrixModule_0, only: CreateBlock, M_Absent, M_Banded, &
      & MatrixElement_T
    use MatrixModule_1, only: Matrix_T
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use Output_M, only: Output
    use Toggles, only: Switches
    use Tree, only: DECORATION, NSONS, SOURCE_REF, SUBTREE
    use VectorsModule, only: M_Tikhonov, Vector_T

    type(Matrix_T), intent(inout) :: A
    integer, intent(in) :: Orders
    integer, intent(in) :: Quants
    integer, intent(in) :: Weights
    type(vector_T), pointer :: WeightVec
    integer, intent(out) :: Rows   ! Last row used; ultimately, number of rows

    integer, parameter :: MaxRegOrd = 56 ! Maximum regularization
    ! order.  26!/(13!)**2 < 1/EPSILON(0.0E0) < 27!/(13!)**2 for 24-bit fraction.
    !         56!/(28!)**2 < 1/EPSILON(0.0D0) < 57!/(28!**2) for 53-bit fraction.

    ! Error message codes
    integer, parameter :: FieldSizes = 1   ! size(regOrders) /= size(regQuants)
    integer, parameter :: OrderTooBig = FieldSizes + 1
    integer, parameter :: RegQuantsReq = OrderTooBig + 1 ! RegQuants required
    integer, parameter :: RegTemplate = RegQuantsReq + 1 ! Weight /= column of J
    integer, parameter :: TooFewRows = RegTemplate + 1   ! Won't fit
    integer, parameter :: Unitless = TooFewRows + 1      ! Orders must be unitless

    integer :: Error               ! non-zero if an error occurs

    integer :: C1, C2              ! Column boundaries, esp. if a mask is used.
    integer :: I, J                ! Subscripts, Loop inductors
    integer :: IB                  ! Which block is being regularized
    integer :: MaxCols             ! Columns in block with with max cols.
    integer :: NB                  ! Number of column blocks of A
    integer :: NCOL                ! Number of columns in a block of A
    integer :: NI, NQ              ! Indices for instance and quantity
    integer :: NROW                ! Number of rows in first row block of A
    integer :: Ord                 ! Order for the current block
    integer :: Type                ! Type of value returned by EXPR
    integer :: Units(2)            ! Units of value returned by EXPR
    double precision :: Value(2)   ! Value returned by EXPR
    logical :: Warn                ! Send warning message to MLSMessage
    real(r8) :: Wt                 ! The weight for the current block
    real(r8), pointer :: WtVec(:)  ! In case a quantity has a vector weight

    error = 0
    nb = a%col%nb
    rows = 0
    warn = .false.
    if ( nsons(orders) /= 2 ) then
      if ( quants == 0 ) then
        call announceError ( regQuantsReq, orders )
      else if ( nsons(orders) /= nsons(quants) ) then
        call announceError ( fieldSizes, orders )
      end if
    end if
    if ( weights == 0 ) then
      wt = 1.0
    else if ( nsons(weights) /= 2 ) then
      if ( quants == 0 ) then
        call announceError ( regQuantsReq, weights )
      else if ( nsons(weights) /= nsons(quants) ) then
        call announceError ( fieldSizes, weights )
      end if
    end if
    if ( associated(weightVec) ) then
      if ( weightVec%template%id /= a%col%vec%template%id ) &
        & call announceError ( regTemplate, orders )
    end if

    if ( error == 0 ) then

      maxCols = 0
      do ib = 1, a%col%nb ! Find number of columns in widest block
        if ( quants == 0 ) then ! Doing all of the blocks
          maxCols = max(maxCols, a%col%nelts(ib))
        else
          do i = 2, nsons(quants)
            if ( decoration(subtree(i,quants)) == &
              & a%col%vec%quantities(a%col%quant(ib))%template%quantityType ) then
                maxCols = max(maxCols, a%col%nelts(ib))
                exit
            end if
          end do
        end if
      end do
      nullify ( wtVec )
      call allocate_test ( wtVec, maxCols, "Weight vector", moduleName )
        if ( index(switches,'reg') /= 0 ) then
          call output ( 'Allocated ' )
          call output ( maxCols )
          call output ( '-element vector for weights.', advance='yes' )
        end if

      do ib = 1, nb             ! Loop over matrix blocks = quantities
        ord = 0
        if ( quants == 0 ) then ! only one order and weight allowed
          call expr ( subtree(2,orders), units, value, type )
          ord = value(1)
          if ( weights /= 0 ) then
            call expr ( subtree(2,weights), units, value, type )
            wt = value(1)
          end if
        else
          do i = 2, nsons(quants)
            if ( decoration(subtree(i,quants)) == &
              & a%col%vec%quantities(a%col%quant(ib))%template%quantityType ) then
              j = min(i,nsons(orders))
              call expr ( subtree(j,orders), units, value, type )
              ord = value(1)
              if ( weights /= 0 ) then
                j = min(i,nsons(weights))
                call expr ( subtree(j,weights), units, value, type )
                wt = value(1)
              end if
              exit
            end if
          end do
        end if
        ncol = a%col%nelts(ib)
        nrow = a%row%nelts(ib)
        if ( ord > ncol-1 ) then
          warn = .true.
          ord = ncol - 1
        end if
        if ( ncol - ord > a%row%nelts(ib) ) then
          print *, 'ncol, ord, nb, a%row%nelts(ib) =', ncol, ord, nb, a%row%nelts(ib)
          call announceError ( tooFewRows, orders )
        end if
        if ( ord > maxRegOrd ) call announceError ( orderTooBig, orders )
        if ( error /= 0 ) warn = .true.
        if ( ord /= 0 .and. wt > 0.0_r8 .and. error == 0 ) then
          call createBlock ( a%block(ib,ib), nrow, ncol, m_banded, &
            & (ncol-ord)*(ord+1) )
          a%block(ib,ib)%r1 = 0         ! in case there's a mask
          a%block(ib,ib)%r2 = 0
          a%block(ib,ib)%values = 0.0_r8

          rows = rows + 1

          ni = a%col%inst(ib)
          nq = a%col%quant(ib)
          if ( .not. associated(a%col%vec%quantities(nq)%mask) ) then
            c1 = 1
            c2 = ncol
            if ( associated(weightVec) ) then
              wtVec(1:weightVec%quantities(nq)%template%instanceLen) = &
                & weightVec%quantities(nq)%values(:,ni)
              call fillBlock ( a%block(ib,ib), ord, rows, 1, ncol, wt, wtVec )
            else
              call fillBlock ( a%block(ib,ib), ord, rows, 1, ncol, wt )
            end if
          else
            ! Scan for blocks of consecutive zero values of M_Tikhonov bits.
            c2 = 1
o:          do while ( c2 <= a%block(ib,ib)%ncols )
              c1 = c2
              do while ( iand(ichar(a%col%vec%quantities(nq)%mask(c1,ni)),M_Tikhonov) &
                  & /= 0 )
                if ( c1 >= a%block(ib,ib)%ncols ) exit o
                c1 = c1 + 1
              end do
              c2 = c1 + 1
              do while ( iand(ichar(a%col%vec%quantities(nq)%mask(c2,ni)),M_Tikhonov) &
                  & == 0 )
                c2 = c2 + 1
                if ( c2 > a%block(ib,ib)%ncols ) exit
              end do
              if ( associated(weightVec) ) then
                wtVec = 0.0_r8
                wtVec ( c1 : c2-1 ) = weightVec%quantities(nq)%values(c1:c2-1,ni)
                call fillBlock ( a%block(ib,ib), ord, rows, c1, c2-1, wt, wtVec )
              else
                call fillBlock ( a%block(ib,ib), ord, rows, c1, c2-1, wt )
              end if
            end do o
          end if

        end if
        error = 0

      end do
      call deallocate_test ( wtVec, "Weight vector", moduleName )
      if ( warn ) call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "Some blocks not regularized, or at lower order than requested" )
    else
      call MLSMessage ( MLSMSG_Error, moduleName, "Regularization failed." )
    end if ! error == 0

  contains

    subroutine AnnounceError ( code, where )
      use Lexer_Core, only: Print_Source

      integer, intent(in) :: Code    ! The message number
      integer, intent(in) :: Where   ! Where in the tree

      error = max(error,1)
      call output ( '***** At or near ' )
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
      case ( regTemplate )
        call output ( "The template for the regularization weights vector is " )
        call output ( "not the same as for the columns of the Jacobian matrix.", &
          & advance='yes' )
      case ( tooFewRows )
        call output ( "Not enough rows in the matrix to do regularization.", &
          & advance="yes" )
      case ( unitless )         ! regOrders must be unitless
        call output ( "The orders shall be unitless.", advance="yes" )
      end select

    end subroutine AnnounceError

    subroutine FillBlock ( B, Ord, Rows, C1, C2, Wt, WtVec )
    ! Fill the block B with regularization coefficients starting at row Rows+1
    ! and columns C1 to C2.  Update Rows to the last row filled.
      type(matrixElement_T), intent(inout) :: B   ! Block to fill
      integer, intent(in) :: Ord        ! Order of regularization operator
      integer, intent(inout) :: Rows    ! Last row used
      integer, intent(in) :: C1, C2     ! Columns to fill
      real(r8), intent(in) :: Wt        ! Scalar weight for coefficients
      real(r8), intent(inout), dimension(:), optional :: WtVec ! Weights vector

      real(r8), dimension(0:maxRegOrd) :: C ! Binomial regularization
      ! coefficients in reverse order.
      integer :: I, J                   ! Subscripts, Loop inductors
      integer :: K                      ! Column being filled
      integer :: MaxRow                 ! Maximum row to be filled = ncol - ord
      integer :: NCOL                   ! Number of columns -- C2 - C1 + 1
      integer :: NV                     ! Next element in VALUES component
      integer :: S                      ! Sign of regularization coefficient.

      ncol = c2 - c1 + 1

      !{ Calculate binomial coefficients with alternating sign,
      ! $(-1)^i C_i^n = (-1)^i \frac{n!}{i! (n-i)!}$
      !  by the recursion
      !  $C_0^n = 1\text{, } C_i^n = -(n-i+1) C_{i-1}^n / i$.
      ! Notice that $C_i^n = -1^n C_{n-i}^n$, so we only need
      ! to go halfway through the array.
      s = 1 - 2*mod(ord,2) ! +1 for even order, -1 for odd order
      c(0) = s
      c(ord) = 1
      do i = 1, ord / 2
        c(i) = ( -(ord-i+1) * c(i-1) ) / i
        c(ord-i) = s * c(i)
      end do
      c(0:ord) = wt * c(0:ord) / 2.0_r8 ** ord

      ! If there is a weight vector, it's the same length as the number
      ! of columns.  But we want to weight the rows.  So construct a new
      ! weight vector that is the same length as the number of rows, i.e.,
      ! (number of columns) - (order), by averaging using the absolute
      ! value of the coefficients.

      if ( present(wtVec) ) then
        do i = 1, ncol
          j = min(i+ord,ncol)
          wtVec(i) = dot_product(wtVec(i:j),abs(c(0:j-i)))
        end do
      end if

      ! Each row has the binomial coefficients.  Therefore, each column
      ! has the binomial coefficients in reverse order (which wouldn't
      ! matter if the signs didn't alternate).  Except the first and
      ! last ord columns have 1, 2, ..., ord and ord, ord-1, ..., 1
      ! elements.  E.g., for order three, the first four rows look like:

      !  1  -3   3  -1
      !      1  -3   3  -1
      !          1  -3   3  -1
      !              1  -3   3  -1

      ! (assuming there are at least seven columns)

      k = c1
      maxRow = ncol - ord
      ! Fill in coefficients from the end of C(:Ord) (but no more than
      ! maxRow-1 of them)
      do i = 1, ord + 1
        nv = b%r2(k-1) + 1
        j = min(i,maxRow) ! Number of coefficients
        b%r1(k) = i
        b%r2(k) = nv + j - 1
        if ( present(wtVec) ) then
          b%values(nv:nv+j-1,1) = c(ord-i+1:ord-i+j) * wtVec(i:i+j-1)
        else
          b%values(nv:nv+j-1,1) = c(ord-i+1:ord-i+j)
        end if
        k = k + 1
      end do
      ! Fill in coefficients from all of C(:Ord)
      do i = ord+2, maxRow
        nv = b%r2(k-1) + 1
        b%r1(k) = i - ord
        b%r2(k) = nv + ord
        if ( present(wtVec) ) then
          b%values(nv:nv+ord,1) = c(0:ord) * wtVec(i:i+ord)
        else
          b%values(nv:nv+ord,1) = c(0:ord)
        end if
        k = k + 1
      end do
      ! Fill in coefficients from the beginning of C(:Ord) (but no more
      ! than maxRow-1 of them)
      j = min(maxrow-1,ord) - 1 ! Index of last coefficient
      do i = max(ord+2,maxRow+1), ncol
        nv = b%r2(k-1) + 1
        b%r1(k) = i - ord
        b%r2(k) = nv + j
        if ( present(wtVec) ) then
          b%values(nv:nv+j,1) = c(0:j) * wtVec(ncol-j:ncol)
        else
          b%values(nv:nv+j,1) = c(0:j)
        end if
        j = j - 1
        k = k + 1
      end do
      rows = rows + maxRow
    end subroutine FillBlock

  end subroutine Regularize

end module Regularization

! $Log$
! Revision 2.19  2002/08/20 19:54:19  vsnyder
! Re-arrange to use a matrix with the same row and column vector definitions.
! Don't try to put all of the regularization into one block.  Instead, fill
! each diagonal block of the matrix with regularization for that quantity
! and instance.
!
! Revision 2.18  2002/08/16 21:33:02  livesey
! Bug fix in tooFewRows message
!
! Revision 2.17  2002/08/15 22:45:05  livesey
! Lots of array bounds/indexing changes (fixes hopefully).
!
! Revision 2.16  2002/08/08 22:02:05  vsnyder
! Implement mask and weight vector.  Make AnnounceError an internal subroutine
! instead of a module subroutine.  Move USE statements into subroutines.
!
! Revision 2.15  2002/08/03 01:15:36  vsnyder
! Fixed the comment about row sizes
!
! Revision 2.14  2002/07/30 23:20:20  vsnyder
! Use order ncol-1 instead of zero for narrow blocks, Warn if ord < ncol-1.
! Use the weight (duh!). Use the row block with the most rows.
!
! Revision 2.13  2002/07/02 01:45:56  vsnyder
! Regularization.f90
!
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
