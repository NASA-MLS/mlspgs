! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module HessianModule_1          ! High-level Hessians in the MLS PGS suite
!=============================================================================

! This module provides the elementary Hessian type.  Blocks of this
! type are used to compose block Hessians.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Deallocate
  use HessianModule_0, only: ClearBlock, CreateBlock, DestroyBlock, HessianElement_T, H_Absent, &
    & H_Sparse, H_Full, OptimizeBlock
  use MLSKinds, only: RM
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use MatrixModule_1, only: DefineRCInfo, DestroyRCInfo, NullifyRCInfo, RC_Info
  use VectorsModule, only: Vector_T

  implicit NONE
  private

  type :: Hessian_T
    integer :: Name = 0  ! Sub-rosa index of hHessian name, if any, else zero
    integer :: Where = 0 ! Source_ref for creation if by L2CF
    type(RC_Info) :: Col, Row  ! Column and row info
    type(HessianElement_T), dimension(:,:,:), pointer :: BLOCK => NULL()
  end type Hessian_T
 
  public :: AddHessianToDatabase, CreateEmptyHessian, DestroyHessian
  public :: DestroyHessianDatabase, Hessian_T
  public :: InsertHessianPlane, Multiply, NullifyHessian, StreamlineHessian

  interface CreateBlock
    module procedure CreateHessianBlock_1
  end interface

  interface InsertHessianPlane
    module procedure InsertHessianPlane_1
  end interface

  interface Multiply
    module procedure Hessian_Vector_Vector_Multiply
  end interface

  interface NullifyHessian
    module procedure NullifyHessian_1
  end interface


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------  AddHessianToDatabase  -----
  integer function AddHessianToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector to a database of such vectors, 
  ! creating the database if necessary.

    ! Dummy arguments
    type (Hessian_T), dimension(:), pointer :: DATABASE
    type (Hessian_T), intent(in) ::            ITEM

    ! Local variables
    type (Hessian_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddHessianToDatabase = newSize
  end function AddHessianToDatabase

  ! ------------------------------ CreateEmptyHessian ---

  type (Hessian_T) function CreateEmptyHessian ( Name, Row, Col, &
    & Row_Quan_First, Col_Quan_First, Text, where ) result ( H )
    use Symbol_Types, only: T_identifier
    use Symbol_Table, only: Enter_Terminal

    integer, intent(in) :: Name         ! Sub-rosa index of its name, or zero
    type (Vector_T), intent(in) :: Row   ! Vector used to define the row
    !                                     space of the matrix.
    type (Vector_T), intent(in) :: Col   ! Vector used to define the column
    !                                     space of the matrix.
    logical, intent(in), optional :: Row_Quan_First    ! True (default false)
      ! means the quantity is the major order of the rows of blocks and the
      ! instance is the minor order.
    logical, intent(in), optional :: Col_Quan_First    ! True (default false)
      ! means the quantity is the major order of the columns of blocks and the
      ! instance is the minor order.
    character(len=*), intent(in), optional :: Text     ! A name to use
      ! instead of "Name."
    integer, intent(in), optional :: Where             ! source_ref

    integer :: I, J, K                  ! Loop counters
    integer :: STATUS                   ! Flag

    call destroyHessian ( H )  ! Avoid a memory leak if it isn't freshly minted
    h%name = name
    if ( present(text) ) h%name = enter_terminal ( text, t_identifier )
    call defineRCInfo ( h%row, row, row_quan_first )
    call defineRCInfo ( h%col, col, col_quan_first )
    if ( present(where) ) h%where = where

    allocate ( h%block ( h%row%nb, h%col%nb, h%col%nb ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "H%Block in CreateEmptyHessian" )

    do i = 1, h%row%nb ! Now create absent blocks with the correct sizes
      do j = 1, h%col%nb
        do k = 1, h%col%nb
          call createBlock ( h, i, j, k, h_absent )
        end do
      end do
    end do
  end function CreateEmptyHessian

  ! ----------------------------------- CreateHessianBlock_1 ----
  subroutine CreateHessianBlock_1 ( H, RowNum, ColNum1, ColNum2, Kind, InitTuples )
  ! Create the hessian block Z%Block(RowNum,ColNum), which sprang into
  ! existence with kind M_Absent.  Create it with the specified Kind.
  ! See HessianModule_0 for a list of the kinds.  If the Kind is
  ! M_Sparse the initial number of tuples is required
    type (Hessian_T), intent(inout) :: H       ! The matrix having the block
    integer, intent(in) :: RowNum, ColNum1, ColNum2    ! Row and column of the block
    integer, intent(in) :: Kind         ! Kind of block, see MatrixModule_0
    integer, intent(in), optional :: InitTuples  ! Number of nonzeros
    call createBlock ( h%block ( rowNum, colNum1, colNum2 ), &
      & h%row%nelts ( rowNum ), &
      & h%col%nelts ( colNum1 ), &
      & h%col%nelts ( colNum2 ), &
      & kind, &
      & initTuples )
  end subroutine CreateHessianBlock_1

  ! --------------------------------------  DestroyHessian  -----
  subroutine DestroyHessian ( hessian )
    ! This subroutine destroys a hessian
    
    ! Dummy argument
    type (Hessian_T), intent (inout) :: hessian
    integer :: I, J, K
    
    integer :: status
    hessian%name = 0
    hessian%where = 0
    do i = 1, hessian%row%nb
      do j = 1, hessian%col%nb
        do k = 1, hessian%col%nb
          call DestroyBlock ( hessian%block(i,j,k) )
        end do
      end do
    end do
    call destroyRCInfo ( hessian%row )
    call destroyRCInfo ( hessian%col )
    if ( associated(hessian%block) ) then
      deallocate ( hessian%block, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // "HESSIAN%H in DestroyHessian" )
    end if
  end subroutine DestroyHessian

  ! --------------------------------------  DestroyHessianDatabase  -----
  subroutine DestroyHessianDatabase ( database )

  ! This subroutine destroys a vector database

    ! Dummy argument
    type (Hessian_T),  dimension(:), pointer :: database

    ! Local variables
    integer :: hessianIndex, Status

    if ( associated(database) ) then
      do hessianIndex = 1, SIZE(database)
        call DestroyHessian ( database(hessianIndex) )
      end do
      deallocate ( database, stat=status )
      call test_deallocate ( status, moduleName, 'database', 0 )
    end if
  end subroutine DestroyHessianDatabase

  ! ------------------------------- Hessian_Vector_Vector_Multiply -----
  subroutine Hessian_Vector_Vector_Multiply ( H, V, P, Update )
  !{ Multiply a Hessian {\tt H} by {\tt V} twice, with a factor of $\frac12$,
  !  giving {\tt P}: $P^k = \frac12 H^k_{ij} V^i V^j$ or $P^k = P^k + \frac12
  !  H^k_{ij} V^i V^j$, depending upon {\tt Update}.  This is the
  !  second-order term of a Taylor series.  {\tt P} is initially set to zero
  !  unless {\tt Update} is present and true.

    use HessianModule_0, only: Multiply
    use VectorsModule, only: Vector_T

    type(Hessian_T), intent(in) :: H
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: P
    logical, intent(in), optional :: Update

    integer :: I, J, K, L
    integer :: IQ, JQ, KQ               ! Quantity indices
    integer :: II, JI, KI               ! Instance indices
    logical :: MYUPDATE                 ! Copy of update

    myUpdate = .false.
    if ( present ( update ) ) myUpdate = update
    ! Error checking on campatibility of H%row with P and H%col with V here.

    do i = 1, h%row%nb
      iq = h%col%quant ( i )
      ii = h%col%inst ( i )
      do j = 1, h%col%nb
        jq = h%col%quant ( j )
        ji = h%col%inst ( j )
        do k = 1, h%col%nb
          kq = h%col%quant ( k )
          ki = h%col%inst ( k )
          do l = 1, h%col%vec%quantities(kq)%template%instanceLen
            call multiply ( h%block(i,j,k), &
              & v%quantities(iq)%values(:,ii), &
              & v%quantities(jq)%values(:,ji), &
              & p%quantities(kq)%values(:,ki), &
              & myUpdate )
          end do
        end do
      end do
    end do

  end subroutine Hessian_Vector_Vector_Multiply

  ! ---------------------------------- InsertHessianPlane_1 ----
  subroutine InsertHessianPlane_1 ( H, M, B, EL, MIRROR )
    ! Insert matrix M as a plane (block B, element EL) of the Hessian H
    ! If set, populate the transpose set also
    use HessianModule_0, only: InsertHessianPlane
    use MatrixModule_1, only: Matrix_T
    ! Dummy arguments
    type(Hessian_T), intent(inout) :: H
    type(Matrix_T), intent(in) :: M
    integer, intent(in) :: B            ! One of the column blocks
    integer, intent(in) :: EL           ! The element of that block
    logical, intent(in), optional :: MIRROR
    ! Local variables
    integer :: CB                       ! Column block
    integer :: RB                       ! Row block
    logical :: MYMIRROR                 ! Copy of mirror

    myMirror = .false.
    if ( present ( mirror ) ) myMirror = mirror

    ! Check that the jacobian and hessian match
    if ( h%col%vec%template%name /= m%col%vec%template%name &
      & .or. h%row%vec%template%name /= m%row%vec%template%name &
      & .or. (h%col%instFirst .neqv. m%col%instFirst) &
      & .or. (h%row%instFirst .neqv. m%row%instFirst) )&
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Hessian and Matrix incompatible in InsertHessianPlane' )

    ! Loop over the rows
    do rb = 1, H%row%nb
      do cb = 1, H%col%nb
        call InsertHessianPlane ( H%block(rb,cb,b), M%block(rb,cb), el )
        if ( myMirror ) &
          & call InsertHessianPlane ( H%block(rb,b,cb), M%block(rb,cb), &
          & el, mirroring=.true. )
      end do
    end do
    
  end subroutine InsertHessianPlane_1

  ! -------------------------------------- NullifyHessian ------
  subroutine NullifyHessian_1 ( H )
    ! Given a matrix, nullify all the pointers associated with it
    type(Hessian_T), intent(inout) :: H

    ! Executable code
    h%name = 0
    call nullifyRCInfo ( h%col )
    call nullifyRCInfo ( h%row )
    nullify ( h%block )
  end subroutine NullifyHessian_1

  ! -------------------------------------- StreamlineHessian ------
  subroutine StreamlineHessian ( H, scaleHeight, geodAngle )
    use Intrinsic, only: L_ZETA
    use MLSKinds, only: R8
    use QuantityTemplates, only: QuantityTemplate_T
    ! Given a Hessian, trim off the elements that are further away than indicated in
    ! scale height (for zeta coordinates) and geodAngle
    type (Hessian_T), intent(inout) :: H
    type (HessianElement_T), pointer :: HB
    real(r8), intent(in) :: scaleHeight
    real(r8), intent(in) :: geodAngle

    type (QuantityTemplate_T), pointer :: Q1, Q2
    integer :: P1, P2, S1, S2           ! Profile and surface for column quantities
    integer :: I, J, K                  ! Loop counters
    integer :: T, II, JJ, KK            ! More loop counters and indices
    logical :: DROPBLOCK

    do i = 1, h%row%nb
      do j = 1, h%col%nb
        p1 = h%col%inst(j)
        q1 => h%col%vec%quantities ( h%col%quant(j) ) % template
        do k = 1, h%col%nb
          p2 = h%col%inst(k)
          q2 => h%col%vec%quantities ( h%col%quant(k) ) % template

          hb => h%block ( i, j, k )
          
          ! Decide whether to drop this on horizontal grounds
          dropBlock = .false.
          if ( q1%stacked .and. q2%stacked .and. ( geodAngle > 0 ) ) &
            & dropBlock = abs ( q1%phi(1,p1) - q2%phi(1,p2) ) > geodAngle
          if ( dropBlock ) then
            call ClearBlock ( hb )
          else
            ! Now clear elements on vertical grounds, but only for simple quantities
            if ( .not. q1%coherent .or. .not. q2%coherent ) continue
            if ( q1%verticalCoordinate /= l_zeta .or. &
              &  q2%verticalCoordinate /= l_zeta ) continue

            select case ( hb%kind )
            case ( h_absent )
            case ( h_sparse )
              do t = 1, hb%tuplesFilled
                s1 = ( hb%tuples(t)%j-1 ) / q1%noChans + 1
                s2 = ( hb%tuples(t)%k-1 ) / q1%noChans + 1
                if ( abs ( q1%surfs(s1,1) - q2%surfs(s2,1) ) > scaleHeight ) &
                  & hb%tuples(t)%h = 0
              end do
            case ( h_full )
              do kk = 1, hb%nCols2
                s2 = ( kk-1 ) / q2%noChans + 1
                do jj = 1, hb%nCols1
                  s1 = ( jj-1 ) / q1%noChans + 1
                  if ( abs ( q1%surfs(s1,1) - q2%surfs(s2,1) ) > scaleHeight ) &
                    & hb%values ( :, jj, kk ) = 0
                end do
              end do
            end select

            ! Now optimize this block
            call OptimizeBlock ( hb )
          end if
        end do
      end do
    end do

  end subroutine StreamlineHessian

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HessianModule_1

! $Log$
! Revision 2.2  2010/02/25 21:14:33  pwagner
! Replaced non-standard component separator '.' with '%'
!
! Revision 2.1  2010/02/25 18:14:02  pwagner
! First commit
!
