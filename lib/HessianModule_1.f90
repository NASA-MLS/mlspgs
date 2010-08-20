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
  use HessianModule_0, only: ClearBlock, CreateBlock, &
    & DestroyBlock, HessianElement_T, H_Absent, &
    & H_Sparse, H_Full, OptimizeBlock
  use LEXER_CORE, only: INIT_LEXER
  use MLSKinds, only: RM
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use MatrixModule_1, only: DefineRCInfo, DestroyRCInfo, NullifyRCInfo, RC_Info
  use Output_M, only: Output, outputNamedValue
  use VectorsModule, only: Vector_T

  implicit NONE
  private

  type :: Hessian_T
    integer :: Name = 0  ! Sub-rosa index of Hessian name, if any, else zero
    integer :: Where = 0 ! Source_ref for creation if by L2CF
    type(RC_Info) :: Col, Row  ! Column and row info
    type(HessianElement_T), dimension(:,:,:), pointer :: BLOCK => NULL()
  end type Hessian_T
 
  public :: AddHessianToDatabase, CreateBlock, CreateEmptyHessian
  public :: DestroyHessian, DestroyHessianDatabase, Diff, Dump, Hessian_T
  public :: InsertHessianPlane, Multiply, NullifyHessian
  public :: OptimizeHessian, StreamlineHessian

  interface CreateBlock
    module procedure CreateHessianBlock_1
  end interface

  interface Diff
    module procedure Diff_Hessians
  end interface

  interface Dump
    module procedure Dump_Hessian, Dump_Hessian_Database
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

  interface OptimizeHessian
    module procedure OptimizeHessian_1
  end interface

  interface StreamlineHessian
    module procedure StreamlineHessian_1
  end interface

  logical, parameter :: DEEBUG = .false.
  logical, parameter :: HIDEBSENTBLOCKS = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------- AddHessianToDatabase -----
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

  ! ------------------------------------------- CreateEmptyHessian -----

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

  ! ----------------------------------------- CreateHessianBlock_1 -----
  subroutine CreateHessianBlock_1 ( H, RowNum, ColNum1, ColNum2, Kind, InitTuples )
  ! Create the hessian block H%Block(RowNum,ColNum), which sprang into
  ! existence with kind M_Absent.  Create it with the specified Kind.
  ! See HessianModule_0 for a list of the kinds.  If the Kind is
  ! M_Sparse the initial number of tuples is required
    type (Hessian_T), intent(inout) :: H ! The matrix having the block
    integer, intent(in) :: RowNum, ColNum1, ColNum2 ! Row and column of the block
    integer, intent(in) :: Kind          ! Kind of block, see HessianModule_0
    integer, intent(in), optional :: InitTuples     ! Number of nonzeros
    call createBlock ( h%block ( rowNum, colNum1, colNum2 ), &
      & h%row%nelts ( rowNum ), &
      & h%col%nelts ( colNum1 ), &
      & h%col%nelts ( colNum2 ), &
      & kind, &
      & initTuples )
  end subroutine CreateHessianBlock_1

  ! ----------------------------------------------- DestroyHessian -----
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

  ! --------------------------------------- DestroyHessianDatabase -----
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

  ! ------------------------------------------------- Diff_Hessians -----
  ! Diff two Hessians, presumed to match somehow
  subroutine Diff_Hessians ( H1, H2, Details, Clean )
    use HessianModule_0, only: Diff
    use Lexer_core, only: Print_Source
    use Output_m, only: Output, NewLine
    use String_Table, only: Display_String

    type (Hessian_T), intent(in) :: H1, H2
    integer, intent(in), optional :: Details ! Print details, default 1
      ! <= -3 => no output
      ! -2 => Just name, size, and where created
      ! -1 => Just row and col info for each block
      ! 0 => Dimensions of each block,
      ! 1 => Values in each block
    logical, intent(in), optional :: CLEAN   ! print \size

    integer :: I, J, K    ! Subscripts, loop inductors
    integer :: My_Details
    integer :: TotalSize  ! of all blocks

    my_details = 1
    if ( present(details) ) my_details = details
    if ( my_details <= -3 ) return
    if ( h1%name > 0 .and. h2%name > 0 ) then
      call output ( 'Diffing ' )
      call display_string ( h1%name )
      call output ( ' and ' )
      call display_string ( h2%name )
    end if
    if ( h1%where > 0 .and. h2%where > 0 ) then
      call output ( ', created at ' )
      call print_source ( h1%where )
      call output ( ' and ' )
      call print_source ( h2%where )
    end if
    call newLine
    if ( .not. associated(h1%block) .or. .not. associated(h2%block) ) then
      call output ( '(the hessians have been destroyed)', advance='yes' )
      return
    end if
    ! Check that the two Hessians match somehow
    if ( h1%col%nb /= h1%col%nb .or. h1%row%nb /= h2%row%nb ) then
      call output ( 'the hessians have different shapes', advance='yes' )
      return
    endif
    if ( my_Details < -1 ) return
    totalSize = 0
    do k = 1, h1%col%nb
      do j = 1, h1%col%nb
        do i = 1, h1%row%nb
          if ( h1%block(i,j,k)%kind /= h2%block(i,j,k)%kind ) then
            call output ( 'the hessians at (i,j,k) ')
            call output( (/ i,j,k /) )
            call output( 'are of different kinds', advance='yes' )
            cycle
          endif
          if ( my_details < 0 .or. &
            & ( HIDEBSENTBLOCKS .and. h1%block(i,j,k)%kind == h_absent ) ) cycle
          call output ( i, before='Block at row ' )
          call output ( j, before=' and columns ' )
          call output ( k, before=', ', after=' (' )
          if ( h1%row%vec%quantities(h1%row%quant(i))%template%name /= 0 ) then
            call display_string ( &
              & h1%row%vec%quantities(h1%row%quant(i))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h1%row%Inst(i) )
          call output (' , ')
          if ( h1%col%vec%quantities(h1%col%quant(j))%template%name /= 0 ) then
            call display_string ( &
              & h1%col%vec%quantities(h1%col%quant(j))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h1%col%Inst(j) )
          call output (' , ')
          if ( h1%col%vec%quantities(h1%col%quant(k))%template%name /= 0 ) then
            call display_string ( &
              & h1%col%vec%quantities(h1%col%quant(k))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h1%col%Inst(k) )
          call output ( ' )' )
          call newLine
          call Diff ( h1%block(i,j,k), h2%block(i,j,k), details=my_details, clean=clean )
        end do
      end do
    end do

  end subroutine Diff_Hessians

  ! ------------------------------------------------- Dump_Hessian -----
  subroutine Dump_Hessian ( H, Name, Details, onlyTheseBlocks, Clean )
    use HessianModule_0, only: Dump
    use Lexer_core, only: Print_Source
    use MatrixModule_1, only: Dump_RC
    use MLSStringLists, only: switchDetail
    use Output_m, only: Output, NewLine
    use String_Table, only: Display_String, GET_STRING

    type (Hessian_T), intent(in) :: H
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: Details ! Print details, default 1
      ! <= -3 => no output
      ! -2 => Just name, size, and where created
      ! -1 => Just row and col info for each block
      ! 0 => Dimensions of each block,
      ! 1 => Values in each block
    character(len=*), dimension(3), intent(in), optional :: onlyTheseBlocks
    logical, intent(in), optional :: CLEAN   ! print \size

    integer :: I, J, K    ! Subscripts, loop inductors
    integer :: My_Details
    character(len=32), dimension(3) :: myBlocks
    integer :: TotalSize  ! of all blocks

    my_details = 1
    if ( present(details) ) my_details = details
    if ( my_details <= -3 ) return
    if ( present(name) ) call output ( name )
    myBlocks = '*'
    if ( present(onlyTheseBlocks) ) myBlocks = onlyTheseBlocks
    if ( h%name > 0 ) then
      if ( present(name) ) call output ( ', ' )
      call output ( 'Name = ' )
      call display_string ( h%name )
    end if
    if ( h%where > 0 ) then
      call output ( ', created at ' )
      call print_source ( h%where )
    end if
    call newLine
    if ( my_details > -2 ) then
      call dump_rc ( h%row, 'row', .true. )
      call dump_rc ( h%col, 'col', .true. )
    end if
    if ( .not. associated(h%block) ) then
      call output ( '      (the hessian has been destroyed)', advance='yes' )
      return
    end if
    totalSize = 0
    do k = 1, h%col%nb
      do j = 1, h%col%nb
        do i = 1, h%row%nb
          if ( HIDEBSENTBLOCKS .and. h%block(i,j,k)%kind == h_absent ) cycle
          select case ( h%block(i,j,k)%kind )
          case ( h_sparse )
            if ( associated(h%block(i,j,k)%tuples) ) &
              & totalSize = totalSize + h%block(i,j,k)%TuplesFilled
          case ( h_full )
            if ( associated(h%block(i,j,k)%values) ) &
              & totalSize = totalSize + size(h%block(i,j,k)%values)
          end select
          if ( my_details < 0 ) cycle
          if ( skipThisBlock ( &
            & 1, h%row%vec%quantities(h%row%quant(i))%template%name &
            &  ) ) cycle
          if ( skipThisBlock ( &
            & 2, h%col%vec%quantities(h%col%quant(j))%template%name &
            &  ) ) cycle
          if ( skipThisBlock ( &
            & 3, h%col%vec%quantities(h%col%quant(k))%template%name &
            &  ) ) cycle
          call output ( i, before='Block at row ' )
          call output ( j, before=' and columns ' )
          call output ( k, before=', ', after=' (' )
          if ( h%row%vec%quantities(h%row%quant(i))%template%name /= 0 ) then
            call display_string ( &
              & h%row%vec%quantities(h%row%quant(i))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h%row%Inst(i) )
          call output (' , ')
          if ( h%col%vec%quantities(h%col%quant(j))%template%name /= 0 ) then
            call display_string ( &
              & h%col%vec%quantities(h%col%quant(j))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h%col%Inst(j) )
          call output (' , ')
          if ( h%col%vec%quantities(h%col%quant(k))%template%name /= 0 ) then
            call display_string ( &
              & h%col%vec%quantities(h%col%quant(k))%template%name )
            call output ( ':' )
          else
            call output ( '<No template>:' )
          end if
          call output ( h%col%Inst(k) )
          call output ( ' )' )
          call newLine
          call dump ( h%block(i,j,k), details=my_details, clean=clean )
        end do
      end do
    end do
    contains
    function skipThisBlock ( s, tid ) result( doWe )
      ! Return TRUE only if we are to skip dumping this block
      integer, intent(in) :: s ! 1 : i; 2 : j; 3 : k
      integer, intent(in) :: tid ! qty template id
      logical :: doWe
      character(len=32)  ::     templateNameStr
      ! Executable
      templateNameStr = '(not found)'
      ! call outputNamedValue( 'Checking whether to skip this block', (/ s, tid /) )
      if ( tid > 0 ) call get_string ( tid, templateNameStr )
      doWe = ( myBlocks(s) /= '*' .and. &
        & SwitchDetail( myBlocks(s), trim(templateNameStr), options='-wcf' ) < 0 &
        & )
    end function skipThisBlock

  end subroutine Dump_Hessian

  ! ---------------------------------------- Dump_Hessian_Database -----
  subroutine Dump_Hessian_Database ( H, Details, Clean )
    type (Hessian_T), intent(in) :: H(:)
    integer, intent(in), optional :: Details ! See Dump_Hessian
    logical, intent(in), optional :: CLEAN   ! print \size

    integer :: I

    do i = 1, size(h)
      call dump ( h(i), details=details, clean=clean )
    end do

  end subroutine Dump_Hessian_Database

  ! ------------------------------- Hessian_Vector_Vector_Multiply -----
  subroutine Hessian_Vector_Vector_Multiply ( H, V, P, Scalar, Update )
  !{ Multiply Hessian {\tt H} by vector {\tt V} twice, with a factor of
  !  $\frac12$, giving {\tt P}: $P^k = P^k + \frac12 H^k_{ij} V^i V^j$.
  !  This is the  second-order term of a Taylor series.  {\tt P} is
  !  initially set to zero unless {\tt Update} is present and true.

    use HessianModule_0, only: Multiply
    use MLSKinds, only: RV
    use VectorsModule, only: Vector_T

    type(Hessian_T), intent(in) :: H
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: P
    real(rv), intent(in) :: Scalar
    logical, intent(in), optional :: Update

    integer :: I, J, K                  ! Loop indices
    integer :: IQ, JQ, KQ               ! Quantity indices
    integer :: II, JI, KI               ! Instance indices
    logical :: MYUPDATE                 ! Copy of update

    myUpdate = .false.
    if ( present ( update ) ) myUpdate = update

    ! Error checking assumes templates were created by L2CF.
    if ( h%row%vec%template%name /= p%template%name ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Row of Hessian (H) incompatible with product (P) vector." )
    if ( h%col%vec%template%name /= v%template%name ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Column of Hessian (H) incompatible with factor (V) vector." )

    if ( DEEBUG ) then
      call output( "Performing H-v-v times", advance="yes" )
      call outputNamedValue( "shape(h%block)", shape(h%block) )
      call outputNamedValue( "max(h%block%1st index)", h%row%nb )
      call outputNamedValue( "max(h%block%2nd index)", h%col%nb )
      call outputNamedValue( "max(h%block%3rd index)", h%col%nb )
    endif
    do i = 1, h%col%nb
      iq = h%col%quant ( i )
      ii = h%col%inst ( i )
      do j = 1, h%col%nb
        jq = h%col%quant ( j )
        ji = h%col%inst ( j )
        do k = 1, h%row%nb
          kq = h%row%quant ( k )
          ki = h%row%inst ( k )
          ! Rows of H%block (first subscript) correspond to quantities
          ! of P; columns (second and third subscripts) correspond to
          ! quantities of V.
          if ( DEEBUG ) call outputNamedValue( "(k,i,j)", (/k,i,j/) )
          call multiply ( h%block(k,i,j), &
            & v%quantities(iq)%values(:,ii), &
            & v%quantities(jq)%values(:,ji), scalar, &
            & p%quantities(kq)%values(:,ki), &
            & update )
        end do
      end do
    end do

  end subroutine Hessian_Vector_Vector_Multiply

  ! ----------------------------------------- InsertHessianPlane_1 -----
  subroutine InsertHessianPlane_1 ( H, M, B, EL, MIRROR )
    ! Insert matrix M as a plane (block B, element EL) of the Hessian H
    ! If Mirror is set, populate the transpose set also
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

    ! Check that the Jacobian and Hessian match
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
        if ( myMirror ) then
          call InsertHessianPlane ( H%block(rb,b,cb), M%block(rb,cb), el, &
          & mirroring=.true. )
        end if
      end do
    end do
    
  end subroutine InsertHessianPlane_1

  ! ----------------------------------------------- NullifyHessian -----
  subroutine NullifyHessian_1 ( H )
    ! Given a matrix, nullify all the pointers associated with it
    type(Hessian_T), intent(inout) :: H

    ! Executable code
    h%name = 0
    call nullifyRCInfo ( h%col )
    call nullifyRCInfo ( h%row )
    nullify ( h%block )
  end subroutine NullifyHessian_1

  ! -------------------------------------------- OptimizeHessian_1 -----
  subroutine OptimizeHessian_1 ( H )
  ! Get rid of duplicates and zeroes
    type (Hessian_T), intent(inout) :: H
    integer :: I, J, K
    do k = lbound(h%block,3), ubound(h%block,3)
      do j = lbound(h%block,2), ubound(h%block,2)
        do i = lbound(h%block,1), ubound(h%block,1)
          call optimizeBlock ( h%block(i,j,k) )
        end do
      end do
    end do
  end subroutine OptimizeHessian_1

  ! ------------------------------------------ StreamlineHessian_1 -----
  subroutine StreamlineHessian_1 ( H, ScaleHeight, GeodAngle, Threshold )

    use HessianModule_0, only: StreamlineHessian
    use Intrinsic, only: L_ZETA
    use MLSKinds, only: R8
    use QuantityTemplates, only: QuantityTemplate_T
    ! Given a Hessian, trim off the elements that are further away than indicated in
    ! scale height (for zeta coordinates) or geodAngle.  In each block, replace
    ! ones that are smaller in magnitude than the maximum in the block by
    ! a factor of threshold by zero.
    type (Hessian_T), intent(inout) :: H
    type (HessianElement_T), pointer :: HBU, HBL ! In Upper and Lower triangle
    real(r8), intent(in) :: ScaleHeight ! Negative if not specified
    real(r8), intent(in) :: GeodAngle   ! Negative if not specified
    real(r8), intent(in) :: Threshold   ! Negative if not specified

    type (QuantityTemplate_T), pointer :: Q1, Q2 ! For 1st, 2nd cols of H
    integer :: P1, P2   ! Profile indices for 1st, 2nd cols of H
    integer :: I, J, K  ! Loop counters
    logical :: DROPBLOCK

    do i = 1, h%row%nb
      do j = 1, h%col%nb - 1
        p1 = h%col%inst(j)
        q1 => h%col%vec%quantities ( h%col%quant(j) ) % template
        do k = j + 1, h%col%nb
          p2 = h%col%inst(k)
          q2 => h%col%vec%quantities ( h%col%quant(k) ) % template

          hbu => h%block ( i, j, k )
          hbl => h%block ( i, k, j )
          
          ! Decide whether to drop these on horizontal grounds
          dropBlock = .false.
          if ( q1%stacked .and. q2%stacked .and. ( geodAngle > 0.0 ) ) &
            & dropBlock = abs ( q1%phi(1,p1) - q2%phi(1,p2) ) > geodAngle
          if ( dropBlock ) then
            call ClearBlock ( hbu )
            call ClearBlock ( hbl )
          else
            ! Now clear elements on vertical grounds, but only for simple quantities
            if ( q1%coherent .and. q1%verticalCoordinate == l_zeta .and. &
               & q2%coherent .and. q2%verticalCoordinate == l_zeta .and. &
               & scaleHeight > 0.0 ) then

              call streamlineHessian ( hbu, q1, q2, scaleHeight, threshold )
              call streamlineHessian ( hbl, q2, q1, scaleHeight, threshold )

            end if
          end if
        end do
      end do
    end do

  end subroutine StreamlineHessian_1

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
! Revision 2.10  2010/08/20 23:17:48  pwagner
! May specify which blocks to dump by name
!
! Revision 2.9  2010/08/13 22:06:22  pwagner
! Added diff; skips mention of absent blocks
!
! Revision 2.8  2010/06/29 19:56:40  vsnyder
! Add SCALAR argument instead of buried 0.5 factor
!
! Revision 2.7  2010/06/28 17:02:28  pwagner
! Fixed a few bugs; added debugging output
!
! Revision 2.6  2010/06/23 22:43:15  vsnyder
! Remove nonsense loop, use h%row where it should be used, in multiply
!
! Revision 2.5  2010/05/14 22:45:33  pwagner
! CreateBlock must be public
!
! Revision 2.4  2010/03/26 23:15:45  vsnyder
! Add Threshold to StreamlineHessian
!
! Revision 2.3  2010/03/24 20:38:14  vsnyder
! Add Dump and Optimize.  Replace 'continue' with 'cycle'.
!
! Revision 2.2  2010/02/25 21:14:33  pwagner
! Replaced non-standard component separator '.' with '%'
!
! Revision 2.1  2010/02/25 18:14:02  pwagner
! First commit
!
