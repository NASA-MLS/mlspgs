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
module HessianModule_0          ! Low-level Hessians in the MLS PGS suite
!=============================================================================

! This module provides the elementary Hessian type.  Blocks of this
! type are used to compose block Hessians.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
    & Test_Allocate, Test_Deallocate
  use HyperSlabs, only: Repopulate
  use, Intrinsic :: Iso_C_Binding, only: C_Intptr_T, C_Loc
  use Diff_1, only: Diff
  use Dump_0, only: Dump
  use HighOutput, only: OutputnamedValue
  use, Intrinsic :: Iso_C_Binding, only: C_Intptr_T, C_Loc
  use MLSKinds, only: Rh=>rm ! Renamed Here To Make It Easier To Change Later
  use Output_M, only: Output
  use Trace_M, only: Trace_Begin, Trace_End

  implicit none
  private

  public :: H_Absent, H_Sparse, H_Full, H_Unknown
  public :: HessianElement_T, Tuple_T
  public :: ClearBlock, CloneBlock, CopyBlock, CreateBlock
  public :: Densify, DestroyBlock, Diff, Dump
  public :: InsertHessianPlane, Multiply, OptimizeBlock, RH
  public :: Sparsify, StreamlineHessian, Unsparsify

  integer, parameter :: H_Absent = 0    ! An absent block -- assumed zero
  integer, parameter :: H_Sparse = 1    ! A 3-way indexed sparse representation
  integer, parameter :: H_Full = 2      ! A regular 3-D array representation
  integer, parameter :: H_Unknown = 3   ! We don't know yet
  ! This is used for reading l2pc files where we don't want to load the block
  ! into memory until we know we need it.  The HessianModule_0 code doesn't
  ! understand such blocks, so use them with care.
  ! Why use them at all? Such trickery leads to fragile results.

  type :: Tuple_T
    real(rh) :: H      ! The value of a Hessian element
    integer :: I, J, K ! Indices of a nonzero element.  I is the "up" index
                       ! and J and K are the "down" indices"
  end type Tuple_T

  type HessianElement_T
    integer :: nRows = 0   ! Extent of first dimensions of H
    integer :: nCols1 = 0  ! Extent of second dimension of H
    integer :: nCols2 = 0  ! Extent of third dimension of H
    integer :: Kind = H_Absent  ! One of H_Absent, H_Sparse, H_Full, H_Unknown
    ! If Kind == H_Sparse the Tuples component is associated with an array of
    ! tuples, each one giving a row index, two column indices, and a value, and
    ! the Values component is null.  All values other than the ones described
    ! by tuples are assumed to be zero.
    ! If Kind == H_Full the Values component is associated with an array of
    ! shape (nRows,nCols1,nCols2), each element containing a value, and the
    ! Tuples component is null.
    ! Otherwise, both Tuples and Values are null, and nRows == nCols1 ==
    ! nCols2 == 0.
    type(tuple_t), pointer :: Tuples(:) => NULL() ! Some may not be filled (yet)
    ! Some explanation of why and when TuplesFilled would differ from
    ! size(tuples) would not be amiss
    ! So far the only time the original code expanded TuplesFilled
    ! beyond zero was during insertHessianPlane operations
    ! Needless to say, this was exceptionally inconvenient when loading
    ! the values(1:TuplesFilled) from a saved dataset
    ! Therefore we also set TuplesFilled when we create a Hessian Block
    integer :: TuplesFilled = 0 ! Number of tuples filled
    real(rh), pointer :: Values(:,:,:) => NULL() ! for full explicit representation
    logical :: optimizedAlready = .false. ! Have we been through OptimizeBlock?
  end type HessianElement_T

  interface ClearBlock
    module procedure ClearHessianBlock_0
  end interface

  interface CreateBlock
    module procedure CreateHessianBlock_0
  end interface

  interface Densify
    module procedure Densify_Hessian
  end interface

  interface DestroyBlock
    module procedure ClearHessianBlock_0
  end interface

  interface Diff
    module procedure Diff_Hessian_Blocks
  end interface

  interface Dump
    module procedure Dump_Hessian_Block
  end interface

  interface InsertHessianPlane
    module procedure InsertHessianPlane_Array
    module procedure InsertHessianPlane_Matrix
  end interface

  interface Multiply
    module procedure Hessian_Vec_Vec_Multiply_D
    module procedure Hessian_Vec_Vec_Multiply_S
  end interface

  interface Sparsify
    module procedure Sparsify_Hessian, Sparsify_Hessian_Array_D
    module procedure Sparsify_Hessian_Array_S
  end interface

  interface StreamlineHessian
    module procedure StreamlineHessian_0
  end interface

  real(rh), parameter:: AUGMENTFACTOR = 1.0
  logical, parameter :: DEEBUG = .false.
  logical, parameter :: DONTOPTIMIZE = .false.
  logical, parameter :: DONTOPTIMIZEFULLS = .false.
  logical, parameter :: DONTSTREAMLINEFULLS = .false.
  logical, parameter :: DUMPASUNSPARSIFIED = .true.
  logical, parameter :: OPTIMIZEMEANSSPARISFY = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------  AugmentHessian  -----
  subroutine AugmentHessian ( H, N, factor )
    ! Make space for N extra elements in the Hessian tuple if it's not
    ! already there If factor is present then don't merely add 'just
    ! enough' values, but increase the number by this fraction
    ! ????????????? Why is that a good idea ?????????????
    ! Among other consequences, h%tuplesFilled will be smaller
    ! than h%tuples, and subsequent elemental procedures passed
    ! h%tuples w/o an array subsection will be operating
    ! with undefined values (paw)
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    
    type(HessianElement_T), intent(inout) :: H
    integer, intent(in) :: N
    real(rh), intent(in), optional :: FACTOR

    ! Local variables
    real(rh) :: myFactor
    type(tuple_t), dimension(:), pointer :: oldTuple  ! Used to expand 'in place'
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: extra                    ! How many to add
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: S                        ! Size in bytes of an object to deallocate
    integer :: stat                     ! Status from allocates etc
    integer :: space                    ! How many tuples are free

    call trace_begin ( me, 'AugmentHessian', cond=.false. )
    myFactor = 1.0_rh
    if ( present(factor) ) myFactor = factor
    if ( h%kind == h_full ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Attempt to agument an already densified Hessian block" )

    if ( h%kind == h_absent .or. .not. associated(h%tuples) ) then
      h%kind = h_sparse
      allocate ( h%tuples ( 0 ), stat=stat )
      call test_allocate ( stat, moduleName, "H%Tuples", (/1/), (/0/) )
      h%tuplesFilled = 0
    end if

    space = size ( h%tuples ) - h%tuplesFilled
    if ( n > space ) then
      ! Save old tuple
      oldTuple => h%tuples

      ! Work out how much to add
      if ( myFactor > 1.0_rh ) then
        extra = max ( N, nint ( size ( h%tuples ) * factor ) )
      else
        extra = N
      end if

      ! Create new tuple
      s = size ( oldTuple ) + extra
      allocate ( H%tuples ( s ), stat=stat )
      addr = 0
      if ( stat == 0 ) addr = transfer(c_loc(H%tuples(1)), addr)
      call test_allocate ( stat, moduleName, "H%Tuple", &
        & (/1/), (/s/), storage_size(oldTuple) / 8, address=addr )

      ! Copy old contents
      if ( associated ( oldTuple ) ) then
        H%tuples ( 1:h%tuplesFilled ) = oldTuple ( 1:h%tuplesFilled )

      ! Destroy old contents
        s = size(oldTuple) * storage_size(oldTuple) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(oldTuple(1)), addr)
        deallocate ( oldTuple, stat=stat )
        call test_deallocate ( stat, moduleName, "oldTuple", s, address=addr )
      end if
    end if
    call trace_end ( 'AugmentHessian', cond=.false. )
  end subroutine AugmentHessian

  ! -------------------------------------------------  CloneBlock  -----
  subroutine CloneBlock ( Z, X, ForWhom ) ! Z = X, except the values
  ! Duplicate a Hessian block, including copying all of its structural
  ! descriptive information, but not its values.
    type(HessianElement_T), intent(inout) :: Z ! intent(inout) so that
      !                            destroyBlock gets a chance to clean up surds
    type(HessianElement_T), intent(in) :: X
    character(len=*), intent(in) :: ForWhom

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: stat

    ! Executable
    call destroyBlock ( z )
    z%nRows  = x%nRows
    z%nCols1 = x%nCols1
    z%nCols2 = x%nCols2
    z%kind   = x%kind
    z%tuplesFilled = 0
    select case ( x%kind )
    case ( h_absent )
      call CreateEmptyBlock ( z )
    case ( h_sparse )
      z%tuplesFilled = x%tuplesFilled
      allocate ( z%tuples ( x%tuplesFilled ), stat=stat )
      addr = 0
      if ( stat==0 .and. x%tuplesFilled>0 ) &
        & addr = transfer(c_loc(z%tuples(1)), addr)
      call test_allocate ( stat, moduleName, 'Z%Tuples', ubounds=x%tuplesFilled, &
        elementSize = storage_size(z%tuples) / 8, address=addr )
    case ( h_full )
      call Allocate_test ( z%values, x%nRows, x%nCols1, x%nCols2, &
        & 'values in loneBlock_0' // ForWhom, &
        & moduleName )
    ! case default ! Nothing more to be done
    end select
  end subroutine CloneBlock

  ! --------------------------------------------------  CopyBlock  -----
  subroutine CopyBlock ( Z, X ) ! Destroy Z, deep Z = X, including the values
    type(HessianElement_T), intent(inout) :: Z ! intent(inout) so that the
      !                            destroyBlock in cloneBlock gets a chance
      !                            to clean up surds
    type(HessianElement_T), intent(in) :: X
    call CloneBlock ( Z, X, "CopyBlock" )
    select case ( x%kind )
    case ( h_absent ) ! Nothing more to be done
    case ( h_sparse )
      z%tuples = x%tuples
    case ( h_full )
      z%values = x%values
    ! case default ! Nothing more to be done
    end select
  end subroutine CopyBlock

  ! ----------------------------------------  ClearHessianBlock_0  -----
  subroutine ClearHessianBlock_0 ( H )
    ! Clear a HessianElement_T structure
    type(HessianElement_T), intent(inout) :: H
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Stat

    if ( associated(h%tuples) ) then
      s = size(h%tuples) * storage_size(h%tuples) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, &
        & "H%Tuples in DestroyHessianBlock", s, address=addr )
    end if

    call deallocate_test ( h%values, moduleName, "H%H in DestroyHessianBlock" )

    h%tuplesFilled = 0
    h%kind = h_absent

  end subroutine ClearHessianBlock_0

  ! ---------------------------------------  Create_Empty_Hessian  -----
  subroutine CreateHessianBlock_0 ( H, nRows, nCols1, nCols2, H_kind, initTuples, &
                                  & Fill )
  ! Create an empty HessianElement_T structure

    type(HessianElement_T), intent(inout) :: H ! inout so we can destroy it before
                                        ! its components are nullified by
                                        ! default initialization
    integer, intent(in) :: nRows, nCols1, nCols2
    integer, intent(in) :: H_kind
    integer, intent(in), optional :: initTuples
    real(rh), intent(in), optional :: Fill ! Fill value if H_kind==h_full

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STAT                     ! Status from allocate

    call trace_begin ( me, 'CreateHessianBlock_0', cond=.false. )
    call clearBlock ( H )

    h%nRows = nRows
    h%nCols1 = nCols1
    h%nCols2 = nCols2
    h%kind = h_kind
    h%tuplesFilled = 0
    if ( DEEBUG ) then
      call output( 'Creating a Hessian Block', advance='yes' )
      call outputNamedValue( 'nRows', nRows )
      call outputNamedValue( 'nCols1', nCols1 )
      call outputNamedValue( 'nCols2', nCols2 )
      call outputNamedValue( 'h_kind', h_kind )
    end if

    if ( h_kind == h_sparse ) then
      h%tuplesFilled = 0
      if ( present ( initTuples ) ) then
        if ( DEEBUG ) then
          call outputNamedValue( 'initTuples', initTuples )
        end if
        allocate ( h%tuples ( initTuples ), stat=stat )
        addr = 0
        if ( stat==0 .and. initTuples>0 ) addr = transfer(c_loc(h%tuples(1)), addr)
        call test_allocate ( stat, moduleName, "H%Tuples", &
          & (/1/), (/initTuples/), elementSize=storage_size(h%tuples)/8, &
          & address=addr )
        h%tuplesFilled = initTuples
      end if
    end if
    if ( h_kind == h_full ) then
      call Allocate_test ( h%values, nRows, nCols1, nCols2, 'values in CreateHessianBlock_0', &
        & moduleName, fill=fill )
    end if
    if ( DEEBUG ) then
      call dump( h%values, 'h%values' )
    end if
    call trace_end ( 'CreateHessianBlock_0', cond=.false. )
  end subroutine CreateHessianBlock_0

  ! --------------------------------------------  Densify_Hessian  -----
  subroutine Densify_Hessian ( H )
  ! Convert a Hessian represented by tuples to an explicit representation

    type(HessianElement_T), intent(inout) :: H

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: n, status

    ! Skip if already done
    select case ( h%kind )
    case ( h_absent )
    case ( h_sparse )
      call allocate_test ( h%values, h%nRows, h%nCols1, h%nCols2, &
        & "H%values in Densify_Hessian", moduleName, fill=0.0_rh )
      
      if ( associated(h%tuples) ) then
        do n = 1, size(h%tuples)
          h%values(h%tuples(n)%i,h%tuples(n)%j,h%tuples(n)%k) = h%tuples(n)%h
        end do

        n = size(h%tuples) * storage_size(h%tuples) / 8
        addr = 0
        if ( n > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
        deallocate ( h%tuples, stat=status )
        call test_deallocate ( status, moduleName, 'H%Tuples', n, address=addr )
      end if
      h%kind = h_full
    case ( h_full )
    end select

  end subroutine Densify_Hessian

  ! ----------------------------------------  Diff_Hessian_Blocks  -----
  subroutine Diff_Hessian_Blocks ( H1, H2, Details, Indices, options, Clean )
    use MLSStringlists, only: Optiondetail, Unquote
    type(HessianElement_T), intent(inout) :: H1, H2
    integer, intent(in), optional :: Details ! Print details, 0 => minimal,
                                             ! 1 => values, default 1
    integer, intent(in), optional :: Indices(:) ! 3 indices of the block
    character(len=*), intent(in), optional :: options
    logical, intent(in), optional :: CLEAN   ! print \size
    ! Internal variables
    integer :: h1kind
    integer :: h2kind
    logical :: My_Clean
    integer :: My_Details
    logical :: myForce
    ! These are in case we need to diff two blocks stored sparsely
    real(rh), pointer :: h1array(:,:,:) => NULL()
    real(rh), pointer :: h2array(:,:,:) => NULL()

    my_Details = 1
    if ( present(details) ) my_Details = details
    my_Clean = .false.
    if ( present(clean) ) my_Clean = clean
    myForce = .false.
    if ( present(options) ) then
      myForce = ( optionDetail( options, 'f' ) == 'yes' )
    end if
    if ( present(indices) ) then
      call output ( indices(1), before = ' (', after=', ' )
      call output ( indices(2), after=', ' )
      call output ( indices(3), after=')' )
    end if
    if ( h1%nRows /= h2%nRows .or. h1%nCols1 /= h2%nCols1 .or. &
      & h1%nCols2 /= h2%nCols2 ) then
      call output ( 'the hessian blocks have different shapes', advance='yes' )
      return
    end if
    call output ( ' has ' )
    call output ( h1%nRows, after=' Rows, ' )
    call output ( h1%nCols1, after = ' First Columns, ' )
    call output ( h1%nCols2, after = ' Second Columns, is ' )
    h1kind = h1%kind
    h2kind = h2%kind
    if ( h1kind /= h2kind ) then
      if ( .not. myForce .or. my_Details < 1 ) then
        call output ( 'the hessian blocks are of different kind', advance='yes' )
        return
      else
        ! OK, we will need to make the two arrays alike somehow
        ! Here's an idea--make both full
        select case ( h1kind )
        case ( h_absent )
          call Unsparsify( h1 )
        case ( h_sparse )
          call Unsparsify( h1 )
        case default
        end select
        select case ( h2kind )
        case ( h_absent )
          call Unsparsify( h2 )
        case ( h_sparse )
          call Unsparsify( h2 )
        case default
        end select
        call output( '(originally unlike)' ) ! return
      end if
    end if
    if ( h1kind /= h2kind ) then
      call outputNamedValue( 'original kind(h1)', h1kind )
      call outputNamedValue( 'original kind(h2)', h2kind )
      call outputNamedValue( 'new kind(h1)', h1%kind )
      call outputNamedValue( 'new kind(h2)', h2%kind )
    end if
    call Diff_like( h1, h2, details, options )
    ! Now did we have to monkey with h1 and k2 becuase they were unkinds?
    if ( h1kind == h2kind ) return
    if ( h1kind /= h_full ) then
      call deallocate_test ( h1%values, moduleName // 'diff', "h1%values" )
      h1%kind = h1kind
    end if
    if ( h2kind /= h_full ) then
      call deallocate_test ( h2%values, moduleName // 'diff', "h2%values" )
      h2%kind = h2kind
    end if
  contains
    subroutine Diff_like( h1, h2, details, options )
    ! Diff two Hessian blocks, confident that they are alike
    type(HessianElement_T), intent(in) :: H1, H2
    integer, intent(in), optional :: Details ! Print details, 0 => minimal,
                                             ! 1 => values, default 1
    character(len=*), intent(in), optional :: options
    ! Internal variables
    character(len=16) :: myoptions ! Keep only array-diffing ones
    ! Executable
    if ( .not. present(options) ) then
      myOptions = merge('c',' ',my_clean)
    else
      myOptions = trim( &
        & unquote( options, quotes='[', cquotes=']', options='-r' ) &
        & ) // merge('c',' ',my_clean)
    end if
    select case ( h1%kind )
    case ( h_absent )
      call output ( 'both absent, therefore identical', advance='yes' )
    case ( h_sparse )
      call output ( size(h1%tuples), before='sparse with ', &
        & after=' nonzero elements', advance='yes' )
      if ( my_details > 0 ) then
        ! We know how to diff full arrays, so 
        ! temporarily convert sparse representation to full
        call output ( '(Diffing as if full)', advance='yes' )
        call allocate_test( h1array, h1%nRows, h1%nCols1, h1%nCols2, &
          & "h1array in Diff_Hessian_Blocks", ModuleName, fill=0.0_rh )
        call allocate_test( h2array, h2%nRows, h2%nCols2, h2%nCols2, &
          & "h2array in Diff_Hessian_Blocks", ModuleName, fill=0.0_rh )
        call Repopulate( h1array, h1%tuples%i, h1%tuples%j, h1%tuples%k, &
          & h1%tuplesFilled, h1%tuples%h )
        call Repopulate( h2array, h2%tuples%i, h2%tuples%j, h2%tuples%k, &
          & h2%tuplesFilled, h2%tuples%h )
        call Diff ( h1array, 'Hessian 1 values', h2array, 'Hessian 2 values' )
        call deallocate_test ( h1array, moduleName, "h1array in Diff_Hessian_Blocks" )
        call deallocate_test ( h2array, moduleName, "h2array in Diff_Hessian_Blocks" )
      end if
    case ( h_full )
      call output ( 'full', advance='yes' )
      if ( details > 0 ) &
        & call Diff ( h1%values, 'Hessian 1 values', h2%values, &
          & 'Hessian 2 values', options=myoptions )
    case ( h_unknown )
      call output ( 'of unknown form, no method to Diff', advance='yes' )
    end select

    end subroutine Diff_like
  end subroutine Diff_Hessian_Blocks

  ! -----------------------------------------  Dump_Hessian_Block  -----
  subroutine Dump_Hessian_Block ( H, Name, Details, Indices, Clean, Options )
    use MLSFillValues, only: Isnan
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Warning
    type(HessianElement_T), intent(in) :: H
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: Details ! Print details, 
                                             !-3, just check for NaNs,
                                             ! 0 => minimal (shape, bandwidth), 
                                             ! 1 => values, 
                                             ! default 1
    integer, intent(in), optional :: Indices(:) ! 3 indices of the block
    logical, intent(in), optional :: CLEAN   ! print \size
    character(len=*), intent(in), optional :: options ! Passed dumping arrays

    integer :: I
    integer :: My_Details
    logical :: My_Clean
    character(len=32) :: myOptions
    ! In case we need to dump as full blocks stored sparsely
    real(rh), pointer :: harray(:,:,:) => NULL()

    my_Details = 1
    if ( present(details) ) my_Details = details
    my_Clean = .false.
    if ( present(clean) ) my_Clean = clean
    myOptions = ''
    if ( present(options) ) myOptions = options

    if ( present(name) ) call output ( name, advance='yes' )
    call outputNamedValue( 'myOptions', myOptions )
    if ( present(indices) ) then
      call output ( indices(1), before = ' (', after=', ' )
      call output ( indices(2), after=', ' )
      call output ( indices(3), after=')' )
    end if
    call output ( ' has ' )
    call output ( h%nRows, after=' Rows, ' )
    call output ( h%nCols1, after = ' First Columns, ' )
    call output ( h%nCols2, after = ' Second Columns, is ' )
    select case ( h%kind )
    case ( h_absent )
      call output ( 'absent', advance='yes' )
    case ( h_sparse )
      if ( h%TuplesFilled < 1 ) then
        call output( '(0 tuplesFilled)', advance='yes' )
        return
      else if ( .not. associated (h%tuples) ) then
        call output( '(tuplesFilled > 0 yet h%tuples not associated; how?)', advance='yes' )
        return
      end if
      call output ( size(h%tuples), before='sparse with ', &
        & after=' nonzero elements', advance='yes' )
      call outputNamedValue( 'Min, max row numbers', &
        & (/ minval(h%tuples(1:h%TuplesFilled)%i), maxval(h%tuples(1:h%TuplesFilled)%i) /) )
      call outputNamedValue( 'Min, max 1st col numbers', &
        & (/ minval(h%tuples(1:h%TuplesFilled)%j), maxval(h%tuples(1:h%TuplesFilled)%j) /) )
      call outputNamedValue( 'Min, max 2nd col numbers', &
        & (/ minval(h%tuples(1:h%TuplesFilled)%k), maxval(h%tuples(1:h%TuplesFilled)%k) /) )
      if ( my_details < -2 ) then
        ! Check for NaNs
        if ( any(isNaN(h%tuples%h)) ) &
          & call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'NaNs found in sparse Hessian block' )
      else if ( my_details > -1 ) then
        if ( DUMPASUNSPARSIFIED .or. my_details == 0 ) then
          ! temporarily convert sparse representation to full
          call output ( '(Dumping as if full)', advance='yes' )
          call allocate_test( harray, h%nRows, h%nCols1, h%nCols2, &
            & "harray in Dump_Hessian_Blocks", ModuleName, fill=0.0_rh )
          call Repopulate( harray, h%tuples%i, h%tuples%j, h%tuples%k, &
            & h%tuplesFilled, h%tuples%h )
          if ( my_details == 0 ) then
            call Dump ( Reshape(harray, (/h%nCols1, h%nCols2, h%nRows/), &
              & order=(/3,2,1/) ), &
              & 'Hessian values', options='-B' ) ! Just bandwidth
          else
            call Dump ( harray, 'Hessian values', options=options )
          endif
          call deallocate_test ( harray, moduleName, "harray in Dump_Hessian_Blocks" )
        else
          do i = 1, size(h%tuples)
            call output ( h%tuples(i)%i, format='(i5)' )
            call output ( h%tuples(i)%j, format='(i5)' )
            call output ( h%tuples(i)%k, format='(i5)' )
            call output ( '   ' )
            call output ( h%tuples(i)%h, advance='yes' )
          end do
        end if
      end if
    case ( h_full )
      call output ( 'full', advance='yes' )
      if ( my_details < -2 ) then
        ! Check for NaNs
        if ( any(isNaN(h%values)) ) &
          & call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'NaNs found in full Hessian block' )
      else if ( details > 0 ) then
        call dump ( h%values, &
          & options=trim(myOptions) // merge('c',' ',my_clean) )
      else if ( details == 0 ) then
        call Dump ( Reshape(h%values, (/h%nCols1, h%nCols2, h%nRows/), &
          & order=(/3,2,1/) ), &
          & 'Hessian values', options='-B' ) ! Just bandwidth
      end if
    case ( h_unknown )
      call output ( 'of unknown form, nothing to dump', advance='yes' )
      call outputNamedValue ( 'h_unknown', h_unknown )
      call outputNamedValue ( 'h_absent', h_absent )
      call outputNamedValue ( 'h_sparse', h_sparse )
      call outputNamedValue ( 'h_full', h_full )
    end select
    call outputNamedValue ( 'already optimized?', h%optimizedAlready )

  end subroutine Dump_Hessian_Block

  ! ---------------------------------  Hessian_Vec_Vec_Multiply_D  -----
  subroutine Hessian_Vec_Vec_Multiply_D ( H, V1, V2, SCALAR, P, Update )
  !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
  !  {\tt SCALAR}, giving {\tt P}:
    !  $P^i = \text{\tt SCALAR}\times H^i_{jk} V_1^j V_2^k$ or
    !  $P^i = P^i + \text{\tt SCALAR}\times H^i_{jk} V_1^j V_2^k$,
    !  depending upon {\tt Update}.  {\tt P} is initially set to zero
    !  unless {\tt Update} is present and true.

    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0d0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(in) :: Scalar
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    ! logical, parameter :: DEEBUG = .true.
    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_D

  ! ---------------------------------  Hessian_Vec_Vec_Multiply_S  -----
  subroutine Hessian_Vec_Vec_Multiply_S ( H, V1, V2, SCALAR, P, Update )
    !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
    !  {\tt SCALAR}, giving {\tt P}:
    !  $P^i = \text{\tt SCALAR}\times H^i_{jk} V_1^j V_2^k$ or
    !  $P^i = P^i + \text{\tt SCALAR}\times H^i_{jk} V_1^j V_2^k$,
    !  depending upon {\tt Update}.  {\tt P} is initially set to zero
    !  unless {\tt Update} is present and true.

    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0e0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(in) :: Scalar
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    ! logical, parameter :: DEEBUG = .true.
    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_S

  ! -----------------------------------  InsertHessianPlane_Array  -----
  subroutine InsertHessianPlane_Array ( H, plane, k, mirroring )
    ! Insert the supplied array in (:,:,k) of the HessianElement_T
    ! or the (:,k,:) plane if mirroring is present and true
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error

    type(HessianElement_T), intent(inout) :: H
    real(rh), intent(in) :: PLANE(:,:)
    integer, intent(in) :: K
    logical, intent(in), optional :: MIRRORING

    integer :: I, J, N
    logical :: MYMIRRORING
    type(tuple_t) :: oneTuple

    myMirroring = .false.
    if ( present ( mirroring ) ) myMirroring = mirroring
    ! Do some checking
    if ( size ( plane, 1 ) /= h%nRows ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Rows dimension of plane does not match Hessian" )
    if ( myMirroring ) then
      if ( size ( plane, 2 ) /= h%nCols2 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    else
      if ( size ( plane, 2 ) /= h%nCols1 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    end if

    ! Try the dense approach first
    if ( h%kind == h_full ) then
      if ( mirroring ) then
        h%values(:,k,:) = plane(:,:)
      else
        h%values(:,:,k) = plane(:,:)
      end if
    else
      ! Otherwise make space for the new entries
      call AugmentHessian ( H, count ( plane /= 0 ), factor=AUGMENTFACTOR )
      
      n = h%tuplesFilled
      if ( myMirroring ) then
        oneTuple%j = k
      else
        oneTuple%k = k
      end if
      do j = 1, size ( plane, 2 )
        if ( myMirroring ) then
          oneTuple%k = j
        else
          oneTuple%j = j
        end if
        do i = 1, size ( plane, 1 )
          if ( plane ( i, j ) /= 0 ) then
            n = n + 1
            oneTuple%h = plane ( i, j )
            oneTuple%i = i
            h%tuples ( n ) = oneTuple
          end if
        end do
      end do
      h%tuplesFilled = n
    end if

    call optimizeBlock ( h )

  end subroutine InsertHessianPlane_Array

  ! ----------------------------------  InsertHessianPlane_Matrix  -----
  subroutine InsertHessianPlane_Matrix ( H, plane, k, mirroring )
    ! Insert the supplied matrix_0 element in (:,:,k) of the HessianElement_T
    ! or the (:,k,:) plane if mirroring is present and true
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use Matrixmodule_0, only: MatrixElement_T, M_Full, M_Absent, M_Column_Sparse, M_Banded
    use MLSStringlists, only: Switchdetail
    use Toggles, only: Switches

    type(HessianElement_T), intent(inout) :: H
    type(MatrixElement_T), intent(in) :: PLANE
    integer, intent(in) :: K
    logical, intent(in), optional :: MIRRORING

    integer :: I, J, P, R, N
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: MYMIRRORING
    type ( tuple_t ) :: oneTuple

    call trace_begin ( me, 'InsertHessianPlane_Matrix', cond=.false. )
    myMirroring = .false.
    if ( present ( mirroring ) ) myMirroring = mirroring

    ! Do some checking
    if ( plane%nRows /= h%nRows ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Rows dimension of plane does not match Hessian" )
    if ( myMirroring ) then
      if ( plane%nCols /= h%nCols2 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    else
      if ( plane%nCols /= h%nCols1 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    end if

    ! Nothing to do if matrix is empty
    if ( plane%kind == M_Absent ) then
      call trace_end ( 'InsertHessianPlane_Matrix', cond=.false. )
      return
    endif

    ! If the Hessian is dense and the matrix is full then get the other routine to do this
    if ( plane%kind == M_Full ) then
      call InsertHessianPlane ( h, plane%values, k, mirroring )
      if ( switchDetail( switches, 'hess' ) > 0 ) then
        call output( 'Full plane of hessian values inserted at ')
        call output( k, advance='yes' )
      end if
    else
      ! Otherwise make space for the new entries
      if ( h%kind == H_Sparse .or. h%kind == H_Absent ) &
        & call AugmentHessian ( H, size ( plane%values, 1 ), factor=AUGMENTFACTOR )
      
      n = h%tuplesFilled
      if ( myMirroring ) then
        oneTuple%j = k
      else
        oneTuple%k = k
      end if

      select case ( plane%kind ) 
      case ( M_Banded )
        do j = 1, plane%nCols
          if ( myMirroring ) then
            oneTuple%k = j
          else
            oneTuple%j = j
          end if

          i = plane%r1(j) ! Row number of first nonzero
          do p = plane%r2(j-1) + 1, plane%r2(j) ! Indices in values of nonzeros
            oneTuple%i = i
            oneTuple%h = plane%values ( p, 1 )
            i = i + 1
            if ( h%kind == H_Sparse ) then
              n = n + 1
              h%tuples(n) = oneTuple
            else
!print*,'Stats: ', associated ( h%values ), oneTuple%i, oneTuple%j, oneTuple%k
              h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
            end if
          end do
        end do
        if ( switchDetail( switches, 'hess' ) > 0 ) then
          call output( 'Banded plane of hessian values inserted at ')
          call output( oneTuple%k, advance='yes' )
        end if
      case ( M_Column_Sparse )
        do j = 1, plane%nCols
          if ( myMirroring ) then
            oneTuple%k = j
          else
            oneTuple%j = j
          end if
          do r = plane%r1(j-1) + 1, plane%r1(j) ! Indices of nonzeros and their row numbers
            oneTuple%i = plane%r2 ( r ) ! Row number of nonzero
            oneTuple%h = plane%values ( r, 1 )
            if ( h%kind == H_Sparse ) then
              n = n + 1
              h%tuples(n) = oneTuple
            else
              h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
            end if
          end do
        end do
        if ( switchDetail( switches, 'hess' ) > 0 ) then
          call output( 'column sparse plane of hessian values inserted at ')
          call output( oneTuple%k, advance='yes' )
        end if
      end select
      if ( h%kind == H_Sparse ) h%tuplesFilled = n
    end if

    call optimizeBlock ( h )

    call trace_end ( 'InsertHessianPlane_Matrix', cond=.false. )
  end subroutine InsertHessianPlane_Matrix

  ! ----------------------------------------------  OptimizeBlock  -----
  ! Rewrite this from scratch before use
  subroutine OptimizeBlock ( H )
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Warning
    use MLSStringlists, only: Switchdetail
    use Toggles, only: Switches
    use Sort_M, only: Sortp
    ! Delete any 'zero' elements and reorder with i as the slowest changing index
    ! (While this is unconventional, it optimizes for invocations of the 2nd order
    ! Taylor series)
    type (HessianElement_T), intent(inout) :: H

    ! Local variables
    type (Tuple_T), dimension(:), pointer :: NEWTUPLES ! Workspace
    integer, dimension(:), pointer :: ORDER ! The new order for the tuples
    integer, dimension(:), pointer :: RANK ! A 'cost' used to generate rank
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: N        ! Number to end up with
    integer :: I,J,K    ! Loop counters
    integer :: MaxRank  ! Used to place a tuple out of the running
    integer :: Me = -1  ! String index for trace cacheing
    integer :: S        ! Size in bytes of an object to deallocate
    integer :: STATUS   ! Flag from allocate
    logical, parameter :: DEEBUG = .true.
    logical :: verbose

    verbose = ( switchDetail(switches, 'hess') > -1 )
    if ( h%optimizedAlready ) return
    h%optimizedAlready = .true.
    call trace_begin ( me, 'OptimizeBlock', cond=.false. )
    if ( DONTOPTIMIZE ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & "Skipping buggy optimizeBlocks (paw)" )
       call trace_end ( 'OptimizeBlock', cond=.false. )
       return
    else if ( DEEBUG .and. any( h%kind == (/ h_full, h_sparse /) ) .and. &
      & h%tuplesFilled > 0 ) then
      call output( 'About to dump h with details=0', advance='yes' )
      call dump( h, details=0, name='Before optimizing Hessian block' )
      if ( h%kind == h_full ) then
        call outputnamedValue( 'count(h%values /= 0.0)', count(h%values /= 0.0) )
      else
        call outputnamedValue( 'h%tuplesFilled', h%tuplesFilled )
        n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0._rh )
        call outputnamedValue( 'count(h%tuplevalues /= 0.0_rh)', n )
        call outputnamedValue( 'is h sparse?', (h%kind == h_Sparse) )
      end if
    end if
    if ( h%kind == h_Absent .or. h%kind == h_Unknown ) then
      call trace_end ( 'OptimizeBlock', cond=.false. )
      return
    endif
    if ( h%kind == h_Full ) then
      n = count(h%values /= 0.0)
      ! Assuming integers take half the space of real(rh), each tuple takes
      ! 2.5 times the space of a single value.
      if ( 2.5 * n > size(h%values) .or. DONTOPTIMIZEFULLS ) then
        h%optimizedAlready = .true.
        call trace_end ( 'OptimizeBlock', cond=.false. )
        return ! sparsifying won't improve things
      else if ( n < 1 ) then
        ! all values 0; so make it absent
        call clearBlock ( h )
      else
        if ( verbose ) call outputNamedValue( 'Sparsifying full Hessian;' // &
          & 'Number of nozero elements', n )
        call Sparsify_Hessian ( H )
        if ( verbose ) call output( 'Returned from Sparsifying full Hesssian', advance='yes' )
      end if
      if ( h%kind /= h_Sparse ) then
        call trace_end ( 'OptimizeBlock', cond=.false. )
        return
      endif
    elseif ( h%kind == h_Sparse ) then
      ! Now we re-sparsify the block
        if ( verbose ) call outputNamedValue( 'Sparsifying sparse Hessian;' // &
          & 'Number of nozero elements', n )
      call Sparsify( H )
    end if
    ! h%optimizedAlready = .true.
    ! OK, so now it must be sparse

    if ( h%tuplesFilled == 0 ) then
      call clearBlock ( h )
      if ( verbose ) call output( 'We optimized the block away', advance='yes' )
      call trace_end ( 'OptimizeBlock', cond=.false. )
      return
    end if

    if ( .not. associated(h%tuples) ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "How do we get h%tuples disassociated but h%tuplesFilled > 0?" )
      call clearBlock ( h )
    end if

    ! How many will I end up with at most
    n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0._rh )
    if ( n == 0 ) then
      call clearBlock ( h )
      if ( verbose ) call output( 'We optimized the block away', advance='yes' )
      call trace_end ( 'OptimizeBlock', cond=.false. )
      return
    end if

    if ( DEEBUG .and. any( h%kind == (/ h_full, h_sparse /) ) .and. &
      & h%tuplesFilled > 0 ) then
      call dump( h, details=0, name='After optimizing Hessian block' )
      if ( h%kind == h_full ) &
        & call outputnamedValue( 'count(h%values /= 0.0)', count(h%values /= 0.0) )
      call outputnamedValue( 'h%tuplesFilled', h%tuplesFilled )
      if ( associated(h%tuples) ) then
        n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0 )
        call outputnamedValue( 'count(h%tuplevalues /= 0.0)', n )
      else if ( h%tuplesFilled > 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName // "optimizeBlocks", &
          & "How is it that h%tuplesFilled > 0 but h%tuples is not associated" )
      end if
    end if

    if ( OPTIMIZEMEANSSPARISFY ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Skipping remainder of optimizeBlocks after sparsifying(paw)" )
      call trace_end ( 'OptimizeBlock', cond=.false. )
      return
    end if
    ! -------------------------------------------------------------
    ! The rest of this is probably useless as well as wrong
    call allocate_Test ( rank,  h%tuplesFilled, "Rank in OptimizeBlock", ModuleName )
    call allocate_Test ( order, h%tuplesFilled, "Order in OptimizeBlock", ModuleName )
    rank = h%tuples(1:h%tuplesFilled)%k + &
      & h%nCols2 * ( h%tuples(1:h%tuplesFilled)%j + &
      &   h%nCols1 * h%tuples(1:h%tuplesFilled)%i )

    ! Put the zeros out of the running
    maxRank = h%nRows * h%nCols1 * h%nCols2 + 1
    where ( h%tuples ( 1:h%tuplesFilled ) % h == 0 )
      rank = maxRank
    end where

    call SortP ( rank, 1, h%tuplesFilled, order )

    ! If the first sorted one has a zero value, they all do, since zeros
    ! were given maxRank.
    ! ?????????? Why ???????????
    ! We did not sort by value, remember, but by "rank"
    ! This is just one of what are many unjustified assumptions here
    ! They pretty much doom this subroutine (paw)
    if ( h%tuples(order(1))%h == 0.0 ) then
      call clearBlock ( h )
    else

      allocate ( newTuples ( n ), stat=status )
      addr = 0
      if ( status == 0 .and. n > 0 ) addr = transfer(c_loc(newTuples(1)), addr)
      call test_allocate ( status, moduleName, "newTuples", (/1/), (/n/), &
        & elementSize=storage_size(newTuples)/8, address=addr )

      ! Store the sorted tuples while eliminating duplicates
      newTuples(1) = h%tuples(order(1))
      i = 1 ! The one to test for duplication
      k = 1 ! How many unique ones so far
o:    do while ( i < n )
        do j = i+1, n ! The one to test whether it's a duplicate
          if ( rank(order(i)) /= rank(order(j)) ) then ! Found next unique one
            i = j
            k = k + 1
            newTuples(k) = h%tuples(order(j)) ! Store the next unique one
            cycle o
          end if
        end do
        exit ! The ones on the end are duplicates
      end do o

      s = size(h%tuples) * storage_size(h%tuples) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
      deallocate ( h%tuples, stat=status )
      call test_deallocate ( status, moduleName, "h%tuples", s, address=addr )

      h%tuples => newTuples
      h%tuplesFilled = k

      ! Densify if that would take less space
      if ( 2.5 * n > h%nRows * h%nCols1 * h%nCols2 ) then
        call densify_Hessian ( h )
      else if ( n > 1.2 * k ) then
        ! We have an inordinately greater number of tuples than needed
        allocate ( newTuples ( k ), stat=status )
        addr = 0
        if ( status == 0 .and. k > 0 ) addr = transfer(c_loc(newTuples(1)), addr)
        call test_allocate ( status, moduleName, "newTuples", (/1/), (/k/), &
          & elementSize=storage_size(newTuples)/8, address=addr )
        newTuples(:k) = h%tuples(:k)
        s = size(h%tuples) * storage_size(h%tuples) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
        deallocate ( h%tuples, stat=status )
        call test_deallocate ( status, moduleName, "h%tuples", s, address=addr )
        h%tuples => newTuples
      end if
    end if

    call Deallocate_Test ( rank, "Rank in OptimizeBlock", ModuleName )
    call Deallocate_Test ( order, "Order in OptimizeBlock", ModuleName )

    if ( DEEBUG .and. any( h%kind == (/ h_full, h_sparse /) ) .and. &
      & h%tuplesFilled > 0 ) then
      call dump( h, details=0, name='After optimizing Hessian block' )
      if ( h%kind == h_full ) &
        & call outputnamedValue( 'count(h%values /= 0.0)', count(h%values /= 0.0) )
      call outputnamedValue( 'h%tuplesFilled', h%tuplesFilled )
      if ( associated(h%tuples) ) then
        n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0 )
        call outputnamedValue( 'count(h%tuplevalues /= 0.0)', n )
      else if ( h%tuplesFilled > 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName // "optimizeBlocks", &
          & "How is it that h%tuplesFilled > 0 but h%tuples is not associated" )
      end if
    end if
    ! -------------------------------------------------------------
    call trace_end ( 'OptimizeBlock', cond=.false. )
  end subroutine OptimizeBlock

  ! --------------------------------------  Sparsify_Full_Hessian  -----
  subroutine Sparsify_Full_Hessian ( H )
    ! Sparsify the representation of a full H

    use MLSStringlists, only: Switchdetail
    use Toggles, only: Switches

    type (HessianElement_T), intent(inout) :: H
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i, j, k, l, n
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: S                        ! Size in bytes of an object to deallocate
    integer :: status
    logical :: verbose
    ! Executable
    if ( h%kind /= h_full .or. .not. associated(h%values) ) return
    verbose = ( switchDetail(switches, 'hess') > -1 ) .or. .true.
    call trace_begin ( me, 'SparsifyFullHessian', cond=.false. )
    n = count(h%values /= 0)
    if ( verbose ) &
      & call outputNamedValue ( 'About to Sparsify Hessian with non-zeros', n )
    h%kind   = H_Sparse
    h%nRows  = size ( h%values, 1 )
    h%nCols1 = size ( h%values, 2 )
    h%nCols2 = size ( h%values, 3 )

    if ( associated(h%tuples) ) then
      call outputNamedValue ( 'Must deallocate old tuples array', h%tuplesFilled )
      s = size(h%tuples) * storage_size(h%tuples) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
      deallocate ( h%tuples, stat=status )
      call test_deallocate ( status, ModuleName // '%SparsifyFullHessian', &
        & "H%Tuples", s, address=addr )
    end if
    call outputNamedValue ( 'About to allocate Hessian tuples', n )
    allocate ( h%tuples(n), stat=status )
    addr = 0
    if ( status == 0 .and. n > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
    call test_allocate ( status, ModuleName // '%SparsifyFullHessian', "H%Tuples", &
      & (/ 1 /), (/ n /), elementSize=storage_size(h%tuples), address=addr )

    l = 0
    do k = 1, h%nCols2
      do j = 1, h%nCols1
        do i = 1, h%nRows
          if ( h%values(i,j,k) /= 0._rh ) then
            l = l + 1
            h%tuples(l) = tuple_t(h%values(i,j,k),i,j,k)
          end if
        end do
      end do
    end do
    h%tuplesFilled = l
    call deallocate_test ( h%values, "H%Values in Sparsify_Hessian", moduleName )
    call trace_end ( 'SparsifyFullHessian', cond=.false. )
  end subroutine Sparsify_Full_Hessian

  ! -------------------------------------------  Sparsify_Hessian  -----
  subroutine Sparsify_Hessian ( H )
    ! (Re-)Sparsify the representation of H
    ! depending on whether its representation is currently Full
    ! or already Sparse

    use MLSStringlists, only: Switchdetail
    use Toggles, only: Switches

    type (HessianElement_T), intent(inout) :: H
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i
    integer :: Me = -1               ! String index for trace cacheing
    integer :: newSize
    integer :: S                     ! Size in bytes of an object to deallocate
    integer :: status
    type(tuple_t), pointer :: Tuples(:) => NULL()
    logical :: verbose

    ! Executable
    ! h%optimizedAlready = .false.
    call trace_begin ( me, 'Sparsify_Hessian', cond=.false. )
    nullify( tuples )
    verbose = ( switchDetail(switches, 'hess') > -1 ) .or. .true.
    select case (h%kind)
    case (h_Full)
      if ( verbose ) call output( 'Sparsifying full hessian', advance='yes' )
      ! if ( associated(h%values) ) call sparsify ( h%values, h )
      ! call deallocate_test ( h%values, "H%Values in Sparsify_Hessian", moduleName )
      if ( associated(h%values) ) call sparsify_full_Hessian( h )
    case (h_Sparse)
      if ( h%tuplesFilled < 1 ) then
        call trace_end ( 'Sparsify_Hessian', cond=.false. )
        return
      end if
      newSize = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0._rh )
      if ( verbose ) call outputNamedValue( 'Sparsifying sparse hessian', newsize )
      if ( newSize ==  h%tuplesFilled ) then
        call trace_end ( 'Sparsify_Hessian', cond=.false. )
        return
      else if ( newSize > 0 ) then
        allocate ( tuples(newSize), stat=status )
        addr = 0
        if ( status==0 ) addr = transfer(c_loc(tuples(1)), addr)
        call test_allocate ( status, moduleName, 'Tuples(NewSize)', &
          & uBounds=[newSize], elementSize=storage_size(tuples), address=addr )
        newSize = 0
        do i = 1, h%tuplesFilled
          if ( h%tuples(i)%h == 0.0 ) cycle
          newSize = newSize + 1
          tuples(newSize) = h%tuples(i)
        end do
      end if
      s = size(h%tuples) * storage_size(h%tuples) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(h%tuples(1)), addr)
      deallocate ( h%tuples, stat=status )
      call test_deallocate ( status, moduleName, 'h%tuples', s, address=addr )
      h%tuplesFilled = newSize
      h%tuples => tuples
    case default
      if ( verbose ) call output( 'Sparsifying ??? hessian', advance='yes' )
      ! We don't need to do anything for these cases, do we?
    end select
    call trace_end ( 'Sparsify_Hessian', cond=.false. )
  end subroutine Sparsify_Hessian

  ! -----------------------------------  Sparsify_Hessian_Array_D  -----
  subroutine Sparsify_Hessian_Array_D ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_D

  ! -----------------------------------  Sparsify_Hessian_Array_S  -----
  subroutine Sparsify_Hessian_Array_S ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_S

  !  ---------------------------------------  StreamlineHessian_0  -----
  subroutine StreamlineHessian_0 ( H, Q1, Q2, Surface, ScaleHeight, Threshold )
    ! Remove elements 
    ! (1) between levels separated vertically by more than ScaleHeight (units); or
    ! (2) between levels separated vertically by more than Surface (index); or
    ! (3) are smaller in magnitude than the maximum in the block by a factor
    ! of threshold
    ! This will make the Hessians, and the containing L2PC files much smaller
    ! However, some of the assumptions are dubious
    ! On balance, StreamlineHessian_0 might not be worth the risk (paw)
    ! Breaking news: running retrievals with analytical, full hessian blocks
    ! takes too long; we must try streamlining after all
    use MLSKinds, only: R8
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Warning
    use MLSStringlists, only: Switchdetail
    use QuantityTemplates, only: QuantityTemplate_T
    use Toggles, only: Switches
    ! Args
    type (HessianElement_T), intent(inout) :: H
    type (QuantityTemplate_T), intent(in) :: Q1, Q2 ! For 1st, 2nd cols of H
    integer, intent(in) ::  Surface     ! Diagonal horizon (set to 0 outside)
    ! Eventually we'll use the surfaces field of the template to
    ! establish a diagonal horizon if we supply scaleheight
    ! For now, however, we treat it just like surface
    real(r8), intent(in) :: ScaleHeight ! Diagonal horizon (set to 0 outside)
    real(r8), intent(in) :: Threshold
    ! Internal variables
    real(r8) :: Cutoff ! Threshold * maxval(abs(...)) if threshold > 0.
    integer :: I, J, K    ! Subscripts
    integer :: nt      ! Num of nTuples
    integer :: nnz     ! Num of non-zero values
    integer :: S1, S2  ! Surface indices for 1st, 2nd cols of H
    logical :: verbose

    call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & "Streamlining Hessians is probably still buggy--use with caution (paw)" )

    call outputNamedValue ( 'switches', trim(switches) )
    verbose = ( switchDetail(switches, 'hess') > -1 )
    h%optimizedAlready = .false.
    if ( verbose ) call outputNamedValue( 'h%kind', h%kind )
    select case ( h%kind )
    case ( h_absent )
      return
    case ( h_sparse )
      ! Note that this cutoff won't be right if h%tuples is
      ! longer than h%tuplesFilled (paw)
      nt = h%tuplesFilled
      if ( nt < 1 .or. .not. associated(h%tuples) ) then
        h%kind = h_absent
        if ( verbose ) &
          & call output( 'h_sparse, but number of tuplees 0' // &
          & 'changing to h_absent', advance='yes' )
        return
      endif
      if ( nt < 2 ) return
      if ( verbose ) call outputNamedValue( 'Before streamlining sparse Hessian;' // &
        & 'Number of nozero elements', count(h%tuples%h /= 0.) )
      cutoff = threshold * maxval(abs(h%tuples(1:nt)%h))
      if ( surface < 0 .and. scaleHeight < 0._r8 ) then
        ! where ( abs( h%tuples%h ) < cutoff ) h%tuples%h = 0.0
        do i=1, nt
          if ( abs( h%tuples(i)%h ) < cutoff ) h%tuples(i)%h = 0.
        enddo
      else if ( surface > 0 ) then
        ! where ( abs( h%tuples%h ) < cutoff .or. &
         ! &     abs( (h%tuples%j-1 ) / q1%noChans + 1 -   &
         ! &          (h%tuples%k-1 ) / q2%noChans + 1 ) > &
         ! &     surface  ) h%tuples%h = 0.0
        do i=1, nt
          s1 = (h%tuples(i)%j-1 ) / q1%noChans + 1
          s2 = (h%tuples(i)%k-1 ) / q2%noChans + 1
          if ( abs( h%tuples(i)%h ) < cutoff .or. &
            &     abs(s1-s2) > surface ) h%tuples(i)%h = 0.
          ! if ( abs( h%tuples(i)%h ) < cutoff .or. &
          !   &     abs( (h%tuples(i)%j-1 ) / q1%noChans + 1 -   &
          !   &          (h%tuples(i)%k-1 ) / q2%noChans + 1 ) > &
          !   &     surface ) h%tuples(i)%h = 0.
        enddo
      else if ( scaleheight > 0._r8 ) then
        ! where ( abs( h%tuples%h ) < cutoff .or. &
        ! &     abs( q1%surfs( (h%tuples%j-1 ) / q1%noChans + 1, 1) -   &
        ! &          q2%surfs( (h%tuples%k-1 ) / q2%noChans + 1, 1) ) > &
        ! &     scaleHeight  ) h%tuples%h = 0.0
        do i=1, nt
          s1 = (h%tuples(i)%j-1 ) / q1%noChans + 1
          s2 = (h%tuples(i)%k-1 ) / q2%noChans + 1
          if ( abs( h%tuples(i)%h ) < cutoff .or. &
            &     abs(s1-s2) > scaleHeight ) h%tuples(i)%h = 0.
          ! if ( abs( h%tuples(i)%h ) < cutoff .or. &
          !   &     abs( (h%tuples(i)%j-1 ) / q1%noChans + 1 -   &
          !   &          (h%tuples(i)%k-1 ) / q2%noChans + 1 ) > &
          !   &     scaleHeight ) h%tuples(i)%h = 0.
        enddo
      end if
      if ( verbose ) call outputNamedValue( 'After streamlining sparse Hessian;' // &
        & 'Number of nozero elements', count(h%tuples%h /= 0.) )
    case ( h_full )
      if ( DONTSTREAMLINEFULLS ) return
      if ( .not. associated(h%values) ) then
        h%kind = h_absent
        if ( verbose ) &
          & call output( 'h_full, but values not associated' // &
          & 'changing to h_absent', advance='yes' )
        return
      endif
      nnz = count(h%values /= 0.)
      if ( nnz < 1 ) then
        h%kind = h_absent
        if ( verbose ) &
          & call output( 'h_full, but all values 0' // &
          & 'changing to h_absent', advance='yes' )
        return
      endif
      if ( verbose .and. nnz > 0) &
        & call outputNamedValue( 'Before streamlining full Hessian;' // &
        & 'Number of nozero elements', nnz )
      do k = 1, h%nCols2
        s2 = ( k-1 ) / q2%noChans + 1
        do j = 1, h%nCols1
          s1 = ( j-1 ) / q1%noChans + 1
          if ( surface < 0 .and. scaleheight < 0._r8 ) then
            ! no op; some compilers may complain
            cutoff = 0._r8
          else if ( surface > 0 .and. abs(s1-s2) > surface ) then
            h%values ( :, j, k ) = 0
          else if ( scaleHeight > 0._r8 .and. abs(s1-s2) > scaleHeight ) then
            h%values ( :, j, k ) = 0
          ! else if ( scaleHeight > 0._r8 .and. &
          !   & abs ( q1%surfs(s1,1) - q2%surfs(s2,1) ) > scaleHeight ) then
          !   h%values ( :, j, k ) = 0
          end if
        end do
      end do
      if ( threshold > 0._r8 ) then
        cutoff = threshold * maxval(abs(h%values))
        ! where ( abs(h%values) < cutoff ) h%values = 0.0
        do k = 1, h%nCols2
          do j = 1, h%nCols1
            do i = 1, h%nRows
              if ( abs( h%values(i,j,k) ) < cutoff ) h%values(i,j,k) = 0.
            enddo
          enddo
        enddo
      end if
      if ( verbose .and. nnz > 0 ) &
        & call outputNamedValue( 'After streamlining full Hessian;' // &
        & 'Number of nozero elements', count(h%values /= 0.) )
    end select

    call optimizeBlock ( h )

  end subroutine StreamlineHessian_0

! =====     Private Procedures     =====================================

  ! -------------------------------------------  CreateEmptyBlock  -----
  subroutine CreateEmptyBlock ( EmptyBlock )
    type(HessianElement_T), intent(inout) :: EmptyBlock
    ! EmptyBlock is intent(inout) so that destroyBlock will have a chance
    ! to clean up surds.  Default initialization for intent(out) would
    ! nullify the pointers before destroyBlock had a chance to deallocate
    ! them.
    call destroyBlock ( emptyBlock )
  end subroutine CreateEmptyBlock

  !  ------------------------------------------------  Unsparsify  -----
  ! Fill out h%values based on its sparse nTuples (or else zeros)
  ! Should you combine this with Repopulate from MLSFillValues module?
  subroutine Unsparsify ( H )
    ! Argument
    type(HessianElement_T), intent(inout) :: H
    ! Internal variables
    integer :: i
    ! Executable
    select case (h%kind)
    case ( h_absent )
      if ( associated(h%values) ) &
        & call Deallocate_Test ( h%values, moduleName // 'Unsparsify', &
        &  'h%values' )
      call allocate_test( h%values, h%nRows, h%nCols1, h%nCols2, &
          & "h%values", ModuleName // 'Unsparsify', fill=0.0_rh )
    case ( h_sparse )
      if ( associated(h%values) ) &
        & call Deallocate_Test ( h%values, moduleName // 'Unsparsify', &
        &  'h%values' )
      call allocate_test( h%values, h%nRows, h%nCols1, h%nCols2, &
          & "h%values", ModuleName // 'Unsparsify', fill=0.0_rh )
      do i=1, h%TuplesFilled
        h%values(h%tuples(i)%i, h%tuples(i)%j, h%tuples(i)%k) = h%tuples(i)%h
      end do
    case default
    end select
  end subroutine Unsparsify

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HessianModule_0

! $Log$
! Revision 2.32  2017/11/03 20:00:45  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.31  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.30  2015/03/28 01:04:29  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.29  2014/09/04 23:42:20  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.28  2014/01/09 00:25:06  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.27  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.26  2012/02/02 01:11:04  pwagner
! Added Procedures to Clone, Copy blocks just like with matrices
!
! Revision 2.25  2012/01/13 01:08:18  pwagner
! Fixed two bugs in StreamlineHessian
!
! Revision 2.24  2011/12/07 01:19:56  pwagner
! Details=0 now dumps Bandwidth, too
!
! Revision 2.23  2011/11/10 16:20:00  pwagner
! Fixed bug in streamline
!
! Revision 2.22  2011/10/14 00:35:06  pwagner
! Further attempts to stop hanging during Streamline
!
! Revision 2.21  2011/10/07 00:01:55  pwagner
! Some improvements to speed; still hangs though in Streamline
!
! Revision 2.20  2011/09/20 22:34:06  pwagner
! Repaired most obvious bugs in Sparsify, Streamline
!
! Revision 2.19  2011/04/02 01:20:57  vsnyder
! Don't look at h%values for sparse Hessians
!
! Revision 2.18  2011/04/01 22:06:34  vsnyder
! Don't look at h%values for sparse Hessians, or h%tuplesFilled for full ones
!
! Revision 2.17  2011/03/02 02:01:35  vsnyder
! Declare local variables of Hessian_Vec_Vec_Multiply_* in the include file
!
! Revision 2.16  2011/02/18 17:52:18  pwagner
! myOptions needed to be longer; long enough?
!
! Revision 2.15  2011/02/05 01:37:31  pwagner
! Passes options to dump routines
!
! Revision 2.14  2010/11/25 01:16:32  pwagner
! Fixed bug in diffing Hessian blocks with options
!
! Revision 2.13  2010/11/19 23:54:19  pwagner
! Turned optimizeBlock back on, hopefully skipping bauggy parts
!
! Revision 2.12  2010/11/05 20:25:10  vsnyder
! Rename RM as RH to make it easier to change later.  Rename KIND argument
! of CreateHessianBlock_0 as H_Kind so the KIND intrinsic is available.  Add
! a Fill optional argument.  Set h%kind = h_full in Densify_Hessian.  Delete
! unused declarations.
!
! Revision 2.11  2010/09/16 23:54:43  pwagner
! dump with details=-3 warns of NaNs
!
! Revision 2.10  2010/08/24 18:04:22  yanovsky
! Add default initialization to the components of HessianElement_T
!
! Revision 2.9  2010/08/20 23:06:05  pwagner
! Sets Augmenting factor to 1.0; skips optimizing blocks
!
! Revision 2.8  2010/08/13 22:05:32  pwagner
! Added diff; dumps sparse as if full
!
! Revision 2.7  2010/06/29 19:56:40  vsnyder
! Add SCALAR argument instead of buried 0.5 factor
!
! Revision 2.6  2010/06/28 17:01:56  pwagner
! Fixed a few bugs; added debugging output
!
! Revision 2.5  2010/06/22 22:45:38  vsnyder
! Correct LaTeX comments for multiply
!
! Revision 2.4  2010/03/26 23:15:45  vsnyder
! Add Threshold to StreamlineHessian
!
! Revision 2.3  2010/03/24 20:36:24  vsnyder
! Add Dump_Hessian_Block.  Replace specific DestroyHessianBlock_0 with
! ClearHessianBlock_0 but keep generic DestroyBlock.  Add OptimizeBlock
! at the end of InsertHessianPlane_Array.  Simplify InsertHessianPlane_Matrix.
! Simplify OptimizeBlock.
!
! Revision 2.2  2010/02/25 18:34:11  pwagner
! Fixed bug in first commit
!
! Revision 2.1  2010/02/25 18:13:35  pwagner
! First commit
!
