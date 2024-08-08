! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ChunkDivide_m

  ! This module replaces the old ScanDivide, and is a new approach to dividing
  ! the data into chunks.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use ChunkDivideConfig_M, only: ChunkDivideConfig_T, ChunkDivideConfig, Dump, &
    & Dump_CriticalSignals
  use Chunks_M, only: MLSChunk_T, Dump
  use Dump_0, only: Dump
  use HighOutput, only: OutputnamedValue
  use Init_Tables_Module, only: F_Crashifphinotmono, &
    & F_Criticalbands, F_Criticalmodules, F_Criticalsignals, &
    & F_Excludepostoverlaps, F_Excludeprioroverlaps, &
    & F_Homemodule, F_Homegeodangle, &
    & F_Maxlength, F_Maxorby, F_Method, F_Module, F_Nochunks, F_Noslaves, &
    & F_Overlap, F_Loweroverlap, &
    & F_Saveobstructions, F_Scanlowerlimit, F_Scanupperlimit, &
    & F_SkipL1Bcheck, F_Upperoverlap, &
    & Field_First, Field_Last, L_Both, L_Either, L_Even, &
    & L_Fixed, L_Polygon, F_Maxgap, L_Orbital, L_Pe, S_Chunkdivide, S_Dump
  use Intrinsic, only: L_None, Field_Indices, Lit_Indices, Phyq_Angle, &
    & Phyq_Dimensionless, Phyq_Length, Phyq_Mafs, Phyq_Time
  use L1bData, only: L1bData_T, Namelen, Precisionsuffix, &
    & AssembleL1Bqtyname, Dump, DeallocateL1BData, GetL1BFile, ReadL1BData
  use MLSCommon, only: MLSFile_T, Range_T, Tai93_Range_T, InRange
  use MLSFinds, only: FindFirst, FindLast, FindLongestRange
  use MLSFiles, only: Dump, GetMLSFilebytype, MLS_OpenFile
  use MLSKinds, only: R8, Rp
  use MLSMessageModule, only: MLSMessage, MLSMsg_Error, MLSMsg_Warning
  use MLSNumerics, only: Hunt
  use MLSFinds, only: FindFirst
  use MLSSignals_M, only: Dump, Modules, GetModuleName
  use MLSStringlists, only: SwitchDetail
  use Output_M, only: Blanks, Output, &
    & RevertOutput, SwitchOutput
  use PCFHdr, only: GlobalAttributes
  use String_Table, only: Get_String, Display_String
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Polygon_M, only: Polygon_Inside, Polygon_Vertices

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
! any_good_signaldata      Return TRUE if and only if any good data in signal
! ChunkDivide              Divide MAFs in processing range among chunks
! DestroyChunkDatabase     Deallocate memory taken by chunk database
! ReduceChunkDatabase      Reduce chunk database to [first, last] chunks
! === (end of toc) ===
! === (start of api) ===
! log any_good_signaldata ( int signal, int sideband,
!  *MLSFile_T filedatabase(:), [int maf], [int maf2], &
!    & [log good_buffer(:)], [int mafRange(2)] )
! ChunkDivide (int root, TAI93_Range_T processingRange,
!    *MLSFile_T filedatabase(:), *mlSChunk_T Chunks(:) )
! DestroyChunkDatabase (*mlSChunk_T Chunks(:) )
! === (end of api) ===
  public :: CHUNKDIVIDE, CFM_CHUNKDIVIDE, &
    & DESTROYCHUNKDATABASE, REDUCECHUNKDATABASE

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

  public :: CHUNKDIVIDECONFIG_T, CHUNKDIVIDECONFIG

  ! The chunk divide methods are:
  !
  ! Fixed - Ignore the L1B file, just give a fixed set of chunks as described
  !
  ! Polygon - Used for QTM, ASMLS
  !
  ! PE - One chunk centred on the given orbit geodetic angle.
  !
  ! Orbital - Chunks are some ideal fraction of an orbit.  The algorithm
  !           desires to keep the chunk boundaries at the same positions each
  !           orbit where possible.
  !
  ! Even - Hope to have chunks all about the same length, as quoted.

  type MAFRange_T
    integer, dimension(2) :: L1BCover ! Range found in L1B files
    integer, dimension(2) :: L2Cover  ! Range covered by L2 ProcessingRange
    integer, dimension(2) :: Expanded ! L2Cover plus prior/post overlaps
  end  type MAFRange_T

  public :: OBSTRUCTION_T
  ! This type describes obstructions in the Level 1 data which will affect the
  ! selection of chunk divisions.
  type Obstruction_T
    logical :: range                    ! If set is a range not a wall
    integer, dimension(2) :: MAFS       ! Affected MAF or MAF range
  end type Obstruction_T

  type (Obstruction_T), dimension(:), public, save, pointer :: OBSTRUCTIONS => null()

  interface dump
    module procedure DUMP_OBSTRUCTIONS
    ! module procedure DUMP_CRITICALSIGNALS
  end interface

  logical, parameter :: CHECKFORSHAREDMAFS = .true.
  integer, parameter :: SEVERITYIFPHINOTMONO = MLSMSG_Error
  logical, parameter :: DONTPAD = .true.
  integer :: swLevel ! How much extra debugging info to print (-1 means none)
  integer, parameter :: VERBOSETHRESHOLD = 1 ! was 0; turn on extra debugging

contains ! ===================================== Public Procedures =====

  !----------------------------------------  DestroyChunkDatabase  -----
  subroutine DestroyChunkDatabase ( Chunks, Shallow )
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS
    logical, intent(in), optional :: Shallow ! Don't deallocate components

    integer(c_intptr_t) :: Addr    ! For tracing
    integer :: CHUNK               ! Index
    logical :: Deep                ! .not. shallow, else .false.
    integer :: Me = -1             ! String index for trace
    integer :: S                   ! Size in bytes of an object to deallocate
    integer :: STATUS              ! From deallocate

    if ( .not. associated ( chunks ) ) return

    call trace_begin ( me, "DestroyChunkDatabase", cond=toggle(gen) )
    deep = .true.
    if ( present(shallow) ) deep = .not. shallow
    if ( deep ) then
      do chunk = 1, size ( chunks )
        call Deallocate_test ( chunks(chunk)%hGridOffsets, &
          & 'chunks(?)%hGridOffsets', ModuleName )
        call Deallocate_test ( chunks(chunk)%hGridTotals, &
          & 'chunks(?)%hGridTotals', ModuleName )
      end do
    end if
    s = size(chunks) * storage_size(chunks) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(chunks(1)), addr)
    deallocate ( chunks, stat=status )
    call test_deallocate ( status, moduleName, "chunks", s, address=addr )
    call trace_end ( "DestroyChunkDatabase", cond=toggle(gen) )
  end subroutine DestroyChunkDatabase

  !----------------------------------------  ReduceChunkDatabase  -----
  subroutine ReduceChunkDatabase ( chunks, firstChunk, lastChunk )

    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    type (MLSChunk_T), dimension(:), pointer :: chunks
    integer, intent(in) :: firstChunk, lastChunk
    ! Local variables
    integer(c_intptr_t) :: Addr    ! For tracing
    type (MLSChunk_T), dimension(:), pointer :: TEMPDATABASE
    integer :: Me = -1             ! String index for trace
    integer :: newSize, status
    ! Executable
    if ( .not. associated ( chunks ) ) return
    call trace_begin ( me, "ReduceChunkDatabase", cond=toggle(gen) )
    newSize = lastChunk - firstChunk + 1
    if ( newSize < 1 ) return
    if ( lastChunk > size(chunks) ) return
    allocate(tempDatabase(newSize), STAT=status)
    addr = 0
    if ( status == 0 .and. newSize>0 ) addr = transfer(c_loc(tempDatabase(1)), addr )
    call test_allocate ( status, moduleName, "tempDatabase", uBounds = newSize, &
      elementSize = storage_size(tempDatabase) / 8, address=addr )

    tempDatabase(1:newSize) = chunks(firstChunk:lastChunk)
    call DestroyChunkDatabase ( chunks, shallow=.true. )
    chunks => tempDatabase
    call trace_end ( "ReduceChunkDatabase", cond=toggle(gen) )

  end subroutine ReduceChunkDatabase

  ! ------------------------------------------------  Chunk Divide -----
  subroutine ChunkDivide ( root, processingRange, filedatabase, chunks )

    use Allocate_Deallocate, only: Test_Deallocate
    use DumpCommand_m, only: DumpCommand
    use Expr_m, only: Expr
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Lexer_Core, only: Print_Source
    use MLSCommon, only: TAI93_Range_T
    use MLSHDF5, only: GetHDF5Attribute
    use MLSL2Options, only: Need_L1BFiles, SpecialDumpfile
    use MLSL2Timings, only: Section_Times, Total_Times
    use MoreTree, only: Get_Field_ID, Get_Spec_ID
    use Next_Tree_Node_m, only: Next_Tree_Node, Next_Tree_Node_State
    use Time_m, only: Time_Now
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Node_ID, NSons, Subtree, Sub_Rosa, Where_At=>Where
    use Tree_Types, only: N_Named

    integer, intent(in) :: ROOT    ! Root of the L2CF tree for ChunkDivide
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type( TAI93_Range_T ), intent(in) :: PROCESSINGRANGE
    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS, TempChunks

    ! Local variables
    integer(c_intptr_t) :: Addr    ! For tracing
    type (MLSFile_T), pointer :: DACSFile
    type (MAFRange_T) :: MAFRange
    integer :: Me = -1                  ! String index for trace
    integer :: S                        ! Size in bytes of an object to deallocate
    integer :: Son                      ! of section root
    type(next_tree_node_state) :: State ! of tree traverser
    integer :: STATUS                   ! From deallocate

    ! For announce_error:
    logical :: DidOne                   ! Found a ChunkDivide spec
    integer :: ERROR                    ! Error level
    integer :: CHUNK                    ! Loop counter
    integer :: NumberOlds               ! Loop counter
    integer :: OLDCHUNK                 ! Loop counter

    integer, parameter :: BadUnits = 1
    integer, parameter :: NotSpecified = BadUnits + 1
    integer, parameter :: Unnecessary = NotSpecified + 1
    integer, parameter :: invalid = Unnecessary + 1
    integer :: instrumentModule
    character (len=NameLen) :: InstrumentModuleName
    type(Range_T), dimension(:), allocatable :: MAFRanges

    ! For timing
    logical :: Timing
    real :: T1

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    nullify(ChunkDivideConfig%criticalSignals)    ! Just for Sun's compiler

    call trace_begin ( me, "ChunkDivide", root, cond=toggle(gen) )

    timing = section_times
    if ( timing ) call time_now ( t1 )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )

    ! First decode the l2cf instructions

    ! Eventually the ChunkDivide command will be free floating, in the meantime
    ! find it within the section
    ! WE CAN GET RID OF THIS BIT WHEN THE COMMAND FLOATS FREE LATER
    didOne = .false.
    do
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit
      if ( node_id(son) == n_named ) son = subtree(2,son) ! Ignore label
      select case ( get_spec_id(son) )
      case ( s_dump )
        call dumpCommand ( son )
      case ( s_chunkDivide )
        call chunkDivideL2CF ( son )
        call dump(ChunkDivideConfig)
        didOne = .true.
      end select
    end do
    if ( .not. didOne ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'ChunkDivide section must include a ChunkDivide Spec' )

    ! For methods other than fixed, we want to survey the L1 data and note the
    ! location of obstructions
    nullify ( obstructions )
    if ( .not. NEED_L1BFILES ) then
      call output( 'Entered ChunkDivide w/o an l1boa file--hope thats OK', &
        & advance='yes' )
    elseif ( ChunkDivideConfig%method /= l_fixed ) then
      call output( 'We must survey the l1b data for this Chunk Divide method', &
        & advance='yes' )
      call SurveyL1BData ( processingRange, filedatabase, mafRange )
      if ( swLevel > -1 ) then
        call output ( 'Requested time range ' )
        call output ( processingRange%startTime )
        call output ( ' : ' )
        call output ( processingRange%endTime, advance='yes' )
        call output ( 'Corresponding MAF range ' )
        call output ( mafRange%L2Cover(1) )
        call output ( ' : ' )
        call output ( mafRange%L2Cover(2), advance='yes' )
        call output ( 'MAF range in L1B files ' )
        call output ( mafRange%L1BCover(1) )
        call output ( ' : ' )
        call output ( mafRange%L1BCover(2), advance='yes' )
        call output ( 'including overlapped days ' )
        call output ( mafRange%Expanded(1) )
        call output ( ' : ' )
        call output ( mafRange%Expanded(2), advance='yes' )
        call output ( 'Method: ' )
        call display_string ( lit_indices(ChunkDivideConfig%method), &
          &             strip=.true., advance='yes' )
      endif
    endif

    ! Now go place the chunks.
    select case ( ChunkDivideConfig%method )
    case ( l_fixed )
      call ChunkDivide_Fixed ( chunks, root )
    case ( l_polygon )
      if ( .not. NEED_L1BFILES ) then
        call AnnounceError ( root, invalid, ChunkDivideConfig%method )
      else
        call ChunkDivide_Polygon ( mafRange, filedatabase, chunks, &
                                 & InstrumentModuleName, root )
      endif
    case ( l_PE )
      if ( .not. NEED_L1BFILES ) then
        call AnnounceError ( root, invalid, ChunkDivideConfig%method )
      else
        call ChunkDivide_PE ( mafRange%Expanded, filedatabase, chunks )
      endif
    case ( l_orbital )
      if ( .not. NEED_L1BFILES ) then
        call AnnounceError ( root, invalid, ChunkDivideConfig%method )
      else
        call ChunkDivide_Orbital ( mafRange, filedatabase, chunks, root )
      endif
    case ( l_even )
      if ( .not. NEED_L1BFILES ) then
        call AnnounceError ( root, invalid, ChunkDivideConfig%method )
      else
        call ChunkDivide_Even ( mafRange%Expanded, filedatabase, chunks )
      endif
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unexpected problem with ChunkDivide' )
    end select

    ! Check that level 1 has done the DACS deconvolution (so we don't have to)
    ! Note that this would be an attribute, a global one, of the DACS
    ! radiance file attached to the group object '/'
    if ( NEED_L1BFILES ) &
      & DACSFile => GetL1BFile ( filedatabase, 'DACsDeconvolved', '-/a' )
    if ( .not. NEED_L1BFILES ) then
      ChunkDivideConfig%DACSDeconvolved = .false.
    elseif ( associated(DACSFile) ) then
      if ( .not. DACSFile%stillOpen ) call mls_OpenFile(DACSFile)
      DACSFile%fileID%sd_id = 0 ! So we don't look here for the attribute
      ! We still need to open the root group '/'
      call h5gopen_f ( DACSFile%fileID%f_id, '/', DACSFile%fileID%grp_id, status )
      call GetHDF5Attribute( DACSFile, 'DACsDeconvolved', &
      & ChunkDivideConfig%DACSDeconvolved )
      call h5gclose_f ( DACSFile%fileID%grp_id, status )
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Attribute DACsDeconvolved not found in probable DACS file' )
      ChunkDivideConfig%DACSDeconvolved = .false.
    endif

    if ( .not. ChunkDivideConfig%DACSDeconvolved ) &
      & call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Failed to confirm that level 1 performed DACS deconvolution--hope' // &
      & ' you are OK with that' )

    ! Tidy up
    if ( swLevel > -1 .or. &
      & switchDetail(switches, 'opt') > -1 ) call dump(ChunkDivideConfig)
    if ( associated(obstructions) ) then
      if ( swLevel > -1 ) call Dump_Obstructions ( obstructions )
      if ( .not. ChunkDivideConfig%saveObstructions ) then
        s = size(obstructions) * storage_size(obstructions) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(obstructions(1)), addr)
        deallocate ( obstructions, stat=status )
        call test_deallocate ( status, moduleName, 'obstructions', s, address=addr )
      end if
    end if
    if ( associated(ChunkDivideConfig%criticalSignals) ) then
      if ( swLevel > -1 ) &
        & call Dump_criticalSignals(ChunkDivideConfig%criticalSignals)
      s = size(ChunkDivideConfig%criticalSignals) * &
        & storage_size(ChunkDivideConfig%criticalSignals) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(ChunkDivideConfig%criticalSignals(1)(1:1)), addr)
      deallocate ( ChunkDivideConfig%criticalSignals, stat=status )
      call test_deallocate ( status, moduleName, &
        & 'ChunkDivideConfig%criticalSignals', s, address=addr )
    end if
    ! Check that no 2 chunks share non-overlapped MAFs
    ! That mafrange(1) <= all_mafs <= mafrange(2)
    ! And that all overlaps are >= 0
    if ( ChunkDivideConfig%method /= l_fixed ) &
      & call CheckChunkSanity ( chunks, mafRange%Expanded )
    if ( .not. associated(chunks) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'ChunkDivide failed to associate the chunks pointer with a target' )
    else if ( size(chunks) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'ChunkDivide failed to produce any chunks' )
    endif
    
    !! Check that no chunk is wholly contained inside another
    NumberOlds = size ( chunks )
    allocate( MAFRanges( NumberOlds ) )
    do chunk = 1, NumberOlds
      MAFRanges(chunk)%Bottom = chunks(chunk)%firstMAFIndex
      MAFRanges(chunk)%Top    = chunks(chunk)%lastMAFIndex
    enddo
    
    call OutputNamedValue('Num Chunks before checking inRange', NumberOlds )
    
    nullify (TempChunks)
    do OldChunk = 1, NumberOlds
      if ( OldChunk == 1 ) then
        if ( any ( inRange( MAFRanges(OldChunk), MAFRanges(2:NumberOlds) ) ) ) cycle
      elseif ( OldChunk < size ( chunks ) ) then
        if ( &
          & any ( inRange( MAFRanges(OldChunk), MAFRanges(1:OldChunk-1) ) ) &
          & .or. &
          & any ( inRange( MAFRanges(OldChunk), MAFRanges(OldChunk+1:NumberOlds) ) ) &
          & ) cycle
      else
        if ( any ( inRange( MAFRanges(OldChunk), MAFRanges(1:OldChunk-1) ) ) ) cycle
      endif
      call AddChunkToDatabase ( TempChunks, Chunks(OldChunk) )
    enddo
    call DestroyChunkDatabase ( chunks, shallow=.true. )
    chunks => TempChunks
    
    call OutputNamedValue('Num Chunks after removing inRange', size ( chunks ) )

    ! Now go through and number the chunks
    do chunk = 1, size ( chunks )
      chunks(chunk)%chunkNumber = chunk
    end do

    if ( swLevel > 0 .or. .true. ) then
      call dump ( chunks )
    else
      call OutputNamedValue( 'Number of chunks', size(chunks) )
    endif
    if ( specialDumpFile /= ' ' ) call revertOutput
    call trace_end ( "ChunkDivide", cond=toggle(gen) )

  contains

  ! ========================================== Internal Procedures =====

    ! ---------------------------------------------- AnnounceError -----
    subroutine AnnounceError ( where, Code, field )
      integer, intent(in) :: where, Code, field

      error = max(error,1)
      call print_source ( where_at(where) )
      call output ( ', ChunkDivide complained: ' )
      select case ( code )
      case ( BadUnits )
        call output ( ' The field ' )
        call display_string ( field_indices(field) )
        call output (' has inappropriate units', advance='yes' )
      case ( notSpecified )
        call output ( ' The parameter ' )
        call display_string ( field_indices(field) )
        call output ( ' is required but not specified.', advance='yes' )
      case ( unnecessary )
        call output ( ' The parameter ' )
        call display_string ( field_indices(field) )
        call output ( ' is specified but not appropriate.', advance='yes' )
      case ( invalid )
        call output ( ' The method ' )
        call display_string ( lit_indices(field) )
        call output ( ' is invalid w/o an l1boa file.', advance='yes' )
      end select
    end subroutine AnnounceError

    !-------------------------------------------- ChunkDivide_Even -----
    subroutine ChunkDivide_Even ( mafRange, filedatabase, chunks, root )
      integer, dimension(2), intent(in) :: MAFRANGE
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS
      integer, intent(in), optional :: Root ! in l2cf tree, for tracing

      ! Local variables
      integer :: Me = -1             ! String index for trace

      ! Exectuable code
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "I'm sorry, we haven't written the code for even chunk divides yet" )
      if ( present(root) ) call trace_begin ( me, 'ChunkDivide_even', root, &
          & cond=toggle(gen) .and. levels(gen) > 2 )
      if ( present(root) ) call trace_end ( &
        & 'ChunkDivide_even', cond=toggle(gen) .and. levels(gen) > 2 )
    end subroutine ChunkDivide_Even

    !------------------------------------------- ChunkDivide_Fixed -----
    subroutine ChunkDivide_Fixed ( chunks, root )
      ! type (ChunkDivideConfig_T), intent(in) :: CONFIG
      use Allocate_Deallocate, only: Test_Allocate
      use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS
      integer, intent(in), optional :: Root ! in l2cf tree, for tracing

      ! Local variables
      integer(c_intptr_t) :: Addr         ! For tracing
      integer :: I                        ! Loop inductor
      integer :: Me = -1                  ! String index for trace
      integer :: STATUS                   ! From allocate
      integer :: MAXLENGTH                ! nint(config%maxLength)
      integer :: LOWEROVERLAP             ! nint(config%lowerOverlap)
      integer :: UPPEROVERLAP             ! nint(config%upperOverlap)

      ! Executable code
      if ( present(root) ) call trace_begin ( me, 'ChunkDivide_fixed', root, &
          & cond=toggle(gen) .and. levels(gen) > 2 )
      swlevel = switchDetail(switches, 'chu' )
      allocate ( chunks(ChunkDivideConfig%noChunks), stat=status )
      addr = 0
      if ( status == 0 .and. ChunkDivideConfig%noChunks > 0 ) &
        & addr = transfer(c_loc(chunks(1)), addr)
      call test_allocate ( status, moduleName, 'chunks', &
        & uBounds = ChunkDivideConfig%noChunks, &
        & elementSize = storage_size(chunks) / 8, address=addr )

      maxLength = nint ( ChunkDivideConfig%maxLength )
      lowerOverlap = nint ( ChunkDivideConfig%lowerOverlap )
      upperOverlap = nint ( ChunkDivideConfig%upperOverlap )
      do i = 1, ChunkDivideConfig%noChunks
        chunks(i)%firstMAFIndex = max ( (i-1)*maxLength - lowerOverlap, 0 )
        chunks(i)%lastMAFIndex = i*maxLength + upperOverlap - 1
        chunks(i)%noMAFsUpperOverlap = upperOverlap
        if ( i > 1 ) then
          chunks(i)%noMAFsLowerOverlap = lowerOverlap
        else
          chunks(i)%noMAFsLowerOverlap = 0
        end if
      end do
      if ( present(root) ) call trace_end ( &
        & 'ChunkDivide_fixed', cond=toggle(gen) .and. levels(gen) > 2 )

    end subroutine ChunkDivide_Fixed

    !------------------------------------------- ChunkDivide_Polygon -----
    subroutine ChunkDivide_Polygon ( mafRange, filedatabase, chunks, &
                                   & InstrumentModuleName, root )
      ! integer, dimension(2), intent(in) :: MAFRANGE
      use Allocate_Deallocate, only: Test_Allocate
      use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
      use Dump_Geolocation_m, only: Dump_H_t
      use Geometry, only: To_XYZ
      use PnPoly_m, only: PnPoly

      type (MAFRange_T) :: MAFRange
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS
      character (len=NameLen) :: InstrumentModuleName
      integer, intent(in), optional :: Root ! in l2cf tree, for tracing


      ! Local variables
      integer(c_intptr_t) :: Addr         ! For tracing
      integer :: I                        ! Loop inductor
      integer :: MAF
      integer :: Me = -1                  ! String index for trace
      integer :: STATUS                   ! From allocate
      integer :: MAXLENGTH                ! nint(config%maxLength)
      integer :: LOWEROVERLAP             ! nint(config%lowerOverlap)
      integer :: UPPEROVERLAP             ! nint(config%upperOverlap)
      type (L1BData_T) :: lats            ! lats data
      type (L1BData_T) :: lons            ! lons data
      type (MLSFile_T), pointer             :: L1BFile
      integer :: L1BFLAG                  ! error Flag
      integer :: NoMAFs
      character (len=NameLen) :: L1BItemName
      logical, dimension(:), pointer :: ptsInPolygon => null()
      real(r8), dimension(3)         :: xyz
      real(r8), dimension(:), allocatable :: Vertices_XX, Vertices_YY
      real(r8), dimension(3)         :: PolyVert_XYZ
      integer :: n
      integer, dimension(2) :: range
      ! integer :: pointTest

      ! Executable code
      if ( present(root) ) call trace_begin ( me, 'ChunkDivide_polygon', root, &
          & cond=toggle(gen) .and. levels(gen) > 2 )
      ! 1st-- read lats and lons from the l1b file
      L1BFile => GetMLSFileByType( filedatabase, content='l1boa' )
      if ( .not. associated(L1BFile) ) &
        & call MLSMessage  ( MLSMSG_Error, ModuleName, &
          & "Can't make progress in ChunkDivide_polygon without L1BOA files" )
      l1bItemName = AssembleL1BQtyName (  instrumentModuleName//".tpGeodLat", L1BFile%HDFVersion, .false. )
      call output (trim(l1bItemName), advance ='yes')
      call ReadL1BData ( L1BFile, l1bItemName, lats, noMAFs, &
        & l1bFlag )
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Where is ' // l1bItemName )
      call output ('lats is ', advance='yes' )
      call output( lats%dpField(1,1,1::NoMAFS), advance='yes')
      call dump ( lats, 0)
      l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpLon", L1BFile%HDFVersion, .false. )
      call output (trim(l1bItemName), advance ='yes')
      call ReadL1BData ( L1BFile, l1bItemName, lons, noMAFs, &
        & l1bFlag )
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Where is ' // l1bItemName )
      call output ('lons is ' )
      call dump ( lons, 0 )
      call allocate_test( ptsInPolygon, noMAFs, 'ptsInPolygon', &
        & ModuleName )
      n = size(polygon_vertices)
      allocate(Vertices_XX(n))
      allocate(Vertices_YY(n))



      do i = 1, n
        PolyVert_XYZ = to_xyz(polygon_vertices(i)%lat, polygon_vertices(i)%lon%d)
        call output ( polygon_vertices(i)%lon%d )
        Vertices_XX(i)= PolyVert_XYZ(1)
        call output ('   ')
        call output ( polygon_vertices(i)%lat, advance='yes' )
        Vertices_YY(i) = PolyVert_XYZ(2)
        call output ('Vertices_XX,  Vertices_YY for polygon is ')
        call output (Vertices_XX(i))
        call output ('   ')
        call output (Vertices_YY(i), advance='yes')
      end do

      ! Now we'll check for whether the points are inside or outside the PolyGon
      call output ('before do MAF loop', advance='yes')
      !pointTest = 0
      do MAF=1, NoMAFs
         call output('lat = ')
         call output( lats%dpField(1,1,MAF) )
         call output ('  lon = ')
         call output( lons%dpField(1,1,MAF), advance='yes')
         xyz = to_xyz ( lats%dpField(1,1,MAF), lons%dpField(1,1,MAF))
         call output ('xyz is ')
         call output (xyz, advance='yes')
         call output(PnPoly(xyz(1), xyz(2), Vertices_XX, Vertices_YY ), advance='yes')
         if ((PnPoly(xyz(1), xyz(2), Vertices_XX, Vertices_YY )) == 1)  then
             ptsInPolygon(MAF) = .true.
             !pointTest = pointTest + 1
           !  if (MAF == 70) then
           !    ptsInPolygon(MAF) = .false.
           !  end if
           !  if (MAF == 100) then
            !    ptsInPolygon(MAF) = .false.
           !  end if
           !  if (MAF == 200) then
           !     ptsInPolygon(MAF) = .false.
           !  end if
           !  if (MAF == 240) then
           !     ptsInPolygon(MAF) = .false.
           !  end if
           !  if (MAF == 300) then
           !     ptsInPolygon(MAF) = .false.
          !   end if
         else
             !ptsInPolygon(MAF) = pointTest
             !pointTest = 0
             ptsInPolygon(MAF) = .false.
         end if

        ! Use To_XYZ to get XX, YY, ZZ
        ! The send XX, YY through Pn_Poly
        ! Set to true or false depending on retrun value (-1 outside, +1 inside)
        ! store in ptsInPolygon
      enddo
      ! Now dump ptsInPolygon
      call dump(ptsInPolygon, "ptsInPolygon")

      call output ( 'ChunckDivide_Polygon is usign Fixed method exactly.\n' )
      call output ( 'vertices \n')
      call dump_h_t ( polygon_vertices, &
      & 'Polygon boundary vertices (lon,lat):', format='(f8.3)' )

      call output ( 'Polygon_Inside \n')
      call output ( polygon_inside%lon%d, before='Point inside polygon (lon,lat): (' )
      call output ( polygon_inside%lat, before=',', after=')', advance='yes' )

      call output ( 'mafRange ', advance='yes')
       call output ( 'Corresponding MAF range ')
        call output ( mafRange%L2Cover(1) )
        call output ( ' : ' )
        call output ( mafRange%L2Cover(2), advance='yes' )
        call output ( 'MAF range in L1B files ' )
        call output ( mafRange%L1BCover(1) )
        call output ( ' : ' )
        call output ( mafRange%L1BCover(2), advance='yes' )
        call output ( 'including overlapped days ' )
        call output ( mafRange%Expanded(1) )
        call output ( ' : ' )
        call output ( mafRange%Expanded(2), advance='yes' )

        call output ( 'before  FindLongestRange  call ', advance='yes')
        call FindLongestRange (ptsInPolygon, .true., range)
        call output ('range is ')
        call output (range, advance='yes')

        call output ( 'after  FindLongestRange  call ')

      swlevel = switchDetail(switches, 'chu' )
      allocate ( chunks(ChunkDivideConfig%noChunks), stat=status )
      addr = 0
      if ( status == 0 .and. ChunkDivideConfig%noChunks > 0 ) &
        & addr = transfer(c_loc(chunks(1)), addr)
      call test_allocate ( status, moduleName, 'chunks', &
        & uBounds = ChunkDivideConfig%noChunks, &
        & elementSize = storage_size(chunks) / 8, address=addr )

      maxLength = nint ( ChunkDivideConfig%maxLength )
      lowerOverlap = nint ( ChunkDivideConfig%lowerOverlap )
      upperOverlap = nint ( ChunkDivideConfig%upperOverlap )
      do i = 1, ChunkDivideConfig%noChunks
        chunks(i)%firstMAFIndex = max ( (i-1)*maxLength - lowerOverlap, 0 )
        chunks(i)%lastMAFIndex = i*maxLength + upperOverlap - 1
        chunks(i)%noMAFsUpperOverlap = upperOverlap
        if ( i > 1 ) then
          chunks(i)%noMAFsLowerOverlap = lowerOverlap
        else
          chunks(i)%noMAFsLowerOverlap = 0
        end if
      end do
      if ( present(root) ) call trace_end ( &
        & 'ChunkDivide_Polygon', cond=toggle(gen) .and. levels(gen) > 2 )

    end subroutine ChunkDivide_Polygon

    !---------------------------------------------- ChunkDivide_PE -----
    subroutine ChunkDivide_PE ( mafRange, filedatabase, chunks, root )
      ! type (ChunkDivideConfig_T), intent(in) :: CONFIG
      use Allocate_Deallocate, only: Test_Allocate
      use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
      integer, dimension(2), intent(in) :: MAFRANGE
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS
      integer, intent(in), optional :: Root ! in l2cf tree, for tracing

      ! Local parameters
      real(r8), parameter :: HOMEACCURACY = 3.0 ! Try to hit homeGeodAngle within this

      ! Local variables
      type (L1BData_T) :: TAITIME         ! From L1BOA
      type (L1BData_T) :: TPGEODANGLE     ! From L1BOA

      character(len=10) :: MODNAMESTR     ! Home module name as string

      integer(c_intptr_t) :: Addr         ! For tracing
      integer :: FLAG                     ! From ReadL1B
      integer :: HOMEMAF                  ! first MAF after homeGeodAngle
      integer :: HOME                     ! Index of home MAF in array
      integer :: Me = -1                  ! String index for trace
      integer :: M1, M2                   ! MafRange + 1
      integer :: NOMAFS                   ! Number of MAFs to consider
      integer :: NOMAFSREAD               ! From ReadL1B
      integer :: ORBIT                    ! Used to locate homeMAF
      integer :: STATUS                   ! From allocate etc.
      integer :: NOCHUNKS                 ! Number of chunks
      integer :: MAXLENGTH                ! Max length as integer (MAFs)

      real(r8) :: ANGLEINCREMENT          ! Increment in hunt for homeMAF
      real(r8) :: MAXANGLE                ! Of range in data
      real(r8) :: MINANGLE                ! Of range in data
      real(r8) :: TESTANGLE               ! Angle to check for

      integer   ::                 l1b_hdf_version
      character(len=namelen) ::    MAF_start, tp_angle
      type (MLSFile_T), pointer :: L1BFile

      ! Executable code
      if ( present(root) ) call trace_begin ( me, 'ChunkDivide_PE', root, &
          & cond=toggle(gen) .and. levels(gen) > 2 )
      swlevel = switchDetail(switches, 'chu' )

      ! Read in the data we're going to need
      call get_string ( lit_indices(ChunkDivideConfig%homeModule), modNameStr, &
        & strip=.true. )
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      if ( .not. associated(L1BFile) ) &
        & call MLSMessage  ( MLSMSG_Error, ModuleName, &
          & "Can't make progress in ChunkDivide_PE without L1BOA files" )
      ! call dump(L1BFile)
      l1b_hdf_version = L1BFile%HDFVersion
      MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
        .false. )
      tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
        & l1b_hdf_version, &
        & .false. )
      call ReadL1BData ( L1BFile, trim(tp_angle), &
        & tpGeodAngle, noMAFsRead, flag, &
        & dontPad=DONTPAD )
      ! If you needed this, should you be using this chunkDivide?
      call smoothOutDroppedMAFs(tpGeodAngle%dpField)
      call smoothOutDroppedMAFs(tpGeodAngle%dpField, monotonize=.true.)
      call ReadL1BData ( L1BFile, trim(MAF_start), &
        & taiTime, noMAFsRead, flag, &
        & dontPad=DONTPAD )
      call smoothOutDroppedMAFs(taiTime%dpField)
      noMAFs = mafRange(2) - mafRange(1) + 1
      m1 = mafRange(1) + 1
      m2 = mafRange(2) + 1

      minAngle = minval ( tpGeodAngle%dpField(1,1,m1:m2) )
      maxAngle = maxval ( tpGeodAngle%dpField(1,1,m1:m2) )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! First try to locate the last MAF before the homeGeodAngle
      orbit = int ( tpGeodAngle%dpField(1,1,m1)/360.0 )
      if ( tpGeodAngle%dpField(1,1,m1) < 0.0 ) orbit = orbit - 1
      testAngle = ChunkDivideConfig%homeGeodAngle + orbit*360.0
      if ( ChunkDivideConfig%maxLengthFamily == PHYQ_Angle ) then
        angleIncrement = ChunkDivideConfig%maxLength
      else
        angleIncrement = 360.0
      end if

      maxLength = nint ( ChunkDivideConfig%maxLength )

      homeHuntLoop: do
        if ( testAngle < minAngle ) then
          testAngle = testAngle + angleIncrement
          cycle
        endif
        if ( testAngle > maxAngle ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Unable to establish a home major frame, using the first' )
          home = 1
          exit homeHuntLoop
        end if
        ! Find MAF which starts before this test angle
        call Hunt ( tpGeodAngle%dpField(1,1,:), testAngle, home, nearest=.true.,&
          & allowTopValue = .true. )
        ! Now if this is close enough, and has enough MAFs around it, accept it
        if ( ( abs ( tpGeodAngle%dpField(1,1,home) - &
          & testAngle ) < HomeAccuracy ) &
          & .and. ( home - maxLength/2 >= 1 ) &
          & .and. ( home + maxLength/2 <= noMAFs ) ) exit homeHuntLoop
        ! Otherwise, keep looking
        testAngle = testAngle + angleIncrement
      end do homeHuntLoop
      homeMAF = home + m1 - 1

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OK, now we have a home MAF, get a first cut for the chunks
      ! We work out the chunk ends for each chunk according to how the
      ! maxLength field is specified.
      ! ChunkDivideConfig%maxLengthFamily == PHYQ_MAFs
      noChunks = 1

      ! Allocate the chunk
      allocate ( chunks(noChunks), stat=status )
      addr = 0
      if ( status == 0 .and. noChunks > 0 ) addr = transfer(c_loc(chunks(1)), addr )
      call test_allocate ( status, ModuleName, 'chunks', uBounds = noChunks, &
        & elementSize = storage_size(chunks) / 8, address=addr )

      ! Work out its position (in the file)
      chunks(1)%firstMAFIndex = homeMAF - maxLength/2 - 1
      chunks(1)%lastMAFIndex = homeMAF + maxLength/2 - 1

      chunks(1)%noMAFsUpperOverlap = 0
      chunks(1)%noMAFsLowerOverlap = 0

      call DeallocateL1BData ( tpGeodAngle )
      call DeallocateL1BData ( taiTime )
      if ( present(root) ) call trace_end ( &
        & 'ChunkDivide_PE', cond=toggle(gen) .and. levels(gen) > 2 )

    end subroutine ChunkDivide_PE

    !--------------------------------------------- ChunkDivideL2CF -----
    subroutine ChunkDivideL2CF ( root )
      ! This subroutine identifies, separates, and checks values from the section
      ! of the MLSCF (ChunkDivide) passed to Scan/Divide.
      use HighOutput, only: LetsDebug, OutputNamedValue
      use Init_Tables_Module, only: Field_Indices
      use moreTree, only: get_boolean
      integer, intent(in) :: ROOT    ! Root of the ChunkDivide command
                                     ! MLSCF abstract syntax tree
      ! type (ChunkDivideConfig_T), intent(out) :: CONFIG ! Result of operation

      integer, target, dimension(2) :: NeededForFixed = &
        & (/ f_maxLength, f_overlap /)
        ! & (/ f_noChunks, f_maxLength, f_overlap /)
      integer, target, dimension(7) :: NotWantedForFixed = &
        & (/ f_noSlaves, f_homeModule, f_homeGeodAngle, f_scanLowerLimit, &
        &    f_scanUpperLimit, f_criticalModules, f_maxGap /)

      integer, target, dimension(3) :: NeededForPE = &
        & (/ f_maxLength, f_homeModule, f_homeGeodAngle /)
      integer, target, dimension(7) :: NotWantedForPE = &
        & (/ f_noChunks, f_overlap, f_noSlaves, f_scanLowerLimit, &
        &    f_scanUpperLimit, f_criticalModules, f_maxGap /)

      integer, target, dimension(5) :: NeededForOrbital = &
        & (/ f_maxLength, f_overlap, f_homeModule, f_homeGeodAngle, &
        &    f_maxGap /)
        ! &    f_criticalModules, f_maxGap /)
      integer, target, dimension(1) :: NotWantedForOrbital = &
        & (/ f_noSlaves /)

      integer, target, dimension(6) :: NeededForEven = &
        & (/ f_maxLength, f_overlap, f_maxLength, f_noSlaves, f_maxGap, &
        &    f_criticalModules /)
      integer, target, dimension(3) :: NotWantedForEven = &
        & (/ f_noChunks, f_homeModule, f_homeGeodAngle /)

      ! Local variables
      logical :: debug
      integer :: FieldIndex               ! Of son of the ChunkDivide command
      integer :: FieldValue
      integer :: GSON                     ! A son of son the ChunkDivide command node
      integer :: I                        ! Loop inductor
      integer :: J                        ! Another loop inductor
      logical :: Log_Value                ! For boolean fields
      integer :: Me = -1                  ! String index for trace cacheing
      integer :: numCriticalSignals       ! How many crit. Sigs
      integer :: signalsNode              ! Node where crit. Sig. begin
      integer :: SON                      ! A son of the ChunkDivide command
      integer :: UNITS(2)                 ! Units of expression
      integer, dimension(:), pointer :: NEEDED ! Which fields are needed
      integer, dimension(:), pointer :: NOTWANTED ! Which fields are not wanted
      logical :: GOT(field_first:field_last) = .false.
      real(rp) :: VALUE(2)                ! Value of expression

      ! Executable code
      call trace_begin ( me, 'ChunkDivideL2CF', root, &
        & cond=toggle(gen) .and. levels(gen) > 0 )
      swlevel = switchDetail(switches, 'chu' )
      debug = LetsDebug ( 'chu', 2 )
      error = 0
      ChunkDivideConfig%where = root

      got = .false.

      ! Loop through the command identifying fields and storing values.
      do i = 2, nsons(root) ! Skip the command
        son = subtree(i,root)
        if ( nsons(son) > 1 ) then
          fieldValue = decoration(subtree(2,son)) ! The field's value
        else
          fieldValue = son
        end if
        fieldIndex = get_field_id(son)
        got(fieldIndex) = .true.
        call trace_begin ( 'ChunkDivideL2CF.loop', root, &
          & string=field_indices(fieldIndex),&
          & cond=toggle(gen) .and. levels(gen) > 1 )
        if ( &
          & any(fieldIndex == &
          & (/ f_excludePostOverlaps, f_excludePriorOverlaps, f_skipL1BCheck, &
          & f_crashIfPhiNotMono, f_saveObstructions /) &
          & ) &
          & ) then
          log_value = get_Boolean(son)
        else
          log_value = .false.
        endif
        if ( debug ) then
        call OutputNamedValue ( 'fieldIndex ', fieldIndex )
        call OutputNamedValue ( 'nsons(son) ', nsons(son) )
        call OutputNamedValue ( 'log_value ', log_value )
        endif
        if (nsons(son) > 1 ) then
          gson = subtree(2,son)
          call expr ( gson, units, value )
          ! log_value = get_boolean ( fieldValue )
        elseif (nsons(son) > 0 ) then
          value = 0.0
          ! log_value = get_boolean ( fieldValue )
        else
          value = 0.0
          ! log_value = .false.
        end if
        ! Get value for this field if appropriate
        select case ( fieldIndex )
        case ( f_method )
          ChunkDivideConfig%method = decoration ( gson )
        case ( f_module )
          call output ('Inside case f_module', advance='yes')
          !ChunkDivideConfig%module = decoration ( gson )
          instrumentModule = decoration(decoration(subtree(2,son)))
          call GetModuleName ( instrumentModule , instrumentModuleName )
          call output ('insturmentModule name is   ')
          call output (instrumentModule , advance='yes')
          call output ('insturmentModuleNAME name is   ')
          call output (instrumentModuleName , advance='yes')
          ChunkDivideConfig%module = instrumentModule
        case ( f_noChunks )
          ChunkDivideConfig%noChunks = nint ( value(1) )
        case ( f_maxLength )
          ChunkDivideConfig%maxLength = value(1)
          ChunkDivideConfig%maxLengthFamily = units(1)
        case ( f_maxOrbY )
          ChunkDivideConfig%maxOrbY = value(1)
          if ( units(1) /= PHYQ_Length ) &
            & call AnnounceError ( root, BadUnits, fieldIndex )
        case ( f_overlap )
          ChunkDivideConfig%overlap = value(1)
          ChunkDivideConfig%overlapFamily = units(1)
        case ( f_lowerOverlap )
          ChunkDivideConfig%lowerOverlap = value(1)
          ChunkDivideConfig%lowerOverlapFamily = units(1)
        case ( f_upperOverlap )
          ChunkDivideConfig%upperOverlap = value(1)
          ChunkDivideConfig%upperOverlapFamily = units(1)
        case ( f_noSlaves )
          ChunkDivideConfig%noSlaves = value(1)
        case ( f_homeModule )
          ChunkDivideConfig%homeModule = decoration ( gson )
          !instrumentModule = decoration(decoration(subtree(2,son)))
          !call GetModuleName ( instrumentModule , instrumentModuleName )
          !ChunkDivideConfig%homeModule = instrumentModule
        case ( f_homeGeodAngle )
          ChunkDivideConfig%homeGeodAngle = value(1)
        case ( f_scanLowerLimit )
          if ( any ( units /= PHYQ_Dimensionless .and. units /= PHYQ_Length ) &
            & .or. .not. any ( units == PHYQ_Length ) ) &
            & call AnnounceError ( root, BadUnits, fieldIndex )
          ChunkDivideConfig%scanLowerLimit = value
          ChunkDivideConfig%scanLLSet = .true.
        case ( f_scanUpperLimit )
          if ( any ( units /= PHYQ_Dimensionless .and. units /= PHYQ_Length ) &
            & .or. .not. any ( units == PHYQ_Length ) ) &
            & call AnnounceError ( root, BadUnits, fieldIndex )
          ChunkDivideConfig%scanUpperLimit = value
          ChunkDivideConfig%scanULSet = .true.
        case ( f_criticalBands )
          call get_string ( sub_rosa ( gson ), &
              & ChunkDivideConfig%criticalBands, strip=.true. )
          call output ( 'reading critical bands ' )
          call output ( trim(ChunkDivideConfig%criticalBands), advance='yes' )
        case ( f_criticalModules )
          ChunkDivideConfig%criticalModules = decoration ( gson )
        case ( f_criticalSignals )
          signalsNode = son
          numCriticalSignals = nsons(signalsNode) - 1
          call allocate_test(ChunkDivideConfig%criticalSignals, numCriticalSignals, &
            & 'critical signals', ModuleName)
          do j = 2, nsons(signalsNode)    ! Skip name
            call get_string ( sub_rosa ( subtree (j, signalsNode) ), &
              & ChunkDivideConfig%criticalSignals(j-1), strip=.true. )
          end do
        case ( f_excludePostOverlaps )
          ChunkDivideConfig%allowPostOverlaps = .not. log_value ! log_value
          if ( .not. ChunkDivideConfig%allowPostOverlaps ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to exclude MAFs after time range' )
        case ( f_excludePriorOverlaps )
          ChunkDivideConfig%allowPriorOverlaps = .not. log_value ! log_value
          if ( .not. ChunkDivideConfig%allowPriorOverlaps ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to exclude MAFs prior to time range' )
        case ( f_maxGap )
          ChunkDivideConfig%maxGap = value(1)
          ChunkDivideConfig%maxGapFamily = units(1)
        case ( f_skipL1BCheck )
          ! print *, 'processing f_skipL1BCheck ', log_value
          ChunkDivideConfig%skipL1BCheck = log_value ! log_value
          if ( ChunkDivideConfig%skipL1BCheck ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to skip checking the l1b data for problems' )
        case ( f_crashIfPhiNotMono )
          ChunkDivideConfig%crashIfPhiNotMono = log_value ! log_value
          if ( ChunkDivideConfig%crashIfPhiNotMono ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to crash if phi, the master geodetic angle, is' // &
            & ' not monotonic' )
        case ( f_saveObstructions )
          ChunkDivideConfig%saveObstructions = log_value ! log_value
          if ( ChunkDivideConfig%saveObstructions ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to save obstructions (possibly for Output_Close)' )
        end select
        call trace_end ( 'ChunkDivideL2CF.loop', &
          & cond=toggle(gen) .and. levels(gen) > 1 )
      end do

      ! Now check the sanity of what we've been given, this varies a bit
      ! depending on the method
      if ( ChunkDivideConfig%skipL1BCheck .and. &
        & ChunkDivideConfig%crashIfPhiNotMono ) &
        & call MLSMessage(MLSMSG_Error, ModuleName, &
            & 'How can we check if phi, the master geodetic angle, is' // &
            & ' not monotonic when we skip checking the l1B file?' )
      select case ( ChunkDivideConfig%method )
      case ( l_fixed )
        needed => NeededForFixed
        notWanted => NotWantedForFixed
        ! We made NoChunks an optional parameter instead of a required one
        ! So if it's not present, the default value, 0, will be
        ! replaced with the number of MAFs in the level 1 files
        ! If that's not what you want, then you had better overide the default
        ! 0 with the actual number
        if ( ChunkDivideConfig%NoChunks < 1 ) &
          & ChunkDivideConfig%NoChunks = 1 + &
          & (GlobalAttributes%LastMAFCtr - GlobalAttributes%FirstMAFCtr) / &
          & ChunkDivideConfig%maxLength
      case ( l_polygon )
        needed => NeededForFixed
        notWanted => NotWantedForFixed
      case ( l_PE )
        needed => NeededForPE
        notWanted => NotWantedForPE
      case ( l_orbital )
        needed => NeededForOrbital
        notWanted => NotWantedForOrbital
      case ( l_even )
        needed => NeededForEven
        notWanted => NotWantedForEven
      end select

      ! Some special thought for overlap cases
      if ( any ( got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) ) then
        if ( .not. all ( got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) ) &
          & call AnnounceError ( root, notSpecified, &
          & FindFirst ( .not. got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) )
        if ( got ( f_overlap ) ) &
          & call AnnounceError ( root, unnecessary, f_overlap )
        got ( f_overlap ) = .true.
        ! Insist that upper/lower overlaps be specified in the same units
        if ( ChunkDivideConfig%lowerOverlapFamily /= ChunkDivideConfig%upperOverlapFamily ) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      else
        ! If single overlap specified, copy it into upper/lower
        if ( got ( f_overlap ) ) then
          ChunkDivideConfig%lowerOverlap = ChunkDivideConfig%overlap
          ChunkDivideConfig%lowerOverlapFamily = ChunkDivideConfig%overlapFamily
          ChunkDivideConfig%upperOverlap = ChunkDivideConfig%overlap
          ChunkDivideConfig%upperOverlapFamily = ChunkDivideConfig%overlapFamily
        end if
      end if

      ! Check we've got all the arguments we need
      do i = 1, size(needed)
        if ( .not. got(needed(i) ) ) &
          & call AnnounceError ( root, notSpecified, needed(i) )
      end do
      ! Check we don't have unnecessary ones
      do i = 1, size(notWanted)
        if ( got(notWanted(i) ) ) &
          & call AnnounceError ( root, unnecessary, notWanted(i) )
      end do

      ! Make other checks of parameters
      if ( ChunkDivideConfig%criticalModules /= l_none ) then
        if ( .not. got(f_scanLowerLimit) ) &
          & call AnnounceError ( root, notSpecified, f_scanLowerLimit )
        if ( .not. got (f_scanUpperLimit) ) &
          & call AnnounceError ( root, notSpecified, f_scanUpperLimit )
      else
        if ( got(f_scanLowerLimit) ) &
          & call AnnounceError ( root, unnecessary, f_scanLowerLimit )
        if ( got (f_scanUpperLimit) ) &
          & call AnnounceError ( root, unnecessary, f_scanUpperLimit )
      end if

      ! Now check the units for various cases
      if ( got(f_maxgap) .and. all ( ChunkDivideConfig%maxGapFamily /= &
        & (/ PHYQ_MAFs, PHYQ_Angle, PHYQ_Time /))) &
        & call AnnounceError ( root, badUnits, f_maxGap )
      select case ( ChunkDivideConfig%method )
      case ( l_PE )
        if ( ChunkDivideConfig%maxLengthFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_maxLength )
      case ( l_orbital )
        if (all(ChunkDivideConfig%maxLengthFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_maxLength )
        if (all(ChunkDivideConfig%lowerOverlapFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_lowerOverlap )
        if (all(ChunkDivideConfig%upperOverlapFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      case ( l_fixed, l_even )
        if ( ChunkDivideConfig%maxLengthFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_maxLength )
        if ( ChunkDivideConfig%lowerOverlapFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_lowerOverlap )
        if ( ChunkDivideConfig%upperOverlapFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      end select

      if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Problem with ChunkDivide command' )

      ! That's it, we're all valid now.
      call trace_end ( 'ChunkDivideL2CF', &
        & cond=toggle(gen) .and. levels(gen) > 0 )

    end subroutine ChunkDivideL2CF

    ! ---------------------------------------------------- SayTime -----
    subroutine SayTime
      real :: T2
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( 'Timing for ChunkDivide = ' )
      call output ( dble(t2-t1), advance='yes' )
      timing = .false.
    end subroutine SayTime

    ! ---------------------------------------------------- CheckChunkSanity -----
    subroutine CheckChunkSanity ( chunks, mafRange )
      type (MLSChunk_T), dimension(:), intent(inout) :: CHUNKS
      integer, dimension(2), intent(in) :: MAFRANGE   ! Processing range in MAFs

      ! Local variables
      integer :: i
      logical :: sharing, outrange, negovlps

      ! Executable statements
      sharing = .false.
      outrange = .false.
      negovlps = .false.
      if ( CHECKFORSHAREDMAFS ) then
        do i = 1, size(chunks) - 1
          sharing = (chunks(i)%lastMAFIndex - chunks(i)%noMAFsUpperOverlap) >= &
            & chunks(i+1)%firstMAFIndex + chunks(i+1)%noMAFsLowerOverlap
          if ( sharing ) exit
        enddo
      endif
      outrange = any(chunks%firstMAFIndex < mafrange(1)) .or. &
        &        any(chunks%lastMAFIndex > mafrange(2))
      negovlps = any(chunks%noMAFsLowerOverlap < 0) .or. &
        &        any(chunks%noMAFsUpperOverlap < 0)
      if ( sharing ) then
        call output ( "Shared non-overlaps MAFs detected ", advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Shared non-overlapped MAFs' )
      elseif ( outrange ) then
        call output ( " MAFs detected outside range: " )
        call output ( mafRange, advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'MAFs outside range' )
      elseif ( negovlps ) then
        call output ( "Negative overlaps detected ", advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Negative overlaps' )
      endif
    end subroutine CheckChunkSanity
  end subroutine ChunkDivide

  ! -------------------------------------------- Dump Obstructions -----
  subroutine Dump_Obstructions ( obstructions )

    use Output_M, only: OUTPUT

    type (Obstruction_T), dimension(:), pointer :: obstructions
    ! Local variables
    integer :: i                        ! Loop counter
    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    if ( associated ( obstructions ) ) then
      if ( size(obstructions) == 0 ) then
        call output ( 'Obstructions is a zero size array.', advance='yes' )
      else
        call output ( 'Dumping ' )
        call output ( size(obstructions) )
        call output ( ' obstructions:', advance='yes' )
        do i = 1, size(obstructions)
          call output ( i )
          if ( obstructions(i)%range ) then
            call output ( ' : Range [ ' )
            call output ( obstructions(i)%mafs(1) )
            call output ( ' : ' )
            call output ( obstructions(i)%mafs(2) )
            call output ( ' ]', advance='yes' )
          else
            call output ( ' : Wall [ ' )
            call output ( obstructions(i)%mafs(1) )
            call output ( ' ]', advance='yes' )
          end if
        end do
      end if
    else
      call output ( 'Obstructions is not associated.', advance='yes')
    end if
  end subroutine Dump_Obstructions

! -----------------------------------------------  ANY_GOOD_SIGNALDATA  -----
  function ANY_GOOD_SIGNALDATA ( signal, sideband, filedatabase, maf, maf2, &
    & good_buffer, mafRange )  &
    & result (answer)
  ! Scalar use:
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  !
  ! Array use
  ! Return one logical value in good_buffer
  ! for each of the_mafs between maf and maf2
  ! Arguments

    use MLSKinds, only: Rk => R8
    use MLSSignals_m, only: GetSignalName

    integer, intent(in)                            :: signal
    integer, intent(in)                            :: sideband
    logical                                        :: answer
    integer, intent(in)                            :: maf
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer, optional, intent(in)                  :: maf2
    logical, dimension(:), optional, intent(inout) :: good_buffer
    integer, dimension(2), optional, intent(in)    :: MAFRANGE   ! Processing
    ! Private                                                   range in MAFs
    integer :: flag
    type(MLSFile_T), pointer :: L1BFile
    integer :: maf_index   ! 1 <= maf_index <= mafrange(2)-mafrange(1)+1
    type (l1bData_T) :: MY_L1BDATA
    integer :: mymaf2
    character(len=127)  :: namestring
    integer :: noMAFs
    integer :: the_maf

    ! Executable
    swlevel = switchDetail(switches, 'chu' )
    answer = .false.
    mymaf2 = maf
    if ( present(maf2) ) mymaf2 = maf2

    ! Set defaults for the good_buffer if present
    if ( present(good_buffer) ) then
      if ( .not. present(mafrange) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'mafRange must be supplied to any_good_signaldata' )
      endif
      good_buffer ( maf - mafRange(1) + 1 : myMAF2 - mafRange(1) + 1 ) = .false.
    end if

    ! OK, try to find this item in an l1brad file
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) return
    call GetSignalName ( signal, nameString, &
      & sideband=sideband, noChannels=.TRUE. )
    nameString = AssembleL1BQtyName ( nameString, L1BFile%hdfVersion, .false. )
    nameString = trim(nameString) // PRECISIONSUFFIX
    ! fileID = FindL1BData (filedatabase, nameString, hdfVersion )
    L1BFile => GetL1bFile(filedatabase, namestring)
    ! If not found, exit appropriately
    if ( .not. associated(L1BFile) ) then
      answer = .false.
      return
    end if

    ! OK, we've found it.  Read the data in.
    call ReadL1BData ( L1BFile, nameString, my_l1bData, noMAFs, flag, &
      & firstMAF=maf, lastMAF=mymaf2, &
      & NeverFail= .true., dontPad = .false. )
    ! Quit if the reading failed.
    if ( flag /= 0 ) then
      answer = .false.
      return
    end if

    ! Give detailed or curt response
    if ( present(good_buffer) ) then
      answer = .true.    ! This value should be ignored by the caller
      do the_maf = maf, mymaf2
        maf_index = the_maf - mafRange(1) + 1
        good_buffer(maf_index) = .not. &
          & all (my_l1bData%DpField(:,:, the_maf+1-maf) < 0.0_rk)
      end do
    else
      answer = .not. all (my_l1bData%DpField < 0._rk)
    end if
    call DeallocateL1bData ( my_l1bData )

  end function ANY_GOOD_SIGNALDATA

  logical function isNotACriticalBand( band )
    use MLSStringLists, only: NumStringElements, StringElement
    use MLSStrings, only: ReadIntsFromChars
    ! Args
    integer, intent(in) :: band
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: i
    integer :: n
    ! Executable
    isNotACriticalBand = .true.
    do i=1, NumStringElements( trim(ChunkDivideConfig%criticalBands), countEmpty )
      call readIntsFromChars( trim(&
        & StringElement( ChunkDivideConfig%criticalBands, i, countEmpty ) &
        & ), n, ignore='*' )
      isNotACriticalBand = isNotACriticalBand .and. (band /= n)
    enddo
  end function isNotACriticalBand

  subroutine smoothOutDroppedMAFs(field, wasSmoothed, monotonize)
    ! detect any fillValues--replace them with nearest neighbor values
    ! or, optionally, detect and correct any departures from monotone growth
    use MLSFillValues, only: IsFillValue
    ! Args
    real(r8), intent(inout)                      :: field(:,:,:)
    logical, dimension(:), optional, intent(out) :: wasSmoothed
    logical, optional, intent(in)                :: monotonize
    ! Internal variables
    integer :: maf, nearest
    logical :: myMonotonize
    real(r8):: lastValue
    ! Executable
    swlevel = switchDetail(switches, 'chu' )
    myMonotonize = .false.
    if ( present(monotonize) ) mymonotonize = monotonize
    if ( present(wasSmoothed) ) wasSmoothed = .false.
    lastValue = field(1,1,1)
    do maf=1, size(field, 3)
      if ( myMonotonize ) then
        nearest = max(maf-1, 1)
        if ( field(1,1,maf) < lastValue ) then
          if ( present(wasSmoothed) ) wasSmoothed(maf) = .true.
          field(:,:,maf) = field(:,:,nearest)
        else
          lastValue = field(1,1,maf)
        endif
      elseif ( any( isFillValue(field(:,:,maf)) ) ) then
        if ( present(wasSmoothed) ) wasSmoothed(maf) = .true.
        if ( maf == 1 ) then
          nearest = findfirst(.not. isFillValue(field(1,1,:)) )
        else
          nearest = maf - 1
        endif
        field(:,:,maf) = field(:,:,nearest)
      endif
    enddo
  end subroutine smoothOutDroppedMAFs

  ! ----------------------------------------- DealWithObstructions -----
  subroutine DealWithObstructions ( chunks, obstructions )
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    ! This routine modifies the chunks according to the information
    ! given in the obstructions.

    ! Local variables
    type (MLSChunk_T) :: NEWCHUNK       ! A chunk to create
    integer :: CHUNK                    ! Index of chunk under some consideration
    integer :: DEADCHUNK                ! Index of chunk to kill
    integer :: FIRSTMAF                 ! Index of first MAF in range
    integer :: LASTMAF                  ! Index of last MAF in range
    integer :: MAF                      ! Index of MAF for wall
    integer :: OBSTRUCTION              ! Loop counter

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )

    do obstruction = 1, size ( obstructions )
      ! Find chunk where this obstruction is/starts
      if ( obstructions(obstruction)%range ) then
        ! A range obstruction

        ! first identify the chunks whose non overlapped portion falls completely
        ! within the range and delete them
        firstMAF = obstructions(obstruction)%mafs(1)
        lastMAF = obstructions(obstruction)%mafs(2)
        insideRange: do
          deadChunk = FindFirst ( &
            & chunks%firstMAFIndex + chunks%noMAFsLowerOverlap >= firstMAF .and. &
            & chunks%lastMAFIndex - chunks%noMAFsUpperOverlap <= lastMAF )
          if ( deadChunk == 0 ) exit insideRange
          call DeleteChunk ( chunks, deadChunk )
        end do insideRange

        ! Now think about chunks whose non overlapped part completely
        ! encompass the range, and split them.
        ! (we'll deal with the issue of ranges spilling into overlap regions below)
        chunk = FindFirst ( &
          & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap <= firstMAF ) .and. &
          & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap >= lastMAF ) )
        if ( chunk /= 0 ) then
          ! Create two new chunks with the wall between, first the lower portion
          newChunk = chunks ( chunk )
          newChunk%lastMAFIndex = firstMAF - 1
          newChunk%noMAFsUpperOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now the upper portion
          newChunk = chunks ( chunk )
          newChunk%firstMAFIndex = lastMAF + 1
          newChunk%noMAFsLowerOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Delete the old chunk
          call DeleteChunk ( chunks, chunk )
        endif

        ! So the cases we have left do not involve deleting or creating new
        ! chunks, just modifying existing ones.  Note that this includes
        ! cases where a range starts end ends in a chunk but one or other boundary
        ! is in the overlap regions.

        ! Look for chunks where the range starts in the chunk
        moveEnd: do
          chunk = FindFirst ( &
            & chunks%firstMAFIndex <= firstMAF .and. &
            & chunks%lastMAFIndex >= firstMAF )
          if ( chunk == 0 ) exit moveEnd
          chunks(chunk)%noMAFsUpperOverlap = max ( 0, &
            & chunks(chunk)%noMAFsUpperOverlap - &
            &    ( chunks(chunk)%lastMAFIndex - (firstMAF-1) ) )
          chunks(chunk)%lastMAFIndex = firstMAF - 1
        end do moveEnd

        ! Look for chunks where the range ends in the chunk
        moveStart: do
          chunk = FindFirst ( &
            & chunks%firstMAFIndex <= lastMAF .and. &
            & chunks%lastMAFIndex >= lastMAF )
          if ( chunk == 0 ) exit moveStart
          chunks(chunk)%noMAFsLowerOverlap = max ( 0, &
            & chunks(chunk)%noMAFsLowerOverlap - &
            &    ( (lastMAF+1) - chunks(chunk)%firstMafIndex ) )
          chunks(chunk)%firstMAFIndex = lastMAF + 1
        end do moveStart

      else ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! A wall obstruction
        maf = obstructions(obstruction)%mafs(1)

        ! First look for where the wall occurs in the non overlapped part
        chunk = FindFirst ( &
          & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap <= maf ) .and. &
          & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap >= maf ) )

        if ( chunk /= 0 ) then
          ! Create two new chunks with the wall between, first the lower portion
          newChunk = chunks ( chunk )
          newChunk%lastMAFIndex = maf - 1
          newChunk%noMAFsUpperOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now the upper portion
          newChunk = chunks ( chunk )
          newChunk%firstMAFIndex = maf
          newChunk%noMAFsLowerOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now delete the original chunk
          call DeleteChunk ( chunks, chunk )
        endif

        ! For chunks where the wall is in the overlap, just make the
        ! overlap shorter
        ! First the lower overlap
        wallInLower: do
          chunk = FindFirst ( &
            & ( chunks%firstMAFIndex < maf ) .and. &
            & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap > maf ) )
          if ( chunk == 0 ) exit wallInLower
          chunks(chunk)%noMAFsLowerOverlap = &
          & max ( chunks(chunk)%noMAFsLowerOverlap - &
          & ( maf - chunks(chunk)%firstMAFIndex ), 0 )
          chunks(chunk)%firstMAFIndex = maf
        end do wallInLower
        ! Now the upper overlap
        wallInUpper: do
          chunk = FindFirst ( &
            & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap < maf ) .and. &
            & ( chunks%lastMAFIndex > maf ) )
          if ( chunk == 0 ) exit wallInUpper
          chunks(chunk)%noMAFsUpperOverlap = &
            & max ( chunks(chunk)%noMAFsUpperOverlap - &
            & ( chunks(chunk)%lastMAFIndex - (maf-1) ), 0 )
          chunks(chunk)%lastMAFIndex = maf - 1
        end do wallInUpper
      end if  ! - - - - - - - - - - - - ! Wall obstructions

    end do                              ! Loop over obstructions

    ! Sort the chunks back into order
    call SortChunks ( chunks )

  end subroutine DealWithObstructions

  ! --------------------------------------------- SortObstructions -----
  subroutine SortObstructions ( obstructions )
    ! Sort the obstructions into order of increasing
    ! mafs(1) (start/wall MAF index)
    type (Obstruction_T), dimension(:), intent(inout) :: OBSTRUCTIONS

    ! Local variables
    type (Obstruction_T) :: TEMP
    integer :: I                        ! Loop counters
    integer, dimension(1) :: TOSWAP     ! Index

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    do i = 1, size(obstructions) - 1
      toSwap = minloc ( obstructions(i:)%mafs(1) ) + (/ i-1 /)
      if ( toSwap(1) /= i ) then
        temp = obstructions(i)
        obstructions(i) = obstructions(toSwap(1))
        obstructions(toSwap(1)) = temp
      end if
    end do
  end subroutine SortObstructions

  ! -------------------------------------------- DeleteObstruction -----
  subroutine DeleteObstruction ( obstructions, index )
    ! Dummy arguments
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (Obstruction_T), pointer, dimension(:) :: OBSTRUCTIONS
    integer, intent(in) :: INDEX

    ! Local variables
    type (Obstruction_T), pointer, dimension(:) :: TEMP
    integer(c_intptr_t) :: Addr    ! For tracing
    integer :: S                   ! Size in bytes of object to deallocate
    integer :: STATUS              ! From allocate

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    allocate ( temp ( size(obstructions) - 1 ), stat=status )
    addr = 0
    if ( status == 0 .and. size(obstructions) > 0 ) addr = transfer(c_loc(temp(1)), addr )
    call test_allocate ( status, ModuleName, 'temp', &
      & uBounds = size(obstructions) - 1, &
      & elementSize = storage_size(temp) / 8, address=addr )

    if ( index > 1 ) temp(1:index-1) = obstructions(1:index-1)
    if ( index < size(obstructions) .and. size(obstructions) > 1 ) &
      & temp(index:) = obstructions(index+1:)

    s = size(obstructions) * storage_size(obstructions) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(obstructions(1)), addr )
    deallocate ( obstructions, stat=status )
    call test_deallocate ( status, moduleName, 'obstructions', s, address=addr )

    obstructions => temp

  end subroutine DeleteObstruction

  ! ---------------------------------------- ChooseCriticalSignals -----
  subroutine ChooseCriticalSignals ( criticalSignals )
    use MLSSignals_m, only: Dump_Signals, GetSignalName, Signals
    use MLSStringLists, only: CatLists, List2Array
    use MLSStrings, only: Lowercase

    ! Args
    character(len=160), dimension(:), pointer :: criticalSignals
    ! Internal variables
    character(len=40) :: critical_module_str
    integer, parameter :: MAXCRITICALSIGNALLENGTH = 40000
    character(len=MAXCRITICALSIGNALLENGTH) :: criticalSignalStr
    character(len=40) :: module_str
    integer :: numCriticalSignals
    character(len=40) :: signal_full
    integer :: signalIndex
    logical, dimension(:), pointer :: isSignalcritical => null()
    ! Executable
    swlevel = switchDetail(switches, 'chu' )
    if ( ChunkDivideConfig%criticalModules == l_none .and. &
      & len_trim(ChunkDivideConfig%criticalBands) < 1 ) return
    critical_module_str = ' '
    if ( ChunkDivideConfig%criticalModules /= l_none ) then
      call get_string( lit_indices(ChunkDivideConfig%criticalModules), signal_full, &
        & strip=.true. )
      critical_module_str = lowercase(signal_full)
      if ( critical_module_str /= 'ghz' ) return
    endif
    criticalSignalStr = ' '
    numCriticalSignals = 0
    if ( swLevel > -1 ) &
      & call outputNamedValue( 'critical module', critical_module_str )
    if ( swLevel > -1 ) &
      & call outputNamedValue( 'critical bands', ChunkDivideConfig%criticalBands )
    if ( swLevel > -1 ) &
      & call outputNamedValue( 'size(signals)', size(signals) )
    call allocate_test( isSignalCritical, size(signals), 'isSignalCritical', &
      & ModuleName )
    do signalIndex=1, size(signals)
      isSignalCritical(signalIndex) = &
        & .not. isNotACriticalBand( signals(signalIndex)%Band )
      if ( swLevel > -1 .and. .false. ) then
        call outputNamedValue( 'signal index', signalIndex )
        call GetSignalName( signalIndex, signal_full )
        call outputNamedValue( 'signal', signal_full )
        call outputNamedValue( 'band', signals(signalIndex)%Band )
        call outputNamedValue( 'critical?', &
          & .not. isNotACriticalBand( signals(signalIndex)%Band ) )
        call get_string( modules(signals(signalIndex)%instrumentModule)%name, &
          &  module_str, strip=.true. )
        call outputNamedValue( 'module', module_str )
      endif
      if ( len_trim(ChunkDivideConfig%criticalBands) > 0 ) then
        ! See if the band is one of the critical bands
        if ( isNotACriticalBand( signals(signalIndex)%Band ) ) cycle
      elseif ( len_trim(critical_module_str) > 0 ) then
        call get_string( modules(signals(signalIndex)%instrumentModule)%name, &
          & signal_full, strip=.true. )
        module_str = lowercase(signal_full)
        if ( module_str /= critical_module_str ) cycle
      endif
      numCriticalSignals = numCriticalSignals + 1
      call GetSignalName( signalIndex, signal_full )
      criticalSignalStr = catLists( criticalSignalStr, signal_full )
    enddo
    if ( swLevel > -1 ) &
      call Dump_Signals( signals, isSignalCritical, 1 )
    call deallocate_test( isSignalCritical,'isSignalCritical', &
      & ModuleName )
    if ( numCriticalSignals < 1 ) return
    call allocate_test( criticalSignals, numCriticalSignals, 'criticalSignals', &
      & ModuleName )
    call List2Array (criticalSignalStr, criticalSignals,  countEmpty=.true. )
    if ( swLevel > -1 ) &
      & call dump( criticalSignals, 'critical Signals chosen', &
      & options=what_options(trim=.true.) )
  end subroutine ChooseCriticalSignals

  ! -------------------------------------------- NoteL1BRADChanges -----
  subroutine NoteL1BRADChanges ( obstructions, mafRange, filedatabase )
    use MLSSignals_M, only: Dumpsignals=>dump, Getsignalindex, Getsignalname, &
      & GetNameofsignal, Signals
    use MLSStringlists, only: NumstringElements, GetstringElement
    use MLSStrings, only: Lowercase
    use Parse_Signal_M, only: Parse_Signal
    ! This routine notes any lack of good data for one of the
    ! signals, and, depending on sensitivity,
    ! augments the database of obstructions
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    type (Obstruction_T) :: NEWOBSTRUCTION ! In progrss
    type (MAFRange_T) :: MAFRange
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    ! Local variables
    ! This next sets the threshold of what it takes to be an obstruction
    ! At the most sensitive, hair-trigger level, any change in goodness
    ! of any signal would suffice
    ! Otherwise, a dead-zone in any signal must persist for at least a
    ! duration of config%maxGap before we declare it an obstruction
    ! and add it to the database
    logical, parameter :: ANYCHANGEISOBSTRUCTION = .false.
    logical :: choseCriticalSignals
    integer :: critical_index
    character(len=40) :: critical_module_str
    integer, dimension(:), pointer            :: criticalIndices => null() ! taken from config or criticalModules
    character(len=160), dimension(:), pointer :: criticalSignals => null() ! taken from config or criticalModules
    integer :: critical_sub_index
    integer, pointer, dimension(:) :: goods_after_gap
    integer, pointer, dimension(:) :: goodness_changes
    logical, dimension(size(signals)) :: good_after_maxgap
    logical, dimension(size(signals)) :: good_signals_now
    logical, dimension(size(signals)) :: good_signals_last
    integer, dimension(size(signals)) :: howlong_nogood
    integer :: i
    integer :: maf              ! mafrange(1) <= maf <= mafrange(2)
    integer :: maf_index   ! 1 <= maf_index <= mafrange(2)-mafrange(1)+1
    integer :: mafset
    integer :: mafset_end
    integer :: mafset_start
    integer, parameter :: MAXMAFSINSET = 250
    character(len=40) :: module_str
    integer :: nmafs
    integer :: nmafsets
    integer :: num_goods_after_gap
    integer :: num_goodness_changes
    logical, dimension(:), pointer  :: or_valids_buffer
    integer :: SignalIndex
    character(len=40) :: signal_full
    character(len=40) :: signal_str
    integer, pointer :: Signal_Indices(:)         ! Indices in the signals
    logical, dimension(:,:), pointer  :: signals_buffer
    logical :: this_maf_valid
    logical, dimension(:), pointer  :: valids_buffer
    ! Won't check if there are no files at all
    if ( .not. associated(filedatabase) ) return
    ! Won't check if there are no radiance files
    if ( .not. any(filedatabase%content == 'l1brad') ) return
    swlevel = SwitchDetail(switches, 'chu')
    nmafs = mafRange%Expanded(2) - mafRange%Expanded(1) + 1
    choseCriticalSignals = .false.
    ! Here we will loop over the signals database
    ! Searching for
    ! (1) mafs where there is no good data for any of the signals
    ! (2) mafs where the good data is switched from one signal to another
    ! (See Construct QuantityTemplates below line 253)

    ! (Possibly) time-consuming step:
    ! Read through the l1b radiance files for all signals,
    ! noting which ones are good, which not

    nullify ( signals_buffer, goods_after_gap, goodness_changes )
    call allocate_test( &
      & signals_buffer, nmafs , size(signals), &
      & 'signals_buffer', ModuleName )
    call allocate_test( goods_after_gap, nmafs,&
      & 'goods_after_gap', ModuleName )
    call allocate_test( goodness_changes, nmafs,&
      & 'goodness_changes', ModuleName )
    good_signals_now = .false.   ! Initializing
    signals_buffer = .false.
    if ( swLevel > -1 ) then
      call outputNamedValue ( 'shape(signals_buffer)', shape(signals_buffer) )
      call outputNamedValue ( 'shape(goods_after_gap)', shape(goods_after_gap) )
      call outputNamedValue ( 'shape(goodness_changes)', shape(goodness_changes) )
    endif
    if ( associated(ChunkDivideConfig%criticalSignals) ) then
      criticalSignals => ChunkDivideConfig%criticalSignals
      if ( swLevel > -1 ) then
        call output('Critical Signals from config', advance='yes' )
        call Dump (  criticalSignals )
      endif
    elseif ( ChunkDivideConfig%chooseCriticalSignals ) then
      ! We can choose as critical signals all the signals belonging
      ! to a critical module
      call chooseCriticalSignals( criticalSignals )
      choseCriticalSignals = .true.
      if ( swLevel > -1 ) then
        call Dump (  criticalSignals, 'Critical Signals you chose', &
          & options=what_options(trim=.true.) )
      endif
    endif
    ! We'll need the index number for each of the critical signals
    call allocate_test( criticalIndices, size(criticalSignals),&
      & 'critical indices', ModuleName )
    do critical_index = 1, size(criticalSignals)
      call GetSignalIndex( criticalSignals(critical_index), &
        & criticalIndices(critical_index) )
      if ( swLevel < VERBOSETHRESHOLD ) cycle
      call outputnamedValue( 'critical signal string', criticalSignals(critical_index) )
      call outputnamedValue( 'its index', criticalIndices(critical_index) )
    enddo
    if ( swLevel >= VERBOSETHRESHOLD ) then
      ! call dumpSignals( signals )
      call GetNameOfSignal ( signals(1), signal_str )
      call outputnamedValue( 'signal(1)%name', signal_str )
      call dump( criticalSignals, 'criticalSignals' )
      call dump( criticalIndices, 'criticalIndices' )
    endif
    do signalIndex=1, size(signals)
      call get_string( lit_indices(ChunkDivideConfig%criticalModules), &
        & signal_full, &
        & strip=.true. )
      critical_module_str = lowercase(signal_full)
      call get_string( modules(signals(signalIndex)%instrumentModule)%name, &
        & signal_full, &
        & strip=.true. )
      module_str = lowercase(signal_full)
      if ( swLevel >= VERBOSETHRESHOLD ) &
        & call outputNamedValue( 'module', module_str )
      if ( ChunkDivideConfig%criticalModules == l_none ) then
        ! No module is critical for signal data being good
      elseif ( module_str /= critical_module_str ) then
        cycle
      elseif ( .not. any(signalIndex == criticalIndices) ) then
        cycle
      endif
      if ( swLevel >= VERBOSETHRESHOLD ) &
        & call dumpSignals( signals(signalIndex) )
      if ( swLevel >= VERBOSETHRESHOLD ) &
        & call outputNamedValue( 'critical module', critical_module_str )
      if ( nmafs <= MAXMAFSINSET ) then
        good_signals_now(signalIndex) = &
          & any_good_signaldata ( signalIndex, signals(signalIndex)%sideband, &
          &   filedatabase, mafRange%Expanded(1), mafRange%Expanded(2), &
          &   signals_buffer(:,signalIndex), mafRange%Expanded )
      else
        nmafsets = (nmafs-1)/MAXMAFSINSET + 1   ! A trick--mafset_start
        mafset_end = mafRange%Expanded(1) - 1   ! is mafRange(1)
        do mafset=1, nmafsets
          mafset_start = mafset_end + 1
          mafset_end = min ( mafRange%Expanded(2), mafset_end + MAXMAFSINSET )
          good_signals_now(signalIndex) = &
            &   any_good_signaldata ( signalIndex, &
            &   signals(signalIndex)%sideband, &
            &   filedatabase, mafset_start, mafset_end, &
            &   signals_buffer(:,signalIndex), mafRange%Expanded )
        enddo
      endif
    enddo
    if ( swlevel > -1 ) &
      & call outputNamedValue ( 'shape(signals_buffer)', shape(signals_buffer) )
    if ( swlevel >= VERBOSETHRESHOLD ) &
      & call dump ( signals_buffer, 'signals_buffer' )

    ! Task (1a): Find mafs where there is at least one signal which
    ! changes from either nogood to good or from good to nogood
    ! compared with the last maf
    ! Task (1b): Find mafs where there is at least one signal which
    ! changes from nogood to good after a dead zone
    ! lasting at least maxGap mafs
    num_goodness_changes = 0
    num_goods_after_gap = 0
    howlong_nogood       = 0
    good_after_maxgap = .false.
    do maf = mafRange%Expanded(1), mafRange%Expanded(2)
      maf_index = maf - mafRange%Expanded(1) + 1
      do signalIndex=1, size(signals)
        if ( .not. any(signalIndex == criticalIndices) ) cycle
        good_signals_now(signalIndex) = signals_buffer(maf_index, signalIndex)
        if ( .not. good_signals_now(signalIndex) ) &
          & howlong_nogood(signalIndex) = howlong_nogood(signalIndex) + 1
        good_after_maxgap(signalIndex) = &
          & good_signals_now(signalIndex) &
          & .and. &
          & ( howlong_nogood(signalIndex) > ChunkDivideConfig%maxGap )
        if ( good_signals_now(signalIndex) ) howlong_nogood(signalIndex) = 0
      enddo
      if ( maf /= mafRange%Expanded(1) ) then
        if ( Any(good_signals_now .neqv. good_signals_last) ) then
          num_goodness_changes = num_goodness_changes + 1
          goodness_changes(num_goodness_changes) = maf
        endif
      endif
      if ( Any(good_after_maxgap) ) then
        num_goods_after_gap = num_goods_after_gap + 1
        goods_after_gap(num_goods_after_gap) = maf
      endif
      good_signals_last = good_signals_now
    enddo
    call outputNamedValue( 'num_goods_after_gap', num_goods_after_gap )

    ! Task (2): Find regions where there is no signal among at least one
    ! of the critical signals
    if ( associated(criticalSignals) ) then
      ! What we're about to do is:
      ! Assume we're given an array of char strings
      ! array = [str_1, str_2, ..]
      ! and each string contains comma-separated signals
      ! str_1 = 'sig_1,sig_2,..'
      ! Then the rule is that for a maf to pass the critical signals test,
      ! for each array element at least one of the signals within that
      ! string must be turned on for that maf
      ! Thus the shorthand formula that we're and-ing the array elements
      ! and or-ing the list elements
      ! You might prefer the formula that a maf is invalid if for any of the
      ! array elements, none of its signals is turned on
      !
      ! Unfortunately, if we chose our critical signals based on a critical module
      ! the role the string containing comma-separated signals is played by
      ! by the array so in that case we'll need to be or-ing the array elements,
      ! not and-ing them
      nullify ( valids_buffer )
      call allocate_test( &
        & valids_buffer, nmafs , &
        & 'valids_buffer', ModuleName)
      nullify ( or_valids_buffer )
      call allocate_test( &
        & or_valids_buffer, nmafs , &
        & 'or_valids_buffer', ModuleName)
      !
      valids_buffer = .not. choseCriticalSignals ! .true.
      if ( SwLevel >= VERBOSETHRESHOLD ) then
        call output ( 'Checking for critical signals: ')
        do critical_index=1, size(criticalSignals)
          call output ( trim(criticalSignals(critical_index)), &
            & advance='yes')
        enddo
        call output ( &
          & 'Which signals match the incomplete signal string ' &
          &  // 'weve been given', advance='yes')
      endif
      do critical_index=1, size(criticalSignals)
        or_valids_buffer = .false.
        do critical_sub_index=1, NumStringElements( &
          & trim(criticalSignals(critical_index)), .FALSE.)
          call GetStringElement( &
            & trim(criticalSignals(critical_index)), signal_str, &
            & critical_sub_index, .FALSE. )
          ! Which signals match the incomplete signal string we've been given?
          nullify(Signal_Indices)
          call Parse_signal(signal_str, signal_indices)
          if ( .not. associated(Signal_Indices) ) exit
          if ( SwLevel >= VERBOSETHRESHOLD ) then
            call output ( 'Signal_Indices: ')
            call output ( Signal_Indices, advance='yes')
            do i=1, size(Signal_Indices)
              call GetSignalName(Signal_Indices(i), signal_full)
              call output ( 'Full Signal Name: ')
              call output ( trim(signal_full), advance='yes')
            enddo
          endif
          do maf = mafRange%Expanded(1), mafRange%Expanded(2)
            maf_index = maf - mafRange%Expanded(1) + 1
            this_maf_valid = .false.   ! not valid unless at least one signal
            do i=1, size(Signal_Indices)
              signalIndex = Signal_Indices(i)
              this_maf_valid = &
                & this_maf_valid .or. signals_buffer(maf_index, signalIndex)
            enddo
            or_valids_buffer(maf_index) = &
              & or_valids_buffer(maf_index) .or. this_maf_valid
          enddo
          call deallocate_test(Signal_Indices, 'Signal_Indices', ModuleName)
        enddo
        ! An ad-hoc filter to ignore bands that are simply switched off
        ! if ( swLevel > -1 ) call outputNamedValue ( 'count(or_valids_buffer)', count(or_valids_buffer), advance='yes')
        ! if ( count(or_valids_buffer) < &
        !  & ( mafRange%Expanded(2) - mafRange%Expanded(1) ) / 50 ) cycle
        ! call dump ( or_valids_buffer, 'or_valids_buffer' )
        if ( choseCriticalSignals ) then
          valids_buffer = &
            & valids_buffer .or. or_valids_buffer
        else
          valids_buffer = &
            & valids_buffer .and. or_valids_buffer
        endif
      enddo
      if ( swLevel >= VERBOSETHRESHOLD ) call dump ( valids_buffer, 'valids_buffer' )
      if ( swLevel > -1 ) then
        call outputNamedValue ( 'count(valids_buffer)', count(valids_buffer), advance='yes')
        call output ( 'Before converting valids to obstructions', advance='yes' )
        call dump(obstructions)
      endif
      call ConvertFlagsToObstructions ( valids_buffer, obstructions, &
        & mafRange%Expanded, obstructionType='critical radiances' )
      if ( swLevel > -1 ) then
        call output ( 'After converting valids to obstructions', advance='yes' )
        call dump(obstructions)
      endif
      call deallocate_test(valids_buffer, 'valids_buffer', ModuleName)
      call deallocate_test(or_valids_buffer, 'or_valids_buffer', ModuleName)
      if ( .not. associated(ChunkDivideConfig%criticalSignals) )  &
        & call deallocate_test( criticalSignals, 'criticalSignals', &
        & ModuleName )
    endif

    ! Depending on sensitivity, add these to Obstructions database
    if ( ANYCHANGEISOBSTRUCTION .and. num_goodness_changes > 0 ) then
      do mafset = 1, num_goodness_changes
        newObstruction%range = .false.
        newObstruction%mafs(1) = goodness_changes(mafset)
        newObstruction%mafs(2) = 0    ! For overzealous Lahey uninitialized checking
        call AddObstructionToDatabase ( obstructions, newObstruction )
        if ( swLevel > -1 ) &
          call outputNamedValue( &
          & 'maf where any change made obstruction', goodness_changes(mafset) )
      enddo
    elseif ( num_goods_after_gap > 0 ) then
      do mafset = 1, num_goods_after_gap
        newObstruction%range = .false.
        newObstruction%mafs(1) = goods_after_gap(mafset)
        newObstruction%mafs(2) = 0    ! For overzealous Lahey uninitialized checking
        call AddObstructionToDatabase ( obstructions, newObstruction )
        if ( swLevel > -1 ) &
          call outputNamedValue( &
          & 'maf where num goods after gap made obstruction', goods_after_gap(mafset) )
      enddo
    endif
    ! OK, we have the mafs where the goodness changes, now what?
    call deallocate_test( signals_buffer, 'signals_buffer', ModuleName )
    call deallocate_test( goods_after_gap, 'goods_after_gap', ModuleName )
    call deallocate_test( goodness_changes, 'goodness_changes', ModuleName )
    call deallocate_test( criticalIndices, 'critical indices', ModuleName )
  end subroutine NoteL1BRADChanges

  ! ----------------------------------- ConvertFlagsToObstructions -----
  subroutine ConvertFlagsToObstructions ( valid, obstructions, &
    & mafRange, obstructionType )
    ! This routine takes an array of logicals indicating good/bad data
    ! and converts it into obstruction information.
    logical, dimension(:), intent(in)           :: VALID
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    integer, dimension(:), intent(in), optional :: MAFRANGE
    character(len=*), intent(in), optional      :: obstructionType

    ! Local variables
    logical :: LASTONEVALID           ! Flag
    integer :: MAF                    ! Loop counter
    type (Obstruction_T) :: NEWOBSTRUCTION ! In progrss
    character(len=64)    :: obstructionTrigger
    integer :: OFFSET                 ! MAF index offset

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    lastOneValid = .true.
    offset = 0
    if ( present(mafRange) ) offset = mafRange(1)
    obstructionTrigger = &
      & 'maf where transition from bad to good made obstruction'
    if ( present(obstructionType) ) obstructionTrigger = &
      & 'maf where transition from bad to good ' // trim(obstructionType) // &
      & 'made obstruction'

    do maf = 1, size(valid)
      if ( valid(maf) .neqv. lastOneValid ) then
        ! A transition either from good to bad or bad to good
        if ( .not. valid(maf) ) then
          ! From good to bad
          newObstruction%range = .true.
          newObstruction%mafs(1) = maf - 1 + offset
        else
          newObstruction%mafs(2) = maf - 1 + offset
          call AddObstructionToDatabase ( obstructions, newObstruction )
          if ( swLevel > -1 ) &
            call outputNamedValue( &
            & trim(obstructionTrigger), maf )
        end if
      end if
      lastOneValid = valid(maf)
    end do
    ! Check that the two mafs after the obstruction end are valid
!     maf = newObstruction%mafs(2) + 1
!     if ( maf <= size(valid) ) then
!       if ( .not. valid(maf) ) print *, 'Warning not valid at maf ', maf
!       newObstruction%mafs(2) = newObstruction%mafs(2) + 1
!     endif
!     maf = newObstruction%mafs(2) + 1
!     if ( maf <= size(valid) ) then
!       if ( .not. valid(maf) ) print *, 'Warning not valid at maf ', maf
!       newObstruction%mafs(2) = newObstruction%mafs(2) + 1
!     endif

    ! Make sure any range at the end gets added
    if ( .not. lastOneValid ) then
      newObstruction%mafs(2) = size(valid) - 1 + offset
      call AddObstructionToDatabase ( obstructions, newObstruction )
      if ( swLevel > -1 ) &
        call outputNamedValue( &
        & trim(obstructionTrigger), size(valid)-1 )
    end if
  end subroutine ConvertFlagsToObstructions

  ! ------------------------------------------------ SurveyL1BData -----
  subroutine SurveyL1BData ( processingRange, filedatabase, mafRange )
    ! This goes through the L1B o/a data file and tries to spot possible
    ! obstructions.
    ! If any breaks in valid data are revealed, they will be added to 
    ! the obsstructions database
    use Monotone, only: isMonotonic
    type (TAI93_Range_T), intent(in)        :: PROCESSINGRANGE
    type (MLSFile_T), dimension(:), pointer :: FILEDATABASE
    type (MAFRange_T), intent(out)          :: MAFRange

    ! Local variables
    ! The following 5 l1b data sets will be examined for possible breaks
    ! in valid data
    type (L1BData_T) :: SCVEL           ! Read from L1BOA file
    type (L1BData_T) :: TAITIME         ! Read from L1BOA file
    type (L1BData_T) :: TPGEODALT       ! Read from L1BOA file
    type (L1BData_T) :: TPGEODANGLE     ! Read from L1BOA file
    type (L1BData_T) :: TPORBY          ! Read from L1BOA file

    character(len=10) :: MODNAMESTR     ! Module name

    integer :: FLAG                     ! From L1B
    integer :: MAF                      ! Loop inductor
    integer :: MOD                      ! Loop inductor
    integer :: NOMAFS                   ! Number of MAFs in processing range
    integer :: NOMAFSREAD               ! From L1B

    logical :: THISONEVALID             ! To go into valid
    logical, dimension(:), pointer :: ANGLEWASSMOOTHED ! Flag for each MAF
    logical, dimension(:), pointer :: VALID ! Flag for each MAF
    logical, dimension(:), pointer :: WASSMOOTHED ! Flag for each MAF

    real(r8) :: ORBYMAX                 ! Maximum value of orbY each maf
    real(r8) :: SCANMAX                 ! Range of scan each maf
    real(r8) :: SCANMIN                 ! Range of scan each maf
    real(r8) :: SCVELMAX                ! Maximum value of abs(sc/VelECR)

    integer   ::                       l1b_hdf_version
    character(len=namelen) ::          MAF_start
    character(len=namelen) ::          sc_vel
    character(len=namelen) ::          tp_alt
    character(len=namelen) ::          tp_orbY
    character(len=namelen) ::          tp_angle
    type (MLSFile_T), pointer       :: L1BFile
    ! Executable code
    swlevel = switchDetail(switches, 'chu' )

    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) &
      & call MLSMessage  ( MLSMSG_Error, ModuleName, &
        & "Can't make progress in SurveyL1BData without L1BOA files" )
    if ( swLevel > -1 ) then
      call output( 'here is the L1BOA', advance='yes' )
      call dump( L1BFile )
    endif
    l1b_hdf_version = L1BFile%HDFVersion
    MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
      .false. )
    ! tp_angle = &
    !  & AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
    !  & l1b_hdf_version, &
    !  & .false. )
    ! Read time from the L1BOA file
    call ReadL1BData ( L1BFile, trim(MAF_start), taiTime, &
      & noMAFsRead, flag, &
      & dontPad=DONTPAD )
    call smoothOutDroppedMAFs(taiTime%dpField)

    ! We shall assume that all l1b files cover this range
    mafRange%L1BCover(1) = taiTime%FirstMAF
    mafRange%L1BCover(2) = taiTime%NoMAFs - 1 ! Because it starts with 0

    ! Deduce the first and last MAFs covered the L2 Processing Range
    call Hunt ( taiTime%dpField(1,1,:), &
      & (/ processingRange%startTime, processingRange%endTime /), &
      & mafRange%L2Cover, allowTopValue=.true., allowBelowValue=.true. )

    ! Check the validity of the MAF range returned
    if ( mafRange%L2Cover(1) == taiTime%noMAFs ) then
      call outputNamedValue ( 'processingRange', &
        & (/ processingRange%startTime, processingRange%endTime /) )
      call outputNamedValue ( 'noMAFsRead', noMAFsRead )
      call dump ( taiTime%dpField(1,1,:), 'MAF_start' )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L1B data ends before requested processing range' )
    endif
    if ( mafRange%L2Cover(1) == taiTime%noMAFs ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L1B data ends before requested processing range' )
    mafRange%L2Cover = min ( noMAFsRead-1, max ( 0, mafRange%L2Cover - 1 ) )          ! Index from zero
    if ( mafRange%L2Cover(2) < mafRange%L2Cover(1) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L2 mafRange(2) < L2 mafRange(1)' )
    elseif ( mafRange%L2Cover(2) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L2 mafRange(2) < 1' )
    endif

    mafRange%Expanded(1) = mafRange%L2Cover(1)
    mafRange%Expanded(2) = mafRange%L2Cover(2)
    ! Expand the range covered by L2 with possible overlaps
    if ( ChunkDivideConfig%allowPriorOverlaps ) then
      mafRange%Expanded(1) = max( mafRange%L1BCover(1), &
        & mafRange%L2Cover(1) - nint(ChunkDivideConfig%lowerOverlap) )
    endif
    if ( ChunkDivideConfig%allowPostOverlaps ) then
      mafRange%Expanded(2) = min( mafRange%L1BCover(2), &
        & mafRange%L2Cover(2) + nint(ChunkDivideConfig%upperOverlap) )
    endif
    ! Why was this here? if ( sharedPCF ) return
    noMAFS = mafRange%Expanded(2) - mafRange%Expanded(1) + 1
    ! Now look through the L1B data, first look for scan problems
    if ( ChunkDivideConfig%criticalModules /= l_none ) then
      call outputNamedValue( 'Looking for scan problems; swLevel ', swlevel )
      nullify ( valid )
      nullify ( anglewasSmoothed )
      nullify ( wasSmoothed )
      call Allocate_test ( valid, noMAFs, 'valid', ModuleName )
      call Allocate_test ( anglewasSmoothed, noMAFs, 'angleWasSmoothed', ModuleName )
      call Allocate_test ( wasSmoothed, noMAFs, 'wasSmoothed', ModuleName )
      ! How to choose initial values for valid:
      ! It depends on whether we insist "both" modules are critical
      ! or not.
      ! In practice we will never insist.
      if ( ChunkDivideConfig%criticalModules == l_both ) then
        valid = .true.
      else
        valid = .false.
      endif
      do mod = 1, size(modules)
        call get_string ( modules(mod)%name, modNameStr, strip=.true. )
        sc_vel = 'sc/VelECR'
        tp_alt = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAlt', &
          & l1b_hdf_version, .false. )
        tp_orby = AssembleL1BQtyName ( trim(modNameStr)//'.tpOrbY', &
          & l1b_hdf_version, .false. )
        tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
          & l1b_hdf_version, .false. )
        if ( .not. modules(mod)%spacecraft .and. &
          & ( any ( ChunkDivideConfig%criticalModules == (/ l_either, l_both /) ) .or. &
          &   lit_indices(ChunkDivideConfig%criticalModules) == modules(mod)%name ) ) then
          ! Read the s/c velocity
          call ReadL1BData ( L1BFile, trim(sc_vel), &
            & scVel, noMAFsRead, flag, &
            & firstMAF=mafRange%Expanded(1), lastMAF=mafRange%Expanded(2), &
            & dontPad=DONTPAD )
          ! Read the tangent point altitude
          call ReadL1BData ( L1BFile, trim(tp_alt), &
            & tpGeodAlt, noMAFsRead, flag, &
            & firstMAF=mafRange%Expanded(1), lastMAF=mafRange%Expanded(2), &
            & dontPad=DONTPAD )
          ! Read the out of plane distance
          call ReadL1BData ( L1BFile, trim(tp_orby), &
            & tpOrbY, noMAFsRead, flag, &
            & firstMAF=mafRange%Expanded(1), lastMAF=mafRange%Expanded(2), &
            & dontPad=DONTPAD )
          call ReadL1BData ( L1BFile, trim(tp_angle), &
            & tpGeodAngle, noMAFsRead, flag, &
            & firstMAF=mafRange%Expanded(1), lastMAF=mafRange%Expanded(2), &
            & dontPad=DONTPAD )
          if ( .not. ChunkDivideConfig%skipL1BCheck ) then
            if ( .not. isMonotonic(tpGeodAngle%dpField(1,1,1:noMAFsRead) ) &
              & ) then
              call MLSMessage ( PhiNotMonotonicFun(), ModuleName, &
                & 'tpGeodAngle is not monotonic--may quit', MLSFile=L1BFile )
            endif
          endif
          ! For a reason forgotten in the fog of antiquity we, and mlsl2,
          ! smooth out discontinuities in the tp quantities if possible.
          call smoothOutDroppedMAFs(tpGeodAngle%dpField, angleWasSmoothed, &
            & monotonize=.true.)
          call smoothOutDroppedMAFs(tpGeodAlt%dpField, wasSmoothed)
          call smoothOutDroppedMAFs(tpOrbY%dpField)
          if ( swLevel > 1 ) then
            call dump( angleWasSmoothed, 'smoothed because of geodetic angle' )
            call dump( wasSmoothed, 'smoothed because of geodetic altitude' )
          endif
          wasSmoothed = (wasSmoothed .or. angleWasSmoothed)
          ! Consider the scan range in each MAF in turn
          do maf = 1, noMAFs
            scanMax = maxval ( tpGeodAlt%dpField(1,:,maf) )
            scanMin = minval ( tpGeodAlt%dpField(1,:,maf) )
            orbyMax = maxval ( tpOrbY%dpField(1,:,maf) )
            scVelMax = maxval ( abs(scVel%dpField(1,:,maf)) )

            thisOneValid = ( &
              &   scanMin >= ChunkDivideConfig%scanLowerLimit(1) .and. &
              &   scanMin <= ChunkDivideConfig%scanLowerLimit(2) ) .and. ( &
              &   scanMax >= ChunkDivideConfig%scanUpperLimit(1) .and. &
              &   scanMax <= ChunkDivideConfig%scanUpperLimit(2) )
            if ( ChunkDivideConfig%maxOrbY > 0.0 ) then
              thisOneValid = thisOneValid .and. orbYMax < ChunkDivideConfig%maxOrbY
            end if
            if ( ChunkDivideConfig%maxscVel > 0.0 ) then
              thisOneValid = thisOneValid .and. scVelMax < ChunkDivideConfig%maxscVel
            end if
            ! this one is not valid if it's valid only by virtue
            ! of having been smoothed
            thisOneValid = thisOneValid .and. .not. wasSmoothed(maf)
            if ( ChunkDivideConfig%criticalModules == l_both ) then
              valid(maf) = valid(maf) .and. thisOneValid
            else
              valid(maf) = valid(maf) .or. thisOneValid
            end if
          end do                        ! Maf loop
          call DeallocateL1BData ( taiTime )
          call DeallocateL1BData ( tpgeodalt )
          call DeallocateL1BData ( tpgeodangle )
          call DeallocateL1BData ( tporby )
          call DeallocateL1BData ( scVel )
        end if                          ! Consider this module
      end do                            ! Module Loop
      ! Convert this information into obstructions and tidy up.
      if ( swLevel > 1 ) call Dump( valid, 'valid array passed to ConvertFlagsToObstructions' )
      maf = FindFirst ( .not. valid )
      call OutputNamedValue ( 'First not valid', maf )
      maf = FindLast ( .not. valid )
      call OutputNamedValue ( 'Last not valid', maf )
      call ConvertFlagsToObstructions ( valid, obstructions, &
        & obstructionType='scan' )
      call Deallocate_test ( valid, 'valid', ModuleName )
      call Deallocate_test ( wasSmoothed, 'wasSmoothed', ModuleName )
      call Deallocate_test ( angleWasSmoothed, 'angleWasSmoothed', ModuleName )
    end if                              ! Consider scan issues

    ! Here we look at radiances and switch changes.
    if ( swLevel > -1 ) call output ( 'NoteL1BRADChanges', advance='yes' )
    if ( .not. ChunkDivideConfig%skipL1BCheck ) &
      call NoteL1BRADChanges ( obstructions, mafRange, filedatabase )

    ! Sort the obstructions into order; prune them of repeats, overlaps etc.
    call PruneObstructions ( obstructions )

    ! Tidy up
    call DeallocateL1BData ( taiTime )

  contains
    function PhiNotMonotonicFun() result( SEVERITY )
    ! Default severity (probably ERROR) incurred
    ! if config options set the crashIfPhiNotMono flag
    ! Args
    integer :: SEVERITY
    ! if ( switchDetail(switches,'nmono') > -1 ) then ! be lenient
    if ( ChunkDivideConfig%crashIfPhiNotMono ) then ! be lenient
      severity = severityifphinotmono
    else
      severity = min( MLSMSG_Warning, severityifphinotmono )
    endif
    end function PhiNotMonotonicFun
  end subroutine SurveyL1BData

  ! -------------------------------------------- PruneObstructions -----
  subroutine PruneObstructions ( obstructions )
    ! This routine merges overlapping range obstructions and deletes
    ! wall obstructions inside ranges.  The job is made easier
    ! by sorting the obstructions into order
    use Allocate_Deallocate, only: Test_Allocate
    type(Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

    ! Local variables
    integer :: I,J                  ! Loop counters
    type (Obstruction_T) :: newObs  ! New Obstruction
    logical :: FOUNDONE             ! Found at least one
    integer :: STATUS               ! Flag from allocate

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    ! If no obstructions make sure allocate to size zero, not just unassociated pointer
    if ( .not. associated(obstructions) ) then
      allocate ( obstructions(0), stat=status )
      call test_allocate ( status, ModuleName, 'obstructions', &
        & uBounds = 0, elementSize = storage_size(obstructions) / 8 )
      return
    end if

    ! Otherwise, do the tidying up
    outerLoop: do
      foundOne = .false.
      i = 0
      call SortObstructions(obstructions)
      middleLoop: do
        i = i + 1
        if ( i >= size(obstructions) ) exit middleLoop
        j = i
        innerLoop: do
          j = j + 1
          if ( j > size(obstructions) ) exit innerLoop
          if ( all ( obstructions((/i,j/))%range ) ) then
            ! --------------------------- ( Range, range )
            if ( obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
              ! Combine overlapping range obstructions
              newObs%range = .true.
              newObs%mafs(1) = obstructions(i)%mafs(1)
              newObs%mafs(2) = &
                & max ( obstructions(i)%mafs(2), obstructions(j)%mafs(2) )
              ! Must delete these in order: otherwise
              ! if deleted i first where i < j, index would
              ! no longer be "j" afterwards
              call DeleteObstruction ( obstructions, j )
              call DeleteObstruction ( obstructions, i )
              call AddObstructionToDatabase ( obstructions, newObs )
              call SortObstructions ( obstructions )
              foundOne = .true.
              exit middleLoop
            end if
          else if ( obstructions(i)%range .and. .not. obstructions(j)%range ) then
            ! --------------------------- ( Range, wall )
            if ( obstructions(j)%mafs(1) >= obstructions(i)%mafs(1) .and. &
              &  obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
              ! Delete wall obstruction inside range
              call DeleteObstruction ( obstructions, j )
              foundOne = .true.
              exit middleLoop
            end if
          else
            ! --------------------------- ( Wall, range ) or ( Wall, wall )
            ! Becuase the obstructions are in order, we know in the wall, range
            ! case that the wall must be at the start of the range, not inside it.
            if ( obstructions(i)%mafs(1) == obstructions(j)%mafs(1) ) then
              ! Delete wall obstruction at start of a range or at another wall
              call DeleteObstruction ( obstructions, i )
              foundOne = .true.
              exit middleLoop
            end if
          end if
          ! I'm pretty sure this covers all the possibilities.  It might seem
          ! not at first glance, but I think the fact that I always re-sort the
          ! obstructions into order means that the above code does catch everything.

        end do innerLoop
      end do middleLoop
      if ( .not. foundOne ) exit outerLoop
    end do outerLoop

  end subroutine PruneObstructions

  !------------------------------------------- ChunkDivide_Orbital -----
  subroutine ChunkDivide_Orbital ( mafRange, filedatabase, chunks, root )
    ! integer, dimension(2), intent(in) :: MAFRANGE
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Trace_m, only: Trace_Begin, Trace_End
    type (MAFRange_T) :: MAFRange
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS
    integer, intent(in), optional :: Root ! of the l2cf tree, for tracing

    ! Local parameters
    real(r8), parameter :: HOMEACCURACY = 3.0 ! Try to hit homeGeodAngle within this
    ! (see homeHuntLoop warning below)
    ! Local variables
    type (L1BData_T) :: TAITIME         ! From L1BOA
    type (L1BData_T) :: TPGEODANGLE     ! From L1BOA

    character(len=10) :: MODNAMESTR     ! Home module name as string

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: CHUNK                    ! Loop counter
    integer :: FLAG                     ! From ReadL1B
    integer :: HOME                     ! Index of home MAF in array
    integer :: M1, M2                   ! MafRange%L2Cover + 1
    integer :: Me = -1                  ! String index for trace
    integer :: MEXP1, MEXP2             ! MafRange%Expanded + 1
    integer :: NOCHUNKSBELOWHOME        ! Used for placing chunks
    integer :: NOMAFSATORABOVEHOME      ! Fairly self descriptive
    integer :: NOMAFSBELOWHOME          ! Fairly self descriptive
    ! integer :: NOMAFS                   ! Number of MAFs to consider
    integer :: NOMAFSREAD               ! From ReadL1B
    integer :: ORBIT                    ! Used to locate home
    integer :: STATUS                   ! From allocate etc.
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: MAXLENGTH                ! Max length as integer (MAFs)

    integer, dimension(:), pointer :: NEWFIRSTMAFS ! For thinking about overlaps
    integer, dimension(:), pointer :: NEWLASTMAFS ! For thinking about overlaps

    real(r8) :: ANGLEINCREMENT          ! Increment in hunt for home
    real(r8) :: MAXANGLE                ! Of range in data
    real(r8) :: MAXTIME                 ! Time range in data
    real(r8) :: MAXV                    ! Either minTime or minAngle
    real(r8) :: MINANGLE                ! Of range in data
    real(r8) :: MINTIME                 ! Time range in data
    real(r8) :: MINV                    ! Either minTime or minAngle
    real(r8) :: TESTANGLE               ! Angle to check for
    real(r8) :: HOMEV                   ! Value of angle/time at home

    real(r8), dimension(:), pointer :: BOUNDARIES ! Used in placing chunks
    real(r8), dimension(:), pointer :: FIELD ! Used in placing chunks

    integer   ::                       l1b_hdf_version
    character(len=namelen) ::         MAF_start, tp_angle
    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code
    if ( present(root) ) call trace_begin ( me, 'ChunkDivide_orbital', root, &
        & cond=toggle(gen) .and. levels(gen) > 2 )

    swlevel = switchDetail(switches, 'chu' )
    if ( swLevel > -1 ) then
      call output('Entering Orbital Chunk Divide', advance='yes')
      call dump( obstructions )
    endif
    ! Read in the data we're going to need
    call get_string ( lit_indices(ChunkDivideConfig%homeModule), modNameStr, strip=.true. )
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) &
      & call MLSMessage  ( MLSMSG_Error, ModuleName, &
        & "Can't make progress in ChunkDivide_Orbital without L1BOA files" )
    ! call dump(L1BFile)
    l1b_hdf_version = L1BFile%HDFVersion
    MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
      .false. )
    tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
      & l1b_hdf_version, &
      & .false. )
    if ( swLevel > -1 ) &
      & call output('Reading Geod Angle', advance='yes')
    call ReadL1BData ( L1BFile, trim(tp_angle), &
      & tpGeodAngle, noMAFsRead, flag, &
      & dontPad=DONTPAD )
    if ( swLevel > -1 ) &
      & call output('1st smoothing', advance='yes')
    call smoothOutDroppedMAFs(tpGeodAngle%dpField)
    if ( swLevel > -1 ) &
      & call output('2nd smoothing', advance='yes')
    call smoothOutDroppedMAFs(tpGeodAngle%dpField, monotonize=.true.)
    if ( swLevel > -1 ) &
      & call output('Reading tai Time', advance='yes')
    call ReadL1BData ( L1BFile, trim(MAF_start), &
      & taiTime, noMAFsRead, flag, &
      & dontPad=DONTPAD )
    call smoothOutDroppedMAFs(taiTime%dpField)
    ! noMAFs = mafRange%L2Cover(2) - mafRange%L2Cover(1) + 1
    m1 = mafRange%L2Cover(1) + 1
    m2 = mafRange%L2Cover(2) + 1
    mexp1 = mafRange%Expanded(1) + 1
    mexp2 = mafRange%Expanded(2) + 1

    minAngle = minval ( tpGeodAngle%dpField(1,1,m1:m2) )
    maxAngle = maxval ( tpGeodAngle%dpField(1,1,m1:m2) )
    minTime = minval ( taiTime%dpField(1,1,m1:m2) )
    maxTime = maxval ( taiTime%dpField(1,1,m1:m2) )
    if ( swLevel > -1 ) then
      call output ( 'Num MAFs in file: ' )
      call output ( noMAFsRead, advance='yes' )
      call output ( 'MAF time range: ' )
      call output ( minTime )
      call output ( ' : ' )
      call output ( maxTime, advance='yes' )
      call output ( 'Angle range: ' )
      call output ( minAngle )
      call output ( ' : ' )
      call output ( maxAngle, advance='yes' )
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! First try to locate the last MAF before the homeGeodAngle
    orbit = int ( tpGeodAngle%dpField(1,1,m1)/360.0 )
    if ( tpGeodAngle%dpField(1,1,m1) < 0.0 ) orbit = orbit - 1
    testAngle = ChunkDivideConfig%homeGeodAngle + orbit*360.0
    if ( ChunkDivideConfig%maxLengthFamily == PHYQ_Angle ) then
      angleIncrement = ChunkDivideConfig%maxLength
    else
      angleIncrement = 360.0
    end if

    if ( swLevel > -1 ) then
      call output ( ' orbit  ', advance='no' )
      call output ( orbit , advance='no' )
      call output ( '    testAngle  ', advance='no' )
      call output ( testAngle , advance='no' )
      call output ( '    angleIncrement  ', advance='no' )
      call output ( angleIncrement , advance='yes' )
    end if
    ! In my opinion (paw) here's what the following loop should do:
    ! Find the 1st MAF within HOMEACCURACY of home_angle
    ! where home_angle has been corrected for the starting orbit number
    ! Afterwards, the preceding MAFs must be divided among one or more
    ! chunks, and the same done with subsequent MAFs
    !
    ! Instead what it actually does is
    ! Find the 1st MAF within HOMEACCURACY of (home_angle + n*angleIncrement)
    ! where home_angle has been corrected for the starting orbit number
    ! In effect the home_angle is set only within an unknown number
    ! of angleIncrements
    ! While there may be few cases in which they don't do about as well
    ! let this be a warning
    homeHuntLoop: do
      if ( testAngle < minAngle ) then
        testAngle = testAngle + angleIncrement
        cycle
      endif
      if ( testAngle > maxAngle ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to establish a home major frame, using the first in your range' )
        home = m1
        exit homeHuntLoop
      end if
      ! Find MAF which starts before this test angle
      call Hunt ( tpGeodAngle%dpField(1,1,:), testAngle, home, nearest=.true.,&
        & allowTopValue = .true. )
      ! Now if this is close enough, accept it
      if ( abs ( tpGeodAngle%dpField(1,1,home) - &
        & testAngle ) < HomeAccuracy ) exit homeHuntLoop
      ! Otherwise, keep looking
      testAngle = testAngle + angleIncrement
    end do homeHuntLoop
    if ( swLevel > -1 ) then
      call output ( 'Test Angle  ' )
      call output ( testAngle , advance='yes' )
      call output ( 'Angle(home)  ' )
      call output ( tpGeodAngle%dpField(1,1,home) , advance='yes' )
      call output ( 'Home  ' )
      call output ( home  )
      call output ( 'Difference  ' )
      call output (  abs ( tpGeodAngle%dpField(1,1,home) - &
        & testAngle ) , advance='yes' )
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! OK, now we have a home MAF, get a first cut for the chunks
    ! We work out the chunk ends for each chunk according to how the
    ! maxLength field is specified.
    if ( ChunkDivideConfig%maxLengthFamily == PHYQ_MAFs ) then
      maxLength = nint ( ChunkDivideConfig%maxLength )
      noMAFsBelowHome = home - m1
      noChunksBelowHome = noMAFsBelowHome / maxLength
      if ( mod ( noMAFsBelowHome, maxLength ) /= 0 ) noChunksBelowHome = noChunksBelowHome + 1
      noMAFsAtOrAboveHome = m2 - home + 1
      if ( ChunkDivideConfig%noChunks == 0 ) then
         ! If user did not request specific number of chunks choose them
         noChunks = noChunksBelowHome + noMAFsAtOrAboveHome / maxLength
         if ( mod ( noMAFsAtOrAboveHome, maxLength ) /= 0 ) &
            & noChunks = noChunks + 1
      else
        ! User requested specific number of chunks
        noChunks = ChunkDivideConfig%noChunks
      end if

      ! Allocate the chunks
      allocate ( chunks(noChunks), stat=status )
      addr = 0
      if ( status == 0 .and. noChunks > 0 ) addr = transfer(c_loc(chunks(1)), addr)
      call test_allocate ( status, ModuleName, 'chunks', &
        & uBounds = noChunks, elementSize = storage_size(chunks) / 8, address=addr )

      ! Work out their positions
      do chunk = 1, noChunks
        chunks(chunk)%lastMAFIndex = home + &
          & ( chunk - noChunksBelowHome ) * maxLength - 1
        ! Subtract one to convert from index in array to index in file
      end do
    else
      ! For angle and time, they are similar enough we'll just do some stuff
      ! with pointers to allow us to use common code to sort them out
      select case ( ChunkDivideConfig%maxLengthFamily )
      case ( PHYQ_Angle )
        field => tpGeodAngle%dpField(1,1,:)
        minV = minAngle
        maxV = maxAngle
      case ( PHYQ_Time )
        field => taiTime%dpField(1,1,:)
        minV = minTime
        maxV = maxTime
      case ( PHYQ_MAFs)
      end select
      homeV = field(home)

      noMAFsBelowHome = -999
      noMAFsAtOrAboveHome = -999
      noChunksBelowHome = int ( &
        & ( homeV - minV ) / ChunkDivideConfig%maxLength )
      if ( homeV > minV ) noChunksBelowHome = noChunksBelowHome + 1
      if ( ChunkDivideConfig%noChunks == 0 ) then
        ! Choose the number of chunks ourselves
        noChunks = noChunksBelowHome + int ( &
          & ( maxV - homeV ) / ChunkDivideConfig%maxLength )
        if ( homeV + ChunkDivideConfig%maxLength * &
          & ( noChunks - noChunksBelowHome ) < maxV ) &
          & noChunks = noChunks + 1
        if ( (mexp1 < m1) .and. ( &
          & homeV + (ChunkDivideConfig%maxLength-1) + &
          & ChunkDivideConfig%maxLength * &
          & ( noChunks - 1 - noChunksBelowHome ) < maxV ) ) &
          & noChunks = noChunks + 1
      else
        noChunks = ChunkDivideConfig%noChunks
      end if

      ! Allocate the chunks
      allocate ( chunks(noChunks), stat=status )
      addr = 0
      if ( status == 0 .and. noChunks > 0 ) addr = transfer(c_loc(chunks(1)), addr )
      call test_allocate ( status, ModuleName, 'chunks', &
        & uBounds = noChunks, elementSize = storage_size(chunks) / 8, address=addr )

      ! Work out their positions
      ! Boundaries are the angles/times at the end of the chunks
      nullify ( boundaries )
      call Allocate_test ( boundaries, noChunks, 'boundaries', ModuleName )
      ! When we allow prior overlaps, the first chunk
      ! sometimes has 1 too many MAFs unless we take extra care
      if ( mexp1 == m1 ) then
        do chunk = 1, noChunks
          boundaries(chunk) = homeV + &
            & ( chunk - noChunksBelowHome ) * ChunkDivideConfig%maxLength
        end do
      else
        boundaries(1) = homeV + &
          & ( 1 - noChunksBelowHome ) * (ChunkDivideConfig%maxLength - 1)
        do chunk = 2, noChunks
          boundaries(chunk) = boundaries(chunk-1) + &
            & ChunkDivideConfig%maxLength
        end do
      endif
      boundaries = min ( boundaries, maxV )
      boundaries = max ( boundaries, minV )
      if ( ChunkDivideConfig%maxLengthFamily == PHYQ_Angle ) then
        chunks(1)%phiStart = homeV
        chunks%phiEnd = boundaries
        chunks(2:noChunks)%phiStart = chunks(1:noChunks-1)%phiEnd
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Regular HGrids will be unable to exploit chunk Geodetic Range' )
      endif

      ! Do some dumping
      if ( swLevel > -1 ) then
        call output ( ' minV: ' )
        call output ( minV  )
        call output ( ' maxV: ' )
        call output ( maxV , advance='yes' )
        call output ( ' homeV: ' )
        call output ( homeV , advance='yes' )
        call output ( ' noMAFSBelowHome: ' )
        call output ( noMAFsBelowHome, advance='yes' )
        call output ( ' noMAFSAtOrAboveHome: ' )
        call output ( noMAFsAtOrAboveHome, advance='yes' )
        call output ( ' NoChunks: ' )
        call output ( NoChunks, advance='yes' )
        call output ( ' noChunksBelowHome: ' )
        call output ( noChunksBelowHome , advance='yes' )
        call dump ( boundaries , 'boundaries' )
        call dump ( chunks%phiStart , 'phiStarts' )
        call dump ( chunks%phiEnd ,   'phiEnds' )
        call dump ( field, 'field' )
      end if

      call Hunt ( field, boundaries, chunks%lastMAFIndex, start=m1, &
        & allowTopValue=.true., nearest=.true. )
      call Deallocate_test ( boundaries, 'boundaries', ModuleName )
    end if

    ! Now deduce the chunk starts from the ends of their predecessors
    if ( noChunks > 1 ) &
      & chunks(2:noChunks)%firstMAFIndex = &
      & chunks(1:noChunks-1)%lastMAFIndex + 1
    chunks(1)%firstMAFIndex = m1

    ! Now bound the chunks to be within the processing range
    chunks%firstMAFIndex = min ( max ( chunks%firstMAFIndex, m1 ), m2 )
    chunks%lastMAFIndex = min ( max ( chunks%lastMAFIndex, m1 ), m2 )
    
    ! Each chunk remembers its own time range
    ! in case we readGriddedDatat during the loop of chunks
    ! chunks%Startime = taiTime%dpField(1,1,chunks%firstMAFIndex)
    ! chunks%Endtime =  taiTime%dpField(1,1,chunks%lastMAFIndex)

    ! Now offset these to the index in the file not the array
    ! chunks%firstMAFIndex = chunks%firstMAFIndex + mafRange(1) - 1
    ! chunks%lastMAFIndex = chunks%lastMAFIndex + mafRange(1) - 1
    chunks%firstMAFIndex = chunks%firstMAFIndex - 1
    chunks%lastMAFIndex = chunks%lastMAFIndex - 1

    ! If at this point the last two chunks end in the same place, this is
    ! a subtle defect in our chunking algorithm, lets avoid it
    if ( noChunks > 1 ) then
      if ( chunks(noChunks-1)%lastMAFIndex == chunks(noChunks)%lastMAFIndex ) then
        call DeleteChunk ( chunks, noChunks )
        noChunks = noChunks - 1
      end if
    end if

    ! Do some dumping
    if ( swLevel > -1 ) then
      call dump ( chunks%lastMAFIndex , 'chunks%lastMAFIndex' )
      call dump ( chunks%firstMAFIndex , 'chunks%firstMAFIndex' )
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Think about overlaps
    nullify ( newFirstMAFs, newLastMAFs )
    call Allocate_test ( newFirstMAFs, noChunks, 'newFirstMAFs', ModuleName )
    call Allocate_test ( newLastMAFs, noChunks, 'newLastMAFs', ModuleName )
    ! We could split things out to deal with mixed unit, but lets make life easier
    ! for ourselves.  ChunkDivideL2PC has already insisted that
    ! lowerOverlapFamily == upperOverlapFamily.
    if ( ChunkDivideConfig%lowerOverlapFamily == PHYQ_MAFs ) then
      newFirstMAFs = max(chunks%firstMAFIndex - nint(ChunkDivideConfig%lowerOverlap), &
        & m1 - 1 )
      newLastMAFs = min(chunks%lastMAFIndex + nint(ChunkDivideConfig%upperOverlap), &
        & m2 - 1 )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The bit of code that deals with non-MAF overlaps is probably broken' )
    end if
    chunks%noMAFsLowerOverlap = chunks%firstMAFIndex - newFirstMAFs
    chunks%noMAFsUpperOverlap = newLastMAFs - chunks%lastMAFIndex
    chunks%firstMAFIndex = newFirstMAFs
    chunks%lastMAFIndex = newLastMAFs
    if ( swLevel > -1 ) then
      call dump ( newFirstMAFs , 'newFirstMAFs' )
      call dump ( newLastMAFs , 'newLastMAFs' )
      call dump ( chunks%noMAFsLowerOverlap , 'chunks%noMAFsLowerOverlap' )
      call dump ( chunks%noMAFsUpperOverlap , 'chunks%noMAFsUpperOverlap' )
    endif
    call Deallocate_test ( newFirstMAFs, 'newFirstMAFs', ModuleName )
    call Deallocate_test ( newLastMAFs, 'newLastMAFs', ModuleName )

    ! Delete any zero length or all overlapped chunks
    call PruneChunks ( chunks )
    if ( .not. associated(chunks) ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No chunks remaining after we pruned them for bad radiances; bad day?' )
    else if ( size(chunks) < 1 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No chunks remaining after we pruned them for bad radiances; bad day?' )
    end if

    if ( swLevel > -1 ) then
      call output ( 'Before dealing with obstructions', advance='yes' )
      call Dump ( chunks )
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Now think about the obstructions
    if ( associated(obstructions) ) &
      & call DealWithObstructions ( chunks, obstructions )

    ! Delete any zero length or all overlapped chunks
    call PruneChunks ( chunks )
    if ( .not. associated(chunks) ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No chunks remaining after we pruned them for zero length; bad day?' )
    else if ( size(chunks) < 1 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No chunks remaining after we pruned them for zero length; bad day?' )
    end if

!     ! Forcibly zero out number of lower (upper) overlaps on 1st (last) chunks
    noChunks = size ( chunks )
!     chunks(1)%noMAFsLowerOverlap = 0
!     chunks(noChunks)%noMAFsUpperOverlap = 0

    ! Add lower overlap to first chunk if allowed
    if ( mexp1 < m1 ) then
      chunks(1)%firstMAFIndex = chunks(1)%firstMAFIndex + mexp1 - m1
      chunks(1)%noMAFsLowerOverlap = m1 - mexp1
    endif

    ! Add upper overlap to last chunk if allowed
    if ( mexp2 > m2 ) then
      chunks(noChunks)%lastMAFIndex = chunks(noChunks)%lastMAFIndex + mexp2 - m2
      chunks(noChunks)%noMAFsUpperOverlap = mexp2 - m2
    endif

    ! Now we assign phiStarts and phiEnds to each chunk
    do chunk = 1, noChunks
      m1 = max( chunks(chunk)%firstMAFIndex + 1, mafRange%L2Cover(1) + 1 )
      m2 = min( chunks(chunk)%lastMAFIndex + 1,  mafRange%L2Cover(2) + 1 )
      chunks(chunk)%phiStart = tpGeodAngle%dpField(1,1,m1)
      chunks(chunk)%phiEnd = tpGeodAngle%dpField(1,1,m2)
      chunks(chunk)%StartTime = taiTime%dpField(1,1,m1)
      chunks(chunk)%EndTime   = taiTime%dpField(1,1,m2)
    enddo
    if ( swLevel > -1 ) then
      call output ( 'After dealing with obstructions', advance='no' )
      call output ( ', poss. overlaps outside proc. range', advance='yes' )
      call Dump ( chunks )
      call dump ( chunks%phiStart  , 'phiStarts' )
      call dump ( chunks%phiEnd    , 'phiEnds  ' )
      call dump ( chunks%StartTime , 'StartTime' )
      call dump ( chunks%EndTime   , 'EndTime  ' )
    endif
    ! Tidy up
    call DeallocateL1BData ( tpGeodAngle )
    call DeallocateL1BData ( taiTime )

    if ( present(root) ) call trace_end ( &
      & 'ChunkDivide_orbital', cond=toggle(gen) .and. levels(gen) > 2 )

  end subroutine ChunkDivide_Orbital

  ! --------------------------------------------------- SortChunks -----
  subroutine SortChunks ( chunks )
    ! Sort the chunks into order of increasing firstMAFIndex
    type (MLSChunk_T), dimension(:), intent(inout) :: CHUNKS

    ! Local variables
    type (MLSChunk_T) :: TEMP
    integer :: I                        ! Loop counters
    integer, dimension(1) :: TOSWAP     ! Index

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    do i = 1, size(chunks) - 1
      toSwap = minloc ( chunks(i:)%firstMAFIndex + &
        & chunks(i:)%noMAFsLowerOverlap ) + (/ i-1 /)
      if ( toSwap(1) /= i ) then
        temp = chunks(i)
        chunks(i) = chunks(toSwap(1))
        chunks(toSwap(1)) = temp
      end if
    end do
  end subroutine SortChunks

  ! -------------------------------------------------- PruneChunks -----
  subroutine PruneChunks ( chunks )
     type (MLSChunk_T), dimension(:), pointer :: CHUNKS
     integer :: CHUNK
     ! Now delete chunks that either:
     !  1 - Are nothing but overlap
     !  2 - Have <=0 MAFs
     pruneChunksLoop: do
        chunk = FindFirst ( &
           & ( chunks%noMAFsLowerOverlap + chunks%noMAFsUpperOverlap ) >= &
           & ( chunks%lastMAFIndex - chunks%firstMAFIndex + 1 ) &
           & .or.&
           & ( chunks%firstMAFIndex > chunks%lastMAFIndex ) )
        if ( chunk == 0 ) exit pruneChunksLoop
        call DeleteChunk ( chunks, chunk )
     end do pruneChunksLoop
  end subroutine PruneChunks

  ! -------------------------------------------------- DeleteChunk -----
  subroutine DeleteChunk ( chunks, index )
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    ! Dummy arguments
    type (MLSChunk_T), pointer, dimension(:) :: CHUNKS
    integer, intent(in) :: INDEX

    ! Local variables
    type (MLSChunk_T), pointer, dimension(:) :: TEMP
    integer(c_intptr_t) :: Addr ! For tracing
    integer :: S                ! Size in bytes of object to deallocate
    integer :: STATUS           ! From allocate

    ! Executable code
    swlevel = switchDetail(switches, 'chu' )
    allocate ( temp ( size(chunks) - 1 ), stat=status )
    addr = 0
    if ( status == 0 .and. size(chunks) > 1 ) addr = transfer(c_loc(temp(1)), addr)
    call test_allocate ( status, ModuleName, 'temp', &
      & uBounds = size(chunks) - 1, elementSize = storage_size(temp) / 8, &
      & address=addr )

    if ( index > 1 ) temp(1:index-1) = chunks(1:index-1)
    if ( index < size(chunks) .and. size(chunks) > 1 ) &
       & temp(index:) = chunks(index+1:)

    s = size(chunks) * storage_size(chunks) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(chunks(1)), addr)
    deallocate ( chunks, stat=status )
    call test_deallocate ( status, moduleName, 'chunks', s, address=addr )

    chunks => temp

  end subroutine DeleteChunk

  !----------------------------------- Add obstruction to database -----
  subroutine AddObstructionToDatabase ( database, item )

     use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

     ! Dummy arguments
     type (Obstruction_T), dimension(:), pointer :: DATABASE
     type (Obstruction_T), intent(in) :: ITEM

     ! Local variables
     type (Obstruction_T), dimension(:), pointer :: TEMPDATABASE

     include "addItemToDatabase.f9h"

  end subroutine AddObstructionToDatabase

  !----------------------------------------- Add chunk to database -----
  subroutine AddChunkToDatabase ( database, item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (MLSChunk_T), dimension(:), pointer :: DATABASE
    type (MLSChunk_T), intent(in) :: ITEM

    ! Local variables
    type (MLSChunk_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"

  end subroutine AddChunkToDatabase

  ! --------------------------------------------  CFM_ChunkDivide  -----
  type(MLSChunk_T) function CFM_ChunkDivide (processingRange, &
     filedatabase, config) result(chunk)

     type (TAI93_Range_T), intent(in) :: processingRange
     type (MLSFile_T), dimension(:), pointer :: filedatabase
     type(ChunkDivideConfig_T), intent(in) :: config

     type (MAFRange_T) :: MAFRange
     type(MLSChunk_T), dimension(:), pointer :: chunks => null()

     ! Executables
     chunkDivideConfig = config

     nullify(chunkDivideConfig%criticalSignals)   ! Just for Sun's compiler
     nullify(obstructions)

     ! not sure how to deal with obstruction, ignore it
     if (chunkDivideConfig%method /= l_fixed) then
        call SurveyL1BData ( processingRange, filedatabase, mafRange)
     end if

     if (chunkdivideconfig%method == l_orbital) then
        call ChunkDivide_Orbital (mafRange, filedatabase, chunks)
        chunk = chunks(1)
     else
        print *, "chunk divide method is unsupported"
     end if
  end function CFM_ChunkDivide

  !------------------------------------------------  what_options  -----
  function what_options( clean, transpose, trim ) result( options )
    use MLSStrings, only: Trim_Safe
    logical, optional, intent(in) :: clean
    logical, optional, intent(in) :: transpose
    logical, optional, intent(in) :: trim
    character(len=8) :: options
    options = ' '
    if ( present(clean) ) then
      if ( clean ) options = trim_safe(options) // 'c'
    endif
    if ( present(transpose) ) then
      if ( transpose ) options = trim_safe(options) // 'p'
    endif
    if ( present(trim) ) then
      if ( trim ) options = trim_safe(options) // 't'
    endif
  end function what_options

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ChunkDivide_m

! $Log$
! Revision 2.139  2024/08/08 20:42:22  pwagner
! Fixed two bugs in handling level 1 obstructions
!
! Revision 2.138  2020/07/22 22:50:04  pwagner
! Added chunks%Start, EndTimes to Dumps ending ChunkDivide_Orbital
!
! Revision 2.137  2020/07/09 23:52:25  pwagner
! Added Start,EndTime components to MLSChunk_T
!
! Revision 2.136  2019/10/16 20:57:37  pwagner
! Fixed a strange bug involving get_boolean
!
! Revision 2.135  2018/10/05 20:42:21  pwagner
! Improved appearance of how we Dump Critical Signals
!
! Revision 2.134  2018/08/03 23:43:42  vsnyder
! Add tracing to several routines.  Some cannonball polishing
!
! Revision 2.133  2018/06/26 00:12:32  pwagner
! Avoids accessing index=1 of 0-sized temp in DeleteChunk
!
! Revision 2.132  2018/06/22 23:22:08  pwagner
! Say why we stop if none remain after 2nd PruneChunks
!
! Revision 2.131  2018/05/14 23:23:58  vsnyder
! Remove tab formatting (again) to eliminate compiler warnings
!
! Revision 2.130  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.129  2018/03/05 19:51:56  pwagner
! Now uses Dump_Signals
!
! Revision 2.128  2018/02/23 22:14:26  mmadatya
! Updated for polygon for ASMLS
!
! Revision 2.127  2018/02/09 00:58:27  pwagner
! Removed unused variables; wider use of outputNamedValue
!
! Revision 2.126  2018/01/03 01:15:22  pwagner
! Prints a little more debugging info
!
! Revision 2.125  2017/09/15 22:15:39  pwagner
! Correct bugs in evaluating excludeOverlap fields
!
! Revision 2.124  2017/09/14 23:19:36  pwagner
! Fixed some errors in ChunkDivide_Polygon
!
! Revision 2.123  2017/09/14 18:36:13  vsnyder
! Remove tab formatting to eliminate compiler warnings
!
! Revision 2.122  2017/03/06 19:55:48  pwagner
! Corrected remaining mistakes in evaluating logical values
!
! Revision 2.121  2017/02/15 00:48:26  pwagner
! Repaired bugs in processing Boolean fields /exclude..
!
! Revision 2.120  2017/02/10 21:57:47  pwagner
! Added the polygon method for ChunkDivide; NoChunks an optional parameter for fixed method
!
! Revision 2.119  2015/10/03 00:31:01  pwagner
! Will now show number of chunks
!
! Revision 2.118  2015/06/24 18:01:56  pwagner
! Halt with useful error message if no radiances, instead of vanishing in puff of smoke
!
! Revision 2.117  2015/03/28 02:19:01  vsnyder
! Added shallow destruction to DestroyChunkDatabase.  Got IsMonotonic from
! Monotone instead of MLSFillValues.
! Added stuff to trace allocate/deallocate addresses.
!
! Revision 2.116  2015/02/27 23:18:47  pwagner
! Require -Schu1 or greater to dump chunks
!
! Revision 2.115  2014/09/05 00:36:15  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Add some
! tracing.
!
! Revision 2.114  2014/08/07 22:47:12  vsnyder
! Delete local declaration of Obstructions from AddChunkToDatabase so that
! when it's nullified the one at module scope gets nullified.  It was
! not otherwise used in AddChunkToDatabase.
!
! Revision 2.113  2014/03/07 19:20:05  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.112  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.111  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.110  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.109  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.108  2013/10/09 23:40:34  vsnyder
! Add Evaluate_Variable
!
! Revision 2.107  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.106  2013/09/21 00:39:43  vsnyder
! Repair some tree examination errors.  Move ChunkDivideConfig_t,
! ChunkDivideConfig, and some dumpers to ChunkDivideConfig_m.  Call
! DumpCommand.
!
! Revision 2.105  2013/08/30 02:45:34  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.104  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.103  2013/06/28 19:02:58  pwagner
! criticalSignals may be used in place of criticalModules; it will be effective everywhere in marking obstructions
!
! Revision 2.102  2013/06/21 17:37:38  pwagner
! /crashIfPhiNotMono flag added to ChunkDivide config; default is to just warn; removed -Snmono switch
!
! Revision 2.101  2013/06/20 17:56:16  pwagner
! -Snmono lets us warn, not quit if phi not monotonic
!
! Revision 2.100  2013/05/21 22:50:24  pwagner
! Workaround for ifort13 bug
!
! Revision 2.99  2012/06/21 00:41:08  pwagner
! Added phi start and end to be used someday by HGrid
!
! Revision 2.98  2012/03/14 16:56:26  pwagner
! NAG-debug required this change
!
! Revision 2.97  2011/11/30 21:33:13  pwagner
! Quits with message if GeodAngle not monotonic
!
! Revision 2.96  2011/06/29 21:50:20  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.95  2011/02/05 01:41:44  pwagner
! Define and use consistently swLevel to control dump verboseness
!
! Revision 2.94  2010/04/20 17:32:23  honghanh
! Abandon attempt to remove the requirement to have a maxMafLength
! in ChunkDivide_Orbital
!
! Revision 2.92  2010/03/23 23:53:08  honghanh
! Move most subroutine inside ChunkDivide out, including
! ChunkDivide_Orbital, create CFM_ChunkDivide, which call
! ChunkDivide_Orbital
!
! Revision 2.91  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.90  2009/06/16 17:42:59  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.89  2009/04/02 18:09:26  pwagner
! Fixed bug where status might be undefined
!
! Revision 2.88  2009/04/01 23:34:26  pwagner
! By default saves obstructions db to written to l2aux file
!
! Revision 2.87  2008/07/12 00:12:07  pwagner
! dump_criticalsignals not generic to appease sun studio 12
!
! Revision 2.86  2008/05/28 21:51:43  pwagner
! May choose critical band(s)
!
! Revision 2.85  2007/12/14 01:55:28  pwagner
! Must delete obstructions observing proper order for indexing
!
! Revision 2.84  2007/11/01 23:30:48  pwagner
! Should permit us to make sids files, i.e. omit l1brad
!
! Revision 2.83  2007/10/24 00:14:58  pwagner
! Removed unused declarations
!
! Revision 2.82  2007/09/06 22:31:50  pwagner
! Raise switch threshold for dumpings signals
!
! Revision 2.81  2007/06/08 22:00:23  vsnyder
! Quit gracefully if no L1BOA files
!
! Revision 2.80  2007/03/23 00:16:57  pwagner
! Prevent crashing while printing extra-long signals names
!
! Revision 2.79  2007/02/09 01:10:00  pwagner
! Fixed bug where criticalSignals not nullified
!
! Revision 2.78  2007/02/06 23:13:40  pwagner
! Now can chooseCriticalSignals based on criticalModule
!
! Revision 2.77  2007/01/25 19:04:59  pwagner
! Warns if DACsDeconvolved attribute not found
!
! Revision 2.76  2006/06/20 00:12:30  pwagner
! Improved two printed formats
!
! Revision 2.75  2006/05/09 23:41:41  pwagner
! Warn if not assured DACS deconvolution performed by level 1
!
! Revision 2.74  2006/04/20 23:22:54  pwagner
! Show both kinds of allowed extra-range overlaps
!
! Revision 2.73  2006/04/10 23:45:18  pwagner
! Reset defaults in read_apriori, not ChunkDivide
!
! Revision 2.72  2006/04/03 20:26:08  pwagner
! More verbose notice when using l2gp a priori exclude outside overlaps
!
! Revision 2.71  2006/03/17 21:57:43  pwagner
! Adjust default behavior to exclude overlaps outside processing if sids run
!
! Revision 2.70  2006/03/17 00:06:31  pwagner
! Change default to allowing overlaps outside processingRange
!
! Revision 2.69  2006/02/10 21:19:24  pwagner
! dumps may go to special dumpfile
!
! Revision 2.68  2006/02/07 00:55:47  pwagner
! Now allows overlaps after data end time
!
! Revision 2.67  2006/01/26 00:34:50  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.66  2005/12/16 00:06:51  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.65  2005/10/22 00:43:43  pwagner
! Gets DACSDeconvolved attribute from l1b file if there
!
! Revision 2.64  2005/09/21 23:25:42  pwagner
! Obstructions public now; optionally saved; should add deallocate in tree_walker
!
! Revision 2.63  2005/09/14 00:10:37  pwagner
! ChunkDivideConfig public so allowPriorOverlaps visible
!
! Revision 2.62  2005/08/09 00:02:09  pwagner
! hdfVersion not left undefined in ANY_GOOD_SIGNALDATA
!
! Revision 2.61  2005/06/03 02:02:17  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.60  2005/06/01 17:39:26  pwagner
! Dont read L1bFile if unassocated
!
! Revision 2.59  2005/05/31 17:51:16  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.58  2005/05/27 19:53:43  vsnyder
! Produce meaningful error message instead of crash for empty section
!
! Revision 2.57  2004/12/14 21:43:33  pwagner
! Repaired bug in reading mafRange rather than all mafs
!
! Revision 2.56  2004/11/03 17:19:09  livesey
! Bug fix in case where only one chunk
!
! Revision 2.55  2004/10/05 23:09:54  pwagner
! Can handle dropped MAFs, maneuvers that disrupt monotonic geodAngle
!
! Revision 2.54  2004/08/23 22:00:39  pwagner
! Made most readl1bData dontpad=.true.
!
! Revision 2.53  2004/08/16 17:10:04  pwagner
! Passes dontPad option to readL1BData
!
! Revision 2.52  2004/08/09 21:43:10  livesey
! Bug fixes and added the maxOrbY argument.
!
! Revision 2.51  2004/08/04 23:19:57  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.50  2004/08/02 23:40:29  livesey
! Bug fixes in the orbital case when chunk length is expressed in MAFs
!
! Revision 2.49  2004/07/31 19:58:28  livesey
! Various bug fixes and clean ups in the light of real data.
!
! Revision 2.47  2004/06/10 00:58:44  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.46  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.45  2003/08/27 20:08:04  livesey
! Removed print statements
!
! Revision 2.44  2003/08/26 19:49:22  livesey
! Removed another print statement.
!
! Revision 2.43  2003/08/26 19:43:42  livesey
! Commented out a print statement
!
! Revision 2.42  2003/08/26 18:05:14  livesey
! Various fixes to the obstruction stuff.  More to come probably.
!
! Revision 2.41  2003/08/25 23:43:08  pwagner
! Added ReduceChunkDatabase
!
! Revision 2.40  2003/08/21 22:51:18  livesey
! Removed yet more print statements!
!
! Revision 2.39  2003/08/21 22:49:22  livesey
! Removed more print statements etc.
!
! Revision 2.38  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.37  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.36  2003/06/09 22:55:23  pwagner
! Improved criticalSignals; some bug fixes
!
! Revision 2.35  2003/06/05 23:37:14  pwagner
! 1st version of criticalSignals in Chunk Divide Orbital
!
! Revision 2.34  2003/05/17 00:06:06  pwagner
! Wont check sidebands for missing radiances
!
! Revision 2.33  2003/05/09 16:43:05  pwagner
! Speedup of L1B radiance check
!
! Revision 2.32  2003/05/07 23:43:05  pwagner
! Optionally may skipL1BCheck
!
! Revision 2.31  2003/04/30 22:06:04  pwagner
! Shouldnt check for good signals if there arent any
!
! Revision 2.30  2003/04/28 23:07:00  pwagner
! Fleshed out notel1brad_changes; not yet tested where needed
!
! Revision 2.29  2003/01/06 20:13:09  livesey
! New split upper/lower overlaps
!
! Revision 2.28  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.27  2002/12/06 01:08:33  pwagner
! Less likely to bomb on single chunks
!
! Revision 2.26  2002/11/13 01:03:11  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.25  2002/11/06 00:21:16  pwagner
! Fixed non-zero starttime prob; some extra checks, printing
!
! Revision 2.24  2002/10/07 23:49:49  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.23  2002/10/07 18:00:10  livesey
! Added call to trace_end, whoops!
!
! Revision 2.22  2002/08/22 01:23:52  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.21  2002/08/04 15:59:45  mjf
! New method "PE" for a single chunk centred on (or near) a given phi.
!
! Revision 2.20  2002/05/24 20:57:34  livesey
! More revisions associated with being able to have a specific number of chunks
!
! Revision 2.19  2002/05/24 16:47:30  livesey
! Allowed user to request specific number of chunks for orbital
!
! Revision 2.18  2001/12/17 23:09:18  livesey
! Now deletes chunks that are nothing but overlap
!
! Revision 2.17  2001/12/16 00:56:43  livesey
! Added dump option
!
! Revision 2.16  2001/11/20 00:10:29  livesey
! Another bug fix
!
! Revision 2.15  2001/11/19 23:53:50  livesey
! Works better, fixed bug with cases with all missing data, also handles
! ends of processing range properly.
!
! Revision 2.14  2001/11/19 23:33:45  livesey
! Interim version, some bugs to track down
!
! Revision 2.13  2001/11/16 20:38:42  livesey
! Nullified some pointers I missed
!
! Revision 2.12  2001/11/15 17:43:46  livesey
! Tidied up, more comments etc.
!
! Revision 2.11  2001/11/14 22:33:40  livesey
! This version seems pretty good.
!
! Revision 2.10  2001/11/14 01:49:12  livesey
! Improvements, getting closer
!
! Revision 2.9  2001/11/12 21:15:34  livesey
! More bug fixes
!
! Revision 2.8  2001/11/12 20:31:23  livesey
! Bug fix check in
!
! Revision 2.7  2001/11/10 01:09:41  livesey
! Bug fixes
!
! Revision 2.6  2001/11/10 01:01:03  livesey
! Just tidied up
!
! Revision 2.5  2001/11/10 00:56:24  livesey
! OK, this is getting close to being ready for testing.
!
! Revision 2.4  2001/11/10 00:03:44  livesey
! More work, still got issues!
!
! Revision 2.3  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.2  2001/11/09 06:34:39  livesey
! Minor bug fix, and added Log stuff
!
