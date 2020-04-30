! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TREE_WALKER

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions of the MLS L2 processing in the order indicated.

  implicit none
  private

  public :: WALK_TREE_TO_DO_MLS_L2
  
  logical, parameter :: ComplainIfSkippedEveryChunk = .true.
  logical, parameter :: MustCheckForCorruptFileDatabase = .false.

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------  WALK_TREE_TO_DO_MLS_L2  -----
  subroutine WALK_TREE_TO_DO_MLS_L2 ( ROOT, ERROR_FLAG, FIRST_SECTION, &
    & COUNTCHUNKS, FILEDATABASE )

    use Algebra_M, only: Algebra
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, FinalMemoryReport
    use AntennaPatterns_M, only: Destroy_Ant_Patterns_Database
    use ChunkDivide_M, only: ChunkDivide, DestroyChunkDatabase
    use Chunks_M, only: Dump, MLSChunk_T
    use Construct, only: MLSL2Construct, MLSL2Deconstruct, &
      & ConstructMIFGeolocation, DestroyMIFGeolocation
    use DirectWrite_M, only: DirectData_T, DestroydirectDatabase
    use Dump_0, only: Dump
    use EmpiricalGeometry, only: ForgetOptimumLon0
    use FGrid, only: Fgrid_T, DestroyFGridDatabase, Dump
    use Fill, only: MLSL2Fill
    use ForwardModelConfig, only: ForwardModelConfig_T, &
      & DestroyFwmConfigDatabase, &
      & StripForwardModelConfigDatabase
    use ForwardModelSupport, only: PrintForwardModelTiming
    use Global_Settings, only: Set_Global_Settings
    use GriddedData, only: GriddedData_T, DestroyGriddedDataDatabase, Dump
    use HessianModule_1, only: DestroyHessianDatabase, Hessian_T
    use HGridsDatabase, only: HGrids_T
    use HGrid, only: ComputeAllHGridOffsets, DestroyHGridGeoLocations
    use HighOutput, only: BeVerbose, GetStamp, HeadLine, OutputNamedValue, &
      & SetStamp
    use Init_Tables_Module, only: L_Chisqchan, L_Chisqmmaf, L_Chisqmmif, &
      & Section_First, Section_Last, &
      & Z_Algebra, Z_Chunkdivide, Z_Construct, Z_Fill, Z_Globalsettings, &
      & Z_Join, Z_Mergegrids, Z_MLSsignals, Z_Output, Z_Readapriori, &
      & Z_Retrieve, Z_Spectroscopy
    use Intrinsic, only: Section_Indices
    use Join, only: MLSL2Join
    use L1BData, only: CheckForCorruptFileDatabase
    use L2AUXData, only: DestroyL2AUXDatabase, L2AUXData_T, Dump
    use L2FWMParallel, only: L2FWMSlaveTask, LaunchFWMSlaves
    use L2GPData, only: DestroyL2GPDatabase, L2GPData_T, Dump
    use L2Parallel, only: GetChunkInfoFromMaster, L2MasterTask
    use L2Parinfo, only: Parallel, CloseParallel
    use L2PC_M, only: DestroyL2PCDatabase, DestroyBinSelectorDatabase
    use L2PCBins_M, only: FlushLockedBins
    use MatrixModule_1, only: DestroyMatrixDatabase, Matrix_Database_T
    use MergeGridsModule, only: MergeGrids
    use MLSCommon, only: TAI93_Range_T, MLSFile_T
    use MLSL2Options, only: Aura_L1BFiles, &
      & CheckPaths, L2Options, ExitToNextChunk, &
      & L2CFNode, MLSL2Message, Need_L1BFiles, PhasesToSkip, &
      & SectionsToSkip, SkipDirectWrites, SkipDirectWritesOriginal, &
      & SlavesCleanUpSelves, SpecialDumpFile, StopAfterSection, &
      & Toolkit
    use MLSMessageModule, only: MLSMSG_Allocate, MLSMSG_Info, &
      & MLSMSG_Error, SummarizeWarnings
    use MLSPCF2, only: MLSPCF_Spectroscopy_End
    use MLSFinds, only: FindFirst
    use MLSSignals_M, only: Bands, DestroyBandDatabase, DestroyModuleDatabase, &
      & DestroyRadiometerDatabase, DestroySignalDatabase, &
      & DestroySpectrometerTypeDatabase, IsSpacecraftAura, &
      & MLSSignals, Modules, Radiometers, &
      & Signals, SpectrometerTypes
    use MLSStringLists, only: ExpandStringRange, IsInList, SwitchDetail
    use MLSStrings, only: Lowercase
    use MLSL2Timings, only: Add_To_Section_Timing, RestartTimings
    use Next_Tree_Node_M, only: Dump, &
      & Next_Tree_Node, Next_Tree_Node_State
    use Open_Init, only: OpenAndInitialize
    use OutputAndClose, only: Output_Close
    use Output_M, only: Output, ResumeOutput, RevertOutput, SwitchOutput
    use PointingGrid_M, only: Destroy_Pointing_Grid_Database
    use QuantityTemplates, only: QuantityTemplate_T ! , &
!       & DestroyQuantityTemplateDatabase
    use ReadApriori, only: Read_Apriori
    use RetrievalModule, only: Retrieve
    use SpectroscopyCatalog_M, only: Destroy_Line_Database, &
      & Destroy_Spectcat_Database, Spectroscopy
    use String_Table, only: Get_String
    use Time_M, only: SayTime, Time_Now
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Subtree
    use VectorsModule, only: DestroyVectorDatabase, Dump_Vectors, &
      & Vector_T, VectorTemplate_T
    use VGridsDatabase, only: DestroyVGridDatabase, VGrids

    integer, intent(in) ::     ROOT        ! Root of the abstract syntax tree
    integer, intent(out) ::    ERROR_FLAG  ! Nonzero means failure
    integer, intent(in) ::     FIRST_SECTION! Index of son of root of first n_cf
    logical, intent(in) ::     COUNTCHUNKS ! Just count, print, and quit

    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    ! Internal variables
    logical ::                                   canWriteL2PC
    integer ::                                   chunkNo ! Index of Chunks
    type (MLSChunk_T), dimension(:), pointer ::  Chunks  ! of data
    logical, dimension(:), pointer :: CHUNKSSKIPPED => null() ! Don't do these
    integer                                   :: details
    logical                                   :: doneWithChunks
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (FGrid_T), dimension(:), pointer ::     FGrids
    ! Forward model configurations:
    integer ::                                   FIRSTCHUNK ! For chunk loop
    type (ForwardModelConfig_T), dimension(:), &
      & pointer ::                               ForwardModelConfigDatabase
    type (GriddedData_T), dimension(:), &
      & pointer ::                               GriddedDataBase
    type (Hessian_T), dimension(:), pointer ::   Hessians
    type (HGrids_T), dimension(:), pointer ::    HGrids
    integer ::                                   I       ! Loop inductors
    integer, dimension(1) :: ICHUNKS
    integer ::                                   LASTCHUNK ! For chunk loop
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (L2GPData_T), dimension(:), pointer  :: L2GPDatabase
    type (Matrix_Database_T), dimension(:), &
      & pointer ::                               Matrices
    integer :: Me = -1                           ! String index for trace
    logical ::                                   now_stop
    type (TAI93_Range_T) ::                      ProcessingRange  ! Data processing range
    integer ::                                   section_index
    character(len=32) ::                         section_name
    logical ::                                   showTime
    logical, dimension(SECTION_FIRST:SECTION_LAST) :: skipSections
    integer ::                                   SON     ! Son of Root
    type(next_tree_node_state) ::                State, Save1, Save2 ! of tree traverser
    logical ::                                   STOPBEFORECHUNKLOOP
    real    ::                                   t1, tChunk
    character(len=24) ::                         textCode
    type (Vector_T), dimension(:), pointer ::    Vectors
    logical :: verbose
    logical :: verboser
    logical :: warnOnDestroy        ! Whether to warn before destroying dbs

    ! Arguments for Construct not declared above:
    type (QuantityTemplate_T), dimension(:), pointer :: MifGeolocation
    type (QuantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (VectorTemplate_T), dimension(:), pointer :: VectorTemplates

    ! Executable
    verbose = BeVerbose ( 'walker', -1 )
    verboser = BeVerbose ( 'walker', 0 )
    call trace_begin ( me, 'WALK_TREE_TO_DO_MLS_L2', &
      & subtree(first_section,root), index=first_section, cond=toggle(gen) )
    call time_now ( t1 )

    nullify ( chunks, forwardModelConfigDatabase, griddedDataBase, hessians, &
      & directDatabase, hGrids, l2auxDatabase, l2gpDatabase, matrices, mifGeolocation, &
      & qtyTemplates, vectorTemplates, fGrids, vGrids )

    ! Allocate Vectors with size zero so nobody has to check whether it's
    ! associated.  This allows to remove the pointer attribute in several
    ! places.
    allocate ( vectors(0), stat=error_flag )
    if ( error_flag /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // 'vectors' )
    allocate ( directDatabase(0), stat=error_flag )
    if ( error_flag /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // 'direct write files' )

    nullify ( chunksSkipped )
    warnOnDestroy = ( switchDetail(switches, 'destroy' ) > -1 )
    canWriteL2PC = .false.
    stopBeforeChunkLoop          = ( &
      & index('global_setting,chunk_divide', &
      & lowercase( trim(stopAfterSection) ) ) > 0 )
    stopBeforeChunkLoop = ( stopBeforeChunkLoop .and. stopAfterSection /= ' ' )
    skipSections = .false.
    do i = section_first, section_last
      call get_string ( section_indices(i), section_name, strip=.true. )
      skipSections(i) = isInList( sectionstoskip, section_name, '-fc' )
    end do
    call OpenAndInitialize ( processingRange, filedatabase )
    call add_to_section_timing ( 'open_init', t1, now_stop )
    if ( now_stop ) then
      call finishUp ( .true. )
      return
    end if

    ! -----------------------------------------------------
    ! ----------------------------------------------------- Loop over tree

    chunkNo = 0
    doneWithChunks = .false.
    ! Now loop over the sections in the tree
    do
      save1 = state
      son = next_tree_node(root,state,start=first_section)
      if ( son == 0 ) exit
      L2CFNODE = son
      section_index = decoration(subtree(1,son))
      call get_string ( section_indices(section_index), section_name, &
        & strip=.true. )
      if ( section_index > SECTION_FIRST - 1 .and. &
        & section_index < SECTION_LAST + 1 ) then
        if( skipSections(section_index) ) &
          & section_index = SECTION_FIRST - 1 ! skip
      end if
      ! First those involved in 'preprocessing'
      select case ( section_index )

        ! -------------------------------------------------------- Init sections
      case ( z_globalsettings )
        if ( verboser ) call Dump ( state, 'at globalsettings' )
        call set_global_settings ( son, forwardModelConfigDatabase, &
          & filedatabase, fGrids, l2gpDatabase, DirectDatabase, processingRange )
        call add_to_section_timing ( 'global_settings', t1, now_stop )
        if ( now_stop .and. .not. parallel%slave ) then
          call finishUp(.true.)
          return
        end if
        if ( verbose .and. MustCheckForCorruptFileDatabase ) then
          call output( 'Done with global settings', advance='yes' )
          call CheckForCorruptFileDatabase( filedatabase )
        endif
      case ( z_mlsSignals )
        if ( verboser ) call Dump ( state, 'at mlsSignals' )
        call MLSSignals ( son )
        if ( switchDetail(switches,'tps') > -1 ) then
          ! call test_parse_signals
          call MLSL2Message ( MLSMSG_Info, ModuleName, &
            & 'Go back and uncomment the previous line in tree_walker' )
        end if
        ! Here's one way for the l2cf to set Aura to .false.
        ! Put a line like
        !  SC: Module, Aura=false, /spaceCraft
        ! in your MLSSignals section
        AURA_L1BFILES = AURA_L1BFILES .and. IsSpaceCraftAura()
        call add_to_section_timing ( 'signals', t1, now_stop )
        if ( now_stop ) then
          call finishUp(.true.)
          return
        end if
      case ( z_spectroscopy )
        if ( verboser ) call Dump ( state, 'at spectroscopy' )
        call spectroscopy ( son, &
          & TOOLKIT, mlspcf_spectroscopy_end, fileDatabase )
        call add_to_section_timing ( 'spectroscopy', t1, now_stop )
        if ( now_stop ) then
          call finishUp(.true.)
          return
        end if
      case ( z_readapriori )
        if ( verboser ) call Dump ( state, 'at readapriori' )
        if ( .not. ( stopBeforeChunkLoop .or. parallel%master ) ) &
          & call read_apriori ( son , &
          & l2gpDatabase, l2auxDatabase, griddedDataBase, fileDataBase )
        call add_to_section_timing ( 'read_apriori', t1, now_stop )
        if ( now_stop ) then
          call finishUp(.true.)
          return
        end if
      case ( z_mergeGrids )
        if ( .not. &
          & ( stopBeforeChunkLoop .or. checkPaths .or. parallel%master ) &
          & ) call mergeGrids ( son, l2gpDatabase, l2auxDatabase, &
          & griddedDataBase, fileDataBase )

        ! ------------------------------------------------------- Chunk divide
        ! Chunk divide can be a special one, in slave mode, we just listen out
        ! for instructions.
        call add_to_section_timing ( 'merge_grids', t1, now_stop )
        if ( now_stop ) then
          call finishUp(.true.)
          return
        end if
      case ( z_chunkdivide )
        if ( verboser ) call Dump ( state, 'at chunkdivide' )
        if ( parallel%slave .and. .not. parallel%fwmParallel ) then
          call GetChunkInfoFromMaster ( chunks, chunkNo )
          firstChunk = chunkNo
          lastChunk = chunkNo
          parallel%ChunkNo = chunkNo
          L2Options%CurrentChunkNumber = chunkNo
        else
          if ( (.not. checkPaths .or. parallel%chunkRange /= '') .and. &
            & NEED_L1BFILES ) then
            if ( verbose .and. MustCheckForCorruptFileDatabase ) then
              call output( 'before chunk divide', advance='yes' )
              call CheckForCorruptFileDatabase( filedatabase )
            endif
            call ChunkDivide ( son, processingRange, filedatabase, chunks )
            if ( verbose .and. MustCheckForCorruptFileDatabase ) then
              call output( 'after chunk divide', advance='yes' )
              call CheckForCorruptFileDatabase( filedatabase )
            endif
            call ComputeAllHGridOffsets ( root, first_section, chunks, &
              & filedatabase, &
              & l2gpDatabase, processingRange )
            if ( verbose .and. MustCheckForCorruptFileDatabase ) then
              call output( 'after computing hgrid offsets', advance='yes' )
              call CheckForCorruptFileDatabase( filedatabase )
            endif
          else
            allocate ( chunks(1), stat=error_flag )
            if ( error_flag /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
              & MLSMSG_Allocate // 'chunks' )
            if ( .not. NEED_L1BFILES ) then
              call output( 'Creating artificial chunk', advance='yes' )
              if ( switchDetail(switches, 'chu' ) > -1 ) &
                & call dump( chunks(1) )
            end if
          end if
          firstChunk = 1
          lastChunk = size(chunks)
          if ( countChunks ) then
            error_flag = 0
            call output ( size(chunks) )
            return
          end if
        end if
        call add_to_section_timing ( 'chunk_divide', t1, now_stop )
        if ( now_stop ) then
          call finishUp(.true.)
          return
        end if
        if ( verbose ) &
          & call outputNamedValue ( 'First, last chunks', (/ firstChunk, lastChunk /) )
      case ( z_algebra )
        call algebra ( son, vectors, matrices, chunks(1), &
          & forwardModelConfigDatabase )
        ! ---------------------------------------------------- Chunk processing
        ! Now construct, fill, join and retrieve live inside the 'chunk loop'
      case ( z_construct, z_fill, z_join, z_retrieve )
        if ( doneWithChunks ) cycle
        if ( verboser ) call Dump ( state, 'at other' )
        ! Do special stuff in some parallel cases, or where there are
        ! no chunks.
        if ( .not. associated(chunks) ) then
          allocate(chunks(1), stat=error_flag)
          if ( error_flag /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
            & 'unable to allocate chunks' )
          firstChunk = 1
          lastChunk = size(chunks)
        end if
        if ( ( size(chunks) < 1 ) .or. &
          & ( parallel%master .and. .not. parallel%fwmParallel ) .or. &
          & ( parallel%slave .and. parallel%fwmParallel ) ) then
          if ( parallel%master .and. .not. parallel%fwmParallel ) then
            call outputNamedValue( 'Calling L2MasterTask with num chunks ', size(chunks) )
            call L2MasterTask ( chunks )
            doneWithChunks = .true.
            call add_to_section_timing ( 'master', t1, now_stop )
            cycle
          end if
          if ( parallel%slave .and. parallel%fwmParallel ) then
            call ExpandStringRange( parallel%chunkRange, iChunks )
            call ConstructMIFGeolocation ( mifGeolocation, filedatabase, &
              & chunks(iChunks(1)) ) 
            call L2FWMSlaveTask ( mifGeolocation )
          end if
          ! Sort out the timings
          select case ( section_index )
          case ( z_algebra )
            call add_to_section_timing ( 'algebra', t1, now_stop )
          case ( z_construct )
            call add_to_section_timing ( 'construct', t1, now_stop )
          case ( z_fill )
            call add_to_section_timing ( 'fill', t1, now_stop )
          case ( z_join )
            call add_to_section_timing ( 'join', t1, now_stop )
          case ( z_retrieve )
            call add_to_section_timing ( 'retrieve', t1, now_stop )
          case default
            ! This is one of the skipped sections
          end select
          ! print the timing for FullForwardModel, the following return
          ! if fmt1 or fmt2 is true
          if ( switchDetail(switches, 'fmt') > -1 .and. &
            & associated(forwardModelConfigDatabase)) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
          end if
          canWriteL2PC = .true.
        else
        ! Otherwise, this is the 'standard' work for these sections.
          call allocate_test( chunksSkipped, size(chunks), 'chunksSkipped', ModuleName )
          chunksSkipped = .false.
          if ( parallel%chunkRange /= '' ) then
            chunksSkipped = .true.
            call ExpandStringRange(trim(parallel%chunkRange), chunksSkipped, &
              & sense=.false.)
          end if  
          canWriteL2PC = ( count(.not. chunksSkipped) < 2 )
          if ( .not. associated(chunks) ) then
            allocate(chunks(1), stat=error_flag)
            if ( error_flag /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
              & 'unable to allocate chunks' )
            firstChunk = 1
            lastChunk = size(chunks)
          end if
          if ( FindFirst( .not. chunksSkipped(firstChunk:lastChunk) ) < 1 .and. &
            & COMPLAINIFSKIPPEDEVERYCHUNK ) &
              & call MLSL2Message ( MLSMSG_Error, ModuleName, &
              & 'We have skipped every chunk' )
          if ( verboser ) call Dump( save1, 'save1' )
          do chunkNo = firstChunk, lastChunk ! --------------------- Chunk loop
            state = save1 ! Back up so first repeated section is next
            call resumeOutput ! In case the last phase was  silent
            call restartTimings( 'flags' ) ! reset only the flags
            if ( chunksSkipped(chunkNo) ) cycle
            call time_now ( tChunk )
            if ( BeVerbose( (/ 'chu', 'pro' /), 0) .or. verbose ) then
              call output ( " ================ Starting processing for chunk " )
              call output ( chunkNo )
              call output ( " ================ ", advance='yes' )
            end if
            if ( parallel%master .and. parallel%fwmParallel ) then
              call LaunchFWMSlaves ( chunks ( chunkNo ) )
            elseif( .not. parallel%slave ) then
              L2Options%CurrentChunkNumber = chunkNo  ! Stored for dumping
            endif
            if ( verbose .and. MustCheckForCorruptFileDatabase ) then
              call output( 'Now in loop of chunks', advance='yes' )
              call CheckForCorruptFileDatabase( filedatabase )
            endif
            exitToNextChunk = .false.
            call FinalMemoryReport
subtrees:   do
              ! in the outer loop
              save2 = state ! before advancing to section below
              son = next_tree_node(root,state)
              if ( verboser ) call Dump ( state, 'inner state' )
              if ( exitToNextChunk ) &
                & call output ( &
                & 'Commanded to Skip remaining sections of this chunk', &
                & advance='yes' )
              if ( son == 0 .or. exitToNextChunk ) exit subtrees
              L2CFNODE = son
              section_index = decoration(subtree(1,son))
              call get_string ( section_indices(section_index), section_name, &
                & strip=.true. )
              ! Are we to skip this section? Reasons could be
              ! (1) It's in the list of sections to Skip
              ! (2) The phase is one we've been asked to skip
              if( skipSections(section_index) ) &
                & section_index = SECTION_FIRST - 1 ! skip
              if ( isInList( PhasesToSkip, trim(L2Options%CurrentPhaseName), '-fc' ) ) &
                & section_index = SECTION_FIRST - 2 ! skip
              if ( verboser ) call MLSL2Message ( MLSMSG_Info, ModuleName, &
                & 'Innermost loop ' // trim(section_name) )
              ! Start inner loop for one chunk at the current position
              select case ( section_index ) ! section index
              case ( z_algebra )
                call algebra ( son, vectors, matrices, chunks(chunkNo), forwardModelConfigDatabase )
              case ( z_construct )
                if ( .not. checkPaths ) &
                & call MLSL2Construct ( son, filedatabase, processingRange, &
                  & chunks(chunkNo), qtyTemplates, vectorTemplates, vectors, &
                  & fGrids, hGrids, l2gpDatabase, forwardModelConfigDatabase, &
                  & griddedDataBase, mifGeolocation )
                call add_to_section_timing ( 'construct', t1, now_stop )
                if ( associated(hGrids) ) then
                  if ( HGrids(1)%the_hGrid%noProfs < 1 ) then
                    call headLine ( 'No profiles, so skipping this chunk' )
                    exit subtrees
                  endif
                endif
              case ( z_fill )
                if ( .not. checkPaths ) then 
                  call MLSL2Fill ( son, filedatabase, griddedDataBase, &
                  & vectorTemplates, vectors, qtyTemplates, matrices, hessians, &
                  & l2gpDatabase, l2auxDatabase, forwardModelConfigDatabase, &
                  & chunks, chunkNo, hgrids )
                end if
                call add_to_section_timing ( 'fill', t1, now_stop )
              case ( z_join )
                call MLSL2Join ( son, vectors, l2gpDatabase, &
                  & l2auxDatabase, DirectDatabase, chunkNo, chunks, &
                  & forwardModelConfigDatabase, fileDatabase, HGrids, &
                  & matrices, Hessians )
                call add_to_section_timing ( 'join', t1, now_stop )
              case ( z_retrieve )
                if ( .not. checkPaths) &
                & call retrieve ( son, vectors, matrices, hessians, &
                  & forwardModelConfigDatabase, &
                  & chunks(chunkNo), fileDataBase )
                call add_to_section_timing ( 'retrieve', t1, now_stop )
              case ( z_output )
                exit subtrees
              case default
                ! exit subtrees
                ! This may simply be a skipped section
                if ( verbose ) then
                  if ( section_index == SECTION_FIRST - 1 ) then
                    call output ( 'Skipping this section', advance='yes' )
                  else
                    call output ( 'Skipping ' // trim(L2Options%CurrentPhaseName), advance='yes' )
                  endif
                endif
              end select
            end do subtrees

            if ( switchDetail(switches,'chi') > 0 .and. chunkNo > 1 ) then
              ! Dumps nothing after 1st chunk
            else if ( switchDetail(switches,'chi') > -1 ) then
              if ( specialDumpFile /= ' ' ) &
                & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
              call output('Here are our diagnostics for chunk ', advance='no')
              call output(chunkNo, advance='yes')
              call dump_vectors( vectors, details=1, &
              & quantityTypes = (/l_chisqchan, l_chisqmmaf, l_chisqmmif/) )
              if ( specialDumpFile /= ' ' ) &
                & call revertOutput
            end if

            ! Now, if we're dealing with more than one chunk destroy stuff
            ! Otherwise, we'll save them as we may need to output them as l2pc files.
            if ( .not. canWriteL2PC ) then
              if ( warnOnDestroy ) call output('About to MLSL2Deconstruct', advance='yes' )
              call MLSL2DeConstruct ( vectorTemplates, hGrids )
              if ( warnOnDestroy ) call output('About to destroy hessian db', advance='yes' )
              call DestroyHessianDatabase ( hessians )
              if ( warnOnDestroy ) call output('About to destroy matrix db', advance='yes' )
              call DestroyMatrixDatabase ( matrices )
              if ( warnOnDestroy ) call output('About to destroy vector db', advance='yes' )
              call DestroyVectorDatabase ( vectors )
            end if
            if ( warnOnDestroy ) call output('About to forget optimum Lon0', advance='yes' )
            call ForgetOptimumLon0
            ! print the timing for FullForwardModel
            ! fmt2: at each chunk, fmt1: at last chunk
            if ( switchDetail(switches, 'fmt') > 1 ) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
            end if
            if ( switchDetail(switches, 'fmt') == 1 .and. &
              & chunkNo == lastChunk ) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
            end if
            if ( warnOnDestroy ) call output('About to strip fwmconfig db', advance='yes' )
            call StripForwardModelConfigDatabase ( forwardModelConfigDatabase )
            if ( switchDetail(switches,'chu') > -1 ) then
              call sayTime( 'processing this chunk', t1=tChunk )
            end if
            if ( now_stop ) then
              call finishUp(.true.)
              return
            end if
            call DestroyMIFGeolocation( mifGeolocation ) 
          end do ! ---------------------------------- End of chunk loop
          doneWithChunks = .true.
          state = save2
          ! Continue the outer loop after the last section processed in
          ! the chunk loop
          ! Clear any locked l2pc bins out.
          if ( warnOnDestroy ) call output('About to flushl2pc bins', advance='yes' )
          call FlushLockedBins
          ! Done with the chunksSkipped array
          call deallocate_test( chunksSkipped, 'chunksSkipped', ModuleName )
        end if

        ! If we're stamping stdout with phase, times, etc.
        ! we'll want to show the last phase is over
        call getStamp( textcode=textcode, showTime=showTime )
        if ( showTime .or. textCode /= ' ' ) then
          call setStamp( textCode="Output" )
        end if
        ! And resume directwrites
        skipDirectwrites = skipDirectwritesOriginal

        ! ------------------------------------------- Output section
      case ( z_output ) ! Write out the data
        call resumeOutput ! In case the last phase was  silent
        if ( parallel%slave .and. .not. slavesCleanUpSelves ) then
          exit
        else if ( .not. parallel%slave ) then
          call Output_Close ( son, l2gpDatabase, l2auxDatabase, DirectDatabase, &
            & matrices, hessians, vectors, fileDataBase, griddedDataBase, &
            & chunks, processingRange, &
            & canWriteL2PC )
        end if

        call add_to_section_timing ( 'output', t1)

      case default
        ! This is one of the skipped sections
        if ( .not. parallel%master ) then
          call output('Skipping ', advance='no')
          call output(trim(section_name), advance='yes')
        end if
      end select
    end do

    ! Now finish up
    call finishUp

  contains

    ! Housekeeping
    ! Fortran compilers are good at garbage collection
    ! Nonetheless, we make a stab at it ourselves, hoping to
    ! avoid memory leaks
    subroutine FinishUp ( Early )
      ! Deallocate and utterly destroy all the arrays we allocated
      ! Some are known to us, ebing arrays of user-defined datatypes
      ! These we try to deallocate first
      use FilterShapes_M, only: Destroy_Filter_Shapes_Database, &
        & Destroy_Dacs_Filter_Database
      logical, intent(in), optional :: Early ! If so, can't deallocate everything
      ! Local variables
      logical :: myEarly
      integer :: numChunks
      ! Executable
      myEarly = .false.
      if ( present(early) ) myEarly = early
      numChunks = 0
      if ( associated(Chunks) ) numChunks = size(Chunks)
      ! This used to be done at the end of the Output section
      ! For case where there was one chunk, destroy vectors etc.
      ! This is to guard against destroying stuff needed by l2pc writing
      if ( canWriteL2PC ) then
        if ( warnOnDestroy ) call output('About to MLSL2Deconstruct', advance='yes' )
        call MLSL2DeConstruct ( vectorTemplates, hGrids )
        if ( warnOnDestroy ) call output('About to destroy hessian db', advance='yes' )
        call DestroyHessianDatabase ( hessians )
        if ( warnOnDestroy ) call output('About to destroy matrix db', advance='yes' )
        call DestroyMatrixDatabase ( matrices )
        if ( warnOnDestroy ) call output('About to destroy vector db', advance='yes' )
        call DestroyVectorDatabase ( vectors )
        if ( warnOnDestroy ) call output('Unable to destroy qty template db', advance='yes' )
        ! call destroyQuantityTemplateDatabase ( QtyTemplates )
      end if

      if ( specialDumpFile /= ' ' ) &
        & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
      ! Now tidy up any remaining `pointer' data.
      ! processingRange needs no deallocation
      details = SwitchDetail(switches, 'gridd')
      if ( details > -1 .and. .not. parallel%slave &
       & .and. associated(griddedDataBase) ) then
        call Dump( griddedDataBase, details )
      end if
      if ( warnOnDestroy ) call output('About to destroy gridded db', advance='yes' )
      call DestroyGriddedDataDatabase ( griddedDataBase )
      details = SwitchDetail(switches, 'l2gp')
      if ( details > -1 .and. .not. parallel%slave) then
        call Dump( l2gpDatabase, details=details )
      else if ( switchDetail(switches,'cab') > -1 .and. .not. parallel%slave) then
        call Dump( l2gpDatabase, ColumnsOnly=.true. )
      end if
      if ( warnOnDestroy ) call output('About to destroy l2gp db', advance='yes' )
      call DestroyL2GPDatabase ( l2gpDatabase )
      details = SwitchDetail(switches, 'l2aux')
      if ( details > -1 .and. .not. parallel%slave) then
        call Dump( l2auxDatabase, details=details )
      end if
      if ( warnOnDestroy ) call output('About to destroy l2aux db', advance='yes' )
      call DestroyL2AUXDatabase ( l2auxDatabase )
      if ( warnOnDestroy ) call output('About to destroy direct db', advance='yes' )
      call DestroyDirectDatabase ( DirectDatabase )
      ! vectors, vectorTemplates and qtyTemplates destroyed at the
      ! end of each chunk
      ! fileDataBase is deallocated in MLSL2
      details = SwitchDetail(switches, 'pro') - 2 ! 'pro' prints only size(DB)
      if ( details > -3 &
        & .and. associated(fileDataBase) ) then
        call Dump( fileDataBase, details=details, table=.true. )
      end if
      if ( specialDumpFile /= ' ' ) &
        & call revertOutput

      call SummarizeWarnings
      call CloseParallel(numChunks, early)
      ! Let's try to deallocate all the module-scope pointers
      ! trusting each module to know what they are
      if ( parallel%slave .and. .not. slavesCleanUpSelves ) return
      if ( .not. (myEarly .or. L2Options%SkipRetrieval .or. checkPaths) ) then
        call trace_begin ( 'Destroying databases', cond=toggle(gen) )
        if ( parallel%slave ) then
          ! Some things a slave can safely destroy
          call destroyChunkDatabase ( chunks )
          ! call destroy_Ant_Patterns_database
          call destroy_DACS_Filter_Database
          ! call destroy_Filter_Shapes_Database
          call destroyBinSelectorDatabase
          ! call destroyL2PCDatabase
          call destroyFWMConfigDatabase ( forwardModelConfigDatabase )
          call destroy_line_database
          ! call destroy_pointing_grid_database
          call destroy_spectcat_database
          call destroyBandDatabase ( Bands )
          call destroyModuleDatabase ( Modules )
          call destroyRadiometerDatabase ( Radiometers )
          call destroySpectrometerTypeDatabase ( SpectrometerTypes )
          call destroySignalDatabase ( Signals )
          call destroyVGridDatabase ( vGrids )
          call destroyFGridDatabase ( fGrids )
        else
          ! Others only a master can safely destroy
          call destroyChunkDatabase ( chunks )
          call DestroyHGridGeoLocations
          call destroy_Ant_Patterns_database
          call destroy_DACS_Filter_Database
          call destroy_Filter_Shapes_Database
          call destroyBinSelectorDatabase
          call destroyL2PCDatabase
          if ( switchDetail ( switches, 'l2pc' ) > -1 ) &
            & call output('Destroyed l2pc db', advance='yes')
          call destroyFWMConfigDatabase ( forwardModelConfigDatabase )
          call destroy_line_database
          call destroy_pointing_grid_database
          call destroy_spectcat_database
          call destroyBandDatabase ( Bands )
          call destroyModuleDatabase ( Modules )
          call destroyRadiometerDatabase ( Radiometers )
          call destroySpectrometerTypeDatabase ( SpectrometerTypes )
          call destroySignalDatabase ( Signals )
          call destroyVGridDatabase ( vGrids )
          call destroyFGridDatabase ( fGrids )
        end if
        call trace_end ( 'Destroying databases', cond=toggle(gen) )
      end if
      error_flag = 0
      call trace_end ( 'WALK_TREE_TO_DO_MLS_L2', cond=toggle(gen) )
    end subroutine FinishUp

  end subroutine WALK_TREE_TO_DO_MLS_L2

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module TREE_WALKER

! $Log$
! Revision 2.213  2019/09/19 16:13:01  pwagner
! Fixed bug in timings when processing multiple chunks serially
!
! Revision 2.212  2018/09/13 20:23:23  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.211  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.210  2018/04/19 23:44:36  pwagner
! Skip may take /nextChunk flag
!
! Revision 2.209  2018/04/19 00:49:44  vsnyder
! Remove USE statements and declarations for unused names
!
! Revision 2.208  2018/03/14 22:08:36  pwagner
! May print section in innemost loop if verbose enough
!
! Revision 2.207  2017/01/25 17:24:22  pwagner
! May skip certain Phases named in phasesToSkip cmdline opt
!
! Revision 2.206  2016/11/08 17:31:26  pwagner
! Use SayTime subroutine from time_m module
!
! Revision 2.205  2016/09/21 00:40:17  pwagner
! Usually dump FileDB as a table
!
! Revision 2.204  2016/08/09 22:08:33  pwagner
! May check for corrupt file database
!
! Revision 2.203  2016/05/18 01:37:30  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.202  2016/05/06 17:50:17  pwagner
! Still unable to reliably destroyQuantityTemplateDatabase
!
! Revision 2.201  2016/05/04 18:22:59  pwagner
! Restored ability of serial runs to process multiple chunks; updated api for MLSL2DeConstruct; deallocated qty templates before end
!
! Revision 2.200  2015/10/13 23:51:53  pwagner
! Will skip a chunk w/o profiles instead of crashing
!
! Revision 2.199  2015/09/03 20:29:47  pwagner
! We were calling L2MasterTask more than once; fixed
!
! Revision 2.198  2015/08/03 21:46:59  pwagner
! May debug Next_tree_node_m calls
!
! Revision 2.197  2015/06/19 20:40:04  pwagner
! Explicitly destroy HGrid gelocations if a master
!
! Revision 2.196  2014/10/07 00:22:32  pwagner
! -gn doesnt automatically dump chunks
!
! Revision 2.195  2014/09/05 01:28:14  vsnyder
! Add some tracing, remove unreferenced USSE name
!
! Revision 2.194  2014/06/30 23:30:45  pwagner
! Renamed global setting to slavesCleanUpSelves
!
! Revision 2.193  2014/04/10 00:43:11  pwagner
! Moved currentChunkNumber, currentPhaseName from MLSL2Timings to MLSL2Options
!
! Revision 2.192  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.191  2014/01/08 21:04:16  vsnyder
! Add more info to entry trace
!
! Revision 2.190  2013/12/13 21:25:30  vsnyder
! Use the iterator correctly so the correct number of chunks are run, and
! the output section is run
!
! Revision 2.189  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.188  2013/10/09 00:24:30  pwagner
! May have multiple Output sections
!
! Revision 2.187  2013/10/01 22:18:56  pwagner
! Passes hgrids to MLSL2Fill
!
! Revision 2.186  2013/08/30 02:45:52  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.185  2013/08/23 23:32:55  pwagner
! May use l2cf to set l1b file type to non-Aura
!
! Revision 2.184  2013/08/17 02:54:32  vsnyder
! Remove references to DEPTH from trace_m
!
! Revision 2.183  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.182  2013/04/05 23:25:41  pwagner
! Made 'master' a 'section' for timings summary
!
! Revision 2.181  2013/02/14 19:04:29  pwagner
! Consistent with changed L2MasterTask api
!
! Revision 2.180  2012/08/16 17:56:54  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.179  2012/04/26 23:15:28  pwagner
! Now tracks currentPhaseName and currentChunkNumber
!
! Revision 2.178  2011/10/07 00:06:02  pwagner
! May dump Matrices, Hessians from Fill, Join
!
! Revision 2.177  2011/06/29 21:51:17  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.176  2011/06/16 23:18:01  pwagner
! Pass details from switches when dumping dbs
!
! Revision 2.175  2011/05/09 18:27:56  pwagner
! Converted to using switchDetail
!
! Revision 2.174  2011/04/20 16:53:43  pwagner
! Added new flexibility to l2cf control flow by run-time booleans affecting gridded data
!
! Revision 2.173  2011/04/08 00:10:00  pwagner
! Raise error if every chunk was skipped
!
! Revision 2.172  2010/04/13 01:43:09  vsnyder
! Move FlushLockedBins from LinearizedForwardModel_m to L2PCBins_m
!
! Revision 2.171  2010/03/17 20:57:35  pwagner
! Print about destroyed l2pc db only if -S'l2pc'
!
! Revision 2.170  2010/02/25 18:20:44  pwagner
! Adds support for new Hessian database
!
! Revision 2.169  2010/02/04 19:10:58  pwagner
! Allocate directwrite db with zero size
!
! Revision 2.168  2009/10/01 19:58:00  vsnyder
! Pass file database to set_global_settings
!
! Revision 2.167  2009/06/23 18:46:19  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.166  2008/12/18 21:12:41  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.165  2008/09/19 23:56:28  pwagner
! May now Destroy GriddedData
!
! Revision 2.164  2008/05/20 00:29:19  vsnyder
! Don't leave J undefined if no chunks processed
!
! Revision 2.163  2008/04/01 17:00:15  pwagner
! Can get hdf5 spectroscopy file path/name from PCF
!
! Revision 2.162  2007/12/07 01:15:01  pwagner
! Removed unused dummp varaibles, etc.
!
! Revision 2.161  2007/10/09 00:34:21  pwagner
! ifort was crashing in destroy_DACS_Filter_Database during checkPaths; we now skip
!
! Revision 2.160  2007/09/13 21:51:39  pwagner
! Resumed destroying dbs (unless slave)
!
! Revision 2.159  2007/09/06 23:35:16  pwagner
! slaves need not do own cleanup
!
! Revision 2.158  2007/08/31 00:03:26  pwagner
! Summarizes warnings at end
!
! Revision 2.157  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.156  2007/06/07 20:41:07  pwagner
! Avoid read_ and Merge_apriori if master task
!
! Revision 2.155  2007/04/03 17:36:44  vsnyder
! Allocate Vectors with zero size so others don't need to check whether
! it's associated; indeed, they don't need the pointer attribute.
!
! Revision 2.154  2007/03/23 00:29:41  pwagner
! Switch destroy warns when destroying dbs
!
! Revision 2.153  2007/02/14 17:31:00  pwagner
! Commented-out destroyL2PCDatabase (was it killing slaves?)
!
! Revision 2.152  2007/01/12 00:34:04  pwagner
! Renamed routine outputNamedValue
!
! Revision 2.151  2006/10/09 18:39:35  pwagner
! Fixed bug preventing a call to Output_CLose
!
! Revision 2.150  2006/10/05 23:32:43  pwagner
! skipSections can skip named sections
!
! Revision 2.149  2006/08/14 16:22:28  pwagner
! Should not double-deallocate if canWriteL2PC
!
! Revision 2.148  2006/08/10 21:46:50  pwagner
! --chunk commandline option now synonym for --chunkRange
!
! Revision 2.147  2006/08/05 02:10:46  vsnyder
! Delete some debugging print I shouldn't have left at the last check-in
!
! Revision 2.146  2006/08/05 00:41:50  vsnyder
! Comment out filter database destruction -- causes memory trouble
!
! Revision 2.145  2006/08/02 19:53:26  vsnyder
! Send Vectors database to outputClose for destroy command
!
! Revision 2.144  2006/07/24 20:35:04  pwagner
! Fixed bug when stopAfterSection is blank
!
! Revision 2.143  2006/07/21 20:12:26  pwagner
! Can select what section to stop after
!
! Revision 2.142  2006/06/12 16:28:25  pwagner
! Added ability to dump Gridded Data
!
! Revision 2.141  2006/05/05 16:49:23  pwagner
! May convertEtaToP and create a VGrid in MergeGrids section
!
! Revision 2.140  2006/04/11 23:30:00  pwagner
! chunkRange option effective in serial runs, too
!
! Revision 2.139  2006/03/04 00:22:04  pwagner
! Restore original skipDirectWrites after chunkLoop ends
!
! Revision 2.138  2006/02/16 00:16:32  pwagner
! switchDetail instead of index
!
! Revision 2.137  2006/02/10 21:15:26  pwagner
! dumps may go to special dumpfile
!
! Revision 2.136  2006/01/11 17:02:16  pwagner
! Repaired erroneous report when single chunk > size(chunks)
!
! Revision 2.135  2006/01/06 01:15:01  pwagner
! Resumes output at start of each chunk, and OutputClose
!
! Revision 2.134  2005/08/19 23:35:01  pwagner
! Allow Output to repair l2gp with HGrid while copying files
!
! Revision 2.133  2005/07/12 17:37:20  pwagner
! Removed unused DestroyL1BInfo
!
! Revision 2.132  2005/06/14 20:44:36  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.131  2005/06/03 02:05:29  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of passing as an argument.
!
! Revision 2.130  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.129  2004/12/14 21:56:50  pwagner
! Changes related to stopping early
!
! Revision 2.128  2004/07/22 20:49:58  cvuu
! Add forwardModelConfigDatabase to the call MLSL2Join and MLSL2Fill
!
! Revision 2.127  2004/05/19 19:16:12  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.126  2004/04/28 23:07:44  livesey
! Now passes more stuff to Algebra
!
! Revision 2.125  2004/02/10 19:29:36  pwagner
! Prints time for processing each chunk at chunks end
!
! Revision 2.124  2004/01/17 00:28:09  vsnyder
! Provide for Algebra section
!
! Revision 2.123  2004/01/14 18:49:58  vsnyder
! Stuff to support the Algebra section
!
! Revision 2.122  2003/12/16 01:28:56  livesey
! Moved the destruction of the chunk database to after closeParallel.
!
! Revision 2.121  2003/12/11 22:59:32  pwagner
! May fill DirectWriteDatabase in global settings
!
! Revision 2.120  2003/12/05 00:42:05  pwagner
! Removed last vestige of explicit garbage collection
!
! Revision 2.119  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.118  2003/11/05 21:27:54  pwagner
! Can enter range of chunks to be processed instead of single
!
! Revision 2.117  2003/10/09 23:31:26  pwagner
! Changed grid switch to gridd
!
! Revision 2.116  2003/09/03 15:56:18  cvuu
! Move the do loop over the forwardModel inside subroutine PrintForwardModelTiming
!
! Revision 2.115  2003/09/02 18:04:38  pwagner
! Can do a singleChunk even if master task
!
! Revision 2.114  2003/08/21 21:24:39  cvuu
! Change the output format for fullForwardModel Timing
!
! Revision 2.113  2003/08/21 16:07:34  livesey
! Now calls FlushLockedBins automatically at the end of each chunk.
!
! Revision 2.112  2003/08/11 23:24:02  pwagner
! Stores ChunkNo as component of L2ParallelInfo_T
!
! Revision 2.111  2003/07/07 23:51:06  pwagner
! Need not pass around l2pc as L2pcf now a saved variable in WriteMetaData
!
! Revision 2.110  2003/06/30 22:56:16  cvuu
! Print mean, std dev for fullForwardModel timing
!
! Revision 2.109  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.108  2003/06/24 23:00:53  livesey
! Nullified directDatabase
!
! Revision 2.107  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.106  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.105  2003/06/09 22:51:35  pwagner
! Renamed scan_divide to chunk_divide in timings table
!
! Revision 2.104  2003/03/08 00:46:08  pwagner
! Per njl, checks for illegal single chunk runs
!
! Revision 2.103  2003/02/21 21:03:39  pwagner
! Disabled tps switch to eliminate need for Test_Parse_Signals
!
! Revision 2.102  2003/01/13 20:59:14  livesey
! Removed a print statement
!
! Revision 2.101  2002/12/05 19:45:20  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.100  2002/12/04 01:18:21  pwagner
! First halting steps toward using filedatabase
!
! Revision 2.99  2002/11/22 12:24:25  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.98  2002/11/22 01:14:06  vsnyder
! Remove USE'd but unreferenced symbols and two unused local variables
!
! Revision 2.97  2002/11/21 18:45:04  livesey
! Changes to way the destroy stuff is called based on chunks
!
! Revision 2.96  2002/11/21 18:38:34  livesey
! Bug fix in canWriteL2PC flag passed to output_close
!
! Revision 2.95  2002/10/08 17:41:50  livesey
! Various bug fixes associated with FWMParallel
!
! Revision 2.94  2002/10/08 17:36:23  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.93  2002/10/05 00:44:29  livesey
! Included the FWMParallel stuff
!
! Revision 2.92  2002/09/25 20:09:32  livesey
! Changes to allow chunk based forwardModel configs
!
! Revision 2.91  2002/09/24 22:17:50  pwagner
! Defines t1 before call to add_to_section_timing
!
! Revision 2.90  2002/09/24 18:18:47  pwagner
! Consistent with add_to_section_timing now calling time_now at its end
!
! Revision 2.89  2002/08/28 22:29:19  pwagner
! Moved DestroyL1BInfo to after end of loop of chunks
!
! Revision 2.88  2002/08/22 01:26:00  vsnyder
! Cosmetic changes
!
! Revision 2.87  2002/08/21 19:05:04  livesey
! Removed calls to H5Open/close as they are now done in MLSL2.
!
! Revision 2.86  2002/08/21 00:55:35  livesey
! Changed error to warning when fail to close hdf5 (at least for the
! moment).
!
! Revision 2.85  2002/08/20 22:10:50  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.84  2002/08/07 00:05:27  livesey
! Added calls to H5Open_F and H5Close_F
!
! Revision 2.83  2002/08/04 16:10:43  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.82  2002/03/20 00:49:44  pwagner
! chi1 dumps only 1st chunk, unlike chi
!
! Revision 2.81  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.80  2002/01/24 00:58:28  livesey
! Now calls MergeGrids at the appropriate time
!
! Revision 2.79  2002/01/22 18:14:47  livesey
! Fixed typo
!
! Revision 2.78  2002/01/21 23:11:06  livesey
! Added call to DestroyBinSelectorsDatabase etc.
!
! Revision 2.77  2002/01/18 18:55:37  livesey
! Code to support the --chunk option
!
! Revision 2.76  2002/01/09 22:56:17  livesey
! Now sends slaves all chunks as regular HGrids need them.
!
! Revision 2.75  2001/12/16 00:57:12  livesey
! Now passes all chunks to construct
!
! Revision 2.74  2001/12/14 01:43:20  livesey
! Passes processingRange to HGrid via Construct
!
! Revision 2.73  2001/12/13 23:21:26  livesey
! Added countChunks option
!
! Revision 2.72  2001/12/10 20:22:09  livesey
! Added code for EmpiricalGeometry.
!
! Revision 2.71  2001/11/20 00:48:15  livesey
! Fixed problem when no chunks to process
!
! Revision 2.70  2001/11/16 17:24:13  livesey
! Now calls new ChunkDivide routine.
!
! Revision 2.69  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.68  2001/10/31 19:07:52  livesey
! Hooked fGrids into quantity templates
!
! Revision 2.67  2001/10/31 18:36:58  livesey
! Added fGrids stuff
!
! Revision 2.66  2001/10/26 23:16:15  pwagner
! Similar dump interfaces for l2gp, l2aux, Griddeddata databases
!
! Revision 2.65  2001/10/12 23:14:22  pwagner
! Debugging when dumping diagnostics; may remove later
!
! Revision 2.64  2001/10/09 23:43:42  pwagner
! Some further improvements in dumping vectors
!
! Revision 2.63  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models
!
! Revision 2.62  2001/09/29 00:01:00  pwagner
! Fixed various timing problems
!
! Revision 2.61  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.60  2001/09/21 17:40:57  pwagner
! May dump diagnostic quantities chi..
!
! Revision 2.59  2001/09/10 23:37:32  livesey
! New GriddedData etc.
!
! Revision 2.58  2001/08/06 18:34:59  pwagner
! Now dumps l2gp and l2aux databases when asked
!
! Revision 2.57  2001/07/11 21:40:21  livesey
! Added -Schu option
!
! Revision 2.56  2001/06/13 20:44:08  livesey
! Moved the CloseParallel higher up, to work around the memory management
! problem (somethine [HDF?] seems to stamp on PointingFrequencyDatabase)
!
! Revision 2.55  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.54  2001/05/23 22:00:16  livesey
! Interim version
!
! Revision 2.53  2001/05/23 01:44:35  livesey
! Parallel stuff taking shape
!
! Revision 2.52  2001/05/10 00:43:23  livesey
! Tree walker now owns hGrids
!
! Revision 2.51  2001/05/04 17:12:25  pwagner
! Passes necessary args to global_settings
!
! Revision 2.50  2001/05/03 20:34:08  vsnyder
! Cosmetic changes
!
! Revision 2.49  2001/05/02 23:23:00  livesey
! Added some of the parallel stuff
!
! Revision 2.48  2001/04/28 01:31:46  livesey
! Changes for new l2pc / matrix handling.
!
! Revision 2.47  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.46  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.45  2001/04/26 00:08:02  livesey
! Stuff to support reading of l2pc files
!
! Revision 2.44  2001/04/25 21:54:16  livesey
! Added canDoL2PC flag to join
!
! Revision 2.43  2001/04/25 21:52:13  livesey
! Moved DeConstruct to after output for 1 chunk cases.
! This is to protect vectors and matrices stored in l2pcs.
!
! Revision 2.42  2001/04/25 20:34:36  livesey
! Now supports writing of l2pc files
!
! Revision 2.41  2001/04/25 19:31:13  livesey
! Fixed bug, now nullifies l2pcDatabase
!
! Revision 2.40  2001/04/24 23:05:54  vsnyder
! Make 'test_parse_signals' depend on 'switches' containing 'tps'
!
! Revision 2.39  2001/04/24 20:20:02  livesey
! L2PC moved to lib, and renamed
!
! Revision 2.38  2001/04/24 20:05:50  livesey
! New stuff to support joining of l2pc's
!
! Revision 2.37  2001/04/21 01:41:35  vsnyder
! Fix memory leaks
!
! Revision 2.36  2001/04/21 01:26:37  livesey
! Now passes l2gpDatabase to more people
!
! Revision 2.35  2001/04/20 17:12:38  livesey
! Add vGrids argument to fill to support fill from vGrid
!
! Revision 2.34  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.33  2001/04/11 17:47:47  pwagner
! presets anText to null
!
! Revision 2.32  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.31  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.30  2001/04/10 00:02:19  vsnyder
! Implement 'matrix' spec in Fill section
!
! Revision 2.29  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.28  2001/04/06 20:12:59  vsnyder
! Make 'call test_parse_signal' depend on 'emit' toggle
!
! Revision 2.27  2001/04/04 02:15:12  vsnyder
! Add Spectroscopy section
!
! Revision 2.26  2001/04/03 20:50:45  pwagner
! Added anText to hold PCF file contents
!
! Revision 2.25  2001/04/02 23:41:09  pwagner
! Now keeps l2pcf and transmits as needed
!
! Revision 2.24  2001/03/29 19:13:03  livesey
! Renamed apriorDatabase to griddedData
!
! Revision 2.23  2001/03/28 23:47:48  livesey
! Added arguments to sids etc.
!
! Revision 2.22  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.21  2001/03/17 03:30:25  vsnyder
! Remove FMI and TFMI from the call to set_global_settings, since it no
! longer uses them.
!
! Revision 2.20  2001/03/17 00:45:53  livesey
! Added forwardModelConfigDatabase
!
! Revision 2.19  2001/03/15 23:26:09  livesey
! Added chunks to call to MLSL2Join
!
! Revision 2.18  2001/03/15 21:20:54  pwagner
! Split between GriddedData and ncep_dao modules
!
! Revision 2.17  2001/03/15 21:05:23  vsnyder
! Set up to test Parse_Signal
!
! Revision 2.16  2001/03/14 02:04:53  vsnyder
! Moved MLSSignals_m to mlspgs/lib
!
! Revision 2.15  2001/03/09 19:57:48  vsnyder
! Add 'z_retrieve' to a case branch
!
! Revision 2.14  2001/03/08 18:21:11  vsnyder
! Even more stuff for L2_Load
!
! Revision 2.13  2001/03/08 03:23:10  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.12  2001/03/07 22:46:05  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.11  2001/03/03 00:13:30  pwagner
! Gets read_apriori from ReadAPriori module
!
! Revision 2.10  2001/02/28 01:17:57  livesey
! Removed obtain_ncep etc. These will later be in lower down modules
!
! Revision 2.9  2001/02/27 17:38:07  livesey
! Tidied up arguments to MLSL2Join
!
! Revision 2.8  2001/02/23 02:51:44  vsnyder
! Improve progress messages triggered by -g option
!
! Revision 2.7  2001/02/08 01:40:53  vsnyder
! Don't know what I have done, but "cvs update" said M instead of U
!
! Revision 2.6  2001/01/03 17:48:43  pwagner
! added chunk args to call to Fill
!
! Revision 2.5  2000/12/05 00:41:50  pwagner
! Added L2AUXDatabase arg in call to MLSL2Fill
!
! Revision 2.4  2000/11/30 00:21:08  pwagner
! passes l2*databses to read a priori
!
! Revision 2.3  2000/11/29 00:27:54  pwagner
! Began changes to open old l2gp
!
! Revision 2.2  2000/10/10 00:37:46  vsnyder
! Added $Log for CVS at end.
!
