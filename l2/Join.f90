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
module Join                     ! Join together chunk based data.
!=============================================================================

  use MLSStringLists, only: SWITCHDETAIL
  ! This module performs the 'join' task in the MLS level 2 software.

  implicit none
  private
  public :: MLSL2Join, JoinL2GPQuantities, JoinL2AuxQuantities

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! logical, parameter, private :: FORCEDIRWRITEREOPEN = .false. ! Usually FALSE
  logical, parameter, private :: SKIPMETADATA = .false. ! Usually FALSE
  real, parameter             :: timeReasonable = 500.

  ! Parameters for Announce_Error

  integer :: ERROR
  logical, parameter :: CATENATESPLITS  = .true.
  integer, parameter :: NO_ERROR_CODE   = 0
  integer, parameter :: NotAllowed      = 1

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  MLSL2Join  -----

  ! This is the main routine for the Join block.  This one just goes
  ! through the tree and dispatches work to other routines.

  subroutine MLSL2Join ( root, vectors, l2gpDatabase, l2auxDatabase, &
    & DirectDataBase, chunkNo, chunks, FWModelConfig, filedatabase, HGrids, &
    & Matrices, Hessians )
    ! Imports
    use ALLOCATE_DEALLOCATE, only: TEST_ALLOCATE
    use CHUNKS_M, only: MLSCHUNK_T
    use DIRECTWRITE_M, only: DIRECTDATA_T
    use DUMPCOMMAND_M, only: DUMPCOMMAND, &
      & MLSCASE, MLSENDSELECT, MLSSELECT, MLSSELECTING, SKIP
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use HESSIANMODULE_1, only: HESSIAN_T
    use HGRIDSDATABASE, only: HGRID_T
    use INIT_TABLES_MODULE, only: S_L2GP, S_L2AUX, S_TIME, S_DIRECTWRITE, &
      & S_ENDSELECT, S_CASE, S_DIFF, S_DUMP, S_LABEL, S_SELECT, S_SKIP
    use L2GPDATA, only: L2GPDATA_T
    use L2AUXDATA, only: L2AUXDATA_T
    use L2PARINFO, only: PARALLEL, WAITFORDIRECTWRITEPERMISSION
    use LEXER_CORE, only: PRINT_SOURCE
    use MatrixModule_1, only: MATRIX_DATABASE_T
    use MLSCOMMON, only: MLSFILE_T
    use MLSL2OPTIONS, only: CHECKPATHS, L2CFNODE, &
      & SKIPDIRECTWRITES, SPECIALDUMPFILE, MLSMESSAGE
    use MLSL2TIMINGS, only: SECTION_TIMES, TOTAL_TIMES, &
      & ADD_TO_DIRECTWRITE_TIMING, ADD_TO_SECTION_TIMING
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR
    use MORETREE, only: Get_Label_And_Spec, GET_SPEC_ID
    use Next_Tree_Node_m, only: Init_Next_Tree_Node, Next_Tree_Node, &
      & Next_Tree_Node_State
    use OUTPUT_M, only: OUTPUT, REVERTOUTPUT, SWITCHOUTPUT
    use TOGGLES, only: GEN, TOGGLE, SWITCHES
    use TREE, only: Where_At => Where
    use TIME_M, only: TIME_NOW
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VECTORSMODULE, only: VECTOR_T
    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the JOIN section in the AST
    type (Vector_T), dimension(:), pointer :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (HGrid_T), dimension(:), pointer ::     HGrids
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (Hessian_T), dimension(:), pointer :: Hessians

    ! Local parameters
    ! integer, parameter :: DELAY = 500000  ! For Usleep, no. microsecs
    ! External (C) function
    external :: Usleep
    ! Local variables
    logical :: AUTODIRECTWRITE
    logical :: CREATEFILE               ! Flag
    integer :: DBINDEX
    logical :: DIDTHEWRITE
    integer :: DIRECTWRITENODEGRANTED   ! Which request was granted
    real :: DWT1                        ! Time we started
    real :: DWT2                        ! Time we finished
    real :: DWT22                       ! Time we finished, too
    integer :: KEY                      ! Tree node
    integer :: Label                    ! Not actually used
    integer :: Me = -1                  ! String index for trace
    integer :: NODIRECTWRITES           ! Array size
    integer :: NODIRECTWRITESCOMPLETED  ! Counter
    integer :: NOEXTRAWRITES            ! Correction to NODIRECTWRITES
    integer :: PASS                     ! Loop counter
    integer :: SON                      ! Tree node
    integer :: SPECID                   ! Type of l2cf line this is
    type(next_tree_node_state) :: State ! of tree traverser
    integer :: STATUS
    integer :: TICKET                   ! Direct write permission ticket
    ! integer :: THEFILE                  ! Direct write permission on file
    character(len=4096) :: THEFILE  ! Which file permission granted for
    logical :: TIMING                   ! Flag
    real :: T1                          ! Time we started
    real :: T2                          ! Time we finished
    logical :: namedFile                ! set true if DirectWrite named file
    logical :: DEEBUG
    
    ! Executable code
    DEEBUG = ( switchDetail(switches, 'direct') > -1 )!  .or. .true.
    call trace_begin ( me, "MLSL2Join", root, cond=toggle(gen) )
    timing = section_times
    if ( timing ) call time_now ( t1 )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ! Will we (perhaps) be automatically assigning Direct writes to
    ! files declared in global settings section?
    autoDirectWrite = .false.
    if ( associated(DirectDataBase) ) then
      if ( size(DirectDataBase) > 0 ) then
        autoDirectWrite = any(DirectDataBase%autoType > 0)
      end if
    end if
    ! This is going to be somewhat atypical, as the code may run in 'passes'
    ! In pass 1 we do all the regular join statements and count the direct writes
    ! In pass 2 we log all our direct write requests.
    ! In the later passes (as many as there are direct writes) we do the
    ! direct writes we've been given permission for.

    ! In the non-parallel slave mode, one pass is sufficient.

    error = 0
    pass = 1
    noDirectWrites = 0
    noDirectWritesCompleted = 0
    ticket = 0                          ! Default value for serial case
    didthewrite = .true.
    passLoop: do
      ! For the later passes, we wait for permission to do one of the direct writes
      if ( pass > 2 ) then
        if ( skipDirectWrites ) exit passLoop
        call time_now ( dwt2 )
        call WaitForDirectWritePermission ( directWriteNodeGranted, ticket, &
          & theFile, createFile )
        call time_now ( dwt22 )
        if ( dwt22-dwt2 > timeReasonable ) then
          call output('Unreasonable time waiting for permission', advance='yes')
        endif
        if ( parallel%verbosity > 0 ) then
          call output ( "Got permission for ticket " )
          call output ( ticket )
          call output ( " node " )
          call output ( directWriteNodeGranted, advance='no' )
          call output ( " file " )
          call output ( trim(theFile), advance='yes' )
        endif
        call add_to_directwrite_timing ( 'waiting', dwt2)
        didthewrite = .false.
      end if
      
      ! Simply loop over lines in the l2cf
      call init_next_tree_node ( state )
      do 
        son = next_tree_node ( root, state )
        if ( son == 0 ) exit
        call get_label_and_spec ( son, label, key )
        L2CFNODE = key
        if ( MLSSelecting .and. &
          & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) cycle
        specId = get_spec_id ( key )
        select case ( specId )
        case ( s_diff, s_dump ) ! ======================= Diff, Dump ==========
          ! Handle disassociated pointers by allocating them with zero size
          status = 0
          if ( CHECKPATHS ) cycle
          if ( .not. associated(vectors) ) allocate ( vectors(0), stat=status )
          call test_allocate ( status, moduleName, 'Vectors', (/0/), (/0/) )
          call dumpCommand ( key, vectors=vectors, HGrids=HGrids, &
            & ForwardModelConfigs=FWModelConfig, FileDataBase=FileDataBase, &
            & MatrixDatabase=matrices, HessianDatabase=Hessians )
        case ( s_skip ) ! ============================== Skip ==========
          ! We'll skip the rest of the section if the Boolean cond'n is TRUE
          if ( Skip(key) ) exit
        case ( s_select ) ! ============ Start of select .. case ==========
          ! We'll start seeking a matching case
          call MLSSelect (key)
        case ( s_case ) ! ============ seeking matching case ==========
          ! We'll continue seeking a match unless the case is TRUE
          call MLSCase (key)
        case ( s_endSelect ) ! ============ End of select .. case ==========
          ! We'done with seeking a match
          call MLSEndSelect (key)
        case ( s_time )
          ! Only say the time the first time round
          if ( pass == 1 ) then
            if ( timing ) then
              call sayTime
            else
              call time_now ( t1 )
              timing = .true.
            end if
          end if
        case ( s_l2gp, s_l2aux )
          ! Only do these the first time round
          if ( pass == 1 .and. .not. checkpaths ) then
            call JoinQuantities ( key, vectors, l2gpDatabase, l2auxDatabase, &
              & chunkNo, chunks )
          end if
        case ( s_label )
          ! Only do these the first time round
          if ( pass == 1 .and. .not. checkpaths ) then
            call LabelVectorQuantity ( son, vectors )
          end if
        case ( s_directWrite )
          ! if ( SKIPDIRECTWRITES ) exit
          call time_now ( dwt1 )
          if ( pass == 1 ) then
            ! On the first pass just count the number of direct writes
            noDirectWrites = noDirectWrites + 1
            ! Unless we're not a slave in which case just get on with it.
            if ( .not. parallel%slave ) then
              if(DEEBUG)print*,'Calling DirectWrite, not slave'
              call time_now ( dwt2 )
              if ( autoDirectWrite .and. .not. explicitFile(son) ) then
                ! Pretend we're given permission to write to each of the files
                do dbIndex=1, size(DirectdataBase)
                  if ( DirectDataBase(dbIndex)%autoType < 1 ) cycle
                  if(DEEBUG)print*,'DirectWrite to ', &
                    & trim(DirectDataBase(dbIndex)%fileNameBase)
                  call DirectWriteCommand ( son, ticket, vectors, &
                    & DirectdataBase, filedatabase, &
                    & chunkNo, chunks, FWModelConfig, HGrids, &
                    & theFile=DirectDataBase(dbIndex)%fileNameBase, &
                    & namedFile=namedFile )
                  if ( namedFile ) exit
                end do
              else
                call DirectWriteCommand ( son, ticket, vectors, &
                  & DirectdataBase, filedatabase, &
                  & chunkNo, chunks, FWModelConfig, HGrids )
              end if
              call add_to_directwrite_timing ( 'writing', dwt2)
            end if
          else if ( pass == 2 ) then
            ! On the second pass, log all our direct write requests.
            ! (and correct the count of Direct Writes if we're automatically
            !  distributing them, instead of relying on one line/one write)
            if(DEEBUG)call print_source ( where_at(son) )
            if(DEEBUG)print*,'Calling direct write to do a setup'
            call DirectWriteCommand ( son, ticket, vectors, &
              & DirectdataBase, fileDatabase,  &
              & chunkNo, chunks, FWModelConfig, HGrids, makeRequest=.true., &
              & NoExtraWrites=noExtraWrites )
            noDirectWrites = noDirectWrites + noExtraWrites
          else
            ! On the later passes we do the actual direct write we've been
            ! given permission for.
            if ( son == directWriteNodeGranted ) then
              didTheWrite = .true.
              call time_now ( dwt2 )
              if(DEEBUG)call print_source ( where_at(son) )
              if(DEEBUG)print*,'Calling direct write to do the write'
              if(DEEBUG)print*,'Asked to create file? ', createFile
              if(DEEBUG)print*,'the file ', trim(theFile)
              call DirectWriteCommand ( son, ticket, vectors, &
                & DirectdataBase, filedatabase, &
                & chunkNo, chunks, FWModelConfig, HGrids, create=createFile, &
                & theFile=theFile )
              call time_now(dwt22)
              if ( dwt22-dwt2 > timeReasonable .and. &
                & switchDetail(switches,'dwreq') > -1 ) then
                call output('Unreasonable time for directwritecommand', advance='yes')
                call output('Chunk: ', advance='no')
                call output(ChunkNo, advance='yes')
                call output('File: ', advance='no')
                call output(trim(theFile), advance='yes')
              endif
              call add_to_directwrite_timing ( 'writing', dwt2)
              noDirectWritesCompleted = noDirectWritesCompleted + 1
              ! If that was the last one then bail out
              if(DEEBUG)print*,'noDirectWritesCompleted: ', noDirectWritesCompleted, &
                & 'noDirectWrites: ', noDirectWrites
              if ( noDirectWritesCompleted == noDirectWrites ) exit passLoop
            end if
          end if
          call add_to_section_timing ( 'directwrite', dwt1)
        end select
      end do                            ! End loop over l2cf lines

      ! If we're not in parallel mode then one pass is enough
      if ( .not. parallel%slave ) exit passLoop

      ! Bail out of pass loop if there are no direct writes, or there was
      ! an error.
      ! if ( noDirectWrites == 0 .or. SKIPDIRECTWRITES .or. error /= 0 ) exit passLoop
      if ( noDirectWrites == 0 .or. error /= 0 ) exit passLoop
      pass = pass + 1
      ! Did we receive permission to write to a "ghost node"
      if ( .not. didTheWrite ) then
        call announce_error(root, no_error_code, &
            & 'Ghost node' ,ticket, 0 )
        print *, 'Ticket: ', ticket
        print *, 'ghost node: ', directWriteNodeGranted
        print*,'noDirectWritesCompleted: ', noDirectWritesCompleted, &
                & 'noDirectWrites: ', noDirectWrites        
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'We received permission to write, but could not find node' )
      end if
    end do passLoop                     ! End loop over passes

    if ( specialDumpFile /= ' ' ) &
      & call revertOutput
    ! Check for errors
    if ( error /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Error in Join section' )

    call trace_end ( "MLSL2Join", cond=toggle(gen) )
    if ( timing ) call sayTime

  contains
    ! Private procedure
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), after='    ', advance = 'no' )
      end if
      call output ( "Timing for MLSL2Join =" )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine MLSL2Join

  ! -----------------------------------------  DirectWriteCommand  -----
  ! 
  ! Direct write is probably going to become the predominant form
  ! of output in the software, as the other forms have become a
  ! little too intensive.

  ! This routine may be called by a standalone run, or by a slave.
  ! If called by a standalone run, just do the write and return.

  ! For a slave it will be called many times for each direct write.  The first
  ! is a pass through to check the syntax etc.
  ! If it is successful, in a second pass (makeRequests=.true.)
  ! a request will be made to the master.

  ! Subsequent calls will arrive with permission to actually do the writing,
  ! each time to a specific named file.

  ! There are two modes:
  ! (1) The DirectWrite command named the file explicitly--then write the
  !     quantities to that file
  ! (2) The DirectWrite command did not name any file
  !     then distribute the quantities among the appropriate DirectFiles
  !     i.e., swaths go to DGG files, hdf datasets to DGM files
  subroutine DirectWriteCommand ( node, ticket, vectors, &
    & DirectDataBase, fileDatabase, &
    & chunkNo, chunks, FWModelConfig, HGrids, makeRequest, create, theFile, &
    & noExtraWrites, namedFile )
    ! Imports
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use CHUNKS_M, only: MLSCHUNK_T
    use DIRECTWRITE_M, only: DIRECTDATA_T, &
      & DIRECTWRITE, DUMP, &
      & EXPANDDIRECTDB, EXPANDSDNAMES, FILENAMETOID
    use DUMP_0, only: DUMP
    use EXPR_M, only: EXPR
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use HDF, only: DFACC_CREATE, DFACC_RDWR
    use HGRIDSDATABASE, only: HGRID_T
    use HIGHOUTPUT, only: OUTPUTNAMEDVALUE
    use INIT_TABLES_MODULE, only: F_CONVERGENCE, F_FILE, F_HDFVERSION, &
      & F_LOWEROVERLAP, F_OPTIONS, F_AscDescMode, F_PRECISION, F_QUALITY, &
      & F_SINGLE, F_SOURCE, F_STATUS, F_TYPE, F_UPPEROVERLAP, F_VECTOR, &
      & FIELD_FIRST, FIELD_LAST
    use INIT_TABLES_MODULE, only: L_L2GP, L_L2AUX, L_L2DGG, L_L2FWM, &
      & L_PRESSURE, L_QUANTITY, L_ZETA
    use INTRINSIC, only: L_NONE, L_HDF, L_SWATH, LIT_INDICES, PHYQ_DIMENSIONLESS
    use L2PARINFO, only: PARALLEL, LOGDIRECTWRITEREQUEST, FINISHEDDIRECTWRITE
    use MANIPULATEVECTORQUANTITIES, only: DOHGRIDSMATCH
    use MLSCOMMON, only: FILENAMELEN, MLSFILE_T, L2METADATA_T
    use MLSFILES, only: HDFVERSION_5, &
      & ADDINITIALIZEMLSFILE, DUMP, GETMLSFILEBYNAME, GETPCFROMREF, &
      & MLS_CLOSEFILE, MLS_OPENFILE, SPLIT_PATH_NAME
    use MLSHDFEOS, only: MLS_SWATH_IN_FILE
    use MLSKINDS, only: R8
    use MLSL2OPTIONS, only: CHECKPATHS, &
      & DEFAULT_HDFVERSION_WRITE, PATCH, SKIPDIRECTWRITES, TOOLKIT
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSPCF2, only: MLSPCF_L2GP_START, MLSPCF_L2GP_END, &
      & MLSPCF_L2DGM_START, MLSPCF_L2DGM_END, MLSPCF_L2FWM_FULL_START, &
      & MLSPCF_L2FWM_FULL_END, &
      & MLSPCF_L2DGG_START, MLSPCF_L2DGG_END
    use MLSFINDS, only: FINDFIRST, FINDNEXT
    use MORETREE, only: GET_FIELD_ID, GET_BOOLEAN
    use OUTPUT_M, only: BLANKS, OUTPUT
    use OUTPUTANDCLOSE, only: ADD_METADATA
    use STRING_TABLE, only: DISPLAY_STRING, GET_STRING
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: GEN, TOGGLE, SWITCHES
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T, VALIDATEVECTORQUANTITY, &
      & GETVECTORQTYBYTEMPLATEINDEX
    ! Dummy arguments
    integer, intent(in) :: NODE         ! Of the JOIN section in the AST
    integer, intent(in) :: TICKET       ! Ticket number from master
    type (Vector_T), dimension(:), pointer :: VECTORS
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (HGrid_T), dimension(:), pointer ::     HGrids
    integer, intent(in) :: CHUNKNO
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    ! The next 3 args are used in the multi-pass followed by each slave:
    ! 1st pass just counts up requests to write
    ! 2nd pass just logs requests to write
    ! later passes come with permission granted to write
    logical, intent(in), optional :: MAKEREQUEST  ! Set on first pass through
    logical, intent(in), optional :: CREATE       ! Set if slave is to create file.
    character(len=*), intent(in), optional :: THEFILE  ! Which file permission granted for
    integer, intent(out), optional :: NOEXTRAWRITES ! How many extras if distributing
    logical, intent(out), optional :: namedFile ! Did it name the file?
    ! Local parameters
    integer, parameter :: MAXFILES = 1000 ! Set for an internal array
    ! Saved variable - used to work out the createFile flag in the serial case
    integer, dimension(maxFiles), save :: CREATEDFILENAMES = 0
    integer, save :: NOCREATEDFILES=0   ! Number of files created
    ! Local variables
    integer :: AFILE
    logical :: CREATEFILEFLAG           ! Flag (often copy of create)
    logical, dimension(:), pointer :: CREATETHISSOURCE
    logical :: CREATETHISSWATH
    logical :: DEEBUG
    type(MLSFile_T), pointer :: directFile
    logical :: DISTRIBUTINGSOURCES      ! No field 'file=...' in l2cf line
    logical :: DUMMY
    integer :: ERRORTYPE
    integer :: EXPECTEDTYPE             ! l2gp/l2aux
    integer :: EXPRUNITS(2)             ! From expr
    real (r8) :: EXPRVALUE(2)           ! From expr
    integer :: FILE                     ! File name string index
    integer :: FILEACCESS               ! DFACC_CREATE or DFACC_RDWR
    integer :: FIELDINDEX               ! Type of field in l2cf line
    character(len=1024) :: FILENAME     ! Output full filename
    character(len=1024) :: FILE_BASE    ! made up of
    integer :: FILETYPE
    logical, dimension(field_first:field_last) :: GOT
    logical :: GOTSOURCE                ! TRUE if already had a source field
    integer :: GSON                     ! Son of son
    integer :: HANDLE                   ! File handle from hdf/hdf-eos
    character(len=1024) :: HDFNAME      ! Output swath/sd name
    integer :: HDFNAMEINDEX             ! String index for output name
    integer :: HDFVERSION               ! 4 or 5
    logical :: ISNEWDIRECT              ! TRUE if not already in database
    integer :: KEYNO                    ! Loop counter, field in l2cf line
    integer :: l2gp_Version
    logical :: lowerOverlap
    ! integer :: LASTFIELDINDEX           ! Type of previous field in l2cf line
    integer :: Me = -1                  ! String index for trace
    integer :: MYFILE                   ! File permission granted for
    logical :: MYMAKEREQUEST            ! Copy of makeRequest
    character(len=FileNameLen) :: MYTHEFILE  ! File permission granted for
    character(len=256), dimension(:), pointer :: NAMEBUFFER
    integer :: NEXTFILE
    integer :: NOSOURCES                ! No. things to output
    integer :: NumPermitted             ! No. things permitted to output
    integer :: NumOutput                ! No. things actually output
    logical, parameter :: OPENHERE = .false.
    character(len=16) :: OPTIONS
    integer :: AscDescModeVECTOR
    integer :: AscDescModeQTYINDX
    integer :: OUTPUTTYPE               ! l_l2gp, l_l2aux, l_l2fwm, l_l2dgg
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character(len=1024) :: PATH         ! path/file_base
    integer :: PCBottom
    integer :: PCTop
    integer :: RETURNSTATUS
    logical :: SINGLE
    logical :: SKIPDGG
    logical :: SKIPDGM
    integer :: SON                      ! A tree node
    integer :: SOURCE                   ! Loop counter
    logical :: upperOverlap

    integer, dimension(:), pointer :: CONVERGVECTORS ! Indices
    integer, dimension(:), pointer :: CONVERGQUANTITIES ! Indices
    integer, dimension(:), pointer :: SOURCEVECTORS ! Indices
    integer, dimension(:), pointer :: SOURCEQUANTITIES ! Indices
    integer, dimension(:), pointer :: PRECISIONVECTORS ! Indices
    integer, dimension(:), pointer :: PRECISIONQUANTITIES ! Indices
    integer, dimension(:), pointer :: QUALITYVECTORS ! Indices
    integer, dimension(:), pointer :: QUALITYQUANTITIES ! Indices
    integer, dimension(:), pointer :: STATUSVECTORS ! Indices
    integer, dimension(:), pointer :: STATUSQUANTITIES ! Indices
    integer, dimension(:), pointer :: DIRECTFILES ! Indices
    type(VectorValue_T), pointer :: CONVERGQTY ! The quantities convergence ratio
    type(VectorValue_T), pointer :: QTY ! The quantity
    type(VectorValue_T), pointer :: PRECQTY ! The quantities precision
    type(VectorValue_T), pointer :: QUALITYQTY ! The quantities quality
    type(VectorValue_T), pointer :: STATUSQTY ! The quantities status
    type(VectorValue_T), pointer :: AscDescModeQTY ! The shared quantity AscDescMode
    type(DirectData_T), pointer  :: thisDirect ! => null()
    type(Vector_T), pointer      :: Vector ! => null()
    real :: TimeIn, TimeSetUp, TimeWriting, timeToClose, TimeOut
    type(L2Metadata_T) :: l2metaData

    ! Executable code
    DEEBUG = ( switchDetail(switches, 'direct') > -1 ) ! .or. SKIPDIRECTWRITES
    SKIPDGG = ( switchDetail(switches, 'skipdgg') > -1 )
    SKIPDGM = ( switchDetail(switches, 'skipdgm') > -1 )
    nullify(thisDirect)

    call trace_begin ( me, "DirectWriteCommand", node, &
      & cond=toggle(gen) .and. switchDetail(switches,'dwreq') > -1 )
    call time_now ( timeIn )
    TimeSetUp = TimeIn
    TimeWriting = TimeIn
    timeToClose = TimeIn
    timeOut = TimeIn
    myMakeRequest = .false.
    if ( present ( makeRequest ) ) myMakeRequest = makeRequest
    myTheFile = 'undefined'   ! 0
    if ( present ( theFile ) ) myTheFile = theFile
    if ( present(noExtraWrites) ) noExtraWrites = 0

    ! lastFieldIndex = 0
    noSources = 0
    hdfVersion = DEFAULT_HDFVERSION_WRITE
    error = 0
    file = 0
    filename = 'undefined'
    lowerOverlap = .false.
    got = .false.
    options = ' '
    outputType=0
    outputtypestr = 'unknown'
    if ( present(namedFile) ) namedFile = .false.
    single = .false.
    upperOverlap = .false.
    gotsource = .false.
    do keyNo = 2, nsons(node)           ! Skip DirectWrite command
      son = subtree ( keyNo, node )
      ! if ( keyNo > 2 ) lastFieldIndex = fieldIndex
      fieldIndex = get_field_id ( son )
      got ( fieldIndex ) = .true.
      select case ( fieldIndex )
      case ( f_source, f_vector )
        noSources = noSources + 1
        gotsource = .true.
      case ( f_AscDescMode )
        gson = subtree(2,son)
        AscDescModeVector = decoration(decoration(subtree(1,gson)))
        AscDescModeQtyIndx = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_precision )
        if ( .not. gotsource ) &
          & call Announce_Error ( son, no_error_code, &
            & 'A precision can only be given following a source' )
      case ( f_quality )
        if ( .not. gotsource ) &
          & call Announce_Error ( son, no_error_code, &
            & 'A quality can only be given following a source' )
      case ( f_status )
        if ( .not. gotsource ) &
          & call Announce_Error ( son, no_error_code, &
            & 'A status can only be given following a source' )
      case ( f_hdfVersion )
        call expr ( subtree(2,son), exprUnits, exprValue )
        if ( exprUnits(1) /= phyq_dimensionless ) &
          & call Announce_error ( son, NO_ERROR_CODE, &
          & 'No units allowed for hdfVersion: just integer 4 or 5')
        hdfVersion = exprValue(1)
      case ( f_options )
        call get_string ( sub_rosa(subtree(2,son)), options, strip=.true. )
      case ( f_file )
        file = sub_rosa(subtree(2,son))
        if ( present(namedFile) ) namedFile = .true.
      case ( f_type )
        outputType = decoration(subtree(2,son))
        call get_string ( lit_indices(outputType), outputTypeStr, strip=.true. )
      case ( f_lowerOverlap )
        lowerOverlap = get_boolean ( son )
      case ( f_upperOverlap )
        upperOverlap = get_boolean ( son )
      case ( f_single )
        single = get_boolean ( son )
      end select
    end do

    if ( file > 0 ) then
      call get_string ( file, filename, strip=.true. )
    else if ( myTheFile /= 'undefined' ) then
      filename = myTheFile
    end if
    myFile =  FileNameToID(trim(filename), DirectDataBase ) 
    distributingSources = (file < 1)
    ! Now identify the quantities we're after
    nullify ( sourceVectors, sourceQuantities, &
      & qualityVectors, qualityQuantities, &
      & convergQuantities, convergVectors, &
      & statusVectors, statusQuantities, &
      & precisionVectors, precisionQuantities, directFiles )
    call Allocate_test ( sourceVectors, noSources, 'sourceVectors', ModuleName )
    call Allocate_test ( sourceQuantities, noSources, 'sourceQuantities', &
      &  ModuleName, Fill=-1 )
    call Allocate_test ( precisionVectors, noSources, 'precisionVectors', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( precisionQuantities, noSources, 'precisionQuantities', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( qualityVectors, noSources, 'qualityVectors', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( qualityQuantities, noSources, 'qualityQuantities', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( statusVectors, noSources, 'statusVectors', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( statusQuantities, noSources, 'statusQuantities', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( convergVectors, noSources, 'convergVectors', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( convergQuantities, noSources, 'convergQuantities', &
      &  ModuleName, Fill=0 )
    call Allocate_test ( directFiles, noSources, 'directFiles', &
      &  ModuleName, Fill=0 )
    ! Go round again and identify each quantity, work out what kind of file
    ! we're talking about
    source = 0
    do keyNo = 2, nsons(node)
      l2gp_Version = 1
      son = subtree ( keyNo, node )
      fieldIndex = get_field_id ( son )
      select case ( fieldIndex )
      case ( f_source )
        source = source + 1
        gson = subtree(2,son)
        sourceVectors(source) = decoration(decoration(subtree(1,gson)))
        sourceQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_convergence )
        if ( all ( outputType /= (/ l_l2gp, l_l2dgg /) ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Convergence only appropriate for l2gp files" )
        gson = subtree(2,son)
        convergVectors(source) = decoration(decoration(subtree(1,gson)))
        convergQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_precision )
        if ( all ( outputType /= (/ l_l2gp, l_l2dgg /) ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Precision only appropriate for l2gp files" )
        gson = subtree(2,son)
        precisionVectors(source) = decoration(decoration(subtree(1,gson)))
        precisionQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_quality )
        if ( all ( outputType /= (/ l_l2gp, l_l2dgg /) ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Quality only appropriate for l2gp files" )
        gson = subtree(2,son)
        qualityVectors(source) = decoration(decoration(subtree(1,gson)))
        qualityQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_status )
        if ( all ( outputType /= (/ l_l2gp, l_l2dgg /) ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Status only appropriate for l2gp files" )
        gson = subtree(2,son)
        statusVectors(source) = decoration(decoration(subtree(1,gson)))
        statusQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_vector )
        source = source + 1
        sourceVectors(source) = decoration(decoration(subtree(2,son)))
      case default
      end select
    end do

    ! if ( .not. (SKIPDIRECTWRITES .or. checkpaths) ) then
    if ( .not. checkpaths ) then
    ! Now go through and do some sanity checking
    nullify ( AscDescModeQty )
    if ( got(f_AscDescMode) ) AscDescModeQty => GetVectorQtyByTemplateIndex ( vectors(AscDescModevector), &
    & AscDescModeQtyIndx )
    do source = 1, noSources
      if ( sourceQuantities(source) < 1 ) cycle
      qty => GetVectorQtyByTemplateIndex ( vectors(sourceVectors(source)), &
      & sourceQuantities(source) )
      if ( qty%label == 0 ) call Announce_Error ( son, no_error_code, &
      & "Quantity does not have a label" )
      if ( precisionVectors(source) /= 0 ) then
        precQty => &
          & GetVectorQtyByTemplateIndex ( vectors(precisionVectors(source)), &
          & precisionQuantities(source) )
        ! Check that this is compatible with its value quantity
        if ( qty%template%name /= precQty%template%name ) &
        & call Announce_Error ( son, no_error_code, &
        & "Precision and quantity do not match" )
      else
        precQty => NULL()
      end if
      if ( qualityVectors(source) /= 0 ) then
        qualityQty => &
          & GetVectorQtyByTemplateIndex ( vectors(qualityVectors(source)), &
          & qualityQuantities(source) )
        ! Check that value and quality share same HGrid
        if ( .not. DoHgridsMatch( qty, qualityQty ) ) &
        & call Announce_Error ( son, no_error_code, &
        & "Source and quality not on matching HGrids" )
      else
        qualityQty => NULL()
      end if
      if ( statusVectors(source) /= 0 ) then
        statusQty => &
          & GetVectorQtyByTemplateIndex ( vectors(statusVectors(source)), &
          & statusQuantities(source) )
        ! Check that value and quality share same HGrid
        if ( .not. DoHgridsMatch( qty, statusQty ) ) &
        & call Announce_Error ( son, no_error_code, &
        & "Source and status not on matching HGrids" )
      else
        statusQty => NULL()
      end if
      if ( convergVectors(source) /= 0 ) then
        convergQty => &
          & GetVectorQtyByTemplateIndex ( vectors(convergVectors(source)), &
          & convergQuantities(source) )
        ! Check that value and convergence share same HGrid
        if ( .not. DoHgridsMatch( qty, convergQty ) ) &
        & call Announce_Error ( son, no_error_code, &
        & "Source and convergence not on matching HGrids" )
      else
        convergQty => NULL()
      end if
      ! Now check that things make sense
      if ( ValidateVectorQuantity ( qty, &
        & coherent=.true., stacked=.true., regular=.true., &
        & verticalCoordinate = (/ l_pressure, l_zeta, l_none/), &
        & minorFrame=.false., majorFrame=.false. ) ) then
        expectedType = l_l2gp
      else
        expectedType = l_l2aux
      end if
      if ( outputType /= expectedType .and. .not. &
        & ( outputType == l_l2dgg .and. expectedType == l_l2gp ) &
        &                           .and. .not.  &
        & ( outputType == l_quantity ) &
        &                           .and. .not.  &
        & ( outputType == l_l2fwm .and. expectedType == l_l2aux ) ) then
        call output ( "Offending quantity " )
        call display_string ( qty%template%name, strip=.true., advance='yes' )
        call output ( "Expected type " )
        call display_string ( lit_indices(expectedType), advance='yes' )
        call output ( "Output file type " )
        call display_string ( lit_indices(outputType), advance='yes' )
        call Announce_Error ( son, no_error_code, &
          & "Inappropriate quantity for this file type in direct write", &
          & Penalty=0 )
      end if
    end do
    end if
    
    ! Bail out at this stage if there is some kind of error.
    if ( error /= 0 ) then
      call trace_end ( "DirectWriteCommand", &
        & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
      return
    endif

    ! Distribute sources among available DirectWrite files if filename undefined
    if ( distributingSources ) then
       call DistributeSources
    end if
    if ( DeeBUG ) then
      call output('Direct write to file', advance='yes')
      call output('File name: ', advance='no')
      call output(trim(filename), advance='yes')
      call output('hdfVersion: ', advance='no')
      call output(hdfVersion, advance='yes')
      call output('Num sources: ', advance='no')
      call output(noSources, advance='yes')
      if ( present(theFile) ) then
        call output('my(theFile): ', advance='no')
        call output(trim(myTheFile), advance='yes')
      end if
      call output('size(DirectDB): ', advance='no')
      call output(size(DirectDatabase), advance='yes')
    end if

    ! If we're skipping all directwrites, let's deallocate and return
    if ( SKIPDIRECTWRITES ) then
      call DeallocateStuff
      call trace_end ( "DirectWriteCommand", &
        & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
      return
    endif

    ! If this is the first pass through, then we just log our request
    ! with the master
    if ( myMakeRequest ) then
      if ( file < 1 ) then
        ! Here's a crude idea: loop until we've wrapped
        ! (A better idea would pick out only unique integers from an array)
        nextFile = directFiles(1)
        call LogDirectWriteRequest ( DirectDatabase(nextfile)%fileIndex, node )
        if ( DEEBug ) then
          call output('Logged directwrite rq for: ', advance='no')
          call output(DirectDatabase(nextfile)%fileIndex, advance='no')
          call output('   node: ', advance='no')
          call output(node, advance='no')
          call output('  (named): ', advance='no')
          call display_string(DirectDatabase(nextfile)%fileIndex, advance='yes')
        end if
        if ( noSources > 1 ) then
          do source = 2, noSources
            aFile = nextFile
            nextFile = directFiles(source)
            if ( nextFile < aFile ) exit
            call LogDirectWriteRequest ( DirectDatabase(nextfile)%fileIndex, node )
            if ( present(noExtraWrites) ) noExtraWrites = noExtraWrites + 1
            if ( DEEBug ) then
              call output('Logged directwrite rq for: ', advance='no')
              call output(DirectDatabase(nextfile)%fileIndex, advance='no')
              call output('   node: ', advance='no')
              call output(node, advance='no')
              call output('  (named): ', advance='no')
              call display_string(DirectDatabase(nextfile)%fileIndex, advance='yes')
            end if
          end do
        end if
        if ( DEEBug .and. present(noExtraWrites) ) then
          call output('noExtraWrites: ', advance='no')
          call output(noExtraWrites, advance='no')
        end if
      else
        call LogDirectWriteRequest ( file, node )
      end if
    else if ( skipdgg .and. any ( outputType == (/ l_l2dgg /) ) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'DirectWriteCommand skipping all dgg writes ' // trim(filename) )
      call DeallocateStuff
      call trace_end ( "DirectWriteCommand", &
        & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
      return
    else if ( skipdgm .and. any ( outputType == (/ l_l2fwm, l_l2aux, l_hdf /) ) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'DirectWriteCommand skipping all dgm/fwm writes ' // trim(filename) )
      call DeallocateStuff
      call trace_end ( "DirectWriteCommand", &
        & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
      return
    else
      ! OK, it's time to write this bit of the file
      if ( parallel%slave ) then
        createFileFlag = .false.
        if ( present ( create ) ) createFileFlag = create
      else if ( distributingSources ) then
        ! Short-circuit this direct write if outputType fails to match
        if ( myFile < 1 ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to auto-direct write ' // trim(filename) )
        else if ( outputType /= DirectDataBase(myFile)%type ) then
          call DeallocateStuff
          if ( DeeBUG ) print *, 'Short-circuiting ' // trim(filename)
          call trace_end ( "DirectWriteCommand", &
            & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
          return
        end if
        createFileFlag = .not. any ( createdFilenames == myFile )
        if ( createFileFlag ) then
          noCreatedFiles = noCreatedFiles + 1
          if ( noCreatedFiles > maxFiles ) call MLSMessage ( &
            & MLSMSG_Error, ModuleName, 'Too many direct write files' )
          createdFilenames ( noCreatedFiles ) = myFile
        end if
      else
        createFileFlag = .not. any ( createdFilenames == file )
        if ( createFileFlag ) then
          noCreatedFiles = noCreatedFiles + 1
          if ( noCreatedFiles > maxFiles ) call MLSMessage ( &
            & MLSMSG_Error, ModuleName, 'Too many direct write files' )
          createdFilenames ( noCreatedFiles ) = file
        end if
      end if
      
      returnStatus = 0
      ! Open/create the file of interest
      isnewdirect = .true.  ! Just so it has a value even for toolkitless runs
      call split_path_name(filename, path, file_base)
      if ( distributingSources .or. TOOLKIT .or. CATENATESPLITS ) then
        call ExpandDirectDB ( DirectDatabase, file_base, thisDirect, &
        & isnewdirect )
        if ( DeeBUG ) then
          call output('Did we need to expand DB for ' // trim(file_base), &
            & advance='yes')
          call output(isnewdirect, advance='yes')
        end if
        if ( .not. associated(thisDirect) ) then
          call Announce_Error ( son, NO_ERROR_CODE, &
              & 'ExpandDirectDB returned unassociated thisDirect' )
          call dump(DirectDatabase)
          call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'ExpandDirectDB returned unassociated thisDirect' )
        end if
      end if
      if ( .not. (TOOLKIT .or. CATENATESPLITS) ) then
        handle = 0
      else if ( .not. isnewdirect ) then
        Filename = thisDirect%fileName
        Handle = thisDirect%Handle
      else if ( .not. TOOLKIT ) then
        Handle = -1
      else if ( any ( outputType == (/ l_l2gp /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
          & mlspcf_l2gp_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      else if ( any ( outputType == (/ l_l2dgg /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2dgg_start, &
          & mlspcf_l2dgg_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      else if ( any ( outputType == (/ l_l2fwm /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2fwm_full_start, &
          & mlspcf_l2fwm_full_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      else
        Handle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
          & mlspcf_l2dgm_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      end if
      if ( returnStatus /= 0 ) then
          call Announce_Error ( node, NO_ERROR_CODE, &
          & 'Uh-oh' )
         call MLSMessage ( &
         & MLSMSG_Error, ModuleName, &
         & 'Failed in GetPCFromRef for ' // trim(filename) )
      end if
      if ( isnewdirect .and. (TOOLKIT .or. CATENATESPLITS) ) then
        thisDirect%Handle = Handle
        thisDirect%FileName = FileName
      end if
      if ( isnewdirect .and. distributingSources ) then
          call Announce_Error ( node, NO_ERROR_CODE, &
          & 'Uh-oh' )
        call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'ExpandDirectDB thinks we need a new directFile; ' // &
          & 'did you enter any in global settings?' )
      end if
      ! Done what we wished to do if just checking paths or SKIPDIRECTWRITES
      if ( checkPaths ) then
      ! if ( SKIPDIRECTWRITES .or. checkPaths ) then
        call trace_end ( "DirectWriteCommand", &
          & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
        return
      endif
      
      if ( createFileFlag .and. .not. patch ) then
        fileaccess = DFACC_CREATE
      else
        fileaccess = DFACC_RDWR
      end if
      if ( DeeBUG ) then
        print *, 'Trying to open ', trim(FileName)
        print *, 'FileAccess ', FileAccess
        print *, 'hdfVersion ', hdfVersion
        print *, 'outputType ', outputType
        print *, 'myFile ', myFile
        print *, 'createFileFlag ', createFileFlag
      end if
      if ( .not. isnewdirect ) then
        ! if ( outputType /= thisDirect%type ) then
        if ( .not. doOutputTypesMatch( outputType, thisDirect%type ) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'DirectWriteCommand mismatched outputTypes for ' // trim(filename) )
        end if
      end if

      ! Some indicators of where to find the file in the PCF
      ! and what type of file it will be (fileType will be redefined later)
      select case ( outputType )
      case ( l_l2gp )
        PCBottom = mlspcf_l2gp_start
        PCTop    = mlspcf_l2gp_end
        fileType = l_swath
      case ( l_l2dgg )
        PCBottom = mlspcf_l2dgg_start
        PCTop    = mlspcf_l2dgg_end
        fileType = l_swath
      case ( l_l2fwm )
        PCBottom = mlspcf_l2fwm_full_start
        PCTop    = mlspcf_l2fwm_full_end
        fileType = l_hdf
      case ( l_l2aux, l_hdf, l_quantity )
        PCBottom = mlspcf_l2dgm_start
        PCTop    = mlspcf_l2dgm_end
        fileType = l_hdf
      end select

      ! Call the l2gp open/create routine.  Filename is 'filename'
      ! file id should go into 'handle'
      directFile => GetMLSFileByName(filedatabase, filename, ignore_paths=.true.)
      if ( .not. associated(directFile) ) then
        if(DEEBUG) call MLSMessage(MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(filename) )
        directFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=Filename, shortName=file_base, &
          & type=fileType, access=DFACC_CREATE, HDFVersion=HDFVERSION_5, &
          & PCBottom=PCBottom, PCTop=PCTop)
      endif
      directFile%access = FileAccess
      if ( OPENHERE ) call mls_openFile(directFile, ErrorType)
      if(DEEBUG) call dump(directFile)

      errorType = 0
      select case ( outputType )
      case ( l_l2gp, l_l2dgg )
        ! Before opening file, see which swaths are already there
        ! and which ones need to be created
        if ( DeeBUG ) print *, 'Allocating ', noSources
        nullify(createThisSource, nameBuffer)
        call Allocate_test ( createThisSource, noSources, 'createThisSource', &
          & ModuleName )
        call Allocate_test ( nameBuffer, noSources, 'nameBuffer', &
          & ModuleName )
        if(DEEBUG)print *, 'Must we create this file?', createFileFlag
        if ( .not. createFileFlag ) then
          do source = 1, noSources
            if(DEEBUG)print*,'Source:', source
            if ( sourceQuantities(source) > 0 ) then
              qty => GetVectorQtyByTemplateIndex ( &
                & vectors(sourceVectors(source)), sourceQuantities(source) )
            else
              qty => GetVectorQtyByTemplateIndex ( &
                & vectors(sourceVectors(source)), 1 )
            endif
            hdfNameIndex = qty%label
            if(DEEBUG) &
              & call display_string ( hdfNameIndex, strip=.true., advance='yes' )
            call get_string ( hdfNameIndex, nameBuffer(source), strip=.true. )
            if(DEEBUG)print*,'Done'
          end do
          if(DEEBUG)print *, '************** Checking for swaths in file ***************'
          if(DEEBUG)call dump(directFile)
          dummy = MLS_SWATH_IN_FILE(directFile%Name, nameBuffer, HdfVersion, &
            & createThisSource, returnStatus )
          if(DEEBUG)print*,'Got out of MLS_SWATH_IN_FILE'
          if ( returnStatus /= 0 ) then
            call MLSMessage(MLSMSG_Warning, ModuleName, &
              & 'Unable to check on swath in ' // trim(filename) )
          endif
          if(DEEBUG)call dump( createThisSource, 'createThisSource' )
          source = findFirst( createThisSource )
          if(DEEBUG)call outputNamedValue ( 'source number of T', source )
        else
          createThisSource = .false.
        end if
      case ( l_l2aux, l_l2fwm, l_hdf, l_quantity )
        ! Nothing special
        ! (Why don't we need to know which SDs are there and which aren't?)
      case default
          call MLSMessage(MLSMSG_Warning, ModuleName, &
          & 'Unrecognized output type ' // trim(directFile%name), &
          & MLSFile=directFile )
      end select
      
      if ( ErrorType /= 0 .and. OPENHERE ) then
        print *, 'Tried to open ', trim(directFile%Name)
        print *, 'ErrorType ', ErrorType
        print *, 'Output Type ', OutputType
        print *, 'Output Type ', OutputTypeStr
        print *, 'FileAccess ', directFile%Access
        print *, 'hdfVersion ', directFile%hdfVersion
        print *, 'myFile ', myFile
        print *, 'createFileFlag ', createFileFlag
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to open ' // trim(filename), &
          & MLSFIle=directFile )
      end if
      ! Loop over the quantities to output
      NumPermitted = 0
      NumOutput = 0
      if(DEEBUG)print *, '************** Loop over quantities to output ***************'
      do source = 1, noSources
        ! At this point we're ready to write the data
        ! Make certain that the sources match the permitted file
        ! if the sources are being distributed, however
        if ( distributingSources ) then
          if ( DirectDataBase(directfiles(source))%filename /= myTheFile &
            & .and. &
            & DirectDataBase(directfiles(source))%filenameBase /= myTheFile &
            & ) cycle
        end if
        NumPermitted = NumPermitted + 1
        if ( sourceQuantities(source) < 1 ) then
          ! This is the case where we write an entire vector
          vector => vectors(sourceVectors(source))
          select case ( outputType )
          case ( l_l2gp, l_l2dgg )
            call DirectWrite ( directFile, vector, &
              & chunkNo, chunks, FWModelConfig, &
              & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
              & single=single, options=options )
          case ( l_l2aux, l_l2fwm, l_hdf )
            call DirectWrite ( directFile, vector, &
              & chunkNo, chunks, FWModelConfig, &
              & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
              & single=single, options=options )
          case ( l_quantity )
            call DirectWrite ( directFile, &
              & vector, &
              & chunkNo, options=options )
          case default
          end select
          cycle
        endif
        qty => GetVectorQtyByTemplateIndex ( vectors(sourceVectors(source)), &
          & sourceQuantities(source) )
        hdfNameIndex = qty%label
        call get_string ( hdfNameIndex, hdfName, strip=.true. )
        if ( TOOLKIT .or. CATENATESPLITS ) &
          & call ExpandSDNames(thisDirect, trim(hdfName))
        if ( (TOOLKIT .or. CATENATESPLITS) .and. DEEBUG ) call dump(thisDirect)
        if ( precisionVectors(source) /= 0 ) then
          precQty => GetVectorQtyByTemplateIndex &
            & ( vectors(precisionVectors(source)), precisionQuantities(source) )
          ! Check that this is compatible with its value quantitiy
          if ( qty%template%name /= precQty%template%name ) &
            & call Announce_Error ( son, no_error_code, &
            & "Precision and quantity do not match" )
        else
          precQty => NULL()
        end if
        if ( qualityVectors(source) /= 0 ) then
          qualityQty => &
            & GetVectorQtyByTemplateIndex ( vectors(qualityVectors(source)), &
            & qualityQuantities(source) )
          ! Check that value and quality share same HGrid
          if ( .not. DoHgridsMatch( qty, qualityQty ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Source and quality not on matching HGrids" )
        else
          qualityQty => NULL()
        end if
        if ( statusVectors(source) /= 0 ) then
          statusQty => &
            & GetVectorQtyByTemplateIndex ( vectors(statusVectors(source)), &
            & statusQuantities(source) )
          ! Check that value and quality share same HGrid
          if ( .not. DoHgridsMatch( qty, statusQty ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Source and status not on matching HGrids" )
        else
          statusQty => NULL()
        end if
        if ( convergVectors(source) /= 0 ) then
          convergQty => &
            & GetVectorQtyByTemplateIndex ( vectors(convergVectors(source)), &
            & convergQuantities(source) )
          ! Check that value and convergence share same HGrid
          if ( .not. DoHgridsMatch( qty, convergQty ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Source and convergence not on matching HGrids" )
        else
          convergQty => NULL()
        end if
        
        if ( DeeBUG ) then
          call output('CreateFileFlag: ', advance='no')
          call output(createFileFlag, advance='yes')
          call output('file access: ', advance='no')
          call output(fileaccess, advance='yes')
          call output('file handle: ', advance='no')
          call output(handle, advance='yes')
          call output('outputType: ', advance='no')
          call output(outputType, advance='yes')
          call output('sd name: ', advance='no')
          call output(trim(hdfname), advance='yes')
        end if
        call time_now ( timeSetup )
        if ( timeSetup-timeIn > timeReasonable .and. &
          & switchDetail(switches,'dwreq') > -1 ) then
          call output('Unreasonable set up time for ' // trim(hdfname), &
            & advance='yes')
        endif

        ! Do the actual DirectWrite
        ! (Why do we need to redefine fileType? Is add_metadata so stupid?)
        select case ( outputType )
        case ( l_l2gp, l_l2dgg )
          ! Call the l2gp swath write routine.  This should write the 
          ! non-overlapped portion of qty (with possibly precision in precQty)
          ! into the l2gp swath named 'hdfName' starting at profile 
          ! qty%template%instanceOffset + 1
          ! May optionally supply first, last profiles
          if ( DEEBUG) then
            call dump(directFile, details=1)
            call output('createSwath: ', advance='no')
            call output(.not. createThisSource(source), advance='yes')
            call outputNamedValue ( 'source number of DW', source )
          endif
          createthisswath = (.not. createThisSource(source))
          ! We had a bug somewhere in hdfeos
          ! When we created the first swath in an hdfeos file
          ! mls_swath_in_file didn't find it
          ! What we'll try is to create the swath, then close the file
          ! and reset its access to read/write--could be that it was confused 
          ! by the DFACC_CREATE
          call DirectWrite ( directFile, &
            & qty, precQty, qualityQty, statusQty, convergQty, AscDescModeQty, &
            & hdfName, chunkNo, HGrids, &
            & createSwath=createthisswath, &
            & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
          if ( createthisswath ) then
            if ( directFile%stillOpen ) &
              & call mls_closeFile(directFile, errorType)
            directFile%access = DFACC_RDWR
          endif
          if ( fileaccess == DFACC_CREATE ) then
            ! OK, because the bug is still there (!), we'll repeat
            call DirectWrite ( directFile, &
              & qty, precQty, qualityQty, statusQty, convergQty, AscDescModeQty, &
              & hdfName, chunkNo, HGrids, &
              & createSwath=createthisswath, &
              & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
          endif
          NumOutput= NumOutput + 1
          if ( outputType == l_l2dgg ) then
            filetype=l_l2dgg
          else
            filetype=l_swath
          end if
        case ( l_l2aux, l_l2fwm, l_hdf )
          ! Call the l2aux sd write routine.  This should write the 
          ! non-overlapped portion of qty (with possibly precision in precQty)
          ! into the l2aux sd named 'hdfName' starting at profile 
          ! qty%template%instanceOffset ( + 1 ? )
          ! Note sure about the +1 in this case, probably depends whether it's a
          ! minor frame quantity or not.  This mixed zero/one indexing is becoming
          ! a real pain.  I wish I never went down that road!
          call DirectWrite ( directFile, qty, precQty, hdfName, &
            & chunkNo, chunks, FWModelConfig, &
            & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
            & single=single, options=options )
          NumOutput = NumOutput + 1
          filetype=l_hdf
        case ( l_quantity )
          ! Write the quantity with all its geolocations
          call DirectWrite ( directFile, qty, hdfName, &
            & chunkNo, options=options )
          NumOutput = NumOutput + 1
          filetype=l_quantity
        case default
          call output('outputType: ', advance='no')
          call output(outputType, advance='yes')
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unrecognized OutputType in ' // trim(filename) )
        end select
      end do ! End loop over swaths/sds
      call time_now ( timeWriting )
      if ( timeWriting-timeSetup > timeReasonable .and. &
        & switchDetail(switches,'dwreq') > -1 ) then
        call output('Unreasonable writing time for ' //trim(hdfname), advance='yes')
      endif
      
      if ( DEEBUG ) then
        print *, 'Num permitted to ', trim(FileName), ' ', NumPermitted
        print *, 'Num actually output ', NumOutput
      end if
      if ( NumPermitted < 1 ) then
        if ( parallel%slave) then
          call Announce_Error ( son, no_error_code, &
            & "NumPermitted=0", penalty=0 )
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No sources permitted for writing to  ' // trim(filename) )
          call FinishedDirectWrite ( ticket )
        else
          ! But did we claim we were created this file? If so, decrement
          if ( createFileFlag ) noCreatedFiles = noCreatedFiles - 1
          if ( DEEBUG ) then
            print *, 'No sources written to ', trim(FileName)
            print *, 'FileAccess ', FileAccess
            print *, 'hdfVersion ', hdfVersion
            print *, 'myFile ', myFile
            print *, 'createFileFlag ', createFileFlag
            print *, 'noCreatedFiles ', noCreatedFiles
          end if
        end if
        call DeallocateStuff
        if ( any ( outputType == (/ l_l2gp, l_l2dgg /) ) ) then
          call Deallocate_test ( createThisSource, 'createThisSource', ModuleName )
          call Deallocate_test ( nameBuffer, 'nameBuffer', ModuleName )
        end if

        ! Don't forget to close file
        select case ( outputType )
        case ( l_l2gp, l_l2dgg, l_l2aux, l_l2fwm, l_hdf )
          if ( OPENHERE ) call mls_closeFile(directFile, ErrorType)
        case default
          call output('outputType: ', advance='no')
          call output(outputType, advance='yes')
          call output('ErrorType: ', advance='no')
          call output(ErrorType, advance='yes')
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Tried to close Unrecognized OutputType in ' // trim(filename), &
            & MLSFile=directFile )
        end select
        call time_now ( timeToClose )
        if ( timeToClose-timeWriting > timeReasonable .and. &
          & switchDetail(switches,'dwreq') > -1 ) then
          call output('Unreasonable closing time for ' // trim(hdfname), &
            & advance='yes')
        endif
        if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to close (1)' // trim(filename), &
          & MLSFile=directFile )
        call trace_end ( "DirectWriteCommand", &
          & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
        return
      end if
      if ( isnewdirect .and. (TOOLKIT .or. CATENATESPLITS) ) then
        thisDirect%type = outputType
        thisDirect%fileNameBase = file_base
      end if
      
      if ( DEEBug ) then
        call output('outputType: ', advance='no')
        call DISPLAY_STRING(lit_indices(outputType), advance='yes')
        call output('fileType: ', advance='no')
        call DISPLAY_STRING(lit_indices(fileType), advance='yes')
      end if

      ! Close the output file of interest (does this need to be split like this?)
      select case ( outputType )
      case ( l_l2gp, l_l2dgg )
        if ( DeeBUG ) print *, 'Deallocating ', noSources
        call Deallocate_test ( createThisSource, 'createThisSource', ModuleName )
        call Deallocate_test ( nameBuffer, 'nameBuffer', ModuleName )
        ! Call the l2gp close routine
        if ( OPENHERE ) call mls_closeFile(directFile, errorType)
        if ( DeeBUG ) then
          print *, 'Tried to close ', trim(FIleName)
          print *, 'Handle ', Handle
          print *, 'hdfVersion ', hdfVersion
          print *, 'errortype ', errortype
        end if
      case ( l_l2aux, l_l2fwm, l_hdf, l_quantity )
        ! Call the l2aux close routine
        if ( OPENHERE ) call mls_closeFile(directFile, errorType)
      case default
        call DISPLAY_STRING(lit_indices(outputType), before='outputType: ', advance='yes')
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Tried to closeUnrecognized output type ' // trim(filename), &
        & MLSFile=directFile )
      end select
      call time_now ( timeToClose )
      if ( timeToClose-timeWriting > timeReasonable .and. &
        & switchDetail(switches,'dwreq') > -1 ) then
        call output('Unreasonable closing time for ' // trim(hdfname), &
          & advance='yes')
      endif
      if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'DirectWriteCommand unable to close (2)' // trim(filename), &
        & MLSFile=directFile )
      if ( createFileFlag .and. TOOLKIT .and. .not. SKIPMETADATA .and. &
        & .not. any( outputType == (/l_l2fwm, l_hdf /)) .and. &
        & .not. (distributingSources .and. CATENATESPLITS) ) then
        call add_metadata ( node, file_base, l2metaData, &
          & hdfVersion, filetype, errortype, NumPermitted, thisDirect%sdNames )
        if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to addmetadata to ' // trim(filename), &
          & MLSFile=directFile )
      end if
      
      ! Tell the master we're done
      if ( parallel%slave ) call FinishedDirectWrite ( ticket )
      call time_now ( timeOut )
      if ( timeOut-timeToClose > timeReasonable .and. &
        & switchDetail(switches,'dwreq') > -1 ) then
        call output('Unreasonable time out ' //trim(hdfname), advance='yes')
      endif
    end if

    call trace_end ( "DirectWriteCommand", &
      & cond=toggle(gen) .and. switchDetail(switches, 'dwreq') > -1 )
    call DeallocateStuff

  contains
    subroutine DistributeSources
      ! Distribute noSources among the files in the directWrite db
      ! Of course, the types must match
      ! Local variables
      integer :: afile
      integer :: nextfile
      integer :: NumSuitablesFiles   ! How many with matching types
      integer :: source
      ! Executable
      if ( size(DirectDataBase) < 1 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, &
        & 'No files in directDatabase to distribute sources among' )
      NumSuitablesFiles = 0
      do afile = 1, size(DirectDataBase)
        if ( outputType == DirectDataBase(afile)%autoType ) &
          & NumSuitablesFiles = NumSuitablesFiles + 1
      end do
      if ( NumSuitablesFiles < 1 ) then
        call Announce_Error ( son, NO_ERROR_CODE, &
              & 'No suitable files in directDatabaset' )
        call output('outputType: ')
        call output(outputType, advance='yes')
        call display_string ( lit_indices(outputType), advance='yes' )
        call output(DirectDataBase%autoType)
        call MLSMessage ( &
        & MLSMSG_Error, ModuleName, &
        & 'No suitable files in directDatabase to distribute sources among' )
      end if
      nextfile = FindFirst(DirectDataBase%autoType, outputType)
      do source = 1, noSources
        afile = nextFile
        directfiles(source) = afile
        nextFile = FindNext(DirectDataBase%autoType, outputType, afile, &
          & wrap=.true.)
      end do
      if ( DeeBug ) then
       call output ( "NumSources = " )
       call output ( NoSources, advance='no' )
       call output ( "    NumSuitablesFiles = " )
       call output ( NumSuitablesFiles, advance='yes' )
       do source=1, noSources
         aFile= directfiles(source)
         call output ( source, after='    ', advance='no' )
         call output ( aFile, after='    ', advance='no' )
         call output ( trim(DirectDatabase(aFile)%fileName), advance='no' )
         call blanks ( 4, advance='no' )
         call display_string(DirectDatabase(aFile)%fileIndex, advance='yes')
       end do
      end if
    end subroutine DistributeSources
    subroutine DeallocateStuff
      call Deallocate_test ( sourceVectors, 'sourceVectors', ModuleName )
      call Deallocate_test ( sourceQuantities, 'sourceQuantities', ModuleName )
      call Deallocate_test ( precisionVectors, 'precisionVectors', ModuleName )
      call Deallocate_test ( precisionQuantities, 'precisionQuantities', ModuleName )
      call Deallocate_test ( qualityVectors, 'qualityVectors', ModuleName )
      call Deallocate_test ( qualityQuantities, 'qualityQuantities', ModuleName )
      call Deallocate_test ( statusVectors, 'statusVectors', ModuleName )
      call Deallocate_test ( statusQuantities, 'statusQuantities', ModuleName )
      call Deallocate_test ( convergVectors, 'convergVectors', ModuleName )
      call Deallocate_test ( convergQuantities, 'convergQuantities', ModuleName )
      call Deallocate_test ( directFiles, 'directFiles', ModuleName )
    end subroutine DeallocateStuff
    function doOutputTypesMatch( type1, type2 ) result( match )
      ! Args
      integer, intent(in) :: type1, type2
      logical :: match
      select case ( type1 )
      case ( l_l2gp, l_l2dgg )
        match = any ( type2 == (/ l_l2gp, l_l2dgg /) )
      case ( l_l2aux, l_hdf )
        match = any ( type2 == (/ l_l2aux, l_hdf /) )
      case ( l_l2fwm )
        match = any ( type2 == (/ l_l2fwm, l_hdf /) )
      case default
        match = .false.
      end select
    end function doOutputTypesMatch
  end subroutine DirectWriteCommand

  ! ------------------------------------------------ LabelVectorQuantity -----
  subroutine LabelVectorQuantity ( node, vectors )
    use INIT_TABLES_MODULE, only: F_LABEL, F_PREFIXSIGNAL, &
      & F_QUANTITY, F_VECTOR
    use MLSSIGNALS_M, only: GETSIGNALNAME
    use MORETREE, only: GET_FIELD_ID, GET_BOOLEAN
    use SYMBOL_TABLE, only: ENTER_TERMINAL
    use SYMBOL_TYPES, only: T_STRING
    use STRING_TABLE, only: GET_STRING
    use TREE, only: NSONS, SUBTREE, SUB_ROSA, DECORATION
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T, &
      & GETVECTORQUANTITY, GETVECTORQTYBYTEMPLATEINDEX
    ! Dummy arguments
    integer, intent(in) :: NODE          ! Tree node for l2cf line
    type (Vector_T), dimension(:), pointer :: VECTORS ! Vectors database
    ! Local variables
    integer :: FIELDINDEX               ! Type of field
    integer :: KEYNO                    ! Field index
    integer :: LABEL                    ! String index
    integer :: QUANTITYINDEX            ! Index into quantities database
    integer :: SON                      ! Tree node
    integer :: SOURCE                   ! Tree node
    integer :: VECTORINDEX              ! Index into database
    logical :: PREFIXSIGNAL             ! From l2cf
    ! logical :: SUFFIXLABEL              ! From l2cf
    integer :: VLABEL                   ! String index
    character(len=8) :: whatToLabel     ! 'quantity' or 'vector'
    ! Executable code
    type (VectorValue_T), pointer :: QTY ! The quantity
    character(len=1024) :: LABELSTR     ! The label itself
    ! character(len=128)  :: QTYSTR       ! The template name

    prefixSignal = .false.
    ! Loop over the fields of the mlscf line
    do keyNo = 2, nsons(node) ! Skip spec name
      son = subtree(keyNo,node)
      fieldIndex = get_field_id(son)
      select case ( fieldIndex )
      case ( f_quantity )
        source = subtree(2,son) ! required to be an n_dot vertex
        vectorIndex = decoration(decoration(subtree(1,source)))
        quantityIndex = decoration(decoration(decoration(subtree(2,source))))
        whatToLabel = 'quantity'
      case ( f_vector )
        vectorIndex = decoration(decoration(subtree(2,son)))
        whatToLabel = 'vector'
      case ( f_prefixSignal )
        prefixSignal = get_boolean ( son )
      ! case ( f_suffixLabel )
        ! suffixLabel = get_boolean ( son )
      case ( f_label )
        label = sub_rosa(subtree(2,son))
      case default ! Can't get here if tree_checker worked properly
      end select
    end do

    select case ( whatToLabel )
    case ( 'quantity' )
      ! Get the quantity
      qty => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )

      ! Adapt the label if the prefix signal flag is set.
      if ( prefixSignal ) then
        if ( qty%template%signal == 0 ) then
          call Announce_Error ( node, no_error_code, &
            & 'The quantity has no signal so prefixSignal is not appropriate' )
          return
        end if
        call GetSignalName ( qty%template%signal, labelStr, &
          & sideband=qty%template%sideband )
        call Get_String( label, labelStr(len_trim(labelStr)+1:), strip=.true. )
        ! Now get an index for this possibly new name which may include the signal
        label = enter_terminal ( trim(labelStr), t_string, caseSensitive=.true. )
      end if

      ! Attach the label
      qty%label = label
    case ( 'vector' )
      do quantityIndex=1, size( vectors(vectorIndex)%quantities )
        qty => GetVectorQuantity( vectors(vectorIndex), quantityIndex )
        labelStr = ' '
        if ( qty%template%name /= 0 ) then
          call get_string( qty%template%name, labelStr, strip=.true. )
        endif
        if ( len_trim(labelStr) == 0 ) then
          write( labelStr, '(a9,i3.3,a1)' ) 'quantity[', quantityIndex, ']'
        endif
        ! call outputNamedValue( 'label string start', labelStr, advance='yes' )
        call Get_String( label, labelStr(len_trim(labelStr)+1:), strip=.true. )
        ! Attach the label
        vlabel = enter_terminal ( trim(labelStr), t_string, caseSensitive=.true. )
        qty%label = vlabel
      enddo
    case default
      call Announce_Error ( node, NO_ERROR_CODE, &
        & 'Must specify quantity or vector in Label' )
    end select
  end subroutine LabelVectorQuantity

  ! --------------------------------------------------  JoinQuantities  -----
  ! This routine parses a line of the l2cf that is designed to join
  ! quantities together into l2gp/l2aux files
  subroutine JoinQuantities ( node, vectors, l2gpDatabase, l2auxDatabase, &
    & chunkNo, chunks )

    use CHUNKS_M, only: MLSCHUNK_T
    ! use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: &
      & F_PRECISION, F_PREFIXSIGNAL, F_SOURCE, F_SDNAME, F_SWATH, FIELD_FIRST, &
      & FIELD_LAST
    use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA, S_L2AUX, S_L2GP
    use INTRINSIC, only: L_NONE! , PHYQ_DIMENSIONLESS
    use L2AUXDATA, only: L2AUXDATA_T
    use L2GPDATA, only: L2GPDATA_T
    use L2PARINFO, only: PARALLEL, SLAVEJOIN
    ! use MLSKINDS, only: R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSSIGNALS_M, only: GETSIGNALNAME
    use MORETREE, only: GET_BOOLEAN, GET_FIELD_ID, GET_SPEC_ID
    use STRING_TABLE, only: GET_STRING
    use SYMBOL_TABLE, only: ENTER_TERMINAL
    use SYMBOL_TYPES, only: T_STRING
    use TREE, only: DECORATION, NSONS, NULL_TREE, SUB_ROSA, SUBTREE
    use VECTORSMODULE, only: GETVECTORQTYBYTEMPLATEINDEX, &
      & VALIDATEVECTORQUANTITY, VECTOR_T, VECTORVALUE_T

    ! Dummy arguments
    integer, intent(in) :: NODE         ! The start of the l2cf line
    type (Vector_T), dimension(:), pointer :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks

    ! Local variables
    ! integer :: EXPRUNITS(2)                 ! From expr
    integer :: FIELDINDEX               ! F_..., see Init_Tables_Module
    ! integer :: FILE                     ! Name of output file for direct write
    integer :: HDFNAMEINDEX             ! Name of swath/sd
    ! integer :: HDFVERSION               ! Version of hdf for directwrite
    integer :: KEY                      ! Index of an L2GP or L2AUX tree
    integer :: KEYNO                    ! Index of subtree of KEY
    integer :: NAME                     ! Sub-rosa index of name of L2GP or L2AUX
    integer :: PRECQTYINDEX             ! Index for precision qty (in database not vector)
    integer :: PRECVECTORINDEX          ! Index for precision vector
    integer :: QUANTITYINDEX            ! ind in qty tmpl database, not vector
    integer :: SON                      ! Son of Key
    integer :: SOURCE                   ! Index in AST
    integer :: VECTORINDEX              ! Index for vector to join
    ! logical :: COMPAREOVERLAPS
    ! logical :: OutputOverlaps
    logical :: PREFIXSIGNAL             ! Prefix (i.e. make) the sd name the signal
    ! real (r8) :: EXPRVALUE(2)               ! From expr

    character(len=132) :: HDFNAME          ! Name for swath/sd
    logical :: GOT(field_first:field_last)
    type (VectorValue_T), pointer :: Quantity
    type (VectorValue_T), pointer :: PrecisionQuantity

    ! We know this node is named
    key = subtree(2,node)
    name = sub_rosa(subtree(1,node))

    got = .false.
    source = null_tree
    ! compareOverlaps = .false.
    ! outputOverlaps = .false.
    hdfNameIndex=name
    prefixSignal = .false.
    ! hdfVersion = 4

    ! Loop over the fields of the mlscf line
    do keyNo = 2, nsons(key) ! Skip spec name
      son = subtree(keyNo,key)
      fieldIndex = get_field_id(son)
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      case ( f_source )
        source = subtree(2,son) ! required to be an n_dot vertex
        vectorIndex = decoration(decoration(subtree(1,source)))
        quantityIndex = decoration(decoration(decoration(subtree(2,source))))
      case ( f_precision )
        source = subtree(2,son) ! required to be an n_dot vertex
        precVectorIndex = decoration(decoration(subtree(1,source)))
        precQtyIndex = decoration(decoration(decoration(subtree(2,source))))
      ! case ( f_hdfVersion )
        ! call expr ( subtree(2,son), exprUnits, exprValue )
        ! if ( exprUnits(1) /= phyq_dimensionless ) &
        !   & call Announce_error ( son, NO_ERROR_CODE, &
        !  & 'No units allowed for hdfVersion: just integer 4 or 5')
        ! hdfVersion = exprValue(1)
      case ( f_prefixSignal )
        prefixSignal = get_boolean(son)
      ! case ( f_compareoverlaps )
        ! compareOverlaps = get_boolean(son)
      ! case ( f_outputoverlaps )
        ! outputOverlaps = get_boolean(son)
      case ( f_swath )
        hdfNameIndex = sub_rosa(subtree(2,son))
      case ( f_sdName )
        hdfNameIndex = sub_rosa(subtree(2,son))
      ! case ( f_file )
      !   file = sub_rosa(subtree(2,son))
      case default ! Can't get here if tree_checker worked properly
      end select
    end do
    
      ! Some final checks
    ! if ( any ( got ( (/ f_file, f_hdfVersion /) ) ) ) &
    !   & call Announce_Error ( key, NO_ERROR_CODE, &
    !  & 'File or hdfVersion not appropriate arguments for output l2aux/l2gp' )

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Errors in configuration prevent proceeding" )

    ! Identify the quantity
    quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex),quantityIndex)
    ! Get the precision quantity too perhaps
    if ( got ( f_precision ) ) then
      precisionQuantity => &
        & GetVectorQtyByTemplateIndex(vectors(precVectorIndex),precQtyIndex)
      if ( quantity%template%name /= precisionQuantity%template%name ) &
        & call announce_error(key, NO_ERROR_CODE, &
        & 'Quantity and precision quantity do not match')
    else
      precisionQuantity => NULL()
    end if
    
    ! Establish a swath/sd name for this quantity.
    hdfName = ''
    if ( prefixSignal ) &
      & call GetSignalName ( quantity%template%signal, hdfName, &
      &   sideband=quantity%template%sideband )
    call Get_String( hdfNameIndex, hdfName(len_trim(hdfName)+1:), strip=.true. )
    ! Now get an index for this possibly new name which may include the signal
    hdfNameIndex = enter_terminal ( trim(hdfName), t_string, caseSensitive=.true. )
    
    ! Now do the join, perhaps as a parallel slave, perhaps more directly.
    if ( parallel%slave ) then
      ! For slave tasks in a PVM system, simply ship this vector off to the master
      call SlaveJoin ( quantity, precisionQuantity, hdfName, key )
    else
      ! Now, depending on the properties of the source we deal with the
      ! vector quantity appropriately.
      if ( ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true., regular=.true., &
        & minorFrame=.false., majorFrame=.false., &
        & verticalCoordinate = (/ l_pressure, l_zeta, l_none/) ) ) then 
        ! Coherent, stacked, regular quantities on pressure surfaces, or
        ! with no vertical coordinate system go in l2gp files.
        if ( get_spec_id(key) /= s_l2gp ) then
          call Announce_Error ( key, NO_ERROR_CODE, &
            & 'This quantity should be joined as an l2gp' )
          call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2gp')
        end if
        call JoinL2GPQuantities ( key, hdfNameIndex, quantity, &
          & precisionQuantity, l2gpDatabase, chunkNo )
      else
        ! All others go in l2aux files.
        if ( get_spec_id(key) /= s_l2aux ) then
          call Announce_Error ( key, NO_ERROR_CODE, &
            & 'This quantity should be joined as an l2aux' )
          call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2aux')
        end if
        call JoinL2AUXQuantities ( key, hdfNameIndex, quantity, &
          & l2auxDatabase, chunks )
      end if
    end if

  end subroutine JoinQuantities

  ! -----------------------------------------  JoinL2GPQuantities  -----

  ! Warning: outmoded--may not work properly
  ! We have not kept this subroutine up-to-date
  ! as we have shifted almost entirely to doing things
  ! by DirectWrite instead of the older Join then Output
  
  ! This routine joins an l2gp line quantity to a database of such quantities.
  ! If this is the first time through, the database is created.

  ! The firstInstance and lastInstance arguments give an optional range of
  ! the instances that we wish to store in the l2gp quantity.  Otherwise, it
  ! defaults to the non overlapped region.

  subroutine JoinL2GPQuantities ( key, name, quantity, &
    & precision, l2gpDatabase, chunkNo, &
    & firstInstance, lastInstance, nameString )

    use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA
    use INTRINSIC, only: L_NONE
    use L2GPDATA, only: ADDL2GPTODATABASE, EXPANDL2GPDATAINPLACE, &
      & L2GPDATA_T, SETUPNEWL2GPRECORD, RGP
    use MLSKINDS, only: RV
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, TOGGLE, LEVELS
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION
    use VECTORSMODULE, only: VECTORVALUE_T

    ! Dummy arguments
    integer, intent(in) :: KEY          ! spec_args to Decorate with the L2GP index
    integer, intent(in) :: NAME         ! For the swath
    type (VectorValue_T), intent(in) :: QUANTITY ! Vector quantity
    type (VectorValue_T), pointer :: precision ! Optional vector quantity
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    integer, intent(in) :: CHUNKNO
    integer, intent(in), optional :: FIRSTINSTANCE, LASTINSTANCE
    character(len=*), intent(in), optional :: nameString
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2gp data.

    ! Local variables
    type (L2GPData_T) :: NewL2GP
    type (L2GPData_T), pointer :: ThisL2GP
    integer :: Index
    integer :: FirstProfile, LastProfile ! Profile range in the l2gp to output to
    real(rv) :: HUGERGP
    integer :: Me = -1                   ! String index for trace
    integer :: NoSurfsInL2GP, NoFreqsInL2GP
    integer :: UseFirstInstance, UseLastInstance, NoOutputInstances
    logical :: L2gpDataIsNew
   
    hugeRgp = real ( huge(0.0_rgp), rv )
    call trace_begin ( me, "JoinL2GPQuantities", key, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    ! If this is the first chunk, we have to setup the l2gp quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.
    l2gpDataIsNew = (.not. associated(l2gpDatabase))
    if ( .not. l2gpDataIsNew ) then
      index = decoration(key)
      l2gpDataIsNew = (index>=0)
    end if

    ! Work out what to do with the first and last instance information
    
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    noOutputInstances = useLastInstance-useFirstInstance+1
    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    if ( l2gpDataIsNew ) then
      ! Now create an empty L2GP record with this dimension

      if (any(quantity%template%verticalCoordinate == (/l_Pressure, l_Zeta /) )) then
        noSurfsInL2GP = quantity%template%noSurfs
      else
        noSurfsInL2GP = 0
      end if

      if ( quantity%template%frequencyCoordinate == l_None) then
         noFreqsInL2GP=0
      else
         noFreqsInL2GP=quantity%template%noChans
      end if

      call SetupNewL2GPRecord ( newL2GP, noFreqsInL2GP, noSurfsInL2GP )
      ! Setup the standard stuff, only pressure as it turns out.
      if ( quantity%template%verticalCoordinate == l_Pressure ) &
        & newL2GP%pressures = quantity%template%surfs(:,1)
      if ( quantity%template%verticalCoordinate == l_Zeta ) &
        & newL2GP%pressures = 10.0**(-quantity%template%surfs(:,1))
      ! It inherits its quantity type from the quantity template
      newL2GP%quantityType=quantity%template%quantityType
      ! Do something about frequency
      if ( associated ( quantity%template%frequencies ) ) then
        newL2GP%frequency = quantity%template%frequencies
      else
        newL2GP%frequency = 0.0
      end if

      ! Add it to the database of l2gp quantities
      index = AddL2GPToDatabase ( l2gpDatabase, newL2GP )
      ! Setup the pointer and index to be used later
      call decorate ( key, -index ) ! Remember where it is
      thisL2GP => l2gpDatabase(index)

    else
      ! Setup the index and pointer
      thisL2GP => l2gpDatabase(-index)
    end if

    ! Expand l2gp (initially all zero-size arrays) to take the new information
    call ExpandL2GPDataInPlace ( thisL2GP, &
      &thisL2GP%nTimes+noOutputInstances )
    thisL2GP%nTimesTotal = quantity%template%grandTotalInstances
    ! Now copy the information from the quantity to the l2gpData

    ! name is an integer, but L2GP%name is Character data
    thisL2GP%nameIndex = name
    if ( present(nameString) ) then
      thisL2GP%name = nameString
    else
      call Get_String( name, thisL2GP%name, strip=.true.)
    end if
    lastProfile=thisL2GP%nTimes
    firstProfile=lastProfile-noOutputInstances+1

    ! Now fill the data, first the geolocation
    thisL2GP%latitude(firstProfile:lastProfile) = &
      & quantity%template%geodLat(1,useFirstInstance:useLastInstance)
    thisL2GP%longitude(firstProfile:lastProfile) = &
      & quantity%template%lon(1,useFirstInstance:useLastInstance)
    thisL2GP%solarTime(firstProfile:lastProfile) = &
      & quantity%template%solarTime(1,useFirstInstance:useLastInstance)
    thisL2GP%solarZenith(firstProfile:lastProfile) = &
      & quantity%template%solarZenith(1,useFirstInstance:useLastInstance)
    thisL2GP%losAngle(firstProfile:lastProfile) = &
      & quantity%template%losAngle(1,useFirstInstance:useLastInstance)
    thisL2GP%geodAngle(firstProfile:lastProfile) = &
      & quantity%template%phi(1,useFirstInstance:useLastInstance)
    thisL2GP%time(firstProfile:lastProfile) = &
      & quantity%template%time(1,useFirstInstance:useLastInstance)
    thisL2GP%chunkNumber(firstProfile:lastProfile)=chunkNo

    ! Now the various data quantities.

    ! For v0.1 we're only thinking about value.  The precision will
    ! come from matrices later in 0.5, and the diagnostics such as status
    ! and quality will come later too (probably 0.5, but maybe 1.0)

    thisL2GP%l2gpValue(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeRgp, min ( hugeRgp, &
      &   quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &  (/max(thisL2GP%nFreqs,1),max(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    if (associated(precision)) then
      thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = &
        & reshape ( max ( -hugeRgp, min ( hugeRgp, &
        &   precision%values(:,useFirstInstance:useLastInstance) ) ), &
        &  (/max(thisL2GP%nFreqs,1),max(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    else
      thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0
    end if
    thisL2GP%status(firstProfile:lastProfile)=0
    thisL2GP%quality(firstProfile:lastProfile)=0.0

    call trace_end ( "JoinL2GPQuantities", &
      & cond=toggle(gen) .and. levels(gen) > 0 )
  end subroutine JoinL2GPQuantities

  ! ----------------------------------------  JoinL2AUXQuantities  -----

  ! This subroutine is like the one above, except that the quantities it joins
  ! are destined to go in L2AUX quantities.

  subroutine JoinL2AUXQuantities ( key, name, quantity, l2auxDatabase, &
   & chunks, firstInstance, lastInstance )

    use CHUNKS_M, only: MLSCHUNK_T
    use INTRINSIC, only: L_GEODANGLE, L_MAF
    use L2AUXDATA, only: ADDL2AUXTODATABASE, RESIZEL2AUXDATA, &
      & L2AUXDATA_T, L2AUXRANK, SETUPNEWL2AUXRECORD
    use MLSKINDS, only: R4, R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    use TOGGLES, only: GEN, TOGGLE, LEVELS, SWITCHES
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION
    use VECTORSMODULE, only: VECTORVALUE_T

    ! Dummy arguments
    integer, intent(in) :: KEY     ! spec_args to decorate with the L2AUX index
    integer, intent(in) :: NAME    ! for the sd
    type (VectorValue_T), intent(in) :: quantity
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    integer, intent(in), optional :: firstInstance, lastInstance
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2aux data.

    ! Local variables
    logical :: DEEBUG
    integer :: FIRSTMAF
    integer :: FIRSTPROFILE
    integer :: DB_INDEX
    integer :: LASTMAF
    integer :: LASTPROFILE
    integer :: MAF
    integer :: Me = -1             ! String index for trace
    integer :: NOMAFS
    integer :: NOOUTPUTINSTANCES
    integer :: USEFIRSTINSTANCE
    integer :: USELASTINSTANCE
    logical :: L2AUXDATAISNEW

    real(r8) :: HUGER4
    type (L2AUXData_T) :: NEWL2AUX
    type (L2AUXData_T), pointer :: THISL2AUX

    ! Executable code

    DEEBUG = ( switchDetail(switches, 'join') > -1 )
    call trace_begin ( me, "JoinL2AUXQuantities", key, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    hugeR4 = real ( huge(0.0_r4), r8 )

    if ( DEEBUG ) then
      call output('Joining vector quantity to L2AUX quantities', advance='yes')
      call output('minor frame? ', advance='no')
      call output(quantity%template%minorFrame, advance='no')
      call output('   major frame? ', advance='no')
      call output(quantity%template%majorFrame, advance='no')
      call output('   template  name ', advance='no')
      if ( quantity%template%name < 1 ) then
        call output('   (unnamed) ', advance='yes')
      else
        call display_string(quantity%template%name, strip=.true., advance='yes' )
      end if
    end if

    ! If this is the first chunk, we have to setup the l2aux quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.
    l2auxDataIsNew = (.not. associated(l2auxDatabase))
    if ( .not. l2auxDataIsNew ) then
      db_index = decoration(key)
      l2auxDataIsNew = (db_index>=0)
    end if

    ! Work out what to do with the first and last Instance information
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if
    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if
    noOutputInstances = useLastInstance - useFirstInstance + 1

    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    if ( DEEBUG ) then
      call output('Joining L2Aux quantity with ', advance='no')
      call output(noOutputInstances, advance='no')
      call output(' instances ', advance='yes')
    end if

    ! Now if this is a new l2aux quantity, we need to setup an l2aux data type
    ! for it.
    if ( l2auxDataIsNew ) then
      ! We need to setup the quantity.  In some cases (minor/major frame)
      ! we can tell how big it is going to be.  Otherwise, create it empty
      if ( any ((/ quantity%template%minorFrame, quantity%template%majorFrame /)) ) then
        firstMAF = minval ( chunks%firstMAFIndex )
        lastMAF = maxval ( chunks%lastMAFIndex )
        noMAFs = lastMAF - firstMAF + 1
        if ( DEEBUG ) then
          call output('  firstMAF ', advance='no')
          call output(firstMAF, advance='no')
          call output('  noMAFs ', advance='no')
          call output(noMAFs, advance='yes')
        end if
      else
        ! Otherwise, we don't know how big it will be (at least in the Join
        ! scenario), so create it empty to begin with.
        firstMAF = 1
        noMAFs = 0
      end if
      ! Create the record accordingly
      call SetupNewL2AUXRecord ( newL2AUX, quantity%template, firstMAF, noMAFs )
      newL2AUX%instrumentModule=quantity%template%instrumentModule
      newL2AUX%quantityType=quantity%template%quantityType

      ! Add this l2aux to the database
      db_index = AddL2AUXToDatabase ( l2auxDatabase, newL2AUX )

      ! Setup the pointer and the index to be used later
      call decorate ( key, -db_index ) ! Remember where it is
      thisL2AUX => l2auxDatabase(db_index)
      thisL2AUX%name = name
    else
      ! Not a new l2aux, so just point ourselves to the old one.
      thisL2AUX => l2auxDatabase(-db_index)
    end if

    ! OK, now thisL2AUX points to an appropriate l2aux to fill, be it newly created,
    ! or old, and be it big enough or not.
    ! First, do we need to expand this?
    if ( .not. any ((/ quantity%template%minorFrame, quantity%template%majorFrame /)) ) then
      ! We need to expand this L2AUX to fit in the latest data.
      call ResizeL2AUXData ( thisL2AUX, noOutputInstances )
    end if

    if ( quantity%template%minorFrame .or. quantity%template%majorFrame ) then
      ! Don't forget instanceOffset is for the first non-overlapped instance (ie MAF)
      ! Also remember the L2AUX data is already indexed from zero! Great!
      if (DEEBUG) call output ( "Doing the special calculation of first/last profile", advance='yes' )
      lastProfile = quantity%template%instanceOffset + quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap - &
        & quantity%template%noInstancesLowerOverlap - 1
    else
      lastProfile = thisL2AUX%dimensions(L2AUXRank)%noValues
    end if
    firstProfile = lastProfile - noOutputInstances + 1

    if ( DEEBUG ) then
      call output('  instance offset ' )
      call output( quantity%template%instanceOffset, advance='yes' )
      call output('  firstProfile ', advance='no')
      call output(firstProfile, advance='no')
      call output('   lastProfile ', advance='no')
      call output(lastProfile, advance='yes')
      call output('  FirstInstance ', advance='no')
      call output(useFirstInstance, advance='no')
      call output('   LastInstance ', advance='no')
      call output(useLastInstance, advance='yes')
      call output('  L2AUX%dimensions ', advance='no')
      call output(trim(thisL2AUX%dim_names), advance='yes')
      call output('  L2AUX%dim_units ', advance='no')
      call output(trim(thisL2AUX%dim_units), advance='yes')
      call output('  L2AUX%value_units ', advance='no')
      call output(trim(thisL2AUX%value_units), advance='yes')
      call output('  L2AUX%dimensions(1)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(1)%noValues, advance='no')
      call output('  L2AUX%dimensions(2)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(2)%noValues, advance='no')
      call output('  L2AUX%dimensions(3)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(3)%noValues, advance='yes')
      call output('shape(l2aux values) ', advance='no')
      call output(shape(thisL2AUX%values), advance='yes')
      if ( any ( thisL2AUX%dimensions(L2AUXRank)%dimensionFamily &
        & == (/ L_GeodAngle, L_MAF /) ) ) then
        call output('   dimensions ', advance='no')
        call output(size(thisL2AUX%dimensions(L2AUXRank)%values), advance='no')
      else
        call output(' (dimensions unassociated)', advance='no')
      end if
      call output('   values 3rd coord ', advance='no')
      call output(size(thisL2AUX%values(1,1,:)), advance='yes')
    end if

    select case (thisL2AUX%dimensions(L2AUXRank)%dimensionFamily)
    case ( L_GeodAngle )
      thisL2AUX%dimensions(L2AUXRank)%values(firstProfile:lastProfile)=&
        & quantity%template%phi(1,useFirstInstance:useLastInstance)
    case ( L_MAF )
      do maf = firstProfile, lastProfile
        thisL2AUX%dimensions(L2AUXRank)%values(maf) = maf
      end do
    case default
    end select
    
    ! Check that reshape has a prayer of succeeding
    if ( DEEBUG ) then
      call output('  num l2aux values/profile ', advance='no')
      call output(size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2), &
       & advance='yes')
      call output('   num dim values/profile ', advance='no')
      call output(thisL2AUX%dimensions(1)%noValues* &
       &              thisL2AUX%dimensions(2)%noValues, advance='yes')
      call output('  num l2aux values (total) ', advance='no')
      call output(size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2)* &
       &              (lastProfile-firstProfile+1), advance='yes')
      call output('   num qty values ', advance='no')
      call output(size(quantity%values, 1)* &
       &              (useLastInstance-useFirstInstance+1), advance='yes')
    end if
    if ( size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2) &
     & /= thisL2AUX%dimensions(1)%noValues*thisL2AUX%dimensions(2)%noValues ) &
     & call MLSMessage ( MLSMSG_Error, &
     & ModuleName, "Reshape fails: size mismatch betw. dims and values" )
    if ( size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2)*(lastProfile-firstProfile+1) &
     & /= size(quantity%values, 1)*(useLastInstance-useFirstInstance+1) ) &
     & call MLSMessage ( MLSMSG_Error, &
     & ModuleName, "Reshape fails: size mismatch betw. quantity and l2aux" )
    thisL2AUX%values(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeR4, min ( hugeR4, &
      & quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &   (/ thisL2AUX%dimensions(1)%noValues, &
      &      thisL2AUX%dimensions(2)%noValues, &
      &      lastProfile-firstProfile+1/) )
    
    call trace_end ( "JoinL2AUXQuantities", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  end subroutine JoinL2AUXQuantities

! =====     Private Procedures     =====================================

  function explicitFile ( node ) result ( gotTheFile )
    ! Loop overs fields in DirectWrite command to discover whether
    ! the file is explicitly named
    use INIT_TABLES_MODULE, only: F_FILE
    use MORETREE, only: GET_FIELD_ID
    use TREE, only: NSONS, SUBTREE
    ! Args
    integer, intent(in) :: node
    logical :: gotTheFile
    ! Private variables
    integer :: keyNo
    integer :: fieldIndex
    integer :: son
    ! Executable
    gotTheFile = .false.
    do keyNo = 2, nsons(node)           ! Skip DirectWrite command
      son = subtree ( keyNo, node )
      fieldIndex = get_field_id ( son )
      select case ( fieldIndex )
      case ( f_file )
        gotTheFile = .true.
      end select
    end do
  end function explicitFile

  ! ---------------------------------------------  Announce_Error  -----
  subroutine ANNOUNCE_ERROR ( where, CODE, ExtraMessage, FIELDINDEX, Penalty )

    use INTRINSIC, only: FIELD_INDICES
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    use TREE, only: Where_At => Where

    integer, intent(in) :: where   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: FIELDINDEX ! Extra information for msg
    character (LEN=*), intent(in), optional :: ExtraMessage

    integer, intent(in), optional :: Penalty
    integer :: myPenalty

    myPenalty = 1
    if ( present(penalty) ) myPenalty = penalty
    error = max(error,myPenalty)

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( where_at(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    select case ( code )
      case ( NotAllowed )
        call output('Field ')
        call display_string(field_indices(fieldIndex))
        call output(' is not allowed in this context',advance='yes')
      case default
        call output ( " command caused an unrecognized programming error", advance='yes' )
    end select
    if ( present(ExtraMessage) ) then
      call output(ExtraMessage, advance='yes')
    end if
  end subroutine ANNOUNCE_ERROR

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Join

!
! $Log$
! Revision 2.162  2015/03/28 02:47:14  vsnyder
! Paul added Quantity type to Direct write.
!
! Revision 2.161  2014/04/07 18:03:28  pwagner
! May specify AscDescMode when DirectWrite-ing swaths
!
! Revision 2.160  2014/03/31 23:43:49  pwagner
! Commented-out unused stuff; renamed procedure ResizeL2AUXData, generalizing it to expand or contract
!
! Revision 2.159  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.158  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.157  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.156  2013/11/01 00:15:15  pwagner
! Match trace_begins and _ends scrupulously
!
! Revision 2.155  2013/10/09 23:41:36  vsnyder
! Add Evaluate_Variable
!
! Revision 2.154  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.153  2013/09/04 17:35:23  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.152  2013/08/30 02:45:42  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.151  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.150  2012/08/16 17:58:24  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.149  2012/07/02 20:40:43  pwagner
! Once-routine output now requires verbosity > 0
!
! Revision 2.148  2012/05/08 17:51:11  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.147  2012/05/01 23:16:35  pwagner
! Partly retreated from last optimistic commit
!
! Revision 2.146  2012/03/12 17:22:20  pwagner
! Believe we can drop odious workaround for hdfeos 'bug'--may have been our mistake
!
! Revision 2.145  2012/02/24 21:19:44  pwagner
! May DirectWrite a /single instance only
!
! Revision 2.144  2011/11/04 00:09:53  pwagner
! Fixed bug that prevented writing matched output types
!
! Revision 2.143  2011/10/07 00:06:02  pwagner
! May dump Matrices, Hessians from Fill, Join
!
! Revision 2.142  2011/05/09 18:18:45  pwagner
! Consistent with new api for DirectWrite
!
! Revision 2.141  2009/10/26 17:11:53  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.140  2009/09/29 23:40:45  pwagner
! Unjam error message when DirectWrite type is unexpected
!
! Revision 2.139  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.138  2009/06/02 17:53:15  cvuu
! Add NRT Lat and Lon bounding to metadata
!
! Revision 2.137  2009/04/23 23:02:25  pwagner
! May specify upperOverlap or lowerOverlap in DirectWrites
!
! Revision 2.136  2009/04/01 23:35:59  pwagner
! First attempt at fixing bug that wrote missing values to Temperature-APriori
!
! Revision 2.135  2008/12/18 21:14:24  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.134  2008/12/02 23:27:09  pwagner
! May automatically label every quantity in a vector now
!
! Revision 2.133  2007/12/07 01:50:52  pwagner
! Removed unused dummy variables, etc.
!
! Revision 2.132  2007/11/05 18:37:19  pwagner
! May Skip remaining lines in Fill, Join, Retrieve sections depending on Boolean
!
! Revision 2.131  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.130  2006/10/11 22:58:00  pwagner
! Will write convergence ratio as another quality-like field of l2gp
!
! Revision 2.129  2006/09/21 18:55:04  pwagner
! Fixed bug causing freezes when band13 off; cosmetic changes too
!
! Revision 2.128  2006/06/21 22:06:56  pwagner
! Downgraded inappropriate quantity output to warning mesg
!
! Revision 2.127  2006/04/11 23:34:27  pwagner
! Fixed bug which added excess profiles
!
! Revision 2.126  2006/03/04 00:20:13  pwagner
! Account for directdatabase even if skipping directWrites
!
! Revision 2.125  2006/02/10 21:17:33  pwagner
! dumps may go to special dumpfile
!
! Revision 2.124  2005/11/11 21:47:08  pwagner
! May have fixed bugs Besetting Alyn and Dong
!
! Revision 2.123  2005/11/04 18:54:04  pwagner
! Accommodates changed interface to add_metadata
!
! Revision 2.122  2005/10/28 23:17:37  pwagner
! Print missing outputType when unable to find suitable for distributing sources
!
! Revision 2.121  2005/10/18 23:11:36  pwagner
! Fixed bug causing standalone, toolkitless run to put split dgg/dgm files in cwd
!
! Revision 2.120  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.119  2005/06/16 18:43:01  pwagner
! Should not bomb if catenating split files w/o toolkit
!
! Revision 2.118  2005/06/14 20:43:19  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.117  2005/03/24 21:31:55  pwagner
! May warn of unreasonable directWrite waiting, writing times
!
! Revision 2.116  2004/12/14 22:51:35  pwagner
! Changes related to stopping early
!
! Revision 2.115  2004/07/22 20:49:58  cvuu
! Add forwardModelConfigDatabase to the call MLSL2Join and MLSL2Fill
!
! Revision 2.114  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.113  2004/05/19 20:16:29  vsnyder
! Remove declarations and uses for unreferenced symbols, polish some cannonballs
!
! Revision 2.112  2004/05/19 19:16:10  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.111  2004/04/27 23:56:17  pwagner
! SKIPDIRECTWRITES acts much like checkpaaths
!
! Revision 2.110  2004/04/24 00:21:33  pwagner
! Bombs less readily when not a slave
!
! Revision 2.109  2004/04/17 07:00:54  livesey
! Bug fix possibly?
!
! Revision 2.108  2004/03/03 19:28:04  pwagner
! Should not add metadata if dgg/dgm to be unsplit
!
! Revision 2.107  2004/02/19 23:56:23  pwagner
! Fixed tiny memory leak; skipdgg/m skips directwrite of filetype
!
! Revision 2.106  2004/02/11 23:11:38  livesey
! Logic wrong in calls to DoHGridsMatch
!
! Revision 2.105  2004/02/11 17:21:51  pwagner
! May DirectWrite l2gp status and quality quantities
!
! Revision 2.104  2004/02/10 19:28:36  pwagner
! Prevents named DirectWrites from being looped over
!
! Revision 2.103  2004/02/05 23:40:35  pwagner
! More bugs fixed in automatic directwrites
!
! Revision 2.102  2004/01/23 01:09:48  pwagner
! Only directwrite files entered in global settings eligible to be auto-sourced
!
! Revision 2.101  2004/01/22 00:56:35  pwagner
! Fixed many bugs in auto-distribution of DirectWrites
!
! Revision 2.100  2004/01/02 23:36:00  pwagner
! DirectWrites may choose files automatically from db
!
! Revision 2.99  2003/12/05 00:41:14  pwagner
! patch option avoids deleting existing data in file
!
! Revision 2.98  2003/12/03 17:50:54  pwagner
! L2GP tracks both nTimes (for this slave) and nTimesTotal (done by all)
!
! Revision 2.97  2003/11/14 23:38:45  pwagner
! Uses DirectWrite databse in preference to repeated calls to toolkit
!
! Revision 2.96  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.95  2003/10/20 18:21:45  pwagner
! Timings breakdown added for directWrite
!
! Revision 2.94  2003/10/10 00:00:24  pwagner
! Should quit properly if SIPS and no PCFid match for file name
!
! Revision 2.93  2003/09/12 21:45:52  pwagner
! Only prints l2gp label during DirectWrite if DEEBUG
!
! Revision 2.92  2003/09/04 22:40:04  pwagner
! Gets dgg file name from dgg PCFid when DirectWrite usingPCF
!
! Revision 2.91  2003/09/03 23:05:49  livesey
! More problems with precision in DirectWrite, hope I've got them all now.
!
! Revision 2.90  2003/09/03 00:53:50  livesey
! Bug fix on directWrite, it was storing the value field in the precision
! too!
!
! Revision 2.89  2003/08/28 23:51:58  livesey
! Renamed some variables to make them more obvious
!
! Revision 2.88  2003/08/14 20:11:30  pwagner
! DirectWrite may take l2fwm types for fwm radiances
!
! Revision 2.87  2003/08/01 20:38:31  pwagner
! Distinguishes between l2dgg and l2gp when writing metadata as part of directwrite
!
! Revision 2.86  2003/07/25 00:51:06  livesey
! Added file type l2dgg to support metadata.
!
! Revision 2.85  2003/07/23 18:30:35  cvuu
! reduce routine printing
!
! Revision 2.84  2003/07/15 23:39:01  pwagner
! Disabled most printing
!
! Revision 2.83  2003/07/11 01:24:20  livesey
! More changes trying to get the direct write going.
!
! Revision 2.82  2003/07/09 21:49:53  pwagner
! Tries to figure out in advance whether to create swath
!
! Revision 2.81  2003/07/08 00:15:51  livesey
! Various tidy ups and reworks
!
! Revision 2.80  2003/07/07 23:52:13  pwagner
! Slave that creates DirectWrite file may also add_metadata
!
! Revision 2.79  2003/07/07 20:29:43  livesey
! Mainly cosmetic changes
!
! Revision 2.78  2003/07/07 17:31:11  livesey
! Various things to get DirectWrite working
!
! Revision 2.77  2003/07/02 00:55:27  pwagner
! Some improvements in DirectWrites of l2aux, l2gp
!
! Revision 2.76  2003/06/26 23:13:52  pwagner
! New debugging output; distinguishes between l2gp/l2aux quantities better
!
! Revision 2.75  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.74  2003/06/24 23:30:00  livesey
! Finished LabelVectorQuantity and made some other bug fixes.
!
! Revision 2.73  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.72  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.71  2003/05/30 00:09:27  livesey
! Can now directWrite major frame quantities
!
! Revision 2.70  2003/05/12 02:06:23  livesey
! Bug fix for prefixSignal L2GPs and also bound r8->r4 conversion
!
! Revision 2.69  2003/02/08 00:31:31  pwagner
! Now saves quantityType in newl2gp
!
! Revision 2.68  2003/01/30 01:03:24  pwagner
! Stores quantity type taken from source vector in l2aux
!
! Revision 2.67  2003/01/17 23:11:26  pwagner
! Moved most ops out of LoinL2AUXData to SetupL2AUXData
!
! Revision 2.66  2002/12/19 15:53:47  livesey
! Allowed verticalCoordinate=l_none quantities back into the l2gp fold.
!
! Revision 2.65  2002/11/26 23:38:01  livesey
! Better joining of major frame quantities
!
! Revision 2.64  2002/10/29 21:54:21  livesey
! Made join less verbose.
!
! Revision 2.63  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.62  2002/08/20 22:10:50  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.61  2002/08/20 20:10:30  livesey
! Dealt with frequency in l2gps
!
! Revision 2.60  2002/05/28 17:09:57  livesey
! Removed print statements
!
! Revision 2.59  2002/05/22 00:48:52  livesey
! Added direct write stuff
!
! Revision 2.58  2002/05/16 22:36:46  livesey
! Fixed a bug with joining minor frame quantities with overlaps.
!
! Revision 2.57  2002/04/08 20:49:17  pwagner
! Swath name optionally passed to JoinL2GPQuantities
!
! Revision 2.56  2002/04/06 00:35:21  pwagner
! Should accept actual case of swathname for l2gp
!
! Revision 2.55  2002/03/20 00:46:47  pwagner
! Removed 2 unused lits
!
! Revision 2.54  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.53  2001/10/30 01:45:21  livesey
! Some modifications/fixes to parallel join
!
! Revision 2.52  2001/10/08 23:38:58  pwagner
! Tiny fixes; not perfect yet
!
! Revision 2.51  2001/10/06 00:27:42  pwagner
! Still some problems with diagnostics
!
! Revision 2.50  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.49  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.48  2001/09/08 00:21:44  pwagner
! Revised to work for new column Abundance in lone swaths
!
! Revision 2.47  2001/09/05 20:34:56  pwagner
! Reverted to pre-columnAbundance state
!
! Revision 2.46  2001/08/03 23:13:52  pwagner
! Began testing; at least now exits normally again
!
! Revision 2.45  2001/08/02 23:58:31  pwagner
! More complete treatment of column abundance(s)
!
! Revision 2.44  2001/08/02 00:18:55  pwagner
! Began adding column quantities; incomplete
!
! Revision 2.43  2001/07/31 23:25:32  pwagner
! Able to accept 2 new fields for join of column; does nothing yet
!
! Revision 2.42  2001/06/19 22:52:31  pwagner
! l_none  no longer got from init_tables_module
!
! Revision 2.41  2001/05/23 21:59:43  livesey
! Interim version, almost there
!
! Revision 2.40  2001/05/23 01:43:19  livesey
! New parallel version in progress
!
! Revision 2.39  2001/05/19 01:19:58  livesey
! Fills precision correctly for l2gp!
!
! Revision 2.38  2001/05/14 22:23:53  livesey
! Embarassing bug fix, to do with renumbering of minor frame L2AUX quantities.
!
! Revision 2.37  2001/05/12 00:18:17  livesey
! Tidied up array bounds for L2AUX/MAF.
!
! Revision 2.36  2001/05/10 16:31:24  livesey
! Added prefix signal option for swath/sd name
!
! Revision 2.35  2001/05/08 23:25:32  livesey
! Added the precision stuff for l2gp's
!
! Revision 2.34  2001/05/08 21:51:02  livesey
! Removed some old xStar, yStar, kStar stuff.
!
! Revision 2.33  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.32  2001/05/02 22:22:43  pwagner
! Removed SDPToolkit use
!
! Revision 2.31  2001/04/28 01:30:14  livesey
! Basically gone back to an earlier version.  As l2pc's now output
! directly as matrices there is no need for Join to think about them.
!
! Revision 2.30  2001/04/27 21:52:39  livesey
! Removed l2pc stuff
!
! Revision 2.29  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.28  2001/04/26 15:59:30  livesey
! Tidied up uses
!
! Revision 2.27  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.26  2001/04/26 00:07:16  livesey
! Insulate vector is gone
!
! Revision 2.25  2001/04/25 21:54:22  livesey
! Added candol2pc flag
!
! Revision 2.24  2001/04/25 21:51:46  livesey
! Tidied up Join for l2pcs
!
! Revision 2.23  2001/04/25 20:33:28  livesey
! Minor improvements to Join l2pc stuff
!
! Revision 2.22  2001/04/24 20:20:27  livesey
! L2PC moved to lib and word bin dropped from types etc.
!
! Revision 2.21  2001/04/24 20:04:54  livesey
! Added l2pc joining
!
! Revision 2.20  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.19  2001/03/15 23:26:56  livesey
! Avoid calling ExpandL2AUXQuantitiesInPlace for minor frame quantities.
! Really saves on memory thrashing.
!
! Revision 2.18  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.17  2001/03/06 22:40:41  livesey
! New L2AUX stuff
!
! Revision 2.16  2001/03/05 20:46:41  livesey
! Removed a debugging statement left behind
!
! Revision 2.15  2001/03/05 01:19:45  livesey
! Removed a print statement
!
! Revision 2.14  2001/03/05 01:01:12  livesey
! Bug fix, now uses GetVectorQtyFromTemplateIndex
!
! Revision 2.13  2001/03/01 18:38:27  livesey
! Fixed bug with verticalCoordinate==l_Zeta
!
! Revision 2.12  2001/02/27 17:38:21  livesey
! Tidied things up, removed unnecessary arguments
!
! Revision 2.11  2001/02/27 00:50:31  livesey
! Added ability to Join verticalCoordinate=l_zeta quantities into l2gp entities.
!
! Revision 2.10  2001/02/16 00:50:17  livesey
! Added error to avoid confusion with L2GP in ReadApriori section
!
! Revision 2.9  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.8  2001/02/09 18:01:46  livesey
! Various further updates, set default values for status and quality
!
! Revision 2.7  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.6  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.5  2000/11/16 02:19:01  vsnyder
! Implement timing.
!
! Revision 2.4  2000/11/13 23:02:21  pwagner
! Adapted for rank2 vectorsModule
!
! Revision 2.3  2000/10/05 16:37:19  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:34:35  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

