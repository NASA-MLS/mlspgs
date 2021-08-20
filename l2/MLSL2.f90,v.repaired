! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program MLSL2
  use Allocate_DeAllocate, only: AllocateLogUnit, FinalMemoryReport, &
    & Set_Garbage_Collection, TrackAllocates
  use Call_Stack_M, only: Show_Final_Summary, Show_Sys_Memory, &
    & Sys_Memory_Ch, Sys_Memory_Convert
  use ChunkDivide_M, only: ChunkDivideConfig
  use Declaration_Table, only: Allocate_Decl, DeAllocate_Decl, Dump_Decl
  use Dump_1, only: Dump
  use EmpiricalGeometry, only: DestroyEmpiricalGeometry
  use HDF, only: Dfacc_Rdonly
  use HighOutput, only: Dump, HeadLine, OutputNamedValue
  use Init_Tables_Module, only: Init_Tables
  use Io_Stuff, only: Write_TextFile
  use Intrinsic, only: Get_Type, L_Ascii, L_Tkgen, Lit_Indices
  use L2GPData, only: AvoidUnlimitedDims
  use L2ParInfo, only: Parallel, InitParallel, AccumulateSlaveArguments, &
    & TransmitSlaveArguments
  use LeakCheck_M, only: LeakCheck
  use Lexer_Core, only: Init_Lexer
  use Machine, only: Getarg, Hp, Io_Error, USleep
  use MLSCommon, only: MLSFile_T, MLSNamesAreDebug, MLSNamesAreVerbose
  use MLSFiles, only: FileStringTable, &
    & AddFileToDatabase, DeAllocate_FileDatabase, Dump, &
    & InitializeMLSFile, MLS_OpenFile, MLS_CloseFile
  use MLSHDF5, only: MLS_H5Open, MLS_H5Close
  use MLSL2Options, only: AllocFile, Aura_L1BFiles, &
    & CheckL2CF, CheckLeak, CheckPaths, CountChunks, Current_Version_Id, &
    & Default_HDFVersion_Read, Default_HDFVersion_Write, Do_Dump, Dump_Tree, &
    & DumpOptions, L2cf_Unit, Level1_HDFVersion, MaxChunkSize, MLSL2Message, &
    & Need_L1BFiles, Normal_Exit_Status, NoteFile, NumSwitches, &
    & L2Options, OriginalOptions, &
    & Patch, PhasesToSkip, ProcessOptions, Quit_Error_Threshold, &
    & Recl, RestartWarnings, RunTimeValues, &
    & SectionsToSkip, SectionTimes, SectionTimingUnits, &
    & SharedPCF, ShowDefaults, & ! Sips_Version, &
    & SkipDirectWrites, SkipDirectWritesOriginal, &
    & SlavesCleanUpSelves, SlaveMAF, &
    & SpecialDumpFile, StateFilledBySkippedRetrievals, &
    & StopAfterSection, StopWithError, &
    & Timing, Toolkit, TotalTimes, UniqueID
  use MLSL2Timings, only: Run_Start_Time, Section_Times, Total_Times, &
    & Add_To_Section_Timing, Dump_Section_Timings
  use MLSMessageModule, only: MLSMSG_Debug, &
    & MLSMessageConfig, MLSMSG_Error, MLSMSG_Severity_To_Quit, &
    & MLSMSG_Success, MLSMSG_Warning, DumpConfig, MLSMessage, MLSMessageExit
  use MLSPCF2 ! Everything
  use MLSStrings, only: Trim_Safe
  use MLSStringLists, only: ExpandStringRange, PutHashElement, SwitchDetail
  use Output_M, only: Blanks, FlushStdout, Output, &
    & BothPrUnit, InvalidPrUnit, MSGLogPrUnit, OutputOptions, PrUnitName, &
    & SetOutputStatus, StampOptions, StdoutPrUnit, SwitchOutput, RevertOutput
  use Parser, only: Clean_Up_Parser, Configuration
  use Parser_Table_M, only: Destroy_Parser_Table, Parser_Table_T
  use Parser_Tables_L2cf, only: Init_Parser_Table
  ! use PCFHdr, only: DumpGlobalAttributes
  use Printit_M, only: Set_Config, StdoutLogUnit
  use PVM, only: ClearPVMArgs, FreePVMArgs
  use SDPToolkit, only: PGSD_DEM_30arc, PGSD_DEM_90arc, &
    & PGSD_DEM_Elev, PGSD_DEM_Water_Land, UseSDPToolkit, PGS_DEM_Close
  use String_Table, only: Destroy_Char_Table, Destroy_Hash_Table, &
    & Destroy_String_Table, Get_String, AddinUnit
  use Symbol_Table, only: Destroy_Symbol_Table
  use Time_M, only: SayTime_Config, Time_Config, &
    & Begin, Dump, Finish, SayTime, Time_Now
  use Toggles, only: Levels, Syn, Switches, Toggle
  use Track_M, only: ReportLeaks
  use Tree, only: Allocate_Tree, DeAllocate_Tree, NSons, SubTree
  use Tree_Checker, only: Check_Tree
  use Tree_Walker, only: Walk_Tree_To_Do_MLS_L2

  ! === (start of toc) ===
  !     c o n t e n t s
  !     - - - - - - - -

  ! Main program for level 2 processing
  ! === (end of toc) ===
  ! (It is assumed that mlsl1 has already been run successfully)
  ! Usage:
  ! mlsl2 [options] [<] [l2cf]
  ! where l2cf is an ascii file which comes from one of
  ! (i)   a file named by a line in the pcf
  !        (if and only if toolkit is TRUE)
  ! (ii)  a file named on the command line w/o the '<' redirection
  ! (iii) stdin or a file redirected as stdin using '<'
  ! and where we expand 'mlsl2 [options]' below
  ! mlsl2 -slo1 -slo2 .. --mlo1 --mlo2 ..
  !    where slon are single-letter options and mlon are multiletter options
  !    E.g., mlsl2 -m -p --nmeta -Sglo,jac
  !    For a list of available options enter 'mlsl2 --help'
  !    For a list of available switches enter 'mlsl2 -S"?"'
  ! In case the l2cf is named, say, "file_name" by (ii), we will try:
  ! First, to find file_name in the current directory
  ! If that fails, we will try file_name.l2cf in the same directory
  ! If that fails, woe is us and we exit with an error

  ! Overview
  ! Level 2 must accomplish some operational tasks (opening and reading files),
  ! prior to any scientific tasks (processing data)
  ! The operational tasks are under the control of any and all of
  ! (a) hard-wired options (see MLSL2Options module)
  ! (b) command-line options (see /mlspgs/notes/options and
  !         /mlspgs/notes/switches)
  ! (c) the l2cf
  ! (d) the pcf (unless the variable 'pcf' is FALSE)
  ! Tasks
  ! (1) Accept hard-wired options (see module MLSL2OPtions.f90)
  ! (2) Initialize parser structures
  ! (3) Read command-line options
  ! (4) Open the l2cf (if not stdin)
  ! (5) Parse the l2cf
  ! (6) Close the l2cf (if not stdin)
  ! (7) Call walk_tree_to_do_MLS_L2 to do the (mostly) scientific tasks
  ! (8) Close down some remaining things and exit with an appropriate status

  ! Parallel processing
  ! In parallel with this, there is a possible cooperation where
  ! a single master process starts up, then tries to pass messages
  ! to slave processes which may be running. Master and slaves, all,
  ! will be running mlsl2
  ! The interprocess communication is handled by pvm
  !
  ! Three alternatives are available to manage these tasks
  ! (1) "submit", where you specify --submit "command" as
  !       a command-line option, and the master uses "command"
  !       with a generic, non-mlsl2-aware batch queue system;
  !       (this has been successfully tested via --submit mlssubmit)
  ! (2) l2q, where you specify --submit l2q as the command-line option;
  !      l2q is a special-purpose mlsl2-aware executable
  !      (Different from (1), but option line looks like (1))
  ! (3) Direct control, where the master has sole responsibilty for
  !      inquiring from pvm as to available hosts, and using them
  !
  ! Notes on parallel processing:
  ! (a) (1) and (2) may permit more than one master task to run
  !     simultaneously. They require a queue manager to be running already.
  !     In contrast, (3) requires only the pvm demon, but you are limited
  !     to one master at a time
  ! (b) The l2cf must be specified; you can't use stdin (the slaves wouldn't
  !     see the master's stdin)
  implicit none
  ! ********************************************************
  ! Should we move these next 4 to MLSL2Options?
  character(len=*), parameter :: L2CFNameExtension = ".l2cf"
  logical, parameter          :: WaitForscript     = .true.
  logical, parameter          :: Wrap              = .false.
  integer, parameter          :: WrapPastColumn    = 0 ! 140 ! 192
  ! ********************************************************

  integer :: Error                 ! Error flag from check_tree
  integer :: First_Section         ! Index of son of root of first n_cf node
  logical :: garbage_collection_by_dt = .false. ! Collect garbage after each deallocate_test?
  integer :: I                     ! counter for command line arguments
  integer, dimension(1) :: IChunks
  integer :: J                     ! index within option
  character(len=2048) :: Line      ! Into which is read the command args
  character(len=1) :: null
  integer :: Numfiles
  type(Parser_Table_t) :: Parser_Table
  integer :: Root                  ! of the abstract syntax tree
  integer :: Status                ! From OPEN
  real :: T0, T1, T2               ! For timing
  integer :: inunit = -1

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  type (MLSFile_T), dimension(:), pointer ::     FileDatabase
  type (MLSFile_T)                        ::     MLSL2CF
  ! DEM nonsense
  integer, dimension(2) :: resolutionList
  integer   :: numResolutions
  integer, dimension(2) :: layerList
  integer   :: numLayers
  ! Executable
  nullify (FileDatabase)
  null = achar(0)
  !---------------- Task (1) ------------------

  call time_now ( t0 )
  call MLS_H5Open ( error )
  if (error /= 0) then
      call MLSL2Message ( MLSMSG_Error, moduleName, &
        & "Unable to MLS_H5Open" )
      if ( not_used_here() ) print *, "This never gets executed"
  end if
! Initialize the lexer, symbol table, and tree checker's tables:
!  ( Under some circumstances, you may need to increase these )
  call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957, &
    & DEBUG=0 )
  call allocate_decl ( ndecls=8000 )
  call allocate_tree ( n_tree=2000000 )
! call init_tables ! Postpone this until after toggles get set
  FileStringTable = .true.
  !---------------- Task (2) ------------------
! Where to send output, how severe an error to quit
! call SwitchOutput ( 'stdout' )
! call Dump( OutputOptions )
! call DumpConfig
! call RevertOutput
  outputOptions%prunit = L2Options%Output_print_unit
  MLSMSG_Severity_to_quit = MAX(Quit_Error_Threshold, MLSMSG_Debug+1)
  call set_config ( severity_to_quit = MLSMSG_Severity_to_quit )

! Clear the command line arguments we're going to accumulate to pass
! to slave tasks
  call ClearPVMArgs

  ! We set up a mirror command line for launching slaves
  call getarg ( hp, parallel%executable )

  !---------------- Task (3) ------------------
  !---------------- Task (3a) ------------------
  ! Command line options
  i = 1+hp
  L2Options%Originalcmds = ' '
  do ! Process Lahey/Fujitsu run-time options; they begin with "-Wl,"
    call getarg ( i, line )
    if ( line(1:4) /= '-Wl,' ) exit
    call AccumulateSlaveArguments(line) ! pass them to slave processes
    i = i + 1
  end do
  do ! Process the command line options to set toggles
    call getarg ( i, line )
    ! print *, trim(line)
    if ( len_trim(line) < 1 ) exit
    i = i + 1
    L2Options%Originalcmds = trim(L2Options%Originalcmds) // null // trim(line)
  enddo
  line = processOptions( trim(L2Options%Originalcmds), null )
  OriginalOptions = L2Options
! call SwitchOutput ( 'stdout' )
! call Dump( OutputOptions )
! call DumpConfig
! call RevertOutput
! stop
  if ( line == 'help' ) then
    call option_usage
  elseif( switchDetail(switches, '?') > -1 .or. &
    & switchDetail(switches, 'help') > -1 ) then
    call switch_usage
  end if
  MLSMSG_Severity_to_quit = MAX(Quit_Error_Threshold, MLSMSG_Debug+1)
  ! print *, 'slaveArguments: ', trim(slaveArguments)
  ! The following no longer does anything
  call Set_garbage_collection(garbage_collection_by_dt)

! Done with command-line parameters; enforce cascading negative options
! (waited til here in case any were (re)set on command line)
! print *, 'Enforcing output destinations'! 
! print *, 'tookit: ', toolkit
! print *, 'showDefaults: ', showDefaults
! print *, 'switches: ', trim(switches)

  !---------------- Task (3b) ------------------
  ! Redirecting output and optional logging by the Toolkit
  if ( L2Options%Overridden .and. .not. parallel%slave ) then
    ! Getting these things from L2Options
     call MLSMessage( MLSMSG_Warning, ModuleName, &
      & 'Setting output destination based on L2Options' )
  elseif ( .not. toolkit .or. showDefaults ) then
     outputOptions%prunit = max( StdoutPrUnit, &
       & outputOptions%prunit)   ! stdout or Fortran unit
  else if ( outputOptions%prunit == INVALIDPRUNIT ) then
     call MLSMessage( MLSMSG_Warning, ModuleName, &
      & 'Avoding all output except possibly MLSMessages' )
  else if ( parallel%master ) then
     outputOptions%prunit = BothPrUnit ! MSGLOGPRUNIT   ! output both logged, not sent to stdout
     call set_config( LogFileUnit=BothPrUnit )
  else if ( parallel%slave ) then
     outputOptions%prunit = StdoutPrUnit   ! output sent only to stdout, not logged
     call set_config( useToolkit=.false. )
  end if
  i = SwitchDetail(switches, 'log')
! print *, 'i: ', i
  if( .not. L2Options%Overridden .and. (i == 0 .or. i > 5 .or. .not. toolkit) ) then
     call set_config ( LogFileUnit = StdoutLogUnit ) ! -1
  end if
  if ( i > 9 ) then
    MLSMessageConfig%MaxModuleNameLength   = i - 10
    MLSMessageConfig%MaxSeverityNameLength = i - 10
  endif

  ! Dont log parent module names unless asked by -Sparent
  outputOptions%logparent = ( SwitchDetail(switches, 'parent') > -1 ) 
  ! How to prefix logged messages?
  ! But don't stomp on effect of laconic
  i = SwitchDetail(switches, 'mess')
  select case (i)
  case (0)
    ! Skip writing module names, severity to log file for routine output
    MLSMessageConfig%skipModuleNamesThr = MLSMSG_Warning
    MLSMessageConfig%skipSeverityThr    = MLSMSG_Warning
  case (1)
    ! Skip writing module names to log file for routine output
    MLSMessageConfig%skipModuleNamesThr = MLSMSG_Warning
  case (2)
    ! Skip writing severity to log file for routine output
    MLSMessageConfig%skipSeverityThr    = MLSMSG_Warning
  case default
    ! Write both module names, severity to log file for routine output
    ! unless countermanded by --laconic commandline option
    if ( MLSMessageConfig%skipModuleNamesThr < MLSMSG_Warning) then
      MLSMessageConfig%skipModuleNamesThr = MLSMSG_Success
      MLSMessageConfig%skipSeverityThr    = MLSMSG_Success
    endif
  end select
  
  ! Report on memory usage during traces?
  Show_Final_Summary = ( SwitchDetail(switches, 'alloc') > -1 )
  section_times = sectiontimes
  total_times = totaltimes

  UseSDPToolkit = toolkit    ! Redundant, but may be needed in lib
  ! Try to prevent calls to output being lost 
  ! between output_m and MLSMessage modules
  if ( MLSMessageConfig%LogFileUnit == StdoutLogUnit ) &
    & outputOptions%prunit = StdoutPrUnit

  if ( time_config%use_wall_clock ) call time_now( run_start_time )
  ! If checking paths, run as a single-chunk case in serial mode
  if ( checkPaths ) then
    parallel%master = .false.
    parallel%slave = .false.
    parallel%chunkRange = '1'
    stopWithError = .false.
    ! Nio longer issue warning about l2pc files
    if ( .false. ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
    & 'checkPaths will fail if l2pc files only on local disks but master runs' &
    & // ' on front end' )
  end if
  ! If doing a range of chunks, the avoidance of unlimited dimensions
  ! in directwrites of l2gp files currently fails
  ! (when will this be fixed? By whom? HDF people or us?)
  if ( parallel%chunkRange /= ' ' ) then
    ! avoidUnlimitedDims = .false.
    call ExpandStringRange( parallel%chunkRange, iChunks )
  else
    iChunks = 0
  endif

  outputOptions%parentName = 'MLSL2'
  if ( Wrap ) then
    call SetOutputStatus ( 'wrap', 1 )
    call SetOutputStatus ( 'wrappast', WrapPastColumn )
  endif
  ! stop

  !---------------- Task (3c) ------------------
  ! Are we sharing our PCF with level 1?
  ! Are we part of a parallel run? If so, master or slave?
  ! If sharing a single PCF with level 1, we need to move any "mobile" PCF ids
  if ( sharedPCF ) call MovePCFIDs
  ! Setup the parallel stuff.  Register our presence with the master if we're a
  ! slave.
  if ( parallel%master ) call TransmitSlaveArguments
  if ( parallel%master .and. parallel%myTid <= 0 ) &
    & call MLSL2Message ( MLSMSG_Error, ModuleName, &
    & 'master Tid <= 0; probably pvm trouble' )
  if ( parallel%fwmParallel .and. parallel%master .and. parallel%chunkRange == ' ' ) &
    & call MLSL2Message ( MLSMSG_Error, ModuleName, &
    & 'fwmParallel mode can only be run for a single chunk' )
  if ( parallel%fwmParallel ) &
    & call MLSL2Message ( MLSMSG_Error, ModuleName, &
    & 'The fwmParallel option is currently inoperative, it needs significant work to fix - NJL' )
  call init_tables
  status = 0
  if ( parallel%slave ) then
    if ( parallel%masterTid <= 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'masterTid of this slave <= 0' )
    call dump_settings
    call InitParallel ( iChunks(1), slaveMAF )
    if ( parallel%myTid <= 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'slave Tid <= 0; probably pvm trouble' )
  end if
  !---------------- Task (4) ------------------
  ! Open the L2CF
  ! First, set up type checking for it
  status = InitializeMLSFile(MLSL2CF, content = 'l2cf', name='<STDIN>', &
      & type=l_ascii, access=DFACC_RDONLY, recordLength=recl, &
      & PCBottom=MLSPCF_L2CF_Start, PCTop=MLSPCF_L2CF_Start)
  MLSL2CF%FileID%f_id = l2cf_unit
  ! print *, 'l2cf from line? ', trim(line)
  ! print *, 'len_trim(line) ', len_trim(line)
  if ( line /= ' ' ) then
    MLSL2CF%name = line
    call mls_openFile( MLSL2CF, status )
    if ( status /= 0 ) then
      MLSL2CF%name = trim(line) // L2CFNAMEEXTENSION
      call MLS_OpenFile(MLSL2CF, status)
    end if
    if ( status /= 0 ) then
      call io_error ( "While opening L2CF", status, line )
      call MLSL2Message ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file: " // trim(line), MLSFile=MLSL2CF )
    else if(switchDetail(switches, 'pro') >= 0) then
      call announce_success( MLSL2CF%name, l2cf_unit )
    end if
    inunit = l2cf_unit
  else if ( TOOLKIT .and. .not. showDefaults ) then
    MLSL2CF%name = '' ! To force reference to PCF entry
    MLSL2CF%type = l_tkgen
    call mls_openFile(MLSL2CF, status)
    ! call dump(MLSL2CF)
    inunit = MLSL2CF%FileID%f_id
    if(status /= 0) then
      call output( 'Non-zero status returned from open_MLSCF: ', &
      & advance='no')
      call output(status, advance='yes')
      call MLSL2Message ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file named in pcf", MLSFile=MLSL2CF )
    else if(switchDetail(switches, 'pro') >= 0) then
      call announce_success( MLSL2CF%name, inunit )
    end if
  end if
  error = status
  if ( MLSL2CF%name == '<STDIN>' ) then
    MLSL2CF%FileID%f_id = -1
  else
    inunit = MLSL2CF%FileID%f_id
    ! We will store the l2cf file name in the Boolean 'l2cf'
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'l2cf', MLSL2CF%name, &
      & countEmpty=.true., &
      & inseparator=runTimeValues%sep )
  end if
  if (inunit /= -1) call AddInUnit(inunit)
  numfiles = AddFileToDataBase(FileDatabase, MLSL2CF)

  Show_Sys_Memory = Show_Sys_Memory .or. switchDetail(switches, 'mem') > -1
  if ( switchDetail(switches, 'memkb', options='-fc') > -1 ) then
    Show_Sys_Memory = .true.
    sys_memory_ch = 'kB'
    sys_memory_convert = 1.0
  elseif ( switchDetail(switches, 'memmb', options='-fc') > -1 ) then
    Show_Sys_Memory = .true.
    sys_memory_ch = 'MB'
    sys_memory_convert = 1.0e-3
  elseif ( switchDetail(switches, 'memgb', options='-fc') > -1 ) then
    Show_Sys_Memory = .true.
    sys_memory_ch = 'GB'
    sys_memory_convert = 1.0e-6
  endif
  call time_now ( t1 )

  if( switchDetail(switches, 'opt') > -1 .or. showDefaults ) then
    do j=1, size(current_version_id)
      call output(trim(current_version_id(j)), advance='yes')
    end do
    ! Don't dump_settings more than once
    if ( .not. parallel%slave ) call dump_settings
  end if
  if ( showDefaults ) stop
  ! In case we set skip retrievals on the commandline
  ! we won't want to accidentally override it later
  L2Options%SkipRetrievalOriginal = L2Options%SkipRetrieval
  call begin('starting mlsl2')
  !---------------- Task (5) ------------------
  if (error == 0) then
    ! Parse the L2CF, producing an abstract syntax tree
    call init_parser_table ( parser_table )
    call configuration ( root, parser_table )
    call destroy_parser_table ( parser_table )
    call clean_up_parser
  else
    root = -1
  end if
  ! call configureSayTime ( Coda=' ** (in s) **' )
  call dump ( Time_config )
  call dump ( sayTime_config )
  if ( timing ) call SayTime ( 'Parsing the L2CF', cumulative=.false. )

  !---------------- Task (6) ------------------
  if ( TOOLKIT .and. error==0) then
    call mls_closeFile ( FileDatabase(numfiles), error )
  else
    if ( inunit >= 0 ) close ( inunit )  ! Don't worry about the status
  end if
  if ( error /= 0) then
    call output ( &
      'An io error occurred with the l2cf -- there is no abstract syntax tree', &
      advance='yes' )
    call dump(MLSL2CF)
  else if ( root <= 0 ) then
    call output ( &
      'A syntax error occurred -- there is no abstract syntax tree', &
      advance='yes' )
      error = 1
  else
    if ( dump_tree >= 0 ) then
      call output ( 'Begin un-type-checked abstract syntax tree:', &
        & advance='yes' )
      call print_interesting_subtrees ( dump_tree )
      call output ( 'End un-type-checked abstract syntax tree', &
        & advance='yes' )
    end if

    ! Check that supra-syntactic conditions are met, e.g. correct
    ! types for fields of commands, correct command order, etc.
    call time_now ( t1 )
    call check_tree ( root, error, first_section )
    if ( timing ) call SayTime ( 'Type checking the L2CF', cumulative=.false. )
    if ( do_dump > 0 ) call dump_decl ( do_dump )
    if ( toggle(syn) ) then
      call output ( 'Begin type-checked abstract syntax tree:', advance='yes' )
      call print_interesting_subtrees ( levels(syn), type_name=get_type )
      call output ( 'End type-checked abstract syntax tree', advance='yes' )
    end if

    t2 = t0
    call add_to_section_timing( 'main', t2 )

    if(error /= 0) then
       call MLSL2Message( MLSMSG_Error, ModuleName, &
       & 'error in check_tree: probably need to repair l2cf ', MLSFile=MLSL2CF )
    end if

    if ( checkLeak ) then
      call time_now ( t1 )
      call leakCheck ( root )
      t2 = t0
      call add_to_section_timing( 'main', t2 )
    end if

! call SwitchOutput ( 'stdout' )
! call Dump( OutputOptions )
! call DumpConfig
! call RevertOutput

    !---------------- Task (7) ------------------
    if ( error == 0 .and. first_section /= 0 .and. .not. checkl2cf ) then
      ! Now do the L2 processing.
      ! stop-early flags => no writing, no retrieval
      skipDirectwrites = ( skipDirectwrites .or. stopAfterSection /= ' ' )
      L2Options%SkipRetrieval = ( L2Options%SkipRetrieval .or.  stopAfterSection /= ' ' )
      skipDirectwritesOriginal = skipDirectwrites
      OriginalOptions = L2Options
      call time_now ( t1 )
      if ( timing ) &
        & call output ( "-------- Processing Begun ------ ", advance='yes' )
      call walk_tree_to_do_MLS_L2 ( root, error, first_section, countChunks, &
        & FileDatabase )
      if ( timing ) then
        call output ( "-------- Processing Ended ------ ", advance='yes' )
        call SayTime ( 'Processing mlsl2', t1=0., cumulative=.false. )
      end if
    end if
  end if

  call time_now ( t0 )
  t1 = t0
  ! Moved up here because we're having chunks die during Task (8)
  ! (Why don't you find out why and fix it?)
  call add_to_section_timing( 'main', t0 )
  if ( switchDetail(switches, 'time') >= 0 ) then
    call output('(Now for the timings summary)', advance='yes')
    call dump_section_timings
  endif
  
  ! Tell wrapper script we Finished by way of noteFile
  if ( parallel%slave .and. len_trim(noteFile) > 0 ) then
    call output('Telling wrapper script we Finished', advance='yes')
    call write_textFile ( noteFile, 'Finished' )
  endif
  !---------------- Task (8) ------------------
  if ( .not. parallel%slave .or. slavesCleanUpSelves ) then
    call destroyEmpiricalGeometry
    call destroy_char_table
    call output('Destroyed char table', advance='yes')
    call destroy_hash_table
    call output('Destroyed hash table', advance='yes')
    call destroy_string_table
    call output('Destroyed string table', advance='yes')
    call destroy_symbol_table
    call output('Destroyed symbol table', advance='yes')
    call deallocate_decl
    call output('Deallocated decl', advance='yes')
    call deallocate_tree
    call output('Deallocated tree', advance='yes')
    call FreePVMArgs
    call output('Freed PVM args', advance='yes')
    if ( parallel%slave .and. &
      & (SKIPDIRECTWRITES .or. L2Options%SkipRetrieval .or. sectionsToSkip /= ' ' ) ) then
      ! call mls_h5close(error)
      ! call MLSMessageExit
    else if ( error == 0 ) then
      if ( .not. parallel%slave ) then
        call Deallocate_FileDatabase(FileDatabase)
        call output('Deallocated FileDatabase', advance='yes')
        call mls_h5close(error)
        call output('Closed hdf5 library', advance='yes')
      else
        call MLSMessage ( MLSMSG_Warning, moduleName, &
          & "We are a slave, unable to mls_*close in this version" )
      endif
      if (error /= 0) then
         call MLSL2Message ( MLSMSG_Error, moduleName, &
          & "Unable to MLS_Close" )
      end if
    end if
  end if
  if ( Toolkit ) then
    ! Coda
    numResolutions = 2
    numLayers = 2
    resolutionList(1) = PGSD_DEM_30ARC
    resolutionList(2) = PGSD_DEM_90ARC
    layerList(1) = PGSD_DEM_ELEV
    layerList(2) = PGSD_DEM_WATER_LAND
    status = PGS_DEM_Close ( resolutionList, numResolutions, &
      & layerList, numLayers )
    call outputNamedValue( 'PGS_DEM_Close status', status )
  end if
  if ( AllocateLogUnit > 0 ) close( Unit=AllocateLogUnit )
  if ( timing ) call SayTime ( 'Closing and deallocating', cumulative=.false. )
  call add_to_section_timing( 'main', t0 )
  if ( trackAllocates > 0 ) call ReportLeaks ( "At end of program execution..." )
  call FinalMemoryReport ( isFinal=.true. )
  call finish ( 'ending mlsl2' )
  if ( parallel%slave .and. len_trim(noteFile) > 0 .and. WAITFORSCRIPT ) then
    ! cycle endlessly until wrapper script kills us
    call output( ' Waiting for wrapper script to kill us', advance='yes' )
    call flushStdout
    do
      call usleep ( 10000*parallel%delay )
    enddo
  endif
  if( error /= 0 .or. STOPWITHERROR ) then
     call MLSMessageExit(1)
  else if( NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave ) then
     call MLSMessageExit( NORMAL_EXIT_STATUS )
  else
     call MLSMessageExit
  end if

contains

  subroutine Print_Interesting_Subtrees ( Level, Type_Name )
    use Tree, only: Dump_Tree_Node, Node_ID, Print_Subtree
    use Tree_Types, only: N_CF, N_If, N_Select
    integer, intent(in) :: Level ! If > 0, print the entire tree rooted at
                               ! Root, else print only subtrees with roots
                               ! having tree indices of n_cf, n_if and n_select
    optional :: Type_Name
    interface
      integer function Type_Name ( Decor )
      ! Return the string index to print for the decoration
        integer, intent(in) :: Decor ! Tree(where)%Decor
      end function Type_Name
    end interface
    integer :: I ! Subtree index
    if ( level > 0 ) then
      call print_subtree ( root, 0, type_name=type_name )
    else
      call output ( root, 5 )
      call output ( ':' )
      call dump_tree_node ( root, 0, type_name=type_name, advance='yes' )
      do i = 1, nsons(root)
        select case ( node_id(subtree(i,root)) )
        case ( N_CF, N_If, N_Select )
          call print_subtree ( subtree(i,root), 1, type_name=type_name )
        end select
      end do
    end if
  end subroutine Print_Interesting_Subtrees

  ! The following offer some information on the options
  ! available either on the command-line or via the PCF
  ! Note the unashamed use of 'print' statements which
  ! are officially discouraged in favor of calls to MLSMessage.
  ! Unfortunately, we have not yet decided which method to use
  ! until *after* processing all the options.

  subroutine Switch_usage
    print *, 'Switch usage: -S"sw1,sw2,..,swn" or -Ssw1 -Ssw2 ..'
    print *, ' where each of the swk may be one of the following'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   !  See mlspgs/notes/switches
   ! === (start of automatic usage lines) ===
    print *, '  A => AntennaPatterns'
   ! === (end of automatic usage lines) ===
    stop
  end subroutine Switch_usage

  subroutine Option_usage
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [L2CF-name]'
    print *, ' Options:'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   !  See mlspgs/notes/options
   ! === (start of automatic option lines) ===
    print *, '  -A: Dump the un-decorated abstract syntax tree.'
   ! === (end of automatic option lines) ===
    stop
  end subroutine Option_usage

  subroutine Dump_settings
  ! Show current run-time settings resulting from
  ! command-line, MLSL2Options, etc.
  ! Intel's ifort does not yet support the following USE statement
  !     use, intrinsic :: ISO_FORTRAN_ENV, only: Compiler_Options, Compiler_Version
    use Printit_m, only: Get_Config
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
    character(len=255) :: Command ! Command that executed the program
    integer :: LogFileUnit
    character(len=8) :: string
    ! Executable
    ! Always show how we were invoked and what the l2cf is
    call getarg ( 0, command )
    if ( command /= '' ) then
      call output ( ' mlsl2 invoked as: ', advance='no' )
      call output ( trim(command), advance='yes' )
    end if
    call output ( ' mlsl2 called with command line options: ', advance='no' )
    call output ( trim(L2Options%Originalcmds), advance='yes' )
    call output ( ' l2cf file:', advance='no' )
    call blanks ( 4, advance='no' )
    call output ( trim(MLSL2CF%name), advance='yes' )
!     call output ( 'Compiler options ' )
!     call output ( compiler_options(), advance='yes' )
!     call output ( 'Compiler version ' )
!     call output ( compiler_version(), advance='yes' )
    ! Now show a nicely-formatted list of all run-time options
    if( SwitchDetail( switches, 'opt') > 0 .or. showDefaults ) then
      call blanks( 80, fillChar='-', advance='yes' )
      call headline( 'Summary of run time options', fillChar='-', before='*', after='*' )
      call outputNamedValue ( 'Use toolkit panoply', toolkit, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Error threshold before halting', Quit_Error_Threshold, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Status on normal exit', normal_exit_status, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'PCF shared with level 1?', SHAREDPCF, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Range of chunks', trim_safe(parallel%chunkRange), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Unique job ID', trim(uniqueID), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Show system memory usage?', &
        & Show_Sys_Memory, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Is this the master task in pvm?', parallel%master, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      if ( parallel%master .or. showDefaults ) then
        call outputNamedValue ( 'Master task number', parallel%myTid, advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        call outputNamedValue ( 'Command line sent to slaves', trim_safe(parallel%pgeName), advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        call outputNamedValue ( 'Command to queue slave tasks', trim_safe(parallel%submit), advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        call outputNamedValue ( 'Maximum failures per chunk', parallel%maxFailuresPerChunk, advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        call outputNamedValue ( 'Maximum failures per machine', parallel%maxFailuresPerMachine, advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        call outputNamedValue ( 'Sleep time in masterLoop (mus)', parallel%delay, advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      end if
      call outputNamedValue ( 'Is this a slave task in pvm?', parallel%slave, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      if ( parallel%slave ) then
        call outputNamedValue ( 'Master task number', parallel%masterTid, advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
        if ( len_trim(noteFile) > 0 ) then
          call output ( '*  Note file:', advance='yes' )
          call output ( '*    ' // trim(noteFile), advance='yes' )
        end if
      end if
      if ( len_trim(stopAfterSection) > 0 ) then
        call outputNamedValue ( 'Stop after section', trim_safe(stopAfterSection), advance='yes', &
          & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      endif
      call outputNamedValue ( 'Preflight check paths?', checkPaths, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Need L1B files?', NEED_L1BFILES, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Reading Aura L1B files?', AURA_L1BFILES, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Skip all direct writes?', SKIPDIRECTWRITES, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Skip all retrievals?', L2Options%SkipRetrieval, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Filling unretrieved states with', STATEFILLEDBYSKIPPEDRETRIEVALS, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Using wall clock instead of cpu time?', time_config%use_wall_clock, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call get_string ( lit_indices(sectionTimingUnits), string, strip=.true. )
      call outputNamedValue ( 'Summarize time in what units', trim_safe(string), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Number of switches set', numSwitches, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Standard output unit', PrUnitName(outputOptions%prunit), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call get_config ( logFileUnit = logFileUnit )
      call outputNamedValue ( 'Log file unit', PrUnitname(logFileUnit), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Crash on any error?', MLSMessageConfig%crashOnAnyError, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Suppress identical warnings after', MLSMessageConfig%limitWarnings, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Restart warnings count at each phase?', restartWarnings, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Set error before stopping?', StopWithError, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      if ( specialDumpFile /= ' ' ) then
        call outputNamedValue ( 'Save special dumps to', trim(specialDumpFile), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Default hdf version for l1b files', level1_hdfversion, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Default hdfeos version on reads', default_hdfversion_read, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Default hdfeos version on writes', default_hdfversion_write, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Max Chunk size Appended to l2gp files', maxChunkSize, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      end if
      if ( len_trim(AllocFile%Name) > 0 ) then
        call outputNamedValue ( 'Log Allocates to', trim(AllocFile%Name), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      endif
      call outputNamedValue ( 'Avoiding unlimited dimensions in directwrites?', &
        & avoidUnlimitedDims, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Allow overlaps outside proc. range?', &
        & (/ChunkDivideConfig%allowPriorOverlaps, ChunkDivideConfig%allowPostOverlaps/), advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Time some processes?', timing, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Is this run in forward model parallel?', parallel%fwmParallel, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call outputNamedValue ( 'Avoid creating file on first directWrite?', patch, advance='yes', &
        & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
      call Blanks( 80, fillChar='-', advance='yes' )
      ! Dump special settings and configurations
      if ( len_trim(switches) > 0 ) &
        & call Dump( switches, 'All Switches' )
      if ( len_trim(MLSNamesAreDebug) > 0 ) &
        & call Dump( MLSNamesAreDebug, 'Modules to debug' )
      if ( len_trim(MLSNamesAreVerbose) > 0 ) &
        & call Dump( MLSNamesAreDebug, 'Modules are verbose' )
      if ( len_trim(SectionsToSkip) > 0 ) &
        & call Dump( SectionsToSkip, 'Sections to Skip' )
      if ( len_trim(PhasesToSkip) > 0 ) &
        & call Dump( PhasesToSkip, 'Phases to Skip' )
      if ( L2Options%Overridden ) call DumpOptions
      call SwitchOutput ( 'stdout' )
      call Dump( OutputOptions )
      call Dump( StampOptions )
      call DumpConfig
      call RevertOutput
      call Blanks( 80, fillChar='-', advance='yes' )
      ! Many of these are defined in the PCF or in the global settings
      ! Neither of these sections have been completed yet
      ! call DumpGlobalAttributes
    end if
  end subroutine Dump_settings

  ! ---------------------------------------------  MovePCFIDs  -----
  subroutine MovePCFIDs
    ! We need to move any "mobile" PCFids
    include 'sharedpcf2.f9h'
    print *, 'Moved PCFIds'
    print *, 'mlspcf_l2cf_start: ', mlspcf_l2cf_start

  end subroutine MovePCFIDs

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, unit_number )
    character(len=*), intent(in)   :: Name
    integer, intent(in)            :: unit_number

    call output ( 'Level 2 configuration file ' )
    call output ( 'unit number : ' )
    call blanks(4)
    call output ( unit_number, advance='no')
    call blanks(10)
    call output ( 'name : ' )
    call blanks(4)
    call output ( trim(Name), advance='yes')
  end subroutine announce_success

!-----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
!---------------------------------------------------------------------------

end program MLSL2

! $Log$
! Revision 2.236  2021/08/20 15:58:17  pwagner
! -Salloc prints Memory summary while tracing
!
! Revision 2.235  2021/07/08 23:34:08  pwagner
! Wider use of MLSL2Message
!
! Revision 2.234  2020/04/30 23:29:42  pwagner
! Add Final call to FinalMemoryReport at end
!
! Revision 2.233  2020/04/09 23:17:46  pwagner
! Improved wording in Dump_settings
!
! Revision 2.232  2019/08/01 23:46:50  pwagner
! Some Housekeeping
!
! Revision 2.231  2019/04/11 23:43:43  pwagner
! cmdline options --help,--version, and -S'?' now print only relevant stuff
!
! Revision 2.230  2019/03/08 17:15:49  pwagner
! Still correcting mis-directed output
!
! Revision 2.229  2019/02/13 17:31:27  pwagner
! change outputOptions%prunit to L2Options only if default value is Overriden
!
! Revision 2.228  2019/01/31 19:23:05  pwagner
! Removed more unused stuff
!
! Revision 2.227  2018/12/07 00:20:47  pwagner
! If cmdline says to skip retrievals, must skip even if l2cf says otherwise
!
! Revision 2.226  2018/09/13 20:25:19  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.225  2018/05/22 23:08:40  pwagner
! Dumped settings include GlobalAttributes
!
! Revision 2.224  2017/12/22 00:32:00  pwagner
! Try harder to make slaves flush stdout buffer
!
! Revision 2.223  2017/03/23 16:57:38  pwagner
! Improved Dump_settings appearance
!
! Revision 2.222  2017/01/25 18:10:22  pwagner
! May skip certain Phases named in phasesToSkip cmdline opt
!
! Revision 2.221  2017/01/11 23:56:53  pwagner
! Dumps both time configs; Dump MLSFile type if l2cf error
!
! Revision 2.220  2016/11/08 17:30:57  pwagner
! Use SayTime subroutine from time_m module
!
! Revision 2.219  2016/11/04 19:32:11  pwagner
! Locate calculation of MLSMSG_Severity_to_quit appropriately
!
! Revision 2.218  2016/05/27 00:05:22  pwagner
! Should now correctly process options containing an embedded space
!
! Revision 2.217  2016/03/18 18:00:15  pwagner
! Dont let -Smess override effect of --loconic option
!
! Revision 2.216  2016/02/29 19:51:07  pwagner
! Usleep got from machine module instead of being an external
!
! Revision 2.215  2015/10/06 00:22:09  pwagner
! Automatically store l2cf name in runtime Boolean 'l2cf'
!
! Revision 2.214  2015/09/24 22:07:57  pwagner
! Dump_settings shows max chunk size appended to l2gp files
!
! Revision 2.213  2015/07/16 22:13:15  pwagner
! Dump of settings says whether it will stop after a named section
!
! Revision 2.212  2014/12/10 23:04:54  pwagner
! Removed unused stuff
!
! Revision 2.211  2014/09/11 18:28:35  pwagner
! Added -Smem[units] switch to show sys memory
!
! Revision 2.210  2014/09/05 01:07:49  vsnyder
! Remove ProcessID stuff.  Destroy Empirical Geometry explicitly
!
! Revision 2.209  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.208  2014/09/02 18:18:23  pwagner
! Uses noteFile mechanism for telling our wrapper script we finished
!
! Revision 2.207  2014/06/30 23:26:58  pwagner
! Can log allocations/deallocations to separate file
!
! Revision 2.206  2014/05/20 23:56:53  vsnyder
! New parser gets its tables from an argument instead of an include
!
! Revision 2.205  2014/03/25 18:19:00  pwagner
! Moved init_tables before dump_settings to avoid crash
!
! Revision 2.204  2014/03/20 01:42:32  vsnyder
! Get Get_Type from Intrinsic, do init_tables later
!
! Revision 2.203  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.202  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.201  2013/12/12 02:10:58  vsnyder
! Add 'debug' to init_lexer, move syntax tree printing into internal subroutine
!
! Revision 2.200  2013/11/26 22:40:51  vsnyder
! Change -A to -A[n] with n>0 meaning dump the entire tree, including the
! type-checking stuff, and n==0 or absent meaning dump only the parser output.
!
! Revision 2.199  2013/11/04 22:56:59  pwagner
! Print which modules are berbose, debugged
!
! Revision 2.198  2013/10/09 23:46:38  vsnyder
! Pass Get_Type to tree dumper
!
! Revision 2.197  2013/09/25 18:51:55  pwagner
! Closes DEM stuff before exiting
!
! Revision 2.196  2013/09/06 20:58:18  pwagner
! Trying again to prevent slaves from logging to toolkit
!
! Revision 2.195  2013/09/04 17:34:52  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.194  2013/08/23 02:52:13  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 2.193  2013/08/17 00:22:14  pwagner
! New cmdline arg relaxes some for non-Aura l1b datasets
!
! Revision 2.192  2013/06/28 19:18:42  pwagner
! -Smess[n] and -Sparent switches; other tweaks to output
!
! Revision 2.191  2013/02/12 18:14:38  pwagner
! Removed SIPS_VERSION
!
! Revision 2.190  2012/12/04 00:16:09  pwagner
! Removed confisuion-causing OUTSIDEOVERLAPS and its cmdline option
!
! Revision 2.189  2012/08/21 23:54:01  pwagner
! Fixed certain problems when showing settings
!
! Revision 2.188  2012/08/16 17:48:47  pwagner
! Exploit new features to fleexibly change how and when to log output
!
! Revision 2.187  2012/07/18 00:39:34  pwagner
! Consistent with module parameters for prUnit, LogFileUnit
!
! Revision 2.186  2012/07/10 15:22:25  pwagner
! Moved stuff to MLSL2Options module so it works right
!
! Revision 2.185  2012/07/02 20:34:41  pwagner
! Avoid dumping settings more than one time
!
! Revision 2.184  2012/06/27 18:02:59  pwagner
! May overwrite command line options with options field to phase spec
!
! Revision 2.183  2012/06/06 20:36:37  vsnyder
! Use Set_Toggles for -aefglpt options
!
! Revision 2.182  2012/04/26 23:12:45  pwagner
! --crash turns off neverCrash; maybe we should change default for neverCrash to .false.
!
! Revision 2.181  2012/04/20 01:30:38  vsnyder
! Use Begin, Finish instead of Output_Date_and_Time
!
! Revision 2.180  2012/03/28 20:06:44  pwagner
! Removed stgmem option--slave tasks lost ability to join quantities
!
! Revision 2.179  2012/02/09 02:45:26  vsnyder
! Remove -C option -- I didn't notice --crash was already there
!
! Revision 2.178  2012/02/09 01:48:11  vsnyder
! Add -C option to crash on any error, to get traceback for debugging
!
! Revision 2.177  2011/06/29 21:43:59  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.176  2010/11/05 22:36:16  pwagner
! Consistent with new unquote api
!
! Revision 2.175  2010/09/17 00:11:19  pwagner
! Fixed misuse of switchDetail function
!
! Revision 2.174  2010/08/06 23:08:48  pwagner
! Pass Hessians, matrices to DumpCommand
!
! Revision 2.173  2010/05/23 03:25:37  honghanh
! Use AddInUnit instead of inunit to adapt with change in the string_table
!
! Revision 2.172  2010/02/02 01:41:03  vsnyder
! Move declaration of Id to a place where it is more difficult for a
! compiler to notice it's not actually referenced.
!
! Revision 2.171  2009/10/21 16:59:10  pwagner
! No longer insist on unlimited dims when given chunkrange
!
! Revision 2.170  2008/01/08 00:22:02  pwagner
! Levels 1 and 2 can use same shared PCF now if --shared option set
!
! Revision 2.169  2007/09/13 21:52:22  pwagner
! Reduced sips limitwarnings to 4
!
! Revision 2.168  2007/09/06 23:35:16  pwagner
! slaves need not do own cleanup
!
! Revision 2.167  2007/08/23 22:18:35  pwagner
! 'walk' switch prints walkback fo calling stack even for warnings
!
! Revision 2.166  2007/08/17 00:36:29  pwagner
! Removed now-redundant [n]everCrash option
!
! Revision 2.165  2007/07/27 00:18:20  vsnyder
! Print the command that invoked MLSL2
!
! Revision 2.164  2007/06/07 20:38:45  pwagner
! Should prevent unit collisions (not as good as get_lun)
!
! Revision 2.163  2007/02/14 20:48:59  pwagner
! Slaves bypass mls_h5close to avoid another bomb
!
! Revision 2.162  2007/02/14 17:31:45  pwagner
! Moved timing summary before statements possibly killing slaves
!
! Revision 2.161  2007/02/07 20:58:03  pwagner
! Permit slaves to write unbuffered stdout when master does
!
! Revision 2.160  2007/01/12 00:35:08  pwagner
! Renamed routine outputNamedValue; may set host name and unbuffered stdout
!
! Revision 2.159  2006/10/09 18:38:41  pwagner
! Trims parallel%chunkRange before outputting
!
! Revision 2.158  2006/10/05 23:32:14  pwagner
! skipSections can skip named sections
!
! Revision 2.157  2006/09/21 18:46:49  pwagner
! Reduce level of dumps in SIDS version
!
! Revision 2.156  2006/08/14 16:21:22  pwagner
! pass chunk range correctly to slave tasks
!
! Revision 2.155  2006/08/10 21:46:49  pwagner
! --chunk commandline option now synonym for --chunkRange
!
! Revision 2.154  2006/08/05 02:13:33  vsnyder
! Add 'where' argument for ReportLeaks
!
! Revision 2.153  2006/07/29 03:42:09  vsnyder
! New --memtrack interpretation
!
! Revision 2.152  2006/07/27 03:49:12  vsnyder
! Detect --leak option, attach leak checker
!
! Revision 2.151  2006/07/21 20:10:37  pwagner
! Can fill state even if skipping retrievals; select what section to stop after
!
! Revision 2.150  2006/06/28 00:00:18  pwagner
! Uses new outputNamedValue
!
! Revision 2.149  2006/06/24 23:11:29  pwagner
! prunit now a component of outputOptions
!
! Revision 2.148  2006/04/20 23:22:39  pwagner
! Show both kinds of allowed extra-range overlaps
!
! Revision 2.147  2006/04/11 23:28:49  pwagner
! Whether to allow overlaps outside of processing range
!
! Revision 2.146  2006/03/04 00:17:20  pwagner
! May skip retrieval, directWrites depending on runtime Booleans
!
! Revision 2.145  2006/02/21 19:15:17  pwagner
! Uses switchDetail only now
!
! Revision 2.144  2006/02/10 21:14:32  pwagner
! May specify skipRetrivel for particular Phases; dumps may go to special dumpfile
!
! Revision 2.143  2005/11/17 20:11:46  pwagner
! Was printing wrong thing for parallel sleeptime
!
! Revision 2.142  2005/09/22 23:38:54  pwagner
! time_config and retry_config now hold config settings
!
! Revision 2.141  2005/08/19 23:29:04  pwagner
! Wider use of SwitchDetail function
!
! Revision 2.140  2005/07/21 23:39:49  pwagner
! use_wall_clock tied to SIPS_VERSION
!
! Revision 2.139  2005/06/29 17:57:40  pwagner
! FILESTRINGTABLE set to TRUE
!
! Revision 2.138  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.137  2005/06/14 20:45:22  pwagner
! Many changes to accommodate the new fields in MLSFile_T
!
! Revision 2.136  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.135  2005/04/12 18:12:30  pwagner
! SIPS version limits warnings to 50; reads scalar intfromchars
!
! Revision 2.134  2005/04/01 00:15:37  pwagner
! Automatcially remove slv switch from SIPS version
!
! Revision 2.133  2005/03/15 23:58:42  pwagner
! Sets MLSMessageConfig appropriately for slaves
!
! Revision 2.132  2005/03/12 00:49:18  pwagner
! -w option added
!
! Revision 2.131  2005/03/03 00:23:42  pwagner
! Added -Rs1,s2,.. Removeswitches option
!
! Revision 2.130  2005/01/22 00:39:08  pwagner
! Reversed buggy evercrash option
!
! Revision 2.129  2004/12/27 23:05:27  pwagner
! Commented out useless prints; warns on checkPaths option
!
! Revision 2.128  2004/12/14 21:55:37  pwagner
! May skip sections, stop early
!
! Revision 2.127  2004/10/30 00:28:00  vsnyder
! Commented out some unused stuff to make NAG stop nagging
!
! Revision 2.126  2004/10/11 16:57:03  pwagner
! Bug fix for options; slaves output not logged; prints proc start, end date_and_time
!
! Revision 2.125  2004/08/19 00:20:21  pwagner
! New crash-related, defaults options
!
! Revision 2.124  2004/08/16 17:12:20  pwagner
! Fixed clearonallocate setting
!
! Revision 2.123  2004/08/05 22:47:47  pwagner
! New --chunkRange option to run selected chunks in parallel mode
!
! Revision 2.122  2004/08/04 23:19:57  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.121  2004/07/08 22:48:44  pwagner
! Made SIPS_VERSION public
!
! Revision 2.120  2004/06/29 00:10:17  pwagner
! Exploit catlist function
!
! Revision 2.119  2004/04/27 23:50:24  pwagner
! Added SKIPDIRECTWRITES option
!
! Revision 2.118  2004/04/15 22:48:30  pwagner
! Multiword options like maxFailuresPerChunk made case-insensitive
!
! Revision 2.117  2004/04/06 23:50:09  livesey
! Added clearOnAllocate flag
!
! Revision 2.116  2004/04/03 05:44:16  livesey
! Added memTrack option
!
! Revision 2.115  2004/03/24 23:54:06  pwagner
! Switched from h5_open/close_f to mls_open/close
!
! Revision 2.114  2004/03/24 01:03:34  livesey
! Increased tree size.
!
! Revision 2.113  2004/02/05 23:26:22  pwagner
! Added --cat option
!
! Revision 2.112  2004/01/09 00:22:12  pwagner
! Unsets avoidUnlimitedDims to bypass bug directWriting range of chunks
!
! Revision 2.111  2004/01/07 23:50:16  livesey
! Added error message about fwmParallel being broken at the moment
!
! Revision 2.110  2003/12/11 22:56:04  pwagner
! Transmits --idents option to slaves; unquotes switches
!
! Revision 2.109  2003/12/07 23:06:29  pwagner
! Should not dump all the chunks agagin and again for each slaves chunk
!
! Revision 2.108  2003/12/05 00:40:11  pwagner
! Added patch option, section timing units
!
! Revision 2.107  2003/11/15 00:33:37  pwagner
! New commandline opts: maxfailuresperchunk, maxfailurespermachine
!
! Revision 2.106  2003/11/14 23:37:13  pwagner
! Lets user change masterLoop delay via commandline option
!
! Revision 2.105  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.104  2003/11/05 21:27:54  pwagner
! Can enter range of chunks to be processed instead of single
!
! Revision 2.103  2003/10/21 00:02:32  pwagner
! Revert to writing just once if toolkit and parallel
!
! Revision 2.102  2003/10/14 18:17:02  pwagner
! Fixed problem with reducing switches to unique list
!
! Revision 2.101  2003/10/09 23:57:35  pwagner
! A few more SIPS-related tweaks
!
! Revision 2.100  2003/10/09 23:35:21  pwagner
! Treats SIPS version special; switches ,-separated
!
! Revision 2.99  2003/09/05 23:22:08  pwagner
! Takes in new --skipRetrieval option
!
! Revision 2.98  2003/08/01 20:26:01  pwagner
! gets slave pge name from command line
!
! Revision 2.97  2003/06/13 20:02:50  vsnyder
! Put snoopMAF before snoop, so it can be found, futzing
!
! Revision 2.96  2003/06/09 22:49:33  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.95  2003/05/14 00:59:40  livesey
! Increased hash table size.
!
! Revision 2.94  2003/05/13 04:48:20  livesey
! Added stgmem option
!
! Revision 2.93  2003/02/27 18:40:29  pwagner
! recl passed to let NAG open l2cf with long lines
!
! Revision 2.92  2003/02/08 00:30:57  pwagner
! Increased default RECL due to new l2cfs in lib
!
! Revision 2.91  2002/12/19 11:54:45  livesey
! Upped tree size
!
! Revision 2.90  2002/12/10 00:39:06  pwagner
! Announces success a la read_apriori
!
! Revision 2.89  2002/12/06 22:33:41  livesey
! Added the snoop name stuff
!
! Revision 2.88  2002/12/05 19:45:20  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.87  2002/12/04 01:18:21  pwagner
! First halting steps toward using filedatabase
!
! Revision 2.86  2002/11/13 01:08:40  pwagner
! Prints more info when called asked to print opts
!
! Revision 2.85  2002/10/30 00:56:51  livesey
! Increased size of hash table
!
! Revision 2.84  2002/10/08 17:41:20  livesey
! Now passes chunk argument to slaves if present.  Needed for
! fwmParallel mode.
!
! Revision 2.83  2002/10/05 00:44:14  livesey
! Included the FMWParallel stuff
!
! Revision 2.82  2002/10/03 23:00:03  pwagner
! You can set l1b, l2gp hdfversions on command line
!
! Revision 2.81  2002/09/24 18:17:20  pwagner
! Consistent with add_to_section_timing now calling time_now at its end
!
! Revision 2.80  2002/08/29 21:47:51  livesey
! Added subblock option
!
! Revision 2.79  2002/08/21 18:22:26  vsnyder
! Revise processing of Lahey/Fujitsy runtime (-Wl) options
!
! Revision 2.78  2002/08/15 21:48:51  pwagner
! h5open(close)_f now called at start (end)
!
! Revision 2.77  2002/07/23 23:15:05  pwagner
! Moved dump_settings call after learning name of l2cf file
!
! Revision 2.76  2002/07/18 21:59:28  vsnyder
! Cosmetic changes, move some stuff around so PRUNIT is stdout in init_tables
!
! Revision 2.75  2002/05/28 22:34:59  livesey
! Increased tree size
!
! Revision 2.74  2002/05/24 17:27:18  pwagner
! Added recl=recl to second-chance open of l2cf
!
! Revision 2.73  2002/05/23 20:57:57  vsnyder
! Add --recl option, some cosmetic changes
!
! Revision 2.72  2002/05/14 00:26:49  livesey
! Larger table sizes (may be unnecessary)
!
! Revision 2.71  2002/04/24 20:21:13  livesey
! Upped some string sizes.
!
! Revision 2.70  2002/04/24 16:54:02  livesey
! Changes for submit option
!
! Revision 2.69  2002/03/20 00:48:29  pwagner
! Option -check just checks l2cf then quits
!
! Revision 2.68  2002/02/20 00:29:04  pwagner
! Retries FN+.l2cf; tracks successful l2cf file name
!
! Revision 2.67  2002/02/12 00:25:00  pwagner
! New switch -opt[n] and new --version option
!
! Revision 2.66  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.65  2002/01/18 18:55:25  livesey
! Added the --chunk option
!
! Revision 2.64  2002/01/09 22:56:46  livesey
! Tidied up parsing of master option.
!
! Revision 2.63  2002/01/09 00:00:48  pwagner
! Fixed small comment; added others explaining unavoidable use of print
!
! Revision 2.62  2001/12/13 23:21:20  livesey
! Added countChunks option
!
! Revision 2.61  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.60  2001/11/09 18:12:38  livesey
! Added --nckbk option.
!
! Revision 2.59  2001/10/12 23:11:28  pwagner
! Clarified what Processing time means?
!
! Revision 2.58  2001/10/09 22:37:55  livesey
! Increased tree size and hash table size to accomodate new
! spectroscopy database from Bill
!
! Revision 2.57  2001/10/04 23:50:25  livesey
! Added the ckbk option
!
! Revision 2.56  2001/10/04 00:16:45  pwagner
! Increased hash table size; added note that size(s) may need to grow
!
! Revision 2.55  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.54  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.53  2001/09/19 23:43:49  livesey
! Added --snoop option
!
! Revision 2.52  2001/07/19 22:00:09  pwagner
! Better behaved when l2cf is stdin
!
! Revision 2.51  2001/07/18 23:56:24  pwagner
! Gets error from close_mlscf
!
! Revision 2.50  2001/07/18 00:16:54  pwagner
! Better control over when to exit with status
!
! Revision 2.49  2001/07/16 23:43:15  pwagner
! With settable NORMAL_EXIT_STATUS
!
! Revision 2.48  2001/05/25 01:03:47  livesey
! Working parallel version
!
! Revision 2.47  2001/05/23 23:21:51  pwagner
! Did the same for options, switches
!
! Revision 2.46  2001/05/23 22:31:30  pwagner
! Automatically updates based on notes/switches
!
! Revision 2.45  2001/05/23 21:59:43  livesey
! Interim version, almost there
!
! Revision 2.44  2001/05/23 01:44:24  livesey
! Parallel code starting to fit into place
!
! Revision 2.43  2001/05/18 01:14:21  vsnyder
! Add a 'stop' in 'switch_usage', plus cosmetic changes
!
! Revision 2.42  2001/05/17 22:34:55  pwagner
! output and MLSMessage modules cooperate better w/o toolkit; switch_usage
!
! Revision 2.41  2001/05/15 23:46:07  pwagner
! Removed 2 settings from MLSL2Opts; now in switches
!
! Revision 2.40  2001/05/11 23:47:00  pwagner
! (Re)Sets prunit depending on toolkit
!
! Revision 2.39  2001/05/11 17:34:31  vsnyder
! Improve built-in usage display
!
! Revision 2.38  2001/05/09 23:33:00  pwagner
! Sets new MLSL2Options
!
! Revision 2.37  2001/05/08 20:33:41  vsnyder
! Moved usage display into a subroutine
!
! Revision 2.36  2001/05/07 23:30:51  pwagner
! Sets USESDPToolkit
!
! Revision 2.35  2001/05/07 21:53:28  vsnyder
! Improve built-in usage display
!
! Revision 2.34  2001/05/07 21:05:03  vsnyder
! Separated '[n]pcf' and [n]cfpcf'
!
! Revision 2.33  2001/05/07 18:16:11  vsnyder
! Added "do [not] output MLSMessage messages through the toolkit" option.
! Resolved conflicts on merge.
!
! Revision 2.32  2001/05/07 17:17:31  pwagner
! Calls MLSMessage to exit with error if check_tree fails
!
! Revision 2.31  2001/05/04 22:55:36  pwagner
! Added cascading negatives; new command line options
!
! Revision 2.30  2001/05/03 01:58:52  vsnyder
! Add error checking for --slave option
!
! Revision 2.29  2001/05/02 23:22:48  livesey
! Added parallel stuff
!
! Revision 2.28  2001/05/01 17:51:47  vsnyder
! Ignore Lahey/Fujitsu run-time library's command-line arguments
!
! Revision 2.27  2001/04/28 01:44:47  vsnyder
! Provide to set levels(emit)
!
! Revision 2.26  2001/04/27 20:59:16  vsnyder
! Change -G option to -S
!
! Revision 2.25  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.24  2001/04/24 23:04:42  vsnyder
! Add -f[digit] to set toggle(emit) and levels(emit)
!
! Revision 2.23  2001/04/21 01:44:05  vsnyder
! Make the timing message prettier
!
! Revision 2.22  2001/04/21 01:42:21  vsnyder
! Add -T option for timing
!
! Revision 2.21  2001/04/20 20:44:18  pwagner
! Sets MLSMSG_Severity_to_quit
!
! Revision 2.20  2001/04/17 22:08:56  pwagner
! Sets prunit according to OUTPUT_PRINT_UNIT
!
! Revision 2.19  2001/04/16 23:46:58  pwagner
! Gets PCF_FOR_INPUT flag from MLSL2Options
!
! Revision 2.18  2001/04/06 20:11:53  vsnyder
! Add -e option
!
! Revision 2.17  2001/04/05 01:33:46  vsnyder
! Increase initial sizes for sevaral paraer tables
!
! Revision 2.16  2001/03/28 01:29:48  vsnyder
! Add description of -G to -h output
!
! Revision 2.15  2001/03/16 21:47:57  vsnyder
! Add -G option
!
! Revision 2.14  2001/03/14 18:59:03  vsnyder
! Add K and k options to control whether the lexer capitalizes identifiers
!
! Revision 2.13  2001/03/08 00:39:37  vsnyder
! Improve some debugging output
!
! Revision 2.12  2001/03/02 02:38:17  vsnyder
! Expand LINE, alphabetize USEs
!
! Revision 2.11  2001/02/28 03:01:48  vsnyder
! Make presence of L2CF-name on command line take precedence over --[n]pcf
!
! Revision 2.10  2001/02/28 02:52:32  vsnyder
! Improve usage description
!
! Revision 2.9  2001/02/28 02:44:24  vsnyder
! Identify abstract syntax tree dumps, show default --[n]pcf in usage
!
! Revision 2.8  2001/02/28 01:59:29  vsnyder
! Access Open_MLSCF and Close_MLSCF from Obtain_MLSCF instead of Open_Init
!
! Revision 2.7  2001/02/23 02:39:56  vsnyder
! Add description of --pcf option to usage instructions.
!
! Revision 2.6  2001/02/23 02:38:34  vsnyder
! Open L2CF either by PCF or by Fortran OPEN or expect it on stdin
!
! Revision 2.5  2001/02/22 23:51:00  vsnyder
! Improved usage messages
!
! Revision 2.4  2001/02/22 23:05:12  vsnyder
! Display usage if -h, -H or -? option is present.
!
! Revision 2.3  2000/10/12 00:33:47  vsnyder
! Insert CVS variables and copyright notice
!
