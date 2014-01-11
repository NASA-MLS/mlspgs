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
MODULE MLSL2Options              !  Options and Settings for the MLSL2 program
!=============================================================================

  use HIGHOUTPUT, only: OUTPUTNAMEDVALUE
  use INTRINSIC, only: L_HOURS, L_MINUTES, L_SECONDS
  use MLSCOMMON, only: MLSFILE_T, MLSNAMESAREDEBUG, MLSNAMESAREVERBOSE
  use MLSFILES, only: WILDCARDHDFVERSION, HDFVERSION_4, HDFVERSION_5
  use MLSMESSAGEMODULE, only: MLSMESSAGECONFIG, &
    & MLSMSG_ERROR, MLSMSG_INFO, MLSMSG_TESTWARNING, &
    & MLSMSG_SEVERITY_TO_WALKBACK, MLSMSG_WARNING, &
    & SAYMESSAGE => MLSMESSAGE
  use MLSPCF2, only: MLSPCF_L1B_RAD_END, MLSPCF_L1B_RAD_START
  use OUTPUT_M, only: OUTPUTOPTIONS, &
    & INVALIDPRUNIT, STDOUTPRUNIT, MSGLOGPRUNIT, BOTHPRUNIT, &
    & OUTPUT
  use PRINTIT_M, only: DEFAULTLOGUNIT, GET_CONFIG, STDOUTLOGUNIT

  implicit none
  public
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains initial or permanent settings. Values
  ! are chosen according to what is most suitable for the environment.
  ! For example
  ! certain settings may be appropriate during development but not
  ! for production use, i.e. sips. Therefore, this is a convenient place
  ! to hold everything that needs to be changed before delivery.
  
  ! See also MLSL2PCF, L2ParInfo.parallel, lib/toggles.switches

  ! --------------------------------------------------------------------------
  ! The following should be adjusted before delivery to sips

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! The following should be TRUE if run with level 1 as a single PGE
  ! sharing a single PCF
  ! (in which case we need to move some of the "mobile" PCF ids)
  logical :: SHAREDPCF                               =  .false. 

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! Update these lines before delivery to sips
  ! id to print out in response to "--version" command-line option       
  character(LEN=*), dimension(2), parameter :: CURRENT_VERSION_ID = (/ &
    & 'v4.00 swdev team               ' , & 
    & 'See license terms for copyright'/)
     
  ! Set the following to MSGLOGPRUNIT before delivering to sips;
  ! (its possible values and their effects on normal output:
  ! INVALIDPRUNIT  on MUTE, no output
  ! STDOUTPRUNIT   sent to stdout (via print *, '...')
  ! MSGLOGPRUNIT   sent to Log file (via MLSMessage)
  ! BOTHPRUNIT     both stdout and Log file
  ! > 0            Fortran 'unit=OUTPUT_PRINT_UNIT')
  integer            :: OUTPUT_PRINT_UNIT             = MSGLOGPRUNIT ! -2

  ! Set the following to MLSMSG_Error before delivering to sips;
  ! when set higher, it allows program keep going despite errors
  ! when set lower, the program would quit even on warnings
  integer, parameter :: QUIT_ERROR_THRESHOLD          = MLSMSG_Error

  ! Set the following to 2 before delivering to sips;
  ! If 0, you won't be able to distinguish normal termination
  ! from some abnormal ones (e.g. in parser) (a bad thing)
  ! if 2, status will be 2 only if run complete         
  ! and without error (a good thing)
  integer, parameter :: NORMAL_EXIT_STATUS            = 2

  ! ---------------------------------------------------------------
  ! None of the following need to be changed before delivery to sips
  
  ! Assume hdf files w/o explicit hdfVersion field are this 
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc.        
  integer            :: DEFAULT_HDFVERSION_WRITE      = HDFVERSION_5
  ! Set to WILDCARDHDFVERSION if you wish to autodetect such files  
  ! on input
  integer            :: DEFAULT_HDFVERSION_READ       = WILDCARDHDFVERSION  
  integer            :: LEVEL1_HDFVERSION             = WILDCARDHDFVERSION  

  ! The following is FALSE only for runs that don't need orbit/attitude info
  logical            :: NEED_L1BFILES                 = .true.
  ! The following is FALSE only for runs that use non-Aura satellite data
  logical            :: AURA_L1BFILES                 = .true.
  ! Set if run must not create file, instead just append to it
  logical            :: PATCH                         = .false. 
  ! Whether to restart printing identical warnings at each new phase
  logical            :: RESTARTWARNINGS               = .true.
  ! What units to use in summarizing timings at end of run
  integer            :: SECTIONTIMINGUNITS            = L_SECONDS
  ! Whether to skip doing the direct writes--quicker when snooping
  logical            :: SKIPDIRECTWRITES              = .false.    
  logical            :: SKIPDIRECTWRITESORIGINAL      = .false.    
  ! Whether to skip doing the retrieval--a pre-flight checkout of paths, etc.
  logical            :: SKIPRETRIEVAL                 = .false.        
  ! A holding place for the above, allowing us to skip for some phases only
  logical            :: SKIPRETRIEVALORIGINAL         = .false. 
  ! Whether each slave deallocates all its arrays, pointers, etc.
  ! Sometimes slaves die or take too long to finish cleaning up
  ! But if system fails to reclaim memory properly, subsequent slaves
  ! may not find enough available and therefore crash
  ! FALSE means let operating system do it automatically
  logical            :: SLAVESDOOWNCLEANUP            = .true.
  ! In case special dumps are to go to a special dumpfile
  character(len=255) :: SPECIALDUMPFILE               = ' '
  ! What to fill state, outputSD with if skipping retrieval
  real               :: STATEFILLEDBYSKIPPEDRETRIEVALS = 0.
  ! Whether to stop after a certain section: which section it is
  character(len=16)  :: STOPAFTERSECTION              = ' ' ! Blank means 'no'

  ! Whether to exit with status 1 no matter what
  logical            :: STOPWITHERROR                 = .false.         
  ! Whether to do only a pre-flight checkout of paths
  logical            :: CHECKPATHS                    = .false.         

  logical            :: TOOLKIT                       =  .true. ! SIPS_VERSION
  logical            :: MLSL2DEBUG                    =  .false.
  ! --------------------------------------------------------------------------

  ! This is the type to store runtime Booleans set and used by the l2cf
  integer, parameter :: RTVSTRINGLENGTH               = 1024
  integer, parameter :: RTVARRAYLENGTH                = 128
  
  character(len=2048) :: command_line ! All the opts
  character(len=2048) :: ORIGINALCMDS ! As set when executed
  ! The following will be used only by MLSL2
  logical :: CHECKL2CF = .false.   ! Just check the l2cf and quit
  logical :: CHECKLEAK = .false.   ! Check parse tree for potential memory leaks
  logical :: COUNTCHUNKS = .false. ! Just count the chunks and quit
  integer :: DO_DUMP = 0           ! Dump declaration table if > 0
  integer :: DUMP_TREE = -1        ! Dump tree after parsing
  ! Wouldn't it be better to use get_lun at the moment we open the l2cf?
  integer, parameter :: L2CF_UNIT = 20  ! Unit # if L2CF is opened by Fortran
  integer :: L2CFNODE        = 0        ! Line #, Col # of L2CF being executed
  integer :: NUMSWITCHES
  integer :: RECL            = 20000    ! Record length for l2cf (but see --recl opt)
  character(len=128) :: SECTIONSTOSKIP = ''
  logical :: SECTIONTIMES    = .false.  ! Show times in each section
  logical :: TOTALTIMES      = .false.  ! Show total times from start
  logical :: SHOWDEFAULTS    = .false.  ! Just print default opts and quit
  integer :: SLAVEMAF        = 0        ! Slave MAF for fwmParallel mode
  logical :: TIMING          = .false.  ! -T option is set
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! The following list of public procedures is for convenience only
  public :: DUMPMACROS
  public :: MLSMESSAGE
  public :: PROCESSOPTIONS
  public :: REMOVERUNTIMEBOOLEAN
  public :: RESTOREDEFAULTS
  integer, private :: i ! For loop constructor below

  type :: runTimeValues_T
    character(len=1)                   :: sep = achar(0)
    ! Two arrays bound as a {keys=>values} hash
    character(len=RTVSTRINGLENGTH)     :: lkeys       = &
      & 'true' // achar(0) // 'false' // achar(0) // 'count' 
    character(len=RTVSTRINGLENGTH)     :: lvalues     = &
      & 'true' // achar(0) // 'false' // achar(0) // 'count' 
  end type runTimeValues_T
    
  type(runTimeValues_T), save :: runTimeValues

  type :: dbItem_T
    ! This stores a generic level 2 datatype indexed by the database indx
    ! and with the type specified by dbType
    integer                            :: indx
    integer                            :: dbType  ! One of DB_* types below
  end type dbitem_T
  ! --------------------------------------------------------------------------
  !   For the various databases/data types
  !   We do not support mixing data types
  !   but for certain abstract operations, e.g. adding an item to
  !   a database, the use of an include file is a workaround
  !   For other cases, e.g. in an error message, we may wish to 
  !   communicate a database entry as an index, datatype pair 
  !   and have the calling procedure decide what to do about it
  !   (Here's how we created the list of types that are pasted in
  !    the table below:
  !    grep 'type (' tree_walker.f90 | \
  !     sed -n 's/type (//; s/_T),.*//p' tree_walker.f90

  integer, parameter :: DB_Direct             = 0
  integer, parameter :: DB_File               = DB_Direct             + 1
  integer, parameter :: DB_Chunk              = DB_File               + 1
  integer, parameter :: DB_DirectData         = DB_Chunk              + 1
  integer, parameter :: DB_FGrid              = DB_DirectData         + 1
  integer, parameter :: DB_ForwardModelConfig = DB_FGrid              + 1
  integer, parameter :: DB_GriddedData        = DB_ForwardModelConfig + 1
  integer, parameter :: DB_Hessian            = DB_GriddedData        + 1
  integer, parameter :: DB_HGrid              = DB_Hessian            + 1
  integer, parameter :: DB_L2AUXData          = DB_HGrid              + 1
  integer, parameter :: DB_L2GPData           = DB_L2AUXData          + 1
  integer, parameter :: DB_Matrix             = DB_L2GPData           + 1
  integer, parameter :: DB_Vector             = DB_Matrix             + 1
  integer, parameter :: DB_QuantityTemplate   = DB_Vector             + 1
  integer, parameter :: DB_VectorTemplate     = DB_QuantityTemplate   + 1
  
  logical, private, parameter :: countEmpty = .true. ! Except where overriden locally
  
!=============================================================================
contains 
  ! -------------------------------------------------  MLSMessage  -----
  ! Process the level 2-savvy MLSMessage call
  ! Optionally give l2cf line #, dump an item, and so on before summoning
  ! regular MLSMessage
  !
  ! For some cases we skip calling MLSMessage, either calling output instead
  ! or remaining mute
  !
  ! In certain other cases we must repeat printing the message via the
  ! output module's commands
  subroutine MLSMessage ( SEVERITY, MODULENAMEIN, MESSAGE, &
    & ADVANCE, MLSFILE, STATUS, ITEM )
    use LEXER_CORE, only: GET_WHERE
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use TOGGLES, only: SWITCHES
    use TREE, only: WHERE
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    type(MLSFile_T), intent(in), optional :: MLSFile
    integer, intent(out), optional :: status ! 0 if msg printed, 1 if suppressed
    type(dbItem_T), intent(in), optional  :: item
    ! Internal variables
    integer :: logFileUnit  ! Where we will log
    integer :: LogThreshold ! Severity at which we will log
    logical :: mustRepeat   ! if we will repeat via output
    character(len=256) :: myMessage
    integer :: myStatus
    logical :: outputInstead   ! if we will call output instead
    ! Executable
    call get_config( logFileUnit=logFileUnit )
    myMessage = Message
    mustRepeat = ( logFileUnit == DEFAULTLOGUNIT .and. &
      & OutputOptions%PrUnit == STDOUTPRUNIT )
    LogThreshold = SwitchDetail( switches, 'log' )
    ! Treat log as if it were 'log6'; i.e. output instead every severity
    if ( LogThreshold == 0 ) LogThreshold = 6
    ! Should we call output instead?
    outputInstead = .false.
    if ( LogThreshold > -1 .and. LogThreshold < 6 ) then
      outputInstead = ( severity < LogThreshold .and. &
        & mustRepeat )
    endif
    if ( MLSL2DEBUG ) then
      print *, 'mustRepeat ', mustRepeat
      print *, 'outputInstead ', outputInstead
      print *, 'STDOUTPRUNIT ', STDOUTPRUNIT
      print *, 'OutputOptions%PrUnit ', OutputOptions%PrUnit
      print *, 'DEFAULTLOGUNIT ', DEFAULTLOGUNIT
      print *, 'Config%logFileUnit ', logFileUnit
    endif
    ! For severity "info" just do the call and return
    ! For severity "warn" do the call and 
    ! check status to see if we need to skip printing
    if ( severity == MLSMSG_INFO ) then
      call SayIt ( Message )
      return
    elseif ( severity == MLSMSG_Warning ) then
      ! This trickery is to determine whether this warning would be suppressed
      call SAYMESSAGE ( MLSMSG_TestWarning, ModuleNameIn, MyMessage, &
        & Advance, MLSFile, myStatus )
      if ( present(status) ) status = myStatus
      if ( myStatus /= 0 ) return
    endif
    ! Do we have an l2cf node we were processing?
    if ( L2CFNode /= 0 ) then
      call get_where ( where(L2CFNode), myMessage, &
        & before='***** At ', after=': ' // Message )
    else
      myMessage = '(no tree node): ' // message
    end if
    if ( present(item) ) then
      call output ( '(Not ready to dump arbitrary datatypes yet) ', advance='yes' )
      select case ( item%dbType )
      case ( DB_Direct             )
      case ( DB_File               )
      case ( DB_Chunk              )
      case ( DB_DirectData         )
      case ( DB_FGrid              )
      case ( DB_ForwardModelConfig )
      case ( DB_GriddedData        )
      case ( DB_Hessian            )
      case ( DB_HGrid              )
      case ( DB_L2AUXData          )
      case ( DB_L2GPData           )
      case ( DB_Matrix             )
      case ( DB_Vector             )
      case ( DB_QuantityTemplate   )
      case ( DB_VectorTemplate     )
      case default
      end select
    endif
    call SayIt ( myMessage )
  contains
    subroutine SayIt ( It )
      ! Say it with MLSMessage
      ! and possibly rpeat it with output
      character (len=*), intent(in) :: It
      if ( .not. outputInstead ) call SAYMESSAGE ( severity, ModuleNameIn, It, &
        & Advance, MLSFile, status )
      if ( mustRepeat ) call output( trim(It), advance='yes' )
    end subroutine SayIt
  end subroutine MLSMessage

  ! ---------------------------------------------  ProcessOptions  -----
  ! Process the command line options; either from
  ! (1) the command line directly, by calls to getNextArg; or
  ! (2) parsing cmdline, if supplied as an arg
  ! The return value will be the filename (if any)
  function ProcessOptions ( cmdLine ) result ( fileName )
    use ALLOCATE_DEALLOCATE, only: TRACKALLOCATES, &
      & CLEARONALLOCATE
    use IO_STUFF, only: GET_LUN
    use L2PARINFO, only: PARALLEL, INITPARALLEL, ACCUMULATESLAVEARGUMENTS, &
      & SNIPLASTSLAVEARGUMENT
    use LEXER_M, only: CAPIDENTIFIERS
    use MACHINE, only: GETARG, HP, IO_ERROR, NEVERCRASH
    use MATRIXMODULE_0, only: CHECKBLOCKS, SUBBLOCKLENGTH
    use MLSCOMMON, only: FILENAMELEN
    use MLSMESSAGEMODULE, only: SetConfig
    use MLSSTRINGLISTS, only: CATLISTS, &
      & GETSTRINGELEMENT, GETUNIQUELIST, &
      & NUMSTRINGELEMENTS, REMOVESWITCHFROMLIST, &
      & SORTLIST, STRINGELEMENT, SWITCHDETAIL, UNQUOTE
    use MLSSTRINGS, only: LOWERCASE, READINTSFROMCHARS
    use PCFHDR, only: GLOBALATTRIBUTES
    use PRINTIT_M, only: SET_CONFIG
    use SET_TOGGLES_M, only: SET_TOGGLES
    use SNOOPMLSL2, only: SNOOPINGACTIVE, SNOOPNAME
    use STRING_TABLE, only: ADD_INCLUDE_DIRECTORY, DO_LISTING
    use TIME_M, only: TIME_CONFIG
    use TOGGLES, only: SWITCHES
    ! Args
    character(len=*), intent(in), optional :: cmdLine
    character(len=FileNameLen)             :: fileName
    ! Internal variables
    character(len=1) :: arg_rhs      ! 'n' part of 'arg=n'
    character(len=16) :: aSwitch
    integer :: DEGREE                ! index affecting degree of option
    logical, parameter :: DEEBUG = .false.
    logical :: EXIST
    logical :: Field_Is_Include      ! Field is include file path, not L2CF name
    integer, dimension(100)           :: iarray
    integer :: J
    character(len=2048) :: LINE      ! Into which is read the command args
    integer :: N
    logical :: OPENED
    character(len=2) :: quotes
    integer :: RECL = 20000          ! Record length for l2cf (but see --recl opt)
    character(len=len(switches)) :: removeSwitches = ''
    integer :: STATUS
    logical :: SWITCH                ! "First letter after -- was not n"
    character(len=len(switches)) :: tempSwitches
    integer :: V
    character(len=2048) :: WORD      ! Some text
    ! Executable
    quotes = char(34) // char(39)   ! {'"}
    filename = 'help' ! This means abnormal options--should dump help mesg
    if ( present( cmdline ) .and. DEEBUG ) then
      print *, 'cmdline: ', trim(cmdline)
    endif
  ! Before looking at command-line options, TOOLKIT is set to SIPS_VERSION
  ! So here's a good place to put any SIPS-specific settings overriding defaults
  if ( .true.  ) then ! SIPS_VERSION
    ! SIPS_VERSION
    parallel%maxFailuresPerMachine = 2
    parallel%maxFailuresPerChunk = 1
    removeSwitches='slv' ! Since slave output already saved to separate files
    ! switches='red'  ! No longer a good idea
    DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
    MLSMessageConfig%limitWarnings = 4 ! 50 ! Why print all that stuff?
    time_config%use_wall_clock = .true. ! SIPS_VERSION
  else
    ! SCF_VERSION
    switches='0sl'
    time_config%use_wall_clock = .false. ! SIPS_VERSION
  end if
    i = 1+hp
    do ! Process Lahey/Fujitsu run-time options; they begin with "-Wl,"
      call getNextArg ( i, line )
      if ( line(1:4) /= '-Wl,' ) then
        call SnipLastSlaveArgument ! Don't want slaves to see this
        exit
      endif
      call AccumulateSlaveArguments(line) ! pass them to slave processes
      i = i + 1
    end do
    ! Now process the other options
    command_line = ' '
    field_is_include = .false.
    cmds:    do
      call getNextArg ( i, line )
      if ( DEEBUG ) print *, i, trim(line)
      if ( len_trim( line ) < 1 ) then
        exit
      endif
      if ( line(1:2) == '--' ) then       ! "word" options
        line = lowercase(line)
        n = 0
        switch = .true.
        if ( line(3:3) == 'n' ) then
          switch = .false.
          n = 1
        end if
        if ( line(3+n:8+n) == 'check ' ) then
          checkl2cf = switch
        ! Using lowercase so either --checkPaths or --checkpaths work
        ! Perhaps we should do this for all multiletter options
        ! (single-letter options are case-sensitive)
        else if ( line(3+n:8+n) == 'checkp' ) then
          checkPaths = switch
        else if ( line(3+n:7) == 'chunk' ) then
          i = i + 1
          call getNextArg ( i, line )
          parallel%chunkRange = line
        else if ( line(3+n:18+n) == 'clearonallocate ' ) then
          clearonallocate = switch
        else if ( line(3+n:7+n) == 'ckbk ' ) then
          checkBlocks = switch
        else if ( line(3+n:14+n) == 'countchunks ' ) then
          countChunks = switch
        else if ( line(3+n:7+n) == 'crash' ) then
          MLSMessageConfig%crashOnAnyError = switch
          neverCrash = .not. switch
        else if ( line(3+n:8+n) == 'defaul' ) then
          showDefaults = switch
        else if ( line(3+n:7+n) == 'delay' ) then
          if ( line(8+n:) /= ' ' ) then
            line(:7+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) parallel%Delay
          if ( status /= 0 ) then
            call io_error ( "After --delay option", status, line )
            stop
          end if
        else if ( line(3+n:7+n) == 'dump ' ) then
          i = i + 1
          call getNextArg ( i, line )
          specialDumpFile = trim(line)
        else if ( line(3+n:14+n) == 'fwmparallel ' ) then
          parallel%fwmParallel = .true.
        else if ( line(3+n:7+n) == 'host ' ) then
          i = i + 1
          call getNextArg ( i, line )
          GlobalAttributes%hostName = trim(line)
        else if ( line(3+n:9+n) == 'idents ' ) then
          i = i + 1
          call getNextArg ( i, line )
        ! else if ( line(3+n:6+n) == 'kit ' ) then
          ! MLSMessageConfig%useToolkit = switch
        else if ( line(3+n:6+n) == 'l1b=' ) then
          arg_rhs = line(7+n:7+n)
          select case(arg_rhs)
          case('4')
            LEVEL1_HDFVERSION = HDFVERSION_4
          case('5')
            LEVEL1_HDFVERSION = HDFVERSION_5
          case default
            LEVEL1_HDFVERSION = WILDCARDHDFVERSION
          end select
        else if ( line(3+n:8+n) == 'l2gpr=' ) then
          arg_rhs = line(9+n:9+n)
          select case(arg_rhs)
          case('4')
            DEFAULT_HDFVERSION_READ = HDFVERSION_4
          case('5')
            DEFAULT_HDFVERSION_READ = HDFVERSION_5
          case default
            DEFAULT_HDFVERSION_READ = WILDCARDHDFVERSION
          end select
        else if ( line(3+n:8+n) == 'l2gpw=' ) then
          arg_rhs = line(9+n:9+n)
          select case(arg_rhs)
          case('4')
            DEFAULT_HDFVERSION_WRITE = HDFVERSION_4
          case('5')
            DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
          case default
            DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
          end select
        else if ( line(3+n:5+n) == 'lac' ) then
          ! print *, 'Got the laconic option'
          ! The laconic command-line option trims logged output
          ! resulting from calls to MLSMessageModule which normally prefix
          ! each line with severity and module name; e.g.,
          ! Info (output_m):         ..Exit ReadCompleteHDF5L2PC at 11:37:19.932 
          ! Usage: laconic n where if n is
          ! 0  just abbreviate module names and severities
          ! 1  omit module names and severity unless Warning or worse
          ! 2  omit module names and severity unless Error
          ! 3  omit module names and severity completely
          ! 11 skip printing anything unless Warning or worse
          ! 12 skip printing anything unless Error
          i = i + 1
          if ( .not. switch ) then
            MLSMessageConfig%MaxModuleNameLength     = 10
            MLSMessageConfig%MaxSeverityNameLength   = 8
            cycle
          endif
          call getNextArg ( i, line )
          read ( line, *, iostat=status ) degree
          if ( status /= 0 ) then
            call io_error ( "After --lac[onic] option", status, line )
            stop
          end if
          MLSMessageConfig%suppressDebugs = (degree > 0)
          ! MLSMessageConfig%AbbreviateModSevNames = (degree == 0)
          if ( degree == 0 ) then
            MLSMessageConfig%MaxModuleNameLength     = 3
            MLSMessageConfig%MaxSeverityNameLength   = 1
          elseif ( degree == 3 ) then
            MLSMessageConfig%MaxModuleNameLength     = 0
            MLSMessageConfig%MaxSeverityNameLength   = 0
          endif
          MLSMessageConfig%skipModuleNamesThr = mod(degree, 10) + 2
          MLSMessageConfig%skipSeverityThr = mod(degree, 10) + 2
          MLSMessageConfig%skipMessageThr = degree - 10 + 2
          ! print *, ' Processing lac option: degree ', degree
          if ( degree > 9 ) then
            removeSwitches = catLists(trim(removeSwitches), 'log' )
            outputOptions%prunit = INVALIDPRUNIT
          endif
        else if ( line(3+n:7+n) == 'leak' ) then
          checkLeak = .true.
        else if ( line(3+n:9+n) == 'master ' ) then
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%master = .true.
          i = i + 1
          call getNextArg ( i, line )
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%slaveFilename = trim ( line )
          call InitParallel ( 0, 0 )
          word = '--slave'
          write ( word(len_trim(word)+1:), * ) parallel%myTid
          call AccumulateSlaveArguments(word)
        else if ( line(3+n:21+n) == 'maxfailuresperchunk' ) then
          if ( line(22+n:) /= ' ' ) then
            line(:21+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) parallel%maxFailuresPerChunk
          if ( status /= 0 ) then
            call io_error ( "After --maxFailuresPerChunk option", status, line )
            stop
          end if
        else if ( line(3+n:23+n) == 'maxfailurespermachine' ) then
          if ( line(24+n:) /= ' ' ) then
            line(:23+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) parallel%maxFailuresPerMachine
          if ( status /= 0 ) then
            call io_error ( "After --maxFailuresPerMachine option", status, line )
            stop
          end if
        else if ( line(3+n:10+n) == 'memtrack' ) then
          v = 1
          if ( line(11+n:) /= ' ' ) then
            read ( line(11+n:), *, iostat=status ) v
            if ( status /= 0 ) then
              call io_error ( "After --memtrack option", status, line )
              stop
            end if
          else
            call getNextArg ( i+1, line )
            read ( line, *, iostat=status ) j
            if ( status == 0 ) then
              i = i + 1
              v = j
            end if
          end if
          if ( switch ) then
            trackAllocates = v
          else
            trackAllocates = 0
          end if
        else if ( line(3+n:9+n) == 'msgconf' ) then
          if ( line(10+n:10+n) /= '=' ) then
            i = i + 1
            call getNextArg ( i, line(11:) )
          end if
          call setConfig ( [ line(11:) ] )
        else if ( line(3+n:4+n) == 'oa' ) then
          NEED_L1BFILES = switch
        else if ( line(3+n:6+n) == 'aura' ) then
          AURA_L1BFILES = switch
        else if ( line(3+n:8+n) == 'patch ' ) then
          patch = switch
        else if ( line(3+n:5+n) == 'pge ' ) then
          call SnipLastSlaveArgument ! Don't want slaves to see this
          i = i + 1
          call getNextArg ( i, line )
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%pgeName = trim(line)
        else if ( line(3+n:6+n) == 'recl' ) then
          if ( line(7+n:) /= ' ' ) then
            line(:6+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) recl
          if ( status /= 0 ) then
            call io_error ( "After --recl option", status, line )
            stop
          end if
        else if ( line(3+n:6+n) == 'shar' ) then
          SHAREDPCF = switch
        else if ( line(3+n:9+n)  == 'skipdir' ) then
          SKIPDIRECTWRITES = switch
        else if ( line(3+n:10+n) == 'skipretr' ) then
          SKIPRETRIEVAL = switch
        else if ( line(3+n:9+n) == 'skipsec' ) then
          i = i + 1
          call getNextArg ( i, line )
          sectionsToSkip = lowercase(line)
        else if ( line(3+n:10+n) == 'slavemaf' ) then
          if ( line(11+n:) /= ' ' ) then
            line(:10+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) slaveMAF
          if ( status /= 0 ) then
            call io_error ( "After --slaveMAF option", status, line )
            stop
          end if
        else if ( line(3+n:7+n) == 'slave' ) then
          parallel%slave = .true.
          if ( line(8+n:) /= ' ' .and. .false.) then
            line(:7+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) parallel%masterTid
          if ( status /= 0 ) then
            call io_error ( "After --slave option", status, line )
            stop
          end if
          MLSMessageConfig%SendErrMsgToMaster = .true.
          MLSMessageConfig%masterTID = parallel%masterTid
          if ( DEEBUG ) print *, 'masterTid: ', parallel%masterTid
        else if ( line(3+n:8+n) == 'snoop ' ) then
          snoopingActive = .true.
        else if ( line(3+n:12+n) == 'snoopname' ) then
          i = i + 1
          call getNextArg ( i, line )
          snoopName = line
        else if ( line(3+n:7+n) == 'state' ) then
          if ( line(8+n:) /= ' ' ) then
            line(:7+n) = ' '
          else
            i = i + 1
            call getNextArg ( i, line )
          end if
          read ( line, *, iostat=status ) stateFilledBySkippedRetrievals
          if ( status /= 0 ) then
            call io_error ( "After --state option", status, line )
            stop
          end if
        else if ( line(3+n:9+n) == 'stdout ' ) then
          if ( .not. switch ) then
            OUTPUT_PRINT_UNIT = INVALIDPRUNIT
            return
          endif
          i = i + 1
          call getNextArg ( i, line )
          select case ( lowercase(line) )
          case ( 'log' )
            OUTPUT_PRINT_UNIT = MSGLOGPRUNIT
          case ( 'out' )
            OUTPUT_PRINT_UNIT = STDOUTPRUNIT
          case ( 'both' )
            OUTPUT_PRINT_UNIT = BOTHPRUNIT
            ! MLSMessageConfig%logFileUnit = STDOUTLOGUNIT
          case ( 'unbuffered' )
            OutputOptions%buffered = .false.
          case default
            OutputOptions%name = trim(line)
            OutputOptions%buffered = .false.
            ! outputOptions%debugUnit = 32
            ! Make certain prUnit won't be l2cf_unit
            open( unit=l2cf_unit, status='unknown' )
            call get_lun ( OutputOptions%prUnit, msg=.false. )
            close( unit=l2cf_unit )
            inquire( unit=OutputOptions%prUnit, exist=exist, opened=opened )
            OUTPUT_PRINT_UNIT = OutputOptions%prUnit
          end select
        else if ( line(3+n:12+n) ==  'stopafter ' ) then
          i = i + 1
          call getNextArg ( i, stopAfterSection )
        else if ( line(3+n:12+n) ==  'stopwither' ) then
          stopWithError = switch
        else if ( line(3+n:11+n) == 'subblock ' ) then
          i = i + 1
          call getNextArg ( i, line )
          read ( line, *, iostat=status ) subBlockLength
          if ( status /= 0 ) then
            call io_error ( "After --subblock option", status, line )
            stop
          end if
        else if ( line(3+n:9+n) == 'submit ' ) then
          i = i + 1
          call getNextArg ( i, line )
          parallel%submit = trim ( line )
        else if ( line(3+n:5+n) == 'tk ' ) then
          toolkit = switch
          call set_config ( useToolkit = switch )
        else if ( line(3+n:9+n) == 'verbose' ) then
          switches = catLists( trim(switches), &
            & 'l2q,glob,mas,bool,opt1,log,pro1,time,apr,phase' )
        else if ( line(3+n:10+n) == 'version ' ) then
          do j=1, size(current_version_id)
            print *, current_version_id(j)
          end do
          stop
        else if ( line(3+n:7+n) == 'wall ' ) then
          time_config%use_wall_clock = switch
        else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
          i = i + 1
          call getNextArg ( i, line )
          exit
        else
          print *, 'unrecognized option ', trim(line), ' ignored.'
          ! call option_usage
          filename = 'help'
          return ! will dump help mesg
        end if
      else if ( line(1:1) == '-' ) then   ! "letter" options
        j = 1
jloop:do while ( j < len_trim(line) )
          j = j + 1
          select case ( line(j:j) )
          case ( ' ' )
            exit
          case ( 'A' )
            dump_tree = 0
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              dump_tree = ichar(line(j+1:j+1)) - ichar('0')
              j = j + 1
            end if
          case ( 'a', 'c', 'f', 'g', 'l', 'p', 't' )
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              call set_toggles ( line(j:j+1) )
              j = j + 1
            else
              call set_toggles ( line(j:j) )
            end if
          case ( 'D' ) ! This turns debugging on for some modules
            MLSNamesAreDebug = catLists(trim(MLSNamesAreDebug), line(j+1:))
            exit ! Took the rest of the string, so there can't be more options
          case ( 'd' ); do_dump = 1
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              do_dump = ichar(line(j+1:j+1)) - ichar('0')
              j = j + 1
            end if
          case ( 'h', 'H', '?' )     ! Describe command line usage
            ! call option_usage
            filename = 'help'
            return ! will dump help mesg
          case ( 'I' ); field_is_include = .true. ! Next field is include path
          case ( 'K' ); capIdentifiers = .true.
          case ( 'k' ); capIdentifiers = .false.
          case ( 'M' ); outputOptions%prunit = MSGLOGPRUNIT ! -2
          case ( 'm' ); outputOptions%prunit = STDOUTPRUNIT ! -1
          case ( 'R' ) ! This does the opposite of what S does
            removeSwitches = catLists(trim(removeSwitches), line(j+1:))
            exit ! Took the rest of the string, so there can't be more options
          case ( 'S' )
            switches = catLists(trim(switches), line(j+1:))
            exit ! Took the rest of the string, so there can't be more options
          case ( 'T' )
            timing = .true.
            do
              if ( j >= len(line) ) exit
              if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
                if( line(j+1:j+1) /= '0' ) &
                  & switches = catLists(trim(switches), 'time')
                sectiontimes = &
                  & switchDetail(switches, 'time') > -1 &
                  & .and. &
                  & switchDetail(switches, 'time') /= 1
                totaltimes = sectiontimes .and. switchDetail(switches, 'time') /= 2
              else if ( lowercase(line(j+1:j+1)) == 's' ) then
                sectionTimingUnits = l_seconds
              else if ( lowercase(line(j+1:j+1)) == 'm' ) then
                sectionTimingUnits = l_minutes
              else if ( lowercase(line(j+1:j+1)) == 'h' ) then
                sectionTimingUnits = l_hours
              end if
              j = j + 1
            end do
          case ( 'V' ) ! This sets verbose to TRUE some modules
            MLSNamesAreVerbose = catLists(trim(MLSNamesAreVerbose), line(j+1:))
            exit ! Took the rest of the string, so there can't be more options
          case ( 'v' ); do_listing = .true.
          case ( 'w' )
            MLSMessageConfig%limitWarnings = 1
            if ( line(j+1:j+1) == 'p' ) then
              RESTARTWARNINGS = .false.
              j = j + 1
            end if
            if ( j < len_trim(line) ) then
              call readIntsFromChars(line(j+1:), MLSMessageConfig%limitWarnings)
              j = j + len_trim(line(j+1:))
            end if
          case default
            print *, 'Unrecognized option -', line(j:j), ' ignored.'
            ! call option_usage
            filename = 'help'
            return ! will dump help mesg
          end select
          ! i = i + 1
        end do jLoop
      else if ( field_is_include ) then
        call add_include_directory ( line )
        field_is_include = .false.
      else
        filename = line
        exit ! This must be the l2cf filename
      end if
      i = i + 1
    enddo cmds
    ! Did we somehow miss the filename among the args?
    if ( filename == 'help' ) then
      i = i + 1
      call getNextArg ( i, line )
      filename = line
    endif
  ! Are any switches inappropriate for master or for slave?
    if ( parallel%master ) &
      & removeSwitches = catLists( trim(removeSwitches), 'bool,walk' )
    if ( parallel%slave ) &
      & removeSwitches = catLists( trim(removeSwitches), 'chu,l2q,mas,slv' )
    ! Remove any quote marks from RemoveSwitches array
    tempSwitches = unquote(removeSwitches, quotes=quotes, options='-p')
    call GetUniqueList( tempSwitches, removeSwitches, numSwitches, &
          & ignoreLeadingSpaces=.true., options='-eSL' )
    ! Remove any quote marks from switches array
    tempSwitches = unquote(switches, quotes=quotes, options='-p')
    ! Now we want to keep only the swich with the highest details level
    call sortList( tempSwitches, iarray, ',', switches )
    tempSwitches = switches
    call GetUniqueList( tempSwitches, Switches, numSwitches, &
          & ignoreLeadingSpaces=.true., options='-eSL' )
    ! Remove any switches embedded in the removeSwitches option 'R'
    ! call outputNamedValue( 'starting List', trim(switches) )
    ! call outputNamedValue( 'to remove', trim(removeSwitches) )
    do i=1, NumStringElements(removeSwitches, countEmpty=.true.)
      call GetStringElement(trim(removeSwitches), aSwitch, i, countEmpty=.true.)
      call RemoveSwitchFromList(switches, tempSwitches, trim(aSwitch))
      switches = tempSwitches
    end do
    ! call outputNamedValue( 'result List', trim(switches) )

    ! If we like, we could move these next few statements to a standalone
    ! subroutine named something like processSwitches
    parallel%verbosity = switchDetail(switches, 'mas') + 1
    if ( switchDetail(switches, 'walk') > -1 ) &
      & MLSMSG_Severity_to_walkback = MLSMSG_Warning
    if ( outputOptions%prunit /= INVALIDPRUNIT ) &
      & outputOptions%prunit = OUTPUT_PRINT_UNIT
    ! print *, 'Ended processing options'
    contains
    subroutine getNextArg( i, line )
      ! Args
      integer, intent(in)           :: i
      character(len=*), intent(out) :: line
      logical, parameter :: countEmpty = .false.
      ! Executable
      if ( present( cmdline ) ) then
        line = StringElement ( cmdline, i, countEmpty, inseparator=' ' )
      else
        call getArg ( i, line )
      endif
      command_line = trim(command_line) // ' ' // trim(line)
      call AccumulateSlaveArguments(line) ! pass them to slave processes
    end subroutine getNextArg
  end function ProcessOptions
  
  ! --------------------------------------------  RestoreDefaults  -----
  ! Restore the options to their default values
  ! Now some things it makes no sense to overwrite, so it makes
  ! no sense to restore them either; e.g., CHECKPATHS, parallel, etc.
  subroutine restoreDefaults
  use MLSMESSAGEMODULE, only: RESTORECONFIG
  use TOGGLES, only: INIT_TOGGLE
    OUTPUT_PRINT_UNIT             = -2
    DEFAULT_HDFVERSION_WRITE      = HDFVERSION_5
    DEFAULT_HDFVERSION_READ       = WILDCARDHDFVERSION
    LEVEL1_HDFVERSION             = WILDCARDHDFVERSION
    NEED_L1BFILES                 = .true.
    AURA_L1BFILES                 = .true.
    PATCH                         = .false. 
    RESTARTWARNINGS               = .true.
    SECTIONTIMINGUNITS            = L_SECONDS
    SKIPDIRECTWRITES              = .false.    
    SKIPDIRECTWRITESORIGINAL      = .false.    
    SKIPRETRIEVAL                 = .false.        
    SKIPRETRIEVALORIGINAL         = .false. 
    SLAVESDOOWNCLEANUP            = .false.
    SPECIALDUMPFILE               = ' '
    STATEFILLEDBYSKIPPEDRETRIEVALS = 0.
    STOPAFTERSECTION              = ' ' ! Blank means 
    STOPWITHERROR                 = .false.         
    CHECKPATHS                    = .false.         
    TOOLKIT                       =  .true. ! SIPS_VERSION
    call init_toggle
    call restoreConfig
  end subroutine restoreDefaults

  ! -------------------------------------------------  DumpMacros  -----
  ! Dump the runtime macros
  subroutine DUMPMACROS
  use DUMP_0, only: DUMP
  call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time macros', separator=runTimeValues%sep )
  end subroutine DumpMacros

  ! ---------------------------------------  RemoveRuntimeBoolean  -----
  ! Dump the runtime macros
  subroutine REMOVERUNTIMEBOOLEAN ( NAME )
  use MLSSTRINGLISTS, only: GETHASHELEMENT, &
    & REMOVEHASHARRAY, REMOVEHASHELEMENT
    ! Dummy args
    character(len=*), intent(in) :: NAME
    ! Internal variables
    character (len=16)                            :: keyString
    character (len=8)                             :: nCh
    ! Executable
    ! First: is the name an array-valued r/t Boolean?
    keyString = trim(name) // 'n'
    call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, keyString, nCh, &
      & countEmpty, runTimeValues%sep )
    if ( nCh /= runTimeValues%sep ) then
      ! Yes, it is array-valued; remove it
      call RemoveHashArray( runTimeValues%lkeys, runTimeValues%lvalues, name, &
        & countEmpty, runTimeValues%sep )
    else
      ! No, so is it scalar-valued?
      call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, name, nCh, &
        & countEmpty, runTimeValues%sep )
      if ( nCh == runTimeValues%sep ) return
      ! Yes, it is scalar-valued; remove it
      call RemoveHashElement( runTimeValues%lkeys, runTimeValues%lvalues, name, &
        & countEmpty, runTimeValues%sep )
    endif
  end subroutine removeRuntimeBoolean

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE MLSL2Options
!=============================================================================

!
! $Log$
! Revision 2.79  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.78  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.77  2013/12/12 02:10:07  vsnyder
! Change type of do_dump from logical to integer
!
! Revision 2.76  2013/12/05 23:52:21  vsnyder
! Process # in -A# option correctly
!
! Revision 2.75  2013/11/26 22:40:51  vsnyder
! Change -A to -A[n] with n>0 meaning dump the entire tree, including the
! type-checking stuff, and n==0 or absent meaning dump only the parser output.
!
! Revision 2.74  2013/11/20 01:00:45  pwagner
! slaves were dumping chunkdivide data for all chunks; fixed
!
! Revision 2.73  2013/11/04 22:56:02  pwagner
! Added -Vmodules and -Dmodules to turn on verbose and debug module-wide
!
! Revision 2.72  2013/09/24 23:29:45  vsnyder
! Add -I option, Use Where instead of Source_Ref for messages
!
! Revision 2.71  2013/09/14 01:22:02  vsnyder
! Add MSGConf option
!
! Revision 2.70  2013/09/06 20:49:17  pwagner
! Solve another case where we repeated warnings
!
! Revision 2.69  2013/09/04 17:34:03  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.68  2013/08/23 02:52:13  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 2.67  2013/08/17 00:22:14  pwagner
! New cmdline arg relaxes some for non-Aura l1b datasets
!
! Revision 2.66  2013/08/13 23:05:48  pwagner
! Removd some fugitive debug printing
!
! Revision 2.65  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.64  2013/06/19 00:40:34  pwagner
! No longer automatically reduce Details of l2cf-born Dumps
!
! Revision 2.63  2013/06/14 01:26:37  vsnyder
! handle --stdout unbuffered
!
! Revision 2.62  2013/06/12 02:37:36  vsnyder
! Cruft removal
!
! Revision 2.61  2013/05/22 20:21:33  pwagner
! Added removeRuntimeBoolean
!
! Revision 2.60  2013/05/17 00:52:58  pwagner
! runtime value sep now achar(0)
!
! Revision 2.59  2013/02/12 18:14:15  pwagner
! Removed SIPS_VERSION
!
! Revision 2.58  2013/02/04 22:01:02  pwagner
! Added '--verbose' option; '--lac' more so
!
! Revision 2.57  2012/12/04 00:15:49  pwagner
! Removed confisuion-causing OUTSIDEOVERLAPS and its cmdline option
!
! Revision 2.56  2012/08/30 20:54:08  pwagner
! Improved adding, removing switches
!
! Revision 2.55  2012/08/16 17:46:17  pwagner
! Added a level 2-savvy MLSMessage to interpose between level 2 procedures and lib version
!
! Revision 2.54  2012/07/18 00:38:00  pwagner
! Consistent with module parameters for prUnit
!
! Revision 2.53  2012/07/10 15:23:42  pwagner
! Works properly now; api adjusted for GetUniqueList
!
! Revision 2.52  2012/07/02 20:29:32  pwagner
! Improve RestoreDefaults, some housekeeping
!
! Revision 2.51  2012/06/27 18:02:09  pwagner
! May overwrite command line options with options field to phase spec
!
! Revision 2.50  2012/05/08 17:53:51  pwagner
! Converted runtimes to character-valued; added DB types
!
! Revision 2.49  2011/06/29 21:43:23  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.48  2010/04/12 22:20:23  pwagner
! Changed vers. id to conform with v3.30 sips id
!
! Revision 2.47  2009/11/05 00:31:21  pwagner
! Updated version id
!
! Revision 2.46  2009/07/24 23:22:47  pwagner
! Updated version id, copyright statement
!
! Revision 2.45  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.44  2009/04/13 21:00:44  pwagner
! update current version id to pre-v3
!
! Revision 2.43  2008/01/08 00:16:29  pwagner
! Added SHAREDPCF so levels 1 and 2 can use same PCF
!
! Revision 2.42  2007/09/06 22:42:28  pwagner
! Added SLAVESDOOWNCLEANUP
!
! Revision 2.41  2007/06/21 22:35:22  pwagner
! Updated version string to v2.22
!
! Revision 2.40  2007/02/06 23:15:48  pwagner
! CURRENT_VERSION_ID now v2.21
!
! Revision 2.39  2006/11/01 20:37:38  pwagner
! CURRENT_VERSION_ID now v2.20
!
! Revision 2.38  2006/07/27 23:06:17  pwagner
! update CURRENT_VERSION_ID
!
! Revision 2.37  2006/07/21 20:09:56  pwagner
! Can fill state even if skipping retrievals; select what section to stop after
!
! Revision 2.36  2006/06/13 00:16:27  pwagner
! catenating split dgg/dgm files now on by default
!
! Revision 2.35  2006/02/21 19:19:27  pwagner
! New things to create, refer to run time booleans in l2cf
!
! Revision 2.34  2006/02/10 21:13:30  pwagner
! May specify skipRetrivel for particular Phases; dumps may go to special dumpfile
!
! Revision 2.33  2005/07/21 23:40:54  pwagner
! Removed unneeded ILLEGALL1BRADID, MAXNUML1BRADIDS
!
! Revision 2.32  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.31  2005/03/12 00:48:01  pwagner
! Added RESTARTWARNINGS; corrected vsn id
!
! Revision 2.30  2004/12/14 00:04:24  pwagner
! New early stop options added for quicker debugging
!
! Revision 2.29  2004/07/08 22:48:44  pwagner
! Made SIPS_VERSION public
!
! Revision 2.28  2004/04/27 23:49:51  pwagner
! Added SKIPDIRECTWRITES option
!
! Revision 2.27  2004/03/12 00:28:56  pwagner
! At last hdf version at output increased to 5
!
! Revision 2.26  2004/01/23 01:06:39  pwagner
! Added CATENATESPLITS
!
! Revision 2.25  2003/12/05 00:39:35  pwagner
! Added patch option, section timing units
!
! Revision 2.24  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.23  2003/10/09 23:58:34  pwagner
! Updated CURRENT_VERSION_ID to 1.4
!
! Revision 2.22  2003/09/05 23:22:52  pwagner
! Has new SKIPRETRIEVAL option
!
! Revision 2.21  2003/06/09 22:49:32  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.20  2003/05/02 20:53:19  pwagner
! Reordered to make SIPS-dependent section clearer; default_hdfversion at read now wildcard
!
! Revision 2.19  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/10/03 23:00:03  pwagner
! You can set l1b, l2gp hdfversions on command line
!
! Revision 2.17  2002/08/28 22:25:42  pwagner
! Moved LEVEL1_HDFVERSION, ILLEGALL1BRADID, MAXNUML1BRADIDS here from global_settings
!
! Revision 2.16  2002/03/14 23:38:28  pwagner
! Gets HDFVERSION_4 and 5 from MLSFiles module
!
! Revision 2.15  2002/02/12 00:25:25  pwagner
! New current_version_id parameter
!
! Revision 2.14  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.13  2002/01/29 23:49:38  pwagner
! Separate DEFAULT_HDFVERSION_(READ)(WRITE)
!
! Revision 2.12  2002/01/23 21:48:16  pwagner
! Added DEFAULT_HDFVERSION
!
! Revision 2.11  2001/09/28 17:57:47  pwagner
! SIPS_VERSION controls other logical options
!
! Revision 2.10  2001/07/16 23:43:15  pwagner
! With settable NORMAL_EXIT_STATUS
!
! Revision 2.9  2001/05/30 22:56:48  pwagner
! Moved PCFL2CFSAMECASE here from OutputAndClose
!
! Revision 2.8  2001/05/15 23:46:08  pwagner
! Removed 2 settings from MLSL2Opts; now in switches
!
! Revision 2.7  2001/05/11 23:48:23  pwagner
! Changed to not echo globals; added note on SIPS
!
! Revision 2.6  2001/05/09 23:34:13  pwagner
! Added ECHO_GLOBAL_STNGS LOG_TO_STDOUT
!
! Revision 2.5  2001/05/06 20:54:40  pwagner
! Default settings should work for most jpl users
!
! Revision 2.4  2001/05/04 22:54:31  pwagner
! Added TOOLKIT, CREATEMETADATA, PCF_FOR_INPUT
!
! Revision 2.3  2001/04/20 20:41:52  pwagner
! Added QUIT_ERROR_THRESHOLD
!
! Revision 2.2  2001/04/17 20:26:28  pwagner
! Added OUTPUT_PRINT_UNIT
!
! Revision 2.1  2001/04/16 23:53:10  pwagner
! First commit
!
