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
module MLSL2Options              !  Options and Settings for the MLSL2 program
!=============================================================================
  ! Subroutines have also found their way into this module. While very handy
  ! they go against the clean-sheet design which reserved it for
  ! data only. We're probably too far down the path now to retreat.
  ! A similar trespass was lib/MLSCommon which began as a data-only module
  ! but has grown to hold lots of indispensible procedures.

  use Dump_0, only: Dump
  use Dump_1, only: Dump
  use HDF, only: Dfacc_Rdwr
  use HighOutput, only: Banner, OutputNamedValue, StyledOutput
  use Intrinsic, only: L_Ascii, L_HDFEOS, L_Hours, L_Minutes, L_NetCDF4, L_Seconds
  use, Intrinsic :: ISO_Fortran_Env, only: Error_Unit
  use MLSCommon, only: MLSFile_T, MLSNamesAreDebug, MLSNamesAreVerbose
  use MLSFiles, only: WildCardHDFVersion, HDFVersion_4, HDFVersion_5, Dump
  use MLSMessageModule, only: MLSMessageConfig, &
    & MLSMsg_Error, MLSMsg_Info, MLSMsg_TestWarning, &
    & MLSMSG_Severity_To_Quit, MLSMsg_Severity_To_Walkback, MLSMsg_Warning, &
    & Bummer, SayMessage => MLSMessage
  use MLSPCF2, only: MLSPCF_L1b_Rad_End, MLSPCF_L1b_Rad_Start
  use MLSStringLists, only: EvaluateFormula
  use MLSStrings, only: IsComment, IsDigits, LowerCase, &
    & ReadIntsFromChars, Replace, WriteIntsToChars
  use PCFHdr, only: GlobalAttributes, SomeGlobalAttributes, &
    & SomeToGlobalAttributes
  use Output_M, only: AdvancedOptions, OutputOptions, StampOptions, &
    & TimeStampOptions, InvalidPrUnit, StdoutPrUnit, MSGLogPrUnit, BothPrUnit, &
    & Output, PrUnitName !, SetTruthPattern
  use Printit_M, only: DefaultLogUnit, Get_Config, StdoutLogUnit

  implicit none
  public
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains initial or permanent settings. 
  ! Choose values according to what is most suitable for the environment.
  ! For example
  ! certain settings may be appropriate during development but not
  ! for production use, i.e. sips. Therefore, this is a convenient place
  ! to hold everything that needs to be changed before delivery.
  
  ! See also MLSL2PCF, L2ParInfo.parallel, lib/toggles.switches
  
  ! --------------------------------------------------------------------------
  ! The following should be adjusted before delivery to sips

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! The following should be TRUE if run with level 1 as a single PGE
  ! sharing a single PCF; i.e., for the near real-time (nrt)
  ! (in which case we need to move some of the "mobile" PCF ids)
  logical :: SharedPCF                               =  .false. 

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! We used to update these lines before each new delivery to the sips
  ! Now we prefer to use the --versId cmdline option or else implement it
  ! by putting "versId=xxxxx" in the optsfile.
  ! id to print out in response to "--version" command-line option  
  ! If no toolkit, then also the PGEVersion saved as a Global Attribute  
  !
  ! Beginning with v4.24/v5.02 we deliver the same compiled programs
  ! for both v4.x and v5.x versions.
  integer, parameter                        :: versIDLen = 32
  character(len=versIDLen), dimension(2)    :: current_version_id = (/ &
    & 'v4/v5 swdev team               ' , & 
    & 'See license terms for copyright'/)
  character(len=32)    :: UniqueID                    = ' '
     
  ! Set the following to MLSMSG_Error before delivering to sips;
  ! when set higher, it allows program keep going despite errors
  ! when set lower, the program would quit even on warnings
  integer            :: Quit_error_threshold          = MLSMSG_Error

  ! Set the following to MLSMSG_Error before delivering to sips;
  ! when set lower, it prints chunknum and phasename for, e.g. warnings, too.
  integer            :: Print_chunk_threshold         = MLSMSG_Warning

  ! Set the following to 2 before delivering to sips;
  ! If 0, you won't be able to distinguish normal termination
  ! from some abnormal ones (e.g. in parser) (a bad thing)
  ! if 2, status will be 2 only if run complete         
  ! and without error (a good thing)
  integer, parameter :: Normal_exit_status            = 2

  ! ---------------------------------------------------------------
  ! None of the following need to be changed before delivery to sips
  
  ! Assume HDF files w/o explicit HDFVersion field are this 
  ! 4 corresponds to HDF4, 5 to HDF5 in L2GP, L2AUX, etc.        
  integer            :: Default_HDFversion_write      = HDFversion_5
  ! Assume std. product l2gp and dgg files are the following format
  integer            :: Default_L2GPFormat_write      = L_HDFEOS ! l_netCDF4
  ! Set to WildcardHDFversion if you wish to autodetect such files  
  ! on input
  integer            :: Default_HDFversion_read       = WildcardHDFversion  
  integer            :: Level1_HDFversion             = WildcardHDFversion  

  ! The following is false only for runs that don't need orbit/attitude info
  logical            :: Need_L1BFiles                 = .true.
  ! The following is false only for runs that use non-Aura satellite data
  logical            :: Aura_L1BFiles                 = .true.
  ! Set if run must not create file, instead just append to it
  logical            :: Patch                         = .false. 
  ! Whether to restart printing identical warnings at each new phase
  logical            :: RestartWarnings               = .true.
  ! Whether to run slave tasks in the background
  ! (needed if they are to allow the wrapper script to forcibly terminate
  ! hanging slave tasks)
  logical            :: RunInBackground               = .false.
  ! File name to which to write "finished" after slave sends sig_finished
  character(len=255) :: NoteFile                      = ' '
  ! What units to use in summarizing timings at end of run
  integer            :: SectionTimingUnits            = L_seconds
  
  ! Whether to skip doing the direct writes--quicker when snooping
  logical            :: SkipDirectWrites              = .false.    
  logical            :: SkipDirectwritesOriginal      = .false.    
  ! Whether each slave deallocates all its arrays, pointers, etc.
  ! Sometimes slaves die or take too long to finish cleaning up
  ! But if system fails to reclaim memory properly, subsequent slaves
  ! may not find enough available and therefore crash
  ! FALSE means let operating system do it automatically
  logical            :: slavesCleanUpSelves           = .true.
  ! Will we be using a pvm channel to pipe data to an
  ! external running process, e.g. idl?
  logical            :: SnoopingActive                = .false.
  character(len=132) :: SnoopName                     = ''

  ! In case special dumps are to go to a special dumpfile
  character(len=255) :: SpecialDumpfile               = ' '
  ! What to fill state, outputsd with if skipping retrieval
  real               :: StateFilledBySkippedretrievals = 0.
  ! Whether to stop after a certain section: which section it is
  character(len=16)  :: StopAfterSection              = ' ' ! Blank means 'no'

  ! Whether to exit with status 1 no matter what
  logical            :: StopWithError                 = .false.         
  ! Whether to do only a pre-flight checkout of paths
  logical            :: CheckPaths                    = .false.         

  logical            :: Toolkit                       = .true. ! SIPS_VERSION 
  logical, parameter :: WriteFileAttributes           = .false.               
  ! --------------------------------------------------------------------------

  ! The following will be used only by MLSL2
  logical :: CHECKL2CF = .false.   ! Just check the l2cf and quit
  logical :: CHECKLEAK = .false.   ! Check parse tree for potential memory leaks
  logical :: COUNTCHUNKS = .false. ! Just count the chunks and quit
  integer :: DO_DUMP = 0           ! Dump declaration table if > 0
  integer :: DUMP_TREE = -1        ! Dump tree after parsing
  logical :: ExitToNextChunk = .false.   ! Skip rest of current chunk
  ! Wouldn't it be better to use get_lun at the moment we open the l2cf?
  integer, parameter :: L2CF_Unit = 20  ! Unit # if L2CF is opened by Fortran
  integer :: L2CFNode        = 0        ! Line #, Col # of L2CF being executed
  integer :: L2CFErrornode   = 0        ! Line #, Col # of L2CF at 1st error
  integer :: Numswitches
  integer :: Recl            = 20000    ! Record length for l2cf (but see --recl opt)
  integer :: MaxChunkSize    = 21     ! Max chunk size for l2gp DirectWrites
  character(len=128) :: phasesToSkip   = ''
  character(len=128) :: SectionsToSkip = ''
  logical :: SectionTimes    = .false.  ! Show times in each section
  logical :: TotalTimes      = .false.  ! Show total times from start
  logical :: ShowDefaults    = .false.  ! Just print default opts and quit
  integer :: SlaveMAF        = 0        ! Slave MAF for fwmParallel mode
  logical :: Timing          = .false.  ! -T option is set
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! This is the type to store the level 2 options that can be overridden
  ! by the l2cf commands Phase or ChangeSettings
  type :: L2Options_T
    character(len=2048) :: Command_Line ! All the opts
    character(len=2048) :: Originalcmds ! As set at launch
    character(len=32)  :: CurrentPhaseName              = ' '
    integer            :: CurrentChunkNumber            = 0
    ! Whether to skip doing the retrieval--a pre-flight checkout of paths, etc.
    logical            :: SkipRetrieval                 = .false.        
    ! A holding place for the above, allowing us to skip for some phases only
    logical            :: SkipRetrievalOriginal         = .false. 
    logical            :: MLSL2Debug                    = .false. 
    logical            :: Overridden                    = .false. ! Did we 
    ! Set the following to MSGLOGPRUNIT before delivering to sips;  override?
    ! (its possible values and their effects on normal output:
    ! INVALIDPRUNIT  on MUTE, no output
    ! STDOUTPRUNIT   sent to stdout (via print *, '...')
    ! MSGLOGPRUNIT   sent to Log file (via MLSMessage)
    ! BOTHPRUNIT     both stdout and Log file
    ! > 0            Fortran 'unit=OUTPUT_PRINT_UNIT')
    integer            :: Output_print_unit             = MSGLOGPRUNIT ! -2
    real               :: GPH_MissingValue              = -1.e15 ! -999.99

  end type L2Options_T
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! This is the type to store runtime Booleans set and used by the l2cf
  ! We sometimes call these runtime values or runtime macros

  ! Suggestion: settle on just one name and use it exclusively
  integer, parameter :: RTVStringLength               = 1024
  integer, parameter :: RTVArrayLength                = 128
  
  type :: RunTimeValues_T
    character(len=1)                   :: sep = achar(0)
    ! Two arrays bound as a {keys=>values} hash
    character(len=RTVSTRINGLENGTH)     :: lkeys       = &
      & 'true' // achar(0) // 'false' // achar(0) // 'count' 
    character(len=RTVSTRINGLENGTH)     :: lvalues     = &
      & 'true' // achar(0) // 'false' // achar(0) // 'count' 
  end type runTimeValues_T
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
  type(L2Options_T), save :: L2Options, OriginalOptions
  type(runTimeValues_T), save :: runTimeValues
  type (MLSFile_T), save      :: AllocFile

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Not used yet; the polymorphic type introduced in Fortran 2003
  ! may be a better solution.
  ! The basic idea is to be able to pass an object of any of the datatypes
  ! level 2 knows about. We might do this during a call to MLSMessage, say,
  ! if an anomalous instance is detected
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
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
  ! The following public procedures are listed here for convenience only
  ! since they are public by default anyway
  public :: DumpOptions
  public :: DumpMacros
  public :: MLSL2Message
  public :: ProcessOptions
  public :: RemoveRuntimeBoolean
  public :: RestoreDefaults

  logical, private, parameter :: countEmpty = .true. ! Except where Overridden locally
  logical, private, parameter :: AllowEnvInOpts = .true. ! Like "USEOPTSENV"

  ! The following is part of a mechanism to prevent L2Options in a .opts file
  ! from overriding either
  ! --SkipRetrieval on the cmdline, or 
  ! [n]SkipRetrieval=true
  logical, private, save      :: L2OptsCantResetSkipRetrvl = .false.
  
  type :: SomeL2Options_T
    ! Whether to skip doing the retrieval--a pre-flight checkout of paths, etc.
    logical            :: SkipRetrieval                 = .false.        
    ! A holding place for the above, allowing us to skip for some phases only
    logical            :: SkipRetrievalOriginal         = .false. 
    logical            :: MLSL2Debug                    = .false.               
    ! Set the following to MSGLOGPRUNIT before delivering to sips;
    ! (its possible values and their effects on normal output:
    ! INVALIDPRUNIT  on MUTE, no output
    ! STDOUTPRUNIT   sent to stdout (via print *, '...')
    ! MSGLOGPRUNIT   sent to Log file (via MLSMessage)
    ! BOTHPRUNIT     both stdout and Log file
    ! > 0            Fortran 'unit=OUTPUT_PRINT_UNIT')
    integer            :: Output_print_unit             = MSGLOGPRUNIT ! -2
    ! This next will be used as Missing values for the L2GPValue field
    ! in any l2gp types whose names or qty template inlude 'GPH'
    real               :: GPH_MissingValue              = -1.e15 ! -999.99

  end type SomeL2Options_T

  type(SomeL2Options_T), save :: SomeL2Options

!=============================================================================
contains 
  ! -------------------------------------------------  MLSL2Message  -----
  ! Process the level 2-savvy MLSMessage call
  ! Optionally give l2cf line #, dump an item, and so on before summoning
  ! regular MLSMessage
  !
  ! For some cases we skip calling MLSMessage, either calling output instead
  ! or remaining mute
  !
  ! In certain other cases we must repeat printing the message via the
  ! output module's commands
  
  ! The above is too vague. What determines these cases? Specifically
  ! (a) When  do we call output instead of MLSL2Message? When are we mute?
  ! (b) When must we call output and MLSL2Message both?
  subroutine MLSL2Message ( Severity, ModuleNameIn, Message, &
    & Advance, MLSFile, Status, Item )
    use Lexer_Core, only: Get_Where
    use MLSStringLists, only: SwitchDetail
    use Toggles, only: Switches
    use Tree, only: Where
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
    character(len=16) :: chunkChars
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
    end if
    if ( L2Options%MLSL2Debug ) then
      print *, 'mustRepeat ', mustRepeat
      print *, 'outputInstead ', outputInstead
      print *, 'STDOUTPRUNIT ', STDOUTPRUNIT
      print *, 'OutputOptions%PrUnit ', OutputOptions%PrUnit
      print *, 'DEFAULTLOGUNIT ', DEFAULTLOGUNIT
      print *, 'Config%logFileUnit ', logFileUnit
    end if
    ! Do we have a current phase name?
    if ( len_trim(L2Options%currentPhaseName) > 0 &
      & .and. severity >= Print_chunk_threshold ) then
      myMessage = ' (' // trim(L2Options%currentPhaseName) // ') ' // myMessage
    end if
    ! Do we have a current chunk number?
    if ( severity >= Print_chunk_threshold ) then
      if ( L2Options%CurrentChunkNumber > 0 ) then
        call writeIntsToChars( L2Options%CurrentChunkNumber, chunkChars )
        myMessage = ' (chunk ' // trim(chunkChars) // ') ' // myMessage
      else
        myMessage = ' (no chunks) ' //  myMessage
      end if
    end if
    ! For severity "info" just do the call and return
    ! For severity "warn" do the call and 
    ! check status to see if we need to skip printing
    if ( severity == MLSMSG_INFO ) then
      call SayIt ( Message )
      return
    else if ( severity == MLSMSG_Warning ) then
      ! This trickery is to determine whether this warning would be suppressed
      call SayMessage ( MLSMSG_TestWarning, ModuleNameIn, MyMessage, &
        & Advance, MLSFile, myStatus )
      if ( present(status) ) status = myStatus
      if ( myStatus /= 0 ) return
    end if
    ! If severe enough an error, we will want to know more
    ! in hopes of figuring out what went wrong, where, and why
    ! Do we have an l2cf node we were processing?
    if ( L2CFErrorNode /= 0 ) then
      call get_where ( where(L2CFErrorNode), myMessage, &
        & before='***** At ', after=': ' // myMessage )
    elseif ( L2CFNode /= 0 ) then
      call get_where ( where(L2CFNode), myMessage, &
        & before='***** At ', after=': ' // myMessage )
    else
      myMessage = '(no tree node): ' // myMessage
    end if
    ! Were we passed an item of a recognized MLS data type to Dump?
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
    end if
    call SayIt ( myMessage )
  contains
    subroutine SayIt ( It )
      ! Say it with MLSMessage
      ! and possibly repeat it with output
      character (len=*), intent(in) :: It
      character(len=256)            :: mesg
      integer, parameter            :: LineLength = 40
      if ( Severity > MLSMSG_Warning ) then
        mesg = trim(it) // ' (in) ' // trim(ModuleNameIn)
        if ( present(MLSFile) ) call dump( MLSFile )
        if ( mustRepeat ) call banner( trim(mesg), LineLength=LineLength )
        if ( .not. outputInstead ) call Bummer ( mesg, &
          & LineLength=LineLength, severity=severity )
      else
        if ( .not. outputInstead ) call SayMessage ( severity, ModuleNameIn, It, &
          & Advance, MLSFile, status )
        if ( mustRepeat ) call output( trim(It), advance='yes' )
      endif
    end subroutine SayIt
  end subroutine MLSL2Message

  ! ---------------------------------------------  ProcessOptions  -----
  ! Process the command line options; either from
  ! (1) the command line directly, by calls to getNextArg; or
  ! (2) parsing cmdline, if supplied as an arg
  ! The return value will be the filename (if any)
  function ProcessOptions ( cmdLine, separator ) result ( fileName )
    use Allocate_Deallocate, only: AllocateLogUnit, TrackAllocates, &
      & ClearOnAllocate, Allocate_Test, Deallocate_Test
    use Evaluate_Variable_M, only: Define_Variable_As_String
    use Io_Stuff, only: Get_Lun, Get_NLines, Read_TextFile
    use L2ParInfo, only: Parallel, InitParallel, AccumulateSlaveArguments, &
      & SnipLastSlaveArgument
    use Lexer_M, only: CapIdentifiers
    use Machine, only: Getarg, Hp, Io_Error, NeverCrash
    use Matrixmodule_0, only: CheckBlocks, SubblockLength
    use MLSCommon, only: FileNameLen
    use MLSMessageModule, only: Setconfig
    use MLSStringLists, only: Catlists, &
      & Getstringelement, Getuniquelist, &
      & Numstringelements, Removeswitchfromlist, &
      & Stringelement, Switchdetail, Unquote
    use PCFHdr, only: GlobalAttributes
    use Printit_M, only: Set_Config
    use Set_Toggles_M, only: Set_Toggles
    use String_Table, only: Add_Include_Directory, Do_Listing
    use Time_M, only: Time_Config
    use Toggles, only: Switches
    ! Args
    character(len=*), intent(in), optional :: cmdLine
    character(len=1), intent(in), optional :: separator
    character(len=FileNameLen)             :: fileName
    ! Internal variables
    character(len=1) :: arg_rhs      ! 'n' part of 'arg=n'
    character(len=16) :: aSwitch
    integer :: Degree                ! index affecting degree of option
    logical, parameter :: DEEBUG = .false.
    logical :: Exist
    logical :: exitLoop
    logical :: Field_Is_Include      ! Field is include file path, not L2CF name
    integer :: I, J
    character(len=2048) :: LINE      ! Into which is read the command args
    integer :: N
    integer :: nLines
    logical :: Opened
    character(len=FileNameLen)             :: optsFile
    character(len=FileNameLen), dimension(:), pointer  :: optLines => null()
    character(len=2) :: quotes
    integer :: Recl = 256          ! Record length for allocation log
    character(len=len(switches)) :: removeSwitches = ''
    integer :: Status
    logical :: Switch                ! "First letter after -- was not n"
    character(len=len(switches)) :: tempSwitches
    integer :: V
    character(len=2048) :: Word      ! Some text
    ! Executable
    quotes = char(34) // char(39)   ! {'"}
    filename = 'helphelphelp' ! This means abnormal options--should dump help mesg
    if ( present( cmdline ) .and. DEEBUG ) then
      print *, 'cmdline: ', trim(cmdline)
    end if
    ! Before looking at command-line options, TOOLKIT is set to true
    ! So here's a good place to put any SIPS-specific settings
    ! --- SIPS_VERSION ---
    parallel%maxFailuresPerMachine = 2
    parallel%maxFailuresPerChunk = 1
    removeSwitches='slv' ! Since slave output already saved to separate files
    ! switches='red'  ! No longer a good idea
    DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
    MLSMessageConfig%limitWarnings = 4 ! Print less
    time_config%use_wall_clock = .true. ! SIPS_VERSION
    ! call SetTruthPattern ( (/ '+ ', '- ' /) ) ! Print these instead of 'T' 'F'
    i = 1 + hp
    ! ----------------- deprecated and deplored ----------------------------
    ! Man, we couldn't use the Lahey compiler even if we wanted to!
    ! Next time, let's remove these commented-out lines
    ! do ! Process Lahey/Fujitsu run-time options; they begin with "-Wl,"
    !   call getNextArg ( i, line )
    !   if ( line(1:4) /= '-Wl,' ) then
    !     call SnipLastSlaveArgument ! Don't want slaves to see this
    !     exit
    !   end if
    !   call AccumulateSlaveArguments(line) ! pass them to slave processes
    !   i = i + 1
    ! end do
    ! Now process the other options
    ! -------------------------------------------------------------
    L2Options%command_line = ' '
    field_is_include = .false.
cmds: do
      call getNextArg ( i, line )
      if ( DEEBUG ) print *, i, trim(line)
      if ( len_trim( line ) < 1 ) then
        exit
      end if
      call processLine( line, filename, exitLoop, entireLine=.false. )
      if ( DEEBUG ) print *, 'filename: ', trim(filename)
      if ( exitLoop ) exit
      i = i + 1
    enddo cmds
    ! Did we somehow miss the filename among the args?
    if ( filename == 'helphelphelp' ) then
      if ( DEEBUG ) print *, 'did we miss the filename?'
      i = i + 1
      call getNextArg ( i, line )
      filename = line
      if ( DEEBUG ) print *, 'filename: ', trim(filename)
    end if
  ! Are any switches inappropriate for master or for slave?
    if ( parallel%master ) &
      & removeSwitches = catLists( trim(removeSwitches), 'bool,walk' )
    if ( parallel%slave ) &
      & removeSwitches = catLists( trim(removeSwitches), 'chu,l2q,mas,slv' )
    ! Remove any quote marks from RemoveSwitches array
    tempSwitches = unquote( removeSwitches, quotes=quotes, options='-p' )
    call GetUniqueList( tempSwitches, removeSwitches, numSwitches, &
          & ignoreLeadingSpaces=.true., options='-eSL' )
    ! Remove any quote marks from switches array
    tempSwitches = switches
    Switches = unquote( tempSwitches, quotes=quotes, options='-p' )
    ! Remove any switches embedded in the removeSwitches option 'R'
    if ( DEEBUG ) print *, 'starting List ', trim(switches) 
    if ( DEEBUG ) print *, 'to remove ', trim(removeSwitches) 
    do i=1, NumStringElements( removeSwitches, countEmpty=.true. )
      call GetStringElement( trim(removeSwitches), aSwitch, i, &
        & countEmpty=.true. )
      call RemoveSwitchFromList( switches, tempSwitches, trim(aSwitch) )
      switches = tempSwitches
    end do
    if ( DEEBUG ) print *, 'result List ', trim(switches) 

    ! If we like, we could move these next few statements to a standalone
    ! subroutine named something like processSwitches
    parallel%verbosity = switchDetail(switches, 'mas') + 1
    if ( switchDetail(switches, 'walk') > -1 ) &
      & MLSMSG_Severity_to_walkback = MLSMSG_Warning
    if (  L2Options%Overridden ) then
      outputOptions%prunit =  L2Options%Output_print_unit
      MLSMessageConfig%logFileUnit =  L2Options%Output_print_unit
      if ( DEEBUG ) &
        & print *, 'Setting logFileUnit ', MLSMessageConfig%logFileUnit
      if ( any( L2Options%Output_print_unit == &
        & (/ MSGLOGPRUNIT, BOTHPRUNIT /) ) ) MLSMessageConfig%useToolkit = .true.
    endif
    ! print *, 'Ended processing options'
  contains
    subroutine getNextArg( i, line )
      ! Args
      integer, intent(in)           :: i
      character(len=*), intent(out) :: line
      logical, parameter :: countEmpty = .false.
      ! Executable
      if ( present( cmdline ) ) then
        line = StringElement ( cmdline, i, countEmpty, inseparator=separator )
        ! print *, 'Using cmdline:', trim(cmdline)
        ! print *, trim(line)
      else
        call getArg ( i, line )
        ! print *, 'arg #:', i
        ! print *, trim(line)
      end if
      ! Check for embedded spaces
      if ( index(trim(line), ' ') > 0 ) then
        line = "'" // trim(line) // "'"
      endif
      L2Options%command_line = trim(L2Options%command_line) // ' ' // trim(line)
      call AccumulateSlaveArguments( line ) ! pass them to slave processes
    end subroutine getNextArg

    ! process the portion of the whole cmdline beginning with the 
    ! current arg's start
    ! recursive because some options, --set and --setf,
    ! provoke more calls to this subroutine
    recursive subroutine processLine( inLine, filename, exitLoop, entireLine  )
      character(len=*), intent(in) :: inLine
      character(len=*)             :: filename
      logical, intent(out)         :: exitLoop
      logical, intent(in)          :: entireLine ! is inLine entire cmdline?
      ! Internal variables
      integer :: iActual, K
      character(len=2048) :: LINE      ! Into which is read the command args
      ! Executable
      line = inLine
      exitLoop = .false.
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
        else if ( line(3+n:7+n) == 'backg' ) then
          runinbackground = switch
        else if ( line(3+n:7) == 'chunk' ) then
          i = i + 1
          ! print *, 'About to read chunk num'
          call myNextArgument( i, inLine, entireLine, line )
          parallel%chunkRange = line
          ! print *, 'parallel%chunkRange: ', trim(parallel%chunkRange)
        else if ( line(3+n:18+n) == 'clearonallocate ' ) then
          clearonallocate = switch
        else if ( line(3+n:7+n) == 'ckbk ' ) then
          checkBlocks = switch
        else if ( line(3+n:14+n) == 'countchunks ' ) then
          countChunks = switch
        else if ( line(3+n:7+n) == 'crash' ) then
          MLSMessageConfig%crashOnAnyError = switch
          neverCrash = .not. switch
        else if ( line(3+n:8+n) == 'curver' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
        else if ( line(3+n:8+n) == 'defaul' ) then
          showDefaults = switch
        else if ( line(3+n:7+n) == 'delay' ) then
          if ( line(8+n:) /= ' ' ) then
            line(:7+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
          end if
          read ( line, *, iostat=status ) parallel%Delay
          if ( status /= 0 ) then
            call io_error ( "After --delay option", status, line )
            stop
          end if
        else if ( line(3+n:7+n) == 'dump ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          specialDumpFile = trim(line)
        else if ( line(3+n:14+n) == 'fwmparallel ' ) then
          parallel%fwmParallel = .true.
        else if ( line(3+n:7+n) == 'host ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          GlobalAttributes%hostName = trim(line)
        else if ( line(3+n:9+n) == 'idents ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
        else if ( line(3+n:7+n) == 'inst ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          GlobalAttributes%InstrumentName = trim(line)
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
            return
          end if
          call myNextArgument( i, inLine, entireLine, line )
          read ( line, *, iostat=status ) degree
          if ( status /= 0 ) then
            call io_error ( "After --lac[onic] option", status, line )
            stop
          end if
          ! MLSMessageConfig%suppressDebugs = (degree > 0) ! Why do this?
          ! MLSMessageConfig%AbbreviateModSevNames = (degree == 0)
          if ( degree == 0 ) then
            MLSMessageConfig%MaxModuleNameLength     = 3
            MLSMessageConfig%MaxSeverityNameLength   = 1
          else if ( degree == 3 ) then
            MLSMessageConfig%MaxModuleNameLength     = 0
            MLSMessageConfig%MaxSeverityNameLength   = 0
          end if
          MLSMessageConfig%skipModuleNamesThr = mod(degree, 10) + 2
          MLSMessageConfig%skipSeverityThr = mod(degree, 10) + 2
          MLSMessageConfig%skipMessageThr = degree - 10 + 2
          if ( degree > 9 ) then
            removeSwitches = catLists(trim(removeSwitches), 'log' )
            outputOptions%prunit = INVALIDPRUNIT
          end if
        else if ( line(3+n:7+n) == 'leak' ) then
          checkLeak = .true.
        else if ( line(3+n:5+n) == 'loc' ) then
          call SnipLastSlaveArgument ! Don't want slaves to see this
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          call SnipLastSlaveArgument ! Don't want slaves to see this
          GlobalAttributes%productionLoc = trim(line)
          if ( len_trim(GlobalAttributes%hostName) < 1 ) &
            & GlobalAttributes%hostName = trim(line)
        else if ( line(3+n:10+n) ==  'logalloc' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, filename )
          call get_lun ( AllocateLogUnit, msg=.true. )
          ! We can't call InitializeMLSFile yet because
          ! the (string, character) tables are still empty and we would seg fault
          ! (can't we fix InitializeMLSFile somehow?)
          
          ! So instead let's do all that "by hand"
          AllocFile%name          = filename
          ! AllocFile%type = l_ascii
          AllocFile%content       = 'logAlloc'
          AllocFile%access        = DFACC_RDWR
          AllocFile%FileId%f_id  = AllocateLogUnit
          open ( unit=AllocateLogUnit, access='sequential', action='readwrite', &
            & form='formatted', &
            & status='unknown', file=trim(fileName), recl=recl, iostat=status )
          if ( status /= 0 ) then
            call io_error ( "Failed to open file to log allocates", &
              & status, filename )
            stop
          else
            print *, 'opened allocations log file ' // trim(fileName) &
              &  // ' unit ', AllocateLogUnit
          endif
        else if ( line(3+n:9+n) == 'master ' ) then
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%master = .true.
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%slaveFilename = trim ( line )
          call InitParallel ( 0, 0 )
          word = '--slave'
          write ( word(len_trim(word)+1:), * ) parallel%myTid
          call AccumulateSlaveArguments(word)
        else if ( line(3+n:14+n) == 'maxchunksize' ) then
          if ( line(15+n:) /= ' ' ) then
            line(:14+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
          end if
          read ( line, *, iostat=status ) maxChunkSize
          if ( status /= 0 ) then
            call io_error ( "After --maxChunkSize option", status, line )
            stop
          end if
        else if ( line(3+n:21+n) == 'maxfailuresperchunk' ) then
          if ( line(22+n:) /= ' ' ) then
            line(:21+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
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
            call myNextArgument( i, inLine, entireLine, line )
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
            call myNextArgument( i, inLine, entireLine, line )
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
            call myNextArgument( i, inLine, entireLine, line(11:) )
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
          call myNextArgument( i, inLine, entireLine, line )
          call SnipLastSlaveArgument ! Don't want slaves to see this
          parallel%pgeName = trim(line)
        else if ( line(3+n:6) == 'pidf' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, noteFile )
        else if ( line(3+n:6+n) == 'quit' ) then
          if ( .not. switch ) &
            & QUIT_ERROR_THRESHOLD = MLSMSG_Error + 10 ! Try not to quit
        else if ( line(3+n:6+n) == 'recl' ) then
          if ( line(7+n:) /= ' ' ) then
            line(:6+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
          end if
          read ( line, *, iostat=status ) recl
          if ( status /= 0 ) then
            call io_error ( "After --recl option", status, line )
            stop
          end if
        else if ( line(3+n:6+n) == 'set ' ) then
          ! While it may seem odd, there may be cases when it is convenient
          ! e.g., when reading a set of name=value pairs from a file
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          ! Must save i, which will otherwise be incremented
          iActual = i
          call parseNameValue ( line )
        else if ( line(3+n:6+n) == 'setf' ) then
          i = i + 1
          call getNextArg ( i, optsFile )
          call get_nLines ( optsFile, nLines )
          if ( nLines < 0 ) then
            print *, 'Sorry, unable to open opts file ' // trim(optsFile)
            call MLSL2Message( MLSMSG_Error, ModuleName, &
              & 'Sorry, unable to open opts file ' // trim(optsFile) )
          elseif ( nLines < 1 ) then
            print *, '0 lines in opts file ' // trim(optsFile)
            call MLSL2Message( MLSMSG_Warning, ModuleName, &
              & 'Unexpectedly empty opts file ' // trim(optsFile) )
          endif
          call allocate_test ( optLines, nLines, 'optLines', &
            & trim(ModuleName) // 'processLine' )
          optLines = ' '
          call read_textFile( optsFile, optLines )
          print *, 'options read from ' // trim(optsFile)
          print *, 'nLines', nLines 
          do k=1, nLines
            print *,  trim(optLines(k))
          enddo
          ! Must save i, which will otherwise be incremented
          iActual = i
          do k=1, nLines
            ! Ignore comments and blank lines
            if ( DEEBUG ) print *, 'len(line)', len_trim(optLines(k)) 
            if ( DEEBUG ) print *, 'is Comment', isComment( optLines(k), '#' ) 
            if ( len_trim(optLines(k)) < 1 .or. &
              & isComment( optLines(k), '#' ) ) cycle
            !
            ! Do we ignore lines w/o the "=" sign?
            ! They could be used to turn on special flags to pre-process
            ! the .opts file
            ! e.g.
            !   USEOPTSENV
            !   nightday=$NIGHTDAY
            ! and if you set the global environmental parameter
            !    NIGHTDAY=day
            ! then mlsl2 would see
            !   nightday=day
            !
            ! If we don't ignore such lines, they will cause a parse error
            !
            ! The parameter ALLOWENVINOPTS, if TRUE.,
            ! allows such lines
            if ( index( optLines(k), '=' ) < 1 ) then
              if ( ALLOWENVINOPTS ) then
                cycle
              else
                print *, 'Sorry, unable to parse this opts line'
                print *, '(expected something like lhs=rhs)'
                print *, 'To allow environment flags in the opts file'
                print *, 'MLSL2Options.f90 must be recompiled '
                print *, 'after setting ALLOWENVINOPTS = .true. '
                print *, trim(optLines(k))
                call MLSL2Message( MLSMSG_Error, ModuleName, &
                  & 'Sorry, Parse error in opts file ' // trim(optsFile) )
              endif
            endif
            if ( DEEBUG ) print *, 'optLines(k)', trim(optLines(k)) 
            call parseNameValue( optLines(k) )
          enddo
          call deallocate_test ( optLines, 'optLines', &
            & trim(ModuleName) // 'processLine' )
          i = iActual ! Restore i
        else if ( line(3+n:6+n) == 'shar' ) then
          SHAREDPCF = switch
        else if ( line(3+n:9+n)  == 'skipdir' ) then
          SKIPDIRECTWRITES = switch
        else if ( line(3+n:10+n) == 'skipretr' ) then
          L2Options%SkipRetrieval = switch .or. L2Options%SkipRetrievalOriginal
          L2OptsCantResetSkipRetrvl = .true.
        else if ( line(3+n:8+n) == 'skipph' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          phasesToSkip = lowercase(line)
        else if ( line(3+n:9+n) == 'skipsec' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          sectionsToSkip = lowercase(line)
        else if ( line(3+n:10+n) == 'slavemaf' ) then
          if ( line(11+n:) /= ' ' ) then
            line(:10+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
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
            call myNextArgument( i, inLine, entireLine, line )
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
          call myNextArgument( i, inLine, entireLine, line )
          snoopName = line
        else if ( line(3+n:7+n) == 'state' ) then
          if ( line(8+n:) /= ' ' ) then
            line(:7+n) = ' '
          else
            i = i + 1
            call myNextArgument( i, inLine, entireLine, line )
          end if
          read ( line, *, iostat=status ) stateFilledBySkippedRetrievals
          if ( status /= 0 ) then
            call io_error ( "After --state option", status, line )
            stop
          end if
        else if ( line(3+n:9+n) == 'stdout ' ) then
          ! print *, 'Got --stdout option'
          if ( .not. switch ) then
            L2Options%Output_print_unit = INVALIDPRUNIT
            return
          end if
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          ! print *, 'Its arg: ', trim(line)
          select case ( lowercase(line) )
          case ( 'log' )
            L2Options%Output_print_unit = MSGLOGPRUNIT
          case ( 'out' )
            L2Options%Output_print_unit = STDOUTPRUNIT
          case ( 'both' )
            L2Options%Output_print_unit = BOTHPRUNIT
            ! MLSMessageConfig%logFileUnit = STDOUTLOGUNIT
          case ( 'error' ) ! Send Output to stderr instead of stdout
            outputOptions%prUnit = error_unit
            outputOptions%prUnitLiteral = .true. ! In case Error_Unit <= 0
            OutputOptions%buffered = .false.
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
            L2Options%Output_print_unit = OutputOptions%prUnit
          end select
          if (  L2Options%Overridden ) then
            outputOptions%prunit =  L2Options%Output_print_unit
            MLSMessageConfig%logFileUnit =  L2Options%Output_print_unit
            if ( any( L2Options%Output_print_unit == &
              & (/ MSGLOGPRUNIT, BOTHPRUNIT /) ) ) &
              & MLSMessageConfig%useToolkit = .true.
            if ( DeeBug ) &
              & print *, 'Setting logFileUnit ', MLSMessageConfig%logFileUnit
          endif
        else if ( line(3+n:12+n) ==  'stopafter ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, stopAfterSection )
        else if ( line(3+n:12+n) ==  'stopwither' ) then
          stopWithError = switch
        else if ( line(3+n:11+n) == 'subblock ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          read ( line, *, iostat=status ) subBlockLength
          if ( status /= 0 ) then
            call io_error ( "After --subblock option", status, line )
            stop
          end if
        else if ( line(3+n:9+n) == 'submit ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          parallel%submit = trim ( line )
        else if ( line(3+n:5+n) == 'tk ' ) then
          toolkit = switch
          call set_config ( useToolkit = switch )
        else if ( line(3+n:5+n) == 'uid' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, uniqueId )
        else if ( line(3+n:6+n) == 'var ' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line ) ! The variable
          i = i + 1
          call myNextArgument( i, inLine, entireLine, word ) ! Its value as a string
          call define_variable_as_string ( line, word )
        else if ( line(3+n:10+n) == 'verbose ' ) then
          switches = catLists( trim(switches), &
            & 'l2q,glob,mas,bool1,opt1,log,pro1,time,apr,phase' )
          StampOptions%showTime = .true.
          StampOptions%dateFormat = 'yyyy-mm-dd'
        else if ( line(3+n:8+n) == 'versid' ) then
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          ! Undo the substitution of '%' for space done during parseNameValue
          current_version_id(1) = Replace( trim(line), '%', ' ' )
          GlobalAttributes%PGEVersion = current_version_id(1)
          ! print *, 'Changed current version id'
          ! do j=1, size(current_version_id)
          !   print *, current_version_id(j)
          ! end do
        else if ( line(3+n:10+n) == 'version ' ) then
          print '(a)', current_version_id
          stop
        else if ( line(3+n:7+n) == 'wall ' ) then
          time_config%use_wall_clock = switch
        else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
          i = i + 1
          call myNextArgument( i, inLine, entireLine, line )
          exitLoop = .true.
          return
        else
          if ( index( line, '-help' ) < 1 ) &
            & print *, 'unrecognized option ', trim(line), ' ignored.'
          ! call option_usage
          filename = 'help'
          exitLoop = .true.
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
          case ( 'E' ) ! Output usually sent to stdout goes to stderr
            outputOptions%prUnit = error_unit
            outputOptions%prUnitLiteral = .true. ! In case Error_Unit <= 0
            OutputOptions%buffered = .false.
          case ( 'h', 'H', '?' )     ! Describe command line usage
            ! call option_usage
            filename = 'help'
            exitLoop = .true.
            return ! will dump help mesg
          case ( 'I' ); field_is_include = .true. ! Next field is include path
          case ( 'K' ); capIdentifiers = .true.
          case ( 'k' ); capIdentifiers = .false.
          case ( 'M' ); outputOptions%prunit = MSGLOGPRUNIT ! -2
          case ( 'm' ); outputOptions%prunit = STDOUTPRUNIT ! -1
          case ( 'R' ) ! This does the opposite of what S does
            removeSwitches = catLists(trim(removeSwitches), line(j+1:))
            if ( DeeBug ) call Dump ( removeSwitches, 'switches to be removed' )
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
        exitLoop = .true. ! This must be the l2cf filename
        print *, 'This must be the l2cf filename ' // trim(line)
        return
      end if
    end subroutine ProcessLine

    ! Separate lhs from rhs, then treat 
    ! as cmdline options depending on whether
    ! lhs is one character in length or two
    ! rhs is true, false, or a different value
    recursive subroutine parseNameValue( inLine )
      character(len=*), intent(in) :: inLine
      ! Local variables
      logical, parameter :: DEEBUG = .false.
      logical :: exitLoop
      character(len=fileNameLen)      :: filename
      character(len=len(inLine)+32)   :: line
      character(len=fileNameLen)      :: name
      character(len=fileNameLen)      :: valu
      logical, parameter              :: countEmpty = .true.
      character(len=1), parameter     :: eqls = '='
      ! Executable
      call SubstituteRuntimeBoolean ( inLine, line )
      call getStringElement( line, name, 1, countEmpty, eqls )
      call getStringElement( line, valu, 2, countEmpty, eqls )

      ! -------------------------------------------------------------
      ! This mechanism optionally sets level 2 to crash with
      ! a walkback
      ! if it logs a message containing a fatal string
      ! E.g., put the next line in your opts file
      ! CrashIfMsgSays=Drop. Dead.
      ! and your run will automatically crash at the point where it logs
      ! any message containing the string "Drop. Dead."
      if ( lowercase(name) == 'crashifmsgsays' ) then
        NeverCrash = .false.
        MLSMessageConfig%CrashIfMsgSays = valu
        return
      endif
      ! See also open_init for the same mechanism implemented in the PCF
      ! -------------------------------------------------------------
      
      ! ---------------------------------------------------------------
      ! Special means for setting user-defined types:
      ! OutputOptions, AdvancedOptions, 
      ! StampOptions, or TimeStampOptions
      ! Of course, here is where you could
      ! add similar means for setting other user-defined types
      ! ---------------------------------------------------------------
      if ( lowercase(name) == 'outputoptions' ) then
        read(valu,*) OutputOptions
        return
      elseif ( lowercase(name) == 'advancedoptions' ) then
        read(valu,*) AdvancedOptions
        return
      elseif ( lowercase(name) == 'stampoptions' ) then
        read(valu,*) StampOptions
        return
      elseif ( lowercase(name) == 'timestampoptions' ) then
        read(valu,*) TimeStampOptions
        return
      elseif ( lowercase(name) == 'globalattributes' ) then
        read(valu,*) SomeGlobalAttributes
        call SomeToGlobalAttributes
        return
      elseif ( lowercase(name) == 'l2options' ) then
        if ( DEEBUG ) print *, 'Reading SomeL2Options'
        read(valu,*) SomeL2Options
        call SomeToL2Options
        call DumpOptions
        return
      endif
      ! ---------------------------------------------------------------

      ! Beware of cases where valu contains an embedded space
      ! Replace such spaces with '%'
      if ( len_trim(name) < 1 .or. name == eqls ) return
      if ( len_trim(valu) > 0 ) valu = Replace ( trim(valu), ' ', '%' )

      ! Special cases:
      ! We prefer to store these in valu
      ! because it is long enough to hold, e.g., mutiple switches
      ! print *, trim(name) // ' set = to ' // trim(valu)
      if ( lowercase(name) == 'switches' ) then
        valu = '-S' // valu
        call ProcessLine( trim(valu), filename, exitLoop, entireLine=.true. )
        return
      elseif ( lowercase(name) == 'remove' ) then
        valu = '-R' // valu
        call ProcessLine( trim(valu), filename, exitLoop, entireLine=.true. )
        if ( DEEBUG ) print *, 'encountered a remove name/value pair ' // trim(valu)
        return
      elseif ( lowercase(name) == 'debugging' ) then
        valu = '-D' // valu
        call ProcessLine( trim(valu), filename, exitLoop, entireLine=.true. )
        return
      elseif ( lowercase(name) == 'verboseness' ) then
        valu = '-V' // valu
        call ProcessLine( trim(valu), filename, exitLoop, entireLine=.true. )
        return
      endif

      ! Standard cases: either single-character or multi-character
      ! If false, prefix name with 'n'
      if ( valu == 'false' ) name = 'n' // name
      if ( len_trim(name ) < 2 ) then
        ! single-character options
        name = '-' // name
      elseif ( isDigits(name(2:) ) ) then
        ! single-character options plus "degree"; e.g. -f4
        name = '-' // name
      else
        ! multi-character options
        name = '--' // name
      endif
      if ( any( valu == (/ 'true ', 'false' /) ) ) then
        ! print *, 'processing cmdline option ', trim(name)
        call ProcessLine( trim(name), filename, exitLoop, entireLine=.true. )
      else
        ! print *, 'processing cmdline option ', trim(name) // ' ' // trim(valu)
        call ProcessLine( trim(name) // ' ' // trim(valu), filename, exitLoop, &
          & entireLine=.true. )
      endif
    end subroutine parseNameValue
    
    ! Return either 
    ! (1) The ith command line argument (if entireLine is false)
    ! (2) The remainder of inLine       (if entireLine is true)
    subroutine myNextArgument( i, inLine, entireLine, line )
      ! Args
      integer, intent(in)           :: i
      character(len=*), intent(in)  :: inLine
      logical, intent(in)           :: entireLine
      character(len=*), intent(out) :: line
      ! Executable
      if ( .not. entireLine ) then
        call getNextArg ( i, line )
        ! print *, 'not entireLine'
        ! print *, trim(line)
      else
        line = stringElement( adjustl( inLine ), 2, &
          & countEmpty=.false., inseparator=' ' )
        ! print *, 'entireLine'
        ! print *, trim(line)
      endif
      ! print *, 'my next arg: ', trim(line)
    end subroutine myNextArgument
    
  end function ProcessOptions
  
  ! --------------------------------------------  RestoreDefaults  -----
  ! Restore the options to their default values
  ! Now some things it makes no sense to overwrite, so it makes
  ! no sense to restore them either; e.g., CHECKPATHS, parallel, etc.
  subroutine RestoreDefaults ( complete )
  use MLSMessageModule, only: RestoreConfig
  use Toggles, only: Init_Toggle
    logical, intent(in), optional            :: complete ! Restore even quit, crash!
    ! Internal variables
    logical, parameter :: DEEBUG = .false.
    logical                                  :: myComplete
    ! Executable
    myComplete                               = present(complete)
    if ( myComplete )   myComplete           = complete
    if ( DEEBUG ) print *, 'Entered RestoreDefaults; myComplete ', myComplete
    L2Options%Output_print_unit             = -2
    Default_hdfversion_write      = HDFVERSION_5
    Default_hdfversion_read       = WILDCARDHDFVERSION
    Level1_hdfversion             = WILDCARDHDFVERSION
    Restartwarnings               = .true.
    Sectiontimingunits            = L_SECONDS
    Skipdirectwrites              = .false.    
    Skipdirectwritesoriginal      = .false.    
    ! Skipretrieval                 = .false.        
    ! Skipretrievaloriginal         = .false. 
    SlavesCleanUpSelves           = .true.
    Specialdumpfile               = ' '
    Statefilledbyskippedretrievals = 0.
    Stopaftersection              = ' ' ! Blank means 
    Stopwitherror                 = .false.         
    call init_toggle

    if ( .not. myComplete ) return
    Need_L1BFiles                 = .true.
    Aura_L1BFiles                 = .true.
    Patch                         = .false. 
    Checkpaths                    = .false.         
    Toolkit                       =  .true. ! SIPS_VERSION
    call restoreConfig ( complete )
  end subroutine RestoreDefaults

  ! --------------------------------------------  SomeToL2Options  -----
  ! Restore the options to their default values
  ! Now some things it makes no sense to overwrite, so it makes
  ! no sense to restore them either; e.g., CHECKPATHS, parallel, etc.
  subroutine SomeToL2Options
    logical, parameter :: DEEBUG = .false.
    ! Executable
    if ( DEEBUG ) print *, 'Entered SomeToL2Options'
    ! Here's the mechanism to prevent L2Options from resetting either
    ! --SkipRetrieval on the cmdline, or
    ! [n]SkipRetrieval=true
    if ( .not. L2OptsCantResetSkipRetrvl ) then
    L2Options%SkipRetrieval           =  SomeL2Options%SkipRetrieval        
    else
      call output( 'L2Options not permitted to reset SkipRetrieval', &
        & advance='yes' )
    endif
    L2Options%SkipRetrievalOriginal   =  SomeL2Options%SkipRetrievalOriginal
    L2Options%MLSL2Debug              =  SomeL2Options%MLSL2Debug           
    L2Options%Output_print_unit       =  SomeL2Options%Output_print_unit    
    L2Options%GPH_MissingValue        =  SomeL2Options%GPH_MissingValue     
    L2Options%Overridden              =   .true.
  end subroutine SomeToL2Options

  ! -------------------------------------------------  DumpOptions  -----
  ! Dump the L2 Options
  ! Print the cmd line in blocs of BlocLength characters
  ! We have not yet decided how to use the details arg
  subroutine DumpOptions ( details )
    use HighOutput, only: AddRow, AddRow_Divider, AddRow_Header, &
      & OutputTable, StartTable
    integer, optional, intent(in)                :: Details ! Not used at present

    call startTable
    call addRow_header ( 'Current Level 2 Options', 'c' )
    call addRow_divider ( '-' )
    call addRow ('Phase            ', trim( L2Options%CurrentPhaseName        ) )
    call addRow ('Chunk            ', L2Options%CurrentChunkNumber         )
    call addRow ('SkipRetrieval    ', L2Options%Skipretrieval         )
    call addRow ('Send output to   ', trim( PrUnitname( L2Options%Output_print_unit )   ) )
    call addRow ('MLSL2Debug       ', L2Options%MLSL2Debug         )
    call addRow ('Overridden       ', L2Options%Overridden         )
    call addRow ('GPH MissingValue ', L2Options%GPH_MissingValue     )
    call addRow ('Cmdline          ', trim( L2Options%Command_line        ), &
      & BlocLen=40, options='-w' )
    call outputTable ( sep='|', border='-' )
  end subroutine DumpOptions

  ! -------------------------------------------------  DumpMacros  -----
  ! Dump the runtime macros
  ! Either simply, or in a table
  subroutine DumpMacros ( details )
    use Dump_1, only: Dump
    use HighOutput, only: OutputTable
    use MLSStringLists, only: List2Array, NumStringElements, SwitchDetail
    use Toggles, only: Switches
    integer, optional, intent(in)           :: Details
    character(len=256), dimension(1024, 2)  :: KeysValues
    integer                                 :: MyDetails
    integer                                 :: NValues
    ! Executable
    myDetails = SwitchDetail( switches, 'bool' )
    if ( present(details) ) myDetails = details
    
    if ( myDetails < -1 ) then
      call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
        & 'Run-time macros', separator=runTimeValues%sep )
    else
      call List2Array( runtimevalues%lkeys, keysValues(2:,1), countEmpty, &
        & inseparator=runTimeValues%sep, &
        & ignoreLeadingSpaces=.true. )
      call List2Array( runtimevalues%lvalues, keysValues(2:,2), countEmpty, &
        & inseparator=runTimeValues%sep, &
        & ignoreLeadingSpaces=.true. )
      nValues = NumStringElements( runtimevalues%lkeys, countEmpty, &
        & runTimeValues%sep )
      ! The first line will be the header
      keysValues(1,1) = 'names'
      keysValues(1,2) = 'Level 2 run-time macro values'
      nValues = nValues + 1
      call outputTable( keysValues(1:nValues, :), border='-', headliner='-' )
    endif
  end subroutine DumpMacros

  ! -----------------------------------  SubstituteRuntimeBoolean  -----
  ! Substitute values for named runtime macros in the original
  subroutine SubstituteRuntimeBoolean ( original, substitute )
    ! Args:
    character(len=*), intent(in)  :: original
    character(len=*), intent(out) :: substitute
    ! Internal variables
    substitute = EvaluateFormula ( original, &
      & runTimeValues%lvalues, runTimeValues%lkeys, runTimeValues%sep )
  end subroutine SubstituteRuntimeboolean

  ! ---------------------------------------  RemoveRuntimeBoolean  -----
  ! Remove the named runtime macros
  subroutine RemoveRuntimeBoolean ( NAME )
  use MLSStringLists, only: GetHashElement, &
    & RemoveHashArray, RemoveHashElement
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
    end if
  end subroutine RemoveRuntimeBoolean

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSL2Options
!=============================================================================

!
! $Log$
! Revision 2.135  2023/05/25 22:26:32  pwagner
! Note that v4 and v5 share a single code base from now on
!
! Revision 2.134  2021/07/22 23:14:53  pwagner
! May print chunknum, phase on warning mesgs
!
! Revision 2.133  2020/04/09 23:19:20  pwagner
! Prevent L2Options from overriding --SkipRetrieval on cmdline or in opts file
!
! Revision 2.132  2020/02/21 21:46:47  pwagner
! Added Default_L2GPFormat_write so we can choose NetCDF4 format
!
! Revision 2.131  2020/02/07 01:16:51  pwagner
! Update versId string and how to re-set it at runtime
!
! Revision 2.130  2019/08/01 23:48:13  pwagner
! Some light Housekeeping
!
! Revision 2.129  2019/07/22 23:22:36  pwagner
! Some light housekeeping, got rid of DumpOptions_old
!
! Revision 2.128  2019/07/09 20:53:22  pwagner
! Use Table ccells to DumpOptions
!
! Revision 2.127  2019/05/15 17:28:39  pwagner
! Warns if opts file is empty
!
! Revision 2.126  2019/05/13 23:34:52  pwagner
! lengthened KeysValues; can usually hold full l2cf path_name now
!
! Revision 2.125  2019/04/18 16:25:40  pwagner
! Allow env lines in opts file (which have no '=' sign)
!
! Revision 2.124  2019/03/08 00:08:55  pwagner
! Various improvements; logging now works properly more reliably
!
! Revision 2.123  2019/02/13 17:33:22  pwagner
! May override L2Options at runtime; added fields GPH_MissingValue and Overriden
!
! Revision 2.122  2019/01/10 21:45:59  pwagner
! SwitchDetail behaviot means no need to discard multiple switches
!
! Revision 2.121  2018/12/07 00:20:26  pwagner
! If cmdline says to skip retrievals, must skip even if l2cf says otherwise
!
! Revision 2.120  2018/10/05 20:49:21  pwagner
! Dont mention switches to be removed unless dbugging
!
! Revision 2.119  2018/09/13 20:20:18  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.118  2018/07/27 23:15:40  pwagner
! Renamed MLSMessage MLSL2Message
!
! Revision 2.117  2018/05/31 18:04:29  pwagner
! The opts file can now override globalattributes
!
! Revision 2.116  2018/05/22 23:07:29  pwagner
! May set InstrumentName by cmdline option
!
! Revision 2.115  2018/04/19 23:44:55  pwagner
! Skip may take /nextChunk flag
!
! Revision 2.114  2018/02/09 01:05:40  pwagner
! Improved comments
!
! Revision 2.113  2017/12/07 01:01:23  vsnyder
! Don't use host-associated variable as a DO index
!
! Revision 2.112  2017/11/30 20:57:10  pwagner
! opts file mechanism may now set OutputOptions, AdvancedOptions, etc.
!
! Revision 2.111  2017/03/24 22:59:23  pwagner
! Made new opts file name CrashIfMsgSays that tells level 2 to crash with walkback if special msg logged
!
! Revision 2.110  2017/01/25 18:05:59  pwagner
! May skip certain Phases named in phasesToSkip cmdline opt
!
! Revision 2.109  2016/11/04 19:30:54  pwagner
! restoreDefaults less complete by default
!
! Revision 2.108  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.107  2016/05/27 00:05:43  pwagner
! Should now correctly process options containing an embedded space
!
! Revision 2.106  2016/03/18 17:57:35  pwagner
! Make certain the L2CF line cited is actal error, not end of section
!
! Revision 2.105  2016/01/12 00:51:35  pwagner
! Repair error in treating 'name=false' option
!
! Revision 2.104  2015/09/24 22:08:34  pwagner
! Added --maxChunkSize option
!
! Revision 2.103  2015/09/17 23:22:46  pwagner
! Now give default value to MaxChunkSize for l2gp DirectWrites
!
! Revision 2.102  2015/07/16 20:51:08  pwagner
! Should quit if cant open optsFile
!
! Revision 2.101  2015/03/05 18:11:04  pwagner
! Some commandline options were being truncated in parseNameValue; fixed
!
! Revision 2.100  2015/03/04 18:33:32  pwagner
! Commandline option --loc also sets hostName
!
! Revision 2.99  2015/02/13 00:19:54  pwagner
! Nay dump macros more nicely as a Table
!
! Revision 2.98  2014/09/11 18:30:19  pwagner
! Removed unused code; corrected parsing of, e.g., f1=true
!
! Revision 2.97  2014/09/05 01:10:03  vsnyder
! Get Error_Unit from intrinsic ISO_Fortran_Env module.  Delete PID stuff.
! Add -E and --stdout error to send output to stderr.
!
! Revision 2.96  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.95  2014/09/02 18:03:12  pwagner
! Correctly distinguish --pid and --pidf (though former may disappear)
!
! Revision 2.94  2014/08/12 23:31:35  pwagner
! commandline options --backg and --pidf noteFile added
!
! Revision 2.93  2014/08/06 23:32:34  vsnyder
! Comment out USE for InitializeMLSFile MLS_Openfile, which are referenced
! only in commented-out code.  Remove USE for MLS_CloseFile, which is not
! referenced.
!
! Revision 2.92  2014/08/05 00:21:29  pwagner
! Add --pId and --uId to set id_strings for slave task
!
! Revision 2.91  2014/08/01 01:46:23  vsnyder
! AllocFile needs SAVE attribute because it has default initialization
!
! Revision 2.90  2014/07/18 23:19:19  pwagner
! Added record length for allocation log
!
! Revision 2.89  2014/06/30 23:27:47  pwagner
! Can log allocations/deallocations to separate file
!
! Revision 2.88  2014/06/25 20:45:44  pwagner
! Show options read from optsFile, if any
!
! Revision 2.87  2014/06/20 20:28:31  pwagner
! Added --set, --setf, and -versId
!
! Revision 2.86  2014/06/16 20:29:05  pwagner
! Updated version vsn id; --recl now affects global setting
!
! Revision 2.85  2014/06/03 22:41:33  pwagner
! Prints phaseName, chunkNum on severe error mesgs
!
! Revision 2.84  2014/04/22 16:33:36  pwagner
! MLSMessage now prints error message as eye-catching banner
!
! Revision 2.83  2014/04/10 00:44:58  pwagner
! Moved more stuff here
!
! Revision 2.82  2014/04/09 00:44:49  vsnyder
! Add --var variable value option.  Cycle instead of returning with
! --nstdout option.  Set print unit immediately at --stdout option.
!
! Revision 2.81  2014/04/07 18:08:57  pwagner
! Stop Writing MLSFile_T attributes by default; they confuse users
!
! Revision 2.80  2014/03/26 17:45:07  pwagner
! added cmdline option --loc to set ProductionLocation
!
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
! Added slavesCleanUpSelves
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
