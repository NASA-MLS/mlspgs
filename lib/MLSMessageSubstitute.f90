! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
module MLSMessageModule         ! Basic messaging for the MLSPGS suite
!==============================================================================

  use MACHINE, only: CRASH_BURN, EXIT_WITH_STATUS, NEVERCRASH
  use MLSCommon, only: MLSFile_T
  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! A low-weight substitute for the full module MLSMessageModule.f90.
  ! which provides low level messaging for the MLSPGS suite.  The main
  ! routine is MLSMessage, which generates log messages as directed by the
  ! user. In the high-fat module, which needs the toolkit,
  ! MLSMessage routine logs a message using the SDPToolkit routine
  ! PGS_SMF_GenerateStatusReport.  This writes a string to the `LogReport'
  ! file (PCF# 10101) in the toolkit.  
  
  ! Here, however, the Toolkit `substitute' just does a simple print.
  
  ! Alternate entries for special circumstances are PVMErrorMessage
  ! (to log a PVM Error)
  ! and ReportTKStatus
  
  ! Another choice to report an error is StopWithErrorMsg
  ! which lets you dump a calling stack you create using MLSMessageCalls
  ! (to report on the severity of a PGS Toolkit return status)
  
  ! Yet another mode is MLSMessage_, useful if you overload MLSMessage with
  ! another module's MLSMessage subroutine accepting extra args
  ! which can then turn around and call MLSMessage_

  ! We dispense with most of the toolkit panoply, needing only the modules
  ! Machine
  ! MLSKinds
  ! MLSCommon

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! MLSMSG_SUCCESS           status returned when all went well
! MLSMSG_DEBUG             should print only if debugging turned on
! MLSMSG_INFO              fyi only
! MLSMSG_WARNING           not fatal, but deserving of attention
! MLSMSG_ERROR             quits after printing
! MLSMSG_CRASH             should give traceback before quitting
! MLSMSG_Severity_to_quit  severity level needed to quit
! MLSMSG_Allocate          mesg prefix for this type of error
! MLSMSG_Fileopen          mesg prefix for this type of error
! MLSMSG_Keyword           mesg prefix for this type of error
! MLSMSG_L1BRead           mesg prefix for this type of error
! MLSMSG_Duplicate         mesg prefix for this type of error
! MLSMSG_DeAllocate        mesg prefix for this type of error
! MLSMSG_PVM               mesg prefix for this type of error
! MLSMessageConfig         configuration controlling where to print, etc.

!     (subroutines and functions)
! MLSMessage               main messaging routine
! MLSMessageSetup          routine interface to change some parts of MLSMessageConfig
! MLSMessageCalls          manage calling stack 
! MLSMessageClose          close MLSMessage log file; but see MLSMessageExit
! MLSMessageExit           recommended way to finish main program
! MLSMessageInternalFile   Returns the complete text that would be printed
! MLSMessageReset          reset flags, counters, etc. during runtime
! PVMErrorMessage          log a PVM error
! ReportTKStatus           converts SDP status to severity, prints if needed
! StopWithErrorMsg         report error msg, dump calling stack, stop

! === (end of toc) ===

! === (start of api) ===
! MLSMessage ( int Severity, char* ModuleNameIn, char* Message, 
!      [char* Advance], [MLSFile_T MLSFile] )
! MLSMessageCalls ( char* command, [char* name] )
! MLSMessageSetup ( [log SuppressDebugs], [int LogFileUnit], [char* Prefix],
!      [log useToolkit], [log CrashOnAnyError] )
! MLSMessageExit ( [int status], [char* farewell] )
! char* MLSMessageInternalFile ( int Severity, char* ModuleNameIn, char* Message, 
!      [char* Advance], [MLSFile_T MLSFile] ) 
! MLSMessageReset ( [int logFileUnit], [log CrashOnAnyError], [log Warnings] )
! PVMErrorMessage ( int INFO, char* PLACE  )
! ReportTKStatus( int status, char* ModuleNameIn, char* Message, 
!      [int Threshold] )
! StopWithErrorMsg ( char* Message, [MLSFile_T MLSFile] )
! === (end of api) ===
  ! ---------------------------------------------------------------------------

  ! Define some low level parameters.  These are used by the calling code to
  ! indicate the severity or otherwise of the messages.

  integer, public, parameter :: MLSMSG_Success = 0
  integer, public, parameter :: MLSMSG_Debug   = MLSMSG_Success + 1
  integer, public, parameter :: MLSMSG_Info    = MLSMSG_Debug + 1
  integer, public, parameter :: MLSMSG_Warning = MLSMSG_Info + 1
  integer, public, parameter :: MLSMSG_Error   = MLSMSG_Warning + 1
  ! Warning--a Crash may not properly close files opened by your run
  ! Use it only for specific debugging where you need a walkback
  ! See also MLSMessageConfig%crashOnAnyError
  integer, public, parameter :: MLSMSG_Crash   = MLSMSG_Error + 1

  ! MLSMSG_Severity_to_quit can be reset in a main program to cause us
  ! to become more lenient (set it higher) or strict (set it lower )
  integer, public            :: MLSMSG_Severity_to_quit = MLSMSG_Error

  private :: SeverityNames
  character (len=*), dimension(MLSMSG_Success:MLSMSG_Crash), parameter :: &
     & SeverityNames = (/&
     & "Success", &
     & "Debug  ", &
     & "Info   ", &
     & "Warning", &
     & "Error  ", &
     & "Crash  " &
     /)

  ! So that we may limit the number of times warnings printed, messagewise
  character(len=*), parameter :: WARNINGSSUPPRESSED = '(No more warnings of this)'
  integer, parameter :: MAXNUMWARNINGS = 40 ! was 80, but most tests < 10
  integer, parameter :: WARNINGMESSLENGTH = 80
  character(len=WARNINGMESSLENGTH), dimension(MAXNUMWARNINGS), save :: &
    &                   warningmessages = ' '
  integer, dimension(MAXNUMWARNINGS), save :: timeswarned = 0
  integer, save :: numwarnings = 0
  ! This set of parameters are simple prefixes for common messages

  character (len=*), public, parameter :: MLSMSG_Allocate = &
     & "Allocation failed: "
  character (len=*), public, parameter :: MLSMSG_Fileopen = &
     & "Failed to open file: "
  character (len=*), public, parameter :: MLSMSG_Keyword = &
     & "Unrecognized configuration file keyword: "
  character (len=*), public, parameter :: MLSMSG_L1BRead = &
     & "Unable to read L1B data item: "
  character (len=*), public, parameter :: MLSMSG_Duplicate = &
     & "There is already an entry with the name "
  character (len=*), public, parameter :: MLSMSG_DeAllocate = &
     & "Deallocation failed: "
  character (len=*), public, parameter :: MLSMSG_PVM = &
     & "PVM Error: "
  ! This datatype describes the configuration of the messaging suite

  integer, private, parameter :: MLSMSG_PrefixLen = 32

   ! May get some of these from MLSLibOptions? 
   ! May get some of these from MLSLibOptions? 
  type, public :: MLSMessageConfig_T
    ! We log messages by toolkit (if useToolkit and UseSDPToolkit are TRUE )
    ! In the following, values would have the effect of adding logged messages:
    ! -2: none added
    ! -1: to stdout
    !  n: to ftn unit n
    integer :: logFileUnit                     = -2
    ! In the following, values would have the effect on identical warnings of:
    ! -1: Print every one without suppression
    !  0: Suppress every one
    !  1: Print every one only once
    integer :: limitWarnings                   = 1000 ! Max number each warning
    integer :: masterTID                       = -1 ! Where to send error msg
    character (len=MLSMSG_PrefixLen) :: prefix = ''   ! Prefix to every msg
    ! Instead of simply calling them Info, you could use something more
    ! informative, like Phase names
    character (len=MLSMSG_PrefixLen) :: Info   = 'Info' ! What to call Info
    logical :: suppressDebugs                  = .false.
    logical :: useToolkit                      = .true.
    logical :: CrashOnAnyError                 = .false. ! See crash warning
    logical :: SendErrMsgToMaster              = .false. ! Whether to send last
    ! last file we were reading/writing if an error occurs and the file
    ! isn't passed in the call statement
    type(MLSFile_T) :: MLSFile ! which file we were reading/writing last
  end type MLSMessageConfig_T

  ! This variable describes the configuration

  type (MLSMessageConfig_T), public, save :: MLSMessageConfig
  
  ! The following can be used to help trace a sequence of calls that led
  ! to an error; it will be dumped (if non-blank) on calling StopWithErrorMsg
  ! You may push a name onto it, pop a name off, or reset it by
  ! appropriate commands sent with subroutine MLSMessageCalls
  
  ! Note the following limitations:
  ! Each name must be shorter than 33 chars, may not contain a '?' character 
  ! (such a name will be split at the '?')
  ! strung together, all the names must not exceed 2048 characters in length
  character(len=2048), public, save :: MLSCallStack = ' '
  
  ! Public procedures
  public :: MLSMessage, MLSMessage_, MLSMessageSetup, MLSMessageClose
  public :: MLSMessageExit, MLSMessageInternalFile
  public :: MLSMessageReset, PVMErrorMessage
  public :: ReportTKStatus
  public :: MLSMessageCalls, StopWithErrorMsg

  interface MLSMessage
    module procedure MLSMessage_
  end interface

contains

  ! ------------------------------------------------  MLSMessage_  -----

  ! This first routine is the main `messaging' code.

  subroutine MLSMessage_( severity, ModuleNameIn,  Message, Advance, MLSFile )
    ! A wraparound subroutine so we can intercept calls 
    ! when the severity is MLSMSG_Error
    ! allowing us to push ModuleNameIn onto the calling stack
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    type(MLSFile_T), intent(in), optional :: MLSFile
    ! Executable
    if ( severity /= MLSMSG_Error ) then
      ! For warnings and so on, just pass args to MLSMessageStd
      call MLSMessageStd( severity, ModuleNameIn,  Message, Advance )
      return
    endif
    call MLSMessageCalls( 'push', constantName=ModuleNameIn )
    call StopWithErrorMsg( Message, MLSFile )
  end subroutine MLSMessage_

  ! --------------------------------------------  MLSMessageCalls  -----

  ! Manage the calling stack MLSCallStack
  ! It will be dumped on calling StopWithErrorMsg
  ! possible commands are 
  ! 'push'    push a new name onto MLSCallStack
  ! 'pop'     pop the last name off
  ! 'rpush'   push a new name underneath
  ! 'rpop'    pop the first name out from underneath
  ! 'clear'   clear the stack
  ! 'print'   print its contents as a single line
  ! 'dump'    print a walkback, one name per line, top to bottom
  ! 'rdump'   print a walkback, one name per line, bottom to top
  
  ! Every command except 'print' and 'clear' change the stack
  ! '[r]push' requires name as an input arg
  ! '[r]pop', and '[r]dump' produce name as an output arg

  ! If name is omitted, it will be lost or assumed blank, as appropriate

  subroutine MLSMessageCalls ( command, name, constantName )
    ! Args
    character(len=*), intent(in)              :: command
    ! Because name is (inout) you cannot call this with a constant
    ! so if you wish to use a constant use constantName instead
    character(len=*), optional, intent(inout) :: name
    character(len=*), optional, intent(in)    :: constantName
    ! Internal variables
    character(len=1), parameter :: comma = '?' ! ','
    integer :: ind
    integer :: m
    character(len=32) :: myName
    ! Executable
    myName = ' '
    if ( index( command, 'push' ) > 0 ) then
      if ( present(name) ) myName = name
      if ( present(constantName) ) myName = constantName
      myName = snipRCSFrom ( myName )
    endif
    ind = index( MLSCallStack, comma, back=.true. )
    m = len_trim(MLSCallStack)
    select case( command )
    case ( 'push' )
      if ( m < 1 ) then
        MLSCallStack = myName
      else
        MLSCallStack = trim(myName) // comma // MLSCallStack
      endif
    case ( 'pop' )
      ind = index( MLSCallStack, comma )
      if ( ind < 1 ) then
        myName = MLSCallStack
        MLSCallStack = ' '
      else
        myName = MLSCallStack( 1 : ind-1 )
        MLSCallStack = MLSCallStack( ind+1 : )
      endif
      if ( present(name) ) name = myName
    case ( 'rpush' )
      if ( m < 1 ) then
        MLSCallStack = myName
      else
        MLSCallStack = trim(MLSCallStack) // comma // myName
      endif
    case ( 'rpop' )
      if ( ind < 1 ) then
        myName = MLSCallStack
        MLSCallStack = ' '
      else
        myName = MLSCallStack( ind+1 : m )
        MLSCallStack = MLSCallStack( 1 : ind-1 )
      endif
      if ( present(name) ) name = myName
    case ( 'clear' )
      MLSCallStack = ' '
    case ( 'rdump' )
      call MLSMessage ( MLSMSG_Info, ModuleName, 'Calling stack (top-down)' )
      do
        if ( m <= ind ) exit
        myName = MLSCallStack( ind+1 : m )
        call MLSMessage ( MLSMSG_Info, ModuleName, myName )
        m = ind - 1
        ind = index( MLSCallStack(1:max(m,1)), comma, back=.true. )
      end do
      if ( present(name) ) name = myName
    case ( 'dump' )
      m = 0
      ind = index( MLSCallStack, comma )
      m =   index( MLSCallStack(ind+1:), comma )
      call MLSMessage ( MLSMSG_Info, ModuleName, 'Calling stack (bottom-up)' )
      do
        ind = m + 1
        m =   ind - 1 + index( MLSCallStack(ind:), comma )
        if ( m < ind + 1 ) m = max( len_trim(MLSCallStack), ind ) + 1
        myName = MLSCallStack( ind : m - 1 )
        call MLSMessage ( MLSMSG_Info, ModuleName, myName )
        if ( present(name) .and. ind == 1 ) name = myName
        if ( m > len_trim(MLSCallStack) - 1 ) exit
      end do
    case ( 'print' )
      if ( len_trim(MLSCallStack) > 0 ) &
        & call MLSMessage ( MLSMSG_Info, ModuleName, trim(MLSCallStack) )
    end select
  end subroutine MLSMessageCalls

  !-----------------------------------------  MLSMessageInternalFile  -----
  function MLSMessageInternalFile( Severity, ModuleNameIn, Message ) result(line)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=512)           :: Line   ! Line to output, should be long enough
    integer                       :: Line_len
    ! Internal variables
    ! Executable
    Line_len = 0
    line = ' '
    call assembleFullLine( Severity, ModuleNameIn, Message, line, line_len, &
      & .false. )
  end function MLSMessageInternalFile

  ! --------------------------------------------  MLSMessageSetup  -----

  ! This routine sets up the MLSMessage suite.  The defaults are of course
  ! sensible, but the user may wish to change things.

  subroutine MLSMessageSetup ( SuppressDebugs, LogFileUnit, Prefix, useToolkit, &
    & CrashOnAnyError  )

    ! Dummy arguments
    logical, optional, intent(in) :: SuppressDebugs
    integer, optional, intent(in) :: LogFileUnit
    character (len=*), optional, intent(in) :: Prefix
    logical, optional, intent(in) :: useToolkit
    logical, optional, intent(in) :: CrashOnAnyError

    ! Local variables

    ! Executable code

    if ( present(suppressDebugs) ) &
      & MLSMessageConfig%suppressDebugs=suppressDebugs

    if ( present(prefix) ) &
      & MLSMessageConfig%prefix=prefix

    if ( present(logFileUnit) ) then
      if ( MLSMessageConfig%logFileUnit /= -1 ) call MLSMessage ( &
        & MLSMSG_Warning, ModuleName,"Already writing to a log file" )
      MLSMessageConfig%logFileUnit = logFileUnit
    end if

    if ( present(useToolkit) ) &
      & MLSMessageConfig%useToolkit=useToolkit
    if ( present(CrashOnAnyError) ) &
      & MLSMessageConfig%CrashOnAnyError=CrashOnAnyError

  end subroutine MLSMessageSetup

  ! --------------------------------------------  MLSMessageClose  -----

  ! This routine simply closes the MLSMessage log file if there is one.

  subroutine MLSMessageClose
    ! Executable code
    MLSMessageConfig%logFileUnit=-1
  end subroutine MLSMessageClose

  ! --------------------------------------------  MLSMessageExit  -----

  ! This routine (optionally) logs farewell, advances
  ! (hopefully) gracefully ends logging, and exits 
  ! (optionally with status )
  ! if farewell present, and non-blank, logs it
  ! if farewell present,  but blank, logs default message
  ! if farewell absent, does not log
  subroutine MLSMessageExit ( status, farewell )
  integer, optional, intent(in) :: STATUS
  character(LEN=*), optional, intent(in) :: FAREWELL
  CHARACTER(LEN=36) :: mesg

    ! Executable code
    if(present(status)) then
      if(present(farewell)) then
        if(farewell == ' ') then
          write(mesg, '(A29, I2, A1)') 'Exiting with status (', &
          status, ')'
          call MLSMessage ( MLSMSG_Info, ModuleName, mesg, advance='y' )
        else
          call MLSMessage ( MLSMSG_Info, ModuleName, farewell, advance='y' )
        end if
      end if
      call MLSMessageClose
      call exit_with_status ( status  )
    else
      if(present(farewell)) then
        if(farewell == ' ') then
          mesg='Exiting normally with "stop"'
          call MLSMessage ( MLSMSG_Info, ModuleName, mesg, advance='y' )
        else
          call MLSMessage ( MLSMSG_Info, ModuleName, farewell, advance='y' )
        end if
      end if
      call MLSMessageClose
      stop
    end if
  end subroutine MLSMessageExit

  ! --------------------------------------------  MLSMessageReset  -----

  ! This routine allows you to reset flags, counters, etc. during runtime

  subroutine MLSMessageReset ( logFileUnit, CrashOnAnyError, Warnings )
    ! Args
    integer, intent(in), optional :: logFileUnit
    logical, intent(in), optional :: CrashOnAnyError
    logical, intent(in), optional :: Warnings
    character(len=6) :: logname
    ! Executable code
    if ( present(logFileUnit) ) then
      if ( logFileUnit /= MLSMessageConfig%logFileUnit ) then
        write(logname, '(i6)') MLSMessageConfig%logFileUnit
        call MLSMessage ( MLSMSG_Info, ModuleName, &
          & 'Closing output on' // logname )
        call MLSMessageClose
        MLSMessageConfig%logFileUnit = logFileUnit
        write(logname, '(i6)') MLSMessageConfig%logFileUnit
        call MLSMessage ( MLSMSG_Info, ModuleName, &
          & 'Opening output on' // logname )
      end if
    end if
    if ( present(CrashOnAnyError) ) MLSMessageConfig%CrashOnAnyError = CrashOnAnyError
    if ( present(Warnings) ) then
      numwarnings = 0
      timeswarned = 0
      warningmessages = ' '
    end if
  end subroutine MLSMessageReset

  ! --------------------------------------------  PVMERRORMESSAGE  -----
  subroutine PVMErrorMessage ( INFO, PLACE  )
    ! This routine is called to log a PVM error
    integer, intent(in) :: INFO
    character (LEN=*) :: PLACE

    character (LEN=132) :: LINE

    write (line, * ) info
    call MLSMessage ( MLSMSG_Error, Place, MLSMSG_PVM // &
      & ' Info='//trim(adjustl(line)))
  end subroutine PVMErrorMessage

  ! --------------------------------------------  ReportTKStatus  -----

  ! Report on severity represented by value returned by toolkit function
  ! Some functions may return non-zero values that are mere warnings or notices
  ! This routine converts these values to a severity level
  ! By default, it then prints only if severity is MLSMSG_Error or higher
  ! but by setting threshold this can be changed
  ! E.g., setting threshold to 0 will print every status, even success

  subroutine ReportTKStatus( status, ModuleNameIn, Message, Threshold )
    ! Dummy arguments
    integer, intent(in) :: status ! e.g. PGS_TD_NOLEAPSECFILE
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in), optional :: Threshold ! Min severity to log message

    ! Internal variables
    integer :: levelmask
    integer :: myThreshold
    integer :: severity

    ! Executable code
    call MLSMessage( MLSMSG_Info, ModuleNameIn, Message )
  end subroutine ReportTKStatus

  ! ------------ StopWithErrorMsg ------------
  subroutine StopWithErrorMsg ( Message, MLSFile )
    ! Print Message, dump calling stack (if any) and stop
    character (len=*), intent(in) :: Message ! Line of text
    type(MLSFile_T), intent(in), optional :: MLSFile
    ! Internal variables
    character(len=32) :: name
    ! Executable
    if ( len_trim(MLSCallStack) > 0 ) call MLSMessageCalls( 'dump', name )
    if ( len_trim(name) < 1 ) name = ModuleName
    call MLSMessageStd( MLSMSG_Error, name, Message, MLSFile=MLSFile )
  end subroutine StopWithErrorMsg

  ! Private procedures
  !-----------------------------------------  accessDFACCToStr  -----
  function accessDFACCToStr ( dfacc ) result(str)

    ! This routine converts an hdf access type
    ! like DFACC_RDONLY into a string like 'rdonly'
    ! If access type is unrecognized, returns 'unknown'
    ! Args
    integer, intent(in)           :: dfacc
    character(len=8)              :: str
    ! Executable
    select case (dfacc)
    case default
      str = 'unknown' ! Why not ' '? Or '?'
    end select
  end function accessDFACCToStr

  ! --------------------------------------------  dumpFile  -----
  subroutine dumpFile ( MLSFile  )
    ! Show everything about it
    type(MLSFile_T) :: MLSFile
    integer, parameter :: SEVERE = MLSMSG_Error
    ! Executable code
    call printitout ( 'MLS File Info: ', MLSMSG_Error )                                  
    call dump ( '(name) ', charValue=trim(MLSFile%Name))                                  
    call dump ( 'short name ', charValue=trim(MLSFile%shortName))                                  
    call dump ( '    Type (int)   : ', MLSFile%Type)
    call dump ( '    Type         : ', charValue=trim(MLSFile%TypeStr))
    call dump ( '    Access       : ', charValue=trim(accessDFACCToStr(MLSFile%access)))
    call dump ( '    content      : ', charValue=trim(MLSFile%content))
    call dump ( '    last Operatn : ', charValue=trim(MLSFile%lastOperation))
    call dump ( '    File ID      : ', MLSFile%FileId%f_id)
    call dump ( '    Group ID     : ', MLSFile%FileId%grp_id)
    call dump ( '    DataSet ID   : ', MLSFile%FileId%sd_id)
    call dump ( '    PCF ID       : ', MLSFile%PCFId)
    call dump ( '    PCF Range    : ', MLSFile%PCFidRange%Bottom)
    call dump ( '                 : ', MLSFile%PCFidRange%Top)
    call dump ( '    hdf version  : ', MLSFile%HDFVersion)
    call dump ( '    record length: ', MLSFile%recordLength)
    call dump ( '    Open?        : ', logValue= MLSFile%StillOpen)
    call dump ( '    error code   : ', MLSFile%errorCode)
  end subroutine dumpFile

  ! -------------------- LastGasp -------------------
  ! In the full module this is called when:
  ! We're a slave and we're about to expire
  ! Before we do, however, try to tell the master why
  subroutine LastGasp ( ModuleNameIn, Message )
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
  end subroutine LastGasp

  !-----------------------------------------  assembleFullLine  -----
  subroutine assembleFullLine( Severity, ModuleNameIn, Message, &
    & line, line_len, nosubsequentwarnings )
    integer, intent(in)           :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character(len=*) ::              Line
    integer                       :: line_len
    logical, intent(in)           :: nosubsequentwarnings
    ! Assemble a full message line

    if ( line_len == 0 ) then
      if ( severity == MLSMSG_Info ) then
        line = MLSMessageConfig%Info
      elseif ( severity > MLSMSG_Success-1 .and. severity < MLSMSG_Crash+1 ) then
        line = SeverityNames(severity)
        if ( MLSMessageConfig%Info /= SeverityNames(MLSMSG_Info) ) then
          line = SeverityNames(severity) // ' ' // &
            & MLSMessageConfig%Info
        endif
      else
        line = 'Unknown'
      end if
      line_len = len_trim(line)
      line(line_len+1:line_len+2) = ' ('
      line(line_len+3:) = snipRCSFrom ( moduleNameIn )
      line_len = len_trim(line) + 3
      line(line_len-2:line_len-1) = '):'
    end if
    if ( nosubsequentwarnings ) then
      line(line_len+1:) = WARNINGSSUPPRESSED // message
      line_len = line_len + len(WARNINGSSUPPRESSED) + len_trim(message)
    else
      line(line_len+1:) = message
      line_len = line_len + len(message) ! Not len-trim, so we can get
      ! trailing blanks into a part of a message.  If there are trailing
      ! blanks remaining when my_adv is true, they'll be trimmed off.
    end if
  end subroutine assembleFullLine

  !-----------------------------------------  level2severity  -----
  function level2severity ( level ) result(severity)

    ! In the full module:
    ! This routine converts a toolkit levelmask to an mls severity
    ! Args
    integer, intent(in)           :: level
    integer :: severity
    severity = MLSMSG_Info
  end function level2severity

  ! --------------------------------------------  dump  -----
  subroutine dump ( name, intValue, charValue, logValue  )
    ! In any way we're asked
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: intValue
    character(len=*), intent(in), optional :: charValue
    logical, intent(in), optional :: logValue
    !
    character(len=132) :: line
    character(len=32) :: value
    value = ''
    if ( present(intValue) ) then
      write(value, '(i10)') intValue
    elseif ( present(logValue) ) then
      write(value, '(l10)') logValue
    endif
    if ( present(charValue) ) then
      line = trim(name) // ' : ' // trim(charvalue)
    else
      line = trim(name) // ' : ' // trim(value)
    endif
    call printitout(trim(line), MLSMSG_Error)
  end subroutine dump

  subroutine MLSMessageStd ( Severity, ModuleNameIn, Message, Advance, MLSFile )

    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    type(MLSFile_T), intent(in), optional :: MLSFile

    ! Local variables
    character (len=512), save :: Line   ! Line to output, should be long enough
    integer, save :: Line_len=0         ! Number of saved characters in line.
    !                                     If nonzero, do not insert prefix.
    integer :: msgLength                  
    logical :: My_adv
    logical :: nosubsequentwarnings
    logical :: newwarning
    integer :: warning_index

    ! Executable code

    my_adv = .true.
    if ( present(advance) ) &
      & my_adv = advance(1:1) /= 'n' .and. advance(1:1) /= 'N'
    ! This is the smaller of the actual length and what we can check for repeats
    msgLength = min(len(message), WARNINGMESSLENGTH)

    ! Here's where we suppress warning messages beyond a limit for each
    nosubsequentwarnings = .false.
    if ( severity == MLSMSG_Warning .and. MLSMessageConfig%limitWarnings > -1 &
      & .and. numwarnings <= MAXNUMWARNINGS .and. message /= ' ' ) then
      ! See if we have seen this message before
      ! newwarning = .not. any(warningmessages == trim(message))
      newwarning = .true.
      do warning_index = 1, numwarnings
        newwarning = newwarning .and. &
          & ( warningmessages(warning_index) /= trim(message(1:msgLength)) )
      enddo
      if ( newwarning .and. numwarnings >= MAXNUMWARNINGS ) then
      else if ( newwarning .or. &
        & numwarnings < 1 ) then
        numwarnings = numwarnings + 1
        warningmessages(numwarnings) = message
        timeswarned(numwarnings) = timeswarned(numwarnings) + 1
        if ( timeswarned(numwarnings) > MLSMessageConfig%limitWarnings ) return
        timeswarned(numwarnings) = min(timeswarned(numwarnings) + 1, &
          & MLSMessageConfig%limitWarnings + 1 )
        nosubsequentwarnings = &
          & (timeswarned(numwarnings) >= MLSMessageConfig%limitWarnings)
      else
        do warning_index = 1, numwarnings
          if ( warningmessages(warning_index) == message(1:msgLength) ) exit
        end do
        if ( warning_index > numwarnings ) return
        if ( timeswarned(warning_index) > MLSMessageConfig%limitWarnings ) return
        timeswarned(warning_index) = min(timeswarned(warning_index) + 1, &
          & MLSMessageConfig%limitWarnings + 1 )
        nosubsequentwarnings = &
          & (timeswarned(warning_index) >= MLSMessageConfig%limitWarnings)
      end if
    end if
      
    if ( (.not. MLSMessageConfig%suppressDebugs).OR. &
         & (severity /= MLSMSG_Debug) ) then
       
      call assembleFullLine( Severity, ModuleNameIn, Message, line, line_len, &
        & nosubsequentwarnings )

       ! Log the message using the toolkit routine
       ! (or its substitute )
       ! if either using toolkit or severity is sufficient to
       ! quit (which means we might have been called directly
       ! rather than from output module )

       if ( my_adv ) then
         call printitout( line, severity, line_len )
         line_len = 0
         line = ' '
       end if

    end if

    ! Now if it's an error, and the message is complete, then try to close
    ! log file if any and quit (or crash)

    if ( my_adv .and. severity >= MLSMSG_Severity_to_quit ) then
      ! Here's a chance to dump facts about last file we were reading/writing
      if ( present(MLSFile) ) then
        call dumpFile(MLSFile)
      elseif( MLSMessageConfig%MLSFile%name /= ' ' ) then
        call dumpFile( MLSMessageConfig%MLSFile )
      endif
      if ( MLSMessageConfig%SendErrMsgToMaster .and. &
        & MLSMessageConfig%masterTID > 0 ) call LastGasp(ModulenameIn, Message )
      if ( MLSMessageConfig%logFileUnit > 0 ) &
        & close ( MLSMessageConfig%logFileUnit  )
      if ( severity >= MLSMSG_Crash .or. MLSMessageConfig%CrashOnAnyError ) then
        NEVERCRASH = .false.
        call crash_burn
      endif
      call exit_with_status ( 1  )
    end if
  end subroutine MLSMessageStd

  ! --------------------------------------------  PRINTITOUT  -----
  subroutine PRINTITOUT ( LINE, SEVERITY, LINE_LEN  )
    ! In any way we're asked
    character(len=*), intent(in) :: LINE
    integer, intent(in) :: SEVERITY
    integer, optional, intent(in) :: LINE_LEN
    ! Now, if we're also logging to a file then write to that too.

    select case ( MLSMessageConfig%logFileUnit  )
    case ( 0 :  )
      write ( UNIT=max(MLSMessageConfig%logFileUnit,1), FMT=* ) TRIM(line)
    case ( -1  )
      write ( UNIT=*, FMT=* ) TRIM(line)
    case default
    end select

  end subroutine PRINTITOUT

  function snipRCSFrom ( with ) result ( without )
    ! Trim nonsense involving RCS system from input "with"
    ! (if present)
    ! Args
    character(len=*), intent(in) :: with
    character(len=len(with))     :: without
      if ( with(1:1) == '$' ) then
      ! The with is <dollar>RCSFile: <filename>,v <dollar>
        without = with(11:(LEN_TRIM(with)-8))
      else
        without = with
      end if
  end function snipRCSFrom
  
!=======================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSMessageModule
!=======================================================================

!
! $Log$
! Revision 2.3  2007/08/13 17:11:07  pwagner
! Implement MLSCallStack for printing walkback
!
! Revision 2.2  2007/01/23 17:18:04  pwagner
! Fixed an obvious bug; now compiles successfully
!
! Revision 2.1  2007/01/12 00:25:26  pwagner
! First commit
!
