! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.


module PrintIt_m

  ! use ISO_FORTRAN_ENV, only: ERROR_UNIT, OUTPUT_UNIT
  use Machine, only: Crash_Burn, Exit_With_Status, NeverCrash
  use MLSCommon, only: MLSFile_T
  use SDPToolkit, only: UseSDPToolkit, Pgs_Smf_GenerateStatusReport

  implicit none
  private

  public :: assembleFullLine, get_config, logUnitName, printItOut
  public :: set_config, snipRCSFrom

  ! These apply if we don't log messages to a Fortran unit number
  ! other than Error_Unit or Output_Unit
  integer, parameter, public :: InvalidLogUnit      = 0 ! max(0,stdoutLogUnit+1)
  integer, parameter, public :: StdoutLogUnit       = InvalidLogUnit - 1 ! Output_Unit
  integer, parameter, public :: DefaultLogUnit      = StdoutLogUnit - 1 ! Error_Unit

  integer, parameter, public :: PGS_S_SUCCESS = 0
  integer, parameter, public :: PrefixLen = 32

  ! May get some of these from MLSCommon? 
  ! Define some low level parameters.  These are used by the calling code to
  ! indicate the severity or otherwise of the messages.
  ! Normally, we treat any severity of Error or worse as reason to stop.
  ! Any Warning is worth recording, and suppressed when too numerous.
  ! Info may be customized to show, phase name, chunk number, etc.
  ! Be advised, Crash may not properly close files opened by your run.
  ! Use it only for specific debugging where you need a walkback.
  ! See also MLSMessageConfig%crashOnAnyError

  integer, public, parameter :: MLSMSG_Success     = PGS_S_SUCCESS ! == 0
  integer, public, parameter :: MLSMSG_Debug       = MLSMSG_Success + 1
  integer, public, parameter :: MLSMSG_Info        = MLSMSG_Debug + 1
  integer, public, parameter :: MLSMSG_Warning     = MLSMSG_Info + 1
  integer, public, parameter :: MLSMSG_Error       = MLSMSG_Warning + 1
  integer, public, parameter :: MLSMSG_Crash       = MLSMSG_Error + 1
  integer, public, parameter :: MLSMSG_TestWarning = MLSMSG_Crash + 1

  character (len=*), public, parameter :: MLSMSG_Allocate = &
     & "Allocation failed: "
  character (len=*), public, parameter :: MLSMSG_DeAllocate = &
     & "Deallocation failed: "
  type, public :: MLSMessageConfig_T
    ! We log messages by toolkit (if useToolkit and UseSDPToolkit are TRUE )
    ! --- ------- ------ ------- ------- ------- -------
    ! In the following, values would have the effect of adding logged messages:
    ! DEFAULTLOGUNIT: none added
    ! STDOUTLOGUNIT:  to stdout
    !  n > 0:         to ftn unit n
    integer :: logFileUnit             = DEFAULTLOGUNIT ! -2
    ! --- ------- ------ ------- ------- ------- -------

    ! In the following, values would have the effect on identical warnings of:
    ! -1: Print every one without suppression
    !  0: Suppress every one
    !  1: Print every one only once
    integer :: limitWarnings           = 1000 ! Max number each warning
    integer :: masterTID               = -1 ! Where to send error msg
    character (len=prefixLen) :: CrashIfMsgSays &
      &                                = ''   ! Crash if any msg has this string
    character (len=prefixLen) :: prefix &
      &                                = ''   ! Prefix to every msg
    ! Instead of showing both module names and severity for every message
    ! you can control thresholds
    ! (1) the severity below which to skip showing module names
    integer :: skipModuleNamesThr      = MLSMSG_Success ! Always show module
    ! (2) the severity below which to skip showing severity
    integer :: skipSeverityThr         = MLSMSG_Success ! Always show severity
    ! (3) the severity below which to skip messages entirely
    integer :: skipMessageThr          = MLSMSG_Success ! Always show messages
    ! (4) simply skip every debug
    logical :: suppressDebugs          = .false.

    ! Instead of simply calling them Info, you could use something more
    ! informative, like Phase names
    character (len=prefixLen) :: Info   &
      &                                = ''      ! What (else) to call Info
    ! Anything else you would like to prefix Warning or Error messages with,
    ! like Phase names?
    character (len=prefixLen) :: Warning   &
      &                                = ''      ! What (else) to call Info
    logical :: useToolkit              = .true.
    integer :: MaxModuleNameLength     = 32      ! Abbreviate longer name
    integer :: MaxSeverityNameLength   = 8       ! Abbreviate longer severity
    logical :: CrashOnAnyError         = .false. ! See crash warning
    logical :: SendErrMsgToMaster      = .false. ! send last gasp to master?
    logical :: ShowCumulativeSeverity  = .false. ! print severity_so_far?
    logical :: StackTrace              = .false. ! Trace via MLSMessageCalls?
    ! Track the last file we were reading/writing if an error occurs and
    ! that file isn't passed in the call statement
    type(MLSFile_T) :: MLSFile ! = MLSFile_T() (crashes under intel v12)

    logical :: AsciifyMessages = .true.
    integer :: Severity_To_Quit = 0
    logical :: UseDefaultFormatStdout = .false.
    
    ! Temporarily skip using Toolkit (if we ever were)
    ! logical :: adHocPrintToStdout = .false.
  end type MLSMessageConfig_T

  ! This variable describes the configuration
  type (MLSMessageConfig_T), public, save :: MLSMESSAGECONFIG
 
  ! MLSMSG_Severity_to_* can be reset in a main program to cause us
  ! to become more lenient (set it higher) or strict (set it lower )
  integer, public, save      :: MLSMSG_Severity_to_quit     = MLSMSG_Error
  integer, public, save      :: MLSMSG_Severity_to_walkback = MLSMSG_Error
  integer, public, save      :: MLSMSG_Severity_so_far      = PGS_S_SUCCESS

  private :: SeverityNames
  character (len=*), dimension(MLSMSG_Success:MLSMSG_Crash), parameter :: &
     & SeverityNames = (/&
     & "Success", &
     & "Debug  ", &
     & "Info   ", &
     & "Warning", &
     & "Error  ", &
     & "Crash  "  &
     /)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  !--------------------------------------------  AssembleFullLine  -----
  ! Assemble the full line out of 
  ! (1) A severity level
  ! (2) The module name
  ! (3) Whatever message we're asked to repeat (%Info? %Warning?)
  subroutine AssembleFullLine( Severity, ModuleNameIn, Message, &
    & line, line_len )
    integer, intent(in)           :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character(len=*)              :: Line
    integer                       :: line_len
    ! Assemble a full message line

    if ( severity < MLSMessageConfig%skipMessageThr ) then
      return
    elseif ( line_len == 0 ) then
      if ( severity < MLSMessageConfig%skipSeverityThr ) then
        line = ' '
      elseif ( severity > MLSMSG_Success-1 .and. severity < MLSMSG_Crash+1 ) then
        line = trim(SeverityNamesFun(severity))
        if ( MLSMessageConfig%ShowCumulativeSeverity .and. len_trim(line) > 0 ) then
          line = trim(line) // ':' // &
            & trim(SeverityNamesFun(MLSMSG_Severity_so_far))
        endif
        ! Do we prefix with MLSMessageConfig%Info
        if ( severity == MLSMSG_Info .and. &
          & len_trim(MLSMessageConfig%Info) > 0 ) &
          & line = trim(line) // ' ' // MLSMessageConfig%Info
        ! Do we prefix with MLSMessageConfig%Warning
        if ( severity > MLSMSG_Info .and. &
          & len_trim(MLSMessageConfig%Warning) > 0 ) &
          & line = trim(line) // ' ' // MLSMessageConfig%Warning
      else
        line = 'Unknown'
      end if
      line_len = len_trim(line)
      if ( severity >= MLSMessageConfig%skipModuleNamesThr .and. &
        & len_trim(snipRCSFrom ( moduleNameIn )) > 0 ) then
        line(line_len+1:line_len+2) = ' ('
        line(line_len+3:) = snipRCSFrom ( moduleNameIn )
        line_len = len_trim(line) + 3
        line(line_len-2:line_len-1) = '):'
      elseif ( severity >= MLSMessageConfig%skipSeverityThr ) then
        line_len = len_trim(line) + 1
        line(line_len:line_len) = ':'
      end if
    end if
    ! Make sure that we don't start prefixing with a blank
    if ( any( severity >= &
      & (/MLSMessageConfig%skipModuleNamesThr, MLSMessageConfig%skipSeverityThr/) &
      & ) .and. line(1:1) == ' ' ) then
      line_len = line_len - 1
      line(1:line_len) = line(2:line_len+1)
    endif
    line(line_len+1:) = message
    line_len = line_len + len(message) ! Not len_trim, so we can get
    ! trailing blanks into a part of a message.  If there are trailing
    ! blanks remaining when my_adv is true, they'll be trimmed off.
  end subroutine AssembleFullLine

  ! -------------------------------------------------  Get_Config  -----
  subroutine Get_Config ( Asciify, LogFileUnit, Prefix, &
    & Severity_to_Quit, UseDefaultFormatStdout, UseToolkit )
    logical, intent(out), optional :: Asciify, UseDefaultFormatStdout, UseToolkit
    integer, intent(out), optional :: LogFileUnit, Severity_to_Quit
    character(len=*), intent(out), optional :: Prefix
    if ( present(asciify) ) asciify = MLSMessageConfig%asciifyMessages 
    if ( present(logFileUnit) ) logFileUnit = MLSMessageConfig%logFileUnit 
    if ( present(prefix) ) prefix = MLSMessageConfig%prefix
    if ( present(severity_to_quit) ) severity_to_quit = MLSMessageConfig%severity_to_quit 
    if ( present(useDefaultFormatStdout) ) useDefaultFormatStdout = MLSMessageConfig%useDefaultFormatStdout
    if ( present(useToolkit) ) useToolkit = MLSMessageConfig%useToolkit
  end subroutine Get_Config

  ! ------------------------------------------------  LogUnitName  -----
  function LogUnitName ( LogUnit ) result( name )
    ! Return an appropriate name for the LogUnit number
    ! Args
    integer, intent(in) :: LogUnit
    character(len=12) :: name
    ! Executable
    select case ( LogUnit )
    case ( stdoutLogUnit )
      name = 'stdout'
    case ( defaultLogUnit )
      name = 'mls LogUnit'
    case ( invalidLogUnit )
      name = 'invalid'
    case default ! > 0
      name = 'Fortran unit'
    end select
  end function LogUnitName

  ! -------------------------------------------------  PrintItOut  -----
  subroutine PrintItOut ( inLine, severity, &
    & line_len, noPrefix, exitStatus, noExit  )
    ! In any way we're asked .. print inLine
    ! After that, maybe exit with status or worse
    ! Args
    character(len=*), intent(in)  :: inLine      ! What to print
    integer, intent(in)           :: severity    ! Tell me Doc, how bad is it?
    integer, optional, intent(in) :: line_len    ! If different from (len(inLine)
    logical, optional, intent(in) :: noPrefix    ! Don't add any prefix
    integer, intent(in), optional :: exitStatus  ! Exit with this status
    logical, optional, intent(in) :: noExit      ! No, just return no matter what
    ! Local variables
    character(len=len(inline)) :: Line
    logical :: log_it
    integer :: loggedLength
    character(len=len(inline)+len(MLSMessageConfig%prefix)) :: loggedLine
    integer :: ioerror
    integer :: maxLineLength
    logical :: myNoExit
    logical :: myNoPrefix
    logical, parameter :: DEEBUG = .false.
    ! Executable
    if ( MLSMessageConfig%AsciifyMessages ) then
      line = asciify(inLine)
    else
      line = inLine
    end if
    loggedLength = len_trim(line)
    if ( present(line_len) ) loggedLength = max( line_len, loggedLength )
    myNoExit = .false.
    if ( present(noExit) ) myNoExit = noExit
    myNoPrefix = .false.
    if ( present(noPrefix) ) myNoPrefix = noPrefix
    loggedLine = line
    if ( trim(MLSMessageConfig%prefix) /= ' ' .and. .not. myNoPrefix ) then
      loggedLength = loggedLength + len_trim(MLSMessageConfig%prefix)
      loggedLine = trim(MLSMessageConfig%prefix) // &
           & trim(line)
    end if
    maxLineLength = min( loggedLength, len(loggedLine) )
    log_it = (MLSMessageConfig%useToolkit .and. UseSDPToolkit) .or. &
      & severity >= MLSMessageConfig%severity_to_quit
    if ( DEEBUG .and. log_it ) then
      print *, 'trim(loggedLine) ', trim(loggedLine)
      print *, 'maxLineLength ', maxLineLength
    endif
    if( log_it .and. maxLineLength > 0 .and. MLSMessageConfig%useToolkit ) then
      ioerror = PGS_SMF_GenerateStatusReport ( loggedLine(1:maxLineLength) )
    end if

    ! Now, if we're also logging to a file then write to that too.
    select case ( MLSMessageConfig%logFileUnit  )
    case ( StdoutLogUnit  )
      if ( MLSMessageConfig%useDefaultFormatStdout ) then
        write ( unit=*, fmt=* ) trim(line)
      else
        write ( unit=*, fmt='(a)' ) trim(line)
      end if
    case ( defaultLogUnit )
    case default
      write ( UNIT=max(MLSMessageConfig%logFileUnit,1), FMT=* ) trim(line)
    end select
    ! Have we been tasked with something more?
    ! Exit? Crash?
    if ( myNoExit ) return
    if ( severity >= MLSMSG_Crash .or. &
      & (severity > MLSMSG_Warning .and. MLSMessageConfig%CrashOnAnyError) &
      & ) then
      NEVERCRASH = .false.
      call crash_burn
    endif
    if ( present( exitStatus ) )  call exit_with_status ( exitStatus  )
    if ( severity > MLSMSG_Warning ) call exit_with_status ( 1  )

  end subroutine PrintItOut

  ! ----------------------------------------------  ModuleNameFun  -----
  function ModuleNameFun ( moduleName ) result (name)
    ! Return name of module unless asked to abbreviate
    character(len=*), intent(in)    :: moduleName
    character(len=len(ModuleName))  :: name
    name = moduleName
    if ( len_trim(moduleName) < 1 ) name = 'Unknown'
    if ( MLSMessageConfig%MaxModuleNameLength < 1 ) then
      name = ' '
    elseif ( MLSMessageConfig%MaxModuleNameLength < len(ModuleName) ) then
      name = name(1:MLSMessageConfig%MaxModuleNameLength) // ' '
    endif
  end function ModuleNameFun

  ! -------------------------------------------  SeverityNamesFun  -----
  function SeverityNamesFun ( severity ) result (name)
    ! Return name of level corresponding to severity, if recognized
    ! If not recignized, return  'Unknown'
    ! Full name unless asked to abbreviate
    integer, intent(in)                   :: severity
    character(len=len(SeverityNames(1)))  :: name
    if ( severity < MLSMSG_Success .or. severity > MLSMSG_Crash ) then
      name = 'Unknown'
    else
      name = SeverityNames( severity )
    endif
    if ( MLSMessageConfig%MaxSeverityNameLength < 1 ) then
      name = ' '
    elseif ( MLSMessageConfig%MaxSeverityNameLength < len(SeverityNames(1)) ) then
      name = name(1:MLSMessageConfig%MaxSeverityNameLength) // ' '
    endif
  end function SeverityNamesFun

  ! ------------------------------------------------  SnipRCSFrom  -----
  function SnipRCSFrom ( with ) result ( without )
    ! Trim nonsense involving RCS system from input "with"
    ! (if present)
    ! Args
    character(len=*), intent(in) :: with
    character(len=len(with))     :: without
    integer :: secondBuck
      if ( with(1:1) == '$' ) then
      ! The with is <dollar>RCSFile: <filename>,v <dollar><otherstuff>
      ! The without is <filename><otherstuff>
        secondBuck = 1 + index( with(2:),'$')
        if ( secondBuck-3 > 11 .and. secondBuck < len_trim(with) ) then
          ! without = with(11:(LEN_TRIM(with)-8))
          without = with(11:secondBuck-4) // with(secondBuck+1:)
        else
          without = with(11:(LEN_TRIM(with)-8))
        endif
      else
        without = with
      end if
      without = ModuleNameFun( without )
  end function SnipRCSFrom
  
  ! -------------------------------------------------  Set_Config  -----
  subroutine Set_Config ( Asciify, LogFileUnit, Prefix, &
    & Severity_to_Quit, UseDefaultFormatStdout, UseToolkit )
    logical, intent(in), optional :: Asciify, UseDefaultFormatStdout, UseToolkit
    integer, intent(in), optional :: LogFileUnit, Severity_to_Quit
    character(len=*), intent(in), optional :: Prefix
    if ( present(asciify) ) MLSMessageConfig%asciifyMessages  = asciify
    if ( present(logFileUnit) ) MLSMessageConfig%logFileUnit  = logFileUnit
    if ( present(prefix) ) MLSMessageConfig%prefix = prefix
    if ( present(severity_to_quit) ) MLSMessageConfig%severity_to_quit  = severity_to_quit
    if ( present(useToolkit) ) MLSMessageConfig%useToolkit = useToolkit
    if ( present(useDefaultFormatStdout) ) MLSMessageConfig%useDefaultFormatStdout = useDefaultFormatStdout
  end subroutine Set_Config

! *****  Private Procedures     ****************************************

  ! -------------------------------------------------  ASCIIFY  -----
  ! takes input string and replaces any non-printing characters
  ! with a substitute '@'
  function ASCIIFY (STR) result (OUTSTR)
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=len(str))      :: OUTSTR

    !----------Local vars----------!
    integer :: I
    !----------Executable part----------!
    outstr=str
    do i=1, len(str)
      if ( .not. isAscii(str(i:i)) ) outstr(i:i) = '@'
    end do
  end function ASCIIFY

  ! ---------------------------------------------------  isAscii  -----
  elemental function isAscii(arg) result(itIs)
    ! Returns TRUE if arg is in range of printing chars [32,126]
    ! Args
    character(len=1), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    integer, parameter :: pcMin = iachar(' ')
    integer, parameter :: pcMax = iachar('~')
    ! Executable
    itis = iachar(arg) >= pcMin .and. iachar(arg) <= pcMax
  end function isAscii

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PrintIt_m

! $Log$
! Revision 2.9  2017/03/23 16:22:03  pwagner
! Programs may optionally crash when MLSMessage logs fatal string
!
! Revision 2.8  2015/05/06 20:40:49  pwagner
! May prefix Warnings or worse with, e.g., phase and chunk num
!
! Revision 2.7  2013/09/09 18:36:58  pwagner
! Workaround for ifort v12 internal compiler error
!
! Revision 2.6  2013/09/06 20:41:38  pwagner
! Remove config in favor of using MLSMessageConfig
!
! Revision 2.5  2013/08/30 23:11:27  pwagner
! NoExit option prevents unwanted stop
!
! Revision 2.4  2013/08/30 03:56:02  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.3  2013/08/29 19:34:52  pwagner
! Fixed some bugs affecting logging via toolkit
!
! Revision 2.2  2013/08/28 00:35:19  pwagner
! Moved more stuff from MLSMessage down to PrintIt module
!
! Revision 2.1  2013/08/23 02:48:07  vsnyder
! Initial commit
!
