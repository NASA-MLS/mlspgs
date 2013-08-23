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

  use ISO_Fortran_Env, only: Error_Unit, Output_Unit

  implicit none
  private

  public Get_Config, LogUnitName, PrintItOut, Set_Config

  ! These apply if we don't log messages to a Fortran unit number
  ! other than Error_Unit or Output_Unit
  integer, parameter, public :: StdoutLogUnit       = Output_Unit
  integer, parameter, public :: DefaultLogUnit      = Error_Unit
  integer, parameter, public :: InvalidLogUnit      = max(0,stdoutLogUnit+1)

  integer, parameter, public :: PrefixLen = 32

  type, private :: Config_t
    logical :: AsciifyMessages = .true.
    integer :: LogFileUnit = DefaultLogUnit
    integer :: Severity_To_Quit = 0
    logical :: UseDefaultFormatStdout = .false.
    logical :: UseToolkit = .true.
    character(len=prefixLen) :: Prefix
  end type Config_t

  type(config_t), private, save :: Config

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  Get_Config  -----
  subroutine Get_Config ( Asciify, LogFileUnit, Prefix, &
    & Severity_to_Quit, UseDefaultFormatStdout, UseToolkit )
    logical, intent(out), optional :: Asciify, UseDefaultFormatStdout, UseToolkit
    integer, intent(out), optional :: LogFileUnit, Severity_to_Quit
    character(len=*), intent(out), optional :: Prefix
    if ( present(asciify) ) asciify = config%asciifyMessages 
    if ( present(logFileUnit) ) logFileUnit = config%logFileUnit 
    if ( present(prefix) ) prefix = config%prefix
    if ( present(severity_to_quit) ) severity_to_quit = config%severity_to_quit 
    if ( present(useDefaultFormatStdout) ) useDefaultFormatStdout = config%useDefaultFormatStdout
    if ( present(useToolkit) ) useToolkit = config%useToolkit
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
  subroutine PrintItOut ( INLINE, SEVERITY, LINE_LEN, NOPREFIX  )
    ! In any way we're asked
    use SDPToolkit, only: UseSDPToolkit, PGS_SMF_GenerateStatusReport
    ! Args
    character(len=*), intent(in) :: INLINE
    integer, intent(in) :: SEVERITY
    integer, optional, intent(in) :: LINE_LEN
    logical, optional, intent(in) :: NOPREFIX
    ! Local variables
    character(len=len(inline)) :: Line
    logical :: log_it
    integer :: loggedLength
    character(len=len(inline)+len(config%prefix)) :: loggedLine
    integer :: ioerror
    integer :: maxLineLength
    logical :: myNoPrefix
    ! Executable
    if ( config%AsciifyMessages ) then
      line = asciify(inLine)
    else
      line = inLine
    end if
    loggedLength = len_trim(line)
    if ( present(line_len) ) loggedLength = line_len
    myNoPrefix = .false.
    if ( present(noPrefix) ) myNoPrefix = noPrefix
    loggedLine = line
    if ( trim(config%prefix) /= ' ' .and. .not. myNoPrefix ) then
      loggedLength = loggedLength + len_trim(config%prefix)
      loggedLine = trim(config%prefix) // &
           & trim(line)
    end if
    maxLineLength = min( loggedLength, len(loggedLine) )
    log_it = (config%useToolkit .and. UseSDPToolkit) .or. &
           & severity >= config%severity_to_quit
    if( log_it .and. loggedLength > 0 .and. config%useToolkit) then
      ioerror = PGS_SMF_GenerateStatusReport ( loggedLine(1:maxLineLength) )
    end if

    ! Now, if we're also logging to a file then write to that too.
    select case ( config%logFileUnit  )
    case ( StdoutLogUnit  )
      if ( config%useDefaultFormatStdout ) then
        write ( unit=*, fmt=* ) trim(line)
      else
        write ( unit=*, fmt='(a)' ) trim(line)
      end if
    case ( defaultLogUnit )
    case default
      write ( UNIT=max(config%logFileUnit,1), FMT=* ) trim(line)
    end select

  end subroutine PrintItOut

  ! -------------------------------------------------  Set_Config  -----
  subroutine Set_Config ( Asciify, LogFileUnit, Prefix, &
    & Severity_to_Quit, UseDefaultFormatStdout, UseToolkit )
    logical, intent(in), optional :: Asciify, UseDefaultFormatStdout, UseToolkit
    integer, intent(in), optional :: LogFileUnit, Severity_to_Quit
    character(len=*), intent(in), optional :: Prefix
    if ( present(asciify) ) config%asciifyMessages  = asciify
    if ( present(logFileUnit) ) config%logFileUnit  = logFileUnit
    if ( present(prefix) ) config%prefix = prefix
    if ( present(severity_to_quit) ) config%severity_to_quit  = severity_to_quit
    if ( present(useToolkit) ) config%useToolkit = useToolkit
    if ( present(useDefaultFormatStdout) ) config%useDefaultFormatStdout = useDefaultFormatStdout
  end subroutine Set_Config

! *****  Private Procedures     ****************************************

  ! -------------------------------------------------  ASCIIFY  -----
  ! takes input string and replaces any non-printing characters
  ! with corresponding ones in range [32,126]
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
! Revision 2.1  2013/08/23 02:48:07  vsnyder
! Initial commit
!
