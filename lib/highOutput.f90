! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HighOutput

  ! Very high level printing and formatting
  
  ! See also dump_0 and output_m
  
  use Dates_Module, only: BuildCalendar, DaysInMonth, &
    & ReformatDate, ReformatTime, Utc_To_Yyyymmdd
  use Machine, only: Crash_Burn, Exit_With_Status, NeverCrash
  use MLSCommon, only: MLSDebug, MLSVerbose
  use MLSFinds, only: FindFirst
  use MLSStringLists, only: ExpandStringRange, GetStringElement, &
    & List2Array, NCharsInFormat, NumStringElements, SwitchDetail, Wrap
  use MLSStrings, only: Asciify, Lowercase, Ncopies, &
    & Trim_Safe, WriteIntsToChars
  use Output_M, only: Advance_Is_Yes_Or_No, Blanks, GetOutputStatus, &
    & Newline, &
    & Output, Output_ => Output_Char_Nocr, &
    & RestoreOutputSettings => RestoreSettings, &
    & OutputOptions, OutputOptions_T, StampOptions, StampOptions_T, &
    & TimeStampOptions, TimeStampOptions_T, &
    & BothPrUnit, InvalidPrUnit, MSGLogPrUnit, &
    & OutputLines, OutputLinesPrUnit, SetOutputStatus, StdoutPrUnit
  use PrintIt_M, only: AssembleFullLine, Get_Config, &
    & MLSMSG_Crash, MLSMSG_Debug, &
    & MLSMSG_Severity_To_Quit, &
    & MLSMSG_Warning, &
    & PrintItOut, MLSMessageConfig
  use Toggles, only: Switches
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (subroutines and functions)
! AddRow                   add name, value row to a 2-d table of cells
! AddRow_header            add a single line stretched across an entire row
! AddRow_divider           add a row composed of a single, repeated character
! AlignToFit               align printed argument to fit column range
! Banner                   surround message with stars and stripes; e.g.,
!                            *-----------------------------------------------*
!                            *            Your message here                  *
!                            *-----------------------------------------------*
! BeVerbose                Should we do extra printing?
! BlanksToColumn           print blanks [or fill chars] out to specified column
! BlanksToTab              print blanks [or fill chars] out to next tab stop
! Dump                     dump output or stamp options
! Dumpsize                 print a nicely-formatted memory size 
! Dumptabs                 print the current tab stop positions
! GetStamp                 get stamp being added to every output
! HeadLine                 print a line with eye-catching features
!                           e.g., '*-------  Your message here   -------*'
! LetsDebug                Should we do debug printing?
! NumNeedsFormat           return what format is needed to output num
! NumToChars               return what string would be printed by output
! OutputCalendar           output nicely-formatted calendar page
! Output_date_and_time     print nicely formatted date and time
! OutputList               output array as comma-separated list; e.g. '(1,2,..)'
! OutputNamedValue         print nicely formatted name and value
! OutputTable              output 2-d array as cells in table
! Resettabs                restore tab stops to what was in effect at start
! RestoreSettings          restore default settings for output, stamps, tabs
! SetStamp                 set stamp to be automatically printed on every line
! SetTabs                  set tab stops (to be used by tab)
! StartTable               initialize a 2-d table of cells to be output later
! StyledOutput             output a line according to options; e.g. "--Banner"
! Tab                      move to next tab stop
! Timestamp                print argument with a timestamp manually
!                            (both stdout and logged output)
! === (end of toc) ===

! === (start of api) ===
! addRow ( char* name, value )
! addRow_header ( char* name, char alignment )
! addRow_divider ( char char )
! alignToFit ( char* chars, int columnRange(2), char alignment, [int skips] )
! banner ( char* chars, [int columnRange(2)], [char alignment], [int skips], 
!    [int lineLength], [char mode], [char pattern] )
! log BeVerbose ( char* switch, threshold )
! blanksToColumn ( int column, [char fillChar], [char* advance] )
! blanksToTab ( [int tabn], [char* fillChar] )
! Dump ( options )
! DumpSize ( n, [char* advance], [units] )
!       where n can be an int or a real, and 
!       units is a scalar of the same type, if present
! DumpTabs ( [int tabs(:)] )
! getStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! headLine ( char* chars, 
!          [char fillChar], [char* Before], [char* After], 
!          [int columnRange(2)], [char alignment], [int skips] )
! log LetsDebug ( char* switch, threshold )
! char* numNeedsFormat ( value )
! char* numToChars ( value, [char* format] )
! output_date_and_time ( [log date], [log time], [char* from_where], 
!          [char* msg], [char* dateFormat], [char* timeFormat], 
!          [double CPU_seconds], [int wallClock_seconds], [char* advance] )
! outputCalendar ( [char* date], [char* datenote], [char* notes(:)], 
!          [log dontwrap, [log moonPhases] ] )
! outputList ( values(:), [char* sep], [char* delims] )
! outputNamedValue ( char* name, value, [char* advance],
!          [char colon], [char fillChar], [char* Before], [char* After], 
!          [integer tabn], [integer tabc], [integer taba], log dont_stamp],
!          [char* options] )
! outputTable ( [array(:,:)], [char sep], [char border], [int cellWidth],
!          [char interior], [char headliner], [char alignment] )
! resetTabs ( [int tabs(:)] )
! restoreSettings ( [log useToolkit] )
! setStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! setTabs ( [char* Range], [int tabs(:)] )
! startTable
! styledOutput ( char* chars, [char* options] )
! tab ( [int tabn], [char* fillChar] )
! timeStamp ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [char*8 style], [log date] )
! timeStamp ( log value, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [char*8 style], [log date] )
! timeStamp ( int int, [int places], [char* advance],
!          [log fill], [char* format], [char* Before], [char* After],
!          [char*8 style], [log date] )
! === (end of api) ===
!
! Note:
! By calling appropriate functions and procedures you can adjust aspects of
! behavior of output, and others can be changed by setting various
! public global parameters directly
! (in OO-speak they are class-level rather than instance-level)
! Sometimes there is more than one way to accomplish the same thing
! E.g., calling timeStamp or using setStamp before calling output
!
! To understand the codes for dateformat and timeFormat, see the dates_module
! 
! To use this module to build a 2d table of names and values
! (1) Call startTable
! (2) Optionally call addRow_header
! (3) Optionally call addRow_divider
! (4) For each name, value pair
!     (a) call AddRow
!     (b) call AddRow_Divider
! (5) Call OutputTable
  public :: addRow, addRow_divider, addRow_header, alignToFit, &
    & banner, beVerbose, blanksToColumn, blanksToTab, &
    & dump, dumpSize, dumpTabs, getStamp, headLine, &
    & letsDebug, nextColumn, nextTab, numNeedsFormat, numToChars, &
    & output_date_and_time, outputCalendar, outputList, outputTable, &
    & outputAnyNamedValue, outputNamedValue, &
    & resetTabs, restoreSettings, &
    & setStamp, setTabs, startTable, styledOutput, tab, timeStamp

  ! These types must be made public because the class instances are public
  public :: outputOptions_T
  public :: stampOptions_T
  public :: timeStampOptions_T

  interface addRow
    module procedure addRow_character
    module procedure addRow_complex
    module procedure addRow_dbl_array, addRow_double
    module procedure addRow_int_array, addRow_integer
    module procedure addRow_log_array, addRow_logical
    module procedure addRow_sngl_array, addRow_single
  end interface

  interface ALIGNTOFIT
    module procedure aligntofit_chars, aligntofit_double, aligntofit_single
    module procedure aligntofit_integer
  end interface

  interface BANNER
    module procedure banner_chars
    module procedure banner_chararray
  end interface

  interface DUMP
    module procedure DUMPOUTPUTOPTIONS, DUMPSTAMPOPTIONS, DUMPTIMESTAMPOPTIONS
  end interface

  interface DUMPSIZE
    module procedure DUMPSIZE_DOUBLE, DUMPSIZE_INTEGER, DUMPSIZE_REAL
  end interface

  interface GETOPTION
    module procedure getOption_char, getOption_log
  end interface

  interface NUMNEEDSFORMAT
    module procedure numNeedsFormat_double, numNeedsFormat_integer, numNeedsFormat_single
    module procedure numNeedsFormat_complex, numNeedsFormat_dcomplx
  end interface

  interface NUMTOCHARS
    module procedure numtochars_double, numtochars_integer, numtochars_single
  end interface

  interface OUTPUTLIST
    module procedure OUTPUTLIST_INTS, OUTPUTLIST_CHARS
  end interface

  interface OUTPUTNAMEDVALUE
    module procedure output_nvp_character
    module procedure output_nvp_complex
    module procedure output_nvp_dbl_array, output_nvp_double
    module procedure output_nvp_int_array, output_nvp_integer
    module procedure output_nvp_log_array, output_nvp_logical
    module procedure output_nvp_sngl_array, output_nvp_single
  end interface

  interface OUTPUTANYNAMEDVALUE
    module procedure output_nvp_whatever
  end interface

  interface TAB
    module procedure blanksToTab
  end interface
  
  interface TIMESTAMP
    module procedure timestamp_char, timestamp_integer, timestamp_logical
  end interface
  
  ! When Calling OutputNamedValue with character values, should we trim them?
  logical, public                      :: TrimCharacterValues = .true.
  
  ! Used for automatic assembly of a table to be neatly formatted and output
  ! The table holds two columns:
  ! names and values  
  ! It might look something like the following
  ! ------------------------------------------------------------
  ! - names  values                                            -
  ! ------------------------------------------------------------
  ! - true   true                                              -
  ! - false  false                                             -
  ! - count  count                                             -
  ! - l2cf   /users/mmadatya/l2tests/ASMLS/AS-01-009-QTM.l2cf  -
  ! ------------------------------------------------------------
  integer, parameter :: MAXCELLSIZE = 128 ! How many chars can 1 cell hold
  character(len=MAXCELLSIZE), dimension(:,:), pointer :: cellDatabase => null()

  integer, private, parameter :: MAXNUMTABSTOPS = 24
  integer, save, private :: WRAPPASTCOLNUM = 0  ! Don't print beyond (if > 0)
  integer, save, private :: OLDWRAPPASTCOLNUM = 0

  ! These next tab stops can be reset using the procedure setTabs
  ! the default values correspond to range coded '5-120+5'
  ! (read as from 5 to 120 in intervals of 5)
  character(len=*), parameter :: INITTABRANGE = '5-120+5'
  integer, dimension(MAXNUMTABSTOPS), save, private :: TABSTOPS = &
    & (/ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, &
    &   65, 70, 75, 80, 85, 90, 95,100,105,110,115,120 /)

  ! For certain numerical values we will use list directed '*' format
  ! unless optional FORMAT specifier supplied
  double precision, parameter, dimension(3) :: DPREFERDEFAULTFORMAT = &
    & (/ -1.d0, 0.d0, 1.d0 /)  ! For which values to use default format '*'
  character(len=12), private :: sdNeedsFormat = '(1pg14.6)'
  character(len=12), private :: sdNeedsFragment = '(1pg14'

  ! This is the type for configuring how to automatically style 
  ! special output formats; e.g., Banner
  ! Note the effect on the "bars" part of "stars and bars"
  ! of choosing different HeadlineFill or BannerPattern characters:
!------------------------------------------------------------------------------*
!                            Test Banner Pattern: -                            *
!------------------------------------------------------------------------------*
!*******************************************************************************
!                            Test Banner Pattern: *                            *
!*******************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!                            Test Banner Pattern: +                            *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!##############################################################################*
!                            Test Banner Pattern: #                            *
!##############################################################################*
!                                                                              *
!                            Test Banner Pattern: 0                            *
!                                                                              *
!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   *
!                            Test Banner Pattern: 1                            *
!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   *
! . .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .  *
!                            Test Banner Pattern: 2                            *
! . .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .  *
! .  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  . *
!                            Test Banner Pattern: 3                            *
! .  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  . *
! .   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   .  *
!                            Test Banner Pattern: 4                            *
! .   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   .  *
! .. .... .... .... .... .... .... .... .... .... .... .... .... .... .... ..  *
!                            Test Banner Pattern: 5                            *
! .. .... .... .... .... .... .... .... .... .... .... .... .... .... .... ..  *
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   *
!                            Test Banner Pattern: 6                            *
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   *
! - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -  *
!                            Test Banner Pattern: 7                            *
! - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -  *
! -  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  - *
!                            Test Banner Pattern: 8                            *
! -  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  - *
! - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - ..   *
!                            Test Banner Pattern: 9                            *
! - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - ..   *
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   *
!                            Test Banner Pattern: A                            *
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   *
!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   *
!                            Test Banner Pattern: B                            *
!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   *
! = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~  *
!                            Test Banner Pattern: C                            *
! = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~  *

  type StyleOptions_T
    ! Headline
    character(len=1) :: HeadLineAlignment        = 'C'
    character(len=1) :: HeadLineFill             = ' '
    character(len=8) :: HeadLineBefore           = ' '
    character(len=8) :: HeadLineAfter            = ' '
    integer, dimension(2) :: HeadLineColumnrange = (/ 1, 80 /)
    integer          :: HeadLineSkips            = 0
    ! Banner
    character(len=1) :: BannerAlignment          = 'C'
    character(len=1) :: BannerPattern            = '-'
    integer, dimension(2) :: BannerColumnrange   = (/ 1, 80 /)
    integer          :: BannerSkips              = 0
    integer          :: BannerLength             = 0
  end type
  type(StyleOptions_T), private, save  :: DefaultStyleOptions
  type(StyleOptions_T), public, save   :: StyleOptions

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------  addRow  -----
  ! This family of routines adds a paired name and value
  ! as a new row to the cellDatabase
  
  ! later to be printed as a neatly-formatted table to stdout
  ! By means of optional args you can create a line like
  ! *   name                   value   *
  subroutine addRow_character ( name, value, format )
    character(len=*), intent(in)          :: name
    character(len=*), intent(in)          :: value
    include 'addRow.f9h'
  end subroutine addRow_character

  subroutine addRow_complex ( name, value, format )
    character(len=*), intent(in)          :: name
    complex, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine addRow_complex

  subroutine addRow_double ( name, value, format )
    character(len=*), intent(in)          :: name
    double precision, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine addRow_double

  subroutine addRow_dbl_array ( name, value, format )
    character(len=*), intent(in)          :: name
    double precision, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine addRow_dbl_array

  subroutine addRow_int_array ( name, value, format )
    character(len=*), intent(in)          :: name
    integer, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine addRow_int_array

  subroutine addRow_integer ( name, value, format )
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine addRow_integer

  subroutine addRow_log_array ( name, value, format )
    character(len=*), intent(in)          :: name
    logical, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine addRow_log_array

  subroutine addRow_logical ( name, value, format )
    character(len=*), intent(in)          :: name
    logical, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine addRow_logical

  subroutine addRow_single ( name, value, format )
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    include 'addRow.f9h'
  end subroutine addRow_single

  subroutine addRow_sngl_array ( name, value, format )
    character(len=*), intent(in)          :: name
    real, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine addRow_sngl_array

  subroutine addRow_header ( name, alignment )
    character(len=*), intent(in)          :: name
    character(len=1), intent(in)          :: alignment !: 'l(eft)', 'c', or 'r'
    character(len=MAXCELLSIZE), dimension(2)   :: item = ' '
    integer                                    :: newSize
    item(1) = '<<' // alignment // '>>' // name
    newSize = addCellRowToDatabase( cellDatabase, item )
  end subroutine addRow_header

  subroutine addRow_divider ( char )
    character(len=1), intent(in)          :: char
    character(len=MAXCELLSIZE), dimension(2)   :: item = ' '
    integer                                    :: newSize
    item(1) = '<<' // 'd' // '>>' // char
    newSize = addCellRowToDatabase( cellDatabase, item )
  end subroutine addRow_divider

  ! -----------------------------------------------------  ALIGNTOFIT  -----
  ! Align chars to fit within column range
  ! Alignment controls whether the chars are
  ! L    Flushed left
  ! R    Flushed right
  ! C    Centered
  ! J    Justified (padding spaces to any existing spaces)
  subroutine ALIGNTOFIT_CHARS ( CHARS, COLUMNRANGE, ALIGNMENT, SKIPS )
    character(len=*), intent(in)      :: CHARS
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), intent(in) :: COLUMNRANGE
    character(len=1), intent(in)      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)     :: SKIPS ! How many spaces between chars
    !
    ! Internal variables
    character(len=max(len(chars), abs(columnRange(2)-columnRange(1)))) :: &
      & ALLCHARS
    integer :: char1
    integer :: char2
    integer :: firstSpace
    integer :: m
    integer :: nc
    integer :: padLeft
    integer :: padRight
    integer :: spaces
    ! Executable
    allChars = chars
    if ( present(skips) ) then
      if ( skips > 0 .and. len_trim(chars) > 0 ) then
        allChars = stretch(chars, skips)
      end if
    end if
    if ( columnRange(1) > 0 ) then
      spaces = columnRange(2) - max( columnRange(1), getOutputStatus( 'column' ) )
      if ( spaces < 1 ) return
      if ( columnRange(1) > getOutputStatus( 'column' ) ) &
        & call blanks( columnRange(1) - getOutputStatus( 'column' ) )
    else
      spaces = columnRange(2) - columnRange(1)
    end if
    firstSpace = 0
    nc = max( len_trim(allchars), 1 )
    select case (lowercase(alignment))
    case ('l')
      char1    = 1
      padLeft  = 0
      char2    = min( nc, spaces )
      padRight = spaces - char2
    case ('r')
      char1    = max(1, nc-spaces+1)
      char2    = nc
      padLeft  = spaces - (char2-char1+1)
      padRight = 0
    case ('c', 'j')
      m = (spaces - nc) / 2
      padLeft  = max( m, 0 )
      padRight = max( spaces - nc - m, 0 )
      char1 = max(1-m, 1)
      char2 = min(nc+m, nc)
      if ( lowercase(alignment) == 'j' .and. padRight > 0 ) &
        & firstSpace = index( allChars, ' ' )
    end select
    if ( firstSpace > 1 ) then
      call output_( allChars(char1:firstSpace-1) )
      call blanks( padRight+padLeft+1 )
      if ( firstSpace+1 < char2 ) call output_( allChars(firstSpace+1:char2) )
    else
      call blanks( padLeft )
      call output_( allChars(char1:char2) )
      call blanks( padRight )
    end if
  end subroutine ALIGNTOFIT_CHARS

  subroutine ALIGNTOFIT_DOUBLE ( value, COLUMNRANGE, ALIGNMENT, FORMAT )
    double precision, intent(in)      :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), intent(in) :: COLUMNRANGE
    character(len=1), intent(in)      :: ALIGNMENT ! L, R, C, or J
    character(len=*), optional, intent(in)     :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call alignToFit( trim(line), columnRange, alignment )
  end subroutine ALIGNTOFIT_DOUBLE

  subroutine ALIGNTOFIT_INTEGER ( value, COLUMNRANGE, ALIGNMENT, FORMAT )
    integer, intent(in)               :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), intent(in) :: COLUMNRANGE
    character(len=1), intent(in)      :: ALIGNMENT ! L, R, C, or J
    character(len=*), optional, intent(in)     :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call alignToFit( trim(line), columnRange, alignment )
  end subroutine ALIGNTOFIT_INTEGER

  subroutine ALIGNTOFIT_SINGLE ( value, COLUMNRANGE, ALIGNMENT, FORMAT )
    real, intent(in)      :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), intent(in) :: COLUMNRANGE
    character(len=1), intent(in)      :: ALIGNMENT ! L, R, C, or J
    character(len=*), optional, intent(in)     :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call alignToFit( trim(line), columnRange, alignment )
  end subroutine ALIGNTOFIT_SINGLE

  ! -----------------------------------------------------  BANNER  -----
  ! Surround your message with stars and stripes; e.g.,
  ! *-----------------------------------------------*
  ! *            Your message here                  *
  ! *-----------------------------------------------*
  ! proclaiming its great importance to an uncaring world.
  ! For multiline messages, you may divide them into elements of
  ! a character array, or else a longer character scalar and
  ! supply LineLength asking the routine to wrap at word boundaries
  subroutine BANNER_CHARS ( chars, &
    & columnRange, alignment, skips, lineLength, mode, pattern )
    character(len=*), intent(in)                :: CHARS
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    integer, optional, intent(in)               :: LINELENGTH
    character (len=*), optional, intent(in)     :: mode ! if not 'hard'
    character (len=1), optional, intent(in)     :: pattern ! if not stripes
    !
    ! Internal variables
    integer :: addedLines
    character(len=1)      :: myAlignment
    character(len=1)      :: myFillChar
    integer, dimension(2) :: myColumnRange
    integer :: lineLen, mySkips, padding
    character(len=2*len(chars))      :: wrappedChars
    character(len=160), dimension(:), allocatable :: lines
    logical, parameter :: DEBUG = .false.
    ! Executable
    myAlignment = StyleOptions%BannerAlignment ! 'C'
    if ( present(alignment) ) myAlignment = alignment
    mySkips = StyleOptions%BannerSkips ! 0
    if ( present(skips) ) mySkips = skips
    myFillChar = StyleOptions%BannerPattern ! '-'
    if ( present(pattern) ) myFillChar = pattern
    lineLen = StyleOptions%BannerLength ! 0
    if ( present(LineLength) ) lineLen = LineLength
    if ( lineLen > 4 ) then
      ! We will wrap the input to fit within LineLength, but remembering
      ! the stars and padding
      call wrap( chars, wrappedChars, width=lineLen-4, &
        & inseparator=achar(0), addedLines=addedLines, mode=mode )
      addedLines = addedLines + 1
      allocate( lines(addedLines) )
      lines = ' '
      call List2Array( wrappedChars, lines, &
        & countEmpty=.true., inseparator=achar(0) )
      call banner( lines, alignment=alignment, pattern=pattern )
      deallocate(lines)
      return
    else if ( present(columnRange) ) then
      myColumnRange = columnRange
    else
      lineLen = max( 80, 4 + len_trim(chars)*(1+mySkips) )
      padding = ( lineLen - len_trim(chars)*(1+mySkips) ) / 2
      myColumnRange(1) = 1 + padding
      myColumnRange(2) = lineLen - padding
    end if
    
    ! define padding as the larger of columnrange(1) and 1
    padding = max( 1, myColumnRange(1) )
    LineLen = padding + myColumnRange(2) - 1
    if ( DEBUG ) then
      call outputnamedValue( 'padding', padding )
      call outputnamedValue( 'LineLen', LineLen )
      call outputnamedValue( 'myColumnRange', myColumnRange )
    end if
    ! Top border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
    ! Left star, then message, then right star
    call output( '*' )
    call alignToFit( chars, myColumnRange, myAlignment, skips )
    call blanksToColumn( lineLen )
    call output( '*', advance = 'yes' )
    ! Bottom border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
  end subroutine BANNER_CHARS

  subroutine BANNER_CHARARRAY ( charArray, &
    & columnRange, alignment, skips, pattern )
    character(len=*), dimension(:), intent(in)  :: CHARARRAY
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    character (len=1), optional, intent(in)     :: pattern ! if not stripes
    !
    ! Internal variables
    integer :: i
    ! Internal variables
    character(len=1)      :: myAlignment
    integer, dimension(2) :: myColumnRange
    character(len=1)      :: myFillChar
    integer :: lineLen, mySkips,  padding
    ! Executable
    myAlignment = 'C'
    if ( present(alignment) ) myAlignment = alignment
    mySkips = 0
    if ( present(skips) ) mySkips = skips
    myFillChar = '-'
    if ( present(pattern) ) myFillChar = pattern
    if ( present(columnRange) ) then
      myColumnRange = columnRange
    else
      lineLen = 80
      padding = LineLen
      do i = 1, size(chararray)
        lineLen = max( lineLen, 4 + len_trim(chararray(i))*(1+mySkips) )
        padding = min( padding, &
          & ( lineLen - len_trim(chararray(i))*(1+mySkips) ) / 2 )
      end do
      myColumnRange(1) = 1 + padding
      myColumnRange(2) = lineLen - padding
    end if
    
    ! define padding as the larger of columnrange(1) and 1
    padding = max( 1, myColumnRange(1) )
    LineLen = padding + myColumnRange(2) - 1
    ! Top border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
    do i = 1, size(chararray)
      ! Left star, then message, then right star
      call output( '*' )
      call alignToFit( chararray(i), myColumnRange, myAlignment, skips )
      call blanksToColumn( lineLen )
      call output( '*', advance = 'yes' )
    end do
    ! Bottom border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
  end subroutine BANNER_CHARARRAY

  ! -----------------------------------------------------  BEVERBOSE  -----
  logical function BEVERBOSE ( SWITCH, THRESHOLD )
    ! Args
    character(len=*), intent(in) :: SWITCH
    integer, intent(in)          :: THRESHOLD
    ! Executable
    BeVerbose = switchDetail( switches, switch ) > threshold .or. MLSVerbose
  end function BEVERBOSE

  ! -----------------------------------------------------  BLANKSTOCOLUMN  -----
  subroutine BLANKSTOCOLUMN ( column, fillchar, advance, dont_stamp )
  ! Output N_BLANKS blanks to PRUNIT out to column COLUMN.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: COLUMN
    character(len=*), intent(in), optional :: advance
    character(len=*), intent(in), optional :: fillchar  ! default is ' '
    logical, intent(in), optional          :: dont_stamp ! Prevent double-stamping
    ! Internal variables
    integer :: nblanks
    ! Executable
    if ( getOutputStatus( 'physicalcolumn' ) == 1 ) then
      ! Don't double indent the line's first field
      nblanks = column - getOutputStatus( 'physicalcolumn' )
    else
      nblanks = column - getOutputStatus( 'column' )
    endif
    if ( nblanks < 1 ) return
    call blanks( nblanks, fillChar, advance, dont_stamp )
  end subroutine BLANKSTOCOLUMN

  ! ------------------------------------------------  blanksToTab  -----
  ! Print blanks out to next tabstop
  ! (or else to tabstop number tabn)
  subroutine blanksToTab ( tabn, fillChar )
    ! Args
    integer, optional, intent(in) :: TABN
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    ! Internal variables
    integer :: nColumn
    integer :: nTab
    ! Executable
    if ( getOutputStatus( 'physicalcolumn' ) == 1 ) then
      ! Don't double indent the line's first field
      nColumn = getOutputStatus( 'physicalcolumn' )
    else
      nColumn = getOutputStatus( 'column' )
    endif
    if ( present(tabn) ) then
      if ( tabn < 1 .or. tabn > MAXNUMTABSTOPS ) return
      if ( nColumn < tabStops(tabn) ) &
        & call blanksToColumn( tabStops(tabn), fillChar )
    else
      nTab = findFirst( tabStops > nColumn )
      if ( nTab > 0 ) &
        & call blanksToColumn( tabStops(nTab), fillChar )
    end if
  end subroutine blanksToTab

  ! ---------------------------------------------- DumpOuputOptions -----
  subroutine DumpOutputOptions( options )
    ! Show output options
    type(outputOptions_T), intent(in) :: options
    ! Internal variables
    logical, parameter :: checkingTabbing = .false.
    character(len=10), parameter :: decade = '1234567890'
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
    integer :: i
    ! Executable
    call blanks(80, fillChar='-', advance='yes')
    call headline( 'Summary of output options', &
      & fillChar='-', before='*', after='*' )
    call outputNamedValue ( 'unit number', options%prUnit, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    if ( options%prUnit < 0 ) then
      call outputNamedValue ( 'meaning', prunitname(options ), &
        & advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    end if
    call outputNamedValue ( 'file name', trim(options%name), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'logging level', options%MLSMSG_Level, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'buffered?', options%buffered, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'skip MLSMSG logging?', options%SKIPMLSMSGLOGGING, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'log Parent Name?', options%logParent, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'use patterned blanks?', options%usePatternedBlanks, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'special fills', trim(options%specialFillChars), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'lineup fills', trim(options%lineupFillChars), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'tab stops', tabstops, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    do i=1, MAXNUMTABSTOPS
      call tab( fillChar=fillChar )
      call output_( '^', advance='no' )
    end do
    call newline
    if ( CHECKINGTABBING ) then
      do i=1, MAXNUMTABSTOPS
        call blanksToColumn( tabStops(i), fillChar=fillChar )
        call output_( '^', advance='no' )
      end do
      call newline
    end if
    do i = 1, 14
      call output_( decade, advance='no' )
      if ( getOutputStatus( 'column' ) > 132 ) exit
    end do
    call newline
    call blanks(80, fillChar='-', advance='yes')
  contains
    function PRUnitName( options ) result( name )
      ! Return an appropriate name for the prUnit number
      ! Args
      type(outputOptions_T), intent(in) :: options
      character(len=12) :: name
      ! Executable
      if ( options%prUnitLiteral ) then
        name = 'literal'
      else
        select case ( options%prUnit )
        case ( STDOUTPRUNIT )
          name = 'stdout'
        case ( MSGLOGPRUNIT )
          name = 'mls logfile'
        case ( BOTHPRUNIT )
          name = 'stdout+log'
        case ( OUTPUTLINESPRUNIT )
          name = 'outputLines'
        case ( INVALIDPRUNIT )
          name = 'invalid'
        case default ! > 0
          name = 'Fortran unit'
        end select
      end if
    end function PRUnitName
  end subroutine DumpOutputOptions

  ! ---------------------------------------------- DumpStampOptions -----
  subroutine DumpStampOptions( options )
    ! Show output options
    type(StampOptions_T), intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     call blanks(80, fillChar='-', advance='yes')
    call headline( 'Summary of automatic stamp options', &
      & fillChar='-', before='*', after='*' )
     call outputNamedValue ( 'never stamp', options%neverStamp, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'stamp end of line', options%post, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'show time', options%showTime, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'extra text', trim_safe(options%textCode), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'date format', trim_safe(options%dateFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'time format', trim_safe(options%timeFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'interval', options%interval, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'style of timeStamps', trim_safe(options%timestampstyle), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call blanks(80, fillChar='-', advance='yes')
  end subroutine DumpStampOptions

  ! ---------------------------------------------- DUMPTIMESTAMPOPTIONS -----
  subroutine DUMPTIMESTAMPOPTIONS( options )
    ! Show output options
    type(TimeStampOptions_T), intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     call blanks(80, fillChar='-', advance='yes')
    call headline( 'Summary of time stamp options', &
      & fillChar='-', before='*', after='*' )
     call outputNamedValue ( 'stamp end of line', options%post, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'show date', options%showDate, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'extra text', trim_safe(options%textCode), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'date format', trim_safe(options%dateFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'time format', trim_safe(options%timeFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'style of timeStamps', trim_safe(options%timestampstyle), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call blanks(80, fillChar='-', advance='yes')
  end subroutine DUMPTIMESTAMPOPTIONS

  ! ---------------------------------------------- DumpSize_double -----
  subroutine DumpSize_double ( n, advance, units, Before, After, Signed )
    double precision, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    real, intent(in), optional :: units
    character(len=*), intent(in), optional :: Before, After
    logical, intent(in), optional :: Signed ! Force "+" if positive
    ! Local parameters
    real, parameter :: KB = 1024.0
    real, parameter :: MB = KB * 1024.0
    real, parameter :: GB = MB * 1024.0
    real, parameter :: TB = GB * 1024.0
    double precision :: Amount ! N * MyUnits
    character(len=8) :: HowMuch
    character        :: mySigned ! blank or "+"
    real             :: myUnits
    character(len=6) :: Suffix
    ! Make a 'nice' output
    if ( present(before) ) call output_ ( before )
    mySigned = ''
    if ( present(signed) ) then
      if ( signed ) mySigned = merge('+', ' ', n >= 0)
    end if
    myUnits = 1.0
    if ( present(units) ) myUnits = units
    if ( myUnits == 0.0 ) then
      call output ( n, format='(e12.1)', before=trim(mySigned) )
      call output_ ( ' (illegal units)', advance=advance )
      return
    end if
    amount = n*myUnits
    if ( abs(n) < kb/myUnits ) then
      suffix = ' bytes'
    else if ( abs(n) < Mb/myUnits ) then
      amount = amount/kb
      suffix = ' kB'
    else if ( abs(n) < Gb/myUnits ) then
      amount = amount/Mb
      suffix = ' MB'
    else if ( abs(n) < Tb/myUnits ) then
      amount = amount/Gb
      suffix = ' GB'
    else
      amount = amount/Tb
      suffix = ' TB'
    end if
    if ( amount < -99999 ) then     ! I6 format limits this
      call output_( '(-HUGE)' )
    else if ( amount > 999999 ) then ! I6 format limits this
      call output_( '(HUGE)' )
    else if ( amount == int(amount) ) then
      write ( howMuch, '(i6)' ) int(amount)
    else
      write ( howMuch, '(f8.1)' ) amount
    end if
    call output_ ( trim(mySigned) )
    call output_ ( trim(adjustl(howMuch)) )
    call output_ ( trim(suffix) )
    if ( present(after) ) call output_ ( after )
    call output_ ( '', advance=advance )
  end subroutine DumpSize_double

  ! --------------------------------------------- DumpSize_integer -----
  subroutine DumpSize_integer ( n, advance, units, Before, After, Signed )
    integer, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    integer, intent(in), optional :: units ! E.g., 1024 for kB
    character(len=*), intent(in), optional :: Before, After
    logical, intent(in), optional :: Signed ! Force "+" if positive
    ! Executable
    if ( present(units) ) then
      call dumpSize ( dble(n), advance=advance, units=real(units), &
        & before=before, after=after, signed=signed )
    else
      call dumpSize ( dble(n), advance=advance, before=before, after=after, &
        & signed=signed )
    end if
  end subroutine DumpSize_integer

  ! ------------------------------------------------ DumpSize_real -----
  subroutine DumpSize_real ( n, advance, units, Before, After, Signed )
    real, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    real, intent(in), optional :: units
    character(len=*), intent(in), optional :: Before, After
    logical, intent(in), optional :: Signed ! Force "+" if positive
    ! Make a 'nice' output
    call dumpsize ( dble(n), advance, units, before, after, signed )
  end subroutine DumpSize_real

  ! ----------------------------------------------  dumpTabs  -----
  ! Show tab stops in effect
  ! Optionally returning them as an integer array
  subroutine dumpTabs ( tabs )
    ! Args
    integer, dimension(:), optional, intent(out) :: tabs
    ! Internal variables
    integer :: n
    ! Executable
    call output( 'Current tab stops', advance='yes' )
    call output( TABSTOPS, advance='yes' )
    if ( present(tabs) ) then
      n = min(size(tabs), MAXNUMTABSTOPS)
      tabs = 0
      tabs(1:n) = TABSTOPS(1:n)
    end if
  end subroutine dumpTabs

  ! ----------------------------------------------  getStamp  -----
  subroutine getStamp ( textCode, showTime, dateFormat, timeFormat, &
    & post, interval )
  ! get stamp being added to every output to PRUNIT.
    character(len=*), optional, intent(out) :: textCode
    logical, optional, intent(out)          :: showTime
    character(len=*), optional, intent(out) :: dateFormat
    character(len=*), optional, intent(out) :: timeFormat
    logical, optional, intent(out)          :: post
    integer, optional, intent(out)          :: interval
    if ( present(showTime) )   showTime   = stampOptions%showTime
    if ( present(textCode) )   textCode   = stampOptions%textCode
    if ( present(dateFormat) ) dateFormat = stampOptions%dateFormat
    if ( present(timeFormat) ) timeFormat = stampOptions%timeFormat
    if ( present(post) )       post       = stampOptions%post
    if ( present(interval) )   interval   = stampOptions%interval
  end subroutine getStamp

  ! -----------------------------------------------------  HEADLINE  -----
  ! Print your message with extra formatting features; e.g.,
  ! *----------------  Your message here   ----------------*
  ! See also banner
  subroutine HEADLINE ( CHARS, fillChar, Before, After, &
    & ColumnRange, Alignment, Skips )
    character(len=*), intent(in)                :: CHARS
    character(len=1), intent(in), optional      :: fillChar      ! For padding
    character(len=*), intent(in), optional      :: Before, After ! text to print
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    !
    ! Internal variables
    character(len=1)      :: myAlignment
    integer, dimension(2) :: myColumnRange      ! To fit chars
    integer, dimension(2) :: myFullColumnRange  ! Must fit before and after, too
    integer :: mySkips,  rightpadding
    character(len=1) :: myFillChar
    character(len=8) :: myBefore, myAfter
    ! Executable
    rightpadding = 0
    if ( present(columnRange) ) then
      myFullColumnRange = columnRange
    else
      myFullColumnRange = StyleOptions%HeadlineColumnRange ! (/ 1, 80 /)
    end if
    ! print *, 'myFullColumnRange ', myFullColumnRange 
    myColumnRange = myFullColumnRange
    mySkips = StyleOptions%HeadlineSkips ! 0
    if ( present(skips) ) mySkips = skips
    myFillChar = StyleOptions%Headlinefill ! ' '
    if ( present(fillChar) ) myFillChar = fillChar
    myAlignment = StyleOptions%HeadlineAlignment ! 'C'
    if ( present(alignment) ) myAlignment = alignment
    myBefore = StyleOptions%HeadlineBefore
    myAfter = StyleOptions%HeadlineAfter
    ! call outputNamedValue ( 'myFullColumnRange', myFullColumnRange )
    ! call outputNamedValue ( 'myColumnRange', myColumnRange )
    ! call outputNamedValue ( 'mySkips', mySkips )
    ! call outputNamedValue ( 'myFillChar', myFillChar )
    ! Adjust for lengths of before, after
    if ( len_trim(myBefore) > 0 ) &
      & myColumnRange(1) = myFullColumnRange(1) + len_trim(myBefore)
    if ( len_trim(myAfter) > 0 ) then
      myColumnRange(2) = myFullColumnRange(2) - len_trim(myAfter)
      rightpadding = len_trim(myAfter) - 1
    end if
    call blanksToColumn( myFullColumnRange(1), advance='no' )
    if ( len_trim(myBefore) > 0 ) call output( trim(myBefore), advance='no' )
    if ( mySkips == 0 .and. myAlignment == 'C' .and. myFillChar /= ' ' ) then
      ! OK, final adjustments of myColumnRange
      myColumnRange(1) = &
        & ( myColumnRange(1) + myColumnRange(2) + 1 - len(chars) )/2
      myColumnRange(2) =  myColumnRange(1) - 1 + len(chars)
      ! print *, 'myColumnRange ', myColumnRange
      call blanksToColumn( myColumnRange(1), fillChar=myFillChar, advance='no' )
      call aligntofit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn( myFullColumnRange(2)-rightpadding, &
        & fillChar=myFillChar, advance='no' )
      if ( len_trim(myAfter) > 0 ) call output( trim(myAfter), advance='no' )
    else
      call aligntofit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn( myFullColumnRange(2)-rightpadding, advance='no' )
      if ( len_trim(myAfter) > 0 ) call output( trim(myAfter), advance='no' )
    end if
    call newLine
  end subroutine HEADLINE

  ! -----------------------------------------------------  LETSDEBUG  -----
  logical function LETSDEBUG ( SWITCH, THRESHOLD )
    ! Args
    character(len=*), intent(in) :: SWITCH
    integer, intent(in)          :: THRESHOLD
    ! Executable
    letsDebug = switchDetail( switches, switch ) > threshold .or. MLSDebug
  end function LETSDEBUG

  ! ----------------------------------------------------  NextColumn  -----
  function NextColumn() result(Column)
    ! Args
    integer :: Column
    Column = getOutputStatus( 'column' )
  end function NextColumn

  ! ----------------------------------------------------  NextTab  -----
  function NextTab() result(Column)
    ! Args
    integer :: Column
    ! Internal variables
    integer :: nTab
    ! Executable
    Column = 0
    nTab = findFirst( tabStops > getOutputStatus( 'column' ) )
    if ( nTab > 0 ) Column = max( tabStops(nTab), getOutputStatus( 'column' ) )
  end function NextTab

  ! ----------------------------------------------------  numNeedsFormat  -----
  ! This family of functions return what format is needed to be printed by output
  function numNeedsFormat_double( value, inFormat ) result ( format )
    ! Args
    double precision, intent(in) :: VALUE
    character(len=*), optional, intent(in)  :: inFormat
    character(len=30) :: format
    ! Internal variables
    character(len=30) :: charValue
    character(len=2)  :: dotm
    character(len=30) :: ndotm
    integer :: I
    ! Executable
    call whatSDNeedsFormat( ndotm, dotm, inFormat )
    charValue = adjustl(numToChars ( value, format=ndotm ))
    ! call outputNamedValue( 'ndotm', ndotm )
    ! call outputNamedValue( 'charValue', charValue )
    i = len_trim(charValue)
    write(charValue, *) i+5
    format = '(1pg' // trim(adjustl(charValue)) // dotm // ')'
  end function numNeedsFormat_double

  function numNeedsFormat_integer( value, inFormat ) result ( format )
    ! Args
    integer, intent(in) :: VALUE
    character(len=*), optional, intent(in)  :: inFormat
    character(len=30) :: format
    ! Internal variables
    character(len=30) :: charValue
    integer :: I
    ! Executable
    charValue = numToChars(value)
    i = len_trim(charValue)
    write(charValue, *) i+5
    format = '(i' // trim(adjustl(charValue)) // ')'
  end function numNeedsFormat_integer

  function numNeedsFormat_single( value, inFormat ) result ( format )
    ! Args
    real, intent(in) :: VALUE
    character(len=*), optional, intent(in)  :: inFormat
    character(len=30) :: format
    ! Internal variables
    character(len=30) :: charValue
    character(len=2)  :: dotm
    character(len=30) :: ndotm
    integer :: I
    ! Executable
    call whatSDNeedsFormat( ndotm, dotm, inFormat )
    charValue = adjustl(numToChars ( value, format=ndotm ))
    i = len_trim(charValue)
    write(charValue, *) i+5
    format = '(1pg' // trim(adjustl(charValue)) // dotm // ')'
  end function numNeedsFormat_single

  function numNeedsFormat_complex( value, inFormat ) result ( format )
    ! Args
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: VALUE
    character(len=*), optional, intent(in)  :: inFormat
    character(len=45) :: format
    ! Internal variables
    character(len=30) :: charValue
    character(len=2)  :: dotm
    character(len=30) :: ndotm
    integer :: I
    ! Executable
    call whatSDNeedsFormat( ndotm, dotm, inFormat )
    charValue = adjustl(numToChars ( abs(value), format=ndotm ))
    i = len_trim(charValue)
    write(charValue, *) i+5
    format = '(1x,"(",1pg' // trim(adjustl(charValue)) // dotm // ',",",1pg' &
      & // trim(adjustl(charValue)) // dotm // ',")")'
  end function numNeedsFormat_complex

  function numNeedsFormat_dcomplx( value, inFormat ) result ( format )
    ! Args
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: VALUE
    character(len=*), optional, intent(in)  :: inFormat
    character(len=45) :: format
    ! Internal variables
    character(len=30) :: charValue
    character(len=2)  :: dotm
    character(len=30) :: ndotm
    integer :: I
    ! Executable
    call whatSDNeedsFormat( ndotm, dotm, inFormat )
    charValue = adjustl(numToChars ( abs(value), format=ndotm ))
    i = len_trim(charValue)
    write(charValue, *) i+5
    format = '(1x,"(",1pg' // trim(adjustl(charValue)) // dotm // ',",",1pg' &
      & // trim(adjustl(charValue)) // dotm // ',")")'
  end function numNeedsFormat_dcomplx

  ! ----------------------------------------------------  numToChars  -----
  ! This family of functions return what would otherwise be printed by output
  function numToChars_double( value, format ) result ( line )
    ! Args
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=30) :: line
    ! Internal variables
    character(len=30) :: FormatSpec
    integer :: I, J, K
    ! Executable
    FormatSpec = outputOptions%sdFormatDefault
    if ( any( value == DPREFERDEFAULTFORMAT ) ) FormatSpec = '*'
    if ( present(Format)  ) then
      if ( format /= '*' ) FormatSpec = Format
    end if
    include 'numToChars.f9h'
  end function numToChars_double

  function numToChars_integer( value, format ) result ( line )
    ! Args
    integer, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=30) :: line
    ! Executable
    if ( present(Format)  ) then
      write( line, Format ) value
    else
      write( line, * ) value
    end if
    line = adjustl(line)
  end function numToChars_integer

  function numToChars_single( value, format ) result ( line )
    ! Args
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=30) :: line
    ! Internal variables
    character(len=30) :: FormatSpec
    integer :: I, J, K
    ! Executable
    FormatSpec = outputOptions%sdFormatDefault
    if ( any( value == DPREFERDEFAULTFORMAT ) ) FormatSpec = '*'
    if ( present(Format)  ) then
      if ( format /= '*' ) FormatSpec = Format
    end if
    include 'numToChars.f9h'
  end function numToChars_single

  ! ---------------------------------------  OUTPUTCALENDAR  -----
  subroutine OUTPUTCALENDAR ( date, datenote, notes, dontWrap, moonPhases )
    use Dates_Module, only: NextMoon
    use MLSStringLists, only: CatLists
    ! output a nicely-formatted calendar of the current month with
    ! today's date marked in "bold"
    ! Args
    character(len=*), intent(in), optional :: date ! use date instead of today
    ! dateNote, (notes), if present, is (an array of)
    ! stringLists, (one per day in the month,)
    character(len=*), optional :: dateNote ! Note for the current date
    ! Each string list contains either a blank for a date, meaning
    ! nothing will be printed in the calendar square for that date,
    ! or else it contains '/'-separated lines of text, each of
    ! which will be printed on a separate line within the square
    character(len=*), dimension(:), optional :: notes
    logical, optional                        :: dontWrap ! Dont wrap notes to fit
    logical, optional                        :: moonPhases ! Show new, full moons
    ! Internal variables
    integer, parameter :: MAXNOTELENGTH = 256
    ! This should be modified for internationalization; e.g. with
    ! an include statement or suchlike
    ! The next two arrays may be customized for non-English users
    character(len=*), dimension(12), parameter :: MONTHNAME = (/ &
      & 'January  ', 'February ', 'March    ', 'April    ', 'May      ', &
      & 'June     ', 'July     ', 'August   ', 'September', 'October  ', &
      & 'November ', 'December '/)

    character(len=*), dimension(7), parameter :: DAYSOFWEEK = (/ &
      & 'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', 'Thursday ', &
      & 'Friday   ', 'Saturday '/)
    logical, parameter :: countEmpty = .true.
    character(len=1), parameter :: inseparator = '/'
    character(len=*), parameter :: utcformat = 'yyyy-mm-dd' ! 'yyyy-Doy'
    integer :: aday
    integer :: col1
    integer :: col2
    character(len=16) :: date2, dateString
    character(len=16) :: fullDate, newDate
    integer :: day
    integer, dimension(6,7) :: days, daysOfYear
    integer :: ErrTyp
    integer :: FullMoonDay
    integer :: iwk
    integer :: month
    logical :: myDontWrap
    integer :: NewMoonDay
    character(len=10) :: noteString
    character(len=128), dimension(42) :: myNotes
    integer :: numRows
    integer :: numWeeks
    integer :: row
    logical :: today
    integer :: wkdy
    character(len=MAXNOTELENGTH) :: wrappedNote
    integer :: year
    ! Executable
    myDontWrap = .false.
    if ( present(dontWrap) ) myDontWrap = dontWrap
    if ( present(date) ) then
      dateString = date
    else
      dateString = '' ! Intel 12 and earlier doesn't fill with blanks
      call date_and_time ( date=dateString )
    end if
    myNotes = ' '
    if ( present(notes) ) myNotes(1:size(notes) ) = notes
    col1 = index(lowercase(dateString), 't')
    if ( col1 > 0 ) then
      date2 = dateString(1:col1-1)
       if ( nCopies(dateString(:col1-1), '-') < 2 ) &
        & date2 = reformatDate( dateString(1:col1-1), fromForm='yyyy-Doy', toForm=utcformat )
    else
      date2 = reformatDate( dateString, fromForm='*', toForm=utcformat )
    end if
    call utc_to_yyyymmdd( date2, ErrTyp, year, month, day )
    if ( month < 0 ) then
    end if
    FullMoonDay = 0
    NewMoonDay = 0
    if ( present(moonPhases) ) then
      if ( moonPhases) then
        fullDate = nextMoon( date2(1:8) // '01', 'full' )
        newDate = nextMoon( date2(1:8) // '01', 'new' )
        read ( fullDate(9:10), * ) FullMoonDay
        read ( newDate(9:10), * ) NewMoonDay
        ! Add these to notes if present
        if ( present(notes) ) then
          do aday=1, min( size(notes), daysInMonth( month, year ) )
            if ( aday == FullMoonday ) then
              myNotes(aday) = catLists( 'Full moon', myNotes(aday), inseparator )
              ! call outputnamedValue ( 'fullmoonday', trim(myNotes(aday)) )
            elseif ( aday == newMoonday ) then
              myNotes(aday) = catLists( 'New moon', myNotes(aday), inseparator )
              ! call outputnamedValue ( 'newmoonday', trim(myNotes(aday)) )
            endif
          enddo
        endif
      endif
    endif
    call buildCalendar( year, month, days, daysOfYear )
    ! Temporary use of   w i d e  tabstops
    call settabs( '14-210+14' )
    call newline
    call alignToFit( trim(monthName(month)), (/ 1, 100 /), 'c', skips=1 )
    call newline
    col2 = 0
    do wkdy=1, 7
      col1 = col2 + 1
      col2 = tabStops(wkdy)
      call alignToFit( trim(daysOfWeek(wkdy)), (/ col1, col2 /), 'c' )
    end do
    call newline
    numWeeks = 4
    if ( any( days(5,:) /= 0 ) ) numWeeks = 5
    if ( any( days(6,:) /= 0 ) ) numWeeks = 6
    ! How many rows will we need?
    numRows = 4
    if ( present(dateNote) ) then
      if ( myDontWrap ) then
        wrappedNote = dateNote
      else
        call wrap( dateNote, wrappedNote, 10, '/' )
      end if
      numRows = max( numRows, &
        & NumStringElements( wrappedNote, countEmpty, inseparator ) + 2 &
        & )
    end if
    if ( present(notes) ) then
      do aday=1, min( size(notes), daysInMonth( month, year ) )
        if ( myDontWrap ) then
          wrappedNote = myNotes(aday)
        else
          call wrap( myNotes(aday), wrappedNote, 10, '/' )
        end if
        numRows = max( numRows, &
          & NumStringElements( wrappedNote, countEmpty, inseparator ) + 2 &
          & )
      end do
    end if
    do iwk = 1, numWeeks
      ! Start with horizontal line
      call blanksToTab( 7, fillChar='-' )
      call newline
      do row=1, numRows
        col2 = 0
        do wkdy=1, 7
          col1 = col2 + 1
          col2 = tabStops(wkdy)
          today = ( days(iwk, wkdy) == day )
          if ( today ) then
            call output_('||')
          else
            call output_('|')
          end if
          if ( days(iwk, wkdy) < 1 ) then
            ! Don't write notes or anything else in "empty" days
          else if ( row == 1 ) then
            call writeIntsToChars( days(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call alignToFit( trim(dateString), (/ col1, col2-1 /), 'r' )
          else if( row == numRows ) then
            call writeIntsToChars( daysOfYear(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call alignToFit( 'd' // trim(dateString), (/ col1, col2-1 /), 'r' )
          else if( present(dateNote) .and. today ) then
            if ( myDontWrap ) then
              wrappedNote = dateNote
            else
              call wrap( dateNote, wrappedNote, 10, '/' )
            end if
            call GetStringElement ( wrappedNote, noteString, &
              & row-1, countEmpty, inseparator )
            if ( noteString == inseparator ) noteString = ' '
            call output_( noteString )
          else if( present(notes) ) then
            if ( days(iwk, wkdy) <= size(notes) ) then
              if ( myDontWrap ) then
                wrappedNote = myNotes(days(iwk, wkdy))
              else
                call wrap( myNotes(days(iwk, wkdy)), wrappedNote, 10, '/' )
              end if
              call GetStringElement ( wrappedNote, noteString, &
                & row-1, countEmpty, inseparator )
              if ( noteString == inseparator ) noteString = ' '
              call output_( noteString )
            end if
          else if ( days(iwk, wkdy) == FullMoonday .and. row < 3 ) then
            call output_( 'Full moon' )
          else if ( days(iwk, wkdy) == NewMoonday .and. row < 3 ) then
            call output_( 'New moon' )
          end if
          if ( today ) then
            call blanksToColumn(col2-1)
            call output_('|')
          else
            call blanksToTab
          end if
        end do ! wkdy
        call output_('|')
        call newline
      end do ! row
      ! begin with 
    end do ! week
    call blanksToTab( 7, fillChar='-' )
    call newline
    ! Restore tabstops
    call settabs( '5-120+5' )
  end subroutine OUTPUTCALENDAR

  ! ---------------------------------------  OUTPUT_DATE_AND_TIME  -----
  subroutine OUTPUT_DATE_AND_TIME ( date, time, &
    & from_where, msg, dateFormat, timeFormat, &
    & CPU_Seconds, wallClock_Seconds, advance )
    ! Output nicely-formatted date, time, and extra message
    ! We'll assume we won't want this line stamped with date and time
    ! (for fear of being redundant, which we fear)
    ! Optionally print CPU and wall clock usage, too.
    logical, intent(in), optional :: date ! output date as character string
    logical, intent(in), optional :: time ! output time as character string
    character(len=*), intent(in), optional :: FROM_WHERE
    character(len=*), intent(in), optional :: MSG
    character(len=*), intent(in), optional :: DATEFORMAT
    character(len=*), intent(in), optional :: TIMEFORMAT
    double precision, intent(in), optional :: CPU_Seconds
    integer, intent(in), optional          :: wallClock_Seconds
    character(len=*), intent(in), optional :: ADVANCE

    character(len=16) :: dateString
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    integer :: HH, MM, MS, SS
    logical :: myDate
    logical :: myTime
    character(len=3) :: MY_ADV
    real :: My_CPU, Seconds
    character(len=16) :: timeString

    myDate = .true.
    if ( present(date) ) myDate = date
    myTime = .true.
    if ( present(time) ) myTime = time
    if ( .not. (myDate .or. myTime) ) return ! Why call if won't print?
    my_adv = 'no'
    if ( .not. present(msg) ) then
      my_adv = 'yes'
      if ( present(advance) ) my_adv = advance
    end if
    dateString = '' ! Intel 12 and earlier doesn't fill with blanks
    timeString = '' ! Intel 12 and earlier doesn't fill with blanks
    call date_and_time ( date=dateString, time=timeString )
    dateString = reFormatDate(trim(dateString), toForm=dateFormat)
    timeString = reFormatTime(trim(timeString), timeFormat)
    if ( myDate .and. myTime ) then
      call output_ ( trim(dateString), from_where=from_where, advance='no', &
        & DONT_STAMP=DONT_STAMP )
      call blanks(3)
      call output_ ( trim(timeString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    else if ( myDate ) then
      call output_ ( trim(dateString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    else if ( myTime ) then
      call output_ ( trim(TimeString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    end if
    if ( .not. present(msg) ) return
    my_adv = 'yes'
    if ( present(advance) ) my_adv = advance
    call blanks ( 3 )
    call output_ ( trim(msg), from_where=from_where, &
      & advance=merge ( 'no ', my_adv, present(CPU_seconds) .or. &
      & present(WallClock_seconds) ), &
      & DONT_STAMP=DONT_STAMP )
    if ( present(CPU_seconds) ) then
      hh = CPU_seconds / 3600
      my_cpu = CPU_seconds - hh * 3600
      mm = my_cpu / 60
      seconds = my_cpu - mm * 60
      ss = seconds
      ms = nint(1000 * (seconds-ss))
      write ( timeString, '(2(i2.2,":"),i2.2,".",i3.3)' ) hh, mm, ss, ms
      call output_ ( '   CPU time ' // trim(timeString), from_where=from_where, &
        & advance=merge ( 'no ', my_adv, present(WallClock_seconds)), &
        & DONT_STAMP=DONT_STAMP )
    end if
    if ( present(WallClock_seconds) ) then
      hh = WallClock_seconds / 3600
      my_cpu = WallClock_seconds - hh * 3600
      mm = my_cpu / 60
      seconds = my_cpu - mm * 60
      ss = seconds
      write ( timeString, '(2(i2.2,":"),i2.2)' ) hh, mm, ss
      call output_ ( '   wall clock time ' // trim(timeString), from_where=from_where, &
        & advance=my_adv, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_DATE_AND_TIME

  ! ----------------------------------------------  OUTPUTLIST  -----
  ! This family of routines outputs an array as a comma-separated list
  ! E.g., given the array (/ 1, 2, 3, .. /) outputs
  ! '(1, 2, 3, .. )'
  ! optionally using sep instead of ',' and delims instead of '()'
  subroutine OUTPUTLIST_CHARS ( array, sep, delims )
    ! Args
    character(len=*), dimension(:), intent(in)      :: array
    character(len=*), optional, intent(in) :: sep
    character(len=*), optional, intent(in) :: delims
    ! Local variables
    character(len=1) :: comma
    integer          :: i
    character(len=2) :: parens
    ! Executable
    if ( size(array) < 1 ) return
    comma = ','
    if ( present(sep) ) comma = sep
    parens = '()'
    if ( present(delims) ) parens = delims
    call output( parens(1:1) )
    do i=1, size(array)
      call output( trim_safe(array(i)) )
      if ( i < size(array) ) call output( comma )
    end do
    call output( parens(2:2) )
  end subroutine OUTPUTLIST_CHARS

  subroutine OUTPUTLIST_INTS ( array, sep, delims )
    ! Args
    integer, dimension(:), intent(in)      :: array
    character(len=*), optional, intent(in) :: sep
    character(len=*), optional, intent(in) :: delims
    ! Local variables
    character(len=1) :: comma
    integer          :: i
    character(len=2) :: parens
    ! Executable
    if ( size(array) < 1 ) return
    comma = ','
    if ( present(sep) ) comma = sep
    parens = '()'
    if ( present(delims) ) parens = delims
    call output( parens(1:1) )
    do i=1, size(array)
      call output( array(i) )
      if ( i < size(array) ) call output( comma )
    end do
    call output( parens(2:2) )
  end subroutine OUTPUTLIST_INTS

  ! ----------------------------------------------  outputNamedValue  -----
  ! This family of routines outputs a paired name and value
  ! (Basically saving you a few lines over the idiom
  !  call output ( trim(name), advance='no' )
  !  call output ( ': ', advance='no' )
  !  call output ( value, advance='yes' )
  
  ! to print following line to stdout
  !  name: value
  ! Optional args control
  ! before: what extra to print at start of each line
  ! after: what extra to print at end of each line
  ! colon: what to print instead of ':'
  ! fillChar: instead of spaces if you use tabs to align name, value
  ! tabn: column number where name begins
  ! tabc: column number where colon occurs
  ! taba: column number where after begins
  ! advance: whether to advance after printing pair (by default we WILL advance)
  ! dont_stamp: override setting to stamp end of each line
  ! By means of optional args you can create a line like
  ! *   name                   value   *
  ! See also startTable, addRow, outputTable
  subroutine output_nvp_whatever ( name, &
   & chvalue, ivalue, cmvalue, dbvalue, snvalue, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    character(len=*), intent(in), optional:: chvalue
    complex, intent(in), optional         :: cmvalue
    double precision, intent(in), optional:: dbvalue
    integer, intent(in), optional         :: ivalue
    real, intent(in), optional            :: snvalue
    character(len=*), intent(in), optional :: ADVANCE
    character(len=1), intent(in), optional :: COLON
    character(len=1), intent(in), optional :: fillChar
    integer, intent(in), optional :: TABN
    integer, intent(in), optional :: TABC
    integer, intent(in), optional :: TABA
    logical, intent(in), optional :: DONT_STAMP
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=*), intent(in), optional :: options
    ! Local variables
    if ( present(chvalue) ) then
      call output_nvp_character ( name, chvalue, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    elseif ( present(cmvalue) ) then
      call output_nvp_complex ( name, cmvalue, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    elseif ( present(dbvalue) ) then
      call output_nvp_double ( name, dbvalue, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    elseif ( present(ivalue) ) then
      call output_nvp_integer ( name, ivalue, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    elseif ( present(snvalue) ) then
      call output_nvp_single ( name, snvalue, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    endif
  end subroutine output_nvp_whatever

  subroutine output_nvp_character ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    character(len=*), intent(in)          :: value
    character(len=*), intent(in), optional :: ADVANCE
    character(len=1), intent(in), optional :: COLON
    character(len=1), intent(in), optional :: fillChar
    integer, intent(in), optional :: TABN
    integer, intent(in), optional :: TABC
    integer, intent(in), optional :: TABA
    logical, intent(in), optional :: DONT_STAMP
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=*), intent(in), optional :: options
    if ( TrimCharacterValues ) then
      call possiblyTrimmedvalue ( name, trim(value), &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    else
      call possiblyTrimmedvalue ( name, value, &
        & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    endif
  contains
    subroutine possiblyTrimmedvalue ( name, value, &
      & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
      character(len=*), intent(in)          :: name
      character(len=*), intent(in)          :: value
      include 'output_name_value_pair.f9h'
    end subroutine possiblyTrimmedvalue
  end subroutine output_nvp_character

  subroutine output_nvp_complex ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    complex, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_complex

  subroutine output_nvp_double ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    double precision, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_double

  subroutine output_nvp_dbl_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    double precision, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_dbl_array

  subroutine output_nvp_int_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    integer, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_int_array

  subroutine output_nvp_integer ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_integer

  subroutine output_nvp_log_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    logical, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_log_array

  subroutine output_nvp_logical ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    logical, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_logical

  subroutine output_nvp_single ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_single

  subroutine output_nvp_sngl_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP, options )
    character(len=*), intent(in)          :: name
    real, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_sngl_array

  ! ----------------------------------------------  outputTable  -----
  ! Outputs a 2d character array as a table
  ! Optionally, 
  ! (1)  you may supply the arry; otherwise, dellDatabase will be used
  ! (2)  the table can have a character separating cells, and another 
  !      marking its outer borders
  ! (3)  the minimum cell width can be set; otherwise
  !      it is computed based on the trimmed lengths of each column
  ! (4)  the alignment within each cell can be set; otherwise
  !      it is flushed left, i.e. 'L'
  ! (5)  each row can be separated by an interior wall of characters;
  !      by default they are consecutive
  !      a special value of interior, null (achar(0)) inserts an empty line
  ! (6)  the first row can be treated as special, and separated from the second
  !      by a wall of special headliner characters
  subroutine outputTable ( array, sep, border, cellWidth, &
    & interior, headliner, alignment )
    ! Args
    character(len=*), dimension(:,:), optional, intent(in)   &
      &                                            :: array
    character(len=1), optional, intent(in)         :: sep       ! between cols
    character(len=1), optional, intent(in)         :: border    ! outside
    integer, optional, intent(in)                  :: cellWidth
    character(len=1), optional, intent(in)         :: interior  ! between rows
    character(len=1), optional, intent(in)         :: headliner ! 1st row are headers
    character(len=1), optional, intent(in)         :: alignment ! L, R, or C
    ! Local variables
    integer                                        :: status
    ! Executable
    if ( present( array ) ) then
      call outputTableArray ( array, sep, border, cellWidth, &
        & interior, headliner, alignment )
    elseif ( .not. associated ( cellDatabase ) ) then
      call banner ( 'Empty table' )
    else
      call outputTableArray ( cellDatabase, sep, border, cellWidth, &
        & interior, headliner, alignment )
      deallocate( cellDatabase, stat=status )
      nullify( cellDatabase )
    endif
  end subroutine outputTable

  ! ----------------------------------------------  resetTabs  -----
  ! Restore tab stops to what was in effect at start
  ! Optionally returning them as an integer array
  subroutine resetTabs ( tabs )
    ! Args
    integer, dimension(:), optional, intent(out) :: tabs
    ! Internal variables
    integer :: n
    ! Executable
    call setTabs( range=INITTABRANGE )
    if ( present(tabs) ) then
      n = min(size(tabs), MAXNUMTABSTOPS)
      tabs = 0
      tabs(1:n) = TABSTOPS(1:n)
    end if
  end subroutine resetTabs

  ! ----------------------------------------------  RestoreSettings  -----
  ! Restore tab stops to what was in effect at start
  ! Optionally returning them as an integer array
  subroutine RestoreSettings ( settings )
    ! Args
    character(len=*), optional, intent(in) :: settings
    ! Local variables
    character(len=*), parameter            :: allSettings = &
      & 'style'
    character(len=64)                      :: mySettings 
    ! Executable
    call RestoreOutputSettings ( settings )
    mySettings = ' '
    if ( present(settings) ) mySettings = settings
    if ( index(mySettings, '*') > 0 ) mySettings = allSettings
    mySettings = lowercase(mySettings)
    if ( index(mySettings, 'style') > 0 ) StyleOptions = DefaultStyleOptions
  end subroutine RestoreSettings 

  ! ----------------------------------------------  setStamp  -----
  subroutine setStamp ( textCode, showTime, dateFormat, timeFormat, &
    & post, interval )
  ! set stamp to be added to every output to PRUNIT.
    character(len=*), optional, intent(in) :: textCode
    logical, optional, intent(in)          :: showTime
    character(len=*), optional, intent(in) :: dateFormat
    character(len=*), optional, intent(in) :: timeFormat
    logical, optional, intent(in)          :: post
    integer, optional, intent(in)          :: interval
    if ( present(showTime) )   stampOptions%showTime   = showTime  
    if ( present(textCode) )   stampOptions%textCode   = textCode  
    if ( present(dateFormat) ) stampOptions%dateFormat = dateFormat
    if ( present(timeFormat) ) stampOptions%timeFormat = timeFormat
    if ( present(post) )       stampOptions%post       = post
    if ( present(interval) )   stampOptions%interval   = interval
  end subroutine setStamp

  ! ----------------------------------------------  setTabs  -----
  subroutine setTabs ( RANGE, TABS )
    ! Set tabstops
    ! Methods:
    ! (1) a string range; e.g., "8, 32-100+8"
    !     is converterd to "8, 32, 40, 48, 56, 64, 72, 80, 88, 96"
    ! (2) an arrays of ints; e.g., (/ 4, 9, 12, 18, 22, 30, 35, 40 /)
    ! (3) reset back to the defaults (equiv to "5-120+5")
    ! Args
    character(len=*), optional, intent(in)         :: Range
    integer, dimension(:), optional, intent(in)    :: Tabs
    ! Internal variables
    integer :: n
    ! Executable
    if ( present(range) ) then
      call ExpandStringRange ( range, TABSTOPS )
    else if ( present(tabs) ) then
      n = min( MAXNUMTABSTOPS, size(tabs) )
      tabStops(1:n) = tabs(1:n) 
    else
      call ExpandStringRange ( '5-120+5', TABSTOPS )
    end if
  end subroutine setTabs

  ! ----------------------------------------------  startTable  -----
  ! Set up a table, initially empty. Subsequent calls to AddRow will
  ! add a new row, consisting of two columns: 
  ! the name field and the value field
  ! Optionally, 
  ! (1)  the table can have a character separating cells, and another 
  !      marking its outer borders
  ! (2)  the minimum cell width can be set; otherwise
  !      it is computed based on the trimmed lengths of each column
  ! (3)  the alignment within each cell can be set; otherwise
  !      it is flushed left, i.e. 'L'
  ! (4)  each row can be separated by an interior wall of characters;
  !      by default they are consecutive
  !      a special value of interior, null (achar(0)) inserts an empty line
  ! (5)  the first row can be treated as special, and separated from the second
  !      by a wall of special headliner characters
  subroutine startTable 
    ! Internal variables
    integer :: status
    ! Executable
    status = 0
    if ( associated( cellDatabase ) ) &
      & deallocate( cellDatabase, stat=status )
    if ( status /= 0 ) call myMessage ( MLSMSG_Warning, 'startTable', &
      & 'Unable to deallocate celldatabase' )
    nullify ( cellDatabase )
    outputLines = ' '
  end subroutine startTable

  ! -----------------------------------------------------  StyledOutput  -----
  ! Print your message according to the style set by options
  ! e.g., "--Banner"
  ! See also banner
  !
  ! We could dig more deeply into options to allow\
  ! it to pass Alignment, Fill, After, etc.
  subroutine StyledOutput ( chars, options )
    character(len=*), intent(in)                :: chars
    character(len=*), intent(in), optional      :: options
    ! Internal variables
    logical :: asBanner
    logical :: asHeadline
    ! Executable
    asBanner   = .false.
    asHeadline = .false.
    if ( .not. present(options ) ) then
      call output( chars, advance='yes' )
      return
    endif
    asBanner   = index( options, 'B' ) > 0
    asHeadline = index( options, 'H' ) > 0
    if ( asBanner ) then
      call Banner( chars ) 
    elseif ( asHeadline ) then
      StyleOptions%HeadLineFill = '-'
      StyleOptions%HeadLineBefore = '*'
      StyleOptions%HeadLineAfter = '*'
      call Headline( chars ) 
      StyleOptions = DefaultStyleOptions
    else
      call output( chars, advance='yes' )
    endif
  end subroutine StyledOutput

  ! ------------------------------------------------  timeStamp  -----
  ! time-stamp output on demand, not automatic:
  ! Either in style pre or post
  ! (pre) '(HH:MM:SS) chars'
  ! (post) 'chars (HH:MM:SS)'
  ! Note that in pre-style, the time will be printed only if getOutputStatus( 'start' ) == 1 true
  ! in post-style, the time will be printed only if MY_ADV is 'yes'
  subroutine timeStamp_char ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, STYLE, DATE )
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    character(len=8), intent(in), optional :: STYLE ! pre or post
    logical, intent(in), optional          :: DATE  ! Include date with time
    !
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    character(len=8) :: my_style
    character(len=3) :: MY_ADV
    logical  ::         myDate
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_style = timeStampOptions%Timestampstyle
    if ( present(style) ) my_style = lowercase(style)
    myDate = timeStampOptions%showDate
    if ( present(date) ) myDate = date
    if ( my_style == 'post' ) then
      call output_( CHARS, &
        & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
        & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK, DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_(' (', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME( date=myDate, &
          & dateFormat=timeStampOptions%dateFormat, &
          & timeFormat=timeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='yes', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      end if
    else
      if ( getOutputStatus( 'start' ) == 1 ) then
        call output_('(', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME( date=myDate, &
          & dateFormat=timeStampOptions%dateFormat, &
          & timeFormat=timeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      end if
      call output_( CHARS, &
        & ADVANCE, FROM_WHERE, DONT_LOG, &
        & LOG_CHARS, INSTEADOFBLANK, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine timeStamp_char

  subroutine timeStamp_integer ( INT, &
    & PLACES, ADVANCE, FILL, FORMAT, Before, After, style, date )
    integer, intent(in) :: INT
    integer, intent(in), optional :: PLACES
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: FILL
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=*), intent(in), optional :: style ! pre or post
    logical, intent(in), optional          :: DATE  ! Include date with time
    !
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    character(len=8) :: my_style
    character(len=3) :: MY_ADV
    logical  ::         myDate
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_style = timeStampOptions%Timestampstyle
    if ( present(style) ) my_style = lowercase(style)
    myDate = timeStampOptions%showDate
    if ( present(date) ) myDate = date
    if ( my_style == 'post' ) then
      call output( INT, PLACES, &
        & ADVANCE='no', FILL=FILL, FORMAT=FORMAT, BEFORE=BEFORE, AFTER=AFTER, &
        & DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_(' (', ADVANCE='no', DONT_STAMP=DONT_STAMP )
        call OUTPUT_DATE_AND_TIME( date=myDate, &
          & dateFormat=timeStampOptions%dateFormat, &
          & timeFormat=timeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='yes', DONT_STAMP=DONT_STAMP)
      end if
    else
      if ( getOutputStatus( 'start' ) == 1 ) then
        call output_('(', ADVANCE='no', DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME( date=myDate, &
          & dateFormat=timeStampOptions%dateFormat, &
          & timeFormat=timeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='no', DONT_STAMP=DONT_STAMP)
      end if
      call output( INT, PLACES, &
        & ADVANCE, FILL, FORMAT, BEFORE, AFTER, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine timeStamp_integer

  subroutine timeStamp_logical ( value, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, STYLE, DATE )
    logical, intent(in) ::                    value
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    character(len=8), intent(in), optional :: STYLE ! pre or post
    logical, intent(in), optional          :: DATE  ! Include date with time
    ! Internal variables
    character(len=1) :: str
    ! Executable
    if ( value ) then
      str = 'T'
    else
      str = 'F'
    end if
    call timeStamp_char(str, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, STYLE, DATE )
  end subroutine timeStamp_logical

  ! ------------------ Private procedures -------------------------
  ! .............................................  addCellRowToDatabase  .....
  ! Add an additional row to end of cellDatabase
  integer function addCellRowToDatabase ( database, item )
    ! Args
    character(len=*), dimension(:,:), pointer  :: database
    character(len=*), dimension(2)             :: item
    ! Local variables
    integer :: newSize
    integer :: status
    character(len=MAXCELLSIZE), dimension(:,:), pointer  &
      &                                        :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    if ( associated(database) ) then ! tree_checker prevents duplicate names
      newSize = SIZE(database,1) + 1
    else
      newSize = 1
    end if
    allocate ( tempDatabase(newSize,2), STAT=status )
    if ( newSize>1 ) tempDatabase(1:newSize-1,:) = database
    if ( associated(database) ) then
      deallocate ( database, stat=status )
    end if
    database => tempDatabase
    database(newSize,:) = ' '
    database(newSize,:) = adjustl(item)

    ! include "addCellRowToDatabase.f9h" 

    addCellRowToDatabase = newSize
  end function addCellRowToDatabase

  ! .............................................  getOption  .....
  ! This family of subroutines parses a multipart advance arg into
  ! its components, returning an appropriate value
  ! Example, say the component is marked by the '-S' flag
  ! value type     component   returned value
  !  logical          -S         true
  ! character       -Sxyz        xyz
  subroutine getOption_char ( arg, flag, val )
    ! Args
    character(len=*), intent(in)  :: arg
    character(len=*), intent(in)  :: flag
    character(len=*), intent(out) :: val
    ! Local variables
    integer :: kFlag, kNext, flagLen
    ! Executable
    val = ' '
    kFlag = index( arg, trim(flag) )
    if ( kFlag < 1 ) return
    flagLen = len(flag)
    ! Find start of next component
    kNext = index( arg(kFlag+flagLen:), ' ' )
    if ( kNext < 1 ) then
      val = arg(kFlag+flagLen:)
    else
      val = arg(kFlag+flagLen:kFlag+flagLen+kNext)
    end if
  end subroutine getOption_char

  subroutine getOption_log ( arg, flag, val )
    ! Args
    character(len=*), intent(in) :: arg
    character(len=*), intent(in) :: flag
    logical, intent(out)         :: val
    ! Local variables
    integer :: kFlag
    ! Executable
    kFlag = index( arg, trim(flag) )
    val = ( kFlag > 0 )
  end subroutine getOption_log

  ! ------------------------------------  myMessage  -----
  subroutine myMessage ( Severity, ModuleNameIn, Message, &
    & Advance )

    ! Print a message (unless printing is suppressed).  If it has %[Nn]
    ! in it, replace that with newline.

    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'

    ! Local variables
    !                                     If nonzero, do not insert prefix.
    logical :: AllOfIt                  ! Print all of it (no %n or %N remains)
    integer :: L1, L2                   ! How far in the line have we printed?
    character (len=512), save :: Line   ! Line to output, should be long enough
    integer, save :: Line_len=0         ! Number of saved characters in line.
    integer :: LogFileUnit
    logical :: My_adv

    ! Executable code
    my_adv = .true.
    if ( present(advance) ) &
      & my_adv = advance(1:1) /= 'n' .and. advance(1:1) /= 'N'

    my_adv = my_adv .and. ( severity >= MLSMessageConfig%skipMessageThr )
    if ( (.not. MLSMessageConfig%suppressDebugs).OR. &
         & (severity /= MLSMSG_Debug) ) then
      l1 = 0
      do
        l2 = index(Message(l1+1:),'%n' )
        if ( l2 == 0 ) l2 = index(Message(l1+1:),'%N')
        allOfIt = l2 == 0
        l2 = l2 + l1 - 1 ! Last character before %n or %N, if any
        if ( allOfIt ) l2 = len(Message) ! no %n or %N
        call assembleFullLine( Severity, ModuleNameIn, Message(l1+1:l2), &
          & line, line_len )
        l1 = l2 + 2 ! "n" or "N" of %n or %N
         ! Log the message using the toolkit routine
         ! (or its substitute )
         ! if either using toolkit or severity is sufficient to
         ! quit (which means we might have been called directly
         ! rather than from output module )

        if ( my_adv .or. .not. allOfIt ) then
          call printitout( line, severity, line_len )
          line_len = 0
          line = ' '
        end if
        if ( allOfIt ) exit
      end do

    end if

    ! Now if it's an error, and the message is complete, then try to close
    ! log file if any and quit (or crash)

    if ( my_adv .and. severity >= MLSMSG_Severity_to_quit ) then
      call get_config ( logFileUnit = logFileUnit )
      if ( logFileUnit > 0 ) close ( logFileUnit  )
      if ( severity >= MLSMSG_Crash .or. MLSMessageConfig%CrashOnAnyError ) then
        NEVERCRASH = .false.
        call crash_burn
      end if
      call exit_with_status ( 1  )
    end if
  end subroutine myMessage

  ! ----------------------------------------------  outputTableArray  -----
  ! outputs a 2d character array as a table
  ! Optionally, 
  ! (1)  the table can have a character separating cells, and another 
  !      marking its outer borders
  ! (2)  the minimum cell width can be set; otherwise
  !      it is computed based on the trimmed lengths of each column
  ! (3)  the alignment within each cell can be set; otherwise
  !      it is flushed left, i.e. 'L'
  ! (4)  each row can be separated by an interior wall of characters;
  !      by default they are consecutive
  !      a special value of interior, null (achar(0)) inserts an empty line
  ! (5)  the first row can be treated as special, and separated from the second
  !      by a wall of special headliner characters
  ! (6)  If any row begins with the special formatting sequence '<<', then
  !      one of the following occurs:
  !      special format             action
  !          <<l>>xxxx       merge:  left-aligned xxx stretched across table
  !          <<c>>xxxx       merge:  centered xxx stretched across table
  !          <<r>>xxxx       merge:  right-aligned xxx stretched across table
  !          <<d>>x          divider: inserts a wall of x's
  subroutine outputTableArray ( array, sep, border, cellWidth, &
    & interior, headliner, alignment )
    ! use Dump_0, only: Dump
    ! Args
    character(len=*), dimension(:,:),intent(in)    :: array
    character(len=1), optional, intent(in)         :: sep       ! between cols
    character(len=1), optional, intent(in)         :: border    ! outside
    integer, optional, intent(in)                  :: cellWidth
    character(len=1), optional, intent(in)         :: interior  ! between rows
    character(len=1), optional, intent(in)         :: headliner ! 1st row are headers
    character(len=1), optional, intent(in)         :: alignment ! L, R, or C
    ! Local variables
    character(len=1)                               :: align
    character(len=MAXCELLSIZE)                     :: cell
    integer                                        :: i
    integer                                        :: j
    integer                                        :: k
    integer                                        :: left
    integer                                        :: minWidth
    character(len=1)                               :: myAlignment
    character(len=1)                               :: myBorder
    character(len=1)                               :: mySep
    character(len=1)                               :: myHeadliner
    character(len=1)                               :: myInterior
    integer                                        :: right
    integer, dimension(size(array,2))              :: widths
    integer, parameter                             :: leftPadding = 1
    integer, parameter                             :: rightPadding = 1
    integer                                        :: totalWidth
    logical, parameter                             :: debug=.false.
    ! Executable
    minWidth = 3 ! Don't know why, but this works
    if ( present(cellWidth) ) minWidth = cellWidth

    mySep = ' '
    if ( present(sep) ) mySep = sep
    myBorder = ' '
    if ( present(border) ) myBorder = border
    myInterior = ' '
    if ( present(Interior) ) myInterior = Interior
    
    myHeadliner = myInterior
    if ( present(headliner) ) myHeadliner = headliner
    
    myAlignment = 'L'
    if ( present(alignment) ) myAlignment = alignment

    ! 1st, compute total table width
    widths = minWidth
    totalWidth = 3
    do j=1, size(array,2)
      widths(j) = maxval( len_trim(array(:,j)) )
      widths(j) = max( widths(j), minWidth )
      totalWidth = totalWidth + widths(j) + leftPadding + rightPadding
      if ( j > 1 .and. j < size(array,2) ) totalWidth = totalWidth + 1
    enddo
    if ( debug ) then
      ! call dump( widths, 'widths' )
      call outputNamedValue( 'totalWidth', totalWidth )
      do i=1, size(array,1)
        call output ( i, advance='no' )
        call blanks(2)
        call output( trim(array(i, 1)), advance='yes' )
      enddo
      do i=1, size(array,1)
        call output ( i, advance='no' )
        call blanks(2)
        call output( trim(array(i, 2)), advance='yes' )
      enddo
    endif
    if ( len_trim(myBorder) > 0 ) &
      & call output( repeat( myBorder, totalWidth ), advance='yes' )
    do i=1, size(array,1)
      right = 0
      if ( len_trim(myBorder) > 0 ) then
        call output( myBorder, advance='no' )
        right = right + 1
      endif
      ! Special formatting?
      k = index(array(i,1), '<<' )
      if ( k > 0 ) then
        align = array(i,1)(k+2:k+2)
        select case ( lowercase(align))
        case ('l') ! left-aligned merge
          call blanks ( leftPadding )
          call output( trim(array(i,1)(k+5:)), advance='no' )
        case ('c') ! centered merge
          call aligntofit ( trim(array(i,1)(k+5:)), &
            & (/ right, totalWidth /), 'c' )
        case ('r') ! right-aligned merge
          call aligntofit ( trim(array(i,1)(k+5:)), &
            & (/ right, totalWidth /), 'r' )
        case ('d') ! wall of chars
          call blanksToColumn( totalWidth, fillChar=array(i,1)(k+5:k+5) )
        end select
        if ( len_trim(myBorder) > 0 ) then
          call blanksToColumn( totalWidth )
          call output( myBorder, advance='yes' )
        else
          call newLine
        endif
        cycle
      endif
      do j=1, size(array,2)
        left = right + 1 + leftPadding
        right = left + widths(j) ! Don't know why, but this works
        call alignToFit ( trim(array(i,j)), (/ left, right /), myAlignment )
        call blanks ( rightPadding )
        if ( len_trim(mySep) > 0 .and. j < size(array,2) ) then
          call output( mySep, advance='no' )
          right = right + 1
        endif
      enddo
      if ( len_trim(myBorder) > 0 ) then
        call blanksToColumn( totalWidth )
        call output( myBorder, advance='yes' )
      else
        call newLine
      endif
      ! Interior cell walls or headliners
      if ( len_trim(myheadliner) > 0 .and. i == 1 .and. &
        & myheadliner /= achar(0) ) then
        call output( repeat( myheadliner, totalWidth ), advance='yes' )
      elseif ( myInterior == achar(0) .and. i < size(array,1) ) then
        call output( myBorder, advance='no' )
        call blanksToColumn( totalWidth )
        call output( myBorder, advance='yes' )
      elseif ( len_trim(myInterior) > 0 .and. i < size(array,1) ) then
        call output( repeat( myInterior, totalWidth ), advance='yes' )
      endif
    enddo
    if ( len_trim(myBorder) > 0 ) &
      & call output( repeat( myBorder, totalWidth ), advance='yes' )
  end subroutine outputTableArray

  ! ------------------------------------  SeparateElements  -----
  ! insert blanks or separator between consecutive elements while outputting
  subroutine SeparateElements (i, n )
    ! Args
    integer, intent(in) :: i ! Element number
    integer, intent(in) :: n ! Number of elements
    ! Executable
    if ( i >= n ) return
    if ( wrappastcolnum > 0 .and. getOutputStatus( 'column' ) >= wrappastcolnum ) then
      call newLine
      return
    end if
    if ( wrappastcolnum == 0 .and. &
      & mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
      call output_ ( '', advance='yes', DONT_STAMP=.true. )
      return
    end if
    if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) then
      call output_( outputOptions%arrayElmntSeparator, advance='no' )
    else
      call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
    end if
  end subroutine SeparateElements

  ! ----------------------------------------------  stamp  -----
  function stamp( chars )
  ! stamp input chars before outputting to PRUNIT.
  ! Args
    character(len=*), intent(in) :: chars
    character(len=len(chars)+64) :: stamp
    character(len=16) :: dateString
    character(len=16) :: timeString
    ! Executable
    stamp = chars
    if ( stampOptions%showTime ) then
      dateString = '' ! Intel 12 and earlier doesn't fill with blanks
      timeString = '' ! Intel 12 and earlier doesn't fill with blanks
      call date_and_time ( date=dateString, time=timeString )
      dateString = reFormatDate(trim(dateString), toForm=stampOptions%dateFormat)
      timeString = reFormatTime(trim(timeString), stampOptions%timeFormat)
      if ( stampOptions%dateFormat /= ' ' ) &
      & stamp = catStrings( stamp, dateString )
      if ( stampOptions%timeFormat /= ' ' ) &
      & stamp = catStrings( stamp, timeString )
    end if
    if ( stampOptions%textCode /= ' ' ) &
        & stamp = catStrings( stamp, stampOptions%textCode )
  contains
    function catStrings(a, b) result(c)
      ! Catenates strings a and b with intervening space
      ! if post then a before b
      ! otherwise b before a
      character(len=*), intent(in) :: a
      character(len=*), intent(in) :: b
      character(len = (len(a)+len(b)+1) ) :: c
      if ( stampOptions%post ) then
        c = trim(a) // ' ' // b
      else
        c = trim(b) // ' ' // a
      end if
    end function catStrings
    
  end function stamp 

  ! ----------------------------------------------  stretch  -----
  function stretch( arg, skips ) result(chars)
  ! stretch input arg by inserting skips number of spaces
  ! between each pair of consecutive characters
  ! Args
    character(len=*), intent(in)      :: arg
    integer, intent(in)               :: skips
    character(len=(1+skips)*len(arg)) :: chars
    ! Internal variables
    integer :: i, k
    ! Executable
    chars = ' '
    if ( len_trim(arg) < 1 ) return
    do i=1, len_trim(arg)
      ! E.g., if skips==1, k ~ 1 3 5 7 ..
      k = 1 + (skips+1)*(i-1)
      chars(k:k) = arg(i:i)
    end do
  end function stretch

  ! ........................................  whatSDNeedsFormat
  ! parse inFormat which might be
  ! (1) absent, in which case format=sdNeedsFormat and dotm='.6'
  ! (2) '(*)', in which case format=sdNeedsFormat and dotm='.6'
  ! (3) '(*.m)', in which case format=(sdNeedsFragment //'.m') and dotm='.m'
  subroutine whatSDNeedsFormat ( format, dotm, inFormat )
    character(len=*), optional, intent(in)  :: inFormat
    character(len=*), intent(out)           :: format
    character(len=*), intent(out)           :: dotm
    integer :: dot
    if ( .not. present(inFormat) ) then
      format = sdNeedsFormat
      dotm = '.6'
    else if ( index(inFormat, '.') < 1 ) then
      format = sdNeedsFormat
      dotm = '.6'
    else
      ! Must find integer after '.'
      dot = index( inFormat, '.' )
      dotm = inFormat(dot:dot+1)
      format = trim(sdNeedsFragment) // trim(dotm) // ')'
    end if
  end subroutine whatSDNeedsFormat

  ! ..............................................  not_used_here  .....
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HIGHOUTPUT

! $Log$
! Revision 2.22  2017/12/14 23:20:48  pwagner
! Added TrimCharacterValues
!
! Revision 2.21  2017/11/30 20:48:59  pwagner
! RestoreSettings may now restore all or just some
!
! Revision 2.20  2017/11/15 00:01:51  pwagner
! Print options%neverStamp as part of Dump
!
! Revision 2.19  2017/10/03 21:44:11  pwagner
! restored bars to Banner; improved comments showing effect of BannerPattern
!
! Revision 2.18  2017/09/29 00:20:27  pwagner
! Added Styled output and options; options as an arg to OutputNamedValue
!
! Revision 2.17  2017/09/07 23:43:40  pwagner
! Improved comments; removed unused myMessage_old
!
! Revision 2.16  2017/01/19 23:33:23  pwagner
! New procedures to add an internal wall or merge cells across a cell table
!
! Revision 2.15  2017/01/13 01:28:32  pwagner
! Rename internal procedure to addCellRowToDatabase
!
! Revision 2.14  2016/12/14 01:22:40  pwagner
! outputCalendar prints new, full moons if moonphases present and ttrue
!
! Revision 2.13  2016/11/15 19:27:19  pwagner
! May print elapsed WallClock_seconds at Finish
!
! Revision 2.12  2016/09/22 22:21:16  pwagner
! May specify format in call to addRow
!
! Revision 2.11  2016/03/25 00:37:06  pwagner
! Added OUTPUTANYNAMEDVALUE
!
! Revision 2.10  2015/09/24 18:50:41  pwagner
! May choose different pattern for stripes in banner
!
! Revision 2.9  2015/05/18 17:42:50  pwagner
! addRow and startTable maintains an internal Table for outputTable to output
!
! Revision 2.8  2015/02/24 23:32:22  pwagner
! Make sure rightpadding defined in headLine
!
! Revision 2.7  2015/02/13 00:17:49  pwagner
! Added procedure to output 2d array as Table
!
! Revision 2.6  2015/02/10 01:00:58  pwagner
! Avoid another double-indent error
!
! Revision 2.5  2015/02/06 01:08:20  pwagner
! Can now print to virtual page indented w.r.t. physical page
!
! Revision 2.4  2014/10/06 23:06:26  vsnyder
! Add Signed argument to Dump_Size
!
! Revision 2.3  2014/09/05 00:23:57  vsnyder
! Better handling of literal output unit.  Convert a local pointer temp to
! allocatable.
!
! Revision 2.2  2014/04/22 16:30:57  pwagner
! Banner accepts mode as optional arg
!
! Revision 2.1  2014/01/09 00:22:32  pwagner
! Split output module procedures between it and new highOutput
!
