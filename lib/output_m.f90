! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module OUTPUT_M

  ! Very high level printing and formatting
  ! see also dump_0 and MLSMessageModule
  
  use DATES_MODULE, only:  BUILDCALENDAR, DAYSINMONTH, &
    & REFORMATDATE, REFORMATTIME, UTC_TO_YYYYMMDD
  use MLSCOMMON, only: FILENAMELEN, FINITE_SIGNAL, IS_WHAT_IEEE
  use MLSMESSAGEMODULE, only: MLSMESSAGE, &
    & MLSMSG_INFO, MLSMSG_ERROR
  use MLSSETS, only: FINDFIRST
  use MLSSTRINGLISTS, only: EXPANDSTRINGRANGE, GETSTRINGELEMENT, &
    & LIST2ARRAY, NUMSTRINGELEMENTS, WRAP
  use MLSSTRINGS, only: REPLACENONASCII, LOWERCASE, NCOPIES, &
    & READINTSFROMCHARS, TRIM_SAFE, WRITEINTSTOCHARS
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! outputLines              If PrUnit = OUTPUTLINESPRUNIT, holds output until flushed
! outputOptions            where to send output and how to format it
! (some components)
! MLSMSG_Level             MLSMessage level if so logged
! PrUnit                   How to direct output:
!                          OUTPUTLINESPRUNIT :: held in outputLines
!                          MSGLOGPRUNIT      :: logged via MLSMessage
!                          STDOUTPRUNIT      :: print stdout
!                          BOTHPRUNIT        :: both stdout and logged
!                          > 0 :: print to Fortran unit number PrUnit
! skipMLSMSGLogging        whether to skip MLSMessage by default
! stampOptions             whether and how to stamp each output automatically
! timeStampOptions         how to stamp when calling timeStamp

!     (subroutines and functions)
! alignToFit               align printed argument to fit column range
! banner                   surround message with stars and stripes
! blanks                   print specified number of blanks [or fill chars]
! blanksToColumn           print blanks [or fill chars] out to specified column
! blanksToTab              print blanks [or fill chars] out to next tab stop
! dump                     dump output or stamp options
! dumpsize                 print a nicely-formatted memory size 
! dumptabs                 print the current tab stop positions
! flushOutputLines         print the current outputLines; then reset to ''
! getStamp                 get stamp being added to every output
! headLine                 print a line with extra formatting features
! isOutputSuspended         returns TRUE if output is suspended
! newline                  print a newline
! numNeedsFormat           return what format is need to output num
! numToChars               return what would be printed by output
! output                   print argument
! outputCalendar           output nicely-formatted calendar page
! output_date_and_time     print nicely formatted date and time
! outputList               output array as comma-separated list; e.g. '(1,2,..)'
! outputNamedValue         print nicely formatted name and value
! resettabs                restore tab stops to what was in effect at start
! restoreSettings          restore default settings for output, stamps, tabs
! revertOutput             revert output to file used before switchOutput
!                           if you will revert, keepOldUnitOpen when switching
! resumeOutput             resume output
! setStamp                 set stamp to be added to every output automatically
! setTabs                  set tab stops (to be used by tab)
! suspendOutput            suspend output
! switchOutput             switch output to a new named file
! tab                      move to next tab stop
! timestamp                print argument with a timestamp manually
!                            (both stdout and logged output)
! === (end of toc) ===

! === (start of api) ===
! alignToFit ( char* chars, int columnRange(2), char alignment, [int skips] )
! banner ( char* chars, [int columnRange(2)], [char alignment], [int skips] )
! blanks ( int n_blanks, [char fillChar], [char* advance] )
! blanksToColumn ( int column, [char fillChar], [char* advance] )
! blanksToTab ( [int tabn], [char* fillChar] )
! Dump ( options )
! DumpSize ( n, [char* advance], [units] )
!       where n can be an int or a real, and 
!       units is a scalar of the same type, if present
! DumpTabs ( [int tabs(:)] )
! flushOutputLines ( [int prUnit] )
! getStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! headLine ( char* chars, 
!          [char fillChar], [char* Before], [char* After], 
!          [int columnRange(2)], [char alignment], [int skips] )
! log isOutputSuspended ()
! NewLine
! char* numNeedsFormat ( value )
! char* numToChars ( value, [char* format] )
! output ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [log dont_stamp], [int newlineval] )
! output ( char* chars(:), [char* advance],
!          [char* insteadOfBlank], [int newlineval] )
! output ( value, [char* format], [char* advance],
!          [char* Before], [char* After] )
!       where value can be any numerical type, either scalar or 1-d array
! output_date_and_time ( [log date], [log time], [char* from_where], 
!          [char* msg], [char* dateFormat], [char* timeFormat], 
!          [double CPU_seconds], [char* advance] )
! outputCalendar ( [char* date], [char* datenote], [char* notes(:)], 
!          [dontwrap] )
! outputList ( values(:), [char* sep], [char* delims] )
! outputNamedValue ( char* name, value, [char* advance],
!          [char colon], [char fillChar], [char* Before], [char* After], 
!          [integer tabn], [integer tabc], [integer taba], log dont_stamp] )
! resetTabs ( [int tabs(:)] )
! resumeOutput
! revertOutput
! restoreSettings
! setStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! setTabs ( [char* Range], [int tabs(:)] )
! suspendOutput
! switchOutput ( char* filename, [int unit] )
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
  integer, save, public :: LINE_WIDTH = 120 ! Not used here, but a convenient
                                            ! place to store it
  ! Where to output?
  ! These apply if we don't output to a fortran unit number (which is > 0)
  integer, parameter, public :: INVALIDPRUNIT      = 0
  integer, parameter, public :: STDOUTPRUNIT       = INVALIDPRUNIT - 1
  integer, parameter, public :: MSGLOGPRUNIT       = STDOUTPRUNIT - 1
  integer, parameter, public :: BOTHPRUNIT         = MSGLOGPRUNIT - 1
  integer, parameter, public :: OUTPUTLINESPRUNIT  = BOTHPRUNIT - 1

  integer, save, private :: OLDUNIT = -1 ! Previous Unit for output.
  logical, save, private :: OLDUNITSTILLOPEN = .TRUE.

  public :: ALIGNTOFIT, BANNER, BLANKS, BLANKSTOCOLUMN, BLANKSTOTAB, &
    & DUMP, DUMPSIZE, DUMPTABS, FLUSHOUTPUTLINES, GETSTAMP, HEADLINE, &
    & NEXTCOLUMN, NEXTTAB, NEWLINE, NUMNEEDSFORMAT, NUMTOCHARS, &
    & OUTPUT, OUTPUT_DATE_AND_TIME, OUTPUTCALENDAR, OUTPUTLIST, &
    & OUTPUTNAMEDVALUE, &
    & RESETTABS, RESTORESETTINGS, RESUMEOUTPUT, REVERTOUTPUT, &
    & SETSTAMP, SETTABS, SUSPENDOUTPUT, SWITCHOUTPUT, TAB, TIMESTAMP

  ! These types made public because the class instances are public
  public :: outputOptions_T
  public :: stampOptions_T
  public :: timeStampOptions_T

  interface aligntofit
    module procedure aligntofit_chars, aligntofit_double, aligntofit_single
    module procedure aligntofit_integer
  end interface

  interface banner
    module procedure banner_chars
    module procedure banner_chararray
  end interface

  interface DUMP
    module procedure DUMPOUTPUTOPTIONS, DUMPSTAMPOPTIONS, DUMPTIMESTAMPOPTIONS
  end interface

  interface DUMPSIZE
    module procedure DUMPSIZE_DOUBLE, DUMPSIZE_INTEGER, DUMPSIZE_REAL
  end interface

  interface getOption
    module procedure getOption_char, getOption_log
  end interface

  interface NUMNEEDSFORMAT
    module procedure numNeedsFormat_double, numNeedsFormat_integer, numNeedsFormat_single
    module procedure numNeedsFormat_complex, numNeedsFormat_dcomplx
  end interface

  interface NUMTOCHARS
    module procedure numtochars_double, numtochars_integer, numtochars_single
  end interface

  ! Embeddded <cr> print multiple lines
  interface OUTPUT
    module procedure output_char, output_char_array, output_complex
    module procedure output_dcomplex, output_double
    module procedure output_integer, output_integer_array
    module procedure output_logical, output_logical_array
    module procedure output_single, output_double_array, output_single_array
    module procedure output_string
  end interface

  ! Don't filter for <cr>
  interface OUTPUT_
    module procedure OUTPUT_CHAR_NOCR
  end interface

  interface OUTPUTLIST
    module procedure OUTPUTLIST_INTS, OUTPUTLIST_CHARS
  end interface

  interface outputNamedValue
    module procedure output_nvp_character
    module procedure output_nvp_complex
    module procedure output_nvp_dbl_array, output_nvp_double
    module procedure output_nvp_int_array, output_nvp_integer
    module procedure output_nvp_log_array, output_nvp_logical
    module procedure output_nvp_sngl_array, output_nvp_single
  end interface

  interface TAB
    module procedure blanksToTab
  end interface
  
  interface TIMESTAMP
    module procedure timestamp_char, timestamp_integer, timestamp_logical
  end interface
  
  ! We can use the OutputLines mechanism for user-controlled
  ! buffering, filtering, grep-ing, or whatever
  integer, parameter :: MAXOUTPUTLINESLEN = 2048 ! How many chars it can hold
  character(len=MAXOUTPUTLINESLEN), public, save     :: outputLines = ' '

  ! This is the type for configuring how to automatically format
  ! lines and whether they should be sent to stdout or elsewhere
  type outputOptions_T
    integer :: PRUNIT = STDOUTPRUNIT    ! Unit for output (see comments above).  
    integer :: MLSMSG_Level        = MLSMSG_Info ! What level if logging
    integer :: newLineVal          = 10 ! 13 means <cr> becomes new line; -999 means ignore
    integer :: nArrayElmntsPerLine = 7
    integer :: nBlanksBtwnElmnts   = 3
    logical :: BUFFERED            = .true.
    logical :: OPENED              = .false.
    logical :: SKIPMLSMSGLOGGING   = .false.
    logical :: usePatternedBlanks  = .true. ! Use patterns for special fillChars
    character(len=9) :: specialFillChars = '123456789'
    character(len=9) :: lineupFillChars =  'ynnnnynnn' ! whether they line up
    character(len=16), dimension(9) :: patterns = (/ & ! on consecutive lines
      &  '(. )            ' , &
      &  '(. .)           ' , &
      &  '(.  .)          ' , &
      &  '(.   .)         ' , &
      &  '(.. ..)         ' , &
      &  '(- )            ' , &
      &  '(- -)           ' , &
      &  '(-  -)          ' , &
      &  '(- .. )         ' /)
      !   12345678901234567890
    character(len=FileNameLen) :: name = 'stdout'
    character(len=3)  :: advanceDefault = 'no' ! if advance=.. missing
    character(len=12) :: sdFormatDefault = '*' ! * means default format spec
    character(len=1)  :: arrayElmntSeparator = ' '
  end type
  
  type(outputOptions_T), public, save :: outputOptions

  ! This is the type for configuring whether and how to automatically stamp
  ! lines sent to stdout
  ! If interval > 1 only a fraction of lines will be stamped
  ! If interval > 10, the stamps will not be in-line, instead
  ! they will appear alone as page headers
  ! (As an alternative, use timeStamp to stamp only individual lines)
  type stampOptions_T
    logical :: neverStamp = .false.  ! if true, forget about automatic stamping
    logical :: post       = .true.      ! Put stamp at end of line?
    logical :: showTime   = .false. ! Don't show date or time unless TRUE
    character(len=24) :: textCode = ' '
    ! Don't show date unless dateFormat is non-blank
    character(len=16) :: dateFormat = ' '
    character(len=16) :: timeFormat = 'hh:mm'
    integer :: interval = 1 ! 1 means stamp every line; n means every nth line
    character(len=8) :: TIMESTAMPSTYLE = 'post' ! 'pre' or 'post'
  end type
  
  type(stampOptions_T), public, save :: stampOptions ! Could leave this private

  ! This is the type for configuring how timeStamp stamps only individual lines)
  type timeStampOptions_T
    logical :: post = .true.      ! Put stamp at end of line?
    logical :: showDate = .false. ! Don't show date unless TRUE
    character(len=24) :: textCode = ' '
    ! Don't show date unless dateFormat is non-blank
    character(len=16) :: dateFormat = 'yyyy-mm-dd'
    character(len=16) :: timeFormat = 'hh:mm:ss'
    character(len=8) :: TIMESTAMPSTYLE = 'post' ! 'pre' or 'post'
  end type
  
  type(timeStampOptions_T), public, save :: timeStampOptions ! Could leave this private

  ! Private parameters
  logical, save, private :: SILENTRUNNING = .false. ! Suspend all further output
  integer, save, private :: ATCOLUMNNUMBER = 1  ! Where we'll print next
  logical, save, private :: ATLINESTART = .true.  ! Whether to stamp if notpost
  integer, save, private :: LINESSINCELASTSTAMP = 0
  logical, private, parameter :: LOGEXTRABLANKS = .false.
  integer, private, parameter :: MAXNUMTABSTOPS = 24
  integer, private, parameter :: RECLMAX = 1024  ! This is NAG's limit
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
  real, parameter, dimension(3) :: RPREFERDEFAULTFORMAT = &
    & (/ -1., 0., 1. /)  ! For which values to use default format '*'
  character(len=16), private, save :: NONFINITEFORMAT = '(1pg14.6)' ! 'NaN, Inf'
  character(len=12), private :: sdNeedsFormat = '(1pg14.6)'
  character(len=12), private :: sdNeedsFragment = '(1pg14'

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

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
      endif
    endif
    if ( columnRange(1) > 0 ) then
      spaces = columnRange(2) - max( columnRange(1), atColumnNumber )
      if ( spaces < 1 ) return
      if ( columnRange(1) > atColumnNumber ) &
        & call blanks( columnRange(1) - atColumnNumber )
    else
      spaces = columnRange(2) - columnRange(1)
    endif
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
    endif
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
  ! For multiline messages, you may divide them into elements of
  ! a character array, or else a longer character scalar and
  ! supply LineLength for the routine to wrap at word boundaries
  subroutine BANNER_CHARS ( CHARS, COLUMNRANGE, ALIGNMENT, SKIPS, LINELENGTH )
    character(len=*), intent(in)                :: CHARS
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    integer, optional, intent(in)               :: LINELENGTH
    !
    ! Internal variables
    integer :: addedLines
    character(len=1)      :: myAlignment
    integer, dimension(2) :: myColumnRange
    integer :: lineLen, mySkips,  padding
    character(len=2*len(chars))      :: wrappedChars
    character(len=160), dimension(:), pointer :: lines
    logical, parameter :: DEBUG = .false.
    ! Executable
    myAlignment = 'C'
    if ( present(alignment) ) myAlignment = alignment
    mySkips = 0
    if ( present(skips) ) mySkips = skips
    if ( present(LineLength) ) then
      ! We will wrap the input to fit within LineLength, but remembering
      ! the stars and padding
      call wrap( chars, wrappedChars, width=LineLength-4, &
        & inseparator=achar(0), addedLines=addedLines )
      addedLines = addedLines + 1
      allocate(lines(addedLines))
      lines = ' '
      call List2Array( wrappedChars, lines, &
        & countEmpty=.true., inseparator=achar(0) )
      call banner( lines, alignment=alignment )
      deallocate(lines)
      return
    elseif ( present(columnRange) ) then
      myColumnRange = columnRange
    else
      lineLen = max( 80, 4 + len_trim(chars)*(1+mySkips) )
      padding = ( lineLen - len_trim(chars)*(1+mySkips) ) / 2
      myColumnRange(1) = 1 + padding
      myColumnRange(2) = lineLen - padding
    endif
    
    ! define padding as the larger of columnrange(1) and 1
    padding = max( 1, myColumnRange(1) )
    LineLen = padding + myColumnRange(2) - 1
    if ( DEBUG ) then
      call outputnamedValue( 'padding', padding )
      call outputnamedValue( 'LineLen', LineLen )
      call outputnamedValue( 'myColumnRange', myColumnRange )
    endif
    ! Top border
    call output( '*' )
    call blanks ( lineLen-2, FillChar='-' )
    call output( '*', advance = 'yes' )
    ! Left star, then message, then right star
    call output( '*' )
    call alignToFit( chars, myColumnRange, myAlignment, skips )
    call blanksToColumn( lineLen )
    call output( '*', advance = 'yes' )
    ! Bottom border
    call output( '*' )
    call blanks ( lineLen-2, FillChar='-' )
    call output( '*', advance = 'yes' )
  end subroutine BANNER_CHARS

  subroutine BANNER_CHARARRAY ( CHARARRAY, COLUMNRANGE, ALIGNMENT, SKIPS )
    character(len=*), dimension(:), intent(in)  :: CHARARRAY
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: SKIPS ! How many spaces between chars
    !
    ! Internal variables
    integer :: i
    ! Internal variables
    character(len=1)      :: myAlignment
    integer, dimension(2) :: myColumnRange
    integer :: lineLen, mySkips,  padding
    ! Executable
    myAlignment = 'C'
    if ( present(alignment) ) myAlignment = alignment
    mySkips = 0
    if ( present(skips) ) mySkips = skips
    if ( present(columnRange) ) then
      myColumnRange = columnRange
    else
      lineLen = 80
      padding = LineLen
      do i = 1, size(chararray)
        lineLen = max( lineLen, 4 + len_trim(chararray(i))*(1+mySkips) )
        padding = min( padding, &
          & ( lineLen - len_trim(chararray(i))*(1+mySkips) ) / 2 )
      enddo
      myColumnRange(1) = 1 + padding
      myColumnRange(2) = lineLen - padding
    endif
    
    ! define padding as the larger of columnrange(1) and 1
    padding = max( 1, myColumnRange(1) )
    LineLen = padding + myColumnRange(2) - 1
    ! Top border
    call output( '*' )
    call blanks ( lineLen-2, FillChar='-' )
    call output( '*', advance = 'yes' )
    do i = 1, size(chararray)
      ! Left star, then message, then right star
      call output( '*' )
      call alignToFit( chararray(i), myColumnRange, myAlignment, skips )
      call blanksToColumn( lineLen )
      call output( '*', advance = 'yes' )
    enddo
    ! Bottom border
    call output( '*' )
    call blanks ( lineLen-2, FillChar='-' )
    call output( '*', advance = 'yes' )
  end subroutine BANNER_CHARARRAY

  ! -----------------------------------------------------  BLANKS  -----
  subroutine BLANKS ( N_BLANKS, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    integer :: I    ! Blanks to write in next WRITE statement
    logical :: lineup
    integer :: ntimes
    integer :: numSoFar
    character(len=16) :: pattern
    integer :: patternLength
    integer :: patternNum
    integer :: theRest
    ! Executable
    if ( present(fillChar) ) then
      if ( outputOptions%usePatternedBlanks .and. &
        & index(outputOptions%specialFillChars, FILLCHAR) > 0 ) then
        ! We need to try to fit our called-for pattern into n_blanks
        ! The 1st question is, how many times could it be done?
        read(fillChar, *) patternNum
        pattern = outputOptions%patterns(patternNum)
        ! The pattern length (adjusted for enclosing parentheses)
        patternLength = len_trim(pattern) - 2
        ! Now we assume we'll always want the first and blanks of n_blanks to be
        ! purely blank
        if ( patternLength > n_blanks - 2 ) then
          ! n_blanks too short--just print blanks
          call pr_blanks ( n_blanks, advance=advance, dont_stamp=dont_stamp )
          return
        end if
        ntimes = (n_blanks-2)/patternLength
        ! In case we want latterns on consecutive lines to line up; viz
        ! a: . . . . . . Something
        ! xy:. . . . . . Something else
        lineup = ( outputOptions%lineupFillChars(patternNum:patternNum) == 'y' )
        if ( lineup ) then
          numSoFar = 0
          ! Make sure that we always begin on an even-numbered column
          ! (This only works for patterns like '. . . ' or '- - - '
          if ( mod(atColumnNumber, 2) /= 0 ) then
            call pr_blanks ( 1, advance='no' )
            numSoFar = 1
          end if
        else
          call pr_blanks ( 1, advance='no' )
          numSoFar = 1
        end if
        do i=1, ntimes
          call output_ ( pattern(2:patternLength+1), advance='no' )
          ! if ( xtraBlanks > 0 ) call pr_blanks ( xtraBlanks, advance='no' )
          numSoFar = numSoFar + patternLength
        enddo
        theRest = n_blanks - numSoFar
        if ( theRest > 0 ) call pr_blanks ( theRest, advance=advance, dont_stamp=dont_stamp )
        return
      end if
    end if
    call pr_blanks ( n_blanks, fillChar=fillChar, advance=advance, dont_stamp=dont_stamp )
  end subroutine BLANKS

  ! -----------------------------------------------------  BLANKSTOCOLUMN  -----
  subroutine BLANKSTOCOLUMN ( COLUMN, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT out to column COLUMN.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: COLUMN
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    ! Executable
    if ( ATCOLUMNNUMBER >= COLUMN ) return
    call blanks( COLUMN-ATCOLUMNNUMBER, fillChar, advance, dont_stamp )
  end subroutine BLANKSTOCOLUMN

  ! ------------------------------------------------  blanksToTab  -----
  ! Print blanks out to next tabstop
  ! (or else to tabstop number tabn)
  subroutine blanksToTab ( tabn, fillChar )
    ! Args
    integer, optional, intent(in) :: TABN
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    ! Internal variables
    integer :: nTab
    ! Executable
    if ( present(tabn) ) then
      if ( tabn < 1 .or. tabn > MAXNUMTABSTOPS ) return
      if ( atColumnNumber < tabStops(tabn) ) &
        & call blanksToColumn( tabStops(tabn), fillChar )
    else
      nTab = findFirst( tabStops > atColumnNumber )
      if ( nTab > 0 ) &
        & call blanksToColumn( tabStops(nTab), fillChar )
    endif
  end subroutine blanksToTab

  ! ---------------------------------------------- DumpOuputOptions -----
  subroutine DumpOutputOptions(options)
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
      call outputNamedValue ( 'meaning', prunitname(options%prUnit), &
        & advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    endif
    call outputNamedValue ( 'file name', trim(options%name), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'logging level', options%MLSMSG_Level, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'buffered?', options%buffered, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'skip MLSMsg logging?', options%SKIPMLSMSGLOGGING, advance='yes', &
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
    enddo
    call newline
    if ( CHECKINGTABBING ) then
      do i=1, MAXNUMTABSTOPS
        call blanksToColumn( tabStops(i), fillChar=fillChar )
        call output_( '^', advance='no' )
      enddo
      call newline
    endif
    do
      call output_( decade, advance='no' )
      if ( atColumnNumber > 132 ) exit
    enddo
    call newline
    call blanks(80, fillChar='-', advance='yes')
  contains
    function PRUnitName( prUnit ) result( name )
      ! Return an appropriate name for the prUnit number
      ! Args
      integer, intent(in) :: prUnit
      character(len=12) :: name
      ! Executable
      select case ( prUnit )
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
    end function PRUnitName
  end subroutine DumpOutputOptions

  ! ---------------------------------------------- DumpStampOptions -----
  subroutine DumpStampOptions(options)
    ! Show output options
    type(StampOptions_T), intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     call blanks(80, fillChar='-', advance='yes')
    call headline( 'Summary of automatic stamp options', &
      & fillChar='-', before='*', after='*' )
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
  subroutine DUMPTIMESTAMPOPTIONS(options)
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
  subroutine DumpSize_double ( n, advance, units, Before, After )
    double precision, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    real, intent(in), optional :: units
    character(len=*), intent(in), optional :: Before, After
    ! Local parameters
    real, parameter :: KB = 1024.0
    real, parameter :: MB = KB * 1024.0
    real, parameter :: GB = MB * 1024.0
    real, parameter :: TB = GB * 1024.0
    double precision :: Amount ! N * MyUnits
    character(len=8) :: HowMuch
    real             :: myUnits
    character(len=6) :: Suffix
    ! Make a 'nice' output
    if ( present(before) ) call output_ ( before )
    myUnits = 1.0
    if ( present(units) ) myUnits = units
    if ( myUnits == 0.0 ) then
      call output ( n, format='(e12.1)' )
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
    elseif ( amount > 999999 ) then ! I6 format limits this
      call output_( '(HUGE)' )
    elseif ( amount == int(amount) ) then
      write ( howMuch, '(i6)' ) int(amount)
    else
      write ( howMuch, '(f6.1)' ) amount
    end if
    call output_ ( trim(adjustl(howMuch)) )
    call output_ ( trim(suffix) )
    if ( present(after) ) call output_ ( after )
    call output_ ( '', advance=advance )
  end subroutine DumpSize_double

  ! --------------------------------------------- DumpSize_integer -----
  subroutine DumpSize_integer ( n, advance, units, Before, After )
    integer, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    integer, intent(in), optional :: units ! E.g., 1024 for kB
    character(len=*), intent(in), optional :: Before, After
    ! Executable
    if ( present(units) ) then
      call dumpSize ( dble(n), advance=advance, units=real(units), &
        & before=before, after=after )
    else
      call dumpSize ( dble(n), advance=advance, before=before, after=after )
    end if
  end subroutine DumpSize_integer

  ! ------------------------------------------------ DumpSize_real -----
  subroutine DumpSize_real ( n, advance, units, Before, After )
    real, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    real, intent(in), optional :: units
    character(len=*), intent(in), optional :: Before, After
    ! Make a 'nice' output
    call dumpsize ( dble(n), advance, units, before, after )
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
    endif
  end subroutine dumpTabs

  ! ----------------------------------------------  flushOutputLines  -----
  ! print or log OutputLines
  ! then reset to ''
  subroutine flushOutputLines ( prUnit )
    ! Args
    integer, optional, intent(in) :: prUnit ! How do you want 'em?
    ! Local arguments
    integer :: kNull  ! 1 past where the end of str occurs
    integer :: myPrUnit
    integer :: oldPrUnit
    ! Executable
    myPrUnit = -1 ! By default, just print to stdout
    if ( present(prUnit) ) myprUnit = prUnit
    oldprUnit = outputOptions%prUnit
    outputOptions%prUnit = myprUnit
    kNull = index( OutputLines, achar(0), back=.true. )
    if ( kNull < 2 ) then
      ! OutputLines hasn't been defined yet
      return
    elseif ( kNull > len(OutputLines) ) then
      ! Something very wrong
      return
    else
      call output( OutputLines(1:kNull-1) )
    endif
    outputOptions%prUnit = oldPrUnit
    OutputLines = ' '
  end subroutine flushOutputLines

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
    & COLUMNRANGE, ALIGNMENT, SKIPS )
    character(len=*), intent(in)                :: CHARS
    character(len=1), intent(in), optional :: fillChar      ! For padding
    character(len=*), intent(in), optional :: Before, After ! text to print
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
    ! Executable
    if ( present(columnRange) ) then
      myFullColumnRange = columnRange
    else
      myFullColumnRange(1) = 1
      myFullColumnRange(2) = 80
    endif
    myColumnRange = myFullColumnRange
    mySkips = 0
    if ( present(skips) ) mySkips = skips
    myFillChar = ' '
    if ( present(fillChar) ) myFillChar = fillChar
    myAlignment = 'C'
    if ( present(alignment) ) myAlignment = alignment
    ! call outputNamedValue ( 'myFullColumnRange', myFullColumnRange )
    ! call outputNamedValue ( 'myColumnRange', myColumnRange )
    ! call outputNamedValue ( 'mySkips', mySkips )
    ! call outputNamedValue ( 'myFillChar', myFillChar )
    ! Adjust for lengths of before, after
    if ( present(before) ) myColumnRange(1) = myFullColumnRange(1) + len(before)
    if ( present(after) ) then
      myColumnRange(2) = myFullColumnRange(2) - len(after)
      rightpadding = len(after) - 1
    endif
    call blanksToColumn( myFullColumnRange(1), advance='no' )
    if ( present(before) ) call output(before, advance='no' )
    if ( mySkips == 0 .and. myAlignment == 'C' .and. myFillChar /= ' ' ) then
      ! OK, final adjustments of myColumnRange
      myColumnRange(1) = &
        & ( myColumnRange(1) + myColumnRange(2) + 1 - len(chars) )/2
      myColumnRange(2) =  myColumnRange(1) - 1 + len(chars)
      ! call outputNamedValue ( 'myColumnRange', myColumnRange )
      call blanksToColumn( myColumnRange(1), fillChar=fillChar, advance='no' )
      call aligntofit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn( myFullColumnRange(2)-rightpadding, &
        & fillChar=fillChar, advance='no' )
      if ( present(after) ) call output(after, advance='no' )
    else
      call aligntofit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn(myFullColumnRange(2)-rightpadding, advance='no' )
      if ( present(after) ) call output(after, advance='no' )
    endif
    call newLine
  end subroutine HEADLINE

  ! ----------------------------------------------  isOutputSuspended  -----
  logical function isOutputSuspended ()
  ! Have we suspended outputting to PRUNIT?
    isOutputSuspended = silentRunning
  end function isOutputSuspended

  ! ----------------------------------------------------  NewLine  -----
  subroutine NewLine
    call output_ ( '', advance='yes' )
  end subroutine NewLine

  ! ----------------------------------------------------  NextColumn  -----
  function NextColumn() result(Column)
    ! Args
    integer :: Column
    Column = atColumnNumber
  end function NextColumn

  ! ----------------------------------------------------  NextTab  -----
  function NextTab() result(Column)
    ! Args
    integer :: Column
    ! Internal variables
    integer :: nTab
    ! Executable
    Column = 0
    nTab = findFirst( tabStops > atColumnNumber )
    if ( nTab > 0 ) Column = max( tabStops(nTab), atColumnNumber )
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
    endif
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
    endif
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
    endif
    include 'numToChars.f9h'
  end function numToChars_single

  ! ---------------------------------------  OUTPUTCALENDAR  -----
  subroutine OUTPUTCALENDAR ( date, datenote, notes, dontWrap )
    ! output a nicely-formatted calendar of current month with
    ! today's date in "bold"
    ! Args
    character(len=*), intent(in), optional :: date ! date instead of current one
    ! dateNote, (notes), if present, is (an array of)
    ! stringLists, (one per day in the month,)
    character(len=*), optional :: dateNote ! Note for the current date
    ! Each string list contains either a blank for a date, meaning
    ! nothing will be printed in the calendar square for that date,
    ! or else it contains '/'-separated lines of text, each of
    ! which will be printed on a separate line within the square
    character(len=*), dimension(:), optional :: notes
    logical, optional                        :: dontWrap ! Dont wrap notes to fit
    ! Internal variables
    integer, parameter :: MAXNOTELENGTH = 256
    ! This should be modified for internationalization; e.g. with
    ! an include statement or suchlike
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
    integer :: day
    integer, dimension(6,7) :: days, daysOfYear
    integer :: ErrTyp
    integer :: iwk
    integer :: month
    logical :: myDontWrap
    character(len=10) :: noteString
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
    endif
    col1 = index(lowercase(dateString), 't')
    if ( col1 > 0 ) then
      date2 = dateString(1:col1-1)
       if ( nCopies(dateString(:col1-1), '-') < 2 ) &
        & date2 = reformatDate( dateString(1:col1-1), fromForm='yyyy-Doy', toForm=utcformat )
    else
      date2 = reformatDate( dateString, fromForm='*', toForm=utcformat )
    endif
    call utc_to_yyyymmdd( date2, ErrTyp, year, month, day )
    if ( month < 0 ) then
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
    enddo
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
      endif
      numRows = max( numRows, &
        & NumStringElements( wrappedNote, countEmpty, inseparator ) + 2 &
        & )
    endif
    if ( present(notes) ) then
      do aday=1, min( size(notes), daysInMonth( month, year ) )
        if ( myDontWrap ) then
          wrappedNote = notes(aday)
        else
          call wrap( notes(aday), wrappedNote, 10, '/' )
        endif
        numRows = max( numRows, &
          & NumStringElements( wrappedNote, countEmpty, inseparator ) + 2 &
          & )
      enddo
    endif
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
          endif
          if ( days(iwk, wkdy) < 1 ) then
            ! Don't write notes or anything else in "empty" days
          elseif ( row == 1 ) then
            call writeIntsToChars( days(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call alignToFit( trim(dateString), (/ col1, col2-1 /), 'r' )
          elseif( row == numRows ) then
            call writeIntsToChars( daysOfYear(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call alignToFit( 'd' // trim(dateString), (/ col1, col2-1 /), 'r' )
          elseif( present(dateNote) .and. today ) then
            if ( myDontWrap ) then
              wrappedNote = dateNote
            else
              call wrap( dateNote, wrappedNote, 10, '/' )
            endif
            call GetStringElement ( wrappedNote, noteString, &
              & row-1, countEmpty, inseparator )
            if ( noteString == inseparator ) noteString = ' '
            call output_( noteString )
          elseif( present(notes) ) then
            if ( days(iwk, wkdy) <= size(notes) ) then
              if ( myDontWrap ) then
                wrappedNote = notes(days(iwk, wkdy))
              else
                call wrap( notes(days(iwk, wkdy)), wrappedNote, 10, '/' )
              endif
              call GetStringElement ( wrappedNote, noteString, &
                & row-1, countEmpty, inseparator )
              if ( noteString == inseparator ) noteString = ' '
              call output_( noteString )
            endif
          endif
          if ( today ) then
            call blanksToColumn(col2-1)
            call output_('|')
          else
            call blanksToTab
          endif
        enddo ! wkdy
        call output_('|')
        call newline
      enddo ! row
      ! begin with 
    enddo ! week
    call blanksToTab( 7, fillChar='-' )
    call newline
    ! Restore tabstops
    call settabs( '5-120+5' )
  end subroutine OUTPUTCALENDAR

  ! ------------------------------------------------  OUTPUT_CHAR  -----
  ! Output CHARS to PRUNIT.
  subroutine OUTPUT_CHAR ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, DONT_STAMP, &
    & NEWLINEVAL, DONT_ASCIIFY )
    ! We will 1st check to see whether any internal characters are
    ! codes for newlines
    ! If any are, we will call newLine in place of printing
    ! them
    ! (This is a new default behavior; you can restore
    ! the old by passing an impossible value for NewLineVal, e.g. -999)
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    integer, intent(in), optional :: NEWLINEVAL ! What char val to treat as <cr>
    logical, intent(in), optional          :: DONT_ASCIIFY ! output binary
    ! Internal variables
    integer :: I ! loop inductor
    integer :: myNewLineVal
    logical :: myAsciify
    ! Executable
    myNewLineVal = outputOptions%newlineVal
    if ( present(newLineVal) ) myNewLineVal = newLineVal
    myAsciify = .true.
    if ( present(dont_asciify) ) myAsciify = .not. dont_asciify
    i = index( chars, achar(myNewLineVal) )
    if ( i < 1 ) then
      if ( myAsciify) then
        call OUTPUT_CHAR_NOCR ( ReplaceNonAscii(CHARS, '@', exceptions=achar(9)), &
          & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, DONT_STAMP )
      else
        call OUTPUT_CHAR_NOCR ( CHARS, &
          & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, DONT_STAMP )
      endif
    else
      do i=1, len(chars)
        if ( chars(i:i) /= achar(myNewLineVal) ) then
          if ( myAsciify ) then
            call OUTPUT_CHAR_NOCR ( ReplaceNonAscii(CHARS(i:i), '@', exceptions=achar(9)), &
              & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
              & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK, &
              & DONT_STAMP=DONT_STAMP )
          else
            call OUTPUT_CHAR_NOCR ( CHARS(i:i), &
              & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
              & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK, &
              & DONT_STAMP=DONT_STAMP )
          endif
        else
          call newLine
        endif
      enddo
      if ( Advance_is_yes_or_no(advance) == 'yes' ) call newLine
    endif
  end subroutine OUTPUT_CHAR

  subroutine OUTPUT_CHAR_NOCR ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, DONT_STAMP )
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    !
    integer :: i1, i2
    integer :: IOBloc
    integer :: nIOBlocs
    character(len=max(16,len(chars)+1)) :: my_chars
    character(len=len(chars)+64) :: stamped_chars ! What to print to stdout
    character(len=max(16,len(chars)+1)) :: the_chars
    logical :: my_dont_log
    logical :: my_dont_stamp
    character(len=3) :: MY_ADV
    integer :: n_chars
    integer :: n_stamp ! How much of stamped_chars to print
    logical :: stamp_header
    logical :: stamped
    integer :: status
    !
    if ( SILENTRUNNING ) return
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_stamp = stampOptions%neverStamp ! .false.
    if ( present(dont_stamp) ) my_dont_stamp = dont_stamp
    my_dont_stamp = ( my_dont_stamp .or. &
      & linesSinceLastStamp < (stampOptions%interval - 1) )
    stamped = .false.
    stamp_header = .false.
    stamped_chars = chars
    ! Are we trying to avoid buffered output?
    ! If so, note the following decision and its effects:
    ! we will close and reopen for append our output file
    ! not at the beginning of a write to a fresh line
    ! but at the end of a write to the previous line.

    ! If, in the middle, we switch output files, say from file1 to file2,
    ! it's possible that a line that should have been the last written
    ! to file1 will instead be written to file2. 
    ! That's the downside.

    ! The good part is that if we crash hard we won't lose that
    ! last line.
    if ( (.not. outputOptions%buffered) .and. &
      & outputOptions%prUnit > 0 .and. &
      & (.not. outputOptions%opened) .and. &
      & ( atcolumnnumber == 1 ) ) then
      open (outputOptions%prUnit, file=trim(outputOptions%name), status='replace', &
        & form='formatted', access='sequential', iostat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to reopen prUnit for replace ' // trim(outputOptions%name) )
      outputOptions%opened = .true.
    endif
    ! Do we need to stamp this line? If so, at beginning or at end?
    if ( my_dont_stamp ) then
    elseif ( stampoptions%interval > 9 ) then
      stamp_header = (my_adv == 'yes')
    elseif ( stampoptions%post ) then
      if ( my_adv == 'yes' ) then
        stamped_chars = stamp(chars)
        stamped = .true.
      end if
    elseif( ATLINESTART ) then
      stamped_chars = stamp(chars)
      stamped = .true.
    end if

    n_chars = len(chars)
    if (LOGEXTRABLANKS) n_chars = max(len(chars), 1)
    my_dont_log = outputOptions%skipmlsmsglogging ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    n_stamp = len_trim(stamped_chars)
    if ( my_adv == 'no' ) n_stamp = n_stamp + len(chars) - len_trim(chars)
    ! Special case: if chars is blank (chars are blank?)
    ! we'll want to print anyway
    if ( len_trim(chars) < 1 ) n_stamp = max(n_stamp, 1)
    if ( any(outputOptions%prunit == (/STDOUTPRUNIT, BOTHPRUNIT/)) .and. &
      &  n_stamp > RECLMAX ) then
      nIOBlocs = 1 + (n_stamp-1)/RECLMAX
      i2 = 0
      do IOBloc=1, nIOBlocs
        i1 = i2 + 1
        i2 = min(i2+RECLMAX, n_stamp)
        write ( *, '(a)', advance='no' ) stamped_chars(i1:i2)
      enddo
      if ( my_adv == 'yes' ) write ( *, '(a)', advance=my_adv ) ' '
    elseif ( any(outputOptions%prunit == (/STDOUTPRUNIT, BOTHPRUNIT/)) .and. &
      &  len(chars) < 1 .and. my_adv == 'yes' ) then
      write ( *, '(a)', advance=my_adv )
    elseif ( any(outputOptions%prunit == (/STDOUTPRUNIT, BOTHPRUNIT/)) .and. &
      &  n_stamp > 0 ) then
      write ( *, '(a)', advance=my_adv ) stamped_chars(1:n_stamp)
    endif
    if ( any(outputOptions%prunit == (/MSGLOGPRUNIT, BOTHPRUNIT/)) .and. .not. my_dont_log  ) then
      ! We must use MLSMessage to log the chars
      the_chars = chars // ' '
      if (LOGEXTRABLANKS) n_chars = max(len(chars), 1)
      if ( present(log_chars) ) then
        if ( len_trim(log_chars) > 0 ) then
          n_chars = len_trim(log_chars) + 1
          the_chars = log_chars(:n_chars-1) // ' '
        end if
      end if
      if ( the_chars == ' ' .and. present(insteadofblank)  ) then
        my_chars = trim(insteadofblank) // ' '
        if (LOGEXTRABLANKS) n_chars = max(len(insteadofblank), 1)
      else
        my_chars = the_chars
      end if
      if ( my_adv == 'no' .and. LOGEXTRABLANKS ) n_chars = n_chars+1
      n_chars = min(n_chars, len(my_chars))
      if ( my_adv == 'yes' ) n_chars = max(n_chars, 1)
      if ( n_chars < 1 ) then
      elseif ( present(from_where)  ) then
        call MLSMessage ( outputOptions%MLSMSG_Level, from_where, my_chars(1:n_chars), &
          & advance=my_adv )
      else
        call MLSMessage ( outputOptions%MLSMSG_Level, ModuleName, my_chars(1:n_chars), &
          & advance=my_adv )
      end if
    elseif ( outputOptions%prunit ==  OUTPUTLINESPRUNIT ) then
      ! Append to outputLines; maybe print later on
      !   if ( len_trim(chars) > 0 ) then
      !     if ( len_trim(outputLines) > 0 ) then
      !       outputLines = trim(outputLines) // trim(chars)
      !     else
      !       outputLines = chars
      !     endif
      call append_chars( outputLines, chars )
        if ( my_adv == 'yes' ) &
          call append_chars( outputLines, achar(outputOptions%NewLineVal) )
      !    & outputLines = trim(outputLines) // achar(outputOptions%NewLineVal) ! add <cr>
      ! endif
    end if
    
    if ( outputOptions%prunit < 0  ) then
      ! Already logged; no output to stdout
    else if ( stamped_chars == ' ' .and. present(insteadofblank)  ) then
      write ( outputOptions%prunit, '(a)', advance=my_adv ) trim_safe(insteadofblank)
    else if ( stamped_chars == ' ' .and. my_adv == 'no' ) then
      write ( outputOptions%prunit, '(a)', advance='no' ) ' '
    elseif ( len_trim(chars) < 1 .and. n_stamp == 1 ) then
      write ( outputOptions%prunit, '(a)', advance=my_adv )
    else
      write ( outputOptions%prunit, '(a)', advance=my_adv ) stamped_chars(1:n_stamp)
    end if
    atLineStart = (my_adv == 'yes')
    if ( atLineStart ) then
      ! Are we trying to avoid buffered output?
      if ( (.not. outputOptions%buffered) .and. &
        & outputOptions%prUnit > 0 .and. &
        & outputOptions%opened ) then
        close(outputOptions%prUnit)
        open (outputOptions%prUnit, file=trim(outputOptions%name), status='old', &
          & position='append', &
          & form='formatted', access='sequential', iostat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to reopen prUnit for append ' // trim(outputOptions%name) )
      endif
      if ( stamp_header ) then
        ! Time to add stamp as a page header on its own line
        stamped_chars = stamp(' ')
        stamped = .true.
        if ( any(outputOptions%prunit == (/STDOUTPRUNIT, BOTHPRUNIT/)) ) then
          write ( *, '(a)', advance='yes' ) trim(stamped_chars)
        elseif( outputOptions%prunit > 0 ) then
          write ( outputOptions%prunit, '(a)', advance='yes' ) trim(stamped_chars)
        end if
      end if
      if ( stamped ) then
        linesSinceLastStamp = 0
      else
        linesSinceLastStamp = linesSinceLastStamp + 1
      end if
    end if
    atColumnNumber = atColumnNumber + n_chars
    if ( atLineStart ) atColumnNumber = 1
  contains
    subroutine append_chars ( str, chars )
      ! Append chars to end of str
      ! where end is marked by the last null character (achar(0))
      ! Args
      character(len=*), intent(inout) :: str
      character(len=*), intent(in)    :: chars
      ! Local variables
      integer :: kNull  ! 1 past where the end of str occurs
      ! Executable
      kNull = index( str, achar(0), back=.true. ) + 1
      if ( kNull < 2 ) then
        ! str hasn't been defined yet
        str = chars // achar(0)
      elseif ( kNull > len(str) + 1) then
        ! Something very wrong
        str = achar(0)
      else
        str(kNull-1:) = chars // achar(0)
      endif
    end subroutine append_chars
  end subroutine OUTPUT_CHAR_NOCR

  ! ------------------------------------------  OUTPUT_CHAR_ARRAY  -----
  subroutine OUTPUT_CHAR_ARRAY ( CHARS, ADVANCE, INSTEADOFBLANK, NEWLINEVAL )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    integer, intent(in), optional :: NEWLINEVAL ! What char val to treat as <cr>
    ! Internal variables
    integer :: I ! loop inductor
    ! Executable
    do i = 1, size(chars)
      call output ( chars(i), &
        & insteadofblank=insteadofblank, newLineVal=newLineVal )
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance )
    end if
  end subroutine OUTPUT_CHAR_ARRAY

  ! ---------------------------------------------  OUTPUT_COMPLEX  -----
  subroutine OUTPUT_COMPLEX ( VALUE, Format, ADVANCE, Before, After, dont_stamp )
    complex, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: DONT_STAMP
    character(len=60) :: LINE

    if ( present(Format)  ) then
      write ( line, format ) value
    else
      write ( line, '("(",1pg15.7,",",1pg15.7,")")' ) value
    end if
    if ( present(before) ) call output_ ( before, dont_log = .true. )
    if ( present(after)  ) then
      call output_ ( trim(line), dont_log = .true. )
      call output_ ( after, advance=advance, dont_log = .true., &
        & dont_stamp=dont_stamp )
    else
      call output_ ( trim(line), advance=advance, dont_log = .true., &
        & dont_stamp=dont_stamp )
    end if
  end subroutine OUTPUT_COMPLEX

  ! ---------------------------------------  OUTPUT_DATE_AND_TIME  -----
  subroutine OUTPUT_DATE_AND_TIME ( date, time, &
    & from_where, msg, dateFormat, timeFormat, CPU_Seconds, advance )
    ! Output nicely-formatted date, time, and extra message
    ! We'll assume we won't want this line stamped with date and time
    ! (for fear of being redundant, which we fear)
    logical, intent(in), optional :: date ! output date as character string
    logical, intent(in), optional :: time ! output time as character string
    character(len=*), intent(in), optional :: FROM_WHERE
    character(len=*), intent(in), optional :: MSG
    character(len=*), intent(in), optional :: DATEFORMAT
    character(len=*), intent(in), optional :: TIMEFORMAT
    double precision, intent(in), optional :: CPU_Seconds
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
      & advance=merge ( 'no ', my_adv, present(CPU_seconds) ), &
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
        & advance=my_adv, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_DATE_AND_TIME

  ! --------------------------------------------  OUTPUT_DCOMPLEX  -----
  subroutine OUTPUT_DCOMPLEX ( VALUE, Format, ADVANCE, Before, After )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=60) :: LINE

    if ( present(Format)  ) then
      write ( line, format ) value
    else
      write ( line, '("(",1pg22.14,",",1pg22.14,")")' ) value
    end if
    if ( present(before) ) call output_ ( before, dont_log = .true. )
    if ( present(after)  ) then
      call output_ ( trim(line), dont_log = .true. )
      call output_ ( after, advance=advance, dont_log = .true. )
    else
      call output_ ( trim(line), advance=advance, dont_log = .true. )
    end if
  end subroutine OUTPUT_DCOMPLEX

  ! ----------------------------------------------  OUTPUT_DOUBLE  -----
  subroutine OUTPUT_DOUBLE ( VALUE, Format, LogFormat, ADVANCE, &
    & Before, After, dont_stamp )
  ! Output "double" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: LogFormat ! How to post to Log
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: DONT_STAMP
    integer :: I, J, K
    character(len=30) :: LINE, LOG_CHARS, FormatSpec

    FormatSpec = outputOptions%sdFormatDefault
    if ( any( value == DPREFERDEFAULTFORMAT ) )  FormatSpec = '*'
    if ( .not. is_what_ieee(finite_signal, value) ) FormatSpec = NONFINITEFORMAT
    if ( present(Format)  ) FormatSpec = Format
    include 'numToChars.f9h'
    log_chars = ' ' ! line
    if ( present(LogFormat)  ) then
      write ( log_chars, LogFormat ) value
    end if
    if ( present(before) ) call output_ ( before )
    if ( present(after)  ) then
      call output_ ( line(:k), log_chars=log_chars )
      call output_ ( after, advance=advance, dont_stamp=dont_stamp )
    else
      call output_ ( line(:k), advance=advance, log_chars=log_chars, &
        & dont_stamp=dont_stamp )
    end if

  end subroutine OUTPUT_DOUBLE

  ! ----------------------------------------  OUTPUT_DOUBLE_ARRAY  -----
  subroutine OUTPUT_DOUBLE_ARRAY ( values, FORMAT, LogFormat, ADVANCE, DONT_STAMP )
  ! Output double-precision values to PRUNIT.
    double precision, intent(in) :: values(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=3) :: MY_ADV
    integer :: I ! loop inductor
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      if ( i < size(values) ) then
        if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) &
          & call output_( outputOptions%arrayElmntSeparator, advance='no' )
        if ( mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
          call output_ ( '', advance='yes', DONT_STAMP=.true. )
        else
          call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
        endif
      endif
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_DOUBLE_ARRAY

  ! ---------------------------------------------  OUTPUT_INTEGER  -----
  subroutine OUTPUT_INTEGER ( INT, PLACES, ADVANCE, FILL, FORMAT, &
    & Before, After, DONT_STAMP )
  ! Output INT to PRUNIT using at most PLACES (default zero) places
  ! If 'fill' is present and true, fill leading blanks with zeroes (only
  ! makes sense if 'places' is specified).

  ! Note that we never have places and format specified simultaneously
  ! we therefore institute the following trick:
  ! if format is present and contains the string 'PLACES=', the integer
  ! following the '=' will be taken to be places
    integer, intent(in) :: INT
    integer, intent(in), optional :: PLACES
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: FILL
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: DONT_STAMP
    !
    logical :: My_Fill
    integer :: I, J
    character(len=12) :: LINE
    integer :: MY_PLACES
    logical :: formatEncodesPlaces ! Have we used format to encode places?
    my_places = 0
    if ( present(places)  ) then; my_places = places; end if
    my_fill = .false.
    if ( present(places) .and. present(fill) ) my_fill = fill
    formatEncodesPlaces = .false.
    if ( present(format) ) then
      if ( index( lowercase(format), 'places=' ) > 0 ) then
        formatEncodesPlaces = .true.
        call readIntsFromChars ( format, my_places, ignore='*=' )
      endif
    endif
    if ( present(format) .and. .not. formatEncodesPlaces ) then
      line = ' '
      write ( line, format ) int
      i = 1
      j = len_trim(line)
    else if ( my_fill  ) then
      write ( line, '(i6.6)' ) int
      i = 1
      j = 6
    else
      write ( line, '(i12)' ) int
      i = max( 1, min(len(line)+1-my_places, index(line,' ',back=.true.)+1) )
      j = len(line)
    end if
    if ( present(before) ) call output_ ( before )
    if ( present(after)  ) then
      call output_ ( line(i:j), DONT_STAMP=DONT_STAMP )
      call output_ ( after, advance=advance, DONT_STAMP=DONT_STAMP )
    else
      call output_ ( line(i:j), advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_INTEGER

  ! ---------------------------------------  OUTPUT_INTEGER_ARRAY  -----
  subroutine OUTPUT_INTEGER_ARRAY ( INTEGERS, ADVANCE, FORMAT, DONT_STAMP )
  ! Output INTEGERS to PRUNIT.
    integer, intent(in) :: INTEGERS(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    logical, optional, intent(in) :: DONT_STAMP
    character(len=3) :: MY_ADV
    integer :: I ! loop inductor
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(integers)
      call output ( integers(i), advance='no', format=format )
      if ( i < size(integers) ) then
        if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) &
          & call output_( outputOptions%arrayElmntSeparator, advance='no' )
        if ( mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
          call output_ ( '', advance='yes', DONT_STAMP=.true. )
        else
          call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
        endif
      endif
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_INTEGER_ARRAY

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL ( LOG, ADVANCE, BEFORE, DONT_STAMP )
  ! Output LOG to PRUNIT using at most PLACES (default zero) places
    logical, intent(in) :: LOG
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    logical, optional, intent(in) :: DONT_STAMP
    character(len=2) :: LINE
    if ( log ) then
      line=' T'
    else
      line=' F'
    end if
    if ( present(before) ) call output_ ( before, DONT_STAMP=DONT_STAMP )
    call output_ ( line, advance=advance, DONT_STAMP=DONT_STAMP )
  end subroutine OUTPUT_LOGICAL

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL_ARRAY ( logs, ADVANCE, BEFORE, DONT_STAMP, ONLYIF )
  ! Output LOG to PRUNIT using at most PLACES (default zero) places
    logical, dimension(:), intent(in) :: logs
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    logical, optional, intent(in) :: DONT_STAMP
    logical, optional, intent(in) :: ONLYIF ! Print only if true (false)
    ! Internal variables
    character(len=1), dimension(size(logs)) :: clogs
    character(len=size(logs)) :: logChars
    integer :: I ! loop inductor
    character(len=1) :: ifonlyWhat
    ! Executable
    if ( present(before) ) call output_ ( before, DONT_STAMP=DONT_STAMP )
    if ( present(onlyif) ) then
      if ( onlyif ) then
        ifonlyWhat = 'T'
      else
        ifonlyWhat = 'F'
      endif
      clogs = ' '
      where (logs .eqv. onlyif)
        clogs = ifonlyWhat
      end where
      logChars = transfer( clogs, logChars )
      call output( logChars, advance='no' )
    else
      do i = 1, size(logs)
        call output ( logs(i), advance='no' )
      if ( i < size(logs) ) then
        if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) &
          & call output_( outputOptions%arrayElmntSeparator, advance='no' )
        if ( mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
          call output_ ( '', advance='yes', DONT_STAMP=.true. )
        else
          call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
        endif
      endif
      end do
    endif
    if ( present(advance) ) call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
  end subroutine OUTPUT_LOGICAL_ARRAY

  ! ----------------------------------------------  OUTPUT_SINGLE  -----
  subroutine OUTPUT_SINGLE ( VALUE, FORMAT, LogFormat, ADVANCE, &
    & Before, After, DONT_STAMP )
  ! Output "SINGLE" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format  ! How to print
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, optional, intent(in) :: DONT_STAMP
    integer :: I, J, K
    character(len=30) :: LINE, LOG_CHARS, FormatSpec

    FormatSpec = outputOptions%sdFormatDefault
    if ( any( value == RPREFERDEFAULTFORMAT ) )  FormatSpec = '*'
    if ( .not. is_what_ieee(finite_signal, value) ) FormatSpec = NONFINITEFORMAT
    if ( present(Format)  ) FormatSpec = Format
    include 'numToChars.f9h'
    log_chars = ' ' ! line
    if ( present(LogFormat)  ) then
      write ( log_chars, LogFormat ) value
    end if
    if ( present(before) ) call output_ ( before, DONT_STAMP=DONT_STAMP )
    if ( present(after)  ) then
      call output_ ( line(:k), log_chars=log_chars, DONT_STAMP=DONT_STAMP )
      call output_ ( after, advance=advance, DONT_STAMP=DONT_STAMP )
    else
      call output_ ( line(:k), advance=advance, log_chars=log_chars, &
        & DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_SINGLE

  ! ----------------------------------------  OUTPUT_SINGLE_ARRAY  -----
  subroutine OUTPUT_SINGLE_ARRAY ( values, FORMAT, LogFormat, ADVANCE, DONT_STAMP )
  ! Output single-precision values to PRUNIT.
    real, intent(in) :: values(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    integer :: I ! loop inductor
    character(len=3) :: MY_ADV
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      if ( i < size(values) ) then
        if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) &
          & call output_( outputOptions%arrayElmntSeparator, advance='no' )
        if ( mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
          call output_ ( '', advance='yes', DONT_STAMP=.true. )
        else
          call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
        endif
      endif
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_SINGLE_ARRAY

  ! ----------------------------------------------  OUTPUT_STRING  -----
  subroutine OUTPUT_STRING ( STRING, LENSTRING, ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS )
  ! Output STRING to PRUNIT.
    character(len=*), intent(in) :: STRING
    integer, intent(in) :: LENSTRING
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    integer :: n_chars
    !
    n_chars = min(len(string), lenstring)
    if ( len(string) < 1  ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Bad string arg in OUTPUT_STRING' )
    else if ( len_trim(string) < 1 .or. LENSTRING < 1  ) then
      call output_ ( '', advance )
    else
      call output_ ( string(:n_chars), advance, from_where, dont_log, log_chars )
    end if
  end subroutine OUTPUT_STRING

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
    enddo
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
    enddo
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
  ! tabn: column number where name begins
  ! tabc: column number where colon occurs
  ! taba: column number where after begins
  ! advance: whether to advance after printing pair (by default we WILL advance)
  ! dont_stamp: override setting to stamp end of each line
  subroutine output_nvp_character ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    character(len=*), intent(in)          :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_character

  subroutine output_nvp_complex ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    complex, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_complex

  subroutine output_nvp_double ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    double precision, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_double

  subroutine output_nvp_dbl_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    double precision, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_dbl_array

  subroutine output_nvp_int_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    integer, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_int_array

  subroutine output_nvp_integer ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_integer

  subroutine output_nvp_log_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    logical, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_log_array

  subroutine output_nvp_logical ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    logical, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_logical

  subroutine output_nvp_single ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_single

  subroutine output_nvp_sngl_array ( name, value, &
   & ADVANCE, colon, fillChar, Before, After, TABN, TABC, TABA, DONT_STAMP )
    character(len=*), intent(in)          :: name
    real, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine output_nvp_sngl_array

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
    endif
  end subroutine resetTabs

  ! ----------------------------------------------  restoreSettings  -----
  subroutine restoreSettings 
  ! resume outputting to PRUNIT.
    outputOptions%PRUNIT               = STDOUTPRUNIT
    outputOptions%MLSMSG_Level         = MLSMSG_Info
    outputOptions%newLineVal           = 10 ! 13
    outputOptions%nArrayElmntsPerLine  = 7
    outputOptions%nBlanksBtwnElmnts    = 3
    outputOptions%BUFFERED             = .true.
    outputOptions%OPENED               = .false.
    outputOptions%SKIPMLSMSGLOGGING    = .false.
    outputOptions%usePatternedBlanks   = .true. 
    outputOptions%specialFillChars     = '123456789'
    outputOptions%patterns             = (/ & ! on consecutive lines
                                            &  '(. )            ' , &
                                            &  '(. .)           ' , &
                                            &  '(.  .)          ' , &
                                            &  '(.   .)         ' , &
                                            &  '(.. ..)         ' , &
                                            &  '(- )            ' , &
                                            &  '(- -)           ' , &
                                            &  '(-  -)          ' , &
                                            &  '(- .. )         ' /)
                                            !   12345678901234567890
    outputOptions%name                 = 'stdout'
    outputOptions%advanceDefault       = 'no'
    outputOptions%sdFormatDefault      = '*'
    outputOptions%arrayElmntSeparator  = ' '

    stampOptions%neverStamp            = .false.
    stampOptions%post                  = .true.
    stampOptions%showTime              = .false.
    stampOptions%textCode              = ' '
    stampOptions%dateFormat            = ' '
    stampOptions%timeFormat            = 'hh:mm'
    stampOptions%interval              = 1
    stampOptions%TIMESTAMPSTYLE        = 'post'

    timeStampOptions%post              = .true.
    timeStampOptions%showDate          = .false.
    timeStampOptions%textCode          = ' '
    timeStampOptions%dateFormat        = 'yyyy-mm-dd'
    timeStampOptions%timeFormat        = 'hh:mm:ss'
    timeStampOptions%TIMESTAMPSTYLE    = 'post'
  end subroutine restoreSettings

  ! ----------------------------------------------  resumeOutput  -----
  subroutine resumeOutput 
  ! resume outputting to PRUNIT.
    silentRunning = .false.
  end subroutine resumeOutput

  ! ----------------------------------------------  revertOutput  -----
  subroutine revertOutput
  ! revert to outputting to OLDUNIT. Close current PRUNIT (if > 0 and open)
    ! Local variables
    logical :: itsOpen
    ! Executable
    call resumeOutput
    if ( .not. OLDUNITSTILLOPEN ) then
      call output_( 'Unable to Revert output--old unit not open', advance='yes' )
      return
    end if
    call output_( 'Reverting output to unit: ', advance='no' )
    call output( OLDUNIT, advance='yes' )
    if ( outputOptions%prunit > 0 ) then
      inquire( unit=outputOptions%prunit, opened=itsOpen )
      if ( itsOpen ) then
        close(outputOptions%prunit)
      end if
    end if
    outputOptions%prunit = OLDUNIT    

  end subroutine revertOutput

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
    elseif ( present(tabs) ) then
      n = min( MAXNUMTABSTOPS, size(tabs) )
      tabStops(1:n) = tabs(1:n) 
    else
      call ExpandStringRange ( '5-120+5', TABSTOPS )
    endif
  end subroutine setTabs

  ! ----------------------------------------------  suspendOutput  -----
  subroutine suspendOutput 
  ! suspend outputting to PRUNIT.
    silentRunning = .true.
  end subroutine suspendOutput

  ! ----------------------------------------------  switchOutput  -----
  subroutine switchOutput ( filename, unit, keepOldUnitOpen )
  ! stop outputting to PRUNIT. Switch to filename [using unit if supplied]
    ! Args
    character(len=*), intent(in)    :: filename
    integer, optional, intent(in)   :: unit
    logical, optional, intent(in)   :: keepOldUnitOpen
    ! Internal variables
    integer, parameter :: DEFAULTSWITCHUNIT = 59
    logical :: dontCloseOldUnit
    logical, save :: NeedToAppend = .false. ! 1st time here overwrite; append later
    integer :: switchUnit
    ! Executable
    dontCloseOldUnit = .false.
    if ( present(keepOldUnitOpen) ) dontCloseOldUnit = keepOldUnitOpen
    oldUnit = outputOptions%prunit
    switchUnit = DEFAULTSWITCHUNIT
    call resumeOutput
    if ( present(unit) ) then
      call output_('Switching further output to: ', advance='no')
      call output_(trim(filename), advance='yes')
      call output_('using unit number: ', advance='no')
      call output(unit, advance='yes')
      switchUnit = unit
    else
      call output_('Switching further output to: ', advance='no')
      call output_(trim(filename), advance='yes')
      ! if ( PRUNIT > 0 ) switchUnit = PRUNIT
    end if
    if ( outputOptions%prunit > 0 .and. .not. dontCloseOldUnit ) then
      close(outputOptions%prunit)
      OLDUNIT = outputOptions%prunit
      OLDUNITSTILLOPEN = .false.
    end if
    if ( NeedToAppend ) then
      open( unit=switchUnit, file=filename, status='old', position='append' )
    else
      open( unit=switchUnit, file=filename, status='replace' )
    end if
    outputOptions%prunit = SwitchUnit
    NeedToAppend = .true.
  end subroutine switchOutput

  ! ------------------------------------------------  timeStamp  -----
  ! time-stamp output on demand, not automatic:
  ! Either in style pre or post
  ! (pre) '(HH:MM:SS) chars'
  ! (post) 'chars (HH:MM:SS)'
  ! Note that in pre-style, the time will be printed only if ATLINESTART true
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
    logical :: my_dont_log
    character(len=8) :: my_style
    character(len=3) :: MY_ADV
    logical  ::         myDate
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_log = outputOptions%skipmlsmsglogging ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
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
      if ( ATLINESTART ) then
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
      call output_integer( INT, PLACES, &
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
      if ( ATLINESTART ) then
        call output_('(', ADVANCE='no', DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME( date=myDate, &
          & dateFormat=timeStampOptions%dateFormat, &
          & timeFormat=timeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='no', DONT_STAMP=DONT_STAMP)
      end if
      call output_integer( INT, PLACES, &
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
    endif
    call timeStamp_char(str, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, STYLE, DATE )
  end subroutine timeStamp_logical

  ! Internal procedures
  
  ! .......................................  Advance_is_yes_or_no  .....
  function Advance_is_yes_or_no ( str ) result ( outstr )
    ! takes '[Yy]...' or '[Nn..] and returns 'yes' or 'no' respectively
    ! also does the same with '[Tt]..' and '[Ff]..'
    ! leaves all other patterns unchanged, but truncated to three
    ! characters.  Returns 'no' if the argument is absent.
    
    ! We are allowing its argument str to do multiple duties by being
    ! composed of multiple space-separated arguments, e.g. 
    ! 'arg1 [arg2] .. [argn]'
    ! (1) the first arg1 is treated as before, basically 'yes' or 'no'
    ! (2) arg2 and beyond will introduce extra options to the
    ! output command, eventually simplifying it to just
    !   call output( something, [advance='arg1 [arg2] .. [argn]' )
    !--------Argument--------!
    character (len=*), intent(in), optional :: Str
    character (len=3) :: Outstr

    !----------Local vars----------!
    character (len=*), parameter :: yeses = 'YyTt'
    character (len=*), parameter :: nose = 'NnFf'
    integer :: kSpace

    if ( .not. present(str)  ) then
      outstr = outputoptions%advanceDefault ! 'no'
      return
    end if

    outstr = adjustl(str)
    kSpace = index( outstr, ' ' )
    if ( kSpace > 1 ) outstr = outstr(:kSpace) ! To snip off arg2 ..
    if ( index( yeses, outstr(:1)) > 0  ) then
      outstr = 'yes'
    else if ( index( nose, outstr(:1)) > 0  ) then
      outstr = 'no'
    else
      outstr = str
    end if
  end function Advance_is_yes_or_no

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
    endif
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

  ! .............................................  nCharsinFormat  .....
  function nCharsinFormat ( Format ) result(nplusm)
    ! Utility to calculate how many characters in a format spec:         
    ! [n{xX}][,]{DEFGdefg}m.b                                             
    ! where n, m, and b are digits (we care only about n and m)           
    ! return (n+m)
    ! Tested for specs: sci. format esm.b and eng. format enm.b
    ! Also for min. width spec: 'f0.b' it will silently return 0
    ! (It's up to you to handle that correctly)
    ! Args                                                                
    character(len=*), intent(in) ::  Format                               
    integer :: nplusm                                                     
    ! Local variables                                                     
    character(len=20) :: kChar, myFormat                                  
    integer :: n, m
    ! Executable                                                          
    nplusm = 0                                                            
    kChar=lowerCase(Format)
    call ourReplaceSubString(kChar, myFormat, 'es', 'f')                   
    call ourReplaceSubString(myFormat, kChar, 'en', 'f')                   
    call ourReplaceSubString(kChar, myFormat, 'g', 'f')                   
    call ourReplaceSubString(myFormat, kChar, 'e', 'f')                   
    call ourReplaceSubString(kChar, myFormat, 'd', 'f')                   
    call ourExtractSubString(TRIM(myFormat), kChar, 'f', '.')             
    if ( kChar == '0' ) return ! Special case of e.g. 'f0.3'
    read (kChar, '(i2)') m                                                
    if (m < 1) call MLSMessage ( MLSMSG_Error, ModuleName, &              
      & 'Bad conversion to m in OUTPUT_xxxLE (format not "{defg}"' )      
    if ( index(TRIM(myFormat), 'x' ) == 0 ) then                          
      n = 0                                                               
    else                                                                  
      call ourExtractSubString(TRIM(myFormat), kChar, '(', 'x')           
      read (kChar, '(i2)') n                                              
      if (n < 1) then                                                     
        print *, trim(kChar)                                              
        print *, trim(myFormat)                                           
        call MLSMessage ( MLSMSG_Error, ModuleName, &                     
          & 'Bad conversion to n in OUTPUT_xxxLE (format not "{defg}"' )  
      end if                                                              
    end if                                                                 
    nplusm = n + m                                                        
  end function nCharsinFormat

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
    elseif ( index(inFormat, '.') < 1 ) then
      format = sdNeedsFormat
      dotm = '.6'
    else
      ! Must find integer after '.'
      dot = index( inFormat, '.' )
      dotm = inFormat(dot:dot+1)
      format = trim(sdNeedsFragment) // trim(dotm) // ')'
    endif
  end subroutine whatSDNeedsFormat
  ! ........................................  ourExtractSubString  .....
  subroutine ourExtractSubString ( instr, outstr, sub1, sub2 )
    ! Extract portion of instr between sub1 and sub2 and return as outstr
    ! Args
    character (len=*), intent(in) :: instr
    character (len=*), intent(out) :: outstr
    character (len=1), intent(in) :: sub1
    character (len=1), intent(in) :: sub2
    ! Internal variables
    integer :: pos1
    integer :: pos2
    ! Begin executable
    outstr = ''
    pos1 = index(instr, sub1) 
    if ( pos1 < 1 ) return
    pos2 = index(instr, sub2)
    if ( pos2-1 < pos1+1 ) return
    outstr = instr(pos1+1:pos2-1)
  end subroutine ourExtractSubString

  ! ........................................  ourReplaceSubString  .....
  subroutine ourReplaceSubString ( instr, outstr, sub1, sub2 )
    ! Swap a single instance in instr of sub1 with sub2 and return as outstr
    ! Args
    character (len=*), intent(in) :: instr
    character (len=*), intent(out) :: outstr
    character (len=*), intent(in) :: sub1
    character (len=*), intent(in) :: sub2
    ! Internal variables
    integer :: d
    integer :: pos
    integer :: pos1
    integer :: pos2
    integer :: pose
    ! Begin executable
    outstr = instr
    pos =index(instr, sub1) 
    if ( pos < 1 .or. pos > len_trim(outstr)) return
    pos1 = pos
    if ( len(sub1) == len(sub2) ) then
      pos2 = pos1 + len(sub1) - 1
      outstr(pos1:pos2) = sub2
    elseif ( len(sub1) > len(sub2) ) then
      d = len(sub1) - len(sub2)
      pos2 = pos1 + len(sub2) - 1
      outstr(pos1:pos2) = sub2
      pose = min(len(instr), len(outstr))
      if ( pos2 == pose ) return
      outstr(pos2+1:) = ' '   ! To fill to the end with blanks
      outstr(pos2+1:pose) = instr(pos2+1+d:pose+d)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &              
      & 'Not yet able to replace shorter substring with longer' ) 
    end if
  end subroutine ourReplaceSubString

  ! -----------------------------------------------------  PR_BLANKS  -----
  subroutine PR_BLANKS ( N_BLANKS, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=3) :: ADV
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    character(len=len(BlankSpace)) :: b
    integer :: I    ! Blanks to write in next WRITE statement
    integer :: N    ! Blanks remaining to write
    character(len=3) :: MY_ADV
    ! Executable
    my_adv = Advance_is_yes_or_no(advance)
    if ( n_blanks < 1 ) then
      if ( my_adv == 'yes' ) &
        & call output_ ( '', advance='yes', dont_stamp=dont_stamp )
      return
    end if
    n = max(n_blanks, 1)
    if ( present(fillChar)  ) then
      do i=1, min(n, len(BlankSpace))
        b(i:i) = fillChar
      end do
    else
      b = BLANKSPACE
    end if
    adv = 'no'
    do
      i = min(n,len(b))
      n = n - i
      if ( n == 0 ) adv = my_adv
      call output_ ( b(:i), advance=adv )
      if ( n < 1 ) exit   ! was if n == 0, but this should be safer
    end do
  end subroutine PR_BLANKS

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
    enddo
  end function stretch

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

end module OUTPUT_M

! $Log$
! Revision 2.97  2012/09/11 18:52:26  pwagner
! Added isOutputSuspended
!
! Revision 2.96  2012/08/16 17:36:30  pwagner
! Removed more unused stuff
!
! Revision 2.95  2012/08/07 18:01:22  pwagner
! output simply prints tab character as is instead of as '@'
!
! Revision 2.94  2012/08/02 21:09:53  pwagner
! Added RestoreSettings
!
! Revision 2.93  2012/08/01 00:08:21  pwagner
! Uses same . . pattern in option dumps
!
! Revision 2.92  2012/07/18 00:33:45  pwagner
! Fixed bug causing double printing
!
! Revision 2.91  2012/07/17 16:38:01  pwagner
! OutputLines mechanism introduced to defer printing; new HeadLine subroutine
!
! Revision 2.90  2012/06/22 00:04:02  pwagner
! May now change default advance option to 'yes'
!
! Revision 2.89  2012/04/20 01:27:14  vsnyder
! Add CPU_Seconds to Output_Date_and_Time
!
! Revision 2.88  2011/08/11 22:24:59  pwagner
! Added banner to set off message with stars and stripes
!
! Revision 2.87  2011/07/12 00:12:25  pwagner
! Added numNeedsFormat
!
! Revision 2.86  2011/05/26 20:38:32  pwagner
! By default, outputting chars substitutes for non-ascii, except for newlines
!
! Revision 2.85  2011/04/29 02:16:32  vsnyder
! Prefill dateString and timeString with blanks to compensate for Intel bug
!
! Revision 2.84  2011/03/12 00:39:31  vsnyder
! Change len=1 to len=* to avoid Intel checking problem
!
! Revision 2.83  2010/10/14 18:43:02  pwagner
! Can now dump and reset tabs; also can outputlists
!
! Revision 2.82  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.81  2010/01/26 17:49:42  pwagner
! Fixed bug that added space before newlines; simplified output_char
!
! Revision 2.80  2009/06/24 22:35:44  pwagner
! Trick to pass places arg into output via format arg
!
! Revision 2.79  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.78  2009/06/16 17:22:42  pwagner
! With outputOptions may make outputted arrays look like text Climatology file
!
! Revision 2.77  2009/05/14 22:00:28  pwagner
! New optional arg onlyif prints logicals only if true (false)
!
! Revision 2.76  2008/11/24 19:29:43  pwagner
! Added a print to not_used_here
!
! Revision 2.75  2008/10/17 00:04:22  pwagner
! NewLine should not add space at line end when output to named file
!
! Revision 2.74  2008/06/17 00:01:53  pwagner
! Separate options for auto and manual time stamping
!
! Revision 2.73  2008/05/02 00:07:20  pwagner
! Correctly handles some rare cases, newLine wont add extra space
!
! Revision 2.72  2008/04/18 16:34:37  pwagner
! achar(13) among chars by default now triggers newLine
!
! Revision 2.71  2008/03/07 01:34:57  pwagner
! Added f.p. array versions of generic outputNamedValue
!
! Revision 2.70  2008/02/22 21:24:14  pwagner
! Lets NAG-built tools dump PCF, L2CF
!
! Revision 2.69  2008/01/09 20:52:03  pwagner
! call output(NaN) now prints 'NaN'; same with Inf
!
! Revision 2.68  2007/12/19 01:29:05  pwagner
! Removed unused variables
!
! Revision 2.67  2007/11/30 18:19:48  pwagner
! outputCalendar handles yyyy-mm-dd formatted date
!
! Revision 2.66  2007/10/18 23:39:46  pwagner
! Added numToChars and alignToFit intercaes for numeric types
!
! Revision 2.65  2007/09/24 20:22:08  pwagner
! Improved outputCalendar
!
! Revision 2.64  2007/09/20 17:38:09  pwagner
! improved outputCalendar; neverStamp field added to stampOptions
!
! Revision 2.63  2007/09/14 00:15:42  pwagner
! Added alignToFit and outputCalendar
!
! Revision 2.62  2007/09/06 22:27:06  pwagner
! Renamed TAB to blanksToTab to avoid conflict with TOGGLES constant
!
! Revision 2.61  2007/08/27 23:55:01  pwagner
! Added many tabstop-related procedures
!
! Revision 2.60  2007/07/27 00:21:59  vsnyder
! Spiff up printing in DumpSize
!
! Revision 2.59  2007/07/17 00:24:18  pwagner
! Treat certain numbers with default list-directed format
!
! Revision 2.58  2007/06/14 18:40:04  pwagner
! Allow sdFormatDefault to be set at class level
!
! Revision 2.57  2007/04/14 00:37:16  vsnyder
! Correction dumpSize to avoid asterisks in I6 output
!
! Revision 2.56  2007/03/23 00:08:21  pwagner
! Guard against negative args confusing dumpSize
!
! Revision 2.55  2007/01/13 01:49:48  pwagner
! Repaired long-standing bug blighting logged output
!
! Revision 2.54  2007/01/12 00:31:00  pwagner
! May use unbuffered output ;renamed routine outputNamedValue
!
! Revision 2.53  2006/07/28 01:58:53  vsnyder
! Cannonball polishing in dumpSize routines
!
! Revision 2.52  2006/07/19 22:25:14  vsnyder
! Add Dumpsize_Double, plus some cannonball polishing
!
! Revision 2.51  2006/06/27 23:58:01  pwagner
! name_v_value works much better
!
! Revision 2.50  2006/06/24 23:05:07  pwagner
! Added outputNamedValue, special blank fills
!
! Revision 2.49  2006/06/03 00:17:29  vsnyder
! Eliminate trailing blanks sometimes
!
! Revision 2.48  2006/02/21 19:08:36  pwagner
! Removed two unused declarations
!
! Revision 2.47  2006/02/15 18:10:44  pwagner
! Fixed bug preventing CR from being printed sometimes
!
! Revision 2.46  2006/02/15 00:00:07  pwagner
! Added automatic stamping features
!
! Revision 2.45  2006/02/10 21:23:05  pwagner
! Added switchOutput, revertOutput
!
! Revision 2.44  2006/01/04 20:28:51  pwagner
! Added suspend- and resumeOutput procedures
!
! Revision 2.43  2005/12/16 23:25:13  pwagner
! dumpSize moved from dump0 to output_m
!
! Revision 2.42  2005/09/22 23:34:56  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.41  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.40  2005/03/19 01:14:41  pwagner
! Fixed error in formatting output_date_and_time
!
! Revision 2.39  2005/01/20 23:17:51  pwagner
! Prevent n_char being gt len(my_chars)
!
! Revision 2.38  2005/01/19 01:09:49  pwagner
! New timeStamp interface to certain output procedures
!
! Revision 2.37  2005/01/07 01:26:04  pwagner
! Advance now an optional arg to OUTPUT_DATE_AND_TIME so it can time-stamp
!
! Revision 2.36  2004/12/31 02:39:51  vsnyder
! Simplified computing My_Adv, simplified Output_Char, added Before argument
! to Output_Logical, some cannonball-polishing.
!
! Revision 2.35  2004/12/28 19:29:57  pwagner
! Changes to handle formats like f0.3, en10.2, es8.2
!
! Revision 2.34  2004/12/14 00:00:50  pwagner
! Optional arg insteadofblank added to char outputs
!
! Revision 2.33  2004/12/13 20:30:19  vsnyder
! Cosmetic cannonball polishing
!
! Revision 2.32  2004/09/23 22:57:36  pwagner
! Added output_date_and_time
!
! Revision 2.31  2004/08/04 23:19:02  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.30  2004/06/10 23:59:29  pwagner
! blanks may take optional fillchar
!
! Revision 2.29  2004/02/26 21:51:15  pwagner
! Added output_string--although it is almost useless
!
! Revision 2.28  2003/10/07 01:12:59  vsnyder
! Add NewLine subroutine, and Before and After text args
!
! Revision 2.27  2003/09/15 23:08:44  vsnyder
! Remove five unused local variables
!
! Revision 2.26  2003/09/08 17:43:25  pwagner
! Fixed bug in nCharsinFormat when no 'x' in Format
!
! Revision 2.25  2003/09/06 01:35:55  pwagner
! Can account for (nx,{defg}m.b} in f.p. format
!
! Revision 2.24  2003/08/25 17:48:37  pwagner
! Remembered formats may be gx.y
!
! Revision 2.23  2003/08/25 17:06:50  pwagner
! Remembered that formats may use ex.y as well as [fd]x.y
!
! Revision 2.22  2003/08/23 00:11:46  pwagner
! Tried to fix prob with fix to output_single; also output_double
!
! Revision 2.21  2003/08/21 21:20:35  cvuu
! Change output of format in OUTPUT_SINGLE
!
! Revision 2.20  2003/07/02 01:07:27  vsnyder
! Add complex output
!
! Revision 2.19  2003/03/20 19:20:17  pwagner
! Changes to prevent double-logging when using MLSMessage
!
! Revision 2.18  2003/02/27 18:35:30  pwagner
! Appends trailing spaces to improve appearance with MLSMessage
!
! Revision 2.17  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.16  2001/10/19 22:31:36  pwagner
! Now can output (small-sized) s.p., d.p. arrays
!
! Revision 2.15  2001/10/08 23:43:28  pwagner
! Allows wider range of advance(s); my_adv implemented uniforml
!
! Revision 2.14  2001/09/26 02:16:22  vsnyder
! Simplify by using output_char internally
!
! Revision 2.13  2001/05/24 22:39:07  vsnyder
! Make output_single work like output_double; cosmetic changes
!
! Revision 2.12  2001/05/24 22:22:48  vsnyder
! Add Output_Single to the generic Output interface
!
! Revision 2.11  2001/05/10 22:52:03  vsnyder
! Increase maximum integer width
!
! Revision 2.10  2001/05/10 18:22:00  pwagner
! Added LogFOrmat to output_double
!
! Revision 2.9  2001/05/08 20:27:24  vsnyder
! Added an optional 'format' argument in a few more places
!
! Revision 2.8  2001/04/25 00:08:01  vsnyder
! Add 'fill' argument to 'output_integer'
!
! Revision 2.7  2001/04/18 23:28:10  pwagner
! Added output_integer_array
!
! Revision 2.6  2001/04/07 01:53:28  vsnyder
! Output 0 instead 0.0e+00
!
! Revision 2.5  2001/03/16 23:14:16  vsnyder
! Don't trim off the last nonzero digit
!
! Revision 2.4  2001/02/28 21:35:34  livesey
! Added output logical
!
! Revision 2.3  2001/02/22 23:54:27  vsnyder
! Added optional "from_where" argument to "output_char"
!
! Revision 2.2  2001/02/22 23:27:16  vsnyder
! Correct routing of output through MLSMessage
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.2  2000/10/09 23:29:54  vsnyder
! Must have updated something -- permissions weren't r--r--r--.
!
! Revision 2.1  2000/10/04 18:07:04  vsnyder
! Added capability to output through MLSMessage
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
