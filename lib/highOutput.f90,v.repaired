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
  
  ! See also Dump_0 and Output_M
  
  use Dates_Module, only: BuildCalendar, DaysInMonth, &
    & ReformatDate, ReformatTime, Utc_To_Yyyymmdd
  use Machine, only: Crash_Burn, Exit_With_Status, NeverCrash
  use MLSCommon, only: LineLen, MLSDebug, MLSVerbose
  use MLSFinds, only: FindFirst, FindNext
  use MLSStringLists, only: ExpandStringRange, GetStringElement, &
    & List2Array, NumStringElements, SwitchDetail, Wrap
  use MLSStrings, only: Asciify, Indexes, Justify, Lowercase, &
    & NCharsInFormat, Ncopies, &
    & Replace, Stretch, Trim_Safe, WriteIntsToChars
  use Output_M, only: Advance_Is_Yes_Or_No, Blanks, GetOutputStatus, &
    & MaxOutputLineslen, Newline, &
    & Output, Output_ => Output_Char_Nocr, &
    & RestoreOutputSettings => RestoreSettings, &
    & OutputOptions, OutputOptions_T, PatternOptions_T, &
    & StampOptions, StampOptions_T, &
    & TimeStampOptions, TimeStampOptions_T, &
    & BothPrUnit, InvalidPrUnit, MSGLogPrUnit, &
    & OutputLines, OutputLinesPrUnit, SetOutputStatus, StdoutPrUnit
  use PrintIt_M, only: AssembleFullLine, Get_Config, &
    & MLSMSG_Crash, MLSMSG_Debug, &
    & MLSMSG_Severity_To_Quit, &
    & MLSMSG_Warning, &
    & PrintItOut, MLSMessageConfig, SeverityNamesFun
  use Toggles, only: Switches
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (datatypes)
! StyleOptions             How Banners and Headlines will appear
!     (subroutines and functions)
! AddRow                   Add a {name, value} row to a 2-d table of cells
! AddRow_header            Add a single line stretched across an entire row
! AddRow_divider           Add a row composed of a single, repeated character
! AlignToFit               Align printed argument to fit column range
! Banner                   Surround message with stars and stripes; e.g.,
!                            *-----------------------------------------------*
!                            *            Your message here                  *
!                            *-----------------------------------------------*
! BeVerbose                Should we print extra, optional data?
! BlanksToColumn           Print blanks [or fill chars] out to specified column
! BlanksToTab              Print blanks [or fill chars] out to next tab stop
! Dump                     Dump output, pattern, or stamp options
! Dumpsize                 Print a nicely-formatted memory size 
! Dumptabs                 Print the current tab stop positions
! FinalMemoryReport        Print a summary of memory allocated/deallocated
! GetStamp                 Get stamp being added to every output
! HeadLine                 Print a line with eye-catching features
!                           e.g., '*-------  Your message here   -------*'
! LetsDebug                Should we print extra debugging data?
! NextColumn               Return next column number that would be printed
! NextTab                  Return next column number that Tab would move to
! NumNeedsFormat           Return what format is needed to output num
! NumToChars               Return what string would be printed by output
! OutputCalendar           Output nicely-formatted calendar page
! Output_Date_And_Time     Print nicely formatted date and time
! OutputList               Output array as comma-separated list; e.g. '(1,2,..)'
! OutputNamedValue         Print nicely formatted name and value
! OutputParagraph          Print text formatted as a paragraph with line breaks, 
!                            indent, etc.
! OutputTable              Output 2-d array as cells in table,; or else
!                          Output the 2d table of cells constucted by AddRow(s)
! ResetTabs                Restore tab stops to what was in effect at start
! RestoreSettings          Restore default settings for output, styles, tabs
! SetStamp                 Set stamp to be automatically printed on every line
! SetTabs                  Set tab stops (to be used by tab)
! StartTable               Initialize a 2-d table of cells to be output later
! StyledOutput             Output a line according to options; e.g. "--Banner"
! Tab                      Move to next tab stop
! Timestamp                Print argument with a timestamp manually
!                            (both stdout and logged output)
! === (end of toc) ===

! === (start of api) ===
! AddRow ( char* name, value, [int BlocLen], [char* options], [char* format] )
! AddRow_header ( char* name, char alignment )
! AddRow_divider ( char char )
! AlignToFit ( char* chars, int columnRange(2), char alignment, [int skips] )
! Banner ( char* chars, [int columnRange(2)], [char alignment], [int skips], 
!    [int lineLength], [char mode], [char pattern], [log underline] )
! log BeVerbose ( char* switch[(:)], threshold )
! BlanksToColumn ( int column, [char fillChar], [char* advance] )
! BlanksToTab ( [int tabn], [char* fillChar] )
! Dump ( options )
! DumpSize ( n, [char* advance], [units] )
!       where n can be an int or a real, and 
!       units is a scalar of the same type, if present
! DumpTabs ( [int tabs(:)] )
! FinalMemoryReport ( [log IsFinal] )
! GetStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! HeadLine ( char* chars, 
!          [char fillChar], [char* Before], [char* After], 
!          [int columnRange(2)], [char alignment], [int skips], [log underline] )
! log LetsDebug ( char* switch[(:)], threshold )
! int NextColumn ( )
! int NextTab ( )
! char* NumNeedsFormat ( value )
! char* NumToChars ( value, [char* format] )
! Output_Date_And_Time ( [log date], [log time], [char* from_where], 
!          [char* msg], [char* dateFormat], [char* timeFormat], 
!          [double CPU_seconds], [int wallClock_seconds], [char* advance] )
! OutputCalendar ( [char* date], [char* datenote], [char* notes(:)], 
!          [log dontwrap, [log moonPhases] ] )
! OutputList ( values(:), [char* sep], [char* delims] )
! OutputNamedValue ( char* name, value, [char* advance],
!          [char colon], [char fillChar], [char* Before], [char* After], 
!          [integer tabn], [integer tabc], [integer taba], log dont_stamp],
!          [char* options] )
! OutputParagraph ( char* text, int ColumnRange(2), [char* alignment], &
!          [int indent] )
! OutputTable ( [array(:,:)], [char sep], [char border], [int cellWidth],
!          [char interior], [char headliner], [char alignment] )
! ResetTabs ( [int tabs(:)] )
! RestoreSettings ( [char* settings] )
! SetStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! SetTabs ( [char* Range], [int tabs(:)] )
! StartTable
! StyledOutput ( char* chars, [char* options] )
! Tab ( [int tabn], [char* fillChar] )
! TimeStamp ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [char*8 style], [log date] )
! TimeStamp ( log value, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [char*8 style], [log date] )
! TimeStamp ( int int, [int places], [char* advance],
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
! E.g., calling TimeStamp or else using SetStamp before calling output
!
! To understand the codes for DateFormat and TimeFormat, see the Dates_Module
! 
! Here's how to build and print a 2d table of names and values
! (1) Call startTable
! (2) Optionally call AddRow_header
! (3) Optionally call AddRow_divider
! (4) For each name, value pair
!     (a) call AddRow
!     (b) call AddRow_Divider
! (5) Call OutputTable
!
! (see CellDatabase)
!
! The aligment arg in AlignToFit, etc.  can be explained best with an example,
! in fact 4 different examples (showing L, R, C, and J in that order)
! ------------------------------------------------------------------------------
! The first line is L                                                               
!                                                           The second line is R
!                               The third line is C                                
! The                                 final                                 line
! ------------------------------------------------------------------------------
! The same 4 styles can be applied to Banner and StyledOutput by
! suitable choices of the options arg.
!
! A better programmer would make a '--help' option
! available for most of these procedures
! and maybe the whole module, too.

! A general rule to be remembered always:
!    * There's no software that can't be improved *
! which is a corollary of the more biting and disheartening rule
!    * There's no software that is free of bugs *


  public :: AddRow, AddRow_Divider, AddRow_Header, AlignToFit, &
    & Banner, BeVerbose, BlanksToColumn, BlanksToTab, &
    & Dump, DumpSize, DumpTabs, FinalMemoryReport, GetStamp, HeadLine, &
    & LetsDebug, NextColumn, NextTab, NumNeedsFormat, NumToChars, &
    & Output_Date_And_Time, OutputCalendar, OutputList, OutputTable, &
    & OutputAnyNamedValue, OutputNamedValue, OutputParagraph, &
    & ResetTabs, RestoreSettings, &
    & SetStamp, SetTabs, StartTable, StyledOutput, Tab, TimeStamp

  ! These types must be made public because the class instances are public
  public :: OutputOptions_T
  public :: StampOptions_T
  public :: TimeStampOptions_T

  interface AddRow
    module procedure AddRow_Character, AddRow_Character_Blocs
    module procedure AddRow_Complex
    module procedure AddRow_Dbl_Array, AddRow_Double
    module procedure AddRow_Int_Array, AddRow_Integer
    module procedure AddRow_Log_Array, AddRow_Logical
    module procedure AddRow_Sngl_Array, AddRow_Single
  end interface

  interface AlignToFit
    module procedure AlignToFit_Chars, AlignToFit_Double, AlignToFit_Single
    module procedure AlignToFit_Integer
  end interface

  interface Banner
    module procedure Banner_Chars
    module procedure Banner_Chararray
  end interface

! Other modules may interrogate Beverbose or LetsDebug to decide whether or not
! to print some intermediate results or messages.
! Do these really belong here, or somewhere else?
  interface Beverbose
    module procedure Beverbose_Chars
    module procedure Beverbose_Chararray
  end interface

  interface Letsdebug
    module procedure Letsdebug_Chars
    module procedure Letsdebug_Chararray
  end interface

  interface Dump
    module procedure DumpOutputOptions, DumpPatternOptions, &
      & DumpStampOptions, DumpTimeStampOptions
  end interface

  interface Dumpsize
    module procedure Dumpsize_Double, Dumpsize_Integer, Dumpsize_Real
  end interface

  interface Getoption
    module procedure Getoption_Char, Getoption_Log
  end interface

  interface NumNeedsFormat
    module procedure NumNeedsFormat_Double, NumNeedsFormat_Integer, NumNeedsFormat_Single
    module procedure NumNeedsFormat_Complex, NumNeedsFormat_Dcomplx
  end interface

  interface NumToChars
    module procedure NumToChars_Double, NumToChars_Integer, NumToChars_Single
  end interface

  interface OutputList
    module procedure OutputList_Ints, OutputList_Chars
  end interface

  interface OutputnamedValue
    module procedure Output_Nvp_Character
    module procedure Output_Nvp_Complex
    module procedure Output_Nvp_Dbl_Array, Output_Nvp_Double
    module procedure Output_Nvp_Int_Array, Output_Nvp_Integer
    module procedure Output_Nvp_Log_Array, Output_Nvp_Logical
    module procedure Output_Nvp_Sngl_Array, Output_Nvp_Single
  end interface

  interface OutputanynamedValue
    module procedure Output_Nvp_Whatever
  end interface

  interface Tab
    module procedure Blankstotab
  end interface
  
  interface TimeStamp
    module procedure TimeStamp_Char, TimeStamp_Integer, TimeStamp_Logical
  end interface
  
  ! -------------------------------------------------------------------------
  ! *      Module settings, parameters, and data types                      *
  ! When Calling OutputNamedValue with character values, should we trim them?
  logical, public        :: TrimCharacterValues = .true.
  logical, save, private :: OldNeverStamp
  integer, save, private :: OLDWRAPPASTCOLNUM = 0
  integer, save, private :: WRAPPASTCOLNUM = 0  ! Don't print beyond (if > 0)
  
  ! -------------------------------------------------------------------------
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
  character(len=MAXCELLSIZE), dimension(:,:), pointer :: CellDatabase => null()
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ! * Tabs                                                                    *
  integer, private, parameter :: MAXNUMTABSTOPS = 24
  ! These next tab stops can be reset using the procedure setTabs
  ! the default values correspond to range coded '5-120+5'
  ! (read as from 5 to 120 in intervals of 5)
  character(len=*), parameter :: INITTABRANGE = '5-120+5'
  integer, dimension(MAXNUMTABSTOPS), save, private :: TABSTOPS = &
    & (/ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, &
    &   65, 70, 75, 80, 85, 90, 95,100,105,110,115,120 /)

  ! -------------------------------------------------------------------------
  ! For certain numerical values we will use list directed '*' format
  ! unless optional FORMAT specifier supplied
  double precision, parameter, dimension(3) :: DPREFERDEFAULTFORMAT = &
    & (/ -1.d0, 0.d0, 1.d0 /)  ! For which values to use default format '*'
  character(len=12), private :: sdNeedsFormat = '(1pg14.6)'
  character(len=12), private :: sdNeedsFragment = '(1pg14'

  ! -------------------------------------------------------------------------
  ! This is the type for configuring how to automatically style 
  ! special output formats; e.g., Banner
  ! Note the effect on the "bars" part of "stars and bars"
  ! of choosing different HeadLineFill or BannerPattern characters:
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
    logical                            :: Underline           = .false.
    ! HeadLine
    integer, dimension(2)              :: HeadLineColumnrange = (/ 1, 80 /)
    integer                            :: HeadLineSkips       = 0
    character(len=1)                   :: HeadLineAlignment   = 'C'
    character(len=1)                   :: HeadLineFill        = ' '
    character(len=8)                   :: HeadLineBefore      = ' '
    character(len=8)                   :: HeadLineAfter       = ' '
    logical                            :: HeadlineStretch     = .false.
    ! Banner
    integer, dimension(2)              :: BannerColumnrange   = (/ 1, 80 /)
    integer                            :: BannerSkips         = 0
    character(len=1)                   :: BannerAlignment     = 'C'
    character(len=1)                   :: BannerPattern       = '-'
    integer                            :: BannerLength        = 0
    logical                            :: BannerStretch     = .false.
  end type
  type(StyleOptions_T), private, save  :: DefaultStyleOptions
  type(StyleOptions_T), public, save   :: StyleOptions

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------  AddRow  -----
  ! This family of routines pairs each added name with its value
  ! by inserting a new row into the CellDatabase
  
  ! later to be printed as a neatly-formatted table to stdout
  ! By means of optional args you can create a line like
  ! *   name                   value   *
  subroutine AddRow_character_blocs ( name, value, BlocLen, &
    & options, format, wrappingChar )
    ! For character values of great length, the values may
    ! span multiple lines, which we'll call "blocs"
    ! Method: divide value up into separate blocs, each BlocLen long
    ! The 1st bloc will have its name in the name cell,
    ! later blocs won't.
    ! If options is present and contains the character 'w', the blocs
    ! will be divided at spaces (or wrappingChar if present) 
    ! instead of arbitrarily
    character(len=*), intent(in)             :: name
    character(len=*), intent(in)             :: value
    integer, intent(in)                      :: BlocLen ! How long is each bloc?
    character(len=*), intent(in), optional   :: options
    character(len=*), intent(in), optional   :: format
    character(len=*), intent(in), optional   :: wrappingChar
    ! Local variables
    integer                                  :: BlocLength
    character(len=1)                         :: break
    integer, parameter                       :: MaxNBlocs = 50
    integer                                  :: c1, c2 ! 1st and last positions
    integer                                  :: i
    character(len=1)                         :: null
    character(len=len(name))                 :: itsName
    integer                                  :: NBlocs ! How many blocs?
    character(len=len(value)+51)             :: wrapped
    logical                                  :: wrapValue
    logical, parameter                       :: DeeBug = .false.
    ! Executable
    BlocLength = max( BlocLen, 1)
    BlocLength = min(MAXCELLSIZE, BlocLength)
    wrapValue = .false.
    if ( present(options) ) wrapValue = (index( options, 'w' ) > 0)
    if ( wrapValue .and. any( indexes( trim(value), (/'"',  "'"/) ) > 0 ) ) then
      call myMessage ( MLSMSG_Warning, 'AddRow_character_blocs', &
        & 'Value contains quoted material--may not wrap properly' )
      ! wrapValue = .false.
    endif
    nBlocs = ( len_trim(value ) - 1)/BlocLength + 1
    if ( wrapvalue ) then
      ! Instead of wrapping hard at fixed bloc sizes, 
      ! try to wrap soft at spaces
      null = achar(10)
      break = ' '
      if ( present(wrappingChar) ) break = wrappingChar
      if ( DeeBug ) call output( trim(value), advance='yes' )
      call wrap( trim(value), wrapped, BlocLength, inseparator=null, &
        & mode='soft', break=break, addedLines=nBlocs )
      wrapped = adjustl(wrapped) ! Still unsure why this is necessary
      if ( DeeBug ) call output( trim(wrapped), advance='yes' )
      nBlocs = nBlocs + 1
      if ( DeeBug ) call outputNamedValue ( 'nBlocs', nBlocs )
      itsName = name
      c2 = 0
      do i=1, nBlocs
        c1 = c2 + 1
        if ( i < 2 ) then
          c2 = FindFirst( wrapped, null )
        elseif ( c1 > len_trim(wrapped) ) then
          exit
        else
          c2 = FindNext( wrapped, null, current=c1 )
        endif
        c2 = c2 - 1 ! Back up to last non-breaking character
        if ( wrapped(c1:c1) == null ) c1 = c1 + 1
        if ( c2 < c1 ) then
          c2 = len_trim(wrapped)
        endif
        call AddRow_character ( itsName, wrapped(c1:c2), format )
        if ( DeeBug ) call output( (/c1, c2/), advance='yes' )
        if ( DeeBug ) call output( wrapped(c1:c2), advance='yes' )
        itsName = ' ' ! A trick so blocs after 1st are nameless
      enddo
      return
    endif
    nBlocs = min( MaxNBlocs, NBlocs )
    itsName = name
    c2 = 0 ! A trick making 1st bloc start at position 1
    do i=1, nBlocs
      c1 = c2 + 1
      c2 = min(c2 + BlocLength, len_trim(value) )
      call AddRow_character ( itsName, value(c1:c2), format )
      itsName = ' ' ! A trick so blocs after 1st are nameless
    enddo
  end subroutine AddRow_character_blocs

  subroutine AddRow_character ( name, value, format )
    character(len=*), intent(in)          :: name
    character(len=*), intent(in)          :: value
    include 'addRow.f9h'
  end subroutine AddRow_character

  subroutine AddRow_complex ( name, value, format )
    character(len=*), intent(in)          :: name
    complex, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine AddRow_complex

  subroutine AddRow_double ( name, value, format )
    character(len=*), intent(in)          :: name
    double precision, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine AddRow_double

  subroutine AddRow_dbl_array ( name, value, format )
    character(len=*), intent(in)          :: name
    double precision, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine AddRow_dbl_array

  subroutine AddRow_int_array ( name, value, format )
    character(len=*), intent(in)          :: name
    integer, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine AddRow_int_array

  subroutine AddRow_integer ( name, value, format )
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine AddRow_integer

  subroutine AddRow_log_array ( name, value, format )
    character(len=*), intent(in)          :: name
    logical, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine AddRow_log_array

  subroutine AddRow_logical ( name, value, format )
    character(len=*), intent(in)          :: name
    logical, intent(in)                   :: value
    include 'addRow.f9h'
  end subroutine AddRow_logical

  subroutine AddRow_single ( name, value, format )
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    include 'addRow.f9h'
  end subroutine AddRow_single

  subroutine AddRow_sngl_array ( name, value, format )
    character(len=*), intent(in)          :: name
    real, dimension(:), intent(in)     :: value
    include 'addRow.f9h'
  end subroutine AddRow_sngl_array

  subroutine AddRow_header ( name, alignment )
    character(len=*), intent(in)          :: name
    character(len=1), intent(in)          :: alignment !: 'l(eft)', 'c', or 'r'
    character(len=MAXCELLSIZE), dimension(2)   :: item = ' '
    integer                                    :: newSize
    item(1) = '<<' // alignment // '>>' // name
    newSize = addCellRowToDatabase( cellDatabase, item )
  end subroutine AddRow_header

  subroutine AddRow_divider ( char )
    character(len=1), intent(in)          :: char
    character(len=MAXCELLSIZE), dimension(2)   :: item = ' '
    integer                                    :: newSize
    item(1) = '<<' // 'd' // '>>' // char
    newSize = addCellRowToDatabase( cellDatabase, item )
  end subroutine AddRow_divider

  ! -----------------------------------------------------  AlignToFit  -----
  ! Align chars to fit within column range
  ! Alignment controls whether the chars are
  ! L    Flushed left
  ! R                                              Flushed right
  ! C                 Centered
  ! J    Justified       (padding spaces to any existing spaces)
  subroutine AlignToFit_Chars ( chars, columnrange, alignment, skips )
    character(len=*), intent(in)                :: Chars
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), intent(in)           :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    integer, optional, intent(in)               :: Skips ! How many spaces between chars
    !
    ! Internal variables
    character(len=max(len(chars), abs(columnRange(2)-columnRange(1)))) :: &
      & Allchars
    character(len=max(len(chars), abs(columnRange(2)-columnRange(1)))) :: &
      & Justified
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
        allChars = stretch( chars, skips, options='a' )
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
    ! Why was this necessary? Are there nulls?
    ! if ( nc > 1 .and. allchars(nc:nc) == ' ' ) nc = nc - 1
    select case (lowercase(alignment))
    case ('l')
      char1    = 1
      padLeft  = 0
      char2    = min( nc, spaces )
      padRight = spaces - char2
    case ('r')
      char1    = max(1, nc-spaces+1)
      char2    = nc
      padLeft  = spaces - (char2-char1)
      padRight = 0
    case ('j')
      ! print *, 'columnRange: ', columnRange
      nc = abs(columnRange(2) - columnRange(1)) + 1
      Justified = Justify( allChars, nc )
      call output_( Justified(1:nc) )
      return
    case ('c')
    ! case ('c', 'j')
      m = (spaces - nc) / 2
      padLeft  = max( m, 0 )
      padRight = max( spaces - nc - m, 0 )
      char1 = max(1-m, 1)
      char2 = min(nc+m, nc)
      if ( lowercase(alignment) == 'j' .and. padRight > 0 ) &
        & firstSpace = index( allChars, ' ' )
    end select
    ! print *, 'char1, char2, padLeft, padRight, firstSpace ', &
    !   & char1, char2, padLeft, padRight, firstSpace
    ! print *, 'nc, alignment, char(char2): ', &
    !   & nc, lowercase(alignment), allChars(char2:char2)
    if ( firstSpace > 1 ) then
      call output_( allChars(char1:firstSpace-1) )
      call blanks( padRight+padLeft+1 )
      if ( firstSpace+1 < char2 ) call output_( allChars(firstSpace+1:char2) )
    else
      call blanks( padLeft )
      call output_( allChars(char1:char2) )
      call blanks( padRight )
    end if
  end subroutine AlignToFit_Chars

  subroutine AlignToFit_Double ( value, columnrange, alignment, format )
    double precision, intent(in)                :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    character(len=*), optional, intent(in)      :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call AlignToFit( trim(line), columnRange, alignment )
  end subroutine AlignToFit_Double

  subroutine AlignToFit_Integer ( value, columnrange, alignment, format )
    integer, intent(in)                         :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    character(len=*), optional, intent(in)      :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call AlignToFit( trim(line), columnRange, alignment )
  end subroutine AlignToFit_Integer

  subroutine AlignToFit_Single ( value, columnrange, alignment, format )
    real, intent(in)                            :: value
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    character(len=*), optional, intent(in)      :: FORMAT
    !
    ! Internal variables
    character(len=30) :: line
    ! Executable
    line = numToChars( value, format )
    call AlignToFit( trim(line), columnRange, alignment )
  end subroutine AlignToFit_Single

  ! -----------------------------------------------------  Banner  -----
  ! Surround your message with stars and stripes; e.g.,
  ! *-----------------------------------------------*
  ! *            Your message here                  *
  ! *-----------------------------------------------*
  ! proclaiming its great importance to an uncaring world.
  ! For multiline messages, you may divide them into elements of
  ! a character array, or else a longer character scalar and
  ! supply LineLength asking the routine to wrap at word boundaries
  !
  ! Another way to combine multiple mesgs is like this
  !   call Banner ( first, Bottomless=.true. )
  !   call Banner ( second, Bottomless=.true., Topless=.true. )
  !     .   .   .
  !   call Banner ( last, Topless=.true. )
  !
  ! See also HeadLine, StyledOutput
  subroutine Banner_Chars ( inChars, &
    & columnRange, alignment, skips, lineLength, &
      & mode, pattern, underline, Stretched, Topless, Bottomless )
    character(len=*), intent(in)                :: inChars
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    integer, optional, intent(in)               :: Skips ! How many spaces between chars
    integer, optional, intent(in)               :: Linelength
    character (len=*), optional, intent(in)     :: mode ! if not 'hard'
    character (len=1), optional, intent(in)     :: pattern ! if not stripes
    logical, optional, intent(in)               :: underline ! beneath non-blank chars
    logical, optional, intent(in)               :: Stretched   ! s t r e t c h
    logical, optional, intent(in)               :: Topless     ! skip top stripe
    logical, optional, intent(in)               :: Bottomless  ! skip bottom
    !
    ! Internal variables
    integer                          :: addedLines
    character(len=160), dimension(:), allocatable :: lines
    character(len=2*len(inchars))    :: Chars
    integer                          :: lineLen, mySkips, padding
    character(len=1)                 :: myAlignment              
    character(len=1)                 :: myFillChar               
    integer, dimension(2)            :: myColumnRange            
    logical                          :: myStretched
    logical                          :: myUnderline              
    logical, parameter               :: DEBUG = .false.          
    logical                          :: skipTop
    logical                          :: skipBottom
    character(len=len(chars))        :: underlineChars
    character(len=2*len(chars))      :: wrappedChars
    ! Executable
    Chars = ' '
    myAlignment = StyleOptions%BannerAlignment ! 'C'
    if ( present(alignment) ) myAlignment = alignment
    mySkips = StyleOptions%BannerSkips ! 0
    if ( present(skips) ) mySkips = skips
    myFillChar = StyleOptions%BannerPattern ! '-'
    if ( present(pattern) ) myFillChar = pattern
    lineLen = StyleOptions%BannerLength ! 0
    myStretched = .false.
    if ( present(Stretched) ) myStretched = Stretched
    myUnderline = .false.
    if ( present(underline) ) myUnderline = underline

    if ( present(LineLength) ) lineLen = LineLength
    if ( myStretched ) then
      Chars = Stretch ( inChars, 1, options='a' )
    else
      Chars = inChars
    endif
    skipTop = .false.
    if ( present(Topless) ) skipTop = Topless
    skipBottom = .false.
    if ( present(Bottomless) ) skipBottom = Bottomless
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
    elseif ( myAlignment == 'C' ) then
      lineLen = max( 80, 4 + len_trim(chars)*(1+mySkips) )
      padding = ( lineLen - len_trim(chars)*(1+mySkips) ) / 2
      myColumnRange(1) = 1 + padding
      myColumnRange(2) = lineLen - padding
    else
      lineLen = max( 80, 4 + len_trim(chars)*(1+mySkips) )
      myColumnRange(1) = 1
      myColumnRange(2) = LineLen
    end if
    
    ! define padding as the larger of columnrange(1) and 1
    padding = max( 1, myColumnRange(1) )
    LineLen = padding + myColumnRange(2) - 1
    if ( DEBUG ) then
      call outputnamedValue( 'padding', padding )
      call outputnamedValue( 'LineLen', LineLen )
      call outputnamedValue( 'myColumnRange', myColumnRange )
    end if
    ! Here we begin printing
    ! Temporarily stop Stamping lines
    OldNeverStamp = stampOptions%neverStamp
    stampOptions%neverStamp = .true.
    ! Top border
    if ( .not. SkipTop ) then
      call output( '*' )
      call blanks ( lineLen-2, FillChar=myFillChar )
      call output( '*', advance = 'yes' )
    endif
    ! Left star, then message, then right star
    call output( '*' )
    call AlignToFit( chars, myColumnRange, myAlignment, skips )
    call BlanksToColumn( lineLen )
    call output( '*', advance = 'yes' )
    if ( myUnderline .and. len_trim(chars) > 0 ) then
      ! The next trick replaces everything with a dash except spaces
      underlineChars = Replace( chars, ' ', '-', reverse=.true. )
      call output( '*' )
      call AlignToFit( underlineChars, myColumnRange, myAlignment, skips )
      call BlanksToColumn( lineLen )
      call output( '*', advance = 'yes' )
    endif
    ! Bottom border
    if ( .not. SkipBottom ) then
      call output( '*' )
      call blanks ( lineLen-2, FillChar=myFillChar )
      call output( '*', advance = 'yes' )
    endif
    ! Restore Stamping
    stampOptions%neverStamp = OldNeverStamp
  end subroutine Banner_Chars

  subroutine Banner_Chararray ( charArray, &
    & columnRange, alignment, skips, pattern )
    character(len=*), dimension(:), intent(in)  :: CHARARRAY
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: COLUMNRANGE
    character(len=1), intent(in), optional      :: ALIGNMENT ! L, R, C, or J
    integer, optional, intent(in)               :: Skips ! How many spaces between chars
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
    ! Here we begin printing
    ! Temporarily stop Stamping lines
    OldNeverStamp = stampOptions%neverStamp
    stampOptions%neverStamp = .true.
    ! Top border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
    do i = 1, size(chararray)
      ! Left star, then message, then right star
      call output( '*' )
      call AlignToFit( chararray(i), myColumnRange, myAlignment, skips )
      call BlanksToColumn( lineLen )
      call output( '*', advance = 'yes' )
    end do
    ! Bottom border
    call output( '*' )
    call blanks ( lineLen-2, FillChar=myFillChar )
    call output( '*', advance = 'yes' )
    ! Restore Stamping
    stampOptions%neverStamp = OldNeverStamp
  end subroutine Banner_Chararray

  ! -----------------------------------------------------  Beverbose_Chars  -----
  logical function Beverbose_Chars ( switch, threshold )
    ! Args
    character(len=*), intent(in) :: SWITCH
    integer, intent(in)          :: THRESHOLD
    ! Executable
    Beverbose_Chars = switchDetail( switches, switch ) > threshold .or. MLSVerbose
  end function Beverbose_Chars

  logical function Beverbose_CharArray ( switcharray, threshold )
    ! Args
    character(len=*), dimension(:), intent(in) :: switcharray
    integer, intent(in)                        :: THRESHOLD
    ! Internal variables
    integer                                    :: i
    ! Executable
    Beverbose_CharArray = MLSVerbose
    if ( Beverbose_CharArray ) return
    do i=1, size(switcharray)
      Beverbose_CharArray = Beverbose_CharArray .or. &
        & switchDetail( switches, switcharray(i) ) > threshold
    enddo
  end function Beverbose_CharArray

  ! -----------------------------------------------------  BlanksToColumn  -----
  subroutine BlanksToColumn ( column, fillchar, advance, dont_stamp )
  ! Output blanks to PRUNIT out to column COLUMN.
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
  end subroutine BlanksToColumn

  ! ------------------------------------------------  blanksToTab  -----
  ! Print blanks out to next tabstop
  ! (or else to tabstop number tabn)
  ! (or optionally that many copies of fillChar)
  subroutine blanksToTab ( tabn, fillChar, MayStayIfAtTab )
    ! Args
    integer, optional, intent(in)          :: TABN
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: MayStayIfAtTab
    ! Internal variables
    integer :: nColumn
    integer :: nTab
    logical :: MayStay
    ! Executable
    MayStay = .false.
    if ( present(MayStayIfAtTab) ) MayStay = MayStayIfAtTab
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
      if ( MayStay ) then
        nTab = findFirst( tabStops >= nColumn )
      else
        nTab = findFirst( tabStops > nColumn )
      endif
      if ( nTab > 0 ) &
        & call blanksToColumn( tabStops(nTab), fillChar )
    end if
  end subroutine blanksToTab

  ! ---------------------------------------------- DumpOutputOptions -----
  subroutine DumpOutputOptions( options )
    ! Show currently-active output options
    ! We put this here because it uses other routines from HighOutput
    ! Why aren't we using the CellDatabase?
    type(OutputOptions_T), intent(in) :: options
    ! Internal variables
    logical                      :: AlwaysWrap
    logical, parameter           :: checkingTabbing = .false.
    character(len=10), parameter :: decade = '1234567890'
    character(len=1), parameter  :: fillChar = '1' ! fill blanks with '. .'
    integer                      :: i
    integer                      :: IndentBy
    integer                      :: status
    integer                      :: WrapPastColumn
    ! Executable
    AlwaysWrap = options%AlwaysWrap
    WrapPastColumn = options%WrapPastColumn
    IndentBy = GetOutputStatus ( 'indent' )
    call blanks(80, fillChar='-', advance='yes')
    call HeadLine( 'Summary of output options', &
      & fillChar='-', before='*', after='*' )
    call outputNamedValue ( 'Unit number', options%prUnit, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    if ( options%prUnit < 0 ) then
      call outputNamedValue ( 'Meaning', prunitname(options ), &
        & advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    end if
    call outputNamedValue ( 'File name', trim(options%name), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Logging level', &
      & SeverityNamesFun(options%MLSMSG_Level), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Buffered?', options%buffered, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Skip MLSMSG logging?', options%SKIPMLSMSGLOGGING, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Log Parent Name?', options%logParent, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'AdvanceDefault', options%advanceDefault, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'SdFormatDefault', options%sdFormatDefault, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Indent By', IndentBy, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Always Wrap', AlwaysWrap, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    if ( AlwaysWrap ) &
  & call outputNamedValue ( 'Wrap past column', WrapPastColumn, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'Tab stops', tabstops, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    do i=1, MAXNUMTABSTOPS
      call tab( fillChar='.' )
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
      type(OutputOptions_T), intent(in) :: options
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

  ! ---------------------------------------------- DumpPatternOptions -----
  subroutine DumpPatternOptions( options )
    ! Show currently-active pattern options
    type(PatternOptions_T), intent(in) :: options
    ! Internal variables
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
    ! Executable
    call blanks(80, fillChar='-', advance='yes')
    call HeadLine( 'Summary of pattern options', &
      & fillChar='-', before='*', after='*' )
    call outputNamedValue ( 'use patterned blanks?', options%usePatternedBlanks, advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'special fills', trim(options%specialFillChars), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call outputNamedValue ( 'lineup fills', trim(options%lineupFillChars), advance='yes', &
      & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call newline
    call blanks(80, fillChar='-', advance='yes')
  end subroutine DumpPatternOptions

  ! ---------------------------------------------- DumpStampOptions -----
  subroutine DumpStampOptions( options )
    ! Show currently-active stamp options
    type(StampOptions_T), intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     call blanks(80, fillChar='-', advance='yes')
    call HeadLine( 'Summary of automatic stamp options', &
      & fillChar='-', before='*', after='*' )
     call outputNamedValue ( 'Never stamp', options%neverStamp, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Stamp end of line', options%post, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Show time', options%showTime, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Extra text', trim_safe(options%textCode), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Date format', trim_safe(options%dateFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Time format', trim_safe(options%timeFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Interval', options%interval, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Style of TimeStamps', trim_safe(options%TimeStampstyle), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call blanks(80, fillChar='-', advance='yes')
  end subroutine DumpStampOptions

  ! ---------------------------------------------- DumpTimeStampOptions -----
  subroutine DumpTimeStampOptions( options )
    ! Show output options
    type(TimeStampOptions_T), intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     call blanks(80, fillChar='-', advance='yes')
    call HeadLine( 'Summary of time stamp options', &
      & fillChar='-', before='*', after='*' )
     call outputNamedValue ( 'Stamp end of line', options%post, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Show date', options%showDate, advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Extra text', trim_safe(options%textCode), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Date format', trim_safe(options%dateFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Time format', trim_safe(options%timeFormat), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
     call outputNamedValue ( 'Style of TimeStamps', trim_safe(options%TimeStampstyle), advance='yes', &
       & fillChar=fillChar, before='* ', after='*', tabn=4, tabc=62, taba=80 )
    call blanks(80, fillChar='-', advance='yes')
  end subroutine DumpTimeStampOptions

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

  !-----------------------------------   FinalMemoryReport  -----
  subroutine FinalMemoryReport ( IsFinal )
    use MLSCommon, only: NoBlocksAllocated, NoBlocksDeAllocated, &
      & NoBytesAllocated, TotalAllocated, TotalDeAllocated
    use Memory_M, only: Memory_Used
    use Optional_M, only: Default
    use Output_M, only: Output
    ! Print Final or interim report summarizing allocates, deallocates
    logical, optional, intent(in)  :: IsFinal
    integer, save                  :: Previous_used = 0
    integer                        :: Total_used     ! Memory used in kilobytes (1024)

    ! Executable code
    call Output ( '* ', advance='no' )
    if ( .not. Default( IsFinal, .false. ) ) then
      call Output ( 'Interim', advance='no' )
    else
      call Output ( 'Final', advance='no' )
    endif
    call Output( ' report on allocates/deallocates *', advance='yes' )
    call OutputNamedValue ( 'Number of calls to _allocate_', &
      & NoBlocksAllocated )
    call OutputNamedValue ( 'Number of calls to _deallocate_', &
      & NoBlocksDeAllocated )
    call OutputNamedValue ( 'Total GB allocated', &
      & TotalAllocated*1.e-9 )
    call OutputNamedValue ( 'Total GB deallocated', &
      & TotalDeAllocated*1.e-9 )
    call OutputNamedValue ( 'Net GB remaining allocated', &
      & NoBytesAllocated*1.e-9 )
    call memory_used ( total=Total_used )
    call OutputNamedValue ( 'Compared with total memory used', &
      & Total_used*1.e-6 )
    call OutputNamedValue ( 'Change from last report', &
      & (Total_used-Previous_used)*1.e-6 )
    Previous_used = Total_used ! To remember for the next call


  end subroutine FinalMemoryReport

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

  ! -----------------------------------------------------  HeadLine  -----
  ! Print your message with extra formatting features; e.g.,
  ! *----------------  Your message here   ----------------*
  ! See also Banner, StyledOutput
  subroutine HeadLine ( inChars, fillChar, Before, After, &
    & ColumnRange, Alignment, Skips, Underline, Stretched )
    character(len=*), intent(in)                ::inChars
    character(len=1), intent(in), optional      :: fillChar      ! For padding
    character(len=*), intent(in), optional      :: Before, After ! text to print
    ! If columnRange(1) < 1, just use starting columns; otherwise move to
    integer, dimension(2), optional, intent(in) :: Columnrange
    character(len=1), intent(in), optional      :: Alignment ! L, R, C, or J
    integer, optional, intent(in)               :: Skips ! How many spaces between chars
    logical, optional, intent(in)               :: Underline
    logical, optional, intent(in)               :: Stretched   ! s t r e t c h
    !
    ! Internal variables
    character(len=2*len(inchars))    :: Chars
    character(len=1)      :: myAlignment
    integer, dimension(2) :: myColumnRange      ! To fit chars
    integer, dimension(2) :: myFullColumnRange  ! Must fit before and after, too
    integer               :: mySkips,  rightpadding
    character(len=1)      :: myFillChar
    character(len=8)      :: myBefore, myAfter
    logical                          :: myStretched
    logical               :: myUnderline
    character(len=len(chars))        :: underlineChars
    ! Executable
    Chars = ' '
    rightpadding = 0
    if ( present(columnRange) ) then
      myFullColumnRange = columnRange
    else
      myFullColumnRange = StyleOptions%HeadLineColumnRange ! (/ 1, 80 /)
    end if
    ! print *, 'myFullColumnRange ', myFullColumnRange 
    myColumnRange = myFullColumnRange
    mySkips = StyleOptions%HeadLineSkips ! 0
    if ( present(skips) ) mySkips = skips
    myFillChar = StyleOptions%HeadLinefill ! ' '
    if ( present(fillChar) ) myFillChar = fillChar
    myAlignment = StyleOptions%HeadLineAlignment ! 'C'
    if ( present(alignment) ) myAlignment = alignment
    myStretched = .false.
    if ( present(Stretched) ) myStretched = Stretched
    myUnderline = .false.
    if ( present(Underline) ) myUnderline = Underline
    if ( myStretched ) then
      Chars = Stretch ( inChars, 1, options='a' )
    else
      Chars = inChars
    endif
    myBefore = StyleOptions%HeadLineBefore
    myAfter = StyleOptions%HeadLineAfter
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
      call AlignToFit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn( myFullColumnRange(2)-rightpadding, &
        & fillChar=myFillChar, advance='no' )
      if ( len_trim(myAfter) > 0 ) call output( trim(myAfter), advance='no' )
    else
      call AlignToFit( chars, myColumnRange, myAlignment, skips )
      call blanksToColumn( myFullColumnRange(2)-rightpadding, advance='no' )
      if ( len_trim(myAfter) > 0 ) call output( trim(myAfter), advance='no' )
    end if
    call newLine
    if ( myUnderline .and. len_trim(chars) > 0 ) then
      ! The next trick replaces everything with a dash except spaces
      underlineChars = Replace( chars, ' ', '-', reverse=.true. )
      call AlignToFit( underlineChars, myColumnRange, myAlignment, skips )
      call newLine
    endif
  end subroutine HeadLine

  ! -----------------------------------------------------  Letsdebug  -----
  logical function Letsdebug_Chars ( switch, threshold )
    ! Args
    character(len=*), intent(in) :: SWITCH
    integer, intent(in)          :: THRESHOLD
    ! Executable
    Letsdebug_Chars = switchDetail( switches, switch ) > threshold .or. MLSDebug
  end function Letsdebug_Chars

  logical function Letsdebug_CharArray ( switcharray, threshold )
    ! Args
    character(len=*), dimension(:), intent(in) :: switcharray
    integer, intent(in)                        :: THRESHOLD
    ! Internal variables
    integer                                    :: i
    ! Executable
    Letsdebug_CharArray = MLSVerbose
    if ( Letsdebug_CharArray ) return
    do i=1, size(SwitchArray)
      Letsdebug_CharArray = Letsdebug_CharArray .or. &
        & switchDetail( switches, switcharray(i) ) > threshold
    enddo
  end function Letsdebug_CharArray

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

  function NumNeedsFormat_dcomplx( value, inFormat ) result ( format )
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
  function numToChars_double( value, format, trim ) result ( line )
    ! Args
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    logical, intent(in), optional :: Trim ! Trim blanks even if Format present
    logical :: Do_Trim
    character(len=30) :: line
    ! Internal variables
    character(len=30) :: FormatSpec
    integer :: I, J, K
    ! Executable
    FormatSpec = OutputOptions%sdFormatDefault
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

  function numToChars_single( value, format, trim ) result ( line )
    ! Args
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    logical, intent(in), optional :: Trim ! Trim blanks even if Format present
    logical :: Do_Trim
    character(len=30) :: line
    ! Internal variables
    character(len=30) :: FormatSpec
    integer :: I, J, K
    ! Executable
    FormatSpec = OutputOptions%sdFormatDefault
    if ( any( value == DPREFERDEFAULTFORMAT ) ) FormatSpec = '*'
    if ( present(Format)  ) then
      if ( format /= '*' ) FormatSpec = Format
    end if
    include 'numToChars.f9h'
  end function numToChars_single

  ! ---------------------------------------  OutputCalendar  -----
  subroutine OutputCalendar ( date, datenote, notes, dontWrap, moonPhases )
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
    logical :: PrintedDate
    integer :: row
    logical :: today ! Print extra '|' around today's date
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
    ! Here we begin printing
    ! Temporarily stop Stamping lines
    OldNeverStamp = stampOptions%neverStamp
    stampOptions%neverStamp = .true.
    call newline
    call AlignToFit( trim(monthName(month)), (/ 1, 100 /), 'c', skips=1 )
    call newline
    col2 = 0
    do wkdy=1, 7
      col1 = col2 + 1
      col2 = tabStops(wkdy)
      call AlignToFit( trim(daysOfWeek(wkdy)), (/ col1, col2 /), 'c' )
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
          PrintedDate = .false.
          col1 = col2 + 1
          col2 = tabStops(wkdy)
          today = ( days(iwk, wkdy) == day )
          if ( today ) then
            call output_('||')
            col2 = col2 - 1 ! Must scoot string one space to the lefts
          else
            call output_('|')
          end if
          if ( days(iwk, wkdy) < 1 ) then
            ! Don't write notes or anything else in "empty" days
          else if ( row == 1 ) then
            call writeIntsToChars( days(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call AlignToFit( trim(dateString), (/ col1, col2-1 /), 'r' )
            PrintedDate = .true.
          else if( row == numRows ) then
            call writeIntsToChars( daysOfYear(iwk, wkdy), dateString )
            dateString = adjustl(dateString)
            call AlignToFit( 'd' // trim(dateString), (/ col1, col2-1 /), 'r' )
            PrintedDate = .true.
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
            if ( .not. PrintedDate ) call Blanks( 1 )
            call output_('|')
          else
            call blanksToTab ( MayStayIfAtTab=.true. )
          end if
        end do ! wkdy
        call output_('|')
        call newline
      end do ! row
      ! begin with 
    end do ! week
    call blanksToTab( 7, fillChar='-' )
    call newline
    ! Restore tabstops, Stamping
    call settabs( '5-120+5' )
    stampOptions%neverStamp = OldNeverStamp
  end subroutine OutputCalendar

  ! ---------------------------------------  Output_Date_And_Time  -----
  subroutine Output_Date_And_Time ( date, time, &
    & from_where, msg, dateFormat, timeFormat, &
    & CPU_Seconds, wallClock_Seconds, advance )
    ! Output nicely-formatted date, time, and extra message
    ! We'll assume we won't want this line stamped with date and time
    ! (for fear of being redundant, which we fear)
    ! Optionally print CPU and wall clock usage, too.
    ! The first two args optionally override the defaults, which are to output
    ! both date and time
    logical, intent(in), optional          :: Date ! output date?
    logical, intent(in), optional          :: Time ! output time?
    ! Anything else to output?
    character(len=*), intent(in), optional :: From_where
    character(len=*), intent(in), optional :: Msg
    character(len=*), intent(in), optional :: Dateformat
    character(len=*), intent(in), optional :: Timeformat
    double precision, intent(in), optional :: Cpu_seconds
    integer, intent(in), optional          :: Wallclock_seconds
    character(len=*), intent(in), optional :: Advance
    ! Internal variables
    character(len=16) :: dateString
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    integer :: HH, MM, MS, SS
    logical :: myDate
    logical :: myTime
    character(len=3) :: MY_ADV
    real :: My_CPU, Seconds
    character(len=16) :: timeString
    ! Executable
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
  end subroutine Output_Date_And_Time

  ! ----------------------------------------------  OutputList  -----
  ! This family of routines outputs an array as a comma-separated list
  ! E.g., given the array (/ 1, 2, 3, .. /) outputs
  ! '(1, 2, 3, .. )'
  ! optionally using sep instead of ',' and delims instead of '()'
  subroutine OutputList_Chars ( array, sep, delims )
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
  end subroutine OutputList_Chars

  subroutine OutputList_Ints ( array, sep, delims )
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
  end subroutine OutputList_Ints

  ! ----------------------------------------------  OutputNamedValue  -----
  ! This family of routines outputs a paired name and value
  ! (Basically saving you a few lines over the idiom
  !  call output ( trim(name), advance='no' )
  !  call output ( ': ', advance='no' )
  !  call output ( value, advance='yes' )
  
  ! to print following line to stdout
  !  name: value
  ! Optional args control
  ! Before:     what extra to print at start of each line
  ! After:      what extra to print at end of each line
  ! Colon:      what to print instead of ':'
  ! FillChar:   instead of spaces if you use tabs to align name, value
  ! Tabn:       column number where name begins
  ! Tabc:       column number where colon occurs
  ! Taba:       column number where after begins
  ! Advance:    whether to advance after printing pair (by default we will advance)
  ! Dont_stamp: override setting to stamp end of each line
  ! By means of optional args you can create a line like the following
  ! *   name                   value   *
  ! See also startTable, AddRow, outputTable
  subroutine Output_Nvp_whatever ( name, &
   & chvalue, ivalue, cmvalue, dbvalue, snvalue, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)           :: name
    character(len=*), intent(in), optional :: chvalue
    complex, intent(in), optional          :: cmvalue
    double precision, intent(in), optional :: dbvalue
    integer, intent(in), optional          :: ivalue
    real, intent(in), optional             :: snvalue
    character(len=*), intent(in), optional :: Advance
    character(len=1), intent(in), optional :: Colon
    character(len=1), intent(in), optional :: Fillchar
    integer, intent(in), optional          :: Tabn
    integer, intent(in), optional          :: Tabc
    integer, intent(in), optional          :: Taba
    logical, intent(in), optional          :: Dont_stamp
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=*), intent(in), optional :: options
    ! Local variables
    if ( present(chvalue) ) then
      call Output_Nvp_character ( name, chvalue, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    elseif ( present(cmvalue) ) then
      call Output_Nvp_complex ( name, cmvalue, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    elseif ( present(dbvalue) ) then
      call Output_Nvp_double ( name, dbvalue, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    elseif ( present(ivalue) ) then
      call Output_Nvp_integer ( name, ivalue, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    elseif ( present(snvalue) ) then
      call Output_Nvp_single ( name, snvalue, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    endif
  end subroutine Output_Nvp_whatever

  subroutine Output_Nvp_character ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
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
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    else
      call possiblyTrimmedvalue ( name, value, &
        & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    endif
  contains
    subroutine possiblyTrimmedvalue ( name, value, &
      & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
      character(len=*), intent(in)          :: name
      character(len=*), intent(in)          :: value
      include 'output_name_value_pair.f9h'
    end subroutine possiblyTrimmedvalue
  end subroutine Output_Nvp_character

  subroutine Output_Nvp_complex ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    complex, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_complex

  subroutine Output_Nvp_double ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    double precision, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_double

  subroutine Output_Nvp_dbl_array ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    double precision, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_dbl_array

  subroutine Output_Nvp_int_array ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    integer, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_int_array

  subroutine Output_Nvp_integer ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_integer

  subroutine Output_Nvp_log_array ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    logical, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_log_array

  subroutine Output_Nvp_logical ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    logical, intent(in)                   :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_logical

  subroutine Output_Nvp_single ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_single

  subroutine Output_Nvp_sngl_array ( name, value, &
   & Advance, colon, fillChar, Before, After, Tabn, Tabc, Taba, Dont_stamp, options )
    character(len=*), intent(in)          :: name
    real, dimension(:), intent(in)     :: value
    include 'output_name_value_pair.f9h'
  end subroutine Output_Nvp_sngl_array

  ! ----------------------------------------------  OutputParagraph  -----
  ! Outputs text formatted as a paragraph where each line begins on
  ! columnrange(1) and ends on columnrange(2).
  ! Optionally, 
  ! (1)  you may specify an alignment, one of L R C or J
  ! (2)  you may specify an ident
  !      (a) of the first line, if indent > 0
  !      (b) of every line except the first
  !          if indent < 0 and columnrange(1) > 1
  subroutine OutputParagraph ( text, columnrange, alignment, indent )
    ! Args
    character(len=*), intent(in)                   :: text ! text to print
    integer, dimension(2), intent(in)              :: columnrange ! left, right columns
    character(len=1), optional, intent(in)         :: alignment ! L R C or J
    integer, optional, intent(in)                  :: indent ! num of spaces
    ! Internal variables
    integer, parameter                             :: maxindent = 15
    integer, parameter                             :: maxLines = 256
    integer, parameter                             :: MaxStrElementLength = LineLen
    character (len=MaxStrElementLength), dimension(:), allocatable    &
      &                                            :: array
    logical, parameter                             :: countEmpty = .true.
    integer                                        :: i
    logical                                        :: ignoreFirstIndent
    character(len=maxindent+len(text))             :: intext ! indented text
    character                                      :: firstAlignment
    character                                      :: myAlignment
    integer                                        :: nElems
    character, parameter                           :: null = achar(0)
    integer                                        :: status
    character(len=maxindent+maxLines+len(text))    :: wrappedtext ! indented text
    ! Executable
    ! Reasonable choices of args?
    if ( len_trim (text) < 1 ) return
    if ( columnrange(2) <= columnrange(1) ) return
    myAlignment = 'l'
    if ( present(Alignment) ) myAlignment = Alignment
    ignoreFirstIndent = .false.
    FirstAlignment = myAlignment
    if ( present(indent) ) then
      if ( indent > 0 ) then
        intext = repeat( ' ', indent ) // adjustl(text)
        if ( Lowercase(myAlignment) == 'j' ) firstAlignment = 'r'
      elseif ( indent < 0 ) then
        ignoreFirstIndent = .true.
        intext = adjustl(text)
      else
        intext = text
      endif
    else
      intext = text
    endif
    ! print *, trim(inText)
    call wrap( intext, wrappedText, columnrange(2)-columnrange(1)+1, &
      & inseparator=null, dontSqueeze=.true. )
    ! print *, trim(wrappedText)
    ! call output( trim(wrappedText), advance='yes' )
    nElems = NumStringElements( wrappedText, countEmpty, null )
    ! print *, 'nElems: ', nElems
    allocate ( Array(nElems), STAT=status )
    ! call List2Array( trim(wrappedText), array, countEmpty, null )
    call wrap( trim(wrappedText), array, columnrange(2)-columnrange(1)+1, &
      & inseparator=null )
    ! print *, 1, ' ', trim(array(1))
    ! print *, 2, ' ', trim(array(2))
    ! print *, NElems, ' ', trim(array(NElems))
    do i=1, NElems
    ! Convert every null into a space
      array(i) = Replace( array(i), null, ' ' )
      if ( i == 1 .and. ignoreFirstIndent ) then
        call AlignToFit( array(i), (/ 1, ColumnRange(2) /) , Alignment )
      elseif ( i == NElems .and. Lowercase(myAlignment) == 'j' ) then
        ! We won't attempt to left-right justify the last line in the paragraph
        call AlignToFit( array(i), ColumnRange, Alignment='l' )
      elseif ( i == 1 ) then
        ! We won't attempt to left-right justify the first line in the paragraph
        ! if it's indented
        call AlignToFit( array(i), ColumnRange, FirstAlignment )
      else
        call AlignToFit( array(i), ColumnRange, Alignment )
      endif
      call newLine
    enddo
  end subroutine OutputParagraph

  ! ----------------------------------------------  OutputTable  -----
  ! Outputs a 2d character array as a table
  ! Optionally, 
  ! (1)  you may supply the array; otherwise, cellDatabase will be used
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
  !      by a wall of special HeadLiner characters
  subroutine outputTable ( array, sep, border, cellWidth, &
    & interior, HeadLiner, alignment )
    ! Args
    character(len=*), dimension(:,:), optional, intent(in)   &
      &                                            :: array
    character(len=1), optional, intent(in)         :: sep       ! between cols
    character(len=1), optional, intent(in)         :: border    ! outside
    integer, optional, intent(in)                  :: cellWidth
    character(len=1), optional, intent(in)         :: interior  ! between rows
    character(len=1), optional, intent(in)         :: HeadLiner ! 1st row are headers
    character(len=1), optional, intent(in)         :: alignment ! L, R, or C
    ! Local variables
    integer                                        :: status
    ! Executable
    if ( present( array ) ) then
      call outputTableArray ( array, sep, border, cellWidth, &
        & interior, HeadLiner, alignment )
    elseif ( .not. associated ( cellDatabase ) ) then
      call banner ( 'Empty table' )
    else
      call outputTableArray ( cellDatabase, sep, border, cellWidth, &
        & interior, HeadLiner, alignment )
      deallocate( cellDatabase, stat=status )
      nullify( cellDatabase )
    endif
  end subroutine outputTable

  ! ----------------------------------------------  ResetTabs  -----
  ! Restore tab stops to what was in effect at start
  ! Optionally returning them as an integer array
  subroutine ResetTabs ( tabs )
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
  end subroutine ResetTabs

  ! ----------------------------------------------  RestoreSettings  -----
  ! Restore output settings by call to RestoreOutputSettings
  ! Optionally restore StyleOptions and Tabs
  subroutine RestoreSettings ( settings )
    ! Args
    character(len=*), optional, intent(in) :: settings
    ! Local variables
    character(len=*), parameter            :: allSettings = 'style'
    character(len=64)                      :: mySettings 
    ! Executable
    call RestoreOutputSettings ( settings )
    mySettings = ' '
    if ( present(settings) ) mySettings = settings
    if ( index(mySettings, '*') > 0 ) mySettings = allSettings
    mySettings = lowercase(mySettings)
    if ( index(mySettings, 'tabs') > 0 ) call ResetTabs
    if ( index(mySettings, 'style') > 0 ) StyleOptions = DefaultStyleOptions
  end subroutine RestoreSettings 

  ! ----------------------------------------------  SetStamp  -----
  subroutine SetStamp ( textCode, showTime, dateFormat, timeFormat, &
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
  end subroutine SetStamp

  ! ----------------------------------------------  setTabs  -----
  subroutine setTabs ( range, tabs )
    ! Set tabstops
    ! Methods:
    ! (1) a string range; e.g., "8, 32-100+8"
    !     is converted to its expanded form 
    !      "8, 32, 40, 48, 56, 64, 72, 80, 88, 96"
    ! (2) an array of ints; e.g., (/ 4, 9, 12, 18, 22, 30, 35, 40 /)
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
  !      by a wall of special HeadLiner characters
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
  ! e.g., "--Banner" or "-B"
  ! See also banner
  !
  ! options contains       style
  !    ----                -----
  !      B               Banner
  !      H               Headline
  !      U               Underline
  !      S               s t r e t c h e d
  !            sub-options of Banner
  !      L, C, R, or J   Alignment (if not default)
  !      w               Don't print top stripe of banner
  !      m               Don't print bottom stripe of banner
  !
  ! Some, though not all, of these options may be combined; e.g. "-USB"
  !
  ! We could dig more deeply into options to allow
  ! it to pass Alignment, Fill, After, etc.
  subroutine StyledOutput ( chars, options )
    character(len=*), intent(in)                :: chars
    character(len=*), intent(in), optional      :: options
    ! Internal variables
    character(len=1)                            :: Alignment
    logical                                     :: AsBanner
    logical                                     :: AsHeadLine
    logical                                     :: Bottomless
    logical                                     :: Stretched
    logical                                     :: Topless
    logical                                     :: Underline
    character(len=4*len(chars))                 :: UnderlineChars
    ! Executable
    asBanner   = .false.
    asHeadLine = .false.
    if ( .not. present(options ) ) then
      call output( chars, advance='yes' )
      return
    endif
    asBanner   = index( options, 'B' ) > 0
    asHeadLine = index( options, 'H' ) > 0
    Stretched  = index( options, 'S' ) > 0
    Underline  = index( options, 'U' ) > 0 .or. StyleOptions%Underline
    ! Sub-options for Banner
    Topless  = index( options, 'w' ) > 0
    Bottomless  = index( options, 'm' ) > 0
    if ( asBanner ) then
      Alignment = StyleOptions%BannerAlignment
      if ( index(options, 'L' ) > 0 ) then
        Alignment = 'L'
      elseif ( index(options, 'R' ) > 0 ) then
        Alignment = 'R'
      elseif ( index(options, 'J' ) > 0 ) then
        Alignment = 'J'
      else
        Alignment = 'C'
      endif
      Stretched = Stretched .or. StyleOptions%BannerStretch
      call Banner( chars, Underline=Underline, Stretched=Stretched, &
        & Alignment=Alignment, &
        & Topless=Topless, BottomLess=BottomLess )
    elseif ( asHeadLine ) then
      StyleOptions%HeadLineFill = '-'
      StyleOptions%HeadLineBefore = '*'
      StyleOptions%HeadLineAfter = '*'
      Stretched = Stretched .or. StyleOptions%HeadlineStretch
      call HeadLine( chars, Underline=Underline, Stretched=Stretched ) 
      StyleOptions = DefaultStyleOptions
    elseif ( Stretched .and. Underline ) then
      underlineChars = Stretch( chars, 1, options='a' )
      call output( underlineChars, advance='yes' )
      underlineChars = Replace( underlineChars, ' ', '-', reverse=.true. )
      call output( underlineChars, advance='yes' )
    elseif ( Underline ) then
      call output( chars, advance='yes' )
      ! The next trick replaces everything with a dash except spaces
      underlineChars = Replace( chars, ' ', '-', reverse=.true. )
      call output( underlineChars, advance='yes' )
    elseif ( Stretched ) then
      underlineChars = Stretch( chars, 1, options='a' )
      call output( underlineChars, advance='yes' )
    else
      call output( chars, advance='yes' )
    endif
  end subroutine StyledOutput

  ! ------------------------------------------------  TimeStamp  -----
  ! time-stamp output on demand, not automatic:
  ! Either in style pre or post
  ! (pre) '(HH:MM:SS) chars'
  ! (post) 'chars (HH:MM:SS)'
  ! Note that in pre-style, the time will be printed only if getOutputStatus( 'start' ) == 1 true
  ! in post-style, the time will be printed only if MY_ADV is 'yes'
  subroutine TimeStamp_char ( chars, &
    & advance, from_where, dont_log, log_chars, insteadofblank, style, date )
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
    my_style = TimeStampOptions%TimeStampstyle
    if ( present(style) ) my_style = lowercase(style)
    myDate = TimeStampOptions%showDate
    if ( present(date) ) myDate = date
    if ( my_style == 'post' ) then
      call output_( CHARS, &
        & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
        & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK, DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_(' (', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call Output_Date_And_Time( date=myDate, &
          & dateFormat=TimeStampOptions%dateFormat, &
          & timeFormat=TimeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='yes', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      end if
    else
      if ( getOutputStatus( 'start' ) == 1 ) then
        call output_('(', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call Output_Date_And_Time( date=myDate, &
          & dateFormat=TimeStampOptions%dateFormat, &
          & timeFormat=TimeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      end if
      call output_( chars, &
        & advance, from_where, dont_log, &
        & log_chars, insteadofblank, dont_stamp=dont_stamp )
    end if
  end subroutine TimeStamp_char

  subroutine TimeStamp_integer ( int, &
    & places, advance, fill, format, Before, After, style, date )
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
    my_style = TimeStampOptions%TimeStampstyle
    if ( present(style) ) my_style = lowercase(style)
    myDate = TimeStampOptions%showDate
    if ( present(date) ) myDate = date
    if ( my_style == 'post' ) then
      call output( INT, PLACES, &
        & ADVANCE='no', FILL=FILL, FORMAT=FORMAT, BEFORE=BEFORE, AFTER=AFTER, &
        & DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_(' (', ADVANCE='no', DONT_STAMP=DONT_STAMP )
        call Output_Date_And_Time( date=myDate, &
          & dateFormat=TimeStampOptions%dateFormat, &
          & timeFormat=TimeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='yes', DONT_STAMP=DONT_STAMP)
      end if
    else
      if ( getOutputStatus( 'start' ) == 1 ) then
        call output_('(', ADVANCE='no', DONT_STAMP=DONT_STAMP)
        call Output_Date_And_Time( date=myDate, &
          & dateFormat=TimeStampOptions%dateFormat, &
          & timeFormat=TimeStampOptions%timeFormat, &
          & advance='no')
        call output_(')', ADVANCE='no', DONT_STAMP=DONT_STAMP)
      end if
      call output( int, places, &
        & advance, fill, format, before, after, dont_stamp=dont_stamp )
    end if
  end subroutine TimeStamp_integer

  subroutine TimeStamp_logical ( value, &
    & advance, from_where, dont_log, log_chars, insteadofblank, style, date )
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
    call TimeStamp_char(str, &
    & advance, from_where, dont_log, log_chars, insteadofblank, style, date )
  end subroutine TimeStamp_logical

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
    allocate ( tempDatabase(newSize,2), stat=status )
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
  !      by a wall of special HeadLiner characters
  ! (6)  If any row begins with the special formatting sequence '<<', then
  !      one of the following occurs:
  !      special format             action
  !          <<l>>xxxx       merge:  left-aligned xxx stretched across table
  !          <<c>>xxxx       merge:  centered xxx stretched across table
  !          <<r>>xxxx       merge:  right-aligned xxx stretched across table
  !          <<d>>x          divider: inserts a wall of x's
  subroutine outputTableArray ( array, sep, border, cellWidth, &
    & interior, HeadLiner, alignment )
    ! use Dump_0, only: Dump
    ! Args
    character(len=*), dimension(:,:),intent(in)    :: array
    character(len=1), optional, intent(in)         :: sep       ! between cols
    character(len=1), optional, intent(in)         :: border    ! outside
    integer, optional, intent(in)                  :: cellWidth
    character(len=1), optional, intent(in)         :: interior  ! between rows
    character(len=1), optional, intent(in)         :: HeadLiner ! 1st row are headers
    character(len=1), optional, intent(in)         :: alignment ! L, R, or C
    ! Local variables
    character(len=1)                               :: align
    integer                                        :: i
    integer                                        :: j
    integer                                        :: k
    integer                                        :: left
    integer                                        :: minWidth
    character(len=1)                               :: myAlignment
    character(len=1)                               :: myBorder
    character(len=1)                               :: mySep
    character(len=1)                               :: myHeadLiner
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
    
    myHeadLiner = myInterior
    if ( present(HeadLiner) ) myHeadLiner = HeadLiner
    
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
    ! Here we begin printing
    ! Temporarily stop Stamping lines
    OldNeverStamp = stampOptions%neverStamp
    stampOptions%neverStamp = .true.
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
          call AlignToFit ( trim(array(i,1)(k+5:)), &
            & (/ right, totalWidth /), 'c' )
        case ('r') ! right-aligned merge
          call AlignToFit ( trim(array(i,1)(k+5:)), &
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
        call AlignToFit ( trim(array(i,j)), (/ left, right /), myAlignment )
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
      ! Interior cell walls or HeadLiners
      if ( len_trim(myHeadLiner) > 0 .and. i == 1 .and. &
        & myHeadLiner /= achar(0) ) then
        call output( repeat( myHeadLiner, totalWidth ), advance='yes' )
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
    ! Restore Stamping
    stampOptions%neverStamp = OldNeverStamp
  end subroutine outputTableArray

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

end module HighOutput

! $Log$
! Revision 2.42  2021/09/09 23:23:15  pwagner
! Show memory_used and its change since last memory report
!
! Revision 2.41  2021/07/22 23:20:14  pwagner
! Repair printing errors in Calendar
!
! Revision 2.40  2020/06/30 23:18:07  pwagner
! Improved comments among toc and api blocs
!
! Revision 2.39  2020/06/09 21:56:32  pwagner
! report on allocates/deallocates now stands out with stars
!
! Revision 2.38  2020/04/30 23:09:54  pwagner
! Added optional arg IsFinal to FinalMemoryReport
!
! Revision 2.37  2020/04/27 21:32:07  pwagner
! Added procedure to print FinalMemoryReport
!
! Revision 2.36  2019/11/11 23:08:27  pwagner
! Added OutputParagraph
!
! Revision 2.35  2019/10/30 20:07:18  pwagner
! Banner and styledOutput may align output as one of {LRCJ}
!
! Revision 2.34  2019/10/01 23:40:51  vsnyder
! Add Trim optional argument to floating-point output
!
! Revision 2.33  2019/08/01 23:44:25  pwagner
! Removed unused stuff; numerous other changes
!
! Revision 2.32  2019/07/17 20:19:21  pwagner
! StyledOutput can now underline and/or stretch, too
!
! Revision 2.31  2019/07/09 23:52:17  pwagner
! The values added by AddRow may now span several lines if needed
!
! Revision 2.30  2019/05/15 23:20:43  pwagner
! Non-essential housekeeping
!
! Revision 2.29  2019/02/21 22:35:26  pwagner
! Improved DumpOutputOptions
!
! Revision 2.28  2019/01/24 18:33:20  pwagner
! Reorganized modules that print to simplify toolkit-free builds
!
! Revision 2.27  2018/05/11 20:33:36  pwagner
! Make sure myChars long enough; reinitialize outputLines after use
!
! Revision 2.26  2018/04/19 23:42:20  pwagner
! Switch arg may be array in BeVerbose and LetsDebug
!
! Revision 2.25  2018/01/05 01:21:16  pwagner
! Corrected data type for OldNeverStamp
!
! Revision 2.24  2018/01/03 01:13:51  pwagner
! Prevent time stamps from interrupting tables, banners
!
! Revision 2.23  2017/12/22 00:25:24  pwagner
! Add move some items from DumpOuputOptions to new DumpPatternOptions
!
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
! May specify format in call to AddRow
!
! Revision 2.11  2016/03/25 00:37:06  pwagner
! Added OUTPUTANYNAMEDVALUE
!
! Revision 2.10  2015/09/24 18:50:41  pwagner
! May choose different pattern for stripes in banner
!
! Revision 2.9  2015/05/18 17:42:50  pwagner
! AddRow and startTable maintains an internal Table for outputTable to output
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
