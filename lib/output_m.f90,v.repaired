! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Output_M

  ! Normal level printing and formatting
  ! Directs output to either stddout or PrUnit
  ! See also dump_0 and Printit_m
  ! For higher-level procedures, see HighOutput
  
  ! Overview:
  ! Text, scalars, and arrays of numbers may be output
  ! in a controlled way according a number of settings
  ! contained in ..Options  datatypes or by arguments
  ! passed during the call
  ! Major features include whether or when to
  !   (*) direct the output to stdout or elsewhere @
  !   (*) time stamp each line #
  !   (*) indent each line
  !   (*) apply a special format
  !   (*) advance to the next line after printing
  !   (*) defer and store up output in a temporary buffer
  !   (*) pause, exit, or crash hard after printing
  !   (*) go 'silent' until commanded to resume printing
  !   (@) If we don't print to stdout, we may use the MLSMessaging facility
  !       or we may do both
  !   (#) Stamping may be applied automatically or on command

  use Machine, only: Crash_Burn, Exit_With_Status, NeverCrash
  use MLSCommon, only: Filenamelen, Finite_Signal, &
    & Is_What_Ieee
  use MLSStrings_0, only: NCharsInFormat, ReplaceNonAscii, Lowercase, &
    & Readintsfromchars, Stretch, Trim_Safe
  use PrintIt_M, only: Assemblefullline, Get_Config, &
    & MLSMSG_Crash, MLSMSG_Debug, MLSMSG_Info, MLSMSG_Error, &
    & MLSMSG_Severity_To_Quit, &
    & InvalidPrUnit          => InvalidLogUnit, &  
    & StdoutPrUnit           => StdoutLogUnit, &   
    & MsgLogPrUnit           => DefaultLogUnit, &  
    & BothPrUnit             => BothLogUnit, &     
    & OutputLinesPrUnit      => BufferedLogUnit, & 
    & PrUnitName             => LogUnitName, &
    & PrintItOut, MLSMessageConfig
  use IO_Stuff, only: Pause

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! OutputLines              If PrUnit = OUTPUTLINESPRUNIT, 
!                            this is where we store the output until flushed
! OutputOptions            where to send output and how to format it
! (some components)
! MLSMSG_Level             MLSMessage level if so logged
! PrUnit                   How to direct output:
!                          OUTPUTLINESPRUNIT :: held in outputLines
!                          MSGLOGPRUNIT      :: logged via MLSMessage
!                          STDOUTPRUNIT      :: print stdout
!                          BOTHPRUNIT        :: both stdout and logged
!                          > 0 :: print to Fortran unit number PrUnit
! SkipMLSMSGLogging        whether to skip MLSMessage by default
! StampOptions             whether and how to stamp each output automatically
! TimeStampOptions         how to stamp when calling timeStamp

!     (subroutines and functions)
! Advance_is_yes_or_no     parse the advance=.. field
!                            to say whether to advance to the next line
!                            and optionally to set advanced formatting
! AddToIndent              add to the number of blanks indented; subtract if < 0
! Beep                     print message to error_unit
! Blanks                   print specified number of blanks [or fill chars]
! FlushOutputLines         print the current outputLines; then reset to ''
! FlushStdout              flush any buffered lines to stdout
! GetOutputStatus          returns normally private data
! IsOutputSuspended        returns TRUE if output is suspended; i.e. running silent
! Newline                  print a newline
! Output                   print argument
! PrintOutputStatus        print normally private settings and options
! PrUnitName               Where output is directed
! ResetIndent              set indenting back to 0
! RestoreSettings          restore default settings for output, stamps, tabs
! RevertOutput             revert output to file used before switchOutput
!                           if you will revert, keepOldUnitOpen when switching
! ResumeOutput             resume suspended output
! SetFillPattern           set a special Fill pattern 
!                           (that can be used in call to blanks, banner, etc.)
! SetTruthPattern          set the chars that will be printed in place of 'T' 'F'
! SetOutputStatus          sets normally private data
! SetAdvancedOption        sets advanced or miscellaneous option values
! SuspendOutput            suspend output; run silent
! SwitchOutput             switch output to a new named file, 
!                            or else to 'stdout'
! === (end of toc) ===

! === (start of api) ===
! char* Advance_is_yes_or_no ( [char* str] )
! Beep ( [char* chars] )
! Blanks ( int n_blanks, [char fillChar], [char* advance] )
! FlushStdout
! FlushOutputLines ( [int prUnit] )
! int GetOutputStatus( char* name, [char* value] )
! log IsOutputSuspended ()
! NewLine ( [log dont_make_blank_line] )
! Output ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [log dont_stamp], [int newlineval] )
! Output ( char* chars(:), [char* advance],
!          [char* insteadOfBlank], [int newlineval] )
! Output ( value, [char* format], [char* advance],
!          [char* Before], [char* After] )
!       where value can be any numerical type, either scalar or 1-d array
!       advance may be
!       'y[es]', 't[rue]'  advance to next line after printing
!       'n[o]', 'f[alse]'  don't advance to next line after printing
!       'stderr'  print to stderr instead of default
! PrintOutputStatus ( [char* keywords] )
! ResumeOutput
! RevertOutput
! RestoreSettings ( [char* settings] )
! SetFillPattern ( pattern, [fillChar] )
! SetTruthPattern ( char* TrueFalse[2] )
! SetAdvancedOption ( [char* str] )
! SetOutputStatus( char* name, int value )
! SuspendOutput
! SwitchOutput ( char* filename, [int unit] )
! === (end of api) ===
!
! Notes:
! (1) By calling appropriate functions and procedures you can adjust aspects of
! how and where we output, and others can be changed by setting various
! public global parameters directly
! (in OO-speak they are class-level rather than instance-level)
! (2) Sometimes there is more than one way to accomplish the same thing
! E.g., calling timeStamp or using setStamp before calling output
! (Should we remove one of these?)
! (3) Another example is the advance argument which may be composed of multiple
! sub-arguments, each to set an advanced feature
! (see advancedOptions below)
! (4) To understand the codes for dateformat and timeFormat, see dates_module
! 
  ! Where to output?
  ! These apply if we don't output to a fortran unit number (which is > 0)
  ! See also Beep command, and advance='stderr' or advance='arg1 ..'
  ! Possible values for PrUnit
  public :: InvalidPrUnit    
  public :: StdoutPrUnit     
  public :: MsglogPrUnit     
  public :: BothPrUnit       
  public :: OutputlinesPrUnit

  ! Which of the above would involve stdout
  logical, parameter, private :: UseStdout(OutputlinesPrUnit:InvalidPrUnit) = &
    !    OUTPUTLINESPRUNIT  BOTHPRUNIT  MSGLOGPRUNIT  STDOUTPRUNIT
    & (/ .false.,           .true.,     .false.,      .true., &
    !    INVALIDPRUNIT
    &    .false. /)

  integer, save, private :: OLDUNIT = -1 ! Previous PrUnit for output.
  logical, save, private :: OLDUNITSTILLOPEN = .TRUE.

  public :: AddToIndent, Advance_Is_Yes_Or_No, Beep, Blanks, &
    & FlushOutputLines, FlushStdout, GetOutputStatus, IsOutputSuspended, &
    & Newline, Output, Output_Char_NoCR, PrintOutputStatus, PrUnitName, &
    & ResetIndent, RestoreSettings, ResumeOutput, RevertOutput, &
    & SetFillPattern, SetAdvancedOption, SetOutputStatus, SetTruthPattern, &
    & SuspendOutput, SwitchOutput

  ! These types made public because the class instances are public
  public :: OutputOptions_T
  public :: PatternOptions_T
  public :: StampOptions_T
  public :: TimeStampOptions_T

  ! We can use the advance=.. mechanism to convey extra formatting info
  interface getOption
    module procedure getOption_char, getOption_log
  end interface

  ! Embeddded <cr> print multiple lines
  interface Output
    module procedure Output_Char, Output_Char_array, output_complex
    module procedure Output_Dcomplex, output_double
    module procedure Output_Integer, output_integer_array
    module procedure Output_Logical, output_logical_array
    module procedure Output_Single, output_double_array, output_single_array
    module procedure Output_String
  end interface

  ! This won't filter for <cr>
  interface Output_
    module procedure Output_Char_NoCR
  end interface

  ! We can use the OutputLines mechanism for user-controlled
  ! buffering, filtering, grep-ing, or whatever
  integer, parameter :: numPatterns = 13 ! How many special patterns may we store
  integer, public, parameter :: MaxOutputlinesLen = 2048 ! How many chars it can hold
  character(len=MAXOUTPUTLINESLEN), public, save     :: OutputLines = ' '

  ! This is the type for configuring how to automatically format
  ! lines and whether they should be sent to stdout or elsewhere
  ! The default values work well; be cautious in overriding them
  type OutputOptions_T
    integer :: PrUnit = StdoutPrUnit    ! Unit for output (see comments above).  
    integer :: MLSMSG_Level        = MLSMSG_Info ! What level if logging
    integer :: NewLineVal          = 10 ! 13 means <cr> becomes new line; -999 means ignore
    integer :: NArrayElmntsPerLine = 7
    integer :: NBlanksBtwnElmnts   = 3
    integer :: WrapPastColumn      = 0
    logical :: AlwaysWrap          = .false. ! Lines longer than WrapPastColumn
    logical :: Buffered            = .true.
    logical :: LogParent           = .false. ! Show who called output, not output
    logical :: PrUnitLiteral       = .false. ! output to prUnit even if < 0
    logical :: SkipMLSMsgLogging   = .false.
    character(len=FileNameLen) :: name = 'stdout'
    character(len=3)  :: AdvanceDefault = 'no' ! if advance=.. missing
    character(len=12) :: SdFormatDefault = '*' ! * means default format spec
    character(len=1)  :: ArrayElmntSeparator = ' '
    character(len=27) :: ParentName = "$RCSfile$"
  end type

  ! This is the type for advanced formatting options from the optional arg
  !    advance=..
  ! Note:
  ! Pausing to wait for user input may hang that job running at the sips
  ! Remember to strip out any debugging use of 
  !    advance='.. pause ..'
  ! before delivery
  type AdvancedOptions_T
    logical :: stretch             = .false. ! p r i n t  l i k e  t h i s   ?
    logical :: bannered            = .false. ! print as a banner
    logical :: headered            = .false. ! print as a headline
    logical :: pause               = .false. ! pause and wait for user input
    type(OutputOptions_T) :: originalOptions
  end type

  ! This is the type for showing special patterns in blanks or other occasions
  ! where something interesting and eye-catching is desired.
  ! These predefined settings should normally suit us well; don't change them
  ! without testing!
  type patternOptions_T
    logical :: usePatternedBlanks  = .true. ! Use patterns for special fillChars
    character(len=numPatterns) :: specialFillChars = '0123456789ABC'
    character(len=numPatterns) :: lineupFillChars =  'yynnnnynnnyyn' ! whether they line up
    character(len=16), dimension(numPatterns) :: patterns = (/ & ! on consecutive lines
                                            &  '                ' , &
                                            &  '(. )            ' , &
                                            &  '(. .)           ' , &
                                            &  '(.  .)          ' , &
                                            &  '(.   .)         ' , &
                                            &  '(.. ..)         ' , &
                                            &  '(- )            ' , &
                                            &  '(- -)           ' , &
                                            &  '(-  -)          ' , &
                                            &  '(- .. )         ' , &
                                            &  '(= )            ' , &
                                            &  '(~ )            ' , &
                                            &  '(= ~ )          ' /)
                                            !   12345678901234567890
    ! Here are samples of the special pre-defined patterns above
    !       Special patterns used in blanks
    ! pattern
    !  0                                                                 
    !  1 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   
    !  2  . .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .   
    !  3  .  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  ..  .   
    !  4  .   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   ..   .   
    !  5  .. .... .... .... .... .... .... .... .... .... .... .... ..   
    !  6 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    !  7  - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -   
    !  8  -  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -   
    !  9  - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - .. - ..    
    !  A = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
    !  B ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   
    !  C  = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~ = ~    
  end type

  type(AdvancedOptions_T), public, save :: AdvancedOptions
  type(AdvancedOptions_T), private, save :: DefaultAdvancedOptions
  type(OutputOptions_T), public, save   :: OutputOptions
  type(PatternOptions_T), public, save  :: PatternOptions
  type(OutputOptions_T), private, save  :: DefaultOutputOptions

  ! This is the type for configuring whether and how to automatically stamp
  ! lines sent to stdout
  ! If interval > 1 only a fraction of lines will be stamped
  ! If interval > 10, the stamps will not be in-line, instead
  ! they will appear alone as page headers
  ! (As an alternative, use timeStamp to stamp only individual lines)
  type stampOptions_T
    logical :: neverStamp = .true.  ! if true, forget about automatic stamping
    logical :: post       = .true.      ! Put stamp at end of line?
    logical :: showTime   = .false. ! Don't show date or time unless TRUE
    character(len=24) :: textCode = ' '
    ! Don't show date unless dateFormat is non-blank
    character(len=16) :: dateFormat = ' '
    character(len=16) :: timeFormat = 'hh:mm'
    integer :: interval = 1 ! 1 means stamp every line; n means every nth line
    character(len=8) :: timeStampStyle = 'post' ! 'pre' or 'post'
  end type

  type(stampOptions_T), public, save :: stampOptions ! Could leave this private
  type(stampOptions_T), private, save :: DefaultStampOptions

  ! This is the type for configuring how the timeStamp stamps its lines)
  type timeStampOptions_T
    logical           :: post = .true.      ! Put stamp at end of line?
    logical           :: showDate = .false. ! Don't show date unless TRUE
    character(len=24) :: textCode = ' '
    ! Don't show date unless dateFormat is non-blank
    character(len=16) :: DateFormat     = 'yyyy-mm-dd'
    character(len=16) :: TimeFormat     = 'hh:mm:ss'
    character(len=8)  :: TimeStampStyle = 'post' ! 'pre' or 'post'
  end type

  type(timeStampOptions_T), public, save :: TimeStampOptions ! Could leave this private
  type(timeStampOptions_T), private, save :: DefaultTimeStampOptions
  logical, public, save                   :: MustRestoreAdvOpts = .false.

  ! Private parameters
  logical :: alreadyLogged ! Would we need to print again?
  logical, save :: DeeBug = .false.
  character(len=2), parameter :: defaultNewLineCode = achar(0) // 'n' ! not '%n'
  logical, save, private :: SWITCHTOSTDOUT = .false.! Temp'ly all to stdout
  logical, save, private :: SILENTRUNNING  = .false. ! Suspend all further output
  integer, save, private :: ATCOLUMNNUMBER = 1  ! Where we'll print next
  ! See below for uses of indentBy
  integer, save, private :: INDENTBY       = 0  ! How many spaces to indent
  logical, save, private :: ATLINESTART = .true.! Are we at the line's start?
  integer, save, private :: LINESSINCELASTSTAMP = 0
  logical, private, parameter :: LOGEXTRABLANKS = .false.
  integer, private, parameter :: RECLMAX = 1024  ! This is NAG's limit
  ! logical, save, private :: AlwaysWrap          = .false. ! once past WrapPastColumn
  ! integer, save, private :: WrapPastColumn = 0  ! Don't print beyond (if > 0)

  ! For certain numerical values we will use list directed '*' format
  ! unless optional FORMAT specifier supplied
  double precision, parameter, dimension(3) :: DPREFERDEFAULTFORMAT = &
    & (/ -1.d0, 0.d0, 1.d0 /)  ! For which values to use default format '*'
  real, parameter, dimension(3) :: RPREFERDEFAULTFORMAT = &
    & (/ -1., 0., 1. /)  ! For which values to use default format '*'
  character(len=16), private, save :: NONFINITEFORMAT = '(1pg14.6)' ! 'NaN, Inf'
  ! For printing logical-valued arguments
  integer, private, parameter          :: true  = 1 ! Index in TruthValues
  integer, private, parameter          :: false = 2
  character(len=2), dimension(2), save :: TruthValues = &
    & (/ 'T ', 'F ' /)
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------- Indenting procedures --------
  ! The idea is to have a virtual printed page offset 
  ! from the actual physical page by an amount, indentBy
  ! So if we print something to column k, it shows up
  ! printed in column (k+indentBy)
  !
  ! The indents will affect all calls to:
  !   blanks
  !   output
  ! and subroutines in other modules that USE them.
  ! Indents will not affect logged output or 'print *,'
  ! or calls made directly to PrintItOut
  ! -----------------------------------------------------  addToIndent  -----
  subroutine addToIndent ( n )
  ! add n blanks To Indent
    integer, intent(in) :: n
    indentBy = max( indentBy + n, 0 )
  end subroutine addToIndent

  ! -----------------------------------------------------  resetIndent  -----
  subroutine resetIndent
  ! Reset Indent back to 0
    indentBy = 0
  end subroutine resetIndent

  ! .......................................  Advance_is_yes_or_no  .....
  function Advance_is_yes_or_no ( str ) result ( outstr )
    ! Process the advance argument in
    !   call output( something, [advance='arg1 [arg2] .. [argn]' )
    ! takes '[Yy]...' or '[Nn..] and returns 'yes' or 'no' respectively
    ! also does the same with '[Tt]..' and '[Ff]..'
    ! Returns OutputOptions%advanceDefault if the argument is absent.
    !        A d v a n c e d   u s a g e
    ! All the other usagess require a longer explanation:
    ! (1)
    ! 'beep', 'stderr' or just 'err' is special case--we return 'beep'
    ! trusting that the caller will print to stderr instead of stdout
    ! (2)
    ! if the initial character isn't one of [YyNn], or (1), we return 
    ! OutputOptions%advanceDefault (just as if it were missing)
    ! (3)
    ! We are allowing the argument str to do multiple duties by being
    ! composed of multiple space-separated sub-arguments, e.g. 
    !   'arg1 [arg2] .. [argn]'
    ! the first arg1 is treated as before, basically 'yes' or 'no' or 'beep'.
    ! arg2 and beyond set AdvancedOptions to the
    ! output command:
    !   sub-arg                   meaning
    !   -------                   -------
    !      (These options take effect permanently)
    !     save         save original OutputOptions to be restored later
    !    restore       restore original OutputOptions
    !    unit n        set print unit to n
    !    level k       set severity level to k if calling MLSMessage
    !   newline m      use achar(m) as character for newLine
    !    wrap          Wrap lines longer than WrapPast Column
    !    wrappast col  Set WrapPast Column to col
    !
    !      (These options take effect only temporarily)
    !     pause        pause and wait for user input
    !    stretch       s t r e t c h  characters before printing
    !                  * -------------------------------------- *
    !    banner        *   surround characters with a banner    *
    !                  * -------------------------------------- *
    !    header        *-- surround characters with a header ---*
    !--------Argument--------!
    character (len=*), intent(in), optional :: Str
    character (len=4) :: Outstr ! 'yes', 'no', or 'beep'

    !----------Local vars----------!
    character (len=*), parameter :: yeses = 'YyTt'
    character (len=*), parameter :: nose = 'NnFf'
    integer                      :: kSpace
    ! Executable
    if ( .not. present(str)  ) then
      outstr = outputoptions%advanceDefault ! 'no'
      if ( OutputOptions%AlwaysWrap .and. AtColumnNumber > OutputOptions%WrapPastColumn .and. &
        & OutputOptions%WrapPastColumn > 0 ) &
        & outstr = 'yes'
      return
    end if
    outstr = adjustl(str)
    kSpace = index( outstr, ',' ) ! Did we separate args with a ',' ?
    if ( kspace < 1 ) kSpace = index( outstr, ' ' )
    if ( kSpace > 1 ) outstr = outstr(:kSpace) ! To snip off arg2 ..
    if ( index( yeses, outstr(:1)) > 0  ) then
      outstr = 'yes'
    else if ( index( nose, outstr(:1)) > 0  ) then
      outstr = 'no'
    else if ( index( lowercase(outstr), 'beep' ) > 0  ) then
      outstr = 'beep'
    else if ( index( lowercase(outstr), 'err' ) > 0  ) then
      outstr = 'beep'
    else if ( index( lowercase(outstr), 'stde' ) > 0  ) then
      outstr = 'beep'
    else
      outstr = outputoptions%advanceDefault ! str
    end if
    if ( OutputOptions%AlwaysWrap .and. AtColumnNumber > OutputOptions%WrapPastColumn .and. &
      & OutputOptions%WrapPastColumn > 0 ) &
      & outstr = 'yes'
    if ( len_trim(str) < 4 ) return
    ! write(*,*) 'str       = ', str
    ! Now set more advanced options based on str
    call SetAdvancedOption( str )
    ! Remember to restore previous advanced options
    MustRestoreAdvOpts = .true.
  end function Advance_is_yes_or_no

  ! -----------------------------------------------------  Beep  -----
  subroutine Beep ( chars )
    use, intrinsic :: ISO_Fortran_Env, only: Error_Unit
  ! Print chars to stderr
    character(len=*), intent(in), optional :: chars
    if ( present(chars) ) then
      write( Error_Unit, * ) trim(chars)
    else
      write( Error_Unit, * ) 'Beep'
    endif
  end subroutine Beep

  ! -----------------------------------------------------  Blanks  -----
  subroutine Blanks ( N_Blanks, Fillchar, Advance, Dont_Stamp )
  ! Output N_Blanks blanks to PRUNIT.
  ! or optionally that many copies of fillChar.
  ! If FillChar is one of the SpecialFillChars maintained in outputOptions
  ! we print the corresponding pattern instead of blanks.
    integer, intent(in) :: N_BlankS
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: Fillchar  ! default is ' '
    logical, intent(in), optional          :: Dont_stamp ! Prevent double-stamping
    integer :: I    ! Blanks to write in next Write statement
    logical :: lineup
    integer :: ntimes
    integer :: numSoFar
    character(len=16) :: pattern
    integer :: patternLength
    integer :: patternNum
    integer :: theRest
    ! Executable
    if ( present(fillChar) ) then
      if ( patternOptions%usePatternedBlanks .and. &
        & index(patternOptions%specialFillChars, FILLCHAR) > 0 ) then
        ! We need to try to fit our called-for pattern into n_Blanks
        ! The 1st question is, how many times could it be done?
        patternNum = index( patternOptions%specialFillChars, FillChar )
        pattern = patternOptions%patterns(patternNum)
        ! The pattern length (adjusted for enclosing parentheses)
        patternLength = len_trim(pattern) - 2
        ! Now we assume we'll always want the first and Blanks of n_Blanks to be
        ! purely Blank
        if ( patternLength > n_blanks - 2 ) then
          ! n_blanks too short--just print blanks
          call Pr_Blanks ( n_blanks, advance=advance, dont_stamp=dont_stamp )
          return
        end if
        ntimes = (n_blanks-2)/patternLength
        ! In case we want latterns on consecutive lines to line up; viz
        ! a: . . . . . . Something
        ! xy:. . . . . . Something else
        lineup = ( patternOptions%lineupFillChars(patternNum:patternNum) == 'y' )
        if ( lineup ) then
          numSoFar = 0
          ! Make sure that we always begin on an even-numbered column
          ! (This only works for patterns like '. . . ' or '- - - '
          if ( mod(atColumnNumber, 2) /= 0 ) then
            call Pr_Blanks ( 1, advance='no' )
            numSoFar = 1
          end if
        else
          call Pr_Blanks ( 1, advance='no' )
          numSoFar = 1
        end if
        do i=1, ntimes
          call output_ ( pattern(2:patternLength+1), advance='no' )
          ! if ( xtraBlanks > 0 ) call Pr_Blanks ( xtraBlanks, advance='no' )
          numSoFar = numSoFar + patternLength
        enddo
        theRest = n_blanks - numSoFar
        if ( theRest > 0 ) call Pr_Blanks ( theRest, advance=advance, dont_stamp=dont_stamp )
        return
      end if
    end if
    call Pr_Blanks ( n_blanks, fillChar=fillChar, advance=advance, dont_stamp=dont_stamp )
  end subroutine Blanks

  ! ----------------------------------------------  flushOutputLines  -----
  ! print or log OutputLines
  ! then reset to ''
  subroutine flushOutputLines ( prUnit )
    ! use, intrinsic :: ISO_Fortran_Env, only: Output_Unit
    ! Args
    integer, optional, intent(in) :: prUnit ! How do you want 'em?
    ! Local arguments
    integer :: kNull  ! 1 past where the end of str occurs
    integer :: myPrUnit
    integer :: oldPrUnit
    ! Executable
    myPrUnit = -1 ! By default, just print to stdout
    if ( present(prUnit) ) myPrUnit = prUnit
    oldprUnit = outputOptions%prUnit
    outputOptions%prUnit = myPrUnit
    kNull = index( OutputLines, achar(0), back=.true. )
    if ( kNull >= 2 .and. kNull <= len(OutputLines) ) &
      & call output( OutputLines(1:kNull-1) )
    OutputLines = ' '
    ! flush ( merge ( output_Unit, outputOptions%prUnit, outputOptions%prUnit < 0 ) )
    outputOptions%prUnit = oldPrUnit
  end subroutine flushOutputLines

  ! ----------------------------------------------  FlushStdout  -----
  ! Flush Whichever output unit is being used
  subroutine FlushStdout
    use, intrinsic :: ISO_Fortran_Env, only: Output_Unit
    ! Executable
    flush( Output_Unit )
  end subroutine FlushStdout

  ! ---------------------------------------------- GetOutputStatus
  ! Returns certain normally private data
  ! intended for modules like highOutput and maybe some others
  ! result will be an integer
  ! equal to the value of integer-valued data
  ! or to 1 if the logical-valued data is TRUE, 0 if FALSE
  function GetOutputStatus( name, value ) result( status )
    ! Args
    character(len=*), intent(in)            :: name
    character(len=*), intent(out), optional :: value
    integer                                 :: status
    ! Executable
    status = -999 ! meaning name was not recognized
    if ( present(value) ) value = 'value not relevant for this name'
    if ( index(lowercase(name), 'physicalcolumn' ) > 0 ) then
      status = atColumnNumber            ! This is the "physical" column
    elseif ( index(lowercase(name), 'column' ) > 0 ) then
      status = atColumnNumber - indentBy ! This is the "virtual" column
    elseif( index(lowercase(name), 'indent' ) > 0 ) then
      status = indentBy
    elseif( index(lowercase(name), 'start' ) > 0 ) then
      status = merge(1, 0, atLineStart)
    elseif( index(lowercase(name), 'lines' ) > 0 ) then
      status = linesSincelastStamp
    elseif( index(lowercase(name), 'silent' ) > 0 ) then
      status = merge(1, 0, silentRunning)
    elseif( index(lowercase(name), 'wrappast' ) > 0 ) then
      status = OutputOptions%WrapPastColumn
    elseif( index(lowercase(name), 'wrap' ) > 0 ) then
      status = merge(1, 0, OutputOptions%AlwaysWrap)
    elseif( index(lowercase(name), 'true' ) > 0 .and. present(value)) then
      status = true
      value = TruthValues (true)
    elseif( index(lowercase(name), 'false' ) > 0 .and. present(value)) then
      status = false
      value = TruthValues (false)
    endif
  end function GetOutputStatus

  ! ---------------------------------------------- printOutputStatus
  ! Prints certain normally private data
  ! revealing what settings and options are in force
  subroutine printOutputStatus ( keywords )
    ! Args
    character(len=*), optional, intent(in) :: keywords
    ! Internal variables
    integer                                :: i
    ! Executable
    print *,  'atColumnNumber ', atColumnNumber            ! This is the "physical" column
    print *,  'atColumnNumber - indentBy ', atColumnNumber - indentBy ! This is the "virtual" column
    print *,  'indentBy ', indentBy
    print *,  'atLineStart? ', merge(1, 0, atLineStart)
    print *,  'linesSincelastStamp ', linesSincelastStamp
    print *,  'silentRunning ', merge(1, 0, silentRunning)
    if ( .not. present(keywords) ) return
    ! keyword tells which Output Options to print
    if ( ok ( 'unit'   ) )  print *, 'prUnit          ', outputOptions%prUnit
    if ( ok ( 'name'   ) )  print *, 'name            ', trim(outputOptions%name)
    if ( ok ( 'newlin' ) )  print *, 'newLineVal      ', outputOptions%newLineVal
    if ( ok ( 'sddef ' ) )  print *, 'sdFormatDefault ', outputOptions%sdFormatDefault
    if ( ok ( 'parent' ) )  print *, 'parentName      ', outputOptions%parentName
    if ( ok ( 'pattern' ))  then
      call output( 'Special patterns used in blanks ', advance='yes' )
      do i=1, numPatterns
        call output( 'pattern: ', advance='no' )
        call output( patternOptions%specialFillChars( i:i ), advance='no' )
        call blanks(1)
        call blanks( 64, FillChar=patternOptions%specialFillChars( i:i ) )
        call NewLine
      enddo
    endif
  contains
    logical function ok ( arg )
      ! Tell us whether keyword means it's ok to print arg
      ! '*' is the wildcard
      character(len=*), intent(in)   :: arg
      ok = ( keywords == '*' .or. index( lowercase(keywords), arg ) > 0 )
    end function ok
  end subroutine printOutputStatus

  ! ----------------------------------------------  IsOutputSuspended  -----
  logical function IsOutputSuspended ()
  ! Have we suspended outputting to PRUNIT?
    IsOutputSuspended = silentRunning
  end function IsOutputSuspended

  ! ----------------------------------------------------  NewLine  -----
  subroutine NewLine ( dont_make_blank_line )
    ! Go on to the next printed line
    ! unless already at column 1 and dont_make_blank_line is true
    ! Args
    ! If the following is TRUE, avoid
    ! adding a blank line; that means don't add
    ! a new line if at column 1
    logical, optional, intent(in) :: dont_make_blank_line
    ! Executable
    if ( present( dont_make_blank_line ) ) then
      if ( dont_make_blank_line .and. ATCOLUMNNUMBER == 1 ) return
    endif
    call output_ ( '', advance='yes' )
  end subroutine NewLine

  ! ------------------------------------------------  Output_Char  -----
  ! Output CHARS to PRUNIT.
  subroutine Output_Char ( Chars, &
    & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp, &
    & NewLineVal, dont_asciify, format )
    ! We will 1st check to see whether any internal characters are
    ! codes for newlines
    ! If any are, we will call newLine in place of printing
    ! them
    ! (This is a new default behavior; you can restore
    ! the old by passing an impossible value for NewLineVal, e.g. -999)
    !
    ! If NewLineVal is present, we'll use it as the code for new lines
    ! Otherwise, we default to outputOptions%newlineVal
    character(len=*), intent(in)           :: Chars
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: From_Where
    logical, intent(in), optional          :: Dont_Log ! prevent double-logging
    character(len=*), intent(in), optional :: Log_Chars
    character(len=*), intent(in), optional :: InsteadOfBlank ! What to output
    logical, intent(in), optional          :: Dont_Stamp ! prevent double-stamping
    integer, intent(in), optional          :: Newlineval ! What char val to treat as <cr>
    logical, intent(in), optional          :: Dont_asciify ! output binary
    character(len=*), intent(in), optional :: Format ! consistent with generic
    ! Internal variables
    integer :: I ! loop inductor
    integer :: BannerLen
    integer :: indent ! Should have chosen different name
    integer :: LineLen ! How many chars to print
    character(len=4) :: MY_ADV ! 'yes' if advance
    integer :: myNewLineVal
    logical :: myAsciify ! Convert non-printing chars to ascii?
    character(len=max(90, 2*(len(chars)+15))) :: newChars ! What to print
    ! Executable
    my_adv = Advance_is_yes_or_no(advance)
    ! call Dump_AdvancedOptions
    myNewLineVal = outputOptions%newlineVal
    if ( present(newLineVal) ) myNewLineVal = newLineVal
    myAsciify = .true.
    if ( present(dont_asciify) ) myAsciify = .not. dont_asciify
    ! If bannered, how large
    BannerLen = max( 80, 4 + len_trim(chars) )
    LineLen = len(chars)
    indent = ( BannerLen - LineLen ) / 2 ! spaces beetween '*' and start of chars
    newChars = chars
    if ( advancedOptions%bannered ) then
      call Output_Char_NoCR ( '*', advance='no' )
      call Output_Char_NoCR ( Repeat('-', BannerLen-2 ), advance='no' )
      call Output_Char_NoCR ( '*', advance='yes' )
    endif
    if ( advancedOptions%headered ) then
      newChars = '* ---- ' // newChars(:LineLen) // ' ---- *'
      LineLen = LineLen + 2*7
    endif
    if ( advancedOptions%stretch ) then
      newChars = stretch( newChars, options='-a' )
      LineLen = 2*LineLen - 1
      BannerLen = max( 80, 4 + LineLen )
      indent = ( BannerLen - LineLen ) / 2
    endif
    ! write (*,*) 'bannered            =  ', advancedOptions%bannered 
    if ( advancedOptions%bannered ) then
      newChars = '*' // Repeat( ' ', indent-1 ) // newChars
      LineLen = Linelen + indent
      newChars(BannerLen:BannerLen) = '*'
      LineLen = Bannerlen
      ! print *, 'newChars ', newChars
    endif
    ! Are any internal chars new line vals?
    i = index( chars, achar(myNewLineVal) )
    if ( i < 1 ) then
      ! No internal new lines
      if ( myAsciify ) newChars = ReplaceNonAscii( newChars, '@', exceptions=achar(9))
      if ( myAsciify ) then
        call Output_Char_NoCR ( newChars(:LineLen), &
          & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
      else
        call Output_Char_NoCR ( newChars(:LineLen), &
          & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
      endif
    else
      ! Print every character one-by-one except internal new lines, at
      ! which we'll call NewLine instead
      do i=1, len(chars)
        if ( chars(i:i) /= achar(myNewLineVal) ) then
          if ( myAsciify ) then
            call Output_Char_NoCR ( ReplaceNonAscii(CHARS(i:i), '@', exceptions=achar(9)), &
              & advance='no', from_where=from_where, dont_log=dont_log, &
              & log_chars=log_chars, insteadofblank=insteadofblank, &
              & dont_stamp=dont_stamp )
          else
            call Output_Char_NoCR ( CHARS(i:i), &
              & advance='no', from_where=from_where, dont_log=dont_log, &
              & log_chars=log_chars, insteadofblank=insteadofblank, &
              & dont_stamp=dont_stamp )
          endif
        else
          call newLine
        endif
      enddo
      if ( my_adv == 'yes' ) call newLine
    endif
    ! write (*,*) 'bannered            =  ', advancedOptions%bannered 
    if ( advancedOptions%bannered ) then
      call Output_Char_NoCR ( '*', advance='no' )
      call Output_Char_NoCR ( Repeat('-', BannerLen-2 ), advance='no' )
      call Output_Char_NoCR ( '*', advance='yes' )
    endif
    if ( MustRestoreAdvOpts ) then
      AdvancedOptions = DefaultAdvancedOptions
      MustRestoreAdvOpts = .false.
    endif
  end subroutine Output_Char

  subroutine Output_Char_NoCR ( Chars, &
    & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
    character(len=*), intent(in)           :: Chars
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: From_Where
    logical, intent(in), optional          :: Dont_Log ! prevent double-logging
    character(len=*), intent(in), optional :: Log_Chars
    character(len=*), intent(in), optional :: InsteadOfBlank ! What to output
    logical, intent(in), optional          :: Dont_Stamp ! prevent double-stamping
    !
    character(len=4) :: MY_ADV
    ! Executable
    my_adv = Advance_is_yes_or_no(advance)
    atLineStart = (my_adv == 'yes')
    if ( indentBy > 0 .and. atColumnNumber == 1 ) then
      call Output_Char_NoCR_indented ( repeat( ' ', indentby ) // chars, &
        & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
    else
      call Output_Char_NoCR_indented ( chars, &
        & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
    endif
    ! Were we asked to Pause?
    ! If so read 1 char from stdin
    ! Note: this may hang a job running under sips control so strip
    ! out any debugging use of this before delivery
    ! print *, 'pause: ', advancedOptions%pause
    if ( advancedOptions%pause ) then
      ! print *, '(P a u s e d .. e n t e r   o k   t o   r e s u m e, p   t o   s t e p)'
      ! read (*,'(a)') my_adv(1:1)
      ! Reset--unless told to ste'p', must explicitly request next pause
      call Pause ( my_adv, Prompts=(/ &
        & '(P a u s e d .. e n t e r   o k   t o   r e s u m e, p   t o   s t e p)' &
        & /) )
      if ( index('pP', my_adv(1:1)) < 1 ) advancedOptions%pause = .false. 
    endif
    if ( MustRestoreAdvOpts ) then
      AdvancedOptions = DefaultAdvancedOptions
      MustRestoreAdvOpts = .false.
    endif
  end subroutine Output_Char_NoCR

  subroutine Output_Char_NoCR_Indented ( chars, &
    & advance, from_where, dont_log, log_chars, insteadofblank, dont_stamp )
    ! -------------------------------------------------------------------
    ! We have arrived at
    !   T h e   w o r k h o r s e
    ! We have taken care of linefeeds, indents, and numeric conversions
    ! so now we do one of
    ! (1) print to stderr if advance = 'stderr' or 'beep'
    ! (2) go silent if we are skipping all output
    ! (3) write to a file unit if PrUnit is to be taken literally
    ! (4) append to our accumulated OutputLines if we are deferring output
    ! (5) print immediately if we have been switched to stdout
    ! (6) wade through a thicket of sub-choices involving possibly
    !     (a) time stamps, either individually or grouped
    !     (b) logging, either instead of printing or in addition
    ! -------------------------------------------------------------------
    ! Should we rename this subroutine something less opaque, less obscure?
    ! Not precisely the lowest level, because we may call 
    !   Beep
    !   MyMesssage
    !   PrintItOut
    use, intrinsic :: ISO_Fortran_Env, only: Output_Unit
    use SDPToolkit, only: Stamp
    character(len=*), intent(in)           :: Chars
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: From_Where
    logical, intent(in), optional          :: Dont_Log ! prevent double-logging
    character(len=*), intent(in), optional :: Log_Chars
    character(len=*), intent(in), optional :: InsteadOfBlank ! What to output
    logical, intent(in), optional          :: Dont_Stamp ! prevent double-stamping
    !
    ! logical :: alreadyLogged
    logical :: DoIt    ! TheUnit /= 0 .or. outputOptions%prUnitLiteral
    integer :: i1, i2
    integer :: IOBloc
    integer :: nIOBlocs
    character(len=max(16,len(chars)+1)) :: my_chars
    character(len=len(chars)+64) :: stamped_chars ! What to print to stdout
    character(len=max(16,len(chars)+1)) :: the_chars
    logical :: my_dont_log
    logical :: my_dont_stamp
    character(len=4) :: MY_ADV
    integer :: n_chars
    integer :: n_stamp ! How much of stamped_chars to print
    logical :: stamp_header
    logical :: stamped
    integer :: status
    integer :: TheUnit ! zero for PrUnit == MSGLOGPRUNIT, else unit number
    ! Executable
    n_chars = len(chars)
    alreadyLogged = .false.
    my_adv = Advance_is_yes_or_no( advance )
    ! Print to stderr instead?
    if ( my_adv == 'beep' ) then
      call Beep( chars )
      if ( MustRestoreAdvOpts ) then
        AdvancedOptions = DefaultAdvancedOptions
        MustRestoreAdvOpts = .false.
      endif
      return
    endif
    if ( SILENTRUNNING ) go to 9 ! When we skip all output
    ! print *, 'outputOptions%prUnit: ', outputOptions%prUnit
    ! print *, 'outputOptions%prUnitLiteral: ', outputOptions%prUnitLiteral
    ! Do any special orders apply to this output?
    if ( outputOptions%prunit ==  OutputLinesPrUnit ) then
      ! Append to OutputLines; maybe print later on
      call append_chars( OutputLines, chars )
      if ( my_adv == 'yes' ) &
        &  call append_chars( OutputLines, achar(outputOptions%NewLineVal) )
      ! print *, 'Appending to OutputLines; now ', trim(OutputLines)
      go to 9
    elseif ( outputOptions%prUnitLiteral ) then
      write( outputOptions%prUnit, '(a)', advance=my_adv ) chars
      ! print *, chars
      go to 9
    elseif ( SWITCHTOSTDOUT ) then
      write( *, '(a)', advance=my_adv ) chars
      ! print *, chars
      go to 9
    end if
    ! If we're not advancing and we've been a string of length 0 to print ..
    if ( my_adv == 'no' .and. len(chars) < 1 ) return
    my_dont_stamp = stampOptions%neverStamp ! .false.
    if ( present(dont_stamp) ) my_dont_stamp = dont_stamp
    my_dont_stamp = ( my_dont_stamp .or. &
      & linesSinceLastStamp < (stampOptions%interval - 1) )
    stamped = .false.
    stamp_header = .false.
    stamped_chars = chars
    theUnit = 0
    if ( outputOptions%prUnit > 0 .or. outputOptions%prUnitLiteral ) then
      theUnit = outputOptions%prUnit
    else if ( useStdout(outputOptions%prUnit) ) then
      theUnit = output_unit
    end if
    doIt = TheUnit /= 0 .or. outputOptions%prUnitLiteral
    ! Are we trying to avoid buffered output?
    ! The good part is that if we crash hard we won't lose that
    ! last line.
    if ( (.not. outputOptions%buffered) .and. doIt .and. &
      & ( atcolumnnumber == 1 ) ) then
      flush ( theUnit, iostat=status )
      if ( status /= 0 ) call PrintItOut ( &
        & trim('Unable to flush prUnit ' // outputOptions%name), &
        & 1, exitStatus=1 )
    endif
    ! Do we need to stamp this line? If so, at beginning or at end?
    if ( my_dont_stamp ) then
    elseif ( stampoptions%interval > 9 ) then
      stamp_header = (my_adv == 'yes')
    elseif ( stampoptions%post ) then
      if ( my_adv == 'yes' ) then
        stamped_chars = stamp ( chars, stampOptions%dateFormat, &
          & stampOptions%timeFormat, stampOptions%textCode, &
          & stampOptions%post, stampOptions%showTime )
        stamped = .true.
      end if
    elseif( ATLINESTART ) then
      stamped_chars = stamp ( chars, stampOptions%dateFormat, &
          & stampOptions%timeFormat, stampOptions%textCode, &
          & stampOptions%post, stampOptions%showTime )
      stamped = .true.
    end if

    if (LOGEXTRABLANKS) n_chars = max(len(chars), 1)
    my_dont_log = outputOptions%skipmlsmsglogging ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    n_stamp = len_trim(stamped_chars)
    if ( my_adv == 'no' ) n_stamp = n_stamp + len(chars) - len_trim(chars)
    ! Special case: if chars is blank (chars are blank?)
    ! we'll want to print anyway
    if ( len_trim(chars) < 1 ) n_stamp = max(n_stamp, 1)
    ! print *, 'n_stamp: ', n_stamp
    ! print *, 'len(chars): ', len(chars)
    alreadyLogged = .true.
    if ( doIt .and. n_stamp > RECLMAX ) then
      nIOBlocs = 1 + (n_stamp-1)/RECLMAX
      i2 = 0
      do IOBloc=1, nIOBlocs
        i1 = i2 + 1
        i2 = min(i2+RECLMAX, n_stamp)
        write ( theUnit, '(a)', advance='no' ) stamped_chars(i1:i2)
      enddo
      if ( my_adv == 'yes' ) write ( theUnit, '(a)' )
    elseif ( doIt .and. len(chars) < 1 .and. my_adv == 'yes' ) then
      write ( theUnit, '(a)', advance=my_adv )
      ! print *, 'wrote: ', stamped_chars(1:n_stamp), ' to ', theUnit
    elseif ( doIt .and. n_stamp > 0 ) then
      write ( theUnit, '(a)', advance=my_adv ) stamped_chars(1:n_stamp)
      ! print *, 'wrote: ', stamped_chars(1:n_stamp), ' to ', theUnit
    else
      alreadyLogged = .false.
    endif
    if ( any(outputOptions%prunit == (/MSGLOGPRUNIT, BOTHPRUNIT/)) .and. &
      & .not. my_dont_log .and. .not. outputOptions%prUnitLiteral ) then
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
        call myMessage ( outputOptions%MLSMSG_Level, from_where, &
          & my_chars(1:n_chars), &
          & advance=my_adv )
      elseif ( outputOptions%logParent ) then
        call myMessage ( outputOptions%MLSMSG_Level, &
          & outputOptions%parentName, &
          & my_chars(1:n_chars), &
          & advance=my_adv )
      else
        call myMessage ( outputOptions%MLSMSG_Level, &
          & ModuleName, &
          & my_chars(1:n_chars), &
          & advance=my_adv )
      end if
    end if

    ! print *, 'alreadyLogged: ', alreadyLogged
    ! print *, 'theUnit: ', theUnit
    if ( (outputOptions%prunit <= 0 .and. .not. outputOptions%prunitLiteral) &
      & .or. alreadyLogged ) then
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
    if ( alreadyLogged .and. outputOptions%prunitLiteral ) &
      & flush ( theUnit, iostat=status )
    atLineStart = (my_adv == 'yes')
    if ( atLineStart ) then
      ! Are we trying to avoid buffered output?
      if ( (.not. outputOptions%buffered) .and. doIt ) then
        flush ( theUnit, iostat=status )
        if ( status /= 0 ) call myMessage ( MLSMSG_Error, ModuleName, &
          & trim('Unable to flush prUnit ' // outputOptions%name) )
      endif
      if ( stamp_header ) then
        ! Time to add stamp as a page header on its own line
        stamped_chars = stamp ( ' ', stampOptions%dateFormat, &
          & stampOptions%timeFormat, stampOptions%textCode, &
          & stampOptions%post, stampOptions%showTime )
        stamped = .true.
        if ( doIt .and. len_trim(stamped_chars) > 0 ) then
          write ( theUnit, '(a)', advance='yes' ) trim(stamped_chars)
        end if
      end if
      if ( stamped ) then
        linesSinceLastStamp = 0
      else
        linesSinceLastStamp = linesSinceLastStamp + 1
      end if
    end if
9   continue
    atColumnNumber = atColumnNumber + n_chars
    if ( atLineStart ) atColumnNumber = 1
    if ( MustRestoreAdvOpts ) then
      AdvancedOptions = DefaultAdvancedOptions
      MustRestoreAdvOpts = .false.
    endif
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
  end subroutine Output_Char_NoCR_INDENTED

  ! ------------------------------------------  Output_Char_ARRAY  -----
  subroutine Output_Char_ARRAY ( CHARS, ADVANCE_AFTER_EACH, ADVANCE, &
    & INSTEADOFBLANK, NEWLINEVAL, format )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS(:)
    character(len=*), intent(in), optional :: Advance_after_each
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: Insteadofblank ! what to output
    integer, intent(in), optional          :: Newlineval ! what char val to treat as <cr>
    character(len=*), intent(in), optional :: Format ! consistent with generic
    ! Internal variables
    integer :: I ! loop inductor
    ! Executable
    do i = 1, size(chars)
      call output ( chars(i), &
        & insteadofblank=insteadofblank, &
        & newLineVal=newLineVal, advance=ADVANCE_AFTER_EACH )
      if ( len(chars(1)) > 1 .and. .not. present(ADVANCE_AFTER_EACH) ) &
        & call SeparateElements( i, size(chars) )
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance )
    end if
  end subroutine Output_Char_ARRAY

  ! ---------------------------------------------  OUTPUT_COMPLEX  -----
  subroutine OUTPUT_COMPLEX ( VALUE, Format, Advance, Before, After, dont_stamp )
    complex, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: Dont_Stamp
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

  ! --------------------------------------------  OUTPUT_DCOMPLEX  -----
  subroutine OUTPUT_DCOMPLEX ( VALUE, Format, Advance, Before, After )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: Advance
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
  subroutine OUTPUT_DOUBLE ( VALUE, Format, LogFormat, Advance, &
    & Before, After, dont_stamp, Trim )
  ! Output "double" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: LogFormat ! How to post to Log
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: Dont_Stamp
    logical, intent(in), optional :: Trim ! Trim blanks even if Format present
    logical :: Do_Trim
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
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    integer :: I ! loop inductor
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      call SeparateElements( i, size(values) )
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
    character(len=*), intent(in), optional :: Advance
    logical, intent(in), optional :: FILL
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, intent(in), optional :: Dont_Stamp
    !
    logical :: My_Fill
    integer :: I, J
    character(len=12) :: LINE
    character(len=12) :: myFormat
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
      myFormat = '(i7)'
      if ( len_trim(format) > 1 ) myFormat = format
      write ( line, myFormat ) int
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
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: FORMAT
    logical, optional, intent(in) :: DONT_STAMP
    integer :: I ! loop inductor
    do i = 1, size(integers)
      call output ( integers(i), advance='no', format=format )
      call SeparateElements( i, size(integers) )
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_INTEGER_ARRAY

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL ( LOG, Advance, Before, DONT_STAMP, format )
  ! Output LOG to PRUNIT using at most PLACES (default zero) places
    logical, intent(in) :: LOG
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: BEFORE
    logical, optional, intent(in) :: DONT_STAMP
    character(len=*), intent(in), optional :: format ! consistent with generic
    character(len=2) :: LINE
    if ( log ) then
      line = TruthValues( true ) ! ' T'
    else
      line = TruthValues( false ) ! ' F'
    end if
    if ( present(before) ) call output_ ( before, DONT_STAMP=DONT_STAMP )
    call output_ ( line, advance=advance, DONT_STAMP=DONT_STAMP )
  end subroutine OUTPUT_LOGICAL

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL_ARRAY ( logs, &
    & Advance, Before, DONT_STAMP, ONLYIF, format )
    ! Output LOG to PRUNIT using at most PLACES (default zero) places
    ! Optionally, print non-blank only if T (or F)
    logical, dimension(:), intent(in) :: logs
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: BEFORE
    logical, optional, intent(in) :: DONT_STAMP
    logical, optional, intent(in) :: ONLYIF ! Print only if true (false)
    character(len=*), intent(in), optional :: format ! consistent with generic
    ! Internal variables
    character(len=1), dimension(size(logs)) :: clogs
    character(len=size(logs)) :: logChars
    integer :: I ! loop inductor
    character(len=1) :: ifonlyWhat
    ! Executable
    if ( present(before) ) call output_ ( before, DONT_STAMP=DONT_STAMP )
    if ( present(onlyif) ) then
      if ( onlyif ) then
        ifonlyWhat = TruthValues( true ) ! ' T'
      else
        ifonlyWhat = TruthValues( false ) ! ' F'
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
        call SeparateElements( i, size(logs) )
      end do
    endif
    if ( present(advance) ) &
      & call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
  end subroutine OUTPUT_LOGICAL_ARRAY

  ! ----------------------------------------------  OUTPUT_SINGLE  -----
  subroutine OUTPUT_SINGLE ( VALUE, FORMAT, LogFormat, ADVANCE, &
    & Before, After, DONT_STAMP, Trim )
  ! Output "SINGLE" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format  ! How to print
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: Before, After ! text to print
    logical, optional, intent(in) :: DONT_STAMP
    logical, intent(in), optional :: Trim ! Trim blanks even if Format present
    logical :: Do_Trim
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
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    integer :: I ! loop inductor
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      call SeparateElements( i, size(values) )
    end do
    if ( present(advance)  ) then
      call output_ ( '', advance=advance, DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_SINGLE_ARRAY

  ! ----------------------------------------------  OUTPUT_STRING  -----
  subroutine OUTPUT_STRING ( STRING, LENSTRING, ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, format )
  ! Output STRING to PRUNIT.
    character(len=*), intent(in) :: STRING
    integer, intent(in) :: LENSTRING
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: format ! consistent with generic
    integer :: n_chars
    !
    n_chars = min(len(string), lenstring)
    if ( len(string) < 1  ) then
      call myMessage ( MLSMSG_Error, ModuleName, &
        & 'Bad string arg in OUTPUT_STRING' )
    else if ( len_trim(string) < 1 .or. LENSTRING < 1  ) then
      call output_ ( '', advance )
    else
      call output_ ( string(:n_chars), advance, from_where, dont_log, log_chars )
    end if
  end subroutine OUTPUT_STRING

  ! ----------------------------------------------  restoreSettings  -----
  subroutine restoreSettings ( settings )
  ! resume outputting to PRUNIT.
  ! optionally set to use or ignore Toolkit
    use PrintIt_m, only: DefaultLogUnit, Set_Config, StdoutLogUnit
    character(len=*), optional, intent(in) :: settings
    ! Local variables
    character(len=*), parameter            :: allSettings = &
      & 'toolkit,stamp,output,time'
    character(len=len(allSettings))        :: mySettings 
    logical                                :: useToolkit
    ! Executable
    mySettings = ' '
    if ( present(settings) ) mySettings = settings
    if ( index(mySettings, '*') > 0 ) mySettings = allSettings
    mySettings = lowercase(mySettings)
    if ( index(mySettings, 'output') > 0 ) outputOptions = DefaultOutputOptions

    if ( index(mySettings, 'stamp') > 0 ) stampOptions = DefaultStampOptions

    if ( index(mySettings, 'time') > 0 ) timeStampOptions = DefaultTimeStampOptions
    useToolkit = ( index(mySettings, 'toolkit') > 0 )
    if ( .not. useToolkit ) return
    call set_config ( useToolkit = useToolkit, &
      & logFileUnit=merge(defaultLogUnit, stdoutLogUnit, useToolkit) )
  end subroutine restoreSettings

  ! ----------------------------------------------  resumeOutput  -----
  subroutine resumeOutput 
  ! resume outputting to PRUNIT.
  ! Reverses effect of suspendOutput
    silentRunning = .false.
  end subroutine resumeOutput

  ! ----------------------------------------------  revertOutput  -----
  subroutine revertOutput
  ! revert to outputting to OLDUNIT. Close current PRUNIT (if > 0 and open)
  ! Unless switchOutput was called with filename == 'stdout'
    ! Local variables
    logical :: itsOpen
    ! Executable
    if ( SWITCHTOSTDOUT ) then
      SWITCHTOSTDOUT = .false.
      return
    endif
    call resumeOutput
    if ( .not. OLDUNITSTILLOPEN ) then
      call output_( 'Unable to Revert output--old unit not open', advance='yes' )
      return
    end if
    if ( outputOptions%prunit > 0 ) then
      inquire( unit=outputOptions%prunit, opened=itsOpen )
      if ( itsOpen ) then
        close(outputOptions%prunit)
      end if
    end if
    outputOptions%prunit = OLDUNIT    
    call output_( 'Reverting output to unit: ', advance='no' )
    call output( OLDUNIT, advance='yes' )

  end subroutine revertOutput

  ! ----------------------------------------------  setFillPattern  -----
  subroutine SetFillPattern ( pattern, fillChar )
  ! Override and set a special Fill pattern that can be used in a call to blanks
  ! e.g., by call blanks, FillChar='0' )
    character(len=1), optional, intent(in) :: fillChar
    character(len=*), optional, intent(in) :: pattern
    ! Local variables
    integer :: patternNum
    ! Executable
    patternNum = 1
    if ( present(fillChar) ) &
      & patternNum = index( patternOptions%specialFillChars, fillChar )
    patternOptions%patterns(patternNum) = pattern
    if ( index( pattern, '(' ) < 1 ) &
      & patternOptions%patterns(patternNum) = '(' // pattern // ')'
  end subroutine SetFillPattern

  ! -----------------------------------------------------  SetAdvancedOption  -----
  subroutine SetAdvancedOption ( str )
  ! Set advanced or miscellaneous options
  
  ! Note: these may be more conveniently done via setOutputStatus
  ! The way we do it is to painfully try picking str apart
  ! See how the options explained under Advance_is_yes_or_no
    character (len=*), intent(in), optional :: Str
    ! Local variables
    ! character (len=3)                       :: adv
    character (len=32)                      :: val
    logical                                 :: logval
    ! Executable
    ! adv =  Advance_is_yes_or_no ( str )
    ! AdvancedOptions%originalOptions = outputOptions
    ! Now check for changing the advanced options
    call getOption ( lowercase(str), 'save', logval, initialize=.true. )
    if ( logval ) then
      advancedOptions%originalOptions = outputOptions
    endif
    call getOption ( lowercase(str), 'restore', logval, initialize=.true. )
    if ( logval ) then
      outputOptions = advancedOptions%originalOptions
    endif
    call getOption ( lowercase(str), 'unit', val, initialize=.true. )
    if ( len_trim(val) > 0 ) read ( val, * ) outputOptions%prUnit
    call getOption ( lowercase(str), 'level', val, initialize=.true. )
    if ( len_trim(val) > 0 ) read ( val, * ) outputOptions%MLSMSG_Level
    call getOption ( lowercase(str), 'newline', val, initialize=.true. )
    if ( len_trim(val) > 0 ) read ( val, * ) outputOptions%newLineVal
    call getOption ( lowercase(str), 'wrappast', val, initialize=.true. )
    if ( len_trim(val) > 0 ) read ( val, * ) OutputOptions%WrapPastColumn
    call getOption ( lowercase(str), 'stretch', advancedOptions%stretch )
    call getOption ( lowercase(str), 'banner', advancedOptions%bannered )
    call getOption ( lowercase(str), 'header', advancedOptions%headered )
    call getOption ( lowercase(str), 'wrap', OutputOptions%AlwaysWrap )
    if ( .not. advancedOptions%pause ) &
      & call getOption ( str, 'pause', advancedOptions%pause )
    ! print *, 'pause: ', advancedOptions%pause
    if ( DeeBug ) print *, 'WrapPastColumn: ', OutputOptions%WrapPastColumn
    if ( DeeBug ) print *, 'AlwaysWrap: ', OutputOptions%AlwaysWrap
  end subroutine SetAdvancedOption

  ! ---------------------------------------------- setOutputStatus
  ! Sets certain normally private data
  ! Sets for modules like highOutput and maybe some others
  ! Effect will be an integer
  ! equal to value if integer-valued data
  ! or to TRUE if value is 1
  subroutine SetOutputStatus( name, value )
    ! Args
    character(len=*), intent(in) :: name
    integer, intent(in)          :: value
    ! Executable
    if ( index(lowercase(name), 'physicalcolumn' ) > 0 ) then
      atColumnNumber = value           ! This is the "physical" column
    elseif( index(lowercase(name), 'indent' ) > 0 ) then
      indentBy = value
    elseif( index(lowercase(name), 'start' ) > 0 ) then
      atLineStart = ( value == 1 )
    elseif( index(lowercase(name), 'lines' ) > 0 ) then
      linesSincelastStamp = value
    elseif( index(lowercase(name), 'silent' ) > 0 ) then
      silentRunning = ( value == 1 )
    elseif( index(lowercase(name), 'wrappast' ) > 0 ) then
      OutputOptions%WrapPastColumn = value
    elseif( index(lowercase(name), 'wrap' ) > 0 ) then
      OutputOptions%AlwaysWrap = ( value == 1 )
    endif
  end subroutine SetOutputStatus

  ! ----------------------------------------------  setTruthPattern  -----
  subroutine SetTruthPattern ( TrueFalse )
  ! Override and set a special Truth pattern that can be used in calls to
  ! output logical-valued scalars and arrays e.g., 
  !  call output( logs, ..)
    character(len=2), dimension(2), intent(in) :: TrueFalse ! (/ 'T ', 'F ' /)
    TruthValues = TrueFalse
  end subroutine SetTruthPattern

  ! ----------------------------------------------  suspendOutput  -----
  subroutine SuspendOutput 
  ! suspend outputting to PRUNIT. Run silent.
  ! Reversible by calling resumeOutput
    silentRunning = .true.
  end subroutine SuspendOutput

  ! ----------------------------------------------  switchOutput  -----
  subroutine SwitchOutput ( filename, unit, keepOldUnitOpen )
  ! stop outputting to PRUNIT. Switch to filename [using unit if supplied]
  ! Special use: if filename == 'stdout', just temporarily print to stdout
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
    ! Check for filename == 'stdout'
    if ( lowercase(filename) == 'stdout' ) then
      SWITCHTOSTDOUT = .true.
      return
    endif
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
  end subroutine SwitchOutput

  ! ------------------ Private procedures -------------------------
  ! ----------------- Dump_AdvancedOptions
  subroutine Dump_AdvancedOptions
    write (*,*) 'stretch             =  ', advancedOptions%stretch  
    write (*,*) 'bannered            =  ', advancedOptions%bannered 
    write (*,*) 'headered            =  ', advancedOptions%headered 
    ! write (*,*) '                       ', advancedOptions
  end subroutine Dump_AdvancedOptions
  
  ! .............................................  getOption  .....
  ! This family of subroutines parses a multipart advance arg into
  ! its components, returning an appropriate value
  ! Example, say the component is marked by the '-S' flag
  ! value type     component   returned value
  !  logical          -S         true
  ! character       -Sxyz        xyz
  !
  ! You may optionally insert a space between S and its value, e.g.
  !                 "-S xyz" also returns "xyz"
  subroutine getOption_char ( arg, flag, val, initialize )
    ! Args
    character(len=*), intent(in)      :: arg
    character(len=*), intent(in)      :: flag
    character(len=*), intent(inout)   :: val
    logical, optional, intent(in)     :: initialize
    ! Local variables
    integer :: kFlag, kNext, flagLen
    ! Executable
    if ( present(initialize) ) val = ' '
    kFlag = index( arg, trim(flag) )
    ! DeeBug = ( index(lowercase(arg), 'wrappast' ) > 0 )
    if ( DeeBug ) print *, 'arg: ', trim(arg)
    if ( DeeBug ) print *, 'flag: ', trim(flag)
    if ( DeeBug ) print *, 'kFlag: ', kFlag
    if ( kFlag < 1 ) return
    val = ' '
    flagLen = len(flag)
    if ( DeeBug ) print *, 'flagLen: ', flagLen
    ! Find start of next component
    kNext = index( arg(kFlag+flagLen+1:), ' ' )
    if ( kNext < 1 ) then
      val = arg(kFlag+flagLen+1:)
    else
      val = arg(kFlag+flagLen+1:kFlag+flagLen+kNext)
    endif
    if ( index(arg, 'wrappast' ) < 1 ) return
    if ( DeeBug ) print *, trim(arg), kflag, flagLen, kNext
    if ( DeeBug ) print *, trim(flag)
    if ( DeeBug ) print *, trim(val)
    if ( DeeBug ) print *, 'Returning'
  end subroutine getOption_char

  subroutine getOption_log ( arg, flag, val, initialize )
    ! Args
    character(len=*), intent(in)      :: arg  
    character(len=*), intent(in)      :: flag 
    logical, intent(inout)            :: val  
    logical, optional, intent(in)     :: initialize
    ! Local variables
    integer :: kFlag
    ! Executable
    if ( present(initialize) ) val = .false.
    kFlag = index( arg, trim(flag) )
    if ( kFlag > 0 ) val = .true.
  end subroutine getOption_log

  ! ------------------------------------  SeparateElements  -----
  ! insert blanks or separator between consecutive elements while outputting
  subroutine SeparateElements (i, n )
    ! Args
    integer, intent(in) :: i ! Element number
    integer, intent(in) :: n ! Number of elements
    ! Executable
    if ( i >= n ) return
    if ( OutputOptions%WrapPastColumn > 0 .and. AtColumnNumber >= OutputOptions%WrapPastColumn ) then
      call newLine
      return
    endif
    if ( OutputOptions%WrapPastColumn == 0 .and. &
      & mod(i, outputOptions%nArrayElmntsPerLine) == 0 ) then
      call output_ ( '', advance='yes', DONT_STAMP=.true. )
      return
    endif
    if ( len_trim(outputOptions%arrayElmntSeparator) > 0 ) then
      call output_( outputOptions%arrayElmntSeparator, advance='no' )
    else
      call blanks ( outputOptions%nBlanksBtwnElmnts, advance='no' )
    endif
  end subroutine SeparateElements

  ! ------------------------------------  myMessage  -----
  subroutine myMessage ( Severity, ModuleNameIn, Message, &
    & Advance, newLineCode )

    ! Print a message (unless printing is suppressed).  If it has %[Nn]
    ! in it, replace that with newline.

    ! We use defaultNewLineCode instead of '%n' because there were some
    ! circumstances where we actally printed '%n'; e.g., 
    ! 'hGrid%noProfsUpperOverlap'

    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in), optional :: Advance ! Do not advance
    character (len=*), intent(in), optional :: newLineCode ! Instead of %n
    !                                 if present and the first character is 'N'
    !                                 or 'n'

    ! Local variables
    logical :: AllOfIt                  ! Print all of it (no %n or %N remains)
    integer :: L1, L2                   ! How far in the line have we printed?
    character (len=512), save :: Line   ! Line to output, should be long enough
    integer, save :: Line_len=0         ! Number of saved characters in line.
    integer :: LogFileUnit
    logical :: My_adv
    character(len=2) :: myNewLineCode

    ! Executable code
    my_adv = .true.
    if ( present(advance) ) &
      & my_adv = advance(1:1) /= 'n' .and. advance(1:1) /= 'N'
    myNewLineCode = defaultNewLineCode
    if ( present(newLineCode) ) myNewLineCode = newLineCode

    my_adv = my_adv .and. ( severity >= MLSMessageConfig%skipMessageThr )
    if ( (.not. MLSMessageConfig%suppressDebugs).OR. &
         & (severity /= MLSMSG_Debug) ) then
      l1 = 0
      do
        l2 = index( Message(l1+1:),myNewLineCode )
        if ( l2 == 0 ) l2 = index( Message(l1+1:),myNewLineCode )
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
          call PrintItOut( line, severity, line_len, &
            & alreadyLogged=alreadyLogged )
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
      endif
      call exit_with_status ( 1  )
    end if
    if ( MustRestoreAdvOpts ) then
      AdvancedOptions = DefaultAdvancedOptions
      MustRestoreAdvOpts = .false.
    endif
  end subroutine myMessage

  ! -----------------------------------------------------  Pr_Blanks  -----
  subroutine Pr_Blanks ( N_BLANKS, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: Advance
    character(len=*), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=3) :: ADV
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    character(len=len(BlankSpace)) :: b
    integer :: I    ! Blanks to write in next WRITE statement
    integer :: N    ! Blanks remaining to write
    character(len=4) :: MY_ADV
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
    if ( MustRestoreAdvOpts ) then
      AdvancedOptions = DefaultAdvancedOptions
      MustRestoreAdvOpts = .false.
    endif
  end subroutine Pr_Blanks

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

end module Output_M

! $Log$
! Revision 2.147  2019/10/01 23:40:52  vsnyder
! Add Trim optional argument to floating-point output
!
! Revision 2.146  2019/08/01 23:42:33  pwagner
! Added SetAdvancedOption, new components to OutputOptions, numerous other changes
!
! Revision 2.145  2019/07/22 22:12:43  pwagner
! Can now setTruthPattern to something other than T and F
!
! Revision 2.144  2019/07/17 20:16:47  pwagner
! Light housekeping
!
! Revision 2.143  2019/04/09 20:30:40  pwagner
! Moved some procedures from MLSStrings to new MLSStrings_0
!
! Revision 2.142  2019/03/18 22:05:12  pwagner
! Dont print again if alreadylogged
!
! Revision 2.141  2019/01/24 18:38:05  pwagner
! Reorganized modules that print to simplify toolkit-free builds
!
! Revision 2.140  2018/10/25 23:25:20  pwagner
! Uses Pause from Io_Stuff
!
! Revision 2.139  2018/10/17 23:03:10  pwagner
! advance=.. can be used to make program pause, e.g. for debugging
!
! Revision 2.138  2018/09/13 20:18:20  pwagner
! Now gets PrUnits from similarly-named units of PrintIt_m; PrUnitname, too
!
! Revision 2.137  2018/05/11 20:33:05  pwagner
! make MAXOUTPUTLINESLEN public
!
! Revision 2.136  2018/04/05 16:55:05  pwagner
! Corrected comments; hopefully made them clearer, too
!
! Revision 2.135  2017/12/22 00:23:23  pwagner
! Move some items from output to new patternOptions; add flushOutputLines
!
! Revision 2.134  2017/11/30 20:50:13  pwagner
! RestoreSettings may now restore all or just some
!
! Revision 2.133  2017/11/15 00:00:17  pwagner
! Avoid adding unwanted blank lines when stamping
!
! Revision 2.132  2017/10/03 21:42:54  pwagner
! Added printOutputStatus; improved comments showing patterns
!
! Revision 2.131  2017/09/29 00:19:01  pwagner
! Added setOutputStatus
!
! Revision 2.130  2017/09/07 20:58:29  pwagner
! Added printOutputStatus
!
! Revision 2.129  2017/07/31 23:01:22  pwagner
! NewLine can be asked not to make a blank line
!
! Revision 2.128  2017/01/25 17:13:44  pwagner
! Output logicals so they line up with integers
!
! Revision 2.127  2016/10/18 17:46:28  pwagner
! Added advancedOptions; may insert extra args in advance='..'
!
! Revision 2.126  2016/09/22 22:51:48  pwagner
! optional format arg now in generic output api
!
! Revision 2.125  2015/09/24 18:52:28  pwagner
! Expanded special Fill patterns; allow user to set own special pattern '0'
!
! Revision 2.124  2015/08/26 23:25:11  pwagner
! Added Beep and advance='stderr' to print to stderr
!
! Revision 2.123  2015/08/25 18:36:56  vsnyder
! Add a FLUSH statement to FlushOutputLines
!
! Revision 2.122  2015/07/14 23:26:54  pwagner
! New advance_after_each arg to output of char arry; note that Reverting output now appears in old unit, e.g. stdout
!
! Revision 2.121  2015/05/18 17:40:03  pwagner
! Made Advance_is_yes_or_no public; reordered where to print
!
! Revision 2.120  2015/03/06 21:37:42  pwagner
! "%n" in output string could trigger unintended new line; fixed
!
! Revision 2.119  2015/02/13 00:16:24  pwagner
! Reordered tests in Output_Char_NoCR_INDENTED more understandably, we hope
!
! Revision 2.118  2015/02/10 00:59:36  pwagner
! Repaired error introduced by last change
!
! Revision 2.117  2015/02/06 00:45:54  pwagner
! Can now print to virtual page indented w.r.t. physical page
!
! Revision 2.116  2015/01/12 22:20:55  pwagner
! swichOutput can switch to 'stdout'
!
! Revision 2.115  2014/09/05 00:28:05  vsnyder
! Better handling of literal output unit
!
! Revision 2.114  2014/08/06 19:26:49  pwagner
! Bugfixes plus one workaround for an ifort v13 bug
!
! Revision 2.113  2014/08/05 18:23:24  pwagner
! prUnitLiteral field lets prUnit go negative
!
! Revision 2.112  2014/07/21 20:56:47  pwagner
! Should not bomb so easily
!
! Revision 2.111  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.110  2014/01/09 00:22:18  pwagner
! Split output module procedures between it and new highOutput
!
! Revision 2.109  2013/11/21 21:21:41  pwagner
! Wrap lines at the right-hand border when outputting named arrays with borders
!
! Revision 2.108  2013/11/04 22:53:51  pwagner
! Added beVerbose, letsDebug
!
! Revision 2.107  2013/09/12 01:56:50  vsnyder
! Change f6.1 format to f8.1 format in DumpSize_Double
!
! Revision 2.106  2013/08/28 00:35:39  pwagner
! Moved more stuff from MLSMessage down to PrintIt module
!
! Revision 2.105  2013/08/23 02:51:04  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 2.104  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.103  2013/07/18 22:34:31  pwagner
! Avoid double-printing when prUnit > 0
!
! Revision 2.102  2013/07/13 00:01:09  vsnyder
! Remove old comments about how unbuffering was done
!
! Revision 2.101  2013/06/28 17:53:44  pwagner
! logParent, parentName added to show who called output module
!
! Revision 2.100  2013/06/14 01:26:11  vsnyder
! Use FLUSH to unbuffer output
!
! Revision 2.99  2013/04/17 00:03:44  pwagner
! Removed LINE_WIDTH; comments note possible solution to unwieldy module
!
! Revision 2.98  2013/02/04 21:57:06  pwagner
! Fixed bug sending invalidPRUnit output to stderr
!
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
! Fixed bug that added space before newlines; simplified Output_Char
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
! Simplify by using Output_Char internally
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
! Added optional "from_where" argument to "Output_Char"
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
