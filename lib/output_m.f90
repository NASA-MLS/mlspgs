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

  use dates_module, only:  reformatDate, reformatTime
  use MLSMessageModule, only: MLSMessage, MLSMSG_Info, MLSMSG_Error
  use MLSStrings, only: lowerCase
  implicit NONE
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! MLSMSG_Level             MLSMessage level if so logged
! outputOptions            where to send output and how to format it
! PrUnit                   How to direct output:
!                          -2 :: logged via MLSMessage
!                          -1 :: print stdout
!                          < -2 :: both
!                          > 0 :: print to Fortran unit number PrUnit
! skipMLSMSGLogging        whether to skip MLSMessage by default
! stampOptions             whether and how to stamp each output to stdout
! timeStampStyle           'pre' (at linestart) or 'post' (at end-of-line)
!                            (only for timeStamp)

!     (subroutines and functions)
! blanks                   "print" specified number of blanks [or fill chars]
! dumpsize                 print a nicely-formatted memory size 
! getStamp                 get stamp being added to every output
! newline                  print a newline
! output                   print argument
! output_date_and_time     print nicely formatted date and time
! output_name_v_pair       print nicely formatted name and value
! revertOutput             revert output to file used before switchOutput
!                           if you will revert, keepOldUnitOpen when switching
! resumeOutput             resume output
! setStamp                 set stamp to be added to every output
! suspendOutput            suspend output
! switchOutput             switch output to a new named file
! timestamp                print argument with a timestamp
!                            (both stdout and logged output)
! === (end of toc) ===

! === (start of api) ===
! blanks ( int n_blanks, [char fillChar], [char* advance] )
! DumpSize ( n, [char* advance], [units] )
!       where n can be an int or a real, and 
!       units is a scalar of the same type, if present
! getStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! NewLine
! output ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank] )
! output ( char* chars(:), [char* advance],
!          [char* insteadOfBlank] )
! output ( value, [char* format], [char* advance],
!          [char* Before], [char* After] )
!       where value can be any numerical type, either scalar or 1-d array
! output_date_and_time ( [log date], [log time], [char* from_where], 
!          [char* msg], [char* dateFormat], [char* timeFormat], [char* advance] )
! output_name_v_pair ( char* name, value, [char* advance],
!          [char colon], [char fillChar], [char* Before], [char* After], 
!          [integer tabn], [integer tabc], [integer taba], log dont_stamp] )
! resumeOutput
! revertOutput
! setStamp ( [char* textCode], [log post], [int interval],
!          [log showTime], [char* dateFormat], [char* timeFormat] )
! suspendOutput
! switchOutput( char* filename, [int unit] )
! timeStamp ( char* chars, [char* advance], [char* from_where], 
!          [log dont_log], [char* log_chars], [char* insteadOfBlank],
!          [char*8 style] )
! timeStamp ( int int, [int places], [char* advance],
!          [log fill], [char* format], [char* Before], [char* After],
!          [char*8 style] )
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
!  integer, save, public :: PRUNIT = -1   ! Unit for output.  
                                         ! -1 means "printer" unit, *
                                         ! -2 means MLSMessage if 
                                         ! < -2, both printer and MLSMSG
                                         ! > 0, actual unit number
                                         !
  integer, save, private :: OLDUNIT = -1 ! Previous Unit for output.
  logical, save, private :: OLDUNITSTILLOPEN = .TRUE.

  public :: BLANKS, DUMPSIZE, GETSTAMP, NEWLINE, &
    & OUTPUT, OUTPUT_DATE_AND_TIME, OUTPUT_NAME_V_PAIR, &
    & RESUMEOUTPUT, revertOutput, &
    & SETSTAMP, SUSPENDOUTPUT, switchOutput, TIMESTAMP

  public :: outputOptions_T
  public :: stampOptions_T

  interface DUMPSIZE
    module procedure DUMPSIZE_INTEGER, DUMPSIZE_REAL
  end interface

  interface OUTPUT
    module procedure output_char, output_char_array, output_complex
    module procedure output_dcomplex, output_double
    module procedure output_integer, output_integer_array
    module procedure output_logical, output_logical_array
    module procedure output_single, output_double_array, output_single_array
    module procedure output_string
  end interface

  interface OUTPUT_NAME_V_PAIR
    module procedure output_nvp_character
    module procedure output_nvp_complex
    module procedure output_nvp_double
    module procedure output_nvp_integer
    module procedure output_nvp_int_array
    module procedure output_nvp_logical
    module procedure output_nvp_log_array
    module procedure output_nvp_single
  end interface

  interface TIMESTAMP
    module procedure timestamp_char, timestamp_integer
  end interface
  
  ! This is the type for configuring how to automatically format
  ! lines and whether they should be sent to stdout or elsewhere
  type outputOptions_T
    integer :: PRUNIT = -1
    integer :: MLSMSG_Level = MLSMSG_Info ! What level if logging
    logical :: SKIPMLSMSGLOGGING = .false.
    logical :: usePatternedBlanks = .true. ! Use patterns for special fillChars
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
  end type
  
  type(outputOptions_T), public, save :: outputOptions

  ! This is the type for configuring whether and how to automatically stamp
  ! lines sent to stdout
  ! If interval > 1 only a fraction of lines will be stamped
  ! If interval > 10, the stamps will not be in-line, instead
  ! they will appear alone as page headers
  ! (As an alternative, use timeStamp to stamp only individual lines)
  type stampOptions_T
    logical :: post = .true.      ! Put stamp at end of line?
    logical :: showTime = .false. ! Don't show date or time unless TRUE
    character(len=24) :: textCode = ' '
    ! Don't show date unless dateFormat is non-blank
    character(len=16) :: dateFormat = ' '
    character(len=16) :: timeFormat = 'hh:mm'
    integer :: interval = 1 ! 1 means stamp every line; n means every nth line
    character(len=8) :: TIMESTAMPSTYLE = 'post' ! 'pre' or 'post'
  end type
  
  type(stampOptions_T), public, save :: stampOptions ! Could leave this private

  ! integer, save, public :: MLSMSG_Level = MLSMSG_Info ! What level if logging
  logical, save, private:: SILENTRUNNING = .false. ! Suspend all further output
  ! logical, save, public :: SKIPMLSMSGLOGGING = .false.
  ! character(len=8), save, public :: TIMESTAMPSTYLE = 'post' ! 'pre' or 'post'
  ! Private parameters
  integer, save, private :: ATCOLUMNNUMBER = 1  ! Where we'll print next
  logical, save, private :: ATLINESTART = .true.  ! Whether to stamp if notpost
  integer, save, private :: LINESSINCELASTSTAMP = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -----------------------------------------------------  BLANKS  -----
  subroutine BLANKS ( N_BLANKS, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=1), intent(in), optional :: FILLCHAR  ! default is ' '
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    character(len=3) :: ADV
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    character(len=len(BlankSpace)) :: b
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
        endif
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
          endif
        else
          call pr_blanks ( 1, advance='no' )
          numSoFar = 1
        endif
        do i=1, ntimes
          call output ( pattern(2:patternLength+1), advance='no' )
          ! if ( xtraBlanks > 0 ) call pr_blanks ( xtraBlanks, advance='no' )
          numSoFar = numSoFar + patternLength
        enddo
        theRest = n_blanks - numSoFar
        if ( theRest > 0 ) call pr_blanks ( theRest, advance=advance, dont_stamp=dont_stamp )
        return
      endif
    endif
    call pr_blanks ( n_blanks, fillChar=fillChar, advance=advance, dont_stamp=dont_stamp )
  end subroutine BLANKS

  ! -----------------------------------------------------  PR_BLANKS  -----
  subroutine PR_BLANKS ( N_BLANKS, FILLCHAR, ADVANCE, DONT_STAMP )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=1), intent(in), optional :: FILLCHAR  ! default is ' '
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
    if ( n_blanks < 1 .and. my_adv == 'yes' ) then
      ! This does not yet prevent us from printing a space before newline
      ! but someday we may achieve that goal
      call output ( '', advance='yes', dont_stamp=dont_stamp )
      return
    endif
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
      call output ( b(:i), advance=adv )
      if ( n < 1 ) exit   ! was if n == 0, but this should be safer
    end do
  end subroutine PR_BLANKS

  ! --------------------------------------------- DumpSize_integer -----
  subroutine DumpSize_integer ( n, advance, units )
    integer, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    integer, intent(in), optional :: units ! E.g., 1024 for kB
    ! Executable
    if ( present(units) ) then
      call dumpSize ( n*1.0, advance=advance, units=units*1.0 )
    else
      call dumpSize ( n*1.0, advance=advance )
    endif
  end subroutine DumpSize_integer

  ! ------------------------------------------------ DumpSize_real -----
  subroutine DumpSize_real ( n, advance, units )
    real, intent(in) :: N
    character(len=*), intent(in), optional :: ADVANCE
    real, intent(in), optional :: units
    ! Local parameters
    real, parameter :: KB = 1024.0
    real, parameter :: MB = KB * 1024.0
    real, parameter :: GB = MB * 1024.0
    real, parameter :: TB = GB * 1024.0
    real            :: myUnits
    ! Make a 'nice' output
    myUnits = 1.0
    if ( present(units) ) myUnits = units
    if ( myUnits == 0.0 ) then
      call output ( n, format='(e12.1)' )
      call output ( ' (illegal units)', advance=advance )
      return
    endif
    if ( n < kb/myUnits ) then
      call output ( n*myUnits, format='(f6.1)' )
      call output ( ' bytes', advance=advance )
    else if ( n < Mb/myUnits ) then
      call output ( n*myUnits/kb, format='(f6.1)' )
      call output ( ' kb', advance=advance )
    else if ( n < Gb/myUnits ) then
      call output ( n*myUnits/Mb, format='(f6.1)' )
      call output ( ' Mb', advance=advance )
    else if ( n < Tb/myUnits ) then
      call output ( n*myUnits/Gb, format='(f6.1)' )
      call output ( ' Gb', advance=advance )
    else
      call output ( n*myUnits/Tb, format='(f6.1)' )
      call output ( ' Tb', advance=advance )
    end if
  end subroutine DumpSize_real

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

  ! ----------------------------------------------------  NewLine  -----
  subroutine NewLine
    call output ( '', advance='yes' )
  end subroutine NewLine

  ! ------------------------------------------------  OUTPUT_CHAR  -----
  subroutine OUTPUT_CHAR ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, DONT_STAMP )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    logical, intent(in), optional          :: DONT_STAMP ! Prevent double-stamping
    !
    character(len=max(16,len(chars)+1)) :: my_chars
    character(len=len(chars)+64) :: stamped_chars ! What to print to stdout
    character(len=max(16,len(chars)+1)) :: the_chars
    logical :: my_dont_log
    logical :: my_dont_stamp
    integer :: n_chars
    character(len=3) :: MY_ADV
    integer :: n_stamp ! How much of stamped_chars to print
    logical :: stamped
    logical :: stamp_header
    !
    if ( SILENTRUNNING ) return
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_stamp = .false.
    if ( present(dont_stamp) ) my_dont_stamp = dont_stamp
    my_dont_stamp = ( my_dont_stamp .or. &
      & linesSinceLastStamp < (stampOptions%interval - 1) )
    stamped = .false.
    stamp_header = .false.
    stamped_chars = chars
    ! Do we need to stamp this line? If so, at beginning or at end?
    if ( my_dont_stamp ) then
    elseif ( stampoptions%interval > 9 ) then
      stamp_header = (my_adv == 'yes')
    elseif ( stampoptions%post ) then
      if ( my_adv == 'yes' ) then
        stamped_chars = stamp(chars)
        stamped = .true.
      endif
    elseif( ATLINESTART ) then
      stamped_chars = stamp(chars)
      stamped = .true.
    endif

    n_chars = max(len(chars), 1)
    my_dont_log = outputOptions%skipmlsmsglogging ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    n_stamp = len_trim(stamped_chars)
    if ( my_adv == 'no' ) n_stamp = n_stamp + len(chars) - len_trim(chars)
    ! Special case: if chars is blank (chars are blank?)
    ! we'll want to print anyway
    if ( len_trim(chars) < 1 ) n_stamp = max(n_stamp, 1)
    if ( (outputOptions%prunit == -1 .or. outputOptions%prunit < -2) .and. n_stamp > 0 ) &
      & write ( *, '(a)', advance=my_adv ) stamped_chars(1:n_stamp)
    if ( outputOptions%prunit < -1 .and. .not. my_dont_log  ) then
      the_chars = chars // ' '
      n_chars = max(len(chars), 1)
      if ( present(log_chars) ) then
        n_chars = len_trim(log_chars) + 1
        the_chars = log_chars(:n_chars-1) // ' '
      end if
      if ( the_chars == ' ' .and. present(insteadofblank)  ) then
        my_chars = trim(insteadofblank) // ' '
        n_chars = max(len(insteadofblank), 1)
      else
        my_chars = the_chars
      end if
      if ( my_adv == 'no' ) n_chars = n_chars+1
      n_chars = min(n_chars, len(my_chars))
      if ( present(from_where)  ) then
        call MLSMessage ( outputOptions%MLSMSG_Level, from_where, my_chars(1:n_chars), &
          & advance=my_adv )
      else
        call MLSMessage ( outputOptions%MLSMSG_Level, ModuleName, my_chars(1:n_chars), &
          & advance=my_adv )
      end if
    end if
    
    if ( outputOptions%prunit < 0  ) then
      ! Already logged; no output to stdout
    else if ( stamped_chars == ' ' .and. present(insteadofblank)  ) then
      write ( outputOptions%prunit, '(a)', advance=my_adv ) trim(insteadofblank)
    else
      write ( outputOptions%prunit, '(a)', advance=my_adv ) stamped_chars(1:n_stamp)
    end if
    atLineStart = (my_adv == 'yes')
    if ( atLineStart ) then
      if ( stamp_header ) then
        ! Time to add stamp as a page header on its own line
        stamped_chars = stamp(' ')
        stamped = .true.
        if ( outputOptions%prunit == -1 .or. outputOptions%prunit < -2 ) then
          write ( *, '(a)', advance='yes' ) trim(stamped_chars)
        elseif( outputOptions%prunit > 0 ) then
          write ( outputOptions%prunit, '(a)', advance='yes' ) trim(stamped_chars)
        endif
      endif
      if ( stamped ) then
        linesSinceLastStamp = 0
      else
        linesSinceLastStamp = linesSinceLastStamp + 1
      endif
    endif
    atColumnNumber = atColumnNumber + n_chars
    if ( atLineStart ) atColumnNumber = 1
  end subroutine OUTPUT_CHAR

  ! ------------------------------------------  OUTPUT_CHAR_ARRAY  -----
  subroutine OUTPUT_CHAR_ARRAY ( CHARS, ADVANCE, INSTEADOFBLANK )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    integer :: I ! loop inductor
    do i = 1, size(chars)
      call output ( chars(i), insteadofblank=insteadofblank )
    end do
    if ( present(advance)  ) then
      call output ( '', advance=advance )
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
    if ( present(before) ) call output ( before, dont_log = .true. )
    if ( present(after)  ) then
      call output ( trim(line), dont_log = .true. )
      call output ( after, advance=advance, dont_log = .true., &
        & dont_stamp=dont_stamp )
    else
      call output ( trim(line), advance=advance, dont_log = .true., &
        & dont_stamp=dont_stamp )
    end if
  end subroutine OUTPUT_COMPLEX

  ! ---------------------------------------  OUTPUT_DATE_AND_TIME  -----
  subroutine OUTPUT_DATE_AND_TIME ( date, time, &
    & from_where, msg, dateFormat, timeFormat, advance )
    ! Output nicely-formatted date, time, and extra message
    ! We'll assume we won't want this line stamped with date and time
    ! (for fear of being redundant, which we fear)
    logical, intent(in), optional :: date ! output date as character string
    logical, intent(in), optional :: time ! output time as character string
    character(len=*), intent(in), optional :: FROM_WHERE
    character(len=*), intent(in), optional :: MSG
    character(len=*), intent(in), optional :: DATEFORMAT
    character(len=*), intent(in), optional :: TIMEFORMAT
    character(len=*), intent(in), optional :: ADVANCE
    !
    character(len=16) :: dateString
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    character(len=16) :: timeString
    logical :: myDate
    logical :: myTime
    character(len=3) :: MY_ADV
    !
    myDate = .true.
    if ( present(date) ) myDate = date
    myTime = .true.
    if ( present(time) ) myTime = time
    if ( .not. (myDate .or. myTime) ) return ! Why call if won't print?
    my_adv = 'no'
    if ( .not. present(msg) ) then
      my_adv = 'yes'
      if ( present(advance) ) my_adv = advance
    endif
    call date_and_time ( date=dateString, time=timeString )
    dateString = reFormatDate(trim(dateString), toForm=dateFormat)
    timeString = reFormatTime(trim(timeString), timeFormat)
    if ( myDate .and. myTime ) then
      call output ( trim(dateString), from_where=from_where, advance='no', &
        & DONT_STAMP=DONT_STAMP )
      call blanks(3)
      call output ( trim(timeString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    else if ( myDate ) then
      call output ( trim(dateString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    else if ( myTime ) then
      call output ( trim(TimeString), from_where=from_where, advance=my_adv, &
        & DONT_STAMP=DONT_STAMP )
    end if
    if ( .not. present(msg) ) return
    my_adv = 'yes'
    if ( present(advance) ) my_adv = advance
    call blanks ( 3 )
    call output ( trim(msg), from_where=from_where, advance=my_adv, &
      & DONT_STAMP=DONT_STAMP )
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
    if ( present(before) ) call output ( before, dont_log = .true. )
    if ( present(after)  ) then
      call output ( trim(line), dont_log = .true. )
      call output ( after, advance=advance, dont_log = .true. )
    else
      call output ( trim(line), advance=advance, dont_log = .true. )
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
    character(len=30) :: LINE, LOG_CHARS

    line = ' '
    if ( .not. present(Format)  ) then
   ! No optional formats: use default char-by-char accretion
      write ( line, * ) value
      if ( scan(line,'123456789') == 0  ) then
        line = '0'
      else
        i = index(line,'.')
        j = scan(line(i:),'DdEe ') + i - 1
        if ( i /= 0  ) then
          if ( j == i ) j = len(line)
          i = i + 1
          k = j
          do while ( j > i )
            j = j - 1
            if ( line(j:j) /= '0' .and. line(j:j) /= ' ') exit
          end do
          line(j+1:) = line(k:)
        end if
        line = adjustl(line)
      end if
      k = len_trim(line)
    ! Use one or both optional formats
    else
      line = ' '
      write ( line, Format ) value
      k = nCharsinFormat(Format)
      if ( k==0 ) k = len_trim(line)
    end if

    log_chars = line
    if ( present(LogFormat)  ) then
      write ( log_chars, LogFormat ) value
    end if
    if ( present(before) ) call output ( before )
    if ( present(after)  ) then
      call output ( line(:k), log_chars=log_chars )
      call output ( after, advance=advance, dont_stamp=dont_stamp )
    else
      call output ( line(:k), advance=advance, log_chars=log_chars, &
        & dont_stamp=dont_stamp )
    end if

  end subroutine OUTPUT_DOUBLE

  ! ----------------------------------------  OUTPUT_DOUBLE_ARRAY  -----
  subroutine OUTPUT_DOUBLE_ARRAY ( values, FORMAT, LogFormat, ADVANCE )
  ! Output double-precision values to PRUNIT.
    double precision, intent(in) :: values(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    character(len=3) :: MY_ADV
    integer :: I ! loop inductor
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      call blanks ( 3, advance='no' )
    end do
    if ( present(advance)  ) then
      if ( outputOptions%prunit == -1 .or. outputOptions%prunit < -2 ) &
        & write ( *, '(a)', advance=my_adv )
      if ( outputOptions%prunit < -1 ) &
        & call MLSMessage ( outputOptions%MLSMSG_Level, ModuleName, '', advance=my_adv )
      if ( outputOptions%prunit >= 0 ) &
        & write ( outputOptions%prunit, '(a)', advance=my_adv )
    end if
  end subroutine OUTPUT_DOUBLE_ARRAY

  ! ---------------------------------------------  OUTPUT_INTEGER  -----
  subroutine OUTPUT_INTEGER ( INT, PLACES, ADVANCE, FILL, FORMAT, &
    & Before, After, DONT_STAMP )
  ! Output INT to PRUNIT using at most PLACES (default zero) places
  ! If 'fill' is present and true, fill leading blanks with zeroes (only
  ! makes sense if 'places' is specified).
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
    my_places = 0
    if ( present(places)  ) then; my_places = places; end if
    my_fill = .false.
    if ( present(places) .and. present(fill) ) my_fill = fill
    if ( present(format)  ) then
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
    if ( present(before) ) call output ( before )
    if ( present(after)  ) then
      call output ( line(i:j), DONT_STAMP=DONT_STAMP )
      call output ( after, advance=advance, DONT_STAMP=DONT_STAMP )
    else
      call output ( line(i:j), advance=advance, DONT_STAMP=DONT_STAMP )
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
      call blanks ( 3, advance='no' )
    end do
    if ( present(advance)  ) then
      ! Why all this needless complication?
      ! if ( outputOptions%prunit == -1 .or. outputOptions%prunit < -2 ) &
      !   & write ( *, '(a)', advance=my_adv )
      ! if ( outputOptions%prunit < -1 ) &
      !   & call MLSMessage ( outputOptions%MLSMSG_Level, ModuleName, '', advance=my_adv )
      ! if ( outputOptions%prunit >= 0 ) &
      !   & write ( outputOptions%prunit, '(a)', advance=my_adv )
      
      ! Just do it the easy way
      call output ( '', advance=advance, DONT_STAMP=DONT_STAMP )
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
    if ( present(before) ) call output ( before, DONT_STAMP=DONT_STAMP )
    call output ( line, advance=advance, DONT_STAMP=DONT_STAMP )
  end subroutine OUTPUT_LOGICAL

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL_ARRAY ( logs, ADVANCE, BEFORE, DONT_STAMP )
  ! Output LOG to PRUNIT using at most PLACES (default zero) places
    logical, dimension(:), intent(in) :: logs
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    logical, optional, intent(in) :: DONT_STAMP
    integer :: I ! loop inductor
    if ( present(before) ) call output ( before, DONT_STAMP=DONT_STAMP )
    do i = 1, size(logs)
      call output ( logs(i), advance='no' )
      call blanks ( 3, advance='no' )
    end do
    if ( present(advance) ) call output ( '', advance=advance, DONT_STAMP=DONT_STAMP )
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
    character(len=30) :: LINE, LOG_CHARS

    line = ' '
    if ( .not. present(Format)  ) then
   ! No optional formats: use default char-by-char accretion
      write ( line, * ) value
      if ( scan(line,'123456789') == 0  ) then
        line = '0'
      else
        i = index(line,'.')
        j = scan(line(i:),'DdEe ') + i - 1
        if ( i /= 0  ) then
          if ( j == i ) j = len(line)
          i = i + 1
          k = j
          do while ( j > i )
            j = j - 1
            if ( line(j:j) /= '0' .and. line(j:j) /= ' ') exit
          end do
          line(j+1:) = line(k:)
        end if
        line = adjustl(line)
      end if
      k = len_trim(line)
    ! Use one or both optional formats
    else
      line = ' '
      write ( line, Format ) value
      k = nCharsinFormat(Format)
      if ( k==0 ) k = len_trim(line)
    end if

    log_chars = line
    if ( present(LogFormat)  ) then
      write ( log_chars, LogFormat ) value
    end if
    if ( present(before) ) call output ( before, DONT_STAMP=DONT_STAMP )
    if ( present(after)  ) then
      call output ( line(:k), log_chars=log_chars, DONT_STAMP=DONT_STAMP )
      call output ( after, advance=advance, DONT_STAMP=DONT_STAMP )
    else
      call output ( line(:k), advance=advance, log_chars=log_chars, &
        & DONT_STAMP=DONT_STAMP )
    end if
  end subroutine OUTPUT_SINGLE

  ! ----------------------------------------  OUTPUT_SINGLE_ARRAY  -----
  subroutine OUTPUT_SINGLE_ARRAY ( values, FORMAT, LogFormat, ADVANCE )
  ! Output single-precision values to PRUNIT.
    real, intent(in) :: values(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    integer :: I ! loop inductor
    character(len=3) :: MY_ADV
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(values)
      call output ( values(i), advance='no', format=format, logFormat=logFormat )
      call blanks ( 3, advance='no' )
    end do
    if ( present(advance)  ) then
      if ( outputOptions%prunit == -1 .or. outputOptions%prunit < -2 ) &
        & write ( *, '(a)', advance=my_adv )
      if ( outputOptions%prunit < -1 ) &
        & call MLSMessage ( outputOptions%MLSMSG_Level, ModuleName, '', advance=my_adv )
      if ( outputOptions%prunit >= 0 ) &
        & write ( outputOptions%prunit, '(a)', advance=my_adv )
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
      call output ( '', advance )
    else
      call output ( string(:n_chars), advance, from_where, dont_log, log_chars )
    end if
  end subroutine OUTPUT_STRING

  ! ----------------------------------------------  OUTPUT_NAME_V_PAIR  -----
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
      call output( 'Unable to Revert output--old unit not open', advance='yes' )
      return
    endif
    call output( 'Reverting output to unit: ', advance='no' )
    call output( OLDUNIT, advance='yes' )
    if ( outputOptions%prunit > 0 ) then
      inquire( unit=outputOptions%prunit, opened=itsOpen )
      if ( itsOpen ) then
        close(outputOptions%prunit)
      endif
    endif
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
      call output('Switching further output to: ', advance='no')
      call output(trim(filename), advance='yes')
      call output('using unit number: ', advance='no')
      call output(unit, advance='yes')
      switchUnit = unit
    else
      call output('Switching further output to: ', advance='no')
      call output(trim(filename), advance='yes')
      ! if ( PRUNIT > 0 ) switchUnit = PRUNIT
    endif
    if ( outputOptions%prunit > 0 .and. .not. dontCloseOldUnit ) then
      close(outputOptions%prunit)
      OLDUNIT = outputOptions%prunit
      OLDUNITSTILLOPEN = .false.
    endif
    if ( NeedToAppend ) then
      open( unit=switchUnit, file=filename, status='old', position='append' )
    else
      open( unit=switchUnit, file=filename, status='replace' )
    endif
    outputOptions%prunit = SwitchUnit
    NeedToAppend = .true.
  end subroutine switchOutput

  ! ------------------------------------------------  timeStamp_char  -----
  subroutine timeStamp_char ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK, STYLE)
  ! time-stamp output:
  ! Either in style pre or post
  ! (pre) '(HH:MM:SS) chars'
  ! (post) 'chars (HH:MM:SS)'
  ! Note that in pre-style, the time will be printed only if ATLINESTART true
  ! in post-style, the time will be printed only if MY_ADV is 'yes'
  
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    character(len=8), intent(in), optional :: STYLE ! pre or post
    !
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    logical :: my_dont_log
    character(len=8) :: my_style
    character(len=3) :: MY_ADV
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_log = outputOptions%skipmlsmsglogging ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    my_style = stampOptions%Timestampstyle
    if ( present(style) ) my_style = lowercase(style)
    if ( my_style == 'post' ) then
      call output_char( CHARS, &
        & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
        & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK, DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_char(' (', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='yes', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      endif
    else
      if ( ATLINESTART ) then
        call output_char('(', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='no', DONT_LOG=DONT_LOG, DONT_STAMP=DONT_STAMP)
      endif
      call output_char( CHARS, &
        & ADVANCE, FROM_WHERE, DONT_LOG, &
        & LOG_CHARS, INSTEADOFBLANK, DONT_STAMP=DONT_STAMP )
    endif
  end subroutine timeStamp_char

  ! ------------------------------------------------  timeStamp_integer  -----
  subroutine timeStamp_integer ( INT, &
    & PLACES, ADVANCE, FILL, FORMAT, Before, After, style)
  ! time-stamp output:
  ! Either in style pre or post
  ! (pre) '(HH:MM:SS) int'
  ! (post) 'int (HH:MM:SS)'
  ! Note that in pre-style, the time will be printed only if ATLINESTART true
  ! in post-style, the time will be printed only if MY_ADV is 'yes'
  
    integer, intent(in) :: INT
    integer, intent(in), optional :: PLACES
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: FILL
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=*), intent(in), optional :: style ! pre or post
    !
    logical, parameter :: DONT_STAMP = .true. ! Don't double-stamp
    character(len=8) :: my_style
    character(len=3) :: MY_ADV
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_style = stampOptions%Timestampstyle
    if ( present(style) ) my_style = lowercase(style)
    if ( my_style == 'post' ) then
      call output_integer( INT, PLACES, &
        & ADVANCE='no', FILL=FILL, FORMAT=FORMAT, BEFORE=BEFORE, AFTER=AFTER, &
        & DONT_STAMP=DONT_STAMP )
      if ( my_adv=='yes' ) then
        call output_char(' (', ADVANCE='no', DONT_STAMP=DONT_STAMP )
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='yes', DONT_STAMP=DONT_STAMP)
      endif
    else
      if ( ATLINESTART ) then
        call output_char('(', ADVANCE='no', DONT_STAMP=DONT_STAMP)
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='no', DONT_STAMP=DONT_STAMP)
      endif
      call output_integer( INT, PLACES, &
        & ADVANCE, FILL, FORMAT, BEFORE, AFTER, DONT_STAMP=DONT_STAMP )
    endif
  end subroutine timeStamp_integer

  ! Internal procedures
  
  ! .......................................  Advance_is_yes_or_no  .....
  function Advance_is_yes_or_no ( str ) result ( outstr )
    ! takes '[Yy]...' or '[Nn..] and returns 'yes' or 'no' respectively
    ! also does the same with '[Tt]..' and '[Ff]..'
    ! leaves all other patterns unchanged, but truncated to three
    ! characters.  Returns 'no' if the argument is absent.
    !--------Argument--------!
    character (len=*), intent(in), optional :: Str
    character (len=3) :: Outstr

    !----------Local vars----------!
    character (len=*), parameter :: yeses = 'YyTt'
    character (len=*), parameter :: nose = 'NnFf'

    if ( .not. present(str)  ) then
      outstr = 'no'
      return
    end if

    outstr = adjustl(str)
    if ( index( yeses, outstr(:1)) > 0  ) then
      outstr = 'yes'
    else if ( index( nose, outstr(:1)) > 0  ) then
      outstr = 'no'
    else
      outstr = str
    end if
  end function Advance_is_yes_or_no

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
    endif                                                                 
    nplusm = n + m                                                        
  end function nCharsinFormat

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
    endif
  end subroutine ourReplaceSubString

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
      call date_and_time ( date=dateString, time=timeString )
      dateString = reFormatDate(trim(dateString), toForm=stampOptions%dateFormat)
      timeString = reFormatTime(trim(timeString), stampOptions%timeFormat)
      if ( stampOptions%dateFormat /= ' ' ) &
      & stamp = catStrings( stamp, dateString )
      if ( stampOptions%timeFormat /= ' ' ) &
      & stamp = catStrings( stamp, timeString )
    endif
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
      endif
    end function catStrings
    
  end function stamp 

  ! ..............................................  not_used_here  .....
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module OUTPUT_M

! $Log$
! Revision 2.51  2006/06/27 23:58:01  pwagner
! name_v_value works much better
!
! Revision 2.50  2006/06/24 23:05:07  pwagner
! Added output_name_v_pair, special blank fills
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
