! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module OUTPUT_M

  use MLSMessageModule, only: MLSMessage, MLSMSG_Info, MLSMSG_Error
  use MLSStrings, only:  lowercase, reformatDate, reformatTime
  implicit NONE
  private

  integer, save, public :: LINE_WIDTH = 120 ! Not used here, but a convenient
                                        ! place to store it
  integer, save, public :: PRUNIT = -1  ! Unit for output.  "printer" unit, *
                                        ! if -1, MLSMessage if -2, both
                                        ! printer and MLSMSG if < -2.

  public :: BLANKS, NEWLINE, OUTPUT, OUTPUT_DATE_AND_TIME, TIMESTAMP
  interface OUTPUT
    module procedure output_char, output_char_array, output_complex
    module procedure output_dcomplex, output_double
    module procedure output_integer, output_integer_array, output_logical
    module procedure output_single, output_double_array, output_single_array
    module procedure output_string
  end interface
  interface TIMESTAMP
    module procedure timestamp_char, timestamp_integer
  end interface

  integer, save, public :: MLSMSG_Level = MLSMSG_Info
  logical, save, public :: SKIPMLSMSGLOGGING = .false.
  character(len=8), save, public :: TIMESTAMPSTYLE = 'post' ! 'pre' or 'post'
  logical, save, private :: ATLINESTART = .true.  ! Used by timeStamp if notpost

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

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

  ! -----------------------------------------------------  BLANKS  -----
  subroutine BLANKS ( N_BLANKS, FILLCHAR, ADVANCE )
  ! Output N_BLANKS blanks to PRUNIT.
  ! (or optionally that many copies of fillChar)
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=1), intent(in), optional :: FILLCHAR  ! default is ' '
    character(len=3) :: ADV
    character(len=*), parameter :: BLANKSPACE = &
    '                                                                    '
    character(len=len(BlankSpace)) :: b
    integer :: I    ! Blanks to write in next WRITE statement
    integer :: N    ! Blanks remaining to write
    character(len=3) :: MY_ADV
    ! Executable
    my_adv = Advance_is_yes_or_no(advance)
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
  end subroutine BLANKS

  ! ----------------------------------------------------  NewLine  -----
  subroutine NewLine
    call output ( '', advance='yes' )
  end subroutine NewLine

  ! ------------------------------------------------  OUTPUT_CHAR  -----
  subroutine OUTPUT_CHAR ( CHARS, &
    & ADVANCE, FROM_WHERE, DONT_LOG, LOG_CHARS, INSTEADOFBLANK)
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FROM_WHERE
    logical, intent(in), optional          :: DONT_LOG ! Prevent double-logging
    character(len=*), intent(in), optional :: LOG_CHARS
    character(len=*), intent(in), optional :: INSTEADOFBLANK ! What to output
    character(len=max(16,len(chars)+1)) :: my_chars
    character(len=max(16,len(chars)+1)) :: the_chars
    logical :: my_dont_log
    integer :: n_chars
    character(len=3) :: MY_ADV
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_log = SKIPMLSMSGLOGGING ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    if ( prunit == -1 .or. prunit < -2 ) &
      & write ( *, '(a)', advance=my_adv ) chars
    if ( prunit < -1 .and. .not. my_dont_log  ) then
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
        call MLSMessage ( MLSMSG_Level, from_where, my_chars(1:n_chars), &
          & advance=my_adv )
      else
        call MLSMessage ( MLSMSG_Level, ModuleName, my_chars(1:n_chars), &
          & advance=my_adv )
      end if
    end if
    
    if ( prunit < 0  ) then
      ! Already logged; no output to stdout
    else if ( chars == ' ' .and. present(insteadofblank)  ) then
      write ( prunit, '(a)', advance=my_adv ) trim(insteadofblank)
    else
      write ( prunit, '(a)', advance=my_adv ) chars
    end if
    ATLINESTART = (my_adv == 'yes')
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
  subroutine OUTPUT_COMPLEX ( VALUE, Format, ADVANCE, Before, After )
    complex, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
    character(len=60) :: LINE

    if ( present(Format)  ) then
      write ( line, format ) value
    else
      write ( line, '("(",1pg15.7,",",1pg15.7,")")' ) value
    end if
    if ( present(before) ) call output ( before, dont_log = .true. )
    if ( present(after)  ) then
      call output ( trim(line), dont_log = .true. )
      call output ( after, advance=advance, dont_log = .true. )
    else
      call output ( trim(line), advance=advance, dont_log = .true. )
    end if
  end subroutine OUTPUT_COMPLEX

  ! ---------------------------------------  OUTPUT_DATE_AND_TIME  -----
  subroutine OUTPUT_DATE_AND_TIME ( date, time, &
    & from_where, msg, dateFormat, timeFormat, advance )
    logical, intent(in), optional :: date ! output date as character string
    logical, intent(in), optional :: time ! output time as character string
    character(len=*), intent(in), optional :: FROM_WHERE
    character(len=*), intent(in), optional :: MSG
    character(len=*), intent(in), optional :: DATEFORMAT
    character(len=*), intent(in), optional :: TIMEFORMAT
    character(len=*), intent(in), optional :: ADVANCE
    character(len=16) :: dateString
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
    dateString = reFormatDate(trim(dateString), dateFormat)
    timeString = reFormatTime(trim(timeString), timeFormat)
    if ( myDate .and. myTime  ) then
      call output ( trim(dateString), from_where=from_where, advance='no' )
      call blanks(3)
      call output ( trim(timeString), from_where=from_where, advance=my_adv )
    else if ( myDate  ) then
      call output ( trim(dateString), from_where=from_where, advance=my_adv )
    else if ( myTime  ) then
      call output ( trim(TimeString), from_where=from_where, advance=my_adv )
    end if
    if ( .not. present(msg) ) return
    my_adv = 'yes'
    if ( present(advance) ) my_adv = advance
    call blanks ( 3 )
    call output ( trim(msg), from_where=from_where, advance=my_adv )
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
  subroutine OUTPUT_DOUBLE ( VALUE, Format, LogFormat, ADVANCE, Before, After )
  ! Output "double" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format    ! How to print
    character(len=*), intent(in), optional :: LogFormat ! How to post to Log
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
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
      call output ( after, advance=advance )
    else
      call output ( line(:k), advance=advance, log_chars=log_chars )
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
      if ( prunit == -1 .or. prunit < -2 ) &
        & write ( *, '(a)', advance=my_adv )
      if ( prunit < -1 ) &
        & call MLSMessage ( MLSMSG_Level, ModuleName, '', advance=my_adv )
      if ( prunit >= 0 ) &
        & write ( prunit, '(a)', advance=my_adv )
    end if
  end subroutine OUTPUT_DOUBLE_ARRAY

  ! ---------------------------------------------  OUTPUT_INTEGER  -----
  subroutine OUTPUT_INTEGER ( INT, PLACES, ADVANCE, FILL, FORMAT, Before, After )
  ! Output INT to PRUNIT using at most PLACES (default zero) places
  ! If 'fill' is present and true, fill leading blanks with zeroes (only
  ! makes sense if 'places' is specified).
    integer, intent(in) :: INT
    integer, intent(in), optional :: PLACES
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: FILL
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: Before, After ! text to print
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
      call output ( line(i:j) )
      call output ( after, advance=advance )
    else
      call output ( line(i:j), advance=advance )
    end if
  end subroutine OUTPUT_INTEGER

  ! ---------------------------------------  OUTPUT_INTEGER_ARRAY  -----
  subroutine OUTPUT_INTEGER_ARRAY ( INTEGERS, ADVANCE, FORMAT )
  ! Output INTEGERS to PRUNIT.
    integer, intent(in) :: INTEGERS(:)
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: FORMAT
    character(len=3) :: MY_ADV
    integer :: I ! loop inductor
    my_adv = Advance_is_yes_or_no(advance)
    do i = 1, size(integers)
      call output ( integers(i), advance='no', format=format )
      call blanks ( 3, advance='no' )
    end do
    if ( present(advance)  ) then
      if ( prunit == -1 .or. prunit < -2 ) &
        & write ( *, '(a)', advance=my_adv )
      if ( prunit < -1 ) &
        & call MLSMessage ( MLSMSG_Level, ModuleName, '', advance=my_adv )
      if ( prunit >= 0 ) &
        & write ( prunit, '(a)', advance=my_adv )
    end if
  end subroutine OUTPUT_INTEGER_ARRAY

  ! ---------------------------------------------  OUTPUT_LOGICAL  -----
  subroutine OUTPUT_LOGICAL ( LOG, ADVANCE, BEFORE )
  ! Output LOG to PRUNIT using at most PLACES (default zero) places
    logical, intent(in) :: LOG
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    character(len=2) :: LINE
    if ( log ) then
      line=' T'
    else
      line=' F'
    end if
    if ( present(before) ) call output ( before )
    call output ( line, advance=advance )
  end subroutine OUTPUT_LOGICAL

  ! ----------------------------------------------  OUTPUT_SINGLE  -----
  subroutine OUTPUT_SINGLE ( VALUE, FORMAT, LogFormat, ADVANCE, Before, After )
  ! Output "SINGLE" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: Format  ! How to print
    character(len=*), intent(in), optional :: LogFormat     ! How to post to Log
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: Before, After ! text to print
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
      call output ( after, advance=advance )
    else
      call output ( line(:k), advance=advance, log_chars=log_chars )
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
      if ( prunit == -1 .or. prunit < -2 ) &
        & write ( *, '(a)', advance=my_adv )
      if ( prunit < -1 ) &
        & call MLSMessage ( MLSMSG_Level, ModuleName, '', advance=my_adv )
      if ( prunit >= 0 ) &
        & write ( prunit, '(a)', advance=my_adv )
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
    logical :: my_dont_log
    character(len=8) :: my_style
    integer :: n_chars
    character(len=3) :: MY_ADV
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_dont_log = SKIPMLSMSGLOGGING ! .false.
    if ( present(dont_log) ) my_dont_log = dont_log
    my_style = TIMESTAMPSTYLE
    if ( present(style) ) my_style = lowercase(style)
    if ( my_style == 'post' ) then
      call output_char( CHARS, &
        & ADVANCE='no', FROM_WHERE=FROM_WHERE, DONT_LOG=DONT_LOG, &
        & LOG_CHARS=LOG_CHARS, INSTEADOFBLANK=INSTEADOFBLANK )
      if ( my_adv=='yes' ) then
        call output_char(' (', ADVANCE='no', DONT_LOG=DONT_LOG)
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='yes', DONT_LOG=DONT_LOG)
      endif
    else
      if ( ATLINESTART ) then
        call output_char('(', ADVANCE='no', DONT_LOG=DONT_LOG)
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='no', DONT_LOG=DONT_LOG)
      endif
      call output_char( CHARS, &
        & ADVANCE, FROM_WHERE, DONT_LOG, &
        & LOG_CHARS, INSTEADOFBLANK )
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
    character(len=8) :: my_style
    integer :: n_chars
    character(len=3) :: MY_ADV
    !
    my_adv = Advance_is_yes_or_no(advance)
    my_style = TIMESTAMPSTYLE
    if ( present(style) ) my_style = lowercase(style)
    if ( my_style == 'post' ) then
      call output_integer( INT, PLACES, &
        & ADVANCE='no', FILL=FILL, FORMAT=FORMAT, BEFORE=BEFORE, AFTER=AFTER )
      if ( my_adv=='yes' ) then
        call output_char(' (', ADVANCE='no')
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='yes')
      endif
    else
      if ( ATLINESTART ) then
        call output_char('(', ADVANCE='no')
        call OUTPUT_DATE_AND_TIME(date=.false., advance='no')
        call output_char(')', ADVANCE='no')
      endif
      call output_integer( INT, PLACES, &
        & ADVANCE, FILL, FORMAT, BEFORE, AFTER )
    endif
  end subroutine timeStamp_integer

  ! Internal procedures
  
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

  ! ..............................................  not_used_here  .....
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module OUTPUT_M

! $Log$
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
