! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module STRING_TABLE

! String table module for compiler

  use HASH, only: HASH_LOOKUP => LOOKUP_AND_INSERT, HASH_FOUND => FOUND, &
                  INSERTED
  use IO_STUFF, only: GET_LUN
  use MACHINE, only: IO_ERROR
  use OUTPUT_M, only: OUTPUT
  implicit NONE
  private

  ! Public procedures
  public :: ADD_CHAR, ALLOCATE_CHAR_TABLE, ALLOCATE_HASH_TABLE
  public :: ALLOCATE_STRING_TABLE, CLEAR_STRING, COMPARE_STRINGS
  public :: DISPLAY_STRING, ENTER_STRING, FLOAT_VALUE, GET_CHAR, GET_STRING
  public :: HOW_MANY_STRINGS, LOOKUP_AND_INSERT, NEW_LINE, NUMERICAL_VALUE
  public :: OPEN_INPUT, STRING_LENGTH, STRING_TABLE_SIZE, UNGET_CHAR

  ! Public parameters
  character, public, parameter :: EOF = ACHAR(4)  ! ^D
  character, public, parameter :: EOL = ACHAR(10) ! ^J

  ! Public control variables
  logical, save, public :: DO_LISTING = .false.
  integer, save, public :: INUNIT = -1       ! Input unit, * if < 0

  ! Public information variables
  integer, save, public :: SOURCE_LINE = 0   ! Current line number
  integer, save, public :: SOURCE_COLUMN     ! Column number of current
                                             ! character

! =====     End of public declarations     ==================================

  ! Configuration parameter
  integer, parameter :: FILE_NAME_LEN = 127  ! Length of file name buffer
  ! in OPEN_INPUT

  ! State shared by GET_CHAR and OPEN_INPUT
  logical, save :: AT_EOF = .false.     ! Input is at EOF
  integer, save :: CUR_END = 0          ! Position of end of input
  integer, save :: CUR_POS = 1          ! Current position -- last one used

  ! Tables
  character, allocatable, save :: CHAR_TABLE (:)
  integer, save :: NCHARS     ! How full, not how big
  integer, allocatable, save :: HASH_TABLE (:,:)
  integer, allocatable, save :: STRINGS (:)  ! STRINGS(i) is the position
  ! in CHAR_TABLE of the last character of the i'th string
  integer, save :: NSTRING    ! How full, not how big

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
  ! =============================================     ADD_CHAR     =====
  subroutine ADD_CHAR ( CHARS, STAT )
  ! Add a character string to the last string in STRINGS.
    character(len=*), intent(in) :: CHARS
    integer, optional, intent(OUT) :: STAT
    integer :: I
    do i = 1, len(chars)
      if ( strings(nstring+1) >= ubound(char_table,1) ) then
        call double_char_table ( stat )
        if ( present(stat) ) then
          if ( stat /= 0 ) then; return; end if
        end if
      end if
      strings(nstring+1) = strings(nstring+1) + 1
      char_table(strings(nstring+1)) = chars(i:i)
    end do
    return
  end subroutine ADD_CHAR
  ! ==================================     ALLOCATE_CHAR_TABLE     =====
  subroutine ALLOCATE_CHAR_TABLE ( AMOUNT, STATUS )
    integer, intent(in) :: AMOUNT  ! How many characters to allocate
    integer, optional, intent(out) :: STATUS ! Status from ALLOCATE statement

    integer :: STAT

    if ( allocated( char_table ) ) then; deallocate ( char_table ); end if
    allocate ( char_table(amount), stat = stat )
    nchars = 0
    if ( present(status) ) then
      status = stat
      return
    end if
    if ( stat == 0 ) then; return; end if
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_CHAR_TABLE-E- Unable to allocate storage', &
      stat, '' )
    stop
  end subroutine ALLOCATE_CHAR_TABLE
  ! ==================================     ALLOCATE_HASH_TABLE     =====
  subroutine ALLOCATE_HASH_TABLE ( AMOUNT, STATUS )
    integer, intent(in) :: AMOUNT  ! How many hash table cells to allocate
                                   ! More might be allocated
    integer, optional, intent(out) :: STATUS ! Status from ALLOCATE statement

    integer MY_AMOUNT
    integer :: STAT

    if ( allocated( hash_table ) ) then; deallocate ( hash_table ); end if
    my_amount = amount
    ! Increase my_amount until it has no factors of 2, 3, or 5.
    if ( mod(my_amount, 2) == 0 ) then; my_amount = my_amount + 1; end if
    do
      if ( mod(my_amount, 3) == 0 ) then
        my_amount = my_amount + 2
        cycle
      end if
      if ( mod(my_amount, 5) == 0 ) then
        my_amount = my_amount + 2
        cycle
      end if
      exit
    end do
    allocate ( hash_table(2,my_amount+1), stat = stat )
    if ( present(status) ) then
      status = stat
      if ( stat == 0 ) then
        hash_table(1,:) = 0
        hash_table(2,1) = size(hash_table,2)
      end if
      return
    end if
    if ( stat == 0 ) then; return; end if
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_HASH_TABLE-E- Unable to allocate storage', &
      stat, '' )
    stop
  end subroutine ALLOCATE_HASH_TABLE
 ! =================================     ALLOCATE_STRING_TABLE     =====
  subroutine ALLOCATE_STRING_TABLE ( AMOUNT, STATUS )
    integer, intent(in) :: AMOUNT  ! How many strings to allocate
    integer, optional, intent(out) :: STATUS ! Status from ALLOCATE statement

    integer :: STAT

    if ( allocated( strings ) ) then; deallocate ( strings ); end if
    allocate ( strings(0:amount), stat = stat )
    strings(0:1) = 0
    nstring = 0
    if ( present(status) ) then
      status = stat
      return
    end if
    if ( stat == 0 ) then; return; end if
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_STRING_TABLE-E- Unable to allocate storage', &
      stat, '' )
    stop
  end subroutine ALLOCATE_STRING_TABLE
  ! =========================================     CLEAR_STRING     =====
  subroutine CLEAR_STRING
  ! Restart the last string in the string table.  Usually used after
  ! doing a lookup and finding the string is already there.
    strings(nstring+1) = strings(nstring)
    return
  end subroutine CLEAR_STRING
  ! ======================================     COMPARE_STRINGS     =====
  integer function COMPARE_STRINGS ( FIRST, SECOND, CASELESS )
  ! Returns 0 if the strings have the same text, <0 if the FIRST string
  ! should be sorted before the SECOND, and >0 if the  FIRST string
  ! should be sorted after the SECOND.  If CASELESS is present and
  ! TRUE, the case of letters is ignored.
    integer, intent(in) :: FIRST, SECOND     ! Indices in STRINGS
    logical, intent(in), optional :: CASELESS
    integer :: I, J      ! Subscripts, loop inductors
    logical :: NOCASE    ! .false. if CASELESS absent, else CASELESS
    nocase = .false.
    if ( present(caseless) ) nocase = caseless
    compare_strings = 0
    j = strings(second-1)
    if ( nocase ) then
      do i = strings(first-1)+1, strings(first)
        j = j + 1
        if ( j > strings(second) ) then ! FIRST string is longer
          compare_strings = 1
!         return
          go to 10
        end if
        compare_strings = iacap(char_table(i)) - iacap(char_table(j))
!       if ( compare_strings /= 0 ) return
        if ( compare_strings /= 0 ) go to 10
      end do
    else
      do i = strings(first-1)+1, strings(first)
        j = j + 1
        if ( j > strings(second) ) then ! FIRST string is longer
          compare_strings = 1
!         return
          go to 10
        end if
        compare_strings = iachar(char_table(i)) - iachar(char_table(j))
!       if ( compare_strings /= 0 ) return
        if ( compare_strings /= 0 ) go to 10
      end do
    end if
    if ( j < strings(second) ) compare_strings = -1 ! FIRST string is shorter
10  continue
!   print *, 'COMPARE_STRINGS = ', compare_strings
    return
  end function COMPARE_STRINGS
  ! =======================================     DISPLAY_STRING     =====
  subroutine DISPLAY_STRING ( STRING, ADVANCE, STRIP )
  ! Write the string indexed by STRING.
    integer, intent(in) :: STRING
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: STRIP

    integer :: offset
    logical :: myStrip
    character (len=1) :: firstChar, lastChar
    
    myStrip=.false.
    if (present(strip)) myStrip=strip

    offset=0
    if (myStrip) then
      firstChar=char_table(strings(string-1)+1)
      lastChar=char_table(strings(string))
      if ( (firstChar=='"' .and. lastChar=='"') .or. &
        &  (firstChar=="'" .and. lastChar=="'") ) offset=1
    end if

    call output ( char_table(strings(string-1)+1+offset: strings(string)-offset), &
                        advance=advance )
    return
  end subroutine DISPLAY_STRING
  ! =========================================     ENTER_STRING     =====
  integer function ENTER_STRING ()
  ! Commit the string being constructed, and return its index.
  ! The next call to ADD_CHAR will add the first character of a new
  ! string.
    nstring = nstring + 1
    if ( nstring > ubound(strings,1) ) then; call double_strings; end if
    strings(nstring+1) = strings(nstring)
    enter_string = nstring
    return
  end function ENTER_STRING
 ! ===========================================     FLOAT_VALUE     =====
  double precision function FLOAT_VALUE ( STRING )
  ! Return the FLOAT value of a string indexed by STRING
    integer, intent(in) :: STRING
    character(len=30) :: MY_CHAR
    call get_value ( &
      transfer(char_table(strings(string-1)+1:strings(string)), my_char), &
      strings(string) - strings(string-1) )
  contains
    subroutine GET_VALUE ( TEXT, N )
      integer, intent(in) :: N
      character(len=n), intent(in) :: TEXT
      read ( text, * ) float_value
    end subroutine GET_VALUE
  end function FLOAT_VALUE
  ! =============================================     GET_CHAR     =====
  subroutine GET_CHAR ( CHAR )
  ! Get a character from the input buffer (which is CHAR_TABLE).
  ! End-of-line is represented by EOL; end-of-file is represented by EOF.
  ! Once EOF has been read, all subsequent characters are EOF until input
  ! is re-opened.  The input is read ahead a line at a time into the
  ! CHAR_TABLE, as necessary, and listed with a line number if
  ! DO_LISTING is .true.
    character, intent(out) :: CHAR

    if ( at_eof ) then
      char = eof
      return
    end if
    if ( cur_pos >= cur_end ) then ! need to read a new line
      cur_end = strings(nstring+1)
      cur_pos = cur_end
      do
        cur_end = cur_end + 1
        if ( cur_end > ubound(char_table,1) ) then
          call double_char_table
        end if
        ! Non-advancing input is used so that trailing blanks in
        ! the record can be distinguished from padding.
        if ( inunit >= 0 ) then
          read ( inunit, '(a)', advance='no', eor=100, end=200 ) &
            char_table(cur_end)
        else
          read ( *, '(a)', advance='no', eor=100, end=200 ) &
            char_table(cur_end)
        end if
      end do
100   char_table(cur_end) = EOL
      go to 300
200   char_table(cur_end) = EOF
      at_eof = .true.
300   continue
      source_line = source_line + 1
      if ( .not. at_eof ) then
        if ( do_listing ) then
          call output ( source_line, 6 )
          call output ( '. ' )
          call output ( char_table(cur_pos+1:cur_end-1 ), advance='yes' )
        end if
      end if
      source_column = 0
    end if
    source_column = source_column + 1
    cur_pos = cur_pos + 1
    char = char_table(cur_pos)
    return
  end subroutine GET_CHAR
  ! ===========================================     GET_STRING     =====
  subroutine GET_STRING ( STRING, STRING_TEXT, CAP, STRIP )
  ! Put as much as will fit of the string indexed by STRING into STRING_TEXT.
  ! If CAP is present and .TRUE., capitalize STRING_TEXT.
    integer, intent(in) :: STRING
    character(len=*), intent(out) :: STRING_TEXT
    logical, intent(in), optional :: CAP
    logical, intent(in), optional :: STRIP
    integer :: I, J, offset
    logical :: MY_CAP, MY_STRIP
    my_cap = .false.
    my_strip = .false.
    offset = 0
    if ( present(cap) ) my_cap = cap
    if ( present(strip) ) my_strip = strip
    if (my_strip) then
      if ( ( (char_table(strings(string-1)+1) == '"') .and. &
        &    (char_table(strings(string)) == '"') ) .or.&
        &  ( (char_table(strings(string-1)+1) == "'") .and.&
        &    (char_table(strings(string)) == "'") ) ) &
        & offset=1
    endif
    call test_string ( string, 'GET_STRING' )
    string_text = ' '
    j = 0
    if ( my_cap ) then
      do i = strings(string-1)+1+offset, strings(string)-offset
        j = j + 1
        if ( j > len(string_text) ) exit
        string_text(j:j) = char(iacap(char_table(i)))
      end do
    else
      do i = strings(string-1)+1+offset, strings(string)-offset
        j = j + 1
        if ( j > len(string_text) ) exit
        string_text(j:j) = char_table(i)
      end do
!     string_text = transfer(char_table(strings(string-1)+1:strings(string)), &
!                            string_text(:strings(string)-strings(string-1))
    end if
    return
  end subroutine GET_STRING
  ! =====================================     HOW_MANY_STRINGS     =====
  integer function HOW_MANY_STRINGS ()
  ! Returns the number of strings in the string table
    how_many_strings = nstring
    return
  end function HOW_MANY_STRINGS
  ! ====================================     LOOKUP_AND_INSERT     =====
  subroutine LOOKUP_AND_INSERT ( STRING, FOUND, CASELESS )
  ! Look for the string built up by Add_Char.  If it is found return the
  ! position at which it was found in STRING, and FOUND = .true.  If it is
  ! not found, add it, return the position at which it was added in
  ! STRING, and FOUND = .false.  In any case, the next call to Add_Char
  ! will add the first character of a new string. If CASELESS is present
  ! and .true., compare caseless.
    integer, intent(out) :: STRING
    logical, intent(out) :: FOUND
    logical, optional, intent(in) :: CASELESS

    integer :: HASH_KEY  ! Integer derived from characters of STRING
    integer :: I         ! Subscript, loop inductor
    integer :: LOC       ! Where HASH_KEY was found in HASH_TABLE
    logical :: NOCASE    ! .false. or CASELESS
    integer :: STATUS    ! Result, see HASH_LOOKUP

    nocase = .false.
    if ( present(caseless) ) nocase = caseless

    ! Construct a hash_key by adding up the numeric representations
    ! of the characters in NSTRING+1, ignoring overflows
    hash_key = 0
    if ( nocase ) then
      do i = strings(nstring)+1, strings(nstring+1)
        hash_key = hash_key + iacap ( char_table(i) )
      end do
    else
      do i = strings(nstring)+1, strings(nstring+1)
        hash_key = hash_key + iachar ( char_table(i) )
      end do
    end if

    ! Look for hash_key
    loc = 0
    do
      call hash_lookup ( hash_key, hash_table, .true., loc, status )
      if ( status == inserted ) then
        found = .false.
        nstring = nstring + 1
        if ( nstring >= ubound(strings,1) ) then; call double_strings; end if
        hash_table(2,loc) = nstring
        string = nstring
        strings(nstring+1) = strings(nstring)
        return
      end if
      if ( status /= hash_found ) then
        write ( *, * ) 'STRING_TABLE%LOOKUP_AND_INSERT-E- ', &
                       'Either the hash table is full'
        write ( *, * ) 'or the program is using the hash software ', &
                       'incorrectly.'
        write ( *, * ) 'Unfortunately, the program doesn''t know how ', &
                       'to increase the hash table size.'
        stop
      end if
      ! The hash key matches; check whether the string does
!write (*, *) 'Compare', hash_table(2,loc), ': ', &
!char_table(strings(hash_table(2,loc)-1)+1:strings(hash_table(2,loc))), &
!' to ', nstring+1, ': ', char_table(strings(nstring)+1:strings(nstring+1))

!  call output ( 'Compare ', advance='no' )
!  call output ( hash_table(2,loc) ); call output ( ': ', advance='no')
!  call display_string ( hash_table(2,loc) )
!  call output ( ' to ' )
!  call output ( nstring+1 ); call output ( ': ')
!  call display_string ( nstring+1, advance='no' )
!  if ( nocase ) call output ( ' caseless')
!  call output ( '', advance='yes' )
      string = hash_table(2,loc)
      found = compare_strings ( string, nstring+1, caseless ) == 0
      if ( found ) then
        call clear_string
        return
      end if
    end do
  end subroutine LOOKUP_AND_INSERT
 ! ==============================================     NEW_LINE     =====
  subroutine NEW_LINE
  ! Skip the rest of the current line, set up so that get_char will
  ! read a new one.
    cur_pos = cur_end + 1
    return
  end subroutine NEW_LINE
 ! =======================================     NUMERICAL_VALUE     =====
  integer function NUMERICAL_VALUE ( STRING )
  ! Return the numerical value of a string indexed by STRING
    integer, intent(in) :: STRING
    character(len=30) :: MY_CHAR
    call get_value ( &
      transfer(char_table(strings(string-1)+1:strings(string)), my_char), &
      strings(string) - strings(string-1) )
  contains
    subroutine GET_VALUE ( TEXT, N )
      integer, intent(in) :: N
      character(len=n), intent(in) :: TEXT
      read ( text, * ) numerical_value
    end subroutine GET_VALUE
  end function NUMERICAL_VALUE
  ! ===========================================     OPEN_INPUT     =====
  subroutine OPEN_INPUT ( FILE_NAME, STAT )
  ! Open the file given by FILE_NAME for input.  If it can't be opened and
  ! STAT is present, return the status.  If it can't be opened and STAT is
  ! absent, ask the user for it.
    character(len=*), intent(in) :: FILE_NAME
    integer, optional, intent(out) :: STAT

    integer :: IOSTAT
    character(file_name_len) :: MY_FILE

    call get_lun ( inunit )   ! get an unused I/O unit number
    my_file = file_name
    do
      open ( inunit, file=my_file, status='OLD', access='SEQUENTIAL', &
             form='FORMATTED', iostat=iostat )
      if ( iostat == 0 ) then
        at_eof = .false.
        source_line = 0
        return
      end if
      if ( present(stat) ) then
        stat = iostat
        return
      end if
      call io_error ( 'STRING_TABLE%OPEN_INPUT-E- Unable to open input file', &
                      iostat, my_file )
      write ( *, * ) 'Enter input file name: '
      read ( *, '(a)', end=999 ) my_file
    end do
999 stop
  end subroutine OPEN_INPUT
 ! =========================================     STRING_LENGTH     =====
  integer function STRING_LENGTH ( STRING )
  ! Return the length of the string indexed by STRING
    integer, intent(in) :: STRING
    call test_string ( string, 'STRING_LENGTH' )
    string_length = strings(string) - strings(string-1)
    return
  end function STRING_LENGTH
  ! ====================================     STRING_TABLE_SIZE     =====
  integer function STRING_TABLE_SIZE ()
  ! Return the allocate size of the string table.
    string_table_size = ubound(strings, 1)
    return
  end function STRING_TABLE_SIZE
 ! ============================================     UNGET_CHAR     =====
  subroutine UNGET_CHAR
  ! Back up the input, but not before the beginning of the current line
    cur_pos = max( cur_pos-1, strings(nstring+1) )
    return
  end subroutine UNGET_CHAR

! =====     Private procedures    ======================================
  ! ------------------------------------     DOUBLE_CHAR_TABLE     -----
  subroutine DOUBLE_CHAR_TABLE ( STAT )
  ! Double the space for the character table
    integer, optional, intent(out) :: STAT
    integer :: MY_STAT
    character, allocatable :: OLD_CHAR(:)
    allocate ( old_char(ubound(char_table,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'stat_TABLE%DOUBLE_CHARS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    old_char = char_table
    deallocate ( char_table )
    allocate( char_table(2*ubound(old_char,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'stat_TABLE%DOUBLE_CHARS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    char_table(:ubound(old_char,1)) = old_char
    deallocate ( old_char )
    return
  end subroutine DOUBLE_CHAR_TABLE
  ! ---------------------------------------     DOUBLE_STRINGS     -----
  subroutine DOUBLE_STRINGS ( STAT )
  ! Double the space for the string table
    integer, optional, intent(out) :: STAT
    integer :: MY_STAT
    integer, allocatable :: OLD_STRING(:)
    allocate ( old_string(0:ubound(strings,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'STRING_TABLE%DOUBLE_STRINGS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    old_string = strings
    deallocate ( strings )
    allocate( strings(0:2*ubound(old_string,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'STRING_TABLE%DOUBLE_STRINGS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    strings(0:ubound(old_string,1)) = old_string
    deallocate ( old_string )
    return
  end subroutine DOUBLE_STRINGS
 ! -------------------------------------------------     IACAP     -----
  integer function IACAP ( CHAR )
  ! Returns the ASCII index of a character, except if it is a lower
  ! case letter, it returns the upper case index.
    character, intent(in) :: CHAR
    integer, parameter :: SHIFT = iachar('A') - iachar('a')
    iacap = iachar ( char )
    if ( char >= 'a' .and. char <= 'z' ) then
      iacap = iacap + shift
    end if
  end function IACAP
  ! ------------------------------------------     TEST_STRING     -----
  subroutine TEST_STRING ( STRING, ROUTINE )
  ! Test whether STRING is within bounds.  If not, use ROUTINE to emit
  ! an error message.
    integer, intent(in) :: STRING
    character(len=*), intent(in) :: ROUTINE
    if ( string < 1 .or. string > nstring ) then
      write ( *, * ) 'STRING_TABLE%', routine, '-E- String index ', &
        string, ' not in 1 .. ', nstring
      stop
    end if
    return
  end subroutine TEST_STRING
end module STRING_TABLE

! $Log$
! Revision 2.4  2001/04/05 00:54:59  vsnyder
! Correct 'increase table sizes automatically' code
!
! Revision 2.3  2001/03/03 00:07:24  livesey
! Added strip argument to get_string
!
! Revision 2.2  2001/03/02 01:33:40  livesey
! Added strip argument to display_string
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.2  2000/08/08 19:44:51  vsnyder
! Removed an inaccessible RETURN statement at line 406.
!
