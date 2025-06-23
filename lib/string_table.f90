! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module STRING_TABLE

! String table module for compiler

  use hash, only: hash_lookup => lookup_and_insert, hash_found => found, &
                  inserted, hash_full => full, hash_not_key => not_key, &
                  hash_bad_loc => bad_loc, hash_empty => empty
  use io_stuff, only: get_lun
  use machine, only: crash_burn, io_error
  use output_m, only: blanks, NewLine, output
  use printIt_m, only: MLSMSG_Deallocate, MLSMSG_Allocate, MLSMSG_Error, &
                       printItOut
  implicit NONE
  private

  ! Public procedures
  public :: add_char, add_include_directory, addinunit, allocate_char_table
  public :: allocate_hash_table, allocate_string_table, clear_string
  public :: compare_strings, create_string, destroy_char_table
  public :: destroy_hash_table, destroy_string_table, display_string
  public :: display_string_list, dump_char_table, dump_inunit_stack, enter_string
!   ifort v17 gets a seg fault if find_file is public.
!   public :: find_file
  public :: float_value, get_char, get_string, Get_String_Char
  public :: how_many_strings, include_stack_top
  public :: index, init_string_table, index_in_string, isStringInTable, len, lookup
  public :: lookup_and_insert, new_line, numerical_value, open_include, open_input
  public :: string_length, string_table_size, unget_char

  ! Public procedure pointers.  Use this to find whether a file exists, get
  ! its full name, and inquire whether it is already open, using something
  ! other than INQUIRE and OPEN statements directly.
  ! This initialization doesn't work with nagfor 1052 or ifort 15.0.2
  ! procedure(Find_File), pointer :: Find_A_File => Find_File
  procedure(Find_File), pointer, public :: Find_A_File => null()

  ! Public variables and named constants
  public :: Do_Listing, EOF, EOL, Includes, Source_Line, Source_Column

  interface DISPLAY_STRING
    module procedure DISPLAY_STRING, DISPLAY_STRING_LIST
  end interface

  interface INDEX
    module procedure INDEX_IN_STRING
  end interface

  interface LEN
    module procedure STRING_LENGTH
  end interface

  ! Public parameters
  character, parameter :: EOF = ACHAR(4)  ! ^D
  character, parameter :: EOL = ACHAR(10) ! ^J

  ! Public control variables
  logical, save :: DO_LISTING = .false.

  ! Public information variables

  integer, save, protected, allocatable :: Includes(:) ! String table indices
                                              ! of directories to search for
                                              ! files; see Find_File
  integer, save, protected :: SOURCE_LINE = 0 ! Current line number
  integer, save, protected :: SOURCE_COLUMN   ! Column number of current
                                              ! character

! =====     End of public declarations     ==================================

  ! Configuration parameter
  integer, parameter :: FILE_NAME_LEN = 127  ! Length of file name buffer
  ! in OPEN_INPUT

  ! State shared by GET_CHAR and OPEN_INPUT
  logical, save :: AT_EOF = .false.     ! Input is at EOF
  integer, save :: CUR_END = 0          ! Position of end of input
  integer, save :: CUR_POS = 1          ! Current position -- last one used
  integer, save :: inunit_counter = 0
  integer, dimension(:), pointer, save :: inunit_list => NULL() ! Input unit, * if null

  type :: Inunit_Stack_t
    integer :: File = 0      ! string table index of complete file name
    integer :: Unit = 0      ! I/O unit
    integer :: Line = 0      ! in the file
    integer :: Column = 0    ! in the line
  end type Inunit_Stack_t

  type(inunit_stack_t), allocatable :: Inunit_Stack(:)
  integer, parameter :: Init_Stack_Size = 10

  integer :: Inunit_Stack_Ptr = 0

  ! Tables
  character, allocatable, save :: CHAR_TABLE (:)
  integer, allocatable, save   :: HASH_TABLE (:,:)
  integer, allocatable, save   :: STRINGS (:)  ! STRINGS(i) is the position
  ! in CHAR_TABLE of the last character of the i'th string
  integer, save                :: NSTRING = 0   ! How full, not how big (??)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
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
          if ( stat /= 0 ) return
        end if
      end if
      strings(nstring+1) = strings(nstring+1) + 1
      char_table(strings(nstring+1)) = chars(i:i)
    end do
  end subroutine ADD_CHAR
  ! ================================     ADD_INCLUDE_DIRECTORY     =====
  subroutine ADD_INCLUDE_DIRECTORY ( Directory, Stat )
    character(len=*), intent(in) :: Directory
    integer, intent(out), optional :: Stat
    integer :: Dir_String ! String index of directory
    integer :: MyStat
    integer, allocatable :: Temp_List(:)
    dir_string = create_string ( trim(directory) )
    if ( .not. allocated(includes) ) then
      allocate ( temp_list(1), stat=myStat )
    else
      allocate ( temp_list(size(includes)+1), stat=myStat )
      temp_list(:size(includes)) = includes
    end if
    if ( present(stat) ) stat = myStat
    if ( myStat /= 0 ) then
      if ( present(stat) ) return
      call io_error &
      ( 'STRING_TABLE%ALLOCATE_CHAR_TABLE-E- Unable to allocate include directory list', &
      & myStat )
    end if
    temp_list(size(temp_list)) = dir_string
    call move_alloc ( temp_list, includes )
  end subroutine ADD_INCLUDE_DIRECTORY
  ! ==================================     ALLOCATE_CHAR_TABLE     =====
  subroutine ALLOCATE_CHAR_TABLE ( AMOUNT, STATUS )
    integer, intent(in) :: AMOUNT  ! How many characters to allocate
    integer, optional, intent(out) :: STATUS ! Status from ALLOCATE statement

    integer :: STAT

    if ( allocated( char_table ) ) then; deallocate ( char_table ); end if
    allocate ( char_table(amount), stat = stat )
    if ( present(status) ) then
      status = stat
      return
    end if
    if ( stat == 0 ) return
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_CHAR_TABLE-E- Unable to allocate storage', &
      stat )
    call crash_burn
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
    if ( stat == 0 ) return
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_HASH_TABLE-E- Unable to allocate storage', &
      stat )
    call crash_burn
  end subroutine ALLOCATE_HASH_TABLE
 ! =================================     ALLOCATE_STRING_TABLE     =====
  subroutine ALLOCATE_STRING_TABLE ( AMOUNT, STATUS )
    integer, intent(in) :: AMOUNT  ! How many strings to allocate
    integer, optional, intent(out) :: STATUS ! Status from ALLOCATE statement

    integer :: STAT

    if ( allocated( strings ) ) deallocate ( strings )
    allocate ( strings(0:amount), stat = stat )
    strings(0:1) = 0
    nstring = 0
    if ( present(status) ) then
      status = stat
      return
    end if
    if ( stat == 0 ) return
    call io_error &
      ( 'STRING_TABLE%ALLOCATE_STRING_TABLE-E- Unable to allocate storage', &
      stat )
    call crash_burn
  end subroutine ALLOCATE_STRING_TABLE
  ! =========================================     CLEAR_STRING     =====
  subroutine CLEAR_STRING
  ! Restart the last string in the string table.  Usually used after
  ! doing a lookup and finding the string is already there.
    strings(nstring+1) = strings(nstring)
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
  end function COMPARE_STRINGS
  ! ========================================     CREATE_STRING     =====
  integer function CREATE_STRING ( TEXT, CASELESS, DEBUG ) result ( STRING )
  ! Put the characters of TEXT into the end of the string table using
  ! Add_Char.  Look up the string with LOOKUP_AND_INSERT using STRING,
  ! CASELESS and DEBUG arguments.  Don't bother returning whether the
  ! string was found or inserted.
    character(len=*), intent(in) :: TEXT
    logical, optional, intent(in) :: CASELESS
    integer, optional, intent(in) :: DEBUG

    logical :: FOUND

    call add_char ( text )
    call lookup_and_insert ( string, found, caseless, debug )
  end function CREATE_STRING
  ! ===================================     DESTROY_CHAR_TABLE     =====
  subroutine DESTROY_CHAR_TABLE ( STATUS )
    integer, intent(out), optional :: STATUS ! From deallocate
    if ( allocated(char_table) ) then
      if ( present(status) ) then
        deallocate ( char_table, stat=status )
      else
        deallocate ( char_table )
      end if
    end if
  end subroutine DESTROY_CHAR_TABLE
  ! ===================================     DESTROY_HASH_TABLE     =====
  subroutine DESTROY_HASH_TABLE ( STATUS )
    integer, intent(out), optional :: STATUS ! From deallocate
    if ( allocated(hash_table) ) then
      if ( present(status) ) then
        deallocate ( hash_table, stat=status )
      else
        deallocate ( hash_table )
      end if
    end if
  end subroutine DESTROY_HASH_TABLE
  ! =================================     DESTROY_STRING_TABLE     =====
  subroutine DESTROY_STRING_TABLE ( STATUS )
    integer, intent(out), optional :: STATUS ! From deallocate
    if ( allocated(strings) ) then
      if ( present(status) ) then
        deallocate ( strings, stat=status )
      else
        deallocate ( strings )
      end if
    end if
    nstring = 0
  end subroutine DESTROY_STRING_TABLE
  ! =======================================     DISPLAY_STRING     =====
  subroutine DISPLAY_STRING ( String, Advance, Strip, Ierr, Before, Pos, Width )
  ! Write BEFORE (if present), followed by the string indexed by STRING.
  ! If IERR is present return 0 if all went well, 1 otherwise.
  ! If Pos is present, add len(Before) + string_length(String) to it, using zero
  ! for len(Before) if Before is absent.  If Pos and Width are present, and Pos
  ! was > 0, and Pos > Width, advance regardless of the value of Advance.
    integer, intent(in) :: String
    character(len=*), intent(in), optional :: Advance
    logical, intent(in), optional :: Strip
    integer, intent(out), optional :: Ierr
    character(len=*), intent(in), optional :: Before
    integer, intent(inout), optional :: Pos ! Current position in line
    integer, intent(in), optional :: Width  ! Maximum line width

    logical :: Fail     ! String is not in the string table
    character (len=1) :: firstChar, lastChar
    character(len=*), parameter :: Msg = '(not found in string table)'
    logical :: myStrip
    integer :: offset
    integer :: W        ! width of string or msg

    fail = .false.
    call test_string ( string, 'Display_String', ierr )
    if ( present(ierr) ) fail = ierr /= 0

    if ( present(pos) ) then
      w = len(msg)
      if ( .not. fail ) w = string_length(string)
      offset = pos ! Temp to avoid new line if the line is empty
      if ( present(before) ) pos = pos + len(before)
      pos = pos + w
      if ( present(width) .and. offset > 0 ) then
        if ( pos > width ) call newLine
      end if
    end if
    if ( present(before) ) call output ( before )
    if ( fail ) then
      call output ( msg, advance=advance )
      return
    end if

    offset=0

    myStrip=.false.
    if (present(strip)) myStrip=strip
    if ( myStrip ) then
      firstChar=char_table(strings(string-1)+1)
      lastChar=char_table(strings(string))
      if ( (firstChar=='"' .and. lastChar=='"') .or. &
        &  (firstChar=="'" .and. lastChar=="'") ) offset=1
    end if

    call output ( char_table(strings(string-1)+1+offset: strings(string)-offset), &
                        advance=advance )
  end subroutine DISPLAY_STRING
  ! ==================================     DISPLAY_STRING_LIST     =====
  subroutine DISPLAY_STRING_LIST ( STRING, ADVANCE, STRIP, IERR, BEFORE )
  ! Write BEFORE (if present), followed every string indexed by STRING, each
  ! preceded by a blank. ADVANCE applies only to the last one.  See
  ! DISPLAY_STRING for other arguments.
    integer, intent(in) :: STRING(:)
    character(len=*), intent(in), optional :: ADVANCE
    logical, intent(in), optional :: STRIP
    integer, intent(out), optional :: IERR
    character(len=*), intent(in), optional :: BEFORE
    integer :: I, JERR
    jerr = 0
    if ( present(before) ) call output ( before )
    do i = 1, size(string)-1
      call display_string ( string(i), strip=strip, ierr=ierr, before=' ' )
      if ( present(ierr) ) jerr = max(jerr,ierr)
    end do
    call display_string ( string(size(string)), advance, strip, ierr, before=' ')
    if ( present(ierr) )  ierr = max(jerr,ierr)
  end subroutine DISPLAY_STRING_LIST
  ! ======================================     Dump_Char_Table     =====
  subroutine Dump_Char_Table ( Start, End )
    ! This is for debugging, assuming you have somehow gotten hold of
    ! how much of the character table you want to dump, probably by
    ! inserting some dumping stuff elsewhere in this module
    integer, intent(in) :: Start, End
    integer :: I, J, K, MyStart, MyEnd
    myStart = start
    if ( start < 1 ) then
      call output ( start, before='Requested start ' )
      call output ( 1, before=' changed to ', advance='yes' )
      myStart = 1
    end if
    myEnd = end
    if ( myEnd > cur_end ) then
      call output ( end, before='Requested end ' )
      call output ( cur_end, before=' changed to ', advance='yes' )
      myEnd = cur_end
    end if
    do i = myStart, myEnd, 100
      j = min(i+99,myEnd)
      call output ( i, before='Char_Table(' )
      call output ( j, before=':' )
      call output ( ') = "' )
      do k = i, j
        if ( char_table(k) == EOL ) then
          call output ( 'EOL' )
        else if ( char_table(k) == EOF ) then
          call output ( 'EOF' )
        else
          call output ( char_table(k) )
        end if
      end do
      call output ( '"', advance='yes' )
    end do
  end subroutine Dump_Char_Table
  ! ====================================     DUMP_INUNIT_STACK     =====
  subroutine DUMP_INUNIT_STACK
    integer :: I
    if ( .not. allocated(inunit_stack) .or. inunit_stack_ptr <= 0 ) then
      call output ( "There is no inunit stack to dump", advance="yes" )
      return
    end if
    call output ( "Inunit Stack, top down", advance="yes" )
    do i = inunit_stack_ptr, 1, -1
      call output ( i, format='(i3)' )
      call output ( inunit_stack(i)%unit, before=': unit ' )
      call output ( inunit_stack(i)%line, before=', line ' )
      call output ( inunit_stack(i)%column, before=', column ' )
      call output ( inunit_stack(i)%file, before=', file ' )
      call display_string ( inunit_stack(i)%file, before=': ', advance='yes' )
    end do
  end subroutine DUMP_INUNIT_STACK
  ! ====================================     DUMP_STRING_TABLE     =====
  subroutine DUMP_STRING_TABLE
  ! Dump the entire string table
    character(len=7) :: FMT ! '(I..: )'
    integer :: I
    i = log10(real(max(nstring,1))) + 1
    write ( fmt, "('(I',i2.2,': ')" ) i
    do i = 1, nstring
      call output ( i, format=fmt, advance='no' )
      call display_string ( i, advance='yes' )
    end do
  end subroutine DUMP_STRING_TABLE
  ! =========================================     ENTER_STRING     =====
  integer function ENTER_STRING ()
  ! Commit the string being constructed, and return its index.
  ! The next call to ADD_CHAR will add the first character of a new
  ! string.
    nstring = nstring + 1
    if ( nstring > ubound(strings,1) ) then; call double_strings; end if
    strings(nstring+1) = strings(nstring)
    enter_string = nstring
  end function ENTER_STRING
  ! ============================================     Find_File     =====
  subroutine Find_File ( Directories, File_Name, Exist, Full_Text, Opened, &
                       & Unit, Stat )
  ! Look for File_Name in directories.  If it is found, set Exist=true
  ! and put its full text in Full_Text, otherwise set Exist=false.
    integer, intent(in) :: Directories(:)   ! String indices
    integer, intent(in) :: File_Name        ! String index
    logical, intent(out) :: Exist
    character(len=*), intent(out) :: Full_Text
    logical, intent(out) :: Opened          ! "The file is already open"
    integer, intent(out), optional  :: Unit ! An unused unit number
    integer, intent(out), optional :: Stat
    logical :: Change
    integer :: I, L
    integer :: MyUnit

    ! No point in working hard on this if there aren't any available units
    call get_lun ( myUnit )   ! get an unused I/O unit number
    if ( myUnit < 0 ) then    ! no LUNs available
      if ( present(stat) ) then
        stat = myUnit
        return
      end if
      call output ( 'STRING_TABLE%OPEN_INPUT-E- Unable to get LUN', advance='yes' )
      call crash_burn
    end if

    exist = .false.
    call get_string ( file_name, full_text, strip=.true. )
    inquire ( file=full_text, exist=exist )
    if ( .not. exist ) then
      do i = 1, size(directories)
        l = 1
        if ( directories(i) > 0 ) then
          call get_string ( directories(i), full_text, strip=.true. )
          l = string_length(directories(i)) + 1
          if ( full_text(l-1:l-1) /= '/' ) then
            full_text(l:l) = '/'
            l = l + 1
          end if
        end if
        call get_string ( file_name, full_text(l:), strip=.true. )
        inquire ( file=full_text, exist=exist )
        if ( exist ) exit
      end do
    end if
    ! Remove multiple slashes from path name and replace instances of /./
    ! with /; inquire opened seems not to catch circular includes in these cases
    change = .true.
    do while ( change )
      change = .false.
      i = index(full_text, '//')
      if ( i /= 0 ) then
        change = .true.
        full_text(i:) = full_text(i+1:)
      end if
      i = index(full_text, '/./')
      if ( i /= 0 ) then
        change = .true.
        full_text(i:) = full_text(i+2:)
      end if
    end do
    ! Remove prefixes of "./"; inquire opened seems not
    ! to catch circular includes in this case
    do while ( full_text(1:2) == './' )
      full_text = full_text(3:)
    end do
    inquire ( file=full_text, exist=exist )
    opened = .false.
    if ( exist ) inquire ( file=full_text, opened=opened )
    if ( present(unit) ) then
      unit = myUnit
    else
      call AddInUnit(myUnit)
    end if
  end subroutine Find_File
 ! ===========================================     FLOAT_VALUE     =====
  double precision function FLOAT_VALUE ( STRING )
  ! Return the FLOAT value of a string indexed by STRING
    integer, intent(in) :: STRING
    character(len=30) :: MY_CHAR
    call test_string ( string, 'Float_Value' )
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

    integer :: IOSTAT

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
        if ( inunit_stack_ptr > 0 ) then
          read ( inunit_stack(inunit_stack_ptr)%unit, '(a)', advance='no', &
            & eor=100, end=200, err=400, iostat=iostat ) char_table(cur_end)
        else if ( associated(inunit_list) ) then
          read ( inunit_list(inunit_counter), '(a)', advance='no', eor=100, &
            & end=200, err=400, iostat=iostat ) char_table(cur_end)
        else
          read ( *, '(a)', advance='no', eor=100, end=200, err=400, &
            iostat=iostat ) char_table(cur_end)
        end if
      end do
100   char_table(cur_end) = EOL
      go to 300
200   if ( inunit_stack_ptr > 0 ) then
        close ( inunit_stack(inunit_stack_ptr)%unit ) ! Ignore status
        source_line = inunit_stack(inunit_stack_ptr)%line
        source_column = inunit_stack(inunit_stack_ptr)%column
        inunit_stack_ptr = inunit_stack_ptr-1
      else
        call SetNextInUnit
      end if
      if (inunit_counter == 0) then ! has no more file to read from
        char_table(cur_end) = EOF
        at_eof = .true.
      else
        char_table(cur_end) = EOL
      end if
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
400 call io_error ( 'While reading input in String_Table%Get_Char', iostat )
    call crash_burn
  end subroutine GET_CHAR
  ! ===========================================     GET_STRING     =====
  subroutine GET_STRING ( STRING, STRING_TEXT, CAP, STRIP, NOERROR, IERR, &
    & START, END )
  ! Put as much as will fit of the string indexed by STRING into STRING_TEXT.
  ! If CAP is present and .TRUE., capitalize STRING_TEXT.
  ! If STRIP is present and .TRUE., remove quotes if any.
  ! If NOERROR is present and TRUE, return safely no matter what
  ! If IERR is present and error occurs, set it to 1 && return safely
  ! If START is not present start at offset where offset is 2 if STRIP is
  !   present and true and the string is quoted else 1, else start at
  !   max(START,1) + offset.
  ! If END is not present, end at the end of the string, or end-1 if STRIP
  !   is present and true and the string is quoted, else end at
  !   min(length,end)-offset.
    integer, intent(in) :: STRING
    character(len=*), intent(out) :: STRING_TEXT
    logical, intent(in), optional :: CAP
    logical, intent(in), optional :: STRIP
    logical, intent(in), optional :: NOERROR
    integer, intent(out), optional :: IERR
    integer, intent(in), optional :: START
    integer, intent(in), optional :: END
    integer :: I, J, MY_END, MY_IERR, MY_START, Offset
    logical :: MY_CAP, MY_NOERROR, MY_STRIP

    my_noerror = .false.
    if ( present(noError) ) my_noerror = noerror

    if ( my_noerror ) then
      call test_string ( string, 'GET_STRING', my_ierr )
      if (present(ierr)) ierr = my_ierr
    else
      ! Won't return if there's an error and IERR is not present
      call test_string ( string, 'GET_STRING', ierr )
      my_ierr = 0
      if ( present(ierr) ) my_ierr = ierr
    end if

    if ( my_ierr /= 0 ) then ! Don't leave string_text undefined
      write ( string_text, '("String index ",i0, " not in 1:")' ) string
      i = len_trim(string_text)
      write ( string_text(i+1:), '(i0)' ) nstring
      return
    end if

    my_cap = .false.
    my_strip = .false.
    if ( present(cap) ) my_cap = cap
    if ( present(strip) ) my_strip = strip

    offset = 0
    if (my_strip) then
      if ( ( (char_table(strings(string-1)+1) == '"') .and. &
        &    (char_table(strings(string)) == '"') ) .or.&
        &  ( (char_table(strings(string-1)+1) == "'") .and.&
        &    (char_table(strings(string)) == "'") ) ) &
        & offset=1
    end if
    my_start = 1
    if ( present(start) ) my_start = max(start, 1)
    my_start = strings(string-1) + my_start + offset
    my_end = strings(string) - offset
    if ( present(end) ) my_end = min(strings(string-1) + end, my_end)

    j = 0
    if ( my_cap ) then
      do i = my_start, my_end
        j = j + 1
        if ( j > len(string_text) ) exit
        string_text(j:j) = char(iacap(char_table(i)))
      end do
    else
      do i = my_start, my_end
        j = j + 1
        if ( j > len(string_text) ) exit
        string_text(j:j) = char_table(i)
      end do
!     string_text = transfer(char_table(strings(string-1)+1:strings(string)), &
!                            string_text(:strings(string)-strings(string-1))
    end if
!    if ( present(ierr) ) ierr=0 ! Already zeroed by call to test_string
    string_text(j+1:) = '' ! Fill the rest with blanks
  end subroutine GET_STRING
  ! ======================================     GET_STRING_CHAR     =====

  pure character function Get_String_Char ( S, I )

    ! Return the I'th character of the string at S.

    ! If I < 1 or I > String_Length(s) return blank

    integer, intent(in) :: S

    integer, intent(in) :: I

    integer :: J

    j = strings(s-1) + i

    if ( i < 1 .or. j > strings(s) ) then

      get_string_char = ' '

    else

      get_string_char = char(j)

    end if

  end function Get_String_Char
  ! =====================================     HOW_MANY_STRINGS     =====
  integer function HOW_MANY_STRINGS ()
  ! Returns the number of strings in the string table
    how_many_strings = nstring
  end function HOW_MANY_STRINGS
  ! ====================================     INCLUDE_STACK_TOP     =====
  subroutine INCLUDE_STACK_TOP ( File )
  ! Get the directory and file name string indices from the top of the
  ! include stack, if there is anything in it, else zeros.
    integer, intent(out) :: File
    file = 0
    if ( .not. allocated(inunit_stack) ) return
    if ( inunit_stack_ptr == 0 ) return
    file = inunit_stack(inunit_stack_ptr)%file
  end subroutine INCLUDE_STACK_TOP
  ! ======================================     INDEX_IN_STRING     =====
  integer function INDEX_IN_STRING ( STRING, SUBSTRING, CASELESS, STRIP )
  ! Works like intrinsic INDEX, but with integer string indices instead of
  ! characters, without the BACK argument, and with the CASELESS argument.
  ! If STRIP is present and true do not consider quotes or apostrophes at
  ! the ends of the strings.
  ! Use the brute-force method instead of a fancy substring method such as
  ! Knuth-Morris-Pratt or Rabin-Karp.  The result value can be used as the
  ! START argument for GET_STRING, provided GET_STRING is invoked with
  ! STRIP having the same value, or STRING and SUBSTRING are not quoted.

    integer, intent(in) :: String
    integer, intent(in) :: Substring
    logical, intent(in), optional :: Caseless ! Default false
    logical, intent(in), optional :: Strip    ! Default false
    character :: C1, C2
    integer, parameter :: Diff = iachar('A') - iachar('a')
    integer :: I, I1, I2, J, J1, J2
    logical :: MyCaseless, MyStrip
    character(len=*), parameter :: Quotes = '"'//"'"

    myCaseless = .false.
    if ( present(caseless) ) myCaseless = caseless
    myStrip = .false.
    if ( present(strip) ) myStrip = strip

    index_in_string = 0
    j1 = 1
    j2 = len(substring)
    if ( myStrip ) then
      j1 = 1  + scan(char_table(strings(substring-1)+1),quotes)
      j2 = j2 - scan(char_table(strings(substring)),    quotes)
    end if
    i1 = 1
    i2 = len(string) - j2 + j1
    if ( myStrip ) then
      i1 = 1  + scan(char_table(strings(string-1)+1),   quotes)
      i2 = i2 - scan(char_table(strings(string)),       quotes)
    end if
    if ( .not. myCaseless ) then ! case sensitive compare
 o:   do i = i1, i2
        do j = j1, j2
          if ( char_table(strings(string-1)+i+j-1) /= &
            &  char_table(strings(substring-1)+j) ) cycle o
        end do
        index_in_string = i
        return
      end do o
    else ! compare in upper case
 l:   do i = i1, i2
        do j = j1, j2
          c1 = char_table(strings(string-1)+i+j-1)
          if ( c1 >= 'a' .and. c1 <= 'z' ) c1 = achar(iachar(c1)+diff)
          c2 = char_table(strings(substring-1)+j)
          if ( c2 >= 'a' .and. c2 <= 'z' ) c2 = achar(iachar(c2)+diff)
          if ( c1 /= c2 ) cycle l
        end do
        index_in_string = i
        return
      end do l
    end if
  end function INDEX_IN_STRING
  ! ====================================     LOOKUP     =====
  subroutine LOOKUP ( STRING, FOUND, CASELESS, DEBUG )
  ! Look for the string built up by Add_Char.  If it is found return the
  ! position at which it was found in STRING, and FOUND = .true.  If it is
  ! not found, set FOUND = .false.  In any case, the next call to Add_Char
  ! will add the first character of a new string. If CASELESS is present
  ! and .true., compare caseless.
    integer, intent(out) :: STRING
    logical, intent(out) :: FOUND
    logical, optional, intent(in) :: CASELESS
    integer, optional, intent(in) :: DEBUG
    call LOOKUP_AND_INSERT_MAYBE ( STRING, FOUND, CASELESS, DEBUG, &
      & LOOKUPONLY=.true. )
  end subroutine LOOKUP
  ! ====================================     LOOKUP_AND_INSERT     =====
  subroutine LOOKUP_AND_INSERT ( STRING, FOUND, CASELESS, DEBUG )
  ! Look for the string built up by Add_Char.  If it is found return the
  ! position at which it was found in STRING, and FOUND = .true.  If it is
  ! not found, add it, return the position at which it was added in
  ! STRING, and FOUND = .false.  In any case, the next call to Add_Char
  ! will add the first character of a new string. If CASELESS is present
  ! and .true., compare caseless.
    integer, intent(out) :: STRING
    logical, intent(out) :: FOUND
    logical, optional, intent(in) :: CASELESS
    integer, optional, intent(in) :: DEBUG
    call LOOKUP_AND_INSERT_MAYBE ( STRING, FOUND, CASELESS, DEBUG, &
      & LOOKUPONLY=.false. )
  end subroutine LOOKUP_AND_INSERT
 ! ==============================================     NEW_LINE     =====
  subroutine NEW_LINE
  ! Skip the rest of the current line, set up so that get_char will
  ! read a new one.
    cur_pos = cur_end + 1
  end subroutine NEW_LINE
 ! =======================================     NUMERICAL_VALUE     =====
  integer function NUMERICAL_VALUE ( STRING )
  ! Return the numerical value of a string indexed by STRING
    integer, intent(in) :: STRING
    character(len=30) :: MY_CHAR
    call test_string ( string, 'Numerical_Value' )
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
  ! =========================================     OPEN_INCLUDE     =====
  subroutine Open_Include ( File_Name, Source, InFile, Stat )
  ! Check the inunit stack to make sure File_Name isn't there.  If not,
  ! put it on the stack and open the file.

    integer, intent(in) :: File_Name ! String index
    integer, intent(in) :: Source, InFile ! 256*line+column, string index
    integer, optional, intent(out) :: STAT

    logical :: Exist  ! Does file exist?  For use in INQUIRE statement
    integer :: MyFile ! File_Name or string index of directory/file
    ! Can't use String_Length here because it's not pure
    ! Assume maximum include path name is 254.
    character(len=strings(file_name) - strings(file_name-1)+255) :: MyName
    integer :: MyStat
    logical :: Opened ! "The file is already open"
    type(inunit_stack_t), allocatable :: Temp_Stack(:)

    myStat = 0
    if ( .not. allocated(inunit_stack) ) &
      & allocate ( inunit_stack(1:init_stack_size), stat=myStat )
    if ( myStat == 0 ) then
      if ( inunit_stack_ptr >= ubound(inunit_stack,1) ) then
        allocate ( temp_stack(1:2*ubound(inunit_stack,1)), stat=myStat )
        if ( myStat == 0 ) then
          temp_stack(:ubound(inunit_stack,1)) = inunit_stack
          call move_alloc ( temp_stack, inunit_stack )
        end if
      end if
    end if
    if ( myStat /= 0 ) then
      if ( present(stat) ) then
        stat = myStat
        return
      end if
      call io_error ( 'STRING_TABLE%OPEN_INCLUDE-E- Unable to allocate unit stack', &
                    &  myStat )
      call crash_burn
      return
    end if
    ! Look for the file name in the directory list
    inunit_stack_ptr = inunit_stack_ptr + 1
    if ( .not. allocated(includes) ) allocate ( includes(1:0) ) ! Assume it worked
    ! When the procedure pointer initialization works, we won't need this test
    if ( associated(find_a_file) ) then
      call find_a_file ( includes, file_name, exist, myName, &
                       & opened, inunit_stack(inunit_stack_ptr)%unit )
    else
      call find_file ( includes, file_name, exist, myName, &
                     & opened, inunit_stack(inunit_stack_ptr)%unit )
    end if
    if ( .not. exist ) then
      call display_string ( file_name, &
        & before= 'STRING_TABLE%OPEN_INCLUDE-E- the include file "', strip=.true. )
      call output ( '" could not be found.', advance='yes' )
      call crash_burn
      return
    end if
    if ( opened )  then
      ! Not allowed to open a file with more than one unit
      call loop
      return
    end if
    myFile = create_string ( trim(myName) )
    if ( any(inunit_stack(:inunit_stack_ptr)%file == myFile) ) then
      call loop
      return
    end if
    inunit_stack(inunit_stack_ptr)%file = myFile
    inunit_stack(inunit_stack_ptr)%line = source_line
    inunit_stack(inunit_stack_ptr)%column = source_column
    call open_input ( myName, unit=inunit_stack(inunit_stack_ptr)%unit, &
                    & stat=stat )
  contains
    subroutine Loop
      call display_string ( myFile, &
        & before='STRING_TABLE%OPEN_INCLUDE-E- Circular Include involving ', &
        & advance='yes' )
      call output ( source/256, before='Line ' )
      call output ( mod(source,256), before=', column ' )
      if ( inFile /= 0 ) call display_string ( inFile, before = ' in ' )
      call output ( '', advance='yes' )
      call crash_burn
    end subroutine Loop
  end subroutine Open_Include
  ! ===========================================     OPEN_INPUT     =====
  subroutine OPEN_INPUT ( FILE_NAME, UNIT, STAT )
  ! Open the file given by FILE_NAME for input.  If it can't be opened and
  ! STAT is present, return the status.  If it can't be opened and STAT is
  ! absent, ask the user for it.
    character(len=*), intent(in) :: FILE_NAME
    integer, intent(in) :: UNIT
    integer, optional, intent(out) :: STAT

    integer :: IOSTAT
    character(max(file_name_len,len(file_name))) :: MY_FILE

    my_file = file_name
    do
      open ( unit, file=my_file, status='OLD', access='SEQUENTIAL', &
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
      call output( 'Enter input file name: ', advance='yes')
      read ( *, '(a)', end=999 ) my_file
    end do
999 call crash_burn
  end subroutine OPEN_INPUT
 ! =========================================     STRING_LENGTH     =====
  integer function STRING_LENGTH ( STRING ) ! generic LEN
  ! Return the length of the string indexed by STRING
    integer, intent(in) :: STRING
    call test_string ( string, 'STRING_LENGTH' )
    string_length = strings(string) - strings(string-1)
  end function STRING_LENGTH
  ! ====================================     STRING_TABLE_SIZE     =====
  pure integer function STRING_TABLE_SIZE ()
  ! Return the allocated upper bound of the string table, or -1 if not allocated.
    if ( allocated(strings) ) then
      string_table_size = ubound(strings, 1)
    else
      string_table_size = -1
    end if
  end function STRING_TABLE_SIZE
 ! ============================================     UNGET_CHAR     =====
  subroutine UNGET_CHAR
  ! Back up the input, but not before the beginning of the current line
    cur_pos = max( cur_pos-1, strings(nstring+1) )
  end subroutine UNGET_CHAR

! =====     Private procedures    ======================================
  ! ------------------------------------     DOUBLE_CHAR_TABLE     -----
  subroutine DOUBLE_CHAR_TABLE ( STAT )
  ! Double the space for the character table
    integer, optional, intent(out) :: STAT
    integer :: MY_STAT
    character, allocatable :: NEW_CHAR(:)
    allocate ( new_char(2*ubound(char_table,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'stat_TABLE%DOUBLE_CHARS-E- Unable to allocate storage', stat )
      call crash_burn
    end if
    new_char(:ubound(char_table,1)) = char_table
    call move_alloc ( new_char, char_table )
  end subroutine DOUBLE_CHAR_TABLE
  ! ---------------------------------------     DOUBLE_STRINGS     -----
  subroutine DOUBLE_STRINGS ( STAT )
  ! Double the space for the string table
    integer, optional, intent(out) :: STAT
    integer :: MY_STAT
    integer, allocatable :: NEW_STRING(:)
    allocate ( new_string(0:2*ubound(strings,1)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'STRING_TABLE%DOUBLE_STRINGS-E- Unable to allocate storage', stat )
      call crash_burn
    end if
    new_string(0:ubound(strings,1)) = strings
    call move_alloc ( new_string, strings )
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
  ! ====================================     LOOKUP_AND_INSERT_MAYBE     =====
  subroutine LOOKUP_AND_INSERT_MAYBE ( STRING, FOUND, CASELESS, DEBUG, LOOKUPONLY )
  ! Look for the string built up by Add_Char.  If it is found return the
  ! position at which it was found in STRING, and FOUND = .true.  If it is
  ! not found, add it, return the position at which it was added in
  ! STRING, and FOUND = .false.  In any case, the next call to Add_Char
  ! will add the first character of a new string. If CASELESS is present
  ! and .true., compare caseless.
  ! LOOKUPONLY, if TRUE, modifies this a little--don't insert
    integer, intent(out) :: STRING
    logical, intent(out) :: FOUND
    logical, optional, intent(in) :: CASELESS
    integer, optional, intent(in) :: DEBUG
    logical, optional, intent(in) :: LOOKUPONLY ! DONTINSERT

    integer :: HASH_KEY  ! Integer derived from characters of STRING
    integer :: I         ! Subscript, loop inductor
    logical :: INSERT    ! Insert if not found
    integer :: LOC       ! Where HASH_KEY was found in HASH_TABLE
    integer :: myDEBUG   ! zero or DEBUG
    logical :: NOCASE    ! .false. (case sensitive) or CASELESS
    integer :: STATUS    ! Result, see HASH_LOOKUP

    nocase = .false.
    if ( present(caseless) ) nocase = caseless
    myDEBUG = 0
    if ( present(DEBUG) ) myDEBUG = DEBUG
    insert = .true.
    if ( present(LookupOnly) ) insert = .not. lookUpOnly

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
      ! call hash_lookup ( hash_key, hash_table, .true., loc, status )
      call hash_lookup ( hash_key, hash_table, insert, loc, status )
      if ( status == inserted ) then
        found = .false.
        if ( .not. insert ) then
          call clear_string
          return
        endif
        nstring = nstring + 1
        if ( nstring >= ubound(strings,1) ) then; call double_strings; end if
        hash_table(2,loc) = nstring
        string = nstring
        strings(nstring+1) = strings(nstring)
        if ( myDEBUG > 2 ) then
!        write ( *, * ) 'STRING_TABLE%LOOKUP_AND_INSERT-E- ', &
!                       'hash_key was not found in table'
!      write (*, *) 'hash key: ', hash_key
!      write (*, *) 'hash table keys: ', hash_table(2,:)
!      write (*, *) 'Compare', hash_table(2,loc), ': ', &
!      char_table(strings(hash_table(2,loc)-1)+1:strings(hash_table(2,loc))), &
!      ' to ',1, ': ', char_table(1:strings(2))
          call output('STRING_TABLE%LOOKUP_AND_INSERT-E- ', advance='no')
          call blanks(3)
          call output('hash_key was not found in table ', advance='yes')

          call output('hash key: ', advance='no')
          call output(hash_key, advance='yes')

          call output('hash table keys: ', advance='no')
          call output(hash_table(2,:), advance='yes')

          call output('Compare ', advance='no')
          call output(hash_table(2,loc), advance='no')
          call blanks(1)
          call output(':', advance='no')
          call blanks(1)
          call output(char_table(strings( &
           & hash_table(2,loc)-1)+1:strings(hash_table(2,loc) &
           & )), advance='yes')

          call output('to', advance='yes')

          call output('1 : ', advance='no')
          call blanks(1)
          call output(char_table(1:strings(2)), advance='yes')
        end if
        return
      end if
      if ( status /= hash_found ) then
!        write ( *, * ) 'STRING_TABLE%LOOKUP_AND_INSERT-E- ', &
!                       'Either the hash table is full'
!        write ( *, * ) 'or the program is using the hash software ', &
!                       'incorrectly.'
!        write ( *, * ) 'Unfortunately, the program doesn''t know how ', &
!                       'to increase the hash table size.'
!      write (*, *) 'status: ', status
          call output('STRING_TABLE%LOOKUP_AND_INSERT-E- ', advance='no')
          call blanks(3)
          call output('Either the hash table is full', advance='yes')
          call output('or the program is using the hash software' // &
          & 'incorrectly.', advance='yes')
          call output('Unfortunately, the program doesn''t know how' // &
          & 'to increase the hash table size.', advance='yes')

          call output('status: ', advance='no')
          call output(status, advance='yes')

          if(status == HASH_FULL) then
            call output( '(hash full)', advance='yes')
          else if(status == HASH_NOT_KEY) then
            call output( '(hash not key)', advance='yes')
          else if(status == HASH_BAD_LOC) then
            call output( '(hash bad loc)', advance='yes')
          else if(status == HASH_EMPTY) then
            call output( '(hash empty)', advance='yes')
          else
            call output( '(unrecognized hash error)', advance='yes')
          end if
!      write (*, *) 'hash key: ', hash_key
!      write (*, *) 'hash table keys: ', hash_table(2,:)
!      write (*, *) 'hash table size ', size(hash_table(2,:))
          call output('hash key: ', advance='no')
          call output(hash_key, advance='yes')
          call output('hash table keys: ', advance='no')
          call output(hash_table(2,:), advance='yes')
          call output('hash table size: ', advance='no')
          call output(size(hash_table(2,:)), advance='yes')
          if(status == HASH_FULL) then
             call output( 'You can probably fix this problem by ' // &
              & 'increasing hash_table_size,', advance='yes')
             call output( 'the last arg in the call to init_lexer ' // &
              & 'in your main program', advance='yes')
          end if
        call crash_burn
      end if
      ! The hash key matches; check whether the string does
      if ( myDEBUG > 1 ) then
! write (*, *) 'Compare', hash_table(2,loc), ': ', &
! char_table(strings(hash_table(2,loc)-1)+1:strings(hash_table(2,loc))), &
! ' to ', nstring+1, ': ', char_table(strings(nstring)+1:strings(nstring+1))
          call output('Compare ', advance='no')
          call output(hash_table(2,loc), advance='no')
          call blanks(1)
          call output(':', advance='no')
          call blanks(1)
          call output(char_table(strings( &
           & hash_table(2,loc)-1)+1:strings(hash_table(2,loc) &
           & )), advance='yes')

          call output('to', advance='yes')

          call output(nstring+1, advance='no')
          call blanks(1)
          call output(':', advance='no')
          call output(char_table(strings(nstring)+1:strings(nstring+1)), advance='yes')

          call output ( 'Compare ', advance='no' )
          call output ( hash_table(2,loc) ); call output ( ': ', advance='no')
          call display_string ( hash_table(2,loc) )
          call output ( ' to ' )
          call output ( nstring+1 ); call output ( ': ')
        !  call display_string ( nstring+1, advance='no' )
          call output ( char_table(strings(nstring)+1:strings(nstring+1)), advance='no' )
          if ( nocase ) call output ( ' caseless')
          call output ( '', advance='yes' )
      end if
      string = hash_table(2,loc)
      found = compare_strings ( string, nstring+1, caseless ) == 0
      if ( found ) then
        call clear_string
        return
      end if
    end do
  end subroutine LOOKUP_AND_INSERT_MAYBE
  ! ------------------------------------------     TEST_STRING     -----
  subroutine TEST_STRING ( STRING, ROUTINE, IERR )
  ! Test whether STRING is within bounds.  If not, use ROUTINE to emit
  ! an error message.
  ! Unless IERR is present, in which case set IERR and return
    integer, intent(in) :: STRING
    integer, intent(out), optional :: IERR
    character(len=*), intent(in) :: ROUTINE
    if ( .not. isStringInTable(string) ) then
      if ( present(ierr) ) then
        ierr=1
        return
      else
        call output( 'STRING_TABLE%' // trim(routine) // &
         & '-E- ', advance='yes')
        call output( 'String index ', advance='no')
        call output( string, advance='no')
        call output( ' not in 1 .. ', advance='no')
        call output( nstring, advance='yes')
        call crash_burn
      end if
    end if
    if ( present(ierr) ) ierr=0
    return
  end subroutine TEST_STRING

  logical function isStringInTable ( string, indices )
  ! Test whether STRING is within bounds.
    integer, intent(in) :: STRING
    integer, dimension(:), optional, intent(in) :: indices
    ! Executable
    if ( .not. present(indices) ) then
      isStringInTable = .not. ( string < 1 .or. string > nstring ) 
    else
      isStringInTable = .not. ( string < 1 .or. string > size(indices) )
      if ( isStringInTable ) &
        & isStringInTable = .not. &
        & ( indices(string) < 1 .or. indices(string) > nstring ) 
    endif
  end function isStringInTable

  ! Add another unit to read from
  subroutine AddInUnit (inunit)
    integer, intent(in) :: inunit
    integer, dimension(:), pointer :: temp_inunits
    integer :: status

    ! Executables
    if (associated(inunit_list)) then
      allocate(temp_inunits(size(inunit_list) + 1), stat=status)
    else
      allocate(temp_inunits(1), stat=status)
    end if

    if (status /= 0) &
      call PRINTITOUT(  MLSMSG_Allocate // "temp_inunits", MLSMSG_Error, &
        & 1, exitStatus=1 )

    if (associated(inunit_list)) then
      temp_inunits(1:size(inunit_list)) = inunit_list
      deallocate(inunit_list, stat=status)
      if (status /= 0) &
        call PRINTITOUT( MLSMSG_DeAllocate // "inunit_list", &
        & 1, exitStatus=1 )
    end if

    temp_inunits(size(temp_inunits)) = inunit
    inunit_list => temp_inunits
    inunit_counter = 1
  end subroutine AddInUnit

  subroutine SetNextInUnit
    integer :: status

    if (.not. associated(inunit_list)) return

    if (inunit_counter == size(inunit_list)) then
      inunit_counter = 0
      deallocate(inunit_list, stat=status)
      if (status /= 0) &
        call PRINTITOUT( MLSMSG_DeAllocate // "inunit_list", &
        & 1, exitStatus=1 )
      inunit_list => NULL()
    else
      inunit_counter = inunit_counter + 1
    end if

  end subroutine SetNextInUnit

  subroutine Init_String_Table
    source_line = 0
    ! no need to initialize source_column
    at_eof = .false.
    cur_end = 0
    cur_pos = 1
    inunit_counter = 0
    ! inunit_list might not be null here if AddInUnit 
    ! has been used before this call
    if (associated(inunit_list)) deallocate(inunit_list)
    inunit_list => NULL()
  end subroutine Init_String_Table

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module STRING_TABLE

! $Log$
! Revision 2.54  2022/02/18 00:14:37  pwagner
! Fixed a bug in Get_String when end is present
!
! Revision 2.53  2021/11/03 23:48:44  pwagner
! Added Get_String_Char
!
! Revision 2.52  2019/07/31 20:03:24  vsnyder
! Add Dump_Char_Table, which might be useful for debugging
!
! Revision 2.51  2016/10/04 20:55:25  vsnyder
! Make Find_File private because ifort 17 crashes with it public
!
! Revision 2.50  2016/08/16 23:15:13  vsnyder
! Add Pos and Width arguments to Display_String
!
! Revision 2.49  2016/05/17 00:05:53  pwagner
! Added optional arg indices to isStringInTable
!
! Revision 2.48  2016/02/25 00:54:39  vsnyder
! Bump inunit_stack_ptr early enough
!
! Revision 2.47  2016/02/12 21:08:09  pwagner
! Reinitialize nstring in DESTROY_STRING_TABLE
!
! Revision 2.46  2016/01/29 01:10:49  vsnyder
! Add procedure pointer to get full name of include from PCF
!
! Revision 2.45  2015/07/14 23:15:43  pwagner
! Added isStringInTable so outside modules can check first
!
! Revision 2.44  2014/03/15 00:08:03  vsnyder
! Return the allocated upper bound from string_table_size, unless if it's
! not allocated return -1.  Simplify routines to double the character and
! string tables.  Some cannonball polishing.
!
! Revision 2.43  2013/12/12 02:00:15  vsnyder
! Change type of debug from logical to integer
!
! Revision 2.42  2013/10/01 02:12:02  vsnyder
! Check current directory first for includes.  Handle the case of a prefix
! in the directory list not ending with "/" correctly.
!
! Revision 2.41  2013/09/27 22:35:36  pwagner
! Supplied lower bounds on allocates in Open_Include to mollify NAG
!
! Revision 2.40  2013/09/25 02:05:15  vsnyder
! Even more improved include loop detector
!
! Revision 2.39  2013/09/25 01:02:10  vsnyder
! Improved include loop detector
!
! Revision 2.38  2013/09/24 23:07:16  vsnyder
! Add Includes, Open_Include, Find_File
!
! Revision 2.37  2013/09/21 00:20:56  pwagner
! Initialize NSTRING so HOW_MANY_STRINGS can be checked at start
!
! Revision 2.36  2013/08/28 00:36:45  pwagner
! Moved more stuff from MLSMessage down to PrintIt module
!
! Revision 2.35  2013/06/12 02:15:38  vsnyder
! Cruft removal
!
! Revision 2.34  2012/01/05 01:12:43  pwagner
! Added Lookup sub routine; lookup_and_insert loses optional parameter lookuponly
!
! Revision 2.33  2011/10/11 16:57:51  honghanh
! Fix a bug in get_string to return an error code when no_error is .true.
!
! Revision 2.32  2011/08/27 13:25:13  honghanh
! Fix a memory leak bug in string_table regarding inunit_list
!
! Revision 2.31  2011/07/22 18:29:22  vsnyder
! Produce a message in "string_text" instead of an empty string if the
! "string" argument to get_string is out of bounds and noError is present
! and true.
!
! Revision 2.30  2010/08/05 17:45:34  honghanh
! Adding subroutine init_string_table in string_table module
!
! Revision 2.29  2010/05/23 03:05:26  honghanh
! Use an inunitList instead of an inunit to read from multiple l2cf
!
! Revision 2.28  2010/05/14 02:15:47  vsnyder
! Calculate length of STRING correctly if SUBSTRING is quoted
!
! Revision 2.27  2010/05/07 02:23:46  vsnyder
! Add STRIP optional argument to INDEX
!
! Revision 2.26  2010/04/30 22:18:51  vsnyder
! Add CASELESS to INDEX_IN_STRING
!
! Revision 2.25  2010/04/14 03:13:52  vsnyder
! Add INDEX_IN_STRING, with a generic INDEX, working similarly to intrinsic
! INDEX but with string indices instead of strings, and without the BACK
! argument.  Add START and END arguments for GET_STRING.  Add a LEN generic
! for STRING_LENGTH.
!
! Revision 2.24  2009/10/01 19:42:14  vsnyder
! Simplify error testing in Get_String, improve comments
!
! Revision 2.23  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.22  2006/07/27 03:56:27  vsnyder
! Make String_Table_Size pure, so it can be a specification function
!
! Revision 2.21  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.20  2004/10/30 00:19:56  vsnyder
! Revised the CVS stuff
!
! Revision 2.19  2004/10/21 00:37:25  vsnyder
! Add Display_String_List and 'before' argument to Display_String
!
! Revision 2.18  2004/08/19 00:14:03  pwagner
! crash_burn instead of stop
!
! Revision 2.17  2004/05/20 23:53:10  vsnyder
! Handle no-luns-available error condition
!
! Revision 2.16  2003/05/14 01:39:30  vsnyder
! Add Dump_String_Table
!
! Revision 2.15  2002/10/08 00:09:14  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.14  2002/08/21 20:38:00  vsnyder
! Add Create_String
!
! Revision 2.13  2002/05/23 20:58:33  vsnyder
! Detect errors while reading input
!
! Revision 2.12  2002/01/09 23:49:38  pwagner
! Each print became call output
!
! Revision 2.11  2001/10/12 23:08:14  pwagner
! New ierr option to prevent (unwanted) stoppings
!
! Revision 2.10  2001/10/04 00:14:21  pwagner
! Some more messages if lookup_and_insert fails
!
! Revision 2.9  2001/06/06 17:30:15  pwagner
! DEBUG optional arg to lookup..
!
! Revision 2.8  2001/05/15 16:44:41  livesey
! Added noError argument to get_string
!
! Revision 2.7  2001/04/20 17:43:35  vsnyder
! OOPS -- previous commit was premature -- forgot to declare a variable
!
! Revision 2.6  2001/04/20 17:40:49  vsnyder
! Add 'Destroy...' subroutines
!
! Revision 2.5  2001/04/11 21:54:26  vsnyder
! Put in more bounds checking
!
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
