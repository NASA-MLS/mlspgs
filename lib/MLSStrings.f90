! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSStrings               ! Some low level string handling stuff
!=============================================================================

  implicit none
  private

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

!
! This module contains some low level string handling stuff for mls

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! Capitalize         tr[a-z] -> [A-Z]
! CompressString     Removes all leading and embedded blanks
! Count_words        Counts the number of space-separated words in a string
! Depunctuate        Replaces punctuation with blanks
! Hhmmss_value       Converts 'hh:mm:ss' formatted string to a real r8
!                    (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
! Ints2Strings       Converts an array of integers to strings using "char" ftn
! LinearSearchStringArray     
!                    Finds string index of substring in array of strings
! LowerCase          tr[A-Z] -> [a-z]
! ReadCompleteLineWithoutComments     
!                    Knits continuations, snips comments
! ReadIntsFromChars  Converts an array of strings to ints using Fortran read
! ReformatDate       Turns 'yyyymmdd' -> 'yyyy-mm-dd'; or more general format
! ReformatTime       Turns 'hhmmss.sss' -> 'hh:mm:ss'
! Reverse            Turns 'a string' -> 'gnirts a'
! SplitWords         Splits 'first, the, rest, last' -> 'first', 'the, rest', 'last'
! Strings2Ints       Converts an array of strings to ints using "ichar" ftn
! WriteIntsToChars   Converts an array of ints to strings using Fortran write
! === (end of toc) ===

! === (start of api) ===
! char* Capitalize (char* str)
! char* CompressString (char* str)
! int count_words (char* str)
! char* depunctuate (char* str)
! ints2Strings (int ints(:,:), char* strs(:))
! int LinearSearchStringArray (char* list(:), char* string, 
!   [log caseInsensitive, [log testSubstring], [log listInString])
! char* LowerCase (char* str)
! readIntsFromChars (char* strs(:), int ints(:), char* forbiddens)
! ReadCompleteLineWithoutComments (int unit, char* fullLine, [log eof], &
!       & [char commentChar], [char continuationChar])
! char* ReformatDate (char* date, [char* fromForm], [char* toForm])
! char* ReformatTime (char* time, [char* form])
! char* Reverse (char* str)
! SplitWords (char *line, char* first, char* rest, [char* last], &
!       & [log threeWay], [char* delimiter])
! strings2Ints (char* strs(:), int ints(:,:))
! writeIntsToChars (int ints(:), char* strs(:))
! Many of these routines take optional arguments that greatly modify
! their default operation

! Warnings: 
! (1) in the routine LinearSearchStringArray
! the input arguments include an array of strings;
! This array is of assumed-size
! I.e., all elements from array(1:size(array)) are relevant
! Therefore in calling one of these you probably need to use the format
!   call SortArray(myArray(1:mySize), ..
! to avoid operating on undefined array elements
! === (end of api) ===

  public :: Capitalize, CompressString, count_words, &
   & depunctuate, hhmmss_value, &
   & ints2Strings, LinearSearchStringArray, &
   & LowerCase, &
   & ReadCompleteLineWithoutComments, readIntsFromChars, &
   & reformatDate, reformatTime, Reverse, &
   & SplitWords, strings2Ints, &
   & writeIntsToChars

  ! hhmmss_value
  integer, public, parameter :: INVALIDHHMMSSSTRING = 1
  ! strings2Ints
  integer, public, parameter :: LENORSIZETOOSMALL=-999
  
  integer, parameter :: YEARMAX = 4999  ! Conversion invalid after 4999 AD
  ! The following arrays contains the maximum permissible day for each month
  ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
  integer, dimension(-1:12), parameter :: DAYMAXLY = (/ &
    & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)
  integer, dimension(-1:12), parameter :: DAYMAXNY = (/ &
    & 365, 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)
  character(len=*), dimension(12), parameter :: MONTHNAME = (/ &
    & 'January  ', 'February ', 'March    ', 'April    ', 'May      ', &
    & 'June     ', 'July     ', 'August   ', 'September', 'October  ', &
    & 'November ', 'December '/)
contains

  ! -------------------------------------------------  CAPITALIZE  -----
  elemental function Capitalize (STR) result (OUTSTR)
    ! takes a-z and replaces with A-Z 
    ! leaving other chars alone
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=len(str)) :: OUTSTR

    !----------Local vars----------!
    integer :: I, ICODE
    integer, parameter :: OFFSET=iachar("A")-iachar("a")
    !----------Executable part----------!
    outstr=str

    do i=1, len(str)
       icode=iachar(outstr(i:i))
       if ( icode >=iachar("a") .and. icode <= iachar("z")) then
          outstr(i:i)=achar(icode+offset)
       end if
    end do

  end function Capitalize

  ! ---------------------------------------------  CompressString  -----
  FUNCTION CompressString (str) RESULT (outstr)

    ! Removes all leading and embedded blanks from a string
    !--------Argument--------!

    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    INTEGER :: i, n

    !----------Executable part----------!

    outstr = " "
    n = 0
    DO i = 1, LEN(str)
       IF (str(i:i) /= " ") THEN
          n = n + 1
          outstr(n:n) = str(i:i)
       END IF
    END DO

  END FUNCTION CompressString

  ! ------------------------------------------------  COUNT_WORDS  -----
  FUNCTION count_words (str) RESULT (no_of_words)
    ! This counts the number of words in a string 
    ! For our purposes, words consist of any non-space characters
    ! and are separated by one or more spaces
    ! -----Added by HCP-------- 
    !--------Argument--------!
    CHARACTER (len=*), INTENT(in) :: str
    !---------result---------!
    INTEGER::no_of_words
    !-----local-variables------!
    INTEGER::j
    !-------Executable-code----!
    IF (str(1:1) /= " ") THEN
       no_of_words=1
    ELSE
       no_of_words=0
    END IF
    DO j = 2, LEN(str)
       IF(str(j:j) /= " " .AND. str(j-1:j-1) == " ") THEN
          no_of_words=no_of_words+1
       END IF
    END DO
  END FUNCTION count_words

  ! This one converts a string to all upper case (taken from HCP routine
  ! of same name) (Except that HCP can spell capitalise, that is. Fnord.)

  ! ------------------------------------------------  DEPUNCTUATE  -----
  Function Depunctuate ( str ) result ( outstr )
    ! Function that removes punctuation and replaces with blanks
    ! Added by HCP. This depends on the native character set being 
    ! ASCII. 
    !--------Argument--------!
    character(len=*),intent(in) :: str
    character(len=len(str)) :: outstr
    !----------Local vars----------!
    integer :: i, icode
    !----------Executable part----------!
    outstr=str
    do i = 1 ,len(str)
        icode=iachar(str(i:i))
        if(  (icode >= 33 .and. icode <= 47).or.&
             (icode >= 58 .and. icode <=64).or. &
             (icode >= 91 .and. icode <=96).or. &
             (icode >= 123)) then
           outstr(i:i)=" "
        end if
    end do

  end Function Depunctuate

  ! ------------------------------------------------  HHMMSS_value  -----
  function HHMMSS_value ( str, ErrTyp, separator, strict ) result ( value )
    use MLSCommon, only: R8
    ! Function that returns the value in seconds of a string 'hh:mm:ss'
    ! where the field separator ':' divides the string into two
    ! integer-like strings 'hh' and 'mm', as well as one float-like
    ! string 'ss' which may have a decimal point plus fractional part
    ! E.g., ss=59.9999

    ! Requires 0 <= hh <= 24
    ! Requires 0 <= mm < 60
    ! Requires 0. <= ss < 60.

    ! Returns ErrTyp == 0 unless an error occurs

    ! Lenient wrt utc and non-compliant formats:
    ! ignores chars in front of 'hh' and a terminal,
    ! non-numerical char: e.g., '2000-01-01T00:00:00.000000Z'
    ! will be treated the same as '00:00:00.0000000'

    ! If given optional arg strict and it's true, not lenient
    ! i.e., non-compliant str always returns non-zero ErrTyp

    ! If not strict, some of the fields can be null, e.g. 12:00 and 12:00:
    ! and 12::0 and ::12 are allowed.

    ! If given optional arg separator, uses separator as field separator

    ! Useful to allow an added way to input time

    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Arguments--------!
    character(len=*), intent(in) :: Str
    integer, intent(out) :: ErrTyp
    character(len=1), intent(in), optional :: Separator
    logical, intent(in), optional :: Strict
    real(r8) :: value
    !----------Locals----------!
    character(len=1), parameter :: Colon=':'
    character(len=1) :: MyColon
    character(len=*), parameter :: Digits='0123456789.'
    integer :: HValue, MValue
    character(len=len_trim(str)+1) :: MyStr
    logical :: MyStrict
    integer :: H1, Sep1, Sep2, S2  ! Indices in STR
    !----------Executable part----------!

    myColon = colon
    if ( present(separator) ) myColon = separator

    myStrict = .false.
    if ( present(strict) ) myStrict = strict

    s2 = len_trim(str)
    value = 0.0_r8
    errTyp = INVALIDHHMMSSSTRING
    if ( s2 == 0 ) then
      if ( .not. myStrict ) errTyp = 0
      return
    end if
    myStr = str(:s2)
    if ( verify(myStr(s2:s2),digits) /= 0 ) then ! Junk at the end
      if ( myStrict ) return
    else
      s2 = s2 + 1
    end if
    myStr(s2:s2) = '/' ! list-directed I/O terminator, might replace junk at end

    sep1 = index(myStr,myColon)
    if ( sep1 == 0 ) then
      if ( myStrict ) return
      sep1 = s2
    else
      myStr(sep1:sep1) = ',' ! list-directed I/O separator
    end if

    h1 = verify(myStr(:sep1-1),digits(1:10),back=.true.) ! '.' not allowed
    if ( h1 /= 0 .and. myStrict ) return ! Junk before hours
    if ( verify(myStr(h1+1:sep1-1),digits(1:10)) /= 0 ) return ! catch blank/comma

    sep2 = index(myStr(sep1+1:),myColon) + sep1
    if ( sep2 == sep1 ) then
      if ( myStrict ) return
      sep2 = s2
    else
      myStr(sep2:sep2) = ',' ! list-directed I/O separator
    end if
    if ( verify(myStr(sep1+1:sep2-1),digits(1:10)) /= 0 ) return ! catch blank/comma
    if ( verify(myStr(sep2+1:s2-1),digits) /= 0 ) return ! catch blank/comma

    if ( myStrict .and. & ! check for null fields
      & ( h1+1 == sep1 .or. sep1+1 == sep2 ) ) return

    ! myStr(h1+1:sep1-1) is hours, myStr(sep1+1:sep2-1) is minutes,
    ! myStr(sep2+1:s2) is seconds
    hvalue=0; mvalue=0 ! Fortran doesn't update null fields
    read ( myStr(h1+1:s2), *, iostat=errTyp ) hvalue, mvalue, value
    if ( errTyp /= 0 ) return

    errTyp = INVALIDHHMMSSSTRING
    if ( hvalue < 0 .or. hvalue >= 24 .or. &
      &  mvalue < 0 .or. mvalue >= 60 .or. &
      &  value < 0 .or. value >= 60 ) return

    errTyp = 0
    value = value + 60 * (mvalue + 60 * hvalue )

  end function HHMMSS_value

  ! --------------------------------------------------  Int2Strings  -----
  subroutine Ints2Strings ( ints, strs )
    ! takes an array of integers and returns string array
    ! using "char"
	 ! Useful due to bug in toolbox swrfld
	 !
	 ! See also strings2Ints
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    integer, intent(in), dimension(:,:) ::          ints
    character (len=*), intent(out), dimension(:) :: strs

    !----------Local vars----------!
    integer :: i, substr, strLen, arrSize
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   strLen = MIN(len(strs(1)), size(ints, dim=1))
   arrSize = MIN(size(strs), size(ints, dim=2))
   strs = ' '
   if ( strLen <= 0 .or. arrSize <= 0 ) return
   do i=1, arrSize
      do substr=1, strLen
         strs(i)(substr:substr) = achar(ints(substr, i))
      end do
   end do

  end subroutine Ints2Strings

  ! ------------------------------------  LinearSearchStringArray  -----

  ! This routine does a simple linear search for a string in an array.
  ! If the case insensitive flag is set the strings are capitalized first.
  ! If the test substring flag is set, the string is tested as a substring.
  ! If the listInString flag is set, the array list is tested as substrings
  !  against the string; otherwise, the string is tested as a substring
  !  against the array list.
  ! If the string is not found, 0 is returned

  FUNCTION LinearSearchStringArray (list, string, caseInsensitive, &
       & testSubstring, listInString) RESULT (sindex)

    ! Dummy arguments
    CHARACTER (LEN=*), DIMENSION(:) :: list
    CHARACTER (LEN=*) :: string
    LOGICAL, intent (in), OPTIONAL :: caseInsensitive
    LOGICAL, intent (in), OPTIONAL :: testSubstring
    LOGICAL, intent (in), OPTIONAL :: listInString


    ! Function result
    INTEGER :: sindex   ! matching string index (0 = not found)

    ! Local variables
    INTEGER :: i
    LOGICAL :: useCaseInsensitive
    LOGICAL :: testForSubstring
    LOGICAL :: testForList
    LOGICAL :: found

    ! Executable code

    IF (PRESENT(caseInsensitive)) THEN
       useCaseInsensitive = caseInsensitive
    ELSE
       useCaseInsensitive = .FALSE.
    END IF

    IF (PRESENT(testSubstring)) THEN
       testForSubstring = testSubstring
    ELSE
       testForSubstring = .FALSE.
    END IF

    IF (PRESENT(listInString)) THEN
       testForList = listInString
    ELSE
       testForList = .FALSE.
    END IF

    found = .FALSE.
    sindex = 0

    ! Put the conditional outside the loop for speed (not that it will make 
    ! much difference for strings)

    IF (useCaseInsensitive) THEN
       linSearchStringInsens: DO i = 1, SIZE(list)
          IF (testForSubstring) THEN
             IF (testForList) THEN
                found = (INDEX (Capitalize(TRIM(string)), &
                     & Capitalize(TRIM(list(i)))) /= 0)
             ELSE
                found = (INDEX (Capitalize(TRIM(list(i))), &
                     & Capitalize(TRIM(string))) /= 0)
             END IF
          ELSE
             found = (Capitalize(TRIM(list(i))) == Capitalize(TRIM(string)))
          END IF
          IF (found) THEN
             sindex = i
             EXIT linSearchStringInsens
          END IF
       END DO linSearchStringInsens
    ELSE
       linSearchStringSens: DO i = 1, SIZE(list)
          IF (testForSubstring) THEN
             IF (testForList) THEN
                found = (INDEX (TRIM(string), TRIM(list(i))) /= 0)
             ELSE
                found = (INDEX (TRIM(list(i)), TRIM(string)) /= 0)
             END IF
          ELSE
             found = (TRIM(list(i)) == TRIM(string))
          END IF
          IF (found) THEN
             sindex = i
             EXIT linSearchStringSens
          END IF
       END DO linSearchStringSens
    END IF

  END FUNCTION LinearSearchStringArray

  ! --------------------------------------------------  LowerCase  -----
  elemental function LowerCase (str) result (outstr)
    ! takes A-Z  and replaces with a-z
    ! leaving other chars alone
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=len(str))      :: OUTSTR

    !----------Local vars----------!
    integer            :: i, icode
    integer, parameter :: offset=IACHAR("a")-IACHAR("A")
    !----------Executable part----------!
    outstr=str

    DO i = 1, LEN(str)
       icode=IACHAR(outstr(i:i))
       IF ( icode >=IACHAR("A") .AND. icode <= IACHAR("Z")) THEN
          outstr(i:i)=achar(icode+offset)
       END IF
    END DO

  END FUNCTION LowerCase

  ! ----------------------------  ReadCompleteLineWithoutComments  -----

  ! This funtion reads a line or set of lines from a text file and returns a
  ! string giving the full command with continuation lines joined, and comments
  ! removed.

  ! EOF can be returned if requested

  ! Note that this doesn't consider quotation marks, comments within quoted
  ! strings are considered comments, and continuation marks can apply within
  ! quoted strings.  Later versions of this routine may be more intelligent.

  SUBROUTINE ReadCompleteLineWithoutComments(unit,fullLine,eof, &
       & commentChar,continuationChar)

    ! Dummy arguments

    INTEGER, INTENT(IN) :: unit ! Input file unit
    ! fullLine changed to intent InOut by HCP. Some (but not all) 
    ! F90 compilers won't let this be intent(out) because the declaration
    ! of inputLine makes use of the length of fullLine even if its _contents_
    ! are immaterial
    CHARACTER(LEN=*), INTENT(INOUT) :: fullLine ! Output line
    CHARACTER(LEN=*), OPTIONAL :: commentChar
    CHARACTER(LEN=*), OPTIONAL :: continuationChar
    LOGICAL, INTENT(OUT), OPTIONAL :: eof ! Set if got to eof

    ! Local variables

    INTEGER :: ioInfo           ! IOSTAT result
    CHARACTER(LEN=LEN(fullLine)) :: inputLine ! One line from file
    INTEGER :: commentStart     ! Start of a comment in line
    INTEGER :: lastChar         ! Last character position in line
    INTEGER :: gotContinuation  ! 1 if continuation needed, 0 if not
    LOGICAL :: firstLine        ! A correction to be applied

    CHARACTER(LEN=1) :: useCommentChar
    CHARACTER(LEN=1) :: useContinuationChar

    ! Executable code

    ! Set default values for arguments
    
    IF (.NOT. PRESENT(commentChar)) THEN
       useCommentChar=";"
    ELSE
       useCommentChar=commentChar
    END IF

    IF (.NOT. PRESENT(continuationChar)) THEN
       useContinuationChar="$"
    ELSE
       useContinuationChar=continuationChar
    END IF

    ! Set up for loop

    fullLine=""
    firstLine=.TRUE.
    IF (PRESENT(eof)) eof=.FALSE.

    readLoop: DO

       ! Try to read a line of text

       READ (UNIT=unit,FMT="(a)",IOSTAT=ioInfo) inputLine
       IF (ioInfo /= 0) THEN 
          IF (PRESENT(eof)) eof=.TRUE.
          EXIT readLoop
       END IF

       ! Now we look for the start of a comment and remove any following text
       ! from this line.

       commentStart=INDEX(inputLine,useCommentChar)
       IF (commentStart /= 0) inputLine=inputLine(1:commentStart-1)

       ! See if the last non blank character is a contination

       lastChar=LEN_TRIM(inputLine)
       ! if bloc inserted by HCP because inputline(lastchar:lastchar:) 
       ! caused errors with array bounds checking turned on with
       ! some f90 compilers.
       if (lastChar > 0) then 
          gotContinuation=INDEX(inputLine(lastChar:lastChar),&
               useContinuationChar)
       else
          gotContinuation=0
       end if
       ! Concatenate this with what we have so far, make sure there's an extra
       ! space there though (not for first line though)

       ! If block inserted 5 Sept. 2000 by HCP to prevent an out-of-range 
       ! error when the input line has 0 length and you have bounds-checking
       ! turned on
       if(LEN_TRIM(inputLine) > 0) then 
          inputLine=inputLine(1:LEN_TRIM(inputLine)-gotContinuation)
       end if
       IF (firstLine) THEN
          fullLine=inputLine
          firstLine=.FALSE.
       ELSE
          fullLine=fullLine(1:LEN_TRIM(fullLine)+1)//inputLine
       END IF
       
       ! If we have a continuation mark, or a blank line then keep reading
       IF ((gotContinuation==0).AND.(LEN_TRIM(fullLine) /= 0)) EXIT readLoop
    END DO readLoop

    ! Do a final trim and exit

    fullLine=TRIM(ADJUSTL(fullLine))

  END SUBROUTINE ReadCompleteLineWithoutComments

  ! --------------------------------------------------  readIntsFromChars  -----
  SUBROUTINE readIntsFromChars (strs, ints, forbiddens)
    ! takes an array of strings and returns integer array
    ! using Fortran "read"
    ! If any element of string array is blank or contains one of forbiddens
    ! the corresponding element of ints is left undefined
	 ! Not useful yet
	 !
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    CHARACTER (LEN=*), INTENT(in), dimension(:) ::   strs
    integer, intent(out), dimension(:)          ::   ints
    CHARACTER (LEN=*), INTENT(in), optional     ::   forbiddens
    character(len=40)                           ::   myForbiddens

    !----------Local vars----------!
    INTEGER :: i, j, arrSize
    LOGICAL :: leave_undef
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   arrSize = MIN(size(strs), size(ints))
   if ( arrSize <= 0 ) then
     ints = LENORSIZETOOSMALL
     return
   endif
   if ( present(forbiddens) ) then
     myForbiddens = adjustl(forbiddens)
   else
     myForbiddens = ' '
   endif
   do i=1, arrSize
      leave_undef = (strs(i) == ' ')
      if ( myForbiddens /= ' ' ) then
        do j=1, len(trim(myForbiddens))
           leave_undef = leave_undef &
            & .or. &
            & ( &
            &    index(strs(i), myForbiddens(j:j)) > 0 &
            &  .and. &
            &    myForbiddens(j:j) /= ' ' &
            & )
        enddo
      endif
      if ( .not. leave_undef ) read(strs(i), *) ints(i)
   enddo

  END SUBROUTINE readIntsFromChars

  ! --------------------------------------------------  reFormatDate  -----
  function reFormatDate(date, fromForm, toForm) result(reFormat)
    ! Reformat yyyymmdd as yyyy-mm-dd
    ! Wouldn't it be clever to allow an optional string arg defining
    ! the output format; E.g. 'dd M yyyy' for '03 September 2005'
    ! or 'yyyy-doy' for '2005-d245'
    ! And yet another optional string to hold the input format
    ! in case it wasn't yyyymmdd?
    ! Args
    character(len=*), intent(in)            :: date
    character(len=len(date)+24)             :: reFormat
    character(len=*), optional, intent(in)  :: fromform ! output format
    character(len=*), optional, intent(in)  :: toform   ! input format
    ! Internal variables
    character(len=1), parameter             :: ymSpacer = '-'
    character(len=1), parameter             :: mdSpacer = '-'
    integer                                 :: i ! format string index
    integer                                 :: j ! date string index
    integer                                 :: doy
    character(len=4)                        :: doyString
    character(len=4)                        :: yyyyString
    integer                                 :: month
    character(len=len(date)+24)             :: tempFormat
    logical                                 :: inputWasDoy
    character(len=1), parameter             :: SPACESUBSTITUTE = '?'
    ! Executable
    tempFormat = date
    if ( present(fromform) ) then
      ! print *, 'fromForm: ', trim(fromForm)
      if ( len_trim(fromform) > 0 ) then
        inputWasDoy = .false.
        tempFormat = ' '
        i = 0
        j = 0
        do
          if ( i >= len_trim(fromform) ) exit
          i = i + 1
          j = j + 1
          ! print *, fromform(i:i)
          select case (fromform(i:i))
          case ('y')
            tempFormat(1:4) = date(j:j+3)
            i = i + 3
            j = j + 3
          case ('m')
            tempFormat(5:6) = date(j:j+1)
            i = i + 1
            j = j + 1
          case ('M')
            ! print *, j, trim(date), trim(date(j:))
            month = monthNameToNumber(date(j:))
            if ( month < 1 .or. month > 12 ) then
              reFormat = 'month name uncrecognized'
              return
            endif
            j = j + len_trim(MONTHNAME(month)) - 1
            write(tempFormat(5:6),'(i2.2)') month
            ! print *, j, trim(date), trim(date(j:))
          case ('d')
            ! print *, fromform(i:i+1)
            if ( fromform(i:i+1) == 'dd' ) then
              tempFormat(7:8) = date(j:j+1)
              i = i + 1
              j = j + 1
            else
              inputWasDoy = .true.
              doyString = date(j:j+3)
              i = i + 2
              j = j + 2
            endif
          case default
            ! i = i + 1
          end select
        enddo
        if ( inputWasDoy ) then
          ! Need to convert from yyyydoy to yyyymmdd
          ! print *, ' Need to convert from yyyydoy to yyyymmdd'
          yyyyString = tempFormat(1:4)
          call yyyydoy_to_yyyymmdd_str(yyyyString, doyString(2:4), tempFormat)
          ! print *, yyyyString, doyString(2:4), tempFormat
        endif
      endif
    endif
    ! print *, tempFormat
    reFormat = tempFormat(1:4) // ymSpacer // tempFormat(5:6) // mdSpacer // tempFormat(7:8)
    if ( .not. present(toform) ) return
    if ( len_trim(toform) < 1 ) return
    reFormat = ' '
    i = 0
    do
      if ( i >= len_trim(toform) ) exit
      i = i + 1
      select case (toform(i:i))
      case ('y')
        reFormat = trim(reFormat) // tempFormat(1:4)
        i = i + 3
      case ('m')
        reFormat = trim(reFormat) // tempFormat(5:6)
        i = i + 1
      case ('M')
        ! print *, 'Now trying to convert to month name ', trim(tempFormat), ' ', tempFormat(5:6)
        read(tempFormat(5:6), *) month
        reFormat = trim(reFormat) // trim(MONTHNAME(month))
      case ('d')
        if ( toform(i:i+1) == 'dd' ) then
          reFormat = trim(reFormat) // tempFormat(7:8)
          i = i + 1
        else
          call yyyymmdd_to_doy_str(tempFormat, doy)
          write(doyString, '(a1, i3.3)') 'd', doy
          reFormat = trim(reFormat) // doyString
          i = i + 2
        endif
      case (' ')
        ! Foolish me--it can't handle spaces because of all the trims
        reFormat = trim(reFormat) // SPACESUBSTITUTE
      case default
        reFormat = trim(reFormat) // toform(i:i)
      end select
    enddo
    if ( index(reFormat, SPACESUBSTITUTE) < 1 ) return
    ! Need to substitute a space for every occurrence of SPACESUBSTITUTE
    do
      i = index(reFormat, SPACESUBSTITUTE)
      if ( i < 1 ) return
      reFormat(i:i) = ' '
    enddo
  end function reFormatDate

  ! --------------------------------------------------  reFormatTime  -----
  function reFormatTime(time, form) result(reFormat)
    ! Reformat hhmmss.sss as hh:mm:ss
    ! (Note it truncates instead of rounding)
    ! Wouldn't it be clever to allow an optional string arg defining
    ! the output format; E.g. 'hh:mm' for '13:23'
    ! or 'HH:mm' for '01:23 PM'
    ! Args
    character(len=*), intent(in)            :: time
    character(len=len(time)+24)             :: reFormat
    character(len=*), optional, intent(in)  :: form
    ! Internal variables
    character(len=1), parameter             :: hmSpacer = ':'
    character(len=1), parameter             :: msSpacer = ':'
    integer                                 :: hours
    integer                                 :: i
    character(len=2)                        :: ampm
    character(len=2)                        :: hh
    ! Executable
    reFormat = time(1:2) // hmSpacer // time(3:4) // msSpacer // time(5:6)
    if ( .not. present(form) ) return
    if ( len_trim(form) < 1 ) return
    ampm = ' '
    reFormat = ' '
    i = 0
    do
      if ( i >= len_trim(form) ) exit
      i = i + 1
      select case (form(i:i))
      case ('h')
        reFormat = trim(reFormat) // time(1:2)
        i = i + 1
      case ('H')
        read(time(1:2), *) hours
        ampm = 'AM'
        if ( hours > 12 ) then
          hours = hours - 12
          ampm = 'PM'
        endif
        write(hh, '(i2.2)') hours
        reFormat = trim(reFormat) // hh
        i = i + 1
      case ('m')
        reFormat = trim(reFormat) // time(3:4)
        i = i + 1
      case ('s')
        reFormat = trim(reFormat) // time(5:6)
        i = i + 1
      case default
        reFormat = trim(reFormat) // form(i:i)
      end select
    enddo
    if ( ampm /= ' ' ) reFormat = trim(reFormat) // ' ' // ampm
  end function reFormatTime

   ! --------------------------------------------------  Reverse  -----
  FUNCTION Reverse (str) RESULT (outstr)
    ! takes a string and returns one with chars in reversed order
	 ! Useful in certain contexts:
	 ! e.g., to remove leading blanks
	 ! arg = Reverse(TRIM(Reverse(arg)))
	 !
	 ! See also ReverseList
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    INTEGER :: i, istr, irev
	 CHARACTER (LEN=1) :: strChar
    !----------Executable part----------!
    outstr=str
    IF(LEN(str) == 1) RETURN

    DO i = 1, LEN(str)-1, 2
       istr = 1 + (i-1)/2				! 1, 2, ..
       irev = LEN(str) - (i-1)/2		! N, N-1, ..
       strChar = str(istr:istr)
		 outstr(istr:istr) = str(irev:irev)
		 outstr(irev:irev) = strChar
    END DO

! Special case: str contains odd number of chars
    IF(MOD(LEN(str), 2) == 1) THEN
       istr = 1 + (LEN(str)-1)/2				! 1, 2, ..
        outstr(istr:istr) = str(istr:istr)
	ENDIF

  END FUNCTION Reverse

  ! -------------------------------------------------  SplitWords  -----

  ! This subroutine is based on my IDL one of the same name.
  ! A line of input is split into sets of words.  There are two ways in which
  ! this can be invoked.  Typically it is split into `first' and 'rest'
  ! However, if the threeway option is given it is split to first, rest and
  ! last.

  ! Note that there is a slight subtlety here, spaces are treated specially
  ! due to the use of TRIM.  Thus while two commas in a row would count as
  ! two separators, two spaces would count as one. Also if , is the separator
  ! then ,<space> counts as complete separator.
  
  ! Apologies--I have replaced almost all uses of the word "delimiter" with
  ! the word "separator" in a global manner, excepting only
  ! the optional last arg to this subroutine

  SUBROUTINE SplitWords(line,first,rest,last,&
       & threeWay,delimiter)

    ! Dummy arguments

    CHARACTER (LEN=*), INTENT(IN) :: line
    CHARACTER (LEN=*), INTENT(OUT) :: first
    CHARACTER (LEN=*), INTENT(OUT) :: rest
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: last

    LOGICAL, INTENT(IN), OPTIONAL :: threeWay
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: delimiter ! really separator

    ! Local variables

    CHARACTER (LEN=1) :: useseparator
    LOGICAL :: useThreeWay
    CHARACTER (LEN=LEN(line)) useLine ! Line with leading spaces removed

    INTEGER :: firstseparatorPos,lastseparatorPos,trimmedLen

    ! Executable code

    useLine=ADJUSTL(line)
    trimmedLen=LEN_TRIM(useLine)

    IF (PRESENT(delimiter)) THEN
       useseparator=delimiter
    ELSE
       useseparator=","
    END IF

    IF (PRESENT(threeWay)) THEN
       useThreeWay=threeWay 
    ELSE 
       useThreeWay=.FALSE.
    END IF

    ! Clear some results by default

    IF (PRESENT(last)) last=""
    rest=""

    ! Find the first separator

    firstseparatorPos=INDEX(useLine,useseparator)

    IF (firstseparatorPos == 0) THEN
       first=useLine
    ELSE
       first=useLine(1:firstseparatorPos-1)
       IF (useThreeWay) THEN
          ! In three way mode, find the last separator
          lastseparatorPos=INDEX(TRIM(useLine),useseparator,back=.TRUE.)
          IF (PRESENT(last) .AND. &
               & lastseparatorPos /= trimmedLen) THEN
             last=TRIM(useLine(lastseparatorPos+1:))
          END IF
          IF (firstseparatorPos+1 <= lastseparatorPos-1) THEN
             rest=TRIM(useLine(firstseparatorPos+1:lastseparatorPos-1))
          END IF
       ELSE
          IF (firstseparatorPos /= trimmedLen) THEN
             rest=TRIM(useLine(firstseparatorPos+1:))
          END IF
       END IF
    END IF

  END SUBROUTINE SplitWords
       
  ! --------------------------------------------------  strings2Ints  -----
  SUBROUTINE strings2Ints (strs, ints)
    ! takes an array of strings and returns integer array
    ! using "ichar"
	 ! Useful due to bug in toolbox swrfld
	 !
	 ! See also ints2Strings
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    CHARACTER (LEN=*), INTENT(in), dimension(:) ::   strs
    integer, intent(out), dimension(:,:) ::          ints

    !----------Local vars----------!
    INTEGER :: i, substr, strLen, arrSize
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   strLen = MIN(len(strs(1)), size(ints, dim=1))
   arrSize = MIN(size(strs), size(ints, dim=2))
   ints = LENORSIZETOOSMALL
   if(strLen <= 0 .or. arrSize <= 0) return
   ints=iachar(' ')
   do i=1, arrSize
      do substr=1, strLen
         ints(substr, i) = iachar(strs(i)(substr:substr))
      enddo
   enddo

  END SUBROUTINE strings2Ints

  ! --------------------------------------------------  writeIntsToChars  -----
  SUBROUTINE writeIntsToChars (ints, strs)
    ! takes an array of integers and returns string array
    ! using Fortran "write"
    ! If any element of string array is blank or contains one of forbiddens
    ! the corresponding element of ints is left undefined
	 ! Not useful yet
	 !
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    character (LEN=*), intent(out), dimension(:) ::   strs
    integer, intent(in), dimension(:)            ::   ints

    !----------Local vars----------!
    integer :: i, j, arrSize
    logical :: leave_undef
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   arrSize = MIN(size(strs), size(ints))
   if ( arrSize <= 0 ) then
     strs = ' '
     return
   endif
   do i=1, arrSize
      write(strs(i), *) ints(i)
   enddo

  END SUBROUTINE writeIntsToChars

  ! ---------------------------------------------  yyyymmdd_to_doy_ints  -----
  subroutine yyyymmdd_to_doy_ints(year, month, day, doy)
    ! Routine that returns the number of days after the year's start
    ! for year, month, day
    !--------Argument--------!
    integer, intent(in) :: year, month, day
    integer, intent(out) :: doy
    !----------Local vars----------!
    integer :: m
    integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     if ( year < 0 .or. year > YEARMAX ) then
       doy = -1
     endif
     doy = day
     if ( month <= 1 ) then
       return
     endif
     if ( leapyear(year) ) then
       DAYMAX = DAYMAXLY
     else
       DAYMAX = DAYMAXNY
     endif
     do m=1, month-1
       doy = doy + DAYMAX(m)
     enddo
     
  end subroutine yyyymmdd_to_doy_ints

  ! ---------------------------------------------  yyyymmdd_to_doy_str  -----
  subroutine yyyymmdd_to_doy_str(str, doy)
    ! Routine that returns the number of days after the year's start
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*),intent(in) :: str
    integer,intent(out) :: doy
    !----------Local vars----------!
    integer :: year, month, day
    integer :: ErrTyp
    !----------Executable part----------!
     ! call utc_to_yyyymmdd_ints(str, ErrTyp, year, month, day, nodash=.true.)
     doy = -1
     if ( len_trim(str) < 8 ) return 
     read(str(1:4), *) year
     read(str(5:6), *) month
     read(str(7:8), *) day
     call yyyymmdd_to_doy_ints(year, month, day, doy)
  end subroutine yyyymmdd_to_doy_str

  ! ---------------------------------------------  yyyydoy_to_yyyymmdd_str  -----
  subroutine yyyydoy_to_yyyymmdd_str(yyyy, doy, yyyymmdd)
    ! Routine that converts the string encoding the number of days 
    ! after the year's start, input as yyyydoy,
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*),intent(in) :: yyyy
    character(len=*),intent(in) :: doy
    character(len=*),intent(out) :: yyyymmdd
    !----------Local vars----------!
    integer :: year, month, day, doynum
    !----------Executable part----------!
     read(yyyy, *) year
     read(doy, *) doynum
     call yeardoy_to_yyyymmdd_ints(year, doynum, yyyymmdd)
  end subroutine yyyydoy_to_yyyymmdd_str

  ! ---------------------------------------------  yeardoy_to_yyyymmdd_ints  -----
  subroutine yeardoy_to_yyyymmdd_ints(year, doy, yyyymmdd)
    ! Routine that converts the string encoding the number of days 
    ! after the year's start, input as yyyydoy,
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    integer,intent(in)           :: year
    integer,intent(in)           :: doy
    character(len=*),intent(out) :: yyyymmdd
    !----------Local vars----------!
    integer :: day
    integer :: daysum
    integer :: m
    integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     if ( year < 0 .or. year > YEARMAX ) then
       yyyymmdd = 'year not in range'
     endif
     daysum = 0
     if ( leapyear(year) ) then
       DAYMAX = DAYMAXLY
     else
       DAYMAX = DAYMAXNY
     endif
     do m=1, 12
       if ( daysum + DAYMAX(m) >= doy ) exit
       daysum = daysum + DAYMAX(m)
     enddo
     if ( m > 12 ) then
       yyyymmdd = 'doy not in range'
       return
     endif
     day = doy - daysum
     write(yyyymmdd,'(I4.4, 2i2.2)') year, m, day
  end subroutine yeardoy_to_yyyymmdd_ints

  ! ---------------------------------------------------  Leapyear  -----
  logical function leapyear(year)
    integer,intent(in) :: year
     ! This is to capture the rule that centuries are leap only
     ! if divisible by 400
     ! Otherwise, as perhaps you knew, leapyears are those years divisible by 4
     leapyear = mod(year,4) == 0 .and. & ! Processor might short-circuit this
       & ( (mod(year,100) /= 0) .or. (mod(year,400) == 0) )
  end function leapyear

  ! ------------------------------------------  monthNameToNumber  -----
  function monthNameToNumber(name) result(number)
    ! Convert month name to corresponding number
    ! E.g., given 'March', returns 3
    ! As a courtesy, name may be case-insensitive
    ! As a further courtesy, name may be followed by any junk you like
    ! Thus 'March 23, 2004 01:59:59.999' still returns 3
    ! If no such month is found, returns -1
    ! Args
    character(len=*), intent(in)             :: name
    integer                                  :: number
    do number=1, size(MONTHNAME)
      if ( index(lowerCase(name), lowercase(trim(MONTHNAME(number)))) > 0 ) return
    enddo
    number = -1
  end function monthNameToNumber

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSStrings
!=============================================================================

! $Log$
! Revision 2.50  2004/10/13 16:23:03  pwagner
! Moved declaration of hhmmss_value result to after use statement
!
! Revision 2.49  2004/10/13 00:52:20  vsnyder
! Move HHMMSS_value here from MLSStringLists and simplify it.
! Remove AnyForbiddenChars and AllAllowedChars because intrinsic Scan and
! Verify can do them, and they weren't used anyway.
!
! Revision 2.48  2004/10/05 23:08:04  pwagner
! Added AnyForbiddenChars and AllApprovedChars; improved reFormatDate
!
! Revision 2.47  2004/09/23 22:56:38  pwagner
! Added reformats of date, time
!
! Revision 2.46  2004/09/16 00:15:52  pwagner
! Added writeIntsToChars
!
! Revision 2.45  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.44  2004/06/29 00:06:13  pwagner
! Tried to straighten out delimiter-separator terms in comments
!
! Revision 2.43  2004/06/16 22:15:28  pwagner
! Make lowerCase elemental
!
! Revision 2.42  2004/06/16 01:25:08  vsnyder
! Make Capitalize elemental
!
! Revision 2.41  2004/06/10 00:57:47  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.40  2004/06/09 00:02:35  pwagner
! GetUniqueList now accepts optional arg str2 returning str not in str2
!
! Revision 2.39  2004/01/27 21:34:02  pwagner
! Fixed some bugs in ExtractSubString
!
! Revision 2.38  2003/12/11 23:02:35  pwagner
! yyyymmdd_to_dai may take 3 ints or str
!
! Revision 2.37  2003/12/07 23:10:42  pwagner
! Added RemoveElemFromList; bug fixes in ReplaceSubString
!
! Revision 2.36  2003/12/05 00:52:18  pwagner
! Added yyyymmdd_to_dai (though arguably this belongs in time_m)
!
! Revision 2.35  2003/10/28 19:28:32  vsnyder
! Make sure outStr always has a value
!
! Revision 2.34  2003/10/15 00:34:19  pwagner
! Fixed the real bug in NumStringElements
!
! Revision 2.33  2003/10/14 18:17:02  pwagner
! Fixed problem with reducing switches to unique list
!
! Revision 2.32  2003/10/09 23:33:11  pwagner
! Added GetUniqueList
!
! Revision 2.31  2003/09/15 23:04:06  vsnyder
! Remove unused local variable
!
! Revision 2.30  2003/04/11 23:29:30  pwagner
! Fixed bug in ReplaceSubString; added ExtractSubString
!
! Revision 2.29  2003/02/27 18:36:57  pwagner
! utc_to_yyyymmdd optionally returns yy-mm-ddT00:00:00Z
!
! Revision 2.28  2003/02/19 19:08:40  pwagner
! Added ReplaceSubString
!
! Revision 2.27  2003/02/01 00:28:32  pwagner
! Added utc_to_yyyymmdd
!
! Revision 2.26  2003/01/15 21:20:00  pwagner
! Added readIntsFromChars
!
! Revision 2.25  2002/10/29 19:55:39  pwagner
! Fixed mistake in StringElementNum that caused crashes
!
! Revision 2.24  2002/10/29 01:00:05  pwagner
! optional param part_match added to str element routines
!
! Revision 2.23  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.22  2002/04/29 17:39:31  pwagner
! Comments re hhmmss_value mention mls_utctotai
!
! Revision 2.21  2002/02/22 23:35:42  pwagner
! SortList checks on lax elem length, not number of elems
!
! Revision 2.20  2002/02/22 01:19:57  pwagner
! SortArray not limited to array sizes lt MAXELEM
!
! Revision 2.19  2002/02/19 23:12:03  pwagner
! New optional args to Sorting routines
!
! Revision 2.18  2002/02/15 01:06:12  pwagner
! Added new array and sorting routines
!
! Revision 2.17  2002/01/09 23:46:05  pwagner
! Removed debugging stuff
!
! Revision 2.16  2001/08/03 00:03:08  pwagner
! Added ints2Strings and strings2Ints
!
! Revision 2.15  2001/06/20 23:23:39  vsnyder
! Same as last time, but for LowerCase instead of Capitalize.
!
! Revision 2.14  2001/06/20 23:21:49  vsnyder
! Replace ICHAR with IACHAR and CHAR with ACHAR, to improve portability.
! Make "offset" a parameter.
!
! Revision 2.13  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.12  2001/05/29 21:17:03  pwagner
! Removed Downcase; added table of contents
!
! Revision 2.11  2001/05/26 00:21:59  livesey
! Added downcase
!
! Revision 2.10  2001/05/24 23:36:17  pwagner
! Fixed problem with hhmmss_value
!
! Revision 2.9  2001/05/15 23:44:42  pwagner
! Added hhmmss_value
!
! Revision 2.8  2001/05/11 23:41:31  pwagner
! Improved unquote
!
! Revision 2.7  2001/05/11 00:06:54  pwagner
! Added unquote
!
! Revision 2.6  2001/03/14 17:34:00  pwagner
! Removed some of the dross, left all of the gold
!
! Revision 2.5  2001/03/02 19:33:14  pwagner
! Added GetIntHashElement
!
! Revision 2.4  2001/02/24 01:02:45  pwagner
! Added GetStringHashElement; alphabetized entries
!
! Revision 2.3  2001/02/23 00:05:56  pwagner
! Added 4 StringList functions
!
! Revision 2.2  2000/12/01 22:38:00  vsnyder
! Add lowercase function, alphabetize procedures, add "bookmarks"
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.21  2000/09/05 10:59:32  pumphrey
! HCP Fixed an out-of-bounds error
!
