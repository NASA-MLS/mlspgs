! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSStrings               ! Some low level string handling stuff
!=============================================================================

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
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
! CatStrings         Concatenate strings with a specified separator
! CompressString     Removes all leading and embedded blanks
! Count_words        Counts the number of space-separated words in a string
! Depunctuate        Replaces punctuation with blanks
! Hhmmss_value       Converts 'hh:mm:ss' formatted string to a real r8
!                    (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
! Indexes            Indexes an array of substrings of a string into an array
! Ints2Strings       Converts an array of integers to strings using "char" ftn
! IsRepeat           Is a string composed entirely of one substring repeated?
! LinearSearchStringArray     
!                    Finds string index of substring in array of strings
! LowerCase          tr[A-Z] -> [a-z]
! NCopies            How many copies of a substring in a string
! ReadCompleteLineWithoutComments     
!                    Knits continuations, snips comments
! ReadIntsFromChars  Converts an array of strings to ints using Fortran read
! Replace            Replaces every instance of oldChar with newChar
! Reverse            Turns 'a string' -> 'gnirts a'
! Reverse_trim       (Reverses after trimming its argument)
! SplitNest          Splits 'part 1 (part 2) part 3' -> 'part 1', 'part 2', 'part 3'
! SplitWords         Splits 'first, the, rest, last' -> 'first', 'the, rest', 'last'
! streq              Generalized strings "==" (optionally ignoring case,
!                      leading blanks, and allowing wildcard matches)
! Strings2Ints       Converts an array of strings to ints using "ichar" ftn
! trim_safe          trims string down, but never to length 0
! WriteIntsToChars   Converts an array of ints to strings using Fortran write
! === (end of toc) ===

! === (start of api) ===
! char* Capitalize (char* str)
! char* CompressString (char* str)
! int count_words (char* str)
! char* depunctuate (char* str)
! int(:) indexes (char* string, char* substrings, [char* mode])
! ints2Strings (int ints(:,:), char* strs(:))
! int LinearSearchStringArray (char* list(:), char* string, 
!   [log caseInsensitive, [log testSubstring], [log listInString])
! log IsRepeat ( char* str, [char* ssubtring] )
! char* LowerCase (char* str)
! int NCopies (char* str, char* substring, [log overlap])
! readIntsFromChars (char* strs(:), int ints(:), char* forbiddens)
! ReadCompleteLineWithoutComments (int unit, char* fullLine, [log eof], &
!       & [char commentChar], [char continuationChar])
! char* Replace (char* str, char oldChar, char newChar)
! char* Reverse (char* str)
! char* Reverse_trim (char* str)
! SplitNest ( char *str, char* part1, char* part2, char* part3, [char* parens] )
! SplitWords (char *line, char* first, char* rest, [char* last], &
!       & [log threeWay], [char* delimiter])
! log streq (char* str1, char* str2, [char* options])
! log streq (char* str1(:), char* str2, [char* options])
! log streq (char* str1, char* str2(:), [char* options])
! strings2Ints (char* strs(:), int ints(:,:))
! char* trim_safe (char* str)
! writeIntsToChars (int ints(:), char* strs(:))
! Many of these routines take optional arguments that greatly modify
! their default operation

! One standard is the character flag "options" which affects how loosely
! string matches may be interpreted
! it may include any of the following (poss. in combination, e.g. "-wc")
! w    Wildcard * which allows 'a*' to equal 'abcd'
! c    case insensitive which allows 'ABCD' to equal 'abcd'
! f    flush left which allows 'abcd' to equal '  abcd'
! (These are different in streq_array1, however--either redo that function
! to make it conform, or rename the options flag there to prevent
! unnecessary confusion)

! The above is to replace the countEmpty, caseSensitive, etc. that are
! separate optional args to many of the current module procedures

! Warnings: 
! (1) in the routine LinearSearchStringArray
! the input arguments include an array of strings;
! This array is of assumed-size
! I.e., all elements from array(1:size(array)) are relevant
! Therefore in calling one of these you probably need to use the format
!   call SortArray(myArray(1:mySize), ..
! to avoid operating on undefined array elements
! (2) in some routines trailing blanks are ignored,
!   while in others they are significant
! How trailing blanks in argument(s) treated
! Ignored                   Significant
! -------                   -----------
! CatStrings                Capitalize
! Depunctuate               CompressString
! HHMMSS_value              count_words
! indexes                   LowerCase
! isRepeat                  Reverse
! LinearSearchStringArray   strings2Ints
! ncopies                   trim_safe
! readIntsFromChars         
! reFormatDate              
! reFormatTime              
! Reverse_trim              
! SplitWords                
! streq                     
! === (end of api) ===

  public :: Capitalize, CatStrings, CompressString, count_words, &
   & depunctuate, hhmmss_value, &
   & indexes, ints2Strings, IsRepeat, &
   & LinearSearchStringArray, LowerCase, NCopies, &
   & ReadCompleteLineWithoutComments, readIntsFromChars, &
   & Replace, Reverse, Reverse_trim, &
   & SplitNest, SplitWords, streq, strings2Ints, trim_safe, &
   & writeIntsToChars

  interface readIntsFromChars
    module procedure readAnIntFromChars, readIntArrayFromChars
  end interface

  interface streq
    module procedure streq_scalar, streq_array1, streq_array2
  end interface

  interface writeIntsToChars
    module procedure writeAnIntToChars, writeIntArrayToChars
  end interface

  ! hhmmss_value
  integer, public, parameter :: INVALIDHHMMSSSTRING = 1
  ! strings2Ints
  integer, public, parameter :: LENORSIZETOOSMALL=-999
  
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

  ! -------------------------------------------------  CatStrings  -----
  subroutine CatStrings ( Strings, Sep, StringsCat, L )
  ! Concatenate Strings with Sep between them, giving StringsCat(:L-1)
    character(len=*), intent(in) :: Strings(:)
    character(len=*), intent(in) :: Sep
    character(len=*), intent(out) :: StringsCat
    integer, intent(out) :: L
    integer :: I, N, T, W
    w = len(sep)
    l = len_trim(strings(1)) + 1
    stringsCat(:l-1) = strings(1)(:l-1)
    do i = 2, size(strings)
      t = len_trim(strings(i))
      n = l + t + w
      stringsCat(l:n-1) = sep // strings(i)(:t)
      l = n
    end do
  end subroutine CatStrings

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

  ! ---------------------------------------------------  indexes  -----
  function indexes(str, substrings, mode) result(array)
    character(len=*), intent(in) :: str
    character(len=*), intent(in), dimension(:) :: substrings
    character(len=*), optional :: mode
    integer, dimension(size(substrings)) :: array
    character(len=len(str)) :: substr
    ! Returns the array of indexes of each element of substrings in str
    ! The mode determines how consecutive elements of array are ordered
    ! (default is first)
    ! mode          order
    ! ------        -----
    ! first         always find first occurrence of substr(i) in str
    ! last          always find last occurrence of substr(i) in str
    ! left          progressively find next substr(i) after substr(i-1)
    ! right         progressively find next substr(i) before substr(i-1)
    ! wrap          left-right, meeting in middle
    !
    ! E.g., if str='ababababababa' and substrings = (/'a', 'a', 'a', 'a', 'a'/)
    ! mode            result
    ! ------          ------
    ! first        (/1, 1, 1, 1, 1/)
    ! last         (/13, 13, 13, 13, 13/)
    ! left         (/1, 3, 5, 7, 9/)
    ! right        (/13, 11, 9, 7, 5/)
    ! wrap         (/1, 3, 7, 11, 13/) ! Currently, we give (/1,3,9,11,13/)
    !
    ! Notes:
    ! mode='wrap' does not return true middle value yet--do we care?
     integer :: i
     integer :: ipos
     integer :: lpos
     integer :: n
     integer :: rpos
     integer, dimension(size(substrings)) :: left
     integer, dimension(size(substrings)) :: right
     character(len=5) :: myMode
     !
     myMode = 'first'
     if ( present(mode) ) myMode = mode
     n = size(substrings)
     ! Simple modes
     if ( lowercase(myMode(1:2)) == 'fi' ) then
       do i = 1, n
         array(i) = index(str, trim(substrings(i)))
       enddo
       return
     elseif ( lowercase(myMode(1:2)) == 'la' ) then
       do i = 1, n
         array(i) = index(str, trim(substrings(i)), back=.true.)
       enddo
       return
     endif
     ! Progressive Modes (which accumulate *pos)
     lpos = 1
     rpos = max(len_trim(str), 1) ! len(trim_safe(str))
     left = 0
     do i = 1, n
       ipos = index(str(lpos:), trim_safe(substrings(i)))
       if ( ipos < 1 ) exit
       left(i) = lpos + ipos - 1
       lpos = left(i) + max(len_trim(substrings(i)), 1) ! len(trim_safe(substr(i)))
       if ( lpos > len(str) ) exit
     enddo
     right = 0
     do i = 1, n
       ipos = index(str(:rpos), trim_safe(substrings(i)), back=.true.)
       if ( ipos < 1 ) exit
       right(i) = ipos
       rpos = right(i) - 1
       if ( rpos < 1 ) exit
     enddo
     select case (lowercase(myMode(1:2)))
     case ('le')
       array = left
     case ('ri')
       array = right
     case ('wr')
       array(1:n) = right(n:1:-1)
       lpos = n/2
       rpos = lpos + 1
       do i = 1, lpos
         array(i) = left(i)
       enddo
       ! Is n an odd number?
       !if ( 2*lpos < n ) &
       !  & array(rpos) = (left(rpos) + right(rpos) ) / 2 ! No, this won't work
     case default
       array = left
     end select
  end function indexes

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

  ! ---------------------------------------------------  isRepeat  -----
  logical function isRepeat(str, substr)
    character(len=*), intent(in) :: str
    character(len=*), intent(in), optional :: substr
     ! Is str formed purely of repeated blocks of substr?
     ! If substr not supplied, 
     ! then is str any one character repeated over and over?
     ! Special cases:
     ! str = ' ' => TRUE 
     ! str != ' ' and substr = ' ' => FALSE
     ! Note that otherwise we're ignoring trailing blanks
     integer :: strlen
     integer :: substrlen
     character(len=1) :: aChar
     isRepeat = .false.
     if ( len_trim(str) < 1 ) then
       isRepeat = .true.
       return
     endif
     strlen = len_trim(str)
     if ( present(substr) ) then
       substrlen = len_trim(substr)
       if ( substrlen < 1 ) return
       isRepeat = (substrlen*ncopies(trim(str), trim(substr)) >= strlen)
       return
     endif
     aChar = str(1:1)
     isRepeat = (ncopies(trim(str), aChar) >= strlen)
  end function isRepeat

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

  ! ---------------------------------------------------  ncopies  -----
  integer function ncopies(str, substr, overlap)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: substr
    logical, optional, intent(in) :: overlap
     ! How many copies of substr are in str?
     ! The copies may optionally overlap
     ! E.g., str = 'aaaa', substr = 'aa'
     ! overlap = FALSE => ncopies = 2
     ! overlap = TRUE => ncopies = 3
     integer :: ipos
     integer :: next
     logical :: mayOverlap
     !
     mayOverlap = .false.
     if ( present(overlap) ) mayOverlap = overlap
     ncopies = 0
     ipos = 1
     do
       if ( ipos > len_trim(str) ) return
       next = index(str(ipos:), substr)
       if ( next < 1 ) then
         return
       endif
       ncopies = ncopies + 1
       if ( .not. mayOverlap ) then
         ipos = ipos + next + len(substr) - 1
       else
         ipos = ipos + next
       endif
     enddo
  end function ncopies

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

  ! --------------------------------------------------  readAnIntFromChars  -----
  SUBROUTINE readAnIntFromChars (str, int, forbiddens, ignore)
    ! takes a string and returns an integer
    ! using Fortran "read"
    ! (which could cause an io error--that's why this
    ! subroutine exists, to filter out invalid characters)
    ! If the string is blank or contains one of forbiddens
    ! the int is left undefined
    
    ! Then snip away any from the set ignore if present
    ! If ignore is '*', that means ignore all alphabetical chars
    ! If ignore contains '*', that means ignore all alphabetical chars
    ! plus any other chars among ignore
    ! If the string is composed entirely of ignorable chars, int is 0
    
    ! Finally attempt to read as an int what remains

    ! Examples:
    ! (1) if str='band13a' and ignore='*', int will be 13
    ! (2) if str='3 cm' and forbiddens='c', int will be left undefined

    ! Limitation: you're unable to "escape" a * so you'll have to
    ! preprocess the * away if you really want to read a string which has
    ! a * in it somewhere
	 !
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(in) ::   str
    integer, intent(out)          ::   int
    CHARACTER (LEN=*), INTENT(in), optional     ::   forbiddens
    CHARACTER (LEN=*), INTENT(in), optional     ::   ignore

    !----------Local vars----------!
    INTEGER :: j, k
    LOGICAL :: leave_undef
    character(len=40)                           ::   myForbiddens
    character(len=40)                           ::   myIgnore
    character(len=len(str))                     ::   myStr
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   if ( present(forbiddens) ) then
     myForbiddens = adjustl(forbiddens)
   else
     myForbiddens = ' '
   endif
   if ( present(ignore) ) then
     myignore = adjustl(ignore)
   else
     myignore = ' '
   endif
   leave_undef = (str == ' ')
   if ( myForbiddens /= ' ' ) then
     do j=1, len(trim(myForbiddens))
        leave_undef = leave_undef &
         & .or. &
         & ( &
         &    index(str, myForbiddens(j:j)) > 0 &
         &  .and. &
         &    myForbiddens(j:j) /= ' ' &
         & )
     enddo
   endif
   if ( leave_undef ) then
     return
   elseif (  myIgnore == "" ) then
     read(str, *) int
   ! elseif (  myIgnore == "*" ) then
   elseif (  index(myIgnore, "*") /= 0 ) then
     int = 0  ! a str made up entirely of ignorables means "0"
     k = 1
     myStr = ""
     do j = 1, len(str)
       if ( .not. isAlphabet(str(j:j)) .and. &
         & index(myIgnore, str(j:j)) < 1 ) then
         myStr(k:k) = str(j:j)
         k = k + 1
       endif
     enddo
     if ( myStr /= "" ) read(mystr, *) int
   else
     int = 0  ! a str made up entirely of ignorables means "0"
     k = 1
     myStr = ""
     do j = 1, len(str)
       if ( index(myIgnore, str(j:j)) < 1 ) then
         myStr(k:k) = str(j:j)
         k = k + 1
       endif
     enddo
     if ( myStr /= "" ) read(mystr, *) int
   endif

  END SUBROUTINE readAnIntFromChars

  ! --------------------------------------------------  readIntArrayFromChars  -----
  SUBROUTINE readIntArrayFromChars (strs, ints, forbiddens)
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
   do i=1, arrSize
     call readAnIntFromChars(strs(i), ints(i), forbiddens)
   enddo

  END SUBROUTINE readIntArrayFromChars

   ! --------------------------------------------------  Replace  -----
  function Replace (str, oldChar, newchar) RESULT (outstr)
    ! takes a string and returns one with oldChar replaced by newChar
    ! E.g., to replace every char(0), which is the NUL character, with a blank
	 ! arg = Replace( arg, char(0), char(32) )
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: oldChar
    character(len=1), intent(in) :: newChar
    character(len=len(str))      :: outstr
    ! Internal variables
    integer :: i, n
    ! Executable
    outstr = str
    if ( len(str) < 1 ) return
    if ( index(str, oldChar) < 1 ) return
    n = len(str)
    do i=1, n
      if ( str(i:i) == oldChar ) outstr(i:i) = newChar
    enddo
  end function Replace

	 !
   ! --------------------------------------------------  Reverse  -----
  elemental function Reverse (str) RESULT (outstr)
    ! takes a string and returns one with chars in reversed order
	 ! Useful in certain contexts:
	 ! e.g., to remove leading blanks
	 ! arg = Reverse(TRIM(Reverse(arg)))
	 !
	 ! See also Reverse_trim, ReverseList
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

  end function Reverse

   ! --------------------------------------------------  Reverse_trim  -----
  function Reverse_trim (str) RESULT (outstr)
    ! takes a string, trims it then returns one with chars in reversed order
	 ! See also Reverse which omits the trim step
    !
    ! E.g., given 'A string    ' reverse_trim returns 'gnirst A   ' while
    ! a simple Reverse returns '   gnirst A'
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=max(LEN_TRIM(str), 1)) :: outstr

    !----------Local vars----------!
    INTEGER :: i, istr, irev
	 CHARACTER (LEN=1) :: strChar
    !----------Executable part----------!
    outstr=str
    IF(LEN_TRIM(str) <= 1) RETURN

    DO i = 1, LEN_TRIM(str)-1, 2
       istr = 1 + (i-1)/2				! 1, 2, ..
       irev = LEN_TRIM(str) - (i-1)/2		! N, N-1, ..
       strChar = str(istr:istr)
		 outstr(istr:istr) = str(irev:irev)
		 outstr(irev:irev) = strChar
    END DO

! Special case: str contains odd number of chars
    IF(MOD(LEN_TRIM(str), 2) == 1) THEN
       istr = 1 + (LEN_TRIM(str)-1)/2				! 1, 2, ..
        outstr(istr:istr) = str(istr:istr)
	ENDIF

  end function Reverse_trim

  ! ---------------------------------------------  SplitNest  -----

  ! This subroutine splits an input string into 3 parts; it is most
  ! easily described with the following diagram:
  ! Given "part 1 (part 2) part 3"
  ! it returns "part 1" "part 2" "part 3"
  ! Note the paramount role played by the nesting parentheses
  ! If given "No parentheses here", part 2 and part 3 will be empty
  ! If given "(part 2)" part 1 and part 3 will be empty
  ! The cases where only part 1 or only part 3 would be empty are obvious
  ! To permit recursive passes, the split is made at the rightmost pair
  ! of nesting parentheses
  ! Thus given
  ! "(((a))) (b (c))"
  ! it will return "(((a))) (b " "c" and ")"
  ! The motivation is to recursively parse expressions such as
  ! "p or not (q and (r or s))"
  ! by turning this by successive steps of collapsing the nesting parentheses
  ! "p or not (q and t)"
  ! "p or not u"
  ! Since we collapse one nesting level each time, we are assured of
  ! arriving at an expression without nests in a finite number of steps

  subroutine SplitNest( str, part1, part2, part3, parens )
    character(len=*), intent(in)           :: str
    character(len=*), intent(out)          :: part1
    character(len=*), intent(out)          :: part2
    character(len=*), intent(out)          :: part3
    character(len=*), intent(in), optional :: parens
    ! Local variables
    character(len=1) :: closeParen ! Usu. right parenthesis
    integer          :: firstOpen
    integer          :: matchingClose
    character(len=1) :: openParen  ! Usu. left parenthesis
    ! Executable
    if ( present(parens) ) then
      openParen  = parens(1:1)
      closeParen = openParen
      if ( len_trim(parens) > 1 ) closeParen = parens(2:2)
    else
      openParen  = '('
      closeParen = ')'
    endif
    part1 = ' '
    part2 = ' '
    part3 = ' '
    if ( str == ' ' ) return
    firstOpen = index( str, openParen, back=.true. )
    if ( firstOpen < 1 ) then
      part1 = str
      return
    endif
    matchingClose = index( str(firstOpen+1:), closeParen )
    if ( matchingClose < 1 ) then
      ! This probably means an ill-formed expression with an unmatched '('
      part1 = str
      return
    endif
    matchingClose = matchingClose + firstOpen ! Referenced back to original
    if ( firstOpen > 1 ) part1 = str(1:firstOpen-1)
    if ( matchingClose < len_trim(str) ) part3 = str(matchingClose+1:)
    part2 = str(firstOpen+1:matchingClose-1)
  end subroutine SplitNest

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
       
  ! -------------------------------------------------  streq_array1  -----
  function streq_array1 (STR1, STR2, OPTIONS) result (relation)
    ! Array version of streq
    ! May return multiple TRUEs (except see 's', 'l' options)
    ! Extra options
    ! 'p' is partial match (STR2 is replaced by '*' // STR2 // '*'
    ! 's' returns TRUE only in element corresponding to shortest STR1
    ! 'l' returns TRUE only in element corresponding to longest STR1
    use MLSSets, only: findFirst
    character (len=*), dimension(:), intent(in) :: STR1
    character (len=*), intent(in)               :: STR2
    character (len=*), intent(in), optional     :: OPTIONS
    logical, dimension(size(STR1))              :: RELATION
    ! Internal variables
    integer :: candidate
    integer :: candidateLength
    integer :: i
    integer, dimension(size(STR1))              :: lengths
    character(len=8) :: myOptions
    character(len=len(str2)+2) :: mystr2
    ! Executable
    myOptions = ''
    if ( present(options) ) myOptions = lowercase(options)
    mystr2 = str2
    if ( index(myOptions, 'p') > 0 ) mystr2 = '*' // trim_safe(str2) // '*'
    do i=1, size(str1)
      relation(i) = streq(str1(i), mystr2, OPTIONS)
    enddo
    if ( .not. present(options) .or. count(relation) < 2 ) return
    lengths = len_trim(str1)
    if ( index(myOptions, 's') > 0 ) then
      candidate = findFirst(relation)
      candidateLength = lengths(candidate)
      do i=1, size(str1)
        if ( relation(i) .and. lengths(i) < candidateLength ) then
          candidate = i
          candidateLength = lengths(candidate)
        endif
      enddo
      relation = .false.
      relation(candidate) = .true.
    elseif ( index(myOptions, 'l') > 0 ) then
      candidate = findFirst(relation)
      candidateLength = lengths(candidate)
      do i=1, size(str1)
        if ( relation(i) .and. lengths(i) > candidateLength ) then
          candidate = i
          candidateLength = lengths(candidate)
        endif
      enddo
      relation = .false.
      relation(candidate) = .true.
    endif
  end function streq_array1

  ! -------------------------------------------------  streq_array2  -----
  function streq_array2 (STR1, STR2, OPTIONS) result (relation)
    ! Array version of streq
    ! (see streq_array1)
    ! Here str2 is array, not str1
    character (len=*), dimension(:), intent(in) :: STR2
    character (len=*), intent(in)               :: STR1
    character (len=*), intent(in), optional     :: OPTIONS
    logical, dimension(size(STR2))              :: RELATION
    relation = streq_array1(str2, str1, options)
  end function streq_array2

  ! -------------------------------------------------  streq_scalar  -----
  function streq_scalar (STR1, STR2, OPTIONS) result (relation)
    ! Are two strings "equal" where equality is modified by
    ! (1) Wildcard * (off by default) which allows 'a*' to equal 'abcd'
    ! (2) case insensitive (off by default) which allows 'ABCD' to equal 'abcd'
    ! (3) flush left (off by default) which allows 'abcd' to equal '  abcd'
    !
    ! Defaults to standard (str1 == str2), but options broaden cases of TRUE
    ! To turn options on, supply optional arg options which is a character
    ! string containing: 'w' => turns on (1); 'c' => turns on (2); 'f' => (3)
    ! e.g., streq('Ab*', 'abcd', '-wc') is TRUE
    !
    ! Notes:
    ! The '-' character in options is ignored and therefore not necessary
    ! A more powerful version can be imagined that would permit full regexp
    ! Trailing spaces are always ignored; e.g. streq('abcd ', 'abcd') is TRUE
    ! Only one of str1, str2 may contain wildcards
    !--------Argument--------!
    character (len=*), intent(in) :: STR1
    character (len=*), intent(in) :: STR2
    character (len=*), intent(in), optional  :: OPTIONS
    logical                       :: RELATION

    ! Internal variables
    integer, parameter :: MAXNUMWILDCARDS = 10 ! How many '*' in the pattern
    character(len=*), parameter :: star = '*'  ! Should we allow others?
    logical :: flushleft
    integer :: i
    logical :: ignorecase
    integer, dimension(MAXNUMWILDCARDS) :: istars
    character(len=8) :: myOptions
    integer :: nstars
    integer :: spos
    character(len=max(len(str1), len(str2))) :: str
    character(len=len(str)) :: ptrn  ! The one with '*'
    character(len=len(str)), dimension(MAXNUMWILDCARDS+1) :: substrs
    logical :: wildcard
    logical, parameter :: deebug = .false.
    !----------Executable part----------!
    relation = .FALSE.
    myOptions = ' '
    if ( present(options) ) myOptions = lowercase(options)
    wildcard = (index(myoptions, 'w') > 0)
    ignorecase = (index(myoptions, 'c') > 0)
    flushleft  = (index(myoptions, 'f') > 0)
    ! Now the wildcard(s) should be in ptrn
    if ( index(str1, star) > 0 ) then
      str = str2
      ptrn = str1
    elseif ( index(str2, star) > 0 ) then
      str = str1
      ptrn = str2
    else
      ! Whoops--don't need wild card after all
      wildcard = .false.
    endif
    ! Special cases: emptry strings or a bare vanilla match
    if ( len_trim(str1 // str2) < 1 ) then
      relation = .true.
      return
    elseif ( len_trim(str1) < 1 ) then
      relation = ( str2 == star .and. wildcard ) 
      return
    elseif ( len_trim(str2) < 1 ) then
      relation = ( str1 == star .and. wildcard ) 
      return
    elseif ( str2 == str1 ) then
      relation = .true.
      return
    endif

    if ( .not. wildcard ) then
      if ( ignorecase ) then
        if ( flushleft ) then
          relation = (adjustl(lowercase(str1)) == adjustl(lowercase(str1)))
        else
          relation = (lowercase(str1) == lowercase(str1))
        endif
      else
        if ( flushleft ) then
          relation = (adjustl(str1) == adjustl(str1))
        ! else plain vanilla case already handled as special case above
        endif    
      endif    
      return
    endif    

    ! How to handle the wildcard? Use it to split ptrn into sub-patterns
    if ( ignorecase ) then
      str = lowercase(str)
      ptrn = lowercase(ptrn)
    endif
    if ( flushleft ) then
      str = adjustl(str)
      ptrn = adjustl(ptrn)
    endif
    ! if ( deebug ) print *, 'str: ', trim(str), '  ptrn: ', trim(ptrn)
    if ( ptrn == star ) then
      relation = .true.
      return
    elseif ( isRepeat(ptrn, star) ) then
      relation = .true.
      return
    endif

    ! 1st -- how many stars?
    nstars =     ncopies(ptrn, star)
    ! if ( deebug ) print *, 'num of * ', nstars
    ! 2nd -- extract substrings from inbetween the wildcards
    substrs = star
    istars(1:nstars) = indexes(ptrn, substrs(1:nstars), mode='left')
    ! if ( deebug ) print *, 'where? ', istars(1:nstars)
    substrs = ' '
    spos = 1
    do i=1, nstars
      if ( spos > len(ptrn) ) exit
      substrs(i) = firstsubstr(ptrn(spos:), star)
      spos = max(spos, istars(i)) + 1
    enddo
    substrs(nstars+1) = Reverse_trim(firstsubstr(Reverse_trim(ptrn), star))
!     if ( deebug ) then
!       do i=1, nstars+1
!         print *, trim(substrs(i))
!       enddo
!     endif
    ! Deal specifically with empty elements of substrs
    relation = .true.
    if ( substrs(1) /= ' ' ) then
      if ( index(str, trim(substrs(1))) /= 1 ) then
        relation = .false.
        return
      endif
    endif
    ! if ( deebug ) print *, 'passed 1st sub-test'
    ! firstSSindex = 2
    if ( substrs(nstars+1) /= ' ' ) then
      if ( index(Reverse_trim(str), Reverse_trim(substrs(nstars+1))) /= 1 ) then
        relation = .false.
        ! if ( deebug ) print *, 'failed 2nd sub-test: ', Reverse_trim(str), Reverse_trim(substrs(nstars+1))
        return
      endif
    endif
    ! lastSSindex = nstars
    ! if ( deebug ) print *, 'passed 2nd sub-test'
    if ( nstars < 2 ) return
    ! Now find the indexes of these sub-patterns according to mode='wrap'
    istars(1:nstars-1) = indexes(str, substrs(2:nstars), mode='wrap')
    ! What we want to do is to check that no istars < 1
    relation = all ( istars(1:nstars-1) > 0 )
  end function streq_scalar

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

  ! -------------------------------------------------  TRIM_SAFE  -----
  function trim_safe (STR) result (OUTSTR)
    ! trims str returning a string of length no less than 1
    ! similar to trim, but will return a single blank character
    ! Useful in those cases where trim would result in strings of length 0
    ! E.g., MakeHDFAttribute(trim(' ')) fails but
    ! E.g., MakeHDFAttribute(trim_safe(' ')) succeeds
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=max(len_trim(str), 1)) :: OUTSTR

    !----------Executable part----------!
    outstr=' '

    if ( len_trim(str) > 0 ) outstr=trim(str)

  end function trim_safe

  ! --------------------------------------------------  writeAnIntToChars  -----
  SUBROUTINE writeAnIntToChars (int, str, fmt, specialInts, specialChars)
    ! takes an integer and returns a string
    ! using Fortran "write"
    ! Unless integer is one of specialInts, in which case
    ! we return corresponding one of specialChars
    ! (So that we can treat -1 as "unlimited' or -999 as 'FillValue')
	 ! We'll just assume both special arrays are of same size
	 ! Not useful yet
	 !
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    integer, intent(in)                                    ::   int
    character (LEN=*), intent(out)                         ::   str
    character (LEN=*), optional, intent(in)                ::   fmt
    integer, intent(in), dimension(:), optional            ::   specialInts
    character (LEN=*), intent(in), dimension(:), optional  ::   specialChars

    !----------Local vars----------!
    integer :: i
    !----------Executable part----------!

   ! Check that we don't have one of special cases
   if ( present(specialInts) ) then
     do i=1, size(specialInts)
       if ( int == specialInts(i) ) then
         str = specialChars(i)
         return
       endif
     enddo
   endif
   if ( present(fmt) ) then
     write(str, fmt=fmt) int
   else
     write(str, *) int
   endif

  END SUBROUTINE writeAnIntToChars

  ! --------------------------------------------------  writeIntArrayToChars  -----
  SUBROUTINE writeIntArrayToChars (ints, strs, fmt, specialInts, specialChars)
    ! takes an array of integers and returns string array
    ! using Fortran "write"
	 ! Not useful yet
	 !
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    integer, intent(in), dimension(:)            ::   ints
    character (LEN=*), intent(out), dimension(:) ::   strs
    character (LEN=*), intent(in), optional      ::   fmt
    integer, intent(in), dimension(:), optional  ::   specialInts
    character (LEN=*), intent(in), dimension(:), optional  ::   specialChars

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
      call writeAnIntToChars(ints(i), strs(i), fmt, specialInts, specialChars)
   enddo

  END SUBROUTINE writeIntArrayToChars

  ! Private procedures and functions
  !
  ! ---------------------------------------------------  isAlphabet  -----
  function isAlphabet(arg, inputcase) result(itIs)
    ! Returns TRUE if arg alphabetical; 
    ! i.e.is one of {'a', 'b', ..}
    ! Note: to check if input is UPPER  lower, either, set
    ! inputcase          case
    ! ----------         ----
    !   UPPER             u
    !   lower             l
    !   either            e (default)
    ! Args
    character(len=1), intent(in) :: arg
    character(len=1), optional, intent(in) :: inputcase
    logical                      :: itIs
    ! Internal variables
    logical :: itsEither
    logical :: itsLower
    character(len=*), parameter :: list='abcdefghijklmnopqrstuvwxyz'
    character(len=1)            :: myCase
    ! Executable
    myCase = 'e'
    if ( present(inputcase) ) myCase = inputcase
    itsLower = ( index(list, arg) > 0 )
    itsEither = ( index(list, lowercase(arg)) > 0 )
    select case(myCase)
    case ('u')
      itIs = itsEither .and. .not. itsLower
    case ('l')
      itIs = itsLower
    case default
      itis = itsEither
    end select    
  end function isAlphabet

  ! ---------------------------------------------------  isDigit  -----
  function isDigit(arg) result(itIs)
    ! Returns TRUE if arg is one of {'1', '2', ..}
    ! Args
    character(len=1), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    character(len=*), parameter :: list='1234567890'
    ! Executable
    itIs = ( index(list, arg) > 0 )
  end function isDigit

  ! ---------------------------------------------------  lastchar  -----
  character function lastchar(str)
    character(len=*), intent(in) :: str
     ! Returns the last non-blank character of str (unless str itself is blank)
     integer :: strlen
     lastchar = ' '
     strlen = len_trim(str)
     if ( strlen < 1 ) return
     lastchar = str(strlen:strlen)
  end function lastchar

  ! ---------------------------------------------------  firstsubstr  -----
  elemental function firstsubstr(str, star) result(substr)
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: star
    character(len=len(str)) :: substr
    ! Returns the substr between the start of str and the 1st occurrence of star
     integer :: strlen
     integer :: ipos
     substr = ' '
     strlen = len_trim(str)
     if ( strlen < 1 ) return
     ipos = index(str, star)
     if ( ipos < 1 ) then
       substr = str
     elseif ( ipos == 1 ) then
       substr = ' '
     else
       substr = str(:ipos-1)
     endif
  end function firstsubstr

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSStrings
!=============================================================================

! $Log$
! Revision 2.62  2006/05/09 00:14:23  pwagner
! Added Replace; useful to replace null chars with blanks
!
! Revision 2.61  2006/02/24 01:14:02  pwagner
! Added SplitNest (is this the best name)
!
! Revision 2.60  2006/02/16 00:58:12  pwagner
! ignore arg to ReadIntsFromChars may include * plus others
!
! Revision 2.59  2005/09/22 23:33:58  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.58  2005/07/21 23:37:22  pwagner
! Added explanation of to-be-standard character flag options
!
! Revision 2.57  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.56  2005/06/14 18:31:41  pwagner
! readIntsFromChars can now ignore certain chars
!
! Revision 2.55  2005/05/31 17:46:01  pwagner
! Added array  interfaces for streq
!
! Revision 2.54  2005/04/12 17:34:53  pwagner
! isRepeat, streq, ncopies, reverse_trim, scalar versions of read/write intchar
!
! Revision 2.53  2005/03/15 23:45:05  pwagner
! Added trim_safe to stop trimming at length 1
!
! Revision 2.52  2005/01/20 01:29:42  vsnyder
! Add CatStrings
!
! Revision 2.51  2004/10/13 20:25:45  pwagner
! In reFormatDate allow day of year w/o letter d; e.g. 2004-274
!
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
