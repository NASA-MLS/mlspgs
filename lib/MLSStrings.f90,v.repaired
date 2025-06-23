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
module MLSStrings               ! Some low level string handling stuff
!=============================================================================
  use MLSFinds, only: FindFirst, FindNext, FindLast
  use MLSStrings_0, only: Asciify, Capitalize, DecimalCode, &
    & IsAllAscii, IsAlphabet, &
    & NCharsInFormat, ReplaceNonAscii, Lowercase, &
    & MNemonicCode, OctalCode, &
    & ReadIntsFromChars, Stretch, Trim_Safe
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

!     (parameters and data)
! StringOptions      Default options string

!     (subroutines and functions)
! Asciify            purify chars to be within printing range [32,126]
!                      (no binary) (see also ReplaceNonAscii, unAsciify)
! Capitalize         tr[a-z] -> [A-Z]
! CatStrings         Concatenate strings with a specified separator
! CharToInt          Convert a single char to its corresponding integer
!                      tr['1'-'9'] -> [1-9]
! ShiftLRC           Shift string to left, center, or right
! CompressString     Removes all leading and embedded blanks
!                      (should rename to 'unstretch'?)
! Count_quotes       Counts the number of quotes-surrounded phrases in a string
! Count_words        Counts the number of space-separated words in a string
! Delete             Deletes each instance of a char
! Depunctuate        Replaces punctuation with blanks
! Enclosure          Finds outermost indices of substring enclosing an arg
!                      with '(' and ')'; optionally some other brackets
! FlushArrayLeft     Flush character array left over blank elements
! Hhmmss_value       Converts 'hh:mm:ss' formatted string to a real r8
!                    (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
! Indexes            Indexes an array of substrings of a string into an array
! Ints2Strings       Converts an array of integers to strings using "char" ftn
! IsAllAscii         Is a string composed entirely of ascii, i.e. non-binary
! isAlphabet         Is the arg an alphabetical character?
! IsComment          Does the string begin with a Comment character?
! IsDigits           Is the string composed entirely of digits?
! IsRepeat           Is a string composed entirely of one substring repeated?
! Justify            Return the string str Justified to fit the available length
! LenTrimToAscii     len_trim of string ignoring all non-ascii chars
! LinearSearchStringArray     
!                    Finds string index of substring in array of strings
! LowerCase          tr[A-Z] -> [a-z]
! NAppearances       The number of times each substring appears in string
! NCharsInFormat     How many characters in a format spec
! NCopies            How many copies of a substring in a string
! ReadCompleteLineWithoutComments     
!                    Knits continuations, snips comments
! ReadNumFromBaseN   Interprets a string as base n representation; returns float
! ReadIntFromBaseN   Interprets a string as base n representation
! ReadIntsFromChars  Converts an [array of] strings to int[s] using Fortran read
! ReadNumsFromChars  Converts an [array of] strings to num[s] using Fortran read
! ReadRomanNumerals  Converts a Roman numeral (e.g. 'ix') to its integer value
! Remap              Replaces every instance in str of char old[k] with 
!                       corresponding new[k]
! Replace            Replaces every instance of oldChar with newChar
! ReplaceNonAscii    Replaces every non-ascii char with newChar (see also Asciify)
! Reverse            Turns 'a string' -> 'gnirts a'
! Reverse_trim       (Reverses after first trimming its argument)
! Rot13              Like ROT13 but for general integer nn
! Size_trim          Returns len_trim of equivalent character scalar for array
! SplitDetails       Splits 'pro1' into 'pro' and 1
! SplitNest          Splits 'part 1 (part 2) part 3' -> 'part 1', 'part 2', 'part 3'
! SplitWords         Splits 'first, the, rest, last' -> 'first', 'the, rest', 'last'
! Squeeze            Snip excess spaces from between words; optionally snip all
! StartCase          Capitalize first letter of each (space-separated) word
! Streq              Generalized strings "==" (optionally ignoring case,
!                      leading blanks, and allowing wildcard matches)
! Stretch            Insert spaces between words; optionally between letters, too
! Strings2Ints       Converts an array of strings to ints using "ichar" ftn
! Swap               Swaps every instance of c1 with c2 and vice versa
! Trim_safe          trims string down, but never to length 0
! TrueList           describe where elements of an array are true
! UnAsciify          restore non-ascii characters in place of their coded forms
!                       (see Asciify)
! UnWrapLines        undo the splitting of commands across multiple lines
! WriteIntAsBaseN    Converts an integer to base n representation
! WriteIntsToChars   Converts an array of ints to strings using Fortran write
! WriteRomanNumerals Converts an integer to Roman numeral (e.g. 9 to 'ix')
! === (end of toc) ===

! === (start of api) ===
! char* Asciify ( char* str, [char* how] )
! char* Capitalize (char* str)
! CatStrings ( char* strings(:), char* sep, char* stringsCat, &
!      & int L, [char* options] )
! int charToInt (char str)
! char* CompressString (char* str)
! int count_quotes (char* str, char lquote, char rquote)
! int count_words (char* str)
! char* delete (char* str, char ch, [int max])
! char* depunctuate (char* str)
! Enclosure ( char* str, int i1, int i2, [ char brackets(:)] )
! FlushArrayLeft ( char* a(:), char* b(:), [char* options] )
! int(:) indexes (char* string, char* substrings, [char* mode])
! ints2Strings (int ints(:,:), char* strs(:))
! log(:) isAllAscii( char* arg(:) )
! log isAlphabet( char arg )
! log(:) isComment( char* arg(:), [char Comment] )
! log(:) isDigits( char* arg(:) )
! log IsRepeat ( char* str, [char* subtring] )
! char* Justify ( char* str, int length, [char mode] )
! int lenTrimToAscii (char* str)
! int LinearSearchStringArray (char* list(:), char* string, 
!   [log caseInsensitive, [log testSubstring], [log listInString])
! char* LowerCase (char* str)
! int(:) NAppearances (char* string, char* substrings)
! int nCharsinFormat ( char* Format )
! int NCopies (char* str, char* substring, [log overlap])
! ReadCompleteLineWithoutComments (int unit, char* fullLine, [log eof], &
!       & [char commentChar], [char continuationChar])
! ReadNumFromBaseN ( char* strs, real float, int N, [char* options] )
! ReadIntFromBaseN ( char* strs, int int, int N, [char* options] )
! ReadIntsFromChars ( char* strs[(:)], int ints[(:)], &
!       & [char* forbiddens], [char* ignore] )
! ReadNumsFromChars ( char* strs[(:)], num num[(:)], &
!       & [char* forbiddens], [char* ignore] )
! ReadRomanNumerals ( char* strs, int int )
! char* remap ( char* str, char* old, char* new )
! char* Replace ( char* str, char oldChar, char newChar, &
!      & [int max], [log reverse] )
! char* ReplaceNonAscii ( char* str, char newChar, [char* exceptions] )
! char* Reverse ( char* str )
! char* Reverse_trim ( char* str )
! char* Rot13 ( char* str, [int nn], [char* otp], [log inverse] )
! char* shiftLRC ( char* str, [char* position], [char fillChar] )
! int size_trim ( char* str(:) )
! SplitDetails ( char *str, char* switch, int details )
! SplitNest ( char *str, char* part1, char* part2, char* part3, [char* parens] )
! SplitWords (char *line, char* first, char* rest, [char* last], &
!       & [log threeWay], [char* delimiter])
! char* squeeze (char* str, [char* options])
! char* StartCase (char* str, [char separator])
! log streq ( [char* str1], [char* str2], [char* options] )
! log streq ( char* str1(:), char* str2, [char* options] )
! log streq ( char* str1, char* str2(:), [char* options] )
! char* stretch ( char* str, [char* options] )
! strings2Ints ( char* strs(:), int ints(:,:) )
! char* Swap ( char* str, char c1, char c2 )
! char* trim_safe ( char* str )
! TrueList ( logical* list, char* str )
! char* unAsciify ( char* str, [char* how] )
! unWrapLines ( char* inLines(:), char* outLines(:), &
!              [int nout], [char escape], [char comment] )
! writeIntAsBaseN (int int, int N, char* strs, [char* options])
! writeIntsToChars (int ints(:), char* strs(:))
! writeRomanNumerals (int int, char* strs)
! Many of these routines take optional arguments that greatly modify
! their default operation

! One standard is the character flag "options" which affects how loosely
! string matches may be interpreted
! it may include any of the following (poss. in combination, e.g. "-wc")
! w    Wildcard * which allows 'a*' to equal 'abcd'
! c    case insensitive which allows 'ABCD' to equal 'abcd'
! e    permit empty string elements
! f    flush left which allows 'abcd' to equal '  abcd'
! h    check only that the heads match; allows 'oldman' to match 'old'
! n    reverse sense of match (where appropriate)
! (These are different in streq_array1, however--either redo that function
! to make it conform, or rename the options flag there to prevent
! unnecessary confusion)

! We hope eventually that options will replace the countEmpty, caseSensitive, 
! etc. separate optional args to many of the current module procedures

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
! HHMMSS_value              Count_words
! Indexes                   LowerCase
! IsRepeat                  Reverse
! LinearSearchStringArray   Strings2Ints
! Nappearances              Trim_safe
! Ncopies
! ReadIntsFromChars         
! ReFormatDate              
! ReFormatTime              
! Reverse_trim              
! SplitWords                
! Streq                     
! === (end of api) ===

  public :: Asciify, &
    & Capitalize, CatStrings, CharToInt, CompressString, Count_Words, &
    & Delete, Depunctuate, Enclosure, FlushArrayLeft, Hhmmss_Value, &
    & Indexes, Ints2Strings, &
    & IsAllAscii, IsAlphabet, IsComment, IsDigits, IsRepeat, Justify, &
    & LenTrimToAscii, LinearSearchStringArray, Lowercase, &
    & NAppearances, NCharsInFormat, NCopies, &
    & ReadCompleteLineWithoutComments, ReadIntFromBaseN, ReadIntsFromChars, &
    & ReadNumFromBaseN, ReadNumsFromChars, ReadRomanNumerals, &
    & Remap, Replace, ReplaceNonAscii, Reverse, Reverse_Trim, &
    & Rot13, &
    & Shiftlrc, Size_Trim, SplitDetails, SplitNest, SplitWords, Squeeze, &
    & StartCase, Streq, Stretch, Strings2Ints, &
    & Swap, Trim_Safe, TrueList, UnAsciify, UnWrapLines, &
    & WriteIntasBaseN, WriteIntsToChars, WriteRomanNumerals

  interface ReadNumFromBaseN
    module procedure ReadFloatFromBaseN, ReadIntFromBaseN
  end interface

  interface ReadNumsFromChars
    module procedure ReadADoubleFromChars, ReadDoubleArrayFromChars
    module procedure ReadARealFromChars, ReadRealArrayFromChars
    module procedure ReadAnIntFromChars, ReadIntArrayFromChars
  end interface

  interface Streq
    module procedure Streq_Scalar, Streq_Array1, Streq_Array2
  end interface

  interface WriteIntsToChars
    module procedure WriteAnIntToChars, WriteIntsToChars_1d, WriteIntsToChars_2d
  end interface

  ! Public data
  character(len=16), public, save :: StringOptions = ' '

  ! hhmmss_value
  integer, public, parameter      :: Invalidhhmmssstring = 1
  ! readAnIntFromChars
  integer, public, parameter      :: StringContainsForbiddens=-999
  ! strings2Ints
  integer, public, parameter      :: LenOrSizeTooSmall=-999
  ! streq max input str lengths
  ! (We had to resort to this hard-wired limit after we made
  ! the possibly bone-headed decision to let
  ! str1 and str2 be optional args)
  integer, public, parameter      :: MaxStreqLen  = 2048
  integer, public, parameter      :: MaxSubstrLen = 256

  logical, private, save          :: CaseSensitive       
  logical, private, save          :: IgnoreLeadingSpaces

contains

  ! -------------------------------------------------  CatStrings  -----
  subroutine CatStrings ( Strings, Sep, StringsCat, L, options )
  ! Concatenate Strings with Sep between them, giving StringsCat(:L-1)
    character(len=*), intent(in) :: Strings(:)
    character(len=*), intent(in) :: Sep
    character(len=*), intent(out) :: StringsCat
    integer, intent(out) :: L
    character(len=*), optional, intent(in) :: options
    ! Internal variables
    logical :: CountEmpty
    integer :: I, N, T, W
    ! Executable
    !countEmpty = ( index( Default( Options, StringOptions), 'e') > 0 )
    countEmpty = index( StringOptions, 'e' ) > 0
    if ( present(options) ) countEmpty = index( options, 'e' ) > 0
    w = len(sep)
    l = len_trim(strings(1)) + 1
    stringsCat(:l-1) = strings(1)(:l-1)
    do i = 2, size(strings)
      t = len_trim(strings(i))
      if ( t < 1 .and. .not. CountEmpty ) cycle
      n = l + t + w
      n = min( len(StringsCat) + 1, n )
      stringsCat(l:n-1) = sep // strings(i)(:t)
      l = n
      if ( l > len(StringsCat) ) return
    end do
  end subroutine CatStrings

  ! ------------------------------------------------  CharToInt  -----
  elemental function CharToInt (str) result (int)
    ! This converts the input character to its integer value, if
    ! it has one
    ! e.g., '9' returns 9
    ! 'a' or '0' or ' ' , indeed anything else all return 0
    !--------Argument--------!
    character (len=1), intent(in) :: str
    !---------result---------!
    integer :: int
    !-----local-variables------!
    character(len=9), parameter :: allChars = &
      & '123456789'
    !-------executable-code----!
    int = index( allChars, str )
  end function CharToInt

  ! ---------------------------------------------  CompressString  -----
  function CompressString ( str ) result ( outstr )

    ! Removes all leading and embedded blanks from a string.
    ! This undoes what the function stretch does.
    ! Should we rename it 'unstretch'?
    ! Well, a truer implementation of an unstretch function would
    ! (1) turn 'a b' -> 'ab'
    ! (2) turn 'a  b' -> 'a b'
    ! (3) turn 'a  .. b' -> 'a b'
    ! i.e., turn every sequence of 2 or more blanks into a single blank
    ! while deleting any singleton blank.
    ! 
    !--------Argument--------!

    character (len=*), intent(in) :: str
    character (len=len(str)) :: outstr

    !----------Local vars----------!
    integer :: i, n

    !----------Executable part----------!

    outstr = " "
    n = 0
    DO i = 1, LEN(str)
       IF (str(i:i) /= " ") THEN
          n = n + 1
          outstr(n:n) = str(i:i)
       END IF
    END DO

  end function CompressString

  ! ------------------------------------------------  count_quotes  -----
  function count_quotes ( str, lquote, rquote ) result (no_of_quotes)
    ! This counts the number of quotes in a string 
    ! For our purposes, quotes consist of any non-space characters
    ! surrounded by an {lquote. rquote} pair

    ! method:
    ! step through string, each stride positioned between
    ! j1, where leftmost lquote is, and j2, where leftmost rquote is
    ! cut off both stride and evrything leftward, call remainder "remainder"
    ! Contrast this approach with what is is done in MLSStringLists module
    ! which uses string slicing and index arithmetic
    ! Should we go back and simplify the latter?
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=1), intent(in) :: lquote
    character (len=1), intent(in) :: rquote
    !---------result---------!
    integer::no_of_quotes
    !-----local-variables------!
    integer::j1
    integer::j2
    character(len=len(str)) :: remaining
    !-------Executable-code----!
    no_of_quotes = 0
    IF (index(str, lquote) < 1 .or. index(str, rquote) < 1) return
    remaining = str
    DO
       j1 = index( remaining, lquote)
       if ( j1 < 1 ) exit
       j2 = index( remaining, rquote)
       if ( j2 < 1 ) exit
       no_of_quotes = no_of_quotes + 1
       if ( j2 >= len_trim(remaining) ) exit
       remaining = remaining(j2+1:)
    END DO
  end function count_quotes

  ! ------------------------------------------------  COUNT_WORDS  -----
  function count_words (str) result (no_of_words)
    ! This counts the number of words in a string 
    ! For our purposes, words consist of any non-space characters
    ! and are separated by one or more spaces
    ! -----Added by HCP-------- 
    !--------Argument--------!
    character (len=*), intent(in) :: str
    !---------result---------!
    integer::no_of_words
    !-----local-variables------!
    integer::j
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
  end function count_words

  ! ------------------------------------------------  Delete  -----
  function Delete ( str, ch, max ) result ( outstr )
    ! function that removes every instance of a char
    !--------Argument--------!
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: ch
    integer, optional, intent(in) :: max ! up to how many such deletions?
    character(len=len(str)) :: outstr
    !----------Local vars----------!
    integer :: i, iout
    integer :: myMax
    integer :: dels
    !----------Executable part----------!
    outstr = str
    if ( index(str, ch) < 1 ) return
    myMax = Huge(0) / 2
    if ( present(max) ) myMax = max
    dels = 0
    outstr = ' '
    iout = 0
    do i = 1, len_trim(str)
      if ( str(i:i) == ch ) then
        dels = dels + 1
        if ( dels <= myMax ) cycle
      endif
      iout = iout + 1
      outstr(iout:iout) = str(i:i)
    end do

  end function Delete

  ! ------------------------------------------------  DEPUNCTUATE  -----
  function Depunctuate ( str ) result ( outstr )
    ! function that removes punctuation and replaces with blanks
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

  end function Depunctuate

  ! ------------------------------------------------  Enclosure  -----
  subroutine Enclosure ( str, i1, i2, brackets )
    ! Find outermost substring indices enclosing an arg
    ! with "(" and ")"; any failure returns 0 and 0
    ! E.g., if
    !   str = 'abs(a+(b+c)/d)'
    !          123456789012345
    ! returns i1=4 and i2=14
    ! Optionally will use 
    !   brackets(1) instead of "(" and 
    !   brackets(2) instead of ")"
    !
    ! Method:
    ! Each '(' adds 1, each ')' subtracts 1, though never going negative
    ! If 0 when done, we're balanced, not if not
    !
    ! brackets(1) and (2) must not be th same
    ! Args
    character(len=*), intent(in)                  :: str
    integer, intent(out)                          :: i1
    integer, intent(out)                          :: i2
    character, dimension(:), optional, intent(in) :: brackets
    ! Internal variables
    integer                                       :: i
    integer                                       :: balance
    character                                     :: lef ! "("
    character                                     :: rig ! ")"
    ! Executable
    lef = '('
    rig = ')'
    i1 = 0
    i2 = 0
    if ( present(brackets) ) then
      if ( size(brackets) > 1 ) then
        lef = brackets(1)
        rig = brackets(2)
      endif
    endif
    if ( lef == rig ) then
      ! print *, 'The lef and rig brackets must not be the same'
      return
    endif
    i1 = index( str, lef )
    if ( i1 < 1 ) then
      ! Failed right off the bat to find a 1st "("
      i1 = 0
      i2 = 0
      return
    endif
    balance = 1
    do i = index( str, lef ) + 1, len_trim(str)
      select case(str(i:i))
      case ('(')
        balance = balance + 1
      case (')')
        balance = max( 0, balance - 1 )
      ! case default (leaves balance unchanged)
      end select
      if ( balance < 1 ) exit
    enddo
    if ( i > len_trim(str) ) then
      ! Failed to find the matching ")"
      i1 = 0
      i2 = 0
    else
      i2 = i
    endif
  end subroutine Enclosure

  ! ------------------------------------------------  FlushArrayLeft  -----
  subroutine FlushArrayLeft ( a, b, options )
    ! Flush array "a" over by leading blanks, returning result as "b"
    ! according to options
    ! options           meaning
    ! -------           -------
    !  a (default)     array-wise: skip over leading blank array elements
    !  e               element-wise: skip over leading blanks in each element
    ! thus, -ae would both skip over leading blank array elements
    ! and also flush each element left
    !
    !--------Argument--------!
    character (len=*), dimension(:), intent(in)  :: a
    character (len=*), dimension(:), intent(out) :: b
    character (len=*), optional, intent(in)      :: options
    !---------local variables---------!
    integer::i, n
    character(len=8) :: myOptions
    !-------Executable-code----!
    myOptions = '-a'
    if ( present(options) ) myOptions = options
    if ( index( lowercase(myOptions), 'a' ) > 0 ) then
      b = ' '
      n = 0
      do i=1, size(a)
        if ( a(i) == ' ' .and. n < 1 ) cycle
        n = n + 1
        b(n) = a(i)
      enddo
    else
      b = a
    endif
    if ( index( lowercase(myOptions), 'e' ) < 1 ) return
    do i=1, size(a)
      b(i) = adjustl(b(i))
    enddo
  end subroutine FlushArrayLeft

  ! ------------------------------------------------  HHMMSS_value  -----
  function HHMMSS_value ( str, ErrTyp, separator, strict ) result ( value )
    ! function that returns the value in seconds of a string 'hh:mm:ss'
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
    double precision :: value
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
    value = 0.0
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
    ! Returns the array of indexes of each element of substrings in str
    ! The mode determines how consecutive elements of array are ordered
    ! (default is first)
    ! mode          order
    ! ------        -----
    ! first         always find first occurrence of substrings(i) in str
    ! last          always find last occurrence of substrings(i) in str
    ! left          progressively find next substrings(i) after substrings(i-1)
    ! right         progressively find next substrings(i) before substrings(i-1)
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
    ! we trim each element of the sub-strings--should we allow an option not to?
    ! (Probably)
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
       lpos = left(i) + max(len_trim(substrings(i)), 1) ! len(trim_safe(substrings(i)))
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

  ! ------------------------------------------------  Justify  -----
  function Justify ( str, length, mode ) result ( Justified )
    ! Return the string str Justified to fit the supplied length
    ! by padding its interior spaces. 
    ! Thus it turns
!The first line                                                                *
    ! into 
!The                                 final                                 line*
   
    ! Note:
    !      meaning both left-right justified; we aim to distribute
    !      any padding fairly, by which we mean not favoring gaps due their
    !      position but only their occurrence number.
    !      An alternate strategy would be based on the substring index number
    !      instead of the gap occurrence number. Could we implement both
    !      and allow the user to choose? Let's use the optional arg 'mode':
    !      mode           means
    !      ----           -----
    !       g             count gaps
    !       s             count substrings
    !
    ! See also Stretch
    ! Args
    character (len=*), intent(in)                :: str
    integer, intent(in)                          :: length
    character, intent(in), optional              :: mode ! 's' or 'g'
    character (len=length)                       :: justified
    ! Local variables
    character(len=length)                        :: newStr ! A temp
    integer                                      :: i ! substring index of str
    integer                                      :: j ! index of last gap
    integer                                      :: k ! index of justified
    character                                    :: myMode
    integer                                      :: n ! integer part of x
    integer                                      :: nWords
    logical                                      :: stillInGap ! Was last char blank?
    integer                                      :: surplus ! how many extra
    real                                         :: dx ! density of extra
    real                                         :: x  ! accumulated extra
    ! Executable
    justified = str
    ! print *, 'length, len_trim(justified): ', length, len_trim(justified)
    if ( length <= len_trim(justified) ) return
    nWords = count_words ( str )
    ! print *, 'nWords: ', nWords
    if ( nWords < 2 ) return
    myMode = 'g' ! Defaults to counting gaps, not substrings
    if ( present(mode) ) myMode = mode
    newStr = adjustl(str)
    surplus = length - len_trim(newstr)
    if ( myMode == 'g' ) then
      ! dx is the avg number of spaces to add to each gap
      dx = surplus / (nWords - 1.0)
    else
      ! dx is the avg number of spaces to add at each char
      dx = surplus / real(len_trim(newStr))
    endif
    ! print *, 'surplus, dx: ', surplus, dx
    ! We'll accumulate them at each gap, 
    ! acting only when we have enough to add an extra space
    x = 0.
    stillInGap = .false.
    k = 0
    do i=1, len_trim(newStr)
      if ( myMode /= 'g' ) x = x + dx
      if ( newStr(i:i) /= ' ' .or. stillInGap ) then
        k = k + 1
        justified(k:k) = newStr(i:i)
        stillInGap = ( newStr(i:i) == ' ' )
        if ( myMode /= 'g'  .and. x > 1. .and. stillInGap ) then
          ! Yes, we must add more spaces. How many?
          n = x
          ! print *, 'Adding spaces: ', n
          justified(k+1:k+n) = ' '
          x = x - n
          k = k + n
        endif
        cycle
      endif
      ! We're starting a new gap--must we add spaces?
      if ( myMode == 'g' ) x = x + dx
      if ( x > 1. ) then
        ! Yes, we must. How many?
        n = x
        ! print *, 'Adding spaces: ', n
        justified(k+1:k+n) = ' '
        x = x - n
        k = k + n
      endif
      stillInGap = .true.
      k = k + 1
      justified(k:k) = newStr(i:i)
    enddo
    ! Did we end at length?
    ! Probably not, so must add more spaces
    ! print *, 'k, length: ', k, length
    if ( k >= length ) return
    newStr = justified
    i = FindFirst( newStr, ' ' )
    n = length - k
    ! print *, 'i, n: ', i, n
    if ( nWords < 3 .or. n < 2 ) then
      ! crowd all n extra spaces inside the first gap
      justified = ' '
      justified = newStr(1:i) // repeat( ' ', n ) // newStr(i+1:)
    else
      ! Put (n/2) extra spaces inside the first gap
      ! and the other (n/2) inside the last gap
      j = FindLast( trim(newStr), ' ' )
      justified = ' '
      justified = newStr(1:i) // repeat( ' ', n/2 ) // newStr(i+1:j) // &
        & repeat( ' ', n - (n/2) ) // newStr(j+1:)
    endif
  end function Justify

  ! ------------------------------------------------  lenTrimToAscii  -----
  function lenTrimToAscii (str) result (trimmedLength)
    ! Return the len_trim of a string treating all non-Ascii as blanks
    ! -----added by hcp-------- 
    !--------argument--------!
    character (len=*), intent(in) :: str
    !---------result---------!
    integer :: trimmedLength
    !-----local-variables------!
    character(len=len(str)) :: newStr
    !-------executable-code----!
    newStr = ReplaceNonAscii( str, char(32) )
    trimmedLength = len_trim( newStr )
  end function lenTrimToAscii

  ! ------------------------------------  LinearSearchStringArray  -----

  ! This routine does a simple linear search for a string in an array.
  ! If the case insensitive flag is set the strings are capitalized first.
  ! If the test substring flag is set, the string is tested as a substring.
  ! If the listInString flag is set, the array list is tested as substrings
  !  against the string; otherwise, the string is tested as a substring
  !  against the array list.
  ! If the string is not found, 0 is returned

  function LinearSearchStringArray (list, string, caseInsensitive, &
       & testSubstring, listInString) result (sindex)

    ! Dummy arguments
    character (len=*), dimension(:) :: list
    character (len=*) :: string
    logical, intent (in), optional :: caseInsensitive
    logical, intent (in), optional :: testSubstring
    logical, intent (in), optional :: listInString


    ! function result
    integer :: sindex   ! matching string index (0 = not found)

    ! Local variables
    integer :: i
    logical :: useCaseInsensitive
    logical :: testForSubstring
    logical :: testForList
    logical :: found

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

  end function LinearSearchStringArray

  ! ---------------------------------------------------  NAppearances  -----
  function NAppearances( str, substrings, dontTrim ) result(array)
    character(len=*), intent(in) :: str
    character(len=*), intent(in), dimension(:) :: substrings
    logical, optional, intent(in) :: dontTrim
    integer, dimension(size(substrings)) :: array
    ! Returns the array of the number of times each element of substrings
    ! appears in str
    ! E.g., if str='ababababababa' and substrings = (/'a', 'ab', 'abab', 'b'/)
    !      result
    !      ------
    !   (/7, 6, 3, 6/)
    !
    ! Method:
    ! Use indexes function to find successive indexes of a single substring
    !
    ! Note:
    !     these are distinct, non-overlapping occurrences of each sub-string
    !     the dontTrim option not yet passed to indexes
    ! Internal variables
    integer, dimension(len(str)) :: tmpArray
    character(len=len(substrings)), dimension(len(str)) :: tmpSubs
    integer :: i
    logical :: myDontTrim
    ! Executable
    myDontTrim = .false.
    if ( present(dontTrim) ) myDontTrim = dontTrim
    do i=1, size(substrings)
      tmpArray = 0
      tmpSubs = substrings(i)
      if ( myDontTrim ) then
        tmpArray = indexes( str, tmpSubs, 'left' )
      else
        tmpArray = indexes( trim_safe(str), tmpSubs, 'left' )
      endif
      array(i) = count( tmpArray > 0 )
    enddo
  end function NAppearances

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

  subroutine ReadCompleteLineWithoutComments(unit,fullLine,eof, &
       & commentChar,continuationChar)

    ! Dummy arguments

    integer, intent(in) :: unit ! Input file unit
    ! fullLine changed to intent InOut by HCP. Some (but not all) 
    ! F90 compilers won't let this be intent(out) because the declaration
    ! of inputLine makes use of the length of fullLine even if its _contents_
    ! are immaterial
    character(len=*), intent(inout) :: fullLine ! Output line
    character(len=*), optional :: commentChar
    character(len=*), optional :: continuationChar
    logical, intent(out), optional :: eof ! Set if got to eof

    ! Local variables

    integer :: ioInfo           ! IOSTAT result
    character(len=len(fullLine)) :: inputLine ! One line from file
    integer :: commentStart     ! Start of a comment in line
    integer :: lastChar         ! Last character position in line
    integer :: gotContinuation  ! 1 if continuation needed, 0 if not
    logical :: firstLine        ! A correction to be applied

    character(len=1) :: useCommentChar
    character(len=1) :: useContinuationChar

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

  end subroutine ReadCompleteLineWithoutComments

  ! ----------------  readNumFromBaseN  -----
  ! read a number which had been reexpressed in base n from a string
  
  ! Use:
  ! Return e.g. the number 1989 in base 5 from str='3 0 4 2 4'
  ! i.e., 4 + 5*(2 + 5*(4 + 5*(0 + 5*3)))
  ! According to the options string, the intermediate integer coefficients
  ! will be based on
  !  options contain            coefficient string
  !    ---                      ------------------
  !   (default)                 decimal integers, separated by spaces
  !                              e.g., ''13 4 29 6''
  !      x                      "Extended" hexadecimal-style (hex if N = 16)
  !                               e.g., '9 b d 2 f 0 0 0'
  !      c                      compressed str; otherwise like 'x'
  !                               i.e. no spaces separate each digit
  !                               e.g., '9bd2f000'
  ! 0123456789012345678901234567890123456789012345678901234567890123456789
  ! 0123456789abcdefghijklmnopqrstuvwxyz;'[]ABCDEFGHIJKLMNOPQRSTUVWXYZ:"{}
  !                              (good only up to base 70)
  !
  ! The default representation can handle arbitrarily large moduli
  ! a or A is limited to the unique range of achar
  ! x is good only up to a modulus of 70.

  ! Notes
  ! (1) we do not check that the str is a valid number in base N
  !     E.g., you can enter '10 0 1' in base 5 and get back a numerical
  !     result of 251 instead of an error message
  ! (2) we do not check that str would be converted to an integer
  !     outside the largest representable
  ! (3) At first we intended to include the following options; however
  !     it never became clear that there was an urgent need they satisfied
  !      a                      achar(c[m]) in the same order
  !                              e.g., CR // EOT // GS // ACR   
  !      A                      like a but converted to ascii
  !                              e.g., '<CR>' // '<EOT>' // '<GS>' // <ACR>   

  subroutine readIntFromBaseN ( STR, K, N, OPTIONS )
    ! reads an integer from a string
    ! using poorly-tested but hopefully non-critical code
    
    ! Note: beware of integer overflow; if str would be converted
    ! to an int > max integer representable, you must
    ! use the api which reads a float from str instead

    !--------Argument--------!
    character (len=*), intent(in)           ::   str
    integer, intent(out)                    ::   k ! integer result
    integer, intent(in)                     ::   N ! the base
    character(len=*), optional, intent(in)  ::   options
    include "ReadNumFromBaseN.f9h"
  end subroutine readIntFromBaseN

  subroutine readFloatFromBaseN ( STR, k, N, OPTIONS )
    ! reads an f.p. number from a string
    ! using poorly-tested but hopefully non-critical code

    !--------Argument--------!
    character (len=*), intent(in)           ::   str
    real, intent(out)                       ::   k ! integer result
    integer, intent(in)                     ::   N ! the base
    character(len=*), optional, intent(in)  ::   options
    include "ReadNumFromBaseN.f9h"
  end subroutine readFloatFromBaseN

  subroutine readARealFromChars (str, num, forbiddens, ignore)
    !--------Argument--------!
    character (len=*), intent(in) ::   str
    real, intent(out)             ::   num
    include 'ReadANumFromChars.f9h'
  end subroutine readARealFromChars

  subroutine readADoubleFromChars (str, num, forbiddens, ignore)
    !--------Argument--------!
    character (len=*), intent(in) ::   str
    double precision, intent(out)             ::   num
    include 'ReadANumFromChars.f9h'
  end subroutine readADoubleFromChars

  subroutine ReadRealArrayFromChars (strs, nums, forbiddens, ignore)
    character (len=*), intent(in), dimension(:) ::   strs
    real, intent(out), dimension(:)             ::   nums
    character (len=*), intent(in), optional     ::   forbiddens
    character (len=*), optional, intent(in)     :: ignore

    !----------Local vars----------!
    integer :: i, arrSize
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   arrSize = MIN(size(strs), size(nums))
   if ( arrSize <= 0 ) then
     nums = LENORSIZETOOSMALL
     return
   endif
   do i=1, arrSize
     call readARealFromChars(strs(i), nums(i), forbiddens, ignore)
   enddo

  end subroutine ReadRealArrayFromChars

  subroutine ReadDoubleArrayFromChars (strs, nums, forbiddens, ignore)
    character (len=*), intent(in), dimension(:) ::   strs
    double precision, intent(out), dimension(:) ::   nums
    character (len=*), intent(in), optional     ::   forbiddens
    character (len=*), optional, intent(in)     :: ignore

    !----------Local vars----------!
    integer :: i, arrSize
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   arrSize = MIN(size(strs), size(nums))
   if ( arrSize <= 0 ) then
     nums = LENORSIZETOOSMALL
     return
   endif
   do i=1, arrSize
     call readADoubleFromChars(strs(i), nums(i), forbiddens, ignore)
   enddo

  end subroutine ReadDoubleArrayFromChars

  ! ----------------  readRomanNumerals  -----
  ! Read a string composed of roman numerals into an int
  
  ! Use:
  ! Given e.g. str='MCMLXXXIX ' return in num the value 1989
  subroutine readRomanNumerals ( str, num )
    ! takes a string and returns an integer
    ! using poorly-tested but hopefully non-critical code

    ! Method:
    ! Process string left to right looking first for highest-valued char
    ! then working down to lowest valued
    ! When ever a lower-valued char is found to left of a higher valued
    ! one, it and all chars up to higher value are azxsigned negative value
    ! e.g., 'cm' is -100 + 'm', or 'xm' i -10 + 'm'
    ! The second example is non-standard roman numerals; we can read them
    ! anyway

    ! Limitation: we treat upper case and lowercase equivalently
    ! we ignore non-roman strings
    !
    !--------Argument--------!
    character (len=*), intent(in) ::   str
    integer, intent(out)          ::   num
    ! Internal variables
    character(len=1), dimension(7), parameter :: romans = &
      & (/  'm', 'd', 'c', 'l', 'x', 'v', 'i' /)
    integer, dimension(7), parameter :: values = &
      & (/ 1000, 500, 100, 50,  10,   5,   1 /)
    integer :: strpos  ! What str character we're on
    integer :: r_index ! What romans index we're searching for
    ! Executable
    num = 0
    if ( len_trim(str) < 1 ) return
    ! Outer loop: string character number
    do strpos=1, len_trim(str)
      ! Inner loop: romans index
      do r_index=1, 7
        if ( lowercase(str(strpos:strpos)) /= romans(r_index) ) cycle
        ! OK, we've found it--but is there a higher-valued one to the right?
        if ( r_index == 1 .or. strpos == len_trim(str) ) then
          num = num + values(r_index)
        elseif ( any(indexes( lowercase(str(strpos+1:)), romans(1:r_index-1)) > 0 ) ) then
          num = num - values(r_index)
        else
          num = num + values(r_index)
        endif
      enddo
    enddo
  end subroutine readRomanNumerals

   ! --------------------------------------------------  remap  -----
  function remap ( str, old, new ) result (outstr)
    ! Replaces every old in the input str with the corresponding new char
    ! E.g., given 
    ! str = 'Four score and seven years' and
    ! old = 'osa' and new = '0$9', then the result will be
    ! result = 'F0ur $c0re 9nd $even year$'
    ! Notes
    ! (1) characters in str not found in old will be left unchanged
    ! (2) if new is shorter than old, some old characters may be replaced by blanks
    ! (3) See also Rot13 and Replace
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=*), intent(in) :: old
    character (len=*), intent(in) :: new
    character (len=max(len_trim(str), 1)) :: outstr

    !----------local vars----------!
    integer :: i
    integer :: it
    !----------executable part----------!
    outstr = str
    if( len_trim(str) < 1 .or. len(old) < 1 ) return
    do i = 1, len_trim(str)
       it = index( old, str(i:i) )
       if ( it > 0 .and. it < len_trim(new)+1 ) then
         outstr(i:i) = new(it:it)
       elseif ( it > len_trim(new) ) then
         outstr(i:i) = ' '
       endif
    end do

  end function remap

   ! --------------------------------------------------  Replace  -----
  function Replace ( str, oldChar, newchar, max, reverse ) result ( outstr )
    ! takes a string and returns one with oldChar replaced by newChar
    ! E.g., to replace every char(0), which is the NUL character, with a blank
    ! arg = Replace( arg, char(0), char(32) )
    !
    ! If TRUE, the optional arg reverse replaces all the chars /= oldChar
    ! E.g., the following can be used to underline chars
    !  underlineChars = Replace( chars, ' ', '-', reverse=.true. )
    character(len=*), intent(in)  :: str
    character(len=1), intent(in)  :: oldChar
    character(len=1), intent(in)  :: newChar
    integer, optional, intent(in) :: max ! up to how many such replacements?
    logical, optional, intent(in) :: reverse ! Replace if /= oldChar?
    character(len=len(str))       :: outstr
    ! Internal variables
    integer :: i, n
    integer :: myMax
    logical :: myReverse
    integer :: reps
    ! Executable
    outstr = str
    if ( len(str) < 1 ) return
    if ( index(str, oldChar) < 1 ) return
    myMax = Huge(0) / 2
    if ( present(max) ) myMax = max
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    n = len(str)
    reps = 0
    do i=1, n
      if ( .not. myReverse ) then
        if ( str(i:i) == oldChar ) then
          reps = reps + 1
          if ( reps > myMax ) exit
          outstr(i:i) = newChar
        endif
      else
        ! Reverse sense: replace only the chars /= oldChar
        if ( str(i:i) /= oldChar ) then
          reps = reps + 1
          if ( reps > myMax ) exit
          outstr(i:i) = newChar
        endif
      endif
    enddo
  end function Replace

   ! --------------------------------------------------  Reverse  -----
  elemental function Reverse (str) result (outstr)
    ! takes a string and returns one with chars in reversed order
    ! Useful in certain contexts:
    ! e.g., to remove leading blanks
    ! arg = Reverse(TRIM(Reverse(arg)))
    !
    ! See also Reverse_trim, ReverseList
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=len(str)) :: outstr

    !----------Local vars----------!
    integer :: i, istr, irev
    character (len=1) :: strChar
    !----------Executable part----------!
    outstr=str
    IF(LEN(str) == 1) RETURN

    DO i = 1, LEN(str)-1, 2
       istr = 1 + (i-1)/2                       ! 1, 2, ..
       irev = LEN(str) - (i-1)/2                ! N, N-1, ..
       strChar = str(istr:istr)
       outstr(istr:istr) = str(irev:irev)
       outstr(irev:irev) = strChar
    END DO

! Special case: str contains odd number of chars
    IF(MOD(LEN(str), 2) == 1) THEN
       istr = 1 + (LEN(str)-1)/2                ! 1, 2, ..
        outstr(istr:istr) = str(istr:istr)
    ENDIF

  end function Reverse

   ! --------------------------------------------------  Reverse_trim  -----
  function Reverse_trim (str) result (outstr)
    ! takes a string, trims it then returns one with chars in reversed order
    ! See also Reverse which omits the trim step
    !
    ! E.g., given 'A string    ' reverse_trim returns 'gnirst A   ' while
    ! a simple Reverse returns '   gnirst A'
    !
    ! Uses:
    ! To find the last non-blank character in a string str
    !   rev = reverse_trim(str)
    !   it = rev(1:1)
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=max(len_trim(str), 1)) :: outstr

    !----------Local vars----------!
    integer :: i, istr, irev
    character (len=1) :: strchar
    !----------Executable part----------!
    outstr=str
    IF(LEN_TRIM(str) <= 1) RETURN

    DO i = 1, LEN_TRIM(str)-1, 2
       istr = 1 + (i-1)/2                       ! 1, 2, ..
       irev = LEN_TRIM(str) - (i-1)/2           ! N, N-1, ..
       strChar = str(istr:istr)
       outstr(istr:istr) = str(irev:irev)
       outstr(irev:irev) = strChar
    END DO

! Special case: str contains odd number of chars
    IF(MOD(LEN_TRIM(str), 2) == 1) THEN
       istr = 1 + (LEN_TRIM(str)-1)/2           ! 1, 2, ..
        outstr(istr:istr) = str(istr:istr)
    ENDIF

  end function Reverse_trim

   ! --------------------------------------------------  Rot13  -----
  function Rot13 ( str, nn, otp, inverse ) result (outstr)
    ! performs generalized ROT13 on the input str
    ! E.g., given 
    ! 'Look, a cryptic message!' Rot13 returns 'Ybbx, n pelcgvp zrffntr!'
    ! Rot13 is self-inverting
    ! ROT(nn) is more general shift by nn, 0 < nn < 26
    ! If ROT(nn) and ROT(mm) are inverses, then nn + mm = 26
    
    ! The optional arg "inverse" if TRUE performs the inverse shift
    
    ! In principle use of a one-time pad makes an encrypted message
    ! unbreakable. The string of characters in otp will be used one-by-one
    ! to shift the input str
    !--------Argument--------!
    character (len=*), intent(in) :: str
    integer, optional, intent(in)         :: nn
    logical, optional, intent(in)         :: inverse
    character (len=*), optional, intent(in) :: otp
    character (len=max(len_trim(str), 1)) :: outstr

    !----------local vars----------!
    integer :: i, it, nprime, o, oprime
    logical :: inv
    !----------executable part----------!
    outstr = ' '
    if(len_trim(str) < 1) return
    nprime = 13
    if ( present(nn) ) nprime = nn
    inv = .false.
    if ( present(inverse) ) inv = inverse
    if ( inv ) nprime = 26 - nprime
    do i = 1, len_trim(str)
       o = iachar( str(i:i) )
       ! Leave non-printing characters alone
       if ( o < 32 ) then
         oprime = o
       elseif ( present(otp) ) then
         if ( len_trim(otp) < i ) then
           it = mod( i-1, len_trim(otp) ) + 1
           nprime = iachar( otp(it:it) )
         else
           nprime = iachar( otp(i:i) )
         endif
         if ( nprime < 32 ) nprime = 32 + nprime
         if ( .not. inv ) then
           oprime = 32 + mod(o-32 + nprime-32, 95)
         else
           oprime = 32 + mod(o-32 - nprime+32+95, 95)
         endif
       elseif ( 64 < o ) then
         if ( o < 91 - nprime ) then
           oprime = o + nprime
         elseif ( o < 91 ) then
           oprime = o + nprime - 26
         elseif ( o < 97 ) then
           oprime = o
         elseif ( o < 123 - nprime ) then
           oprime = o + nprime
         elseif ( o < 123 ) then
           oprime = o + nprime - 26
         else
           oprime = o
         endif
       else
         oprime = o
       endif
       outstr(i:i) = achar(oprime)
    end do

  end function Rot13

   ! --------------------------------------------------  shiftLRC  -----
  elemental function shiftLRC (str, position, fillChar) result (outstr)
    ! shifts the characters in a string with leading or trailing blanks to:
    ! position
    !  'l[eft]'    all the way to the left
    !  'c[enter]'  so that they are centered within the string
    !  'r[ight]'   all the way to the right (default)
    ! optionally, replaces the surrounding blanks with fillChar
    !--------Argument--------!
    character (len=*), intent(in)           :: str
    character (len=*), optional, intent(in) :: position
    character (len=1), optional, intent(in) :: fillChar
    character (len=len(str))                :: outstr

    !----------Local vars----------!
    integer :: i_average
    integer :: i_leading
    integer :: i_trailing
    character(len=1) :: myPosition
    !----------Executable part----------!
    outstr=str
    if( len_trim(str) < 2 .or. len(str) == len_trim(str) ) return
    myPosition = 'r'
    if ( present(position) ) myPosition = lowerCase(position(1:1))
    select case (myPosition)
    case ('l')
      outstr = adjustl(str)
      if ( present(fillChar) ) then
        i_trailing = len(outstr) - len_trim(outstr)
        outstr = trim(outstr) // repeat(fillChar, i_trailing)
      endif
    case ('r')
      outstr = adjustr(str)
      if ( present(fillChar) ) then
        i_leading = len_trim(outstr) - len_trim(adjustl(outstr))
        outstr = repeat(fillChar, i_leading) // adjustl(outstr)
      endif
    case ('c')
      ! Method: average the number of leading and trailing blanks
      i_leading = len_trim(str) - len_trim(adjustl(str))
      i_trailing = len(str) - len_trim(str)
      i_average = (i_leading+i_trailing)/2
      if ( .not. present(fillChar) ) then
        outstr = ' '
        outstr(i_average+1:) = adjustl(str)
      else
        outstr = repeat(fillChar, i_average) // trim(adjustl(outstr)) // &
          &         repeat(fillChar, i_average)
      endif
    case default
      outstr = adjustr(str)
      if ( present(fillChar) ) then
        i_leading = len_trim(outstr) - len_trim(adjustl(outstr))
        outstr = repeat(fillChar, i_leading) // adjustl(outstr)
      endif
    end select
  end function shiftLRC

  ! ------------------------------------------------  size_trim  -----
  function size_trim ( strs, safe ) result (length)
    ! This performs len_trim of character scalar equivalent to 
    !--------Argument--------!
    character (len=*), dimension(:), intent(in) :: strs
    logical, optional, intent(in)               :: safe
    !---------result---------!
    integer::length
    !-----local-variables------!
    character(len=size(strs)*len(strs)) :: together, ttogether
    !-------Executable-code----!
    together = transfer( strs, ttogether )
    length = len_trim( together )
    if ( .not. present(safe) ) return
    if ( safe ) length = max( length, 1 )
  end function size_trim

  ! ---------------------------------------------  SplitDetails  -----

  ! This subroutine splits an input string into 2 parts:
  ! a string which is the switch
  ! and an integer which is the details
  ! It assumes that they are catenated as in the example
  !     'pro1' => 
  !         {string = 'pro', 
  !          details = 1}
  ! This corresponds directly to how the Switches data is
  ! encoded and used by level 2
  ! Any other usefulness would be purest serendipity
  ! (a tempting name for this subroutine, but for once
  ! we succeeded in resisting temptation--for each 'serendipity' we
  ! would have a thousand 'crudeButUsefulHack's)

  subroutine SplitDetails( str, switch, details )
    character(len=*), intent(in)           :: str
    character(len=*), intent(out)          :: switch
    integer, intent(out)                   :: details
    ! Local variables
    character(len=*), parameter :: digits    = '1234567890'
    character(len=*), parameter :: nondigits = 'abcdefghijklmnopqrstuvwxyz' // &
      &                                        'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer                     :: firstDigitPosition
    integer                     :: lastNonDigitPosition
    ! Executable
    switch = ' '
    details = 0   ! If details is absent, this is its default
    if ( len_trim(str) < 1 ) return
    firstDigitPosition = scan( str, digits )
    if ( firstDigitPosition < 1 ) then
      ! No details present
      switch = str
    elseif ( firstDigitPosition == 1 ) then
      ! No switch present
      call readAnIntFromChars ( str, details )
    else
      ! OK, we can't naively assume what we found marks the switch detail start
      ! because there might be internal digits in the switch, e.g. 'l1bread'
      ! Therefore, we must determine the last non-digit position
      lastNonDigitPosition = scan( str, nondigits, back=.true. )
      firstDigitPosition = lastNonDigitPosition + 1
      if ( firstDigitPosition > len_trim(str) ) then
        ! No details present
        switch = str
      else
        switch = str(:firstDigitPosition-1)
        call readAnIntFromChars ( str(firstDigitPosition:), details )
      endif
    endif
  end subroutine SplitDetails

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

  elemental subroutine SplitWords(line,first,rest,last,&
       & threeWay,delimiter)

    ! Dummy arguments

    character (len=*), intent(in) :: line
    character (len=*), intent(out) :: first
    character (len=*), intent(out) :: rest
    character (len=*), intent(out), optional :: last

    logical, intent(in), optional :: threeWay
    character (len=*), intent(in), optional :: delimiter ! really separator

    ! Local variables

    character (len=1) :: useseparator
    logical :: useThreeWay
    character (len=len(line)) useLine ! Line with leading spaces removed

    integer :: firstseparatorPos,lastseparatorPos,trimmedLen

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

  end subroutine SplitWords
  
  ! ------------------------------------------------- squeeze --------
  ! Snip excess spaces between words
  ! E.g., turns ' a   man   from   Sai- Pan' into ' a man from Sai- Pan'
  ! optionally snips all spaces making it 'amanfromSai-Pan'
  
  ! options, if present, can contain the following characters
  !  character                 effect
  !     a                   snip all spaces (see CompressString)
  function squeeze( str, options ) result( squeezed )
    ! Args
    character(len=*), intent(in)           :: str
    character(len=*), optional, intent(in) :: options
    character(len=len(str))                :: squeezed
    ! Internal variables
    integer                                :: cpos ! current position in str
    integer                                :: cposq ! current position in squeezed
    logical                                :: newWord
    logical                                :: snipAll
    character(len=1)                       :: space
    ! Executable
    snipAll = .false.
    if ( present(options) ) snipAll = ( index(options, 'a') > 0 )
    space = ' '
    squeezed = str
    if ( len_trim(str) < 2 ) return
    squeezed = ' '
    if ( snipAll ) then
      if ( .true. ) then
        cposq = 0
        do cpos = 1, len_trim(str)
          ! This is easy -- snip every space no matter where
          if ( str(cpos:cpos) /= space ) then
            cposq = cposq + 1
            squeezed(cposq:cposq) = str(cpos:cpos)
          endif
        enddo
      else
        squeezed = CompressString( str )
      endif
    else
      squeezed(1:1) = str(1:1)
      cposq = 1
      newWord = ( str(1:1) == space )
      do cpos = 2, len_trim(str)
        if ( newWord ) then
          ! Already have at least one space, so must snip any others
          ! i.e., snip unless not a space
          if ( str(cpos:cpos) /= space ) then
            cposq = cposq + 1
            squeezed(cposq:cposq) = str(cpos:cpos)
            newWord = .false.
          endif
        else
          ! don't snip, even if a space
          cposq = cposq + 1
          squeezed(cposq:cposq) = str(cpos:cpos)
          ! Have we reached a space which divides words?
          newWord = ( str(cpos:cpos) == space )
        endif
      enddo
    endif
  end function squeeze
       
  ! -------------------------------------------------  StartCase  -----
  elemental function StartCase ( STR, SEPARATOR ) result (OUTSTR)
    ! Capitalize first letter of each (space-separated) word
    !--------Argument--------!
    character (len=*), intent(in)           :: STR
    character (len=1), optional, intent(in) :: SEPARATOR
    character (len=len(str))                :: OUTSTR
    ! Internal variables
    character(len=1) :: space
    integer :: i
    ! Executable
    outstr = lowercase(adjustl(str))
    if ( len_trim(str) < 1 ) return
    space = ' '
    if ( present(separator) ) space=separator
    outstr(1:1) = Capitalize(outstr(1:1))
    if ( index(trim(outstr), space) < 1 .or. len_trim(outstr) < 3 ) return
    do i = 3, len_trim(outstr)
       if( outstr(i:i) /= space .and. outstr(i-1:i-1) == space ) then
          outstr(i:i) = Capitalize(outstr(i:i))
       end if
    enddo
  end function StartCase

  ! -------------------------------------------------  streq_array1  -----
  function streq_array1 (STR1, STR2, OPTIONS) result (relation)
    ! Array version of streq
    ! May return multiple TRUEs (except see 's', 'l' options)
    ! Extra options
    ! 'P' is partial match (each STR1(i) is replaced by '*' // STR1(i) // '*')
    ! 'p' is partial match (STR2 is replaced by '*' // STR2 // '*')
    ! 's' returns TRUE only in element corresponding to shortest STR1
    ! 'l' returns TRUE only in element corresponding to longest STR1
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
    character(len=len(str1)+2) :: mystr1
    character(len=len(str2)+2) :: mystr2
    ! Executable
    myOptions = ''
    if ( present(options) ) myOptions = options
    mystr2 = str2
    if ( index(myOptions, 'p') > 0 ) mystr2 = '*' // trim_safe(str2) // '*'
    do i=1, size(str1)
      mystr1 = str1(i)
      if ( index(myOptions, 'P') > 0 ) mystr1 = '*' // trim_safe(str1(i)) // '*'
      relation(i) = streq( mystr1, mystr2, OPTIONS )
      ! print *, trim(mystr1)
      ! print *, trim(mystr2)
      ! print *, relation(i)
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
  function streq_scalar ( str1, str2, options ) result ( relation )
    ! Are two strings "equal" where equality is modified by
    ! (w) Wildcard * (off by default) which allows 'a*' to equal 'abcd'
    ! (c) case insensitive (off by default) which allows 'ABCD' to equal 'abcd'
    ! (f) flush left (off by default) which allows 'abcd' to equal '  abcd'
    ! (h) head match (off by default) which allows 'abcd' to equal 'abc '
    ! (n) reverse sense ("Are two strings unequal?")
    !
    ! Defaults to standard (str1 == str2), but options broaden cases of TRUE
    ! To turn options on, supply optional arg options which is a character
    ! string containing: 
    ! 'w' => turns on (1); 'c' => (2); 'f' => (3); 'n' => (4)
    ! e.g., streq('Ab*', 'abcd', '-wc') is TRUE
    !
    ! Notes:
    ! The '-' character in options is ignored and therefore not necessary
    ! A more powerful version can be imagined that would permit full regexp
    ! Trailing spaces are always ignored; e.g. streq('abcd ', 'abcd') is TRUE
    ! Only one of str1, str2 may contain wildcards
    !--------Argument--------!
    ! str1 and str2 are optional so we can streq on optional args
    ! and always return false
    ! (But should we return false even when the sense is reversed?)
    character (len=*), intent(in), optional  :: str1
    character (len=*), intent(in), optional  :: str2
    character (len=*), intent(in), optional  :: options
    logical                                  :: relation

    ! Internal variables
    integer, parameter :: MAXNUMWILDCARDS = 10 ! How many '*' in the pattern
    character(len=*), parameter :: star = '*'  ! Should we allow others?
    logical :: flushleft
    integer :: i
    logical :: ignorecase
    integer, dimension(MAXNUMWILDCARDS) :: istars
    character(len=8) :: myOptions
    integer :: nstars
    logical :: reverseSense
    integer :: spos
    ! If len(str) is used for ptrn and substrs, Intel Fortran 10.0.023 crashes
    character(len=MaxStreqLen)    :: str
    character(len=MaxStreqLen)    :: ptrn  ! The one with '*'
    character(len=MaxSubstrLen), dimension(MAXNUMWILDCARDS+1) :: substrs
    logical :: head
    logical :: wildcard
    logical, parameter :: deebug = .false.
    !----------Executable part----------!
    ! Check for any absent atr1, str2
    relation = ( present(str1) .and. present(str2) )
    if ( .not. relation ) return
    ! Check for length violations
    relation = ( len(str1) <= MaxStreqLen .and. len(str2) <= MaxStreqLen )
    if ( .not. relation ) return
    !
    relation = .FALSE.
    myOptions = ' '
    if ( present(options) ) myOptions = lowercase(options)
    wildcard = (index(myoptions, 'w') > 0)
    head = (index(myoptions, 'h') > 0)
    ignorecase = (index(myoptions, 'c') > 0)
    flushleft  = (index(myoptions, 'f') > 0)
    reverseSense  = (index(myoptions, 'n') > 0)
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
      goto 90 ! return
    elseif ( len_trim(str1) < 1 ) then
      relation = ( str2 == star .and. wildcard ) 
      goto 90 ! return
    elseif ( len_trim(str2) < 1 ) then
      relation = ( str1 == star .and. wildcard ) 
      goto 90 ! return
    elseif ( str2 == str1 ) then
      relation = .true.
      goto 90 ! return
    endif

    if ( head ) then
      ! Check only that the heads of each string match, i.e.
      ! up to the rightmost non-blank character of the shorter string
      if ( ignorecase ) then
        str  = lowercase(str1)
        ptrn = lowercase(str2)
      else
        str  = str1
        ptrn = str2
      endif
      if ( flushleft ) then
        str  = adjustl(str )
        ptrn = adjustl(ptrn)
      endif
      i = min( len_trim(str), len_trim(ptrn) )
      relation = ( str(:i) == ptrn(:i) )
      goto 90
    elseif ( .not. wildcard ) then
      if ( ignorecase ) then
        if ( flushleft ) then
          relation = (adjustl(lowercase(str1)) == adjustl(lowercase(str2)))
        else
          relation = (lowercase(str1) == lowercase(str2))
        endif
      else
        if ( flushleft ) then
          relation = (adjustl(str1) == adjustl(str2))
        ! else plain vanilla case already handled as special case above
        endif    
      endif    
      goto 90 ! return
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
    if ( deebug ) print *, 'str: ', trim(str), '  ptrn: ', trim(ptrn)
    if ( ptrn == star ) then
      relation = .true.
      goto 90 ! return
    elseif ( isRepeat(ptrn, star) ) then
      relation = .true.
      goto 90 ! return
    endif

    ! 1st -- how many stars?
    nstars =     ncopies(ptrn, star)
    if ( deebug ) print *, 'num of * ', nstars
    ! 2nd -- extract substrings from inbetween the wildcards
    substrs = star
    istars(1:nstars) = indexes(ptrn, substrs(1:nstars), mode='left')
    if ( deebug ) print *, 'where? ', istars(1:nstars)
    substrs = ' '
    spos = 1
    do i=1, nstars
      if ( spos > len(ptrn) ) exit
      substrs(i) = firstsubstr(ptrn(spos:), star)
      spos = max(spos, istars(i)) + 1
    enddo
    substrs(nstars+1) = Reverse_trim(firstsubstr(Reverse_trim(ptrn), star))
    if ( deebug ) then
      do i=1, nstars+1
        print *, trim_safe(substrs(i))
      enddo
    endif
    ! Deal specifically with empty elements of substrs
    if ( deebug ) print *, 'Deal specifically with empty elements of substrs'
    relation = .true.
    if ( deebug ) print *, 'Is substr(1) non-blank? ', substrs(1) /= ' '
    if ( substrs(1) /= ' ' ) then
      if ( deebug ) then
        print *, 'About to check index'
        print *, 'str ', str
        print *, 'trim(substrs(1)) ', trim(substrs(1))
        print *, 'len(str) ', len(str)
        print *, 'len(trim(substrs(1))) ', len(trim(substrs(1)))
        print *, 'index: ', index( str, trim(substrs(1)) )
      endif
      if ( index(str, trim(substrs(1))) /= 1 ) then
        relation = .false.
        goto 90 ! return
      endif
    endif
    if ( deebug ) print *, 'passed 1st sub-test'
    ! firstSSindex = 2
    if ( substrs(nstars+1) /= ' ' ) then
      if ( index(Reverse_trim(str), Reverse_trim(substrs(nstars+1))) /= 1 ) then
        relation = .false.
        ! if ( deebug ) print *, 'failed 2nd sub-test: ', Reverse_trim(str), Reverse_trim(substrs(nstars+1))
        goto 90 ! return
      endif
    endif
    ! lastSSindex = nstars
    if ( deebug ) print *, 'passed 2nd sub-test'
    if ( nstars < 2 ) goto 90 ! return
    ! Now find the indexes of these sub-patterns according to mode='wrap'
    istars(1:nstars-1) = indexes(str, substrs(2:nstars), mode='wrap')
    ! What we want to do is to check that no istars < 1
    relation = all ( istars(1:nstars-1) > 0 )
    ! Do the options instruct us to reverse the sense?
90  if ( reverseSense ) relation = .not. relation
  end function streq_scalar

  ! --------------------------------------------------  strings2Ints  -----
  subroutine strings2Ints (strs, ints)
    ! takes an array of strings and returns integer array
    ! using "ichar"
    ! Useful due to bug in toolbox swrfld
    !
    ! See also ints2Strings
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    character (len=*), intent(in), dimension(:) ::   strs
    integer, intent(out), dimension(:,:) ::          ints

    !----------Local vars----------!
    integer :: i, substr, strLen, arrSize
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

  end subroutine strings2Ints

   ! --------------------------------------------------  Swap  -----
  function Swap ( str, c1, c2 ) result ( outstr )
    ! Takes a string and returns one with c1 c2 swapped for each other
    ! E.g., to Swap every 'p' with a 'q' and every 'q' with a 'p'
    ! would transform 'perquisite' to 'qerpuisite'
    character(len=*), intent(in)  :: str
    character(len=1), intent(in)  :: c1
    character(len=1), intent(in)  :: c2
    character(len=len(str))       :: outstr
    ! Internal variables
    integer :: i, n
    ! Executable
    outstr = str
    if ( len(str) < 1 ) return
    do i=1, len(str)
      if ( str(i:i) == c1 ) then
        outstr(i:i) = c2
      elseif ( str(i:i) == c2 ) then
        outstr(i:i) = c1
      endif
    enddo
  end function Swap

  ! --------------------------------------------------  TrueList  -----
  subroutine TrueList ( Array, Str )
    ! Return Str with integers describing where Array has true values.
    ! Sequences of more than two consecutive true values are represented
    ! by n1:n2.
    logical, intent(in) :: Array(:)
    character(len=*), intent(out) :: Str
    integer :: L, N1, N2
    character :: Before
    before = ''
    l = 1
    n1 = 0
    str = ''
    do while ( n1 < size(array) )
      n1 = n1 + 1
      if ( .not. array(n1) ) cycle
      do n2 = n1, size(array)-1
        if ( .not. array(n2+1) ) exit
      end do
      ! n1 .. n2 are true
      if ( n1 == n2 ) then
        write ( str(l:), '(a,i0,:,":",i0)' ) trim(before), n1
      else
        write ( str(l:), '(a,i0,:,":",i0)' ) trim(before), n1, n2
      end if
      before = ','
      l = len_trim(str) + 1
      n1 = n2 + 1
    end do
  end subroutine TrueList
    
  ! -------------------------------------------------  unAsciify  -----
  function unAsciify ( str, how ) result ( outstr )
  ! takes input string and replaces any coded forms non-printing characters
  ! with corresponding binary characters outside ascii range
  ! leaving other chars alone
  ! This reverses the effect of the Asciify procedure which recoded
  ! non-ascii characters as '<something>'
  !
  ! How the decoding is done is according to the optional arg
  !    how         coding
  !    ---         -------
  !  'decimal'    <nnn> where nnn is the decimal value (e.g. 0)
  !  'octal'      <nnn> where nnn is the octal value (e.g. 000)
  !  'mnemonic'   <ID> where ID is the mnemonic code (e.g. NUL)
  ! Note that we may output a shorter string than the input

  ! method:
  ! (a) step through str, each stride positioned between
  ! j1, where leftmost lquote is, and j2, where leftmost rquote is
  ! (b) decode the current stride
  ! (c) accumulate evrything leftward plus decoded portion
  ! (d) cut off both stride and everything leftward. Call remainder "remainder"
  ! Contrast this approach with what is is done in MLSStringLists module
  ! which uses string slicing and index arithmetic
  ! Should we go back and simplify the latter?

  ! (see also ReplaceNonAscii, Asciify)
    !--------Argument--------!
    character (len=*), intent(in)           :: STR
    character (len=len(str))                :: OUTSTR
    character (len=*), optional, intent(in) :: HOW
    ! Internal variables
    character(len=1) :: bchar  ! A non-ascii character
    character(len=5) :: codedStr
    integer :: J1
    integer :: J2
    character(len=8) :: myHow
    character(len=len(str)) :: remaining
    character(len=1), parameter :: LQUOTE = '<'
    character(len=1), parameter :: RQUOTE = '>'
    ! Executable
    outstr=str
    myHow = 'mnemonic'
    if ( present(how) ) myHow = how
    ! print *, 'count_quotes ', count_quotes(str, LQUOTE, RQUOTE )
    if ( count_quotes(str, LQUOTE, RQUOTE ) < 1 ) return
    outstr = ' '
    remaining = str
    do
       j1 = index( remaining, lquote)
       if ( j1 < 1 ) exit
       j2 = index( remaining, rquote)
       if ( j2 < 1 ) exit
       codedStr = remaining(j1:j2)
       ! print *, 'j1, j2 ', j1, j2
       ! print *, 'codedStr ', codedStr
       bChar = decoder( codedStr, myHow )
       if ( j1 == 1 ) then
         outstr = trim(outstr) // bChar
       else
         outstr = trim(outstr) // remaining(:j1-1) // bChar
       endif
       if ( j2 < len_trim(remaining) ) then
         remaining = remaining(j2+1:)
       else
         remaining = ' '
         exit
       endif
    enddo
    if ( len_trim(remaining) > 1 ) &
      & outstr = trim(outstr) // remaining
    ! print *, 'unAsciify returning char-valued ', iachar(outstr(1:1))
  end function unAsciify

  ! --------------------------------------------------  unWrapLines  -----
  ! undo the splitting of commands across multiple lines by escaping line feeds
  ! i.e, a special escape character at line's end to denote a continuation
  ! to the following line
  ! optionally remove comment lines, i.e. lines beginning with comment character
  subroutine unWrapLines ( inLines, outLines, nOut, escape, comment )
    !
    character (len=*), dimension(:), intent(in)            ::   inLines
    character (len=*), dimension(:), intent(out)           ::   outLines
    integer, optional, intent(out)                         ::   nOut
    character (len=*), optional, intent(in)                ::   escape
    character (len=*), optional, intent(in)                ::   comment

    !----------Local vars----------!
    integer                      :: i, j, k, line
    character(len=len(outLines)) :: inLine
    character(len=1)             :: myComment
    character(len=1)             :: myEscape
    character(len=len(outLines)) :: outLine
    !----------Executable part----------!
    myEscape = '\'
    if ( present(escape) ) myEscape = escape
    myComment = achar(0)
    if ( present(comment) ) myComment = comment
    outlines = ' '
    line = 0
    k = 0
    do i=1, size(inLines)
      if ( line > size(outLines) - 1 ) exit
      inLine = adjustl(inLines(i))
      j = len_trim(inLine)
      if ( j < 1 ) then
        if ( k < 1 ) cycle
        ! We are done with the line
        line = line + 1
        outLines(line) = outLine
        k = 0
      else
        if ( inLine(1:1) == myComment ) cycle
        if ( inLine(j:j) == myEscape ) then
          if ( j < 2 ) then
            cycle
          elseif ( k < 1 ) then
            outline = inLine(:j-1)
            k = j-1
          else
            outLine = outLine(:k) // inLine(:j-1)
            k = k + j - 1
          endif
        elseif ( k < 1 ) then
          line = line + 1
          outLines(line) = inLine
        else
          line = line + 1
          outLines(line) = outLine(:k) // inLine
          k = 0
        endif
      endif
    enddo
    if ( present(nOut) ) nOut = line
  end subroutine unWrapLines

  ! --------------------------------------------------  writeAnIntToChars  -----
  ! takes an integer and returns a string
  ! using Fortran "write"
  ! Unless integer is one of specialInts, in which case
  ! we return corresponding one of specialChars
  ! (So that we can treat -1 as "unlimited' or -999 as 'FillValue')
  ! We'll just assume both special arrays are of same size
  ! and, in case of array versions of 
  subroutine writeAnIntToChars (int, str, fmt, specialInts, specialChars)
    !
    integer, intent(in)                                    ::   int
    character (len=*), intent(out)                         ::   str
    character (len=*), optional, intent(in)                ::   fmt
    integer, intent(in), dimension(:), optional            ::   specialInts
    character (len=*), intent(in), dimension(:), optional  ::   specialChars

    !----------Local vars----------!
    integer :: i
    character(len=16) :: MyStr
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
   else ! Don't use *; RTL can legally insert any number of leading blanks
     write ( myStr,'(i16)' ) int
     str = adjustl(myStr)
   endif

  end subroutine writeAnIntToChars

  subroutine writeIntsToChars_1d (ints, strs, fmt, specialInts, specialChars)
    ! takes an array of integers and returns string array
    integer, intent(in), dimension(:)            ::   ints
    character (len=*), intent(out), dimension(:) ::   strs
    character (len=*), intent(in), optional      ::   fmt
    integer, intent(in), dimension(:), optional  ::   specialInts
    character (len=*), intent(in), dimension(:), optional  ::   specialChars

    !----------Local vars----------!
    integer :: i, arrSize
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

  end subroutine writeIntsToChars_1d

  subroutine writeIntsToChars_2d (ints, strs, fmt, specialInts, specialChars)
    integer, intent(in), dimension(:,:)            ::   ints
    character (len=*), intent(out), dimension(:,:) ::   strs
    character (len=*), intent(in), optional      ::   fmt
    integer, intent(in), dimension(:), optional  ::   specialInts
    character (len=*), intent(in), dimension(:), optional  ::   specialChars

    !----------Local vars----------!
    integer :: i, j, shp(2)
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   shp = shape(ints)
   if ( any(shp <= 0) ) then
     strs = ' '
     return
   endif
   do j=1, shp(2)
     do i=1, shp(1)
       call writeAnIntToChars(ints(i,j), strs(i,j), &
         & fmt, specialInts, specialChars)
     enddo
   enddo

  end subroutine writeIntsToChars_2d

  ! ----------------  writeIntAsBaseN  -----
  ! write an integer reexpressed in base n as a string
  
  ! Use:
  ! Given e.g. the number 1989 in base 5 return str='30424'
  ! i.e., 4 + 5*(2 + 5*(4 + 5*(0 + 5*3)))
  subroutine writeIntAsBaseN ( K, N, STR, OPTIONS )
    ! takes an integer and returns a string
    ! using poorly-tested but hopefully non-critical code

    ! Method:
    ! First determine how many terms in expansion, or equivalently,
    ! highest power of N; thus want N^m < num
    ! m < (log k) / (log N)
    ! For its coefficient c[m], must have
    ! c[m] < k / N^m
    ! FFor succeeding coefficients
    ! c[m-1] < k / N^(m-1) - N c[m]
    ! c[m-2] < k / N^(m-2) - N^2 c[m] - N c[m-1]
    ! ...    ...    ...    ...    ...
    ! c[m-j] < k / N^(m-j) - N^j c[m] - N^(j-1) c[m-1] - .. - N c[m-j+1]
    ! until for j up to m
    
    ! According to the options string, the resulting integer coefficients
    ! will be returned as
    !  options contain            coefficient string
    !    ---                      ------------------
    !   (default)                 decimal integers, separated by spaces
    !                              e.g., ''13 4 29 6''
    !      a                      achar(c[m]) in the same order
    !                              e.g., CR // EOT // GS // ACR   
    !      A                      like a but converted to ascii
    !                              e.g., '<CR>' // '<EOT>' // '<GS>' // <ACR>   
    !      x                      "Extended" hexadecimal
    ! 0123456789012345678901234567890123456789012345678901234567890123456789
    ! 0123456789abcdefghijklmnopqrstuvwxyz;'[]ABCDEFGHIJKLMNOPQRSTUVWXYZ:"{}
    !                              (good only up to base 70)
    !
    ! The default representation can handle arbitrarily large moduli
    ! a or A is limited to the unique range of achar
    ! x is good only up to a modulus of 70.
    
    ! Note that the returned string is in normal order. The coefficient of the
    ! highest power comes first, then the next highest power, and
    ! so on. If you want the coeffieicents in reversed order, call
    ! this module's Reverse function afterward.
    
    ! If you want to remove the embedded spaces from one of the non-default
    ! representations, call CompressString with the result
    
    ! If you want to convert a floating number, say '4.523' to another base,
    ! (i) multiply it by the base as many times as you will want after
    !     that base's "decimal" point. So for base 13 to 3 places, form
    !    13*13*13*4.523
    ! (ii) takes the integer part of the result
    ! (iii) call this subroutine with that integer
    ! (iv) Insert the "decimal" point 3places before the result's last char
    !
    !--------Argument--------!
    integer, intent(in)                     ::   k ! integer to be reexpressed
    integer, intent(in)                     ::   N ! the new base
    character (len=*), intent(out)          ::   str
    character(len=*), optional, intent(in)  ::   options
    ! Internal variables
    integer, dimension(0:128)               :: C ! One for each digit of K
    integer                                 :: j
    integer                                 :: m
    character(len=8)                        :: myOptions
    double precision                        :: S, SHigher
    character(len=8)                        :: substr
    character(len=*), parameter             :: xtended = &
      & '0123456789abcdefghijklmnopqrstuvwxyz;''[]ABCDEFGHIJKLMNOPQRSTUVWXYZ:"{}'
    ! Executable
    str = ' '
    if ( k < 0 .or. N < 2 ) return
    myOptions = ' '
    if ( present(options) ) myOptions = options
    ! What is the highest power needed of N?
    m = 0
    if ( k > 1 ) m = ( log (k*1.d0) / log (N*1.d0) ) + 1.d-3
    if ( k < N ) then
      C(0) = k
    else
      ! Recursively define
      ! S[m]     = k / N^m
      ! S[m-1] = N ( S[m] - c[m] )
      ! S[m-j] = N ( S[m-j+1] - c[m-j+1] )
      ! but we'll restate this for each j as
      ! S = N ( SHigher - c[m-j=1] )
      SHigher = k*1.d0 / N**m
      j = 0
      do
        if ( (m - j) < 0 ) exit
        c(m-j) = SHigher + 1.d-3
        S = N*( SHigher - c(m-j) )
        SHigher = S
        j = j + 1
      end do
    end if
    ! print *, 'm ', m
    ! call dump( c, 'coeffs' )
    ! Now we turn the integers into a character string
    ! in an appropriate way
    if ( index(myOptions, 'a') > 0 .or. index(myOptions, 'A') > 0 ) then
      do j = m, 0, -1
        if ( c(j) < 177 ) str = trim(str) // ' ' // achar(c(j))
        ! print *, 'str ', Asciify(str)
      end do
      str = adjustl(str)
      if ( index(myOptions, 'A') > 0 ) str = Asciify(str, 'mnemonic')
    else if ( index(myOptions, 'x') > 0 ) then
      do j = m, 0, -1
        if ( c(j) < 71 ) str = trim(str) // ' ' // xtended(c(j)+1:c(j)+1)
        ! print *, 'str ', str
      end do
      str = adjustl(str)
    else
      ! default
      do j = m, 0, -1
        call writeIntsToChars ( c(j), substr )
        str = trim(str) // ' ' // adjustl(substr)
        ! print *, 'substr ', substr
        ! print *, 'str ', str
      end do
      str = adjustl(str)
    end if
  end subroutine writeIntAsBaseN

  subroutine writeRomanNumerals ( num, str, options )
    ! takes an integer and returns a string
    ! using poorly-tested but hopefully non-critical code

    ! Method:
    ! Peel off values from num, working from highest to lowest,
    ! each time appending the corresponding romans to the end of str
    
    ! The optional string options may control the following
    ! if options contains  meaning
    ! 4     use 'iiii' for 4 instead of 'iv'; seen in some clock faces
    ! n     non-standard

    ! If non_standard is chosen, check for exceptional abbreviations
    ! like 'im' for 999, or 'cic' for 199
    !
    ! Note: The '4' option mnioned above doesn't seem tto be operative.
    ! Did you change your mind?
    !--------Argument--------!
    integer, intent(in)                     ::   num
    character (len=*), intent(out)          ::   str
    character(len=*), optional, intent(in)  ::   options
    ! Internal variables
    character(len=2), dimension(13), parameter :: romans = &
      & (/  'm ','cm','d ','cd', 'c ','xc','l ','xl','x ','ix','v ','iv','i ' /)
    integer, dimension(13), parameter :: values = &
      & (/ 1000, 900,  500, 400,  100, 90, 50,  40,  10,   9,   5,   4,   1 /)
    integer :: remainder
    integer :: r_index ! What romans index we're searching for
    logical :: non_standard
    ! Executable
    str = ' '
    if ( len(str) < 1 .or. num < 1 ) return
    non_standard = .false.
    if ( present(options) ) non_standard = index(lowercase(options), 'n') > 0
    remainder = num
    ! Outer loop: romans index
    do r_index=1, 13
      ! Will adding an 'i' or 'x' boost the remainder above values(r_index)?
      if ( non_standard .and. remainder < values(r_index) ) then
        select case (r_index)
        case (1,3,5,7) ! m, d, c, l
          if ( remainder + values(13) >= values(r_index) ) then
            str = trim(str) // romans(13)
            remainder = remainder + values(13)
          elseif( remainder + values(9) >= values(r_index) ) then
            str = trim(str) // romans(9)
            remainder = remainder + values(9)
          endif
        end select
      endif
      do while ( remainder >= values(r_index) )
        str = trim(str) // romans(r_index)
        remainder = remainder - values(r_index)
      enddo
    enddo
    ! print *, 'Remainder: ', remainder
  end subroutine writeRomanNumerals

  ! Private procedures and functions
!============================ Private ==============================
  subroutine prepOptions( options )
    ! Process options into separate optional args
    ! You should call this at the start of every procedure
    ! that uses options to set countEmpty, etc.
    ! Args:
    character(len=*), intent(in), optional  :: options
    ! Internal variables
    character(len=16) :: myOptions
    ! Executable
    caseSensitive       = .true.
    ignoreLeadingSpaces = .false.
    myOptions = STRINGOPTIONS
    if ( present(options) ) myOptions = options
    if ( len_trim(myOptions) > 0 ) then
      caseSensitive       = ( index(myOptions, 'c') < 1 )
      ignoreLeadingSpaces = ( index(myOptions, 'f') > 0 )
    endif
  end subroutine prepOptions

  ! ---------------------------------------------------  isComment  -----
  elemental function isComment( arg, Comment ) result( itIs )
    ! Returns TRUE if first non-blank char is comment
    ! Args
    character(len=*), intent(in)    :: arg
    character, optional, intent(in) :: Comment
    logical                         :: itIs
    ! Internal variables
    character :: first
    character :: myComment
    ! Executable
    myComment = '#'
    if ( present(Comment) ) myComment = Comment
    first = adjustl(arg)
    itIs = ( first == myComment )
  end function isComment

  ! ---------------------------------------------------  isDigit  -----
  elemental function isDigit(arg) result(itIs)
    ! Returns TRUE if arg is one of {'1', '2', ..}
    ! Args
    character(len=1), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    character(len=*), parameter :: list='1234567890'
    ! Executable
    itIs = ( index(list, arg) > 0 )
  end function isDigit

  ! ---------------------------------------------------  isDigits  -----
  elemental function isDigits(arg) result(itIs)
    ! Returns TRUE if arg is made up of {'1', '2', ..} and trailing space(s)
    ! Args
    character(len=*), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    integer :: i
    ! Executable
    itIs = .false.
    if ( len_trim(arg) < 1 ) return
    do i=1, len_trim(arg)
      if ( .not. isDigit(arg(i:i)) ) return
    enddo
    itIs = .true.
  end function isDigits

  ! ---------------------------------------------------  lastChar  -----
  character function lastChar(str)
    character(len=*), intent(in) :: str
     ! Returns the last non-blank character of str (unless str itself is blank)
     integer :: strlen
     lastchar = ' '
     strlen = len_trim(str)
     if ( strlen < 1 ) return
     lastchar = str(strlen:strlen)
  end function lastChar

  ! ---------------------------------------------------  decoder  -----
  function decoder( said, how ) result ( meant )
    ! decode what you said, a coded ascii string like this
    ! '<NUL>'
    ! into what you actually meant, i.e. a null char
    character(len=*), intent(in)           :: said
    character(len=*), intent(in)           :: how
    character(len=1)                       :: meant
    ! Executable
    ! print *, 'said ', said
    ! print *, 'how ', how
    if ( How == 'decimal' ) then
      meant = decimalCode( said, invert=.true. )
    elseif ( How == 'octal' ) then
      meant = octalCode( said, invert=.true. )
    elseif ( How == 'mnemonic' ) then
      meant = mnemonicCode( said, invert=.true. )
    else
      ! Uh-oh, that is not a recognized how
      ! print *, 'Uh-oh, that is not a recognized how'
      meant = achar(0) ! return a NUL char; should we crash, instead?
    endif
  end function decoder

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

  subroutine readAnIntFromChars (str, num, forbiddens, ignore)
    !--------Argument--------!
    character (len=*), intent(in)               ::   str
    integer, intent(out)                        ::   num
    character (len=*), intent(in), optional     ::   forbiddens
    character (len=*), optional, intent(in)     ::   ignore
    call ReadIntsFromChars (str, num, forbiddens, ignore)
  end subroutine readAnIntFromChars

  subroutine readIntArrayFromChars (strs, ints, forbiddens, ignore)
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    character (len=*), intent(in), dimension(:) ::   strs
    integer, intent(out), dimension(:)          ::   ints
    character (len=*), intent(in), optional     ::   forbiddens
    character (len=*), optional, intent(in)     ::   ignore

    !----------Local vars----------!
    integer :: i, arrSize
    !----------Executable part----------!
    call ReadIntsFromChars (strs, ints, forbiddens, ignore)
  end subroutine readIntArrayFromChars

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSStrings
!=============================================================================

! $Log$
! Revision 2.112  2021/07/09 21:48:27  pwagner
! Clearer comments in CompressString
!
! Revision 2.111  2019/10/30 20:04:31  pwagner
! Added Justify function
!
! Revision 2.110  2019/10/16 20:50:13  pwagner
! Added Enclosure
!
! Revision 2.109  2019/07/17 20:13:42  pwagner
! Replace may now take the arg reverse to replace non-matching oldChar; added Swap
!
! Revision 2.108  2019/04/09 20:35:17  pwagner
! Moved some procedures from MLSStrings to new MLSStrings_0
!
! Revision 2.107  2019/03/21 23:44:38  pwagner
! Fixed error in nCharsinFormat; cosmetic improvements
!
! Revision 2.106  2019/01/24 18:31:12  pwagner
! Reorganized modules that print to simplify toolkit-free builds
!
! Revision 2.105  2018/03/27 22:06:52  pwagner
! No longer uses optional_m; now one less module needed to build toolkitlessly
!
! Revision 2.104  2017/12/08 00:28:44  pwagner
! Add optional 'options' arg to CatStrings; by default will not cat empty string
!
! Revision 2.103  2017/01/25 21:10:17  pwagner
! str1 and str2 args to streq now optional
!
! Revision 2.102  2016/12/08 00:15:41  pwagner
! readNumsFromChars can now ignore non-numerical stuff
!
! Revision 2.101  2014/09/11 18:23:57  pwagner
! Added isDigits
!
! Revision 2.100  2014/09/05 00:12:11  vsnyder
! Convert pointer temp to explicit shape.
!
! Revision 2.99  2014/07/25 21:42:37  pwagner
! Fixed bugs in readIntFromBaseN; now generic as readNumFromBaseN; added remap
!
! Revision 2.98  2014/02/21 01:24:54  pwagner
! Added CharToInt; made SplitWords elemental
!
! Revision 2.97  2014/02/12 23:59:34  pwagner
! isAlphabet now public
!
! Revision 2.96  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.95  2013/08/09 00:43:06  pwagner
! Added count_quotes and unAsciify
!
! Revision 2.94  2013/07/30 23:26:44  pwagner
! New procedures to read/write base N integers
!
! Revision 2.93  2013/07/24 19:02:43  pwagner
! Added isComment function; streq adds new option 'h' to match string heads
!
! Revision 2.92  2013/06/28 23:56:47  pwagner
! Fixed an old bug in SplitDetails
!
! Revision 2.91  2013/06/17 21:35:26  pwagner
! You may specify better where stretch adds space
!
! Revision 2.90  2013/06/14 03:02:37  vsnyder
! Simplify ReplaceNonAscii
!
! Revision 2.89  2013/06/13 00:39:21  pwagner
! Removed two unused params
!
! Revision 2.88  2012/12/04 00:11:55  pwagner
! Improved comments
!
! Revision 2.87  2012/09/05 21:43:55  pwagner
! streq options now case sensitive; 'P'artial match means str1(i) is the part
!
! Revision 2.86  2012/08/07 18:02:37  pwagner
! ReplaceNonAscii now takes optional arg exceptions which dont get replaced
!
! Revision 2.85  2012/07/10 15:15:33  pwagner
! Added SplitDetails
!
! Revision 2.84  2012/05/01 22:10:26  vsnyder
! Add TrueList subroutine
!
! Revision 2.83  2011/06/23 17:25:50  pwagner
! Added ability to read, write Roman numerals
!
! Revision 2.82  2011/06/16 00:14:51  pwagner
! Added new procedures to unwrap commands and Start Case strings
!
! Revision 2.81  2011/02/18 17:58:10  pwagner
! Replace, Delete no take an optional arg; added shiftLRC
!
! Revision 2.80  2011/02/05 01:34:23  pwagner
! Added Delete function
!
! Revision 2.79  2010/09/24 23:45:45  pwagner
! Removed dependence on MLSCommon by substituting double for r8
!
! Revision 2.78  2010/06/23 20:42:21  honghanh
! Change readAnIntFromChars to only read the first number in the string
! if there are many numbers separated by ignored characters.
! Update comment of the method.
!
! Revision 2.77  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.76  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.75  2009/06/16 17:08:42  pwagner
! Added ReadNumsFromChars, squeeze, stretch functions
!
! Revision 2.74  2008/06/04 21:44:43  pwagner
! Added lenTrimToAscii and ReplaceNonAscii
!
! Revision 2.73  2008/05/09 00:22:37  pwagner
! Nappearances has new optional arg to not trim substrings
!
! Revision 2.72  2008/02/07 18:46:55  pwagner
! isAllAscii, isAscii now public; Asciify generic
!
! Revision 2.71  2007/09/13 21:06:25  pwagner
! Added 2-d array interface for writeIntsToChars
!
! Revision 2.70  2007/08/29 19:52:18  pwagner
! Added Asciify function
!
! Revision 2.69  2007/07/31 22:46:08  pwagner
! undefined status now defined in readAnIntFromChars ;'n' option added to streq
!
! Revision 2.68  2007/07/25 21:58:16  vsnyder
! Replace tabs by spaces because tabs are not standard
!
! Revision 2.67  2007/07/18 00:06:46  pwagner
! Added Rot13 function
!
! Revision 2.66  2007/05/22 20:57:18  vsnyder
! Don't use the length of one automatic variable to specify the length of
! another one: Intel ifort 10.0.023 crashes at run time on this.
! Don't use list-directed output to internal files.
!
! Revision 2.65  2007/04/26 20:32:15  pwagner
! Coded around bug in Intel compiler causing streq to bomb
!
! Revision 2.64  2007/01/03 20:40:25  pwagner
! Added NAppearances
!
! Revision 2.63  2006/10/05 23:34:13  pwagner
! Fixed bug in streq making identity comparisons
!
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
