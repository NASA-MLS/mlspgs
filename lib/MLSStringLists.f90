! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSStringLists               ! Module to treat string lists
!=============================================================================

  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
    & MLSMSG_Allocate, MLSMSG_DeAllocate
  use MLSCommon, only: i4, r8, NameLen, BareFNLen
  use MLSSets, only: FindFirst
  use MLSStrings, only: lowerCase, Capitalize, reverse

  implicit NONE
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

!     (parameters)
! KEYNOTFOUND        key not found among keyList
! KEYBEYONDHASHSIZE  index of key in keyList > size(hash array)
! INVALIDHHMMSSSTRING
!                    (if strict) str not in format '00:00:00.0000000'
! LENORSIZETOOSMALL  Either charsize of strs or size(ints) too small

!     (subroutines and functions)
! Array2List         Converts an array of strings to a single string list
! catLists           cats 2 string lists, taking care if either one is blank
! ExpandStringRange  Turns '1,2-5,7' into '1,2,3,4,5,7'
! ExtractSubString   Extracts portion of string sandwiched between sub1 and sub2
! GetIntHashElement  Returns int from hash array corresponding to key string
! GetStringElement   Returns n'th element of string list
! GetStringHashElement   
!                    Returns string from hash list corresponding to key string
! GetUniqueList      Returns str list of only unique entries from input list
! GetUniqueStrings   Returns array of only unique entries from input array
! hhmmss_value       Converts 'hh:mm:ss' formatted string to a real r8
!                    (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
! List2Array         Converts a single string list to an array of strings
! NumStringElements  Returns number of elements in string list
! RemoveElemFromList removes occurrence(s) of elem in a string list
! ReplaceSubString   replaces occurrence(s) of sub1 with sub2 in a string
! ReverseList        Turns 'abc,def,ghi' -> 'ghi,def,abc'
! SortArray          Turns (/'def','ghi','abc'/) -> (/'abc','def','ghi'/)
! SortList           Turns 'def,ghi,abc' -> 'abc,def,ghi'
! StringElementNum   Returns element number of test string in string list
! unquote            Removes surrounding [quotes]
! utc_to_yyyymmdd    Parses yyyy-mm-ddThh:mm:ss.sss or yyyy-dddThh:mm:ss.sss
! yyyymmdd_to_dai    Converts yyyymmdd to days after Jan 1, 2001
! === (end of toc) ===

! === (start of api) ===
! Array2List (char* inArray(:), char* outList(:), &
!   & [char inseparator], [int ordering], [char leftRight]) 
! char* catLists (char* str1, char* str2)
! ExpandStringRange (char* str, char* outst)
! ExtractSubString (char* str, char* outstr, char* sub1, char* sub2, &
!       & [char* how], [log no_trim])
! GetIntHashElement (strlist keyList, hashArray(:), char* key, 
!   int ErrType, log countEmpty, [char inseparator], [log part_match])
! GetStringElement (strlist inList, char* outElement,
!   i4 nElement, log countEmpty, [char inseparator])
! GetStringHashElement (strlist keyList, strlist hashList, char* key, 
!   char* outElement, log countEmpty, [char inseparator], [log part_match])
! GetUniqueList (char* str, char* outstr(:), int noUnique, &
!   & log countEmpty, [char inseparator], [log IgnoreLeadingSpaces]) 
! GetUniqueStrings (char* inList(:), char* outList(:), int noUnique) 
! r8 hhmmss_value (char* str, int ErrTyp, [char separator], [log strict])
! List2Array (strlist inList, char* outArray(:), log countEmpty, [char inseparator],
!    [log IgnoreLeadingSpaces])
! int NumStringElements(strlist inList, log countEmpty, &
!   & [char inseparator], [int LongestLen])
! ReplaceSubString (char* str, char* outstr, char* sub1, char* sub2, &
!       & [char* which], [log no_trim])
! strlist ReverseList (strlist str, [char inseparator])
! SortArray (char* inStrArray(:), int outIntArray(:), log CaseSensitive, &
!   & [char* sortedArray(:)], [log shorterFirst], [char leftRight])
! SortList (strlist inStrArray, int outIntArray(:), log CaseSensitive, &
!   & log countEmpty, [char inseparator], [log IgnoreLeadingSpaces], 
!     [strlist sortedList], [char leftRight])
! int StringElementNum(strlist inList, char* test_string, log countEmpty, &
!    & [char inseparator], [log part_match])
! char* unquote (char* str, [char* quotes], [char* cquotes], [log strict])

! in the above, a string list is a string of elements (usu. comma-separated)
! e.g., units='cm,m,in,ft'
! an array is a Fortran array of strings or integers
! a hash is a list of key strings and either
! (1) a list of associated strings
! (2) an array of associated integers
! (an idea called a hash in perl or a dictionary in python)
! Many of these routines take optional arguments that greatly modify
! their default operation

! Warnings: 
! (1) in the routines Array2List, and SortArray
! the input arguments include an array of strings;
! This array is of assumed-size
! I.e., all elements from array(1:size(array)) are relevant
! Therefore in calling one of these you probably need to use the format
!   call SortArray(myArray(1:mySize), ..
! to avoid operating on undefined array elements
! (2) In operating on string lists it is sometimes assumed that no
! element is longer than a limit: MAXSTRELEMENTLENGTH
! === (end of api) ===

  public :: catLists, Array2List, ExpandStringRange, ExtractSubString, &
   & GetIntHashElement, GetStringElement, GetStringHashElement, &
   & GetUniqueStrings, GetUniqueList, &
   & hhmmss_value, &
   & List2Array, NumStringElements, &
   & RemoveElemFromList, ReplaceSubString, ReverseList, &
   & SortArray, SortList, StringElementNum, &
   & unquote, utc_to_yyyymmdd, yyyymmdd_to_dai

  interface switch
    module procedure switch_ints
  end interface

  interface utc_to_yyyymmdd
    module procedure utc_to_yyyymmdd_strs, utc_to_yyyymmdd_ints
  end interface

  interface yyyymmdd_to_dai
    module procedure yyyymmdd_to_dai_str, yyyymmdd_to_dai_ints
  end interface

  ! Error return values from:
  ! GetIntHashElement
  integer, public, parameter :: KEYNOTFOUND=-1
  integer, public, parameter :: KEYBEYONDHASHSIZE=KEYNOTFOUND-1
  ! hhmmss_value
  integer, public, parameter :: INVALIDHHMMSSSTRING = 1
  ! utc_to_yyyymmdd
  integer, public, parameter :: INVALIDUTCSTRING = 1
  ! strings2Ints
  integer, public, parameter :: LENORSIZETOOSMALL=-999
  
  ! A limitation among string list operations
  integer, private, parameter :: MAXSTRELEMENTLENGTH = BareFNLen

  integer, parameter :: YEARMAX = 4999  ! Conversion invalid after 4999 AD
  ! The following arrys contains the maximum permissible day for each month
  ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
  integer, dimension(-1:12), parameter :: DAYMAXLY = (/ &
    & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)
  integer, dimension(-1:12), parameter :: DAYMAXNY = (/ &
    & 365, 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)
CONTAINS

  ! ---------------------------------------------  Array2List  -----

  ! This subroutine returns a (usually) comma-separated string list, interpreted it
  ! as a list of individual elements, given an equivalent array of
  ! sub-strings in which the n'th element becomes the n'th element

  ! As an optional arg the separator may supplied, in case it isn't a comma
  ! As an optional arg the ordering in which the array elements are to be
  ! taken may be supplied; e.g. (/4, 1, 3, 2/) means 1st take 4th element,
  ! then 1st, then 3rd, and finally 2nd: list[k] = array[ordering[k]]
  ! (unless the further optional arg leftRight is also supplied and equals
  ! one of {"l", "L"} in which case list[ordering[k]] = array[k])

  SUBROUTINE Array2List ( inArray, outList, inseparator, ordering, leftRight )
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(OUT)                :: outList
    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN)   :: inArray
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: inseparator
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)   :: ordering
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: leftRight

    ! Local variables
    INTEGER(i4) :: listElem, arrayElem, nElems

    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: BLANK = ' '   ! Returned for any element empty
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    CHARACTER (LEN=1)               :: myLeftRight
    ! Executable code

    IF(PRESENT(inseparator)) THEN
      separator = inseparator
    ELSE
      separator = COMMA
    END IF

    IF(PRESENT(leftRight)) THEN
      myleftRight = Capitalize(leftRight)
    ELSE
      myleftRight = "R"
    ENDIF

    if ( len(outList) <= 0 ) return
    outList = BLANK
    nElems = size(inArray)
    if ( nElems <= 0 ) return
	 listElem = 1
    DO
      if (.not. present(ordering) ) then
        arrayElem = ListElem
      elseif (myLeftRight == "R") then
        arrayElem = ordering(ListElem)
      else
        ! Try to invert ordering function
        do arrayElem=1, nElems
          if ( ordering(arrayElem) == listElem ) exit
        enddo
        arrayElem = min(arrayElem, nElems)
      endif
      if ( listElem == 1 ) then
        outList = trim(inArray(arrayElem))
      else
        outList = trim(outList) // separator // trim(inArray(arrayElem))
      endif
      listElem = listElem + 1
      if ( listElem > min(nElems, len(outList)) ) return
    ENDDO

  END SUBROUTINE Array2List

  ! -------------------------------------------------  catLists  -----
  function catLists (STR1, STR2, inseparator) result (OUTSTR)
    ! cats 2 string lists, taking care if either is blank
    ! E.g., given str1 = 'a,b,c' and str2 = 'd,e,f'
    ! returns 'a,b,c,d,e,f'
    ! If either is blank, returns the other
    ! If both blank, returns a blank
    !--------Argument--------!
    character (len=*), intent(in) :: STR1
    character (len=*), intent(in) :: STR2
    character (len=1), optional, intent(in)       :: inseparator
    character (len=len(str1)+len(str2)+1) :: OUTSTR

    !----------Local vars----------!
    character (len=1), parameter    :: COMMA = ','
    character (len=1)               :: separator
    !----------executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = comma
    end if
    if ( len_trim(str2) < 1 ) then
      outstr=str1
    elseif ( len_trim(str1) < 1 ) then
      outstr=str2
    else
      outstr = trim(str1) // separator // trim(str2)
    endif
  end function catLists

  ! --------------------------------------------------  ExpandStringRange  -----
  subroutine ExpandStringRange (instr, outstr)
    ! Takes a string list and expands any ranges found within:
    ! Ranges are marked by patterns
    ! (simple) 'n-m' integers from n to m inclusive, where m > n
    ! (stride) 'n-m+s' integers from n to m with stride s, where s > 0
    ! Examples:
    ! '1-10' becomes '1,2,3,4,5,6,7,8,9,10'
    ! '2-21+2' bcecomes '2,4,6,8,10,12,14,16,18,20' (missing '21')

    ! Errors and limitations:
    ! '5-5' becomes simply '5' (range is inclusive, but not duplicative)
    ! '6-4' becomes '' (we can't go backward)
    ! m, n must be non-negative
    ! The separator between elements must be ','
    ! Although countEmpty set to .true. below, actual behavior
    ! ignores blank elements; e.g. 
    ! '0,1,,2,,4-6,,' becomes '0,1,2,4,5,6'
    ! Args
    character (len=*), intent(in) :: instr
    character (len=*), intent(inout) :: outstr
    ! Internal variables
    character (len=len(outstr)) :: tempstr
    integer :: dashpos
    integer :: elem
    integer :: ErrTyp
    integer :: m
    character (len=8) :: mChar
    integer :: n
    character (len=8) :: nChar
    integer :: nelem
    integer :: pluspos
    integer :: s
    character (len=8) :: sChar
    integer :: t
    character (len=8) :: tChar
    logical, parameter :: countEmpty=.true.
    ! Executable
    outstr = instr
    nelem = NumStringElements(instr, countEmpty)
    dashpos = index(instr, '-')
    if ( nelem < 1 .or. dashpos < 1 ) return
    outstr = ' '
    do elem = 1, nelem
      call GetStringElement (instr, tempstr, elem, countEmpty)
      dashpos = index(trim(tempstr), '-')
      if ( dashpos > 0 ) then
        n = -999
        m = -999
        s = 1
        call GetStringElement (trim(tempstr), nChar, 1, countEmpty, &
          & inseparator='-')
        pluspos = index(trim(tempstr), '+')
        if ( pluspos > 0 ) then
          call ExtractSubString (trim(tempstr), mChar, '-', '+')
          call GetStringElement (trim(tempstr), sChar, 2, countEmpty, &
            & inseparator='+')
          ! non-simple pattern is 'n-m+s'
        else
          ! simple pattern is 'n-m'
          call GetStringElement (trim(tempstr), mChar, 2, countEmpty, &
            & inseparator='-')
          sChar = ''
        endif
        read(nChar, *, iostat=ErrTyp) n
        if ( ErrTyp == 0 ) read(mChar, *, iostat=ErrTyp) m
        if ( ErrTyp == 0 .and. sChar /= '') read(sChar, *, iostat=ErrTyp) s
        ! print *, 'nChar: ', trim(nChar)
        ! print *, 'mChar: ', trim(mChar)
        ! print *, 'sChar: ', trim(sChar)
        ! print *, 'n, m, s: ', n, m, s
        tempstr = ''
        if ( m >= n ) then
          do t=n, m, s
            write(tChar, *) t
            ! print *, 'tChar: ', trim(tChar)
            tempstr = catLists(trim(tempstr), adjustl(tChar))
          enddo
        endif
      endif
      outstr = catLists(trim(outstr), adjustl(tempstr))
    enddo
  end subroutine ExpandStringRange

  ! --------------------------------------------------  ExtractSubString  -----
  SUBROUTINE ExtractSubString (instr, outstr, sub1, sub2, how, no_trim)
    ! Takes a string and extracts what is sandwiched between sub1 and sub2
	 ! Defaults to choosing only the first occurrence of sub1 and sub2
    ! But if how == 'greedy' chooses last occurrence of sub2
    ! or if how == 'stingy' chooses last occurrence of sub1
    ! Note that, depending on how, we extract:
    !    (let sub1='abc' sub2='def' str='abcabc123defdef')
    ! (a) if how == default => 'abc123'
    ! (b) if how == greedy => 'abc123def'
    ! (c) if how == stingy => '123'
    ! If no_trim is TRUE, sub1 and sub2 may have trailing spaces
    ! that will not be trimmed before attempting to match
    ! Method:
    ! Replace substrings sub1 and sub2 with separator character
    ! and then use GetStringElement to get subelement number 2
    !  
    ! A fundamental issue arises if sub2 occurs before sub1 in the string
    ! Do we want to interpret the request such that we
    ! (1) return a blank
    ! (2) look for occurrences of sub2 in the string after sub1
    ! I think we should aim for 2, as it produces a generalization
    ! of picking elements out of a comma-separated list
    !
    ! Will this still work if sub1 has leading or trailing blanks? 
    ! How about sub2?
    ! Do we need an optional arg, no_trim, say, that will leave them?
    ! Tried coding it, but can't say for sure it works
    ! What if sub1 is a substring of sub2, or vice versa?
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: instr
    CHARACTER (LEN=*), INTENT(IN) :: sub1
    CHARACTER (LEN=*), INTENT(IN) :: sub2
    CHARACTER (LEN=*), INTENT(INOUT) :: outstr
    character (len=*), intent(in), optional :: how
    logical, intent(in), optional :: no_trim

    !----------Local vars----------!
    CHARACTER (LEN=len(instr)) :: str
    integer, parameter         :: MAXREPLACEMENTS = 100
    integer, parameter         :: EARLYSUB2INTERPRETATION = 2 ! 2 or 1
    INTEGER :: i, isub1, isub2, strlen, tmpstrlen
    character (len=7) :: my_how
    character(len=1) :: separator
    character(len=*), parameter :: separators =',.$%#{}()'
    character (len=max(len(str), len(outstr))) :: tmpstr
    character (len=max(len(str), len(outstr))) :: tmpstr2
    logical :: my_no_trim, trimming
    !----------Executable part----------!
    my_how = 'first'
    if ( present(how) ) my_how = lowercase(how)
    my_no_trim = .false.
    if ( present(no_trim) ) my_no_trim = no_trim
    trimming = .not. my_no_trim
    outstr = ' '
    strlen = LEN_trim(instr)
    if (strlen < 1 .or. instr == ' ') return
    ! print *, 'instr: ', trim(instr)
    ! print *, 'strlen: ', strlen
    ! Which interpretation of sub2 occurring before sub1 do we make?
    isub1 = index(instr, trim(sub1))
    isub2 = index(instr, trim(sub2))
    ! print *, 'sub1: ', trim(sub1)
    ! print *, 'isub1: ', isub1
    ! print *, 'sub2: ', trim(sub2)
    ! print *, 'isub2: ', isub2
    if ( isub2 < isub1 ) then
      if ( isub2 == 0 ) then
        return
      elseif ( EARLYSUB2INTERPRETATION == 2 ) then
        ! zap every occurrence of sub2 up to position isub1
        ! print *, 'zap every occurrence of sub2 up to position isub1'
        ! print *, 'before zapping: ', instr(1:isub1-1)
        ! print *, 'tail: ', instr(isub1:strlen)
        call ReplaceSubString (instr(1:isub1-1), tmpstr, sub2, '', &
          & which='all', no_trim=.false.)
        ! print *, 'afterzapping: ', trim(tmpstr), '//', instr(isub1:strlen)
        tmpstrlen = len_trim(tmpstr)
        ! print *, 'tmpstrlen: ', tmpstrlen
        tmpstrlen = len(trim(tmpstr))
        ! print *, 'tmpstrlen(2): ', tmpstrlen
        str = ''
        if ( tmpstrlen < 1 ) then
          str = instr(isub1:strlen)
        else
          str = tmpstr(1:tmpstrlen) // instr(isub1:strlen)
          ! str = tmpstr(1:tmpstrlen)
          ! print *, tmpstr(1:tmpstrlen)
          ! print *, str(1:tmpstrlen)
          ! str(tmpstrlen+1:) = instr(isub1:strlen)
          ! print *, instr(isub1:strlen)
          ! print *, trim(str(tmpstrlen+1:))
        endif
      else
        str = instr
      endif
    else
      str = instr
    endif
    ! print *, 'str: ', trim(str)
    if ( trimming ) then
      if (len_trim(sub1) < 1 &
        & .or. &
        & len_trim(sub2) < 1 .or. index(str, trim(sub1)) == 0 &
        & .or. &
        & index(str, trim(sub2)) == 0 ) RETURN
    else
      if (index(str, sub1) == 0 &
        & .or. &
        & index(str, sub2) == 0 ) RETURN
    endif
    do i=1, len(separators)
      if ( index(str, separators(i:i)) == 0 ) exit
    enddo
    if ( i > len(separators) ) return   ! This means our method will fail
    separator = separators(i:i)
    ! print *, 'separator: ', separator
    select case (trim(my_how))
    case ('greedy')
      if ( trimming ) then
        call ReplaceSubString (Reverse(trim(str)), tmpstr, &
          & Reverse(trim(sub2)), separator)
        tmpstr2 = Reverse(trim(tmpstr))
        call ReplaceSubString (tmpstr2, tmpstr, sub1, separator)
      else
        call ReplaceSubString (Reverse(str), tmpstr, &
          & Reverse(sub2), separator, no_trim=.true.)
        tmpstr2 = Reverse(tmpstr)
        call ReplaceSubString (tmpstr2, tmpstr, sub1, separator, no_trim=.true.)
      endif
    case ('stingy')
      if ( trimming ) then
        call ReplaceSubString (Reverse(trim(str)), tmpstr, &
          & Reverse(trim(sub1)), separator)
        tmpstr2 = Reverse(trim(tmpstr))
        call ReplaceSubString (tmpstr2, tmpstr, sub2, separator)
      else
        call ReplaceSubString (Reverse(str), tmpstr, &
          & Reverse(sub1), separator, no_trim=.true.)
        tmpstr2 = Reverse(tmpstr)
        call ReplaceSubString (tmpstr2, tmpstr, sub2, separator, no_trim=.true.)
      endif
    case default
      ! print *, 'Replacing: ', sub1, ' with ', separator
      ! print *, 'in: ', trim(str)
      call ReplaceSubString (str, tmpstr2, sub1, separator, &
        & which='first', no_trim=no_trim)
      ! print *, 'results in: ', trim(tmpstr2)
      ! print *, 'Replacing: ', sub2, ' with ', separator
      call ReplaceSubString (tmpstr2, tmpstr, sub2, separator, &
        & which='first', no_trim=no_trim)
      ! print *, 'results in: ', trim(tmpstr)
    end select
    call GetStringElement (tmpstr, outstr, 2, .true., &
      & inseparator=separator )

  END SUBROUTINE ExtractSubString

  ! ---------------------------------------------  GetIntHashElement  -----

  ! This function takes one (usually) comma-separated string list, interprets it
  ! it as a list of elements, and an array of ints
  ! treating the list as keys and the array as
  ! a hash table, associative array or dictionary
  ! It returns the int from the hash table corresponding to the key
  ! If the key is not found in the array of keys, it sets ErrType=KEYNOTFOUND
  ! otherwise ErrType=0
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns a match for the 
  ! first hash element merely found in the key; e.g.
  ! 'won, to, tree' and key 'protocol.dat' matches 'to'

  ! Basic premise: Use StringElementNum on key in keyList to find index
  ! Use this index for the array of ints
  
  FUNCTION GetIntHashElement(keyList, hashArray, key, ErrType, &
  & countEmpty, inseparator, part_match) RESULT (hashInt)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)             :: keyList
    INTEGER, DIMENSION(:), INTENT(IN)         :: hashArray
    INTEGER                                   :: hashInt
    CHARACTER (LEN=*), INTENT(IN)             :: key
    INTEGER, INTENT(OUT)                      :: ErrType
    LOGICAL, INTENT(IN)                       :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator
    LOGICAL, OPTIONAL, INTENT(IN)             :: part_match

    ! Local variables
	INTEGER :: elem

    ! Executable code

   ErrType = 0
	elem = StringElementNum(keyList, key, countEmpty, inseparator, part_match)
	hashInt = elem
	IF(elem <= 0) THEN
		ErrType = KEYNOTFOUND
	ELSEIF(elem > SIZE(hashArray)) THEN
		ErrType = KEYBEYONDHASHSIZE
	ELSE
		hashInt = hashArray(elem)
	ENDIF

  END FUNCTION GetIntHashElement

  ! ---------------------------------------------  GetStringElement  -----

  ! This subroutine takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the
  ! sub-string which is the n'th element
  ! if n is too large or small, it returns the separator
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! See also SplitWords

  SUBROUTINE GetStringElement(inList, outElement, nElement, countEmpty, inseparator)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: inList
    CHARACTER (LEN=*), INTENT(OUT)  :: outElement
    INTEGER(i4), INTENT(IN)         :: nElement 	! Entry number to return
    LOGICAL, INTENT(IN)   :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator

    ! Local variables
    INTEGER(i4) :: i           ! Loop counters
    INTEGER(i4) :: elem, nextseparator

    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: BLANK = ' '   ! Returned if element empty
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    ! Executable code

    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

	 IF(nElement.LE.0) THEN
	     outElement = separator
	 ELSEIF(LEN(inList) < nElement) THEN
	     outElement = separator
    ENDIF
    i = 1
	 elem = 1
    DO
	     nextseparator = i - 1 + INDEX(inList(i:), separator)

	! No more separators
		  IF(nextseparator == i - 1) THEN
		      IF(elem >= nElement) THEN
				    outElement = inList(i:)
			    ELSE
				    outElement = separator
			    ENDIF
				 RETURN

	! Next separator is the adjacent char
			ELSEIF(nextseparator == i) THEN
				IF(countEmpty) THEN
		     	 IF(elem >= nElement) THEN
				    	outElement = BLANK
						RETURN
			   	 ELSE
					 	elem = elem+1
			    	ENDIF
				ENDIF

	! Until next separator is the next element
			ELSE
		      IF(elem >= nElement) THEN
				    IF(i < nextseparator) THEN
				       outElement = inList(i:nextseparator-1)
						ELSE
				       outElement = separator
						ENDIF
   				 RETURN
			    ELSEIF(nextseparator >= LEN(inList)) THEN
				    outElement = separator
				    RETURN
			    ELSE
					 elem = elem+1
			    ENDIF
			ENDIF
			i = nextseparator+1
	 ENDDO

  END SUBROUTINE GetStringElement

  ! ---------------------------------------------  GetStringHashElement  -----

  ! This subroutine takes two (usually) comma-separated string lists, interprets it
  ! each as a list of elements, treating the first as keys and the second as
  ! a hash table, associative array or dictionary
  ! It returns the sub-string from the hash table corresponding to the key
  ! If the key is not found in the array of keys, it returns the separator
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns a match for the 
  ! first hash element merely found in the key; e.g.
  ! 'won, to, tree' and key 'protocol.dat' matches 'to'

  ! Basic premise: Use StringElementNum on key in keyList to find index
  ! Use this index to GetStringElement from HashList

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  
  SUBROUTINE GetStringHashElement(keyList, hashList, key, outElement, &
  & countEmpty, inseparator, part_match)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: keyList
    CHARACTER (LEN=*), INTENT(IN)   :: hashList
    CHARACTER (LEN=*), INTENT(IN)   :: key
    CHARACTER (LEN=*), INTENT(OUT)  :: outElement
    LOGICAL, INTENT(IN)   :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator
    LOGICAL, OPTIONAL, INTENT(IN)             :: part_match

    ! Local variables
	INTEGER(i4) :: elem
    CHARACTER (LEN=1)                          :: separator
    CHARACTER (LEN=1), PARAMETER               :: COMMA = ','

    ! Executable code

    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

	elem = StringElementNum(keyList, key, countEmpty, inseparator, part_match)
	IF(elem <= 0) THEN
		outElement = separator
	ELSE
		CALL GetStringElement(hashList, outElement, elem, &
        & countEmpty, inseparator)
	ENDIF

  END SUBROUTINE GetStringHashElement

  ! ---------------------------------------------  GetUniqueList  -----

  ! This subroutine takes a string list and returns another containing
  ! only the unique entries. The resulting list is supplied by the caller
  ! (You may safely use the same variable for str and outStr)
  ! E.g., given 'one,two,three,one,four' returns 'one,two,three,four'
  ! If optional string list str2 is supplied, instead
  ! returns list from str that are not also in str2

  SUBROUTINE GetUniqueList(str, outStr, noUnique, countEmpty, &
    & inseparator, IgnoreLeadingSpaces, str2)
    ! Dummy arguments
    CHARACTER (LEN=*), intent(in) :: str
    CHARACTER (LEN=*), intent(out) :: outstr
    INTEGER :: noUnique ! Number of unique entries
    LOGICAL, INTENT(IN)                           :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: inseparator
    LOGICAL, OPTIONAL, INTENT(IN)       :: IgnoreLeadingSpaces
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN)       :: str2

    ! Local variables
    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    CHARACTER (LEN=MAXSTRELEMENTLENGTH), DIMENSION(:), ALLOCATABLE    &
      &                             :: inStringArray, outStringArray, inStrAr2
    integer :: nElems
    integer :: nElems2
    integer :: LongestLen
    integer :: status

    ! Executable code
    IF(PRESENT(inseparator)) THEN
      separator = inseparator
    ELSE
      separator = COMMA
    ENDIF
    if ( len(str) <= 0 .or. len(outstr) <= 0 ) return
    nElems = NumStringElements(str, countEmpty, inseparator, LongestLen)
    noUnique = nElems
    if ( present(str2) ) then
      outStr = ''
      if ( nElems < 1 ) return
    else
      outStr = str
      if ( nElems <= 1 ) return
    endif
    if ( LongestLen > MAXSTRELEMENTLENGTH ) then
      ! print *, 'str: ', trim(str)
      ! print *, 'len(str): ', len(str)
      ! print *, 'LongestLen: ', LongestLen
      ! print *, 'nElems: ', nElems
      call MLSMessage(MLSMSG_Error, ModuleName, &
         & "Element length too long in GetUniqueList")
      return
    end if
    ALLOCATE (inStringArray(nElems), outStringArray(nElems), STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"stringArray in GetUniqueList")
    call list2Array(str, inStringArray, countEmpty, inseparator, &
     & IgnoreLeadingSpaces)
    if ( present(str2) ) then
      nElems2 = NumStringElements(str2, countEmpty, inseparator, LongestLen)
      ALLOCATE (inStrAr2(nElems2), STAT=status)
      IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
           & MLSMSG_Allocate//"stringArray2 in GetUniqueList")
      call list2Array(str2, inStrAr2, countEmpty, inseparator, &
       & IgnoreLeadingSpaces)
      call GetUniqueStrings(inStringArray, outStringArray, noUnique, inStrAr2)
      if ( noUnique > 0 ) then
        call Array2List(outStringArray(1:noUnique), outStr, &
         & inseparator)
      else
        outStr=''
      endif
      DEALLOCATE(inStringArray, outStringArray, inStrAr2)
    else
      call GetUniqueStrings(inStringArray, outStringArray, noUnique)
      if ( noUnique > 0 ) then
        call Array2List(outStringArray(1:noUnique), outStr, &
         & inseparator)
      else
        outStr=''
      endif
      DEALLOCATE(inStringArray, outStringArray)
    endif
  END SUBROUTINE GetUniqueList

  ! ---------------------------------------------  GetUniqueStrings  -----

  ! This subroutine takes an array of strings and returns another containing
  ! only the unique entries. The resulting array is supplied by the caller
  ! If optional extra array is supplied, instead
  ! returns entries from first array not also found in second
  ! Some checking is done to make sure it's appropriate

  SUBROUTINE GetUniqueStrings(inList, outList, noUnique, extra)
    ! Dummy arguments
    CHARACTER (LEN=*), DIMENSION(:) :: inList
    CHARACTER (LEN=*), DIMENSION(:) :: outList
    INTEGER :: noUnique ! Number of unique entries
    CHARACTER (LEN=*), optional, DIMENSION(:) :: extra

    ! Local variables
    INTEGER :: i,j,k           ! Loop counters
    LOGICAL, DIMENSION(:), ALLOCATABLE :: duplicate ! Set if already found
    INTEGER :: status        ! Status from allocate

    INTEGER :: extraSize
    integer :: howManyMax
    INTEGER :: inSize

    ! Executable code, setup arrays
    inSize=SIZE(inList)
    ALLOCATE (duplicate(inSize), STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"duplicate")
    if ( present(extra) ) then
      extraSize=size(extra)
      howManyMax = inSize
      ! print *, 'SIZE(inList) ', inSize
      ! print *, 'SIZE(extra) ', extraSize
    else
      extraSize = -1
      howManyMax = inSize-1 ! Don't bother with last one
    endif
    duplicate = .FALSE.

    ! Go through and find duplicates

    DO i = 1, howManyMax
       IF (.NOT. duplicate(i)) THEN
         if ( extraSize < 1 ) then
          DO j = i+1, inSize
             IF (inList(j)==inList(i)) duplicate(j)=.TRUE.
          END DO
         else
          DO j = 1, extraSize
             IF (extra(j)==inList(i)) duplicate(i)=.TRUE.
          END DO
         endif
       END IF
    END DO

    ! Count how many unique ones there are

    noUnique=count(.NOT. duplicate)

    IF (noUnique>SIZE(outList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "outList too small")
    IF (LEN(outList)<LEN(inList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "outList strings to small")
    outList=""

    if ( noUnique > 0 ) then
      ! do j=1, inSize, 20
      !   print *, (duplicate(j+i), i=0, min(19, inSize-j))
      ! enddo
      j=1
      UniqueLoop: DO i = 1, noUnique
         ! UniqueHuntLoop: DO
         !   IF (.NOT. duplicate(j)) EXIT UniqueHuntLoop
         !   j=j+1
         !   if ( j > inSize ) exit UniqueLoop
         ! END DO UniqueHuntLoop
         k = findFirst(.not. duplicate(j:))
         ! print *, 'j: ', j, '   k: ', k
         if ( k+j-1 > inSize ) then
           call MLSMessage(MLSMSG_Error, ModuleName, &
             & "k goes past array end in GetUniqueStrings")
           outList(i)=inList(inSize)
           return
         elseif ( k > 0 ) then
           outList(i)=inList(k+j-1)  ! was inList(j)
           j = j + k
         else
           exit UniqueLoop
         endif
         ! j=j+1
         if ( j > inSize ) exit UniqueLoop
      END DO UniqueLoop
    endif

    DEALLOCATE(duplicate)
  END SUBROUTINE GetUniqueStrings

  ! ------------------------------------------------  hhmmss_value  -----
  Function hhmmss_value(str, ErrTyp, separator, strict) result (value)
    ! Function that returns the value in seconds of a string 'hh:mm:ss'
    ! where the field separator ':' divides the string into two
    ! integer-like strings 'hh' and 'mm', as well as one float-like
    ! string 'ss' which may have a decimal point plus fractional part
    ! E.g., ss=59.9999
    
    ! Requires 0 <= hh <= 24
    ! Requires 0 <= mm < 60
    ! Requires 0. <= ss < 60.

    ! Returns ErrTyp=0 unless an error occurs

    ! Lenient wrt utc and non-compliant formats:
    ! ignores chars in front of 'hh' and a terminal,
    ! non-numerical char: e.g., '2000-01-01T00:00:00.000000Z'
    ! will be treated the same as '00:00:00.0000000'

    ! If given optional arg strict, not lenient
    ! i.e., non-compliant str always returns non-zero ErrTyp
    
    ! If given optional arg separator, uses separator as field separator
    
    ! Useful to allow an added way to input time
    
    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Argument--------!
    character(len=*),intent(in) :: str
    real(r8) :: value
    integer, intent(out) :: ErrTyp
    character(len=1),intent(in), optional :: separator
    logical,intent(in), optional :: strict
    !----------Local vars----------!
    character(len=1), parameter :: colon=':'
    character(len=1) :: myColon
    character(len=2) :: mm
    character(len=NameLen) :: ss
    character(len=NameLen) :: hh
    character(len=10), parameter :: digits='0123456789'
    character(LEN=*), parameter :: time_conversion='(I2)'
    character(LEN=*), parameter :: real_conversion='(F32.0)'
    integer :: i
    logical :: mystrict
    integer :: hvalue, mvalue
    !----------Executable part----------!

   if(present(separator)) then
      myColon=separator
   else
      myColon=colon
   endif
         
   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(len_trim(str) <= 0) then
      value=0.
      if(mystrict) then
         ErrTyp=INVALIDHHMMSSSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDHHMMSSSTRING
   value=0.
   
   call GetStringElement(str, hh, 1, countEmpty=.true., inseparator=myColon)
   call GetStringElement(str, mm, 2, countEmpty=.true., inseparator=myColon)
   call GetStringElement(str, ss, 3, countEmpty=.true., inseparator=myColon)
   
   ! Check if ss terminates in a non-digit
   ss=Reverse(trim(ss))
   
   if(len_trim(ss) <= 0) then
      if(mystrict) then
         return
      endif      
   elseif( .not. (index(digits, ss(1:1)) > 0) ) then
      if(mystrict) then
         return
      else
         ss=Reverse(trim(ss(2:)))
      endif
   else
      ss=Reverse(trim(ss))
   endif

   do i=1, len_trim(ss)
      if( .not. (index(digits, ss(i:i)) > 0 .or. ss(i:i) == '.') ) return
   enddo

   ! Check if mm complies
   if(len_trim(mm) <= 0) then
      if(mystrict) then
         return
      endif      
   else
      do i=1, len_trim(mm)
        if( .not. (index(digits, mm(i:i)) > 0) ) return
      enddo
   endif      

   ! Check if hh complies
   hh=Reverse(trim(hh))
   
   if(len_trim(hh) <= 0) then
      if(mystrict) then
         return
      endif      
   elseif(mystrict) then
      do i=1, len_trim(hh)
        if( .not. (index(digits, hh(i:i)) > 0) ) return
      enddo
   endif

   hh=Reverse(hh(:2))
   
   ErrTyp=0
   
   ! Convert to value
   if(hh == ' ') then
      hvalue=0
   else
      read(hh(1:2), time_conversion, iostat=ErrTyp) hvalue
   endif
   
   if(ErrTyp /= 0) then
      return
   elseif(hvalue < 0 .or. hvalue > 24) then
      ErrTyp=INVALIDHHMMSSSTRING
      return
   endif

   if(mm == ' ') then
      mvalue=0
   else
      read(mm, time_conversion, iostat=ErrTyp) mvalue
   endif

   if(ErrTyp /= 0) then
      return
   elseif(mvalue < 0 .or. mvalue > 59) then
      ErrTyp=INVALIDHHMMSSSTRING
      return
   endif

   ! Convert to value
   if(ss == ' ') then
      value=0
   else
      read(ss, real_conversion, iostat=ErrTyp) value
   endif
   
   if(ErrTyp /= 0) then
      return
   elseif(value < 0. .or. value > 60.) then
      ErrTyp=INVALIDHHMMSSSTRING
      return
   endif
   
   value = value + 60*(mvalue + 60*hvalue)

  end Function hhmmss_value

  ! ---------------------------------------------  List2Array  -----

  ! This subroutine takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns an equivalent array of
  ! sub-strings in which the n'th element is the n'th element

  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! If the optional arg ignoreLeadingSpaces is TRUE, "a, b, c" is
  ! treated like "a,b,c"; otherwise the leading spaces are retained

  SUBROUTINE List2Array(inList, outArray, countEmpty, inseparator, &
   & IgnoreLeadingSpaces)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)                 :: inList
    CHARACTER (LEN=*), DIMENSION(:), INTENT(OUT)  :: outArray
    LOGICAL, INTENT(IN)                           :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: inseparator
    LOGICAL, OPTIONAL, INTENT(IN)                 :: IgnoreLeadingSpaces

    ! Local variables
    INTEGER(i4) :: elem, nElems

    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: BLANK = ' '   ! Returned for any element empty
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    logical                         :: myIgnoreLeadingSpaces
    ! Executable code

    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

    IF(PRESENT(IgnoreLeadingSpaces)) THEN
	     myIgnoreLeadingSpaces = IgnoreLeadingSpaces
	 ELSE
	     myIgnoreLeadingSpaces = .false.
	 ENDIF

    if ( size(outArray) <= 0 ) return
    outArray = BLANK
	 elem = 1
    nElems = NumStringElements(inList, countEmpty, inseparator)
    if ( nElems <= 0 ) return
    DO
      call GetStringElement(inList, outArray(elem), elem, countEmpty, inseparator)
      if ( myIgnoreLeadingSpaces ) outArray(elem) = adjustl(outArray(elem))
      elem = elem + 1
      if ( elem > min(nElems, size(outArray)) ) return
	 ENDDO

  END SUBROUTINE List2Array

  ! ---------------------------------------------  NumStringElements  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the
  ! number of elements
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists
  !
  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the separator may supplied, in case it isn't comma

  ! See also GetStringElement

  FUNCTION NumStringElements(inList, countEmpty, &
   & inseparator, LongestLen) RESULT (nElements)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)             :: inList
    LOGICAL, INTENT(IN)                       :: countEmpty
	 INTEGER                                   :: nElements
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator
    INTEGER, OPTIONAL, INTENT(OUT)            :: LongestLen  ! Length of longest

    ! Local variables
    INTEGER :: i, sinceLastseparated           ! Loop counters
	 LOGICAL :: lastWasNotseparated

    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    ! Executable code

    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

	! Count the number of separators
   if ( present(LongestLen) ) &
     & LongestLen =0
	! nElements-1 = number of separators
	IF(LEN_TRIM(inList) <= 0) THEN
		nElements=0
      if ( present(LongestLen) ) LongestLen = 0
		RETURN
	ENDIF
	
	lastWasNotseparated = .FALSE.
	nElements = 1
   sinceLastseparated = 0
	DO i=1, LEN_TRIM(inList)
		IF(inList(i:i) == separator) THEN
			IF(countEmpty .OR. lastWasNotseparated) THEN
				nElements = nElements+1
            if ( present(LongestLen) ) &
             & LongestLen = max(LongestLen, sinceLastseparated)
			ENDIF
			lastWasNotseparated = .FALSE.
         sinceLastseparated = 0
		ELSE
			lastWasNotseparated = .TRUE.
         sinceLastseparated = sinceLastseparated + 1
		ENDIF
	ENDDO
   if ( present(LongestLen) ) &
     & LongestLen = max(LongestLen, sinceLastseparated)

  END FUNCTION NumStringElements

  ! --------------------------------------------------  RemoveElemFromList  -----
  SUBROUTINE RemoveElemFromList (inList, outList, elem, inseparator)
    ! Takes a list and removes all occurrence(s) of elem
	 ! E.g., given 'a,b,c,d,..,z' and asked to remove 'c' returns 'a,b,d,..z'
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: inList
    CHARACTER (LEN=*), INTENT(IN) :: elem
    CHARACTER (LEN=*), INTENT(OUT)                :: outList
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: inseparator
    ! Method:
    ! Prepend elem onto start of list, make it unique,
    ! Then snip it back off
    !----------Local vars----------!
    character(len=len(inList)+len(elem)+1) :: temp_list, unique_list
    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: BLANK = ' '   ! Returned for any element empty
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    integer :: numUnique
    !----------Executable part----------!
    IF(PRESENT(inseparator)) THEN
      separator = inseparator
    ELSE
      separator = COMMA
    END IF

    outList = inList
    IF (LEN_trim(elem) < 1 .or. len_trim(inList) < 1 &
      & .or. StringElementNum(inList, elem, countEmpty=.true., &
    & inseparator=inseparator) < 1 ) RETURN
    temp_list = trim(elem) // separator // trim(inList)
    call GetUniqueList(temp_list, unique_list, numUnique, countEmpty=.true., &
    & inseparator=inseparator, ignoreLeadingSpaces=.true.)
    outList = unique_list(len(elem)+1:)
  END SUBROUTINE RemoveElemFromList

  ! --------------------------------------------------  ReplaceSubString  -----
  SUBROUTINE ReplaceSubString (str, outstr, sub1, sub2, which, no_trim)
    ! Takes a string and replaces occurrence(s) of sub1 with sub2
	 ! Defaults to replacing only the first
    ! But if which == 'all' replaces all
    ! or if which == 'last' replaces last
    ! Note that, depending on no_trim, 'all' does the following:
    ! (a) if no_trim == TRUE, multiple passes (up to 100) until no
    !     further replacements are possible
    !     ( which could be bad; e.g., if sub1 is 'sub1' and sub2 is 'sub11'
    !     then (blah)sub1(blah)sub1..' becomes '(blah)sub1111...' )
    ! (b) if no_trim is FALSE or missing, a single pass after chopping
    !     the string up into separate sub1-containing pieces
    !     ( e.g., '(blah)sub1(blah)sub1(blah)..' becomes
    !      '(blah)sub2(blah)sub2(blah)..' )
    !  
    ! Will this still work if sub1 has leading or trailing blanks? 
    ! How about sub2?
    ! Do we need an optional arg, no_trim, say, that will leave them?
    ! Tried coding it, but can't say for sure it works
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=*), INTENT(IN) :: sub1
    CHARACTER (LEN=*), INTENT(IN) :: sub2
    CHARACTER (LEN=*) :: outstr
    character (len=*), intent(in), optional :: which
    logical, intent(in), optional :: no_trim

    !----------Local vars----------!
    integer, parameter         :: MAXREPLACEMENTS = 100
    INTEGER :: i, array_size
    character (len=5) :: my_which
    character(len=max(len(str), len(outstr))) :: head
    character(len=max(len(str), len(outstr))) :: tail
    character(len=max(len(str), len(outstr))) :: sub_str
    character(len=max(len(str), len(outstr))), dimension(MAXREPLACEMENTS) &
      &                                      :: str_array
    logical :: my_no_trim
    !----------Executable part----------!
    head = ''
    tail = ''
    sub_str = ''
    outstr = str
    IF (LEN_trim(str) < 1 .or. len_trim(sub1) < 1) RETURN
    my_which = 'first'
    if ( present(which) ) my_which = lowercase(which)
    my_no_trim = .false.
    if ( present(no_trim) ) my_no_trim = no_trim

    select case (my_no_trim)
    case (.false.)
      if ( index(str, trim(sub1)) < 1 ) return
      select case (my_which)
      case ('first')
        call Replace_me ( str, outstr, .false. )
      case ('last')
        call Replace_me ( str, outstr, .true. )
      case ('all')
        outstr = ' '
        str_array = ' '
        call Split_me
        do i=1, array_size
          ! print *, i, ' ', trim(str_array(i))
          call Replace_me ( trim(str_array(i)), sub_str, .false. )
          outstr = adjustl( trim(outstr) // sub_str )
        enddo
      end select
    case (.true.)
      if ( index(str, sub1) < 1 ) return
      select case (my_which)
      case ('first')
        call Replace_me_no_trim ( str, outstr, .false. )
      case ('last')
        call Replace_me_no_trim ( str, outstr, .true. )
      case ('all')
        ! Originally, I despaired of solving this
        ! CALL MLSMessage(MLSMSG_Error, ModuleName, &
        ! & 'Unable to ReplaceSubStrings with which=all and no_trim=TRUE yet')
        ! Then I had an idea: Why not reinterpret this as multiple passes?
        i = 0
        str_array(1) = str
        do
          i = i + 1
          ! print *, 'i ', i
          ! print *, 'str_array(1) ', str_array(1)
          call Replace_me_no_trim ( str_array(1), str_array(2), .false. )
          ! print *, 'str_array(2) ', str_array(2)
          if ( str_array(2) == str_array(1) .or. i > MAXREPLACEMENTS ) exit
          str_array(1) = str_array(2)
        enddo
        outstr = str_array(2)
      end select
    end select
    
    contains
      subroutine Replace_me ( the_orig, after_sub, back )
        ! This replaces an instance of sub1 with sub2 in
        ! the string the_orig
        ! Either the first instance (if back == FALSE) or the last
        ! Arguments
        character(len=*), intent(in)  :: the_orig
        character(len=*), intent(inout) :: after_sub
        logical, intent(in)           :: back
        ! Local variables
        integer :: istrt1, istrt2, ihead
        if ( index(the_orig, trim(sub1)) == 0 ) then
          after_sub = the_orig
          return
        endif
        istrt1 = index(the_orig, trim(sub1), back=back)
        istrt2 = istrt1 + len_trim(sub1)
        ihead = 1
        head = ' '
        tail = ' '
        if ( istrt1 > 1 ) then
          head = the_orig(1:istrt1-1)
          ihead = istrt1 - 1
        endif
        if ( istrt2 < len_trim(the_orig)+1 ) then
          tail = the_orig(istrt2:)
        endif
        ! print *, 'len(after_sub): ', len(after_sub)
        ! print *, 'len(head): ', len(head)
        ! print *, 'len(tail): ', len(tail)
        ! print *, 'head: ', trim(head), '  ihead: ', ihead
        ! print *, 'tail: ', trim(tail), len_trim(tail)
        ! print *, 'sub1: ', trim(sub1), len_trim(sub1)
        ! print *, 'sub2: ', trim(sub2), len_trim(sub2)
        if ( sub2 /= ' ' ) then
          after_sub = adjustl(head(1:ihead) // trim(sub2) // trim(tail))
        else
          after_sub = adjustl(head(1:ihead) // trim(tail))
        endif
      end subroutine Replace_me
      subroutine Replace_me_no_trim ( the_orig, after_sub, back )
        ! This replaces an instance of sub1 with sub2 in
        ! the string the_orig -- w/o trimming leading or trailing blanks
        ! Either the first instance (if back == FALSE) or the last
        ! Arguments
        character(len=*), intent(in)  :: the_orig
        character(len=*), intent(inout) :: after_sub
        logical, intent(in)           :: back
        ! Local variables
        integer :: istrt1, istrt2, ihead
        if ( index(the_orig, sub1) == 0 ) then
          after_sub = the_orig
          return
        endif
        istrt1 = index(the_orig, sub1, back=back)
        istrt2 = istrt1 + len(sub1)
        ihead = 0
        head = ' '
        tail = ' '
        if ( istrt1 > 1 ) then
          head = the_orig(1:istrt1-1)
          ihead = istrt1 - 1
        endif
        if ( istrt2 < len(the_orig)+1 ) then
          tail = the_orig(istrt2:)
        endif
        ! Now all the possibilities:
        ! (1) the_orig = sub1
        if ( ihead == 0 .and. istrt2 > len(the_orig) ) then
          after_sub = sub2
        ! (2) the_orig = (head)sub1
        elseif ( istrt2 > len(the_orig) ) then
          after_sub = head(1:ihead) // sub2
        ! (3) the_orig = sub1(tail)
        elseif ( ihead == 0 ) then
          after_sub = sub2 // tail
        ! (4) the_orig = (head)sub1(tail)
        else
          after_sub = head(1:ihead) // sub2 // tail
        endif
      end subroutine Replace_me_no_trim
      subroutine Split_me
        ! Will this still work if some of the str_arrays end in one or more ' '?
        ! Arguments (none)
        ! Local variables
        integer :: istrt1, istrt2
        array_size = 0
        istrt2 = 0
        do
          if ( istrt2 > len_trim(str) - 1 ) return
          if ( index(str(istrt2+1:), trim(sub1)) < 1 ) then
            array_size = min(array_size+1, MAXREPLACEMENTS)
            str_array(array_size) = str(istrt2+1:)
            return
          endif
          istrt1 = istrt2 + index(str(istrt2+1:), trim(sub1))
          array_size = min(array_size+1, MAXREPLACEMENTS)
          str_array(array_size) = str(istrt2+1:istrt1 + len_trim(sub1) - 1)
          istrt2 = istrt1 + len_trim(sub1) - 1
        enddo
      end subroutine Split_me
  END SUBROUTINE ReplaceSubString

  ! --------------------------------------------------  ReverseList  -----
  FUNCTION ReverseList (str, inseparator) RESULT (outstr)
    ! takes a string list, usually comma-separated,
	 ! and returns one with elements in reversed order

	 ! E.g., given "alpha, beta, gamma" => "gamma, beta, alpha"

	 ! Limitation:
	 ! No element may be longer than MAXWORDLENGTH
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator

    !----------Local vars----------!
    INTEGER(i4) :: i, istr, irev, elem, iBuf
    INTEGER, PARAMETER :: MAXWORDLENGTH=80
    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    CHARACTER (LEN=1), DIMENSION(:), ALLOCATABLE :: charBuf
    CHARACTER (LEN=MAXWORDLENGTH) :: word
! Treat consecutive separators as if enclosing an empty element
	LOGICAL, PARAMETER :: countEmpty = .TRUE.    

    !----------Executable part----------!
    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

!  Special case--only one element of str
    outstr = str
    IF(LEN(str) == 1 .OR. INDEX(str, separator) == 0) RETURN
	 
! General case
	 ALLOCATE(charBuf(LEN(str)+1), STAT=istr)
    IF (istr /= 0) THEN
	 	CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"charBuf")
		RETURN
	ENDIF

    outstr = ' '

! Loop over elements
	elem = 1
	iBuf=0
	DO
		CALL GetStringElement(str, word, elem, countEmpty, separator)
		IF(word == separator) THEN
			EXIT
		ELSEIF(iBuf > LEN(str)) THEN
			EXIT
		ELSE
			istr = MAX(LEN_TRIM(word), 1)
			word = Reverse(word(:istr))
			DO i=1, istr
				iBuf=iBuf+1
				charBuf(iBuf) = word(i:i)
			ENDDO
			iBuf=iBuf+1
			charBuf(iBuf) = separator
			elem = elem+1
		ENDIF
	ENDDO
	
	IF(charBuf(iBuf) == separator) THEN
		iBuf = iBuf-1
	ENDIF
	
	DO i=1, iBuf
		irev = iBuf - i + 1
		outstr(irev:irev) = charBuf(i)
	ENDDO

	DEALLOCATE(charBuf)

  END FUNCTION ReverseList

  ! ---------------------------------------------  SortArray  -----

  ! This subroutine takes an array of strings
  ! and returns the array of ordered integers
  ! sorting the array; i.e., if ss[n] is the sub-string which is
  ! the n'th element, and ia[k] is the k'th element of the integer array
  ! then {psl[ia[k]]=ss[k], k=1..n} yields the properly sorted array
  ! (unless the further optional arg leftRight is also supplied and equals
  ! one of {"r", "R"} in which case {psl[k]=ss[ia[k]], k=1..n})
  ! Parallel use of ia is how you would normally 
  ! sort any other arrays associated with ss
  
  ! The sorting is ordered by ascii collating sequence:
  ! "0" < "9" < "A" < "Z" < "a" < "z"
  ! unless caseSensitive is FALSE, when "0" < "9" < "A" < "a" < "Z" < "z"

  ! As an optional arg the properly sorted array is returned, too
  ! You may safely supply the same arg for both inStrArray and sortedArray
  ! If the optional arg shorterFirst is TRUE, the sorting is modified
  ! so that shorter strings come first
  ! e.g., (/'abc', 'st', 'Z', '1'/) -> (/'1', 'Z', 'st', 'abc'/)
  
  ! If shorterFirst, leading spaces are always ignored
  ! otherwise they are always significant
  ! (See SortList for contrasting treatment options)
  !  (If you want them ignored, it's easy enough: create a tempArray
  !     tempArray(1:N) = adjustl(strArray(1:N))
  !   and pass it in instead)

  ! Method:
  ! The strings are sifted one character at a time through a series
  ! of ever-finer bins using the selection sort embodied in
  ! subroutine tie_breaker (which see; it surely can be easily improved
  ! upon, but the overall computational gains would be modest)
  ! until each bin is occupied by no more than one string
  ! The bin number is the ranking index of that string which
  ! is returned as outIntArray
  SUBROUTINE SortArray(inStrArray, outIntArray, CaseSensitive, &
   & sortedArray, shorterFirst, leftRight)
    ! Dummy arguments
    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN)   :: inStrArray
    INTEGER, DIMENSION(:), INTENT(OUT)            :: outIntArray
    LOGICAL, INTENT(IN)                           :: caseSensitive
    CHARACTER (LEN=*), DIMENSION(:), OPTIONAL, INTENT(OUT)  &
     &                                            :: sortedArray
    LOGICAL, OPTIONAL, INTENT(IN)                 :: shorterFirst
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: leftRight

    ! Local variables
    INTEGER(i4) :: elem, nElems
    integer, parameter                     :: MAXCHARVALUE = 256
    integer, parameter                     :: MAXELEM = MAXSTRELEMENTLENGTH
!    integer, dimension(MAXELEM)           :: chValue, cvInvBN
!    integer, dimension(MAXELEM)           :: binNumber, invBinNumber
!    integer, dimension(MAXELEM)           :: jsort, inTheBin
    integer, dimension(:), allocatable     :: chValue, cvInvBN
    integer, dimension(:), allocatable     :: binNumber, invBinNumber 
    integer, dimension(:), allocatable     :: jsort, inTheBin
    integer                                :: numBins, oldNumBins
    integer                                :: i, bin, ck, strPos
    integer                                :: status
    integer                                :: maxStrPos
    logical                                :: allTheSameInThisBin
    logical                                :: myShorterFirst
    CHARACTER (LEN=1)                      :: theChar  
    CHARACTER (LEN=1), PARAMETER           :: BLANK = ' '
    CHARACTER (LEN=MAXSTRELEMENTLENGTH), DIMENSION(:), ALLOCATABLE    &
      &                                    :: stringArray
    CHARACTER (LEN=MAXSTRELEMENTLENGTH)    :: theString  
    CHARACTER (LEN=1)                      :: myLeftRight
    logical, parameter                     :: DeeBUG = .false.

    ! Executable code
    IF(PRESENT(shorterFirst)) THEN
	     myshorterFirst = shorterFirst
	 ELSE
	     myshorterFirst = .false.
	 ENDIF
    IF(PRESENT(leftRight)) THEN
	     myleftRight = Capitalize(leftRight)
	 ELSE
	     myleftRight = "L"
	 ENDIF

    nElems = size(inStrArray)
    if ( size(outIntArray) <= 0 .or. nElems <= 0 ) then
      return
!    elseif ( nElems > MAXELEM ) then
!       CALL MLSMessage(MLSMSG_Error, ModuleName, &
!         & 'Too many elements in inStrArray in SortArray')
!       return
    endif
!    ALLOCATE (stringArray(nElems), STAT=status)
    ALLOCATE (stringArray(nElems), chValue(nElems), cvInvBN(nElems), &
     & binNumber(nElems), invBinNumber(nElems), &
     & jsort(nElems), inTheBin(nElems), &
     & STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
         & MLSMSG_Allocate//"stringArray, etc. in SortArray")
    outIntArray = 0
    numBins = 1
    maxStrPos = 1                ! This will hold max string length needed
    do elem = 1, nElems    
      outIntArray(elem) = 1
      if ( myShorterFirst ) then
        maxStrPos = max(maxStrPos, len_trim(adjustl(inStrArray(elem))))
      else
        maxStrPos = max(maxStrPos, len_trim(inStrArray(elem)))
      endif
    enddo                  
    if ( DEEBUG ) then
      do elem = 1, nElems    
        print *, 'Array element ', elem, ' ', trim(inStrArray(elem))
      enddo                  
    endif
    do elem = 1, nElems    
      if ( myshorterFirst ) then
        ! This causes shorter strings to have more leading spaces
        ! and therefore come up first when sorted
        ! (which is why we always ignore leading spaces in inStrArray)
        theString = adjustl(inStrArray(elem))
        stringArray(elem) = adjustr(theString(1:maxStrPos))
      else
        stringArray(elem) = inStrArray(elem)
      endif
    enddo                  
    DO strPos = 1, maxStrPos
      
      if ( DEEBUG ) then
        print *, 'string position: ', strPos
        print *, 'array of bins: ', (outIntArray(elem), elem=1, nElems)
      endif
      do elem = 1, nElems
        theChar = stringArray(elem)(strPos:strPos)
        if ( stringArray(elem) == ' ' ) then
          chValue(elem) = MAXCHARVALUE
        elseif ( theChar == ' ' ) then
          chValue(elem) = 0
        elseif ( CaseSensitive ) then
          chValue(elem) = IACHAR(theChar)
        elseif (IACHAR("a") <= IACHAR(theChar) .and. &
          & IACHAR(theChar) <= IACHAR("z") ) then
          chValue(elem) = 2*IACHAR(Capitalize(theChar)) - IACHAR("A") + 1
        elseif (IACHAR("A") <= IACHAR(theChar) .and. &
          & IACHAR(theChar) <= IACHAR("Z") ) then
          chValue(elem) = 2*IACHAR(theChar) - IACHAR("A")
        elseif (IACHAR("Z") < IACHAR(theChar) ) then
          chValue(elem) = 2*IACHAR(theChar)
        else
          chValue(elem) = IACHAR(theChar)
        endif
      enddo
      if ( DEEBUG ) print *, 'array of chValues: ', (chValue(elem), elem=1, nElems)
      oldNumBins = numBins
      do elem=1, nElems
        binNumber(elem) = outIntArray(elem)
      enddo
      numBins = 0
      ck = 0
      if ( DEEBUG ) print *, 'number of bins: ', oldNumBins
      do bin=1, oldNumBins
        if ( DEEBUG ) print *, 'bin number: ', bin
        call warm_up(bin)
        if ( DEEBUG ) then
          print *, 'number in bin: ', inTheBin(bin)
          print *, 'array of invBinNumber: ', &
           & (invBinNumber(elem), elem=1, inTheBin(bin))
          print *, 'array of cvInvBN: ', (cvInvBN(elem), elem=1, inTheBin(bin))
        endif
        call tie_breaker(bin)
        if ( DEEBUG ) print *, 'array of jsort: ', (jsort(elem), elem=1, inTheBin(bin))
        numBins = numBins + 1
        ck = cvInvBN(jsort(1))
        do i=1, inTheBin(bin)
          if ( ck /= cvInvBN(jsort(i)) .or. &
           & ( &
           &  allTheSameInThisBin .and. i > 1 &
           & ) &
           & ) then
            numBins = numBins + 1
            ck = cvInvBN(jsort(i))
          endif
          outIntArray(invBinNumber(jsort(i))) = numBins
        enddo
      enddo
      if ( numBins >= min(nElems, size(outIntArray)) ) exit
	 ENDDO
    if ( DEEBUG ) then
      print *, 'Final number of bins: ', numBins
      print *, 'Sorting order: ', (outIntArray(i), i=1, nElems)
    endif

    if ( present(sortedArray) ) then
      do elem=1, nElems
        i = max(1, outIntArray(elem))
        i = min(i, nElems, size(sortedArray))
        if ( myShorterFirst ) then
          sortedArray(i) = adjustl(stringArray(elem))
        else
          sortedArray(i) = stringArray(elem)
        endif
      enddo
    endif
    if ( myLeftRight == 'R' ) then
      ! Need to 'invert' outIntArray
      do elem=1, nElems
        do i=1, nElems
          if ( outIntArray(i) == elem ) invBinNumber(elem) = i
        enddo
      enddo
      do elem=1, nElems
        outIntArray(elem) = invBinNumber(elem)
      enddo
    endif
    DEALLOCATE(stringArray, chValue, cvInvBN, binNumber, invBinNumber, &
     & jsort, inTheBin, STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
         & MLSMSG_DeAllocate//"stringArray, etc. in SortArray")

   contains
     subroutine warm_up(theBin)
       ! Form array invBinNumber = {i[j], j=1 .. }
       ! such that binNumber[i] = theBin
       ! Then form cvInvBN = {c[j] = chValue[i[j]], j=1..}
       integer, intent(in) :: theBin
       integer :: j, i
       j=0
       do i=1, nElems
         if ( binNumber(i) == theBin ) then
           j=j+1
           invBinNumber(j) = i
           cvInvBN(j) = chValue(i)
         endif
       enddo
       inTheBin(theBin) = j
     end subroutine warm_up
     subroutine tie_breaker(theBin)
       ! Form array jsort = j[k] = {j_1, j_2, .., j_N}
       ! sorted so that c[j_1} <= c[j_2] <= .. <= c{j_N]
       ! This is a naive selection sort--make any improvements you wish
       ! (Order N^2 sorting algorithms are inefficient)
       integer, intent(in)      :: theBin
       integer                  :: kp, k, ck, jsortie
       CHARACTER (LEN=MAXSTRELEMENTLENGTH)  :: stringElement  
       allTheSameInThisBin = (inTheBin(theBin) /= 1)
       stringElement = stringArray(invBinNumber(1))
       do k=1, inTheBin(theBin)
         jsort(k) = k
         allTheSameInThisBin = allTheSameInThisBin .and. &
           & stringElement == stringArray(invBinNumber(k))
       enddo
       if ( inTheBin(theBin) == 1 .or. allTheSameInThisBin) return
       do k=1, inTheBin(theBin) - 1
         ck = cvInvBN(jsort(k))
         do kp=k+1, inTheBin(theBin)
           if ( cvInvBN(jsort(kp)) < ck ) then
           ! Pull the old switcheroo
             ck = cvInvBN(jsort(kp))
             jsortie = jsort(kp)
             jsort(kp) = jsort(k)
             jsort(k) = jsortie
           endif
         enddo
       enddo
     end subroutine tie_breaker

  END SUBROUTINE SortArray

  ! ---------------------------------------------  SortList  -----

  ! This subroutine takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the array of ordered integers
  ! sorting the list; i.e., if ss[n] is the sub-string which is
  ! the n'th element, and ia[k] is the k'th element of the integer array
  ! then {psl[ia[k]]=ss[k], k=1..n} yields the properly sorted list
  ! (unless the further optional arg leftRight is also supplied and equals
  ! one of {"r", "R"} in which case {psl[k]=ss[ia[k]], k=1..n})
  ! Parallel use of ia is how you would normally 
  ! sort any other arrays associated with ss
  
  ! The sorting is ordered by ascii collating sequence:
  ! "0" < "9" < "A" < "Z" < "a" < "z"
  ! unless caseSensitive is FALSE, when "0" < "9" < "A" < "a" < "Z" < "z"

  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the properly sorted list is returned, too
  ! You may safely supply the same arg for both inList and sortedList
  ! As an optional arg the separator may supplied, in case it isn't comma
  ! If the optional arg ignoreLeadingSpaces is TRUE, "a, b, c" is
  ! sorted like "a,b,c"; otherwise the leading spaces make" b, c,a"

  ! Method:
  ! (see SortArray)
  SUBROUTINE SortList(inList, outArray, CaseSensitive, countEmpty, &
   & ignoreLeadingSpaces, inseparator, sortedList, leftRight)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)                 :: inList
    INTEGER, DIMENSION(:), INTENT(OUT)            :: outArray
    LOGICAL, INTENT(IN)                           :: CaseSensitive
    LOGICAL, INTENT(IN)                           :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: inseparator
    CHARACTER (LEN=*), OPTIONAL, INTENT(OUT)      :: sortedList
    LOGICAL, OPTIONAL, INTENT(IN)                 :: IgnoreLeadingSpaces
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)       :: leftRight

    ! Local variables
    integer, parameter              :: MAXELEM = MAXSTRELEMENTLENGTH
    INTEGER(i4) :: nElems, status, LongestLen

    CHARACTER (LEN=1)               :: separator
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    CHARACTER (LEN=MAXSTRELEMENTLENGTH), DIMENSION(:), ALLOCATABLE    &
      &                             :: stringArray
    CHARACTER (LEN=1)               :: myLeftRight
    logical, parameter              :: DeeBUG = .false.
    ! Executable code
    IF(PRESENT(inseparator)) THEN
	     separator = inseparator
	 ELSE
	     separator = COMMA
	 ENDIF

    IF(PRESENT(leftRight)) THEN
	     myleftRight = Capitalize(leftRight)
	 ELSE
	     myleftRight = "L"
	 ENDIF

    if ( DEEBUG ) then
       print *, 'Entered SortList'
       print *, 'present(inseparator)?: ', PRESENT(inseparator)
       print *, 'separator: ', separator
       print *, 'string: ', trim(inList)
    endif
    if ( size(outArray) <= 0 ) return
    outArray = 0
    nElems = NumStringElements(inList, countEmpty, inseparator, LongestLen)
    if ( nElems <= 0 ) then
      return
    elseif ( LongestLen > MAXSTRELEMENTLENGTH ) then
      call MLSMessage(MLSMSG_Error, ModuleName, &
         & "Element length too long in SortList")
      return
!    elseif ( nElems > MAXELEM ) then
!      call MLSMessage(MLSMSG_Error, ModuleName, &
!         & "Too many elements needed in SortList")
!      return
    endif
    ALLOCATE (stringArray(nElems), STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"stringArray in SortList")
    call list2Array(inList, stringArray, countEmpty, inseparator, &
     & IgnoreLeadingSpaces)
    call SortArray(stringArray(1:nElems), outArray, CaseSensitive, &
     & leftRight=leftRight)
    if ( present(sortedList) ) then
      if ( myLeftRight == 'R' ) then
        call Array2List(stringArray(1:nElems), sortedList, &
         & inseparator, outArray, leftRight='R')
      else
        call Array2List(stringArray(1:nElems), sortedList, &
         & inseparator, outArray, leftRight='L')
      endif
    endif
    DEALLOCATE(stringArray)

  END SUBROUTINE SortList

  ! ---------------------------------------------  StringElementNum  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements, and a test string which may be an element
  ! It returns the element number of the test string in the string list
  ! or, 0 if the test string is not found
  
  ! Any leading blanks are disregarded before making the comparison;
  ! e.g., 'stare' is the same as ' stare' and is the second element of 
  ! the list 'lex, stare, decisis'
  
  ! Note: if there are multiple matches between the test string and elements
  ! of inList we return only the first
  
  ! If you want the last instead, use ReverseList on inList && subtract
  ! the answer from nElements
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists
  
  ! It will be the immediate precursor function in a hash table
  ! == aka associative array == aka dictionary
  !
  ! If countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns the number of the 
  ! first element merely found in the test string; e.g.
  ! 'won, to, tree' and test 'protocol.dat' returns 2

  ! See also GetStringElement, NumStringElements

  FUNCTION StringElementNum(inList, test_string, countEmpty, &
    & inseparator, part_match) RESULT (elem)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)             :: inList
    CHARACTER (LEN=*), INTENT(IN)             :: test_string
    LOGICAL, INTENT(IN)                       :: countEmpty
	 INTEGER                                   :: elem
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inseparator
    LOGICAL, OPTIONAL, INTENT(IN)             :: part_match

    ! Local variables
    INTEGER :: nElements
    INTEGER , PARAMETER :: MAXELEMENTLENGTH = 80

    CHARACTER (LEN=MAXELEMENTLENGTH)           :: listElement
    logical ::                                    match
!    CHARACTER (LEN=1)                          :: separator
!    CHARACTER (LEN=1), PARAMETER               :: COMMA = ','
    ! Executable code

	nElements = NumStringElements(inList, countEmpty, inseparator)
	
	IF(nElements <= 0) THEN
		elem = 0
		RETURN
	ENDIF
   match = .false.
   if ( present(part_match) ) match = part_match

	! Check for matches--snipping off any leading blanks
	DO elem=1, nElements
		CALL GetStringElement(inList, listElement, elem, countEmpty, inseparator)
      if ( match ) then
        if (trim(listElement) /= ' ' .and. &
          & index(trim(test_string), trim(listElement)) > 0) RETURN
      else
  	     IF(adjustl(listElement) == adjustl(test_string)) RETURN
      endif
	ENDDO
	
	elem = 0

  END FUNCTION StringElementNum

  ! ------------------------------------------------  unquote  -----
  Function unquote(str, quotes, cquotes, strict, stripany) result (outstr)
    ! Function that removes a single pair of surrounding quotes from string

    ! E.g., given "Let me see." or 'Let me see.' returns
    !    Let me see.
    ! If no surrounding quotes are found, returns string unchanged; unless
    ! (1) mismatched quotes, e.g. 'Let me see." will:
    !     remove leading quote but leave trailing quote
    ! (2) a single unpaired quote found at beginning or end, will:
    !  (a) remove it if the resulting string is non-empty; or
    !  (b) return the single unpaired quote if that was the entire str
    
    ! If given optional arg strict, options (1) and (2) above disregarded
    ! i.e., surrounding quotes must match, else returns string unchanged
    
    ! If given optional arg stripany, any quotes, surrounding or internal,
    ! will be removed
    
    ! If given optional arg quotes, removes only surrounding pair:
    ! quotes[i:i] for each i=1..len[quotes]
    ! E.g., given /a\ regexp/ with quotes='/' returns
    !    a\ regexp
    
    ! If given optional args quotes & cquotes, removes only surrounding pair:
    ! quotes[i:i] on the left, cquotes[i:i] on the right, i=1..len[quotes]
    ! E.g., given [a particle] with quotes='[' cquotes=']' returns
    !    a particle
    ! (For this case, strict matching is always on)
    
    ! Useful because the parser will return quote-surrounded strings if that's
    ! how they appear in the lcf
    
    ! Calling get_string with "strip=.true." renders this unnecessary.
    ! However, you might find another use for it, especially with
    ! feature of being able to trim other, user-supplied detritus:
    ! e.g., braces, parentheses, extraneous separators
    
    ! (Aside from switches, we haven't found such a use so far;
    ! instead see more powerful ExtractSubString or ReplaceSubString)
    !--------Argument--------!
    character(len=*),intent(in) :: str
    character(len=len(str)) :: outstr, tmpstr
    character(len=*),intent(in), optional :: quotes
    character(len=*),intent(in), optional :: cquotes
    logical,intent(in), optional :: strict
    logical,intent(in), optional :: stripany
    !----------Local vars----------!
    character(len=1), parameter :: sq=''''
    character(len=1), parameter :: dq='"'
    integer :: first, last, ult, prim
    character(len=1) :: quote, cquote
    integer :: i
    logical :: mystrict
    logical :: mystripany
    !----------Executable part----------!

   ult = len_trim(str)    ! Position of last non-blank char
   prim = ult - len_trim(adjustl(str)) + 1    ! Position of 1st non-blank char
   outstr=str
      
   ! length of non-blank portion of string to be trimmed must be at least 2
   if(ult-prim+1 <= 1) then
      outstr=str
      return
   endif

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
   
   if(present(stripany)) then
      mystripany=stripany
   else
      mystripany=.false.
   endif
   
   ! These are initialized so that if no matching quotes found
   ! we will return    outstr = adjustl(str)
   first = prim
   last = ult

   ! trim surrounding user-supplied marks?

   if(present(quotes)) then
      if(len_trim(quotes) <= 0) then
       outstr=str
       return
      endif
      
      ! Loop over char class in string quotes
      do i=1, len_trim(quotes)
      
         quote = quotes(i:i)
         
         ! Stripany option in force?
         if ( mystripany ) then
            ! print *, 'Replacing ', quote, ' in ', trim(outstr)
            call ReplaceSubString(outstr, tmpstr, quote, '', 'all')
            outstr = tmpstr
            ! print *, trim(outstr)
            cycle
         endif

         ! Supplied with paired left and right quotes?
         if(present(cquotes)) then
            cquote=cquotes(i:i)
            mystrict=.true.
         else
            cquote=quote
         endif

         if(mystrict) then
           if(str(prim:prim) == quote .and. str(ult:ult) == cquote) then
               outstr=str(prim+1:ult-1)
               return
          endif
      
        else
         if(str(prim:prim) == quote) then
           first=prim+1
          endif

          if(str(ult:ult) == cquote) then
             last=ult-1
           endif
        endif

      enddo

   ! insist surrounding marks match?
   elseif(present(strict)) then
      if( &
      & str(prim:prim) == str(ult:ult) &
      & .and. &
        & (str(prim:prim) == sq .or. str(prim:prim) == dq) &
        & ) then
            outstr=str(prim+1:ult-1)
          else
            outstr=str
          endif
         return
      
   elseif(str(prim:prim) == sq) then
      first=prim+1
      if(str(ult:ult) == sq) then
         last=ult-1
      endif
      
   elseif(str(prim:prim) == dq) then
      first=prim+1
      if(str(ult:ult) == dq) then
         last=ult-1
      endif

   else
      first=prim
      if(str(ult:ult) == dq .or. str(ult:ult) == sq) then
         last=ult-1
      endif

   endif
   if ( mystripany ) then
       return
   elseif(last >= first) then
       outstr=str(first:last)
   else
       outstr=str
   endif
      
  end Function unquote

  ! ---------------------------------------------  utc_to_yyyymmdd_ints  -----
  subroutine utc_to_yyyymmdd_ints(str, ErrTyp, year, month, day, strict, nodash)
    ! Routine that returns the year, month, and day from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! The date substring in subdivided by the separator '-'
    ! into either two or three fields
    ! In case (A), the 3 fields are year, month, and day of month
    ! in case (B) the two fields are year and day of year
    
    ! For case (A) returns year, month, and day=day of month
    ! For case (B) returns year, month=-1, and day=day of year
    ! Useful to decode utc inputs into attribute values
    
    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Argument--------!
    character(len=*),intent(in) :: str
    integer, intent(out) :: ErrTyp
    integer, intent(out) :: year
    integer, intent(out) :: month
    integer, intent(out) :: day
    logical, intent(in), optional :: strict
    logical, intent(in), optional :: nodash   ! No dash separating date fields
    !----------Local vars----------!
    character(len=1), parameter :: dash='-'
    character(len=NameLen) :: date
    character(len=NameLen) :: yyyy
    character(len=NameLen) :: mm
    character(len=NameLen) :: dd
    character(LEN=*), parameter :: time_conversion='(I4)'
    logical :: mystrict
    logical :: mynodash
    character(len=1) :: utc_format        ! 'a' or 'b'
    ! The following arrys contains the maximum permissible day for each month
    ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
    integer, dimension(-1:12), parameter :: DAYMAX = (/ &
      & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
      & /)
    !----------Executable part----------!

   year = -1
   month = -1
   day = -1
   mm = ' '

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(present(nodash)) then
      mynodash=nodash
   else
      mynodash=.false.
   endif
         
   if(len_trim(str) <= 0) then
      if(mystrict) then
         ErrTyp=INVALIDUTCSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDUTCSTRING
   ! Snip off time fields from date fields
   call GetStringElement(lowercase(str), date, 1, &
     & countEmpty=.true., inseparator='t')
   if ( date == ' ' ) then
     if ( .not. mystrict) Errtyp = 0
     return
   endif
   if ( myNoDash ) then
     yyyy = date(1:4)
     mm = date(5:6)
     dd = date(7:8)
   else
     call GetStringElement(trim(date), yyyy, 1, countEmpty=.true., inseparator=dash)
     if ( &
       & NumStringElements(trim(date), countEmpty=.true., inseparator=dash) == 2) then
       call GetStringElement(trim(date), dd, 2, countEmpty=.true., inseparator=dash)
       utc_format = 'b'
     else
       call GetStringElement(trim(date), mm, 2, countEmpty=.true., inseparator=dash)
       call GetStringElement(trim(date), dd, 3, countEmpty=.true., inseparator=dash)
       utc_format = 'a'
     endif
   endif
   
   ErrTyp=0
   
   ! Convert to value
   if(yyyy /= ' ') then
      read(yyyy, time_conversion, iostat=ErrTyp) year
   endif
   
   if(ErrTyp /= 0) then
      return
   elseif(year < 0 .or. year > YEARMAX) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif

   if(mm /= ' ') then
      read(mm, time_conversion, iostat=ErrTyp) month
   endif

   if(utc_format == 'b') then
     ErrTyp = 0
     month = -1
   elseif(ErrTyp /= 0) then
      return
   elseif(month < 1 .or. month > 12) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif
   ! Coming out of the above, month should be in the interval [-1, 12]

   if(dd /= ' ') then
      read(dd, time_conversion, iostat=ErrTyp) day
   endif

   if(ErrTyp /= 0) then
      return
   elseif(day < 1 .or. day > DAYMAX(month)) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif
  end subroutine utc_to_yyyymmdd_ints

  ! ---------------------------------------------  utc_to_yyyymmdd_strs  -----
  subroutine utc_to_yyyymmdd_strs(str, ErrTyp, year, month, day, &
    & strict, utcAt0z)
    ! Routine that returns the year, month, and day from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! The date substring in subdivided by the separator '-'
    ! into either two or three fields
    ! In case (A), the 3 fields are year, month, and day of month
    ! in case (B) the two fields are year and day of year
    
    ! For case (A) returns year, month, and day=day of month
    ! For case (B) returns year, month=-1, and day=day of year
    ! Useful to decode utc inputs into attribute values
    
    ! Optionally returns the input string in utcAt0z modified so that 
    ! the hh:mm:ss.sss is 00:00:00Z
    
    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Argument--------!
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: ErrTyp
    character(len=*), intent(out) :: year
    character(len=*), intent(out) :: month
    character(len=*), intent(out) :: day
    logical,intent(in), optional  :: strict
    character(len=*),intent(out), optional   :: utcAt0z
    !----------Local vars----------!
    character(len=1), parameter :: dash='-'
    character(len=NameLen) :: date
    logical :: mystrict
    character(len=1) :: utc_format        ! 'a' or 'b'
    character(len=*), parameter :: chars_0z = 'T00:00:00Z'
    !----------Executable part----------!

   year = ' '
   month = ' '
   day = ' '

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(len_trim(str) <= 0) then
      if(mystrict) then
         ErrTyp=INVALIDUTCSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDUTCSTRING
   ! Snip off time fields from date fields
   call GetStringElement(lowercase(str), date, 1, &
     & countEmpty=.true., inseparator='t')
   if ( date == ' ' ) then
     if ( .not. mystrict) Errtyp = 0
     if ( present(utcAt0z) ) utcAt0z = ' '
     return
   endif
   if ( present(utcAt0z) ) utcAt0z = trim(date) // chars_0z
   call GetStringElement(trim(date), year, 1, countEmpty=.true., inseparator=dash)
   if ( &
     & NumStringElements(trim(date), countEmpty=.true., inseparator=dash) == 2) then
     call GetStringElement(trim(date), day, 2, countEmpty=.true., inseparator=dash)
     utc_format = 'b'
   else
     call GetStringElement(trim(date), month, 2, countEmpty=.true., inseparator=dash)
     call GetStringElement(trim(date), day, 3, countEmpty=.true., inseparator=dash)
     utc_format = 'a'
   endif
   
   ErrTyp=0
   
  end subroutine utc_to_yyyymmdd_strs

  ! ---------------------------------------------  yyyymmdd_to_dai_ints  -----
  subroutine yyyymmdd_to_dai_ints(yyyy, mm, dd, dai, startingDate)
    ! Routine that returns the number of days after a starting date
    ! from 3 ints: the form yyyymmdd
    !--------Argument--------!
    integer ,intent(in) :: yyyy
    integer ,intent(in) :: mm
    integer ,intent(in) :: dd
    integer,intent(out) :: dai
    character(len=*),intent(in),optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    character(len=8) :: mystartingDate
    integer :: yyyy1, mm1, dd1, doy1
    integer :: yyyy2, doy2
    integer :: ErrTyp
    logical :: daiNegative
    integer :: y
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate='20010101'
   endif
   call utc_to_yyyymmdd_ints(mystartingDate, ErrTyp, yyyy1, mm1, dd1, nodash=.true.)
   call yyyymmdd_to_doy_str(mystartingDate, doy1)
   call yyyymmdd_to_doy_ints(yyyy, mm, dd, doy2)
   yyyy2 = yyyy
   daiNegative = yyyy1 > yyyy2
   if ( daiNegative ) then
     call switch(yyyy1, yyyy2)
     call switch(doy1, doy2)
   elseif ( yyyy1 == yyyy2 ) then
     dai = doy2 - doy1
     return
   endif
   dai = doy2 - doy1
   do y = yyyy1, yyyy2 - 1
     if ( leapyear(y) ) then
       dai = dai + DAYMAXLY(-1)
     else
       dai = dai + DAYMAXNY(-1)
     endif
   enddo
   if ( daiNegative ) dai = -dai
  end subroutine yyyymmdd_to_dai_ints

  ! ---------------------------------------------  yyyymmdd_to_dai_str  -----
  subroutine yyyymmdd_to_dai_str(str, dai, startingDate)
    ! Routine that returns the number of days after a starting date
    ! from a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*),intent(in) :: str
    integer,intent(out) :: dai
    character(len=*),intent(in),optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    character(len=8) :: mystartingDate
    integer :: yyyy1, mm1, dd1, doy1
    integer :: yyyy2, mm2, dd2, doy2
    integer :: ErrTyp
    logical :: daiNegative
    integer :: y
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate='20010101'
   endif
   call utc_to_yyyymmdd_ints(mystartingDate, ErrTyp, yyyy1, mm1, dd1, nodash=.true.)
   call utc_to_yyyymmdd_ints(str, ErrTyp, yyyy2, mm2, dd2, nodash=.true.)
   call yyyymmdd_to_doy_str(mystartingDate, doy1)
   call yyyymmdd_to_doy_str(str, doy2)
   daiNegative = yyyy1 > yyyy2
   if ( daiNegative ) then
     call switch(yyyy1, yyyy2)
     call switch(doy1, doy2)
   elseif ( yyyy1 == yyyy2 ) then
     dai = doy2 - doy1
     return
   endif
   dai = doy2 - doy1
   do y = yyyy1, yyyy2 - 1
     if ( leapyear(y) ) then
       dai = dai + DAYMAXLY(-1)
     else
       dai = dai + DAYMAXNY(-1)
     endif
   enddo
   if ( daiNegative ) dai = -dai
  end subroutine yyyymmdd_to_dai_str

!=============================================================================
  ! ---------------------------------------------  switch_ints  -----
  subroutine switch_ints(x1, x2)
    ! Switch args x1 <=> x2
    !--------Argument--------!
    integer,intent(inout) :: x1
    integer,intent(inout) :: x2
    !----------Local vars----------!
    integer :: x
    x = x1
    x1 = x2
    x2 = x
  end subroutine switch_ints

  ! ---------------------------------------------  yyyymmdd_to_doy_ints  -----
  subroutine yyyymmdd_to_doy_ints(year, month, day, doy)
    ! Routine that returns the number of days after the year's start
    ! for year, month, day
    !--------Argument--------!
    integer, intent(in) :: year, month, day
    integer, intent(out) :: doy
    !----------Local vars----------!
    integer :: ErrTyp
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
    ! integer :: m
    ! integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     call utc_to_yyyymmdd_ints(str, ErrTyp, year, month, day, nodash=.true.)
     call yyyymmdd_to_doy_ints(year, month, day, doy)
  end subroutine yyyymmdd_to_doy_str

  logical function leapyear(year)
    integer,intent(in) :: year
     ! This is to capture rule that centuries are leap only
     ! if divisible by 400
     ! Otherwise, as prehaps you knew, leapyears are those years divisible by 4
     if ( 100 * (year/100) >= year ) then
       leapyear = ( 400 * (year/400) >= year )
     else
       leapyear = ( 4 * (year/4) >= year )
     endif
  end function leapyear

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSStringLists
!=============================================================================

! $Log$
! Revision 2.1  2004/08/04 23:17:30  pwagner
! First commit
!
