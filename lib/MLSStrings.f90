! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSStrings               ! Some low level string handling stuff
!=============================================================================

  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  USE MLSCommon, only: i4, r8, NameLen

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

!
! This module contains some low level string handling stuff for mls

!     c o n t e n t s
!     - - - - - - - -

! Capitalize         tr[a-z] -> [A-Z]
! CompressString     Removes all leading and embedded blanks
! count_words        Counts the number of words in a string
! depunctuate        Replaces punctuation with blanks
! GetIntHashElement  Returns int from array of hashes corresponding to key string
! GetStringElement   Returns n'th element of string list
! GetStringHashElement   Returns string from hash list corresponding to key string
! GetUniqueStrings   Returns array of only unique entries from input array
! hhmmss_value       Converts 'hh:mm:ss' formatted string to a real r8
! LinearSearchStringArray     Finds string index of substring in array of strings
! LowerCase          tr[A-Z] -> [a-z]
! NumStringElements  Returns number of elements in list of strings
! ReadCompleteLineWithoutComments     Knits continuations, snips comments
! Reverse            turns 'a string' -> 'gnirts a'
! ReverseList        turns 'abc, def, ghi' -> 'ghi, def, abc'
! SplitWords         Splits 'first, the, rest, last' -> 'first', 'the, rest', 'last'
! StringElementNum   Returns element number of test string in string list
! unquote            Removes surrounding [quotes]

! in the above, a list is a string of usu. comma-delimited words
! an array is a Fortran array of strings or integers
! a hash is a list of key strings and either
! (1) a list of associated strings
! (2) an array of associated integers
!

CONTAINS

  ! -------------------------------------------------  CAPITALIZE  -----
  FUNCTION Capitalize (str) RESULT (outstr)
    ! takes a-z and replaces with A-Z 
    ! leaving other chars alone
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    INTEGER :: i, icode
    integer, parameter :: offset=IACHAR("A")-IACHAR("a")
    !----------Executable part----------!
    outstr=str

    DO i=1, LEN(str)
       icode=IACHAR(outstr(i:i))
       IF ( icode >=IACHAR("a") .AND. icode <= IACHAR("z")) THEN
          outstr(i:i)=achar(icode+offset)
       END IF
    END DO

  END FUNCTION Capitalize

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
  Function depunctuate(str) result (outstr)
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

  end Function depunctuate

  ! ---------------------------------------------  GetIntHashElement  -----

  ! This function takes one (usually) comma-delimited string list, interprets it
  ! it as a list of elements, and an array of ints
  ! treating the list as keys and the array as
  ! a hash table, associative array or dictionary
  ! It returns the int from the hash table corresponding to the key
  ! If the key is not found in the array of keys, it sets ErrType=KEYNOTFOUND
  ! otherwise ErrType=0
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-delimited lists

  ! If countEmpty is TRUE, consecutive delimiters, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single delimiter
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the delimiter may supplied, in case it isn't comma

  ! Basic premise: Use StringElementNum on key in keyList to find index
  ! Use this index for the array of ints
  
  FUNCTION GetIntHashElement(keyList, hashArray, key, ErrType, &
  & countEmpty, inDelim) RESULT (hashInt)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: keyList
    INTEGER(i4), DIMENSION(:), INTENT(IN)   :: hashArray
    INTEGER(i4)                             :: hashInt
    CHARACTER (LEN=*), INTENT(IN)   :: key
    INTEGER, INTENT(OUT)  :: ErrType
    LOGICAL, INTENT(IN)   :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    ! Local variables
	INTEGER, PARAMETER :: KEYNOTFOUND=-1
	INTEGER, PARAMETER :: KEYBEYONDHASHSIZE=KEYNOTFOUND-1
	INTEGER :: elem
!    CHARACTER (LEN=1)                          :: Delim
!    CHARACTER (LEN=1), PARAMETER               :: COMMA = ','

    ! Executable code

	elem = StringElementNum(keyList, key, countEmpty, inDelim)
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

  ! This subroutine takes a (usually) comma-delimited string list, interprets it
  ! as a list of individual elements and returns the
  ! sub-string which is the n'th element
  ! if n is too large or small, it returns the delimiter
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-delimited lists

  ! If countEmpty is TRUE, consecutive delimiters, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single delimiter
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the delimiter may supplied, in case it isn't comma
  ! See also SplitWords

  SUBROUTINE GetStringElement(inList, outElement, nElement, countEmpty, inDelim)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: inList
    CHARACTER (LEN=*), INTENT(OUT)  :: outElement
    INTEGER(i4), INTENT(IN)         :: nElement 	! Entry number to return
    LOGICAL, INTENT(IN)   :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    ! Local variables
    INTEGER(i4) :: i           ! Loop counters
    INTEGER(i4) :: elem, nextDelim

    CHARACTER (LEN=1)               :: Delim
    CHARACTER (LEN=1), PARAMETER    :: BLANK = ' '   ! Returned if element empty
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    ! Executable code

    IF(PRESENT(inDelim)) THEN
	     Delim = inDelim
	 ELSE
	     Delim = COMMA
	 ENDIF

	 IF(nElement.LE.0) THEN
	     outElement = Delim
	 ELSEIF(LEN(inList) < nElement) THEN
	     outElement = Delim
    ENDIF
    i = 1
	 elem = 1
    DO
	     nextDelim = i - 1 + INDEX(inList(i:), Delim)

	! No more delimiters
		  IF(nextDelim == i - 1) THEN
		      IF(elem >= nElement) THEN
				    outElement = inList(i:)
			    ELSE
				    outElement = Delim
			    ENDIF
				 RETURN

	! Next delimiter is the adjacent char
			ELSEIF(nextDelim == i) THEN
				IF(countEmpty) THEN
		     	 IF(elem >= nElement) THEN
				    	outElement = BLANK
						RETURN
			   	 ELSE
					 	elem = elem+1
			    	ENDIF
				ENDIF

	! Until next delimiter is the next element
			ELSE
		      IF(elem >= nElement) THEN
				    IF(i < nextDelim) THEN
				       outElement = inList(i:nextDelim-1)
						ELSE
				       outElement = Delim
						ENDIF
   				 RETURN
			    ELSEIF(nextDelim >= LEN(inList)) THEN
				    outElement = Delim
				    RETURN
			    ELSE
					 elem = elem+1
			    ENDIF
			ENDIF
			i = nextDelim+1
	 ENDDO

  END SUBROUTINE GetStringElement

  ! ---------------------------------------------  GetStringHashElement  -----

  ! This subroutine takes two (usually) comma-delimited string lists, interprets it
  ! each as a list of elements, treating the first as keys and the second as
  ! a hash table, associative array or dictionary
  ! It returns the sub-string from the hash table corresponding to the key
  ! If the key is not found in the array of keys, it returns the delimiter
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-delimited lists

  ! If countEmpty is TRUE, consecutive delimiters, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single delimiter
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! If TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the delimiter may supplied, in case it isn't comma

  ! Basic premise: Use StringElementNum on key in keyList to find index
  ! Use this index to GetStringElement from HashList

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  
  SUBROUTINE GetStringHashElement(keyList, hashList, key, outElement, &
  & countEmpty, inDelim)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: keyList
    CHARACTER (LEN=*), INTENT(IN)   :: hashList
    CHARACTER (LEN=*), INTENT(IN)   :: key
    CHARACTER (LEN=*), INTENT(OUT)  :: outElement
    LOGICAL, INTENT(IN)   :: countEmpty
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    ! Local variables
	INTEGER(i4) :: elem
    CHARACTER (LEN=1)                          :: Delim
    CHARACTER (LEN=1), PARAMETER               :: COMMA = ','

    ! Executable code

    IF(PRESENT(inDelim)) THEN
	     Delim = inDelim
	 ELSE
	     Delim = COMMA
	 ENDIF

	elem = StringElementNum(keyList, key, countEmpty, inDelim)
	IF(elem <= 0) THEN
		outElement = Delim
	ELSE
		CALL GetStringElement(hashList, outElement, elem, countEmpty, inDelim)
	ENDIF

  END SUBROUTINE GetStringHashElement

  ! ---------------------------------------------  GetUniqueStrings  -----

  ! This subroutine takes an array of strings and returns another containing
  ! only the unique entries. The resulting array is supplied by the caller
  ! Some checking is done to make sure it's appropriate

  SUBROUTINE GetUniqueStrings(inList,outList,noUnique)
    ! Dummy arguments
    CHARACTER (LEN=*), DIMENSION(:) :: inList
    CHARACTER (LEN=*), DIMENSION(:) :: outList
    INTEGER :: noUnique ! Number of unique entries

    ! Local variables
    INTEGER :: i,j           ! Loop counters
    LOGICAL, DIMENSION(:), ALLOCATABLE :: duplicate ! Set if already found
    INTEGER :: status        ! Status from allocate

    INTEGER :: inSize

    ! Executable code, setup arrays

    inSize=SIZE(inList)
    ALLOCATE (duplicate(inSize), STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"duplicate")
    DO i = 1, inSize
       duplicate(i)=.FALSE.
    END DO

    ! Go through and find duplicates

    DO i = 1, inSize-1 ! Don't bother with last one
       IF (.NOT. duplicate(i)) THEN
          DO j = i+1, inSize
             IF (inList(j)==inList(i)) duplicate(j)=.TRUE.
          END DO
       END IF
    END DO

    ! Count how many unique ones there are

    noUnique=0
    DO i = 1, inSize
       IF (.NOT. duplicate(i)) noUnique=noUnique+1
    END DO

    IF (noUnique>SIZE(outList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "outList too small")
    IF (LEN(outList)<LEN(inList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "outList strings to small")
    outList=""

    j=1
    DO i = 1, noUnique
       UniqueHuntLoop: DO
          IF (.NOT. duplicate(j)) EXIT UniqueHuntLoop
          j=j+1
       END DO UniqueHuntLoop
       outList(i)=inList(j)
       j=j+1
    END DO

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

    ! Lenient wrt non-compilant formats:
    ! ignores chars in front of 'hh' and a terminal,
    ! non-numerical char: e.g., '2000-01-01T00:00:00.000000Z'
    ! will be treated the same as '00:00:00.0000000'

    ! If given optional arg strict, not lenient
    ! i.e., non-complient str always returns non-zero ErrTyp
    
    ! If given optional arg separator, uses separator as field separator
    
    ! Useful to allow an added way to input time
    !--------Argument--------!
    character(len=*),intent(in) :: str
    real(r8) :: value
    integer, intent(out) :: ErrTyp
    character(len=1),intent(in), optional :: separator
    logical,intent(in), optional :: strict
    !----------Local vars----------!
    character(len=1), parameter :: colon=':'
    integer, parameter :: INVALIDHHMMSSSTRING = 1
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
    logical, parameter :: DeeBUG = .false.
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
   
   call GetStringElement(str, hh, 1, countEmpty=.true., inDelim=myColon)
   call GetStringElement(str, mm, 2, countEmpty=.true., inDelim=myColon)
   call GetStringElement(str, ss, 3, countEmpty=.true., inDelim=myColon)
   
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

   if(DeeBUG) then
      print *, 'hh ', hh
      print *, 'mm ', mm
      print *, 'ss ', trim(ss)
   endif
   do i=1, len_trim(ss)
      if( .not. (index(digits, ss(i:i)) > 0 .or. ss(i:i) == '.') ) return
   enddo

   if(DeeBUG) then
      print *, 'ss passed first test'
   endif
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
   if(DeeBUG) then
      print *, 'mm passed first test'
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
   if(DeeBUG) then
      print *, 'hh value conversion error ', ErrTyp
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

   if(DeeBUG) then
      print *, 'mm value conversion error ', ErrTyp
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
   
   if(DeeBUG) then
      print *, 'ss value conversion error ', ErrTyp
   endif
   if(ErrTyp /= 0) then
      return
   elseif(value < 0. .or. value > 60.) then
      ErrTyp=INVALIDHHMMSSSTRING
      return
   endif
   
   value = value + 60*(mvalue + 60*hvalue)

  end Function hhmmss_value

  ! ------------------------------------  LinearSearchStringArray  -----

  ! This routine does a simple linear search for a string in a list.
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
  FUNCTION LowerCase (str) RESULT (outstr)
    ! takes A-Z  and replaces with a-z
    ! leaving other chars alone
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    INTEGER :: i, icode
    INTEGER, parameter :: offset=IACHAR("a")-IACHAR("A")
    !----------Executable part----------!
    outstr=str

    DO i = 1, LEN(str)
       icode=IACHAR(outstr(i:i))
       IF ( icode >=IACHAR("A") .AND. icode <= IACHAR("Z")) THEN
          outstr(i:i)=achar(icode+offset)
       END IF
    END DO

  END FUNCTION LowerCase

  ! ---------------------------------------------  NumStringElements  -----

  ! This function takes a (usually) comma-delimited string list, interprets it
  ! as a list of individual elements and returns the
  ! number of elements
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-delimited lists
  !
  ! If countEmpty is TRUE, consecutive delimiters, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single delimiter
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the delimiter may supplied, in case it isn't comma

  ! See also GetStringElement

  FUNCTION NumStringElements(inList, countEmpty, inDelim) RESULT (nElements)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)             :: inList
    LOGICAL, INTENT(IN)                       :: countEmpty
	 INTEGER                                   :: nElements
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    ! Local variables
    INTEGER :: i           ! Loop counters
	 LOGICAL :: lastWasNotDelim

    CHARACTER (LEN=1)               :: Delim
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    ! Executable code

    IF(PRESENT(inDelim)) THEN
	     Delim = inDelim
	 ELSE
	     Delim = COMMA
	 ENDIF

	! Count the number of delimiters
	! nElements-1 = number of delimiters
	IF(LEN_TRIM(inList) <= 0) THEN
		nElements=0
		RETURN
	ENDIF
	
	lastWasNotDelim = .FALSE.
	nElements=1
	DO i=1, LEN_TRIM(inList)
		IF(inList(i:i) == Delim) THEN
			IF(countEmpty .OR. lastWasNotDelim) THEN
				nElements = nElements+1
			ENDIF
			lastWasNotDelim = .FALSE.
		ELSE
			lastWasNotDelim = .TRUE.
		ENDIF
	ENDDO

  END FUNCTION NumStringElements

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

  ! --------------------------------------------------  ReverseList  -----
  FUNCTION ReverseList (str, inDelim) RESULT (outstr)
    ! takes a string list, usually comma-delimited,
	 ! and returns one with elements in reversed order

	 ! E.g., given "alpha, beta, gamma" => "gamma, beta, alpha"

	 ! Limitation:
	 ! No element may be longer than MAXWORDLENGTH
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    !----------Local vars----------!
    INTEGER(i4) :: i, istr, irev, elem, iBuf
    INTEGER, PARAMETER :: MAXWORDLENGTH=80
    CHARACTER (LEN=1)               :: Delim
    CHARACTER (LEN=1), PARAMETER    :: COMMA = ','
    CHARACTER (LEN=1), DIMENSION(:), ALLOCATABLE :: charBuf
    CHARACTER (LEN=MAXWORDLENGTH) :: word
! Treat consecutive delimiters as if enclosing an empty element
	LOGICAL, PARAMETER :: countEmpty = .TRUE.    

    !----------Executable part----------!
    IF(PRESENT(inDelim)) THEN
	     Delim = inDelim
	 ELSE
	     Delim = COMMA
	 ENDIF

!  Special case--only one element of str
    outstr = str
    IF(LEN(str) == 1 .OR. INDEX(str, Delim) == 0) RETURN
	 
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
		CALL GetStringElement(str, word, elem, countEmpty, Delim)
		IF(word == Delim) THEN
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
			charBuf(iBuf) = Delim
			elem = elem+1
		ENDIF
	ENDDO
	
	IF(charBuf(iBuf) == Delim) THEN
		iBuf = iBuf-1
	ENDIF
	
	DO i=1, iBuf
		irev = iBuf - i + 1
		outstr(irev:irev) = charBuf(i)
	ENDDO

	DEALLOCATE(charBuf)

  END FUNCTION ReverseList

  ! -------------------------------------------------  SplitWords  -----

  ! This subroutine is based on my IDL one of the same name.
  ! A line of input is split into sets of words.  There are two ways in which
  ! this can be invoked.  Typically it is split into `first' and 'rest'
  ! However, if the threeway option is given it is split to first, rest and
  ! last.

  ! Note that there is a slight subtlety here, spaces are treated specially
  ! due to the use of TRIM.  Thus while two commas in a row would count as
  ! two delimters, two spaces would count as one. Also if , is the delimeter
  ! then ,<space> counts as complete delimiter.

  SUBROUTINE SplitWords(line,first,rest,last,&
       & threeWay,delimiter)

    ! Dummy arguments

    CHARACTER (LEN=*), INTENT(IN) :: line
    CHARACTER (LEN=*), INTENT(OUT) :: first
    CHARACTER (LEN=*), INTENT(OUT) :: rest
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: last

    LOGICAL, INTENT(IN), OPTIONAL :: threeWay
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: delimiter

    ! Local variables

    CHARACTER (LEN=1) :: useDelimiter
    LOGICAL :: useThreeWay
    CHARACTER (LEN=LEN(line)) useLine ! Line with leading spaces removed

    INTEGER :: firstDelimiterPos,lastDelimiterPos,trimmedLen

    ! Executable code

    useLine=ADJUSTL(line)
    trimmedLen=LEN_TRIM(useLine)

    IF (PRESENT(delimiter)) THEN
       useDelimiter=delimiter
    ELSE
       useDelimiter=","
    END IF

    IF (PRESENT(threeWay)) THEN
       useThreeWay=threeWay 
    ELSE 
       useThreeWay=.FALSE.
    END IF

    ! Clear some results by default

    IF (PRESENT(last)) last=""
    rest=""

    ! Find the first delimiter

    firstDelimiterPos=INDEX(useLine,useDelimiter)

    IF (firstDelimiterPos == 0) THEN
       first=useLine
    ELSE
       first=useLine(1:firstDelimiterPos-1)
       IF (useThreeWay) THEN
          ! In three way mode, find the last delimiter
          lastDelimiterPos=INDEX(TRIM(useLine),useDelimiter,back=.TRUE.)
          IF (PRESENT(last) .AND. &
               & lastDelimiterPos /= trimmedLen) THEN
             last=TRIM(useLine(lastDelimiterPos+1:))
          END IF
          IF (firstDelimiterPos+1 <= lastDelimiterPos-1) THEN
             rest=TRIM(useLine(firstDelimiterPos+1:lastDelimiterPos-1))
          END IF
       ELSE
          IF (firstDelimiterPos /= trimmedLen) THEN
             rest=TRIM(useLine(firstDelimiterPos+1:))
          END IF
       END IF
    END IF

  END SUBROUTINE SplitWords
       
  ! ---------------------------------------------  StringElementNum  -----

  ! This function takes a (usually) comma-delimited string list, interprets it
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
  ! comma-delimited lists
  
  ! It will be the immediate precursor function in a hash table
  ! == aka associative array == aka dictionary
  !
  ! If countEmpty is TRUE, consecutive delimiters, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single delimiter
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the delimiter may supplied, in case it isn't comma

  ! See also GetStringElement, NumStringElements

  FUNCTION StringElementNum(inList, element, countEmpty, inDelim) RESULT (elem)
    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)             :: inList
    CHARACTER (LEN=*), INTENT(IN)             :: element
    LOGICAL, INTENT(IN)                       :: countEmpty
	 INTEGER                                   :: elem
    CHARACTER (LEN=1), OPTIONAL, INTENT(IN)   :: inDelim

    ! Local variables
    INTEGER :: nElements
    INTEGER , PARAMETER :: MAXELEMENTLENGTH = 80

    CHARACTER (LEN=MAXELEMENTLENGTH)           :: listElement
!    CHARACTER (LEN=1)                          :: Delim
!    CHARACTER (LEN=1), PARAMETER               :: COMMA = ','
    ! Executable code

	nElements = NumStringElements(inList, countEmpty, inDelim)
	
	IF(nElements <= 0) THEN
		elem = 0
		RETURN
	ENDIF

	! Check for matches--snipping off any leading blanks
	DO elem=1, nElements
		CALL GetStringElement(inList, listElement, elem, countEmpty, inDelim)
		IF(adjustl(listElement) == adjustl(element)) RETURN
	ENDDO
	
	elem = 0

  END FUNCTION StringElementNum

  ! ------------------------------------------------  unquote  -----
  Function unquote(str, quotes, cquotes, strict) result (outstr)
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
    
    ! If given optional arg quotes, removes only surrounding pair:
    ! quotes[i:i] for each i=1..len[quotes]
    ! E.g., given /a\ regexp/ with quotes='/' returns
    !    a\ regexp
    
    ! If given optional args quotes & cquotes, removes only surrounding pair:
    ! quotes[i:i] on the left, cquotes[i:i] on the right, i=1..len[quotes]
    ! E.g., given [a particle] with quotes='[' cquotes=']' returns
    !    a particle
    ! (For this case, strict matching is always on)
    
    ! Useful because the parser will return quote-delimited strings if that's
    ! how they appear in the lcf
    
    ! Calling get_string with "strip=.true." renders this unnecessary.
    ! However, you might find another use for it, especially with
    ! feature of being able to trim other, user-supplied detritus:
    ! e.g., braces, parentheses, extraneous delimiters
    !--------Argument--------!
    character(len=*),intent(in) :: str
    character(len=len(str)) :: outstr
    character(len=*),intent(in), optional :: quotes
    character(len=*),intent(in), optional :: cquotes
    logical,intent(in), optional :: strict
    !----------Local vars----------!
    character(len=1), parameter :: sq=''''
    character(len=1), parameter :: dq='"'
    integer :: first, last, ult, prim
    character(len=1) :: quote, cquote
    integer :: i
    logical :: mystrict
    !----------Executable part----------!

   ult = len_trim(str)    ! Position of last non-blank char
   prim = ult - len_trim(adjustl(str)) + 1    ! Position of 1st non-blank char
      
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

   if(last >= first) then
       outstr=str(first:last)
   else
       outstr=str
   endif
      
  end Function unquote

!=============================================================================
END MODULE MLSStrings
!=============================================================================

! $Log$
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
