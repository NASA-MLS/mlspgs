!=============================================================================
MODULE MLSStrings               ! Some low level string handling stuff
!=============================================================================

  USE MLSMessageModule

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: id,ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

!
! This module contains some low level string handling stuff
!

CONTAINS
  FUNCTION count_words(str) RESULT (no_of_words)
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
    ENDIF
    DO j=2,LEN(str)
       IF(str(j:j) /= " " .AND. str(j-1:j-1) == " ") THEN
          no_of_words=no_of_words+1
       ENDIF
    ENDDO
  END FUNCTION count_words

  ! This one converts a string to all upper case (taken from HCP routine
  ! of same name) (Except that HCP can spell capitalise, that is. Fnord.)

  FUNCTION Capitalize(str) RESULT (outstr)
    ! takes a-z and replaces with A-Z 
    ! leaving other chars alone
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    CHARACTER(LEN=LEN(STR))::capstr
    INTEGER::i,icode,offset
    !----------Executable part----------!
    capstr=str
    offset=ICHAR("A")-ICHAR("a")

    DO i=1,LEN(str)
       icode=ICHAR(capstr(i:i))
       IF ( icode >=ICHAR("a") .AND. icode <= ICHAR("z")) THEN
          capstr(i:i)=char(icode+offset)
       ENDIF
    ENDDO

    outstr=capstr
  END FUNCTION Capitalize

!=============================== lowercase ==========================
FUNCTION lowercase(str) RESULT (outstr)
!=============================== lowercase ==========================
    ! takes A-Z and replaces with a-z 
    ! leaving other chars alone
! (Obviously, a crude theft from the above)
    !--------Argument--------!
    CHARACTER (LEN=*), INTENT(IN) :: str
    CHARACTER (LEN=LEN(str)) :: outstr

    !----------Local vars----------!
    CHARACTER(LEN=LEN(STR))::capstr
    INTEGER::i,icode,offset
    !----------Executable part----------!
    capstr=str
    offset=ICHAR("a")-ICHAR("A")

    DO i=1,LEN(str)
       icode=ICHAR(capstr(i:i))
       IF ( icode >=ICHAR("A") .AND. icode <= ICHAR("Z")) THEN
          capstr(i:i)=char(icode+offset)
       ENDIF
    ENDDO

    outstr=capstr
END FUNCTION lowercase

  Function depunctuate(str) result (outstr)
    ! Function that removes punctuation and replaces with blanks
    ! Added by HCP. This depends on the native character set being 
    ! ASCII. 
    !--------Argument--------!
    character(len=*),intent(in)::str
    character(len=len(str))::outstr
    !----------Local vars----------!
    integer::i,icode
    !----------Executable part----------!
    outstr=str
    do i=1 ,len(str)
        icode=ichar(str(i:i))
        if(  (icode >= 33 .and. icode <= 47).or.&
             (icode >= 58 .and. icode <=64).or. &
             (icode >= 91 .and. icode <=96).or. &
             (icode >= 123)) then
           outstr(i:i)=" "
        endif
    enddo

  end Function depunctuate

  !---------------------------------------------------------------------------

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
    ENDIF

    IF (.NOT. PRESENT(continuationChar)) THEN
       useContinuationChar="$"
    ELSE
       useContinuationChar=continuationChar
    ENDIF

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
       ENDIF

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
       endif
       ! Concatenate this with what we have so far, make sure there's an extra
       ! space there though (not for first line though)

       ! If block inserted 5 Sept. 2000 by HCP to prevent an out-of-range 
       ! error when the input line has 0 length and you have bounds-checking
       ! turned on
       if(LEN_TRIM(inputLine) > 0) then 
          inputLine=inputLine(1:LEN_TRIM(inputLine)-gotContinuation)
       endif
       IF (firstLine) THEN
          fullLine=inputLine
          firstLine=.FALSE.
       ELSE
          fullLine=fullLine(1:LEN_TRIM(fullLine)+1)//inputLine
       ENDIF
       
       ! If we have a continuation mark, or a blank line then keep reading
       IF ((gotContinuation==0).AND.(LEN_TRIM(fullLine) /= 0)) EXIT readLoop
    END DO readLoop

    ! Do a final trim and exit

    fullLine=TRIM(ADJUSTL(fullLine))

    END SUBROUTINE ReadCompleteLineWithoutComments

    !---------------------------------------------------------------------------

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
      ENDIF

      IF (PRESENT(threeWay)) THEN
         useThreeWay=threeWay 
      ELSE 
         useThreeWay=.FALSE.
      ENDIF

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
            ENDIF
            IF (firstDelimiterPos+1 <= lastDelimiterPos-1) THEN
               rest=TRIM(useLine(firstDelimiterPos+1:lastDelimiterPos-1))
            ENDIF
         ELSE
            IF (firstDelimiterPos /= trimmedLen) THEN
               rest=TRIM(useLine(firstDelimiterPos+1:))
            ENDIF
         ENDIF
      ENDIF

     END SUBROUTINE SplitWords

     ! ------------------------------------------------------------------------
     
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
       DO i=1,inSize
          duplicate(i)=.FALSE.
       END DO

       ! Go through and find duplicates

       DO i=1,inSize-1 ! Don't bother with last one
          IF (.NOT. duplicate(i)) THEN
             DO j=i+1,inSize
                IF (inList(j)==inList(i)) duplicate(j)=.TRUE.
             END DO
          ENDIF
       END DO
       
       ! Count how many unique ones there are

       noUnique=0
       DO i=1,inSize
          IF (.NOT. duplicate(i)) noUnique=noUnique+1
       END DO
       
       IF (noUnique>SIZE(outList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "outList too small")
       IF (LEN(outList)<LEN(inList)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "outList strings to small")
       outList=""

       j=1
       DO i=1,noUnique
          UniqueHuntLoop: DO
             IF (.NOT. duplicate(j)) EXIT UniqueHuntLoop
             j=j+1
          END DO UniqueHuntLoop
          outList(i)=inList(j)
          j=j+1
       ENDDO

       DEALLOCATE(duplicate)
     END SUBROUTINE GetUniqueStrings

     ! ------------------------------------------------------------------------

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
       ENDIF

       IF (PRESENT(testSubstring)) THEN
          testForSubstring = testSubstring
       ELSE
          testForSubstring = .FALSE.
       ENDIF

       IF (PRESENT(listInString)) THEN
          testForList = listInString
       ELSE
          testForList = .FALSE.
       ENDIF

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
             ENDIF
             IF (found) THEN
                sindex = i
                EXIT linSearchStringInsens
             ENDIF
          END DO linSearchStringInsens
       ELSE
          linSearchStringSens: DO i = 1, SIZE(list)
             IF (testForSubstring) THEN
                IF (testForList) THEN
                   found = (INDEX (TRIM(string), TRIM(list(i))) /= 0)
                ELSE
                   found = (INDEX (TRIM(list(i)), TRIM(string)) /= 0)
                ENDIF
             ELSE
                found = (TRIM(list(i)) == TRIM(string))
             ENDIF
             IF (found) THEN
                sindex = i
                EXIT linSearchStringSens
             ENDIF
          END DO linSearchStringSens
       ENDIF

     END FUNCTION LinearSearchStringArray

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
       
!=============================================================================
END MODULE MLSStrings
!=============================================================================

! $Log$
! Revision 2.1  2000/11/30 00:24:49  pwagner
! lowercase function moved here from Fill
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.21  2000/09/05 10:59:32  pumphrey
! HCP Fixed an out-of-bounds error
!
