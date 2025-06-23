!=================================
PROGRAM WordListTest ! tests MLSStrings module
!=================================

   USE MLSStrings , ONLY: GetStringElement, NumStringElements, Reverse, &
	& ReverseList, StringElementNum, SortList, SortArray, Capitalize

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the GetStringElement subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering comma-separated lists of words; a blank line terminates

! Variables

   INTEGER, PARAMETER                :: MAXLISTLENGTH=80
   INTEGER, PARAMETER                :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH)     :: theList, sortedList
	CHARACTER (LEN=MAXWORDLENGTH)     :: aWord
	CHARACTER (LEN=MAXWORDLENGTH), DIMENSION(MAXLISTLENGTH) &
    &                                :: words, sortedWords
   INTEGER                           :: listLength, NumWords, wordNum
	LOGICAL, PARAMETER                :: kountEmpty = .TRUE.
	LOGICAL, PARAMETER                :: CaseSensitive = .false.
	LOGICAL, PARAMETER                :: ignoreLeadingSpaces = .false.
	LOGICAL, PARAMETER                :: shorterFirst = .true.
	CHARACTER (LEN=1), PARAMETER      :: listArray = 'l'  ! Sort as list or array?
	CHARACTER (LEN=1), PARAMETER      :: leftRight = 'r'  ! How to interpret iArray
   integer, dimension(MAXLISTLENGTH) :: iArray

	DO

! Prompt for input

 	  print *, 'Enter comma-separated list of words'
		read(*, '(A80)') theList
      print *, 'You entered ', trim(theList)
		IF(theList == ' ' .or. Capitalize(theList(1:4)) == 'STOP') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
	
!	Process theList
		NumWords = NumStringElements(theList, kountEmpty)
		listLength = LEN_TRIM(theList)
	
   	print *, 'Number of chars read = ', listLength
   	print *, 'Number of words read = ', NumWords
		IF(NumWords <= 0) CYCLE
	
		DO wordNum=1, NumWords
			CALL GetStringElement(theList, aWord, wordNum, kountEmpty)
         words(wordNum) = aWord
! Print the output
			print *, 'Word = ', aWord, ' ;  droW = ', Reverse(aWord)
 	  		print *, 'Enter comma-separated list of words'
		ENDDO
		print *, 'tsiLeht = ', ReverseList(theList)
		if( listArray == 'l' ) then
        call SortList(theList, iArray, CaseSensitive, kountEmpty, &
         & ignoreLeadingSpaces=ignoreLeadingSpaces, inDelim=',', &
         & sortedList=theList, leftRight=leftRight)
        print *, 'Sorted list'
        print *, theList
      else
        call SortArray(words(1:NumWords), iArray, CaseSensitive, &
         & sortedArray=words, shorterFirst=shorterFirst, leftRight=leftRight)
        print *, 'Sorted Array'
        print *, (trim(words(wordNum)) //' ', wordNum=1, NumWords)
      endif
      print *, 'Sorted word order (l/r?) ', leftRight
      print *, (iArray(wordNum), wordNum=1, NumWords)
		DO
			read(*, '(A16)') aWord
			IF(aWord == ' ' .or. Capitalize(aWord(1:3)) == 'END') then
           EXIT
         elseif(Capitalize(aWord(1:4)) == 'STOP') then
           stop
         endif
   		print *, 'Element number of word in list = ', &
			& StringElementNum(theList, aWord, kountEmpty)
		ENDDO
	ENDDO
!==================
END PROGRAM WordListTest
!==================

!# $Log$
