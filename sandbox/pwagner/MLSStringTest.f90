!=================================
PROGRAM MLSStringTest ! tests subroutine
!=================================

   USE MLSStrings , ONLY: GetStringElement, NumStringElements, Reverse, &
	& ReverseList, StringElementNum

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
! then enter "make -f MakeFC depends" followed by "make -f MakeFC"


! Then run it, entering comma-separated lists of words; a blank line terminates

! Variables

   INTEGER, PARAMETER :: MAXLISTLENGTH=80
   INTEGER, PARAMETER :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH) :: theList
	CHARACTER (LEN=MAXWORDLENGTH) :: aWord
   INTEGER :: listLength, NumWords, wordNum
	LOGICAL, PARAMETER :: kountEmpty = .TRUE.

	DO

! Prompt for input

 	  print *, 'Enter comma-separated list of words'
		read(*, '(A80)') theList
			IF(theList == ' ') THEN
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
! Print the output

			print *, 'Word = ', aWord, ' ;  droW = ', Reverse(aWord)

 	  		print *, 'Enter comma-separated list of words'
		ENDDO
		print *, 'tsiLeht = ', ReverseList(theList)
		
		DO
			read(*, '(A16)') aWord
			IF(aWord == ' ') EXIT
   		print *, 'Element number of word in list = ', &
			& StringElementNum(theList, aWord, kountEmpty)
		ENDDO
	ENDDO


!==================
END PROGRAM MLSStringTest
!==================

!# $Log$
!#
