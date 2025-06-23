!=================================
PROGRAM ReplaceSubStringTest ! tests MLSStrings module
!=================================

   USE MLSStrings , ONLY: GetStringElement, NumStringElements, Reverse, &
	& ReverseList, StringElementNum, SortList, SortArray, Capitalize, &
   & ReplaceSubString

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the ReplaceSubString subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering comma-separated lists of words; a blank line terminates

! Variables

   INTEGER, PARAMETER                :: MAXLISTLENGTH=80
   INTEGER, PARAMETER                :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH)     :: thesub1List, thesub2List, theString
	CHARACTER (LEN=MAXLISTLENGTH)     :: anotherString
	CHARACTER (LEN=MAXWORDLENGTH)     :: sub1, sub2
   INTEGER                           :: listLength, NumWords, wordNum
	LOGICAL, PARAMETER                :: kountEmpty = .TRUE.
	LOGICAL, PARAMETER                :: CaseSensitive = .false.

   ! 1st--just a simple test of no_trim option
   theString = 'Thearg1 the arg1  the  arg1the  arg1   thearg1    the    arg1'
   print *, 'theString ', trim(theString)
   call ReplaceSubString(theString, anotherString, ' arg1', ' arg2', &
     & which='all', no_trim=.true.)
   print *, '(after replacing ( arg1) with ( arg2)', trim(anotherString)
   call ReplaceSubString(theString, anotherString, 'arg1 ', 'arg2 ', &
     & which='all', no_trim=.true.)
   print *, '(after replacing (arg1 ) with (arg2 )', trim(anotherString)
   call ReplaceSubString(theString, anotherString, ' arg1 ', ' arg2 ', &
     & which='all', no_trim=.true.)
   print *, '(after replacing ( arg1 ) with ( arg2 )', trim(anotherString)
   call ReplaceSubString(theString, anotherString, 'arg1  ', 'arg2  ', &
     & which='all', no_trim=.true.)
   print *, '(after replacing (arg1  ) with (arg2  )', trim(anotherString)
   call ReplaceSubString(theString, anotherString, '  arg1', '  arg2', &
     & which='all', no_trim=.true.)
   print *, '(after replacing (  arg1) with (  arg2)', trim(anotherString)
   call ReplaceSubString(theString, anotherString, '   arg1', '   arg2', &
     & which='all', no_trim=.true.)
   print *, '(after replacing (   arg1) with (   arg2)', trim(anotherString)
   call ReplaceSubString(theString, anotherString, '   arg1 ', '   arg2 ', &
     & which='all', no_trim=.true.)
   print *, '(after replacing (   arg1 ) with (   arg2 )', trim(anotherString)
	DO

      ! Prompt for input

 	   print *, 'Enter comma-separated list of sub1 (originals)'
		read(*, '(A80)') thesub1List
      print *, 'You entered ', trim(thesub1List)
		IF(thesub1List == ' ' .or. Capitalize(thesub1List(1:4)) == 'STOP') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
 	   print *, 'Enter comma-separated list of sub2 (replacements)'
		read(*, '(A80)') thesub2List
      print *, 'You entered ', trim(thesub2List)
	
      !	Process theList
		NumWords = NumStringElements(thesub1List, kountEmpty)
		listLength = LEN_TRIM(thesub1List)
   	print *, '----- sub1 list (originals) ----'
   	print *, 'Number of chars read = ', listLength
   	print *, 'Number of words read = ', NumWords
		IF(NumWords <= 0) CYCLE
		NumWords = NumStringElements(thesub1List, kountEmpty)

		listLength = LEN_TRIM(thesub2List)
   	print *, '----- sub2 list (replacements) ----'
   	print *, 'Number of chars read = ', listLength
   	print *, 'Number of words read = ', NumWords
		IF(NumWords <= 0) CYCLE
	
      do
 	     print *, 'Enter string to be processed'
		  read(*, '(A80)') theString
        print *, 'You entered ', trim(theString)
		  IF(theString == ' ') THEN
   	     CYCLE
        elseif(Capitalize(theString) == 'STOP') then
           stop
		  ENDIF

		  DO wordNum=1, NumWords
		     CALL GetStringElement(thesub1List, sub1, wordNum, kountEmpty)
		     CALL GetStringElement(thesub2List, sub2, wordNum, kountEmpty)
           print *, 'sub1: ', trim(sub1), ' sub2: ', trim(sub2)
           call ReplaceSubString(theString, theString, sub1, sub2)
           ! theString = anotherString
		  ENDDO
		  print *, '(after processing) theString = ', theString
		ENDDO
	ENDDO
!==================
END PROGRAM ReplaceSubStringTest
!==================

!# $Log$
