!=================================
PROGRAM ExtractSubStringTest ! tests MLSStrings module
!=================================

   USE MLSStrings , ONLY: Capitalize, ExtractSubString, GetStringElement, &
     & NumStringElements

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the ExtractSubString subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering comma-separated lists of words; a blank line terminates

! Variables

   INTEGER, PARAMETER                :: MAXLISTLENGTH=80
   INTEGER, PARAMETER                :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH)     :: theString
	CHARACTER (LEN=MAXLISTLENGTH)     :: anotherString
	CHARACTER (LEN=4)                 :: sub1, sub2
   INTEGER                           :: listLength, NumWords, wordNum
	LOGICAL, PARAMETER                :: kountEmpty = .TRUE.
	LOGICAL, PARAMETER                :: CaseSensitive = .false.

   ! 1st--just a simple test of how option
   theString = 'abcabc123defdef'
   sub1 = 'abc'
   sub2 = 'def'
   print *, 'theString ', trim(theString)
   print *, 'sub1 ', trim(sub1)
   print *, 'sub2 ', trim(sub2)
   call ExtractSubString(theString, anotherString, sub1, sub2)
   print *, '(after extracting with default how', trim(anotherString)
   call ExtractSubString(theString, anotherString, sub1, sub2, how='stingy')
   print *, '(after extracting with stingy how', trim(anotherString)
   call ExtractSubString(theString, anotherString, sub1, sub2, how='greedy')
   print *, '(after extracting with greedy how', trim(anotherString)
   ! 1st--just a simple test of no_trim option
   theString = 'abc abc 123def def'
   sub1 = 'abc '
   sub2 = 'def '
   print *, 'theString ', trim(theString)
   print *, 'sub1 "', sub1, '"'
   print *, 'sub2 "', sub2, '"'
   call ExtractSubString(theString, anotherString, sub1, sub2, &
     & no_trim=.true.)
   print *, '(after extracting with default how', trim(anotherString)
   call ExtractSubString(theString, anotherString, sub1, sub2, how='stingy', &
     & no_trim=.true.)
   print *, '(after extracting with stingy how', trim(anotherString)
   call ExtractSubString(theString, anotherString, sub1, sub2, how='greedy', &
     & no_trim=.true.)
   print *, '(after extracting with greedy how', trim(anotherString)
	DO

      ! Prompt for input

 	   print *, 'Enter sub1'
		read(*, '(A)') sub1
      print *, 'You entered ', trim(sub1)
		IF(sub1 == ' ' .or. Capitalize(sub1) == 'STOP') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
 	   print *, 'Enter sub2'
		read(*, '(A)') sub2
      print *, 'You entered ', trim(sub2)
	
      do
 	     print *, 'Enter string to be processed'
		  read(*, '(A)') theString
        print *, 'You entered ', trim(theString)
		  IF(theString == ' ') THEN
   	     CYCLE
        elseif(Capitalize(theString) == 'STOP') then
           stop
		  ENDIF

        call ExtractSubString(theString, anotherString, sub1, sub2)
		  print *, '(after default extracting) theString = ', anotherString
        call ExtractSubString(theString, anotherString, sub1, sub2, how='stingy')
		  print *, '(after stingy extracting) theString = ', anotherString
        call ExtractSubString(theString, anotherString, sub1, sub2, how='greedy')
		  print *, '(after greedy extracting) theString = ', anotherString
		ENDDO
	ENDDO
!==================
END PROGRAM ExtractSubStringTest
!==================

! $Log$
