!=================================
PROGRAM UnquoteTest ! tests subroutine
!=================================

   USE MLSStrings , ONLY: Unquote

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the unquote subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering quote-surrounded strings of chars; a blank line terminates

! Variables

   INTEGER, PARAMETER :: MAXLISTLENGTH=80
   INTEGER, PARAMETER :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH) :: theString
	CHARACTER (LEN=MAXWORDLENGTH) :: quotes, cquotes
	LOGICAL, PARAMETER :: strict = .TRUE.

 	  print *, 'Enter the types of open "quotes" you will be using'
 	  print *, 'E.g., if (my string) then quotes = ('
 	  print *, '(leave empty if you wish to use default single && double-quotes)'
		read(*, '(A16)') quotes
 	  print *, 'Enter the types of closing "quotes" you will be using'
 	  print *, '(if the same as quotes, just leave blank)'
		read(*, '(A16)') cquotes
	DO

! Prompt for input

 	  print *, 'Enter a string surrounded by quotes'
		read(*, '(A80)') theString
			IF(theString == ' ') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
	
!	Process theString
    if( quotes == ' ') then
      theString = unquote(theString, strict=strict)
    elseif( cquotes == ' ') then
      theString = unquote(theString, quotes, strict=strict)
    else
      theString = unquote(theString, quotes, cquotes=cquotes, strict=strict)
    endif

      print *, 'Unquoted, your string = ' // trim(theString)

	ENDDO


!==================
END PROGRAM UnquoteTest
!==================

!# $Log$
!#
