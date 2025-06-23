!=================================
PROGRAM hhmmss_value_test ! tests subroutine
!=================================

   USE MLSCommon , ONLY: r8, FileNameLen
   USE MLSStrings , ONLY: hhmmss_value, utc_to_yyyymmdd
   use SDPTOOLKIT, only: mls_utctotai

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the hhmmss_value subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering quote-surrounded strings of chars; a blank line terminates

! Variables

   INTEGER, PARAMETER :: MAXLISTLENGTH=80
   INTEGER, PARAMETER :: MAXWORDLENGTH=16
	CHARACTER (LEN=MAXLISTLENGTH) :: theString
	CHARACTER (LEN=MAXWORDLENGTH) :: quotes, cquotes
   character(LEN=FileNameLen) :: LeapSecFileName
	LOGICAL, PARAMETER :: strict = .false.
   real(r8) :: theValue
   integer :: ErrTyp

 	print *, 'Enter the leap second file name (with path)'
	read(*, '(A126)') LeapSecFileName
	print *, 'leap second file: ', trim(LeapSecFileName)
	DO

! Prompt for input

 	  print *, 'Enter a string hh:mm:ss'
		read(*, '(A80)') theString
			IF(theString == ' ') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
	
!	Process theString
    theValue = hhmmss_value(theString, ErrTyp, strict=strict)

      if(ErrTyp==0) then
        print *, 'The value of your string = ', theValue
      else
        print *, 'Sorry--can''t grok your string as hh:mm:ss'
      endif
! Convert to tai
      if ( LeapSecFileName /= ' ' ) then
        ErrTyp = mls_utctotai(trim(LeapSecFileName), trim(theString), &
        & theValue)
        print *, '           (tai) ', theValue
      endif
	ENDDO


!==================
END PROGRAM hhmmss_value_test
!==================

!# $Log$
!# Revision 1.1  2002/02/19 23:28:19  pwagner
!# First commit
!#
!#
