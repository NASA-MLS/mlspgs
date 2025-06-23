!=================================
PROGRAM streqTest ! tests MLSStrings module
!=================================

   USE MLSStrings, ONLY: indexes, isRepeat, lowercase, NCopies, streq
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the streq and other new functions.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering comma-separated lists of words; a blank line terminates

! Variables

   INTEGER, PARAMETER                :: MAXARRAYLENGTH=100
   INTEGER, PARAMETER                :: MAXLISTLENGTH=800
	CHARACTER (LEN=MAXLISTLENGTH)     :: theString
	CHARACTER (LEN=MAXLISTLENGTH)     :: pattern
	CHARACTER (LEN=8            )     :: options
	CHARACTER (LEN=8), dimension(MAXARRAYLENGTH)     :: substrings
	CHARACTER (LEN=1            )     :: test
   integer, dimension(MAXARRAYLENGTH):: ints
   logical, dimension(MAXARRAYLENGTH):: logs
   logical                           :: isit
   integer                           :: highfit
   integer                           :: nsubs
   logical, parameter :: countEmpty=.true.

 	 print *, 'What test do you wish to run?'
 	 print *, 'streq (s), indexes (x), ncopies (n), isRepeat (r)'
	 read(*, '(A)') test
    select case (lowercase(test))
    case ('s')

 	   print *, 'Enter pattern to be matched against'
		read(*, '(A)') pattern
      print *, 'You entered ', trim(pattern)
 	   print *, 'Enter options (e.g., -wcf)'
		read(*, '(A)') options
      print *, 'You entered ', trim(options)
      do
 	     print *, 'Enter string to be processed'
		  read(*, '(A)') theString
		  IF(theString == ' ') exit
        print *, 'You entered ', trim(theString)
        isit = streq(theString, pattern, options)
        print *, 'matches? ', isit
		ENDDO
    case ('x')

 	   print *, 'Enter string'
		read(*, '(A)') theString
      print *, 'You entered ', trim(theString)
 	   ! print *, 'Enter mode (left, right, first,last, wrap)'
		! read(*, '(A)') options
      ! print *, 'You entered ', trim(options)
      nsubs = 0
      do
 	     print *, 'Enter substring to be processed (blank after last one)'
        nsubs = nsubs+1
		  read(*, '(A)') substrings(nsubs)
		  IF(substrings(nsubs) == ' ') exit
		ENDDO
      nsubs = nsubs - 1
      if ( nsubs < 1 ) stop
      ! ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode=options)
      ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode='first')
      print *, 'indexes (mode=first): ', ints(1:nsubs)
      ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode='last')
      print *, 'indexes (mode=last): ', ints(1:nsubs)
      ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode='left')
      print *, 'indexes (mode=left): ', ints(1:nsubs)
      ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode='right')
      print *, 'indexes (mode=right): ', ints(1:nsubs)
      ints(1:nsubs) = indexes(theString, substrings(1:nsubs), mode='wrap')
      print *, 'indexes (mode=wrap): ', ints(1:nsubs)
    case ('n')

 	   print *, 'Enter string'
		read(*, '(A)') theString
      print *, 'You entered ', trim(theString)
 	   print *, 'Enter substring to be processed'
		read(*, '(A)') substrings(1)
		IF(substrings(1) == ' ') then
        substrings(1) = theString(1:1)
      endif
        nsubs = 1
 	   ! print *, 'May substrings overlap?'
		! read(*, '(L)') isit
      nsubs = ncopies(theString, trim(substrings(1)), overlap=.false.)
      print *, 'ncopies(non-overlapping): ', nsubs
      nsubs = ncopies(theString, trim(substrings(1)), overlap=.true.)
      print *, 'ncopies(overlapping): ', nsubs
    case ('r')

 	   print *, 'Enter string'
		read(*, '(A)') theString
      print *, 'You entered ', trim(theString)
 	   print *, 'Enter substring to be processed (or blank if none)'
		read(*, '(A)') substrings(1)
		IF(substrings(1) == ' ') then
        nsubs = 0
      else
        nsubs = 1
      endif
      if ( nsubs < 1 ) then
        isit = isRepeat(theString)
      else
        isit = isRepeat(theString, substrings(1))
      endif
      print *, 'is repeat?: ', isit
    case default
      print *, 'Unrecognized test: ', test
    end select
!==================
END PROGRAM streqTest
!==================

! $Log$
