!=================================
PROGRAM timetest ! tests subroutine
!=================================

   USE MLSCommon , ONLY: r8
   USE MLSStrings , ONLY: hhmmss_value, RemoveElemFromList, unquote, &
     & utc_to_yyyymmdd, yyyymmdd_to_dai
   use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
   use time_m, only: SET_STARTTIME, TIME_NOW, USE_WALL_CLOCK

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

   ! The following arrys contains the maximum permissible day for each month
   ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
   integer, dimension(-1:12), parameter :: DAYMAXLY = (/ &
     & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
     & /)
   integer, dimension(-1:12), parameter :: DAYMAXNY = (/ &
     & 365, 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
     & /)

   INTEGER, PARAMETER :: MAXLISTLENGTH=8
   INTEGER, PARAMETER :: ntimesteps=50
   INTEGER, PARAMETER :: delay=500000
	CHARACTER (LEN=MAXLISTLENGTH) :: date
	LOGICAL, PARAMETER :: strict = .false.
   real :: t1, t2
   integer :: mydaysoff
   integer :: timesteps
   logical, parameter :: total_times = .true.
   integer, dimension(8) :: values
   integer, dimension(8) :: dvalues
   character(len=80) :: switches, tempSwitches
   character(len=8) :: badSwitch
   character(len=2) :: quotes
   quotes = char(34) // char(39)
	DO

! Prompt for input

 	  print *, 'Enter a string yyyymmdd'
		read(*, '(A8)') date
			IF(date == ' ') THEN
   		print *, 'GAME OVER'
			exit
		ENDIF
	
!	Process date
    call yyyymmdd_to_dai(date, mydaysoff)

    print *, 'The days past Jan 2001 your string = ', mydaysoff
	ENDDO

	DO

! Prompt for input

 	  print *, 'Enter a string ,-separated switches'
		read(*, '(A80)') switches
			IF(switches == ' ') THEN
   		print *, 'GAME OVER'
			exit
		ENDIF
      tempSwitches = unquote(switches, quotes=quotes, stripany=.true.)
      switches=tempswitches
	
     print *, '(After unquoting)'
     print *, trim(Switches)
 	  print *, 'Now the switch to remove'
		read(*, '(A8)') badSwitch
	
!	Remove badSwitch
    call RemoveElemFromList(switches, tempSwitches, trim(badSwitch))

    print *, trim(tempSwitches)
	ENDDO

! Now test time_now function
  USE_WALL_CLOCK = .true.
  call time_now(t1)
  call output ( "Starting time  = " )
  call output ( dble(t1), advance = 'yes' )
  do timesteps=1, ntimesteps
    call usleep(delay)
    call saytime
    call time_now(t1)
  enddo
  values =  (/2003, 12, 5, 0, 23, 50, 0, 0/)
  dvalues = (/0,     0, 0, 0, 0,  5,  0, 0/)
  call set_startTime( values )
  call time_now(t1, values)
  do timesteps=1, ntimesteps
    call advance(values, dvalues)
    ! call output('daysoff, startdaysoff: ', advance='no')
    ! call output((/daysoff, startdaysoff/), advance='no')
    ! call output(' values: ', advance='no')
    ! call output(values, advance='yes')
    call saytime(values)
    call time_now(t1, values)
  enddo
contains
  subroutine advance ( values, dvalues )
    integer, dimension(8), intent(inout) :: values
    integer, dimension(8), intent(in)    :: dvalues
    integer, dimension(8), parameter     :: lolimits = &
      & (/1,     1,  1, 0,  0, 0,   0,  0/)
    integer, dimension(8), parameter     :: uplimits = &
      & (/2000, 12, -1, 0, 23, 59, 59, 999/)
    integer :: hindx
    integer :: indx
    integer :: limit
    integer, dimension(-1:12) :: DAYMAX
    if ( leapyear(values(1)) ) then
      DAYMAX = DAYMAXLY
    else
      DAYMAX = DAYMAXNY
    endif
    values = values + dvalues
    do indx = 8, 2, -1
      limit = uplimits(indx)
      if ( limit == 0 ) cycle
      if ( uplimits(indx-1) > 0 ) then
        hindx = indx-1
      else
        hindx = indx-2
      endif
      if ( limit < 0 ) then
        limit = DAYMAX(values(2))
      endif
      if ( values(indx) > limit ) then
        values(hindx) = values(hindx) + 1
        values(indx) = values(indx) - (limit+1) + lolimits(indx)
      endif
    enddo
  end subroutine advance
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

  subroutine SayTime ( values )
    integer, dimension(8), intent(in), optional :: values
    call time_now ( t2, values )
    if ( total_times ) then
      call output ( "Total time = " )
      call output ( dble(t2), advance = 'no' )
      call blanks ( 4, advance = 'no' )
    endif
    call output ( "increment  = " )
    call output ( dble(t2 - t1), advance = 'yes' )
  end subroutine SayTime
!==================
END PROGRAM timetest
!==================

! $Log$
