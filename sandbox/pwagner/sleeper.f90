!=================================
PROGRAM sleeper ! tests how date_and_time behaves crossing midnight
!=================================

   USE MLSCommon , ONLY: r8
   use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
   use time_m, only: SET_STARTTIME, TIME_NOW, USE_WALL_CLOCK

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program sleeps for 10 minutes, prints date_and_time data, then repeats

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
   INTEGER, PARAMETER :: ntimesteps=60  ! How many 10 minute intervals this will do that
   INTEGER, PARAMETER :: udelay=500000  ! corresponds to 1/2 second
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
! Now test time_now function
  USE_WALL_CLOCK = .true.
  call time_now(t1)
  call output ( "Starting time  = " )
  call output ( dble(t1), advance = 'yes' )
  
  do timesteps=1, ntimesteps
    call take_ten
    call saytime
    call time_now(t1)
  enddo

contains

  subroutine take_ten 
  ! Sleep for 10 minutes
  integer :: timesteps
  integer, parameter :: ntimesteps = 20*60
  do timesteps=1, ntimesteps
    call usleep(udelay)
  enddo
  end subroutine take_ten

  subroutine SayTime ( values )
    integer, dimension(8), intent(in), optional :: values
    integer, dimension(8) :: myvalues
    call time_now ( t2, values )
    if ( present(values) ) then
      call output ( "(values: " )
      call output ( values, advance = 'no' )
      call output ( ")" )
      call blanks ( 2, advance = 'no' )
    endif
    if ( total_times ) then
      call output ( "Total time = " )
      call output ( dble(t2), advance = 'no' )
      call blanks ( 4, advance = 'no' )
    endif
    call output ( "increment  = " )
    call output ( dble(t2 - t1), advance = 'yes' )
    call date_and_time(date=date, values=myvalues)
    call output ( "(values: " )
    call output ( myvalues, advance = 'no' )
    call output ( ")" )
    call blanks ( 2, advance = 'no' )
    call output ( 'date: ' // date, advance = 'yes' )
  end subroutine SayTime
!==================
END PROGRAM sleeper
!==================

! $Log$
