! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Calendar

  ! Several procedures and defined operators to deal with time.
  ! Doesn't account for leap seconds.

  implicit NONE
  private

  public :: Caldat, Caldat_MJD_YMD, Caldat_MJD_YMDHMS, Caldat_T, Caldat_T_MJD
  public :: Caldat_YMD, Caldat_YMDHMS, Duration
  public :: Duration_Plus_Time, Duration_Plus_Time_I
  public :: Julday, MJD, T_Julday, T_MJD
  public :: Time_Minus_Duration, Time_Minus_Duration_I
  public :: Time_Plus_Duration, Time_Plus_Duration_I
  public :: Time_T, YMD_Julday, YMD_MJD, YMDHMS_MJD
  public :: Operator(+), Operator(-)

  integer, public, parameter :: TK = selected_real_kind(12) ! Kind for Julday

  type Time_T
    integer :: Year, Month, Day, Hours=0, Minutes=0
    real :: Seconds=0.0
  end type Time_T

  interface Caldat
    module procedure Caldat_YMD, Caldat_YMDHMS
  end interface

  interface Caldat_MJD
    module procedure Caldat_MJD_YMD, Caldat_MJD_YMDHMS
  end interface

  interface Julday
    module procedure T_Julday, YMD_Julday, YMDHMS_Julday
  end interface

  interface MJD
    module procedure T_MJD, YMD_MJD, YMDHMS_MJD
  end interface

  interface OPERATOR(-)
    module procedure Duration, Time_Minus_Duration, Time_Minus_Duration_I
  end interface

  interface OPERATOR(+)
    module procedure Duration_Plus_Time, Duration_Plus_Time_I
    module procedure Time_Plus_Duration, Time_Plus_Duration_I
  end interface

! From http://tycho.usno.navy.mil/mjd.html:

! The Modified Julian Day (MJD) is an abbreviated version of the old Julian Day
! (JD) dating method which has been in use for centuries by astronomers,
! geophysicists, chronologers, and others who needed to have an unambiguous
! dating system based on continuing day counts.
! 
! The JD counts have very little to do with the Julian calendar, which was
! introduced by Julius Caesar (46 BC) and in force until 1582 when Pope Gregory
! XIII directed the use of an improved calendar, now known as the Gregorian
! Calendar. In the case of the Julian day count, the name was given because at
! the time, the Julian calendar was in use and, therefore, the epoch of the day
! count was fixed in respect to it. The JD counts days within one Julian Period
! of exactly 7980 Julian years of 365.25 days.
! 
! Start of the JD count is from 0 at 12 noon 1 JAN -4712 (4713 BC), Julian
! proleptic calendar. Note that this day count conforms with the astronomical
! convention starting the day at noon, in contrast with the civil practice where
! the day starts with midnight (in popular use the belief is widespread that the
! day ends with midnight, but this is not the proper scientific use).
! 
! The Julian Period is given by the time it takes from one coincidence to the
! next of a solar cycle (28), a lunar cycle (19), and the ancient Roman
! Indiction (a tax cycle of 15 years). At any rate, this period is of interest
! only in regard to the adoption of the start, at which time all periods counted
! backwards were in coincidence.
! 
! The Modified Julian Day, on the other hand, was introduced by space scientists
! in the late 1950's. It is defined as
! 
! MJD = JD - 2400000.5, i.e., time since midnight 17 November 1858.
! 
! The half day is subtracted so that the day starts at midnight in conformance
! with civil time reckoning. This MJD has been sanctioned by various
! international commissions such as IAU, CCIR, and others who recommend it as a
! decimal day count which is independent of the civil calendar in use. To give
! dates in this system is convenient in all cases where data are collected over
! long periods of time. Examples are double star and variable star observations,
! the computation of time differences over long periods of time such as in the
! computation of small rate differences of atomic clocks, etc.
! 
! The MJD is a convenient dating system with only 5 digits, sufficient for most
! modern purposes. The days of the week can easily be computed because the same
! weekday is obtained for the same remainder of the MJD after division by 7.
! 
! EXAMPLE: MJD 49987 = WED., 27 SEPT, 1995
! 
! Division of the MJD by 7 gives a remainder of 0. All Wednesdays have
! remainder 0 mod 7 of the MJD.  All Mondays have remainder mod 0 of the JD.

! The MJD (and even more so the JD) has to be well distinguished from this day
! of the year (DOY). This is also often but erroneously called Julian Date,
! when in fact it is a Gregorian Date expressed as number of days in the year.
! This is a grossly misleading practice that was introduced by some who were
! simply ignorant and too careless to learn the proper terminology. It creates
! a confusion which should not be taken lightly. Moreover, a continuation of
! the use of expressions "Julian" or "J" day in the sense of a Gregorian Date
! will make matters even worse. It will inevitably lead to dangerous mistakes,
! increased confusion, and it will eventually destroy whatever standard
! practices exist.
! 
! The MJD has been officially recognized by the International Astronomical
! Union (IAU), and by the Consultative Committee for Radio (CCIR), the advisory
! committee to the International Telecommunications Union (ITU). The pertinent
! document is
! 
! CCIR RECOMMENDATION 457-1, USE OF THE MODIFIED JULIAN DATE BY THE STANDARD
! FREQUENCY AND TIME-SIGNAL SERVICES.
! 
! This document is contained in the CCIR "Green Book," Volume VII. Additional,
! extensive documentation regarding the JD is contained in the Explanatory
! Supplement to the Astronomical Ephemeris and the Nautical Almanac , and in
! the yearbooks themselves, now called The Astronomical Almanac . The Almanac
! for Computers also provides information on such matters.
! 
! NOTE: The MJD is always referred to as a time reckoned in Universal Time (UT)
! or the closely related Coordinated Universal Time (UTC) , International
! Atomic Time (TAI), or Terrestrial Dynamic Time (TDT). The same is not true
! for the DOY. This is usually meant in a local time sense, but in all data
! which are given here at the observatory, we refer the DOY to UT also, except
! where specifically noted. One could call it then something like Universal Day
! of the Year to emphasize the point. However, this would introduce a
! completely new term, not authorized by any convention. Moreover, it is not
! really necessary to use a different term because we simply follow logically
! the same practice of extending a time and date measure to the UT reference as
! we do when we give any date or hour.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  CALDAT_T  -----
  elemental &
  type(time_t) function CALDAT_T ( JULDAY ) result ( TIME )

    ! Returns the calendar date corresponding to a given Julian day.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.

    real(tk), intent(in) :: JULDAY

    call caldat ( julday, time%year, time%month, time%day, &
      &                   time%hours, time%minutes, time%seconds )

  end function CALDAT_T

  ! -----------------------------------------------  CALDAT_T_MJD  -----
  elemental &
  type(time_t) function CALDAT_T_MJD ( MJD ) result ( TIME )

    ! Returns the calendar date corresponding to a given Julian day.
    ! The Modified Julian day is the number of days since Midnight 17 November
    ! 1858.  MJD = JULDAY - 2400000.5

    real(tk), intent(in) :: MJD

    call caldat_mjd ( mjd, time%year, time%month, time%day, &
      &                    time%hours, time%minutes, time%seconds )

  end function CALDAT_T_MJD

  ! ---------------------------------------------  CALDAT_MJD_YMD  -----
  elemental &
  subroutine CALDAT_MJD_YMD ( MJD, YEAR, MONTH, DAY )

    ! Returns the calendar date corresponding to a given Modified Julian day.
    ! The Modified Julian day is the number of days since Midnight 17 November
    ! 1858.  MJD = JULDAY - 2400000.5

    real(tk), intent(in) :: MJD
    integer, intent(out) :: YEAR, MONTH, DAY

    call caldat ( mjd + 2400000.5, year, month, day )

  end subroutine CALDAT_MJD_YMD

  ! -------------------------------------------  CALDAT_MJD_YMDHMS  -----
  elemental &
  subroutine CALDAT_MJD_YMDHMS ( MJD, YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS )

    ! Returns the calendar date corresponding to a given Modified Julian day.
    ! Modified Julian day zero extends from Midnight UTC 17 November 1858
    ! to Midnight UTC 18 November 1858.  MJD = Julian day - 2400000.5

    real(tk), intent(in) :: MJD
    integer, intent(out) :: YEAR, MONTH, DAY, HOURS, MINUTES
    real, intent(out) :: SECONDS
    real :: F

    call caldat ( mjd + 2400000.5, year, month, day )

    ! Get the HMS from the MJD instead of the JD to avoid loss of precision.
    f = mjd - int(mjd)
    hours = 24.0_tk * f
    f = f - hours / 24.0_tk
    minutes = 1440.0_tk * f
    seconds = (f-minutes/1440.0_tk) * 86400.0_tk

  end subroutine CALDAT_MJD_YMDHMS

  ! -------------------------------------------------  CALDAT_YMD  -----
  elemental &
  subroutine CALDAT_YMD ( JULDAY, YEAR, MONTH, DAY )

    ! Returns the calendar date corresponding to a given Julian day.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.

    ! Adapted from IDL procedure of the same name.

    real(tk), intent(in) :: JULDAY
    integer, intent(out) :: YEAR, MONTH, DAY

    integer, parameter :: GREG = 2299161 ! Julian date for 15 Oct 1582 on the
      ! Gregorian calendar.  The Gregorian calendar was adopted on this day,
      ! which is 5 Oct 1582 on the Julian calendar.

    integer :: A, B, C, D, E, J

    j = int(julday + 0.5_tk)
    if ( j >= greg ) then
      a = ((j - 1867216) - 0.25_tk) / 36524.25_tk
      a = j + 1 + a - int(0.25_tk*a)
    else
      a = j
    end if

    b = a + 1524
    c = 6680.0_tk + ((b-2439870)-122.1_tk)/365.25_tk
    d = 365.0_tk * c + (0.25_tk * c)
    e = (b - d) / 30.6001_tk

    day = b - d - int(30.6001_tk * e)
    month = e - 1
    month = mod(month-1,12) + 1
    year = c - 4715
    if ( month > 2 ) year = year - 1
    if ( year <= 0 ) year = year - 1

  end subroutine CALDAT_YMD

  ! -------------------------------------------------  CALDAT_YMDHMS  -----
  elemental &
  subroutine CALDAT_YMDHMS ( JULDAY, YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS )

    ! Returns the calendar date corresponding to a given Julian day.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.

    ! Adapted from IDL function of the same name.

    real(tk), intent(in) :: JULDAY
    integer, intent(out) :: YEAR, MONTH, DAY
    integer, intent(out) :: HOURS, MINUTES
    real, intent(out) :: SECONDS

    integer, parameter :: GREG = 2299161 ! Julian date for 15 Oct 1582 on the
      ! Gregorian calendar.  The Gregorian calendar was adopted on this day,
      ! which is therefore 5 Oct 1582 on the Julian calendar.

    real :: F
    real(tk) :: JH

    call caldat ( julday, year, month, day )

    jh = julday + 0.5_tk
    f = jh - int(jh)
    hours = 24.0_tk * f
    f = f - hours / 24.0_tk
    minutes = 1440.0_tk * f
    seconds = (f-minutes/1440.0_tk) * 86400.0_tk

  end subroutine CALDAT_YMDHMS

  ! ---------------------------------------------------  DURATION  -----
  elemental &
  real(tk) function DURATION ( TIME2, TIME1 )
  ! The number of days between Time1 and Time2.
    type(time_t), intent(in) :: Time2, TIme1
    duration = t_julday(time2) - t_julday(time1)
  end function DURATION

  ! -----------------------------------------  DURATION_PLUS_TIME  -----
  elemental &
  type(time_t) function DURATION_PLUS_TIME ( DURATION, TIME )
  ! The time at TIME + DURATION
    real(tk), intent(in) :: DURATION
    type (time_t), intent(in) :: TIME
    if ( duration - int(duration) == 0.0_tk ) then
      ! Avoid round-off problems for whole-day increments.
      duration_plus_time = &
        & caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
      duration_plus_time%hours = time%hours
      duration_plus_time%minutes = time%minutes
      duration_plus_time%seconds = time%seconds
    else
      duration_plus_time = caldat_t ( t_julday(time) + duration )
    end if
  end function DURATION_PLUS_TIME

  ! ---------------------------------------  DURATION_PLUS_TIME_I  -----
  elemental &
  type(time_t) function DURATION_PLUS_TIME_I ( DURATION, TIME ) result ( R )
  ! The time at TIME + DURATION
    integer, intent(in) :: DURATION
    type (time_t), intent(in) :: TIME
    ! Avoid round-off problems for whole-day increments.
    r = caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
    r%hours = time%hours
    r%minutes = time%minutes
    r%seconds = time%seconds
  end function DURATION_PLUS_TIME_I

  ! ---------------------------------------------------  T_JULDAY  -----
  elemental &
  real(tk) function T_JULDAY ( TIME )

    ! Returns the Julian day for the date given by TIME.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.

    type(time_t), intent(in) ::  Time

    t_julday = julday ( time%year, time%month, time%day, &
      &                 time%hours, time%minutes, time%seconds )

  end function T_JULDAY

  ! ------------------------------------------------------  T_MJD  -----
  elemental &
  real(tk) function T_MJD ( TIME )

    ! Returns the Modified Julian day for the date given by TIME.
    ! Modified Julian day zero extends from Midnight UTC 17 November 1858
    ! to Midnight UTC 18 November 1858.  MJD = Julian day - 2400000.5

    type(time_t), intent(in) ::  Time

    t_mjd = mjd ( time%year, time%month, time%day, &
      &           time%hours, time%minutes, time%seconds )

  end function T_MJD

  ! ----------------------------------------  TIME_MINUS_DURATION  -----
  elemental &
  type(time_t) function TIME_MINUS_DURATION ( TIME, DURATION ) result ( R )
  ! The time at TIME - DURATION
    type (time_t), intent(in) :: TIME
    real(tk), intent(in) :: DURATION
    if ( duration - int(duration) == 0.0_tk ) then
      ! Avoid round-off problems for whole-day increments.
      r = caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
      r%hours = time%hours
      r%minutes = time%minutes
      r%seconds = time%seconds
    else
      r = caldat_t ( t_julday(time) - duration )
    end if
  end function TIME_MINUS_DURATION

  ! --------------------------------------  TIME_MINUS_DURATION_I  -----
  elemental &
  type(time_t) function TIME_MINUS_DURATION_I ( TIME, DURATION ) result ( R )
  ! The time at TIME - DURATION
    type (time_t), intent(in) :: TIME
    integer, intent(in) :: DURATION
    ! Avoid round-off problems for whole-day increments.
    r = caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
    r%hours = time%hours
    r%minutes = time%minutes
    r%seconds = time%seconds
  end function TIME_MINUS_DURATION_I

  ! -----------------------------------------  TIME_PLUS_DURATION  -----
  elemental &
  type(time_t) function TIME_PLUS_DURATION ( TIME, DURATION ) result ( R )
  ! The time at TIME + DURATION
    type (time_t), intent(in) :: TIME
    real(tk), intent(in) :: DURATION
    if ( duration - int(duration) == 0.0_tk ) then
      ! Avoid round-off problems for whole-day increments.
      r = caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
      r%hours = time%hours
      r%minutes = time%minutes
      r%seconds = time%seconds
    else
      r = caldat_t ( t_julday(time) + duration )
    end if
  end function TIME_PLUS_DURATION

  ! ---------------------------------------  TIME_PLUS_DURATION_I  -----
  elemental &
  type(time_t) function TIME_PLUS_DURATION_I ( TIME, DURATION ) result ( R )
  ! The time at TIME + DURATION
    type (time_t), intent(in) :: TIME
    integer, intent(in) :: DURATION
    ! Avoid round-off problems for whole-day increments.
    r = caldat_t ( t_julday(time_t(time%year,time%month,time%day,0,0,0)) + duration )
    r%hours = time%hours
    r%minutes = time%minutes
    r%seconds = time%seconds

  end function TIME_PLUS_DURATION_I

  ! -------------------------------------------------  YMD_JULDAY  -----
  elemental &
  integer function YMD_JULDAY ( YEAR, MONTH, DAY ) result ( J )

    ! Returns the Julian day for the date given by MONTH, DAY, YEAR.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.

    ! The civil calendar doesn't have a year zero.  Instead of issuing an
    ! error message, this function simply assumes years zero and -1 are
    ! the same.

    ! Algorithm from Wikipedia.

    integer, intent(in) ::  MONTH, DAY, YEAR

    ! Gregorian Calendar was adopted on Oct. 15, 1582
    integer, parameter :: GREG = 2299171 ! Incorrect Julian day for 25 Oct 1592
    integer A, F, M, Y

    a = (14-month)/12
    y = year + 4800 - a
    m = month + 12*a - 3
    f = y/100 - y/400

    j = day + (153*m+2)/5 + 365*y + y/4 - f - 32045

    if ( j < greg ) j = j + f - 38

  end function YMD_JULDAY

  ! ----------------------------------------------  YMDHMS_JULDAY  -----
  elemental &
  real(tk) function &
    & YMDHMS_JULDAY ( YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS ) result ( R )

    ! Returns the Julian day for the date given by MONTH, DAY, YEAR.
    ! Julian day zero extends from Noon UTC 1 January 4713 BC to Noon UTC
    ! 2 January 4713 BC.
    
    ! The civil calendar doesn't have a year zero.  Instead of issuing an
    ! error message, this function simply assumes years zero and -1 are
    ! the same.

    ! Algorithm from Wikipedia.

    integer, intent(in) ::  MONTH, DAY, YEAR, HOURS, MINUTES
    real, intent(in) :: SECONDS

    r = julday(year,month,day) - 0.5_tk + &
      & ( ( seconds / 60.0 + minutes ) / 60.0 + hours ) / 24.0

  end function YMDHMS_JULDAY

  ! ----------------------------------------------------  YMD_MJD  -----
  elemental &
  real(tk) function YMD_MJD ( YEAR, MONTH, DAY ) result ( R )

    ! Returns the Modified Julian day for the date given by MONTH, DAY, YEAR.
    ! Modified Julian day zero extends from Midnight UTC 17 November 1858
    ! to Midnight UTC 18 November 1858.  MJD = Julian day - 2400000.5

    integer, intent(in) ::  MONTH, DAY, YEAR

    r = julday(year,month,day,0,0,0.0) - 2400000.5_tk

  end function YMD_MJD

  ! -------------------------------------------------  YMDHMS_MJD  -----
  elemental &
  real(tk) function YMDHMS_MJD ( YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS ) &
    & result ( R )

    ! Returns the Modified Julian day for the date given by MONTH, DAY, YEAR.
    ! Modified Julian day zero extends from Midnight UTC 17 November 1858
    ! to Midnight UTC 18 November 1858.  MJD = Julian day - 2400000.5

    integer, intent(in) ::  MONTH, DAY, YEAR, HOURS, MINUTES
    real, intent(in) :: SECONDS

    ! Add HMS onto YMD MJD instead of getting YMDHMS_JULDAY and subtracting
    ! 2400000.5, to minimize loss of precision from the SECONDS.
    r = mjd ( year, month, day ) + &
      & hours/24.0_tk + minutes/1440.0_tk + seconds/86400.0_tk

  end function YMDHMS_MJD

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Calendar

! $Log$
! Revision 2.1  2004/10/18 20:49:52  vsnyder
! Initial commit
!
