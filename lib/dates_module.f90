! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module dates_module

  ! Converts dates, times between various formats. The resolution of everything
  ! on a scale finer than a day may differ from results obtained using
  ! toolkit-supplied functions because the toolkit compulsively
  ! consults a leap-second file before converting to internal
  ! time representations (TAI93) while this module maintains
  ! its own internal database of leap seconds.
  
  ! There are two other date-handling and time-handling modules:
  !   Calendar.f90
  !   CCSDS_Time.f90
  ! No other module uses them currently.

  use IO_Stuff, only: PrintMessage
  use MLSCommon, only: NameLen, MLSMSG_Warning
  use MLSFinds, only: FindFirst
  use MLSStringlists, only: GetStringElement, NumStringElements
  use MLSStrings, only: Capitalize, CharToInt, Depunctuate, &
    & Indexes, IsAlphabet, IsDigits, Lowercase, &
    & Reverse_Trim, ReadIntsFromChars, SplitWords, WriteIntsToChars

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Cal is calendar date as in "25 Jan 1998" or "25 1 1998"
  ! Separators can be any non-alphanumeric character. Leading 0s no problem
  ! If you want the three items in an order that is not day/month/year then
  ! you can supply the optional second argument perm. This must be a 3-element
  ! integer array with elements 1,2,3 in any order
  ! EUDTF is Extended UDTF -- Y2K compliant variant of UDTF
  ! EUDTF is an integer of form yyyyddd with yyyy==year and ddd==day of year
  ! (I (HCP) invented this, it isn't a standard)
  
  ! utc is a combined dateTime string in one of two forms also
  ! described below under CCSDS
  ! (a) yyy-doyThh:mm:ss.ssssZ
  ! (b) yyy-mm-ddThh:mm:ss.ssssZ

  ! See reformatDate below for extra options to describe
  ! input and output formats

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)

!     (subroutines and functions)
! adddaystoutc       add specified number of days to a utc date-time format
! addhourstoutc      add specified number of hours to a utc date-time format
! addsecondstoutc    add specified number of seconds to a utc date-time format
! buildCalendar      builds a calendar page; 6x7 array of ints or chars
!                      columns are Sunday-Saturday, rows are individual weeks
! cal2eudtf          cal date -> yyyyddd
! ccsds2tai          ccsds -> tai (days, not s)
! ccsds2eudtf        ccsds -> eudtf
! ccsdsa2b           yyyy-mm-dd -> yyyy-DDD
! ccsdsb2a           yyyy-DDD -> yyyy-mm-dd
! dai_to_yyyymmdd    Converts days after Jan 1, 2001 to yyyymmdd
! dateForm           Determines what format a date is in; e.g. 'yyyy-doy'
! datesbetween2utcs  Return the dates between 2 date-times
! daysbetween2utcs   How many days between 2 date-times
! dayOfWeek          Determines what day of the week a date is; e.g. 'Monday'
! daysInMonth        Determines how many days in given month in given year
! daysince2eudtf     days since starting date -> eudtf
! days_in_year       how many days in (leap, normal) year
! eudtf2cal          yyyyddd -> cal date
! eudtf2daysince     eudtf -> days since starting date
! FromUARSDate       Converts a uars date, e.g. 'd0007' to yyyy-mm-dd
! gethid             like tai93s2hid but neither pure nor a function (what?)
! GetStartingDate    Return starting date for the tai format
! hoursbetween2utcs  How many hours between 2 date-times
! hoursinday         How many hours since the start of the day
! isUTCInRange       Is the utc a date within allowable range of dates?
! lastday            How many days in normal month (does anyone use this?)
! nextMoon           utc date, time of next (new) moon 
!                     (or next after input date); less accurate than toolkit AA
! precedesutc        Did utcdate1 precede utcdate2?
! ReformatDate       Turns 'yyyymmdd' -> 'yyyy-mm-dd'; or more general format
! ReformatTime       Turns 'hhmmss.sss' -> 'hh:mm:ss'
! ResetStartingDate  Choose a different starting date for the tai format
! RestoreStartingDate
!                    Revert to Jan 01 1993 for the tai starting date
! secondsbetween2utcs 
!                    How many seconds between 2 date-times
! secondsinday       How many seconds since the start of the day
! splitDateTime      Splits dateTtime into date, time
! tai2ccsds          tai (days, not s) -> ccsds (in "B" format)
! tai93s2hid         tai (s, not days) -> hours-in-day
! tai93s2utc         tai (s, not days) -> yyyy-mm-ddThh:mm:ss.sss (i.e. "B" format)
! timeForm           Determines what format time is in
! toUARSDate         Converts yyyy-mm-dd to a uars date, e.g. 'd0007'
! utc_to_date        Returns date portion from dateTtime; e.g. yyyy-dddThh:mm:ss
! utc_to_time        Returns time portion from dateTtime; e.g. yyyy-dddThh:mm:ss
! utc_to_yyyymmdd    Parses yyyy-mm-ddThh:mm:ss.sss or yyyy-dddThh:mm:ss.sss
! utcForm            Determines what format ('a', 'b', 'n', or 'u') a utc is in
!                     (a) yyyy-mm-ddThh:mm:ss.sss
!                     (b) yyyy-dddThh:mm:ss.sss
!                     (n) yyyydddThh:mm:ss.sss (no-dash)
! utc2tai93s         yyyy-mm-ddThh:mm:ss.sss (i.e. "B" format) -> tai (s, not days)
!                     (see also secondsbetween2utcs)
! yyyyDoy_to_mmdd    Converts yyyy and day-of-year to month and day
! yyyymmdd_to_dai    Converts yyyymmdd to days after Jan 1, 2001
! yyyymmdd_to_Doy    Converts yyyy, mont, and day to day-of-year

! Bugs and limitations:
! These procedures are valid only for dates on or after the tai starting date,
! by default Jan 1 1993. For dates prior to that date, you must first call
! ResetStartingDate.
! An internal database is created of leap seconds (leapSecDates below). 
! You must extend it for use after the last date shown if any subsequent
! leap seconds are declared.

! === (end of toc) ===

! === (start of api) ===
! char* adddaystoutc (char* utc, int days)
! char* addhourstoutc (char* utc, int hours)
! char* addsecondstoutc (char* utc, dble seconds)
! buildCalendar ( int year, int month, int days(6, 7), [int daysOfYear(6,7)] )
! dai_to_yyyymmdd (int dai, int yyyy, int mm, int dd, [char* startingDate])
! dai_to_yyyymmdd (int dai, char* str, [char* startingDate])
! char* dateForm( char* date )
! char* dayOfWeek (char* date, [char* fromForm])
! datesbetween2utcs ( char* utc1, char* utc2, char* dates(:) )
! int daysbetween2utcs (char* utc1, char* utc2)
! int daysInMonth (int month, int year)
! gethid( dble tai, [log leapsec] )
! GetStartingDate ( char* Date )
! int hoursbetween2utcs (char* utc1, char* utc2)
! dble hoursinday (char* utc1)
! log IsUTCInRange ( char* utc,  [char* utc1], [char* utc2] )
! char* nextMoon ([char* date], [phase])
! log precedesutc ( char* date1,  char* date2 )
! char* ReformatDate (char* date, [char* fromForm], [char* toForm])
! char* ReformatTime (char* time, [char* form])
! resetStartingDate ( char* newDate )
! restoreStartingDate
! dble secondsbetween2utcs ( char* utc1, char* utc2, [log leapsec] )
! dble secondsinday (char* utc1)
! splitDateTime(char* utc, int ErrTyp, char* date, char* time, [log strict])
! char* tai2ccsds( int tai )
! dble tai93s2hid( dble tai, [log leapsec] )
! char* tai93s2utc( dble tai, [log leapsec] )
! char* timeForm( char* time )
! char utcForm (char* utc)
! utc_to_date(char* utc, int ErrTyp, char* date, [log strict], [char* utcAt0z])
! utc_to_time(char* utc, int ErrTyp, char* time, [log strict])
! dble utc2tai93s( char*  utc, [log leapsec] )
! yyyyDoy_to_mmdd (int yyyy, int mm, int dd, int Doy)
! yyyymmdd_to_dai (int yyyy, int mm, int dd, int dai, [char* startingDate])
! yyyymmdd_to_dai (char* str, int dai, [char* startingDate])
! yyyymmdd_to_Doy (int yyyy, int mm, int dd, int Doy)

! Note: in fromForm and in toForm the following rules are in effect
! non-alphabetic symbols like '-', '/', and ' ' take their own values
! otherwise
!                  date formats
!    key       is replaced by                        example
!    ---       --------------                       ---------
!    dd        2-digits day of month                07     
!   doy        3-digits day of year with d-prefix   d270     
!   Doy        3-digits day of year w/o d-prefix    270     
!    mm        2-digits month number of year        01     
!   yyyy       4-digits year number                 2005     
!     M        full month name                      September

! Putting these together we would have for the same date
!   "M dd yyyy" => "September 07 2005"
!   "yyyy-doy " => "2005-d250"

!                  time formats
!    key       is replaced by                        example
!    ---       --------------                       ---------
!    hh        2-digits hours in day (of 24)           17     
!    HH        2-digits hours in day (of 12)           05 (adds 'am' or 'pm'     
!    mm        2-digits minutes in hour (of 60)        30     
!    ss        seconds in minute (of 60)               30.000000     
!    SS        2-digits seconds in minute (of 60)      30

! Putting these together we would have for the same time
!   "hh:mm:ss" => "17:30:30.0000000"
!   "HH:mm   " => "5:30 pm"

! Special values of fromForm (either ' ' or '*') turn on an
! "automatic recognition" feature which automatically recognizes
! any of the most common date formats

! time format codings are governed by the following rules
!    key       is replaced by             
!    ---       --------------             
!    hh        2-digits in range 00 .. 23 (hours)
!    mm        2-digits in range 00 .. 59 (minutes)
!    ss        2-digits in range 00 .. 59 (seconds)
!    HH        2-digits in range 00 .. 12 (hours) (followed by "AM" or "PM"

  ! Since starting this, I (HCP again) discover that the standard text 
  ! format for EOS will be CCSDS format. This comes in two sorts: 
  ! Form A:  yyyy-mm-dd   where mm= month, dd=day in month
  ! Form B:  yyyy-DDD     where DDD=day of year. No spaces. 
  ! These are clearly special cases of calendar format, and eudtf format
  ! but we provide convenience calls to change from form A to form B

  ! The SDP Toolkit has some sort of date format that fits in a double
  ! precision number. It is called the TAI93 format and is the number 
  ! of seconds since midnight, 1 Jan 1993. We provide a format called
  ! daysince , which is the integer number of days since a particular 
  ! midnight. This can be used to get UARS day or TAI93 time to an accuracy
  ! of 1 Day. 

  ! Be warned that the format we call TAI here is in DAYS 
  ! (unless specified otherwise)
  ! while genuine TAI93 is in SECONDS. 
  
  ! The toolkit's TAI93 accounts for leap seconds
  
  ! PAW chose a different starting date for dai (might be better called dai01)
  ! which is January 1 2001
  ! so try to keep distinct our 3 numerical date types
  ! name             meaning                             numerical type
  ! ----             -------                             --------------
  ! tai93s (toolkit) seconds since midnight 1993 Jan 1        d.p.
  ! tai93  (here) days since midnight 1993 Jan 1              int
  ! dai01  (here) days since midnight 2001 Jan 1              int
  !
  ! Also we allow several ways of encoding dates as strings
  ! which is perhaps the reason this module has grown so large
  
  ! There are two internal date formats:
  ! eudtf          -  yyyddd expressed as an integer, e.g. 1993001 is '1993-001'
  ! MLSDate_Time_T - a user-defined type containing 3 fields: 
  !   dai
  !   seconds
  !   leapseconds
  
  ! Further notes and Limitations:
  ! It would be useful for this module to supply functions converting among
  ! our 3 numerical date types
  ! E.g.,
  !    tai93 (in days) = dai01 + nDaysOffset
  ! where we can 
  !   nDaysOffset = daysbetween2utcs( '1993-01-01', '2001-01-01' )
  ! and, ignoring the effect of leap seconds, 
  !   tai93s (toolkit) = SECONDSINADAY*tai93 (in days)
  !
  ! It would not be that difficult to enable us to take advantage
  ! of a leapsec file, and, after having read it, fill out our leapsec database.
  ! Difficult or easy, we have not yet bestirred ourselves. Instead we
  ! define a character-valued array of leapsecond dates
  
  ! Our implementation of leapseconds is not strictly accurate:
  ! (a) Our database includes dates prior to 1972, contrary to the usual definition
  ! (b) Our implementation assumes that each leapsecond adds exactly 1 second
  !     the definition allows for possible negative leap seconds;
  !     early leap seconds were sometimes even fractional
! === (end of api) ===

! Further notes:
! (1) We maintain a private datatype, MLSDate_Time_T
!     it is an intermediary in various operations and conversions
!     Would it be useful to uncloak it so callers might access it directly?
! (2) We call eudtf an internal data format, and yet procedures
!     to convert to and from it are public. Contradictory?
! (3) This module has grown rather long. Should we consider
!     Splitting it into two? If so, what would be the criteria for
!     grouping procedures into a common module?

  ! Here are the functions this module provides
  public :: Adddaystoutc, Addhourstoutc, Addsecondstoutc, Buildcalendar, &
   & Cal2eudtf, Ccsds2tai, Ccsds2eudtf, Ccsdsa2b, Ccsdsb2a, &
   & Dai_To_Yyyymmdd, Dateform, Datesbetween2utcs, Dayofweek, &
   & Daysbetween2utcs, Daysinmonth, Daysince2eudtf, Days_In_Year, &
   & Eudtf2cal, Eudtf2daysince, Fromuarsdate, Gethid, Getstartingdate, &
   & Hoursbetween2utcs, Hoursinday, Isutcinrange, &
   & Lastday, Monthname, Nextmoon, Precedesutc, Reformatdate, Reformattime, &
   & Resetstartingdate, Restorestartingdate, &
   & Secondsbetween2utcs, Secondsinday, Splitdatetime, &
   & Tai2ccsds, Tai93s2hid, Tai93s2utc, Timeform, Touarsdate, &
   & Utcform, Utc_To_Date, Utc_To_Time, Utc_To_Yyyymmdd, Utc2tai93s, &
   & Yyyydoy_To_Mmdd, Yyyymmdd_To_Dai, Yyyymmdd_To_Doy

 interface Buildcalendar
    module procedure  Buildcalendar_Ints, Buildcalendar_Str
  end interface

  interface Dai_To_Yyyymmdd
    module procedure Dai_To_Yyyymmdd_Str, Dai_To_Yyyymmdd_Ints
  end interface

  interface Dayofweek
    module procedure Dayofweek_Int, Dayofweek_Str
  end interface

  interface Dump
    module procedure Dumpdatetime
  end interface

 interface Secondsinday
    module procedure Secondsindaydble
  end interface

  interface Switch
    module procedure Switch_ints
  end interface

  interface Utc_To_Yyyymmdd
    module procedure Utc_To_Yyyymmdd_Strs, Utc_To_Yyyymmdd_Ints
  end interface

  interface Yyyydoy_To_Mmdd
    module procedure Yyyydoy_To_Mmdd_Ints
  end interface

  interface Yyyymmdd_To_Dai
    module procedure Yyyymmdd_To_Dai_Str, Yyyymmdd_To_Dai_Ints
  end interface

  interface Yyyymmdd_To_Doy
    module procedure Yyyymmdd_To_Doy_Ints, Yyyymmdd_To_Doy_Str
  end interface

  ! utc_to_yyyymmdd
  integer, public, parameter :: Invalidutcstring = 1
  integer, public, parameter :: Maxutcstrlength = 32

  integer, private, parameter :: Secondsinaday = 24*60*60
  
  ! Should we read the following from a file, like the Toolkit does?
  ! Could we even automatically generate that include file from the
  ! toolkit's leapsec.dat?
  ! Let's give it a try.
!  character(len=*), dimension(41), parameter :: leapSecDates = (/ &
!     '1961 JAN 1', & 
!     '1961 AUG 1', & 
!     '1962 JAN 1', & 
!     '1963 NOV 1', & 
!     '1964 JAN 1', & 
!     '1964 APR 1', & 
!     '1964 SEP 1', & 
!     '1965 JAN 1', & 
!     '1965 MAR 1', & 
!     '1965 JUL 1', & 
!     '1965 SEP 1', & 
!     '1966 JAN 1', & 
!     '1968 FEB 1', & 
!     '1972 JAN 1', & 
!     '1972 JUL 1', & 
!     '1973 JAN 1', & 
!     '1974 JAN 1', & 
!     '1975 JAN 1', & 
!     '1976 JAN 1', & 
!     '1977 JAN 1', & 
!     '1978 JAN 1', & 
!     '1979 JAN 1', & 
!     '1980 JAN 1', & 
!     '1981 JUL 1', & 
!     '1982 JUL 1', & 
!     '1983 JUL 1', & 
!     '1985 JUL 1', & 
!     '1988 JAN 1', & 
!     '1990 JAN 1', & 
!     '1991 JAN 1', & 
!     '1992 JUL 1', & 
!     '1993 JUL 1', & 
!     '1994 JUL 1', & 
!     '1996 JAN 1', & 
!     '1997 JUL 1', & 
!     '1999 JAN 1', & 
!     '2006 JAN 1', & 
!     '2009 JAN 1', & 
!     '2012 JUL 1', & 
!     '2015 JUL 1', & 
!     '2017 JAN 1'  & 
!    /) 
  include 'LeapSecDates.f9h'

  ! These somewhat similar parameters are used in the separately-coded
  ! but redundant functions and procedures moved here from their
  ! original slots in MLSStrings
  integer, parameter :: Yearmax = 4999  ! Conversions valid until 4999 AD

  ! The following arrays contains the maximum permissible day for each month
  ! indexed by month number
  ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
  !        Months in leap-years
  integer, dimension(-1:12), parameter :: DayMaxLY = (/ &
    & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)
  !        Months in normal-years
  integer, dimension(-1:12), parameter :: DayMaxNY = (/ &
    & 365, 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
    & /)

  include 'DayMonthNames.f9h'
  ! This is the utc for the first (new) moon of 2001
  ! We assume that the moon's phase repeats perfectly; it would be 
  ! more accurate to use the toolkit's AA (Astronimcal Almanac) components
  ! As it is we may be off by as much as one day.
  ! Initializing
  character(len=*), parameter :: Firstnewmoon = &
    & '2001-024T13:07:00.0000Z'
  double precision, parameter :: Lunarperiod = 60.d0*(44 + 60.d0*( &
    & 12 + 24.d0*29 ) ) ! 29d 12h 44m

  ! These are the starting dates for the two date formats we drag around:
  ! the first is for the tai format
  character(len=16)            :: TAIStartingDate   = '1993-01-01'
  ! the second is for our internal MLSdate_Time format
  character(len=16)            :: MLSStartingDate   = '2001-01-01'
  
  ! This is a private type used only internally
  ! Note we don't bother with leap seconds
  ! which rather limits its accuracy and usefulness
  ! One easy improvement would be to incorporate the starting
  ! date--that would allow us to unify the different numerical date types
  type MLSDate_Time_T
    integer :: dai = 0                  ! days after 1 Jan 2001
    double precision :: seconds = 0.00  ! seconds after midnight
    integer :: LeapSeconds = 0          ! optional leap seconds to be added
    ! character(len=16) :: startingDate = MLSStartingDate
  end type MLSDate_Time_T
  
  ! This will be the starting date for all MLSDate_Time_T
  ! We could switch to letting each MLSDate_Time_T keep its
  ! own StartingDate.
  character(len=16), save :: DTStartingDates = '2001-01-01'
  ! This accounts for the number of days between the dates
  ! '1993-01-01' and '2001-01-01':
  ! 6 non-leap years and 2 leap years
  integer, parameter :: DAI93TODAI01    = 6*365 + 2*366
  integer            :: DAISTARTTODAI01 = DAI93TODAI01
  
contains
  ! ----------- Public funtions and subroutines ------------
  ! Many of these merely convert from one format to another
  ! Often they work by first converting to MLSDate_Time_T
  ! and then from MLSDate_Time_T to whatever end format you wanted
  ! ---------------------------------------------  adddaystoutc  -----
  function adddaystoutc( utc, days ) result(after)
    ! Given a utc return a date later (or earlier) by days
    ! Args
    character(len=*), intent(in) :: utc
    integer, intent(in)          :: days
    character(len=len(utc))      :: after
    ! Internal
    type(MLSDate_time_T)         :: datetime
    ! Executable
    datetime = utc2datetime(utc)
    datetime%dai = datetime%dai + days
    after = datetime2utc(datetime)
  end function adddaystoutc

  ! ---------------------------------------------  addhourstoutc  -----
  function addhourstoutc( utc, hours ) result(after)
    ! Given a utc return a date later (or earlier) by hours
    ! Args
    character(len=*), intent(in) :: utc
    integer, intent(in)          :: hours
    character(len=len(utc))      :: after
    ! Internal
    type(MLSDate_Time_T)         :: datetime
    ! Executable
    datetime = utc2datetime(utc)
    datetime%seconds = datetime%seconds + 60*60*hours
    ! call dumpDateTime(datetime)
    call reducedatetime(datetime)
    ! call dumpDateTime(datetime)
    after = datetime2utc(datetime)
  end function addhourstoutc

  ! ---------------------------------------------  addsecondstoutc  -----
  function addsecondstoutc( utc, seconds ) result(after)
    ! Given a utc return a date later (or earlier) by seconds
    ! Args
    character(len=*), intent(in) :: utc
    double precision, intent(in) :: seconds
    character(len=len(utc))      :: after
    ! Internal
    type(MLSDate_Time_T)         :: datetime
    ! Executable
    datetime = utc2datetime(utc)
    ! call dump ( datetime )
    datetime%seconds = datetime%seconds + seconds
    ! call dump ( datetime )
    call reducedatetime(datetime)
    ! call dump ( datetime )
    after = datetime2utc(datetime)
  end function addsecondstoutc

  function ccsds2eudtf(ccsds) result (eudtf)
    !Converts CCSDS dates to eudtf
    !----args----!
    character(len=*),intent(in)::ccsds
    !---function--result---!
    integer::eudtf
    !---local vars-----!
    character(len=30)::ccsdsb,ccsdsi
    integer::year,day
    !----Executable----!
    ccsdsi=adjustl(ccsds)
    ! Check if it looks like CCSDSA format and, if so,  bludgeon into B format
    if (ccsdsi(8:8)=="-") then
       ccsdsb=ccsdsa2b(ccsdsi)
    else
       ccsdsb=ccsdsi
    endif
    ! Get year and day out of B format
    read(unit=ccsdsb(1:4),fmt="(i4)")year
    read(unit=ccsdsb(6:8),fmt="(i3)")day
    eudtf=year*1000+day
  end function ccsds2eudtf

  ! Converts CCSDS date to TAI Date. This is DAYS since the TAI93=0
  ! time _not_ genuine TAI93 which is in seconds. 
  function ccsds2tai(ccsds) result (tai)
    character(len=*),intent(in)::ccsds
    !---function--result---!
    integer::tai
    !----local -----!
    integer:: eudtf
    eudtf=ccsds2eudtf(ccsds)
    tai=eudtf2daysince(eudtf,1993001)
  end function ccsds2tai

  function eudtf2ccsds(eudtf) result (ccsds)
    ! Converts eudtf to CCSDS dates (in "B" format)
    ! Who concocted the eudtf?
    !----args----!
    integer, intent(in) :: eudtf
    !---function--result---!
    character(len=8) :: ccsds
    !---local vars-----!
    integer::year,day
    !----Executable----!
    year = eudtf / 1000
    day = eudtf - 1000*year
    write(unit=ccsds(1:4),fmt="(i4)")year
    write(unit=ccsds(6:8),fmt="(i3.3)")day
    ccsds(5:5) = '-'
  end function eudtf2ccsds

  logical function isUTCInRange ( UTC, lower, upper ) result( itIs )
    ! Is the UTC within an allowable range?
    ! Args
    character(len=*), intent(in)             :: UTC
    character(len=*), intent(in), optional   :: lower
    character(len=*), intent(in), optional   :: upper
    ! Internal variables
    type(MLSDate_Time_T)                     :: dateTime
    type(MLSDate_Time_T)                     :: lowerDT
    type(MLSDate_Time_T)                     :: upperDT
    ! Executable
    datetime = utc2datetime( utc )
    if ( present(lower) .and. present(upper) ) then
      lowerDT = utc2datetime( lower )
      upperDT = utc2datetime( upper )
      itIs = isDateTimeInRange( dateTime, lowerDT, upperDT )
    else
      itIs = isDateTimeInRange( dateTime )
    endif
  end function isUTCInRange

  ! Converts TAI date to CCSDS Date. This is DAYS since the TAI93=0
  ! time _not_ genuine TAI93 which is in seconds. 
  function tai2ccsds(tai) result (ccsds)
    integer, intent(in) :: tai
    !---function--result---!
    character(len=8)::ccsds
    !----local -----!
    integer:: eudtf
    eudtf = daysince2eudtf(tai,1993001)
    ccsds = eudtf2ccsds(eudtf)
  end function tai2ccsds

  ! Converts TAI in seconds to hours-in-day.
  ! Unless, option leapsec is present and TRUE, ignore leap seconds 
  elemental function tai93s2hid( tai93s, leapsec ) result (hid)
    double precision, intent(in)   :: tai93s
    logical, optional, intent(in) :: leapsec
    !---function--result---!
    double precision :: hid
    !----local -----!
    type(MLSDate_Time_T)           :: datetime
    ! integer :: THEDAY
    ! Executable
    dateTime = tai93s2datetime( tai93s, leapsec )
    hid = dateTime%seconds/3600
    if ( .not. present(leapsec) ) return
    if ( .not. leapsec ) return
    hid = (dateTime%seconds - datetime%leapseconds)/3600
    ! theDay = (tai93s+10.d0) / (24*3600.d0)
    ! hid = (tai93s - theDay*24*3600.d0) / 3600
  end function tai93s2hid

  ! Converts TAI in seconds to hours-in-day.
  ! Unless, option leapsec is present and TRUE, ignore leap seconds 
  ! Here so we may dump if the mood strikes
  subroutine gethid( tai93s, hid, leapsec )
    double precision, intent(in)   :: tai93s
    double precision, intent(out)  :: hid
    logical, optional, intent(in)  :: leapsec
    !----local -----!
    type(MLSDate_Time_T)           :: datetime
    ! integer :: THEDAY
    ! Executable
    dateTime = tai93s2datetime( tai93s, leapsec )
    ! call dump ( dateTime )
    hid = dateTime%seconds/3600
    if ( .not. present(leapsec) ) return
    if ( .not. leapsec ) return
    hid = (dateTime%seconds - datetime%leapseconds)/3600
    ! theDay = (tai93s+10.d0) / (24*3600.d0)
    ! hid = (tai93s - theDay*24*3600.d0) / 3600
  end subroutine gethid

  ! Converts TAI in seconds to utc Date.
  ! Unless, option leapsec is present and TRUE, ignore leap seconds 
  function tai93s2utc( tai93s, leapsec ) result (utc)
    double precision, intent(in)   :: tai93s
    logical, optional, intent(in) :: leapsec
    !---function--result---!
    character(len=MAXUTCSTRLENGTH) :: utc
    !----local -----!
    type(MLSDate_Time_T)           :: datetime
    ! Executable
    dateTime = tai93s2datetime( tai93s, leapsec )
    ! print *, 'StartingDate (TAI,MLS) ', TAIStartingDate, MLSStartingDate
    ! print *, 'dateTime&dai,seconds, leap ', &
    !   & dateTime%dai, dateTime%seconds, dateTime%leapseconds
    utc = datetime2utc( datetime, leapsec )
  end function tai93s2utc

  ! ---------------------------------------------  datetime2utc  -----
  function datetime2utc( datetime, leapsec ) result(utc)
    ! Given an MLSDate_Time_T return a utc
    ! Args
    type(MLSDate_Time_T), intent(in)  :: datetime
    logical, optional, intent(in) :: leapsec
    character(len=MAXUTCSTRLENGTH)    :: utc
    ! Internal
    type(MLSDate_Time_T)         :: datetimeRdcd ! So we don't clobber datetime
    integer                      :: dai93
    character(len=16)            :: hhmmss
    character(len=16)            :: yyyymmdd
    ! Executable
    datetimeRdcd = datetime
    call reduceDatetime( datetimeRdcd, leapsec )
    ! print *, 'After reduction, dateTime&seconds, leap ', dateTime%seconds, dateTime%leapseconds
    ! call dump ( datetimeRdcd )
    ! Now convert to our internal representations
    if ( datetimeRdcd%dai < 0 ) then
      dai93 = datetimeRdcd%dai + daysbetween2utcs( TAIStartingDate, MLSStartingDate )
      call dai_to_yyyymmdd( dai93, yyyymmdd, startingDate=TAIStartingDate )
      ! print *, 'dai93 ', dai93
      ! print *, 'yyyymmdd ', yyyymmdd
    else
      call dai_to_yyyymmdd( datetimeRdcd%dai, yyyymmdd )
    endif
    hhmmss = '00:00:00'
    utc = yyyymmdd(1:4) // '-' // yyyymmdd(5:6) // '-' // yyyymmdd(7:8)
    if ( datetimeRdcd%seconds /= 0. ) then
      hhmmss = s2hhmmss(datetimeRdcd%seconds)
      ! print *, 'hhmmss: ', hhmmss
      ! print *, 'utc: ', utc
    endif
    utc = trim(utc) // 'T' // adjustl(hhmmss)
  end function datetime2utc

  ! Get the starting date
  subroutine GetStartingDate ( TAIStartDate, MLSStartDate, DateOffset )
    character(len=*), optional, intent(out) :: TAIStartDate
    character(len=*), optional, intent(out) :: MLSStartDate
    integer, optional, intent(out)          :: DateOffset
    if ( present(TAIStartDate) ) TAIStartDate = TAIStartingDate
    if ( present(MLSStartDate) ) MLSStartDate = MLSStartingDate
    if ( present(dateOffset) )   DateOffset   = DAISTARTTODAI01
  end subroutine GetStartingDate

  ! Reset the starting date
  subroutine resetStartingDate ( newTAIDate, newMLSDate )
    character(len=*), optional, intent(in) :: newTAIDate
    character(len=*), optional, intent(in) :: newMLSDate
    ! print *, 'Resetting starting date to ' // newdate
    if ( present(newTAIDate) ) TAIStartingDate = newTAIDate
    if ( present(newMLSDate) ) MLSStartingDate = newMLSDate
    DAISTARTTODAI01 = daysbetween2utcs( TAIStartingDate, MLSStartingDate )
  end subroutine resetStartingDate

  ! Restore the starting date
  subroutine RestoreStartingDate
    TAIStartingDate = '19930101'
    MLSStartingDate = '20010101'
    DAISTARTTODAI01 = DAI93TODAI01
  end subroutine RestoreStartingDate

  ! Converts TAI93 in seconds to an MLS DateTime datatype
  elemental function tai93s2datetime( tai93s, leapsec ) &
    & result ( datetime )
    double precision, intent(in)           :: tai93s
    logical, optional, intent(in)          :: leapsec
    !---function--result---!
    type(MLSDate_Time_T)           :: datetime
    !----local -----!
    integer          :: dai93         ! Number of days after 1 Jan 1993
    logical :: myLeapSec
    ! Executable
    myLeapSec = .false.
    if ( present(leapsec) ) myLeapSec = leapsec
    dai93 = int(tai93s / 86400 + 0.5d-8)
    if ( myLeapSec ) then
      dateTime%leapSeconds = HowManyleapSeconds( dai93, start='19930101' )
      datetime%seconds = tai93s - 86400*dai93 !  - dateTime%leapSeconds
    else
      datetime%seconds = tai93s - 86400*dai93
    endif
    ! Now apply offset converting dai93 to dai (2001)
    datetime%dai = dai93 - DAISTARTTODAI01 ! DAI93TODAI01
    call reducedatetime( datetime, leapsec )
  end function tai93s2datetime

  ! --------------------------------------------- buildCalendar ---
  ! Build a calendar page (a 6x7 array)
  ! each entry is a day-of-month, or else 0
  ! E.g., for Sept 2007
  !( S  M  T  W  T  F  S )
  !  0  0  0  0  0  0  1
  !  2  3  4  5  6  7  8
  !    .   .   .
  ! 23 24 25 26 27 28 29
  ! 30  0  0  0  0  0  0

  ! So each column represents a day-of-week, 
  ! and each row represents a week
  subroutine buildCalendar_ints ( year, month, days, daysOfYear )
    ! Args
    integer, intent(in)                  :: year
    integer, intent(in)                  :: month
    integer, dimension(:,:), intent(out) :: days
    integer, dimension(:,:), optional, intent(out) :: daysOfYear
    ! Internal variables
    integer :: daiFirst, daiLast
    integer :: day
    integer :: dayOfYear
    integer :: wkdy
    integer :: last
    integer :: row
    ! Executable
    days = 0
    if ( present(daysOfYear) ) daysOfYear = 0
    if ( size(days,1) < 6 ) return ! Too few rows
    if ( size(days,2) < 7 ) return ! Too few columns
    if ( month < 1 .or. month > 12 ) return ! Illegal month
    ! What is last day of the month?
    if ( leapyear(year) ) then
      last = DAYMAXLY(month)
    else
      last = DAYMAXNY(month)
    endif
    call yyyymmdd_to_dai_ints(year, month, 1, daiFirst )
    call yyyymmdd_to_dai_ints(year, month, last, daiLast )
    call yyyymmdd_to_doy_ints( year, month, 1, dayOfYear )
    row = 1
    ! print *, 'last day of month ', last
    ! print *, 'week day of 1st   ', DayNumberOfWeek( daiFirst )
    ! print *, 'week day of last  ', DayNumberOfWeek( daiLast )
    do day=1, last
      wkdy = DayNumberOfWeek( daiFirst + day - 1 )
      days(row, wkdy) = day
      if ( present(daysOfYear) ) DaysOfYear(row, wkdy) = dayOfYear
      dayOfYear = dayOfYear + 1
      if ( wkdy > 6 ) row = row + 1 ! Tomorrow will begin a new week
    enddo    
  end subroutine buildCalendar_ints

  subroutine buildCalendar_str ( year, month, strdays )
    ! Args
    integer, intent(in)                    :: year
    integer, intent(in)                    :: month
    character, dimension(:,:), intent(out) :: strdays
    ! Internal variables
    integer, dimension(size(strdays,1),size(strdays,1)) :: days
    ! Executable
    strdays = ' '
    call buildCalendar( year, month, days)
    call writeIntsToChars( days, strdays, &
      & specialInts = (/ 0 /), specialChars = (/ ' ' /) )
  end subroutine buildCalendar_str

  ! ---------------------------------------------  datesbetween2utcs  -----
  subroutine datesbetween2utcs( utc1, utc2, dates )
    ! Given two utcs return an array of dates between them
    ! Args
    character(len=*), intent(in)   :: utc1, utc2
    character(len=*), dimension(:) :: dates
    ! Internal
    integer                      :: dai, inDates
    type(MLSDate_Time_T)         :: datetime1, datetime2
    ! Executable
    datetime1 = utc2datetime(utc1)
    datetime2 = utc2datetime(utc2)
    ! call dump(datetime1)
    ! call dump(datetime2)
    dates = ' '
    do dai = datetime1%dai, datetime2%dai
      inDates = dai - datetime1%dai + 1
      if ( inDates > size(dates) ) return
      call dai_to_yyyymmdd ( dai, dates( inDates ) )
    enddo
  end subroutine datesbetween2utcs

  ! ---------------------------------------------  dayOfWeek  -----
  function dayOfWeek_int(dai) result(day)
    ! Given a date as a dai return day of week
    ! Args
    integer, intent(in)  :: dai
    character(len=9)     :: day
    ! Internal
    integer :: daynum
    ! Executable
    day = ' '
    daynum = DayNumberOfWeek(dai)
    if ( daynum > 0 .and. daynum < 8 ) day = DAYSOFWEEK(daynum)
  end function dayOfWeek_int

  function dayOfWeek_str(date, fromForm) result(day)
    ! Given a date as a date string return day of week
    ! Args
    character(len=*), intent(in)            :: date
    character(len=*), optional, intent(in)  :: fromform ! input format
    character(len=9)     :: day
    ! Internal
    integer :: dai
    integer :: daynum
    character(len=16) :: myFromForm
    character(len=16) :: yyyymmdd
    ! Executable
    day = ' '
    myFromForm = ' '
    if ( present(fromForm) ) myFromForm = fromForm
    if ( len_trim(myFromForm) < 1 ) myFromForm = dateForm(date)
    yyyymmdd = reFormatDate(date, &
        & fromForm=trim(myFromForm), toForm='yyyymmdd')
    call yyyymmdd_to_dai(yyyymmdd, dai)
    daynum = DayNumberOfWeek(dai)
    if ( daynum > 0 .and. daynum < 8 ) day = DAYSOFWEEK(daynum)
  end function dayOfWeek_str

  ! ---------------------------------------------  daysbetween2utcs  -----
  function daysbetween2utcs( utc1, utc2 ) result(days)
    ! Given two utcs return the number days between them
    ! Args
    character(len=*), intent(in) :: utc1, utc2
    integer                      :: days
    ! Internal
    type(MLSDate_Time_T)         :: datetime1, datetime2
    ! Executable
    datetime1 = utc2datetime(utc1)
    datetime2 = utc2datetime(utc2)
    days = datetime2%dai - datetime1%dai
  end function daysbetween2utcs

  ! ---------------------------------------------  daysInMonth  -----
  elemental function daysInMonth( month, year ) result(days)
    ! Given a month, year, return the number of days in the month
    ! Return 0 if the month is not a valid number
    ! Args
    integer, intent(in)          :: month
    integer, intent(in)          :: year
    integer                      :: days
    ! Executable
    ! Valid month? If not, return 0
    if ( month < lbound(DayMaxLY, 1) .or. month > ubound(DayMaxLY, 1) ) then
      days = 0
    elseif ( leapyear(year) ) then
       days = DAYMAXLY(month)
    else
       days = DAYMAXNY(month)
    endif
  end function daysInMonth

  ! -----------------------------------------------  DumpDateTime  -----
  ! Not yet a public procedure
  subroutine dumpDateTime( dateTime, name )
    ! Dump an MLSDate_Time_T
    ! Args
    type(MLSDate_Time_T), intent(in)         :: datetime
    character(len=*), optional, intent(in)   :: name
    ! Internal variables
    character(len=32)                         :: daiStr
    character(len=32)                         :: secondsStr
    ! Executable
    print *, "MLSDate_Time: (initialized at " // trim(MLSStartingDate) // ') '
    if ( present(name) ) then
      print *, trim(name)
    endif
    write( daiStr, * ) datetime%dai
    write( secondsStr, * ) datetime%seconds
    print *,  'dai: ' // trim(daiStr)
    print *, 'seconds: ' // trim(secondsStr)

  end subroutine dumpDateTime

  ! ---------------------------------------------  hoursbetween2utcs  -----
  function hoursbetween2utcs( utc1, utc2 ) result(hours)
    ! Given a utc return a date later (or earlier) by days
    ! Args
    character(len=*), intent(in) :: utc1, utc2
    integer                      :: hours
    ! Executable
    hours = secondsbetween2utcs( utc1, utc2 ) / (60*60)
  end function hoursbetween2utcs

  ! ---------------------------------------------  hoursinday  -----
  function hoursinday( utc1 ) result(hours)
    ! Given a utc return a date later (or earlier) by days
    ! Args
    character(len=*), intent(in) :: utc1
    double precision             :: hours
    ! Internal
    type(MLSDate_Time_T)         :: datetime1
    ! Executable
    datetime1 = utc2datetime(utc1)
    hours = datetime1%seconds/3600
  end function hoursinday

  function lastday(imonth) result (day)
    ! See also daysInMonth
    integer,intent(in)::imonth
    integer::day
    if(imonth < 1 .or. imonth >12) then
       call PrintMessage( MLSMSG_Warning, ModuleName,&
       "in function lastday: month is out of range")
       day=31
    else
       day=DAYMAXNY(imonth)
    endif
  end function lastday

  elemental function days_in_year(year) result(days) 
    !Public function returns no. of days in a year
    integer,intent(in)::year
    integer::days
    if ( leapyear(year) ) then
       days=366
    else
       days=365
    endif

  end function days_in_year

  function eudtf2daysince(eudtf,eudtf0) result(daysince)
    ! converts eudtf to days since a fixed day ( also given as EUDTF format)
    ! set eudtf0 to 11 Sept 1991 to get UARS days
    ! Set eudtf0 to 1 Jan 1993 for TAI day 
    ! Note that this gives you DAYS since TAI93=0, not TAI93 itself, (which
    ! is _seconds_ since TAI93=0)
    !-----args------------!
    integer,intent(in) :: eudtf,eudtf0
    !----result-----------!
    integer::daysince
    !------local vars ----!
    integer::year,year0,day,day0,daysinyear,direction
    !----------Executable----!

    year0=eudtf0/1000
    year=eudtf/1000
    day0=modulo(eudtf0,1000)
    day=modulo(eudtf,1000)
    daysince=0

    !    print*,"Start",year0,year,day0,day,daysince
    if (year < year0) then 
       
       direction=year
       year=year0
       year0=direction
       
       direction=day
       day=day0
       day0=direction

       direction=-1
    else
       direction=1
    endif

    yrsloop:do
        daysinyear=days_in_year(year0)
        if (year0 == year) then 
           daysince=daysince+day-day0
           exit yrsloop
        endif
        daysince=daysince+daysinyear-day0+1
        day0=1
        year0=year0+1
        !           print*,"In loop",year0,year,day0,day,daysince

    end do yrsloop
    !        print*,"Final",year0,year,day0,day,daysince

    daysince=daysince*direction
  end function eudtf2daysince

  function daysince2eudtf(daysince,eudtf0) result (eudtf)
    ! Converts days since a given eudtf date to an eudtf date
    !----args-----!
    integer,intent(in)::daysince,eudtf0
    !---result---!
    integer::eudtf
    !---locals---!
    integer::year,days,year0,day0
    !---Executable-!
    year0=eudtf0/1000
    day0=modulo(eudtf0,1000)
    year=year0
    days=daysince+day0 ! days is now days since beginning of year0
    yrsloop:do 
!        print*,"In yearsloop: days=",days," year=",year
        if (days > days_in_year(year)) then 
           days=days-days_in_year(year)
           year=year+1
        else if (days <= 0)then
           year=year-1
           days=days+days_in_year(year)
        else
           exit yrsloop 
        endif 
    enddo yrsloop
    eudtf=1000*year+days
  end function daysince2eudtf

  function eudtf2cal(eudtf,perm,num,sep) result (cal)
    ! Public function to convert eudtf (yyyyddd) to cal date
    ! Returns a character(len=11) 

    !---------args-----------
    integer,intent(in) :: eudtf
    integer,dimension(:),intent(in),optional::perm
    logical,intent(in),optional::num
    character(len=*),intent(in),optional::sep
    !-------- Result---------
    character(len=11)::cal
    !------Local vars--------
    integer:: year,dayofyear,month,dayofmonth,daysinyear,j,cumul_days,i
    integer,dimension(12)::days_in_month
    character(len=20),dimension(3)::tmpstring
    character(len=1)::sep_char
    character(len=5)::str1,str2
    integer,dimension(3)::order
    logical::num_month
    !--------Executable bit-------------
    !Check for optional args!
    if(present(num)) then
       num_month=num
    else
       num_month=.false.
    endif
    if(present(sep)) then
       sep_char=sep(1:1)
    else
       sep_char=" "
    endif


    ! This is necessary for consistency if perm=3,1,2 or 2,3,1
    ! just doing order=perm is the Wrong Thing here.
    if(present(perm)) then
       do i=1,3
           do j=1,3
               if(perm(j)==i) then
                  order(i)=j
               endif
           enddo
       enddo
    else
       order=(/ 1,2,3 /)
    endif
    year=eudtf/1000
    dayofyear=modulo(eudtf,1000)
    if(year <1) then ! Trap bad year
       call PrintMessage( MLSMSG_Warning,ModuleName,&
            "Module dates_module,function eudtf2cal: year <1" // &
            "I Can not do BC dates. Why on earth do you want one?" ) 
       cal="01 Jan 0001"
       return
    endif
    days_in_month=DAYMAXNY(1:12)
    if ( leapyear(year) ) then
       daysinyear=366
       days_in_month(2)=29 !Sodding February
    else
       daysinyear=365
    endif
    if (dayofyear < 1 .or. dayofyear > daysinyear) then
       write(str1,fmt="(i5)")dayofyear
       dayofyear=max(1,dayofyear)
       dayofyear=min(dayofyear,daysinyear)
       write(str2,fmt="(i5)")dayofyear
       call PrintMessage( MLSMSG_Warning,ModuleName,&
       " in function eudtf2cal: day "//str1//" is out of range."//&
       "Setting it to "//str2 )
    endif
    ! Year and day are now in range
    ! Work out which month this day is in
    cumul_days=0
    do j=1,12
        cumul_days=cumul_days+days_in_month(j)
        if (dayofyear<=cumul_days) then
           cumul_days=cumul_days-days_in_month(j)
           month=j
           exit
        endif
    enddo
    dayofmonth=dayofyear-cumul_days
    ! put the three pieces of the date into 3 elements of tmpstring
    ! in order indicated by order (i.e. by the perm arg if present)
    write(unit=tmpstring(order(1)),fmt="(i2.2)")dayofmonth
    if(num_month)then
       write(unit=tmpstring(order(2)),fmt="(i2.2)")month
    else
       tmpstring(order(2))=monthname(month)(1:3)
    endif
    write(unit=tmpstring(order(3)),fmt=*)year
    ! Stick the three bits into a string to be returned
    cal=trim(adjustl(tmpstring(1)))//sep_char//&
         trim(adjustl(tmpstring(2)))//sep_char//&
         trim(adjustl(tmpstring(3)))
    return
  end function eudtf2cal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !Public function to convert cal dates (dd mmm yyyy) to eudtf (yyyyddd)
  function cal2eudtf(caldate,perm) result (eudtf)
    !------- Arguments ---------!
    character(len=*),intent(in)::caldate
    integer,dimension(:),intent(in),optional::perm
    !------- function result --------!
    integer::eudtf
    !---------Locals---------------!
    character(len=len(caldate))::dpdate ! auto-length string ???OK???
    integer::year,dayofyear,month,dayofmonth,i,j,iw
    integer,dimension(3)::order
    character(len=5)::errstr
    character(len=3)::monthstring
    character(len=50)::instring
    character(len=20),dimension(3)::tmpstring
    character(len=1)::tchar
    integer,dimension(12)::days_in_month
    !--------Executable statements--------!
    ! Get year, month and day out of string
    ! To be robust, we want the spacing to not matter and to allow
    ! things other than blank to be separators
    ! We also need to be able to handle DMY, MDY and YMD dates
    if (present(perm)) then
       order=perm
    else
       order=(/1,2,3/)
    endif
    do j=1,3
        tmpstring(j)=" "
    enddo
    dpdate=depunctuate(caldate) !get rid of punctuation
    instring=" "//adjustl(dpdate)
    iw=0
    l1:do j=1,len_trim(instring)-1
        tchar=instring(j+1:j+1)
        if(tchar /= " " .and. instring(j:j) == " ") then
           iw=iw+1
        endif
        if (iw > 3) then
           call PrintMessage( MLSMSG_Warning,ModuleName,&
           "in fn cal2eudtf: Warning: date"//caldate//" contains >3 words" )
           exit l1
        endif
        tmpstring(order(iw))(j:j)=tchar
    enddo l1
    do j=1,3
        tmpstring(j)=adjustl(tmpstring(j))
    enddo
    read(unit=tmpstring(1),fmt=*)dayofmonth
    read(unit=tmpstring(2),fmt="(a)")monthstring
    read(unit=tmpstring(3),fmt=*)year
    !    print*,"dayofmonth=",dayofmonth,"Month= -->",monthstring,"<--","yr=",year
    if(year <1) then ! Trap bad year
       call PrintMessage( MLSMSG_Warning,ModuleName,&
            "Module dates_module,function cal2eudtf: year <1" // &
            "I Can not do BC dates. Why on earth do you want one?" ) 
       year=0000
       return
    endif

    days_in_month=DAYMAXNY(1:12)
    if ( leapyear(year) ) then
       days_in_month(2)=29 !Sodding February
    endif

    !Get number of this month
    j=verify(monthstring," 0123456789")
    if(j==0) then !month was provided as number
       read(unit=monthstring,fmt=*)month
    else! Month provided as letters
       monthstring=capitalize(monthstring)
       month=0
       mcloop:do i=1,12
           if ( monthstring == capitalize(monthname(i)(1:3)) ) then
              month=i
              exit mcloop
           endif
       enddo mcloop
       if(month==0) then 
          call PrintMessage( MLSMSG_Warning,ModuleName,&
               "in function cal2eudtf: Cannot interpret month "//monthstring )
       endif
    endif
    if (month < 1 .or. month > 12) then
       write(unit=errstr,fmt="(i5)")month
        call PrintMessage( MLSMSG_Warning,ModuleName,&
             "in function cal2eudtf Month "//errstr//" Not valid. "//&
             "Setting it to 1 (==Jan)" )
        month=1
    endif


    dayofyear=sum(days_in_month(1:month-1)) + dayofmonth
    eudtf=year*1000 + dayofyear

  end function cal2eudtf

  ! Converts CCSDS dates from A format (yyyy-mm-dd) to B format (yyyy-DDD)
  function ccsdsa2b(a) result(b)
    !-----Arg--------!
    character(len=*),intent(in)::a
    !-----Function result-----!
    character(len=8)::b
    !------locals--------!
    integer::eudtf
    eudtf=cal2eudtf(a,perm=(/3,2,1/))
    write(unit=b(1:4),fmt="(i4.4)")eudtf/1000
    write(unit=b(6:8),fmt="(i3.3)")modulo(eudtf,1000)
    b(5:5)="-"
  end function ccsdsa2b

  ! Converts CCSDS dates from  B format (yyyy-DDD) to A format (yyyy-mm-dd) 
  function ccsdsb2a(b) result(a)
    !-----Arg--------!
    character(len=*),intent(in)::b
    !-----Function result-----!
    character(len=10)::a
    !------locals--------!
    integer::eudtf,year
    character(len=20)::btr

    btr=adjustl(b)
    read(unit=btr(1:4),fmt="(i4)")year
    read(unit=btr(6:8),fmt="(i3)")eudtf
    eudtf=eudtf+1000*year
    a=eudtf2cal(eudtf,perm=(/3,2,1/),sep="-",num=.true.)

  end function ccsdsb2a

  ! ---------------------------------------------  dai_to_yyyymmdd_ints  -----
  subroutine dai_to_yyyymmdd_ints(dai, yyyy, mm, dd, startingDate)
    ! Routine that given the number of days after a starting date
    ! returns 3 ints: the form yyyymmdd
    !--------Argument--------!
    integer, intent(in)  :: dai
    integer, intent(out) :: yyyy
    integer, intent(out) :: mm
    integer, intent(out) :: dd
    character(len=*), intent(in), optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    integer :: doy1
    integer :: ErrTyp
    integer :: loss
    integer :: mydai
    integer :: gain
    character(len=16) :: mystartingDate
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate=MLSStartingDate
   endif
   call utc_to_yyyymmdd_ints( mystartingDate, ErrTyp, yyyy, mm, dd, nodash=.true. )
   call yyyymmdd_to_doy_str(mystartingDate, doy1)
   if ( dai >= 0 ) then
     mydai = dai
     ! Here's what we do:
     ! Given doy1 (the day-of-year of the starting date)
     ! we keep trying to add dai to it
     ! If the result is greater than the number of days in that year (yyyy),
     ! we increment the starting date's year counter (yyyy), its doy1,
     ! and reduce the dai accordingly and try again
     ! If the result is less than the number of days in that year
     ! Then compute mm and dd for yyyy-(doy1+dai)
     ! print *, 'mydai, yyyy at start ', mydai, yyyy
     do
       if ( mydai + doy1 <= days_in_year(yyyy) ) exit
       ! What we said we'd do
       loss = days_in_year(yyyy) - doy1 + 1
       yyyy = yyyy + 1
       doy1 = 1
       mydai = mydai - loss
     enddo
     ! print *, 'mydai, yyyy at finish ', mydai, yyyy
   else
     ! Since dai is negative, we'll reduce yyyy to compensate
     ! This is the reverse of what we do above, incidentally
     mydai = dai
     ! print *, 'mydai, yyyy at start ', mydai, yyyy
     do
       if ( mydai > 0 ) exit
       gain = days_in_year(yyyy-1) - doy1 + 1
       yyyy = yyyy - 1
       doy1 = 1
       mydai = mydai + gain
     enddo
     ! print *, 'mydai, yyyy at finish ', mydai, yyyy
   endif
   ! Now convert from doy to mmdd
   doy1 = 0 ! How many days into year yyyy added by prior months
   do mm=1, 12
     if ( leapyear(yyyy) ) then
       if ( doy1 + DAYMAXLY(mm) > mydai ) exit
       doy1 = doy1 + DAYMAXLY(mm)
     else
       if ( doy1 + DAYMAXNY(mm) > mydai ) exit
       doy1 = doy1 + DAYMAXNY(mm)
     endif
   enddo
   dd = mydai - doy1 + 1
  end subroutine dai_to_yyyymmdd_ints

  ! ---------------------------------------------  dai_to_yyyymmdd_str  -----
  subroutine dai_to_yyyymmdd_str(dai, str, startingDate)
    ! Routine that given the number of days after a starting date
    ! returns a string: the form yyyymmdd
    !--------Argument--------!
    integer, intent(in)           :: dai
    character(len=*), intent(out) :: str
    character(len=*), intent(in), optional :: startingDate  ! If not Jan 1 2001
    ! Internal variables
    integer :: yyyy, mm, dd
    character(len=8) :: year, month, day
    ! executable
    str = ' '
    call dai_to_yyyymmdd(dai, yyyy, mm, dd, startingDate)
    ! print *, 'dai: ', dai
    ! print *, 'yyyy: ', yyyy
    ! print *, 'mm: ', mm
    ! print *, 'dd: ', dd
    call writeIntsToChars(yyyy , year )
    call writeIntsToChars(mm   , month, fmt='(i2.2)')
    call writeIntsToChars(dd   , day  , fmt='(i2.2)')
    str(1:4) = adjustl(year )
    str(5:6) = adjustl(month)
    str(7:8) = adjustl(day  )
  end subroutine dai_to_yyyymmdd_str

  ! ------------------ dateForm -------------------
  function dateForm(date) result(form)
    ! Determine what format the date is in
    ! E.g., given '2004-d271' returns 'yyyy-doy'
    ! Args
    character(len=*), intent(in) :: date
    character(len=len(date)+8) :: form
    ! Internal variables
    integer :: i
    integer :: j
    integer :: month
    character(len=1)            :: s  ! The expected date field
    ! Executable
    form = 'unknown format'
    if ( len_trim(date) < 1 ) return
    form = ' '
    s = 'y'
    i = 0
    j = 0
    do
      if ( j >= len_trim(date) ) exit
      i = i + 1
      j = j + 1
      select case (date(j:j))
      case ('d')
        form(i:i+2) = 'doy'
        i = i + 3
        j = j + 3
        ! print *, 'After d field: ', form
      case ('J', 'F', 'M', 'A', 'S', 'O', 'N', 'D')
        month = monthNameToNumber(date(j:))
        if ( month < 1 .or. month > 12 ) then
          form = 'month name uncrecognized in ' // trim(date(j:))
          return
        endif
        j = j + len_trim(MONTHNAME(month)) - 1
        form(i:i) = 'M'
        s = 'd'
        ! write(tempFormat(5:6),'(i2.2)') month
        ! print *, 'After M field: ', form
      case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        select case (s)
        case ('m')
          ! Was yyyy, now mm
          form(i:i+1) = 'mm'
          s = 'd'
          i = i + 1
          j = j + 1
          ! print *, 'After 0-9  m field: ', form
        case ('d')
          ! Was mm, now dd
          form(i:i+1) = 'dd'
          s = ' '
          i = i + 1
          j = j + 1
          ! print *, 'After 0-9  d field: ', form
        case ('y')
          ! yyyy
          form(i:i+3) = 'yyyy'
          s = 'm'
          i = i + 3
          j = j + 3
          ! print *, 'After 0-9  y field: ', form
        case default
          ! Huh? Already finished with dd
          ! if ( options%verbose ) print *, 'Unexpected digit in dateForm'
        end select
      case default
        form(i:i) = date(j:j)
        ! print *, 'After default field: ', form
      end select
    enddo
  end function dateForm

  ! ------------------ nextMoon -------------------
  function nextMoon( date, phase ) result( next )
    ! Determine when the next (new) moon will occur
    ! or the next one after optional date
    ! or the occurrence of optional phase where phase may be
    ! 'full'    full moon
    ! 'new'     new moon (default)
    ! 'first'   first quarter
    ! 'last'    last quarter
    
    ! Note:
    ! We assume that the phases repeat every LUNARPERIOD
    ! This is inaccurate
    ! The toolkit has AA component functions/procedures that could
    ! be used to give this information much better

    ! Args
    character(len=*), optional, intent(in) :: date
    character(len=*), optional, intent(in) :: phase
    character(len=32) :: next
    ! Internal variables
    character(len=*), parameter :: dateFormat = 'yyyy-doy'
    character(len=*), parameter :: timeFormat = 'hh:mm:ss'
    character(len=MAXUTCSTRLENGTH) :: myDate
    character(len=8)  :: myPhase
    character(len=16) :: dateString
    type(MLSDate_Time_T)  :: datetime
    double precision :: s, sfirst, snext, sPrev
    character(len=16) :: timeString
    ! Executable
    if ( .not. present(date) ) then
      call date_and_time ( date=dateString, time=timeString )
      dateString = reFormatDate( trim(dateString), toForm=dateFormat )
      timeString = reFormatTime( trim(timeString), timeFormat )
      myDate = trim(dateString) // 'T' // timeString
    else
      myDate = date
    endif
    myPhase = 'new'
    if ( present(phase) ) myPhase = phase
    sFirst = abs(secondsbetween2utcs( '2001-001T00:00:00Z', FIRSTNEWMOON ))
    ! How many seconds since first moon of 2001?
    s = abs(secondsbetween2utcs( FIRSTNEWMOON, myDate ))
    ! When was last new moon? (They recur every LUNARPERIOD)
    sPrev = LUNARPERIOD * int( (s-1.d-3)/LUNARPERIOD )
    select case (myPhase)
    case ('first')
      sNext = sPrev + 0.25*LUNARPERIOD
    case ('full')
      sNext = sPrev + 0.5*LUNARPERIOD
    case ('last')
      sNext = sPrev + 0.75*LUNARPERIOD
    case default
      sNext = sPrev + LUNARPERIOD
    end select
    if ( sNext < s ) sNext = sNext + LUNARPERIOD
    datetime%dai = 0
    datetime%seconds = sFirst + sNext
    ! call dumpDateTime( dateTime, 'Before reducing' )
    call reducedatetime(datetime)
    ! call dumpDateTime( dateTime, 'After reducing' )
    next = datetime2utc(datetime)
    ! print *, 'days before next phase ', daysbetween2utcs( date, next )
    ! print *, 'hours before next phase ', hoursbetween2utcs( date, next )
  end function nextMoon

  ! ------------------ precedesutc -------------------
  function precedesutc( date1, date2 ) result( earlier )
    ! Determine whether date1 came before date2
    ! Notes: 
    ! (1) we assume both date1 and date2 are in utc format
    ! (2): we take pains not convert our args into MLS_Date_Time
    ! because we may have called this function to check if
    ! we must reset TAIStarting

    ! Args
    character(len=*),  intent(in) :: date1
    character(len=*),  intent(in) :: date2
    logical                       :: earlier
    ! Internal variables
    character(len=16) :: hhmmss1
    character(len=16) :: hhmmss2
    integer :: error
    integer :: year1, month1, day1, Doy1
    integer :: year2, month2, day2, Doy2
    ! Executable
    earlier = .true.
    call utc_to_yyyymmdd( date1, error, year1, month1, day1 )
    call utc_to_yyyymmdd( date2, error, year2, month2, day2 )
    ! Obviously if the years are different ..
    if ( year1 < year2 ) then
      return
    elseif ( year2 < year1 ) then
      earlier = .false.
      return
    endif
    ! The years are equal, so we check days-of-year
    call yyyymmdd_to_doy( year1, month1, day1, Doy1 )
    call yyyymmdd_to_doy( year2, month2, day2, Doy2 )
    if ( Doy1 < Doy2 ) then
      return
    elseif ( Doy2 < Doy1 ) then
      earlier = .false.
      return
    endif
    ! The days are equal, so we're down to the times; ugh
    call utc_to_time( date1, error, hhmmss1 )
    call utc_to_time( date2, error, hhmmss2 )
    earlier = ( hhmmss2seconds(hhmmss1) < hhmmss2seconds(hhmmss2) )
  end function precedesutc

  ! ---------------------------------------------  reduceDatetime  -----
  pure subroutine reduceDatetime( datetime, leapsec )
    ! Reduces the seconds field of a datetime to a permissible value
    ! possibly adjusting the dai field in compensation
    ! If leapsec is present, and TRUE, account for leapseconds
    ! Args
    type(MLSDate_Time_T), intent(inout) :: datetime
    logical, optional, intent(in)       :: leapsec
    ! Internal
    integer                             :: dai
    integer                             :: extra  ! Extra days adjustment to dai
    logical                             :: myLeapSeconds
    real                                :: seconds
    ! Executable
    seconds = datetime%seconds
    dai = datetime%dai
    myLeapSeconds = .false.
    if ( present(leapsec) ) myLeapSeconds = leapsec
    ! Leap seconds?
    if ( myLeapSeconds ) then
      datetime%seconds = datetime%seconds - datetime%leapseconds
      datetime%leapseconds = 0 ! (Lest we accidentally do this again)
    endif
    if ( seconds < 0. ) then
      extra = -1 + seconds/SECONDSINADAY
      dai = dai + extra
      seconds = seconds - extra*SECONDSINADAY
    elseif ( seconds >= real(SECONDSINADAY) ) then
      extra = seconds/SECONDSINADAY + 1.e-8
      dai = dai + extra
      seconds = seconds - extra*SECONDSINADAY
      seconds = max( seconds, 0. )
    else
      ! No adjustment necessary
      return
    endif
    datetime%seconds = seconds
    datetime%dai = dai
  end subroutine reduceDatetime

  ! --------------------------------------------------  reFormatDate  -----
  function reFormatDate(date, fromForm, toForm, options) result(reFormat)
    ! Reformat yyyymmdd as yyyy-mm-dd
    ! Allows an optional string toForm defining
    ! the output format; E.g. 'dd M yyyy' for '03 September 2005'
    ! (Thus 'M' expands into the full month name)
    ! or 'yyyy-doy' for '2005-d245' (note the inclusion of the letter 'd')
    ! or 'yyyy-Doy' for '2005-245' (note the absence of the letter 'd')
    ! Another optional string fromForm to hold the input format
    ! in case it wasn't yyyymmdd?
    ! Args
    character(len=*), intent(in)            :: date
    character(len=len(date)+24)             :: reFormat
    character(len=*), optional, intent(in)  :: fromform ! input format
    character(len=*), optional, intent(in)  :: toform   ! output format
    character(len=*), optional, intent(in)  :: options
    ! Internal variables
    character(len=1), dimension(len(date))  :: dateChars
    character(len=1), parameter             :: ymSpacer = '-'
    character(len=1), parameter             :: mdSpacer = '-'
    integer                                 :: i ! format string index
    integer                                 :: j ! date string index
    integer                                 :: k ! aux date string index
    integer                                 :: doy
    character(len=4)                        :: doyString
    character(len=4)                        :: yyyyString
    integer                                 :: month
    character(len=8)                        :: myOptions
    character(len=16)                       :: myFromForm
    character(len=len(date)+24)             :: tempFormat
    logical                                 :: inputWasDoy
    character(len=1), parameter             :: SPACESUBSTITUTE = '?'
    logical                                 :: verbose
    ! Executable
    ! print *, 'date: ', trim(date)
    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = index( myOptions, 'v' ) > 0
    tempFormat = date
    if ( present(fromform) ) then
      ! print *, 'fromForm: ', trim(fromForm)
      ! Treat ' ' and '*' a hints to auto-recognize input format
      if ( len_trim(fromForm) < 1 .or. fromForm == '*' ) then
        myFromForm = dateForm(date)
      else
        myFromForm = fromForm
      endif
      if ( len_trim(myFromForm) > 0 ) then
        inputWasDoy = .false.
        tempFormat = ' '
        i = 0
        j = 0
        do
          if ( i >= len_trim(myFromForm) ) exit
          i = i + 1
          j = j + 1
          ! print *, myFromForm(i:i)
          ! print *, 'i, j', i, j
          select case (myFromForm(i:i))
          case ('y')
            tempFormat(1:4) = date(j:j+3)
            i = i + 3
            j = j + 3
          case ('m')
            tempFormat(5:6) = date(j:j+1)
            i = i + 1
            j = j + 1
          case ('M')
            month = monthNameToNumber( date(j:), abbrev=.true. )
            if ( month < 1 .or. month > 12 ) then
              reFormat = 'month name uncrecognized'
              ! print *, 'month name uncrecognized'
              return
            endif
            ! print *, 'date(j:) ', date(j:)
            ! print *, 'month ', month
            ! No, this next is wrong if we abbreviate the month's name
            ! Instead we need to find the next non-alphabetical char number
            dateChars = transfer( date, dateChars )
            k = FindFirst( .not. isAlphabet( dateChars(j:) ) )
            j = j + k - 2
            ! print *, 'k ', k
            write(tempFormat(5:6),'(i2.2)') month
          case ('D')
            inputWasDoy = .true.
            doyString = 'd' // date(j:j+2)
            i = i + 2
            j = j + 2
          case ('d')
            if ( myFromForm(i:i+1) == 'dd' ) then
              tempFormat(7:8) = date(j:j+1)
              if ( date(j+1:j+1) == ' ' ) tempFormat(7:8) = ' ' // date(j:j)
              ! print *, 'date(j:j+1) ', date(j:j+1)
              i = i + 1
              j = j + 1
            else
              inputWasDoy = .true.
              doyString = date(j:j+3)
              i = i + 2
              j = j + 2
            endif
          case default
            ! i = i + 1
          end select
        enddo
        if ( inputWasDoy ) then
          ! Need to convert from yyyydoy to yyyymmdd
          ! print *, ' Need to convert from yyyydoy to yyyymmdd'
          yyyyString = tempFormat(1:4)
          call yyyydoy_to_yyyymmdd_str(yyyyString, doyString(2:4), tempFormat)
          ! print *, yyyyString, doyString(2:4), tempFormat
        endif
      endif
    endif
    ! print *, 'tempFormat: ', tempFormat
    reFormat = tempFormat(1:4) // ymSpacer // tempFormat(5:6) // mdSpacer // tempFormat(7:8)
    if ( verbose ) then
      print *,  'tempFormat: ' // trim(tempFormat)
      print *, ymSpacer
      print *, mdSpacer
    endif
    if ( .not. present(toform) ) return
    if ( len_trim(toform) < 1 ) return
    reFormat = ' '
    i = 0
    do
      if ( i >= len_trim(toform) ) exit
      i = i + 1
      select case (toform(i:i))
      case ('y')
        reFormat = trim(reFormat) // tempFormat(1:4)
        i = i + 3
      case ('m')
        reFormat = trim(reFormat) // tempFormat(5:6)
        i = i + 1
      case ('M')
        read(tempFormat(5:6), *) month
        reFormat = trim(reFormat) // trim(MONTHNAME(month))
      case ('D')
        call yyyymmdd_to_doy_str(tempFormat, doy)
        write(doyString, '(i3.3)') doy
        reFormat = trim(reFormat) // doyString
        i = i + 2
      case ('d')
        if ( toform(i:i+1) == 'dd' ) then
          reFormat = trim(reFormat) // tempFormat(7:8)
          i = i + 1
        else
          ! print *, 'tempFormat: ', tempFormat
          call yyyymmdd_to_doy_str(tempFormat, doy)
          write(doyString, '(a1, i3.3)') 'd', doy
          reFormat = trim(reFormat) // doyString
          i = i + 2
        endif
      case (' ')
        ! Foolish me--it can't handle spaces because of all the trims
        reFormat = trim(reFormat) // SPACESUBSTITUTE
      case default
        reFormat = trim(reFormat) // toform(i:i)
      end select
    enddo
    if ( index(reFormat, SPACESUBSTITUTE) < 1 ) return
    ! Need to substitute a space for every occurrence of SPACESUBSTITUTE
    do
      i = index(reFormat, SPACESUBSTITUTE)
      if ( i < 1 ) return
      reFormat(i:i) = ' '
    enddo
  end function reFormatDate

  ! --------------------------------------------------  reFormatTime  -----
  function reFormatTime( time, form ) result( reFormat )
    ! Reformat hhmmss.sss as hh:mm:ss
    ! (Note it truncates instead of rounding)
    ! "form" is an optional string arg defining
    ! the output format; E.g. 'hh:mm' for '13:23' (note 24-hour time)
    ! or 'HH:mm' for '01:23 PM' (note AM/PM marking)
    ! 's' fills the remaining chars with the seconds plus any decimal portion
    ! 'S' fills just one character at a time, so 'hh:mm:SS' would be '12:30:00'
    ! while 'hh:mm:ss' could be '12:30:00.0000000'
    ! Args
    character(len=*), intent(in)            :: time
    character(len=len(time)+24)             :: reFormat
    character(len=*), optional, intent(in)  :: form
    ! Internal variables
    character(len=1), parameter             :: hmSpacer = ':'
    character(len=1), parameter             :: msSpacer = ':'
    integer                                 :: hours
    integer                                 :: i, j
    character(len=2)                        :: ampm
    character(len=2)                        :: hh
    character(len=2)                        :: hc, mc
    character(len=len(time))                :: sc
    ! Executable
    ! In case we came here with ':' separators
    if ( index(time(1:4), ':') > 0 ) then
       reFormat = time
    else
      reFormat = time(1:2) // hmSpacer // time(3:4) // msSpacer // time(5:)
    endif
    if ( .not. present(form) ) return
    if ( len_trim(form) < 1 ) return
    hc = reFormat(1:2)
    mc = reFormat(4:5)
    sc = adjustl(reFormat(7:))
    ampm = ' '
    reFormat = ' '
    i = 0
    j = 0
    do
      if ( i >= len_trim(form) ) exit
      i = i + 1
      select case (form(i:i))
      case ('h')
        reFormat = trim(reFormat) // hc
        i = i + 1
      case ('H')
        read(hc, *) hours
        ampm = 'AM'
        if ( hours > 12 ) then
          hours = hours - 12
          ampm = 'PM'
        endif
        write(hh, '(i2.2)') hours
        reFormat = trim(reFormat) // hh
        i = i + 1
      case ('m')
        reFormat = trim(reFormat) // mc
        i = i + 1
      case ('S')
        j = j + 1
        reFormat = trim(reFormat) // sc(j:j)
      case ('s')
        reFormat = trim(reFormat) // sc
        i = i + 1
      case default
        reFormat = trim(reFormat) // form(i:i)
      end select
    enddo
    if ( ampm /= ' ' ) reFormat = trim(reFormat) // ' ' // ampm
  end function reFormatTime

  ! ------------------ timeForm -------------------
  function timeForm(time) result(form)
    ! Determine what format the time is in
    ! E.g., given '01:23 PM' returns 'HH:mm'
    ! Note: may be confused if any field has fewer than 2 didits
    ! E.g., will wrongly identify format of  "5:18 AM" as "HHmm"
    ! So must use "05:18 AM"
    ! Args
    character(len=*), intent(in) :: time
    character(len=len(time)+8) :: form
    ! Internal variables
    character(len=2) :: hh
    integer          :: i, j
    character(len=2) :: s
    ! Executable
    form = ' '
    if ( len_trim(time) < 1 ) return
    if ( any ( indexes(lowercase(time), (/'am', 'pm'/) ) > 0 ) ) then
      hh = 'HH'
    else
      hh = 'hh'
    endif
    i = 0
    j = 0
    s = 'h' ! We always begin with the hours field, then minutes, then seconds
    do
      if ( j >= len_trim(time) ) exit ! So we don't loop forever
      i = i + 1
      j = j + 1
      select case ( lowercase(time(j:j)) )
      case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.')
        select case (s)
        case ('h')
          ! 'hh' or 'HH'
          form(i:i+1) = hh
          s = 'm'
          i = i + 1
          j = j + 1
          ! print *, 'After 0-9  m field: ', form
        case ('m')
          ! Was hh, now mm
          form(i:i+1) = 'mm'
          s = 's'
          i = i + 1
          j = j + 1
          ! print *, 'After 0-9  m field: ', form
        case ('s')
          ! Was mm, now dd
          form(i:i+1) = 'ss'
          s = ' '
          i = i + 1
          j = j + 1
          ! print *, 'After 0-9  d field: ', form
        case default
          ! Probably just continuing with ss.sss field
        end select
      case ('a', 'p', 'm', ' ')
        ! Probably just the 'AM/PM' mark
      case ('z')
        ! An optional terminator
      case default
        form(i:i) = time(j:j) ! Expecting ':' between fields
      end select
    end do
  end function timeForm

  ! ---------------------------------------------  secondsbetween2utcs  -----
  function secondsbetween2utcs( utc1, utc2, leapsec ) result(seconds)
    ! Find how many seconds separate two utcs
    ! Optionally accound for leap seconds
    ! Args
    character(len=*), intent(in) :: utc1, utc2
    logical, optional, intent(in) :: leapsec
    double precision             :: seconds
    ! Internal
    type(MLSDate_Time_T)         :: datetime1, datetime2
    integer                      :: days
    ! Executable
    datetime1 = utc2datetime( utc1 )
    datetime2 = utc2datetime( utc2 )
    ! print *, 'utc1, utc2 ', utc1, ' ', utc2
    ! call dump (datetime1)
    ! call dump (datetime2)
    days = datetime2%dai - datetime1%dai
    ! print *, 'days ', days
    seconds = 24.d0*60*60*days + ( datetime2%seconds - datetime1%seconds )
    if ( .not. present(leapsec) ) return
    if ( leapsec ) seconds = seconds + ( datetime2%leapseconds - datetime1%leapseconds )
  end function secondsbetween2utcs

  ! ---------------------------------------------  secondsindaydble  -----
  function secondsindaydble( utc1 ) result(seconds)
    ! Given a utc return how many seconds since its day began
    ! Args
    character(len=*), intent(in) :: utc1
    double precision             :: seconds
    ! Internal
    type(MLSDate_Time_T)         :: datetime1
    ! Executable
    datetime1 = utc2datetime(utc1)
    seconds = datetime1%seconds
  end function secondsindaydble

  ! ---------------------------------------------  splitDateTime  -----
  subroutine splitDateTime(str, ErrTyp, date, time, strict)
    ! Routine that returns the date, time portions from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! (See also utc_to_date and utc_to_time)
    !--------Argument--------!
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: ErrTyp
    character(len=*), intent(out) :: date
    character(len=*), intent(out) :: time
    logical,intent(in), optional  :: strict
    !----------Local vars----------!
    logical :: mystrict
    !----------Executable part----------!

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   date = ' '
   time = ' '

   if(len_trim(str) <= 0) then
      if(mystrict) then
         ErrTyp=INVALIDUTCSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDUTCSTRING
   ! Snip off time fields from date fields
   call GetStringElement(lowercase(str), date, 1, &
     & countEmpty=.true., inseparator='t')
   if ( date == 't' ) then
     if ( .not. mystrict) Errtyp = 0
     date = ' '
     return
   endif
   call GetStringElement(lowercase(str), time, 2, &
     & countEmpty=.true., inseparator='t')
   if ( time == 't' ) then
     if ( .not. mystrict) Errtyp = 0
     time = ' '
     return
   endif
   ErrTyp=0
   
  end subroutine splitDateTime

  ! ---------------------------------------------  FromUARSdate  -----
  ! Converts yyyy-mm-dd From uars date format, e.g. 'd0007'
  subroutine FromUARSdate( uars, utc, leapsec )
    character(len=*), intent(in)   :: uars
    character(len=*), intent(out)  :: utc
    logical, optional, intent(in)  :: leapsec
    !----local -----!
    type(MLSDate_Time_T)           :: datetime
    type(MLSDate_Time_T)           :: fiducialdatetime
    integer                        :: dai
    character(len=16)              :: daiStr
    integer                        :: ErrTyp
    character(len=64)              :: temp
    character(len=16)              :: timeStr
    ! Executable
    call splitDateTime( uars, ErrTyp, daiStr, timeStr )
    ! Now we must find the offset in days, knowing that 
    ! 1991d261 was 'd0007'
    fiducialdatetime = utc2datetime( '1991d261' )
    ! call Dump ( fiducialdatetime, '1991d261 fiducial' )
    
    ! Check that we can invert utc2datetime
    temp = datetime2utc( fiducialdatetime )
    ! print *, 'Inverting, the utc for fiducialdatetime is ', trim(temp)
    call readIntsFromChars ( daiStr, dai, ignore='*' )
    ! print *, 'daiStr, dai ', daiStr, dai
    datetime%dai = dai - 7 + fiducialdatetime%dai
    ! call Dump( dateTime, 'daiStr convrted' )
    temp = datetime2utc(datetime)
    call splitDateTime( temp, ErrTyp, utc, daiStr )
    ! print *, 'utc: ', utc
    ! print *, 'daiStr: ', daiStr
    if ( len_trim(timeStr) > 0 ) utc = trim(utc) // 'T' // trim(timeStr)
  end subroutine FromUARSdate

  ! ---------------------------------------------  toUARSdate  -----
  ! Converts yyyy-mm-dd to uars date format, e.g. 'd0007'
  subroutine toUARSdate( utc, uars, leapsec )
    character(len=*), intent(in)   :: utc
    character(len=*), intent(out)  :: uars
    logical, optional, intent(in)  :: leapsec
    !----local -----!
    type(MLSDate_Time_T)           :: datetime
    type(MLSDate_Time_T)           :: fiducialdatetime
    integer                        :: dai
    integer                        :: ErrTyp
    character(len=16)              :: daiStr
    character(len=16)              :: timeStr
    character(len=16)              :: rtSiad
    ! integer :: THEDAY
    ! Executable
    dateTime = utc2datetime( utc )
    call splitDateTime( utc, ErrTyp, daiStr, timeStr )
    ! call dump ( dateTime )
    ! Now we must find the offset in days, knowing that 
    ! 1991d261 was 'd0007'
    fiducialdatetime = utc2datetime( '1991d261' )
    dai = 7 + datetime%dai - fiducialdatetime%dai
    ! The following trick should result in a string exactly 5 chars long
    ! 'd' + '0..0' + dai
    call WriteIntsToChars ( dai, daiStr )
    daiStr = '00000000' // trim(daiStr)
    rtSiad = reverse_trim( daiStr )
    uars = 'd' // reverse_trim( rtSiad(1:4) )
    if ( len_trim(timeStr) > 0 ) uars = trim(uars) // 'T' // trim(timeStr)
  end subroutine toUARSdate

  ! ---------------------------------------------  utc_to_date  -----
  subroutine utc_to_date(str, ErrTyp, date, &
    & strict, utcAt0z)
    ! Routine that returns the date portion from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! The date substring in subdivided by the separator '-'
    ! into either two or three fields
    ! In case (A), the 3 fields are year, month, and day of month
    ! in case (B) the two fields are year and day of year
    
    ! For case (A) returns year, month, and day=day of month
    ! For case (B) returns year, month=-1, and day=day of year
    ! Useful to decode utc inputs into attribute values
    
    ! Optionally returns the input string in utcAt0z modified so that 
    ! the hh:mm:ss.sss is 00:00:00Z
    
    ! (See also utc_to_time and utc_to_yyymmdd)
    !--------Argument--------!
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: ErrTyp
    character(len=*), intent(out) :: date
    logical,intent(in), optional  :: strict
    character(len=*),intent(out), optional   :: utcAt0z
    !----------Local vars----------!
    logical :: mystrict
    character(len=*), parameter :: chars_0z = 'T00:00:00Z'
    !----------Executable part----------!

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(len_trim(str) <= 0) then
      if(mystrict) then
         ErrTyp=INVALIDUTCSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDUTCSTRING
   ! Snip off time fields from date fields
   call GetStringElement(lowercase(str), date, 1, &
     & countEmpty=.true., inseparator='t')
   if ( date == ' ' ) then
     if ( .not. mystrict) Errtyp = 0
     if ( present(utcAt0z) ) utcAt0z = ' '
     return
   endif
   if ( present(utcAt0z) ) utcAt0z = trim(date) // chars_0z
   ErrTyp=0
   
  end subroutine utc_to_date

  ! ---------------------------------------------  utc_to_time  -----
  subroutine utc_to_time(str, ErrTyp, time, strict)
    ! Routine that returns the time portion from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! (See also utc_to_date and utc_to_yyyymmdd)
    !--------Argument--------!
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: ErrTyp
    character(len=*), intent(out) :: time
    logical,intent(in), optional  :: strict
    !----------Local vars----------!
    logical :: mystrict
    integer :: zpos
    !----------Executable part----------!

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(len_trim(str) <= 0) then
      if(mystrict) then
         ErrTyp=INVALIDUTCSTRING
      else
         ErrTyp=0
      endif
      return
   endif
   
   ErrTyp=INVALIDUTCSTRING
   ! Separate off time fields from date fields
   call GetStringElement(lowercase(str), time, 2, &
     & countEmpty=.true., inseparator='t')
   if ( time == 't' ) then
     if ( .not. mystrict) Errtyp = 0
     time = ' '
     return
   endif

   ! Snip off terminal 'Z' (and everything beyond)
   zpos = index(time, 'z')
   if ( zpos > 0 ) then
     time(zpos:) = ' '
   endif
   ErrTyp=0
   
  end subroutine utc_to_time

  ! ---------------------------------------------  utc_to_yyyymmdd_ints  -----
  subroutine utc_to_yyyymmdd_ints(str, ErrTyp, year, month, day, strict, nodash)
    ! Routine that returns the year, month, and day from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! (N) yyyydddThh:mm:ss.sss (no-dash)
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! The date substring in subdivided by the separator '-'
    ! into either two or three fields
    ! In case (A), the 3 fields are year, month, and day of month
    ! in case (B) the two fields are year and day of year
    
    ! For case (A) returns year, month, and day=day of month
    ! For case (B) returns year, month=-1, and day=day of year
    ! Useful to decode utc inputs into attribute values
    
    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Argument--------!
    character(len=*),intent(in) :: str
    integer, intent(out) :: ErrTyp
    integer, intent(out) :: year
    integer, intent(out) :: month
    integer, intent(out) :: day
    logical, intent(in), optional :: strict
    logical, intent(in), optional :: nodash   ! No dash separating date fields
    !----------Local vars----------!
    logical, parameter :: countEmpty = .false.
    character(len=1), parameter :: dash='-'
    character(len=NameLen) :: date
    character(len=NameLen) :: yyyy
    character(len=NameLen) :: mm
    character(len=NameLen) :: dd
    character(len=NameLen) :: doy
    character(LEN=*), parameter :: time_conversion='(I4)'
    logical :: mystrict
    logical :: mynodash
    character(len=1) :: utc_format        ! 'a' or 'b' or 'n' (no-dash)
    ! The following arrys contains the maximum permissible day for each month
    ! where month=-1 means the whole year, month=1..12 means Jan, .., Dec
    integer, dimension(-1:12), parameter :: DAYMAX = (/ &
      & 366, 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 &
      & /)
    !----------Executable part----------!

   year = -1
   month = -1
   day = -1
   mm = ' '

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   if(present(nodash)) then
      mynodash=nodash
   else
      mynodash=.false.
   endif
         
   call utc_to_date(str, ErrTyp, date, strict= .true.)
   ! mynodash = mynodash .or. index( date, ' -' ) < 1
   
   ! A trick! If date is non-numerical, e.g. '(undefined)', we
   ! will process it as if it was yyyy-mm-dd
   if ( .not. isDigits(str(1:4) ) ) mynodash = .false.
   if ( mynodash ) mynodash = index( date, '-' ) < 1
   if ( str(5:5) == 'd' ) mynodash = .true.

   ! print *, ' '
   ! print *, 'str, mystrict, mynodash, date, ErrTyp ', &
   !   & trim(str), ' ', mystrict, mynodash, ' ', trim(date), ' ', ErrTyp
   
   if ( ErrTyp /= 0 ) then
     if ( .not. mystrict ) ErrTyp = 0
     return
   endif
   if ( myNoDash ) then
     utc_format = 'n'
     yyyy = date(1:4)
     if ( index(date, 'd') < 1 ) then
       mm = date(5:6)
       dd = date(7:8)
     else
       doy = date(6:8)
       call yyyydoy_to_yyyymmdd_str( yyyy, doy, date )
       mm = date(5:6)
       dd = date(7:8)
     endif
   else
     call GetStringElement(trim(date), yyyy, 1, countEmpty=countEmpty, inseparator=dash)
     if ( &
       & NumStringElements(trim(date), countEmpty=countEmpty, inseparator=dash) == 2) then
       call GetStringElement(trim(date), dd, 2, countEmpty=countEmpty, inseparator=dash)
       utc_format = 'b'
     else
       call GetStringElement(trim(date), mm, 2, countEmpty=countEmpty, inseparator=dash)
       call GetStringElement(trim(date), dd, 3, countEmpty=countEmpty, inseparator=dash)
       utc_format = 'a'
     endif
   endif
   
   ErrTyp=0
   
   ! print *, 'utc_format, yyyy, mm, dd ', &
   !   & trim(utc_format), ' ', trim(yyyy), ' ', trim(mm), ' ', trim(dd)
   ! Convert to value
   if(yyyy /= ' ') then
      read(yyyy, time_conversion, iostat=ErrTyp) year
   endif
   
   if(ErrTyp /= 0) then
      return
   elseif(year < 0 .or. year > YEARMAX) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif

   if(mm /= ' ') then
      read(mm, time_conversion, iostat=ErrTyp) month
   endif

   if(utc_format == 'b') then
     ErrTyp = 0
     month = -1
   elseif(ErrTyp /= 0) then
      return
   elseif(month < 1 .or. month > 12) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif
   ! Coming out of the above, month should be in the interval [-1, 12]

   if(dd /= ' ') then
      read(dd, time_conversion, iostat=ErrTyp) day
   endif
   ! print *, 'month, day', month, day

   if(ErrTyp /= 0) then
      return
   elseif(day < 1 .or. day > DAYMAX(month)) then
      ErrTyp=INVALIDUTCSTRING
      return
   endif
  end subroutine utc_to_yyyymmdd_ints

  ! ------------ utc2tai93s ------
  ! Function returns time since midnight Jan 1 1993 in s
  ! Optionally accounts for leap seconds
  ! (see also secondsbetween2utcs)
  function utc2tai93s ( utc, leapsec ) result ( tai93s )
    character(len=*), intent(in) :: utc
    logical, optional, intent(in) :: leapsec
    double precision             :: tai93s
    type(MLSDate_Time_T)         :: datetime
    integer :: eudtf
    datetime = utc2datetime( utc )
    ! call dump( datetime )
    ! However, our internal datetime object stores days since 2001,
    ! so we need to offset it
    eudtf = daysince2eudtf( datetime%dai, 2001001 )
    tai93s = datetime%seconds + 24*3600*eudtf2daysince( eudtf, 1993001 )
    if ( .not. present(leapsec) ) return
    if ( leapsec ) &
      & tai93s = tai93s + HowManyleapSeconds( &
      & eudtf2daysince( eudtf, 1993001 ), start='19930101' &
      & )
  end function utc2tai93s
  
  ! ---------------------------------------------  utc_to_yyyymmdd_strs  -----
  subroutine utc_to_yyyymmdd_strs(str, ErrTyp, year, month, day, &
    & strict)
    ! Routine that returns the year, month, and day from a string of the form
    ! (A) yyyy-mm-ddThh:mm:ss.sss
    ! (B) yyyy-dddThh:mm:ss.sss
    ! where the field separator 'T' divides the string into two
    ! sub-strings encoding the date and time
    ! The date substring in subdivided by the separator '-'
    ! into either two or three fields
    ! In case (A), the 3 fields are year, month, and day of month
    ! in case (B) the two fields are year and day of year
    
    ! For case (A) returns year, month, and day=day of month
    ! For case (B) returns year, month=-1, and day=day of year
    ! Useful to decode utc inputs into attribute values
    
    ! (See also PGS_TD_UTCtoTAI and mls_UTCtoTAI)
    !--------Argument--------!
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: ErrTyp
    character(len=*), intent(out) :: year
    character(len=*), intent(out) :: month
    character(len=*), intent(out) :: day
    logical,intent(in), optional  :: strict
    !----------Local vars----------!
    character(len=1), parameter :: dash='-'
    character(len=NameLen) :: date
    logical :: mystrict
    !----------Executable part----------!

   year = ' '
   month = '-1'
   day = ' '

   if(present(strict)) then
      mystrict=strict
   else
      mystrict=.false.
   endif
         
   call utc_to_date(str, ErrTyp, date, strict= .true.)
   ! print *, 'date: ', trim(date)
   ! print *, 'ErrTyp: ', ErrTyp
   if ( ErrTyp /= 0 ) then
     if ( .not. mystrict ) ErrTyp = 0
     return
   endif
   call GetStringElement(trim(date), year, 1, countEmpty=.true., inseparator=dash)
   if ( &
     & NumStringElements(trim(date), countEmpty=.true., inseparator=dash) == 2) then
     call GetStringElement(trim(date), day, 2, countEmpty=.true., inseparator=dash)
   else
     call GetStringElement(trim(date), month, 2, countEmpty=.true., inseparator=dash)
     call GetStringElement(trim(date), day, 3, countEmpty=.true., inseparator=dash)
   endif
   ! print *, 'num: ', NumStringElements(trim(date), countEmpty=.true., inseparator=dash)
   ! print *, 'utc_format: ', utc_format
   
   ErrTyp=0
   
  end subroutine utc_to_yyyymmdd_strs

  ! ------------------ utcForm ---------------------
  function utcForm(utc) result(which)
    ! Determine which of recognized formats is utc
    ! 'a', 'b', 'n' (no dashes) or 'u' (unrecognized)
    ! Arg
    character(len=*), intent(in) :: utc
    character                    :: which
    ! Internal variables
    character(len=1), parameter :: dash='-'
    character(len=NameLen) :: date
    integer                :: ErrTyp, numDashes
    ! Executable
    which = 'u'
    call utc_to_date(utc, ErrTyp, date, strict= .true.)
    ! print *, 'ErrTyp: ', ErrTyp
    if ( ErrTyp /= 0 ) return
    numDashes = -1 + &
      &  NumStringElements( trim(date), countEmpty=.true., inseparator=dash )
    ! print *, 'numDashes: ', numDashes
    select case (numDashes)
    case (0)
      which = 'n'
    case (1)
      which = 'b'
    case (2)
      which = 'a'
    case default
      which = 'u'
    end select
  end function utcForm

  ! ---------------------------------------------  yyyymmdd_to_dai_pure  -----
  elemental subroutine yyyymmdd_to_dai_pure( yyyymmdd, dai, startingDate )
    ! Routine that returns the number of days after a starting date
    ! from yyyymmdd taking great pains to be elemental
    !--------Argument--------!
    character(len=*), intent(in) :: yyyymmdd
    integer,intent(out) :: dai
    character(len=*),intent(in), optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    character(len=16) :: mystartingDate
    integer :: yyyy1, mm1, dd1, doy1 ! Based on input yyyymmdd
    integer :: yyyy2 !, doy2
    integer :: yr
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate=MLSStartingDate
   endif
   yyyy1 = 1000*chartoint(yyyymmdd(1:1)) + 100*chartoint(yyyymmdd(2:2)) + &
     & 10*chartoint(yyyymmdd(3:3)) + chartoint(yyyymmdd(4:4))
   mm1 = 10*chartoint(yyyymmdd(5:5)) + chartoint(yyyymmdd(6:6))
   if ( len_trim(yyyymmdd) > 7 ) then
     dd1 = 10*chartoint(yyyymmdd(7:7)) + chartoint(yyyymmdd(8:8))
   else
     dd1 = chartoint(yyyymmdd(7:7))
   endif
   doy1 = DayOfYear( yyyy1, mm1, dd1 )

   yyyy2 = 1000*chartoint(mystartingDate(1:1)) + 100*chartoint(mystartingDate(2:2)) + &
     & 10*chartoint(mystartingDate(3:3)) + chartoint(mystartingDate(4:4))
   mm1 = 10*chartoint(mystartingDate(5:5)) + chartoint(mystartingDate(6:6))
   if ( len_trim(mystartingDate) > 7 ) then
     dd1 = 10*chartoint(mystartingDate(7:7)) + chartoint(mystartingDate(8:8))
   else
     dd1 = chartoint(mystartingDate(7:7))
   endif
   ! doy2 = DayOfYear( yyyy2, mm1, dd1 )
   dai = doy1
   if ( yyyy2 < yyyy1 ) then
    do yr = yyyy2, yyyy1 - 1
      dai = dai + days_in_year ( yr )
    enddo
   elseif ( yyyy2 > yyyy1 ) then
     ! Uh-oh, looks like dai will have to be negative
    do yr = yyyy1 - 1, yyyy2
      dai = dai - days_in_year ( yr )
    enddo
   ! elseif ( yyyy2 == yyyy1 ) then
   !   dai = doy1 
   endif
  end subroutine yyyymmdd_to_dai_pure

  ! ---------------------------------------------  yyyymmdd_to_dai_ints  -----
  subroutine yyyymmdd_to_dai_ints(yyyy, mm, dd, dai, startingDate)
    ! Routine that returns the number of days after a starting date
    ! from 3 ints: the form yyyymmdd
    !--------Argument--------!
    integer ,intent(in) :: yyyy
    integer ,intent(in) :: mm
    integer ,intent(in) :: dd
    integer,intent(out) :: dai
    character(len=*),intent(in),optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    character(len=16) :: mystartingDate
    integer :: yyyy1, mm1, dd1, doy1
    integer :: yyyy2, doy2
    integer :: ErrTyp
    logical :: daiNegative
    integer :: y
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate=MLSStartingDate
   endif
   call utc_to_yyyymmdd_ints(mystartingDate, ErrTyp, yyyy1, mm1, dd1, nodash=.true.)
   call yyyymmdd_to_doy_str(mystartingDate, doy1)
   call yyyymmdd_to_doy_ints(yyyy, mm, dd, doy2)
   ! print *, 'doy1, doy2, yyyy1, yyyy ', doy1, doy2, yyyy1, yyyy
   ! print *, 'days between ', yyyy1, doy1, ' and ', yyyy, mm, dd, doy2
   yyyy2 = yyyy
   daiNegative = yyyy1 > yyyy2
   if ( daiNegative ) then
     call switch(yyyy1, yyyy2)
     call switch(doy1, doy2)
   elseif ( yyyy1 == yyyy2 ) then
     dai = doy2 - doy1
     return
   endif
   dai = doy2 - doy1
   do y = yyyy1, yyyy2 - 1
     if ( leapyear(y) ) then
       dai = dai + DAYMAXLY(-1)
     else
       dai = dai + DAYMAXNY(-1)
     endif
   enddo
   if ( daiNegative ) dai = -dai
  end subroutine yyyymmdd_to_dai_ints

  ! ---------------------------------------------  yyyymmdd_to_dai_str  -----
  subroutine yyyymmdd_to_dai_str(str, dai, startingDate)
    ! Routine that returns the number of days after a starting date
    ! from a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*),intent(in) :: str
    integer, intent(out) :: dai
    character(len=*),intent(in),optional :: startingDate  ! If not Jan 1 2001
    !----------Local vars----------!
    character(len=16) :: mystartingDate
    integer :: yyyy1, mm1, dd1, doy1
    integer :: yyyy2, mm2, dd2, doy2
    integer :: ErrTyp
    logical :: daiNegative
    integer :: y
    !----------Executable part----------!
   if(present(startingDate)) then
      mystartingDate=startingDate
   else
      mystartingDate=MLSStartingDate
   endif
   call utc_to_yyyymmdd_ints(mystartingDate, ErrTyp, yyyy1, mm1, dd1, nodash=.true.)
   call utc_to_yyyymmdd_ints(str, ErrTyp, yyyy2, mm2, dd2, nodash=.true.)
   call yyyymmdd_to_doy_str(mystartingDate, doy1)
   call yyyymmdd_to_doy_str(str, doy2)
   daiNegative = yyyy1 > yyyy2
   if ( daiNegative ) then
     call switch(yyyy1, yyyy2)
     call switch(doy1, doy2)
   elseif ( yyyy1 == yyyy2 ) then
     dai = doy2 - doy1
     return
   endif
   dai = doy2 - doy1
   do y = yyyy1, yyyy2 - 1
     if ( leapyear(y) ) then
       dai = dai + DAYMAXLY(-1)
     else
       dai = dai + DAYMAXNY(-1)
     endif
   enddo
   if ( daiNegative ) dai = -dai
  end subroutine yyyymmdd_to_dai_str
!! ------------------------------------------------------------------- !!
  ! Private procedures
  logical function isDateTimeInRange ( dateTime, lower, upper ) result( itIs )
    ! Is the date part of dateTime within an allowable range?
    ! Args
    type(MLSDate_Time_T), intent(in)             :: datetime
    type(MLSDate_Time_T), intent(in), optional   :: lower
    type(MLSDate_Time_T), intent(in), optional   :: upper
    ! Executable
    if ( present(lower) .and. present(upper) ) then
      itIs = ( dateTime%dai >= lower%dai .and. dateTime%dai <= upper%dai )
      itis = itIs .and. dateTime%seconds >= lower%seconds .and. &
        & dateTime%seconds <= upper%seconds
    else
      ! How far back in years do we think we should be allowed to retreat?
      ! And how far in the future do we think we will ever need to advance?
      ! Quite arbitrarily, let's say no more than 500 years each way.
      ! which, by the way, is more than 1.5 *10^10 seconds; log_2 of
      ! that is more than 33, meaning we will have blown past short
      ! 32-bit intger reprresentation before that point is reached
      itIs = abs( dateTime%dai ) < 125*1461  ! 1461 is 4*365.25
    endif
  end function isDateTimeInRange

  elemental subroutine yyyyMondd_yyyymmdd_str( yyyyMondd, yyyymmdd )
    ! Routine that returns the string yyyymmdd
    ! from a string of the form yyyyMondd
    !--------Argument--------!
    character(len=*), intent(in)  :: yyyyMondd
    character(len=*), intent(out) :: yyyymmdd
    !----------Local vars----------!
    character(len=16) :: yyyy, Mon, dd, mm
    !----------Executable part----------!
    call SplitWords( yyyyMondd, yyyy, Mon, dd, &
       & threeWay=.true., delimiter=' ' )
    if ( len_trim(dd) < 2 ) then
      dd = '0' // dd
    elseif ( len_trim(dd) < 1 ) then
      dd = '01'
    endif
    select case ( lowercase(Mon(1:3)) )
    case ( 'jan' )
      mm = '01'
    case ( 'feb' )
      mm = '02'
    case ( 'mar' )
      mm = '03'
    case ( 'apr' )
      mm = '04'
    case ( 'may' )
      mm = '05'
    case ( 'jun' )
      mm = '06'
    case ( 'jul' )
      mm = '07'
    case ( 'aug' )
      mm = '08'
    case ( 'sep' )
      mm = '09'
    case ( 'oct' )
      mm = '10'
    case ( 'nov' )
      mm = '11'
    case ( 'dec' )
      mm = '12'
    case default
      mm = '01'
    end select
    yyyymmdd = yyyy(1:4) // mm(1:2) // dd(1:2)
  end subroutine yyyyMondd_yyyymmdd_str

  ! ---------------------------------------------  yyyymmdd_to_doy_ints  -----
  subroutine yyyymmdd_to_doy_ints(year, month, day, doy)
    ! Routine that returns the number of days after the year's start
    ! for year, month, day
    !--------Argument--------!
    integer, intent(in) :: year, month, day
    integer, intent(out) :: doy
    !----------Local vars----------!
    integer :: m
    integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     if ( year < 0 .or. year > YEARMAX ) then
       doy = -1
     endif
     doy = day
     if ( month <= 1 ) then
       return
     endif
     if ( leapyear(year) ) then
       DAYMAX = DAYMAXLY
     else
       DAYMAX = DAYMAXNY
     endif
     do m=1, month-1
       doy = doy + DAYMAX(m)
     enddo
     
  end subroutine yyyymmdd_to_doy_ints

  ! ---------------------------------------------  yyyymmdd_to_doy_str  -----
  subroutine yyyymmdd_to_doy_str( str, doy )
    ! Routine that returns the number of days after the year's start
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*), intent(in) :: str
    integer, intent(out) :: doy
    !----------Local vars----------!
    integer :: year, month, day
    character(len=2) :: dayCh, monCh
    character(len=4) :: yearCh
    !----------Executable part----------!
     doy = -1
     if ( len_trim(str) < 8 ) return
     ! print *, 'str: ', str
     ! In case we came here with any illegal characters in the fields
     if ( index(str(1:8), '*') > 0 ) return
     ! In case we came here with '-' separators between the fields
     if ( index(str(1:8), '-') > 0 ) then
       call GetStringElement( str, yearCh, 1, .true., '-' )
       call GetStringElement( str, monCh, 2, .true., '-' )
       call GetStringElement( str, dayCh, 3, .true., '-' )
       read(yearCh, *) year
       read(monCh, *) month
       read(dayCh, *) day
     else
       read(str(1:4), *) year
       read(str(5:6), *) month
       read(str(7:8), *) day
     endif
     call yyyymmdd_to_doy_ints(year, month, day, doy)
  end subroutine yyyymmdd_to_doy_str

  ! ---------------------------------------------  yyyydoy_to_yyyymmdd_str  -----
  subroutine yyyydoy_to_yyyymmdd_str(yyyy, doy, yyyymmdd)
    ! Routine that converts the string encoding the number of days 
    ! after the year's start, input as yyyydoy,
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    character(len=*),intent(in) :: yyyy
    character(len=*),intent(in) :: doy ! w/o the letter 'd'
    character(len=*),intent(out) :: yyyymmdd
    !----------Local vars----------!
    integer :: year, doynum
    !----------Executable part----------!
     read(yyyy, *) year
     read(doy, *) doynum
     call yeardoy_to_yyyymmdd_ints(year, doynum, yyyymmdd)
  end subroutine yyyydoy_to_yyyymmdd_str

  ! ---------------------------------------------  yeardoy_to_yyyymmdd_ints  -----
  subroutine yeardoy_to_yyyymmdd_ints(year, doy, yyyymmdd)
    ! Routine that converts the string encoding the number of days 
    ! after the year's start, input as yyyydoy,
    ! for a string of the form yyyymmdd
    !--------Argument--------!
    integer,intent(in)           :: year
    integer,intent(in)           :: doy
    character(len=*),intent(out) :: yyyymmdd
    !----------Local vars----------!
    integer :: day
    integer :: daysum
    integer :: m
    integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     if ( year < 0 .or. year > YEARMAX ) then
       yyyymmdd = 'year not in range'
     endif
     daysum = 0
     if ( leapyear(year) ) then
       DAYMAX = DAYMAXLY
     else
       DAYMAX = DAYMAXNY
     endif
     do m=1, 12
       if ( daysum + DAYMAX(m) >= doy ) exit
       daysum = daysum + DAYMAX(m)
     enddo
     if ( m > 12 ) then
       yyyymmdd = 'doy not in range'
       return
     endif
     day = doy - daysum
     write(yyyymmdd,'(I4.4, 2i2.2)') year, m, day
  end subroutine yeardoy_to_yyyymmdd_ints

  ! ---------------------------------------------  yyyyDoy_to_mmdd_ints  -----
  subroutine yyyyDoy_to_mmdd_ints(year, month, day, doy)
    ! Routine that converts the number of days 
    ! after the year's start, input as yyyy and doy,
    ! into the month and day of that month
    !--------Argument--------!
    integer,intent(in)           :: year
    integer,intent(in)           :: doy
    integer,intent(out)          :: month
    integer,intent(out)          :: day
    !----------Local vars----------!
    integer :: daysum
    integer, dimension(-1:12) :: DAYMAX
    !----------Executable part----------!
     month = 0
     day = 0
     daysum = 0
     if ( leapyear(year) ) then
       DAYMAX = DAYMAXLY
     else
       DAYMAX = DAYMAXNY
     endif
     do month=1, 12
       if ( daysum + DAYMAX(month) >= doy ) exit
       daysum = daysum + DAYMAX(month)
     enddo
     if ( month > 12 ) then
       return
     endif
     day = doy - daysum
  end subroutine yyyyDoy_to_mmdd_ints

  ! ---------------------------------------------------  DayNumberOfWeek  -----
  integer function DayNumberOfWeek(dai)
    integer,intent(in) :: dai
     ! Days in the week are numbered American-style
     ! I.e., 'Sunday' is 1, .., Saturday is 7
     !
     ! Now we need to know what day of the week January 1 2001 was
     ! Let's assume it was 'Monday", which corresponds to 2
     integer, parameter :: STARTINGDAYOFWEEK = 2
     DayNumberOfWeek = mod( STARTINGDAYOFWEEK + dai - 1, 7 ) + 1
  end function DayNumberOfWeek

  ! ---------------------------------------------------  DayOfYear  -----
  elemental integer function DayOfYear( year, month, day )
    integer, intent(in) :: year, month, day
     ! Given year, month, and day return day-of-year
     integer :: mon
     DayOfYear = day
     do mon = 1, 12
       if ( mon >= month ) return
       DayOfYear = DayOfYear + daysInMonth( mon, year )
     enddo
  end function DayOfYear

  ! ---------------------------------------------------  HowManyleapSeconds  -----
  elemental integer function HowManyleapSeconds( dai, Start )
    integer,intent(in) :: dai
    character(len=*), optional, intent(in) :: Start
    ! How many leap seconds have been added since January 1 2001 until dai
    ! If Start present, use it instead of January 1 2001
    
    ! Bugs and limitations:
    ! If dai is before January 1 2001 (or whatever starting date we use) returns 0
    ! If Start is present, we don't check for its validity
    ! If there are more leap second dates declared after the last leapSecDates
    ! you must add it to the array and recompile
    ! Variables
    integer :: i
    integer, dimension(size(leapSecDates)) :: leapSecDais
    character(len=16) :: yyyymmdd
    ! First, calculate dai for each leap second date
    do i =1, size(leapSecDates)
      ! Note that reFormatDate can't be made elemental due to its use of 
      ! read/write to internal files
      ! Someday we must repair that
      ! yyyymmdd = reFormatDate( leapSecDates(i), &
      ! & fromForm='yyyy MMM dd', toForm='yyyymmdd' )
     call yyyyMondd_yyyymmdd_str( leapSecDates(i), yyyymmdd )
     ! For the same reason, yyyymmdd_to_dai can't be made elemental, either
     !  (Sigh--all this hackery-quackery because of something so simple?)
     call yyyymmdd_to_dai_pure( yyyymmdd, leapSecDais(i), start )
    enddo
    ! Method: for each +tive leap second date, compare it with dai
    ! If it was prior to dai, add another leap second
    ! If it was after, the current leap second total is the final value
    ! If the leap second date is < 0, then it occurred before January 1 2001
    ! (or whatever starting time we are using), so we must not add
    ! another leap second
    HowManyleapSeconds = 0
    do i =1, size(leapSecDates)
      if ( leapSecDais(i) > dai ) return
      if ( leapSecDais(i) > 0 ) HowManyleapSeconds = HowManyleapSeconds + 1
    enddo
  end function HowManyleapSeconds

  ! ---------------------------------------------------  Leapyear  -----
  elemental logical function leapyear(year)
    integer,intent(in) :: year
     ! This is to capture the rule that centuries are leap only
     ! if divisible by 400
     ! Otherwise, as perhaps you knew, leapyears are those years divisible by 4
     leapyear = mod(year,4) == 0 .and. & ! Processor might short-circuit this
       & ( (mod(year,100) /= 0) .or. (mod(year,400) == 0) )
  end function leapyear

  ! ------------------------------------------  monthNameToNumber  -----
  function monthNameToNumber( name, abbrev ) result(number)
    ! Convert month name to corresponding number
    ! E.g., given 'March', returns 3
    ! As a courtesy, name may be case-insensitive
    ! As a further courtesy, name may be followed by any junk you like
    ! Thus 'March 23, 2004 01:59:59.999' still returns 3
    ! If no such month is found, returns -1
    ! Options:
    ! If abbrev is present and TRUE, check just the first 3 characters of name
    ! for a match
    ! Args
    character(len=*), intent(in)             :: name
    integer                                  :: number
    logical, optional :: abbrev
    ! Local variables
    logical :: myAbbrev
    ! Executable
    myAbbrev = .false.
    if ( present(abbrev) ) myAbbrev = abbrev
    do number=1, size(MONTHNAME)
      if ( index(lowerCase(name), lowercase(trim(MONTHNAME(number)))) > 0 ) return
      if ( myAbbrev .and. &
        & index(lowerCase(name), lowercase(trim(MONTHNAME(number)(1:3)))) > 0 ) return
    enddo
    number = -1
  end function monthNameToNumber

  ! ---------------------------------------------  switch_ints  -----
  subroutine switch_ints(x1, x2)
    ! Switch args x1 <=> x2
    !--------Argument--------!
    integer,intent(inout) :: x1
    integer,intent(inout) :: x2
    !----------Local vars----------!
    integer :: x
    x = x1
    x1 = x2
    x2 = x
  end subroutine switch_ints

  ! ---------------------------------------------  hhmmss2seconds  -----
  function hhmmss2seconds(time) result(seconds)
    ! Given a time return the number of seconds after midnight
    ! Args
    character(len=*), intent(in)        :: time ! E.g., "11:20:33.000"
    real                                :: seconds
    ! Internal
    integer                             :: hours
    integer                             :: minutes
    real                                :: sdotss
    ! Executable
    hours = 0
    minutes = 0
    sdotss = 0.
    read( time(1:2), '(i2)' ) hours
    if ( len_trim(time) > 3 ) read( time(4:5), '(i2)' ) minutes
    if ( len_trim(time) > 6 ) read( time(7:), * ) sdotss
    seconds = sdotss + 60.*( minutes + 60*hours )
  end function hhmmss2seconds

  ! ---------------------------------------------  s2hhmmss  -----
  function s2hhmmss(s) result(hhmmss)
    ! Args
    double precision, intent(in)        :: s
    character(len=16)       :: hhmmss
    ! Internal
    integer                 :: hours
    integer                 :: minutes
    real                    :: seconds
    character(len=2)        :: hh, mm
    character(len=16)       :: ss
    ! Executable
    if ( s <= 0. ) then
      hhmmss = '00:00:00'
      return
    endif
    hours = s / (60*60)
    minutes = (s - 60*60*hours) / 60
    seconds = s - 60.*(minutes + 60*hours)
    seconds = max(seconds, 0.)
    write(hh, '(i2.2)') hours
    write(mm, '(i2.2)') minutes
    write(ss, '(f12.9)') seconds
    if ( seconds < 10.) then
      hhmmss = hh // ':' // mm // ':0' // adjustl(ss)
    else
      hhmmss = hh // ':' // mm // ':' // adjustl(ss)
    endif
  end function s2hhmmss

  ! ---------------------------------------------  utc2datetime  -----
  function utc2datetime(arg) result(datetime)
    ! Given a utc return an MLSDate_Time_T
    ! Args
    character(len=*), intent(in) :: arg
    type(MLSDate_Time_T)         :: datetime
    ! Internal
    integer                      :: day
    integer                      :: ErrTyp
    character(len=16)            :: hhmmss
    character(len=32)            :: utc
    integer                      :: month
    integer                      :: year
    ! Executable
    ! print *, 'utc: ', trim(arg)
    ! Beware of formats like yyyymmdd
    utc = arg
    if ( index(arg, '-') < 1 ) then
      if ( isDigits(arg) ) utc = arg(1:4) // '-' // arg(5:6) // '-' // arg(7:8)
    endif
    call utc_to_yyyymmdd_ints( utc, ErrTyp, year, month, day )
    call utc_to_time( utc, ErrTyp, hhmmss )
    ! print *, 'year, month, day: ', year, month, day
    ! print *, 'hhmmss: ', trim(hhmmss)
    ! Now convert to our internal representations
    call yyyymmdd_to_dai_ints( year, month, day, datetime%dai )
    datetime%seconds = hhmmss2seconds(hhmmss)
    ! call dump(datetime, 'From utc2datetime')
    dateTime%leapSeconds = HowManyleapSeconds ( datetime%dai )
  end function utc2datetime

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module dates_module
! $Log$
! Revision 2.46  2019/09/19 15:59:47  pwagner
! daysInMonth returns 0 if month out of range
!
! Revision 2.45  2019/07/09 23:12:55  pwagner
! Made MonthName public
!
! Revision 2.44  2018/12/11 01:21:14  pwagner
! No longer uses Printit_M
!
! Revision 2.43  2017/07/10 18:37:08  pwagner
! Correct error in dai_to_yyyymmdd_ints that arose when dai < 0; added GetStartingDate
!
! Revision 2.42  2016/09/14 00:23:14  pwagner
! Added isUTCInRange so we can check for invalid utc values
!
! Revision 2.41  2016/08/08 22:56:11  pwagner
! Added leap second to 2016
!
! Revision 2.40  2016/07/28 01:37:04  vsnyder
! Cannonball polishing
!
! Revision 2.39  2015/06/25 23:17:16  pwagner
! One more go at fixing this (is uars so important?)
!
! Revision 2.38  2015/06/24 18:00:22  pwagner
! A workaround, fixing gold brick, until better solution found
!
! Revision 2.37  2015/06/23 23:53:42  pwagner
! Can canvert to and from uars dates
!
! Revision 2.36  2015/02/27 23:10:06  pwagner
! Corrected tai93s2datetime; added gethid
!
! Revision 2.35  2015/01/21 00:49:47  pwagner
! Added 2015 leap second
!
! Revision 2.34  2014/10/14 21:38:57  pwagner
! Added datesbetween2utcs
!
! Revision 2.33  2014/03/05 20:14:50  pwagner
! Commented-out stray prints; improved some remarks
!
! Revision 2.32  2014/02/21 01:26:54  pwagner
! Added leapseconds; corrected bugs
!
! Revision 2.31  2014/02/13 00:01:24  pwagner
! Repaired bugs in Reformatting dates containing month names
!
! Revision 2.30  2013/08/28 00:38:17  pwagner
! Added a local version of PrintMessage to evade possible circular dependency
!
! Revision 2.29  2013/08/21 00:21:58  pwagner
! Added a function to tell whether one utcdate precedes another
!
! Revision 2.28  2013/08/14 17:25:15  pwagner
! Added reset and restore StartingDates to handle dates before Jan 1 1993
!
! Revision 2.27  2013/06/18 23:02:22  pwagner
! Removed more unused stuff
!
! Revision 2.26  2013/04/05 23:19:22  pwagner
! Added tai93s2hid
!
! Revision 2.25  2011/08/02 16:50:00  honghanh
! Add function utc2tai93s
!
! Revision 2.24  2011/04/26 20:56:16  pwagner
! Can convert tai93s to other date formats
!
! Revision 2.23  2011/01/04 00:46:49  pwagner
! hoursinday and secondsInDay now public functions
!
! Revision 2.22  2010/06/28 16:59:35  pwagner
! Added more comments about numerical date types
!
! Revision 2.21  2009/08/26 16:20:13  pwagner
! Corrected utcForm function--was reversing 'a' and 'b'
!
! Revision 2.20  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.19  2009/06/16 17:23:38  pwagner
! Added tai2ccsds
!
! Revision 2.18  2007/12/19 01:28:46  pwagner
! Removed unused variables
!
! Revision 2.17  2007/09/24 20:23:32  pwagner
! Added new procedures converting to/from yyyyDoy
!
! Revision 2.16  2007/09/20 17:39:15  pwagner
! Added daysInMonth
!
! Revision 2.15  2007/09/14 00:14:50  pwagner
! Added buildCalendar
!
! Revision 2.14  2007/09/06 22:27:59  pwagner
! Added nextMoon
!
! Revision 2.13  2007/07/23 23:19:53  pwagner
! Added many new procedures
!
! Revision 2.12  2007/01/18 19:37:38  pwagner
! New addDaysToUTC and addHoursToUTC functions
!
! Revision 2.11  2006/09/29 00:27:20  pwagner
! Added dateForm to auto-recognize date formats
!
! Revision 2.10  2005/09/23 20:43:45  pwagner
! Removed a few redundancies, a few
!
! Revision 2.9  2005/09/22 23:33:58  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.8  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2003/04/04 13:05:24  hcp
! Added a number of comments. In particular, made it clearer that the TAI returned
! by routines in here is in DAYS not SECONDS
!
! Revision 2.6  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.5  2002/10/01 22:04:19  pwagner
! Renamed moduleNameIn to ModuleName
!
! Revision 2.4  2002/10/01 18:40:35  bwknosp
! Added RCS Ident Info Block
!
