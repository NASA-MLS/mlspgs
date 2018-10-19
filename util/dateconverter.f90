! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program dateconverter
!=================================

   use Dates_Module, only: AddDaysToUTC, AddHoursToUTC, AddSecondsToUTC, &
     & Dateform, DayOfWeek, FromUARSDate, HoursInDay, PrecedesUTC, &
     & ReformatDate, ResetStartingDate, SecondsInDay, SplitDateTime, &
     & Tai93s2utc, ToUARSDate
   use Machine, only: Hp, Getarg
   use MLSStringLists, only: ExpandStringRange
   use MLSStrings, only: Lowercase, Ncopies, ReadNumsFromChars

   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! This program converts an input date to a different format

! E.g., given "2004d274" it will output "2004 September 30"
! Useful only because possibly one will use doy 
! while the other uses month-and-day
! Useful also in automatically computing gmao entries in PCF

  type options_T
    integer     :: offset = 0                  ! How many days to add/subtract
    integer     :: hoursOffset = 0             ! How many hours to add/subtract
    integer     :: weekdayLength= 0     ! how many chars to output day-of-week
    double precision :: secondsOffset = 0.d0   ! How many seconds to add/subtract
    logical     :: leapsec = .false.           ! account for leap seconds?
    logical     :: utcFormat = .false.         ! output date+time?
    logical     :: compare = .false.           ! compare pairs of dates
    logical     :: debug = .false.
    logical     :: verbose = .false.
    character(len=255) :: outputFormat= ' '    ! output format
    character(len=255) :: inputFormat= ' '     ! input format
    character(len=255) :: inputTime= ' '       ! input time
    character(len=255) :: argRange= ' '        ! which date numbers to do
  end type options_T

  type ( options_T ) :: options

! Variables
   integer, parameter :: MAXLISTLENGTH=24
   integer, parameter ::          MAXDATES = 100
   character (len=2)             :: comparison
   character (len=MAXLISTLENGTH) :: converted_date
   character (len=MAXLISTLENGTH) :: converted_time
   character (len=MAXLISTLENGTH) :: date
   character(len=MAXLISTLENGTH), dimension(MAXDATES) :: dates
   logical, dimension(MAXDATES)  :: doThisDate
   character(len=*), parameter   :: DOYFORMAT = 'yyyy-doy'
   integer                       :: ErrTyp
   character (len=MAXLISTLENGTH) :: fromForm
   double precision              :: hours
   integer                       :: i
   character (len=MAXLISTLENGTH) :: intermediate_date
   ! character (len=*), parameter  :: intermediateForm = 'yyyymmdd'
   character(len=*), parameter   :: MFORMAT = 'yyyy M dd'
   integer                       :: n_dates = 0
   double precision              :: seconds
   double precision              :: secondsperday = 24*3600.
   double precision              :: tai
   character (len=MAXLISTLENGTH) :: time
   character (len=MAXLISTLENGTH) :: toForm
   character (len=MAXLISTLENGTH) :: uars_date
   character (len=16           ) :: weekday
  ! Executable
  do      ! Loop over options
     call get_date(date, n_dates, options)
     if ( index('-+', date(1:1)) > 0 ) cycle
     if ( date == ' ' ) exit
     if ( options%inputFormat /= ' ' ) then
       ! Because the internal dateForm will not be needed
     elseif ( index(dateForm(date), 'yyyy') == 0 ) then
       print *, 'Sorry--date format not found: ', trim(date)
       cycle
     endif
     n_dates = n_dates + 1
     dates(n_dates) = date
  enddo
  doThisDate = .true.
  time = options%inputTime
  if ( options%inputTime == ' ' ) time = '00:00:00'
  converted_time = time
  if ( options%argRange /= ' ' ) &
    & call ExpandStringRange ( options%argRange, doThisDate )
  do i=1, n_dates
    if ( .not. doThisDate(i)) cycle
    date = dates(i)
    if ( options%verbose ) print *, 'Input was ', trim(date)
    fromForm = options%inputFormat
    if ( len_trim(options%inputFormat) < 1 ) fromForm = dateForm(date)
    if ( index(lowercase(options%inputFormat), 'uars' ) > 0 ) then
      call resetStartingdate( newMLSDate='1980-01-01' )
      call splitDateTime ( date, ErrTyp, uars_date, converted_time )
      call FromUARSDate( uars_date, intermediate_date )
      date = intermediate_date
      if ( options%debug ) print *, 'uars_date: ', uars_date
      if ( options%debug ) print *, 'date: ', date
    elseif ( options%inputFormat == 'tai' ) then
      call readNumsFromChars ( date, tai )
      date = tai93s2utc( tai, options%leapsec )
      fromForm = 'yyyy-mm-dd'
      options%utcFormat = .true.
      call splitDateTime ( date, ErrTyp, intermediate_date, converted_time )
    endif
    ! Figure out logical format to convert it to
    if ( index(lowercase(options%outputFormat), 'utc' ) > 0 ) then
      toForm = 'yyyy-Doy'
      if ( index(lowercase(options%outputFormat), 'b' ) > 0 ) toForm = 'yyyy-mm-dd'
      options%utcFormat = .true.
    elseif ( index(lowercase(options%outputFormat), 'uars' ) > 0 ) then
      ! Convert the date to its uars format; e.g. 'd0007'
      toForm = 'uars'
    elseif ( options%outputFormat /= ' ' ) then
      toForm = options%outputFormat
    elseif ( index(lowercase(fromForm), 'doy') > 0 ) then
      toForm = MFORMAT
    else
      toForm = DOYFORMAT
    endif
    if ( options%verbose ) then
      print *, 'Input format was ', trim(fromForm)
      print *, 'We parsed it as ', trim(date)
      print *, 'Output format is ', trim(toForm)
      print *, 'input time is ', trim(options%inputTime)
      print *, 'Offset days is ', options%offset
      print *, 'Offset hours is ', options%hoursoffset
      print *, 'Offset seconds is ', options%secondsoffset
    endif
    dates(i) = date    
  enddo
  ! Are we comparing pairs of dates?
  if ( options%compare ) then
    do i=1, n_dates, 2
      ! A slight complication:
      ! if the two dates are equal, then neither precedes the other
      ! +: date1 < date2
      ! -: date2 < date1
      ! A hint of hackery-quackery now follows
      ! We set the comparison string to one of 4 values
      ! value           meaning
      ! -----           -------
      ! (blank)      date1 = date2
      !    +         date1 < date2
      !    -         date2 < date1
      !    +-        date1 = date2
      comparison = ' '
      if ( PrecedesUTC( dates(i), dates(i+1) ) ) comparison = '+'
      if ( PrecedesUTC( dates(i+1), dates(i) ) ) comparison = trim(comparison) // '-'
      if ( len_trim(comparison) /= 1 ) comparison = '='
      call print_string( trim(comparison) )
    enddo
    stop
  endif
  do i=1, n_dates
    if ( .not. doThisDate(i)) cycle
    date = dates(i)
    ! Process date
    ! First: three special codes
    if ( index(options%outputFormat, 'sec') > 0 ) then
      if ( index(date, 'T') < 1 ) date = 'T' // date
      seconds = secondsinday(date)
      write(*,'(f9.1, " s")') seconds
      cycle
    elseif ( index(lowercase(options%inputFormat), 'uars' ) > 0 ) then
      ! Convert the date from its uars format; e.g. 'd0007'
      ! Have you tested this yet?
      call resetStartingdate( newMLSDate='1980-01-01' )
      converted_date = reFormatDate(date, &
        & fromForm='yyyy-mm-dd', toForm=trim(toForm))
    elseif ( index(lowercase(options%outputFormat), 'uars' ) > 0 ) then
      ! Convert the date to its uars format; e.g. 'd0007'
      call resetStartingdate( newMLSDate='1980-01-01' )
      call splitDateTime ( date, ErrTyp, converted_date, converted_time )
      intermediate_date = reFormatDate( converted_date, &
        & fromForm=trim(fromForm), toForm='yyyymmdd' )
      call ToUARSDate( intermediate_date, converted_date )
      if ( len_trim(converted_time) > 0 ) &
        & converted_date = trim(converted_date) // 'T' // adjustl(Time)
    elseif ( index(options%outputFormat, 'hour') > 0 ) then
      if ( index(date, 'T') < 1 ) date = 'T' // date
      hours = hoursinday(date)
      write(*,'(f9.4, " h")') hours
      cycle
    elseif ( options%offset /= 0 ) then
    ! We will always have an intermediate_date in yyyy-mm-dd format
      converted_date = reFormatDate( date, &
        & fromForm=trim(fromForm), toForm='yyyy-mm-dd' )
      intermediate_date = trim(converted_date) // 'T' // &
        & adjustl(Time)
      if ( options%debug ) print *, 'intermediate_date', intermediate_date
      converted_date = adddaystoutc( intermediate_date, options%offset )
      if ( options%debug ) print *, 'intermediate_date (advanced)', converted_date
      call splitDateTime( converted_date, ErrTyp, intermediate_date, converted_time )
      converted_date = reFormatDate(intermediate_date, &
        fromForm='yyyy-mm-dd', toForm=trim(toForm))
    elseif ( options%hoursoffset /= 0 ) then
      converted_date = reFormatDate( date, &
        & fromForm=trim(fromForm), toForm='yyyy-mm-dd' )
      if ( options%debug ) print *, '1st date: ', trim(converted_date)
      ! if ( options%inputTime == ' ' ) options%inputTime = '00:00:00'
      ! Mash date and time together in unholy union
      intermediate_date = trim(converted_date) // 'T' // &
        & adjustl(Time)
      if ( options%debug ) print *, 'mash date: ', trim(intermediate_date)
      converted_date = addhourstoutc( intermediate_date, options%hoursoffset )
      if ( options%verbose ) print *, 'advanced date: ', trim(converted_date)
      ! Now split them asunder
      call splitDateTime( converted_date, ErrTyp, intermediate_date, converted_time )
      converted_date = reFormatDate(intermediate_date, &
        fromForm='yyyy-mm-dd', toForm=trim(toForm))
    elseif ( options%secondsoffset /= 0 ) then
      converted_date = reFormatDate( date, &
        & fromForm=trim(fromForm), toForm='yyyy-mm-dd' )
      if ( options%debug ) print *, '1st date: ', trim(converted_date)
      ! if ( options%inputTime == ' ' ) options%inputTime = '00:00:00'
      ! Mash date and time together in unholy union
      intermediate_date = trim(converted_date) // 'T' // &
        & adjustl(Time)
      if ( options%debug ) print *, 'mash date: ', trim(intermediate_date)
      converted_date = addsecondstoutc( intermediate_date, options%secondsoffset )
      if ( options%verbose ) print *, 'advanced date: ', trim(converted_date)
      ! Now split them asunder
      call splitDateTime( converted_date, ErrTyp, intermediate_date, converted_time )
      converted_date = reFormatDate(intermediate_date, &
        fromForm='yyyy-mm-dd', toForm=trim(toForm))
    else
      intermediate_date = reFormatDate( date, &
        & fromForm=trim(fromForm), toForm='yyyy-mm-dd' )
      if ( options%debug ) print *, 'intermediate: ' // trim(intermediate_date)
      converted_date = reFormatDate(intermediate_date, &
        & fromForm='yyyy-mm-dd', toForm=trim(toForm))
    endif
    if ( options%weekdayLength > 0 ) then
      weekday = dayOfWeek( intermediate_date, fromForm='yyyy-mm-dd' )
      converted_date = weekday(1:options%weekdayLength) // ' ' // converted_date
    endif
    if ( options%utcFormat ) then
      call print_string( trim(converted_date) // 'T' // trim(converted_time) // 'Z' )
    elseif ( options%inputTime /= ' ' ) then
      call print_string( trim(converted_date) // ' ' // trim(converted_time) )
    else
      call print_string( trim(converted_date) )
    endif
   enddo
contains
!------------------------- get_date ---------------------
    subroutine get_date(date, n_dates, options)
    ! Added for command-line processing
     character(len=*), intent(out) :: date          ! date
     integer, intent(in)             :: n_dates
     type ( options_T ), intent(inout) :: options
     ! Local variables
     character(len=255) ::              arg
     integer ::                         error = 1
     integer, save ::                   i = 1
  ! Get date, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, date )
      ! print *, i, ' th Arg: ', trim(date)
      error = 0
      if ( index('-+', date(1:1)) < 1 ) exit
      if ( date(1:3) == '-h ' ) then
        call print_help
      elseif ( date(1:3) == '-H ' ) then
        call getarg ( i+1+hp, arg )
        call readNumsFromChars(arg, options%hoursOffset)
        options%hoursOffset = - options%hoursOffset ! Because we will subtract
        i = i + 1
        exit
      elseif ( date(1:3) == '+H ' ) then
        call getarg ( i+1+hp, arg )
        call readNumsFromChars(arg, options%hoursOffset)
        i = i + 1
        exit
      elseif ( date(1:3) == '-R ' ) then
        call getarg ( i+1+hp, arg )
        call resetStartingdate( arg )
        i = i + 1
        exit
      elseif ( date(1:3) == '-S ' ) then
        call getarg ( i+1+hp, arg )
        read(arg, * ) options%secondsOffset
        options%secondsOffset = - options%secondsOffset ! Because we will subtract
        i = i + 1
        exit
      elseif ( date(1:3) == '+S ' ) then
        call getarg ( i+1+hp, arg )
        read(arg, * ) options%secondsOffset
        i = i + 1
        exit
      elseif ( date(1:4) == '-arg' ) then
        call getarg ( i+1+hp, options%argRange )
        i = i + 1
        exit
      elseif ( date(1:3) == '-i ' ) then
        call getarg ( i+1+hp, options%inputFormat )
        i = i + 1
        exit
      elseif ( date(1:3) == '-t ' ) then
        call getarg ( i+1+hp, options%inputTime )
        i = i + 1
        exit
      elseif ( date(1:3) == '-o ' ) then
        call getarg ( i+1+hp, options%outputFormat )
        i = i + 1
        exit
      elseif ( date(1:5) == '-leap' ) then
        options%leapsec = .true.
        exit
      elseif ( date(1:3) == '-c ' ) then
        options%compare = .true.
        exit
      elseif ( date(1:3) == '-d ' ) then
        options%debug = .true.
        options%verbose = .true.
        exit
      elseif ( date(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( date(1:3) == '-w ' ) then
        options%weekdayLength = 2
        exit
      elseif ( date(1:2) == '-W' ) then
        options%weekdayLength = ncopies( date, 'W' )
        exit
      elseif ( index('0123456789', date(2:2)) > 0 ) then
        call readNumsFromChars(date, options%offset)
        exit
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    if (trim(date) == ' ' .and. n_dates == 0) then

    ! Last chance to enter date
      ! print *,  "Enter the date you wish to convert."
      read(*,'(a)') date
    endif
    
  end subroutine get_date
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:dateconverter [options] [dates]'
      write (*,*) &
      & ' If no dates supplied, you will be prompted to supply one'
      write (*,*) &
      & ' If no time-of-day supplied, none will be output'
      write (*,*) ' Options:'
      write (*,*) ' -o format   => output format to use (e.g. yyyymmdd)'
      write (*,*) '               by default output will complement input'
      write (*,*) '               e.g., "2004 October 01" <=> 2004-d275'
      write (*,*) '               special codes:'
      write (*,*) '               uars means uars dates'
      write (*,*) '                    e.g., 1991-11-01 becomes d0051'
      write (*,*) '               utc fuses date and time'
      write (*,*) '                    e.g., 2007-274T23:59:59.9999Z'
      write (*,*) '               sec prints seconds-in-day'
      write (*,*) '                    e.g., 2007-274T00:01:59.99Z'
      write (*,*) '                    prints 119.99'
      write (*,*) '               hour prints hours-in-day'
      write (*,*) '                    e.g., 2007-274T14:40:00Z'
      write (*,*) '                    prints 14.66667'
      write (*,*) ' -i format   => input format'
      write (*,*) '               by default attempt to auto-recognize'
      write (*,*) '               special codes:'
      write (*,*) '               if format is "tai", treat input as '
      write (*,*) '               (double-precision) seconds since'
      write (*,*) '               1993-01-01T00:00:00'
      write (*,*) ' -R date     => Reset starting date to "date"'
      write (*,*) '               (needed if any dates prior to Jan 1 1993)'
      write (*,*) ' -t time     => optional time-of-day (military-style)'
      write (*,*) '               e.g., "06:53:10"'
      write (*,*) ' -number     => subtract "number" days'
      write (*,*) ' +number     => add "number" days'
      write (*,*) ' -H number   => subtract "number" hours'
      write (*,*) ' +H number   => add "number" hours'
      write (*,*) ' -S number   => subtract "number" seconds'
      write (*,*) ' +S number   => add "number" seconds'
      write (*,*) ' -arg  range =>'
      write (*,*) '       Run just the dates indexed by the expression range'
      write (*,*) '       e.g., 7 means run only the 7th of the dates, '
      write (*,*) '           1,3-9+2,12 means run dates 1,3,5,7,9,12'
      write (*,*) ' -leapsec    => account for leap seconds; e.g. for "tai"'
      write (*,*) ' -c          => compare date1 with date2:'
      write (*,*) '                  condition         returns'
      write (*,*) '                date1 < date2          +'
      write (*,*) '                date2 < date1          -'
      write (*,*) '                date1 = date2          ='
      write (*,*) ' -d          => switch on debug mode'
      write (*,*) ' -v          => switch on verbose mode'
      write (*,*) ' -w          => print day-of-week in 2 characters'
      write (*,*) ' -WW..       => print day-of-week'
      write (*,*) '                in count[W] characters'
      write (*,*) ' -h          => print brief help'
      stop
  end subroutine print_help
!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!==================
end program dateconverter
!==================

! $Log$
! Revision 1.13  2018/05/22 23:25:33  pwagner
! Last update broke all cases except -c usage; ffixed
!
! Revision 1.12  2018/05/04 16:38:11  pwagner
! Added new commandline option -c comparing 2 dates
!
! Revision 1.11  2017/03/30 23:35:29  pwagner
! Add explanation of special uars formatting to help page
!
! Revision 1.10  2015/06/24 18:03:01  pwagner
! Fix some bugs related to uars input
!
! Revision 1.9  2015/06/23 23:54:16  pwagner
! Can canvert to and from uars dates
!
! Revision 1.8  2014/03/05 20:16:04  pwagner
! New -leapsec option accounts for leap seconds
!
! Revision 1.7  2013/08/14 17:26:38  pwagner
! Added -R comdline option to handle dates before Jan 1 1993
!
! Revision 1.6  2011/01/04 00:53:04  pwagner
! Among other improvements, can now convert times to {hours,seconds}-in-day
!
! Revision 1.5  2010/06/28 17:04:13  pwagner
! Added 'tai' format to convert l2gp%time field
!
! Revision 1.4  2010/06/03 23:36:53  pwagner
! Added option to print day-of-week
!
! Revision 1.3  2007/08/17 00:44:30  pwagner
! May add seconds offset; rely more completely on dates module
!
! Revision 1.2  2007/05/10 23:41:42  pwagner
! Fixed error in declared char size of date; now assumed shape
!
! Revision 1.1  2007/01/18 23:36:04  pwagner
! Moved here from sandbox
!
