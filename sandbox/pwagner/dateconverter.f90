!=================================
PROGRAM dateconverter
!=================================

   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   USE MLSStrings , ONLY: lowerCase, reFormatDate

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program converts an input date to a different format

! E.g., given "2004d274" it will output "2004 September 30"
! Useful only because possibly one will use doy 
! while the other uses month-and-day

  type options_T
    logical     :: verbose = .false.
    character(len=255) :: outputFormat= ' '        ! output format
  end type options_T

  type ( options_T ) :: options
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
  character(len=*), dimension(12), parameter :: MONTHNAME = (/ &
    & 'January  ', 'February ', 'March    ', 'April    ', 'May      ', &
    & 'June     ', 'July     ', 'August   ', 'September', 'October  ', &
    & 'November ', 'December '/)

   INTEGER, PARAMETER :: MAXLISTLENGTH=24
   integer, parameter ::          MAXDATES = 100
	CHARACTER (LEN=MAXLISTLENGTH) :: converted_date
	CHARACTER (LEN=MAXLISTLENGTH) :: date
   character(len=MAXLISTLENGTH), dimension(MAXDATES) :: dates
	CHARACTER (LEN=MAXLISTLENGTH) :: toForm
	CHARACTER (LEN=MAXLISTLENGTH) :: fromForm
   character(len=*), parameter   :: MFORMAT = 'yyyy M dd'
   character(len=*), parameter   :: DOYFORMAT = 'yyyy-doy'
   integer                       :: i
   integer                       :: n_dates = 0
  ! Executable
  do      ! Loop over options
     call get_date(date, n_dates, options)
     if ( date(1:1) == '-' ) cycle
     if ( date == ' ' ) exit
     if ( index(dateForm(date), 'yyyy') == 0 ) then
       print *, 'Sorry--date format not found: ', trim(date)
       cycle
     endif
     n_dates = n_dates + 1
     dates(n_dates) = date
  enddo
  do i=1, n_dates
    date = dates(i)
    fromForm = trim(dateForm(date))
    ! Figure out logical format to convert it to
    if ( options%outputFormat /= ' ' ) then
      toForm = options%outputFormat
    elseif ( index(fromForm, 'doy') > 0 ) then
      toForm = MFORMAT
    else
      toForm = DOYFORMAT
    endif
    if ( options%verbose) then
      print *, 'Input was ', trim(date)
      print *, 'Input format was ', trim(fromForm)
      print *, 'Output format is ', trim(toForm)
    endif
    
    ! Process date
    converted_date = reFormatDate(date, &
      & fromForm=trim(fromForm), toForm=trim(toForm))

    print *, trim(converted_date)
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
    character(len=1), parameter :: y = 'y'
    character(len=1), parameter :: m = 'm'
    character(len=1), parameter :: d = 'd'
    ! Executable
    form = 'unknown format'
    if ( len_trim(date) < 1 ) return
    form = ' '
    s = 'y'
    i = 0
    j = 0
    do
      if ( i >= len_trim(date) ) exit
      i = i + 1
      j = j + 1
      select case (date(j:j))
      case ('d')
        form(i:i+2) = 'doy'
        i = i + 3
        j = j + 3
      case ('J', 'F', 'M', 'A', 'S', 'O', 'N', 'D')
        month = monthNameToNumber(date(i:))
        if ( month < 1 .or. month > 12 ) then
          form = 'month name uncrecognized in ' // trim(date(i:))
          return
        endif
        j = j + len_trim(MONTHNAME(month)) - 1
        form(i:i) = 'M'
        s = 'd'
        ! write(tempFormat(5:6),'(i2.2)') month
      case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        select case (s)
        case ('m')
          ! Was yyyy, now mm
          form(i:i+1) = 'mm'
          s = 'd'
          i = i + 1
          j = j + 1
        case ('d')
          ! Was mm, now dd
          form(i:i+1) = 'dd'
          s = ' '
          i = i + 1
          j = j + 1
        case ('y')
          ! yyyy
          form(i:i+3) = 'yyyy'
          s = 'm'
          i = i + 3
          j = j + 3
        case default
          ! Huh? Already finished with dd
          print *, 'Unexpected digit in dateForm'
        end select
      case default
        form(i:i) = date(j:j)
      end select
    enddo
  end function dateForm
  function monthNameToNumber(name) result(number)
    ! Convert month name to corresponding number
    ! E.g., given 'March', returns 3
    ! As a courtesy, name may be case-insensitive
    ! As a further courtesy, name may be followed by any junk you like
    ! Thus 'March 23, 2004 01:59:59.999' still returns 3
    ! If no such month is found, returns -1
    ! Args
    character(len=*), intent(in)             :: name
    integer                                  :: number
    do number=1, size(MONTHNAME)
      if ( index(lowerCase(name), lowercase(trim(MONTHNAME(number)))) > 0 ) return
    enddo
    number = -1
  end function monthNameToNumber

!------------------------- get_date ---------------------
    subroutine get_date(date, n_dates, options)
    ! Added for command-line processing
     character(LEN=255), intent(out) :: date          ! date
     integer, intent(in)             :: n_dates
     type ( options_T ), intent(inout) :: options
     ! character(LEN=*), intent(inout) :: outputFile        ! output date
     ! logical, intent(inout)          :: verbose
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, date )
      ! print *, i, ' th Arg: ', trim(date)
      error = 0
      if ( date(1:1) /= '-' ) exit
      if ( date(1:3) == '-h ' ) then
        call print_help
      elseif ( date(1:3) == '-o ' ) then
        call getarg ( i+1+hp, options%outputFormat )
        i = i + 1
        exit
      elseif ( date(1:3) == '-v ' ) then
        options%verbose = .true.
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
      print *,  "Enter the date you wish to convert."
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
      write (*,*) ' Options: -o format   => output format to use (e.g. yyyymmdd)'
      write (*,*) '          -v          => switch on verbose mode'
      write (*,*) '          -h          => print brief help'
      stop
  end subroutine print_help

!==================
END PROGRAM dateconverter
!==================

! $Log$
