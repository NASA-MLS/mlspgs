! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CCSDS_Time

  ! Several procedures to convert between CCSDS time representations
  ! (yyyy-doyThh:mm:ss.sss or yyyy-mm-ddThh:mm:ss.sss) and the Time_T
  ! structure from the Calendar module.

  implicit NONE
  private
  public :: CCSDS_TO_TIME, CHECK_CCSDS_TIME, TIME_TO_CCSDS

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------  CCSDS_TO_TIME  -----
  type(time_t) function CCSDS_TO_TIME ( CCSDS ) result ( T )
    ! If t%day == 0, there's an error code from Check_CCSDS_Time in t%year
    use Calendar, only: Operator(+), Time_T
    character(len=*), intent(in) :: CCSDS
    character(len=len_trim(ccsds)+1) :: C
    integer :: STAT
    c = ccsds
    stat = check_ccsds_time(c)
    select case ( stat )
    case ( :-1 ) ! yyyy-dddThh:mm:ss.sss format
      read ( c, * ) t%year, t%day, t%hours, t%minutes, t%seconds
      t = time_t(t%year, 1, 0, t%hours, t%minutes, t%seconds) + t%day
    case ( 0 )   ! yyyy-mm-ddThh:mm:ss.sss format
      read ( c, * ) t%year, t%month, t%day, t%hours, t%minutes, t%seconds
    case ( 1: )  ! An error, return zeros
      t = time_t(stat,0,0,0,0,0)
    end select
    if ( t%year == 0 ) t%year = -1 ! No year zero, assume 1 BC.
  end function CCSDS_TO_TIME

  ! -------------------------------------------  CHECK_CCSDS_TIME  -----
  integer function CHECK_CCSDS_TIME ( CCSDS )

    ! Returns -1 if CCSDS is yyyy-dddThh:mm:ss.sss[Z], 0 if CCSDS is
    ! yyyy-mm-ddThh:mm:ss.sss[Z], else either the position of the error or
    ! a processor's I/O error code.  ss.sss might have more digits, or
    ! no decimal point.  Checks values of hh, mm (both), ss.  Checks for
    ! dd[d] <= 0.  Puts commas into CCSDS in place of -, T, :.
    ! CCSDS allows either the date or time and the T separator to be omitted.
    ! We allow only the time to be omitted.  The :ss.sss or :mm:ss.sss can
    ! also be omitted, with zeros assumed in those cases.  If the string
    ! does not end with a non-digit, the last digit will be lost (covered
    ! up by a slash).
    ! This is a little bit lenient in that it does not require exactly two
    ! digits for hh, mm, ss, dd, four for yyyy, or three for ddd.  It also
    ! allows negative years (for BC dates).

    character(len=*), intent(inout) :: CCSDS
    character(len=*), parameter :: DIGITS = '0123456789.T-:Z'
    logical :: DOY
    integer :: H1, M2, SEP1, SEP2, Y2 ! locations of separators
    integer :: L
    character(len=len_trim(ccsds)+1) :: T
    integer :: YY, MO, DD, HH, MM
    real :: SS

    l = len_trim(ccsds)
    check_ccsds_time = verify(ccsds(:l),digits) ! catch most junk, including blanks
    if ( check_ccsds_time /= 0 ) return
    check_ccsds_time = l
    h1 = scan(ccsds(:l),'Tt')
    if ( h1 /= 0 ) then ! Got the T
      ccsds(h1:h1) = ','
      check_ccsds_time = h1
    else
      h1 = l + 1
    end if
    y2 = scan(ccsds(2:h1-1),'-') + 1 ! Allow leading hyphen -- for BC years
    if ( y2 == 1 ) return ! No hyphen-after-year at all
    ccsds(y2:y2) = ','
    doy = .true.
    m2 = scan(ccsds(y2+1:h1-1),'-') + y2
    if ( m2 /= y2 ) then ! Two hyphens
      doy = .false.
      ccsds(m2:m2) = ','
    end if
    if ( h1 /= 0 ) then ! Got the T
      if ( scan(ccsds(l:l),'Zz') /= 0 ) l = l - 1 ! Final Z is optional
      sep1 = scan(ccsds(h1+1:l),':') + h1
      check_ccsds_time = l
      if ( sep1 /= h1 ) then ! Got the first colon
        ccsds(sep1:sep1) = ','
        sep2 = scan(ccsds(h1+1:l),':') + h1
        check_ccsds_time = l
        if ( sep2 /= h1 ) ccsds(sep2:sep2) = ','
        ss = 0
      else
        mm = 0
      end if
    else
      hh = 0
    end if
    l = min(l+1,len(ccsds))
    ccsds(l:l) = '/' ! So that Fortran I/O will skip missing fields.
    do l = l-1, 1, -1
      if ( ccsds(l:l) /= ',' ) exit
      ccsds(l:l) = '/' ! Replace trailing commas by /
    end do

    if ( doy ) then
      read ( ccsds, *, iostat=l ) yy, dd, hh, mm, ss
    else
      read ( ccsds, *, iostat=l ) yy, mo, dd, hh, mm, ss
      check_ccsds_time = y2 + 1
      if ( mo <= 0 .or. mo > 12 ) return
    end if
    check_ccsds_time = h1 - 1
    if ( dd <= 0 ) return
    check_ccsds_time = h1 + 1
    if ( hh < 0 .or. hh >= 24 ) return
    check_ccsds_time = sep1 + 1
    if ( mm < 0 .or. mm >= 60 ) return
    check_ccsds_time = sep2 + 1
    if ( ss < 0.0d0 .or. ss >= 60.0 ) return
    check_ccsds_time = 0
    if ( doy ) check_ccsds_time = -1

  end function CHECK_CCSDS_TIME

  ! ----------------------------------------------  TIME_TO_CCSDS  -----
  subroutine TIME_TO_CCSDS ( TIME, CCSDS, DOY, SECDIG )
  ! Convert TIME to CCSDS.  yyyy-dddThh:mm:ss if DOY, else yyyy-mm-ddThh:mm:ss
  ! SECDIG specifies the number of digits after the decimal for the seconds.
  ! <0 => no decimal (integer format).  Default -1.
  ! It's your responsibility to make sure CCSDS is long enough!
    use Calendar, only: Operator(-), Time_T
    type(time_t), intent(in) :: TIME
    character(len=*), intent(out) :: CCSDS
    logical, intent(in), optional :: DOY ! Default .false.
    integer, intent(in), optional :: SECDIG
    integer :: Dig, L
    logical :: MyDOY
    character(10) :: FMT
    myDoy = .false.
    if ( present(doy) ) myDoy = doy
    dig = -1
    if ( present(secdig) ) dig = secdig
    if ( myDoy ) then
      write ( ccsds, '(i4.4,"-",i3.3,"T",2(i2.2,":"))' ) time%year, &
        & nint(time_t(time%year,time%month,time%day,0,0,0)-time_t(time%year,1,0,0,0,0)), &
        & time%hours, time%minutes
      l = 16
    else
      write ( ccsds, '(i4.4,2("-",i2.2),"T",2(i2.2,":"))' ) time%year, &
        & time%month, time%day, time%hours, time%minutes
      l = 18
    end if
    if ( dig < 0 ) then
      write ( ccsds(l:), '(i2)' ) nint(time%seconds)
    else
      write ( fmt, '("(f",i3,".",i3,")")' ) dig+3,dig
      write ( ccsds(l:), fmt ) time%seconds
    end if
    if ( ccsds(l:l) == ' ' ) ccsds(l:l) = '0'
  end subroutine TIME_TO_CCSDS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CCSDS_Time

! $Log$
! Revision 2.1  2004/10/18 20:49:53  vsnyder
! Initial commit
!
