module TIME_MOD
!     .  Copyright (C) 1989-1999, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.

! TK is "time kind," the kind parameter for real numbers for times in
! seconds since YBASE (an argument to TEXT_TIME and TIME_TEXT).  14 digits
! is adequate for 316 years with millisecond resolution.

  implicit none
  public
  integer, parameter :: TK = selected_real_kind(14)

  character(len=*), private, parameter :: FILL = '1996-260T00:00:00.000000'

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

!================================================  DIFF_TEXT_TIME  =====
  real(kind=tk) function DIFF_TEXT_TIME ( TEXT1, TEXT2 )

! Return the difference in seconds from the time given by TEXT1 (see
! TEXT_TIME) to the time given by TEXT2.

    character(len=*), intent(in) :: TEXT1, TEXT2
    diff_text_time = text_time(text2,1992) - text_time(text1,1992)
    return
  end function DIFF_TEXT_TIME

!=================================================  DURATION_TEXT  =====
  subroutine DURATION_TEXT ( DURATION, TEXT )

! Given a duration, convert it to hh:mm:ss.ssssss
! The "hh" part will be "**" if the duration is 100 hours or more.

    real(tk), intent(in) :: DURATION
    character(len=*), intent(out) :: TEXT

    integer :: HH, L, MM, MMS, MS, SS

    hh = duration / 3600.0_tk
    mm = mod(int(duration/60.0_tk), 60)
    ss = mod(int(duration), 60)
    ms = 1000 * mod(duration, 1.0_tk)
    mms = 1000 * mod(duration*1000.0_tk, 1.0_tk)
    text = '00:00:00.000000'
    l = len(text)
    if ( l >= 2 ) then
      write ( text(1:2), '(i2.2)' ) hh             ! hours
      if ( l >= 5 ) then
        write ( text(4:5), '(i2.2)' ) mm           ! minutes
        if ( l >= 8 ) then
          write ( text(7:8), '(i2.2)' ) ss         ! seconds
          if ( l >= 12) then
            write ( text(10:12), '(i3.3)' ) ms     ! milliseconds
            if ( l >= 15) &
              write ( text(13:15), '(i3.3)' ) mms  ! microseconds
          end if
        end if
      end if
    end if
  end subroutine DURATION_TEXT
!=================================================  INC_TEXT_TIME  =====
  function INC_TEXT_TIME ( TEXT, INC ) result ( TEXT_INCD )

! Increment the time represented by TEXT (see TEXT_TIME) by INC
! seconds, and return the result as text.

    character(len=*), intent(in) :: TEXT
    real(kind=tk), intent(in) :: INC
    character(len=len(fill)) :: TEXT_INCD

    text_incd = time_text( text_time( text, 1992 ) + inc, 1992 )
    return
  end function INC_TEXT_TIME

!================================================  INC_YMDHS_TIME  =====
  subroutine INC_YMDHS_TIME ( YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, INC )

! Increment the time given by YEAR, MONTH, DAY, HOUR, MINUTE, SECOND by
! INC seconds.

    integer, intent(inout) :: YEAR, MONTH, DAY, HOUR, MINUTE
    real(kind=tk), intent(inout) :: SECOND
    real(kind=tk), intent(in) :: INC

    call time_ydhms ( ydhms_time(year, month, day, hour, second, year), &
                           year, year, month, day, hour, second )
    return
  end subroutine INC_YMDHS_TIME

!========================================================  JULDAY  =====
  integer function JULDAY (MONTH, DAY, YEAR)

! Return the Julian day for the date given by MONTH, DAY, YEAR.

    integer, intent(in) :: MONTH, DAY, YEAR

! Gregorian Calander was adopted on Oct. 15, 1582
    integer, parameter :: GREG = 15 + 31 * (10 + 12 * 1582)
    integer :: JA, JM, JY

    if (month .gt. 2) then
      jy = year
      jm = month + 1
    else
      jy = year - 1
      jm = month + 13
    end if

    julday = INT(365.25_TK * jy) + INT(30.6001_TK * jm) + &
                   day  + 1720995

! Test whether to change to Gregorian Calendar
    if ((day + 31 * (month + 12 * year)) .ge. greg) then
      ja = 0.01 * jy
      julday = julday + 2 - ja + 0.25 * ja
    endif

    return
  end function JULDAY

!=====================================================  NOW_CCSDS  =====

  subroutine NOW_CCSDS ( TEXT )

! Get the current time and put it in CCSDS format: yyyy-dddThh:mm:ss.ssssss
! The Fortran standard DATE_AND_TIME intrinsic has only millisecond
! resolution, so the last three digits will be zero.

    character(len=*), intent(out) :: TEXT

    integer :: L              ! LEN(TEXT)
    integer :: VALUES(8)      ! From DATE_AND_TIME intrinsic

    text = '0000-000T00:00:00.000000' ! yyyy-dddThh:mm:ss.ssssss

    call date_and_time ( values=values )
    l = len(text)
    if ( l >= 4 ) then
      write ( text(1:4), '(i4.4)' ) values(1)             ! yyyy
      if ( l >= 8 ) then
        write ( text(6:8), '(i3.3)' ) &                   ! ddd
                julday(values(2),values(3),values(1)) - julday(1,0,values(1))
        if ( l >= 11 ) then
          write ( text(10:11), '(i2.2)' ) values(5)       ! hh
          if ( l >= 14 ) then
            write ( text(13:14), '(i2.2)' ) values(6)     ! mm
            if ( l >= 17 ) then
              write ( text(16:17), '(i2.2)' ) values(7)   ! ss
              if ( l >= 21 ) &
                write ( text(19:21), '(i3.3)' ) values(8) ! .sss
            end if
          end if
        end if
      end if
    end if
    return
  end subroutine NOW_CCSDS

!=====================================================  TEXT_TIME  =====
  function TEXT_TIME (TEXT, YBASE, MBASE, DBASE) result ( TIME )

! Convert a time in TEXT given as yyyy-dddThh:mm:ss.ssssss to
! seconds since midnight of year YBASE, month MBASE, day DBASE.
! MBASE and DBASE are optional; their defaults are both 1.

! Returns -1.0 if TEXT isn't of the correct format.

    character(len=*), intent(in) :: TEXT
    integer, intent(in) :: YBASE
    integer, intent(in), optional :: MBASE, DBASE
    real(kind=tk) :: TIME

    character(len=len(fill)) :: MY_TEXT

    integer :: BASE, YEAR, DAY, HOUR, MINUTE, MY_DBASE, MY_MBASE
    real(kind=tk) :: SECOND

    my_dbase = 1
    if (present(dbase)) my_dbase = dbase
    my_mbase = 1
    if (present(mbase)) my_mbase = mbase
    base = len_trim(text) + 1
    my_text(:base-1) = text(:base-1)
    my_text(base:) = fill(base:)
    if (my_text(5:5)   .ne. '-' .or. my_text(9:9)   .ne. 'T' &
   .or.my_text(12:12) .ne. ':' .or. my_text(15:15) .ne. ':' &
   .or.my_text(18:18) .ne. '.') go to 200

    read (my_text,100,err=200) year, day, hour, minute, second
100 format (i4,1x,i3,1x,2(i2,1x),f9.6)

    base = julday (my_mbase, my_dbase, ybase)
    time = real(julday(1,1,year) - base + day - 1, tk) * 86400.0_TK + &
                second + 60.0_TK * (minute + 60.0_TK * hour)
    return

200 time = -1.0_TK
    return

  end function TEXT_TIME

!=====================================================  TIME_TEXT  =====
  character*(21) function TIME_TEXT ( SINCE, YBASE )

! Return a text string of the form YYYY-DDDTHH:MM:SS.SSS giving the
! date and time corresponding to time in seconds SINCE 0:00 1 January
! YBASE.

    real(KIND=TK), intent(in) :: SINCE  ! Time in seconds after YBASE
    integer, intent(in) :: YBASE

    real(KIND=TK) :: SECOND
    integer :: YEAR, DAY, HOUR, MINUTE, SEC, MS

    call time_ydhms ( since, ybase, year, day, hour, minute, second )
    sec = second
    ms = 1000 * (second - sec)
    write (time_text, 200) year, day, hour, minute, sec, ms
200 format (i4.4,'-',i3.3,'T',2(i2.2,':'),i2.2,'.',i3.3)
    return
  end function TIME_TEXT

!====================================================  TIME_YDMHS  =====

  subroutine TIME_YDHMS ( SINCE, YBASE, YEAR, DAY, HOUR, MINUTE, SECOND )

! Return year, day, hour, minute, second corresponding to time in
! seconds SINCE 0:00 January 1 YBASE.

    real(KIND=TK), intent(in) :: SINCE  ! Time in seconds after YBASE
    integer, intent(in) :: YBASE
    integer, intent(out) :: YEAR, DAY, HOUR, MINUTE
    real(KIND=TK), intent(out) :: SECOND

    real(KIND=TK) :: DATE, DHOUR, DMIN
    integer :: LEAP

    date = since / 86400.0_tk + 1
    year = ybase
    do
      leap = 0
      if (mod(year,4) == 0 &
       .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) leap = 1
      if (date <= 365 + leap) exit
      year = year + 1
      date = date - 365 - leap
    end do
    day = date
    dhour = 24.0 * (date - day)
    hour = dhour
    dmin = (24.0*(date - day) - hour) * 60.0
    minute = dmin
    second = (dmin - minute) * 60.0
    return
  end subroutine TIME_YDHMS

!====================================================  YDHMS_TIME  =====
  function YDHMS_TIME (YEAR, DAY, HOUR, MINUTE, SECOND, YBASE, MBASE, DBASE) &
  result ( TIME )

! Convert YEAR, DAY, HOUR, MINUTE, SECOND to seconds since 0:00 of year
! YBASE, month MBASE, day DBASE.
! HOUR, MINUTE and SECOND are optional; their defaults are zero.
! MBASE and DBASE are optional; their defaults are both 1.

    integer, intent(in) :: YEAR, DAY
    integer, intent(in), optional :: HOUR, MINUTE
    real(kind=tk), intent(in), optional :: SECOND
    integer, intent(in) :: YBASE
    integer, intent(in), optional :: MBASE, DBASE
    real(kind=tk) :: TIME

    integer :: BASE, MY_HOUR, MY_MINUTE, MY_DBASE, MY_MBASE
    real(kind=tk) :: MY_SECOND

    my_hour = 0
    if (present(hour)) my_hour = hour
    my_minute = 0
    if (present(minute)) my_minute = minute
    my_second = 0.0_tk
    if (present(second)) my_second = second
    my_dbase = 1
    if (present(dbase)) my_dbase = dbase
    my_mbase = 1
    if (present(mbase)) my_mbase = mbase

    base = julday (my_mbase, my_dbase, ybase)
    time = real(julday(1,1,year) - base + day - 1, tk) * 86400.0_TK + &
                my_second + 60.0_TK * (my_minute + 60.0_TK * my_hour)
    return
  end function YDHMS_TIME

end module TIME_MOD

! $Log$
! Revision 1.3  2000/05/04 23:45:37  vsnyder
! Initial entry into CVS
!
