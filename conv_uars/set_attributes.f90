! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Set_Attributes_m

  implicit NONE
  private

  public :: Set_Attributes

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public procedures ===================================

  subroutine Set_Attributes ( sdId, YrDoy, Start_Time, End_Time )

    use Calendar, only: Calend
    use PCFHdr, ONLY: h5_writeglobalattr, GlobalAttributes, FillTAI93Attribute
    use Time_m, only: MS_to_HMS

    integer, intent (in) :: sdId, YrDoy, Start_Time(2), End_Time(2)

    integer :: day, doy, month, year, yr_range(2), doy_range(2)
    integer :: hrs, mins, secs

    character (len=17) :: start_utc, end_utc
    character (len=*), parameter :: &
         DateFmt = "(I4, '-', I3.3, 'T', I2.2, ':', I2.2, ':', I2.2)"
    character (len=*), parameter :: inst_name = 'MLS UARS'

    ! year and doy of file:

    year = yrdoy / 1000 + 1900
    doy = mod ( yrdoy, 1000 )

    ! convert DOY to Month, Day:

    call calend ( year, doy, month, day )
    print '(a,4(1x,i0))', 'year, doy, month, day: ', year, doy, month, day

    ! convert start millisecs to HMS:

    call ms_to_hms ( start_time(2), hrs, mins, secs )

    yr_range(1) = start_time(1) / 1000 + 1900
    doy_range(1) = mod (start_time(1), 1000)
    write ( start_utc, fmt=DateFmt ) yr_range(1), doy_range(1), hrs, mins, secs

    ! convert stop millisecs to HMS:

    call ms_to_hms (end_time(2), hrs, mins, secs)

    yr_range(2) = end_time(1) / 1000 + 1900
    doy_range(2) = mod (end_time(1), 1000)
    write ( end_utc, fmt=DateFmt ) yr_range(2), doy_range(2), hrs, mins, secs

    ! Store global attributes

    GlobalAttributes%StartUTC = Start_UTC
    GlobalAttributes%EndUTC = End_UTC
    GlobalAttributes%ProcessLevel = '1'
    GlobalAttributes%InstrumentName = 'MLS UARS'
    GlobalAttributes%PGEVersion = 'Test'
    GlobalAttributes%GranuleYear =  year
    GlobalAttributes%GranuleMonth = month
    GlobalAttributes%GranuleDay = day

    call FillTAI93Attribute

    ! Write global attributes:

    call h5_writeglobalattr ( sdId, skip_if_already_there=.false. )

  end subroutine Set_Attributes

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Set_Attributes_m

! $Log$
! Revision 1.2  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
