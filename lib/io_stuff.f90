! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module IO_STUFF

! Useful stuff for I/O

  implicit NONE

  private
  public GET_LUN

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

! ================================================     GET_LUN     =====

  subroutine GET_LUN ( LUN )
  ! Find a Fortran logical unit number that's not in use.
    integer, intent(out) :: LUN    ! The logical unit number
    logical :: OPENED              ! Used to inquire about the unit
    do lun = 1, 100
      inquire ( unit=lun, opened=opened )
      if (.not. opened) return
    end do
    write(*,*) 'IO_STUFF%GET_LUN-E- Unable to get a logical unit number'
    lun = -1
    return
  end subroutine GET_LUN
end module IO_STUFF

! $Log$
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
