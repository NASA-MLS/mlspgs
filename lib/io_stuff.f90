! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module IO_STUFF

! Useful stuff for I/O

  implicit NONE

  private
  public GET_LUN

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ================================================     GET_LUN     =====

  subroutine GET_LUN ( LUN, MSG )
  ! Find a Fortran logical unit number that's not in use.
    integer, intent(out) :: LUN          ! The logical unit number
    logical, intent(in), optional :: MSG ! Print failure message if absent or .true.
    logical :: EXIST, OPENED             ! Used to inquire about the unit
    do lun = 20, 100
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) return
    end do
    lun = -1
    if ( present(msg) ) then
      if ( .not. msg ) return
    end if
    write(*,*) 'IO_STUFF%GET_LUN-E- Unable to get a logical unit number'
    return
  end subroutine GET_LUN

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module IO_STUFF

! $Log$
! Revision 2.4  2004/05/19 23:00:18  vsnyder
! Add optional MSG argument
!
! Revision 2.3  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2001/04/26 02:39:11  vsnyder
! Fix up CVS stuff
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
