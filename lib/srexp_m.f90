! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SREXP_M

! ----------------------------------------------------------------------
!            EVALUATE THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------

  implicit NONE
  private
  public :: REXP, SREXP

  interface REXP; module procedure SREXP; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

  function SREXP ( X ) result ( REXP )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X
    real(rk) :: REXP
    include 'rexp.f9h'
  end function SREXP

end module SREXP_M

! $Log$
! Revision 2.3  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.2.2.1  2003/04/16 20:12:33  vsnyder
! Move argument declarations from rexp.f9h to here
!
! Revision 2.2  2002/10/10 20:36:50  vsnyder
! Darn, I forgot the CVS log comment
!
