module D_GET_ONE_ETA_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: D_GET_ONE_ETA, GET_ONE_ETA
  interface GET_ONE_ETA; module procedure D_GET_ONE_ETA; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = r8

contains

! This subroutine gets the eta function for temperature for one height & one
! coefficient:
!
  Subroutine D_GET_ONE_ETA ( h, peaks, no_peaks, iq, eta )
    include 'get_one_eta.f9h'
  End Subroutine D_GET_ONE_ETA
end module D_GET_ONE_ETA_M
