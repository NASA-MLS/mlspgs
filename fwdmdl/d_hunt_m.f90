module D_HUNT_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: D_HUNT, HUNT
  interface HUNT; module procedure D_HUNT; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = r8

contains

! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.

  subroutine D_HUNT ( ELEMENT, ARRAY, N, JLO, JHI )
    include 'hunt.f9h'
  end subroutine D_HUNT
end module D_HUNT_M

! $Log$
