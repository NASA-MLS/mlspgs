module S_HUNT_M
  use MLSCommon, only: I4, R4
  implicit NONE
  private
  public :: S_HUNT, HUNT
  interface HUNT; module procedure S_HUNT; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = r4

contains

! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.

  subroutine S_HUNT ( ELEMENT, ARRAY, N, JLO, JHI )
    include 'hunt.f9h'
  end subroutine S_HUNT
end module S_HUNT_M

! $Log$
