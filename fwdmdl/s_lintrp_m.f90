module S_LINTRP_M
  use S_HUNT_M, only: HUNT
  implicit NONE
  private
  public S_LINTRP, LINTRP
  interface LINTRP; module procedure S_LINTRP; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = kind(0.0e0)

contains

  subroutine S_LINTRP ( XIN, XOUT, YIN, YOUT, NIN, NOUT )
    include 'lintrp.f9h'
  end subroutine S_LINTRP
end module S_LINTRP_M

! $Log$
