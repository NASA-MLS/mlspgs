module S_CSPLINE_M
  use S_HUNT_M, only: HUNT
  use S_PCSPL_M, only: PCSPL
  implicit NONE
  private
  public S_CSPLINE, CSPLINE
  interface CSPLINE; module procedure S_CSPLINE; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = kind(0.0e0)

contains
!
  subroutine S_CSPLINE ( XIN, XOUT, YIN, YOUT, NIN, NOUT )
    include 'cspline.f9h'
  end subroutine S_CSPLINE
end module S_CSPLINE_M

! $Log$
