module D_CSPLINE_M
  use D_HUNT_M, only: HUNT
  use D_PCSPL_M, only: PCSPL
  implicit NONE
  private
  public D_CSPLINE, CSPLINE
  interface CSPLINE; module procedure D_CSPLINE; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
!
  subroutine D_CSPLINE (XIN, XOUT, YIN, YOUT, NIN, NOUT)
    include 'cspline.f9h'
  end subroutine D_CSPLINE
end module D_CSPLINE_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
