module SCSPLINE_DER_M
  use S_HUNT_M, only: HUNT
  use S_PCSPL_M, only: PCSPL
  implicit NONE
  private
  public :: SCSPLINE_DER, CSPLINE_DER
  interface CSPLINE_DER; module procedure SCSPLINE_DER; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, parameter :: RK = kind(0.0e0)
!
contains

  subroutine SCSPLINE_DER ( XIN, XOUT, YIN, YOUT, DYOUT, NIN, NOUT )
    include 'cspline_der.f9h'
  end subroutine SCSPLINE_DER
end module SCSPLINE_DER_M

! $Log$
