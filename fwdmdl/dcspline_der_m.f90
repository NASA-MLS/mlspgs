module DCSPLINE_DER_M
  use D_HUNT_M, only: HUNT
  use D_PCSPL_M, only: PCSPL
  implicit NONE
  private
  public :: DCSPLINE_DER, CSPLINE_DER
  interface CSPLINE_DER; module procedure DCSPLINE_DER; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, parameter :: RK = kind(0.0d0)
!
contains
  subroutine DCSPLINE_DER ( XIN, XOUT, YIN, YOUT, DYOUT, NIN, NOUT )
    include 'cspline_der.f9h'
  end subroutine DCSPLINE_DER
end module DCSPLINE_DER_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
