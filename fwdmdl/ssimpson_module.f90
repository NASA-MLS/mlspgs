module SSIMPSON_MODULE
  use MLSCommon, only: I4, R4
  use S_CSPLINE_M, only: CSPLINE
  implicit NONE
  private
  public SSIMPS, SIMPS
  interface SIMPS; module procedure SSIMPS; end interface
  integer, parameter :: RK = r4

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains

  Subroutine SSIMps ( F, DX, N, R )

!  Simpson's Integration of discrete equal spacing

    include 'simpson.inc'

  End Subroutine SSIMps
end module SSIMPSON_MODULE

! $Log$
