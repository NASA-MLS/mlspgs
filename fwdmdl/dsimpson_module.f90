module DSIMPSON_MODULE
  use D_CSPLINE_M, only: CSPLINE
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public DSIMPS, SIMPS
  interface SIMPS; module procedure DSIMPS; end interface
  integer, parameter :: RK = r8
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Subroutine DSimps (F,DX,N,R)
!  Simpson's Integration of discrete equal spacing
    include 'simpson.f9h'
  End Subroutine DSimps
end module DSIMPSON_MODULE
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
