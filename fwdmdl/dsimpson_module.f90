! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
