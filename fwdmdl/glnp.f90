! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GLNP
! module GL3P
! **************  Gauss-Legendre 3 point formula ***************
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: NG, GX, GW
  Integer(i4), parameter :: Ng = 3
!
! These are the 3-point-Gauss-Legendre abscissa (X-axis) values in [-1,1]:
!
  Real(r8), parameter :: Gx(Ng) = (/ &
     & -7.74596669241483377036d-1,  0.00000000000000000000d+0, &
     &  7.74596669241483377036d-1 /)
!
! These are the corresponding 3-point-Gauss-Legendre Weights values:
!
  Real(r8), parameter :: Gw(Ng) = (/ &
     &  5.55555555555555555556d-1,  8.88888888888888888889d-1, &
     &  5.55555555555555555556d-1 /)
!
!---------------------------- RCS Ident Info -------------------------------
!
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
! end module GL3P
end module GLNP
! $Log$
! Revision 1.1  2001/06/06 00:00:00  Z.Shippony
! Initial conversion to Fortran 90
