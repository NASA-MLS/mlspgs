! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GLNP

! **************  Gauss-Legendre 3 point formula ***************
  use MLSCommon, only: R8
  implicit NONE
  private
  public :: NG, GX, GW
  integer, parameter :: Ng = 3
!
! These are the 3-point-Gauss-Legendre abscissae (X-axis) values in [-1,1]:
!
  real(r8), parameter :: Gx(Ng) = (/ & ! sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)
     & -7.74596669241483377036e-1_r8,  0.00000000000000000000e+0_r8, &
     &  7.74596669241483377036e-1_r8 /)
!
! These are the corresponding 3-point-Gauss-Legendre Weights values:
!
  real(r8), parameter :: Gw(Ng) = (/ & ! 5.0/9.0, 8.0/9.0, 5.0/9.0
     &  5.55555555555555555556e-1_r8,  8.88888888888888888889e-1_r8, &
     &  5.55555555555555555556e-1_r8 /)
!
!---------------------------- RCS Ident Info -------------------------------
!
  character (len=*), parameter :: IdParm = &
    & "$id: glnp.f90,v 2.0 2001/09/17 20:26:27 livesey Exp $"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

end module GLNP

! $Log$
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1  2001/06/21 13:07:09  zvi
! Speed enhancement MAJOR update
!
! Revision 1.1  2001/06/06 00:00:00  Z.Shippony
! Initial conversion to Fortran 90
