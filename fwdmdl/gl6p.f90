module GL6P
! **************  Gauss-Legendre 6 point formula ***************
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: NG, GX, GW
  Integer(i4), parameter :: Ng = 6
!
! These are the 6-point-Gauss-Legendre abscissa (X-axis) values in [-1,1]:
!
  Real(r8), parameter :: Gx(Ng) = (/ &
 &  -9.32469514203152027812d-1, -6.61209386466264513661d-1,   &
 &  -2.38619186083196908631d-1,  2.38619186083196908631d-1,   &
 &   6.61209386466264513661d-1,  9.32469514203152027812d-1 /)
!
! These are the corresponding 6-point-Gauss-Legendre Weights values:
!
  Real(r8), parameter :: Gw(Ng) = (/ &
 &   1.71324492379170345040d-1,  3.60761573048138607569d-1,   &
 &   4.67913934572691047390d-1,  4.67913934572691047390d-1,   &
 &   3.60761573048138607569d-1,  1.71324492379170345040d-1 /)
!
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
end module GL6P
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
