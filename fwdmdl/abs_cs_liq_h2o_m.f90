! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ABS_CS_LIQ_H2O_M
  implicit NONE
  private
  public :: ABS_CS_LIQ_H2O
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
! This function computes the liquid water correction
!
! pure real(rk) function ABS_CS_LIQ_H2O ( frequency, temperature )
  real(rk) function ABS_CS_LIQ_H2O ( frequency, temperature )

    real(rk), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
!
! This function when multiplied by mass density (gm/m^3) of liquid droplet
! water gives absorption in Km^-1. Function comes from Liebe 1985 radio
! science paper and others.
!
    real(rk) :: TAU, EPSILON, THETA
!
! Begin calculation
!
    theta = 300.0_rk / temperature
    tau = 4.17e-8_rk * frequency * theta * exp(7.13_rk * theta)
    epsilon = (185.0_rk - 113.0_rk/theta) / (1.0_rk + tau * tau)
    abs_cs_liq_h2o = 1.886e-4_rk * frequency * tau * epsilon &
                   / ((6.9_rk + epsilon)**2 + (tau*epsilon)**2)
!
  End function ABS_CS_LIQ_H2O
end module ABS_CS_LIQ_H2O_M
! $Log$
! Revision 1.5  2001/06/07 23:30:33  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
