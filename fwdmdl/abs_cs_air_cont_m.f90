! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ABS_CS_AIR_CONT_M
  implicit NONE
  private
  public :: ABS_CS_AIR_CONT
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
! This function computes the water continuum contribution to the dry air

! pure real(rk) Function ABS_CS_AIR_CONT ( temperature, pressure, frequency )
  real(rk) Function ABS_CS_AIR_CONT ( temperature, pressure, frequency )

    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
    real(rk), intent(in) :: PRESSURE    ! in mbar
    real(rk), intent(in) :: FREQUENCY   ! in MegaHertz
!
! Set numbers to the LAb numbers (fits for EOS for now ..)
!
    real(rk), Parameter :: BETA=1.07E-19_RK, FPWR=-1.85E-12_RK, &
                           PPWR=2.0_RK, TPWR=3.63_RK
!
    real(rk) :: THETA, FSQR
!
    theta = 300.0_rk / temperature
    fsqr = frequency * frequency
    abs_cs_air_cont = beta * fsqr * (pressure ** ppwr) *    &
                      (theta ** tpwr) * dexp(fpwr * fsqr)
!
  End function abs_cs_air_cont
end module ABS_CS_AIR_CONT_M
! $Log$
! Revision 1.5  2001/04/20 23:29:45  zvi
! Setting the correct (Lab) numbers for N2, fits for EOS for now ..
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
