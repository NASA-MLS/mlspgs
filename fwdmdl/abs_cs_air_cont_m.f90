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

  pure real(rk) Function ABS_CS_AIR_CONT ( temperature, pressure, frequency )
    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
    real(rk), intent(in) :: PRESSURE    ! in mbar
    real(rk), intent(in) :: FREQUENCY   ! in MegaHertz
!
    real(rk), Parameter :: BETA=1.90E-19_RK, FPWR=-1.85E-12_RK, &
                           PPWR=2.0_RK, TPWR=2.79_RK
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
