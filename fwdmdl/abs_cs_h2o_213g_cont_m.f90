module ABS_CS_H2O_213G_CONT_M
  implicit NONE
  private
  public :: ABS_CS_H2O_213G_CONT

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = kind(0.0d0)

contains

! This function computes the water continuum contribution to the h2o
! cross section.
!
  pure real(rk) Function ABS_CS_H2O_213G_CONT ( temperature, pressure, h2o_mr )
    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
    real(rk), intent(in) :: PRESSURE    ! in mbar
    real(rk), intent(in) :: H2O_MR      ! approximate mixing ratio (unitless)
!
    real(rk), Parameter :: BETA=5.68E-5_RK, HCOEF=18.1_RK, &
                           TPWR1=3.67_RK, TPWR2=1.98_RK
!
    real(rk) :: THETA
!
    theta = 300.0_rk / temperature
    abs_cs_h2o_213g_cont = beta * pressure * pressure *  &
           (theta**tpwr1) * (1.0_rk + hcoef * h2o_mr * (theta**tpwr2))
!
  End Function ABS_CS_H2O_213G_CONT
end module ABS_CS_H2O_213G_CONT_M

! $Log$
