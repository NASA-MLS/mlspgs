module ABS_CS_N2_CONT_M
  implicit NONE
  private
  public :: ABS_CS_N2_CONT
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)

contains

! This function computes the N2 continuum contribution

  real(rk) Function ABS_CS_N2_CONT (cont, temperature, pressure, frequency )

    real(rk), intent(in) :: CONT(:)     ! continuum parameters
    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
    real(rk), intent(in) :: PRESSURE    ! in mbar
    real(rk), intent(in) :: FREQUENCY   ! in MegaHertz
!
    real(rk) :: THETA, FSQR, FSXT
!
    theta = 300.0_rk / temperature
    fsqr = frequency * frequency
    fsxt = fsqr * theta
    abs_cs_n2_cont = pressure * pressure * fsqr * (theta**cont(2)) * &
                   & ( cont(1) * exp(cont(3) * fsxt) +               &
                   &   cont(4) * exp(cont(5) * fsxt) * (cont(6)**2 + fsqr))
!
  End function abs_cs_n2_cont
end module ABS_CS_N2_CONT_M
! $Log$
! Revision 2.2  2001/10/17 17:03:17  zvi
! Fix a bug
!
! Revision 2.1  2001/10/16 14:45:02  zvi
! New Dry air function (N2) with continuum parameters passed in
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
