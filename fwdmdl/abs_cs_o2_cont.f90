module ABS_CS_O2_CONT_M
  implicit NONE
  private
  public :: ABS_CS_O2_CONT
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)

contains

! This function computes the O2 continuum contribution

  real(rk) Function ABS_CS_O2_CONT (cont, temperature, pressure, frequency )

    real(rk), intent(in) :: CONT(:)     ! continuum parameters
    real(rk), intent(in) :: TEMPERATURE ! in Kelvin
    real(rk), intent(in) :: PRESSURE    ! in mbar
    real(rk), intent(in) :: FREQUENCY   ! in MegaHertz
!
    real(rk) :: THETA, FSQR
!
    theta = 300.0_rk / temperature
    fsqr = frequency * frequency
    abs_cs_o2_cont = cont(1) * pressure * pressure * fsqr * (theta**cont(2)) &
                   & / (fsqr + (cont(3) * pressure * (theta**cont(4)) )**2 )
!
  End function abs_cs_o2_cont
end module ABS_CS_O2_CONT_M
! $Log$
! Revision 2.1  2001/10/16 22:24:25  zvi
! O2 contiuum function
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
