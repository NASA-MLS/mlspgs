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

! pure real(rk) Function ABS_CS_O2_CONT (p, t, frq )
  real(rk) Function ABS_CS_O2_CONT (p, t, frq )

    real(rk), intent(in) :: P    ! pressure, mbar
    real(rk), intent(in) :: T    ! temperature, Kelvins
    real(rk), intent(in) :: FRQ  ! frequency, MegaHertz
!
    real(rk), Parameter :: DEBYE=1.10816E-3_RK, W_NONRES=0.48_RK, &
                           NO2=0.8_RK
!
    real(rk) :: F2, QT, SLABS_NONRES, Y_NONRES, Y2
!
    f2 = frq * frq
    qt = (300.0_rk / t) ** no2
    slabs_nonres = debye * p / (t * t)
    y_nonres = p * w_nonres * qt
    y2 = y_nonres * y_nonres
    abs_cs_o2_cont = slabs_nonres * f2 * y_nonres / (f2 + y2)
!
  End Function ABS_CS_O2_CONT
end module ABS_CS_O2_CONT_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
