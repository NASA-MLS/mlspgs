module AITKEN_INT_M
  use D_CSPLINE_M, only: CSPLINE
  use DAITKEN_MODULE, only: AITKEN
  use DSIMPSON_MODULE, only: SIMPS
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: Aitken_int
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Function Aitken_int (F_GRID, F_GRID_FLTR, FLTR_FUNC, RAD, JP, NFP, IER)
    real(r8) :: Aitken_int
    real(r8), intent(in) :: F_GRID(:)
    real(r8) :: F_GRID_FLTR(*)
    real(r8), intent(in) :: FLTR_FUNC(:)
    real(r8), intent(in) :: RAD(:)
    integer(i4), intent(in) :: JP
    integer(i4), intent(in) :: NFP
    integer(i4), intent(out) :: IER
    Integer(i4) :: N, KZ, JUMP
    Integer(i4), Parameter :: M = 161       ! For 161 maxfiltpts
    Real(r8) :: A(3)
    Real(r8) :: RXF(m), TMPARY(m), TMPF(m), DF
    ier = 2000
    if ( nfp > m ) then
      Print *, '** Error in subroutine Aitken_Int !'
      Print *, '   Number of filter points too large:', nfp
      Print *, '   Maximum allowed:', m
      Return
    endif
    ier = 0
    jump = 1
    df = f_grid_fltr(2) - f_grid_fltr(1)
    do kz = 1, 3
      n = jp/jump
      tmpf(1:n) = f_grid(1:jp:jump)
      rxf(1:n) = rad(1:jp:jump)
      Call Cspline ( tmpf, f_grid_fltr, rxf, tmpary, n, nfp )
      rxf = tmpary * fltr_func
      Call Simps ( rxf, df, nfp, a(kz) )
      jump = 2 * jump
    end do
    Aitken_int = Aitken(a(3),a(2),a(1))
    Return
  End Function Aitken_int
end module AITKEN_INT_M
! $Log$
! Revision 1.1  2000/06/21 21:56:09  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
