module FREQ_AVG_M
  use D_CSPLINE_M, only: CSPLINE
  use DSIMPSON_MODULE, only: SIMPS
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: Freq_Avg
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!
  Subroutine Freq_Avg(F_grid,F_grid_fltr,Fltr_func,Rad,N,Nfp,Avg,Ier)
!
    Real(r8), intent(in) :: FLTR_FUNC(:)
    Real(r8), intent(in) :: F_GRID(*), RAD(*), F_GRID_FLTR(*)

    Integer(i4), intent(in) :: N, NFP
!
    Integer(i4), intent(out) :: IER
    Real(r8), intent(out)    :: Avg
!
    Integer(i4), Parameter :: M = 161       ! For 161 maxfiltpts
!
    Real(r8) :: RXF(m), TMPARY(m), DF
!
    ier = 2000
    if ( nfp > m ) then
      Print *, '** Error in subroutine Freq_Avg routine !'
      Print *, '   Number of filter points too large:', nfp
      Print *, '   Maximum allowed:', m
      Return
    endif

    ier = 0
    df = f_grid_fltr(2) - f_grid_fltr(1)
    Call Cspline (f_grid, f_grid_fltr, rad, tmpary, n, nfp)
    rxf = tmpary * fltr_func
    Call Simps (rxf, df, nfp, Avg)

    Return

  End Subroutine Freq_Avg

end module FREQ_AVG_M
! $Log$
! Initial conversion to Fortran 90
! 2000/11/20 21:56:09  zvi
! First version D.P.
