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
  Subroutine Freq_Avg(F_grid,F_grid_fltr,Fltr_func,Rad,n,nfp,Avg)
!
    Real(r8), intent(in) :: Fltr_func(:)
    Real(r8), intent(in) :: F_grid(:), Rad(:), F_grid_fltr(:)

    Integer(i4), intent(in) :: n, nfp
!
    Real(r8), intent(out)    :: Avg
!
    Real(r8) :: rxf(nfp), tmpary(nfp), df
!
    df = abs(F_grid_fltr(2) - F_grid_fltr(1))
    Call Cspline(F_grid, F_grid_fltr, Rad, tmpary, n, nfp)
    rxf(1:nfp) = tmpary(1:nfp) * Fltr_func(1:nfp)
    Call Simps (rxf, df, nfp, Avg)

    RETURN

  End Subroutine Freq_Avg

end module FREQ_AVG_M
! $Log$
! Revision 1.4  2001/03/24 01:17:36  livesey
! Bug fix.
!
! Revision 1.3  2001/02/19 22:20:40  zvi
! Latest modification: Conv/NoConv
!
! Revision 1.2  2001/02/19 22:14:21  zvi
!
! Initial conversion to Fortran 90
! 2000/11/20 21:56:09  zvi
! First version D.P.
