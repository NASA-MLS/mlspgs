! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FREQ_AVG_M
  use D_CSPLINE_M, only: CSPLINE
  use DSIMPSON_MODULE, only: SIMPS
  use D_HUNT_M, only: HUNT
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

    Integer(i4), intent(IN) :: n, nfp
!
    Real(r8), intent(OUT)   :: Avg
!
    Integer(i4) :: klo,khi,i
    Real(r8) :: rxf(nfp), tmpary(nfp), Fmin, Fmax, Rmin, Rmax, dF
!
    i = nfp / 2
    dF = F_grid_fltr(i+1) - F_grid_fltr(i)
    if(dF > 0.0_r8) then
      Fmin = F_grid_fltr(001)
      Fmax = F_grid_fltr(nfp)
    else
      dF = -dF
      Fmin = F_grid_fltr(nfp)
      Fmax = F_grid_fltr(001)
    endif

    klo = -1
    Call Hunt(Fmin,F_grid,n,klo,i)
    Call Hunt(Fmax,F_grid,n,i,khi)
!
    Rmin = MINVAL(Rad(klo:khi))
    Rmax = MAXVAL(Rad(klo:khi))
!
    Call Cspline(F_grid, F_grid_fltr, Rad, tmpary, n, nfp, Rmin, Rmax)

    rxf(1:nfp) = tmpary(1:nfp) * Fltr_func(1:nfp)
    Call Simps (rxf, dF, nfp, Avg)

    RETURN

  End Subroutine Freq_Avg

end module FREQ_AVG_M
! $Log$
! Revision 2.1  2002/04/18 10:46:26  zvi
! Better spline use
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.6  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.5  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
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
