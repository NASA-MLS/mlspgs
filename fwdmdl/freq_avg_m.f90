! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Freq_Avg_m

  implicit NONE
  private
  public :: Freq_Avg

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Freq_Avg ( F_grid, F_grid_fltr, Fltr_func, Rad, N, Nfp, Avg )

    use D_CSPLINE_M, only: CSPLINE
    use DSIMPSON_MODULE, only: SIMPS
    use D_HUNT_M, only: HUNT
    use MLSCommon, only: I4, R8, RP

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(rp), intent(in) :: Rad(:)

    integer(i4), intent(in) :: N, Nfp

    real(rp), intent(out)   :: Avg

    integer(i4) :: Klo, Khi, I
    real(r8) :: Rxf(nfp), Tmpary(nfp), Fmin, Fmax, Rmin, Rmax, dF

    i = nfp / 2
    dF = F_grid_fltr(i+1) - F_grid_fltr(i)
    if ( dF > 0.0_r8 ) then
      Fmin = F_grid_fltr(001)
      Fmax = F_grid_fltr(nfp)
    else
      dF = -dF
      Fmin = F_grid_fltr(nfp)
      Fmax = F_grid_fltr(001)
    end if

    klo = -1
    call Hunt ( Fmin, F_grid, n, klo, i )
    call Hunt ( Fmax, F_grid, n, i, khi )

    rmin = minval(Rad(klo:khi))
    rmax = maxval(Rad(klo:khi))

    call Cspline ( F_grid, F_grid_fltr, Rad, tmpary, n, nfp, Rmin, Rmax )

    rxf(1:nfp) = tmpary(1:nfp) * Fltr_func(1:nfp)
    call Simps ( rxf, dF, nfp, Avg )

    return

  end subroutine Freq_Avg

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Freq_Avg_m

! $Log$
! Revision 2.4  2002/09/07 02:20:27  vsnyder
! Fix a type
!
! Revision 2.3  2002/09/07 02:18:37  vsnyder
! Move USEs from module scope to procedure scope, cosmetic changes
!
! Revision 2.2  2002/05/08 08:53:46  zvi
! Modify to accomodate cspline.f9h
!
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
