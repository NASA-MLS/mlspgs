! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Freq_Avg_m

  implicit NONE
  private
  public :: Freq_Avg, Freq_Avg_DACS

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  Freq_Avg  -----
  subroutine Freq_Avg ( F_grid, F_grid_fltr, Fltr_func, Rad, Avg )

    use D_CSPLINE_M, only: CSPLINE
    use DSIMPSON_MODULE, only: SIMPS
    use D_HUNT_M, only: HUNT
    use MLSCommon, only: I4, R8, RP

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(rp), intent(in) :: Rad(:)

    real(rp), intent(out)   :: Avg

    integer(i4) :: I, Khi, Klo, N, Nfp
    real(r8) :: dF, Fmax, Fmin, Rmax, Rmin, Tmpary(size(f_grid_fltr))

    n = size(f_grid)
    nfp = size(f_grid_fltr)
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

    tmpary = tmpary * Fltr_func(1:nfp)
    call Simps ( tmpary, dF, nfp, Avg )

    return

  end subroutine Freq_Avg

  ! ----------------------------------------------  Freq_Avg_DACS  -----
  subroutine Freq_Avg_DACS ( F_grid, F_grid_fltr, Fltr_func, LO_Apod, CH_Norm, &
    & Rad, Avg )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use D_CSPLINE_M, only: CSPLINE
    use D_HUNT_M, only: HUNT
    use DFFT_M, only: DTCST
    use MLSCommon, only: I4, R8, RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(r8), intent(in) :: LO_Apod(:), CH_Norm(:)
    real(rp), intent(in) :: Rad(:)

    real(rp), intent(out)   :: Avg(:) ! Radiances for all channels

    ! Saved stuff for DTCST
    integer, save :: L_Long, L_Short                 ! Log_2 of lengths
    integer, save :: MS_Long = 0, MS_Short = 0       ! Sine table lengths
    real(r8), save, pointer :: S_Long(:), S_Short(:) ! Sine tables

    integer(i4) :: I, Khi, Klo, N, Nfp, Nshort
    real(r8) :: Fmax, Fmin, Rmax, Rmin, Tmpary(size(f_grid_fltr))

    n = size(f_grid)
    nfp = size(f_grid_fltr)
    nShort = size(lo_apod)
    i = nfp / 2

    ! Allocate sine tables
    if ( ms_long == 0 ) then
      nullify ( s_long )
    else if ( ms_long /= size(s_long) ) then
      call deallocate_test ( s_long, 'S_Long', moduleName )
      ms_long = 0
    end if
    if ( ms_long == 0 ) then
      call allocate_test ( s_long, nfp, 'S_Long', moduleName )
      l_long = 0
      do
        l_long = l_long + 1
        if ( 2**l_long + 1 == nfp ) exit
        if ( 2**l_long > nfp ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Size(F_Grid_Fltr) /= 2**k + 1' )
      end do
    end if
    if ( ms_short == 0 ) then
      nullify ( s_short )
    else if ( ms_short /= size(lo_apod) ) then
      call deallocate_test ( s_short, 'S_short', moduleName )
      ms_short = 0
    end if
    if ( ms_short == 0 ) then
      call allocate_test ( s_short, nShort, 'S_short', moduleName )
      l_short = 0
      do
        l_short = l_short + 1
        if ( 2**l_short + 1 == nShort ) exit
        if ( 2**l_short > nShort ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Size(LO_Apod) /= 2**k + 1' )
      end do
    end if

    if ( F_grid_fltr(i+1) > F_grid_fltr(i) ) then
      Fmin = F_grid_fltr(001)
      Fmax = F_grid_fltr(nfp)
    else
      Fmin = F_grid_fltr(nfp)
      Fmax = F_grid_fltr(001)
    end if

    klo = -1
    call Hunt ( Fmin, F_grid, n, klo, i )
    call Hunt ( Fmax, F_grid, n, i, khi )

    rmin = minval(Rad(klo:khi))
    rmax = maxval(Rad(klo:khi))

    call Cspline ( F_grid, F_grid_fltr, Rad, tmpary, n, nfp, Rmin, Rmax )

    tmpary = tmpary * Fltr_func(1:nfp)
    call dtcst ( tmpary, 'C', 'A', (/ l_long /), 1, ms_long, s_long )
    tmpary(:nShort) = tmpary(:nShort) * lo_apod
    call dtcst ( tmpary(:nShort), 'C', 'S', (/ l_short /), 1, ms_short, s_short )
    avg = tmpary(:nShort) / ch_norm

    return

  end subroutine Freq_Avg_DACS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Freq_Avg_m

! $Log$
! Revision 2.5  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
! Changing the (*) to (:) everywhere
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
