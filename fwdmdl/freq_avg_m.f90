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
  subroutine Freq_Avg_DACS ( F_grid, DACSFilter, Rad, Avg )

    use D_CSPLINE_M, only: CSPLINE
    use D_HUNT_M, only: HUNT
    use DFFT_M, only: DTCST
    use FilterShapes_m, only: DACSFilterShape_T
    use MLSCommon, only: I4, R8, RP
    use SineTables_m, only: CreateSineTable, n_sine => LogSize_SineTable_R8, &
      & sines => SineTable_R8

    real(r8), intent(in) :: F_grid(:) ! Frequency grid
    type(DACSFilterShape_T), intent(in) :: DACSFilter
    real(rp), intent(in) :: Rad(:)    ! Radiances on F_Grid

    real(rp), intent(out)   :: Avg(:) ! Radiances for all channels

    integer :: L_Apod, L_Filter, L_Norm ! Logarithms base 2 of transform lengths
    integer :: L_Sines                  ! Log_2 of sine table length

    integer(i4) :: I, Khi, Klo, N, Nfilter, Nnorm
    real(r8) :: Fmax, Fmin, Rmax, Rmin, Tmpary(size(DACSFilter%filterGrid))

    nFilter = size(DACSFilter%filterGrid)
    nNorm = size(DACSFilter%ch_norm)
    i = nFilter / 2

    l_apod = DACSFilter%logApod
    l_filter = DACSFilter%logFilter
    l_norm = DACSFilter%logNorm
    l_sines = max(l_filter, l_norm)

    ! Allocate sine table.  CreateSineTable stores l_sines + 1 into n_sine.
    if ( l_sines + 1 > n_sine ) call createSineTable ( l_sines - 1 )

    if ( DACSFilter%filterGrid(i+1) > DACSFilter%filterGrid(i) ) then
      Fmin = DACSFilter%filterGrid(001)
      Fmax = DACSFilter%filterGrid(nFilter)
    else
      Fmin = DACSFilter%filterGrid(nFilter)
      Fmax = DACSFilter%filterGrid(001)
    end if

    n = size(f_grid)
    klo = -1
    call Hunt ( Fmin, F_grid, n, klo, i )
    call Hunt ( Fmax, F_grid, n, i, khi )

    rmin = minval(Rad(klo:khi))
    rmax = maxval(Rad(klo:khi))

    ! Interpolate from (x=F_grid, y=Rad) to (x=DACSFilter%filterGrid, y=tmpary)
    call Cspline ( F_grid, DACSFilter%filterGrid, Rad, tmpary, n, nFilter, Rmin, Rmax )

    if ( l_apod == l_filter ) then
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      tmpary = tmpary * DACSFilter%lo_apod
      call dtcst ( tmpary, 'C', 'S', (/ l_filter /), 1, n_sine, sines )
      tmpary = tmpary * DACSFilter%filterShape
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      call dtcst ( tmpary(:nNorm), 'C', 'S', (/ l_norm /), 1, n_sine, sines )
    else
      tmpary = tmpary * DACSFilter%filterShape
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      tmpary(:nNorm) = tmpary(:nNorm) * DACSFilter%lo_apod
      call dtcst ( tmpary(:nNorm), 'C', 'S', (/ l_norm /), 1, n_sine, sines )
    end if
    ! The factor of real(nFilter-1)/(nNorm-1) accounts for different normalization
    ! for forward and inverse transforms of different lengths.
    avg = real(nFilter-1)/(nNorm-1) * tmpary(:nNorm) / DACSFilter%ch_norm

  end subroutine Freq_Avg_DACS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Freq_Avg_m

! $Log$
! Revision 2.8  2004/02/12 02:20:22  vsnyder
! Use SineTables_m for sine tables for FFTs
!
! Revision 2.7  2004/02/03 02:49:23  vsnyder
! Work on DACs frequency convolution
!
! Revision 2.6  2003/07/15 23:07:05  vsnyder
! Simplify Freq_Avg, implement Freq_Avg_DACS
!
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
