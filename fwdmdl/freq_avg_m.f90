! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Freq_Avg_m

  implicit NONE
  private
  public :: Freq_Avg, Freq_Avg_Avg, Freq_Avg_DACS, Freq_Avg_Setup

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  Freq_Avg  -----
  subroutine Freq_Avg ( F_grid, F_grid_fltr, Fltr_func, Rad, Avg )

    use MLSKinds, only: R8, RP

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(rp), intent(in) :: Rad(:)

    real(rp), intent(out)   :: Avg

    integer :: Khi, Klo
    real(r8) :: dF

    call freq_avg_setup ( f_grid, f_grid_fltr, klo, khi, dF )
    call freq_avg_avg ( f_grid, f_grid_fltr, fltr_func, rad, klo, khi, dF, avg )

  end subroutine Freq_Avg

  ! -----------------------------------------------  Freq_Avg_Avg  -----
  subroutine Freq_Avg_Avg ( F_grid, F_grid_fltr, Fltr_func, Rad, Klo, Khi, dF, &
                          & Avg )

    ! use D_CSPLINE_M, only: CSPLINE
    use MLSNumerics, only: InterpolateValues, Simps => SimpsonsSub
    use MLSKinds, only: R8, RP

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(rp), intent(in) :: Rad(:)
    integer, intent(in) :: KLo, KHi
    real(r8), intent(in) :: dF

    real(rp), intent(out)   :: Avg

    integer :: Nfp
    real(r8) :: Rmax, Rmin, Tmpary(size(f_grid_fltr))

    nfp = size(f_grid_fltr)

    rmin = minval(Rad(klo:khi))
    rmax = maxval(Rad(klo:khi))

    ! InterpolateValues args are:
    ! oldX, oldY, newX, newY, method
    call InterpolateValues ( F_grid, Rad, F_grid_fltr, tmpary, method='C', &
      & YMIN=Rmin, YMAX=Rmax )

    ! Since the filter function and its derivatives are zero at the ends
    ! of the interval, or at least outside of it, use the Euler-Maclaurin
    ! summation formula.  Why doesn't this work?

!     avg = dF * dot_product(tmpary, fltr_func(:nfp) )

    tmpary = tmpary * Fltr_func(1:nfp)
    call Simps ( tmpary, dF, nfp, Avg )

  end subroutine Freq_Avg_Avg

  ! ----------------------------------------------  Freq_Avg_DACS  -----
  subroutine Freq_Avg_DACS ( F_grid, DACSFilter, Rad, Avg )

    ! use D_CSPLINE_M, only: CSPLINE
    use DFFT_M, only: DTCST
    use FilterShapes_m, only: DACSFilterShape_T
    use MLSKinds, only: I4, R8, RP
    use MLSNumerics, only: InterpolateValues
    use Pure_Hunt_m, only: PureHunt
    use SineTables_m, only: CreateSineTable, N_Sine => Logsize_SineTable_R8, &
      & Sines => SineTable_R8

    real(r8), intent(in) :: F_grid(:) ! Frequency grid
    type(DACSFilterShape_T), intent(in) :: DACSFilter
    real(rp), intent(in) :: Rad(:)    ! Radiances on F_Grid

    real(rp), intent(out)   :: Avg(:) ! Radiances for all channels

    integer :: L_Apod, L_Filter, L_Norm ! Logarithms base 2 of transform lengths
    integer :: L_Sines                  ! Log_2 of sine table length

    integer(i4) :: I, Khi, Klo, N, Nfilter, Nnorm
    real(r8) :: Fmax, Fmin, Rmax, Rmin, Tmpary(size(DACSFilter%filter%filterGrid))

    nFilter = size(DACSFilter%filter%filterGrid)
    nNorm = size(DACSFilter%ch_norm)
    i = nFilter / 2

    l_apod = DACSFilter%logApod
    l_filter = DACSFilter%logFilter
    l_norm = DACSFilter%logNorm
    l_sines = max(l_filter, l_norm)

    ! Allocate sine table.  CreateSineTable stores l_sines + 1 into n_sine.
    if ( l_sines + 1 > n_sine ) call createSineTable ( l_sines - 1 )

    if ( DACSFilter%filter%filterGrid(i+1) > DACSFilter%filter%filterGrid(i) ) then
      Fmin = DACSFilter%filter%filterGrid(001)
      Fmax = DACSFilter%filter%filterGrid(nFilter)
    else
      Fmin = DACSFilter%filter%filterGrid(nFilter)
      Fmax = DACSFilter%filter%filterGrid(001)
    end if

    n = size(f_grid)
    klo = -1
    call purehunt ( Fmin, F_grid, n, klo, i )
    call purehunt ( Fmax, F_grid, n, i, khi )

    rmin = minval(Rad(klo:khi))
    rmax = maxval(Rad(klo:khi))

    ! Interpolate from (x=F_grid, y=Rad) to (x=DACSFilter%filter%filterGrid, y=tmpary)
    ! call Cspline ( F_grid, DACSFilter%filter%filterGrid, Rad, tmpary, n, &
    !   & nFilter, Rmin, Rmax )
    call InterpolateValues ( F_grid, Rad, DACSFilter%filter%filterGrid, tmpary, method='C', &
      & YMIN=Rmin, YMAX=Rmax )

    if ( l_apod == l_filter ) then
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      tmpary = tmpary * DACSFilter%lo_apod
      call dtcst ( tmpary, 'C', 'S', (/ l_filter /), 1, n_sine, sines )
      tmpary = tmpary * DACSFilter%filter%filterShape
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      call dtcst ( tmpary(:nNorm), 'C', 'S', (/ l_norm /), 1, n_sine, sines )
    else
      tmpary = tmpary * DACSFilter%filter%filterShape
      call dtcst ( tmpary, 'C', 'A', (/ l_filter /), 1, n_sine, sines )
      tmpary(:nNorm) = tmpary(:nNorm) * DACSFilter%lo_apod
      call dtcst ( tmpary(:nNorm), 'C', 'S', (/ l_norm /), 1, n_sine, sines )
    end if
    ! The factor of real(nFilter-1)/(nNorm-1) accounts for different normalization
    ! for forward and inverse transforms of different lengths.
    avg = real(nFilter-1)/(nNorm-1) * tmpary(:nNorm) / DACSFilter%ch_norm

  end subroutine Freq_Avg_DACS

  ! ---------------------------------------------  Freq_Avg_Setup  -----
  subroutine Freq_Avg_Setup ( F_grid, F_grid_fltr, Klo, Khi, dF )

    ! Determine which frequencies from F_Grid to use to span F_Grid_Fltr

    use Pure_Hunt_m, only: PureHunt
    use MLSKinds, only: R8

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:)

    integer, intent(out)   :: Klo, Khi
    real(r8), intent(out) :: dF

    integer :: I, N, Nfp
    real(r8) :: Fmax, Fmin

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
    call purehunt ( Fmin, F_grid, n, klo, i )
    call purehunt ( Fmax, F_grid, n, i, khi )

  end subroutine Freq_Avg_Setup

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Freq_Avg_m

! $Log$
! Revision 2.20  2013/06/12 02:24:58  vsnyder
! Cruft removal
!
! Revision 2.19  2013/05/22 00:09:12  vsnyder
! Remove unreferenced use names
!
! Revision 2.18  2011/08/26 17:53:53  pwagner
! purehunt recovers optimized functionality of fwdmdls own hunt
!
! Revision 2.17  2011/08/26 01:21:16  pwagner
! Fixed obvious bugs in call to Hunt
!
! Revision 2.16  2011/08/26 00:31:09  pwagner
! CSpline and Hunt now USE MLSNumerics
!
! Revision 2.15  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.14  2006/04/25 23:25:37  vsnyder
! Revise DACS filter shape data structure
!
! Revision 2.13  2005/10/24 20:24:55  vsnyder
! Insert Euler-Maclaurin idea as a comment.  Why doesn't it work?
!
! Revision 2.12  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.11  2004/12/28 00:26:26  vsnyder
! Remove unreferenced use names
!
! Revision 2.10  2004/10/06 21:25:45  vsnyder
! Split Freq_Avg into 'setup' and 'do it' steps
!
! Revision 2.9  2004/02/14 00:23:48  vsnyder
! New DACS convolution algorithm
!
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
