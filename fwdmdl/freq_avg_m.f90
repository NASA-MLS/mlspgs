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

    use D_CSPLINE_M, only: CSPLINE
    use D_HUNT_M, only: HUNT
    use DFFT_M, only: DTCST
    use MLSCommon, only: I4, R8, RP
    use SineTables_m, only: CreateSineTable, n_sine => LogSize_SineTable_R8, &
      & sines => SineTable_R8

    real(r8), intent(in) :: F_grid(:), F_grid_fltr(:), Fltr_func(:)
    real(r8), intent(in) :: LO_Apod(:), CH_Norm(:)
    real(rp), intent(in) :: Rad(:)

    real(rp), intent(out)   :: Avg(:) ! Radiances for all channels

    ! Saved stuff for DTCST
    integer, save :: L_Long, L_Short  ! Log_2 of lengths

    integer(i4) :: I, Khi, Klo, N, Nfp, Nshort, NSines
    real(r8) :: Fmax, Fmin, Rmax, Rmin, Tmpary(size(f_grid_fltr))

    n = size(f_grid)
    nfp = size(f_grid_fltr)
    nShort = size(lo_apod)
    i = nfp / 2

    ! Allocate sine table
    nSines = max(nfp,nshort) / 2 - 1
    if ( nSines > 2 ** n_sine - 1 ) then
      call check_size ( nfp, l_long, 'F_Grid_Fltr' )
      call check_size ( nShort, l_short, 'LO_Apod' )
      call createSineTable ( l_long - 1 )
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

    ! Interpolate from (x=F_grid, y=Rad) to (x=F_grid_fltr, y=tmpary)
    call Cspline ( F_grid, F_grid_fltr, Rad, tmpary, n, nfp, Rmin, Rmax )

    tmpary = tmpary * Fltr_func(1:nfp)
    call dtcst ( tmpary, 'C', 'A', (/ l_long /), 1, n_sine, sines )
    tmpary(:nShort) = tmpary(:nShort) * lo_apod
    call dtcst ( tmpary(:nShort), 'C', 'S', (/ l_short /), 1, n_sine, sines )
    ! The factor of real(nfp-1)/(nshort-1) accounts for different normalization
    ! for forward and inverse transforms of different lengths.
    avg = real(nfp-1)/(nshort-1) * tmpary(:nShort) / ch_norm

  contains

    subroutine Check_Size ( ArrSize, LogSize, Name )
      ! Make sure ArrSize is 2**k+1 for some k.  Return that k in LogSize
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      integer, intent(in) :: ArrSize     ! Must be 2**k+1 for some k
      integer, intent(out) :: LogSize    ! Base-2 logarithm ( ArrSize ) = K.
      character(len=*), intent(in) :: Name
      logSize = 0
      do
        logSize = logSize + 1
        if ( 2**logSize + 1 == arrSize ) exit
        if ( 2**logSize > arrSize ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Size(' // name // ') /= 2**k + 1' )
      end do
    end subroutine Check_Size

  end subroutine Freq_Avg_DACS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Freq_Avg_m

! $Log$
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
