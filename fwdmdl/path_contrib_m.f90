! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Path_Contrib_M

  implicit NONE
  private
  public :: Get_GL_Inds, Path_Contrib

  interface Get_GL_Inds
    module procedure Get_GL_Inds_and_Where, Get_GL_Inds_Etc, Get_GL_Inds_Only
  end interface

  interface Path_Contrib
    module procedure Path_Contrib_Scalar, Path_Contrib_Polarized
  end interface

  ! If dTaudn < Black_Out, we haven't been using GL, even if Tol < 0.  To do
  ! GL everywhere, even where delta < Black_Out (dTaudn is at first used for
  ! delta), except where dTaudn = exp(delta) is zero, set GL_Everywhere true.
  ! I tried this on a few cases, and it didn't make a noticeable difference in
  ! radiances or derivatives.
  logical, parameter, private :: GL_Everywhere = .false.

  ! The dummy argument I_End could be updated to stop worrying about
  ! incremental optical depth after black out.  At some time in the distant
  ! past, a comment was added above the call in the full forward model to
  ! Tau_m%Get_Tau that calculates Tau from incremental optical depth: "this
  ! breaks the gold brick."  Tau_m%Get_Tau re-calculates delta, and might black
  ! out at a different point on the path.  Set this .true. to update I_End
  ! anyway.
  logical, parameter, private :: Update_I_End = .false.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!-------------------------------------------  Path_Contrib_Scalar  -----
! Estimate the error (along the path) of the integral in each interval
! of the (coarse) pre-selected integration grid.  Use that estimate to
! select where to do Gauss-Legendre quadrature.  Then specify where to
! do the Gauss-Legendre quadratures.

  subroutine Path_Contrib_Scalar ( incoptdepth, tan_pt, i_start, i_end, &
    &                              e_rflty, tol, do_gl )

    use MLSKinds, only: RK => RP, IP
    use Tau_m, only: Black_Out

  ! inputs

    real(rk), intent(in) :: IncOptDepth(:) ! layer optical depth
    integer, intent(in) :: tan_pt          ! Tangent point index in IncOptDepth
    integer, intent(in) :: i_start         ! Start of path to worry about
    integer, intent(inout) :: i_end        ! End of path to worry about
    real(rk), intent(in) :: e_rflty        ! earth reflectivity
    real(rk), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical, intent(inout) :: do_gl(:)     ! set true for indices in coarse
                                           ! path to do gl computation

  ! Internal stuff

    real(rk) :: dtaudn(size(incoptdepth))  ! optical depth (delta) at first,
                                           ! then path derivative of the
                                           ! transmission function (tau)
    integer(ip) :: last

    real(rk) :: MyTol
    real(rk), parameter :: temp = 250.0_rk
    real(rk), parameter :: TolScale = 2.0_rk / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

    real(rk), parameter :: My_Black_Out = &
      & merge(-log(huge(1.0_rk)),black_out,GL_Everywhere)

  ! start code

    myTol = - tolScale * tol ! Negative because we're summing -incoptdepth

    ! Compute dtaudn = delta = indefinite sum of (-incoptdepth).

o:  block
      dtaudn(i_start) = 0.0_rk
      do last = i_start+1, min(tan_pt, i_end)
        dtaudn(last) = dtaudn(last-1) - incoptdepth(last)
        if ( dtaudn(last) < my_black_out ) exit o
      end do

      if ( last == tan_pt+1 ) dtaudn(tan_pt+1) = dtaudn(tan_pt)

      do last = max(i_start+1,tan_pt+2), i_end
        dtaudn(last) = dtaudn(last-1) - incoptdepth(last-1)
        if ( dtaudn(last) < my_black_out ) exit o
      end do
      last = i_end
    end block o

    ! compute the tau path derivative dTau/ds ~ exp(delta) d delta/ds.
    dtaudn(i_start:last) = &
      & (eoshift(dtaudn(i_start:last),1,dtaudn(last)) -             &
      &  eoshift(dtaudn(i_start:last),-1,dtaudn(i_start))) * &
      & exp(dtaudn(i_start:last))

    dtaudn(tan_pt+1:last) = dtaudn(tan_pt+1:last) * e_rflty

    ! find where the tau derivative is large.  Remember, tau is monotone
    ! decreasing, so "large" means "large and negative."

    do_gl(:i_start) = .false. ! irrelevant
    do_gl(i_start+1:last) = dtaudn(i_start+1:last) < myTol .or. GL_Everywhere
    do_gl(last+1:) = GL_Everywhere  ! Tau is blacked out, so no point in doing GL
    if ( update_i_end ) i_end = last

  end subroutine Path_Contrib_Scalar

!----------------------------------------  Path_Contrib_Polarized  -----
! Estimate the contributions (along the path) of each interval of the
! (coarse) pre-selected integration grid.  Use that estimate to select
! where to do Gauss-Legendre quadrature.  Then allocate and fill
! arrays that control the Gauss-Legendre quadratures.

  subroutine Path_Contrib_Polarized ( deltau, tan_pt, e_rflty, tol, do_gl )

!   use CS_Expmat_M, only: CS_Expmat
    use MLSCommon, only: RK => RP, IP

  ! inputs

!    complex(rk), intent(in) :: incoptdepth(:,:,:) ! layer optical depth
    complex(rk), intent(in) :: deltau(:,:,:)    ! = E == exp(-incoptdepth)
                                           ! (2,2,:)
    integer, intent(in) :: tan_pt          ! Tangent point index in deltau
    real(rk), intent(in) :: e_rflty        ! earth reflectivity
    real(rk), intent(in) :: tol            ! accuracy target in K

  ! outputs

    logical, intent(inout) :: do_gl(:)     ! set true for indices in coarse
                                           ! path to do gl computation

  ! Internal stuff

    complex(rk) :: dtaudn(2,2)             ! path derivative of the
                                           ! transmission function
    complex(rk), parameter :: Ident(2,2) = reshape( (/ 1.0, 0.0, &
      &                                                0.0, 1.0 /), (/ 2,2 /) )
    real(rk) :: MyTol
    complex(rk) :: P(2,2,size(deltau,3)), Tau(2,2,size(deltau,3))

    integer(ip) :: i, n_path

    real(rk), parameter :: temp = 250.0_rk
    real(rk), parameter :: TolScale = 2.0_rk / temp ! 2.0 comes from centered
                                           ! difference used to compute dtaudn

  ! start code

    n_path = size(deltau,3)
    myTol = tolScale * tol

  ! Compute exp(incoptdepth) for all but the last level
  ! (now done outside)

  !  do i = 1, n_path - 1
  !    call cs_expmat ( incoptdepth(:,:,i), deltau(:,:,i) )
  !  end do

    P(:,:,1) = ident
    Tau(:,:,1) = ident

  ! Multiply the exp(incoptdepth) matrices together.

  !{ $\mathbf{P}_i = \prod_{j=1}^{i-1} \mathbf{E}_i$;
  !  $\mathbf{\tau}_i = \mathbf{P}_i \mathbf{P}_i^\dagger$.

    do i = 2, tan_pt
  !Note: Indexing of deltau changes at the tangent point because it is a "layer quantity"
  !      while P, Tau are defined at the boundary of a layer closest to the spacecraft  
!      P(:,:,i) =  matmul ( P(1:2,1:2,i-1),  deltau(1:2,1:2,i-1) )
      P(:,:,i) =  matmul ( P(1:2,1:2,i-1),  deltau(1:2,1:2,i) )
      Tau(:,:,i) = matmul ( P(1:2,1:2,i), conjg(transpose(P(1:2,1:2,i))) )
    end do

    P(:,:,tan_pt+1) = P(:,:,tan_pt) * sqrt(e_rflty)
    Tau(:,:,tan_pt+1) = Tau(:,:,tan_pt) * e_rflty

    do i = tan_pt+2, n_path
      P(:,:,i) =  matmul ( P(1:2,1:2,i-1),  deltau(1:2,1:2,i-1) )
      Tau(:,:,i) = matmul ( P(1:2,1:2,i), conjg(transpose(P(1:2,1:2,i))) )
    end do

  ! Where is the derivative of that product large?

    do i = 2, n_path-1
      dtaudn = tau(:,:,i+1) - tau(:,:,i-1)
      if ( any(abs(real(dtaudn))  >= myTol ) .or. &
           any(abs(aimag(dtaudn)) >= myTol ) ) do_gl(i) = .true.
    end do

  end subroutine Path_Contrib_Polarized

  ! --------------------------------------  Get_GL_Inds_and_Where  -----
  subroutine Get_GL_Inds_and_Where ( Do_GL, Tan_pt, GL_Inds, NGL, CG_Inds, NCG, Where_GL )
  ! Fill the arrays that control application of GL

    use GLnp, only: NG, NGP1
    use MLSCommon, only: IP

    logical, intent(in) :: DO_GL(:)         ! Set true for indices in coarse
                                            ! path to do gl computation. 
                                            ! First and last are set false
                                            ! here.
    integer, intent(in) :: Tan_Pt           ! Index of tangent point in Do_GL
    integer(ip), intent(out) :: GL_Inds(:)  ! Indices of GL points within
                                            ! coarse & fine path.
    integer(ip), intent(out) :: NGL         ! How much of GL_Inds to use
    integer(ip), intent(out) :: CG_Inds(:)  ! If K=CG_Inds(i) <= Tan_Pt then
                                            ! panel(k-1:k) on the coarse path
                                            ! needs GL, else panel(k:k+1) needs
                                            ! GL.
    integer(ip), intent(out) :: NCG         ! How much of CG_Inds to use
    integer(ip), intent(out) :: Where_GL(:) ! K = Where_GL(i) indicates
                                            ! GL_Inds(k:k+ng-1) are GL_Inds
                                            ! for panel(i:i+1) in the coarse
                                            ! path.

    integer :: I, I1, I2, J, N_Path, Offset

    integer, parameter :: GLIX(ng) = (/ (i ,i = 1-ng, 0) /)

    n_path = size(do_gl)

    ngl = 0
    ncg = 0
    ! The complication here arises from two sources.  First, we never do
    ! GL between the two tangent points, so we never insert GL points
    ! between tan_pt and tan_pt+1.  Second, before the tangent point, Do_GL(k)
    ! indicates that we need GL on panel(i-1:i), while after the tangent
    ! point, DO_GL indicates that we need GL on panel(i:i+1).
    where_gl(tan_pt) = 0 ! panel(tan_pt:tan_pt+1) never needs GL
    i1 = 2
    i2 = tan_pt
    offset = 1
    do j = 0, ngp1, ngp1
      do i = i1, i2
        if ( do_gl(i) ) then
          where_gl(i - offset) = ngl + 1
          ngl = ngl + ng
          gl_inds(ngl-ng+1:ngl) = Ngp1 * (i - 1) + glix + j
          ncg = ncg + 1
          cg_inds(ncg) = i
        else
          where_gl(i - offset) = 0
        end if
      end do
      i1 = tan_pt+1
      i2 = n_path - 1
      offset = 0
    end do

  end subroutine Get_GL_Inds_and_Where

  ! --------------------------------------------  Get_GL_Inds_Etc  -----
  subroutine Get_GL_Inds_Etc ( Do_GL, Tan_pt, GL_Inds, NGL, CG_Inds, NCG )
  ! Fill the arrays that control application of GL

    use GLnp, only: NG, NGP1
    use MLSCommon, only: IP

    logical, intent(in) :: DO_GL(:)         ! Set true for indices in coarse
                                            ! path to do gl computation. 
                                            ! First and last are set false
                                            ! here.
    integer, intent(in) :: Tan_Pt           ! Index of tangent point in Do_GL
    integer(ip), intent(out) :: GL_Inds(:)  ! Indices of GL points within
                                            ! coarse & fine path.
    integer(ip), intent(out) :: NGL         ! How much of GL_Inds to use
    integer(ip), intent(out) :: CG_Inds(:)  ! If K=CG_Inds(i) <= Tan_Pt then
                                            ! panel(k-1:k) on the coarse path
                                            ! needs GL, else panel(k:k+1) needs
                                            ! GL.
    integer(ip), intent(out) :: NCG         ! How much of CG_Inds to use

    integer :: I, I1, I2, J, N_Path

    integer, parameter :: GLIX(ng) = (/ (i ,i = 1-ng, 0) /)

    n_path = size(do_gl)

    ngl = 0
    ncg = 0
    ! The complication here arises from two sources.  First, we never do
    ! GL between the two tangent points, so we never insert GL points
    ! between tan_pt and tan_pt+1.  Second, before the tangent point, Do_GL(k)
    ! indicates that we need GL on panel(i-1:i), while after the tangent
    ! point, DO_GL indicates that we need GL on panel(i:i+1).
    i1 = 2
    i2 = tan_pt
    do j = 0, ngp1, ngp1
      do i = i1, i2
        if ( do_gl(i) ) then
          ngl = ngl + ng
          gl_inds(ngl-ng+1:ngl) = Ngp1 * (i - 1) + glix + j
          ncg = ncg + 1
          cg_inds(ncg) = i
        end if
      end do
      i1 = tan_pt+1
      i2 = n_path - 1
    end do

  end subroutine Get_GL_Inds_Etc

  ! -------------------------------------------  Get_GL_Inds_Only  -----
  subroutine Get_GL_Inds_Only ( Do_GL, Tan_pt, GL_Inds, NGL )
  ! Fill the arrays that control application of GL

    use GLnp, only: NG, NGP1
    use MLSCommon, only: IP

    logical, intent(in) :: DO_GL(:)            ! Set true for indices in coarse
                                               ! path to do gl computation. 
                                               ! First and last are set false
                                               ! here.
    integer, intent(in) :: Tan_Pt              ! Index of tangent point in Do_GL
    integer(ip), intent(out) :: GL_Inds(:)     ! Indices of GL points within
                                               ! coarse & fine path.
    integer(ip), intent(out) :: NGL            ! How much of GL_INDS to use

    integer :: I, I1, I2, J, N_Path

    integer, parameter :: GLIX(ng) = (/ (i ,i = 1-ng, 0) /)

    n_path = size(do_gl)

    ngl = 0
    ! The complication here arises from two sources.  First, we never do
    ! GL between the two tangent points, so we never insert GL points
    ! between tan_pt and tan_pt+1.  Second, before the tangent point, Do_GL(k)
    ! indicates that we need GL on panel(i-1:i), while after the tangent
    ! point, DO_GL indicates that we need GL on panel(i:i+1).
    i1 = 2
    i2 = tan_pt
    do j = 0, ngp1, ngp1
      do i = i1, i2
        if ( do_gl(i) ) then
          ngl = ngl + ng
          gl_inds(ngl-ng+1:ngl) = Ngp1 * (i - 1) + glix + j
        end if
      end do
      i1 = tan_pt+1
      i2 = n_path - 1
    end do

  end subroutine Get_GL_Inds_Only

!-----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Path_Contrib_M

! $Log$
! Revision 2.30  2018/10/26 22:04:09  vsnyder
! Make I_End INOUT, add a switch to control whether to update it.  Get
! kinds from MLSKinds instead of MLSCommon.
!
! Revision 2.29  2018/08/28 20:25:58  vsnyder
! Add a named constant to change black out to -log(huge(1.0_rk))
!
! Revision 2.28  2018/05/14 23:31:54  vsnyder
! Add several more Get_GL_Inds routines, make generic for them
!
! Revision 2.27  2017/08/09 20:37:43  vsnyder
! Use a BLOCK with EXIT to eliminate some GO TO statements and a label.
! Hoist test for optional argument out of a loop, making two versions of it.
!
! Revision 2.26  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.25  2011/07/29 01:58:56  vsnyder
! Make CG_INDS, NCG optional
!
! Revision 2.24  2010/02/02 01:29:20  vsnyder
! Don't reference undefined parts of dtaudn
!
! Revision 2.23  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.22  2009/06/13 01:11:07  vsnyder
! Specify start and end of path
!
! Revision 2.21  2007/12/04 01:56:41  vsnyder
! Don't bother with earth reflectivity after black out
!
! Revision 2.20  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.19  2006/06/16 20:32:31  vsnyder
! Define NGP1 in glnp
!
! Revision 2.18  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.17  2003/11/04 02:49:13  vsnyder
! Calculate coarse-path indices where GL is needed
!
! Revision 2.16  2003/10/30 20:32:28  vsnyder
! Simplify gl_inds computation
!
! Revision 2.15  2003/10/09 19:30:36  vsnyder
! Simplify computation of gl_inds
!
! Revision 2.13  2003/08/13 22:19:11  michael
! Fixed indexing bug at line 150 in path_contrib_polarized
!
! Revision 2.12  2003/08/12 23:05:30  vsnyder
! Fix a bug at line 162
!
! Revision 2.11  2003/08/12 19:34:53  michael
! Some futzing by Van.
!
! Revision 2.10  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.9  2003/05/15 23:57:21  michael
! Fixed Path_Contrib_Polarized to properly calculate Tau
! and renamed variables to agree with polarized ATBD
!
! Revision 2.8  2003/02/07 03:26:24  vsnyder
! Compute deltau instead of incoptdepth in Path_Contrib_Polarized
!
! Revision 2.7  2003/02/06 21:36:18  vsnyder
! Cosmic -- errr -- cosmetic changes
!
! Revision 2.6  2003/02/03 22:56:20  vsnyder
! Get rid of an array temp
!
! Revision 2.5  2003/01/31 01:53:28  vsnyder
! Calculate where to do GL with one less array temp
!
! Revision 2.4  2003/01/30 19:31:18  vsnyder
! Undo change that didn't work -- tried to compute gl_inds without array temp
!
! Revision 2.3  2003/01/18 02:22:58  vsnyder
! IMAG should have been AIMAG
!
! Revision 2.2  2003/01/18 01:42:42  vsnyder
! Added complex 3-quantity path_contrib subroutine.
! Separated Get_GL_Inds into a separate subroutine.
!
! Revision 2.1  2003/01/08 00:09:56  vsnyder
! Moved from rad_tran, where it didn't seem to belong
!
