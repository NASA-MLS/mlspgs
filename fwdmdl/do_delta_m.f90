! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DO_DELTA_M

  implicit NONE
  private
  public :: PATH_OPACITY, HYD_OPACITY, POLARIZED_PATH_OPACITY

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!--------------------------------------------------  PATH_OPACITY  -----

  subroutine PATH_OPACITY ( DEL_ZETA, SINGULARITY, FUNCT, DS_DH_GL, DH_DZ_GL, &
                     &      INTEGRAL, C_Inds, F_Inds )

! Assume FUNCT, DS_DH_GL and DH_DZ_GL have been evaluated at NG Gauss
! abscissae in each DEL_ZETA interval.  Then for each of those intervals,
! estimate the integral of (FUNCT-SINGULARITY) * DS_DH_GL * DH_DZ_GL
! using the NG-point Gauss quadrature.

    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: Del_zeta(:) ! difference in integration boundary
                                        ! in -log(p) units
    real(rp), intent(in) :: Singularity(:) ! value of function at lower boundary
    real(rp), intent(in) :: Funct(:)    ! function evaluated on gl integration
                                        ! grid
    real(rp), intent(in) :: Ds_dh_gl(:) ! path length derivative wrt height on
                                        ! gl grid
    real(rp), intent(in) :: Dh_dz_gl(:) ! height derivative wrt zeta on gl grid

! Output

    real(rp), intent(out) :: Integral(:) ! result from integration

! Optional

    integer(ip), intent(in) :: C_inds(:) ! Coarse path inds, for Singularity
    integer(ip), intent(in) :: F_Inds(:) ! Subset of ds_dh_gl, dh_dz_gl

! Internals

    integer(ip) a, aa, i

! Start calculation

    !{ Apply Gauss-Legendre quadrature to the panels indicated by
    !  {\tt C\_inds}.  We remove a singularity (which actually only
    !  occurs at the tangent point) by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
    !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !   \text{d}\zeta$.  The first integral is easy -- it's just
    !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  We don't use it here.
    !  In the second integral, $G(\zeta)$ is {\tt funct} -- which has
    !  already been evaluated at the appropriate abscissae -- and
    !  $G(\zeta_i)$ is {\tt singularity}.  The weights are {\tt gw}.

    a = 1
    do i = 1, size(c_inds)
      aa = f_inds(a)
      integral(i) = 0.5_rp * del_zeta(i) * &
               &  sum((funct(a:a+ng-1) - singularity(c_inds(i))) * &
               &     ds_dh_gl(aa:aa+ng-1) * dh_dz_gl(aa:aa+ng-1) * gw)
      a = a + ng
    end do

  end subroutine PATH_OPACITY

! ---------------------------------------  POLARIZED_PATH_OPACITY  -----

  subroutine POLARIZED_PATH_OPACITY ( DEL_ZETA, SINGULARITY, &
                     &      FUNCT, DS_DH_GL, DH_DZ_GL, &
                     &      INTEGRAL, C_Inds, F_Inds )
    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: del_zeta(:) ! difference in integration boundary
                                        ! in -log(p) units
!   logical, intent(in) :: do_gl(:)  ! Where INTEGRAL needs to be evaluated
    complex(rp), intent(in) :: singularity(-1:,:) ! value of function at lower
    complex(rp), intent(in) :: funct(-1:,:)    ! function evaluated on gl integration
                                        ! grid
    real(rp), intent(in) :: ds_dh_gl(:) ! path length derivative wrt height on
                                        ! gl grid
    real(rp), intent(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid

! Output

    complex(rp), intent(out) :: integral(-1:,:) ! result from integration

! Optional

    integer(ip), intent(in) :: C_inds(:) ! Coarse path inds, for
                                         ! Singularity
    integer(ip), intent(in) :: F_Inds(:) ! Subset of ds_dh_gl, dh_dz_gl

! Internals

    integer(ip) a, aa, i, j
    complex(rp) :: ds_dh_dh_dz_gw(ng) ! ds_dh_gl * dh_dz_gl * gw is the same
      !                                 for all j = -1..1

! Start calculation

    !{ Apply Gauss-Legendre quadrature to the panels indicated by
    !  {\tt C\_inds}.  We remove a singularity (which actually only
    !  occurs at the tangent point) by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
    !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !   \text{d}\zeta$.  The first integral is easy -- it's just
    !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  We don't use it here.
    !  In the second integral, $G(\zeta)$ is {\tt funct} -- which has
    !  already been evaluated at the appropriate abscissae -- and
    !  $G(\zeta_i)$ is {\tt singularity}.  The weights are {\tt gw}.

    a = 1
    do i = 1, size(c_inds)
      aa = f_inds(a)
      do j = -1, 1
        integral(j,i) = 0.5_rp * del_zeta(c_inds(i)) * &
               &  sum((funct(j,a:a+ng-1) - singularity(j,c_inds(i))) * &
               &  ds_dh_gl(aa:aa+ng-1) * dh_dz_gl(aa:aa+ng-1) * gw)
      end do 
      a = a + ng
    end do

  end subroutine POLARIZED_PATH_OPACITY


!---------------------------------------------------  HYD_OPACITY  -----

  subroutine HYD_OPACITY ( DEL_ZETA, SINGULARITY, ALPHA_PATH, H_PATH,     &
                     &     DH_DT_PATH, T_PATH, H_TAN, DH_DT_TAN, ETA_ZXP, &
                     &     DS_DH_GL, DH_DZ_GL, C_inds, F_inds, INTEGRAL )
    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: del_zeta(:) ! difference in integration boundary
                                        ! in -log(p) units
    real(rp), intent(in) :: singularity(:) ! value of function at lower
                                        ! integration boundary.
    real(rp), intent(in) :: alpha_path(:) ! absorption coefficient on gl grid.
    real(rp), intent(in) :: h_path(:)   ! path heights + req on gl grid.
    real(rp), intent(in) :: dh_dt_path(:) ! path height derivative wrt
                                        ! temperature coefficient on gl grid.
    real(rp), intent(in) :: t_path(:)   ! path temperature on gl grid.
    real(rp), intent(in) :: h_tan       ! tangent height + req.
    real(rp), intent(in) :: dh_dt_tan   ! height derivative wrt temperature
                                        ! coefficient at the tangent.
    real(rp), intent(in) :: eta_zxp(:)  ! basis function for temperature
                                        ! coefficient on gl grid.
    real(rp), intent(in) :: ds_dh_gl(:) ! path length derivative wrt height on
                                        ! gl grid on gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid.
    integer(ip), intent(in) :: C_inds(:) ! Coarse path inds, for Singularity
    integer(ip), intent(in) :: F_inds(:) ! The first GL path ind for the
                                        ! corresponding C_inds element, for
                                        ! other than Del_zeta and Singularity.
                                        ! Something like all_inds(::ng)
! Output

    real(rp), intent(out) :: integral(:) ! result from integration

! Internals

    integer(ip) a, b, i

! Start calculation

    do i = 1 , size(c_inds)
      a = f_inds(i)
      b = a + ng - 1
      integral(i) = 0.5_rp * del_zeta(c_inds(i)) *                             &
        & sum( ( alpha_path(a:b) - singularity(c_inds(i)) )            &
        &      * (((2.0_rp*h_path(a:b)**2 - 3.0_rp*h_tan**2)           &
        &         * dh_dt_path(a:b) + h_path(a:b) * h_tan * dh_dt_tan) &
        &         / (sqrt(h_path(a:b)**2 - h_tan**2))**3               &
        &         + eta_zxp(a:b)*ds_dh_gl(a:b)/t_path(a:b)) *          &
        &         dh_dz_gl(a:b) * gw )
    end do

  end subroutine HYD_OPACITY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DO_DELTA_M
!---------------------------------------------------
! $Log$
! Revision 2.9  2003/09/26 01:25:56  vsnyder
! Restore a deleted CVS log comment
!
! Revision 2.8  2003/09/25 20:04:48  vsnyder
! Insert TeXnicalities
!
! Revision 2.7  2003/09/24 22:19:55  vsnyder
! Get rid of some array temps
!
! Revision 2.6  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.5  2003/05/15 03:25:20  vsnyder
! Added some comments
!
! Revision 2.4  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.3.2.2  2003/03/06 21:52:39  vsnyder
! Use the correct size
!
! Revision 2.3.2.1  2003/03/05 03:40:35  vsnyder
! More polarized work
!
! Revision 2.3  2003/02/07 00:39:58  michael
! add polarized_path_opacity
!
! Revision 2.2  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/10/02 20:08:16  vsnyder
! Insert copyright notice, other cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.1.2.2  2001/09/12 21:38:45  zvi
! Added CVS stuff
!
!
