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
                     &      INTEGRAL )

! Assume FUNCT, DS_DH_GL and DH_DZ_GL have been evaluated at NG Gauss
! abscissae in each DEL_ZETA interval.  Then for each of those intervals,
! estimate the integral of (FUNCT-SINGULARITY) * DS_DH_GL * DH_DZ_GL
! using the NG-point Gauss quadrature.

    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: del_zeta(:) ! difference in integration boundary
!                                         in -log(p) units
    real(rp), intent(in) :: singularity(:) ! value of function at lower boundary
    real(rp), intent(in) :: funct(:)    ! function evaluated on gl integration
!                                         grid
    real(rp), intent(in) :: ds_dh_gl(:) ! path length derivative wrt height on
!                                         gl grid
    real(rp), intent(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid

! Output

    real(rp), intent(out) :: integral(:) ! result from integration

! Internals

    integer(ip) a, b, i

! Start calculation

    a = 1
    b = ng
    do i = 1, size(singularity)
      integral(i) = 0.5_rp * del_zeta(i) * sum((funct(a:b) - singularity(i))  &
               &  * ds_dh_gl(a:b) * dh_dz_gl(a:b) * gw)
      a = a + ng
      b = b + ng
    end do

  end subroutine PATH_OPACITY
!--------------------------------------------------  PATH_OPACITY  -----

  subroutine POLARIZED_PATH_OPACITY ( DEL_ZETA, DO_GL, SINGULARITY, &
                     &      FUNCT, DS_DH_GL, DH_DZ_GL, &
                     &      INTEGRAL )
    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: del_zeta(:) ! difference in integration boundary
!                                         in -log(p) units
    logical, intent(in) :: do_gl(:)  ! Where INTEGRAL needs to be evaluated
    complex(rp), intent(in) :: singularity(-1:,:) ! value of function at lower
    complex(rp), intent(in) :: funct(-1:,:)    ! function evaluated on gl integration
!                                         grid
    real(rp), intent(in) :: ds_dh_gl(:) ! path length derivative wrt height on
!                                         gl grid
    real(rp), intent(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid

! Output

    complex(rp), intent(out) :: integral(-1:,:) ! result from integration

! Internals

    integer(ip) a, b, i, j

! Start calculation

    a = 1
    b = ng
    do i = 1, size(singularity,2)
      if ( do_gl(i) ) then
        do j = -1, 1
          integral(j,i) = 0.5_rp * del_zeta(i) * sum((funct(j,a:b) - singularity(j,i))  &
                 &  * ds_dh_gl(a:b) * dh_dz_gl(a:b) * gw)
        end do 
        a = a + ng
        b = b + ng
      end if
    end do

  end subroutine POLARIZED_PATH_OPACITY


!---------------------------------------------------  HYD_OPACITY  -----

  subroutine HYD_OPACITY ( DEL_ZETA, SINGULARITY, ALPHA_PATH, H_PATH,     &
                     &     DH_DT_PATH, T_PATH, H_TAN, DH_DT_TAN, ETA_ZXP, &
                     &     DS_DH_GL, DH_DZ_GL, INTEGRAL )
    use GLNP, only: NG, GW
    use MLSCommon, only: IP, RP

! Inputs

    real(rp), intent(in) :: del_zeta(:) ! difference in integration boundary
!                                         in -log(p) units
    real(rp), intent(in) :: singularity(:) ! value of function at lower
!                                        integration boundary.
    real(rp), intent(in) :: alpha_path(:) ! absorption coefficient on gl grid.
    real(rp), intent(in) :: h_path(:)   ! path heights + req on gl grid.
    real(rp), intent(in) :: dh_dt_path(:) ! path height derivative wrt
!                                       temperature coefficient on gl grid.
    real(rp), intent(in) :: t_path(:)   ! path temperature on gl grid.
    real(rp), intent(in) :: h_tan ! tangent height + req.
    real(rp), intent(in) :: dh_dt_tan   ! height derivative wrt temperature
!                           coefficient at the tangent.
    real(rp), intent(in) :: eta_zxp(:)  ! basis function for temperature
!                                         coefficient.
    real(rp), intent(in) :: ds_dh_gl(:) ! path length derivative wrt height on
!                                         gl grid
    real(rp), intent(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid

! Output

    real(rp), intent(out) :: integral(:) ! result from integration

! Internals

    integer(ip) a, b, i

! Start calculation

    a = 1
    b = ng
    do i = 1 , size(singularity)
      integral(i) = 0.5_rp*del_zeta(i)*sum((alpha_path(a:b) - singularity(i)) &
                  * (((2.0_rp*h_path(a:b)**2 - 3.0_rp*h_tan**2) &
                  * dh_dt_path(a:b) + h_path(a:b) * h_tan * dh_dt_tan) &
                  / (sqrt(h_path(a:b)**2 - h_tan**2))**3   &
                  + eta_zxp(a:b)*ds_dh_gl(a:b)/t_path(a:b))*dh_dz_gl(a:b)*gw)
      a = a + ng
      b = b + ng
    end do

  end subroutine HYD_OPACITY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DO_DELTA_M
!---------------------------------------------------
! $Log$
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
