module DO_DELTA_M
  USE GLNP, ONLY: NG, GW
  use MLSCommon, only: IP, RP
  implicit NONE
  private
  PUBLIC :: PATH_OPACITY, HYD_OPACITY
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------------
!
 SUBROUTINE path_opacity(del_zeta,singularity,funct,ds_dh_gl,dh_dz_gl, &
                     &   integral)
!
! Inputs
!
  REAL(rp), INTENT(in) :: del_zeta(:) ! difference in integration boundary
!                                       in -log(p) units
  REAL(rp), INTENT(in) :: singularity(:) ! value of function at lower
  REAL(rp), INTENT(in) :: funct(:)    ! function evaluated on gl integration
!                                       grid
  REAL(rp), INTENT(in) :: ds_dh_gl(:) ! path length derivative wrt height on
!                                       gl grid
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid
!
! Output
!
  REAL(rp), INTENT(out) :: integral(:) ! result from integration
! Internals
!
  INTEGER(ip) a,b,i
!
! Start calculation
!
  a = 1
  b = ng
  DO i = 1 , SIZE(singularity)
    integral(i) = 0.5_rp * del_zeta(i) * SUM((funct(a:b) - singularity(i))  &
             &  * ds_dh_gl(a:b) * dh_dz_gl(a:b) * Gw)
    a = a + ng
    b = b + ng
  ENDDO
!
 END SUBROUTINE path_opacity
!--------------------------------------------------------------------
!
 SUBROUTINE hyd_opacity(del_zeta,singularity,alpha_path,h_path,dh_dt_path, &
                     &  t_path,h_tan,dh_dt_tan,eta_zxp,ds_dh_gl,dh_dz_gl, &
                     &  integral)
!
! Inputs
!
  REAL(rp), INTENT(in) :: del_zeta(:) ! difference in integration boundary
!                                       in -log(p) units
  REAL(rp), INTENT(in) :: singularity(:) ! value of function at lower
!                                      integration boundary.
  REAL(rp), INTENT(in) :: alpha_path(:) ! absorption coefficient on gl grid.
  REAL(rp), INTENT(in) :: h_path(:) ! path heights + req on gl grid.
  REAL(rp), INTENT(in) :: dh_dt_path(:) ! path height derivative wrt
!                                     temperature coefficient on gl grid.
  REAL(rp), INTENT(in) :: t_path(:) ! path temperature on gl grid.
  REAL(rp), INTENT(in) :: h_tan ! tangent height + req.
  REAL(rp), INTENT(in) :: dh_dt_tan ! height derivative wrt temperature
!                         coefficient at the tangent.
  REAL(rp), INTENT(in) :: eta_zxp(:) ! basis function for temperature
!                                      coefficient.
  REAL(rp), INTENT(in) :: ds_dh_gl(:) ! path length derivative wrt height on
!                                       gl grid
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! height derivative wrt zeta on gl grid
!
! Output
!
  REAL(rp), INTENT(out) :: integral(:) ! result from integration
!
! Internals
!
  INTEGER(ip) a,b,i
!
! Start calculation
!
  a = 1
  b = ng
  DO i = 1 , SIZE(singularity)
    integral(i) = 0.5_rp*del_zeta(i)*SUM((alpha_path(a:b) - singularity(i)) &
                * (((2.0_rp*h_path(a:b)**2 - 3.0_rp*h_tan**2) &
                * dh_dt_path(a:b) + h_path(a:b) * h_tan * dh_dt_tan) &
                / (SQRT(h_path(a:b)**2 - h_tan**2))**3   &
                + eta_zxp(a:b)*ds_dh_gl(a:b)/t_path(a:b))*dh_dz_gl(a:b)*Gw)
    a = a + ng
    b = b + ng
  ENDDO
!
 END SUBROUTINE hyd_opacity
!
END module DO_DELTA_M
!---------------------------------------------------
! $Log$
! Revision 1.1.2.2  2001/09/12 21:38:45  zvi
! Added CVS stuff
!
!
