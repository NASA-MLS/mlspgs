module SCRT_DN_M
  use MLSCommon, only: IP, RP
  implicit NONE
  private
  PUBLIC :: SCRT_DN, GET_DSCRT_NO_T_DN, GET_DSCRT_DN
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
!
 SUBROUTINE scrt_dn(t_script,e_rflty,incoptdepth,tau,radiance,i_stop)
!
! inputs:
!
  REAL(rp), INTENT(in) :: t_script(:) ! differential temperatures (K).
  REAL(rp), INTENT(in) :: e_rflty ! Earth surface reflectivity (0--1).
  REAL(rp), INTENT(in) :: incoptdepth(:) ! layer incremental optical depth,
!                         this comes in as a positive number.
! outputs:
!
  REAL(rp), INTENT(out) :: tau(:) ! transmission function
  REAL(rp), INTENT(out) :: radiance    ! radiance(K).
  INTEGER(ip), INTENT(out) :: i_stop ! integration stop index
!
! internals
!
  INTEGER(ip) :: n_path, half_path
  REAL(rp) :: total_opacity
  REAL(rp), PARAMETER :: black_out = -15.0_rp
!
! begin code
!
  n_path = SIZE(t_script)
  half_path = n_path/2
  tau(1) = 1.0_rp
  total_opacity = 0.0_rp
  radiance = t_script(1)
  i_stop = 2
!
  DO
    IF(total_opacity < black_out .or. i_stop > half_path) EXIT
    total_opacity = total_opacity - incoptdepth(i_stop)
    tau(i_stop) = EXP(total_opacity)
    radiance = radiance + t_script(i_stop) * tau(i_stop)
    i_stop = i_stop + 1
  ENDDO
!
  IF(i_stop <= half_path) THEN
    tau(i_stop:n_path) = 0.0_rp
    i_stop = i_stop - 1
    RETURN
  ENDIF
!
! the tangent
!
  tau(i_stop) = e_rflty * tau(i_stop-1)
  radiance = radiance + t_script(i_stop) * tau(i_stop)
!
  DO
    IF(total_opacity < black_out .or. i_stop == n_path) EXIT
    total_opacity = total_opacity - incoptdepth(i_stop)
    i_stop = i_stop + 1
    tau(i_stop) = EXP(total_opacity)
    radiance = radiance + t_script(i_stop) * tau(i_stop)
  ENDDO
!
  IF(i_stop < n_path) THEN
    tau(i_stop:n_path) = 0.0_rp
    RETURN
  ENDIF
!
 END SUBROUTINE scrt_dn
!-------------------------------------------------------------------------
! This subroutine computes the scalarized condensed radiative transfer
! derivatives, without Temperature derivatives
!
  SUBROUTINE get_dscrt_no_t_dn(d_delta_dx,t_script,tau,i_start,i_end,drad_dx)
!
! Inputs
!
  REAL(rp), INTENT(in) :: d_delta_dx(:) ! path opacity derivatives
  REAL(rp), INTENT(in) :: t_script(:) ! differential temparatures
  REAL(rp), INTENT(in) :: tau(:) ! path transmission path
!
  INTEGER(ip), INTENT(in) :: i_start ! where non-zeros on the path begin
  INTEGER(ip), INTENT(in) :: i_end   ! where non-zeros on path end
!
! Output
!
  REAL(rp), INTENT(out) :: drad_dx ! radiance derivative wrt x
!
! internals
!
  Integer(ip) :: i, n_path, half_path
  REAL(rp) :: w
!
  w = 0.0_rp
  drad_dx = 0.0_rp
  n_path = SIZE(t_script)
  half_path = n_path/2
!
  DO i = MAX(2,i_start),MIN(i_end,half_path)
    w = w - d_delta_dx(i)
    drad_dx = drad_dx + t_script(i) * w * tau(i)
  ENDDO
!
  IF(i_end <= half_path) RETURN
!
  drad_dx = drad_dx + t_script(half_path+1) * w * tau(half_path)
!
  DO i = half_path+2 , MIN(n_path,i_end)
    w = w - d_delta_dx(i-1)
    drad_dx = drad_dx + t_script(i) * w * tau(i)
  END DO
!
  RETURN
!
 END  Subroutine get_dscrt_no_t_dn
!---------------------------------------------------------------------------
! This subroutine computes the scalarized condensed radiative transfer
! derivatives
!
  Subroutine get_dscrt_dn(DDER, T_SCRIPT, TAU, DT_SCRIPT_DX, I_start, &
                  &       i_end, drad_dx)
!
    Integer(ip), intent(in) :: i_start, i_end
    Real(rp), intent(in) :: DDER(:), T_SCRIPT(:), TAU(:), DT_SCRIPT_DX(:)
    Real(rp), intent(out) :: drad_dx

    Integer(ip) :: i, n_path, half_path
    Real(rp) :: w
!
    w = 0.0_rp
    drad_dx = dt_script_dx(1)

    n_path = SIZE(t_script)
    half_path = n_path/2
!
    DO i = MAX(2,i_start),MIN(i_end,half_path)
      w = w - dder(i)
      drad_dx = drad_dx + (dt_script_dx(i) + t_script(i) * w) * tau(i)
    end do
!
    IF(i_end <= half_path) RETURN

    drad_dx = drad_dx + (dt_script_dx(half_path+1)  &
            + t_script(half_path+1) * w) * tau(half_path)
!
    DO i = half_path+2 , MIN(n_path,i_end)
      w = w - dder(i-1)
      drad_dx = drad_dx + (dt_script_dx(i) + t_script(i) * w) * tau(i)
    end do
!
    Return
!
  END Subroutine get_dscrt_dn
!
END module SCRT_DN_M
! $Log$
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.6.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
