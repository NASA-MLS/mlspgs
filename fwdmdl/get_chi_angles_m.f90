! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_CHI_ANGLES_M

  implicit NONE

  private
  public :: GET_CHI_ANGLES

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
  ! Set up array of pointing angles

  subroutine GET_CHI_ANGLES ( sc_geoc_alt, tan_index_refr, tan_ht, &
             & phi_tan, Req, elev_offset, ptg_angle, dx_dh, dh_dz, tan_dh_dt, &
             & tan_d2h_dhdt, dx_dt, d2x_dxdt )

    use Units, only: Ln10
    use L2PC_FILE_PARAMETERS, only: DEG2RAD
    use MLSCommon, only: RP

  !  ===============================================================
  !  Declaration of variables for sub-program: get_chi_angles
  !  ===============================================================
  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------
  ! inputs

    real(rp), intent(in) :: sc_geoc_alt ! geocentric spacecraft(observer)
                                        ! altitude in km
    real(rp), intent(in) :: tan_index_refr ! tangent index of refraction
    real(rp), intent(in) :: tan_ht ! tangent height relative to Req
    real(rp), intent(in) :: phi_tan ! tangent orbit plane projected
                                    ! geodetic angle in radians
    real(rp), intent(in) :: Req ! equivalent earth radius in km
    real(rp), intent(in) :: elev_offset ! radiometer pointing offset in radians
    real(rp), intent(in) :: dh_dz ! dh/dz  at the tangent point

  ! outputs

    real(rp), intent(out) :: ptg_angle ! pointing angle in radians
    real(rp), intent(out) :: dx_dh     ! derivative of pointing angle
                                       ! wrt height

  ! keywords

    real(rp), optional, intent(in) :: tan_dh_dt(:) ! derivative of tangent
                                                   ! height wrt temperature
    real(rp), optional, intent(in) :: tan_d2h_dhdt(:) ! 2nd derivative of
                                      ! tangent height wrt temperature & height
    real(rp), optional, intent(out) :: dx_dt(:) ! derivative of pointing angle
                                                ! wrt temperature
    real(rp), optional, intent(out) :: d2x_dxdt(:) ! second derivative of
                                     ! tangent wrt temperature, pointing angle
  !  ----------------
  !  Local variables:
  !  ----------------

    real(rp), parameter :: ampl = 38.9014
    real(rp), parameter :: phas = 51.6814 * Deg2Rad

    real(rp) :: ht, tp, hs, x, q, Np1

  ! Start code:

    ht = Req + tan_ht
    Np1 = 1.0_rp + tan_index_refr
    hs = sc_geoc_alt + ampl * sin(2.0*(phi_tan-phas))
    x = Np1 * ht / hs
    ptg_angle = elev_offset + asin(x * min(ht,Req)/Req)

    if ( tan_ht >= 0.0_rp ) then
      q = Np1 - ht * tan_index_refr * Ln10 / dh_dz
    else
      q = 2.0_rp * ht * Np1 / Req
    end if
    dx_dh = q / (hs * Cos(x))

  ! Do temperature stuff if user requests it
  ! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
  ! (NOTE: These entities has NO PHI dimension, so take the center Phi in dh_dt)

    if ( present(tan_dh_dt) ) then
      tp = tan(ptg_angle)
      dx_dt = tp * tan_dh_dt / ht
      d2x_dxdt = tp * tp * tan_dh_dt / ht + tan_d2h_dhdt
    end if

    return

  end subroutine GET_CHI_ANGLES

end module GET_CHI_ANGLES_M
! $Log$
! Revision 2.10  2002/09/26 00:48:45  vsnyder
! Insert copyright notice.  Move USEs from module scope to procedure scope.
! Cosmetic changes.  Get Ln10 from Geometry module.
!
! Revision 2.9  2002/07/05 07:52:47  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.8  2002/06/28 11:06:50  zvi
! computes dx_dh on output grid as well
!
! Revision 2.7  2002/06/24 21:11:24  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.3  2002/06/11 22:20:45  bill
! eliminate zero-out feature--wgr
!
! Revision 2.2  2002/06/07 04:50:36  bill
! fixes and improvements--wgr
!
! Revision 2.1  2001/10/16 22:32:11  zvi
! Correcting for Earth intersection case
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.10.2.3  2001/09/13 11:18:21  zvi
! use the correct eta
!
! Revision 1.10.2.2  2001/09/12 21:38:47  zvi
! Added CVS stuff
!
! Revision 1.10.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
