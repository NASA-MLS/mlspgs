module GET_CHI_ANGLES_M
  use MLSCommon, only: RP, IP
  use L2PC_FILE_PARAMETERS, only: DEG2RAD
  Implicit NONE
  Private
  Public :: get_chi_angles
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! Set up array of pointing angles

SUBROUTINE get_chi_angles(sc_geoc_alt,tan_index_refr,tan_ht, &
           phi_tan,Req,elev_offset,ptg_angle,dx_dh,dh_dz,tan_dh_dt,&
           tan_d2h_dhdt,dx_dt,d2x_dxdt)

!  ===============================================================
!  Declaration of variables for sub-program: get_chi_angles
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
! inputs
!
  REAL(rp), INTENT(IN) :: sc_geoc_alt ! geocentric spacecraft(observer)
                                      ! altitude in km
  REAL(rp), INTENT(IN) :: tan_index_refr ! tangent index of refraction
  REAL(rp), INTENT(IN) :: tan_ht ! tangent height relative to Req
  REAL(rp), INTENT(IN) :: phi_tan ! tangent orbit plane projected
                                  ! geodetic angle in radians
  REAL(rp), INTENT(IN) :: Req ! equivalent earth radius in km
  REAL(rp), INTENT(IN) :: elev_offset ! radiometer pointing offset in radians
  REAL(rp), INTENT(IN) :: dh_dz ! dh/dz  at the tangent point
!
! output
!
  REAL(rp), INTENT(OUT) :: ptg_angle ! pointing angle in radians
  REAL(rp), INTENT(OUT) :: dx_dh     ! derivative of pointing angle
                                     ! wrt height
!
! keywords
!
  REAL(rp), OPTIONAL, INTENT(IN) :: tan_dh_dt(:,:) ! derivative of tangent
!                                                  height wrt temperature
  REAL(rp), OPTIONAL, INTENT(IN) :: tan_d2h_dhdt(:,:) ! 2nd derivative of
! tangent height wrt temperature & height
  REAL(rp), OPTIONAL, INTENT(OUT) :: dx_dt(:,:) ! derivative of pointing angle
!                                       wrt temperature
  REAL(rp), OPTIONAL, INTENT(OUT) :: d2x_dxdt(:,:) ! second derivative of
! tangent wrt temperature, pointing angle
!  ----------------
!  Local variables:
!  ----------------
!
  Real(rp), PARAMETER :: ampl = 38.9014
  Real(rp), PARAMETER :: phas = 51.6814 * Deg2Rad
  Real(rp), PARAMETER :: Ln10 = 2.302585092994045684_rp
!
  REAL(rp) :: ht, tp, hs, x, q, Np1
!
! Start code:
!
  ht = Req + tan_ht
  Np1 = 1.0_rp + tan_index_refr
  hs = sc_geoc_alt + ampl * SIN(2.0*(phi_tan-phas))
  x = Np1 * ht / hs
  ptg_angle = elev_offset + ASIN(x * min(ht,Req)/Req)
!
  if(tan_ht >= 0.0_rp) then
    q = Np1 - ht * tan_index_refr * Ln10 / dh_dz
  else
    q = 2.0_rp * ht * Np1 / Req
  endif
  dx_dh = q / (hs * Cos(x))
!
! Do temperature stuff if user requests it
! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
! (NOTE: These entities has NO PHI dimension, so take the center Phi in dh_dt)
!
  IF(PRESENT(tan_dh_dt)) THEN
    tp = TAN(ptg_angle)
    dx_dt = tp * tan_dh_dt / ht
    d2x_dxdt = tp*tp*tan_dh_dt/ht + tan_d2h_dhdt
  ENDIF

  RETURN

END SUBROUTINE get_chi_angles

end module GET_CHI_ANGLES_M
! $Log$
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
