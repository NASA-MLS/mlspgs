module GET_CHI_ANGLES_M
  use MLSCommon, only: RP, IP
  use L2PC_FILE_PARAMETERS, only: DEG2RAD
  use GET_ETA_M, only: GET_ETA
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
           phi_tan,req,elev_offset,ptg_angle,tan_dh_dt,   &
           zeta_basis,phi_basis,tan_temp,tan_press,dx_dt,d2x_dxdt)

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
  REAL(rp), INTENT(IN) :: tan_ht ! tangent height relative to req
  REAL(rp), INTENT(IN) :: phi_tan ! tangent orbit plane projected
                                  ! geodetic angle in radians
  REAL(rp), INTENT(IN) :: req ! equivalent earth radius in km
  REAL(rp), INTENT(IN) :: elev_offset ! radiometer pointing offset in radians
!
! output
!
  REAL(rp), INTENT(OUT) :: ptg_angle ! pointing angle in radians
!
! keywords
!
  REAL(rp), OPTIONAL, INTENT(IN) :: tan_dh_dt(:) ! derivative of tangent
!                                                  height wrt temperature
  REAL(rp), OPTIONAL, INTENT(IN) :: zeta_basis(:) ! temperature zeta basis
  REAL(rp), OPTIONAL, INTENT(IN) :: phi_basis(:) ! temperature phi basis
  REAL(rp), OPTIONAL, INTENT(IN) :: tan_temp ! tangent temperature
  REAL(rp), OPTIONAL, INTENT(IN) :: tan_press ! tangent pressure
  REAL(rp), OPTIONAL, INTENT(OUT) :: dx_dt(:) ! derivative of pointing angle
!                                       wrt temperature
  REAL(rp), OPTIONAL, INTENT(OUT) :: d2x_dxdt(:) ! second derivative of tangent
!                                       wrt temperature, pointing angle
!  ----------------
!  Local variables:
!  ----------------
!
  Real(rp), PARAMETER :: ampl = 38.9014
  Real(rp), PARAMETER :: phas = 51.6814 * deg2rad
!
  INTEGER(ip) :: p_coeffs, z_coeffs
  REAL(rp) :: ht, tp
  REAL(rp), ALLOCATABLE :: Eta(:,:)
!
! Start code:
!
  ht = req + tan_ht
  ptg_angle = elev_offset + ASIN((1.0_rp + tan_index_refr)*ht &
            / (sc_geoc_alt + ampl * SIN(2.0*(phi_tan-phas))))
!
! do temperature stuff if user requests it
!
! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
! (NOTE: These entities has NO PHI dimension, so take the center Phi in dh_dt)
!
!  First: Get table of temperature basis functions
!
  IF(PRESENT(tan_dh_dt)) THEN
    IF(tan_ht > 0.0_rp) THEN
      p_coeffs = SIZE(phi_basis)
      z_coeffs = SIZE(zeta_basis)
      ALLOCATE(eta(1:1,1:z_coeffs))
      CALL get_eta((/tan_press/),zeta_basis,1,z_coeffs,Eta)
      tp = TAN(ptg_angle)
      dx_dt = tp * tan_dh_dt / ht
      d2x_dxdt = (2.0_rp+tp*tp)*tan_dh_dt/ht +  &
              &  RESHAPE(eta,(/z_coeffs/)) / tan_temp
      DEALLOCATE(eta)
    ELSE
      dx_dt = 0.0_rp
      d2x_dxdt = 0.0_rp
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE get_chi_angles
end module GET_CHI_ANGLES_M
! $Log$
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
