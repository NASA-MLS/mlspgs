module SPECTRO_DERIVATIVE_M
  use GET_DRAD_M, only: GET_DRAD
  use L2PCdim, only: N2LVL
  use L2PC_PFA_STRUCTURES, only: MAXAITKENPTS, MAXPFACH, SPECTRO_PARAM
  use MLSCommon, only: I4, R4, R8
  use PFA_DB_DELTA_M, only: PFA_DB_DELTA
  implicit NONE
  private
  public :: SPECTRO_DERIVATIVE
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Subroutine SPECTRO_DERIVATIVE ( z_path, t_path, h_path, phi_path,      &
 &           DHDZ_PATH, N_lvls, jf, jch, nf, no_coeffs_f, mxco, f_basis, &
 &           mr_f, ref_corr, mnp, no_phi_f, phi_basis_f, IndxR, IndxL,   &
 &           path_brkpt, beta_t_power, pfa_dbeta_s, tau, t_script,       &
 &           Ary_Zero, ilo, ihi, spectro, Rad_s, Ier )
!
    real(r8), intent(in) :: Z_PATH(*), T_PATH(*), H_PATH(*), PHI_PATH(*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    integer(i4), intent(in) :: N_LVLS, JF, JCH, NF, NO_COEFFS_f(*), MXCO
    real(r8), intent(in) :: F_BASIS(mxco,*), MR_F(mxco,mnp,*), REF_CORR(*)
    integer(i4), intent(in) :: MNP, NO_PHI_F(*)
    real(r8), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: INDXR, INDXL, PATH_BRKPT(*)
    real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,2,*)
!   real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,maxpfach,*)
    real(r8), intent(in) :: PFA_DBETA_S(N2lvl,maxaitkenpts,2,*)
!   real(r8), intent(in) :: PFA_DBETA_S(N2lvl,maxaitkenpts,maxpfach,*)
    real(r8), intent(in) :: TAU(*), T_SCRIPT(*), ARY_ZERO(*)
    integer(i4), intent(in) :: ILO, IHI
    type(spectro_param), intent(in) :: SPECTRO
    real(r8), intent(out) :: Rad_s(maxaitkenpts,mxco,*)
    integer(i4), intent(out) :: ier
    Real(r8) ::  delta_s(N2lvl)
    Real(r8) :: zeta_basis(mxco), phi_basis(MNP)
    Integer(i4) :: s_nz,s_np,ip,iz
!
    s_np = spectro%no_phi_values
    s_nz = spectro%no_zeta_values
    if(s_nz < 1 .or. s_np < 1) Return
!
    phi_basis(1:s_np)  = dble(spectro%phi_basis(1:s_np))
    zeta_basis(1:s_nz) = dble(spectro%zeta_basis(1:s_nz))
!
    do iz = 1, s_nz
!
      do ip = 1, s_np
!
        Call pfa_db_delta ( z_path, t_path, h_path, phi_path, dHdz_path,    &
   &         N_lvls, jf, jch, nf, no_coeffs_f, mxco, N2lvl, maxaitkenpts,   &
   &         maxpfach, f_basis, mr_f, ref_corr, mnp, no_phi_f, phi_basis_f, &
   &         IndxR, IndxL, path_brkpt, beta_t_power, pfa_dbeta_s,           &
   &         zeta_basis, phi_basis, s_nz, s_np, iz, ip, delta_s, Ier )
        if (Ier /= 0) Return
!
! Now assemble the spectral derivatives for this frequency:
!
        Call get_drad ( delta_s, t_script, tau, Ary_Zero, IndxR,            &
   &                    IndxL, ilo, ihi, Rad_s(jf, iz, ip))
!
      end do                   ! ip loop
!
    end do                     ! iz loop
!
    Return
  End Subroutine SPECTRO_DERIVATIVE
end module SPECTRO_DERIVATIVE_M
! $Log$
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
