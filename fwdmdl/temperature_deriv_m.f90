module TEMPERATURE_DERIV_M
  use D_DELTA_DT_M, only: D_DELTA_DT
  use GET_DRAD_M, only: GET_DRAD
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PC_PFA_STRUCTURES, only: GEOPHYS_PARAM
  use L2PCdim, only: N2LVL, NLVL, NSPS
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: TEMPERATURE_DERIV

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
!
  Subroutine TEMPERATURE_DERIV ( geophysic, tndx, band, cs, Frq,       & 
 &           f_basis, z_path, t_path, h_path, phi_path, dHdz_path,     &
 &           dh_dt_path, Kgp, mnp, mdb_freq, mdb_pres, mdb_temp, mr_f, &
 &           mr_g, no_coeffs_f, no_freqs_f, no_phi_f, no_phi_t, no_t,  &
 &           N_lvls,  n_sps, phi_basis_f, ptg_i, ref_corr, sps_tbl,    &
 &           t_phi_basis,  tau, t_script, dt_script_dc, mnf, ilo, ihi, &
 &           IndxR, IndxL, path_brkpt, K_star_geophys, Ier )
!
    integer(i4), intent(in) :: KGP
    integer(i4), intent(in) :: MNF
    integer(i4), intent(in) :: NO_PHI_T

    type(geophys_param), intent(in) :: GEOPHYSIC(*)
    integer(i4), intent(in) :: TNDX
    integer(i4), intent(in) :: BAND
    real(r4), intent(inout) :: CS(Nlvl,no_phi_t,mnf,*)
    real(r8), intent(in) :: FRQ
    real(r4), intent(in) :: F_BASIS(mxco,*)
    real(r4), intent(in) :: Z_PATH(*)
    real(r4), intent(in) :: T_PATH(*)
    real(r4), intent(in) :: H_PATH(*)
    real(r4), intent(in) :: PHI_PATH(*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    real(r4), intent(in) :: DH_DT_PATH(Kgp,mnp,*)
    integer(i4), intent(in) :: MNP
    real(r8), intent(in) :: MDB_FREQ(mnf,*)
    real(r4), intent(in) :: MDB_PRES(*)
    real(r4), intent(in) :: MDB_TEMP(*)
    real(r4), intent(in) :: MR_F(mxco,mnp,*)
    real(r4), intent(in) :: MR_G(mxco,mnp,*)
    integer(i4), intent(in) :: NO_COEFFS_F(*)
    integer(i4), intent(in) :: NO_FREQS_F(*)
    integer(i4), intent(in) :: NO_PHI_F(*)
    integer(i4), intent(in) :: NO_T
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: N_SPS
    real(r4), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: PTG_I
    real(r4), intent(in) :: REF_CORR(*)
    integer(i4), intent(in) :: SPS_TBL(Nsps,*)
    real(r4), intent(in) :: T_PHI_BASIS(*)
    real(r4), intent(in) :: TAU(*)
    real(r4), intent(in) :: T_SCRIPT(*)
    real(r4), intent(in) :: DT_SCRIPT_DC(N2lvl,mxco,*)
    integer(i4), intent(in) :: ILO, IHI
    integer(i4), intent(in) :: INDXR, INDXL
    integer(i4), intent(in) :: PATH_BRKPT(*)
    real(r4), intent(out) :: K_STAR_GEOPHYS(Nlvl,mxco,mnp,*)
    integer(i4), intent(out) :: IER
!
    integer(i4) :: IN, IP
    real(r4) :: D_DELTA_DTNP(N2lvl)
!
! Begin sweep through elements to see if user wants derivatives for this
! channel at this frequency
!
    do in = 1, no_t              ! Loop over the Temp. zeta coeffs.
!
      do ip = 1, no_phi_t        ! Loop over the Temp. phi's coeffs.
!
! Compute the temperature derivative of delta:
!
        Call d_delta_dt(z_path,t_path,h_path,phi_path,dHdz_path,       &
   &         dh_dt_path,Kgp,N_lvls,cs,Frq,n_sps,no_coeffs_f,           &
   &         sps_tbl(1,band),mxco,Nlvl,f_basis,                        &
   &         ref_corr,geophysic(tndx)%basis_peaks,mr_g(1,1,tndx),      &
   &         no_t,t_phi_basis,no_phi_t,mdb_pres,mdb_temp,mdb_freq,     &
   &         mnp,mnf,no_freqs_f,no_phi_f,phi_basis_f,mr_f,in,ip,       &
   &         IndxR,IndxL,path_brkpt,d_delta_dtnp,Ier)
        if (Ier /= 0) Return
!
! Now assemble the derivative:
!
          Call get_drad(d_delta_dtnp,t_script,tau,                     &
     &                  dt_script_dc(1,in,ip),IndxR,IndxL,ilo,ihi,     &
     &                  K_star_geophys(ptg_i,in,ip,tndx))
!
        end do
!
      end do
!
    Return
  End Subroutine TEMPERATURE_DERIV
end module TEMPERATURE_DERIV_M

! $Log$
