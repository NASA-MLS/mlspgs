module RAD_TRAN_WD_M
  use GL6P, only: NG
  use L2PCDim, only: NLVL, NSPS, N2LVL, MNP => max_no_phi
  use MLSCommon, only: I4, R8
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, SPECTRO_PARAM
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA, PATH_DERIVATIVE
  use D_T_SCRIPT_DTNP_M, only: D_T_SCRIPT_DTNP
  use GET_DELTA_M, only: GET_DELTA
  use SPECTRO_DERIVATIVE_M, only: SPECTRO_DERIVATIVE
  use TEMPERATURE_DERIV_M, only: TEMPERATURE_DERIV
  use ZATMOS_DERIV_M, only: ZATMOS_DERIV
  implicit NONE
  private
  public :: RAD_TRAN_WD
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This is the radiative transfer model with derivatives

Subroutine Rad_Tran_WD(frq_i,band,Frq,N_lvls,n_sps,temp_der,atmos_der,&
      &    spect_der,z_path,h_path,t_path,phi_path,dHdz_path,atmospheric,&
      &    beta_path,spsfunc_path,t_z_basis,f_basis,no_coeffs_f,mr_f, &
      &    no_t, ref_corr, no_phi_f, phi_basis_f, no_phi_t, &
      &    t_phi_basis, dh_dt_path, spect_atmos, spectroscopic, k_temp,&
      &    k_atmos, k_spect_dw, k_spect_dn, k_spect_dnu,is_f_log,brkpt,&
      &    no_ele, mid, ilo, ihi, t_script, tau, Ier)
!
    Logical, intent(in) :: temp_der,atmos_der,spect_der
    Logical, intent(in) :: IS_F_LOG(:)

    Integer(i4), intent(in) :: FRQ_I,N_LVLS,NO_PHI_T,NO_T, BAND, &
   &                           N_SPS, BRKPT, NO_ELE, MID, ILO, IHI

    Integer(i4), intent(in) :: NO_COEFFS_F(:), NO_PHI_F(:), SPECT_ATMOS(:)

    Integer(i4), intent(out) :: IER

    Real(r8), intent(in) :: FRQ

    Real(r8), intent(in) :: T_Z_BASIS(:), T_PHI_BASIS(:)

    Real(r8), intent(in) :: TAU(:)
    Real(r8), intent(in) :: T_SCRIPT(:)
    Real(r8), intent(in) :: REF_CORR(:)
    Real(r8), intent(in) :: F_BASIS(:,:)
    Real(r8), intent(in) :: PHI_BASIS_F(:,:)
    Real(r8), intent(in) :: MR_F(:,:,:)

    Type(atmos_comp), intent(in) :: ATMOSPHERIC(:)

    Type(spectro_param), intent(in) :: SPECTROSCOPIC(:)

    Type(path_beta), intent(in) :: BETA_PATH(:)   ! (Nsps)

    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH

    Real(r8), intent(in) :: dh_dt_path(:,:,:)

    Type(path_derivative), intent(in out) :: k_temp
    Type(path_derivative), intent(in out) :: k_atmos(:), k_spect_dw(:), &
                                             k_spect_dn(:), k_spect_dnu(:)
!
    CHARACTER (LEN=01) :: CA

    Integer(i4) :: i, j, Ngp1, Spectag
    Integer(i4) :: nf, sa, s_np, s_nz

    Real(r8) :: dt_scrpt_dnp(N2lvl,mxco,mnp), delta(N2lvl,mxco,mnp,Nsps)

    Real(r8), Parameter, Dimension(N2lvl) :: ary_zero = 0.0
!
    Ngp1 = Ng + 1

    if(atmos_der) then
!
!  Atmospheric derivatives:
!
      Call GET_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path,   &
   &       beta_path,dHdz_path,n_sps,N_lvls,mxco,no_coeffs_f,   &
   &       Nlvl,f_basis,ref_corr,mnp,no_phi_f, &
   &       phi_basis_f,spsfunc_path,mr_f,is_f_log,delta,Ier)
      if (Ier /= 0) Return
!
! Compute atmosperic derivatives for this channel
!
      Call zatmos_deriv(atmospheric,n_sps,band,frq_i, &
   &              no_coeffs_f,mid,delta,t_script,tau,ilo,ihi, &
   &              no_phi_f,k_atmos,Ier)
      if (Ier /= 0) Return
!
    endif
!
! Compute temperature derivative for this channel (if requested)
!
    if (temp_der) then
!
! Create the dt_scrpt_dnp arrays for all coefficients:
!
      do i = 1, no_phi_t
        do j = 1, no_t
          CALL d_t_script_dtnp(Frq, t_z_basis, t_phi_basis,    &
         &     brkpt, no_ele, z_path, t_path, phi_path,        &
         &     Ng, j, i, no_t, no_phi_t, dt_scrpt_dnp(1:,j,i))
        end do
      end do
!
      Call temperature_deriv(mid,brkpt,no_ele,t_z_basis,z_path,t_path,   &
     &     h_path,phi_path,beta_path,dHdz_path,dh_dt_path,no_phi_t,no_t, &
     &     N_lvls,n_sps,ref_corr,t_phi_basis,tau,t_script, &
     &     dt_scrpt_dnp,spsfunc_path,ilo,ihi,frq_i,k_temp)
!
    end if
!
    if (spect_der) then
!
!  ** Spectroscopic derivatives here:
!
!  Get the dI/d(Spectral parameters) (w, n & nu0) FOR EACH SPECIE,
!  and for each Zeta/Phi coefficient.
!
      do nf = 1, n_sps
!
        sa = spect_atmos(nf)
        if (sa < 1) CYCLE
!
        if(.not.spectroscopic(sa)%DER_CALC(band)) CYCLE
!
        Spectag = spectroscopic(sa)%spectag
!
        DO
!
          if(spectroscopic(sa)%Spectag /= Spectag) EXIT
!
          s_np = spectroscopic(sa)%no_phi_values
          s_nz = spectroscopic(sa)%no_zeta_values
          if(s_nz < 1 .or. s_np < 1) EXIT
!
          CA = spectroscopic(sa)%type

          if (CA == 'W') then
            Call spectro_derivative(mid, brkpt, no_ele, z_path,         &
 &               h_path, phi_path, DHDZ_PATH, N_lvls,mxco, ref_corr,mnp,&
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dw, tau,  &
 &               t_script,ary_zero,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,k_spect_dw(nf), Ier )
!
          else if (CA == 'N') then
            Call spectro_derivative(mid, brkpt, no_ele, z_path,         &
 &               h_path, phi_path, DHDZ_PATH, N_lvls,mxco, ref_corr,mnp,&
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dn, tau,  &
 &               t_script,ary_zero,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,k_spect_dn(nf), Ier )
!
          else if (CA == 'V') then
            Call spectro_derivative(mid, brkpt, no_ele, z_path,         &
 &               h_path, phi_path, DHDZ_PATH, N_lvls,mxco, ref_corr,mnp,&
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dnu, tau, &
 &               t_script,ary_zero,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,k_spect_dnu(nf), Ier )
!
          end if
!
          if (Ier /= 0) Return
!
          sa = sa + 1
          if(sa > 3 * n_sps) EXIT
!
        END DO
!
      end do                      ! On nf (Specie loop)
!
    end if                        ! on Derivatives 'if'
!
    Return
  End Subroutine RAD_TRAN_WD
end module RAD_TRAN_WD_M
! $Log$
! Revision 1.5  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.4  2001/03/08 00:11:16  zvi
! *** empty log message ***
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
