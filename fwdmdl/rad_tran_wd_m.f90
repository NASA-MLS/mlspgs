module RAD_TRAN_WD_M
  use GL6P, only: NG
  use ELLIPSE_M, only: ELLIPSE
  use MLSCommon, only: I4, R8
  use ForwardModelConfig, only: ForwardModelConfig_T
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use L2PC_PFA_STRUCTURES, only: SPECTRO_PARAM
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

Subroutine Rad_Tran_WD(ForwardModelConfig, FwdModelExtra, FwdModelIn, &
       &   elvar,frq_i,band,Frq,n_sps,z_path,h_path,t_path,phi_path,  &
       &   dHdz_path,beta_path,spsfunc_path,t_z_basis,no_t,ref_corr,  &
       &   no_phi_t,t_phi_basis,dh_dt_path,spect_atmos,spectroscopic, &
       &   k_temp,k_atmos,k_spect_dw,k_spect_dn,k_spect_dnu,brkpt,    &
       &   no_ele,mid,ilo,ihi,t_script,tau,max_zeta_dim,max_phi_dim,Ier)
!
    type(forwardModelConfig_T), intent(in) :: forwardModelConfig
    type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra
!
    Type(ELLIPSE), INTENT(IN OUT) :: elvar

    Integer(i4), intent(in) :: FRQ_I,NO_PHI_T,NO_T, BAND, &
   &                           N_SPS, BRKPT, NO_ELE, MID, ILO, IHI, &
   &                           max_zeta_dim, max_phi_dim

    Integer(i4), intent(in) :: SPECT_ATMOS(:)

    Integer(i4), intent(out) :: IER

    Real(r8), intent(in) :: FRQ

    Real(r8), intent(in) :: T_Z_BASIS(:), T_PHI_BASIS(:)

    Real(r8), intent(in) :: TAU(:)
    Real(r8), intent(in) :: T_SCRIPT(:)
    Real(r8), intent(in) :: REF_CORR(:)

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

    Integer(i4) :: i, j, Ngp1, Spectag, N_lvls
    Integer(i4) :: nf, sa, s_np, s_nz

    Real(r8) :: dt_scrpt_dnp(size(tau),no_t,no_phi_t)
    Real(r8) :: delta(size(tau),max_zeta_dim,max_phi_dim,n_sps)
!
!  Begin code:
!
    ier = 0
    Ngp1 = Ng + 1
    N_lvls = ForwardModelConfig%integrationGrid%noSurfs

    if(forwardModelConfig%atmos_der) then
!
!  Atmospheric derivatives:
!
      Call GET_DELTA(ForwardModelConfig, FwdModelExtra, FwdModelIn, &
   &       mid,brkpt,no_ele,z_path,h_path,phi_path,beta_path,dHdz_path, &
   &       n_sps, N_lvls, ref_corr,spsfunc_path,elvar,delta,Ier)
      if (Ier /= 0) Return
!
! Compute atmosperic derivatives for this channel
!
      Call ZATMOS_DERIV(ForwardModelConfig, FwdModelExtra, FwdModelIn, &
     &       frq_i, mid, delta, t_script, tau, ilo, ihi, k_atmos,Ier)
      if (Ier /= 0) Return
!
    endif
!
! Compute temperature derivative for this channel (if requested)
!
    if (ForwardModelConfig%temp_der) then
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
     &     dt_scrpt_dnp,spsfunc_path,ilo,ihi,frq_i,elvar,k_temp)
!
    end if
!
    if (ForwardModelConfig%spect_der) then
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
 &               h_path, phi_path, DHDZ_PATH, N_lvls,ref_corr,    &
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dw, tau,  &
 &               t_script,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,elvar,k_spect_dw(nf), Ier )
!
          else if (CA == 'N') then
            Call spectro_derivative(mid, brkpt, no_ele, z_path,         &
 &               h_path, phi_path, DHDZ_PATH, N_lvls,ref_corr,    &
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dn, tau,  &
 &               t_script,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,elvar,k_spect_dn(nf), Ier )
!
          else if (CA == 'V') then
            Call spectro_derivative(mid, brkpt, no_ele, z_path,         &
 &               h_path, phi_path, DHDZ_PATH, N_lvls,ref_corr,    &
 &               spsfunc_path(nf)%values, beta_path(nf)%dbeta_dnu, tau, &
 &               t_script,s_np,s_nz,ilo,ihi,spectroscopic(sa), &
 &               frq_i,elvar,k_spect_dnu(nf), Ier )
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
! Revision 1.9  2001/04/09 20:52:07  zvi
! Debugging Derivatives version
!
! Revision 1.8  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.7  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.6  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.5  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.4  2001/03/08 00:11:16  zvi
! *** empty log message ***
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
