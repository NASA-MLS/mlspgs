module TEMPERATURE_DERIV_M
  use ELLIPSE_M, only: ELLIPSE
  use D_DELTA_DT_M, only: D_DELTA_DT
  use GET_DRAD_M, only: GET_DRAD
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA, PATH_DERIVATIVE
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
    Subroutine temperature_deriv(mid,brkpt,no_ele,t_z_basis,z_path,  &
   &           t_path,h_path,phi_path,beta_path,dHdz_path,dh_dt_path, &
   &           no_phi_t,no_t,N_lvls,n_sps,ref_corr,t_phi_basis, tau, &
   &           t_script,dt_script_dc,spsfunc_path,ilo,ihi,frq_i,     &
   &           elvar,k_temp)
!
    Integer(i4), intent(in) :: MID,BRKPT,NO_ELE,NO_PHI_T,NO_T, &
                               N_LVLS,N_SPS,ILO,IHI,FRQ_I

    Type(ELLIPSE), intent(in out) :: elvar

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: z_path, t_path, h_path, phi_path, &
   &                                 dhdz_path

    Real(r8), intent(in) :: dh_dt_path(:,:,:)
!
    Real(r8), intent(in) :: T_Z_BASIS(:), T_PHI_BASIS(:)
    Real(r8), intent(in) :: TAU(:), T_SCRIPT(:), REF_CORR(:)

    Real(r8), intent(in) :: DT_SCRIPT_DC(:,:,:)

    Type(path_derivative), INTENT(in out) :: k_temp
!
    Integer(i4) :: IN, IP
    Real(r8) :: D_DELTA_DTNP(Size(tau)), r
!
! Begin sweep through elements to see if user wants derivatives for this
! channel at this frequency
!
    do in = 1, no_t              ! Loop over the Temp. zeta coeffs.
!
      do ip = lbound(k_temp%values,3), ubound(k_temp%values,3)
        ! Loop over the Temp. phi's coeffs.
!
! Compute the temperature derivative of delta:
!
        Call d_delta_dt(mid,brkpt,no_ele,z_path,t_path,h_path,phi_path, &
       &     beta_path,dHdz_path,dh_dt_path(:,ip,in),N_lvls,n_sps,      &
       &     ref_corr,t_z_basis,no_t,t_phi_basis,no_phi_t,spsfunc_path, &
       &     in,ip,elvar,d_delta_dtnp)
!
! Now assemble the derivative:
!
        Call get_drad(d_delta_dtnp,t_script,tau,dt_script_dc(1:,in,ip), &
       &              mid,ilo,ihi,r)
        k_temp%values(frq_i,in,ip) = r
!
      end do
!
    end do
!
    Return
  End Subroutine TEMPERATURE_DERIV
end module TEMPERATURE_DERIV_M
! $Log$
! Revision 1.9  2001/03/31 23:40:56  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.8  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.7  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.6  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.5  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1  2000/06/21 21:56:17  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
