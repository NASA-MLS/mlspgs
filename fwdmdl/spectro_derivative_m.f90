module SPECTRO_DERIVATIVE_M
  use ELLIPSE_M, only: ELLIPSE
  use GET_DRAD_M, only: GET_DRAD
  use L2PCdim, only: N2LVL
  use L2PC_PFA_STRUCTURES, only: SPECTRO_PARAM
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_DERIVATIVE
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
  Subroutine SPECTRO_DERIVATIVE(mid,brkpt,no_ele,z_path,h_path,phi_path, &
 &           DHDZ_PATH,N_lvls,mxco,ref_corr,mnp,spsfunc_s,pfa_dbeta_s,   &
 &           tau,t_script,Ary_Zero,s_np,s_nz,ilo,ihi,spectro,frq_i,elvar,&
 &           k_spect,ier)
!
    Integer(i4), intent(in) :: MNP, N_LVLS, MXCO, MID, BRKPT, ILO,   &
   &                           IHI, S_NP, S_NZ, NO_ELE, frq_i
!
    Type(path_vector), intent(in) :: z_path, h_path, phi_path, dhdz_path

    Type(ELLIPSE), INTENT(IN OUT) :: elvar

    Real(r8), intent(in) :: REF_CORR(:), TAU(:)
    Real(r8), intent(in) :: T_SCRIPT(:), ARY_ZERO(:), PFA_DBETA_S(:), &
   &                        SPSFUNC_S(:)
!
    Type(spectro_param), intent(in) :: SPECTRO
!
    Type(path_derivative), intent(in out) :: k_spect
!
    Integer(i4), intent(out) :: ier
!
    Integer(i4) :: ip, iz
    Real(r8)    :: r, delta_s(N2lvl), zeta_basis(mxco), phi_basis(MNP)
!
    phi_basis(1:s_np)  = REAL(spectro%phi_basis(1:s_np),r8)
    zeta_basis(1:s_nz) = REAL(spectro%zeta_basis(1:s_nz),r8)
!
    do iz = 1, s_nz
!
      do ip = 1, s_np
!
        Call pfa_db_delta (mid, brkpt, no_ele, z_path, h_path, phi_path, &
 &           DHDZ_PATH, N_lvls, ref_corr, spsfunc_s, pfa_dbeta_s,        &
 &           zeta_basis,phi_basis,s_nz,s_np,iz,ip,elvar,delta_s,Ier )
        if (Ier /= 0) Return
!
! Now assemble the spectral derivatives for this frequency:
!
        Call get_drad (delta_s, t_script, tau, Ary_Zero, mid, ilo, ihi, r)
        k_spect%values(frq_i,iz,ip) = r
!
      end do                   ! ip loop
!
    end do                     ! iz loop
!
    Return
  End Subroutine SPECTRO_DERIVATIVE
end module SPECTRO_DERIVATIVE_M
! $Log$
! Revision 1.6  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.5  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1  2000/06/21 21:56:17  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
