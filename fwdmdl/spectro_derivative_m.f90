! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SPECTRO_DERIVATIVE_M
  use ELLIPSE_M, only: ELLIPSE
  use GET_DRAD_NOTDER_M, only: GET_DRAD_NOTDER
  use L2PC_PFA_STRUCTURES, only: SPECTRO_PARAM
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR,PATH_DERIVATIVE,PATH_INT_VECTOR_2D
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
 &           DHDZ_PATH,N_lvls,ref_corr,spsfunc_s,pfa_dbeta_s,    &
 &           tau,t_script,s_np,s_nz,ilo,ihi,spectro,frq_i,elvar, &
 &           midval_ndx,no_midval_ndx,gl_ndx,no_gl_ndx,midval_delta, &
 &           Sps_zeta_loop, Sps_phi_loop, k_spect,ier)
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, ILO,   &
   &                           IHI, S_NP, S_NZ, NO_ELE, frq_i

    Integer(i4), intent(in) :: no_midval_ndx,no_gl_ndx
    Integer(i4), intent(in) :: midval_ndx(:,:), gl_ndx(:,:)
!
    Type(path_int_vector_2d), intent(in) :: Sps_zeta_loop,Sps_phi_loop
    Type(path_vector), intent(in) :: z_path, h_path, phi_path, dhdz_path

    Type(ELLIPSE), INTENT(IN OUT) :: elvar

    Real(r8), intent(in) :: midval_delta(:)    ! (N2lvl)

    Real(r8), intent(in) :: REF_CORR(:), TAU(:)
    Real(r8), intent(in) :: T_SCRIPT(:), PFA_DBETA_S(:), SPSFUNC_S(:)
!
    Type(spectro_param), intent(in) :: SPECTRO
!
    Type(path_derivative), intent(in out) :: k_spect
!
    Integer(i4), intent(out) :: ier
!
    Integer(i4) :: ip, iz
    Real(r8) :: r, delta_s(2*(N_lvls+1))
    Real(r8) :: zeta_basis(s_nz), phi_basis(s_np)
!
    phi_basis(1:s_np)  = REAL(spectro%phi_basis(1:s_np),r8)
    zeta_basis(1:s_nz) = REAL(spectro%zeta_basis(1:s_nz),r8)
!
    do iz = 1, s_nz
!
      do ip = 1, s_np
!
        Call PFA_DB_DELTA (mid, brkpt, no_ele, z_path, h_path, phi_path, &
 &           dhdz_path, N_lvls, ref_corr, spsfunc_s, pfa_dbeta_s,        &
 &           zeta_basis, phi_basis, s_nz, s_np, iz, ip, elvar,           &
 &           midval_delta, midval_ndx, no_midval_ndx, gl_ndx, no_gl_ndx, &
 &           Sps_zeta_loop, Sps_phi_loop, delta_s, Ier)
        if (Ier /= 0) Return
!
! Now assemble the spectral derivatives for this frequency:
!
        Call get_drad_notder (delta_s, t_script, tau, mid, ilo, ihi, r)
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
! Revision 1.10  2001/06/07 23:39:32  pwagner
! Added Copyright statement
!
! Revision 1.9  2001/04/07 23:51:33  zvi
! *** empty log message ***
!
! Revision 1.8  2001/03/31 23:40:56  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.7  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
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
