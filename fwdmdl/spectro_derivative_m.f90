! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SPECTRO_DERIVATIVE_M
  use ELLIPSE_M, only: ELLIPSE
  use GET_DRAD_NOTDER_M, only: GET_DRAD_NOTDER
  use L2PC_PFA_STRUCTURES, only: SPECTRO_PARAM
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR,PATH_DERIVATIVE
  use SPECTRO_DELTA_INTEGRAL_M, only: SPECTRO_DELTA_INTEGRAL
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
 &           k_spect,ier)
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, ILO,   &
   &                           IHI, S_NP, S_NZ, NO_ELE, frq_i

    Integer(i4), intent(in) :: no_midval_ndx,no_gl_ndx
    Integer(i4), intent(in) :: midval_ndx(:,:), gl_ndx(:,:)
!
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
    Real(r8), ALLOCATABLE, DIMENSION(:) :: Integrand
!
    phi_basis(1:s_np)  = REAL(spectro%phi_basis(1:s_np),r8)
    zeta_basis(1:s_nz) = REAL(spectro%zeta_basis(1:s_nz),r8)
!
    ALLOCATE(Integrand(no_ele), STAT=ier)
    IF(ier /= 0) THEN
      Ier = 1
      Print *,'** Error: ALLOCATION error in SPECTRO_DERIVATIVE ..'
      Return
    endif
!
!  Define the integrand array:
!
    integrand(1:no_ele) = pfa_dbeta_s(1:no_ele) * spsfunc_s(1:no_ele)
!
    do iz = 1, s_nz
!
      do ip = 1, s_np
!
        Call spectro_delta_integral(mid, brkpt, no_ele, z_path, h_path, &
          &  phi_path,dhdz_path,N_lvls,ref_corr,zeta_basis,phi_basis,   &
          &  s_nz,s_np,iz,ip,integrand,elvar,midval_ndx,no_midval_ndx,  &
          &  gl_ndx,no_gl_ndx,midval_delta,delta_s)
!
! Now assemble the spectral derivatives for this frequency:
!
        Call get_drad_notder (delta_s, t_script, tau, mid, ilo, ihi, r)
        k_spect%values(frq_i,iz,ip) = r
!
      end do                   ! ip loop
!
    end do                     ! iz loop

    DEALLOCATE(Integrand, STAT=iz)
!
    Return
  End Subroutine SPECTRO_DERIVATIVE
end module SPECTRO_DERIVATIVE_M
! $Log$
! Revision 1.11  2001/06/21 13:07:09  zvi
! Speed enhancement MAJOR update
!
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
