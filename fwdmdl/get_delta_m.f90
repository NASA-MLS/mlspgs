module GET_DELTA_M
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  use GENERIC_DELTA_INTEGRAL_M, only: GENERIC_DELTA_INTEGRAL
  implicit NONE
  private
  public :: GET_DELTA
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrate in ZETA Space !
!
  Subroutine GET_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path,         &
 &           beta_path,dHdz_path,n_sps,N_lvls,Nc,ncoeffs,Nlvl,  &
 &           z_basis,ref_corr,mnp,no_phi_f,phi_basis,spsfunc_path,mr_f, &
 &           is_f_log,delta,Ier)

    Logical,     intent(in) :: IS_F_LOG(*)
    Integer(i4), intent(in) :: NLVL, NC, N_SPS, N_LVLS, MNP
    Integer(i4), intent(in) :: NCOEFFS(*), NO_PHI_F(*)
    Integer(i4), intent(in) :: mid, brkpt, no_ele

    Real(r8), intent(in) :: REF_CORR(*)
    Real(r8), intent(in) :: Z_BASIS(Nc,*)
    Real(r8), intent(in) :: MR_F(Nc,mnp,*)
    Real(r8), intent(in) :: PHI_BASIS(mnp,*)

    Real(r8), intent(inout) :: DELTA(2*Nlvl,Nc,mnp,*)

    Integer(i4), intent(out) :: IER

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, DHDZ_PATH
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: J, K, IP, IZ, NCO, NPF

    Real(r8) :: Q
!
    Real(r8), ALLOCATABLE, DIMENSION(:) :: Integrand
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
    DEALLOCATE(Integrand, STAT=k)
!
    ALLOCATE(Integrand(no_ele), STAT=ier)
    IF(ier /= 0) THEN
      Ier = 1
      Print *,'** Error: ALLOCATION error in GET_DELTA routine ..'
      goto 99
    endif
!
!  Initialize all arrays:
!
    k = 2 * N_lvls
    do j = 1, n_sps
      delta(1:k,1:ncoeffs(j),1:no_phi_f(j),j) = 0.0
    end do
!
    do j = 1, n_sps
!
      npf = no_phi_f(j)
      nco = ncoeffs(j)
!
      integrand(1:no_ele) = beta_path(j)%values(1:no_ele)
!
      if(is_f_log(j)) then
        integrand(1:no_ele) = integrand(1:no_ele) *  &
       &                      spsfunc_path(j)%values(1:no_ele)
      endif
!
! Loop over the specie's Phi's
!
      do ip = 1, npf
!
! Loop over the specie's zeta's
!
        do iz = 1, nco
!
          q = 1.0
          if(is_f_log(j)) q = 1.0 / mr_f(iz,ip,j)
!
          Call generic_delta_integral(mid, brkpt, no_ele, z_path,   &
         &     h_path, phi_path, dhdz_path, N_lvls, ref_corr,       &
         &     integrand, z_basis(1:,j), phi_basis(1:,j), nco, npf, &
         &     iz, ip, q, delta(1:,iz,ip,j), Ier)
          IF(ier /= 0) goto 99
!
        end do
!
      end do
!
    end do
!
 99  DEALLOCATE(Integrand, STAT=k)
!
    Return
!
  End Subroutine GET_DELTA
!
end module GET_DELTA_M
! $Log$
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
