! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RAD_TRAN_M
  use GL6P, only: NG
  use MLSCommon, only: I4, R8
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE_M, only: ELLIPSE
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA
  use DO_T_SCRIPT_M, only: DO_T_SCRIPT
  use FAST_DELTA_M, only: FAST_DELTA
  use FAST_ZOPACITY_M, only: FAST_ZOPACITY
  use SCRT_DN_M, only: SCRT_DN
  implicit NONE
  private
  public :: RAD_TRAN
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This is the radiative transfer model, radiances only !

    Subroutine Rad_Tran(elvar,Frq,N_lvls,h_tan,n_sps, &
      &    ndx_path, z_path, h_path, t_path, phi_path, dHdz_path,     &
      &    earth_ref,beta_path, spsfunc_path, ref_corr, s_temp,brkpt, &
      &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
!
    Integer(i4), intent(in) :: N_LVLS, N_SPS

    Integer(i4), intent(out) :: BRKPT, NO_ELE, MID, ILO, IHI, IER

    Real(r8), intent(in) :: FRQ, H_TAN, EARTH_REF, S_TEMP

    Real(r8), intent(in) :: REF_CORR(:)

    Type(ELLIPSE), intent(in out) :: elvar

    Type(path_beta), intent(in) :: BETA_PATH(:)   ! (Nsps)

    Type(path_index), intent(in)  :: NDX_PATH
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)

    Real(r8), intent(out) :: T_SCRIPT(:), TAU(:)
    Real(r8), intent(out) :: RAD
!
    Integer(i4) :: Ngp1, i

    Real(r8) :: CSE, RS

    Real(r8) :: del_opacity(2*(n_lvls+1))
    Real(r8) :: delta(2*(n_lvls+1),n_sps)
!
!  Begin code
!
! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in ?_path%values(1...no_ele)
   
    Ier = 0
    Ngp1 = Ng + 1
    brkpt = ndx_path%break_point_index
    no_ele = ndx_path%total_number_of_elements

    elvar%EarthX = .false.
    elvar%ht = h_tan
    elvar%Rr = elvar%ht + elvar%RoC
    elvar%ht2 = elvar%Rr * elvar%Rr
    if (h_tan < -0.01) then
      elvar%Rr = elvar%Rr / elvar%RoC
      cse = earth_ref
      Call Earth_Intersection(elvar,Rs)
    else
      cse = 1.0
      elvar%Rr = 0.0d0
      elvar%Phi_s = elvar%Phi_tan
      elvar%NPhi_s = elvar%NPhi_tan
    end if
!
!  Compute the appropriate t_script & dt_scrpt_dnp along the integration
!  path of this tanget:
!
    CALL do_t_script(Ngp1, Frq, s_temp, brkpt, no_ele, t_path, &
   &                 mid, t_script)
!
    Call FAST_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path,beta_path, &
 &       dHdz_path,spsfunc_path,n_sps,N_lvls,ref_corr,elvar,delta,Ier)
    if (Ier /= 0) Return
!
! Initialize the tau & del_opacity arrays:
!
    tau(:) = 0.0
    del_opacity(:) = 0.0
!
    CALL FAST_ZOPACITY(n_sps, Ngp1, brkpt, no_ele, delta, del_opacity)
!
    Call Scrt_dn(t_script, N_lvls, cse, del_opacity, tau, Rad, mid, &
   &             ilo, ihi)
!
    Return

  End Subroutine RAD_TRAN

end module RAD_TRAN_M
! $Log$
! Revision 1.8  2001/04/09 23:33:41  zvi
! Initialize error flag
!
! Revision 1.7  2001/04/09 20:52:07  zvi
! Debugging Derivatives version
!
! Revision 1.6  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.5  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.4  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.3  2001/02/26 09:01:16  zvi
! New version - Using "Super-Structures"
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
