module RAD_TRAN_M
  use GL6P, only: NG
  use L2PCDim, only: NLVL, NSPS, N2LVL
  use MLSCommon, only: I4, R8
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE, only: EARTHX, HT, HT2, NPHI_TAN, NPHI_S, PHI_S, PHI_TAN, &
                     ROC, RR
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

Subroutine Rad_Tran(Frq,N_lvls,h_tan,band,n_sps, sps_tbl, &
      &    ndx_path, z_path, h_path, t_path, phi_path, dHdz_path,     &
      &    earth_ref,beta_path, spsfunc_path, ref_corr, s_temp,brkpt, &
      &    no_ele, mid, ilo, ihi, cse, Rs, t_script, tau, Rad, Ier)
!
    Integer(i4), intent(in) :: N_LVLS, BAND, N_SPS

    Integer(i4), intent(in) :: SPS_TBL(Nsps,*)

    Integer(i4), intent(out) :: BRKPT, NO_ELE, MID, ILO, IHI, IER

    Real(r8), intent(in) :: FRQ, H_TAN, EARTH_REF, S_TEMP

    Real(r8), intent(in) :: REF_CORR(*)

    Type(path_beta), intent(in) :: BETA_PATH(:)   ! (Nsps)

    Type(path_index), intent(in)  :: NDX_PATH
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)

    Real(r8), intent(out) :: T_SCRIPT(*), TAU(*)
    Real(r8), intent(out) :: CSE, RS, RAD
!
    Integer(i4) :: Ngp1

    Real(r8) :: del_opacity(N2lvl), delta(N2lvl,Nsps)
!
! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in ?_path%values(1...no_ele)

    Ngp1 = Ng + 1
    brkpt = ndx_path%break_point_index
    no_ele = ndx_path%total_number_of_elements

    EarthX = .false.
    ht = h_tan
    Rr = ht + RoC
    ht2 = Rr * Rr
    if (h_tan < -0.01) then
      Rr = Rr / RoC
      cse = earth_ref
      Call Earth_Intersection(Rs)
    else
      cse = 1.0
      Rr = 0.0d0
      Phi_s = Phi_tan
      NPhi_s = NPhi_tan
    end if
!
!  Compute the appropriate t_script & dt_scrpt_dnp along the integration
!  path of this tanget:
!
    CALL do_t_script(Ngp1, Frq, s_temp, brkpt, no_ele, t_path, &
   &                 mid, t_script)
!
    Call FAST_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path,beta_path, &
 &       dHdz_path,spsfunc_path,n_sps,N_lvls,sps_tbl(1:,band),Nlvl, &
 &       ref_corr,delta,Ier)
    if (Ier /= 0) Return
!
! Initialize the tau & del_opacity arrays:
!
    tau(1:N2lvl) = 0.0
    del_opacity(1:N2lvl) = 0.0
!
    CALL FAST_ZOPACITY(sps_tbl(1:,band), n_sps, Ngp1, N2lvl, brkpt, &
   &                   no_ele, delta, del_opacity)
!
    Call Scrt_dn(t_script, N_lvls, cse, del_opacity, tau, Rad, mid, &
   &             ilo, ihi)
!
    Return
  End Subroutine RAD_TRAN
end module RAD_TRAN_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
