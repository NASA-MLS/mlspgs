module D_DELTA_DT_M
  use ELLIPSE_M, only: ELLIPSE
  use GL6P, only: GW, NG
  use MLSCommon, only: I4, R8
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  implicit NONE
  private
  public :: D_DELTA_DT
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
!  This routine computes the derivative of delta w.r.t. to Temperature in
!  both dimensions (zeta & phi). The Integration is done using the
!  Gauss-Legendre method.  THIS IS THE SLOW VERSION !
!
!  ** NOTE: This routine integrates in ZETA Space !!
!
  Subroutine d_delta_dt(mid,brkpt,no_ele,z_path,t_path,h_path,phi_path, &
 &           beta_path,dHdz_path,dh_dt_path,N_lvls,n_sps,Nlvl,  &
 &           ref_corr,t_z_basis,no_t,t_phi_basis,no_phi_t,spsfunc_path, &
 &           in,ip,elvar,d_delta_dtnp)
!
    Integer(i4), intent(in) :: NLVL,NO_PHI_T,N_LVLS,N_SPS, &
   &             NO_T,IN,IP,MID,BRKPT,NO_ELE
!
    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)
!
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH

    real(r8), intent(in) :: T_Z_BASIS(:), T_PHI_BASIS(:),DH_DT_PATH(:)
    real(r8), intent(in) :: REF_CORR(:)

    Type(ELLIPSE), intent(in out) :: elvar

    real(r8), intent(out) :: D_DELTA_DTNP(:)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: H_I, MP, NGP1, kk

    Real(r8) :: TP_GL(Ng,20) = 0.0, TP_ZS(20) = 0.0
    Real(r8) :: BETAxF_GL(Ng,20) = 0.0, BETAxF_ZS(20) = 0.0

    Real(r8) :: H_GL(Ng), PHI_GL(Ng), T_GL(Ng), Z_GL(Ng)
    Real(r8) :: DHDT_GL(Ng), GW_DHDZ(Ng), VETAP(Ng)

    Real(r8) :: ETANP_SING, DHDTH, DHDTL, HH, HL, HD, HTAN, HTAN2, &
   &            HTXDHT,PH, PL, RC, SA, SB, SING1, SING2, SUM1, SUM2, &
   &            ZH, ZL, TH, TL
!
! -----     Executable statements     ----------------------------------
!
!  Define htan and its square, also the dhdt at the tangent:
!
    htan = elvar%ht + elvar%RoC
    htan2 = htan * htan
!
    elvar%ps = -1.0
    Ngp1 = Ng + 1
!
    htxdht = htan * dh_dt_path(brkpt)
!
!  Initialize all the arrays:
!
    d_delta_dtnp(1:2*Nlvl) = 0.0
!
! First, do the right hand side of the ray path:
!
    mp = 1
    zh = z_path%values(mp)
    th = t_path%values(mp)
    hh = h_path%values(mp)
    ph = phi_path%values(mp)

    hd = hh + elvar%RoC
    sb = Sqrt(abs(hd*hd-htan2))
    dhdth = dh_dt_path(mp)
!
    do h_i = 1, mid
!
      mp = mp + Ngp1
      if(mp > brkpt) EXIT
!
      hl = hh
      hh = h_path%values(mp)
      if (hh < elvar%ht) EXIT
!
      zl = zh
      zh = z_path%values(mp)
!
      tl = th
      th = t_path%values(mp)
!
      pl = ph
      ph = phi_path%values(mp)
!
      sa = sb
      hd = hh + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      if (abs(sa-sb) < 0.05) EXIT
!
      dhdtl = dhdth
      dhdth = dh_dt_path(mp)
!
      rc = ref_corr(h_i+1)
!
      call define_gl_grid_entities
!
      kk = mp
      call compute_beta_and_temp_power (ph,zh)
!
      call thermal_and_singularities (th)
!
      call gauss_legendre
!
      call add_singularities_back_in

    end do
!
! Second, do the left hand side of the ray path:
!
    elvar%ps = 1.0
    h_i = mid
    mp = brkpt + 1
    hh = h_path%values(mp)
    do while (hh < elvar%ht)
      h_i = h_i + 1
      mp = mp + Ngp1
      hh = h_path%values(mp)
    end do
!
    zh = z_path%values(mp)
    th = t_path%values(mp)
    ph = phi_path%values(mp)

    hd = hh + elvar%RoC
    sb = Sqrt(abs(hd*hd-htan2))
    dhdth = dh_dt_path(mp)
!
    do while (h_i < 2 * N_lvls - 1)
!
      mp = mp + Ngp1
      if(mp > no_ele) EXIT
!
      hl = hh
      hh = h_path%values(mp)
!
      sa = sb
      hd = hh + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      zl = zh
      zh = z_path%values(mp)
!
      tl = th
      th = t_path%values(mp)
!
      pl = ph
      ph = phi_path%values(mp)
!
      dhdtl = dhdth
      dhdth = dh_dt_path(mp)
!
      rc = ref_corr(h_i)
!
      call define_gl_grid_entities
!
      kk = mp - Ngp1
      call compute_beta_and_temp_power (pl,zl)
!
      call thermal_and_singularities (tl)
!
      call gauss_legendre
!
      call add_singularities_back_in
!
      h_i = h_i + 1
!
    end do
!
    Return
  ! *****     Internal procedures     **********************************
  contains
! --------------------------------     DEFINE_GL_GRID_ENTITIES     -----
    subroutine DEFINE_GL_GRID_ENTITIES
!
! Define the various GL grid entities for this sub-interval:
!
      integer :: I, J
      real(r8) :: AYM, r
!
      aym = 0.5_r8 * abs(zh - zl)
!
      j = mp - Ng
      do i = 1, Ng
        z_GL(i) = z_path%values(j)
        t_GL(i) = t_path%values(j)
        h_GL(i) = h_path%values(j)
        phi_GL(i) = phi_path%values(j)
        Gw_dHdZ(i) = Gw(i) * dhdz_path%values(j) * aym
        dhdt_GL(i) = dh_dt_path(j)
        Call get_one_eta(phi_GL(i),t_phi_basis,no_phi_t,ip,r)
        vetap(i) = r
        j = j + 1
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! ----------------------------     COMPUTE_BETA_AND_TEMP_POWER     -----
    subroutine COMPUTE_BETA_AND_TEMP_POWER ( FS, ZS )
!
      real(r8), intent(in) :: FS, ZS
!
      real(r8) :: ETAP, ETAZ
      integer :: I, J, K
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do j = 1, n_sps
        tp_zs(j) = beta_path(j)%t_power(kk)
        betaxf_zs(j) = beta_path(j)%values(kk) *  &
       &               spsfunc_path(j)%values(kk)
        k = mp - Ngp1
        do i = 1, Ng
          k = k + 1
          tp_GL(i,j) = beta_path(j)%t_power(k)
          betaxf_GL(i,j) = beta_path(j)%values(k) *  &
         &                 spsfunc_path(j)%values(k)
        end do
      end do
!
      etanp_sing = 0.0
      Call get_one_eta(zs,t_z_basis,no_t,in,etaz)
      if (etaz > 0.0) then
        etanp_sing = etaz
        Call get_one_eta(fs,t_phi_basis,no_phi_t,ip,etap)
        etanp_sing = etanp_sing * etap
      end if
!
    end subroutine COMPUTE_BETA_AND_TEMP_POWER
!
! -----------------------------------------     GAUSS_LEGENDRE     -----
!
    subroutine GAUSS_LEGENDRE
!
! Compute the Gauss-Legendre quadrature, subtracting the 'singularities
! factors':
!
      Integer :: I, J
      Real(r8) :: DS, ETANP, ETAZ, FV1, FV2, HD, HD2, HYD, Q, &
     &            PHI,R,RZ

      sum1 = 0.0
      sum2 = 0.0
!
      do i = 1, Ng
!
        rz = z_GL(i)
        phi = phi_GL(i)
!
        Call get_one_eta(rz,t_z_basis,no_t,in,etaz)
        etanp = (etaz * vetap(i)) / t_GL(i)
!
! Compute the "Hydrostatic" contribution to the derivative:
!
        hd = h_GL(i) + elvar%RoC
        hd2 = hd * hd
        hyd = dhdt_GL(i) * (2.0*hd2 - 3.0*htan2) + hd*htxdht
        r = hd2 - htan2
        q = Sqrt(r)
        ds = hd / q
        hyd = (hyd / r + hd * etanp) / q
!
        fv1 = -sing1
        fv2 = -sing2
        do j = 1, n_sps
          q = betaxf_GL(i,j)
          if (etanp > 0.0) fv1 = fv1 + q * etanp * tp_GL(i,j)
          fv2 = fv2 + q
        end do
!
! The final integrand:
!
        sum1 = sum1 + Gw_dHdz(i) * fv1 * ds
        sum2 = sum2 + Gw_dHdz(i) * fv2 * hyd
!
      end do
    end subroutine GAUSS_LEGENDRE
!
! ------------------------------     THERMAL_AND_SINGULARITIES     -----
!
    Subroutine THERMAL_AND_SINGULARITIES ( TS )
!
! Compute the contribution of the thermal sensitivity of the absorption:
!
      real(r8), intent(in) :: TS
!
      Real(r8) :: Q
      integer :: J
!
      sing1 = 0.0
      sing2 = 0.0

      do j = 1, n_sps
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
        q = betaxf_zs(j)
        sing1 = sing1 + q * etanp_sing * tp_zs(j) / ts
        sing2 = sing2 + q
!
      end do
    end subroutine THERMAL_AND_SINGULARITIES
! ------------------------------     ADD_SINGULARITIES_BACK_IN     -----
    subroutine ADD_SINGULARITIES_BACK_IN
!
! Add the 'singularities factors' back in, multiplied by their
! respective analytical integral:
!
      real(r8) :: DSDT, HD, DS, Q, R

      dsdt = 0.0
      hd = hh + elvar%RoC
      ds = hd * hd - htan2
      if (ds > 0.0) dsdt = (hd*dhdth-htxdht)/Sqrt(ds)
!
      hd = hl + elvar%RoC
      ds = hd * hd - htan2
      if (ds > 0.0) dsdt = dsdt - (hd*dhdtl-htxdht)/Sqrt(ds)
!
      q = sum1 + sing1 * abs(sb-sa)
      r = sum2 + sing2 * dsdt * elvar%ps
!
      d_delta_dtnp(h_i) = (q + r) * rc       ! for (in,ip)
!
    end subroutine ADD_SINGULARITIES_BACK_IN
!
  End Subroutine D_DELTA_DT
!
end module D_DELTA_DT_M
! $Log$
! Revision 1.6  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.5  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
