! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_DELTA_DT_M
  use ELLIPSE_M, only: ELLIPSE
  use GLNP, only: GW, NG
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
 &           beta_path,dHdz_path,dh_dt_path,N_lvls,n_sps,ref_corr,      &
 &           t_z_basis,no_t,t_phi_basis,no_phi_t,spsfunc_path,iz,ip,    &
 &           elvar,no_midval_ndx,midval_ndx,no_gl_ndx,gl_ndx,d_delta_dtnp)
!
    Integer(i4), intent(in) :: NO_PHI_T,N_LVLS,N_SPS,NO_T,IZ,IP,MID, &
                           &   BRKPT,NO_ELE,no_midval_ndx,no_gl_ndx

    Integer(i4), intent(in) :: midval_ndx(:,:),gl_ndx(:,:)
!
    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)
!
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH

    real(r8), intent(in) :: T_Z_BASIS(:),T_PHI_BASIS(:),DH_DT_PATH(:)
    real(r8), intent(in) :: REF_CORR(:)

    Type(ELLIPSE), intent(in out) :: elvar

    real(r8), intent(out) :: D_DELTA_DTNP(:)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: j, k, l, n, h_i, mp, Ngp1, sps_i, hend

    Real(r8) :: TP_GL(Ng,20) = 0.0, TP_ZS(20) = 0.0
    Real(r8) :: BETAxF_GL(Ng,20) = 0.0, BETAxF_ZS(20) = 0.0

    Real(r8) :: H_GL(Ng), PHI_GL(Ng), T_GL(Ng), Z_GL(Ng)
    Real(r8) :: DHDT_GL(Ng), GW_DHDZ(Ng), VETAP(Ng)

    Real(r8) :: ETANP_SING, DHDTH, DHDTL, HH, HL, HD, HTAN, HTAN2, &
   &            HTXDHT,PH, PL, RC, SA, SB, SING1, SING2, SUM1, SUM2, &
   &            ZH, ZL, TH, TL, ds, tpl, tph, fl, fh, dsdt, q ,r
!
! -----     Executable statements     ----------------------------------
!
!  Define htan and its square, also the dhdt at the tangent:
!
    Ngp1 = Ng + 1
    htan = elvar%ht + elvar%RoC
    htan2 = htan * htan
!
    htxdht = htan * dh_dt_path(brkpt)
!
!  Initialize all the arrays:
!
    d_delta_dtnp(1:) = 0.0
!
    if(no_midval_ndx > 0) then
!
      elvar%ps = -1.0
      do j = 1, no_midval_ndx

        mp = midval_ndx(j,2)
        h_i = midval_ndx(j,1)
        if(mp > brkpt) elvar%ps = 1.0
!
        sa = 0.0
        sb = 0.0
        dsdt = 0.0
!
        zl = z_path%values(mp)
        tl = t_path%values(mp)
        hl = h_path%values(mp)
        pl = phi_path%values(mp)
        dhdtl = dh_dt_path(mp)

        hd = hl + elvar%RoC
        q = hd * hd - htan2
        if (q > 0.0) then
          sa = Sqrt(q)
          dsdt = dsdt - (hd*dhdtl-htxdht)/sa
        endif
!
        l = mp + Ngp1
        zh = z_path%values(l)
        th = t_path%values(l)
        hh = h_path%values(l)
        ph = phi_path%values(l)
        dhdth = dh_dt_path(l)
!
        hd = hh + elvar%RoC
        q = hd * hd - htan2
        if (q > 0.0) then
          sb = Sqrt(q)
          dsdt = dsdt + (hd*dhdth-htxdht)/sb
        endif

        ds = abs(sa-sb)
        if (ds < 0.05) EXIT
!
        rc = ref_corr(h_i)
!
        Sum1 = 0.0
        Sum2 = 0.0
        do sps_i = 1, n_sps
          tpl = beta_path(sps_i)%t_power(mp)
          tph = beta_path(sps_i)%t_power(l)
          fl  = beta_path(sps_i)%values(mp) *  &
       &               spsfunc_path(sps_i)%values(mp)
          fh  = beta_path(sps_i)%values(l) *  &
       &               spsfunc_path(sps_i)%values(l)
          Sum1 = Sum1 + 0.5 * (fl*tpl/tl + fh*tph/th)
          Sum2 = Sum2 + 0.5 * (fl + fh)
        end do
!
        q = Sum1 * abs(sb-sa)
        r = Sum2 * dsdt * elvar%ps
        d_delta_dtnp(h_i) = (q + r) * rc       ! for (iz,ip)
!
      end do
!
    endif

    if(no_gl_ndx < 1) Return
!
! Now, do the GL deltas:
! First, do the right hand side of the ray path:
!
    elvar%ps = -1.0
!
    do j = 1, no_gl_ndx
!
      mp = gl_ndx(j,2)
      if (mp >= brkpt) EXIT
!
      h_i = gl_ndx(j,1)
      if(h_i >= mid) EXIT
!
      zl = z_path%values(mp)
      tl = t_path%values(mp)
      hl = h_path%values(mp)
      pl = phi_path%values(mp)
      dhdtl = dh_dt_path(mp)

      sa = 0.0
      hd = hl + elvar%RoC
      q = hd*hd - htan2
      if(q > 0.0) sa = Sqrt(q)
!
      l = mp + Ngp1
      zh = z_path%values(l)
      th = t_path%values(l)
      hh = h_path%values(l)
      ph = phi_path%values(l)
      dhdth = dh_dt_path(l)
!
      sb = 0.0
      hd = hh + elvar%RoC
      q = hd*hd - htan2
      if(q > 0.0) sb = Sqrt(q)
!
      ds = abs(sa-sb)
      if (ds < 0.05) EXIT
!
      rc = ref_corr(h_i+1)
!
      call define_gl_grid_entities
!
      n = l
      call compute_beta_and_temp_power (ph,zh)
!
      call thermal_and_singularities (th)
!
      call gauss_legendre
!
      k = j
      call add_singularities_back_in
!
    end do

    if(k == no_gl_ndx) Return
!
    j = k
    elvar%ps = 1.0
    do
      j = j + 1
      if(j >= no_gl_ndx) EXIT
      mp = gl_ndx(j,2)
      if(mp - l == 1) EXIT
      k = j
    end do
!
! Second, do the left hand side of the ray path:
!
    hend = 2 * N_lvls - 1
    do j = k+1, no_gl_ndx
!
      mp = gl_ndx(j,2)
      l = mp + Ngp1
      if (l > no_ele) EXIT
!
      h_i = gl_ndx(j,1)
      if(h_i >= hend) Return
!
      zl = z_path%values(mp)
      tl = t_path%values(mp)
      hl = h_path%values(mp)
      pl = phi_path%values(mp)
      dhdtl = dh_dt_path(mp)
!
      sa = 0.0
      hd = hl + elvar%RoC
      q = hd*hd - htan2
      if(q > 0.0) sa = Sqrt(q)
!
      zh = z_path%values(l)
      th = t_path%values(l)
      hh = h_path%values(l)
      ph = phi_path%values(l)
      dhdth = dh_dt_path(l)
!
      sb = 0.0
      hd = hh + elvar%RoC
      q = hd*hd - htan2
      if(q > 0.0) sb = Sqrt(q)

      ds = abs(sa-sb)
      if (ds < 0.05) EXIT
!
      rc = ref_corr(h_i)
!
      call define_gl_grid_entities
!
      n = mp
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
      j = mp
      do i = 1, Ng
        j = j + 1
        z_GL(i) = z_path%values(j)
        t_GL(i) = t_path%values(j)
        h_GL(i) = h_path%values(j)
        phi_GL(i) = phi_path%values(j)
        Gw_dHdZ(i) = Gw(i) * dhdz_path%values(j) * aym
        dhdt_GL(i) = dh_dt_path(j)
        Call get_one_eta(phi_GL(i),t_phi_basis,no_phi_t,ip,r)
        vetap(i) = r
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! ----------------------------     COMPUTE_BETA_AND_TEMP_POWER     -----
    subroutine COMPUTE_BETA_AND_TEMP_POWER ( FS, ZS )
!
      real(r8), intent(in) :: FS, ZS
!
      real(r8) :: etap, etaz
      integer :: i, j
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do j = 1, n_sps
        tp_zs(j) = beta_path(j)%t_power(n)
        betaxf_zs(j) = beta_path(j)%values(n) *  &
                   &   spsfunc_path(j)%values(n)
        tp_GL(1:Ng,j) = beta_path(j)%t_power(mp+1:mp+Ng)
        betaxf_GL(1:Ng,j) = beta_path(j)%values(mp+1:mp+Ng) *  &
                        &   spsfunc_path(j)%values(mp+1:mp+Ng)
      end do
!
      etanp_sing = 0.0
      Call get_one_eta(zs,t_z_basis,no_t,iz,etaz)
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
        Call get_one_eta(rz,t_z_basis,no_t,iz,etaz)
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
      real(r8) :: DSDT, HD, Q, R

      dsdt = 0.0
      hd = hh + elvar%RoC
      if(sb > 0.0) dsdt = (hd*dhdth-htxdht)/sb
!
      hd = hl + elvar%RoC
      if (sa > 0.0) dsdt = dsdt - (hd*dhdtl-htxdht)/sa
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
! Revision 1.9  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.8  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
