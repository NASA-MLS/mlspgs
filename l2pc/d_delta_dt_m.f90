module D_DELTA_DT_M
  use CS_INTRP_M, only: CS_INTRP
  use ELLIPSE, only: HT, PS, ROC
  use GL6P, only: GW, NG
  use MLSCommon, only: I4, R4, R8
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use TWO_D_POLATE_M, only: TWO_D_POLATE
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
  Subroutine d_delta_dt ( z_path, t_path, h_path, phi_path, dHdz_path,      &
 &           dh_dt_path, Kgp, N_lvls, cs, Frq, n_sps, ncoeffs, sps_tbl, Nc, &
 &           Nlvl, f_basis, ref_corr, t_basis, t_profile,                   &
 &           no_t, t_phi_basis, no_phi_t, mdb_pres, mdb_temp, mdb_freq,     &
 &           mnp, mnf, cs_mnf, no_freqs_f, no_phi_f, phi_basis_f, mr_f, in, &
 &           ip, IndxR, IndxL, path_brkpt, d_delta_dtnp, Ier )
!
    integer(i4), intent(in) :: KGP, MNF, CS_MNF, MNP, NLVL
    integer(i4), intent(in) :: NO_PHI_T
!
    real(r8), intent(in) :: Z_PATH(*), T_PATH(*), H_PATH(*), PHI_PATH(*)
    real(r4), intent(in) :: DH_DT_PATH(Kgp,mnp,*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    integer(i4), intent(in) :: N_LVLS
    real(r8), intent(inout) :: CS(Nlvl,no_phi_t,cs_mnf,*)
    real(r8), intent(in) :: FRQ
    integer(i4), intent(in) :: N_SPS, NCOEFFS(*), SPS_TBL(*)
    integer(i4), intent(in) :: NC
    real(r8), intent(in) :: F_BASIS(Nc,*)
    real(r8), intent(in) :: REF_CORR(*)
    real(r8), intent(in) :: T_BASIS(*)
    real(r8), intent(in) :: T_PROFILE(Nc,*) ! not used
    integer(i4), intent(in) :: NO_T
    real(r8), intent(in) :: T_PHI_BASIS(*)
    real(r8), intent(in) :: MDB_PRES(*)
    real(r8), intent(in) :: MDB_TEMP(*)
    real(r8), intent(in) :: MDB_FREQ(cs_mnf,*)
    integer(i4), intent(in) :: NO_FREQS_F(*) ! not used
    integer(i4), intent(in) :: NO_PHI_F(*)
    real(r8), intent(in) :: PHI_BASIS_F(mnp,*)
    real(r8), intent(in) :: MR_F(Nc,mnp,*)
    integer(i4), intent(in) :: IN, IP
    integer(i4), intent(in) :: INDXR, INDXL, PATH_BRKPT(*)
    real(r8), intent(out) :: D_DELTA_DTNP(*)
    integer(i4), intent(in) :: IER ! not used
!
! -----     Local variables     ----------------------------------------
!
    real(r8) :: H_GL(Ng), PHI_GL(Ng), T_GL(Ng), Z_GL(Ng)
    real(r8) :: BETA_GL(Ng,20) = 0.0, BETA_ZS(20) = 0.0
    Real(r8) :: DHDT_GL(Ng), DHDTH, DHDTL
    real(r8) :: ETANP_SING
    real(r8) :: GW_DHDZ(Ng)
    Integer(i4) :: H_I
    Real(r8) :: HH, HL
    real(r8) :: HD, HTAN, HTAN2, HTXDHT
    Integer(i4) :: MP, NGP1
    Real(r8) :: PH, PL
    Real(r8) :: RC
    real(r8) :: SA, SB, SING1, SING2
    real(r8) :: SUM1, SUM2
    Real(r8) :: TH, TL, TP_GL(Ng,20) = 0.0, TP_ZS(20) = 0.0, VETAP(Ng)
    Real(r8) :: ZH, ZL
!
! -----     Executable statements     ----------------------------------
!
!  Define htan and its square, also the dhdt at the tangent:
!
    htan = ht + RoC
    htan2 = htan * htan
!
    ps = -1.0
    Ngp1 = Ng + 1
!
    htxdht = htan * dble(DH_DT_PATH(path_brkpt(1),ip,in))
!
!  Initialize all the arrays:
!
    d_delta_dtnp(1:2*Nlvl) = 0.0
!
! First, do the right hand side of the ray path:
!
    mp = 1
    zh = z_path(mp)
    th = t_path(mp)
    hh = h_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
    dhdth = dble(DH_DT_PATH(mp,ip,in))
!
    do h_i = 1, IndxR
!
      mp = mp + Ngp1
!
      hl = hh
      hh = h_path(mp)
      if (hh < ht) exit
!
      zl = zh
      zh = z_path(mp)
!
      tl = th
      th = t_path(mp)
!
      pl = ph
      ph = phi_path(mp)
!
      sa = sb
      hd = dble(hh) + RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      if (abs(sa-sb) < 0.05) exit
!
      dhdtl = dhdth
      dhdth = dble(DH_DT_PATH(mp,ip,in))
!
      rc = ref_corr(max(1,N_lvls-h_i))
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power ( ph, th, zh )
!
      call thermal_and_singularities ( ph, th, zh )
!
      call gauss_legendre
!
      call add_singularities_back_in
    end do
!
! Second, do the left hand side of the ray path:
!
    ps = 1.0
!
    h_i = IndxL - 1
    mp = path_brkpt(2)
    hh = h_path(mp)
    do while (hh < ht)
      h_i = h_i + 1
      mp = mp + Ngp1
      hh = h_path(mp)
    end do
!
    zh = z_path(mp)
    th = t_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
    dhdth = dble(DH_DT_PATH(mp,ip,in))
!
    do while (h_i < 2 * N_lvls - 1)
!
      mp = mp + Ngp1
!
      hl = hh
      hh = h_path(mp)
!
      sa = sb
      hd = dble(hh) + RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      zl = zh
      zh = z_path(mp)
!
      tl = th
      th = t_path(mp)
!
      pl = ph
      ph = phi_path(mp)
!
      dhdtl = dhdth
      dhdth = dble(DH_DT_PATH(mp,ip,in))
!
      rc = ref_corr(max(1,h_i-N_lvls+1))
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power ( pl, tl, zl )
!
      call thermal_and_singularities ( pl, tl, zl )
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
! ------------------------------     ADD_SINGULARITIES_BACK_IN     -----
    subroutine ADD_SINGULARITIES_BACK_IN
!
! Add the 'singularities factors' back in, multiplied by their
! respective analytical integral:
!
      real(r8) :: DSDT, HD, DS, Q, R
      dsdt = 0.0
      hd = dble(hh) + RoC
      ds = hd * hd - htan2
      if (ds > 0.0) dsdt = (hd*dble(dhdth)-htxdht)/Sqrt(ds)
!
      hd = dble(hl) + RoC
      ds = hd * hd - htan2
      if (ds > 0.0) dsdt = dsdt - (hd*dble(dhdtl)-htxdht)/Sqrt(ds)
!
      q = sum1 + sing1 * abs(sb-sa)
      r = sum2 + sing2 * dsdt * ps
!
      d_delta_dtnp(h_i) = (q + r) * rc       ! for (in,ip)
    end subroutine ADD_SINGULARITIES_BACK_IN
! ----------------------------     COMPUTE_BETA_AND_TEMP_POWER     -----
    subroutine COMPUTE_BETA_AND_TEMP_POWER ( FS, TS, ZS )
      real(r8), intent(in) :: FS
      real(r8), intent(in) :: TS
      real(r8), intent(in) :: ZS
!
      real(r8) :: ETAP, ETAZ
      integer :: I, J, SPS_I
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do sps_i = 1, n_sps
        j = sps_tbl(sps_i)
        Call cs_intrp(zs,ts,Frq,N_lvls,mdb_pres,mdb_temp,         &
   &              mdb_freq(1:cs_mnf,j),cs(:,:,:,j),beta_zs(j),tp_zs(j))
        do i = 1, Ng
          Call cs_intrp(z_GL(i),t_GL(i),Frq,N_lvls,mdb_pres,      &
   &           mdb_temp,mdb_freq(1:cs_mnf,j),cs(:,:,:,j),beta_GL(i,j),   &
   &           tp_GL(i,j))
        end do
      end do
!
      etanp_sing = 0.0
      Call get_one_eta(zs,t_basis,no_t,in,etaz)
      if (etaz > 0.0) then
        etanp_sing = dble(etaz)
        Call get_one_eta(fs,t_phi_basis,no_phi_t,ip,etap)
        etanp_sing = etanp_sing * dble(etap)
      end if
    end subroutine COMPUTE_BETA_AND_TEMP_POWER
! --------------------------------     DEFINE_GL_GRID_ENTITIES     -----
    subroutine DEFINE_GL_GRID_ENTITIES
!
! Define the various GL grid entities for this sub-interval:
!
      real(r8) :: AYM
      integer :: I, J
      real(r8) :: XM, YM, Z1, Z2
!
      z1 = dble(zl)
      z2 = dble(zh)
      xm = 0.5 * (z2 + z1)
      ym = 0.5 * (z2 - z1)
      aym = abs(ym)
!
      j = mp - Ng
      do i = 1, Ng
        z_GL(i) = z_path(j)
        t_GL(i) = t_path(j)
        h_GL(i) = h_path(j)
        phi_GL(i) = phi_path(j)
        dhdt_GL(i) = dble(DH_DT_PATH(j,ip,in))
        Gw_dHdZ(i) = Gw(i) * dble(DHDZ_PATH(j)) * aym
        Call get_one_eta(phi_GL(i),t_phi_basis,no_phi_t,ip,vetap(i))
        j = j + 1
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    subroutine GAUSS_LEGENDRE
!
! Compute the Gauss-Legendre quadrature, subtracting the 'singularities
! factors':
!
      real(r8) :: DS, ETANP
      real(r8) :: ETAZ, F
      real(r8) :: FV1, FV2, HD, HD2, HYD
      integer :: I, J
      real(r8) :: Q
      real(r8) :: PHI
      real(r8) :: R
      real(r8) :: RZ
      integer :: SPS_I
      sum1 = 0.0
      sum2 = 0.0
!
      do i = 1, Ng
!
        rz = z_GL(i)
        phi = phi_GL(i)
!
        Call get_one_eta(rz,t_basis,no_t,in,etaz)
        etanp = (etaz * vetap(i)) / dble(t_GL(i))
!
! Compute the "Hydrostatic" contribution to the derivative:
!
        hd = dble(h_GL(i)) + RoC
        hd2 = hd * hd
        hyd = dble(dhdt_GL(i)) * (2.0*hd2 - 3.0*htan2) + hd*htxdht
        r = hd2 - htan2
        q = Sqrt(r)
        ds = hd / q
        hyd = (hyd / r + hd * etanp) / q
!
        fv1 = -sing1
        fv2 = -sing2
        do sps_i = 1, n_sps
          j = sps_tbl(sps_i)
          Call Two_d_polate(f_basis(1,j),mr_f(1,1,j),Nc,ncoeffs(j),     &
   &                        phi_basis_f(1,j),no_phi_f(j),rz,phi,f)
          if (f /= 0.0) then
            q = dble(f) * beta_GL(i,j)
            if (etanp > 0.0) fv1 = fv1 + q * etanp * tp_GL(i,j)
            fv2 = fv2 + q
          end if
        end do
!
! The final integrand:
!
        sum1 = sum1 + Gw_dHdz(i) * fv1 * ds
        sum2 = sum2 + Gw_dHdz(i) * fv2 * hyd
!
      end do
    end subroutine GAUSS_LEGENDRE
! ------------------------------     THERMAL_AND_SINGULARITIES     -----
    subroutine THERMAL_AND_SINGULARITIES ( FS, TS, ZS )
!
! Compute the contribution of the thermal sensitivity of the absorption:
!
      real(r8), intent(in) :: FS
      real(r8), intent(in) :: TS
      real(r8), intent(in) :: ZS
!
      Real(r8) :: F
      integer :: J
      real(r8) :: Q
      integer :: SPS_I
!
      sing1 = 0.0
      sing2 = 0.0
      do sps_i = 1, n_sps
!
        j = sps_tbl(sps_i)
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
! First: Two-Dimensional interpolation of the mixing ratio (in: zeta & phi)
!
        Call Two_d_polate(f_basis(1,j),mr_f(1,1,j),Nc,ncoeffs(j),      &
   &                      phi_basis_f(1,j),no_phi_f(j),zs,fs,f)
!
! Now compute the 'singularities factors':
!
        q = dble(f) * beta_zs(j)
        sing1 = sing1 + q * etanp_sing * tp_zs(j) / ts
        sing2 = sing2 + q
!
      end do
    end subroutine THERMAL_AND_SINGULARITIES
  End Subroutine D_DELTA_DT
end module D_DELTA_DT_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
