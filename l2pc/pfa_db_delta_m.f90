module PFA_DB_DELTA_M
  use ELLIPSE, only: HT, PS, ROC
  use GL6P, only: GW, NG
  use MLSCOmmon, only: I4, R4, R8
  use PFA_DBETA_INTRP_M, only: PFA_DBETA_INTRP
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  implicit NONE
  private
  public :: PFA_DB_DELTA
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta_db function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrates in ZETA Space !
!
  Subroutine PFA_DB_DELTA ( z_path, t_path, h_path, phi_path, DHDZ_PATH, &
 &           N_lvls, jf, ch, jsps, ncoeffs, Nc, N2lvl, maxaitkenpts,     &
 &           maxpfach, f_basis, mr_f, ref_corr, mnp, no_phi_f,           &
 &           phi_basis_f, IndxR, IndxL, path_brkpt, beta_t_power,        &
 &           pfa_dbeta_s, s_z_basis, s_phi_basis, s_nz, s_np, iz, ip,    &
 &           delta_s, Ier )
!
    real(r8), intent(in) :: Z_PATH(*)
    real(r8), intent(in) :: T_PATH(*)
    real(r8), intent(in) :: H_PATH(*)
    real(r8), intent(in) :: PHI_PATH(*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: JSPS
    integer(i4), intent(in) :: NCOEFFS(*)
    integer(i4), intent(in) :: NC
    integer(i4), intent(in) :: N2LVL
    integer(i4), intent(in) :: MAXAITKENPTS
    integer(i4), intent(in) :: MAXPFACH
    real(r8), intent(in) :: F_BASIS(Nc,*)
    real(r8), intent(in) :: MR_F(Nc,mnp,*)
    real(r8), intent(in) :: REF_CORR(*)
    integer(i4), intent(in) :: MNP
    integer(i4), intent(in) :: NO_PHI_F(*)
    real(r8), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: INDXR, INDXL
    integer(i4), intent(in) :: PATH_BRKPT(*)
!   Real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,2,*)
!   Real(r8), intent(in) :: PFA_DBETA_S(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r8), intent(in) :: PFA_DBETA_S(N2lvl,maxaitkenpts,2,*)
    real(r8), intent(in) :: S_Z_BASIS(*)
    real(r8), intent(in) :: S_PHI_BASIS(*)
    integer(i4), intent(in) :: S_NZ
    integer(i4), intent(in) :: S_NP
    integer(i4), intent(in) :: IZ
    integer(i4), intent(in) :: IP
    Real(r8), intent(out) :: DELTA_S(*)
    integer(i4), intent(out) :: IER
    real(r8) :: AYM
    integer(i4) :: CH
    real(r8) :: H_GL(Ng), PHI_GL(Ng), T_GL(Ng), Z_GL(Ng)
    real(r8) :: DBETA_S_GL(Ng)
    real(r8) :: DBETA_S_ZS, ETANP_SING
    real(r8) :: DS, ETAP, ETAZ, F, FS, FV_S, GW_DHDZ(Ng)
    integer(i4) :: H_I
    real(r8) :: HD
    integer(i4) :: HEND
    real(r8) :: HH, HL
    real(r8) :: HTAN2
    integer(i4) :: I, J, JF, K, KP1, KP2, M, MP, NGP1, NCO, NPF
    real(r8) :: PH, PHI, PL, RC, RZ
    real(r8) :: SA, SB, SING_S, SUM_S
    real(r8) :: TH, TL, VETAP(Ng), VETAZ(Ng)
    real(r8) :: YM, Z1, Z2
    real(r8) :: ZH, ZL, ZS
!
    Ier = 0
    hd = ht + RoC
    htan2 = hd * hd
!
!  Initialize all arrays:
!
    dbeta_s_GL = 0.0
!
    hend = 2 * N_lvls - 1
    Kp2 = hend + 1
    delta_s(1:Kp2) = 0.0
!
    npf = no_phi_f(jsps)
    nco = ncoeffs(jsps)
!
! First, do the right hand side of the ray path:
!
    ps = -1.0
    Ngp1 = Ng + 1
!
    Kp1 = 1
    Kp2 = IndxR * Ngp1 - Ng
!
    mp = 1
    zh = z_path(mp)
    th = t_path(mp)
    hh = h_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do h_i = 1, IndxR
!
      mp = mp + Ngp1
      if (mp > Kp2) exit
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
      zs = zh
      fs = ph
!
      k = min(IndxR,h_i + 1)
      dbeta_s_zs = dble(pfa_dbeta_s(k,jf,ch,jsps))
!
      etanp_sing = 0.0
      Call get_one_eta(zs,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        etanp_sing = dble(etaz)
        Call get_one_eta(fs,s_phi_basis,s_np,ip,etap)
        etanp_sing = etanp_sing * dble(etap)
      end if
!
      m = max(1,N_lvls-h_i)
      rc = ref_corr(m)
      call define_gl_grid_entities
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do i = 1, Ng
        Call pfa_dbeta_intrp(z_GL(i),t_GL(i),z_path,t_path,1,IndxR,    &
   &         Kp1,Kp2,Ngp1,beta_t_power(1,jf,ch,jsps),                  &
   &         pfa_dbeta_s(1,jf,ch,jsps),dbeta_s_GL(i))
      end do
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
! First: Two-Dimensional interpolation of the mixing ratio (in: zeta & phi)
!
      Call Two_d_polate(f_basis(1,jsps),mr_f(1,1,jsps),Nc,nco,         &
   &                    phi_basis_f(1,jsps),npf,zs,fs,f)
!
! Now compute the 'singularities factors':
!
      sing_s = dble(f) * dbeta_s_zs * etanp_sing
!
      call gauss_legendre(0)
!
    end do
!
! Second, do the left hand side of the ray path:
!
    ps = 1.0
!
    h_i = IndxL
    mp = path_brkpt(2)
    hh = h_path(mp)
    do while (hh < ht)
      h_i = h_i + 1
      mp = mp + Ngp1
      hh = h_path(mp)
    end do
!
    Kp1 = mp
    Kp2 = hend * Ngp1 - Ng + 1
!
    zh = z_path(mp)
    th = t_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do while (h_i <= hend)
!
      mp = mp + Ngp1
      if (mp > Kp2) Return
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
      hl = hh
      hh = h_path(mp)
!
      sa = sb
      hd = dble(hh) + RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      zs = zl
      fs = pl
      dbeta_s_zs = dble(pfa_dbeta_s(h_i,jf,ch,jsps))
!
      etanp_sing = 0.0
      Call get_one_eta(zs,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        etanp_sing = dble(etaz)
        Call get_one_eta(fs,s_phi_basis,s_np,ip,etap)
        etanp_sing = etanp_sing * dble(etap)
      end if
!
      m = max(1,h_i-N_lvls)
      rc = ref_corr(m)
      call define_gl_grid_entities
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do i = 1, Ng
        Call pfa_dbeta_intrp(z_GL(i),t_GL(i),z_path,t_path,IndxL,      &
   &         hend,Kp1,Kp2,Ngp1,beta_t_power(1,jf,ch,jsps),             &
   &         pfa_dbeta_s(1,jf,ch,jsps),dbeta_s_GL(i))
      end do
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
! First: Two-Dimensional interpolation of the mixing ratio (in: zeta & phi)
!
      Call Two_d_polate(f_basis(1,jsps),mr_f(1,1,jsps),Nc,nco,         &
   &                    phi_basis_f(1,jsps),npf,zs,fs,f)
!
! Now compute the 'singularities factors':
!
      sing_s = dble(f) * dbeta_s_zs * etanp_sing
!
      call gauss_legendre(-1)
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
      z1 = dble(zl)
      z2 = dble(zh)
      ym = 0.5 * (z2 - z1)
      aym = abs(ym)
!
      j = mp - Ngp1
      do i = 1, Ng
        j = j + 1
        z_GL(i) = z_path(j)
        t_GL(i) = t_path(j)
        h_GL(i) = h_path(j)
        phi_GL(i) = phi_path(j)
        Gw_dHdZ(i) = Gw(i) * dble(DHDZ_PATH(j)) * aym
        rz = z_GL(i)
        phi = phi_GL(i)
        Call get_one_eta(rz,s_z_basis,s_nz,iz,vetaz(i))
        Call get_one_eta(phi,s_phi_basis,s_np,ip,vetap(i))
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    subroutine GAUSS_LEGENDRE(SHFT)
!
    integer(i4), intent(in) :: SHFT
!
! Now compute the Gauss-Legendre quadrature, subtructing the 'singularities
! factors':
!
      Sum_s = 0.0
!
      do i = 1, Ng
!
! Compute the "Hydrostatic" contribution to the derivative:
!
        hd = dble(h_GL(i)) + RoC
        ds = hd / Sqrt(hd*hd-htan2)
!
        fv_s = -sing_s
!
        rz = z_GL(i)
        phi = phi_GL(i)
!
        Call Two_d_polate(f_basis(1,jsps),mr_f(1,1,jsps),Nc,nco,       &
   &                      phi_basis_f(1,jsps),npf,rz,phi,f)
        if (f /= 0.0) then
          fv_s = fv_s + dble(dbeta_s_GL(i)) * dble(f)
        end if
!
! The final integrand:
!
        Sum_s = Sum_s + fv_s * Gw_dHdz(i) * ds
!
      end do
!
! Now add the 'singularities factors' back in, multiplied by their
! respective analytical integral:
!
      fv_s = Sum_s + sing_s * abs(sb-sa)
!
! And Finally - define the delta:
!
      delta_s(h_i+shft) = fv_s * rc          ! for (iz,ip)
!
    end subroutine GAUSS_LEGENDRE
  End Subroutine PFA_DB_DELTA
end module PFA_DB_DELTA_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
