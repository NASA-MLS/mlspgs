module GET_PFA_DELTA_M
  use ELLIPSE, only: HT, PS, ROC
  use GL6P, only: GW, NG
  use MLSCommon, only: I4, R4, R8
  use PFA_BETA_INTRP_M, only: PFA_BETA_INTRP
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  implicit NONE
  private
  public :: GET_PFA_DELTA
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrates in ZETA Space !
!
  Subroutine GET_PFA_DELTA(z_path, t_path, h_path, phi_path, DHDZ_PATH,   &
 &           N_lvls, jf, ch, n_sps, pfa_beta_coeff, ncoeffs,              &
 &           sps_tbl, Nc, Nlvl, N2lvl, maxaitkenpts, maxpfach, Nsps,      &
 &           no_phi_t, f_basis, SCLoop, ref_corr,                         &
 &           mnp, no_freqs_f, no_phi_f, phi_basis_f,                      &
 &           IndxR, IndxL, path_brkpt, beta_t_power, delta, Ier )
!
    integer(i4), intent(in) :: MAXAITKENPTS
    integer(i4), intent(in) :: MAXPFACH
    real(r8), intent(in) :: Z_PATH(*)
    real(r8), intent(in) :: T_PATH(*)
    real(r8), intent(in) :: H_PATH(*)
    real(r8), intent(in) :: PHI_PATH(*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: JF
    integer(i4), intent(in) :: CH
    integer(i4), intent(in) :: N_SPS
!   Real(r8), intent(in) :: PFA_BETA_COEFF(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r8), intent(in) :: PFA_BETA_COEFF(N2lvl,maxaitkenpts,2,*)
    integer(i4), intent(in) :: NCOEFFS(*)
    integer(i4), intent(in) :: SPS_TBL(*)
    integer(i4), intent(in) :: NC
    integer(i4), intent(in) :: NLVL
    integer(i4), intent(in) :: N2LVL
    integer(i4), intent(in) :: NSPS          ! not used
    integer(i4), intent(in) :: NO_PHI_T      ! not used
    real(r8), intent(in) :: F_BASIS(Nc,*)
    integer(i4), intent(in) :: SCLOOP(2,Nlvl,*)
    real(r8), intent(in) :: REF_CORR(*)
    integer(i4), intent(in) :: MNP
    integer(i4), intent(in) :: NO_FREQS_F(*) ! not used
    integer(i4), intent(in) :: NO_PHI_F(*)
    real(r8), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: INDXR, INDXL
    integer(i4), intent(in) :: PATH_BRKPT(*)
!   Real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r8), intent(in) :: BETA_T_POWER(N2lvl,maxaitkenpts,2,*)
    Real(r8), intent(out) :: DELTA(2*Nlvl,Nc,mnp,*)
    integer(i4), intent(out) :: IER
    real(r8) :: H_GL(Ng), PHI_GL(Ng), T_GL(Ng), Z_GL(Ng)
    real(r8) :: BETA_GL(Ng,20), BETA_ZS(20)
    Real(r8) :: GW_DHDZ(Ng)
    integer(i4) :: H_I
    real(r8) :: HD
    integer(i4) :: HEND
    real(r8) :: HH, HL
    Real(r8) :: HTAN2
    integer(i4) :: I, J, K, KP1, KP2, M, MP, NGP1
    real(r8) :: PH, PL, RC
    Real(r8) :: SA, SB
    integer(i4) :: SPS_I
    real(r8) :: TH, TL, TP_GL(Ng,20), TP_ZS(20)
    real(r8) :: ZH, ZL
!
    Ier = 0
    hd = ht + RoC
    htan2 = hd * hd
!
!  Initialize all arrays:
!
    tp_zs = 0.0
    beta_zs = 0.0
    tp_GL = 0.0
    beta_GL = 0.0
!
    hend = 2 * N_lvls - 1
!
    do sps_i = 1, n_sps
      j = sps_tbl(sps_i)
      delta(1:hend + 1,1:ncoeffs(j),1:no_phi_f(j),j) = 0.0
    end do
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
      k = min(IndxR,h_i+1)
      do sps_i = 1, n_sps
        j = sps_tbl(sps_i)
        tp_zs(j) = beta_t_power(k,jf,ch,j)
        beta_zs(j) = pfa_beta_coeff(k,jf,ch,j)
      end do
!
      m = max(1,N_lvls-h_i)
      rc = ref_corr(m)
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power(1,IndxR)
!
      call gauss_legendre(ph,zh)
!
!
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
    Kp1 = mp
    Kp2 = hend * Ngp1 - Ng + 1
!
    zh = z_path(mp)
    th = t_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do while (h_i < hend)
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
      k = max(h_i,IndxL)
      do sps_i = 1, n_sps
        j = sps_tbl(sps_i)
        tp_zs(j) = beta_t_power(h_i+1,jf,ch,j)
        beta_zs(j) = pfa_beta_coeff(h_i+1,jf,ch,j)
      end do
!
      m = max(1,h_i-N_lvls+1)
      rc = ref_corr(m)
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power(IndxL,hend)
!
      call gauss_legendre(pl,zl)
!
      h_i = h_i + 1
!
    end do
!
    Return
  ! *****     Internal procedures     **********************************
  contains
! ----------------------------     COMPUTE_BETA_AND_TEMP_POWER     -----
    subroutine COMPUTE_BETA_AND_TEMP_POWER(I1,I2)
      integer(i4) :: I1, I2
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      integer :: I, J, SPS_I
!
      do sps_i = 1, n_sps
        j = sps_tbl(sps_i)
        do i = 1, Ng
          Call pfa_beta_intrp(z_GL(i),t_GL(i),z_path,t_path,I1, &
   &           I2,Kp1,Kp2,Ngp1,pfa_beta_coeff(1,jf,ch,j),        &
   &           beta_t_power(1,jf,ch,j),beta_GL(i,j),tp_GL(i,j))
        end do
      end do
    end subroutine COMPUTE_BETA_AND_TEMP_POWER
! --------------------------------     DEFINE_GL_GRID_ENTITIES     -----
    subroutine DEFINE_GL_GRID_ENTITIES
!
! Define the various GL grid entities for this sub-interval:
!
      real(r8) :: AYM
      integer :: I, J
      Real(r8) :: YM, Z1, Z2
!
      z1 = dble(zl)
      z2 = dble(zh)
      ym = 0.5_r8 * (z2 - z1)
      aym = abs(ym)
!
      j = mp - Ng
      do i = 1, Ng
        z_GL(i) = z_path(j)
        t_GL(i) = t_path(j)
        h_GL(i) = h_path(j)
        phi_GL(i) = phi_path(j)
        Gw_dHdZ(i) = Gw(i) * dble(DHDZ_PATH(j)) * aym
        j = j + 1
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    subroutine GAUSS_LEGENDRE(FS,ZS)
      real(r8), intent(in) :: FS
      real(r8), intent(in) :: ZS
!
      real(r8) :: BZS, DS
      real(r8) :: ETAPHS, ETAZ, ETFI
      real(r8) :: FV
      integer :: FI, I, J, K, LC1, LC2, NCO, NPF
      real(r8) :: PHI, RZ
      real(r8) :: SING
      integer :: SPS_I
      real(r8) :: SUM
!
      do sps_i = 1, n_sps
!
        j = sps_tbl(sps_i)
        npf = no_phi_f(j)
        nco = ncoeffs(j)
!
        bzs = dble(beta_zs(j))
!
        lc1 = SCLoop(1,m,sps_i)
        lc2 = SCLoop(2,m,sps_i)
!
! Loop over the specie's Phi's
!
        do fi = 1, npf
!
          Call get_one_eta(fs,phi_basis_f(1:npf,j),npf,fi,etaphs)
!
! Compute this specie contribution to the integral:
!
          do k = lc1, lc2
!
! Now k is a coefficient value. Note that f_basis(k,i) is the weighting
! function peak for species i coefficient k .
!
            Sum = 0.0
!
            sing = 0.0
            if (etaphs > 0.0) then
              Call get_one_eta(zs,f_basis(1:nco,j),nco,k,etaz)
              if (etaz > 0.0) sing = bzs * etaz * etaphs
            end if
!
! Now compute the Gauss-Legendre quadrature:
!
            do i = 1, Ng
              fv = -sing
              rz = z_GL(i)
              phi = phi_GL(i)
              Call get_one_eta(rz,f_basis(1:nco,j),nco,k,etaz)
              if (etaz > 0.0) then
                Call get_one_eta(phi,phi_basis_f(1:npf,j),npf,fi,etfi)
                if (etfi > 0.0) then
                  fv = fv + dble(beta_GL(i,j))*etaz*etfi
                end if
              end if
              hd = dble(h_GL(i)) + RoC
              ds = hd / Sqrt(hd*hd-htan2)
              Sum = Sum + Gw_dHdz(i) * fv * ds
            end do
!
            fv = Sum + sing * abs(sb-sa)
            delta(h_i,k,fi,j) = fv * rc
!
          end do
!
        end do
!
      end do
    end subroutine GAUSS_LEGENDRE
  End Subroutine GET_PFA_DELTA
end module GET_PFA_DELTA_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
