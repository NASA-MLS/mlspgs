module GET_DELTA_M
  use CS_INTRP_M, only: CS_INTRP
  use ELLIPSE, only: HT, PS, ROC
  use GL6P, only: GW, H_GL, NG, PHI_GL, T_GL, Z_GL
  use MDBETA, only: MAX_NO_FREQ
  use MLSCommon, only: I4, R4, R8
  use S_GET_ONE_ETA_M, only: GET_ONE_ETA
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
  Subroutine GET_DELTA ( z_path, t_path, h_path, phi_path, dHdz_path,      &
 &           N_lvls, cs, Frq, n_sps, ncoeffs, sps_tbl, Nc, Nlvl, no_phi_t, &
 &           f_basis, SCLoop, ref_corr, mdb_pres, mdb_temp, mdb_freq,      &
 &           mnp, mnf, no_freqs_f, no_phi_f, phi_basis_f, IndxR, IndxL,    &
 &           path_brkpt, delta, Ier )
    integer(i4), intent(in) :: NO_PHI_T
    integer(i4), intent(in) :: MNF

    real(r4), intent(in) :: Z_PATH(*)
    real(r4), intent(in) :: T_PATH(*)
    real(r4), intent(in) :: H_PATH(*)
    real(r4), intent(in) :: PHI_PATH(*)
    real(r4), intent(in) :: DHDZ_PATH(*)
    integer(i4), intent(in) :: N_LVLS
    real(r4), intent(in) :: CS(Nlvl,no_phi_t,mnf,*)
    real(r8), intent(in) :: FRQ
    integer(i4), intent(in) :: N_SPS
    integer(i4), intent(in) :: NCOEFFS(*)
    integer(i4), intent(in) :: SPS_TBL(*)
    integer(i4), intent(in) :: NC
    integer(i4), intent(in) :: NLVL
    real(r4), intent(in) :: F_BASIS(Nc,*)
    integer(i4), intent(in) :: SCLOOP(2,Nlvl,*)
    real(r4), intent(in) :: REF_CORR(*)
    real(r4), intent(in) :: MDB_PRES(*)
    real(r4), intent(in) :: MDB_TEMP(*)
    real(r8), intent(in) :: MDB_FREQ(mnf,*)
    integer(i4), intent(in) :: MNP
    integer(i4), intent(in) :: NO_FREQS_F(*)   ! Not used
    integer(i4), intent(in) :: NO_PHI_F(*)
    real(r4), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: INDXR, INDXL
    integer(i4), intent(in) :: PATH_BRKPT(*)
    real(r4), intent(inout) :: DELTA(2*Nlvl,Nc,mnp,*)
    integer(i4), intent(out) :: IER
!
! -----     Local variables     ----------------------------------------
!
    real(r4) :: BETA_GL(Ng,20), BETA_ZS(20)
    real(r8) :: GW_DHDZ(Ng)
    integer(i4) :: H_I
    real(r8) :: HD
    real(r4) :: HH, HL
    real(r8) :: HTAN2
    integer(i4) :: J, LC1, M, MP, NGP1
    real(r4) :: PH, PL, RC
    real(r8) :: SA, SB
    integer(i4) :: SPS_I
    real(r4) :: TH, TL, TP_GL(Ng,20), TP_ZS(20)
    real(r4) :: ZH, ZL
!
! -----     Executable statements     ----------------------------------
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
    do sps_i = 1, n_sps
      j = sps_tbl(sps_i)
      delta(1:2 * N_lvls,1:ncoeffs(j),1:no_phi_f(j),j) = 0.0
    end do
!
! First, do the right hand side of the ray path:
!
    ps = -1.0
    Ngp1 = Ng + 1
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
      m = max(1,N_lvls-h_i)
      rc = ref_corr(m)
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power ( th, zh )
!
      call gauss_legendre ( ph, zh )
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
    zh = z_path(mp)
    th = t_path(mp)
    ph = phi_path(mp)
    hd = dble(hh) + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do while (h_i < 2 * N_lvls - 1)
!
      mp = mp + Ngp1
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
      m = max(1,h_i-N_lvls+1)
      rc = ref_corr(m)
!
      call define_gl_grid_entities
!
      call compute_beta_and_temp_power ( tl, zl )
!
      call gauss_legendre ( pl, zl )
!
      h_i = h_i + 1
!
    end do
!
    Return
  ! *****     Internal procedures     **********************************
  contains
! ----------------------------     COMPUTE_BETA_AND_TEMP_POWER     -----
    subroutine COMPUTE_BETA_AND_TEMP_POWER ( TS, ZS )
      real(r4), intent(in) :: TS
      real(r4), intent(in) :: ZS
!
      integer :: I, J, SPS_I
!
! Compute beta and temp. power on the Gauss-Legendre grid for the current
! sub-interval:
!
      do sps_i = 1, n_sps
        j = sps_tbl(sps_i)
        Call cs_intrp(zs,ts,Frq,N_lvls,mdb_pres,mdb_temp,        &
   &              mdb_freq(1:max_no_freq,j),cs(:,:,:,j),beta_zs(j),tp_zs(j))
        do i = 1, Ng
          Call cs_intrp(z_GL(i),t_GL(i),Frq,N_lvls,mdb_pres,     &
   &           mdb_temp,mdb_freq(1:max_no_freq,j),cs(:,:,:,j),beta_GL(i,j),  &
   &           tp_GL(i,j))
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
      real(r8) :: XM, YM, Z1, Z2
!
      z1 = dble(zl)
      z2 = dble(zh)
      xm = 0.5d0 * (z2 + z1)
      ym = 0.5d0 * (z2 - z1)
      aym = abs(ym)
!
      j = mp - Ng
      do i = 1, Ng
        z_GL(i) = z_path(j)
        t_GL(i) = t_path(j)
        h_GL(i) = h_path(j)
        phi_GL(i) = phi_path(j)
        Gw_dHdZ(i) = Gw(i) * dhdz_path(j) * aym
        j = j + 1
      end do
    end subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    subroutine GAUSS_LEGENDRE ( FS, ZS )
      real(r4), intent(in) :: FS
      real(r4), intent(in) :: ZS
!
      real(r8) :: BZS
      real(r4) :: ETAPHS, ETAZ, ETFI
      real(r8) :: FV
      integer :: FI, I, J, K, LC1, LC2, NCO, NPF
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
          Call Get_one_eta(fs,phi_basis_f(1:,j),npf,fi,etaphs)
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
              Call Get_one_eta(zs,f_basis(1:,j),nco,k,etaz)
              if (etaz > 0.0) sing = bzs * etaz * etaphs
            end if
!
! Now compute the Gauss-Legendre quadrature:
!
            do i = 1, Ng
              fv = -sing
              Call Get_one_eta(z_GL(i),f_basis(1:,j),nco,k,etaz)
              if (etaz > 0.0) then
                Call Get_one_eta(phi_GL(i),phi_basis_f(1:,j),npf,fi,etfi)
                if (etfi > 0.0) then
                  fv = fv + dble(beta_GL(i,j))*etaz*etfi
                end if
              end if
              hd = dble(h_GL(i)) + RoC
              Sum = Sum + Gw_dHdz(i) * fv * (hd / Sqrt(hd*hd-htan2))
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
  End Subroutine GET_DELTA
end module GET_DELTA_M

! $Log$
