module CREATE_PFA_BETA_COEFF_M
  use CREATE_BETA_COEFF_M, only: CREATE_BETA_COEFF
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES
  use Gl6P, only: NG
  use L2PC_PFA_STRUCTURES, only: MAXAITKENPTS, PFA_SLAB
  use L2PCDIM, only: N2LVL
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: Create_pfa_beta_coeff

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  Subroutine CREATE_PFA_BETA_COEFF ( jch, band, IndxL, IndxR, N_lvls,   &
 &           sps_tbl, z_path, t_path, path_brkpt, no_int_frqs, f_grid,  &
 &           pfa_spectrum, mdb_hdr, mdb_rec, pfa_beta_coeff,            &
 &           pfa_dbeta_dw, pfa_dbeta_dn, pfa_dbeta_dnu, beta_t_power )

    integer(i4), intent(in) :: JCH, BAND, INDXL, INDXR, N_LVLS
    integer(i4), intent(in) :: SPS_TBL(*)
    real(r4), intent(in) :: Z_PATH(*), T_PATH(*)
    integer(i4), intent(in) :: PATH_BRKPT(*), NO_INT_FRQS
    real(r8), intent(in) :: F_GRID(*)
    type(pfa_slab), intent(in) :: PFA_SPECTRUM(6,*)
    type(eos_mdb_hdr), intent(in) :: MDB_HDR(*)
    type(eos_mdb_rec), intent(in) :: MDB_REC(max_no_lines,*)
!   real(r4), intent(out) :: PFA_BETA_COEFF(N2lvl,maxaitkenpts,maxpfach,*)
    real(r4), intent(out) :: PFA_BETA_COEFF(N2lvl,maxaitkenpts,2,*)
!   Real(r4), intent(out) :: PFA_DBETA_DW(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r4), intent(out) :: PFA_DBETA_DW(N2lvl,maxaitkenpts,2,*)
!   Real(r4), intent(out) :: PFA_DBETA_DN(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r4), intent(out) :: PFA_DBETA_DN(N2lvl,maxaitkenpts,2,*)
!   Real(r4), intent(out) :: PFA_DBETA_DNU(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r4), intent(out) :: PFA_DBETA_DNU(N2lvl,maxaitkenpts,2,*)
!   Real(r4), intent(out) :: BETA_T_POWER(N2lvl,maxaitkenpts,maxpfach,*)
    Real(r4), intent(out) :: BETA_T_POWER(N2lvl,maxaitkenpts,2,*)
!
! -----     Local variables     ----------------------------------------
!
    Real(r8) :: Freq
    Integer(i4) :: frq_i, h_i, hend, j, mp, Ngp1, no_sps, sps_i
    Real(r4) :: pn, pnu, pw, q, t, t_p, z
!
! -----     Executable statements     ----------------------------------
!
    Ngp1 = Ng + 1
    hend = 2 * N_lvls - 1
    no_sps = pfa_spectrum(band,1)%no_sps
!
    do frq_i = 1, no_int_frqs
      Freq = f_grid(frq_i)
      do sps_i = 1, no_sps
        j = sps_tbl(sps_i)
        mp = 1 - Ngp1
        do h_i = 1, IndxR
          q = 0.0
          mp = mp + Ngp1
          z = z_path(mp)
          if(z > -4.5) then
            t = t_path(mp)
            Call Create_beta_coeff(band,sps_i,z,t,Freq,pfa_spectrum, &
   &                               mdb_hdr,mdb_rec,q,t_p,pw,pn,pnu)
          end if
          pfa_beta_coeff(h_i,frq_i,jch,j) = q
          pfa_dbeta_dw(h_i,frq_i,jch,j) = pw
          pfa_dbeta_dn(h_i,frq_i,jch,j) = pn
          pfa_dbeta_dnu(h_i,frq_i,jch,j) = pnu
          beta_t_power(h_i,frq_i,jch,j) = t_p
        end do
        mp = path_brkpt(2) - Ngp1
        do h_i = IndxL, hend
          q = 0.0
          mp = mp + Ngp1
          z = z_path(mp)
          if(z > -4.5) then
            t = t_path(mp)
            Call Create_beta_coeff(band,sps_i,z,t,Freq,pfa_spectrum, &
   &                               mdb_hdr,mdb_rec,q,t_p,pw,pn,pnu)
          end if
          pfa_beta_coeff(h_i,frq_i,jch,j) = q
          pfa_dbeta_dw(h_i,frq_i,jch,j) = pw
          pfa_dbeta_dn(h_i,frq_i,jch,j) = pn
          pfa_dbeta_dnu(h_i,frq_i,jch,j) = pnu
          beta_t_power(h_i,frq_i,jch,j) = t_p
        end do
      end do
    end do
!
    Return
  end Subroutine CREATE_PFA_BETA_COEFF
end module CREATE_PFA_BETA_COEFF_M
! $Log$
