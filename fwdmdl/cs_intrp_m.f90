module CS_INTRP_M
  use L2PCdim, only: NLVL
  use MDBETA, only: MAX_NO_FREQ, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use S_HUNT_M, only: HUNT
  implicit NONE
  private
  public :: CS_INTRP

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

! Interpolate the beta using the Power interpolation model
!
  Subroutine CS_INTRP ( z, t, f, N_lvls, mdb_pres, mdb_temp, mdb_freq, &
 &                      Cs, Cs_Ntrp, t_power )
    real(r4), intent(in) :: Z
    real(r4), intent(in) :: T
    real(r8), intent(in) :: F
    integer(i4), intent(in) :: N_LVLS
    real(r4), intent(in) :: MDB_PRES(*)
    real(r4), intent(in) :: MDB_TEMP(*)
    real(r8), intent(in) :: MDB_FREQ(*)
    real(r4), intent(in) :: CS(Nlvl,no_t_phi,*)
    real(r4), intent(out) :: CS_NTRP, T_POWER
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: FKLO, J, K, M, NDXP, NDXM, TKLO, ZKLO, ZKHI
    Real(r8) :: Q, V
    Real(r4) :: F_POWER, P_POWER, TMP
!
! -----     Executable statements     ----------------------------------
!
!  Find the closest zeta (Log(pressure)) index
!
    Call Hunt(z,mdb_pres,N_lvls,zklo,zkhi)
    if (abs(mdb_pres(zklo)-z) > abs(mdb_pres(zkhi)-z)) zklo = zkhi
!
!  Find the closest temperature index (In the Phi dimension)
!
    m = (no_t_phi+1) / 2
    v = 1500.0
    tmp = mdb_temp(zklo) + k * 10.0
    do k = -(m-1), m-1
      tmp = tmp + 10.0
      q = abs(t-tmp)
      if (q < v) then
        v = q
        tklo = k + m
        if (v == 0.0) exit
      end if
    end do
!
!  Find the closest frequency index
!
    k = 0
    fklo = 1
    v = 1.0e20
    do k = 1, max_no_freq
      q = abs(f-mdb_freq(k))
      if (q < v) then
        v = q
        fklo = k
        if (v == 0.0) exit
      end if
    end do
!
    Cs_Ntrp = Cs(zklo,tklo,fklo)
!
!  Interpolation in zeta
!
    j = min(zklo+1,N_lvls)
    if (j == zklo) j = max(zklo-1,1)
    v = Cs(j,tklo,fklo) / Cs_Ntrp
    if (abs(v-1.0) > 1.0e-7_r8) then
      p_power = Log(v) / (mdb_pres(j) - mdb_pres(zklo))
      q = p_power * (z - mdb_pres(zklo))
      Cs_Ntrp = Cs_Ntrp * exp(q)
    end if
!
!  Interpolation in temperature
!
    ndxm = tklo
    ndxp = tklo
    t_power = 0.0
    tmp = mdb_temp(zklo) + (tklo - m) * 10.0
    if (tklo < no_t_phi) then
      if (tklo > 1 .and. t < tmp) then
        ndxm = tklo - 1
        q = tmp / (tmp - 10.0)
      else
        ndxp = tklo + 1
        q = (tmp + 10.0) / tmp
      end if
    else
      ndxm = tklo - 1
      q = tmp / (tmp - 10.0)
    end if
!
    if (ndxp > ndxm) then
      v = Cs(zklo,ndxp,fklo) / Cs(zklo,ndxm,fklo)
      if (abs(v-1.0) > 1.0e-7_r8) then
        t_power = Log(v) / Log(q)
        q = t / tmp
        Cs_Ntrp = Cs_Ntrp * (q**t_power)
      end if
    end if
!
!  Interpolation in Frequency
!
    ndxm = fklo
    ndxp = fklo
    q = mdb_freq(fklo)
    if (fklo < max_no_freq) then
      if (fklo > 1 .and. f <= q) then
        ndxm = fklo - 1
      else if (f > q .and. mdb_freq(fklo+1) > 0.1) then
        ndxp = fklo + 1
      end if
    else
      ndxm = fklo - 1
    end if
!
    if (ndxp > ndxm) then
      v = Cs(zklo,tklo,ndxp) / Cs(zklo,tklo,ndxm)
      q = mdb_freq(ndxp) / mdb_freq(ndxm)
      if (abs(v-1.0) > 1.0e-7_r8) then
        f_power = Log(v) / Log(q)
        q = f / mdb_freq(fklo)
        Cs_Ntrp = Cs_Ntrp * (q**f_power)
      end if
    end if
!
    Return
  End Subroutine CS_INTRP
end module CS_INTRP_M

! $Log$
