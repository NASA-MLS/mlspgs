module CREATE_BETA_M
  use ABS_CS_AIR_CONT_M, only: ABS_CS_AIR_CONT
  use ABS_CS_H2O_213G_CONT_M, only: ABS_CS_H2O_213G_CONT
  use ABS_CS_LIQ_H2O_M, only: ABS_CS_LIQ_H2O
  use D_HUNT_M, only: HUNT
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_TEMP, &
                     MAX_ZETA, MAX_FREQ
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: Create_beta
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
! *****     Public Subroutine     **************************************
! --------------------------------------     Create_beta     -----
!
  Subroutine Create_beta (Spectag, z, t, Fgr, mdb_hdr, mdb_rec, &
             beta_value, t_power, dbeta_dw, dbeta_dn, dbeta_dnu0, Ier)
!
!  For a given channel, frequency and height, compute beta_value function.
!  This routine should be called for primary and image seperately.
!
    Integer(i4), intent(in) :: SPECTAG
    Real(r8), intent(in) :: Z, T, Fgr

    Type(eos_mdb_hdr), intent(in) :: MDB_HDR
    Type(eos_mdb_rec), intent(in) :: MDB_REC(*)

    Integer(i4), intent(out) :: Ier
    Real(r8), intent(out) :: BETA_VALUE, T_POWER, DBETA_DW
    Real(r8), intent(out) :: DBETA_DN, DBETA_DNU0
!
! -----  Parameters Declaration ----------------------------------------
!
    Integer(i4), parameter :: MZ=44, NZ=max_zeta, NT=max_temp
!
    Real(r8), parameter :: F213 = 2.15e5_r8, DF213 = 1.5e4_r8
    Real(r8), Parameter :: TINY = epsilon(F213)
!
! -----     Local variables     ----------------------------------------
!
    Real(r8) :: DELV, FREQS(MZ)
    Real(r8) :: B_NTRP, DB_DN_NTRP, DB_DNU0_NTRP, DB_DW_NTRP
    Real(r8) :: FREQ, PRESSURE, Q, Q_LOG(3), RA, RN, RNU0, RW

    Integer(i4) :: LN_I, NO_FREQ, NO_PFA_LINES
!
    Real(r8), SAVE :: mdb_pres(NZ), mdb_temp(NT), x_grid(MAX_FREQ)
!
    Real(r8) :: S, T_P, SPS_MASS, SPS_N, SPS_PS, SPS_V0, SPS_W, WD
    Real(r8) :: TEMP, V0S, X2NU(MZ), X2NU_P(MZ), X2NU_M(MZ)
!
!  Setup absorption coefficients function
!
    Ier = 0
    dbeta_dw = 0.0
    dbeta_dn = 0.0
    dbeta_dnu0 = 0.0
    beta_value = tiny
!
!  Setup frequency independent cross-section data
!
    q = -z
    Temp = t
    Pressure = 10.0_r8 ** q
!
    no_pfa_lines = 0
!
! Check for Water with frequency in (2.0d5, 2.3d5) interval
!
    if ( spectag == 18003 .and. abs(Fgr-f213) < df213 ) then
!
!  Water with frequency in (2.0d5, 2.3d5) interval
!
!  NOTE:
!   This is a problem, we do not have the H2O Mixing Ratio at this point,
!   Thus we use the value of zero...
!
      s = 0.0d0
!
      beta_value = abs_cs_h2o_213g_cont(Temp,Pressure,s)
!     t_power = -3.67           ! See code ...
      wd = T + 10.0d0
      ra = abs_cs_h2o_213g_cont(wd,Pressure,s)
      t_power = Log(ra/beta_value)/Log(wd/T)
      Return
!
    end if
!
!  Now get the beta_value:
!
    if (spectag == 18999) then
!
!  Liquid water
!
      beta_value = abs_cs_liq_h2o(Fgr,Temp)
      wd = Temp + 10.0d0
      ra = abs_cs_liq_h2o(Fgr,wd)
      t_power = Log(ra/beta_value)/Log(wd/Temp)
      Return
!
    else if (spectag == 28964) then
!
!  Dry air contribution
!
      beta_value = abs_cs_air_cont(Temp,Pressure,Fgr)
!     t_power = -2.79           ! See code ...
      wd = Temp + 10.0d0
      ra = abs_cs_air_cont(wd,Pressure,Fgr)
      t_power = Log(ra/beta_value)/Log(wd/Temp)
      Return
!
    else if (spectag == 28965) then
!
!  EXTINCTN molecule
!
      beta_value = 1.0
      t_power = 0.0
      Return
!
    end if
!
!  Check for anything but liquid water and dry air:
!
    if(mdb_hdr%Spectag /= Spectag) then
      Ier = 1
      Print 901,Spectag, mdb_hdr%Spectag
 901  format('** Error in subroutine: Create_beta ..',/,          &
     &       '   For some reason mdb_hdr%spectag /= Spectag !',/, &
     &       '   Spectag:',i6,'  mdb_hdr%spectag:',i6)
      Return
    endif

    x_grid(1:max_freq) = 0.0
    mdb_pres(1:nz) = mdb_hdr%Zeta(1:nz)
    mdb_temp(1:nt) = mdb_hdr%Log_Temp(1:nt)
!
    q_log(1) = mdb_hdr%q_log(1)
    q_log(2) = mdb_hdr%q_log(2)
    q_log(3) = mdb_hdr%q_log(3)
!
    sps_mass = Real(Spectag / 1000)
    no_pfa_lines = mdb_hdr%no_lines
!
    do ln_i = 1, no_pfa_lines
!
      sps_n = mdb_hdr%n(ln_i)
      sps_w = mdb_hdr%w(ln_i)
      sps_ps = mdb_hdr%ps(ln_i)
      sps_v0 = mdb_hdr%v0(ln_i)
!
      Call Get_V0s_X2N(Pressure, t, sps_n, sps_ps, sps_w, sps_mass, &
 &                     sps_v0, v0s, q)
      freqs(ln_i) = v0s
      x2nu(ln_i) = q
!
      ra = t + 10.0
      Call Get_V0s_X2N(Pressure, ra, sps_n, sps_ps, sps_w, sps_mass, &
 &                     sps_v0, v0s, q)
      x2nu_p(ln_i) = q
!
      ra = t - 10.0
      Call Get_V0s_X2N(Pressure, ra, sps_n, sps_ps, sps_w, sps_mass, &
 &                     sps_v0, v0s, q)
      x2nu_m(ln_i) = q
!
    end do
!
!  Now get the beta_value for anything but Liquid water, water with
!  frequency in the (2.0d5,2.3d5) region or Dry Air continuum:
!
    ra = 0.0d0
    rw = 0.0d0
    rn = 0.0d0
    rNu0 = 0.0d0
    t_power = 0.0
!
    do ln_i = 1, no_pfa_lines
!
      Freq = freqs(ln_i)
      delv = Fgr - Freq
!
      no_freq = mdb_hdr%no_f_grid(ln_i)
      x_grid(1:no_freq) = mdb_hdr%x_grid(1:no_freq,ln_i)
!
      Call beta_intrp(z, t, delv, nz, nt, no_freq, mdb_pres, mdb_temp, &
   &       x_grid, x2nu(ln_i), x2nu_m(ln_i), x2nu_p(ln_i),             &
   &       mdb_rec(ln_i)%Log_beta, mdb_rec(ln_i)%dLog_beta_dw,         &
   &       mdb_rec(ln_i)%dLog_beta_dn, mdb_rec(ln_i)%dLog_beta_dNu0,   &
   &       B_Ntrp, dB_dw_Ntrp, dB_dn_Ntrp, dB_dNu0_Ntrp,t_p)
!
      ra = ra + B_Ntrp
      rw = rw + dB_dw_Ntrp
      rn = rn + dB_dn_Ntrp
      rNu0 = rNu0 + dB_dNu0_Ntrp
      t_power = t_power + t_p
!
    end do
!
    dbeta_dw = rw
    dbeta_dn = rn
    beta_value = ra
    dbeta_dnu0 = rNu0
    t_power = t_power / no_pfa_lines
!
    Return
  End Subroutine Create_beta
! *****     Private Procedures     *************************************
! ---------------------------------------------     Get_V0S_X2N    -----
! Compute the shifted Centered frequency & the x-to-Nu factor
!
  Subroutine Get_V0s_X2N(p, t, sps_n, sps_ps, sps_w, sps_mass, sps_v0, &
 &                       v0s, x2nu)
!
    Real(r8), intent(in)  :: p, t, sps_n, sps_ps, sps_w, sps_mass, sps_v0
    Real(r8), intent(out) :: v0s, x2nu
!
    Real(r8), parameter :: WF = 3.58117369e-7_r8
    Real(r8), parameter :: FOUR_LN2 = 2.77258872223978_r8
    Real(r8), parameter :: SQRTLN2 = 8.32554611157698e-1_r8
!
    Real(r8) :: onedt, t3t, ns, wd, x1, y, r
!
! Begin code:
!
    onedt = 1.0d0 / t
    t3t = 300.0d0 * onedt
    ns = 0.25d0 + 1.5d0 * sps_n
    v0s = sps_v0 + sps_ps * p * (t3t**ns)
!
    Wd = v0s * Sqrt(T/sps_mass) * wf
    x1 = Sqrtln2 / Wd
    y = x1 * sps_w * p * (t3t**sps_n)
    r = y + Sqrt(y*y + Four_ln2)
    x2nu = 0.5d0 * wd * r / Sqrtln2
!
    Return
  End Subroutine Get_V0s_X2N
! ---------------------------------------------     BETA_INTRP     -----
! Interpolate the beta using the Power interpolation model
!
  Subroutine BETA_INTRP (z, t, dNu, nz, nt, no_freq, mdb_pres, mdb_temp, &
 &           x_grid, x2nu, x2nu_m, x2nu_p, LogB, dLogB_dw, dLogB_dn,     &
 &           dLogB_dNu0,B_Ntrp,dB_dw_Ntrp,dB_dn_Ntrp,dB_dNu0_Ntrp,t_power)
!
    Real(r8), intent(in) :: Z, T, DNU, X2NU, X2NU_M, X2NU_P
!
    Integer(i4), intent(in) :: NZ, NT, NO_FREQ
!
    Real(r8), intent(in) :: MDB_PRES(*), MDB_TEMP(*), X_GRID(*)
!
!=====================================================================
! The following are real*4 since they are mdb_rec entities...
!
    Real(r4),intent(in) :: LogB(nz,nt,*), dLogB_dw(nz,nt,*)
    Real(r4),intent(in) :: dLogB_dn(nz,nt,*), dLogB_dNu0(nz,nt,*)
!
!=====================================================================
!
    Real(r8), intent(out) :: B_Ntrp, dB_dw_Ntrp, dB_dn_Ntrp, dB_dNu0_Ntrp
    Real(r8), intent(out) :: t_power
!
    Integer(i4) :: jz1,jz2,jt1,jt2,jf1,jf2,pjt1
    Integer(i4) :: i,j,k,m,Ndx(5),ONdx(5)
!
    Real(r8) :: LNT, PV, TEMP, X2N
    Real(r8) :: DELT, DELZ, DZ, DT, Q, R, U, V, Y
    Real(r8) :: VB(5), VF(5), VN(5), VNU0(5), VW(5), OVF(5)

    Real(r8), parameter :: One = 1.0_r8
!
    jz1 = -1
    Call Hunt(z,mdb_pres,nz,jz1,jz2)
!
    delz = z - mdb_pres(jz1)
    dz = mdb_pres(jz2) - mdb_pres(jz1)
    u = delz / dz
    if(u < 0.0) then
      u = 0.0
    else if(u > 1.0) then
      u = 1.0
    endif
!
    Temp = t
    x2n = x2nu
!
    jt1 = -1
    jf1 = -1
!
 10 lnt = Log(Temp)
    Call Hunt(lnt,mdb_temp,nt,jt1,jt2)
!
    delt = lnt - mdb_temp(jt1)
    dt = mdb_temp(jt2) - mdb_temp(jt1)
    v = delt / dt
    if(v < 0.0) then
      v = 0.0
    else if(v > 1.0) then
      v = 1.0
    endif
!
    q = dNu / x2n
    Call Hunt(q,x_grid,no_freq,jf1,jf2)
!
    m = 5
    i = m / 2
    j = max(1,jf1-i)
    k = min(j,no_freq-m+1)
    do i = 1, m
      j = k + i - 1
      Ndx(i) = j
      vf(i) = x_grid(j) * x2n
    end do
!
    do i = 1, m
      j = Ndx(i)
      y = (One-u) * (One-v) * LogB(jz1,jt1,j) + &
   &         u    * (One-v) * LogB(jz2,jt1,j) + &
   &      (One-u) *    v    * LogB(jz1,jt2,j) + &
   &         u    *    v    * LogB(jz2,jt2,j)
      vb(i) = y
    end do
!
    Call Rational_Interp(vf,vb,m,dNu,y)
!
    if(Temp == t) then
      B_Ntrp = Exp(y)               ! Beta(z,t,f)
      pv = v
      pjt1 = jt1
      Ovf(1:m) = vf(1:m)
      ONdx(1:m) = Ndx(1:m)
      Temp = t - 10.0
      x2n = x2nu_m
      goto 10
    else if(Temp == t-10.0) then
      vw(1) = y
      Temp = t + 10.0
      x2n = x2nu_p
      goto 10
    else if(Temp == t+10.0) then
      v = Temp / (t - 10.0)
      t_power = (y - vw(1)) / Log(v)
    endif
!
    v = pv
    jt1 = pjt1
    jt2 = jt1 + 1
!
    do i = 1, m
      j = ONdx(i)
      y = (One-u) * (One-v) * dLogB_dw(jz1,jt1,j) + &
   &         u    * (One-v) * dLogB_dw(jz2,jt1,j) + &
   &      (One-u) *    v    * dLogB_dw(jz1,jt2,j) + &
   &         u    *    v    * dLogB_dw(jz2,jt2,j)
      vw(i) = y
      y = (One-u) * (One-v) * dLogB_dn(jz1,jt1,j) + &
   &         u    * (One-v) * dLogB_dn(jz2,jt1,j) + &
   &      (One-u) *    v    * dLogB_dn(jz1,jt2,j) + &
   &         u    *    v    * dLogB_dn(jz2,jt2,j)
      vn(i) = y
      y = (One-u) * (One-v) * dLogB_dNu0(jz1,jt1,j) + &
   &         u    * (One-v) * dLogB_dNu0(jz2,jt1,j) + &
   &      (One-u) *    v    * dLogB_dNu0(jz1,jt2,j) + &
   &         u    *    v    * dLogB_dNu0(jz2,jt2,j)
      vNu0(i) = y
    end do
!
    Call Rational_Interp(Ovf,vw,m,dNu,r)
    dB_dw_Ntrp = r * B_Ntrp
!
    Call Rational_Interp(Ovf,vn,m,dNu,r)
    dB_dn_Ntrp = r * B_Ntrp
!
    Call Rational_Interp(Ovf,vNu0,m,dNu,r)
    dB_dNu0_Ntrp = r * B_Ntrp
!
    Return
  End Subroutine BETA_INTRP
!
!-----------------------------------------     Rational_Interp     -----
!
  Subroutine Rational_Interp ( Xa, Ya, N, x, y)
!
    Integer(i4), intent(in) :: N
!
    Real(r8), intent(in) :: Xa(*), Ya(*), x
!
    Real(r8), intent(out) :: Y
!
    Integer, Parameter :: Nmax=7
    Integer :: I, KHI, KLO, M, NS
!
    Real(r8), Parameter :: Tiny=1.0d-14
!
    Real(r8) :: C(nmax), D(nmax), DD, DY, DYDX, H, HH, T, W, YLIN
!
    y = Ya(1)
    h = abs(x-Xa(1))
    if (h <= Tiny) Return
!
    klo = -1
    Call Hunt(X,Xa,N,klo,khi)
!
    dydx = (Ya(khi)-Ya(klo))/(Xa(khi)-Xa(klo))
    ylin = Ya(klo) + dydx * (x - Xa(klo))
!
    c(1) = Ya(1)
    d(1) = Ya(1) + Tiny
!
    ns = 1
    hh = abs(x-Xa(1))
!
    do i = 2, n
      h = abs(x-Xa(i))
      if (h <= Tiny) then
        y = Ya(i)
        Return
      end if
      c(i) = Ya(i)
      d(i) = Ya(i) + Tiny
      if(h < hh) then
        ns = i
        hh = h
      endif
    end do
!
    y = Ya(ns)
    ns = ns - 1
!
    do m = 1, n - 1
!
      do i = 1, n - m
!
        w = c(i+1) - d(i)
        h = Xa(i+m) - x
        t = (Xa(i) - x) * d(i) / h
        dd = t - c(i+1)
        if (dd == 0.0) then
          y = ylin
          Return
        end if
        dd = w / dd
        d(i) = c(i+1) * dd
        c(i) = t * dd
!
      end do
!
      if (2*ns < n-m) then
        dy = c(ns+1)
      else
        dy = d(ns)
        ns = ns - 1
      end if
!
      y = y + dy
!
    end do
!
    Return
!
  End Subroutine Rational_Interp
end module CREATE_BETA_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
