module CREATE_BETA_COEFF_M
  use ABS_CS_AIR_CONT_M, only: ABS_CS_AIR_CONT
  use ABS_CS_H2O_213G_CONT_M, only: ABS_CS_H2O_213G_CONT
  use ABS_CS_LIQ_H2O_M, only: ABS_CS_LIQ_H2O
  use D_HUNT_M, only: HUNT
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES, MAX_TEMP, MAX_ZETA
  use L2PC_PFA_STRUCTURES, only: PFA_SLAB
  use MLSCommon, only: I4, R4, R8
  use S_HUNT_M, only: HUNT
  use SLABS, only: SLABS_PREP
  implicit NONE
  private
  public :: Create_beta_coeff

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
! *****     Public Subroutine     **************************************
! --------------------------------------     Create_beta_coeff     -----
!
  Subroutine Create_beta_coeff ( band, sps_i, z, t, Fgr, pfa_spectrum, &
 &           mdb_hdr, mdb_rec, beta_coeff, t_power, dbeta_coef_dw,     &
 &           dbeta_coef_dn, dbeta_coef_dnu0 )
!
!  For a given channel, frequency and height, compute beta_coeff function.
!  This routine should be called for primary and image seperately.
!
    integer(i4), intent(in) :: BAND, SPS_I
    real(r4), intent(in) :: Z
    real(r4), intent(in) :: T
    real(r8), intent(in) :: Fgr
    type(pfa_slab), intent(in) :: PFA_SPECTRUM(6,*)
    type(eos_mdb_hdr), intent(in) :: MDB_HDR(*)
    type(eos_mdb_rec), intent(in) :: MDB_REC(max_no_lines,*)
    real(r4), intent(out) :: BETA_COEFF, T_POWER, DBETA_COEF_DW
    real(r4), intent(out) :: DBETA_COEF_DN, DBETA_COEF_DNU0
!
! -----     Local variables     ----------------------------------------
!
    real(r8), parameter :: F213 = 2.15e5_r8, DF213 = 1.5e4_r8
    Real(r8), parameter :: FOUR_LN2 = 2.77258872223978_r8
    Integer(i4), parameter :: MZ=44, NZ=max_zeta, NT=max_temp
    Real(r8), parameter :: SQRTLN2 = 8.32554611157698e-1_r8
    Real(r4), Parameter :: TINY = 1.0e-14
    Real(r8), parameter :: WF = 3.58117369e-7_r8
!
    real(r4) :: B_NTRP, DB_DN_NTRP, DB_DNU0_NTRP, DB_DW_NTRP
    real(r8) :: DELV, DSLABS1_DV0(MZ)
    Real(r8) :: FREQ, FREQS(MZ), PRESSURE, Q, Q_LOG(3), R, RA, RN, RNU0, RW
    Integer(i4) :: I, LN_I, NO_FREQ, NO_PFA_LINES, SPECTAG
    Integer(i4), save :: MS = -1, PREV_SPEC = -1
    Real(r8) :: S, SLABS1(MZ), SPS_DELTA, SPS_EL, SPS_GAMMA, SPS_MASS, SPS_N
    Real(r8) :: SPS_N1, SPS_N2, SPS_PS, SPS_STR, SPS_V0, SPS_W
    Real(r4) :: T_P
    Real(r8) :: TEMP, V0S, WD, X1(MZ), X2NU(MZ), Y(MZ), YI(MZ)
!
!  Setup absorption coefficients function
!
    beta_coeff = tiny
    dbeta_coef_dw = 0.0
    dbeta_coef_dn = 0.0
    dbeta_coef_dnu0 = 0.0
!
!  Setup frequency independent cross-section data
!
    Temp = Dble(t)
    q = -Dble(z)
    Pressure = 10.0_r8 ** q
!
    no_pfa_lines = 0
    Spectag = pfa_spectrum(band,sps_i)%sps_spectag
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
      beta_coeff = abs_cs_h2o_213g_cont(Temp,Pressure,s)
!     t_power = 3.67           ! See code ...
      wd = Temp + 10.0d0
      ra = abs_cs_h2o_213g_cont(wd,Pressure,s)
      t_power = -Log(ra/beta_coeff)/Log(wd/Temp)
      Return
!
    end if
!
!  Now get the beta_coeff:
!
    if (spectag == 18999) then
!
!  Liquid water
!
      beta_coeff = abs_cs_liq_h2o(Fgr,Temp)
      wd = Temp + 10.0d0
      ra = abs_cs_liq_h2o(Fgr,wd)
      t_power = -Log(ra/beta_coeff)/Log(wd/Temp)
      Return
!
    else if (spectag == 28964) then
!
!  Dry air contribution
!
      beta_coeff = abs_cs_air_cont(Temp,Pressure,Fgr)
!     t_power = 2.79           ! See code ...
      wd = Temp + 10.0d0
      ra = abs_cs_air_cont(wd,Pressure,Fgr)
      t_power = -Log(ra/beta_coeff)/Log(wd/Temp)
      Return
!
    else if (spectag == 28965) then
!
!  EXTINCTN molecule
!
      beta_coeff = 1.0
      t_power = 0.0
      Return
!
    end if
!
!  Check for anything but liquid water and dry air:
!
    if (Spectag /= prev_spec .or. ms < 1) then
      i = 0
      ms = 0
      do while(ms < 1.and.i < 20)
        i = i + 1
        if (mdb_hdr(i)%Spectag == Spectag) then
          ms = i
          i = 22
          prev_spec = Spectag
        end if
      end do
    end if
!
    q_log(1) = mdb_hdr(ms)%q_log(1)
    q_log(2) = mdb_hdr(ms)%q_log(2)
    q_log(3) = mdb_hdr(ms)%q_log(3)
!
    sps_mass = real(spectag / 1000)
    no_pfa_lines = pfa_spectrum(band,sps_i)%no_lines
    do ln_i = 1, no_pfa_lines
      sps_n = mdb_hdr(ms)%n(ln_i)
      sps_w = mdb_hdr(ms)%w(ln_i)
      sps_el = mdb_hdr(ms)%el(ln_i)
      sps_ps = mdb_hdr(ms)%ps(ln_i)
      sps_n1 = mdb_hdr(ms)%n1(ln_i)
      sps_n2 = mdb_hdr(ms)%n2(ln_i)
      sps_v0 = mdb_hdr(ms)%v0(ln_i)
      sps_str = mdb_hdr(ms)%log_i(ln_i)
      sps_delta = mdb_hdr(ms)%delta(ln_i)
      sps_gamma = mdb_hdr(ms)%gamma(ln_i)
!
      Call slabs_prep(Temp,sps_mass,sps_v0,sps_el,sps_w,sps_ps, &
   &       Pressure,sps_n,sps_str,q_log,sps_delta,sps_gamma,    &
   &       sps_n1,sps_n2,v0s,x1(ln_i),q,yi(ln_i),slabs1(ln_i),r)
!
      y(ln_i) = q
      dslabs1_dv0(ln_i) = r
      freqs(ln_i) = v0s
!
      wd = wf * sps_v0 * Sqrt(Temp/sps_mass)
      r = q + Sqrt(q*q + Four_ln2)
      x2nu(ln_i) = 0.5d0 * wd * r / Sqrtln2
    end do
!
!  Now get the beta_coeff for anything but Liquid water, water with
!  frequency in the (2.0d5,2.3d5) region or Dry Air continuum:
!
    ra = 0.0d0
    rw = 0.0d0
    rn = 0.0d0
    rNu0 = 0.0d0
    t_power = 0.0
!
    ln_i = 0
    do ln_i = 1, no_pfa_lines
      Freq = freqs(ln_i)
      delv = Fgr - Freq
!
      no_freq = mdb_hdr(ms)%no_f_grid(ln_i)
      Call beta_intrp(z,t,delv,nz,nt,no_freq,                  &
   &       mdb_hdr(ms)%Zeta,mdb_hdr(ms)%Log_Temp,              &
   &       mdb_hdr(ms)%x_grid(1,ln_i),x2nu(ln_i),              &
   &       mdb_rec(ln_i,ms)%Log_beta,                          &
   &       mdb_rec(ln_i,ms)%dLog_beta_dw,                      &
   &       mdb_rec(ln_i,ms)%dLog_beta_dn,                      &
   &       mdb_rec(ln_i,ms)%dLog_beta_dNu0,                    &
   &       mdb_rec(ln_i,ms)%Log_beta_intrp,B_Ntrp,dB_dw_Ntrp,  &
   &       dB_dn_Ntrp,dB_dNu0_Ntrp,t_p,nt)
!
      ra = ra + B_Ntrp
      rw = rw + dB_dw_Ntrp
      rn = rn + dB_dn_Ntrp
      rNu0 = rNu0 + dB_dNu0_Ntrp
      t_power = t_power + t_p
!
    end do
!
    beta_coeff = ra
    dbeta_coef_dw = rw
    dbeta_coef_dn = rn
    dbeta_coef_dnu0 = rNu0
    t_power = t_power / no_pfa_lines
!
    Return
  End Subroutine Create_beta_coeff
! *****     Private Procedures     *************************************
! ---------------------------------------------     BETA_INTRP     -----
! Interpolate the beta using the Power interpolation model
!
  Subroutine BETA_INTRP ( z, t, dNu, no_zeta, no_temp, no_freq,    &
 &           mdb_pres, mdb_temp, x_grid, x2nu, LogB, dLogB_dw,     &
 &           dLogB_dn, dLogB_dNu0, LogB_Intrp, B_Ntrp, dB_dw_Ntrp, &
 &           dB_dn_Ntrp, dB_dNu0_Ntrp, t_power, nt)

    integer(i4), intent(in) :: NT
    real(r4), intent(in) :: Z
    real(r4), intent(in) :: T
    real(r8), intent(in) :: DNU
    integer(i4), intent(in) :: NO_ZETA
    integer(i4), intent(in) :: NO_TEMP, NO_FREQ
    real(r4), intent(in) :: MDB_PRES(*), MDB_TEMP(*)
    real(r8), intent(in) :: X_GRID(*)
    Real(r8), intent(in) :: X2NU
    real(r4), intent(in) :: LogB(no_zeta,nt,*), dLogB_dw(no_zeta,nt,*)
    real(r4), intent(in) :: dLogB_dn(no_zeta,nt,*), dLogB_dNu0(no_zeta,nt,*)
    character, intent(in) :: LogB_Intrp(no_zeta,nt,*)
    real(r4), intent(out) :: B_Ntrp, dB_dw_Ntrp, dB_dn_Ntrp, dB_dNu0_Ntrp
    real(r4), intent(out) :: t_power

    Integer(i4) :: fklo,i,j,k,ik,m,tklo,zklo
    Character :: IType
    Real(r4) :: LNT, P_POWER
    Real(r8) :: DELT, DELZ, DZ, DT, Q, R
    Real(r8) :: V, VB(5), VBE(5), VF(5), VN(5), VNU0(5), VW(5), Y
!
    Call Hunt(z,mdb_pres,no_zeta,zklo,i)
    if (abs(mdb_pres(i)-z) < abs(mdb_pres(zklo)-z)) zklo=i
!
    lnt = Log(t)
    Call Hunt(lnt,mdb_temp,no_temp,tklo,i)
    if (abs(mdb_temp(i)-lnt) < abs(mdb_temp(tklo)-lnt)) tklo=i
!
    q = dNu / x2nu
    Call Hunt(q,x_grid,no_freq,fklo,i)
    if (abs(x_grid(i)-q) < abs(x_grid(fklo)-q)) fklo=i
!
    y = LogB(zklo,tklo,fklo)
    IType = LogB_Intrp(zklo,tklo,fklo)
!
!  Interpolation in zeta
!
    delz = z - mdb_pres(zklo)
    if (abs(delz) >= 1.0e-3) then
      dz = mdb_pres(zklo+1) - mdb_pres(zklo)
      v = LogB(zklo+1,tklo,fklo) - LogB(zklo,tklo,fklo)
      p_power = v / dz
      y = y + p_power * delz
    end if
!
!  Interpolation in temperature
!
    delt = lnt - mdb_temp(tklo)
    dt = mdb_temp(tklo+1) - mdb_temp(tklo)
    v = LogB(zklo,tklo+1,fklo) - LogB(zklo,tklo,fklo)
    t_power = v / dt
    if (delt /= 0.0) y = y + t_power * delt
!
!  Interpolation in Frequency
!
    m = 5
    ik = 1
    i = m/2
    j = max(1,fklo-i)
    k = min(j,no_freq-m+1)
    do i = 1, m
      j = k+i-1
      vb(i) = LogB(zklo,tklo,j)
      vbe(i) = Exp(vb(i))
      vf(i) = x2nu * x_grid(j)
      vw(i) = dLogB_dw(zklo,tklo,j)
      vn(i) = dLogB_dn(zklo,tklo,j)
      vNu0(i) = dLogB_dNu0(zklo,tklo,j)
      if (j == fklo) ik = i
    end do
!
    j = 2
    if (ik < 3) j = 1
    if (IType == 'L') then
      Call Freq_Inverse_Quad_Intrp(vf(j),vbe(j),dNu,r)
      if (r <= 0.0) Call Freq_Gauss_Intrp(vf(j),vb(j),dNu,r)
    else if (IType == 'Q') then
      Call Freq_Quad_Intrp(vf(j),vbe(j),dNu,r)
      if (r <= 0.0) Call Freq_Gauss_Intrp(vf(j),vb(j),dNu,r)
    else
      Call Freq_Gauss_Intrp(vf(j),vb(j),dNu,r)
    end if
!
    q = r / vbe(ik)
!
!  ** Finaly - The actual Beta:
!
    B_Ntrp = Exp(y) * q
!
!  **** Now we interpolate d_beta_dw:
!
    y = dLogB_dw(zklo,tklo,fklo)
!
!  Interpolation in zeta
!
    if (abs(delz) >= 1.0e-3) then
      v = dLogB_dw(zklo+1,tklo,fklo) - dLogB_dw(zklo,tklo,fklo)
      q = v / dz
      y = y + q * delz
    end if
!
!  Interpolation in temperature
!
    if (abs(delt) >= 1.0e-3) then
      v = dLogB_dw(zklo,tklo+1,fklo) - dLogB_dw(zklo,tklo,fklo)
      q = v / dt
      y = y + q * delt
    end if
!
!  Interpolation in Frequency
!
!   Call Freq_Quad_Intrp(vf(j),vw(j),dNu,r)
!
    Call Rational_Interp(vf,vw,m,ik,dNu,r)
    q = r / vw(ik)
!
!  ** Finaly - The actual dLog_Beta_dw:
!
    dB_dw_Ntrp = y * q * B_Ntrp
!
!  **** Now we interpolate d_beta_dw:
!
    y = dLogB_dn(zklo,tklo,fklo)
!
!  Interpolation in zeta
!
    if (abs(delz) >= 1.0e-3) then
      v = dLogB_dn(zklo+1,tklo,fklo) - dLogB_dn(zklo,tklo,fklo)
      q = v / dz
      y = y + q * delz
    end if
!
!  Interpolation in temperature
!
    if (abs(delt) >= 1.0e-3) then
      v = dLogB_dn(zklo,tklo+1,fklo) - dLogB_dn(zklo,tklo,fklo)
      q = v / dt
      y = y + q * delt
    end if
!
!  Interpolation in Frequency
!
!   Call Freq_Quad_Intrp(vf(j),vn(j),dNu,r)
!
    Call Rational_Interp(vf,vn,m,ik,dNu,r)
    q = r / vn(ik)
!
!  ** Finaly - The actual dLog_Beta_dn:
!
    dB_dn_Ntrp = y * q * B_Ntrp
!
!  **** Now we interpolate d_beta_dNu0:
!
    y = dLogB_dNu0(zklo,tklo,fklo)
!
!  Interpolation in zeta
!
    if (abs(delz) >= 1.0e-3) then
      v = dLogB_dNu0(zklo+1,tklo,fklo) - dLogB_dNu0(zklo,tklo,fklo)
      q = v / dz
      y = y + q * delz
    end if
!
!  Interpolation in temperature
!
    if (abs(delt) >= 1.0e-3) then
      v = dLogB_dNu0(zklo,tklo+1,fklo) - dLogB_dNu0(zklo,tklo,fklo)
      q = v / dt
      y = y + q * delt
    end if
!
!  Interpolation in Frequency
!
    Call Rational_Interp(vf,vNu0,m,ik,dNu,r)
    q = r / vNu0(ik)
!
!  ** Finaly - The actual dLog_Beta_dNu0:
!
    dB_dNu0_Ntrp = y * q * B_Ntrp
!
    Return
  End Subroutine BETA_INTRP
!
!---------------------------------     Freq_Inverse_Quad_Intrp     -----
! Interpolate the beta for Frequency (Using Lorentzian, or the equivalent
! Inverse Quadrature model)
!
!  y = A*A/(B*B+(x-C)^2)  or:  y = 1.0d0 / (a+x*(b+x*c))
!
  Subroutine Freq_Inverse_Quad_Intrp ( vf, vb, f, beta )
!
    Real(r8), intent(in) :: VF(*), VB(*), F
    Real(r8), intent(out) :: BETA
    Real(r8) :: A, B, C, Y, Y1, Y2, Y3, ZA, ZB
!
    if (abs(f-vf(1)) <= 0.001) then
      beta = vb(1)
      Return
    else if (abs(f-vf(2)) <= 0.001) then
      beta = vb(2)
      Return
    else if (abs(f-vf(3)) <= 0.001) then
      beta = vb(3)
      Return
    end if
!
    y1 = 1.0 / vb(1)
    y2 = 1.0 / vb(2)
    y3 = 1.0 / vb(3)
!
    za = (y2 - y1) / (vf(2) - vf(1))
    zb = (y3 - y2) / (vf(3) - vf(2))
!
    c = (zb - za) / (vf(3) - vf(1))
    b = za - c * (vf(2) + vf(1))
    a = y2 - vf(2) * (b + c * vf(2))
!
    y = a + f * (b + c * f)
    beta = 1.0 / y
!
    Return
  End Subroutine Freq_Inverse_Quad_Intrp
!
!----------------------------------------     Freq_Gauss_Intrp     -----
! Interpolate the beta for Frequency (Using Gaussian model)
!  y = A * exp(Alfa*(x-Beta)^2)
!
  Subroutine Freq_Gauss_Intrp ( vf, vb, x, y )
    Real(r8), intent(in) :: VF(*), VB(*), X
    Real(r8), intent(out) :: Y
    Real(r8) :: A, AL, ALFA, BETA, BL, CL, R, Y1, Y2, Y3, ZA, ZB
!
    if (abs(x-vf(1)) <= 0.001) then
      y = Exp(vb(1))
      Return
    else if (abs(x-vf(2)) <= 0.001) then
      y = Exp(vb(2))
      Return
    else if (abs(x-vf(3)) <= 0.001) then
      y = Exp(vb(3))
      Return
    end if
!
    y1 = vb(1)
    y2 = vb(2)
    y3 = vb(3)
!
    za = (y2 - y1) / (vf(2) - vf(1))
    zb = (y3 - y2) / (vf(3) - vf(2))
!
    cl = (zb - za) / (vf(3) - vf(1))
    bl = za - cl * (vf(2) + vf(1))
    al = y2 - vf(2) * (bl + cl * vf(2))
!
    Alfa = -cl
    Beta = bl / (2.0 * Alfa)
    r = vf(2) - Beta
    A = Exp(vb(2) + Alfa*r*r)
!
    r = x - Beta
    al = -Alfa * r * r
    if (abs(al) > 85.0) al = sign(1.0_r8,al) * 85.0
    y = A * exp(al)
!
    Return
  End Subroutine Freq_Gauss_Intrp
!
!-----------------------------------------     Freq_Quad_Intrp     -----
! Quadratic interpolation routine: y = a + x * (b + x * c)
!
  Subroutine Freq_Quad_Intrp ( vf, vb, f, r )
    Real(r8), intent(in) :: VF(*), VB(*), F
    Real(r8), intent(out) :: R
    Real(r8) :: A, B, C, Y1, Y2, Y3, ZA, ZB
!
    if (abs(f-vf(1)) <= 0.001) then
      r = vb(1)
      Return
    else if (abs(f-vf(2)) <= 0.001) then
      r = vb(2)
      Return
    else if (abs(f-vf(3)) <= 0.001) then
      r = vb(3)
      Return
    else if (f < vf(1)) then              ! Linear interpolation
      y1 = (f-vf(1)) / (vf(2)-vf(1))
      r = vb(1) + (vb(2)-vb(1)) * y1
      Return
    else if (f > vf(3)) then              ! Linear interpolation
      y3 = (f-vf(2)) / (vf(3)-vf(2))
      r = vb(2) + (vb(3)-vb(2)) * y3
      Return
    end if
!
    y1 = vb(1)
    y2 = vb(2)
    y3 = vb(3)
!
    za = (y2 - y1) / (vf(2) - vf(1))
    zb = (y3 - y2) / (vf(3) - vf(2))
!
    c = (zb - za) / (vf(3) - vf(1))
    b = za - c * (vf(2) + vf(1))
    a = y2 - vf(2) * (b + c * vf(2))
!
    r = a + f * (b + c * f)
!
    Return
  End Subroutine Freq_Quad_Intrp
!
!-----------------------------------------     Rational_Interp     -----
!
  Subroutine Rational_Interp ( Xa, Ya, n, ik, x, y)
    Real(r8), intent(in) :: Xa(*), Ya(*)
    Integer(i4), intent(in) :: N, IK
    Real(r8), intent(in) :: X
    Real(r8), intent(out) :: Y

    Integer, Parameter :: Nmax=7
    Real(r8), Parameter :: Tiny=1.0d-14
    Real(r8) :: C(nmax), D(nmax), DD, DY, DYDX, H
    Integer :: I, KHI, KLO, M, NS
    Real(r8) :: T, W, YLIN
!
    y = Ya(1)
    h = abs(x-Xa(1))
    if (h <= Tiny) Return
!
    klo = min(ik,n-1)
    khi = klo + 1
!
    dydx = (Ya(khi)-Ya(klo))/(Xa(khi)-Xa(klo))
    ylin = Ya(klo) + dydx * (x - Xa(klo))
!
    c(1) = Ya(1)
    d(1) = Ya(1) + Tiny
!
    do i = 2, n
      h = abs(x-Xa(i))
      if (h <= Tiny) then
        y = Ya(i)
        Return
      end if
      c(i) = Ya(i)
      d(i) = Ya(i) + Tiny
    end do
!
    ns = ik
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
end module CREATE_BETA_COEFF_M

! $Log$
