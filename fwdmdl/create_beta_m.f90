module CREATE_BETA_M
  use MLSCommon, only: R8, RP, IP
  use ABS_CS_CONT_M,    only: ABS_CS_CONT
  use ABS_CS_N2_CONT_M, only: ABS_CS_N2_CONT
  use ABS_CS_O2_CONT_M, only: ABS_CS_O2_CONT
  use ABS_CS_LIQ_H2O_M, only: ABS_CS_LIQ_H2O
  use SLABS_SW_M, only: DVOIGT_SPECTRAL, VOIGT_LORENTZ, SLABSWINT, SLABS
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
  Subroutine Create_beta (Spectag, cont, pressure, Temp, Fgr, nl, pfaw,  &
         &   v0s, x1,y, yi, slabs1, beta_value, dslabs1_dv0, v0sp, x1p,  &
         &   yp, yip, slabs1p, v0sm, x1m, ym, yim, slabs1m,              &
         &   t_power, dbeta_dw, dbeta_dn, dbeta_dv)
!
!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image seperately.
!
! Inputs:
  INTEGER(ip), INTENT(in) :: SPECTAG ! molecule id tag
  REAL(rp), INTENT(in) :: cont(:) ! continuum parameters
  REAL(rp), INTENT(in) :: pressure ! pressure in hPa
  REAL(rp), INTENT(in) :: temp ! temperature in K
  REAL(rp), INTENT(in) :: fgr ! frequency in MHz
  INTEGER(ip), INTENT(in) :: nl ! no of lines
  REAL(rp), INTENT(in) :: pfaw(:) ! line widths
  REAL(r8), INTENT(in) :: v0s(:) ! pressure shifted line centers
  REAL(rp), INTENT(in) :: x1(:) ! Doppler width
  REAL(rp), INTENT(in) :: y(:) ! ratio Pressure to Doppler widths
  REAL(rp), INTENT(in) :: yi(:) ! Interference coefficients
  REAL(rp), INTENT(in) :: slabs1(:) ! strengths
! optional inputs for spectral derivatives
  REAL(rp), OPTIONAL, INTENT(in) :: dslabs1_dv0(:) ! strength derivative
!                                wrt line position
! optional inputs for temperature derivatives
  REAL(r8), OPTIONAL, INTENT(in) :: v0sp(:) ! pressure shifted line centers
  REAL(rp), OPTIONAL, INTENT(in) :: x1p(:)! Doppler width
  REAL(rp), OPTIONAL, INTENT(in) :: yp(:) ! ratio Pressure to Doppler widths
  REAL(rp), OPTIONAL, INTENT(in) :: yip(:) ! Interference coefficients
  REAL(rp), OPTIONAL, INTENT(in) :: slabs1p(:) ! strengths
  REAL(r8), OPTIONAL, INTENT(in) :: v0sm(:) ! pressure shifted line centers
  REAL(rp), OPTIONAL, INTENT(in) :: x1m(:)! Doppler width
  REAL(rp), OPTIONAL, INTENT(in) :: ym(:) ! ratio Pressure to Doppler widths
  REAL(rp), OPTIONAL, INTENT(in) :: yim(:) ! Interference coefficients
  REAL(rp), OPTIONAL, INTENT(in) :: slabs1m(:) ! strengths
! outputs
  REAL(rp), INTENT(out) :: beta_value
! optional outputs
  REAL(rp), OPTIONAL, INTENT(out) :: T_POWER ! for temperature derivative
  REAL(rp), OPTIONAL, INTENT(OUT) :: DBETA_DW ! line width derivative
  REAL(rp), OPTIONAL, INTENT(OUT) :: DBETA_DN ! temperature dependence deriv
  REAL(rp), OPTIONAL, INTENT(OUT) :: DBETA_DV ! line position derivative
!
! -----     Local variables     ----------------------------------------
!
    Integer(ip) :: LN_I

    Real(rp) :: w,ra,dNu,tp,bp,tm,bm,bv,dw,dn,ds,dbdw,dbdn,dbdv
!
    bv = 0.0_rp
    bp = 0.0_rp
    bm = 0.0_rp
    beta_value = 0.0_rp
    tp = Temp + 10.0_rp
    tm = Temp - 10.0_rp
!
!  Setup absorption coefficients function
! NEED TO ADD THE PARDO WATER VAPOR CONTINUUM FOR 18003
!  Now get the beta_value:
!
    if (spectag == 18999) then
!
!  Liquid water
!
      bv = abs_cs_liq_h2o(Fgr,Temp)
      IF (PRESENT(t_power)) THEN
        bm = abs_cs_liq_h2o(Fgr,tm)
        bp = abs_cs_liq_h2o(Fgr,tp)
      ENDIF
!
    else if (spectag == 28964) then
!
!  Dry air contribution (N2)
!
      bv = abs_cs_n2_cont(cont,Temp,Pressure,Fgr)
      IF (PRESENT(t_power)) THEN
        bm = abs_cs_n2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_n2_cont(cont,tp,Pressure,Fgr)
      ENDIF
!
    else if (spectag == 28965) then
!
!  EXTINCTN molecule
!
      beta_value = 1.0_rp
      IF (PRESENT(t_power)) t_power = 0.0_rp
      Return 
!
    else if (spectag == 32001) then

      bv = abs_cs_o2_cont(cont,Temp,Pressure,Fgr)
      IF (PRESENT(t_power)) THEN
        bm = abs_cs_o2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_o2_cont(cont,tp,Pressure,Fgr)
      ENDIF

    else

      bv = abs_cs_cont(cont,Temp,Pressure,Fgr)
      IF (PRESENT(t_power)) THEN
        bm = abs_cs_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_cont(cont,tp,Pressure,Fgr)
      ENDIF

    end if
!
    beta_value = bv
    IF(nl < 1) THEN
      IF(PRESENT(t_power)) THEN
        ds = Log(bp/bv)/Log(tp/Temp)     ! Estimate over [temp+10,temp]
        ra = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
        dw = Log(bv/bm)/Log(Temp/tm)     ! Estimate over [temp,temp-10]
        t_power = (ds + 2.0 * ra + dw) / 4.0  ! Weighted Average
      ENDIF
      Return
    ENDIF
!
    IF(PRESENT(DBETA_DW).OR.PRESENT(DBETA_DN).OR.PRESENT(DBETA_DV)) THEN
!
      dbdw = 0.0_rp
      dbdn = 0.0_rp
      dbdv = 0.0_rp

      do ln_i = 1, nl
!
        dNu = Fgr - v0s(ln_i)
!
! If too far from line center, skip it (to fit /mlspgs/ code).
! Use criterion of 3000.0 (NOT 2000.0) so the l2_gridding program works !!
!
        if(abs(dNu) > 3000.0_rp) CYCLE             ! To fit /mlspgs/ code
!
        w = pfaw(ln_i)
        IF(abs(y(ln_i))+0.666666_rp*abs(x1(ln_i)*dNu) > 100.0_rp) THEN
          Call Voigt_Lorentz(dNu,v0s(ln_i),x1(ln_i),yi(ln_i), &
            &  y(ln_i),w,Temp,slabs1(ln_i),bv,dslabs1_dv0(ln_i),dw,dn,ds)
        ELSE
          Call dvoigt_spectral(dNu,v0s(ln_i),x1(ln_i),yi(ln_i),y(ln_i), &
         &     w,Temp,slabs1(ln_i),bv,dslabs1_dv0(ln_i),dw,dn,ds)
        ENDIF
!
        beta_value = beta_value + bv
        dbdw = dbdw + dw
        dbdn = dbdn + dn
        dbdv = dbdv + ds
!
      end do

      IF(PRESENT(DBETA_DW)) dbeta_dw = dbdw
      IF(PRESENT(DBETA_DN)) dbeta_dn = dbdn
      IF(PRESENT(DBETA_DV)) dbeta_dv = dbdv
!
    ELSE                ! No derivatives required
!
      IF(MAXVAL(ABS(yi)) < 1.0e-06_rp) THEN
        do ln_i = 1, nl
          dNu = Fgr - v0s(ln_i)
! To fit /mlspgs/ code
          if(abs(dNu) <= 3000.0_rp) beta_value = beta_value &
            + Slabs(dNu,v0s(ln_i),x1(ln_i),slabs1(ln_i),y(ln_i))
        end do
      ELSE
        do ln_i = 1, nl
          dNu = Fgr - v0s(ln_i)
! To fit /mlspgs/ code
          if(abs(dNu) <= 3000.0_rp) beta_value = beta_value &
             + Slabswint(dNu,v0s(ln_i),x1(ln_i),slabs1(ln_i),y(ln_i),yi(ln_i))
        enddo
      ENDIF
!
    ENDIF

    IF(PRESENT(t_power)) THEN
!
!  Find the temperatue power dependency now:
!
      IF(MAXVAL(ABS(yi)) < 1.0e-06_rp) THEN
        do ln_i = 1, nl
          ds = Fgr - v0s(ln_i)
          IF(abs(ds) <= 3.0e3) THEN
            dNu = Fgr - v0sp(ln_i)
            bp = bp + Slabs(dNu,v0sp(ln_i),x1p(ln_i),slabs1p(ln_i),yp(ln_i))
            dNu = Fgr - v0sm(ln_i)
            bm = bm + Slabs(dNu,v0sm(ln_i),x1m(ln_i),slabs1m(ln_i),ym(ln_i))
          ENDIF
        end do
      ELSE
        do ln_i = 1, nl
          ds = Fgr - v0s(ln_i)
          IF(abs(ds) <= 3.0e3) THEN
            dNu = Fgr - v0sp(ln_i)
            bp = bp + Slabswint(dNu,v0sp(ln_i),x1p(ln_i),slabs1p(ln_i), &
                             &  yp(ln_i),yip(ln_i))
            dNu = Fgr - v0sm(ln_i)
            bm = bm + Slabswint(dNu,v0sm(ln_i),x1m(ln_i),slabs1m(ln_i), &
                             &  ym(ln_i),yim(ln_i))
          ENDIF
        end do
      ENDIF

      bv = beta_value
      ds = Log(bp/bv)/Log(tp/Temp)     ! Estimate over [temp+10,temp]
      ra = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
      dw = Log(bv/bm)/Log(Temp/tm)     ! Estimate over [temp,temp-10]
      t_power = (ds + 2.0 * ra + dw) / 4.0  ! Weighted Average
!
    ENDIF
!
  End Subroutine Create_beta
end module CREATE_BETA_M
! $Log$
! Revision 2.4  2001/10/18 15:59:26  zvi
! Modification for speed
!
! Revision 2.3  2001/10/18 07:13:04  zvi
! Make routine more efficient for nl=0
!
! Revision 2.2  2001/10/17 18:19:10  zvi
! Fixing bug in computing t_power
!
! Revision 2.1  2001/10/16 15:07:08  zvi
! Continuum parameters are now part of Catalog
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.14.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2001/02/01 18:12:04  vsnyder
! Initial conversion to Fortran 90
