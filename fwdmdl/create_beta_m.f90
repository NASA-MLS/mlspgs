module CREATE_BETA_M
  use MLSCommon, only: I4, R8
  use ABS_CS_AIR_CONT_M, only: ABS_CS_AIR_CONT
  use ABS_CS_H2O_213G_CONT_M, only: ABS_CS_H2O_213G_CONT
  use ABS_CS_LIQ_H2O_M, only: ABS_CS_LIQ_H2O
  use SpectroscopyCatalog_m, only: Catalog_T, Lines
  use SLABS_SW_M, only: SLABSWINT, DVOIGT_SPECTRAL
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
  Subroutine Create_beta (Spectag, pressure, Temp, Fgr, nl, Catalog,     &
         &   v0s, x1,y, yi, slabs1, dx1_dv0, dy_dv0, dslabs1_dv0, v0sp,  &
         &   x1p, yp, yip, slabs1p, v0sm, x1m, ym,yim,slabs1m,beta_value,&
         &   t_power, dbeta_dw, dbeta_dn, dbeta_dnu0, Frq_Gap, temp_der, &
         &   spect_der,Ier)
!
!  For a given channel, frequency and height, compute beta_value function.
!  This routine should be called for primary and image seperately.
!
    Integer(i4), intent(in) :: SPECTAG, nl
    Logical, intent(in) :: temp_der, spect_der

    Real(r8), intent(in) :: Pressure, Temp, Fgr, Frq_Gap

    Real(r8), intent(in) :: x1(:),y(:),yi(:),slabs1(:),slabs1m(:), &
   &          dx1_dv0(:),dy_dv0(:),dslabs1_dv0(:),v0sp(:),x1p(:),  &
   &          v0s(:),yp(:),yip(:),slabs1p(:),v0sm(:),x1m(:),ym(:), &
   &          yim(:)
!
    Type(Catalog_T), INTENT(IN) :: Catalog

    Integer(i4), intent(out) :: Ier
    Real(r8), intent(out) :: BETA_VALUE, T_POWER, DBETA_DW
    Real(r8), intent(out) :: DBETA_DN, DBETA_DNU0
!
! -----  Parameters Declaration ----------------------------------------
!
    Real(r8), parameter :: F213 = 2.15e5_r8, DF213 = 1.5e4_r8
    Real(r8), Parameter :: TINY = epsilon(F213)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: LN_I, j

    Real(r8) :: w,s,wd,ra,dNu,tp,bp,tm,bm,bb,dw,dn,ds
!
!  Setup absorption coefficients function
!
    Ier = 0
!
    t_power = 0.0
    dbeta_dw = 0.0
    dbeta_dn = 0.0
    dbeta_dnu0 = 0.0
    beta_value = tiny
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
      beta_value = abs_cs_h2o_213g_cont(Temp,Pressure,s)
!     t_power = -3.67           ! See code ...
      wd = Temp + 10.0d0
      ra = abs_cs_h2o_213g_cont(wd,Pressure,s)
      t_power = Log(ra/beta_value)/Log(wd/Temp)
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
    do ln_i = 1, nl
!
      j = Catalog%lines(ln_i)
      w = Lines(j)%W
!
! Prepare the temperature weighted coefficients:
!
      dNu = Fgr - v0s(ln_i)

      if(Frq_Gap > 0.0  .and. abs(dNu) >= Frq_Gap) CYCLE

      Call dvoigt_spectral(dNu,v0s(ln_i),x1(ln_i),yi(ln_i),y(ln_i), &
     &       w,temp,slabs1(ln_i),dx1_dv0(ln_i),dy_dv0(ln_i),          &
     &       dslabs1_dv0(ln_i),bb,spect_der,dw,dn,ds)
!
      beta_value = beta_value + bb
!
      if(spect_der) then
        dbeta_dw = dbeta_dw + dw
        dbeta_dn = dbeta_dn + dn
        dbeta_dnu0 = dbeta_dnu0 + ds
      endif

      if(.not. temp_der) CYCLE
!
!  Find the temperatue power dependency now:
!
      tp = temp + 10.0
      dNu = Fgr - v0sp(ln_i)
      bp = Slabswint(dNu,v0sp(ln_i),x1p(ln_i),slabs1p(ln_i),yp(ln_i), &
     &               yip(ln_i))
!
      tm = temp - 10.0
      dNu = Fgr - v0sm(ln_i)
      bm = Slabswint(dNu,v0sm(ln_i),x1m(ln_i),slabs1m(ln_i),ym(ln_i), &
     &               yim(ln_i))
!
      ds = Log(bp/bb)/Log(tp/temp)     ! Estimate over [temp+10,temp]
      ra = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
      wd = Log(bb/bm)/Log(temp/tm)     ! Estimate over [temp,temp-10]
!
      t_power = t_power + (ds+wd+2.0*ra)/4.0  ! Weighted Average
!
    end do
!
    t_power = t_power / nl
!
    Return
  End Subroutine Create_beta
end module CREATE_BETA_M
! $Log$
! Revision 1.10  2001/05/14 23:16:31  zvi
! Added Freq. Gap test..
!
! Revision 1.9  2001/05/14 23:14:54  zvi
! Added Freq. Gap test..
!
! Revision 1.8  2001/04/05 21:58:47  zvi
! Implementing l2cf inputs for FilterShape & Spectroscopy instead of FMI
!
! Revision 1.7  2001/04/03 07:32:45  zvi
! Modify the spectral structure - eliminating sps_ from the names
!
! Revision 1.6  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.5  2001/02/19 22:20:40  zvi
! Latest modification: Conv/NoConv
!
! Revision 1.4  2001/02/19 22:14:21  zvi
!
! Revision 1.1  2001/02/01 18:12:04  vsnyder
! Initial conversion to Fortran 90
