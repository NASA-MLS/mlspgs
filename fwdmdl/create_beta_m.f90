module CREATE_BETA_M

  implicit NONE
  private
  public :: Create_beta

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$id: create_beta_m.f90,v 2.11 2002/07/29 21:41:48 bill Exp $"
  character ( len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------------
contains
! *****     Public Subroutine     **************************************
! ----------------------------------------------  Create_beta  ---------

  subroutine Create_beta ( Spectag, cont, pressure, Temp, Fgr, nl, pfaw, &
         &   v0s, x1,y, yi, slabs1, beta_value, dslabs1_dv0, v0sp, x1p,  &
         &   yp, yip, slabs1p, v0sm, x1m, ym, yim, slabs1m, t_power,     &
         &   dbeta_dw, dbeta_dn, dbeta_dv )

!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image seperately.

    use ABS_CS_CONT_M,    only: ABS_CS_CONT
    use ABS_CS_LIQ_H2O_M, only: ABS_CS_LIQ_H2O
    use ABS_CS_N2_CONT_M, only: ABS_CS_N2_CONT
    use ABS_CS_O2_CONT_M, only: ABS_CS_O2_CONT
    use MLSCommon, only: R8, RP, IP
    use Molecules, only: L_Air_Cont, L_Extinction, L_Liq_H2O, L_O2, Spec_tags
    use SLABS_SW_M, only: DVOIGT_SPECTRAL, VOIGT_LORENTZ, SLABSWINT, SLABS

! Inputs:
    integer(ip), intent(in) :: SPECTAG ! molecule id tag
    real(rp), intent(in) :: cont(:) ! continuum parameters
    real(rp), intent(in) :: pressure ! pressure in hPa
    real(rp), intent(in) :: temp ! temperature in K
    real(rp), intent(in) :: fgr ! frequency in MHz
    integer(ip), intent(in) :: nl ! no of lines
    real(rp), intent(in) :: pfaw(:) ! line widths
    real(r8), intent(in) :: v0s(:) ! pressure shifted line centers
    real(rp), intent(in) :: x1(:) ! Doppler width
    real(rp), intent(in) :: y(:) ! ratio Pressure to Doppler widths
    real(rp), intent(in) :: yi(:) ! Interference coefficients
    real(rp), intent(in) :: slabs1(:) ! strengths
! optional inputs for spectral derivatives
    real(rp), optional, intent(in) :: dslabs1_dv0(:) ! strength derivative
!                                wrt line position
! optional inputs for temperature derivatives
    real(r8), optional, intent(in) :: v0sp(:) ! pressure shifted line centers
    real(rp), optional, intent(in) :: x1p(:)! Doppler width
    real(rp), optional, intent(in) :: yp(:) ! ratio Pressure to Doppler widths
    real(rp), optional, intent(in) :: yip(:) ! Interference coefficients
    real(rp), optional, intent(in) :: slabs1p(:) ! strengths
    real(r8), optional, intent(in) :: v0sm(:) ! pressure shifted line centers
    real(rp), optional, intent(in) :: x1m(:)! Doppler width
    real(rp), optional, intent(in) :: ym(:) ! ratio Pressure to Doppler widths
    real(rp), optional, intent(in) :: yim(:) ! Interference coefficients
    real(rp), optional, intent(in) :: slabs1m(:) ! strengths
! outputs
    real(rp), intent(out) :: beta_value
! optional outputs
    real(rp), optional, intent(out) :: T_POWER ! for temperature derivative
    real(rp), optional, intent(out) :: DBETA_DW ! line width derivative
    real(rp), optional, intent(out) :: DBETA_DN ! temperature dependence deriv
    real(rp), optional, intent(out) :: DBETA_DV ! line position derivative

! -----     Local variables     ----------------------------------------

    integer(ip) :: LN_I

    real(rp) :: w, ra, dNu, tp, bp, tm, bm, bv, dw, dn, ds, dbdw, dbdn, dbdv

    bv = 0.0_rp
    bp = 0.0_rp
    bm = 0.0_rp
    beta_value = 0.0_rp
    tp = Temp + 10.0_rp
    tm = Temp - 10.0_rp

!  Setup absorption coefficients function
!  Now get the beta_value:

    if ( spectag == spec_tags(l_liq_h2o) ) then

!  Liquid water

      bv = abs_cs_liq_h2o(Fgr,Temp)
      if ( present(t_power) ) then
        bm = abs_cs_liq_h2o(Fgr,tm)
        bp = abs_cs_liq_h2o(Fgr,tp)
      end if

    else if ( spectag == spec_tags(l_air_cont) ) then

!  Dry air contribution (N2)

      bv = abs_cs_n2_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_n2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_n2_cont(cont,tp,Pressure,Fgr)
      end if

    else if ( spectag == spec_tags(l_extinction) ) then

!  EXTINCTN molecule

      beta_value = 1.0_rp
      if ( present(t_power)) t_power = 0.0_rp
      Return

    else if ( spectag == spec_tags(l_o2) ) then ! O2

      bv = abs_cs_o2_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_o2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_o2_cont(cont,tp,Pressure,Fgr)
      end if

    else

      bv = abs_cs_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_cont(cont,tp,Pressure,Fgr)
      end if

    end if

    beta_value = bv
    if ( nl < 1 ) then
      if ( present(t_power) ) then
        ds = Log(bp/bv)/Log(tp/Temp)     ! Estimate over [temp+10,temp]
        ra = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
        dw = Log(bv/bm)/Log(Temp/tm)     ! Estimate over [temp,temp-10]
        t_power = (ds + 2.0 * ra + dw) / 4.0  ! Weighted Average
      end if
      return
    end if

    if ( present(dbeta_dw) .or. present(dbeta_dn) .or. present(dbeta_dv) ) then

      dbdw = 0.0_rp
      dbdn = 0.0_rp
      dbdv = 0.0_rp

      do ln_i = 1, nl

        dNu = Fgr - v0s(ln_i)

        w = pfaw(ln_i)
        if ( abs(y(ln_i))+0.666666_rp*abs(x1(ln_i)*dNu) > 100.0_rp ) then
          call Voigt_Lorentz ( dNu, v0s(ln_i), x1(ln_i), yi(ln_i), &
            &  y(ln_i), w, Temp,slabs1(ln_i), bv, dslabs1_dv0(ln_i), dw, dn, ds )
        else
          call DVoigt_Spectral ( dNu, v0s(ln_i), x1(ln_i), yi(ln_i), y(ln_i), &
         &     w, Temp, slabs1(ln_i), bv, dslabs1_dv0(ln_i), dw, dn, ds )
        end if

        beta_value = beta_value + bv
        dbdw = dbdw + dw
        dbdn = dbdn + dn
        dbdv = dbdv + ds

      end do

      if ( present(dbeta_dw)) dbeta_dw = dbdw
      if ( present(dbeta_dn)) dbeta_dn = dbdn
      if ( present(dbeta_dv)) dbeta_dv = dbdv

    else                ! No derivatives required

      if ( maxval(ABS(yi)) < 1.0e-06_rp ) then
        do ln_i = 1, nl
          dNu = Fgr - v0s(ln_i)
          beta_value = beta_value  + &
            &  Slabs(dNu,v0s(ln_i),x1(ln_i),slabs1(ln_i),y(ln_i))
        end do
      else
        do ln_i = 1, nl
          dNu = Fgr - v0s(ln_i)
          beta_value = beta_value + &
            &  Slabswint(dNu,v0s(ln_i),x1(ln_i),slabs1(ln_i),y(ln_i),yi(ln_i))
        end do
      end if

    end if

    if ( present(t_power) ) then

!  Find the temperatue power dependency now:

      if ( maxval(abs(yi)) < 1.0e-06_rp ) then
        do ln_i = 1, nl
          dNu = Fgr - v0sp(ln_i)
          bp = bp + Slabs(dNu,v0sp(ln_i),x1p(ln_i),slabs1p(ln_i),yp(ln_i))
          dNu = Fgr - v0sm(ln_i)
          bm = bm + Slabs(dNu,v0sm(ln_i),x1m(ln_i),slabs1m(ln_i),ym(ln_i))
        end do
      else
        do ln_i = 1, nl
          dNu = Fgr - v0sp(ln_i)
          bp = bp + Slabswint(dNu,v0sp(ln_i),x1p(ln_i),slabs1p(ln_i), &
                           &  yp(ln_i),yip(ln_i))
          dNu = Fgr - v0sm(ln_i)
          bm = bm + Slabswint(dNu,v0sm(ln_i),x1m(ln_i),slabs1m(ln_i), &
                           &  ym(ln_i),yim(ln_i))
        end do
      end if

      bv = beta_value
      ds = Log(bp/bv)/Log(tp/Temp)     ! Estimate over [temp+10,temp]
      ra = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
      dw = Log(bv/bm)/Log(Temp/tm)     ! Estimate over [temp,temp-10]
      t_power = (ds + 2.0 * ra + dw) / 4.0  ! Weighted Average

    end if

  end Subroutine Create_beta
end module CREATE_BETA_M

! $Log$
! Revision 2.11  2002/07/29 21:41:48  bill
! no changes, just debugging
!
! Revision 2.10  2002/03/06 02:29:28  zvi
! Removing more limits on large dNu
!
! Revision 2.9  2002/02/28 07:16:41  zvi
! Removing limit on large dNu
!
! Revision 2.8  2002/02/27 07:04:10  zvi
! Fixing limit on large dNu
!
! Revision 2.7  2001/11/15 01:22:00  zvi
! Remove Extiction debug
!
! Revision 2.6  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.5  2001/10/18 16:01:37  zvi
! Fix a small bug
!
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
