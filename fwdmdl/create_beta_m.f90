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
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
! *****     Public Subroutine     **************************************
! ----------------------------------------------  Create_beta  ---------

  subroutine Create_beta ( Spectag, cont, pressure, Temp, Fgr, pfaw,             &
         &   slabs_0, beta_value, Incl_Cld, cld_beta, slabs_p, slabs_m, t_power, &
         &   dbeta_dw, dbeta_dn, dbeta_dv )

!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image separately.

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use Molecules, only: SP_Air_Cont, SP_Extinction, SP_Liq_H2O, SP_O2
    use SLABS_SW_M, only: DVOIGT_SPECTRAL, VOIGT_LORENTZ, SLABSWINT, SLABS
    use Cloud_Extinction, only: get_beta_cloud

! Inputs:
    integer(ip), intent(in) :: SPECTAG ! molecule id tag
    real(rp), intent(in) :: cont(:) ! continuum parameters
    real(rp), intent(in) :: pressure ! pressure in hPa
    real(rp), intent(in) :: temp ! temperature in K
    real(rp), intent(in) :: fgr ! frequency in MHz
    real(rp), intent(in) :: pfaw(:) ! line widths
    type(slabs_struct), intent(in) :: slabs_0 ! contains, among others:
    Logical :: Incl_Cld


!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths
!    dslabs1_dv0(:) ! strength derivative wrt line position
! optional inputs for temperature derivatives
    type(slabs_struct), intent(in), optional :: slabs_p, slabs_m
! outputs
    real(rp), intent(out) :: beta_value, cld_beta
! optional outputs
    real(rp), optional, intent(out) :: T_POWER ! for temperature derivative
    real(rp), optional, intent(out) :: DBETA_DW ! line width derivative
    real(rp), optional, intent(out) :: DBETA_DN ! temperature dependence deriv
    real(rp), optional, intent(out) :: DBETA_DV ! line position derivative

! -----     Local variables     ----------------------------------------

    integer(ip) :: LN_I
    integer(ip) :: NL ! no of lines

    real(rp) :: ra, dNu, tp, bp, tm, bm, bv, dw, dn, ds, dbdw, dbdn, dbdv
    include 'constants.f9h'
    REAL(rp) :: WC(N)
    REAL(rp) :: W0(N)       ! SINGLE SCATTERING ALBEDO
    REAL(rp) :: PHH(N,NU)   ! PHASE FUNCTION
!----------------------------------------------------------------------------

    nl = size(pfaw)

    tp = Temp + 10.0_rp
    tm = Temp - 10.0_rp

!  Setup absorption coefficients function
!  Now get the beta_value:

    IF ( .not. Incl_Cld ) THEN

       cld_beta =0._rp   
  
    ELSE
       !  WC, IPSD, N, NU, NUA, NAB, NR needed to be inputed from l2cf
       !  will input WC later      
       call get_beta_cloud (Fgr, Temp, Pressure,               &
                          &  WC, IPSD, N, NU, NUA, NAB, NR,    &
                          &  cld_beta, W0, PHH                 )      
    ENDIF

    if ( spectag == sp_liq_h2o ) then ! ..................  Liquid Water

      beta_value = abs_cs_liq_h2o(Fgr,Temp)
      if ( present(t_power) ) then
        bm = abs_cs_liq_h2o(Fgr,tm)
        bp = abs_cs_liq_h2o(Fgr,tp)
      end if

    else if ( spectag == sp_air_cont ) then ! .................  Dry Air

      beta_value = abs_cs_n2_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_n2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_n2_cont(cont,tp,Pressure,Fgr)
      end if

    else if ( spectag == sp_extinction ) then ! ............  Extinction

      beta_value = 1.0_rp
      if ( present(t_power)) t_power = 0.0_rp
      return

    else if ( spectag == sp_o2 ) then ! ............................  O2

      beta_value = abs_cs_o2_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_o2_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_o2_cont(cont,tp,Pressure,Fgr)
      end if

    else ! ......................................................  Other

      beta_value = abs_cs_cont(cont,Temp,Pressure,Fgr)
      if ( present(t_power) ) then
        bm = abs_cs_cont(cont,tm,Pressure,Fgr)
        bp = abs_cs_cont(cont,tp,Pressure,Fgr)
      end if

    end if

    if ( nl < 1 ) then
      if ( present(t_power) ) then
        ds = log(bp/beta_value)/log(tp/temp)  ! Estimate over [temp+10,temp]
        ra = log(bp/bm)/        log(tp/tm)    ! Estimate over [temp+10,temp-10]
        dw = log(beta_value/bm)/log(temp/tm)  ! Estimate over [temp,temp-10]
        t_power = 0.25 * (ds + 2.0 * ra + dw) ! Weighted Average
      end if
      return
    end if

    if ( present(dbeta_dw) .or. present(dbeta_dn) .or. present(dbeta_dv) ) then

      dbdw = 0.0_rp
      dbdn = 0.0_rp
      dbdv = 0.0_rp

      do ln_i = 1, nl

        dNu = Fgr - slabs_0%v0s(ln_i)

        if ( abs(slabs_0%y(ln_i))+0.666666_rp*abs(slabs_0%x1(ln_i)*dNu) > 100.0_rp ) then
          call Voigt_Lorentz ( dNu, slabs_0%v0s(ln_i), slabs_0%x1(ln_i), &
            &  slabs_0%yi(ln_i), slabs_0%y(ln_i), pfaw(ln_i), Temp, &
            &  slabs_0%slabs1(ln_i), bv, slabs_0%dslabs1_dv0(ln_i), dw, dn, ds )
        else
          call DVoigt_Spectral ( dNu, slabs_0%v0s(ln_i), slabs_0%x1(ln_i), &
            &  slabs_0%yi(ln_i), slabs_0%y(ln_i), pfaw(ln_i), Temp, &
            &  slabs_0%slabs1(ln_i), bv, slabs_0%dslabs1_dv0(ln_i), dw, dn, ds )
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

      if ( maxval(ABS(slabs_0%yi)) < 1.0e-06_rp ) then
        do ln_i = 1, nl
          beta_value = beta_value + &
            &  Slabs(Fgr - slabs_0%v0s(ln_i), slabs_0%v0s(ln_i), &
            &        slabs_0%x1(ln_i), slabs_0%slabs1(ln_i), slabs_0%y(ln_i))
        end do
      else
        do ln_i = 1, nl
          beta_value = beta_value + &
            &  Slabswint(Fgr - slabs_0%v0s(ln_i), slabs_0%v0s(ln_i), &
            &            slabs_0%x1(ln_i), slabs_0%slabs1(ln_i), &
            &            slabs_0%y(ln_i), slabs_0%yi(ln_i))
        end do
      end if

    end if

    if ( present(t_power) ) then

!  Find the temperature power dependency now:

      if ( maxval(abs(slabs_0%yi)) < 1.0e-6_rp ) then
        do ln_i = 1, nl
          bp = bp + Slabs(Fgr - slabs_p%v0s(ln_i), slabs_p%v0s(ln_i), &
            &             slabs_p%x1(ln_i), slabs_p%slabs1(ln_i),slabs_p%y(ln_i))
          bm = bm + Slabs(Fgr - slabs_m%v0s(ln_i), slabs_m%v0s(ln_i), &
            &             slabs_m%x1(ln_i), slabs_m%slabs1(ln_i),slabs_m%y(ln_i))
        end do
      else
        do ln_i = 1, nl
          bp = bp + Slabswint(Fgr - slabs_p%v0s(ln_i), slabs_p%v0s(ln_i), &
            &                 slabs_p%x1(ln_i), slabs_p%slabs1(ln_i), &
            &                 slabs_p%y(ln_i), slabs_p%yi(ln_i))
          bm = bm + Slabswint(Fgr - slabs_m%v0s(ln_i), slabs_m%v0s(ln_i), &
            &                 slabs_m%x1(ln_i), slabs_m%slabs1(ln_i), &
            &                 slabs_m%y(ln_i), slabs_m%yi(ln_i))
        end do
      end if

      ds = Log(bp/beta_value)/Log(tp/Temp)  ! Estimate over [temp+10,temp]
      ra = Log(bp/bm)/        Log(tp/tm)    ! Estimate over [temp+10,temp-10]
      dw = Log(beta_value/bm)/Log(Temp/tm)  ! Estimate over [temp,temp-10]
      t_power = 0.25 * (ds + 2.0 * ra + dw) ! Weighted Average

    end if

  contains ! ===============================  Internal procedures  =====

    ! ----------------------------------------------  Abs_CS_Cont  -----

    ! Compute the general continuum contribution
    pure real(rp) function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency )
    ! real(rp) function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency )

      real(rp), intent(in) :: CONT(:)     ! continuum parameters
      real(rp), intent(in) :: TEMPERATURE ! in Kelvin
      real(rp), intent(in) :: PRESSURE    ! in mbar
      real(rp), intent(in) :: FREQUENCY   ! in MegaHertz

      abs_cs_cont = cont(1) * pressure * pressure * frequency * frequency * &
        & ( (300.0_rp / temperature)**cont(2) )

    end function Abs_CS_Cont

    ! -------------------------------------------  Abs_CS_Liq_H2O  -----

    ! Compute the liquid water correction
    pure real(rp) function Abs_CS_Liq_H2O ( Frequency, Temperature )
    ! real(rp) function ABS_CS_LIQ_H2O ( Frequency, Temperature )

      real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
      real(rp), intent(in) :: TEMPERATURE ! in Kelvin

    ! This function when multiplied by mass density (gm/m^3) of liquid droplet
    ! water gives absorption in Km^-1. Function comes from Liebe 1985 radio
    ! science paper and others.

      real(rp) :: TAU, EPSILON, THETA

      theta = 300.0_rp / temperature
      tau = 4.17e-8_rp * frequency * theta * exp(7.13_rp * theta)
      epsilon = (185.0_rp - 113.0_rp/theta) / (1.0_rp + tau * tau)
      abs_cs_liq_h2o = 1.886e-4_rp * frequency * tau * epsilon / &
                     & ((6.9_rp + epsilon)**2 + (tau*epsilon)**2)

    end function Abs_CS_Liq_H2O

    ! -------------------------------------------  Abs_CS_N2_Cont  -----

    ! Compute the N2 continuum contribution
    pure real(rp) function Abs_CS_N2_Cont ( Cont, Temperature, Pressure, Frequency )
    ! real(rp) Function Abs_CS_N2_cont ( Cont, Temperature, Pressure, Frequency )

      real(rp), intent(in) :: CONT(:)     ! continuum parameters
      real(rp), intent(in) :: TEMPERATURE ! in Kelvin
      real(rp), intent(in) :: PRESSURE    ! in mbar
      real(rp), intent(in) :: FREQUENCY   ! in MegaHertz

      real(rp) :: THETA, FSQR, FSXT

      theta = 300.0_rp / temperature
      fsqr = frequency * frequency
      fsxt = fsqr * theta
      abs_cs_n2_cont = pressure * pressure * fsqr * (theta**cont(2)) * &
                     & ( cont(1) * exp(-cont(3) * fsxt * theta) + &
                     &   cont(4) * exp(-cont(5) * fsxt * theta) * &
                     & (cont(6)**2 + fsqr))

    end function Abs_CS_N2_Cont

    ! -------------------------------------------  Abs_CS_O2_Cont  -----

    ! Compute the O2 continuum contribution
    pure real(rp) function Abs_CS_O2_Cont ( Cont, Temperature, Pressure, Frequency )
    ! real(rp) Function ABS_CS_O2_CONT ( Cont, Temperature, Pressure, Frequency )

      real(rp), intent(in) :: CONT(:)     ! continuum parameters
      real(rp), intent(in) :: TEMPERATURE ! in Kelvin
      real(rp), intent(in) :: PRESSURE    ! in mbar
      real(rp), intent(in) :: FREQUENCY   ! in MegaHertz

      real(rp) :: THETA, FSQR

      theta = 300.0_rp / temperature
      fsqr = frequency * frequency
      abs_cs_o2_cont = cont(1) * pressure * pressure * fsqr * (theta**cont(2)) &
                     & / (fsqr + (cont(3) * pressure * (theta**cont(4)) )**2 )

    end function Abs_CS_O2_Cont

  end Subroutine Create_beta

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CREATE_BETA_M

! $Log$
! Revision 2.17  2003/01/31 17:16:28  jonathan
! add Inc_Cld, and cld_ext
!
! Revision 2.16  2002/12/13 02:06:50  vsnyder
! Use a SLABS structure for the slabs quantities
!
! Revision 2.15  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.14  2002/09/24 23:16:48  vsnyder
! Fix up some comments
!
! Revision 2.13  2002/09/24 00:49:19  vsnyder
! Move Abs_CS_... to be internal procedures, some cosmetic changes
!
! Revision 2.12  2002/09/12 23:00:04  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
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
