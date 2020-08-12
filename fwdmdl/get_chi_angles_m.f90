! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Chi_Angles_m

  implicit NONE

  private
  public :: Get_Chi_Angles

  interface Get_Chi_Angles
    module procedure Get_Chi_Angles_All, Get_Chi_Angles_Simple
    module procedure Get_Chi_Angles_Simple_Deriv
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! -----------------------------------------  Get_Chi_Angles_All  -----
  subroutine Get_Chi_Angles_All ( SC_Geoc_Alt, Tan_Refr_Index, Inst_Refr_Index, &
             & Tan_Ht, Phi_Tan, Req, Elev_Offset, Ptg_Angle, dX_dH, dH_dZ, &
             & Tan_dH_dT, Tan_d2H_dHdT, dX_dT, d2X_dXdT )

  ! Compute pointing angle and its derivative with respect to tangent
  ! height. Optionally compute derivatives with respect to temperature.

    use Constants, only: Deg2Rad, Ln10
    use MLSKinds, only: RP

  ! inputs

    real(rp), intent(in) :: SC_geoc_alt ! geocentric spacecraft(observer)
                                        ! altitude in km
    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req in km
    real(rp), intent(in) :: Phi_tan     ! tangent orbit plane projected
                                        ! geodetic angle in radians
    real(rp), intent(in) :: Req         ! equivalent earth radius in km at tangent
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
                                        ! positive is towards the earth,
                                        ! negative is towards space.
    real(rp), intent(in) :: dh_dz       ! dh/dz at the tangent point

  ! outputs

    real(rp), intent(out) :: Ptg_angle  ! pointing angle in radians
    real(rp), intent(out) :: dx_dh      ! derivative of pointing angle
                                        ! wrt height

  ! keywords

    real(rp), optional, intent(in) :: Tan_dh_dt(:) ! derivative of tangent
                                                   ! height wrt temperature
    real(rp), optional, intent(in) :: Tan_d2h_dhdt(:) ! 2nd derivative of
                                        ! tangent height wrt temperature & height
    real(rp), optional, intent(out) :: dx_dt(:) ! derivative of pointing angle
                                                ! wrt temperature
    real(rp), optional, intent(out) :: d2x_dxdt(:) ! second derivative of
                                        ! tangent wrt temperature, pointing angle
  !  ----------------
  !  Local variables:
  !  ----------------

    real(rp), parameter :: ampl = 38.9014
    real(rp), parameter :: phas = 51.6814 * Deg2Rad

    real(rp) :: hs, ht, Np1, q, sinChi, tp, x

  ! Start code:

    ht = Req + tan_ht
    Np1 = ( 1.0_rp + tan_refr_index )         ! N_t

    !{ Empirical formula (8.5): $H_s = R_s + 38.9014\,
    !  \sin 2(\phi_t-51^\circ\hspace{-4.2pt}.\hspace{1pt}6814 )$\\
    !  This only works for EMLS at its current orbit inclination!

    hs = sc_geoc_alt + ampl * sin(2.0*(phi_tan-phas))
    x = min ( Np1 / ( 1.0_rp + inst_refr_index ) * ht / hs, 1.0_rp )

    !{ $\sin \chi^{\text{refr}}_{\text{eq}} =
    !    \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t}{H_s}
    !    \frac{ \text{min} ( H_t, H^{\oplus}_t ) }{H^{\oplus}_t}$,
    !  where $\mathcal{N}_s$ is the index of refraction at the spacecraft. 
    !  This is a bit different from Equation (8.3) in the 19 Aug 2004 ATBD
    !  because that equation does not consider the case of $H_t <
    !  H^{\oplus}_t$, and assumes the instrument is in orbit, and therefore
    !  $N_s = 1.0$.

    ! ptg_angle = asin(x * min(ht,Req)/Req) - elev_offset

    if ( tan_ht >= 0.0_rp ) then      ! min(ht,Req)/Req = 1

      !{ $H_t \geq H^{\oplus}_t: ~
      !  \sin \chi^{\text{refr}}_{\text{eq}} =
      !    \frac{N_t}{N_s} \frac{ H_t }{ H_s }; $ \\
      !  $ \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }{ \text{d} H_t}
      !    \cos \chi^{\text{refr}}_{\text{eq}} =
      !    \frac{1}{N_s H_s} \left ( N_t +
      !     \frac{\text{d} N_t}{\text{d} H_t } H_t \right) =
      !    \frac{1}{N_s H_s} \left ( N_t +
      !     \frac{ \text{d} N_t }{ \text{d} \zeta_t }
      !      \frac{ \text{d} \zeta_t }{ \text{d} H_t } H_t \right )$ \\ ~\\
      !  $N_t = 1 + a\, e^{-\zeta_t \ln 10} ~\Rightarrow~ 
      !    \frac{ \text{d} N_t }{ \text{d} \zeta_t } = - \ln 10 ( N_t - 1 )
      !  \Rightarrow $ \\
      !  $\frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }
      !                   { \text{d} H_t}
      !    \cos \chi^{\text{refr}}_{\text{eq}} =
      !  \sin \chi^{\text{refr}}_{\text{eq}} 
      !   \left[ \frac1{H_t} - \frac1{N_t} ( \ln 10 ( N_t-1 ) )
      !    \frac{ \text{d} \zeta_t }{ \text{d} H_t} \right] $

      sinChi = x
      q = 1.0_rp / ht - 1.0_rp/np1 * tan_refr_index * Ln10 / dh_dz
    else                              ! min(ht,Req)/Req = ht/Req

      !{ $H_t < H^{\oplus}_t: ~~
      !  \sin \chi^{\text{refr}}_{\text{eq}} = 
      !    \frac{N_t}{N_s} \frac{ H_t^2 }{ H_s H^{\oplus}_t }; ~~
      !  \frac{ \text{d} N_t }{ \text{d} H_t} = 0 ~~
      !  \Rightarrow ~~
      !  \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }{ \text{d} H_t}
      !    \cos \chi^{\text{refr}}_{\text{eq}} =
      !    2 \frac{N_t}{N_s} \frac{H_t }{H_s H^{\oplus}_t } =
      !    \frac2{H_t} \sin \chi^{\text{refr}}_{\text{eq}} $

      sinChi = x*ht/Req
      q = 2.0_rp / ht
    end if

    ptg_angle = ASIN(sinChi) - elev_offset
    if ( sinChi < 1.0_rp ) then
      dx_dh = q * sinChi / sqrt(1.0_rp - sinChi**2)
    else
      dx_dh = 0.0_rp
    end if

  ! Do temperature stuff if user requests it
  ! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
  ! (NOTE: These entities have NO PHI dimension, so take the center Phi in dh_dt)

    if ( present(tan_dh_dt) ) then
      tp = tan(ptg_angle)
      dx_dt = tp * tan_dh_dt / ht
      d2x_dxdt = tp * dx_dt + tan_d2h_dhdt
    end if

  end subroutine Get_Chi_Angles_All

  ! --------------------------------------  Get_Chi_Angles_Simple  -----
  subroutine Get_Chi_Angles_Simple ( SC_Geoc_Alt, Tan_Refr_Index, &
             & Inst_Refr_Index, Tan_Ht, Phi_Tan, Req, Elev_Offset, Ptg_Angle )

  ! Compute pointing angle.

    use Constants, only: Deg2Rad
    use MLSKinds, only: RP

  ! inputs

    real(rp), intent(in) :: SC_geoc_alt ! geocentric spacecraft(observer)
                                        ! altitude in km
    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req
    real(rp), intent(in) :: Phi_tan     ! tangent orbit plane projected
                                        ! geodetic angle in radians
    real(rp), intent(in) :: Req         ! equivalent earth radius in km
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
                                        ! positive is towards the earth,
                                        ! negative is towards space.

  ! outputs

    real(rp), intent(out) :: Ptg_angle  ! pointing angle in radians

  !  ----------------
  !  Local variables:
  !  ----------------

    real(rp), parameter :: ampl = 38.9014
    real(rp), parameter :: phas = 51.6814 * Deg2Rad

    real(rp) :: hs, ht, Np1, sinChi, x

  ! Start code:

    ht = Req + tan_ht
    Np1 = ( 1.0_rp + tan_refr_index ) / ( 1.0_rp + inst_refr_index )

    !{ Empirical formula (8.5): $H_s = R_s + 38.9014\,
    !  \sin 2(\phi_t-51^\circ\hspace{-4.2pt}.\hspace{1pt}6814 )$\\
    !  This only works for EMLS at its current orbit inclination!

    hs = sc_geoc_alt + ampl * sin(2.0*(phi_tan-phas))
    x = min ( Np1 * ht / hs, 1.0_rp )

    !{ $\sin \chi^{\text{refr}}_{\text{eq}} =
    !    \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t}{H_s}
    !    \frac{ \text{min} ( H_t, H^{\oplus}_t ) }{H^{\oplus}_t}$,
    !  where $\mathcal{N}_s$ is the index of refraction at the spacecraft. 
    !  This is a bit different from Equation (8.3) in the 19 Aug 2004 ATBD
    !  because that equation does not consider the case of $H_t <
    !  H^{\oplus}_t$, and assumes the instrument is orbit, and therefore
    !  $N_s = 1.0$.

    ! ptg_angle = asin(x * min(ht,Req)/Req) - elev_offset

    if ( tan_ht >= 0.0_rp ) then      ! min(ht,Req)/Req = 1

      !{ $H_t \geq H^{\oplus}_t: ~~
      !  \sin \chi^{\text{refr}}_{\text{eq}} = 
      !  \frac{N_t}{N_s} \frac{ H_t }{ H_s }$

      sinChi = x
    else                              ! min(ht,Req)/Req = ht/Req

      !{ $H_t < H^{\oplus}_t: ~~
      !  \sin \chi^{\text{refr}}_{\text{eq}} = 
      !    \frac{N_t}{N_s} \frac{ H_t^2 }{ H_s H^{\oplus}_t }$

      sinChi = x*ht/Req
    end if

    ptg_angle = ASIN(sinChi) - elev_offset

  end subroutine Get_Chi_Angles_Simple

  ! --------------------------------  Get_Chi_Angles_Simple_Deriv  -----
  subroutine Get_Chi_Angles_Simple_Deriv ( SC_Geoc_Alt, Tan_Refr_Index, &
             & Inst_Refr_Index, Tan_Ht, Phi_Tan, Req, Elev_Offset, Ptg_Angle, &
             & Tan_dH_dT, Tan_d2H_dHdT, dX_dT, d2X_dXdT )

  ! Compute Pointing angle and its derivative with respect to temperature.

    use Constants, only: Deg2Rad
    use MLSCommon, only: RP

  ! inputs

    real(rp), intent(in) :: SC_geoc_alt ! geocentric spacecraft(observer)
                                        ! altitude in km
    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req
    real(rp), intent(in) :: Phi_tan     ! tangent orbit plane projected
                                        ! geodetic angle in radians
    real(rp), intent(in) :: Req         ! equivalent earth radius in km
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
!                                         positive is towards the earth,
!                                         negative is towards space.
    real(rp), intent(in) :: Tan_dh_dt(:) ! derivative of tangent
                                         ! height wrt temperature
    real(rp), intent(in) :: Tan_d2h_dhdt(:) ! 2nd derivative of
                                        ! tangent height wrt temperature & height
  ! outputs

    real(rp), intent(out) :: Ptg_angle  ! pointing angle in radians
    real(rp), intent(out) :: dx_dt(:)   ! derivative of pointing angle
                                        ! wrt temperature
    real(rp), intent(out) :: d2x_dxdt(:) ! second derivative of
                                        ! tangent wrt temperature, pointing angle

  !  ----------------
  !  Local variables:
  !  ----------------

    real(rp) :: tp

  ! Start code:

   call Get_Chi_Angles_Simple ( SC_Geoc_Alt, Tan_Refr_Index, &
                              & Inst_Refr_Index, Tan_Ht, Phi_Tan, Req, &
                              & Elev_Offset, Ptg_Angle )

  ! Do temperature stuff
  ! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
  ! (NOTE: These entities have NO PHI dimension, so take the center Phi in dh_dt)

    tp = tan(ptg_angle)
    dx_dt = tp * tan_dh_dt /  ( Req + tan_ht )
    d2x_dxdt = tp * dx_dt + tan_d2h_dhdt

  end subroutine Get_Chi_Angles_Simple_Deriv

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Chi_Angles_m
! $Log$
! Revision 2.26  2020/04/22 01:59:43  vsnyder
! Some work on QTM chi angles
!
! Revision 2.25  2019/06/24 23:28:17  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.24  2019/02/28 01:58:18  vsnyder
! More TeXnicalities only
!
! Revision 2.23  2019/02/28 01:45:14  vsnyder
! TeXnicalities only
!
! Revision 2.22  2016/05/27 01:25:28  vsnyder
! TeXnicalities
!
! Revision 2.21  2013/06/12 02:25:43  vsnyder
! Cruft removal
!
! Revision 2.20  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.19  2009/05/13 20:03:01  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.18  2008/05/20 00:18:46  vsnyder
! Separate angles from angles-and-derivatives
!
! Revision 2.17  2005/12/07 00:32:58  vsnyder
! Add some TeXnicalities
!
! Revision 2.16  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.15  2004/05/17 22:05:40  livesey
! Minor changes to avoid explosions due to asin(>1)
!
! Revision 2.14  2003/05/14 22:23:10  bill
! corrected sign convention for elev offset angle
!
! Revision 2.13  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/10/03 22:23:03  vsnyder
! Fix a bug in derivative calculation.  Get Deg2Rad from Units.  Cosmetic
! changes, including LaTeX equations.
!
! Revision 2.11  2002/09/26 20:15:26  vsnyder
! Get Ln10 from Units module
!
! Revision 2.10  2002/09/26 00:48:45  vsnyder
! Insert copyright notice.  Move USEs from module scope to procedure scope.
! Cosmetic changes.  Get Ln10 from Geometry module.
!
! Revision 2.9  2002/07/05 07:52:47  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.8  2002/06/28 11:06:50  zvi
! computes dx_dh on output grid as well
!
! Revision 2.7  2002/06/24 21:11:24  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.3  2002/06/11 22:20:45  bill
! eliminate zero-out feature--wgr
!
! Revision 2.2  2002/06/07 04:50:36  bill
! fixes and improvements--wgr
!
! Revision 2.1  2001/10/16 22:32:11  zvi
! Correcting for Earth intersection case
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.10.2.3  2001/09/13 11:18:21  zvi
! use the correct eta
!
! Revision 1.10.2.2  2001/09/12 21:38:47  zvi
! Added CVS stuff
!
! Revision 1.10.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
