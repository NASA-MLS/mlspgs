!{Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Chi_Angles_m

  implicit NONE

  private
  public :: Get_Chi_Angles

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ---------------------------------------------  Get_Chi_Angles  -----
  subroutine Get_Chi_Angles ( sc_geoc_alt, tan_index_refr, tan_ht, &
             & phi_tan, Req, elev_offset, ptg_angle, dx_dh, dh_dz, tan_dh_dt, &
             & tan_d2h_dhdt, dx_dt, d2x_dxdt )

  ! Set up array of pointing angles

    use Units, only: Deg2Rad, Ln10
    use MLSCommon, only: RP

  ! inputs

    real(rp), intent(in) :: SC_geoc_alt ! geocentric spacecraft(observer)
                                        ! altitude in km
    real(rp), intent(in) :: Tan_index_refr ! tangent index of refraction
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req
    real(rp), intent(in) :: Phi_tan     ! tangent orbit plane projected
                                        ! geodetic angle in radians
    real(rp), intent(in) :: Req         ! equivalent earth radius in km
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
!                                         positive is towards the earth,
!                                         negative is towards space.
    real(rp), intent(in) :: dh_dz       ! dh/dz  at the tangent point

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
    Np1 = 1.0_rp + tan_index_refr

    !{ Empirical formula: $H_s = R_s + 38.9014\, \sin 2(\phi_t-51^{\circ}.6814 )$

    hs = sc_geoc_alt + ampl * sin(2.0*(phi_tan-phas))
    x = min ( Np1 * ht / hs, 1.0_rp )

    !{ $\sin \chi^{\text{refr}}_{\text{eq}} = N_t \frac{H_t}{H_s}
    !    \frac{ \text{min} ( H_t, H^{\oplus} ) }{H^{\oplus}}$

    ! ptg_angle = asin(x * min(ht,Req)/Req) - elev_offset

    if ( tan_ht >= 0.0_rp ) then      ! min(ht,Req)/Req = 1

      !{ $H_t \geq H^{\oplus}: ~~
      !  \sin \chi^{\text{refr}}_{\text{eq}} = N_t \frac{ H_t }{ H_s }; ~~
      !  \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }{ \text{d} H_t}
      !    \cos \chi^{\text{refr}}_{\text{eq}} =\frac{1}{H_s} \left ( N_t +
      !    \frac{ \text{d} N_t }{ \text{d} \zeta_t }
      !      \frac{ \text{d} \zeta_t }{ \text{d} H_t } H_t \right )$ \\ ~\\
      !  $N_t = 1 + a\, e^{-\zeta_t \ln 10} ~\Rightarrow~ 
      !    \frac{ \text{d} N_t }{ \text{d} \zeta_t } = - \ln 10 ( N_t - 1 )$

      sinChi = x
      q = Np1 - ht * tan_index_refr * Ln10 / dh_dz
    else                              ! min(ht,Req)/Req = ht/Req

      !{ $H_t < H^{\oplus}: ~~
      !  \sin \chi^{\text{refr}}_{\text{eq}} = 
      !    N_t \frac{ H_t^2 }{ H_s H^{\oplus} }; ~~
      !  \frac{ \text{d} N_t }{ \text{d} H_t} = 0 ~~
      !  \Rightarrow ~~
      !  \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }{ \text{d} H_t}
      !    \cos \chi^{\text{refr}}_{\text{eq}} =
      !    2 \frac{ N_t H_t }{ H_s H^{\oplus} }$

      sinChi = x*ht/Req
      q = 2.0_rp * ht * Np1 / Req
    end if

    ptg_angle = ASIN(sinChi) - elev_offset
    if ( sinChi < 1.0_rp ) then
      dx_dh = q / (hs * sqrt(1.0_rp - sinChi**2))
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

    return

  end subroutine Get_Chi_Angles

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Chi_Angles_m
! $Log$
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
