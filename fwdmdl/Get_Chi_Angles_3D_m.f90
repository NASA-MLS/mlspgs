! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Chi_Angles_3D_m

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
  subroutine Get_Chi_Angles_All ( Tan_Refr_Index, Inst_Refr_Index, &
             & InstECR, LOS, Tan_Ht, Req, Elev_Offset, Ptg_Angle, &
             & dX_dH, dH_dZ, Tan_dH_dT, Tan_d2H_dHdT, dX_dT, d2X_dXdT )

  ! Compute pointing angle and its derivative with respect to tangent
  ! height. Optionally compute derivatives with respect to temperature.

    use Constants, only: Deg2Rad, Ln10
    use Geolocation_0, only: RG ! Kind for ECR components
    use Geometry, only: XYZ_to_Geod
    use MLSKinds, only: RP
    use VectorsModule, only: RV ! Kind for LOS, ScECR

  ! inputs

    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rv), intent(in) :: InstECR(3)  ! Instrument position in ECR
    real(rv), intent(in) :: LOS(3)      ! Line-of-sight unit vector in ECR
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req in km
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
                                        ! tangent wrt temperature, 
                                        ! pointing angle

  ! Local variables

    real(rp) :: Geod(3)      ! Geodetic coordinates of instrument (lat,lon,ht)
    real(rp) :: Grad(3)      ! Geodetic normal at instrument position
    real(rp) :: Grad_Los     ! dot_product(grad,los)**2
    real(rp) :: Ht           ! Tangent height from center of equivalent
                             ! spherical earth
    real(rp) :: Ns, Nt       ! Index of refraction at instrument, tangent
    real(rp) :: SinChi       ! Sin \chi_eq^refr
    real(rp) :: TanChi       ! Tan \chi_eq^refr

    geod = xyz_to_geod ( instECR )
    grad = [ cos(geod(1)) * cos(geod(2)), cos(geod(1)) * sin(geod(2)), &
           & sin(geod(1)) ]  ! Equation (18) in wvs-146.
    ht = tan_ht + req
    ns = Inst_Refr_Index + 1.0
    nt = Tan_Refr_Index + 1

  !{ The refraction-corrected pointing angle $\chi_\text{eq}^\text{refr}$ is
  !  given by Equation (14) in wvs-157:
  !
  !  $H_t \geq R_{\text{eq}}$:
  !  \begin{equation*}
  !    \sin \chi_\text{eq}^\text{refr}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t}{H_s}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \sin \chi_\text{eq}
  !  \end{equation*}
  !
  !  where $\chi_\text{eq}$ is the angle between the geodetic normal at the
  !  instrument position and the line-of-sight =
  !  $\sin^{-1} \frac{H_t}{H_s} =$
  !  $\cos^{-1}$( {\tt grad} $\cdot$ {\tt LOS} ).
  !  The latter does not require converting the instrument ECR position to
  !  geodetic co\"ordinates and height, then computing the radius of curvature
  !  at the geodetic co\"ordinates of the instrument, to get $H_s$.
  !
  !  $H_t < R_{\text{eq}}$:
  !  \begin{equation*}
  !    \sin \chi_\text{eq}^\text{refr}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t^2}{H_s R_{\text{eq}}}
  !      = \frac{H_t}{R_{\text{eq}}} \frac{\mathcal{N}_t}{\mathcal{N}_s}
  !        \sin \chi_\text{eq}
  !  \end{equation*}

    grad_los = dot_product(grad,los) ** 2 ! = cos^2 chi_eq (unrefracted)
    sinChi = ( nt / ns ) * sqrt( 1 - grad_los ) ! refracted
    if ( tan_ht < 0.0 ) sinChi = sinChi * ht/req
    ptg_angle = asin(sinChi) - elev_offset

  !{ $H_t \geq R_{\text{eq}}$:
  !  The derivative of the refraction-corrected pointing angle with respect to
  !  tangent height is given by Equation (17) in wvs-157:
  !
  !  \begin{equation*}
  !  \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }
  !                   { \text{d} H_t}
  !    \cos \chi^{\text{refr}}_{\text{eq}} =
  !  \sin \chi^{\text{refr}}_{\text{eq}} 
  !   \left[ \frac1{H_t} - \ln 10 \, \frac{\mathcal{N}_t-1}{\mathcal{N}_t}
  !    \frac{ \text{d} \zeta_t }{ \text{d} H_t} \right].
  !  \end{equation*}
  !
  !  $H_t < R_{\text{eq}}$:
  !  $\frac{\text{d} \mathcal{N}_t}{\text{d} H_t} = 0
  !  \Rightarrow
  !  \frac{ \text{d} \chi^{\text{refr}}_{\text{eq}} }
  !      { \text{d} H_t}\cos \chi^{\text{refr}}_{\text{eq}} =
  !  2 \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t}{H_s R_{\text{eq}}} =
  !  \frac2{H_t} \sin \chi^{\text{refr}}_{\text{eq}}$
  
    if ( sinChi < 1.0 ) then
      tanChi = sinChi / sqrt( ( 1.0 - sinChi ) * ( 1 + sinChi ) )
      if ( tan_ht >= 0.0 ) then
        dx_dh = tanChi * &
              & ( grad_los / ht - ln10 * tan_refr_index / nt / &
                           & dh_dz ) ! Equation (16) in wvs-157
      else
        dx_dh = tanChi * 2.0 / ht
      end if
    else
      dx_dh = 0.0
    end if

  ! Do temperature stuff if user requests it
  ! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
  ! (NOTE: These entities have NO PHI dimension, so take the center Phi in dh_dt)

    if ( present(tan_dh_dt) ) then
      block
        real(rp) :: TP
        tp = tan(ptg_angle)
        dx_dt = tp * tan_dh_dt / ht
        d2x_dxdt = tp * dx_dt + tan_d2h_dhdt
      end block
    end if

  end subroutine Get_Chi_Angles_All

  ! --------------------------------------  Get_Chi_Angles_Simple  -----
  subroutine Get_Chi_Angles_Simple ( Tan_Refr_Index, Inst_Refr_Index, &
             & InstECR, LOS, Tan_Ht, Req, Elev_Offset, Ptg_Angle )

  ! Compute pointing angle.

    use Constants, only: Deg2Rad, Ln10
    use Geolocation_0, only: RG ! Kind for ECR components
    use Geometry, only: XYZ_to_Geod
    use MLSKinds, only: RP
    use VectorsModule, only: RV ! Kind for LOS, ScECR

  ! inputs

    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rv), intent(in) :: InstECR(3)  ! Instrument position in ECR
    real(rv), intent(in) :: LOS(3)      ! Line-of-sight unit vector in ECR
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req in km
    real(rp), intent(in) :: Req         ! equivalent earth radius in km at tangent
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
                                        ! positive is towards the earth,
                                        ! negative is towards space.

  ! outputs

    real(rp), intent(out) :: Ptg_angle  ! pointing angle in radians

  ! Local variables

    real(rp) :: Geod(3)      ! Geodetic coordinates of instrument (lat,lon,ht)
    real(rp) :: Grad(3)      ! Geodetic normal at instrument position
    real(rp) :: Ht           ! Tangent height from center of equivalent
                             ! spherical earth
    real(rp) :: Ns, Nt       ! Index of refraction at instrument, tangent
    real(rp) :: SinChi       ! Sin \chi_eq^refr

    geod = xyz_to_geod ( instECR )
    grad = [ cos(geod(1)) * cos(geod(2)), cos(geod(1)) * sin(geod(2)), &
           & sin(geod(1)) ]  ! Equation (18) in wvs-146.
    ht = tan_ht + req
    ns = Inst_Refr_Index + 1.0
    nt = Tan_Refr_Index + 1

  !{ The refraction-corrected pointing angle $\chi_\text{eq}^\text{refr}$ is
  !  given by Equation (14) in wvs-157:
  !
  !  $H_t \geq R_{\text{eq}}$:
  !  \begin{equation*}
  !    \sin \chi_\text{eq}^\text{refr}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t}{H_s}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \sin \chi_\text{eq}
  !  \end{equation*}
  !
  !  where $\chi_\text{eq}$ is the angle between the geodetic normal at the
  !  instrument position and the line-of-sight =
  !  $\sin^{-1} \frac{H_t}{H_s} =$
  !  $\cos^{-1}$( {\tt grad} $\cdot$ {\tt LOS} ).
  !  The latter does not require converting the instrument ECR position to
  !  geodetic co\"ordinates and height, then computing the radius of curvature
  !  at the geodetic co\"ordinates of the instrument, to get $H_s$.
  !
  !  $H_t < R_{\text{eq}}$:
  !  \begin{equation*}
  !    \sin \chi_\text{eq}^\text{refr}
  !      = \frac{\mathcal{N}_t}{\mathcal{N}_s} \frac{H_t^2}{H_s R_{\text{eq}}}
  !      = \frac{H_t}{R_{\text{eq}}} \frac{\mathcal{N}_t}{\mathcal{N}_s}
  !        \sin \chi_\text{eq}
  !  \end{equation*}

    sinChi = ( ( nt / ns ) * &
             & sqrt( 1 - dot_product(grad,los) ** 2 ) ) ! 1 - cos^2 \chi_eq
    if ( tan_ht < 0.0 ) sinChi = sinChi * ht/req
    ptg_angle = asin(sinChi) - elev_offset

  end subroutine Get_Chi_Angles_Simple

  ! --------------------------------  Get_Chi_Angles_Simple_Deriv  -----
  subroutine Get_Chi_Angles_Simple_Deriv ( Tan_Refr_Index, Inst_Refr_Index, &
             & InstECR, LOS, Tan_Ht, Req, Elev_Offset, Ptg_Angle, &
             & Tan_dH_dT, Tan_d2H_dHdT, dX_dT, d2X_dXdT )

  ! Compute pointing angle and its derivatives with respect to temperature

    use Constants, only: Deg2Rad, Ln10
    use Geolocation_0, only: RG ! Kind for ECR components
    use Geometry, only: XYZ_to_Geod
    use MLSKinds, only: RP
    use VectorsModule, only: RV ! Kind for LOS, ScECR

  ! inputs

    real(rp), intent(in) :: Tan_Refr_Index ! Index of refraction - 1 at tangent
    real(rp), intent(in) :: Inst_Refr_Index ! Index of refraction - 1 at instrument
    real(rv), intent(in) :: InstECR(3)  ! Instrument position in ECR
    real(rv), intent(in) :: LOS(3)      ! Line-of-sight unit vector in ECR
    real(rp), intent(in) :: Tan_ht      ! tangent height relative to Req in km
    real(rp), intent(in) :: Req         ! equivalent earth radius in km at tangent
    real(rp), intent(in) :: Elev_offset ! radiometer pointing offset in radians
                                        ! positive is towards the earth,
                                        ! negative is towards space.
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
  !  Local variables:

    real(rp) :: tp

    call Get_Chi_Angles_Simple ( Tan_Refr_Index, Inst_Refr_Index, &
             & InstECR, LOS, Tan_Ht, Req, Elev_Offset, Ptg_Angle )

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

end module Get_Chi_Angles_3D_m
! $Log$
! Revision 2.3  2020/08/21 19:56:09  vsnyder
! Account for d (f.g) / d H_t
!
! Revision 2.2  2020/08/17 22:45:11  vsnyder
! Correct some mistakes. Correct some comments. Handle subsurface tangent.
!
! Revision 2.1  2020/08/12 00:04:10  vsnyder
! Initial commit
!
