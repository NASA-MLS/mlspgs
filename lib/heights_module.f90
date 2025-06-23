module heights_module

! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

  ! This module provides conversions from geodetic (geometric height above
  ! earth surface) to geopotential height (above geoid).
  ! Formulae taken from  "Theoretical basis for EOS-MLS geopotential
  ! and scan calculations" by N. J. Livesey


  use Constants, only: Deg2Rad, Rad2Deg

  use Geometry, only: Earth_Major_Axis => EarthRadA, &
                    & Earth_Minor_Axis => EarthRadB, EarthSurfaceGPH, GM, G0, &
                    & GeodToGeocLat, J2, J4, Omega => W

  use MLSKinds, only: R8

  use MLSMessageModule, only: MLSMessage, MLSMSG_Warning

  implicit NONE
  private
  public:: Geom_to_GPH, GPH_to_geom

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------------  geom_to_gph  -----
  function geom_to_gph ( geom_height, lat_geoc, lat_geod ) result (gph)
    ! Converts geometric height in meters above the geoid to geopotential
    ! height in meters above the geoid.
    ! This applies Equation 40 in "Theoretical basis for EOS-MLS geopotential
    ! and scan calculations" by N. J. Livesey
    ! If I have got this right, the difference between the two heights is
    ! piffling at ground level, going to about 1km at 80km.

    !----Argument----!
    real(kind=r8), intent(in) :: geom_height ! above geoid
    real(kind=r8), intent(in), optional :: lat_geoc, lat_geod
    !----Function result----!
    real(kind=r8)::gph
    !---Local variabules---!
    real(kind=r8) :: p2, p4, sinlat, h_inf, ahrat, coslat, h0, hc, lat
    !---------Executable statements-------!

    if ( present(lat_geoc) .and. present(lat_geod) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
           "in function geom_to_gph,"// &
           "Both geometric and geodetic lats supplied: using geometric" )
      lat = lat_geoc
    else if ( present(lat_geoc) ) then
      lat = lat_geoc
    else if ( present(lat_geod) ) then
      lat = rad2deg * geodToGeocLat(lat_geod)
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
           "in function geom_to_gph, Called with no latitude - assuming 45deg N" )
      lat = 45.0_r8
    end if

    sinlat = sin(lat*deg2rad)**2
    coslat = 1.0_r8 - sinlat

    !{ Calculate height of Earth's surface at given geocentric latitude $\beta$
    !  \begin{equation*}
    !   H_0 = \sqrt{R_a^2 \cos^2\beta + R_b^2 \sin^2\beta}
    !  \end{equation*}
    !  where $R_a$ is the Earth's semi-major axis and $R_b$ is the
    !  semi-minor axis.

    h0 = sqrt( earth_major_axis**2*coslat + earth_minor_axis**2*sinlat )
    hc = h0 + geom_height ! requested height, from center of geoid

    !{ Coefficients in the spherical harmonic expansion of the Earth's
    !  figure: $P_2 = ( 3 \sin^2 \beta - 1 ) /2$ and
    !          $P_4 = ( 35 \sin^4 \beta - 30 \sin^2 \beta + 3 ) / 8$.
    p2 = (3*sinlat -1)/2
    p4 = (35*(sinlat**2) - 30*sinlat +3)/8

    !{ Let $r$ be the ratio of the Earth's major axis to the height
    !  measured from the center of the geoid, $r = R_a / H_c$
    ahrat = earth_major_axis / hc
    ahrat = ahrat*ahrat

    !{ Compute geopotential height, accounting for latitude, the Earth's
    !  figure up to fourth-order zonal harmonics, and centripetal acceleration:
    !  \begin{equation*}
    !   H_\infty = \frac{Gm}{g_0 H_c} \left(
    !    1 - J_2 P_2 r^2 - J_4 P_4 r^4 \right) +
    !     \frac{\omega^2 r^2 \cos^2 \beta}{2 g_0}
    !  \end{equation*}
    !  where $J_2$ and $J_4$ are zonal spherical harmonics.
    h_inf = (gm/(g0*hc))*(1-j2*p2*ahrat - j4*p4*ahrat*ahrat) + &
      &  (omega*omega*hc*hc*coslat)/(2*g0)

    gph = earthSurfaceGPH - h_inf

  end function geom_to_gph

  ! ------------------------------------------------  gph_to_geom  -----
  function gph_to_geom ( gph, lat_geoc, lat_geod ) result (geom_height)
    ! This converts Geopotential height (GPH) in meters above the geoid
    ! back to geometric height in meters above the geoid.
    ! We use a yucky secant-method iteration mainly out of idleness --
    ! we can't be bothered to invert the function analytically. Of course,
    ! the present hack has the advantage that gph_to_geom and geom_to_gph
    ! are always consistent with each other.
    ! Tests suggest that the secant method converges after either 3 or 4
    ! iterations over the 0-100km range.

    !----Argument----!
    real(kind=r8), intent(in) :: gph
    real(kind=r8), intent(in), optional :: lat_geoc, lat_geod
    !----Function result----!
    real(kind=r8)::geom_height
    !---Local variabules---!
    real(kind=r8),parameter:: accuracy=0.001
    real(kind=r8)::geom1,geom2,newgeom,dgph1,dgph2,lat
    !---------Executable statements-------!

    if ( present(lat_geoc) .and. present(lat_geod) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
           "in function gph_to_geom,"// &
           "Both geometric and geodetic lats supplied: using geometric" )
      lat = lat_geoc
    else if ( present(lat_geoc) ) then
      lat = lat_geoc
    else if ( present(lat_geod) ) then
      lat = rad2deg * geodToGeocLat(lat_geod)
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
           "in function gph_to_geom, heights_module"// &
           "Called with no latitude - assuming 45deg N" )
      lat = 45.0_r8
    end if

    geom1 = gph
    geom2 = gph+1.0
    dgph1 = geom_to_gph(geom1,lat)-gph
  secantloop: do
      dgph2 = geom_to_gph(geom2,lat)-gph
      newgeom = geom1 - dgph1*(geom2-geom1)/(dgph2-dgph1)
      if ( abs(newgeom-geom2) < accuracy ) then
         exit secantloop
      end if
      geom1 = geom2
      dgph1 = dgph2
      geom2 = newgeom
    end do secantloop
    geom_height = newgeom

  end function gph_to_geom

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module heights_module

! $Log$
! Revision 2.10  2015/03/28 01:57:03  vsnyder
! Deleted Lat_geod_to_geoc (its functionality is in Geometry).
!
! Revision 2.9  2013/05/02 00:17:48  vsnyder
! LaTeX corrections, but it doesn't matter because nobody uses this module
!
! Revision 2.8  2013/05/01 23:18:52  vsnyder
! LaTeX cannonball polishing
!
! Revision 2.7  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2009/05/13 20:40:05  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.5  2006/09/28 20:39:22  vsnyder
! Correct singularity in lat_geod_to_geoc
!
! Revision 2.4  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/09/26 20:56:53  vsnyder
! Move Earth_Axis_Ratio_Squared from heights_module to Geometry
!
! Revision 2.1  2002/09/26 20:54:43  vsnyder
! Get constants from Geometry and Units, cosmetic changes
!
