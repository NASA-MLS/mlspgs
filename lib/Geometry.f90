! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Geometry

  ! This module contains some geometry routines and constants common to the
  ! forward model and the scan model.

  use MLSKinds, only: R8, RP

  implicit NONE
  private

  ! Constants
  public :: Earth_Axis_Ratio_Squared, Earth_Axis_Ratio_Squared_m1 ! a^2/b^2, a^2/b^2-1
  public :: EarthRadA, EarthRadB, EarthSurfaceGPH
  public :: GM, G0, J2, J4, SecPerYear, W, MaxRefraction

  ! Procedures
  public :: GeocToGeodLat, GeodToGeocLat, Get_R_Eq, To_Cart

  ! Earth dimensions.

  real(r8), parameter :: EarthRadA = 6378137.0_r8    ! Major axis radius in m 
  real(r8), parameter :: EarthRadB = 6356752.3141_r8 ! Minor axis radius in m 
  real(r8), parameter :: Earth_Axis_Ratio_Squared = EarthRadA**2 / EarthRadB**2
  real(r8), parameter :: Earth_Axis_Ratio_Squared_m1 = &
    & (EarthRadA / EarthRadB - 1) * (EarthRadA / EarthRadB + 1)

  ! Gravity-related terms.

  ! IUGG (International Union of Geodesy and Geodynamics) 1999 values from
  ! http://www.gfy.ku.dk/~iag/Travaux_99/sc3.htm (part III)

  ! Geocentric Earth gravitational constant, including Earth's atmosphere

! real(rp), parameter :: GM = 3.986004418e14_rp ! m^3/sec^2 (IUGG 1999 value)
  real(rp), parameter :: GM = 3.98600436e14_rp  ! m^3/sec^2

  ! Mean equatorial gravity

  real(r8), parameter :: G0 = 9.80665           ! Nominal little g ms-2
! real(r8), parameter :: G0 = 9.7803278         ! IUGG 1999 value

  ! Stokes second- and fourth-degree zonal harmonics

  real(rp), parameter :: J2 = 0.0010826256_rp   ! 1980 reference geoid value
! real(rp), parameter :: J2 = 0.0010826359_rp   ! IUGG 1999 value
  real(rp), parameter :: J4 = -.0000023709122_rp  ! 1980 reference geoid value

  ! earth rotational velocity.

  real(rp), parameter :: W = 7.292115e-05_rp    ! rad/sec

  ! Earth surface geopotential height.

  real(r8), parameter :: EarthSurfaceGPH = 6387182.265_r8 ! meters

  ! Seconds per tropical year, 1994-1998, based on orbital elements by
  ! Laskar.  See http://scienceworld.wolfram.com/astronomy/TropicalYear.html

  real(r8), parameter :: SecPerYear = 365.242190_r8 * 86400.0_r8

  ! This is the maximum amount of refraction allowed
  real(r8), parameter :: MaxRefraction = 0.0004 ! Add one to get refractive index

  interface GeocToGeodLat
    module procedure GeocToGeodLat_D, GeocToGeodLat_S
  end interface

  interface GeodToGeocLat
    module procedure GeodToGeocLat_D, GeodToGeocLat_S
  end interface

  interface To_Cart
    module procedure To_Cart_D, To_Cart_S
  end interface

! *****     Private Constants     **************************************

! AQUAD   Square of major half axis for earth ellipsoid
! BQUAD   Square of minor half axis for earth ellipsoid

  real, parameter :: AQUAD = (EarthRadA/1000.0)**2 ! Km**2
  real, parameter :: BQUAD = (EarthRadB/1000.0)**2 ! Km**2

! ERAD    Earth radius for normalization of Cartesian
!         coordinates (6371.2 Km) as recommended by the International
!         Astronomical Union

  real, parameter :: ERAD = 6371.2

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------  GeocToGeodLat_D  -----
  double precision elemental function GeocToGeodLat_D ( geocLat ) result ( GeodLat )

  !{ Convert a geocentric latitude $\lambda$ to a geodetic one $\mu$, both
  !  in degrees.
  !  Use the relation $\mu = \tan^{-1} \frac{f^2 \sin\lambda}{\cos\lambda}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    use Constants, only: Deg2Rad, Rad2Deg

    double precision, intent(in) :: geocLat

    double precision :: Lat ! Geocentric latitude in radians
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( f2 * sin(lat), cos(lat) )

  end function GeocToGeodLat_D

  ! --------------------------------------------  GeocToGeodLat_S  -----
  real elemental function GeocToGeodLat_S ( geocLat ) result ( GeodLat )

  !{ Convert a geocentric latitude $\lambda$ to a geodetic one $\mu$. both
  !  in degrees.
  !  Use the relation $\mu = \tan^{-1} \frac{f^2 \sin\lambda}{\cos\lambda}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    use Constants, only: Deg2Rad, Rad2Deg

    real, intent(in) :: geocLat

    real :: Lat ! Geocentric latitude in radians
    real, parameter :: F2 = Earth_Axis_Ratio_Squared

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( f2 * sin(lat), cos(lat) )

  end function GeocToGeodLat_S

  ! --------------------------------------------  GeodToGeocLat_D  -----
  double precision elemental function GeodToGeocLat_D ( GeodLat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ (IN DEGREES!) into a geocentric one
  !  $\lambda$ (IN RADIANS!)
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    use Constants, only: Deg2Rad

    ! Arguments
    double precision, intent(in) :: GeodLat

    double precision :: Lat ! Geodetic latitude in radians
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared

    lat = geodLat * deg2rad
    geocLat = atan2 ( sin(lat), f2*cos(lat) )

  end function GeodToGeocLat_D

  ! --------------------------------------------  GeodToGeocLat_S  -----
  real elemental function GeodToGeocLat_S ( GeodLat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ (IN DEGREES!) into a geocentric one
  !  $\lambda$ (IN RADIANS!)
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    use Constants, only: Deg2Rad

    ! Arguments
    real, intent(in) :: GeodLat

    real :: Lat ! Geodetic latitude in radians
    real, parameter :: F2 = Earth_Axis_Ratio_Squared

    lat = geodLat * deg2rad
    geocLat = atan2 ( sin(lat), f2*cos(lat) )

  end function GeodToGeocLat_S

  ! ----------------------------------------------------  Get_R_eq -----
  real(rp) elemental function Get_R_Eq ( Phi, Csq ) result ( R_eq )

  !{ Given the orbit geodetic longitude {\tt Phi} = $\phi$ in radians and the
  !  square of the minor axis of the orbit plane projected Earth ellipse in
  !  meters {\tt Csq} = $R_c^2$ compute the radius in kilometers of an
  !  equivalent circular Earth tangent to the elliptical Earth and having the
  !  same radius of curvature as the elliptical Earth at $\phi$.
  !%
  ! \begin{equation*}
  ! R_{eq} = \sqrt \frac{R_a^4 \sin^2 \phi + R_c^4 \cos^2 \phi}
  !                    {R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}
  !        = \sqrt \frac{R_a^4 - (R_a^2+R_c^2)(R_a^2-R_c^2) \cos^2 \phi}
  !                       {R_c^2 +              (R_a^2-R_c^2) \cos^2 \phi}
  ! \end{equation*}
  !%
  ! This is Equation (5.21) in the 19 August 2004 ATBD JPL D-18130.

  real(rp), intent(in) :: Phi
  real(rp), intent(in) :: Csq

  real(rp), parameter :: Earthrada_sq = earthrada ** 2
  real(rp), parameter :: Earthrada_4 = earthrada_sq ** 2

  r_eq = (earthrada_sq - csq) * COS(phi)**2
  ! Earthrad[abc] are in meters, but r_eq needs to be in km.
  r_eq = 0.001_rp * SQRT( &
    & ( earthrada_4 -(earthrada_sq + csq) * r_eq ) / &
    & ( csq + r_eq ) )

  end function Get_R_Eq

  !{ Converting geodetic Latitude (DEGREES!), Longitude (DEGREES!), and Height
  !  (km) above the mean Earth ellipsoid (mean sea level), to Cartesian
  !  coordinates in an Earth-Centered-Rotating frame:
  !  %
  !  \begin{equation*}\begin{split}
  !  X =\,& (N(\theta)+h) \cos\theta \cos\phi \\
  !  Y =\,& (N(\theta)+h) \cos\theta \sin\phi \\
  !  Z =\,& (N(\theta)(1-e^2)+h) \sin\theta \\
  !  \end{split}\end{equation*}
  !  %
  !  where $\theta$ is geodetic Latitude (DEGREES!), $\phi$ is Longitude
  !  (DEGREES!), $h$ is height above mean sea level (km), $e^2 = 1 -
  !  \frac{b^2}{f^2}$ is the eccentricity, $b$ is the semi-minor axis
  !  (polar radius), $f = 1 - \frac{b}a$ is the inverse of flattening, $a$
  !  is the semi-major exis (equatorial radius), and
  !  %
  !  \begin{equation*}
  !  N(\theta) = \frac{a}{\sqrt{1-e^2 \sin^2 \theta}} =
  !              \frac{a^2}{\sqrt{a^2 \cos^2 \theta + b^2 \sin^2 \theta}} =
  !              \frac{a^2}d \,.
  !  \end{equation*}
  !  is the distance along the normal from the surface to the polar axis.
  !
  !  The units of the results here ($X$, $Y$, and $Z$) are average Earth
  !  radii, not km; see the constant {\tt ERAD} above.

  ! --------------------------------------------------  To_Cart_D  -----
  subroutine To_Cart_D ( WHERE, CART, CT, ST, CP, SP )
  ! Convert Geodetic latitude and longitude (degrees), and
  ! altitude (km above sea level), to Cartesian
    use Constants, only: Deg2Rad
    double precision, intent(in) :: WHERE(3) ! Latitude (degrees north),
                                             ! Longitude (degrees east),
                                             ! Altitude (km above sea level)
    double precision, intent(out) :: CART(3) ! X, Y, Z, in average Earth radii
                                             ! (See the constant ERAD above)
    double precision, intent(out), optional :: CT, ST, CP, SP ! Cosine, Sine
                                             ! of Lat, Lon
    double precision :: D, MyCt, MySt, MyCp, MySp, RHO, RLAT, RLON

    rlat = where(1)*deg2Rad
    myst = sin(rlat)
    myct = cos(rlat)
    d = sqrt(aquad-(aquad-bquad)*myst*myst)
    rlon = where(2)*deg2Rad
    mycp = cos(rlon)
    mysp = sin(rlon)
    cart(3) = (where(3)+bquad/d)*myst/erad
    rho = (where(3)+aquad/d)*myct/erad
    cart(1) = rho*mycp
    cart(2) = rho*mysp
    if ( present(ct) ) ct = myst
    if ( present(st) ) st = myct
    if ( present(cp) ) cp = mycp
    if ( present(sp) ) sp = mysp

  end subroutine To_Cart_D

  ! --------------------------------------------------  To_Cart_S  -----
  subroutine To_Cart_S ( WHERE, CART, CT, ST, CP, SP )
  ! Convert Geodetic latitude and longitude (degrees), and
  ! altitude (km above sea level), to Cartesian
    use Constants, only: Deg2Rad
    real, intent(in) :: WHERE(3)        ! Latitude (degrees north),
                                        ! Longitude (degrees east),
                                        ! Altitude (km above sea level)
    real, intent(out) :: CART(3)        ! X, Y, Z, in average Earth radii
                                        ! (See the constant ERAD above)
    real, intent(out), optional :: CT, ST, CP, SP ! Cosine, Sine
                                        ! of Lat, Lon
    real :: D, MyCt, MySt, MyCp, MySp, RHO, RLAT, RLON

    rlat = where(1)*deg2Rad
    myst = sin(rlat)
    myct = cos(rlat)
    d = sqrt(aquad-(aquad-bquad)*myst*myst)
    rlon = where(2)*deg2Rad
    mycp = cos(rlon)
    mysp = sin(rlon)
    cart(3) = (where(3)+bquad/d)*myst/erad
    rho = (where(3)+aquad/d)*myct/erad
    cart(1) = rho*mycp
    cart(2) = rho*mysp
    if ( present(ct) ) ct = myst
    if ( present(st) ) st = myct
    if ( present(cp) ) cp = mycp
    if ( present(sp) ) sp = mysp

  end subroutine To_Cart_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Geometry

! $Log$
! Revision 2.21  2014/10/29 21:03:11  vsnyder
! Replaced confusing "myst=cos(rlat); myct=sin(rlat)" in To_Cart.  Described
! units of arguments.  Added a LaTeX comment to describe the mathematics.
!
! Revision 2.20  2013/08/16 02:27:56  vsnyder
! Add GeocToGeod, move To_Cart here from igrf_int
!
! Revision 2.19  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.18  2009/05/13 20:13:58  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.17  2008/10/08 01:11:37  vsnyder
! Add PRINT statement in not_used_here to prevent compilation cascades
!
! Revision 2.16  2008/10/08 01:07:59  vsnyder
! Add Get_R_Eq function
!
! Revision 2.15  2006/09/28 20:51:46  vsnyder
! Add Earth_Axis_Ratio_Squared_m1
!
! Revision 2.14  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.13  2003/01/16 23:13:10  livesey
! Added maxRefraction
!
! Revision 2.12  2003/01/15 02:45:50  vsnyder
! Make SecPerYear public (oops!)
!
! Revision 2.11  2003/01/15 02:35:08  vsnyder
! Add SecPerYear, move a USE to procedure scope
!
! Revision 2.10  2003/01/10 21:55:12  vsnyder
! Move SpeedOfLight from Geometry ot Units
!
! Revision 2.9  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/10/02 21:05:39  vsnyder
! Add SpeedOfLight constant
!
! Revision 2.7  2002/09/26 20:56:53  vsnyder
! Move Earth_Axis_Ratio_Squared from heights_module to Geometry
!
! Revision 2.6  2002/09/26 20:33:19  vsnyder
! Remove PI and LN10 -- they're in Units
!
! Revision 2.5  2002/09/26 16:27:14  livesey
! Changes from Van, new constants etc.
!
! Revision 2.4  2001/03/28 19:50:19  vsnyder
! Change constants from d0 to _r8
!
! Revision 2.3  2001/03/27 18:09:11  vsnyder
! Revised CVS stuff to use CHARACTER(len=*) parameter
!
! Revision 2.2  2001/03/26 13:23:34  livesey
! Added CVS stuff
!
