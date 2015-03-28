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

  use Constants, only: Pi
  use MLSKinds, only: R8, RP

  implicit NONE
  private

  ! Constants
  public :: Earth_Axis_Ratio, Earth_Axis_Ratio_Squared ! a^2/b^2
  public :: EarthRadA, EarthRadB, EarthSurfaceGPH, Eccentricity_Sq, ERad
  public :: GM, G0, J2, J4, SecPerYear, W, MaxRefraction

  ! Procedures
  public :: GeocToGeodLat, GeodToGeocLat, Get_R_Eq
  public :: Orbit_Plane_Minor_Axis_sq, To_Cart, To_XYZ
  public :: XYZ_to_Geod, XYZ_to_Geod_Bowring, XYZ_to_Geod_Fukushima

  ! Earth dimensions.

  real(r8), parameter :: EarthRadA = 6378137.0_r8    ! Semi-Major axis in m
  real(r8), parameter :: EarthRadB = 6356752.3141_r8 ! Semi-Minor axis in m
  real(r8), parameter :: Earth_Axis_Ratio = EarthRadB / EarthRadA
  real(r8), parameter :: Earth_Axis_Ratio_Squared = Earth_Axis_Ratio**2
  real(r8), parameter :: Eccentricity_Sq = &
    & 1.0_r8 - Earth_Axis_Ratio_Squared ! 1 - b**2 / a**2

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

  interface Orbit_Plane_Minor_Axis_sq
    module procedure Orbit_Plane_Minor_Axis_sq_D, Orbit_Plane_Minor_Axis_sq_S
  end interface

  interface To_Cart
    module procedure To_Cart_D, To_Cart_S
  end interface

  interface To_XYZ
    module procedure To_XYZ_D, To_XYZ_S
  end interface

  interface XYZ_to_Geod_Fukushima
    module procedure XYZ_to_Geod_Fukushima_D, XYZ_to_Geod_Fukushima_S
  end interface

  interface XYZ_to_Geod_Bowring
    module procedure XYZ_to_Geod_Bowring_D, XYZ_to_Geod_Bowring_S
  end interface

  interface XYZ_to_Geod
    module procedure XYZ_to_Geod_Fukushima_D, XYZ_to_Geod_Fukushima_S
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
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( sin(lat), f2 * cos(lat) )

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
    real, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( sin(lat), f2 * cos(lat) )

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
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = geodLat * deg2rad
    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

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
    real, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = geodLat * deg2rad
    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

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

  !{ Compute the square of the minor axis of orbit plane projected Earth
  !  ellipse $c$, where
  !  $c^2 = \frac{a^2\,b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta} =
  !         \frac{a^2}{\left(\frac{a^2}{b^2}-1\right) \sin^2 \beta + 1} =
  !         \frac{b^2}{1 - e^2 \cos^2 \beta}$ where $e^2$ is the square of
  !         the eccentricity, given by $1 - \frac{b^2}{a^2}$.
  !  This is Equation (5.3) in the 19 August 2004 ATBD JPL D-18130.

! ----------------------------------  Orbit_Plane_Minor_Axis_sq_D  -----

  pure function Orbit_Plane_Minor_Axis_sq_D ( Beta ) result ( Csq )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Beta  ! Orbit inclination, radians
    real(rk) :: Csq               ! Square of the minor axis of the
                                  ! orbit-plane projected ellipse.
    csq = EarthRadB**2 / ( 1.0_rk - Eccentricity_sq * cos(beta)**2 )
  end function Orbit_Plane_Minor_Axis_sq_D

  pure function Orbit_Plane_Minor_Axis_sq_S ( Beta ) result ( Csq )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Beta  ! Orbit inclination, radians
    real(rk) :: Csq               ! Square of the minor axis of the
                                  ! orbit-plane projected ellipse.
    csq = EarthRadB**2 / ( 1.0_rk - Eccentricity_sq * cos(beta)**2 )
  end function Orbit_Plane_Minor_Axis_sq_S

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
  !  \frac{b^2}{a^2} = 1 - ( 1 - f )^2$ is the eccentricity, $b$ is the
  !  semi-minor axis (polar radius), $f = 1 - \frac{b}a$ is the inverse of
  !  flattening, $a$ is the semi-major exis (equatorial radius), and
  !  %
  !  \begin{equation*}
  !  N(\theta) = \frac{a}{\sqrt{1-e^2 \sin^2 \theta}} =
  !              \frac{a^2}{\sqrt{a^2 \cos^2 \theta + b^2 \sin^2 \theta}} =
  !              \frac{a^2}d \,.
  !  \end{equation*}
  !  is the distance along the normal from the surface to the polar axis.
  !
  !  The units of the results here ($X$, $Y$, and $Z$) are average Earth
  !  radii (see the constant {\tt ERAD} above), or km.

  ! --------------------------------------------------  To_Cart_D  -----
  subroutine To_Cart_D ( WHERE, CART, CT, ST, CP, SP, KM )
  ! Convert Geodetic latitude and longitude (degrees), and altitude (km
  ! above sea level), to Cartesian coordinates.  The units are ERAD (see
  ! above) unless KM is present and true, in which case the units are
  ! kilometers.
    use Constants, only: Deg2Rad
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: WHERE(3) ! Latitude (degrees north),
                                     ! Longitude (degrees east),
                                     ! Altitude (km above sea level)
    real(rk), intent(out) :: CART(3) ! X, Y, Z, in average Earth radii
                                     ! (See the constant ERAD above)
    real(rk), intent(out), optional :: CT, ST, CP, SP ! Cosine, Sine
                                     ! of Lat, Lon
    logical, intent(in), optional :: KM ! Output units are kilometers

    include "To_Cart.f9h"

  end subroutine To_Cart_D

  ! --------------------------------------------------  To_Cart_S  -----
  subroutine To_Cart_S ( WHERE, CART, CT, ST, CP, SP, KM )
  ! Convert Geodetic latitude and longitude (degrees), and altitude (km
  ! above sea level), to Cartesian coordinates.  The units are ERAD (see
  ! above) unless KM is present and true, in which case the units are
  ! kilometers.
    use Constants, only: Deg2Rad
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: WHERE(3) ! Latitude (degrees north),
                                     ! Longitude (degrees east),
                                     ! Altitude (km above sea level)
    real(rk), intent(out) :: CART(3) ! X, Y, Z, in average Earth radii
                                     ! (See the constant ERAD above)
    real(rk), intent(out), optional :: CT, ST, CP, SP ! Cosine, Sine
                                     ! of Lat, Lon
    logical, intent(in), optional :: KM ! Output units are kilometers

    include "To_Cart.f9h"

  end subroutine To_Cart_S

  ! ---------------------------------------------------  To_XYZ_D  -----
  pure function To_XYZ_D ( Lat, Lon, Radians ) result ( XYZ )
    ! Convert north geocentric Lat (default degrees) and east Lon
    ! (default degrees) to a unit vector in ECR.
    use Constants, only: Deg2Rad
    double precision, intent(in) :: Lat, Lon
    logical, intent(in), optional :: Radians
    double precision :: XYZ(3)
    double precision :: MyLat, MyLon
    logical :: MyRad
    myRad = .false.
    if ( present(radians) ) myRad = radians
    if ( myRad ) then
      myLat = lat
      myLon = lon
    else
      myLat = lat * deg2rad
      myLon = lon * deg2rad
    end if
    xyz(1) = cos(myLat) * cos(myLon)
    xyz(2) = cos(myLat) * sin(myLon)
    xyz(3) = sin(myLat)
  end function To_XYZ_D

  ! ---------------------------------------------------  To_XYZ_S  -----
  pure function To_XYZ_S ( Lat, Lon, Radians ) result ( XYZ )
    ! Convert north geocentric Lat (default degrees) and east Lon
    ! (default degrees) to a unit vector in ECR.
    use Constants, only: Deg2Rad
    real, intent(in) :: Lat, Lon
    logical, intent(in), optional :: Radians
    real :: XYZ(3)
    real :: MyLat, MyLon
    logical :: MyRad
    myRad = .false.
    if ( present(radians) ) myRad = radians
    if ( myRad ) then
      myLat = lat
      myLon = lon
    else
      myLat = lat * deg2rad
      myLon = lon * deg2rad
    end if
    xyz(1) = cos(myLat) * cos(myLon)
    xyz(2) = cos(myLat) * sin(myLon)
    xyz(3) = sin(myLat)
  end function To_XYZ_S

! ------------------------------------------  XYZ_to_Geod_Bowring  -----

!{ Convert Geocentric Cartesian coordinates to Longitude, Geodetic Latitude,
!  and Geodetic Height using Bowring's method.
!  Start by estimating the geodetic latitude
!
! \begin{equation*}
! \beta = \tan^{-1} \left(\frac{az}{bs} \right)
! \end{equation*}
!
! where $s = \sqrt{x^2+y^2}$.  Then iterate the following until $\beta$
! converges:
!
! \begin{equation*}\begin{split}
! t = \,& \frac{b z + ( a^2-b^2 ) \sin^3 \beta}
!              {a s - ( a^2-b^2 ) \cos^3 \beta} \\
! \,& \\
! \beta = \,& \tan^{-1} \left( t \right) \\
! \end{split}\end{equation*}
!
! The geodetic height is then
!
! \begin{equation*}\begin{split}
! h = \,& s \cos \mu + ( z + e^2 N \sin \mu ) \sin \mu - N \\
!   = \,& s \cos \mu + z \sin \mu - N ( 1 - e^2 \sin^2 \mu ) \\
!   = \,& s \cos \mu + z \sin \mu - d \\
! \end{split}\end{equation*}
!
! where $\mu = \tan^{-1} \left( \frac{a}b t \right)$,
! $e^2 = 1 - \frac{b^2}{a^2}$ is the square of the first eccentricity,
!
! \begin{equation*}
! N = \frac{a}{\sqrt{1-e^2 \sin^2 \mu}} = \frac{a^2}d
! \end{equation*}
!
! is the radius of curvature in the meridian, and
! $d = a \sqrt{1 - e^2 \sin^2 \mu} = \sqrt{a^2 cos^2 \mu + b^2 \sin^2 \mu}$.
!
! \vspace*{5pt}
! The longitude is $\phi = \tan^{-1}\left( \frac{y}x \right)$.

  pure function XYZ_to_Geod_Bowring_D ( XYZ ) result ( Geod )
  ! Convert Geocentric Cartesian coordinates (XYZ) in meters to Longitude
  ! (Geod(2)) in radians, Geodetic Latitude (Geod(1)) in radians, and
  ! Geodetic Height (Geod(3)) in meters.
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: XYZ(3)
    real(rk) :: Geod(3) ! Geodetic latitude, longitude, geodetic height
    include 'XYZ_to_Geod_Bowring.f9h'
  end function XYZ_to_Geod_Bowring_D

  pure function XYZ_to_Geod_Bowring_S ( XYZ ) result ( Geod )
  ! Convert Geocentric Cartesian coordinates (XYZ) in meters to Longitude
  ! (Geod(2)) in radians, Geodetic Latitude (Geod(1)) in radians, and
  ! Geodetic Height (Geod(3)) in meters.
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: XYZ(3)
    real(rk) :: Geod(3) ! Geodetic latitude, longitude, geodetic height
    include 'XYZ_to_Geod_Bowring.f9h'
  end function XYZ_to_Geod_Bowring_S

! ----------------------------------------  XYZ_to_Geod_Fukushima  -----

!{In the Journal of Geodesy (2006) 79: 689-693, Toshio Fukushima presented 
! a method based on Halley's third-order iteration for the tangent of the
! reduced latitude.  Rather than iterating on the tangent, which requires
! divisions in each iteration, Fukushima iterates on denormalized quantities
! related to the sine and cosine, thereby not requiring any divisions or
! evaluations of trigonometric functions during the iteration.
! 
! \begin{equation*}\begin{split}
! s = \,& \sqrt{x^2 + y^2} \\
! e_c = \,& \frac{b}a \\
! E = \,& e^2 = 1 - \frac{b^2}{a^2} = 1 - e_c^2\\
! P = \,& \frac{s}a \\
! S_0 = \,& z \\
! C_0 = \,& e_c P \\
! A_n = \,& \sqrt{S_n^2+C_n^2} \\
! B_n = \,&  1.5 E S_n C_n^2 [ ( P S_n - z C_n ) A_n - E S_n C_n ] \\
! B_0 = \,& 1.5 E^2 P S_0^2 C_0^2 ( A_0 - e_c ) \\
! D_n = \,& z A_n^3 + E S_n^3 \\
! F_n = \,& P A_n^3 - E C_n^3 \\
! S_{n+1} = \,& D_n F_n - B_n S_n \\
! C_{n+1} = \,& F_n^2 - B_n C_n \\
! C_c = \,& e_c C_{n+1} \\
! \lambda = \,& \text{signum}(z) \, \tan^{-1} \frac{S_{n+1}}{C_c} \\
! h = \,& \frac{s C_c + |z| S_{n+1} - b A_{n+1}}
!              {\sqrt{C_c^2 + S_{n+1}^2}} \\
! \end{split}\end{equation*}
! 
! Fukushima's tests show that the method is accurate to within a few micro
! arcseconds for heights $<$ 30,000 km, and faster than other methods, with
! only one iteration.

  pure function XYZ_to_Geod_Fukushima_D ( XYZ ) result ( Geod )
  ! Convert Geocentric Cartesian coordinates (XYZ) in meters to Longitude
  ! (Geod(2)) in radians, Geodetic Latitude (Geod(1)) in radians, and
  ! Geodetic Height (Geod(3)) in meters.
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: XYZ(3)
    real(rk) :: Geod(3) ! Geodetic latitude, longitude, geodetic height
    include 'XYZ_to_Geod_Fukushima.f9h'
  end function XYZ_to_Geod_Fukushima_D

  pure function XYZ_to_Geod_Fukushima_S ( XYZ ) result ( Geod )
  ! Convert Geocentric Cartesian coordinates (XYZ) in meters to Longitude
  ! (Geod(2)) in radians, Geodetic Latitude (Geod(1)) in radians, and
  ! Geodetic Height (Geod(3)) in meters.
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: XYZ(3)
    real(rk) :: Geod(3) ! Geodetic latitude, longitude, geodetic height
    include 'XYZ_to_Geod_Fukushima.f9h'
  end function XYZ_to_Geod_Fukushima_S

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
! Revision 2.23  2015/03/28 02:58:05  vsnyder
! Added Earth_Axis_Ratio
!
! Revision 2.22  2015/03/28 01:55:20  vsnyder
! Deleted Earth_Axis_Ratio_Squared_m1.  Added Eccentricity_Sq, ERad,
! Orbit_Plane_Minor_Axis_sq, To_XYZ, XYZ_to_Geod, XYZ_to_Geod_Bowring,
! XYZ_to_Geod_Fukushima, Earth_Axis_Ratio.  Corrected errors in
! GeodToGeocLat and GeocToGeodLat.
!
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
