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

  use Constants, only: Deg2Rad, Pi, Rad2Deg
  use Earth_Constants, only: Earth_Axis_Ratio, Earth_Axis_Ratio_Squared, &
    & EarthRadA, EarthRadB, EarthSurfaceGPH, Eccentricity_Sq, &
    & GM, G0, J2, J4, SecPerYear, W
  use MLSKinds, only: RP

  implicit NONE
  private

  ! Constants gotten from Earth_Constants
  public :: Earth_Axis_Ratio, Earth_Axis_Ratio_Squared ! b^2/b^a
  public :: EarthRadA, EarthRadB, EarthSurfaceGPH, Eccentricity_Sq, ERad
  public :: GM, G0, J2, J4, SecPerYear, W

  ! Procedures
  public :: GeocToECRu, GeocToGeodLat, GeodToECRm, GeodToGeocAlt
  public :: GeodToGeocLat, GeodToGeocLatRad, Get_R_Eq, Great_Circle_Points
  public :: Phi_To_Lat_Deg, Phi_To_Lat_Rad
  public :: Orbit_Plane_Minor_Axis_sq, To_Cart, To_XYZ, XYZ_to_Geod
  public :: XYZ_to_Geod_Bowring, XYZ_to_Geod_Fukushima
  public :: XZ_to_Geod, XZ_to_Geod_Fukushima

  interface GeocToECRu ! Convert longitude and geocentric latitude
                       ! (both in degrees ) to an unit vector in ECR.
    module procedure To_XYZ_D, To_XYZ_S
  end interface

  interface GeocToGeodLat
    module procedure GeocToGeodLat_D, GeocToGeodLat_S
  end interface

  ! Convert geodetic coordinates to ECR in meters
  interface GeodToECRm
    module procedure GeodToECRm_D, GeodToECRm_S
  end interface GeodToECRm

  interface GeodToGeocAlt
    module procedure GeodToGeocAlt_D, GeodToGeocAlt_S
  end interface

  interface GeodToGeocLat ! Degrees -> Radians
    module procedure GeodToGeocLat_D, GeodToGeocLat_S
  end interface

  interface GeodToGeocLatRad ! Radians -> Radians
    module procedure GeodToGeocLatRad_D, GeodToGeocLatRad_S
  end interface

  interface Great_Circle_Points
    module procedure Great_Circle_Points_D, Great_Circle_Points_S
  end interface

  interface Phi_To_Lat_Deg
    module procedure Phi_To_Lat_Deg_D, Phi_To_Lat_Deg_S
  end interface

  interface Phi_To_Lat_Rad
    module procedure Phi_To_Lat_Rad_D, Phi_To_Lat_Rad_S
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

  interface XZ_to_Geod_Fukushima ! ( X_in, Z_in, A, B, GeodLat, GeodHt )
    ! This can be used to convert ECR to geodetic latitude and geodetic height
    ! with respect to the orbit-plane projected ellipse by using that ellipse's
    ! minor axis for B and sqrt(ECR%xyz(1)**2 + ECR%xyz(2)**2) for Z_in.
    module procedure XZ_to_Geod_Fukushima_D, XZ_to_Geod_Fukushima_S
  end interface

  interface XZ_to_Geod ! ( X_in, Z_in, A, B, GeodLat, GeodHt )
    ! This can be used to convert ECR to geodetic latitude and geodetic height
    ! with respect to the orbit-plane projected ellipse by using that ellipse's
    ! minor axis for B and sqrt(ECR%xyz(1)**2 + ECR%xyz(2)**2) for Z_in.
    module procedure XZ_to_Geod_Fukushima_D, XZ_to_Geod_Fukushima_S
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

  ! ----------------------------------------------  GeocToGeodLat  -----
  double precision elemental function GeocToGeodLat_D ( geocLat ) result ( GeodLat )

  !{ Convert a geocentric latitude $\lambda$ to a geodetic one $\mu$, both
  !  in degrees.
  !  Use the relation $\mu = \tan^{-1} \frac{f^2 \sin\lambda}{\cos\lambda}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    double precision, intent(in) :: geocLat

    double precision :: Lat ! Geocentric latitude in radians
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( sin(lat), f2 * cos(lat) )

  end function GeocToGeodLat_D

  real elemental function GeocToGeodLat_S ( geocLat ) result ( GeodLat )

  !{ Convert a geocentric latitude $\lambda$ to a geodetic one $\mu$, both
  !  in degrees.
  !  Use the relation $\mu = \tan^{-1} \frac{f^2 \sin\lambda}{\cos\lambda}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    real, intent(in) :: geocLat

    real :: Lat ! Geocentric latitude in radians
    real, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = deg2Rad * geocLat
    geodLat = rad2Deg * atan2 ( sin(lat), f2 * cos(lat) )

  end function GeocToGeodLat_S

  ! -----------------------------------------------  GeodToECRm_D  -----

  !{ Convert geodetic coordinates Latitude (DEGREES!), Longitude (DEGREES!)
  !  and meters above the mean Earth ellipsoid (mean sea level), to ECR
  !  coordinates in meters:
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

  pure function GeodToECRm_D ( Geod ) result ( ECR )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Geod(3) ! Latitude (degrees north),
                                    ! Longitude (degrees east),
                                    ! Altitude (meters above sea level)
    real(rk) :: ECR(3)              ! X, Y, Z in meters from the Earth center
    real(rk), parameter :: AQUAD = EarthRadA**2 ! m**2
    real(rk), parameter :: BQUAD = EarthRadB**2 ! m**2
    real(rk) :: D
    real(rk) :: MyCT ! Cos(theta) (lat)
    real(rk) :: MyLat, MyLon
    real(rk) :: MyST ! Sin(theta) (lat)
    real(rk) :: Rho
    myLat = geod(1) * deg2rad
    myLon = geod(2) * deg2rad
    myst = sin(myLat)
    myct = cos(myLat)
    d = sqrt(aquad*myct*myct + bquad*myst*myst)
!   d = sqrt(aquad-(aquad-bquad)*myst*myst)
    ECR(3) = ( geod(3) + bquad/d ) * mySt
    rho = ( geod(3) + aquad/d ) * myCt
    ECR(1) = rho * cos(myLon)
    ECR(2) = rho * sin(myLon)
  end function GeodToECRm_D

  pure function GeodToECRm_S ( Geod ) result ( ECR )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Geod(3) ! Latitude (degrees north),
                                    ! Longitude (degrees east),
                                    ! Altitude (meters above sea level)
    real(rk) :: ECR(3)              ! X, Y, Z in meters from the Earth center
    real(rk), parameter :: AQUAD = EarthRadA**2 ! m**2
    real(rk), parameter :: BQUAD = EarthRadB**2 ! m**2
    real(rk) :: D
    real(rk) :: MyCT ! Cos(theta) (lat)
    real(rk) :: MyLat, MyLon
    real(rk) :: MyST ! Sin(theta) (lat)
    real(rk) :: Rho
    myLat = geod(1) * deg2rad
    myLon = geod(2) * deg2rad
    myst = sin(myLat)
    myct = cos(myLat)
    d = sqrt(aquad*myct*myct + bquad*myst*myst)
!   d = sqrt(aquad-(aquad-bquad)*myst*myst)
    ECR(3) = ( geod(3) + bquad/d ) * mySt
    rho = ( geod(3) + aquad/d ) * myCt
    ECR(1) = rho * cos(myLon)
    ECR(2) = rho * sin(myLon)
  end function GeodToECRm_S

  ! ----------------------------------------------  GeodToGeocAlt  -----
  double precision function GeodToGeocAlt_D ( Where ) result ( GeocAlt )

  ! Convert geodetic altitude in meters above sea level to geocentric height in meters

    ! Arguments
    double precision, intent(in) :: Where(3) ! Geod. Lat(degrees), Lon(Degrees),
                                             ! Geod. height (Meters) above sea level

    ! Get length of ECR vector in meters
    geocAlt = norm2( geodToECRm ( [ where(1), where(2), where(3) ] ) )

  end function GeodToGeocAlt_D

  double precision function GeodToGeocAlt_S ( Where ) result ( GeocAlt )

  ! Convert geodetic altitude in meters above sea level to geocentric height in meters

    ! Arguments
    real, intent(in) :: Where(3)             ! Geod. Lat(degrees), Lon(Degrees),
                                             ! Geod. height (Meters) above sea level

    ! Get length of ECR vector in meters
    geocAlt = norm2( geodToECRm ( [ where(1), where(2), where(3) ] ) )

  end function GeodToGeocAlt_S

  ! ----------------------------------------------  GeodToGeocLat  -----
  double precision elemental function GeodToGeocLat_D ( GeodLat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ (IN DEGREES!) into a geocentric one
  !  $\lambda$ (IN RADIANS!)
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    ! Arguments
    double precision, intent(in) :: GeodLat

    double precision :: Lat ! Geodetic latitude in radians
    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = geodLat * deg2rad
    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

  end function GeodToGeocLat_D

  real elemental function GeodToGeocLat_S ( GeodLat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ (IN DEGREES!) into a geocentric one
  !  $\lambda$ (IN RADIANS!)
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    ! Arguments
    real, intent(in) :: GeodLat

    real :: Lat ! Geodetic latitude in radians
    real, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    lat = geodLat * deg2rad
    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

  end function GeodToGeocLat_S

  double precision elemental function GeodToGeocLatRad_D ( Lat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ in radians into a geocentric one
  !  $\lambda$ in radians
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    ! Arguments
    double precision, intent(in) :: Lat ! Geodetic latitude in radians

    double precision, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

  end function GeodToGeocLatRad_D

  real elemental function GeodToGeocLatRad_S ( Lat ) result ( GeocLat )

  !{ Convert a geodetic latitude $\mu$ in radians into a geocentric one
  !  $\lambda$ in radians
  !  Use the relation $\lambda = \tan^{-1} \frac{\sin\mu}{f^2 \cos\mu}$,
  !  where $f$ is the ratio of the equatorial to polar Earth radii.

    ! Arguments
    real, intent(in) :: Lat ! Geodetic latitude in radians

    real, parameter :: F2 = Earth_Axis_Ratio_Squared ! b^2/a^2

    geocLat = atan2 ( f2 * sin(lat), cos(lat) )

  end function GeodToGeocLatRad_S

  ! ------------------------------------------------------  Get_N  -----

  pure real(rp) elemental function Get_N ( Phi, Csq ) result ( N )

  !{ Given the orbit geodetic angle {\tt Phi} = $\phi$ in radians and
  !  the square of the minor axis of the orbit plane projected Earth
  !  ellipse in meters {\tt Csq} = $R_c^2$ compute the function $N(\phi)$
  !  defined immediately before Equation (5.11) in the 19 August 2004 ATBD
  !  JPL D-18130, in meters.  This would be the radius of curvature in the
  !  prime vertical of the Earth reference ellipsoid if Csq were $R_{b}^2$.
  !%
  !  \begin{equation*}
  !  N(\phi) = \frac{R_a^2}
  !                   {\sqrt{R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}}
  !          = \frac{R_a^2}
  !                 {\sqrt{(R_a^2 - R_c^2) \cos^2 \phi + R_c^2}}
  !          = \frac{R_a}{\sqrt{1 - e^2 \sin^2 \phi}} \,,
  ! \end{equation*}
  ! where eccentricity $e = 1 - \frac{R_c^2}{R_a^2}$.

    real(rp), intent(in) :: Phi ! Radians
    real(rp), intent(in) :: Csq ! Meters

    real(rp), parameter :: Earthrada_sq = earthrada ** 2

    N = earthrada_sq / sqrt ( ( earthrada_sq - csq ) * cos(phi)**2 + csq )

  end function Get_N

  ! ---------------------------------------------------  Get_R_eq  -----
  pure real(rp) elemental function Get_R_Eq ( Phi, Csq ) result ( R_eq )

  !{ Given the orbit geodetic angle {\tt Phi} = $\phi$ in radians and the
  !  square of the minor axis of the orbit plane projected Earth ellipse in
  !  meters {\tt Csq} = $R_c^2$ compute the radius in kilometers of an
  !  equivalent circular Earth tangent to the elliptical Earth and having the
  !  same radius of curvature as the elliptical Earth at $\phi$.
  !%
  ! \begin{equation*}\begin{split}
  ! R_{eq} =\,& \sqrt \frac{R_a^4 \sin^2 \phi + R_c^4 \cos^2 \phi}
  !                    {R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}
  !        = \sqrt \frac{R_a^4 - (R_a^2+R_c^2)(R_a^2-R_c^2) \cos^2 \phi}
  !                       {R_c^2 +              (R_a^2-R_c^2) \cos^2 \phi} \\
  !        = \,& N(\phi) \sqrt{ \sin^2 \phi +
  !                               \frac{R_c^4}{R_a^4} \cos^2 \phi}
  !        = N(\phi) \sqrt{ 1 + e^2 ( 1-2 e^2 ) \cos^2 \phi}\,, \\
  ! \end{split}\end{equation*}
  !%
  ! where eccentricity $e^2 = 1 - \frac{R_c^2}{R_a^2}$,
  ! and $N(\phi)$ is given above in function {\tt Get_N}.
  ! This is Equation (5.21) in the 19 August 2004 ATBD JPL D-18130,
  ! and
  !%
  ! \begin{equation*}
  !  R_c^2 = \frac{R_a^2 R_b^2}{R_a^2 \sin^2 \beta + R_b^2 \cos^2 \beta}\,,
  ! \end{equation*}
  ! where $R_a$ is the Earth equatorial radius, and $R_b$ is the Earth polar
  ! radius, which is Equation (5.3) in the ATBD.
  !
  ! This expression for $R_{eq}$ is incorrect.  In its derivation, the
  ! d N(phi) / d phi term was neglected. The correct result is
  !%
  ! \begin{equation*}
  ! R_{eq} = \frac{R_a^2 R_c^2}
  !               {(R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi)^\frac32}
  !        = \frac{R_a^2 R_c^2}
  !               {((R_a^2-R_c^2) \cos^2 \phi + R_c^2)^\frac32}
  !        = R_a \, \frac{1-e^2}{(1-e^2 \sin^2 \phi)^\frac32} \,.
  ! \end{equation*}

    real(rp), intent(in) :: Phi ! Radians
    real(rp), intent(in) :: Csq ! Meters

    real(rp), parameter :: EarthRada_sq = earthRada ** 2
    real(rp), parameter :: EarthRada_4 = earthRada_sq ** 2
    real(rp) :: E2 ! Eccentricity**2 in the orbit plane
    real(rp) :: AxisRatio_sq ! c^2/a^2

    logical, parameter :: Correct = .false.

    if ( correct ) then
      axisRatio_sq = Earthrada_sq / csq
      e2 = 1.0 - axisRatio_sq
      ! Earthrad[abc] are in meters, but r_eq needs to be in km.
      r_eq = 0.001_rp * earthRada * axisRatio_sq / &
                        sqrt( 1.0_rp - e2 * sin(phi)**2 )**3
    else
      r_eq = (earthrada_sq - csq) * cos(phi)**2
      ! Earthrad[abc] are in meters, but r_eq needs to be in km.
      r_eq = 0.001_rp * sqrt( &
        & ( earthrada_4 -(earthrada_sq + csq) * r_eq ) / &
        & ( csq + r_eq ) )
    end if

  end function Get_R_Eq

  ! --------------------------------------------  Get_R_eq_Center  -----
  pure function Get_R_eq_Center ( Phi, Req, Csq ) result ( D )

  !{ Given the orbit geodetic angle {\tt Phi} = $\phi$ in radians, the
  !  equivalent circular earth radius {\tt Req} = $R_{\text{eq}}^\oplus
  !  = H_t^\oplus$ in KILOMETERS(!), and the square of the minor axis of
  !  the orbit plane projected Earth ellipse {\tt Csq} = $R_c^2$ in
  !  METERS$^2$ compute the vector $d$, the vector from the center of the
  !  earth to the center of the equivalent circular earth in the orbit
  !  plane, in KILOMETERS(!), defined in Equation (5.22) in the 19
  !  August 2004 ATBD JPL D-18130.
  !%
  ! \begin{equation*}\begin{split}
  !  x = \,& ( N(\phi) - H_t^\oplus ) \cos \phi \\
  !  y = \,& ( \frac{c^2}{a^2} N(\phi) - H_t^\oplus ) \sin \phi \\
  !  d = \,& [\, x,\, y,\, 0\, ] \text{ where}\\
  !  H_t^\oplus \equiv \,& R_\text{eq}^\oplus \\
  ! \end{split}\end{equation*}

    real(rp), intent(in) :: Phi ! Radians
    real(rp), intent(in) :: Req ! Kilometers
    real(rp), intent(in) :: Csq ! Meters
    real(rp) :: D(3)            ! Kilometers.

    real(rp), parameter :: Earthrada_sq = earthrada ** 2
    real(rp) :: N

    N = 0.001_rp * get_N ( phi, csq ) ! KILOMETERS!
    d = [ ( N - req ) * cos(phi), &
        & ( csq/earthrada_sq * N  - req) * sin(phi), &
        & 0.0_rp ]

  end function Get_R_eq_Center

  ! -------------------------------------------------  Get_H_eq_S  -----
  pure real(rp) function Get_H_eq_S ( SCECR, D )

  !{ Given the instrument position SCECR in ECR (kilometers) and the vector
  !  from the center of the Earth to the center of the equivalent circular
  !  Earth in the orbit plane projected Earth ellipse (kilometers), compute
  !  the instrument height measured from the center of the equivalent
  !  circular Earth, in kilometers. This is $|\vec{H}_s|$ where $\vec{H}_s =
  !  \vec{R}_s - \vec{d}$, $\vec{R}_s$ is given by SCECR, and $\vec{d}_s$ is
  !  given by D.  This appears before Equation (8.5) in the 19 August 2004
  !  ATBD JPL D-18130.

    real(rp), intent(in) :: SCECR(3) ! Kilometers
    real(rp), intent(in) :: D(3)     ! Kilometers

    get_H_eq_S = norm2 ( scecr - d )

  end function Get_H_eq_S

! ------------------------------------------  Great_Circle_Points  -----

  subroutine Great_Circle_Points_D ( Lon1, GeocLat1, Lon2, GeocLat2, Phi, Lon, GeocLat )
  ! Compute longitude and geocentric latitude of points on the great circle
  ! defined by R1 and R2, with Phi(1) at R1, and subsequent points along the
  ! great circle spaced at |Phi(n)-Phi(1)|.
    use Cross_m, only: Cross
    use Rotation_m, only: Rotate_3d
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Lon1, GeocLat1, Lon2, GeocLat2
    real(rk), intent(in) :: Phi(:)         ! Along-track angles, degrees
    real(rk), intent(out) :: Lon(:)        ! Longitude, degrees
    real(rk), intent(out) :: GeocLat(:)    ! Geocentric latitude, degrees
    include "Great_Circle_Points.f9h"
  end subroutine Great_Circle_Points_D

  subroutine Great_Circle_Points_S ( Lon1, GeocLat1, Lon2, GeocLat2, Phi, Lon, GeocLat )
  ! Compute longitude and geocentric latitude of points on the great circle
  ! defined by R1 and R2, with Phi(1) at R1, and subsequent points along the
  ! great circle spaced at |Phi(n)-Phi(1)|.
    use Cross_m, only: Cross
    use Rotation_m, only: Rotate_3d
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Lon1, GeocLat1, Lon2, GeocLat2
    real(rk), intent(in) :: Phi(:)         ! Along-track angles, degrees
    real(rk), intent(out) :: Lon(:)        ! Longitude, degrees
    real(rk), intent(out) :: GeocLat(:)    ! Geocentric latitude, degrees
    include "Great_Circle_Points.f9h"
  end subroutine Great_Circle_Points_S

! ------------------------------------  Orbit_Plane_Minor_Axis_sq  -----

  !{ Compute the square of the ratio of the minor axis to the major axis
  !  lengths of the orbit plane projected Earth ellipse $c$, where
  !  $c^2 = \frac{a^2\,b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta} =
  !         \frac{a^2}{\left(\frac{a^2}{b^2}-1\right) \sin^2 \beta + 1} =
  !         \frac{b^2}{1 - e^2 \cos^2 \beta}$ and $e^2$ is the square of
  !         the eccentricity, given by $1 - \frac{b^2}{a^2}$.
  !  This is Equation (5.3) in the 19 August 2004 ATBD JPL D-18130.

  pure elemental function Orbit_Plane_Minor_Axis_sq_D ( Beta ) result ( Csq )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Beta  ! Orbit inclination, radians
    real(rk) :: Csq               ! Square of the minor axis of the
                                  ! orbit-plane projected ellipse.
    csq = EarthRadB**2 / ( 1.0_rk - Eccentricity_sq * cos(beta)**2 )
  end function Orbit_Plane_Minor_Axis_sq_D

  pure elemental function Orbit_Plane_Minor_Axis_sq_S ( Beta ) result ( Csq )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Beta  ! Orbit inclination, radians
    real(rk) :: Csq               ! Square of the minor axis of the
                                  ! orbit-plane projected ellipse.
    csq = EarthRadB**2 / ( 1.0_rk - Eccentricity_sq * cos(beta)**2 )
  end function Orbit_Plane_Minor_Axis_sq_S

! ---------------------------------------------------  Phi_To_Lat  -----

  !{ Given the orbit angle $\phi$ and inclination $\beta$, compute the
  !  latitude $\lambda$ corresponding to the orbit angle.  If the orbit
  !  angle is geodetic (geocentric), the latitude is geodetic
  !  (geocentric).
  !  \begin{equation}
  !    \lambda = \sin^{-1} \left( \sin \phi \sin \beta \right)
  !  \end{equation}
  !  where $a$ is the Earth's semi-major (equatorial) axis length, and $c$
  !  is the ratio of the minor axis to the major axis lengths of the
  !  Earth-projected orbit-plane ellipse.

  pure elemental function Phi_To_Lat_Deg_D ( Phi, Beta ) result ( Lat )
    integer, parameter :: RK = kind(0.0d0)
    real(rk) :: Lat              ! Latitude, degrees (geocentric or geodetic)
    real(rk), intent(in) :: Phi  ! Orbit angle, degrees (geocentric or geodetic)
    real(rk), intent(in) :: Beta ! Orbit inclination, degrees
    lat = rad2deg * asin( sin(deg2rad * phi) * sin(deg2rad * beta) )
  end function Phi_To_Lat_Deg_D

  pure elemental function Phi_To_Lat_Deg_S ( Phi, Beta ) result ( Lat )
    integer, parameter :: RK = kind(0.0e0)
    real(rk) :: Lat              ! Latitude, degrees (geocentric or geodetic)
    real(rk), intent(in) :: Phi  ! Orbit angle, degrees (geocentric or geodetic)
    real(rk), intent(in) :: Beta ! Orbit inclination, degrees
    lat = rad2deg * asin( sin(deg2rad * phi) * sin(deg2rad * beta) )
  end function Phi_To_Lat_Deg_S

  pure elemental function Phi_To_Lat_Rad_D ( Phi, Beta ) result ( Lat )
    integer, parameter :: RK = kind(0.0d0)
    real(rk) :: Lat              ! Latitude, radians (geocentric or geodetic,
                                 ! depending on Phi)
    real(rk), intent(in) :: Phi  ! Orbit angle, radians (geocentric or geodetic)
    real(rk), intent(in) :: Beta ! Orbit inclination, radians
    lat = asin( sin(phi) * sin(beta) )
  end function Phi_To_Lat_Rad_D

  pure elemental function Phi_To_Lat_Rad_S ( Phi, Beta ) result ( Lat )
    integer, parameter :: RK = kind(0.0e0)
    real(rk) :: Lat              ! Latitude, radians (geocentric or geodetic,
                                 ! depending on Phi)
    real(rk), intent(in) :: Phi  ! Orbit angle, radians (geocentric or geodetic)
    real(rk), intent(in) :: Beta ! Orbit inclination, radians
    lat = asin( sin(phi) * sin(beta) )
  end function Phi_To_Lat_Rad_S

  ! ----------------------------------------------------  To_Cart  -----

  !{ Convert geodetic Latitude (DEGREES!), Longitude (DEGREES!), and Height
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

  subroutine To_Cart_D ( WHERE, CART, CT, ST, CP, SP, KM )
  ! Convert Geodetic latitude and longitude (degrees), and altitude (km
  ! above sea level), to Cartesian coordinates.  The units are ERAD (see
  ! above) unless KM is present and true, in which case the units are
  ! kilometers.

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

  subroutine To_Cart_S ( WHERE, CART, CT, ST, CP, SP, KM )
  ! Convert Geodetic latitude and longitude (degrees), and altitude (km
  ! above sea level), to Cartesian coordinates.  The units are ERAD (see
  ! above) unless KM is present and true, in which case the units are
  ! kilometers.

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

  ! -----------------------------------------------------  To_XYZ  -----
  pure function To_XYZ_D ( Lat, Lon, Radians ) result ( XYZ )
    ! Convert north geocentric Lat (default degrees) and east Lon
    ! (default degrees) to a unit vector in ECR.

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

  pure function To_XYZ_S ( Lat, Lon, Radians ) result ( XYZ )
    ! Convert north geocentric Lat (default degrees) and east Lon
    ! (default degrees) to a unit vector in ECR.

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

  elemental subroutine XZ_to_Geod_Fukushima_D ( X_in, Z_in, A, B, GeodLat, GeodHt )
  ! Convert Cartesian coordinates of a point on an ellipse to geodetic latitude
  ! in radians and geodetic height in the same units as X_in, Z_in, A, and B.
  ! This can be used to convert ECR to geodetic latitude and geodetic height
  ! with respect to the orbit-plane projected ellipse by using that ellipse's
  ! minor axis for B and sqrt(ECR%xyz(1)**2 + ECR%xyz(2)**2) for Z_in.
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: X_in ! X-coordinate of a point on the ellipse
    real(rk), intent(in) :: Z_in ! Z-coordinate of a point on the ellipse
    real(rk), intent(in) :: A    ! Ellipse semi-major axis
    real(rk), intent(in) :: B    ! Ellipse semi-minor axis
    real(rk), intent(out) :: GeodLat ! Geodetic latitude of (X_in, Z_in), radians
    real(rk), intent(out) :: GeodHt  ! Geodetic height of (X_in, Z_in), in same
                                     ! units as X_in, Z_in, A, B
    include 'XZ_to_Geod_Fukushima.f9h'
  end subroutine XZ_to_Geod_Fukushima_D

  elemental subroutine XZ_to_Geod_Fukushima_S ( X_in, Z_in, A, B, GeodLat, GeodHt )
  ! Convert Cartesian coordinates of a point on an ellipse to geodetic latitude
  ! in radians and geodetic height in the same units as X_in, Z_in, A, and B
  ! This can be used to convert ECR to geodetic latitude and geodetic height
  ! with respect to the orbit-plane projected ellipse by using that ellipse's
  ! minor axis for B and sqrt(ECR%xyz(1)**2 + ECR%xyz(2)**2) for Z_in.
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: X_in ! X-coordinate of a point on the ellipse
    real(rk), intent(in) :: Z_in ! Z-coordinate of a point on the ellipse
    real(rk), intent(in) :: A    ! Ellipse semi-major axis
    real(rk), intent(in) :: B    ! Ellipse semi-minor axis
    real(rk), intent(out) :: GeodLat ! Geodetic latitude of (X_in, Z_in), radians
    real(rk), intent(out) :: GeodHt  ! Geodetic height of (X_in, Z_in), in same
                                     ! units as X_in, Z_in, A, B
    include 'XZ_to_Geod_Fukushima.f9h'
  end subroutine XZ_to_Geod_Fukushima_S

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
! Revision 2.37  2020/01/23 01:58:54  vsnyder
! Correct a typo in an equation in TeXnicalities in Get_R_Eq
!
! Revision 2.36  2020/01/08 21:49:54  vsnyder
! Added a "correct" method to compute equivalent-circular-eqrth radius,
! under a named constant currently false.  Corrected and improved some
! TeXnicalities.
!
! Revision 2.35  2019/11/26 19:40:41  vsnyder
! Correct a typo in a comment
!
! Revision 2.34  2016/12/07 23:00:55  vsnyder
! Remove unused use name
!
! Revision 2.33  2016/09/02 00:23:25  vsnyder
! Add GeodToGeocLatRad, correct Phi_To_Lat
!
! Revision 2.32  2016/08/30 20:27:51  vsnyder
! Add Phi_To_Lat functions
!
! Revision 2.31  2016/06/02 02:11:32  vsnyder
! Make Orbit_Plane_Minor_Axis_sq elemental
!
! Revision 2.30  2016/05/27 01:10:52  vsnyder
! Publish XZ_to_Geod
!
! Revision 2.29  2016/05/25 01:52:00  vsnyder
! Add XZ_to_Geod and XZ_to_Geod_Fukushima
!
! Revision 2.28  2016/01/23 02:45:27  vsnyder
! Add GeodToECRm (meters); get constants from Earth_Constants
!
! Revision 2.27  2015/10/22 20:16:11  vsnyder
! Add Great_Circle_Points
!
! Revision 2.26  2015/09/23 22:37:01  vsnyder
! Correct a comment
!
! Revision 2.25  2015/09/22 23:05:06  vsnyder
! Correct some comments about geodetic vs geocentric latitudes
!
! Revision 2.24  2015/04/29 00:54:52  vsnyder
! Add GeodToGeocAlt
!
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
