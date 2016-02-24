! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Geolocation_0
!=============================================================================

  ! Define types used in geolocation modules.

  use MLSKinds, only: RG => R8 ! Kind of REAL geolocation variables
  use Earth_Constants, only: Earth_Axis_Ratio, EarthRadA, EarthRadB, &
    & Eccentricity_Sq

  implicit NONE

  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  type :: H_t           ! For horizontal (lat/lon) grids, geoc/geod unspecified
    real(rg) :: Lon     ! Degrees
    real(rg) :: Lat     ! Degrees
  contains
    procedure :: ECR => H_t_To_ECR_Surf ! at the Earth's surface
  end type H_t

  type, extends(h_t) :: H_Geoc
  !  real(rg) :: Lon    ! Longitude, degrees
  !  real(rg) :: Lat    ! Geocentric latitude, degrees
  contains
    procedure :: Geod => Geoc_To_Geod
  end type H_Geoc

  type, extends(h_t) :: H_Geod
  !  real(rg) :: Lon    ! Longitude, degrees
  !  real(rg) :: Lat    ! Geodetic latitude, degrees
  contains
    procedure :: ECR => Geod_To_ECR_Surf
    procedure :: Geoc => Geod_To_Geoc
  end type H_Geod

  type, abstract :: Lat_t ! For latitude; can't instantiate it
    real(rg) :: D       ! Degrees
  end type Lat_t

  type, extends(lat_t) :: GeocLat_t
  ! real(rg) :: D       ! Geocentric latitude, degrees
  contains
    procedure :: Geod => Geoc_To_Geod_Lat
  end type GeocLat_t

  type, extends(lat_t) :: GeodLat_t
  ! real(rg) :: D       ! Geodetic latitude, degrees
  contains
    procedure :: Geoc => Geod_To_Geoc_Lat
  end type GeodLat_t

  type :: Lon_t         ! Longitude, degrees
    real(rg) :: D       ! Degrees
  end type Lon_t

  type :: V_t           ! For vertical grids of unspecified type
    real(rg) :: V       ! Geocentric (meters), geodetic (meters), or zeta
  end type V_t

  type, extends(v_t) :: V_geoc
  ! real(rg) :: V       ! geocentric height (meters from the Earth center)
  end type V_geoc

  type, extends(v_t) :: V_geod
  ! real(rg) :: V       ! geodetic height (meters above sea level)
  end type V_geod

  type, extends(v_t) :: V_zeta
  ! real(rg) :: V       ! Zeta -- geocentric latitude doesn't make sense if
                        ! the coordinate system is stacked
  end type V_zeta

  ! For irregular grids, or lines of sight, where H and V are expected to
  ! have the same extents
  type, extends(h_t) :: H_V_t
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! Meters or zeta
  end type H_V_t

  type, extends(h_geoc) :: H_V_Geoc
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! geocentric height (meters from the Earth center)
  contains
    procedure :: ECR => Geoc_To_ECR
  ! procedure :: ECR_Unit => Geod_To_ECR_Unit ! Inherited from H_t via H_Geoc
    procedure :: GeodV => GeocV_To_GeodV
  end type H_V_Geoc

  type, extends(h_geod) :: H_V_Geod
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! geodetic height (meters above sea level)
  contains
    procedure :: ECR => Geod_To_ECR
  ! procedure :: ECR_Unit => Geod_To_ECR_Unit ! Inherited from H_t via H_Geod
    procedure :: GeocV => GeodV_To_GeocV
  end type H_V_Geod

  type, extends(h_geod) :: H_V_Zeta
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! Zeta
  end type H_V_Zeta

  type :: ECR_t
    real(rg) :: XYZ(3)  ! Cartesian ECR, meters, or a unit vector
  contains
    procedure :: Add => ECR_add
    generic :: operator(+) => Add
    procedure :: Cross => ECR_Cross
    generic :: operator(.CROSS.) => Cross
    procedure :: Dot => ECR_Dot
    generic :: operator(.DOT.) => Dot
    procedure :: Div => ECR_Divide
    generic :: operator(/) => Div
    procedure :: Geoc => ECR_to_Geoc
    procedure :: Geod => ECR_To_Geod_Fukushima
    procedure, pass(B) :: Scale_L => ECR_Scale_L
    procedure :: Scale_R => ECR_Scale_R
    procedure :: Subtract => ECR_Subtract
    generic :: operator(-) => Subtract
    generic :: operator(*) => Scale_L, Scale_R
  end type ECR_t

  interface Cross
    module procedure ECR_Cross
  end interface

  interface Dot_Product
    module procedure ECR_Dot
  end interface

  interface H_Geoc ! Construct H_Geoc object from an H_t one
    module procedure H_t_to_H_Geoc
  end interface

  interface H_Geod ! Construct H_Geod object from an H_t one
    module procedure H_t_to_H_Geod
  end interface

! *****     Private Constants     **************************************

! AQUAD   Square of major half axis for earth ellipsoid
! BQUAD   Square of minor half axis for earth ellipsoid

  real(rg), parameter :: AQUAD = EarthRadA**2 ! m**2
  real(rg), parameter :: BQUAD = EarthRadB**2 ! m**2

contains

  pure elemental type(ECR_t) function ECR_Add ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_Add = ECR_t ( a%xyz + b%xyz )
  end function ECR_Add

  pure type(ECR_t) function ECR_Cross ( A, B )
    use Cross_m, only: Cross
    class(ECR_t), intent(in) :: A, B
    ECR_Cross%xyz = cross(a%xyz,b%xyz)
  end function ECR_Cross

  pure type(ECR_t) function ECR_Divide ( A, B )
    class(ECR_t), intent(in) :: A
    real(rg), intent(in) :: B
    ECR_Divide = ECR_t ( a%xyz / b )
  end function ECR_Divide

  pure real(rg) function ECR_Dot ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_Dot = dot_product ( a%xyz, b%xyz )
  end function ECR_Dot

  pure type(ECR_t) function ECR_Scale_L ( A, B )
    real(rg), intent(in) :: A
    class(ECR_t), intent(in) :: B
    ECR_Scale_L = ECR_t ( a * b%xyz )
  end function ECR_Scale_L

  pure type(ECR_t) function ECR_Scale_R ( A, B )
    class(ECR_t), intent(in) :: A
    real(rg), intent(in) :: B
    ECR_Scale_R = ECR_t ( a%xyz * b )
  end function ECR_Scale_R

  pure elemental type(ECR_t) function ECR_Subtract ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_subtract = ECR_t ( a%xyz - b%xyz )
  end function ECR_Subtract

  pure elemental type(h_v_geoc) function ECR_To_Geoc ( ECR ) result ( Geoc )
    use Constants, only: rad2deg
    class(ECR_t), intent(in) :: ECR
    real(rg) :: Rho, Rho2
    rho2 = ECR%xyz(1)**2 + ECR%xyz(2)**2
    rho = sqrt(rho2)
    geoc%lon = rad2deg * atan2(ECR%xyz(2), ECR%xyz(1))
    geoc%lat = rad2deg * asin(ECR%xyz(3)/rho)
    geoc%v = sqrt(rho2 + ECR%xyz(3)**2)
  end function ECR_To_Geoc

  pure elemental type(h_v_geod) function ECR_To_Geod_Fukushima ( ECR ) result ( Geod )
    use Constants, only: rad2deg
    class(ECR_t), intent(in) :: ECR

    real(rg), parameter :: A = EarthRadA, B = EarthRadB
    real(rg), parameter :: Ec = Earth_Axis_Ratio ! b/a
    ! Square of first eccentricity = 1 - a^2/b^2
    real(rg), parameter :: E2 = Eccentricity_Sq
    real(rg), parameter :: Tol = tiny(1.0_rg)**(1.0/3.0) / A
    integer, parameter :: MaxIt = 1  ! Maximum number of iterations.  In Journal
                                     ! of Geodesy (2006) 79: 689-693, Fukushima
                                     ! says one iteration gets within a few micro
                                     ! arcseconds for heights below 30,000 km.

    real(rg) :: An, An3, Ap, Bn, Dn, Fn  ! Terms in the iteration
    real(rg) :: Cc                       ! Ec * Cn
    real(rg) :: Cn, Cp                   ! Cosine-like iterate, previous value
    integer :: I
    real(rg) :: P                        ! S / a
    real(rg) :: S                        ! sqrt ( x**2 + y**2 )
    real(rg) :: Sn, Sp                   ! Sine-like iterate, previous value
    real(rg) :: Z                        ! Ec * |z| / a

    ! Computing longitude is easy
    geod%lon = rad2deg * atan2 ( ECR%xyz(2), ECR%xyz(1) )

    s = sqrt ( ECR%xyz(1)**2 + ECR%xyz(2)**2 )
    if ( s < sqrt(tiny(1.0_rg)) ) then ! Too near the poles
      geod%lat = sign(90.0_rg,ECR%xyz(3))
      geod%v = ECR%xyz(3) - b
      return
    end if

    ! Now the Fukushima iteration for latitude
    P = s / a
    Z = abs(ECR%xyz(3)) * Ec / a
    Sp = Z
    Cp = Ec * P
    Ap = sqrt ( Sp**2 + Cp**2 )
    Bn = 1.5_rg * E2**2 * P * Sp**2 * Cp**2 * ( Ap - Ec )
    i = 1
    do
      An3 = Ap ** 3
      Dn = Z * An3 + E2 * Sp**3
      Fn = P * An3 - E2 * Cp**3
      Sn = Dn * Fn - Bn * Sp
      Cn = Fn**2 - Bn * Cp
      An = sqrt ( Sn**2 + Cn**2 )
      if ( i >= maxIt ) exit
      i = i + 1
      if ( abs(Cn-Cp) < tol * abs(p) .or. abs(Sn-Sp) < tol * Z ) exit
      Ap = An
      Bn = 1.5_rg * E2 * Sn * Cn**2 * ( ( P * Sn - Z * Cn ) * Ap - &
         & E2 * Sn * Cn )
      Cp = Cn
      Sp = Sn
    end do
    Cc = Ec * Cn
    ! Compute geodetic latitude
    geod%lat = sign ( rad2deg * atan2(Sn,Cc), ECR%xyz(3) )
    ! Compute geodetic height
    geod%v = ( s * Cc + abs(ECR%xyz(3)) * Sn - b * An ) / sqrt ( Cc**2 + Sn**2 )

  end function ECR_To_Geod_Fukushima

  pure elemental type(ECR_t) function Geoc_To_ECR ( Geo ) result ( ECR )
    ! Convert geocentric coordinates to ECR coordinates in meters
    use Constants, only: deg2rad
    class(h_v_geoc), intent(in) :: Geo
    real(rg) ::Lat, Lon, Rho
    lat = geo%lat * deg2rad
    lon = geo%lon * deg2rad
    rho = geo%v * cos(lat)
    ECR%xyz(1) = rho * cos(lon)
    ECR%xyz(2) = rho * sin(lon)
    ECR%xyz(3) = geo%v * sin(lat)
  end function Geoc_To_ECR

  pure elemental function Geoc_To_Geod ( Geoc ) result ( Geod )
    use Constants, only: deg2rad, rad2deg
    use Earth_Constants, only: F2 => Earth_Axis_Ratio_Squared ! (b/a)**2
    class(h_geoc), intent(in) :: Geoc
    type(h_geod) :: Geod
    real(rg) :: Lat
    geod%lon = geoc%lon
    lat = geoc%lat * deg2rad
    geod%lat = atan2 ( sin(lat), f2 * cos(lat) ) * rad2deg
  end function Geoc_To_Geod

  pure elemental function Geoc_To_Geod_Lat ( Geoc ) result ( Geod )
    use Constants, only: deg2rad, rad2deg
    use Earth_Constants, only: F2 => Earth_Axis_Ratio_Squared ! (b/a)**2
    class(geocLat_t), intent(in) :: Geoc
    type(geodLat_t) :: Geod
    real(rg) :: Lat
    lat = geoc%d * deg2rad
    geod%d = atan2 ( sin(lat), f2 * cos(lat) ) * rad2deg
  end function Geoc_To_Geod_Lat

  pure elemental type(h_v_geod) function GeocV_To_GeodV ( Geoc ) result ( Geod )
    class(h_v_geoc), intent(in) :: Geoc
    type(ECR_t) :: ECR
  ! geod = geoc%ECR()%geod()   ! Prohibited (why?)
    ECR = geoc%ECR()
    geod = ECR%geod()
  end function GeocV_To_GeodV

  pure elemental type(ECR_t) function Geod_To_ECR ( Geo ) result ( ECR )
    ! Convert geodetic coordinates to ECR coordinates in meters
    use Constants, only: deg2rad
    class(h_v_geod), intent(in) :: Geo
    real(rg) :: MyLon
    real(rg) :: Rho
    call geod_to_geoc_support ( geo, ECR%xyz(3), rho )
    myLon = geo%lon * deg2rad
    ECR%xyz(1) = rho * cos(myLon)
    ECR%xyz(2) = rho * sin(myLon)
  end function Geod_To_ECR

  pure elemental type(ECR_t) function Geod_To_ECR_Surf ( Geo ) result ( ECR )
    ! Convert geodetic coordinates at the Earth surface to ECR coordinates
    use Constants, only: deg2rad
    class(h_geod), intent(in) :: Geo
    real(rg) :: CT   ! cos(lat)
    real(rg) :: D    ! sqrt(a^2 ct^2 + b^2 st^2)
    real(rg) :: Lat  ! geodetic latitude in radians
    real(rg) :: Lon  ! longitude in radians
    real(rg) :: Rho  ! a^2/d cos(lat)
    real(rg) :: ST   ! sin(lat)
    lat = geo%lat * deg2rad
    ct = cos(lat)
    st = sin(lat)
    d = 1.0 / sqrt(aquad*ct*ct + bquad*st*st)
    rho = aquad * d * ct
    lon = geo%lon * deg2rad
    ecr%xyz(1) = rho * cos(lon)
    ecr%xyz(2) = rho * sin(lon)
    ecr%xyz(3) = bquad * d * st
  end function Geod_To_ECR_Surf

  pure elemental type(ECR_t) function Geod_To_ECR_Unit ( Geod ) result ( ECR )
    ! Convert geodetic coordinates to unit ECR coordinates
    class(h_v_geod), intent(in) :: Geod
    ECR = geod%ECR()
    ECR%xyz = ECR%xyz / norm2(ECR%xyz)
  end function Geod_To_ECR_Unit

  pure elemental type(h_geoc) function Geod_To_Geoc ( Geod ) result ( Geoc )
    use Constants, only: deg2rad, rad2deg
    use Earth_Constants, only: F2 => Earth_Axis_Ratio_Squared ! (b/a)**2
    class(h_geod), intent(in) :: Geod
    real(rg) :: Lat
    geoc%lon = geod%lon
    lat = geod%lat * deg2rad
    geoc%lat = atan2 ( f2 * sin(lat), cos(lat) ) * rad2deg
  end function Geod_To_Geoc

  pure elemental type(GeocLat_t) function Geod_To_Geoc_Lat ( Geod ) result ( Geoc )
    use Constants, only: deg2rad, rad2deg
    use Earth_Constants, only: F2 => Earth_Axis_Ratio_Squared ! (b/a)**2
    class(GeodLat_t), intent(in) :: Geod
    real(rg) :: Lat
    lat = geod%d * deg2rad
    geoc%d = atan2 ( f2 * sin(lat), cos(lat) ) * rad2deg
  end function Geod_To_Geoc_Lat

  pure elemental type(h_v_geoc) function GeodV_To_GeocV ( Geod ) result ( Geoc )
    ! Convert geodetic coordinates to Geocentric coordinates, including heights
    use Constants, only: rad2deg
    class(h_v_geod), intent(in) :: Geod
    real(rg) :: Rho
    real(rg) :: Z    ! Would be Z coordinate of ECR
    call geod_to_geoc_support ( geod, z, rho )
    geoc%lon = geod%lon
    geoc%lat = atan2 ( z, rho ) * rad2deg
    geoc%v = sqrt ( rho**2 + z**2 ) ! geocentric height
  end function GeodV_To_GeocV

  pure elemental type(ECR_t) function H_t_To_ECR_Surf ( Geo ) result ( ECR )
    ! Convert geocentric coordinates at the Earth surface to ECR coordinates
    use Constants, only: deg2rad
    use Earth_Constants, only: E2 => Eccentricity_Sq
    class(h_t), intent(in) :: Geo
    real(rg) :: Lat, Lon, R, Rho, ST
    lat = geo%lat * deg2rad
    st = sin(lat)
    ! Geocentric radius at geocentric latitude
    ! r = ab / sqrt ( a^2 cos^2(lat) + b^2 sin^2(lat) )
    r = earthRadB / sqrt( 1 - e2 * st * st )
    lon = geo%lon * deg2rad
    rho = r * cos(lat)
    ECR%xyz(1) = rho * cos(lon)
    ECR%xyz(2) = rho * sin(lon)
    ECR%xyz(3) = r * st
  end function H_t_To_ECR_Surf

!  -----     Private Procedures     ------------------------------------

  pure elemental subroutine Geod_To_Geoc_Support ( Geod, Z, Rho )
    use Constants, only: deg2rad
    class(h_v_geod), intent(in) :: Geod
    real(rg), intent(out) :: Z    ! Z coordinate of ECR
    real(rg), intent(out) :: Rho  ! sqrt(X^2+Y^2) of ECR
    real(rg) :: Ct   ! cos(lat)
    real(rg) :: D
    real(rg) :: Lat
    real(rg) :: St   ! Sin(lat)
    lat = geod%lat * deg2rad
    st = sin(lat)
    ct = cos(lat)
    d = 1.0_rg / sqrt(aquad*ct*ct + bquad*st*st)
!   d = 1.0_rg / sqrt(aquad-(aquad-bquad)*st*st)
    z = ( geod%v + bquad * d ) * st
    rho = ( geod%v + aquad * d ) * ct
  end subroutine Geod_To_Geoc_Support

  ! Construct an H_Geoc object from an H_t one
  pure elemental type(h_geoc) function H_t_to_H_Geoc ( H ) result ( G )
    type(h_t), intent(in) :: H
    g = h_geoc(h%lon,h%lat)
  end function H_t_to_H_Geoc

  ! Construct an H_Geod object from an H_t one
  pure elemental type(h_geod) function H_t_to_H_Geod ( H ) result ( G )
    type(h_t), intent(in) :: H
    g = h_geod(h%lon,h%lat)
  end function H_t_to_H_Geod

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Geolocation_0

! $Log$
! Revision 2.5  2016/02/24 01:18:21  vsnyder
! Add more coordinate conversions, add arithmetic for ECR.
!
! Revision 2.4  2016/01/23 02:46:21  vsnyder
! Add type-bound coordinate conversions
!
! Revision 2.3  2015/12/31 00:04:14  vsnyder
! Make V_t and H_V_t nonabstract
!
! Revision 2.2  2015/12/01 21:05:03  vsnyder
! Make H_t non-abstract
!
! Revision 2.1  2015/11/14 18:02:36  vsnyder
! Initial commit
!
