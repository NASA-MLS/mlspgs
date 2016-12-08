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

  use Constants, only: Deg2Rad, Rad2Deg
  use Earth_Constants, only: Earth_Axis_Ratio, Earth_Axis_Ratio_Squared, &
    & EarthRadA, EarthRadB, Eccentricity_Sq, RP
  use MLSKinds, only: RG => R8 ! Kind of REAL geolocation variables

  implicit NONE

  public

  ! Users should get these from their original definitions, not from here
  private :: Deg2Rad, Rad2Deg
  private :: Earth_Axis_Ratio, Earth_Axis_Ratio_Squared, &
    & EarthRadA, EarthRadB, Eccentricity_Sq, RP

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  type, abstract :: Lat_t ! For latitude; can't instantiate it
    real(rg) :: D       ! Degrees
  contains
    procedure :: Geoc => Lat_t_to_GeocLat_t  ! Simple type conversion only
    procedure :: Geod => Lat_t_to_GeodLat_t  ! Simple type conversion only
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

  type :: H_t           ! For horizontal (lon/lat) grids, geoc/geod unspecified
    type(lon_t) :: Lon  ! Degrees
    real(rg) :: Lat     ! Degrees
  contains
    procedure :: ECR => H_t_To_ECR_Surf ! at the Earth's surface
    procedure, pass(geo) :: From_ECR => Surf_H_t_From_ECR ! useful because the
                                        ! desired dynamic result type (H_t) of
                                        ! ECR_to_Geoc cannot be used for
                                        ! dispatch
    procedure :: Geoc => H_t_to_H_Geoc  ! Simple type conversion only
    procedure :: Geod => H_t_to_H_Geod  ! Simple type conversion only
    procedure :: H_t_fields             ! A nonpolymorphic "self" component
                                        ! to avoid a defined output subroutine
    procedure :: Surf_ECR => H_t_To_ECR_Surf ! at the Earth's surface
  end type H_t

  type, extends(h_t) :: H_Geoc
  !  type(lon_t) :: Lon ! Longitude, degrees
  !  real(rg) :: Lat    ! Geocentric latitude, degrees
  contains
    procedure :: Geod => Geoc_To_Geod
  end type H_Geoc

  type, extends(h_t) :: H_Geod
  !  type(lon_t) :: Lon ! Longitude, degrees
  !  real(rg) :: Lat    ! Geodetic latitude, degrees
  contains
    procedure :: ECR => H_Geod_To_ECR_Surf
    procedure :: Geoc => Geod_To_Geoc
    procedure, pass(geo) :: From_ECR => Surf_Geod_From_ECR ! useful because the
                                        ! desired dynamic result type (H_geod)
                                        ! of ECR_to_Geod cannot be used for
                                        ! dispatch
    procedure, pass(geo) :: Surf_ECR => H_Geod_To_ECR_Surf
  end type H_Geod

  type :: V_t           ! For vertical grids of unspecified type
    real(rg) :: V       ! Geocentric (meters), geodetic (meters), or zeta
  end type V_t

  type, extends(v_t) :: V_geoc
  ! real(rg) :: V       ! geocentric height (meters from the Earth center)
  end type V_geoc

  type, extends(v_t) :: V_geod
  ! real(rg) :: V       ! geodetic height (meters above mean geoid)
  end type V_geod

  type, extends(v_t) :: V_zeta
  ! real(rg) :: V       ! Zeta -- geocentric latitude doesn't make sense if
                        ! the coordinate system is stacked
  end type V_zeta

  ! For irregular grids, or lines of sight, where H and V are expected to
  ! have the same extents
  type, extends(h_t) :: H_V_t
  ! type(lon_t) :: Lon  ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! Meters or zeta
  contains
    procedure :: H_v_t_fields           ! A nonpolymorphic "self" component
                                        ! to avoid a defined output subroutine
    procedure, pass(geo) :: Surf_From_ECR => Surf_H_v_t_From_ECR ! useful
                        ! because the desired dynamic result type (H_t) of
                        ! ECR_to_Geoc cannot be used for dispatch
    procedure, pass(geo) :: Surf_ECR => H_V_t_To_ECR_Surf ! At the Earth's surface
  end type H_V_t

  type, extends(h_v_t) :: H_Geoc_V_Geod
  ! type(lon_t) :: Lon  ! Degrees
  ! real(rg) :: Lat     ! Degrees Geocentric
  ! real(rg) :: V       ! geodetic height (meters above mean geoid)
  contains
    procedure :: ECR => Geoc_Geod_To_ECR
  end type H_Geoc_V_Geod

  type, extends(h_v_t) :: H_V_Geoc
  ! type(lon_t) :: Lon  ! Degrees
  ! real(rg) :: Lat     ! Degrees Geocentric
  ! real(rg) :: V       ! geocentric height (meters from the Earth center)
  contains
    procedure :: ECR => Geoc_To_ECR
    procedure :: GeodV => GeocV_To_GeodV  ! Would be nice to call this Geod,
                        ! but the result type of the one it would override is
                        ! H_Geoc and we don't want to make the results
                        ! polymorphic
    procedure, pass(geo) :: From_ECR => Geoc_From_ECR ! useful because the
                        ! desired dynamic result type (H_v_geoc) of ECR_to_Geoc
                        ! cannot be used for dispatch
  end type H_V_Geoc

  type, extends(h_v_t) :: H_V_Geod
  ! type(lon_t) :: Lon  ! Degrees
  ! real(rg) :: Lat     ! Degrees Geodetic
  ! real(rg) :: V       ! geodetic height (meters above mean geoid)
  contains
    procedure :: ECR => Geod_To_ECR
    procedure :: GeocV => GeodV_To_GeocV  ! Would be nice to call this Geoc,
                        ! but the result type of the one it would override is
                        ! H_Geoc and we don't want to make the results
                        ! polymorphic
    procedure, pass(geo) :: From_ECR => Geod_From_ECR ! useful because the
                        ! desired dynamic result type (H_v_geod) of ECR_to_Geod
                        ! cannot be used for dispatch
    procedure, pass(geo) :: Surf_ECR => H_V_Geod_To_ECR_Surf ! ECR at the surface
    procedure, pass(geo) :: Surf_From_ECR => Surf_H_v_Geod_From_ECR ! useful
                        ! because the desired dynamic result type (H_geod) of
                        ! ECR_to_Geod cannot be used for dispatch
  end type H_V_Geod

  type, extends(h_v_t) :: H_V_Zeta
  ! type(lon_t) :: Lon  ! Degrees
  ! real(rg) :: Lat     ! Degrees
  ! real(rg) :: V       ! Zeta
  end type H_V_Zeta

  type :: ECR_t
    real(rg) :: XYZ(3)  ! Cartesian ECR, meters, or a unit vector
  contains
    procedure :: Add => ECR_add
    generic :: operator(+) => Add
    procedure :: Cross_Norm => ECR_Cross_Norm ! Unit length
    generic :: operator(.CROSSNORM.) => Cross_Norm
    procedure :: Cross_Product => ECR_Cross
    generic :: operator(.CROSS.) => Cross_Product
    procedure :: Dot_Product => ECR_Dot
    generic :: operator(.DOT.) => Dot_Product
    procedure :: Div => ECR_Divide ! All components by the same value
    procedure :: Div_Components => ECR_Divide_Components ! Separately
    generic :: operator(/) => Div, Div_Components
    procedure :: Geoc => ECR_to_Geoc
    procedure :: Geod => ECR_To_Geod_Fukushima
    procedure :: Geod_Surf => ECR_To_Geod_Surface
    generic :: Grad => Grad_Ellipsoid ! Unit gradient at surface
    generic :: Grad => Grad_Geoid     ! Unit gradient at surface
    procedure, pass(where) :: Grad_Ellipsoid => Ellipsoid_Gradient_RG
    procedure :: Grad_Geoid => Earth_Geoid_Gradient
    procedure :: Negate => ECR_Negate
    procedure :: Norm2 => ECR_Norm2
    procedure, pass(B) :: Scale_L => ECR_Scale_L
    procedure :: Scale_R => ECR_Scale_R
    procedure :: Subtract => ECR_Subtract
    generic :: operator(-) => Negate, Subtract
    generic :: operator(*) => Scale_L, Scale_R
  end type ECR_t

  type :: S_t        ! Descriptors of points along a line
    real(rg) :: S    ! Distance along Line(1) in the direction of Line(2)
  contains
    generic :: operator (*) => Scaled_Line_L, Scaled_Line_R
    procedure :: Scaled_Line_L            ! S_t%s * ECR_t => ECR_t
    procedure, pass(S) :: Scaled_Line_R   ! ECR_t * S_t%s => ECR_t
  end type S_t

  interface Cross
    module procedure ECR_Cross_Norm_Opt
  end interface

  interface Dot_Product
    module procedure ECR_Dot
  end interface

  interface Norm2
    module procedure ECR_Norm2
  end interface

! Uncomment these when we get past ifort 14 and 15

!   interface H_Geoc ! Construct H_Geoc object from an H_t one
!     module procedure H_t_to_H_Geoc
!   end interface
! 
!   interface H_Geod ! Construct H_Geod object from an H_t one
!     module procedure H_t_to_H_Geod
!   end interface

! *****     Private Constants     **************************************

! AQUAD   Square of major half axis for earth ellipsoid
! BQUAD   Square of minor half axis for earth ellipsoid

  real(rg), parameter, private :: AQUAD = EarthRadA**2 ! m**2
  real(rg), parameter, private :: BQUAD = EarthRadB**2 ! m**2

! *****     Procedures     *********************************************

contains

  pure elemental type(ECR_t) function Earth_Geoid_Gradient ( Where ) result ( Grad )
    class(ECR_t), intent(in) :: Where ! Where on the Geoid the gradient is desired
    real(rp), parameter :: A = EarthRadA, B = EarthRadB
    grad = ECR_t ( where%xyz / [ A, A, B ]**2 )
    grad = grad / grad%norm2()
  end function Earth_Geoid_Gradient

  pure type(ECR_t) function Ellipsoid_Gradient_RG ( Axes, Center, Where ) &
    & result ( Grad )
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(in) :: Center  ! Center of the ellipsoid
    class(ECR_t), intent(in) :: Where  ! Where on the Geoid the gradient is desired
    grad = ECR_t ( ( where%xyz - center%xyz ) / axes**2 )
    grad = grad / grad%norm2()
  end function Ellipsoid_Gradient_RG

  pure elemental type(ECR_t) function ECR_Add ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_Add = ECR_t ( a%xyz + b%xyz )
  end function ECR_Add

  pure elemental type(ECR_t) function ECR_Cross ( A, B )
    use Cross_m, only: Cross
    class(ECR_t), intent(in) :: A, B
    ECR_Cross%xyz = cross(a%xyz,b%xyz)
  end function ECR_Cross

  pure elemental type(ECR_t) function ECR_Negate ( A )
    class(ECR_t), intent(in) :: A
    ECR_Negate%xyz = -A%xyz
  end function ECR_Negate

  pure elemental type(ECR_t) function ECR_Cross_Norm ( A, B ) result ( ECR_Cross )
    use Cross_m, only: Cross
    class(ECR_t), intent(in) :: A, B
    ECR_Cross%xyz = cross(a%xyz,b%xyz,norm=.true.)
  end function ECR_Cross_Norm

  pure elemental type(ECR_t) function ECR_Cross_Norm_Opt ( A, B, Norm ) result ( ECR_Cross )
    use Cross_m, only: Cross
    class(ECR_t), intent(in) :: A, B
    logical, intent(in), optional :: Norm
    ECR_Cross%xyz = cross(a%xyz,b%xyz,norm)
  end function ECR_Cross_Norm_Opt

  ! Divide all the elements of the component by the same value
  pure elemental type(ECR_t) function ECR_Divide ( A, B )
    class(ECR_t), intent(in) :: A
    real(rg), intent(in) :: B
    ECR_Divide = ECR_t ( a%xyz / b )
  end function ECR_Divide

  ! Divide each element of the component by a different value
  pure type(ECR_t) function ECR_Divide_Components ( A, B )
    class(ECR_t), intent(in) :: A
    real(rg), intent(in) :: B(3)
    ECR_Divide_Components = ECR_t ( a%xyz / b )
  end function ECR_Divide_Components

  pure elemental real(rg) function ECR_Dot ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_Dot = dot_product ( a%xyz, b%xyz )
  end function ECR_Dot

  pure elemental real(rg) function ECR_Norm2 ( A )
    class(ECR_t), intent(in) :: A
    ECR_Norm2 = norm2(a%xyz)
  end function ECR_Norm2

  pure elemental type(ECR_t) function ECR_Scale_L ( A, B )
    real(rg), intent(in) :: A
    class(ECR_t), intent(in) :: B
    ECR_Scale_L = ECR_t ( a * b%xyz )
  end function ECR_Scale_L

  pure elemental type(ECR_t) function ECR_Scale_R ( A, B )
    class(ECR_t), intent(in) :: A
    real(rg), intent(in) :: B
    ECR_Scale_R = ECR_t ( a%xyz * b )
  end function ECR_Scale_R

  pure elemental type(ECR_t) function ECR_Subtract ( A, B )
    class(ECR_t), intent(in) :: A, B
    ECR_subtract = ECR_t ( a%xyz - b%xyz )
  end function ECR_Subtract

  pure elemental type(h_v_geoc) function ECR_To_Geoc ( ECR ) result ( Geoc )
    class(ECR_t), intent(in) :: ECR
    real(rg) :: Rho, Rho2
    rho2 = ECR%xyz(1)**2 + ECR%xyz(2)**2
    rho = sqrt(rho2)
    geoc%lon%d = rad2deg * atan2(ECR%xyz(2), ECR%xyz(1))
    geoc%lat = rad2deg * asin(ECR%xyz(3)/rho)
    geoc%v = sqrt(rho2 + ECR%xyz(3)**2)
  end function ECR_To_Geoc

  pure elemental type(h_v_geod) function ECR_To_Geod_Fukushima ( ECR ) result ( Geod )
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
    geod%lon%d = rad2deg * atan2 ( ECR%xyz(2), ECR%xyz(1) )

    s = sqrt ( ECR%xyz(1)**2 + ECR%xyz(2)**2 )
    if ( s < sqrt(tiny(1.0_rg)) ) then ! Too near the poles
      geod%lat = sign(90.0_rg,ECR%xyz(3))
      geod%v = abs(ECR%xyz(3)) - b
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

  pure elemental function ECR_To_Geod_Surface ( ECR ) result ( Geo )
    class(ECR_t), intent(in) :: ECR
    type(h_geod) :: Geo
    type(h_v_geoc) :: Geoc
    geoc = ecr%geoc()
    geo = geoc%geod()
  end function ECR_To_Geod_Surface

  pure elemental subroutine Geoc_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_v_geoc), intent(inout) :: Geo ! intent(out) is prohibited
    type(h_v_geoc) :: MyGeo ! Needed because assignment to polymorphic Geo is
                            ! prohibited
    myGeo = ECR%geoc()
    geo%lon = myGeo%lon
    geo%lat = myGeo%lat
    geo%v = myGeo%v
  end subroutine Geoc_From_ECR

  pure elemental type(ECR_t) function Geoc_Geod_To_ECR ( Geo, Norm ) &
    & result ( ECR )
    ! Convert geocentric coordinates at the Earth surface plus geodetic
    ! height to ECR
    class(h_geoc_v_geod), intent(in) :: Geo
    logical, intent(in), optional :: Norm ! Compute normalized ECR
    type(ECR_t) :: Surf
    surf = geo%h_t%ECR() ! ECR coordinates of geocentric coordinates at the
                         ! Earth geoid surface
    ECR = surf + geo%v * surf%grad()
    if ( present(norm) ) then
      if ( norm ) ECR = ECR/norm2(ECR)
    end if
  end function Geoc_Geod_To_ECR

  pure elemental type(ECR_t) function Geoc_To_ECR ( Geo, Norm ) result ( ECR )
    ! Convert geocentric coordinates to ECR coordinates in meters
    class(h_v_geoc), intent(in) :: Geo
    logical, intent(in), optional :: Norm ! Compute normalized ECR
    real(rg) :: CosLat, Lat, Lon
    logical MyNorm
    myNorm = .false.
    if ( present(norm) ) myNorm = norm
    lat = geo%lat * deg2rad
    lon = geo%lon%d * deg2rad
    cosLat = cos(lat)
    ECR%xyz(1) = cosLat * cos(lon)
    ECR%xyz(2) = cosLat * sin(lon)
    ECR%xyz(3) = sin(lat)
    if ( .not. myNorm ) ECR = ECR * geo%v
  end function Geoc_To_ECR

  pure elemental function Geoc_To_Geod ( Geo ) result ( Geod )
    class(h_geoc), intent(in) :: Geo
    type(h_geod) :: Geod
    real(rp), parameter :: F2 = Earth_Axis_Ratio_Squared ! (b/a)**2
    real(rg) :: Lat
    geod%lon = geo%lon
    lat = geo%lat * deg2rad
    geod%lat = atan2 ( sin(lat), f2 * cos(lat) ) * rad2deg
  end function Geoc_To_Geod

  pure elemental function Geoc_To_Geod_Lat ( Geo ) result ( Geod )
    class(geocLat_t), intent(in) :: Geo
    type(geodLat_t) :: Geod
    real(rp), parameter :: F2 = Earth_Axis_Ratio_Squared ! (b/a)**2
    real(rg) :: Lat
    lat = geo%d * deg2rad
    geod%d = atan2 ( sin(lat), f2 * cos(lat) ) * rad2deg
  end function Geoc_To_Geod_Lat

  pure elemental type(h_v_geod) function GeocV_To_GeodV ( Geo ) result ( Geod )
    class(h_v_geoc), intent(in) :: Geo
    type(ECR_t) :: ECR
  ! geod = geoc%ECR()%geod()   ! Prohibited (why?)
    ECR = geo%ECR()
    geod = ECR%geod()
  end function GeocV_To_GeodV

  pure elemental subroutine Geod_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_v_geod), intent(inout) :: Geo ! intent(out) is prohibited
    type(h_v_geod) :: MyGeo ! Needed because assignment to polymorphic Geo is
                            ! prohibited
    myGeo = ECR%geod()
    geo%lon = myGeo%lon
    geo%lat = myGeo%lat
    geo%v = myGeo%v
  end subroutine Geod_From_ECR

  pure elemental type(ECR_t) function Geod_To_ECR ( Geo, Norm ) result ( ECR )
    ! Convert geodetic coordinates to ECR coordinates in meters
    class(h_v_geod), intent(in) :: Geo
    logical, intent(in), optional :: Norm
    real(rg) :: MyLon
    real(rg) :: Rho
    call h_v_geod_to_geoc_support ( geo, ECR%xyz(3), rho )
    myLon = geo%lon%d * deg2rad
    ECR%xyz(1) = rho * cos(myLon)
    ECR%xyz(2) = rho * sin(myLon)
    if ( present(norm) ) then
      if ( norm ) ECR = ECR / norm2(ECR)
    end if
  end function Geod_To_ECR

  pure elemental type(h_geoc) function Geod_To_Geoc ( Geo ) result ( Geoc )
    class(h_geod), intent(in) :: Geo
    real(rp), parameter :: F2 = Earth_Axis_Ratio_Squared ! (b/a)**2
    real(rg) :: Lat
    geoc%lon = geo%lon
    lat = geo%lat * deg2rad
    geoc%lat = atan2 ( f2 * sin(lat), cos(lat) ) * rad2deg
  end function Geod_To_Geoc

  pure elemental type(GeocLat_t) function Geod_To_Geoc_Lat ( Geo ) result ( Geoc )
    class(GeodLat_t), intent(in) :: Geo
    real(rp), parameter :: F2 = Earth_Axis_Ratio_Squared ! (b/a)**2
    real(rg) :: Lat
    lat = geo%d * deg2rad
    geoc%d = atan2 ( f2 * sin(lat), cos(lat) ) * rad2deg
  end function Geod_To_Geoc_Lat

  pure elemental type(h_v_geoc) function GeodV_To_GeocV ( Geo ) result ( Geoc )
    ! Convert geodetic coordinates to Geocentric coordinates, including heights
    class(h_v_geod), intent(in) :: Geo
    real(rg) :: Rho
    real(rg) :: Z    ! Would be Z coordinate of ECR
    call h_v_geod_to_geoc_support ( geo, z, rho )
    geoc%lon = geo%lon
    geoc%lat = atan2 ( z, rho ) * rad2deg
    geoc%v = sqrt ( rho**2 + z**2 ) ! geocentric height
  end function GeodV_To_GeocV

  pure elemental type(ECR_t) function H_Geod_To_ECR_Surf ( Geo, Norm ) result ( ECR )
    ! Convert geodetic coordinates at the Earth surface to ECR coordinates
    class(h_geod), intent(in) :: Geo
    logical, intent(in), optional :: Norm
    ECR = geod_to_ECR_surf ( geo%lon%d, geo%lat, norm )
  end function H_Geod_To_ECR_Surf

  ! Nonpolymorphic result; useful to avoid writing a defined output routine
  pure elemental type(h_t) function H_t_Fields ( Geo )
    class(h_t), intent(in) :: Geo
    h_t_fields = h_t(lon=geo%lon, lat=geo%lat)
  end function H_t_Fields

  pure elemental type(ECR_t) function H_t_To_ECR_Surf ( Geo, Norm ) result ( ECR )
    ! Convert geocentric coordinates at the Earth surface to ECR coordinates
    class(h_t), intent(in) :: Geo
    logical, intent(in), optional :: Norm
    ECR = geoc_to_ECR_surf ( geo%lon%d, geo%lat, norm )
  end function H_t_To_ECR_Surf

  ! Construct an H_Geoc object from an H_t one
  pure elemental type(h_geoc) function H_t_to_H_Geoc ( Geo ) result ( G )
    class(h_t), intent(in) :: Geo
    g = h_geoc(geo%lon,geo%lat)
  end function H_t_to_H_Geoc

  ! Construct an H_Geod object from an H_t one
  pure elemental type(h_geod) function H_t_to_H_Geod ( Geo ) result ( G )
    class(h_t), intent(in) :: Geo
    g = h_geod(geo%lon,geo%lat)
  end function H_t_to_H_Geod

  pure elemental type(ECR_t) function H_V_Geod_To_ECR_Surf ( Geo, Norm ) result ( ECR )
    ! Convert geodetic coordinates at the Earth surface to ECR coordinates
    class(h_v_geod), intent(in) :: Geo
    logical, intent(in), optional :: Norm
    ECR = geod_to_ECR_surf ( geo%lon%d, geo%lat, norm )
  end function H_V_Geod_To_ECR_Surf

  ! Nonpolymorphic result; useful to avoid writing a defined output routine
  pure elemental type(h_v_t) function H_v_t_Fields ( Geo )
    class(h_v_t), intent(in) :: Geo
    h_v_t_fields = h_v_t(lon=geo%lon, lat=geo%lat, v=geo%v)
  end function H_v_t_Fields

  pure elemental type(ECR_t) function H_V_t_To_ECR_Surf ( Geo, Norm ) result ( ECR )
    ! Convert geodetic coordinates at the Earth surface to ECR coordinates
    class(h_v_t), intent(in) :: Geo
    logical, intent(in), optional :: Norm
    ECR = geoc_to_ECR_surf ( geo%lon%d, geo%lat, norm )
  end function H_V_t_To_ECR_Surf

  pure elemental type (GeocLat_t) function Lat_t_to_GeocLat_t ( Geo ) result ( G )
    class(lat_t), intent(in) :: Geo
    g = GeocLat_t(geo%d)
  end function Lat_t_to_GeocLat_t

  pure elemental type (GeodLat_t) function Lat_t_to_GeodLat_t ( Geo ) result ( G )
    class(lat_t), intent(in) :: Geo
    g = GeodLat_t(geo%d)
  end function Lat_t_to_GeodLat_t

  pure elemental function Scaled_Line_L ( S, Line ) result ( R )
    class(S_t), intent(in) :: S
    type(ECR_t), intent(in) :: Line
    type(ECR_t) :: R ! ifort 16.0.2 rejects type decl in the function header
    r = s%s * line
  end function Scaled_Line_L

  pure elemental function Scaled_Line_R ( Line, S ) result ( R )
    type(ECR_t), intent(in) :: Line
    class(S_t), intent(in) :: S
    type(ECR_t) :: R ! ifort 16.0.2 rejects type decl in the function header
    r = line * s%s
  end function Scaled_Line_R

  pure elemental subroutine Surf_Geod_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_geod), intent(inout) :: Geo ! intent(out) is prohibited
    type(h_v_geod) :: MyGeo
    myGeo = ECR%geod()
    geo%lon = myGeo%lon
    geo%lat = myGeo%lat
  end subroutine Surf_Geod_From_ECR

  pure elemental subroutine Surf_H_t_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_t), intent(inout) :: Geo ! intent(out) is prohibited
    geo%lon%d = rad2deg * atan2(ECR%xyz(2), ECR%xyz(1))
    geo%lat = rad2deg * asin(ECR%xyz(3)/sqrt(ECR%xyz(1)**2 + ECR%xyz(2)**2))
  end subroutine Surf_H_t_From_ECR

  pure elemental subroutine Surf_H_v_Geod_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_v_geod), intent(inout) :: Geo ! intent(out) is prohibited
    type(h_v_geod) :: MyGeo
    myGeo = ECR%geod()
    geo%lon = myGeo%lon
    geo%lat = myGeo%lat
    geo%v = 0
  end subroutine Surf_H_v_Geod_From_ECR

  pure elemental subroutine Surf_H_v_t_From_ECR ( ECR, Geo )
    class(ECR_t), intent(in) :: ECR
    class(h_v_t), intent(inout) :: Geo ! intent(out) is prohibited
    call geo%h_t%from_ECR ( ECR )
    geo%v = 0
  end subroutine Surf_H_v_t_From_ECR

!  =====     Private Procedures     ====================================

  pure elemental type(ECR_t) function Geoc_To_ECR_Surf ( Lon, Lat, Norm ) result ( ECR )
    ! Convert geocentric coordinates at the Earth surface to ECR coordinates
    real(rg), intent(in) :: Lon, Lat ! Degrees
    logical, intent(in), optional :: Norm
    real(rp), parameter :: E2 = Eccentricity_Sq
    real(rg) :: MyLat, MyLon, R, Rho, ST
    myLat = lat * deg2rad
    myLon = lon * deg2rad
    st = sin(myLat)
    ! Geocentric radius at geocentric latitude
    ! r = ab / sqrt ( a^2 cos^2(myLat) + b^2 sin^2(myLat) )
    r = earthRadB / sqrt( 1 - e2 * st * st )
    rho = r * cos(myLat)
    ECR%xyz(1) = rho * cos(myLon)
    ECR%xyz(2) = rho * sin(myLon)
    ECR%xyz(3) = r * st
    if ( present(norm) ) then
      if ( norm ) ECR = ECR / norm2(ECR)
    end if
  end function Geoc_To_ECR_Surf

  pure elemental type(ECR_t) function Geod_To_ECR_Surf ( Lon, Lat, Norm ) result ( ECR )
    ! Convert geodetic coordinates at the Earth surface to ECR coordinates
    real(rg), intent(in) :: Lon, Lat ! Degrees
    logical, intent(in), optional :: Norm
    real(rg) :: D      ! sqrt(a^2 ct^2 + b^2 st^2)
    real(rg) :: MyLon  ! longitude in radians
    real(rg) :: Rho    ! a^2/d cos(lat)
    real(rg) :: ST     ! sin(lat)
    call geod_to_geoc_support ( lat, d, rho, st )
    myLon = lon * deg2rad
    ecr%xyz(1) = rho * cos(myLon)
    ecr%xyz(2) = rho * sin(myLon)
    ecr%xyz(3) = bquad * d * st
    if ( present(norm) ) then
      if ( norm ) ECR = ECR / norm2(ECR)
    end if
  end function Geod_To_ECR_Surf

  pure elemental subroutine Geod_To_Geoc_Support ( Lat, D, Rho, St )
    real(rg), intent(in) :: Lat   ! Geodetic latitude, degrees
    real(rg), intent(out) :: D
    real(rg), intent(out) :: Rho  ! sqrt(X^2+Y^2) of ECR
    real(rg), intent(out) :: St   ! Sin(myLat)
    real(rg) :: Ct    ! cos(myLat)
    real(rg) :: MyLat ! Geodetic latitude, radians
    myLat = lat * deg2rad
    st = sin(myLat)
    ct = cos(myLat)
    d = 1.0_rg / sqrt(aquad*ct*ct + bquad*st*st)
!   d = 1.0_rg / sqrt(aquad-(aquad-bquad)*st*st)
    rho = aquad * d * ct
  end subroutine Geod_To_Geoc_Support

  pure elemental subroutine H_V_Geod_To_Geoc_Support ( Geod, Z, Rho )
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
  end subroutine H_V_Geod_To_Geoc_Support

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
! Revision 2.17  2016/12/08 02:00:57  vsnyder
! Correct F***ingcd ../lib! units mistake (degrees vs. radians)
!
! Revision 2.16  2016/11/04 22:50:53  vsnyder
! Make as many type-bound procedures as possible elemental
!
! Revision 2.15  2016/11/03 20:40:28  vsnyder
! Add ECR_To_Geod_Surface
!
! Revision 2.14  2016/10/17 23:06:56  vsnyder
! Intent(OUT) polymorphic dummy argument of pure subroutine is prohibited
!
! Revision 2.13  2016/09/24 02:08:11  vsnyder
! Add Surf_ECR to H_t, H_V_t, and extensions other than H_V_Zeta
!
! Revision 2.12  2016/09/23 01:34:03  vsnyder
! Add H_v_t_fields to H_v_t, for the same reason as adding H_t_fields
!
! Revision 2.11  2016/09/22 20:32:11  vsnyder
! Add H_t_fields to type H_t.  This is a subterfuge to access essentially
! the entire object from a polymorphic version of it, so we don't need to
! write defined I/O subroutines.  Ugh.
!
! Revision 2.10  2016/04/16 02:02:48  vsnyder
! Add Negate, type to represent position on a line
!
! Revision 2.9  2016/03/25 00:24:29  vsnyder
! Change type of Lon component from real(rg) to type(lon_t).  Add
! Div_Components binding to ECR_t and add it to Operator(/) generic.
! Change the Cross generic to ECR_Cross_Norm_Opt because operator
! generics cannot have optional arguments.  Add Norm option to
! Geoc_Geod_To_ECR, Geoc_to_ECR, Geod_To_ECR, Geod_To_ECR_Surf, and
! H_t_To_ECR_Surf.
!
! Revision 2.8  2016/03/02 21:48:02  vsnyder
! Add H_Geoc_V_Geod and ellipsoid gradient functions
!
! Revision 2.7  2016/03/02 19:23:49  vsnyder
! Add From_ECR bindings to H_t, H_Geod, and H_V_Geod -- H_Geoc inherits from
! H_t.  Add Surf_From_ECR bindings to H_V_t and H_V_Geod.  Add NORM2 binding
! to ECR_t.  Add NORM2 generic interface.  Make ECR_Divide, ECR_Scale_L, and
! ECR_Scale_R elemental.
!
! Revision 2.6  2016/02/25 21:15:55  vsnyder
! Comment out generics for H_Geoc and H_Geod because ifort 15 doesn't like
! them.  Add a bunch of type-conversion bindings.  Move some USE statements.
! Make stuff used from Constants Earth_Constants private.
!
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
