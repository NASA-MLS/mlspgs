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
module Tangent_Quantities_m
!=============================================================================

  ! Given vectors C and U that define lines of sight of the form C + s * U,
  ! and therefore also define plane-projected ellipses that are the
  ! intersections of the Earth with the planes defined by the lines of sight
  ! and the Earth center, calculate the positions of the geodetic tangents
  ! or Earth intersections.

  ! Produce new LOS representations where C' gives the geodetic tangent
  ! location, and U' is an unit vector.  In the LOS representation C' + s * U',
  ! the tangent is at s == 0.

  ! For each tangent or intersection, compute
  ! (1) its geodetic angle in degrees from the Equator in the plane-projected
  !     ellipse,
  ! (2) the inclination of that plane,
  ! (3) the square of the minor axis of the plane-projected ellipse, and
  ! (4) the tangent height, i.e., the geodetic height above the plane-projected
  !     ellipse (zero, not negative, for an intersection).

  implicit NONE

  private

  public :: Tangent_Quantities

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Tangent_Quantities ( MAF, FwdModelConf, FwdModelIn, FwdModelExtra, &
                                & EarthRadC_sq, Incline, LOS, PhiTan, TanHt )

    use Constants, only: Rad2Deg
    use ForwardModelConfig, only: ForwardModelConfig_T
    use Earth_Constants, only: EarthRadA ! major axis, meters
    use Geolocation_0, only: ECR_t, Norm2, RG
    use Geometry, only: Orbit_Plane_Minor_Axis_sq, XZ_to_Geod
    use Intrinsic, only: l_mifLOS, l_mifRadC, l_orbitInclination, l_phiTan, &
                       & l_tngtGeodAlt
    use Line_And_Ellipsoid_m, only: Exact_Line_Nearest_Ellipsoid, &
                                  & Line_And_Ellipsoid
    use VectorsModule, only: CloneVectorQuantity, RV, Vector_T, VectorValue_T

    integer, intent(in) :: MAF
    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    type(vector_t), intent(in) :: FwdModelIn, FwdModelExtra

    type(vectorValue_t), intent(out) :: EarthRadC_sq ! m^2
    type(vectorValue_t), intent(out) :: Incline ! Degrees, geocentric
    type(vectorValue_t), intent(out) :: LOS     ! LOS as C' + s U'
    type(vectorValue_t), intent(out) :: PhiTan  ! Degrees, geodetic, w.r.t.
                                                ! plane-projected ellipse
    type(vectorValue_t), intent(out) :: TanHt   ! Geodetic height, above plane-
                                                ! projected ellipse, meters

    real(rv), parameter :: A = EarthRadA ! In case kind(earthRadA) /= RV
    type(ECR_t) :: Equator, Normal, Pole ! of plane-projected ellipse
    real(rv) :: H ! TanHt proxy
    integer :: I
    type(ECR_t), allocatable :: Intersections(:) ! of LOS with Earth, if any
    type(ECR_t) :: Line(2) ! C vector, U' vector
    real(rv) :: P ! PhiTan proxy
    real(rg) :: R ! sqrt(x**2 + y**2), radius of unit normal vector projected
                  ! into circle of latitude at tangent point
    type(ECR_t) :: Tangent
    real(rv) :: X, Z ! Cartesian coordinates w.r.t. plane-projected ellipse

    ! Get quantities
    call cloneVectorQuantity ( earthRadC_sq, GetQuantity ( l_mifRadC ) )
    call cloneVectorQuantity ( incline, GetQuantity ( l_orbitInclination ) )
    call cloneVectorQuantity ( LOS, GetQuantity ( l_mifLOS ) )
    call cloneVectorQuantity ( phiTan, GetQuantity ( l_phiTan ) )
    call cloneVectorQuantity ( tanHt, GetQuantity ( l_tngtGeodAlt ) )

    ! Verify the quantities have compatible templates to LOS, and if
    ! they are compatible, create values components for them.
    call check_compatible_and_create ( earthRadC_sq, "EarthRadC_sq" )
    call check_compatible_and_create ( incline, "Incline" )
    call check_compatible_and_create ( phiTan, "PhiTan" )
    call check_compatible_and_create ( tanHt, "TanHt" )

    ! Check that LOS has at least 6 channels (two ECR quantities).
    ! We don't care if the others have more than one (maybe we should).
    if ( LOS%template%noChans < 6 ) call wrongShape ( LOS, "LOS", 6 )

    do i = 1, LOS%template%noSurfs
      ! Get input LOS as two ECR vectors, and replace U with U'.
      line(1) = ECR_t(LOS%value3(1:3,i,MAF))      ! C
      line(2) = ECR_t(LOS%value3(4:6,i,MAF)) / &
                norm2(LOS%value3(4:6,i,MAF))      ! U'
      ! Compute intersection of LOS with Earth surface, if any, or tangent point.
      call line_and_ellipsoid ( line, intersections )
      select case ( size(intersections) )
      case ( 1 ) ! One intersection, so the tangent is at the Earth's surface.
        tangent = intersections(1)
      case ( 2 ) ! Two intersections.  Assume Line(2) points from the
        ! instrument toward the tangent point, but Line(1) could be opposite
        ! the tangent point from the instrument, not at the instrument or
        ! nearer to it than the tangent point.  Therefore, first use Tangent
        ! as a temp to represent a point hopefully "behind" the instrument. 
        ! Assume s = -10 * A is big enough to get "behind" the instrument,
        ! and use that instead of -huge(1.0_rg) so there will be some digits
        ! in the difference between that point and the intersections, and
        ! there won't be an overflow while computing norm2.
        tangent = line(1) - ( 10 * A ) * line(2)
        ! Choose the intersection nearest the instrument.
!       ifort 14 didn't like this:
!       tangent = intersections(minloc(norm2 ( tangent - intersections ), 1 ))
        call Choose_Nearest( Tangent, Intersections )
      case default ! zero => No intersection, compute the tangent point.
        ! Hopefully, if we asked to get the tangent height here, it would be
        ! the same value as the one we get from XZ_to_Geod below.
        call exact_line_nearest_ellipsoid ( line, p=tangent )
      end select

      ! Unit normal to plane containing LOS and Earth center is C .crossnorm. U.
      normal = line(1) .crossnorm. line(2)
      r = norm2 ( normal%xyz(1:2) ) ! Radius of unit normal in plane parallel
                                    ! to the Equator.
      ! Inclination of plane containing LOS and Earth center = 90 - atan(Z / R).
      ! Maybe we need to worry here about ascending vs. descending.
      incline%values(i,MAF) = 90.0 - rad2deg * atan2 ( normal%xyz(3), r )
      ! (Minor axis)**2 of ellipse in plane containing LOS and Earth center.
      earthRadC_sq%values(i,MAF) = &
        & orbit_plane_minor_axis_sq ( incline%values(i,MAF) )
      ! Unit vector along the intersection of the plane-projected ellipse and
      ! the Equator = North_Pole .cross. Normal = (0,0,1) .cross. Normal.
      equator = ECR_t ( [ -normal%xyz(2), normal%xyz(1), 0.0_rg ] ) / r
      ! X coordinate of tangent point w.r.t. plane-projected ellipse is the
      ! projection of C' onto its equator.
      x = tangent .dot. equator
      ! Compute the unit "pole" of the plane-projected ellipse, which would be
      ! the north pole (0,0,1) if the inclination were 90 degrees.
      pole = normal .crossnorm. equator
      ! Z coordinate of tangent point w.r.t. the plane-projected ellipse is
      ! the projection of C' onto its pole.
      z = pole .dot. tangent
      ! Compute geodetic tangent coordinates w.r.t. the plane-projected ellipse,
      ! radians and meters above the ellipse.
      call xz_to_geod ( x, z, a, sqrt(earthRadC_sq%values(i,MAF)), p, h )
      if ( size(intersections) > 0 ) h = 0 ! Don't use negative height
      tanHt%values(i,MAF) = h
      p = rad2deg * p ! XZ_to_Geod produced angle in radians
      ! Put PhiTan in the correct quadrant.  The signs of phiTan and the Z
      ! component of the tangent are the same, so we don't need to worry
      ! about North vs. South.  The case that needs correcting is when the
      ! signs of phiTan and the Z component of the LOS direction differ.
      if ( line(2)%xyz(3) * p < 0.0 ) p = sign(180.0_rg - abs(p), p )
      ! Put p in 0 ... 360 instead of -180 ... 180.
      if ( p < 0.0 ) p = 360.0 + p
      phiTan%values(i,MAF) = p
    end do

  contains

    subroutine Choose_Nearest ( Tangent, Intersections )
      ! Choose the intersection nearest the instrument.
      type(ECR_t), intent(out) :: Tangent
      type(ECR_t), intent(in)  :: Intersections(:)
      ! Internal variables
      type(ECR_t), dimension(size(intersections)) :: diff
      real, dimension(size(intersections))        :: n2
      integer, dimension(1)                       :: k
      ! Executable
      diff = tangent - intersections
      n2 = norm2 ( diff )
      k = minloc( n2 )
      tangent = intersections( k(1) )
    end subroutine Choose_Nearest

    subroutine Check_Compatible_and_Create ( V, Name )
      ! Check whether the template of V is compatible to the template of LOS,
      ! and if so create V's values component.
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use MoreMessage, only: MLSMessage
      use QuantityTemplates, only: QuantitiesAreCompatible
      use VectorsModule, only: CreateVectorValue
      type(vectorValue_t), intent(in) :: V
      character(*), intent(in) :: Name
      if ( quantitiesAreCompatible ( LOS%template, PhiTan%template, &
        & differentTypeOK=.true., differentChansOK=.true. ) ) then
        call createVectorValue ( V, Name )
        return
      end if
      if ( LOS%template%name == 0 ) then
        if ( v%template%name == 0 ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Quantities LOS and " // name // " are incompatible." )
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Quantities LOS and %s are incompatible.", v%template%name )
        end if
      else
        if ( v%template%name == 0 ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Quantities %s and " // name // " are incompatible.", &
            & LOS%template%name )
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Quantities %s and %s are incompatible.", &
            & [ LOS%template%name, v%template%name ] )
        end if
      end if
    end subroutine Check_Compatible_and_Create

    function GetQuantity ( What )
      use ForwardModelVectorTools, only: GetQuantityForForwardModel
      type(vectorValue_T), pointer :: GetQuantity
      integer, intent(in) :: What ! L_whatever
      getQuantity => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=what, config=fwdModelConf, &
        & instrumentModule=fwdModelConf%signals(1)%instrumentModule )
    end function GetQuantity

    subroutine WrongShape ( V, Name, Shape )
      use MLSMessageModule, only: MLSMSG_Error
      use MoreMessage, only: MLSMessage
      type(vectorValue_t), intent(in) :: V
      character(*), intent(in) :: Name
      integer, intent(in) :: Shape
      if ( v%template%name == 0 ) then
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Quantity " // name // " does not have %d ECR components.", shape )
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Quantity %s does not have %d ECR components.", &
          & [ v%template%name, shape ] )
      end if
    end subroutine WrongShape

  end subroutine Tangent_Quantities

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

end module Tangent_Quantities_m

! $Log$
! Revision 2.5  2016/11/09 00:34:13  vsnyder
! Don't put first line reference at the tangent point.  Use .crossnorm. for
! cross products that should return unit vectors.  Specify intents for
! arguments of Choose_Nearets.  Correct an error message.
!
! Revision 2.4  2016/08/20 00:54:17  vsnyder
! Correct a typo in a comment
!
! Revision 2.3  2016/07/27 23:47:36  vsnyder
! Remove unreferenced USE
!
! Revision 2.2  2016/06/13 21:03:01  vsnyder
! Polish up some comments
!
! Revision 2.1  2016/06/07 18:45:57  pwagner
! First commit
!
