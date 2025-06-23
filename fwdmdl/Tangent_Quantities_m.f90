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

  ! The minor-frame line of sight (LOS) is represented by two vectors C and U
  ! that define lines of the form C + s * U, and therefore also define
  ! plane-projected ellipses that are the intersections of the Earth with the
  ! planes defined by the lines of sight and the Earth center, calculate the
  ! positions of the geodetic tangents or Earth intersections.

  ! For each tangent or intersection, compute
  ! (1) its geodetic angle in degrees from the Equator in the plane-projected
  !     ellipse (PhiTan),
  ! (2) the square of the minor axis of the plane-projected ellipse, and
  ! (3) the tangent height, i.e., the geodetic height above the plane-projected
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
                                & LOS, EarthRadC_sq, PhiTan, TanHt )

    use Constants, only: Rad2Deg
    use ForwardModelConfig, only: ForwardModelConfig_T
    use Earth_Constants, only: EarthRadA ! major axis, meters
    use Geolocation_0, only: ECR_t, Norm2, RG
    use Geometry, only: Orbit_Plane_Minor_Axis_sq, XZ_to_Geod
    use Intrinsic, only: l_mifRadC, l_phiTan, l_tngtGeodAlt
    use Line_And_Ellipsoid_m, only: Exact_Line_Nearest_Ellipsoid, &
                                  & Line_And_Ellipsoid
    use Path_Representation_m, only: Path_t
    use VectorsModule, only: CloneVectorQuantity, CreateVectorValue, RV, &
                           & Vector_T, VectorValue_T

    integer, intent(in) :: MAF
    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    type(vector_t), intent(in) :: FwdModelIn, FwdModelExtra

    type(Path_t), intent(in) :: LOS(:)          ! LOS as C' + s U'

    type(vectorValue_t), intent(out) :: EarthRadC_sq ! Square of  minor axis of
                                                ! plane-projected ellipse, m^2
    type(vectorValue_t), intent(out) :: PhiTan  ! Degrees from equator,
                                                ! geodetic, w.r.t.
                                                ! plane-projected ellipse
    type(vectorValue_t), intent(out) :: TanHt   ! Geodetic height, above plane-
                                                ! projected ellipse, meters

    real(rv), parameter :: A = EarthRadA ! In case kind(earthRadA) /= RV
    type(ECR_t) :: Equator, Normal, Pole ! of plane-projected ellipse
    real(rv) :: H ! TanHt proxy
    integer :: I
    real(rv) :: Incline
    type(ECR_t), allocatable :: Intersections(:) ! of LOS with Earth, if any
    type(ECR_t) :: Line(2)               ! C and U vectors for a line
    real(rv) :: P ! PhiTan proxy
    real(rg) :: R ! sqrt(x**2 + y**2), radius of unit normal vector projected
                  ! into circle of latitude at tangent point
    type(ECR_t) :: Tangent
    real(rv) :: X, Z ! Cartesian coordinates w.r.t. plane-projected ellipse

    ! Get quantities
    call cloneVectorQuantity ( earthRadC_sq, GetQuantity ( l_mifRadC ) )
    call cloneVectorQuantity ( phiTan, GetQuantity ( l_phiTan ) )
    call cloneVectorQuantity ( tanHt, GetQuantity ( l_tngtGeodAlt ) )

    ! Create values components
    call createVectorValue ( earthRadC_sq, "EarthRadC_sq" )
    call createVectorValue ( phiTan, "PhiTan" )
    call createVectorValue ( tanHt, "TanHt" )

    do i = 1, size(LOS)
      ! Get input LOS as two ECR vectors, and replace U with U'.
      line(1) = LOS(i)%lines(1,1)                              ! C
      line(2) = LOS(i)%lines(2,1) / LOS(i)%lines(2,1)%norm2()  ! U'
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
      incline = 90.0 - rad2deg * atan2 ( normal%xyz(3), r )
      ! (Minor axis)**2 of ellipse in plane containing LOS and Earth center.
      earthRadC_sq%values(i,MAF) = &
        & orbit_plane_minor_axis_sq ( incline )
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

    function GetQuantity ( What )
      use ForwardModelVectorTools, only: GetQuantityForForwardModel
      type(vectorValue_T), pointer :: GetQuantity
      integer, intent(in) :: What ! L_whatever
      getQuantity => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=what, config=fwdModelConf, &
        & instrumentModule=fwdModelConf%signals(1)%instrumentModule )
    end function GetQuantity

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
! Revision 2.7  2017/10/31 17:36:15  vsnyder
! Change QTM path from vector quantity to Path_t
!
! Revision 2.6  2016/11/11 01:57:38  vsnyder
! Eliminate Incline argument because nobody used it.  Don't check LOS chans
! because that's done earlier, which also allows to eliminate Wrong_Shape
! internal subroutine.  Spiff some comments.
!
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
