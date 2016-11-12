! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Path_Representation_m
!=============================================================================

  ! Types to represent information about the path of integration of the
  ! radiative-transfer equation, and information along it.

  use Geolocation_0, only: ECR_t, RG

  implicit NONE
  private

  type, public :: Path_t
    ! A path is represented by two lines, one from the instrument to the
    ! tangent or intersection, and one from there onward. The given line is the
    ! one from the instrument to the tangent or intersection, represented by
    ! Lines(1,1) + s * Lines(2,1). A line defined by
    ! Lines(1,2) + s * Lines(2,2) is produced, which is the continuation
    ! of Lines(:,1) after the tangent point if Lines(:,1) does not
    ! intersect the Earth reference ellipsoid, or the reflection of
    ! Lines(:,1) if it does intersect the Earth reference ellipsoid. 
    ! Lines(2,1) and Lines(2,2) are made to be unit vectors.
    type(ECR_t) :: Lines(2,2)
    real(rg) :: SMax(2), SMin(2) ! Interesting intervals of S along Lines
    logical :: Ready = .false.   ! Lines(:,2), SMin and SMax have been computed.
  contains
    procedure :: Get_Path_Ready
    procedure :: New_Path
  end type Path_t

  type, public :: Value_t
    integer :: N = 0  ! Serial number of QTM vertex to which V is germane.
    integer :: NP = 1 ! Subscript of vertex near the path to which V is germane.
    integer :: NZ = 1 ! Vertical subscript of layer below vertex to which V is germane.
    real(rg) :: V = 0 ! Value at a point on or near the path
  end type Value_t

  type, public :: Flag_t
    integer :: N = 0  ! Serial number of QTM vertex to which F is germane.
    integer :: NP = 1 ! Subscript of vertex near the path to which F is germane.
    integer :: NZ = 1 ! Vertical subscript of layer below vertex to which F is germane.
    logical :: F = .false. ! Flag at a point on or near the path
  end type Flag_t

  type, public :: Facets_and_Vertices_t
    integer, allocatable :: Facets(:)   ! Indices in QTM_Tree%Q
    integer, allocatable :: Vertices(:) ! Indices in QTM_Tree%Geo_In
  end type Facets_and_Vertices_t

  public :: Path_Continuation

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Get_Path_Ready ( Path )

    ! Given a line defined by a point in ECR, and a vector in ECR parallel
    ! to that line, compute the position of the tangent on that line, and
    ! the extents of the interesting parts of the line.

    ! The given line is Path%Lines(1,1) + s * Path%Lines(2,1).
    ! Path%Lines(2,1) is made an unit vector here.  A line defined by
    ! Path%Lines(1,2) + s * Path%Lines(2,2) is produced, which is the
    ! continuation of Path%Lines(:,1) after the tangent point if
    ! Path%Lines(:,1) does not intersect the Earth reference ellipsoid, or
    ! the reflection of Path%Lines(:,1) if it does intersect the Earth
    ! reference ellipsoid.  Path%Lines(2,2) is an unit vector.

    ! The tangent or surface-intersection point is Path%Lines(1,2).

    class(path_t), intent(inout) :: Path

    real(rg) :: Tangent   ! S-value of tangent point
    real(rg) :: Tan_Dir   ! +/- 1; Tan_Dir * Lines(2,i) is directed from
                          ! Lines(1,i) toward the tangent point.

    if ( path%ready ) return

    ! Make Path%Lines(1:2) an unit vector to simplify later calculations.
    ! Get the tangent or intersection point of Path%Lines(:,1) with the Earth
    ! reference ellipsoid.
    ! Path%Lines(1,1) + s * Path%Lines(2,1) describes the line before the
    ! tangent point.
    ! Path%Lines(1,2) + s * Path%Lines(2,2) describes the line after the
    ! tangent point.
    ! We could, in principle, divide the line into an arbitrary number of
    ! segments, to account (approximately) for refraction.
    call path_continuation ( path%lines, tangent )
    tan_dir = sign(1.0_rg,tangent)
    if ( tan_dir >= 0 ) then
      ! Direction of Path%Lines(2,1) is from Path%Lines(1,1) toward tangent,
      ! therefore direction of Path%Lines(2,2) is from tangent toward +infinity
      path%sMin(1) = -sqrt(huge(0.0_rg)) ! Allow intersections before Lines(1,1)
      path%sMax(1) = tangent             ! Disallow intersections after tangent
      path%sMin(2) = 0                   ! Disallow intersections before tangent
      path%sMax(2) = sqrt(huge(0.0_rg))  ! Allow intersections arbitrarily far away
    else ! sMax(1) < 0
      ! Direction of Lines(2,1) is from Lines(1,1) away from tangent, therefore
      ! direction of Lines(2,2) is from -infinity toward tangent
      path%sMin(1) = tangent             ! Disallow intersections before tangent
      path%sMax(1) = sqrt(huge(0.0_rg))  ! Allow intersections arbitrarily far away
      path%sMin(2) = -sqrt(huge(0.0_rg)) ! Allow intersections arbitrarily far away
      path%sMax(2) = 0                   ! Disallow intersections after tangent
    end if

  end subroutine Get_Path_Ready

  subroutine New_Path ( Default_Initialized )

    class(path_t), intent(out) :: Default_Initialized

  end subroutine New_Path

  subroutine Path_Continuation ( Lines, Tangent )

    ! Calculate the path continuation after the tangent point or intersection.
    ! Lines(1,1) + s * Lines(2,1) describes the line before the tangent point.
    ! Lines(1,2) + s * Lines(2,2) describes the line after the tangent point.
    ! We could, in principle, divide the line into an arbitrary number of
    ! segments, to account (approximately) for refraction.

    use Geolocation_0, only: Norm2
    use Line_And_Ellipsoid_m, only: Line_And_Ellipsoid, &
      & Exact_Line_Nearest_Ellipsoid
    use Line_And_Plane_m, only: Line_Reflection
    type(ECR_t), intent(inout) :: Lines(2,2)
    real(rg), intent(out) :: Tangent   ! S-value of tangent point

    type(ECR_t) :: Grad           ! Gradient to Earth reference ellipsoid
                                  ! at tangent point
    real(rg) :: H                 ! Tangent height above Earth surface,
                                  ! Lines is an intersection if H < 0.
    integer :: I
    type(ECR_t), allocatable :: Ints(:) ! Intersections with Earth reference ellipsoid
    real(rg), allocatable :: W(:) ! Where line and ellipsoid intersect

    ! Make Lines(2,1) an unit vector, to simplify later calculations.
    lines(2,1) = lines(2,1) / lines(2,1)%norm2()

    ! Get the tangent or intersection point of Lines(:,1) with the Earth
    ! reference ellipsoid.
    call exact_line_nearest_ellipsoid ( lines(:,1), tangent, h=h )
    if ( h >= 0 ) then ! No intersection, Lines(:,2) is colinear with Lines(:,1)
      lines(1,2) = lines(1,1) + tangent * lines(2,1)
      lines(2,2) = lines(2,1) ! Parallel to incident line
    else               ! Compute reflection direction
      ! Assume lines(1,1) is not inside the ellipsoid
      call line_and_ellipsoid ( lines(:,1), ints, w )
      ! If lines(:,1) is not tangent to the Earth's surface, there are
      ! two intersections.  Choose the one closest to lines(1,1).
      i = minloc(norm2(ints - lines(1,1)),1)
      lines(1,2) = ints(i) ! The continuation starts at that intersection
      tangent = w(i)
      grad = lines(1,2)%grad() ! Gradient to Earth reference ellipsoid at
      ! the intersection Compute Lines(2,2) such that Lines(2,2) is at the
      ! same angle from Grad as Lines(2,1), but on the opposite side of the
      ! tangent from Lines(2,1).
      call line_reflection ( lines(2,1), grad, lines(2,2) )
      deallocate ( ints, w )
    end if
  end subroutine Path_Continuation

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

end module Path_Representation_m

! $Log$
! Revision 2.7  2016/11/12 01:32:23  vsnyder
! Add NP and NZ components to Value_t and Flag_t, and default initialize
!
! Revision 2.6  2016/11/11 01:46:41  vsnyder
! Add Facets_and_Vertices_t
!
! Revision 2.5  2016/11/08 02:04:56  vsnyder
! Use Exact_Line_Nearest_Ellipsoid
!
! Revision 2.4  2016/11/07 23:50:55  vsnyder
! Remove unused variable declaration
!
! Revision 2.3  2016/11/07 23:48:52  vsnyder
! Make Path_Continuation public
!
! Revision 2.2  2016/11/04 01:26:32  vsnyder
! Spiff some comments
!
! Revision 2.1  2016/10/26 01:20:19  vsnyder
! Initial commit
!
