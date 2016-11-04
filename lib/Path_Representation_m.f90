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
    ! tangent, and one from the tangent onward. The given line is the one
    ! from the instrument to the tangent, represented by
    ! Lines(1,1) + s * Lines(2,1). A line defined by
    ! Lines(1,2) + s * Lines(2,2) is produced, which is the continuation
    ! of Lines(:,1) after the tangent point if Lines(:,1) does not
    ! intersect the Earth reference ellipsoid, or the reflection of
    ! Lines(:,1) if it does intersect the Earth reference ellipsoid. 
    ! Lines(2,1) and Lines(2,2) are unit vectors.
    type(ECR_t) :: Lines(2,2)
    real(rg) :: SMax(2), SMin(2) ! Interesting intervals of S along Lines
    logical :: Ready = .false.   ! Lines(:,2), SMin and SMax have been computed.
  contains
    procedure :: Get_Path_Ready
    procedure :: New_Path
  end type Path_t

  type, public :: Value_t
    integer :: N   ! Position in array-element order to which V is germane. 
                   ! If it's a position in a rank-2 array, its array-element
                   ! order is computed by Element_Position in Array_Stuff. 
                   ! If there is none (which shouldn't happen), N is one and
                   ! V is zero.
    real(rg) :: V  ! Value at a point on or near the path
  end type Value_t

  type, public :: Flag_t
    integer :: N   ! Position in array-element order to which F is germane. 
                   ! If it's a position in a rank-2 array, its array-element
                   ! order is computed by Element_Position in Array_Stuff. 
                   ! If there is none (which shouldn't happen), N is one and
                   ! F is zero.
    logical :: F   ! Flag at a point on or near the path
  end type Flag_t

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

    use Geolocation_0, only: Norm2
    use Line_And_Ellipsoid_m, only: Line_And_Ellipsoid, Line_Nearest_Ellipsoid
    use Line_And_Plane_m, only: Line_Reflection

    class(path_t), intent(inout) :: Path

    type(ECR_t) :: Grad   ! Gradient to Earth reference ellipsoid at tangent point
    integer :: I
    type(ECR_t), allocatable :: Ints(:) ! Intersections with Earth reference ellipsoid
    real(rg) :: R         ! R = 1 => intersection, >= 1 => tangent
    real(rg) :: Tangent   ! S-value of tangent point
    real(rg) :: Tan_Dir   ! +/- 1; Tan_Dir * Lines(2,i) is directed from
                          ! Lines(1,i) toward the tangent point.
    real(rg), allocatable :: W(:) ! Where line and ellipsoid intersect

    if ( path%ready ) return

    ! Make Lines(2,1) an unit vector, to simplify later calculations.
    path%lines(2,1) = path%lines(2,1) / norm2(path%lines(2,1))

    ! Get the tangent or intersection point of Lines(:,1) with the Earth
    ! reference ellipsoid.
    ! Lines(1,1) + s * Lines(2,1) describes the line before the tangent point.
    ! Lines(1,2) + s * Lines(2,2) describes the line after the tangent point.
    ! We could, in principle, divide the line into an arbitrary number of
    ! segments, to account (approximately) for refraction.
    call line_nearest_ellipsoid ( path%lines(:,1), tangent, r )
    if ( r >= 1 ) then ! No intersection, Lines(:,2) is colinear with Lines(:,1)
      path%lines(1,2) = path%lines(1,1) + tangent * path%lines(2,1)
      path%lines(2,2) = path%lines(2,1) ! Parallel to incident line
    else               ! Compute reflection direction
      ! Assume lines(1,1) is not inside the ellipsoid
      call line_and_ellipsoid ( path%lines(:,1), ints, w )
      ! If path%lines(:,1) is not tangent to the Earth's surface, there are
      ! two intersections.  Choose the one closest to path%lines(1,1).
      i = minloc(norm2(ints - path%lines(1,1)),1)
      path%lines(1,2) = ints(i) ! The continuation starts at that intersection
      tangent = w(i)
      grad = path%lines(1,2)%grad() ! Gradient to Earth reference ellipsoid at
      ! the intersection Compute Lines(2,2) such that Lines(2,2) is at the
      ! same angle from Grad as Lines(2,1), but on the opposite side of the
      ! tangent from Lines(2,1).
      call line_reflection ( path%lines(2,1), grad, path%lines(2,2) )
      deallocate ( w )
    end if
    tan_dir = sign(1.0_rg,tangent)
    if ( tan_dir >= 0 ) then
      ! Direction of Lines(2,1) is from Lines(1,1) toward tangent, therefore
      ! direction of Lines(2,2) is from tangent toward +infinity
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
! Revision 2.2  2016/11/04 01:26:32  vsnyder
! Spiff some comments
!
! Revision 2.1  2016/10/26 01:20:19  vsnyder
! Initial commit
!
