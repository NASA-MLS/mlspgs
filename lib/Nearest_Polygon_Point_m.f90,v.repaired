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
module Nearest_Polygon_Point_m
!=============================================================================

  implicit NONE
  private

  ! Find the point on a polygon's boundary that is nearest to a specified
  ! point, and report interpolating coefficients to its intersection.

  public :: Nearest_Polygon_Point, Nearest_Polygon_Point_Geo
  public :: Nearest_Polygon_Point_ZOT

  interface Nearest_Polygon_Point
    module procedure &
      & Nearest_Polygon_Point_Geo, &
      & Nearest_Polygon_Point_ZOT
  end interface Nearest_Polygon_Point

  ! Find the vertex of a polygon that is nearest to a specified point.

  public :: Nearest_Polygon_Vertex, Nearest_Polygon_Vertex_Geoc
  public :: Nearest_Polygon_Vertex_Geod, Nearest_Polygon_Vertex_ZOT

  interface Nearest_Polygon_Vertex
    module procedure Nearest_Polygon_Vertex_Geoc, &
      & Nearest_Polygon_Vertex_Geod, Nearest_Polygon_Vertex_H, &
      & Nearest_Polygon_Vertex_ZOT
  end interface Nearest_Polygon_Vertex

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  ! Find the point on the boundary of a polygon that is nearest to a point
  ! specified using longitude and latitude.  We don't report which polygon
  ! vertices are involved.
  subroutine Nearest_Polygon_Point_Geo ( Point, Polygon, P )
    use Geolocation_0, only: H_Geod, H_t, RG
    use Line_and_Polygon, only: Perpendicular_Intersection
    class(h_t), intent(in) :: Point       ! If the dynamic types of Point and
    class(h_t), intent(in) :: Polygon(:)  ! Polygon are not H_t, we assume
                                          ! they are both geodetic or both
                                          ! geocentric.
    class(h_t), intent(out) :: P          ! Nearest point

    logical :: Between     ! Perpendicular intersection is between Polygon(Which)
    real(rg) :: Cross(2)   ! Point of perpendicular intersection
    real(rg) :: D, D_min   ! Distance**2 from Point to edges (N), minimum so far
    integer :: I, J
    integer :: N(2)        ! N(1) is nearest Polygon vertex to Point,
                           ! N(2) is adjacent Polygon vertex
    integer :: S           ! Size(Polygon)

    select type ( point )
    type is ( h_geod )
      n(1) = nearest_polygon_vertex ( point, polygon )
    class default ! Either H_Geoc or H_t
      n(1) = nearest_polygon_vertex ( point, polygon )
    end select
    s = size(polygon)
    j = 0
    p%lon = polygon(n(1))%lon ! Assume no perpendicular intersection within 
    p%lat = polygon(n(1))%lat ! either incident edge
    d_min = huge(d_min)
    do i = 1, 2
      n(2) = 1 + mod(n(1)+j,s)     ! Adjacent vertex
      call perpendicular_intersection ( &
        & [ polygon(n)%lon%d, point%lon%d ], &
        & [ polygon(n)%lat, point%lat ], cross, between )
      if ( between ) then
        d = (point%lon%d-cross(1))**2 + (point%lat-cross(2))**2
        if ( d < d_min ) then
          p%lon%d = cross(1)
          p%lat = cross(2)
        end if
      end if
      j = s - 2
    end do

  end subroutine Nearest_Polygon_Point_Geo

  ! Find the point on the boundary of a polygon that is nearest to a point
  ! specified using ZOT coordinates.  We don't report which polygon
  ! vertices are involved.
  subroutine Nearest_Polygon_Point_ZOT ( Point, Polygon, P )
    use Line_and_Polygon, only: Perpendicular_Intersection
    use QTM_m, only: ZOT_t, RG
    type(ZOT_t), intent(in) :: Point
    class(ZOT_t), intent(in) :: Polygon(:)
    type(ZOT_t), intent(out) :: P        ! Nearest point

    logical :: Between     ! Perpendicular intersection is between Polygon(Which)
    real(rg) :: Cross(2)   ! Point of perpendicular intersection
    real(rg) :: D, D_min   ! Distance**2 from Point to edges (N), minimum so far
    integer :: I, J
    integer :: N(2)        ! N(1) is nearest Polygon vertex to Point,
                           ! N(2) is adjacent Polygon vertex
    integer :: S           ! Size(Polygon)

    n(1) = nearest_polygon_vertex ( point, polygon )
    s = size(polygon)
    j = 0
    p = polygon(n(1)) ! Assume no perpendicular intersection within either
                      ! incident edge
    d_min = huge(d_min)
    do i = 1, 2
      n(2) = 1 + mod(n(1)+j,s)     ! Adjacent vertex
      call perpendicular_intersection ( &
        & [ polygon(n)%x, point%x ], &
        & [ polygon(n)%y, point%y ], cross, between )
      if ( between ) then
        d = (point%x-cross(1))**2 + (point%y-cross(2))**2
        if ( d < d_min ) p = ZOT_t(cross(1),cross(2))
      end if
      j = s - 2
    end do

  end subroutine Nearest_Polygon_Point_ZOT

  ! Find the vertex of a polygon that is nearest to a point specified using
  ! geocentric longitude and latitude.  If there are equidistant ones, report
  ! the first one.
  integer function Nearest_Polygon_Vertex_Geoc ( Point, Polygon ) result ( N )
    use Geolocation_0, only: H_Geoc, H_t
    type(h_geoc), intent(in) :: Point
    class(h_t), intent(in) :: Polygon(:) ! Assumed also to be geocentric

    n = nearest_polygon_vertex ( point%h_t, polygon )

  end function Nearest_Polygon_Vertex_Geoc

  ! Find the vertex of a polygon that is nearest to a point specified using
  ! geodetic longitude and latitude.  If there are equidistant ones, report
  ! the first one.
  integer function Nearest_Polygon_Vertex_Geod ( Point, Polygon ) result ( N )
    use Geolocation_0, only: H_Geod, H_t, RG
    use Geometry, only: GeodToECRm
    type(h_geod), intent(in) :: Point
    class(h_t), intent(in) :: Polygon(:) ! Assumed also to be geodetic
    real(rg) :: D         ! Dot product of Pt_ECR and V_ECR
    real(rg) :: D_Max     ! Maximum value of D
    integer :: I
    real(rg) :: Pt_ECR(3) ! Cartesian (ECR) coordinates corresponding to Point
    real(rg) :: V_ECR(3)  ! Cartesian (ECR) coordinates corresponding to a vertex
    pt_ECR = GeodToECRm ( [ point%lat, point%lon%d, 0.0_rg ] )
    pt_ECR = pt_ECR / norm2(pt_ECR) ! Make it a unit vector
    d_max = -2
    do i = 1, size(polygon)
      v_ECR = GeodToECRm ( [ polygon(i)%lat, polygon(i)%lon%d, 0.0_rg ] )
      v_ECR = v_ECR / norm2(v_ECR)
      d = dot_product ( v_ECR, pt_ECR )
      if ( d > d_max ) then
        d_max = d
        n = i
      end if
    end do
  end function Nearest_Polygon_Vertex_Geod

  ! Find the vertex of a polygon that is nearest to a point specified using
  ! geocentric longitude and latitude.  If there are equidistant ones, report
  ! the first one.
  integer function Nearest_Polygon_Vertex_H ( Point, Polygon ) result ( N )
    use Geolocation_0, only: H_t, RG
    use Geometry, only: To_XYZ
    type(h_t), intent(in) :: Point      ! Assumed to be geocentric
    type(h_t), intent(in) :: Polygon(:) ! Assumed also to be geocentric
    real(rg) :: D         ! Dot product of Pt_ECR and V_ECR
    real(rg) :: D_Max     ! Maximum value of D
    integer :: I
    real(rg) :: Pt_ECR(3) ! Cartesian (ECR) coordinates corresponding to Point
    real(rg) :: V_ECR(3)  ! Cartesian (ECR) coordinates corresponding to a vertex
    pt_ECR = to_xyz ( point%lat, point%lon%d ) ! Result is a unit vector
    d_max = -2
    do i = 1, size(polygon)
      v_ecr = to_xyz ( polygon(i)%lat, polygon(i)%lon%d ) ! Result is a unit vector
      d = dot_product ( v_ECR, pt_ECR )
      if ( d > d_max ) then
        d_max = d
        n = i
      end if
    end do
  end function Nearest_Polygon_Vertex_H

  ! Find the vertex of a polygon that is nearest to a point specified using
  ! ZOT coordinates.  Distance is computed in the L-1 or "taxicab" norm.  If
  ! there are equidistant ones, report the first one.
  integer function Nearest_Polygon_Vertex_ZOT ( Point, Polygon ) result ( N )
    use QTM_m, only: ZOT_t, RG
    type(ZOT_t), intent(in) :: Point
    class(ZOT_t), intent(in) :: Polygon(:)
    real(rg) :: D         ! L-1 distance between Point and Polygon(i)
    real(rg) :: D_Max     ! Maximum value of D
    integer :: I
    d_max = 99
    do i = 1, size(polygon)
      d = abs(point%x - polygon(i)%x) + abs(point%y - polygon(i)%y)
      if ( d < d_max ) then
        d_max = d
        n = i
      end if
    end do
  end function Nearest_Polygon_Vertex_ZOT

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

end module Nearest_Polygon_Point_m

! $Log$
! Revision 2.3  2016/03/25 00:46:44  vsnyder
! Lon component now needs to acces its %d component
!
! Revision 2.2  2016/01/23 02:48:09  vsnyder
! Replace To_Cart with GeodToECRm
!
! Revision 2.1  2015/12/31 00:03:23  vsnyder
! Initial Commit
!
