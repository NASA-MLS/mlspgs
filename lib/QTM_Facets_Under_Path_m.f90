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
module QTM_Facets_Under_Path_m
!=============================================================================

  ! Find the facets of a QTM that are under a path.

  ! This is useful to avoid computations above the surface, such as
  ! interpolating temperature to the entire Zeta Gauss-Legendre grid, and
  ! hydrostatic equilibrium, at vertices of other facets, that are not under
  ! the path, and for which therefore the coefficients to interpolate from
  ! their vertices to the path will be zero.

  implicit NONE
  private

  public :: QTM_Facets_Under_Path

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine QTM_Facets_Under_Path ( Path, QTM, Facets, Vertices )
    ! Find the facets of QTM that are under Path, and list their vertices.
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_t, H_Geoc, H_Geod, H_V_Geoc, Lon_t, RG
    use Nearest_Lines_m, only: Nearest_Lines
    use Path_Representation_m, only: Path_t
    use QTM_m, only: Stack_t
    type(path_t), intent(inout) :: Path
    type(QTM_tree_t), intent(in) :: QTM
    integer, allocatable, intent(out) :: Facets(:)   ! Indices of QTM%Q(:)
    integer, allocatable, intent(out) :: Vertices(:) ! Indices of QTM%Geo_In(:)

    type(h_t) :: Centroid   ! of a facet
    integer :: Facet        ! Facet found to be under the path
    type(h_t) :: Geo        ! lon/lat coordinates of a surface point
    type(h_v_geoc) :: Geoc3 ! Geocentric coordinates of an ECR point
    type(ECR_t) :: Grad(2)  ! Normal to the centroid of the polygon
    integer :: I, J
    integer, allocatable :: I_Temp(:) ! To reallocate Facets and Vertices
    logical :: Keep
    integer :: NF           ! Number of facets the path crosses
    type(ECR_t) :: NP       ! Unit normal to the plane containing the path and
                            ! the center of the Earth
    integer :: NV           ! Number of vertices adjacent to the path
    integer :: Q(QTM%n)     ! Queue of facets to investigate
    integer :: QH, QT       ! Queue head, tail
    type(stack_t) :: Stack  ! For faster QTM searches
    logical :: Parallel     ! Path and normal to centroid of polygon are parallel
    real(rg) :: S1, S2      ! Positions of points on path and normal to centroid
                            ! of the polygon that are nearest to each other.
    integer :: V            ! Serial number of a vertex of a facet
    logical :: Visited(QTM%n) ! "Facet(i) has been visited"

    call path%get_path_ready ! Compute tangent or intersection point

    grad(1) = QTM%Polygon_Geo_Centroid%ECR()
    grad(2) = grad(1)%grad() ! Unit gradient

    call nearest_lines ( path%lines(:,1), grad, s1, s2, parallel )

    if ( parallel ) then ! Must be a nadir view.  Use the facet below
                         ! the path reference point, if there is one
      facet = qtm%find_facet ( path%lines(1,1), stack )
      if ( facet <= 0 ) then ! Path reference point is not above the polygon
        allocate ( facets(0), vertices(0) )
      else
        allocate ( facets(1), vertices(3) )
        facets(1) = facet               ! Only the facet facet
        vertices(:3) = QTM%Q(facet)%ser ! Serial numbers of the only facet
      end if
      return
    end if

    facet = qtm%find_facet ( path%lines(1,1) + s1 * path%lines(2,1), stack )

    if ( facet <= 0 ) then
      ! Point on path nearest to gradient at polygon centroid is not above the
      ! polygon. Use a brute-force method: find the point on the path that is
      ! nearest to the gradient from the centroid of each facet until a point
      ! on the path is found that is above a facet.
      do i = 1, QTM%n
        associate ( Q => QTM%Q(i) )
          if ( Q%depth == QTM%level ) then
            ! Search tree node represents a finest-refinement facet
            centroid = h_t( lon_t(sum(Q%geo%lon%d) / 3), sum(Q%geo%lat) / 3 )
            grad(1) = centroid%ECR()
            grad(2) = grad(1)%grad()
            call nearest_lines ( path%lines(:,1), grad, s1, s2, parallel )
            ! If parallel, s1 == 0
            facet = qtm%find_facet ( path%lines(1,1) + s1 * path%lines(2,1), stack )
            if ( facet > 0 ) exit
           end if
        end associate
      end do
    end if

    ! ifort 17.0.0.098 doesn't like this:
!   np = path%lines(1,1) .crossnorm. path%lines(2,1)
    np = path%lines(1,1)%cross_norm(path%lines(2,1))

b:  if ( facet <= 0 ) then ! The "hunt for a centroid" method didn't work.
      ! There is no facet such that the point on the path nearest to the
      ! gradient at the centroid of that facet is above some facet.  Search
      ! for facets for which the path crosses an edge.
      do facet = 1, QTM%n
        if ( path_crosses_facet(facet) ) exit b
      end do
      facet = 0
    end if b

    if ( facet <= 0 ) then ! Path apparently does not cross the polygon
      allocate ( facets(0), vertices(0) )
      return
    end if

    keep = .true. ! Keep the facet
    visited = .false.
    qh = 0
    qt = 0
    allocate ( facets(QTM%n), vertices(QTM%n_in) )
    nf = 0
    vertices(1:3) = QTM%q(facet)%ser ! Serial numbers of vertices
    nv = 0

    do
      visited(facet) = .true.
      if ( keep ) then
        nf = nf + 1
        facets(nf) = facet
        call insert_vertices ( facet )
        do j = 1, 3
          v = QTM%q(facet)%ser(j)
          ! Put facets adjacent to V on the queue
          do i = 1, QTM%n_adjacent_in(v)
            facet = QTM%adjacent_in(i,v) ! facet is a temp here
            if ( .not. visited(facet) ) then ! put facet on the queue
              qt = qt + 1
              q(qt) = facet
            end if
          end do
        end do
      end if
      qh = qh + 1
      if ( qh > qt ) exit ! Queue is empty
      facet = q(qh)
      keep = path_crosses_facet ( facet )
    end do

    ! Adjust Facets and Vertices arrays to the correct sizes
    if ( nf /= size(facets) ) then
      allocate ( i_temp(nf) )
      i_temp(:nf) = facets(:nf)
      call move_alloc ( i_temp, facets )
    end if
    if ( nv /= size(vertices) ) then
      allocate ( i_temp(nv) )
      i_temp(:nv) = vertices(:nv)
      call move_alloc ( i_temp, vertices )
    end if

  contains

    subroutine Insert_Vertices ( F )
      integer, intent(in) :: F ! Facet index
      integer :: I
      do i = 1, 3
        if ( any(vertices(:nv) == QTM%q(f)%ser(i)) ) cycle
        nv = nv + 1
        vertices(nv) = QTM%q(f)%ser(i)
      end do
    end subroutine Insert_Vertices

    logical function Path_Crosses_Facet ( F )
      integer, intent(in) :: F ! Facet index in QTM%Q(:)
      integer :: I
      type(h_geoc) :: Geoc     ! Coordinates of a vertex
      type(h_geod) :: Geod     ! Coordinates of a vertex
      type(ECR_t) :: V(3)      ! Coordinates of facet vertices

      ! Get unit vectors to facet vertices
      if ( QTM%geodetic ) then
        do i = 1, 3
          geod = QTM%q(f)%geo(i)%geod()
          v(i) = geod%ECR()
          v(i) = v(i) / v(i)%norm2()
        end do
      else
        do i = 1, 3
          geoc = QTM%q(f)%geo(i)%geoc()
          v(i) = geoc%ECR()
          v(i) = v(i) / v(i)%norm2()
        end do
      end if

      path_crosses_facet = .true.
      if ( path_intersects_facet_edge ( f, v(1), v(2) ) ) return
      if ( path_intersects_facet_edge ( f, v(1), v(3) ) ) return
      if ( path_intersects_facet_edge ( f, v(2), v(3) ) ) return
      path_crosses_facet = .false.

    end function Path_Crosses_Facet

    logical function Path_Intersects_Facet_Edge ( F, V1, V2 )
      integer, intent(in) :: F          ! Facet index in QTM%Q(:)
      type(ECR_t), intent(in) :: V1, V2 ! Unit vectors to vertices of facet F

      real(rg) :: D            ! V1 .dot. V2
      type(ECR_t) :: NE        ! Unit normal to the plane containing edges
                               ! (V1,V2) of facet F and the center of the Earth
      type(ECR_t) :: N3        ! Normal to NP and NE, and therefore in the
                               ! intersection of the planes normal to them.

      ne = v1 .crossnorm. v2
      n3 = np .cross. ne ! NP and NE are unit vectors here, so n3%norm2() is
                         ! the sine of the angle between them
      path_intersects_facet_edge = n3%norm2() < sqrt(huge(1.0_rg))
      if ( path_intersects_facet_edge ) return ! path and edge are in same plane

      d = v1 .dot. v2 ! Cosine of angle between facet vertices
      ! Use ABS in the check here because N3 might point to the other side
      ! of the Earth from the facet, depending on the ordering of vertices,
      ! which oughtn't make a difference.
      path_intersects_facet_edge = abs(n3 .dot. v1) >= d .and. &
                                 & abs(n3 .dot. v2) >= d

    end function Path_Intersects_Facet_Edge

  end subroutine QTM_Facets_Under_Path

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

end module QTM_Facets_Under_Path_m

! $Log$
! Revision 2.2  2016/11/04 01:25:33  vsnyder
! Detect intersection at V2.  Don't care about vertex order
!
! Revision 2.1  2016/11/04 01:19:28  vsnyder
! Initial commit
!
