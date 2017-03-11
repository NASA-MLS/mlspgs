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

  subroutine QTM_Facets_Under_Path ( Path, QTM, F_and_V )
    ! Find the facets of QTM that are under Path, and list their vertices.
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_t, H_Geoc, H_Geod, Lon_t, RG
    use Nearest_Lines_m, only: Nearest_Lines
    use Path_Representation_m, only: Facets_and_Vertices_t, Path_t
    use QTM_m, only: Stack_t
    use Where_m, only: Where

    type(path_t), intent(inout) :: Path
    type(QTM_tree_t), intent(in) :: QTM
    type(facets_and_vertices_t), intent(out) :: F_and_V

    logical :: Adjacent(QTM%n_in) ! "Vertex I is adjacent to the path"
    type(h_t) :: Centroid    ! of a facet
    integer :: Facet         ! Facet found to be under the path
    type(ECR_t) :: Grad(2)   ! Normal to the centroid of the polygon
    integer :: I, J
    logical :: Kept(QTM%n)   ! List of kept facets
    type(ECR_t) :: NP        ! Unit normal to the plane containing the path and
                             ! the center of the Earth; it's the same for both
                             ! halves of the path
    integer :: Q(QTM%n)      ! Queue of facets to investigate
    integer :: QH, QT        ! Queue head, tail
    type(stack_t) :: Stack   ! For faster QTM searches
    logical :: Parallel      ! Path and normal to centroid of polygon are parallel
    real(rg) :: S1, S2       ! Positions of points on path and normal to centroid
                             ! of the polygon that are nearest to each other.
    integer :: V             ! Serial number of a vertex of a facet
    logical :: Visited(QTM%n) ! "Facet(i) has been visited"

    call path%get_path_ready ! Compute tangent or intersection point, and the
                             ! continuation of the path therefrom.

    ! Compute the unit normal to the plane containing the path and the
    ! center of the Earth.  Both parts of the path are in the same plane.

! ifort 17.0.0.098 doesn't like this:
!   np = path%lines(1,1) .crossnorm. path%lines(2,1)
    np = path%lines(1,1)%cross_norm(path%lines(2,1))

    ! Is the tangent point or intersection above the polygon?
    ! We expect this frequently to be true.

    facet = qtm%find_facet ( path%lines(1,2), stack )

    if ( facet <= 0 ) then ! Tangent or intersection is not above the polygon.
      ! Check whether a point, on either part of the path, that is nearest to
      ! the gradient at the centroid of the polygon, is above the polygon.  If
      ! it's a nadir view (i.e., parallel to grad), which we expect to be rare
      ! (or nonexistent), we won't find any adjacent facets, so don't bother
      ! making that a special case.
      grad(1) = QTM%Polygon_Geo_Centroid%ECR()
      grad(2) = grad(1)%grad() ! Unit gradient
      do j = 1, 2
        call nearest_lines ( path%lines(:,j), grad, s1, s2, parallel )
        facet = qtm%find_facet ( path%lines(1,j) + s1 * path%lines(2,j), stack )
        if ( facet > 0 ) exit
      end do
    end if

    if ( facet <= 0 ) then
      ! Neither the tangent or intersection point, nor a point on the path
      ! nearest to the gradient at the polygon centroid, is above the polygon. 
      ! Check whether the point on the path that is nearest to the gradient
      ! from the centroid of any facet above a facet.  It's not obvious that
      ! this method saves anything compared to the next method.
o:    do j = 1, 2 ! Both parts of the path
        do i = 1, QTM%n
          associate ( Q => QTM%Q(i) )
            if ( Q%depth == QTM%level ) then
              ! Search tree node represents a finest-refinement facet
              centroid = h_t( lon_t(sum(Q%geo%lon%d) / 3), sum(Q%geo%lat) / 3 )
              grad(1) = centroid%ECR()
              grad(2) = grad(1)%grad()
              call nearest_lines ( path%lines(:,j), grad, s1, s2, parallel )
              ! If parallel, then s1 == 0
              facet = qtm%find_facet ( path%lines(1,j) + s1 * path%lines(2,j), stack )
              if ( facet > 0 ) exit o
             end if
          end associate
        end do
      end do o
    end if

    ! There is no facet such that the point on the path nearest to the
    ! gradient at the centroid of that facet is above some facet.  Search
    ! for facets for which the path crosses an edge.  We don't need to check
    ! both parts of the path separately because their tracks on the surface
    ! of the Earth are on the same great circle.
    do i = 1, QTM%n
      if ( facet > 0 ) exit
      if ( path_crosses_facet(i) ) facet = i
    end do

    if ( facet <= 0 ) then ! Path apparently does not cross the polygon
      allocate ( f_and_v%facets(0), f_and_v%vertices(0) )
      return
    end if

    kept = .false.     ! We haven't kept any facets yet
    kept(facet) = .true.  ! except the starting facet
    visited = .false.  ! We haven't visited any facets yet
    adjacent = .false. ! We don't know of any adjacent vertices yet
    qh = 0             ! Queue of unvisited facets is empty
    qt = 0

    do
      if ( .not. visited(facet) ) then ! Don't process a facet more than once.
        visited(facet) = .true.
        if ( kept(facet) ) then
          do j = 1, 3 ! All three vertices of the facet are adjacent to the path
            v = QTM%q(facet)%ser(j)
            adjacent(v) = .true.
            ! Put facets adjacent to V on the queue
            do i = 1, QTM%n_adjacent_in(v)
              facet = QTM%adjacent_in(i,v) ! facet is a temp here
              if ( .not. visited(facet) ) then ! put unvisited facet on the queue
                qt = qt + 1
                q(qt) = facet
              end if
            end do
          end do
        end if
      end if
      qh = qh + 1
      if ( qh > qt ) exit ! Queue is empty
      facet = q(qh)
      ! We don't need to check both parts of the path because their tracks
      ! on the surface of the Earth are on the same great circle.
      kept(facet) = path_crosses_facet ( facet )
    end do

    ! Fill Facets and Vertices arrays
    f_and_v%nf = count(kept)
    f_and_v%nv = count(adjacent)
    allocate ( f_and_v%facets(f_and_v%nf), f_and_v%vertices(f_and_v%nv) )
    call where ( kept, f_and_v%facets )
    call where ( adjacent, f_and_v%vertices )

  contains

    logical function Path_Crosses_Facet ( F )

      integer, intent(in) :: F       ! Facet index in QTM%Q(:)

      integer :: I
      type(h_geoc) :: Geoc           ! Coordinates of a vertex
      type(h_geod) :: Geod           ! Coordinates of a vertex
      type(ECR_t) :: V(3)            ! Coordinates of facet vertices

      ! Get unit vectors to vertices of facet F, which are represented in
      ! QTM%Q by longitude and latitude.
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
      if ( path_intersects_facet_edge ( v(1), v(2) ) ) return
      if ( path_intersects_facet_edge ( v(1), v(3) ) ) return
      if ( path_intersects_facet_edge ( v(2), v(3) ) ) return
      path_crosses_facet = .false.

    end function Path_Crosses_Facet

    logical function Path_Intersects_Facet_Edge ( V1, V2 )

      type(ECR_t), intent(in) :: V1, V2 ! Unit vectors to vertices of facet

      real(rg) :: D            ! First |NP X NE|, then V1 .dot. V2
      type(ECR_t) :: NE        ! Unit normal to the plane containing edges
                               ! (V1,V2) of facet F and the center of the Earth
      type(ECR_t) :: N3        ! Normal to NP and NE, and therefore in the
                               ! intersection of the planes normal to them, and
                               ! therefore the vector to one of the
                               ! intersections of the path's great circle with
                               ! the edge's great circle.

      ! We can't use the triple-cross-product formula here because we want NE
      ! to be the unit vector in the direction of V1 X V2.
      ne = v1 .crossnorm. v2   ! Unit normal to the plane of the edge
      n3 = np .cross. ne       ! NP and NE are unit vectors here.
      d = n3%norm2()           ! |sine of the angle between NP and NE|.
      path_intersects_facet_edge = d > sqrt(huge(1.0_rg))
      ! If path and edge are too close to parallel, i.e., if the angle between
      ! NP and NE is too small, let another edge decide.  Maybe the path
      ! crosses an adjacent facet instead of the one currently of interest.
      if ( .not. path_intersects_facet_edge ) return

      n3 = n3 / d     ! Make N3 an unit vector
      d = v1 .dot. v2 ! Cosine of angle between facet vertices
      ! Use ABS in the check here because N3 might point to the other side
      ! of the Earth from the facet, depending on the ordering of vertices,
      ! which oughtn't make a difference.  If the angle to the point where
      ! the great circle of the path crosses the great circle of the facet
      ! edge, to both vertices of the edge, is less than the angle between
      ! the vertices (i.e., the cosine is greater), the path crosses the edge
      ! between the vertices.
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
! Revision 2.5  2017/03/11 00:49:06  vsnyder
! Assign values to NF and NV components of F_and_V
!
! Revision 2.4  2016/11/11 01:47:24  vsnyder
! Use Facets_and_Vertices_t instead of schlepping them indivicually
!
! Revision 2.3  2016/11/07 21:21:46  vsnyder
! Numerous simplifications and corrections
!
! Revision 2.2  2016/11/04 01:25:33  vsnyder
! Detect intersection at V2.  Don't care about vertex order
!
! Revision 2.1  2016/11/04 01:19:28  vsnyder
! Initial commit
!
