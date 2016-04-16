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
module Generate_QTM_m
!=============================================================================

  ! Generate a QTM inside a specified polygon.  Put the unique ZOT coordinates
  ! into a one-dimensional array, and remember their serial numbers in the QTM
  ! structure.  The ZOT projection is a projection of a sphere.  Therefore, the
  ! definition of "inside" a polygon is ambiguous, so an additional point that
  ! is defined to be inside the polygon is required.

  use Geolocation_0, only: Lon_t, H_t, RG
  use QTM_m, only: QK, ZOT_t

  implicit NONE
  private

  public :: Cross_Meridian, Dump, Dump_QTM, Dump_QTM_Tree, Generate_QTM, &
          & Get_QTM_Lats, Polygon_To_ZOT, QTM_node_t, QTM_Tree_t, ZOT_t

  ! One facet of a QTM:
  type :: QTM_Node_t
    logical :: Leaf = .false.  ! A leaf node
    integer(qk) :: QID = 0     ! QTM ID
    integer :: Son(0:3) = 0    ! Sub facet subscripts in the tree
    integer :: XN = 0, YN = 0  ! Which son is xNode, yNode?  Central node is 0,
                               ! pole node is 6 - ( xn + yn )
    type(ZOT_t) :: Z(3)        ! ZOT coordinates of vertices
    integer :: ZOT_N(3) = 0    ! Serial numbers of ZOT coordinates of vertices;
                               ! zero if the vertex is not in a polygon.
  end type QTM_Node_t

  type :: QTM_Tree_t
    integer :: Level = 6       ! Level of QTM refinement
    integer(qk) :: N           ! Number of used elements in Q
    type(ZOT_t) :: In = ZOT_t(-999,-999) ! A point defined to be inside the
                               ! polygon. This is needed because the concept of
                               ! "inside a polygon" is ambiguous on a sphere.
    type(h_t) :: In_Geo = h_t(lon_t(-999),-999) ! QTM_Tree_t%In as longitude
                               ! and latitude
    logical :: In_or_Out       ! True iff PnPoly reports QTM_Tree_t%In is
                               ! inside the polygon, from the point of view of
                               ! the ZOT projection.  A point is considered to
                               ! be inside if PnPoly returns the same value.
    integer(qk) :: N_In = 0    ! Number of vertices inside the polygon
    type(h_t), allocatable :: Polygon_Geo(:)
    type(ZOT_t), allocatable :: Polygon_ZOT(:)
    logical, allocatable :: Ignore_Edge(:) ! Ignore edge (i,1+mod(i,n)) because
                               ! it joins two aliased points on a southern-
                               ! hemisphere meridian at mod(lon,90)=0 and there
                               ! is an edge antiparallel to it and colinear
                               ! with it.
    type(QTM_node_t), allocatable :: Q(:)
  contains
    procedure :: Find_Facet_Geo
    procedure :: Find_Facet_QID
    procedure :: Find_Facet_ZOT
    generic :: Find_Facet => Find_Facet_Geo
    generic :: Find_Facet => Find_Facet_QID
    generic :: Find_Facet => Find_Facet_ZOT
  end type QTM_Tree_t

  interface Dump
    module procedure Dump_QTM, Dump_QTM_Tree
  end interface Dump

  ! Initial size to allocate for QTM_Tree_t%Q if it's not allocated
  integer, public :: Default_Initial_Size = 100

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  integer function Find_Facet_Geo ( QTM_Trees, Geo, Stack ) result ( F )
    ! Return the subscript in QTM_Trees of the facet containing the longitude
    ! and latitude specified by Geo.  The refinement level is QTM_Trees%Level.
    ! If a leaf is reached before running out of digits, the result is
    ! the negative of the subscript of the last vertex visited.
    use QTM_m, only: Geo_to_ZOT, Stack_t, ZOT_t
    class(QTM_Tree_t), intent(in) :: QTM_Trees
    class(h_t), intent(in) :: Geo
    type(stack_t), intent(inout), optional :: Stack ! to make subsequent search faster
    type(stack_t) :: My_Stack ! in case Stack is not present
    type(ZOT_t) :: Zot
    zot = geo_to_zot ( geo )
    if ( present(stack) ) then
      f = qtm_trees%find_facet ( zot%qtm_encode(qtm_trees%level, stack) )
    else
      f = qtm_trees%find_facet ( zot%qtm_encode(qtm_trees%level, my_stack) )
    end if
  end function Find_Facet_Geo

  integer function Find_Facet_QID ( QTM_Trees, QID ) result ( F )
    ! Return the subscript in QTM_Trees of the facet identified by QID.
    ! If a leaf is reached before running out of digits, the result is
    ! the negative of the subscript of the last vertex visited.
    use QTM_m, only: High_Bit_Index, QK
    class(QTM_Tree_t), intent(in) :: QTM_Trees
    integer(qk), intent(in) :: QID
    integer :: H, Next  ! High bit index, next tree node to visit
    integer :: QD       ! One two-bit digit of QID

    f = 1
    if ( qid < 8 ) return ! Bogus QID
    h = high_bit_index(QID)
    ! The number of iterations is bounded by the refinement level.
    do while ( h > 0 )
      h = h - 2
      qd = iand(int(shiftr(qid,h)),3)
      next = QTM_Trees%Q(f)%son(qd)
      if ( next == 0 ) exit
      f = next
    end do
    if ( QTM_Trees%Q(f)%QID /= QID ) f = -f
  end function Find_Facet_QID

  integer function Find_Facet_ZOT ( QTM_Trees, ZOT, Stack ) result ( F )
    ! Return the subscript in QTM_Trees of the facet containing ZOT.
    ! The refinement level is QTM_Trees%Level.
    ! If a leaf is reached before running out of digits, the result is
    ! the negative of the subscript of the last vertex visited.
    use QTM_m, only: Stack_t, ZOT_t
    class(QTM_Tree_t), intent(in) :: QTM_Trees
    type(ZOT_t), intent(in) :: ZOT
    type(stack_t), intent(inout), optional :: Stack ! to make subsequent search faster
    type(stack_t) :: My_Stack ! in case Stack is not present
    if ( present(stack) ) then
      f = qtm_trees%find_facet ( zot%qtm_encode(qtm_trees%level, stack) )
    else
      f = qtm_trees%find_facet ( zot%qtm_encode(qtm_trees%level, my_stack) )
    end if
  end function Find_Facet_ZOT

  subroutine Generate_QTM ( QTM_Trees, ZOT )
    ! Generate a QTM refined to QTM_Trees%Level (in 1..QTM_Depth) within a
    ! Polygon defined by QTM_Trees%Polygon_Geo or QTM_Trees%Polygon_ZOT.  A
    ! point inside the polygon, specified by QTM_Trees%In_Geo or QTM_Trees%In,
    ! is needed because the concept of "inside a polygon" is ambiguous on a
    ! sphere.

    ! Whether a QTM vertex is within a polygon is determined within the ZOT
    ! projection.  Depending upon how edges of facets are defined on the
    ! Earth, a point determined to be within (outwith) a polygon in the ZOT
    ! projection might be outwith (within) the polygon on the surface of the
    ! Earth if it's near an edge.

    ! The ZOT coordinates of the vertices are in ZOT.

    use PnPoly_m, only: PnPoly
    use QTM_m, only: Geo_to_ZOT, QTM_Decode, QTM_Depth, Stack_t, ZOT_t

    type(QTM_Tree_t), intent(inout) :: QTM_Trees
    type(ZOT_t), intent(out), allocatable :: ZOT(:) ! ZOT coordinates of QTM
                                   ! vertices that are within a polygon

    integer :: H(2,hash_size()) ! Assume QK = kind(0) for now
    integer :: Hemisphere ! 2 = north, 3 = south
    integer :: L          ! min(QTM_Trees%Level,QTM_Depth)
    integer :: Octant     ! Octant being refined.
    integer :: Quadrant   ! Quadrant in a hemisphere, 0..3.
    type(QTM_node_t), allocatable :: Q_Temp(:)
    type(stack_t) :: S    ! Put here so we don't need a new one in every
                          ! instance of Add_QTM_Vertex_To_Tree.  Speeds up
                          ! QTM_Decode in Add_QTM_Vertex_To_Tree if it's
                          ! not new for every instance.
    real(rg) :: X, Y      ! Disambiguated ZOT coordinates
    type(ZOT_t), allocatable :: Z_Temp(:) ! For reallocating ZOT

    QTM_Trees%n = 0 ! In case setup is wrong
    if ( .not. allocated(QTM_Trees%polygon_ZOT) ) then
      if ( .not. allocated(QTM_Trees%polygon_geo) ) return ! without doing anything
      call polygon_to_ZOT ( QTM_Trees )
    end if

    if ( QTM_Trees%in%x < -10 .or. QTM_Trees%in%y < -10 ) then
      if ( QTM_Trees%in_geo%lon%d < -500 .or. QTM_Trees%in_geo%lat < -500 ) return
      QTM_Trees%in = geo_to_ZOT(QTM_Trees%in_geo)
    end if
    l = min(QTM_Trees%level,QTM_Depth)
    h = 0

    x = merge(QTM_Trees%in%x,abs(QTM_Trees%in%x),abs(QTM_Trees%in%y)/=1.0_rg)
    y = merge(QTM_Trees%in%y,abs(QTM_Trees%in%y),abs(QTM_Trees%in%x)/=1.0_rg)

    QTM_Trees%in_or_out = pnPoly ( x, y, QTM_Trees%polygon_ZOT%x, &
                                 & QTM_Trees%polygon_ZOT%y ) >= 0

    if ( .not. allocated(QTM_Trees%Q) ) allocate ( QTM_Trees%Q(default_initial_size) )
    allocate ( zot(default_initial_size) )
    ! Create the top eight octants in the tree
    QTM_Trees%n = 3
    QTM_Trees%Q(1)%son = [ 0, 0, 2, 3 ]
    QTM_Trees%Q(2)%son = [ 4, 5,  6,  7 ]
    QTM_Trees%Q(3)%son = [ 8, 9, 10, 11 ]
    do hemisphere = 2, 3
      do quadrant = 0, 3 ! Octant number = 4*hemisphere + octant
        octant = 4*hemisphere + quadrant
        QTM_Trees%Q(hemisphere)%son(quadrant) = Add_QTM_Vertex_To_Tree ( octant )
      end do
    end do

    ! Re-size QTM_Trees%Q to the number actually used
    allocate ( Q_Temp(QTM_Trees%n) )
    Q_temp(:) = QTM_Trees%Q(1:QTM_Trees%n)
    call move_alloc ( Q_temp, QTM_Trees%Q )

    if ( QTM_Trees%n_in /= size(ZOT) ) then ! Reallocate ZOT to the correct size
      allocate ( z_temp(QTM_Trees%n_in) )
      z_temp = ZOT(1:QTM_Trees%n_in)
      call move_alloc ( z_temp, ZOT )
    end if

  contains

    integer recursive function Add_QTM_Vertex_To_Tree ( QID ) result( Root )

      use Line_and_Polygon, only: Line_Intersects_Any_Edge
      use Triangle_Interpolate_m, only: Triangle_Interpolate

      integer(qk), intent(in) :: QID

      integer :: F      ! Which facet is being created
      integer :: I
      logical :: In     ! Facet vertex is in a polygon, or polygon vertex is
                        ! in the facet
      real(rg) :: W(3)  ! Barycentric coordinates of polygon vertex in facet,
                        ! to determine whether any of its vertices are within it.

      QTM_trees%n = QTM_trees%n + 1
      root = QTM_trees%n

      if ( QTM_trees%n > size(QTM_Trees%Q) ) then
        ! If QTM_Trees%Q overflows, double its size.
        allocate ( Q_Temp(2*size(QTM_Trees%Q)) )
        q_temp(1:size(QTM_Trees%Q)) = QTM_Trees%Q
        call move_alloc ( Q_Temp, QTM_Trees%Q )
      end if

      QTM_Trees%Q(root)%leaf = .true.

      call QTM_Decode ( QID, S ) ! Get ZOT coordinates of QID into top of S

      QTM_Trees%Q(root)%qid = qid
      QTM_Trees%Q(root)%xn = s%xNode(s%top(octant),octant)  ! xNode number
      QTM_Trees%Q(root)%yn = s%yNode(s%top(octant),octant)  ! yNode number
      QTM_Trees%Q(root)%z = s%z(:,s%top(octant),octant)     ! All the vertices

      ! Determine which vertices, if any, are within the polygon
      in = .false.

      do i = 1, 3
        if ( is_it_in ( s%z(i,s%top(octant),octant) ) ) then
          ! In the polygon; give it a serial number
          in = .true.
          call add_one ( QTM_Trees%Q(root), i )
        end if
      end do

      ! Is any vertex of the polygon within the facet?
      do i = 1, size(QTM_Trees%polygon_ZOT)
        if ( in ) exit
        call triangle_interpolate ( QTM_Trees%Q(root)%z%x, QTM_Trees%Q(root)%z%y, &
          & QTM_Trees%polygon_ZOT(i)%x, QTM_Trees%polygon_ZOT(i)%y, w )
        in = all(w >= 0)
      end do

      ! Does an edge of the facet intersect an edge of the polygon?
      do i = 1, 3
        if ( in ) exit
        in = Line_Intersects_Any_Edge ( &
          & [ QTM_Trees%Q(root)%z(i)%x, QTM_Trees%Q(root)%z(mod(i,3)+1)%x ], &
          & [ QTM_Trees%Q(root)%z(i)%y, QTM_Trees%Q(root)%z(mod(i,3)+1)%y ], &
          & QTM_Trees%polygon_ZOT%x, QTM_Trees%polygon_ZOT%y, &
          & QTM_Trees%ignore_edge ) > 0
      end do

      ! Is the mesh fine enough in this facet?
      if ( s%top(octant) == l ) then
        if ( in ) then
          ! Make sure all vertices of final facets have serial numbers. 
          ! Thereby, a vertex outside the polygon that is a vertex of a final
          ! facet that has any vertex inside the polygon, or any vertex of the
          ! polygon inside it, or any edge of the polygon that intersects one
          ! of its edges, has a serial number.
          do i = 1, 3
            call add_one ( QTM_Trees%Q(root), i )
          end do
        end if
        in = .false. ! So we don't refine further
      end if

      ! If any vertices are within a polygon, or any polygon vertex is
      ! within the facet, refine the facet.
      ! What we really want here is to know whether there's any overlap
      ! between the facet and a polygon.  This can happen even if there are
      ! no vertices of the facet within the polygon, and no vertices of
      ! any polygon within the facet.

      if ( in .and. s%top(octant) <= l) then ! s%top(octant) <= l should
                                             ! always be true
        ! Is the mesh fine enough in this facet?
        do f = 0, 3
          ! The assignment needs to be done in two stages because of a
          ! restriction in the standard: "The evaluation of a function
          ! reference shall neither affect nor be affected by the evaluation
          ! of any other entity within the statement."  Recursive reference to
          ! Add_QTM_Vertex_To_Tree will reallocate QTM_Trees%Q if it runs out
          ! of space.
          QTM_Trees%Q(root)%leaf = .false.
          i = Add_QTM_Vertex_To_Tree ( 4*QID+f )
          QTM_Trees%Q(root)%son(f) = i
        end do
      end if

      ! Don't keep the tree node if it's not interesting
      if ( all(QTM_Trees%Q(root)%son == 0) .and. &
         & all(QTM_Trees%Q(root)%zot_n == 0) ) then
        root = 0
        QTM_Trees%n = QTM_Trees%n - 1
      end if

    end function Add_QTM_Vertex_To_Tree

    subroutine Add_One ( QTM, Node )
      ! Add one vertex to the QTM tree.
      ! Add its coordinates to ZOT if they're not there already.  The ZOT
      ! coordinates are disambiguated as explained below.  These coordinates
      ! should be used to represent the state vector.  The original
      ! coordinates, not the disambiguated ones, are stored in the tree, to be
      ! used for interpolation.  The original coordinates, not the
      ! disambiguated ones, should be used for plotting in the ZOT projection.

      use Hash, only: Found, Inserted, Lookup_And_Insert

      type(QTM_node_t), intent(inout) :: QTM ! Vertex of QTM tree
      integer, intent(in) :: Node  ! Which node 1..3

      ! ZOT coordinates are not unique because the outer edges of the ZOT
      ! projection represent each edge of the southern hemisphere of the
      ! original octahedron twice.  Therefore, if |x| == 1, use |y|, and if
      ! |y| == 1, use |x|.
      ! Observe that 2**L + 1 ZOT coordinates are equally spaced in [-1,1]. 
      ! Therefore, 2**L*(1+[xy]) is in the range [0,2**(L+1)].  Therefore,
      ! 2**(L-1) * ( 1 + x + (2**L + 1) * ( 1 + y ) ) is an (almost)
      ! unique integer identifier in [0, 2**(2*L+2)], which could be rapidly
      ! located in a hash table.  By using |y| if |x| == 1 and |x| if
      ! |y| == 1, 2**(2*L)+2 unique identifiers can be produced.  If it is
      ! desired to represent ZOT coordinates in this compressed form instead
      ! of as two real values, inverting this representation to recover ZOT
      ! coordinates is easier if 2**(L-1) * ( 1 + x + 2**(L+1) * ( 1 + y ) )
      ! is used.  This creates an identifier in the range [2**(L-1),2**(2*L+1)].
      ! The upper bound is the same as the maximum number of facets at level L. 
      ! Since a QID represents a facet, and it fits in an integer, so does this
      ! identifier.

      integer(qk) :: CZ ! z%condense_ZOT()
      integer :: Loc    ! Where was it put or found in the hash table
      integer :: Stat   ! "Inserted" if Z%condense_ZOT() inserted in H
      real(rg) :: X, Y

      if ( QTM%zot_n(node) /= 0 ) return ! Already has a serial number
      cz = qtm%z(node)%condense_ZOT()
      loc = 0 ! Start lookup at z%condense_ZOT()
      do
        call lookup_and_insert ( cz, h, .true., loc, stat )
        if ( stat == inserted ) then
          if ( QTM_Trees%n_in >= size(zot) ) then
            allocate ( z_temp ( 2*size(zot) ) )
            z_temp(1:QTM_Trees%n_in) = zot(1:QTM_Trees%n_in)
            call move_alloc ( z_temp, zot )
          end if
          QTM_Trees%n_in = QTM_Trees%n_in + 1
          ! Disambiguate ZOT coordinates, because the outer edges of the ZOT
          ! projection represent each edge of the southern hemisphere of the
          ! original octahedron twice.
          x = merge(qtm%z(node)%x,abs(qtm%z(node)%x),abs(qtm%z(node)%y)/=1.0_rg)
          y = merge(qtm%z(node)%y,abs(qtm%z(node)%y),abs(qtm%z(node)%x)/=1.0_rg)
          zot(QTM_Trees%n_in) = zot_t(x,y)
          h(2,loc) = QTM_Trees%n_in ! Remember coordinates' serial number
          QTM%zot_n(node) = h(2,loc) ! Coordinates' serial number
          return
        end if
        ! Hash got a match; make sure it's the same, else go around again.
        if ( zot(h(2,loc))%condense_ZOT() == cz ) exit
      end do
      if ( stat == found ) QTM%zot_n(node) = h(2,loc) ! Coordinates' serial number

    end subroutine Add_One

    logical function Is_It_In ( Z )
      ! Returns true if one of the vertices is in the polygon (if
      ! QTM_Trees%in_or_out is true), or outside (if QTM_Trees%in_or_out is false).
      type(ZOT_t), intent(in) :: Z ! Coordinates
      ! Assume Z is not in any polygon
      is_it_in = QTM_Trees%in_or_out .eqv. &
               & pnPoly ( z%x, z%y, QTM_Trees%polygon_ZOT%x, &
                        & QTM_Trees%polygon_ZOT%y ) >= 0
    end function Is_It_In

  end subroutine Generate_QTM

  subroutine Get_QTM_Lats ( QTM_Geo, QTM_Lats )
    ! Get the unique set of QTM latitudes.
    use Geolocation_0, only: GeocLat_t, GeodLat_t, H_Geod, H_t, Lat_t
    class(h_t), intent(in) :: QTM_Geo(:)
    class(lat_t), intent(out), allocatable :: QTM_Lats(:)

    integer :: I, N
    class(lat_t), allocatable :: T(:) ! Temp for QTM_Lats

    select type ( QTM_geo )
    class is ( h_geod )
      allocate ( geodLat_t :: t(size(QTM_geo)) )
    class default
      allocate ( geocLat_t :: t(size(QTM_geo)) )
    end select

    n = 0
    do i = 1, size(QTM_geo)
      if ( .not. any(QTM_geo(i)%lat == t(:n)%d) ) then
        n = n + 1
        t(n)%d = QTM_geo(i)%lat
      end if
    end do

    allocate ( QTM_lats(n), mold=t )
    QTM_lats(:n)%d = t(:n)%d

  end subroutine Get_QTM_Lats

  integer function Cross_Meridian ( G, Z ) result ( T )
    ! Returns the quadrant, i.e., int(longitude/90), of a meridian in the
    ! southern hemisphere for which mod(longitude,90) == 0, that is crossed
    ! by the great circle containing G1 and G2, if it crosses only one,
    ! in which case Z is the ZOT coordinate of the crossing point.  If
    ! -1 there aren't any.  If -3 it crosses two such meridians, and we
    ! can't handle it because we don't know which way around the earth the
    ! edge goes.  If it appears to cross three meridians, i.e., if one end
    ! is in quadrant 3 and the other in quadrant zero, we assume it goes the
    ! short way around and crosses only one.
    use Constants, only: Rad2Deg
    use Geolocation_0, only: Cross, Dot_Product, ECR_t
    use QTM_m, only: Geo_To_ZOT, RG, ZOT_t
    type(h_t), intent(in) :: G(2)
    type(ZOT_t), intent(out) :: Z

    real(rg) :: D            ! dot_product(v1,v2)
    real(rg) :: Lat
    real(rg) :: Lon(2)       ! Longitudes of end points of edge
    integer :: M             ! Longitude of crossed meridian / 90
    type(ECR_t) :: N         ! Normal to great circle containing an edge
    integer :: Q(2)          ! Quadrant number
    type(ECR_t) :: V1, V2    ! Cartesian coordinates of ends of polygon
                             ! edge, then test vectors
    type(ECR_t) :: VX        ! Cartesian coordinates of edge's intersection
                             ! with a meridian having mod(longitude,90) = 0.
    ! X(m) and Y(m) coordinates of normal to meridian M
    real, parameter :: X(0:3) = [  0, -1,  0,  1 ]
    real, parameter :: Y(0:3) = [  1,  0, -1,  0 ]
    t = -1
    if ( g(1)%lat >= 0 .and. g(2)%lat >= 0 ) return ! both ends in the northern
                                                    ! hemisphere
    lon = g%lon%d
    ! Put longitude in [0, 360).
    lon = mod(g%lon%d,360.0_rg)
    where ( lon < 0.0 ) lon = lon + 360.0
    q = int(lon/90.0_rg) ! Each Q in [0,3]
    t = -abs(q(1)-q(2)) - 1
    ! Does the edge have ends in the same quadrant, or two quadrants apart?
    if ( mod(t,2) /= 0 ) return ! t = -1 or -3 => No crossing or two crossings
    v1 = g(1)%ECR()
    v2 = g(2)%ECR()
    n = cross ( v1, v2 ) ! Normal to circle containing G(1) and G(2)
    if ( t == -4 ) then  ! Three crossings only if Q = (3,0) or (0,3),
      m = 0              ! therefore crossing meridian 0.
      if ( g(2)%lon%d > g(1)%lon%d ) n = -n
    else                 ! One crossing, meridian # is max of quadrant #'s
      m = max(q(1),q(2))
      if ( g(1)%lon%d > g(2)%lon%d ) n = -n
    end if
    t = -1               ! No crossings for now.
    ! [ x(m), y(m), 0 ] is normal to the crossed meridian
    vx = ECR_t([ -n%xyz(3)*y(m), n%xyz(3)*x(m), n%xyz(1)*y(m)- n%xyz(2)*x(m) ]) ! N X E
    lat = atan2( vx%xyz(3), abs(vx%xyz(1+mod(m,2))) )
    ! If VX isn't between V1 and V2, i.e., if V1.VX < V1.V2 or V2.VX < V1.V2,
    ! assume the other crossing.
    d = dot_product(v1,v2)
    if ( dot_product(v1,vx) < d .or. dot_product(v2,vx) < d ) lat = -lat
    if ( lat < 0 ) then ! Crosses southern hemisphere meridian M
      z = geo_to_ZOT(h_t(lon_t(90.0_rg*m), lat * rad2deg))
      t = m
    end if
  end function Cross_Meridian

  recursive subroutine Dump_QTM ( QTM, Which, Depth, LatLon, Only, Before, &
                                & Sons, Format )
    ! Dump the frame of the QTM tree specified by Which, preceded by Depth
    ! dots. Then dump its sons with Depth + 1.  You should start this with
    ! Depth == 0.
    use Output_m, only: Output
    use QTM_m, only: Dump_QID, RG

    type(QTM_node_t), intent(in) :: QTM(:)
    integer, intent(in) :: Which
    integer, intent(in) :: Depth
    logical, intent(in), optional :: LatLon ! instead of ZOT, default .false.
    logical, intent(in), optional :: Only   ! Only dump Which frame
    character(*), intent(in), optional :: Before
    logical, intent(in), optional :: Sons   ! Dump sons' indices
    character(*), intent(in), optional :: Format ! For coordinates

    character(6) :: DefaultFmt(2) = [ '(f8.3)', '(f8.4)' ]
    logical :: DoLatLon
    integer :: I
    character(31) :: Fmt
    integer :: NX, NY
    real(rg) :: X, Y

    doLatLon = .false.
    if ( present(latLon) ) doLatLon = latLon
    if ( present(before) ) call output ( before )
    if ( present(format) ) then
      fmt = format
    else if ( doLatLon ) then
      fmt = defaultFmt(1)
    else
      fmt = defaultFmt(2)
    end if
    call output ( which, format='(i5)', after=':' )
    call output ( repeat('.',depth) )
    call dump_qid ( QTM(which)%QID, before=' QID ' )
    nx = QTM(which)%xn; ny = QTM(which)%yn
    call output ( nx, before=' XN ' )
    call output ( ny, before=' YN ' )
    do i = 1, 3
      call get_xy ( QTM(which)%z(i), x, y, doLatLon )
      call output ( x, before=' (', format=fmt )
      call output ( y, before=',', after=')', format=fmt )
    end do
    call output ( QTM(which)%ZOT_n(1), before=' Serial # ' )
    call output ( QTM(which)%ZOT_n(2), before=' ' )
    call output ( QTM(which)%ZOT_n(3), before=' ' )
    if ( present(sons) ) then
      if ( sons ) then
        call output ( ' SONS' )
        do i = 0, 3
          call output ( QTM(which)%son(i), before=' ' )
        end do
      end if
    end if
    call output ( "", advance='yes' )
    if ( present(only) ) then
      if ( only ) return
    end if
    do i = 0, 3
      if ( QTM(which)%son(i) /= 0 ) &
        & call dump_QTM ( QTM, QTM(which)%son(i), depth+1, latLon, sons=sons )
    end do
  contains
    subroutine Get_XY ( Z, X, Y, LatLon )
      type(zot_t), intent(in) :: Z
      real(rg), intent(out) :: X, Y
      logical, intent(in) :: LatLon
      type(h_t) :: Geo
      if ( latLon ) then
        geo = z%zot_to_geo()
        x = geo%lon%d
        y = geo%lat
      else
        x = z%x
        y = z%y
      end if
    end subroutine Get_XY
  end subroutine Dump_QTM

  subroutine Dump_QTM_Tree ( QTM_Trees, Format, LatLon, Sons )
    use Output_m, only: NewLine, Output

    type(QTM_tree_t), intent(in) :: QTM_Trees
    character(*), intent(in), optional :: Format(2) ! to override default,
                                            ! Format(1) is Geo, Format(2) is ZOT
    logical, intent(in), optional :: LatLon ! instead of ZOT, default .false.
    logical, intent(in), optional :: Sons   ! Dump sons' indices

    character(31) :: Fmt(2)
    integer :: I, J

    fmt = [ '(f8.3)', '(f8.4)' ]
    if ( present(format) ) fmt = format

    ! Dump the "inside" point
    call output ( QTM_Trees%in_geo%lon%d, before='Inside Geo (Lon,Lat) (' )
    call output ( QTM_Trees%in_geo%lat, before=',' )
    call output ( QTM_Trees%in%x, before=') Inside ZOT (X,Y) (', format=fmt(2) )
    call output ( QTM_Trees%in%y, before=',', format=fmt(2) )
    call output ( ')', advance='yes' )
    ! Dump the polygons
    call output ( 'Polygon_Geo:', advance='yes' )
    do i = 1, size(QTM_Trees%polygon_geo)
      if ( mod(i,5) == 1 ) call output ( i, format='(i4,"#")' )
      call output ( QTM_Trees%polygon_geo(i)%lon%d, before=' (', format=fmt(1) )
      call output ( QTM_Trees%polygon_geo(i)%lat, before=',', format=fmt(1) )
      call output ( ')' )
      if ( mod(i,5) == 0 .or. i == size(QTM_Trees%polygon_geo) ) call newLine
    end do
    call output ( 'Polygon_ZOT:', advance='yes' )
    do i = 1, size(QTM_Trees%polygon_ZOT)
      if ( mod(i,5) == 1 ) call output ( i, format='(i4,"#")' )
      call output ( QTM_Trees%polygon_ZOT(i)%x, before=' (', format=fmt(2) )
      call output ( QTM_Trees%polygon_ZOT(i)%y, before=',', format=fmt(2) )
      call output ( ')' )
      if ( mod(i,5) == 0 .or. i == size(QTM_Trees%polygon_ZOT) ) call newLine
    end do
    if ( allocated(QTM_Trees%ignore_edge) ) then
      if ( size(QTM_Trees%ignore_edge) > 0 ) then
        call output ( 'Edges on 90*n degree meridians that have antiparallel partners:', &
          & advance='yes' )
        do i = 1, size(QTM_Trees%ignore_edge), 20
          call output ( i, 4, after=": " )
          do j = i, min(i+20,size(QTM_Trees%ignore_edge))
            call output ( QTM_Trees%ignore_edge(j), before=" " )
          end do
          call newLine
        end do
      end if
    end if

    call output ( QTM_Trees%n, before='QTM tree has ' )
    call output ( ' vertices:', advance='yes' )
    if ( QTM_Trees%n > 0 ) &
      & call dump ( QTM_Trees%q, 1, 0, latLon=latLon, sons=sons )
    call output ( QTM_Trees%level, before='Mesh refined to level ' )
    call output ( QTM_Trees%n_in, before=' has ' )
    call output ( ' vertices within or adjacent to the polygon', advance='yes' )

  end subroutine Dump_QTM_Tree

  subroutine Polygon_To_ZOT ( QTM_Trees )
    ! Convert the coordinates that define a polygon on the surface of the
    ! Earth using geographic coordinates (longitude, latitude) to ones in
    ! the ZOT projection.  If an edge of the polygon crosses a meridian in
    ! the southern hemisphere with mod(longitude,90) == 0, break it into
    ! two edges, and insert a third edge that joins the two aliased ZOT
    ! coordinates of the crossing latitude.
    use QTM_m, only: Geo_to_ZOT, Get_Octant, ZOT_t

    type(QTM_Tree_t), intent(inout) :: QTM_Trees

    integer :: IG, IZ ! Subscripts for GEO and ZOT coordinates
    integer :: I, IP, J, JP
    integer :: M      ! Which meridian is crossed, if any; M < 0 means none
    integer :: NG     ! Number of Geo coordinates
    integer :: QN, QP ! Octant of next, previous vertices
    type(ZOT_t), allocatable :: T(:) ! Temp for ZOT coordinates.  Polygon_ZOT
                      ! is allocated with the right size at the end, and copied
                      ! from T.
    integer :: XM     ! What happens to the sign of X as we cross meridian M
                      ! (the change of the sign of  Y is the negative of the
                      ! change of the sign of X).
    type(ZOT_t) :: ZN ! ZOT coordinate of next point; either ZOT coordinate
                      ! of geo coordinate, or ZX.
    type(ZOT_t) :: ZX ! ZOT coordinate where an edge crosses a
                      ! meridian having mod(longitude,90) = 0.

    ng = size(QTM_Trees%polygon_geo)
    if ( ng == 0 ) then
      allocate ( QTM_Trees%polygon_ZOT(1:0) )
      return
    end if

    allocate ( t(3*ng) ) ! Adjust size later
    iz = 1
    t(1) = geo_to_ZOT(QTM_Trees%polygon_geo(1))
    qn = get_octant(t(1)) ! Southern hemisphere octants are 12..15
    do ig = 1, ng-1
      qp = qn
      zn = geo_to_ZOT(QTM_Trees%polygon_geo(ig+1))
      qn = get_octant(zn) ! Southern hemisphere octants are 12..15
      m = -1 ! Assume no crossing
      if ( qp /= qn .and. ( qp >= 12 .or. qn >= 12 ) ) &
        & m = cross_meridian ( QTM_Trees%polygon_geo(ig:ig+1), ZX )
      ! Does edge (ig,ig+1) crosses meridian M from octant QP to QN?
      if ( m >= 0 ) call add_edges
      iz = iz + 1
      t(iz) = zn
    end do
    m = -1 ! Assume no crossing
    qp = qn
    qn = get_octant(t(1))
    if ( qp /= qn .and. ( qp >= 12 .or. qn >= 12 ) ) &
      & m = cross_meridian ( [ QTM_Trees%polygon_geo(ng), QTM_Trees%polygon_geo(1) ], ZX )
    ! Does edge (ng,1) crosses meridian M from octant QP to QN?
    if ( m >= 0 ) call add_edges

    allocate ( QTM_Trees%polygon_ZOT(iz) )
    QTM_Trees%polygon_ZOT(:iz) = t(:iz)

    ! Mark antiparallel edges on southern-hemisphere meridians at mod(lon,90)=0
    ! to be ignored.
    allocate ( QTM_Trees%ignore_edge(iz), source=.false. )
    do i = 1, iz
      ip = 1 + mod(i,iz)
      if ( abs(QTM_Trees%polygon_ZOT(i)%x) == 1.0 .and. &
         & abs(QTM_Trees%polygon_ZOT(ip)%x) == 1.0 .and. &
         & abs(QTM_Trees%polygon_ZOT(i)%y) == abs(QTM_Trees%polygon_ZOT(ip)%y) ) then
        ! Edge I connects two aliased points on a southern-hemisphere meridian
        do j = i, iz
          jp = 1 + mod(j,iz)
          if ( QTM_Trees%polygon_ZOT(i)%x == QTM_Trees%polygon_ZOT(j)%x .and. &
             & QTM_Trees%polygon_ZOT(ip)%x == QTM_Trees%polygon_ZOT(jp)%x .and. &
             & ( QTM_Trees%polygon_ZOT(i)%y - QTM_Trees%polygon_ZOT(jp)%y ) * &
             & ( QTM_Trees%polygon_ZOT(ip)%y - QTM_Trees%polygon_ZOT(j)%y ) <= 0 ) then
            ! Edge J is antiparallel to I and overlaps it.
            QTM_Trees%ignore_edge(i) = .true.
            QTM_Trees%ignore_edge(j) = .true.
          end if
        end do
      else if ( abs(QTM_Trees%polygon_ZOT(i)%y) == 1.0 .and. &
              & abs(QTM_Trees%polygon_ZOT(ip)%y) == 1.0 .and. &
              & abs(QTM_Trees%polygon_ZOT(i)%x) == abs(QTM_Trees%polygon_ZOT(ip)%x) ) then
        ! Edge connects two aliased points on a southern-hemisphere meridian
        do j = i, iz
          jp = 1 + mod(j,iz)
          if ( QTM_Trees%polygon_ZOT(i)%y == QTM_Trees%polygon_ZOT(j)%y .and. &
             & QTM_Trees%polygon_ZOT(ip)%y == QTM_Trees%polygon_ZOT(jp)%y .and. &
             & ( QTM_Trees%polygon_ZOT(i)%x - QTM_Trees%polygon_ZOT(jp)%x ) * &
             & ( QTM_Trees%polygon_ZOT(ip)%x - QTM_Trees%polygon_ZOT(j)%x ) <= 0 ) then
            QTM_Trees%ignore_edge(i) = .true.
            QTM_Trees%ignore_edge(j) = .true.
          end if
        end do
      end if
    end do

  contains

    subroutine Add_Edges
      ! Introduce an edge from T(iz) to ZX, then one from ZX to its alias.
      zx%x = sign(zx%x,t(iz)%x) ! Keep the new edge aimed
      zx%y = sign(zx%y,t(iz)%y) !   in the same direction
      iz = iz + 2
      t(iz-1) = zx
      xm = 2 * mod(m,2) - 1
      t(iz) = ZOT_t(zx%x*xm,-zx%y*xm) ! Reverse the sign of X or Y
    end subroutine Add_Edges

  end subroutine Polygon_To_ZOT

  ! -----     Private Procedures     -----------------------------------

  ! Compute a hash table size that is at least 2**QTM_Depth + 2, and that
  ! has no prime factors less than 11.
  pure integer function Hash_Size ( ) result ( NH )
    use QTM_m, only: QTM_Depth
    nh = 2**(2*QTM_depth)+3
    do
      if ( mod(nh,2) /= 0 .and. mod(nh,3) /= 0 .and. mod(nh,5) /= 0 .and. &
         & mod(nh,7) /= 0 ) exit
      nh = nh + 1
    end do
  end function Hash_Size

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

end module Generate_QTM_m

! $Log$
! Revision 2.5  2016/04/16 02:01:21  vsnyder
! Use ECR_t instead of arrays, remove dependence upon Dump_0
!
! Revision 2.4  2016/03/25 00:27:35  vsnyder
! Lon component now needs to access its %d component.  Make Geo argument
! of Find_Facet_Geo polymorphic.  Simplify call to Triangle_Interpolate.
!
! Revision 2.3  2016/02/26 02:03:37  vsnyder
! Add Get_QTM_Lats subroutine and Leaf component.
!
! Revision 2.2  2016/01/26 19:55:08  vsnyder
! Cannonball polishing
!
! Revision 2.1  2015/12/31 00:01:36  vsnyder
! Initial commit
!
