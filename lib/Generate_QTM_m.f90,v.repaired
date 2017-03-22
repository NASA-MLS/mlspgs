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

  use Geolocation_0, only: ECR_t, Lat_t, Lon_t, H_Geod, H_t, H_V_Geoc, RG
  use QTM_m, only: QK, ZOT_t

  implicit NONE
  private

  public :: Cross_Meridian, Generate_QTM, Get_QTM_Lats, Polygon_To_ZOT, &
          & QTM_node_t, QTM_Tree_t, ZOT_t

  ! One facet of a QTM, at the finest refinement if %depth == QTM_Tree_t%Level.
  ! All(%son==0) means a facet is not refined, but it might not be at the
  ! finest refinement because it's entirely outwith the polygon.
  type :: QTM_Node_t
    integer :: Depth = 0       ! A leaf node if Depth == QTM_Tree_t%Level
    integer(qk) :: QID = 0     ! QTM ID
    integer :: Son(0:3) = 0    ! Sub facet subscripts in the tree
    integer :: XN = 0, YN = 0  ! Which son is xNode, yNode?  Central node is 0,
                               ! pole node is 6 - ( xn + yn )
    type(h_t) :: Geo(3)        ! (lon,lat) coordinates of vertices, degrees
    type(ZOT_t) :: Z(3)        ! ZOT coordinates of vertices
    integer :: Ser(3) = 0      ! Serial numbers of coordinates of vertices in 
                               ! QTM_Tree_t%ZOT_In and QTM_Tree_t%Geo_In if
                               ! nonzero, else the vertex is not either in or
                               ! adjacent to the polygon.  A leaf node always
                               ! has three nonzero vertex serial numbers, and
                               ! four zero son indices.
  end type QTM_Node_t

  type :: QTM_Tree_t
    integer :: Level = 7       ! Level of QTM refinement, 20000 / 2^level km
                               ! or 180 / 2^level degrees latitude spacing. 
                               ! Not greater than QTM_Depth from QTM_m
    integer(qk) :: N           ! Number of used elements in Q
    logical :: Geodetic = .true. ! "Latitudes are geodetic" or not
    type(ZOT_t) :: In = ZOT_t(-999,-999) ! A point defined to be inside the
                               ! polygon. This is needed because the concept of
                               ! "inside a polygon" is ambiguous on a sphere.
    type(h_t) :: In_Geo = h_t(lon_t(-999),-999) ! QTM_Tree_t%In as longitude
                               ! and latitude, degrees, either geodetic or
                               ! geocentric.
    logical :: In_or_Out       ! True iff PnPoly reports QTM_Tree_t%In is
                               ! inside the polygon, from the point of view of
                               ! the ZOT projection.  A point is considered to
                               ! be inside if PnPoly returns the same value.
    integer(qk) :: N_Facets = 0 ! Number of finest-refinement facets of the QTM
    integer(qk) :: N_In = 0    ! Number of vertices inside the polygon or
                               ! immediately adjacent to it.
    type(h_t), allocatable :: Polygon_Geo(:) ! (lon,lat) degrees, either
                               ! geodetic or geocentric.
    type(h_t) :: Polygon_Geo_Centroid        ! (lon,lat) degrees, either
                               ! geodetic or geocentric.
    type(ZOT_t), allocatable :: Polygon_ZOT(:)
    type(ZOT_t):: Polygon_ZOT_Centroid
    logical, allocatable :: Ignore_Edge(:)   ! Ignore edge (i,1+mod(i,n))
                               ! because it joins two aliased points on a
                               ! southern-hemisphere meridian at mod(lon,90)=0
                               ! and there is an edge antiparallel to it and
                               ! colinear with it.
    type(QTM_node_t), allocatable :: Q(:)    ! The nodes of the search tree
    integer, allocatable :: Adjacent_in(:,:) ! Search tree node indices of
                               ! facets adjacent to each vertex of Geo_In or
                               ! ZOT_In.  First dimension is 6.
    integer, allocatable :: N_Adjacent_in(:) ! Effective first extents of
                               ! Adjacent_in.
    class(h_t), allocatable :: Geo_In(:)     ! H_T coordinates corresponding to
                               ! ZOT_In, (lon,lat) degrees, either geodetic or
                               ! geocentric.
    type(ZOT_t), allocatable :: ZOT_In(:)    ! ZOT coordinates of vertices of
                               ! QTM that are inside or adjacent to the
                               ! polygon. Indexed by Q(.)%Ser(.).
    class(lat_t), allocatable :: QTM_Lats(:) ! Unique latitudes within QTM.
                               ! They're geocentric or geodetic, depending upon
                               ! the dynamic type of Geo_In.
    integer, allocatable :: Path_Vertices(:) ! Mapping from QTM vertex indices
                               ! to indices adjacent to a line of sight.
                               ! Filled separately for each path.  Inverse of
                               ! Facets_and_Vertices_t%vertices.
  contains
    procedure :: Find_Facet_ECR
    procedure :: Find_Facet_Geo
    procedure :: Find_Facet_QID
    procedure :: Find_Facet_ZOT
    generic :: Find_Facet => Find_Facet_ECR
    generic :: Find_Facet => Find_Facet_Geo
    generic :: Find_Facet => Find_Facet_QID
    generic :: Find_Facet => Find_Facet_ZOT
  end type QTM_Tree_t

  ! Initial size to allocate for QTM_Tree_t%Q if it's not allocated
  integer, public :: Default_Initial_Size = 100

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  integer function Find_Facet_ECR ( QTM_Trees, ECR, Stack ) result ( F )
    ! Return the subscript in QTM_Trees of the facet containing the position
    ! where the line from the Earth's center to ECR pierces the surface.
    use QTM_m, only: Stack_t
    class(QTM_Tree_t), intent(in) :: QTM_Trees
    type(ECR_t), intent(in) :: ECR
    type(stack_t), intent(inout), optional :: Stack ! to make subsequent search faster
    type(h_v_geoc) :: Geoc
    type(h_geod) :: Geod
    if ( QTM_Trees%geodetic ) then
      geod = ecr%geod_surf()
      f = qtm_trees%find_facet ( geod, stack )
    else
      geoc = ecr%geoc()
      f = qtm_trees%find_facet ( geoc, stack )
    end if
  end function Find_Facet_ECR

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

    ! Start at the root (f == 1).  The highest-order two nonzero bits of the
    ! QID are either 2 (indexing the Northern hemisphere frame) or 3 (indexing
    ! the Southern hemisphere frame).  Thereafter, they are subdivision
    ! indices, according to how basis numbers are reassigned as the QTM is
    ! refined.  QTM_Trees%Q(1)%son(2) == 2 and QTM_Trees%Q(1)%son(3) == 3.

    f = -1
    if ( qid < 8 ) return ! Bogus QID
    h = high_bit_index(QID)
    ! The number of iterations is bounded by the refinement level.
    do while ( h > 0 )
      h = h - 2
      qd = iand(int(shiftr(qid,h)),3) ! Peel off the next two bits of the QID
      next = QTM_Trees%Q(f)%son(qd)   ! Use those bits to select the next frame
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

  subroutine Generate_QTM ( QTM_Trees )

    ! Generate a QTM refined to QTM_Trees%Level (in 1..QTM_Depth) within a
    ! Polygon defined by QTM_Trees%Polygon_Geo or QTM_Trees%Polygon_ZOT.  A
    ! point inside the polygon, specified by QTM_Trees%In_Geo or QTM_Trees%In,
    ! is needed because the concept of "inside a polygon" is ambiguous on a
    ! sphere.

    ! Whether a QTM vertex is within a polygon is determined within the ZOT
    ! projection.  Depending upon how edges of facets are defined on the
    ! Earth, and the resolution of the polygon, a point determined to be
    ! within (outwith) a polygon in the ZOT projection might be outwith
    ! (within) the polygon on the surface of the Earth if it's near an edge.

    ! The ZOT coordinates of the vertices that are within or adjacent to the
    ! polygon are in QTM_Trees%ZOT_In.  The (lon,lat) coordinates are in
    ! QTM_Trees%Geo_In.  For each facet, they are also in QTM_Trees%Q(F)%Z
    ! and QTM_Trees%Q(F)%Geo, which could be gotten from
    ! QTM_Trees%ZOT_In(QTM_Trees%Q(F)%Ser) and
    ! QTM_Trees%Geo_In(QTM_Trees%Q(F)%Ser), provided
    ! QTM_Trees%Q(F)%Ser /= 0, which is true for facets at the finest
    ! refinement but not necessarily at coarser refinement.  If a vertex V of
    ! a coarse-refinement facet C is not also a vertex of a finest-refinement
    ! facet, QTM_Trees%Q(C)%Ser(v) == 0.  The search tree indices of facets
    ! adjacent to QTM_Trees%Geo_In(V) or QTM_Trees%ZOT_In(V) are in
    ! QTM_Trees%Adjacent_in(:QTM_Trees%N_Adjacent_in(V),V).

    use Geolocation_0, only: H_Geoc, H_Geod
    use PnPoly_m, only: PnPoly
    use QTM_m, only: Geo_to_ZOT, QTM_Decode, QTM_Depth, Stack_t, ZOT_t

    type(QTM_Tree_t), intent(inout), target :: QTM_Trees

    integer, allocatable :: H(:,:) ! Assume QK = kind(0) for now, until
                          ! we have a generic Hash module
    type(h_t), pointer :: Geo(:) ! Nonpolymorphic handle for QTM_Trees%geo_in
    integer :: Hemisphere ! 2 = north, 3 = south
    integer, allocatable :: I1_Temp(:)   ! 1-dimensional integer temp
    integer, allocatable :: I2_Temp(:,:) ! 2-dimensional integer temp
    integer :: L          ! min(QTM_Trees%Level,QTM_Depth)
    integer :: Octant     ! Octant being refined.
    integer :: Quadrant   ! Quadrant in a hemisphere, 0..3.
    type(QTM_node_t), allocatable :: Q_Temp(:) ! For reallocating Q
    type(stack_t) :: S    ! Put here so we don't need a new one in every
                          ! instance of Add_QTM_Vertex_To_Tree.  Speeds up
                          ! QTM_Decode in Add_QTM_Vertex_To_Tree if it's
                          ! not new for every instance.
    integer :: Son        ! A temp
    real(rg) :: X, Y      ! Disambiguated ZOT coordinates
    type(ZOT_t), allocatable :: Z_Temp(:) ! For reallocating ZOT

    QTM_Trees%n = 0 ! In case setup is wrong
    if ( .not. allocated(QTM_Trees%polygon_ZOT) ) then
      if ( .not. allocated(QTM_Trees%polygon_geo) ) return ! without doing anything
      call polygon_to_ZOT ( QTM_Trees )
    end if

    ! The polygon centroid is later useful for finding a point on a line of
    ! sight that is nearest to it, and from that finding the facet below that
    ! point, and from that finding all facets crossed by that line.
    QTM_Trees%polygon_geo_centroid = &
      & h_t(lon_t(sum(QTM_Trees%polygon_geo%lon%d)/size(QTM_Trees%polygon_geo)), &
         &  sum(QTM_Trees%polygon_geo%lat)/size(QTM_Trees%polygon_geo) )
    QTM_Trees%polygon_ZOT_centroid = &
      & ZOT_t(sum(QTM_Trees%polygon_ZOT%x)/size(QTM_Trees%polygon_ZOT), &
           &  sum(QTM_Trees%polygon_ZOT%y)/size(QTM_Trees%polygon_ZOT) )

    ! Now that we know the centroid of the polygon, we can calculate the
    ! fraction of the Earth's surface occupied by a spherical cap whose angular
    ! extent is the maximum angle between any polygon vertex and the centroid,
    ! and use that fraction to calculate the hash size.

    allocate ( h(2,hash_size(QTM_Trees)) )

    if ( QTM_Trees%in%x < -10 .or. QTM_Trees%in%y < -10 ) then
      if ( QTM_Trees%in_geo%lon%d < -500 .or. QTM_Trees%in_geo%lat < -500 ) return
      QTM_Trees%in = geo_to_ZOT(QTM_Trees%in_geo)
    end if
    l = min(QTM_Trees%level,QTM_Depth)
    QTM_Trees%level = l
    h = 0

    x = merge(QTM_Trees%in%x,abs(QTM_Trees%in%x),abs(QTM_Trees%in%y)/=1.0_rg)
    y = merge(QTM_Trees%in%y,abs(QTM_Trees%in%y),abs(QTM_Trees%in%x)/=1.0_rg)

    QTM_Trees%in_or_out = pnPoly ( x, y, QTM_Trees%polygon_ZOT%x, &
                                 & QTM_Trees%polygon_ZOT%y ) >= 0

    if ( .not. allocated(QTM_Trees%Q) ) allocate ( QTM_Trees%Q(default_initial_size) )
    allocate ( QTM_Trees%ZOT_In(default_initial_size) )
    allocate ( QTM_Trees%adjacent_in(6,default_initial_size) )
    allocate ( QTM_Trees%n_adjacent_in(default_initial_size), source=0 )
    ! Create the root and top two hemispheres in the tree.
    QTM_Trees%n = 3
    QTM_Trees%Q(1)%son = [ 0, 0,  2,  3 ]; QTM_Trees%Q(1)%depth=0 ! Root
    QTM_Trees%Q(1)%geo = h_t(lon_t(-90),90); QTM_Trees%Q(1)%Z = ZOT_t(-1,1)
    QTM_Trees%Q(2)%depth=0 ! Northern hemisphere
    QTM_Trees%Q(2)%geo = h_t(lon_t(-90),90); QTM_Trees%Q(2)%Z = ZOT_t(-1,1)
    QTM_Trees%Q(3)%depth=0 ! Southern hemisphere
    QTM_Trees%Q(3)%geo = h_t(lon_t(-90),-90); QTM_Trees%Q(3)%Z = ZOT_t(-1,1)
    do hemisphere = 2, 3
      do quadrant = 0, 3 ! Octant number = 4*hemisphere + octant
        octant = 4*hemisphere + quadrant
        ! The assignment needs to be done in two stages because of a
        ! restriction in the standard: "The evaluation of a function
        ! reference shall neither affect nor be affected by the evaluation
        ! of any other entity within the statement."  Recursive reference to
        ! Add_QTM_Vertex_To_Tree will reallocate QTM_Trees%Q if it runs out
        ! of space.
        son = Add_QTM_Vertex_To_Tree ( octant, 0 )
        QTM_Trees%Q(hemisphere)%son(quadrant) = son
      end do
    end do

    if ( QTM_Trees%n /= size(QTM_Trees%Q) ) then
      ! Reallocate QTM_Trees%Q to the number actually used
      allocate ( Q_Temp(QTM_Trees%n) )
      Q_temp(:) = QTM_Trees%Q(1:QTM_Trees%n)
      call move_alloc ( Q_temp, QTM_Trees%Q )
    end if

    ! Allocate Path_Vertices so Get_Lines_of_Sight has a place to
    ! store the mapping from QTM vertices to vertices adjacent to
    ! line of sight path.
    allocate ( QTM_Trees%path_vertices(QTM_Trees%n_in) )

    if ( QTM_Trees%n_in /= size(QTM_Trees%ZOT_In) ) then
      ! Reallocate ZOT etc. to the correct size
      allocate ( z_temp(QTM_Trees%n_in) )
      z_temp = QTM_Trees%ZOT_In(1:QTM_Trees%n_in)
      call move_alloc ( z_temp, QTM_Trees%ZOT_In )
      allocate ( i1_temp(QTM_Trees%n_in) )
      i1_temp = QTM_Trees%n_adjacent_in(1:QTM_Trees%n_in)
      call move_alloc ( i1_temp, QTM_Trees%n_adjacent_in )
      allocate ( i2_temp(6,QTM_Trees%n_in) )
      i2_temp = QTM_Trees%adjacent_in(1:6,1:QTM_Trees%n_in)
      call move_alloc ( i2_temp, QTM_Trees%adjacent_in )
    end if

    ! Compute (lon,lat) coordinates corresponding to ZOT coordinates
    if ( QTM_Trees%Geodetic ) then
      allocate ( h_geod :: QTM_Trees%geo_in(1:QTM_Trees%n_in) )
    else
      allocate ( h_geoc :: QTM_Trees%geo_in(1:QTM_Trees%n_in) )
    end if

    geo => QTM_Trees%geo_in ! Get a nonpolymorphic handle
    geo = QTM_Trees%ZOT_In%zot_to_geo()

    ! Get the unique latitudes
    call Get_QTM_Lats ( QTM_Trees%geo_in, QTM_Trees%QTM_Lats )

  contains

    integer recursive function Add_QTM_Vertex_To_Tree ( QID, Depth ) result( Root )

      use Line_and_Polygon, only: Line_Intersects_Any_Edge
      use Triangle_Interpolate_m, only: Triangle_Interpolate

      integer(qk), intent(in) :: QID
      integer, intent(in) :: Depth

      real(rg) :: B(3)  ! Barycentric coordinates of polygon vertex in facet,
                        ! to determine whether any of its vertices are within it.
      integer :: F      ! Which facet is being created
      integer :: I
      logical :: In     ! Facet vertex is in a polygon, or polygon vertex is
                        ! in the facet, or an edge of the polygon crosses an
                        ! edge of the facet.
      integer :: V      ! Serial number of a vertex.

      QTM_trees%n = QTM_trees%n + 1
      root = QTM_trees%n

      if ( QTM_trees%n > size(QTM_Trees%Q) ) then
        ! If QTM_Trees%Q overflows, double its size.
        allocate ( Q_Temp(2*size(QTM_Trees%Q)) )
        q_temp(1:size(QTM_Trees%Q)) = QTM_Trees%Q
        call move_alloc ( Q_Temp, QTM_Trees%Q )
      end if

      call QTM_Decode ( QID, S ) ! Get ZOT coordinates of QID into top of S

      QTM_Trees%Q(root)%depth = depth + 1
      QTM_Trees%Q(root)%qid = qid
      QTM_Trees%Q(root)%xn = s%xNode(s%top(octant),octant)  ! xNode number
      QTM_Trees%Q(root)%yn = s%yNode(s%top(octant),octant)  ! yNode number
      QTM_Trees%Q(root)%z = s%z(:,s%top(octant),octant) ! All the vertices ZOT
      QTM_Trees%Q(root)%geo = QTM_Trees%Q(root)%z%ZOT_to_geo() ! Vertices Geo

      in = .false.

      ! Determine which vertices, if any, are within the polygon
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
          & QTM_Trees%polygon_ZOT(i)%x, QTM_Trees%polygon_ZOT(i)%y, b )
        in = all(b >= 0)
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

      ! Is this facet of the mesh too coarse?
      if ( QTM_Trees%Q(root)%depth >= l ) then
        ! Should we refine this facet?
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

      ! If any vertices are within the polygon, or any polygon vertex is
      ! within the facet, or any edge of the polygon crosses an edge of
      ! the facet (i.e., there's any overlap between the polygon and
      ! the facet), and the facet is not of sufficiently fine resolution,
      ! refine the facet.

      if ( in ) then
        ! The mesh is not fine enough in this facet and it needs to be refined.
        do f = 0, 3
          ! The assignment needs to be done in two stages because of a
          ! restriction in the standard: "The evaluation of a function
          ! reference shall neither affect nor be affected by the evaluation
          ! of any other entity within the statement."  Recursive reference to
          ! Add_QTM_Vertex_To_Tree will reallocate QTM_Trees%Q if it runs out
          ! of space.  A processor might decide where to store the result before
          ! it invokes the function, which would be the wrong place if
          ! QTM_Trees%Q gets reallocated.
          QTM_Trees%Q(root)%depth = depth + 1
          i = Add_QTM_Vertex_To_Tree ( 4*QID+f, depth + 1 )
          QTM_Trees%Q(root)%son(f) = i
        end do
      end if

      ! Don't keep the tree node if it's not interesting, that is, the facet
      ! is too coarse and not refined.
      if ( all(QTM_Trees%Q(root)%son == 0) .and. &
         & all(QTM_Trees%Q(root)%ser == 0) ) then
        root = 0
        QTM_Trees%n = QTM_Trees%n - 1
      else if ( QTM_Trees%Q(root)%depth == QTM_Trees%level ) then
        QTM_Trees%n_facets = QTM_Trees%n_facets + 1
        ! Add facet to list of facets adjacent to its vertices
        do i = 1, 3
          v = QTM_Trees%Q(root)%ser(i)
          f = QTM_Trees%n_adjacent_in(v)
          if ( any(QTM_Trees%adjacent_in(:f,v) == root) ) cycle
          f = f + 1
          QTM_Trees%n_adjacent_in(v) = f
          QTM_Trees%adjacent_in(f,v) = root
        end do
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

      use Hash, only: Full, Inserted, Lookup_And_Insert
      use MLSMessageModule, only: MLSMSG_Error
      use MoreMessage, only: MLSMessage

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
      integer :: N      ! Size(QTM_Trees%ZOT_in) before adding the vertex
      integer :: Stat   ! "Inserted" if Z%condense_ZOT() inserted in H
      real(rg) :: X, Y

      if ( QTM%ser(node) /= 0 ) return ! Already has a serial number
      cz = qtm%z(node)%condense_ZOT()
      loc = 0 ! Start lookup at mod(z%condense_ZOT(),size(h,2))
      do
        call lookup_and_insert ( cz, h, .true., loc, stat )
        if ( stat == full ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & 'Hash table, allocated with %d elements, is full.', size(h,2) )
        end if
        if ( stat == inserted ) then
          QTM_Trees%n_in = QTM_Trees%n_in + 1
          n = size(QTM_Trees%ZOT_in)
          if ( QTM_Trees%n_in >= n ) then
            ! Make ZOT_in etc twice as big as before
            n = 2 * n
            allocate ( z_temp ( n ) )
            z_temp(1:QTM_Trees%n_in) = QTM_Trees%ZOT_in(1:QTM_Trees%n_in)
            call move_alloc ( z_temp, QTM_Trees%ZOT_in )
            allocate ( i1_temp ( n ) )
            i1_temp(1:QTM_Trees%n_in) = &
              & QTM_Trees%n_adjacent_in(1:QTM_Trees%n_in)
            call move_alloc ( i1_temp, QTM_Trees%n_adjacent_in )
            allocate ( i2_temp(6,n) )
            i2_temp(1:6,1:QTM_Trees%n_in) = &
              & QTM_Trees%adjacent_in(1:6,1:QTM_Trees%n_in)
            call move_alloc ( i2_temp, QTM_Trees%adjacent_in )
          end if
          ! Disambiguate ZOT coordinates, because the outer edges of the ZOT
          ! projection represent each edge of the southern hemisphere of the
          ! original octahedron twice.
          x = merge(qtm%z(node)%x,abs(qtm%z(node)%x),abs(qtm%z(node)%y)/=1.0_rg)
          y = merge(qtm%z(node)%y,abs(qtm%z(node)%y),abs(qtm%z(node)%x)/=1.0_rg)
          QTM_Trees%ZOT_in(QTM_Trees%n_in) = zot_t(x,y)
          QTM_Trees%n_adjacent_in(QTM_Trees%n_in) = 0
          h(2,loc) = QTM_Trees%n_in ! Remember coordinates' serial number
          exit
        end if
        ! Hash got a match; make sure it's the same, else go around again.
        if ( QTM_Trees%ZOT_in(h(2,loc))%condense_ZOT() == cz ) exit
      end do
      QTM%ser(node) = h(2,loc) ! Coordinates' serial number

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

  ! Compute a hash table size that is big enough to represent all the
  ! vertices that might be in the polygon, and that has no prime factors
  ! less than 11.

!   pure &
  integer function Hash_Size ( QTM_Trees ) result ( NH )

    type(QTM_Tree_t), intent(in) :: QTM_Trees
    type(ECR_t) :: Centroid ! Polygon centroid
    integer :: I, L
    real(rg) :: MinCos      ! minimum of cos(between Centroid and Point)
    type(ECR_t) :: Point    ! A polygon vertex

    ! The area of a spherical cap of angular extent A is 2 pi R^2 ( 1 - cos(A)).
    ! A spherical cap that encloses the polygon has an angular extent of the
    ! maximum angle between any vertex and the centroid.  The area of the Earth
    ! (as a sphere) is 4 pi R^2.  So the fraction of the Earth covered by a
    ! spherical cap big enough to enclose the polygon is 0.5 * ( 1 - cos(A) ). 
    ! We also need to allow for vertices adjacent to the polygon by adding the
    ! equivalent of one more parallel of latitude, which is done by multiplying
    ! the area by 1 + 1.0 / 2**level

    minCos = 2
    centroid = QTM_Trees%polygon_geo_centroid%surf_ECR(norm=.true.)
    do i = 1, size(QTM_Trees%polygon_geo)
      point = QTM_Trees%polygon_geo(i)%surf_ECR(norm=.true.)
      minCos = min(minCos, point .dot. centroid)
    end do

    l = 2 ** ( QTM_Trees%level - 1 )

    nh = ceiling(5 * l * l * (1 + 2.0/l) * (1 - minCos)) + 3

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
! Revision 2.26  2017/03/22 23:19:24  vsnyder
! The hash table size calculation is too small for very small polygons.
! Make it five times bigger.
!
! Revision 2.25  2016/12/09 00:44:26  vsnyder
! Extra space in hash table for boundary vertices
!
! Revision 2.24  2016/12/08 02:41:06  vsnyder
! Fill Geo and Z components of first three tree vertices, just so they don't
! have undefined values, in case somebody wants to print them.
!
! Revision 2.23  2016/12/08 02:22:35  vsnyder
! Corrected three mistakes concerning re-allocating components as QTM grew:
! 1. Checked size before increment requirement.  2-3. After computing new
! size and allocating i1_temp and i2_temp, assigned to them without (1:n_in),
! so they were automatically reallocated to the old size.  Calculate the hash
! table size using the fraction of the Earth's surface occupied by a
! spherical cap with a radial extent equal to the maximum angle between any
! polygon vertex and its centroid, instead of the whole Earth.
!
! Revision 2.22  2016/11/12 01:31:28  vsnyder
! Add Path_Vertices, to store mapping from QTM to vertices near the path
!
! Revision 2.21  2016/11/03 20:52:51  vsnyder
! Return negative facet number for point outside QTM
!
! Revision 2.20  2016/11/03 20:40:39  vsnyder
! Add Find_Facet_ECR
!
! Revision 2.19  2016/11/02 22:58:58  vsnyder
! Compute Polygon centroid in Geo and ZOT coordinates.  Compute facets
! adjacent to each vertex.
!
! Revision 2.18  2016/10/18 00:43:15  vsnyder
! Remove declaration for unused variable
!
! Revision 2.17  2016/10/05 23:28:22  vsnyder
! Replace ZOT_n component name with Ser because it's a serial number for
! more than just the ZOT coordinates.
!
! Revision 2.16  2016/09/23 01:56:16  vsnyder
! Add Geodetic and QTM_Lats components.  Fill QTM_Lats.  Make geolocation
! components polymorphic so they encode whether they're geodetic or not.
!
! Revision 2.15  2016/09/15 22:49:42  vsnyder
! Resize QTM_Trees%Q only if it's the wrong size, add and revise a ton of
! comments.
!
! Revision 2.14  2016/09/14 20:08:32  vsnyder
! Correct depth component calculation
!
! Revision 2.13  2016/09/13 20:08:12  vsnyder
! Replace Leaf component with Depth
!
! Revision 2.12  2016/09/10 01:57:46  vsnyder
! Use LatLon argument to print 'Polygon has n vertices' header
!
! Revision 2.11  2016/09/10 01:50:41  vsnyder
! Use MyLatLon instead of optional LatLon in Dump_QTM_Tree
!
! Revision 2.10  2016/09/10 01:48:07  vsnyder
! Indent dump headers so they line up with other dumps
!
! Revision 2.9  2016/08/24 01:35:51  vsnyder
! Add (lon,lat) coordinates of facet vertices in QTM_Node_t
!
! Revision 2.8  2016/08/23 01:25:58  vsnyder
! Allocate QTM_Trees%geo_in explicitly because ifort doesn't do it
! automagigically without a command-line option.
!
! Revision 2.7  2016/08/23 00:40:48  vsnyder
! Put coordinates of vertices within or adjacent to the polygon in components
! within the QTM_Tree_t structure.
!
! Revision 2.6  2016/07/30 00:51:35  vsnyder
! Change default QTM level to 7 (156 km)
!
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
