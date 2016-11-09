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
module Metrics_3D_m
!=============================================================================

  use Geolocation_0, only: ECR_t, RG

  implicit NONE

  public

  real(rg), public :: Height_Tol = 100.0 ! Meters.  How close must the height
                     ! at an extrapolation along the line-of-sight be to the
                     ! height where the line-of-sight meets the region of
                     ! interest?

  interface Metrics_3D
    module procedure Metrics_3D_QTM   ! The entire path
    module procedure Metrics_3D_QTM_1 ! Half of the path, before or after
                                      ! tangent or reflection
  end interface Metrics_3D

  ! Which surface penetrations to put into the path:
  enum, bind(C)
    enumerator :: All        ! Both horizontal and vertical surfaces
    enumerator :: Horizontal ! Horizontal (zeta or height) surfaces
    enumerator :: Vertical   ! Vertical surfaces (planes and latitude cones)
  end enum

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Metrics_3D_QTM ( Path, QTM_Tree, H, S, Tangent_Index, Pad, Facets, &
                            & Which )
    ! Given a line defined by a point in ECR, and a vector in ECR parallel
    ! to that line, compute all intersections of that line with a face of
    ! the grid whose horizontal grid is QTM_Tree.  It is assumed the vertical
    ! grid given by H is stacked, but it needn't be coherent.

    ! Path%Lines(2,1) and Path%Lines(2,2) are made unit vectors here.

    ! The given line is Path%Lines(1,1) + s * Path%Lines(2,1).  A line defined
    ! by Path%Lines(1,2) + s * Path%Lines(2,2) is produced, which is the
    ! continuation of Path%Lines(:,1) after the tangent point if
    ! Path%Lines(:,1) does not intersect the Earth reference ellipsoid, or the
    ! reflection of Path%Lines(:,1) if it does intersect the Earth reference
    ! ellipsoid.  Path%Lines(2,2) is an unit vector.

    ! The tangent point is Path%Lines(2,1).

    ! Values of S(1:Tangent_Index) are in order such that the first one is
    ! farthest from the tangent point in the direction toward Path%Lines(1,1)
    ! and the last one is the tangent point.  Values of S(Tangent_Index:) are
    ! in order such that the first one is the tangent point (again) and the
    ! last one is farthest from Path%Lines(1,1).

    use Generate_QTM_m, only: QTM_Tree_t
    use Path_Representation_m, only: Path_t
    use QTM_Interpolation_Weights_3D_m, only: S_QTM_t

    type(path_t), intent(inout) :: Path ! Path%Lines(1,1) + s * Path%Lines(2,1)
                                   ! defines the initial line, i.e., before
                                   ! the tangent point.
                                   ! Path%Lines(2,1) + s * Path%Lines(2,2)
                                   ! defines the line after the tangent point.
    type(QTM_Tree_t), intent(inout), target :: QTM_Tree
    real(rg), intent(in), contiguous :: H(:,:) ! Heights X size(QTM_Tree%Geo_in)
                                   ! The second subscript is the QTM vertex
                                   ! serial number, also used as a subscript
                                   ! for QTM_Tree%Geo_in.
    type(S_QTM_t), intent(out), allocatable :: S(:) ! S-values, and interpolation
                                   ! coefficients
    integer, intent(out) :: Tangent_Index ! Index in S of tangent or
                                   ! intersection.
                                   ! S(Tangent_Index) == S(Tangent_Index+1)
    integer, intent(in) :: Facets(:) ! Facets under Path
    integer, intent(in) :: Pad     ! Amount of padding to introduce in S between
                                   ! before-the-tangent and after-the-tangent.
    integer, intent(in), optional :: Which ! Which intersections to detect,
                                   ! default all

    type :: Temp
      type(S_QTM_t), allocatable :: S(:)
    end type Temp

    integer :: I
    type(Temp) :: S_Int(2)         ! Intersections from Metrics_3D_QTM_1

    call path%get_path_ready

    ! Do the real work, one half at a time -- before and after the
    ! tangent or Earth-surface intersection.
    do i = 1, 2
      call Metrics_3D_QTM_1 ( path%Lines(:,i), path%sMin(i), path%sMax(i), &
        & QTM_tree, h, s_int(i)%s, facets, which )
    end do
    ! Concatenate the two halves of the path, with Pad copies of the tangent
    ! between.  s(tangent_index:tangent_index+pad+1) should all be the same.
    ! s(tangent_index+1:tangent_index+pad) probably aren't used anywhere, but
    ! just to be safe we make sure they're defined.
    tangent_index = size(s_int(1)%s,1)
    allocate ( s, source=[ s_int(1)%s, spread(s_int(2)%s(1),1,pad), &
                         & s_int(2)%s ] )

  end subroutine Metrics_3D_QTM

  subroutine Metrics_3D_QTM_1 ( Line, SMin, SMax, QTM_Tree, H, Intersections, &
                              & Facets, Which )

    ! Given a line defined by a point in ECR, and a vector in ECR parallel
    ! to that line, compute all intersections of that line with a face of
    ! the grid defined by QTM_Tree and H.  The grid is assumed to be stacked
    ! but is not required to be coherent.  The line is Line(1) + s * Line(2).
    ! Only intersections with sMin <= s <= sMax are reported.

    ! The algorithm proceeds in four stages:
    ! 1. Compute all intersections with horizontal boundary surfaces, which
    !    are approximated by a spherical cap passing through the three points
    !    bounding the surface, with a radius of curvature given by the average
    !    curvature of the Earth at the circumcenter of the QTM facet.  We
    !    assume that the geodetic height is so small that the curvature at
    !    the boundary surface is essentially the same as at the surface.  The
    !    center is on the gradient to the surface at the circumcenter.
    ! 2. Compute all intersections with latitude cones.
    ! 3. Compute all intersections with vertical boundary faces that are not
    !    latitude cones.  These are planes, but not necessarily passing through
    !    the center of the Earth.
    ! 4. Compute all intersections with surfaces of constant height above the
    !    Earth reference ellipsoid that are not within the polygon in which the
    !    QTM is constructed.  The heights are interpolated on the face of the
    !    QTM facet intersected by Line that is nearest to the edge of the
    !    polygon in which the QTM was constructed, at the latitude and longitude
    !    of the intersection.
    ! At the end, the intersections are sorted according to S.

    !!!!!=======================================================!!!!!
    !!!!!                                                       !!!!!
    !!!!!      THIS DOES NOT WORK FOR CONCAVE POLYGONS!         !!!!!
    !!!!!                                                       !!!!!
    !!!!!=======================================================!!!!!

    use Center_of_Sphere_m, only: Center_of_Sphere, Circumcenter
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_Geod, H_V_Geoc, H_V_Geod, H_V_t
    use Line_And_Cone_m, only: Line_And_Cone
    use Line_And_Ellipsoid_m, only: Line_And_Sphere
    use Line_And_Plane_m, only: Line_And_Plane
    use QTM_m, only: Stack_t
    use QTM_Interpolation_Weights_3D_m, only: S_QTM_t, Cone_Face, &
      & Top_Face, X_Face, Y_Face
    use Radius_of_Curvature_m, only: Radius_of_Curvature_Mean
    use Sort_m, only: SortP

    type(ECR_t), intent(inout) :: Line(2)  ! The line is Line(1) + s * Line(2);
                                           ! Line(2) is made an unit vector here
    real(rg), intent(in) :: SMin, SMax     ! Reject intersections outside this range
    type(QTM_tree_t), intent(inout) :: QTM_Tree
    real(rg), intent(in), contiguous :: H(:,:) ! Heights X size(QTM_Tree%Geo_in)
    type(S_QTM_t), intent(out), allocatable :: Intersections(:) ! S-values of
                                           ! intersections of Line with any
                                           ! surface of any prism of the
                                           ! QTM-based 3D grid, provided they
                                           ! are within [SMin, SMax]
    integer, intent(in) :: Facets(:)       ! Facets under Line
    integer, intent(in), optional :: Which ! Which intersections to detect,
                                           ! default all

    type(S_QTM_t), allocatable :: Cone_Int(:) ! All intersections of Line with
                                           ! a latitude cone of the QTM
    integer :: I, J, K
    type(S_QTM_t), allocatable :: Inside(:) ! [ Cone_Int, V_Int, Top_Int ]
    integer :: MyWhich                     ! Which intersections to detect
    integer :: N_Cone                      ! Number of intersections of Line
                                           ! with a latitude cone of the QTM
    integer :: N_Height                    ! Ubound(H,1)
    integer :: N_Top                       ! Number of intersections of Line with
                                           ! a horizontal boundary surface
    integer :: N_Vert                      ! Number of intersections of Line
                                           ! with a vertical face of a prism of
                                           ! the QTM that is not on a latitude
                                           ! cone
    integer, allocatable :: P(:) ! Permutation vector for sorting
    type(stack_t) :: Stack                 ! To make QTM searches faster
    type(S_QTM_t) :: Tangent ! SMin or SMax, whichever has the smaller magnitude
    type(S_QTM_t) :: Top_Int(2*ubound(h,1)) ! All intersections of Line with the
                                           ! sphere defined by a facet, the
                                           ! mean radius of curvature at its
                                           ! circumcenter, and the geodetic
                                           ! heights at its vertices.  We use
                                           ! the mean radius of curvature
                                           ! instead of one in the direction of
                                           ! Line so that it will be the same
                                           ! curvature, and therefore the same
                                           ! boundary surface, for all lines.
    type(S_QTM_t), allocatable :: V_Int(:) ! All intersections of Line with a
                                           ! vertical face of a prism of the
                                           ! QTM that is not on a latitude
                                           ! cone

    myWhich = all
    if ( present(which) ) myWhich = which

    n_height = ubound(h,1)

    ! Get all intersections of Line with horizontal boundary surfaces
    ! of prisms of the QTM.  These are spheres that have the same radius
    ! of curvature as the Earth's surface at the circumcenter of the
    ! QTM facet, and pass through the three points at the height level
    ! above the facet.

    n_top = 0
    if ( myWhich /= vertical ) call Intersect_Line_And_Horizontal_Boundary

    if ( myWhich /= horizontal ) then

      ! Get all intersections of Line with latitude cones of the QTM.
      call Intersect_Line_And_Latitude_Cone

      ! A vertical face of a prism can be intersected by the line only if its
      ! cone face is or one of its horizontal boundaries is.  The intersected
      ! faces might be in the same prism, the one above it, the one below it, or
      ! the ones below and above it on the opposite facet of its latitude cone.
      ! For now, it's simpler just to check all the vertical faces than to check
      ! the eight possible faces while avoiding duplicates.
      call Intersect_Line_And_Vertical_Boundary
    else
      allocate ( cone_int(1:0) )
      allocate ( v_int(1:0) )
    end if

    ! The absolute value of one of sMin or sMax is sqrt(huge(0.0_rg)).
    ! The other one is the value of S at the tangent or Earth-surface
    ! intersection.
    tangent%s = sMin
    if ( abs(sMax) < abs(sMin) ) tangent%s = sMax

    ! Make sure the tangent is in the list of intersections
    if ( any ( tangent%s == cone_int(1:n_cone)%s ) .or. &
       & any ( tangent%s == top_int(1:n_top)%s ) .or. &
       & any ( tangent%s == v_int(1:n_vert)%s ) ) then
      allocate ( inside, source = [ cone_int(1:n_cone), &
                                  & top_int(1:n_top), &
                                  & v_int(1:n_vert) ] )
    else
      allocate ( inside, source = [ cone_int(1:n_cone), &
                                  & top_int(1:n_top), &
                                  & v_int(1:n_vert), tangent ] )
    end if

    ! Sort the intersections found so far according to their positions along
    ! Line
    k = size(inside,1)
    allocate ( p(k+n_height) )
    call sortp ( inside%s, 1, k, p )
    inside = inside(p(:k))
    ! Maybe at this point we want to eliminate points that are too
    ! close together to bother keeping both of them

    ! Extrapolate along Line to surfaces of constant height that are above
    ! the intersection of Line with the outermost face of a prism of the QTM.
    ! The heights used are the ones at the longitude and geodetic latitude
    ! of the intersection.

    ! Then concatenate the outside intersections with the inside ones,
    ! and sort everything on S.

    if ( myWhich /= vertical ) then
      call Intersect_Line_And_Extrapolated_Height
    else
      call move_alloc ( inside, intersections )
    end if

  contains

    subroutine Intersect_Line_And_Extrapolated_Height
      use Geolocation_0, only: Norm2
      use Line_And_Ellipsoid_m, only: Line_And_Ellipsoid
      use QTM_m, only: Geo_To_ZOT, ZOT_t
      use Triangle_Interpolate_m, only: Triangle_Interpolate
      type(ECR_t) :: Edge      ! ECR coordinates of Inside(I_edge)
      real(rg) :: Eta(3)       ! Interpolation coefficients to compute height
      integer :: F             ! Index of facet containing Edge or Geod
      type(H_V_Geod) :: Geod   ! Geodetic coordinates of Inside(I_edge)
      integer :: I_Edge        ! Index in sorted Inside of edge of QTM polygon
      integer :: I_H(3)        ! Second subscripts of H
      integer :: N_H           ! Number of QTM vertices to use for height
                               ! interpolation
      integer :: N_Out         ! How much of Outside is used
      type(S_QTM_t) :: Outside(size(h,1)) ! Intersections with surfaces of
                               ! constant height outside the QTM
      real(rg), allocatable :: S(:) ! Intersections of Line with
                               ! a surface at constant height above the Earrh
                               ! reference ellipsoid.  Size is in 0..2.
      type(ECR_t) :: Surf      ! ECR coordinates at the surface
      type(ECR_t), allocatable :: Test_Int(:) ! Intersections of Line with
                               ! a surface at constant height above the Earrh
                               ! reference ellipsoid.  Size is in 0..2.
      real(rg) :: Want_H       ! The desired height at the extrapolation
      type(ZOT_t) :: Z         ! ZOT coordinates of Geod
      integer :: Ser(3)        ! Vertices of facet to use for height
                               ! interpolation

      if ( tangent%s < 0.0 ) then
        i_edge = maxloc(inside%s,1)
      else
        i_edge = minloc(inside%s,1)
      end if
      edge = line(1) + inside(i_edge)%s * line(2)
      geod = edge%geod()
      f = QTM_tree%find_facet ( geod, stack )
      i = 6 - QTM_tree%q(f)%xn - QTM_tree%q(f)%yn ! Polar vertex index
      if ( inside(i_edge)%face /= top_face ) then
        select case ( inside(i_edge)%face )
        case ( cone_face )
          ser(1) = QTM_tree%q(f)%xn
          ser(2) = QTM_tree%q(f)%yn
        case ( x_face )
          ser(1) = QTM_tree%q(f)%xn
          ser(2) = i
        case ( y_face )
          ser(1) = i
          ser(2) = QTM_tree%q(f)%yn
        end select
        n_h = 2 ! Number of QTM vertices to use for height interpolation
        ser(1:2) = QTM_tree%q(f)%ser(ser(1:2))
        ! The parent type of H_v_geod is H_v_t, not H_Geod
        surf = geod%surf_ECR()
        eta(1) = norm2(QTM_tree%geo_in(ser(2))%ecr() - surf) / &
               & norm2(QTM_tree%geo_in(ser(2))%ecr() - &
                      &QTM_tree%geo_in(ser(1))%ecr())
        eta(2) = 1 - eta(1)
        i_h(1) = min(ser(1),size(h,2))
        i_h(2) = min(ser(2),size(h,2))
      else ! This interpolates in a plane; maybe it should be the same
           ! as in Intersect_Line_And_Horizontal_Boundary
        n_h = 3 ! Number of QTM vertices to use for height interpolation
        ser = QTM_tree%q(f)%ser
        i_h = min(QTM_tree%q(f)%ser,size(h,2))
        ! Compute the interpolation coefficients in ZOT coordinates.
        z = geo_to_ZOT ( geod )
        call triangle_interpolate ( QTM_tree%q(f)%z%x, QTM_tree%q(f)%z%y, &
          & z%x, z%y, eta )
      end if
      n_out = 0
      do i = 1, n_height
        want_h = sum ( h(i,i_h(:n_h)) * eta(:n_h) )
        ! Compute the intersections of Line with a surface at height Want_H
        ! above the Earth reference ellipsoid, +/- Height_Tol.
        call line_and_ellipsoid ( line, want_h, height_tol, test_int, s )
        do j = 1, size(test_int)
          if ( s(j) >= sMin .and. s(j) <= sMax ) then
            n_out = n_out + 1
            outside(n_out) = S_QTM_t( s=s(j), face=-inside(i_edge)%face,  &
                                    & facet=f, h=want_h, h_ind=j, n_coeff=n_h )
            outside(n_out)%coeff(:n_h) = eta(:n_h)
            outside(n_out)%ser(:n_h) = ser(:n_h)
          end if
        end do
      end do
      ! Concatenate the outside intersections with the inside ones
      allocate ( intersections, source = [ inside, outside(:n_out) ] )
      ! Sort intersections according to S
      call sortP ( intersections%s, 1, k+n_out, P )
      intersections = intersections(p(:k+n_out))
    end subroutine Intersect_Line_And_Extrapolated_Height

    subroutine Intersect_Line_And_Horizontal_Boundary
      ! Get all intersections of Line with horizontal boundary surfaces
      ! of prisms of the QTM.  These are spheres that have the same radius
      ! of curvature as the Earth's surface at the circumcenter of the
      ! QTM facet, and pass through the three points at the height level
      ! above the facet.
      type(ECR_t) :: C                  ! Vector from V(3) to circumcenter
                                        ! of facet V (see below)
      type(ECR_t) :: CC                 ! Circumcenter of facet V (see
                                        ! below) = V(3) + C, or its centroid
      type(ECR_t) :: Center             ! Center of sphere containing facet V
      integer :: F                      ! Index of a facet, from Facets(:)
      real(rg) :: G(3)                  ! Normalized geocentric angular
                                        ! distance coordinates
      type(h_v_geod) :: Geod            ! Geodetic coordinates of CC
      type(h_v_geod) :: Geod_f(3)       ! Geodetic coordinates of corners
                                        ! of facet                                        
      integer :: I, J, K, L, M
      logical :: Keep                   ! "Keep the intersection"
      type(ECR_t) :: N                  ! Normal vector to V (see below)
      real(rg) :: N2                    ! |N|**2
      real(rg) :: R                     ! Mean radius of curvature at the
                                        ! geodetic latitude of the
                                        ! circumcenter of a horizontal
                                        ! boundary surface of a facet
      real(rg), allocatable :: S_Int(:) ! Intersections of Line with the
                                        ! horizontal boundary of one facet
                                        ! of the QTM
      type(ECR_t) :: V(3)               ! Vertices of a horizontal boundary
                                        ! surface above a QTM facet

      do i = 1, size(facets)
        f = facets(i)
        ! We know qtm_tree%q(f)%depth == qtm_tree%level, f.e., it's a facet
        ! at the finest level of refinement.
        do j = 1, n_height

          do k = 1, 3
          ! Get geodetic coordinates of vertices of QTM facet V at height
          ! H(j,.) above the Earth's surface.
            l = min(qtm_tree%q(f)%ser(k),size(h,2)) ! Coherent?
            geod_f(k) = h_v_geod ( QTM_tree%geo_in(qtm_tree%q(f)%ser(k))%lon, &
                                 & QTM_tree%geo_in(qtm_tree%q(f)%ser(k))%lat, &
                                 & h(j,l) )
            ! Get the ECR coordinates of geod_f
            v(k) = geod_f(k)%ecr()
          end do
          ! Get the vector C from V(3) to the circumcenter of facet V, a
          ! normal N to the plane containing V, and N2 = |N|^2
          call circumcenter ( v, c, n, n2 )
          ! Get the circumcenter of facet V
          cc = c + v(3)
          ! Get the geodetic coordinates of the circumcenter of facet V
          geod = cc%geod()
          ! Get the radius of curvature at CC.  We assume that the geodetic
          ! height of the boundary surface is so small that the radius of
          ! curvature is essentially the same as at the Earth's surface.
          r = radius_of_curvature_mean ( geod%lat )
          ! Get the center of the sphere having mean radius of curvature R
          ! and including V, and center nearest the Earth's center
          call center_of_sphere ( v(3), cc, n, n2, r, center )
          ! Compute intersections of Line with the sphere at Center
          ! and radius R
          call line_and_sphere ( r, line, s=s_int, center=center )
          ! Eliminate intersections not in [SMin,SMax]
          m = size(s_int)
          do k = 1, m
            keep = s_int(k) >= sMin .and. s_int(k) <= sMax
            if ( keep ) then
              ! Get the geodetic coordinates of the centroid of facet V
              cc = ( v(1) + v(2) + v(3) ) / 3.0_rg
              geod = cc%geod()
              ! Keep the intersection if it's with the current facet.
              keep = f == QTM_tree%find_facet ( geod, stack )
            end if
            if ( keep ) then
              n_top = n_top + 1
              top_int(n_top) = S_QTM_t(s=s_int(k), facet=f, h=geod%v, &
                             & face=top_face, n_coeff=3, h_ind=j,     &
                             & ser=qtm_tree%q(f)%ser )
              do l = 1, 3 ! Compute horizontal interpolation coefficients.
                          ! Use normalized geocentric angular distances.
                          ! C and CC are temporary variables here
                cc = line(1) + s_int(k) * line(2)
                cc = cc / cc%norm2()
                c = QTM_tree%geo_in(qtm_tree%q(f)%ser(l))%ecr(norm=.true.)
                g(f) = acos( cc .dot. c )
              end do
              top_int(n_top)%coeff(1) = 0.5 * ( 1.0 - g(2) - g(3) )
              top_int(n_top)%coeff(2) = 0.5 * ( 1.0 - g(1) - g(3) )
              top_int(n_top)%coeff(3) = 0.5 * ( 1.0 - g(1) - g(2) )
            end if
          end do
        end do
      end do
    end subroutine Intersect_Line_And_Horizontal_Boundary

    subroutine Intersect_Line_And_Latitude_Cone
      ! Get all intersections of Line with latitude cones of the QTM.
      real(rg) :: Eta_h(2)    ! Interpolation coefficient for longitude
      integer :: F            ! Index of a facet, from Facets(:)
      class(h_v_t), allocatable :: Geo     ! H_V_Geoc or H_V_Geod, depending
                              ! upon QTM_tree%geo_in
      integer :: H1, H2       ! Subscripts for second dimension
                              ! of H at non-polar vertices of F
      integer :: I, J
      logical :: Keep         ! "Keep the intersection"
      real(rg), allocatable :: S_Int(:)    ! Intersections of Line with
                              ! one latitude cone of the QTM.
      integer :: S1, S2       ! Serial numbers of X or Y
                              ! vertices of QTM facet

      allocate ( cone_int(2*size(QTM_Tree%QTM_lats)) )
      n_cone = 0
      select type ( Q => QTM_tree%geo_in )
      class is ( h_geod )
        allocate ( h_v_geod :: geo )
      class default
        allocate ( h_v_geoc :: geo )
      end select

      do i = 1, size(QTM_Tree%QTM_lats)
        call line_and_cone ( QTM_Tree%QTM_lats(i), line, s=s_int )
        ! Eliminate intersections outside the QTM or not in [SMin,SMax]
        keep = .true.
        do j = 1, size(s_int) ! size(s_int) is in 0..2 here
          ! Convert ECR coordinates of intersection to lon and geoc/geod lat
          ! because that's what Find_Facet needs
          call geo%from_ECR(line(1) + s_int(j) * line(2))
          f = QTM_tree%find_facet ( geo, stack )
          ! Intersection with cone is on a QTM leaf facet that is within
          ! or intersects the polygon, and is within [SMin,SMax]
          keep = qtm_tree%q(f)%depth == qtm_tree%level .and. &
               & s_int(j) >= SMin .and. s_int(j) <= SMax
          if ( keep ) then
            s1 = qtm_tree%q(f)%ser(qtm_tree%q(f)%xn)
            h1 = min(s1,ubound(h,2))
            s2 = qtm_tree%q(f)%ser(qtm_tree%q(f)%yn)
            h2 = min(s2,ubound(h,2))
            eta_h(1) = ( QTM_tree%geo_in(s2)%lon%d - geo%lon%d ) / &
                     & ( QTM_tree%geo_in(s2)%lon%d - QTM_tree%geo_in(s1)%lon%d )
            eta_h(2) = 1 - eta_h(2)
            ! Keep the intersection if it's not too high or too low
            keep = &
              & geo%v >= ( h(1,h1) * eta_h(1) + h(1,h2) * eta_h(2) ) .and. &
              & geo%v <= ( h(n_height,h1) * eta_h(1) + h(n_height,h2) * eta_h(2) )
          end if
          if ( keep ) then
            n_cone = n_cone + 1
            cone_int(n_cone)%s = s_int(j)
            cone_int(n_cone)%facet = f
            cone_int(n_cone)%h = geo%v
            cone_int(n_cone)%face = cone_face
            cone_int(n_cone)%n_coeff = 2
            cone_int(n_cone)%ser(:2) = qtm_tree%q(f)%ser([s1,s2])
            cone_int(n_cone)%coeff(:2) = eta_h
           end if
        end do
      end do
    end subroutine Intersect_Line_And_Latitude_Cone

    subroutine Intersect_Line_And_Vertical_Boundary
      ! Get all intersections of Line with vertical faces of prisms of the QTM
      ! that are not on latitude cones.  These are faces for which one vertex
      ! of the QTM is a polar vertex and the other one is not.
      real(rg) :: Eta_h(2)    ! Interpolation coefficient for latitude
      integer :: F            ! Index of a facet, from Facets(:)
      integer, parameter :: Faces(2) = [ X_face, Y_face ]
      type(h_v_geod) :: Geod_f(2,2) ! Geodetic coordinates of a point on a plane
      type(h_v_geod) :: Geod_int    ! Geodetic coordinates of intersection
      logical :: Hit(n_height,QTM_tree%n_in) ! Each vertex is adjacent to only
                              ! one polar vertex, so its serial number, along
                              ! with the height, can be used to identify a face
      real(rg) :: H1, H2      ! Heights at Eta on top and bottom boundaries
      integer :: I
      type(ECR_t) :: Int      ! Intersection, if there is one
      integer :: Intersect    ! -1 => no intersection,
                              !  0 => line is in the plane,
                              ! +1 => one intersection.
      integer :: J, K, L, M, N
      type(ECR_t) :: Plane(4) ! Vertices that define a vertical face
      integer :: SN(0:1)      ! Serial numbers of QTM vertices of vertical face
      real(rg) :: S_Int       ! S-value of the intersection, if there is one
      integer :: V(2,2)       ! Indices of QTM vertices of edges not on a
                              ! latitude cone

      allocate ( v_int(2*size(QTM_Tree%QTM_lats)) )
      n_vert = 0
      hit = .false.
   facet: do i = 1, size(facets) ! Examine all the facets under the path
        f = facets(i)
        ! F is the index of a facet of the finest refinement.
        ! Get indices of vertices of edges not on a latitude cone.
        ! One of the vertices of such an edge will be the pole node.
        v(1,:) = 6 - qtm_tree%q(f)%xn - qtm_tree%q(f)%yn ! pole node
        v(2,1) = qtm_tree%q(f)%xn
        v(2,2) = qtm_tree%q(f)%yn
        do j = 1, size(h,1)-1
          do k = 1, 2     ! Which vertical plane
            if ( hit(j,v(2,k)) ) cycle facet ! Already checked this face
            do m = 1, 2   ! 1 => polar, 2 => non-polar vertex of the QTM
              ! Get geodetic and then ECR coordinates of vertices of vertical
              ! plane K at heights H(j:j+1,.) above the Earth's surface.
              do n = 0, 1
                sn(n) = min(qtm_tree%q(f)%ser(v(m,k)),size(h,2))
                geod_f(m,n+1) = h_v_geod ( QTM_tree%geo_in(qtm_tree%q(f)%ser(v(m,k)))%lon, &
                                         & QTM_tree%geo_in(qtm_tree%q(f)%ser(v(m,k)))%lat, &
                                         & h(j+n,sn(n)) )
                plane(2*n+m) = geod_f(m,n+1)%ECR()
              end do ! n
              ! Geod_f(1,:) are on vertical edge at polar vertex
              ! Geod_f(2,:) are on vertical edge at non-polar
              ! Geod_f(:,1) are on bottom surface
              ! Geod_f(:,2) are on top surface
            end do ! m
            call line_and_plane ( plane(1:3), line, intersect, int, s_int )
            if ( s_int < sMin.or. s_int > sMax ) cycle ! Outside range
            if ( intersect >= 0 ) then
              if ( intersect == 0 ) then ! Place intersection at the centroid
                int = 0.25_rg * ( plane(1) + plane(2) + plane(3) + plane(4) )
                l = maxloc(abs(line(2)%xyz),1) ! don't divide by zero; assume
                                               ! line(2) is not the zero vector
                s_int = ( int%xyz(l) - line(1)%xyz(l) ) / line(2)%xyz(l)
              end if
              call geod_int%from_ECR(int) ! Get geodetic coordinates of Int
              ! We're looking at a face that's not on a parallel of latitude.
              ! Therefore, if the line intersects the face, the latitude of
              ! the intersection will be between the minimum and maximum
              if ( geod_int%lat >= minval(geod_f(:,1)%lat) .and. &
                 & geod_int%lat <= minval(geod_f(:,1)%lat) ) then
                ! Intersection is within horizontal range.  Interpolate
                ! heights along top and bottom edges to geod_int%lat.
                eta_h(1) = ( geod_f(2,1)%lat - geod_int%lat ) / &
                         & ( geod_f(2,1)%lat - geod_f(1,1)%lat )
                eta_h(2) = 1 - eta_h(2)
                h1 = eta_h(1) * geod_f(2,1)%v + eta_h(2) * geod_f(1,1)%v
                h2 = eta_h(1) * geod_f(2,2)%v + eta_h(2) * geod_f(1,2)%v
                if ( geod_int%v >= h1 .and. geod_int%v <= h2 ) then
                  ! Intersection is within vertical range too
                  n_vert = n_vert + 1
                  v_int(n_vert) = S_QTM_t( s=s_int, facet=f, h=geod_int%v, &
                                         & face=faces(m), n_coeff=2 )
                  v_int(n_vert)%ser(:2) = qtm_tree%q(f)%ser(v(:,1))
                  ! Interpolation coefficients are in latitude only.
                  v_int(n_vert)%coeff(:2) = eta_h
                  hit(j,v(2,k)) = .true.
                end if
              end if
            end if
          end do ! which plane
        end do ! heights
      end do facet
    end subroutine Intersect_Line_And_Vertical_Boundary

  end subroutine Metrics_3D_QTM_1

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Metrics_3D_m

! $Log$
! Revision 2.8  2016/11/09 00:36:13  vsnyder
! Remove Metrics_3D_Grid, pass in list of facets to use
!
! Revision 2.7  2016/11/03 19:11:47  vsnyder
! Inching toward 3D forward model
!
! Revision 2.6  2016/10/05 23:29:04  vsnyder
! Replace ZOT_n component name with Ser because it's a serial number for
! more than just the ZOT coordinates.
!
! Revision 2.5  2016/09/30 21:47:47  pwagner
! Intel v16 unhappy w/o these changes
!
! Revision 2.4  2016/09/23 18:42:53  vsnyder
! Inching along
!
! Revision 2.3  2016/09/14 18:18:05  vsnyder
! Replace QTM_Node_t%Leaf with Depth
!
! Revision 2.2  2016/05/10 00:13:59  vsnyder
! Remember all three QTM vertex indices used for height interpolation.
!
! Revision 2.1  2016/04/16 02:06:23  vsnyder
! Initial commit
!
