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
module QTM_Interpolation_Weights_3D_m
!=============================================================================

  use Geolocation_0, only: RG, S_t
  use QTM_Interpolation_Weights_m, only: QTM_Interpolation_Weights
  use Path_Representation_m, only: Value_t
  use Weight_1D_m, only: Weight_1D

  implicit NONE
  private

  ! Use QTM_Interpolation_Weights to get weights in 2D, and linear
  ! interpolation to get weights in 3D.  It is assumed that the 3D
  ! grid is stacked and coherent.

  type, extends(S_t), public :: S_QTM_t ! Descriptors of points along a line
                     ! through a stacked but not necessarily coherent 3D grid
                     ! built atop a QTM surface grid.  The grid consists of
                     ! triangular prisms resting upon the QTM.  The horizontal
                     ! boundaries of the prisms are spherical caps.  One
                     ! vertical face is a latitude cone; the other two are
                     ! planes.
  ! real(rg) :: S    ! Distance along Line(1) in the direction of Line(2);
                     ! Inherited from parent type.
    real(rg) :: Coeff(3) = 0.0_rg ! Horizontal interpolation coefficients.
    integer :: Face  ! Cone_Face => latitude cone face, bounded horizontally
                     !      by two non-polar vertices of facet QID.
                     ! X_face => vertical face bounded horizontally by the
                     !      polar vertex and the X-node vertex of facet QID.
                     ! Y_face => vertical face bounded horizontally by the
                     !      polar vertex and the Y-node vertex of facet QID.
                     ! Top_Face => horizontal face of facet QID at height
                     !      indexed by H.
                     ! Inside_Prism => The intersection is inside a prism
                     !      resting on the surface QTM instead of with one
                     !      of its faces.
                     ! If the point is outside the QTM, this will be the
                     ! negative of one of the above values, other than
                     ! Inside_Prism; it will be whatever face the extrapolated
                     ! line intersects.  The coordinates indexed by Ser are
                     ! those where the line intersects the QTM, not where it
                     ! intersects a surface of a height equal to one of the
                     ! heights in the intersection of all heights in the state
                     ! vector.  The H component probably isn't useful.
    integer :: Facet = 0 ! Index in Grid%QTM_Tree of QTM facet intersected at
                     ! the point.
    real(rg) :: H    ! Height at which the point intersects a face.  This is
                     ! the height where the extrapolated line intersects a face
                     ! of the QTM for points outside the QTM.
    integer :: H_ind = 0 ! Index in height grid.  If |Face| /= Top_Face, H_Ind
                     ! is the index of the next lower height surface.  If
                     ! |Face| == Top_Face but the point is below the first
                     ! element of the height grid, H_ind is zero.  This could
                     ! happen if the ray is an Earth-reflecting ray and the
                     ! minimum height in the height-reference array is above
                     ! the Earth surface.
    integer :: N_Coeff = 0 ! How many elements of Coeff and Ser are used.
    integer :: PV(3) = 0   ! Indices of vertices near path, from 
                     ! QTM%path_vertices(ser), second subscript for heights.
    integer :: Ser(3) = 0  ! Serial numbers of facet coordinates, see QTM_Node_t
                     ! in the Generate_QTM_m module.
  end type S_QTM_t

  integer, parameter, public :: MaxCoeff = 6 ! First dimension of Weights arguments

  ! Interpolation coefficients within a 3D grid that has a QTM horizontal grid
  type, public :: Weight_ZQ_t
    integer :: N                  ! count(w%v /= 0), <= MaxCoeff
    type(value_t) :: W(maxCoeff)  ! Weights, and array element positions in
                                  ! the state vector or at points near the path
  end type Weight_ZQ_t

  integer, parameter, public :: Cone_Face = 1
  integer, parameter, public :: X_face = cone_face + 1
  integer, parameter, public :: Y_face = x_face + 1
  integer, parameter, public :: Top_Face = y_face + 1
  integer, parameter, public :: Inside_Prism = top_face + 1

  public :: QTM_Interpolation_Weights

  interface QTM_Interpolation_Weights
    module procedure &
      & QTM_Interpolation_Weights_Geo_3D, &
      & QTM_Interpolation_Weights_Geo_3D_Incoherent, &
      & QTM_Interpolation_Weights_Geo_3D_Incoherent_Line, &
      & QTM_Interpolation_Weights_Geo_3D_Incoherent_List, &
      & QTM_Interpolation_Weights_Geo_3D_Line, &
      & QTM_Interpolation_Weights_Geo_3D_List
  end interface QTM_Interpolation_Weights

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine QTM_Interpolation_Weights_Geo_3D ( QTM_Tree, Heights, Point, &
                                              & Weights, Facet, Stack, Used )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  Weights%N are QTM vertex serial numbers.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t, H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK, Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)   ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Point    !   the same units as for Point%V
    type(weight_ZQ_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    class(h_t), intent(out), optional :: Used ! Horizontal coordinates of the
                                         !   point used for interpolation

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: Index, J
    integer(qk) :: N_QTM ! Number of vertices of QTM within or adjacent to
                         ! the polygon
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    n_QTM = QTM_Tree%n_in

    ! Get horizontal interpolation coefficients and serial numbers
    call QTM_Interpolation_Weights ( QTM_Tree, Point, Weights%w(1:3), Facet, &
      & Stack, Used )

    ! Get vertical interpolation coefficients
    call hunt ( heights, point%v, index, &
              & allowTopValue=.true., allowBelowValue = .true. )

    ! Combine horizontal and vertical interpolation coefficients
    if ( index == 0 ) then
      weights%w(4:6)%n = 0 ! Constant vertical extrapolation outside
                           !   below of Heights(:)
    else if ( index >= size(heights) ) then
      do j = 1, 3
        weights%w(j)%nz = size(heights)
      end do
      weights%w(4:6)%n = 0 ! Constant vertical extrapolation outside
                           !   above of Heights(:)
    else
      dH = heights(index+1) - heights(index)
      if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
        weights%w(4:6)%n = 0 ! Don't use points other than at Height Index
      else
        w = ( heights(index+1) - point%v ) / dh
        if ( w == 1.0_rg ) then
          do j = 1, 3
            weights%w(j)%nz = index
          end do
          weights%w(4:6)%n = 0 ! Don't use points other than at Height Index
        else if ( w == 0.0_rg ) then
          do j = 1, 3
            weights%w(j)%nz = index + 1
          end do
          weights%w(4:6)%n = 0 ! Don't use points other than at Height Index + 1
        else
          weights%w(4:6)%v = weights%w(1:3)%v * ( 1.0_rg - w )
          weights%w(1:3)%v = weights%w(1:3)%v * w
          do j = 1, 3
            weights%w(j+3)%nz = index + 1
            weights%w(j)%nz = index
          end do
        end if
      end if
    end if

    where ( weights%w%v == 0 ) weights%w%n = 0
    weights%n = count ( weights%w%n /= 0 ) 
    weights%w = pack ( weights%w,weights%w%n /= 0 )
    weights%w(weights%n+1:) = value_t(0,0,0,0)

  end subroutine QTM_Interpolation_Weights_Geo_3D

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent ( QTM_Tree, Heights, &
                                     & Point, Weights, Facet, Stack, Used )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  Weights%N are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t, H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK, Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! from which to interpolate, extents
                                         ! are (# heights, # vertices near path)
                                         ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Point    !   the same units as for Point%V
    type(weight_ZQ_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    class(h_t), intent(out), optional :: Used ! Horizontal coordinates of the
                                         !   point used for interpolation

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: Index
    integer(qk) :: N_QTM ! Number of vertices of QTM within or adjacent to
                         ! the polygon
    integer :: P         ! Index of a corner of a facet
    integer :: S         ! Serial number of a point in the QTM
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, Heights(:,1), &
                                     & Point, Weights, Facet, Stack, Used )
      weights%w%np = 1
      return
    end if

    n_QTM = QTM_Tree%n_in

    ! Get horizontal interpolation coefficients and serial numbers
    call QTM_Interpolation_Weights ( QTM_Tree, Point, Weights%w(1:3), Facet, &
      & Stack, Used )

    ! Default is constant vertical extrapolation outside range of Heights(:)
    weights%w(4:6)%n = 0
    ! Get vertical interpolation coefficients
    do p = 1, 3
      index = 0
      s = weights%w(p)%n
      if ( s /= 0 ) then
        s = QTM_tree%path_vertices(s)
        weights%w(p)%np = s
        call hunt ( heights(:,s), point%v, index, &
                  & allowTopValue=.true., allowBelowValue = .true. )
      end if

      ! Combine horizontal and vertical interpolation coefficients
      if ( index /= 0 ) then
        dH = heights(index+1,s) - heights(index,s)
        if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
          weights%w(p+3)%n = 0 ! Don't use points other than at Height Index
          weights%w(p+3)%nz = 0
        else
          w = ( heights(index+1,p) - point%v ) / dh
          if ( w == 1.0_rg ) then
            weights%w(p)%nz = index
          else if ( w == 0.0_rg ) then
            weights%w(p)%nz = index+1
          else
            ! Same horizontal coefficients at both levels
            weights%w(p+3)%n = weights%w(p)%n
            weights%w(p+3)%np = weights%w(p)%np
            weights%w(p+3)%v = weights%w(p)%v * ( 1.0_rg - w )
            weights%w(p)%v = weights%w(p)%v * w
            weights%w(p+3)%nz = index + 1
            weights%w(p)%nz =index
          end if
        end if
      end if

    end do ! p

    where ( weights%w%v == 0 ) weights%w%n = 0 ! Don't use zero weights
    weights%n = count ( weights%w%n /= 0 )
    weights%w = pack ( weights%w,weights%w%n /= 0 )
    weights%w(weights%n+1:) = value_t(0,0,0,0)

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line ( QTM_Tree, &
                         & Heights, Line, Points, Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  Weights%N are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.
    ! The points are described by Line, which consists of a vector in ECR to
    ! a point on a line and an unit vector along the line, and Points, which
    ! describe distances along that line and how the line intersects a prism
    ! of a QTM-based grid (or not).

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_V_Geod

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! from which to interpolate.  Extents
                                         ! are (# heights, # points near path).
    type(ECR_t), intent(in) :: Line(2)   ! Vector to a point on the line,
                                         ! vector along the line.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(weight_ZQ_t), intent(out) :: Weights(size(points))

    type(H_V_Geod) :: G                  ! Geodetic coordinate of Points(i)
    integer :: I, J
    integer :: Jlo, Jhi                  ! Indices of H[pxy] in Heights
    integer :: Nh                        ! Number of heights
    integer :: NP(3)                     ! Second subscripts of Heights
    type(ECR_t) :: P                     ! ECR coordinate of Points(i)
    real(rg) :: W(2,2)                   ! Interpolation weights

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, Heights(:,1), &
                                     & Line, Points, Weights )
      return
    end if

    nh = size(heights,1)

    do i = 1, size(points)
      p = line(1) + points(i)*line(2)
      g = p%geod()
      select case ( abs(points(i)%face) )
      case ( cone_face, X_face, Y_face )
        ! Interpolate in height at the X and Y, X and P, or Y and P nodes
        ! of the QTM.
        jlo = 1; jhi = nh
        np(1:2) = QTM_tree%path_vertices(points(i)%ser(1:2))
        call weight_1d ( heights(:,np(1)), g%v, jlo, jhi, w(:,1) )
        call weight_1d ( heights(:,np(2)), g%v, jlo, jhi, w(:,2) )
        ! Assume that it's extremely rare to hit a vertex, so compute four
        ! weights.
        weights(i)%n = 4
        weights(i)%w(1:maxCoeff)%v = [ w(:,1)*points(i)%coeff(1), &
                                     & w(:,2)*points(i)%coeff(2), &
                                     & 0.0_rg, 0.0_rg ]
        weights(i)%w(1:maxCoeff)%n = &
          & [ points(i)%ser(1), points(i)%ser(1), &
            & points(i)%ser(2), points(i)%ser(2), 0, 0 ]
        weights(i)%w(1:maxCoeff)%np = [ np(1), np(1), np(2), np(2), 0, 0 ]
        weights(i)%w(1:maxCoeff)%nz = [ jlo, jhi, jlo, jhi, 0, 0 ]

      case ( inside_prism ) ! Use the general method
        call QTM_interpolation_weights ( QTM_tree, heights, g, weights(i) )
      case ( top_face )
        ! Get the two vertical interpolation weights
        do j = 1, 3
          jlo = 1; jhi = nh
          np(j) = QTM_tree%path_vertices(points(i)%ser(j))
          call weight_1d ( heights(:,np(j)), points(i)%h, jlo, jhi, w(:,1) )
          weights(i)%n = 6
          ! Compute the outer product of the vertical and horizontal
          ! interpolation coefficients
          weights(i)%w(2*j-1:2*j)%v = w(:,1)*points(i)%coeff(j)
          weights(i)%w(2*j-1:2*j)%n = points(i)%ser(j)
          weights(i)%w(2*j-1:2*j)%np = np(j)
          weights(i)%w(2*j-1:2*j)%nz = [ jlo, jhi ]
        end do
      end select

      where ( weights(i)%w%v == 0 ) weights(i)%w%n = 0 ! Don't use zero weights(i)
      weights(i)%n = count ( weights(i)%w%n /= 0 )
      weights(i)%w = pack ( weights(i)%w,weights(i)%w%n /= 0 )
      weights(i)%w(weights(i)%n+1:) = value_t(0,0,0,0)

    end do

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_List ( QTM_Tree, &
                         & Heights, Points, Path_Heights, Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  Weights%N are array element
    ! order positions, which can be converted to 2D subscripts (height, QTM
    ! serial number) by the Subscripts function in the Array_Stuff module.

    use Generate_QTM_m, only: QTM_Tree_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:)  ! from which to interpolate.  Extents
                                          ! are (heights, QTM_Tree%N_In).
    class(s_QTM_t), intent(in) :: Points(:) ! horizontal interpolation coefficients
                                          ! are here
    real(rg), intent(in) :: Path_Heights(:) ! to which to interpolate; same units
                                          ! as heights, same size as Points
    type(weight_ZQ_t), intent(out) :: Weights(size(points))

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: I, Index
    integer(qk) :: N_QTM ! Number of vertices of QTM within or adjacent to
                         ! the polygon
    integer :: NP        ! Second subscript of Heights
    integer :: P         ! Index of a corner of a facet
    integer :: S         ! Serial number of a point in the QTM
    real(rg) :: V        ! Horizontal interpolation coefficient temp
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, heights(:,1), &
                                     & points, path_heights, weights )
      return
    end if

    n_QTM = QTM_Tree%n_in

    do i = 1, size(points)
      ! Get horizontal interpolation weights; serial numbers are gotten later
      weights(i)%w(1:3)%v = points(i)%coeff
      ! Default is constant vertical extrapolation outside range of Heights(:)
      weights(i)%w(4:6)%n = 0
      do p = 1, 3
        index = 0
        s = weights(i)%w(p)%n
        if ( s /= 0 ) then
          np = QTM_tree%path_vertices(s)
          weights(i)%w(p)%np = np
          weights(i)%w(p+3)%np = np
          ! Get vertical interpolation indices
          call hunt ( heights(:,np), path_heights(i), index, &             
                    & allowTopValue=.true., allowBelowValue = .true. )
        end if
        ! Combine horizontal and vertical interpolation coefficients
        if ( index == 0 ) then !   below range of Heights(:)
          weights(i)%w(p)%nz = 1 ! Constant vertical extrapolation
        else if ( index >= size(heights,1) ) then
          weights(i)%w(p)%nz = size(heights,1) ! Constant vertical extrapolation
        else
          dH = heights(index+1,np) - heights(index,np)
          if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
            weights(i)%w(p)%nz = index
          else
            w = ( heights(index+1,np) - path_heights(i) ) / dh
            if ( w == 1.0_rg ) then
              weights(i)%w(p)%nz = index
            else if ( w == 0.0_rg ) then
              weights(i)%w(p)%nz = index+1
            else
              v = weights(i)%w(p)%v
              weights(i)%w(p) = value_t ( n = points(i)%ser(p), np = np, &
                                        & nz = index, &
                                        & v = v * ( 1.0_rg - w ) )
              weights(i)%w(p+3) = value_t ( n = points(i)%ser(p), np = np, &
                                        & nz = index + 1, v = v * w )
            end if
          end if
        end if

      end do ! p

      where ( weights(i)%w%v == 0 ) weights(i)%w%n = 0
      weights(i)%n = count ( weights(i)%w%n /= 0 ) 
      weights(i)%w = pack ( weights(i)%w,weights(i)%w%n /= 0 )
      weights(i)%w(weights(i)%n+1:6) = value_t(0,0,0,0)

    end do ! i

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_List

  subroutine QTM_Interpolation_Weights_Geo_3D_Line ( QTM_Tree, &
                         & Heights, Line, Points, Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  Weights%N are array element order positions,
    ! which can be converted to 2D subscripts (height, QTM serial number) by
    ! the Subscripts function in the Array_Stuff module. The points are
    ! described by Line, which consists of a vector in ECR to a point on a
    ! line and an unit vector along the line, and Points, which describe
    ! distances along that line and how the line intersects a prism of a
    ! QTM-based grid (or not).

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_V_Geod

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)   ! from which to interpolate.  Coherent.
    type(ECR_t), intent(in) :: Line(2)   ! Vector to a point on the line,
                                         ! vector along the line.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(weight_ZQ_t), intent(out) :: Weights(size(points))

    type(h_v_geod) :: G                  ! Geodetic coordinate of Points(i)
    integer :: I
    integer :: Jlo, Jhi                  ! Indices of H[pxy] in Heights
    integer :: Nh                        ! Number of heights
    type(ECR_t) :: P                     ! ECR coordinate of Points(i)
    real(rg) :: W(2)                     ! Interpolation weights

    nh = size(heights,1)

    do i = 1, size(points)
      p = line(1) + points(i)*line(2)
      g = p%geod()
      select case ( abs(points(i)%face) )
      case ( cone_face, X_face, Y_face )
        ! Interpolate in height at the X and Y, X and P, or Y and P nodes
        ! of the QTM.
        ! Get the two vertical interpolation weights
        jlo = 1; jhi = nh
        call weight_1d ( heights, g%v, jlo, jhi, w )
        ! Assume that it's extremely rare to hit a vertex, so compute four
        ! weights.
        weights(i)%n = 4
        ! Compute the outer product of the vertical and horizontal
        ! interpolation coefficients
        weights(i)%w(1:6)%v = [ w*points(i)%coeff(1), &
                              & w*points(i)%coeff(2), 0.0_rg, 0.0_rg ]
        weights(i)%w(1:6)%n = &
          & [ points(i)%ser(1), points(i)%ser(1), &
            & points(i)%ser(2), points(i)%ser(2), 0, 0 ]
        weights(i)%w(1:6)%np = [ QTM_tree%path_vertices(weights(i)%w(1:4)%n), 0, 0 ]
        weights(i)%w(1:6)%nz = [ jlo, jhi, jlo, jhi, 0, 0 ]
      case ( inside_prism ) ! Use the general method
        call QTM_interpolation_weights ( QTM_tree, heights, g, weights(i) )
      case ( top_face )
        ! Get the two vertical interpolation weights
        jlo = 1; jhi = nh
        call weight_1d ( heights, g%v, jlo, jhi, w )
        weights(i)%n = 6
        ! Compute the outer product of the vertical and horizontal
        ! interpolation coefficients
        weights(i)%w(1:6)%v = [ w*points(i)%coeff(1), &
                              & w*points(i)%coeff(2), &
                              & w*points(i)%coeff(3) ]
        weights(i)%w(1:6)%n = &
          & [ points(i)%ser(1),points(i)%ser(1), &
            & points(i)%ser(2),points(i)%ser(2), &
            & points(i)%ser(3),points(i)%ser(3) ]
        weights(i)%w(1:6)%np = QTM_tree%path_vertices(weights(i)%w(1:6)%n)
        weights(i)%w(1:6)%nz = [ jlo, jhi, jlo, jhi, jlo, jhi ]
      end select

      where ( weights(i)%w%v == 0 ) weights(i)%w%n = 0
      weights(i)%n = count ( weights(i)%w%n /= 0 ) 
      weights(i)%w = pack ( weights(i)%w,weights(i)%w%n /= 0 )
      weights(i)%w(weights(i)%n+1:6) = value_t(0,0,0,0)

    end do

  end subroutine QTM_Interpolation_Weights_Geo_3D_Line

  subroutine QTM_Interpolation_Weights_Geo_3D_List ( QTM_Tree, Heights, Points, &
                                                   & Path_Heights, Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  Weights%W%N are array element order positions,
    ! which can be converted to 2D subscripts (height, QTM serial number)
    ! by the Subscripts function in the Array_Stuff module.

    use Generate_QTM_m, only: QTM_Tree_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)    ! from which to interpolate.  Coherent
    class(s_QTM_t), intent(in) :: Points(:) ! horizontal interpolation coefficients
                                          ! are here
    real(rg), intent(in) :: Path_Heights(:) ! to which to interpolate; same units
                                          ! as heights, same size as Points
    type(weight_ZQ_t), intent(out) :: Weights(size(points))

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: I, Index, Indices(size(points)), J
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    ! Get vertical interpolation indices
    call hunt ( heights, path_heights, indices, &
              & allowTopValue=.true., allowBelowValue = .true. )

    ! Combine horizontal and vertical interpolation coefficients
    do i = 1, size(points)
      ! Get horizontal interpolation weights; serial numbers are gotten later
      weights(i)%w(1:3)%n = points(i)%ser
      where ( weights(i)%w(1:3)%n /= 0 ) &
        & weights(i)%w(1:3)%np = QTM_tree%path_vertices(weights(i)%w(1:3)%n)
      weights(i)%w(1:3)%v = points(i)%coeff
      index = indices(i)
      if ( index == 0 ) then
        weights(i)%w(4:maxCoeff)%n = 0  ! Constant vertical extrapolation
                                        !   below range of Heights(:)
      else if ( index >= size(heights,1) ) then
        do j = 1, 3
          weights(i)%w(j)%n = points(i)%ser(j)
          weights(i)%w(j)%nz = size(heights,1)
        end do
        weights(i)%w(4:maxCoeff)%n = 0  ! Constant vertical extrapolation
                                        !   above range of Heights(:)
      else
        dH = heights(index+1) - heights(index)
        if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
          weights(i)%w(4:6)%n = 0 ! Don't use points other than at Height Index
        else
          w = ( heights(index+1) - path_heights(i) ) / dh
          if ( w == 1.0_rg ) then
            weights(i)%w(1:3)%nz = index
            weights(i)%w(4:maxCoeff)%n = 0 ! Don't use points other than at Height Index
          else if ( w == 0.0_rg ) then
            weights(i)%w(j)%nz = index + 1
            weights(i)%w(4:maxCoeff)%n = 0 ! Don't use points other than at Height Index
          else
            weights(i)%w(4:6)%v = weights(i)%w(1:3)%v * ( 1.0_rg - w )
            weights(i)%w(1:3)%v = weights(i)%w(1:3)%v * w
            weights(i)%w(4:6)%n = weights(i)%w(1:3)%n
            weights(i)%w(4:6)%np = weights(i)%w(1:3)%np
            weights(i)%w(1:6)%nz = [ index, index+1, index, index+1, index, index+1 ]
          end if
        end if
      end if

      where ( weights(i)%w%v == 0 ) weights(i)%w%n = 0
      weights(i)%n = count ( weights(i)%w%n /= 0 ) 
      weights(i)%w = pack ( weights(i)%w,weights(i)%w%n /= 0 )
      weights(i)%w(weights(i)%n+1:maxCoeff) = value_t(0,0,0,0)

    end do ! i

  end subroutine QTM_Interpolation_Weights_Geo_3D_List

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

end module QTM_Interpolation_Weights_3D_m

! $Log$
! Revision 2.11  2016/11/12 01:37:02  vsnyder
! Add PV component in case somebody wants to store the path subscript for a
! QTM serial number.  Add Facet argument to routines with scalar Point, in
! case somebody already knows it.  Compute NP and NZ indices newly added
! to Value_t, instead of conflating 2-D subscript into N.
!
! Revision 2.10  2016/11/02 22:55:53  vsnyder
! Use Value_t from Path_Representation_m instead of Weight_t from
! QTM_Interpolation_Weights_m.
!
! Revision 2.9  2016/10/25 18:23:48  vsnyder
! Inching toward 3D-QTM forward model support
!
! Revision 2.8  2016/10/18 00:44:28  vsnyder
! Correct some comments
!
! Revision 2.7  2016/10/05 23:28:22  vsnyder
! Replace ZOT_n component name with Ser because it's a serial number for
! more than just the ZOT coordinates.
!
! Revision 2.6  2016/09/14 20:46:21  vsnyder
! Add H_Ind component, move Weight_1D to Weight_1D_m
!
! Revision 2.5  2016/04/28 23:09:50  vsnyder
! Interpret negative face number as the index of the face intersected by a
! line segment outside the QTM.  Alphabetize components of S_QTM_t.
!
! Revision 2.4  2016/04/16 02:04:06  vsnyder
! Add S_QTM_t, routines for using positions on a line
!
! Revision 2.3  2016/03/30 01:42:02  vsnyder
! Detect and exploit coherent heights in incoherent routine
!
! Revision 2.2  2016/01/26 02:30:32  vsnyder
! Add QTM_Interpolation_Weights_Geo_3D_Incoherent_List to generic
! QTM_Interpolation_Weights.
!
! Revision 2.1  2015/12/31 00:56:47  vsnyder
! Initial commit
!
