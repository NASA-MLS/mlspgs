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
  use QTM_Interpolation_Weights_m, only: QTM_Interpolation_Weights, Weight_t
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
                     ! Inside_Prism => The intersection is not with a face of a
                     !      prism resting on the surface QTM; it is inside a
                     !      prism resting on the surface QTM.
                     ! If the point is outside the QTM, this will be the
                     ! negative of one of the above values, other than
                     ! Inside_Prism; it will be whatever face the extrapolated
                     ! line intersects.  The coordinates indexed by ZOT_N are
                     ! those where the line intersects the QTM, not where it
                     ! intersects a surface of a height equal to one of the
                     ! heights in the intersection of all heights in the state
                     ! vector.  The H component probably isn't useful.
    integer :: Facet = 0 ! Index in Grid%QTM_Tree of QTM facet intersected at
                     ! the point.
    real(rg) :: H    ! Height at which the point intersects a face.  This is
                     ! the height where the extrapolated line intersects a face
                     ! of the QTM for points outside the QTM.
    integer :: H_ind = 0 ! Index in height grid.  If the intersection is not
                     ! with a horizontal boundary within the QTM, H_Ind is the
                     ! index of the next lower height surface.
    integer :: N_Coeff = 0  ! How many elements of Coeff and ZOT_N are used.
    integer :: ZOT_N(3) = 0 ! Serial numbers of ZOT coordinates, see QTM_Node_t
                     ! in the Generate_QTM_m module.
  end type S_QTM_t

  integer, parameter, public :: Cone_Face = 1
  integer, parameter, public :: X_face = cone_face + 1
  integer, parameter, public :: Y_face = x_face + 1
  integer, parameter, public :: Top_Face = y_face + 1
  integer, parameter, public :: Inside_Prism = top_face + 1

  public :: QTM_Interpolation_Weights, Weight_t

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
                                              & Weights, N_Weights, Stack, Used )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  The positions are array element order positions,
    ! which can be converted to 2D subscripts (QTM serial number, height)
    ! by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t, H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK, Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)   ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Point    !   the same units as for Point%V
    type(weight_t), intent(out) :: Weights(6)
    integer :: N_Weights                 ! Number of nonzero weights
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
    call QTM_Interpolation_Weights ( QTM_Tree, Point, Weights(1:3), Stack, Used )

    ! Get vertical interpolation coefficients
    call hunt ( heights, point%v, index, &
              & allowTopValue=.true., allowBelowValue = .true. )

    ! Combine horizontal and vertical interpolation coefficients
    if ( index == 0 ) then
      weights(4:6)%which = 0 ! Constant vertical extrapolation outside
                             !   range of Heights(:)
    else if ( index >= size(heights) ) then
      do j = 1, 3
        weights(j)%which = element_position ( [ weights(j)%which, size(heights) ], &
                                            & [ n_QTM, size(heights) ] )
      end do
      weights(4:6)%which = 0 ! Constant vertical extrapolation outside
                             !   range of Heights(:)
    else
      dH = heights(index+1) - heights(index)
      if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
        weights(4:6)%which = 0 ! Don't use points other than at Height Index
      else
        w = ( heights(index+1) - point%v ) / dh
        if ( w == 1.0_rg ) then
          do j = 1, 3
            weights(j)%which = element_position ( [ weights(j)%which, index ], &
                                                & [ n_QTM, size(heights) ] )
          end do
          weights(4:6)%which = 0 ! Don't use points other than at Height Index
        else if ( w == 0.0_rg ) then
          do j = 1, 3
            weights(j)%which = element_position ( [ weights(j)%which, index+1 ], &
                                                & [ n_QTM, size(heights) ] )
          end do
          weights(4:6)%which = 0 ! Don't use points other than at Height Index + 1
        else
          weights(4:6)%weight = weights(1:3)%weight * ( 1.0_rg - w )
          weights(1:3)%weight = weights(1:3)%weight * w
          do j = 1, 3
            weights(j+3)%which = element_position ( [ weights(j)%which, index+1 ], &
                                                  & [ n_QTM, size(heights) ] )
            weights(j)%which = element_position ( [ weights(j)%which, index ], &
                                                & [ n_QTM, size(heights) ] )
          end do
        end if
      end if
    end if
    n_weights = count ( weights%weight /= 0 ) 
    where ( weights%weight == 0 ) weights%which = 0
    weights = pack ( weights,weights%which /= 0 )
    weights(n_weights+1:6) = weight_t(0,0)

  end subroutine QTM_Interpolation_Weights_Geo_3D

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent ( QTM_Tree, Heights, &
                                     & Point, Weights, N_Weights, Stack, Used )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  The positions are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t, H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK, Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! Extents are (heights, QTM_Tree%N_In)
                                         ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Point    !   the same units as for Point%V
    type(weight_t), intent(out) :: Weights(6)
    integer :: N_Weights                 ! Number of nonzero weights
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
                                & Point, Weights, N_Weights, Stack, Used )
      return
    end if

    n_QTM = QTM_Tree%n_in

    ! Get horizontal interpolation coefficients and serial numbers
    call QTM_Interpolation_Weights ( QTM_Tree, Point, Weights(1:3), Stack, Used )

    ! Get vertical interpolation coefficients
    do p = 1, 3
      index = 0
      s = weights(p)%which
      if ( s /= 0 ) then
        call hunt ( heights(:,s), point%v, index, &
                  & allowTopValue=.true., allowBelowValue = .true. )
      end if

      ! Combine horizontal and vertical interpolation coefficients
      if ( index == 0 ) then
        weights(4:6)%which = 0 ! Constant vertical extrapolation outside
                               !   range of Heights(:)
      else if ( index >= size(heights) ) then
        weights(p)%which = element_position ( [ weights(p)%which, size(heights,1) ], &
          & [ n_QTM, size(heights,1) ] )
        weights(4:6)%which = 0 ! Constant vertical extrapolation outside
                               !   range of Heights(:)
      else
        dH = heights(index+1,s) - heights(index,s)
        if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
          weights(4:6)%which = 0 ! Don't use points other than at Height Index
        else
          w = ( heights(index+1,p) - point%v ) / dh
          if ( w == 1.0_rg ) then
            weights(p)%which = element_position ( [ weights(p)%which, index ], &
                                                & [ n_QTM, size(heights,1) ] )
            weights(4:6)%which = 0 ! Don't use points other than at Height Index
          else if ( w == 0.0_rg ) then
            weights(p)%which = element_position ( [ weights(p)%which, index+1 ], &
                                                & [ n_QTM, size(heights,1) ] )
            weights(4:6)%which = 0 ! Don't use points other than at Height Index + 1
          else
            weights(4:6)%weight = weights(1:3)%weight * ( 1.0_rg - w )
            weights(1:3)%weight = weights(1:3)%weight * w
            weights(p+3)%which = element_position ( [ weights(p)%which, index+1 ], &
                                                  & [ n_QTM, size(heights,1) ] )
            weights(p)%which = element_position ( [ weights(p)%which, index ], &
                                                & [ n_QTM, size(heights,1) ] )
          end if
        end if
      end if

    end do ! p

    n_weights = count ( weights%weight /= 0 ) 
    where ( weights%weight == 0 ) weights%which = 0
    weights = pack ( weights,weights%which /= 0 )
    weights(n_weights+1:6) = weight_t(0,0)

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line ( QTM_Tree, &
                         & Heights, Line, Points, Weights, N_Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  The positions are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.
    ! The points are described by Line, which consists of a vector in ECR to
    ! a point on a line and an unit vector along the line, and Points, which
    ! describe distances along that line and how the line intersects a prism
    ! of a QTM-based grid (or not).

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_V_Geod
    use Pure_Hunt_m, only: PureHunt

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! Extents are (heights, QTM_Tree%N_In)
                                         ! Assume Heights and QTM_Tree came from
                                         ! the same Geolocation_t structure
                                         ! (see Geolocation_m).
    type(ECR_t), intent(in) :: Line(2)   ! Vector to a point on the line,
                                         ! vector along the line.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(weight_t), intent(out) :: Weights(6,size(points))
    integer, intent(out) :: N_Weights(size(points)) ! Number of nonzero weights

    type(H_V_Geod) :: G                  ! Geodetic coordinate of Points(i)
    integer :: I, J
    integer :: Jlo, Jhi                  ! Indices of H[pxy] in Heights
    integer :: Nh                        ! Number of heights
    integer :: NQ                        ! Number of QTM vertices
    type(ECR_t) :: P                     ! ECR coordinate of Points(i)
    real(rg) :: W(2,2)                   ! Interpolation weights

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, Heights(:,1), &
                                     & Line, Points, Weights, N_Weights )
      return
    end if

    nh = size(heights,1)
    nq = size(QTM_Tree%q)

    do i = 1, size(points)
      p = line(1) + points(i)*line(2)
      g = p%geod()
      select case ( abs(points(i)%face) )
      case ( cone_face, X_face, Y_face )
        ! Interpolate in height at the X and Y, X and P, or Y and P nodes
        ! of the QTM.
        jlo = 1; jhi = nh
        call weight_1d ( heights(:,points(i)%ZOT_n(1)), g%v, jlo, jhi, w(:,1) )
        call weight_1d ( heights(:,points(i)%ZOT_n(2)), g%v, jlo, jhi, w(:,2) )
        ! Assume that it's extremely rare to hit a vertex, so compute four
        ! weights.
        n_weights(i) = 4
        weights(1:4,i)%weight = [ w(:,1)*points(i)%coeff(1), &
                                & w(:,2)*points(i)%coeff(2) ]
        weights(5:6,i)%weight = 0
        weights(1:4,i)%which = &
          & [ element_position([jlo,points(i)%ZOT_n(1)], [nh,nq]), &
            & element_position([jhi,points(i)%ZOT_n(1)], [nh,nq]), &
            & element_position([jlo,points(i)%ZOT_n(2)], [nh,nq]), &
            & element_position([jhi,points(i)%ZOT_n(2)], [nh,nq]) ]
        weights(5:6,i)%which = 0
      case ( inside_prism ) ! Use the general method
        call QTM_interpolation_weights ( QTM_tree, heights, g, weights(:,i), &
                                       & n_weights(i) )
      case ( top_face )
        ! Don't need height interpolation, but we do need to find the index
        ! in the heights array, because the heights array here might be a
        ! subset of the heights array given to Metrics_3D.
        n_weights(i) = 3
        weights(1:3,i)%weight = points(i)%coeff(1:3)
        weights(4:6,i)%weight = 0
        do j = 1, 3
          call purehunt ( points(i)%h, heights(:,points(i)%ZOT_n(j)), &
                        & nh, jlo, jhi )
          weights(j,i)%which = &
            & element_position([jlo,points(i)%ZOT_n(j)], [nh,nq])
        end do
        weights(4:6,i)%which = 0
      end select
    end do

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_Line

  subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_List ( QTM_Tree, &
                         & Heights, Points, Weights, N_Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked but not necessarily coherent.  The positions are array element
    ! order positions, which can be converted to 2D subscripts (QTM serial
    ! number, height) by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:,:) ! Extents are (heights, QTM_Tree%N_In)
                                         ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Points(:) !   the same units as for Points%V
    type(weight_t), intent(out) :: Weights(6,size(points))
    integer :: N_Weights(size(points))   ! Number of nonzero weights

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: I, Index
    integer(qk) :: N_QTM ! Number of vertices of QTM within or adjacent to
                         ! the polygon
    integer :: P         ! Index of a corner of a facet
    integer :: S         ! Serial number of a point in the QTM
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    if ( size(heights,2) == 1 ) then ! Coherent heights
      call QTM_Interpolation_Weights ( QTM_Tree, Heights(:,1), &
                                & Points, Weights, N_Weights )
      return
    end if

    n_QTM = QTM_Tree%n_in

    ! Get horizontal interpolation coefficients and serial numbers
    call QTM_Interpolation_Weights ( QTM_Tree, Points, Weights(1:3,:) )

    ! Get vertical interpolation coefficients
    do i = 1, size(points)
      do p = 1, 3
        index = 0
        s = weights(p,i)%which
        if ( s /= 0 ) then
          call hunt ( heights(:,s), points(i)%v, index, &
                    & allowTopValue=.true., allowBelowValue = .true. )
        end if

        ! Combine horizontal and vertical interpolation coefficients
        if ( index == 0 ) then
          weights(4:6,i)%which = 0 ! Constant vertical extrapolation outside
                                   !   range of Heights(:)
        else if ( index >= size(heights) ) then
          weights(p,i)%which = element_position ( [ weights(p,i)%which, size(heights,1) ], &
                                                & [ n_QTM, size(heights,1) ] )
          weights(4:6,i)%which = 0 ! Constant vertical extrapolation outside
                                   !   range of Heights(:)
        else
          dH = heights(index+1,s) - heights(index,s)
          if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
            weights(4:6,i)%which = 0 ! Don't use points other than at Height Index
          else
            w = ( heights(index+1,p) - points(i)%v ) / dh
            if ( w == 1.0_rg ) then
              weights(p,i)%which = element_position ( [ weights(p,i)%which, index ], &
                                                    & [ n_QTM, size(heights,1) ] )
              weights(4:6,i)%which = 0 ! Don't use points other than at Height Index
            else if ( w == 0.0_rg ) then
              weights(p,i)%which = element_position ( [ weights(p,i)%which, index+1 ], &
                                                    & [ n_QTM, size(heights,1) ] )
              weights(4:6,i)%which = 0 ! Don't use points other than at Height Index + 1
            else
              weights(4:6,i)%weight = weights(1:3,i)%weight * ( 1.0_rg - w )
              weights(1:3,i)%weight = weights(1:3,i)%weight * w
              weights(p+3,i)%which = element_position ( [ weights(p,i)%which, index+1 ], &
                                                      & [ n_QTM, size(heights,1) ] )
              weights(p,i)%which = element_position ( [ weights(p,i)%which, index ], &
                                                    & [ n_QTM, size(heights,1) ] )
            end if
          end if
        end if

      end do ! p

      n_weights(i) = count ( weights(:,i)%weight /= 0 ) 
      where ( weights(:,i)%weight == 0 ) weights(:,i)%which = 0
      weights(:,i) = pack ( weights(:,i),weights(:,i)%which /= 0 )
      weights(n_weights(i)+1:6,i) = weight_t(0,0)
    end do ! i

  end subroutine QTM_Interpolation_Weights_Geo_3D_Incoherent_List

  subroutine QTM_Interpolation_Weights_Geo_3D_Line ( QTM_Tree, &
                         & Heights, Line, Points, Weights, N_Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  The positions are array element order positions,
    ! which can be converted to 2D subscripts (QTM serial number, height) by
    ! the Subscripts function in the Array_Stuff module. The points are
    ! described by Line, which consists of a vector in ECR to a point on a
    ! line and an unit vector along the line, and Points, which describe
    ! distances along that line and how the line intersects a prism of a
    ! QTM-based grid (or not).

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: ECR_t, H_V_Geod
    use Hunt_m, only: PureHunt

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)   ! Assume Heights and QTM_Tree came from
                                         ! the same Geolocation_t structure
                                         ! (see Geolocation_m).
    type(ECR_t), intent(in) :: Line(2)   ! Vector to a point on the line,
                                         ! vector along the line.
    type(S_QTM_t), intent(in) :: Points(:) ! Points along Line
    type(weight_t), intent(out) :: Weights(6,size(points))
    integer, intent(out) :: N_Weights(size(points)) ! Number of nonzero weights

    type(h_v_geod) :: G                  ! Geodetic coordinate of Points(i)
    integer :: I, J
    integer :: Jlo, Jhi                  ! Indices of H[pxy] in Heights
    integer :: Nh                        ! Number of heights
    integer :: NQ                        ! Number of QTM vertices
    type(ECR_t) :: P                     ! ECR coordinate of Points(i)
    real(rg) :: W(2)                     ! Interpolation weights

    nh = size(heights,1)
    nq = size(QTM_Tree%q)

    do i = 1, size(points)
      p = line(1) + points(i)*line(2)
      g = p%geod()
      select case ( abs(points(i)%face) )
      case ( cone_face, X_face, Y_face )
        ! Interpolate in height at the X and Y, X and P, or Y and P nodes
        ! of the QTM.
        jlo = 1; jhi = nh
        call weight_1d ( heights, g%v, jlo, jhi, w )
        ! Assume that it's extremely rare to hit a vertex, so compute four
        ! weights.
        n_weights(i) = 4
        weights(1:4,i)%weight = [ w*points(i)%coeff(1), &
                                & w*points(i)%coeff(2) ]
        weights(5:6,i)%weight = 0
        weights(1:4,i)%which = &
          & [ element_position([jlo,points(i)%ZOT_n(1)], [nh,nq]), &
            & element_position([jhi,points(i)%ZOT_n(1)], [nh,nq]), &
            & element_position([jlo,points(i)%ZOT_n(2)], [nh,nq]), &
            & element_position([jhi,points(i)%ZOT_n(2)], [nh,nq]) ]
        weights(5:6,i)%which = 0
      case ( inside_prism ) ! Use the general method
        call QTM_interpolation_weights ( QTM_tree, heights, g, weights(:,i), &
                                       & n_weights(i) )
      case ( top_face )
        ! Don't need height interpolation, but we do need to find the index
        ! in the heights array, because the heights array here might be a
        ! subset of the heights array given to Metrics_3D.
        n_weights(i) = 3
        weights(1:3,i)%weight = points(i)%coeff(1:3)
        weights(4:6,i)%weight = 0
        do j = 1, 3
          call purehunt ( points(i)%h, heights, nh, jlo, jhi )
          weights(j,i)%which = &
            & element_position([jlo,points(i)%ZOT_n(j)], [nh,nq])
        end do
        weights(4:6,i)%which = 0
      end select
    end do

  end subroutine QTM_Interpolation_Weights_Geo_3D_Line

  subroutine QTM_Interpolation_Weights_Geo_3D_List ( QTM_Tree, Heights, Points, &
                                                   & Weights, N_Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  The positions are array element order positions,
    ! which can be converted to 2D subscripts (QTM serial number, height)
    ! by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_V_t
    use Hunt_m, only: Hunt
    use QTM_m, only: QK

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)    ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Points(:) !   the same units as for Points%V
    type(weight_t), intent(out) :: Weights(6,size(points))
    integer :: N_Weights(size(points))    ! Number of nonzero weights

    real(rg) :: dH       ! Difference in Heights coordinates
    integer :: I, Index, Indices(size(points)), J
    integer(qk) :: N_QTM ! Number of vertices of QTM within or adjacent to
                         ! the polygon
    real(rg) :: W        ! Weight in the vertical direction for the lower level

    n_QTM = QTM_Tree%n_in

    ! Get horizontal interpolation coefficients and serial numbers
    call QTM_Interpolation_Weights ( QTM_Tree, Points, Weights(1:3,:) )

    ! Get vertical interpolation coefficients
    call hunt ( heights, points%v, indices, &
              & allowTopValue=.true., allowBelowValue = .true. )

    ! Combine horizontal and vertical interpolation coefficients
    do i = 1, size(points)
      index = indices(i)
      if ( index == 0 ) then
        weights(4:6,i)%which = 0 ! Constant vertical extrapolation outside
                               !   range of Heights(:)
      else if ( index >= size(heights) ) then
        do j = 1, 3
          weights(j,i)%which = element_position ( [ weights(j,i)%which, size(heights) ], &
            & [ n_QTM, size(heights) ] )
        end do
        weights(4:6,i)%which = 0 ! Constant vertical extrapolation outside
                               !   range of Heights(:)
      else
        dH = heights(index+1) - heights(index)
        if ( dH == 0 ) then ! shouldn't happen if Heights(:) are distinct
          weights(4:6,i)%which = 0 ! Don't use points other than at Height Index
        else
          w = ( heights(index+1) - points(i)%v ) / dh
          if ( w == 1.0_rg ) then
            do j = 1, 3
              weights(j,i)%which = element_position ( [ weights(j,i)%which, index ], &
                                                    & [ n_QTM, size(heights) ] )
            end do
            weights(4:6,i)%which = 0 ! Don't use points other than at Height Index
          else if ( w == 0.0_rg ) then
            do j = 1, 3
              weights(j,i)%which = element_position ( [ weights(j,i)%which, index+1 ], &
                                                    & [ n_QTM, size(heights) ] )
            end do
            weights(4:6,i)%which = 0 ! Don't use points other than at Height Index + 1
          else
            weights(4:6,i)%weight = weights(1:3,i)%weight * ( 1.0_rg - w )
            weights(1:3,i)%weight = weights(1:3,i)%weight * w
            do j = 1, 3
              weights(j,i+3)%which = element_position ( [ weights(j,i)%which, index+1 ], &
                                                      & [ n_QTM, size(heights) ] )
              weights(j,i)%which = element_position ( [ weights(j,i)%which, index ], &
                                                    & [ n_QTM, size(heights) ] )
            end do
          end if
        end if
      end if
      n_weights(i) = count ( weights(:,i)%weight /= 0 ) 
      where ( weights(:,i)%weight == 0 ) weights(:,i)%which = 0
      weights(:,i) = pack ( weights(:,i),weights(:,i)%which /= 0 )
      weights(n_weights(i)+1:6,i) = weight_t(0,0)
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
