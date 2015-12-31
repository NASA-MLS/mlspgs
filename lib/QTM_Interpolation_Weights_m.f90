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
module QTM_Interpolation_Weights_m
!=============================================================================

  use Geolocation_0, only: RG

  implicit NONE
  private

  ! Determine which vertices of QTM to use for interpolation to a specified
  ! point, and calculate their interpolation weights.

  public :: QTM_Interpolation_Weights
  public :: QTM_Interpolation_Weights_Geo, QTM_Interpolation_Weights_Geo_list
  public :: QTM_Interpolation_Weights_Geo_3D, QTM_Interpolation_Weights_Geo_3D_List
  public :: QTM_Interpolation_Weights_ZOT, QTM_Interpolation_Weights_ZOT_list

  interface QTM_Interpolation_Weights
    module procedure &
      & QTM_Interpolation_Weights_Geo, &
      & QTM_Interpolation_Weights_ZOT, &
      & QTM_Interpolation_Weights_Geo_list, &
      & QTM_Interpolation_Weights_Geo_3D, &
      & QTM_Interpolation_Weights_Geo_3D_List, &
      & QTM_Interpolation_Weights_ZOT_list
  end interface QTM_Interpolation_Weights

  public :: Weight_t

  type :: Weight_t
    integer :: Which      ! Position in array element order to which weight is
                          ! germane; zero if not to be used because no facet
                          ! was found (which shouldn't happen) or because the
                          ! Weight is zero.
    real(rg) :: Weight    ! Weight -- interpolation coefficient
  end type Weight_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine QTM_Interpolation_Weights_Geo ( QTM_Tree, Point, Weights, Stack, &
                                           & Used )

    ! Get the interpolation weights for Point in the QTM.  This is the
    ! routine that actually does the work for Point being in some kind of Geo
    ! coordinates.  This is a 2D-only routine; it's used by the 3D routines.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use QTM_m, only: Geo_to_ZOT, Stack_t, ZOT_t
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point
    type(weight_t), intent(out) :: Weights(3)
    type(stack_t), intent(inout), optional :: Stack
    class(h_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(h_t) :: P    ! Coordinates of boundary point nearest to Point.
    type(ZOT_t) :: Z  ! ZOT coordinates of Point.

    f = QTM_tree%find_facet ( point, stack )
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ZOT_n > 0)
    if ( inside ) then
      z = geo_to_ZOT ( point )
      weights%which = QTM_tree%Q(f)%ZOT_n
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & z%x, z%y, weights%weight )
      if ( present(used) ) then
        used%lon = point%lon
        used%lat = point%lat
      end if
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_geo, p )
      z = geo_to_ZOT ( p )
      f = QTM_tree%find_facet ( z, stack )
      ! P should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ZOT_n > 0)
      if ( inside ) then
        weights%which = QTM_tree%Q(f)%ZOT_n
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & z%x, z%y, weights%weight )
      else ! Shouldn't get here
        weights = weight_t(0,0)
      end if
      if ( present(used) ) then
        used%lon = p%lon
        used%lat = p%lat
      end if
    end if

  end subroutine QTM_Interpolation_Weights_Geo

  subroutine QTM_Interpolation_Weights_Geo_list ( QTM_Tree, Point, Weights, Stack )

    ! Get the interpolation weights for Point in the QTM.  This is the
    ! routine that actually does the work for Point being in some kind of Geo
    ! coordinates.  This is a 2D-only routine; it's used by the 3D routines.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t
    use QTM_m, only: Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point(:)          ! Assume size(point) == 
    type(weight_t), intent(out) :: Weights(:,:) ! size(weights,2) and
                                                ! size(weights,1) == 3
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights,2)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(:,i), stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights,2)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(:,i), myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_Geo_list

  subroutine QTM_Interpolation_Weights_Geo_3D ( QTM_Tree, Heights, Point, &
                                              & Weights, N_Weights, Stack, Used )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  The positions are array element order positions,
    ! which can be converted to 2D subscripts (QTM serial number, height)
    ! by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t, H_V_t, RG
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

  subroutine QTM_Interpolation_Weights_Geo_3D_List ( QTM_Tree, Heights, Points, &
                                                   & Weights, N_Weights )

    ! Construct interpolation weights and positions for a point in a 3D grid
    ! for which the horizontal grid is QTM.  The 3D grid is assumed to be
    ! stacked and coherent.  The positions are array element order positions,
    ! which can be converted to 2D subscripts (QTM serial number, height)
    ! by the Subscripts function in the Array_Stuff module.

    use Array_Stuff, only: Element_Position
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_V_t, RG
    use Hunt_m, only: Hunt
    use QTM_m, only: QK

    type(QTM_tree_t), intent(in) :: QTM_Tree
    real(rg), intent(in) :: Heights(:)    ! Assume Heights (meters, zeta) are
    class(h_v_t), intent(in) :: Points(:) !   the same units as for Point%V
    type(weight_t), intent(out) :: Weights(6,size(points))
    integer :: N_Weights                  ! Number of nonzero weights

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
      n_weights = count ( weights(:,i)%weight /= 0 ) 
      where ( weights(:,i)%weight == 0 ) weights(:,i)%which = 0
      weights(:,i) = pack ( weights(:,i),weights(:,i)%which /= 0 )
      weights(n_weights+1:6,i) = weight_t(0,0)
    end do

  end subroutine QTM_Interpolation_Weights_Geo_3D_List

  subroutine QTM_Interpolation_Weights_ZOT ( QTM_Tree, Point, Weights, Stack, &
                                           & Used )

    ! Get the interpolation weights for Point in the QTM.
    ! This is a 2D-only routine.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point
    type(weight_t), intent(out) :: Weights(3)
    type(stack_t), intent(inout), optional :: Stack
    type(ZOT_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(ZOT_t) :: P  ! Coordinates of boundary point nearest to Point.

    f = QTM_tree%find_facet ( point, stack )
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ZOT_n > 0)
    if ( inside ) then
      weights%which = QTM_tree%Q(f)%ZOT_n
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & point%x, point%y, weights%weight )
      if ( present(used) ) used = point
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_zot, p )
      f = QTM_tree%find_facet ( p, stack )
      ! It should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ZOT_n > 0)
      if ( inside ) then
        weights%which = QTM_tree%Q(f)%ZOT_n
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & point%x, point%y, weights%weight )
      else ! Shouldn't get here
        weights = weight_t(0,0)
      end if
      if ( present(used) ) used = p
    end if

  end subroutine QTM_Interpolation_Weights_ZOT

  subroutine QTM_Interpolation_Weights_ZOT_list ( QTM_Tree, Point, Weights, Stack  )

    ! Get the interpolation weights for Point in the QTM.
    ! This is a 2D-only routine.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point(:)         ! Assume size(point) ==
    type(weight_t), intent(out) :: Weights(:,:) ! size(weights,2) and
                                                ! size(weights,1) == 3
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights,2)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(:,i), stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(:,i), myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_ZOT_list

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

end module QTM_Interpolation_Weights_m

! $Log$
! Revision 2.1  2015/12/31 00:05:02  vsnyder
! Initial commit
!
