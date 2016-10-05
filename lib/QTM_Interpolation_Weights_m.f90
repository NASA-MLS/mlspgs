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
  public :: QTM_Interpolation_Weights_ZOT, QTM_Interpolation_Weights_ZOT_list

  interface QTM_Interpolation_Weights
    module procedure &
      & QTM_Interpolation_Weights_Geo, &
      & QTM_Interpolation_Weights_ZOT, &
      & QTM_Interpolation_Weights_Geo_list, &
      & QTM_Interpolation_Weights_ZOT_list
  end interface QTM_Interpolation_Weights

  public :: Weight_t

  type :: Weight_t
    integer :: Which      ! Position in array element order to which weight is
                          ! germane; if no facet was found (which shouldn't
                          ! happen), Which is one and Weight is zero.
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
    ! coordinates.  This is a 2D-only routine.

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
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      z = geo_to_ZOT ( point )
      weights%which = QTM_tree%Q(f)%ser
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
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
      if ( inside ) then
        weights%which = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & z%x, z%y, weights%weight )
      else ! Shouldn't get here
        weights = weight_t(1,0)
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
    ! coordinates.  This is a 2D-only routine.

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
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      weights%which = QTM_tree%Q(f)%ser
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & point%x, point%y, weights%weight )
      if ( present(used) ) used = point
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_zot, p )
      f = QTM_tree%find_facet ( p, stack )
      ! It should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
      if ( inside ) then
        weights%which = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & point%x, point%y, weights%weight )
      else ! Shouldn't get here
        weights = weight_t(1,0)
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
! Revision 2.4  2016/10/05 23:28:22  vsnyder
! Replace ZOT_n component name with Ser because it's a serial number for
! more than just the ZOT coordinates.
!
! Revision 2.3  2016/08/24 22:56:34  vsnyder
! Make Weights%Which == 1 with Weights%Weight == 0 if a facet is not found,
! so Weights can still be used without checking whether a facet was found.
! The result of interpolating with it will be zero.
!
! Revision 2.2  2015/12/31 00:57:32  vsnyder
! Move 3D routines to QTM_Interpolation_Weights_3D_m
!
! Revision 2.1  2015/12/31 00:05:02  vsnyder
! Initial commit
!
