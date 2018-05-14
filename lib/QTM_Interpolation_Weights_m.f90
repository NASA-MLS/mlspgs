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

  use MLSKinds, only: RK => RP

  implicit NONE
  private

  ! Determine which vertices of QTM to use for interpolation to a specified
  ! point, and calculate their horizontal interpolation weights.

  public :: QTM_Interpolation_Weights
  public :: QTM_Interpolation_Weights_Geo, QTM_Interpolation_Weights_Geo_list
  public :: QTM_Interpolation_Weights_ZOT, QTM_Interpolation_Weights_ZOT_list

  interface QTM_Interpolation_Weights
    module procedure &
      & QTM_Interpolation_Weights_Geo, &
      & QTM_Interpolation_Weights_ZOT, &
      & QTM_Interpolation_Weights_Geo_list, &
      & QTM_Interpolation_Weights_ZOT_list, &
      & QTM_Interpolation_Weights_Geo_2d, &
      & QTM_Interpolation_Weights_ZOT_2d, &
      & QTM_Interpolation_Weights_Geo_list_2d, &
      & QTM_Interpolation_Weights_ZOT_list_2d
  end interface QTM_Interpolation_Weights

  public :: Value_QTM_1D_List_t, Value_QTM_2D_List_t, Value_QTM_2D_t

  type :: Value_List ! ( RK ) ! Base type for Value_*List_t
!     integer, kind :: RK
  end type Value_List

  ! For one interpolation weight from a 1D array
  type :: Value_1D_t ! ( RK )
!     integer, kind :: RK
    real(rk) :: V = 0.0    ! Value to be applied at N
    integer :: J = 0       ! Subscript at which to apply V
  end type Value_1D_t

  ! For one interpolation weight from a QTM, or an extract of it onto
  ! an array of profiles adjacent to an integration path, e.g., in Grids_t.
  type, extends(value_1D_t) :: Value_QTM_1D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at N or NP
  ! integer :: J = 0       ! Index in QTM at which to apply V
    integer :: JP = 0      ! Index of path-adjacent part of QTM at which to
                           ! apply V, e.g., in Grids_t.  This is intentionally
                           ! the same name as the second (probably Phi) subscript
                           ! in Value_2D_t so that Interpolate_2D_List.f9h can be
                           ! used for both interpolations.
  end type Value_QTM_1D_t

  ! For interpolating to a list (probably an integration path) from a QTM,
  ! or an extract of it onto an array of profiles adjacent to an integration
  ! path.
  type, extends(Value_List) :: Value_QTM_1D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 3       ! Number of useful elements of V
!     type(value_QTM_1d_t(rk)) :: V(3) = value_QTM_1d_t(rk)()
    type(value_QTM_1d_t) :: V(3) = value_QTM_1d_t()
  end type Value_QTM_1D_List_t

  ! For one interpolation weight from a QTM, or an extract of it onto
  ! an array of profiles adjacent to an integration path, with a zeta basis
  ! at each vertex.
  type, extends(value_QTM_1D_t) :: Value_QTM_2D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at (NZ, N or NP)
  ! integer :: J = 0       ! Index in QTM at which to apply V
  ! integer :: JP = 0      ! Index of path-adjacent part of QTM at which to
                           ! apply V, e.g., in Grids_t.  This is intentionally
                           ! the same name as the second (probably Phi) subscript
                           ! in Value_2D_t so that Interpolate_2D_List.f9h can be
                           ! used for both interpolations.
    integer :: JZ = 0      ! Vertical (probably zeta) index
  end type Value_QTM_2D_t

  ! For interpolating to a list (probably an integration path) from a QTM,
  ! or an extract of it onto an array of profiles adjacent to an integration
  ! path, with a zeta basis at each profile.
  type, extends(Value_List) :: Value_QTM_2D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 6       ! Number of useful elements of V
!     type(value_QTM_2d_t(rk)) :: V(6) = value_QTM_2d_t(rk)()
    type(value_QTM_2d_t) :: V(6) = value_QTM_2d_t()
  end type Value_QTM_2D_List_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine QTM_Interpolation_Weights_Geo ( QTM_Tree, Point, Weights, Facet, &
                                           & Stack, Used )

    ! Get the horizontal interpolation weights for Point in the QTM.  This is
    ! the routine that actually does the work for Point being in some kind of
    ! Geo coordinates.  This is a 2D-only routine.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t! , RG
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use QTM_m, only: Geo_to_ZOT, Stack_t, ZOT_t
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point
!     type(value_QTM_1D_List_t(rg)), intent(out) :: Weights
    type(value_QTM_1D_List_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    class(h_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(h_t) :: P    ! Coordinates of boundary point nearest to Point.
    type(ZOT_t) :: Z  ! ZOT coordinates of Point.

    if ( present(facet) ) then
      f = facet
    else
      f = QTM_tree%find_facet ( point, stack )
    end if
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      z = geo_to_ZOT ( point )
      weights%v%j = QTM_tree%Q(f)%ser
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & z%x, z%y, weights%v%v )
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
        weights%v%j = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & z%x, z%y, weights%v%v )
      else ! Shouldn't get here
!         weights = value_QTM_1D_List_t(rg)(n=1) ! Default weight is zero
        weights = value_QTM_1D_List_t(n=1) ! Default weight is zero
      end if
      if ( present(used) ) then
        used%lon = p%lon
        used%lat = p%lat
      end if
    end if

  end subroutine QTM_Interpolation_Weights_Geo

  subroutine QTM_Interpolation_Weights_Geo_list ( QTM_Tree, Point, Weights, Stack )

    ! Get the horizontal interpolation weights for Point in the QTM.  This is
    ! the routine that actually does the work for Point being in some kind of
    ! Geo coordinates.  This is a 2D-only routine.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t! , RG
    use QTM_m, only: Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point(:)             ! Assume size(point) ==
!     type(value_QTM_1D_list_t(rg)), intent(out) :: Weights(:) ! size(weights), and
    type(value_QTM_1D_list_t), intent(out) :: Weights(:) ! size(weights), and
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_Geo_list

  subroutine QTM_Interpolation_Weights_ZOT ( QTM_Tree, Point, Weights, Facet, &
                                           & Stack, Used )

    ! Get the horizontal interpolation weights for Point in the QTM.
    ! This is a 2D-only routine.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t
!     use Geolocation_0, only: RG
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point
!     type(value_QTM_1D_list_t(rg)), intent(out) :: Weights
    type(value_QTM_1D_list_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    type(ZOT_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(ZOT_t) :: P  ! Coordinates of boundary point nearest to Point.

    if ( present(facet) ) then
      f = facet
    else
      f = QTM_tree%find_facet ( point, stack )
    end if
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      weights%v%j = QTM_tree%Q(f)%ser
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & point%x, point%y, weights%v%v )
      if ( present(used) ) used = point
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_zot, p )
      f = QTM_tree%find_facet ( p, stack )
      ! It should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
      if ( inside ) then
        weights%v%j = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & point%x, point%y, weights%v%v )
      else ! Shouldn't get here
!         weights = value_QTM_1D_List_t(rg)(n=1) ! Default weight is zero
        weights = value_QTM_1D_List_t(n=1) ! Default weight is zero
      end if
      if ( present(used) ) used = p
    end if

  end subroutine QTM_Interpolation_Weights_ZOT

  subroutine QTM_Interpolation_Weights_ZOT_list ( QTM_Tree, Point, Weights, Stack  )

    ! Get the horizontal interpolation weights for Point in the QTM.
    ! This is a 2D-only routine.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t
!     use Geolocation_0, only: RG

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point(:)            ! Assume size(point) ==
!     type(value_QTM_1D_list_t(rg)), intent(out) :: Weights(:) ! size(weights), and
    type(value_QTM_1D_list_t), intent(out) :: Weights(:) ! size(weights), and
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights,2)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_ZOT_list

  subroutine QTM_Interpolation_Weights_Geo_2d ( QTM_Tree, Point, Weights, Facet, &
                                           & Stack, Used )

    ! Get the horizontal interpolation weights for Point in the QTM.  This is
    ! the routine that actually does the work for Point being in some kind of
    ! Geo coordinates.  This is a 2D-only routine, but the Weights argument
    ! is for 3D.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t! , RG
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use QTM_m, only: Geo_to_ZOT, Stack_t, ZOT_t
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point
!     type(value_QTM_2D_List_t(rg)), intent(out) :: Weights
    type(value_QTM_2D_List_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    class(h_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(h_t) :: P    ! Coordinates of boundary point nearest to Point.
    type(ZOT_t) :: Z  ! ZOT coordinates of Point.

    if ( present(facet) ) then
      f = facet
    else
      f = QTM_tree%find_facet ( point, stack )
    end if
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      z = geo_to_ZOT ( point )
      weights%v(1:3)%j = QTM_tree%Q(f)%ser
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & z%x, z%y, weights%v(1:3)%v )
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
        weights%v(1:3)%j = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & z%x, z%y, weights%v(1:3)%v )
      else ! Shouldn't get here
!         weights = value_QTM_2D_List_t(rg)(n=1) ! Default weight is zero
        weights = value_QTM_2D_List_t(n=1) ! Default weight is zero
      end if
      if ( present(used) ) then
        used%lon = p%lon
        used%lat = p%lat
      end if
    end if

  end subroutine QTM_Interpolation_Weights_Geo_2d

  subroutine QTM_Interpolation_Weights_Geo_list_2d ( QTM_Tree, Point, Weights, Stack )

    ! Get the horizontal interpolation weights for Point in the QTM.  This is
    ! the routine that actually does the work for Point being in some kind of
    ! Geo coordinates.  This is a 2D-only routine, but the Weights argument
    ! is for 3D.

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_t! , RG
    use QTM_m, only: Stack_t

    type(QTM_tree_t), intent(in) :: QTM_Tree
    class(h_t), intent(in) :: Point(:)             ! Assume size(point) ==
!     type(value_QTM_2D_List_t(rg)), intent(out) :: Weights(:) ! size(weights), and
    type(value_QTM_2D_List_t), intent(out) :: Weights(:) ! size(weights), and
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_Geo_list_2d

  subroutine QTM_Interpolation_Weights_ZOT_2d ( QTM_Tree, Point, Weights, Facet, &
                                           & Stack, Used )

    ! Get the horizontal interpolation weights for Point in the QTM.
    ! This is a 2D-only routine, but the Weights argument
    ! is for 3D.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t
!     use Geolocation_0, only: RG
    use Nearest_Polygon_Point_m, only: Nearest_Polygon_Point
    use Triangle_Interpolate_m, only: Triangle_Interpolate

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point
!     type(value_QTM_2D_List_t(rg)), intent(out) :: Weights
    type(value_QTM_2D_List_t), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    type(ZOT_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(ZOT_t) :: P  ! Coordinates of boundary point nearest to Point.

    if ( present(facet) ) then
      f = facet
    else
      f = QTM_tree%find_facet ( point, stack )
    end if
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( inside ) then
      weights%v(1:3)%j = QTM_tree%Q(f)%ser
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & point%x, point%y, weights%v(1:3)%v )
      if ( present(used) ) used = point
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_zot, p )
      f = QTM_tree%find_facet ( p, stack )
      ! It should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
      if ( inside ) then
        weights%v(1:3)%j = QTM_tree%Q(f)%ser
        call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                  & point%x, point%y, weights%v(1:3)%v )
      else ! Shouldn't get here
!         weights = value_QTM_2D_List_t(rg)(n=1) ! Default weight is zero
        weights = value_QTM_2D_List_t(n=1) ! Default weight is zero
      end if
      if ( present(used) ) used = p
    end if

  end subroutine QTM_Interpolation_Weights_ZOT_2d

  subroutine QTM_Interpolation_Weights_ZOT_list_2d ( QTM_Tree, Point, Weights, Stack  )

    ! Get the horizontal interpolation weights for Point in the QTM.
    ! This is a 2D-only routine, but the Weights argument
    ! is for 3D.

    use QTM_m, only: Stack_T, ZOT_t
    use Generate_QTM_m, only: QTM_Tree_t
!     use Geolocation_0, only: RG

    type(QTM_tree_t), intent(in) :: QTM_Tree
    type(ZOT_t), intent(in) :: Point(:)            ! Assume size(point) ==
!     type(value_QTM_2D_List_t(rg)), intent(out) :: Weights(:) ! size(weights), and
    type(value_QTM_2D_List_t), intent(out) :: Weights(:) ! size(weights), and
    type(stack_t), intent(inout), optional :: Stack

    integer :: I
    type(stack_t) :: MyStack

    if ( present(stack) ) then
      do i = 1, size(point) ! Assume size(point) == size(weights,2)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=stack )
      end do
    else
      do i = 1, size(point) ! Assume size(point) == size(weights)
        call QTM_Interpolation_Weights ( QTM_Tree, point(i), weights(i), &
          & stack=myStack )
      end do
    end if

  end subroutine QTM_Interpolation_Weights_ZOT_list_2d

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
! Revision 2.11  2018/05/14 23:25:29  vsnyder
! Change to sparse eta representation
!
! Revision 2.10  2017/08/28 20:27:47  livesey
! Changed the n,nf,np,nz elements to j,jf,...
!
! Revision 2.9  2016/11/23 20:08:45  vsnyder
! Comment out Use...RG until processors support KPDT
!
! Revision 2.8  2016/11/23 00:09:36  vsnyder
! Use types from Indexed_Values_m.  Add _2d routines to do horizontal
! interpolation for later 3D interpolation.
!
! Revision 2.7  2016/11/12 01:34:06  vsnyder
! Add Facet argument for routines with scalar Point, in case we know it
!
! Revision 2.6  2016/11/02 22:55:15  vsnyder
! Remove Weight_t in favor of Value_t in Path_Representation_m
!
! Revision 2.5  2016/10/25 18:23:48  vsnyder
! Inching toward 3D-QTM forward model support
!
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
! Move 3D routines to QTM_Interpolation_Weights_2d_m
!
! Revision 2.1  2015/12/31 00:05:02  vsnyder
! Initial commit
!
