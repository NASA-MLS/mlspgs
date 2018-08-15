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
  use Sparse_m, only: Sparse_Element_t

  implicit NONE
  private

  ! Determine which vertices of QTM to use for interpolation to a specified
  ! point, and calculate their horizontal interpolation weights.

  public :: Fill_Weights
  public :: QTM_Interpolation_Weights
  public :: QTM_Interpolation_Weights_Geo
  public :: QTM_Interpolation_Weights_ZOT
  public :: QTM_Weights_t

  interface Fill_Weights
    module procedure Fill_Weights_List
    module procedure Fill_Weights_Sparse
  end interface Fill_Weights

  interface QTM_Interpolation_Weights
    module procedure QTM_Interpolation_Weights_Geo
    module procedure QTM_Interpolation_Weights_ZOT
  end interface QTM_Interpolation_Weights

  public :: Value_QTM_1D_List_t, Value_QTM_2D_List_t, Value_QTM_2D_t

  type :: QTM_Weights_t
    integer :: N = 3               ! How many weights
    ! The weights, 3 for 2D, 6 for 3D:
    type(sparse_element_t) :: W(6) = sparse_element_t(0.0_rk,0,0,0,0)
  contains
    procedure :: Fill_Weights_Sparse
    generic :: Fill_Weights => Fill_Weights_Sparse
    procedure :: Copy => Copy_Weights_Sparse
!????? Crashes Intel 17 17.0.0.098 Build 20160721
!????? where call x%copy ( y ) or x = y appears.
!     generic :: assignment ( = ) => Copy
  end type QTM_Weights_t

  type :: Value_List ! ( RK ) ! Base type for Value_*List_t
!     integer, kind :: RK
  contains
    procedure :: Fill_Weights_List
    generic :: Fill_Weights => Fill_Weights_List
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

  subroutine Fill_Weights_List ( Weights, Ser, V )
    ! Fill the Value_QTM_[12]D_List_t structures using given QTM serial numbers
    ! Ser and weight values V.  This is only a 2D routine, but it allows the
    ! dynamic type of Weights to be for a 3D interpolation.
    class(value_list), intent(out) :: Weights
    integer, intent(in) :: Ser(:)  ! Up to six serial numbers
    real(rk), intent(in) :: V(:)   ! One weight for each Ser
    integer :: N
    n = size(ser)
    select type ( weights )
    type is (value_QTM_1D_List_t)
      weights%n = 3
      weights%v(:n)%j = ser
      weights%v(:n)%v = v
    type is (value_QTM_2D_List_t)
      weights%n = 3
      weights%v(:n)%j = ser
      weights%v(:n)%v = v
    end select
  end subroutine Fill_Weights_List

  subroutine Copy_Weights_Sparse ( To_Weights, From_Weights )
    ! Only copy From_Weights%w(1:From_Weights%n)
    class( QTM_Weights_t ), intent(inout) :: To_Weights
    class( QTM_Weights_t ), intent(in) :: From_Weights
    to_weights%n = from_weights%n
    to_weights%w(1:to_weights%n) = from_weights%w(1:from_weights%n)
  end subroutine Copy_Weights_Sparse

  subroutine Fill_Weights_Sparse ( Weights, Ser, V )
    ! Fill the QTM_Weights_t structure using given QTM serial numbers Ser and
    ! weight values V.
    class( QTM_Weights_t ), intent(out) :: Weights
    integer, intent(in) :: Ser(:)  ! Up to six serial numbers
    real(rk), intent(in) :: V(:)   ! One weight for each Ser
    integer :: N
    n = size(ser)
    weights%n = n
    weights%w(:n)%r = 0            ! Row number (path position) not known yet
    weights%w(:n)%c = ser          ! Column numbers are QTM serial numbers
    weights%w(:n)%v = v
  end subroutine Fill_Weights_Sparse

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
    class(value_list), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    real(rk) :: V(3)
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
    else
      call nearest_polygon_point ( point, QTM_tree%polygon_geo, p )
      z = geo_to_ZOT ( p )
      f = QTM_tree%find_facet ( z, stack )
      ! P should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    end if
    if ( inside ) then
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & z%x, z%y, v )
      call weights%fill_weights ( QTM_tree%Q(f)%ser, v )
    else ! Shouldn't get here
      call weights%fill_weights ( [0], [0.0_rk] )
    end if
    if ( present(used) ) then
      used%lon = p%lon
      used%lat = p%lat
    end if

  end subroutine QTM_Interpolation_Weights_Geo

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
    class(value_list), intent(out) :: Weights
    integer, intent(in), optional :: Facet
    type(stack_t), intent(inout), optional :: Stack
    type(ZOT_t), intent(out), optional :: Used ! The point used for interpolation

    logical :: Inside ! Point is inside a facet with three serial numbers.
    integer :: F      ! Index in QTM_Tree of facet containing Point, if any.
    type(ZOT_t) :: P  ! Coordinates of boundary point nearest to Point.
    real(rk) :: V(3)

    if ( present(facet) ) then
      f = facet
    else
      f = QTM_tree%find_facet ( point, stack )
    end if
    inside = f > 0
    if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    if ( .not. inside ) then
      call nearest_polygon_point ( point, QTM_tree%polygon_zot, p )
      f = QTM_tree%find_facet ( p, stack )
      ! It should be inside!
      inside = f > 0
      if ( inside ) inside = all(QTM_tree%Q(f)%ser > 0)
    end if
    if ( inside ) then
      call triangle_interpolate ( QTM_tree%Q(f)%z%x, QTM_tree%Q(f)%z%y, &
                                & point%x, point%y, v )
      call weights%fill_weights ( QTM_tree%Q(f)%ser, v )
    else ! Shouldn't get here
      call weights%fill_weights ( [0], [0.0_rk] )
    end if
    if ( present(used) ) used = p

  end subroutine QTM_Interpolation_Weights_ZOT

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
! Revision 2.13  2018/08/15 17:48:17  vsnyder
! Repair botched comment from previous check in:
! Remove list interpolators.  Add QTM_Weights_t and Fill_Weights.
!
! Revision 2.12  2018/08/15 17:45:14  vsnyder
! Remove list interpolators.  Add QTM_Weights_t and Fill_Weights.
!
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
