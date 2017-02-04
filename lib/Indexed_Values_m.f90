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
module Indexed_Values_m
!=============================================================================

  ! Types to represent indexed values.  These provide a subscript or
  ! subscripts, and a value that applies at a position in an array
  ! indexed by that (those) subscript(s).

  ! These are usually used for interpolation coefficients.

  ! Using RP for RK is temporary until all processors support parameterized
  ! derived types.  There is an implicit assumption that RP == RG (from
  ! Geolocation_0).

  use MLSKinds, only: RK => RP

  implicit NONE
  public

  ! Inheritance heirarchy:

  ! Value_1D_t -> Value_QTM_1D_t -> Value_QTM_2D_t -> Value_QTM_3D_t
  !            -> Value_2D_t -> Value_3D_t
  ! Value_List (abstract) ->
  !            -> Value_1D_List_t
  !            -> Value_1D_p_t
  !            -> Value_QTM_1D_List_t
  !            -> Value_2D_List_t
  !            -> Value_3D_List_t
  !            -> Value_QTM_2D_List_t
  !            -> Value_QTM_3D_List_t
  ! Value_Lists (abstract) ->
  !            -> Value_1D_Lists_t
  !            -> Value_1D_Lists_t
  !            -> Value_1D_Lists_t
  !            -> Value_2D_Lists_t
  !            -> Value_3D_Lists_t
  !            -> Value_QTM_1D_Lists_t
  !            -> Value_QTM_2D_Lists_t
  !            -> Value_QTM_3D_Lists_t

  type :: Value_List ! ( RK ) ! Base type for Value_*List_t
!     integer, kind :: RK
  end type Value_List

  ! For one interpolation weight from a 1D array
  type :: Value_1D_t ! ( RK )
!     integer, kind :: RK
    real(rk) :: V = 0.0    ! Value to be applied at N
    integer :: N = 0       ! Subscript at which to apply V
  end type Value_1D_t

  ! For one interpolation weight from a 1D array to a 1D array
  type, extends(Value_List) :: Value_1D_p_t ! ( RK )
!     integer, kind :: RK
    real(rk) :: V = 0.0    ! Value to be applied at N
    integer :: N = 0       ! Subscript at which to apply V
    integer :: P = 0       ! Subscript at which interpolated value accumulates
  end type Value_1D_p_t

  ! For interpolating to a list from a 1D array, assuming linear interolation
  type, extends(Value_List) :: Value_1D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 2       ! Number of useful elements of V
!     type(value_1d_t(rk)) :: V(2) = value_1d_t(rk)()
    type(value_1d_t) :: V(2) = value_1d_t()
  end type Value_1D_List_t

  ! For one interpolation weight from a QTM, or an extract of it onto
  ! an array of profiles adjacent to an integration path.
  type, extends(value_1D_t) :: Value_QTM_1D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at N or NP
  ! integer :: N = 0       ! Index in QTM at which to apply V
    integer :: NP = 0      ! Index of path-adjacent part of QTM at which to
                           ! apply V
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

  ! For one interpolation weight from a two-dimensional rectangular array, such
  ! as a state vector.
  type, extends(value_1D_t) :: Value_2D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at N or NP
  ! integer :: N = 0       ! First (probably zeta) subscript at which to apply V
    integer :: NP = 0      ! Second (probably phi) subscript at which to apply V
  end type Value_2D_t

  ! For interpolating to a list (probably an integration path) from a two-
  ! dimensional rectangular array, such as a state vector, assuming bilinear
  ! interpolation.
  type, extends(Value_List) :: Value_2D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 4       ! Number of useful elements of V
!     type(value_2d_t(rk)) :: V(4) = value_2d_t(rk)()
    type(value_2d_t) :: V(4) = value_2d_t()
  end type Value_2D_List_t

  ! For one interpolation weight from a three-dimensional rectangular array,
  ! such as a state vector.
  type, extends(value_2D_t) :: Value_3D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at N or NP
  ! integer :: N = 0       ! Second (probably zeta) subscript at which to apply V
  ! integer :: NP = 0      ! Third (probably phi) subscript at which to apply V
    integer :: NF = 0      ! First (probaly frequency) subscript at which to apply V
  end type Value_3D_t

  ! For interpolating to a list (probably an integration path) from a three-
  ! dimensional rectangular array, such as a state vector, assuming trilinear
  ! interpolation.
  type, extends(Value_List) :: Value_3D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 8       ! Number of useful elements of V
!     type(value_3d_t(rk)) :: V(8) = value_3d_t(rk)()
    type(value_3d_t) :: V(8) = value_3d_t()
  end type Value_3D_List_t

  ! For one interpolation weight from a QTM, or an extract of it onto
  ! an array of profiles adjacent to an integration path, with a zeta basis
  ! at each vertex.
  type, extends(value_QTM_1D_t) :: Value_QTM_2D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at (NZ, N or NP)
  ! integer :: N = 0       ! Index in QTM at which to apply V
  ! integer :: NP = 0      ! Index of path-adjacent part of QTM at which to
                           ! apply V
    integer :: NZ = 0      ! Vertical (probably zeta) index
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

  ! For one interpolation weight from a QTM, or an extract of it onto
  ! an array of profiles adjacent to an integration path, with a zeta basis
  ! and a frequency basis at each vertex.
  type, extends(value_QTM_2D_t) :: Value_QTM_3D_t
  ! real(rk) :: V = 0.0    ! Value to be applied at (NF, NZ, N or NP)
  ! integer :: N = 0       ! Index in QTM at which to apply V
  ! integer :: NP = 0      ! Index of path-adjacent part of QTM at which to
                           ! apply V
  ! integer :: NZ = 0      ! Vertical (probably zeta) index
    integer :: NF = 0      ! Frequency subscript
  end type Value_QTM_3D_t

  ! For interpolating to a list (probably an integration path) from a QTM,
  ! or an extract of it onto an array of profiles adjacent to an integration
  ! path, with a zeta basis and frequency basis at each profile.
  type, extends(Value_List) :: Value_QTM_3D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 12      ! Number of useful elements of V
!     type(value_QTM_3d_t(rk)) :: V(12) = value_QTM_3d_t(rk)()
    type(value_QTM_3d_t) :: V(12) = value_QTM_3d_t()
  end type Value_QTM_3D_List_t

  ! Use this type to create an array of lists, potentially of different sizes.
  type, abstract :: Value_Lists
    integer :: N ! Number of elements of Eta currently in use.
  end type Value_Lists

  type, extends(value_Lists) :: Value_1D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_1D_list_t), allocatable :: Eta(:)
  end type Value_1D_Lists_t

  type, extends(value_Lists) :: Value_2D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_2D_list_t), allocatable :: Eta(:)
  end type Value_2D_Lists_t

  type, extends(value_Lists) :: Value_3D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_3D_list_t), allocatable :: Eta(:)
  end type Value_3D_Lists_t

  type, extends(value_Lists) :: Value_QTM_1D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_QTM_1D_list_t), allocatable :: Eta(:)
  end type Value_QTM_1D_Lists_t

  type, extends(value_Lists) :: Value_QTM_2D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_QTM_2D_list_t), allocatable :: Eta(:)
  end type Value_QTM_2D_Lists_t

  type, extends(value_Lists) :: Value_QTM_3D_Lists_t
!   integer :: N ! Number of elements of Eta currently in use.
    type(value_QTM_3D_list_t), allocatable :: Eta(:)
  end type Value_QTM_3D_Lists_t

  interface Convert
    module procedure Convert_QTM_1D_to_P
  end interface

  interface Dot_Product
    module procedure Dot_Product_1D!_d, Dot_Product_1D_s
  end interface

  interface Dump
    module procedure Dump_Value_1D_List, Dump_Value_1D_p_List
!     module procedure Dump_Value_1D_List_d, Dump_Value_1D_List_s
!     module procedure Dump_Value_1D_p_List_d, Dump_Value_1D_p_List_s
    module procedure Dump_Value_2D_List
!     module procedure Dump_Value_2D_List_d, Dump_Value_2D_List_s
    module procedure Dump_Value_3D_List
!     module procedure Dump_Value_3D_List_d, Dump_Value_3D_List_s
  end interface

  interface Interpolate
    module procedure Interpolate_1D_List!_d, Interpolate_1D_List_s
    module procedure Interpolate_1D_Lists!_d, Interpolate_1D_Lists_s
    module procedure Interpolate_1D_P!_d, Interpolate_1D_P_s
    module procedure Interpolate_QTM_1D_List!_d, Interpolate_QTM_1D_List_s
    module procedure Interpolate_2D_List!_d, Interpolate_2D_List_s
    module procedure Interpolate_2D_Lists!_d, Interpolate_2D_Lists_s
    module procedure Interpolate_2D_P!_d, Interpolate_2D_P_s
    module procedure Interpolate_QTM_2D_List!_d, Interpolate_QTM_2D_List_s
    module procedure Interpolate_3D_List!_d, Interpolate_3D_List_s
    module procedure Interpolate_3D_Lists!_d, Interpolate_3D_Lists_s
    module procedure Interpolate_QTM_3D_List!_d, Interpolate_QTM_3D_List_s
  end interface

  interface Interpolate_Polymorphic
    module procedure Interpolate_Polymorphic_1D
!     module procedure Interpolate_Polymorphic_1D_d, Interpolate_Polymorphic_1D_s
    module procedure Interpolate_Polymorphic_2D
!     module procedure Interpolate_Polymorphic_2D_d, Interpolate_Polymorphic_2D_s
    module procedure Interpolate_Polymorphic_3D
!     module procedure Interpolate_Polymorphic_3D_d, Interpolate_Polymorphic_3D_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Convert_QTM_1D_to_P ( From, To, N, Path )
    type(value_QTM_1D_List_t), intent(in) :: From(:)
    type(value_1D_p_t), intent(out) :: To(:)
    integer, intent(out) :: N             ! Number of elements put into To
    logical, intent(in), optional :: Path ! Use NP component of From
    integer :: I, J
    logical :: MyPath

    myPath = .false.
    if ( present(path) ) myPath = path
    n = 0
    do i = 1, size(from)
      do j = 1, from(i)%n
        n = n + 1
        if ( myPath ) then
          to(n) = value_1D_p_t(v=from(i)%v(j)%v, n=from(i)%v(j)%np, p=i )
        else
          to(n) = value_1D_p_t(v=from(i)%v(j)%v, n=from(i)%v(j)%n, p=i )
        end if
      end do
    end do
  end subroutine Convert_QTM_1D_to_P

  pure function Dot_Product_1D ( Vector, List ) result ( Dot )
    real(rk), intent(in) :: Vector(:)
    type(value_1d_list_t), intent(in) :: List
    real(rk) :: Dot
    dot = sum ( vector(list%v(:list%n)%n) * list%v(:list%n)%v )
  end function Dot_Product_1D!_d

  subroutine Dump_Value_1D_List ( Value, Name, Format )
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: NewLine, Output
!     type(value_1D_list_t(rk)), intent(in) :: Value(:)
    class(value_1D_list_t), intent(in) :: Value(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    include "Dump_Value_1D_List.f9h"
  end subroutine Dump_Value_1D_List

  subroutine Dump_Value_1D_p_List ( Value, Name, Format )
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: Output
!     type(value_1D_list_t(rk)), intent(in) :: Value(:)
    class(value_1D_p_t), intent(in) :: Value(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    include "Dump_Value_1D_p_List.f9h"
  end subroutine Dump_Value_1D_p_List

  subroutine Dump_Value_2D_List ( Value, Name, Format )
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: NewLine, Output
!     type(value_2D_list_t(rk)), intent(in) :: Value(:)
    class(value_2D_list_t), intent(in) :: Value(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    include "Dump_Value_2D_List.f9h"
  end subroutine Dump_Value_2D_List

  subroutine Dump_Value_3D_List ( Value, Name, Format )
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: NewLine, Output
!     type(value_2D_list_t(rk)), intent(in) :: Value(:)
    class(value_3D_list_t), intent(in) :: Value(:)
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Format
    include "Dump_Value_3D_List.f9h"
  end subroutine Dump_Value_3D_List

  subroutine Interpolate_Polymorphic_1D ( Field, Eta, Path )
    ! Select one interpolator, depending upon the type of Eta
    real(rk), intent(in) :: Field(:)
    class(value_list), intent(in) :: Eta(:)     ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    select type ( Eta )
    type is ( value_1D_list_t )
      call interpolate ( field, eta, path )
    type is ( value_1D_p_t )
      call interpolate ( field, eta, path )
    type is ( value_QTM_1D_list_t )
      call interpolate ( field, eta, path )
    end select
  end subroutine Interpolate_Polymorphic_1D

  subroutine Interpolate_Polymorphic_2D ( Field, Eta, Path )
    ! Select one of the interpolators, depending upon the type of Eta
    real(rk), intent(in) :: Field(:,:)
    class(value_list), intent(in) :: Eta(:)     ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    select type ( Eta )
    type is ( value_2D_list_t )
      call interpolate ( field, eta, path )
    type is ( value_QTM_2D_list_t )
      call interpolate ( field, eta, path )
    end select
  end subroutine Interpolate_Polymorphic_2D

  subroutine Interpolate_Polymorphic_3D ( Field, Eta, Path )
    ! Select one interpolator, depending upon the type of Eta
    real(rk), intent(in) :: Field(:,:,:)
    class(value_list), intent(in) :: Eta(:)     ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    select type ( Eta )
    type is ( value_3D_list_t )
      call interpolate ( field, eta, path )
    type is ( value_QTM_3D_list_t )
      call interpolate ( field, eta, path )
    end select
  end subroutine Interpolate_Polymorphic_3D

  subroutine Interpolate_1D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:)
    type(value_1D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    include "Interpolate_1D_List.f9h"
  end subroutine Interpolate_1D_List

  subroutine Interpolate_1D_Lists ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:)
    type(value_1D_lists_t), intent(in) :: Eta   ! Eta%n assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    call interpolate_1d_list ( field, eta%eta(:eta%n), path )
  end subroutine Interpolate_1D_Lists

  subroutine Interpolate_1D_P ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! in Eta.
    real(rk), intent(in) :: Field(:)
    type(value_1D_p_t), intent(in) :: Eta(:)    ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    include "Interpolate_1D_P.f9h"
  end subroutine Interpolate_1D_P

  subroutine Interpolate_QTM_1D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:)
    type(value_QTM_1D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)                ! the same as Size(Path)
    include "Interpolate_1D_List.f9h"
  end subroutine Interpolate_QTM_1D_List

  subroutine Interpolate_2D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:)
    type(value_2D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    include "Interpolate_2D_List.f9h"
  end subroutine Interpolate_2D_List

  subroutine Interpolate_2D_Lists ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:)
    type(value_2D_lists_t), intent(in) :: Eta   ! Eta%n assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    call interpolate_2d_list ( field, eta%eta(:eta%n), path )
  end subroutine Interpolate_2D_Lists

  subroutine Interpolate_2D_P ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the second subscript for Field
    ! and Path is in Eta.  The same interpolation is applied to all of the
    ! first subscripts for Field and Path, which are assumed to have the
    ! same extents.  This is for quantities that have several components that
    ! don't depend upon a basis grid, such as magnetic fields.
    real(rk), intent(in) :: Field(:,:)          ! Size(Field,1) assumed to be
    type(value_1D_p_t), intent(in) :: Eta(:)
    real(rk), intent(out) :: Path(:,:)          ! the same as Size(Path,1)
    include "Interpolate_2D_P.f9h"
  end subroutine Interpolate_2D_P

  subroutine Interpolate_QTM_2D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:)
    type(value_QTM_2D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)                ! the same as Size(Path)
    include "Interpolate_2D_List.f9h"
  end subroutine Interpolate_QTM_2D_List

  subroutine Interpolate_3D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:,:)
    type(value_3D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    include "Interpolate_3D_List.f9h"
  end subroutine Interpolate_3D_List

  subroutine Interpolate_3D_Lists ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:,:)
    type(value_3D_lists_t), intent(in) :: Eta   ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)            ! the same as Size(Path)
    call interpolate_3d_list ( field, eta%eta(:eta%n), path )
  end subroutine Interpolate_3D_Lists

  subroutine Interpolate_QTM_3D_List ( Field, Eta, Path )
    ! Interpolate Field * Eta to Path, where the subscript for Path is
    ! the same as the subscript for Eta.
    real(rk), intent(in) :: Field(:,:,:)
    type(value_QTM_3D_list_t), intent(in) :: Eta(:) ! Size(Eta) assumed to be
    real(rk), intent(out) :: Path(:)                ! the same as Size(Path)
    include "Interpolate_3D_List.f9h"
  end subroutine Interpolate_QTM_3D_List

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

end module Indexed_Values_m

! $Log$
! Revision 2.10  2017/02/04 02:06:41  vsnyder
! Add Dot_Product.  Make the list argument of some dumps polymorphic so
! it will work with extensions (but it only dumps the parent part).
!
! Revision 2.9  2017/01/21 01:59:53  vsnyder
! Add Dump_Value_3D_List
!
! Revision 2.8  2017/01/14 01:46:12  vsnyder
! Add List of Lists, add Format argument to dumps
!
! Revision 2.7  2016/12/15 18:45:09  vsnyder
! Move argument declarations from includes to Indexed_Values
!
! Revision 2.6  2016/12/03 02:43:46  vsnyder
! Add Convert_QTM_1D_to_P
!
! Revision 2.5  2016/12/02 01:58:44  vsnyder
! Move argument declarations here from include files.  Add 2D P interpolator.
! Add Dump_Value_1D_p_List.
!
! Revision 2.4  2016/11/29 00:28:54  vsnyder
! Add interpolators
!
! Revision 2.3  2016/11/23 21:32:59  vsnyder
! Add Value_3D_*, Dump
!
! Revision 2.2  2016/11/23 20:07:43  vsnyder
! Add Dump_Value_1D_List
!
! Revision 2.1  2016/11/23 00:05:48  vsnyder
! Initial commit
!
