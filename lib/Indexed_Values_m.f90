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

  ! For one interpolation weight from a 1D array
  type :: Value_1D_t ! ( RK )
!     integer, kind :: RK
    real(rk) :: V = 0.0    ! Value to be applied at N
    integer :: N = 0       ! Subscript at which to apply V
  end type Value_1D_t

  ! For interpolating to a list from a 1D array, assuming linear interolation
  type :: Value_1D_List_t ! ( RK )
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
  type :: Value_QTM_1D_List_t ! ( RK )
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
  type :: Value_2D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 4       ! Number of useful elements of V
!     type(value_2d_t(rk)) :: V(4) = value_2d_t(rk)()
    type(value_2d_t) :: V(4) = value_2d_t()
  end type Value_2D_List_t

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
  type :: Value_QTM_2D_List_t ! ( RK )
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
  type :: Value_QTM_3D_List_t ! ( RK )
!     integer, kind :: RK
    integer :: N = 12      ! Number of useful elements of V
!     type(value_QTM_3d_t(rk)) :: V(12) = value_QTM_3d_t(rk)()
    type(value_QTM_3d_t) :: V(12) = value_QTM_3d_t()
  end type Value_QTM_3D_List_t

  interface Dump
    module procedure Dump_Value_1D_List
!     module procedure Dump_Value_1D_List_d, Dump_Value_1D_List_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Dump_Value_1D_List ( Value, Name )
    use Output_m, only: NewLine, Output
    include "Dump_Value_1D_List.f9h"
  end subroutine Dump_Value_1D_List

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
! Revision 2.2  2016/11/23 20:07:43  vsnyder
! Add Dump_Value_1D_List
!
! Revision 2.1  2016/11/23 00:05:48  vsnyder
! Initial commit
!
