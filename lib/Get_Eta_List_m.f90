! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Eta_List_m

  use Indexed_Values_m, only: Value_List, Value_1D_List_t, Value_1D_p_t, &
    & Value_1D_t, Value_2D_List_t, Value_3D_List_t, &
    & Value_QTM_1D_List_t, Value_QTM_2D_List_t, Value_QTM_3D_List_t

  implicit NONE

  private

  ! Type
  public :: Eta_Lists_t

  ! Procedures
  public :: Eta_List_1D, Eta_List_1D_1D_2D, Eta_List_1D_2D_3D
  public :: Eta_List_1D_QTM_1D_2D, Eta_List_1D_QTM_2D_3D

  ! Compute Eta to interpolate from Basis to Grid
  interface Eta_List_1D
    module procedure Eta_List_1D_D, Eta_List_1D_S
    module procedure Eta_List_1D_p_D, Eta_List_1D_p_S
    module procedure Eta_List_1D_1D_2D!_d, Eta_List_1D_1D_2D_s
    module procedure Eta_List_1D_2D_3D!_d, Eta_List_1D_2D_3D_s
    module procedure Eta_List_1D_QTM_1D_2D!_d, Eta_List_1D_QTM_1D_2D_s
    module procedure Eta_List_1D_QTM_2D_3D!_d, Eta_List_1D_QTM_2D_3D_s
  end interface

  ! Use this type to create an array of Eta lists, potentially of different
  ! types and sizes.
  type :: Eta_Lists_t
    class(value_list), allocatable :: Eta(:)
  end type Eta_Lists_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid.
  subroutine Eta_List_1D_D ( Basis, Grid, Eta, Sorted )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Basis(:)
    real(rk), intent(in) :: Grid(:)
  !   type(value_1D_list_t(rk)), intent(out) :: Eta(:) ! size(grid)
    type(value_1D_list_t), intent(out) :: Eta(:) ! size(grid)
    logical, intent(in), optional :: Sorted ! "Basis is sorted" -- default true
    include "Eta_List_1D.f9h"
  end subroutine Eta_List_1D_D

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid.
  subroutine Eta_List_1D_S ( Basis, Grid, Eta, Sorted )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Basis(:)
    real(rk), intent(in) :: Grid(:)
  !   type(value_1D_list_t(rk)), intent(out) :: Eta(:) ! size(grid)
    type(value_1D_list_t), intent(out) :: Eta(:) ! size(grid)
    logical, intent(in), optional :: Sorted ! "Basis is sorted" -- default true
    include "Eta_List_1D.f9h"
  end subroutine Eta_List_1D_S

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid, returning also
  ! the amount of Eta actually used.
  subroutine Eta_List_1D_p_D ( Basis, Grid, Eta, N, Sorted )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: Basis(:)
    real(rk), intent(in) :: Grid(:)
  !   type(value_1D_p_t(rk)), intent(out) :: Eta(:) ! 2*size(grid)
    type(value_1D_p_t), intent(out) :: Eta(:) ! 2*size(grid)
    integer, intent(out) :: N
    logical, intent(in), optional :: Sorted ! "Basis is sorted" -- default true
    include "Eta_List_1D_p.f9h"
  end subroutine Eta_List_1D_p_D

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid, returning also
  ! the amount of Eta actually used.
  subroutine Eta_List_1D_p_S ( Basis, Grid, Eta, N, Sorted )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: Basis(:)
    real(rk), intent(in) :: Grid(:)
  !   type(value_1D_p_t(rk)), intent(out) :: Eta(:) ! 2*size(grid)
    type(value_1D_p_t), intent(out) :: Eta(:) ! 2*size(grid)
    integer, intent(out) :: N
    logical, intent(in), optional :: Sorted ! "Basis is sorted" -- default true
    include "Eta_List_1D_p.f9h"
  end subroutine Eta_List_1D_p_S

  ! Compute Eta to interpolate from 2-D basis to 1-D grid.  The matrix C is
  ! the outer product of interpolation matrices A and B, all represented as
  ! sparse lists.  The field for which A X B is the basis is assumed to be
  ! represented by a rank-one array with elements of the "real" rank-two field
  ! F(NA,*) taken in array-element order.  A%N are first subscripts of F. B%N
  ! are second subscripts of F.  "Sorted" (optional, default true) indicates
  ! whether both A and B are sorted on their P component.  NC is the number of
  ! elements put into C, which might be less than size(A) * size(B) if any
  ! elements of A%V or B%V are zero.
  subroutine Eta_List_1D_1D_p_2D ( A, NA, B, C, NC, Sorted )
!   type(value_1D_p_t(rk)), intent(in) :: A(:)
!   type(value_1D_p_t(rk)), intent(in) :: B(:)
!   type(value_1D_p_t(rk)), intent(out) :: C(:)
    type(value_1D_p_t), intent(in) :: A(:)
    type(value_1D_p_t), intent(in) :: B(:)
    type(value_1D_p_t), intent(out) :: C(:) ! at least size(A) * size(B)
    integer, intent(in) :: NA
    integer, intent(out) :: NC ! Number of used elements of C
    logical, intent(in), optional :: Sorted
    include "Eta_List_1D_1D_p_2D.f9h"
  end subroutine Eta_List_1D_1D_p_2D

  ! Compute Eta to interpolate from 2-D basis to 1-D grid.  The matrix C
  ! is the outer product of interpolation matrices A and B, all represented
  ! as sparse lists.
  subroutine Eta_List_1D_1D_2D ( A, B, C )
!   type(value_1D_list_t(rk)), intent(in) :: A(:)
!   type(value_1D_list_t(rk)), intent(in) :: B(:)
!   type(value_2D_list_t(rk)), intent(out) :: C(:)
    type(value_1D_list_t), intent(in) :: A(:)
    type(value_1D_list_t), intent(in) :: B(:)
    type(value_2D_list_t), intent(out) :: C(:)
    include "Eta_List_1D_1D_2D.f9h"
  end subroutine Eta_List_1D_1D_2D

  ! Compute Eta to interpolate from 3-D basis to 1-D grid.  The matrix C
  ! is the outer product of interpolation matrices A and B, all represented
  ! as sparse lists.
  subroutine Eta_List_1D_2D_3D ( A, B, C )
!   type(value_1D_list_t(rk)), intent(in) :: A(:)
!   type(value_2D_list_t(rk)), intent(in) :: B(:)
!   type(value_3D_list_t(rk)), intent(out) :: C(:)
    type(value_1D_list_t), intent(in) :: A(:)
    type(value_2D_list_t), intent(in) :: B(:)
    type(value_3D_list_t), intent(out) :: C(:)
    include "Eta_List_1D_2D_3D.f9h"
  end subroutine Eta_List_1D_2D_3D

  ! Compute Eta to interpolate from a 2-D QTM basis and a 1-D vertical
  ! basis to a 1-D grid.  The matrix C is the outer product of
  ! interpolation matrices A and B, all represented as sparse lists.
  subroutine Eta_List_1D_QTM_1D_2D ( A, B, C )
!   type(value_1D_list_t(rk)), intent(in) :: A(:)
!   type(value_qtm_1d_list_t(rk)), intent(in) :: B(:)
!   type(value_qtm_2d_list_t(rk)), intent(out) :: C(:)
    type(value_1D_list_t), intent(in) :: A(:)
    type(value_qtm_1d_list_t), intent(in) :: B(:)
    type(value_qtm_2d_list_t), intent(out) :: C(:)
    include "Eta_List_1D_QTM_1D_2D.f9h"
  end subroutine Eta_List_1D_QTM_1D_2D

  ! Compute Eta to interpolate from a 2-D QTM basis, a 1-D vertical, and a
  ! 1-D frequency  basis to a 1-D grid.  The matrix C is the outer product
  ! of interpolation matrices A and B, all represented as sparse lists.
  subroutine Eta_List_1D_QTM_2D_3D ( A, B, C )
!   type(value_1D_list_t(rk)), intent(in) :: A(:)
!   type(value_qtm_2d_list_t(rk)), intent(in) :: B(:)
!   type(value_qtm_3d_list_t(rk)), intent(out) :: C(:)
    type(value_1D_list_t), intent(in) :: A(:)
    type(value_qtm_2d_list_t), intent(in) :: B(:)
    type(value_qtm_3d_list_t), intent(out) :: C(:)
    include "Eta_List_1D_QTM_2D_3D.f9h"
  end subroutine Eta_List_1D_QTM_2D_3D

!=========================================================================

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Eta_List_m
!---------------------------------------------------
! $Log$
! Revision 2.3  2016/12/02 02:02:01  vsnyder
! Move dummy argument declarations here from include files.  Add Eta_Lists_t.
! Add 1D Eta P generators.  Add outer product of 1D P Etas.
!
! Revision 2.2  2016/11/23 21:32:25  vsnyder
! Add Eta_List_1D_... routines to compute outer products of lower-
! dimensionality lists.
!
! Revision 2.1  2016/11/23 00:05:47  vsnyder
! Initial commit
!
