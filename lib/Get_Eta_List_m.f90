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

  use Indexed_Values_m, only: Value_1D_List_t, Value_1D_t, &
    & Value_2D_List_t, Value_3D_List_t, &
    & Value_QTM_1D_List_t, Value_QTM_2D_List_t, Value_QTM_3D_List_t

  implicit NONE

  private

  public :: Eta_List_1D, Eta_List_1D_1D_2D, Eta_List_1D_2D_3D
  public :: Eta_List_1D_QTM_1D_2D, Eta_List_1D_QTM_2D_3D

  ! Compute Eta to interpolate from Basis to Grid
  interface Eta_List_1D
    module procedure Eta_List_1D_D, Eta_List_1D_S
    module procedure Eta_List_1D_1D_2D!_d, Eta_List_1D_1D_2D_s
    module procedure Eta_List_1D_2D_3D!_d, Eta_List_1D_2D_3D_s
    module procedure Eta_List_1D_QTM_1D_2D!_d, Eta_List_1D_QTM_1D_2D_s
    module procedure Eta_List_1D_QTM_2D_3D!_d, Eta_List_1D_QTM_2D_3D_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid
  subroutine Eta_List_1D_D ( Basis, Grid, Eta, Sorted )
    integer, parameter :: RK = kind(0.0d0)
    include "Eta_List_1D.f9h"
  end subroutine Eta_List_1D_D

  ! Compute Eta to interpolate from 1-D Basis to 1-D Grid
  subroutine Eta_List_1D_S ( Basis, Grid, Eta, Sorted )
    integer, parameter :: RK = kind(0.0e0)
    include "Eta_List_1D.f9h"
  end subroutine Eta_List_1D_S

  ! Compute Eta to interpolate from 2-D basis to 1-D grid.  The matrix C
  ! is the outer product of interpoation matrices A and B, all represented
  ! as sparse lists.
  subroutine Eta_List_1D_1D_2D ( A, B, C )
    include "Eta_List_1D_1D_2D.f9h"
  end subroutine Eta_List_1D_1D_2D

  ! Compute Eta to interpolate from 3-D basis to 1-D grid.  The matrix C
  ! is the outer product of interpoation matrices A and B, all represented
  ! as sparse lists.
  subroutine Eta_List_1D_2D_3D ( A, B, C )
    include "Eta_List_1D_2D_3D.f9h"
  end subroutine Eta_List_1D_2D_3D

  ! Compute Eta to interpolate from a 2-D QTM basis and a 1-D vertical
  ! basis to a 1-D grid.  The matrix C is the outer product of
  ! interpoation matrices A and B, all represented as sparse lists.
  subroutine Eta_List_1D_QTM_1D_2D ( A, B, C )
    include "Eta_List_1D_QTM_1D_2D.f9h"
  end subroutine Eta_List_1D_QTM_1D_2D

  ! Compute Eta to interpolate from a 2-D QTM basis, a 1-D vertical, and a
  ! 1-D frequency  basis to a 1-D grid.  The matrix C is the outer product
  ! of interpoation matrices A and B, all represented as sparse lists.
  subroutine Eta_List_1D_QTM_2D_3D ( A, B, C )
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
! Revision 2.2  2016/11/23 21:32:25  vsnyder
! Add Eta_List_1D_... routines to compute outer products of lower-
! dimensionality lists.
!
! Revision 2.1  2016/11/23 00:05:47  vsnyder
! Initial commit
!
