! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelIntermediate
!=============================================================================

  ! This module simply defines a data type passed around forward model routines. 
  ! It contains hydrostatic information and so on, and is remembered by the
  ! calling code.

  use MLSCommon, only: R8
  use Path_entities_m, only: PATH_VECTOR, PATH_VECTOR_2D
  use Ellipse_m, only: ELLIPSE

  implicit none
  private

  type, public :: ForwardModelIntermediate_T

    real (r8), dimension(:,:),             pointer :: h_glgrid
    real (r8), dimension(:,:),             pointer :: t_glgrid
    real (r8), dimension(:),               pointer :: z_glgrid

    real (r8), dimension(:,:,:),           pointer :: dh_dt_glgrid
    real (r8), dimension(:,:),             pointer :: dhdz_glgrid

    real (r8), dimension(:,:),             pointer :: tan_hts
    real (r8), dimension(:,:),             pointer :: tan_temp
    real (r8), dimension(:,:,:),           pointer :: tan_dh_dt

    type (path_vector), dimension(:,:),    pointer :: z_path
    type (path_vector), dimension(:,:),    pointer :: h_path
    type (path_vector), dimension(:,:),    pointer :: t_path

    type (path_vector), dimension(:,:),    pointer :: phi_path
    type (path_vector), dimension(:,:),    pointer :: dhdz_path
    type (path_Vector_2d), dimension(:,:), pointer :: eta_phi
    type (ellipse), dimension(:,:),        pointer :: elvar

  end type ForwardModelIntermediate_T

  type, public :: ForwardModelStatus_T
    logical :: newHydros                ! Need to recompute hydrostatic
    integer :: maf                      ! The next maf to process
    logical :: finished                 ! Flag to calling code to indicate completion
  end type ForwardModelStatus_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

end module ForwardModelIntermediate

! $Log$
