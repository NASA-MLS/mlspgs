! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelIntermediate
!=============================================================================

  ! This module simply defines a data type passed around forward model routines. 
  ! It contains hydrostatic information and so on, and is remembered by the
  ! calling code.

  use MLSCommon, only: R8
  use Path_entities_m, only: PATH_VECTOR, PATH_VECTOR_2D, PATH_INDEX
  use Ellipse_m, only: ELLIPSE

  implicit none
  private

  type, public :: ForwardModelIntermediate_T

    real (r8), dimension(:,:),             pointer :: H_GLGRID=>NULL()
    real (r8), dimension(:,:),             pointer :: T_GLGRID=>NULL()
    real (r8), dimension(:),               pointer :: z_glgrid=>NULL()

    real (r8), dimension(:,:,:),           pointer :: dh_dt_glgrid=>NULL()
    real (r8), dimension(:,:),             pointer :: dhdz_glgrid=>NULL()
    integer                                        :: gl_count

    real (r8), dimension(:,:),             pointer :: tan_hts=>NULL()
    real (r8), dimension(:,:),             pointer :: tan_temp=>NULL()
    real (r8), dimension(:,:,:),           pointer :: tan_dh_dt=>NULL()
    real (r8), dimension(:),               pointer :: GEOC_LAT=>NULL()
    real (r8), dimension(:),               pointer :: E_RAD=>NULL()

    type (path_vector), dimension(:,:),    pointer :: z_path=>NULL()
    type (path_vector), dimension(:,:),    pointer :: h_path=>NULL()
    type (path_vector), dimension(:,:),    pointer :: t_path=>NULL()

    type (path_vector), dimension(:,:),    pointer :: phi_path=>NULL()
    type (path_vector), dimension(:,:),    pointer :: dhdz_path=>NULL()
    type (path_Vector_2d), dimension(:,:), pointer :: eta_phi=>NULL()
    type (path_index),  dimension(:,:),    pointer :: NDX_PATH=>NULL()
    type (ellipse), dimension(:),          pointer :: elvar=>NULL()

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
! Revision 1.1  2001/04/10 22:17:05  livesey
! Renamed module
!
! Revision 1.1  2001/04/10 22:06:26  livesey
! Very first version.  Bare bones only.
!
