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
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_Deallocate, MLSMSG_ERROR
  use Allocate_Deallocate, only: DEALLOCATE_TEST

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

    integer, dimension(:),                 pointer :: closestInstances=>NULL()

    ! These ones are for the scan model
    real (r8), dimension(:,:),             pointer :: basisGph=>NULL()
    real (r8), dimension(:,:),             pointer :: RT=>NULL()
    real (r8), dimension(:),               pointer :: R=>NULL()
    integer :: belowRef                 ! T. basis at or below refGPH

  end type ForwardModelIntermediate_T

  type, public :: ForwardModelStatus_T
    logical :: newHydros                ! Need to recompute hydrostatic
    logical :: newScanHydros            ! Scan model needs to recompute hydrostatic
    integer :: maf                      ! The next maf to process
    logical :: finished                 ! Flag to calling code to indicate completion
    logical, dimension(:), pointer :: rows=>NULL() ! Flag to indicate this row has non zeros
  end type ForwardModelStatus_T

  public :: DestroyForwardModelIntermediate

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------------ DestroyForwardModelIntemediate ---
  subroutine DestroyForwardModelIntermediate ( ifm )
    type (ForwardModelIntermediate_T), intent(inout) :: ifm

    ! Local variables
    integer :: i,j                      ! Loop counters
    integer :: noMAFs, no_tan_hts       ! Dimensions
    integer :: status                   ! Flag

    ! Exectuable code

    if ( associated( ifm%z_path) ) then
      no_tan_hts = size ( ifm%z_path, 1)
      noMAFs = size ( ifm%z_path, 2)
      
      deallocate (ifm%ndx_path, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Deallocate//'ndx_path' )
      
      do j = 1, noMAFs
        do i = 1, No_tan_hts
          deallocate (ifm%z_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'z_path%values' )
          deallocate (ifm%h_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'h_path%values' )
          deallocate (ifm%t_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'t_path%values' )
          deallocate (ifm%phi_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'phi_path%values' )
          deallocate (ifm%dhdz_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'dhdz_path%values' )
          deallocate (ifm%eta_phi(i,j)%values, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Deallocate//'eta_phi%values' )
        end do
      end do
      
      deallocate (ifm%z_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Deallocate//'z_path' )
    endif

    deallocate (ifm%h_path, stat=status )
!     if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'h_path' )
    Deallocate (ifm%t_path, stat=status )
!     if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'t_path' )
    deallocate (ifm%phi_path, stat=status )
!     if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'phi_path' )
    deallocate (ifm%dhdz_path, stat=status )
!     if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'dhdz_path' )
    deallocate (ifm%eta_phi, stat=status )
!     if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'eta_phi' )

    deallocate ( ifm%elvar, stat=status )
!     if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
!       & MLSMSG_Deallocate//'elvar' )

    call deallocate_test ( ifm%geoc_lat, 'geoc_lat', ModuleName )
    call deallocate_test ( ifm%e_rad, 'e_rad', ModuleName )

    call deallocate_test ( ifm%h_glgrid, 'h_glgrid', ModuleName )
    call deallocate_test ( ifm%t_glgrid, 't_glgrid', ModuleName )
    call deallocate_test ( ifm%z_glgrid, 'z_glgrid', ModuleName )
    call deallocate_test ( ifm%dh_dt_glgrid, 'dh_dt_glgrid', ModuleName )
    call deallocate_test ( ifm%dhdz_glgrid, 'dhdz_glgrid', ModuleName )
    call deallocate_test ( ifm%tan_hts,'tan_hts', ModuleName )
    call deallocate_test ( ifm%tan_temp,'tan_temp', ModuleName )
    call deallocate_test ( ifm%tan_dh_dt, 'tan_dh_dt', ModuleName )

    call deallocate_test ( ifm%closestInstances, 'closestInstances', ModuleName )
    
  end subroutine DestroyForwardModelIntermediate

end module ForwardModelIntermediate

! $Log$
! Revision 1.9  2001/05/03 23:07:12  livesey
! Added scan model stuff
!
! Revision 1.8  2001/05/01 00:23:00  livesey
! Nullified fmStat%rows by default
!
! Revision 1.7  2001/04/28 07:05:54  livesey
! Tidy up destruction of ForwardModelIntermediate, for cases where
! forward model has not been run.
!
! Revision 1.6  2001/04/23 21:40:46  livesey
! Moved closestInstances into ifm
!
! Revision 1.5  2001/04/19 23:57:00  livesey
! New fmStat
!
! Revision 1.4  2001/04/19 21:04:42  livesey
! Bug fix, sized arrays wrongly
!
! Revision 1.3  2001/04/19 20:29:46  livesey
! Added destroy routine
!
! Revision 1.2  2001/04/10 23:16:14  livesey
! Working version.
!
! Revision 1.1  2001/04/10 22:17:05  livesey
! Renamed module
!
! Revision 1.1  2001/04/10 22:06:26  livesey
! Very first version.  Bare bones only.
!
