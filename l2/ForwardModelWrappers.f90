! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ForwardModelWrappers

  ! This module contains a wrapper routine for calling the various forward
  ! models we have.

  implicit none
  private

  public :: ForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------
  
contains ! ============= Public Procedures ==========================

  !----------------------------------------- ForwardModel -----------
  subroutine ForwardModel ( Config, FwdModelIn, FwdModelExtra, &
    FwdModelOut, Ifm, fmStat, Jacobian )

    use BaselineForwardModel_m, only: BASELINEFORWARDMODEL
    use HybridForwardModel_m, only: HYBRIDFORWARDMODEL
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, &
      & FORWARDMODELSTATUS_T
    use FullCloudForwardModel, only: FULLCLOUDFORWARDMODELWRAPPER
    use FullForwardModel_m, only: FULLFORWARDMODEL
    use Init_tables_module, only: L_LINEAR, L_SCAN, L_SCAN2D, L_FULL, L_CLOUDFULL, &
      & L_SWITCHINGMIRROR, L_HYBRID, L_POLARLINEAR
    use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
    use PolarLinearModel_m, only: POLARLINEARMODEL
    use MatrixModule_1, only: MATRIX_T
    use MLSL2Timings, only: Add_to_retrieval_timing
    use ScanModelModule, only: SCANFORWARDMODEL, TWODSCANFORWARDMODEL
    use SwitchingMirrorModel_m, only: SWITCHINGMIRRORMODEL
    use VectorsModule, only: VECTOR_T
    use Time_M, only: Time_Now
    use String_table, only: GET_STRING
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Toggles, only: Emit, Toggle

    ! Dummy arguments
    type(ForwardModelConfig_T), intent(inout) :: CONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    
    ! Local variables
    real :: time_start, time_end, deltaTime  
    character(len=132) :: THISNAME
    logical :: radianceModel

    ! Executable code
    ! Report we're starting
    if ( toggle(emit) ) then
      if ( config%name /= 0 ) then
        call get_string ( config%name, thisName )
      else
        thisName = '[unnamed]'
      end if
      call trace_begin ( 'Forward model config: ' // trim(thisName) )
    end if

    ! Setup the timing
    call time_now (time_start)

    ! Do the actual forward models
    radianceModel = any ( config%fwmType == &
      & (/ l_full, l_linear, l_polarLinear, l_hybrid, l_cloudFull /) )
    select case (config%fwmType)
    case ( l_full )
      call FullForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'full_fwm' )
    case ( l_linear )
      call LinearizedForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'linear_fwm' )
    case ( l_hybrid )
      call HybridForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'full_fwm' )
    case ( l_polarLinear )
      call PolarLinearModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'full_fwm' )
    case ( l_scan )
      call ScanForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'scan_fwm' )
    case ( l_scan2d )
      call TwoDScanForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'twod_scan_fwm' )
    case ( l_cloudFull )
      call FullCloudForwardModelWrapper ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call add_to_retrieval_timing( 'fullcloud_fwm' )
    case ( l_switchingMirror )
      call SwitchingMirrorModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
    case default ! Shouldn't get here if parser etc. worked
    end select
    
    if ( radianceModel ) then
      call BaselineForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
      call SwitchingMirrorModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
    end if

    ! Report we're finished
    if ( toggle(emit) ) then
      call trace_end ( 'Forward model config: ' // trim(thisName) )
    end if

    ! Do the timing stuff
    call time_now (time_end)
    deltaTime = time_end - time_start
    config%Ntimes = config%Ntimes + 1
    config%sum_DeltaTime = &
      & config%sum_DeltaTime + deltaTime
    config%sum_squareDeltaTime = &
      & config%sum_squareDeltaTime + (deltaTime * deltaTime)
    
  end subroutine ForwardModel

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ForwardModelWrappers

! $Log$
! Revision 2.21  2003/08/13 00:49:56  livesey
! Added PolarLinear model
!
! Revision 2.20  2003/07/15 22:11:12  livesey
! Added hybrid model and slight reorganization
!
! Revision 2.19  2003/07/15 18:18:39  livesey
! Made timing apply to all configs.
!
! Revision 2.18  2003/06/30 22:55:01  cvuu
! Find mean, std dev timing of fullForwardModel calls
!
! Revision 2.17  2003/06/03 19:24:56  livesey
! Added the ability to call the switching mirror model in isolation
!
! Revision 2.16  2003/05/29 16:42:34  livesey
! Added calls to SwitchingMirrorModel
!
! Revision 2.15  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.14  2002/08/21 23:43:33  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.13  2002/07/23 00:06:05  pwagner
! No upper-case allowed in section names
!
! Revision 2.12  2002/07/22 22:51:56  pwagner
! Restored name of 2d scan model in timings
!
! Revision 2.11  2002/06/24 22:14:59  livesey
! Changed name of 2d scan model in timings
!
! Revision 2.10  2002/06/24 18:27:09  livesey
! New 2D scan model
!
! Revision 2.9  2001/11/27 23:34:49  pwagner
! Split forward model timings into four types
!
! Revision 2.8  2001/10/02 16:55:10  livesey
! Bug fix, forgot use statement
!
! Revision 2.7  2001/10/02 16:53:18  livesey
! Added call to BaselineForwardModel for Full and Linearized forward models.
!
! Revision 2.6  2001/07/17 22:36:32  jonathan
! add cloud_width, jonathan/paul
!
! Revision 2.5  2001/05/29 23:22:20  livesey
! FullForwardModel moved, also added (but commented out)
! call to FullCloudForwardModelWrapper
!
! Revision 2.4  2001/05/03 23:42:48  livesey
! Activated scan model.
!
! Revision 2.3  2001/04/28 17:48:48  livesey
! Removed some unnecessary checks
!
! Revision 2.2  2001/04/26 23:54:26  livesey
! Now uses linear forward model
!
! Revision 2.1  2001/04/26 19:47:41  livesey
! First version
!
