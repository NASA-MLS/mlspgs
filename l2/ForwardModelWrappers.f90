! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ForwardModelWrappers

  ! This module contains a wrapper routine for calling the various forward
  ! models we have.

  implicit none
  private

  public :: ForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------
  
contains ! ============= Public Procedures ==========================

  !----------------------------------------- ForwardModel -----------
  subroutine ForwardModel ( Config, FwdModelIn, FwdModelExtra, &
    FwdModelOut, fmStat, Jacobian, vectors )

    use BaselineForwardModel_m, only: BASELINEFORWARDMODEL
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
    use FullCloudForwardModel, only: FULLCLOUDFORWARDMODELWRAPPER
    use FullForwardModel_m, only: FULLFORWARDMODEL
    use HybridForwardModel_m, only: HYBRIDFORWARDMODEL
    use Init_tables_module, only: L_LINEAR, L_SCAN, L_SCAN2D, L_FULL, L_CLOUDFULL, &
      & L_SWITCHINGMIRROR, L_HYBRID, L_POLARLINEAR, L_BASELINE
    use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
    use MatrixModule_1, only: MATRIX_T, CHECKINTEGRITY
    use MLSL2Timings, only: Add_to_retrieval_timing
    use MLSMessageModule, only: MLSMessage, MLSMessageCalls, MLSMSG_Error, MLSMSG_Warning
    use PolarLinearModel_m, only: POLARLINEARMODEL
    use ScanModelModule, only: SCANFORWARDMODEL, TWODSCANFORWARDMODEL
    use String_table, only: Display_String, GET_STRING
    use SwitchingMirrorModel_m, only: SWITCHINGMIRRORMODEL
    use Time_M, only: Time_Now
    use Toggles, only: Emit, Switches, Toggle
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use VectorsModule, only: CheckNaN, Dump, VECTOR_T

    ! Dummy arguments
    type(ForwardModelConfig_T), intent(inout) :: CONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(Vector_t), dimension(:), target, optional :: VECTORS ! Vectors database

    ! Local variables
    integer :: K
    real :: time_start, time_end, deltaTime  
    character(len=132) :: THISNAME
    logical :: radianceModel

    ! Executable code
    ! Report we're starting
    if ( config%name /= 0 ) then
      call get_string ( config%name, thisName )
    else
      thisName = '[unnamed]'
    end if

    if ( toggle(emit) ) then
      call trace_begin ( 'ForwardModel ' // trim(thisName) )
    else
      call MLSMessageCalls( 'push', constantName='ForwardModel ' // trim(thisName) )
    end if
    ! Setup the timing
    call time_now (time_start)

    ! Do the actual forward models
    radianceModel = any ( config%fwmType == &
      & (/ l_full, l_linear, l_polarLinear, l_hybrid, l_cloudFull /) )
    select case (config%fwmType)
    case ( l_baseline )
      call MLSMessageCalls( 'push', constantName='BaselineForwardModel' )
      call BaselineForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'baseline' )
    case ( l_full )
      call MLSMessageCalls( 'push', constantName='FullForwardModel' )
      call FullForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'full_fwm' )
    case ( l_linear )
      call MLSMessageCalls( 'push', constantName='LinearizedForwardModel' )
      call LinearizedForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian, vectors )
      call add_to_retrieval_timing( 'linear_fwm' )
    case ( l_hybrid )
      call MLSMessageCalls( 'push', constantName='HybridForwardModel' )
      call HybridForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian, vectors )
      call add_to_retrieval_timing( 'hybrid' )
    case ( l_polarLinear )
      call MLSMessageCalls( 'push', constantName='PolarForwardModel' )
      call PolarLinearModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian, vectors )
      call add_to_retrieval_timing( 'polar_linear' )
    case ( l_scan )
      call MLSMessageCalls( 'push', constantName='ScanForwardModel' )
      call ScanForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'scan_fwm' )
    case ( l_scan2d )
      call MLSMessageCalls( 'push', constantName='TwoDForwardModel' )
      call TwoDScanForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'twod_scan_fwm' )
    case ( l_cloudFull )
      call MLSMessageCalls( 'push', constantName='CloudForwardModel' )
      call FullCloudForwardModelWrapper ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'fullcloud_fwm' )
    case ( l_switchingMirror )
      call MLSMessageCalls( 'push', constantName='SwitchingForwardModel' )
      call SwitchingMirrorModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call add_to_retrieval_timing( 'switching_mirror' )
    case default ! Shouldn't get here if parser etc. worked
    end select
    call MLSMessageCalls( 'pop' ) ! for all the cases

    if ( radianceModel ) then
      call MLSMessageCalls( 'push', constantName='BaselinForwardModel' )
      call BaselineForwardModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call MLSMessageCalls( 'pop' )
      call add_to_retrieval_timing( 'baseline' )
      call MLSMessageCalls( 'push', constantName='SwitchingForwardModel' )
      call SwitchingMirrorModel ( config, FwdModelIn, FwdModelExtra, &
        FwdModelOut, fmStat, Jacobian )
      call MLSMessageCalls( 'pop' )
      call add_to_retrieval_timing( 'switching_mirror' )
    end if

    k = 0
    if ( index(switches,'FMNAN') > 0 ) then
      k = 3 ! Check, print and stop if any
    else if ( index(switches,'fmNaN') > 0 ) then
      k = 2 ! Check, print if any
    else if ( index(switches,'fmnan') > 0 ) then
      k = 1 ! Check, print name if any
    end if
    if ( k > 0 ) then
      ! Check radiances
      if ( checkNaN(fwdModelOut, k-1, 'ForwardModelOut') ) then
        if ( k > 1 ) then
          call dump ( fwdModelIn, k-1, 'ForwardModelIn' )
          call dump ( fwdModelExtra, k-1, 'ForwardModelExtra' )
        end if
        call display_string ( config%name, &
          & before='Forward model config name: ', advance='yes' )
        if ( k > 2 ) then
          k = MLSMSG_Error
        else
          k = MLSMSG_Warning
        end if
        call MLSMessage ( k, ModuleName, 'NaNs found in forward model output' )
      end if

      ! Check Jacobians if relevant
      if ( present ( Jacobian ) ) then 
        if ( .not. checkIntegrity ( Jacobian, noError=.true. ) ) then
          if ( k > 2 ) then
            k = MLSMSG_Error
          else
            k = MLSMSG_Warning
          end if
          call MLSMessage ( k, ModuleName, 'Problem (NANs?) found in Jacobians' )
        end if
      end if
      
    end if
      
    ! Do the timing stuff
    call time_now (time_end)
    deltaTime = time_end - time_start
    config%Ntimes = config%Ntimes + 1
    config%sum_DeltaTime = &
      & config%sum_DeltaTime + deltaTime
    config%sum_squareDeltaTime = &
      & config%sum_squareDeltaTime + (deltaTime * deltaTime)

    ! Report we're finished
    if ( toggle(emit) ) then
      call trace_end ( 'ForwardModel ' // trim(thisName) )
    else
      call MLSMessageCalls( 'pop' )
    end if

  end subroutine ForwardModel

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ForwardModelWrappers

! $Log$
! Revision 2.32  2009/11/18 22:18:07  livesey
! Added ability to check Jacobians as well as radiances is the various
! fmNAN flags are set
!
! Revision 2.31  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.30  2007/10/04 01:48:30  vsnyder
! Make sure thisName has a value, handle call stack properly
!
! Revision 2.29  2007/10/02 22:38:19  vsnyder
! Add code to check for NaNs in forward models' output
!
! Revision 2.28  2007/08/20 22:05:06  pwagner
! Many procedures now push their names onto MLSCallStack
!
! Revision 2.27  2007/06/29 19:32:07  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.26  2007/04/03 17:46:12  vsnyder
! Replace pointer attribute on VectorDatabase with target attribute
!
! Revision 2.25  2005/06/03 02:08:24  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.24  2003/10/20 18:22:47  pwagner
! New forwardModel types added to RetrievalTimings breakdown
!
! Revision 2.23  2003/09/11 23:15:10  livesey
! Added vectors argument which is handed on to some but not all models.
! This is needed to support the xStar/yStar capability of the linear
! forward model (and by inference all those that call it.)
!
! Revision 2.22  2003/08/16 01:18:29  livesey
! Added baseline forward model on its own.
!
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
