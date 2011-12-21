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
    & FwdModelOut, fmStat, Jacobian, Hessian, Vectors, FwmState, FwmJacobian )

    ! Call the forward model selected by Config.

    ! If FwmState and FwmJacobian are present (both or neither), transform
    ! FwdModelIn to FwmState before calling the forward model, and transform
    ! FwmJacobian to Jacobian after returning from the forward model.

    use BASELINEFORWARDMODEL_M, only: BASELINEFORWARDMODEL
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELINTERMEDIATE, only: FORWARDMODELSTATUS_T
    use FULLCLOUDFORWARDMODEL, only: FULLCLOUDFORWARDMODELWRAPPER
    use FULLFORWARDMODEL_M, only: FULLFORWARDMODEL
    use HESSIANMODULE_1, only: HESSIAN_T
    use HYBRIDFORWARDMODEL_M, only: HYBRIDFORWARDMODEL
    use INIT_TABLES_MODULE, only: L_LINEAR, L_SCAN, L_SCAN2D, L_FULL, &
      & L_CLOUDFULL, L_SWITCHINGMIRROR, L_HYBRID, L_MIFExtinction, &
      & L_MIFExtinctionV2, L_POLARLINEAR, L_BASELINE
    use Intrinsic, only: L_LowestRetrievedPressure, L_VMR
    use LINEARIZEDFORWARDMODEL_M, only: LINEARIZEDFORWARDMODEL
    use MATRIXMODULE_1, only: CHECKINTEGRITY, FindBlock, MATRIX_T
    use MLSL2TIMINGS, only: ADD_TO_RETRIEVAL_TIMING
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMESSAGECALLS, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use Molecules, only: First_Molecule, Last_Molecule, &
      & L_Extinction, L_ExtinctionV2
    use POLARLINEARMODEL_M, only: POLARLINEARMODEL
    use SCANMODELMODULE, only: SCANFORWARDMODEL, TWODSCANFORWARDMODEL
    use STRING_TABLE, only: DISPLAY_STRING, GET_STRING
    use SWITCHINGMIRRORMODEL_M, only: SWITCHINGMIRRORMODEL
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: EMIT, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VECTORSMODULE, only: CHECKNAN, DUMP, GetVectorQuantityByType, &
      & MoveVectorQuantity, VECTOR_T, VECTORVALUE_T

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: CONFIG
    type(vector_T), intent(inout), target :: FWDMODELIN ! The retriever state,
      ! and the forward model state if no transformations are requested
    type(vector_T), intent(in) :: FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional, target :: JACOBIAN ! The retriever
      ! Jacobian, used for the Newton iteration.  The Jacobian produced
      ! by the forward model if no transformations are requested, else
      ! transformed from FwmJacobian
    type(Hessian_T), intent(inout), optional :: HESSIAN ! No transformation
    type(Vector_t), dimension(:), target, optional :: VECTORS ! Vectors database
    type(Vector_t), intent(in), optional, target :: FwmState ! The forward
      ! model state, transformed from FwdModelIn if transformations requested
    type(matrix_T), intent(inout), optional, target :: FwmJacobian ! The
      ! Jacobian produced by the forward model if transformations are requested.

    ! Local variables
    real :: DeltaTime
    integer :: DumpTransform            ! Dump transformed vector quantities
    integer :: FRow, FCol               ! Row, Column of block in fwmJacobian
    integer :: I, K
    integer :: Inst                     ! Instance index
    integer :: JRow, JCol               ! Row, Column of block in Jacobian
    type(vectorValue_t), pointer :: LRP ! Lowest Retrieved Pressure
    type(vectorValue_t), pointer :: Qty ! Transformed quantity in fwmState
    logical :: RadianceModel
    type(vector_t), pointer :: TheState
    type(vectorValue_t), pointer :: T_Qty ! Quantity in State to be transformed
    type(matrix_T), pointer :: TheJacobian
    character(len=132) :: THISNAME
    real :: Time_start, Time_end 

    interface MINLOC_S
      module procedure MINLOC_S_S, MINLOC_S_D
    end interface

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

    dumpTransform = switchDetail(switches,'dxfq')

    ! Do the actual forward models
    if ( present(fwmJacobian) ) then ! Need transformations
      ! Lowest Returned Pressure is needed for extinction transformations
      lrp => getVectorQuantityByType ( fwdModelExtra, &
               & quantityType=l_lowestRetrievedPressure )
      ! Move quantities in State for which no transformation is requested to
      ! FwmState.  Transform those for which transformation is requested. We
      ! move the quantities so as to provide a consistent view of the state
      ! to the Forward model.  The un-transformed quantities in FwdModelIn
      ! are moved instead of copied to avoid aliasing and memory leaks.
      do k = 1, size(fwdModelIn%quantities)
        t_qty => fwdModelIn%quantities(k)    ! Quantity to be transformed
        select case ( t_qty%template%quantityType )
        case ( l_MIFExtinction, l_MIFExtinctionV2 )
          qty => getVectorQuantityByType ( & ! Transformed quantity
            &    fwmState, quantityType=t_qty%template%quantityType, &
            &    radiometer=t_qty%template%radiometer )
          call transform_MIF_extinction &
            & ( fmStat%MAF, T_Qty, lrp%values(1,1), dumpTransform, Qty )
        case default ! Just move the quantity
          qty => getVectorQuantityByType ( fwmState, &
              &  quantityType=t_qty%template%quantityType )
          call moveVectorQuantity ( t_qty, qty )
        end select
      end do

      ! Run the forward model with the transformed quantities
      call doForwardModels ( fwmState, fwmJacobian )

      ! Move quantities in FwmState for which no transformation is requested
      ! back to State.  Move columns in fwmJacobian for which no transformation
      ! is requested back to Jacobian.  Transform fwmJacobian and radiances for
      ! quantities for which transformation is requested.
      do k = 1, size(fwdModelIn%quantities)
        t_qty => fwdModelIn%quantities(k)
        select case ( t_qty%template%quantityType )
        case ( l_MIFExtinction, l_MIFExtinctionV2 )
          qty => getVectorQuantityByType ( & ! Transformed quantity
              &    fwmState, quantityType=t_qty%template%quantityType, &
              &    radiometer=t_qty%template%radiometer )
          call transform_FWM_extinction ( config, &
            & fmStat%MAF, fwdModelOut, qty, t_qty, fwmJacobian, Jacobian )
        case default ! Just move the quantity and column
          qty => getVectorQuantityByType ( fwmState, &
              &  quantityType=t_qty%template%quantityType )
          call moveVectorQuantity ( qty, t_qty )
          call move_Jacobian_Columns ( qty, fwmJacobian, Jacobian )
        end select
      end do

    else ! Run the forward model without transformation
      call doForwardModels ( fwdModelIn, Jacobian )
    end if

    call MLSMessageCalls( 'pop' ) ! for all the cases

    radianceModel = any ( config%fwmType == &
      & (/ l_full, l_linear, l_polarLinear, l_hybrid, l_cloudFull /) )
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
    if ( switchDetail(switches,'FMNAN') > -1 ) then
      k = 3 ! Check, print and stop if any
    else if ( switchDetail(switches,'fmNaN') > -1 ) then
      k = 2 ! Check, print if any
    else if ( switchDetail(switches,'fmnan') > -1 ) then
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
      
      ! Check Hessians if relevant
        if ( present( Hessian ) ) then
          if ( hessian%col%vec%template%name /= jacobian%col%vec%template%name &
            & .or. hessian%row%vec%template%name /= jacobian%row%vec%template%name &
            & .or. (hessian%col%instFirst .neqv. jacobian%col%instFirst) &
            & .or. (hessian%row%instFirst .neqv. jacobian%row%instFirst) )&
            & call MLSMessage( MLSMSG_Error, ModuleName, 'Hessian is not compatible with Jacobian' )
        end if

      else if ( present( Hessian ) ) then      
        call MLSMessage ( MLSMSG_Error, ModuleName, 'Hessian is present without Jacobian' )
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

  contains

    subroutine DoForwardModels ( FwdModelIn, Jacobian )
      type(vector_t), intent(in) :: FwdModelIn
      type(matrix_t), intent(inout), optional :: Jacobian
      select case (config%fwmType)
      case ( l_baseline )
        call MLSMessageCalls( 'push', constantName='BaselineForwardModel' )
        call BaselineForwardModel ( config, FwdModelIn, FwdModelExtra, &
          FwdModelOut, fmStat, Jacobian )
        call add_to_retrieval_timing( 'baseline' )
      case ( l_full )
        call MLSMessageCalls( 'push', constantName='FullForwardModel' )
        call FullForwardModel ( config, FwdModelIn, FwdModelExtra, &
          FwdModelOut, fmStat, Jacobian, Hessian )
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
    end subroutine DoForwardModels

  end subroutine ForwardModel

  pure integer function MINLOC_S_S ( A ) result ( MS )
    ! Return the subscript in A of the minimum absolute value, as a scalar
    real, intent(in) :: A(:)
    integer :: M(1)
    m = minloc(abs(a))
    ms = m(1)
  end function MINLOC_S_S

  pure integer function MINLOC_S_D ( A ) result ( MS )
    ! Return the subscript in A of the minimum absolute value, as a scalar
    double precision, intent(in) :: A(:)
    integer :: M(1)
    m = minloc(abs(a))
    ms = m(1)
  end function MINLOC_S_D

  subroutine Move_Jacobian_Columns ( Qty, FwmJacobian, Jacobian )
    ! Move the blocks in all rows in the columns of FwmJacobian corresponding
    ! to Qty to Jacobian.

    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use MatrixModule_0, only: Move_Block
    use MatrixModule_1, only: FindBlock, Matrix_T
    use VectorsModule, only: Vector_t, VectorValue_t

    type(vectorValue_t), intent(in) :: Qty ! Quantity in fwmState
    type(matrix_T), intent(inout) :: FwmJacobian
    type(matrix_T), intent(inout) :: Jacobian

    integer :: Chan, FCol, Inst, JCol, JRow

    do jRow = 1, jacobian%row%nb
      do inst = 1, qty%template%noInstances
        jCol = findBlock ( Jacobian%col, qty%index, inst )
        fCol = findBlock ( fwmJacobian%col, qty%index, inst )
        call move_block ( fwmJacobian%block(jRow,fCol), Jacobian%block(jRow,jCol) )
      end do
    end do

  end subroutine Move_Jacobian_Columns

  subroutine Transform_FWM_extinction ( Config, MAF, FwdModelOut, Qty, T_Qty, &
    &                                   FwmJacobian, Jacobian )
    ! Transform the forward model FwmJacobian to the retriever Jacobian
    ! for the case of qty%template%quantityType == l_extinction[v2].

    ! See Equations (2) and (3) in wvs-107.

    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use MatrixModule_0, only: DestroyBlock
    use MatrixModule_1, only: FindBlock, GetFullBlock, Matrix_T
    use MLSKinds, only: RM, RV
    use VectorsModule, only: Vector_t, VectorValue_t

    type(forwardModelConfig_T), intent(in) :: Config
    integer, intent(in) :: MAF               ! MAjor Frame number
    type(vector_t), intent(inout) :: FwdModelOut
    type(vectorValue_t), intent(in) :: Qty   ! Transformed quantity in fwmState
    type(vectorValue_t), intent(in) :: T_Qty ! Quantity in State
    type(matrix_T), intent(in) :: FwmJacobian
    type(matrix_T), intent(inout) :: Jacobian

    integer :: Chan    ! c in wvs-107
    integer :: CZ      ! ci, channels X zetas, in wvs-107
    integer :: FCols(qty%template%noInstances) ! of FwmJacobian, j in wvs-107
    integer :: Inst    ! Loop inductor to compute fCols, j in wvs-107
    integer :: JCol    ! of Jacobian
    integer :: JCols(t_qty%template%noInstances) ! of Jacobian, n in wvs-107
    integer :: JRow    ! of both Jacobians, n in wvs-107
    integer :: NChans  ! Number of channels, N_C in wvs-107
    integer :: NSurfs  ! Number of surfaces in MIF, N_m in wvs-107
    real(rv) :: P      ! 10**(-2*(zeta(surf)-zeta(jRow)))
    real(rm) :: RowSum ! sum(FwmJacobian%block(jRow,fCols)%values(cz,surf))
    integer :: Surf    ! column of FwmJacobian%block(jRow,fCol)%values, 
                       ! g in wvs-107
    integer :: VSurf   ! Surface (zeta) index in a MIF, i in wvs-107

    ! Get the block column subscripts for instances of qty.
    do inst = 1, qty%template%noInstances
      fCols(inst) = findBlock ( fwmJacobian%col, qty%index, inst )
    end do

    do inst = 1, t_qty%template%noInstances
      jCols(inst) = findBlock ( Jacobian%col, t_qty%index, inst )
    end do

    nChans = size(config%channels)
    do jRow = 1, jacobian%row%nb
      nSurfs = fwdModelOut%quantities(jRow)%template%noSurfs
      do jCol = 1, size(jCols)
        if ( Jacobian%row%inst(jRow) /= Jacobian%col%inst(jCols(jCol)) ) then
          call destroyBlock ( Jacobian%block(jRow,jCols(jCol)) )
        else ! (jRow,jCol) = nn in wvs-107
          call getFullBlock ( Jacobian, jrow, jCols(jCol), 'Transform_FWM_extinction' )
          do chan = 1, nChans
            do vSurf = 1, nSurfs ! Same for all fCols
              cz = config%channels(chan)%used + 1 - & ! ci, Channel and Zeta
                 & config%channels(chan)%origin + nChans*(vSurf-1)
              do surf = 1, size(FwmJacobian%block(jRow,fCols(inst))%values,2) ! g
                ! Compute the inner sum in Equations (3) and (4) in wvs-107. 
                ! But we can't do
                ! rowSum = sum(FwmJacobian%block(jRow,fCols)%values(cz,surf))
                ! because "values" has the POINTER attribute.
                rowSum = 0.0
                do inst = 1, size(fCols)
                  rowSum = rowSum + &
                    & FwmJacobian%block(jRow,fCols(inst))%values(cz,surf)
                end do

                p = 10.0_rm ** ( -2.0_rm * ( qty%template%surfs(surf,1) - &
                                             t_qty%template%surfs(jRow,maf) ) )

                ! Equation (3) in wvs-107
                Jacobian%block(jRow,jCols(jCol))%values(cz,1) = rowSum * p

                ! Equation (4) in wvs-107
                fwdModelOut%quantities(jRow)%values(cz,maf) = &
                  & fwdModelOut%quantities(jRow)%values(cz,maf) + &
                    & rowSum * ( t_qty%values(vSurf,maf) * p - &
                               & qty%values(surf,maf) )
              end do
            end do
          end do
        end if
      end do
    end do

  end subroutine Transform_FWM_extinction

  subroutine Transform_MIF_Extinction ( MAF, T_Qty, Lrp, DumpTransform, Qty )

    ! Interpolate from a MAF-indexed minor-frame quantity to the first
    ! column of an L2GP quantity.  Then spread it out in orbit geodetic
    ! angle as if it were a 2D quantity.

    ! Equation (1) in wvs-107.

    use MLSKinds, only: RV
    use MLSNumerics, only: InterpolateValues
    use Sort_m, only: Sortp
    use VectorsModule, only: Dump, VectorValue_t

    integer, intent(in) :: MAF                ! MAjor Frame number
    type(vectorValue_t), intent(in) :: T_Qty  ! State quantity to be transformed
    real(rv), intent(in) :: LRP               ! Lowest retrieved pressure
    integer, intent(in) :: DumpTransform      ! Dump T_Qty and Qty at the end
    type(vectorValue_t), intent(inout) :: Qty ! Forward model quantity

    integer :: I       ! Subscript and loop inductor
    integer :: P(t_qty%template%noSurfs)
    integer :: Q_LRP   ! Index in Qty%Surfs of LRP
    integer :: T_LRP   ! Index in T_Qty%Surfs of LRP

    interface MINLOC_S
      module procedure MINLOC_S_S, MINLOC_S_D
    end interface

    ! Interpolate in Zeta only, from t_qty to the first column of qty
    ! PTan might be out of order, so sort t_qty%template%surfs
    call sortp ( t_qty%template%surfs(:,maf), 1, t_qty%template%noSurfs, p )
    call interpolateValues (                              &
      & t_qty%template%surfs(p,maf), t_qty%values(p,maf), &
      & qty%template%surfs(:,1),   qty%values(:,1), 'L', 'C' )
    q_lrp = minloc_s(lrp - qty%template%surfs(:,1))
    t_lrp = minloc_s(lrp - t_qty%template%surfs(:,1))
    ! Apply a P**(-2) dependence, Equation (2) in wvs-107
    do i = 1, q_lrp-1
      ! hGrids of qty and t_qty have same extent and spacing.
      ! Vertical coordinate is zeta, so 10**(-2*zeta) is P**(-2).
      qty%values(i,:) = t_qty%values(t_lrp,maf) * &
        & 10 ** ( - 2.0 * ( qty%template%surfs(i,1)) - &
                          & t_qty%template%surfs(t_lrp,1) )              
    end do
    do i = q_lrp, qty%template%noSurfs
      qty%values(i,2:) = qty%values(i,1)
    end do
    if ( dumpTransform > -1 ) then
      call dump ( t_qty, details=dumpTransform, name='from fwdModelIn' )
      call dump ( qty, details=dumpTransform, name='to fwmState' )
    end if
  end subroutine Transform_MIF_Extinction

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
! Revision 2.35  2011/12/21 01:42:22  vsnyder
! Add MIFExtinction transformation
!
! Revision 2.34  2011/05/09 18:10:02  pwagner
! Converted to using switchDetail
!
! Revision 2.33  2010/08/27 06:25:38  yanovsky
! ForwardModel subroutine has Hessian as dummy argument.
! Actual argument Hessian is passed in a call to FullForwardModel subroutine.
!
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
