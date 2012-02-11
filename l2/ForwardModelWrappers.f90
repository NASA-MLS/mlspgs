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
    use FORWARDMODELCONFIG, only: DERIVEFROMFORWARDMODELCONFIG, &
      & DESTROYFORWARDMODELDERIVED, FORWARDMODELCONFIG_T
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELINTERMEDIATE, only: FORWARDMODELSTATUS_T
    use FULLCLOUDFORWARDMODEL, only: FULLCLOUDFORWARDMODELWRAPPER
    use FULLFORWARDMODEL_M, only: FULLFORWARDMODEL
    use HESSIANMODULE_1, only: HESSIAN_T
    use HYBRIDFORWARDMODEL_M, only: HYBRIDFORWARDMODEL
    use INIT_TABLES_MODULE, only: L_BASELINE, L_CLOUDFULL, &
      & L_Extinction, L_ExtinctionV2, L_FULL, L_HYBRID, L_LINEAR, &
      & L_MIFExtinction, L_MIFExtinctionV2, L_POLARLINEAR, L_Ptan, L_SCAN, &
      & L_SCAN2D, L_SWITCHINGMIRROR
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
    integer :: DumpTransform              ! Dump transformed vector quantities
    integer :: FCol, FRow                 ! Column, Row of block in fwmJacobian
    type(vectorValue_t), pointer :: F_Qty ! Transformed FWM quantity in fwmState
    integer :: I, K
    integer :: Inst                       ! Instance index
    integer :: JRow, JCol                 ! Row, Column of block in Jacobian
    type(vectorValue_t), pointer :: LRP   ! Lowest Retrieved Pressure
    type(vectorValue_t), pointer :: Ptan  ! Tangent pressure
    logical :: RadianceModel
    type(vectorValue_t), pointer :: S_Qty ! MIF-basis Quantity in State
    type(vector_t), pointer :: TheState
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

    call deriveFromForwardModelConfig ( config )

    if ( present(fwmJacobian) ) then ! Need transformations
      ! Lowest Returned Pressure is needed for extinction transformations
      lrp => getVectorQuantityByType ( fwdModelExtra, &
               & quantityType=l_lowestRetrievedPressure )
      ptan => getVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
               & quantityType=l_ptan )
      ! Move quantities in the retriever's state vector (FwdModelIn) for which
      ! no transformation is requested to FwmState.  Transform those for which
      ! transformation is requested. We move the quantities so as to provide a
      ! consistent view of the state to the Forward model.  The un-transformed
      ! quantities in FwdModelIn are moved instead of copied to avoid aliasing
      ! and memory leaks.
      do k = 1, size(fwdModelIn%quantities)
        s_qty => fwdModelIn%quantities(k)    ! State quantity to be transformed
        select case ( s_qty%template%quantityType )
        case ( l_MIFExtinction )
          f_qty => getVectorQuantityByType ( & ! Transformed quantity
            &    fwmState, quantityType=l_vmr, molecule=l_extinction, &
            &    radiometer=s_qty%template%radiometer )
          call transform_MIF_extinction &
            & ( fmStat%MAF, s_qty, ptan, lrp%values(1,1), dumpTransform, f_qty )
        case ( l_MIFExtinctionV2 )
          f_qty => getVectorQuantityByType ( & ! Transformed quantity
            &    fwmState, quantityType=l_vmr, molecule=l_extinctionV2, &
            &    radiometer=s_qty%template%radiometer )
          call transform_MIF_extinction &
            & ( fmStat%MAF, s_qty, ptan, lrp%values(1,1), dumpTransform, f_qty )
        case default ! Just move the quantity
          f_qty => getVectorQuantityByType ( fwmState, &
              &  quantityType=s_qty%template%quantityType )
          call moveVectorQuantity ( s_qty, f_qty )
        end select
      end do

      ! Run the forward model with the transformed quantities
      call doForwardModels ( fwmState, fwmJacobian )

      ! Move quantities in FwmState for which no transformation is requested
      ! back to State.  Move columns in fwmJacobian for which no transformation
      ! is requested back to Jacobian.  Transform fwmJacobian and radiances for
      ! quantities for which transformation is requested.
      do k = 1, size(fwdModelIn%quantities)
        s_qty => fwdModelIn%quantities(k)
        select case ( s_qty%template%quantityType )
        case ( l_MIFExtinction )
          f_qty => getVectorQuantityByType ( & ! Transformed quantity
              &    fwmState, quantityType=l_vmr, molecule=l_extinction, &
              &    radiometer=s_qty%template%radiometer )
          call transform_FWM_extinction ( config, &
            & fmStat%MAF, fwdModelOut, f_qty, s_qty, ptan, &
            & fwmJacobian, Jacobian, dumpTransform )
        case ( l_MIFExtinctionV2 )
          f_qty => getVectorQuantityByType ( & ! Transformed quantity
              &    fwmState, quantityType=l_vmr, molecule=l_extinctionV2, &
              &    radiometer=s_qty%template%radiometer )
          call transform_FWM_extinction ( config, &
            & fmStat%MAF, fwdModelOut, f_qty, s_qty, ptan, &
            & fwmJacobian, Jacobian, dumpTransform )
        case default ! Just move the quantity and column
          f_qty => getVectorQuantityByType ( fwmState, &
              &  quantityType=s_qty%template%quantityType )
          call moveVectorQuantity ( f_qty, s_qty )
          call move_Jacobian_Columns ( f_qty, fwmJacobian, Jacobian )
        end select
      end do

    else ! Run the forward model without transformation
      call doForwardModels ( fwdModelIn, Jacobian )
    end if

    call destroyForwardModelDerived ( config )

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

  subroutine Move_Jacobian_Columns ( F_Qty, FwmJacobian, Jacobian )
    ! Move the blocks in all rows in the columns of FwmJacobian corresponding
    ! to F_Qty to Jacobian.

    use MatrixModule_0, only: Move_Block
    use MatrixModule_1, only: FindBlock, Matrix_T
    use VectorsModule, only: Vector_t, VectorValue_t

    type(vectorValue_t), intent(in) :: F_Qty ! FWM Quantity in fwmState
    type(matrix_T), intent(inout) :: FwmJacobian
    type(matrix_T), intent(inout) :: Jacobian

    integer :: Chan, FCol, Inst, JCol, JRow

    do jRow = 1, jacobian%row%nb
      do inst = 1, f_qty%template%noInstances
        jCol = findBlock ( Jacobian%col, f_qty%index, inst )
        fCol = findBlock ( fwmJacobian%col, f_qty%index, inst )
        call move_block ( fwmJacobian%block(jRow,fCol), Jacobian%block(jRow,jCol) )
      end do
    end do

  end subroutine Move_Jacobian_Columns

  subroutine Transform_FWM_extinction ( Config, MAF, FwdModelOut, F_Qty, S_Qty, &
    &                                   PTan, FwmJacobian, Jacobian, DumpTransform )
    ! Transform the forward model FwmJacobian to the retriever Jacobian
    ! for the case of f_qty%template%quantityType == l_extinction[v2].
    ! Transform the forward model radiances to retriever radiances.

    ! See Equations (3) and (4) in wvs-107.

    use Dump_0, only: Dump
    use ForwardModelConfig, only: ForwardModelConfig_T
    use MatrixModule_0, only: DestroyBlock, Dump, M_Banded
    use MatrixModule_1, only: FindBlock, CreateBlock, Dump, Matrix_T
    use MLSKinds, only: RM, RP, RV
    use VectorsModule, only: Dump, Vector_t, VectorValue_t

    type(forwardModelConfig_T), intent(in) :: Config
    integer, intent(in) :: MAF               ! MAjor Frame number
    type(vector_t), intent(inout) :: FwdModelOut
    type(vectorValue_t), intent(in) :: F_Qty ! Transformed FWM quantity in fwmState
    type(vectorValue_t), intent(in) :: S_Qty ! Quantity in State
    type(vectorValue_t), intent(in) :: PTan  ! Tangent pressure, for S_Qty
    type(matrix_T), intent(in) :: FwmJacobian
    type(matrix_T), intent(inout) :: Jacobian
    integer, intent(in) :: DumpTransform

    integer :: Chan    ! c in wvs-107
    integer :: CG      ! cg, channels X zetas, in wvs-107
    integer :: CV      ! Subscript corresponding to Chan
    integer :: CZ      ! ci, channels X zetas, in wvs-107
    integer :: FCols(f_qty%template%noInstances) ! of FwmJacobian, j in wvs-107
    integer :: Inst    ! Loop inductor to compute fCols, j in wvs-107
    integer :: JCol    ! of Jacobian
    integer :: JCols(s_qty%template%noInstances) ! of Jacobian, n in wvs-107
    integer :: JRow    ! of both Jacobians, n in wvs-107
    integer :: NConfigChans  ! Number of channels in configuration, N_C in wvs-107
    integer :: NVecChans ! Number of channels in radiance
    integer :: NSurfs  ! Number of surfaces in MIF, N_m in wvs-107
    real(rv) :: P      ! 10**(-2*(zeta(surf)-zeta(vSurf)))
    real(rm) :: RowSum ! sum(FwmJacobian%block(jRow,fCols)%values(cz,surf))
    integer :: Surf    ! column of FwmJacobian%block(jRow,fCol)%values, 
                       ! g in wvs-107
    integer :: VSurf   ! Surface (zeta) index in a MIF, i in wvs-107

    ! Get the block column subscripts for instances of f_qty.
    do inst = 1, f_qty%template%noInstances
      fCols(inst) = findBlock ( fwmJacobian%col, f_qty%index, inst )
    end do

    do inst = 1, s_qty%template%noInstances
      jCols(inst) = findBlock ( Jacobian%col, s_qty%index, inst )
    end do

    if ( dumpTransform > -1 ) then
      call dump ( fwdModelOut, 1, 'fwdModelOut before transformation' )
      call dump ( pTan%values, name='PTan zetas' )
    end if
    nConfigChans = size(config%channels)
    do jRow = 1, Jacobian%row%nb
      nSurfs = fwdModelOut%quantities(jRow)%template%noSurfs
      nVecChans = fwdModelOut%quantities(jRow)%template%noChans
      do jCol = 1, size(jCols)
        if ( Jacobian%row%inst(jRow) /= Jacobian%col%inst(jCols(jCol)) ) then
          ! Zero off-diagonal extinction blocks of retriever Jacobian
          call destroyBlock ( Jacobian%block(jRow,jCols(jCol)) )
        else ! (inst(jRow),inst(jCol)) = nn in wvs-107
          call createBlock ( Jacobian, jrow, jCols(jCol), m_banded, &
            & FwmJacobian%row%nelts(jRow), bandHeight=nVecChans, &
            & forWhom='Transform_FWM_extinction' )
          Jacobian%block(jRow,jCols(jCol))%values = 0
          do vSurf = 1, nSurfs ! i in wvs-107, Same for all fCols
            fwdModelOut%quantities(jRow)%values(vSurf,maf) = 0
            ! Jacobian is banded.  The only nonzeros in a column are in
            ! rows for the same MIF as the column, and the maximum number of
            ! nonzeros is the number of channels in the band.  This could be
            ! sharpened to the range from the smallest to largest channel
            ! numbers, but the computations would be messier.
            do chan = 1, nConfigChans
              cv = config%channels(chan)%used + 1 - config%channels(chan)%origin
              cz = cv + nVecChans*(vSurf-1) ! ci in wvs-107
              ! The surf loop runs in reverse order to sum small-to-large
              do surf = f_qty%template%noSurfs, 1, -1 ! g in wvs-107
                rowSum = 0.0
                do inst = 1, size(fCols)
    !{ Compute the inner sum in Equations (3) and (4) in wvs-107.
    ! \begin{equation*}
    !  \sum_{j=1}^{N_\zeta} K^F_{nj,cig}
    ! \end{equation*}
    ! where $n$ is the Jacobian block row number,
    ! $j$ is the Jacobian block column number,
    ! $c$ is the channel number,
    ! $i$ is the MIF number, and
    ! $g$ is the surface ($\zeta$) number.

    ! But we can't do
    ! rowSum = sum(FwmJacobian%block(jRow,fCols)%values(cz,surf))
    ! because "values" has the POINTER attribute.
                  rowSum = rowSum + &
                    & FwmJacobian%block(jRow,fCols(inst))%values(cz,surf)
print '(4(a,i0),1p,2(a,g14.6))', &
'FwmJacobian%block(',jrow,',',fcols(inst),')%values(',cz,',',surf,') =', &
FwmJacobian%block(jRow,fCols(inst))%values(cz,surf),', rowSum =', rowSum
                end do ! inst

                cg = cv + nVecChans*(surf-1) ! cg in wvs-107
                p = 10.0_rm ** ( -2.0_rm * ( f_qty%template%surfs(surf,1) - &
                                             ptan%values(vSurf,maf) ) )
print '(a,i0,1p,a,g14.6,2(a,i0),2(a,g14.6))',&
'f_qty%template%surfs(',surf,',1) =', f_qty%template%surfs(surf,1),&
'ptan%values(',vSurf,',',maf,') =', ptan%values(vSurf,maf), ', p =', p
    !{Compute the retriever's Jacobian using Equation (3) in wvs-107
    ! \begin{equation*}
    ! \begin{array}{lll}
    !  K^R_{nn,cii} = \sum_{g=1}^{N_\zeta}
    !   \left( \sum_{j=1}^{N_\phi} K^F_{nj,cig} \right)
    !   10^{-2(\zeta_g - \zeta_i)} & i = 1, \dots, N_m & c = 1, \dots, N_C\,, \\
    ! \end{array}
    ! \end{equation*}
    ! where the subscripts are as above.

                Jacobian%block(jRow,jCols(jCol))%values(cz,1) = &
                  & Jacobian%block(jRow,jCols(jCol))%values(cz,1) + rowSum * p
print '(1p,a,g14.6,3(a,i0),a,g14.6)', &
'rowSum*p =', rowSum*p, &
', Jacobian%block(',jrow,',',jCols(jCol),')%values(',cz,',1) =', &
Jacobian%block(jRow,jCols(jCol))%values(cz,1)
    !{Compute radiance using Equation (4) in wvs-107
    ! \begin{equation*}
    ! \begin{array}{lll}
    ! I^R_{ci,n} = I^F_{ci,n} + \sum_{g=1}^{N_\zeta}
    !  \left( \sum_{j=1}^{N_\phi} K^F_{nj,cig} \right)
    !  \left( E^R_{i,n} 10^{-2(\zeta_g - \zeta_i)} -
    !   E^F_{g,1} \right) & i = 1, \dots, N_m & c = 1, \dots, N_C \,, \\
    ! \end{array}
    ! \end{equation*}
    ! where the subscripts are as above.

                fwdModelOut%quantities(jRow)%values(vSurf,maf) = &
                  & fwdModelOut%quantities(jRow)%values(vSurf,maf) + &
                    & rowSum * ( s_qty%values(vSurf,maf) * p - &
                               & f_qty%values(surf,1) )
              end do ! surf
            end do ! chan
          end do ! vSurf
          if ( dumpTransform > -1 ) then
            do inst = 1, size(fCols)
              call dump ( fwmJacobian, details=9, name='Forward model Jacobian', &
                        & row=jrow, column=fCols(inst) )
            end do
            call dump ( Jacobian, details=9, name='Transformed Jacobian', &
                      & row=jrow, column=jCols(jCol) )
          end if
        end if
      end do ! jCol
    end do ! jRow
    if ( dumpTransform > -1 ) &
      & call dump ( fwdModelOut, 1, 'fwdModelOut after transformation' )

  end subroutine Transform_FWM_extinction

  subroutine Transform_MIF_Extinction ( MAF, S_Qty, PTan, Lrp, DumpTransform, &
    &                                   F_Qty )

    ! Interpolate from a MAF-indexed minor-frame quantity to the first
    ! column of an L2GP quantity.  Then spread it out in orbit geodetic
    ! angle as if it were a 2D quantity.

    ! Equations (1) and (2) in wvs-107.

    use MLSKinds, only: RV
    use MLSNumerics, only: InterpolateValues
    use Sort_m, only: Sortp
    use VectorsModule, only: Dump, VectorValue_t

    integer, intent(in) :: MAF                ! MAjor Frame number
    type(vectorValue_t), intent(in) :: S_Qty  ! State quantity to be transformed
    type(vectorValue_t), intent(in) :: PTan   ! Zetas, for S_Qty, which has
                                              ! altitudes in meters
    real(rv), intent(in) :: LRP               ! Lowest retrieved pressure
    integer, intent(in) :: DumpTransform      ! Dump S_Qty and F_Qty at the end
    type(vectorValue_t), intent(inout) :: F_Qty ! Forward model quantity

    integer :: F_LRP   ! Index in F_Qty%Surfs of LRP
    integer :: I       ! Subscript and loop inductor
    integer :: P(s_qty%template%noSurfs)
    integer :: S_LRP   ! Index in S_Qty%Surfs of LRP

    interface MINLOC_S
      module procedure MINLOC_S_S, MINLOC_S_D
    end interface

    ! Interpolate in Zeta only, from S_Qty to the first column of F_Qty
    ! PTan might be out of order, so sort ptan%values
    call sortp ( ptan%values(:,maf), 1, ptan%template%noSurfs, p )
    call interpolateValues (                              &
      & ptan%values(p,maf), s_qty%values(p,maf), &
      & f_qty%template%surfs(:,1), f_qty%values(:,1), 'L', 'C' )
    f_lrp = minloc_s(lrp - f_qty%template%surfs(:,1))
    s_lrp = minloc_s(lrp - ptan%values(:,1))
     !{ Apply a $P^{-2}$ dependence, Equation (2) in wvs-107,
     !  to compute extinction for the forward model
     ! \renewcommand{\arraystretch}{2}
     ! \begin{equation*}
     ! E^F_{g,j} = \left\{
     ! \begin{array}{llll}
     !  \overline{E}^R_{n}(\zeta_g)
     !   & \zeta_g \geq \zeta_r & j = 1, \dots, N_\phi & g = 1, \dots, N_\zeta \\
     !  \overline{E}^R_{n}(\zeta_r) \times 10^{-2(\zeta_g - \zeta_r)}
     !   & \zeta_g < \zeta_r    & j = 1, \dots, N_\phi & g = 1, \dots, N_\zeta \\
     ! \end{array} \right.
     ! \end{equation*}
     ! {\tt hGrids} of {\tt F_Qty} and {\tt S_Qty} have same extent and spacing.
     ! Vertical coordinate is $\zeta = \log_{10}P$,
     ! so $10^{-2\zeta}$ is $P^{-2}$).

    do i = 1, f_lrp-1
      f_qty%values(i,:) = s_qty%values(s_lrp,maf) * &
        & 10 ** ( - 2.0 * ( f_qty%template%surfs(i,1)) - &
                          & ptan%values(s_lrp,1) )              
    end do
    do i = f_lrp, f_qty%template%noSurfs
      f_qty%values(i,2:) = f_qty%values(i,1)
    end do
    if ( dumpTransform > -1 ) then
      call dump ( s_qty, details=dumpTransform, name='from fwdModelIn' )
      call dump ( f_qty, details=dumpTransform, name='to fwmState' )
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
! Revision 2.36  2012/02/11 21:28:31  vsnyder
! Interim MIF extinction commit
!
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
