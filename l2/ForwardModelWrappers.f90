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
    & FwdModelOut, fmStat, Jacobian, Hessian, Vectors )

    ! Call the forward model selected by Config.

    use BASELINEFORWARDMODEL_M, only: BASELINEFORWARDMODEL
    use FORWARDMODELCONFIG, only: DERIVEFROMFORWARDMODELCONFIG, &
      & DESTROYFORWARDMODELDERIVED, FORWARDMODELCONFIG_T, QtyStuff_t
    use FORWARDMODELINTERMEDIATE, only: FORWARDMODELSTATUS_T
    use ForwardModelVectorTools, only: GETQUANTITYFORFORWARDMODEL
    use FULLCLOUDFORWARDMODEL, only: FULLCLOUDFORWARDMODELWRAPPER
    use FULLFORWARDMODEL_M, only: FULLFORWARDMODEL
    use HESSIANMODULE_1, only: HESSIAN_T
    use HYBRIDFORWARDMODEL_M, only: HYBRIDFORWARDMODEL
    use INIT_TABLES_MODULE, only: L_BASELINE, L_CLOUDFULL, &
      & L_EXTINCTION, L_EXTINCTIONV2, L_FULL, L_HYBRID, L_LINEAR, &
      & L_MIFEXTINCTION, L_MIFEXTINCTIONV2, L_POLARLINEAR, L_PTAN, L_SCAN, &
      & L_SCAN2D, L_SWITCHINGMIRROR
    use Intrinsic, only: L_LOWESTRETRIEVEDPRESSURE, L_VMR
    use LINEARIZEDFORWARDMODEL_M, only: LINEARIZEDFORWARDMODEL
    use MATRIXMODULE_1, only: CHECKINTEGRITY, MATRIX_T
    use MLSL2OPTIONS, only: MLSMESSAGE
    use MLSL2TIMINGS, only: ADD_TO_RETRIEVAL_TIMING
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use Molecules, only: L_EXTINCTION, L_EXTINCTIONV2
    use POLARLINEARMODEL_M, only: POLARLINEARMODEL
    use SCANMODELMODULE, only: SCANFORWARDMODEL, TWODSCANFORWARDMODEL
    use STRING_TABLE, only: DISPLAY_STRING, GET_STRING
    use SWITCHINGMIRRORMODEL_M, only: SWITCHINGMIRRORMODEL
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: EMIT, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VECTORSMODULE, only: CHECKNAN, DUMP, VECTOR_T, VECTORVALUE_T

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: CONFIG
    type(vector_T), intent(inout), target :: FWDMODELIN ! The state
    type(vector_T), intent(in) :: FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional, target :: JACOBIAN ! The
      ! Jacobian, used for the Newton iteration.
    type(hessian_T), intent(inout), optional :: HESSIAN ! No transformation
    type(vector_t), dimension(:), target, optional :: VECTORS ! Vectors database

    ! Local variables
    logical :: Clean                      ! Dumps are clean, from switch dxfc
    real :: DeltaTime
    integer :: DumpTransform(3)           ! Dump levels for transformed stuff
                                          ! 1 = 1's digit => Input and output
                                          ! 2 = 10's digit => FWM Jacobian
                                          ! 3 = 100's digit => Transformed Jacobian
    integer :: FMNaN                      ! Level of fmnan switch
    integer :: I, K
    type(vectorValue_t), pointer :: LRP   ! Lowest Retrieved Pressure
    integer :: NQty                       ! Number of quantities to transform
    type(vectorValue_t), pointer :: Ptan  ! Tangent pressure
    type(qtyStuff_t) :: &                 ! Pointers to MIF extinction quantities
      & Qtys(count( &                     ! (:,1) are MIF-basis, (:,2) are gridded
          &   fwdModelIn%quantities%template%quantityType == L_MIFExtinction .or. &
        &     fwdModelIn%quantities%template%quantityType == L_MIFExtinctionV2 &
        & ) ,2)
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
    if ( dumpTransform(1) >= 0 ) then
      dumpTransform(3) = dumpTransform(3)/100
      dumpTransform(2) = mod(dumpTransform(2)/10,10)
      dumpTransform(1) = mod(dumpTransform(1),10)
    end if
    clean = switchDetail(switches,'dxfc') > -1

    ! Do the actual forward models

    if ( config%fwmType == l_full) call deriveFromForwardModelConfig ( config )

    if ( config%transformMIFextinction ) then
      ! Transform MIFextinction quantities to extinction molecules before
      ! calling the forward model, and transform extinction molecules and
      ! associated columns of the Jacobian to MIF extinction quantities after
      ! return.  See wvs-107.
      if ( size(qtys) == 0 ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Forward model config ' // trim(thisName) // &
          & ' requests MIF extinction transformation, but there are no' // &
          & ' MIF extinction quantities' )
      ! Find MIF extinction quantities in the state vector, and their
      ! corresponding extinction quantities
      nQty = 0
      do i = 1, size(fwdModelIn%quantities)
        if ( fwdModelIn%quantities(i)%template%radiometer /= &
           & config%signals(1)%radiometer ) cycle
        if ( fwdModelIn%quantities(i)%template%quantityType == l_MIFextinction ) then
          nQty = nQty + 1
          qtys(nQty,1)%qty => fwdModelIn%quantities(i)
          qtys(nQty,2)%qty => GetQuantityForForwardModel ( fwdModelIn, &
               & quantityType=l_vmr, molecule=l_extinction, &
               & radiometer=fwdModelIn%quantities(i)%template%radiometer )
        else if( fwdModelIn%quantities(i)%template%quantityType == l_MIFextinctionv2 ) then
          nQty = nQty + 1
          qtys(nQty,1)%qty => fwdModelIn%quantities(i)
          qtys(nQty,2)%qty => GetQuantityForForwardModel ( fwdModelIn, &
               & quantityType=l_vmr, molecule=l_extinctionv2, &
               & radiometer=fwdModelIn%quantities(i)%template%radiometer )
        end if
      end do
      ! Lowest Returned Pressure is needed for extinction transformations
      lrp => GetQuantityForForwardModel ( fwdModelExtra, &
               & quantityType=l_lowestRetrievedPressure )
      !??? Future upgrade ???:  Use the mask field on MIFExtinction
      !??? to determine the lowest pressure.
      ptan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
               & quantityType=l_ptan, &
               & instrumentModule=config%signals(1)%instrumentModule )
      ! Transform MIF extinction quantities to extinction molecules
      do k = 1, nQty
        call transform_MIF_extinction &
          & ( fmStat%MAF, qtys(k,1)%qty, ptan, lrp%values(1,1), qtys(k,2)%qty, &
            & dumpTransform )
      end do

      ! Run the forward model with the transformed quantities
      call doForwardModels

      ! Transform fwmJacobian and radiances for quantities for which
      ! transformation is requested.
      do k = 1, nQty
        call transform_FWM_extinction ( config, fmStat%MAF, fwdModelOut, &
          & qtys(k,2)%qty, qtys(k,1)%qty, ptan, Jacobian, dumpTransform, clean )
      end do

    else ! Run the forward model without transformation
      call doForwardModels
    end if

    call destroyForwardModelDerived ( config )

    fmnan = switchDetail(switches,'fmnan')
    if ( fmnan > 0 ) then
      ! Check radiances
      if ( checkNaN(fwdModelOut, k-1, 'ForwardModelOut') ) then
        if ( fmnan > 1 ) then
          call dump ( fwdModelIn, k-1, 'ForwardModelIn' )
          call dump ( fwdModelExtra, k-1, 'ForwardModelExtra' )
        end if
        call display_string ( config%name, &
          & before='Forward model config name: ', advance='yes' )
        call MLSMessage ( merge(MLSMSG_Error,MLSMSG_Warning,fmnan>2), &
          & ModuleName, 'NaNs found in forward model output' )
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

    subroutine DoForwardModels

      select case (config%fwmType)
      case ( l_baseline )
        call MLSMessageCalls( 'push', constantName='BaselineForwardModel' )
        call BaselineForwardModel ( config, FwdModelIn, FwdModelExtra, &
          FwdModelOut, fmStat, Jacobian )
        call add_to_retrieval_timing( 'baseline' )
      case ( l_cloudFull )
        call MLSMessageCalls( 'push', constantName='CloudForwardModel' )
        call FullCloudForwardModelWrapper ( config, FwdModelIn, FwdModelExtra, &
          FwdModelOut, fmStat, Jacobian )
        call add_to_retrieval_timing( 'fullcloud_fwm' )
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
      case ( l_switchingMirror )
        call MLSMessageCalls( 'push', constantName='SwitchingForwardModel' )
        call SwitchingMirrorModel ( config, FwdModelIn, FwdModelExtra, &
          FwdModelOut, fmStat, Jacobian )
        call add_to_retrieval_timing( 'switching_mirror' )
      case default ! Shouldn't get here if parser etc. worked
      end select

      if ( config%isRadianceModel ) then
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

      call MLSMessageCalls( 'pop' ) ! for all the cases

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

!{\cleardoublepage
  subroutine Transform_FWM_extinction ( Config, MAF, FwdModelOut, F_Qty, S_Qty, &
    &                                   PTan, Jacobian, DumpTransform, Clean )
    ! Transform blocks of the Jacobian associated with f_qty to blocks
    ! associated with s_qty.
    ! Transform the forward model radiances to retriever radiances.

    ! See Equations (3) and (4) in wvs-107.

    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELVECTORTOOLS, only: GETQUANTITYFORFORWARDMODEL
    use INTRINSIC, only: L_RADIANCE
    use MATRIXMODULE_0, only: DESTROYBLOCK, M_ABSENT, M_BANDED, M_FULL
    use MATRIXMODULE_1, only: FINDBLOCK, CREATEBLOCK, DUMP, MATRIX_T
    use MLSKINDS, only: RM, RV
    use MLSL2OPTIONS, only: MLSMESSAGE
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR
    use OUTPUT_M, only: OUTPUT
    use VECTORSMODULE, only: DUMP, M_IGNORE, M_LINALG, VECTOR_T, VECTORVALUE_T

    type(forwardModelConfig_T), intent(in) :: CONFIG
    integer, intent(in) :: MAF               ! MAjor Frame number
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(vectorValue_t), intent(in) :: F_Qty ! Extinction quantity in fwmState
    type(vectorValue_t), intent(in) :: S_Qty ! MIF Extinction quantity in State
    type(vectorValue_t), intent(in) :: PTan  ! Tangent pressure, for S_Qty
    type(matrix_T), intent(inout) :: Jacobian
    integer, intent(in) :: DumpTransform(3)
    logical, intent(in) :: Clean             ! for the dumps

    integer :: Chan    ! Index of channel in Config
    integer :: CV      ! c in wvs-107
    integer :: CZ      ! ci, channels X zetas, in wvs-107
    integer :: FCols(f_qty%template%noInstances) ! of Jacobian, j in wvs-107
    integer :: Inst    ! Loop inductor to compute fCols, j in wvs-107
    integer :: JCol    ! of Jacobian
    integer :: JCols(s_qty%template%noInstances) ! of Jacobian, n in wvs-107
    integer :: JRow    ! of Jacobian, n in wvs-107
    integer :: NVecChans ! Number of channels in radiance
    type(vectorValue_t), pointer :: O_Qty        ! Qty of output vector
    real(rv) :: P(size(ptan%values,1),f_qty%template%noSurfs) ! 10**(-2*(zeta(surf)-zeta(vSurf)))
    real(rm) :: RowSum ! sum(Jacobian%block(jRow,fCols)%values(cz,:))
    integer :: SB      ! Sideband, zero or from the first signal in config
    integer :: Surf    ! column of Jacobian%block(jRow,fCol)%values, 
                       ! g in wvs-107
    integer :: VSurf   ! Surface (zeta) index in a MIF, i in wvs-107

    sb = merge( 0, config%signals(1)%sideband, config%forceFoldedOutput )

    ! Get the block column subscripts for instances of f_qty.
    do inst = 1, f_qty%template%noInstances
      fCols(inst) = findBlock ( Jacobian%col, f_qty%index, inst )
    end do

    ! Get the block column subscripts for instances of s_qty.
    do inst = 1, s_qty%template%noInstances
      jCols(inst) = findBlock ( Jacobian%col, s_qty%index, inst )
    end do

    !{ Evaluate Equation (2) from wvs-107, \emph{viz.}
    ! \begin{equation*}
    ! 10^{-2(\zeta_g-\zeta_r)} = \left( \frac{P_g}{P_r} \right)^2
    ! \end{equation*}

    forall ( vSurf=1:size(ptan%values,1), surf=1:f_qty%template%noSurfs ) &
      & p(vSurf,surf) = &
        & 10.0_rm ** ( -2.0_rm * ( f_qty%template%surfs(surf,1) - &
                                   ptan%values(vSurf,maf) ) )

    do chan = 1, size(config%channels)
      cv = config%channels(chan)%used + 1 - config%channels(chan)%origin
      o_qty => getQuantityForForwardModel ( fwdModelOut, &
        & quantityType=l_radiance, &
        & signal=config%signals(config%channels(chan)%signal)%index, &
        & sideband=sb )
      jRow = findBlock ( Jacobian%row, o_qty%index, MAF )
      nVecChans = o_qty%template%noChans
      ! Dump the radiance, and Jacobian columns, to be transformed
      if ( dumpTransform(1) >= 1 ) then
        call output ( MAF, before='MAF ' )
        call dump ( o_qty, dumpTransform(1), &
          & ' fwdModelOut before transformation', options=merge('-c','  ',clean) )
      end if
      if ( dumpTransform(2) >= 1 ) then
        do inst = 1, size(fCols)
          call dump ( Jacobian, 'Jacobian from forward model', dumpTransform(2), &
            & row=jRow, column=fCols(inst) )
        end do
      end if
      do jCol = 1, size(jCols)
        if ( Jacobian%col%inst(jCols(jCol)) /= MAF .or. &
           & all(Jacobian%block(jRow,fCols)%kind == m_absent) ) then
          ! Zero MIF extinction blocks of Jacobian that are off-diagonal
          ! or would not be filled because all corresponding extinction
          ! blocks are absent.
          call destroyBlock ( Jacobian%block(jRow,jCols(jCol)) )
        else ! (inst(jRow),inst(jCol)) = nn in wvs-107
          select case ( Jacobian%block(jRow,jCols(jCol))%kind )
          case ( m_absent )
            call createBlock ( Jacobian, jrow, jCols(jCol), m_banded, &
              & Jacobian%row%nelts(jRow), bandHeight=nVecChans, &
              & forWhom='Transform_FWM_extinction' )
              Jacobian%block(jRow,jCols(jCol))%values = 0
          case ( m_banded )
            if ( ubound(Jacobian%block(jRow,jCols(jCol))%values,1) /= &
               & Jacobian%row%nelts(jRow) ) &
                 & call MLSMessage ( MLSMSG_Error, moduleName, &
                   & 'Band structure wrong for Transformed MIF extinction block' )
          case default
            call MLSMessage ( MLSMSG_Error, moduleName, &
              & 'Transformed MIF extinction block neither absent or banded' )
          end select
          do vSurf = 1, size(ptan%values,1) ! i in wvs-107, Same for all fCols
            ! For what rows of the Jacobian and elements of the state vector
            ! do we evaluate Equations (3) and (4)?
            ! ??? Why doesn't this work?  Perhaps because it's  ???
            ! ??? fiddling row (radiance) masks when it should  ???
            ! ??? be fiddling column (MIFExtinction) masks?     ???
!             if ( associated(s_qty%mask) ) then
!               if ( iand(ichar(s_qty%mask(vsurf,1)), m_linAlg) /= 0 ) cycle
!             end if
            ! Jacobian is banded.  The only nonzeros in a column are in
            ! rows for the same MIF as the column, and the maximum number of
            ! nonzeros is the number of channels in the band.  This could be
            ! sharpened to the range from the smallest to largest channel
            ! numbers, but the computations would be messier.
            cz = cv + nVecChans*(vSurf-1) ! ci in wvs-107
            ! The surf loop runs in reverse order to sum rowSum*p in
            ! small-to-large order
            do surf = f_qty%template%noSurfs, 1, -1 ! g in wvs-107

    !{ Compute the inner sum in Equations (3) and (4) in wvs-107.
    ! \begin{equation*}
    !  \sum_{j=1}^{N_\zeta} K^F_{nj,cig}
    ! \end{equation*}
    ! where $n$ is the Jacobian block row number,\\
    ! $j$ is the Jacobian block column number,\\
    ! $c$ is the channel number,\\
    ! $i$ is the MIF number, and\\
    ! $g$ is the surface ($\zeta$) number.

    ! But we can't do
    ! rowSum = sum(Jacobian%block(jRow,fCols)%values(cz,surf))
    ! because "values" has the POINTER attribute.

              rowSum = 0.0
              do inst = 1, size(fCols)
                select case ( Jacobian%block(jRow,fCols(inst))%kind )
                case ( m_full )
                  rowSum = rowSum + &
                    & Jacobian%block(jRow,fCols(inst))%values(cz,surf)
                case ( m_absent )
                case default
                  call MLSMessage ( MLSMSG_Error, moduleName, &
                    & "Transform_FWM_extinction cannot handle sparse or banded blocks" )
                end select
              end do ! inst

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
                & Jacobian%block(jRow,jCols(jCol))%values(cz,1) + &
                  & rowSum * p(vSurf,surf)

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

              o_qty%values(cz,maf) = o_qty%values(cz,maf) + &
                & rowSum * ( s_qty%values(vSurf,maf) * p(vSurf,surf) - &
                           & f_qty%values(surf,1) )
            end do ! Surf
          end do ! vSurf
        end if
      end do ! jCol
      if ( dumpTransform(1) >= 1 ) then
        call output ( MAF, before='MAF ' )
        call dump ( o_qty, dumpTransform(1), &
          & ' fwdModelOut after transformation', options=merge('-c','  ',clean) )
      end if
      if ( dumpTransform(3) >= 1 ) then
        do inst = 1, size(jCols)
          call dump ( Jacobian, 'Transformed Jacobian', dumpTransform(3), &
            & row=jRow, column=jCols(inst) )
        end do
      end if
    end do ! chan

    ! Destroy columns of Jacobian corresponding to f_qty so retriever has
    ! no hope of trying to retrieve extinction.
    do jCol = 1, size(fCols)
      do jRow = 1, Jacobian%row%nb
        call destroyBlock ( Jacobian%block(jRow,fCols(jCol)) )
      end do ! jRow
    end do ! jCol

    ! Make values of f_Qty zero so as not to confuse the calculation of
    ! aj%axmax in the retriever
    f_qty%values = 0.0

  end subroutine Transform_FWM_extinction

!{\cleardoublepage
  subroutine Transform_MIF_Extinction ( MAF, S_Qty, PTan, Lrp, F_Qty, &
    &                                   DumpTransform )

    ! Interpolate from a MAF-indexed minor-frame quantity to the first
    ! column of an L2GP quantity.  Then spread it out in orbit geodetic
    ! angle as if it were a 2D quantity.

    ! Equations (1) and (2) in wvs-107.

    use Dump_0, only: Dump
    use MLSKinds, only: RM, RV
    use MLSNumerics, only: InterpolateValues
    use Sort_m, only: Sortp
    use VectorsModule, only: Dump, VectorValue_t

    integer, intent(in) :: MAF                ! MAjor Frame number
    type(vectorValue_t), intent(in) :: S_Qty  ! State quantity to be transformed
    type(vectorValue_t), intent(in) :: PTan   ! Zetas, for S_Qty, which has
                                              ! altitudes in meters
    real(rv), intent(in) :: LRP               ! Lowest retrieved pressure
    type(vectorValue_t), intent(inout) :: F_Qty ! Forward model quantity
    integer, intent(in) :: DumpTransform(3)   ! Dump S_Qty and F_Qty

    integer :: F_LRP   ! Index in F_Qty%Surfs of LRP
    integer :: I       ! Subscript and loop inductor
    integer :: P(s_qty%template%noSurfs)
    integer :: S_LRP   ! Index in S_Qty%Surfs of LRP

    interface MINLOC_S
      module procedure MINLOC_S_S, MINLOC_S_D
    end interface

    ! PTan might be out of order, so sort ptan%values
    call sortp ( ptan%values(:,maf), 1, ptan%template%noSurfs, p )

    f_lrp = minloc_s(lrp - f_qty%template%surfs(:,1))
    s_lrp = minloc_s(lrp - ptan%values(p,1))

    ! Interpolate in Zeta only, from S_Qty to the first column of F_Qty
    call interpolateValues (                              &
      & ptan%values(p,maf), s_qty%values(p,maf), &
      & f_qty%template%surfs(:,1), f_qty%values(:,1), 'L', 'C' )
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
      f_qty%values(i,:) = f_qty%values(f_lrp,1) * &
        & 10.0_rm ** ( - 2.0 * ( f_qty%template%surfs(i,1) - &
                               & ptan%values(p(s_lrp),1) ) )
    end do

    do i = f_lrp, ubound(f_qty%values,1)
      f_qty%values(i,2:) = f_qty%values(i,1)
    end do

    if ( dumpTransform(1) >= 1 ) then
      call dump ( ptan%values(p,maf), name='PTan zetas' )
      call dump ( s_qty, details=dumpTransform(1), name='from fwdModelIn' )
      call dump ( f_qty, details=dumpTransform(1), name='to forward model' )
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
! Revision 2.54  2013/03/15 20:35:26  vsnyder
! Change debug print level threshold from zero to one
!
! Revision 2.53  2012/08/16 18:19:54  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.52  2012/08/14 00:38:14  vsnyder
! Specify instrument module for PTAN instead of hoping for the best
!
! Revision 2.51  2012/08/10 22:49:26  vsnyder
! Don't test linalg flag
!
! Revision 2.50  2012/07/31 00:49:33  vsnyder
! Correct testing a dump switch
!
! Revision 2.49  2012/07/06 01:54:33  vsnyder
! Transform only rows of Jacobian and radiance that are produced by current
! forward model config.  Revise some dumps.
!
! Revision 2.48  2012/07/06 00:57:20  vsnyder
! Bogus comment: ForwardModelWrappers.f90
!
! Revision 2.47  2012/07/04 02:15:01  vsnyder
! Don't compute masked rows of Jacobian.
!
! Revision 2.46  2012/06/15 23:29:39  vsnyder
! Spread extinction interpolated from MIF extinction, below the lowest
! retrieved pressure, to every profile.  Change dump switch interpretation.
!
! Revision 2.45  2012/06/07 00:47:16  vsnyder
! Handle switchDetail properly
!
! Revision 2.44  2012/06/06 20:41:59  vsnyder
! More and better dumps and messages
! don't clobber Jacobian from prior FWM for same MAF
!
! Revision 2.43  2012/05/01 22:24:06  vsnyder
! Use IsRadianceModel component, some cannonball polishing
!
! Revision 2.42  2012/04/20 01:57:15  vsnyder
! Add some dumps, add clean switch.  Handle MAF selection properly. Some
! cannonball polishing.
!
! Revision 2.41  2012/03/28 00:56:49  vsnyder
! Move check for signals with MIF extinction from Wrappers to Support
!
! Revision 2.40  2012/03/14 21:37:30  wgread
! added mif extinction transform capability vws&wgr
!
! Revision 2.39  2012/02/23 00:59:43  vsnyder
! Forgot to add RM in use for MLSKinds
!
! Revision 2.38  2012/02/23 00:08:08  vsnyder
! Maybe MIF extinction transformations work now
!
! Revision 2.37  2012/02/14 19:01:29  pwagner
! Fixed bug that broke nrt
!
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
