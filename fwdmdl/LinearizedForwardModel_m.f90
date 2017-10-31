
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module LinearizedForwardModel_m

  ! Approximate a forward model with a first-order Taylor's series.

  implicit none
  private
  public :: LinearizedForwardModel, LinearizedForwardModelAuto

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  logical, parameter :: EXPECTDELTAXALLZERO = .false.

contains ! =====     Public Procedures     =============================

  ! -------------------------------------  LinearizedForwardModel  -----
  subroutine LinearizedForwardModel ( fmConf, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, fmStat, Jacobian, vectors )

    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use FORWARDMODELINTERMEDIATE, only: FORWARDMODELSTATUS_T
    use INTRINSIC, only: L_RADIANCE, L_OPTICALDEPTH
    use MATRIXMODULE_1, only: MATRIX_T
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSSIGNALS_M, only: SIGNAL_T
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T
    use FORWARDMODELVECTORTOOLS, only: GETQUANTITYFORFORWARDMODEL

    ! Dummy arguments
    type(forwardModelConfig_T), intent(in) :: FMCONF
    ! target attribute prevents STATEQ from being undefined when
    ! FindMatchForL2PCQ returns:
    type(vector_T), intent(in), target ::  FWDMODELIN
    type(vector_T), intent(in), target ::  FWDMODELEXTRA
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_t), dimension(:), target, optional :: VECTORS ! Vectors database

    type(VectorValue_T), pointer :: RADIANCE ! The radiance quantity to fill
    integer :: SIDEBAND
    type(Signal_T), pointer :: SIGNAL   ! Signal from the configuration

    if ( size ( fmConf%signals ) /= 1 ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Can only have one signal for linearized models')

    signal => fmConf%signals(1)
    ! Probably access the sideband associated with the signal ...
    sideband = signal%sideband
    ! ... but in certain rare circumstances (when we're called by the hybrid
    ! model) we might want to force the folded one.
    if ( fmConf%forceFoldedOutput ) sideband = 0
    radiance => GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
      & signal=signal%index, sideband=sideband, noError=.true., config=fmConf )

    ! Now, it's possible we're really being asked to deal with optical depth,
    ! not radiance.
    if ( .not. associated ( radiance ) ) &
      & radiance => GetQuantityForForwardModel (fwdModelOut, quantityType=l_opticalDepth, &
        & signal=signal%index, sideband=sideband, noError=.true., config=fmConf )

    ! Now, some possible error messages
    if ( .not. associated ( radiance ) ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Unable to find a radiance or optical depth quantity for this signal' )
    if ( radiance%template%quantityType == l_opticalDepth ) then
      if ( present(jacobian)  ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Not appropriate to ask for derivatives for optical depth' )
      if ( fmConf%xStar /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Not appropriate to supply x/yStar for optical depth' )
      if  ( sideband /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Not appropriate to request optical depth from unfolded signal' )
    end if
    if ( fmConf%xStar /= 0 .and. .not. present(vectors) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'x/yStar supplied by vectors database argument not present' )

    call LinearizedForwardModelAuto ( fmConf, FwdModelIn, FwdModelExtra,&
    & fmStat, radiance, noChans=radiance%template%noChans, &
    & noMifs=radiance%template%noSurfs, &
    & noMifsJ=merge( radiance%template%noSurfs,0,present(Jacobian) ), &
    & Jacobian=Jacobian, Vectors=vectors )

  end subroutine LinearizedForwardModel

  ! ---------------------------------  LinearizedForwardModelAuto  -----
  subroutine LinearizedForwardModelAuto ( fmConf, FwdModelIn, FwdModelExtra,&
    & fmStat, radiance, noChans, noMifs, noMifsJ, Jacobian, Vectors )

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use DUMP_0, only: DUMP
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T, DUMP
    use FORWARDMODELINTERMEDIATE, only: FORWARDMODELSTATUS_T
    use FORWARDMODELVECTORTOOLS, only: GETQUANTITYFORFORWARDMODEL
    use HESSIANMODULE_1, only: MULTIPLY
    use INTRINSIC, only: L_LIMBSIDEBANDFRACTION, L_OPTICALDEPTH, &
      & L_PTAN, L_RADIANCE, L_TEMPERATURE, L_VMR, PHYQ_Angle
    use L2PC_M, only: L2PC_T, L2PCDATABASE, POPULATEL2PCBIN
    use L2PCBINS_M, only: FINDMATCHFORL2PCQ, SELECTL2PCBINS
    use MANIPULATEVECTORQUANTITIES, only: FINDONECLOSESTINSTANCE
    use MATRIXMODULE_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
      & MATRIXELEMENT_T, CREATEBLOCK, DENSIFY, CHECKFORSIMPLEBANDEDLAYOUT
    use MATRIXMODULE_1, only: MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP, &
      & FINDBLOCK, CREATEBLOCK
    use MLSKINDS, only: R8, RM, RV
    use MLSSIGNALS_M, only: SIGNAL_T
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSNUMERICS, only: Coefficients, HUNT, INTERPOLATEARRAYSETUP, &
      & INTERPOLATEARRAYTEARDOWN, INTERPOLATEVALUES
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use MOLECULES, only: ISEXTINCTION
    use MOREMESSAGE, only: MLSMESSAGE
    use OUTPUT_M, only: OUTPUT !, OUTPUTNAMEDVALUE
    use STRING_TABLE, only: DISPLAY_STRING
    use TOGGLES, only: EMIT, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VECTORSMODULE, only: ASSIGNMENT(=), OPERATOR(-), OPERATOR(+), &
      & OPERATOR(/=), CLONEVECTOR,  DESTROYVECTORINFO, DUMP, &
      & GETVECTORQUANTITYINDEXBYNAME, GETVECTORQUANTITYBYTYPE, M_LINALG, &
      & RV, VECTOR_T, VECTORVALUE_T
    use SORT_M, only: SORTP

    ! Dummy arguments
    type(forwardModelConfig_T), intent(in) :: FMCONF
    ! target attribute prevents STATEQ from being undefined when
    ! FindMatchForL2PCQ returns:
    type(vector_T), intent(in), target ::  FWDMODELIN
    type(vector_T), intent(in), target ::  FWDMODELEXTRA
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(vectorvalue_t), intent(inout) :: RADIANCE ! The radiance quantity to fill
    integer, intent(in) :: NoChans, NoMifs, NoMifsJ ! Dimensions
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_t), dimension(:), target, optional :: VECTORS ! Vectors database

    ! Local variables
    integer :: CENTER                   ! Center instance of l2pc
    integer :: CHAN                     ! Loop counter
    integer :: CLOSESTINSTANCE          ! A closest profile to this MAF
    integer :: COLJBLOCK                ! Column index in jacobian
    integer :: COLLBLOCK                ! Column index in l2pc
    integer :: DELTAINSTANCE            ! Instance offset between center an current
    integer :: I                        ! Array index
    integer :: INSTANCELEN              ! For the state quantity
    integer :: LOWER                    ! Array index
    integer :: MAF                      ! Major frame to do
    integer :: Me = -1                  ! String index for trace
    integer :: MIF                      ! Minor frame loop counter
    integer :: NOPOINTINGS              ! Number of pointings in the l2pc file
    integer :: QTYIND                   ! Loop index for main loop
    integer :: ROWJBLOCK                ! Row index in jacobian
    integer :: ROWLBLOCK                ! Row index in l2pc
    integer :: SIDEBAND                 ! Loop index
    integer :: SIDEBANDSTART            ! For sideband loop
    integer :: SIDEBANDSTEP             ! For sideband loop
    integer :: SIDEBANDSTOP             ! For sideband loop
    integer :: TOP                      ! For interpolating
    integer :: UPPER                    ! Array index
    integer :: XINSTANCE                ! Instance in x corresponding to xStarInstance
    integer :: XSTARINSTANCE            ! Loop counter

    logical :: DODERIVATIVES            ! Flag
    logical :: DOELEMENT                ! Flag
    logical :: FOUNDINFIRST             ! Flag for state quantities
    logical :: PTANINFIRST              ! PTan was found in the first vector

    integer, dimension(-1:1) :: L2PCBINS ! Which l2pc to use
    integer, dimension(noMIFsJ) :: mifPointingsLower ! Result of a hunt
    integer, dimension(size(mifPointingsLower)) :: mifPointingsUpper ! mifPointingsLower+1

    logical, dimension(noChans) :: doChannel ! Do this channel?

    real (r8) :: DELTAPHI               ! Difference in geod Angle in l2pc

    ! The `prime' quantities are important.
    ! - yPrime (yP) contains one maf of the relevant radiances, but on the
    !   l2pc's vertical coordinates
    ! - xPrime (xP) is like xStar but contains state values
    ! - kPrime is the jacobian for these two.

    real (r8), dimension(size(mifPointingsLower)) :: lowerWeight ! For interpolation
    real (r8), dimension(size(mifPointingsLower)) :: upperWeight ! For interpolation
    real (r8), dimension(merge(noMifs,0,radiance%template%quantityType==l_opticalDepth)) &
      & :: tangentTemperature ! For optical depth
    real (r8), dimension(noChans) :: thisFraction ! Sideband fraction values

    real (r8), dimension(noMIFs,noChans) :: resultMapped ! Remapped values of result
    real (r8), dimension(noMIFs,noChans) :: dyByDX ! Raw dRad/dPtan
    real (rm), dimension(:,:), pointer :: dense  ! Densified matrix from l2pc
    real (rm), dimension(:,:), pointer :: kBit ! Remapped values of l2pc

    type(coefficients(rv)) :: Coeffs    ! For interpolation
    type(vector_T) :: XP                ! Same form as xStar, contents as x
    type(vector_T) :: YP                ! Same form as yStar,=kstar*(xp-xStar)
    type(vector_T) :: DELTAX            ! xp-xStar
    type(Signal_T), pointer :: signal   ! Signal from the configuration

    type(L2PC_T), pointer :: L2PC     ! The l2pc to use
    type(MatrixElement_T), pointer :: L2PCBLOCK  ! A block from the l2pc
    type(MatrixElement_T), pointer :: JBLOCK     ! A block from the jacobian

    type(VectorValue_T), pointer :: RADINL2PC ! The relevant radiance part of yStar
    type(VectorValue_T), pointer :: STATEQ ! A state vector quantity
    type(VectorValue_T), pointer :: PTAN ! Tangent pressure quantity
    type(VectorValue_T), pointer :: TEMP ! Temperature quantity
    type(VectorValue_T), pointer :: XSTARPTAN ! Tangent pressure in l2pc
    type(VectorValue_T), pointer :: L2PCQ ! A quantity in the l2pc
    type(VectorValue_T), pointer :: SIDEBANDFRACTION ! From the state vector
    type(VectorValue_T), pointer :: LOWERSIDEBANDFRACTION ! From the state vector
    type(VectorValue_T), pointer :: UPPERSIDEBANDFRACTION ! From the state vector
    type(VectorValue_T), pointer :: THISXSTARQ ! Quantitiy from supplied XStar vector

    ! Executable code

    ! Identify the band/maf we're looking for

    maf = fmStat%maf

    call trace_begin ( me, 'LinearizedForwardModel, MAF=', index=maf, &
      & cond=toggle(emit) )

    nullify ( dense )

    signal => fmConf%signals(1)

    doChannel = .true.
    if ( associated ( signal%channels ) ) doChannel = signal%channels

    call SelectL2PCBins ( fmConf, FwdModelIn, FwdModelExtra, &
      & radiance, signal%sideband, maf, &
      & l2pcBins, sidebandStart, sidebandStop, sidebandStep )

    ! If we're doing a split calculation then get the relevant
    ! sideband information.
    ! Note, this is *NOT* the same decision process as for the full forward model
    ! and deliberately so.
    if ( sidebandStart /= sidebandStop .or. fmConf%forceSidebandFraction ) then
      sidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, &
        & sideband=0, noError=.true. )
      lowerSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, &
        & sideband=-1, noError=.true. )
      upperSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, &
        & sideband=1, noError=.true. )
      if (.not. associated (sidebandFraction) .and. .not. &
        & ( associated ( lowerSidebandFraction) .and. &
        &   associated ( upperSidebandFraction ) ) ) &
        & call MLSMessage( MLSMSG_Error, ModuleName, &
        & "No sideband fraction supplied" )
    end if

    ! --------- Loop over sidebands ------------------------------------------------
    do sideband = sidebandStart, sidebandStop, sidebandStep
      ! Setup a sideband fraction array
      ! Note, this is *NOT* the same decision process as for the full forward model
      ! and deliberately so.
      if ( sidebandStart /= sidebandStop .or. fmConf%forceSidebandFraction ) then !????
        if ( associated ( sidebandFraction ) ) then
          thisFraction = sidebandFraction%values(:,1)
          if ( sideband == 1 ) thisFraction = 1.0 - thisFraction
        else
          if ( sideband == -1 ) then
            thisFraction = lowerSidebandFraction%values(:,1)
          else
            thisFraction = upperSidebandFraction%values(:,1)
          end if
        end if
      else                  ! Otherwise, want just unfolded signal
        thisFraction = 1.0
      end if

      ! Load this bin if necessary
      ! But don't load Hessians if we've elected to ignore them
      call PopulateL2PCBin ( l2pcBins(sideband), fmConf%ignoreHessian )

      ! Work out which l2pc bin we want
      l2pc => l2pcDatabase(l2pcBins(sideband))

      ! Set a dimension
      radInl2pc => GetVectorQuantityByType ( &
        & l2pc%j%row%vec, quantityType=l_radiance, &
        & signal=signal%index, sideband=sideband )
      noPointings = radInL2PC%template%noSurfs
      if ( radInL2PC%template%noChans /= noChans ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Channel dimension in l2pc not same as in measurements" )

      ! Now we loop over the quantities in the l2pc file and construct an xPrime
      ! for them
      call cloneVector ( xP, l2pc%j%col%vec, vectorNameText='_xP' )
      call cloneVector ( deltaX, xP, vectorNameText='_deltaX' ) ! sets values to 0.0

      ! Set up some other stuff before main loop
      ! Get the two ptans; we'll need these for interpolation
      ptan => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra, &
        & quantityType = l_ptan, &
        & instrumentModule = radiance%template%instrumentModule,&
        & foundInFirst = ptanInFirst )

      xStarPtan => GetVectorQuantityByType ( l2pc%j%col%vec, &
        & quantityType = l_ptan )

      rowLBlock = 1
      if ( present (jacobian) ) rowJBlock = &
        & FindBlock ( jacobian%row, radiance%index, maf )

      ! Now, if we need any derivatives, we need to setup some arrays
      if ( present(jacobian) ) then
        call Hunt ( xStarPtan%values(:,1), ptan%values(:,maf), mifPointingsLower )
        mifPointingsUpper = mifPointingsLower + 1

        upperWeight = &
          & ( ptan%values(:,maf) - xStarPtan%values(mifPointingsLower,1) ) / &
          & ( xStarPtan%values(mifPointingsUpper,1) - &
          &   xStarPtan%values(mifPointingsLower,1) )
        lowerWeight = 1.0 - upperWeight
        upperWeight = min ( 1.0_r8, max (upperWeight, 0.0_r8 ) )
        lowerWeight = min ( 1.0_r8, max (lowerWeight, 0.0_r8 ) )
      end if

      ! -------- Main loop over xStar quantities -------------------------------
      quantityLoop: do qtyInd = 1, size ( l2pc%j%col%vec%quantities )

        ! Identify this quantity in xStar
        l2pcQ => l2pc%j%col%vec%quantities(qtyInd)

        ! Now see if we are wanting to deal with this
        if ( l2pcQ%template%quantityType == l_ptan ) cycle ! Get this from interpolation
        if ( l2pcQ%template%quantityType == l_vmr ) then
          if (.not. associated(fmConf%molecules) ) cycle
          if ( .not. any (l2pcQ%template%molecule == &
            &   fmConf%molecules)) cycle
          if ( isExtinction(l2pcQ%template%molecule) ) cycle
        end if

        ! Identify this quantity in x
        call FindMatchForL2PCQ ( l2pcQ, fmConf, fwdModelIn, fwdModelExtra, &
          & stateQ, foundInFirst )

        ! If it's not in the state vector, and the l2pc does contain
        ! derivative information for this then make a fuss.
        if ( .not. associated(stateQ) ) then
          do xStarInstance = 1, l2pcQ%template%noInstances
            colLBlock = FindBlock ( l2pc%j%col, l2pcQ%index, xStarInstance )
            if ( l2pc%j%block(rowLBlock,colLBlock)%kind /= M_Absent ) then
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & "No quantity in state vectors found to match %s", &
                & datum=l2pcQ%template%name )
            end if
          end do
          cycle quantityLoop            ! Go to next l2pc quantity
        end if

        if ( toggle(emit) .and. levels(emit) > 1 ) then
          call output ( 'Dealing with xStar Quantity named ' )
          call display_string ( l2pc%j%col%vec%quantities(qtyInd)%template%name, &
            & advance='yes' )
        end if

        ! Do we need derivatives for this?
        doDerivatives = present(jacobian) .and. foundInFirst
        if ( doDerivatives ) then
          select case ( l2pcQ%template%quantityType )
          case ( l_temperature )
            doDerivatives = fmConf%temp_der
          case ( l_vmr )
            doDerivatives = associated(fmConf%molecules) .and. fmConf%atmos_der
            if ( doDerivatives ) then
              doDerivatives = any( l2pcQ%template%molecule == &
                & pack(fmConf%molecules, fmConf%moleculeDerivatives) )
            end if
          case default ! Can we get here?
            doDerivatives = .false.
          end select
        end if

        ! OK, we're legit, let's get going.
        instanceLen = l2pcQ%template%instanceLen
        closestInstance = FindOneClosestInstance ( stateQ, radiance, maf )

        ! Loop over profiles.
        ! Assumes l2pcQ is stacked.
        center = l2pcQ%template%noInstances/2 + 1
        do xStarInstance = 1, l2pcQ%template%noInstances
          ! Identify this instance in state
          if ( sum(fmConf%phiWindow) == 0.0 ) then
            ! If 0, user wants 1D, always choose the closest instance
            deltaInstance = 0
          else
            ! Try to fit within phiWindow
            deltaInstance = xStarInstance - center
            if ( fmConf%windowUnits == PHYQ_Angle ) then
              phiWindowLoop: do
                deltaPhi = l2pcQ%template%phi(1,center+deltaInstance) - &
                  & l2pcQ%template%phi(1,center)
                if ( deltaPhi > fmConf%phiWindow(2) ) then
                  deltaInstance = deltaInstance - 1
                else if ( deltaPhi < -fmConf%phiWindow(1) ) then
                  deltaInstance = deltaInstance + 1
                else
                  exit phiWindowLoop
                end if
              end do phiWindowLoop
            else ! fmConf%windowUnits == PHYQ_Profiles
              deltaInstance = max ( -nint(fmConf%phiWindow(1)), &
                                  & min ( nint(fmConf%phiWindow(2)), &
                                        & deltaInstance ) )
            end if
          end if

          xInstance = closestInstance + deltaInstance
          xInstance = max ( 1, min ( xInstance, stateQ%template%noInstances ) )

          ! Fill this part of xP, if we have been supplied an xStar fill xP with
          ! the delta otherwise, xP is just the direct state
          if ( fmConf%xStar /= 0 ) then
            i = GetVectorQuantityIndexByName ( vectors(fmConf%xStar), stateQ%template%name )
            if ( i == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to find quantity in xStar vector supplied' )
            thisXStarQ => vectors(fmConf%xStar)%quantities(i)
            xP%quantities(qtyInd)%values(:,xStarInstance) = &
              & stateQ%values(:,xInstance) - thisXStarQ%values(:,xInstance)
          else
            xP%quantities(qtyInd)%values(:,xStarInstance) = &
              & stateQ%values(:,xInstance)
          end if

          ! Note that we've ignored the mask in stateQ here. We might
          ! in later versions want to apply one mask field?

          ! Interpolate kStar to K if required.
          if ( present ( jacobian ) ) fmStat%rows(rowJBlock) = .true.
          ! We want to set the row flag even if we don't compute derivatives
          ! as we want our contribution to the radiances known.
          if ( doDerivatives ) then
            colLBlock = FindBlock ( l2pc%j%col, l2pcQ%index, xStarInstance )
            colJBlock = FindBlock ( jacobian%col, stateQ%index, xInstance )
            l2pcBlock => l2pc%j%block(rowLBlock,colLBlock)

            if ( l2pcBlock%kind /= M_Absent ) then
              ! OK, before, I was doing transpose / reshapes / 
              ! interpolates / transpose / reshapes.  That was taking too long,
              ! so I'm going to do the interpolations by hand

              ! Get pointer to l2pc information
              if ( any ( l2pcBlock%kind == &
                & (/ M_Column_Sparse, M_Banded /) ) ) then
                call allocate_test ( dense, &
                  & noChans * noPointings, instanceLen, &
                  & 'dense', ModuleName )
                call densify ( dense, l2pcBlock )
                kBit => dense
              else
                kBit => l2pcBlock%values
              end if

              ! Indentify the block in the jacobian
              jBlock => jacobian%block(rowJblock,colJblock)
              select case ( jBlock%kind )
              case ( M_Absent )
                call CreateBlock ( jBlock, &
                  & noChans*noMIFs, instanceLen, &
                  & M_Full, init=0.0_rm )
              case ( M_Banded, M_Column_Sparse )
                call MLSMessage( MLSMSG_Error, ModuleName, &
                  & "Code not written for adding to non full blocks" )
              case default
              end select

              ! Now do the interpolation
!$OMP PARALLEL DO private ( i, lower, upper, chan, doElement )
              do mif = 1, noMIFs
                i = ( mif - 1 ) * noChans + 1
                lower = (mifPointingsLower(mif)-1)*noChans + 1
                upper = (mifPointingsUpper(mif)-1)*noChans + 1
                do chan = 1, noChans
                  doElement = doChannel(chan)
                  if ( doElement .and. associated ( radiance%mask ) ) &
                    & doElement = iand ( ichar ( radiance%mask(i,maf)), m_linAlg ) == 0
                  if ( doElement ) then
                    ! It's tempting to think we need a if sideband==sidebandStart 
                    ! case here where we ensure we zero out jBlock the first time.
                    ! However! The 2D computation makes that a mistake as at the chunk
                    ! edges, all the derivatives need to be piled up on each other.
                    jBlock%values ( i , : ) = &
                      & jBlock%values ( i , : ) + &
                      &   thisFraction(chan) * ( &
                      &     lowerWeight(mif) * kBit( lower, : ) + &
                      &     upperWeight(mif) * kBit( upper, : ) )
                  end if
                  i = i + 1
                  lower = lower + 1
                  upper = upper + 1
                end do
              end do
!$OMP END PARALLEL DO

              call deallocate_test ( dense, 'dense', ModuleName )
            end if                        ! Any matrix here?
          end if                          ! do derivatives?
        end do                            ! End loop over xStar Profiles

        ! Compute this instance of deltaX
        if ( fmConf%xStar /= 0 ) then
          ! The delta has already been done in the xStar supplied case
          deltaX%quantities(qtyInd)%values = xP%quantities(qtyInd)%values
        else
          ! Otherwise, difference it from the l2pc state
          deltaX%quantities(qtyInd)%values = &
            & xP%quantities(qtyInd)%values - &
            & l2pc%j%col%vec%quantities(qtyInd)%values
        end if

      end do quantityLoop               ! End loop over quantities

      ! Now compute yP
      if ( (deltaX /= 0.0_rv) .and. EXPECTDELTAXALLZERO ) then
        call dump ( fmConf )
        call dump ( deltaX, name='deltaX' )
        call MLSMessage( MLSMSG_Error, ModuleName, &
          & "deltaX not all zero" )
      end if

      if ( toggle(emit) .and. levels(emit) > 8 ) then
        call dump ( deltaX, name='deltaX' )
        call dump ( l2pc%j%col%inst, 'l2pc%j%col%inst' )
        call dump ( l2pc%j%col%quant, 'l2pc%j%col%quant' )

        call dump ( l2pc%j%row%inst, 'l2pc%j%row%inst' )
        call dump ( l2pc%j%row%quant, 'l2pc%j%row%quant' )
      end if

      call cloneVector( yp, l2pc%j%row%vec, vectorNameText='_yP' )
      call MultiplyMatrixVectorNoT ( l2pc%j, deltaX, yP, update = .false. )
      if ( .not. fmConf%ignoreHessian .and. l2pc%goth ) &
        ! We have a Hessian and want to use it; add a term to the Taylor series
        & call Multiply ( l2pc%h, deltaX, yP, 0.5_rv, update = .true. )

      if ( toggle(emit) .and. levels(emit) > 8 ) then
        call dump ( yp, name='yP' )
        call dump ( l2pc%j%row%vec, name='l2pc%j%row%vec' )
      end if

      ! Now, if no yStar has been supplied add the one in the l2pc file to yP
      if ( fmConf%yStar == 0 ) yP = yP + l2pc%j%row%vec

      ! Now we interpolate yP to ptan
      ! Do this using setup, loop, teardown instead of using the array
      ! interpolator to avoid allocating and deallocating ypMapped and
      ! ypMapped = transpose ( &
      !   & reshape ( yp%quantities(1)%values(:,1), (/noChans, noPointings/) ) )
      call InterpolateArraySetup ( &
        & xStarPtan%values(:,1), &      ! OldX
        & ptan%values(:,maf), &         ! NewX
        & 'Spline', &                   ! use spline
        & coeffs, &                     ! the coefficients
        & extrapolate='Constant', &     ! no extrapolation
        & dyByDx=.true. )               ! need more coefficients for this

      ! call outputNamedValue( 'shape(xp)', shape(xStarPtan%values) )
      ! call outputNamedValue( 'shape(yp)', shape(yp%quantities(1)%values) )
      ! call outputNamedValue( 'noChans', noChans )
      ! call outputNamedValue( 'noPointings', noPointings )
      do chan = 1, noChans
        top = noChans*noPointings + chan  - noChans
        ! call outputNamedValue( 'chan', chan )
        ! call outputNamedValue( 'top', top )
        ! call outputNamedValue( 'shape(OldX)', shape(xStarPtan%values(:,1)) )
        ! call outputNamedValue( 'shape(OldY)', shape(yp%quantities(1)%values(chan:top:noChans,1)) )
        ! call outputNamedValue( 'shape(NewX)', shape(ptan%values(:,maf)) )
        ! call outputNamedValue( 'shape(NewY)', shape(resultMapped(:,chan)) )
        call InterpolateValues ( coeffs, &
          & xStarPtan%values(:,1), &                       ! OldX
          & yp%quantities(1)%values(chan:top:noChans,1), & ! OldY
          & ptan%values(:,maf), &                          ! NewX
          & resultMapped(:,chan), &                        ! NewY
          & 'Spline', &                                    ! use spline
          & extrapolate='Constant', &                      ! No extrapolation
          & dyByDx=dyByDx(:,chan) )
      end do

      call InterpolateArrayTeardown ( coeffs )

      ! Either place or add, make decision outside loop!
      if ( sideband == sidebandStart ) then
        do chan = 1, noChans
          if ( doChannel ( chan ) ) then
            do mif = 1, noMIFs
              radiance%values( chan + (mif-1)*noChans, maf ) = &
                & thisFraction(chan)*resultMapped ( mif, chan )
            end do
          end if
        end do
      else
        do chan = 1, noChans
          if ( doChannel ( chan ) ) then
            do mif = 1, noMIFs
              radiance%values( chan + (mif-1)*noChans, maf ) = &
                & radiance%values( chan + (mif-1)*noChans, maf ) + &
                &   thisFraction(chan)*resultMapped ( mif, chan )
            end do
          end if
        end do
      end if

      ! Now deal with ptan derivatives
      if ( present( jacobian ) .and. ptanInFirst ) then
        colJBlock = FindBlock ( jacobian%col, ptan%index, maf )
        jBlock => jacobian%block( rowJBlock, colJBlock )
        select case (jBlock%kind)
        case (m_absent)
          call CreateBlock ( jacobian, rowJBlock, colJBlock, m_banded, &
            & noMIFs*noChans, bandHeight=noChans, init=0.0_rm )
        case (m_banded)
          call CheckForSimpleBandedLayout ( jBlock, noChans, 'jBlock in Linearized model' )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Wrong matrix type for ptan derivative')
        end select

        ! Either place or add, make decision outside loop!
        if ( sideband == sidebandStart ) then
          do chan = 1, noChans
            if ( doChannel ( chan ) ) then
              do mif = 1, noMIFs
                jBlock%values( chan + (mif-1)*noChans, 1 ) = &
                  & thisFraction(chan) * dyByDx ( mif, chan )
              end do                    ! Minor frames
            end if                      ! Doing this channel
          end do                        ! Channels
        else                            ! Must be doing second sideband
          do chan = 1, noChans
            if ( doChannel ( chan ) ) then
              do mif = 1, noMIFs
                jBlock%values( chan + (mif-1)*noChans, : ) = &
                  & jBlock%values ( chan + (mif-1)*noChans, : ) + &
                  &   thisFraction(chan) * dyByDx ( mif, chan )
              end do                    ! Minor frames
            end if                      ! Doing this channel
          end do                        ! Channel
        end if                          ! First/second sideband
      end if                            ! Folding

      ! Now think about optical depth
      if ( radiance%template%quantityType == l_opticalDepth ) then
        ! Get temperature and interpolate to tangent points
        temp => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_temperature, config=fmConf )
        closestInstance = FindOneClosestInstance ( temp, radiance, maf )
        call InterpolateValues ( &
          & temp%template%surfs(:,1), &
          & temp%values(:,closestInstance), &
          & ptan%values(:,maf), &
          & tangentTemperature, &
          & 'Spline', extrapolate='Constant' )
        ! Now convert radiance to optical depth
        radiance%values(:,maf) = 1.0 - &
          & radiance%values(:,maf) / reshape ( spread ( &
          & tangentTemperature, 1, noChans ), (/ radiance%template%instanceLen /) )
        where ( radiance%values(:,maf) > 0.0 )
          radiance%values(:,maf) = -log( radiance%values(:,maf) )
        elsewhere
          radiance%values(:,maf) = 1.0e5 ! Choose a suitable large value
        end where
      end if

      call DestroyVectorInfo ( xP )
      call DestroyVectorInfo ( deltaX )
      call DestroyVectorInfo ( yP )

    end do                                ! End of sideband loop

    ! If a yStar has been supplied, we need to add that on.
    ! We also need to add on the impact of any perturbation in tangent pressure
    ! in the xStar/yStar supplied case too
    if ( fmConf%yStar /= 0 ) call Add_y_Star

    if ( switchDetail(switches,'rad') + switchDetail(switches,'RAD') > -2 ) then
      call dump ( radiance, name='Linearized radiances', details=1, options="c" )
      if ( switchDetail(switches,'RAD') > -1 ) stop
    end if

    call trace_end ( 'LinearizedForwardModel', cond=toggle(emit) )

  contains

    ! ======================================== Internal procedures =====

    ! -----------------------------------------------  Add_y_Star  -----
    subroutine Add_y_Star
    ! If a yStar has been supplied, we need to add that on.
    ! We also need to add on the impact of any perturbation in tangent pressure
    ! in the xStar/yStar supplied case too

      real (r8), dimension(noMIFs) :: deltaPtan ! Change of ptan between x and supplied xStar
      real (r8), dimension(noMIFs, noChans) :: dummy ! Workspace
      real (r8), dimension(noMIFs, noChans) :: dYStarByDPtan ! Remapped values of l2pc
      integer, dimension(noMIFs) :: ptanOrder ! Ordering for ptan
      real (r8), dimension(noMIFs) :: sortedPtan ! supplied ptan in ascending order
      real (r8), dimension(noMIFs, noChans) :: sortedDYStarByDPtan ! Remapped values of l2pc
      real (r8), dimension(noMIFs, noChans) :: sortedYStarMapped ! Remapped values of l2pc
      real (r8), dimension(noMIFs, noChans) :: yStarMapped ! Remapped values of l2pc

      type(VectorValue_T), pointer :: THISYSTARQ ! Quantitiy from supplied YStar vector

      ! Identify the radiance in the supplied yStar vector
      thisYStarQ => GetVectorQuantityByType ( vectors(fmConf%yStar), &
        & quantityType=l_radiance, signal=signal%index, sideband=signal%sideband )

      ! We're going to interpolate this in ptan.  Extract it and make
      ! sure ptan is in ascending order.
      yStarMapped = transpose ( reshape ( thisYStarQ%values(:,maf), (/ noChans, noMIFs /) ) )
      call SortP ( ptan%values(:,maf), 1, noMIFs, ptanOrder )
      do i = 1, noMIFs
        sortedPtan ( i ) = ptan%values ( ptanOrder(i), maf )
        sortedyStarMapped ( i, : ) = yStarMapped ( ptanOrder(i), : )
      end do

      ! 'Interpolate' this to get the derivatives
      call InterpolateValues ( &
        & sortedPtan, &                 ! OldX
        & sortedyStarMapped, &          ! OldY
        & sortedPtan, &                 ! NewX (the same as oldX)
        & dummy, &                      ! NewY (don't care about it)
        & 'Spline', &                   ! Use spline
        & extrapolate='Constant', &     ! No extrapolation
        & dyByDx=sortedDyStarByDPtan, & ! This is what we're after
        & skipNewY=.true. )

      ! Now 'unsort' the result
      do i = 1, noMIFs
        dYStarByDptan ( ptanOrder(i), : ) = sortedDyStarByDPtan ( i, : )
      end do

      ! Now compute 'deltaPtan'
      xStarPtan => GetVectorQuantityByType ( vectors ( fmConf%xStar ), quantityType = l_ptan, &
        & instrumentModule = radiance%template%instrumentModule )
      deltaPtan = ptan%values(:,maf) - xStarPtan%values(:,maf)

      ! Add these impacts to our radiances
      do chan = 1, noChans
        if ( doChannel ( chan ) ) then
          do mif = 1, noMIFs
            i = chan + ( mif - 1 ) * noChans
            radiance%values( i, maf ) = radiance%values( i, maf ) + &
              & thisYStarQ%values ( i, maf ) + &
              &   dYStarByDPtan ( mif, chan ) * deltaPtan(mif)
          end do
        end if
      end do

      ! Output the ptan derivatives if they're wanted
      if ( present ( jacobian ) .and. ptanInFirst ) then
        ! We know the block is already of the right type (banded) from our
        ! work above.
        colJBlock = FindBlock ( jacobian%col, ptan%index, maf )
        jBlock => jacobian%block( rowJBlock, colJBlock )

        do chan = 1, noChans
          if ( doChannel ( chan ) ) then
            do mif = 1, noMIFs
              i = chan + (mif-1) * noChans
              jBlock%values( i, 1 ) = jBlock%values( i, 1 ) + &
                & dyStarByDPtan ( mif, chan )
            end do                    ! Minor frames
          end if                      ! Doing this channel
        end do                        ! Channels
      end if

    end subroutine Add_y_Star

  end subroutine LinearizedForwardModelAuto

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LinearizedForwardModel_m

! $Log$
! Revision 2.92  2016/04/13 00:47:57  vsnyder
! Use operator(/=) from VectorsModule instead of AreEqual function
!
! Revision 2.91  2015/08/25 17:24:32  vsnyder
! PhiWindow is a tuple, with the first element specifying the angles or
! number of profiles/MAFs before the tangent point, and the second
! specifying the angles or number after.  Pay attention to the units of
! PhiWindow instead of assuming it's the number of profiles.
!
! Revision 2.90  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.89  2013/08/12 23:48:09  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.88  2012/07/31 00:45:02  vsnyder
! Comment out OUTPUTNAMEDVALUE in USE since its use is commented out
!
! Revision 2.87  2011/11/11 00:42:06  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.86  2011/08/20 02:31:47  vsnyder
! Delete unused variable declaration
!
! Revision 2.85  2011/08/20 02:10:06  vsnyder
! Simplify emitting an error message
!
! Revision 2.84  2011/05/09 17:50:34  pwagner
! Converted to using switchDetail
!
! Revision 2.83  2011/03/28 17:52:08  pwagner
! Fixed bug in calculating top before call to InterpolateValues
!
! Revision 2.82  2010/09/25 01:08:39  vsnyder
! Cannonball polishing
!
! Revision 2.81  2010/06/29 19:57:58  vsnyder
! Add 0.5 factor to Hessian-Vector-Vector multiply
!
! Revision 2.80  2010/06/07 23:25:18  vsnyder
! Revise how failure to find desired blocks is handled
!
! Revision 2.79  2010/05/19 17:52:53  pwagner
! Removed unused stuff
!
! Revision 2.78  2010/05/19 00:33:46  vsnyder
! Get rid of a temp, hopefully interpolate faster, pass ignoreHessian into
! PopulateL2PCBin, cosmetic changes.
!
! Revision 2.77  2010/05/13 23:46:00  pwagner
! Temporary expedients for l2pc files with Hessians; needs more code
!
! Revision 2.76  2010/04/30 23:57:40  vsnyder
! Remove StateQ from SelectL2PCBins call
!
! Revision 2.75  2010/04/13 01:42:30  vsnyder
! Move FindMatchForL2PCQ, FlushLockedBins, SelectL2PCBins to L2PCBins_m
!
! Revision 2.74  2010/03/26 23:13:12  vsnyder
! Add ignoreHessian field to forward model config
!
! Revision 2.73  2010/02/25 18:02:17  pwagner
! Conforms with changed l2pc structure
!
! Revision 2.72  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.71  2009/04/20 18:44:27  pwagner
! Needed changes when identical types with different names allowed in L2PC files
!
! Revision 2.70  2008/10/03 16:39:29  livesey
! Bug fix with EXTINCTIONV2
!
! Revision 2.69  2008/10/03 16:31:13  livesey
! Added EXTINCTIONV2
!
! Revision 2.68  2008/06/06 22:51:44  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.67  2008/05/07 20:55:32  vsnyder
! OOPS, can't test optional args in specification exprs
!
! Revision 2.66  2008/05/02 20:14:09  vsnyder
! Further simplification using MERGE in dimensions
!
! Revision 2.65  2008/05/02 19:44:39  vsnyder
! Cure three memory leaks, plus some reorg
!
! Revision 2.64  2007/06/29 19:32:42  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.63  2007/04/03 17:45:10  vsnyder
! Add target attribute to FwdModelIn and FwdModelExtra, replace pointer
! attribute on Vectors with Target attribute.
!
! Revision 2.62  2006/07/19 22:33:02  vsnyder
! Cannonball polishing
!
! Revision 2.61  2006/06/19 15:53:35  livesey
! My first 'bug fix' was a mistake
!
! Revision 2.60  2006/06/16 23:41:57  livesey
! Bug fix, if called twice on same band species derivatives doubled rather
! than rewritten
!
! Revision 2.59  2006/06/15 17:40:33  pwagner
! Stops with err instead of continuing with warning if xStar qty not in statevectors
!
! Revision 2.58  2005/06/30 22:42:34  livesey
! Bug fix for split sideband case
!
! Revision 2.57  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.56  2004/11/01 20:24:48  vsnyder
! Reorganization of representation for molecules and beta groups
!
! Revision 2.55  2004/08/06 01:09:49  livesey
! Updated to fix Van's new definition of fmConf%molecules
!
! Revision 2.54  2004/07/07 19:42:11  vsnyder
! Use new Init argument of CreateBlock
!
! Revision 2.53  2004/02/07 00:45:36  livesey
! Fixed a bug in the sideband handling.  Was overzelous yesterday.
!
! Revision 2.52  2004/02/05 23:30:39  livesey
! Fixed long standing problem with single sideband radiometers.
!
! Revision 2.51  2003/11/01 18:44:39  livesey
! Bug fixes in bin selection
!
! Revision 2.50  2003/10/29 00:44:53  livesey
! Bug fix in forceSidebandFraction handling
!
! Revision 2.49  2003/10/28 23:44:15  livesey
! Various changes involved in adding the forceFoldedOutput option to
! support the hybrid model.
!
! Revision 2.48  2003/10/09 22:17:10  livesey
! Various bug fixes, and added the call to CheckForSimpleBandedLayout
!
! Revision 2.47  2003/09/15 17:42:52  livesey
! Various bug fixes associated with getting the polar linear model
! working.
!
! Revision 2.46  2003/09/11 23:10:32  livesey
! Added option to linearize around pre-computed state/radiances instead of
! those in the l2pc file.
!
! Revision 2.45  2003/08/20 20:07:22  livesey
! Bug fix in ptan derivatives (not actually a problem yet but would have
! become one).  Also, more commented out if statements for the single
! sideband case, and some cosmetic changes.
!
! Revision 2.44  2003/08/14 20:24:08  livesey
! Added the exact bin criterion
!
! Revision 2.43  2003/08/13 00:48:23  livesey
! Removed the faking of binSelectors, forwardModelSupport now ensures
! there's always one to hand.
!
! Revision 2.42  2003/07/15 22:10:38  livesey
! Added support for hybrid model
!
! Revision 2.41  2003/07/09 20:16:26  livesey
! Anticipative bug fix in sideband stuff (not really used so far)
!
! Revision 2.40  2003/05/29 16:37:45  livesey
! Renamed sideband fraction
!
! Revision 2.39  2003/02/22 00:40:22  livesey
! Now can use OpenMP
!
! Revision 2.38  2003/02/17 00:33:41  livesey
! Added deallocation of possible/cost to fix memory leak.
!
! Revision 2.37  2003/02/06 01:13:28  livesey
! Added the dump stuff
!
! Revision 2.36  2003/02/06 00:46:03  livesey
! New SelectL2PCBins routine.
!
! Revision 2.35  2003/02/05 21:58:31  livesey
! Just a tidy up in preparation for the new bin selectors stuff
!
! Revision 2.34  2003/01/29 22:43:25  livesey
! More intelligent setting of rows flags.
!
! Revision 2.33  2003/01/08 23:50:15  livesey
! Now doesn't bother to create derivatives for masked radiances
!
! Revision 2.32  2002/12/16 15:31:57  mjf
! Use GetQuantityForForwardModel instead of GetVectorQuantityByType.
!
! Revision 2.31  2002/11/14 17:54:06  livesey
! Changed 'no quantity found' warning message to work better.
!
! Revision 2.30  2002/11/14 17:46:08  livesey
! Bypassed a warning message that was starting to irritate me.
!
! Revision 2.29  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.28  2002/10/02 23:19:51  livesey
! Various bug fixes associated with extinction
!
! Revision 2.27  2002/10/02 02:42:55  livesey
! Changes to allow the new mixed extinction stuff
!
! Revision 2.26  2002/09/13 22:53:22  vsnyder
! Move USE statements from module scope to procedure scope.  Cosmetic changes.
! Move some loop-invariant allocate/deallocates out of loops.
!
! Revision 2.25  2002/09/13 18:09:09  pwagner
! May change matrix precision rm from r8
!
! Revision 2.24  2002/09/11 17:43:39  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.23  2002/09/03 23:46:04  livesey
! Moved debug print statement to location where it will be less verbose.
!
! Revision 2.22  2002/07/22 03:25:50  livesey
! Minor bug fix
!
! Revision 2.21  2002/07/17 06:02:01  livesey
! Got HDF5 l2pcs working
!
! Revision 2.20  2002/07/09 17:37:45  livesey
! Removed ptan extrapolation
!
! Revision 2.19  2002/07/02 19:51:21  livesey
! Embarassing bug, was choosing 'center' wrong, lead to lag
! and all sorts of messyness
!
! Revision 2.18  2002/06/12 17:46:21  livesey
! Better treatment of real phiWindow
!
! Revision 2.17  2002/06/12 17:02:25  livesey
! Very TEMPORARY! change to LinearizedForwardModel to skip round real
! value of phiWindow.  Fix it properly later.
!
! Revision 2.16  2002/05/28 22:34:21  livesey
! Added more useful diagnostic message
!
! Revision 2.15  2002/05/14 22:32:53  livesey
! Added single sideband stuff
!
! Revision 2.14  2002/05/03 23:29:31  livesey
! Added split sideband ratio stuff
!
! Revision 2.13  2002/03/15 21:23:09  livesey
! New bin selectors based on 'latitude' allowed
!
! Revision 2.12  2002/02/20 03:00:23  livesey
! Made the verbose dump more unlikely
!
! Revision 2.11  2002/02/12 21:48:21  livesey
! Fixed minor bug for case when no bin selectors
!
! Revision 2.10  2002/02/09 21:35:22  livesey
! Minor bug fixes
!
! Revision 2.9  2002/02/09 19:11:19  livesey
! Added optical depth calculation etc.
!
! Revision 2.8  2002/02/08 22:51:00  livesey
! Working (hopefully) version of bin selection
!
! Revision 2.7  2002/02/06 01:34:20  livesey
! Changed to use FindOneClosestInstance
!
! Revision 2.6  2002/01/21 22:52:21  livesey
! A few more comments
!
! Revision 2.5  2002/01/21 22:46:08  livesey
! Wrote, but not hooked up SelectL2PCBin
!
! Revision 2.4  2002/01/17 02:17:04  livesey
! Somewhat temporary, skip l2pc state vector components we don't like
!
! Revision 2.3  2001/10/22 16:51:40  livesey
! Allow radiance extrapolation (at least for now).
!
! Revision 2.2  2001/10/07 00:23:50  livesey
! Fixed minor bug with dyByDx
!
! Revision 2.1  2001/10/02 16:51:41  livesey
! Removed fmStat%finished and reordered loops in forward models
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.22  2001/06/04 22:43:26  livesey
! Now works when no molecules at all
!
! Revision 1.21  2001/06/01 20:35:55  livesey
! Now obeys the phiWindow parameter to give 1D forward models.
!
! Revision 1.20  2001/06/01 20:09:28  livesey
! Embarassing bug to do with edge effects removed.  This version makes L2 work!
!
! Revision 1.19  2001/05/25 19:12:21  livesey
! Embarassing bug fix.  Was using quantity basis as framework for interpolation,
! not pointings in the l2pc file!
!
! Revision 1.18  2001/05/19 01:21:33  livesey
! Bug fix, was only storing derivatives for second sideband.
!
! Revision 1.17  2001/05/09 19:46:49  vsnyder
! Use new bandHeight argument of createBlock
!
! Revision 1.16  2001/05/08 23:26:36  livesey
! Embarassing typo/bug fixed
!
! Revision 1.15  2001/05/03 22:58:38  livesey
! Removed bug work around
!
! Revision 1.14  2001/05/03 02:02:54  vsnyder
! Use emit toggle.  Add names to cloned vectors
!
! Revision 1.13  2001/05/02 20:26:17  vsnyder
! Give names to cloned vectors.  Use toggle(emit) and levels(emit) to
! control debugging output.
!
! Revision 1.12  2001/05/02 20:22:13  livesey
! Removed some unused variables.
!
! Revision 1.11  2001/05/01 06:55:07  livesey
! Intermediate bug avoiding version.
!
! Revision 1.10  2001/04/30 15:50:29  livesey
! Commented out code associated with mask
!
! Revision 1.9  2001/04/28 22:32:07  livesey
! Nothing much changed apart from comments.
! Need to add code to have it ignore state stuff according to MASK.
!
! Revision 1.8  2001/04/28 21:23:33  livesey
! Working version.  Does folded sidebands (not tested yet), nice and fast
!
! Revision 1.7  2001/04/28 19:40:37  livesey
! Working version. Assembles jacobian in a nice fast manner.
! Next thing is to do sideband folding.
!
! Revision 1.6  2001/04/28 07:05:25  livesey
! Another interim version
!
! Revision 1.5  2001/04/28 04:40:56  livesey
! Another interim version.
!
! Revision 1.4  2001/04/28 01:54:54  livesey
! Interim, non working but compiling version
!
! Revision 1.3  2001/04/27 17:35:42  livesey
! An interim version for some testing.
!
! Revision 1.2  2001/04/26 23:54:32  livesey
! Interim version
!
! Revision 1.1  2001/04/25 20:49:27  vsnyder
! Initial commit
!
