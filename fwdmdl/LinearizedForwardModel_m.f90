! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module LinearizedForwardModel_m

  ! Approximate a forward model with a first-order Taylor's series.

  implicit none
  private
  public :: FlushLockedBins, LinearizedForwardModel

  ! This array is used to keep track of which bins to use for each (side)band.
  integer, dimension(:,:), pointer, private, save :: lockedBins => NULL()

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------  FlushLockedBins  -----
  subroutine FlushLockedBins
    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use L2PC_m, only: FLUSHL2PCBINS

    ! This could set them to zero again, but I think I'll deallocate here, just
    ! because otherwise I'd have to write a separate deallocate routine.
    call deallocate_test ( lockedBins, 'lockedBins', moduleName )
    call FlushL2PCBins
  end subroutine FlushLockedBins

  ! -------------------------------------  LinearizedForwardModel  -----
  subroutine LinearizedForwardModel ( fmConf, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, Ifm, fmStat, Jacobian, vectors )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Dump_0, only: DUMP
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T, &
      & FORWARDMODELINTERMEDIATE_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: LIT_INDICES
    use Intrinsic, only: L_NONE, L_RADIANCE, L_TEMPERATURE, L_PTAN, L_VMR, &
      & L_LIMBSIDEBANDFRACTION, L_ZETA, L_OPTICALDEPTH, L_LATITUDE, L_FIELDSTRENGTH, &
      & L_FIELDELEVATION, L_FIELDAZIMUTH
    use L2PC_m, only: L2PCDATABASE, BINSELECTORS, POPULATEL2PCBIN
    use ManipulateVectorQuantities, only: FINDONECLOSESTINSTANCE, &
      & DOHGRIDSMATCH, DOVGRIDSMATCH, DOFGRIDSMATCH
    use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
      & MATRIXELEMENT_T, CREATEBLOCK, DENSIFY
    use MatrixModule_1, only: MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP, &
      & FINDBLOCK, CREATEBLOCK
    use MLSCommon, only: r8, rm
    use MLSSignals_m, only: Signal_T, GetSidebandLoop, GetSignalName
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
      & MLSMSG_Allocate, MLSMSG_Deallocate, MLSMSG_WARNING
    use MLSNumerics, only: HUNT, INTERPOLATEVALUES
    use Molecules, only: L_EXTINCTION
    use Output_m, only: Output
    use QuantityTemplates, only: QuantityTemplate_T
    use String_Table, only: Display_String, Get_String
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_begin, Trace_end
    use VectorsModule, only: assignment(=), OPERATOR(-), OPERATOR(+), &
      & CLONEVECTOR, CONSTRUCTVECTORTEMPLATE, COPYVECTOR, CREATEVECTOR,&
      & DESTROYVECTORINFO, GETVECTORQUANTITYBYTYPE, VECTOR_T, &
      & VECTORVALUE_T, VECTORTEMPLATE_T, DUMP, &
      & VALIDATEVECTORQUANTITY, M_LINALG, GETVECTORQUANTITYINDEXBYNAME
    use Sort_m, only: SORTP

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FMCONF
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_t), dimension(:), pointer, optional :: VECTORS ! Vectors database

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
    integer :: MIF                      ! Minor frame loop counter
    integer :: NOCHANS                  ! Dimension
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOPOINTINGS              ! Number of pointings in the l2pc file
    integer :: QTYIND                   ! Loop index for main loop
    integer :: ROWJBLOCK                ! Row index in jacobian
    integer :: ROWLBLOCK                ! Row index in l2pc
    integer :: SIDEBAND                 ! Loop index
    integer :: SIDEBANDSTART            ! For sideband loop
    integer :: SIDEBANDSTEP             ! For sideband loop
    integer :: SIDEBANDSTOP             ! For sideband loop
    integer :: UPPER                    ! Array index
    integer :: XINSTANCE                ! Instance in x corresponding to xStarInstance
    integer :: XSTARINSTANCE            ! Loop counter

    logical :: ANYBLOCKS                ! Flag for checking derivatives
    logical :: DODERIVATIVES            ! Flag
    logical :: DOELEMENT                ! Flag
    logical :: FOUNDINFIRST             ! Flag for state quantities
    logical :: PTANINFIRST              ! PTan was found in the first vector

    integer, dimension(-1:1) :: L2PCBINS ! Which l2pc to use
    integer, dimension(:), pointer :: mifPointingsLower ! Result of a hunt
    integer, dimension(:), pointer :: mifPointingsUpper ! mifPointingsLower+1
    integer, dimension(:), pointer :: ptanOrder ! Ordering for ptan

    logical, dimension(:), pointer :: doChannel ! Do this channel?

    real (r8) :: BESTCOST               ! Output from SelectL2PCBin
    real (r8) :: DELTAPHI               ! Difference in geod Angle in l2pc

    character (len=80) :: WORD          ! A word to output

    ! The `prime' quantities are important.
    ! - yPrime (yP) is contains one maf of the relevant radiances, but on the
    !   l2pc's vertical coordinates
    ! - xPrime (xP) is like xStar but contains state values
    ! - kPrime is the jacobian for these two.

    real (r8), dimension(:), pointer :: lowerWeight ! For interpolation
    real (r8), dimension(:), pointer :: upperWeight ! For interpolation
    real (r8), dimension(:), pointer :: tangentTemperature ! For optical depth
    real (r8), dimension(:), pointer :: thisFraction ! Sideband fraction values
    real (r8), dimension(:), pointer :: deltaPtan ! Change of ptan between x and supplied xStar
    real (r8), dimension(:), pointer :: sortedPtan ! supplied ptan in ascending order

    real (r8), dimension(:,:), pointer :: yPmapped ! Remapped values of yP
    real (r8), dimension(:,:), pointer :: resultMapped ! Remapped values of result
    real (r8), dimension(:,:), pointer :: dyByDX ! Raw dRad/dPtan
    real (rm), dimension(:,:), pointer :: dense  ! Densified matrix from l2pc
    real (rm), dimension(:,:), pointer :: kBit ! Remapped values of l2pc
    real (r8), dimension(:,:), pointer :: dYStarByDPtan ! Remapped values of l2pc
    real (r8), dimension(:,:), pointer :: sortedDYStarByDPtan ! Remapped values of l2pc
    real (r8), dimension(:,:), pointer :: yStarMapped ! Remapped values of l2pc
    real (r8), dimension(:,:), pointer :: sortedYStarMapped ! Remapped values of l2pc
    real (r8), dimension(:,:), pointer :: dummy ! Workspace

    type(vector_T) :: XP                ! Same form as xStar, contents as x
    type(vector_T) :: YP                ! Same form as yStar,=kstar*(xp-xStar)
    type(vector_T) :: DELTAX            ! xp-xStar
    type(Signal_T) :: signal            ! A signal

    type(Matrix_T), pointer :: L2PC     ! The l2pc to use
    type(MatrixElement_T), pointer :: L2PCBLOCK  ! A block from the l2pc
    type(MatrixElement_T), pointer :: JBLOCK     ! A block from the jacobian

    type(VectorValue_T), pointer :: RADIANCE ! The radiance quantity to fill
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
    type(VectorValue_T), pointer :: THISYSTARQ ! Quantitiy from supplied YStar vector

    ! Executable code
    if ( toggle(emit) ) call trace_begin ( 'LinearizedForwardModel' )

    nullify ( yPmapped, resultMapped, dyByDX )
    nullify ( dense, mifPointingsLower, mifPointingsUpper )
    nullify ( lowerWeight, upperWeight )
    nullify ( thisFraction, doChannel )
    nullify ( tangentTemperature )

    ! Identify the band/maf we're looking for

    maf = fmStat%maf

    if ( toggle(emit) .and. levels(emit) > 0 ) then
      call output ( 'Linear model doing maf: ' )
      call output ( maf, advance='yes' )
    end if
    if ( size ( fmConf%signals ) /= 1 ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Can only have one signal for linearized models')

    signal = fmConf%signals(1)
    radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal=signal%index, sideband=signal%sideband, noError=.true. )

    ! Now, it's possible we're really being asked to deal with optical depth, not
    ! radiance.
    if ( .not. associated ( radiance ) ) then
      radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_opticalDepth, &
        & signal=signal%index, sideband=signal%sideband, noError=.true. )
    end if

    ! Now, some possible error messages
    if ( .not. associated ( radiance ) ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Unable to find a radiance or optical depth quantity for this signal' )
    if ( radiance%template%quantityType == l_opticalDepth .and. &
      & present(jacobian)  ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Not appropriate to ask for derivatives for optical depth' )
    if ( radiance%template%quantityType == l_opticalDepth .and. &
      & fmConf%xStar /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Not appropriate to supply x/yStar for optical depth' )
    if ( fmConf%xStar /= 0 .and. .not. present(vectors) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'x/yStar supplied by vectors database argument not present' )

    ! Set some dimensions
    noChans = radiance%template%noChans
    noMIFs = radiance%template%noSurfs

    call allocate_test ( thisFraction, noChans, 'thisFraction', ModuleName )
    call allocate_test ( doChannel, noChans, 'doChannel', ModuleName )

    doChannel = .true.
    if ( associated ( signal%channels ) ) doChannel = signal%channels

    call SelectL2PCBins ( radiance, maf, l2pcBins, &
      & sidebandStart, sidebandStop, sidebandStep )

    ! If we're doing a split calculation then get the relevant
    ! sideband information
    if ( sidebandStart /= sidebandStop .or. fmConf%forceSidebandFraction ) then 
      ! Change to this if statement later !????
      ! if ( signal%sideband == 0 .or. fmConf%forceSidebandFraction ) then
      sidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, noError=.true. )
      lowerSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, &
        & sideband=-1, noError=.true. )
      upperSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType = l_limbSidebandFraction, signal=signal%index, &
        & sideband=1, noError=.true. )
      if (.not. associated (sidebandFraction) .and. .not. &
        & ( associated ( lowerSidebandFraction) .and. &
        &   associated ( upperSidebandFraction ) ) ) &
        & call MLSMessage(MLSMSG_Error,ModuleName, &
        & "No sideband fraction supplied")
    end if

    ! Now, if we need any derivatives, we need to setup some arrays
    if ( present(jacobian) ) then
      call allocate_test ( mifPointingsLower, noMIFs, &
        & 'mifPointingsLower', ModuleName )
      call allocate_test ( mifPointingsUpper, noMIFs, &
        & 'mifPointingsUpper', ModuleName )
      call allocate_test ( lowerWeight, noMIFs, &
        & 'lowerWeight', ModuleName )
      call allocate_test ( upperWeight, noMIFs, &
        & 'upperWeight', ModuleName )
    end if

    ! --------- Loop over sidebands ------------------------------------------------
    do sideband = sidebandStart, sidebandStop, sidebandStep
      ! Load this bin if necessary
      call PopulateL2PCBin ( l2pcBins(sideband) )

      ! Setup a sideband fraction array
      if ( sidebandStart /= sidebandStop ) then   ! We're folding
        ! Change to this if statement later
        ! if ( signal%sideband == 0 .and. .not. fmConf%forceSidebandFraction ) then !????
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

      ! Work out which l2pc bin we want
      l2pc => l2pcDatabase(l2pcBins(sideband))

      ! Set a dimension
      radInl2pc => GetVectorQuantityByType ( &
        & l2pc%row%vec, quantityType=l_radiance, &
        & signal=signal%index, sideband=sideband )
      noPointings = radInL2PC%template%noSurfs
      if ( radInL2PC%template%noChans /= noChans ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Channel dimension in l2pc not same as in measurements" )

      ! Now we loop over the quantities in the l2pc file and construct an xPrime
      ! for them
      call cloneVector ( xP, l2pc%col%vec, vectorNameText='_xP' )
      call cloneVector ( deltaX, xP, vectorNameText='_deltaX' ) ! sets values to 0.0

      ! Set up some other stuff before main loop
      ! Get the two ptans; we'll need these for interpolation
      ptan => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra, &
        & quantityType = l_ptan, &
        & instrumentModule = radiance%template%instrumentModule,&
        & foundInFirst = ptanInFirst )

      xStarPtan => GetVectorQuantityByType ( l2pc%col%vec, &
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
      quantityLoop: do qtyInd = 1, size ( l2pc%col%vec%quantities )

        ! Identify this quantity in xStar
        l2pcQ => l2pc%col%vec%quantities(qtyInd)

        ! Now see if we are wanting to deal with this
        if ( l2pcQ%template%quantityType == l_ptan ) cycle ! Get this from interpolation
        if ( l2pcQ%template%quantityType == l_vmr ) then
          if (.not. associated(fmConf%molecules) ) cycle
          if ( .not. any (l2pcQ%template%molecule == &
            &   fmConf%molecules)) cycle
          if ( l2pcQ%template%molecule == l_extinction .and. &
            & l2pcQ%template%radiometer /= signal%radiometer ) cycle
        end if

        ! Identify this quantity in x
        call FindMatchForL2PCQ ( l2pcQ, fmConf, fwdModelIn, fwdModelExtra, &
          & stateQ, foundInFirst )

        ! If it's not in the state vector, and the l2pc does contain
        ! derivative information for this then make a fuss.
        if ( .not. associated(stateQ) ) then
          anyBlocks = .false.
          do xStarInstance = 1, l2pcQ%template%noInstances
            colLBlock = FindBlock ( l2pc%col, l2pcQ%index, xStarInstance )
            if ( l2pc%block(rowLBlock,colLBlock)%kind /= M_Absent ) then
              anyBlocks = .true.
              exit
            end if
          end do
          if ( anyBlocks ) then
            call get_string ( l2pcQ%template%name, word, strip=.true. )
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              &  "No quantity in state vectors found to match "//trim(word) )
          end if
          cycle quantityLoop            ! Go to next l2pc quantity
        end if

        if ( toggle(emit) .and. levels(emit) > 1 ) then
          call output ( 'Dealing with xStar Quantity named ' )
          call display_string ( l2pc%col%vec%quantities(qtyInd)%template%name, &
            & advance='yes' )
        end if

        ! OK, we're legit, lets get going.
        instanceLen = l2pcQ%template%instanceLen
        closestInstance = FindOneClosestInstance ( stateQ, radiance, maf )

        ! Do we need derivatives for this?
        doDerivatives = present(jacobian) .and. foundInFirst
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_temperature) ) &
          & doDerivatives = fmConf%temp_der
        if (doDerivatives .and. &
          & (l2pcQ%template%quantityType==l_vmr) .and. &
          & associated(fmConf%molecules) ) then
          doDerivatives = fmConf%atmos_der
          if ( doDerivatives .and. .not. any (l2pcQ%template%molecule == &
            & pack(fmConf%molecules, fmConf%moleculeDerivatives))) &
            & doDerivatives = .false.
        end if

        ! Loop over profiles
        center = l2pcQ%template%noInstances/2 + 1
        do xStarInstance = 1, l2pcQ%template%noInstances
          ! Identify this instance in state
          if ( fmConf%phiWindow == 0.0 ) then
            ! If 0, user wants 1D, always choose the closest instance
            deltaInstance = 0
          else
            deltaInstance = xStarInstance - center
            ! Try to fit within phiWindow
            phiWindowLoop: do
              deltaPhi = l2pcQ%template%phi(1,center+deltaInstance) - &
                & l2pcQ%template%phi(1,center)
              if ( deltaPhi > fmConf%phiWindow ) then
                deltaInstance = deltaInstance - 1
              else if ( deltaPhi < -fmConf%phiWindow ) then
                deltaInstance = deltaInstance + 1
              else
                exit phiWindowLoop
              end if
            end do phiWindowLoop
          end if

          xInstance = closestInstance + deltaInstance
          xInstance = max ( 1, min ( xInstance, stateQ%template%noInstances ) )

          ! Fill this part of xP, if we have been supplied an xStar fill xP with the delta
          ! otherwise, xP is just the direct state
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
          if (doDerivatives) then
            colLBlock = FindBlock ( l2pc%col, l2pcQ%index, xStarInstance )
            colJBlock = FindBlock ( jacobian%col, stateQ%index, xInstance )
            l2pcBlock => l2pc%block(rowLBlock,colLBlock)

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
                  & M_Full )
                jBlock%values = 0.0_rm
              case ( M_Banded, M_Column_Sparse )
                call MLSMessage( MLSMSG_Error, ModuleName, &
                  & "Not written code for adding to non full blocks" )
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

              if ( any ( l2pcBlock%kind == &
                & (/ M_Column_Sparse, M_Banded /) ) ) &
                & call deallocate_test ( dense, 'dense', ModuleName )
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
            & l2pc%col%vec%quantities(qtyInd)%values
        end if

      end do quantityLoop               ! End loop over quantities

      ! Now compute yP
      if ( toggle(emit) .and. levels(emit) > 8 ) then
        call dump ( (/deltaX/) )

        call dump ( l2pc%col%inst, 'l2pc%col%inst' )
        call dump ( l2pc%col%quant, 'l2pc%col%quant' )

        call dump ( l2pc%row%inst, 'l2pc%row%inst' )
        call dump ( l2pc%row%quant, 'l2pc%row%quant' )
      end if

      call cloneVector( yp, l2pc%row%vec, vectorNameText='_yP' )
      call MultiplyMatrixVectorNoT ( l2pc, deltaX, yP, update = .false. )

      if ( toggle(emit) .and. levels(emit) > 8 ) then
        call dump ( (/yp, l2pc%row%vec/) )
      end if

      ! Now, if no yStar has been supplied add yP to the one in the l2pc file
      if ( fmConf%yStar == 0 ) yP = yP + l2pc%row%vec

      ! Now we interpolate yP to ptan
      call allocate_test( ypMapped, radInL2PC%template%noSurfs, &
        & noChans, 'ypMapped', ModuleName )
      call allocate_test( resultMapped, radiance%template%noSurfs, &
        & noChans, 'resultMapped', ModuleName )
      call allocate_test( dyByDx, radiance%template%noSurfs, &
        & noChans, 'dyByDX', ModuleName )

      yPmapped = transpose ( &
        & reshape ( yp%quantities(1)%values(:,1), &
        &           (/radInL2PC%template%noChans, radInL2PC%template%noSurfs/) ) )

      call InterpolateValues ( &
        & xStarPtan%values(:,1), &      ! OldX
        & yPmapped, &                   ! OldY
        & ptan%values(:,maf), &         ! NewX
        & resultMapped, &               ! NewY
        & 'Spline', &                   ! use spline
        & extrapolate='Constant', &     ! No extrapolation
        & dyByDx=dyByDx )

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
            & noMIFs*noChans, bandHeight=noChans )
          jBlock%values = 0.0_rm
        case (m_banded)
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
        call Allocate_test ( tangentTemperature, noMIFs, &
          & 'tangentTemperature', ModuleName )
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
      call Deallocate_test ( dyByDx, 'dyByDx', ModuleName)
      call Deallocate_test ( resultMapped, 'resultMapped', ModuleName )
      call Deallocate_test ( ypMapped, 'ypMapped', ModuleName )

    end do                                ! End of sideband loop

    ! If a yStar has been supplied, we need to add that on.
    ! We also need to add on the impact of any perturbation in tangent pressure
    ! in the xStar/yStar supplied case too
    if ( fmConf%yStar /= 0 ) then
      ! Setup arrays etc.
      nullify ( deltaPtan, dummy, dyStarByDPtan, ptanOrder, sortedPtan, sorteddYStarByDPtan, &
        & sortedyStarMapped, yStarMapped )
      call Allocate_test ( deltaPtan, noMIFs, 'deltaPtan', ModuleName )
      call Allocate_test ( dummy, noMIFs, noChans, 'dummy', ModuleName )
      call Allocate_test ( dyStarByDptan, noMIFs, noChans, 'dyStarByDPtan', ModuleName )
      call Allocate_test ( ptanOrder, noMIFs, 'ptanOrder', ModuleName )
      call Allocate_test ( sortedPtan, noMIFs, 'sortedPtan', ModuleName )
      call Allocate_test ( sorteddyStarByDptan, noMIFs, noChans, 'sorteddyStarByDPtan', ModuleName )
      call Allocate_test ( sortedyStarMapped, noMIFs, noChans, 'sortedyStarMapped', ModuleName )
      call Allocate_test ( yStarMapped, noMIFs, noChans, 'yStarMapped', ModuleName )

      ! Identify the radiance in the supplied yStar vector
      thisYStarQ => GetVectorQuantityByType ( vectors(fmConf%yStar), quantityType=l_radiance, &
        & signal=signal%index, sideband=signal%sideband )

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

      ! Extra tidying up in the supplied xStar/yStar case
      call Deallocate_test ( deltaPtan, 'deltaPtan', ModuleName )
      call Deallocate_test ( dummy, 'dummy', ModuleName )
      call Deallocate_test ( dyStarByDptan, 'dyStarByDPtan', ModuleName )
      call Deallocate_test ( ptanOrder, 'ptanOrder', ModuleName )
      call Deallocate_test ( sortedPtan, 'sortedPtan', ModuleName )
      call Deallocate_test ( sorteddyStarByDptan, 'sorteddyStarByDPtan', ModuleName )
      call Deallocate_test ( sortedyStarMapped, 'sortedyStarMapped', ModuleName )
      call Deallocate_test ( yStarMapped, 'yStarMapped', ModuleName )
    end if

    if ( present (jacobian) ) then         ! Destroy working arrays
      call deallocate_test ( mifPointingsLower, 'mifPointingsLower', ModuleName )
      call deallocate_test ( mifPointingsUpper, 'mifPointingsUpper', ModuleName )
      call deallocate_test ( lowerWeight, 'lowerWeight', ModuleName )
      call deallocate_test ( upperWeight, 'upperWeight', ModuleName )
    end if

    if ( toggle(emit) ) call trace_end ( 'LinearizedForwardModel' )

  contains
    ! ======================================== Internal procudures =====

    ! ------------------------------------------ FindMatchForL2PCQ ---
    subroutine FindMatchForL2PCQ ( l2pcQ, fmConf, FwdModelIn, FwdModelExtra, &
      & stateQ, foundInFirst )
      type (VectorValue_T), intent(in) :: L2PCQ ! Quantity to search for
      type (ForwardModelConfig_T), intent(in) :: fmConf ! Forward model config
      type (Vector_T), intent(in), target :: FWDMODELIN ! State vector
      type (Vector_T), intent(in), target :: FWDMODELEXTRA ! Extra state vector
      type (VectorValue_T), pointer :: STATEQ ! Result
      logical, intent(out) :: FOUNDINFIRST ! If set, found in first vector
      ! This routine looks through the supplied state vectors and finds a match
      ! for the supplied quantity from the l2pc state vector.  It then goes on to
      ! test that the HGrids and VGrids for the two quantities match appropriately.

      ! Local variables
      integer :: QTY, VEC
      type (Vector_T), pointer :: V

      ! Executable code
      foundInFirst = .false.
      select case ( l2pcQ%template%quantityType )
      case ( l_temperature )
        stateQ => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_temperature, config=fmConf, &
          & foundInFirst = foundInFirst, noError=.true. )
      case ( l_vmr )
        ! Here we may need to be a little more intelligent
        if ( l2pcQ%template%molecule /= l_extinction ) then
          stateQ => GetQuantityForForwardModel ( FwdModelIn, FwdModelExtra,&
            & quantityType = l_vmr, config=fmConf, &
            & molecule = l2pcQ%template%molecule, &
            & foundInFirst = foundInFirst, noError=.true. )
        else
          searchLoop: do vec = 1, 2
            ! Point to appropraite vector
            if ( vec == 1 ) then
              v => FwdModelIn
            else
              v => FwdModelExtra
            end if
            ! Skip this vector if it doesn't exist
            if ( .not. associated ( v ) ) cycle searchLoop
            ! Loop over this vector
            do qty = 1, size ( v%quantities )
              stateQ => v%quantities(qty)
              if ( stateQ%template%quantityType == l_vmr .and. &
                &  stateQ%template%molecule == l_extinction .and. &
                &  stateQ%template%radiometer == l2pcQ%template%radiometer ) then
                if ( DoFGridsMatch ( l2pcQ, stateQ ) ) exit searchLoop
              end if
            end do
          end do searchLoop
          foundInFirst = ( vec == 1 )
        end if
      case default
        ! For the moment, just ignore things we don't understand.
        stateQ => NULL()
      end select

      ! Now check that these match.
      if ( associated ( stateQ ) ) then
        if ( .not. DoVGridsMatch ( stateQ, l2pcQ ) ) stateQ => NULL()
      end if
      if ( associated ( stateQ ) ) then
        if ( .not. DoHGridsMatch ( stateQ, l2pcQ, spacingOnly=.true. ) ) stateQ => NULL()
      end if

      if ( .not. associated ( stateQ ) ) foundInFirst = .false.

    end subroutine FindMatchForL2PCQ

    ! ---------------------------------------------  NormalizePhi  -----
    elemental real (r8) function NormalizePhi ( phi ) result ( lat )
      ! This does an approximate phi to latitude conversion
      real (r8), intent(in) :: PHI      ! Input geodetic angle
      ! Executable code
      lat = modulo ( phi, 360.0_r8 )
      if ( (lat > 90.0) .and. (lat <= 180.0) ) then
        lat = 180.0 - lat
      else if ( ( lat > 180.0 ) .and. ( lat <= 270.0 ) ) then
        lat = 180.0 - lat
      else if ( lat > 270.0 ) then
        lat = lat - 360.0
      end if
    end function NormalizePhi

    ! --------------------------------------------- SelectL2PCBins -----
    subroutine SelectL2PCBins ( radiance, maf, l2pcBins, &
      & sidebandStart, sidebandStop, sidebandStep )
      use MLSSignals_m, only: signals
      use L2PC_M, only: BINSELECTORS, BINSELECTOR_T
      use Intrinsic, only: L_NAMEFRAGMENT, L_SZA
      use Toggles, only: SWITCHES
      use MLSNumerics, only: ESSENTIALLYEQUAL

      type (VectorValue_T), intent(in) :: RADIANCE ! The radiance we're after
      integer, intent(in) :: MAF                      ! MAF index
      integer, dimension(-1:1), intent(out) :: L2PCBINS ! Result
      integer, intent(out) :: SIDEBANDSTART ! Resulting loop indices
      integer, intent(out) :: SIDEBANDSTOP
      integer, intent(out) :: SIDEBANDSTEP

      ! Local variables
      character (len=132) :: BINNAME      ! The name of the bin
      character (len=132) :: NAMEFRAGMENT ! Fragment of name to try to match

      integer :: BIN                    ! Loop counter
      integer :: L2PCINSTANCE           ! Instance index
      integer :: MAF1                   ! Subset limit
      integer :: MAFN                   ! Subset limit
      integer :: MYMAF                  ! A loop counter
      integer :: NOBINS                 ! Number of l2pc bins
      integer :: NOSELECTORS            ! A loop limit
      integer :: SIDEBAND               ! This sideband
      integer :: SELECTOR               ! Bin selector index
      integer :: SIGNAL                 ! Signal index
      integer :: STATEINSTANCE          ! Instance index
      integer :: STATUS                 ! Flag
      logical :: SPLIT                  ! Need a split calculation
      logical :: FOUNDINFIRST           ! Flag, not used
      real(r8) :: THISCOST              ! Interim cost

      integer, dimension(1) :: S1, S2   ! Results of minloc
      logical, dimension(:,:), pointer :: POSSIBLE ! Flags for each bin
      real(r8), dimension(:), pointer :: COST ! Cost of each bin
      real(r8), dimension(-1:1) :: BESTCOST
      type (QuantityTemplate_T), pointer :: BINRAD ! Quantity template
      type (BinSelector_T), pointer :: SEL ! One bin selector

      ! Executable code

      ! Firstly, if we're in locked bins mode, then just return the locked bin
      signal = radiance%template%signal
      sideband = radiance%template%sideband
      noBins = size ( l2pcDatabase )

      ! Now special code for dealing with the locked bins case
      if ( fmConf%lockBins ) then
        if ( .not. associated ( lockedBins ) ) then
          allocate ( lockedBins ( -1:1, size(signals) ), STAT=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Allocate//'lockedBins' )
          lockedBins = 0
        end if
        if ( any ( lockedBins ( :, signal ) /= 0 ) ) then
          l2pcBins = lockedBins ( :, signal )
          ! If got folded, use that by preference (3rd argument below).
          ! Assert that sidebands are present if needed here.
          call GetSidebandLoop ( signal, sideband, &
            & (lockedBins(0,signal)==0), sidebandStart, sidebandStop, &
            & sidebandStep )
          return
        end if
      end if

      ! Code will only get to here if the bins are unlocked, or will be locked
      ! once we've decided on them.

      ! Setup the arrays we need, possible is a flag to indicate (for
      ! each sideband) whether a bin is even worth considering, cost is
      ! the associated cost.
      nullify ( possible, cost )
      allocate ( possible ( -1:1, nobins ), STAT=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'possible' )
      call Allocate_test ( cost, noBins, 'cost', ModuleName )

      ! Get a first cut at the 'possible' flag
      possible = .false.
      do bin = 1, noBins
        binRad => l2pcDatabase(bin)%row%vec%quantities(1)%template
        possible ( binRad%sideband, bin ) = ( binRad%signal == signal )
      end do
      ! Set the cost to zero by default
      cost = 0.0_r8

      ! Setup some stuff we'll need from time to time
      if ( fmConf%lockBins ) then
        maf1 = 1 + radiance%template%noInstancesLowerOverlap
        mafN = radiance%template%noInstances - &
          & radiance%template%noInstancesUpperOverlap
      else
        maf1 = maf
        mafN = maf
      end if

      ! OK, now we have to loop over the bin selectors for each bin and
      ! apply the rules they give.
      noSelectors = size ( fmConf%binSelectors )
      do bin = 1, noBins
        binRad => l2pcDatabase(bin)%row%vec%quantities(1)%template
        do selector = 1, noSelectors
          sel => binSelectors ( fmConf%binSelectors(selector) )
          select case ( sel%selectorType )
          case ( l_nameFragment )
            call get_string ( l2pcDatabase(bin)%name, binName, strip=.true. )
            call get_string ( sel%nameFragment, nameFragment, strip=.true. )
            if ( len_trim(nameFragment) /= 0 ) then
              possible ( :, bin ) = index ( trim(binName), trim(nameFragment) ) /= 0
            end if
          case ( l_latitude )
            ! When we say latitude, we really mean an empirical phi.
            thisCost = sqrt ( sum ( &
              & (NormalizePhi(radiance%template%phi(1,maf1:mafN)) - &  
              &  NormalizePhi(binRad%phi(1,1)) )**2 ) / (mafN-maf1+1) ) / &
              & sel%cost
            if ( sel%exact ) then
              possible ( :, bin ) = EssentiallyEqual ( thisCost, 0.0_r8 )
            else
              cost(bin) = cost(bin) + thisCost
            end if
          case ( l_sza )
            thisCost = sqrt ( sum ( &
              & ( radiance%template%solarZenith(1,maf1:mafN) - &
              &   binRad%solarZenith(1,1) )**2 ) / &
              & (mafN-maf1+1) ) / sel%cost
            if ( sel%exact ) then
              possible ( :, bin ) = EssentiallyEqual ( thisCost, 0.0_r8 )
            else
              cost(bin) = cost(bin) + thisCost
            end if
          case ( l_temperature, l_vmr, l_fieldStrength, l_fieldElevation, l_fieldAzimuth )
            ! This one involves matching elements of xStar with x.
            if ( sel%selectorType == l_vmr ) then
              l2pcQ => GetVectorQuantityByType ( &
                & l2pcDatabase(bin)%col%vec, quantityType=l_vmr, &
                & molecule=sel%molecule, noError=.true. )
            else
              l2pcQ => GetVectorQuantityByType ( &
                & l2pcDatabase(bin)%col%vec, quantityType=sel%selectorType, &
                & noError=.true. )
            end if
            if ( .not. ValidateVectorQuantity ( l2pcQ,&
              & verticalCoordinate = (/ l_zeta /) ) ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Expected zeta coordinate for quantity in binSelector' )
            ! Find the relevant corresponding quantity in the state vector
            call FindMatchForL2PCQ ( l2pcQ, fmConf, fwdModelIn, &
              & fwdModelExtra, stateQ, foundInFirst )
            ! If we've got both of them make sure they match
            if ( associated(stateQ) .and. associated(l2pcQ) ) then
              ! OK, identify height range
              s1 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                & log10(binSelectors(selector)%heightRange(1) ) ) )
              s2 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                & log10(binSelectors(selector)%heightRange(2) ) ) )
              ! Here we'll just compare the central profile in the
              ! l2pc with the state profile closest to each maf.
              l2pcInstance = l2pcQ%template%noInstances/2 + 1
              thisCost = 0.0_r8
              do myMaf = maf1, mafN
                ! Only compare the closest profile to the maf
                stateInstance = FindOneClosestInstance ( stateQ, radiance, myMaf )
                thisCost = thisCost + sum ( &
                  & ( stateQ%values ( s1(1):s2(1), stateInstance ) - &
                  &   l2pcQ%values  ( s1(1):s2(1), l2pcInstance  ) ) **2 )
              end do
              ! If we require an exact match set the possible flag,
              ! otherwise just report our cost.
              if ( sel%exact ) then
                possible ( :, bin ) = EssentiallyEqual ( thisCost, 0.0_r8 )
              else
                cost ( bin ) = cost ( bin ) + sqrt ( thisCost / &
                  &  ( ( s2(1)-s1(1)+1 ) * ( mafN - maf1 + 1 ) ) ) / sel%cost
              end if
            end if
          end select                  ! Bin selector type
        end do                        ! Loop over selectors
      end do                          ! Loop over bins

      if ( index ( switches, 'binsel' ) /= 0 ) then
        call output ( 'Choosing bin for ' )
        call GetSignalName ( signal, namefragment, sideband=sideband )
        call output ( trim(nameFragment), advance='yes' )
        do bin = 1, noBins
          if ( any ( possible ( :, bin ) ) ) then
            call output  ( 'Candidate: ' )
            call display_string ( l2pcDatabase(bin)%name, strip=.true. )
            call output ( ' cost = ' )
            call output ( cost(bin), advance='yes' )
          end if
        end do
      end if

      ! OK, now we've surveyed the bins, let's cut things down.
      ! When computing folded radiances I'll always choose folded bins
      ! over unfolded ones, even if they cost more.  I think it's
      ! unlikely that this will really be an issue anyway.
      split = .false.
      if ( sideband == 0 ) then
        ! If we've got a match for the folded case, forget the others.
        if ( any ( possible(0,:) ) ) then
          where ( possible(0,:) ) 
            possible(-1,:) = .false.
            possible(1,:) = .false.
          end where
        else
          split = .true.
        end if
      end if

      ! Work out the range of the sideband loop
      call GetSidebandLoop ( signal, sideband, split, &
        & sidebandStart, sidebandStop, sidebandStep )

      ! Choose the bin(s)
      bestCost = huge ( cost(1) )
      l2pcBins = 0
      do bin = 1, noBins
        ! Reuse the sideband variable here, we don't need it's old
        ! definition anymore
        do sideband = sidebandStart, sidebandStop
          if ( possible(sideband,bin) .and. &
            & cost(bin) < bestCost(sideband) ) then
            bestCost(sideband) = cost(bin)
            l2pcBins(sideband) = bin
          end if
        end do
      end do

      ! Check that we've got the bins we need
      if (any(l2pcBins ((/sidebandStart,sidebandStop/)) == 0 )) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to find l2pc bins to match request' )

      ! Record this in the 'locked bins' information
      if ( fmConf%lockBins ) lockedBins ( :, signal ) = l2pcBins

      ! Tidy up
      call Deallocate_test ( cost, 'cost', ModuleName )
      deallocate ( possible )

    end subroutine SelectL2PCBins

  end subroutine LinearizedForwardModel

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module LinearizedForwardModel_m

! $Log$
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
