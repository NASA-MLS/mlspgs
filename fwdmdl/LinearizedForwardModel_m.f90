! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module LinearizedForwardModel_m

  ! Approximate a forward model with a first-order Taylor's series.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T, &
    & FORWARDMODELINTERMEDIATE_T
  use Intrinsic, only: LIT_INDICES
  use Intrinsic, only: L_NONE, L_RADIANCE, L_TEMPERATURE, L_PTAN, L_VMR, &
    & L_SIDEBANDRATIO, L_ZETA, L_OPTICALDEPTH, L_LATITUDE
  use L2PC_m, only: L2PCDATABASE, BINSELECTORS
  use ManipulateVectorQuantities, only: FINDONECLOSESTINSTANCE
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T, CREATEBLOCK, DENSIFY
  use MatrixModule_1, only: MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP, &
    & FINDBLOCK, CREATEBLOCK
  use MLSCommon, only: r8
  use MLSSignals_m, only: Signal_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use MLSNumerics, only: HUNT, INTERPOLATEVALUES
  use Molecules, only: L_EXTINCTION
  use Output_m, only: Output
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Display_String
  use Toggles, only: Emit, Levels, Toggle
  use Trace_m, only: Trace_begin, Trace_end
  use VectorsModule, only: assignment(=), OPERATOR(-), OPERATOR(+), &
    & CLONEVECTOR, CONSTRUCTVECTORTEMPLATE, COPYVECTOR, CREATEVECTOR,&
    & DESTROYVECTORINFO, GETVECTORQUANTITYBYTYPE, VECTOR_T, &
    & VECTORVALUE_T, VECTORTEMPLATE_T, DUMP, &
    & VALIDATEVECTORQUANTITY

  implicit none
  private
  public :: LinearizedForwardModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! -------------------------------------  LinearizedForwardModel  -----
  subroutine LinearizedForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, Ifm, fmStat, Jacobian )
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Local variables
    integer :: CHAN                     ! Loop counter
    integer :: CLOSESTINSTANCE          ! A closest profile to this MAF
    integer :: COLJBLOCK                ! Column index in jacobian
    integer :: COLLBLOCK                ! Column index in l2pc
    integer :: DELTAINSTANCE            ! Instance offset between center an current
    integer :: I                        ! Array index
    integer :: INSTANCELEN              ! For the state quantity
    integer :: L2PCBIN                  ! Which l2pc to use
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

    logical :: DODERIVATIVES            ! Flag
    logical :: FOUNDINFIRST             ! Flag for state quantities
    logical :: PTANINFIRST              ! PTan was found in the first vector

    integer, dimension(:), pointer :: closestInstances ! Of qty to radiance
    integer, dimension(:), pointer :: mifPointingsLower ! Result of a hunt
    integer, dimension(:), pointer :: mifPointingsUpper ! mifPointingsLower+1

    logical, dimension(:), pointer :: doChannel ! Do this channel?

    real (r8) :: BESTCOST               ! Output from SelectL2PCBin

    ! The `prime' quantities are important.
    ! - yPrime (yP) is contains one maf of the relevant radiances, but on the
    !   l2pc's vertical coordinates
    ! - xPrime (xP) is like xStar but contains state values
    ! - kPrime is the jacobian for these two.

    real (r8), dimension(:), pointer :: lowerWeight ! For interpolation
    real (r8), dimension(:), pointer :: upperWeight ! For interpolation
    real (r8), dimension(:), pointer :: tangentTemperature ! For optical depth
    real (r8), dimension(:), pointer :: thisRatio ! Sideband ratio values

    real (r8), dimension(:,:), pointer :: yPmapped ! Remapped values of yP
    real (r8), dimension(:,:), pointer :: resultMapped ! Remapped values of result
    real (r8), dimension(:,:), pointer :: dyByDX ! Raw dRad/dPtan
    real (r8), dimension(:,:), pointer :: dense  ! Densified matrix from l2pc
    real (r8), dimension(:,:), pointer :: kBit ! Remapped values of l2pc

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
    type(VectorValue_T), pointer :: SIDEBANDRATIO ! From the state vector
    type(VectorValue_T), pointer :: LOWERSIDEBANDRATIO ! From the state vector
    type(VectorValue_T), pointer :: UPPERSIDEBANDRATIO ! From the state vector

    ! Executable code
    if ( toggle(emit) ) call trace_begin ( 'LinearizedForwardModel' )

    nullify ( yPmapped, resultMapped, dyByDX, closestInstances )
    nullify ( dense, mifPointingsLower, mifPointingsUpper )
    nullify ( lowerWeight, upperWeight )
    nullify ( thisRatio, doChannel )
    nullify ( tangentTemperature )

    ! Identify the band/maf we're looking for

    maf = fmStat%maf

    if ( toggle(emit) .and. levels(emit) > 0 ) then
      call output ( 'Linear model doing maf: ' )
      call output ( maf, advance='yes' )
    end if
    if ( size ( forwardModelConfig%signals ) /= 1 ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Can only have one signal for linearized models')

    signal = forwardModelConfig%signals(1)
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

    ! Set some dimensions
    noChans = radiance%template%noChans
    noMIFs = radiance%template%noSurfs

    call allocate_test ( thisRatio, noChans, 'thisRatio', ModuleName )
    call allocate_test ( doChannel, noChans, 'doChannel', ModuleName )

    doChannel = .true.
    if ( associated ( signal%channels ) ) doChannel = signal%channels

    if ( signal%sideband == 0 ) then
      ! User wants folded.  Can we find this in the l2pc file?
      ! If we have a choice between unfoled or folded l2pc bins,
      ! we'll always choose the folded one for speed, and only think about bin
      ! matching later on.  I imagine most of the time we won't have 'mixed'
      ! l2pc files.
      do l2pcBin = 1, size ( l2pcDatabase ) 
        radInl2pc => GetVectorQuantityByType ( &
          & l2pcDatabase(l2pcBin)%row%vec, quantityType=l_radiance, &
          & signal=signal%index, sideband=signal%sideband, noError=.true. )
        if ( associated ( radInL2PC ) ) exit
      end do
      if ( associated ( radInL2PC ) ) then
        ! We have a folded signal in the l2pc
        thisRatio = 1.0
        sidebandStart = 0
        sidebandStop = 0
        sidebandStep = 1
      else
        ! We have to fold it ourselves
        sidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_sidebandRatio, signal=signal%index, noError=.true. )
        lowerSidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_sidebandRatio, signal=signal%index, &
          & sideband=-1, noError=.true. )
        upperSidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_sidebandRatio, signal=signal%index, &
          & sideband=1, noError=.true. )
        if (.not. associated (sidebandRatio) .and. .not. &
          & ( associated ( lowerSidebandRatio) .and. associated ( upperSidebandRatio ) ) ) &
          & call MLSMessage(MLSMSG_Error,ModuleName, &
          & "No sideband ratio supplied")
        if ( signal%singleSideband == 0 ) then
          ! This is not a single sideband radiometer
          sidebandStart = -1
          sidebandStop = 1
          sidebandStep = 2
        else
          ! This is a single sideband radiometer
          sidebandStart = signal%singleSideband
          sidebandStop = sidebandStart
          sidebandStep = 1
        end if
      end if
    else
      sidebandStart = forwardModelConfig%signals(1)%sideband
      sidebandStop = sideBandStart
      sidebandStep = 1
    endif

    ! --------- Loop over sidebands ------------------------------------------------
    do sideband = sidebandStart, sidebandStop, sidebandStep

      ! Setup a sideband ratio array
        if ( sidebandStart /= sidebandStop ) then   ! We're folding
          if ( associated ( sidebandRatio ) ) then
            thisRatio = sidebandRatio%values(:,1)
            if ( sideband == 1 ) thisRatio = 1.0 - thisRatio
          else
            if ( sideband == -1 ) then
              thisRatio = lowerSidebandRatio%values(:,1)
            else
              thisRatio = upperSidebandRatio%values(:,1)
            end if
          end if
        else                  ! Otherwise, want just unfolded signal
          thisRatio = 1.0
        end if

      ! Work out which l2pc bin we want
      l2pcBin = SelectL2PCBin ( FwdModelIn, FwdModelExtra, radiance, &
        & sideband, maf, bestCost )
      if ( l2pcBin == 0 ) &
        & call MLSMessage( MLSMSG_Error, ModuleName, &
        &  "No appropriate l2pc bin found" )
      l2pc => l2pcDatabase(l2pcBin)
      if ( toggle(emit) .and. levels(emit) > 0 ) then
        call output ( 'Choosing l2pc bin: ' )
        call output ( l2pcBin )
        call output ( ' Cost was: ' )
        call output ( bestCost, advance='yes' )
      end if

      ! Set a dimension
      radInl2pc => GetVectorQuantityByType ( &
        & l2pcDatabase(l2pcBin)%row%vec, quantityType=l_radiance, &
        & signal=signal%index, sideband=sideband, noError=.true. )
      noPointings = radInL2PC%template%noSurfs
      if ( radInL2PC%template%noChans /= noChans ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Channel dimension in l2pc not same as in measurements" )

      ! Now we loop over the quantities in the l2pc file and construct an xPrime
      ! for them
      call cloneVector ( xP, l2pc%col%vec, vectorNameText='_xP' )
      call cloneVector ( deltaX, xP, vectorNameText='_deltaX' ) ! sets values to 0.0

      ! Set up some other stuff before main loop
      call Allocate_test ( closestInstances, radiance%template%noInstances, &
        & 'closestInstances', ModuleName )

      ! Get the two ptans, we'll need these for interpolation
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
        call allocate_test ( mifPointingsLower, noMIFs, &
          & 'mifPointingsLower', ModuleName )
        call allocate_test ( mifPointingsUpper, noMIFs, &
          & 'mifPointingsUpper', ModuleName )
        call allocate_test ( lowerWeight, noMIFs, &
          & 'lowerWeight', ModuleName )
        call allocate_test ( upperWeight, noMIFs, &
          & 'upperWeight', ModuleName )

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
      do qtyInd = 1, size ( l2pc%col%vec%quantities )

        if ( toggle(emit) .and. levels(emit) > 1 ) then
          call output ( 'Dealing with xStar Quantity named ' )
          call display_string ( l2pc%col%vec%quantities(qtyInd)%template%name, &
            & advance='yes' )
        end if

        ! Identify this quantity in xStar
        l2pcQ => l2pc%col%vec%quantities(qtyInd)

        ! Now see if we are wanting to deal with this
        if ( l2pcQ%template%quantityType == l_ptan ) cycle ! Get this from interpolation
        if ( l2pcQ%template%quantityType == l_vmr ) then
          if (.not. associated(forwardModelConfig%molecules) ) cycle
          if ( .not. any (l2pcQ%template%molecule == &
            &   forwardModelConfig%molecules)) cycle
        end if

        ! Identify this quantity in x
        select case ( l2pcQ%template%quantityType )
        case ( l_temperature )
          stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
            & quantityType = l_temperature, &
            & foundInFirst = foundInFirst, noError=.true. )
        case ( l_vmr )
          stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
            & quantityType = l_vmr, &
            & molecule = l2pcQ%template%molecule, &
            & foundInFirst = foundInFirst, noError=.true. )
        case ( l_extinction )
!           stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
!             & quantityType = l_extinction, &
!             & radiometer = l2pcQ%template%radiometer, &
!             & foundInFirst = foundInFirst, noError=.true. )

          ! Temporary
          stateQ => NULL()
        case default
          ! For the moment, just ignore things we don't understand.
          stateQ => NULL()
!           stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
!             & quantityType = l2pcQ%template%quantityType, &
!             & foundInFirst = foundInFirst, noError=.true. )
        end select

        ! If it's not in the state vector, perhaps make a fuss
        if ( .not. associated(stateQ) ) then
          if ( any(l2pcQ%template%quantityType == &
            & (/ l_temperature, l_vmr /))) &
            &   call MLSMessage ( MLSMSG_Error, ModuleName, &
            &  "Temperature or vmr quantity absent from state")
          cycle                         ! Go to next l2pc quantity
        end if

        ! Now check that surfaces are the same
        if (stateQ%template%noSurfs /= l2pcQ%template%noSurfs) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Number of state vector surfaces not same as l2pc" )
        ! The l2pc writer insisted that things were on zeta coords or none
        ! at all, so check surfaces OK.
        if ( any (abs (stateQ%template%surfs-l2pcQ%template%surfs) > 0.01)) &
          & call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'State vector surface values not same as l2pc')

        ! Check that no. chans. is the same (will it ever /=1!!?)
        if (stateQ%template%noChans /= l2pcQ%template%noChans) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "State vector channels not same as l2pc" )

        ! OK, we're legit, lets get going.
        instanceLen = l2pcQ%template%instanceLen
        closestInstance = FindOneClosestInstance ( stateQ, radiance, maf )

        ! Do we need derivatives for this?
        doDerivatives = present(jacobian) .and. foundInFirst
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_temperature) ) &
          & doDerivatives = forwardModelConfig%temp_der
        if (doDerivatives .and. &
          & (l2pcQ%template%quantityType==l_vmr) .and. &
          & associated(forwardModelConfig%molecules) ) then
          doDerivatives = forwardModelConfig%atmos_der
          if ( doDerivatives .and. .not. any (l2pcQ%template%molecule == &
            & pack(ForwardModelConfig%molecules, &
            &      ForwardModelConfig%moleculeDerivatives))) &
            & doDerivatives = .false.
        end if

        ! Loop over profiles
        do xStarInstance = 1, l2pcQ%template%noInstances
          ! Identify this instance in state
          deltaInstance = xStarInstance - l2pcQ%template%noInstances/2 - 1
          deltaInstance = max ( -forwardModelConfig%phiWindow/2, &
            &                   min ( deltaInstance, forwardModelConfig%phiWindow/2) )
          xInstance = closestInstance + deltaInstance
          xInstance = max ( 1, min ( xInstance, stateQ%template%noInstances ) )

          ! Fill this part of xP
          xP%quantities(qtyInd)%values(:,xStarInstance) = &
            & stateQ%values(:,xInstance)

          ! Note that we've ignored the mask in stateQ here. We might
          ! in later versions want to apply one mask field?

          ! Interpolate kStart to K if required.
          if (doDerivatives) then
            fmStat%rows(rowJBlock) = .true.
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
              endif

              ! Indentify the block in the jacobian
              jBlock => jacobian%block(rowJblock,colJblock)
              select case ( jBlock%kind )
              case ( M_Absent )
                call CreateBlock ( jBlock, &
                  & noChans*noMIFs, instanceLen, &
                  & M_Full )
                jBlock%values = 0.0_r8
              case ( M_Banded, M_Column_Sparse )
                call MLSMessage( MLSMSG_Error, ModuleName, &
                  & "Not written code for adding to non full blocks" )
              case default
              end select

              ! Now do the interpolation
              i = 1
              do mif = 1, noMIFs
                lower = (mifPointingsLower(mif)-1)*noChans + 1
                upper = (mifPointingsUpper(mif)-1)*noChans + 1
                  do chan = 1, noChans
                    if ( doChannel(chan) ) then
                      jBlock%values ( i , : ) = &
                        & jBlock%values ( i , : ) + &
                        &   thisRatio(chan) * ( &
                        &     lowerWeight(mif) * kBit( lower, : ) + &
                        &     upperWeight(mif) * kBit( upper, : ) )
                    endif
                    i = i + 1
                    lower = lower + 1
                    upper = upper + 1
                  end do
              end do

              if ( any ( l2pcBlock%kind == &
                & (/ M_Column_Sparse, M_Banded /) ) ) &
                & call deallocate_test ( dense, 'dense', ModuleName )
            end if                        ! Any matrix here?
          end if                          ! do derivatives?
        end do                            ! End loop over xStar Profiles

        ! Compute this instance of deltaX
        deltaX%quantities(qtyInd)%values = &
          & xP%quantities(qtyInd)%values - &
          & l2pc%col%vec%quantities(qtyInd)%values

        ! NOTE FOR FUTURE ******************************************
        ! Now make sure that the deltaX's where mask is set are zero
        ! Write this bit later!!!!! !?????? ***********************
        ! Also, do we need to zero out the corresponding columns of K.
        ! I think not as Van's code skips them already.

      end do                              ! End loop over quantities

      if ( present (jacobian) ) then         ! Destroy working arrays
        call deallocate_test ( mifPointingsLower, &
          & 'mifPointingsLower', ModuleName )
        call deallocate_test ( mifPointingsUpper, &
          & 'mifPointingsUpper', ModuleName )
        call deallocate_test ( lowerWeight, &
          & 'lowerWeight', ModuleName )
        call deallocate_test ( upperWeight, &
          & 'upperWeight', ModuleName )
      end if

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
      yP = yP + l2pc%row%vec

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
        & xStarPtan%values(:,1), &        ! OldX
        & yPmapped, &                     ! OldY
        & ptan%values(:,maf), &           ! NewX
        & resultMapped, &                 ! NewY
        & 'Spline', &                     ! use spline
        & extrapolate='Allow', &          ! Allow extrapolation in radiance
        & dyByDx=dyByDx )

      if ( sidebandStart == sidebandStop ) then
        radiance%values(:,maf) = reshape(transpose(resultMapped),&
          & (/radiance%template%instanceLen/))
      else
        ! Either place or add, make decision outside loop!
        if ( sideband == sidebandStart ) then
          do chan = 1, noChans
            if ( doChannel ( chan ) ) then
              do mif = 1, noMIFs
                radiance%values( chan + (mif-1)*noChans, maf ) = &
                  & thisRatio(chan)*resultMapped ( mif, chan )
              end do
            endif
          enddo
        else
          do chan = 1, noChans
            if ( doChannel ( chan ) ) then
              do mif = 1, noMIFs
                radiance%values( chan + (mif-1)*noChans, maf ) = &
                  & radiance%values( chan + (mif-1)*noChans, maf ) + &
                  &   thisRatio(chan)*resultMapped ( mif, chan )
              end do
            endif
          enddo
        end if
      end if                            ! Doing folded
        
      ! Now deal with ptan derivatives
      if ( present( Jacobian ) .and. ptanInFirst ) then
        colJBlock = FindBlock ( Jacobian%col, ptan%index, maf )
        jBlock => Jacobian%block( rowJBlock, colJBlock )
        select case (jBlock%kind)
        case (m_absent)
          call CreateBlock ( Jacobian, rowJBlock, colJBlock, m_banded, &
            & noMIFs*noChans, bandHeight=noChans )
          jBlock%values = 0.0_r8
        case (m_banded)
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Wrong matrix type for ptan derivative')
        end select
        
        if ( sidebandStart == sidebandStop ) then
          do chan = 1, noChans
            if ( doChannel ( chan ) ) then
              do mif = 1, noMIFs
                jBlock%values( chan + (mif-1)*noChans, 1 ) = &
                  & dyByDx ( mif, chan )
              end do                  ! Minor frames
            endif                     ! Doing this channel
          enddo                       ! Channels
        else
          ! Either place or add, make decision outside loop!
          if ( sideband == sidebandStart ) then
            do chan = 1, noChans
              if ( doChannel ( chan ) ) then
                do mif = 1, noMIFs
                  jBlock%values( chan + (mif-1)*noChans, 1 ) = &
                    & thisRatio(chan) * dyByDx ( mif, chan )
                end do                  ! Minor frames
              endif                     ! Doing this channel
            enddo                       ! Channels
          else                          ! Must be doing second sideband
            do chan = 1, noChans
              if ( doChannel ( chan ) ) then
                do mif = 1, noMIFs
                  jBlock%values( chan + (mif-1)*noChans, : ) = &
                    & jBlock%values ( chan + (mif-1)*noChans, : ) + &
                    &   thisRatio(chan) * dyByDx ( mif, chan )
                end do                  ! Minor frames
              endif                     ! Doing this channel
            enddo                       ! Channel
          end if                        ! First/second sideband
        end if                          ! Folding
      end if                            ! Want ptan derivatives

      ! Now think about optical depth
      if ( radiance%template%quantityType == l_opticalDepth ) then
        ! Get temperature and interpolate to tangent points
        temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_temperature )
        closestInstance = FindOneClosestInstance ( temp, radiance, maf )
        call Allocate_test ( tangentTemperature, noMIFs, &
          & 'tangentTemperature', ModuleName )
        call InterpolateValues ( &
          & temp%template%surfs(:,1), &
          & temp%values(:,closestInstance), &
          & ptan%values(:,maf), &
          & tangentTemperature, &
          & 'Spline', extrapolate='Allow' )
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
      call Deallocate_test ( closestInstances, 'closestInstances', ModuleName )
      call Deallocate_test ( dyByDx, 'dyByDx', ModuleName)
      call Deallocate_test ( resultMapped, 'resultMapped', ModuleName )
      call Deallocate_test ( ypMapped, 'ypMapped', ModuleName )

    end do                                ! End of sideband loop

    if ( toggle(emit) ) call trace_end ( 'LinearizedForwardModel' )

  end subroutine LinearizedForwardModel

  ! ======================================= Private procudures

  ! ------------------------------------- SelectL2PCBin -----
  integer function SelectL2PCBin ( FwdModelIn, FwdModelExtra, &
    & radiance, sideband, maf, bestCost )
    ! Dummy arguments
    type (Vector_T), intent(in) :: FWDMODELIN ! Main state vector
    type (Vector_T), intent(in) :: FWDMODELEXTRA ! Other state vector
    type (VectorValue_T), intent(in) :: RADIANCE ! The radiance we're after
    integer, intent(in) :: SIDEBAND     ! What sideband of that radiance
    integer, intent(in) :: MAF          ! What maf for that radiance
    real (r8), intent(out), optional :: BESTCOST ! Best cost

    ! Local variables
    integer :: L2PCBIN                  ! Loop counter
    integer :: SELECTOR                 ! Loop counter
    integer :: S1(1), S2(1)             ! Surface range
    type (VectorValue_T), pointer :: RADINL2PC ! Radiance vector quantity
    type (VectorValue_T), pointer :: STATEQ ! Radiance vector quantity
    type (VectorValue_T), pointer :: L2PCQ ! Radiance vector quantity
    real(r8) :: COSTS(size(L2PCDatabase)) ! Cost for each bin
    logical :: FLAG                     ! Flag to simplify an if
    integer :: RESULTASARRAY(1)         ! Result
    integer :: L2PCINSTANCE             ! Instance index in l2pc
    integer :: SIGNAL                   ! Signal index
    integer :: STATEINSTANCE            ! Instance index in state vector
    logical :: APPROPRIATE(size(L2PCDatabase)) ! Is this bin appropriate
    real (r8) :: MEANRADPHI             ! Mean value of phi for this maf
    real (r8) :: MEANL2PCPHI            ! Mean value of phi in center scan of bin

    ! Executable code
    costs = 0.0_r8
    appropriate = .true.
    signal = radiance%template%signal
    do l2pcBin = 1, size(l2pcDatabase)
      radInL2PC => GetVectorQuantityByType ( &
        & l2pcDatabase(l2pcBin)%row%vec, quantityType=l_radiance, &
        & signal=signal, sideband=sideband, noError=.true. )
      if ( associated(radInL2PC) ) then
        ! OK this signal is present in the l2pc file.
        ! Now loop over the bin selectors and apply the relevant rules
        if ( associated(binSelectors) ) then
          do selector = 1, size(binSelectors)
            ! Only deal with the selectors for this band etc.
            if ( any ( &
              & ( binSelectors(selector)%signals == signal ) .and. &
              & ( ( binSelectors(selector)%sidebands == 0 ) .or. &
              &   ( binSelectors(selector)%sidebands == sideband ) ) ) ) then

              ! Some bin selectors are based on quantities, other on positions
              ! Do the position based ones first

              if ( binSelectors(selector)%selectorType == l_latitude ) then
                ! --------------------------- Position based selectors
                ! When we say latitude, we really mean phi 
                ! (at least for the moment)
                ! WORRY LATER ABOUT MISSING MINOR FRAMES?
                l2pcInstance = radInL2PC%template%noInstances/2 + 1
                meanRadPhi = sum ( radiance%template%phi(:,maf) ) / &
                  & radiance%template%noSurfs
                meanL2PCPhi = sum ( radInL2PC%template%phi(:,l2pcInstance) ) / &
                  & radInL2PC%template%noSurfs
                costs(l2pcBin) = costs(l2pcBin) + &
                  & binSelectors(selector)%cost * abs ( &
                  &   phiToLat(meanRadPhi) - phiToLat(meanL2PCPhi) )
              else
                ! --------------------------- Value based selectors
                ! Locate the quantities in the state and l2pc xStar vectors
                select case ( binSelectors(selector)%selectorType )
                case ( l_temperature )
                  stateQ => GetVectorQuantityByType ( &
                    & fwdModelIn, fwdModelExtra, quantityType=l_temperature, &
                    & noError=.true. )
                  l2pcQ => GetVectorQuantityByType ( &
                    & l2pcDatabase(l2pcBin)%col%vec, quantityType=l_temperature, &
                    & noError=.true. )
                case ( l_vmr )
                  stateQ => GetVectorQuantityByType ( &
                    & fwdModelIn, fwdModelExtra, quantityType=l_vmr, &
                    & molecule=binSelectors(selector)%molecule, noError=.true. )
                  l2pcQ => GetVectorQuantityByType ( &
                    & l2pcDatabase(l2pcBin)%col%vec, quantityType=l_vmr, &
                    & molecule=binSelectors(selector)%molecule, noError=.true. )
                case default
                  call MLSMessage ( MLSMSG_Error, ModuleName, &
                    & 'Unexpected problem with bin selection' )
                end select
                ! If we've got both of them make sure they match
                if ( associated(stateQ) .and. associated(l2pcQ) ) then
                  flag = &
                    & ValidateVectorQuantity ( stateQ, coherent=.true., &
                    & verticalCoordinate=(/l_zeta/), frequencyCoordinate=(/l_none/) )
                  flag = flag .and. &
                    & ValidateVectorQuantity ( l2pcQ, coherent=.true., &
                    & verticalCoordinate=(/l_zeta/), frequencyCoordinate=(/l_none/) )
                  if ( flag ) flag = flag .and. &
                    & (stateQ%template%noSurfs == l2pcQ%template%noSurfs)
                  if ( flag ) flag = flag .and. &
                    & all (abs (stateQ%template%surfs-l2pcQ%template%surfs) < 0.01)
                  if ( flag ) then
                    ! OK, identify height range
                    s1 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                      & log10(binSelectors(selector)%heightRange(1) ) ) )
                    s2 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                      & log10(binSelectors(selector)%heightRange(2) ) ) )
                    ! We're at last ready to go.  Just compare central profiles in
                    ! each l2pc bin with each state vector profile
                    l2pcInstance = l2pcQ%template%noInstances/2 + 1
                    ! Only compare the closest profile to the maf
                    stateInstance = FindOneClosestInstance ( stateQ, radiance, maf )
                    costs ( l2pcBin ) = costs ( l2pcBin ) + sum ( abs ( &
                      & stateQ%values ( s1(1):s2(1), stateInstance ) - &
                      & l2pcQ%values ( s1(1):s2(1), l2pcInstance ) ) ) / &
                      & ( binSelectors(selector)%cost * ( s2(1)-s1(1)+1 ) )
                  else
                    ! This l2pcBin doesn't match the quantity dimensions etc.
                    appropriate ( l2pcBin ) = .false.
                  end if                ! OK to do the test
                end if                  ! Quantity in both state and l2pc
              end if                    ! A quantity based selector
            end if                      ! This selector relevant
          end do                        ! End loop over selectors
        end if                          ! Any selectors
      else                              ! This l2pc bin not relevant
        appropriate ( l2pcBin ) = .false.
      end if
    end do                              ! End loop over l2pc bins
    if ( any ( appropriate ) ) then 
      where ( .not. appropriate )
        costs = huge ( 0.0_r8)
      end where
      resultAsArray = minloc ( costs, appropriate )
      SelectL2PCBin = resultAsArray(1)
      if ( present(bestCost) ) bestCost = costs( resultAsArray(1) )
    else
      SelectL2PCBin = 0
      if ( present(bestCost) ) bestCost = huge( costs(1) )
    end if

  contains
    
    ! This does an approximate phi to latitude conversion
    real (r8) function PhiToLat ( phi ) result ( lat )
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
    end function Phitolat

  end function SelectL2PCBin
end module LinearizedForwardModel_m

! $Log$
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
