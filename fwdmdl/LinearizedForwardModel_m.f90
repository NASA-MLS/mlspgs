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
  use Intrinsic, only: L_RADIANCE, L_TEMPERATURE, L_PTAN, L_VMR, &
    & L_SIDEBANDRATIO
  use L2PC_m, only: L2PCDATABASE
  use ManipulateVectorQuantities, only: FINDCLOSESTINSTANCES
  use MatrixModule_0, only: M_ABSENT, M_BANDED, M_COLUMN_SPARSE, M_FULL, &
    & MATRIXELEMENT_T, CREATEBLOCK, DENSIFY
  use MatrixModule_1, only: MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP, &
    & FINDBLOCK, CREATEBLOCK
  use MLSCommon, only: r8
  use MLSSignals_m, only: Signal_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use MLSNumerics, only: HUNT, INTERPOLATEVALUES
  use Output_m, only: Output
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Display_String
  use Toggles, only: Emit, Levels, Toggle
  use Trace_m, only: Trace_begin, Trace_end
  use VectorsModule, only: assignment(=), OPERATOR(-), OPERATOR(+), &
    & CLONEVECTOR, CONSTRUCTVECTORTEMPLATE, COPYVECTOR, CREATEVECTOR,&
    & DESTROYVECTORINFO, GETVECTORQUANTITYBYTYPE, VECTOR_T, &
    & VECTORVALUE_T, VECTORTEMPLATE_T, DUMP

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

    ! The `prime' quantities are important.
    ! - yPrime (yP) is contains one maf of the relevant radiances, but on the
    !   l2pc's vertical coordinates
    ! - xPrime (xP) is like xStar but contains state values
    ! - kPrime is the jacobian for these two.

    real (r8), dimension(:), pointer :: lowerWeight ! For interpolation
    real (r8), dimension(:), pointer :: upperWeight ! For interpolation
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
    type(VectorValue_T), pointer :: XSTARPTAN ! Tangent pressure in l2pc
    type(VectorValue_T), pointer :: L2PCQ ! A quantity in the l2pc
    type(VectorValue_T), pointer :: SIDEBANDRATIO ! From the state vector

    ! Executable code
    if ( toggle(emit) ) call trace_begin ( 'LinearizedForwardModel' )

    nullify ( yPmapped, resultMapped, dyByDX, closestInstances )
    nullify ( dense, mifPointingsLower, mifPointingsUpper )
    nullify ( lowerWeight, upperWeight )
    nullify ( thisRatio, doChannel )

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
      & signal= signal%index, sideband=signal%sideband )

    ! Set some dimensions
    noChans = radiance%template%noChans
    noMIFs = radiance%template%noSurfs

    call allocate_test ( thisRatio, noChans, 'thisRatio', ModuleName )
    call allocate_test ( doChannel, noChans, 'doChannel', ModuleName )

    doChannel = .true.
    if ( associated ( signal%channels ) ) doChannel = signal%channels

    if ( signal%sideband == 0 ) then
      ! Use wants folded.  Can we find this in the l2pc file?
      ! Later in multi bin case think more carefully here.
      do l2pcBin = 1, size ( l2pcDatabase ) 
        radInl2pc => GetVectorQuantityByType ( &
          & l2pcDatabase(l2pcBin)%row%vec, quantityType=l_radiance, &
          & signal=signal%index, sideband=signal%sideband, noError=.true. )
        if ( associated ( radInL2PC ) ) exit
      end do
      if ( associated ( radInL2PC ) ) then
        ! We have a folded signal in the l2p
        thisRatio = 1.0
        sidebandStart = 0
        sidebandStop = 0
        sidebandStep = 1
      else
        ! We have to fold it ourselves
        sidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_sidebandRatio, signal=signal%index )
        thisRatio = sidebandRatio%values(:,1)
        sidebandStart = -1
        sidebandStop = 1
        sidebandStep = 2
      end if
    else
      thisRatio = 1.0
      sidebandStart = forwardModelConfig%signals(1)%sideband
      sidebandStop = sideBandStart
      sidebandStep = 1
    endif

    ! --------- Loop over sidebands ------------------------------------------------
    do sideband = sidebandStart, sidebandStop, sidebandStep

      ! Setup a sideband ratio array
      if ( sidebandStart /= sidebandStop ) then
        thisRatio = sidebandRatio%values(:,1)
        if ( sideband == 1 ) thisRatio = 1.0 - thisRatio
      end if

      ! Find an l2pc bin containing this signal.
      ! DO SOMETHING MORE CLEVER LATER!
      do l2pcBin = 1, size ( l2pcDatabase ) 
        radInl2pc => GetVectorQuantityByType ( &
          & l2pcDatabase(l2pcBin)%row%vec, quantityType=l_radiance, &
          & signal=signal%index, sideband=sideband, noError=.true. )
        if ( associated ( radInL2PC ) ) exit
      end do
      if ( .not. associated ( radInL2PC ) ) &
        & call MLSMessage( MLSMSG_Error, ModuleName, &
        &  "No appropriate l2pc bin found" )
      l2pc => l2pcDatabase(l2pcBin)

      ! Set a dimension
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

      ! -------- Main loop over xStar quantities -------------------------------
      do qtyInd = 1, size ( l2pc%col%vec%quantities )

        ! Identify this quantity in xStar
        l2pcQ => l2pc%col%vec%quantities(qtyInd)

        ! Now see if we are wanting to deal with this
        if ( l2pcQ%template%quantityType == l_ptan ) cycle ! Get this from interpolation
        if ( l2pcQ%template%quantityType == l_vmr ) then
          if ( .not. any &
            & (l2pcQ%template%molecule == forwardModelConfig%molecules)) cycle
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
        case default
          stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
            & quantityType = l2pcQ%template%quantityType, &
            & foundInFirst = foundInFirst, noError=.true. )
        end select

        ! If it's not in the state vector, perhaps make a fuss
        if ( .not. associated(stateQ) ) then
          if ( any(l2pcQ%template%quantityType == &
            & (/ l_temperature, l_vmr /))) &
            &   call MLSMessage ( MLSMSG_Error, ModuleName, &
            &  "Temperature or vmr quantity absent from state")
        end if

        ! Now check that surfaces are the same
        if (stateQ%template%noSurfs /= l2pcQ%template%noSurfs) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "State vector surfaces not same as l2pc" )
        ! The l2pc writer insisted that things were on zeta coords or none
        ! at all, so check surfaces OK.
        if ( any (abs (stateQ%template%surfs-l2pcQ%template%surfs) > 0.01)) &
          & call MLSMessage(MLSMSG_Error,ModuleName,&
          & 'State vector surfaces not same as l2pc')

        ! Check that no. chans. is the same (will it ever /=1!!?)
        if (stateQ%template%noChans /= l2pcQ%template%noChans) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "State vector channels not same as l2pc" )

        ! OK, we're legit, lets get going.
        instanceLen = l2pcQ%template%instanceLen
        call FindClosestInstances ( stateQ, radiance, closestInstances)
        closestInstance=closestInstances(maf)

        ! Do we need derivatives for this?
        doDerivatives = present(jacobian) .and. foundInFirst
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_temperature) ) &
          & doDerivatives = forwardModelConfig%temp_der
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_vmr) ) then
          doDerivatives = forwardModelConfig%atmos_der
          if ( doDerivatives .and. .not. any (l2pcQ%template%molecule == &
            & pack(ForwardModelConfig%molecules, &
            &      ForwardModelConfig%moleculeDerivatives))) &
            & doDerivatives = .false.
        end if

        ! If so, set up some arrays
        if (doDerivatives) then
          call allocate_test ( mifPointingsLower, noMIFs, &
            & 'mifPointingsLower', ModuleName )
          call allocate_test ( mifPointingsUpper, noMIFs, &
            & 'mifPointingsUpper', ModuleName )
          call allocate_test ( lowerWeight, noMIFs, &
            & 'lowerWeight', ModuleName )
          call allocate_test ( upperWeight, noMIFs, &
            & 'upperWeight', ModuleName )
          if ( l2pcQ%template%noSurfs /= 1 ) then
            call Hunt ( l2pcQ%template%surfs(:,1), ptan%values(:,maf), mifPointingsLower )
            mifPointingsUpper = mifPointingsLower + 1

            upperWeight = &
              & ( ptan%values(:,maf) - l2pcQ%template%surfs(mifPointingsLower,1) ) / &
              & ( l2pcQ%template%surfs(mifPointingsUpper,1) - &
              &   l2pcQ%template%surfs(mifPointingsLower,1) )
            lowerWeight = 1.0 - lowerWeight
            upperWeight = min ( 1.0_r8, max (upperWeight, 0.0_r8 ) )
            lowerWeight = min ( 1.0_r8, max (lowerWeight, 0.0_r8 ) )
          else
            mifPointingsLower = 1
            mifPointingsUpper = 1
            lowerWeight = 1.0
            upperWeight = 0.0
          end if
        end if

        ! Loop over profiles
        do xStarInstance = 1, l2pcQ%template%noInstances
          ! Identify this instance in state
          xInstance = closestInstance - &
            & l2pcQ%template%noInstances/2 + &
            & xStarInstance - 1
          xInstance = max ( 1, min ( xInstance, stateQ%template%noInstances ) )

          ! Fill this part of xP
          xP%quantities(qtyInd)%values(:,xStarInstance) = &
            & stateQ%values(:,xInstance)

          ! If the state vector part has a mask, copy it to DeltaX
! *************** MORE WORK NEEDED HERE NJL
!           if ( associated ( stateQ%mask ) ) then
!             if ( .not. associated ( deltaX%quantities(qtyInd)%mask ) ) &
!               & CreateMask ( deltaX%quantities(qtyInd) )
!             deltaX%mask(qtyInd)%values(:,xStarInstance) = &
!               stateQ%mask(:,xInstance)
!           endif

          ! If so, interpolate this block of kStar and place in jacobian
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
                      & thisRatio(chan) * ( &
                      &   lowerWeight(mif) * kBit( lower , : ) + &
                      &   upperWeight(mif) * kBit( upper, : ) )
                  endif
                  i = i + 1
                  lower = lower + 1
                  upper = upper
                end do
              end do

              if ( any ( l2pcBlock%kind == &
                & (/ M_Column_Sparse, M_Banded /) ) ) &
                & call deallocate_test ( dense, 'dense', ModuleName )
            end if                        ! Any matrix here?
          end if                          ! do derivatives?
        end do                            ! End loop over xStar Profiles

        if ( doDerivatives ) then         ! Destroy working arrays
          call deallocate_test ( mifPointingsLower, &
            & 'mifPointingsLower', ModuleName )
          call deallocate_test ( mifPointingsUpper, &
            & 'mifPointingsUpper', ModuleName )
          call deallocate_test ( lowerWeight, &
            & 'lowerWeight', ModuleName )
          call deallocate_test ( upperWeight, &
            & 'upperWeight', ModuleName )
        end if
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

      ! Now compute yP

      if ( toggle(emit) .and. levels(emit) > 2 ) then
        call dump ( (/deltaX/) )

        call dump ( l2pc%col%inst, 'l2pc%col%inst' )
        call dump ( l2pc%col%quant, 'l2pc%col%quant' )

        call dump ( l2pc%row%inst, 'l2pc%row%inst' )
        call dump ( l2pc%row%quant, 'l2pc%row%quant' )
      end if

      call cloneVector( yp, l2pc%row%vec, vectorNameText='_yP' )
      call MultiplyMatrixVectorNoT ( l2pc, deltaX, yP, update = .false. )

      if ( toggle(emit) .and. levels(emit) > 2 ) then
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
        & extrapolate='Constant', &       ! dont extrapolate, clamp
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
            &   noMIFs*noChans )
          jBlock%values = 0.0_r8
          do mif = 1, noMIFs
            jBlock%r1(mif) = 1 + noChans * (mif-1)
            jBlock%r2(mif) = noChans * mif
          end do
        case (m_banded)
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Wrong matrix type for ptan derivative')
        end select
        
        if ( sidebandStart == sidebandStop ) then
          jBlock%values(:,1) = reshape(transpose(dyByDx),&
            & (/radiance%template%instanceLen/))
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
      
      call DestroyVectorInfo ( xP )
      call DestroyVectorInfo ( deltaX )
      call Deallocate_test ( closestInstances, 'closestInstances', ModuleName )
      call Deallocate_test ( dyByDx, 'dyByDx', ModuleName)
      call Deallocate_test ( resultMapped, 'resultMapped', ModuleName )
      call Deallocate_test ( ypMapped, 'ypMapped', ModuleName )

    end do                                ! End of sideband loop

    if ( maf == radiance%template%noInstances ) fmStat%finished = .true.

    if ( toggle(emit) ) call trace_end ( 'LinearizedForwardModel' )

  end subroutine LinearizedForwardModel
end module LinearizedForwardModel_m

! $Log$
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
