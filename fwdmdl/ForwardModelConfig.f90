! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelConfig
!=============================================================================

! Set up the forward model configuration, except for actually processing
! the command.


  use MLSCommon, only: R8, RP
  use MLSSignals_M, only: Signal_T
  use VGridsDatabase, only: VGrid_T, DestroyVGridContents

  implicit NONE
  private

  ! Public procedures:
  interface Dump
    module procedure Dump_ForwardModelConfig, Dump_ForwardModelConfigDatabase
  end interface Dump

  public :: AddForwardModelConfigToDatabase, DeriveFromForwardModelConfig
  public :: DestroyFWMConfigDatabase, DestroyForwardModelDerived, Dump
  public :: NullifyForwardModelConfig
  public :: StripForwardModelConfigDatabase, PVMPackFWMConfig, PVMUnpackFWMConfig
 
  ! Public Types:

  ! Quantities derived from forward models, but not carted around by
  ! PVMPackFWMConfig and PVMUnpackFWMConfig.  Rather, they are computed
  ! by ForwardModelDerive when a ForwardModelConfig_T is created, or
  ! when arrives by way of PVMUnpackFWMConfig

  ! Channel information from the signals database
  type, public :: Channels_T
    integer :: Used       ! Which channel is this?
    integer :: Origin     ! Index of first channel (zero or one)
    integer :: Signal     ! Signal index for the channel
    integer :: DACS       ! DACS index if any, else zero
    integer, pointer :: PFAData(:,:) => NULL()    ! Sidebands X numPFA
                                                  ! Indices in PFADataBase%PFAData
    integer, pointer :: PFAMolecules(:) => NULL() ! L_... from PFAData for this channel
  end type Channels_T

  ! Now all of the derived stuff
  type, public :: ForwardModelDerived_T
    real(rp), dimension(:,:), pointer :: DACsStaging => NULL() ! Temporary
                                         ! space for DACS radiances
    integer, dimension(:), pointer :: USEDDACSSIGNALS => NULL() ! Indices in
                                         ! FwdModelConf_T%Signals
                                         ! of signals for our dacs
    type(channels_T), pointer, dimension(:) :: Channels => NULL()
  end type ForwardModelDerived_T
  
  ! These components are sorted in the order they are to make the packing
  ! and unpacking for PVM as easy as possible to maintain
  type, public :: ForwardModelConfig_T
    ! First the lit_indices
    integer :: Name                   ! String index of config name
    integer :: Cloud_der              ! Compute cloud sensitivity in cloud models.
                                      ! l_iwc_low_height, l_iwc_high_height, l_iwp
                                      ! l_none
    integer :: FwmType                ! l_linear, l_full or l_scan
    integer :: I_saturation           ! Flag to determine saturation status
                                      ! l_clear, l_clear_110rh_below_top
                                      ! l_clear_0rh, l_clear_lowest_0_110rh
                                      ! l_clear_110rh_below_tropopause,
                                      ! l_cloudy_110rh_below_top
                                      ! l_cloudy_110rh_in_cloud,
                                      ! l_cloudy_nearside_only
    integer :: InstrumentModule       ! Module for scan model (actually a spec index)
    integer :: LinearSideband         ! For hybrid model, which SB is linear?
    ! Now the other integers
    integer :: FirstPFA               ! Index of first PFA in Molecules
    integer :: No_cloud_species       ! No of Cloud Species '2'
    integer :: No_model_surfs         ! No of Model surfaces '640'
    integer :: Ntimes = 0	      ! Number of times calling FullForwardModel
    integer :: Num_ab_terms           ! No of AB terms '50'
    integer :: Num_azimuth_angles     ! No of azmuth angles '8'
    integer :: Num_scattering_angles  ! No of scattering angles '16'
    integer :: Num_size_bins          ! No of size bins '40'
    integer :: SidebandStart, SidebandStop ! Folded or SSB config?
    integer :: SurfaceTangentIndex    ! Index in Tangentgrid of Earth's surface
    integer :: WindowUnits            ! Either degrees or profiles
    integer :: xStar                  ! Index of specific vector to use for linearized model
    integer :: yStar                  ! Index of specific vector to use for linearized model
    ! Now the logicals
    logical :: AllLinesForRadiometer  ! As opposed to just using lines designated for band.
    logical :: AllLinesInCatalog      ! Use all the lines
    logical :: Atmos_der              ! Do atmospheric derivatives
    logical :: Default_spectroscopy   ! Using Bill's spectroscopy data
    logical :: DifferentialScan       ! Differential scan model
    logical :: Do_1d                  ! Do 1D forward model calculation
    logical :: Do_baseline            ! Do a baseline computation
    logical :: Do_conv                ! Do convolution
    logical :: Do_freq_avg            ! Do Frequency averaging
    logical :: ForceFoldedOutput      ! Output to folded sideband even if signal is other (linear only)
    logical :: ForceSidebandFraction  ! If set mult. by SBfrac even if single sideband
    logical :: GlobalConfig           ! If set is shared between all chunks
    logical :: Incl_cld ! Include cloud extinction calculation in Bill's forward model
    logical :: LockBins               ! Use same l2pc bin for whole chunk
    logical :: Polarized              ! Use polarized model for Zeeman-split lines
    logical :: SkipOverlaps           ! Don't calculate for MAFs in overlap regions
    logical :: Spect_Der              ! Do spectroscopy derivatives
    logical :: SwitchingMirror        ! Model radiance at the switching mirror
    logical :: Temp_Der               ! Do temperature derivatives
    ! Now the reals
    real (r8) :: PhiWindow            ! Window size for examining stuff
    real (r8) :: Tolerance            ! Accuracy desired when choosing approximations
    real :: sum_DeltaTime = 0.0	      ! sum of delta time calling FullForwardModel 
    real :: sum_squareDeltaTime = 0.0 ! sum of the square of delta times calling FullForwardModel 
    ! Now the arrays
    integer, dimension(:), pointer :: BinSelectors=>NULL() ! List of relevant bin selectors
    integer, dimension(:), pointer :: Molecules=>NULL() ! Which molecules to consider
      ! >0 = beginning of a group or a lonesome molecule, <0 = member of a group.
      ! Size is one more than number in config; last one is a huge positive sentinel.
    logical, dimension(:), pointer :: MoleculeDerivatives=>NULL() ! Want Jacobians
    integer, dimension(:), pointer :: SpecificQuantities=>NULL() ! Specific quantities to use
    ! Now the derived types
    type (Signal_T), dimension(:), pointer :: Signals=>NULL()
    type (vGrid_T), pointer :: IntegrationGrid=>NULL() ! Zeta grid for integration
    type (vGrid_T), pointer :: TangentGrid=>NULL()     ! Zeta grid for integration
    ! Finally stuff that PVMPackFWMConfig and PVMUnpackFWMConfig don't cart around
    type (forwardModelDerived_T) :: ForwardModelDerived
  end type ForwardModelConfig_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  ! ----------------------------  AddForwardModelConfigToDatabase  -----
  integer function AddForwardModelConfigToDatabase ( Database, Item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_Deallocate, MLSMSG_Error

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database
    type (ForwardModelConfig_T), intent(in) :: Item

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddForwardModelConfigToDatabase = newSize
  end function AddForwardModelConfigToDatabase

  ! -------------------------------  DeriveFromForwardModelConfig  -----
  subroutine DeriveFromForwardModelConfig ( FwdModelConf )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use FilterShapes_m, only: DACSFilterShapes
    use MLSMessageModule, only: MLSMessage,  MLSMSG_Allocate, MLSMSG_Error
    use MLSSets, only: FindFirst, Intersection, Union
    use MLSSignals_M, only: MatchSignal
    use PFADataBase_m, only: PFAData
    use Toggles, only: Switches

    type (ForwardModelConfig_T), intent(inout) :: FwdModelConf

    integer :: Channel
    integer :: HowManyPFA                  ! How many PFA data records are for the
                                           ! molecules in fwdModelConf%Molecules?
    integer :: I, Ier, J
    integer :: LBoundDACs, UBoundDACs      ! How many channels in a DAC
    integer :: M1, M2                      ! Matched signal indices
    integer :: NoUsedChannels, NoUsedDACS
    integer :: NumPFA                      ! Like HowManyPFA, but for one channel
    integer :: SigInd
    logical :: SignalFlag(size(fwdModelConf%signals))
    integer :: S1, S2                      ! SidebandStart, SidebandStop
    integer, pointer :: T1(:), T2(:)       ! Temps for set operations
    integer, pointer :: WhichPFA(:)        ! Which PFA data records are for the
                                           ! molecules in fwdModelConf%Molecules?


    ! Shorthand pointers into fwdModelConf%forwardModelDerived
    real(rp), dimension(:,:), pointer :: DACsStaging  ! Temporary space for DACS radiances
    integer, pointer :: USEDDACSSIGNALS(:) ! Indices in FwdModelConf_T%Signals
                                           ! of signals for our dacs
    type(channels_T), pointer :: Channels(:)

    s1 = fwdModelConf%sidebandStart
    s2 = fwdModelConf%sidebandStop

    ! Identify which of our signals are DACS and how many unique DACS are involved
    ! Compute NoUsedDACs
    ! Allocate and compute UsedDACSSignals and allocate DACsStaging.

    nullify ( DACsStaging, usedDACSSignals, whichPFA )
    signalFlag = .false.
    lBoundDACs = 0; uBoundDACs = 0
    noUsedDACs = 0
    do sigInd = 1, size(fwdModelConf%signals)
      if ( fwdModelConf%signals(sigInd)%dacs .and. &
        & .not. signalFlag(sigind) ) then
        signalFlag(sigind) = .true.
        noUsedDACs = noUsedDACs + 1
        if ( noUsedDACs == 1 ) then
          lBoundDACs = lbound(fwdModelConf%signals(sigInd)%frequencies,1 )
          uBoundDACs = ubound(fwdModelConf%signals(sigInd)%frequencies,1 )
        else
          if ( lBoundDACs /= lbound ( fwdModelConf%signals(sigInd)%frequencies,1 ) .or. &
            &  uBoundDACs /= ubound ( fwdModelConf%signals(sigInd)%frequencies,1 ) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Two DACS have different number of channels' )
        end if
      end if
    end do
    if ( noUsedDACs > 0 .and. .not. associated(DACsFilterShapes) ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'DACS in use but no filter shapes provided.' )
    call allocate_test ( usedDACSSignals, noUsedDACs, 'usedDACSSignals', ModuleName )
    usedDACSSignals = pack ( (/ (i, i=1, size(signalFlag)) /), signalFlag )
    call allocate_test ( DACsStaging, uBoundDACs, noUsedDACs, &
      & 'DACsStaging', moduleName, low1 = lBoundDACs )

    ! Work out which channels are used.
    noUsedChannels = 0
    do sigInd = 1, size(fwdModelConf%signals)
      noUsedChannels = noUsedChannels + &
        & count( fwdModelConf%signals(sigInd)%channels )
    end do
    allocate ( channels(noUsedChannels), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'fwdModelConf%forwardModelDerived%channels' )

    ! Collect channel information from signals database.
    channel = 1
    do sigInd = 1, size(fwdModelConf%signals)
      do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
        if ( fwdModelConf%signals(sigInd)%channels(i) ) then
          channels(channel)%origin = &
            & lbound ( fwdModelConf%signals(sigInd)%frequencies, 1 )
          channels(channel)%used = i + channels(channel)%origin - 1
          channels(channel)%signal = sigInd
          channels(channel)%dacs = FindFirst ( usedDACSSignals, sigind )
          channel = channel + 1
        end if
      end do
    end do

    if ( associated(pfaData) ) then
      ! Work out which PFA data are germane to fwdModelConf%Molecules
      ! We do this because there are thousands of PFA data, but maybe only
      ! a few that are germane to this fwdModelConf
      howManyPFA = 0
      do i = 1, size(pfaData)
        do j = fwdModelConf%firstPFA, size(fwdModelConf%molecules) - 1
          howManyPFA = howManyPFA + &
            & count(fwdModelConf%molecules(j) == PFAData(i)%molecules)
        end do
      end do
      call allocate_test ( whichPFA, howManyPFA, 'whichPFA', moduleName )
      howManyPFA = 0
      do i = 1, size(pfaData)
        do j = fwdModelConf%firstPFA, size(fwdModelConf%molecules) - 1
          if ( any(fwdModelConf%molecules(j) == PFAData(i)%molecules) ) then
            howManyPFA = howManyPFA + 1
            whichPFA(howManyPFA) = i
          end if
        end do
      end do

      ! Work out PFA abstracts for each channel
      do i = 1, size(channels)
        numPFA = 0
        nullify ( t1, channels(i)%PFAMolecules )
        call allocate_test ( channels(i)%PFAMolecules, 0, &
          & 'channels(i)%PFAMolecules', moduleName ) ! so we can do unions
        do j = 1, howManyPFA
          if ( channels(i)%used < lbound(PFAData(whichPFA(j))%theSignal%channels,1) .or. &
               channels(i)%used > lbound(PFAData(whichPFA(j))%theSignal%channels,1) ) cycle
          if ( .not. PFAData(whichPFA(j))%theSignal%channels(channels(i)%used) ) cycle
          if ( .not. matchSignal ( &
            &  fwdModelConf%signals(channels(i)%signal:channels(i)%signal), &
            &  PFAData(whichPFA(j))%theSignal, DSBSSB=.true. ) /= 0) cycle
          numPFA = numPFA + 1
          t2 => intersection(fwdModelConf%molecules(fwdModelConf%firstPFA:), &
            &                PFAData(whichPFA(j))%molecules)
          t1 => union(channels(i)%PFAMolecules,t2)
          call deallocate_test ( t2, 't2', moduleName )
          call deallocate_test ( channels(i)%PFAMolecules, &
            & 'channels(i)%PFAMolecules', moduleName )
          channels(i)%PFAMolecules => t1
        end do
        allocate ( channels(i)%PFAData(s1:s2,numPFA), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'channels(i)%PFAData' )
        channels(i)%PFAData = 0
        numPFA = 0
        do j = 1, howManyPFA
          if ( channels(i)%used < lbound(PFAData(whichPFA(j))%theSignal%channels,1) .or. &
               channels(i)%used > lbound(PFAData(whichPFA(j))%theSignal%channels,1) ) cycle
          if ( .not. PFAData(whichPFA(j))%theSignal%channels(channels(i)%used) ) cycle
          m1 = matchSignal ( &
            &  fwdModelConf%signals(channels(i)%signal:channels(i)%signal), &
            &  PFAData(whichPFA(j))%theSignal, DSBSSB=.true., sideband=s1 )
          m2 = matchSignal ( &
            &  fwdModelConf%signals(channels(i)%signal:channels(i)%signal), &
            &  PFAData(whichPFA(j))%theSignal, DSBSSB=.true., sideband=s2 )
          if ( m1 + m2 /= 0 ) then
            numPFA = numPFA + 1
            if ( m1 /= 0 ) channels(i)%PFAData(s1,numPFA) = whichPFA(j)
            if ( m2 /= 0 ) channels(i)%PFAData(s2,numPFA) = whichPFA(j)
          end if
        end do
      end do
      call deallocate_test ( whichPFA, 'HowManyPfa', moduleName )
    end if ! associated(pfaData)

    ! Hook the shortcuts into the structure
    fwdModelConf%forwardModelDerived%channels => channels
    fwdModelConf%forwardModelDerived%DACsStaging => DACsStaging
    fwdModelConf%forwardModelDerived%usedDACSSignals => usedDACSSignals

    if ( index(switches,'fwmd') /= 0 ) call dump ( fwdModelConf )

  end subroutine DeriveFromForwardModelConfig

  ! --------------------------  DestroyForwardModelConfigDatabase  -----
  subroutine DestroyFWMConfigDatabase ( Database, Deep )

    use MLSMessageModule, only: MLSMessage,  MLSMSG_Deallocate, MLSMSG_Error

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database
    logical, optional, intent(in) :: DEEP

    ! Local variables
    integer :: Config                   ! Loop counter
    integer :: Status                   ! Flag

    if ( associated(database) ) then
      do config = 1, size(database)
        call DestroyOneForwardModelConfig ( database(config), deep=deep )
      end do

      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "Database" )
    end if
  end subroutine DestroyFWMConfigDatabase

  ! ---------------------------------  DestroyForwardModelDerived  -----
  subroutine DestroyForwardModelDerived ( FwdModelConf )
    ! Destroy FwdModelConf%ForwardModelDerived

    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type ( ForwardModelConfig_T ), intent(inout) :: FwdModelConf

    integer :: I, Ier

    if ( associated(fwdModelConf%forwardModelDerived%channels) ) then
      do i = 1, size(fwdModelConf%forwardModelDerived%channels)
        call deallocate_test ( &
          & fwdModelConf%forwardModelDerived%channels(i)%PFAData, &
          & 'fwdModelConf%forwardModelDerived%channels(i)%PFAData', moduleName )
        call deallocate_test ( &
          & fwdModelConf%forwardModelDerived%channels(i)%PFAMolecules, &
          & 'fwdModelConf%forwardModelDerived%channels(i)%PFAMolecules', moduleName )
      end do
      deallocate ( fwdModelConf%forwardModelDerived%channels, stat = ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'fwdModelConf%forwardModelDerived%channels' )
    end if

    call deallocate_test ( fwdModelConf%forwardModelDerived%DACSStaging, &
      & 'fwdModelConf%forwardModelDerived%DACSStaging', moduleName )

    call deallocate_test ( fwdModelConf%forwardModelDerived%usedDACSSignals, &
      & 'fwdModelConf%forwardModelDerived%usedDACSSignals', moduleName )

  end subroutine DestroyForwardModelDerived
  ! ------------------------------------ NullifyForwardModelConfig -----
  subroutine NullifyForwardModelConfig ( F )
    ! Given a forward model config, nullify all the pointers associated with it
    type ( ForwardModelConfig_T ), intent(out) :: F

    ! Executable code not needed since => NULL() initializes pointer
    ! components of intent(out) dummy arguments.
  end subroutine NullifyForwardModelConfig

  ! ------------------------------------------- PVMPackFwmConfig --------
  subroutine PVMPackFWMConfig ( config )
    use PVMIDL, only: PVMIDLPack
    use PVM, only: PVMErrorMessage
    use MorePVM, only: PVMPackLitIndex, PVMPackStringIndex
    use MLSSignals_m, only: PVMPackSignal
    use VGridsDatabase, only: PVMPackVGrid
    ! Dummy arguments
    type ( ForwardModelConfig_T ), intent(in) :: CONFIG
    ! Local variables
    integer :: INFO                     ! Flag from PVM
    integer :: I                        ! Loop counter

    ! Executable code
    ! First pack the lit indices
    call PVMPackStringIndex ( config%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig name" )
    call PVMPackLitIndex ( config%cloud_der, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig cloud_der" )
    call PVMPackLitIndex ( config%fwmType, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig fwmType" )
    call PVMPackLitIndex ( config%i_saturation, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig i_saturation" )
    call PVMPackLitIndex ( config%instrumentModule, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig instrumentModule" )
    call PVMPackLitIndex ( config%windowUnits, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig windowUnits" )

    ! Now pack the integer scalars
    call PVMIDLPack ( (/ &
      & config%linearSideband, &
      & config%no_cloud_species, config%no_model_surfs, &
      & config%num_ab_terms, config%num_azimuth_angles, &
      & config%num_scattering_angles, config%num_size_bins, &
      & config%sidebandStart, config%sidebandStop, &
      & config%surfaceTangentIndex /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig integers" )

    ! Now the logical scalars
    call PVMIDLPack ( (/ &
      & config%allLinesForRadiometer, config%atmos_der, &
      & config%default_spectroscopy, config%differentialScan,&
      & config%do_1d, config%do_baseline, config%do_conv, &
      & config%do_freq_avg, config%forceFoldedOutput, config%forceSidebandFraction, &
      & config%globalConfig, config%incl_cld, &
      & config%lockBins, config%polarized, config%skipOverlaps, &
      & config%spect_Der, config%temp_Der /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig logicals" )

    ! Now pack the reals
    call PVMIDLPack ( (/ config%phiWindow, config%tolerance /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing fwmConfig reals" )

    ! ------------- The rest are arrays and/or types
    ! Bin selectors
    if ( associated ( config%binSelectors ) ) then
      call PVMIDLPack ( size ( config%binSelectors ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing number of binSelectors" )
      call PVMIDLPack ( config%binSelectors, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing binSelectors" )
    else
      call PVMIDLPack ( 0, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing 0 binSelectors" )
    end if

    ! Molecules / derivatives
    if ( associated ( config%molecules ) ) then
      call PVMIDLPack ( size ( config%molecules ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing number of molecules" )
      call PVMIDLPack ( config%firstPFA, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing index of first PFA molecule" )
      if ( size ( config%molecules ) > 0 ) then
        do i = 1, size(config%molecules) - 1
          call PVMPackLitIndex ( abs ( config%molecules(i) ), info )
          if ( info /= 0 ) call PVMErrorMessage ( info, "Packing a molecule" )
          call PVMIDLPack ( (config%molecules(i) .gt. 0.0), info )
          if ( info /= 0 ) call PVMErrorMessage ( info, "Packing molecule sign" )
        end do
        call PVMIDLPack ( config%moleculeDerivatives, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "Packing molecule derivatives" )
      end if
    else
      call PVMIDLPack ( 0, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing 0 molecules" )
    end if

    ! Specific quantities
    if ( associated ( config%specificQuantities ) ) then
      call PVMIDLPack ( size ( config%specificQuantities ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing number of specificQuantities" )
      call PVMIDLPack ( config%specificQuantities, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing specificQuantities" )
    else
      call PVMIDLPack ( 0, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing 0 specificQuantities" )
    end if

    ! Pack the other structures - signals
    if ( associated ( config%signals ) ) then
      call PVMIDLPack ( size ( config%signals ), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing number of signals" )
      do i = 1, size ( config%signals )
        call PVMPackSignal ( config%signals(i) )
      end do
    else
      call PVMIDLPack ( 0, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Packing 0 signals" )
    end if

    ! Vgrids
    call PVMIDLPack ( (/ associated ( config%integrationGrid ), &
      & associated ( config%tangentGrid ) /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Packing vGrid flags" )
    if ( associated ( config%integrationGrid ) ) &
      & call PVMPackVGrid ( config%integrationGrid )
    if ( associated ( config%tangentGrid ) ) &
      & call PVMPackVGrid ( config%tangentGrid )

  end subroutine PVMPackFWMConfig

  ! ----------------------------------------- PVMUnpackFWMConfig ---------
  subroutine PVMUnpackFWMConfig ( CONFIG )
    use PVMIDL, only: PVMIDLUnpack
    use PVM, only: PVMErrorMessage
    use MorePVM, only: PVMUnpackLitIndex, PVMUnpackStringIndex
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
    use Allocate_Deallocate, only: Allocate_test
    use MLSSignals_m, only: PVMUnpackSignal
    use VGridsDatabase, only: PVMUnpackVGrid
    ! Dummy arguments
    type ( ForwardModelConfig_T ), intent(out) :: CONFIG
    ! Local variables
    integer :: INFO                     ! Flag from PVM
    logical :: FLAG                     ! A flag from the sender
    logical, dimension(17) :: LS        ! Temporary array
    integer, dimension(10) :: IS         ! Temporary array
    real(r8), dimension(2) :: RS        ! Temporary array
    integer :: I                        ! Loop counter
    integer :: N                        ! Array size

    ! Executable code
    ! First unpack the lit indices
    call PVMUnpackStringIndex ( config%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig name" )
    call PVMUnpackLitIndex ( config%cloud_der, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig cloud_der" )
    call PVMUnpackLitIndex ( config%fwmType, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig fwmType" )
    call PVMUnpackLitIndex ( config%i_saturation, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig i_saturation" )
    call PVMUnpackLitIndex ( config%instrumentModule, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig instrumentModule" )
    call PVMUnpackLitIndex ( config%windowUnits, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig windowUnits" )

    ! Now the integer scalars
    call PVMIDLUnpack ( is, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig integers" )
    i = 1
    config%linearsideband        = is(i) ; i = i + 1
    config%no_cloud_species      = is(i) ; i = i + 1
    config%no_model_surfs        = is(i) ; i = i + 1
    config%num_ab_terms          = is(i) ; i = i + 1
    config%num_azimuth_angles    = is(i) ; i = i + 1
    config%num_scattering_angles = is(i) ; i = i + 1
    config%num_size_bins         = is(i) ; i = i + 1
    config%sideBandStart         = is(i) ; i = i + 1
    config%sideBandStop          = is(i) ; i = i + 1
    config%surfaceTangentIndex   = is(i) ; i = i + 1

    ! Now the logical scalars. Array LS has to be long enough for this.
    ! If you add any items here, don't forget to make LS longer at the top
    ! of the program unit.
    call PVMIDLUnpack ( ls, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig logicals" )
    i = 1
    config%allLinesForRadiometer = ls(i) ; i = i + 1
    config%atmos_der             = ls(i) ; i = i + 1
    config%default_spectroscopy  = ls(i) ; i = i + 1
    config%differentialScan      = ls(i) ; i = i + 1
    config%do_1d                 = ls(i) ; i = i + 1
    config%do_baseline           = ls(i) ; i = i + 1
    config%do_conv               = ls(i) ; i = i + 1
    config%do_freq_avg           = ls(i) ; i = i + 1
    config%forceFoldedOutput     = ls(i) ; i = i + 1
    config%forceSidebandFraction = ls(i) ; i = i + 1
    config%globalConfig          = ls(i) ; i = i + 1
    config%incl_cld              = ls(i) ; i = i + 1
    config%lockBins              = ls(i) ; i = i + 1
    config%polarized             = ls(i) ; i = i + 1
    config%skipOverlaps          = ls(i) ; i = i + 1
    config%spect_der             = ls(i) ; i = i + 1
    config%temp_der              = ls(i) ; i = i + 1

    ! Now the real scalars
    call PVMIDLUnpack ( rs, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking fwmConfig reals" )
    i = 1
    config%phiWindow = rs(i) ; i = i + 1
    config%tolerance = rs(i) ; i = i + 1

    ! ------- The rest are arrays and/or types
    ! Bin selectors
    call PVMIDLUnpack ( n, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking number of specific quantities" )
    if ( n > 0 ) then
      call Allocate_test ( config%binSelectors, n, &
        & 'config%binSelectors', ModuleName )
      call PVMIDLUnpack ( config%binSelectors, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking binSelectors" )
    end if

    ! Molecules / derivatives
    call PVMIDLUnpack ( n, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking number of molecules" )
    call PVMIDLUnpack ( config%firstPFA, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking index of PFA molecule" )
    if ( n > 0 ) then
      call Allocate_test ( config%molecules, n, 'config%molecules', ModuleName )
      call Allocate_test ( config%moleculeDerivatives, &
        & n, 'config%moleculeDerivatives', ModuleName )
      do i = 1, n - 1
        call PVMUnpackLitIndex ( config%molecules(i), info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking a molecule" )
        call PVMIDLUnpack ( flag, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking a molecule sign flag" )
        if ( .not. flag ) config%molecules(i) = - config%molecules(i)
      end do
      config%molecules(n) = huge(config%molecules(n)) ! Sentinel
      call PVMIDLUnpack ( config%moleculeDerivatives, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking moleculeDerivatives" )
    end if

    ! Specific quantities
    call PVMIDLUnpack ( n, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking number of specific quantities" )
    if ( n > 0 ) then
      call Allocate_test ( config%specificQuantities, n, &
        & 'config%specificQuantities', ModuleName )
      call PVMIDLUnpack ( config%specificQuantities, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking specific quantities" )
    end if

    ! Unpack other structures - signals
    call PVMIDLUnpack ( n, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking number of signals" )
    if ( n > 0 ) then
      allocate ( config%signals(n), STAT=info )
      if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'config%signals' )
      do i = 1, n
        call PVMUnpackSignal ( config%signals(i) )
      end do
    end if

    ! Vgrids
    call PVMIDLUnpack ( ls(1:2), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "Unpacking vGrid flags" )
    if ( ls(1) ) then
      allocate ( config%integrationGrid, STAT=info )
      if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'config%integrationGrid' )
      call PVMUnpackVGrid ( config%integrationGrid )
    end if
    if ( ls(2) ) then
      allocate ( config%tangentGrid, STAT=info )
      if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'config%tangentGrid' )
      call PVMUnpackVGrid ( config%tangentGrid )
    end if

  end subroutine PVMUnpackFWMConfig

  ! --------------------------  StripForwardModelConfigDatabase --------
  subroutine StripForwardModelConfigDatabase ( database )
    ! This routine removes the non-global forward model configs from the database
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_Deallocate, MLSMSG_Error

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: DATABASE

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: TMPDATABASE
    integer :: CONFIG                   ! Loop counter
    integer :: STATUS

    ! Executable code
    ! Clear out dying configs
    if ( .not. associated ( database ) ) return
    do config = 1, size ( database )
      if ( .not. database(config)%globalConfig ) &
        & call DestroyOneForwardModelConfig ( database(config) )
    end do

    ! Create new database in tmp space, pack old one into
    allocate ( tmpDatabase ( count ( database%globalConfig ) ), STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // 'tmpDatabase' )
    tmpDatabase = pack ( database, database%globalConfig )

    ! Destroy old database, then point to new one
    deallocate ( database, STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate // 'database' )
    database => tmpDatabase
  end subroutine StripForwardModelConfigDatabase

  ! =====     Private Procedures     =====================================

  ! ------------------------------------ DestroyOneForwardModelConfig --
  subroutine DestroyOneForwardModelConfig ( Config, Deep )
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSSignals_M, only: DestroySignalDatabase
    use MLSMessageModule, only: MLSMSG_Deallocate, MLSMSG_Error, MLSMessage

    ! Dummy arguments
    type ( ForwardModelConfig_T), intent(inout) :: config
    logical, optional, intent(in) :: DEEP ! Do a really deep destroy

    ! Local variables
    integer :: STATUS                   ! Flag from allocate etc.
    logical :: MYDEEP                   ! Copy of deep

    ! Executable code
    myDeep = .false.
    if ( present ( deep ) ) myDeep = deep
    if ( associated(config%signals) ) &
      & call DestroySignalDatabase ( config%signals, justChannels=.not. myDeep )
    if ( myDeep ) then
      if ( associated ( config%integrationGrid ) ) then
        call DestroyVGridContents ( config%integrationGrid )
        deallocate ( config%integrationGrid, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate // "config%integrationGrid" )
      end if
      if ( associated ( config%tangentGrid ) ) then
        call DestroyVGridContents ( config%tangentGrid )
        deallocate ( config%tangentGrid, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate // "config%tangentGrid" )
      end if
    end if
    ! Otherwise don't destroy integrationGrid and tangentGrid.  Assume they will
    ! be (or already are) destroyed by destroyVGridDatabase.
    call deallocate_test ( config%molecules, &
      & "config%molecules", moduleName )
    call deallocate_test ( config%moleculeDerivatives, &
      & "config%moleculeDerivatives", moduleName )
    call Deallocate_test ( config%specificQuantities, &
      & "config%specificQuantities", ModuleName )
    call Deallocate_test ( config%binSelectors, &
      & "config%binSelectors", ModuleName )
    call destroyForwardModelDerived ( config )
  end subroutine DestroyOneForwardModelConfig

  ! ----------------------------  Dump_ForwardModelConfigDatabase  -----
  subroutine Dump_ForwardModelConfigDatabase ( Database, Where )

    use MoreTree, only: StartErrorMessage
    use Output_M, only: Output

    type (ForwardModelConfig_T), pointer, dimension(:) :: Database
    integer, optional, intent(in) :: Where ! Tree node index

    ! Local variables
    integer :: I                         ! Loop counters

    ! executable code
    if ( associated(database) ) then
      do i = 1, size(database)
        call Dump_ForwardModelConfig( database(i) )
      end do
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No forward model database to dump.', advance='yes' )
    end if
  end subroutine Dump_ForwardModelConfigDatabase

  ! -----------------------------------  Dump_ForwardModelConfig  -----
  subroutine Dump_ForwardModelConfig ( Config )

    use Dump_0, only: DUMP
    use Intrinsic, only: Lit_indices
    use MLSSignals_M, only: GetNameOfSignal, MaxSigLen
    use Output_M, only: NewLine, Output
    use String_Table, only: Display_String

    type (ForwardModelConfig_T):: Config

    ! Local variables
    integer ::  I, J                         ! Loop counters
    character (len=MaxSigLen) :: SignalName  ! A line of text

    ! executable code

    call output ( '  Forward Model Config Name: ' )
    call display_string ( Config%name, advance='yes' )
    call output ( '  Atmos_der:' )
    call output ( Config%atmos_der, advance='yes' )
    call output ( '  Do_conv:' )
    call output ( Config%do_conv, advance='yes' )
    call output ( '  Do_Baseline:' )
    call output ( Config%do_Baseline, advance='yes' )
    call output ( '  Default_spectroscopy:' )
    call output ( Config%Default_spectroscopy, advance='yes' )
    call output ( '  Do_freq_avg:' )
    call output ( Config%do_freq_avg, advance='yes' )
    call output ( '  Do_1D:' )
    call output ( Config%do_1d, advance='yes' )
    call output ( '  Incl_Cld:' )
    call output ( Config%incl_cld, advance='yes' )
    call output ( '  SkipOverlaps:' )
    call output ( Config%skipOverlaps, advance='yes' )
    call output ( '  Spect_der:' )
    call output ( Config%spect_der, advance='yes' )
    call output ( '  Temp_der:' )
    call output ( Config%temp_der, advance='yes' )
    call output ( '  Molecules: ', advance='yes' )
    if ( associated(Config%molecules) ) then
      do j = 1, size(Config%molecules) - 1
        call output ( '    ' )
        if ( j == config%firstPFA ) call output ( 'PFA: ' )
        if ( Config%molecules(j) < 0 ) call output ( '-' )
        call display_string(lit_indices(abs(Config%molecules(j))))
        if (Config%moleculeDerivatives(j)) then
          call output (' compute derivatives', advance='yes')
        else
          call output (' no derivatives', advance='yes')
        end if
      end do
    end if
    call output ( '  Signals:', advance='yes')
    do j = 1, size(Config%signals)
      call getNameOfSignal ( Config%signals(j), signalName)
      call output ( '  ' )
      call output ( trim(signalName)//' channelIncluded:', advance='yes')
      call dump ( Config%signals(j)%channels )
    end do
    ! Dump ForwardModelDerived
    call output ( '  ForwardModelDerived:', advance='yes' )
    if ( associated(Config%forwardModelDerived%usedDACSSignals) ) &
      call dump  ( Config%forwardModelDerived%usedDACSSignals, &
        & name='   Used DACS Signals' )
    if ( associated(Config%forwardModelDerived%channels) ) then
      call output ( '   Channel info:', advance='yes' )
      do j = 1, size(Config%forwardModelDerived%channels)
        call output ( Config%forwardModelDerived%channels(j)%used, before='    Used: ' )
        call output ( Config%forwardModelDerived%channels(j)%origin, before='    Origin: ' )
        call output ( Config%forwardModelDerived%channels(j)%signal, before='    Signal: ' )
        call output ( Config%forwardModelDerived%channels(j)%DACS, before='    DACS: ' )
        if ( associated(Config%forwardModelDerived%channels(j)%PFAData) ) &
          & call dump ( Config%forwardModelDerived%channels(j)%PFAData, &
            & name='    PFAData' )
        if ( associated(Config%forwardModelDerived%channels(j)%PFAMolecules) ) then
          call output ( '    PFA Molecules:', advance='no' )
          if ( size(Config%forwardModelDerived%channels(j)%PFAMolecules) == 0 ) &
            call output ( ' Empty' )
          do i = 1, size(Config%forwardModelDerived%channels(j)%PFAMolecules)
            call output ( ' ', advance='no' )
            call display_string ( &
              & lit_indices(Config%forwardModelDerived%channels(j)%PFAMolecules(i)), &
              & advance='no' )
          end do
          call newLine
        end if
      end do
    end if
  end subroutine Dump_ForwardModelConfig

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ForwardModelConfig

! $Log$
! Revision 2.55  2004/08/03 21:40:06  vsnyder
! Inching further toward PFA
!
! Revision 2.54  2004/07/16 19:13:18  vsnyder
! Correct problem in dump routine
!
! Revision 2.53  2004/07/08 02:35:29  vsnyder
! Put all line-by-line molecules before PFA ones
!
! Revision 2.52  2004/06/23 02:14:06  vsnyder
! Added PFA stuff, some cannonball polishing
!
! Revision 2.51  2004/06/11 01:33:29  vsnyder
! Declare all pointer components to be initially NULL
!
! Revision 2.50  2004/06/10 00:59:55  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.49  2004/05/26 23:54:14  vsnyder
! Don't dump the database if it's not allocated
!
! Revision 2.48  2004/05/01 04:00:59  vsnyder
! Added pfaMolecules field
!
! Revision 2.47  2004/03/24 13:50:56  hcp
! Made array LS 17 elements instead of 16. 17 things are taken out of it,
! so presumably it should be that long. NAG f95 flagged this as an error.
!
! Revision 2.46  2004/03/22 18:23:20  livesey
! Added AllLinesInCatalog flag
!
! Revision 2.45  2003/10/28 23:43:47  livesey
! Added forceFoldedOutput
!
! Revision 2.44  2003/10/18 01:15:58  livesey
! Various changes to the pack/unpack stuff.  This currently needs more
! attention.
!
! Revision 2.43  2003/09/15 23:45:00  vsnyder
! Remove unused local variables and USEs
!
! Revision 2.42  2003/09/11 23:10:04  livesey
! Added xStar and yStar
!
! Revision 2.41  2003/07/22 22:43:39  michael
! Added a Dump_ForwardModelConfig subroutine for a single configuration and made
! Dump_ForwardModelConfigDatabase call Dump_ForwardModelConfig to dump an array
! of configurations.  Either is called with the generic interface "dump".
! Formerly, only pointers to arrays of forwardmodelconsigurations could be dumped.
!
! Revision 2.40  2003/07/15 22:09:59  livesey
! Added forceSidebandFraction and linearSideband
!
! Revision 2.39  2003/07/15 18:16:26  livesey
! Added name to configuration
!
! Revision 2.38  2003/06/30 22:55:01  cvuu
! Find mean, std dev timing of fullForwardModel calls
!
! Revision 2.37  2003/06/18 01:58:01  vsnyder
! Add SidebandStart and SidebandStop fields
!
! Revision 2.36  2003/05/29 16:37:02  livesey
! Added switchingMirror
!
! Revision 2.35  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.34  2003/05/05 23:00:24  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.33  2003/04/04 00:26:44  jonathan
! change dimension(12) :: i11 TO dimension(11) :: i11
!
! Revision 2.32  2003/04/02 21:46:51  jonathan
! remove cloud_fov, changed i12 to i11
!
! Revision 2.31  2003/03/07 03:16:00  livesey
! Changed use of DestroySignal
!
! Revision 2.30.2.2  2003/04/08 23:40:11  jonathan
! remove cloud_fov
!
! Revision 2.30.2.1  2003/02/22 00:48:08  vsnyder
! Delete moleculesPol, moleculeDerivativesPol, add Polarized to ForwardModelConfig
!
! Revision 2.30  2003/02/06 22:04:25  vsnyder
! Add f_moleculesPol, f_moleculeDerivativesPol, delete f_polarized
!
! Revision 2.29  2003/02/06 20:16:23  livesey
! Minor bug fix in the database deallocation
!
! Revision 2.28  2003/02/05 21:57:27  livesey
! Tidy up and added binSelectors, removed nameFragment
!
! Revision 2.27  2003/01/30 22:01:30  livesey
! Tidy up of the logical array packing.
!
! Revision 2.26  2003/01/30 18:29:40  jonathan
! change dimension l13 to 14
!
! Revision 2.25  2003/01/30 17:28:01  jonathan
! add logical incl_cld
!
! Revision 2.24  2003/01/29 01:48:52  vsnyder
! Add 'polarized' field to forwardModel
!
! Revision 2.23  2003/01/26 04:42:42  livesey
! Added profiles/angle options for phiWindow
!
! Revision 2.22  2003/01/17 00:01:44  livesey
! Another bug fix in the packing/unpacking
!
! Revision 2.21  2003/01/16 05:53:19  livesey
! Bug fix to Jonathans new configs
!
! Revision 2.20  2003/01/16 05:50:59  livesey
! Bug fix in do_1d handling
!
! Revision 2.19  2003/01/16 00:55:27  jonathan
! add do_1d, also fix bug of reversed  do_freq_avg do_baseline order
!
! Revision 2.18  2003/01/13 17:16:23  jonathan
! chane cloud_width to i_saturation
!
! Revision 2.17  2002/12/04 21:55:22  livesey
! Added the name fragment packing
!
! Revision 2.16  2002/11/22 12:14:31  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.15  2002/11/15 01:32:53  livesey
! Added allLinesForRadiometer
!
! Revision 2.14  2002/10/08 17:40:01  livesey
! Various bug fixes in the pack/unpack routines
!
! Revision 2.13  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/10/06 01:09:22  livesey
! Made the pvm routines public
!
! Revision 2.11  2002/10/05 00:41:01  livesey
! Added pvm pack and unpack config and deep option on destroy
!
! Revision 2.10  2002/09/25 20:06:42  livesey
! Added specificQuantities, which necessitated globalConfig to allow for
! some configs inside construct.  This in turn required
! StripForwardModelConfigDatabase
!
! Revision 2.9  2002/08/21 23:53:57  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.8  2002/07/17 06:02:13  livesey
! New config elements for hdf5 l2pcs
!
! Revision 2.7  2002/06/12 17:00:49  livesey
! Changed phiWindow to float
!
! Revision 2.6  2002/03/07 17:17:49  livesey
! Removed frqGap
!
! Revision 2.5  2002/02/14 23:01:35  livesey
! Added justChannels in call to destroySignal
!
! Revision 2.4  2002/02/13 00:09:24  livesey
! Added differential Scan model
!
! Revision 2.3  2001/11/15 23:50:11  jonathan
! rename DF_spectroscopy to default_spectroscopy
!
! Revision 2.2  2001/11/15 20:56:38  jonathan
! add df_spectroscopy
!
! Revision 2.1  2001/10/02 20:37:09  livesey
! Added do_baseline
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.15  2001/09/04 15:59:01  jonathan
! add cloud_fov, jonathan
!
! Revision 1.14  2001/07/17 22:38:05  jonathan
! add cloud_width, jonathan/paul
!
! Revision 1.13  2001/07/16 22:07:57  jonathan
! change cloud_der to integer-type, jonathan
!
! Revision 1.12  2001/07/06 18:55:13  jonathan
! Modified for cloud model, Paul/Jonathan
!
! Revision 1.11  2001/06/21 15:05:53  livesey
! Added tolerance
!
! Revision 1.10  2001/05/31 23:07:45  livesey
! Added cloud_der
!
! Revision 1.9  2001/05/25 20:26:09  livesey
! Added skipOverlaps option
!
! Revision 1.8  2001/05/14 23:17:35  livesey
! Added frqGap parameter
!
! Revision 1.7  2001/05/03 23:07:02  livesey
! Added scan model stuff
!
! Revision 1.6  2001/05/02 20:30:36  livesey
! Removed frequency from config
!
! Revision 1.5  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.4  2001/04/21 01:08:57  vsnyder
! Deallocate Molecules and MoleculeDerivatives in DestroyFWMConfigDatabase
!
! Revision 1.3  2001/04/12 17:00:08  vsnyder
! Comment out a line with an undefined variable on it
!
! Revision 1.2  2001/04/10 22:17:05  livesey
! Renamed module
!
! Revision 1.1  2001/04/07 01:56:25  vsnyder
! Initial commit
!
