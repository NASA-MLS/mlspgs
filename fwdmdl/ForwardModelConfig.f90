! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelConfig
!=============================================================================

! Set up the forward model configuration, except for actually processing
! the command.

  use MLSCommon, only: R8, RP
  use MLSSignals_M, only: Signal_T
  use SpectroscopyCatalog_m, only: Catalog_t
  use VectorsModule, only: VectorValue_T
  use VGridsDatabase, only: VGrid_T, DestroyVGridContents

  implicit NONE
  private

  ! Public procedures:
  interface Dump
    module procedure Dump_Beta_Group, Dump_ForwardModelConfig
    module procedure Dump_ForwardModelConfigDatabase, Dump_Qty_Stuff
  end interface Dump

  public :: AddForwardModelConfigToDatabase, DeriveFromForwardModelConfig
  public :: DestroyFWMConfigDatabase, DestroyForwardModelDerived
  public :: Dump, NullifyForwardModelConfig
  public :: StripForwardModelConfigDatabase, PVMPackFWMConfig, PVMUnpackFWMConfig
 
  ! Public Types:

  ! Quantities derived from forward models, but not carted around by
  ! PVMPackFWMConfig and PVMUnpackFWMConfig.  Rather, they are computed by
  ! DeriveFromForwardModelConfig, either when a ForwardModelConfig_T is
  ! created or when one arrives by way of PVMUnpackFWMConfig.

  type, public :: QtyStuff_T ! So we can have an array of pointers to QTY's
    type (VectorValue_T), pointer :: QTY => NULL()
    logical :: FoundInFirst
  end type QtyStuff_T

  ! Beta group type declaration.  Each entry in the Molecules list of the form
  ! "m" has one of these with n_elements == 1, referring to "m".  Each entry
  ! of the form "[m,m1,...,mn]" has one of these with n_elements == n, referring
  ! to m1,...mn (but not m).  The Cat_Index, Qty and Ratio components are
  ! filled in Get_Species_Data.
  type, public :: Beta_Group_T
    ! For the group as a whole:
    logical :: Derivatives = .false.  ! "Compute derivatives w.r.t. mixing ratio"
    logical :: Group = .false.        ! "Molecule group", i.e., [m,m1,...,mn]
    integer :: Molecule               ! Group name, i.e., "m".
    type(qtyStuff_t) :: Qty           ! The Qty's vector and foundInFirst
    ! For LBL molecules in the group:
    integer, pointer  :: Cat_Index(:) => NULL() ! 1:size(LBL_Molecules).  Indices
                                      ! for config%catalog and gl_slabs.
                                      ! Allocated and filled in Get_Species_Data.
    integer, pointer :: LBL_Molecules(:) => NULL() ! LBL molecules in the group
                                      ! if a group, i.e., "m1...mn", else "m" if
                                      ! "m" is LBL, else zero size.
    real(rp), pointer :: LBL_Ratio(:) => NULL() ! 1:size(LBL_Molecules).  Isotope
                                      ! ratio.  Allocated in ForwardModelSupport
                                      ! with value 1.0, but could be filled in
                                      ! Get_Species_Data.
    ! For PFA molecules in the group:
    integer, pointer :: PFA_Indices(:,:,:) => NULL() ! Sidebands x Channels
                                      ! x 1:size(PFA_Molecules).  Indices in
                                      ! PFADataBase%PFAData.  Allocated and
                                      ! filled in Get_Species_Data.
    integer, pointer :: PFA_Molecules(:) => NULL() ! PFA molecules in the group
                                      ! if a group, i.e., "m1...mn", else "m" if
                                      ! "m" is PFA, else zero size.
    real(rp), pointer :: PFA_Ratio(:) => NULL() ! 1:size(PFA_Molecules).  Isotope
                                      ! ratio.  Allocated in ForwardModelSupport
                                      ! with value 1.0, but could be filled in
                                      ! Get_Species_Data.
  end type Beta_Group_T

  ! Channel information from the signals database
  type, public :: Channels_T
    integer :: Used       ! Which channel is this?
    integer :: Origin     ! Index of first channel (zero or one)
    integer :: Signal     ! Signal index for the channel
    integer :: DACS       ! DACS index if any, else zero
    integer, pointer :: PFAIndex(:,:) => NULL()   ! Sidebands X numPFA
                          ! Indices in PFADataBase%PFAData
    integer, pointer :: PFAMolecules(:) => NULL() ! L_... from Molecules(firstPFA:)
                          ! for this channel
    integer, pointer :: BetaIndex(:) => NULL()    ! Allocated here but filled
                          ! later in Get_Species_Data; BetaIndex(i) is second
                          ! subscript for beta_path for PFAMolecules(i).
  end type Channels_T

  ! The scalar components are sorted in the order they are to make the packing
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
    logical :: AnyPFA                 ! Set if there are any PFA molecules
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
    type (beta_group_t), dimension(:), pointer :: Beta_Group => NULL() ! q.v. above
      ! Some of this is filled in each time the forward model is invoked.  Those
      ! parts aren't dragged around by PVM.
    type (Signal_T), dimension(:), pointer :: Signals=>NULL()
    type (vGrid_T), pointer :: IntegrationGrid=>NULL() ! Zeta grid for integration
    type (vGrid_T), pointer :: TangentGrid=>NULL()     ! Zeta grid for integration
    ! Finally stuff that PVMPackFWMConfig and PVMUnpackFWMConfig don't cart around
    ! This stuff is filled in each time the forward model is invoked, partly by
    ! DeriveFromForwardModel and partly by Get_Species_Data.
    type (catalog_t), pointer :: Catalog(:,:) => NULL()     ! sidebands,1:noNonPFA
      ! Catalog's second subscript comes from Beta_Group%Cat_Index for each
      ! element of the beta group.  We do this indirectly so that the whole
      ! catalog can be turned into gl_slabs data at once, and then the gl_slabs
      ! can be indexed by Beta_Group%Cat_Index as well.
    real(rp), dimension(:,:), pointer :: DACsStaging => NULL() ! Temporary
      ! space for DACS radiances
    integer, dimension(:), pointer :: USEDDACSSIGNALS => NULL() ! Indices in
      ! FwdModelConf_T%Signals of signals for our dacs
    type(channels_T), pointer, dimension(:) :: Channels => NULL()
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
    use MLSMessageModule, only: MLSMessage,  MLSMSG_Allocate, MLSMSG_Error, &
      & MLSMSG_Warning
    use MLSSets, only: FindFirst
    use MLSSignals_m, only: MatchSignal
!   use Output_m, only: Output
    use PFADataBase_m, only: Dump, PFAData, PFA_By_Molecule, SortPFAData, &
      & Sort_PFADataBase
    use Toggles, only: Switches

    type (ForwardModelConfig_T), intent(inout) :: FwdModelConf

    integer :: DumpFwm = -1                ! -1 = not called yet, 0 = no dumps,
                                           ! 1 = dump, 2 = dump and stop
    logical :: Error
    integer :: S1, S2                      ! SidebandStart, SidebandStop

    if ( dumpFwm < 0 ) then ! done only once
      dumpFwm = 0
      if ( index(switches,'fwmd') /= 0 )  dumpFwm = 1
      if ( index(switches,'fwmD') /= 0 )  dumpFwm = 2
    end if

    call sort_PFADatabase ! Only does anything once

    error = .false.

    s1 = fwdModelConf%sidebandStart
    s2 = fwdModelConf%sidebandStop

    ! Identify which of our signals are DACS and how many unique DACS are involved
    ! Allocate and compute UsedDACSSignals and allocate DACsStaging.
    call DACS_Stuff ( fwdModelConf%DACsStaging, &
                    & fwdModelConf%usedDACSSignals ) ! Below

    ! Work out which channels are used.
    call channel_stuff ( fwdModelConf%channels, fwdModelConf%usedDACSSignals ) ! Below

    ! Work out the spectroscopy we're going to need.
    call SpectroscopyCatalogExtract ! Below

    ! Work out the PFA stuff
    call PFA_Stuff ! Below

    if ( dumpFwm > 0 .or. error ) then
      call dump ( fwdModelConf, 'DeriveFromForwardModelConfig' )
      if ( dumpFwm > 1 .and. .not. error ) stop ! error message will stop later
    end if

    if ( error ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unrecoverable errors in forward model configuration' )

  contains

    ! ............................................  Channel_Stuff  .....
    subroutine Channel_Stuff ( Channels, UsedDACSSignals )

      ! Work out which channels are used.

      type(channels_T), pointer :: Channels(:)
      integer :: UsedDACSSignals(:)          ! Indices in FwdModelConf_T%Signals
                                             ! of signals for our dacs

      integer :: Channel
      integer :: I, Ier
      integer :: NoUsedChannels
      integer :: SigInd

      noUsedChannels = 0
      do sigInd = 1, size(fwdModelConf%signals)
        noUsedChannels = noUsedChannels + &
          & count( fwdModelConf%signals(sigInd)%channels )
      end do
      allocate ( channels(noUsedChannels), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'fwdModelConf%channels' )

      ! Collect channel information from signals database.
      channel = 0
      do sigInd = 1, size(fwdModelConf%signals)
        do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
          if ( fwdModelConf%signals(sigInd)%channels(i) ) then
            channel = channel + 1
            channels(channel)%origin = &
              & lbound ( fwdModelConf%signals(sigInd)%frequencies, 1 )
            channels(channel)%used = i + channels(channel)%origin - 1
            channels(channel)%signal = sigInd
            channels(channel)%dacs = FindFirst ( usedDACSSignals, sigind )
          end if
        end do
      end do

    end subroutine Channel_Stuff

    ! ...............................................  DACS_Stuff  .....
    subroutine DACS_Stuff ( DACsStaging, UsedDACSSignals )

      use FilterShapes_m, only: DACSFilterShapes

      ! Identify which of our signals are DACS and how many unique DACS are involved
      ! Allocate and compute UsedDACSSignals and allocate DACsStaging.

      real(rp), pointer :: DACsStaging(:,:)  ! DACS radiances
      integer, pointer :: UsedDACSSignals(:) ! Indices in FwdModelConf_T%Signals
                                             ! of signals for our dacs

      integer :: I
      integer :: LBoundDACs, UBoundDACs      ! How many channels in a DAC
      integer :: NoUsedDACS
      integer :: SigInd
      logical :: SignalFlag(size(fwdModelConf%signals))

      nullify ( DACsStaging, usedDACSSignals )
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
    end subroutine DACS_Stuff

    ! ................................................  PFA_Stuff  .....
    subroutine PFA_Stuff

      integer :: B                           ! Index for beta groups
      integer :: Channel
!     logical :: Hit                         ! Hit a molecule to be used
!     integer :: HowManyPFA                  ! How many PFA data records are for the
                                             ! molecules in fwdModelConf%Molecules?
      integer :: I
!     integer :: J, K
!     integer :: NumPFA                      ! Like HowManyPFA, but for one channel
      integer :: P                           ! Index for PFA molecules in a beta group
!     integer, pointer :: PFAWork(:,:)       ! PFA Indices for a channel.
!                                            ! Sideband X Molecules.
      integer :: SB                          ! Sideband index
!     integer, pointer :: T1(:), T2(:)       ! Temps for set operations
!     integer, pointer :: WhichMolecule(:)   ! To which element of config%molecules
!                                            ! does WhichPfa(i) correspond?
!     integer, pointer :: WhichPFA(:)        ! Which PFA data records are for the
!                                            ! molecules in fwdModelConf%Molecules?

      if ( associated(pfaData) ) then
        do b = 1, size(fwdModelConf%beta_group)
          call allocate_test ( fwdModelConf%beta_group(b)%pfa_indices, &
            & fwdModelConf%sidebandStop, &
            & size(fwdModelConf%channels), &
            & size(fwdModelConf%beta_group(b)%pfa_molecules), &
            & 'Beta_group(b)%PFA_indices', moduleName, &
            & lowBound_1=fwdModelConf%sidebandStart )
          fwdModelConf%beta_group(b)%pfa_indices = 0 ! OK if some missing, but no junk
          do p = 1, size(fwdModelConf%beta_group(b)%PFA_Molecules)
            do i = PFA_by_molecule(fwdModelConf%beta_group(b)%PFA_Molecules(p)-1)+1, &
              &    PFA_by_molecule(fwdModelConf%beta_group(b)%PFA_Molecules(p))
              do sb = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
                do channel = 1, size(fwdModelConf%channels)
                  if ( matchSignal ( PFAData(sortPFAdata(i))%theSignal, &
                    &  fwdModelConf%signals(fwdModelConf%channels(channel)%signal), sideband=sb, &
                    &  channel=fwdModelConf%channels(channel)%used ) > 0 ) &
                      & fwdModelConf%beta_group(b)%pfa_indices(sb,channel,p) = &
                        & sortPFAdata(i)
                end do ! channel
              end do ! sb
            end do ! i (molecules)
          end do ! p
        end do ! b
      else
        do b = 1, size(fwdModelConf%beta_group)
          nullify ( fwdModelConf%beta_group(b)%pfa_indices ) ! Why is this needed?
        end do
      end if ! associated(pfaData)

!      if ( associated(pfaData) ) then
!
!        call allocate_test ( PFAWork, s2, size(fwdModelConf%molecules), &
!          & 'PFAWork', moduleName, lowBound_1=s1, lowBound_2=fwdModelConf%firstPFA )
!
!        ! Work out which PFA data are germane to fwdModelConf%Molecules
!        ! We do this because there are thousands of PFA data, but maybe only
!        ! a few that are germane to this fwdModelConf
!        howManyPFA = 0
!        do i = 1, size(pfaData)
!          hit = .false.
!          do j = fwdModelConf%firstPFA, size(fwdModelConf%molecules) - 1
!            if ( any(fwdModelConf%molecules(j) == PFAData(i)%molecules) ) then
!              if ( matchSignal(fwdModelConf%signals, PFAData(i)%theSignal, &
!                & DSBSSB=.true.) /= 0 ) then
!                if ( hit ) then
!                  error = .true.
!                  call MLSMessage ( MLSMSG_Warning, moduleName, &
!                    & 'Two molecules in same PFA Datum selected' )
!                  call dump ( PFAData(whichPFA(j)), details=0 )
!                end if
!                hit = .true.
!              end if
!            end if
!          end do
!          if ( hit ) howManyPFA = howManyPFA + 1
!        end do
!        call allocate_test ( whichPFA, howManyPFA, 'whichPFA', moduleName )
!        call allocate_test ( whichMolecule, howManyPFA, 'whichMolecule', moduleName )
!        howManyPFA = 0
!        do i = 1, size(pfaData)
!          do j = fwdModelConf%firstPFA, size(fwdModelConf%molecules) - 1
!            if ( any(fwdModelConf%molecules(j) == PFAData(i)%molecules) ) then
!              if ( matchSignal(fwdModelConf%signals, PFAData(i)%theSignal, &
!                & DSBSSB=.true.) /= 0 ) then
!                howManyPFA = howManyPFA + 1
!                whichPFA(howManyPFA) = i
!                whichMolecule(howManyPFA) = j
!              end if
!            end if
!          end do
!        end do
!        ! Work out PFA abstracts for each channel
!        if ( .not. error ) then
!          do i = 1, size(channels)
!            PFAWork = 0
!            nullify ( t1, channels(i)%PFAMolecules, channels(i)%betaIndex )
!            call allocate_test ( channels(i)%PFAMolecules, 0, &
!              & 'channels(i)%PFAMolecules', moduleName ) ! so we can do unions
!            ! Where there is a signal and a molecule, make whichPFA negative.
!            ! Make sure that only one molecule is selected from each PFA Datum
!            do j = 1, howManyPFA
!              if ( channels(i)%used < lbound(PFAData(whichPFA(j))%theSignal%channels,1) .or. &
!                   channels(i)%used > ubound(PFAData(whichPFA(j))%theSignal%channels,1) ) cycle
!              if ( .not. PFAData(whichPFA(j))%theSignal%channels(channels(i)%used) ) cycle
!              if ( matchSignal ( fwdModelConf%signals(channels(i)%signal), &
!                &  PFAData(whichPFA(j))%theSignal, DSBSSB=.true. ) == 0 ) cycle
!              t2 => intersection(fwdModelConf%molecules(fwdModelConf%firstPFA:), &
!                &                PFAData(whichPFA(j))%molecules)
!              if ( size(t2) == 0 ) cycle
!              whichPFA(j) = -whichPFA(j) ! Indicate above tests succeeded
!              t1 => union(channels(i)%PFAMolecules,t2)
!              call deallocate_test ( t2, 't2', moduleName )
!              call deallocate_test ( channels(i)%PFAMolecules, &
!                & 'channels(i)%PFAMolecules', moduleName )
!              channels(i)%PFAMolecules => t1
!            end do ! j = 1, howManyPFA
!            numPFA = size(channels(i)%PFAMolecules)
!            call allocate_test ( channels(i)%PFAIndex, s2, numPFA, &
!              & 'channels(i)%PFAIndex', moduleName, lowBound_1=s1 )
!            call allocate_test ( channels(i)%betaIndex, numPFA, &
!              & 'channels(i)%BetaIndex', moduleName )
!            channels(i)%betaIndex = 0 ! in case a dump is requested
!            ! Arrange the PFAData indices by sideband and molecule in PFAWork
!            do sb = s1, s2, 2
!              do j = 1, howManyPFA
!                if ( whichPFA(j) > 0 ) cycle ! some test in previous loop failed
!                if ( matchSignal ( PFAData(-whichPFA(j))%theSignal, &
!                  &  fwdModelConf%signals(channels(i)%signal), sideband=sb, &
!                  &  channel=channels(i)%used ) > 0 ) &
!                    & PFAWork(sb,whichMolecule(j)) = -whichPFA(j)
!              end do ! j = 1, howManyPFA
!            end do ! sb = s1, s2
!            ! Copy nonzero columns of PFAWork to channels(i)%PFAIndex, which
!            ! has exactly the same number of columns as the number of PFAWork's
!            ! nonzero columns
!            k = 0
!            do j = lbound(pfawork,2), ubound(pfawork,2)
!              if ( pfawork(-1,j)+pfawork(+1,j) /= 0 ) then
!                k = k + 1
!                channels(i)%PFAIndex(s1:s2:2,k) = pfawork(s1:s2:2,j)
!              end if
!            end do
!            ! Check whether we have PFA data for both sidebands if LBL is DSB
!            if ( s1 /= s2 ) then
!              do j = 1, numPFA
!                if ( channels(i)%PFAIndex(s1,j) * channels(i)%PFAIndex(s2,j) == 0 ) then
!    !             error = .true.
!                  call MLSMessage ( MLSMSG_warning, moduleName, &
!                    & 'LBL signal is DSB but both PFAs are not available' )
!                  call output ( channels(i)%used, before=' Channel ' )
!                  call output ( ', LBL signal: ' )
!                  call displaySignalName ( fwdModelConf%signals(channels(i)%signal), &
!                    & advance='yes' )
!                  call output ( ' Available ' )
!                  if ( channels(i)%PFAIndex(s1,j) + channels(i)%PFAIndex(s2,j) == 0 ) then
!                    call output ( 'PFA Datum: None', advance='yes' )
!                  else
!                    call dump ( PFAData( &
!                      & channels(i)%PFAIndex(s1,j) + channels(i)%PFAIndex(s2,j)), &
!                      & index=channels(i)%PFAIndex(s1,j) + channels(i)%PFAIndex(s2,j), &
!                      & details=0 )
!                  end if
!                end if
!              end do ! j = 1, numPFA
!            end if
!            whichPFA = abs(whichPFA)
!          end do ! i = 1, size(channels)
!        end if
!        call deallocate_test ( PFAWork, 'PFAWork', moduleName )
!        call deallocate_test ( whichPFA, 'whichPFA', moduleName )
!        call deallocate_test ( whichMolecule, 'whichMolecule', moduleName )
!     end if ! associated(pfaData)
    end subroutine PFA_Stuff

    ! ...............................  SpectroscopyCatalogExtract  .....
    subroutine SpectroscopyCatalogExtract
      use Allocate_Deallocate, only: Allocate_Test
      use Intrinsic, only: LIT_INDICES, L_NONE
      use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
        & MLSMSG_Warning
      use MLSSignals_m, only: GetRadiometerFromSignal
      use SpectroscopyCatalog_m, only: Catalog, Empty_Cat, Line_t, &
        & Lines, MostLines
      use String_table, only: GET_STRING
      use Toggles, only: Switches

      integer :: B         ! Beta_group index
      integer :: C         ! Spectroscopy catalog extract size/index
      logical :: DoThis    ! Flag for lines in catalog item
      integer :: I
      integer :: L         ! Index for lines, or number of lines
      integer, pointer :: LINEFLAG(:) ! /= 0 => Use this line
      integer :: M         ! Index for beta_group's molecule-size stuff
      character(len=32) :: MolName ! For messages
      integer :: N         ! Molecule index, L_... from Intrinsic
      integer :: NoLinesMsg = -1 ! From switches
      integer, target :: MaxLineFlag(mostLines)
      integer :: Polarized ! -1 => One of the selected lines is Zeeman split
                           ! +1 => None of the selected lines is Zeeman split
      integer :: S         ! Index for sidebands                
      integer :: STAT      ! Status from allocate or deallocate 
      type (line_T), pointer :: ThisLine
      integer :: Z         ! Index for fwdModelConf%Signals

      if ( noLinesMsg < 0 ) noLinesMsg = index(switches, '0sl') ! Done once

      ! Allocate the spectroscopy catalog extract
      c = 0
      do b = 1, size(fwdModelConf%beta_group) ! Get total catalog size
        c = c + size(fwdModelConf%beta_group(b)%lbl_molecules)
      end do

      allocate ( fwdModelConf%catalog(fwdModelConf%sidebandStart:fwdModelConf%sidebandStop,c), &
        & stat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'fwdModelConf%catalog' )

      ! Work out the spectroscopy we're going to need.
      fwdModelConf%catalog = empty_cat

      do s = s1, s2, 2
        c = 0
        do b = 1, size(fwdModelConf%beta_group)
          do m = 1, size(fwdModelConf%beta_group(b)%lbl_molecules)
            c = c + 1
            fwdModelConf%beta_group(b)%cat_index(m) = c
            n = fwdModelConf%beta_group(b)%lbl_molecules(m)
            if ( catalog(n)%molecule == l_none ) then
              call get_string ( lit_indices(n), molName )
              call MLSMessage ( MLSMSG_Error, moduleName, &
                & 'No spectroscopy catalog for ' // molName )
            end if
            fwdModelConf%catalog(s,c) = catalog(n)
            ! Don't deallocate them by mistake -- fwdModelConf%catalog is a shallow copy
            nullify ( fwdModelConf%catalog(s,c)%lines, fwdModelConf%catalog(s,c)%polarized )
            if ( associated ( catalog(n)%lines ) ) then
              ! Now subset the lines according to the signal we're using
              lineFlag => MaxLineFlag(:size(catalog(n)%lines))
              lineFlag = 0
              if ( fwdModelConf%allLinesInCatalog ) then
                ! NOTE: If allLinesInCatalog is set, then no lines can be polarized;
                ! this is checked for in ForwardModelSupport.
                lineFlag = 1
              else
                do l = 1, size ( catalog(n)%lines )
                  thisLine => lines(catalog(n)%lines(l))
                  if ( associated(thisLine%signals) ) then
                    polarized = 1 ! not polarized
                    ! Work out whether to do this line
                    do z = 1, size(fwdModelConf%signals)
                      if ( fwdModelConf%allLinesForRadiometer ) then
                        doThis = .false.
                        do i = 1, size(thisLine%signals)
                          ! Tried to make GetRadiometerFromSignal elemental, but compile time
                          ! in LF95 (optimized) for Construct.f90 went nuts! :-(
                          if ( GetRadiometerFromSignal ( thisLine%signals(i) ) == &
                            & fwdModelConf%signals(z)%radiometer ) then
                            doThis = .true.
                            if ( .not. fwdModelConf%polarized ) &
                              exit   ! loop over signals for line -- no need to check for
                            ! polarized lines
                            if ( associated(thisLine%polarized) ) then
                              if ( thisLine%polarized(i) ) then
                                polarized = -1 ! polarized
                                exit   ! loop over signals for line -- one signal
                                ! that sees a polarized line is enough to turn on
                                ! the polarized method
                              end if
                            end if
                          end if
                        end do ! End loop over signals for line
                      else
                        ! Not doing all lines for radiometer, be more selective
                        doThis = any ( &
                          & ( thisLine%signals == fwdModelConf%signals(z)%index ) .and. &
                          & ( ( thisLine%sidebands == 0 ) .or. ( thisLine%sidebands == s ) ) )
                        if ( fwdModelConf%polarized .and. doThis .and. &
                          & associated(thisLine%polarized) ) then
                          if ( any(thisLine%polarized) ) polarized = -1 ! polarized
                        end if
                      end if

                      if ( fwdModelConf%sidebandStart == fwdModelConf%sidebandStop ) &
                        & doThis = doThis .and. &
                        & any( ( thisLine%sidebands == fwdModelConf%sidebandStart ) &
                        & .or. ( thisLine%sidebands == 0 ) )
                      if ( doThis ) then
                        lineFlag(l) = polarized
                        if ( polarized < 0 .or. .not. fwdModelConf%polarized ) &
                          exit   ! loop over signals requested in fwm
                      end if
                    end do ! z End loop over signals requested in fwm
                  end if
                end do     ! l End loop over lines
              end if       ! End case where allLinesInCatalog not set

              ! Check we have at least one line for this specie

              l = count(lineFlag /= 0)
              if ( l == 0 .and. all ( fwdModelConf%catalog(s,c)%continuum == 0.0 ) &
                & .and. noLinesMsg > 0 ) then
                call get_string ( lit_indices(n), molName )
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & 'No relevant lines or continuum for '//trim(molName) )
              end if
              call allocate_test ( fwdModelConf%catalog(s,c)%lines, l, &
                & 'fwdModelConf%catalog(?,?)%lines', moduleName )
              call allocate_test ( fwdModelConf%catalog(s,c)%polarized, l, &
                & 'fwdModelConf%catalog(?,?)%polarized', moduleName )
              fwdModelConf%catalog(s,c)%lines = pack ( catalog(n)%lines, lineFlag /= 0 )
              fwdModelConf%catalog(s,c)%polarized = pack ( lineFlag < 0, lineFlag /= 0 )

            else

              ! No lines for this specie.  However, its continuum is still
              ! valid so don't set it to empty.
              ! Don't bother checking that continuum /= 0 as if it were then
              ! presumably having no continuum and no lines it wouldn't be in
              ! the catalog!
              call allocate_test ( fwdModelConf%catalog(s,c)%lines, 0, &
                & 'fwdModelConf%catalog(?,?)%lines(0)', moduleName )
              call allocate_test ( fwdModelConf%catalog(s,c)%polarized, 0, &
                & 'fwdModelConf%catalog(?,?)%polarized(0)', moduleName )
            end if
          end do ! m Molecules in fwdModelConf
        end do ! b Beta groups
      end do ! s Sidebands

    end subroutine SpectroscopyCatalogExtract

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
    ! Destroy stuff in FwdModelConf derived for one forward model run

    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type ( ForwardModelConfig_T ), intent(inout) :: FwdModelConf

    integer :: C, B, I, Ier, S

    if ( associated(fwdModelConf%channels) ) then
      do i = 1, size(fwdModelConf%channels)
        call deallocate_test ( &
          & fwdModelConf%channels(i)%PFAIndex, &
          & 'fwdModelConf%channels(i)%PFAIndex', moduleName )
        call deallocate_test ( &
          & fwdModelConf%channels(i)%PFAMolecules, &
          & 'fwdModelConf%channels(i)%PFAMolecules', moduleName )
        call deallocate_test ( &
          & fwdModelConf%channels(i)%betaIndex, &
          & 'fwdModelConf%channels(i)%BetaIndex', moduleName )
      end do
      deallocate ( fwdModelConf%channels, stat = ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'fwdModelConf%channels' )
  ! else
  !   It was already deallocated at the end of FullForwardModel
    end if

    do b = 1, size(fwdModelConf%beta_group)
      call deallocate_test ( fwdModelConf%beta_group(b)%PFA_indices, &
        'Beta_group(b)%PFA_indices', moduleName )
    end do ! b

    call deallocate_test ( fwdModelConf%DACSStaging, &
      & 'fwdModelConf%DACSStaging', moduleName )

    call deallocate_test ( fwdModelConf%usedDACSSignals, &
      & 'fwdModelConf%usedDACSSignals', moduleName )

    if ( associated(fwdModelConf%catalog) ) then
      do s = lbound(fwdModelConf%catalog,1), ubound(fwdModelConf%catalog,1), 2
        do c = 1, size(fwdModelConf%catalog,2)
          ! We don't deallocate the signals/sidebands stuff for each line because
          ! they're shallow copies of the main spectroscopy catalog stuff
          call deallocate_test ( fwdModelConf%catalog(s,c)%lines, &
            & 'fwdModelConf%catalog(?,?)%lines', moduleName )
          call deallocate_test ( fwdModelConf%catalog(s,c)%polarized, &
            & 'fwdModelConf%catalog(?,?)%polarized', moduleName )
        end do
      end do

      deallocate ( fwdModelConf%catalog, stat=s )
      if ( s /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Deallocate // 'fwdModelConf%catalog' )
  ! else
  !   It was already deallocated at the end of FullForwardModel
    end if

  end subroutine DestroyForwardModelDerived

  ! ------------------------------------ NullifyForwardModelConfig -----
  subroutine NullifyForwardModelConfig ( F )
    ! Given a forward model config, nullify all the pointers associated with it
    type ( ForwardModelConfig_T ), intent(out) :: F

    ! Executable code not needed since => NULL() initializes pointer
    ! components of intent(out) dummy arguments.
  end subroutine NullifyForwardModelConfig

  ! --------------------------------------------- PVMPackFwmConfig -----
  subroutine PVMPackFWMConfig ( config )
    use PVMIDL, only: PVMIDLPack
    use MorePVM, only: PVMPackLitIndex, PVMPackStringIndex
    use MLSSignals_m, only: PVMPackSignal
    use VGridsDatabase, only: PVMPackVGrid
    ! Dummy arguments
    type ( ForwardModelConfig_T ), intent(in) :: CONFIG
    ! Local variables
    integer :: I                        ! Loop counter

    ! Executable code
    ! First pack the lit indices
    call PVMPackStringIndex ( config%name, msg = "Packing fwmConfig name" )
    call PVMPackLitIndex ( config%cloud_der, msg = "Packing fwmConfig cloud_der" )
    call PVMPackLitIndex ( config%fwmType, msg = "Packing fwmConfig fwmType" )
    call PVMPackLitIndex ( config%i_saturation, msg = "Packing fwmConfig i_saturation" )
    call PVMPackLitIndex ( config%instrumentModule, msg = "Packing fwmConfig instrumentModule" )
    call PVMPackLitIndex ( config%windowUnits, msg = "Packing fwmConfig windowUnits" )

    ! Now pack the integer scalars
    call PVMIDLPack ( (/ &
      & config%linearSideband, &
      & config%no_cloud_species, config%no_model_surfs, &
      & config%num_ab_terms, config%num_azimuth_angles, &
      & config%num_scattering_angles, config%num_size_bins, &
      & config%sidebandStart, config%sidebandStop, &
      & config%surfaceTangentIndex /), msg = "Packing fwmConfig integers" )

    ! Now the logical scalars
    call PVMIDLPack ( (/ &
      & config%allLinesForRadiometer, config%allLinesInCatalog, config%anyPFA, &
      & config%atmos_der, config%default_spectroscopy, config%differentialScan,&
      & config%do_1d, config%do_baseline, config%do_conv, config%do_freq_avg, &
      & config%forceFoldedOutput, config%forceSidebandFraction, &
      & config%globalConfig, config%incl_cld, config%lockBins, config%polarized, &
      & config%skipOverlaps, config%spect_Der, config%switchingMirror, &
      & config%temp_Der /), msg = "Packing fwmConfig logicals" )

    ! Now pack the reals
    call PVMIDLPack ( (/ config%phiWindow, config%tolerance /), &
      & msg = "Packing fwmConfig reals" )

    ! ------------- The rest are arrays and/or types
    ! Bin selectors
    if ( associated ( config%binSelectors ) ) then
      call PVMIDLPack ( size ( config%binSelectors ), msg = "Packing number of binSelectors" )
      call PVMIDLPack ( config%binSelectors, msg = "Packing binSelectors" )
    else
      call PVMIDLPack ( 0, msg = "Packing 0 binSelectors" )
    end if

    ! Molecules / derivatives
    if ( associated ( config%molecules ) ) then
      call PVMIDLPack ( size ( config%molecules ), msg = "Packing number of molecules" )
      if ( size ( config%molecules ) > 0 ) then
        do i = 1, size(config%molecules) - 1
          call PVMPackLitIndex ( abs ( config%molecules(i) ), &
            & msg = "Packing a molecule" )
          call PVMIDLPack ( (config%molecules(i) > 0.0), msg = "Packing molecule sign" )
        end do
        call PVMIDLPack ( config%moleculeDerivatives, msg = "Packing molecule derivatives" )
      end if
    else
      call PVMIDLPack ( 0, msg = "Packing 0 molecules" )
    end if

    ! Specific quantities
    if ( associated ( config%specificQuantities ) ) then
      call PVMIDLPack ( size ( config%specificQuantities ), &
        & msg = "Packing number of specificQuantities" )
      call PVMIDLPack ( config%specificQuantities, msg = "Packing specificQuantities" )
    else
      call PVMIDLPack ( 0, msg = "Packing 0 specificQuantities" )
    end if

    ! Pack the other structures - signals
    if ( associated ( config%signals ) ) then
      call PVMIDLPack ( size ( config%signals ), msg = "Packing number of signals" )
      do i = 1, size ( config%signals )
        call PVMPackSignal ( config%signals(i) )
      end do
    else
      call PVMIDLPack ( 0, msg = "Packing 0 signals" )
    end if

    ! Vgrids
    call PVMIDLPack ( (/ associated ( config%integrationGrid ), &
      & associated ( config%tangentGrid ) /), msg = "Packing vGrid flags" )
    if ( associated ( config%integrationGrid ) ) &
      & call PVMPackVGrid ( config%integrationGrid )
    if ( associated ( config%tangentGrid ) ) &
      & call PVMPackVGrid ( config%tangentGrid )

  end subroutine PVMPackFWMConfig

  ! ----------------------------------------- PVMUnpackFWMConfig ---------
  subroutine PVMUnpackFWMConfig ( CONFIG )
    use PVMIDL, only: PVMIDLUnpack
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
    logical, dimension(20) :: LS        ! Temporary array, for logical scalars
    integer, dimension(10) :: IS        ! Temporary array, for integer scalars
    real(r8), dimension(2) :: RS        ! Temporary array, for real scalars
    integer :: I                        ! Loop counter
    integer :: N                        ! Array size

    ! Executable code
    ! First unpack the lit indices
    call PVMUnpackStringIndex ( config%name, msg = "Unpacking fwmConfig name" )
    call PVMUnpackLitIndex ( config%cloud_der, msg = "Unpacking fwmConfig cloud_der" )
    call PVMUnpackLitIndex ( config%fwmType, msg = "Unpacking fwmConfig fwmType" )
    call PVMUnpackLitIndex ( config%i_saturation, msg = "Unpacking fwmConfig i_saturation" )
    call PVMUnpackLitIndex ( config%instrumentModule, &
      msg = "Unpacking fwmConfig instrumentModule" )
    call PVMUnpackLitIndex ( config%windowUnits, msg = "Unpacking fwmConfig windowUnits" )

    ! Now the integer scalars
    call PVMIDLUnpack ( is, msg = "Unpacking fwmConfig integers" )
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
    call PVMIDLUnpack ( ls, msg = "Unpacking fwmConfig logicals" )
    i = 1
    config%allLinesForRadiometer = ls(i) ; i = i + 1
    config%allLinesInCatalog     = ls(i) ; i = i + 1
    config%anyPFA                = ls(i) ; i = i + 1
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
    config%switchingMirror       = ls(i) ; i = i + 1
    config%temp_der              = ls(i) ; i = i + 1

    ! Now the real scalars
    call PVMIDLUnpack ( rs, msg = "Unpacking fwmConfig reals" )
    i = 1
    config%phiWindow = rs(i) ; i = i + 1
    config%tolerance = rs(i) ; i = i + 1

    ! ------- The rest are arrays and/or types
    ! Bin selectors
    call PVMIDLUnpack ( n, msg = "Unpacking number of specific quantities" )
    if ( n > 0 ) then
      call Allocate_test ( config%binSelectors, n, &
        & 'config%binSelectors', ModuleName )
      call PVMIDLUnpack ( config%binSelectors, msg = "Unpacking binSelectors" )
    end if

    ! Molecules / derivatives
    call PVMIDLUnpack ( n, msg = "Unpacking number of molecules" )
    if ( n > 0 ) then
      call Allocate_test ( config%molecules, n, 'config%molecules', ModuleName )
      call Allocate_test ( config%moleculeDerivatives, &
        & n, 'config%moleculeDerivatives', ModuleName )
      do i = 1, n - 1
        call PVMUnpackLitIndex ( config%molecules(i), msg = "Unpacking a molecule" )
        call PVMIDLUnpack ( flag, msg = "Unpacking a molecule sign flag" )
        if ( .not. flag ) config%molecules(i) = - config%molecules(i)
      end do
      config%molecules(n) = huge(config%molecules(n)) ! Sentinel
      call PVMIDLUnpack ( config%moleculeDerivatives, msg = "Unpacking moleculeDerivatives" )
    end if

    ! Specific quantities
    call PVMIDLUnpack ( n, msg = "Unpacking number of specific quantities" )
    if ( n > 0 ) then
      call Allocate_test ( config%specificQuantities, n, &
        & 'config%specificQuantities', ModuleName )
      call PVMIDLUnpack ( config%specificQuantities, msg = "Unpacking specific quantities" )
    end if

    ! Unpack other structures - signals
    call PVMIDLUnpack ( n, msg = "Unpacking number of signals" )
    if ( n > 0 ) then
      allocate ( config%signals(n), STAT=info )
      if ( info /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'config%signals' )
      do i = 1, n
        call PVMUnpackSignal ( config%signals(i) )
      end do
    end if

    ! Vgrids
    call PVMIDLUnpack ( ls(1:2), msg = "Unpacking vGrid flags" )
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
    integer :: B                        ! Subscript for Beta_Group
    integer :: STATUS                   ! Flag from allocate etc.
    logical :: MYDEEP                   ! Copy of deep

    ! Executable code
    myDeep = .false.
    if ( present ( deep ) ) myDeep = deep
    call destroyForwardModelDerived ( config )
    ! Destroy the beta groups
    do b = 1, size(config%beta_group)
      call deallocate_test ( config%beta_group(b)%cat_index, 'Cat_Index', moduleName )
      call deallocate_test ( config%beta_group(b)%lbl_molecules, 'LBL_Molecules', moduleName )
      call deallocate_test ( config%beta_group(b)%lbl_ratio, 'LBL_Ratio', moduleName )
      call deallocate_test ( config%beta_group(b)%pfa_molecules, 'PFA_Molecules', moduleName )
      call deallocate_test ( config%beta_group(b)%pfa_ratio, 'PFA_Ratio', moduleName )
      ! PFA_Indices are created and destroyed on each call to FullForwardModel.
    end do
    deallocate ( config%beta_group, stat=b )
    if ( b /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'Config%Beta_group' )
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
    ! Don't deallocate config%molecules because it's a pointer into beta_group
    call Deallocate_test ( config%specificQuantities, &
      & "config%specificQuantities", ModuleName )
    call Deallocate_test ( config%binSelectors, &
      & "config%binSelectors", ModuleName )
  end subroutine DestroyOneForwardModelConfig

  ! --------------------------------------------  Dump_Beta_Group  -----
  subroutine Dump_Beta_Group ( Beta_Group, Name )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_indices
    use Output_m, only: NewLine, Output
    use String_Table, only: Display_String

    type(beta_group_t), intent(in) :: Beta_Group(:)
    character(len=*), intent(in), optional :: Name

    integer :: I, J

    call output ( '  Beta groups' )
    if ( present(name) ) call output ( ' ' // trim(name) )
    call output ( size(beta_group), before=', SIZE = ', advance='yes' )
    do i = 1, size(beta_group)
      call output ( i, before='  Beta group ', after=': ' )
      call display_string ( lit_indices(beta_group(i)%molecule) )
      if ( beta_group(i)%derivatives ) call output ( ' with derivative' )
      if ( size(beta_group(i)%lbl_molecules) > 0 ) call display_string ( &
          & lit_indices(beta_group(i)%lbl_molecules), before=', LBL:' )
      if ( size(beta_group(i)%pfa_molecules) > 0 ) call display_string ( &
          & lit_indices(beta_group(i)%pfa_molecules), before=', PFA:' )
      call newLine
      if ( size(beta_group(i)%lbl_molecules) > 0 ) then
        call dump ( beta_group(i)%lbl_ratio, name='   LBL Ratio' )
        call dump ( beta_group(i)%cat_index, name='   Cat_Index' )
      end if
      if ( size(beta_group(i)%pfa_molecules) > 0 ) then
        call dump ( beta_group(i)%pfa_ratio, name='   PFA Ratio' )
        if ( associated(beta_group(i)%PFA_indices) ) then
          do j = lbound(beta_group(i)%PFA_indices,1), &
                 ubound(beta_group(i)%PFA_indices,1), 2
            call output ( j, before='   PFA Indices for sideband ', advance='yes' )
            call dump ( beta_group(i)%PFA_indices(j,:,:) )
          end do
        end if
      end if
      if ( associated ( beta_group(i)%qty%qty ) ) call dump ( beta_group(i)%qty )
    end do
  end subroutine Dump_Beta_Group

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
  subroutine Dump_ForwardModelConfig ( Config, Where )

    use Dump_0, only: DUMP
    use Intrinsic, only: Lit_indices
    use MLSSignals_M, only: GetNameOfSignal, MaxSigLen
    use Output_M, only: NewLine, Output
    use String_Table, only: Display_String

    type (ForwardModelConfig_T), intent(in) :: Config
    character(len=*), optional, intent(in) :: Where

    ! Local variables
    integer ::  I, J                         ! Loop counters
    character (len=MaxSigLen) :: SignalName  ! A line of text

    ! executable code

    call output ( '  Forward Model Config Name: ' )
    call display_string ( Config%name )
    if ( present(where) ) then
      call output ( ' from ' )
      call output ( where )
    end if
    call newLine
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
    if ( associated(Config%Beta_group) ) call dump ( Config%Beta_group )
    call output ( '  Molecules: ', advance='yes' )
    if ( associated(Config%molecules) ) then
      i = 0
      do j = 1, size(Config%molecules) - 1
        call output ( '    ' )
        if ( Config%molecules(j) < 0 ) call output ( '-' )
        call display_string(lit_indices(abs(Config%molecules(j))))
        if ( Config%molecules(j) > 0 ) then
          i = i + 1
          call output ( i, before=':' )
        end if
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
    if ( associated(Config%usedDACSSignals) ) &
      call dump  ( Config%usedDACSSignals, name='   Used DACS Signals' )
    if ( associated(Config%channels) ) then
      call output ( '   Channel info:', advance='yes' )
      do j = 1, size(Config%channels)
        call output ( Config%channels(j)%used, before='    Used: ' )
        call output ( Config%channels(j)%origin, before='    Origin: ' )
        call output ( Config%channels(j)%signal, before='    Signal: ' )
        call output ( Config%channels(j)%DACS, before='    DACS: ', &
          & advance='yes' )
        if ( associated(Config%channels(j)%PFAIndex) ) &
          & call dump ( Config%channels(j)%PFAIndex(-1:1:2,:), &
            & name='    PFAData' )
        if ( associated(Config%channels(j)%PFAMolecules) ) then
          call output ( '    PFA Molecules:', advance='no' )
          if ( size(Config%channels(j)%PFAMolecules) == 0 ) &
            call output ( ' Empty' )
          do i = 1, size(Config%channels(j)%PFAMolecules)
            call output ( ' ', advance='no' )
            call display_string ( &
              & lit_indices(Config%channels(j)%PFAMolecules(i)), &
              & advance='no' )
            if ( Config%channels(j)%betaIndex(i) /= 0 ) &
              & call output ( Config%channels(j)%betaIndex(i), &
                & before=':' )
          end do
          call newLine
        end if
      end do
    end if
  end subroutine Dump_ForwardModelConfig

  ! ---------------------------------------------  Dump_Qty_Stuff  -----
  subroutine Dump_Qty_Stuff ( Qty )
    use Output_m, only: NewLine, Output
    use VectorsModule, only: Dump
    type(qtyStuff_t), intent(in) :: Qty
    call dump ( qty%qty, details=-2 )
    if ( qty%foundInFirst ) call output ( ', Found in first' )
    call newLine
  end subroutine Dump_Qty_Stuff

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ForwardModelConfig

! $Log$
! Revision 2.64  2004/12/13 20:35:22  vsnyder
! Moved a bunch of stuff that doesn't depend on the state vector here from
! get_species data.  Some support for PFA, too.
!
! Revision 2.63  2004/11/05 19:38:11  vsnyder
! Moved some stuff here from Get_Species_Data, more work on PFA, some rearranging
!
! Revision 2.62  2004/11/04 03:42:08  vsnyder
! Provide for both LBL_Ratio and PFA_Ratio in beta_group
!
! Revision 2.61  2004/11/03 01:25:30  vsnyder
! Don't deallocate config%molecules -- it's a pointer into beta_group
!
! Revision 2.60  2004/11/01 20:18:23  vsnyder
! Reorganization of representation for molecules and beta groups
!
! Revision 2.59  2004/10/06 21:23:50  vsnyder
! More work on PFA data structures
!
! Revision 2.58  2004/09/01 00:32:40  vsnyder
! Add PFAIndex field, polish up the dump routines
!
! Revision 2.57  2004/08/07 01:17:26  vsnyder
! Forgot an advance='yes' in a dump routine
!
! Revision 2.56  2004/08/05 20:57:52  vsnyder
! Put a sentinel at the end of %molecules
!
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
