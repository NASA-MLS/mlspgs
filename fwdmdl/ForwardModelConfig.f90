! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject bto U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module ForwardModelConfig
!=============================================================================

! Set up the forward model configuration, except for actually processing
! the command.

  use MLSKinds, only: R8, RP
  use MLSSignals_M, only: Signal_T
  use SpectroscopyCatalog_M, only: Catalog_T
  use VectorsModule, only: RV, VectorValue_T
  use VGridsDatabase, only: VGrid_T, DestroyVGridContents, Dump

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
  public :: PVMPackFWMConfig, PVMUnpackFWMConfig
  public :: StripForwardModelConfigDatabase
 
  ! Public Types:

  ! Quantities derived from forward models, but not carted around by
  ! PVMPackFWMConfig and PVMUnpackFWMConfig.  Rather, they are computed by
  ! DeriveFromForwardModelConfig, either when a ForwardModelConfig_T is
  ! created or when one arrives by way of PVMUnpackFWMConfig.

  type, public :: QtyStuff_T ! So we can have an array of pointers to QTY's
    type (VectorValue_T), pointer :: Qty => NULL()
    real(rv), pointer :: Value1(:) => NULL()      ! For logBasis quantities
                                                  ! Freq * Vert * Horiz
    real(rv), pointer :: Values(:,:) => NULL()    ! Associated with Value1
                                                  ! Freq * Vert x Horiz
    real(rv), pointer :: Values3(:,:,:) => NULL() ! Associated with Value1
                                                  ! Freq x Vert x Horiz
    logical :: FoundInFirst = .false.
    logical :: WasSpecific = .false.
    logical :: DerivOK = .false. ! There is a place for the forward model to
                                 ! put the derivatives.  If .not. FoundInFirst,
                                 ! then the forward model must have gotten
                                 ! ExtraJacobian.
  end type QtyStuff_T

  ! Data structures to indicate which spectral parameters are being solved.
  type, public :: SpectroParam_T
    integer :: Beta(2)       ! Beta group LBL indices
    integer :: Molecule
    type (QtyStuff_T) :: Qty ! The Qty's vector
  end type SpectroParam_T

  ! Subscripts for spectral parameter stuff
  integer, parameter, public :: LineCenter = 1, LineWidth = 2, LineWidth_TDep = 3

  ! Beta group type declaration.  Each entry in the Molecules list of the form
  ! "m" has one of these with n_elements == 1, referring to "m".  Each entry
  ! of the form "[m,m1,...,mn]" has one of these with n_elements == n, referring
  ! to m1,...mn (but not m).  The Ratio components are filled in
  ! Get_Species_Data.
  type, public :: LBL_T ! For LBL molecules in the group.  All the same
    ! size, if associated at all.
    integer, pointer  :: Cat_Index(:) => NULL() ! Indices for config%catalog
                                      ! and gl_slabs. Allocated in
                                      ! ConstructForwardModelConfig. Filled in
                                      ! DeriveFromForwardModelConfig.
    integer, pointer :: Molecules(:) => NULL() ! LBL molecules in the group
                                      ! if a group, i.e., "m1...mn", else "m" if
                                      ! "m" is LBL, else zero size.
    real(rp), pointer :: Ratio(:) => NULL() ! Isotope ratio.  Allocated in
                                      ! ForwardModelSupport with value 1.0, but
                                      ! could be filled in Get_Species_Data.
    integer :: Spect_der_ix(lineCenter:lineWidth_tDep) = 0 ! Do deriv w.r.t.
      ! spectral param if nonzero. If Spect_der_ix(lineCenter) is nonzero,
      ! fwdModelConf%%LineCenter(Spect_der_ix_ix(lineCenter)) is added onto the
      ! line center in the gl_slabs array for the i'th line-by-line molecule in
      ! sideband s.  Similarly for %LineWidth and %LineWidth_Tdep.
  end type LBL_T

  type, public :: PFA_T ! For PFA molecules in the group:
    integer, pointer :: Molecules(:) => NULL() ! PFA molecules in the group if a
                                      ! group, i.e., "m1...mn", else "m" if "m"
                                      ! is PFA, else zero size.
    integer, pointer :: Data(:,:) => NULL() ! Channels x 1:size(Molecules).
                                      ! Indices in PFAData.  Allocated and
                                      ! filled in DeriveFromForwardModelConfig.
    real(rp), pointer :: Ratio(:) => NULL() ! 1:size(Molecules).  Isotope
                                      ! ratio.  Allocated in ForwardModelSupport
                                      ! with value 1.0, but could be filled in
                                      ! Get_Species_Data.
  end type PFA_T

  type, public :: Beta_Group_T
    ! For the group as a whole:
    logical :: Derivatives = .false.        ! "Compute derivatives w.r.t. mixing ratio"
    logical :: SecondDerivatives = .false.  ! "Compute second derivatives w.r.t. mixing ratio"
    logical :: Group = .false.        ! "Molecule group", i.e., [m,m1,...,mn]
    integer :: Molecule               ! Group name, i.e., "m".
    type(qtyStuff_t) :: Qty           ! The Qty's vector and foundInFirst, filled
                                      ! in Get_Species_Data.
    type(lbl_t) :: LBL(2)             ! LSB, USB for Line-By-Line stuff
    type(pfa_t) :: PFA(2)             ! LSB, USB for Pre-Frequency-Averaged stuff
  end type Beta_Group_T

  ! Channel information from the signals database.
  type, public :: Channels_T
    integer :: Used       ! Which channel is this?
    integer :: Origin     ! Index of first channel (zero or one)
    integer :: Signal     ! Signal index for the channel
    integer :: DACS       ! DACS index if any, else zero
    integer :: ShapeInds(2) ! Filter shape indices, by sideband, 1 => LSB, 2 => USB
  end type Channels_T

  ! The scalar components are sorted in the order they are to make the packing
  ! and unpacking for PVM as easy as possible to maintain
  type, public :: ForwardModelConfig_T
    ! First the lit_indices
    integer :: Name                   ! String index of config name
    integer :: Where                  ! Tree node index of config (for messages)
    integer :: Cloud_der              ! Compute cloud sensitivity in cloud models.
                                      ! l_iwc_low_height, l_iwc_high_height, l_iwp
                                      ! l_none
    integer :: FwmType                ! l_linear, l_full, l_scan, ....
    integer :: I_saturation           ! Flag to determine saturation status
                                      ! l_clear, l_clear_110rh_below_top
                                      ! l_clear_0rh, l_clear_lowest_0_110rh
                                      ! l_clear_110rh_below_tropopause,
                                      ! l_cloudy_110rh_below_top
                                      ! l_cloudy_110rh_in_cloud,
                                      ! l_cloudy_nearside_only
    integer :: InstrumentModule       ! Module for scan model (actually a spec index)
    integer :: Trapezoid              ! "Wrong" or "Correct"
    integer :: WindowUnits            ! Either degrees or profiles
    ! Now the other integers
    integer :: Cat_Size(2)            ! Catalog size, by sideband, 1 = LSB, 2 = USB
    integer :: LinearSideband         ! For hybrid model, which SB is linear?
    integer :: MIFTangent             ! L_ECRtoFOV or L_PTAN
    integer :: No_cloud_species       ! No of Cloud Species '2'
    integer :: No_model_surfs         ! No of Model surfaces '640'
    integer :: NoUsedChannels         ! Total in all signals
    integer :: Ntimes = 0             ! Number of times calling FullForwardModel
    integer :: Num_ab_terms           ! No of AB terms '50'
    integer :: Num_azimuth_angles     ! No of azmuth angles '8'
    integer :: Num_scattering_angles  ! No of scattering angles '16'
    integer :: Num_size_bins          ! No of size bins '40'
    integer :: ReferenceMIF = 1       ! MIF number to use for MAF geolocation
    integer :: SidebandStart, SidebandStop ! Folded or SSB config?
    integer :: SurfaceTangentIndex    ! Index in Tangentgrid of Earth's surface
    integer :: TScatMIF               ! Which MIF to use for TScat LOS VEL and PHITAN
    integer :: xStar                  ! Index of specific vector to use for linearized model
    integer :: yStar                  ! Index of specific vector to use for linearized model
    ! Now the logicals
    logical :: AllLinesForRadiometer  ! As opposed to just using lines designated for band.
    logical :: AllLinesInCatalog      ! Use all the lines
    logical :: AnyLBL(2)              ! "there are LBL molecules in the sideband"
    logical :: AnyPFA(2)              ! "there are PFA molecules in the sideband"
    logical :: Atmos_der              ! Do atmospheric derivatives
    logical :: Atmos_second_der       ! Do atmospheric second derivatives
    logical :: Default_spectroscopy   ! Using Bill's spectroscopy data
    logical :: DifferentialScan       ! Differential scan model
    logical :: Do_1d                  ! Do 1D forward model calculation
    logical :: Do_baseline            ! Do a baseline computation
    logical :: Do_conv                ! Do convolution
    logical :: Do_freq_avg            ! Do Frequency averaging
    logical :: Do_Path_Norm           ! Do path normalization
    logical :: ForceFoldedOutput      ! Output to folded sideband even if signal is other (linear only)
    logical :: ForceSidebandFraction  ! If set mult. by SBfrac even if single sideband
    logical :: GenerateTScat          ! Generate TScat tables
    LOGICAL :: GlobalConfig           ! If set is shared between all chunks
    LOGICAL :: hmag_der               ! magnetic field magnitude derivative
    LOGICAL :: htheta_der             ! magnetic field theta angle derivative
    logical :: hphi_der               ! magnetic field phi angle derivative
    logical :: IgnoreHessian          ! Don't do 2nd-order Taylor series
                                      ! in quasi-linear model even if L2PC has
                                      ! a Hessian
    logical :: Incl_cld ! Include cloud extinction calculation in Bill's forward model
    logical :: IsRadianceModel        ! The forward model is a radiance model
    logical :: LockBins               ! Use same l2pc bin for whole chunk
    logical :: No_Magnetic_Field      ! Set magnetic field to zero for testing
    logical :: Polarized              ! Use polarized model for Zeeman-split lines
    logical :: Refract                ! Compute refractive correction for PhiTan
    logical :: ScanAverage            ! Average scan over MIF
    logical :: SkipOverlaps           ! Don't calculate for MAFs in overlap regions
    logical :: Spect_Der              ! Do spectroscopy derivatives
    logical :: SwitchingMirror        ! Model radiance at the switching mirror
    logical :: Temp_Der               ! Do temperature derivatives
    logical :: TransformMIFExtinction ! Transform MIF extinction, see wvs-107
    logical :: TransformMIFRHI        ! Transform MIF RHI
    logical :: UseTScat               ! Use TScat tables + linear model in full model
    ! Now the reals
    real (r8) :: FrqTol               ! MHz, how close to desired frequency must
                                      ! Mie table be?
    real (r8) :: PhiWindow(2)         ! Window for examining stuff; (1) before
                                      ! tangent point, (2) after tangent point
    real (r8) :: Tolerance            ! Accuracy desired when choosing approximations
    real :: sum_DeltaTime = 0.0       ! sum of delta time calling FullForwardModel 
    real :: sum_squareDeltaTime = 0.0 ! sum of the square of delta times calling FullForwardModel 
    ! Now the arrays
    integer, dimension(:), pointer :: BinSelectors=>NULL() ! List of relevant bin selectors
    integer, dimension(:), pointer :: Molecules=>NULL() ! Which molecules to consider
    logical, dimension(:), pointer :: MoleculeDerivatives=>NULL() ! Want Jacobians
    logical, dimension(:), pointer :: MoleculeSecondDerivatives=>NULL() ! Want Hessians
    integer, dimension(:), pointer :: SpecificQuantities=>NULL() ! Specific quantities to use
    ! Now the derived types
    type (beta_group_t), dimension(:), pointer :: Beta_Group => NULL() ! q.v. above
    type (spectroParam_t), dimension(:), pointer :: LineCenter => NULL()
    type (spectroParam_t), dimension(:), pointer :: LineWidth => NULL()
    type (spectroParam_t), dimension(:), pointer :: LineWidth_TDep => NULL()
      ! Some of this is filled in each time the forward model is invoked.  Those
      ! parts aren't dragged around by PVM.
    type (Signal_T), dimension(:), pointer :: Signals=>NULL()
    integer, dimension(:), pointer :: SignalIndices=>NULL() ! in signals database
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
    type(qtyStuff_T) :: Temp ! Temperature stuff
  end type ForwardModelConfig_T

  !------------- RCS Ident Info (more below in not_used_here) ----------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile: ForwardModelConfig.f90,v $"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  ! ----------------------------  AddForwardModelConfigToDatabase  -----
  integer function AddForwardModelConfigToDatabase ( Database, Item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

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

    use Allocate_Deallocate, only: Allocate_Test
    use MLSMessageModule, only: MLSMessage,  MLSMsg_Error, MLSMsg_Warning
    use MLSStringLists, only: SwitchDetail
    use PFADatabase_M, only: Dump
    use Read_Mie_M, only: F_S
    use Output_M, only: Output
    use SpectroscopyCatalog_M, only: Dump
    use Toggles, only: Switches

    type (ForwardModelConfig_T), intent(inout) :: FwdModelConf
    integer :: DumpFwm = -2                ! -2 = not called yet, -1 = no dumps,
                                           ! low-order digit: catalog dump level
                                           ! high-order digit: 1 => stop
    logical :: Error
    integer :: S1, S2                      ! SidebandStart, SidebandStop

    if ( dumpFwm < -1 ) dumpFwm = switchDetail(switches,'fwmd')

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

    ! Work out the PFA stuff.  The PFA stuff is done here instead of in
    ! ForwardModelSupport because the PFA stuff might be large (at least
    ! compared to the LBL stuff), so it is useful to allocate and destroy
    ! it separately for each forward model run.
    call PFA_Stuff ! Below

    if ( fwdModelConf%useTScat .and. .not. allocated(f_s) ) then
      call output ( 'UseTScat requested but no Mie tables loaded', advance='yes' )
      error = .true.
    end if

    if ( dumpFwm > -1 .or. error ) then
      call dump ( fwdModelConf, 'DeriveFromForwardModelConfig' )
      call dump ( fwdModelConf%catalog, details=mod(dumpFwm,10) )
      if ( dumpFwm > 9 .and. .not. error ) stop ! error message will stop later
    end if

    if ( error ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unrecoverable errors in forward model configuration' )

  contains

    ! ............................................  Channel_Stuff  .....
    subroutine Channel_Stuff ( Channels, UsedDACSSignals )

      ! Work out which channels are used.

      use Allocate_Deallocate, only: Test_Allocate
      use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
      use FilterShapes_M, only: DACSFilterShapes, FilterShapes
      use MLSFinds, only: FindFirst
      use MLSSignals_M, only: MatchSignal

      type(channels_T), pointer :: Channels(:)
      integer :: UsedDACSSignals(:)          ! Indices in FwdModelConf_T%Signals
                                             ! of signals for our dacs

      integer(c_intptr_t) :: Addr         ! For tracing
      integer :: Channel
      integer :: I, Ier
      integer :: SigInd
      integer :: SX, ThisSideband ! Sideband indices

      allocate ( channels(fwdModelConf%noUsedChannels), stat=ier )
      if ( ier == 0 .and. fwdModelConf%noUsedChannels > 0 ) &
        & addr = transfer(c_loc(channels(1)), addr)
      call test_allocate ( ier, ModuleName, 'info%channels', &
        & ubounds=(/ fwdModelConf%noUsedChannels /), &
        & elementSize = storage_size(channels) / 8, address=addr )

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
            ! Look for channel shape information, but we only need to
            ! do this if we've been asked to perform frequeny averaging
            if ( fwdModelConf%do_freq_avg ) then
              if ( channels(channel)%dacs == 0 ) then
                do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
                  sx = (thisSideband +3) / 2
                  channels(channel)%shapeInds(sx) = MatchSignal ( &
                    & filterShapes%signal, fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channels(channel)%used )
                  if ( channels(channel)%shapeInds(sx) == 0) &
                    & call MLSMessage ( MLSMSG_Error, ModuleName, &
                    &    "No matching channel shape information" )
                end do
              else
                do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
                  sx = (thisSideband +3) / 2
                  channels(channel)%shapeInds(sx) = MatchSignal ( &
                    & DACSFilterShapes%filter%signal, fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channels(channel)%used )
                  if ( channels(channel)%shapeInds(sx) == 0 ) &
                    & call MLSMessage ( MLSMSG_Error, ModuleName, &
                    &    "No matching DACS channel shape information" )
                end do
              end if ! filter bank or DACS
            end if ! Need filter shape information
          end if
        end do
      end do

      if ( fwdModelConf%do_freq_avg ) then
        ! Check again whether we have channel shapes.  We can't just check in
        ! ConstructForwardModelConfig whether frequency averaging was requested
        ! and both FilterShapes and DACSFilterShapes are associated, because if
        ! there are no DACS (or only DACS) then one of them might be harmlessly
        ! disassociated.  To check in ConstructForwardModelConfig, it would be
        ! necessary to determine whether any channels are DACS channels, which
        ! is only done here.
        ier = 0
        do i = 1, channel
          do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
            sx = (thisSideband + 3) / 2
            if ( channels(channel)%shapeInds(sx) == 0 ) ier = 1
          end do
        end do
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
           & "Frequency averaging requested but there is no channel " // &
           & "shape information for some channel" )
      end if

    end subroutine Channel_Stuff

    ! ...............................................  DACS_Stuff  .....
    subroutine DACS_Stuff ( DACsStaging, UsedDACSSignals )

      use FilterShapes_M, only: DACSFilterShapes

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

      use Allocate_Deallocate, only: Allocate_Test
      use Intrinsic, only: Lit_Indices
      use MLSMessageModule, only: MLSMessage, MLSMsg_Error
      use MLSSignals_M, only: DisplaySignalName
      use Molecules, only: L_CloudIce
      use MoreTree, only: StartErrorMessage
      use PFADatabase_M, only: Test_And_Fetch_PFA
      use Read_Mie_M, only: Beta_C_A, Beta_C_S
      use String_Table, only: Display_String

      integer :: B                    ! Index for beta groups
      integer :: Channel              ! Index in fwdModelConf%channels
      integer :: P                    ! Index for PFA molecules in a beta group
      integer :: SB                   ! Sideband index, -1 .. +1
      integer :: SX                   ! Sideband index, 1 .. 2

      ! Fill fwdModelConf%beta_group%pfa%data

      do sb = s1, s2, 2
        sx = ( sb + 3 ) / 2
        do b = 1, size(fwdModelConf%beta_group)
          call allocate_test ( fwdModelConf%beta_group(b)%pfa(sx)%data, &
            & size(fwdModelConf%channels), &
            & size(fwdModelConf%beta_group(b)%pfa(sx)%molecules), &
            & 'Beta_group(b)%PFA(sx)%data', moduleName, fill=0 )
          do p = 1, size(fwdModelConf%beta_group(b)%pfa(sx)%molecules)
            if ( fwdModelConf%beta_group(b)%pfa(sx)%molecules(p) == l_cloudIce ) then
              if ( .not. allocated(beta_c_a) ) &
                call MLSMessage ( MLSMSG_Error, moduleName, &
                  'No Mie tables for Cloud_A beta' )
              cycle
              if ( .not. allocated(beta_c_s) )  &
                call MLSMessage ( MLSMSG_Error, moduleName, &
                  'No Mie tables for Cloud_S beta' )
              cycle
            end if
            do channel = 1, size(fwdModelConf%channels)
              ! Look up PFA data and read it if necessary
              fwdModelConf%beta_group(b)%pfa(sx)%data(channel,p) = &
                test_and_fetch_PFA(fwdModelConf%beta_group(b)%pfa(sx)%molecules(p), &
                  & fwdModelConf%signalIndices(fwdModelConf%channels(channel)%signal), &
                  & sb, fwdModelConf%channels(channel)%used, fwdModelConf%spect_der)
              if ( fwdModelConf%beta_group(b)%pfa(sx)%data(channel,p) == 0 ) then
                call startErrorMessage ( fwdModelConf%where )
                call display_string ( &
                  & lit_indices(fwdModelConf%beta_group(b)%molecule), &
                  & before=' PFA table not found for ' )
                call displaySignalName ( &
                  & fwdModelConf%signals(fwdModelConf%channels(channel)%signal), &
                  & advance='yes', before=' and ', sideband=sb, &
                  & channel=fwdModelConf%channels(channel)%used )
                error = .true.
              end if
            end do ! channel
          end do ! p
        end do ! b
      end do ! sb

    end subroutine PFA_Stuff

    ! ...............................  SpectroscopyCatalogExtract  .....
    subroutine SpectroscopyCatalogExtract
      use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
      use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
      use Intrinsic, only: Lit_Indices, L_None
      use MLSSignals_M, only: GetRadiometerFromSignal
      use MLSStringLists, only: SwitchDetail
      use MoreTRee, only: StartErrorMessage
      use SpectroscopyCatalog_M, only: Catalog, Empty_Cat, Lines, MostLines
      use String_Table, only: Display_String ! , Get_String
      use Toggles, only: Switches

      integer(c_intptr_t) :: Addr         ! For tracing
      integer :: B         ! Beta_group index
      integer :: C         ! Spectroscopy catalog extract size/index
      logical :: DoThis    ! Flag for lines in catalog item
      integer :: I
      integer :: L         ! Index for lines, or number of lines
      integer, pointer :: LineFlag(:) ! /= 0 => Use this line
      integer :: M         ! Index for beta_group's molecule-size stuff
      integer :: N         ! Molecule index, L_... from Intrinsic
      integer :: NoLinesMsg = -1 ! From switches
      integer, target :: MaxLineFlag(mostLines)
!     character(len=32) :: MoleculeName ! for error message
      integer :: Polarized ! -1 => One of the selected lines is Zeeman split
                           ! +1 => None of the selected lines is Zeeman split
      integer :: S, SX     ! Indices for sidebands
      logical :: SawNoLines ! Saw a species without lines or continuum
      integer :: STAT      ! Status from allocate or deallocate
      integer :: Z         ! Index for fwdModelConf%Signals

      if ( noLinesMsg < 0 ) noLinesMsg = switchDetail(switches, '0sl') ! Done once
      sawNoLines = .false.

      ! Allocate the spectroscopy catalog extract
      c = maxval(fwdModelConf%cat_size)
      allocate ( fwdModelConf%catalog(s1:s2,c), stat=stat )
      addr = 0
      if ( stat == 0 ) then
        if ( size(fwdModelConf%catalog) > 0 ) &
          & addr = transfer(c_loc(fwdModelConf%catalog(s1,1)), addr)
      end if
      call test_allocate ( stat, moduleName, 'fwdModelConf%catalog', &
        & (/s1,1/), (/s2,c/), storage_size(fwdModelConf%catalog) / 8, address=addr )

      ! Work out the spectroscopy we're going to need.
      fwdModelConf%catalog = empty_cat

      do s = s1, s2, 2
        sx = (s + 3) / 2 ! 1 or 2, instead of -1 or 1.
        c = 0
        do b = 1, size(fwdModelConf%beta_group)
          do m = 1, size(fwdModelConf%beta_group(b)%lbl(sx)%molecules)
            c = c + 1
            fwdModelConf%beta_group(b)%lbl(sx)%cat_index(m) = c
            n = fwdModelConf%beta_group(b)%lbl(sx)%molecules(m)
            if ( n > ubound(catalog,1) .or. n < lbound(catalog,1) ) then ! Probably RHi
              call startErrorMessage ( fwdModelConf%where )
              call display_string ( lit_indices(n), &
                & before=' No spectroscopy catalog for ', advance='yes' )
              error = .true.
              cycle
            else if ( catalog(n)%molecule == l_none ) then
              call startErrorMessage ( fwdModelConf%where )
              call display_string ( lit_indices(n), &
                & before=' No spectroscopy catalog for ', advance='yes' )
              error = .true.
              cycle
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
                  associate ( thisLine => lines(catalog(n)%lines(l)) )
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
                  end associate
                end do     ! l End loop over lines
              end if       ! End case where allLinesInCatalog not set

              ! Check we have at least one line for this specie

              l = count(lineFlag /= 0)
              if ( l == 0 .and. all ( fwdModelConf%catalog(s,c)%continuum == 0.0 ) &
                & .and. noLinesMsg > 0 ) then
                sawNoLines = .true.
                call startErrorMessage ( fwdModelConf%where )
                call display_string ( lit_indices(n), &
                  & before='No relevant lines or continuum for ', advance='yes' )
              end if
              call allocate_test ( fwdModelConf%catalog(s,c)%lines, l, &
                & 'fwdModelConf%catalog(?,?)%lines', moduleName )
              call allocate_test ( fwdModelConf%catalog(s,c)%polarized, l, &
                & 'fwdModelConf%catalog(?,?)%polarized', moduleName )
              fwdModelConf%catalog(s,c)%lines = pack ( catalog(n)%lines, lineFlag /= 0 )
              fwdModelConf%catalog(s,c)%polarized = pack ( lineFlag < 0, lineFlag /= 0 )

            else

              ! No lines for this specie.  However, its continuum may still
              ! be valid so don't set it to empty.
              call allocate_test ( fwdModelConf%catalog(s,c)%lines, 0, &
                & 'fwdModelConf%catalog(?,?)%lines(0)', moduleName )
              call allocate_test ( fwdModelConf%catalog(s,c)%polarized, 0, &
                & 'fwdModelConf%catalog(?,?)%polarized(0)', moduleName )
              if ( all(catalog(n)%continuum == 0.0) .and. noLinesMsg >= 0 ) then
                sawNoLines = .true.
                call startErrorMessage ( fwdModelConf%where )
                call display_string ( lit_indices(n), &
                  & before='No lines or continuum for ', advance='yes' )
!               ! DON'T DO THIS! it catches radiometer-dependent species
!               ! such as O3_R1A that intentionally have no lines or continuum
!                 call get_string ( lit_indices(n), moleculeName )
!                 call MLSMessage ( MLSMSG_Error, moduleName, &
!                   & 'No lines or continuum for ' // trim(moleculeName) )
              end if
            end if
          end do ! m Molecules in fwdModelConf
        end do ! b Beta groups
      end do ! s Sidebands

      if ( sawNoLines .and. noLinesMsg >= 0 ) &
        & call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'At least one species has no lines or continuum' )

    end subroutine SpectroscopyCatalogExtract

  end subroutine DeriveFromForwardModelConfig

  ! --------------------------  DestroyForwardModelConfigDatabase  -----
  subroutine DestroyFWMConfigDatabase ( Database, Deep )

    use ALLOCATE_DEALLOCATE, only: TEST_DEALLOCATE
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database
    logical, optional, intent(in) :: DEEP

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Config                   ! Loop counter
    integer :: S                        ! Size in bytes of object to deallocate
    integer :: Status                   ! Flag

    if ( associated(database) ) then
      do config = 1, size(database)
        call DestroyOneForwardModelConfig ( database(config), deep=deep )
      end do

      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, ModuleName, 'Database', s, address=addr )
    end if
  end subroutine DestroyFWMConfigDatabase

  ! ---------------------------------  DestroyForwardModelDerived  -----
  subroutine DestroyForwardModelDerived ( FwdModelConf )
    ! Destroy stuff in FwdModelConf derived for one forward model run

    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    type ( ForwardModelConfig_T ), intent(inout) :: FwdModelConf

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: C, B, Ier, S

    if ( associated(fwdModelConf%channels) ) then
      s = size(fwdModelConf%channels) * storage_size(fwdModelConf%channels) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(fwdModelConf%channels(1)), addr)
      deallocate ( fwdModelConf%channels, stat = ier )
      call test_deallocate ( ier, ModuleName, 'fwdModelConf%channels', s, address=addr )
  ! else
  !   It was already deallocated at the end of FullForwardModel
    end if

    do b = 1, size(fwdModelConf%beta_group)
      do s = 1, 2
        call deallocate_test ( fwdModelConf%beta_group(b)%PFA(s)%data, &
          & 'Beta_group(b)%PFA(s)%data', moduleName )
      end do ! s
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

      s = size(fwdModelConf%catalog) * storage_size(fwdModelConf%catalog) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc( &
        & fwdModelConf%catalog(lbound(fwdModelConf%catalog,1),1)), addr)
      deallocate ( fwdModelConf%catalog, stat=ier )
      call test_deallocate ( ier, ModuleName, 'fwdModelConf%catalog', s, address=addr )
  ! else
  !   It was already deallocated at the end of FullForwardModel
    end if

  end subroutine DestroyForwardModelDerived

  ! ------------------------------------ NullifyForwardModelConfig -----
  subroutine NullifyForwardModelConfig ( IntentionallyNotUsed )
    ! Given a forward model config, nullify all the pointers associated with it
    type ( ForwardModelConfig_T ), intent(out) :: IntentionallyNotUsed

    ! Executable code not needed since => NULL() initializes pointer
    ! components of intent(out) dummy arguments.
  end subroutine NullifyForwardModelConfig

  ! --------------------------------------------- PVMPackFwmConfig -----
  subroutine PVMPackFWMConfig ( Config )
    use PVMIDL, only: PVMIDLPack
    use MorePVM, only: PVMPackLitIndex, PVMpackStringIndex
    use MLSSignals_M, only: PVMPackSignal
    use VGridsDatabase, only: PVMPackVGrid
    ! Dummy arguments
    type ( ForwardModelConfig_T ), intent(in) :: Config
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
      & config%referenceMIF, config%sidebandStart, config%sidebandStop, &
      & config%surfaceTangentIndex /), msg = "Packing fwmConfig integers" )

    ! Now the logical scalars
    call PVMIDLPack ( (/ &
      & config%allLinesForRadiometer, config%allLinesInCatalog, config%anyLBL, &
      & config%anyPFA, config%atmos_der, config%atmos_second_der, config%default_spectroscopy, &
      & config%differentialScan, config%do_1d, config%do_baseline, &
      & config%do_conv, config%do_freq_avg, config%forceFoldedOutput, &
      & config%forceSidebandFraction, config%generateTScat, &
      & config%globalConfig, config%hmag_der, config%htheta_der, &
      & config%hphi_der, config%ignoreHessian, config%incl_cld, &
      & config%isRadianceModel, config%lockBins,  &
      & config%no_magnetic_field, config%polarized, &
      & config%refract, config%scanAverage, config%skipOverlaps, config%spect_Der, &
      & config%switchingMirror, config%temp_Der, config%transformMIFextinction, &
      & config%useTScat /), &
      & msg ="Packing fwmConfig logicals" )

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
        call PVMIDLPack ( config%moleculeSecondDerivatives, msg = "Packing molecule second derivatives" )
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
  subroutine PVMUnpackFWMConfig ( Config )
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use MLSMessageModule, only: MLSMessage, MLSMsg_Error
    use MLSSignals_M, only: PVMUnpackSignal
    use MorePVM, only: PVMUnpackLitIndex, PVMUnpackStringIndex
    use PVMIDL, only: PVMIDLUnpack
    use VGridsDatabase, only: PVMUnpackVGrid
    ! Dummy arguments
    type ( ForwardModelConfig_T ), intent(out) :: Config
    ! Local variables
    integer, parameter     :: ISMax = 11 ! Number of integers to unpack
    integer, parameter     :: LSMax = 32 ! Number of logicals to unpack
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: INFO                     ! Flag from PVM
    logical :: FLAG                     ! A flag from the sender
    integer, dimension(ISMax) :: IS     ! Temporary array, for integer scalars
    logical, dimension(LSMax) :: LS     ! Temporary array, for logical scalars
    real(r8), dimension(3) :: RS        ! Temporary array, for reals
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

    ! Now the integer scalars. Array IS has to be long enough for this.
    ! If you add any items here, don't forget to make IS longer at the top
    ! of the program unit.
    call PVMIDLUnpack ( is, msg = "Unpacking fwmConfig integers" )
    i = 1
    config%linearsideband         = is(i) ; i = i + 1
    config%no_cloud_species       = is(i) ; i = i + 1
    config%no_model_surfs         = is(i) ; i = i + 1
    config%num_ab_terms           = is(i) ; i = i + 1
    config%num_azimuth_angles     = is(i) ; i = i + 1
    config%num_scattering_angles  = is(i) ; i = i + 1
    config%num_size_bins          = is(i) ; i = i + 1
    config%referenceMIF           = is(i) ; i = i + 1
    config%sideBandStart          = is(i) ; i = i + 1
    config%sideBandStop           = is(i) ; i = i + 1
    config%surfaceTangentIndex    = is(i) !  ; i = i + 1
    if ( i > LSMAX ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // 'PVMUnpackFWMConfig', &
      & 'programming error--ISMAX too small' )

    ! Now the logical scalars. Array LS has to be long enough for this.
    ! If you add any items here, don't forget to make LS longer at the top
    ! of the program unit.
    call PVMIDLUnpack ( ls, msg = "Unpacking fwmConfig logicals" )
    i = 1
    config%allLinesForRadiometer  = ls(i) ; i = i + 1
    config%allLinesInCatalog      = ls(i) ; i = i + 1
    config%anyLBL                 = ls(i:i+1) ; i = i + 2
    config%anyPFA                 = ls(i:i+1) ; i = i + 2
    config%atmos_der              = ls(i) ; i = i + 1
    config%atmos_second_der       = ls(i) ; i = i + 1
    config%default_spectroscopy   = ls(i) ; i = i + 1
    config%differentialScan       = ls(i) ; i = i + 1
    config%do_1d                  = ls(i) ; i = i + 1
    config%do_baseline            = ls(i) ; i = i + 1
    config%do_conv                = ls(i) ; i = i + 1
    config%do_freq_avg            = ls(i) ; i = i + 1
    config%forceFoldedOutput      = ls(i) ; i = i + 1
    config%forceSidebandFraction  = ls(i) ; i = i + 1
    config%globalConfig           = ls(i) ; i = i + 1
    config%generateTScat          = ls(i) ; i = i + 1
    config%hmag_der               = ls(i) ; i = i + 1
    config%htheta_der             = ls(i) ; i = i + 1
    config%hphi_der               = ls(i) ; i = i + 1
    config%ignoreHessian          = ls(i) ; i = i + 1
    config%incl_cld               = ls(i) ; i = i + 1
    config%isRadianceModel        = ls(i) ; i = i + 1
    config%lockBins               = ls(i) ; i = i + 1
    config%no_magnetic_field      = ls(i) ; i = i + 1
    config%polarized              = ls(i) ; i = i + 1
    config%refract                = ls(i) ; i = i + 1
    config%scanAverage            = ls(i) ; i = i + 1
    config%skipOverlaps           = ls(i) ; i = i + 1
    config%spect_der              = ls(i) ; i = i + 1
    config%switchingMirror        = ls(i) ; i = i + 1
    config%temp_der               = ls(i) ; i = i + 1
    config%transformMIFextinction = ls(i) ; i = i + 1
    config%useTScat               = ls(i) !  ; i = i + 1
    if ( i > LSMAX ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // 'PVMUnpackFWMConfig', &
      & 'programming error--LSMAX too small' )

    ! Now the real scalars
    call PVMIDLUnpack ( rs, msg = "Unpacking fwmConfig reals" )
    i = 1
    config%phiWindow(1) = rs(i) ; i = i + 1
    config%phiWindow(2) = rs(i) ; i = i + 1
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
      call Allocate_test ( config%moleculeSecondDerivatives, &
        & n, 'config%moleculeSecondDerivatives', ModuleName )
      do i = 1, n - 1
        call PVMUnpackLitIndex ( config%molecules(i), msg = "Unpacking a molecule" )
        call PVMIDLUnpack ( flag, msg = "Unpacking a molecule sign flag" )
        if ( .not. flag ) config%molecules(i) = - config%molecules(i)
      end do
      config%molecules(n) = huge(config%molecules(n)) ! Sentinel
      call PVMIDLUnpack ( config%moleculeDerivatives, msg = "Unpacking moleculeDerivatives" )
      call PVMIDLUnpack ( config%moleculeSecondDerivatives, msg = "Unpacking moleculeSecondDerivatives" )
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
      addr = 0
      if ( info == 0 ) addr = transfer(c_loc(config%signals(1)), addr)
      call test_allocate ( info, ModuleName, 'config%signals', ubounds = [n], &
        & elementSize = storage_size(config%signals) / 8, address=addr )
      do i = 1, n
        call PVMUnpackSignal ( config%signals(i) )
      end do
    end if

    ! Vgrids
    call PVMIDLUnpack ( ls(1:2), msg = "Unpacking vGrid flags" )
    if ( ls(1) ) then
      allocate ( config%integrationGrid, STAT=info )
      addr = 0
      if ( info == 0 ) addr = transfer(c_loc(config%integrationGrid), addr)
      call test_allocate ( info, ModuleName, 'config%integrationGrid', &
        & uBounds = [1], elementSize = storage_size(config%integrationGrid) / 8, &
        & address=addr )
      call PVMUnpackVGrid ( config%integrationGrid )
    end if
    if ( ls(2) ) then
      allocate ( config%tangentGrid, STAT=info )
      addr = 0
      if ( info == 0 ) addr = transfer(c_loc(config%tangentGrid), addr)
      call test_allocate ( info, ModuleName, 'config%tangentGrid', &
        & uBounds = [1], elementSize = storage_size(config%tangentGrid) / 8, &
        & address=addr )
      call PVMUnpackVGrid ( config%tangentGrid )
    end if

  end subroutine PVMUnpackFWMConfig

  ! --------------------------  StripForwardModelConfigDatabase --------
  subroutine StripForwardModelConfigDatabase ( Database )
    ! This routine removes the non-global forward model configs from the database
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: TmpDatabase
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Config                   ! Loop counter
    integer :: N
    integer :: Status

    ! Executable code
    ! Clear out dying configs
    if ( .not. associated ( database ) ) return
    do config = 1, size ( database )
      if ( .not. database(config)%globalConfig ) &
        & call DestroyOneForwardModelConfig ( database(config) )
    end do

    ! Create new database in tmp space, pack old one into
    n = count ( database%globalConfig )
    allocate ( tmpDatabase ( n ), STAT=status )
    addr = 0
    if ( status == 0 .and. n > 0 ) addr = transfer(c_loc(tmpDatabase(1)), addr)
    call test_allocate ( status, ModuleName, 'tmpDatabase', &
      & uBounds = [n], elementSize = storage_size(tmpDatabase) / 8, address=addr )
    tmpDatabase = pack ( database, database%globalConfig )

    ! Destroy old database, then point to new one
    n = size(database) * storage_size(database) / 8
    addr = 0
    if ( n > 0 ) addr = transfer(c_loc(database(1)), addr)
    deallocate ( database, STAT=status )
    call test_deallocate ( status, ModuleName, 'database', n, address=addr )

    database => tmpDatabase
  end subroutine StripForwardModelConfigDatabase

  ! =====     Private Procedures     =====================================

  ! ------------------------------------ DestroyOneForwardModelConfig --
  subroutine DestroyOneForwardModelConfig ( Config, Deep )
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use MLSSignals_M, only: DestroySignalDatabase

    ! Dummy arguments
    type ( ForwardModelConfig_T), intent(inout) :: config
    logical, optional, intent(in) :: Deep ! Do a really deep destroy

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: B                        ! Subscript for Beta_Group
    integer :: S                        ! 1 = LSB, 2 = USB, or
                                        ! size in bytes of an object to deallocate
    integer :: Status                   ! Flag from allocate etc.
    logical :: MyDeep                   ! Copy of deep

    ! Executable code
    myDeep = .false.
    if ( present ( deep ) ) myDeep = deep
    call destroyForwardModelDerived ( config )
    ! Destroy the beta groups
    do b = 1, size(config%beta_group)
      do s = 1, 2
        call deallocate_test ( config%beta_group(b)%LBL(s)%cat_index, 'Cat_Index', moduleName )
        call deallocate_test ( config%beta_group(b)%LBL(s)%molecules, 'LBL Molecules', moduleName )
        call deallocate_test ( config%beta_group(b)%LBL(s)%ratio, 'LBL Ratio', moduleName )
        call deallocate_test ( config%beta_group(b)%PFA(s)%molecules, 'PFA molecules', moduleName )
        call deallocate_test ( config%beta_group(b)%PFA(s)%ratio, 'PFA ratio', moduleName )
      end do ! s
    end do ! b = 1, size(config%beta_group)

    if ( associated(config%lineCenter) ) then
      s = size(config%lineCenter) * storage_size(config%lineCenter) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(config%lineCenter(1)), addr)
      deallocate ( config%lineCenter, stat=status )
      call test_deallocate ( status, moduleName, 'LineCenter', s, address=addr )
    end if
    if ( associated(config%lineWidth) ) then
      s = size(config%lineWidth) * storage_size(config%lineWidth) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(config%lineWidth(1)), addr)
      deallocate ( config%lineWidth, stat=status )
      call test_deallocate ( status, moduleName, 'LineWidth', s, address=addr )
    end if
    if ( associated(config%lineWidth_TDep) ) then
      s = size(config%lineWidth_TDep) * storage_size(config%lineWidth_TDep) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(config%lineWidth_TDep(1)), addr)
      deallocate ( config%lineWidth_TDep, stat=status )
      call test_deallocate ( status, moduleName, 'LineWidth_TDep', s, address=addr )
    end if

    if ( associated(config%beta_group) ) then
      s = size(config%beta_group) * storage_size(config%beta_group) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(config%beta_group(1)), addr)
      deallocate ( config%beta_group, stat=status )
      call test_deallocate ( status, moduleName, 'config%Beta_group', s, address=addr )
    end if

    if ( associated(config%signals) ) &
      & call DestroySignalDatabase ( config%signals, justChannels=.not. myDeep )
    call deallocate_test ( config%signalIndices, 'SignalIndices', moduleName )
    if ( myDeep ) then
      if ( associated ( config%integrationGrid ) ) then
        call DestroyVGridContents ( config%integrationGrid )
        s = storage_size(config%integrationGrid) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(config%integrationGrid), addr)
        deallocate ( config%integrationGrid, stat=status )
        call test_deallocate ( status, moduleName, 'Config%integrationGrid', s, &
          & address=addr )
      end if
      if ( associated ( config%tangentGrid ) ) then
        call DestroyVGridContents ( config%tangentGrid )
        s = storage_size(config%tangentGrid) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(config%tangentGrid), addr)
        deallocate ( config%tangentGrid, stat=status )
        call test_deallocate ( status, moduleName, 'Config%tangentGrid', s, &
          & address=addr )
      end if
    end if
    ! Otherwise don't destroy integrationGrid and tangentGrid.  Assume they
    ! will be (or already are) destroyed by destroyVGridDatabase.
    ! Don't deallocate config%molecules because it's a pointer into beta_group
    call deallocate_test ( config%specificQuantities, &
      & "config%specificQuantities", ModuleName )
    call deallocate_test ( config%binSelectors, &
      & "config%binSelectors", ModuleName )
  end subroutine DestroyOneForwardModelConfig

  ! --------------------------------------------  Dump_Beta_Group  -----
  subroutine Dump_Beta_Group ( Beta_Group, Name, Sidebands, Details )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use Output_M, only: Blanks, NewLine, Output
    use PFADatabase_M, only: Dump, PFAData
    use String_Table, only: Display_String

    type(beta_group_t), intent(in) :: Beta_Group(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Sidebands(2)
    integer, intent(in), optional :: Details ! for Dump_PFADatum (default 0)

    integer :: B, C, M, S, S1, S2
    logical :: Missing ! PFA data
    integer :: MyDetails
    character(3), parameter :: SB(2) = (/ 'Low', 'Upp' /)
    character(*), parameter :: ParamNames(lineCenter:lineWidth_tDep) = &
      & (/ ' LineCenter=    ', &
      &    ' LineWidth=     ', &
      &    ' LineWidth_tDep=' /)

    s1 = 1; s2 = 2
    if ( present(sidebands) ) then
      s1 = (sidebands(1)+3)/2
      s2 = (sidebands(2)+3)/2
    end if
    myDetails = 0
    if ( present(details) ) myDetails = details
    call output ( '  Beta groups' )
    if ( present(name) ) call output ( ' ' // trim(name) )
    call output ( size(beta_group), before=', SIZE = ', advance='yes' )
    do b = 1, size(beta_group)
      call output ( b, before='  Beta group ', after=': ' )
      call display_string ( lit_indices(beta_group(b)%molecule) )
      if ( beta_group(b)%derivatives ) call output ( ' with derivative' )
      if ( beta_group(b)%secondDerivatives ) call output ( ' with second derivative' )
      call newLine
      do s = s1, s2
        call output ( '  ' // sb(s) // 'er sideband:', advance='yes' )
        if ( size(beta_group(b)%lbl(s)%molecules) > 0 ) then
          call display_string ( &
            & lit_indices(beta_group(b)%lbl(s)%molecules), before='   LBL:', &
            & advance='yes' )
          call dump ( beta_group(b)%lbl(s)%ratio, name='    Isotope ratio' )
          call dump ( beta_group(b)%lbl(s)%cat_index, &
            & name='    Spectroscopy catalog extract index' )
          if ( any(beta_group(b)%lbl(s)%spect_der_ix /= 0) ) then
            call output ( '    Spectrosopy parameter derivatives:' )
            do m = lineCenter, lineWidth_tDep
              if ( beta_group(b)%lbl(s)%spect_der_ix(m) /= 0 ) &
                & call output ( beta_group(b)%lbl(s)%spect_der_ix(m), &
                  & before=trim(paramNames(m)) )
            end do
            call newLine
          end if
        end if
        if ( size(beta_group(b)%pfa(s)%molecules) > 0 ) then
          call display_string ( &
            & lit_indices(beta_group(b)%pfa(s)%molecules), before='   PFA:', &
            & advance='yes' )
          call dump ( beta_group(b)%pfa(s)%ratio, name='    Isotope ratio' )
          if ( associated(beta_group(b)%pfa(s)%data) ) then
            do m = 1, size(beta_group(b)%pfa(s)%data,2)
              missing = .false.
              do c = 1, size(beta_group(b)%pfa(s)%data,1)
                if ( beta_group(b)%pfa(s)%data(c,m) /= 0 ) then
                  call blanks ( 4 )
                  call dump ( PFAData(beta_group(b)%pfa(s)%data(c,m)), details=myDetails )
                else
                  missing = .true.
                end if
              end do ! c
              if ( missing ) call display_string ( &
                & lit_indices(beta_group(b)%pfa(s)%molecules(m)), &
                & before='    Some PFA data missing for ', advance='yes' )
            end do ! m
          end if
        end if
      end do ! s
      call Dump_Qty_Stuff ( beta_group(b)%qty )
    end do
  end subroutine Dump_Beta_Group

  ! ----------------------------  Dump_ForwardModelConfigDatabase  -----
  subroutine Dump_ForwardModelConfigDatabase ( Database, &
    & Where, Details, SkipPFA, QuantityTemplatesDB )

    use MoreTree, only: StartErrorMessage
    use Output_M, only: Output
    use QuantityTemplates, only: QuantityTemplate_T

    type (forwardModelConfig_T), pointer, dimension(:) :: Database
    integer, optional, intent(in) :: Where ! Tree node index
    integer, intent(in), optional :: Details ! for Dump_Beta_Group
    logical, optional, intent(in) :: SkipPFA
    type (quantityTemplate_t), pointer, optional, dimension(:) :: QuantityTemplatesDB

    ! Local variables
    integer :: I                         ! Loop counters

    ! executable code
    if ( associated(database) ) then
      do i = 1, size(database)
        call Dump_ForwardModelConfig( database(i), details=details, &
          & skipPFA=skipPFA, quantityTemplatesDB=quantityTemplatesDB )
      end do
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No forward model database to dump.', advance='yes' )
    end if
  end subroutine Dump_ForwardModelConfigDatabase

  ! -----------------------------------  Dump_ForwardModelConfig  -----
  subroutine Dump_ForwardModelConfig ( Config, Where, &
    & Details, SkipPFA, QuantityTemplatesDB )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices, PHYQ_Indices
    use Lexer_Core, only: Print_Source
    use MLSSignals_M, only: GetNameOfSignal, MaxSigLen, Modules
    use Output_M, only: NewLine, Output
    use QuantityTemplates, only: QuantityTemplate_T
    use String_Table, only: Display_String
    use Tree, only: Null_Tree, Where_At => Where

    type (forwardModelConfig_T), intent(in) :: Config
    character(len=*), intent(in), optional :: Where
    integer, intent(in), optional :: Details ! for Dump_Beta_Group
    logical, optional, intent(in) :: SkipPFA
    type (quantityTemplate_t), pointer, optional, dimension(:) :: QuantityTemplatesDB

    ! Local variables
    logical :: dumpPFA
    integer :: J, S                          ! Loop counters
    integer :: MyDetails
    integer :: S1, S2                        ! Sideband limits
    character (len=MaxSigLen) :: SignalName  ! A line of text

    ! executable code
    dumpPFA = .true.
    if ( present(skipPFA) ) dumpPFA = .not. skipPFA
    myDetails = 0
    if ( present(details) ) myDetails = details

    s1 = (config%sidebandStart+3)/2; s2 = (config%sidebandStop+3)/2
    call display_string ( config%name, before='  Forward Model Config Name: ' )
    if ( config%where /= null_tree ) then
      call output ( ' defined at ' )
      call print_source ( where_at(config%where) )
    end if
    if ( present(where) ) then
      call output ( ' from ' )
      call output ( where )
    end if
    call newLine
    if ( myDetails < 0 ) return
    ! Logical scalars
    call output ( config%allLinesForRadiometer, before='  AllLinesForRadiometer: ', advance='yes' )
    call output ( config%allLinesInCatalog, before='  AllLinesInCatalog: ', advance='yes' )
    call output ( config%anyLBL(1), before='  AnyLBL: ' )
    call output ( config%anyLBL(2) )
    call output ( config%anyPFA(1), before='  AnyPFA: ' )
    call output ( config%anyPFA(2), advance='yes' )
    call output ( config%atmos_der, before='  Atmos_der: ', advance='yes' )
    call output ( config%atmos_second_der, before='  Atmos_second_der: ', advance='yes' )
    call output ( config%default_spectroscopy, before='  Default_spectroscopy: ', advance='yes' )
    call output ( config%DifferentialScan, before='  DifferentialScan: ', advance='yes' )
    call output ( config%do_1d, before='  Do_1D: ', advance='yes' )
    call output ( config%do_Baseline, before='  Do_Baseline: ', advance='yes' )
    call output ( config%do_conv, before='  Do_conv: ', advance='yes' )
    call output ( config%do_freq_avg, before='  Do_freq_avg: ', advance='yes' )
    call output ( config%do_path_norm, before='  Do_path_norm: ', advance='yes' )
    call output ( config%forceFoldedOutput, before='  ForceFoldedOutput: ', advance='yes' )
    call output ( config%forceSidebandFraction, before='  ForceSidebandFraction: ', advance='yes' )
    call output ( config%globalConfig, before='  GlobalConfig: ', advance='yes' )
    CALL output ( config%hmag_der, before='  Hmag_der: ', advance='yes')
    CALL output ( config%htheta_der, before='  Htheta_der: ', advance='yes')
    CALL output ( config%hphi_der, before='  Hphi_der: ', advance='yes')
    call output ( config%incl_cld, before='  Incl_Cld: ', advance='yes' )
    CALL output ( config%lockBins, before='  LockBins: ', advance='yes' )
    call output ( config%no_magnetic_field, before='  No_Magnetic_Field: ', advance='yes' )
    call output ( config%polarized, before='  Polarized: ', advance='yes' )
    call output ( config%refract, before='  Refract: ', advance='yes' )
    call output ( config%scanAverage, before='  ScanAverage: ', advance='yes' )
    call output ( config%skipOverlaps, before='  SkipOverlaps: ',advance='yes' )
    call output ( config%spect_der, before='  Spect_der: ', advance='yes' )
    call output ( config%switchingMirror, before='  SwitchingMirror: ', advance='yes' )
    call output ( config%temp_der, before='  Temp_der: ', advance='yes' )
    call output ( config%transformMIFextinction, before='  TransformMIFextinction: ', advance='yes' )
    ! Strings
    call display_string ( lit_indices(config%cloud_der), before='  Cloud_der: ', advance='yes' )
    call display_string ( lit_indices(config%fwmType), before='  FwmType: ', advance='yes' )
    call display_string ( lit_indices(config%i_saturation), before='  I_saturation: ', advance='yes' )
    call output ( config%instrumentModule, before='  InstrumentModule: ' )
    if ( config%instrumentModule /= 0 ) &
      & call display_string ( modules(config%instrumentModule)%name, &
        & before=' - ' )
    call newline
    call display_string ( lit_indices(config%MIFTangent), before=' MIFTangent: ', advance='yes' )
    ! Integer scalars
    call output ( config%LinearSideband, before='  LinearSideband: ', advance='yes' )
    call output ( config%No_cloud_species, before='  No_cloud_species: ', advance='yes' )
    call output ( config%No_model_surfs, before='  No_model_surfs: ', advance='yes' )
    call output ( config%NoUsedChannels, before='  NoUsedChannels: ', advance='yes' )
    call output ( config%Num_ab_terms, before='  Num_ab_terms: ', advance='yes' )
    call output ( config%Num_azimuth_angles, before='  Num_azimuth_angles: ', advance='yes' )
    call output ( config%Num_scattering_angles, before='  Num_scattering_angles: ', advance='yes' )
    call output ( config%Num_size_bins, before='  Num_size_bins: ', advance='yes' )
    call output ( config%ReferenceMIF, before='  ReferenceMIF: ', advance='yes' )
    call output ( config%sidebandStart, before='  Sidebands: ' )
    call output ( config%sidebandStop, before=' ', advance='yes' )
    call output ( config%SurfaceTangentIndex, before='  SurfaceTangentIndex: ', advance='yes' )
    call output ( config%TScatMIF, before='  TScatMIF: ', advance='yes' )
    call output ( config%xStar, before='  xStar: ', advance='yes' )
    call output ( config%yStar, before='  yStar: ', advance='yes' )
    ! Integer arrays
    if ( associated(config%specificQuantities) ) then
      if ( present(quantityTemplatesDB) ) then
        call output ( '  Specific quantities:', advance='yes' )
        do j = 1, size(config%specificQuantities)
          call display_string ( quantityTemplatesDB(config%specificQuantities(j))%name, &
            before='    ', advance='yes' )
        end do
      else
        call dump ( config%specificQuantities, name='  SpecificQuantities' )
      end if
    end if
    ! Real scalars
    call output ( config%PhiWindow(1), before='  PhiWindow: ' )
    call output ( config%PhiWindow(2), before=' ' )
    call display_string ( phyq_indices(config%windowUnits), before=' ', advance='yes' )
    call output ( config%Tolerance, before='  Tolerance: ', advance='yes' )
    ! Bin selectors
    if ( associated ( config%binSelectors ) ) &
      & call dump ( config%binSelectors, name=  '  BinSelectors: ' )

    if ( associated(config%Beta_group) .and. dumpPFA ) &
      & call dump ( config%Beta_group, &
        & sidebands=(/config%sidebandStart,config%sidebandStop/), details=details )
    call output ( '  Molecules: ', advance='yes' )
    if ( associated(config%molecules) ) then
      do j = 1, size(config%molecules)
        call display_string ( lit_indices(config%molecules(j)), before='    ' )
        call output ( j, before=':' )

        if (config%moleculeDerivatives(j)) then
          call output (' compute derivatives' )
        else
          call output (' no derivatives' )
        end if

        if (config%moleculeSecondDerivatives(j)) then
          call output (' compute second derivatives' )
        else
          call output (' no second derivatives' )
        end if  

        call newline

      end do
    end if
    if ( size(config%lineCenter) > 0 .or. size(config%lineWidth) > 0 .or. &
      &  size(config%lineWidth_Tdep) > 0 ) then
      call output ( '  Spectroscopy parameters:', advance='yes')
      if ( size(config%lineCenter) > 0 ) &
        & call display_string ( lit_indices(config%lineCenter%molecule), &
          & advance='yes', before='    Line centers:' )
      if ( size(config%lineWidth) > 0 ) &
        & call display_string ( lit_indices(config%lineWidth%molecule), &
          & advance='yes', before='    Line widths:' )
      if ( size(config%lineWidth_Tdep) > 0 ) &
        & call display_string ( lit_indices(config%lineWidth%molecule), &
          & advance='yes', before='    Line width temperature dependencies:' )
    end if
    call output ( '  Signals:', advance='yes')
    if ( associated(config%signals) ) then
      do j = 1, size(config%signals)
        if ( config%signals(j)%index < 0 ) then
          call output ( j, before='unable to parse signal ', advance='yes' )
        else
          call getNameOfSignal ( config%signals(j), signalName)
          call output ( '  '//trim(signalName)//' channels Included:' )
          if ( associated(config%signals(j)%channels) ) then
            call newLine
            call dump ( config%signals(j)%channels )
          else
            call output ( ' None', advance='yes' )
          end if
        end if
      end do
    else
      call output ( '  (none associated yet)', advance='yes')
    endif
    ! Dump ForwardModelDerived
    call output ( '  ForwardModelDerived:', advance='yes' )
    if ( associated(config%usedDACSSignals) ) &
      call dump  ( config%usedDACSSignals, name='   Used DACS Signals' )
    if ( associated(config%channels) ) then
      call output ( '   Channel info:', advance='yes' )
      do j = 1, size(config%channels)
        call output ( config%channels(j)%used, before='    Used: ' )
        call output ( config%channels(j)%origin, before='    Origin: ' )
        call output ( config%channels(j)%signal, before='    Signal: ' )
        call output ( config%channels(j)%DACS, before='    DACS: ' )
        call output ( '    Shape Inds:' )
        do s = s1, s2
          call output ( config%channels(j)%shapeInds(s), before=' ' )
        end do
        call newLine
      end do ! j = 1, size(config%channels)
    end if
    if ( associated(config%IntegrationGrid) ) then
      call dump ( config%IntegrationGrid, details=details, what='Integration grid' )
    else
      call output ( '  no IntegrationGrid', advance='yes' )
    end if
    if ( associated(config%tangentGrid) ) then
      call dump ( config%tangentGrid, details=details, what='Tangent grid' )
    else
      call output ( '  no tangentGrid', advance='yes' )
    end if
    if ( associated(config%beta_Group) .and. .not. dumpPFA ) &
      & call dump ( config%beta_Group, &
        & sidebands=(/config%sidebandStart,config%sidebandStop/), &
        & details=details, name='Beta group' )
    if ( .not. associated(config%beta_Group) ) &
      & call output ( ' no Beta group (surprised?)', advance='yes' )
  end subroutine Dump_ForwardModelConfig

  ! ---------------------------------------------  Dump_Qty_Stuff  -----
  subroutine Dump_Qty_Stuff ( Qty, Details )
    use HighOutput, only: OutputNamedValue
    use Output_M, only: NewLine, Output
    use VectorsModule, only: Dump
    type(QtyStuff_t), intent(in) :: Qty
    integer, optional, intent(in) :: Details
    if ( associated(qty%qty) ) then
      call dump ( qty%qty, details=details )
    else
      call output( '  Quantity not assocated for this qtyStuff type', advance='yes' )
    endif
    call outputNamedValue ( '  Found in first', qty%foundInFirst )
    call outputNamedValue ( '  Was Specific', qty%wasSpecific )
    call newLine
  end subroutine Dump_Qty_Stuff

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id: ForwardModelConfig.f90,v 2.145 2019/10/07 20:04:12 vsnyder Exp $"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ForwardModelConfig

! $Log: ForwardModelConfig.f90,v $
! Revision 2.146  2024/08/15 15:50:02  wgread
! Add hmag_der, htheta_der, hphi_der
!
! Revision 2.145  2019/10/07 20:04:12  vsnyder
! Add trapezoid field for quadrature in FullForwardModel
!
! Revision 2.144  2019/04/24 19:17:24  vsnyder
! Add MIFTangent field to forwardModel
!
! Revision 2.143  2018/08/06 19:58:24  vsnyder
! Use ASSOCIATE construct to avoid necessity for Lines database to have the
! TARGET attribute.  Some cannonball polishing.
!
! Revision 2.142  2018/05/15 03:26:25  vsnyder
! Change Mie tables from pointer to allocatable
!
! Revision 2.141  2018/04/11 22:25:23  vsnyder
! Remove USE for unused names
!
! Revision 2.140  2017/09/20 01:07:57  vsnyder
! Add check that filter shapes are provided if frequency averaging is
! selected.
!
! Revision 2.139  2017/09/15 15:45:18  livesey
! Modified to allow omission of filter shapes file under appropriate circumstances
!
! Revision 2.138  2017/08/17 16:29:54  livesey
! Allowed there to be a missing filter shapes file if frequency
! averaging is not being called for
!
! Revision 2.137  2017/02/04 02:17:40  vsnyder
! Undo all-caps switch in USE statements and some declarations.  Add Values1,
! Values (rank 2) and Values3 pointers to QtyStuff_t, but they're not used
! yet, and might never be used.
!
! Revision 2.136  2016/05/02 23:31:50  vsnyder
! Add QtyStuff component for temperature
!
! Revision 2.135  2015/08/25 17:21:47  vsnyder
! PhiWindow is a tuple, with the first element specifying the angles or
! number of profiles/MAFs before the tangent point, and the second
! specifying the angles or number after.
!
! Revision 2.134  2015/03/28 01:59:22  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.133  2014/09/29 20:25:58  vsnyder
! Add No_Magnetic_Field to config
!
! Revision 2.132  2014/09/05 20:48:44  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.131  2014/08/01 01:01:55  vsnyder
! Change noLinesMsg level from 1 to zero
!
! Revision 2.130  2014/01/09 00:26:39  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.129  2013/09/24 23:28:17  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.128  2013/08/16 02:32:14  vsnyder
! Remove ModelPlaneMIF
!
! Revision 2.127  2013/08/12 23:48:08  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.126  2013/08/09 01:02:58  vsnyder
! Add ReferenceMIF component
!
! Revision 2.125  2013/08/08 02:34:54  vsnyder
! Add derivOK component to Qty_Stuff
!
! Revision 2.124  2013/07/25 00:22:09  vsnyder
! Replace TransformRHI with TransformMIFRHI
!
! Revision 2.123  2013/07/19 01:18:45  vsnyder
! Add TransformRHI component
!
! Revision 2.122  2013/07/13 00:07:28  vsnyder
! Add Model_Plane_MIF component
!
! Revision 2.121  2013/05/15 03:08:55  vsnyder
! Revise processing of dump switch
!
! Revision 2.120  2013/03/30 00:11:39  vsnyder
! Add Quantity database to dump, dump more stuff in config, add wasSpecific
! to QtyStuff_t.
!
! Revision 2.119  2013/03/20 22:46:42  vsnyder
! Give default value to QtyStuff_t%foundInFirst
!
! Revision 2.118  2013/03/01 01:08:11  pwagner
! Snip all but name of dumped config if details < 0
!
! Revision 2.117  2012/05/01 22:18:50  vsnyder
! Add IsRadianceModel component
!
! Revision 2.116  2012/03/07 00:45:28  vsnyder
! Add TransformMIFExtinction component
!
! Revision 2.115  2012/01/25 00:07:07  vsnyder
! Use Test_Allocate, Test_Deallocate
!
! Revision 2.114  2011/07/29 01:51:38  vsnyder
! Remove TScatMolecules and TScatMoleculeDerivatives fields.  Make CloudIce
! a molecule.  Look for CloudIce instead of Cloud_A and Cloud_S
!
! Revision 2.113  2011/05/09 17:45:38  pwagner
! Converted to using switchDetail
!
! Revision 2.112  2011/03/31 19:50:29  vsnyder
! Don't dump beta group twice
!
! Revision 2.111  2010/12/06 19:15:42  pwagner
! More detailed dump of config to aid debugging
!
! Revision 2.110  2010/09/25 01:08:39  vsnyder
! Cannonball polishing
!
! Revision 2.109  2010/08/27 06:13:37  yanovsky
! Add atmos_second_der, MoleculeSecondDerivatives.
! Add SecondDerivatives component to Beta_Group_T.
!
! Revision 2.108  2010/06/09 16:33:59  pwagner
! Complain if IS, LS array bounds exceeded
!
! Revision 2.107  2010/06/07 23:20:51  vsnyder
! Add UseTScat, TScatMolecules, TScatMoleculeDerivatives, change name
! of PhaseFrqTol to FrqTol, check that Mie tables are loaded of UseTScat
! is selected.
!
! Revision 2.106  2010/05/14 02:18:16  vsnyder
! Dump BinSelectors
!
! Revision 2.105  2010/03/26 23:12:51  vsnyder
! Add ignoreHessian field for quasi-linear model
!
! Revision 2.104  2010/01/23 01:04:26  vsnyder
! Make sure Mie tables have been read if needed
!
! Revision 2.103  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.102  2008/08/27 19:56:51  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.101  2008/07/31 17:39:34  vsnyder
! Add PhaseFrqTol to the config
!
! Revision 2.100  2008/05/20 00:27:32  vsnyder
! Add GenerateTScat field
!
! Revision 2.99  2007/11/07 03:07:42  vsnyder
! Add Do_Path_Norm switch
!
! Revision 2.98  2007/10/02 22:34:49  vsnyder
! Cannonball polishing
!
! Revision 2.97  2007/06/25 20:34:46  vsnyder
! Replace tabs by spaces, since tabs are nonstandard
!
! Revision 2.96  2006/11/30 23:32:17  vsnyder
! Use SIZE() instead of ASSOCIATED for spectroscopy derivative stuff
!
! Revision 2.95  2006/07/21 00:17:39  vsnyder
! Remove unused declarations and USEs
!
! Revision 2.94  2006/06/03 01:46:10  vsnyder
! Remove no_dup_mol flag from config structure
!
! Revision 2.93  2006/06/01 02:55:21  vsnyder
! Make sure config%lines is associated
!
! Revision 2.92  2006/05/31 22:02:13  vsnyder
! Display 'no lines or continuum' msg if none in the catalog
!
! Revision 2.91  2006/05/11 19:36:14  pwagner
! Added option to disallow duplicate molecules
!
! Revision 2.90  2006/04/25 23:25:36  vsnyder
! Revise DACS filter shape data structure
!
! Revision 2.89  2006/04/11 18:37:21  vsnyder
! Check for DACS channel information
!
! Revision 2.88  2006/02/23 00:58:58  vsnyder
! Don't crash while dumping config if the signal can't be parsed
!
! Revision 2.87  2006/02/08 21:37:36  vsnyder
! Delay halting on error until all possible error messages are emitted.
! Add descriptive titles to integration grid and tangent grid dumps.
!
! Revision 2.86  2006/02/08 01:02:01  vsnyder
! More stuff for spectroscopy derivatives
!
! Revision 2.85  2005/12/29 01:12:04  vsnyder
! Add 'refract' field, polished up dump
!
! Revision 2.84  2005/11/15 00:23:21  pwagner
! Must not attempt to dump config%signals if not associated
!
! Revision 2.83  2005/11/02 21:38:24  vsnyder
! Hoist some stuff out of MAF loop, delete debugging code
!
! Revision 2.82  2005/11/01 23:01:36  vsnyder
! Precompute ShapeInds and stash in config
!
! Revision 2.81  2005/09/17 00:48:42  vsnyder
! Cannonball polishing
!
! Revision 2.80  2005/09/03 01:21:33  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.79  2005/08/19 23:32:06  pwagner
! option to allow skipping voluminous PFA DB dump when dumping FwdMdl
!
! Revision 2.78  2005/08/03 18:04:09  vsnyder
! Some spectroscopy derivative stuff
!
! Revision 2.77  2005/06/03 01:58:53  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! Revise PFA data structures.
!
! Revision 2.76  2005/05/28 03:27:21  vsnyder
! Simplify PFAStuff
!
! Revision 2.75  2005/05/26 20:12:16  vsnyder
! Delete some debugging code
!
! Revision 2.74  2005/05/26 20:11:31  vsnyder
! Don't delete PFA molecules and ratio in DestroyForwardModelDerived
!
! Revision 2.73  2005/05/26 02:15:14  vsnyder
! Use molecule.signal.sideband.channel structure for PFA
!
! Revision 2.72  2005/05/24 01:54:35  vsnyder
! Delete unused symbols
!
! Revision 2.71  2005/05/05 20:48:02  vsnyder
! Don't check IER if deallocate isn't done
!
! Revision 2.70  2005/05/05 01:14:22  vsnyder
! Make sure fields of PFA_t are nullified.
! Don't try to deallocate channels once for each sideband -- there's only one
! Don't try to deallocate fwdModelConf%beta_group(b)%PFA(s)%data if it's not
! associated.
!
! Revision 2.69  2005/05/02 23:04:03  vsnyder
! Stuff for PFA Cacheing
!
! Revision 2.68  2005/03/28 20:27:51  vsnyder
! Lots of PFA stuff
!
! Revision 2.67  2005/02/17 02:35:13  vsnyder
! Remove PFA stuff from Channels part of config
!
! Revision 2.66  2005/02/16 23:16:49  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.65  2004/12/28 00:27:29  vsnyder
! Remove unreferenced use names
!
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
