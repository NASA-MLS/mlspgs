! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullCloudForwardModel

! ===========================================================================
! THIS MODULE CONTAINS THE FULL CLOUD FORWARD MODEL  
! ===========================================================================

  ! This one is here instead of inside FullCloudForwardModelWrapper because
  ! it has a dummy argument of the same name.
  use ForwardModelConfig,         only: FORWARDMODELCONFIG_T   

  implicit none
  private

  public :: FullCloudForwardModelWrapper

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------

         ! ---------------------------------------------------------------------
contains ! THIS SUBPROGRAM CONTAINS THE WRAPPER ROUTINE FOR CALLING THE FULL
         ! CLOUD FORWARD MODEL
         ! ---------------------------------------------------------------------

  !----------------------------------------  FullCloudForwardModelWrapper  -----
  subroutine FullCloudForwardModelWrapper ( ForwardModelConfig, FwdModelIn,  &
                                            FwdModelExtra, FwdModelOut, Ifm, &
                                            fmStat, Jacobian                 )  

    use Allocate_deallocate,        only: Allocate_test, Deallocate_test
    use AntennaPatterns_m,          only: ANTENNAPATTERNS
    use CloudySkyModule,            only: CLOUD_MODEL
    use CloudySkyRadianceModel,     only: CloudForwardModel
    use ForwardModelIntermediate,   only: FORWARDMODELINTERMEDIATE_T, &
                                        & FORWARDMODELSTATUS_T
    use ManipulateVectorQuantities, only: FindClosestInstances
    use MatrixModule_0,             only: MATRIXELEMENT_T, &
                                        & M_FULL, CREATEBLOCK
    use MatrixModule_1,             only: MATRIX_T, FINDBLOCK
    use MLSCommon,                  only: r8, rm, rp
    use MLSMessageModule,           only: MLSMessage, MLSMSG_Allocate, &
                                        & MLSMSG_Error, MLSMSG_Warning, &
                                        & MLSMSG_Deallocate
    use Molecules,                  only: L_H2O, L_N2O, L_O3
    use MLSSignals_m,               only: ARESIGNALSSUPERSET, GetNameOfSignal, &
                                        & GetSidebandStartStop, SIGNAL_T
    use MLSNumerics,                only: InterpolateValues
!   use Output_m,                   only: OUTPUT
    use Physics,                    only: C => SpeedOfLight
    use SpectroscopyCatalog_m,      only: CATALOG, CATALOG_T, EMPTY_CAT, &
                                        & LINE_T, LINES, MostLines
    use String_table,               only: GET_STRING
    use Toggles,                    only: Emit, Levels, Toggle
    use Trace_M,                    only: Trace_begin, Trace_end
    use Units,                      only: Deg2Rad                         
    use VectorsModule,              only: GETVECTORQUANTITYBYTYPE, VECTOR_T, &
                                        & VECTORVALUE_T, VALIDATEVECTORQUANTITY

! ----------------------------------------------------------
! DEFINE INTRINSIC CONSTANTS NEEDED FROM Init_Tables_Module
! ----------------------------------------------------------

    use Intrinsic, only: &
                       & L_CLOUDEXTINCTION,                                    &
                       & L_CLOUDICE,                                           &
                       & L_CLOUDINDUCEDRADIANCE,                               &
                       & L_CLOUDRADSENSITIVITY,                                &
                       & L_CLOUDWATER,                                         &
                       & L_EARTHRADIUS,                                        &
                       & L_EFFECTIVEOPTICALDEPTH,                              &
                       & L_ELEVOFFSET,                                         &
                       & L_GPH,                                                &
                       & L_IWC_HIGH_HEIGHT,                                    &
                       & L_IWC_LOW_HEIGHT,                                     &
                       & L_LIMBSIDEBANDFRACTION,                               &
                       & L_LOSTRANSFUNC,                                       &
                       & L_LOSVEL,                                             &
                       & L_MASSMEANDIAMETERICE,                                & 
                       & L_MASSMEANDIAMETERWATER,                              &
                       & L_NONE,                                               &
                       & L_PTAN,                                               &
                       & L_RADIANCE,                                           &
                       & L_SCGEOCALT,                                          &
                       & L_SIZEDISTRIBUTION,                                   &
                       & L_SURFACETYPE,                                        &
                       & L_TEMPERATURE,                                        &
                       & L_TOTALEXTINCTION,                                    &
                       & L_VMR,                                                &
                       & LIT_INDICES


    ! Dummy arguments
    type(forwardModelConfig_T),       intent(inout) :: FORWARDMODELCONFIG
    type(vector_T),                   intent(in)    :: FWDMODELIN, FwdModelExtra
    type(vector_T),                   intent(inout) :: FWDMODELOUT             ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM                     ! Workspace
    type(forwardModelStatus_t),       intent(inout) :: FMSTAT                  ! Reverse comm. stuff
    type(matrix_T),                   intent(inout), optional :: JACOBIAN


    ! Local parameters ---------------------------------------------------------
    character(len=*), parameter :: INVALIDQUANTITY = &
     & "Invalid vector quantity for "

    ! Local variables
    type (VectorValue_T), pointer :: CLOUDICE                   ! Profiles
    type (VectorValue_T), pointer :: CLOUDWATER                 ! Profiles
    type (VectorValue_T), pointer :: CLOUDEXTINCTION            ! Profiles
    type (VectorValue_T), pointer :: modelCLOUDRADIANCE         ! modelled cloud radiance
!   type (VectorValue_T), pointer :: constrainCldRad            ! observed cloud radiance
    type (VectorValue_T), pointer :: CLOUDRADSENSITIVITY        ! Like radiance
    type (VectorValue_T), pointer :: EFFECTIVEOPTICALDEPTH      ! Quantity
    type (VectorValue_T), pointer :: GPH                        ! Geop height
    type (VectorValue_T), pointer :: MASSMEANDIAMETERICE        ! Quantity
    type (VectorValue_T), pointer :: MASSMEANDIAMETERWATER      ! Quantity
    type (VectorValue_T), pointer :: PTAN                       ! Tgt pressure
    type (VectorValue_T), pointer :: RADIANCE                   ! Quantity
    type (VectorValue_T), pointer :: SIZEDISTRIBUTION           ! Integer really
    type (VectorValue_T), pointer :: EARTHRADIUS                ! Scalar 
    type (VectorValue_T), pointer :: SURFACETYPE                ! Integer really
    type (VectorValue_T), pointer :: TEMP                       ! Temperature 
    type (VectorValue_T), pointer :: TOTALEXTINCTION            ! Profile
    type (VectorValue_T), pointer :: VMR                        ! Quantity
    type (VectorValue_T), pointer :: SCGEOCALT                  ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: ELEVOFFSET                 ! Elevation offset quantity
    type (VectorValue_T), pointer :: LOSVEL                     ! Line of sight velocity
    type (Signal_T)               :: signal                     ! A signal
    type (VectorValue_T), pointer :: STATE_ext                  ! A state vector quantity
    type (VectorValue_T), pointer :: STATE_los                  ! A state vector quantity
    type (VectorValue_T), pointer :: SIDEBANDFRACTION              ! The sideband ratio to use
    type (VectorValue_T), pointer :: LOWERSIDEBANDFRACTION         ! From the state vector
    type (VectorValue_T), pointer :: UPPERSIDEBANDFRACTION         ! From the state vector

    type (catalog_T), dimension(:,:), pointer :: MY_CATALOG 
    type (catalog_T), pointer :: thisCatalogEntry
    type (line_T),    pointer :: thisLine

    ! for jacobian
    type(MatrixElement_T), pointer :: JBLOCK       ! A block from the jacobian
    integer :: COLJBLOCK                           ! Column index in jacobian
    integer :: ROWJBLOCK                           ! Row index in jacobian
    integer :: noInstances                         ! no of instance
    integer :: noMIFs                              ! Number of minor frames
    integer :: noSgrid                             ! no of elements in S grid
    integer :: noSurf                              ! Number of pressure levels
    integer :: novmrSurf                           ! Number of vmr levels
    integer :: noExtSurf                           ! Number of cloud ext levels
    integer :: noIWCSurf                           ! Number of cloud IWC levels
    integer :: noLWCSurf                           ! Number of cloud LWC levels
    integer :: noPSDSurf                           ! Number of size dist. levels
    integer :: NOFREQS                             ! Number of frequencies
    integer :: DIRECTION                           ! Direction of channel numbering

    integer :: i                                   ! Loop counter
    integer :: j                                   ! Loop counter
    integer :: k                                   ! Loop counter
    integer :: IER                                 ! Status flag from allocates
    integer :: mif
    integer :: MAF                                 ! major frame counter
    integer :: INSTANCE                            ! Relevant instance for temperature
    integer :: MinInst                             ! lower bound of instance
    integer :: MaxInst                             ! upper bound of instance
    integer :: nfine                               ! no of fine resolution grids
    integer :: nNear                               ! no of nearest profiles
    integer :: status                              ! allocation status 
    integer :: SIDEBANDSTART                       ! For sideband loop
    integer :: SIDEBANDSTOP                        ! For sideband loop
    integer :: THISSIDEBAND                        ! Loop counter for sidebands
    integer :: SIGIND                              ! Signal index, loop counter
    integer :: nspec                               ! no of species in cloud fw model
    integer :: ispec                               ! species index in cloud fw model

    integer :: iCloudHeight                        ! Index for Cloud Top Height

    integer, dimension(:), pointer :: closestInstances 

    integer :: WHICHChannel                        ! which single channel is used
    integer :: WHICHPATTERN                        ! Index of antenna pattern
    integer :: MAXSUPERSET                         ! Max. value of superset
    integer, dimension(:), pointer :: SUPERSET     ! Result of AreSignalsSuperset
    integer, dimension(1) :: WHICHPATTERNASARRAY   ! Result of minloc

    real(r8), dimension(:,:), pointer :: A_CLEARSKYRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDINDUCEDRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDEXTINCTION
    real(r8), dimension(:,:), pointer :: A_CLOUDRADSENSITIVITY
    real(r8), dimension(:,:), pointer :: A_EFFECTIVEOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_MASSMEANDIAMETER
    real(r8), dimension(:,:), pointer :: A_TOTALEXTINCTION
    real(r8), dimension(:,:), pointer :: VMRARRAY  ! The VMRs
    real(r8), dimension(:), allocatable :: Slevl   ! S grid
    real(r8), dimension(:), allocatable :: Zt      ! tangent height
    real(r8), dimension(:), pointer :: thisFraction ! Sideband ratio values
    real(rp) :: Vel_Cor                 ! Velocity correction due to Vel_z

    real(r8), dimension(:), pointer :: phi_fine    ! Fine resolution for phi 
    real(r8), dimension(:), pointer :: z_fine      ! Fine resolution for z
    real(r8), dimension(:), pointer :: zp_fine     ! Fine resolution for zp
    real(r8), dimension(:), pointer :: s_fine      ! Fine resolution for s
    real(r8), dimension(:), pointer :: ds_fine     ! Fine resolution for ds
    real(r8), dimension(:), pointer :: w_fine      ! weight along s_fine

    real(r8), dimension(:,:,:), allocatable  :: A_TRANS
    real(r8), dimension(:), pointer :: FREQUENCIES 
    real(r8), dimension(:,:), allocatable :: WC
    real(r8), dimension(:), allocatable :: PSD

    real(r8) :: dz                                 ! thickness of state quantity
    real(r8) :: dphi                               ! phi interval of state quantity
    real(r8) :: tLat                               ! temperature 'window' latitude
    real(r8) :: CloudHeight                        ! Cloud Top Height

    logical, dimension(:), pointer :: doChannel    ! Do this channel?
    logical :: DoThis                              ! Flag for lines
    logical :: FOUNDINFIRST                        ! Flag to indicate derivatives
    logical, dimension(:), pointer :: LINEFLAG     ! Use this line (noLines per species)
    logical, target :: MaxLineFlag(mostLines)
    logical :: prt_log = .false.

    character :: Cloudtype                         ! cloud profile type
    character (len=32) :: MolName                  ! Name of a molecule
    character (len=32) :: SignalName               ! Name of a signal

    !---------------------------------------------------------------------------
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>> Executable code  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !---------------------------------------------------------------------------

    if ( toggle(emit) ) call trace_begin ( 'FullCloudForwardModel' )

    !--------------------------
    ! Nullify all the pointers
    !--------------------------

    nullify( a_clearSkyradiance, a_cloudExtinction, a_cloudInducedradiance,   &
             a_cloudRadsensitivity, a_effectiveOpticaldepth,                  &
             a_massMeandiameter, a_totalExtinction, cloudExtinction,          &
             cloudIce, cloudRadsensitivity, closestInstances, cloudWater,     &
             doChannel, earthradius, effectiveOpticaldepth, frequencies, gph, &
             jblock, losvel, massMeandiameterice,                             &
             massMeandiameterwater, modelCloudradiance, my_catalog,           &
!            constrainCldRad, &
             & ptan,  radiance, sizeDistribution, state_ext,                  &
             state_los, superset, surfaceType, temp, thisLine, thisFraction,  &
             totalExtinction, vmr, vmrarray )
 
    ! ------------------------------------
    ! Find which maf is called at present
    ! ------------------------------------
    maf = fmStat%maf

    ! --------------------------------------------
    ! Get the quantities we need and don't depend on
    ! the signal from the vectors
    ! --------------------------------------------

    ! --------
    ! Outputs:
    ! --------
    cloudExtinction => GetVectorQuantityByType ( fwdModelOut,            & 
        & quantityType=l_cloudExtinction, noerror=.true.)
    totalExtinction => GetVectorQuantityByType ( fwdModelOut,            &
      & quantityType=l_totalExtinction, noerror=.true.)
    massMeanDiameterIce => GetVectorQuantityByType ( fwdModelOut,        &
      & quantityType=l_massMeanDiameterIce, noerror=.true. )
    massMeanDiameterWater => GetVectorQuantityByType ( fwdModelOut,      &
      & quantityType=l_massMeanDiameterWater, noerror=.true. )

    ! -------
    ! Inputs:
    ! -------
    temp => GetVectorQuantityByType ( fwdModelExtra,                     &
      & quantityType=l_temperature )
    gph => GetVectorQuantityByType ( fwdModelExtra,                      &
      & quantityType=l_gph )
    cloudIce => GetVectorQuantityByType ( fwdModelExtra,                 &
      & quantityType=l_cloudIce )
    cloudWater => GetVectorQuantityByType ( fwdModelExtra,               &
      & quantityType=l_cloudWater )
    surfaceType => GetVectorQuantityByType ( fwdModelExtra,              &
      & quantityType=l_surfaceType )
    sizeDistribution=>GetVectorQuantityByType(fwdModelExtra,             &
      & quantityType=l_sizeDistribution )
    earthradius=>GetVectorQuantityByType ( fwdModelExtra,                &
      & quantityType=l_earthradius ) 
    scGeocAlt => GetVectorQuantityByType ( fwdModelExtra,                &
      & quantityType=l_scGeocAlt )

    !-----------------------------------------
    ! Make sure the quantities we need are got and with correct format
    !-----------------------------------------
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(gph, stacked=.true., coherent=.true.,  &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )

    ! ----------------------------
    ! Get some basic dimensions
    ! ----------------------------
    noSurf  = temp%template%noSurfs     ! Number of model layers
    if ( noSurf /= GPH%template%nosurfs ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'number of levels in gph does not match no of levels in temp' )

    !----------------------------
    ! get state quantity type (need both ext and los quantities for Jacobians)
    !----------------------------
    state_ext => GetVectorQuantityByType(FwdModelIn,quantityType=l_cloudExtinction)
    state_los => GetVectorQuantityByType(FwdModelIn,quantityType=l_LosTransFunc)

    ! Get s dimension
    noSgrid = state_los%template%noChans

    allocate ( Slevl(noSgrid), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'Slevl' )
    Slevl = state_los%template%frequencies

    ! Get number of cloud surfaces for retrieval
    noPSDSurf=sizeDistribution%template%noSurfs
    noIWCSurf=cloudice%template%noSurfs
    noLWCSurf=cloudwater%template%noSurfs
    noExtSurf=state_ext%template%noSurfs

    ! use Temp vGrid for all cloud variables
    allocate ( WC(ForwardModelConfig%no_cloud_species,NoSurf), STAT=status )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'WC' )
    allocate ( PSD(NoSurf), STAT=status )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'PSD' )

    ! check if these dimensions are same
    if ( noIWCSurf /= noLWCSurf .or. noIWCSurf /= noExtSurf ) then
      call MLSMessage(MLSMSG_Error, ModuleName, &
      &'IWC, LWC and Extinction profiles have different vertical grids')
    end if

    nspec = size(forwardModelConfig%molecules)
    if ( nspec < 2 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Not enough molecules' )

    call allocate_test ( vmrArray, nspec, noSurf, 'vmrArray', ModuleName )
    vmrarray = 0._r8

    allocate ( My_Catalog(-1:1,nspec), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )   

    call allocate_test ( superset, size(antennaPatterns), &
      & 'superset', ModuleName )
    call Allocate_test ( a_massMeanDiameter, 2, noSurf, &
      & 'a_massMeanDiameter', ModuleName )

    !--------------------------------------------
    ! Loop over signals
    !--------------------------------------------
    do sigInd = 1, size(forwardModelConfig%signals)

    ! -------------------------------------
    ! Identify the signal (band)
    ! -------------------------------------

      signal = forwardModelConfig%signals(sigInd)

      if ( prt_log ) print*,'signal%index', signal%index

      ! find which channel is used
      noFreqs = size (signal%frequencies)
      do i=1, noFreqs
        if ( signal%channels(i) ) whichChannel=i    
      end do

      !-----------------------
      ! Think about sidebands
      !-----------------------
      call getSidebandStartStop ( signal, sidebandStart, sidebandStop ) 

      ! --------------------------------------------
      ! Get the quantities we need and that do depend
      ! on the signal from the vectors
      ! --------------------------------------------

      ! --------
      ! Outputs:
      ! --------
      radiance => GetVectorQuantityByType ( fwdModelOut,                   &
        & quantityType=l_radiance,                                         &
        & signal=signal%index, sideband=signal%sideband )
      modelCloudRadiance => GetVectorQuantityByType ( fwdModelOut,         &
        & quantityType=l_cloudInducedRadiance,                             &
        & signal=signal%index, sideband=signal%sideband )
      cloudRADSensitivity => GetVectorQuantityByType ( fwdModelOut,        &
        & quantityType=l_cloudRADSensitivity, noerror=.true.,              &
        & signal=signal%index, sideband=signal%sideband )
      effectiveOpticalDepth => GetVectorQuantityByType ( fwdModelOut,      &
        & quantityType=l_effectiveOpticalDepth, noerror=.true.,            &
        & signal=signal%index, sideband=signal%sideband )

      ! -------
      ! Inputs:
      ! -------
      ptan => GetVectorQuantityByType ( fwdModelExtra,                     &
        & quantityType=l_ptan, instrumentModule = Signal%instrumentModule )
      losVel => GetVectorQuantityByType ( fwdModelExtra,                   &
        & quantityType=l_losVel, instrumentModule=Signal%instrumentModule )

      !-----------------------------------------
      ! Make sure the quantities we need are got and with correct format
      !-----------------------------------------
      if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true.,               &
         & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
         & ModuleName, InvalidQuantity//'ptan' )

      !----------------------------------
      ! Set up some temporary quantities
      !----------------------------------

      call Allocate_test ( closestInstances, radiance%template%noInstances,    &
        & 'closestInstances', ModuleName )      

      !------------------------
      ! Assemble the vmr array
      !------------------------
      !---------------------------------------------------------
      ! Work out the closest instances from temperature
      !---------------------------------------------------------
      call FindClosestInstances ( temp, radiance, closestInstances )
      instance = closestInstances(maf)
      tLat = temp%template%geodLat(1,instance)    ! get latitude for each instance

      ! ----------------------------
      ! Get some basic dimensions
      ! ----------------------------
      noMIFs  = radiance%template%noSurfs

      ! Interpolate PSD surfaces onto temperature vGrid if necessary
      if ( sizeDistribution%template%vGridIndex /= temp%template%vGridIndex ) then
        call InterpolateValues ( &
         & reshape(sizeDistribution%template%surfs(:,1),(/noPSDSurf/)), & ! Old X
         & reshape(sizeDistribution%values(:,instance),(/noPSDSurf/)),  & ! Old Y
         & reshape(temp%template%surfs(:,1),(/noSurf/)),                & ! New X
         & PSD, 'Linear')                                                 ! New Y
        call MLSMessage(MLSMSG_Warning, ModuleName, &
        &'Particle size dist vGrid is NOT same as Temperature and interpolation is enforced')
      else
          PSD = sizeDistribution%values(:,instance)
      end if

      if ( radiance%template%nosurfs /= ptan%template%nosurfs ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'number of levels in radiance does not match number of levels in ptan' )

      allocate ( a_Trans(noSgrid, noMIFs, noFreqs), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'a_Trans' )

      ! now checking spectroscopy

      ! When using Bill's Spectral data, work out which spectroscopy we're going to need
      if ( forwardModelConfig%default_spectroscopy ) then  !Bill's clear-sky spectroscopy

        my_catalog = empty_cat
        do thisSideband = sidebandStart, sidebandStop, 2
          do j = 1, nspec      ! Loop over species
            ! ForwardModelSupport prevents beta grouping for cloud models.
            k = forwardModelConfig%molecules(j)
            if ( catalog(k)%molecule == l_none ) then
              call get_string ( lit_indices(k), molName )
              Call MLSMessage ( MLSMSG_Warning, ModuleName, &
                & 'No spectroscopy catalog for '//trim(molName) )
            end if

            if ( associated ( catalog(k)%lines ) ) then
              ! Now subset the lines according to the signal we're using
              lineFlag => maxLineFlag(:size(catalog(k)%lines))
              lineFlag = .FALSE.
              do k = 1, size ( catalog(k)%lines )
                thisLine => lines(catalog(k)%lines(k))
                if ( associated(thisLine%signals) ) then
                  doThis = any ( ( thisLine%signals == signal%index ) .and. &
                      & ( ( thisLine%sidebands == 0 ) .or. &
                      &   ( thisLine%sidebands == thisSideband ) ) )
                    ! If we're only doing one sideband, maybe we can remove some more lines
                    if ( sidebandStart==sidebandStop ) doThis = doThis .and. &
                    & any( ( thisLine%sidebands == sidebandStart ) .or. &
                    & ( thisLine%sidebands == 0 ) )
                  lineFlag(k) = lineFlag(k) .or. doThis
                end if
              end do               ! End loop over lines

              My_Catalog(thisSideband,j) = catalog(k)
              nullify ( my_catalog(thisSideband,j)%lines ) ! Don't deallocate it by mistake 

              ! Check we have at least one line for this

              call allocate_test ( my_catalog(thisSideband,j)%lines, count(lineFlag),&
                & 'my_catalog(?,?)%lines', ModuleName )
              my_catalog(thisSideband,j)%lines = pack ( catalog(k)%lines, lineFlag )
            else
              ! No lines for this species
              lineFlag => maxLineFlag(1:0)
              my_catalog(thisSideband,j) = catalog(k)
              nullify ( my_catalog(thisSideband,j)%lines ) ! Don't deallocate it by mistake
              call Allocate_test ( my_catalog(thisSideband,j)%lines, 0, &
                & 'my_catalog(?,?)%lines(0)', ModuleName )
            end if

            if ( count(lineFlag) == 0 ) then
              call get_string ( lit_indices( abs(forwardModelConfig%molecules(j)) ), &
                & molName )
              call getNameOfSignal ( signal, SignalName, sideband=thisSideband )
              call MLSMessage ( MLSMSG_Warning, ModuleName, &
                & 'No lines for '//trim(molName)//' and '//trim(signalName) )
            end if

          end do ! End of Loop over species
        end do ! Loop over sidebands

      else          ! cloudy-sky spectroscopy

        do j = 1, nspec      ! Loop over species
          select case ( forwardModelConfig%molecules(j) )
          case ( l_h2o )
            ispec = 1
          case ( l_o3 )
            ispec = 2
          case ( l_n2o )
            ispec = 3
          case default
            cycle
          end select

          vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,            &
            & quantityType=l_vmr, molecule=forwardModelConfig%molecules(j) )

          if ( .not.ValidateVectorQuantity( vmr, stacked=.true., coherent=.true., &
             & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,  &
             & ModuleName, InvalidQuantity//'vmr' )

          novmrSurf = vmr%template%nosurfs

          ! Interpolate vmr onto same grid as temperature if necessary
          if ( vmr%template%vGridIndex /= temp%template%vGridIndex ) then
            call InterpolateValues ( &
            & reshape(vmr%template%surfs(:,1),(/novmrSurf/)), &    ! Old X
            & reshape(vmr%values(:,instance),(/novmrSurf/)),  &    ! Old Y
            & reshape(temp%template%surfs(:,1),(/noSurf/)),   &    ! New X
            & vmrArray(ispec,:),                              &    ! New Y
            & 'Linear', extrapolate='Clamp' )
          else
            vmrArray(ispec,:) = vmr%values(:,instance)
          end if

        end do ! End of Loop over species
      end if

      call allocate_test ( doChannel, noFreqs, 'doChannel', ModuleName )

      allocate ( zt(noMifs) )

      doChannel = .true.
      if ( associated ( signal%channels ) ) doChannel = signal%channels

      !-----------------------------------------------------
      ! Make temporary arrays for the cloud forward model
      !-----------------------------------------------------
      call Allocate_test ( a_clearSkyRadiance,                                 &
        & radiance%template%noSurfs, noFreqs,                                  &
        & 'a_clearSkyRadiance', ModuleName )
      call Allocate_test ( a_cloudInducedRadiance,                             &
        & radiance%template%noSurfs, noFreqs,                                  &
        & 'a_cloudInducedRadiance', ModuleName )
      call Allocate_test ( a_effectiveOpticalDepth,                            &
        & radiance%template%noSurfs, noFreqs,                                  &
        & 'a_effectiveOpticalDepth', ModuleName )
      call Allocate_test ( a_cloudRADSensitivity,                              &
        & radiance%template%noSurfs, noFreqs,                                  &
        & 'a_cloudRADSensitivity', ModuleName )
      call Allocate_test ( a_totalExtinction,                                  &
        & temp%template%noSurfs, noFreqs,                                      &
        & 'a_totalExtinction', ModuleName )
      call Allocate_test ( a_cloudExtinction,                                  &
        & temp%template%noSurfs, noFreqs,                                      &
        & 'a_cloudExtinction', ModuleName )

      do i = 1, ForwardModelConfig%no_cloud_species
        if ( noSurf == noIWCSurf ) then
          if ( i == 1 ) WC (i,:) = CloudIce%values(:,instance)
          if ( i == 2 ) WC (i,:) = CloudWater%values(:,instance)
        else
          if ( i == 1 ) &
            & call InterpolateValues ( &
            & reshape(CloudIce%template%surfs(:,1),(/noIWCSurf/)), &  ! Old X
            & reshape(CloudIce%values(:,instance),(/noIWCSurf/)),  &  ! Old Y
            & reshape(temp%template%surfs(:,1),(/noSurf/)),      &    ! New X
            & WC (i,:), 'Linear' )                                    ! New Y
          if ( i == 2 ) &
            & call InterpolateValues ( &
            & reshape(CloudWater%template%surfs(:,1),(/noLWCSurf/)), & ! Old X
            & reshape(CloudWater%values(:,instance),(/noLWCSurf/)),  & ! Old Y
            & reshape(temp%template%surfs(:,1),(/noSurf/)),          & ! New X
            & WC (i,:), 'Linear' )                                     ! New Y
        end if      
      end do

      if ( prt_log ) print*, 'jacobian is true'

      ! get tangent height from tangent pressure
      call InterpolateValues ( &
        & reshape(gph%template%surfs(:,1),(/noSurf/)), &    ! Old X
        & reshape(gph%values(:,instance),(/noSurf/)),  &    ! Old Y
        & reshape(ptan%values(:,maf),(/noMifs/)),      &    ! New X
        & zt, &                                             ! New Y
        & 'Linear' )

      if ( present(jacobian) ) then
        ! transmission functions for retrieval are calculated using artificial
        ! cloud profiles depending on where the measurement is taken in latitude
        if ( tLat >= -30._r8 .and. tLat <= 30._r8 ) then
          CloudType='Convective'
        else
          CloudType='Frontal'
        end if

        ! find cloud top index from observed Tcir, threshold to be determined
        !     for high Zt, use Tcir(maf)
        !     for low Zt, use Tcir(maf-2)
        ! --------------------------------------------------------------------
        ! find the cloud top using the cloudRadiance in fwdModelExtra 
        ! (should be only one cloud radiance quantity there)
        ! --------------------------------------------------------------------
!       constrainCldRad => GetVectorQuantityByType ( fwdModelExtra,     &           
!           & quantityType=l_cloudInducedRadiance)                                  

!       if ( .not. associated(constrainCldRad) ) then
!          call MLSMessage( MLSMSG_Error, ModuleName,                             & 
!                       'Need cloud radiances to estimate cloud top in retrieval' ) 
!       else                                                                        
        iCloudHeight = 0
!       if ( forwardModelConfig%cloud_der == l_iwc_high_height ) then
!         do mif = 1, noMifs
!           if ( constrainCldRad%values(whichchannel+(mif-1)*noFreqs,maf) /= 0.0_r8) &
!                & iCloudHeight = mif
!         end do
!       else
!         do mif = 1, noMifs
!           if ( constrainCldRad%values(whichchannel+(mif-1)*noFreqs,maf) /= 0.0_r8) & 
!                & iCloudHeight = mif
!         end do
!       end if
        ! if no cloud is found, cloud top is 20 km
        CloudHeight = 20.e3_r8     ! meters
        if ( iCloudHeight /= 0 ) CloudHeight = min(zt(iCloudHeight),CloudHeight)

        ! set up artificial cloud profile for retrieval use only
        call CLOUD_MODEL (CloudType, CloudHeight, gph%values(:,instance), &
          & noSurf, WC, ForwardModelConfig%no_cloud_species)
!       end if
      
      end if

      call Allocate_test ( frequencies, noFreqs, 'frequencies', ModuleName )

      ! It is equavelent to shift either spectral line frequencies or filter frequencies
      ! but the sign is opposite.
      ! Vel_Cor = 1.0_rp - losVel%values(1,maf)/c  ! for shift spectral lines
      Vel_Cor = 1.0_rp + losVel%values(1,maf)/c    ! for shift spectral channnels

      !--------------------------------------------
      ! Loop over sidebands 
      !--------------------------------------------

      call allocate_test ( thisFraction, noFreqs, 'thisFraction', ModuleName )

      if ( signal%sideband == 0 .or. forwardModelConfig%forceSidebandFraction ) then
        sidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_limbSidebandFraction, signal=signal%index, noError=.true. )
        lowerSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_limbSidebandFraction, signal=signal%index, &
          & sideband=-1, noError=.true. )
        upperSidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType = l_limbSidebandFraction, signal=signal%index, &
          & sideband=1, noError=.true. )
        if ( .not. associated (sidebandFraction) .and. .not. &
          & ( associated ( lowerSidebandFraction) .and. associated ( upperSidebandFraction ) ) ) &
          & call MLSMessage(MLSMSG_Error,ModuleName, &
          & "No sideband ratio supplied")
      end if

      if ( toggle(emit) .and. levels(emit) > 0 ) then
        call Trace_Begin ( 'FullCloudForwardModel.SidebandLoop' )
      end if

      do thisSideband = sidebandStart, sidebandStop, 2

        if ( toggle(emit) .and. levels(emit) > 1 ) &
          & call Trace_Begin ( 'FullCloudForwardModel.Sideband ', index=thisSideband )

        !-----------------------------
        ! Setup a sideband ratio array
        !------------------------------

        elevOffset => GetVectorQuantityByType ( fwdModelExtra,               &
          & quantityType=l_elevOffset, signal=Signal%index, sideband=thisSideband )

        if ( signal%sideband == 0 .or. forwardModelConfig%forceSidebandFraction ) then
            if ( thisSideband == -1 ) then
              thisFraction = lowerSidebandFraction%values(:,1)
            else
              thisFraction = upperSidebandFraction%values(:,1)
            end if
        else                  ! Otherwise, want just unfolded signal
          thisFraction = 1.0
        end if

        direction = signal%direction

        frequencies = signal%centerFrequency + direction*signal%frequencies 

        if ( signal%sideband == 0 ) &
        frequencies = signal%lo + thisSideband * frequencies    ! double sideband cases

        if ( signal%sideband /= 0 ) &
        frequencies = signal%lo + signal%sideband * frequencies    ! signal sideband cases

        !--------------------------------------------
        ! VELOCITY shift correction to frequencies
        !--------------------------------------------

        frequencies =  Vel_Cor * frequencies

        if ( prt_log ) then
           print*, ' '
           print*,'No. of Frequencies:', noFreqs 
           print*,'Frequency=', frequencies/1e3_r8
        end if

        do j = 1, size(antennaPatterns)
          superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
               & ForwardModelConfig%signals( (/sigInd/) ), sideband=thisSideband, channel =1 )
        end do

        whichPattern = 1
        maxSuperset = maxval ( superset )
        where ( superset < 0 )
              superset = maxSuperset + 1
        end where
        whichPatternAsArray = minloc ( superset )
        whichPattern = whichPatternAsArray(1)
!        if ( toggle(emit) .and. levels(emit) > 2 ) then
!           call output ( 'Using antenna pattern: ' )    
!           call output ( whichPattern, advance='yes' )  
!        end if                                          

        !---------------------------------------------
        ! Now call the Full CloudForwardModel routine
        !---------------------------------------------
        !print*,instance,maxval(wc)

        call CloudForwardModel ( doChannel,                               &
          & noFreqs,                                                      &
          & noSurf,                                                       & 
          & noMifs,                                                       &
          & nspec,                                                        &
          & ForwardModelConfig%no_cloud_species,                          &
          & ForwardModelConfig%no_model_surfs,                            &
          & frequencies/1e3_r8,                                           &
          & 10.0**(-temp%template%surfs),                                 &
          & gph%values(:, instance),                                      &
          & temp%values(:,instance),                                      &
          & vmrArray,                                                     &
          & WC,                                                           &
          & int(PSD),                                                     &
          & zt,                                                           &
          & earthradius%values(1,1),                                      &
          & int(surfaceType%values(1, instance)),                         &
          & forwardModelConfig%cloud_der,                                 &
          & forwardModelConfig%i_saturation,                              &
          & forwardModelConfig%do_conv,                                   &
          & forwardModelConfig%default_spectroscopy,                      &
          & scGeocAlt%values(1,1),                                        &
          & elevOffset%values(1,1),                                       &
          & antennaPatterns(whichPattern),                                &
          & a_clearSkyRadiance,                                           &
          & a_cloudInducedRadiance,                                       &
          & a_trans,                                                      &
          & a_totalExtinction,                                            &
          & a_cloudExtinction,                                            &
          & a_massMeanDiameter,                                           &
          & a_effectiveOpticalDepth,                                      &
          & a_cloudRADSensitivity,                                        &
          & forwardModelConfig%NUM_SCATTERING_ANGLES,                     &  
          & forwardModelConfig%NUM_AZIMUTH_ANGLES,                        &
          & forwardModelConfig%NUM_AB_TERMS,                              &
          & forwardModelConfig%NUM_SIZE_BINS,                             &
          & Slevl*1000._r8, noSgrid,                                      &
          & My_Catalog(thisSideband,:),                                   &
          & losVel%values(1,1) )   
!          & losVel%values(1,maf) )
                                
        if ( prt_log ) print*, 'Successfully done with Full Cloud Foward Model ! '

        ! Output minor frame quantities if they are asked: check channel and type

        if ( sidebandStart == sidebandStop ) then
          radiance%values ( :, maf) =                                              &
            & reshape ( transpose(a_clearSkyRadiance),                             &
            & (/radiance%template%instanceLen/) )
          modelCloudRadiance%values ( :, maf ) =                                   &
            & reshape ( transpose(a_cloudInducedRadiance),                         &
            & (/modelCloudRadiance%template%instanceLen/) )
          if ( associated(effectiveOpticalDepth) )                                 &
            & effectiveOpticalDepth%values ( :, maf ) =                            &
            & reshape ( transpose(a_effectiveOpticalDepth),                        &
            & (/effectiveOpticalDepth%template%instanceLen/) )
          if ( associated(cloudRADSensitivity) )                                   &
            & cloudRADSensitivity%values ( :, maf ) =                              &
            & reshape ( transpose(a_cloudRADSensitivity),                          &
            & (/cloudRADSensitivity%template%instanceLen/) )
!          print*,maf,thissideband,a_clearSkyRadiance(:,25)
        else

          do i =1 , noFreqs
            if ( doChannel(i) ) then
              do mif=1, noMIFs

                radiance%values (i+(mif-1)*noFreqs, maf) =                         &
                  &               radiance%values (i+(mif-1)*noFreqs, maf)         &
                  &             + thisFraction(i)*a_clearSkyRadiance(mif,i) 

                modelCloudRadiance%values (i+(mif-1)*noFreqs, maf ) =              &
                  &        modelCloudRadiance%values (i+(mif-1)*noFreqs, maf )     &
                  &      + thisFraction(i)*a_cloudInducedRadiance(mif,i)

                if ( associated(effectiveOpticalDepth) )                              &
                  & effectiveOpticalDepth%values (i+(mif-1)*noFreqs, maf ) =       &
                  &        effectiveOpticalDepth%values (i+(mif-1)*noFreqs, maf )  &
                  &      + thisFraction(i)*a_effectiveOpticalDepth(mif,i)

                if ( associated(cloudRADSensitivity) )                                &
                  & cloudRADSensitivity%values (i+(mif-1)*noFreqs, maf ) =         &
                  &        cloudRADSensitivity%values (i+(mif-1)*noFreqs, maf )    &  
                  &      + thisFraction(i)*a_cloudRADSensitivity(mif,i)

              end do
!   print*,maf,thissideband,i,thisFraction(i),a_cloudRADSensitivity(1,i),a_clearSkyRadiance(1,i)
            end if
          end do
! print*,maf,thissideband,a_clearSkyRadiance(:,25)
        end if

        if ( toggle(emit) .and. levels(emit) > 1 ) &
          & call Trace_End ( 'FullCloudForwardModel.Sideband ',index=thisSideband )

      !--------------------------------------------
      ! End of sideband loop 
      !--------------------------------------------
      end do

      if ( toggle(emit) .and. levels(emit) > 0 ) &      
        & call Trace_End ( 'FullCloudForwardModel.SidebandLoop' )

      ! -----------------------------------------------------------------------------
      ! output L2GP quantities
      ! -----------------------------------------------------------------------------

      ! To save disk space, we output one channel per band (first true doChannel)
      ! in the last signal (band) will over write all the previous signals (bands)

      FOUNDINFIRST = .true.
      do i=1, noFreqs
        if ( doChannel(i) .and. FOUNDINFIRST ) then
          FOUNDINFIRST = .false.
          if ( associated(cloudExtinction) ) &
            & cloudExtinction%values (:, instance )    = a_cloudExtinction(:,i)
          if ( associated(totalExtinction) ) &
            & totalExtinction%values (:, instance )    = a_totalExtinction (:,i)
        end if
      end do

      if ( associated(massMeanDiameterIce) ) &
        & massMeanDiameterIce%values (:,instance)   = a_massMeanDiameter(1,:)
      if ( associated(massMeanDiameterWater) ) &
        & massMeanDiameterWater%values(:, instance) =  a_massMeanDiameter(2,:)

      !-----------------------
      ! Start output Jacobian
      !-----------------------

        noInstances = state_ext%template%noInstances

      !--------------------------------------------
      ! Jacobian for high tangent height retrieval (only one frequency)
      !--------------------------------------------

      if ( forwardModelConfig%cloud_der == l_iwc_high_height &
        & .and. present(jacobian) ) then
      ! do not handle multiple signals and channels  
        if ( size(forwardModelConfig%signals) > 1 .and. size(doChannel) > 1 ) then
           print*,'only one frequency and one signal is allowed'
           stop
        end if

        colJBlock = FindBlock ( Jacobian%col, state_ext%index, maf)
        rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
        fmStat%rows(rowJBlock) = .true.

        jBlock => jacobian%block(rowJblock,colJblock)

          call CreateBlock ( jBlock, noMIFs, noExtSurf*noInstances, M_Full )
          jBlock%values = 0._rm

        !-------------------------------------------------------------------
        ! we use 100 times better resolution to compute weighting functions
        !-------------------------------------------------------------------
        nfine = 100
        nNear = ForwardModelConfig%phiWindow         ! default = 5
        ! only nearest instances are mattered
          minInst = instance - (nNear-1)/2
          maxInst = instance + (nNear-1)/2
           if ( minInst < 1 ) minInst = 1
           if ( maxInst > noInstances ) maxInst = noInstances

        allocate( phi_fine(nfine*nNear), stat=status )
        allocate( z_fine(nfine*nNear), stat=status )
        allocate( zp_fine(nfine*nNear), stat=status )
        allocate( s_fine(nfine*nNear), stat=status )
        allocate( w_fine(nfine*nNear), stat=status )
        allocate( ds_fine(nfine*nNear), stat=status )

        do i=1,nNear*nfine
         phi_fine(i) = state_ext%template%phi(1,minInst) + &
           & (i-1._r8)/nfine/nNear*(state_ext%template%phi(1,maxInst) - &
           & state_ext%template%phi(1,minInst))
        end do

        !----------------------------------
        ! find vertical and horizontal intervals of stateQ at maf
        !----------------------------------
        ! vertical
        dz = abs(state_ext%template%surfs(2,1)-state_ext%template%surfs(1,1))
        ! horiozontal: in case of the last instance, use the previous one
        if ( instance < noInstances ) then
           dphi = abs(state_ext%template%phi(1,instance+1)-state_ext%template%phi(1,instance))
        else
           dphi = abs(state_ext%template%phi(1,instance-1)-state_ext%template%phi(1,instance))
        end if

        do mif = 1, noMIFs
          ! Jacobians at only tangent heights less than the top level of retrieval are calculated
          if ( state_ext%template%surfs(noExtSurf,1) > ptan%values(mif,maf) .and. &
             & state_ext%template%surfs(1,1) < ptan%values(mif,maf) ) then
            !------------------------------
            ! find z for given phi_fine
            !------------------------------
            z_fine = (earthradius%values(1,maf)+ zt(mif)) / &
             & cos((phi_fine - radiance%template%phi(mif,maf))*Deg2Rad) - &
             & earthradius%values(1,maf)
             ! convert back to log tangent pressure
             call InterpolateValues ( &
                & reshape(gph%values(:,instance),(/noSurf/)), &     ! Old X
                & reshape(gph%template%surfs(:,1),(/noSurf/)), &    ! Old Y
                & z_fine, &                                         ! New X
                & zp_fine, &                                        ! New Y
                & 'Linear' )

            !--------------------------------
            ! find ds and weight for each (z,phi) pair
            !--------------------------------
            s_fine = (earthradius%values(1,maf)+ zt(mif)) * &
             & sin((phi_fine - radiance%template%phi(mif,maf))*Deg2Rad)
            ds_fine = 0._r8    ! initialize it
            do i=1,nfine*nNear-1
             ds_fine(i) = s_fine(i+1) - s_fine(i)
            end do

             call InterpolateValues ( &
                & sLevl, &                                               ! Old X
                & reshape(a_trans(:,mif,whichChannel),(/noSgrid/)), &    ! Old Y
                & s_fine, &                                              ! New X
                & w_fine, &                                              ! New Y
                & 'Linear' )

            ! ds needs to be weighted by transmission function

            ds_fine = ds_fine*w_fine    ! ds is in meters

            !----------------------------------------------------------
            ! determine weights by the length inside each state domain
            !----------------------------------------------------------

             do i = minInst,maxInst             ! loop over closer profiles
             do j = 1,noExtSurf                 ! loop over cloudQty surface
             do k = 1,nfine*nNear               ! sum up all the lengths
               if ( abs(zp_fine(k) - state_ext%template%surfs(j,1)) < dz/2._r8 &
               & .AND. abs(phi_fine(k) - state_ext%template%phi(1,i)) < dphi/2._r8 ) &
               & jBlock%values(mif,j+(i-1)*noExtSurf) = &
               & jBlock%values(mif,j+(i-1)*noExtSurf) + ds_fine(k)/1.e3_r8
             end do
             end do
             end do
          end if
        end do         ! mif

        Deallocate( phi_fine, stat=status )
        Deallocate( z_fine, stat=status )
        Deallocate( zp_fine, stat=status )
        Deallocate( ds_fine, stat=status )
        Deallocate( w_fine, stat=status )
        Deallocate( s_fine, stat=status )

      end if     ! high tangent height case

      !--------------------------------------------
      ! Jacobian for low tangent height retrieval
      !--------------------------------------------
      if ( forwardModelConfig%cloud_der == l_iwc_low_height &
        & .and. present(jacobian) ) then

        colJBlock = FindBlock ( Jacobian%col, state_los%index, maf )
        rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
        fmStat%rows(rowJBlock) = .true.

        jBlock => jacobian%block(rowJblock,colJblock)

      ! to save space, the jacobian is packed in a full rectangle matrix
        call CreateBlock ( jBlock, noFreqs, noSgrid*noMIFs, M_Full )
        jBlock%values = 0.0_rm

        do j = 1, noFreqs
          if ( doChannel(j) ) then
             do mif = 1, noMIFs
               do i=1,noSgrid
               ! now we normalize cloud extinction weighting functions at 200GHz 
               ! and output the transmission functions via Jacobian
                 jBlock%values(j,i+(mif-1)*noSgrid)= a_trans(i,mif,j)* &
                    & (frequencies(j)/200000._r8)**4
               end do
             end do
          end if  ! doChannel
        end do    ! channel

      end if      ! low tangent height case

      !------------------------------
      ! End of output jacobian
      !------------------------------

      !------------------------------
      ! Remove temporary quantities
      !------------------------------

      call Deallocate_test ( a_cloudExtinction,       'a_cloudExtinction',       ModuleName )
      call Deallocate_test ( a_totalExtinction,       'a_totalExtinction',       ModuleName )
      call Deallocate_test ( a_cloudRADSensitivity,   'a_cloudRADSensitivity',   ModuleName )
      call Deallocate_test ( a_effectiveOpticalDepth, 'a_effectiveOpticalDepth', ModuleName )
      call Deallocate_test ( a_cloudInducedRadiance,  'a_cloudInducedRadiance',  ModuleName )
      call Deallocate_test ( a_clearSkyRadiance,      'a_clearSkyRadiance',      ModuleName )
      call Deallocate_test ( closestInstances,        'closestInstances',        ModuleName )
      call Deallocate_test ( doChannel,               'doChannel',               ModuleName )
      call Deallocate_test ( thisFraction,            'thisFraction',            ModuleName )
      call Deallocate_test ( frequencies,             'frequencies',             ModuleName )
      deallocate (a_trans, Zt, stat=ier )
      if ( ier /= 0 ) Call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'a_trans or Zt' )

      if ( prt_log ) then
        print*, ' '
        print*, 'Time Instance: ', instance
      end if

    end do  ! End of signals

    call deallocate_test ( a_massMeanDiameter, 'a_massMeanDiameter', ModuleName )
    call deallocate_test ( superset,           'superset',           ModuleName )
    call deallocate_test ( vmrArray,           'vmrArray',           ModuleName )

    deallocate (PSD, Slevl, WC, stat=ier )
    if ( ier /= 0 ) Call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'PSD, Slevl, WC, or Zt' )

    do i = 1, size(my_catalog,2)
      do thisSideband = sidebandStart, sidebandStop, 2
        call deallocate_test ( my_catalog(thisSideband,i)%lines, &
          & 'my_catalog(?,?)%lines', ModuleName )
        ! Note that we don't deallocate the signals/sidebands stuff for each line
        ! as these are shallow copies of the main spectroscopy catalog stuff
      end do
    end do
    deallocate ( my_catalog, stat=ier )
    if ( ier /= 0 ) Call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'my_catalog' )

    if ( toggle(emit) ) call trace_end ( 'FullCloudForwardModel' )

    if ( prt_log ) print*, 'Successfully done with full cloud forward wrapper !'

  end subroutine FullCloudForwardModelWrapper

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module FullCloudForwardModel


! $Log$
! Revision 1.134  2004/11/03 01:26:32  vsnyder
! Missed one size(molecules)-1 that should be size(molecules)
!
! Revision 1.133  2004/11/01 20:26:57  vsnyder
! Reorganization of representation for molecules and beta groups
!
! Revision 1.132  2004/10/21 01:33:35  vsnyder
! Move fetching of elevOffset to after trace_begin, some cannonball polishing
!
! Revision 1.131  2004/10/20 23:44:13  pwagner
! elevOffset Vector now found by signal, sideband
!
! Revision 1.130  2004/10/13 01:08:05  vsnyder
! Moved some checking to ForwardModelSupport
!
! Revision 1.129  2004/08/05 21:03:07  vsnyder
! Add sentinel at end of %molecules
!
! Revision 1.128  2004/06/10 01:00:14  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 1.127  2004/02/07 00:45:49  livesey
! Minor typo
!
! Revision 1.126  2004/02/05 23:30:20  livesey
! Finally turned on code to do correct handing of single sideband
! radiometers and their SB fractions.
!
! Revision 1.125  2004/01/08 00:52:43  jonathan
! delet output Using antenna pattern: etc
!
! Revision 1.124  2004/01/08 00:24:32  jonathan
! add tracing signals
!
! Revision 1.123  2003/11/21 19:17:40  dwu
! fix a bug from the last modification
!
! Revision 1.122  2003/11/21 17:55:57  dwu
! add vGrid check for size distribution
!
! Revision 1.121  2003/11/20 20:30:58  dwu
! fit the problem when IWC and Temp vGrid mismatch
!
! Revision 1.120  2003/07/09 20:17:22  livesey
! Anticipative bug fix in sideband fraction
!
! Revision 1.119  2003/05/29 16:38:03  livesey
! Renamed sideband fraction
!
! Revision 1.118  2003/05/16 23:53:56  livesey
! Removed reference to spectags
!
! Revision 1.117  2003/05/15 22:48:28  dwu
! changes
!
! Revision 1.116  2003/05/07 22:49:03  jonathan
! some clean-up and cosmetic changes
!
! Revision 1.115  2003/04/10 20:40:59  dwu
! make i_saturation and cloud_der as a verbal argument
!
! Revision 1.114  2003/04/08 20:03:55  dwu
! fix a bug in handling no of cloud species. now this number is meanful
!
! Revision 1.113  2003/04/05 17:30:43  dwu
! clean up
!
! Revision 1.112  2003/04/03 22:38:05  dwu
! change the way vmr and molecule are handled
!
! Revision 1.111  2003/04/03 01:18:41  dwu
! allow molecules in random order
!
! Revision 1.110  2003/04/03 01:15:51  dwu
! allow molecules in random order
!
! Revision 1.109  2003/04/02 20:00:12  dwu
! some clearup and replace ifov with do_conv
!
! Revision 1.108  2003/01/23 00:19:09  pwagner
! Some cosmetic only (or so I hope) changes
!
! Revision 1.107  2003/01/17 07:19:33  dwu
! properly initialize vmrArray
!
! Revision 1.106  2003/01/17 01:08:11  jonathan
! add vmrArray=0. after allocate
!
! Revision 1.105  2003/01/17 00:52:05  jonathan
! if specices missing, set VMR=0 accordingly
!
! Revision 1.104  2003/01/16 18:39:41  pwagner
! Removed some unused variables
!
! Revision 1.103  2003/01/13 17:59:52  jonathan
!  change cloud_width to i_saturation
!
! Revision 1.102  2003/01/09 21:08:09  dwu
! drop orbI
!
! Revision 1.101  2002/12/18 16:09:24  jonathan
! add phi_tan
!
! Revision 1.100  2002/11/30 21:31:46  dwu
! fix signal loop and move constrainCldRad inside jacobian loop
!
! Revision 1.99  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.98  2002/10/04 00:02:15  vsnyder
! Change handling of GOT variable, some cosmetic changes
!
! Revision 1.97  2002/10/03 23:25:52  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 1.96  2002/09/11 17:43:39  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 1.95  2002/08/22 00:14:17  jonathan
! upgrade to include more molecules
!
! Revision 1.94  2002/08/19 22:22:03  jonathan
! debug stuff
!
! Revision 1.93  2002/08/08 22:46:30  jonathan
! newly improved version
!
! Revision 1.92  2002/05/08 17:03:39  jonathan
! use earthradius(1,1) for now
!
! Revision 1.91  2002/01/14 19:30:16  jonathan
! minor changes
!
! Revision 1.90  2001/11/16 00:49:05  jonathan
! clean log
!
! Revision 1.89  2001/11/16 00:47:13  jonathan
! change ptan from radiance%template%instrumentModule to Signal%instrumentModule
!
! Revision 1.88  2001/11/16 00:41:00  jonathan
! add losVel
!
! Revision 1.87  2001/11/15 23:52:12  jonathan
! rename DF_spectroscopy to default_spectroscopy
!
! Revision 1.86  2001/11/15 23:50:21  jonathan
! add DF_spectroscopy
!
! Revision 1.85  2001/11/09 18:12:13  jonathan
! add deallocate my_catalog
!
! Revision 1.84  2001/11/09 18:05:24  jonathan
! pass spectra catalog into CloudySkyRadianceModel
!
! Revision 1.83  2001/11/08 21:36:13  jonathan
! add SpectroscopyCatalog
!
! Revision 1.82  2001/11/07 23:47:52  dwu
! some minor changes
!
! Revision 1.81  2001/11/07 05:22:06  dwu
! fixed a bug in computing dphi
!
! Revision 1.80  2001/11/06 21:54:47  dwu
! use phiWindow to save time
!
! Revision 1.79  2001/11/06 20:06:31  dwu
! speed up Jacobian calculation for high tangent heights
!
! Revision 1.78  2001/11/06 19:52:42  dwu
! speed up Jacobian calculation for high tangent heights
!
! Revision 1.77  2001/11/06 18:28:06  dwu
! some cleanups
!
! Revision 1.76  2001/11/06 00:54:11  dwu
! add two cloud radiances: modelled and observed
!
! Revision 1.75  2001/11/06 00:29:38  dwu
! set up cloud height estimation using DTcir
!
! Revision 1.74  2001/11/05 22:39:57  dwu
! high Zt Jacobian
!
! Revision 1.73  2001/11/05 20:25:14  dwu
! *** empty log message ***
!
! Revision 1.72  2001/11/02 01:14:17  dwu
! correction in high cloud Jacobian
!
! Revision 1.71  2001/11/02 01:00:13  jonathan
! add IWC1
!
! Revision 1.70  2001/11/02 00:47:34  dwu
! correction in high cloud Jacobian
!
! Revision 1.69  2001/11/02 00:45:33  dwu
! correction in high cloud Jacobian
!
! Revision 1.68  2001/11/02 00:42:05  dwu
! correction in high cloud Jacobian
!
! Revision 1.67  2001/11/02 00:25:58  dwu
! correction in high cloud Jacobian
!
! Revision 1.66  2001/10/30 05:25:48  dwu
! assign whichchannle
!
! Revision 1.65  2001/10/19 19:29:50  dwu
! initialize output quantities
!
! Revision 1.64  2001/10/19 16:26:09  dwu
! some minors
!
! Revision 1.63  2001/10/18 22:17:02  dwu
! pretection for sensitivity=0
!
! Revision 1.62  2001/10/18 06:06:24  dwu
! minor
!
! Revision 1.61  2001/10/12 16:58:31  dwu
! distinguish number surfaces between model temperature grid and retrieval ext grid
!
! Revision 1.60  2001/10/11 22:44:01  dwu
! modify high zt Jacobian and use cloud_der as the switch between high and low Zt
!
! Revision 1.59  2001/10/11 17:01:11  jonathan
! for (cloud/total)extinction/output one channel per band
!
! Revision 1.58  2001/10/10 18:55:17  dwu
! why elevOffset is not maf-dependent?
!
! Revision 1.57  2001/10/10 18:33:53  dwu
! normalize cloud extinction weighting function to 200GHz
!
! Revision 1.56  2001/10/09 22:11:54  jonathan
! *** empty log message ***
!
! Revision 1.55  2001/10/09 17:50:17  jonathan
! *** empty log message ***
!
! Revision 1.54  2001/10/09 17:47:33  jonathan
! some changes
!
! Revision 1.53  2001/10/08 23:43:00  dwu
! fix jBlock%kind initialization
!
! Revision 1.52  2001/10/08 21:46:39  jonathan
! add CloudySkyModule
!
! Revision 1.50  2001/10/08 21:42:29  dwu
! *** empty log message ***
!
! Revision 1.49  2001/10/08 21:35:33  dwu
! *** empty log message ***
!
! Revision 1.48  2001/10/08 20:57:05  dwu
! change coljBlock finder
!
! Revision 1.47  2001/10/08 20:46:26  dwu
! add default cloudheight as 18km
!
! Revision 1.46  2001/10/08 20:34:18  dwu
! *** empty log message ***
!
! Revision 1.45  2001/10/08 20:23:54  jonathan
! some changes
!
! Revision 1.44  2001/10/08 19:25:48  jonathan
! delet vmrINS
!
! Revision 1.43  2001/10/07 23:42:17  jonathan
! add CloudProfile module
!
! Revision 1.42  2001/10/05 22:25:40  dwu
! allow multiple signals
!
! Revision 1.41  2001/10/05 20:46:39  dwu
! clean up input statements
!
! Revision 1.40  2001/10/05 20:26:12  dwu
! make sure the model output fields are associated
!
! Revision 1.39  2001/10/04 23:34:19  dwu
! *** empty log message ***
!
! Revision 1.38  2001/10/04 16:27:12  jonathan
! added framework for double sideband calculation, unfinished
!
! Revision 1.37  2001/10/04 00:29:36  dwu
! fix coljBlock
!
! Revision 1.36  2001/10/02 17:08:02  jonathan
! some adjustment due to construction
!
! Revision 1.35  2001/10/02 16:27:36  livesey
! Removed reference to fmStat%finished
!
! Revision 1.34  2001/10/01 23:40:26  jonathan
! construct codes for double sideband
!
! Revision 1.33  2001/09/28 21:45:50  dwu
! modify low cloud Jacobian output format
!
! Revision 1.32  2001/09/28 15:54:39  jonathan
! minor
!
! Revision 1.31  2001/09/26 19:17:02  dwu
! normalize weights for high tangent retrieval
!
! Revision 1.30  2001/09/24 23:16:40  dwu
! add derivatives for high tangent height retrievals
!
! Revision 1.29  2001/09/24 23:12:53  dwu
! add derivatives for high tangent height retrievals
!
! Revision 1.28  2001/09/21 15:51:37  jonathan
! modified F95 version
!
! Revision 1.27  2001/09/19 16:46:22  dwu
! some minor
!
! Revision 1.26  2001/09/19 00:25:59  dwu
! add M_banded to Jacobian
!
! Revision 1.25  2001/09/04 15:59:44  jonathan
! add cloud_fov, jonathan
!
! Revision 1.24  2001/08/17 21:48:41  jonathan
! Added FOV average, Jonathan
!
! Revision 1.23  2001/08/07 17:17:50  jonathan
! add radiance%template%instrumentModule to ptan
!
! Revision 1.22  2001/08/02 01:03:16  dwu
! add doChannel to frequency loop
!
! Revision 1.21  2001/08/01 20:51:30  dwu
! add delTau100
!
! Revision 1.20  2001/08/01 17:24:29  jonathan
! updated version
!
! Revision 1.19  2001/08/01 00:20:22  dwu
! add Jacobian -Jonathan/Wu
!
! Revision 1.17  2001/07/27 22:12:13  jonathan
! fixed bug in output extinction
!
! Revision 1.16  2001/07/27 20:26:24  jonathan
! jonathan
!
! Revision 1.15  2001/07/27 15:17:58  jonathan
! First Successful f90 runs
!



