
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullCloudForwardModel

! -------------------------------------------------------------------------
! THIS MODULE CONTAINS THE FULL CLOUD FORWARD MODEL  
! Jonathan Jiang,  Paul Wagner, July 16, 2001 
! Jonathan Jiang,  Dong Wu, add Jacobian, August 3, 2001
! Jonathan Jiang,  add Field Of View Convolution, August 16, 2001  
! Jonathan Jiang,  add sideband ratio, October 4, 2001 
! Jonathan Jiang,  add CloudProfile module, October 5, 2001 
! -------------------------------------------------------------------------

  use Allocate_deallocate, only: Allocate_test, Deallocate_test
  use AntennaPatterns_m, only: ANTENNAPATTERNS
  use CloudSkyModule, only: CLOUD_MODEL
  use CloudySkyRadianceModel, only: CloudForwardModel
  use Hdf, only: DFACC_READ, DFACC_CREATE
  use HDFEOS, only: SWOPEN,     SWCLOSE
  use L2GPData, only: L2GPData_T, ReadL2GPData, WriteL2GPData
  use MLSCommon,only: NameLen,    FileNameLen, r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSSignals_m, only: SIGNAL_T, ARESIGNALSSUPERSET
  use MatrixModule_0, only: M_Absent, M_BANDED, MATRIXELEMENT_T, M_BANDED,   &
                          & M_COLUMN_SPARSE, CREATEBLOCK, M_FULL
  use MatrixModule_1, only: MATRIX_T, FINDBLOCK
  use ManipulateVectorQuantities, only: FindClosestInstances
  use MLSNumerics, only: InterpolateValues
  use Molecules, only: L_H2O, L_O3
  use Output_m, only: OUTPUT
  use PointingGrid_m, only: POINTINGGRIDS
  use Toggles, only: Emit, Levels, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Units
  use VectorsModule, only: GETVECTORQUANTITYBYTYPE,                          &
                         & VECTOR_T, VECTORVALUE_T,                          &
                         & VALIDATEVECTORQUANTITY
  
! -----------------------------------------------------------------------
! THE FOLLOWING IS MODIFICATIONS FOR THE CLOUD FORWARD MODEL PARAMETERS
! -----------------------------------------------------------------------

  use ForwardModelConfig, only: FORWARDMODELCONFIG_T   
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T,            &
                                    & FORWARDMODELSTATUS_T

! ----------------------------------------------------------
! DEFINE INTRINSIC CONSTANTS NEEDED BY Init_Tables_Module
! ----------------------------------------------------------

  use Intrinsic, only: L_TEMPERATURE,L_PTAN,L_VMR,L_GPH,L_RADIANCE,L_NONE,   &
                     & L_CLOUDINDUCEDRADIANCE,                               &
                     & L_EFFECTIVEOPTICALDEPTH,                              &
                     & L_CLOUDRADSENSITIVITY,                                &
                     & L_TOTALEXTINCTION,                                    &
                     & L_CLOUDEXTINCTION,                                    &
                     & L_MASSMEANDIAMETERICE,                                & 
                     & L_MASSMEANDIAMETERWATER,                              &
                     & L_SURFACETYPE,                                        &
                     & L_SIZEDISTRIBUTION,                                   &
                     & L_TNGTGEOCALT,                                        &
                     & L_EARTHRADIUS,                                        &
                     & L_CLOUDICE,                                           &
                     & L_CLOUDWATER,                                         &
                     & L_LOSTRANSFUNC,                                       &
   		     & L_SCGEOCALT,                                          &
		     & L_ELEVOFFSET,                                         &
                     & L_SIDEBANDRATIO,                                      &
                     & L_CHANNEL,                                            &
                     & L_NONE

  implicit none
  private

  public :: FullCloudForwardModelWrapper

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------

 ! Local parameters ---------------------------------------------------------

  character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

         ! ---------------------------------------------------------------------
contains ! THIS SUBPROGRAM CONTAINS THE WRAPPER ROUTINE FOR CALLING THE FULL
         ! CLOUD FORWARD MODEL
         ! ---------------------------------------------------------------------

  subroutine FullCloudForwardModelWrapper ( ForwardModelConfig, FwdModelIn,  &
                                            FwdModelExtra, FwdModelOut, Ifm, &
                                            fmStat, Jacobian                 )  
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT                ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM      ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT         ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Local variables
    type (VectorValue_T), pointer :: CLOUDICE                   ! Profiles
    type (VectorValue_T), pointer :: CLOUDWATER                 ! Profiles
    type (VectorValue_T), pointer :: CLOUDEXTINCTION            ! Profiles
    type (VectorValue_T), pointer :: CLOUDINDUCEDRADIANCE       ! Like radiance
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
    type (VectorValue_T), pointer :: SCGEOCALT     ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: ELEVOFFSET    ! Elevation offset quantity
    type (Signal_T) :: signal                      ! A signal
    type(MatrixElement_T), pointer :: JBLOCK       ! A block from the jacobian
    !type(VectorValue_T), pointer :: STATEQ        ! A state vector quantity
    type(VectorValue_T), pointer :: STATE_ext      ! A state vector quantity
    type(VectorValue_T), pointer :: STATE_los      ! A state vector quantity

    type(VectorValue_T), pointer :: SIDEBANDRATIO  ! From the state vector

    ! for jacobian
    integer :: COLJBLOCK                ! Column index in jacobian
    integer :: ROWJBLOCK                ! Row index in jacobian
    integer :: XINSTANCE                ! Instance in x corresponding to xStarInstance
    integer :: noInstances              ! no of instance
    integer :: noChans                  ! Dimension
    integer :: noMIFs                   ! Number of minor frames
    integer :: noSgrid                  ! no of elements in S grid
    integer :: noSurf                   ! Number of pressure levels
    integer :: NQ1
    integer :: NQ2                      ! no of quantities in extraModelIn
    integer :: NOFREQS                  ! Number of frequencies
    integer :: noNonZero                ! Number of nonzero values in Jacobian

    integer :: chan
    integer :: i                        ! Loop counter
    integer :: j                        ! Loop counter
    integer :: k                        ! Loop counter
    integer :: ivmr
    integer :: mif
    integer :: MAF                      ! major frame counter
!    integer :: VMRINST                  ! Instance index
    integer :: INSTANCE                 ! Relevant instance for temperature
    integer :: GPHINST                  ! Relevant instance for GPH
    integer :: NOLAYERS                 ! temp.noSurfs - 1
    integer :: nfine                    ! no of fine resolution grids
    integer :: status                   ! allocation status 
    integer :: SIDEBAND                 ! Loop index
    integer :: SIDEBANDSTART            ! For sideband loop
    integer :: SIDEBANDSTEP             ! For sideband loop
    integer :: SIDEBANDSTOP             ! For sideband loop
    integer :: THISSIDEBAND             ! Loop counter for sidebands
    integer :: NOUSEDCHANNELS           ! How many channels are we considering
    integer :: SIGIND                   ! Signal index, loop counter

    integer :: quantity_type, L_quantity_type       ! added on Jul 13
    integer :: L_state_type

    integer :: iCloudHeight                          ! Index for Cloud Top Height

    integer, dimension(:), pointer :: closestInstances 

    integer :: WHICHPATTERN                           ! Index of antenna pattern
    integer :: MAXSUPERSET                            ! Max. value of superset
    integer, dimension(:), pointer :: SUPERSET        ! Result of AreSignalsSuperset
    integer, dimension(1) :: WHICHPATTERNASARRAY      ! Result of minloc
    integer, dimension(:), pointer :: USEDCHANNELS    ! Which channel is this
    integer, dimension(:), pointer :: USEDSIGNALS     ! Which signal is this channel from

    real(r8), dimension(:,:), pointer :: A_CLEARSKYRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDINDUCEDRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDEXTINCTION
    real(r8), dimension(:,:), pointer :: A_CLOUDRADSENSITIVITY
    real(r8), dimension(:,:), pointer :: A_EFFECTIVEOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_MASSMEANDIAMETER
    real(r8), dimension(:,:), pointer :: A_TOTALEXTINCTION
    real(r8), dimension(:,:), pointer :: VMRARRAY     ! The VMRs

    real (r8), dimension(:), pointer :: thisRatio     ! Sideband ratio values

    real(r8), dimension(:), pointer :: phi_fine       ! Fine resolution for phi 
    real(r8), dimension(:), pointer :: z_fine         ! Fine resolution for z
    real(r8), dimension(:), pointer :: ds_fine        ! Fine resolution for ds
    real(r8) :: ds_tot                                ! Total length of all ds

    real(r8), dimension(:,:), pointer :: A_TRANS
    real(r8), dimension(:), pointer :: FREQUENCIES
    real(r8), dimension(:,:), allocatable :: WC
    real(r8), dimension(:,:), allocatable :: TransOnS
    real(r8) :: phi_tan
    real(r8) :: dz                                    ! thickness of state quantity
    real(r8) :: dphi                                  ! phi interval of state quantity
    real(r8) :: tLat                                  ! temperature 'window' latitude
    real(r8) :: CloudHeight                           ! Cloud Top Height

    logical, dimension(:), pointer :: doChannel       ! Do this channel?
    logical :: DoHighZt                               ! Flag
    logical :: DoLowZt                                ! Flag
    logical :: Got(2)  = .false.  
    logical :: QGot(8) = .false.  
    logical :: dee_bug = .true.  
    logical :: prt_log = .false.
    logical :: FOUNDINFIRST                           ! Flag to indicate derivatives

    character :: cloudtype                            ! cloud profile type

    !---------------------------------------------------------------------------
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>> Executable code  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !---------------------------------------------------------------------------

    if ( toggle(emit) ) call trace_begin ( 'FullCloudForwardModel' )

    !--------------------------
    ! Nullify all the pointers
    !--------------------------

    nullify( CLOUDICE, CLOUDWATER, CLOUDEXTINCTION, CLOUDINDUCEDRADIANCE,    &
             CLOUDRADSENSITIVITY, EFFECTIVEOPTICALDEPTH, GPH,                &
             MASSMEANDIAMETERICE, MASSMEANDIAMETERWATER, PTAN,               &
             RADIANCE, SIZEDISTRIBUTION, EARTHRADIUS, SURFACETYPE,           &
             TEMP, TOTALEXTINCTION, VMR, VMRARRAY,closestInstances,          &
             A_CLEARSKYRADIANCE, A_CLOUDINDUCEDRADIANCE,                     &
             A_CLOUDEXTINCTION, A_CLOUDRADSENSITIVITY,                       &
             A_EFFECTIVEOPTICALDEPTH, A_MASSMEANDIAMETER,                    &
             A_TOTALEXTINCTION, A_TRANS,FREQUENCIES,                         &
             superset, usedchannels, usedsignals, thisRatio,                 &
             JBLOCK, state_ext, state_los )
             
    nullify ( doChannel )
    
    ! ------------------------------------
    ! Find which maf is called at present
    ! ------------------------------------

    maf = fmStat%maf
  
    !------------------------------------------------------------------------
    ! Current version can only have one signal for FullCloudForwardModel
    ! This will be updated soon
    !------------------------------------------------------------------------

    do sigInd = 1, size(forwardModelConfig%signals)

    !    if ( size ( forwardModelConfig%signals ) /= 1 )                        &
    !       & call MLSMessage ( MLSMSG_Error, ModuleName,                       &
    !       & 'Cannot call the full cloud forward model with multiple signals' )

    ! -------------------------------------
    ! Identify the signal (band)
    ! -------------------------------------

    !signal = forwardModelConfig%signals(1)
      
     signal = forwardModelConfig%signals(sigInd)

     if (prt_log) print*,signal%index

    ! --------------------------------------------
    ! Get the quantities we need from the vectors
    ! --------------------------------------------

    ! --------
    ! Outputs:
    ! --------
        radiance => GetVectorQuantityByType ( fwdModelOut,                 &
          & quantityType=l_radiance,                                       &
          & signal=signal%index, sideband=signal%sideband )
        cloudInducedRadiance => GetVectorQuantityByType ( fwdModelOut,     &
          & quantityType=l_cloudInducedRadiance,                           &
          & signal=signal%index, sideband=signal%sideband )
        cloudExtinction => GetVectorQuantityByType ( fwdModelOut,          & 
            & quantityType=l_cloudExtinction, noerror=.true. )
        cloudRADSensitivity => GetVectorQuantityByType ( fwdModelOut,      &
          & quantityType=l_cloudRADSensitivity, noerror=.true.,            &
          & signal=signal%index, sideband=signal%sideband )
        totalExtinction => GetVectorQuantityByType ( fwdModelOut,          &
          & quantityType=l_totalExtinction, noerror=.true. )
        effectiveOpticalDepth => GetVectorQuantityByType ( fwdModelOut,    &
          & quantityType=l_effectiveOpticalDepth, noerror=.true.,          &
          & signal=signal%index, sideband=signal%sideband )
        massMeanDiameterIce => GetVectorQuantityByType ( fwdModelOut,      &
          & quantityType=l_massMeanDiameterIce, noerror=.true. )
        massMeanDiameterWater => GetVectorQuantityByType ( fwdModelOut,    &
          & quantityType=l_massMeanDiameterWater, noerror=.true. )

    ! -------
    ! Inputs:
    ! -------
        ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,       &
          & quantityType=l_ptan, instrumentModule = radiance%template%instrumentModule)
        temp => GetVectorQuantityByType ( fwdModelIn,  fwdModelExtra,      &
          & quantityType=l_temperature )
        gph => GetVectorQuantityByType ( fwdModelIn,  fwdModelExtra,       &
          & quantityType=l_gph )
        cloudIce => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,   &
          & quantityType=l_cloudIce )
        cloudWater => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_cloudWater )
        surfaceType => GetVectorQuantityByType ( fwdModelIn,fwdModelExtra, &
          & quantityType=l_surfaceType )
        sizeDistribution=>GetVectorQuantityByType(fwdModelIn,fwdModelExtra,&
          & quantityType=l_sizeDistribution )
        earthradius=>GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,  &
          & quantityType=l_earthradius ) 
	     scGeocAlt => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_scGeocAlt )
        elevOffset => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_elevOffset, radiometer=Signal%radiometer )	

    !-----------------------------------------
    ! Make sure the quantities we need are got and with correct format
    !-----------------------------------------

        if(.not. (associated(ptan) .and. associated(temp) .and. associated(gph) &
           & .and. associated(cloudIce) .and. associated(cloudWater) &
           & .and. associated(surfaceType) .and. associated(sizeDistribution) &
           & .and. associated(earthradius) .and. associated(scGeocAlt) &
           & .and. associated(elevOffset) ) ) &
          call MLSMessage ( MLSMSG_Error, ModuleName,                        &
                            'Cloud Model inputs are not sufficient')

    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(gph, stacked=.true., coherent=.true.,  &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true.,               &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'ptan' )

    !----------------------------------
    ! Set up some temporary quantities
    !----------------------------------

    call Allocate_test ( closestInstances, radiance%template%noInstances,   &
      & 'closestInstances', ModuleName )      

    !------------------------
    ! Assemble the vmr array
    !------------------------

    if ( size(forwardModelConfig%molecules) .lt. 2 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, 'Not enough molecules' )
    endif
    call allocate_test ( vmrArray,                                           &
      & size(forwardModelConfig%molecules), temp%template%noSurfs,           &
      & 'vmrArray', ModuleName )

    !---------------------------------------------------------
    ! Work out the closest instances from temperature
    !---------------------------------------------------------
    call FindClosestInstances ( temp, radiance, closestInstances )
    instance = closestInstances(maf)


    ivmr=0
    do i = 1, size(forwardModelConfig%molecules)
      select case (forwardModelConfig%molecules(i))
        case(L_H2O)
          ivmr=1
        case(L_O3)
          ivmr=2
        case default
          ivmr=0
      end select
      if(ivmr==0) then
        cycle
      endif

      vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,            &
        & quantityType=l_vmr, molecule=forwardModelConfig%molecules(i) )
      if (.not.ValidateVectorQuantity( vmr, stacked=.true., coherent=.true., &
        & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,  &
        & ModuleName, InvalidQuantity//'vmr' )

!      call FindClosestInstances ( vmr, radiance, closestInstances )
!      vmrInst = closestInstances(maf)
!      instance = closestInstances(maf)

      call InterpolateValues ( &
        & vmr%template%surfs(:,1), &    ! Old X
        & vmr%values(:,instance), &            ! Old Y
        & temp%template%surfs(:,1), &   ! New X
        & vmrArray(i,:), &                     ! New Y
        & 'Linear', extrapolate='Clamp' )
      Got(ivmr)=.true.
    end do

    if ( .not. got(1) .or. .not. got(2) ) then
      call MLSMessage( MLSMSG_Error, ModuleName,                             &
                      'Missing the required molecules' )
    endif

    ! ----------------------------
    ! Get some basic dimensions
    ! ----------------------------
    noChans = radiance%template%noChans
    noMIFs  = radiance%template%noSurfs
    noFreqs = size (signal%frequencies) ! Note: noFreq and noChans should be the same
    noSurf  = temp%template%noSurfs     ! Number of model layers

    call allocate_test ( thisRatio, noChans, 'thisRatio', ModuleName )
    call allocate_test ( doChannel, noFreqs, 'doChannel', ModuleName )

    doChannel = .true.
    if ( associated ( signal%channels ) ) doChannel = signal%channels

    !----------------------------------------------------------
    ! Old version has only single sideband
    ! Now changed to double sideband, Oct 2, 2001
    !----------------------------------------------------------
    ! if ( signal%sideband == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,  &
    !   & 'Only single sidebands allowed in FullForwardCloudModel for now' )
    !----------------------------------------------------------

    !-------------------------------------------------
    ! Allow folded double sideband
    !-------------------------------------------------

    if ( signal%sideband == 0 ) then
       sidebandStart = -1
       sidebandStop  = 1
       sidebandStep  = 2
       sidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
                     & quantityType = l_sidebandRatio, signal=signal%index )
       thisRatio = sidebandRatio%values(:,1)
    else
       sidebandStart = 0
       sidebandStop  = 0
       sidebandStep  = 1
    endif

    !--------------------------------------------
    ! Loop over sidebands 
    !--------------------------------------------

    do sideband = sidebandStart, sidebandStop, sidebandStep

    !-----------------------------
    ! Setup a sideband ratio array
    !------------------------------

    if ( sidebandStart /= sidebandStop ) then
       thisRatio = sidebandRatio%values(:,1)
       if ( sideband == 1 ) thisRatio = 1.0 - thisRatio
    end if

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
    call Allocate_test ( a_massMeanDiameter,                                 &
      & 2, temp%template%noSurfs,                                            &
      & 'a_massMeanDiameter', ModuleName )
    call Allocate_test ( a_Trans,                                            &
      & temp%template%noSurfs, noFreqs,                                      &
      & 'a_trans', ModuleName )
    
    if (noSurf /= GPH%template%nosurfs) then
      call MLSMessage ( MLSMSG_Error, ModuleName,                            &
      & 'number of levels in gph does not match no of levels in temp' )
     else if (radiance%template%nosurfs /= ptan%template%nosurfs) then
        call MLSMessage ( MLSMSG_Error, ModuleName,                          &
        & 'number of levels in radiance does not match no of levels in ptan' )
    endif

    allocate ( WC(2,NOsurf), STAT=status )

    WC (1,:) = CloudIce%values(:,instance)
    WC (2,:) = CloudWater%values(:,instance)

    phi_tan = Deg2Rad * temp%template%phi(1,instance)

    call allocate_test ( superset, size(antennaPatterns), &
         & 'superset', ModuleName )

        do j = 1, size(antennaPatterns)
          superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
           & ForwardModelConfig%signals, sideband=signal%sideband , channel=1 )
        end do

    !    if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
    !           & "No matching antenna patterns." )

    whichPattern = 1
    maxSuperset = maxval ( superset )
    where ( superset < 0 )
          superset = maxSuperset + 1
    end where
    whichPatternAsArray = minloc ( superset )
    whichPattern = whichPatternAsArray(1)
    if ( toggle(emit) .and. levels(emit) > 2 ) then
       call output ( 'Using antenna pattern: ' )
       call output ( whichPattern, advance='yes' )
    end if


    call Allocate_test ( frequencies, noFreqs, 'frequencies', ModuleName )

    frequencies = signal%sideband * ( signal%centerFrequency +   &
                  signal%frequencies)
   
        select case ( sideband )
        case ( -1 )
          frequencies = signal%lo - frequencies
        case ( +1 )
          frequencies = signal%lo + frequencies
        case ( 0 )
          frequencies = signal%lo + frequencies
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Bad value of signal%sideband' )
        end select

        ! frequencies = signal%lo + signal%sideband * ( signal%centerFrequency +   &
        !               signal%frequencies)

    if (prt_log) then
       print*, ' '
       print*,'No. of Frequencies:', noFreqs 
       print*,'Frequency=', frequencies/1e3_r8
    endif

    tLat = temp%template%geodLat(1,instance)

    IF ( present(jacobian) ) THEN

    if (prt_log) print*, 'jacobian is true'

    !----------------------------
    ! get state quantity type 
    !----------------------------
       state_ext => GetVectorQuantityByType(FwdModelIn,quantityType=l_cloudExtinction)
       state_los => GetVectorQuantityByType(FwdModelIn,quantityType=l_LosTransFunc)

       if (tLat .ge. -30._r8 .and. tLat .le. 30._r8) then
          CloudType='Convective'
       else
          CloudType='Frontal'
       endif

       iCloudHeight = 0
       do i = 1, noSurf  
          if(state_ext%values(i,instance) .ne. 0.) then
            iCloudHeight = i                    ! FIND INDEX FOR CLOUD-TOP              
          endif
        enddo

        CloudHeight = 18.e3_r8     ! meters  as a default
        if(iCloudHeight .ne. 0) CloudHeight = gph%values(iCloudHeight, instance)

        call CLOUD_MODEL ( CloudType, CloudHeight, gph%values(:,instance),   &
	     &            noSurf, WC )

    ENDIF

    !------------------------------------------
    ! Now call the full CloudForwardModel code
    !------------------------------------------

    call CloudForwardModel ( doChannel,                                      &
      & noFreqs,                                                             &
      & noSurf,                                                              & 
      & radiance%template%noSurfs,                                           &
      & size(ForwardModelConfig%molecules),                                  &
      & ForwardModelConfig%no_cloud_species,                                 &
      & ForwardModelConfig%no_model_surfs,                                   &
      & frequencies/1e3_r8,                                                  &
      & 10.0**(-temp%template%surfs),                                        &
      & gph%values(:, instance),                                             &
      & temp%values(:,instance),                                             &
      & vmrArray,                                                            &
      & WC,                                                                  &
      & int(sizeDistribution%values(:,instance)),                            &
      & 10.0**(-ptan%values(:,maf)),                                         &
      & earthradius%values(1,maf),                                           &
      & int(surfaceType%values(1, instance)),                                &
      & forwardModelConfig%cloud_der,                                        &
      & forwardModelConfig%cloud_width,                                      &
      & forwardModelConfig%cloud_fov,                                        &
      & phi_tan,                                                             &
      & scGeocAlt%values(1,1),                                               &
      & elevOffset%values(1,1),                                              &
      & antennaPatterns(whichPattern),                                       &
      & a_clearSkyRadiance,                                                  &
      & a_cloudInducedRadiance,                                              &
      & a_trans,                                                             &
      & a_totalExtinction,                                                   &
      & a_cloudExtinction,                                                   &
      & a_massMeanDiameter,                                                  &
      & a_effectiveOpticalDepth,                                             &
      & a_cloudRADSensitivity,                                               &
      & forwardModelConfig%NUM_SCATTERING_ANGLES,                            &  
      & forwardModelConfig%NUM_AZIMUTH_ANGLES,                               &
      & forwardModelConfig%NUM_AB_TERMS,                                     &
      & forwardModelConfig%NUM_SIZE_BINS )

    if (prt_log) print*, 'Successfully done with Full Cloud Foward Model ! '

    deallocate (WC, stat=status)
     
    !------------------------------------------------------------------------
    ! Now store results in relevant vectors
    ! Vectors are stored (noChannels*noSurfaces, noInstances), so transpose
    ! and reshape temporary variables to be in the right form.
    !------------------------------------------------------------------------
    
    !------------------------------
    ! First the minor frame stuff
    !------------------------------

    radiance%values ( :, maf) =                                              &
      & reshape ( transpose(a_clearSkyRadiance),                             &
      & (/radiance%template%instanceLen/) )

    if(associated(cloudInducedRadiance)) &
      & cloudInducedRadiance%values ( :, maf ) =                             &
      & reshape ( transpose(a_cloudInducedRadiance),                         &
      & (/cloudInducedRadiance%template%instanceLen/) )

    if(associated(effectiveOpticalDepth)) &
      & effectiveOpticalDepth%values ( :, maf ) =                            &
      & reshape ( transpose(a_effectiveOpticalDepth),                        &
      & (/effectiveOpticalDepth%template%instanceLen/) )

    if(associated(cloudRADSensitivity)) &
      & cloudRADSensitivity%values ( :, maf ) =                              &
      & reshape ( transpose(a_cloudRADSensitivity),                          &
      & (/cloudRADSensitivity%template%instanceLen/) )

    ! -----------------------------------------------------------------------------
    ! For layer(noTempSurfs-1) stuff make sure all are zero to start, then do rest
    ! -----------------------------------------------------------------------------

    if(associated(cloudExtinction)) &
      & cloudExtinction%values ( :, instance )    = a_cloudExtinction(:,1)
    if(associated(massMeanDiameterIce)) &
      & massMeanDiameterIce%values (:,instance)   = a_massMeanDiameter(1,:)
    if(associated(massMeanDiameterWater)) &
      & massMeanDiameterWater%values(:, instance) =  a_massMeanDiameter(2,:)
    if(associated(totalExtinction)) &
      & totalExtinction%values ( :, instance )    = a_totalExtinction (:,1)

    !-----------------------
    ! Start output Jacobian
    !-----------------------
    IF ( present(jacobian) ) THEN

    !-------------------------------------------------------
    ! Determine which retrieval is to be used
    !-------------------------------------------------------
      noInstances = state_ext%template%noInstances

    !--------------------------------------------
    ! Jacobian for high tangent height retrieval
    !--------------------------------------------
    !doHighZt = present(jacobian) .and. (l_stateQ_type == l_cloudExtinction)

    doHighZt = .false.
    if (doHighZt) then      

    ! Get dimension
      noSurf = state_ext%template%noSurfs
      
      colJBlock = FindBlock ( Jacobian%col, state_ext%index, maf )
      rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
      fmStat%rows(rowJBlock) = .true.
!      colJBlock = 0
!      do while (colJBlock <= jacobian%col%nb .and. &
!           jacobian%col%inst(colJBlock) /= maf)
!           colJBlock = colJBlock +1 
!      end do

      jBlock => jacobian%block(rowJblock,colJblock)

      select case ( jBlock%kind )
      case ( M_Absent )
        call CreateBlock ( jBlock, noMIFs, noSurf*noInstances, M_Full )
        jBlock%values = 0.0_r8

      !-------------------------------------------------------------------
      ! we use 100 times better resolution to compute weighting functions
      !-------------------------------------------------------------------
      nfine=100
      allocate( phi_fine(nfine*noInstances), stat=status )
      allocate( z_fine(nfine*noInstances), stat=status )
      allocate( ds_fine(nfine*noInstances), stat=status )
      
      do i=1,noInstances*nfine
       phi_fine(i) = minval(state_ext%template%phi(1,:)) + &
         & 1._r8*(i-1)/nfine/noInstances * &
         & (maxval(state_ext%template%phi(1,:)) - minval(state_ext%template%phi(1,:)))
      end do 
      
      do mif = 1, noMIFs
        
        !----------------------------------
        ! find intervals of stateQ at maf
        !----------------------------------
        dz = abs(state_ext%template%surfs(2,maf)-state_ext%template%surfs(1,maf))
        dphi = abs(state_ext%template%phi(2,maf)-state_ext%template%phi(1,maf))

        !------------------------------
        ! find z for given phi_fine
        !------------------------------
        z_fine = (earthradius%values(1,maf)+ &
         & (ptan%values(mif,maf)+3.)*16.) / &
         & cos((phi_fine - radiance%template%phi(mif,maf))*pi/180._r8) - &
         & earthradius%values(1,maf)
        z_fine = z_fine/16._r8 - 3._r8    ! convert back to log pressure

        !--------------------------------
        ! find ds for each (z,phi) pair
        !--------------------------------
        ds_fine = (earthradius%values(1,maf)+ &
         & (ptan%values(mif,maf)+3.)*16.) / &
         & cos((phi_fine - radiance%template%phi(mif,maf))*pi/180._r8)**2

        !-----------------------------------------------
        ! find the total length for this tangent height
        !-----------------------------------------------
        ds_tot = sum(ds_fine)

        !----------------------------------------------------------
        ! determine weights by the length inside each state domain
        !----------------------------------------------------------
         do i = 1,noInstances             ! loop over profile
         do j = 1,noSurf                  ! loop over surface
         do k = 1, nfine*noInstances      ! sum up all the lengths
           if(abs(z_fine(k) - state_ext%template%surfs(j,i)) < dz/2. &
           & .AND. abs(phi_fine(k) - state_ext%template%phi(j,i)) < dphi/2.) &
           & jBlock%values(mif,j+(i-1)*noInstances) = &
           & jBlock%values(mif,j+(i-1)*noInstances) + ds_fine(k)/ds_tot
         end do
         end do
         end do 
      end do
      
      Deallocate( phi_fine, stat=status )
      Deallocate( z_fine, stat=status )
      Deallocate( ds_fine, stat=status )

      case default
         call MLSMessage(MLSMSG_Error, ModuleName, &
                 & "Invalid matrix block kind in CreateBlock")              
      end select
    end if

    !--------------------------------------------
    ! Jacobian for low tangent height retrieval
    !--------------------------------------------
    ! doLowZt = present(jacobian) .and. (l_stateQ_type == l_LosTransFunc)

    doLowZt = .true.
    if (doLowZt) then
    
        ! Get dimension
        noSgrid=state_los%template%noChans

        colJBlock = FindBlock ( Jacobian%col, state_los%index, maf )
        rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
        fmStat%rows(rowJBlock) = .true.

        jBlock => jacobian%block(rowJblock,colJblock)

        select case ( jBlock%kind )

        case ( M_Absent )

        ! ---------------------------------------------------------------
        ! In the absent case, the jacobian is stored in a special format
        !----------------------------------------------------------------

        call CreateBlock ( jBlock, noChans, &
            & noSgrid*noMIFs, M_Full )
        jBlock%values = 0.0_r8

        allocate( TransOnS(noSgrid, noMIFs), stat=status )

        !------------------------------------------
        ! Now fill the jacobian, case ( M_Absent )
        !------------------------------------------
        do chan = 1, noChans
          if ( doChannel(chan) ) then
          !-----------------------------------------------------------------
          ! now, we define beta as transmission function in Sensitivity.f90
          ! and we interpolate it onto Sgrid
          !-----------------------------------------------------------------
          call FindTransForSgrid (                                   &
                      &     ptan%values(:,maf),                      &
                      &     earthradius%values(1,maf)*1.e-3_r8,      &
                      &     noMIFs,                                  &
                      &     temp%template%noSurfs,                   &
                      &     noSgrid,                                 &
                      &     gph%template%Surfs,                      &
                      &     a_trans(:,chan),               &
                      &     state_los%template%frequencies,             &
                      &     TRANSonS )                

                    do mif = 1, noMIFs
                    do i=1,noSgrid
                     jBlock%values(chan, & 
                       & i+(mif-1)*noSgrid)= TransOnS(i,mif)
                    end do
                    end do
          end if
        end do

        Deallocate(TransOnS,stat=status)

        case ( M_Banded )     
        
        noNonZero = noSgrid*noMIFs*noInstances
        
        call CreateBlock ( jBlock, noChans*noMIFs, &
                         & noSgrid*noMIFs*noInstances, &
                         & M_Banded, noNonZero )

        allocate( TransOnS(noSgrid, noMIFs), stat=status )

       !------------------------------------------
       ! Now fill the jacobian, case ( M_Banded )
       !------------------------------------------
        do chan = 1, noChans
          if ( doChannel(chan) ) then
          !-----------------------------------------------------------------
          ! now, we define beta as transmission function in Sensitivity.f90
          ! and we interpolate it onto Sgrid
          !-----------------------------------------------------------------
          call FindTransForSgrid (                                   &
                      &     ptan%values(:,maf),                      &
                      &     earthradius%values(1,maf)*1.e-3_r8,      &
                      &     noMIFs,                                  &
                      &     temp%template%noSurfs,                   &
                      &     noSgrid,                                 &
                      &     gph%template%Surfs,                      &
                      &     a_trans(:,chan),                         &
                      &     state_los%template%frequencies,             &
                      &     TRANSonS )                

            jBlock%values = 0.0_r8
            jBlock%r2(0) = 0
            
            do i=1,noSgrid
            do mif=1,noMIFs
            jBlock%r1(i+(mif-1)*noSgrid) = noChans*(mif-1)
            jBlock%r2(i+(mif-1)*noSgrid) = noChans* & 
               & (i+(mif-1)*noSgrid+(maf-1)*noInstances)
            jBlock%values(i+(mif-1)*noSgrid,1) = TransOnS(i,mif)
            enddo
            enddo

         end if
         end do

              Deallocate(TransOnS,stat=status)
         case default
               call MLSMessage(MLSMSG_Error, ModuleName, &
                 & "Invalid matrix block kind in CreateBlock")              
         end select
               
    endif

    ENDIF   
    !------------------------------
    ! End of output jacobian
    !------------------------------

    !------------------------------
    ! Remove temporary quantities
    !------------------------------
    call deallocate_test ( superset, 'superset',          ModuleName )

    call Deallocate_test ( a_massMeanDiameter,                               &
                          'a_massMeanDiameter',           ModuleName )
    call Deallocate_test ( a_cloudExtinction,                                &
                          'a_cloudExtinction',            ModuleName )
    call Deallocate_test ( a_trans,                                          &
                          'a_trans',                      ModuleName )
    call Deallocate_test ( a_totalExtinction,                                &
                          'a_totalExtinction',            ModuleName )
    call Deallocate_test ( a_cloudRADSensitivity,                            &
                          'a_cloudRADSensitivity',        ModuleName )
    call Deallocate_test ( a_effectiveOpticalDepth,                          &
                          'a_effectiveOpticalDepth',      ModuleName )
    call Deallocate_test ( a_cloudInducedRadiance,                           &
                          'a_cloudInducedRadiance',       ModuleName )
    call Deallocate_test ( a_clearSkyRadiance,                               &
                          'a_clearSkyRadiance',           ModuleName )
    call Deallocate_test ( vmrArray,                                         &
                          'vmrArray',                     ModuleName )
    call Deallocate_test ( closestInstances,                                 &
                          'closestInstances',             ModuleName )

    call Deallocate_test ( doChannel, 'doChannel',        ModuleName )

    !--------------------------------------------
    ! End of sideband loop 
    !--------------------------------------------
    enddo

    if (prt_log) then
      print*, ' '
      print*, 'Time Instance: ', instance
    endif

    end do  ! End of signals
    
    if ( toggle(emit) ) call trace_end ( 'FullCloudForwardModel' )

    if(prt_log) print*, 'Successful done with full cloud forward wapper !'

  end subroutine FullCloudForwardModelWrapper

end module FullCloudForwardModel


subroutine FindTransForSgrid ( PT, Re, NT, NZ, NS, Zlevel, TRANSonZ, Slevel, TRANSonS)

  use MLSCommon,only: r8
  use MLSNumerics, only: INTERPOLATEVALUES
 
  integer :: NT                 ! No. of tangent pressures
  integer :: NZ                 ! No. of pressure levels
  integer :: NS                 ! No. of S levels
  integer :: mif
  
  real(r8) :: PT(NT)            ! this -log10 pressure
  real(r8) :: Re                ! in km
  real(r8) :: Zlevel(NZ)
  real(r8) :: Slevel(NS)
  real(r8) :: TRANSonZ(NZ)
  real(r8) :: TRANSonS(NS,NT)

  real(r8) :: zt(nt),x_out(ns)
    TransOns = 0._r8
    zt = (pt+3.)*16.                      ! converted to height in km
    do mif=1,nt

      ! find altitude of each s grid
      x_out = Slevel**2/2./(re + zt(mif))
      ! converted to zeta
      x_out = x_out/16. + pt(mif)

      CALL INTERPOLATEVALUES(Zlevel,TransOnZ,x_out,TransOnS(:,mif),method='Linear')

     enddo

end subroutine FindTransForSgrid

! $Log$
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
