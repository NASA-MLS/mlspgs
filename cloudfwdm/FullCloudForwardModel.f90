
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
  use CloudySkyModule, only: CLOUD_MODEL
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
    type (VectorValue_T), pointer :: modelCLOUDRADIANCE         ! modelled cloud radiance
    type (VectorValue_T), pointer :: obsCLOUDRADIANCE           ! observed cloud radiance
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
    type(VectorValue_T), pointer :: STATE_ext      ! A state vector quantity
    type(VectorValue_T), pointer :: STATE_los      ! A state vector quantity

    type(VectorValue_T), pointer :: SIDEBANDRATIO  ! From the state vector

    ! for jacobian
    type(MatrixElement_T), pointer :: JBLOCK       ! A block from the jacobian
    integer :: COLJBLOCK                ! Column index in jacobian
    integer :: ROWJBLOCK                ! Row index in jacobian
    integer :: noInstances              ! no of instance
    integer :: noMIFs                   ! Number of minor frames
    integer :: noSgrid                  ! no of elements in S grid
    integer :: noSurf                   ! Number of pressure levels
    integer :: novmrSurf                   ! Number of vmr levels
    integer :: noCldSurf                   ! Number of cloud ext levels
    integer :: NOFREQS                  ! Number of frequencies

    integer :: i                        ! Loop counter
    integer :: j                        ! Loop counter
    integer :: k                        ! Loop counter
    integer :: ivmr
    integer :: mif
    integer :: MAF                      ! major frame counter
    integer :: INSTANCE                 ! Relevant instance for temperature
    integer :: MinInst                  ! lower bound of instance
    integer :: MaxInst                  ! upper bound of instance
    integer :: nfine                    ! no of fine resolution grids
    integer :: nNear                    ! no of nearest profiles
    integer :: status                   ! allocation status 
    integer :: SIDEBAND                 ! Loop index
    integer :: SIDEBANDSTART            ! For sideband loop
    integer :: SIDEBANDSTEP             ! For sideband loop
    integer :: SIDEBANDSTOP             ! For sideband loop
    integer :: THISSIDEBAND             ! Loop counter for sidebands
    integer :: SIGIND                   ! Signal index, loop counter

    integer :: iCloudHeight                          ! Index for Cloud Top Height

    integer, dimension(:), pointer :: closestInstances 

    integer :: WHICHChannel                           ! which single channel is used
    integer :: WHICHPATTERN                           ! Index of antenna pattern
    integer :: MAXSUPERSET                            ! Max. value of superset
    integer, dimension(:), pointer :: SUPERSET        ! Result of AreSignalsSuperset
    integer, dimension(1) :: WHICHPATTERNASARRAY      ! Result of minloc

    real(r8), dimension(:,:), pointer :: A_CLEARSKYRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDINDUCEDRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDEXTINCTION
    real(r8), dimension(:,:), pointer :: A_CLOUDRADSENSITIVITY
    real(r8), dimension(:,:), pointer :: A_EFFECTIVEOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_MASSMEANDIAMETER
    real(r8), dimension(:,:), pointer :: A_TOTALEXTINCTION
    real(r8), dimension(:,:), pointer :: VMRARRAY     ! The VMRs
    real(r8), dimension(:), allocatable :: Slevl      ! S grid
    real(r8), dimension(:), allocatable :: Zt         ! tangent height
    real (r8), dimension(:), pointer :: thisRatio     ! Sideband ratio values

    real(r8), dimension(:), pointer :: phi_fine       ! Fine resolution for phi 
    real(r8), dimension(:), pointer :: z_fine         ! Fine resolution for z
    real(r8), dimension(:), pointer :: zp_fine         ! Fine resolution for zp
    real(r8), dimension(:), pointer :: s_fine        ! Fine resolution for s
    real(r8), dimension(:), pointer :: ds_fine        ! Fine resolution for ds
    real(r8), dimension(:), pointer :: w_fine        ! weight along s_fine

    real(r8), dimension(:,:,:), allocatable  :: A_TRANS
    real(r8), dimension(:), pointer :: FREQUENCIES
    real(r8), dimension(:,:), allocatable :: WC

    real(r8) :: phi_tan
    real(r8) :: dz                                    ! thickness of state quantity
    real(r8) :: dphi                                  ! phi interval of state quantity
    real(r8) :: tLat                                  ! temperature 'window' latitude
    real(r8) :: CloudHeight                           ! Cloud Top Height

    logical, dimension(:), pointer :: doChannel       ! Do this channel?
    logical :: DoHighZt                               ! Flag
    logical :: DoLowZt                                ! Flag
    logical :: Got(2)  = .false.  
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

    nullify( CLOUDICE, CLOUDWATER, CLOUDEXTINCTION, modelCLOUDRADIANCE,    &
             CLOUDRADSENSITIVITY, EFFECTIVEOPTICALDEPTH, obsCLOUDRADIANCE,GPH, &
             MASSMEANDIAMETERICE, MASSMEANDIAMETERWATER, PTAN,               &
             RADIANCE, SIZEDISTRIBUTION, EARTHRADIUS, SURFACETYPE,           &
             TEMP, TOTALEXTINCTION, VMR, VMRARRAY,closestInstances,          &
             A_CLEARSKYRADIANCE, A_CLOUDINDUCEDRADIANCE,                     &
             A_CLOUDEXTINCTION, A_CLOUDRADSENSITIVITY,                       &
             A_EFFECTIVEOPTICALDEPTH, A_MASSMEANDIAMETER,                    &
             A_TOTALEXTINCTION,FREQUENCIES,                         &
             superset, thisRatio, JBLOCK, state_ext, state_los )
             
    nullify ( doChannel )
    
    ! ------------------------------------
    ! Find which maf is called at present
    ! ------------------------------------

    maf = fmStat%maf
  
    !-------------------------------------------------------
    ! Determine which retrieval is to be used
    !    cloud_der = 0     high tangent height
    !    cloud_der = 1     low tangent height
    !-------------------------------------------------------
    doHighZt = .false.
    doLowZt = .false.
    if(forwardModelConfig%cloud_der == 0) doHighZt = .true.
    if(forwardModelConfig%cloud_der == 1) doLowZt = .true.

    !------------------------------------------------------------------------
    ! Current version can only have one signal for FullCloudForwardModel
    ! This will be updated soon
    !------------------------------------------------------------------------

    do sigInd = 1, size(forwardModelConfig%signals)

    ! -------------------------------------
    ! Identify the signal (band)
    ! -------------------------------------
      
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
        modelCloudRadiance => GetVectorQuantityByType ( fwdModelOut,     &
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
        ptan => GetVectorQuantityByType ( fwdModelExtra,       &
          & quantityType=l_ptan, instrumentModule = radiance%template%instrumentModule)
        temp => GetVectorQuantityByType ( fwdModelExtra,      &
          & quantityType=l_temperature )
        gph => GetVectorQuantityByType ( fwdModelExtra,       &
          & quantityType=l_gph )
        obsCloudRadiance => GetVectorQuantityByType ( fwdModelExtra,     &
          & quantityType=l_cloudInducedRadiance, noerror=.true.,  &
          & signal=signal%index, sideband=signal%sideband )
        cloudIce => GetVectorQuantityByType ( fwdModelExtra,   &
          & quantityType=l_cloudIce )
        cloudWater => GetVectorQuantityByType ( fwdModelExtra, &
          & quantityType=l_cloudWater )
        surfaceType => GetVectorQuantityByType ( fwdModelExtra, &
          & quantityType=l_surfaceType )
        sizeDistribution=>GetVectorQuantityByType(fwdModelExtra,&
          & quantityType=l_sizeDistribution )
        earthradius=>GetVectorQuantityByType ( fwdModelExtra,  &
          & quantityType=l_earthradius ) 
        scGeocAlt => GetVectorQuantityByType ( fwdModelExtra, &
          & quantityType=l_scGeocAlt )
        elevOffset => GetVectorQuantityByType ( fwdModelExtra, &
          & quantityType=l_elevOffset, radiometer=Signal%radiometer )	

    !-----------------------------------------
    ! Make sure the quantities we need are got and with correct format
    !-----------------------------------------

    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(gph, stacked=.true., coherent=.true.,  &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true.,               &
       & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,   &
       & ModuleName, InvalidQuantity//'ptan' )

    ! ----------------------------
    ! Get some basic dimensions
    ! ----------------------------
    noMIFs  = radiance%template%noSurfs
    noFreqs = size (signal%frequencies) ! Note: noFreq and noChans should be the same
    noSurf  = temp%template%noSurfs     ! Number of model layers

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
      & size(forwardModelConfig%molecules), noSurf,           &
      & 'vmrArray', ModuleName )

    !---------------------------------------------------------
    ! Work out the closest instances from temperature
    !---------------------------------------------------------
    call FindClosestInstances ( temp, radiance, closestInstances )
    instance = closestInstances(maf)
    tLat = temp%template%geodLat(1,instance)    ! get latitude for each instance
!   print*,'Lat=',tLat


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

	   novmrSurf = vmr%template%nosurfs
      call InterpolateValues ( &
        & reshape(vmr%template%surfs(:,1),(/novmrSurf/)), &    ! Old X
        & reshape(vmr%values(:,instance),(/novmrSurf/)), &     ! Old Y
        & reshape(temp%template%surfs(:,1),(/noSurf/)), &   ! New X
        & vmrArray(i,:), &                     ! New Y
        & 'Linear', extrapolate='Clamp' )
      Got(ivmr)=.true.
    end do

    if ( .not. got(1) .or. .not. got(2) ) then
      call MLSMessage( MLSMSG_Error, ModuleName,                             &
                      'Missing the required molecules' )
    endif

    call allocate_test ( thisRatio, noFreqs, 'thisRatio', ModuleName )
    call allocate_test ( doChannel, noFreqs, 'doChannel', ModuleName )
    allocate(zt(noMifs))

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
       sidebandRatio => GetVectorQuantityByType ( fwdModelExtra, &
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


    if (prt_log) print*, 'jacobian is true'

    !----------------------------
    ! get state quantity type (need both ext and los quantities for Jacobians)
    !----------------------------
       state_ext => GetVectorQuantityByType(FwdModelIn,quantityType=l_cloudExtinction)
       state_los => GetVectorQuantityByType(FwdModelIn,quantityType=l_LosTransFunc)

        ! Get number of cloud surfaces for retrieval
        noCldSurf=state_ext%template%noSurfs
        ! Get s dimension
        noSgrid=state_los%template%noChans
        
        Allocate( a_Trans(noSgrid, noMIFs, noFreqs))
        Allocate( Slevl(noSgrid))

        Slevl = state_los%template%frequencies

     ! get tangent height from tangent pressure
      call InterpolateValues ( &
        & reshape(temp%template%surfs,(/noSurf/)), &    ! Old X
        & reshape(gph%values(:,instance),(/noSurf/)), &            ! Old Y
        & reshape(ptan%values(:,maf),(/noMifs/)), &   ! New X
        & zt, &                     ! New Y
        & 'Linear' )

     IF ( present(jacobian) ) THEN
      ! transmission functions for retrieval are calculated using artificial
      ! cloud profiles depending on where the measurement is taken in latitude
       if (tLat .ge. -30._r8 .and. tLat .le. 30._r8) then
          CloudType='Convective'
       else
          CloudType='Frontal'
       endif

      ! find cloud top index from observed Tcir, threshold to be determined
      !     for high Zt, use Tcir(maf)
      !     for low Zt, use Tcir(maf-2)
      if(.not. associated(obsCloudRadiance)) then
         call MLSMessage( MLSMSG_Error, ModuleName,                             &
                      'Need cloud radiances to estimate cloud top in retrieval' )
      else
       iCloudHeight = 0
       if(doHighZt) then
         do i = 1, noMifs
            if(obsCloudRadiance%values(i,maf) .ne. 0.0_r8) iCloudHeight = i
         enddo
       else
         do i = 1, noMifs
            if(obsCloudRadiance%values(i,maf) .ne. 0.0_r8) iCloudHeight = i
         enddo
       end if

       ! if no cloud is found, cloud top is 18 km
       CloudHeight = 18.e3_r8     ! meters
       if(iCloudHeight .ne. 0) CloudHeight = zt(iCloudHeight)

      ! set up artificial cloud profile for retrieval use only
       call CLOUD_MODEL (CloudType, CloudHeight, gph%values(:,instance), noSurf, WC)
      end if
      
    ENDIF

    !------------------------------------------
    ! Now call the full CloudForwardModel code
    !------------------------------------------

    call CloudForwardModel ( doChannel,                                      &
      & noFreqs,                                                             &
      & noSurf,                                                              & 
      & noMifs,                                           &
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
      & zt, &
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
      & forwardModelConfig%NUM_SIZE_BINS,                                    &
      & Slevl*1000._r8, noSgrid )

    if (prt_log) print*, 'Successfully done with Full Cloud Foward Model ! '

    ! Output minor frame quantities if they are asked: check channel and type
    do i =1 , noFreqs
    if (doChannel(i)) then
    do mif=1, noMIFs
    radiance%values (i+(mif-1)*noFreqs, maf) = a_clearSkyRadiance(mif,i) 

    modelCloudRadiance%values (i+(mif-1)*noFreqs, maf ) = &
      & a_cloudInducedRadiance(mif,i)

    if(associated(effectiveOpticalDepth)) &
      & effectiveOpticalDepth%values (i+(mif-1)*noFreqs, maf ) = &
      & a_effectiveOpticalDepth(mif,i)

    if(associated(cloudRADSensitivity)) &
      & cloudRADSensitivity%values (i+(mif-1)*noFreqs, maf ) =   &
      & a_cloudRADSensitivity(mif,i)
    enddo
    endif
    enddo

    ! -----------------------------------------------------------------------------
    ! output L2GP quantities
    ! -----------------------------------------------------------------------------

    ! To save disk space, we output one channel per band (first true doChannel)
    ! in the last signal (band) will over write all the previous signals (bands)
    FOUNDINFIRST = .true.
    do i=1, noFreqs
    if( doChannel(i) .and. FOUNDINFIRST ) then 
      FOUNDINFIRST = .false.
      if(associated(cloudExtinction)) &
        & cloudExtinction%values (:, instance )    = a_cloudExtinction(:,i)
      if(associated(totalExtinction)) &
        & totalExtinction%values (:, instance )    = a_totalExtinction (:,i)
    endif
    enddo

    if(associated(massMeanDiameterIce)) &
      & massMeanDiameterIce%values (:,instance)   = a_massMeanDiameter(1,:)
    if(associated(massMeanDiameterWater)) &
      & massMeanDiameterWater%values(:, instance) =  a_massMeanDiameter(2,:)


    !-----------------------
    ! Start output Jacobian
    !-----------------------

      noInstances = state_ext%template%noInstances

    !--------------------------------------------
    ! Jacobian for high tangent height retrieval (only one frequency)
    !--------------------------------------------

    if (doHighZt .and. present(jacobian)) then      
    ! do not handle multiple signals and channels  
      if(size(forwardModelConfig%signals) > 1 .and. size(doChannel) > 1) then
         print*,'only one frequency and one signal is allowed'
         stop
      end if
      do i=1, noFreqs
        if (doChannel(i)) whichChannel=i
      end do

      colJBlock = FindBlock ( Jacobian%col, state_ext%index, maf)
      rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
      fmStat%rows(rowJBlock) = .true.

      jBlock => jacobian%block(rowJblock,colJblock)

        call CreateBlock ( jBlock, noMIFs, noCldSurf*noInstances, M_Full )
        jBlock%values = 0._r8

      !-------------------------------------------------------------------
      ! we use 100 times better resolution to compute weighting functions
      !-------------------------------------------------------------------
      nfine = 100
      nNear = ForwardModelConfig%phiWindow         ! default = 5
      ! only nearest instances are mattered
        minInst = instance - (nNear-1)/2
        maxInst = instance + (nNear-1)/2
         if(minInst < 1) minInst = 1
         if(maxInst > noInstances) maxInst = noInstances
         
      allocate( phi_fine(nfine*nNear), stat=status )
      allocate( z_fine(nfine*nNear), stat=status )
      allocate( zp_fine(nfine*nNear), stat=status )
      allocate( s_fine(nfine*nNear), stat=status )
      allocate( w_fine(nfine*nNear), stat=status )
      allocate( ds_fine(nfine*nNear), stat=status )
      
      do i=1,nNear*nfine
       phi_fine(i) = minval(state_ext%template%phi(1,minInst:maxInst)) + &
         & (i-1._r8)/nfine/nNear*(maxval(state_ext%template%phi(1,minInst:maxInst)) - &
         & minval(state_ext%template%phi(1,minInst:maxInst)))
      end do

      !----------------------------------
      ! find vertical and horizontal intervals of stateQ at maf
      !----------------------------------
      ! vertical
      dz = abs(state_ext%template%surfs(2,1)-state_ext%template%surfs(1,1))
      ! horiozontal: in case of the last instance, use the previous one
      if(instance < noInstances) &
         & dphi = abs(state_ext%template%phi(1,instance+1)-state_ext%template%phi(1,instance))
      
      do mif = 1, noMIFs
      ! Jacobians at only tangent heights less than the top level of retrieval are calculated
      if(gph%values(noCldSurf,instance) > zt(mif)) then
        !------------------------------
        ! find z for given phi_fine
        !------------------------------
        z_fine = (earthradius%values(1,maf)+ zt(mif)) / &
         & cos((phi_fine - radiance%template%phi(mif,maf))*Deg2Rad) - &
         & earthradius%values(1,maf)
         ! convert back to log tangent pressure
         call InterpolateValues ( &
            & reshape(gph%values(:,instance),(/noSurf/)), &            ! Old X
            & reshape(temp%template%surfs,(/noSurf/)), &    ! Old Y
            & z_fine, &   ! New X
            & zp_fine, &                     ! New Y
            & 'Linear' )

        !--------------------------------
        ! find ds and weight for each (z,phi) pair
        !--------------------------------
        s_fine = (earthradius%values(1,maf)+ zt(mif)) * &
         & sin((phi_fine - radiance%template%phi(mif,maf))*Deg2Rad)
        ds_fine = 0._r8    ! initialize it
        ds_fine(1:nfine*nNear-1) = s_fine(2:nfine*nNear) - &
         & s_fine(1:nfine*nNear-1)

         call InterpolateValues ( &
            & sLevl, &            ! Old X
            & reshape(a_trans(:,mif,whichChannel),(/noSgrid/)), &    ! Old Y
            & s_fine, &   ! New X
            & w_fine, &                     ! New Y
            & 'Linear' )
            
        ! ds needs to be weighted by transmission function
        
        ds_fine = ds_fine*w_fine    ! ds is in meters

        !----------------------------------------------------------
        ! determine weights by the length inside each state domain
        !----------------------------------------------------------
        
         do i = minInst,maxInst             ! loop over closer profiles
         do j = 1,noCldSurf               ! loop over cloudQty surface
         do k = 1,nfine*nNear             ! sum up all the lengths
           if(abs(zp_fine(k) - state_ext%template%surfs(j,1)) < dz/2._r8 &
           & .AND. abs(phi_fine(k) - state_ext%template%phi(1,i)) < dphi/2._r8) &
           & jBlock%values(mif,j+(i-1)*noCldSurf) = &
           & jBlock%values(mif,j+(i-1)*noCldSurf) + ds_fine(k)/1.e3_r8
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
      
    end if     ! doHighZt

    !--------------------------------------------
    ! Jacobian for low tangent height retrieval
    !--------------------------------------------
    if (doLowZt .and. present(jacobian)) then
    
        colJBlock = FindBlock ( Jacobian%col, state_los%index, maf )
        rowJBlock = FindBlock ( jacobian%row, radiance%index, maf)
        fmStat%rows(rowJBlock) = .true.

        jBlock => jacobian%block(rowJblock,colJblock)

      ! to save space, the jacobian is packed in a full rectangle matrix
        call CreateBlock ( jBlock, noFreqs, noSgrid*noMIFs, M_Full )
        jBlock%values = 0.0_r8

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

    endif      ! doLowZt

    !------------------------------
    ! End of output jacobian
    !------------------------------

    !------------------------------
    ! Remove temporary quantities
    !------------------------------
    Deallocate (WC, a_trans, Slevl, Zt )
     
    call deallocate_test ( superset, 'superset',          ModuleName )
    call Deallocate_test ( a_massMeanDiameter,'a_massMeanDiameter',ModuleName )
    call Deallocate_test ( a_cloudExtinction,'a_cloudExtinction',ModuleName )
    call Deallocate_test ( a_totalExtinction,'a_totalExtinction',ModuleName )
    call Deallocate_test ( a_cloudRADSensitivity,'a_cloudRADSensitivity',ModuleName )
    call Deallocate_test ( a_effectiveOpticalDepth,'a_effectiveOpticalDepth',ModuleName )
    call Deallocate_test ( a_cloudInducedRadiance,'a_cloudInducedRadiance',ModuleName )
    call Deallocate_test ( a_clearSkyRadiance,'a_clearSkyRadiance',ModuleName )
    call Deallocate_test ( vmrArray,'vmrArray',ModuleName )
    call Deallocate_test ( closestInstances,'closestInstances',ModuleName )
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


! $Log$
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
