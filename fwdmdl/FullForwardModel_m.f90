! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullForwardModel_m

  use GLNP, only: NG
  use MLSCommon, only: I4, R4, R8
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
  use ELLIPSE_M, only: ELLIPSE
  use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
  use GET_PATH_SPSFUNC_NGRID_M, only: GET_PATH_SPSFUNC_NGRID
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
      PATH_DERIVATIVE, PATH_VECTOR_2D, PATH_INT_VECTOR_2D
  use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_BETA_PATH_M, only: GET_COARSE_BETA_PATH, GET_GLBETA_PATH
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use PATH_CONTRIB_M, only: PATH_CONTRIB
  use COARSE_DELTA_M, only: COARSE_DELTA
  use RAD_TRAN_M, only: RAD_TRAN
  use RAD_TRAN_WD_M, only: RAD_TRAN_WD
  use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
  use FREQ_AVG_M, only: FREQ_AVG
  use CONVOLVE_ALL_M, only: CONVOLVE_ALL
  use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
  use D_LINTRP_M, only: LINTRP
  use D_HUNT_M, only: hunt_zvi => HUNT
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VALIDATEVECTORQUANTITY, &
                       &   GETVECTORQUANTITYBYTYPE
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, &
                                  &   FORWARDMODELSTATUS_T
  use AntennaPatterns_m, only: ANTENNAPATTERNS
  use FilterShapes_m, only: FILTERSHAPES
  use PointingGrid_m, only: POINTINGGRIDS
  use MatrixModule_1, only: MATRIX_T
  use Trace_M, only: Trace_begin, Trace_end
  use MLSSignals_m, only: SIGNAL_T, MATCHSIGNAL, ARESIGNALSSUPERSET
  use SpectroscopyCatalog_m, only: CATALOG_T, LINES, CATALOG
  use Intrinsic, only: L_TEMPERATURE, L_RADIANCE, L_PTAN, L_ELEVOFFSET, &
    & L_ORBITINCLINATION, L_SPACERADIANCE, L_EARTHREFL, L_LOSVEL,       &
    & L_SCGEOCALT, L_SIDEBANDRATIO, L_NONE, L_CHANNEL, L_VMR, L_REFGPH
  use Units, only: Deg2Rad
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error
  use MLSNumerics, only: HUNT
  use Toggles, only: Emit, Gen, Levels, Switches, Toggle
  use Molecules, only: spec_tags
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Output_m, only: OUTPUT
  use ManipulateVectorQuantities, only: FindClosestInstances

  ! This module contains the `full' forward model.

  implicit none
  private
  public :: FullForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! ================================ FullForwardModel routine ======

  ! -----------------------------------------------  ForwardModel  -----
  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    &                       FwdModelOut, Ifm, FmStat, Jacobian )


    ! Dummy arguments --------------------------------------------------------

    ! From ForwardModelSetup
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: Ifm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Local parameters ---------------------------------------------------------

    character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

    ! Local variables ----------------------------------------------------------

    ! First the old stuff which we hope to get rid of or redefine
    integer(i4) :: brkpt, ch, frq_i, i, ier, ihi, ilo, j, k, lmax, m, maf, &
      max_phi_dim, max_zeta_dim, mid, n, no_ele, no_tan_hts, &
      ptg_i, si, Spectag

    type(path_derivative) :: K_temp_frq
    type(path_derivative), allocatable, dimension(:) :: K_atmos_frq

    type(path_beta), dimension(:,:), pointer :: Beta_path

    real(r8) :: Frq, Geod_lat, H_tan, Phi_tan, Rad

    ! This is the `legit stuff' we hope will stay; they are all pointers to
    ! VectorValue_T's containing vector quantities.
    type (VectorValue_T), pointer :: EARTHREFL     ! Earth reflectivity
    type (VectorValue_T), pointer :: ELEVOFFSET    ! Elevation offset quantity
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: LOSVEL        ! Line of sight velocity
    type (VectorValue_T), pointer :: ORBINCLINE    ! Orbital inclination (beta)
    type (VectorValue_T), pointer :: PTAN          ! PTAN quantity
    type (VectorValue_T), pointer :: FIRSTRADIANCE ! One radiance quantity to be filled
    type (VectorValue_T), pointer :: THISRADIANCE ! One radiance quantity to be filled
    type (VectorValue_T), pointer :: REFGPH        ! Reference GPH, (zRef and hRef)
    type (VectorValue_T), pointer :: SCGEOCALT     ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: SIDEBANDRATIO ! Sideband ratio for radiance
    type (VectorValue_T), pointer :: SPACERADIANCE ! Space radiance
    type (VectorValue_T), pointer :: TEMP          ! Temperature quantity

    integer :: CHANNEL                  ! Loop counter
    integer :: INSTANCE                 ! Loop counter
    integer :: MAFTINSTANCE             ! Temperature instance closest to this MAF
    integer :: MAXNOFREQS               ! Used for sizing arrays
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXSUPERSET              ! Max. value of superset
    integer :: MAXVERT                  ! Number of points in gl grid
    integer :: N2LVL                    ! Twice size of tangent grid
    integer :: NLVL                     ! Size of tangent grid
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOSPECIES                ! Number of molecules we're considering
    integer :: NOUSEDCHANNELS           ! Number of channels to output
    integer :: NO_PHI_T                 ! No. of Temp. profiles in the chunk
    integer :: PHIWINDOW                ! Copy of forward model config%phiWindow
    integer :: SHAPEIND                 ! Index into filter shapes
    integer :: SIDEBANDSTART            ! Loop limit
    integer :: SIDEBANDSTEP             ! Loop step
    integer :: SIDEBANDSTOP             ! Loop limit
    integer :: SIGIND                   ! Loop counter
    integer :: SPECIE                   ! Loop counter
    integer :: STATUS                   ! From allocates etc.
    integer :: SURFACE                  ! Loop counter
    integer :: THISSIDEBAND             ! Loop counter
    integer :: TOTALSIGNALS             ! Used when hunting for pointing grids
    integer :: WHICHPATTERN             ! Index of antenna pattern
    integer :: WHICHPOINTINGGRID        ! Index of pointing grid
    integer :: WINDOWFINISH             ! Range of window
    integer :: WINDOWSTART              ! Range of window

    real (r8) :: CENTERFREQ             ! Of band
    real (r8) :: CENTER_ANGLE           ! For angles
    real (r8) :: R                      ! To convert the kind of output from
    !                                     Freq_Avg
    real (r8) :: THISRATIO              ! A sideband ratio

    integer, dimension(:), pointer :: CHANNELINDEX ! E.g. 1..25
    integer, dimension(:), pointer :: GRIDS ! Frq grid for each tan_press
    integer, dimension(:), pointer :: SUPERSET ! Result of AreSignalsSuperset
    integer, dimension(:), pointer :: USEDCHANNELS ! Array of indices used
    integer, dimension(:), pointer :: USEDSIGNALS ! Array of indices used

    logical :: FOUNDINFIRST                     ! Flag to indicate derivatives

    real(r8), dimension(:,:),   pointer :: D2X_DXDT    ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:,:,:), pointer :: DH_DT_PATH  ! (pathSize, Tsurfs, Tinstance)
    real(r8), dimension(:), allocatable :: Dum
    real(r8), dimension(:,:),   pointer :: DX_DT       ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:),     pointer :: FREQUENCIES ! Frequency points
    real(r8), dimension(:,:),   pointer :: I_STAR_ALL    ! (noMIFs,noChans)
    real(r4), dimension(:,:,:,:,:), pointer :: K_ATMOS ! (channel,Nptg,mxco,mnp,Nsps)
    real(r4), dimension(:,:,:,:), pointer :: K_TEMP    ! (channel,Nptg,mxco,mnp)
    real(r8), dimension(:),     pointer :: PTG_ANGLES  ! (no_tan_hts)
    real(r8), dimension(:,:),   pointer :: RADIANCES     ! (Nptg,noChans)
    real(r8), dimension(:),     pointer :: RadV
    real(r8), dimension(:,:),   pointer :: REF_CORR    ! (n2lvl, no_tan_hts)
    real(kind(k_temp_frq%values)), &
      & dimension(:), pointer :: TOAVG ! Stuff to be passed to frq.avg.
    real(r8), dimension(:),     pointer :: T_SCRIPT    ! (n2lvl)
    real(r8), dimension(:),     pointer :: TAU         ! (n2lvl)

    integer, dimension(1) :: WHICHPOINTINGGRIDASARRAY ! Result of minloc
    integer, dimension(1) :: WHICHPATTERNASARRAY      ! Result of minloc

    type(path_vector), dimension(:), allocatable :: N_PATH    ! (No_tan_hts)

    ! dimensions of SPSFUNC_PATH are: (Nsps,No_tan_hts)
    type(path_vector), allocatable, dimension(:,:) :: SPSFUNC_PATH

    ! dimensions of SPS_PHI_LOOP  are: (Nsps,No_tan_hts)
    ! dimensions of SPS_ZETA_LOOP are: (Nsps,No_tan_hts)
    type(path_int_vector_2d), allocatable, dimension(:,:) :: SPS_PHI_LOOP
    type(path_int_vector_2d), allocatable, dimension(:,:) :: SPS_ZETA_LOOP

    real(r8) :: tol
    Integer :: no_midval_ndx, no_gl_ndx, max_nl, max_num_frq
    Integer, DIMENSION(:,:), ALLOCATABLE :: gl_ndx, midval_ndx
    real(r8), dimension(:,:),  pointer :: midval_delta   ! (N2lvl,Nsps)

    Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

    type(signal_t) :: FirstSignal
    type(signal_t) :: ThisSignal
    type(catalog_T), dimension(:), pointer :: My_Catalog

    ! Executable code --------------------------------------------------------

    if ( toggle(emit) ) call trace_begin ( 'ForwardModel' )

    ! Nullify a bunch of pointers so that Allocate_Test doesn't try to
    ! deallocate them.  We don't want them to be initialized NULL()
    ! because that makes them SAVEd.

    nullify (beta_path, channelIndex, d2x_dxdt, dh_dt_path, dx_dt, &
      & frequencies, grids, i_star_all, k_atmos, k_temp, my_Catalog, &
      & ptg_angles, radiances, radV, ref_corr, superset, t_script, &
      & tau, usedChannels, usedSignals, midval_delta )

    ! Identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    firstSignal = fwdModelConf%signals(1)

    ! Now make sure all the signals we're dealing with are same module,
    ! radiometer and sideband.
    if ( any( fwdModelConf%signals%sideband .ne. &
      & firstSignal%sideband ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed sidebands in forward model config")
    if ( any( fwdModelConf%signals%radiometer .ne. &
      & firstSignal%radiometer ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed radiometers in forward model config")

    ! Now from that we identify the radiance quantity we'll be outputting
    firstRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal=firstSignal%index, sideband=firstSignal%sideband )

    ! Identify the appropriate state vector components, save vmrs for later
    temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature )
    ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=firstSignal%instrumentModule )
    elevOffset => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_elevOffset, radiometer=firstSignal%radiometer )
    orbIncline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_orbitInclination )
    spaceRadiance => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_spaceRadiance )
    earthRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_earthRefl )
    refGPH => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_refGPH )
    losVel => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_losVel, instrumentModule=firstSignal%instrumentModule )
    scGeocAlt => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scGeocAlt )
    sidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_sidebandRatio, signal= firstSignal%index, noError=.true. )

    ! We won't seek for molecules here as we can't have an array of pointers.
    ! When we do want molecule i we would do something like
    ! vmr => GetVectorQuantityBytype (fwdModelIn, fwdModelExtra, &
    !   quantityType=l_vmr, molecule=fwdModelConf.molecules(i))

    ! Now we're going to validate the quantities we've been given, don't forget
    ! we already know what their quantityType's are as that's how we found them
    !, so we don't need to check that.
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )
    if ( .not. ValidateVectorQuantity(elevOffset, verticalCoordinate=(/l_none/), &
      & frequencyCoordinate=(/l_none/), noInstances=(/1/)) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'elevOffset' )
    ! There will be more to come here.

    noSpecies = size(fwdModelConf%molecules)

    !  Create a subset of the catalog composed only of those molecules to be
    !  used for this run

    maxNoFSurfs = 0
    do specie = 1, noSpecies
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=fwdModelConf%molecules(specie) )
      maxNoFSurfs = max(maxNoFSurfs, f%template%noSurfs)
    end do

    allocate ( My_Catalog(noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )

    max_nl = 1 
    do j = 1, noSpecies
      Spectag = spec_tags(fwdModelConf%molecules(j))
      do i = 1, Size(Catalog)
        if ( Catalog(i)%Spec_Tag == Spectag ) then
          My_Catalog(j) = Catalog(i)
          m = Size(Catalog(i)%Lines)
          if(m > max_nl) max_nl = m
          EXIT
        end if
      end do
    end do

    ! Get the max. dimension in zeta coeff. space and phi coeff. space
    ! (To be used later in rad_tran_wd, for automatic arrays asignement)
    max_phi_dim = 1
    max_zeta_dim = 1
    if ( FwdModelConf%temp_der ) then
      max_zeta_dim = temp%template%noSurfs
      max_phi_dim = temp%template%noInstances
    end if

    if ( fwdModelConf%atmos_der ) then
      do k = 1, noSpecies
        if ( fwdModelConf%moleculeDerivatives(k) ) then
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            &     quantityType=l_vmr, molecule=fwdModelConf%molecules(k))
          j = f%template%noInstances
          max_phi_dim = max(max_phi_dim,j)
          j = f%template%noSurfs
          max_zeta_dim = max(max_zeta_dim,j)
        end if
      end do
    end if

    ! Deal with fmStat%rows
    if ( present(Jacobian) .and. ( .not. associated (fmStat%rows) ) ) then
      call Allocate_test ( fmStat%rows, Jacobian%row%nb, 'fmStat%rows', &
        & ModuleName)
      fmStat%rows = .false.
    endif

    ! Get some dimensions that we'll use a lot
    noMAFs = firstRadiance%template%noInstances
    noMIFs = firstRadiance%template%noSurfs
    no_phi_t = temp%template%noInstances
    no_tan_hts = FwdModelConf%TangentGrid%nosurfs
    nlvl=size(FwdModelConf%integrationGrid%surfs)
    n2lvl=2*nlvl
!   maxVert = n2lvl * (NG+1)
    maxVert = 2 * ( (NG+1) * (nlvl-1) + 1)
    phiWindow = FwdModelConf%phiWindow

    tol = +0.2
!   tol = -0.2

    if ( toggle(emit) ) then
      print*,'Dimensions:'
      print*,'noMAFs:',noMAFs
      print*,'no_phi_t:',no_phi_t
      print*,'no_tan_hts:',no_tan_hts
      print*,'maxVert:',maxVert
      print*,'nlvl:',nlvl
      print*,'n2lvl:',n2lvl
      print*,'phiWindow:',phiWindow
      print*,'noSpecies:',noSpecies
      print*,'maxNoFSurfs:',maxNoFSurfs
      print*,'MAF:',fmStat%maf
      print*,'tol:',tol
    end if

    ! Work out which channels are used, also check we have radiances for them.
    noUsedChannels = 0
    do sigInd = 1, size(fwdModelConf%signals)
      thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
        & signal=fwdModelConf%signals(sigInd)%index, sideband=firstSignal%sideband )
      if ( .not. ValidateVectorQuantity(thisRadiance, minorFrame=.true.,&
        & frequencyCoordinate=(/l_channel/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'radiance' )
      noUsedChannels = noUsedChannels + &
        & count( fwdModelConf%signals(sigInd)%channels )
    end do
    call allocate_test ( usedChannels, noUsedChannels, &
      & 'usedChannels', ModuleName )
    call allocate_test ( usedSignals, noUsedChannels, &
      & 'usedSignals', ModuleName )
    channel = 1
    do sigInd = 1, size(fwdModelConf%signals)
      do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
        if (fwdModelConf%signals(sigInd)%channels(i)) then
          usedChannels(channel) = i
          usedSignals(channel) = sigInd
          channel = channel + 1
        end if
      end do
    end do

    ! --------------- Hydrostatic stuff ---------------------------------
    if ( fmStat%newHydros ) then

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.hydrostatic' )

      ! Now we're going to create the many temporary arrays we need
      allocate ( ifm%ndx_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'ndx_path' )
      allocate ( ifm%dhdz_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'dhdz_path' )
      allocate ( ifm%h_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'h_path' )
      allocate ( ifm%phi_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'phi_path' )
      allocate ( ifm%t_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'t_path' )
      allocate ( ifm%z_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'z_path' )

      allocate ( ifm%eta_phi(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'eta_phi' )

      allocate ( ifm%elvar(no_phi_t), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'elvar' )

      call allocate_test ( ifm%geoc_lat, no_phi_t, 'geoc_lat', ModuleName )
      call allocate_test ( ifm%e_rad, no_phi_t, 'e_rad', ModuleName )

      call allocate_test ( ifm%z_glgrid, maxVert/2, 'z_glgrid', ModuleName )
      call allocate_test ( ifm%h_glgrid, maxVert, no_phi_t, 'h_glgrid', &
        &  ModuleName )
      call allocate_test ( ifm%t_glgrid, maxVert, no_phi_t, 't_glgrid', &
        &  ModuleName )
      call allocate_test ( ifm%dh_dt_glgrid, maxVert, no_phi_t, &
        & temp%template%noSurfs,'dh_dt_glgrid', ModuleName )
      call allocate_test ( ifm%dhdz_glgrid, maxVert, no_phi_t, &
        &  'dhdz_glgrid', ModuleName )
      call allocate_test ( ifm%tan_hts, no_tan_hts, no_phi_t, 'tan_hts', &
        &  ModuleName )
      call allocate_test ( ifm%tan_temp, no_tan_hts, no_phi_t, 'tan_hts', &
        &  ModuleName )
      call allocate_test ( ifm%tan_dh_dt, no_tan_hts, no_phi_t, &
        & temp%template%noSurfs, 'tan_dh_dt', ModuleName )
      call Allocate_test( ifm%closestInstances, noMAFs, 'closestInstances', ModuleName)

      ! Setup for hydrostatic calculation
      call FindClosestInstances ( temp, firstRadiance, ifm%closestInstances )

      do i = 1, no_phi_t
        phi_tan = Deg2Rad*temp%template%phi(1,i)
        geod_lat= Deg2Rad*temp%template%geodLat(1,i)
        call geoc_geod_conv ( ifm%elvar(i), orbIncline%values(1,1), &
          &  phi_tan, geod_lat, ifm%geoc_lat(i), ifm%E_rad(i) )
      end do

      ! Now compute a hydrostatic grid given the temperature and refGPH
      ! information.
      call hydrostatic_model ( FwdModelConf%SurfaceTangentIndex, &
        &  no_phi_t, ifm%geoc_lat, 0.001*refGPH%values(1,:), &
        &  refGPH%template%surfs(1,1), &
        &  FwdModelConf%integrationGrid%surfs, &
        &  temp%template%surfs(:,1), temp%values, &
        &  ifm%z_glgrid, ifm%h_glgrid, ifm%t_glgrid, &
        &  ifm%dhdz_glgrid, ifm%dh_dt_glgrid, &
        &  FwdModelConf%TangentGrid%surfs, &
        &  ifm%tan_hts, ifm%tan_temp, ifm%tan_dh_dt, &
        &  ifm%gl_count, Ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Hydrostatic model failed' )

      ! Now compute stuff along the path given this hydrostatic grid.
      call comp_path_entities ( temp, ifm%closestInstances, &
        &  FwdModelConf%integrationGrid%noSurfs, &
        &  temp%template%noSurfs, ifm%gl_count, ifm%ndx_path, ifm%z_glgrid, &
        &  ifm%t_glgrid, ifm%h_glgrid, ifm%dhdz_glgrid, ifm%tan_hts,        &
        &  no_tan_hts, ifm%z_path, ifm%h_path, ifm%t_path, ifm%phi_path,    &
        &  ifm%dhdz_path, ifm%eta_phi, no_phi_t,                            &
        &  temp%template%phi(1,:)*Deg2Rad, noMAFs, phiWindow, ifm%elvar, Ier)
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Hydrostatic model failed' )

      fmStat%newHydros = .false.

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.hydrostatic' )
    end if

    ! ------ End of hydrostatic setup stuff --------------------------

    ! Skip this MAF if it's in an overlap region possibly
    maf=fmStat%maf
    if ( (.not. fwdModelConf%skipOverlaps) .or. &
      & ( maf > firstRadiance%template%noInstancesLowerOverlap .and. &
      &   maf <= noMAFs - firstRadiance%template%noInstancesUpperOverlap ) ) then

      ! ------ Begin main MAF Specific stuff ---------------------------

      ! Now allocate other stuff
      call allocate_test ( t_script, n2lvl, 't_srcipt', ModuleName )
      call allocate_test ( ref_corr, n2lvl, no_tan_hts, 'ref_corr', ModuleName )
      call allocate_test ( tau, n2lvl, 'tau', ModuleName )
      call allocate_test ( ptg_angles, no_tan_hts, 'ptg_angles', ModuleName )

      call allocate_test ( midval_delta, n2lvl, noSpecies, &
                     &    'midval_delta', ModuleName )

      allocate ( k_atmos_frq(noSpecies), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'k_atmos_frq' )

      allocate ( n_path(No_tan_hts), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'n_path' )
      allocate ( spsfunc_path(noSpecies,No_tan_hts), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'spsfunc_path' )
      call allocate_test ( dx_dt, No_tan_hts, temp%template%noSurfs, &
        & 'dx_dt', ModuleName )
      call allocate_test ( d2x_dxdt, No_tan_hts, temp%template%noSurfs, &
        & 'd2x_dxdt', ModuleName )

      allocate ( sps_phi_loop(noSpecies,No_tan_hts), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'sps_phi_loop' )
      allocate ( sps_zeta_loop(noSpecies,No_tan_hts), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'sps_zeta_loop' )

      call allocate_test ( radiances, no_tan_hts, noUsedChannels, &
        & 'Radiances', ModuleName )
      call allocate_test ( i_star_all, noUsedChannels, noMIFs, &
        & 'i_star_all', ModuleName )

      ! Now set radiances to zero, forward model just adds in terms
      do i = 1, noUsedChannels
        thisRadiance = GetVectorQuantityByType (fwdModelOut, &
          & quantityType=l_radiance, &
          & signal=fwdModelConf%signals(usedSignals(i))%index, &
          & sideband=firstSignal%sideband )
        ch = usedChannels(i)
        do j = 1, noMIFs
          thisRadiance%values( ch + (j-1)*thisRadiance%template%noChans, fmStat%maf) = 0.0
        end do
      end do

      if ( fwdModelConf%signals(1)%sideband == 0 ) then
        if (.not. associated (sidebandRatio) ) &
          & call MLSMessage(MLSMSG_Error,ModuleName, &
          & "No sideband ratio supplied")
        sidebandStart = -1
        sidebandStop = 1
        sidebandStep = 2
      else
        sidebandStart = fwdModelConf%signals(1)%sideband
        sidebandStop = sideBandStart
        sidebandStep = 1
      endif

      ! ----------------- Begin loop over sidebands -----------------------
      do thisSideband = sidebandStart, sidebandStop, sidebandStep


        if ( toggle(emit) .and. levels(emit) > 0 ) then
          call trace_begin ( 'ForwardModel.sideband' )
          call output ( ' Doing sideband ' )
          call output ( thisSideband )
          call output ( ' (' ); call output ( sidebandStart )
          call output ( ', ' ); call output ( sidebandStop )
          call output ( ')', advance='yes' )
        end if

        ! Now code splits into two sections, one for when we're doing frequency
        ! averaging, and one when we're not.
        if ( fwdModelConf%do_freq_avg ) then ! --- Doing freq. avg. ---
          if ( toggle(emit) .and. levels(emit) > 0 ) &
            & call trace_begin ( 'ForwardModel.FreqAvg' )

          call allocate_test ( superset, size(pointingGrids), &
            & 'superset', ModuleName )
          do i = 1, size(pointingGrids)
            superset(i) = AreSignalsSuperset ( pointingGrids(i)%signals, &
              & fwdModelConf%signals, sideband=thisSideband )
          end do
          if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching pointing frequency grids." )

          maxSuperset = maxval ( superset )
          where ( superset < 0 )
            superset = maxSuperset + 1
          end where
          whichPointingGridAsArray = minloc ( superset )
          whichPointingGrid = whichPointingGridAsArray(1)
          call deallocate_test ( superset, 'superset', ModuleName )

          if ( toggle(emit) ) then
            call output ( 'Using pointing frequency grid: ' )
            call output ( whichPointingGrid, advance='yes' )
          end if

          ! Now we've identified the pointing grids.  Locate the tangent grid
          ! within it.
          call allocate_test ( grids, FwdModelConf%TangentGrid%nosurfs, &
            "Grids", ModuleName )
          call Hunt ( PointingGrids(whichPointingGrid)%oneGrid%height, &
            & FwdModelConf%TangentGrid%surfs, grids, allowTopValue=.true. )
          if ( toggle(emit) .and. levels(emit) > 0 ) &
            & call trace_end ( 'ForwardModel.FreqAvg' )

        else ! ------------------------- Not frequency averaging ---------

          if ( toggle(emit) .and. levels(emit) > 0 ) &
            & call trace_begin ( 'ForwardModel.NotFreqAvg' )

          call allocate_test ( frequencies,noUsedChannels, "frequencies", ModuleName )
          do channel = 1, noUsedchannels
            frequencies(channel) = &
              & fwdModelConf%signals(usedSignals(channel))%centerFrequency + &
              & fwdModelConf%signals(usedSignals(channel))% &
              &  frequencies(usedChannels(channel))
          end do
          select case ( thisSideband )
          case ( -1 )
            frequencies = firstSignal%lo - frequencies
          case ( +1 )
            frequencies = firstSignal%lo + frequencies
          case ( 0 )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Folded signal requested in forward model' )
          case default
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Bad value of signal%sideband' )
          end select
          noFreqs = noUsedChannels
          if ( toggle(emit) .and. levels(emit) > 0 ) &
            & call trace_end ( 'ForwardModel.NotFreqAvg' )
        end if

        ! ----------- Done the gnarly frequency stuff ----------

        ! Now work out what `window' we're inside.  This will need to be changed
        ! a bit in later versions to avoid the noMAFS==noTemp/f instances
        ! assertion

        mafTInstance = ifm%closestInstances(maf)

        windowStart  = max(1, mafTInstance - phiWindow/2)
        windowFinish = min(mafTInstance + phiWindow/2, no_phi_t)

        if ( toggle(emit) .and. levels(emit) > 0 ) then
          print *, 'Doing MAF: ', maf
          Print *, 'mafTInstance:',mafTInstance
          print *, 'WindowStart:',WindowStart
          print *, 'WindowFinish:',WindowFinish
        end if

        allocate ( k_temp(noUsedChannels, no_tan_hts, temp%template%noSurfs, &
          & windowStart:windowFinish), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate//'k_temp' )
        allocate ( k_atmos(noUsedChannels, no_tan_hts, maxNoFSurfs, &
          & windowStart:windowFinish, noSpecies), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate//'k_atmos' )

        ! Compute the specie function (spsfunc) and the refraction along
        ! all the paths for the current maf

        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.get_path_spsfunc_ngrid' )
        Call get_path_spsfunc_ngrid ( fwdModelIn, fwdModelExtra, &
          &  fwdModelConf%molecules, ifm%ndx_path(:,maf), no_tan_hts, &
          &  ifm%z_path(:,maf), ifm%t_path(:,maf), ifm%phi_path(:,maf), &
          &  n_path, spsfunc_path, sps_zeta_loop, sps_phi_loop, ier)
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.get_path_spsfunc_ngrid' )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'get_path_sps_fun_ngrid failed' )

        !??? Choose better value for phi_tan later
        phi_tan = Deg2Rad * temp%template%phi(1,mafTInstance)

        ! Compute the ptg_angles (chi) for Antenna convolution, also the
        ! derivatives of chi w.r.t to T and other parameters
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.get_chi_angles' )
        call get_chi_angles ( ifm%ndx_path(:,maf), n_path, &
          &  fwdModelConf%tangentGrid%surfs, &
          &  ifm%tan_hts(:,mafTInstance),ifm%tan_temp(:,mafTInstance),&
          &  phi_tan,ifm%elvar(maf)%Roc,&
          &  0.001*scGeocAlt%values(1,1),  &
          &  elevOffset%values(1,1), &
          &  ifm%tan_dh_dt(:,mafTInstance,:), no_tan_hts, &
          &  temp%template%noSurfs, &
          &  temp%template%surfs(:,1), &
          &  fwdModelConf%SurfaceTangentIndex, &
          &  center_angle, ptg_angles, dx_dt, d2x_dxdt, ier )
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.get_chi_angles' )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'get_chi_angles failed' )

        ! Compute the refraction correction scaling matrix for this mmaf:
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.refraction_correction' )
        call refraction_correction ( no_tan_hts, ifm%tan_hts(:,mafTInstance), &
          &  ifm%h_path(:,maf), n_path, ifm%ndx_path(:,maf),      &
          &  ifm%E_rad(mafTInstance), ref_corr )
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.refraction_correction' )

        Radiances = 0.0

        ! If we're not doing frequency averaging, instead outputting radiances
        ! corresponding to delta function responses, we can set up the frequency
        ! information here.  In the more common case where we are doing the
        ! averaging, the frequency grid varies from pointing to pointing, and is
        ! allocated inside the pointing loop.

        ! First we have a mini loop over pointings to work out an upper limit
        ! for the number of frequencies we're going to be dealing with
        if ( fwdModelConf%do_freq_avg ) then
          maxNoFreqs = size(PointingGrids(whichPointingGrid)%OneGrid(grids(1))%Frequencies)
          do ptg_i = 2, no_tan_hts - 1
            maxNoFreqs = max ( maxNoFreqs, size(PointingGrids(whichPointingGrid) &
              & %OneGrid(grids(ptg_i))%Frequencies) )
          end do
        else
          maxNoFreqs = noFreqs
        end if

        ! Now allocate arrays this size
        if ( fwdModelConf%temp_der ) then
          allocate ( k_temp_frq%values( maxNoFreqs, temp%template%noSurfs, &
            & windowStart:windowFinish), stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//'k_temp_frq' )
          k_temp_frq%values = 0.0_r8
        end if

        call allocate_test ( radV,maxNoFreqs, 'radV', ModuleName )

        do specie = 1, noSpecies
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_vmr, molecule=fwdModelConf%molecules(specie) )

          ! Allocate intermediate space for vmr derivatives
          if ( fwdModelConf%moleculeDerivatives(specie) ) then
            allocate ( k_atmos_frq(specie)%values(maxNoFreqs,f%template%noSurfs,&
              & windowStart:windowFinish), stat=status )
            if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
              & MLSMSG_Allocate//'k_atmos_frq' )
          end if

        end do ! End loop over species

!  Find the maximum number of frequencies over all heights

        max_num_frq = 1
        i = whichPointingGrid
        if ( FwdModelConf%do_freq_avg ) then
          do j = 1, no_tan_hts
            m = Size(PointingGrids(i)%oneGrid(grids(j))%frequencies)
            if ( m > max_num_frq) max_num_frq = m
          end do
        end if
!
        ! Now we can go ahead and loop over pointings
        ! ------------------------------ Begin loop over pointings --------
        do ptg_i = 1, no_tan_hts - 1

          if ( toggle(emit) .and. levels(emit) > 1 ) then
            call trace_begin ( 'ForwardModel.Pointing' )
            call output ( 'Ptg = ' ); call output ( ptg_i, advance='yes' )
          end if

          k = ptg_i
          h_tan = ifm%tan_hts(k,mafTInstance)
          lmax = ubound(ifm%eta_phi(ptg_i,maf)%values,2)
!
          ifm%elvar(maf)%ht = h_tan
          ifm%elvar(maf)%Rr = ifm%elvar(maf)%ht + ifm%elvar(maf)%RoC
          ifm%elvar(maf)%ht2 = ifm%elvar(maf)%Rr * ifm%elvar(maf)%Rr

 ! Compute the beta's along the path, for this tanget hight and this mmaf:

          brkpt = ifm%ndx_path(ptg_i,maf)%break_point_index
          no_ele = ifm%ndx_path(ptg_i,maf)%total_number_of_elements
!
          if(ptg_i == 1) then
            j = no_ele / Ng
            allocate(midval_ndx(j,2),gl_ndx(j,2),stat=status)
            if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
              & MLSMSG_Allocate//'midval_ndx/gl_ndx' )
          endif

          ! If we're doing frequency averaging, get the frequencies we need for
          ! this pointing.
          if ( FwdModelConf%do_freq_avg ) then
            frequencies => PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
            noFreqs = size(frequencies)
          end if ! If not, we dealt with this outside the loop

          if(noFreqs > 1) tau(1:) = 1.0
!
! Compute ALL the slabs_prep entities over the path's GL grid for this
! pointing & mmaf:
!
          Call get_gl_slabs_arrays(my_Catalog,ifm%z_path(ptg_i,maf),&
             & ifm%t_path(ptg_i,maf),0.001*losVel%values(1,maf), &
             & gl_slabs,ptg_i,no_ele,max_nl,ier)
          if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'get_gl_slabs_arrays' )
!
! Compute the beta's along the path on the COARSE grid only, for this pointing
! and this mmaf (recall that k=ptg_i):
!
          Call get_coarse_beta_path(k,my_Catalog,ifm%ndx_path(k,maf), &
             & gl_slabs,noFreqs,max_num_frq,h_tan,frequencies,        &
             & ifm%z_path(k,maf),ifm%h_path(k,maf),ifm%t_path(k,maf), &
             & fwdModelConf%frqGap,fwdModelConf%temp_der,             &
             & fwdModelConf%spect_der, beta_path, Ier)
          if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'get_coarse_beta_path' )
!
          ! Define the dh_dt_path for this pointing and this MAF:

          ! Need to allocate this even if no derivatives as we pass it

          call allocate_test ( dh_dt_path, no_ele, phiWindow, &
            & temp%template%noSurfs, "dh_dt_path", ModuleName )

          if ( fwdModelConf%temp_der ) then
            allocate ( dum(no_ele), stat=ier )
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
              & MLSMSG_Allocate // 'dum' )
            do j = 1, temp%template%noSurfs
              do i = 1, phiWindow
                m = min(lmax,i+windowStart-1)
                call Lintrp ( ifm%z_glgrid, ifm%z_path(ptg_i,maf)%values,&
                  &           ifm%dh_dt_glgrid(:,m,j), dum, ifm%gl_count,&
                  &           no_ele )
                dh_dt_path(:,i,j) = dum(:) * ifm%eta_phi(ptg_i,maf)%values(:,m)
              end do
            end do
            deallocate ( dum, stat=ier )
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
              & MLSMSG_DeAllocate // 'dum' )
          end if

          ! ------------------------------- Begin loop over frequencies ------
          do frq_i = 1, noFreqs
!
            Call path_contrib(tau, brkpt, no_ele, tol, mid, midval_ndx, &
              &  no_midval_ndx, gl_ndx, no_gl_ndx, Ier)
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Path_contrib failed' )
!
            Call coarse_delta(mid,brkpt,no_ele,ifm%h_path(k,maf), &
              &  beta_path(:,frq_i),spsfunc_path(:,k),        &
              &  noSpecies,Nlvl,ref_corr(:,k),ifm%elvar(maf),     &
              &  midval_delta,Ier)
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'coarse_delta failed' )
!
            Frq = frequencies(frq_i)
            if ( toggle(emit) .and. levels(emit) > 2 ) then
              call trace_begin ( 'ForwardModel.Frequencies' )
              call output ( 'Frq = ' ); call output ( frq_i, advance='yes' )
            end if
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
            Call get_glbeta_path(k,frq_i,my_Catalog,gl_ndx,no_gl_ndx, &
               & gl_slabs,Frq,ifm%z_path(k,maf),ifm%t_path(k,maf),    &
               & fwdModelConf%frqGap,fwdModelConf%temp_der,           &
               & fwdModelConf%spect_der, beta_path, Ier)
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'get_glbeta_path failed' )
!
            Call Rad_Tran ( ifm%elvar(maf), Frq,                         &
              &  fwdModelConf%integrationGrid%noSurfs,h_tan,noSpecies,   &
              &  ifm%z_path(k,maf),ifm%h_path(k,maf),&
              &  ifm%t_path(k,maf),ifm%phi_path(k,maf),                  &
              &  ifm%dHdz_path(k,maf),earthRefl%values(1,1),             &
              &  beta_path(:,frq_i),spsfunc_path(:,k),ref_corr(:,k),     &
              &  spaceRadiance%values(1,1),brkpt,no_ele,mid,ilo,ihi,     &
              &  t_script,tau,midval_ndx,no_midval_ndx,gl_ndx,no_gl_ndx, &
              &  midval_delta,Sps_zeta_loop(:,k),Sps_phi_loop(:,k),Rad,Ier)
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,  &
              & 'rad_tran failed' )

            RadV(frq_i) = Rad

            ! Now, Compute the radiance derivatives:

            !??? Do we need to do this if there's no Jacobian or no derivatives requested ???

            Call Rad_Tran_WD ( FwdModelConf, FwdModelExtra, FwdModelIn,    &
              &  ifm%elvar(maf),frq_i,Frq,noSpecies,ifm%z_path(k,maf),     &
              &  ifm%h_path(k,maf),ifm%t_path(k,maf),ifm%phi_path(k,maf),  &
              &  ifm%dHdz_path(k,maf),beta_path(:,frq_i),spsfunc_path(:,k),&
              &  temp%template%surfs(:,1),temp%template%noSurfs,           &
              &  ref_corr(:,k),temp%template%noInstances,                  &
              &  temp%template%phi(1,:)*Deg2Rad,dh_dt_path,k_temp_frq,     &
              &  k_atmos_frq,brkpt,no_ele,mid,ilo,ihi,t_script,tau,        &
              &  max_zeta_dim,max_phi_dim,midval_ndx,no_midval_ndx,        &
              &  gl_ndx,no_gl_ndx,midval_delta,Sps_zeta_loop(:,k),         &
              &  Sps_phi_loop(:,k),Ier)
            if( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,     &
              & 'rad_tran_wd failed' )

            if ( toggle(emit) .and. levels(emit) > 2 ) &
              & call trace_end ( 'ForwardModel.Frequencies' )

          end do                          ! Frequency loop

          ! ----------------------------- End loop over frequencies ----

          ! Here we either frequency average to get the unconvolved radiances, or
          ! we just store what we have as we're using delta funciton channels

          if ( toggle(emit) .and. levels(emit) > 1 ) &
            & call trace_begin ( 'ForwardModel.FrequencyAvg' )
          if ( fwdModelConf%do_freq_avg ) then
            do i = 1, noUsedChannels
              sigInd = usedSignals(i)
              ch = usedChannels(i)
              if ( toggle(emit) .and. levels(emit) > 2 ) then
                call output ( 'Channel = ' )
                call output ( i )
                call output ( ' ( ' )
                call output ( sigInd )
                call output ( ':' )
                call output ( ch )
                call output ( ' )', advance='yes' )
              end if
              centerFreq = firstSignal%lo + &
                & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), sideband = thisSideband )
              if ( shapeInd == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & "No matching channel shape information" )
              if ( toggle(emit) .and. levels(emit) > 2 ) then
                call output ( 'Using filter shape:' )
                call output ( shapeInd, advance='yes' )
              endif

              call Freq_Avg ( frequencies, &
                & centerFreq+thisSideband * &
                & FilterShapes(shapeInd)%FilterGrid(ch,:), &
                & FilterShapes(shapeInd)%FilterShape(ch,:), RadV, noFreqs,  &
                & Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), Radiances(ptg_i,i) )
            end do
          else
            Radiances(ptg_i,:) = RadV(1:noFreqs)
          end if

          ! Frequency Average the temperature derivatives with the appropriate
          ! filter shapes
          !??? Do we need to do this if there's no Jacobian ???
          if ( fwdModelConf%temp_der ) then
            if ( fwdModelConf%do_freq_avg ) then
              do i = 1, noUsedChannels
                sigInd = usedSignals(i)
                ch = usedChannels(i)
                centerFreq = firstSignal%lo + &
                  & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
                shapeInd = MatchSignal ( filterShapes%signal, &
                  & fwdModelConf%signals(sigInd), &
                  & sideband = thisSideband, channel=ch )
                do instance = lbound(k_temp_frq%values,3), &
                  & ubound(k_temp_frq%values,3)
                  do surface = 1, temp%template%noSurfs
                    ToAvg => k_temp_frq%values(1:noFreqs,surface,instance)
                    call Freq_Avg ( frequencies,                        &
                      & centerFreq+thisSideband* &
                      & FilterShapes(shapeInd)%FilterGrid(ch,:), &
                      & FilterShapes(shapeInd)%FilterShape(ch,:),&
                      & real(ToAvg,r8), noFreqs, &
                      & Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), r )
                    k_temp(i,ptg_i,surface,instance) = r
                  end do                  ! Surface loop
                end do                    ! Instance loop
              end do                      ! Channel loop
            else
              do i = 1, noUsedChannels
                k_temp(i,ptg_i,1:temp%template%noSurfs,:) = &
                  &  k_temp_frq%values(i,1:temp%template%noSurfs,:)
              end do
            end if
          end if

          ! Frequency Average the atmospheric derivatives with the appropriate
          ! filter shapes
          !??? Do we need to do this if there's no Jacobian ???
          do specie = 1, noSpecies
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
                &  quantityType=l_vmr, molecule=fwdModelConf%molecules(specie))
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  ch = usedChannels(i)
                  centerFreq = firstSignal%lo + &
                    & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=ch )
                  do instance = lbound(k_atmos_frq(specie)%values,3),&
                    & ubound(k_atmos_frq(specie)%values,3)
                    do surface = 1, f%template%noSurfs
                      ToAvg => k_atmos_frq(specie)%values(1:noFreqs,surface,instance)
                      call Freq_Avg ( frequencies,                      &
                        & centerFreq+thisSideband * &
                        & FilterShapes(shapeInd)%FilterGrid(ch,:), &
                        & FilterShapes(shapeInd)%FilterShape(ch,:), &
                        & real(ToAvg,r8),  &
                        & noFreqs, Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), r )
                      k_atmos(i,ptg_i,surface,instance,specie) = r
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! Else not frequency averaging
                surface = f%template%noSurfs
                do i = 1, noUsedChannels
                  k_atmos(i,ptg_i,1:surface,:,specie) = &
                    &  k_atmos_frq(specie)%values(i,1:surface,:)
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over species

          if ( toggle(emit) .and. levels(emit) > 1 ) &
            & call trace_end ( 'ForwardModel.FrequencyAvg' )

          call deallocate_test ( dh_dt_path, 'dh_dt_path', ModuleName )

          if ( toggle(emit) .and. levels(emit) > 1 ) &
            & call trace_end ( 'ForwardModel.Pointing' )

        end do                            ! Pointing Loop
        ! ---------------------------------- End of Pointing Loop ---------------

        ! Complete the radiance's last location; also complete k_temp last
        ! location.

        do i = 1, noUsedChannels
          ch = usedChannels(i)
          Radiances(no_tan_hts,i) = Radiances(no_tan_hts-1,i)
          if ( FwdModelConf%temp_der ) then
            n = temp%template%noSurfs
            k_temp(i,no_tan_hts,1:n,:) = k_temp(i,no_tan_hts-1,1:n,:)
          end if
          if ( FwdModelConf%atmos_der ) then
            do m = 1, noSpecies
              f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
                & quantityType=l_vmr, molecule=fwdModelConf%molecules(m),&
                & foundInFirst=foundInFirst )
              if ( foundInFirst ) then
                n = f%template%noSurfs
                k_atmos(i,no_tan_hts,1:n,:,m)= k_atmos(i,no_tan_hts-1,1:n,:,m)
              end if
            end do
          end if
        end do

        !  Here comes the Convolution code
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.Convolution' )
        call allocate_test ( superset, size(antennaPatterns), &
          & 'superset', ModuleName )
        do i = 1, noUsedChannels

          ch = usedChannels(i)
          sigInd = usedSignals(i)
          thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
            & signal=fwdModelConf%signals(sigInd)%index, &
            & sideband=firstSignal%sideband )

          if ( sidebandStart /= sidebandStop ) then ! We're folding
            thisRatio = sidebandRatio%values(ch,1)
            if ( thisSideband == 1 ) thisRatio = 1.0 - thisRatio
          else ! Otherwise, want just unfolded signal
            thisRatio = 1.0
          end if

          if ( FwdModelConf%do_conv ) then

            do j = 1, size(antennaPatterns)
              superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
                & fwdModelConf%signals, sideband=thisSideband, channel=ch )
            end do
            if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "No matching antenna patterns." )

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

            call convolve_all ( fwdModelConf, fwdModelIn, maf, ch, &
              &     windowStart, windowFinish, mafTInstance-windowStart+1, &
              &     temp, ptan, thisRadiance, &
              &     FwdModelConf%tangentGrid%surfs, ptg_angles, &
              &     ifm%tan_temp(:,maf), dx_dt, d2x_dxdt, si, center_angle, &
              &     Radiances(:, i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
              &     thisRatio, Jacobian, fmStat%rows,  &
              &     antennaPatterns(whichPattern), ier )
            !??? Need to choose some index other than 1 for AntennaPatterns ???
            if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'convolve_all failed' )
          else
            call no_conv_at_all ( fwdModelConf, fwdModelIn, maf, ch,  &
              &     windowStart, windowFinish, temp, ptan, thisRadiance, &
              &     FwdModelConf%tangentGrid%surfs, &
              &     Radiances(:,i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
              &     thisRatio, Jacobian, fmStat%rows )

          end if

        end do                            ! Channel loop
        call deallocate_test ( superset, 'superset', ModuleName )
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.Convolution' )

        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.sideband' )
      end do
      ! ---------------------------- End of loop over sideband ------------------

      ! ------------------------------ End of Major Frame Specific stuff --------

      if ( associated(beta_path) ) then
        do i = 1, size(beta_path,1)
          do j = 1, size(beta_path,2)
            deallocate ( beta_path(i,j)%values, beta_path(i,j)%t_power, &
              & beta_path(i,j)%dbeta_dw, beta_path(i,j)%dbeta_dn, &
              & beta_path(i,j)%dbeta_dnu, STAT=k )
          end do
        end do
      end if

      deallocate ( beta_path, STAT=k )

      if ( FwdModelConf%temp_der ) call deallocate_test &
        & ( k_temp_frq%values, "k_temp_frq%values", ModuleName )
      if ( FwdModelConf%atmos_der ) then
        do j = 1, noSpecies
          call deallocate_test(k_atmos_frq(j)%values,"k_atmos_frq(j)%values",&
            & ModuleName )
        end do
      end if

      if ( index(switches,'rad') /= 0 ) then
        ! *** DEBUG Print
        if ( FwdModelConf%do_conv ) then
          print *,'Convolution: ON'
        else
          print *,'Convolution: OFF'
        end if

        if ( FwdModelConf%do_freq_avg ) then
          print *,'Frequency Averaging: ON'
        else
          Frq = Frequencies(1)
          print *,'Frequency Averaging: OFF'
          print '(A,f12.4,a)', ' (All computations done at Frq =',Frq,')'
        end if

        print *
        k=ptan%template%noSurfs
        print 901, k
901     format ( 'ptan\',i3.3)
        Print 902,Ptan%values(1:k,maf)
902     format ( 4(4x, f10.7))

        print *
        do i = 1, noUsedChannels
          ch = usedChannels(i)
          print 903, ch, char(92), ptan%template%noSurfs
903       format ( 'ch', i2.2, '_pfa_rad', a1, i3.3 )
          print 905, ( firstRadiance%values(ch+(k-1)*firstRadiance%template%noChans,maf),&
            & k = 1, ptan%template%noSurfs )
905       format ( 4(2x, 1pg15.8) )
        end do
      end if

      do j = 1, No_tan_hts
        deallocate ( n_path(j)%values, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate//'n_path%values' )
        do i = 1, noSpecies
          deallocate ( spsfunc_path(i,j)%values, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Deallocate//'spsfunc_path%values' )
          deallocate ( sps_phi_loop(i,j)%values, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Deallocate//'sps_phi_loop%values' )
          deallocate ( sps_zeta_loop(i,j)%values, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Deallocate//'sps_zeta_loop%values' )
        end do
      end do

      deallocate ( n_path, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'n_path' )

      deallocate ( spsfunc_path, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'spsfunc_path' )

      deallocate ( sps_phi_loop, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'sps_phi_loop' )
      deallocate ( sps_zeta_loop, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'sps_zeta_loop' )

      call deallocate_test ( dx_dt, 'dx_dt', ModuleName )
      call deallocate_test ( d2x_dxdt, 'd2x_dxdt', ModuleName )

      call deallocate_test ( t_script, 't_srcipt', ModuleName )
      call deallocate_test ( ref_corr, 'ref_corr', ModuleName )
      call deallocate_test ( tau, 'tau', ModuleName )
      call deallocate_test ( ptg_angles, 'ptg_angles', ModuleName )

      deallocate ( k_atmos_frq, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'k_atmos_frq' )
      deallocate ( k_temp, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'k_temp' )
      deallocate ( k_atmos, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'k_atmos' )
      deallocate ( My_Catalog, stat=status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'My_catalog' )

      call deallocate_test ( radiances, 'Radiances', ModuleName )
      call deallocate_test ( i_star_all, 'i_star_all', ModuleName )
      call deallocate_test ( radv, 'rad_v', ModuleName )
      call deallocate_test ( grids, 'grids',  ModuleName )

    else
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call output ( 'This MAF skipped as in overlap region', advance='yes' )
    endif ! --------------------------- Possible skip of this major frame

    if ( maf == noMAFs ) fmStat%finished = .true.

    call Deallocate_test ( usedChannels, 'usedChannels', ModuleName )
    call Deallocate_test ( usedSignals, 'usedSignals', ModuleName )

    if ( .not. fwdModelConf%do_freq_avg ) call deallocate_test ( &
      & frequencies, "frequencies", ModuleName )

    if ( toggle(emit) ) call trace_end ( 'ForwardModel' )

  end subroutine FullForwardModel

end module FullForwardModel_m

! $Log$
! Revision 1.1  2001/05/29 22:53:51  livesey
! First version, taken from old ForwardModelInterface.f90
!
