! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullForwardModel_m

  use GLNP, only: NG, GX
  use MLSCommon, only: I4, R4, R8, RP, IP, FINDFIRST
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATEONESLABS, &
                              &  DESTROYCOMPLETESLABS
  use COMP_ETA_DOCALC_NO_FRQ_M, only: comp_eta_docalc_no_frq
  use COMP_SPS_PATH_FRQ_M, only: comp_sps_path_frq
  use EVAL_SPECT_PATH_M, only: EVAL_SPECT_PATH
  use REFRACTION_M, only: REFRACTIVE_INDEX, COMP_REFCOR, PATH_DS_DH
  use TWO_D_HYDROSTATIC_M, only: two_d_hydrostatic
  use METRICS_M
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use RAD_TRAN_M, only: PATH_CONTRIB, RAD_TRAN, DRAD_TRAN_DF,DRAD_TRAN_DT, &
                       &  DRAD_TRAN_DX
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
  use Load_sps_data_m, only: LOAD_SPS_DATA, Grids_T
  use MatrixModule_1, only: MATRIX_T
  use Trace_M, only: Trace_begin, Trace_end
  use MLSSignals_m, only: SIGNAL_T, MATCHSIGNAL, ARESIGNALSSUPERSET
  use String_table, only: GET_STRING, DISPLAY_STRING
  use SpectroscopyCatalog_m, only: CATALOG_T, LINE_T, LINES, CATALOG
  use intrinsic, only: L_TEMPERATURE, L_RADIANCE, L_PTAN, L_ELEVOFFSET, &
    & L_ORBITINCLINATION, L_SPACERADIANCE, L_EARTHREFL, L_LOSVEL,       &
    & L_SCGEOCALT, L_SIDEBANDRATIO, L_NONE, L_CHANNEL, L_VMR, L_REFGPH, LIT_INDICES
  use Units, only: Deg2Rad
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error, MLSMSG_Warning
  use MLSNumerics, only: HUNT
  use Toggles, only: Emit, Gen, Levels, Switches, Toggle
  use Molecules, only: spec_tags
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Output_m, only: OUTPUT
  use ManipulateVectorQuantities, only: FindOneClosestInstance
  use Trace_M, only: Trace_begin, Trace_end

  use Dump_0, only: DUMP

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
  Subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                             &  FwdModelOut, oldIfm, FmStat, Jacobian )
    ! This is the full radiative transfer forward model, the workhorse
    ! code
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Define local parameters
    character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "
    integer, parameter :: NGP1=NG+1     ! NG + 1

    ! Now define local variables, group by type and then
    ! alphabetically

    integer :: BRKPT                    ! Index of midpoint of path
    integer :: CATINDEX                 ! Index for molecule
    integer :: CHANNEL                  ! A Loop counter
    integer :: FRQ_I                    ! Frequency loop index
    integer :: F_LEN                    ! Total number of f's
    integer :: P_LEN                    ! Partial number of f's (No freq.)
    integer :: H2O_IND                  ! Index of h2o inside f array
    integer :: I                        ! Loop index and other uses .
    integer :: IER                      ! Status flag from allocates
    integer :: I_STOP                   ! Upper index for radiance comp.
    integer :: INSTANCE                 ! Loop counter
    integer :: J                        ! Loop index and other uses ..
    integer :: K                        ! Loop index and other uses ..
    integer :: MAF                      ! MAF under consideration
    integer :: MAFTINSTANCE             ! Temperature instance closest to this MAF
    integer :: MAXNOFREQS               ! Used for sizing arrays
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXSUPERSET              ! Max. value of superset
    integer :: MAXVERT                  ! Number of points in gl grid
    integer :: NL                       ! Number of lines
    integer :: NLM1                     ! NLVL - 1
    integer :: NLVL                     ! Size of integration grid
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOSPECIES                ! Number of molecules under consideration
    integer :: NOUSEDCHANNELS           ! How many channels are we considering
    integer :: NO_ELE                   ! Length of a gl path
    integer :: NO_GL_NDX                ! Number of GL points to do
    integer :: NO_TAN_HTS               ! Number of tangent heights
    integer :: NPC                      ! Length of coarse path
    integer :: N_T_PHI                  ! Number of phis for temperature
    integer :: N_T_ZETA                 ! Number of zetas for temperature
    integer :: PHIWINDOW                ! From fwdMdlConfig
    integer :: PTG_I                    ! Loop counter for the pointings
    integer :: SHAPEIND                 ! Index into filter shapes
    integer :: SIDEBANDSTART            ! Loop limit
    integer :: SIDEBANDSTEP             ! Loop step
    integer :: SIDEBANDSTOP             ! Loop limit
    integer :: SIGIND                   ! Signal index, loop counter
    integer :: SPECTAG                  ! A single spectag
    integer :: SURFACETANGENTINDEX      ! Index in tangent grid of earth's surface
    integer :: THISSIDEBAND             ! Loop counter for sidebands
    integer :: WHICHPOINTINGGRID        ! Index into the pointing grids
    integer :: WINDOWFINISH             ! End of temperature `window'
    integer :: WINDOWSTART              ! Start of temperature `window'
    integer :: SPECIE                   ! Loop counter
    integer :: SV_I                     ! Loop index and other uses .
    integer :: SV_DW_LEN                ! Length of DW in vector
    integer :: SV_DN_LEN                ! Length of DN in vector
    integer :: SV_DV_LEN                ! Length of DV in vector
    integer :: SV_START                 ! Temporary sv_i
    integer :: SV_T_LEN                 ! Number of t_phi*t_zeta in the window
    integer :: SURFACE                  ! Loop counter
    integer :: WHICHPATTERN             ! Index of antenna pattern

    logical :: DOTHIS                   ! Flag for lines
    character (len=32) :: molName       ! Name of a molecule

    integer, dimension(1) :: WHICHPOINTINGGRIDASARRAY ! Result of minloc
    integer, dimension(1) :: WHICHPATTERNASARRAY      ! Result of minloc

    integer, dimension(:), pointer :: GRIDS ! Heights in ptgGrid for each tangent
    integer, dimension(:), pointer :: USEDCHANNELS ! Which channel is this
    integer, dimension(:), pointer :: USEDSIGNALS ! Which signal is this channel from
    integer, dimension(:), pointer :: SUPERSET ! Used for matching signals
    integer, dimension(:), pointer :: EXT_IND_C ! Indecies on coarse grid
    integer, dimension(:), pointer :: TAN_INDS ! Index of tangent grid into gl grid
    integer, dimension(:), pointer :: GL_INDS ! Index of GL indecies

    integer, dimension(:,:), pointer :: GL_NDX ! Packed Index array of GL intervals
    integer, dimension(:,:), pointer :: GL_INDGEN ! Temp. array of indecies

    logical, dimension(:), pointer :: DO_GL ! GL indicator
    logical, dimension(:), pointer :: LIN_LOG ! Is this vmr on a log basis? (noSpecies)
    logical, dimension(:), pointer :: LINEFLAG ! Use this line (noLines per species)

    logical, dimension(:,:), pointer :: DO_CALC_ZP  ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_FZP ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_DN ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DV ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DW ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_HYD ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_T ! 'Avoid zeros'

    real(r8) :: FRQ                     ! Frequency
    real(r8) :: CENTERFREQ              ! Of band

    real(rp) :: CENTER_ANGLE            ! For angles
    real(rp) :: DEL_TEMP   ! Temp. step-size in evaluation of Temp. power dep.
    real(rp) :: ELEV_OFFSET             ! Elevation offset
    real(rp) :: E_RFLTY                 ! Earth reflectivity at given tan. point
    real(rp) :: NEG_TAN_HT              ! GP Height (in KM.) of tan. press. below surface
    real(rp) :: PHI_TAN                 ! Phi at tanget point for given pointing
    real(rp) :: R                       ! real variable for various uses
    real(rp) :: RAD                     ! Radiance
    real(rp) :: REQ                     ! Equivalent Earth Radius
    real(rp) :: THISRATIO               ! A sideband ratio

    real(r8), dimension(:), pointer :: FREQUENCIES ! What frequencies to compute for
    real(rp), dimension(:), pointer :: ALPHA_PATH_C ! coarse grid Sing.
    real(rp), dimension(:), pointer :: DEL_S ! Integration lengths along the path
    real(rp), dimension(:), pointer :: DHDZ_PATH ! dH/dZ on path
    real(rp), dimension(:), pointer :: DRAD_DF ! dI/dVmr
    real(rp), dimension(:), pointer :: DRAD_DN ! dI/dN
    real(rp), dimension(:), pointer :: DRAD_DT ! dI/dT
    real(rp), dimension(:), pointer :: DRAD_DV ! dI/dV
    real(rp), dimension(:), pointer :: DRAD_DW ! dI/dW
    real(rp), dimension(:), pointer :: H_PATH ! Heights on path
    real(rp), dimension(:), pointer :: INCOPTDEPTH ! Incremental Optical depth
    real(rp), dimension(:), pointer :: N_PATH ! Refractivity on path
    real(rp), dimension(:), pointer :: ONE_TAN_HT ! ***
    real(rp), dimension(:), pointer :: ONE_TAN_TEMP ! ***
    real(rp), dimension(:), pointer :: PATH_DSDH ! dS/dH on path
    real(rp), dimension(:), pointer :: PHI_BASIS ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DN ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DV ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DW ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_PATH ! Phi's on path
    real(rp), dimension(:), pointer :: P_GLGRID ! Pressure on glGrid surfs
    real(rp), dimension(:), pointer :: P_PATH ! Pressure on path
    real(rp), dimension(:), pointer :: RADV ! Radiances for 1 pointing on freqGrid
    real(rp), dimension(:), pointer :: REF_CORR ! Refraction correction
    real(rp), dimension(:), pointer :: SPS_VALUES ! all vmrs (f_len)
    real(rp), dimension(:), pointer :: TAU ! Optical depth
    real(rp), dimension(:), pointer :: T_PATH ! Temperatures on path
    real(rp), dimension(:), pointer :: T_SCRIPT ! ********
    real(rp), dimension(:), pointer :: XM ! Midpoint of integration grid intervals
    real(rp), dimension(:), pointer :: YM ! Half length of integration grid intervals
    real(rp), dimension(:), pointer :: ZGX ! Gauss weights (with -1 at front)
    real(rp), dimension(:), pointer :: Z_BASIS ! zeta basis per species (n_f_zeta)
    real(rp), dimension(:), pointer :: Z_BASIS_DN ! zeta basis for dw (n_dn_z)
    real(rp), dimension(:), pointer :: Z_BASIS_DV ! zeta basis for dw (n_dv_z)
    real(rp), dimension(:), pointer :: Z_BASIS_DW ! zeta basis for dw (n_dw_z)
    real(rp), dimension(:), pointer :: Z_GLGRID ! Zeta on glGrid surfs
    real(rp), dimension(:), pointer :: Z_PATH ! Zeta on path

    real(rp), dimension(:,:), pointer :: BETA_PATH ! Beta on path
    real(rp), dimension(:,:), pointer :: BETA_PATH_C ! Beta on path coarse
    real(rp), dimension(:,:), pointer :: BETA_PATH_F ! Beta on path fine
    real(rp), dimension(:,:), pointer :: DBETA_DN_PATH_C ! dBeta_dn on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DN_PATH_F ! dBeta_dn on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DT_PATH_C ! dBeta_dt on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DT_PATH_F ! dBeta_dt on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DV_PATH_C ! dBeta_dv on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DV_PATH_F ! dBeta_dv on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DW_PATH_C ! dBeta_dw on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DW_PATH_F ! dBeta_dw on fine grid
    real(rp), dimension(:,:), pointer :: DHDZ_GLGRID ! dH/dZ on glGrid surfs
    real(rp), dimension(:,:), pointer :: DH_DT_PATH ! dH/dT on path
    real(rp), dimension(:,:), pointer :: DX_DT       ! (No_tan_hts, Tsurfs)
    real(rp), dimension(:,:), pointer :: D2X_DXDT    ! (No_tan_hts, Tsurfs)
    real(rp), dimension(:,:), pointer :: ETA_ZP  ! Eta_z x Eta_p
    real(rp), dimension(:,:), pointer :: ETA_FZP ! Eta_z x Eta_p * Eta_f
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DN ! Eta_z x Eta_p for N
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DV ! Eta_z x Eta_p for V
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DW ! Eta_z x Eta_p for W
    real(rp), dimension(:,:), pointer :: ETA_ZXP_T ! Eta_t_z x Eta_t_p
    real(rp), dimension(:,:), pointer :: H_GLGRID ! H on glGrid surfs
    real(rp), dimension(:,:), pointer :: K_ATMOS_FRQ ! Storage for Atmos deriv.
    real(rp), dimension(:,:), pointer :: K_SPECT_DN_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_SPECT_DV_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_SPECT_DW_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_TEMP_FRQ ! Storage for Temp. deriv.
    real(rp), dimension(:,:), pointer :: PTG_ANGLES ! (no_tan_hts,noMaf)
    real(rp), dimension(:,:), pointer :: RADIANCES     ! (Nptg,noChans)
    real(rp), dimension(:,:), pointer :: SPS_PATH ! spsices on path
    real(rp), dimension(:,:), pointer :: TAN_DH_DT ! dH/dT at Tangent
    real(rp), dimension(:,:), pointer :: TAN_TEMP ! ***
    real(rp), dimension(:,:), pointer :: T_GLGRID ! Temp on glGrid surfs

    real(rp), dimension(:,:,:), pointer :: DH_DT_GLGRID ! *****

    type (VectorValue_T), pointer :: EARTHREFL ! Earth reflectivity
    type (VectorValue_T), pointer :: ELEVOFFSET ! Elevation offset
    type (VectorValue_T), pointer :: FIRSTRADIANCE ! Radiance qty for first signal
    type (VectorValue_T), pointer :: LOSVEL ! Line of sight velocity
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: ORBINCLINE ! Orbital inclination
    type (VectorValue_T), pointer :: PTAN ! Tangent pressure component of state vector
    type (VectorValue_T), pointer :: REFGPH ! Reference geopotential height
    type (VectorValue_T), pointer :: SCGEOCALT ! S/C geocentric altitude /m
    type (VectorValue_T), pointer :: SIDEBANDRATIO ! The sideband ratio to use
    type (VectorValue_T), pointer :: SPACERADIANCE ! Emission from space
    type (VectorValue_T), pointer :: TEMP ! Temperature component of state vector
    type (VectorValue_T), pointer :: THISRADIANCE ! A radiance vector quantity

    type (Signal_T) :: FIRSTSIGNAL      ! The first signal we're dealing with

    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS ! ***
    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS_P ! ***
    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS_M ! ***

    type (catalog_T), dimension(:), pointer :: MY_CATALOG ! ***
    type (catalog_T), pointer :: thisCatalogEntry
    type (line_T), pointer :: thisLine

    type (Grids_T) :: Grids_f   ! All the coordinates
    type (Grids_T) :: Grids_dw  ! All the spectroscopy(W) coordinates
    type (Grids_T) :: Grids_dn  ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_dv  ! All the spectroscopy(V) coordinates

    ! ZVI's dumping ground for variables he's too busy to put in the right
    ! place, and doesn't want to write comments for

    ! Local storage places for derivatives..(Temporary..)
    real(r4), dimension(:,:,:,:)  , pointer :: K_TEMP
    real(r4), dimension(:,:,:,:,:), pointer :: K_ATMOS
    real(r4), dimension(:,:,:,:,:), pointer :: K_SPECT_DW
    real(r4), dimension(:,:,:,:,:), pointer :: K_SPECT_DN
    real(r4), dimension(:,:,:,:,:), pointer :: K_SPECT_DV

    ! Executable code --------------------------------------------------------
    ! ------------------------------------------------------------------------

    if ( toggle(emit) ) &
      & call trace_begin ( 'ForwardModel, MAF=', index=fmstat%maf )

    ! Nullify all our pointers!

    nullify ( grids, usedchannels, usedsignals, superset, ext_ind_c, &
      & tan_inds, gl_inds )
    nullify ( gl_ndx, gl_indgen )
    nullify ( do_gl, lin_log )
    nullify ( do_calc_zp, do_calc_dn, do_calc_dv, do_calc_dw, &
      & do_calc_hyd, do_calc_t, do_calc_fzp )
    nullify ( k_temp, k_atmos, k_spect_dw, k_spect_dn, k_spect_dv )
    nullify ( frequencies )
    nullify ( alpha_path_c, del_s, dhdz_path, drad_df, drad_dn, &
      & drad_dt, drad_dv, drad_dw, h_path, incoptdepth, n_path, &
      & one_tan_ht, one_tan_temp, path_dsdh, phi_basis, phi_basis_dn, &
      & phi_basis_dv, phi_basis_dw, phi_path, p_glgrid, p_path, radv, &
      & ref_corr, sps_values, tau, t_path, t_script, xm, ym, zgx,&
      & z_basis, z_basis_dn, z_basis_dv, z_basis_dw, z_glgrid, z_path )

    nullify ( beta_path, beta_path_c, beta_path_f, dbeta_dn_path_c, &
      & dbeta_dn_path_f, dbeta_dt_path_c, dbeta_dt_path_f, &
      & dbeta_dv_path_c, dbeta_dv_path_f, dbeta_dw_path_c, dbeta_dw_path_f, &
      & dhdz_glgrid, dh_dt_path, dx_dt, d2x_dxdt, eta_zp, eta_zxp_dn, &
      & eta_zxp_dv, eta_zxp_dw, eta_zxp_t, h_glgrid, k_atmos_frq, &
      & k_spect_dn_frq, k_spect_dv_frq, k_spect_dw_frq, eta_fzp, &
      & k_temp_frq, ptg_angles, radiances, sps_path, tan_dh_dt, tan_temp, &
      & t_glgrid, dh_dt_glgrid )

    nullify ( lineFlag )

    ! Work out what we've been asked to do -----------------------------------

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

    ! Start sorting out stuff from state vector ------------------------------

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

    ! Think about sidebands
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

    ! Sort out some important dimensions
    noSpecies = size ( fwdModelConf%molecules )
    n_t_phi = temp%template%noInstances
    n_t_zeta = temp%template%noSurfs
    no_tan_hts = fwdModelConf%tangentGrid%noSurfs
    MAF = fmStat%maf
    noMAFs = firstRadiance%template%noInstances
    noMIFs = firstRadiance%template%noSurfs

    !  Get some dimensions that we'll use a lot

    ! Work out the `window' stuff for temperature.
    phiWindow = fwdModelConf%phiWindow
    mafTInstance = FindOneClosestInstance ( temp, firstRadiance, maf )
    windowStart  = max(1, mafTInstance - phiWindow/2)
    windowFinish = min(mafTInstance + phiWindow/2, n_t_phi)

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

    maxNoFSurfs = 0
    do specie = 1, noSpecies
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=fwdModelConf%molecules(specie) )
      maxNoFSurfs = max(maxNoFSurfs, f%template%noSurfs)
    end do

    ! Work out which spectroscopy we're going to need ------------------------
    allocate ( My_Catalog(noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )

    do j = 1, noSpecies
      Spectag = spec_tags(fwdModelConf%molecules(j))
      catIndex = FindFirst(catalog%spec_tag == spectag)

      if ( catIndex < 1 ) then
        call output ( 'No catalog entry for molecule ' )
        call display_string ( lit_indices ( fwdModelConf%molecules(j) ), advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Possible problem with spectroscopy catalog" )
      endif

      thisCatalogEntry => Catalog(catIndex)
      if ( associated ( thisCatalogEntry%lines ) ) then
        ! Now subset the lines according to the signal we're using
        call Allocate_test ( lineFlag, size(thisCatalogEntry%lines), &
          & 'lineFlag', ModuleName )
        lineFlag = .false.
        do k = 1, size ( thisCatalogEntry%lines )
          thisLine => lines(thisCatalogEntry%lines(k))
          do sigInd = 1, size(fwdModelConf%signals)
            if ( associated(thisLine%signals) ) then
              doThis = any ( thisLine%signals == &
                & fwdModelConf%signals(sigInd)%index )
              ! If we're only doing one sideband, maybe we can remove some more lines
              if ( sidebandStart==sidebandStop ) doThis = doThis .and. &
                & any( ( thisLine%sidebands == sidebandStart ) .or. &
                & ( thisLine%sidebands == 0 ) )
            else
              doThis = .true.
            end if
            lineFlag(k) = lineFlag(k) .or. doThis
          end do ! End loop over signals requested in fwm
        end do ! End loop over lines
        My_Catalog(j) = thisCatalogEntry
        nullify ( my_catalog(j)%lines ) ! Don't deallocate it by mistake
        ! Check we have at least one line for this
        if ( count(lineFlag) == 0 ) then
          call get_string ( lit_indices(fwdModelConf%molecules(j)), molName )
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No relevant lines for '//trim(molName) )
        endif
        call Allocate_test ( my_catalog(j)%lines, count(lineFlag),&
          & 'my_catalog(?)%lines', ModuleName )
        my_catalog(j)%lines = pack ( thisCatalogEntry%lines, lineFlag )
        call Deallocate_test ( lineFlag, 'lineFlag', ModuleName )
      else
        ! No lines for this species
        my_catalog(j) = thisCatalogEntry
        call Allocate_test ( my_catalog(j)%lines, 0, 'my_catalog(?)%lines(0)', &
          & ModuleName )
      end if
    end do ! Loop over species

    ! Work out which frequencies we're going to need in non frequency --------
    ! averaging case

    ! Now, allocate other variables we're going to need later ----------------

    call allocate_test ( one_tan_ht, 1, 'one_tan_ht', ModuleName )
    call allocate_test ( one_tan_temp, 1, 'one_tan_temp', ModuleName )

    allocate ( k_temp(noUsedChannels, no_tan_hts, n_t_zeta, &
      & windowStart:windowFinish), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_temp' )
    allocate ( k_atmos(noUsedChannels, no_tan_hts, maxNoFSurfs, &
      & windowStart:windowFinish, noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_atmos' )
    allocate ( k_spect_dw(noUsedChannels, no_tan_hts, maxNoFSurfs, &
      & windowStart:windowFinish, noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dw' )
    allocate ( k_spect_dn(noUsedChannels, no_tan_hts, maxNoFSurfs, &
      & windowStart:windowFinish, noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dn' )
    allocate ( k_spect_dv(noUsedChannels, no_tan_hts, maxNoFSurfs, &
      & windowStart:windowFinish, noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dv' )

    ! Setup our temporary `state vector' like arrays -------------------------
    call load_sps_data ( fwdModelIn, fwdModelExtra, fwdModelConf%molecules, &
     &   p_len, f_len, h2o_ind, lin_log, sps_values, Grids_f, Grids_dw, &
     &   Grids_dn, Grids_dv)

    ! Compute Gauss Legendre (GL) grid ---------------------------------------
    nlvl = size(FwdModelConf%integrationGrid%surfs)
    maxVert = (Nlvl-1) * Ng + Nlvl

    ! Work out the dimensions of it
    NLm1 = nlvl - 1

    ! Allocate GL grid stuff
    call Allocate_test ( xm, NLm1, 'xm', ModuleName )
    call Allocate_test ( ym, NLm1, 'ym', ModuleName )
    call Allocate_test ( zgx, Ngp1, 'zgx', ModuleName )

    call Allocate_test ( z_glGrid, maxVert, 'z_glGrid', ModuleName )
    call Allocate_test ( p_glGrid, maxVert, 'z_glGrid', ModuleName )

    call allocate_test ( h_glgrid, maxVert, n_t_phi, 'h_glgrid', ModuleName )
    call allocate_test ( t_glgrid, maxVert, n_t_phi, 't_glgrid', ModuleName )
    call allocate_test ( dhdz_glgrid, maxVert, n_t_phi, 'dhdz_glgrid', ModuleName )
    call allocate_test ( dh_dt_glgrid, maxVert, n_t_phi, n_t_zeta, &
      &  'dh_dt_glgrid', ModuleName )

    ! From the selected integration grid pressures define the GL pressure
    ! grid:
    zGx(1) = -1.0_rp
    zGx(2:Ngp1) = Gx(1:Ng)
    xm(1:NLm1) = 0.5_rp * ( &
      &   FwdModelConf%integrationGrid%surfs(2:nlvl) + &
      &   FwdModelConf%integrationGrid%surfs(1:NLm1) )
    ym(1:NLm1) = 0.5_rp * ( &
      &   FwdModelConf%integrationGrid%surfs(2:nlvl) - &
      &   FwdModelConf%integrationGrid%surfs(1:NLm1) )
    z_glgrid(1:maxVert-1) = reshape ((spread(xm,1,Ngp1) +          &
      &   spread(ym,1,Ngp1) * spread(zGx,2,NLm1)), &
      &  (/maxVert-1/))
    z_glgrid(maxVert) = FwdModelConf%integrationGrid%surfs(nlvl)
    p_glgrid = 10.0_rp**(-z_glgrid)

    call Deallocate_test ( xm, 'xm', ModuleName )
    call Deallocate_test ( ym, 'ym', ModuleName )
    call Deallocate_test ( zgx,'zgx', ModuleName )

    ! Compute hydrostatic grid -----------------------------------------------

    ! Insert into bill's 2d hydrostatic equation.
    ! The phi input for this program are the orbit plane projected
    ! geodetic locations of the temperature phi basis--not necessarily
    ! the tangent phi's which may be somewhat different.

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_Begin ( 'ForwardModel.Hydrostatic' )
    call two_d_hydrostatic(temp%template%surfs(:,1), &
      &  temp%template%phi(1,:)*Deg2Rad,  &
      &  temp%values, &
      &  spread(refGPH%template%surfs(1,1),1,n_t_phi),&
      &  0.001*refGPH%values(1,:),  &
      &  z_glgrid,orbIncline%values(1,1)*Deg2Rad,    &
      &  t_glgrid, h_glgrid, dhdz_glgrid, &
      &  dh_dt_glgrid )

    call allocate_test ( tan_inds, no_tan_hts, 'tan_inds', ModuleName )
    call Hunt ( z_glgrid-0.0001_rp, &
      & fwdModelConf%tangentGrid%surfs, &
      & tan_inds, allowTopValue=.true. )
    surfaceTangentIndex = count ( tan_inds == 1 )

    elev_offset= 0.0_rp

    call allocate_test ( tan_temp, no_tan_hts, n_t_phi, 'tan_temp', &
      &  ModuleName )

    ! Allocate path quantities -----------------------------------------------

    ! First, Allocate gl_slab arrays....

    no_ele = 2*maxVert     ! maximum possible

    allocate ( gl_slabs ( no_ele, noSpecies), stat=ier )
    if ( ier /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"gl_slabs" )

    do i = 1, NOSPECIES
      nl = size(My_Catalog(i)%Lines)
      gl_slabs(1:no_ele,i)%no_lines = nl
      do j = 1, no_ele
        call AllocateOneSlabs ( gl_slabs(j,i), nl )
      end do
    end do

    if(FwdModelConf%temp_der) then
      allocate(gl_slabs_p(no_ele,NOSPECIES), &
        &  gl_slabs_m(no_ele,NOSPECIES), STAT=ier)
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//"gl_slabs_[pm]" )

      do i = 1, NOSPECIES
        nl = size(My_Catalog(i)%Lines)
        gl_slabs_m(1:no_ele,i)%no_lines = nl
        gl_slabs_p(1:no_ele,i)%no_lines = nl
        do j = 1, no_ele
          call AllocateOneSlabs ( gl_slabs_p(j,i), nl )
          call AllocateOneSlabs ( gl_slabs_m(j,i), nl )
        end do
      end do
    endif

    ! Now allocate all path related... with maximum length..

    brkpt = maxVert
    npc = 2 * (brkpt + Ng) / Ngp1

    ! This can be put outside the mmaf loop

    call Allocate_test ( dhdz_path, no_ele, 'dhdz_path', ModuleName )
    call Allocate_test ( h_path,    no_ele, 'h_path',    ModuleName )
    call Allocate_test ( p_path,    no_ele, 'p_path',    ModuleName )
    call Allocate_test ( path_dsdh, no_ele, 'path_dsdh', ModuleName )
    call Allocate_test ( phi_path,  no_ele, 'phi_path',  ModuleName )
    call Allocate_test ( t_path,    no_ele, 't_path',    ModuleName )
    call Allocate_test ( z_path,    no_ele, 'z_path',    ModuleName )

    call Allocate_Test ( alpha_path_c, npc, 'alpha_path_c', ModuleName )
    call Allocate_Test ( incoptdepth,  npc, 'incoptdept',   ModuleName )
    call Allocate_Test ( t_script,     npc, 't_script',     ModuleName )
    call Allocate_Test ( tau,          npc, 'tau',          ModuleName )
    call Allocate_test ( del_s,        npc, 'del_s',        ModuleName )
    call Allocate_test ( do_gl,        npc, 'do_gl',        ModuleName )
    call Allocate_test ( ext_ind_c,    npc, 'ext_ind_c',    ModuleName )
    call Allocate_test ( n_path,       npc, 'n_path',       ModuleName )
    call Allocate_test ( ref_corr,     npc, 'ref_corr',     ModuleName )

    call Allocate_test ( beta_path_c, npc, noSpecies, 'beta_path_c', ModuleName )
    call Allocate_test ( do_calc_fzp, no_ele, f_len, 'do_calc_fzp', ModuleName )
    call Allocate_test ( do_calc_zp, no_ele, p_len, 'do_calc_zp', ModuleName )
    call Allocate_test ( eta_zp, no_ele, p_len, 'eta_zp', ModuleName )
    call Allocate_test ( eta_fzp, no_ele, f_len, 'eta_zp', ModuleName )
    call Allocate_test ( sps_path, no_ele, noSpecies, 'sps_path', ModuleName )

    if(FwdModelConf%temp_der) then

      ! Allocation for metrics routine when Temp. derivative is needed:

      sv_t_len = n_t_zeta*(windowFinish-windowStart+1)

      call Allocate_test ( dRad_dt, sv_t_len, 'dRad_dt', ModuleName )
      call Allocate_test ( dbeta_dt_path_c, npc, noSpecies, 'dbeta_dt_path_c', ModuleName )
      call Allocate_test ( dh_dt_path, no_ele, sv_t_len, 'dh_dt_path', ModuleName )
      call Allocate_test ( do_calc_hyd, no_ele, sv_t_len, 'do_calc_hyd', ModuleName )
      call Allocate_test ( do_calc_t, no_ele, sv_t_len, 'do_calc_t', ModuleName )
      call Allocate_test ( eta_zxp_t, no_ele, sv_t_len, 'eta_zxp_t', ModuleName )
      call Allocate_test ( tan_dh_dt, 1, n_t_zeta, 'tan_dh_dt', ModuleName )

    endif

    if ( FwdModelConf%atmos_der ) then
      call Allocate_test ( dRad_df, f_len, 'dRad_df', ModuleName )
    end if

    if(FwdModelConf%spect_der) then

      ! Allocation when spectaral derivative are needed:

      call Allocate_test ( dbeta_dw_path_c, npc, noSpecies, &
        & 'dbeta_dw_path_c', ModuleName )
      call Allocate_test ( dbeta_dn_path_c, npc, noSpecies, &
        & 'dbeta_dn_path_c', ModuleName )
      call Allocate_test ( dbeta_dv_path_c, npc, noSpecies, &
        & 'dbeta_dv_path_c', ModuleName )

      sv_dw_len = SUM(Grids_dw%no_z(:) * Grids_dw%no_p(:) * Grids_dw%no_f(:))
      sv_dn_len = SUM(Grids_dn%no_z(:) * Grids_dn%no_p(:) * Grids_dn%no_f(:))
      sv_dv_len = SUM(Grids_dv%no_z(:) * Grids_dv%no_p(:) * Grids_dv%no_f(:))

      call Allocate_test ( do_calc_dw, no_ele, sv_dw_len, 'do_calc_dw', &
                        &  ModuleName )
      call Allocate_test ( do_calc_dn, no_ele, sv_dn_len, 'do_calc_dn', &
                        &  ModuleName )
      call Allocate_test ( do_calc_dv, no_ele, sv_dv_len, 'do_calc_dv', &
                        &  ModuleName )

      call Allocate_test ( eta_zxp_dw, no_ele, sv_dw_len, 'eta_zxp_dw', &
                        &  ModuleName )
      call Allocate_test ( eta_zxp_dn, no_ele, sv_dn_len, 'eta_zxp_dn', &
                        &  ModuleName )
      call Allocate_test ( eta_zxp_dv, no_ele, sv_dv_len, 'eta_zxp_dv', &
                        &  ModuleName )

      call Allocate_test ( drad_dw, sv_dw_len, 'drad_dw', ModuleName )
      call Allocate_test ( drad_dn, sv_dn_len, 'drad_dn', ModuleName )
      call Allocate_test ( drad_dv, sv_dv_len, 'drad_dv', ModuleName )

    endif

    call allocate_test ( ptg_angles, no_tan_hts, noMAFs, 'ptg_angles', &
      &    ModuleName )
    call allocate_test ( radiances, no_tan_hts, noUsedChannels, &
      & 'radiances', ModuleName )
    call allocate_test ( dx_dt, no_tan_hts, n_t_zeta, 'dx_dt', ModuleName )
    call allocate_test ( d2x_dxdt, no_tan_hts, n_t_zeta, 'd2x_dxdt', &
                      &  ModuleName )

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_End ( 'ForwardModel.Hydrostatic' )

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_Begin ( 'ForwardModel.SidebandLoop' )

    ! Loop over sidebands ----------------------------------------------------
    do thisSideband = sidebandStart, sidebandStop, sidebandStep
      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call Trace_Begin ( 'ForwardModel.Sideband ', index=thisSideband )


      ! Work out which pointing frequency grid we're going to need if ----------
      ! frequency averaging

      ! Code splits into two sections, one for when we're doing frequency
      ! averaging, and one when we're not.
      if ( fwdModelConf%do_freq_avg ) then ! --- Doing freq. avg. ---

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

        ! Now we've identified the pointing grids.  Locate the tangent grid
        ! within it.
        call allocate_test ( grids, FwdModelConf%TangentGrid%nosurfs, &
          "Grids", ModuleName )
        call Hunt ( PointingGrids(whichPointingGrid)%oneGrid%height, &
          & FwdModelConf%TangentGrid%surfs, grids, allowTopValue=.true. )
        ! Work out the maximum number of frequencies
        maxNoFreqs = 0
        do ptg_i = 1, no_tan_hts
          maxNoFreqs = max ( maxNoFreqs, &
            & Size(pointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies) )
        end do

      else ! ------------------------- Not frequency averaging ---------

        call allocate_test ( frequencies,noUsedChannels, "frequencies", &
          &   ModuleName )
        do channel = 1, noUsedchannels
          frequencies(channel) = &
            & fwdModelConf%signals(usedSignals(channel))%centerFrequency + &
            & fwdModelConf%signals(usedSignals(channel))% &
            & frequencies(usedChannels(channel))
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
        maxNoFreqs = noUsedChannels
      end if

      call Allocate_test ( radv, maxNoFreqs, 'radV', ModuleName )

      if (fwdModelConf%temp_der) &
        & call Allocate_test (k_temp_frq,maxNoFreqs,sv_t_len,'k_temp_frq', &
        &  ModuleName )

      if (fwdModelConf%atmos_der) &
        & call Allocate_test ( k_atmos_frq,maxNoFreqs,f_len,'k_atmos_frq',&
        &   ModuleName )

      if (fwdModelConf%spect_der) then
        call Allocate_test ( k_spect_dw_frq , maxNoFreqs, sv_dw_len , &
          & 'k_spect_dw_frq', ModuleName )
        call Allocate_test ( k_spect_dn_frq , maxNoFreqs, sv_dn_len , &
          & 'k_spect_dn_frq', ModuleName )
        call Allocate_test ( k_spect_dv_frq , maxNoFreqs, sv_dv_len , &
          & 'k_spect_dv_frq', ModuleName )
      endif

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_Begin ( 'ForwardModel.PointingLoop' )

      ! Loop over pointings --------------------------------------------------
      do ptg_i = 1, no_tan_hts
        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call Trace_Begin ( 'ForwardModel.Pointing ', index=ptg_i )

        ! allocate the path stuff
        brkpt = MaxVert + 1 - tan_inds(ptg_i) ! path tangent index
        no_ele = 2 * brkpt
        npc = 2 * (brkpt + Ng) / Ngp1

        ! This is not pretty but we need some coarse grid extraction indicies
        k = Ngp1
        j = (npc+1)/2
        ext_ind_c(:) = 0
        ext_ind_c(1:npc) = (/(i*k-Ng,i=1,j),((i-1)*k-Ng+1,i=j+1,npc)/)

        ! Compute z_path & p_path
        z_path(1:no_ele) = (/(z_glgrid(i),i=MaxVert,tan_inds(ptg_i),-1), &
          (z_glgrid(i),i=tan_inds(ptg_i),MaxVert)/)
        p_path(1:no_ele) = (/(p_glgrid(i),i=MaxVert,tan_inds(ptg_i),-1), &
          (p_glgrid(i),i=tan_inds(ptg_i),MaxVert)/)

        ! Compute the h_path,t_path,dhdz_path,phi_path,dhdt_path

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_Begin ( 'ForwardModel.MetricsEtc' )

        if (ptg_i < surfaceTangentIndex) then
          neg_tan_ht = temp%values(1,mafTInstance) * &
            &  (fwdModelConf%tangentGrid%surfs(ptg_i) - z_glgrid(1)) / 14.8_rp
          e_rflty = earthRefl%values(1,1)
          phi_tan = firstRadiance%template%phi(1,MAF)

          ! *** This is where we will interpolate Phi_tan

          if(FwdModelConf%temp_der) then
            ! Set up temperature representation basis stuff
            call metrics((/phi_tan/),(/tan_inds(ptg_i)/),                    &
              &  temp%template%geodLat(1,windowStart:windowFinish)*Deg2Rad,  &
              &  z_glgrid,h_glgrid(:,windowStart:windowFinish),              &
              &  t_glgrid(:,windowStart:windowFinish),                       &
              &  dhdz_glgrid(:,windowStart:windowFinish),                    &
              &  orbIncline%values(1,1)*Deg2Rad,h_path(1:no_ele),            &
              &  phi_path(1:no_ele),t_path(1:no_ele),dhdz_path(1:no_ele),    &
              &  Req,TAN_PHI_H_GRID=one_tan_ht,                              &
              &  TAN_PHI_T_GRID=one_tan_temp,                                &
              &  NEG_H_TAN = (/neg_tan_ht/),DHTDTL0=tan_dh_dt,               &
              &  DHIDTLM=dh_dt_glgrid, DHITDTLM=dh_dt_path(1:no_ele,:),      &
              &  Z_BASIS = temp%template%surfs(:,1),                         &
              &  ETA_ZXP=eta_zxp_t(1:no_ele,:),                              &
              &  DO_CALC_T = do_calc_t(1:no_ele,:),                          &
              &  DO_CALC_HYD = do_calc_hyd(1:no_ele,:))
          else
            call metrics((/phi_tan/),(/tan_inds(ptg_i)/),                    &
              &  temp%template%geodLat(1,windowStart:windowFinish)*Deg2Rad,  &
              &  z_glgrid,h_glgrid(:,windowStart:windowFinish),              &
              &  t_glgrid(:,windowStart:windowFinish),                       &
              &  dhdz_glgrid(:,windowStart:windowFinish),                    &
              &  orbIncline%values(1,1)*Deg2Rad,h_path(1:no_ele),            &
              &  phi_path(1:no_ele),t_path(1:no_ele),dhdz_path(1:no_ele),    &
              &  Req,TAN_PHI_H_GRID=one_tan_ht,TAN_PHI_T_GRID=one_tan_temp,  &
              &  NEG_H_TAN = (/neg_tan_ht/))
          endif

          ! Tan heights for a negative tan height from metrics is not correctly working.
        else
          e_rflty = 1.0_rp
          if(FwdModelConf%temp_der) then
            ! Set up temperature representation basis stuff
            call metrics((/phi_tan/),(/tan_inds(ptg_i)/),                    &
              &  temp%template%geodLat(1,windowStart:windowFinish)*Deg2Rad,  &
              &  z_glgrid,h_glgrid(:,windowStart:windowFinish),              &
              &  t_glgrid(:,windowStart:windowFinish),                       &
              &  dhdz_glgrid(:,windowStart:windowFinish),                    &
              &  orbIncline%values(1,1)*Deg2Rad,h_path(1:no_ele),            &
              &  phi_path(1:no_ele),t_path(1:no_ele),dhdz_path(1:no_ele),    &
              &  Req,TAN_PHI_H_GRID=one_tan_ht,TAN_PHI_T_GRID=one_tan_temp,  &
              &  DHTDTL0=tan_dh_dt,DHIDTLM=dh_dt_glgrid,                     &
              &  DHITDTLM=dh_dt_path(1:no_ele,:),                            &
              &  Z_BASIS = temp%template%surfs(:,1),                         &
              &  ETA_ZXP=eta_zxp_t(1:no_ele,:),                              &
              &  DO_CALC_T = do_calc_t(1:no_ele,:),                          &
              &  DO_CALC_HYD = do_calc_hyd(1:no_ele,:))
          else
            call metrics((/phi_tan/),(/tan_inds(ptg_i)/),                    &
              &  temp%template%geodLat(1,windowStart:windowFinish)*Deg2Rad,  &
              &  z_glgrid,h_glgrid(:,windowStart:windowFinish),              &
              &  t_glgrid(:,windowStart:windowFinish),                       &
              &  dhdz_glgrid(:,windowStart:windowFinish),                    &
              &  orbIncline%values(1,1)*Deg2Rad,h_path(1:no_ele),            &
              &  phi_path(1:no_ele),t_path(1:no_ele),dhdz_path(1:no_ele),    &
              &  Req,TAN_PHI_H_GRID=one_tan_ht,TAN_PHI_T_GRID=one_tan_temp)
          endif
        endif

        !  ** Determine the eta_zxp_dw, eta_zxp_dn, eta_zxp_dv
        if(FwdModelConf%spect_der) then
          call eval_spect_path(Grids_dw,z_path(1:no_ele), &
            &  phi_path(1:no_ele),do_calc_dw(1:no_ele,:),         &
            &  eta_zxp_dw(1:no_ele,:))
          call eval_spect_path(Grids_dn,z_path(1:no_ele), &
            &  phi_path(1:no_ele),do_calc_dn(1:no_ele,:),         &
            &  eta_zxp_dn(1:no_ele,:))
          call eval_spect_path(Grids_dv,z_path(1:no_ele), &
            &  phi_path(1:no_ele),do_calc_dv(1:no_ele,:),         &
            &  eta_zxp_dv(1:no_ele,:))
        endif
        tan_temp(ptg_i,maf) = one_tan_temp(1)

        ! Now compute the eta_zp & do_calc_zp (for Zeta & Phi only)
        call comp_eta_docalc_no_frq(Grids_f,z_path(1:no_ele), &
          &  phi_path(1:no_ele),do_calc_zp(1:no_ele,:),eta_zp(1:no_ele,:))

        if (h2o_ind > 0) then
          Frq = 0.0
          call comp_sps_path_frq(Grids_f,Frq,sps_values,eta_zp(1:no_ele,:), &
            &  do_calc_zp(1:no_ele,:),lin_log,sps_path(1:no_ele,:),       &
            &  do_calc_fzp(1:no_ele,:),eta_fzp(1:no_ele,:))
          call refractive_index(p_path(ext_ind_c(1:npc)), &
            &  t_path(ext_ind_c(1:npc)),n_path(1:npc),    &
            &  h2o_path=sps_path((ext_ind_c(1:npc)),h2o_ind))
        else
          call refractive_index(p_path(ext_ind_c(1:npc)), &
            &   t_path(ext_ind_c(1:npc)),n_path(1:npc))
        endif

        if(FwdModelConf%temp_der) then
          call get_chi_angles(0.001*scGeocAlt%values(1,1),n_path(npc/2),    &
          &    one_tan_ht(1),phi_tan,Req,elev_offset,ptg_angles(ptg_i,maf), &
          &    tan_dh_dt(1,:), temp%template%surfs(:,1),                    &
          &    temp%template%geodLat(1,:)*Deg2Rad, one_tan_temp(1),         &
          &    fwdModelConf%tangentGrid%surfs(ptg_i),dx_dt(ptg_i,:),        &
          &    d2x_dxdt(ptg_i,:))
        else
          call get_chi_angles(0.001*scGeocAlt%values(1,1),n_path(npc/2), &
          &    one_tan_ht(1),phi_tan,Req,elev_offset,ptg_angles(ptg_i,maf))
        endif

        call comp_refcor(Req+h_path(ext_ind_c(1:npc)), &
          &   1.0_rp+n_path(1:npc),Req+one_tan_ht(1), &
          &   del_s(1:npc),ref_corr(1:npc))

        ! This only needs to be computed on the gl (not coarse) grid thus there is
        ! some duplication here.
        path_dsdh(2:brkpt-1) = path_ds_dh(Req+h_path(2:brkpt-1), &
          &   Req+one_tan_ht(1))
        path_dsdh(brkpt+2:no_ele-1)=path_ds_dh(Req+h_path(brkpt+2:no_ele-1), &
          &          Req+one_tan_ht(1))

        ! Compute ALL the slabs_prep entities over the path's GL grid for this
        ! pointing & mmaf:

        del_temp = 0.0_rp
        call get_gl_slabs_arrays(my_Catalog,p_path(1:no_ele), &
          &  t_path(1:no_ele),0.001*losVel%values(1,maf),gl_slabs, &
          &  no_ele,del_temp)

        if(FwdModelConf%temp_der) then
          del_temp = 10.0_rp
          call get_gl_slabs_arrays(my_Catalog,p_path(1:no_ele), &
            &  t_path(1:no_ele),0.001*losVel%values(1,maf),gl_slabs_p, &
            &  no_ele,del_temp)
          call get_gl_slabs_arrays(my_Catalog,p_path(1:no_ele), &
            &  t_path(1:no_ele),0.001*losVel%values(1,maf),gl_slabs_m, &
            &  no_ele,-del_temp)
        endif

        ! Work out what frequencies we're using for --------------------------
        ! frequency averaging case for this pointing

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.
        if ( FwdModelConf%do_freq_avg ) then
          frequencies =>  &
            &   PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
          noFreqs = size(frequencies)
        end if ! If not, we dealt with this outside the loop

        ! Loop over frequencies ----------------------------------------------
        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_End ( 'ForwardModel.MetricsEtc' )

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_Begin ( 'ForwardModel.FrequencyLoop' )

        do frq_i = 1, noFreqs
          if ( toggle(emit) .and. levels(emit) > 5 ) &
            & call Trace_Begin ('ForwardModel.Frequency ',index=frq_i)

          Frq = frequencies(frq_i)

          ! Setup path quantities --------------------------------------

          ! Compute the sps_path for this Frequency
          call comp_sps_path_frq(Grids_f,Frq,sps_values,eta_zp(1:no_ele,:), &
               & do_calc_zp(1:no_ele,:),lin_log,sps_path(1:no_ele,:),       &
               & do_calc_fzp(1:no_ele,:),eta_fzp(1:no_ele,:))

          if(FwdModelConf%temp_der  .and. FwdModelConf%spect_der) then

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),     &
              & my_Catalog,gl_slabs,ext_ind_c(1:npc),beta_path_c(1:npc,:), &
              & GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp, &
              & GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp, &
              & DBETA_DT_PATH=dbeta_dt_path_c(1:npc,:),                    &
              & DBETA_DW_PATH=dbeta_dw_path_c(1:npc,:),                    &
              & DBETA_DN_PATH=dbeta_dn_path_c(1:npc,:),                    &
              & DBETA_DV_PATH=dbeta_dv_path_c(1:npc,:) )

          else if(FwdModelConf%temp_der) then

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),     &
              & my_Catalog,gl_slabs,ext_ind_c(1:npc),beta_path_c(1:npc,:), &
              & GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp, &
              & GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp, &
              & DBETA_DT_PATH=dbeta_dt_path_c(1:npc,:))

          else if(FwdModelConf%spect_der) then

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),    &
              & my_Catalog,gl_slabs,ext_ind_c(1:npc),beta_path_c(1:npc,:),&
              & DBETA_DW_PATH=dbeta_dw_path_c(1:npc,:),                   &
              & DBETA_DN_PATH=dbeta_dn_path_c(1:npc,:),                   &
              & DBETA_DV_PATH=dbeta_dv_path_c(1:npc,:) )

          else

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele), &
              &   my_Catalog,gl_slabs,ext_ind_c(1:npc),  &
              &   beta_path_c(1:npc,:))

          endif

          alpha_path_c(1:npc) = SUM(sps_path(ext_ind_c(1:npc),:) &
            &    * beta_path_c(1:npc,:),DIM=2)
          call path_contrib(alpha_path_c(1:npc),del_s(1:npc),e_rflty, &
            &  fwdModelConf%tolerance,tau(1:npc),        &
            &  incoptdepth(1:npc),do_gl(1:npc))

          ! ALLOCATE gl grid beta

          no_gl_ndx = count(do_gl(1:npc))
          j = Ng * no_gl_ndx

          call Allocate_test ( gl_inds, j, 'gl_inds', ModuleName )
          call Allocate_test ( beta_path_f, j, noSpecies, 'beta_path_f', ModuleName )
          call Allocate_test ( gl_indgen, Ng, no_gl_ndx, 'gl_indgen', ModuleName )
          call Allocate_test ( gl_ndx, no_gl_ndx, 2, 'gl_ndx', ModuleName )

          gl_ndx(:,1) = pack((/(i,i=1,npc)/),do_gl(1:npc))

          ! Make (/(j-Ng-1,j=1,Ng)/), (/(j,j=1,Ng)/) parameter variables later on

          do i = 1 , no_gl_ndx
            if(gl_ndx(i,1) > npc/2) then
              gl_ndx(i,2) = 1 - Ng + Ngp1 * (gl_ndx(i,1) - 1)
              gl_indgen(:,i) = (/(j,j=1,Ng)/)
            else
              gl_ndx(i,2) = 1 +      Ngp1 * (gl_ndx(i,1) - 1)
              gl_indgen(:,i) = (/(j-Ng-1,j=1,Ng)/)
            endif
          enddo

          ! compute the gl indicies

          gl_inds = reshape(spread(gl_ndx(1:no_gl_ndx,2),1,Ng) +  &
            &  gl_indgen,(/Ng*no_gl_ndx/))

          j = Ng*no_gl_ndx
          if(FwdModelConf%temp_der  .and. FwdModelConf%spect_der) then

            call Allocate_test ( dbeta_dt_path_f, j, noSpecies, &
              & 'dbeta_dt_path_f', ModuleName )
            call Allocate_test ( dbeta_dw_path_f, j, noSpecies, &
              & 'dbeta_dw_path_f', ModuleName )
            call Allocate_test ( dbeta_dn_path_f, j, noSpecies, &
              & 'dbeta_dn_path_f', ModuleName )
            call Allocate_test ( dbeta_dv_path_f, j, noSpecies, &
              & 'dbeta_dv_path_f', ModuleName )

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),      &
              & my_Catalog,gl_slabs,gl_inds,beta_path_f,                    &
              & GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp,  &
              & GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp,  &
              & DBETA_DT_PATH=dbeta_dt_path_f,DBETA_DW_PATH=dbeta_dw_path_f,&
              & DBETA_DN_PATH=dbeta_dn_path_f,DBETA_DV_PATH=dbeta_dv_path_f)

          else if(FwdModelConf%temp_der) then

            call Allocate_test ( dbeta_dt_path_f, j, noSpecies, &
              & 'dbeta_dt_path_f', ModuleName )
            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),      &
              &   my_Catalog,gl_slabs,gl_inds,beta_path_f,                  &
              &   GL_SLABS_M=gl_slabs_m,T_PATH_M=t_path(1:no_ele)-del_temp, &
              &   GL_SLABS_P=gl_slabs_p,T_PATH_P=t_path(1:no_ele)+del_temp, &
              &   DBETA_DT_PATH=dbeta_dt_path_f)

          else if(FwdModelConf%spect_der) then

            call Allocate_test ( dbeta_dw_path_f, j, noSpecies, &
              & 'dbeta_dw_path_f', ModuleName )
            call Allocate_test ( dbeta_dn_path_f, j, noSpecies, &
              & 'dbeta_dn_path_f', ModuleName )
            call Allocate_test ( dbeta_dv_path_f, j, noSpecies, &
              & 'dbeta_dv_path_f', ModuleName )
            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele),      &
              & my_Catalog,gl_slabs,gl_inds,beta_path_f,                    &
              & DBETA_DW_PATH=dbeta_dw_path_f,DBETA_DN_PATH=dbeta_dn_path_f,&
              & DBETA_DV_PATH=dbeta_dv_path_f)

          else

            call get_beta_path(Frq,p_path(1:no_ele),t_path(1:no_ele), &
              &    my_Catalog,gl_slabs,gl_inds,beta_path_f)

          endif

          ! Compute radiative transfer ---------------------------------------

          call rad_tran(Frq,spaceRadiance%values(1,1),e_rflty,             &
            & z_path(ext_ind_c(1:npc)),t_path(ext_ind_c(1:npc)),            &
            & alpha_path_c(1:npc),ref_corr(1:npc),do_gl(1:npc),             &
            & incoptdepth(1:npc),SUM(sps_path(gl_inds,:)*beta_path_f,DIM=2),&
            & path_dsdh(gl_inds),dhdz_path(gl_inds),t_script(1:npc),        &
            & tau(1:npc),Rad,i_stop)

          RadV(frq_i) = Rad

          ! Compute derivatives if needed ----------------------------------

          if(FwdModelConf%atmos_der) then

            call drad_tran_df(z_path(ext_ind_c(1:npc)),Grids_f%no_z,       &
              &  Grids_f%no_p,lin_log,sps_values,beta_path_c(1:npc,:),     &
              &  eta_fzp(ext_ind_c(1:npc),:),sps_path(ext_ind_c(1:npc),:), &
              &  do_calc_fzp(ext_ind_c(1:npc),:),beta_path_f,              &
              &  eta_fzp(gl_inds,:),sps_path(gl_inds,:),                   &
              &  do_calc_fzp(gl_inds,:),do_gl(1:npc),del_s(1:npc),         &
              &  ref_corr(1:npc),path_dsdh(gl_inds),dhdz_path(gl_inds),    &
              &  t_script(1:npc),tau(1:npc),i_stop,drad_df,ptg_i,frq_i)

            k_atmos_frq(frq_i,:) = drad_df

          endif

          if(FwdModelConf%temp_der) then

            call drad_tran_dt(z_path(ext_ind_c(1:npc)),                      &
              & Req+h_path(ext_ind_c(1:npc)),                                &
              & t_path(ext_ind_c(1:npc)),dh_dt_path(ext_ind_c(1:npc),:),     &
              & alpha_path_c(1:npc),SUM(sps_path(ext_ind_c(1:npc),:)         &
              & * dbeta_dt_path_c(1:npc,:) * beta_path_c(1:npc,:),DIM=2),    &
              & eta_zxp_t(ext_ind_c(1:npc),:),do_calc_t(ext_ind_c(1:npc),:), &
              & do_calc_hyd(ext_ind_c(1:npc),:),del_s(1:npc),ref_corr(1:npc),&
              & Req + one_tan_ht(1),dh_dt_path(brkpt,:),frq,do_gl(1:npc),    &
              & req + h_path(gl_inds),t_path(gl_inds),dh_dt_path(gl_inds,:), &
              & SUM(sps_path(gl_inds,:)*beta_path_f,DIM=2),                  &
              & SUM(sps_path(gl_inds,:)*beta_path_f*dbeta_dt_path_f,DIM=2),  &
              & eta_zxp_t(gl_inds,:),do_calc_t(gl_inds,:),path_dsdh(gl_inds),&
              & dhdz_path(gl_inds),t_script(1:npc),tau(1:npc),i_stop,drad_dt,&
              & ptg_i,frq_i)

            k_temp_frq(frq_i,:) = drad_dt

          endif

          if(FwdModelConf%spect_der) then

            ! Spectroscopic derivative  wrt: W

            call drad_tran_dx(z_path(ext_ind_c(1:npc)),Grids_dw%no_z, &
              &  Grids_dw%no_p,dbeta_dw_path_c(1:npc,:),              &
              &  eta_zxp_dw(ext_ind_c(1:npc),:),sps_path(ext_ind_c(1:npc),:), &
              &  do_calc_dw(ext_ind_c(1:npc),:),dbeta_dw_path_f,      &
              &  eta_zxp_dw(gl_inds,:),sps_path(gl_inds,:),           &
              &  do_calc_dw(gl_inds,:),do_gl(1:npc),del_s(1:npc),     &
              &  ref_corr(1:npc),path_dsdh(gl_inds),dhdz_path(gl_inds), &
              &  t_script(1:npc),tau(1:npc),i_stop,drad_dw,ptg_i,frq_i)

            k_spect_dw_frq(frq_i,:) = drad_dw

            ! Spectroscopic derivative  wrt: N

            call drad_tran_dx(z_path(ext_ind_c(1:npc)),Grids_dn%no_z,   &
              &  Grids_dn%no_p,dbeta_dn_path_c(1:npc,:),                &
              &  eta_zxp_dn(ext_ind_c(1:npc),:),sps_path(ext_ind_c(1:npc),:), &
              &  do_calc_dn(ext_ind_c(1:npc),:),dbeta_dn_path_f,        &
              &  eta_zxp_dn(gl_inds,:),sps_path(gl_inds,:),             &
              &  do_calc_dn(gl_inds,:),do_gl(1:npc),del_s(1:npc),       &
              &  ref_corr(1:npc),path_dsdh(gl_inds),dhdz_path(gl_inds), &
              &  t_script(1:npc),tau(1:npc),i_stop,drad_dn,ptg_i,frq_i)

            k_spect_dn_frq(frq_i,:) = drad_dn

            ! Spectroscopic derivative  wrt: Nu0

            call drad_tran_dx(z_path(ext_ind_c(1:npc)),Grids_dv%no_z,   &
              &  Grids_dv%no_p,dbeta_dv_path_c(1:npc,:),                &
              &  eta_zxp_dv(ext_ind_c(1:npc),:),sps_path(ext_ind_c(1:npc),:), &
              &  do_calc_dv(ext_ind_c(1:npc),:),dbeta_dv_path_f,        &
              &  eta_zxp_dv(gl_inds,:),sps_path(gl_inds,:),             &
              &  do_calc_dv(gl_inds,:),do_gl(1:npc),del_s(1:npc),       &
              &  ref_corr(1:npc),path_dsdh(gl_inds),dhdz_path(gl_inds), &
              &  t_script(1:npc),tau(1:npc),i_stop,drad_dv,ptg_i,frq_i)

            k_spect_dv_frq(frq_i,:) = drad_dv

          endif

          call Deallocate_test ( gl_inds, 'gl_inds', ModuleName )
          call Deallocate_test ( beta_path_f, 'beta_path_f', ModuleName )
          call Deallocate_test ( gl_ndx, 'gl_ndx', ModuleName )
          call Deallocate_Test ( gl_indgen, 'gl_indgen', ModuleName )

          if ( FwdModelConf%temp_der ) &
            & call Deallocate_test (dbeta_dt_path_f, 'dbeta_dt_path_f', &
            & ModuleName )

          if ( FwdModelConf%spect_der ) then
            call Deallocate_test ( dbeta_dw_path_f, 'dbeta_dw_path_f', ModuleName )
            call Deallocate_test ( dbeta_dn_path_f, 'dbeta_dn_path_f', ModuleName )
            call Deallocate_test ( dbeta_dv_path_f, 'dbeta_dv_path_f', ModuleName )
          end if

          ! End of frequency loop ----------------------------------------------

          if ( toggle(emit) .and. levels(emit) > 5 ) &
            & call Trace_End ('ForwardModel.Frequency ',index=frq_i)
        end do            ! End freq. loop

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_End ( 'ForwardModel.FrequencyLoop' )

        ! Work out which channel shape information we're going need ----------


        ! Frequency averaging if needed --------------------------------------

        ! Here we either frequency average to get the unconvolved radiances, or
        ! we just store what we have as we're using delta funciton channels

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call trace_begin ( 'ForwardModel.FrequencyAvg' )
        if ( fwdModelConf%do_freq_avg ) then
          do i = 1, noUsedChannels
            sigInd = usedSignals(i)
            channel = usedChannels(i)
            centerFreq = firstSignal%lo + &
              & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
            shapeInd = MatchSignal ( filterShapes%signal, &
              & fwdModelConf%signals(sigInd), sideband = thisSideband )
            if ( shapeInd == 0 ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              &    "No matching channel shape information" )
            call Freq_Avg ( frequencies, &
              & centerFreq+thisSideband * &
              & FilterShapes(shapeInd)%FilterGrid(channel,:), &
              & FilterShapes(shapeInd)%FilterShape(channel,:),&
              & RadV, noFreqs,  &
              & size(FilterShapes(shapeInd)%FilterGrid(channel,:)), &
              & Radiances(ptg_i,i) )
          end do
        else
          Radiances(ptg_i,1:noFreqs) = RadV(1:noFreqs)
        end if

        ! Frequency averaging of derivatives if needed -----------------------

        ! Frequency Average the temperature derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( fwdModelConf%temp_der ) then
          if ( fwdModelConf%do_freq_avg ) then
            do i = 1, noUsedChannels
              sigInd = usedSignals(i)
              channel = usedChannels(i)
              centerFreq = firstSignal%lo + thisSideband * &
                & fwdModelConf%signals(sigInd)%centerFrequency
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), &
                & sideband = thisSideband, channel=channel )
              sv_i = 1
              do instance = 1, n_t_phi
                do surface = 1, n_t_zeta
                  call Freq_Avg ( frequencies,                        &
                    & centerFreq+thisSideband* &
                    & FilterShapes(shapeInd)%FilterGrid(channel,:), &
                    & FilterShapes(shapeInd)%FilterShape(channel,:),&
                    & k_temp_frq(:,sv_i), noFreqs, &
                    & size(FilterShapes(shapeInd)%FilterGrid(channel,:)),r)
                  k_temp(i,ptg_i,surface,instance) = r
                  sv_i = sv_i + 1
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! Channel loop
          else
            do i = 1, noUsedChannels
              sv_i = 1
              do instance = 1, n_t_phi
                do surface = 1, n_t_zeta
                  k_temp(i,ptg_i,surface,instance) = k_temp_frq(1,sv_i)
                  sv_i = sv_i + 1
                end do
              end do
            end do
          end if
        end if

        ! Frequency Average the atmospheric derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( fwdModelConf%atmos_der ) then

          sv_i = 1
          do specie = 1, noSpecies
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                sv_start = sv_i
                do i = 1, noUsedChannels
                  sv_i = sv_start
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid(channel,:))
                  centerFreq = firstSignal%lo + thisSideband * &
                    & fwdModelConf%signals(sigInd)%centerFrequency
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  do instance = 1, Grids_f%no_p(specie)
                    do surface = 1, Grids_f%no_z(specie)
                      call Freq_Avg ( frequencies,                      &
                        & centerFreq+thisSideband * &
                        & FilterShapes(shapeInd)%FilterGrid(channel,:), &
                        & FilterShapes(shapeInd)%FilterShape(channel,:), &
                        & k_atmos_frq(1:noFreqs,sv_i), noFreqs, j, r )
                      k_atmos(i,ptg_i,surface,instance,specie) = r
                      sv_i = sv_i + 1
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! Else not frequency averaging
                sv_start = sv_i
                do i = 1, noUsedChannels
                  sv_i = sv_start
                  do instance = 1, Grids_f%no_p(specie)
                    do surface = 1, Grids_f%no_z(specie)
                      k_atmos(i,ptg_i,surface,instance,specie) = &
                        k_atmos_frq(1,sv_i)
                      sv_i = sv_i + 1
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over species

        end if                          ! Want derivatives for atmos

        ! Frequency Average the spectroscopic derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( fwdModelConf%spect_der ) then

          !  *** dI/dW

          do specie = 1, noSpecies
            ! ** ZEBUG     if ( fwdModelConf%moleculeSpectDerivatives(specie) ) then
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid(channel,:))
                  centerFreq = firstSignal%lo + thisSideband * &
                    & fwdModelConf%signals(sigInd)%centerFrequency
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dw%no_p(specie)
                    do surface = 1, Grids_dw%no_z(specie)
                      call Freq_Avg ( frequencies,                      &
                        & centerFreq+thisSideband * &
                        & FilterShapes(shapeInd)%FilterGrid(channel,:), &
                        & FilterShapes(shapeInd)%FilterShape(channel,:), &
                        & k_spect_dw_frq(:,sv_i), noFreqs, j, r )
                      k_spect_dw(i,ptg_i,surface,instance,specie) = r
                      sv_i = sv_i + 1
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! Else not frequency averaging
                do i = 1, noUsedChannels
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dw%no_p(specie)
                    do surface = 1, Grids_dw%no_z(specie)
                      k_spect_dw(i,ptg_i,surface,instance,specie) = &
                        k_spect_dw_frq(1,sv_i)
                      sv_i = sv_i + 1
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over species

          !  *** dI/dN

          do specie = 1, noSpecies
            ! ** ZEBUG     if ( fwdModelConf%moleculeSpectDerivatives(specie) ) then
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid(channel,:))
                  centerFreq = firstSignal%lo + thisSideband * &
                    & fwdModelConf%signals(sigInd)%centerFrequency
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dn%no_p(specie)
                    do surface = 1, Grids_dn%no_z(specie)
                      call Freq_Avg ( frequencies,                      &
                        & centerFreq+thisSideband * &
                        & FilterShapes(shapeInd)%FilterGrid(channel,:), &
                        & FilterShapes(shapeInd)%FilterShape(channel,:), &
                        & k_spect_dn_frq(:,sv_i), noFreqs, j, r )
                      k_spect_dn(i,ptg_i,surface,instance,specie) = r
                      sv_i = sv_i + 1
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! Else not frequency averaging
                do i = 1, noUsedChannels
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dn%no_p(specie)
                    do surface = 1, Grids_dn%no_z(specie)
                      k_spect_dn(i,ptg_i,surface,instance,specie) = &
                        k_spect_dn_frq(1,sv_i)
                      sv_i = sv_i + 1
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over species

          !  *** dI/dV

          do specie = 1, noSpecies
            ! ** ZEBUG     if ( fwdModelConf%moleculeSpectDerivatives(specie) ) then
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid(channel,:))
                  centerFreq = firstSignal%lo + thisSideband * &
                    & fwdModelConf%signals(sigInd)%centerFrequency
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dv%no_p(specie)
                    do surface = 1, Grids_dv%no_z(specie)
                      call Freq_Avg ( frequencies,                      &
                        & centerFreq+thisSideband * &
                        & FilterShapes(shapeInd)%FilterGrid(channel,:), &
                        & FilterShapes(shapeInd)%FilterShape(channel,:), &
                        & k_spect_dv_frq(:,sv_i), noFreqs, j, r)
                      k_spect_dv(i,ptg_i,surface,instance,specie) = r
                      sv_i = sv_i + 1
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! Else not frequency averaging
                do i = 1, noUsedChannels
                  if(specie == 1) sv_i = 1
                  do instance = 1, Grids_dv%no_p(specie)
                    do surface = 1, Grids_dv%no_z(specie)
                      k_spect_dv(i,ptg_i,surface,instance,specie) = &
                        k_spect_dv_frq(1,sv_i)
                      sv_i = sv_i + 1
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over species

        end if                        ! Want derivatives for spect

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call trace_end ( 'ForwardModel.FrequencyAvg' )

        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call trace_end ( 'ForwardModel.Pointing ', index=ptg_i )

        ! End of pointing loop -------------------------------------------------
      end do

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_End ( 'ForwardModel.PointingLoop' )

      ! Convolution if needed, or interpolation to ptan ----------------------

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_begin ( 'ForwardModel.Convolution' )

      ! Work out which antenna patterns we're going to need ------------------

      call allocate_test ( superset, size(antennaPatterns), &
        & 'superset', ModuleName )

      do i = 1, noUsedChannels

        channel = usedChannels(i)
        sigInd = usedSignals(i)
        thisRadiance =>  &
          GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=firstSignal%sideband )

        if ( sidebandStart /= sidebandStop ) then   ! We're folding
          thisRatio = sidebandRatio%values(channel,1)
          if ( thisSideband == 1 ) thisRatio = 1.0 - thisRatio
        else                  ! Otherwise, want just unfolded signal
          thisRatio = 1.0
        end if

        ! Here comes the Convolution codes

        if ( FwdModelConf%do_conv ) then

          do j = 1, size(antennaPatterns)
            superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
              & fwdModelConf%signals, sideband=thisSideband, &
              & channel=channel )
          end do

          if ( all( superset < 0 ) ) &
            call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching antenna patterns." )

          maxSuperset = maxval ( superset )
          where ( superset < 0 )
            superset = maxSuperset + 1
          end where
          whichPatternAsArray = minloc ( superset )
          whichPattern = whichPatternAsArray(1)

          center_angle = ptg_angles(surfaceTangentIndex,maf)
          call convolve_all ( fwdModelConf, fwdModelIn, maf, channel, &
            &  windowStart, windowFinish, mafTInstance-windowStart+1, &
            &  temp, ptan, thisRadiance,FwdModelConf%tangentGrid%surfs,&
            &  ptg_angles(:,maf), tan_temp(:,maf), dx_dt, d2x_dxdt,  &
            &  surfaceTangentIndex, center_angle, &
            &  Radiances(:, i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
            &  thisRatio, Jacobian, fmStat%rows,  &
            &  antennaPatterns(whichPattern), ier )
          !??? Need to choose some index other than 1 for AntennaPatterns ???
          if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'convolve_all failed' )

        else          ! No convolution needed ..

          call no_conv_at_all ( fwdModelConf, fwdModelIn, maf, channel,  &
            &     windowStart, windowFinish, temp, ptan, thisRadiance, &
            &     FwdModelConf%tangentGrid%surfs, &
            &     Radiances(:,i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
            &     thisRatio, Jacobian, fmStat%rows )

        end if

      end do                            ! Channel loop

      call deallocate_test ( superset, 'superset', ModuleName )
      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_end ( 'ForwardModel.Convolution' )


      ! Deallocate maxNoFreqs stuff
      call Deallocate_test ( radv, 'radV', ModuleName )

      if ( fwdModelConf%temp_der ) &
        & call Deallocate_test ( k_temp_frq, 'k_temp_frq', ModuleName )

      if (fwdModelConf%atmos_der) &
        & call Deallocate_test ( k_atmos_frq, 'k_atmos_frq', ModuleName )

      if (fwdModelConf%spect_der) then
        call Deallocate_test ( k_spect_dw_frq, 'k_spect_dw_frq', ModuleName )
        call Deallocate_test ( k_spect_dn_frq, 'k_spect_dn_frq', ModuleName )
        call Deallocate_test ( k_spect_dv_frq, 'k_spect_dv_frq', ModuleName )
      endif

      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call trace_end ( 'ForwardModel.sideband ',index=thisSideband )

      ! End of loop over sidebands ---------------------------------------------
    end do
    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_End ( 'ForwardModel.SidebandLoop' )

    !  **** DEBUG Printing cycle ...

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
901   format ( 'ptan\ ',i3.3)
      Print 902,Ptan%values(1:k,maf)
902   format ( 4(4x, f10.7))

      print *
      do i = 1, noUsedChannels
        channel = usedChannels(i)
        print 903, channel, char(92), ptan%template%noSurfs
903     format ( 'ch', i2.2, '_pfa_rad', a1, i3.3 )
        print 905, &
          & firstRadiance%values(i:firstRadiance%template%InstanceLen:25,maf)
905     format ( 4(2x, 1pg15.8) )
      end do

933   format ( 5(2x, 1pg13.6) )
934   format (A,i2.2,A1,i5.5)

    end if

    !  **** End of Printing cycle ...

    ! Now deallocate lots of stuff
    do i = 1, size(my_catalog)
      call Deallocate_test ( my_catalog(i)%lines, 'my_catalog(?)%lines', &
        & ModuleName )
    end do
    deallocate ( my_catalog, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'my_catalog' )

    call Deallocate_test ( usedChannels, 'usedChannels', ModuleName )
    call Deallocate_test ( usedSignals, 'usedSignals', ModuleName )
    call Deallocate_test ( one_tan_ht, 'one_tan_ht', ModuleName )
    call Deallocate_test ( one_tan_temp, 'one_tan_temp', ModuleName )

    deallocate ( k_temp, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_temp' )
    deallocate ( k_atmos, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_atmos' )
    deallocate ( k_spect_dw, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_spect_dw' )
    deallocate ( k_spect_dn, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_spect_dn' )
    deallocate ( k_spect_dv, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_spect_dv' )

    ! DESTROY THE SPS DATA STUFF

    call Deallocate_test ( z_glGrid, 'z_glGrid', ModuleName )
    call Deallocate_test ( p_glGrid, 'z_glGrid', ModuleName )

    call Deallocate_test ( h_glgrid, 'h_glgrid', ModuleName )
    call Deallocate_test ( t_glgrid, 't_glgrid', ModuleName )
    call Deallocate_test ( dhdz_glgrid, 'dhdz_glgrid', ModuleName )
    call Deallocate_test ( dh_dt_glgrid, 'dh_dt_glgrid', ModuleName )

    call Deallocate_test ( tan_inds, 'tan_inds', ModuleName )
    call Deallocate_test ( tan_temp, 'tan_temp', ModuleName )

    call DestroyCompleteSlabs ( gl_slabs )
    if ( fwdModelConf%temp_der ) then
      call DestroyCompleteSlabs ( gl_slabs_p )
      call DestroyCompleteSlabs ( gl_slabs_m )
    end if

    call Deallocate_test ( dhdz_path, 'dhdz_path', ModuleName )
    call Deallocate_test ( h_path,    'h_path',    ModuleName )
    call Deallocate_test ( p_path,    'p_path',    ModuleName )
    call Deallocate_test ( path_dsdh, 'path_dsdh', ModuleName )
    call Deallocate_test ( phi_path,  'phi_path',  ModuleName )
    call Deallocate_test ( t_path,    't_path',    ModuleName )
    call Deallocate_test ( z_path,    'z_path',    ModuleName )

    call Deallocate_test ( alpha_path_c, 'alpha_path_c', ModuleName )
    call Deallocate_test ( incoptdepth,  'incoptdept',   ModuleName )
    call Deallocate_test ( t_script,     't_script',     ModuleName )
    call Deallocate_test ( tau,          'tau',          ModuleName )
    call Deallocate_test ( del_s,        'del_s',        ModuleName )
    call Deallocate_test ( do_gl,        'do_gl',        ModuleName )
    call Deallocate_test ( ext_ind_c,    'ext_ind_c',    ModuleName )
    call Deallocate_test ( n_path,       'n_path',       ModuleName )
    call Deallocate_test ( ref_corr,     'ref_corr',     ModuleName )

    call Deallocate_test ( beta_path_c, 'beta_path_c', ModuleName )
    call Deallocate_test ( do_calc_zp, 'do_calc_zp', ModuleName )
    call Deallocate_test ( do_calc_fzp, 'do_calc_fzp', ModuleName )
    call Deallocate_test ( eta_zp, 'eta_zp', ModuleName )
    call Deallocate_test ( eta_fzp, 'eta_fzp', ModuleName )
    call Deallocate_test ( sps_path, 'sps_path', ModuleName )

    if(FwdModelConf%temp_der) then
      call Deallocate_test ( dRad_dt, 'dRad_dt', ModuleName )
      call Deallocate_test ( dbeta_dt_path_c, 'dbeta_dt_path_c', ModuleName )
      call Deallocate_test ( dh_dt_path, 'dh_dt_path', ModuleName )
      call Deallocate_test ( do_calc_hyd, 'do_calc_hyd', ModuleName )
      call Deallocate_test ( do_calc_t, 'do_calc_t', ModuleName )
      call Deallocate_test ( eta_zxp_t, 'eta_zxp_t', ModuleName )
      call Deallocate_test ( tan_dh_dt, 'tan_dh_dt', ModuleName )
    endif

    if ( FwdModelConf%atmos_der ) then
      call Deallocate_test ( dRad_df, 'dRad_df', ModuleName )
    end if

    if(FwdModelConf%spect_der) then

      call Deallocate_test ( dbeta_dw_path_c, 'dbeta_dw_path_c', ModuleName )
      call Deallocate_test ( dbeta_dn_path_c, 'dbeta_dn_path_c', ModuleName )
      call Deallocate_test ( dbeta_dv_path_c, 'dbeta_dv_path_c', ModuleName )

      call Deallocate_test ( do_calc_dw, 'do_calc_dw', ModuleName )
      call Deallocate_test ( do_calc_dn, 'do_calc_dn', ModuleName )
      call Deallocate_test ( do_calc_dv, 'do_calc_dv', ModuleName )

      call Deallocate_test ( eta_zxp_dw, 'eta_zxp_dw', ModuleName )
      call Deallocate_test ( eta_zxp_dn, 'eta_zxp_dn', ModuleName )
      call Deallocate_test ( eta_zxp_dv, 'eta_zxp_dv', ModuleName )

      call Deallocate_test ( drad_dw, 'drad_dw', ModuleName )
      call Deallocate_test ( drad_dn, 'drad_dn', ModuleName )
      call Deallocate_test ( drad_dv, 'drad_dv', ModuleName )

    endif

    call Deallocate_test ( ptg_angles, 'ptg_angles', ModuleName )
    call Deallocate_test ( radiances, 'radiances', ModuleName )
    call Deallocate_test ( dx_dt, 'dx_dt', ModuleName )
    call Deallocate_test ( d2x_dxdt, 'd2x_dxdt', ModuleName )

    ! Deallocate all variables allocated earlier -----------------------------

    if ( toggle(emit) ) then
      call trace_end ( 'ForwardModel MAF=',fmStat%maf )
    end if

  end subroutine FullForwardModel

 end module FullForwardModel_m

! $Log$
! Revision 2.7  2001/11/02 10:47:57  zvi
! Implementing frequecy grid
!
! Revision 2.6  2001/10/12 20:40:25  livesey
! Moved sideband ratio check
!
! Revision 2.5  2001/10/09 22:39:08  livesey
! Allow for molecules with zero lines.  This may need attention
! from Bill/Zvi later on.
!
! Revision 2.4  2001/10/02 16:51:41  livesey
! Removed fmStat%finished and reordered loops in forward models
!
! Revision 2.3  2001/09/19 04:38:48  livesey
! Lines per band stuff works now
!
! Revision 2.2  2001/09/18 02:04:38  livesey
! Bug fix with signals/spectroscopy interaction
!
! Revision 2.1  2001/09/18 01:23:19  livesey
! Added band discrimination for lines catalog.  Not tested yet.
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.5.2.56  2001/09/14 22:19:39  livesey
! Fixed bug with sv_i in frequency averaging of k_atmos_frq
!
! Revision 1.5.2.55  2001/09/13 19:57:43  livesey
! Fixed a small memory leak
!
! Revision 1.5.2.54  2001/09/13 19:36:27  livesey
! Added some more useful diagnotics/trace statements
!
! Revision 1.5.2.53  2001/09/13 11:18:15  zvi
! Fix temp. derv. bug
!
! Revision 1.5.2.52  2001/09/13 02:03:08  livesey
! Left a print statement in by mistake
!
! Revision 1.5.2.51  2001/09/13 01:48:39  livesey
! Added lots of deallocates, slight problem with temperature derivatives
! and convolution.
!
! Revision 1.5.2.50  2001/09/13 00:42:18  livesey
! Fixed memory leak

! Revision 1.5.2.49  2001/09/12 23:44:18  livesey
! Fixed bug with calling sequence for metrics

! Revision 1.5.2.48  2001/09/12 22:56:05  livesey
! Got rid of print statement.

! Revision 1.5.2.47  2001/09/12 22:46:26  livesey
! Put dimension limit back on dh_dt_path

! Revision 1.5.2.46  2001/09/12 22:38:25  livesey
! Bug fix

! Revision 1.5.2.45  2001/09/12 04:42:32  zvi
! Fixing Conv. bug

! Revision 1.5.2.44  2001/09/12 01:00:34  livesey
! Fixed problem with vmr derivatives

! Revision 1.5.2.43  2001/09/12 00:32:36  livesey
! Fixed printing

! Revision 1.5.2.42  2001/09/12 00:30:00  zvi
! Put printing stmt. back

! Revision 1.5.2.41  2001/09/12 00:12:03  livesey
! Single channel no convolution, derivatives or frequency averaging works

! Revision 1.5.2.40  2001/09/11 21:24:27  livesey
! Interim version

! Revision 1.5.2.39  2001/09/11 20:48:04  livesey
! Added dump and stop statement

! Revision 1.5.2.38  2001/09/11 20:48:22  zvi
! Develop.

! Revision 1.5.2.37  2001/09/11 08:13:33  zvi
! Fixing bugs, adding printing code

! Revision 1.5.2.36  2001/09/11 01:36:54  livesey
! More tidy ups

! Revision 1.5.2.35  2001/09/11 01:27:14  livesey
! It compiles

! Revision 1.5.2.34  2001/09/11 00:50:32  zvi
! Convolution code..done

! Revision 1.5.2.33  2001/09/11 00:01:06  zvi
! adding convolution..incomplete yet

! Revision 1.5.2.32  2001/09/10 23:50:47  livesey
! Added a use statement

! Revision 1.5.2.31  2001/09/10 23:34:51  zvi
! Added freq_avg code..

! Revision 1.5.2.30  2001/09/10 21:07:58  livesey
! Interim

! Revision 1.5.2.29  2001/09/10 21:06:42  livesey
! Interim

! Revision 1.5.2.28  2001/09/10 20:46:42  livesey
! Interim

! Revision 1.5.2.27  2001/09/10 20:24:06  livesey
! More tidying up

! Revision 1.5.2.26  2001/09/10 19:56:36  livesey
! Added call to AllocateOneSlabs

! Revision 1.5.2.25  2001/09/10 19:38:19  livesey
! Tidied up variable definitions a bit.

! Revision 1.5.2.24  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90

! Revision 1.5.2.23  2001/09/09 11:17:43  zvi
! Cleaning up..

! Revision 1.5.2.22  2001/09/09 03:16:48  livesey
! End of working day (ha ha!)

! Revision 1.5.2.21  2001/09/09 03:03:25  livesey
! More

! Revision 1.5.2.20  2001/09/09 03:02:53  zvi
! more definitions ..

! Revision 1.5.2.19  2001/09/09 02:42:52  zvi
! more definitions ..

! Revision 1.5.2.18  2001/09/09 02:16:56  livesey
! More..

! Revision 1.5.2.17  2001/09/09 02:18:03  zvi
! more definirions ..

! Revision 1.5.2.16  2001/09/09 01:57:31  livesey
! Work

! Revision 1.5.2.15  2001/09/09 01:57:41  zvi
! more Allocatios ..

! Revision 1.5.2.14  2001/09/09 01:40:18  livesey
! More work

! Revision 1.5.2.13  2001/09/09 01:41:58  zvi
! Allocatios ..

! Revision 1.5.2.12  2001/09/09 00:43:34  zvi
! more metric work

! Revision 1.5.2.11  2001/09/09 00:09:15  livesey
! Another intermediate

! Revision 1.5.2.10  2001/09/08 23:47:52  livesey
! More updates

! Revision 1.5.2.9  2001/09/08 23:11:23  livesey
! Sent to zvi

! Revision 1.5.2.8  2001/09/08 21:41:14  zvi
! Starting to import codes

! Revision 1.5.2.7  2001/09/07 22:49:10  livesey
! Very intermediate

! Revision 1.5.2.6  2001/09/07 22:34:21  zvi
! change comments..

! Revision 1.5.2.5  2001/09/07 22:18:04  livesey
! More comments

! Revision 1.5.2.4  2001/09/07 20:16:37  livesey
! Changed stuff to lower case

! Revision 1.5.2.3  2001/09/07 20:05:26  livesey
! Cosmetic change again.

! Revision 1.5.2.2  2001/09/07 19:58:49  zvi
! Starting new code developement

! Revision 1.5  2001/07/04 00:34:24  zvi
! Modified & new code(s) for better timing

! Revision 1.4  2001/06/22 02:35:08  zvi
! Fixing some memory leaks..

! Revision 1.3  2001/06/21 15:05:42  livesey
! Gets tolerance from fwdModelConf

! Revision 1.2  2001/06/21 13:07:08  zvi
! Speed enhancement MAJOR update

! Revision 1.1  2001/05/29 22:53:51  livesey
! First version, taken from old ForwardModelInterface.f90
