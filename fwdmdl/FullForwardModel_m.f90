! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullForwardModel_m

  ! This module contains the `full' forward model.

  implicit NONE
  private
  public :: FullForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!-----------------------------------------------------------------------
contains

! ---------------------------------------------  FullForwardModel  -----
  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                             &  FwdModelOut, oldIfm, FmStat, Jacobian )

  ! This is the full radiative transfer forward model, the workhorse
  ! code

    use Allocate_Deallocate, only: ALLOCATE_TEST, deallocate_test
    use AntennaPatterns_m, only: ANTENNAPATTERNS
    use Comp_Eta_Docalc_No_Frq_m, only: Comp_Eta_Docalc_No_Frq
    use Comp_Sps_Path_Frq_m, only: Comp_Sps_Path_Frq
    use Convolve_All_m, only: Convolve_All
    use D_Hunt_m, only: Hunt
    use D_Lintrp_m, only: Lintrp
    use Dump_0, only: Dump
    use Eval_Spect_Path_m, only: Eval_Spect_Path
    use FilterShapes_m, only: FilterShapes
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelIntermediate_t, &
                                    &   ForwardModelStatus_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Freq_Avg_m, only: Freq_Avg
    use Geometry, only: EarthRadA, EarthRadB, SpeedOfLight
    use Get_Beta_Path_m, only: Get_Beta_Path, Beta_Group_T
    use Get_Chi_Angles_m, only: Get_Chi_Angles
    use Get_Chi_Out_m, only: Get_Chi_Out
    use GLnp, only: NG, GX
    use Intrinsic, only: L_TEMPERATURE, L_RADIANCE, L_PHITAN, L_PTAN, &
      & L_ELEVOFFSET, LIT_INDICES, L_ISOTOPERATIO, L_VMR, &
      & L_ORBITINCLINATION, L_SPACERADIANCE, L_EARTHREFL, L_LOSVEL, &
      & L_SCGEOCALT, L_SIDEBANDRATIO, L_NONE, L_CHANNEL, L_REFGPH
    use Load_sps_data_m, ONLY: LOAD_SPS_DATA, Grids_T, destroygrids_t
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATEONESLABS, &
                                &  DESTROYCOMPLETESLABS
    use Make_Z_Grid_M, only: Make_Z_Grid
    use ManipulateVectorQuantities, only: FindInstanceWindow
    use MatrixModule_0, only: M_ABSENT, M_BANDED, M_FULL
    use MatrixModule_1, only: CreateBlock, FindBlock, MATRIX_T
    use Metrics_m, only: Metrics
    use MLSCommon, only: I4, R4, R8, RP, IP, FINDFIRST
    use MLSFiles, only: get_free_lun
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
      & MLSMSG_Error, MLSMSG_Warning
    use MLSNumerics, ONLY: Hunt, InterpolateValues
    use MLSSignals_m, only: Signal_t, MatchSignal, AreSignalsSuperset, Dump, &
                          & GetNameOfSignal
    use Molecules, only: L_EXTINCTION, spec_tags
    use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
    use Output_m, only: OUTPUT
    use PointingGrid_m, only: POINTINGGRIDS
    use RAD_TRAN_M, only: PATH_CONTRIB, RAD_TRAN, DRAD_TRAN_DF,DRAD_TRAN_DT, &
                         &  DRAD_TRAN_DX
    use REFRACTION_M, only: REFRACTIVE_INDEX, COMP_REFCOR, PATH_DS_DH
    use SpectroscopyCatalog_m, only: CATALOG_T, LINE_T, LINES, CATALOG
    use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
    use String_table, only: GET_STRING
    use Toggles, only: Emit, Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use TWO_D_HYDROSTATIC_M, only: Two_D_Hydrostatic
    use Units, only: Deg2Rad
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VALIDATEVECTORQUANTITY, &
                         &   GETVECTORQUANTITYBYTYPE, &
                         &   M_FullDerivatives

    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Define local parameters
    character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "
    integer, parameter :: Ngp1 = Ng+1     ! NG + 1

    ! Now define local variables, group by type and then
    ! alphabetically

    integer :: BRKPT                    ! Index of midpoint of path
    integer :: CHANIND                  ! A 1 based channel index
    integer :: CHANNEL                  ! A Loop counter
    integer :: DIRECTION                ! Direction of channel numbering
    integer :: FRQ_I                    ! Frequency loop index
    integer :: F_LEN                    ! Total number of f's
    integer :: P_LEN                    ! Partial number of f's (No freq.)
    integer :: H2O_IND                  ! Index of h2o inside f array
    integer :: EXT_IND                  ! Index of extinction inside f array
    integer :: IER                      ! Status flag from allocates
    integer :: I_STOP                   ! Upper index for radiance comp.
    integer :: INSTANCE                 ! Loop counter
    integer :: I                        ! Loop index and other uses .
    integer :: J                        ! Loop index and other uses ..
    integer :: K                        ! Loop index and other uses ..
    integer :: L                        ! Loop index and other uses ..
    integer :: M                        ! Loop index and other uses ..
    integer :: JF                       ! Loop index and other uses ..
    integer :: MAF                      ! MAF under consideration
    integer :: MAXNOPTGFREQS            ! Used for sizing arrays
    integer :: MAXNOFFREQS              ! Max. no. frequencies for any molecule
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXSUPERSET              ! Max. value of superset
    integer :: MAXVERT                  ! Number of points in gl grid
    integer :: NL                       ! Number of lines
    integer :: NLM1                     ! Nlvl - 1
    integer :: Nlvl                     ! Size of integration grid
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOSPECIES                ! No. of molecules under consideration
    integer :: NO_MOL                   ! Number of major molecules (NO iso/vib)
    integer :: NOUSEDCHANNELS           ! How many channels are we considering
    integer :: NO_ELE                   ! Length of a gl path
    integer :: NO_GL_NDX                ! Number of GL points to do
    integer :: NO_TAN_HTS               ! Number of tangent heights
    integer :: NPC                      ! Length of coarse path
    integer :: N_T_ZETA                 ! Number of zetas for temperature
    integer :: PTG_I                    ! Loop counter for the pointings
    integer :: SHAPEIND                 ! Index into filter shapes
    integer :: SIDEBANDSTART            ! Loop limit
    integer :: SIDEBANDSTEP             ! Loop step
    integer :: SIDEBANDSTOP             ! Loop limit
    integer :: SIGIND                   ! Signal index, loop counter
    integer :: SPECTAG                  ! A single spectag
    integer :: SURFACETANGENTINDEX      ! Index in tangent grid of earth's
                                        ! surface
    integer :: THISSIDEBAND             ! Loop counter for sidebands
    integer :: WHICHPOINTINGGRID        ! Index into the pointing grids
    integer :: WINDOWFINISH             ! End of temperature `window'
    integer :: WINDOWSTART              ! Start of temperature `window'
    integer :: SPECIE                   ! Loop counter
    integer :: SV_I                     ! Loop index and other uses .
    integer :: F_LEN_DW                 ! Length of DW in vector
    integer :: F_LEN_DN                 ! Length of DN in vector
    integer :: F_LEN_DV                 ! Length of DV in vector
    integer :: SV_START                 ! Temporary sv_i
    integer :: SV_T_LEN                 ! Number of t_phi*t_zeta in the window
    integer :: SURFACE                  ! Loop counter
    integer :: WHICHPATTERN             ! Index of antenna pattern

    logical :: doThis                   ! Flag for lines
    logical :: temp_der, atmos_der, spect_der, ptan_der ! Flags for various derivatives
    logical :: Update                   ! Just update radiances etc.

    character (len=32) :: molName       ! Name of a molecule

    logical :: dummy(2) = (/.FALSE.,.FALSE./)  ! dummy Flag array

    integer, dimension(:), pointer :: GRIDS !Heights in ptgGrid for each tangent
    integer, dimension(:), pointer :: CHANNELORIGINS ! Does this band start at 0 or 1
    integer, dimension(:), pointer :: USEDCHANNELS ! Which channel is this
    integer, dimension(:), pointer :: USEDSIGNALS ! Which signal is this channel from
    integer, dimension(:), pointer :: SUPERSET ! Used for matching signals
    integer, dimension(:), pointer :: INDICES_C ! Indecies on coarse grid
    integer, dimension(:), pointer :: TAN_INDS ! Index of tangent grid into gl grid
    integer, dimension(:), pointer :: GL_INDS ! Index of GL indecies

    integer, dimension(:,:), pointer :: GL_NDX ! Packed Index array of GL intervals
    integer, dimension(:,:), pointer :: GL_INDGEN ! Temp. array of indecies

    logical, dimension(:), pointer :: DO_GL ! GL indicator
    logical, dimension(:), pointer :: LINEFLAG ! Use this line (noLines per species)

    logical, dimension(:,:), pointer :: DO_CALC_ZP    ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_FZP   ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_DN    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DV    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DW    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_HYD   ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_T     ! 'Avoid zeros'

! Array of Flags indicating  which Temp. coefficient to process

    real(r8) :: FRQ                     ! Frequency

    real(rp) :: DEL_TEMP   ! Temp. step-size in evaluation of Temp. power dep.
    real(rp) :: E_RFLTY                 ! Earth reflectivity at given tan. point
    real(rp) :: NEG_TAN_HT              ! GP Height (in KM.) of tan. press.
                                        ! below surface
    real(rp) :: R,R1,R2                 ! real variables for various uses
    real(rp) :: RAD                     ! Radiance
    real(rp) :: REQ                     ! Equivalent Earth Radius
    real(rp) :: Vel_Cor                 ! Velocity correction due to Vel_z
    real(rp) :: THISRATIO               ! A sideband ratio

    real(rp), dimension(1) :: ONE_TAN_HT ! ***
    real(rp), dimension(1) :: ONE_TAN_TEMP ! ***
    real(rp), dimension(:), pointer :: TAN_HTS ! Accumulation of ONE_TAN_HT
    real(rp), dimension(:), pointer :: TAN_TEMPS ! Accumulation of ONE_TAN_TEMP
    real(rp), dimension(:), pointer :: REQS      ! Accumulation of REQ

    real(r8), dimension(:), pointer :: FREQUENCIES ! Frequencies to compute for

    real(rp), dimension(:), pointer :: ALPHA_PATH_C ! coarse grid Sing.
    real(rp), dimension(:), pointer :: DEL_S ! Integration lengths along path
    real(rp), dimension(:), pointer :: DHDZ_PATH ! dH/dZ on path
    real(rp), dimension(:), pointer :: DRAD_DF ! dI/dVmr
    real(rp), dimension(:), pointer :: DRAD_DN ! dI/dN
    real(rp), dimension(:), pointer :: DRAD_DT ! dI/dT
    real(rp), dimension(:), pointer :: DRAD_DV ! dI/dV
    real(rp), dimension(:), pointer :: DRAD_DW ! dI/dW
    real(rp), dimension(:), pointer :: H_PATH ! Heights on path
    real(rp), dimension(:), pointer :: INCOPTDEPTH ! Incremental Optical depth
    real(rp), dimension(:), pointer :: N_PATH ! Refractivity on path
    real(rp), dimension(:), pointer :: PATH_DSDH ! dS/dH on path
    real(rp), dimension(:), pointer :: PHI_BASIS ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DN ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DV ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_BASIS_DW ! phi basis per species
    real(rp), dimension(:), pointer :: PHI_PATH ! Phi's on path
    real(rp), dimension(:), pointer :: P_GLGRID ! Pressure on glGrid surfs
    real(rp), dimension(:), pointer :: P_PATH ! Pressure on path
    real(rp), dimension(:), pointer :: RADV ! Radiances for 1 pointing on
                                            ! Freq_Grid
    real(rp), dimension(:), pointer :: REF_CORR ! Refraction correction
    real(rp), dimension(:), pointer :: TAU ! Optical depth
    real(rp), dimension(:), pointer :: TAN_TEMP ! ***
    real(rp), dimension(:), pointer :: T_PATH ! Temperatures on path
    real(rp), dimension(:), pointer :: T_SCRIPT ! ********
    real(rp), dimension(:), pointer :: Z_BASIS !zeta basis per specie (n_f_zeta)
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
    real(rp), dimension(:,:), pointer :: DX_DT       ! (No_tan_hts, nz*np)
    real(rp), dimension(:,:), pointer :: D2X_DXDT    ! (No_tan_hts, nz*np)
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
    real(rp), dimension(:),   pointer :: PTG_ANGLES ! (no_tan_hts)
    real(rp), dimension(:,:), pointer :: RADIANCES     ! (Nptg,noChans)
    real(rp), dimension(:,:), pointer :: SPS_PATH ! spsices on path
    REAL(rp), DIMENSION(:,:), POINTER :: TAN_DH_DT ! dH/dT at Tangent
    real(rp), dimension(:,:), pointer :: T_GLGRID ! Temp on glGrid surfs

    real(rp), dimension(:,:,:), pointer :: DH_DT_GLGRID ! *****

! Some declarations by bill

    integer(ip) :: sps_i  ! a species counter
    integer(ip) :: no_sv_p_t ! number of phi basis for temperature
    integer(ip) :: beg_ind, end_ind, beg_ind_z, end_ind_z
    integer(ip) :: beg_ind_p, end_ind_p

    real(rp) :: surf_angle(1), one_dhdz(1), one_dxdh(1)
    real(rp) :: earthradc ! minor axis of orbit plane projected Earth ellipse

    integer(ip), dimension(:), pointer :: rec_tan_inds ! recommended tangent
!                        point indecies from make_z_grid
    real(rp), dimension(:), pointer :: z_tmp ! temporary zeta storage
    real(rp), dimension(:), pointer :: z_all ! mass storage of representation
!                                      bases for z_grid determination
    real(rp), dimension(:), pointer :: z_psig(:) ! recommended PSIG for
!                                      radiative transfer calculations
! THIS VARIABLE REPLACES FwdModelConf%integrationGrid%surfs

    real(rp), dimension(:), pointer :: tan_chi_out
    real(rp), dimension(:), pointer :: dx_dh_out
    real(rp), dimension(:), pointer :: dhdz_out
    real(rp), dimension(:), pointer :: req_out
    real(rp), dimension(:), pointer :: tan_press
    real(rp), dimension(:), pointer :: tan_phi
    real(rp), dimension(:), pointer :: est_scgeocalt
    real(rp), dimension(:,:), pointer :: dxdt_tan
    real(rp), dimension(:,:), pointer :: d2xdxdt_tan
    real(rp), dimension(:,:), pointer :: dxdt_surface
    real(rp), dimension(:,:), pointer :: d2xdxdt_surface
    real(rp), dimension(:,:), pointer :: tan_d2h_dhdt
    real(rp), dimension(:,:,:), pointer :: ddhidhidtl0

! THIS VARIABLE REPLACES fwdModelConf%tangentGrid%surfs

    type (VectorValue_T), pointer :: EARTHREFL ! Earth reflectivity
    type (VectorValue_T), pointer :: ELEVOFFSET ! Elevation offset
    type (VectorValue_T), pointer :: FIRSTRADIANCE ! Radiance qty for first signal
    type (VectorValue_T), pointer :: LOSVEL ! Line of sight velocity
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: ORBINCLINE ! Orbital inclination
    type (VectorValue_T), pointer :: PHITAN ! Tangent geodAngle component of state vector
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

    type (Grids_T) :: Grids_f   ! All the coordinates for VMR
    type (Grids_T) :: Grids_dw  ! All the spectroscopy(W) coordinates
    type (Grids_T) :: Grids_dn  ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_dv  ! All the spectroscopy(V) coordinates
    type (Grids_T) :: Grids_tmp ! All the coordinates for TEMP

    ! ZVI's dumping ground for variables he's too busy to put in the right
    ! place, and doesn't want to write comments for

    ! Local storage places for derivatives..(Temporary..)
    real(r4), dimension(:,:,:,:)  , pointer :: K_TEMP
    real(r4), dimension(:,:,:), pointer :: K_ATMOS
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DW
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DN
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DV

!  The 'all_radiometers grid file' approach variables declaration:

    real(rp) :: max_ch_freq_grid, min_ch_freq_grid

! *** Beta & Molecules grouping variables:
    real(rp) :: beta_ratio
    integer, dimension(:), pointer :: mol_cat_index

    type (beta_group_T), dimension(:), pointer :: beta_group

    ! Executable code --------------------------------------------------------
    ! ------------------------------------------------------------------------

!   Print *, '** Enter ForwardModel, MAF =',fmstat%maf   ! ** ZEBUG

    temp_der = present ( jacobian ) .and. FwdModelConf%temp_der
    atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der

! ** Re-instate when appropriate code is done
!   spect_der = FwdModelConf%spect_der
    spect_der = present ( jacobian ) .and. .false.    ! ** ZEBUG

    if ( toggle(emit) ) & ! set by -f command-line switch
      & call trace_begin ( 'ForwardModel, MAF=', index=fmstat%maf )

    ! Nullify all our pointers!

    nullify ( alpha_path_c, beta_path, beta_path_c, beta_path_f, &
      & channelOrigins, d2x_dxdt, d2xdxdt_surface, d2xdxdt_tan, &
      & dbeta_dn_path_c, dbeta_dn_path_f, dbeta_dt_path_c, dbeta_dt_path_f, &
      & dbeta_dv_path_c, dbeta_dv_path_f, dbeta_dw_path_c, dbeta_dw_path_f, &
      & ddhidhidtl0, del_s, dh_dt_glgrid, dh_dt_path, dhdz_glgrid, &
      & dhdz_out, dhdz_path, do_calc_dn, do_calc_dv, do_calc_dw, &
      & do_calc_fzp, do_calc_hyd, do_calc_t, do_calc_zp, do_gl, drad_df, &
      & drad_dn, drad_dt, drad_dv, drad_dw, dx_dh_out, dx_dt, &
      & dxdt_surface, dxdt_tan, est_scgeocalt, eta_fzp, eta_zp, &
      & eta_zxp_dn, eta_zxp_dv, eta_zxp_dw, eta_zxp_t, frequencies, &
      & gl_indgen, gl_inds, gl_ndx, grids, h_glgrid, h_path, incoptdepth, &
      & indices_c, k_atmos, k_atmos_frq, k_spect_dn, k_spect_dn_frq, &
      & k_spect_dv, k_spect_dv_frq, k_spect_dw, k_spect_dw_frq, k_temp, &
      & k_temp_frq, lineFlag, n_path, path_dsdh, p_glgrid, phi_basis, &
      & phi_basis_dn, phi_basis_dv, phi_basis_dw, phi_path, p_path, &
      & ptg_angles, radiances, RadV, ref_corr, req_out, sps_path, &
      & superset, tan_chi_out, tan_d2h_dhdt, tan_dh_dt, tan_inds, &
      & tan_phi, tan_press, tan_temp, tau, t_glgrid, t_path, t_script, &
      & usedchannels, usedsignals, z_all, z_basis, z_basis_dn, &
      & z_basis_dv, z_basis_dw, z_glgrid, z_path, z_tmp, tan_temps, &
      & tan_hts, reqs )

    ! Work out what we've been asked to do -----------------------------------

    ! Identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    firstSignal = fwdModelConf%signals(1)

    ! Now make sure all the signals we're dealing with are same module,
    ! radiometer and sideband.
    if ( any( fwdModelConf%signals%sideband /= firstSignal%sideband ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed sidebands in forward model config" )
    if ( any( fwdModelConf%signals%radiometer /= firstSignal%radiometer ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed radiometers in forward model config" )

    ! Now from that we identify the radiance quantity we'll be outputting
    firstRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal=firstSignal%index, sideband=firstSignal%sideband )

    ! Start sorting out stuff from state vector ------------------------------

    ! Identify the appropriate state vector components, save vmrs for later
    temp => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature, config=fwdModelConf )
    ptan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=firstSignal%instrumentModule, &
      & foundInFirst=ptan_der, config=fwdModelConf )
    phitan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, instrumentModule=firstSignal%instrumentModule, &
      & config=fwdModelConf )
    elevOffset => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_elevOffset, radiometer=firstSignal%radiometer, &
      & config=fwdModelConf )
    orbIncline => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_orbitInclination, config=fwdModelConf )
    spaceRadiance => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_spaceRadiance, config=fwdModelConf )
    earthRefl => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_earthRefl, config=fwdModelConf )
    refGPH => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_refGPH, config=fwdModelConf )
    losVel => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_losVel, instrumentModule=firstSignal%instrumentModule, &
      & config=fwdModelConf )
    scGeocAlt => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scGeocAlt )

    ! We won't seek for molecules here as we can't have an array of pointers.
    ! When we do want molecule i we would do something like
    ! vmr => GetQuantityForForwardModel (fwdModelIn, fwdModelExtra, &
    !   quantityType=l_vmr, molecule=fwdModelConf.molecules(i))

    ! Now we're going to validate the quantities we've been given; don't
    ! forget we already know what their quantityType's are as that's how we
    ! found them, so we don't need to check that.
    if ( .not. ValidateVectorQuantity(temp, stacked=.TRUE., coherent=.TRUE., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.TRUE., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )
    if ( .not. ValidateVectorQuantity(phitan, minorFrame=.TRUE., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'phitan' )
    if ( .not. ValidateVectorQuantity(elevOffset, verticalCoordinate=(/l_none/), &
      & frequencyCoordinate=(/l_none/), noInstances=(/1/)) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'elevOffset' )
    ! There will be more to come here.

    ! Think about sidebands
    if ( ( fwdModelConf%signals(1)%sideband == 0 ) .and.&
      &  ( fwdModelConf%signals(1)%singleSideband == 0 ) ) then
      ! Do a folded measurement
      sidebandStart = -1
      sidebandStop = 1
      sidebandStep = 2
    else
      ! It's either a single sideband radiometer, or the user requested a
      ! specific sideband.
      ! Check sanity, if they are both non zero they should be the same.
      if ( ( fwdModelConf%signals(1)%singleSideband /= 0 ) .and. &
        &  ( fwdModelConf%signals(1)%sideband /= 0 ) .and. &
        &  ( fwdModelConf%signals(1)%singleSideband /= &
        &    fwdModelConf%signals(1)%sideband ) ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, &
        & "User requested a sideband that doesn't exist" )
      ! OK, use whichever one is given
      if ( fwdModelConf%signals(1)%singleSideband /= 0 ) then
        sidebandStart = fwdModelConf%signals(1)%singleSideband
      else
        sidebandStart = fwdModelConf%signals(1)%sideband
      end if
      sidebandStop = sidebandStart
      sidebandStep = 1
    end if

 ! Sort out some important dimensions
    noSpecies = size ( fwdModelConf%molecules )
    no_mol = count ( fwdModelConf%molecules > 0)
    n_t_zeta = temp%template%noSurfs

    MAF = fmStat%maf

    Vel_Cor = 1.0_rp - losvel%values(1,maf) / speedOfLight

! Sort out a remaining flag
    ptan_der = ptan_der .and. present ( jacobian )

! Work out the `window' stuff for temperature. Create the Grids_tmp structure:

    call findInstanceWindow ( temp, phitan, maf, fwdModelConf%phiWindow, &
                            & windowStart, windowFinish )

    no_sv_p_t = windowFinish-windowStart+1
    sv_t_len = n_t_zeta * no_sv_p_t

    call allocate_test ( grids_tmp%no_z, 1, 'Grids_tmp%no_z', moduleName )
    call allocate_test ( grids_tmp%no_p, 1, 'Grids_tmp%no_p', moduleName )
    call allocate_test ( grids_tmp%no_f, 1, 'Grids_tmp%no_f', moduleName )
    call allocate_test ( grids_tmp%windowstart, 1, 'Grids_tmp%windowstart', &
                       & moduleName )
    call allocate_test ( grids_tmp%windowfinish, 1, 'Grids_tmp%windowfinish', &
                       & moduleName )
    call allocate_test ( grids_tmp%lin_log, 1, 'lin_log', moduleName )

    grids_tmp%no_f = 1
    grids_tmp%no_z = n_t_zeta
    grids_tmp%no_p = no_sv_p_t
    grids_tmp%lin_log = .false.
    grids_tmp%windowStart(1) = windowStart
    grids_tmp%windowFinish(1) = windowFinish

! Allocate space for the zeta, phi & freq. basis componenets

    k = sv_t_len
    call allocate_test ( grids_tmp%zet_basis, n_t_zeta, 'Grids_tmp%zet_basis', &
                       & moduleName )
    call allocate_test ( grids_tmp%phi_basis, no_sv_p_t, 'Grids_tmp%phi_basis', &
                       & moduleName )
    call allocate_test ( grids_tmp%frq_basis, 1, 'Grids_tmp%frq_basis', moduleName )
    call allocate_test ( grids_tmp%values, k, 'Grids_tmp%values', moduleName )
    call allocate_test ( grids_tmp%deriv_flags, k, 'Grids_tmp%deriv_flags', &
                       & moduleName )

    grids_tmp%frq_basis = 0.0
    grids_tmp%zet_basis = temp%template%surfs(1:n_t_zeta,1)
    grids_tmp%phi_basis = temp%template%phi(1,windowStart:windowFinish)*Deg2Rad

    grids_tmp%values = reshape(temp%values(:,windowStart:windowFinish),(/k/))

! ** Initialize to ALL derivatives flags to TRUE :
    grids_tmp%deriv_flags(1:k) = .TRUE.

! ** Now Load the Temp. derivative coeff. flag according to the L2CF

    if ( associated(temp%mask) ) Grids_tmp%deriv_flags = RESHAPE((IAND( &
       & M_FullDerivatives,ICHAR(temp%mask(:,WindowStart:WindowFinish)))==0), &
       & (/sv_t_len/))

! Work out which channels are used, also check we have radiances for them.

    noUsedChannels = 0
    do sigInd = 1, size(fwdModelConf%signals)
      thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
        & signal=fwdModelConf%signals(sigInd)%index, sideband=firstSignal%sideband )
      if ( .not. ValidateVectorQuantity(thisRadiance, minorFrame=.TRUE.,&
        & frequencyCoordinate=(/l_channel/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'radiance' )
      noUsedChannels = noUsedChannels + &
        & count( fwdModelConf%signals(sigInd)%channels )
    end do
    call allocate_test ( usedChannels, noUsedChannels, &
      & 'usedChannels', moduleName )
    call allocate_test ( channelOrigins, noUsedChannels, &
      & 'channelOrigins', moduleName )
    call allocate_test ( usedSignals, noUsedChannels, &
      & 'usedSignals', moduleName )
    channel = 1
    do sigInd = 1, size(fwdModelConf%signals)
      do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
        if ( fwdModelConf%signals(sigInd)%channels(i) ) then
          channelOrigins(channel) = &
            & lbound ( fwdModelConf%signals(sigInd)%frequencies, 1 )
          usedChannels(channel) = i + channelOrigins(channel) - 1
          usedSignals(channel) = sigInd
          channel = channel + 1
        end if
      end do
    end do

    maxNoFFreqs = 0
    maxNoFSurfs = 0
    do specie = 1, noSpecies
      l = fwdModelConf%molecules(specie)
      if ( l > 0 ) then
        f => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_vmr, molIndex=specie, config=fwdModelConf, &
          & radiometer=firstSignal%radiometer )
        maxNoFFreqs = max(maxNoFFreqs, f%template%noChans)
        maxNoFSurfs = max(maxNoFSurfs, f%template%noSurfs)
      end if
    end do

    allocate ( mol_cat_index(no_mol), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'mol_cat_index' )

    mol_cat_index(:) = 0
    mol_cat_index = PACK((/(i,i=1,noSpecies)/),fwdModelConf%molecules > 0)

    allocate ( beta_group(no_mol), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'beta_group' )

    k = max(1,noSpecies-no_mol)
    do i = 1, no_mol
      allocate ( beta_group(i)%cat_index(k), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'beta_group%cat_index' )
      allocate ( beta_group(i)%ratio(k), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'beta_group%ratio' )
      beta_group(i)%n_elements = 0
      beta_group(i)%ratio(1:k) = 0.0
      beta_group(i)%cat_index(1:k) = 0
    end do

    if ( all ( fwdModelConf%molecules > 0 ) ) then

      sv_i = 0
      beta_ratio = 1.0_rp   ! Always, for single element (no grouping)
      do j = 1, noSpecies
        l = fwdModelConf%molecules(j)
!        if ( l == l_extinction ) CYCLE
        sv_i = sv_i + 1
        beta_group(sv_i)%n_elements   = 1
        beta_group(sv_i)%cat_index(1) = j
        beta_group(sv_i)%ratio(1)     = beta_ratio
      end do

    else

      k = noSpecies
      call allocate_test ( gl_inds, k+1, 'gl_inds', moduleName )
      gl_inds(1:k) = fwdModelConf%molecules(1:k)
      gl_inds(k+1) = k

      sv_i = 0
      do j = 1, noSpecies
        k = gl_inds(j)
        l = abs(k)
!        if ( l == l_extinction ) CYCLE
        beta_ratio = 1.0_rp
        if ( k > 0 ) then
          if ( gl_inds(j+1) > 0 ) then
            sv_i = sv_i + 1
            beta_group(sv_i)%n_elements   = 1
            beta_group(sv_i)%cat_index(1) = j
            beta_group(sv_i)%ratio(1)     = beta_ratio
          end if
        else
          if ( gl_inds(j-1) > 0) sv_i = sv_i + 1
          f => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_isotoperatio, molecule=l, noError=.TRUE., &
              & config=fwdModelConf )
          if ( associated ( f ) ) beta_ratio = f%values(1,1)
          i = beta_group(sv_i)%n_elements + 1
          beta_group(sv_i)%n_elements   = i
          beta_group(sv_i)%cat_index(i) = j
          beta_group(sv_i)%ratio(i)     = beta_ratio
        end if
      end do

      call deallocate_test ( gl_inds, 'gl_inds', moduleName )

    end if

! Work out which spectroscopy we're going to need ------------------------

    nullify ( my_catalog )
    allocate ( My_Catalog(noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )

    do j = 1, noSpecies
      ! Skip if the next molecule is negative (indicates that this one is a
      ! parent)
      if ( (j < noSpecies) .and. (fwdModelConf%molecules(j) > 0) ) then
        if ( fwdModelConf%molecules(j+1) < 0 ) then
          nullify ( my_catalog(j)%lines ) ! Don't deallocate it by mistake
          call allocate_test ( my_catalog(j)%lines, 0, &
                            & 'my_catalog(?)%lines(0)', moduleName )
          cycle
        end if
      end if
      l=abs(fwdModelConf%molecules(j))
      Spectag = spec_tags(l)
      thisCatalogEntry => Catalog(FindFirst(catalog%spec_tag == spectag ) )
      if ( associated ( thisCatalogEntry%lines ) ) then
        ! Now subset the lines according to the signal we're using
        call allocate_test ( lineFlag, size(thisCatalogEntry%lines), &
                         &  'lineFlag', moduleName )
        lineFlag = .FALSE.
        do k = 1, size ( thisCatalogEntry%lines )
          thisLine => lines(thisCatalogEntry%lines(k))
          if ( associated(thisLine%signals) ) then
            do sigInd = 1, size(fwdModelConf%signals)
              doThis = any ( thisLine%signals == &
                & fwdModelConf%signals(sigInd)%index )

! If we're only doing one sideband, maybe we can remove some more lines

              if ( sidebandStart==sidebandStop ) doThis = doThis .and. &
                & any( ( thisLine%sidebands == sidebandStart ) .or. &
                & ( thisLine%sidebands == 0 ) )
              lineFlag(k) = lineFlag(k) .or. doThis
            end do ! End loop over signals requested in fwm
          end if
        end do               ! End loop over lines

        My_Catalog(j) = thisCatalogEntry
        nullify ( my_catalog(j)%lines ) ! Don't deallocate it by mistake

! Check we have at least one line for this

        if ( count(lineFlag) == 0 .and. all ( my_catalog(j)%continuum == 0.0 ) ) then
          call get_string ( lit_indices(l), molName )
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No relevant lines or continuum for '//trim(molName) )
        end if
        call allocate_test ( my_catalog(j)%lines, count(lineFlag),&
          & 'my_catalog(?)%lines', moduleName )
        my_catalog(j)%lines = pack ( thisCatalogEntry%lines, lineFlag )
        call deallocate_test ( lineFlag, 'lineFlag', moduleName )

      else

        ! No lines for this species
        my_catalog(j) = thisCatalogEntry
        nullify ( my_catalog(j)%lines ) ! Don't deallocate it by mistake
        call allocate_test ( my_catalog(j)%lines, 0, 'my_catalog(?)%lines(0)', &
          & moduleName )
      end if

    end do         ! Loop over species

! Now, allocate other variables we're going to need later ----------------

! Setup our temporary `state vector' like arrays -------------------------

    if ( spect_der ) then
      call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
       &   firstSignal%radiometer, mol_cat_index, p_len, f_len, h2o_ind,     &
       &   ext_ind, Grids_f, f_len_dw, Grids_dw, f_len_dn, Grids_dn,         &
       &   f_len_dv, Grids_dv )
    else
      call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
       &   firstSignal%radiometer, mol_cat_index, p_len, f_len, h2o_ind,     &
       &   ext_ind, Grids_f )
    end if

! set up output pointing angles------------------------------------------
! note we have to compute req !!!!!!!

! compute equivalent earth radius at phi_t(1), nearest surface

    call allocate_test ( req_out, ptan%template%nosurfs, 'req_out', &
                       & moduleName )
    earthradc = earthRadA*earthRadB / &
          & SQRT((earthRadA**2-earthRadB**2) * &
                 &   SIN(orbIncline%values(1,maf)*Deg2Rad)**2 + &
                 & earthRadB**2)
    req_out = COS(phitan%values(:,maf)*Deg2Rad)**2
    req_out = 0.001_rp*SQRT((earthrada**4*(1.0_rp-req_out) &
    & + earthradc**4*req_out) / (earthrada**2*req_out &
    & + earthradc**2*(1.0-req_out)))
    call allocate_test ( tan_chi_out,ptan%template%nosurfs, 'tan_chi_out', &
                       & moduleName )
    call allocate_test ( dx_dh_out,ptan%template%nosurfs, 'dx_dh_out', &
                       & moduleName )
    call allocate_test ( dhdz_out,ptan%template%nosurfs, 'dhdz_out', &
                       & moduleName )
    if ( h2o_ind > 0 .and. .not. temp_der ) then
      end_ind_z = sum(grids_f%no_z(1:h2o_ind))
      beg_ind_z = end_ind_z - grids_f%no_z(h2o_ind) + 1
      end_ind_p = sum(grids_f%no_p(1:h2o_ind))
      beg_ind_p = end_ind_p - grids_f%no_p(h2o_ind) + 1
      end_ind = sum(grids_f%no_z(1:h2o_ind) * grids_f%no_p(1:h2o_ind) * &
                 &  grids_f%no_f(1:h2o_ind))
      beg_ind = end_ind - grids_f%no_z(h2o_ind)*grids_f%no_p(h2o_ind) * &
                       &  grids_f%no_f(h2o_ind) + 1
      call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*Deg2Rad, &
         & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp, &
         & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
         & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
         & req_out, tan_chi_out, dhdz_out, dx_dh_out, &
         & h2o_zeta_basis=grids_f%zet_basis(beg_ind_z:end_ind_z), &
         & h2o_phi_basis=grids_f%phi_basis(beg_ind_p:end_ind_p), &
         & h2o_coeffs=grids_f%values(beg_ind:end_ind), &
         & lin_log=grids_f%lin_log(h2o_ind))
    else if ( h2o_ind == 0 .and. .not. temp_der ) then
      call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad, &
         & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp, &
         & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
         & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
         & req_out, tan_chi_out, dhdz_out, dx_dh_out )
    else if ( h2o_ind > 0 .and.  temp_der ) then
      call allocate_test ( dxdt_tan,ptan%template%nosurfs,sv_t_len, &
                         & 'dxdt_tan',moduleName )
      call allocate_test ( d2xdxdt_tan,ptan%template%nosurfs,sv_t_len, &
                         & 'd2xdxdt_tan',moduleName )
      end_ind_z = SUM(grids_f%no_z(1:h2o_ind))
      beg_ind_z = end_ind_z - grids_f%no_z(h2o_ind) + 1
      end_ind_p = SUM(grids_f%no_p(1:h2o_ind))
      beg_ind_p = end_ind_p - grids_f%no_p(h2o_ind) + 1
      end_ind = SUM(grids_f%no_z(1:h2o_ind)*grids_f%no_p(1:h2o_ind) * &
                  & grids_f%no_f(1:h2o_ind))
      beg_ind = end_ind - grids_f%no_z(h2o_ind)*grids_f%no_p(h2o_ind) * &
                     &  grids_f%no_f(h2o_ind) + 1
      call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad, &
         & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp, &
         & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
         & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
         & req_out, tan_chi_out, dhdz_out, dx_dh_out, &
         & h2o_zeta_basis=grids_f%zet_basis(beg_ind_z:end_ind_z), &
         & h2o_phi_basis=grids_f%phi_basis(beg_ind_p:end_ind_p), &
         & h2o_coeffs=grids_f%values(beg_ind:end_ind), &
         & lin_log=grids_f%lin_log(h2o_ind), &
         & dxdt_tan=dxdt_tan, d2xdxdt_tan=d2xdxdt_tan )
    else
      call allocate_test ( dxdt_tan,ptan%template%nosurfs,sv_t_len, &
                        &  'dxdt_tan',moduleName )
      call allocate_test ( d2xdxdt_tan,ptan%template%nosurfs,sv_t_len, &
                         & 'd2xdxdt_tan',moduleName )
      call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad, &
         & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp, &
         & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
         & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
         & req_out, tan_chi_out, dhdz_out, dx_dh_out, &
         & dxdt_tan=dxdt_tan, d2xdxdt_tan=d2xdxdt_tan )
    end if

! Compute Gauss Legendre (GL) grid ---------------------------------------

! Insert automatic preselected integration gridder here. Need to make a
! large concatenated vector of bases and pointings.

    call allocate_test ( z_all,temp%template%nosurfs+2, 'z_all', moduleName )

! the -3.000 is a designated "surface" value
! initialize step

    z_all = (/-3.000_rp,RESHAPE(temp%template%surfs(:,1), &
         &  (/temp%template%nosurfs/)),4.0_rp/)
! see if pointing grid is associated, if so concatenate it to the
! state vector
    if ( associated(FwdModelConf%tangentGrid) ) then
      call allocate_test ( z_tmp,SIZE(z_ALL) + FwdModelConf%tangentGrid%nosurfs, &
      & 'z_tmp',moduleName )
      z_tmp = (/z_all,FwdModelConf%tangentGrid%surfs/)

! Move z_tmp to z_all

      call deallocate_test ( z_all, 'z_all',moduleName )
      call allocate_test ( z_all,SIZE(z_tmp), 'z_all',moduleName )
      z_all = z_tmp
      call deallocate_test ( z_tmp, 'z_tmp',moduleName )
    end if

    do sps_i = 1 , no_mol

      f => GetQuantityForForwardModel(fwdmodelin,fwdmodelextra, &
          &  quantitytype = l_vmr, molIndex=mol_cat_index(sps_i), config=fwdModelConf, &
          &  radiometer = firstsignal%radiometer )

! Concatenate vector

      j = SIZE(z_all) + f%template%nosurfs
      call allocate_test ( z_tmp, j, 'z_tmp', moduleName )
      z_tmp = (/z_all,f%template%surfs(:,1)/)

! Move z_tmp to z_all

      call allocate_test ( z_all,SIZE(z_tmp), 'z_all',moduleName )
      z_all = z_tmp
      call deallocate_test ( z_tmp, 'z_tmp',moduleName )

    end do

! On top of all of these, add the original Integration Grid:

    j = SIZE(z_all) + Size(FwdModelConf%integrationGrid%surfs)
    call allocate_test ( z_tmp, j, 'z_tmp', moduleName )
    z_tmp = (/z_all,FwdModelConf%integrationGrid%surfs/)

! Move z_tmp to z_all

    call allocate_test ( z_all, SIZE(z_tmp), 'z_all', moduleName )
    z_all = z_tmp
    call deallocate_test ( z_tmp, 'z_tmp', moduleName )

! Now, create the final grid:

    call make_z_grid ( z_all, z_psig, rec_tan_inds )
    call deallocate_test ( z_all, 'z_all', moduleName )

! note that z_psig(1) is the designated surface
    Nlvl = SIZE(z_psig)
    NLm1 = Nlvl - 1
    maxVert = NLm1 * Ngp1 + 1

! Allocate GL grid stuff

    call allocate_test ( z_glGrid, maxVert, 'z_glGrid', moduleName )
    call allocate_test ( p_glGrid, maxVert, 'p_glGrid', moduleName )

    call allocate_test ( h_glgrid, maxVert, no_sv_p_t, 'h_glgrid', moduleName )
    call allocate_test ( t_glgrid, maxVert, no_sv_p_t, 't_glgrid', moduleName )
    call allocate_test ( dhdz_glgrid, maxVert, no_sv_p_t, 'dhdz_glgrid', &
                      &  moduleName )
    call allocate_test ( dh_dt_glgrid, maxVert,  n_t_zeta, no_sv_p_t, &
                      &  'dh_dt_glgrid', moduleName )
    call allocate_test ( ddhidhidtl0, maxVert, n_t_zeta, no_sv_p_t, &
                      &  'ddhidhidtl0', moduleName )

! From the selected integration grid pressures define the GL pressure grid:

    z_glgrid(1:maxVert-1) = reshape ( &
      ! Midpoint of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) + z_psig(1:Nlm1)),1,Ngp1) + &
      ! Half length of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) - z_psig(1:Nlm1)),1,Ngp1) * &
      ! Gauss points (with -1 at front):
      & spread((/-1.0_rp,Gx(1:Ng)/),2,NLm1), (/maxVert-1/))
    z_glgrid(maxVert) = z_psig(Nlvl)
    p_glgrid = 10.0_rp**(-z_glgrid)

    call deallocate_test ( z_psig, 'z_psig', moduleName )

    ! Compute hydrostatic grid -----------------------------------------------

    ! Insert into bill's 2d hydrostatic equation.
    ! The phi input for this program are the orbit plane projected
    ! geodetic locations of the temperature phi basis--not necessarily
    ! the tangent phi's which may be somewhat different.

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_Begin ( 'ForwardModel.Hydrostatic' )

    call two_d_hydrostatic ( Grids_tmp, &
      &  (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
      &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
      &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
      &  dhdz_glgrid, dh_dt_glgrid, DDHDHDTL0=ddhidhidtl0 )

! Compute tan_press from fwdModelConf%tangentGrid%surfs

    j = COUNT(fwdModelConf%tangentGrid%surfs < (z_glgrid(1) - 0.0001_rp))
    no_tan_hts = Nlvl + j
    call allocate_test ( tan_inds, no_tan_hts, 'tan_inds', moduleName )
    call allocate_test ( tan_press, no_tan_hts, 'tan_press', moduleName )
    call allocate_test ( tan_phi, no_tan_hts, 'tan_phi', moduleName )
    call allocate_test ( est_scgeocalt, no_tan_hts, 'est_scgeocalt', moduleName )
    call allocate_test ( tan_hts, no_tan_hts, 'tan_hts', moduleName )
    call allocate_test ( tan_temps, no_tan_hts, 'tan_temps', moduleName )
    call allocate_test ( reqs, no_tan_hts, 'reqs', moduleName )
    tan_inds(1:j) = 1
    tan_inds(j+1:no_tan_hts) = (rec_tan_inds - 1) * Ngp1 + 1
    call deallocate_test ( rec_tan_inds, 'rec_tan_inds', moduleName )

    tan_press(1:j) = fwdModelConf%tangentGrid%surfs(1:j)
    tan_press(j+1:no_tan_hts) = z_glgrid(tan_inds(j+1:no_tan_hts))

    surfaceTangentIndex = COUNT(tan_inds == 1)

! estimate tan_phi

    tan_phi(1:j) = phitan%values(1,MAF)
    est_scgeocalt(1:j) = scGeocAlt%values(1,maf)

! Since the interpolateValues routine needs the OldX array to be sorted
! we have to sort ptan%values and re-arrange phitan%values & scgeocalt%values

    k = ptan%template%noSurfs
    call allocate_test ( z_path, k, 'h_path', moduleName )
    call allocate_test ( p_path, k, 'p_path', moduleName )
    call allocate_test ( t_path, k, 't_path', moduleName )

    z_path(1:k) = ptan%values(1:k,maf)
    p_path(1:k) = phitan%values(1:k,maf)
    t_path(1:k) = scgeocalt%values(1:k,maf)

    do i = 2, k
      m = -1
      do jf = k, i, -1
        if ( z_path(jf) < z_path(jf-1) ) then
          m = jf - 1
          r = z_path(jf)
          z_path(jf) = z_path(m)
          z_path(m) = r
          r = p_path(jf)
          p_path(jf) = p_path(m)
          p_path(m) = r
          r = t_path(jf)
          t_path(jf) = t_path(m)
          t_path(m) = r
        end if
      end do
      if ( m < 1) exit
    end do

    call interpolateValues ( z_path, p_path, tan_press(j+1:no_tan_hts), &
      &  tan_phi(j+1:no_tan_hts), METHOD = 'L' )
    call interpolateValues ( z_path, t_path, tan_press(j+1:no_tan_hts), &
       & est_scgeocalt(j+1:no_tan_hts), METHOD='L' )

    tan_phi(1:no_tan_hts) = tan_phi(1:no_tan_hts) * deg2rad

    call deallocate_test ( z_path, 'h_path', moduleName )
    call deallocate_test ( p_path, 'p_path', moduleName )
    call deallocate_test ( t_path, 't_path', moduleName )

! This is a lazy way to get the surface angle

    if ( temp_der ) then
      call allocate_test ( dxdt_surface,1,sv_t_len, 'dxdt_surface',moduleName )
      call allocate_test ( d2xdxdt_surface,1,sv_t_len, 'd2xdxdt_surface', &
                       & moduleName )
      if ( h2o_ind > 0 ) then
        end_ind_z = SUM(grids_f%no_z(1:h2o_ind))
        beg_ind_z = end_ind_z - grids_f%no_z(h2o_ind) + 1
        end_ind_p = SUM(grids_f%no_p(1:h2o_ind))
        beg_ind_p = end_ind_p - grids_f%no_p(h2o_ind) + 1
        end_ind = SUM(grids_f%no_z(1:h2o_ind)*grids_f%no_p(1:h2o_ind) * &
                & grids_f%no_f(1:h2o_ind))
        beg_ind = end_ind - grids_f%no_z(h2o_ind)*grids_f%no_p(h2o_ind) * &
                & grids_f%no_f(h2o_ind) + 1
        call get_chi_out ( tan_press(1:1), tan_phi(1:1), &
           & 0.001_rp*est_scgeocalt(1:1), Grids_tmp, &
           & SPREAD(refGPH%template%surfs(1,1),1,no_sv_p_t), &
           & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
           & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
           & (/req_out(1)/), surf_angle, one_dhdz, one_dxdh, &
           & h2o_zeta_basis=grids_f%zet_basis(beg_ind_z:end_ind_z), &
           & h2o_phi_basis=grids_f%phi_basis(beg_ind_p:end_ind_p), &
           & h2o_coeffs=grids_f%values(beg_ind:end_ind), &
           & lin_log=grids_f%lin_log(h2o_ind), &
           & dxdt_tan=dxdt_surface, d2xdxdt_tan=d2xdxdt_surface )
      else if ( h2o_ind == 0 ) then
        call get_chi_out ( tan_press(1:1), tan_phi(1:1), &
           & 0.001_rp*est_scgeocalt(1:1), Grids_tmp, &
           & SPREAD(refGPH%template%surfs(1,1),1,no_sv_p_t), &
           & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
           & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
           & (/req_out(1)/), surf_angle, one_dhdz, one_dxdh,  &
           & dxdt_tan=dxdt_surface, d2xdxdt_tan=d2xdxdt_surface )
      end if
      call deallocate_test ( d2xdxdt_surface, 'd2xdxdt_surface',moduleName )
    else
! no temperature derivatives
      if ( h2o_ind > 0 ) then
        end_ind_z = SUM(grids_f%no_z(1:h2o_ind))
        beg_ind_z = end_ind_z - grids_f%no_z(h2o_ind) + 1
        end_ind_p = SUM(grids_f%no_p(1:h2o_ind))
        beg_ind_p = end_ind_p - grids_f%no_p(h2o_ind) + 1
        end_ind = SUM(grids_f%no_z(1:h2o_ind)*grids_f%no_p(1:h2o_ind) * &
                & grids_f%no_f(1:h2o_ind))
        beg_ind = end_ind - grids_f%no_z(h2o_ind)*grids_f%no_p(h2o_ind) * &
                & grids_f%no_f(h2o_ind) + 1
        call get_chi_out ( tan_press(1:1), tan_phi(1:1), &
           & 0.001_rp*est_scgeocalt(1:1), Grids_tmp, &
           & SPREAD(refGPH%template%surfs(1,1),1,no_sv_p_t), &
           & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
           & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad, &
           & (/req_out(1)/), surf_angle, one_dhdz, one_dxdh, &
           & h2o_zeta_basis=grids_f%zet_basis(beg_ind_z:end_ind_z), &
           & h2o_phi_basis=grids_f%phi_basis(beg_ind_p:end_ind_p), &
           & h2o_coeffs=grids_f%values(beg_ind:end_ind), &
           & lin_log=grids_f%lin_log(h2o_ind) )
      else if ( h2o_ind == 0 ) then
        call get_chi_out ( tan_press(1:1), tan_phi(1:1), &
           & 0.001_rp*est_scgeocalt(1:1), Grids_tmp, &
           & SPREAD(refGPH%template%surfs(1,1),1,no_sv_p_t), &
           & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
           & orbIncline%values(1,maf)*Deg2Rad, elevoffset%values(1,1)*Deg2Rad,&
           & (/req_out(1)/), surf_angle, one_dhdz, one_dxdh )
      end if
      call deallocate_test ( d2xdxdt_surface, 'd2xdxdt_surface', moduleName )
    end if
    call deallocate_test ( req_out, 'req_out', moduleName )

 ! Now, allocate other variables we're going to need later ----------------

    allocate ( k_spect_dw(noUsedChannels, no_tan_hts, maxNoFFreqs, &
      & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dw' )
    allocate ( k_spect_dn(noUsedChannels, no_tan_hts, maxNoFFreqs, &
      & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dn' )
    allocate ( k_spect_dv(noUsedChannels, no_tan_hts, maxNoFFreqs, &
      & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_spect_dv' )
    call allocate_test ( tan_temp, no_tan_hts, 'tan_temp', moduleName )

    ! Allocate path quantities -----------------------------------------------

    ! First, Allocate gl_slab arrays....

    no_ele = 2*maxVert     ! maximum possible

    allocate ( gl_slabs ( no_ele, noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"gl_slabs" )

    do i = 1, noSpecies
      nl = size(My_Catalog(i)%Lines)
      gl_slabs(1:no_ele,i)%no_lines = nl
      do j = 1, no_ele
        call AllocateOneSlabs ( gl_slabs(j,i), nl )
      end do
    end do

    if ( temp_der ) then
      allocate ( gl_slabs_p(no_ele,noSpecies), &
        &  gl_slabs_m(no_ele,noSpecies), STAT=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//"gl_slabs_[pm]" )

      do i = 1, noSpecies
        nl = size(My_Catalog(i)%Lines)
        gl_slabs_m(1:no_ele,i)%no_lines = nl
        gl_slabs_p(1:no_ele,i)%no_lines = nl
        do j = 1, no_ele
          call AllocateOneSlabs ( gl_slabs_p(j,i), nl )
          call AllocateOneSlabs ( gl_slabs_m(j,i), nl )
        end do
      end do
    end if

    ! Now allocate all path related... with maximum length..

    brkpt = maxVert
    npc = 2 * (brkpt + Ng) / Ngp1

    ! This can be put outside the mmaf loop

    call allocate_test ( dhdz_path, no_ele, 'dhdz_path', moduleName )
    call allocate_test ( h_path,    no_ele, 'h_path',    moduleName )
    call allocate_test ( p_path,    no_ele, 'p_path',    moduleName )
    call allocate_test ( path_dsdh, no_ele, 'path_dsdh', moduleName )
    call allocate_test ( phi_path,  no_ele, 'phi_path',  moduleName )
    call allocate_test ( t_path,    no_ele, 't_path',    moduleName )
    call allocate_test ( z_path,    no_ele, 'z_path',    moduleName )

    call allocate_test ( alpha_path_c, npc, 'alpha_path_c', moduleName )
    call allocate_test ( incoptdepth,  npc, 'incoptdept',   moduleName )
    call allocate_test ( t_script,     npc, 't_script',     moduleName )
    call allocate_test ( tau,          npc, 'tau',          moduleName )
    call allocate_test ( del_s,        npc, 'del_s',        moduleName )
    call allocate_test ( do_gl,        npc, 'do_gl',        moduleName )
    call allocate_test ( indices_c,    npc, 'indices_c',    moduleName )
    call allocate_test ( n_path,       npc, 'n_path',       moduleName )
    call allocate_test ( ref_corr,     npc, 'ref_corr',     moduleName )

    call allocate_test ( beta_path_c, npc, no_mol, 'beta_path_c', moduleName )
    call allocate_test ( do_calc_fzp, no_ele, f_len, 'do_calc_fzp', moduleName )
    call allocate_test ( do_calc_zp, no_ele, p_len, 'do_calc_zp', moduleName )
    call allocate_test ( eta_zp, no_ele, p_len, 'eta_zp', moduleName )
    call allocate_test ( eta_fzp, no_ele, f_len, 'eta_fzp', moduleName )
    call allocate_test ( sps_path, no_ele, no_mol, 'sps_path', moduleName )

    if ( temp_der ) then

! Allocation for metrics routine when Temp. derivative is needed:

!      call allocate_test ( k_temp, noUsedChannels, no_tan_hts, n_t_zeta, &
!                         & no_sv_p_t, 'k_temp',moduleName )
      allocate ( k_temp(noUsedChannels, no_tan_hts, n_t_zeta, no_sv_p_t) )

      call allocate_test ( dRad_dt, sv_t_len, 'dRad_dt', moduleName )
      call allocate_test ( dbeta_dt_path_c,npc,no_mol, 'dbeta_dt_path_c', &
                         & moduleName )
      call allocate_test ( dh_dt_path, no_ele, sv_t_len, 'dh_dt_path', &
                         & moduleName )
      call allocate_test ( do_calc_hyd, no_ele, sv_t_len, 'do_calc_hyd', &
                         & moduleName )
      call allocate_test ( do_calc_t, no_ele, sv_t_len, 'do_calc_t', &
                         & moduleName )
      call allocate_test ( eta_zxp_t, no_ele, sv_t_len, 'eta_zxp_t', &
                         & moduleName )
      call allocate_test ( tan_dh_dt,1,sv_t_len, 'tan_dh_dt', moduleName )
      call allocate_test ( tan_d2h_dhdt,1,sv_t_len, 'tan_d2h_dhdt', moduleName )

    end if

    if ( atmos_der ) then
      call allocate_test ( dRad_df, f_len, 'dRad_df', moduleName )
      call allocate_test ( k_atmos,noUsedChannels,no_tan_hts,f_len, 'k_atmos',&
                       & moduleName )
      k_atmos = 0.0
    end if

    if ( spect_der ) then

      ! Allocation when spectral derivative are needed:

      call allocate_test ( dbeta_dw_path_c, npc, no_mol, &
        & 'dbeta_dw_path_c', moduleName )
      call allocate_test ( dbeta_dn_path_c, npc, no_mol, &
        & 'dbeta_dn_path_c', moduleName )
      call allocate_test ( dbeta_dv_path_c, npc, no_mol, &
        & 'dbeta_dv_path_c', moduleName )

      f_len_dw = SUM(Grids_dw%no_z(:) * Grids_dw%no_p(:) * Grids_dw%no_f(:))
      f_len_dn = SUM(Grids_dn%no_z(:) * Grids_dn%no_p(:) * Grids_dn%no_f(:))
      f_len_dv = SUM(Grids_dv%no_z(:) * Grids_dv%no_p(:) * Grids_dv%no_f(:))

      call allocate_test ( do_calc_dw, no_ele, f_len_dw, 'do_calc_dw', &
                        &  moduleName )
      call allocate_test ( do_calc_dn, no_ele, f_len_dn, 'do_calc_dn', &
                        &  moduleName )
      call allocate_test ( do_calc_dv, no_ele, f_len_dv, 'do_calc_dv', &
                        &  moduleName )

      call allocate_test ( eta_zxp_dw, no_ele, f_len_dw, 'eta_zxp_dw', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dn, no_ele, f_len_dn, 'eta_zxp_dn', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dv, no_ele, f_len_dv, 'eta_zxp_dv', &
                        &  moduleName )

      call allocate_test ( drad_dw, f_len_dw, 'drad_dw', moduleName )
      call allocate_test ( drad_dn, f_len_dn, 'drad_dn', moduleName )
      call allocate_test ( drad_dv, f_len_dv, 'drad_dv', moduleName )

    end if

    call allocate_test ( ptg_angles,no_tan_hts, 'ptg_angles',moduleName )
    call allocate_test ( radiances, no_tan_hts, noUsedChannels, &
                       & 'radiances', moduleName )

    call allocate_test ( dx_dt, no_tan_hts,sv_t_len, 'dx_dt',moduleName )
    call allocate_test ( d2x_dxdt,no_tan_hts,sv_t_len, 'd2x_dxdt',moduleName )

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
          & 'superset', moduleName )
        do i = 1, size(pointingGrids)
          superset(i) = AreSignalsSuperset ( pointingGrids(i)%signals, &
            & fwdModelConf%signals, sideband=thisSideband )
        end do
        if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error,ModuleName, &
               & "No matching pointing frequency grids." )
        maxSuperset = maxval ( superset )
        where ( superset < 0 ) superset = maxSuperset + 1
        whichPointingGrid = minloc ( superset, 1 )
        call deallocate_test ( superset, 'superset', moduleName )

        ! Now we've identified the pointing grids.  Locate the tangent grid
        ! within it.
        call allocate_test ( grids, no_tan_hts, "Grids", moduleName )
        call Hunt ( PointingGrids(whichPointingGrid)%oneGrid%height, &
                 &  tan_press, grids, allowTopValue=.TRUE., nearest=.TRUE. )
        ! Work out the maximum number of frequencies
        maxNoPtgFreqs = 0
        do ptg_i = 1, no_tan_hts
          k = Size(pointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))% &
                 & frequencies)
          maxNoPtgFreqs = max ( maxNoPtgFreqs, k )
        end do

        min_ch_freq_grid =  huge(min_ch_freq_grid)
        max_ch_freq_grid = -huge(min_ch_freq_grid)
        do i = 1, noUsedChannels
          sigInd = usedSignals(i)
          channel = usedChannels(i)
          shapeInd = MatchSignal ( filterShapes%signal, &
            & fwdModelConf%signals(sigInd), sideband = thisSideband, &
            & channel=channel )
          if ( shapeInd == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                               & "No matching channel shape information" )
          k = Size(FilterShapes(shapeInd)%FilterGrid(:))
          r1 = FilterShapes(shapeInd)%FilterGrid(1)
          r2 = FilterShapes(shapeInd)%FilterGrid(k)
          min_ch_freq_grid = MIN(r1, r2, min_ch_freq_grid)
          max_ch_freq_grid = MAX(r1, r2, max_ch_freq_grid)
        end do

      else ! ------------------------- Not frequency averaging ---------

        call allocate_test ( frequencies, noUsedChannels, "frequencies", &
                           & moduleName )
        do channel = 1, noUsedchannels
          direction = fwdModelConf%signals(usedSignals(channel))%direction
          frequencies(channel) = &
            & fwdModelConf%signals(usedSignals(channel))%centerFrequency + &
            & direction * fwdModelConf%signals(usedSignals(channel))% &
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
        maxNoPtgFreqs = noUsedChannels
      end if

      call allocate_test ( RadV, maxNoPtgFreqs, 'RadV', moduleName )

      if ( temp_der) &
        & call allocate_test ( k_temp_frq, maxNoPtgFreqs, sv_t_len, 'k_temp_frq', &
                             &  moduleName )

      if ( atmos_der ) &
        & call allocate_test ( k_atmos_frq, maxNoPtgFreqs, f_len, 'k_atmos_frq',&
                             & moduleName )

      if ( spect_der ) then
        call allocate_test ( k_spect_dw_frq , maxNoPtgFreqs, f_len_dw , &
                           & 'k_spect_dw_frq', moduleName )
        call allocate_test ( k_spect_dn_frq , maxNoPtgFreqs, f_len_dn , &
                           & 'k_spect_dn_frq', moduleName )
        call allocate_test ( k_spect_dv_frq , maxNoPtgFreqs, f_len_dv , &
                           & 'k_spect_dv_frq', moduleName )
      end if

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

        ! This is not pretty but we need some coarse grid extraction indices
        k = Ngp1
        j = (npc+1)/2
        indices_c(1:npc) = (/(i*k-Ng,i=1,j),((i-1)*k-Ng+1,i=j+1,npc)/)
        indices_c(npc+1:) = 0

        ! Compute z_path & p_path
        z_path(1:no_ele) = (/(z_glgrid(i),i=MaxVert,tan_inds(ptg_i),-1), &
                           & (z_glgrid(i),i=tan_inds(ptg_i),MaxVert)/)
        p_path(1:no_ele) = (/(p_glgrid(i),i=MaxVert,tan_inds(ptg_i),-1), &
                           & (p_glgrid(i),i=tan_inds(ptg_i),MaxVert)/)

        ! Compute the h_path, t_path, dhdz_path, phi_path, dhdt_path

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_Begin ( 'ForwardModel.MetricsEtc' )

        if ( ptg_i < surfaceTangentIndex ) then
          neg_tan_ht = temp%values(1,(windowstart+windowfinish)/2) * &
            &  (tan_press(ptg_i) - z_glgrid(1)) / 14.8_rp
          e_rflty = earthRefl%values(1,1)

          if ( temp_der ) then
            ! Set up temperature representation basis stuff
            call metrics ( tan_phi(ptg_i:ptg_i), tan_inds(ptg_i:ptg_i),      &
              &  Grids_tmp%phi_basis, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid, &
              &  orbIncline%values(1,maf)*Deg2Rad, Grids_tmp%deriv_flags,    &
              &  h_path(1:no_ele), phi_path(1:no_ele),                       &
              &  t_path(1:no_ele), dhdz_path(1:no_ele), Req,                 &
              &  TAN_PHI_H_GRID=one_tan_ht, TAN_PHI_T_GRID=one_tan_temp,     &
              &  NEG_H_TAN = (/neg_tan_ht/), DHTDTL0 = tan_dh_dt,            &
              &  DDHIDHIDTL0 = ddhidhidtl0, DDHTDHTDTL0 = tan_d2h_dhdt,      &
              &  DHIDTLM = dh_dt_glgrid,                                     &
              &  DHITDTLM = dh_dt_path(1:no_ele,:),                          &
              &  Z_BASIS = Grids_tmp%zet_basis,                              &
              &  ETA_ZXP=eta_zxp_t(1:no_ele,:),                              &
              &  DO_CALC_T = do_calc_t(1:no_ele,:),                          &
              &  DO_CALC_HYD = do_calc_hyd(1:no_ele,:) )
          else
            call metrics ( tan_phi(ptg_i:ptg_i), tan_inds(ptg_i:ptg_i),      &
              &  Grids_tmp%phi_basis, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid, &
              &  orbIncline%values(1,maf)*Deg2Rad,                           &
              &  dummy, h_path(1:no_ele), phi_path(1:no_ele),                &
              &  t_path(1:no_ele), dhdz_path(1:no_ele), Req,                 &
              &  TAN_PHI_H_GRID = one_tan_ht, TAN_PHI_T_GRID = one_tan_temp, &
              &  NEG_H_TAN = (/neg_tan_ht/) )
          end if

! Tan heights for negative tan height from metrics is not correctly working.

        else
          e_rflty = 1.0_rp
          if ( temp_der ) then
            ! Set up temperature representation basis stuff
            call metrics ( tan_phi(ptg_i:ptg_i), tan_inds(ptg_i:ptg_i),      &
              &  Grids_tmp%phi_basis, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid, &
              &  orbIncline%values(1,maf)*Deg2Rad, Grids_tmp%deriv_flags,    &
              &  h_path(1:no_ele), phi_path(1:no_ele),                       &
              &  t_path(1:no_ele), dhdz_path(1:no_ele), Req,                 &
              &  TAN_PHI_H_GRID = one_tan_ht, TAN_PHI_T_GRID = one_tan_temp, &
              &  DHTDTL0 = tan_dh_dt, DDHIDHIDTL0 = ddhidhidtl0,             &
              &  DDHTDHTDTL0 = tan_d2h_dhdt, DHIDTLM = dh_dt_glgrid,         &
              &  DHITDTLM = dh_dt_path(1:no_ele,:),                          &
              &  Z_BASIS = Grids_tmp%zet_basis,                              &
              &  ETA_ZXP = eta_zxp_t(1:no_ele,:),                            &
              &  DO_CALC_T = do_calc_t(1:no_ele,:),                          &
              &  DO_CALC_HYD = do_calc_hyd(1:no_ele,:) )
          else
            call metrics ( tan_phi(ptg_i:ptg_i), tan_inds(ptg_i:ptg_i),      &
              &  Grids_tmp%phi_basis, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid, &
              &  orbIncline%values(1,maf)*Deg2Rad,                           &
              &  dummy, h_path(1:no_ele), phi_path(1:no_ele),                &
              &  t_path(1:no_ele), dhdz_path(1:no_ele), Req,                 &
              &  TAN_PHI_H_GRID = one_tan_ht, TAN_PHI_T_GRID = one_tan_temp )
          end if
        end if
        ! Fill the diagnostic arrays
        tan_temps ( ptg_i ) = one_tan_temp ( 1 )
        tan_hts ( ptg_i ) = one_tan_ht ( 1 )
        reqs ( ptg_i ) = req
        !  ** Determine the eta_zxp_dw, eta_zxp_dn, eta_zxp_dv
        if ( spect_der ) then
          call eval_spect_path ( Grids_dw, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dw(1:no_ele,:), eta_zxp_dw(1:no_ele,:) )
          call eval_spect_path ( Grids_dn, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dn(1:no_ele,:), eta_zxp_dn(1:no_ele,:) )
          call eval_spect_path ( Grids_dv, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dv(1:no_ele,:), eta_zxp_dv(1:no_ele,:) )
        end if

        tan_temp(ptg_i) = one_tan_temp(1)

        ! Now compute the eta_zp & do_calc_zp (for Zeta & Phi only)
        call comp_eta_docalc_no_frq ( Grids_f,z_path(1:no_ele), &
          &  phi_path(1:no_ele), do_calc_zp(1:no_ele,:), eta_zp(1:no_ele,:) )

       ! Now compute sps_path with a FAKE frequency, mainly to get the
       ! WATER (H2O) contribution for refraction calculations, but also
       ! to compute sps_path for all those witn no frequency component

        Frq = 0.0
        call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
          & Frq,eta_zp(1:no_ele,:), &
          & do_calc_zp(1:no_ele,:), sps_path(1:no_ele,:),      &
          & do_calc_fzp(1:no_ele,:), eta_fzp(1:no_ele,:) )

        if ( h2o_ind > 0 ) then
          call refractive_index ( p_path(indices_c(1:npc)), &
            &  t_path(indices_c(1:npc)), n_path(1:npc),     &
            &  h2o_path=sps_path((indices_c(1:npc)), h2o_ind) )
        else
          call refractive_index ( p_path(indices_c(1:npc)), &
            &  t_path(indices_c(1:npc)), n_path(1:npc) )
        end if

        if ( temp_der ) then
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path(npc/2),&
             & one_tan_ht(1), tan_phi(ptg_i), Req, 0.0_rp,                &
             & ptg_angles(ptg_i), r, 1.0_rp, tan_dh_dt(1,:),              &
             & tan_d2h_dhdt(1,:), dx_dt(ptg_i,:), d2x_dxdt(ptg_i,:) )
        else
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path(npc/2),&
             & one_tan_ht(1), tan_phi(ptg_i),Req, 0.0_rp,                 &
             & ptg_angles(ptg_i), r, 1.0_rp )
        end if

        call comp_refcor ( Req+h_path(indices_c(1:npc)), 1.0_rp+n_path(1:npc), &
                      &  Req+one_tan_ht(1), del_s(1:npc), ref_corr(1:npc) )

! This only needs to be computed on the gl (not coarse) grid thus there is
! some duplication here.
        path_dsdh(2:brkpt-1) = path_ds_dh(Req+h_path(2:brkpt-1), &
              & Req+one_tan_ht(1))
        path_dsdh(brkpt+2:no_ele-1)=path_ds_dh(Req+h_path(brkpt+2:no_ele-1), &
              & Req+one_tan_ht(1))

        ! Compute ALL the slabs_prep entities over the path's GL grid for this
        ! pointing & mmaf:

        del_temp = 0.0_rp
        call get_gl_slabs_arrays ( my_Catalog, p_path(1:no_ele), &
          &  t_path(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs, &
          &  no_ele, del_temp )

        if ( temp_der ) then
          del_temp = 10.0_rp
          call get_gl_slabs_arrays ( my_Catalog, p_path(1:no_ele), &
            &  t_path(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs_p, &
            &  no_ele, del_temp )
          call get_gl_slabs_arrays ( my_Catalog, p_path(1:no_ele), &
            &  t_path(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs_m, &
            &  no_ele, -del_temp )
        end if

        ! Work out what frequencies we're using for --------------------------
        ! frequency averaging case for this pointing

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.

        if ( FwdModelConf%do_freq_avg ) then
          j = -1
          k = SIZE(PointingGrids(whichPointingGrid)%oneGrid( &
                       & grids(ptg_i))%frequencies)
          call Hunt ( min_ch_freq_grid, PointingGrids(whichPointingGrid)% &
                       & oneGrid(grids(ptg_i))%frequencies, k, j, frq_i )
          call Hunt ( max_ch_freq_grid, PointingGrids(whichPointingGrid)% &
                       & oneGrid(grids(ptg_i))%frequencies, k, frq_i, m )
          noFreqs = m - j + 1
          call allocate_test ( frequencies,noFreqs, "frequencies", moduleName )
          frequencies(1:noFreqs) = PointingGrids(whichPointingGrid)%&
                                   &oneGrid(grids(ptg_i))%frequencies(j:m)

! VELOCITY shift correction to frequency grid

          frequencies =  Vel_Cor * frequencies

        end if

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
          call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
            & Frq,eta_zp(1:no_ele,:), &
            & do_calc_zp(1:no_ele,:), sps_path(1:no_ele,:),      &
            & do_calc_fzp(1:no_ele,:), eta_fzp(1:no_ele,:) )

          if ( temp_der  .and. spect_der ) then

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),  &
              &  my_Catalog, beta_group, gl_slabs, indices_c(1:npc),       &
              &  beta_path_c(1:npc,:), GL_SLABS_M=gl_slabs_m,              &
              &  T_PATH_M=t_path(1:no_ele)-del_temp, GL_SLABS_P=gl_slabs_p,&
              &  T_PATH_P=t_path(1:no_ele)+del_temp,                       &
              &  DBETA_DT_PATH=dbeta_dt_path_c(1:npc,:),                   &
              &  DBETA_DW_PATH=dbeta_dw_path_c(1:npc,:),                   &
              &  DBETA_DN_PATH=dbeta_dn_path_c(1:npc,:),                   &
              &  DBETA_DV_PATH=dbeta_dv_path_c(1:npc,:) )

          else if ( temp_der ) then

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),  &
              &  my_Catalog, beta_group, gl_slabs, indices_c(1:npc),       &
              &  beta_path_c(1:npc,:),                                     &
              &  GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp,&
              &  GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp,&
              &  DBETA_DT_PATH=dbeta_dt_path_c(1:npc,:) )

          else if ( spect_der ) then

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),    &
              &  my_Catalog, beta_group, gl_slabs, indices_c(1:npc),         &
              &  beta_path_c(1:npc,:), DBETA_DW_PATH=dbeta_dw_path_c(1:npc,:),&
              &  DBETA_DN_PATH=dbeta_dn_path_c(1:npc,:),                     &
              &  DBETA_DV_PATH=dbeta_dv_path_c(1:npc,:) )

          else

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),  &
              &  my_Catalog, beta_group, gl_slabs, indices_c(1:npc),       &
              &  beta_path_c(1:npc,:) )

          end if

          alpha_path_c(1:npc) = SUM(sps_path(indices_c(1:npc),:) *  &
                                  & beta_path_c(1:npc,:),DIM=2)

          call path_contrib ( alpha_path_c(1:npc), del_s(1:npc), e_rflty,   &
                  & fwdModelConf%tolerance, tau(1:npc), incoptdepth(1:npc), &
                  & do_gl(1:npc) )

          ! ALLOCATE gl grid beta

          no_gl_ndx = count(do_gl(1:npc))
          j = Ng * no_gl_ndx

          call allocate_test ( gl_inds, j, 'gl_inds', moduleName )
          call allocate_test ( beta_path_f, j, no_mol, 'beta_path_f', &
                             & moduleName )
          call allocate_test ( gl_indgen, Ng, no_gl_ndx, 'gl_indgen', &
                             & moduleName )
          call allocate_test ( gl_ndx, no_gl_ndx, 2, 'gl_ndx', moduleName )

          gl_ndx(:,1) = pack((/(i,i=1,npc)/),do_gl(1:npc))

  ! Make (/(j-Ng-1,j=1,Ng)/), (/(j,j=1,Ng)/) parameter variables later on

          do i = 1 , no_gl_ndx
            if ( gl_ndx(i,1) > npc/2 ) then
              gl_ndx(i,2) = 1 - Ng + Ngp1 * (gl_ndx(i,1) - 1)
              gl_indgen(:,i) = (/(j,j=1,Ng)/)
            else
              gl_ndx(i,2) = 1 +      Ngp1 * (gl_ndx(i,1) - 1)
              gl_indgen(:,i) = (/(j-Ng-1,j=1,Ng)/)
            end if
          end do

          ! compute the gl indicies

          gl_inds = reshape(spread(gl_ndx(1:no_gl_ndx,2),1,Ng) +  &
            &  gl_indgen,(/Ng*no_gl_ndx/))

          j = Ng * no_gl_ndx
          if ( temp_der  .and. spect_der ) then

            call allocate_test ( dbeta_dt_path_f, j, no_mol, &
              & 'dbeta_dt_path_f', moduleName )
            call allocate_test ( dbeta_dw_path_f, j, no_mol, &
              & 'dbeta_dw_path_f', moduleName )
            call allocate_test ( dbeta_dn_path_f, j, no_mol, &
              & 'dbeta_dn_path_f', moduleName )
            call allocate_test ( dbeta_dv_path_f, j, no_mol, &
              & 'dbeta_dv_path_f', moduleName )

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),   &
              & my_Catalog, beta_group, gl_slabs, gl_inds, beta_path_f,     &
              & GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp,  &
              & GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp,  &
              & DBETA_DT_PATH=dbeta_dt_path_f, DBETA_DW_PATH=dbeta_dw_path_f,&
              & DBETA_DN_PATH=dbeta_dn_path_f, DBETA_DV_PATH=dbeta_dv_path_f )

          else if ( temp_der ) then

            call allocate_test ( dbeta_dt_path_f, j, no_mol, &
              & 'dbeta_dt_path_f', moduleName )
            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),    &
              &   my_Catalog, beta_group, gl_slabs, gl_inds, beta_path_f,    &
              &   GL_SLABS_M=gl_slabs_m, T_PATH_M=t_path(1:no_ele)-del_temp, &
              &   GL_SLABS_P=gl_slabs_p, T_PATH_P=t_path(1:no_ele)+del_temp, &
              &   DBETA_DT_PATH=dbeta_dt_path_f )

          else if ( spect_der ) then

            call allocate_test ( dbeta_dw_path_f, j, no_mol, &
                              & 'dbeta_dw_path_f', moduleName )
            call allocate_test ( dbeta_dn_path_f, j, no_mol, &
                              & 'dbeta_dn_path_f', moduleName )
            call allocate_test ( dbeta_dv_path_f, j, no_mol, &
                              & 'dbeta_dv_path_f', moduleName )
            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele),    &
              & my_Catalog, beta_group, gl_slabs, gl_inds, beta_path_f,      &
              & DBETA_DW_PATH=dbeta_dw_path_f, DBETA_DN_PATH=dbeta_dn_path_f,&
              & DBETA_DV_PATH=dbeta_dv_path_f )

          else

            call get_beta_path ( Frq, p_path(1:no_ele), t_path(1:no_ele), &
              &  my_Catalog, beta_group, gl_slabs, gl_inds, beta_path_f )

          end if

          ! Compute radiative transfer ---------------------------------------

          call rad_tran ( Frq, spaceRadiance%values(1,1), e_rflty,           &
            & z_path(indices_c(1:npc)), t_path(indices_c(1:npc)),            &
            & alpha_path_c(1:npc), ref_corr(1:npc), do_gl(1:npc),            &
            & incoptdepth(1:npc), SUM(sps_path(gl_inds,:)*beta_path_f,DIM=2),&
            & path_dsdh(gl_inds), dhdz_path(gl_inds), t_script(1:npc),       &
            & tau(1:npc),Rad,i_stop )

          RadV(frq_i) = Rad

          ! Compute derivatives if needed ----------------------------------

          if ( atmos_der ) then

            call drad_tran_df ( z_path(indices_c(1:npc)), Grids_f,             &
              &  beta_path_c(1:npc,:), eta_fzp(indices_c(1:npc),:),            &
              &  sps_path(indices_c(1:npc),:), do_calc_fzp(indices_c(1:npc),:),&
              &  beta_path_f, eta_fzp(gl_inds,:), sps_path(gl_inds,:),    &
              &  do_calc_fzp(gl_inds,:), do_gl(1:npc), del_s(1:npc),      &
              &  ref_corr(1:npc), path_dsdh(gl_inds), dhdz_path(gl_inds), &
              &  t_script(1:npc), tau(1:npc), i_stop, drad_df, ptg_i, frq_i )

            k_atmos_frq(frq_i,1:f_len) = drad_df(1:f_len)

          end if

          if ( temp_der ) then

            call drad_tran_dt ( z_path(indices_c(1:npc)),                      &
              & Req+h_path(indices_c(1:npc)),                                  &
              & t_path(indices_c(1:npc)), dh_dt_path(indices_c(1:npc),:),      &
              & alpha_path_c(1:npc), SUM(sps_path(indices_c(1:npc),:) *        &
              & dbeta_dt_path_c(1:npc,:) * beta_path_c(1:npc,:),DIM=2),        &
              & eta_zxp_t(indices_c(1:npc),:), do_calc_t(indices_c(1:npc),:),  &
              & do_calc_hyd(indices_c(1:npc),:), del_s(1:npc), ref_corr(1:npc),&
              & Req + one_tan_ht(1), dh_dt_path(brkpt,:), frq, do_gl(1:npc),   &
              & req + h_path(gl_inds), t_path(gl_inds), dh_dt_path(gl_inds,:), &
              & SUM(sps_path(gl_inds,:)*beta_path_f,DIM=2),                    &
              & SUM(sps_path(gl_inds,:)*beta_path_f*dbeta_dt_path_f,DIM=2),    &
              & eta_zxp_t(gl_inds,:), do_calc_t(gl_inds,:), path_dsdh(gl_inds),&
              & dhdz_path(gl_inds), t_script(1:npc), tau(1:npc), i_stop,       &
              & drad_dt, ptg_i, frq_i )

            k_temp_frq(frq_i,:) = drad_dt

          end if

          if ( spect_der ) then

            ! Spectroscopic derivative  wrt: W

            call drad_tran_dx ( z_path(indices_c(1:npc)), Grids_dw,           &
              &  dbeta_dw_path_c(1:npc,:), eta_zxp_dw(indices_c(1:npc),:),    &
              &  sps_path(indices_c(1:npc),:), do_calc_dw(indices_c(1:npc),:),&
              &  dbeta_dw_path_f, eta_zxp_dw(gl_inds,:), sps_path(gl_inds,:), &
              &  do_calc_dw(gl_inds,:), do_gl(1:npc), del_s(1:npc),           &
              &  ref_corr(1:npc), path_dsdh(gl_inds), dhdz_path(gl_inds),     &
              &  t_script(1:npc), tau(1:npc), i_stop, drad_dw, ptg_i, frq_i )

            k_spect_dw_frq(frq_i,1:1:f_len_dw) = drad_dw(1:1:f_len_dw)

            ! Spectroscopic derivative  wrt: N

            call drad_tran_dx ( z_path(indices_c(1:npc)), Grids_dn,           &
              &  dbeta_dn_path_c(1:npc,:), eta_zxp_dn(indices_c(1:npc),:),    &
              &  sps_path(indices_c(1:npc),:), do_calc_dn(indices_c(1:npc),:),&
              &  dbeta_dn_path_f, eta_zxp_dn(gl_inds,:), sps_path(gl_inds,:), &
              &  do_calc_dn(gl_inds,:), do_gl(1:npc), del_s(1:npc),           &
              &  ref_corr(1:npc), path_dsdh(gl_inds), dhdz_path(gl_inds),     &
              &  t_script(1:npc), tau(1:npc), i_stop, drad_dn, ptg_i, frq_i )

            k_spect_dn_frq(frq_i,1:f_len_dn) = drad_dn(1:f_len_dn)

            ! Spectroscopic derivative  wrt: Nu0

            call drad_tran_dx ( z_path(indices_c(1:npc)), Grids_dv,           &
              &  dbeta_dv_path_c(1:npc,:), eta_zxp_dv(indices_c(1:npc),:),    &
              &  sps_path(indices_c(1:npc),:), do_calc_dv(indices_c(1:npc),:),&
              &  dbeta_dv_path_f, eta_zxp_dv(gl_inds,:), sps_path(gl_inds,:), &
              &  do_calc_dv(gl_inds,:), do_gl(1:npc), del_s(1:npc),           &
              &  ref_corr(1:npc), path_dsdh(gl_inds), dhdz_path(gl_inds),     &
              &  t_script(1:npc), tau(1:npc), i_stop, drad_dv, ptg_i, frq_i )

            k_spect_dv_frq(frq_i,1:1:f_len_dv) = drad_dv(1:1:f_len_dv)

          end if

          call deallocate_test ( gl_inds, 'gl_inds', moduleName )
          call deallocate_test ( beta_path_f, 'beta_path_f', moduleName )
          call deallocate_test ( gl_ndx, 'gl_ndx', moduleName )
          call deallocate_test ( gl_indgen, 'gl_indgen', moduleName )

          if ( temp_der ) &
            & call deallocate_test ( dbeta_dt_path_f, 'dbeta_dt_path_f', &
            & moduleName )

          if ( spect_der ) then
            call deallocate_test ( dbeta_dw_path_f, 'dbeta_dw_path_f',moduleName )
            call deallocate_test ( dbeta_dn_path_f, 'dbeta_dn_path_f',moduleName )
            call deallocate_test ( dbeta_dv_path_f, 'dbeta_dv_path_f',moduleName )
          end if

          ! End of frequency loop ----------------------------------------------

          if ( toggle(emit) .and. levels(emit) > 5 ) &
            & call Trace_End ('ForwardModel.Frequency ',index=frq_i)

        end do            ! End freq. loop

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_End ( 'ForwardModel.FrequencyLoop' )

        ! Work out which channel shape information we're going to need -------

        ! Frequency averaging if needed --------------------------------------

        ! Here we either frequency average to get the unconvolved radiances, or
        ! we just store what we have as we're using monochromatic channels

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call trace_begin ( 'ForwardModel.FrequencyAvg' )

        if ( fwdModelConf%do_freq_avg ) then
          do i = 1, noUsedChannels
            sigInd = usedSignals(i)
            channel = usedChannels(i)
            shapeInd = MatchSignal ( filterShapes%signal, &
              & fwdModelConf%signals(sigInd), sideband = thisSideband, &
              & channel=channel )
            if ( shapeInd == 0 ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              &    "No matching channel shape information" )
            k = size(FilterShapes(shapeInd)%FilterGrid)
            call Freq_Avg ( frequencies, &
              &   FilterShapes(shapeInd)%FilterGrid, &
              &   FilterShapes(shapeInd)%FilterShape, &
              &   RadV, noFreqs, k, Radiances(ptg_i,i) )
          end do
        else
          Radiances(ptg_i,1:noUsedChannels) = RadV(1:)
        end if

        ! Frequency averaging of derivatives if needed -----------------------

        ! Frequency Average the temperature derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( temp_der ) then
          if ( fwdModelConf%do_freq_avg ) then
            do i = 1, noUsedChannels
              sigInd = usedSignals(i)
              channel = usedChannels(i)
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), &
                & sideband = thisSideband, channel=channel )
              sv_i = 1
              do instance = 1, no_sv_p_t
                do surface = 1, n_t_zeta
                  call Freq_Avg ( frequencies, &
                    & FilterShapes(shapeInd)%FilterGrid, &
                    & FilterShapes(shapeInd)%FilterShape, &
                    & k_temp_frq(:,sv_i), noFreqs, &
                    & size(FilterShapes(shapeInd)%FilterGrid), r )
                  k_temp(i,ptg_i,surface,instance) = r
                  sv_i = sv_i + 1
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! Channel loop
          else
            do i = 1, noUsedChannels
              sv_i = 1
              do instance = 1, no_sv_p_t
                do surface = 1, n_t_zeta
                  k_temp(i,ptg_i,surface,instance) = k_temp_frq(i,sv_i)
                  sv_i = sv_i + 1
                end do
              end do
            end do
          end if
        end if

        ! Frequency Average the atmospheric derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( atmos_der ) then

          sv_i = 1
          do k = 1, no_mol
            specie = mol_cat_index(k)
            ! Did have a test here for moleculeDerivatives(specie), however
            ! this didn't work as it screwed up sv_1 which is supposed to accumulate
            ! over all molecules.  However, deriv_flags conveys the same information
            ! so not much speed will be lost.  NJL.
            if ( fwdModelConf%do_freq_avg ) then
              sv_start = sv_i
              do i = 1, noUsedChannels
                sv_i = sv_start
                sigInd = usedSignals(i)
                channel = usedChannels(i)
                j = Size(FilterShapes(shapeInd)%FilterGrid)
                shapeInd = MatchSignal ( filterShapes%signal, &
                  & fwdModelConf%signals(sigInd), &
                  & sideband = thisSideband, channel=channel )
                do instance = Grids_f%WindowStart(k), Grids_f%WindowFinish(k)
                  do surface = 1, Grids_f%no_f(k)*Grids_f%no_z(k)
                    if ( grids_f%deriv_flags(sv_i) ) then
                      call Freq_Avg ( frequencies, &
                        & FilterShapes(shapeInd)%FilterGrid, &
                        & FilterShapes(shapeInd)%FilterShape, &
                        & k_atmos_frq(1:noFreqs,sv_i), noFreqs, j, r )
                    else
                      r = 0.0
                    end if
                    k_atmos(i,ptg_i,sv_i) = r
                    sv_i = sv_i + 1
                  end do                ! Surface loop
                end do                  ! Instance loop
              end do                    ! Channel loop
            else                        ! Else not frequency averaging
              sv_start = sv_i
              do i = 1, noUsedChannels
                sv_i = sv_start
                do instance = Grids_f%WindowStart(k), Grids_f%WindowFinish(k)
                  do surface = 1, Grids_f%no_f(k)*Grids_f%no_z(k)
                    if ( grids_f%deriv_flags(sv_i) ) then
                      k_atmos(i,ptg_i,sv_i) = k_atmos_frq(i,sv_i)
                    else
                      k_atmos(i,ptg_i,sv_i) = 0.0
                    end if
                    sv_i = sv_i + 1
                  end do
                end do
              end do
            end if                      ! Frequency averaging or not
          end do                          ! Loop over major molecules
          !
        end if                          ! Want derivatives for atmos
        
        ! Frequency Average the spectroscopic derivatives with the appropriate
        ! filter shapes
        !??? Do we need to do this if there's no Jacobian ???

        if ( spect_der ) then

          !  *** dI/dW

          do k = 1, no_mol
            specie = mol_cat_index(k)
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid)
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dw%no_z(k)
                      do jf = 1, Grids_dw%no_f(k)
                        call Freq_Avg ( frequencies, &
                          & FilterShapes(shapeInd)%FilterGrid, &
                          & FilterShapes(shapeInd)%FilterShape,&
                          & k_spect_dw_frq(:,sv_i), noFreqs, j, r )
                        k_spect_dw(i,ptg_i,jf,surface,instance,k) = r
                        sv_i = sv_i + 1
                      end do              ! Frequencies loop
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! else not frequency averaging
                do i = 1, noUsedChannels
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dw%no_z(k)
                      do jf = 1, Grids_dw%no_f(k)
                        k_spect_dw(i,ptg_i,jf,surface,instance,k) = &
                          k_spect_dw_frq(i,sv_i)
                        sv_i = sv_i + 1
                      end do
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over major molecules

          !  *** dI/dN

          do k = 1, no_mol
            specie = mol_cat_index(k)
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid)
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dn%no_z(k)
                      do jf = 1, Grids_dn%no_f(k)
                        call Freq_Avg ( frequencies, &
                          & FilterShapes(shapeInd)%FilterGrid, &
                          & FilterShapes(shapeInd)%FilterShape,&
                          & k_spect_dn_frq(:,sv_i), noFreqs, j, r )
                        k_spect_dn(i,ptg_i,jf,surface,instance,k) = r
                        sv_i = sv_i + 1
                      end do              ! Frequencies loop
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! else not frequency averaging
                do i = 1, noUsedChannels
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dn%no_z(k)
                      do jf = 1, Grids_dn%no_f(k)
                        k_spect_dn(i,ptg_i,jf,surface,instance,k) = &
                          k_spect_dn_frq(i,sv_i)
                        sv_i = sv_i + 1
                      end do
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over major molecules

          !  *** dI/dV

          do k = 1, no_mol
            specie = mol_cat_index(k)
            if ( fwdModelConf%moleculeDerivatives(specie) ) then
              if ( fwdModelConf%do_freq_avg ) then
                do i = 1, noUsedChannels
                  sigInd = usedSignals(i)
                  channel = usedChannels(i)
                  j = Size(FilterShapes(shapeInd)%FilterGrid)
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd), &
                    & sideband = thisSideband, channel=channel )
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dv%no_z(k)
                      do jf = 1, Grids_dv%no_f(k)
                        call Freq_Avg ( frequencies, &
                          & FilterShapes(shapeInd)%FilterGrid, &
                          & FilterShapes(shapeInd)%FilterShape,&
                          & k_spect_dv_frq(:,sv_i), noFreqs, j, r )
                        k_spect_dv(i,ptg_i,jf,surface,instance,k) = r
                        sv_i = sv_i + 1
                      end do              ! Frequencies loop
                    end do                ! Surface loop
                  end do                  ! Instance loop
                end do                    ! Channel loop
              else                        ! else not frequency averaging
                do i = 1, noUsedChannels
                  if ( k == 1) sv_i = 1
                  do instance = WindowStart, WindowFinish
                    do surface = 1, Grids_dv%no_z(k)
                      do jf = 1, Grids_dv%no_f(k)
                        k_spect_dv(i,ptg_i,jf,surface,instance,k) = &
                          k_spect_dv_frq(i,sv_i)
                        sv_i = sv_i + 1
                      end do
                    end do
                  end do
                end do
              end if                      ! Frequency averaging or not
            end if                        ! Want derivatives for this
          end do                          ! Loop over major molecules

        end if                        ! Want derivatives for spect

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call trace_end ( 'ForwardModel.FrequencyAvg' )

        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call trace_end ( 'ForwardModel.Pointing ', index=ptg_i )

        if ( FwdModelConf%do_freq_avg ) &
           & call deallocate_test ( frequencies, 'frequencies',moduleName )

        ! End of pointing loop -------------------------------------------------
      end do

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_End ( 'ForwardModel.PointingLoop' )

!       ! EXTRA DEBUG FOR NATHANIEL/BILL ********************
!       call dump ( tan_temps, 'tan_temps' )
!       call dump ( tan_press, 'tan_press' )
!       call dump ( ptg_angles, 'ptg_angles' )
!       call dump ( tan_hts, 'tan_hts' )
!       call dump ( reqs, 'reqs' )
!       call dump ( est_scgeocalt, 'est_scgeocalt' )
!       call dump ( grids_f%values, 'grids_f' )

      ! Convolution if needed, or interpolation to ptan ----------------------

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_begin ( 'ForwardModel.Convolution' )

      ! Work out which antenna patterns we're going to need ------------------

      call allocate_test ( superset, size(antennaPatterns), &
        & 'superset', moduleName )

      do i = 1, noUsedChannels

        channel = usedChannels(i)
        chanInd = channel + 1 - channelOrigins(i)
        sigInd = usedSignals(i)
        thisRadiance =>  &
          GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=fwdModelConf%signals(sigInd)%sideband )
        if ( sidebandStart /= sidebandStop ) then   ! We're folding
          sidebandRatio => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_sidebandRatio, signal=fwdModelConf%signals(sigInd)%index, &
            & sideband=thisSideband, config=fwdModelConf )
          thisRatio = sidebandRatio%values(chanInd,1)
        else                  ! Otherwise, want just unfolded signal
          thisRatio = 1.0
        end if

        ! Here comes the Convolution codes
        update = ( thisSideband /= sidebandStart )

        if ( FwdModelConf%do_conv ) then

          do j = 1, size(antennaPatterns)
            superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
              & fwdModelConf%signals( (/sigInd/) ), sideband=thisSideband, &
              & channel=channel )
          end do

          if ( all( superset < 0 ) ) &
            call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching antenna patterns." )

          maxSuperset = maxval ( superset )
          where ( superset < 0 ) superset = maxSuperset + 1
          whichPattern = minloc ( superset, 1 )

          ! Now change channel from starting at 0 or 1 to definately 1
          
          j = sv_t_len
          if ( .not. temp_der .AND. .not. atmos_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, mol_cat_index, temp,ptan, &
               & thisRadiance, update, ptg_angles, Radiances(:,i), tan_chi_out,        &
               & dhdz_out, dx_dh_out, thisRatio,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & PTAN_DER=ptan_der )
          else if ( temp_der .AND. .not. atmos_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra,maf,   &
               & chanInd, windowStart, windowFinish, mol_cat_index, temp,ptan, &
               & thisRadiance, update, ptg_angles, Radiances(:,i), tan_chi_out,        &
               & dhdz_out, dx_dh_out, thisRatio,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),        &
               & DX_DT=dx_dt, D2X_DXDT=d2x_dxdt, DXDT_TAN=dxdt_tan,            &
               & DXDT_SURFACE=dxdt_surface, PTAN_DER=ptan_der )
          else if ( atmos_der .AND. .not. temp_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, mol_cat_index, temp, ptan,&
               & thisRadiance, update, ptg_angles, Radiances(:,i), tan_chi_out,        &
               & dhdz_out, dx_dh_out, thisRatio,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!              & DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,f_len/))) )
          else
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, mol_cat_index, temp, ptan,&
               & thisRadiance, update, ptg_angles, Radiances(:,i), tan_chi_out,        &
               & dhdz_out, dx_dh_out, thisRatio,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),        &
               & DX_DT=dx_dt, D2X_DXDT=d2x_dxdt, DXDT_TAN=dxdt_tan,            &
               & DXDT_SURFACE=dxdt_surface,                                    &
               & DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!              & DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,f_len/))) )
          end if

        else          ! No convolution needed ..

          j = sv_t_len
          if ( .not. temp_der .AND. .not. atmos_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out,          &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisRatio, mol_cat_index, fmStat%rows,   &
              &  Jacobian, PTAN_DER=ptan_der)
          else if ( temp_der .AND. .not. atmos_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out,          &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisRatio, mol_cat_index, fmStat%rows,   &
              &  Jacobian,                                                &
              &  DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),   &
              &  PTAN_DER=ptan_der )
          else if ( atmos_der .AND. .not. temp_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out,          &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisRatio, mol_cat_index, fmStat%rows,   &
              &  Jacobian, DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!             &  DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,f_len/))) )
          else
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out,          &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisRatio, mol_cat_index, fmStat%rows,   &
              &  Jacobian,                                                &
              &  DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),   &
              &  DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!             &  DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,f_len/))) )
          end if

        end if

      end do                            ! Channel loop

      call deallocate_test ( superset, 'superset', moduleName )
      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_end ( 'ForwardModel.Convolution' )

      ! Deallocate maxNoPtgFreqs stuff
      call deallocate_test ( Radv, 'RadV', moduleName )

      if ( temp_der ) &
        & call deallocate_test ( k_temp_frq, 'k_temp_frq', moduleName )

      if ( atmos_der ) &
        & call deallocate_test ( k_atmos_frq, 'k_atmos_frq', moduleName )

      if ( spect_der ) then
        call deallocate_test ( k_spect_dw_frq, 'k_spect_dw_frq', moduleName )
        call deallocate_test ( k_spect_dn_frq, 'k_spect_dn_frq', moduleName )
        call deallocate_test ( k_spect_dv_frq, 'k_spect_dv_frq', moduleName )
      end if

      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call trace_end ( 'ForwardModel.sideband ',index=thisSideband )

    end do            ! End of loop over sidebands -------------------------

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_End ( 'ForwardModel.SidebandLoop' )

    !  **** DEBUG Printing cycle ...

! *** Create *seez* file for "nasty" purposes:

    if ( index(switches, 'seez') /= 0 ) then
      include 'dump_print_code.f9h'
    end if

! *** End of include

    call GetNameOfSignal ( firstSignal, molName )
    i = Index(molName, '.B')
    do j = i+1,32
      if ( molName(j:j) == '.' ) then
        molName(j:)=' '
        exit
      end if
    end do

    if ( index(switches, 'rad') /= 0 ) then
      ! *** DEBUG Print
      if ( FwdModelConf%do_conv ) then
        print *, 'Convolution: ON'
      else
        print *, 'Convolution: OFF'
      end if

      if ( FwdModelConf%do_freq_avg ) then
        print *, 'Frequency Averaging: ON'
      else
        Frq = Frequencies(1)
        print *, 'Frequency Averaging: OFF'
        print '(a,f12.4,a)', ' (All computations done at Frq =',Frq, ')'
      end if

      k = ptan%template%noSurfs
      print "( /'ptan\ ',i3.3)", k
      Print "( 4(3x, f11.7) )", Ptan%values(1:k,maf)

      Print *
      do i = 1, noUsedChannels
        channel = usedChannels(i)
        print "(/, 'ch', i2.2, '_pfa_rad\ ', i3.3 )", channel, k
        j = thisRadiance%template%noChans
        print "( 4(2x, 1pg15.8) )", &
          & thisRadiance%values(channel:channel+j*(k-1):j, maf)
      end do
      Print *

    end if

    !  **** End of Printing cycle ...

    ! Now deallocate lots of stuff
    do i = 1, size(my_catalog)
      if ( associated ( my_catalog(i)%lines ) ) &
        & call deallocate_test ( my_catalog(i)%lines, 'my_catalog(?)%lines', &
        & moduleName )
    end do
    deallocate ( my_catalog, stat=ier )
    ! Note that we don't deallocate the signals/sidebands stuff for each line
    ! as these are shallow copies of the main spectroscopy catalog stuff

    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'my_catalog' )

    do i = 1, no_mol
      deallocate ( beta_group(i)%cat_index, stat=j )
      deallocate ( beta_group(i)%ratio, stat=j )
    end do

    deallocate ( mol_cat_index, stat=j )
    deallocate ( beta_group, stat=j )

    call deallocate_test ( usedChannels, 'usedChannels', moduleName )
    call deallocate_test ( channelOrigins, 'channelOrigins', moduleName )
    call deallocate_test ( usedSignals, 'usedSignals', moduleName )

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

    call deallocate_test ( z_glGrid, 'z_glGrid', moduleName )
    call deallocate_test ( p_glGrid, 'z_glGrid', moduleName )

    call deallocate_test ( h_glgrid, 'h_glgrid', moduleName )
    call deallocate_test ( t_glgrid, 't_glgrid', moduleName )
    call deallocate_test ( dhdz_glgrid, 'dhdz_glgrid', moduleName )
    call deallocate_test ( dh_dt_glgrid, 'dh_dt_glgrid', moduleName )
    call deallocate_test ( ddhidhidtl0, 'ddhidhidtl0', moduleName )

    call deallocate_test ( tan_inds, 'tan_inds', moduleName )
    call deallocate_test ( tan_press, 'tan_press', moduleName )
    call deallocate_test ( tan_phi, 'tan_phi', moduleName )
    call deallocate_test ( tan_hts, 'tan_hts', moduleName )
    call deallocate_test ( tan_temps, 'tan_temps', moduleName )
    call deallocate_test ( reqs, 'reqs', moduleName )
    call deallocate_test ( est_scgeocalt, 'est_scgeocalt', moduleName )
    call deallocate_test ( tan_temp, 'tan_temp', moduleName )

    call DestroyCompleteSlabs ( gl_slabs )
    call destroygrids_t ( grids_f )
    call destroygrids_t ( grids_tmp )

    call deallocate_test ( dhdz_path, 'dhdz_path', moduleName )
    call deallocate_test ( h_path,    'h_path',    moduleName )
    call deallocate_test ( p_path,    'p_path',    moduleName )
    call deallocate_test ( path_dsdh, 'path_dsdh', moduleName )
    call deallocate_test ( phi_path,  'phi_path',  moduleName )
    call deallocate_test ( t_path,    't_path',    moduleName )
    call deallocate_test ( z_path,    'z_path',    moduleName )

    call deallocate_test ( alpha_path_c, 'alpha_path_c', moduleName )
    call deallocate_test ( incoptdepth,  'incoptdept',   moduleName )
    call deallocate_test ( t_script,     't_script',     moduleName )
    call deallocate_test ( tau,          'tau',          moduleName )
    call deallocate_test ( del_s,        'del_s',        moduleName )
    call deallocate_test ( do_gl,        'do_gl',        moduleName )
    call deallocate_test ( indices_c,   'indices_c',   moduleName )
    call deallocate_test ( n_path,       'n_path',       moduleName )
    call deallocate_test ( ref_corr,     'ref_corr',     moduleName )

    call deallocate_test ( beta_path_c, 'beta_path_c', moduleName )
    call deallocate_test ( do_calc_zp, 'do_calc_zp', moduleName )
    call deallocate_test ( do_calc_fzp, 'do_calc_fzp', moduleName )
    call deallocate_test ( eta_zp, 'eta_zp', moduleName )
    call deallocate_test ( eta_fzp, 'eta_fzp', moduleName )
    call deallocate_test ( sps_path, 'sps_path', moduleName )

    call deallocate_test ( tan_chi_out, 'tan_chi_out',moduleName )
    call deallocate_test ( dx_dh_out, 'dx_dh_out',moduleName )
    call deallocate_test ( dhdz_out, 'dhdz_out',moduleName )

    if ( temp_der ) then
      deallocate ( k_temp, STAT=i )
      call deallocate_test ( dRad_dt, 'dRad_dt', moduleName )
      call deallocate_test ( dbeta_dt_path_c, 'dbeta_dt_path_c', moduleName )
      call deallocate_test ( dh_dt_path, 'dh_dt_path', moduleName )
      call deallocate_test ( do_calc_hyd, 'do_calc_hyd', moduleName )
      call deallocate_test ( do_calc_t, 'do_calc_t', moduleName )
      call deallocate_test ( eta_zxp_t, 'eta_zxp_t', moduleName )
      call deallocate_test ( tan_dh_dt, 'tan_dh_dt', moduleName )
      call deallocate_test ( tan_d2h_dhdt, 'tan_d2h_dhdt', moduleName )
      call deallocate_test ( dxdt_tan, 'dxdt_tan',moduleName )
      call deallocate_test ( d2xdxdt_tan, 'd2xdxdt_tan',moduleName )
      call deallocate_test ( dxdt_surface, 'dxdt_surface',moduleName )
      call DestroyCompleteSlabs ( gl_slabs_p )
      call DestroyCompleteSlabs ( gl_slabs_m )
    end if

    if ( atmos_der ) then
      call deallocate_test ( k_atmos, 'k_atmos', moduleName )
      call deallocate_test ( dRad_df, 'dRad_df', moduleName )
    end if

    if ( spect_der ) then

      call deallocate_test ( dbeta_dw_path_c, 'dbeta_dw_path_c', moduleName )
      call deallocate_test ( dbeta_dn_path_c, 'dbeta_dn_path_c', moduleName )
      call deallocate_test ( dbeta_dv_path_c, 'dbeta_dv_path_c', moduleName )

      call deallocate_test ( do_calc_dw, 'do_calc_dw', moduleName )
      call deallocate_test ( do_calc_dn, 'do_calc_dn', moduleName )
      call deallocate_test ( do_calc_dv, 'do_calc_dv', moduleName )

      call deallocate_test ( eta_zxp_dw, 'eta_zxp_dw', moduleName )
      call deallocate_test ( eta_zxp_dn, 'eta_zxp_dn', moduleName )
      call deallocate_test ( eta_zxp_dv, 'eta_zxp_dv', moduleName )

      call deallocate_test ( drad_dw, 'drad_dw', moduleName )
      call deallocate_test ( drad_dn, 'drad_dn', moduleName )
      call deallocate_test ( drad_dv, 'drad_dv', moduleName )

    end if

    call deallocate_test ( ptg_angles, 'ptg_angles', moduleName )
    call deallocate_test ( radiances, 'radiances', moduleName )
    call deallocate_test ( dx_dt, 'dx_dt', moduleName )
    call deallocate_test ( d2x_dxdt, 'd2x_dxdt', moduleName )

    ! Deallocate all variables allocated earlier -----------------------------

    if ( toggle(emit) ) then
      call trace_end ( 'ForwardModel MAF=',fmStat%maf )
    end if

  end subroutine FullForwardModel

! --------------------------------------------------  not_used_here  -----
  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module FullForwardModel_m

! $Log$
! Revision 2.98  2002/10/26 00:13:35  livesey
! Made the warning about lines less common as it checks for continuum
! aswell.
!
! Revision 2.97  2002/10/25 01:12:43  livesey
! Added an array to accumulate the stuff such as one_tan_ht etc.
! Also put in but commented out some code to dump it.
!
! Revision 2.96  2002/10/10 19:38:22  vsnyder
! Mostly cosmetic, some performance improvements
!
! Revision 2.95  2002/10/10 01:46:50  livesey
! Whoops, typo fix
!
! Revision 2.94  2002/10/10 01:28:06  livesey
! Bug fix in k_temp windowing
!
! Revision 2.93  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.92  2002/10/04 23:49:50  vsnyder
! More cosmetic changes
!
! Revision 2.91  2002/10/02 23:06:42  vsnyder
! Add 'seez' switch, instead of uncommenting include, to get Zvi's debug print.
! Get SpeedOfLight from Geometry.  Cosmetic changes.
!
! Revision 2.90  2002/09/26 18:02:49  livesey
! Now uses GetQuantityForForwardModel.
!
! Revision 2.89  2002/09/10 17:05:38  livesey
! New update arguments to convolve/noconvole
!
! Revision 2.88  2002/09/07 02:18:50  vsnyder
! More cosmetic changes
!
! Revision 2.87  2002/09/06 20:23:22  livesey
! Merged in bug fixes with change from Van
!
! Revision 2.86  2002/09/06 18:19:04  vsnyder
! Cosmetic changes
!
! Revision 2.85  2002/09/05 21:00:38  vsnyder
! Get rid of some auxiliary variables
!
! Revision 2.84  2002/08/26 20:02:22  livesey
! Checks whether the jacobian is present when setting internal
! temp/atmos_der
!
! Revision 2.83  2002/08/22 23:13:03  livesey
! New intermediate frequency based frq_bases
!
! Revision 2.82  2002/08/20 23:33:23  livesey
! Fixed bug with handling of extinction
!
! Revision 2.81  2002/08/20 22:37:04  livesey
! Moved uses inside routine
!
! Revision 2.80  2002/08/02 00:12:07  bill
! just testing
!
! Revision 2.79  2002/07/31 23:27:54  bill
! added feature to user supply grid points to the psig and tangent grids
!
! Revision 2.78  2002/07/31 00:03:45  livesey
! Embarassing bug fix, was seeking wrong vector quantity.
!
! Revision 2.77  2002/07/30 20:03:59  livesey
! More sideband fixes
!
! Revision 2.76  2002/07/30 20:03:19  livesey
! Fixed bug which had it using the wrong sideband ratio
!
! Revision 2.75  2002/07/29 21:41:30  bill
! no changes, just debugging
!
! Revision 2.74  2002/07/23 22:26:31  livesey
! Added ptan_der, set k_atmos to zero on creation
!
! Revision 2.73  2002/07/19 23:35:49  bill
! fixed undefined surf angle
!
! Revision 2.72  2002/07/11 20:51:20  bill
! fixed bug regarding req
!
! Revision 2.71  2002/07/09 17:37:10  livesey
! Fixed more DACs problems
!
! Revision 2.70  2002/07/08 17:45:37  zvi
! Make sure spect_der is turned off for now
!
! Revision 2.69  2002/07/05 07:52:45  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.68  2002/06/28 21:41:36  livesey
! Repeat of bug fix with atmos_der being deallocated in wrong place.
!
! Revision 2.67  2002/06/28 11:06:46  zvi
! Now computing dI/dPtan using chain rule
!
! Revision 2.66  2002/06/26 19:58:48  livesey
! Bug fix with DAC channel numbering
!
! Revision 2.65  2002/06/24 21:11:24  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.61  2002/06/17 17:12:15  bill
! fixed yet another bug--wgr
!
! Revision 2.60  2002/06/17 16:30:52  bill
! inc zvis changes--wgr
!
! Revision 2.58  2002/06/13 22:40:38  bill
! some variable name changes--wgr
!
! Revision 2.57  2002/06/11 22:20:10  bill
! work in progress--wgr
!
! Revision 2.56  2002/06/07 23:22:58  bill
! add debug switch--wgr
!
! Revision 2.55  2002/06/07 23:01:25  bill
! work in progress
!
! Revision 2.54  2002/06/07 04:50:03  bill
! fixes and improvements--wgr
!
! Revision 2.53  2002/06/05 17:20:28  livesey
! Fixed tan_temp
!
! Revision 2.52  2002/06/04 23:06:47  livesey
! On the way to having phiTan
!
! Revision 2.51  2002/06/04 10:27:59  zvi
! Encorporate deriv. flag into convolution
!
! Revision 2.50  2002/05/28 17:09:14  livesey
! Removed print statement
!
! Revision 2.49  2002/05/24 17:10:57  livesey
! Fixed bug with my_catalog(?)%lines not being associated for parent
! species.
!
! Revision 2.48  2002/05/23 21:01:11  livesey
! No, that was the wrong thing to do.  Think a bit more.
!
! Revision 2.47  2002/05/23 20:55:20  livesey
! Put more checking around case where a molecule has no lines.
!
! Revision 2.46  2002/05/22 19:44:51  zvi
! Fix a bug in the mol. index loop
!
! Revision 2.45  2002/05/17 22:13:20  livesey
! Bug fix for case where channels start at zero.
!
! Revision 2.44  2002/05/14 22:40:45  livesey
! Bug fix in change in line gathering.  Never got to run it mercifully!
!
! Revision 2.43  2002/05/14 22:32:45  livesey
! Added single sideband stuff.  Also skip line gathering for parent
! molecules.
!
! Revision 2.42  2002/05/14 00:19:10  livesey
! Minor bug fixes
!
! Revision 2.41  2002/05/10 16:18:45  livesey
! Code for dealing with new channel shape information
!
! Revision 2.40  2002/05/08 08:53:42  zvi
! All radiometers grid concept implemented
!
! Revision 2.39  2002/05/03 23:29:18  livesey
! Added direction and split sideband ratio stuff
!
! Revision 2.38  2002/02/22 00:52:06  bill
! fixed units error for light speed--wgr
!
! Revision 2.37  2002/02/20 22:19:45  zvi
! Reversing the subset logic ..
!
! Revision 2.36  2002/02/16 10:32:16  zvi
! Fixing small bug..
!
! Revision 2.35  2002/02/16 06:49:59  zvi
! Retain deriv flag code ..
!
! Revision 2.32  2002/02/14 19:05:01  bill
! Fixed no spectral avg bug--wgr
!
! Revision 2.31  2002/02/13 20:35:47  livesey
! Added some nullifies
!
! Revision 2.30  2002/02/08 00:46:05  zvi
! Some cosmetic changes..
!
! Revision 2.29  2002/02/07 00:36:31  zvi
! Fix a bug - phi_tan non defined..
!
! Revision 2.28  2002/02/05 21:54:29  zvi
! Fix a bug ..
!
! Revision 2.27  2002/02/04 22:44:40  zvi
! Fixing some bugs in the automatic grid selection process
!
! Revision 2.26  2002/02/02 11:20:17  zvi
! Code to overwrite the l2cf integration & tanget grids
!
! Revision 2.25  2002/01/30 01:11:18  zvi
! Fix bug in user selectable coeff. code
!
! Revision 2.24  2002/01/27 08:37:45  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.22  2002/01/08 01:01:19  livesey
! Some changes to my_catalog and one_tan_height and one_tan_temp
!
! Revision 2.21  2001/12/26 04:05:04  zvi
! Convert phi_tan to Radians
!
! Revision 2.20  2001/12/14 23:43:05  zvi
! Modification for Grouping concept
!
! Revision 2.19  2001/11/25 07:57:12  zvi
! Fixing inconsistency in k_xxx instance loops
!
! Revision 2.18  2001/11/20 10:23:27  zvi
! Fixing window bug & diemsion
!
! Revision 2.17  2001/11/20 01:18:59  zvi
! Fixing Shifting Window bug
!
! Revision 2.16  2001/11/15 01:21:57  zvi
! Extiction debug fix
!
! Revision 2.15  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.14  2001/11/08 21:52:23  jonathan
! add spec_tags to Molecules
!
! Revision 2.13  2001/11/08 09:56:59  zvi
! Fixing a bug..
!
! Revision 2.12  2001/11/08 00:11:29  livesey
! Added extinction stuff
!
! Revision 2.11  2001/11/07 21:19:01  livesey
! Put Zvi's change comments back
!
! Revision 2.10  2001/11/07 21:16:56  livesey
! Now defaults to *not* using a line if no bands listed.
!
! Revision 2.9  2001/11/07 09:58:41  zvi
! More effective code for sps_path calculations
!
! Revision 2.8  2001/11/03 01:33:35  livesey
! Add more informative message if no spectroscopy information available
! for a molecule
!
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
