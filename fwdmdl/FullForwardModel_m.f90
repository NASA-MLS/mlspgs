! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module FullForwardModel_m

  ! This module contains the `full' forward model.

  implicit NONE
  private
  public :: FullForwardModel

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  ! -------------------------------------------- FullForwardModel -----

  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                             &  FwdModelOut, oldIfm, FmStat, Jacobian )

  ! This gets a little bit of data and then computes the sizes for quantities
  ! in the full forward model proper, FullForwardModelAuto.

    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use Compute_Z_PSIG_m, only: Compute_Z_PSIG
    use ForwardModelConfig, only: DeriveFromForwardModelConfig, &
      & DestroyForwardModelDerived, ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelIntermediate_t, &
                                        ForwardModelStatus_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Get_Species_Data_M, only:  Get_Species_Data
    use GLnp, only: NGP1
    use Intrinsic, only: L_CLOUDICE, L_MAGNETICFIELD, L_PHITAN, L_PTAN, &
      & L_TEMPERATURE
    use Load_Sps_Data_m, only: DestroyGrids_t, Dump, EmptyGrids_T, Grids_T, &
      & Load_One_Item_Grid, Load_Sps_Data
    use MatrixModule_1, only: MATRIX_T
    use MLSKinds, only: RP
    use Toggles, only: Switches
    use VectorsModule, only: Vector_T, VectorValue_T

    type(forwardModelConfig_T), intent(inout) :: FwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    real(rp), pointer :: Z_PSIG(:)       ! Surfs from Temperature, tangent grid
                                         ! and species grids, sans duplicates.
    real(rp), pointer :: TAN_PRESS(:)    ! Pressures corresponding to Z_PSIG
    type (Grids_T) :: Grids_tmp ! All the coordinates for TEMP
    type (Grids_T) :: Grids_f   ! All the coordinates for VMR
    type (Grids_T) :: Grids_iwc ! All the coordinates for WC
    type (Grids_T) :: Grids_mag ! All the coordinates for Magnetic field
    type (Grids_T) :: Grids_n   ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_v   ! All the spectroscopy(V) coordinates
    type (Grids_T) :: Grids_w   ! All the spectroscopy(W) coordinates
    type (VectorValue_T), pointer :: PHITAN ! Tangent geodAngle component of state vector
    type (VectorValue_T), pointer :: PTAN   ! Tangent pressure component of state vector
    type (VectorValue_T), pointer :: TEMP   ! Temperature component of state vector
    integer :: No_Mol           ! Number of molecules
    integer :: NoUsedChannels   ! Number of channels used
    integer :: No_sv_p_T        ! number of phi basis for temperature
    integer :: N_T_zeta         ! Number of zetas for temperature
    integer :: EXT_IND          ! Index of extinction inside f array
    integer :: H2O_IND          ! Index of h2o inside f array, else zero
    integer :: Sv_T_len         ! Number of t_phi*t_zeta in the window
    integer :: Nlvl             ! Number of levels in coarse zeta grid
    integer :: NO_TAN_HTS       ! Number of tangent heights
    integer :: SURFACETANGENTINDEX  ! Index in tangent grid of earth's
                                    ! surface
    integer :: MAXVERT          ! Number of points in gl-refined vertical grid:
                                ! (NLVL-1) * (NG+1) + 1, i.e., 1 + NG
                                ! per level, except the last,
                                ! where there's no GL space.
    integer :: MAX_C            ! Length of longest possible coarse path,
                                ! Z_PSIG & Min Zeta & surface Zeta
    integer :: MAX_F            ! Length of longest possible fine path
                                ! (all npf<max_f)
    integer :: S_T  ! Multiplier for temp derivative sizes, 0 or 1
    integer :: S_A  ! Multiplier for atmos derivative sizes, 0 or 1
    integer :: S_LC ! Multiplier for line center deriv sizes, 0 or 1
    integer :: S_LW ! Multiplier for line width deriv sizes, 0 or 1
    integer :: S_TD ! Multiplier for temp dependence deriv sizes, 0 or 1
    integer :: S_P  ! Multiplier for polarized sizes, 0 or 1
    integer :: S_PFA ! Multiplier for PFA sizes, 0 or 1
    integer :: S_I  ! Multiplier for ice/cloud sizes, 0 or 1

    logical :: temp_der, atmos_der, spect_der, ptan_der ! Flags for various derivatives
    logical :: Spect_Der_Center, Spect_Der_Width, Spect_Der_Width_TDep

    nullify ( z_psig, tan_press )

    call deriveFromForwardModelConfig ( fwdModelConf )

    ! Create the data structures for the species.  Get the
    ! spectroscopy parameters from the state vector.
    ! This has to be done AFTER deriveFromForwardModelConfig.

    call get_species_data ( fwdModelConf, fwdModelIn, fwdModelExtra )

    no_mol = size(fwdModelConf%beta_group)
    noUsedChannels = size(fwdModelConf%channels)
    phitan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, &
      & instrumentModule=fwdModelConf%signals(1)%instrumentModule, &
      & config=fwdModelConf )
    ptan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, &
      & instrumentModule=fwdModelConf%signals(1)%instrumentModule, &
      & foundInFirst=ptan_der, config=fwdModelConf )
    temp => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature, config=fwdModelConf )
    call load_one_item_grid ( grids_tmp, temp, phitan, fmStat%maf, fwdModelConf, .true. )
    no_sv_p_t = grids_tmp%l_p(1) ! phi == windowFinish - windowStart + 1
    n_t_zeta = grids_tmp%l_z(1)  ! zeta
    sv_t_len = grids_tmp%p_len   ! zeta X phi == n_t_zeta * no_sv_p_t

    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_f, &
      & h2o_ind, ext_ind )

    if ( FwdModelConf%polarized ) then
      call load_one_item_grid ( grids_mag, &
        & GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_magneticField, config=fwdModelConf ), phitan, &
        & fmStat%maf, fwdModelConf, .false. )
    else
      call emptyGrids_t ( grids_mag ) ! Allocate components with zero size
    end if

    if ( FwdModelConf%incl_cld ) then
      call load_one_item_grid ( grids_iwc, &
        & GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra,  &
        & quantityType=l_cloudIce, noError=.true., config=fwdModelConf ), &
        & phitan, fmStat%maf, fwdModelConf, .false., .false. )
    else
      call emptyGrids_t ( grids_iwc ) ! Allocate components with zero size
    end if

  ! Compute the preselected integration grid (all surfs from temperature,
  ! tangent grid and species).  Tan_Press is thrown in for free.

    call compute_Z_PSIG ( fwdModelConf, temp, fwdModelConf%beta_group%qty, nlvl, &
      &                   no_tan_hts, surfaceTangentIndex, z_psig, tan_press )

    max_c = 2*nlvl + 3 ! Maximum coarse path length
    maxVert = (nlvl-1) * ngp1 + 1
    max_f = 2 * maxVert + ngp1 ! Maximum fine path, including minimum Zeta panel

    temp_der = present ( jacobian ) .and. FwdModelConf%temp_der
    atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der
    ptan_der = ptan_der .and. present ( jacobian )

    spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
    spect_der_center = spect_der .and. size(fwdModelConf%lineCenter) > 0
    spect_der_width = spect_der .and. size(fwdModelConf%lineWidth) > 0
    spect_der_width_TDep = spect_der .and. size(fwdModelConf%lineWidth_TDep) > 0

    ! Allocate and fill spectroscopy derivative grids.  They'll be empty
    ! if fwdModelConf%line* has size zero.
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_v, &
      & qtyStuffIn=fwdModelConf%lineCenter%qty )
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_w, &
      & qtyStuffIn=fwdModelConf%lineWidth%qty )
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_n, &
      & qtyStuffIn=fwdModelConf%lineWidth_TDep%qty )

    if ( index(switches,'grids') /= 0 ) then ! dump the grids
      call dump ( grids_f )
      if ( size(fwdModelConf%lineCenter) > 0 ) call dump ( grids_v )
      if ( size(fwdModelConf%lineWidth) > 0 ) call dump ( grids_w )
      if ( size(fwdModelConf%lineWidth_TDep) > 0 ) call dump ( grids_n )
    end if

    s_t = merge(1,0,temp_der)
    s_a = merge(1,0,atmos_der)
    s_lc = merge(1,0,spect_der_center)
    s_lw = merge(1,0,spect_der_width)
    s_td = merge(1,0,spect_der_width_TDep)
    s_p = merge(1,0,FwdModelConf%polarized)
    s_pfa = merge(1,0,FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2))
    s_i = merge(1,0,FwdModelConf%incl_cld)

    call FullForwardModelAuto ( FwdModelConf, FwdModelIn, FwdModelExtra,     &
                              & FwdModelOut, oldIfm, FmStat, z_psig,         &
                              & tan_press, grids_tmp, grids_f, grids_mag,    &
                              & grids_iwc, grids_n, grids_v, grids_w,        &
                              & ptan, phitan, temp,                          &
                              & no_mol, noUsedChannels, no_sv_p_t, n_t_zeta, &
                              & sv_t_len, nlvl, no_tan_hts,                  &
                              & surfaceTangentIndex,                         &
                              & max_c, maxVert, max_f, EXT_ind, H2O_ind,     &
                              & ptan_der,                                    &
                              & s_t, s_a, s_lc, s_lw, s_td, s_p, s_pfa, s_i, &
                              ! Optional:
                              & Jacobian )

    call destroygrids_t ( grids_f )
    call destroygrids_t ( grids_iwc )
    call destroygrids_t ( grids_mag )
    call destroygrids_t ( grids_tmp )
    call destroygrids_t ( grids_n )
    call destroygrids_t ( grids_w )
    call destroygrids_t ( grids_v )

    call deallocate_test ( z_psig,       'z_psig',       moduleName )

    call destroyForwardModelDerived ( fwdModelConf )

  end subroutine FullForwardModel

  ! ---------------------------------------- FullForwardModelAuto  -----

  subroutine FullForwardModelAuto ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                              & FwdModelOut, oldIfm, FmStat, z_psig,         &
                              & tan_press, grids_tmp,  grids_f, grids_mag,   &
                              & grids_iwc, grids_n, grids_v, grids_w,        &
                              & ptan, phitan, temp,                          &
                              & no_mol, noUsedChannels, no_sv_p_t, n_t_zeta, &
                              & sv_t_len, nlvl, no_tan_hts,                  &
                              & surfaceTangentIndex,                         &
                              & max_c, maxVert, max_f, EXT_ind, H2O_ind,     &
                              & ptan_der,                                    &
                              & s_t, s_a, s_lc, s_lw, s_td, s_p, s_pfa, s_i, &
                              ! Optional:
                              & Jacobian )

  ! This is the full radiative transfer forward model, the workhorse code.
  ! It's called FullForwardModelAuto because most of the variable-size
  ! work arrays are automatic arrays, instead of being explicitly allocated
  ! and deallocated.

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use AntennaPatterns_m, only: ANTENNAPATTERNS
    use Comp_Eta_Docalc_No_Frq_m, only: Comp_Eta_Docalc_No_Frq
    use Comp_Sps_Path_Frq_m, only: Comp_Sps_Path, Comp_Sps_Path_Frq, &
      & Comp_Sps_Path_No_Frq, Comp_1_Sps_Path_No_Frq
    use Compute_GL_Grid_m, only: Compute_GL_Grid
    use D_Hunt_m, only: Hunt
    use D_T_SCRIPT_DTNP_M, only: DT_SCRIPT_DT
    use Dump_0, only: Dump
    use FilterShapes_m, only: DACSFilterShapes, FilterShapes
    use ForwardModelConfig, only: Beta_Group_T, Channels_T, &
      & ForwardModelConfig_t, LineCenter, LineWidth, LineWidth_TDep
    use ForwardModelIntermediate, only: ForwardModelIntermediate_t, &
                                    &   ForwardModelStatus_t, &
                                    &   B_Ptg_Angles, B_Refraction
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Freq_Avg_m, only: Freq_Avg, Freq_Avg_DACS
    use Geometry, only: Earth_Axis_Ratio_Squared_m1, EarthRadA, MaxRefraction
    use Get_Chi_Angles_m, only: Get_Chi_Angles
    use GLnp, only: GW, GX, NG, NGP1
    use Intrinsic, only: L_A, L_BOUNDARYPRESSURE, L_CLEAR, &
      & L_CLOUDWATER, L_EARTHREFL, L_ECRtoFOV, &
      & L_ELEVOFFSET, L_GPH, &
      & L_LOSVEL, L_LIMBSIDEBANDFRACTION, &
      & L_ORBITINCLINATION, L_RADIANCE, L_REFGPH, L_SCGEOCALT, &
      & L_SIZEDISTRIBUTION, L_SPACERADIANCE, L_SurfaceHeight, L_VMR
    use Load_Sps_Data_m, only: DestroyGrids_t, Dump, Grids_T, &
      & Load_One_Item_Grid
    use MatrixModule_1, only: MATRIX_T
    use Metrics_m, only: More_Metrics, Height_Metrics, Tangent_Metrics
    use Min_Zeta_m, only: Get_Min_Zeta
    use MLSKinds, only: R4, R8, RP, RV
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSNumerics, only: Hunt, InterpolateValues
    use MLSSignals_m, only: AreSignalsSuperset, GetNameOfSignal, MatchSignal, &
      & Radiometers, Signal_t
    use Molecules, only: L_H2O, L_N2O, L_O3
    use Output_m, only: Output
    use Path_Contrib_M, only: Get_GL_Inds
    use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction
    use Physics, only: H_OVER_K, SpeedOfLight
    use PointingGrid_m, only: POINTINGGRIDS
    use REFRACTION_M, only: REFRACTIVE_INDEX, COMP_REFCOR
    use SLABS_SW_M, only: ALLOCATESLABS, DESTROYCOMPLETESLABS, &
      & GET_GL_SLABS_ARRAYS, SLABS_STRUCT
! use testfield_m
    use Tau_M, only: Destroy_Tau, Dump, Get_Tau, Tau_T
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use TWO_D_HYDROSTATIC_M, only: Two_D_Hydrostatic
    use Units, only: Deg2Rad
    use VectorsModule, only: GETVECTORQUANTITYBYTYPE, VECTOR_T, VECTORVALUE_T
use Get_Eta_Matrix_m, only: Get_Eta_Sparse

    ! Extra space in the ..._R variables:  Replacement for tangent Zeta,
    ! new space for minimum Zeta, plus GL points around them.
    integer, parameter :: NXG = 2 + 4 * NG

    type(forwardModelConfig_T), intent(inout) :: FwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    real(rp), intent(in) :: Z_PSIG(:)       ! Surfs from Temperature, tangent
                                            ! grid and species grids, sans
                                            ! duplicates.
    real(rp), intent(in) :: TAN_PRESS(:)    ! Pressures corresponding to Z_PSIG
    type (Grids_T), intent(inout) :: Grids_tmp ! All the coordinates for TEMP,
                                            ! inout for cloud model to change
                                            ! supersat
    type (Grids_T), intent(inout) :: Grids_f   ! All the coordinates for VMR,
                                            ! inout for cloud model to change
                                            ! supersat
    type (Grids_T), intent(in) :: Grids_iwc ! All the coordinates for WC
    type (Grids_T), intent(in) :: Grids_mag ! All the coordinates for Magnetic field
    type (Grids_T), intent(in) :: Grids_n   ! All the spectroscopy(N) coordinates
    type (Grids_T), intent(in) :: Grids_v   ! All the spectroscopy(V) coordinates
    type (Grids_T), intent(in) :: Grids_w   ! All the spectroscopy(W) coordinates
    type (VectorValue_T), intent(in) :: PTAN ! Tangent pressure component of state vector
    type (VectorValue_T), intent(in) :: PHITAN ! Tangent geodAngle component of state vector
    type (VectorValue_T), intent(in) :: TEMP ! Temperature component of state vector
    integer, intent(in) :: No_Mol           ! Number of molecules
    integer, intent(in) :: NoUsedChannels   ! Number of channels used
    integer, intent(in) :: No_sv_p_T        ! number of phi basis for temperature
    integer, intent(in) :: N_T_zeta         ! Number of zetas for temperature
    integer, intent(in) :: Sv_T_len         ! Number of t_phi*t_zeta in the window
    integer, intent(in) :: Nlvl             ! Number of levels in coarse zeta grid
    integer, intent(in) :: NO_TAN_HTS       ! Number of tangent heights
    integer, intent(in) :: SURFACETANGENTINDEX ! Index in tangent grid of
                                            ! earth's surface
    integer, intent(in) :: MAX_C            ! Length of longest possible coarse path
    integer, intent(in) :: MAXVERT          ! Number of points in gl-refined vertical grid:
                                            ! (NLVL-1) * (NG+1) + 1, i.e., 1 + NG
                                            ! per level, except the last,
                                            ! where there's no GL space.
    integer, intent(in) :: MAX_F            ! Length of longest possible path (all npf<max_f)
    integer, intent(in) :: EXT_IND          ! Index of extinction inside f array
    integer, intent(in) :: H2O_IND          ! Index of h2o inside f array, else zero
    logical, intent(in) :: PTAN_Der
    integer, intent(in) :: S_T  ! Multiplier for temp derivative sizes, 0 or 1
    integer, intent(in) :: S_A  ! Multiplier for atmos derivative sizes, 0 or 1
    integer, intent(in) :: S_LC ! Multiplier for line center deriv sizes, 0 or 1
    integer, intent(in) :: S_LW ! Multiplier for line width deriv sizes, 0 or 1
    integer, intent(in) :: S_TD ! Multiplier for temp dependence deriv sizes, 0 or 1
    integer, intent(in) :: S_P  ! Multiplier for polarized sizes, 0 or 1
    integer, intent(in) :: S_PFA ! Multiplier for PFA sizes, 0 or 1
    integer, intent(in) :: S_I  ! Multiplier for ice/cloud sizes, 0 or 1

    type(matrix_T), intent(inout), optional :: Jacobian

    ! Now define local variables, group by type and then
    ! alphabetically

    integer :: CHANNEL            ! A Loop counter
    integer :: Frq_Avg_Sel        ! Summarizes combinations of PFA, LBL,
                                  ! Frequency averaging and derivatives.
                                  ! See Frequency_Average below.
    integer :: IER                ! Status flag from allocates
    integer :: I                  ! Loop index and other uses .
    integer :: INST               ! Relevant instance for temperature
    integer :: J                  ! Loop index and other uses ..
    integer :: K                  ! Loop index and other uses ..
    integer :: MAF                ! MAF under consideration
    integer :: MIF                ! MIF number for tan_press(ptg_i)
    integer :: Min_Index          ! If > 0, P_Path(min_index) <= Min_Phi <=
                                  ! P_Path(min_index+1), else min zeta is at
                                  ! or too close to the tangent
    integer :: Min_Index_c        ! If Min_Index > 0, P_Path(min_index_c) <=
                                  ! Min_Phi <= P_Path(min_index_c+ng+1) and
                                  ! min_index_c is a coarse path index.
    integer :: NCG                ! Number of panels needing GL = Size(cg_inds)
    integer :: NGLMAX             ! NGL if all panels need GL
    integer :: NGL                ! Total # of GL points = Size(gl_inds)
    integer :: NOFREQS            ! Number of frequencies for a pointing
    integer :: NOUSEDDACS         ! Number of different DACS in this run.
    integer :: NPC                ! Length of coarse path
    integer :: NPF                ! Length of a gl path
    integer :: NZ_IF              ! Effective size of Z_GLgrid_o and cohorts
    integer :: NZ_IG              ! Effective size of Z_IG, size(z_psig) <=
                                  ! nz_ig <= size(z_psig)+2
    integer :: PTG_I              ! Loop counter for the pointings
    integer :: SIDEBAND           ! Either zero or from firstSignal
    integer :: SIGIND             ! Signal index, loop counter
    integer :: SV_I               ! Loop index and other uses .
    integer :: SX                 ! 1 = LSB, 2 = USB
    integer :: TAN_IND_C          ! Index of tangent point in coarse zeta grid
    integer :: TAN_IND_F          ! Index of tangent point in fine grid
    integer :: TAN_PT_C           ! Index of tangent point in coarse path
    integer :: TAN_PT_F           ! Index of tangent point in fine path
    integer :: THISSIDEBAND       ! Loop counter for sidebands, -1 = LSB, +1 = USB
    integer :: WHICHPOINTINGGRID  ! Index into the pointing grids
    integer :: WINDOWFINISH       ! End of temperature `window'
    integer :: WINDOWSTART        ! Start of temperature `window'

    integer :: NoSurf             ! Number of pressure levels
    integer :: NovmrSurf          ! Number of vmr levels
    integer :: Nspec              ! No of species for cloud model
    integer :: Ispec              ! Species index in cloud model

    logical :: Any_Der            ! temp_der .or. atmos_der .or. spect_der
    logical :: cld_fine = .false.
    logical :: Clean              ! Used for dumping
    logical :: Do_Zmin            ! "Do minimum Zeta calculation"
    logical, parameter :: PFAFalse = .false.
    logical, parameter :: PFATrue = .true.
    logical :: Print_Grids        ! For debugging
    logical :: Print_Mag          ! For debugging
    logical :: Print_Min_Zeta     ! For debugging
    logical :: Print_Ptg          ! For debugging
    logical :: Print_Rad          ! For debugging
    logical :: Print_Seez         ! For debugging
    logical :: Print_TauL         ! For debugging
    logical :: Print_TauP         ! For debugging
    logical :: temp_der, atmos_der, spect_der ! Flags for various derivatives
    logical :: Spect_Der_Center, Spect_Der_Width, Spect_Der_Width_TDep

    character (len=32) :: SigName       ! Name of a Signal

    integer, target :: C_INDS_B(max_c)  ! Base array for C_INDS
    integer :: CG_INDS(max_c)     ! Indices on coarse grid where GL needed
    integer :: F_INDS(max_f)      ! Indices on fine grid
    integer, target :: GL_INDS_B(max_f) ! Base array for GL_INDS
    integer :: GRIDS(no_tan_hts)  ! Heights in ptgGrid for each tangent
    integer :: IPSD(s_i*max_f)

    logical :: DO_GL(max_c)       ! GL indicator
    logical :: T_der_path_flags(s_t*max_f) ! a flag that tells where an
      ! absorption coefficient is needed for a temperature derivative.
      ! Only useful when subsetting temperature derivatives.

    ! 'Avoid zeros' indicators
    logical :: DO_CALC_FZP(max_f, size(grids_f%values))
    logical :: DO_CALC_IWC(max_f, size(grids_iwc%values))
    logical :: DO_CALC_HYD(max_f, sv_t_len)
    logical :: DO_CALC_HYD_C(max_c, sv_t_len)  ! DO_CALC_HYD on coarse grid
    logical :: DO_CALC_N(max_f, size(grids_n%values) ) ! on entire grid
    logical :: DO_CALC_T(max_f, sv_t_len)
    logical :: DO_CALC_T_C(max_c, sv_t_len)    ! DO_CALC_T on coarse grid
    logical :: DO_CALC_T_F(max_f, sv_t_len)    ! DO_CALC_T on fine grid
    logical :: DO_CALC_V(max_f, size(grids_v%values) ) ! on entire grid
    logical :: DO_CALC_W(max_f, size(grids_w%values) ) ! on entire grid
    logical :: DO_CALC_W_C(max_c, size(grids_w%values) ) ! on coarse grid
    logical :: DO_CALC_ZP(max_f, grids_f%p_len)

    logical, pointer :: DO_CALC_Tscat(:,:)    ! 'Avoid zeros' indicator
    logical, pointer :: DO_CALC_Salb(:,:)     ! 'Avoid zeros' indicator
    logical, pointer :: DO_CALC_cext(:,:)     ! 'Avoid zeros' indicator
    logical, pointer :: DO_CALC_Tscat_ZP(:,:) ! 'Avoid zeros' indicator
    logical, pointer :: DO_CALC_Salb_ZP(:,:)  ! 'Avoid zeros' indicator
    logical, pointer :: DO_CALC_Cext_ZP(:,:)  ! 'Avoid zeros' indicator

    real(r8) :: WC(s_i*fwdModelConf%no_cloud_species, max_f)
    real(r8) :: Scat_ang(s_i*fwdModelConf%num_scattering_angles)

    real(rp) :: ALPHA_PATH_C(max_c)   ! coarse grid abs coeff.
    real(rp) :: ALPHA_PATH_F(max_f)   ! fine grid abs coeff.
    real(rp) :: B(max_c)              ! Planck radiation function
    real(rp) :: BETA_PATH_cloud_C(s_i*max_c) ! Beta on path coarse
    real(r8), target :: ChannelCenters(noUsedChannels) ! for PFA or non-frequency-averaging
    real(rp) :: DALPHA_DT_PATH_C(max_c)   ! dAlpha/dT on coarse grid
    real(rp) :: DALPHA_DT_PATH_F(max_f)   ! dAlpha/dT on fine grid
    real(rp) :: DEL_S(max_c)          ! Integration lengths along path
    real(rp) :: DEL_ZETA(max_c)       ! Integration lengths in Zeta coords
    real(rp) :: DHDZ_PATH(max_f)      ! dH/dZ on path
    real(rp) :: DHDZ_GW_PATH(max_f)   ! dH/dZ * GW on path
    real(rp) :: DSDZ_C(max_c)         ! ds/dH * dH/dZ on coarse path
    real(rp) :: DSDZ_GW_PATH(max_f)   ! ds/dH * dH/dZ * GW on path
    real(rp) :: DTanh_DT_C(max_c)     ! 1/tanh1_c d/dT tanh1_c
    real(rp) :: DTanh_DT_F(max_f)     ! 1/tanh1_f d/dT tanh1_f
    real(rp) :: H_PATH(max_f)         ! Heights on path
    real(rp) :: H_PATH_C(max_c)       ! H_PATH on coarse grid
    real(rp) :: H_PATH_F(max_f)       ! H_PATH on fine grid
    real(rp) :: INCOPTDEPTH(max_c)    ! Incremental Optical depth on coarse grid
    real(rp) :: N_PATH_C(max_c)       ! Refractive index - 1 on coarse path
    real(rp) :: N_PATH_F(max_f)       ! Refractive index - 1 on fine path
    real(rp) :: PATH_DSDH(max_f)      ! dS/dH on path
    real(rp) :: PHI_PATH(max_f)       ! Phi's on path, Radians
    real(rp), target :: P_GLGRID_O(maxvert) ! Pressure on glGrid surfs, original
    real(rp), target :: P_GLGRID_R(maxvert+nxg) ! Pressure on glGrid surfs, revised
    real(rp) :: P_PATH(max_f)         ! Pressure on path
    real(rp) :: PTG_ANGLES(no_tan_hts)
    real(rp) :: REF_CORR(max_c)       ! Refraction correction
    real(rp) :: TAN_DH_DT(s_t*sv_t_len) ! dH/dT at Tangent
    real(rp) :: TANH1_C(max_c)        ! tanh(0.5 h nu / k T)
    real(rp) :: TANH1_F(max_f)        ! tanh1 on fine grid
    real(rp) :: T_PATH(max_f)         ! Temperatures on path
    real(rp) :: T_PATH_C(max_c)       ! T_PATH on coarse grid
    real(rp) :: T_PATH_F(max_f)       ! T_PATH on fine grid
    real(rp) :: TT_PATH_C(s_i*max_c)  ! tscat on path coarse
    real(rp) :: W0_PATH_C(s_i*max_c)  ! w0 on path coarse
    real(rp) :: Z_COARSE(max_c)       ! Z_PSIG & Z_min & surface zeta
    real(rp), target :: Z_GLGRID_O(maxvert) ! Zeta on initial glGrid surfs, original
    real(rp), target :: Z_GLGRID_R(maxvert+nxg) ! Zeta on initial glGrid surfs, revised
    real(rp) :: Z_IG(size(Z_PSIG)+2)  ! Z_PSIG + min Zeta + Earth intersection
    real(rp) :: Z_PATH(max_f)         ! Zeta on fine grid path tangent grid and
                                      ! species grids, sans duplicates.
    real(rp) :: BETA_PATH_C(max_c,no_mol)  ! Beta on path coarse
    real(rp) :: BETA_PATH_F(max_f,no_mol)  ! Beta on path fine
    real(rp) :: D_DELTA_DF(s_a*max_c,size(grids_f%values)) ! Incremental 
                                      ! opacity derivative schlep from drad_tran_dt
                                      ! to get_d_deltau_pol_df.  Path x SVE.
    real(rp) :: D_T_SCR_dT(max_c,s_t*sv_t_len)  ! D Delta_B in some notes
                                      ! path x state-vector-components
    real(rp) :: D2X_DXDT(no_tan_hts,s_t*sv_t_len)  ! (No_tan_hts, nz*np)
    real(rp) :: DBETA_DN_PATH_C(max_c,s_td*size(fwdModelConf%lineWidth_TDep)) ! dBeta_dn on coarse grid
    real(rp) :: DBETA_DN_PATH_F(max_f,s_td*size(fwdModelConf%lineWidth_TDep)) ! dBeta_dn on fine grid
    real(rp) :: DBETA_DT_PATH_C(max_c,s_t*no_mol)  ! dBeta_dT on coarse grid
    real(rp) :: DBETA_DT_PATH_F(max_f,s_t*no_mol)  ! dBeta_dT on fine grid
    real(rp) :: DBETA_DV_PATH_C(max_c,s_lc*size(fwdModelConf%lineCenter)) ! dBeta_dv on coarse grid
    real(rp) :: DBETA_DV_PATH_F(max_f,s_lc*size(fwdModelConf%lineCenter)) ! dBeta_dv on fine grid
    real(rp) :: DBETA_DW_PATH_C(max_c,s_lw*size(fwdModelConf%lineWidth)) ! dBeta_dw on coarse grid
    real(rp) :: DBETA_DW_PATH_F(max_f,s_lw*size(fwdModelConf%lineWidth)) ! dBeta_dw on fine grid
    real(rp) :: DH_DT_PATH(max_f,s_t*sv_t_len)     ! dH/dT on path
    real(rp) :: DH_DT_PATH_C(max_c,s_t*sv_t_len)   ! DH_DT_PATH on coarse grid
    real(rp) :: DH_DT_PATH_F(max_f,s_t*sv_t_len)   ! DH_DT_PATH on fine grid
    real(rp), target :: DHDZ_GLGRID_O(maxVert,no_sv_p_t) ! dH/dZ on glGrid surfs, original
    real(rp), target :: DHDZ_GLGRID_R(maxVert+nxg,no_sv_p_t) ! dH/dZ on glGrid surfs, revised
    real(rp) :: DX_DT(no_tan_hts,s_t*sv_t_len)     ! (No_tan_hts, nz*np)
    real(rp) :: ETA_FZP(max_f,size(grids_f%values)) ! Eta_z x Eta_p * Eta_f
    real(rp) :: ETA_IWC_ZP(max_f,grids_iwc%p_len)
    real(rp) :: ETA_Mag_ZP(max_f,grids_mag%p_len)  ! Eta_z x Eta_p
    real(rp) :: ETA_ZP(max_f,grids_f%p_len)        ! Eta_z x Eta_p
    real(rp) :: ETA_ZXP_N(max_f,size(grids_n%values)) ! Eta_z x Eta_p for N
    real(rp) :: ETA_ZXP_T(max_f,s_t*sv_t_len)      ! Eta_t_z x Eta_t_p
    real(rp) :: ETA_ZXP_T_C(max_c,s_t*sv_t_len)    ! ETA_ZXP_T on coarse grid
    real(rp) :: ETA_ZXP_T_F(max_f,s_t*sv_t_len)    ! ETA_ZXP_T on fine grid
    real(rp) :: ETA_ZXP_V(max_f,size(grids_v%values)) ! Eta_z x Eta_p for V
    real(rp) :: ETA_ZXP_W(max_f,size(grids_w%values)) ! Eta_z x Eta_p for W
    real(rp), target :: H_GLGRID_O(maxVert,no_sv_p_t) ! H on glGrid surfs, original
    real(rp), target :: H_GLGRID_R(maxVert+nxg,no_sv_p_t) ! H on glGrid surfs, revised
    real(rp) :: IWC_PATH(max_f,s_i)                ! IWC on path
    real(rp), target :: MAG_PATH(s_p*max_f,4)      ! Magnetic field on path
    real(rp), target :: RAD_AVG_PATH(max_c,s_pfa*noUsedChannels) ! Freq. Avgd.
                                                   ! LBL radiance along the path
    real(rp) :: RADIANCES(no_tan_hts,noUsedChannels) ! (Nptg,noChans)
    real(rp) :: SPECT_N_PATH(max_f,size(fwdModelConf%lineWidth_TDep)) ! Line Width Temperature Dependence
    real(rp) :: SPECT_V_PATH(max_f,size(fwdModelConf%lineCenter)) ! Line Center
    real(rp) :: SPECT_W_PATH(max_f,size(fwdModelConf%lineWidth)) ! Line Width
    real(rp) :: SPS_PATH(max_f,no_mol)             ! species on path
    real(rp), target :: T_GLGRID_O(maxVert,no_sv_p_t) ! Temp on glGrid surfs, original
    real(rp), target :: T_GLGRID_R(maxVert+nxg,no_sv_p_t) ! Temp on glGrid surfs, revised
    real(rp) :: T_SCRIPT_PFA(max_c,s_pfa*noUsedChannels) ! Delta_B in some notes
    real(rp) :: TT_PATH(max_f,s_i)                 ! TScat on path along the LOS
    real(r8) :: VMRARRAY(no_mol,s_i*n_t_zeta)      ! The VMRs

    ! Temporary space for DACS radiances if we're doing frequency averaging
    ! and there are any LBL molecules and any derivatives are calculated
    real(rp) :: DACsStaging2(lbound(fwdModelConf%DACsStaging,1): &
      &                      ubound(fwdModelConf%DACsStaging,1), &
      & merge(1,0,fwdModelConf%do_freq_avg .and. &
      &           any(fwdModelConf%anyLBL((fwdModelConf%sidebandStart+3)/2: &
      &                                   (fwdModelConf%sidebandStop+3)/2))) * &
      &   max(s_t,s_a,s_lc,s_lw,s_td) * & ! merge(1,0,any_der)
      &   max(sv_t_len,size(grids_f%values), &
      &       size(grids_w%values), size(grids_n%values), size(grids_v%values)), &
      & size(fwdModelConf%usedDACSSignals) )

    real(rp), target :: DH_DT_GLGRID_O(maxVert,n_t_zeta,s_t*no_sv_p_t)
    real(rp), target :: DH_DT_GLGRID_R(maxVert+nxg,n_t_zeta,s_t*no_sv_p_t)

    complex(rp) :: D_RAD_POL_DF(2,2,s_p*s_a*size(grids_f%values)) ! From mcrt_der
    complex(rp) :: D_RAD_POL_DT(2,2,s_p*s_t*sv_t_len) ! From mcrt_der
    complex(rp) :: RAD_POL(2,2)  ! polarized radiance output of mcrt for one freq and pointing
      ! (-1,:,:) are Sigma_-, (0,:,:) are Pi, (+1,:,:) are Sigma_+
    complex(rp) :: ALPHA_PATH_POLARIZED(-1:1,s_p*max_c)
    complex(rp) :: ALPHA_PATH_POLARIZED_F(-1:1,s_p*max_f)
    complex(rp) :: dALPHA_dT_POLARIZED_PATH_C(-1:1,s_p*s_t*max_c)
    complex(rp) :: dALPHA_dT_POLARIZED_PATH_F(-1:1,s_p*s_t*max_f)
    complex(rp) :: BETA_PATH_POLARIZED(-1:1,s_p*max_c,no_mol)
    complex(rp) :: BETA_PATH_POLARIZED_F(-1:1,s_p*max_f,no_mol)
    complex(rp) :: dBETA_dT_POLARIZED_PATH_C(-1:1,s_p*s_t*max_c,no_mol)
    complex(rp) :: dBETA_dT_POLARIZED_PATH_F(-1:1,s_p*s_t*2*max_f,no_mol)
    complex(rp) :: DE_DF(2,2,s_p*s_a*max_c,size(grids_f%values)) ! DE/Df in Michael's notes
    complex(rp) :: DE_DT(2,2,s_p*s_t*max_c,sv_t_len) ! DE/DT in Michael's notes
    complex(rp) :: DELTAU_POL(2,2,s_p*max_c) ! E in Michael's notes
!   complex(rp) :: DINCOPTDEPTH_POL_DT(2,2,s_p*s_t*max_c) ! D Incoptdepth_Pol / DT
!   complex(rp) :: GL_DELTA_POLARIZED(-1:1,s_p*max_f)
    complex(rp) :: INCOPTDEPTH_POL(2,2,s_p*max_c)
    complex(rp) :: PROD_POL(2,2,s_p*max_c)   ! P in Michael's notes
    complex(rp) :: TAU_POL(2,2,s_p*max_c)    ! Tau in Michael's notes

    real(rp) :: EST_SCGEOCALT(no_tan_hts) ! Est S/C geocentric altitude /m
    real(rp) :: EST_LOS_VEL(no_tan_hts)   ! Est S/C line-of-sight velocity M/S
    real(rp) :: TAN_D2H_DHDT(s_t*sv_t_len)
    real(rp) :: TAN_PHI(no_tan_hts)
    real(rp), target :: DDHIDHIDTL0_O(maxVert,n_t_zeta,s_t*no_sv_p_t) ! Original
    real(rp), target :: DDHIDHIDTL0_R(maxVert+nxg,n_t_zeta,s_t*no_sv_p_t) ! Revised
    real(r4) :: K_ATMOS(noUsedChannels,no_tan_hts,s_a*size(grids_f%values))
    ! Channels x pointings x grid values == frequencies x surfaces x instances x molecules:
    real(r4) :: K_SPECT_DN(noUsedChannels,no_tan_hts,s_td*size(grids_n%values))
    real(r4) :: K_SPECT_DV(noUsedChannels,no_tan_hts,s_lc*size(grids_v%values))
    real(r4) :: K_SPECT_DW(noUsedChannels,no_tan_hts,s_lw*size(grids_w%values))
    real(r4) :: K_TEMP(noUsedChannels,no_tan_hts,s_t*sv_t_len)

    integer, pointer :: C_INDS(:)   ! Indices on coarse grid
    integer, pointer :: GL_INDS(:)  ! Index of GL points -- subset of f_inds
    integer, pointer :: LineCenter_IX(:) ! Where are line center offsets?
    integer, pointer :: LineWidth_IX(:)  ! Where are line width offsets?
    integer, pointer :: LineWidth_TDep_IX(:)  ! Where are line width TDep offsets?
    integer, pointer :: USEDDACSSIGNALS(:) ! Indices in FwdModelConf of signals
                                    ! for our dacs

    real(rp) :: E_RFLTY       ! Earth reflectivity at given tan. point
    real(rp), save :: E_Stop  = 1.0_rp ! X for which Exp(X) is too small to worry
    real(rp) :: H_Surf, H_Tan ! Height at surface and tangent
    real(rp) :: Min_Zeta      ! Minimum zeta along the path
    real(rp) :: Min_Phi       ! Phi at which minimum zeta occurs
    real(rp), parameter :: Min_Phi_Tol = 0.25 * gx(1)**2 ! First GL point
    real(rp) :: TAN_HT        ! Height at the tangent, from equivalent Earth center
    real(rp) :: R             ! real variable for various uses
    real(rp) :: REQ           ! Equivalent Earth Radius
    real(rp) :: ROT(3,3)      ! ECR-to-FOV rotation matrix
    real(rp) :: Vel_Cor       ! Velocity correction due to Vel_z, 1 - Vel_z/c
    real(rp) :: Vel_Rel       ! Vel_z / c

    real(rp), pointer :: CT(:)           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(r8), pointer :: FREQUENCIES(:)  ! Frequencies to compute for
    real(rp), pointer :: H(:)            ! Magnetic field on path, in
                                         ! IFOVPP
    real(rp), pointer :: RADV(:)         ! Radiances for 1 pointing on
                                         ! Freq_Grid
    real(rp), pointer :: STCP(:)         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), pointer :: STSP(:)         ! Sin(Theta) Sin(Phi)
    real(rp), pointer :: Cext_PATH(:,:)  ! Cloud extinction on path
    real(rp), pointer :: DACsStaging(:,:) ! Temporary space for DACS radiances

    ! GL Grid quantities to/from Hydrostatic
    real(rp), pointer :: DDHIDHIDTL0(:,:,:) ! Either DDHIDHIDTL0_O or DDHIDHIDTL0_R
    real(rp), pointer :: DH_DT_GLGRID(:,:,:) ! Either DH_DT_GLgrid_O or DH_DT_GLgrid_R
    real(rp), pointer :: DHDZ_GLGRID(:,:) ! Either DHDZ_GLgrid_O or DHDZ_GLgrid_R
    real(rp), pointer :: P_GLGRID(:)     ! Either P_GLgrid_O or P_GLgrid_R
    real(rp), pointer :: H_GLGRID(:,:)   ! Either H_GLgrid_O or H_GLgrid_R
    real(rp), pointer :: T_GLGRID(:,:)   ! Either T_GLgrid_O or T_GLgrid_R
    real(rp), pointer :: Z_GLGRID(:)     ! Either Z_GLgrid_O or Z_GLgrid_R

    ! Incremental opacity derivatives, Path X SVE:
    real(rp), pointer :: ETA_Tscat(:,:)    !
    real(rp), pointer :: ETA_Tscat_ZP(:,:) !
    real(rp), pointer :: ETA_Salb(:,:)     !
    real(rp), pointer :: ETA_Salb_ZP(:,:)  !
    real(rp), pointer :: ETA_Cext(:,:)     !
    real(rp), pointer :: ETA_Cext_ZP(:,:)  !
    real(rp), pointer :: INC_RAD_PATH(:,:) ! Incremental radiance along the path
    real(rp), pointer :: K_ATMOS_FRQ(:,:)  ! dI/dVMR, ptg.frq X vmr-SV
    real(rp), pointer :: K_SPECT_DN_FRQ(:,:) ! ****
    real(rp), pointer :: K_SPECT_DV_FRQ(:,:) ! ****
    real(rp), pointer :: K_SPECT_DW_FRQ(:,:) ! ****
    real(rp), pointer :: K_TEMP_FRQ(:,:)   ! dI/dT, ptg.frq X T-SV
    real(rp), pointer :: Salb_PATH(:,:)    ! Single Scattering Albedo on path
    real(rp), pointer :: T_SCRIPT_LBL(:,:) ! Delta_B in some notes
    real(rp), pointer :: Tscat_PATH(:,:)   ! TScat on path

    ! Used only to schlep from Both_Sidebands_Setup to Convolution
    real(rp) :: DH_DZ_OUT(ptan%template%nosurfs)
    real(rp) :: DX_DH_OUT(ptan%template%nosurfs)
    real(rp) :: DXDT_SURFACE(1,s_t*sv_t_len)
    real(rp) :: DXDT_TAN(ptan%template%nosurfs,sv_t_len)
    real(rv), pointer :: L1BMIF_TAI(:,:)   ! MIF Times
    real(rv), pointer :: MIFDEADTIME(:,:)  ! Not collecting data
    real(rp) :: surf_angle(1)
    real(rp) :: TAN_CHI_OUT(ptan%template%nosurfs)

    real(rp) :: earthradc ! minor axis of orbit plane projected Earth ellipse
    real(rp) :: earthradc_sq ! earthradc**2

    type (VectorValue_T), pointer :: BOUNDARYPRESSURE
    type (VectorValue_T), pointer :: CLOUDWATER    ! Profiles
    type (VectorValue_T), pointer :: EARTHREFL     ! Earth reflectivity
    type (VectorValue_T), pointer :: ECRtoFOV      ! Rotation matrices
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: GPH           ! Geopotential height
    type (VectorValue_T), pointer :: LOSVEL        ! Line of sight velocity
    type (VectorValue_T), pointer :: ORBINCLINE    ! Orbital inclination
    type (VectorValue_T), pointer :: REFGPH        ! Reference geopotential height
    type (VectorValue_T), pointer :: SCGEOCALT     ! S/C geocentric altitude /m
    type (VectorValue_T), pointer :: SIZEDISTRIBUTION ! Integer really
    type (VectorValue_T), pointer :: SPACERADIANCE ! Emission from space
    type (VectorValue_T), pointer :: SurfaceHeight ! km above mean sea level
    type (VectorValue_T), pointer :: THISRADIANCE  ! A radiance vector quantity
    type (VectorValue_T), pointer :: VMR           ! Quantity

    type (Signal_T), pointer :: FIRSTSIGNAL        ! The first signal we're dealing with

    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS ! Freq. indep. stuff

    type (Grids_T) :: Grids_Tscat ! All the coordinates for scaterring source function
    type (Grids_T) :: Grids_Salb ! All the coordinates for single scaterring albedo
    type (Grids_T) :: Grids_Cext ! All the coordinates for cloud extinction

!  The 'all_radiometers grid file' approach variables declaration:

    real(rp) :: max_ch_freq_grid, min_ch_freq_grid
    real(r8) :: TOL = 1.D-190

! *** Beta & Molecules grouping variables:
    type(beta_group_t), pointer :: Beta_Group(:) ! from FwdModelConf%Beta_Group

    type(tau_t) :: Tau_LBL, Tau_PFA

! Channel information from the signals database as specified by fwdModelConf
    type(channels_T), pointer, dimension(:) :: Channels 

! Scattering source function for each temperature surface
    type (VectorValue_T) :: scat_src
    type (VectorValue_T) :: scat_alb
    type (VectorValue_T) :: cld_ext

    ! Executable code --------------------------------------------------------
    ! ------------------------------------------------------------------------
!   Print *, '** Enter ForwardModel, MAF =',fmstat%maf   ! ** ZEBUG

    if ( toggle(emit) ) & ! set by -f command-line switch
      & call trace_begin ( 'Full ForwardModel, MAF=', index=fmstat%maf )

    ! Set flags from command-line switches
    clean = index(switches, 'clean') /= 0
    do_zmin = index(switches, 'nozm') == 0 ! Do minimum zeta unless told otherwise
    print_Grids = index(switches, 'grids') /= 0
    print_Mag = index(switches, 'mag') /= 0
    print_Min_Zeta = index(switches, 'zmin') /= 0
    print_Ptg = index(switches,'ptg') /= 0
    print_Rad = index(switches, 'rad') /= 0
    print_Seez = index(switches, 'seez') /= 0
    print_TauL = index(switches, 'taul') /= 0
    print_TauP = index(switches, 'taup') /= 0

    ! Nullify all our pointers that are allocated because the first thing
    ! Allocate_Test does is ask if they're associated.  If we don't nullify
    ! them, they're undefined, i.e., junk that might be mistaken for
    ! associated.

    nullify ( cext_path, do_calc_Cext, do_calc_Cext_zp, do_calc_Salb, &
      & do_calc_Salb_zp, do_calc_tscat, do_calc_tscat_zp, eta_cext,   &
      & eta_cext_zp, eta_salb, eta_salb_zp, eta_tscat, eta_tscat_zp,  &
      & frequencies, inc_rad_path, k_atmos_frq, k_spect_dn_frq,       &
      & k_spect_dv_frq, k_spect_dw_frq, k_temp_frq, RadV, salb_path,  &
      & tscat_path, t_script_lbl, vmr )

    ! Nullify pointers that are used to control whether calculations get done
    nullify ( linecenter_ix, linewidth_ix, linewidth_tdep_ix )

    call both_sidebands_setup

    ! Compute hydrostatic grid -----------------------------------------------
    ! If this is moved into BOTH_SIDEBANDS_SETUP the run time increases by
    ! a large unexplainable amount, at least in LF95.

    ! Insert into bill's 2d hydrostatic equation.
    ! The phi input for this program are the orbit plane projected
    ! geodetic locations of the temperature phi basis -- not necessarily
    ! the tangent phi's, which may be somewhat different.

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_Begin ( 'ForwardModel.Hydrostatic' )

    ! Temperature's windowStart:windowFinish are correct here.
    ! RefGPH and temperature have the same horizontal basis.
    ! RefGPH is in meters, but Two_D_Hydrostatic wants it in km.

    if ( temp_der ) then
      call two_d_hydrostatic ( Grids_tmp, &
        &  (/ (refGPH%template%surfs(1,1), j=windowStart,windowFinish) /), &
        &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
        &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
        &  dhdz_glgrid, dh_dt_glgrid, DDHDHDTL0=ddhidhidtl0 )
    else
      call two_d_hydrostatic ( Grids_tmp, &
        &  (/ (refGPH%template%surfs(1,1), j=windowStart,windowFinish) /), &
        &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
        &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
        &  dhdz_glgrid )
    end if

    if ( toggle(emit) .and. levels(emit) > 0 ) then
      call Trace_End ( 'ForwardModel.Hydrostatic' )
      call Trace_Begin ( 'ForwardModel.SidebandLoop' )
    end if

    ! Loop over sidebands ----------------------------------------------------
    do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2

      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call Trace_Begin ( 'ForwardModel.Sideband ', index=thisSideband )

      sx = ( thisSideband + 3 ) / 2 ! [-1,+1] => [1,2]

      ! Work out Frequency averaging / LBL / PFA / Derivatives steering.
      frq_avg_sel = 0
      if ( fwdModelConf%do_freq_avg ) frq_avg_sel = ior(frq_avg_sel, 8)
      if ( fwdModelConf%anyPFA(sx) .and. .not. &
        &  fwdModelConf%anyLBL(sx) )  frq_avg_sel = ior(frq_avg_sel, 8)
      if ( fwdModelConf%anyPFA(sx) )  frq_avg_sel = ior(frq_avg_sel, 4)
      if ( fwdModelConf%anyLBL(sx) )  frq_avg_sel = ior(frq_avg_sel, 2)
      if ( any_der )                  frq_avg_sel = ior(frq_avg_sel, 1)

      ! Now, allocate gl_slabs arrays
      call allocateSlabs ( gl_slabs, max_f, &
        & fwdModelConf%catalog(thisSideband,:fwdModelConf%cat_size(sx)), &
        & moduleName, temp_der )

      call frequency_setup_1

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_Begin ( 'ForwardModel.PointingLoop' )

      ! Loop over pointings --------------------------------------------------
      do ptg_i = 1, no_tan_hts
        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call Trace_Begin ( 'ForwardModel.Pointing ', index=ptg_i )

        if ( FwdModelConf%polarized ) &
          & mif = minloc(abs(tan_press(ptg_i) - &
          &                  ptan%values(:ptan%template%nosurfs,maf)),1)

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_Begin ( 'ForwardModel.MetricsEtc' )

        ! Compute where the tangent is, the equivalent Earth radius, tangent
        ! height and surface height, and determine whether the ray
        ! intersects the Earth surface
        tan_ind_c = max(1,ptg_i-surfaceTangentIndex+1) ! On coarse grid
        nz_ig = nlvl
        if ( associated(surfaceHeight) ) then
          call tangent_metrics ( tan_phi(ptg_i), Grids_tmp%phi_basis, z_psig, &
            &                    h_glgrid, earthradc_sq,     & ! in
            &                    tan_ind_c, nz_ig,           & ! inout
            &                    req, h_surf, h_tan, z_ig,   & ! output
            &                    surf_height=surfaceHeight%values(1,:) ) ! opt
        else if ( ptg_i < surfaceTangentIndex ) then
          call tangent_metrics ( tan_phi(ptg_i), Grids_tmp%phi_basis, z_psig, &
            &                    h_glgrid, earthradc_sq,     & ! in
            &                    tan_ind_c, nz_ig,           & ! inout
            &                    req, h_surf, h_tan, z_ig,   & ! output
            &                    Tan_Press=tan_press(ptg_i), & ! optional
            &                    Surf_Temp=temp%values(1,windowstart:windowfinish) )
        else
          call tangent_metrics ( tan_phi(ptg_i), Grids_tmp%phi_basis, z_psig, &
            &                    h_glgrid, earthradc_sq,     & ! in
            &                    tan_ind_c, nz_ig,           & ! inout
            &                    req, h_surf, h_tan, z_ig )    ! output
        end if

        nz_if = (nz_ig-1) * ngp1 + 1                ! On Z_GLgrid
        tan_ind_f = (tan_ind_c-1) * ngp1 + 1        ! On Z_GLgrid
        tan_pt_f = (nz_ig-1) * ngp1 + 2 - tan_ind_f ! fine path tangent index
        tan_pt_c = (tan_pt_f + Ng) / Ngp1           ! coarse path tangent index
        npc = 2 * tan_pt_c
        npf = 2 * tan_pt_f

        if ( nz_ig == nlvl ) then
          z_coarse(:tan_pt_c) = z_psig(nlvl:tan_ind_c:-1)
          z_coarse(tan_pt_c+1:npc) = z_psig(tan_ind_c:nlvl)
        end if

        if ( h_tan < 0.0 ) then ! Handle Earth-intersecting ray
          e_rflty = earthRefl%values(1,1)
          if ( nz_ig > nlvl ) then ! Added a new zeta
            z_coarse(:tan_pt_c) = z_ig(nz_ig:tan_ind_c:-1)
            z_coarse(tan_pt_c+1:npc) = z_ig(tan_ind_c:nz_ig)
            ! Use the "revised" GL-grid arrays, which depend upon the new zeta.
            ddhidhidtl0 => ddhidhidtl0_r(:nz_if,:,:)
            dh_dt_glgrid => dh_dt_glgrid_r(:nz_if,:,:)
            dhdz_glgrid => dhdz_glgrid_r(:nz_if,:)
            h_glgrid => h_glgrid_r(:nz_if,:)
            p_glgrid => p_glgrid_r(:nz_if)
            t_glgrid => t_glgrid_r(:nz_if,:)
            z_glgrid => z_glgrid_r(:nz_if)
            call compute_GL_grid ( z_ig(:nz_ig), p_glgrid, z_glgrid )
            if ( temp_der ) then
              call two_d_hydrostatic ( Grids_tmp, &
                &  (/ (refGPH%template%surfs(1,1), j=windowStart,windowFinish) /), &
                &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
                &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
                &  dhdz_glgrid, dh_dt_glgrid, DDHDHDTL0=ddhidhidtl0 )
            else
              call two_d_hydrostatic ( Grids_tmp, &
                &  (/ (refGPH%template%surfs(1,1), j=windowStart,windowFinish) /), &
                &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
                &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
                &  dhdz_glgrid )
            end if
          end if
        else ! Not an Earth-intersecting ray
          e_rflty = 1.0
          ! Use the "original" GL-grid arrays -- the ones filled before the
          ! sideband loop.  These depend only upon the original Zetas
          ! (Z_PSIG) derived from the L2CF using make_z_grid.
          ddhidhidtl0 => ddhidhidtl0_o
          dh_dt_glgrid => dh_dt_glgrid_o
          dhdz_glgrid => dhdz_glgrid_o
          h_glgrid => h_glgrid_o
          p_glgrid => p_glgrid_o
          t_glgrid => t_glgrid_o
          z_glgrid => z_glgrid_o
        end if

        ! Get H_Path and Phi_Path on the fine grid.
        call Height_Metrics ( tan_phi(ptg_i), tan_ind_f, Grids_tmp%phi_basis, &
          &  z_glgrid, h_glgrid, req, h_surf, h_tan, & ! in
          &  h_path(1:npf), phi_path(1:npf) )          ! out

        tan_ht = h_path(tan_pt_f) ! Includes Earth radius

        if ( do_zmin ) then
          ! Get minimum zeta on the path
          call Get_Min_Zeta ( Grids_tmp%phi_basis, h_glgrid(tan_ind_f,:), &
                            & t_glgrid(tan_ind_f,:), z_glgrid(tan_ind_f), &
                            & phi_path, tan_ind_f, tan_ht,                &
                            & min_zeta, min_phi, min_index )

          ! Add minimum zeta to the path
          if ( min_index > 0 ) then ! minimum zeta not at or near tangent point
            min_index_c = min_index - mod(min_index-1,ngp1)
            if ( min_index > tan_pt_f ) min_index_c = min_index_c + 1
            if ( min(abs(min_phi-phi_path(min_index)), &
              &      abs(min_phi-phi_path(min_index+1))) <= &
              &  min_phi_tol * (phi_path(min_index_c+ngp1)-phi_path(min_index_c)) ) then
              ! Min zeta very close to an existing point
              if ( abs(min_phi-phi_path(min_index)) > abs(min_phi-phi_path(min_index+1)) ) &
                & min_index = min_index+1
              ! Min_index is now the index of the point
              if ( min_index == min_index_c .or. min_index == min_index_c+ngp1 ) then
                ! Min_index is at a coarse grid point.
                ! All we need to do is change the zeta.
                min_index_c = min_index
                i = min_index / ngp1 + 1
                if ( print_Min_Zeta ) then
                  call output ( i, before='Replacing Z_Coarse(' )
                  call output ( z_coarse(i), before=') = ' )
                  call output ( min_zeta, before=' with minimum Zeta = ', advance='yes' )
                end if
                z_coarse(i) = min_zeta
                min_index_c = 0 ! Indicate nothing more to do
              end if
            end if
            if ( min_index_c > 0 ) then
              ! Min zeta not near an existing coarse point: add one to the path
              npc = npc + 1
              npf = npf + ngp1
              i = min_index_c / ngp1 + 1
              z_coarse(i+1:npc) = z_coarse(i:npc-1) ! Make room
              z_coarse(i) = min_zeta
              if ( print_Min_Zeta ) then
                call output ( min_zeta, before='Added minimum Zeta = ' )
                call output ( i, before=' after ', advance='yes' )
              end if
              if ( min_index_c < tan_pt_f ) then
                tan_pt_c = tan_pt_c + 1
                tan_pt_f = tan_pt_f + ngp1
              end if
            end if
          end if
        end if

        ! Compute Gauss Legendre (GL) grid ----------------------------------
        call compute_GL_grid ( z_coarse(:tan_pt_c), p_path(:tan_pt_f), &
          &                    z_path(:tan_pt_f) )
        call compute_GL_grid ( z_coarse(tan_pt_c+1:npc), p_path(tan_pt_f+1:npf), &
          &                    z_path(tan_pt_f+1:npf) )

        ! This is not pretty but we need some coarse path extraction indices
        c_inds => c_inds_b(:npc)
        c_inds = (/(i*Ngp1-Ng,i=1,tan_pt_c),((i-1)*Ngp1-Ng+1,i=tan_pt_c+1,npc)/)
        ! And some fine path extraction indices
        do_gl(1:npc:npc-1) = .false.; do_gl(2:npc-1) = .true.
        call get_gl_inds ( do_gl(:npc), tan_pt_c, f_inds, cg_inds, nglMax, ncg )

        del_zeta(1:npc:npc-1) = 0.0_rp ! First and last ones
        del_zeta(2:tan_pt_c) = 0.5_rp * ( z_path(c_inds(1:tan_pt_c-1)) - &
          &                               z_path(c_inds(2:tan_pt_c)) )
        del_zeta(tan_pt_c+1:npc-1) = 0.5_rp * ( z_path(c_inds(tan_pt_c+2:npc)) - &
          &                                     z_path(c_inds(tan_pt_c+1:npc-1)) )

        ! Do phi refractive correction
        if ( h2o_ind == 0 ) then
          call refractive_index ( p_path(1:npf), t_path(1:npf), n_path_f(1:npf) )
          n_path_c(:npc) = n_path_f(c_inds)
        end if
        if ( FwdModelConf%refract ) then
          ! Get t_path (and dhdz_path, which we don't need yet)
          call more_metrics ( tan_phi(ptg_i), tan_ind_f, tan_pt_f,          &
            &  Grids_tmp%phi_basis, z_glgrid, t_glgrid, h_path(1:npf),      &
            &  dhdz_glgrid, phi_path(1:npf),                                &
            &  t_path(1:npf), dhdz_path(1:npf) )
          ! Compute refractive index on the path.
          if ( h2o_ind > 0 ) then
            ! Compute eta_zp & do_calc_zp (Zeta & Phi only) for water
            call comp_eta_docalc_no_frq ( Grids_f, z_path(1:npf), &
              &  phi_path(1:npf), eta_zp(1:npf,:), do_calc_zp(1:npf,:) )
            call comp_1_sps_path_no_frq ( Grids_f, h2o_ind, eta_zp(1:npf,:), &
              & sps_path(1:npf,h2o_ind) )
            call refractive_index ( p_path(1:npf), t_path(1:npf), &
              &  n_path_f(1:npf),  &
              &  h2o_path=sps_path(1:npf, h2o_ind) )
          end if
          ! Do the refractive correction.  Use t_path to store the correction,
          ! since we're going to recompute t_path right away.
          n_path_f(1:npf) = min ( n_path_f(1:npf), MaxRefraction )
          call phi_refractive_correction ( tan_pt_f, n_path_f(1:npf), &
            & h_path(1:npf), t_path(1:npf) )
          phi_path(:tan_pt_f) = phi_path(:tan_pt_f) - t_path(:tan_pt_f)
          phi_path(tan_pt_f+1:npf) = phi_path(tan_pt_f+1:npf) + t_path(tan_pt_f+1:npf)
        end if

        ! Now get other metrics-related quantities, t_path, dhdz_path, dhdt_path
        if ( temp_der ) then
          call more_metrics ( tan_phi(ptg_i), tan_ind_f, tan_pt_f,          &
            &  Grids_tmp%phi_basis, z_glgrid, t_glgrid, h_path(1:npf),      &
            &  dhdz_glgrid, phi_path(1:npf),                                &
            &  t_path(1:npf), dhdz_path(1:npf),                             &
            !  Stuff for temperature derivatives:
            &  DHTDTL0 = tan_dh_dt, DDHIDHIDTL0 = ddhidhidtl0,              &
            &  DDHTDHTDTL0 = tan_d2h_dhdt, DHIDTLM = dh_dt_glgrid,          &
            &  DHITDTLM = dh_dt_path(1:npf,:),                              &
            &  T_DERIV_FLAG = Grids_tmp%deriv_flags,                        &
            &  Z_BASIS = Grids_tmp%zet_basis,                               &
            &  ETA_ZXP = eta_zxp_t(1:npf,:),                                &
            &  DO_CALC_T = do_calc_t(1:npf,:),                              &
            &  DO_CALC_HYD = do_calc_hyd(1:npf,:) )
          dh_dt_path_c(1:npc,:) = dh_dt_path(c_inds,:)
          do_calc_hyd_c(1:npc,:) = do_calc_hyd(c_inds,:)
          do_calc_t_c(1:npc,:) = do_calc_t(c_inds,:)
          eta_zxp_t_c(1:npc,:) = eta_zxp_t(c_inds,:)
          t_der_path_flags(1:npf) = any(do_calc_t(1:npf,:),2)
        else
          call more_metrics ( tan_phi(ptg_i), tan_ind_f, tan_pt_f,          &
            &  Grids_tmp%phi_basis, z_glgrid, t_glgrid, h_path(1:npf),      &
            &  dhdz_glgrid, phi_path(1:npf),                                &
            &  t_path(1:npf), dhdz_path(1:npf) )
        end if

!       dhdz_gw_path(f_inds) = dhdz_path(f_inds) * (/ ( gw, i = 1, nglMax/ng ) /)
        do i = 1, nglMax, ng ! Avoid a temp for (/ ( gw, i = 1, nglMax/ng ) /)
          dhdz_gw_path(f_inds(i:i+ng-1)) = dhdz_path(f_inds(i:i+ng-1)) * gw
        end do
        h_path_c(1:npc) = h_path(c_inds)
        t_path_c(1:npc) = t_path(c_inds)

        ! Compute the eta_zp & do_calc_zp (for Zeta & Phi only)

        call comp_eta_docalc_no_frq ( Grids_f, z_path(1:npf), &
          &  phi_path(1:npf), eta_zp(1:npf,:), do_calc_zp(1:npf,:) )

        ! Compute sps_path with a FAKE frequency, mainly to get the
        ! WATER (H2O) contribution for refraction calculations, but also
        ! to compute sps_path for all those with no frequency component

      ! Frq = 0.0_r8
        call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
          & 0.0_r8, eta_zp(1:npf,:),                &
          & do_calc_zp(1:npf,:), sps_path(1:npf,:), &
          & do_calc_fzp(1:npf,:), eta_fzp(1:npf,:) )

        if ( h2o_ind > 0 ) then
          ! Even if we did the refractive correction we need to do this,
          ! because the refractive correction changes phi_path.
          call refractive_index ( p_path(c_inds), &
            &  t_path_c(1:npc), n_path_c(1:npc),  &
            &  h2o_path=sps_path(c_inds, h2o_ind) )
        end if

        n_path_c(1:npc) = min ( n_path_c(1:npc), MaxRefraction )

        if ( size(fwdModelConf%lineCenter) > 0 ) then
          call comp_eta_docalc_no_frq ( grids_v, z_path(1:npf), &
            & phi_path(1:npf), eta_zxp_v(1:npf,:), do_calc_v(1:npf,:) )
          call comp_sps_path_no_frq ( grids_v, eta_zxp_v(1:npf,:), &
            & spect_v_path(1:npf,:) )
          lineCenter_ix => beta_group%lbl(sx)%spect_der_ix(lineCenter)
        end if
        if ( size(fwdModelConf%lineWidth) > 0 ) then
          call comp_eta_docalc_no_frq ( grids_w, z_path(1:npf), &
            & phi_path(1:npf), eta_zxp_w(1:npf,:), do_calc_w(1:npf,:) )
          do_calc_w_c(1:npc,:) = do_calc_w(c_inds,:)
          call comp_sps_path_no_frq ( grids_w, eta_zxp_w(1:npf,:), &
            & spect_w_path(1:npf,:) )
          lineWidth_ix => beta_group%lbl(sx)%spect_der_ix(lineWidth)
        end if
        if ( size(fwdModelConf%lineWidth_TDep) > 0 ) then
          call comp_eta_docalc_no_frq ( grids_n, z_path(1:npf), &
            & phi_path(1:npf), eta_zxp_n(1:npf,:), do_calc_n(1:npf,:) )
          call comp_sps_path_no_frq ( grids_n, eta_zxp_n(1:npf,:), &
            & spect_n_path(1:npf,:) )
          lineWidth_TDep_ix => beta_group%lbl(sx)%spect_der_ix(lineWidth_TDep)
        end if

        ! Special path quantities for cloud model
        if ( fwdModelConf%Incl_Cld ) then ! s_i == 1 here

          !set cloud parameters to zero
          iwc_path(1:npf,1) = 0.
          WC(1,1:npf)=iwc_path(1:npf,1)
          WC(2,1:npf)=0.
          IPSD(1:npf)=1000

          call comp_eta_docalc_no_frq ( Grids_Iwc, z_path(1:npf), &
            &  phi_path(1:npf), eta_iwc_zp(1:npf,:) )

          ! Compute IWC_PATH
          call comp_sps_path_no_frq ( Grids_iwc, eta_iwc_zp(1:npf,:), &
            & iwc_path(1:npf,:) )
          WC(1,1:npf)=iwc_path(1:npf,1)
        end if

        ! Special path quantities for Polarized (magnetic) model
        if ( FwdModelConf%polarized ) then

          call comp_eta_docalc_no_frq ( Grids_Mag, z_path(1:npf), &
            &  phi_path(1:npf), eta_mag_zp(1:npf,:) )

          ! Compute the first three components of MAG_PATH
          call comp_sps_path ( Grids_mag, 1, eta_mag_zp(1:npf,:), &
            & mag_path(1:npf,1:3) )

          rot = reshape(ECRtoFOV%values(9*mif-8:9*mif,maf), (/3,3/))

          do j = 1, npf
            ! Rotate mag_path from ECR to IFOVPP (R1A) coordinates.  Use
            ! the rotation matrix for the MIF nearest to the current
            ! pointing angle instead of interpolating.  They are nearly
            ! identical anyway.
            mag_path(j,1:3) = matmul ( rot, mag_path(j,1:3) )
            ! Put the magnitude of mag_path(j,1:3) in mag_path(j,4)
            mag_path(j,4) = sqrt(sum(mag_path(j,1:3)**2))
            ! Normalize mag_path(j,1:3).
            if ( mag_path(j,4) /= 0.0_rp ) then
              mag_path(j,1:3) = mag_path(j,1:3) / mag_path(j,4)
            else
              mag_path(j,1:3) = 0.0_rp
              mag_path(j,3) = 1.0_rp !arbitrarily, theta=0 for zero field
            end if
          end do

          ct => mag_path(1:npf,3)   ! cos(theta)
          stcp => mag_path(1:npf,1) ! sin(theta) cos(phi)
          stsp => mag_path(1:npf,2) ! sin(theta) sin(phi)
          h => mag_path(1:npf,4)    ! magnitude of magnetic field

          if ( print_Mag ) then
            call dump ( h, 'H', clean=clean )
            call dump ( ct, 'Cos(theta)', clean=clean )
            call dump ( stcp, 'Sin(theta) Cos(phi)', clean=clean )
            call dump ( stsp, 'Sin(theta) Sin(phi)', clean=clean )
          end if

        end if

        if ( temp_der ) then
         ! Ext_SCgeocAlt is in meters, but Get_Chi_Angles wants it in km.
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path_c(tan_pt_c),&
             & tan_ht-Req, tan_phi(ptg_i), Req, 0.0_rp, ptg_angles(ptg_i),      &
             & r, 1.0_rp, tan_dh_dt, tan_d2h_dhdt,       &
             & dx_dt(ptg_i,:), d2x_dxdt(ptg_i,:) )
        else
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path_c(tan_pt_c),&
             & tan_ht-Req, tan_phi(ptg_i), Req, 0.0_rp, ptg_angles(ptg_i),      &
             & r, 1.0_rp )
        end if

        n_path_c(1:npc) = n_path_c(1:npc) + 1.0_rp

        call comp_refcor ( tan_pt_c, h_path_c(:npc), n_path_c(:npc), &
                      &    tan_ht, del_s(:npc), ref_corr(:npc), ier )
        if ( ier /= 0 ) fmStat%flags = ior(fmStat%flags,b_refraction)

        ! We need path_dsdh on the fine grid for Gauss-Legendre or Gauss-
        ! Lobatto quadrature, and on the coarse grid except at the tangent
        ! point for trapezoidal quadrature and Gauss-Lobatto quadrature, so
        ! compute it everywhere except at the tangent point.  Besides, it's
        ! probably faster not to use a vector subscript to restrict it to
        ! the fine grid.

        path_dsdh(:tan_pt_f-1) = h_path(:tan_pt_f-1) / &
          & ( sqrt(h_path(:tan_pt_f-1)**2 - tan_ht**2 ) )
        path_dsdh(tan_pt_f+2:npf) = h_path(tan_pt_f+2:npf) / &
          & ( sqrt(h_path(tan_pt_f+2:npf)**2 - tan_ht**2 ) )
        path_dsdh(tan_pt_f:tan_pt_f+1) = 0.0

        dsdz_gw_path(f_inds(:nglMax)) = path_dsdh(f_inds(:nglMax)) * &
          & dhdz_gw_path(f_inds(:nglMax))

        ! We need dsdz = ds/dh * dh/dz, not multiplied by GW, for
        ! trapezoidal quadrature on the coarse grid.
        dsdz_c(:npc) = path_dsdh(c_inds) * dhdz_path(c_inds)

        ! Compute ALL the slabs_prep entities over the path's GL grid for this
        ! pointing & mmaf:
        if ( temp_der ) then
          call get_gl_slabs_arrays ( p_path(1:npf), t_path(1:npf), &
            &  est_los_vel(ptg_i), gl_slabs(1:npf,:), fwdModelConf%Do_1D, &
            &  spect_v_path, lineCenter_ix, &
            &  spect_w_path, lineWidth_ix, &
            &  spect_n_path, lineWidth_TDep_ix, &
            &  t_der_path_flags(1:npf) )
        else
          call get_gl_slabs_arrays ( p_path(1:npf), t_path(1:npf), &
            &  est_los_vel(ptg_i), gl_slabs(1:npf,:), fwdModelConf%Do_1D, &
            &  spect_v_path, lineCenter_ix, &
            &  spect_w_path, lineWidth_ix, &
            &  spect_n_path, lineWidth_TDep_ix )
        end if

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.

        if ( FwdModelConf%do_freq_avg .and. fwdModelConf%anyLBL(sx) ) &
          & call frequency_setup_2 ( vel_cor * &
          & PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies )

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_End ( 'ForwardModel.MetricsEtc' )

!{{\bfseries Notations}
! 
! $N_f$ is the number of points in the pointing frequency grid.\\
! $n$ is an index in the pointing frequency grid.\\
! $N_p$ is the number of points in the line-of-sight path.\\
! $i$ is an index in the line-of-sight path.  $c$ is a channel index.\\
! $s$ indicates a strong-line (LBL) result.
!   $w$ indicates a weak-line (PFA) result.\\
! $\delta I^\sigma_{iq} = \Delta B_{iq} \tau^\sigma_{iq}$ is the incremental
!  radiance contribution at the $i^{\text{th}}$ point along the line-of-sight
!  path, for $\sigma$ either $s$  or $w$ and $q$ either $c$ or $n$.
!
! There are four possible combinations of LBL, PFA and
! frequency-averaging.

!{{\bfseries Frequency averaged, LBL only}
!      \begin{equation*}\begin{split}
!       I^s_c = \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc}
!             \sum_{i=1}^{N_p} \delta I^s_{in}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I^s_n}{\partial x_k} =& 
!        \sum_{i=1}^{N_p}
!         \tau^s_{in} \frac{\partial \Delta B_{in}}{\partial x_k}
!         - \delta I^s_{in}
!          \sum_{j=1}^i \frac{\partial \delta^s_{jn}}{\partial x_k}\\
! %
!       \frac{\partial I_c}{\partial x_k} =&
!       \frac{\partial I^s_c}{\partial x_k} =
!        \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc}
!        \frac{\partial I^s_n}{\partial x_k}
!      \end{split}\end{equation*}

!{{\bfseries Monochromatic, LBL only}
!      \begin{equation*}
!       I^s_c = \sum_{i=1}^{N_p} \delta I^s_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I^s_c}{\partial x_k} =
!        \sum_{i=1}^{N_p}
!         \tau^s_{ic} \frac{\partial \Delta B_{ic}}{\partial x_k}
!         - \delta I^s_{ic}
!          \sum_{j=1}^i \frac{\partial \delta^s_{jc}}{\partial x_k}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!       \frac{\partial I^s_c}{\partial x_k}
!      \end{equation*}

!{{\bfseries Frequency averaged, PFA only}
!      \begin{equation*}
!       I_c = \sum_{i=1}^{N_p} \delta I^w_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!        \sum_{i=1}^{N_p}
!         \tau^w_{ic} \frac{\partial \Delta B_{ic}}{\partial x_k}
!           - \delta I^w_{ic}
!            \sum_{j=1}^i \frac{\partial \delta^w_{jc}}{\partial x_k}
!       \text{\phantom{xxxx}}
!      \end{equation*}

!{{\bfseries Frequency averaged, LBL and PFA}
!      \begin{equation*}
! %
!       \overline{\delta I^s_{ic}} =
!        \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc} \delta I^s_{in}
!       \text{\phantom{xxxx}}
! %
!       I_{ic} = \overline{\delta I^s_{ic}} \tau^w_{ic}
!       \text{\phantom{xxxx}}
! %
!       I_c = \sum_{i=1}^{N_p} I_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!       \frac{\partial I^s_c}{\partial x_k}
!        - \sum_{i=1}^{N_p} I_{ic}
!          \sum_{j=1}^i \frac{\partial \delta^w_{jc}}{\partial x_k}
!      \end{equation*}

!{{\bfseries Program variables}
! 
! \begin{tabular}{llll}
! $\Delta B_{in}$ is {\tt T\_Script\_LBL} &
! $\Delta B_{ic}$ is {\tt T\_Script\_PFA} &
! $\tau^s_{in}$ is {\tt Tau\_LBL} &
! $\tau^w_{ic}$ is {\tt Tau\_PFA}
! \\
! $\delta I^\sigma_{iq}$ is {\tt Inc\_Rad\_Path} &
! $\overline{\delta I^s_{ic}}$ is {\tt Rad\_Avg\_Path} & 
! $I_{ic}$ is also {\tt Rad\_Avg\_Path} &
! \\
! $\sum_{i=1}^{N_p} \delta I^\sigma_{iq}$ is {\tt RadV} & 
! $I_c$ or $I^s_c$ is {\tt Radiances} &
! $\frac{\partial I^\sigma_q}{\partial x_k}$ is {\tt K\_}$x${\tt\_FRQ} &
! $\frac{\partial I_c}{\partial x_k}$ is {\tt K\_}$x$
! \\
! \end{tabular}

        if ( FwdModelConf%anyLBL(sx) ) then
          call frequency_loop ( alpha_path_c(:npc), beta_path_c(:npc,:), c_inds,  &
            & del_s(:npc), del_zeta(:npc),do_calc_fzp(:npf,:),                    &
            & do_calc_zp(:npf,:), do_GL(:npc), eta_fzp(:npf,:),                   &
            & eta_zp(:npf,:), frequencies, h_path_c, incoptdepth(:npc),           &
            & p_path(:npf), pfaFalse, ref_corr(:npc), sps_path(:npf,:),           &
            & tau_lbl, t_path_c(:npc), t_script_lbl(:npc,:), tanh1_c(:npc),       &
            & tt_path_c(:s_i*npc), w0_path_c(:s_i*npc), z_path(:npf) )
          if ( print_TauL ) then
            call output ( thisSideband, before='Sideband ' )
            call output ( ptg_i, before=' Pointing ' )
            call dump ( tau_lbl, noFreqs, ' Tau_LBL:' )
            call dump ( t_script_lbl(:npc,:noFreqs), 'T_Script_LBL' )
          end if
        end if

        ! Handle PFA molecules
        if ( FwdModelConf%anyPFA(sx) ) then
          if ( frq_avg_sel == 15 ) then ! FRQ_avg + LBL + PFA + Derivs
            ! For every channel, frequency average the incremental radiance at
            ! every point along the path, giving Rad_Avg_Path for every channel
            ! and every point along the path.  Multiply by Tau_PFA to combine
            ! PFA contribution in Frequency_Loop.
            call frequency_avg_rad_path
            call frequency_average_derivatives ( .false. )
          end if
          call frequency_loop ( alpha_path_c(:npc), beta_path_c(:npc,:), c_inds, &
            & del_s(:npc), del_zeta(:npc), do_calc_fzp(:npf,:),                  &
            & do_calc_zp(:npf,:), do_GL(:npc), eta_fzp(:npf,:),                  &
            & eta_zp(:npf,:), channelCenters, h_path_c, incoptdepth(:npc),       &
            & p_path(:npf), pfaTrue, ref_corr(:npc), sps_path(:npf,:),           &
            & tau_pfa, t_path_c(:npc), t_script_pfa(:npc,:), tanh1_c(:npc),      &
            & tt_path_c(:s_i*npc), w0_path_c(:s_i*npc), z_path(:npf) )
          if ( print_TauP ) then
            call output ( thisSideband, before='Sideband ' )
            call output ( ptg_i, before=' Pointing ' )
            call dump ( tau_pfa, noUsedChannels, ' Tau_PFA:' )
            call dump ( t_script_pfa(:npc,:), 'T_Script_PFA' )
          end if

        end if

        call frequency_average ! or maybe just store

        ! If we're doing frequency averaging, there's a different frequency
        ! grid for each pointing, but we don't need to deallocate it here
        ! because the allocate_test in frequency_setup_2 will deallocate it.

        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call trace_end ( 'ForwardModel.Pointing ', index=ptg_i )

        ! End of pointing loop -------------------------------------------------
      end do

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_End ( 'ForwardModel.PointingLoop' )

      ! Frequency averaging and any LBL:
      if ( .not. associated(frequencies,channelCenters) ) &
        & call deallocate_test ( frequencies, 'frequencies', moduleName )
      
      call deallocate_test ( inc_rad_path,   'Inc_Rad_Path',   moduleName )

      call convolution ! or interpolation to ptan

      ! Deallocate maxNoPtgFreqs-sized stuff
      call deallocate_test ( radv, 'RadV', moduleName )
      if ( FwdModelConf%anyLBL(sx) ) then
        call deallocate_test ( t_script_LBL, 'T_Script_LBL', moduleName )
        call destroy_tau ( tau_LBL, "Tau_LBL", moduleName )
      end if

      call DestroyCompleteSlabs ( gl_slabs )
      if ( temp_der ) &
        & call deallocate_test ( k_temp_frq,   'k_temp_frq',     moduleName )

      if ( atmos_der ) &
        & call deallocate_test ( k_atmos_frq,  'k_atmos_frq',    moduleName )

      if ( spect_der ) then
        call deallocate_test ( k_spect_dw_frq, 'k_spect_dw_frq', moduleName )
        call deallocate_test ( k_spect_dn_frq, 'k_spect_dn_frq', moduleName )
        call deallocate_test ( k_spect_dv_frq, 'k_spect_dv_frq', moduleName )
      end if

      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call trace_end ( 'ForwardModel.Sideband ',index=thisSideband )

    end do            ! End of loop over sidebands -------------------------

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_End ( 'ForwardModel.SidebandLoop' )

    if ( FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2) ) &
      & call destroy_tau ( tau_PFA, "Tau_PFA", moduleName )

    !  **** DEBUG Printing cycle ...

! *** Create *seez* file for "nasty" purposes:

    if ( print_Seez ) call dump_print_code

! *** End of include

    call GetNameOfSignal ( firstSignal, sigName )
    i = index(sigName, '.B')
    j = index(sigName(i+2:), '.' )
    if ( j /= 0 ) sigName(i+j+1:) = ''

    if ( print_Rad ) then
      if ( FwdModelConf%do_conv ) then
        print *, 'Convolution: ON'
      else
        print *, 'Convolution: OFF'
      end if

      if ( FwdModelConf%do_freq_avg .and. any(fwdModelConf%anyLBL) ) then
        print *, 'Frequency Averaging: ON'
      else
        print *, 'Frequency Averaging: OFF'
!       print '(a,f12.4,a)', ' (All computations done at Frq =',Frequencies(1), ')'
      end if

      k = ptan%template%noSurfs
      print "( /'ptan\ ',i3.3)", k
      print "( 4(3x, f11.7) )", Ptan%values(1:k,maf)

      do i = 1, noUsedChannels
        print "(/, 'ch', i3.3, '_pfa_rad\ ', i3.3 )", channels(i)%used, k
        thisRadiance =>  &
          GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(channels(i)%signal)%index, sideband=sideband )
        j = thisRadiance%template%noChans
        channel = channels(i)%used - channels(i)%origin + 1
        print "( 4(2x, 1pg15.8) )", &
          & thisRadiance%values(channel:channel+j*(k-1):j, maf)
      end do
      print *

    end if

    ! **** End of Printing cycle ...

    call destroygrids_t ( grids_tscat )
    call destroygrids_t ( grids_salb )
    call destroygrids_t ( grids_cext )

    if ( FwdModelConf%incl_cld ) then
      call deallocate_test ( scat_src%values, 'scat_src%values', moduleName )
      call deallocate_test ( scat_alb%values, 'scat_alb%values', moduleName )
      call deallocate_test ( cld_ext%values,  'cld_ext%values',  moduleName )
    end if

    if ( toggle(emit) ) call trace_end ( 'Full ForwardModel MAF=', fmStat%maf )

  contains

  ! .......................................  Both_Sidebands_Setup  .....
    subroutine Both_Sidebands_Setup
    ! All of the setup stuff done for both sidebands.

      use Get_Chi_Out_m, only: Get_Chi_Out
      use Load_Sps_Data_m, only: Modify_values_for_supersat
      use Intrinsic, only: L_MIFDEADTIME, L_L1BMIF_TAI
      use ManipulateVectorQuantities, only: DoHGridsMatch, FindOneClosestInstance

      real(rp) :: D2XDXDT_SURFACE(1,sv_t_len) ! Would s_t*sv_t_len work?
      real(rp) :: D2XDXDT_TAN(ptan%template%nosurfs,sv_t_len) ! Would s_t*sv_t_len work?
      real(rp) :: One_dhdz(1), One_dxdh(1)
      real(rp) :: REQ_OUT(phitan%template%nosurfs)
      type (VectorValue_T), pointer :: WORK   ! Temporary stuff

      if ( toggle(emit)  .and. levels(emit) > 0 ) &
      &  call trace_begin ( 'ForwardModel.Both_Sidebands_Setup' )

      if ( e_stop > 0.0_rp ) e_stop = log(epsilon(0.0_rp)) ! only once

      fmStat%flags = 0 ! Assume no errors

      temp_der = present ( jacobian ) .and. FwdModelConf%temp_der
      atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der

      spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
      spect_der_center = spect_der .and. size(fwdModelConf%lineCenter) > 0
      spect_der_width = spect_der .and. size(fwdModelConf%lineWidth) > 0
      spect_der_width_TDep = spect_der .and. size(fwdModelConf%lineWidth_TDep) > 0
      spect_der = spect_der .and. &
        & ( spect_der_center .or. spect_der_width .or. spect_der_width_TDep )

      any_der = temp_der .or. atmos_der .or. spect_der

      ! Work out what we've been asked to do -----------------------------------

      beta_group => fwdModelConf%beta_group
      channels => fwdModelConf%channels
      DACsStaging => fwdModelConf%DACsStaging
      usedDACSSignals => fwdModelConf%usedDACSSignals
      noUsedDacs = size(usedDACSSignals)

      ! Identify the vector quantities we're going to need.
      ! The key is to identify the signal we'll be working with first
      firstSignal => fwdModelConf%signals(1) ! Config has verified that signals
        ! are all for same radiometer (actually LO), module and sideband
      sideband = merge ( 0, firstSignal%sideband, fwdModelConf%forceFoldedOutput )

      ! Start sorting out stuff from state vector ------------------------------

      ! Identify the appropriate state vector components
      ! VMRS are in beta_group%qty, gotten by get_species_data
      gph => GetVectorQuantityByType ( fwdModelExtra, quantityType=l_gph )
      earthRefl => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_earthRefl, config=fwdModelConf )
      losVel => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_losVel, instrumentModule=firstSignal%instrumentModule, &
        & config=fwdModelConf )
      orbIncline => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_orbitInclination, config=fwdModelConf )
      refGPH => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_refGPH, config=fwdModelConf )
      scGeocAlt => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_scGeocAlt )
      spaceRadiance => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_spaceRadiance, config=fwdModelConf )
      surfaceHeight => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_surfaceHeight, config=fwdModelConf, noError=.true. )
      if ( FwdModelConf%polarized ) then
        ECRtoFOV => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_ECRtoFOV, config=fwdModelConf )
      end if
      if ( FwdModelConf%incl_cld ) then
        cloudWater => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_cloudWater, noError=.true., config=fwdModelConf )
        sizeDistribution => GetQuantityForForwardModel( fwdModelIn, fwdModelExtra, &
          & quantityType=l_sizeDistribution, noError=.true., config=fwdModelConf )
      end if
      if ( FwdModelConf%i_saturation /= l_clear ) then
        boundaryPressure => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_boundaryPressure, config=fwdModelConf )
      end if

      ! Check that RefGPH and Temp have the same hGrid.  This is not checked in
      ! Construct or when the config is created.
      if ( .not. doHGridsMatch ( refGPH, temp ) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Different horizontal grids for refGPH and temperature' )

      ! Check that we have radiances for the channels that are used
      do sigInd = 1, size(fwdModelConf%signals)
        ! This just emits an error message and stops if we don't have a radiance.
        ! We don't use the vector quantity -- at least not right away.  We get
        ! it again later.
        thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, sideband=sideband )
      end do

      MAF = fmStat%maf

      Vel_Rel = losvel%values(1,maf) / speedOfLight ! Needed for PFA
      Vel_Cor = 1.0_rp - Vel_Rel

      windowStart = grids_tmp%windowStart(1)
      windowFinish = grids_tmp%windowFinish(1)

      ! Stuff for clouds

      if ( FwdModelConf%incl_cld ) then  ! Do this block only if incl_cld is true

        !??? Is thisRadiance appropriate here ???
        inst = FindOneClosestInstance ( temp, thisRadiance, MAF )

        ! checking done in ForwardModelSupport%ConstructForwardModelConfig
        nspec = no_mol ! Will be at least 3 if l_n2o is included, because
                       ! l_h2o and l_o3 are required
        noSurf  = temp%template%noSurfs
        vmrarray = 0.0_r8

        do j = 1, nspec      ! Loop over species

          if ( fwdModelConf%molecules(j) == l_h2o ) then
            ispec = 1
          else if ( fwdModelConf%molecules(j) == l_o3 ) then
            ispec = 2
          else if ( fwdModelConf%molecules(j) == l_n2o ) then
            ispec = 3
          else
            cycle
          end if

          vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,            &
            & quantityType=l_vmr, molecule=fwdModelConf%molecules(j) )

          novmrSurf = vmr%template%nosurfs

          call InterpolateValues ( &
          & reshape(vmr%template%surfs(:,1),(/novmrSurf/)), &    ! Old X
          & reshape(vmr%values(:,inst),(/novmrSurf/)),      &    ! Old Y
          & reshape(temp%template%surfs(:,1),(/noSurf/)),   &    ! New X
          & vmrArray(ispec,:),                              &    ! New Y
          & 'Linear', extrapolate='Clamp' )

        end do ! End of Loop over species

        ! This can be put outside the mmaf loop

        call allocate_test ( scat_src%values, n_t_zeta, fwdModelConf%num_scattering_angles, &
                             &'scat_src', moduleName )
        call allocate_test ( scat_alb%values, n_t_zeta, 2, 'scat_alb', moduleName )
        call allocate_test (  cld_ext%values, n_t_zeta, 2, 'cld_ext', moduleName )

      end if

  ! Set up our temporary `state vector' like arrays ------------------------

  ! modify h2o mixing ratio if a special supersaturation is requested
      if ( fwdModelConf%i_saturation /= l_clear ) &
        & call modify_values_for_supersat ( fwdModelConf, grids_f, h2o_ind, &
          & grids_tmp, boundaryPressure )

  ! set up output pointing angles ------------------------------------------
  ! note we have to compute req !!!!!!!

  !{ Compute equivalent earth radius $c$ at phi_t(1), nearest surface, where
  !  $c^2 = \frac{a^2\,b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta} =
  !         \frac{a^2}{\left(\frac{a^2}{b^2}-1\right) \sin^2 \beta + 1}$

      earthradc_sq = earthRadA**2 / &
            &     (Earth_Axis_Ratio_Squared_m1 * &
                   &   SIN(orbIncline%values(1,maf)*Deg2Rad)**2 + 1)
      earthradc = sqrt(earthradc_sq)

  !{\begin{equation*}\begin{split}
  ! R_{eq} =\;& \sqrt \frac{R_a^4 \sin^2 \phi + R_c^4 \cos^2 \phi}
  !                       {R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}\\
  !        =\;& \sqrt \frac{R_a^4 - (R_a^2+R_c^2)(R_a^2-R_c^2) \cos^2 \phi}
  !                       {R_c^2 +              (R_a^2-R_c^2) \cos^2 \phi}
  ! \end{split}\end{equation*}
  !
  ! Although this is the same mathematical formula as used in {\tt metrics},
  ! the $\phi$ used here is different.  Therefore, we can't use these
  ! values in {\tt metrics}, or where its output {\tt req} value is used.

      req_out = (earthrada-earthradc)*(earthrada+earthradc) * &
        & COS(phitan%values(:,maf)*Deg2Rad)**2
      ! Earthrad[abc] are in meters, but Req_Out needs to be in km.
      req_out = 0.001_rp * SQRT( &
        & ( earthrada**4 - (earthrada**2+earthradc_sq) * req_out ) / &
        & ( earthradc_sq + req_out ) )

  ! Compute reference Gauss Legendre (GL) grid ------------------------------

      ! Fill the "original" GL-grid arrays.  These depend only upon the
      ! original Zetas (Z_PSIG) derived from the L2CF using make_z_grid.
      ddhidhidtl0 => ddhidhidtl0_o
      dh_dt_glgrid => dh_dt_glgrid_o
      dhdz_glgrid => dhdz_glgrid_o
      h_glgrid => h_glgrid_o
      p_glgrid => p_glgrid_o
      t_glgrid => t_glgrid_o
      z_glgrid => z_glgrid_o

      call compute_GL_grid ( z_psig(:nlvl), p_glgrid, z_glgrid )

      ! estimate tan_phi and scgeocalt
      call estimate_tan_phi ( no_tan_hts, nlvl, maf, phitan, ptan, &
                            & scgeocalt, losvel, tan_press, tan_phi, est_scgeocalt, &
                            & est_los_vel )

   ! Now, allocate other variables we're going to need later ----------------

      if ( FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2) ) then
        call allocate_test ( tau_PFA%tau, max_c, noUsedChannels, 'Tau_PFA%Tau', &
          & moduleName )
        call allocate_test ( tau_PFA%i_stop, noUsedChannels, 'Tau_PFA%I_Stop', &
          & moduleName )
      end if

      ! This is only used for convolution, which is done for both sidebands.
      if ( fwdModelConf%scanAverage ) then
        work => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_l1bMIF_TAI )
        l1bMIF_TAI => work%values
        work => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_MIFDeadTime )
        MIFDeadTime => work%values ! Only the (1,1) element is used.
      else
        nullify ( l1bMIF_TAI, MIFDeadTime )
      end if

      ! Temperature's windowStart:windowFinish are correct here.
      ! RefGPH and Temperature have the same horizontal basis.
      ! Grids_F is only needed for H2O, for calculating refractive index.
      ! SCgeocAlt and RefGPH are in meters, but Get_Chi_Out wants them in km.
      ! This is only used for convolution, which is done for both sidebands.
      call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad, &
         & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp,                    &
         & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /),              &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish),             &
         & orbIncline%values(1,maf)*Deg2Rad, 0.0_rp,                       &
         & req_out, grids_f, h2o_ind, tan_chi_out, dh_dz_out, dx_dh_out,   &
         & dxdt_tan=dxdt_tan, d2xdxdt_tan=d2xdxdt_tan )

      ! This is a lazy way to get the surface angle
      ! Temperature's windowStart:windowFinish are correct here.
      ! refGPH and temperature have the same horizontal basis.
      ! Grids_tmp is only needed for H2O, for calculating refractive index.
      ! Est_scgeocalt and RefGPH are in meters, but Get_Chi_Out wants them in km.
      ! This is only used for convolution, which is done for both sidebands.
      call get_chi_out ( tan_press(1:1), tan_phi(1:1),                     &
         & 0.001_rp*est_scgeocalt(1:1), Grids_tmp,                         &
         & (/ (refGPH%template%surfs(1,1), j=windowStart,windowFinish) /), &
         & 0.001_rp*refGPH%values(1,windowStart:windowFinish),             &
         & orbIncline%values(1,maf)*Deg2Rad, 0.0_rp,                       &
         & req_out(1:1), grids_f, h2o_ind, surf_angle, one_dhdz, one_dxdh, &
         & dxdt_tan=dxdt_surface, d2xdxdt_tan=d2xdxdt_surface )

      if ( toggle(emit) .and. levels(emit) > 0 ) &
      &  call trace_end ( 'ForwardModel.Both_Sidebands_Setup' )

    end subroutine Both_Sidebands_Setup

  ! ................................................  Convolution  .....
    subroutine Convolution ! or simply interpolation to output grid

      ! Convolution if needed, or interpolation to ptan ----------------

      use Convolve_All_m, only: Convolve_Radiance, Convolve_Temperature_Deriv, &
        & Convolve_Other_Deriv, Interpolate_Radiance, &
        & Interpolate_Temperature_Deriv, Interpolate_Other_Deriv
      use FOV_Convolve_m, only: Convolve_Support_T, FOV_Convolve_Setup, &
        & FOV_Convolve_Teardown, No_FFT
      use MLSNumerics, only: Coefficients => Coefficients_r8, &
        & InterpolateArraySetup, InterpolateArrayTeardown

      logical, parameter :: OLDPATCHER = .false. ! Old one was buggy

      integer ChanInd, Channel, I, J, SigInd, Ptg_I
      real(rp) :: DELTAPTG         ! Used for patching the pointings
      integer :: MINSUPERSET       ! Min. value of superset > 0
      logical :: PATCHEDAPTG       ! Used in patching the pointings
      integer :: PTG_J, PTG_K      ! Loop counters for patching the pointings
      real(r8) :: RAD_FFT(s_t*no_fft) ! Convolved radiance on FFT grid
      integer :: SUPERSET          ! Output from AreSignalsSuperset
      real(rp) :: THISELEV         ! An elevation offset
      real(rp) :: THISFRACTION     ! A sideband fraction
      logical :: Update            ! Just update radiances etc.
      integer :: WHICHPATTERN      ! Index of antenna pattern
      type (Coefficients) :: Coeffs ! For interpolation
      type (Convolve_Support_T) :: Convolve_Support
      type (VectorValue_T), pointer :: ELEVOFFSET       ! Elevation offset
      type (VectorValue_T), pointer :: SIDEBANDFRACTION ! The sideband fraction to use
      type (VectorValue_T), pointer :: THISRADIANCE     ! A radiance vector quantity

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_begin ( 'ForwardModel.Convolution' )

      ! Check that the angles are in the correct order.  If they
      ! are not it means (give or take some approximations in the
      ! horizontal according to Bill), that the rays crossed over
      ! between the tangent point and the spacecraft.  One could dream
      ! up all sorts of elegant schemes to get around that problem, but
      ! it's simplest just to bail out (and is certainly preferable to
      ! the infinite loop in the convolution (Hunt on angles) that
      ! results otherwise).

      if ( print_Ptg ) &
        & call Dump ( ptg_angles, 'ptg_angles (before any patch)', format='(1PG22.17)' )
      
      ! This code is needed to ensure that the ptg_angles are monotonic
      ! (and not flat even)
      ptg_i = 2
      deltaPtg = 1e-3              ! Some starting value
      patchedAPtg = .false.
      do ptg_i = 2, no_tan_hts
        if ( ptg_angles(ptg_i) <= ptg_angles(ptg_i-1) .and. OLDPATCHER) then
          patchedAPtg = .true.
          ! This one is at or below its predecessor, find the next one above
          ptg_j = ptg_i + 2
          patchPtgInnerLoop: do 
            if ( ptg_j > no_tan_hts ) exit patchPtgInnerLoop
            if ( ptg_angles(ptg_j) > ptg_angles(ptg_i) ) exit patchPtgInnerLoop
            ptg_j = ptg_j + 1
          end do patchPtgInnerLoop
          ! Work out the spacing to fill in with
          if ( ptg_j > no_tan_hts ) then
            ! Fell off the end of the list, just use previous spacing
            ptg_j = no_tan_hts
          else
            ! Didn't fall off, so work out spacing
            deltaPtg = ( ptg_angles(ptg_j) - ptg_angles(ptg_i) ) / ( ptg_j - ptg_i )
          end if
          do ptg_k = ptg_i, ptg_j - 1
            ptg_angles(ptg_k) = ptg_angles(ptg_i-1) + ( ptg_k - ptg_i + 1 ) * deltaPtg
          end do
          ! Don't worry about missing the last one here, it will get caught by
          ! the next iteration of the outer loop
        else if ( ptg_angles(ptg_i) <= ptg_angles(ptg_i-1) ) then
          patchedAPtg = .true.
          ! This one is at or below its predecessor, find the next one above
          ptg_j = ptg_i + 1
          patchPtgInnerLoop2: do 
            if ( ptg_j > no_tan_hts ) exit patchPtgInnerLoop2
            if ( ptg_angles(ptg_j) > ptg_angles(ptg_i-1) ) exit patchPtgInnerLoop2
            ptg_j = ptg_j + 1
          end do patchPtgInnerLoop2
          ! Work out the spacing to fill in with
          if ( ptg_j > no_tan_hts ) then
            ! Fell off the end of the list, just use previous spacing
            ptg_j = no_tan_hts
          else
            ! Didn't fall off, so work out spacing
            deltaPtg = ( ptg_angles(ptg_j) - ptg_angles(ptg_i-1) ) / ( ptg_j - ptg_i + 1 )
          end if
          ptg_angles(ptg_i) = ptg_angles(ptg_i-1) + deltaPtg
        else
          ! This value us above the previous one so compute a delta from it 
          ! to use if needed later
          deltaPtg = ptg_angles(ptg_i) - ptg_angles(ptg_i-1)
        end if
      end do
      
      if ( patchedAPtg ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Had to patch some out-of-order ptg_angles' )
        fmStat%flags = ior(fmStat%flags,b_ptg_angles)
        if ( print_Ptg ) &
          & call Dump ( ptg_angles, 'ptg_angles (after patching)', format='(1PG22.17)' )
      end if

      ! Work out which antenna patterns we're going to need ------------------
      do i = 1, noUsedChannels
        channel = channels(i)%used
        chanInd = channel + 1 - channels(i)%origin
        sigInd = channels(i)%signal
        ! Get the radiance
        thisRadiance =>  &
          GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, sideband=sideband )
        ! Get the sideband fraction if we need to
        if ( firstSignal%sideband == 0 .or. fwdModelConf%forceSidebandFraction ) then
          sidebandFraction => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_limbSidebandFraction, &
            & signal=fwdModelConf%signals(sigInd)%index, &
            & sideband=thisSideband, config=fwdModelConf )
          thisFraction = sidebandFraction%values(chanInd,1)
        else                  ! Otherwise, want just unfolded signal
          thisFraction = 1.0
        end if
        ! Get the elevation offset
        elevOffset => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_elevOffset, signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=thisSideband, config=fwdModelConf )
        thisElev = elevOffset%values(chanInd,1) * deg2Rad

        ! Here comes the Convolution (or not) codes
        update = ( thisSideband /= fwdModelConf%sidebandStart )

        if ( FwdModelConf%do_conv ) then

          whichPattern = -1
          minSuperset = huge(0)
          do j = 1, size(antennaPatterns)
            superset = AreSignalsSuperset ( antennaPatterns(j)%signals, &
              & (/ fwdModelConf%signals(sigInd) /), sideband=thisSideband, &
              & channel=channel )
            if ( superset >= 0 .and. superset <= minSuperset ) then
              minSuperset = superset
              whichPattern = j
            end if
          end do
          if ( whichPattern < 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching antenna patterns." )

          call fov_convolve_setup ( antennaPatterns(whichPattern), ptg_angles, &
            & tan_chi_out-thisElev, convolve_support, &
            & do_dRad_dx=ptan_der, do_Scan_Avg=fwdModelConf%scanAverage )

          if ( temp_der ) then
            call convolve_radiance ( convolve_support, maf, chanInd, &
              & radiances(:,i), thisFraction, update, ptan, thisRadiance, &
              & L1BMIF_TAI, MIFDeadTime, &
              & Jacobian, fmStat%rows, dh_dz_out, dx_dh_out, ptan_der, rad_FFT )
            call convolve_temperature_deriv ( convolve_support, maf, chanInd, &
              & radiances(:,i), rad_fft, thisFraction, update, thisRadiance, &
              & temp, grids_tmp, surf_angle(1), L1BMIF_TAI, MIFDeadTime, &
              & real(k_temp(i,:,:),kind=rp), &
              & dx_dt, d2x_dxdt, dxdt_tan, dxdt_surface, Jacobian, fmStat%rows )
          else ! No temperature derivative
            call convolve_radiance ( convolve_support, maf, chanInd, &
              & radiances(:,i), thisFraction, update, ptan, thisRadiance, &
              & L1BMIF_TAI, MIFDeadTime, &
              & Jacobian, fmStat%rows, dh_dz_out, dx_dh_out, ptan_der )
          end if

          if ( atmos_der ) &
            & call convolve_other_deriv ( convolve_support, maf, chanInd, &
              & thisFraction, update, thisRadiance, beta_group%qty, Grids_f, &
              & L1BMIF_TAI, MIFDeadTime, real(k_atmos(i,:,:),kind=rp), &
              & Jacobian, fmStat%rows )

          if ( spect_der_center ) &
            & call convolve_other_deriv ( convolve_support, maf, chanInd, &
              & thisFraction, update, thisRadiance, &
              & fwdModelConf%lineCenter%qty, grids_v, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dv(i,:,:),kind=rp), Jacobian, fmStat%rows )
          if ( spect_der_Width ) &
            & call convolve_other_deriv ( convolve_support, maf, chanInd, &
              & thisFraction, update, thisRadiance, &
              & fwdModelConf%lineWidth%qty, grids_w, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dw(i,:,:),kind=rp), Jacobian, fmStat%rows )
          if ( spect_der_Width_TDep ) &
            & call convolve_other_deriv ( convolve_support, maf, chanInd, &
              & thisFraction, update, thisRadiance, &
              & fwdModelConf%lineWidth_TDep%qty, grids_n, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dn(i,:,:),kind=rp), Jacobian, fmStat%rows )

          call fov_convolve_teardown ( convolve_support )

        else          ! No convolution needed ..

          call interpolateArraySetup ( ptg_angles, tan_chi_out-thisElev, &
            & method='S', extrapolate='C', coeffs=coeffs, &
            & dyByDx=ptan_der.or.fwdModelConf%scanAverage )

          call interpolate_radiance ( coeffs, maf, chanInd, ptg_angles, &
            & radiances(:,i), thisFraction, update, ptan, tan_chi_out-thisElev, &
            & thisRadiance, L1BMIF_TAI, MIFDeadTime, Jacobian, fmStat%rows, &
            & dh_dz_out, dx_dh_out, ptan_der )

          if ( temp_der ) &
            & call interpolate_temperature_deriv ( coeffs, maf, chanInd, &
              & ptg_angles, thisFraction, update, tan_chi_out-thisElev,  &
              & thisRadiance, temp, grids_tmp, L1BMIF_TAI, MIFDeadTime, &
              & real(k_temp(i,:,:),kind=rp), Jacobian, fmStat%rows )

          if ( atmos_der ) &
            & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
              & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
              & beta_group%qty, grids_f, L1BMIF_TAI, MIFDeadTime, &
              & real(k_atmos(i,:,:),kind=rp), Jacobian, fmStat%rows, linear=.true. )
 
          if ( spect_der_center ) &
            & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
              & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
              & fwdModelConf%lineCenter%qty, grids_v, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dv(i,:,:),kind=rp), Jacobian, fmStat%rows, linear=.true. )

          if ( spect_der_Width ) &
            & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
              & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
              & fwdModelConf%lineWidth%qty, grids_w, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dw(i,:,:),kind=rp), Jacobian, fmStat%rows, linear=.true. )

          if ( spect_der_Width_TDep ) &
            & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
              & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
              & fwdModelConf%lineWidth_TDep%qty, grids_n, L1BMIF_TAI, MIFDeadTime, &
              & real(k_spect_dn(i,:,:),kind=rp), Jacobian, fmStat%rows, linear=.true. )

          call interpolateArrayTeardown ( coeffs )

        end if ! Convolve or interpolate

      end do                            ! Channel loop

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_end ( 'ForwardModel.Convolution' )

    end subroutine Convolution

  ! ............................................  Dump_Print_Code  .....
    subroutine Dump_Print_Code
      include "dump_print_code.f9h"
    end subroutine Dump_Print_Code

  ! ..........................................  Frequency_Average  .....
    subroutine Frequency_Average

      ! Here we either frequency average to get the unconvolved radiances, or
      ! we just store what we have as we're using monochromatic channels

      use Freq_Avg_m, only: Freq_Avg
      use SCRT_dn_m, only: SCRT_PFA

      integer :: C, ShapeInd

      if ( toggle(emit) .and. levels(emit) > 4 ) &
        & call trace_begin ( 'ForwardModel.FrequencyAvg' )

!          8 Frequency averaging?  N N Y Y - - Y Y - N
!          4 PFA?                  N N N N Y Y Y Y N Y
!          2 LBL?                  Y Y Y Y N N Y Y N Y
!          1 Derivatives?          N Y N Y N Y N Y - -
!            Frq_Avg_Sel value     2 3 A B C D E F
!
! =====================================================
!
! Impossible?                                      x x
! Radiances = RadV                 x x     x x
! K = K_frq                          x       x
! Frq Avg path integrated rad          x x     x
! Combine total path radiances                 x
! Frq Avg path integrated LBL deriv      x       x
! Frq Avg rad Along Path                         x
! Combine radiances along path                   x
! Combine LBL and PFA derivs                     x

      ! Do DACs stuff for all DACs channels first
      select case ( frq_avg_sel )
      case ( 10 : 15 )
        do c = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%filter%signal, &
            & fwdModelConf%signals(usedDacsSignals(c)), sideband = thisSideband )
          call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
            & RadV(:noFreqs), DACsStaging(:,c) )
        end do
      end select

      select case ( frq_avg_sel )
      case ( 2, 3, 12, 13 ) ! Just copy radiance. PFA + DACS - LBL impossible
        radiances(ptg_i,:) = radV(:noFreqs)
      case ( 10, 11, 14 )   ! Frq Avg path integrated radiance
        ! Now go through channel by channel
        do c = 1, noUsedChannels
          if ( channels(c)%dacs == 0 ) then
            shapeInd = channels(c)%shapeInds(sx)
            if ( fwdModelConf%anyPFA(sx) ) then ! frq_avg_sel = 14
            ! Combine LBL and PFA Tau's to get radiances.
            ! It's OK to combine the Tau's before doing the frequency
            ! averaging because the filter function is normalized.
              call SCRT_PFA (c, tau_LBL, tau_PFA, t_script_pfa, radV(:noFreqs) )
            end if
              call freq_Avg ( frequencies, &
                &   filterShapes(shapeInd)%filterGrid,  &
                &   filterShapes(shapeInd)%filterShape, &
                &   radV(:noFreqs), radiances(ptg_i,c) )
          else
            radiances(ptg_i,c) = DACsStaging(channels(c)%used,channels(c)%dacs)
          end if
        end do
      case ( 15 )   ! Frq Avg rad along path
        ! For every channel, we've frequency averaged the incremental radiance
        ! at every point along the path, giving Rad_Avg_Path for every channel
        ! and every point along the path.  Rad_Avg_Path was multiplied by
        ! Tau_PFA in Frequency_Loop to combine LBL and PFA contributions.
        do c = 1, noUsedChannels
          if ( channels(c)%dacs == 0 ) then
            radiances(ptg_i,c) = radV(c) ! Computed in Frequency_Loop
          else
            radiances(ptg_i,c) = DACsStaging(channels(c)%used,channels(c)%dacs)
          end if
        end do
      case default          ! Impossible cases: No LBL or PFA,
                            ! PFA and no frequency averaging
      end select

      if ( any_der ) call frequency_average_derivatives ( frq_avg_sel == 15 )

      if ( toggle(emit) .and. levels(emit) > 4 ) &
        & call trace_end ( 'ForwardModel.FrequencyAvg' )
    end subroutine Frequency_Average

  ! ...............................  Frequency_Average_Derivative  .....
    subroutine Frequency_Average_Derivative ( Grids, K_Frq, K, Mol, Combine )

      ! Frequency average or simply copy K_Frq to give K, the final
      ! Jacobian.

      type(grids_T), intent(in) :: Grids
      real(rp), intent(inout) :: K_FRQ(:,:) ! To be averaged  Frq X Grid
      real(r4), intent(out) :: K(:,:)       ! Averaged        Chan X Grid
      integer, intent(in) :: Mol            ! Which molecule
      logical, intent(in) :: Combine        ! "Combine LBL and PFA"

      integer :: C, ShapeInd
      real(rp) :: R    ! Frequency-averaged value
      integer :: SV_I  ! State-vector index

      if ( combine ) then
        ! Simply add newly-computed PFA derivatives in K_frq to
        ! previously-averaged LBL derivatives in K.  Remember that
        ! for PFA, the frequency dimension has extent noUsedChannels,
        ! not maxNoPtgFreqs.  The first dimension of K_frq is at least
        ! noUsedChannels, so we're guaranteed this will fit.
        do c = 1, noUsedChannels
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            k(c,sv_i) = k(c,sv_i) + k_frq(c,sv_i)
          end do
        end do
        return
      end if

      ! Only possible values for Frq_avg_sel here are 3, 11, 13, 15.
      ! See Frequency_Average for definition of Frq_avg_sel.

      select case ( frq_avg_sel )
      case ( 11, 15 ) ! See Frequency_Average.
        ! Do DACs stuff for all DACs channels first
        do c = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%filter%signal, &
            & fwdModelConf%signals(usedDacsSignals(c)), sideband = thisSideband )
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
              & k_frq(:,sv_i), DACsStaging2(:,sv_i,c) )
          end do                  ! Surface loop X Instance loop
        end do
        ! Now go through channel by channel
        do c = 1, noUsedChannels
          shapeInd = channels(c)%shapeInds(sx)
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            if ( grids%deriv_flags(sv_i) ) then
              if ( channels(c)%dacs == 0 ) then
                call Freq_Avg ( frequencies,            &
                  & FilterShapes(shapeInd)%FilterGrid,  &
                  & FilterShapes(shapeInd)%FilterShape, &
                  & k_frq(:,sv_i), r )
              else
                r = DACsStaging2 ( channels(c)%used, sv_i, channels(c)%dacs )
              end if
            else
              r = 0.0
            end if
            k(c,sv_i) = r
          end do                ! Grid loop
        end do                  ! Channel loop
      case ( 3, 13 )            ! Not frequency averaging, or PFA alone; copy.
        do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
          if ( grids%deriv_flags(sv_i) ) then
            k(:,sv_i) = min ( max ( k_frq(:,sv_i), &
                          &   real(-huge(0.0_r4), rp ) ), &
                          &   real( huge(0.0_r4), rp ) )
          else
            k(:,sv_i) = 0.0
          end if
        end do
      case default             ! Impossible
      end select               ! Frequency averaging or not
    end subroutine Frequency_Average_Derivative

  ! ..............................  Frequency_Average_Derivatives  .....
    subroutine Frequency_Average_Derivatives ( Combine )
      logical, intent(in) :: Combine        ! "Combine LBL and PFA"

      ! Frequency Average the temperature derivatives with the appropriate
      ! filter shapes

      integer :: UB ! Upper bound for first dimension of k_..._frq

      ub = noFreqs
      if ( combine ) ub = max(ub,noUsedChannels)

      if ( temp_der ) call frequency_average_derivative ( grids_tmp, &
        &               k_temp_frq(:ub,:), k_temp(:,ptg_i,:), 1, combine )

      ! Frequency Average the atmospheric derivatives with the appropriate
      ! filter shapes

      if ( atmos_der ) then
        do k = 1, no_mol
          if ( fwdModelConf%moleculeDerivatives(k) ) &
            & call frequency_average_derivative ( grids_f, &
              & k_atmos_frq(:ub,:), k_atmos(:,ptg_i,:), k, combine )
        end do                        ! Loop over major molecules
      end if                          ! Want derivatives for atmos

      ! Frequency Average the spectroscopic derivatives with the appropriate
      ! filter shapes

      do k = 1, size(fwdModelConf%lineCenter)
        call frequency_average_derivative &
          & ( grids_v, k_spect_dv_frq(:ub,:), k_spect_dv(:,ptg_i,:), k, &
          & combine )
      end do
      do k = 1, size(fwdModelConf%lineWidth)
        call frequency_average_derivative &
          & ( grids_w, k_spect_dw_frq(:ub,:), k_spect_dw(:,ptg_i,:), k, &
          & combine )
      end do
      do k = 1, size(fwdModelConf%lineWidth_TDep)
        call frequency_average_derivative &
          & ( grids_n, k_spect_dn_frq(:ub,:), k_spect_dn(:,ptg_i,:), k, &
          & combine )
      end do

    end subroutine Frequency_Average_Derivatives

  ! .....................................  Frequency_Avg_Rad_Path  .....
    subroutine Frequency_Avg_Rad_Path
      ! For every channel, frequency average the incremental radiance at
      ! every point along the path, giving Rad_Avg_Path for every channel
      ! and every point along the path.

      use Freq_Avg_m, only: Freq_Avg

      integer :: C, P, ShapeInd

      do c = 1, noUsedChannels
        shapeInd = channels(c)%shapeInds(sx)
        if ( channels(c)%dacs == 0 ) then
          do p = 1, npc
            call Freq_Avg ( frequencies,            &
                    & FilterShapes(shapeInd)%FilterGrid,  &
                    & FilterShapes(shapeInd)%FilterShape, &
                    & inc_rad_path(p,:), rad_avg_path(p,c) )
          end do
        end if
      end do
    end subroutine Frequency_Avg_Rad_Path

  ! .............................................  Frequency_Loop  .....
    subroutine Frequency_Loop ( Alpha_Path_c, Beta_Path_c, C_Inds, Del_S,    &
      & Del_Zeta, Do_Calc_fzp, Do_Calc_zp, Do_GL, Eta_fzp, Eta_zp,           &
      & Frequencies, H_Path_C, IncOptDepth, P_Path, PFA, Ref_Corr, Sps_Path, &
      & Tau, T_Path_c, T_Script, Tanh1_c, TT_Path_c, W0_Path_c, Z_Path )

      ! Having arguments instead of using host association serves two
      ! purposes:  The array sizes are implicit, so we don't need explicitly
      ! to mention them, and the pointer attribute gets stripped during the
      ! trip through the CALL statement -- hopefully thereby helping optimizers.

      use CS_Expmat_m, only: CS_Expmat
      use DO_T_SCRIPT_M, only: TWO_D_T_SCRIPT, TWO_D_T_SCRIPT_CLOUD
      use Get_Beta_Path_m, only: Get_Beta_Path, Get_Beta_Path_Cloud, &
        & Get_Beta_Path_PFA, Get_Beta_Path_Polarized
      use Get_d_Deltau_pol_m, only: Get_d_Deltau_pol_df, Get_d_Deltau_pol_dT
      use Mcrt_m, only: Mcrt_der
      use Opacity_m, only: Opacity
      use Path_Contrib_M, only: Path_Contrib
      use RAD_TRAN_M, only: RAD_TRAN_POL, DRAD_TRAN_DF, &
        & DRAD_TRAN_DT, DRAD_TRAN_DX
      use ScatSourceFunc, only: T_SCAT, Interp_Tscat, Convert_Grid

      real(rp), intent(out) :: Alpha_Path_c(:) ! Beta_Path * mixing ratio
      real(rp), intent(out) :: Beta_Path_c(:,:) ! path x species
      real(rp), intent(in) :: Del_S(:)    ! Integration lengths along path
      real(rp), intent(in) :: Del_Zeta(:) ! Integration lengths in Zeta coords
      integer, intent(in) :: C_Inds(:)    ! Selectors from complete path to coarse path
      logical, intent(inout) :: Do_Calc_fzp(:,:) ! 'Avoid zeros' indicator
      logical, intent(in) :: Do_Calc_zp(:,:) ! 'Avoid zeros' indicator
      logical, intent(out) :: Do_GL(:)    ! Where to do GL correction
      real(rp), intent(inout) :: Eta_fzp(:,:) ! path x (Eta_f x Eta_z x Eta_p)
      real(rp), intent(in) :: Eta_zp(:,:) ! path x (Eta_z x Eta_p)
      real(r8), intent(in) :: Frequencies(:)  ! The frequency grid
      real(rp), intent(in) :: H_Path_C(:) ! Heights on coarse path
      real(rp), intent(out) :: IncOptDepth(:)  ! Incremental optical depth
      real(rp), intent(in) :: P_Path(:)   ! Pressures along complete path
      logical, intent(in) :: PFA          ! Are we doing PFA or not?
      real(rp), intent(in) :: Ref_Corr(:) ! Refraction correction
      real(rp), intent(inout) :: Sps_Path(:,:) ! Species on path
      type(tau_t), intent(inout) :: Tau   ! Optical depth, inout so as not to
                                          ! undefine components' association status
      real(rp), intent(in) :: T_Path_c(:) ! Temperature on coarse path
      real(rp), intent(out) :: T_Script(:,:)! Delta_B in some notes
      real(rp), intent(out) :: TT_Path_c(:) ! tscat on coarse path
      real(rp), intent(out) :: Tanh1_c(:) ! tanh(frqhk/t_path_c)
      real(rp), intent(out) :: W0_Path_c(:) ! w0 on coarse path
      real(rp), intent(in) :: Z_Path(:)   ! -Log10(Pressures) along complete path

      integer :: FRQ_I                    ! Frequency loop index
      real(r8) :: FRQ                     ! Frequency
      real(r8) :: FRQHK                   ! 0.5 * Frq * H_Over_K
      integer :: I, J, L                  ! Loop inductor and subscript
      integer :: I_STOP                   ! Upper index for radiance comp.
      integer :: P_Stop                   ! Where to stop in polarized case
      logical :: PFA_or_not_pol           ! PFA .or. .not. fwdModelConf%polarized

      ! The following kludge is to work around a defect in the Fujitsu
      ! compiler: It appears to want to take a copy of inc_rad_path when
      ! a slice of it is used as an actual argument, thereby increasing
      ! run time by 10-20%.
      real(rp), pointer :: Inc_Rad_Path_Slice(:)

      ! Loop over frequencies ----------------------------------------------

      if ( toggle(emit) .and. levels(emit) > 4 ) &
        & call Trace_Begin ( 'ForwardModel.FrequencyLoop' )

      pfa_or_not_pol = pfa .or. .not. fwdModelConf%polarized

      do frq_i = 1, size(frequencies)

        if ( toggle(emit) .and. levels(emit) > 5 ) &
          & call Trace_Begin ( 'ForwardModel.Frequency ', index=frq_i )

        Frq = frequencies(frq_i)

        do_gl = .false.

        ! The following kludge is to work around a defect in the Fujitsu
        ! compiler: It appears to want to take a copy of inc_rad_path when
        ! a slice of it is used as an actual argument, thereby increasing
        ! run time by 10-20%.
        inc_rad_path_slice => inc_rad_path(:npc,frq_i)

        ! Set up path quantities --------------------------------------

        ! Compute the sps_path for this Frequency
        call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
          & Frq, eta_zp, do_calc_zp, sps_path, do_calc_fzp, eta_fzp )

        if ( pfa ) then
          call get_beta_path_PFA ( frq, frq_i, z_path, c_inds, t_path_c,  &
            & beta_group, sx, vel_rel, beta_path_c, t_der_path_flags,     &
            & dbeta_dT_path_c, dbeta_dw_path_c, dbeta_dn_path_c, dbeta_dv_path_c )
        else
          frqhk = 0.5_r8 * frq * h_over_k    ! h nu / 2 k
          tanh1_c = tanh( frqhk / t_path_c ) ! tanh ( h nu / 2 k T )
          ! dTanh_dT = -h nu / (2 k T**2) 1/tanh1 d(tanh1)/dT
          if ( temp_der ) dTanh_dT_c(:npc) = &
              & frqhk / t_path_c**2 * ( tanh1_c - 1.0_rp / tanh1_c )
          call get_beta_path ( Frq, p_path, t_path_c, tanh1_c,                &
            &  beta_group, sx, fwdModelConf%polarized, gl_slabs, c_inds,      &
            &  beta_path_c, t_der_path_flags, dTanh_dT_c, vel_rel,            &
            &  dbeta_dT_path_c, dbeta_dw_path_c, dbeta_dn_path_c, dbeta_dv_path_c )
        end if

        if ( FwdModelConf%incl_cld .and. .not. pfa ) then
          ! Compute Scattering source function based on temp prof at all
          ! angles U for each temperature layer assuming a plane parallel
          ! atmosphere.

          if ( ptg_i == 1 ) then
          ! ??? Can this work?  On all pointings after the first one, ???
          ! ??? it uses the result for the last frequency.            ???
          ! ??? We have a different frequency grid for each pointing, ???
          ! ??? so is this the wrong idea in the first place?         ???

            call T_scat ( temp%values(:,inst), Frq, GPH%values(:,inst),    & 
               & 10.0**(-temp%template%surfs), vmrArray, nspec,            &
               & fwdModelConf%num_scattering_angles,                       &
               & fwdModelConf%num_azimuth_angles,                          &
               & fwdModelConf%num_ab_terms, fwdModelConf%num_size_bins,    &
               & fwdModelConf%no_cloud_species,                            &
               & scat_src%values, scat_alb%values, cld_ext%values, Scat_ang)

          end if

          scat_src%template = temp%template

          call load_one_item_grid ( grids_tscat, scat_src, phitan, maf, &
            & fwdModelConf, .false., .true. )

          call allocate_test ( do_calc_tscat, npf, size(grids_tscat%values),           &
                             & 'do_calc_tscat', moduleName )
          call allocate_test ( do_calc_tscat_zp, npf, grids_tscat%p_len,               &
                             & 'do_calc_tscat_zp', moduleName )
          call allocate_test ( eta_tscat,     npf, size(grids_tscat%values),           &
                             & 'eta_tscat',     moduleName )
          call allocate_test ( eta_tscat_zp,  npf, grids_tscat%p_len,                  &
                             & 'eta_tscat_zp',  moduleName )
          call allocate_test ( tscat_path,    npf, fwdModelConf%num_scattering_angles, &
                               & 'tscat_path',  moduleName )

          call comp_eta_docalc_no_frq ( Grids_Tscat, z_path(1:npf), &
            &  phi_path(1:npf), eta_tscat_zp(1:npf,:),              &
            &  do_calc_tscat_zp(1:npf,:) )

        ! Frq=0.0
          call comp_sps_path_frq ( Grids_tscat, firstSignal%lo, thisSideband, &
            & 0.0_r8, eta_tscat_zp(1:npf,:),                               &
            & do_calc_tscat_zp(1:npf,:), tscat_path(1:npf,:),              &
            & do_calc_tscat(1:npf,:), eta_tscat(1:npf,:) )

          ! project Tscat onto LOS
          call interp_tscat ( tscat_path(1:npf,:), Scat_ang(:), &
            & phi_path(1:npf), tt_path(1:npf,:) )

          if ( .not. cld_fine ) then                 ! interpolate onto gl grids along the LOS

            scat_alb%template = temp%template
            cld_ext%template  = temp%template

            call load_one_item_grid ( grids_salb,  scat_alb, phitan, maf, fwdModelConf, .false.)
            call load_one_item_grid ( grids_cext,  cld_ext,  phitan, maf, fwdModelConf, .false.)             

            do i = 1, size(grids_salb%values)
               if ( abs(grids_salb%values(i)) < TOL) then
                       grids_salb%values(i) = 0.0
               end if
               if ( abs(grids_cext%values(i)) < TOL) then
                       grids_cext%values(i) = 0.0
               end if
            end do

            call allocate_test (do_calc_salb, npf, size(grids_salb%values),'do_calc_salb',moduleName)
            call allocate_test (do_calc_salb_zp, npf, grids_salb%p_len,               &
                               & 'do_calc_salb_zp', moduleName )
            call allocate_test (eta_salb,    npf, size(grids_salb%values), 'eta_salb', moduleName)
            call allocate_test (eta_salb_zp, npf, grids_salb%p_len, 'eta_salb_zp', moduleName)
            call allocate_test ( salb_path,  npf, 1, 'salb_path', moduleName )

            call allocate_test (do_calc_cext,npf,size(grids_cext%values),'do_calc_cext',moduleName)
            call allocate_test (do_calc_cext_zp, npf, grids_cext%p_len,               &
                               & 'do_calc_cext_zp', moduleName )
            call allocate_test (eta_cext,    npf, size(grids_cext%values), 'eta_cext', moduleName)
            call allocate_test (eta_cext_zp, npf, grids_cext%p_len,  'eta_cext_zp', moduleName)
            call allocate_test (cext_path,   npf, 1, 'cext_path', moduleName)

            call comp_eta_docalc_no_frq ( Grids_salb, z_path(1:npf), &
              &  phi_path(1:npf), eta_salb_zp(1:npf,:),              &
              &  do_calc_salb_zp(1:npf,:) )

          ! Frq=0.0
            call comp_sps_path_frq ( Grids_salb, firstSignal%lo, thisSideband, &
              & 0.0_r8, eta_salb_zp(1:npf,:),                                  &
              & do_calc_salb_zp(1:npf,:), salb_path(1:npf,:),                  &
              & do_calc_salb(1:npf,:), eta_salb(1:npf,:) )

            call comp_eta_docalc_no_frq ( Grids_cext, z_path(1:npf), &
              &  phi_path(1:npf), eta_cext_zp(1:npf,:),              &
              &  do_calc_cext_zp(1:npf,:) )

          ! Frq=0.0
            call comp_sps_path_frq ( Grids_cext, firstSignal%lo, thisSideband, &
              & 0.0_r8, eta_cext_zp(1:npf,:),                                  &
              & do_calc_cext_zp(1:npf,:), cext_path(1:npf,:),                  &
              & do_calc_cext(1:npf,:), eta_cext(1:npf,:) )

            call convert_grid ( salb_path(1:npf,:), cext_path(1:npf,:),  & 
                              & tt_path(1:npf,:), c_inds,                & 
                              & beta_path_cloud_c(1:npc), w0_path_c,           & 
                              & tt_path_c )

          else ! Not cld_fine              re-compute cext and w0 along the LOS

            call get_beta_path_cloud ( Frq, p_path, t_path(1:npf),      &
              &  tt_path(1:npf,:), c_inds,                              &
              &  beta_path_cloud_c(1:npc), w0_path_c,  tt_path_c,       &
              &  IPSD(1:npf),  WC(:,1:npf), fwdModelConf )

          end if

          call deallocate_test ( do_calc_tscat,    'do_calc_tscat',    moduleName )
          call deallocate_test ( do_calc_tscat_zp, 'do_calc_tscat_zp', moduleName )
          call deallocate_test ( eta_tscat,        'eta_tscat',        moduleName )
          call deallocate_test ( eta_tscat_zp,     'eta_tscat_zp',     moduleName )           
          call deallocate_test ( tscat_path,       'tscat_path',       moduleName )

          call deallocate_test ( do_calc_salb,     'do_calc_salb',     moduleName )
          call deallocate_test ( do_calc_salb_zp,  'do_calc_salb_zp',  moduleName )
          call deallocate_test ( eta_salb,         'eta_salb',         moduleName )
          call deallocate_test ( eta_salb_zp,      'eta_salb_zp',      moduleName )
          call deallocate_test ( salb_path,        'salb_path',        moduleName )

          call deallocate_test ( do_calc_cext,     'do_calc_salb',     moduleName )
          call deallocate_test ( do_calc_cext_zp,  'do_calc_salb_zp',  moduleName )
          call deallocate_test ( eta_cext,         'eta_cext',         moduleName )
          call deallocate_test ( eta_cext_zp,      'eta_cext_zp',      moduleName )
          call deallocate_test ( cext_path,        'cext_path',        moduleName )

          do j = 1, npc ! Don't trust compilers to fuse loops
            alpha_path_c(j) = dot_product( sps_path(c_inds(j),:), &
                                  &        beta_path_c(j,:) )     &
                                  & + beta_path_cloud_c(j)
            incoptdepth(j) = alpha_path_c(j) * del_s(j)
          end do

          ! Needed to compute inc_rad_path_slice and by rad_tran_pol
          call two_d_t_script_cloud ( t_path_c, tt_path_c, w0_path_c, &  
            & spaceRadiance%values(1,1), frq, t_script(:,frq_i), B(:npc) )

        else ! Not cloud model

          !{ {\tt incoptdepth} is $\Delta \delta_{s\rightarrow s-1} =
          !  \int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
          !    \frac{\text{d}h}{\text{d}\zeta} \text{d} \zeta$ where
          !  $G(\zeta)$ is {\tt alpha\_path\_c}, which is approximated
          !  here by the rectangle rule, \emph{viz.}
          !  $\Delta \delta_{i\rightarrow i-1} \approx G(\zeta_i) \delta s_i$.

          do j = 1, npc ! Don't trust compilers to fuse loops
            alpha_path_c(j) = dot_product( sps_path(c_inds(j),:), &
                                         & beta_path_c(j,:) )
            incoptdepth(j) = alpha_path_c(j) * del_s(j)
          end do

          ! Needed to compute inc_rad_path_slice and by rad_tran_pol
          call two_d_t_script ( t_path_c, &  
            & spaceRadiance%values(1,1), frq, t_script(:,frq_i), B(:npc) )

        end if ! end of check cld 

        if ( .not. fwdModelConf%polarized ) then
          ! Determine where to use Gauss-Legendre for scalar instead of a trapezoid.

          call path_contrib ( incoptdepth, tan_pt_c, e_rflty, fwdModelConf%tolerance, &
            &                 do_gl )

        else ! extra stuff for polarized case

          call get_beta_path_polarized ( frq, h, beta_group%lbl(sx), gl_slabs, &
            & c_inds, beta_path_polarized, dBeta_dT_polarized_path_c )

          ! We put an explicit extent of -1:1 for the first dimension in
          ! the hope a clever compiler will do better optimization with
          ! a constant extent.
          ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
          ! to alpha here

          do j = 1, npc
            alpha_path_polarized(-1:1,j) = &
              & matmul( beta_path_polarized(-1:1,j,:),       &
              &         sps_path(c_inds(j),:) ) * tanh1_c(j) &
              & + 0.25 * alpha_path_c(j)
            alpha_path_polarized(0,j) = alpha_path_polarized(0,j) + &
              & 0.25 * alpha_path_c(j)
          end do

          ! Turn sigma-, pi, sigma+ into 2X2 matrix incoptdepth_pol
          call opacity ( ct(1:npc), stcp(1:npc), stsp(1:npc), &
            & alpha_path_polarized(:,1:npc), incoptdepth_pol(:,:,1:npc) )

          ! We don't add unpolarized incremental optical depth to diagonal
          ! of polarized incremental optical depth because we added the
          ! scalar alpha_path to the sigma-, pi and sigma+ parts of
          ! alpha_path_polarized above.  If we did add it here, we would
          ! need 0.5 factors to scale unpolarized "power absorption" to
          ! get "field absorption"

          do j = 2, npc-1
            incoptdepth_pol(1,1,j) = - incoptdepth_pol(1,1,j) * del_s(j)
            incoptdepth_pol(2,1,j) = - incoptdepth_pol(2,1,j) * del_s(j)
            incoptdepth_pol(1,2,j) = - incoptdepth_pol(1,2,j) * del_s(j)
            incoptdepth_pol(2,2,j) = - incoptdepth_pol(2,2,j) * del_s(j)

            ! deltau_pol = exp(incoptdepth_pol)
            call cs_expmat ( incoptdepth_pol(:,:,j), deltau_pol(:,:,j) )
          end do

          ! Determine where to do GL
          call path_contrib ( deltau_pol(:,:,1:npc), tan_pt_c, e_rflty, &
             & fwdModelConf%tolerance, do_gl )

        end if

        !{ We want $\Delta \delta_{s\rightarrow s-1} =
        ! \int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
        !   \frac{\text{d}h}{\text{d}\zeta} \text{d} \zeta$ where
        ! $G(\zeta)$ is {\tt alpha\_path\_c}, but $\frac{\text{d}s}{\text{d}h}$
        ! is singular at the tangent point, so we compute
        ! $\Delta \delta_{s\rightarrow s-1} = G(\zeta_i) \delta s_i +
        ! \int_{\zeta_i}^{\zeta_{i-1}} \left(G(\zeta)-G(\zeta_i)\right)
        ! \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
        ! \text{d} \zeta$.
        !
        ! Where we do GL, the second integral is approximated using GL (in
        ! {\tt Get\_Tau}).  Where we don't do GL, approximate it using the
        ! trapezoid rule (here).

        do j = 2, tan_pt_c
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = incoptdepth(j) + &
              & 0.5 * ( alpha_path_c(j-1) - alpha_path_c(j) ) * dsdz_c(j-1)*del_zeta(j)
        end do
        do j = tan_pt_c+1, npc-1
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = incoptdepth(j) + &
              & 0.5 * ( alpha_path_c(j+1) - alpha_path_c(j) ) * dsdz_c(j+1)*del_zeta(j)
        end do

        call get_GL_inds ( do_gl, tan_pt_c, gl_inds_b, cg_inds, ngl, ncg )
        gl_inds => gl_inds_b(:ngl)
        ! ngl is ng * count(do_gl)
        t_path_f(:ngl) = t_path(gl_inds)

        if ( pfa ) then

          call get_beta_path_PFA ( frq, frq_i, z_path, gl_inds, t_path_f(:ngl),   &
            & beta_group, sx, vel_rel, beta_path_f(:ngl,:), t_der_path_flags,     &
            & dbeta_dT_path_f, dbeta_dw_path_f, dbeta_dn_path_f, dbeta_dv_path_f )

        else

          tanh1_f(1:ngl) = tanh( frqhk / t_path_f(:ngl) )
          ! dTanh_dT = -h nu / (2 k T**2) 1/tanh1 d(tanh1)/dT
          if ( temp_der ) &
            & dTanh_dT_f(1:ngl) = frqhk / t_path_f(1:ngl)**2 * &
              & ( tanh1_f(1:ngl) - 1.0_rp / tanh1_f(1:ngl) )

          ! The derivatives that get_beta_path computes depend on which
          ! derivative arrays are allocated, not which ones are present.
          ! This avoids having four paths through the code, each with a
          ! different set of optional arguments.

          call get_beta_path ( Frq, p_path, t_path_f(:ngl), tanh1_f(1:ngl),    &
            & beta_group, sx, fwdModelConf%polarized, gl_slabs,                &
            & gl_inds, beta_path_f(:ngl,:), t_der_path_flags, dTanh_dT_f,      &
            & vel_rel, &
            & dbeta_dT_path_f, dbeta_dw_path_f, dbeta_dn_path_f, dbeta_dv_path_f )

        end if ! .not. pfa

        do j = 1, ngl ! loop around dot_product instead of doing sum(a*b,2)
                      ! to avoid path-length array temps
          alpha_path_f(j) = dot_product( sps_path(gl_inds(j),:),  &
                                       & beta_path_f(j,:) )
        end do

        if ( .not. fwdModelConf%polarized ) then

        ! Compute SCALAR radiative transfer --------------------------

          call get_tau ( frq_i, gl_inds, cg_inds(1:ncg), e_rflty, del_zeta, &
            & alpha_path_c, ref_corr, incoptdepth, tan_pt_c,                &
            & alpha_path_f(1:ngl), dsdz_gw_path, tau )
            i_stop = tau%i_stop(frq_i)

          ! Get incremental radiance and radiance from Tau and T_Script
          ! Don't clobber them if doing PFA and already did LBL.
          if ( .not. pfa .or. .not. fwdModelConf%anyLBL(sx) ) then
            radV(frq_i) = 0.0
            do j = 1, i_stop
              inc_rad_path_slice(j) = t_script(j,frq_i) * tau%tau(j,frq_i)
              radV(frq_i) = radV(frq_i) + inc_rad_path_slice(j)
            end do ! j
            inc_rad_path_slice(i_stop+1:) = 0
          end if

          if ( pfa .and. frq_avg_sel == 15 ) then ! See Frequency_Average.
            ! Multiply Rad_Avg_Path by Tau to combine LBL and PFA.  Then
            ! sum to give RadV.  Remember, if PFA, Frq_I is a channel number.
            radV(frq_i) = 0.0
            do j = 1, i_stop
              Rad_Avg_Path(j,frq_i) = Rad_Avg_Path(j,frq_i) * tau%tau(j,frq_i)
              radV(frq_i) = radV(frq_i) + Rad_Avg_Path(j,frq_i)
            end do ! j
            inc_rad_path_slice => Rad_Avg_Path(:npc,frq_i)
          end if
        else ! Polarized model; can't combine with PFA

        ! Compute POLARIZED radiative transfer -----------------------------

          i_stop = npc ! needed by drad_tran_df

          ! Get the corrections to integrals for layers that need GL for
          ! the polarized species.
          call get_beta_path_polarized ( frq, h, beta_group%lbl(sx), gl_slabs, &
            & gl_inds, beta_path_polarized_f, dBeta_dT_polarized_path_f )

          ! We put an explicit extent of -1:1 for the first dimension in
          ! the hope a clever compiler will do better optimization with
          ! a constant extent.
          ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
          ! to alpha here.
          do j = 1, ngl
            alpha_path_polarized_f(-1:1,j) = &
              & matmul( beta_path_polarized_f(-1:1,j,:), &
              &         sps_path(gl_inds(j),:) ) * tanh1_f(j) &
              & + 0.25 * alpha_path_f(j)
            alpha_path_polarized_f(0,j) = alpha_path_polarized_f(0,j) +               &
              & 0.25 * alpha_path_f(j)
          end do

          call rad_tran_pol ( tan_pt_c, gl_inds, cg_inds(1:ncg), e_rflty, del_zeta, &
            & alpha_path_polarized(:,1:npc), ref_corr, incoptdepth_pol(:,:,1:npc),  &
            & deltau_pol(:,:,1:npc), alpha_path_polarized_f(:,1:ngl), dsdz_gw_path, &
            & ct, stcp, stsp, t_script(:,frq_i), prod_pol(:,:,1:npc),               &
            & tau_pol(:,:,1:npc), rad_pol, p_stop )

          if ( p_stop < 0 ) then ! exp(incoptdepth_pol(:,:,-p_stop)) failed
            call output ( 'Exp(incoptdepth_pol(:,:,' )
            call output ( -p_stop )
            call output ( ') failed.  Value is', advance='yes' )
            call dump ( incoptdepth_pol(:,:,-p_stop), clean=.true. )
            call output ( 'thisSideband = ' ); call output ( thisSideband )
            call output ( ', ptg_i = ' ); call output ( ptg_i )
            call output ( ', frq_i = ' ); call output ( frq_i, advance='true' )
            call MLSMessage ( MLSMSG_Error, moduleName, &
              & 'exp(incoptdepth_pol) failed' )
          end if

          if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
            RadV(frq_i) = real(rad_pol(1,1))
          else
            RadV(frq_i) = real(rad_pol(2,2))
          end if

        end if

        ! Compute derivatives if needed ----------------------------------

        if ( atmos_der ) then

          call drad_tran_df ( c_inds, gl_inds, del_zeta, Grids_f, &
            &  beta_path_c, eta_fzp, sps_path, do_calc_fzp,       &
            &  beta_path_f, do_gl, del_s, ref_corr, dsdz_gw_path, &
            &  inc_rad_path_slice, tan_pt_c, i_stop,              &
            &  d_delta_df(1:npc,:), k_atmos_frq(frq_i,:) )
          if ( FwdModelConf%anyPFA(sx) ) then

          else if ( fwdModelConf%polarized ) then

            ! VMR derivatives for polarized radiance.
            ! Compute DE / Df from D Incoptdepth_pol / Df and put
            ! into DE_DF.
            call Get_D_Deltau_Pol_DF ( ct, stcp, stsp, c_inds(1:p_stop),   &
              &  del_zeta, Grids_f, beta_path_polarized(:,1:p_stop,:),     &
              &  tanh1_c(:npc), eta_fzp, do_calc_fzp, sps_path, del_s,     &
              &  incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),        &
              &  d_delta_df(1:npc,:), de_df(:,:,1:p_stop,:) )

            ! Compute D radiance / Df from Tau, Prod, T_Script
            ! and DE / Df.
            call mcrt_der ( t_script(:,frq_i), sqrt(e_rflty),    &
              & deltau_pol(:,:,1:npc), de_df(:,:,1:npc,:),       &
              & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
              & tan_pt_c, d_rad_pol_df )

            if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
              k_atmos_frq(frq_i,:) = real(d_rad_pol_df(1,1,:))
            else
              k_atmos_frq(frq_i,:) = real(d_rad_pol_df(2,2,:))
            end if

          end if ! PFA

        end if

        if ( temp_der ) then

          ! get d Delta B / d T * d T / eta
          call dt_script_dt ( t_path_c, B(:npc), eta_zxp_t_c(1:npc,:), &
                            & frq, d_t_scr_dt(1:npc,:) )

          dh_dt_path_f(:ngl,:) = dh_dt_path(gl_inds,:)
          do_calc_t_f(:ngl,:) = do_calc_t(gl_inds,:)
          eta_zxp_t_f(:ngl,:) = eta_zxp_t(gl_inds,:)
          h_path_f(:ngl) = h_path(gl_inds)
          dAlpha_dT_path_c(:npc) = sum( sps_path(c_inds,:) *  &
                                        dBeta_dT_path_c(1:npc,:),dim=2 )
          dAlpha_dT_path_f(:ngl) = sum( sps_path(gl_inds,:) * &
                                        dBeta_dT_path_f(1:ngl,:),dim=2 )

          if ( pfa_or_not_pol ) then

            call drad_tran_dt ( del_zeta, h_path_c,                          &
              & dh_dt_path_c(1:npc,:), alpha_path_c,                         &
              & dAlpha_dT_path_c(:npc), eta_zxp_t_c(1:npc,:),                &
              & do_calc_t_c(1:npc,:), do_calc_hyd_c(1:npc,:), del_s,         &
              & ref_corr, tan_ht, dh_dt_path(tan_pt_f,:),                    &
              & do_gl, gl_inds, h_path_f(:ngl), t_path_f(:ngl),              &
              & dh_dt_path_f(:ngl,:), alpha_path_f(1:ngl),                   &
              & dAlpha_dT_path_f(:ngl), eta_zxp_t_f(:ngl,:),                 &
              & do_calc_t_f(:ngl,:), path_dsdh, dhdz_gw_path, dsdz_gw_path,  &
              & d_t_scr_dt(1:npc,:), tau%tau(:npc,frq_i),                    &
              & inc_rad_path_slice, tan_pt_c, i_stop, grids_tmp%deriv_flags, &
              & pfa .and. frq_avg_sel == 15, k_temp_frq(frq_i,:) )

          else ! pol and not pfa

            ! Temperature derivatives for polarized radiance
            ! Compute DE / DT from D Incoptdepth_Pol / DT and put
            ! into DE_DT.

            dAlpha_dT_polarized_path_c(:,1:npc) = 0.0
            dAlpha_dT_polarized_path_f(:,1:ngl) = 0.0
            do j = 1, no_mol
              do l = -1, 1
                dAlpha_dT_polarized_path_c(l,1:npc) = &
              & dAlpha_dT_polarized_path_c(l,1:npc) + &
                  & sps_path(c_inds,j) * dBeta_dT_polarized_path_c(l,1:npc,j)
                dAlpha_dT_polarized_path_f(l,1:ngl) = &
              & dAlpha_dT_polarized_path_f(l,1:ngl) + &
                  & sps_path(gl_inds,j) * dBeta_dT_polarized_path_f(l,1:ngl,j)
              end do ! l
            end do ! j
            do l = -1, 1
              dAlpha_dT_polarized_path_c(l,1:npc) =               &
                & dAlpha_dT_polarized_path_c(l,1:npc) * tanh1_c + &
                & alpha_path_polarized(l,:npc) * dTanh_dT_c(:npc)
              dAlpha_dT_polarized_path_f(l,1:ngl) =                     &
                & dAlpha_dT_polarized_path_f(l,1:ngl) * tanh1_f(:ngl) + &
                & alpha_path_polarized_f(l,:ngl) * dTanh_dT_f(:ngl)
            end do

            call get_d_deltau_pol_dT ( ct, stcp, stsp, tan_pt_c,            &
              & t_path_f(:ngl), alpha_path_polarized(:,1:p_stop),           &
              & alpha_path_polarized_f(:,1:ngl),                            &
              & dAlpha_dT_path_c(:npc), dAlpha_dT_path_f(:ngl),             &
              & dAlpha_dT_polarized_path_c, dAlpha_dT_polarized_path_f,     &
              & eta_zxp_t_c(1:p_stop,:), eta_zxp_t_f(:ngl,:), del_s,        &
              & gl_inds, del_zeta, do_calc_t_c(1:p_stop,:),                 &
              & do_calc_t_f(:ngl,:), do_gl(1:p_stop), path_dsdh,            &
              & dhdz_gw_path, dsdz_gw_path, incoptdepth_pol(:,:,1:p_stop),  &
              & ref_corr(1:p_stop), h_path_c, h_path_f(:ngl),               &
              & dh_dt_path_c(1:p_stop,:),dh_dt_path_f(:ngl,:),              &
              & tan_ht, dh_dt_path(tan_pt_f,:),                             &
              & do_calc_hyd_c(1:p_stop,:), grids_tmp%deriv_flags,           &
              & de_dt(:,:,1:p_stop,:) )

            ! Compute D radiance / DT from Tau, Prod, T_Script, D_T_Scr_dT
            ! and DE / DT.

            call mcrt_der ( t_script(:,frq_i), sqrt(e_rflty),    &      
              & deltau_pol(:,:,1:npc), de_dt(:,:,1:npc,:),       &      
              & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &      
              & tan_pt_c, d_rad_pol_dt, d_t_scr_dt(1:npc,:) )

            if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
              k_temp_frq(frq_i,:) = real(d_rad_pol_dt(1,1,:))
            else
              k_temp_frq(frq_i,:) = real(d_rad_pol_dt(2,2,:))
            end if

          end if

        end if

        if ( spect_der ) then

          if ( fwdModelConf%polarized ) then
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "Spectroscopic derivatives for Polarized species not implemented yet." )
          else if ( pfa ) then
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "Spectroscopic derivatives for PFA not implemented yet." )
          else

            ! Spectroscopic derivative  wrt: W

            if ( spect_der_width ) &
              & call drad_tran_dx ( c_inds, gl_inds, del_zeta, grids_w,      &
                &  eta_zxp_w, sps_path, fwdModelConf%lineWidth%beta(sx),     &
                &  do_calc_w, dbeta_dw_path_c, dbeta_dw_path_f, do_gl, del_s,&
                &  ref_corr, dhdz_gw_path, inc_rad_path_slice, tan_pt_c,     &
                &  i_stop, k_spect_dw_frq(frq_i,:) )

            ! Spectroscopic derivative  wrt: N

            if ( spect_der_width_TDep ) &
              & call drad_tran_dx ( c_inds, gl_inds, del_zeta, grids_n,      &
                &  eta_zxp_n, sps_path, fwdModelConf%lineWidth_tDep%beta(sx),&
                &  do_calc_n, dbeta_dn_path_c, dbeta_dn_path_f, do_gl, del_s,&
                &  ref_corr, dhdz_gw_path, inc_rad_path_slice, tan_pt_c,     &
                &  i_stop, k_spect_dn_frq(frq_i,:) )

            ! Spectroscopic derivative  wrt: Nu0

            if ( spect_der_center ) &
              & call drad_tran_dx ( c_inds, gl_inds, del_zeta, grids_v,      &
                &  eta_zxp_v, sps_path, fwdModelConf%lineCenter%beta(sx),    &
                &  do_calc_v, dbeta_dv_path_c, dbeta_dv_path_f, do_gl, del_s,&
                &  ref_corr, dhdz_gw_path, inc_rad_path_slice, tan_pt_c,     &
                &  i_stop, k_spect_dv_frq(frq_i,:) )

          end if

        end if

        ! End of frequency loop ----------------------------------------------

        if ( toggle(emit) .and. levels(emit) > 5 ) &
          & call Trace_End ( 'ForwardModel.Frequency ', index=frq_i )

      end do            ! End freq. loop

      if ( toggle(emit) .and. levels(emit) > 4 ) &
        & call Trace_End ( 'ForwardModel.FrequencyLoop' )
    end subroutine Frequency_Loop

  ! ..........................................  Frequency_Setup_1  .....
    subroutine Frequency_Setup_1

      ! Work out which pointing frequency grid we're going to need

      ! Code splits into two sections, one for when we're doing frequency
      ! averaging, and one when we're not.

      integer :: I, K, PTG_I, ShapeInd
      integer :: MAXNOPTGFREQS     ! Used for sizing arrays
      integer :: MINSUPERSET       ! Min. value of superset > 0
      real(rp) :: R1, R2           ! real variables for various uses
      integer :: SUPERSET          ! Output from AreSignalsSuperset

      if ( fwdModelConf%do_freq_avg .and. fwdModelConf%anyLBL(sx) ) then

        whichPointingGrid = -1
        minSuperset = huge(0)
        do i = 1, size(pointingGrids)
          superset = AreSignalsSuperset ( pointingGrids(i)%signals, &
            & fwdModelConf%signals, sideband=thisSideband )
          if ( superset >= 0 .and. superset <= minSuperset ) then
            minSuperset = superset
            whichPointingGrid = i
          end if
        end do
        if ( whichPointingGrid < 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
               & "No matching pointing frequency grids." )

        ! Now we've identified the pointing grids.  Locate the tangent grid
        ! within it.
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
          shapeInd = channels(i)%shapeInds(sx)
          if ( channels(i)%dacs == 0 ) then
            k = Size(FilterShapes(shapeInd)%FilterGrid)
            r1 = FilterShapes(shapeInd)%FilterGrid(1)
            r2 = FilterShapes(shapeInd)%FilterGrid(k)
          else
            k = Size(dacsFilterShapes(shapeInd)%filter%FilterGrid)
            r1 = dacsFilterShapes(shapeInd)%filter%FilterGrid(1)
            r2 = dacsFilterShapes(shapeInd)%filter%FilterGrid(k)
          end if
          min_ch_freq_grid = MIN(r1, r2, min_ch_freq_grid)
          max_ch_freq_grid = MAX(r1, r2, max_ch_freq_grid)
        end do

        if ( FwdModelConf%anyPFA(sx) ) call get_channel_centers ( channelCenters )

      else ! ----------------------------- Not frequency averaging -----

        noFreqs = noUsedChannels
        maxNoPtgFreqs = noUsedChannels
        call get_channel_centers ( channelCenters )
        frequencies => channelCenters

      end if ! ----------------- Either frequency averaging or not -----

      call allocate_test ( RadV, maxNoPtgFreqs, 'RadV', moduleName )

      if ( FwdModelConf%anyLBL(sx) ) then
        call allocate_test ( T_Script_LBL, max_c, maxNoPtgFreqs, 'T_Script_LBL', moduleName )
        call allocate_test ( tau_LBL%tau,  max_c, maxNoPtgFreqs, 'Tau_LBL%Tau', &
          & moduleName )
        call allocate_test ( tau_LBL%i_stop, maxNoPtgFreqs, 'Tau_LBL%I_Stop', &
          & moduleName )
      end if

      call allocate_test ( inc_rad_path, max_c, maxNoPtgFreqs, 'Inc_Rad_path', &
        & moduleName )

      if ( temp_der ) &
        & call allocate_test ( k_temp_frq, max(maxNoPtgFreqs,noUsedChannels), &
                             & sv_t_len, 'k_temp_frq', moduleName )

      if ( atmos_der ) &
        & call allocate_test ( k_atmos_frq, max(maxNoPtgFreqs,noUsedChannels), &
                             & size(grids_f%values), 'k_atmos_frq', moduleName )

      if ( spect_der_width ) &
        & call allocate_test ( k_spect_dw_frq, max(maxNoPtgFreqs,noUsedChannels), &
                             & size(grids_w%values), 'k_spect_dw_frq', moduleName )
      if ( spect_der_width_TDep ) &
        & call allocate_test ( k_spect_dn_frq, max(maxNoPtgFreqs,noUsedChannels), &
                             & size(grids_n%values), 'k_spect_dn_frq', moduleName )
      if ( spect_der_center ) &
        & call allocate_test ( k_spect_dv_frq, max(maxNoPtgFreqs,noUsedChannels), &
                             & size(grids_v%values), 'k_spect_dv_frq', moduleName )

    end subroutine Frequency_Setup_1

  ! ..........................................  Frequency_Setup_2  .....
    subroutine Frequency_Setup_2 ( GridFrequencies )

      ! Work out what frequencies we're using for
      ! frequency averaging case for this pointing

      real(r8), intent(in) :: GridFrequencies(:) ! from PointingGrids
      ! Include the VELOCITY shift correction in GridFrequencies!

      integer :: J, K, L, M

      j = -1
      k = SIZE(GridFrequencies)
      call Hunt ( min_ch_freq_grid, GridFrequencies, k, j, l )
      call Hunt ( max_ch_freq_grid, GridFrequencies, k, l, m )
      noFreqs = m - j + 1
      call allocate_test ( frequencies, noFreqs, "frequencies", moduleName )

      frequencies = GridFrequencies(j:m)

    end subroutine Frequency_Setup_2

  ! ........................................  Get_Channel_Centers  .....
    subroutine Get_Channel_Centers ( channelCenters )
      real(r8) :: ChannelCenters(:)

      integer :: Channel           ! Loop inductor and subscript
      real(r8) :: Dir              ! +1 for USB, -1 for LSB
      integer :: Sig               ! Subscript for fwdModelConf%signals

      select case ( thisSideband )
      case ( -1, +1 ) ! OK
      case ( 0 )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Folded signal requested in non-frequency-averaged forward model' )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Bad value of signal%sideband' )
      end select

      dir = thisSideband
      do channel = 1, noUsedChannels
        sig = channels(channel)%signal
        channelCenters(channel) = firstSignal%lo + dir * &
          & ( fwdModelConf%signals(sig)%centerFrequency + &
          &   fwdModelConf%signals(sig)%direction * &
          &   fwdModelConf%signals(sig)%frequencies(channels(channel)%used) )
      end do

    end subroutine Get_Channel_Centers

  end subroutine FullForwardModelAuto

! =====     Private procedures     =====================================

  ! -------------------------------------------  Estimate_Tan_Phi  -----
  subroutine Estimate_Tan_Phi ( no_tan_hts, nlvl, maf, phitan, ptan, &
                              & scgeocalt, losvel, tan_press, &
                              & tan_phi, est_scgeocalt, est_los_vel )

  ! Estimate Tan_Phi and SC_geoc_alt.

    use MLSKinds, only: RP
    use MLSNumerics, ONLY: InterpolateValues
    use Units, only: Deg2Rad
    use VectorsModule, only: VectorValue_T

    implicit NONE

  ! Inputs
    integer, intent(in) :: No_Tan_Hts              ! Number of tangent heights
    integer, intent(in) :: NLVL                    ! Size of integration grid
    integer, intent(in) :: MAF                     ! MAF under consideration
    type (VectorValue_T), intent(in) :: PHITAN     ! Tangent geodAngle component of state vector
    type (VectorValue_T), intent(in) :: PTAN       ! Tangent pressure component of state vector
    type (VectorValue_T), intent(in) :: SCGEOCALT  ! S/C geocentric altitude /m
    type (VectorValue_T), intent(in) :: losvel     ! line of sight velocity by mif and maf
    real(rp), dimension(:), intent(in) :: Tan_press

  ! Outputs
    real(rp), dimension(:), intent(out) :: Tan_phi
    real(rp), dimension(:), intent(out) :: Est_scgeocalt ! Est S/C geocentric altitude /m
    real(rp), dimension(:), intent(out) :: est_los_vel

  ! Local variables
    integer :: I, JF, K, SUB
    real(rp) :: R, R1, R2, r3       ! real variables for various uses
    real(rp), dimension(ptan%template%noSurfs) :: &
      & P_PATH, &               ! Pressure on path
      & T_PATH, &               ! Temperatures on path
      & mif_vel,&               ! LOS velocities
      & Z_PATH                  ! Zeta on path

    sub = no_tan_hts - nlvl ! # subsurface levels = SurfaceTangentIndex-1

    tan_phi(1:sub) = phitan%values(1,MAF)
    est_scgeocalt(1:sub) = scGeocAlt%values(1,maf)
    est_los_vel(1:sub) = losvel%values(1,maf)

! Since the interpolateValues routine needs the OldX array to be sorted
! we have to sort ptan%values and re-arrange phitan%values, scgeocalt%values
! and losvel%values accordingly

    k = ptan%template%noSurfs

    z_path = ptan%values(1:k,maf)
    p_path = phitan%values(1:k,maf)
    t_path = scgeocalt%values(1:k,maf)
    mif_vel = losvel%values(1:k,maf)

    ! Sort z_path.  Permute p_path, t_path and mif_vel the same way.
    ! Use insertion sort since things may be nearly in order already.
    do i = 2, k ! Invariant: z_path(1:i-1) are sorted.
      r = z_path(i)
      if ( r < z_path(i-1) ) then
        r1 = p_path(i)
        r2 = t_path(i)
        r3 = mif_vel(i)
        jf = i
        do ! Find where to insert R.  Make room as we go.
          z_path(jf) = z_path(jf-1)
          p_path(jf) = p_path(jf-1)
          t_path(jf) = t_path(jf-1)
          mif_vel(jf) = mif_vel(jf-1)
          jf = jf - 1
          if ( jf == 1 ) exit
          if ( r >= z_path(jf-1) ) exit
        end do
        z_path(jf) = r
        p_path(jf) = r1
        t_path(jf) = r2
        mif_vel(jf) = r3
      end if
    end do

    call interpolateValues ( z_path, p_path, tan_press(sub+1:no_tan_hts), &
      &  tan_phi(sub+1:no_tan_hts), METHOD = 'L', EXTRAPOLATE='C' )
    call interpolateValues ( z_path, t_path, tan_press(sub+1:no_tan_hts), &
       & est_scgeocalt(sub+1:no_tan_hts), METHOD='L', EXTRAPOLATE='C' )
    call interpolateValues ( z_path, mif_vel, tan_press(sub+1:no_tan_hts), &
       & est_los_vel(sub+1:no_tan_hts), METHOD='L', EXTRAPOLATE='C' )
    tan_phi = tan_phi * deg2rad

  end subroutine Estimate_Tan_Phi

! ------------------------------------------------  not_used_here  -----
  logical function NOT_USED_HERE()
!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
!------------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module FullForwardModel_m

! $Log$
! Revision 2.277  2007/01/18 00:27:10  vsnyder
! Split Pure_Metrics into Tangent_Metrics and Height_Metrics, insert Earth
! intersection into Zeta grid.
!
! Revision 2.276  2006/12/21 22:59:00  vsnyder
! Get rid of some unused variables
!
! Revision 2.275  2006/12/21 01:34:53  vsnyder
! Finally implemented minimum Zeta
!
! Revision 2.274  2006/12/20 21:22:16  vsnyder
! Split metrics into pure H-Phi calculation, and everything else, in
! preparation for inserting the minimum-Zeta point into the path.
!
! Revision 2.273  2006/12/19 02:53:15  vsnyder
! Change some names, send max coarse path from FullForwardModel to
! FullForwardModelAuto instead of using 2*NLVL, get rid of STATUS
! argument from Metrics, reference H_Path from Earth center instead
! of surface.
!
! Revision 2.272  2006/12/13 02:32:02  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.271  2006/12/08 23:57:08  vsnyder
! Revise earth-intersecting metrics
!
! Revision 2.270  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.267  2006/08/25 19:42:31  vsnyder
! Recommited to get correct comment into the log: Compute earthradc_sq for
! metrics, since we both need it, more accurate tracing, cannonball polishing.
!
! Revision 2.266  2006/08/25 19:37:07  vsnyder
! Earthradc_sq, etc.
!
! Revision 2.265  2006/07/06 23:16:19  vsnyder
! Need to allocate dh_dt_glgrid even if no temperature derivatives
!
! Revision 2.264  2006/06/29 19:34:43  vsnyder
! Changes due to metrics handling only one tangent
!
! Revision 2.263  2006/06/19 19:26:56  vsnyder
! OOPS, out of bounds subscript possible in path_dsdh
!
! Revision 2.262  2006/06/16 20:32:30  vsnyder
! Define NGP1 in glnp
!
! Revision 2.261  2006/06/16 00:49:10  vsnyder
! Improved non-GL correction, add TeXnicalities
!
! Revision 2.260  2006/04/25 23:25:36  vsnyder
! Revise DACS filter shape data structure
!
! Revision 2.259  2006/04/21 22:25:20  vsnyder
! Cannonball polishing
!
! Revision 2.258  2006/04/11 18:34:37  vsnyder
! Add tanh(h nu / k T) to Get_D_Deltau_Pol.  Use DACS frequencies in
! Frequency_Setup_1.
!
! Revision 2.257  2006/03/06 20:44:06  vsnyder
! Avoid appearance of out-of-bounds subscript during frequency averaging
!
! Revision 2.256  2006/02/08 21:38:18  vsnyder
! Add relative humidity (RHi) calculation
!
! Revision 2.255  2006/02/08 01:02:01  vsnyder
! More stuff for spectroscopy derivatives
!
! Revision 2.254  2006/01/05 00:03:52  vsnyder
! Implement refractive correction for Phi
!
! Revision 2.253  2005/12/10 01:51:54  vsnyder
! Cannonball polishing
!
! Revision 2.252  2005/12/07 19:43:32  vsnyder
! Mistakenly committed needing Phi_Refractive_Correction
!
! Revision 2.251  2005/12/07 01:30:04  vsnyder
! More on getting correct size for RadV
!
! Revision 2.250  2005/12/07 00:35:04  vsnyder
! Allocate RadV with correct size for PFA and no LBL
!
! Revision 2.249  2005/11/21 22:57:41  vsnyder
! PFA derivatives stuff, plug some memory leaks
!
! Revision 2.248  2005/11/05 03:38:13  vsnyder
! Frequency_Average_Derivative doesn't need Tau, cannonball polishing
!
! Revision 2.247  2005/11/03 03:57:45  vsnyder
! Don't try to look at filter shapes for DACS
!
! Revision 2.246  2005/11/02 21:38:48  vsnyder
! Repair a broken tracing message
!
! Revision 2.245  2005/11/01 23:02:08  vsnyder
! PFA Derivatives, use precomputed ShapeInds
!
! Revision 2.244  2005/09/17 01:02:38  vsnyder
! Out of bounds subscript
!
! Revision 2.243  2005/09/17 00:49:53  vsnyder
! Revise arrays for spectroscopic derivatives, plus some cannonball polishing
!
! Revision 2.242  2005/09/03 01:21:33  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.241  2005/08/03 18:03:42  vsnyder
! Scan averaging, some spectroscopy derivative stuff
!
! Revision 2.240  2005/07/08 19:40:51  vsnyder
! OOPS, forgot to nullify Rad_FFT
!
! Revision 2.239  2005/07/08 00:12:11  vsnyder
! Get Rad_FFT from Convolve_Radiance to Convolve_Temperature_Deriv
!
! Revision 2.238  2005/07/06 02:17:20  vsnyder
! Revisions for spectral parameter derivatives
!
! Revision 2.237  2005/06/09 02:34:16  vsnyder
! Move stuff from l2pc_pfa_structures to slabs_sw_m
!
! Revision 2.236  2005/06/03 01:59:25  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.235  2005/04/19 20:16:45  livesey
! Added a couple of nullifys
!
! Revision 2.234  2005/03/28 20:26:25  vsnyder
! Lots of PFA stuff
!
! Revision 2.233  2005/03/10 00:28:09  pwagner
! Better patching of ptg_angles avoids list-out-of-order in Hunt
!
! Revision 2.232  2005/02/17 02:35:29  vsnyder
! Do PFA on fine path if necessary
!
! Revision 2.231  2005/02/16 23:16:49  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.230  2004/12/28 00:28:02  vsnyder
! Remove unreferenced use name
!
! Revision 2.229  2004/12/13 20:37:23  vsnyder
! Moved stuff from get_species_data to ForwardModelConfig, some cannonball polishing
!
! Revision 2.228  2004/11/05 19:38:39  vsnyder
! Got rid of DerivedFromForwardModel component of config
!
! Revision 2.227  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
! Revision 2.226  2004/10/13 01:08:27  vsnyder
! Moved some checking to ForwardModelSupport
!
! Revision 2.225  2004/10/07 23:26:09  vsnyder
! Changes in Beta_Group structure
!
! Revision 2.224  2004/10/06 21:27:41  vsnyder
! More work on PFA -- seems to work now
!
! Revision 2.223  2004/09/14 17:20:45  bill
! add mif dependent los vel correction-wgr
!
! Revision 2.222  2004/09/04 01:50:30  vsnyder
! get_beta_path_m.f90
!
! Revision 2.221  2004/09/01 01:48:27  vsnyder
! Status flags, more work on PFA
!
! Revision 2.220  2004/08/06 22:40:17  livesey
! Better patch for ptg_angles
!
! Revision 2.219  2004/08/06 01:40:39  livesey
! Upgraded the precision of the ptg dump.
!
! Revision 2.218  2004/08/06 01:24:55  livesey
! Typo!
!
! Revision 2.217  2004/08/06 01:24:22  livesey
! Minor bug fix in ptg_angles dump
!
! Revision 2.216  2004/08/05 20:53:50  vsnyder
! More PFA preparations, some cannonball polishing
!
! Revision 2.215  2004/08/05 20:24:00  livesey
! Added ptg switch to dump ptg_angles
!
! Revision 2.214  2004/08/03 22:07:10  vsnyder
! Inching further toward PFA
!
! Revision 2.213  2004/07/30 19:53:18  livesey
! Bug fix, forbid extrapolation in Estimate_tan_phi
!
! Revision 2.212  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.211  2004/06/10 00:59:56  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.210  2004/05/27 23:23:50  pwagner
! named parameter clean= to dump procedures
!
! Revision 2.209  2004/05/17 22:05:11  livesey
! A change to refraction and to k_atmos handling to avoid explosions.
!
! Revision 2.208  2004/04/24 02:27:05  vsnyder
! Cosmetic changes
!
! Revision 2.207  2004/04/19 21:01:37  vsnyder
! Put size of gl_slabs in call to get_gl_slabs_arrays
!
! Revision 2.206  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.205  2004/04/05 21:08:42  jonathan
! delet temp_prof
!
! Revision 2.204  2004/04/02 00:59:24  vsnyder
! Get catalog from slabs structure
!
! Revision 2.203  2004/03/31 20:35:26  jonathan
! bug fix in handling clouds
!
! Revision 2.202  2004/03/30 02:26:17  livesey
! Bug fix in Jonathan's w0 handling
!
! Revision 2.201  2004/03/30 02:00:36  vsnyder
! Remove USE for unreferenced symbol.  Don't try to fill tpath_m and
! tpath_p if they're not allocated.
!
! Revision 2.200  2004/03/27 03:35:27  vsnyder
! Add pointer to catalog in slabs_struct.  Use it so as not to need to drag
! line centers and line widths around.  Write slabs_lines and slabswint_lines
! to get sum of beta over all lines; put slabs_struct instead of its components
! in the calling sequence.
!
! Revision 2.199  2004/03/20 04:05:50  vsnyder
! Moved SpeedOfLight from units to physics
!
! Revision 2.198  2004/03/20 01:15:16  jonathan
! add in scattering correction term in two t_script
!
! Revision 2.196  2004/03/01 19:22:14  jonathan
! following the changes made to load_one_item
!
! Revision 2.195  2004/02/14 00:23:48  vsnyder
! New DACS convolution algorithm
!
! Revision 2.194  2004/02/05 23:30:01  livesey
! Finally implemented code to do correct handing of sideband fraction in
! single sideband radiometers.
!
! Revision 2.193  2004/02/03 02:48:35  vsnyder
! Progress (hopefully) on polarized temperature derivatives.
! Implement DACs frequency convolution.
!
! Revision 2.192  2004/01/23 19:13:42  jonathan
! add an extra-flag for tscat in load_one_item
!
! Revision 2.191  2003/12/08 21:38:33  jonathan
! some minor changes
!
! Revision 2.190  2003/12/08 17:52:02  jonathan
! update for 2d cldfwm
!
! Revision 2.189  2003/12/07 19:45:46  jonathan
! update for 2D cloud FWM
!
! Revision 2.188  2003/12/01 17:24:05  jonathan
! add scat_alb
!
! Revision 2.187  2003/11/24 22:10:14  vsnyder
! Multiply beta_path_polarized_f by tanh1_f if derivatives needed
!
! Revision 2.186  2003/11/20 17:33:50  pwagner
! Nullify some things otherwise left unassociated
!
! Revision 2.185  2003/11/19 22:21:34  jonathan
! interpolate scat_src to tscat_path
!
! Revision 2.184  2003/11/17 18:04:15  jonathan
! scat_src output from T_scat is correct
!
! Revision 2.183  2003/11/12 00:10:55  jonathan
! some changes due to cloud construction
!
! Revision 2.182  2003/11/07 03:18:49  vsnyder
! Cosmetic changes
!
! Revision 2.181  2003/11/06 01:13:14  bill
! fixed DACS freq extrapolation problem
!
! Revision 2.180  2003/11/05 19:26:01  jonathan
! add stuff for use later in cld model
!
! Revision 2.179  2003/11/04 02:49:13  vsnyder
! Calculate coarse-path indices where GL is needed
!
! Revision 2.178  2003/11/04 01:57:15  vsnyder
! Move trapezoid correction to a better place, cosmetics
!
! Revision 2.177  2003/11/03 23:15:16  vsnyder
! Get rid of path_ds_dh procedure -- a one-liner used in one place
!
! Revision 2.176  2003/11/01 03:02:57  vsnyder
! Compute del_zeta, ds_dz_gw and dh_dz_gw for [d]rad_tran*; change
! indices_c to c_inds for consistency with usage in rad_tran_m.
!
! Revision 2.175  2003/10/30 20:37:00  vsnyder
! Compute del_zeta here for *rad_tran_*
!
! Revision 2.174  2003/10/29 00:43:51  livesey
! Added support for the forceFoldedOutput option
!
! Revision 2.173  2003/10/28 22:05:53  jonathan
! add l_gph for use in cloud model
!
! Revision 2.172  2003/10/09 19:30:18  vsnyder
! Cosmetic changes
!
! Revision 2.171  2003/09/24 02:53:48  vsnyder
! Cosmetic changes
!
! Revision 2.170  2003/09/09 00:04:27  vsnyder
! Supply E and Sqrt_Earth_Rflty to mcrt_der
!
! Revision 2.169  2003/08/16 01:14:58  vsnyder
! Use radiometers%polarization to choose 1,1 or 2,2 element
!
! Revision 2.168  2003/08/15 20:29:26  vsnyder
! Implement polarized VMR derivatives
!
! Revision 2.167  2003/08/15 18:50:21  vsnyder
! Preparing the way for polarized vmr derivatives
!
! Revision 2.166  2003/08/12 23:07:32  vsnyder
! Futzing with comments
!
! Revision 2.165  2003/08/12 21:58:37  vsnyder
! Use trapezoid instead of rectangle to integrate non-GL panels
!
! Revision 2.164  2003/08/12 18:22:10  michael
! Contribution of scalar molecules to polarized absorption is now added to
! coarse grid alpha_path_polarized (1/4 1/2 1/4) instead of to diagonal of
! incoptdepth_pol.  This make alpha_path_polarized correct for later use
! in gl corrections.
!
! Revision 2.163  2003/07/15 23:07:21  vsnyder
! Simplify Freq_Avg
!
! Revision 2.162  2003/07/15 22:10:09  livesey
! Added support for hybridModel
!
! Revision 2.161  2003/07/15 18:16:48  livesey
! Catalog now split by sideband, also changed no_ele to max_ele in
! allocates
!
! Revision 2.160  2003/07/09 23:40:13  vsnyder
! Use new AllocateSlabs routine
!
! Revision 2.159  2003/07/09 22:27:42  vsnyder
! More futzing
!
! Revision 2.158  2003/07/09 20:23:50  vsnyder
! Futzing
!
! Revision 2.157  2003/07/09 20:14:19  livesey
! Anticipative bug fix commented out.
!
! Revision 2.156  2003/07/04 03:40:13  vsnyder
! Simplify dump in case exp(incoptdepth_pol) fails
!
! Revision 2.155  2003/07/04 02:50:15  vsnyder
! Simplify interface to Get_GL_Slabs_Arrays, correct a blunder introduced around July 17
!
! Revision 2.154  2003/06/27 23:43:34  vsnyder
! Remove unreferenced USE names
!
! Revision 2.153  2003/06/27 22:09:48  vsnyder
! Check status from rad_tran_pol and report an error if overflow occurred
!
! Revision 2.152  2003/06/27 00:59:53  vsnyder
! Simplify interface to Get_Species_Data
!
! Revision 2.151  2003/06/25 02:41:37  vsnyder
! Futzing
!
! Revision 2.150  2003/06/18 22:26:41  vsnyder
! Restored the check to the wrong place at 2.149
!
! Revision 2.149  2003/06/18 19:29:30  vsnyder
! Replace a check inadvertently deleted in rev 2.146
!
! Revision 2.148  2003/06/18 17:23:16  bill
! fixed NAG associated bug
!
! Revision 2.147  2003/06/18 14:54:04  bill
! added subsetting feature for T-ders
!
! Revision 2.146  2003/06/18 01:59:20  vsnyder
! Move checking that all signals in config are for same radiometer, module
! and sideband, and computation of sidebandStart and sidebandStop, to
! ForwardModelSupport.  Remove vector quantity validation, as that's now
! done in Construct.  Futzing.
!
! Revision 2.145  2003/06/13 00:00:27  vsnyder
! Move multiplication of beta_path by tanh into FullForwardModel
!
! Revision 2.144  2003/06/10 15:06:54  bill
! fixed polarized t-derivs
!
! Revision 2.143  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.142  2003/06/09 17:38:47  livesey
! Use GetQuantityForForwardModel in more places
!
! Revision 2.141  2003/05/29 16:37:38  livesey
! Renamed sideband fraction
!
! Revision 2.140  2003/05/26 01:42:50  michael
! Two temporary fixes only relevant to the polarized model.
! Added a bug-fix removing scalar contribution to magnetic GL.
! I don't know why this works but it makes agreement with scalar model
! much better.  I also added an undocumented switch -Scrosspol to allow
! the use of the (2,2) element of the radiance tensor rather than the (1,1).
! This code should be removed when we set up the l2cf to do this properly.
!
! Revision 2.139  2003/05/22 20:01:17  vsnyder
! Cosmetic changes
!
! Revision 2.138  2003/05/22 04:03:41  livesey
! Now handles elevation offset as a channel by channel quantity
!
! Revision 2.137  2003/05/20 00:05:39  vsnyder
! Move some stuff to subroutines
!
! Revision 2.136  2003/05/15 03:29:44  vsnyder
! Implement polarized model's temperature derivatives
!
! Revision 2.135  2003/05/09 19:26:36  vsnyder
! Expect T+DT instead of T and DT separately in Get_GL_Slabs_Arrays,
! initial stuff for temperature derivatives of polarized radiance.
!
! Revision 2.134  2003/05/05 23:00:24  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.126.2.40  2003/04/24 21:57:05  vsnyder
! Check FwdModelConf%incl_cld instead of associated(cloudIce)
!
! Revision 2.126.2.39  2003/04/24 19:01:19  vsnyder
! Compute MIF correctly
!
! Revision 2.126.2.38  2003/04/16 20:44:20  vsnyder
! Move allocate for scat_src%values out of the loops
!
! Revision 2.126.2.37  2003/04/16 19:32:02  vsnyder
! Move working storage for T_Scat into T_Scat
!
! Revision 2.126.2.36  2003/04/15 23:45:45  vsnyder
! Correct choice of mag field coordinate indices, futzing
!
! Revision 2.126.2.35  2003/04/15 17:55:19  jonathan
! FullForwardModel_m.f90 change scat_src as non-pointer
!
! Revision 2.126.2.34  2003/04/15 14:18:22  jonathan
! interpolate scat_src%values to LOS
!
! Revision 2.126.2.33  2003/04/10 16:16:06  jonathan
! fix bug
!
! Revision 2.126.2.32  2003/04/08 21:44:23  jonathan
! remove if condition for Modify_values_for_supersat
!
! Revision 2.126.2.31  2003/04/07 17:13:52  jonathan
! modified cal to T_scat
!
! Revision 2.126.2.30  2003/04/04 20:43:53  jonathan
! update call to T_scat
!
! Revision 2.126.2.29  2003/03/28 18:14:19  jonathan
! add new scaterring source module
!
! Revision 2.126.2.28  2003/03/27 23:32:19  vsnyder
! Houst a reshape out of a loop
!
! Revision 2.126.2.27  2003/03/27 00:50:35  vsnyder
! Get ECRtoFOV, use it to rotate mag_path
!
! Revision 2.126.2.26  2003/03/22 04:19:13  vsnyder
! Cosmetic changes
!
! Revision 2.126.2.25  2003/03/22 04:03:04  vsnyder
! Move Beta_Group_T and Dump_Beta_Group from get_beta_path to Get_Species_Data
!
! Revision 2.126.2.24  2003/03/22 02:48:53  vsnyder
! Interpolate the magnetic field onto the path
!
! Revision 2.126.2.23  2003/03/21 02:49:55  vsnyder
! Use an array of pointers to quantities instead of searching several times
! using GetQuantityForForwardModel.  Move stuff around.  Cosmetic changes.
!
! Revision 2.126.2.22  2003/03/20 19:21:05  vsnyder
! More futzing with grids_t and stuff that uses it
!
! Revision 2.126.2.21  2003/03/20 01:42:25  vsnyder
! Revise Grids_T structure
!
! Revision 2.126.2.19  2003/03/07 23:17:05  vsnyder
! Use ASSOCIATED instead of PRESENT for "optional" arguments GET_CHI_OUT
!
! Revision 2.126.2.18  2003/03/06 22:47:29  vsnyder
! Polarized radiative transfer sort-of works
!
! Revision 2.126.2.17  2003/03/05 03:40:07  vsnyder
! Do tanh correction for polarized, more polarized work, simplifications
!
! Revision 2.126.2.16  2003/03/04 20:24:44  dwu
! add a temporay fix for the tangent height crossover problem
!
! Revision 2.126.2.15  2003/03/01 03:19:28  vsnyder
! More work on polarized radiative transfer.  Still doesn't work, but at
! least it doesn't break the nonpolarized case.
!
! Revision 2.126.2.14  2003/02/28 00:15:59  vsnyder
! Print 3 digits in channel number in 'rad' diagnostic output
!
! Revision 2.126.2.13  2003/02/27 23:35:01  vsnyder
! Revise polarized processing
!
! Revision 2.126.2.12  2003/02/25 00:53:09  jonathan
! add grids_tscat
!
! Revision 2.126.2.11  2003/02/24 23:40:31  jonathan
! change input for get_beta_path_cloud
!
! Revision 2.126.2.10  2003/02/24 23:26:43  jonathan
! change temp_supersat to temp_prof
!
! Revision 2.126.2.9  2003/02/22 00:49:26  vsnyder
! Delete moleculesPol, moleculeDerivativesPol, add Polarized to ForwardModelConfig
!
! Revision 2.126.2.8  2003/02/21 21:34:32  michael
! Made index for dumping dacs radiances be 1:129 for dacs channels 0:128
!
! Revision 2.126.2.7  2003/02/15 00:29:11  vsnyder
! Don't exp(incoptdepth) in mcrt -- it's done here
!
! Revision 2.126.2.6  2003/02/14 23:31:14  vsnyder
! Delete unreferenced names
!
! Revision 2.126.2.5  2003/02/14 23:27:31  vsnyder
! Move stuff to Get_Species_Data
!
! Revision 2.126.2.4  2003/02/14 03:53:41  vsnyder
! Initial commit of polarized model
!
! Revision 2.126.2.3  2003/02/14 00:21:42  jonathan
! add singl. scat. albedo W0, ph funct PHH
!
! Revision 2.126.2.2  2003/02/13 22:26:18  jonathan
! changes dimension for beta_path_cloud also delocate it
!
! Revision 2.126.2.1  2003/02/13 17:35:01  bill
! fixes gl_ind bug and interfaces to get_beta
!
! Revision 2.126  2003/02/11 00:48:18  jonathan
! changes made after adding get_beta_path_cloud
!
! Revision 2.125  2003/02/08 01:03:00  livesey
! Bug fix in call to rad_tran
!
! Revision 2.124  2003/02/07 01:07:41  jonathan
! add in option to compute dry and super-saturation case in load_sps
!
! Revision 2.123  2003/02/06 22:04:25  vsnyder
! Add f_moleculesPol, f_moleculeDerivativesPol, delete f_polarized
!
! Revision 2.122  2003/02/06 19:15:47  jonathan
!  fix a bug
!
! Revision 2.121  2003/02/06 19:06:49  jonathan
! add eta_iwc, eta_iwc_zp, do_calc_iwc, do_cala_iwc_zp and also not passing
! through comp_eta_docalc and comp_sps_path_frq if fwdModelConf%Incl_Cld is
! false
!
! Revision 2.120  2003/02/06 05:55:47  livesey
! Fix to sort of fix Jonathan's cloud ice stuff.
!
! Revision 2.119  2003/02/06 00:20:08  jonathan
! Add in many stuff to deal with clouds CloudIce, iwc_path, etc
!
! Revision 2.118  2003/02/03 23:18:43  vsnyder
! Squash a bug in deallocating beta_path_polarized
!
! Revision 2.117  2003/02/03 22:58:17  vsnyder
! Plug a memory leak, delete gl_ndx, some polarized stuff
!
! Revision 2.116  2003/02/03 19:00:36  bill
! changed interface to rad tran to speed up program
!
! Revision 2.115  2003/02/01 02:33:22  vsnyder
! Plug a bunch of memory leaks
!
! Revision 2.114  2003/01/31 17:53:39  jonathan
! change z_path to z_path_c in passing to get_beta_path
!
! Revision 2.113  2003/01/31 17:15:49  jonathan
! add Inc_Cld to get_beta_path
!
! Revision 2.112  2003/01/31 01:53:01  vsnyder
! Move array temps to arrays explicitly allocated outside the loop
!
! Revision 2.111  2003/01/30 00:16:35  jonathan
! add z_path to get_beta_path
!
! Revision 2.110  2003/01/26 04:42:42  livesey
! Added profiles/angle options for phiWindow
!
! Revision 2.109  2003/01/21 18:20:32  vsnyder
! Put dimensions back onto actual arguments to path_contrib
!
! Revision 2.108  2003/01/18 03:36:09  vsnyder
! Undo ill-advised cosmetic changes -- that weren't cosmetic
!
! Revision 2.106  2003/01/16 23:13:03  livesey
! Added MaxRefraction stuff
!
! Revision 2.105  2003/01/16 18:04:01  jonathan
! add Do_1D option to get_gl_slabs_arrays
!
! Revision 2.104  2003/01/14 21:48:58  jonathan
! add i_saturation
!
! Revision 2.103  2003/01/10 21:55:26  vsnyder
! Move SpeedOfLight from Geometry ot Units
!
! Revision 2.102  2003/01/08 00:16:39  vsnyder
! Use "associated" instead of "present" to control optional computations.
! Cosmetic changes, too.
!
! Revision 2.101  2002/12/12 01:12:47  vsnyder
! Let InvalidQuantity have a length > 1
!
! Revision 2.100  2002/11/15 01:33:08  livesey
! Added allLinesForRadiometer functionality.
!
! Revision 2.99  2002/11/13 17:07:44  livesey
! Passes FwdModelExtra into convolve/no_conv
!
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
