! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullForwardModel_m

  ! This module contains the `full' forward model.

  implicit NONE
  private
  public :: FullForwardModel

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------
contains

!---------------------------------------------  FullForwardModel  -------------
  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                             &  FwdModelOut, oldIfm, FmStat, Jacobian )

  ! This is the full radiative transfer forward model, the workhorse code

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use AntennaPatterns_m, only: ANTENNAPATTERNS
    use Comp_Eta_Docalc_No_Frq_m, only: Comp_Eta_Docalc_No_Frq
    use Comp_Sps_Path_Frq_m, only: Comp_Sps_Path, Comp_Sps_Path_Frq
    use Compute_GL_Grid_m, only: Compute_GL_Grid
    use Convolve_All_m, only: Convolve_All
    use CS_Expmat_m, only: CS_Expmat
    use D_Hunt_m, only: Hunt
    use D_T_SCRIPT_DTNP_M, only: DT_SCRIPT_DT
    use DO_T_SCRIPT_M, ONLY: TWO_D_T_SCRIPT
    use Dump_0, only: Dump
    use Eval_Spect_Path_m, only: Eval_Spect_Path
    use FilterShapes_m, only: FilterShapes, DACSFilterShapes
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelIntermediate_t, &
                                    &   ForwardModelStatus_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel, QtyStuff_T
    use Freq_Avg_m, only: Freq_Avg, Freq_Avg_DACS
    use Geometry, only: EarthRadA, EarthRadB, MaxRefraction
    use Get_Beta_Path_m, only: Get_Beta_Path, Get_Beta_Path_Cloud, &
      & Get_Beta_Path_Polarized
    use Get_Chi_Angles_m, only: Get_Chi_Angles
    use Get_Chi_Out_m, only: Get_Chi_Out
    use Get_d_Deltau_pol_m, only: Get_d_Deltau_pol_df, Get_d_Deltau_pol_dT
    use Get_Species_Data_M, only: Beta_Group_T, Get_Species_Data, &
      Destroy_Species_Data, Destroy_Beta_Group
    use GLnp, only: GW, NG
    use Intrinsic, only: L_A, L_BOUNDARYPRESSURE, L_CLEAR, L_CLOUDICE, &
      & L_CLOUDWATER, L_DN, L_DV, L_DW, L_EARTHREFL, L_ECRtoFOV, &
      & L_ELEVOFFSET, L_GPH, L_LOSVEL, L_LIMBSIDEBANDFRACTION, L_MAGNETICFIELD, &
      & L_ORBITINCLINATION, L_PHITAN, L_PTAN, L_RADIANCE, L_REFGPH, L_SCGEOCALT, &
      & L_SIZEDISTRIBUTION, L_SPACERADIANCE, L_TEMPERATURE, L_VMR, &
      & LIT_INDICES
    use Load_Sps_Data_m, only: DestroyGrids_t, Grids_T, Load_One_Item_Grid, &
      & Load_Sps_Data, Modify_values_for_supersat, create_grids_1, fill_grids_1, create_grids_2, fill_grids_2
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATESLABS, &
                                &  DESTROYCOMPLETESLABS
    use ManipulateVectorQuantities, only: DoHGridsMatch, FindClosestInstances
    use MatrixModule_1, only: MATRIX_T
    use Mcrt_m, only: Mcrt_der
    use Metrics_m, only: Metrics
    use MLSCommon, only: R4, R8, RP, IP, FindFirst
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
      & MLSMSG_Error, MLSMSG_Warning
    use MLSNumerics, only: Hunt, InterpolateValues
    use MLSSignals_m, only: AreSignalsSuperset, GetNameOfSignal, MatchSignal, &
      & Radiometers, Signal_t
    use Molecules, only: FIRST_MOLECULE, LAST_MOLECULE,                      &
                       & L_Extinction, & ! Used in include dump_print_code.f9h
                       & L_H2O, L_H2O_18, L_N2, L_N2O, L_O_18_O, L_O2, L_O3
    use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
    use Opacity_m, only: Opacity
    use Output_m, only: Output
    use Path_Contrib_M, only: Get_GL_Inds, Path_Contrib
    use Physics, only: H_OVER_K
    use PointingGrid_m, only: POINTINGGRIDS
    use RAD_TRAN_M, only: RAD_TRAN, RAD_TRAN_POL, DRAD_TRAN_DF, DRAD_TRAN_DT, &
      & DRAD_TRAN_DX
    use REFRACTION_M, only: REFRACTIVE_INDEX, COMP_REFCOR
    use ScatSourceFunc, only: T_SCAT,  Interp_Tscat,  Convert_Grid
    use SpectroscopyCatalog_m, only: CATALOG_T
    use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
    use String_table, only: GET_STRING ! Used in include dump_print_code.f9h
! use testfield_m
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use TWO_D_HYDROSTATIC_M, only: Two_D_Hydrostatic
    use Units, only: Deg2Rad, SpeedOfLight
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQUANTITYBYTYPE

    type(forwardModelConfig_T), intent(in) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    !--------------------------------------------
    real(r8), dimension(:,:), pointer :: WC
    real(r8), dimension(:), pointer :: Scat_ang
    integer, dimension(:), pointer :: IPSD
    !--------------------------------------------

    ! Define local parameters
    real(rp), parameter :: DEL_TEMP = 10.0_rp ! Temp. step size for derivatives
    integer, parameter :: Ngp1 = Ng+1         ! NG + 1

    ! Now define local variables, group by type and then
    ! alphabetically

    integer :: BRKPT                    ! Index of midpoint of path
    integer :: CHANIND                  ! A 1 based channel index
    integer :: CHANNEL                  ! A Loop counter
    integer :: DIRECTION                ! Direction of channel numbering
    integer :: EXT_IND                  ! Index of extinction inside f array
    integer :: F_LEN_DN                 ! Length of DN in vector
    integer :: F_LEN_DV                 ! Length of DV in vector
    integer :: F_LEN_DW                 ! Length of DW in vector
    integer :: FRQ_I                    ! Frequency loop index
    integer :: H2O_IND                  ! Index of h2o inside f array, else zero
    integer :: IER                      ! Status flag from allocates
    integer :: I                        ! Loop index and other uses .
    integer :: INSTANCE                 ! Loop counter
    integer :: INST                     ! Relevant instance for temperature
    integer :: I_STOP                   ! Upper index for radiance comp.
    integer :: J                        ! Loop index and other uses ..
    integer :: K                        ! Loop index and other uses ..
    integer :: L                        ! Loop index and other uses ..
    integer :: M                        ! Loop index and other uses ..
    integer :: MAF                      ! MAF under consideration
    integer :: MAX_ELE                  ! Length of longest possible path (all no_ele<max_ele)
    integer :: MAXNOFFREQS              ! Max. no. frequencies for any molecule
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXNOPTGFREQS            ! Used for sizing arrays
    integer :: MAXVERT                  ! Total number of points in gl grid
                                        ! NLVL * (NG+1) - NG, i.e., 1 + NG
                                        ! per level, except the last, where
                                        ! there's no GL space.
    integer :: MID                      ! NPC / 2
    integer :: MIF                      ! MIF number for tan_press(ptg_i)
    integer :: MINSUPERSET              ! Min. value of superset > 0
    integer :: NCG                      ! Number of panels needing GL = Size(cg_inds)
    integer :: NGL                      ! Total # of GL points = Size(gl_inds)
    integer :: NGLMAX                   ! NGL if all panels need GL
    integer :: Nlvl                     ! Size of integration grid
    integer :: NO_ELE                   ! Length of a gl path
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NO_MOL                   ! Number of major molecules (NO iso/vib)
    integer :: NOSPECIES                ! No. of molecules under consideration
    integer :: NO_TAN_HTS               ! Number of tangent heights
    integer :: NOUSEDDACS               ! Number of different DACS in this run.
    integer :: NOUSEDCHANNELS           ! How many channels are we considering
    integer :: NPC                      ! Length of coarse path
    integer :: N_T_ZETA                 ! Number of zetas for temperature
    integer :: P_Stop                   ! Where to stop in polarized case
    integer :: PTG_I                    ! Loop counter for the pointings
    integer :: SHAPEIND                 ! Index into filter shapes
    integer :: SIGIND                   ! Signal index, loop counter
    integer :: SPECIE                   ! Loop counter
    integer :: SUPERSET                 ! Output from AreSignalsSuperset
    integer :: SURFACE                  ! Loop counter
    integer :: SURFACETANGENTINDEX      ! Index in tangent grid of earth's
                                        ! surface
    integer :: SV_I                     ! Loop index and other uses .
    integer :: SV_T_LEN                 ! Number of t_phi*t_zeta in the window
    integer :: THISSIDEBAND             ! Loop counter for sidebands
    integer :: WHICHPATTERN             ! Index of antenna pattern
    integer :: WHICHPOINTINGGRID        ! Index into the pointing grids
    integer :: WINDOWFINISH             ! End of temperature `window'
    integer :: WINDOWSTART              ! Start of temperature `window'

    integer :: noSurf                   ! Number of pressure levels
    integer :: novmrSurf                ! Number of vmr levels
    integer :: nspec                    ! No of species for cloud model
    integer :: ispec                    ! Species index in cloud model
    integer, dimension(:), pointer :: closestInstances

    logical :: Any_Der                  ! temp_der .or. atmos_der .or. spect_der
    logical :: cld_fine = .false.
    logical :: do_cld = .false.
    logical :: temp_der, atmos_der, spect_der, ptan_der ! Flags for various derivatives
    logical :: Update                   ! Just update radiances etc.

    character (len=32) :: SigName       ! Name of a Signal
    character (len=32) :: molName       ! Name of a molecule

    logical :: Clean                    ! Used for dumping
    logical :: dummy(2) = (/.FALSE.,.FALSE./)  ! dummy Flag array

    integer, dimension(:), pointer :: C_INDS  ! Indices on coarse grid
    integer, dimension(:), pointer :: CG_INDS ! Indices on coarse grid where GL needed
    integer, dimension(:), pointer :: F_INDS  ! Indices on fine grid
    integer, dimension(:), pointer :: GL_INDS ! Index of GL points -- subset of f_inds
    integer, dimension(:), pointer :: GRIDS ! Heights in ptgGrid for each tangent
    integer, dimension(:), pointer :: TAN_INDS ! Index of tangent grid into gl grid
    integer, dimension(:), pointer :: USEDDACSSIGNALS ! Indices in FwdModelConf
                                         ! of signals for our dacs

    logical, dimension(:), pointer :: DO_GL ! GL indicator

    logical, dimension(:,:), pointer :: DO_CALC_DN    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DN_C  ! DO_CALC_DN on coarse grid
    logical, dimension(:,:), pointer :: DO_CALC_DN_F  ! DO_CALC_DN on fine grid
    logical, dimension(:,:), pointer :: DO_CALC_DV    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DV_C  ! DO_CALC_DV on coarse grid
    logical, dimension(:,:), pointer :: DO_CALC_DV_F  ! DO_CALC_DV on fine grid
    logical, dimension(:,:), pointer :: DO_CALC_DW    ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_DW_C  ! DO_CALC_DW  on coarse grid
    logical, dimension(:,:), pointer :: DO_CALC_DW_F  ! DO_CALC_DW  on fine grid
    logical, dimension(:,:), pointer :: DO_CALC_FZP   ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_IWC   ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_Tscat ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_Salb  ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_cext  ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_HYD   ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_HYD_C ! DO_CALC_HYD on coarse grid
    logical, dimension(:,:), pointer :: DO_CALC_T     ! 'Avoid zeros'
    logical, dimension(:,:), pointer :: DO_CALC_T_C   ! DO_CALC_T on coarse grid
    logical, dimension(:,:), pointer :: DO_CALC_T_F   ! DO_CALC_T on fine grid
    logical, dimension(:,:), pointer :: DO_CALC_ZP    ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_IWC_ZP   ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_Tscat_ZP ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_Salb_ZP  ! 'Avoid zeros' indicator
    logical, dimension(:,:), pointer :: DO_CALC_Cext_ZP  ! 'Avoid zeros' indicator
    logical :: Got( FIRST_MOLECULE : LAST_MOLECULE )  !
    LOGICAL, DIMENSION(:), pointer :: true_path_flags ! array of trues
    LOGICAL, DIMENSION(:), pointer :: t_der_path_flags! a flag that tells the
! where an absorption coefficient is needed for a temperature derivative.
! Only useful when subsetting temperature derivatives.

! Array of Flags indicating  which Temp. coefficient to process

    real(rp) :: E_RFLTY       ! Earth reflectivity at given tan. point
    real(rp), save :: E_Stop  = 1.0_rp ! X for which Exp(X) is too small to worry
    real(r8) :: FRQ           ! Frequency
    real(r8) :: FRQHK         ! 0.5 * Frq * H_Over_K
    real(rp) :: NEG_TAN_HT    ! GP Height (in KM.) of tan. press.
                              ! below surface
    real(rp) :: R,R1,R2       ! real variables for various uses
    real(rp) :: REQ           ! Equivalent Earth Radius
    real(rp) :: ROT(3,3)      ! ECR-to-FOV rotation matrix
    real(rp) :: THISFRACTION     ! A sideband fraction
    real(rp) :: THISELEV      ! An elevation offset
    real(rp) :: Vel_Cor       ! Velocity correction due to Vel_z

    real(rp), dimension(1) :: ONE_TAN_HT ! ***
    real(rp), dimension(1) :: ONE_TAN_TEMP ! ***
    real(rp), dimension(:), pointer :: ALPHA_PATH_C ! coarse grid abs coeff.
    real(rp), dimension(:), pointer :: ALPHA_PATH_F ! fine grid abs coeff.
    real(rp), dimension(:), pointer :: BETA_PATH_cloud_C ! Beta on path coarse
    real(rp), dimension(:), pointer :: CT           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), dimension(:), pointer :: DEL_S        ! Integration lengths along path
    real(rp), dimension(:), pointer :: DEL_ZETA     ! Integration lengths in Zeta coords
    real(rp), dimension(:), pointer :: DHDZ_PATH    ! dH/dZ on path
    real(rp), dimension(:), pointer :: DHDZ_GW_PATH ! dH/dZ * GW on path
    real(rp), dimension(:), pointer :: DRAD_DN      ! dI/dN
    real(rp), dimension(:), pointer :: DRAD_DT      ! dI/dT
    real(rp), dimension(:), pointer :: DRAD_DV      ! dI/dV
    real(rp), dimension(:), pointer :: DRAD_DW      ! dI/dW
    real(rp), dimension(:), pointer :: DSDZ_GW_PATH ! ds/dH * dH/dZ * GW on path
    real(r8), dimension(:), pointer :: FREQUENCIES  ! Frequencies to compute for
    real(rp), dimension(:), pointer :: H            ! Magnetic field on path, in
                                                    ! IFOVPP
    real(rp), dimension(:), pointer :: H_PATH       ! Heights on path
    real(rp), dimension(:), pointer :: H_PATH_C     ! H_PATH on coarse grid
    real(rp), dimension(:), pointer :: H_PATH_F     ! H_PATH on fine grid
    real(rp), dimension(:), pointer :: INCOPTDEPTH  ! Incremental Optical depth
    real(rp), dimension(:), pointer :: N_PATH       ! Refractivity on path
    real(rp), dimension(:), pointer :: PATH_DSDH    ! dS/dH on path
    real(rp), dimension(:), pointer :: P_GLGRID     ! Pressure on glGrid surfs
    real(rp), dimension(:), pointer :: PHI_PATH     ! Phi's on path
    real(rp), dimension(:), pointer :: P_PATH       ! Pressure on path
    real(rp), dimension(:), pointer :: P_PATH_C     ! P_PATH on coarse grid
    real(rp), dimension(:), pointer :: PTG_ANGLES   ! (no_tan_hts)
    real(rp), dimension(:), pointer :: RADV         ! Radiances for 1 pointing on
                                                    ! Freq_Grid
    real(rp), dimension(:), pointer :: REF_CORR     ! Refraction correction
    real(rp), dimension(:), pointer :: SPS_BETA_DBETA_C ! SUM(sps*beta*dbeta)
    real(rp), dimension(:), pointer :: SPS_BETA_DBETA_F ! SUM(sps*beta*dbeta)
    real(rp), dimension(:), pointer :: STCP         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), dimension(:), pointer :: STSP         ! Sin(Theta) Sin(Phi)
    real(rp), dimension(:), pointer :: TAN_TEMP     ! ***
    real(rp), dimension(:), pointer :: TANH1_C      ! tanh(0.5 h nu / k T)
    real(rp), dimension(:), pointer :: TANH1_F      ! tanh1_c on fine grid
    real(rp), dimension(:), pointer :: TAU          ! Optical depth
    real(rp), dimension(:), pointer :: T_PATH       ! Temperatures on path
    real(rp), dimension(:), pointer :: T_PATH_C     ! T_PATH on coarse grid
    real(rp), dimension(:), pointer :: T_PATH_F     ! T_PATH on fine grid
    real(rp), dimension(:), pointer :: T_PATH_M     ! T_PATH - del_temp
    real(rp), dimension(:), pointer :: T_PATH_P     ! T_PATH + del_temp
    real(rp), dimension(:), pointer :: T_SCRIPT     ! Delta_B in some notes
    real(rp), dimension(:), pointer :: Z_GLGRID     ! Zeta on glGrid surfs
    real(rp), dimension(:), pointer :: Z_PATH       ! Zeta on path

    real(rp), dimension(:,:), pointer :: BETA_PATH_C ! Beta on path coarse
    real(rp), dimension(:),   pointer :: W0_PATH_C   ! w0 on path coarse
    real(rp), dimension(:),   pointer :: TT_PATH_C   ! tscat on path coarse
    real(r8), dimension(:,:), pointer :: VMRARRAY          ! The VMRs
    real(rp), dimension(:,:), pointer :: BETA_PATH_F ! Beta on path fine
    real(rp), dimension(:,:), pointer :: D_DELTA_DF ! Incremental opacity derivative
                                           ! schlep from drad_tran_dt to
                                           ! get_d_deltau_pol_df.  Path x SVE.
!    real(rp), dimension(:,:), pointer :: D_DELTA_DT ! Incremental opacity derivative
!                                           ! schlep from drad_tran_dt to
!                                           ! get_d_deltau_pol_dt.  Path x SVE.
    real(rp), dimension(:,:), pointer :: D_T_SCR_dT  ! D Delta_B in some notes
                                           ! path x state-vector-components
    real(rp), dimension(:,:), pointer :: D2X_DXDT    ! (No_tan_hts, nz*np)
    real(rp), dimension(:,:), pointer :: DACsStaging ! Temporary space for DACS radiances
    real(rp), dimension(:,:), pointer :: DBETA_DN_PATH_C ! dBeta_dn on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DN_PATH_F ! dBeta_dn on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DT_PATH_C ! n in beta = beta_0 (T/T_0)**n
                                           ! to compute dBeta_dt on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DT_PATH_F ! n in beta = beta_0 (T/T_0)**n
                                           ! to compute dBeta_dt on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DV_PATH_C ! dBeta_dv on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DV_PATH_F ! dBeta_dv on fine grid
    real(rp), dimension(:,:), pointer :: DBETA_DW_PATH_C ! dBeta_dw on coarse grid
    real(rp), dimension(:,:), pointer :: DBETA_DW_PATH_F ! dBeta_dw on fine grid
    real(rp), dimension(:,:), pointer :: DH_DT_PATH_C ! DH_DT_PATH on coarse grid
    real(rp), dimension(:,:), pointer :: DH_DT_PATH   ! dH/dT on path
    real(rp), dimension(:,:), pointer :: DH_DT_PATH_F ! DH_DT_PATH on fine grid
    real(rp), dimension(:,:), pointer :: DHDZ_GLGRID  ! dH/dZ on glGrid surfs
    real(rp), dimension(:,:), pointer :: DX_DT        ! (No_tan_hts, nz*np)
    real(rp), dimension(:,:), pointer :: ETA_FZP      ! Eta_z x Eta_p * Eta_f
    real(rp), dimension(:,:), pointer :: ETA_IWC      !
    real(rp), dimension(:,:), pointer :: ETA_Iwc_ZP   !
    real(rp), dimension(:,:), pointer :: ETA_Tscat    !
    real(rp), dimension(:,:), pointer :: ETA_Tscat_ZP !
    real(rp), dimension(:,:), pointer :: ETA_Salb     !
    real(rp), dimension(:,:), pointer :: ETA_Salb_ZP  !
    real(rp), dimension(:,:), pointer :: ETA_Cext     !
    real(rp), dimension(:,:), pointer :: ETA_Cext_ZP  !
    real(rp), dimension(:,:), pointer :: ETA_Mag_ZP   ! Eta_z x Eta_p
    real(rp), dimension(:,:), pointer :: ETA_ZP       ! Eta_z x Eta_p
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DN_C ! ETA_ZXP_DN on coarse grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DN   ! Eta_z x Eta_p for N
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DN_F ! ETA_ZXP_DN on fine grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DV_C ! ETA_ZXP_DV on coarse grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DV   ! Eta_z x Eta_p for V
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DV_F ! ETA_ZXP_DV on fine grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DW_C ! ETA_ZXP_DW on coarse grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DW   ! Eta_z x Eta_p for W
    real(rp), dimension(:,:), pointer :: ETA_ZXP_DW_F ! ETA_ZXP_DW on fine grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_T_C  ! ETA_ZXP_T on coarse grid
    real(rp), dimension(:,:), pointer :: ETA_ZXP_T    ! Eta_t_z x Eta_t_p
    real(rp), dimension(:,:), pointer :: ETA_ZXP_T_F  ! ETA_ZXP_T on fine grid
    real(rp), dimension(:,:), pointer :: H_GLGRID     ! H on glGrid surfs
    real(rp), dimension(:,:), pointer :: IWC_PATH     ! IWC on path
    real(rp), dimension(:,:), pointer :: Tscat_PATH   ! TScat on path
    real(rp), dimension(:,:), pointer :: TT_PATH      ! TScat on path along the LOS
    real(rp), dimension(:,:), pointer :: Salb_PATH    ! Single Scattering Albedo on path
    real(rp), dimension(:,:), pointer :: Cext_PATH    ! Cloud extinction on path
    real(rp), dimension(:,:), pointer :: K_ATMOS_FRQ  ! dI/dVMR, frq X vmr
    real(rp), dimension(:,:), pointer :: K_SPECT_DN_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_SPECT_DV_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_SPECT_DW_FRQ ! ****
    real(rp), dimension(:,:), pointer :: K_TEMP_FRQ   ! Storage for Temp. deriv.
    real(rp), dimension(:,:), pointer :: MAG_PATH     ! Magnetic field on path
    real(rp), dimension(:,:), pointer :: RADIANCES    ! (Nptg,noChans)
    real(rp), dimension(:,:), pointer :: SPS_PATH     ! species on path
    real(rp), dimension(:,:), pointer :: TAN_DH_DT    ! dH/dT at Tangent
    real(rp), dimension(:,:), pointer :: T_GLGRID     ! Temp on glGrid surfs

    real(rp), dimension(:,:,:), pointer :: DH_DT_GLGRID ! *****

    real(rp), dimension(:,:,:), pointer :: DACsStaging2 ! Temporary space for DACS radiances

    complex(rp), pointer :: D_RAD_POL_DF(:,:,:) ! From mcrt_der
    complex(rp), pointer :: D_RAD_POL_DT(:,:,:) ! From mcrt_der
    complex(rp) :: RAD_POL(2,2)  ! polarized radiance output of mcrt for one freq and pointing

    complex(rp), dimension(:,:), pointer :: ALPHA_PATH_POLARIZED
    complex(rp), dimension(:,:), pointer :: ALPHA_PATH_POLARIZED_F
      ! (-1,:,:) are Sigma_-, (0,:,:) are Pi, (+1,:,:) are Sigma_+
    complex(rp), dimension(:,:,:), pointer :: BETA_PATH_POLARIZED
    complex(rp), dimension(:,:,:), pointer :: BETA_PATH_POLARIZED_F
    complex(rp), dimension(:,:,:,:), pointer :: DE_DF    ! DE/Df in Michael's notes
    complex(rp), dimension(:,:,:,:), pointer :: DE_DT    ! DE/DT in Michael's notes
    complex(rp), dimension(:,:,:), pointer :: DELTAU_POL ! E in Michael's notes
    complex(rp), dimension(:,:,:), pointer :: DINCOPTDEPTH_POL_DT ! D Incoptdepth_Pol / DT
    complex(rp), dimension(:,:),   pointer :: GL_DELTA_POLARIZED
    complex(rp), dimension(:,:,:), pointer :: INCOPTDEPTH_POL ! 2 x 2 x path
    complex(rp), dimension(:,:,:), pointer :: INCOPTDEPTH_POL_GL ! Corrections to INCOPTDEPTH_POL
    complex(rp), dimension(:,:,:), pointer :: PROD_POL ! P in Michael's notes
    complex(rp), dimension(:,:,:), pointer :: TAU_POL  ! Tau in Michael's notes

! Some declarations by bill

    integer(ip) :: no_sv_p_t ! number of phi basis for temperature
    integer(ip) :: sps_i  ! a species counter

    real(rp) :: earthradc ! minor axis of orbit plane projected Earth ellipse
    real(rp) :: surf_angle(1), one_dhdz(1), one_dxdh(1)

    real(rp), dimension(:), pointer :: dhdz_out
    real(rp), dimension(:), pointer :: dx_dh_out
    real(rp), dimension(:), pointer :: est_scgeocalt
    real(rp), dimension(:), pointer :: req_out
    real(rp), dimension(:), pointer :: tan_chi_out
    real(rp), dimension(:), pointer :: tan_phi
    real(rp), dimension(:), pointer :: tan_press
    real(rp), dimension(:,:), pointer :: d2xdxdt_surface
    real(rp), dimension(:,:), pointer :: d2xdxdt_tan
    real(rp), dimension(:,:), pointer :: dxdt_surface
    real(rp), dimension(:,:), pointer :: dxdt_tan
    real(rp), dimension(:,:), pointer :: tan_d2h_dhdt
    real(rp), dimension(:,:,:), pointer :: ddhidhidtl0

    type (qtyStuff_t), pointer :: Qtys(:)      ! Array of pointers to Qty's.

! THIS VARIABLE REPLACES fwdModelConf%tangentGrid%surfs

    type (VectorValue_T), pointer :: BOUNDARYPRESSURE
    type (VectorValue_T), pointer :: CLOUDICE      ! Profiles
    type (VectorValue_T), pointer :: CLOUDWATER    ! Profiles
    type (VectorValue_T), pointer :: EARTHREFL     ! Earth reflectivity
    type (VectorValue_T), pointer :: ECRtoFOV      ! Rotation matrices
    type (VectorValue_T), pointer :: ELEVOFFSET    ! Elevation offset
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: GPH           ! Geopotential height
    type (VectorValue_T), pointer :: LOSVEL        ! Line of sight velocity
    type (VectorValue_T), pointer :: MAGFIELD      ! Profiles
    type (VectorValue_T), pointer :: ORBINCLINE    ! Orbital inclination
    type (VectorValue_T), pointer :: PHITAN        ! Tangent geodAngle component of state vector
    type (VectorValue_T), pointer :: PTAN          ! Tangent pressure component of state vector
    type (VectorValue_T), pointer :: REFGPH        ! Reference geopotential height
    type (VectorValue_T), pointer :: SCGEOCALT     ! S/C geocentric altitude /m
    type (VectorValue_T), pointer :: SIDEBANDFRACTION ! The sideband fraction to use
    type (VectorValue_T), pointer :: SIZEDISTRIBUTION ! Integer really
    type (VectorValue_T), pointer :: SPACERADIANCE ! Emission from space
    type (VectorValue_T), pointer :: TEMP          ! Temperature component of state vector
    type (VectorValue_T), pointer :: THISRADIANCE  ! A radiance vector quantity
    type (VectorValue_T), pointer :: VMR           ! Quantity

    type (Signal_T), pointer :: FIRSTSIGNAL        ! The first signal we're dealing with

    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS ! ***
    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS_M ! ***
    type (slabs_struct), dimension(:,:), pointer :: GL_SLABS_P ! ***

    type (catalog_T), dimension(:,:), pointer :: MY_CATALOG ! ***

    type (Grids_T) :: Grids_dn  ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_dv  ! All the spectroscopy(V) coordinates
    type (Grids_T) :: Grids_dw  ! All the spectroscopy(W) coordinates
    type (Grids_T) :: Grids_f   ! All the coordinates for VMR
    type (Grids_T) :: Grids_iwc ! All the coordinates for WC
    type (Grids_T) :: Grids_mag ! All the coordinates for Magnetic field
    type (Grids_T) :: Grids_tmp ! All the coordinates for TEMP
    type (Grids_T) :: Grids_Tscat ! All the coordinates for scaterring source function
    type (Grids_T) :: Grids_Salb ! All the coordinates for single scaterring albedo
    type (Grids_T) :: Grids_Cext ! All the coordinates for cloud extinction

!   Extra DEBUG for Nathaniel and Bill
!   real(rp), dimension(:), pointer :: REQS         ! Accumulation of REQ
!   real(rp), dimension(:), pointer :: TAN_HTS      ! Accumulation of ONE_TAN_HT
!   real(rp), dimension(:), pointer :: TAN_TEMPS    ! Accumulation of ONE_TAN_TEMP

    ! ZVI's dumping ground for variables he's too busy to put in the right
    ! place, and doesn't want to write comments for

    ! Local storage places for derivatives..(Temporary..)
    real(r4), dimension(:,:,:,:)  , pointer :: K_TEMP
    real(r4), dimension(:,:,:), pointer :: K_ATMOS
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DN
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DV
    real(r4), dimension(:,:,:,:,:,:), pointer :: K_SPECT_DW

!  The 'all_radiometers grid file' approach variables declaration:

    real(rp) :: max_ch_freq_grid, min_ch_freq_grid
    real(r8) :: TOL = 1.D-190

! *** Beta & Molecules grouping variables:
    integer, dimension(:), pointer :: mol_cat_index

    type (beta_group_T), dimension(:), pointer :: beta_group

! scatering source function for each temperature surface
    type (VectorValue_T) :: scat_src
    type (VectorValue_T) :: scat_alb
    type (VectorValue_T) :: cld_ext

! Channel information from the signals database
    type :: Channels_T
      integer :: Used       ! Which channel is this?
      integer :: Origin     ! Index of first channel (zero or one)
      integer :: Signal     ! Signal index for the channel
      integer :: DACS       ! DACS index if any, else zero
    end type Channels_T

    type(channels_T), allocatable, dimension(:) :: Channels 
     
    ! Executable code --------------------------------------------------------
    ! ------------------------------------------------------------------------
!   Print *, '** Enter ForwardModel, MAF =',fmstat%maf   ! ** ZEBUG

    if ( e_stop > 0.0_rp ) e_stop = log(epsilon(0.0_rp)) ! only once

    temp_der = present ( jacobian ) .and. FwdModelConf%temp_der
    atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der

! ** Re-instate when appropriate code is done
!   spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
    spect_der = present ( jacobian ) .and. .false.    ! ** ZEBUG

    any_der = temp_der .or. atmos_der .or. spect_der

    if ( toggle(emit) ) & ! set by -f command-line switch
      & call trace_begin ( 'ForwardModel, MAF=', index=fmstat%maf )

    ! Nullify all our pointers that are allocated!

    nullify ( alpha_path_c, alpha_path_f, alpha_path_polarized, &
      & alpha_path_polarized_f, beta_path_c, beta_path_cloud_c, beta_path_f, &
      & beta_path_polarized, cext_path, cg_inds, c_inds, &
      & cld_ext%values, closestInstances, d2x_dxdt, d2xdxdt_surface, &
      & d2xdxdt_tan, DACsStaging, DACsStaging2, dbeta_dn_path_c, &
      & dbeta_dn_path_f, dbeta_dt_path_c, dbeta_dt_path_f, dbeta_dv_path_c, &
      & dbeta_dv_path_f, dbeta_dw_path_c, dbeta_dw_path_f, d_delta_df, &
      & ddhidhidtl0, de_df, de_dt, del_s, deltau_pol, del_zeta, &
      & dh_dt_glgrid, dh_dt_path, dh_dt_path_c, dh_dt_path_f, dhdz_glgrid, &
      & dhdz_gw_path, dhdz_out, dhdz_path, dincoptdepth_pol_dt, &
      & do_calc_Cext, do_calc_Cext_zp, do_calc_dn, do_calc_dn_c, &
      & do_calc_dn_F, do_calc_dv, do_calc_dv_c, do_calc_dv_f, do_calc_dw, &
      & do_calc_dw_c, do_calc_dw_f, do_calc_fzp, do_calc_hyd, do_calc_hyd_c, &
      & do_calc_iwc, do_calc_iwc_zp, do_calc_Salb, do_calc_Salb_zp, &
      & do_calc_t, do_calc_t_c, do_calc_t_f, do_calc_tscat, &
      & do_calc_tscat_zp, do_calc_zp, do_gl, drad_dn, drad_dt, drad_dv, &
      & drad_dw, d_rad_pol_df, d_rad_pol_dt, dsdz_gw_path, d_t_scr_dt, &
      & dx_dh_out, dx_dt, dxdt_surface, dxdt_tan, eta_cext, eta_cext_zp, &
      & eta_fzp, eta_iwc, eta_iwc_zp, eta_mag_zp, eta_salb, eta_salb_zp, &
      & eta_tscat, eta_tscat_zp, eta_zp, eta_zxp_dn, eta_zxp_dn_c, &
      & eta_zxp_dn_f, eta_zxp_dv, eta_zxp_dv_c, eta_zxp_dv_f, eta_zxp_dw, &
      & eta_zxp_dw_c, eta_zxp_dw_f, eta_zxp_t, eta_zxp_t_c, eta_zxp_t_f, &
      & f_inds, frequencies, gl_delta_polarized, gl_inds, gph, grids, &
      & h_glgrid, h_path, h_path_c, h_path_f, incoptdepth, incoptdepth_pol, &
      & incoptdepth_pol_gl, ipsd, iwc_path, k_atmos, k_atmos_frq, &
      & k_spect_dn, k_spect_dn_frq, k_spect_dv, k_spect_dv_frq, k_spect_dw, &
      & k_spect_dw_frq, k_temp, k_temp_frq, mag_path, mol_cat_index, n_path, &
      & path_dsdh, phi_path, p_path, p_path_c, prod_pol, ptg_angles, &
      & radiances, RadV, ref_corr, req_out, salb_path, scat_alb%values, &
      & scat_ang, scat_src%values, sps_beta_dbeta_c, sps_beta_dbeta_f, &
      & sps_path, tan_chi_out, tan_d2h_dhdt, tan_dh_dt, tanh1_c, tanh1_f, &
      & tan_phi, tan_temp, tau, tau_pol, t_der_path_flags, t_glgrid, t_path, &
      & t_path_c, t_path_f, t_path_m, t_path_p, true_path_flags, tscat_path, &
      & t_script, tt_path, tt_path_c, &
      & usedDacsSignals, vmr, vmrarray, w0_path_c, wc, z_path )

    ! Extra DEBUG for Nathaniel and Bill
!   nullify ( reqs, tan_hts, tan_temps )

    ! Work out what we've been asked to do -----------------------------------

    ! Identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    firstSignal => fwdModelConf%signals(1) ! Config has verified that signals
      ! are all for same radiometer, module and sideband

    ! Create the data structures for the species

    call get_species_data ( fwdModelConf, fwdModelIn, fwdModelExtra, &
      & noSpecies, no_mol, beta_group, my_catalog )

    call allocate_test ( mol_cat_index, no_mol, 'mol_cat_index', moduleName )

    mol_cat_index = PACK((/(i,i=1,noSpecies)/), fwdModelConf%molecules > 0)

    allocate ( qtys(no_mol), stat = ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'QTYS' )

    ! Get state vector quantities for species
    do sps_i = 1 , no_mol
      qtys(sps_i)%qty => GetQuantityForForwardModel(fwdmodelin, fwdmodelextra, &
        &  quantitytype=l_vmr, molIndex=mol_cat_index(sps_i), config=fwdModelConf, &
        &  radiometer=firstsignal%radiometer, foundInFirst=qtys(sps_i)%foundInFirst )
    end do

    ! Start sorting out stuff from state vector ------------------------------

    ! Identify the appropriate state vector components, save vmrs for later
    gph => GetVectorQuantityByType ( fwdModelExtra, quantityType=l_gph )
    temp => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature, config=fwdModelConf )
    ptan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=firstSignal%instrumentModule, &
      & foundInFirst=ptan_der, config=fwdModelConf )
    phitan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, instrumentModule=firstSignal%instrumentModule, &
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
    if ( FwdModelConf%polarized ) then
      magField => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_magneticField, config=fwdModelConf )
      ECRtoFOV => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_ECRtoFOV, config=fwdModelConf )
    end if
    if ( FwdModelConf%incl_cld ) then
      cloudIce => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra,  &
        & quantityType=l_cloudIce, noError=.true., config=fwdModelConf )
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

    MAF = fmStat%maf

    Vel_Cor = 1.0_rp - losvel%values(1,maf) / speedOfLight

! Sort out a remaining flag
    ptan_der = ptan_der .and. present ( jacobian )

! Create the Grids_tmp structure:

    !??? Should .true. be temp_der ???
    call load_one_item_grid ( grids_tmp, temp, phitan, maf, fwdModelConf, .true. )

    windowStart = grids_tmp%windowStart(1)
    windowFinish = grids_tmp%windowFinish(1)
    no_sv_p_t = windowFinish - windowStart + 1
    n_t_zeta = grids_tmp%l_z(1)
    sv_t_len = grids_tmp%p_len ! zeta X phi

    if ( FwdModelConf%incl_cld ) &
      & call load_one_item_grid ( grids_iwc, cloudIce, phitan, maf, &
        & fwdModelConf, .false., .false. )

    if ( FwdModelConf%polarized ) &
      & call load_one_item_grid ( grids_mag, magfield, phitan, maf, &
        & fwdModelConf, .false. )

    ! Identify which of our signals are DACS and how many unique DACS are involved
    call identifyDACS

    ! Work out which channels are used; also check we have radiances for them.
    noUsedChannels = 0
    do sigInd = 1, size(fwdModelConf%signals)
      ! This just emits an error message and stops if we don't have a radiance.
      ! We don't use the vector quantity -- at least not right away.  We get
      ! it again later.
      thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
        & signal=fwdModelConf%signals(sigInd)%index, &
        & sideband=merge ( 0, firstSignal%sideband, fwdModelConf%forceFoldedOutput ) )
      noUsedChannels = noUsedChannels + &
        & count( fwdModelConf%signals(sigInd)%channels )
    end do
    allocate ( channels(noUsedChannels), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'channels' )

    ! Collect channel information from signals database.
    channel = 1
    do sigInd = 1, size(fwdModelConf%signals)
      do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
        if ( fwdModelConf%signals(sigInd)%channels(i) ) then
          channels(channel)%origin = &
            & lbound ( fwdModelConf%signals(sigInd)%frequencies, 1 )
          channels(channel)%used = i + channels(channel)%origin - 1
          channels(channel)%signal = sigInd
          channels(channel)%dacs = FindFirst ( sigind == usedDACSSignals )
          channel = channel + 1
        end if
      end do
    end do

    maxNoFFreqs = 0
    maxNoFSurfs = 0
    do sps_i = 1 , no_mol
      maxNoFFreqs = max(maxNoFFreqs, qtys(sps_i)%qty%template%noChans)
      maxNoFSurfs = max(maxNoFSurfs, qtys(sps_i)%qty%template%noSurfs)
    end do

! Now, allocate other variables we're going to need later ----------------

    ! Stuff for clouds
    if ( FwdModelConf%incl_cld ) then  ! Do this block only if incl_cld is true

      call Allocate_test ( closestInstances, thisRadiance%template%noInstances,    &
      & 'closestInstances', ModuleName )
      call FindClosestInstances ( temp, thisRadiance, closestInstances )
      inst = closestInstances(MAF)

      call allocate_test ( scat_src%values, n_t_zeta, &
      & fwdModelConf%num_scattering_angles, 'scat_src', moduleName )

      call allocate_test ( scat_alb%values, n_t_zeta, 2, 'scat_alb', moduleName )
      call allocate_test (  cld_ext%values, n_t_zeta, 2, 'cld_ext', moduleName )

      if ( size(FwdModelConf%molecules) .lt. 2 ) then
         !   make sure we have enough molecules
         call MLSMessage ( MLSMSG_Error, ModuleName, 'Not enough molecules' )
      end if

      ! now checking spectroscopy
      got = .false.
      nspec = size(FwdModelConf%molecules)
      noSurf  = temp%template%noSurfs
      call allocate_test ( vmrArray, nspec, n_t_zeta, 'vmrArray', ModuleName )
      vmrarray = 0.0_r8

      do j = 1, nspec      ! Loop over species
        call get_string ( lit_indices( abs(fwdModelConf%molecules(j)) ), molName )

        select case (fwdModelConf%molecules(j))
        case ( L_H2O, L_O3, L_N2O )
          ispec = 0
          if(fwdModelConf%molecules(j) == l_h2o) ispec = 1
          if(fwdModelConf%molecules(j) == l_o3) ispec = 2
          if(fwdModelConf%molecules(j) == l_n2o) ispec = 3
          if(fwdModelConf%molecules(j) == l_h2o) got(L_H2O) = .true.
          if(fwdModelConf%molecules(j) == l_o3) got(L_O3) = .true.

          vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,            &
            & quantityType=l_vmr, molecule=fwdModelConf%molecules(j) )

          novmrSurf = vmr%template%nosurfs

          call InterpolateValues ( &
          & reshape(vmr%template%surfs(:,1),(/novmrSurf/)), &    ! Old X
          & reshape(vmr%values(:,inst),(/novmrSurf/)),      &    ! Old Y
          & reshape(temp%template%surfs(:,1),(/noSurf/)),   &    ! New X
          & vmrArray(ispec,:),                              &    ! New Y
          & 'Linear', extrapolate='Clamp' )

        case ( L_N2, L_O2, L_H2O_18, L_O_18_O)
          call MLSMessage(MLSMSG_Warning, ModuleName, &
          &'cloud fwd model internally has this molecule: '//trim(molName))
        case default
          call MLSMessage(MLSMSG_Error, ModuleName, &
          &'cloud fwd model currently cannot accept this molecule: '//trim(molName))
        end select
      end do ! End of Loop over species

      if ( .not. got(l_h2o) .or. .not. got(l_o3) ) then
      !make sure we have at least two molecules h2o and o3.
        call MLSMessage(MLSMSG_Error, ModuleName,'Missing molecules H2O or O3 in cloud FM')
      end if

    end if

! Set up our temporary `state vector' like arrays ------------------------

    call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
      & firstSignal%radiometer, mol_cat_index, l_vmr, &
      & grids_f, h2o_ind, ext_ind )
    if ( spect_der ) then
      call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
        & firstSignal%radiometer, mol_cat_index, l_dn, grids_dn )
      call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
        & firstSignal%radiometer, mol_cat_index, l_dv, grids_dv )
      call load_sps_data ( FwdModelConf,  FwdModelIn, FwdModelExtra, FmStat, &
        & firstSignal%radiometer, mol_cat_index, l_dw, grids_dw )
    end if

! modify h2o mixing ratio if a special supersaturation is requested
    if ( fwdModelConf%i_saturation /= l_clear ) then
      call modify_values_for_supersat ( fwdModelConf, grids_f, h2o_ind, &
        & grids_tmp, boundaryPressure )
    end if

! set up output pointing angles ------------------------------------------
! note we have to compute req !!!!!!!

! compute equivalent earth radius at phi_t(1), nearest surface

    call allocate_test ( req_out, ptan%template%nosurfs, 'req_out', &
                       & moduleName )
    earthradc = earthRadA*earthRadB / &
          & SQRT((earthRadA**2-earthRadB**2) * &
                 &   SIN(orbIncline%values(1,maf)*Deg2Rad)**2 + &
                 & earthRadB**2)

!{\begin{equation*}\begin{split}
! R_{eq} =\;& \sqrt \frac{R_a^4 \sin^2 \phi + R_c^4 \cos^2 \phi}
!                       {R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}\\
!        =\;& \sqrt \frac{R_a^4 - (R_a^2+R_c^2)(R_a^2-R_c^2) \cos^2 \phi}
!                       {R_c^2 +              (R_a^2-R_c^2) \cos^2 \phi}
! \end{split}\end{equation*}

    req_out = (earthrada-earthradc)*(earthrada+earthradc) * &
      & COS(phitan%values(:,maf)*Deg2Rad)**2
    req_out = 0.001_rp * SQRT( &
      & ( earthrada**4 - (earthrada**2+earthradc**2) * req_out ) / &
      & ( earthradc**2 + req_out ) )

    call allocate_test ( tan_chi_out, ptan%template%nosurfs, 'tan_chi_out', &
                       & moduleName )
    call allocate_test ( dx_dh_out, ptan%template%nosurfs, 'dx_dh_out', &
                       & moduleName )
    call allocate_test ( dhdz_out, ptan%template%nosurfs, 'dhdz_out', &
                       & moduleName )

    if ( temp_der ) then
      call allocate_test ( dxdt_tan, ptan%template%nosurfs, sv_t_len, &
                         & 'dxdt_tan',moduleName )
      call allocate_test ( d2xdxdt_tan, ptan%template%nosurfs, sv_t_len, &
                         & 'd2xdxdt_tan',moduleName )
    end if

    ! Temperature's windowStart:windowFinish are correct here.
    ! RefGPH and Temperature have the same horizontal basis.
    call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad, &
       & 0.001_rp*scGeocAlt%values(:,maf), Grids_tmp, &
       & (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
       & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
       & orbIncline%values(1,maf)*Deg2Rad, 0.0_rp, &
       & req_out, grids_f, h2o_ind, tan_chi_out, dhdz_out, dx_dh_out, &
       & dxdt_tan=dxdt_tan, d2xdxdt_tan=d2xdxdt_tan )

! Compute Gauss Legendre (GL) grid ---------------------------------------

    call compute_GL_grid ( fwdModelConf, temp, qtys, nlvl, maxVert, &
      &                    p_glgrid, z_glgrid, tan_inds, tan_press )

! Allocate more GL grid stuff

    call allocate_test ( h_glgrid,    maxVert, no_sv_p_t, 'h_glgrid', moduleName )
    call allocate_test ( t_glgrid,    maxVert, no_sv_p_t, 't_glgrid', moduleName )
    call allocate_test ( dhdz_glgrid, maxVert, no_sv_p_t, 'dhdz_glgrid', &
                      &  moduleName )
    call allocate_test ( dh_dt_glgrid, maxVert,  n_t_zeta, no_sv_p_t, &
                      &  'dh_dt_glgrid', moduleName )
    call allocate_test ( ddhidhidtl0, maxVert, n_t_zeta, no_sv_p_t, &
                      &  'ddhidhidtl0', moduleName )

    surfaceTangentIndex = COUNT(tan_inds == 1)

    no_tan_hts = size(tan_inds)

    ! Extra DEBUG for Nathaniel and Bill
!   call allocate_test ( tan_hts,       no_tan_hts, 'tan_hts',       moduleName )
!   call allocate_test ( tan_temps,     no_tan_hts, 'tan_temps',     moduleName )
!   call allocate_test ( reqs,          no_tan_hts, 'reqs',          moduleName )

    ! estimate tan_phi and scgeocalt
    call estimate_tan_phi ( no_tan_hts, nlvl, maf, phitan, ptan, &
                          & scgeocalt, tan_press, &
                          & tan_phi, est_scgeocalt )


    ! Compute hydrostatic grid -----------------------------------------------

    ! Insert into bill's 2d hydrostatic equation.
    ! The phi input for this program are the orbit plane projected
    ! geodetic locations of the temperature phi basis--not necessarily
    ! the tangent phi's which may be somewhat different.

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_Begin ( 'ForwardModel.Hydrostatic' )

    ! Temperature's windowStart:windowFinish are correct here.
    ! RefGPH and temperature have the same horizontal basis.
    call two_d_hydrostatic ( Grids_tmp, &
      &  (/ (refGPH%template%surfs(1,1), j=1,no_sv_p_t) /), &
      &  0.001*refGPH%values(1,windowStart:windowFinish), z_glgrid, &
      &  orbIncline%values(1,maf)*Deg2Rad, t_glgrid, h_glgrid, &
      &  dhdz_glgrid, dh_dt_glgrid, DDHDHDTL0=ddhidhidtl0 )

! This is a lazy way to get the surface angle

    if ( temp_der ) then
      call allocate_test ( dxdt_surface, 1, sv_t_len,  'dxdt_surface', moduleName )
      call allocate_test ( d2xdxdt_surface, 1, sv_t_len,  'd2xdxdt_surface', &
                       & moduleName )
    end if

    ! Temperature's windowStart:windowFinish are correct here.
    ! refGPH and temperature have the same horizontal basis.
    call get_chi_out ( tan_press(1:1), tan_phi(1:1), &
       & 0.001_rp*est_scgeocalt(1:1), Grids_tmp, &
       & SPREAD(refGPH%template%surfs(1,1),1,no_sv_p_t), &
       & 0.001_rp*refGPH%values(1,windowStart:windowFinish), &
       & orbIncline%values(1,maf)*Deg2Rad, 0.0_rp, &
       & (/req_out(1)/), grids_f, h2o_ind, surf_angle, one_dhdz, one_dxdh,  &
       & dxdt_tan=dxdt_surface, d2xdxdt_tan=d2xdxdt_surface )

    call deallocate_test ( d2xdxdt_surface, 'd2xdxdt_surface',moduleName )
    call deallocate_test ( req_out, 'req_out', moduleName )

 ! Now, allocate other variables we're going to need later ----------------

    if ( toggle(emit) .and. levels(emit) > 0 ) then
      call Trace_End ( 'ForwardModel.Hydrostatic' )
      call Trace_Begin ( 'ForwardModel.Allocate' )
    end if

    call allocate_test ( tan_temp, no_tan_hts, 'tan_temp', moduleName )

    ! Allocate path quantities -----------------------------------------------

    max_ele = 2*maxVert     ! maximum possible

    ! Now allocate all path related... with maximum length..

    brkpt = maxVert
    npc = 2 * (brkpt + Ng) / Ngp1
    ! MJS says that this is the same as npc = 2 * Nlvl

    ! This can be put outside the mmaf loop

    call allocate_test ( alpha_path_c,        npc, 'alpha_path_c',     moduleName )
    call allocate_test ( alpha_path_f,    max_ele, 'alpha_path_f',     moduleName )
    call allocate_test ( beta_path_cloud_c,   npc, 'beta_path_cloud_c',moduleName )
    call allocate_test ( w0_path_c,           npc, 'w0_path_c',        moduleName )
    call allocate_test ( tt_path_c,           npc, 'tt_path_c',        moduleName )
    call allocate_test ( c_inds,              npc, 'c_inds',           moduleName )
    call allocate_test ( cg_inds,             npc, 'cg_inds',          moduleName )
    call allocate_test ( del_s,               npc, 'del_s',            moduleName )
    call allocate_test ( del_zeta,            npc, 'del_zeta',         moduleName )
    call allocate_test ( dhdz_path,       max_ele, 'dhdz_path',        moduleName )
    call allocate_test ( dhdz_gw_path,    max_ele, 'dhdz_gw_path',     moduleName )
    call allocate_test ( do_gl,               npc, 'do_gl',            moduleName )
    call allocate_test ( dsdz_gw_path,    max_ele, 'dsdz_gw_path',     moduleName )
    call allocate_test ( f_inds,         ng * npc, 'f_inds',           moduleName )
    call allocate_test ( gl_inds,         max_ele, 'gl_inds',          moduleName )
    call allocate_test ( h_path_c,            npc, 'h_path_c',         moduleName )
    call allocate_test ( h_path_f,        max_ele, 'h_path_f',         moduleName )
    call allocate_test ( h_path,          max_ele, 'h_path',           moduleName )
    call allocate_test ( incoptdepth,         npc, 'incoptdept',       moduleName )
    call allocate_test ( n_path,              npc, 'n_path',           moduleName )
    call allocate_test ( path_dsdh,       max_ele, 'path_dsdh',        moduleName )
    call allocate_test ( phi_path,        max_ele, 'phi_path',         moduleName )
    call allocate_test ( p_path_c,            npc, 'p_path_c',         moduleName )
    call allocate_test ( p_path,          max_ele, 'p_path',           moduleName )
    call allocate_test ( ref_corr,            npc, 'ref_corr',         moduleName )
    call allocate_test ( sps_beta_dbeta_c,    npc, 'sps_beta_dbeta_c', moduleName )
    call allocate_test ( sps_beta_dbeta_f, max_ele, 'sps_beta_dbeta_f', moduleName )
    call allocate_test ( tanh1_c,             npc, 'tanh1_c',          moduleName )
    call allocate_test ( tanh1_f,         max_ele, 'tanh1_f',          moduleName )
    call allocate_test ( tau,                 npc, 'tau',              moduleName )
    call allocate_test ( t_path_c,            npc, 't_path_c',         moduleName )
    call allocate_test ( t_path_f,        max_ele, 't_path_f',         moduleName )
    call allocate_test ( t_path_m,        max_ele, 't_path_m',         moduleName )
    call allocate_test ( t_path_p,        max_ele, 't_path_p',         moduleName )
    call allocate_test ( t_path,          max_ele, 't_path',           moduleName )
    call allocate_test ( t_script,            npc, 't_script',         moduleName )
    call allocate_test ( z_path,          max_ele, 'z_path',           moduleName )

    call allocate_test ( beta_path_c,      npc, no_mol, 'beta_path_c',   moduleName )
    call allocate_test ( beta_path_f,   max_ele, no_mol, 'beta_path_f',   moduleName )
    call allocate_test ( do_calc_fzp,   max_ele, size(grids_f%values),  'do_calc_fzp',   moduleName )
    call allocate_test ( do_calc_zp,    max_ele, grids_f%p_len,  'do_calc_zp',     moduleName )
    call allocate_test ( eta_fzp,       max_ele, size(grids_f%values),  'eta_fzp', moduleName )
    call allocate_test ( eta_zp,        max_ele, grids_f%p_len,  'eta_zp', moduleName )
    call allocate_test ( sps_path,      max_ele, no_mol, 'sps_path',       moduleName )
    call allocate_test ( true_path_flags, max_ele, 'true_path_flags',moduleName)
    true_path_flags = .true.
    if ( fwdModelConf%Incl_Cld ) then
      call allocate_test ( do_calc_iwc,    max_ele, size(grids_iwc%values),  'do_calc_iwc',  moduleName )
      call allocate_test ( eta_iwc,        max_ele, size(grids_iwc%values),  'eta_iwc',      moduleName )
      call allocate_test ( eta_iwc_zp,     max_ele, grids_iwc%p_len,  'eta_iwc_zp',    moduleName )
      call allocate_test ( iwc_path,       max_ele, 1, 'iwc_path',           moduleName )
      call allocate_test ( ipsd,           max_ele, 'IPSD', moduleName )
      call allocate_test ( wc,         fwdModelConf%no_cloud_species, max_ele, 'WC', moduleName )
    end if
    if ( temp_der ) then

! Allocation for metrics routine when Temp. derivative is needed:

!      call allocate_test ( k_temp, noUsedChannels, no_tan_hts, n_t_zeta, &
!                         & no_sv_p_t, 'k_temp',moduleName )
      allocate ( k_temp(noUsedChannels, no_tan_hts, n_t_zeta, no_sv_p_t) )

      call allocate_test ( dRad_dt, sv_t_len, 'dRad_dt', moduleName )
      call allocate_test ( d_t_scr_dt,         npc, sv_t_len, 'd_t_scr_dt', &
                                                              & moduleName )
      call allocate_test ( dbeta_dt_path_c,    npc, no_mol,   'dbeta_dt_path_c', &
                                                              & moduleName )
      call allocate_test ( dbeta_dt_path_f, max_ele, no_mol,   'dbeta_dt_path_f', &
                                                              & moduleName )
      call allocate_test ( dh_dt_path,      max_ele, sv_t_len, 'dh_dt_path', &
                                                              & moduleName )
      call allocate_test ( dh_dt_path_c,       npc, sv_t_len, 'dh_dt_path_c', &
                                                              & moduleName )
      call allocate_test ( dh_dt_path_f,    max_ele, sv_t_len, 'dh_dt_path_f', &
                                                              & moduleName )
      call allocate_test ( do_calc_hyd,     max_ele, sv_t_len, 'do_calc_hyd', &
                                                              & moduleName )
      call allocate_test ( do_calc_hyd_c,      npc, sv_t_len, 'do_calc_hyd_c', &
                                                              & moduleName )
      call allocate_test ( do_calc_t,       max_ele, sv_t_len, 'do_calc_t', &
                                                              & moduleName )
      call allocate_test ( do_calc_t_c,        npc, sv_t_len, 'do_calc_t_c', &
                                                              & moduleName )
      call allocate_test ( do_calc_t_f,     max_ele, sv_t_len, 'do_calc_t_f', &
                                                              & moduleName )
      call allocate_test ( eta_zxp_t,       max_ele, sv_t_len, 'eta_zxp_t', &
                                                              & moduleName )
      call allocate_test ( eta_zxp_t_c,        npc, sv_t_len, 'eta_zxp_t_c', &
                                                              & moduleName )
      call allocate_test ( eta_zxp_t_f,     max_ele, sv_t_len, 'eta_zxp_t_f', &
                                                              & moduleName )
      call allocate_test ( tan_dh_dt,    1, sv_t_len, 'tan_dh_dt',    moduleName )
      call allocate_test ( tan_d2h_dhdt, 1, sv_t_len, 'tan_d2h_dhdt', moduleName )
      call allocate_test ( t_der_path_flags, max_ele,          't_der_path_flags', &
                                                              & moduleName )

    end if ! temp_der

    if ( atmos_der ) then
      call allocate_test ( d_delta_df, npc, size(grids_f%values), 'd_delta_df', &
                                                              & moduleName )
      call allocate_test ( k_atmos, noUsedChannels, no_tan_hts, size(grids_f%values), &
        & 'k_atmos', moduleName )
      k_atmos = 0.0
    end if ! atmos_der

    if ( spect_der ) then

      !??? Are temperature's windowStart:windowFinish correct here.  Should ???
      !??? we use minval(grids_d?%windowStart):maxval(grids_d?%windowFinish) ???
      allocate ( k_spect_dw(noUsedChannels, no_tan_hts, maxNoFFreqs, &
        & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'k_spect_dw' )
      allocate ( k_spect_dn(noUsedChannels, no_tan_hts, maxNoFFreqs, &
        & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'k_spect_dn' )
      allocate ( k_spect_dv(noUsedChannels, no_tan_hts, maxNoFFreqs, &
        & maxNoFSurfs, windowStart:windowFinish, no_mol), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'k_spect_dv' )

      call allocate_test ( dbeta_dw_path_c,    npc, no_mol, &
        & 'dbeta_dw_path_c', moduleName )
      call allocate_test ( dbeta_dw_path_f, max_ele, no_mol, &
        & 'dbeta_dw_path_f', moduleName )
      call allocate_test ( dbeta_dn_path_c,    npc, no_mol, &
        & 'dbeta_dn_path_c', moduleName )
      call allocate_test ( dbeta_dn_path_f, max_ele, no_mol, &
        & 'dbeta_dn_path_f', moduleName )
      call allocate_test ( dbeta_dv_path_c,    npc, no_mol, &
        & 'dbeta_dv_path_c', moduleName )
      call allocate_test ( dbeta_dv_path_f, max_ele, no_mol, &
        & 'dbeta_dv_path_f', moduleName )

      f_len_dw = grids_dw%l_v(ubound(grids_dw%l_v,1))
      f_len_dn = grids_dn%l_v(ubound(grids_dn%l_v,1))
      f_len_dv = grids_dv%l_v(ubound(grids_dv%l_v,1))

      call allocate_test ( do_calc_dw, max_ele, f_len_dw, 'do_calc_dw', &
                        &  moduleName )
      call allocate_test ( do_calc_dw_c, max_ele, f_len_dw, 'do_calc_dw_c', &
                        &  moduleName )
      call allocate_test ( do_calc_dw_f, max_ele, f_len_dw, 'do_calc_dw_f', &
                        &  moduleName )
      call allocate_test ( do_calc_dn, max_ele, f_len_dn, 'do_calc_dn', &
                        &  moduleName )
      call allocate_test ( do_calc_dn_c, max_ele, f_len_dn, 'do_calc_dn_c', &
                        &  moduleName )
      call allocate_test ( do_calc_dn_f, max_ele, f_len_dn, 'do_calc_dn_f', &
                        &  moduleName )
      call allocate_test ( do_calc_dv, max_ele, f_len_dv, 'do_calc_dv', &
                        &  moduleName )
      call allocate_test ( do_calc_dv_c, max_ele, f_len_dv, 'do_calc_dv_c', &
                        &  moduleName )
      call allocate_test ( do_calc_dv_f, max_ele, f_len_dv, 'do_calc_dv_f', &
                        &  moduleName )

      call allocate_test ( eta_zxp_dw, max_ele, f_len_dw, 'eta_zxp_dw', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dw_c, max_ele, f_len_dw, 'eta_zxp_dw_c', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dw_f, max_ele, f_len_dw, 'eta_zxp_dw_f', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dn, max_ele, f_len_dn, 'eta_zxp_dn', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dn_c, max_ele, f_len_dn, 'eta_zxp_dn_c', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dn_f, max_ele, f_len_dn, 'eta_zxp_dn_f', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dv, max_ele, f_len_dv, 'eta_zxp_dv', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dv_c, max_ele, f_len_dv, 'eta_zxp_dv_c', &
                        &  moduleName )
      call allocate_test ( eta_zxp_dv_f, max_ele, f_len_dv, 'eta_zxp_dv_f', &
                        &  moduleName )

      call allocate_test ( drad_dw, f_len_dw, 'drad_dw', moduleName )
      call allocate_test ( drad_dn, f_len_dn, 'drad_dn', moduleName )
      call allocate_test ( drad_dv, f_len_dv, 'drad_dv', moduleName )

    else

      f_len_dw = 0
      f_len_dn = 0
      f_len_dv = 0

    end if

    if ( FwdModelConf%polarized ) then
      call allocate_test ( eta_mag_zp,     max_ele, grids_mag%p_len,        'eta_mag_zp',     moduleName )
      call allocate_test ( mag_path,       max_ele, magfield%template%noChans+1, 'mag_path', &
        & moduleName )
      allocate ( alpha_path_polarized(-1:1,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'alpha_path_polarized' )
      allocate ( beta_path_polarized(-1:1,npc,no_mol), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'beta_path_polarized' )
      allocate ( alpha_path_polarized_f(-1:1,max_ele), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'alpha_path_polarized_f' )
      allocate ( beta_path_polarized_f(-1:1,max_ele,no_mol), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'beta_path_polarized_f' )
      allocate ( gl_delta_polarized(-1:1,max_ele), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'gl_delta_polarized' )
      allocate ( incoptdepth_pol(2,2,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'incoptdepth_pol' )
      allocate ( incoptdepth_pol_gl(2,2,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'incoptdepth_pol_gl' )
      allocate ( prod_pol(2,2,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'prod_pol' )
      allocate ( tau_pol(2,2,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'tau_pol' )
      allocate ( deltau_pol(2,2,npc), stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'deltau_pol' )
      if ( atmos_der ) then
        allocate ( d_rad_pol_df(2,2,size(grids_f%values)), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'d_rad_pol_df' )
        allocate ( de_df(2,2,npc,size(grids_f%values)), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'de_df' )
      end if
      if ( temp_der ) then
        allocate ( d_rad_pol_dt(2,2,sv_t_len), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'d_rad_pol_dt' )
        allocate ( de_dt(2,2,npc,sv_t_len), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'de_dt' )
        allocate ( dincoptdepth_pol_dt(2,2,npc), stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate//'dincoptdepth_pol_dt' )
      end if
    end if

    call allocate_test ( ptg_angles,no_tan_hts, 'ptg_angles',moduleName )
    call allocate_test ( radiances, no_tan_hts, noUsedChannels, &
                       & 'radiances', moduleName )

    call allocate_test ( dx_dt, no_tan_hts,sv_t_len, 'dx_dt',moduleName )
    call allocate_test ( d2x_dxdt,no_tan_hts,sv_t_len, 'd2x_dxdt',moduleName )

    if ( toggle(emit) .and. levels(emit) > 0 ) then
      call Trace_End ( 'ForwardModel.Allocate' )
      call Trace_Begin ( 'ForwardModel.SidebandLoop' )
    end if

    ! Loop over sidebands ----------------------------------------------------
    do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call Trace_Begin ( 'ForwardModel.Sideband ', index=thisSideband )

      ! Now, allocate gl_slabs arrays
      call allocateSlabs ( gl_slabs, max_ele, my_catalog(thisSideband,:), moduleName )
      if ( temp_der ) then
        call allocateSlabs ( gl_slabs_p, max_ele, my_catalog(thisSideband,:), moduleName )
        call allocateSlabs ( gl_slabs_m, max_ele, my_catalog(thisSideband,:), moduleName )
      end if

      ! Work out which pointing frequency grid we're going to need if ----------
      ! frequency averaging

      ! Code splits into two sections, one for when we're doing frequency
      ! averaging, and one when we're not.
      if ( fwdModelConf%do_freq_avg ) then ! --- Doing freq. avg. ---

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
          sigInd = channels(i)%signal
          channel = channels(i)%used
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
          direction = fwdModelConf%signals(channels(channel)%signal)%direction
          frequencies(channel) = &
            & fwdModelConf%signals(channels(channel)%signal)%centerFrequency + &
            & direction * fwdModelConf%signals(channels(channel)%signal)% &
            & frequencies(channels(channel)%used)
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
                             & moduleName )

      if ( atmos_der ) &
        & call allocate_test ( k_atmos_frq, maxNoPtgFreqs, size(grids_f%values), 'k_atmos_frq',&
                             & moduleName )

      if ( spect_der ) then
        call allocate_test ( k_spect_dw_frq , maxNoPtgFreqs, f_len_dw, &
                           & 'k_spect_dw_frq', moduleName )
        call allocate_test ( k_spect_dn_frq , maxNoPtgFreqs, f_len_dn, &
                           & 'k_spect_dn_frq', moduleName )
        call allocate_test ( k_spect_dv_frq , maxNoPtgFreqs, f_len_dv, &
                           & 'k_spect_dv_frq', moduleName )
      end if

      if ( any_der ) &
        & call allocate_test ( DACsStaging2, ubound(DACsStaging,1), &
                             &               max(sv_t_len,size(grids_f%values), &
                                             &   f_len_dw, f_len_dn, f_len_dv), &
                             &               noUsedDACs, &
                             & 'DACsStaging2', moduleName, &
                             & low1=lbound(DACsStaging,1) )

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call Trace_Begin ( 'ForwardModel.PointingLoop' )

      ! Loop over pointings --------------------------------------------------
      do ptg_i = 1, no_tan_hts
        if ( toggle(emit) .and. levels(emit) > 3 ) &
          & call Trace_Begin ( 'ForwardModel.Pointing ', index=ptg_i )

        if ( FwdModelConf%polarized ) &
          & mif = minloc(abs(tan_press(ptg_i) - &
          &                  ptan%values(:ptan%template%nosurfs,maf)),1)

        ! allocate the path stuff
        brkpt = MaxVert + 1 - tan_inds(ptg_i) ! path tangent index
        no_ele = 2 * brkpt
        mid = (brkpt + Ng) / Ngp1
        npc = 2 * mid

        ! This is not pretty but we need some coarse grid extraction indices
        c_inds(1:npc) = (/(i*Ngp1-Ng,i=1,mid),((i-1)*Ngp1-Ng+1,i=mid+1,npc)/)
        c_inds(npc+1:) = 0
        ! And some fine grid extraction indices
        do_gl(:npc) = (/ .false., (.true., i=2,npc-1), .false. /)
        call get_gl_inds ( do_gl(:npc), f_inds, cg_inds, nglMax, ncg )

        ! Compute z_path & p_path
        z_path(1:no_ele) = (/z_glgrid(MaxVert:tan_inds(ptg_i):-1), &
                           & z_glgrid(tan_inds(ptg_i):MaxVert)/)
        del_zeta(2:mid) = 0.5_rp * ( z_path(c_inds(1:mid-1)) - &
          &                          z_path(c_inds(2:mid)) )
        del_zeta(mid+1:npc-1) = 0.5_rp * ( z_path(c_inds(mid+2:npc)) - &
          &                                z_path(c_inds(mid+1:npc-1)) )
        del_zeta((/1,npc/)) = 0.0_rp

        p_path(1:no_ele) = (/(p_glgrid(i),i=MaxVert,tan_inds(ptg_i),-1), &
                           & (p_glgrid(i),i=tan_inds(ptg_i),MaxVert)/)
        p_path_c(1:npc) = p_path(c_inds(1:npc))

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
              &  ETA_ZXP = eta_zxp_t(1:no_ele,:),                            &
              &  DO_CALC_T = do_calc_t(1:no_ele,:),                          &
              &  DO_CALC_HYD = do_calc_hyd(1:no_ele,:) )
            dh_dt_path_c(1:npc,:) = dh_dt_path(c_inds(1:npc),:)
            do_calc_hyd_c(1:npc,:) = do_calc_hyd(c_inds(1:npc),:)
            do_calc_t_c(1:npc,:) = do_calc_t(c_inds(1:npc),:)
            eta_zxp_t_c(1:npc,:) = eta_zxp_t(c_inds(1:npc),:)
            t_der_path_flags(1:no_ele) = any(do_calc_t(1:no_ele,:),2)
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
            dh_dt_path_c(1:npc,:) = dh_dt_path(c_inds(1:npc),:)
            do_calc_hyd_c(1:npc,:) = do_calc_hyd(c_inds(1:npc),:)
            do_calc_t_c(1:npc,:) = do_calc_t(c_inds(1:npc),:)
            eta_zxp_t_c(1:npc,:) = eta_zxp_t(c_inds(1:npc),:)
            t_der_path_flags(1:no_ele) = any(do_calc_t(1:no_ele,:),2)
          else
            call metrics ( tan_phi(ptg_i:ptg_i), tan_inds(ptg_i:ptg_i),      &
              &  Grids_tmp%phi_basis, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid, &
              &  orbIncline%values(1,maf)*Deg2Rad,                           &
              &  dummy, h_path(1:no_ele), phi_path(1:no_ele),                &
              &  t_path(1:no_ele), dhdz_path(1:no_ele), Req,                 &
              &  TAN_PHI_H_GRID = one_tan_ht, TAN_PHI_T_GRID = one_tan_temp )
          end if
        end if
!       dhdz_gw_path(f_inds) = dhdz_path(f_inds) * (/ ( gw, i = 1, nglMax/ng ) /)
        do i = 1, nglMax, ng ! Avoid a temp for (/ ( gw, i = 1, nglMax/ng ) /)
          dhdz_gw_path(f_inds(i:i+ng-1)) = dhdz_path(f_inds(i:i+ng-1)) * gw
        end do
        h_path(1:no_ele) = req + h_path(1:no_ele)
        h_path_c(1:npc) = h_path(c_inds(1:npc))
        t_path_m = t_path - del_temp ! for computing temperature derivatives
        t_path_p = t_path + del_temp ! for computing temperature derivatives
        t_path_c(1:npc) = t_path(c_inds(1:npc))
        ! Fill the diagnostic arrays for Nathaniel and Bill
!       tan_temps ( ptg_i ) = one_tan_temp ( 1 )
!       tan_hts ( ptg_i ) = one_tan_ht ( 1 )
!       reqs ( ptg_i ) = req
        !  ** Determine the eta_zxp_dw, eta_zxp_dn, eta_zxp_dv
        if ( spect_der ) then
          call eval_spect_path ( Grids_dw, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dw(1:no_ele,:), eta_zxp_dw(1:no_ele,:) )
          do_calc_dw_c(1:npc,:) = do_calc_dw(c_inds(1:npc),:)
          eta_zxp_dw_c(1:npc,:) = eta_zxp_dw(c_inds(1:npc),:)
          call eval_spect_path ( Grids_dn, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dn(1:no_ele,:), eta_zxp_dn(1:no_ele,:) )
          do_calc_dn_c(1:npc,:) = do_calc_dn(c_inds(1:npc),:)
          eta_zxp_dn_c(1:npc,:) = eta_zxp_dn(c_inds(1:npc),:)
          call eval_spect_path ( Grids_dv, firstSignal%lo, thisSideband, &
            & z_path(1:no_ele), phi_path(1:no_ele), &
            & do_calc_dv(1:no_ele,:), eta_zxp_dv(1:no_ele,:) )
          do_calc_dv_c(1:npc,:) = do_calc_dv(c_inds(1:npc),:)
          eta_zxp_dv_c(1:npc,:) = eta_zxp_dv(c_inds(1:npc),:)
        end if

        !??? This is never used for anything ???
!       tan_temp(ptg_i) = one_tan_temp(1)

        ! Now compute the eta_zp & do_calc_zp (for Zeta & Phi only)

        ! Things you do whether or not you are doing magnetic or cloud
        call comp_eta_docalc_no_frq ( Grids_f, z_path(1:no_ele), &
          &  phi_path(1:no_ele), eta_zp(1:no_ele,:), do_calc_zp(1:no_ele,:) )

       ! Now compute sps_path with a FAKE frequency, mainly to get the
       ! WATER (H2O) contribution for refraction calculations, but also
       ! to compute sps_path for all those with no frequency component

        Frq = 0.0
        call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
          & Frq, eta_zp(1:no_ele,:), &
          & do_calc_zp(1:no_ele,:), sps_path(1:no_ele,:),      &
          & do_calc_fzp(1:no_ele,:), eta_fzp(1:no_ele,:) )

        ! Special path quantities for cloud model
        if ( fwdModelConf%Incl_Cld ) then

          !set cloud parameters to zero
          iwc_path(1:no_ele,1) = 0.
          WC(1,1:no_ele)=iwc_path(1:no_ele,1)
          WC(2,1:no_ele)=0.
          IPSD(1:no_ele)=1000

          call comp_eta_docalc_no_frq ( Grids_Iwc, z_path(1:no_ele), &
            &  phi_path(1:no_ele), eta_iwc_zp(1:no_ele,:), &
            &  do_calc_iwc_zp(1:no_ele,:) )

          ! Compute IWC_PATH
          Frq = 0.0
          call comp_sps_path_frq ( Grids_iwc, firstSignal%lo, thisSideband, &
            & Frq, eta_zp(1:no_ele,:), &
            & do_calc_zp(1:no_ele,:), iwc_path(1:no_ele,:),      &
            & do_calc_iwc(1:no_ele,:), eta_iwc(1:no_ele,:) )
          WC(1,1:no_ele)=iwc_path(1:no_ele,1)
        end if

        ! Special path quantities for Polarized (magnetic) model
        if ( FwdModelConf%polarized ) then

          call comp_eta_docalc_no_frq ( Grids_Mag, z_path(1:no_ele), &
            &  phi_path(1:no_ele), eta_mag_zp(1:no_ele,:) )

          ! Compute the three components of MAG_PATH
          call comp_sps_path ( Grids_mag, 1, eta_mag_zp(1:no_ele,:), &
            & mag_path(1:no_ele,1:3) )

          rot = reshape(ECRtoFOV%values(9*mif-8:9*mif,maf), (/3,3/))

          do j = 1, no_ele
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

          ct => mag_path(1:no_ele,3)   ! cos(theta)
          stcp => mag_path(1:no_ele,1) ! sin(theta) cos(phi)
          stsp => mag_path(1:no_ele,2) ! sin(theta) sin(phi)
          h => mag_path(1:no_ele,4)    ! magnitude of magnetic field

          if ( index(switches,'mag') /= 0 ) then
            clean = index(switches,'clean') /= 0
            call dump ( h, 'H', clean )
            call dump ( ct, 'Cos(theta)', clean )
            call dump ( stcp, 'Sin(theta) Cos(phi)', clean )
            call dump ( stsp, 'Sin(theta) Sin(phi)', clean )
          end if

        end if

        if ( h2o_ind > 0 ) then
          call refractive_index ( p_path_c(1:npc), &
            &  t_path_c(1:npc), n_path(1:npc),     &
            &  h2o_path=sps_path(c_inds(1:npc), h2o_ind) )
        else
          call refractive_index ( p_path_c(1:npc), &
            &  t_path_c(1:npc), n_path(1:npc) )
        end if

        if ( temp_der ) then
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path(npc/2),&
             & one_tan_ht(1), tan_phi(ptg_i), Req, 0.0_rp,                &
             & ptg_angles(ptg_i), r, 1.0_rp, tan_dh_dt(1,:),              &
             & tan_d2h_dhdt(1,:), dx_dt(ptg_i,:), d2x_dxdt(ptg_i,:) )
        else
          call get_chi_angles ( 0.001*est_scGeocAlt(ptg_i), n_path(npc/2),&
             & one_tan_ht(1), tan_phi(ptg_i), Req, 0.0_rp,                &
             & ptg_angles(ptg_i), r, 1.0_rp )
        end if

        n_path(1:npc) = min ( n_path(1:npc), MaxRefraction ) + 1.0_rp

        call comp_refcor ( h_path_c(1:npc), n_path(1:npc), &
                      &  Req+one_tan_ht(1), del_s(1:npc), ref_corr(1:npc) )

        path_dsdh(f_inds(:nglMax)) = h_path(f_inds(:nglMax)) / &
          & ( sqrt( h_path(f_inds(:nglMax))**2 - (Req+one_tan_ht(1))**2 ) )
        dsdz_gw_path(f_inds(:nglMax)) = path_dsdh(f_inds(:nglMax)) * &
          & dhdz_gw_path(f_inds(:nglMax))

        ! Compute ALL the slabs_prep entities over the path's GL grid for this
        ! pointing & mmaf:

        call get_gl_slabs_arrays ( my_Catalog(thisSideband,:), p_path(1:no_ele), &
          &  t_path(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs, &
          &  fwdModelConf%Do_1D, true_path_flags(1:no_ele) )

        if ( temp_der ) then
          call get_gl_slabs_arrays ( my_Catalog(thisSideband,:), p_path(1:no_ele), &
            &  t_path_p(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs_p, &
            &  fwdModelConf%Do_1D, t_der_path_flags(1:no_ele) )
          call get_gl_slabs_arrays ( my_Catalog(thisSideband,:), p_path(1:no_ele), &
            &  t_path_m(1:no_ele), 0.001*losVel%values(1,maf), gl_slabs_m, &
            &  fwdModelConf%Do_1D, t_der_path_flags(1:no_ele) )
        end if

        ! Work out what frequencies we're using for --------------------------
        ! frequency averaging case for this pointing

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.

        if ( FwdModelConf%do_freq_avg ) then
          j = -1
          k = SIZE(PointingGrids(whichPointingGrid)% &
                       & oneGrid(grids(ptg_i))%frequencies)
          call Hunt ( min_ch_freq_grid, vel_cor &
                       & * PointingGrids(whichPointingGrid)% &
                       & oneGrid(grids(ptg_i))%frequencies, k, j, frq_i )
          call Hunt ( max_ch_freq_grid,vel_cor &
                       & * PointingGrids(whichPointingGrid)% &
                       & oneGrid(grids(ptg_i))%frequencies, k, frq_i, m )
          noFreqs = m - j + 1
          call allocate_test ( frequencies, noFreqs, "frequencies", moduleName )
          frequencies = PointingGrids(whichPointingGrid)% &
                      & oneGrid(grids(ptg_i))%frequencies(j:m)

          ! VELOCITY shift correction to frequency grid

          frequencies =  Vel_Cor * frequencies

        end if

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_End ( 'ForwardModel.MetricsEtc' )

        ! Loop over frequencies ----------------------------------------------

        if ( toggle(emit) .and. levels(emit) > 4 ) &
          & call Trace_Begin ( 'ForwardModel.FrequencyLoop' )

        do frq_i = 1, noFreqs
          if ( toggle(emit) .and. levels(emit) > 5 ) &
            & call Trace_Begin ('ForwardModel.Frequency ',index=frq_i)

          Frq = frequencies(frq_i)
          frqhk = 0.5_r8 * frq * h_over_k ! h nu / 2 k T

          ! Set up path quantities --------------------------------------

          tanh1_c(1:npc) = tanh( frqhk / t_path_c(1:npc) )

          ! Compute the sps_path for this Frequency
          call comp_sps_path_frq ( Grids_f, firstSignal%lo, thisSideband, &
            & Frq, eta_zp(1:no_ele,:),                            &
            & do_calc_zp(1:no_ele,:), sps_path(1:no_ele,:),       &
            & do_calc_fzp(1:no_ele,:), eta_fzp(1:no_ele,:) )

          call get_beta_path ( Frq,                               &
            &  p_path(1:no_ele), t_path_c(1:npc), tanh1_c(1:npc), &
            &  my_Catalog(thisSideband,:), beta_group,            &
            &  FwdModelConf%polarized,                            &
            &  gl_slabs, c_inds(1:npc), beta_path_c(1:npc,:),     &
            &  gl_slabs_m, t_path_m(1:no_ele),                    &
            &  gl_slabs_p, t_path_p(1:no_ele),                    &
            &  t_der_path_flags,                                  &
            &  dbeta_dt_path_c, dbeta_dw_path_c,                  &
            &  dbeta_dn_path_c, dbeta_dv_path_c )

          do_gl = .false.

          do_cld = .true. !JJ !for Jonathan use Only

          if ( FwdModelConf%incl_cld ) then
!          if ( do_cld ) then !JJ
            ! Compute Scattering source function based on temp_prof at all
            ! angles U for each temperature layer assuming a plane parallel
            ! atmosphere.

            call allocate_test ( Scat_ang, fwdModelConf%num_scattering_angles, 'Scat_ang', moduleName )

            call T_scat ( temp%values(:,inst), Frq, GPH%values(:,inst), &
            & 10.0**(-temp%template%surfs), vmrArray, nspec,            &
            & fwdModelConf%num_scattering_angles,                       &
            & fwdModelConf%num_azimuth_angles,                          &
            & fwdModelConf%num_ab_terms, fwdModelConf%num_size_bins,    &
            & fwdModelConf%no_cloud_species,                            &
            & scat_src%values, scat_alb%values, cld_ext%values, Scat_ang)

            call Deallocate_test ( vmrArray,'vmrArray',ModuleName )

            scat_src%template = temp%template

            call load_one_item_grid ( grids_tscat, scat_src, phitan, maf, &
              & fwdModelConf, .false., .true. )

            call allocate_test ( do_calc_tscat, no_ele, size(grids_tscat%values),           &
                               & 'do_calc_tscat', moduleName )

            call allocate_test ( do_calc_tscat_zp, no_ele, grids_tscat%p_len,               &
                               & 'do_calc_tscat_zp', moduleName )

            call allocate_test ( eta_tscat,     no_ele, size(grids_tscat%values),           &
                               & 'eta_tscat',     moduleName )
                        
            call allocate_test ( eta_tscat_zp,  no_ele, grids_tscat%p_len,                  &
                               & 'eta_tscat_zp',  moduleName )
            call allocate_test ( tscat_path,    no_ele, fwdModelConf%num_scattering_angles, &
                                 & 'tscat_path',  moduleName )

            call comp_eta_docalc_no_frq ( Grids_Tscat, z_path(1:no_ele), &
              &  phi_path(1:no_ele), eta_tscat_zp(1:no_ele,:), &
              &  do_calc_tscat_zp(1:no_ele,:) )

            Frq=0.0
            call comp_sps_path_frq ( Grids_tscat, firstSignal%lo, thisSideband, &
              & Frq, eta_tscat_zp(1:no_ele,:), &
              & do_calc_tscat_zp(1:no_ele,:), tscat_path(1:no_ele,:),      &
              & do_calc_tscat(1:no_ele,:), eta_tscat(1:no_ele,:) )

            call allocate_test ( tt_path, no_ele, 1, 'tt_path', moduleName )

            ! project Tscat onto LOS
            call interp_tscat ( tscat_path(1:no_ele,:), Scat_ang(:), phi_path(1:no_ele), tt_path )

            if ( .not. cld_fine ) then                 ! interpolate onto gl grids along the LOS

              scat_alb%template = temp%template
              cld_ext%template  = temp%template
                            
              call load_one_item_grid ( grids_salb,  scat_alb, phitan, maf, fwdModelConf, .false.)
              call load_one_item_grid ( grids_cext,  cld_ext,  phitan, maf, fwdModelConf, .false.)             

              do i=1, size(grids_salb%values)
                 if ( abs(grids_salb%values(i)) < TOL) then
                         grids_salb%values(i) = 0.0
                 end if
                 if ( abs(grids_cext%values(i)) < TOL) then
                         grids_cext%values(i) = 0.0
                 end if
              enddo

              call allocate_test (do_calc_salb, no_ele, size(grids_salb%values),'do_calc_salb',moduleName)

              call allocate_test (do_calc_salb_zp, no_ele, grids_salb%p_len,               &
                                 & 'do_calc_salb_zp', moduleName )

              call allocate_test (eta_salb,    no_ele, size(grids_salb%values), 'eta_salb', moduleName)

              call allocate_test (eta_salb_zp, no_ele, grids_salb%p_len, 'eta_salb_zp', moduleName)

              call allocate_test ( salb_path,  no_ele, 1, 'salb_path', moduleName )

              call allocate_test (do_calc_cext,no_ele,size(grids_cext%values),'do_calc_cext',moduleName)
              call allocate_test (do_calc_cext_zp, no_ele, grids_cext%p_len,               &
                                 & 'do_calc_cext_zp', moduleName )
              call allocate_test (eta_cext,    no_ele, size(grids_cext%values), 'eta_cext', moduleName)
              call allocate_test (eta_cext_zp, no_ele, grids_cext%p_len,  'eta_cext_zp', moduleName)
              call allocate_test (cext_path,   no_ele, 1, 'cext_path', moduleName)

              call comp_eta_docalc_no_frq ( Grids_salb, z_path(1:no_ele), &
                &  phi_path(1:no_ele), eta_salb_zp(1:no_ele,:), &
                &  do_calc_salb_zp(1:no_ele,:) )

              Frq=0.0
              call comp_sps_path_frq ( Grids_salb, firstSignal%lo, thisSideband, &
                & Frq, eta_salb_zp(1:no_ele,:), &
                & do_calc_salb_zp(1:no_ele,:), salb_path(1:no_ele,:),          &
                & do_calc_salb(1:no_ele,:), eta_salb(1:no_ele,:) )

              call comp_eta_docalc_no_frq ( Grids_cext, z_path(1:no_ele), &
                &  phi_path(1:no_ele), eta_cext_zp(1:no_ele,:), &
                &  do_calc_cext_zp(1:no_ele,:) )

              Frq=0.0
              call comp_sps_path_frq ( Grids_cext, firstSignal%lo, thisSideband, &
                & Frq, eta_cext_zp(1:no_ele,:), &
                & do_calc_cext_zp(1:no_ele,:), cext_path(1:no_ele,:),      &
                & do_calc_cext(1:no_ele,:), eta_cext(1:no_ele,:) )

              call convert_grid ( salb_path(1:no_ele,:), cext_path(1:no_ele,:),   &
                                & tt_path(1:no_ele,:), c_inds(1:npc),             &
                                & beta_path_cloud_c(1:npc), w0_path_c(1:npc),     &
                                & tt_path_c(1:npc) )

            else                           ! re-compute cext and w0 along the LOS

              call get_beta_path_cloud ( Frq,                               &
                &  p_path(1:no_ele), t_path(1:no_ele), tt_path(1:no_ele,:), &
                &  beta_group, c_inds(1:npc), beta_path_cloud_c(1:npc),     &
                &  w0_path_c(1:npc),  tt_path_c(1:npc),      &
                &  IPSD(1:no_ele),  WC(:,1:no_ele), fwdModelConf )

            end if

            do j = 1, npc
              alpha_path_c(j) = dot_product( sps_path(c_inds(j),:), &
                                    & beta_path_c(j,:) ) + &
                                    & beta_path_cloud_c(j)

              incoptdepth(j) = alpha_path_c(j) * del_s(j)
            end do

            ! Determine where to use Gauss-Legendre instead of a trapezoid.

            call path_contrib ( incoptdepth(1:npc), e_rflty, &
              & fwdModelConf%tolerance, do_gl(1:npc) )
          
          else ! Not cloud model

            do j = 1, npc ! Don't trust compilers to fuse loops
              alpha_path_c(j) = dot_product( sps_path(c_inds(j),:), &
                                    & beta_path_c(j,:) )

              incoptdepth(j) = alpha_path_c(j) * del_s(j)

              tt_path_c(j) =0.0
              w0_path_c(j) =0.0

            end do

            if ( .not. FwdModelConf%polarized ) then
              ! Determine where to use Gauss-Legendre for scalar instead of a trapezoid.

              call path_contrib ( incoptdepth(1:npc), e_rflty, &
                    & fwdModelConf%tolerance, do_gl(1:npc) )

            else ! extra stuff for polarized case

              call get_beta_path_polarized ( frq, h, my_Catalog(thisSideband,:), beta_group, &
                & gl_slabs, c_inds(1:npc), beta_path_polarized )

              ! We put an explicit extent of -1:1 for the first dimension in
              ! the hope a clever compiler will do better optimization with
              ! a constant extent.
              if ( any_der ) then
                ! Will need beta_path_polarized * tanh1_c
                ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
                ! to alpha here
                do j = 1, npc
                  beta_path_polarized(-1:1,j,:) = beta_path_polarized(-1:1,j,:) * tanh1_c(j)
                  alpha_path_polarized(-1:1,j) = matmul( beta_path_polarized(-1:1,j,:), &
                    & sps_path(c_inds(j),:) ) + 0.25 * alpha_path_c(j)
                  alpha_path_polarized(0,j) = alpha_path_polarized(0,j) + &
                    & 0.25 * alpha_path_c(j)
                end do
              else
                ! Won't need beta_path_polarized * tanh1_c
                ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
                ! to alpha here
                do j = 1, npc
                  alpha_path_polarized(-1:1,j) = matmul( beta_path_polarized(-1:1,j,:), &
                    & sps_path(c_inds(j),:) ) * tanh1_c(j) + &
                    & 0.25 * alpha_path_c(j)
                  alpha_path_polarized(0,j) = alpha_path_polarized(0,j) + &
                    & 0.25 * alpha_path_c(j)
                end do
              end if

              ! Turn sigma-, pi, sigma+ into 2X2 matrix incoptdepth_pol
              call opacity ( ct(1:npc), stcp(1:npc), stsp(1:npc), &
                & alpha_path_polarized(:,1:npc), incoptdepth_pol(:,:,1:npc) )

              ! We don't add unpolarized incremental optical depth to diagonal
              ! of polarized incremental optical depth because we added the
              ! scalar alpha_path to the sigma-, pi and sigma+ parts of
              ! alpha_path_polarized above.  If we did add it here, we would
              ! need 0.5 factors to scale unpolarized "power absorption" to
              ! get "field absorption"

              ! Do not trust the compiler to fuse loops
              do j = 1, npc
                incoptdepth_pol(1,1,j) = - incoptdepth_pol(1,1,j) * del_s(j)
                incoptdepth_pol(2,1,j) = - incoptdepth_pol(2,1,j) * del_s(j)
                incoptdepth_pol(1,2,j) = - incoptdepth_pol(1,2,j) * del_s(j)
                incoptdepth_pol(2,2,j) = - incoptdepth_pol(2,2,j) * del_s(j)

                ! deltau_pol = exp(incoptdepth_pol)
                call cs_expmat ( incoptdepth_pol(:,:,j), deltau_pol(:,:,j) )
              end do

              ! Determine where to do GL
              call path_contrib ( deltau_pol(:,:,1:npc), e_rflty, &
                 & fwdModelConf%tolerance, do_gl(1:npc) )

            end if

          end if   ! end of check cld 

          ! Where we don't do GL, replace the rectangle rule by the
          ! trapezoid rule.
          do j = 1, npc - 1
            if ( .not. do_gl(j) ) &
              & incoptdepth(j) = &
                & 0.5 * ( alpha_path_c(j) + alpha_path_c(j+1) ) * del_s(j)
          end do

          call get_GL_inds ( do_gl(1:npc), gl_inds, cg_inds, ngl, ncg )
          ! ngl is ng * count(do_gl)

          t_path_f(:ngl) = t_path(gl_inds(:ngl))
          tanh1_f(1:ngl) = tanh( frqhk / t_path_f(:ngl) )

          ! The derivatives that get_beta_path computes depend on which
          ! derivative arrays are allocated, not which ones are present.
          ! This avoids having four paths through the code, each with a
          ! different set of optional arguments.

          call get_beta_path ( Frq,                                   &
            & p_path(1:no_ele), t_path_f(:ngl), tanh1_f(1:ngl),       &
            & my_Catalog(thisSideband,:), beta_group,                 &
            & FwdModelConf%polarized,                                 &
            & gl_slabs, gl_inds(:ngl), beta_path_f(:ngl,:),           &
            & gl_slabs_m, t_path_m(1:no_ele),                         &
            & gl_slabs_p, t_path_p(1:no_ele),                         &
            & t_der_path_flags,                                       &
            & dbeta_dt_path_f, dbeta_dw_path_f,                       &
            & dbeta_dn_path_f, dbeta_dv_path_f )

          do j = 1, ngl ! loop around dot_product instead of doing sum(a*b,2)
                        ! to avoid path-length array temps
            alpha_path_f(j) = dot_product( sps_path(gl_inds(j),:),  &
                                  & beta_path_f(j,:) )
          end do

          ! Needed by both rad_tran and rad_tran_pol
          call two_d_t_script ( t_path_c(1:npc), tt_path_c(1:npc),  &  
            & w0_path_c, spaceRadiance%values(1,1), frq,            &
            & t_script(1:npc) )

          if ( .not. FwdModelConf%polarized ) then


          ! Compute SCALAR radiative transfer ---------------------------------------

          ! Compute radiative transfer ---------------------------------------

              call rad_tran ( gl_inds(1:ngl), cg_inds(1:ncg), e_rflty,     &
                & del_zeta(1:npc), alpha_path_c(1:npc), ref_corr(1:npc),   &
                & do_gl(1:npc), incoptdepth(1:npc), alpha_path_f(1:ngl),   &
                & dsdz_gw_path, t_script(1:npc),  &
                & tau(1:npc), RadV(frq_i), i_stop )

          else ! Polarized model

            i_stop = npc ! needed by drad_tran_df

            ! get the corrections to integrals for layers that need gl for
            ! the polarized species
            call get_beta_path_polarized ( frq, h, my_Catalog(thisSideband,:), &
              & beta_group, gl_slabs, gl_inds(:ngl), beta_path_polarized_f )

            ! The explicit -1:1 is written in the hope that a clever compiler
            ! can exploit it to optimize.
            if ( any_der ) then
              ! Will need beta_path_polarized_f * tanh1_f
              do j = 1, ngl
                beta_path_polarized_f(-1:1,j,:) = beta_path_polarized_f(-1:1,j,:) * tanh1_f(j)
                alpha_path_polarized_f(-1:1,j) = matmul( beta_path_polarized_f(-1:1,j,:), &
                  & sps_path(gl_inds(j),:) ) + 0.25 * alpha_path_f(j)
                alpha_path_polarized_f(0,j) = alpha_path_polarized_f(0,j) + &
                  & 0.25 * alpha_path_f(j)
              end do
            else
              ! Won't need beta_path_polarized_f * tanh1_f
              do j = 1, ngl
                alpha_path_polarized_f(-1:1,j) = matmul( beta_path_polarized_f(-1:1,j,:), &
                  & sps_path(gl_inds(j),:) ) * tanh1_f(j) + 0.25 * alpha_path_f(j)
                alpha_path_polarized_f(0,j) = alpha_path_polarized_f(0,j) + &
                  & 0.25 * alpha_path_f(j)
              end do
            end if

            call rad_tran_pol ( gl_inds(1:ngl), cg_inds(1:ncg), e_rflty,         &
              & del_zeta(1:npc), alpha_path_polarized(:,1:npc), ref_corr(1:npc), &
              & do_gl(1:npc), incoptdepth_pol(:,:,1:npc), deltau_pol(:,:,1:npc), &
              & alpha_path_polarized_f(:,1:ngl), dsdz_gw_path,                   &
              & ct, stcp, stsp, t_script(1:npc),  prod_pol(:,:,1:npc),           &
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

            call drad_tran_df ( c_inds(1:npc), gl_inds(1:ngl), del_zeta, Grids_f, &
              &  beta_path_c(1:npc,:), eta_fzp, sps_path, do_calc_fzp(1:no_ele,:), &
              &  beta_path_f, do_gl(1:npc), del_s(1:npc), ref_corr(1:npc), &
              &  dsdz_gw_path, t_script(1:npc), tau(1:npc), &
              &  i_stop, d_delta_df(1:npc,:), k_atmos_frq(frq_i,:) )

            if ( FwdModelConf%polarized ) then

              ! VMR derivatives for polarized radiance.
              ! Compute DE / Df from D Incoptdepth_pol / Df and put
              ! into DE_DF.
              call Get_D_Deltau_Pol_DF ( ct, stcp, stsp, c_inds(1:p_stop), &
                &  del_zeta, Grids_f, beta_path_polarized(:,1:p_stop,:),      &
                &  eta_fzp, do_calc_fzp(1:no_ele,:), sps_path, del_s(1:npc),  &
                &  incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),         &
                &  d_delta_df(1:npc,:), de_df(:,:,1:p_stop,:) )

              ! Compute D radiance / Df from Tau, Prod, T_Script
              ! and DE / Df.

              call mcrt_der ( t_script(1:npc), sqrt(e_rflty),            &
                & deltau_pol(:,:,1:npc), de_df(:,:,1:npc,:),             &
                & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop,       &
                & d_rad_pol_df )

              if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
                k_atmos_frq(frq_i,:) = real(d_rad_pol_df(1,1,:))
              else
                k_atmos_frq(frq_i,:) = real(d_rad_pol_df(2,2,:))
              end if

            end if

          end if

          if ( temp_der ) then

            ! get d Delta B / d T * d T / eta
            call dt_script_dt ( t_path_c(1:npc), eta_zxp_t_c(1:npc,:), frq, &
                              & d_t_scr_dt(1:npc,:) )

            dh_dt_path_f(:ngl,:) = dh_dt_path(gl_inds(:ngl),:)
            do_calc_t_f(:ngl,:) = do_calc_t(gl_inds(:ngl),:)
            eta_zxp_t_f(:ngl,:) = eta_zxp_t(gl_inds(:ngl),:)
            h_path_f(:ngl) = h_path(gl_inds(:ngl))
            sps_beta_dbeta_c(:npc) = SUM(sps_path(c_inds(1:npc),:) * &
              &                          beta_path_c(1:npc,:) * &
              &                          dbeta_dt_path_c(1:npc,:),DIM=2)
            sps_beta_dbeta_f(:ngl) = SUM(sps_path(gl_inds(1:ngl),:) * &
              &                          beta_path_f(1:ngl,:) * &
              &                          dbeta_dt_path_f(1:ngl,:),DIM=2)

            if ( .not. FwdModelConf%polarized ) then

              call drad_tran_dt ( del_zeta(1:npc ), h_path_c(1:npc),          &
                & t_path_c(1:npc), dh_dt_path_c(1:npc,:),                     &
                & alpha_path_c(1:npc), sps_beta_dbeta_c(:npc),                &
                & eta_zxp_t_c(1:npc,:), do_calc_t_c(1:npc,:),                 &
                & do_calc_hyd_c(1:npc,:), del_s(1:npc), ref_corr(1:npc),      &
                & Req + one_tan_ht(1), dh_dt_path(brkpt,:), do_gl(1:npc),     &
                & gl_inds(1:ngl), h_path_f(:ngl), t_path_f(:ngl),             &
                & dh_dt_path_f(:ngl,:), alpha_path_f(1:ngl),                  &
                & sps_beta_dbeta_f(:ngl), eta_zxp_t_f(:ngl,:),                &
                & do_calc_t_f(:ngl,:), path_dsdh, dhdz_gw_path, dsdz_gw_path, &
                & t_script(1:npc), d_t_scr_dt(1:npc,:), tau(1:npc), i_stop,   &
                & grids_tmp%deriv_flags, drad_dt )
              k_temp_frq(frq_i,:) = drad_dt

            else ! FwdModelConf%polarized

              ! Temperature derivatives for polarized radiance
              ! Compute DE / DT from D Incoptdepth_Pol / DT and put
              ! into DE_DT.

              call get_d_deltau_pol_dT ( frq, h, ct, stcp, stsp,              &
                & my_catalog(thisSideband,:), beta_group, gl_slabs_m, gl_slabs_p, &
                & t_path_c(1:p_stop), t_path_m(1:no_ele), t_path_p(1:no_ele), &
                & t_path_f(:ngl), beta_path_polarized(:,1:p_stop,:),          &
                & beta_path_polarized_f(:,1:ngl,:), sps_path,                 &
                & alpha_path_polarized(:,1:p_stop), &
                & alpha_path_polarized_f(:,1:ngl), &
                & sps_beta_dbeta_c(:npc), sps_beta_dbeta_f(:ngl), &
                & eta_zxp_t_c(1:p_stop,:), eta_zxp_t_f(:ngl,:), del_s(1:npc), &
                & c_inds(1:p_stop), gl_inds(:ngl), del_zeta(1:npc),           &
                & do_calc_t_c(1:p_stop,:), do_calc_t_f(:ngl,:), do_gl(1:p_stop), &
                & path_dsdh, dhdz_gw_path, dsdz_gw_path,                      &
                & incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),          &
                & h_path_c(1:npc), h_path_f(:ngl), dh_dt_path_c(1:p_stop,:),  &
                & dh_dt_path_f(:ngl,:), Req + one_tan_ht(1), dh_dt_path(brkpt,:), &
                & do_calc_hyd_c(1:p_stop,:), grids_tmp%deriv_flags,           &
                & de_dt(:,:,1:p_stop,:) )

              ! Compute D radiance / DT from Tau, Prod, T_Script, D_T_Scr_dT
              ! and DE / DT.

              call mcrt_der ( t_script(1:npc), sqrt(e_rflty),            &
                & deltau_pol(:,:,1:npc), de_dt(:,:,1:npc,:),             &
                & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop,       &
                & d_rad_pol_dt, d_t_scr_dt(1:npc,:) )

              if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
                k_temp_frq(frq_i,:) = real(d_rad_pol_dt(1,1,:))
              else
                k_temp_frq(frq_i,:) = real(d_rad_pol_dt(2,2,:))
              end if

            end if

          end if

          if ( spect_der ) then

            if ( .not. FwdModelConf%polarized ) then
              ! Spectroscopic derivative  wrt: W

              do_calc_dw_f(:ngl,:) = do_calc_dw(gl_inds(:ngl),:)
              eta_zxp_dw_f(:ngl,:) = eta_zxp_dw(gl_inds(:ngl),:)
              call drad_tran_dx ( del_zeta(1:npc), Grids_dw,             &
                &  dbeta_dw_path_c(1:npc,:), eta_zxp_dw_c(1:npc,:),      &
                &  sps_path(c_inds(1:npc),:), do_calc_dw_c(1:npc,:),     &
                &  dbeta_dw_path_f(:ngl,:), eta_zxp_dw_f(:ngl,:),        &
                &  sps_path(gl_inds(1:ngl),:), do_calc_dw_f(:ngl,:),     &
                &  do_gl(1:npc), gl_inds(:ngl), del_s(1:npc),            &
                &  ref_corr(1:npc), dhdz_gw_path,                        &
                &  t_script(1:npc), tau(1:npc), i_stop, drad_dw )

              k_spect_dw_frq(frq_i,1:1:f_len_dw) = drad_dw(1:1:f_len_dw)

              ! Spectroscopic derivative  wrt: N

              do_calc_dn_f(:ngl,:) = do_calc_dn(gl_inds(:ngl),:)
              eta_zxp_dn_f(:ngl,:) = eta_zxp_dn(gl_inds(:ngl),:)
              call drad_tran_dx ( del_zeta(1:npc), Grids_dn,             &
                &  dbeta_dn_path_c(1:npc,:), eta_zxp_dn_c(1:npc,:),      &
                &  sps_path(c_inds(1:npc),:), do_calc_dn_c(1:npc,:),     &
                &  dbeta_dn_path_f(:ngl,:), eta_zxp_dn_f(:ngl,:),        &
                &  sps_path(gl_inds(1:ngl),:), do_calc_dn_f(:ngl,:),     &
                &  do_gl(1:npc), gl_inds(:ngl), del_s(1:npc),            &
                &  ref_corr(1:npc), dhdz_gw_path,                        &
                &  t_script(1:npc), tau(1:npc), i_stop, drad_dn )

              k_spect_dn_frq(frq_i,1:f_len_dn) = drad_dn(1:f_len_dn)

              ! Spectroscopic derivative  wrt: Nu0

              do_calc_dv_f(:ngl,:) = do_calc_dv(gl_inds(:ngl),:)
              eta_zxp_dv_f(:ngl,:) = eta_zxp_dv(gl_inds(:ngl),:)
              call drad_tran_dx ( del_zeta(1:npc), Grids_dv,             &
                &  dbeta_dv_path_c(1:npc,:), eta_zxp_dv_c(1:npc,:),      &
                &  sps_path(c_inds(1:npc),:), do_calc_dv_c(1:npc,:),     &
                &  dbeta_dv_path_f(:ngl,:), eta_zxp_dv_f(:ngl,:),        &
                &  sps_path(gl_inds(1:ngl),:), do_calc_dv_f(:ngl,:),     &
                &  do_gl(1:npc), gl_inds(:ngl), del_s(1:npc),            &
                &  ref_corr(1:npc), dhdz_gw_path,                        &
                &  t_script(1:npc), tau(1:npc), i_stop, drad_dv )

              k_spect_dv_frq(frq_i,1:1:f_len_dv) = drad_dv(1:1:f_len_dv)
            end if

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
          ! Do DACs stuff for all DACs channels first
          do i = 1, noUsedDACS
            shapeInd = MatchSignal ( dacsFilterShapes%signal, &
              & fwdModelConf%signals(usedDacsSignals(i)), sideband = thisSideband )
            call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
              & RadV, DACsStaging(:,i) )
          end do
          ! Now go through channel by channel
          do i = 1, noUsedChannels
            sigInd = channels(i)%signal
            channel = channels(i)%used
            if ( channels(i)%dacs == 0 ) then
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), sideband = thisSideband, &
                & channel=channel )
              if ( shapeInd == 0 ) &
                & call MLSMessage ( MLSMSG_Error, ModuleName, &
                &    "No matching channel shape information" )
              call Freq_Avg ( frequencies, &
                &   FilterShapes(shapeInd)%FilterGrid,  &
                &   FilterShapes(shapeInd)%FilterShape, &
                &   RadV, Radiances(ptg_i,i) )
            else
              radiances(ptg_i,i) = DACsStaging(channel,channels(i)%dacs)
            end if
          end do
        else
          Radiances(ptg_i,1:noUsedChannels) = RadV(1:)
        end if

        ! Frequency averaging of derivatives if needed -----------------------

        ! Frequency Average the temperature derivatives with the appropriate
        ! filter shapes

        if ( temp_der ) then
          if ( fwdModelConf%do_freq_avg ) then
            ! Do DACs stuff for all DACs channels first
            do i = 1, noUsedDACS
              shapeInd = MatchSignal ( dacsFilterShapes%signal, &
                & fwdModelConf%signals(usedDacsSignals(i)), sideband = thisSideband )
              sv_i = 1
              do instance = 1, no_sv_p_t
                do surface = 1, n_t_zeta
                  call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
                    & k_temp_frq(:,sv_i), DACsStaging2(:,sv_i,i) )
                  sv_i = sv_i + 1
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! i -- DACs loop
            ! Now go through channel by channel
            do i = 1, noUsedChannels
              sigInd = channels(i)%signal
              channel = channels(i)%used
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), &
                & sideband = thisSideband, channel=channel )
              sv_i = 1
              do instance = 1, no_sv_p_t
                do surface = 1, n_t_zeta
                  if ( channels(i)%dacs == 0 ) then
                    call Freq_Avg ( frequencies, &
                      & FilterShapes(shapeInd)%FilterGrid, &
                      & FilterShapes(shapeInd)%FilterShape, &
                      & k_temp_frq(:,sv_i), r )
                  else
                    r = DACsStaging2 ( channel, sv_i, channels(i)%dacs )
                  end if
                  k_temp(i,ptg_i,surface,instance) = r
                  sv_i = sv_i + 1
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! Channel loop
          else
            k_temp(1:noUsedChannels,ptg_i,1:n_t_zeta,1:no_sv_p_t) = &
              & reshape( k_temp_frq(1:noUsedChannels,1:no_sv_p_t*n_t_zeta), &
                &        (/ noUsedChannels, n_t_zeta, no_sv_p_t /) )
          end if
        end if

        ! Frequency Average the atmospheric derivatives with the appropriate
        ! filter shapes

        if ( atmos_der ) then

          do k = 1, no_mol
            if ( fwdModelConf%moleculeDerivatives(mol_cat_index(k)) ) then
              if ( fwdModelConf%do_freq_avg ) then
                ! Do DACs stuff for all DACs channels first
                do i = 1, noUsedDACS
                  shapeInd = MatchSignal ( dacsFilterShapes%signal, &
                    & fwdModelConf%signals(usedDacsSignals(i)), sideband = thisSideband )
                  do sv_i = grids_f%l_v(k-1)+1, grids_f%l_v(k)
                    call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
                      & k_atmos_frq(:,sv_i), DACsStaging2(:,sv_i,i) )
                  end do                  ! Surface loop X Instance loop
                end do
                ! Now go through channel by channel
                do i = 1, noUsedChannels
                  sigInd = channels(i)%signal
                  channel = channels(i)%used
                  shapeInd = MatchSignal ( filterShapes%signal, &
                    & fwdModelConf%signals(sigInd),             &
                    & sideband = thisSideband, channel=channel )
                  do sv_i = grids_f%l_v(k-1)+1, grids_f%l_v(k)
                    if ( grids_f%deriv_flags(sv_i) ) then
                      if ( channels(i)%dacs == 0 ) then
                        call Freq_Avg ( frequencies,            &
                          & FilterShapes(shapeInd)%FilterGrid,  &
                          & FilterShapes(shapeInd)%FilterShape, &
                          & k_atmos_frq(1:noFreqs,sv_i), r )
                      else
                        r = DACsStaging2 ( channel, sv_i, channels(i)%dacs )
                      end if
                    else
                      r = 0.0
                    end if
                    k_atmos(i,ptg_i,sv_i) = r
                  end do                ! Surface loop X Instance loop
                end do                  ! Channel loop
              else                      ! Else not frequency averaging
                do sv_i = grids_f%l_v(k-1)+1, grids_f%l_v(k)
                  if ( grids_f%deriv_flags(sv_i) ) then
                    k_atmos(1:noUsedChannels,ptg_i,sv_i) = k_atmos_frq(1:noUsedChannels,sv_i)
                  else
                    k_atmos(1:noUsedChannels,ptg_i,sv_i) = 0.0
                  end if
                end do
              end if                    ! Frequency averaging or not
            end if                      ! Want derivatives for this specie?
          end do                        ! Loop over major molecules
          !
        end if                          ! Want derivatives for atmos

        ! Frequency Average the spectroscopic derivatives with the appropriate
        ! filter shapes

        if ( spect_der ) then

          do k = 1, no_mol
            if ( fwdModelConf%moleculeDerivatives(mol_cat_index(k)) ) then
              !  *** dI/dW
              call dI_dSomething ( grids_dw, k_spect_dw_frq, k_spect_dw, k )
              !  *** dI/dN
              call dI_dSomething ( grids_dn, k_spect_dn_frq, k_spect_dn, k )
              !  *** dI/dV
              call dI_dSomething ( grids_dv, k_spect_dv_frq, k_spect_dv, k )
            end if                        ! Want derivatives for this
          end do                          ! Loop over major molecules

          !??? So now we have k_spect_d?.  What do we do with them ???

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

      ! Now check that the angles are in the correct order.  If they
      ! are not it means (give or take some approximations in the
      ! horizontal according to Bill), that the rays crossed over
      ! between the tangent point and the spacecraft.  One could dream
      ! up all sorts of elegant schemes to get around that problem, but
      ! it's simplest just to bail out (and is certainly preferable to
      ! the infinite loop in the convolution (Hunt on angles) that
      ! results otherwise).

      do ptg_i = 2, no_tan_hts-1
        ! this is a temporary fix
        if ( ptg_angles(ptg_i) < ptg_angles(ptg_i-1) )  &
           & ptg_angles(ptg_i) = (ptg_angles(ptg_i-1) + ptg_angles(ptg_i+1))/2
!        if ( ptg_angles(ptg_i) < ptg_angles(ptg_i-1) )  &
!          & call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & 'Pointing angles in wrong order, too much refraction?' )
      end do

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

      do i = 1, noUsedChannels
        channel = channels(i)%used
        chanInd = channel + 1 - channels(i)%origin
        sigInd = channels(i)%signal
        ! Get the radiance
        thisRadiance =>  &
          GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=merge ( 0, fwdModelConf%signals(sigInd)%sideband, fwdModelConf%forceFoldedOutput ) )
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

        ! Here comes the Convolution codes
        update = ( thisSideband /= fwdModelConf%sidebandStart )

        if ( FwdModelConf%do_conv ) then

          whichPattern = -1
          minSuperset = huge(0)
          do j = 1, size(antennaPatterns)
            superset = AreSignalsSuperset ( antennaPatterns(j)%signals, &
              & fwdModelConf%signals( (/sigInd/) ), sideband=thisSideband, &
              & channel=channel )
            if ( superset >= 0 .and. superset <= minSuperset ) then
              minSuperset = superset
              whichPattern = j
            end if
          end do
          if ( whichPattern < 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching antenna patterns." )

          ! Now change channel from starting at 0 or 1 to definitely 1

          j = sv_t_len
          !??? Which grids_... should give the windowStart:windowFinish to use here ???
          if ( .not. temp_der .AND. .not. atmos_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, qtys, temp,ptan, &
               & thisRadiance, update, ptg_angles, Radiances(:,i), &
               & tan_chi_out-thisElev, &
               & dhdz_out, dx_dh_out, thisFraction,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & PTAN_DER=ptan_der )
          else if ( temp_der .AND. .not. atmos_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, qtys, temp, ptan,&
               & thisRadiance, update, ptg_angles, Radiances(:,i), &
               & tan_chi_out-thisElev, &
               & dhdz_out, dx_dh_out, thisFraction,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),        &
               & DX_DT=dx_dt, D2X_DXDT=d2x_dxdt, DXDT_TAN=dxdt_tan,            &
               & DXDT_SURFACE=dxdt_surface, PTAN_DER=ptan_der )
          else if ( atmos_der .AND. .not. temp_der ) then
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, qtys, temp, ptan,&
               & thisRadiance, update, ptg_angles, Radiances(:,i), &
               & tan_chi_out-thisElev, &
               & dhdz_out, dx_dh_out, thisFraction,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!              & DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,size(grids_f%values)/))) )
          else
            call convolve_all ( FwdModelConf, FwdModelIn, FwdModelExtra, maf,  &
               & chanInd, windowStart, windowFinish, qtys, temp, ptan,&
               & thisRadiance, update, ptg_angles, Radiances(:,i), &
               & tan_chi_out-thisElev, &
               & dhdz_out, dx_dh_out, thisFraction,                               &
               & antennaPatterns(whichPattern), Grids_tmp%deriv_flags,         &
               & Grids_f, Jacobian, fmStat%rows, SURF_ANGLE=surf_angle(1),     &
               & DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),        &
               & DX_DT=dx_dt, D2X_DXDT=d2x_dxdt, DXDT_TAN=dxdt_tan,            &
               & DXDT_SURFACE=dxdt_surface,                                    &
               & DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!              & DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,size(grids_f%values)/))) )
          end if

        else          ! No convolution needed ..

          !??? Which grids_... should give the windowStart:windowFinish to use here ???
          j = sv_t_len
          if ( .not. temp_der .AND. .not. atmos_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out-thisElev, &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisFraction, qtys, fmStat%rows, Jacobian,  &
              &  PTAN_DER=ptan_der)
          else if ( temp_der .AND. .not. atmos_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out-thisElev, &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisFraction, qtys, fmStat%rows, Jacobian,  &
              &  DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),   &
              &  PTAN_DER=ptan_der )
          else if ( atmos_der .AND. .not. temp_der ) then
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out-thisElev, &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisFraction, qtys, fmStat%rows, Jacobian,  &
              &  DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!             &  DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,size(grids_f%values)/))) )
          else
            call no_conv_at_all ( fwdModelConf, fwdModelIn, fwdModelExtra, maf, chanInd, &
              &  windowStart, windowFinish, temp, ptan, thisRadiance, update, &
              &  Grids_tmp%deriv_flags, ptg_angles, tan_chi_out-thisElev, &
              &  dhdz_out, dx_dh_out, Grids_f,                            &
              &  Radiances(:,i), thisFraction, qtys, fmStat%rows, Jacobian,  &
              &  DI_DT=DBLE(RESHAPE(k_temp(i,:,:,:),(/no_tan_hts,j/))),   &
              &  DI_DF=DBLE(k_atmos(i,:,:)), PTAN_DER=ptan_der )
!             &  DI_DF=DBLE(RESHAPE(k_atmos(i,:,:),(/no_tan_hts,size(grids_f%values)/))) )
          end if

        end if

      end do                            ! Channel loop

      if ( toggle(emit) .and. levels(emit) > 2 ) &
        & call trace_end ( 'ForwardModel.Convolution' )

      ! Deallocate maxNoPtgFreqs stuff
      call deallocate_test ( Radv, 'RadV', moduleName )

      call DestroyCompleteSlabs ( gl_slabs )
      if ( temp_der ) then
        call DestroyCompleteSlabs ( gl_slabs_p )
        call DestroyCompleteSlabs ( gl_slabs_m )
        call deallocate_test ( k_temp_frq,     'k_temp_frq',     moduleName )
      end if

      if ( atmos_der ) &
        & call deallocate_test ( k_atmos_frq,  'k_atmos_frq',    moduleName )

      if ( spect_der ) then
        call deallocate_test ( k_spect_dw_frq, 'k_spect_dw_frq', moduleName )
        call deallocate_test ( k_spect_dn_frq, 'k_spect_dn_frq', moduleName )
        call deallocate_test ( k_spect_dv_frq, 'k_spect_dv_frq', moduleName )
      end if

      if ( any_der ) &
        & call deallocate_test ( DACsStaging2, 'DACsStaging2',   moduleName )

      if ( toggle(emit) .and. levels(emit) > 1 ) &
        & call trace_end ( 'ForwardModel.Sideband ',index=thisSideband )

    end do            ! End of loop over sidebands -------------------------

    if ( toggle(emit) .and. levels(emit) > 0 ) &
      & call Trace_End ( 'ForwardModel.SidebandLoop' )

    !  **** DEBUG Printing cycle ...

! *** Create *seez* file for "nasty" purposes:

    if ( index(switches, 'seez') /= 0 ) then
      include 'dump_print_code.f9h'
    end if

! *** End of include

    call GetNameOfSignal ( firstSignal, sigName )
    i = index(sigName, '.B')
    j = index(sigName(i+2:), '.' )
    if ( j /= 0 ) sigName(i+j+1:) = ''

    if ( index(switches, 'rad') /= 0 ) then
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

      do i = 1, noUsedChannels
        channel = channels(i)%used - channels(i)%origin + 1
        print "(/, 'ch', i3.3, '_pfa_rad\ ', i3.3 )", channels(i)%used, k
        j = thisRadiance%template%noChans
        print "( 4(2x, 1pg15.8) )", &
          & thisRadiance%values(channel:channel+j*(k-1):j, maf)
      end do
      Print *

    end if

    !  **** End of Printing cycle ...

    ! Now deallocate lots of stuff

    call destroy_species_data ( my_catalog )
    call destroy_beta_group ( beta_group )

    deallocate ( mol_cat_index, stat=j )

    call Deallocate_Test ( WC, 'WC', moduleName )
    call Deallocate_Test ( Scat_ang, 'Scat_ang', moduleName )
    call Deallocate_Test ( IPSD, 'IPSD', moduleName )

    deallocate ( channels, stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_DeAllocate//'dincoptdepth_pol_dt' )

    call deallocate_test ( usedDACSSignals, 'usedDACSSignals', moduleName )

    ! DESTROY THE SPS DATA STUFF

    call deallocate_test ( z_glGrid, 'z_glGrid', moduleName )
    call deallocate_test ( p_glGrid, 'z_glGrid', moduleName )

    call deallocate_test ( h_glgrid,     'h_glgrid',     moduleName )
    call deallocate_test ( t_glgrid,     't_glgrid',     moduleName )
    call deallocate_test ( dhdz_glgrid,  'dhdz_glgrid',  moduleName )
    call deallocate_test ( dh_dt_glgrid, 'dh_dt_glgrid', moduleName )
    call deallocate_test ( ddhidhidtl0,  'ddhidhidtl0',  moduleName )

    call deallocate_test ( tan_inds,      'tan_inds',      moduleName )
    call deallocate_test ( tan_press,     'tan_press',     moduleName )
    call deallocate_test ( tan_phi,       'tan_phi',       moduleName )
    call deallocate_test ( est_scgeocalt, 'est_scgeocalt', moduleName )
    call deallocate_test ( tan_temp,      'tan_temp',      moduleName )

    ! Extra DEBUG for Nathaniel and Bill
!   call deallocate_test ( tan_hts,       'tan_hts',       moduleName )
!   call deallocate_test ( tan_temps,     'tan_temps',     moduleName )
!   call deallocate_test ( reqs,          'reqs',          moduleName )

    call destroygrids_t ( grids_f )
    call destroygrids_t ( grids_iwc )
    call destroygrids_t ( grids_mag )
    call destroygrids_t ( grids_tmp )
    call destroygrids_t ( grids_tscat )
    call destroygrids_t ( grids_salb )
    call destroygrids_t ( grids_cext )

    call deallocate_test ( dhdz_path,        'dhdz_path',        moduleName )
    call deallocate_test ( dhdz_gw_path,     'dhdz_gw_path',     moduleName )
    call deallocate_test ( dsdz_gw_path,     'dsdz_gw_path',     moduleName )
    call deallocate_test ( h_path,           'h_path',           moduleName )
    call deallocate_test ( h_path_c,         'h_path_c',         moduleName )
    call deallocate_test ( h_path_f,         'h_path_f',         moduleName )
    call deallocate_test ( p_path,           'p_path',           moduleName )
    call deallocate_test ( p_path_c,         'p_path_c',         moduleName )
    call deallocate_test ( path_dsdh,        'path_dsdh',        moduleName )
    call deallocate_test ( phi_path,         'phi_path',         moduleName )
    call deallocate_test ( t_path,           't_path',           moduleName )
    call deallocate_test ( t_path_c,         't_path_c',         moduleName )
    call deallocate_test ( z_path,           'z_path',           moduleName )

    call deallocate_test ( alpha_path_c,     'alpha_path_c',     moduleName )
    call deallocate_test ( alpha_path_f,     'alpha_path_f',     moduleName )
    call deallocate_test ( beta_path_c,      'beta_path_c',      moduleName )
    call deallocate_test ( beta_path_cloud_c,'beta_path_cloud_c',moduleName )
    call deallocate_test ( w0_path_c,        'w0_path_c',        moduleName )
    call deallocate_test ( tt_path_c,        'tt_path_c',        moduleName )
    call deallocate_test ( c_inds,           'c_inds',           moduleName )
    call deallocate_test ( cg_inds,          'cg_inds',          moduleName )
    call deallocate_test ( del_s,            'del_s',            moduleName )
    call deallocate_test ( do_gl,            'do_gl',            moduleName )
    call deallocate_test ( f_inds,           'f_inds',           moduleName )
    call deallocate_test ( gl_inds,          'gl_inds',          moduleName )
    call deallocate_test ( incoptdepth,      'incoptdept',       moduleName )
    call deallocate_test ( n_path,           'n_path',           moduleName )
    call deallocate_test ( ref_corr,         'ref_corr',         moduleName )
    call deallocate_test ( sps_beta_dbeta_c, 'sps_beta_dbeta_c', moduleName )
    call deallocate_test ( sps_beta_dbeta_f, 'sps_beta_dbeta_f', moduleName )
    call deallocate_test ( tau,              'tau',              moduleName )
    call deallocate_test ( tanh1_c,          'tanh1_c',          moduleName )
    call deallocate_test ( tanh1_f,          'tanh1_f',          moduleName )
    call deallocate_test ( t_script,         't_script',         moduleName )

    call deallocate_test ( beta_path_f,      'beta_path_f',      moduleName )
    call deallocate_test ( do_calc_fzp,      'do_calc_fzp',      moduleName )
    call deallocate_test ( do_calc_iwc,      'do_calc_iwc',      moduleName )
    call deallocate_test ( do_calc_iwc_zp,   'do_calc_iwc_zp',   moduleName )
    call deallocate_test ( do_calc_tscat,    'do_calc_tscat',    moduleName )
    call deallocate_test ( do_calc_tscat_zp, 'do_calc_tscat_zp', moduleName )
    call deallocate_test ( do_calc_salb,     'do_calc_salb',     moduleName )
    call deallocate_test ( do_calc_salb_zp,  'do_calc_salb_zp',  moduleName )
    call deallocate_test ( do_calc_cext,     'do_calc_salb',     moduleName )
    call deallocate_test ( do_calc_cext_zp,  'do_calc_salb_zp',  moduleName )
    call deallocate_test ( do_calc_zp,       'do_calc_zp',       moduleName )
    call deallocate_test ( eta_fzp,          'eta_fzp',          moduleName )
    call deallocate_test ( eta_iwc,          'eta_iwc',          moduleName )
    call deallocate_test ( eta_iwc_zp,       'eta_iwc_zp',       moduleName )
    call deallocate_test ( eta_tscat,        'eta_tscat',        moduleName )
    call deallocate_test ( eta_tscat_zp,     'eta_tscat_zp',     moduleName )
    call deallocate_test ( eta_salb,         'eta_salb',         moduleName )
    call deallocate_test ( eta_salb_zp,      'eta_salb_zp',      moduleName )
    call deallocate_test ( eta_cext,         'eta_cext',         moduleName )
    call deallocate_test ( eta_cext_zp,      'eta_cext_zp',      moduleName )
    call deallocate_test ( eta_mag_zp,       'eta_mag_zp',       moduleName )
    call deallocate_test ( eta_zp,           'eta_zp',           moduleName )
    call deallocate_test ( iwc_path,         'iwc_path',         moduleName )
    call deallocate_test ( tscat_path,       'tscat_path',       moduleName )
    call deallocate_test ( tt_path,          'tt_path',          moduleName )
    call deallocate_test ( salb_path,        'salb_path',        moduleName )
    call deallocate_test ( cext_path,        'cext_path',        moduleName )
    call deallocate_test ( mag_path,         'mag_path',         moduleName )
    call deallocate_test ( sps_path,         'sps_path',         moduleName )
    call deallocate_test ( true_path_flags,  'true_path_flags',  moduleName )

    call deallocate_test ( tan_chi_out,      'tan_chi_out',      moduleName )
    call deallocate_test ( dx_dh_out,        'dx_dh_out',        moduleName )
    call deallocate_test ( dhdz_out,         'dhdz_out',         moduleName )

    call deallocate_test ( DACsStaging,      'DACsStaging',      moduleName )

    if ( FwdModelConf%incl_cld ) &
      & call deallocate_test ( scat_src%values, 'scat_src%values', &
        &                      moduleName )

    if ( FwdModelConf%incl_cld ) &
      & call deallocate_test ( scat_alb%values, 'scat_alb%values', &
        &                      moduleName )

    if ( FwdModelConf%incl_cld ) &
      & call deallocate_test ( cld_ext%values, 'cld_ext%values', &
        &                      moduleName )

    if ( temp_der ) then
      deallocate ( k_temp, STAT=i )
      call deallocate_test ( d_t_scr_dt,      'd_t_scr_dt',      moduleName )
      call deallocate_test ( dRad_dt,         'dRad_dt',         moduleName )
      call deallocate_test ( dbeta_dt_path_c, 'dbeta_dt_path_c', moduleName )
      call deallocate_test ( dbeta_dt_path_f, 'dbeta_dt_path_f', moduleName )
      call deallocate_test ( dh_dt_path,      'dh_dt_path',      moduleName )
      call deallocate_test ( dh_dt_path_c,    'dh_dt_path_c',    moduleName )
      call deallocate_test ( dh_dt_path_f,    'dh_dt_path_f',    moduleName )
      call deallocate_test ( do_calc_hyd,     'do_calc_hyd',     moduleName )
      call deallocate_test ( do_calc_hyd_c,   'do_calc_hyd_c',   moduleName )
      call deallocate_test ( do_calc_t,       'do_calc_t',       moduleName )
      call deallocate_test ( do_calc_t_c,     'do_calc_t_c',     moduleName )
      call deallocate_test ( do_calc_t_f,     'do_calc_t_f',     moduleName )
      call deallocate_test ( eta_zxp_t,       'eta_zxp_t',       moduleName )
      call deallocate_test ( eta_zxp_t_c,     'eta_zxp_t_c',     moduleName )
      call deallocate_test ( eta_zxp_t_f,     'eta_zxp_t_f',     moduleName )
      call deallocate_test ( tan_dh_dt,       'tan_dh_dt',       moduleName )
      call deallocate_test ( tan_d2h_dhdt,    'tan_d2h_dhdt',    moduleName )
      call deallocate_test ( dxdt_tan,        'dxdt_tan',        moduleName )
      call deallocate_test ( d2xdxdt_tan,     'd2xdxdt_tan',     moduleName )
      call deallocate_test ( dxdt_surface,    'dxdt_surface',    moduleName )
      call deallocate_test ( t_der_path_flags,'t_der_path_flags',moduleName )
      call deallocate_test ( true_path_flags, 'true_path_flags', moduleName )
    end if

    if ( atmos_der ) then
      call deallocate_test ( d_delta_df,      'd_delta_df',      moduleName )
      call deallocate_test ( k_atmos,         'k_atmos',         moduleName )
    end if

    if ( spect_der ) then

      call deallocate_test ( dbeta_dw_path_c, 'dbeta_dw_path_c', moduleName )
      call deallocate_test ( dbeta_dw_path_f, 'dbeta_dw_path_f', moduleName )
      call deallocate_test ( dbeta_dn_path_c, 'dbeta_dn_path_c', moduleName )
      call deallocate_test ( dbeta_dn_path_f, 'dbeta_dn_path_f', moduleName )
      call deallocate_test ( dbeta_dv_path_c, 'dbeta_dv_path_c', moduleName )
      call deallocate_test ( dbeta_dv_path_f, 'dbeta_dv_path_f', moduleName )

      call deallocate_test ( do_calc_dw,      'do_calc_dw',      moduleName )
      call deallocate_test ( do_calc_dw_c,    'do_calc_dw_c',    moduleName )
      call deallocate_test ( do_calc_dw_f,    'do_calc_dw_f',    moduleName )
      call deallocate_test ( do_calc_dn,      'do_calc_dn',      moduleName )
      call deallocate_test ( do_calc_dn_c,    'do_calc_dn_c',    moduleName )
      call deallocate_test ( do_calc_dn_f,    'do_calc_dn_f',    moduleName )
      call deallocate_test ( do_calc_dv,      'do_calc_dv',      moduleName )
      call deallocate_test ( do_calc_dv_c,    'do_calc_dv_c',    moduleName )
      call deallocate_test ( do_calc_dv_f,    'do_calc_dv_f',    moduleName )

      call deallocate_test ( eta_zxp_dw,      'eta_zxp_dw',      moduleName )
      call deallocate_test ( eta_zxp_dw_c,    'eta_zxp_dw_c',    moduleName )
      call deallocate_test ( eta_zxp_dw_f,    'eta_zxp_dw_f',    moduleName )
      call deallocate_test ( eta_zxp_dn,      'eta_zxp_dn',      moduleName )
      call deallocate_test ( eta_zxp_dn_c,    'eta_zxp_dn_c',    moduleName )
      call deallocate_test ( eta_zxp_dn_f,    'eta_zxp_dn_f',    moduleName )
      call deallocate_test ( eta_zxp_dv,      'eta_zxp_dv',      moduleName )
      call deallocate_test ( eta_zxp_dv_c,    'eta_zxp_dv_c',    moduleName )
      call deallocate_test ( eta_zxp_dv_f,    'eta_zxp_dv_f',    moduleName )

      call deallocate_test ( drad_dw,         'drad_dw',         moduleName )
      call deallocate_test ( drad_dn,         'drad_dn',         moduleName )
      call deallocate_test ( drad_dv,         'drad_dv',         moduleName )

      deallocate ( k_spect_dw, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'k_spect_dw' )
      deallocate ( k_spect_dn, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'k_spect_dn' )
      deallocate ( k_spect_dv, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'k_spect_dv' )

    end if

    if ( FwdModelConf%polarized ) then
      deallocate ( alpha_path_polarized, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'alpha_path_polarized' )
      deallocate ( beta_path_polarized, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'beta_path_polarized' )
      deallocate ( alpha_path_polarized_f, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'alpha_path_polarized_f' )
      deallocate ( beta_path_polarized_f, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'beta_path_polarized_f' )
      deallocate ( gl_delta_polarized, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'gl_delta_polarized' )
      deallocate ( incoptdepth_pol, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'incoptdepth_pol' )
      deallocate ( incoptdepth_pol_gl, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'incoptdepth_pol_gl' )
      if ( atmos_der ) then
        deallocate ( d_rad_pol_df, stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate//'d_rad_pol_df' )
        deallocate ( de_df, stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate//'de_df' )
      end if
      if ( temp_der ) then
        deallocate ( d_rad_pol_dt, stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate//'d_rad_pol_dt' )
        deallocate ( de_dt, stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate//'de_dt' )
        deallocate ( dincoptdepth_pol_dt, stat=ier )
        if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate//'dincoptdepth_pol_dt' )
      end if
      deallocate ( deltau_pol, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'deltau_pol' )
      deallocate ( prod_pol, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'prod_pol' )
      deallocate ( tau_pol, stat=ier )
      if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//'tau_pol' )
    end if

    call deallocate_test ( ptg_angles,       'ptg_angles',     moduleName )
    call deallocate_test ( radiances,        'radiances',      moduleName )
    call deallocate_test ( dx_dt,            'dx_dt',          moduleName )
    call deallocate_test ( d2x_dxdt,         'd2x_dxdt',       moduleName )

    call deallocate_test ( sps_beta_dbeta_c, 'sps_beta_dbeta', moduleName )
    call deallocate_test ( sps_beta_dbeta_f, 'sps_beta_dbeta', moduleName )

    deallocate ( qtys, stat = ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_DeAllocate // 'QTYS' )

    if ( toggle(emit) ) then
      call trace_end ( 'ForwardModel MAF=', fmStat%maf )
    end if

  contains

    subroutine dI_dSomething ( grids_dS, k_spect_dS_frq, k_spect_dS, K )

    ! Do or don't do frequency averaging, to finish up derivative of I w.r.t. S

      type(grids_T), intent(in) :: grids_dS
      real(rp), dimension(:,:), intent(in) :: K_SPECT_DS_FRQ ! ****
      real(r4), dimension(:,:,:,:,:,:), intent(out) :: K_SPECT_DS
      integer, intent(in) :: K ! Which molecule

      integer :: I, Instance, JF, Surface, SV_I

      if ( fwdModelConf%do_freq_avg ) then
        ! Do DACs stuff for all DACs channels first
        do i = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%signal, &
            & fwdModelConf%signals(usedDacsSignals(i)), sideband = thisSideband )
          sv_i = Grids_dS%l_v(k-1)
          do instance = Grids_dS%WindowStart(k), Grids_dS%WindowFinish(k)
            do surface = 1, Grids_dS%l_z(k) - Grids_dS%l_z(k-1)
              do jf = 1, Grids_dS%l_f(k) - Grids_dS%l_f(k-1)
                sv_i = sv_i + 1
                call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
                  & k_spect_dS_frq(:,sv_i), DACsStaging2(:,sv_i,i) )
              end do   ! jf -- Frequencies loop
            end do     ! Surface loop
          end do       ! Instance loop
        end do         ! i -- DACS loop
        do i = 1, noUsedChannels
          sigInd = channels(i)%signal
          channel = channels(i)%used
          shapeInd = MatchSignal ( filterShapes%signal, &
            & fwdModelConf%signals(sigInd),             &
            & sideband = thisSideband, channel=channel )
          sv_i = Grids_dS%l_v(k-1)
          do instance = Grids_dS%WindowStart(k), Grids_dS%WindowFinish(k)
            do surface = 1, Grids_dS%l_z(k) - Grids_dS%l_z(k-1)
              do jf = 1, Grids_dS%l_f(k) - Grids_dS%l_f(k-1)
                sv_i = sv_i + 1
                if ( channels(i)%dacs == 0 ) then
                  call Freq_Avg ( frequencies,           &
                    & FilterShapes(shapeInd)%FilterGrid, &
                    & FilterShapes(shapeInd)%FilterShape,&
                    & k_spect_dS_frq(:,sv_i), r )
                else
                  r = DACsStaging2 ( channel, sv_i, channels(i)%dacs )
                end if
                k_spect_dS(i,ptg_i,jf,surface,instance,k) = r
              end do              ! Frequencies loop
            end do                ! Surface loop
          end do                  ! Instance loop
        end do                    ! Channel loop
      else                        ! else not frequency averaging
        k_spect_dS( 1:noUsedChannels, ptg_i, &
          &         1:Grids_dS%l_f(k) - Grids_dS%l_f(k-1), &
          &         1:Grids_dS%l_z(k) - Grids_dS%l_z(k-1), &
          &         Grids_dS%WindowStart(k):Grids_dS%WindowFinish(k), k ) = &
          & reshape(k_spect_dS_frq( 1:noUsedChannels, &
          &                         Grids_dS%l_v(k-1) + 1: Grids_dS%l_v(k) ), &
            &       (/ noUsedChannels, Grids_dS%l_f(k) - Grids_dS%l_f(k-1), &
            &          Grids_dS%l_z(k) - Grids_dS%l_z(k-1), &
            &          Grids_dS%WindowFinish(k) - Grids_dS%WindowStart(k) + 1 /) )
      end if                      ! Frequency averaging or not
    end subroutine dI_dSomething

    subroutine IdentifyDACS
      ! Compute NoUsedDACs
      ! Allocate and compute UsedDACSSignals and allocate DACsStaging.
      integer :: LBoundDACs, UBoundDACs    ! How many channels in a DAC
      logical :: signalFlag(size(fwdModelConf%signals))

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

    end subroutine IdentifyDACS

  end subroutine FullForwardModel

! =====     Private procedures     =====================================

  ! -------------------------------------------  Estimate_Tan_Phi  -----
  subroutine Estimate_Tan_Phi ( no_tan_hts, nlvl, maf, phitan, ptan, &
                              & scgeocalt, tan_press, &
                              & tan_phi, est_scgeocalt )

  ! Estimate Tan_Phi and SC_geoc_alt.

    use Allocate_Deallocate, only: Allocate_Test
    use MLSCommon, only: RP
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
    real(rp), dimension(:), intent(in) :: Tan_press

  ! Outputs -- would be intent(out) if we could say so.
  ! These are nullified here, so don't expect them to be deallocated if
  ! they were previously allocated.

    real(rp), dimension(:), pointer :: Tan_phi
    real(rp), dimension(:), pointer :: Est_scgeocalt

  ! Local variables
    integer :: I, J, JF, K
    real(rp) :: R, R1, R2       ! real variables for various uses
    real(rp), dimension(ptan%template%noSurfs) :: &
      & P_PATH, &               ! Pressure on path
      & T_PATH, &               ! Temperatures on path
      & Z_PATH                  ! Zeta on path

    nullify ( tan_phi, est_scgeocalt )
    call allocate_test ( tan_phi,       no_tan_hts, 'tan_phi',       moduleName )
    call allocate_test ( est_scgeocalt, no_tan_hts, 'est_scgeocalt', moduleName )
    j = no_tan_hts - nlvl

    tan_phi(1:j) = phitan%values(1,MAF)
    est_scgeocalt(1:j) = scGeocAlt%values(1,maf)

! Since the interpolateValues routine needs the OldX array to be sorted
! we have to sort ptan%values and re-arrange phitan%values & scgeocalt%values

    k = ptan%template%noSurfs

    z_path = ptan%values(1:k,maf)
    p_path = phitan%values(1:k,maf)
    t_path = scgeocalt%values(1:k,maf)

    ! Sort z_path.  Permute p_path and t_path the same way.
    do i = 2, k ! Invariant: z_path(1:i-1) are sorted.
      r = z_path(i)
      if ( r < z_path(i-1) ) then
        r1 = p_path(i)
        r2 = t_path(i)
        jf = i
        do ! Find where to insert R.  Make room as we go.
          z_path(jf) = z_path(jf-1)
          p_path(jf) = p_path(jf-1)
          t_path(jf) = t_path(jf-1)
          jf = jf - 1
          if ( jf == 1 ) exit
          if ( r >= z_path(jf-1) ) exit
        end do
        z_path(jf) = r
        p_path(jf) = r1
        t_path(jf) = r2
      end if
    end do

    call interpolateValues ( z_path, p_path, tan_press(j+1:no_tan_hts), &
      &  tan_phi(j+1:no_tan_hts), METHOD = 'L' )
    call interpolateValues ( z_path, t_path, tan_press(j+1:no_tan_hts), &
       & est_scgeocalt(j+1:no_tan_hts), METHOD='L' )

    tan_phi = tan_phi * deg2rad

  end subroutine Estimate_Tan_Phi

! ------------------------------------------------  not_used_here  -----
  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module FullForwardModel_m

! $Log$
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
