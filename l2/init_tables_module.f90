! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! init_tables_module includes two files: field_parm.f9h and field_add.f9h
! They will normally be created automatically by the auxiliary program
! init_gen based on the contents of the file field_names.txt
!              What This Means
! Don't edit field_parm.f9h and field_add.f9h directly
! nor declare fields and add_idents for fields to init_tables_module
! Instead add or delete any field names in field_names.txt

! It has happened that the allocate_test/deallocate_test caused
! a segment fault of the compiled code when the compiler
! is NAG--if that happens, try turning off optimization
! Another solution would be to insert some time-wasting computations
! ahead of deallocate_test, but that would be an awkward hack

! There is a parameter declared, ID_LAST_MAX, that limits
! how many args you can accumulate by acorn between calls to
! make_tree (see below for usage of acorn)
! If the program quits with the message "Accumulated too many ids in acorn;"
! hunt down its definition below and try doubling its current value

! Declaring the definitions is handled by the tree walker.

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Init_MLSSignals_m ! Everything. Init_MLSSignals, Field_First,
    ! Last_Signal_Field, Spec_First, Last_Signal_Spec, Numerous S_....
  use Init_Spectroscopy_m ! Everything.
  use INTRINSIC ! Everything. ADD_IDENT, BEGIN, D, F, FIRST_LIT,
    ! INIT_INTRINSIC, L, L_<several>, LAST_INTRINSIC_LIT,
    ! N, NADP, ND, NDP, NP, NR, P, S, T, <all>_INDICES,
    ! T_BOOLEAN, T_FIRST, T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE,
    ! T_STRING, S_TIME and Z are used here, but everything is included so
    ! that it can be gotten by USE INIT_TABLES_MODULE.
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MOLECULES ! Everything.
  use Units, only: Init_Units

  implicit NONE
  public ! This would be a MUCH LONGER list than the list of private
  !        names below.
  private :: ADD_IDENT, INIT_INTRINSIC, INIT_MOLECULES, INIT_SPECTROSCOPY

  ! The following will be used by subroutine acorn to build an array 
  ! id_cum(1:id_last) that can be passed directly to make_tree
  ! The advantage is that successive calls to acorn extend id_cum
  ! eliminating the (too) many continuation statements a single call to
  ! make_tree would require
  ! To exploit this, replace a many-line-spanning statement such as:
  !   call make_tree ( (/ &
  !      begin, stuff, .. &
  !      begin, morestuff, .. &
  !        ..
  !      begin, laststuff, .. &
  !      / ) )  
  ! with the following:
  !   id_last = 0
  !   call acorn ( (/begin, stuff, ../) )
  !   call acorn ( (/begin, morestuff, ../) )
  !        ..
  !   call acorn ( (/begin, laststuff, ../) )
  !   call make_tree ( id_cum(1:id_last) )
  !   ! remember to re-initialize id_last to zero before next call to acorn
  !
  ! In fact, you can split the original args in the call to make_tree
  ! into however many pieces you like. In other words, 
  ! you can do it line-by-line or once each occurrence of the token "begin"
  ! or according to some other sytem that makes sense to you.
  
  private :: ID_CUM, ID_LAST, ID_LAST_MAX, ACORN
  integer :: ID_LAST
  integer, parameter :: ID_LAST_MAX = 500     ! You will be told if too small
  integer, pointer, dimension(:) :: ID_CUM => null()

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Enumeration types:
  integer, parameter :: T_CHUNKDIVIDEMETHOD = Last_spectroscopy_type+1
  integer, parameter :: T_CRITICALMODULE = t_chunkDivideMethod+1
  integer, parameter :: T_FGRIDCOORD     = t_criticalmodule+1
  integer, parameter :: T_FILLMETHOD     = t_fGridCoord+1
  integer, parameter :: T_FWMTYPE        = t_fillmethod+1
  integer, parameter :: T_GRIDDEDORIGIN  = t_fwmType+1
  integer, parameter :: T_HGRIDTYPE      = t_griddedOrigin+1
  integer, parameter :: T_MASKS          = t_hgridtype+1
  integer, parameter :: T_MATRIX         = t_masks+1
  integer, parameter :: T_METHOD         = t_matrix+1
  integer, parameter :: T_MODULE         = t_method+1
  integer, parameter :: T_OUTPUTTYPE     = t_module+1
  integer, parameter :: T_QUANTITYTYPE   = t_outputtype+1
  integer, parameter :: T_SCALE          = t_quantitytype+1
  integer, parameter :: T_SPECIES        = t_scale+1
  integer, parameter :: T_UNITS          = t_species+1
  integer, parameter :: T_VGRIDCOORD     = t_units+1
  integer, parameter :: T_VGRIDTYPE      = t_vgridcoord+1
  integer, parameter :: T_LAST           = t_vgridtype
! Field indices:
! Don't edit the following file directly--it is generated automatically
! based on the file field_names.txt
  include 'field_parm.f9h'    
! Enumeration literals (there are more in INTRINSIC and MOLECULES):
! Don't edit the following file directly--it is generated automatically
! based on the file lit_names.txt
  include 'lit_parm.f9h'    
! Section identities.  Indices are in the order the sections are allowed to
! appear.  They're also used to index SECTION_ORDERING, so BE CAREFUL if
! you change them!
  integer, parameter :: Z_CHUNKDIVIDE    = 6
  integer, parameter :: Z_CONSTRUCT      = 7
  integer, parameter :: Z_FILL           = 8
  integer, parameter :: Z_GLOBALSETTINGS = 3
  integer, parameter :: Z_JOIN           = 10
  integer, parameter :: Z_MERGEGRIDS     = 5
  integer, parameter :: Z_MLSSIGNALS     = 1
  integer, parameter :: Z_OUTPUT         = 11
  integer, parameter :: Z_READAPRIORI    = 4
  integer, parameter :: Z_RETRIEVE       = 9
  integer, parameter :: Z_SPECTROSCOPY   = 2
  integer, parameter :: SECTION_FIRST = z_mlsSignals, &
                                SECTION_LAST = z_Output
! Specification indices don't overlap parameter indices, so a section can
! have both parameters and specifications:
  integer, parameter :: S_APRIORI            = last_Spectroscopy_Spec + 1
  integer, parameter :: S_BINSELECTOR        = s_apriori + 1
  integer, parameter :: S_CHUNKDIVIDE        = s_binselector + 1
  integer, parameter :: S_DESTROY            = s_chunkDivide + 1
  integer, parameter :: S_DUMP               = s_destroy + 1
  integer, parameter :: S_DUMPBLOCKS         = s_dump + 1
  integer, parameter :: S_EMPIRICALGEOMETRY  = s_dumpblocks + 1
  integer, parameter :: S_FGRID              = s_empiricalGeometry + 1
  integer, parameter :: S_FILL               = s_fGrid + 1
  integer, parameter :: S_FILLCOVARIANCE     = s_fill + 1
  integer, parameter :: S_FILLDIAGONAL       = s_fillcovariance + 1
  integer, parameter :: S_FORGE              = s_filldiagonal + 1
  integer, parameter :: S_FORWARDMODEL       = s_forge + 1
  integer, parameter :: S_FORWARDMODELGLOBAL = s_forwardModel + 1
  integer, parameter :: S_GRIDDED            = s_forwardModelGlobal + 1
  integer, parameter :: S_HGRID              = s_gridded + 1
  integer, parameter :: S_L1BRAD             = s_hgrid + 1
  integer, parameter :: S_L1BOA              = s_l1brad + 1
  integer, parameter :: S_L2AUX              = s_l1boa + 1
  integer, parameter :: S_L2GP               = s_l2aux + 1
  integer, parameter :: S_MATRIX             = s_l2gp + 1
  integer, parameter :: S_MERGE              = s_matrix + 1
  integer, parameter :: S_OUTPUT             = s_merge + 1
  integer, parameter :: S_QUANTITY           = s_output + 1
  integer, parameter :: S_RETRIEVE           = s_quantity + 1
  integer, parameter :: S_SIDS               = s_retrieve + 1
  integer, parameter :: S_SNOOP              = s_sids + 1
  integer, parameter :: S_SUBSET             = s_snoop + 1
  integer, parameter :: S_TEMPLATE           = s_subset + 1
  integer, parameter :: S_TRANSFER           = s_template + 1
  integer, parameter :: S_VECTOR             = s_transfer + 1
  integer, parameter :: S_VECTORTEMPLATE     = s_vector + 1
  integer, parameter :: S_VGRID              = s_vectortemplate + 1
  integer, parameter :: SPEC_LAST = s_vGrid

! Parameter names:
  ! In GlobalSettings section:
  integer, parameter :: FIRST_PARM = spec_last + 1
  integer, parameter :: P_ALLOW_CLIMATOLOGY_OVERLOADS = first_parm
  integer, parameter :: P_INPUT_VERSION_STRING        = p_allow_climatology_overloads + 1
  integer, parameter :: P_OUTPUT_VERSION_STRING       = p_input_version_string + 1
  integer, parameter :: P_VERSION_COMMENT             = p_output_version_string + 1
  integer, parameter :: P_CYCLE                       = p_version_comment + 1
  integer, parameter :: P_STARTTIME                   = p_cycle + 1
  integer, parameter :: P_ENDTIME                     = p_starttime + 1
  integer, parameter :: P_INSTRUMENT                  = p_endtime + 1
  ! In ChunkDivide section:
  integer, parameter :: P_CRITICAL_BANDS              = p_instrument + 1
  integer, parameter :: P_CRITICAL_SCANNING_MODULES   = p_critical_bands + 1
  integer, parameter :: P_HOME_GEOD_ANGLE             = p_critical_scanning_modules + 1
  integer, parameter :: P_HOME_MODULE                 = p_home_geod_angle + 1
  integer, parameter :: P_IDEAL_LENGTH                = p_home_module + 1
  integer, parameter :: P_IGNOREL1B                   = p_ideal_length + 1
  integer, parameter :: P_MAX_GAP                     = p_ignorel1b + 1
  integer, parameter :: P_NOCHUNKS                    = p_max_gap + 1
  integer, parameter :: P_OVERLAP                     = p_noChunks + 1
  integer, parameter :: P_SCAN_LOWER_LIMIT            = p_overlap + 1
  integer, parameter :: P_SCAN_UPPER_LIMIT            = p_scan_lower_limit + 1
  integer, parameter :: LAST_PARM = p_scan_upper_limit

! Table for section ordering:
  integer, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: |             globalSettings    chunkDivide       retrieve             |
!     | mlsSignals        readApriori       construct          join          |
!     |       spectroscopy      mergeGrids         fill             output   |
! ====|======================================================================|== From: ==
        (/OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! Start
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! mlsSignals
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! spectroscopy
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! globalSettings
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! readApriori
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,    0,  & ! mergeGrids
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,    0,  & ! chunkDivide
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,    0,  & ! Construct
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,    0,  & ! Fill
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,    0,  & ! Retrieve
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   OK,  & ! Join
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
    use TREE_TYPES, only: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
                          N_NAME_DEF, N_SECTION, N_SPEC_DEF

    call allocate_test(id_cum, id_last_max, &
      & 'id_cum', ModuleName)

  ! Put intrinsic predefined identifiers into the symbol table.
    call init_Spectroscopy ( t_last, field_last, last_lit, &
    & first_parm, last_parm, section_last, spec_last )

  ! Put nonintrinsic predefined identifiers into the symbol table.
    ! Put enumeration type names into the symbol table
    data_type_indices(t_chunkDivideMethod) = add_ident ( 'chunkDivideMethod' )
    data_type_indices(t_criticalmodule) =  add_ident ( 'criticalModule' )
    data_type_indices(t_fillmethod) =      add_ident ( 'fillMethod' )
    data_type_indices(t_fgridcoord) =      add_ident ( 'fGridCoord' )
    data_type_indices(t_fwmType) =         add_ident ( 'fwmType' )
    data_type_indices(t_griddedOrigin) =   add_ident ( 'griddedOrigin' )
    data_type_indices(t_hgridtype) =       add_ident ( 'hGridType' )
    data_type_indices(t_matrix) =          add_ident ( 'matrixType' )
    data_type_indices(t_method) =          add_ident ( 'method' )
    data_type_indices(t_module) =          add_ident ( 'module' )
    data_type_indices(t_outputtype) =      add_ident ( 'outputType' )
    data_type_indices(t_quantitytype) =    add_ident ( 'quantityType' )
    data_type_indices(t_scale) =           add_ident ( 'scale' )
    data_type_indices(t_species) =         add_ident ( 'species' )
    data_type_indices(t_units) =           add_ident ( 'units' )
    data_type_indices(t_vgridcoord) =      add_ident ( 'vGridCoord' )
    data_type_indices(t_vgridtype) =       add_ident ( 'vGridType' )
    ! Put field names into the symbol table.  Don't add ones that are
    ! put in by init_MLSSignals.
    ! Don't edit the following file directly--it is generated automatically
    ! based on the file field_names.txt
    include 'field_add.f9h'    
    ! Put enumeration literals into the symbol table.  Don't add ones
    ! that are already put in by init_intrinsic or init_molecules.
    ! Don't edit the following file directly--it is generated automatically
    ! based on the file lit_names.txt
    include 'lit_add.f9h'    
    ! Put parameter names into the symbol table:
    parm_indices(p_allow_climatology_overloads) = &
                                           add_ident ( 'AllowClimatologyOverloads' )
    parm_indices(p_input_version_string) = add_ident ( 'InputVersionString' )
    parm_indices(p_ignoreL1B) =            add_ident ( 'IgnoreL1B' )
    parm_indices(p_output_version_string) =add_ident ( 'OutputVersionString' )
    parm_indices(p_version_comment) =      add_ident ( 'VersionComment' )
    parm_indices(p_cycle) =                add_ident ( 'Cycle' )
    parm_indices(p_starttime) =            add_ident ( 'StartTime' )
    parm_indices(p_endtime) =              add_ident ( 'EndTime' )
    parm_indices(p_instrument) =           add_ident ( 'Instrument' )
    parm_indices(p_critical_bands) =       add_ident ( 'CriticalBands' )
    parm_indices(p_critical_scanning_modules) = &
                                           add_ident ( 'CriticalScanningModules' )
    parm_indices(p_home_geod_angle) =      add_ident ( 'HomeGeodAngle' )
    parm_indices(p_home_module) =          add_ident ( 'HomeModule' )
    parm_indices(p_ideal_length) =         add_ident ( 'IdealLength' )
    parm_indices(p_max_gap) =              add_ident ( 'MaxGap' )
    parm_indices(p_noChunks) =             add_ident ( 'noChunks' )
    parm_indices(p_overlap) =              add_ident ( 'Overlap' )
    parm_indices(p_scan_lower_limit) =     add_ident ( 'ScanLowerLimit' )
    parm_indices(p_scan_upper_limit) =     add_ident ( 'ScanUpperLimit' )
    ! Put section names into the symbol table:
    section_indices(z_chunkDivide) =       add_ident ( 'chunkDivide' )
    section_indices(z_construct) =         add_ident ( 'construct' )
    section_indices(z_fill) =              add_ident ( 'fill' )
    section_indices(z_globalSettings) =    add_ident ( 'globalSettings' )
    section_indices(z_join) =              add_ident ( 'join' )
    section_indices(z_mergeGrids) =        add_ident ( 'mergeGrids' )
    section_indices(z_mlsSignals) =        add_ident ( 'mlsSignals' )
    section_indices(z_output) =            add_ident ( 'output' )
    section_indices(z_readApriori) =       add_ident ( 'readApriori' )
    section_indices(z_retrieve) =          add_ident ( 'retrieve' )
    section_indices(z_spectroscopy) =      add_ident ( 'spectroscopy' )
    ! Put spec names into the symbol table.  Don't add ones that are
    ! put in by init_MLSSignals.
    spec_indices(s_apriori) =              add_ident ( 'apriori' )
    spec_indices(s_binSelector) =          add_ident ( 'binSelector' )
    spec_indices(s_chunkDivide) =          add_ident ( 'chunkDivide' )
    spec_indices(s_empiricalGeometry) =    add_ident ( 'EmpiricalGeometry' )
    spec_indices(s_destroy) =              add_ident ( 'destroy' )
    spec_indices(s_dump) =                 add_ident ( 'dump' )
    spec_indices(s_dumpblocks) =           add_ident ( 'dumpblocks' )
    spec_indices(s_fGrid) =                add_ident ( 'fGrid' )
    spec_indices(s_fill) =                 add_ident ( 'fill' )
    spec_indices(s_fillCovariance) =       add_ident ( 'fillCovariance' )
    spec_indices(s_fillDiagonal)   =       add_ident ( 'fillDiagonal' )
    spec_indices(s_forge) =                add_ident ( 'forge' )
    spec_indices(s_forwardModel) =         add_ident ( 'forwardModel' )
    spec_indices(s_forwardModelGlobal) =   add_ident ( 'forwardModelGlobal' )
    spec_indices(s_gridded) =              add_ident ( 'gridded' )
    spec_indices(s_hgrid) =                add_ident ( 'hgrid' )
    spec_indices(s_l1brad) =               add_ident ( 'l1brad' )
    spec_indices(s_l1boa) =                add_ident ( 'l1boa' )
    spec_indices(s_l2aux) =                add_ident ( 'l2aux' )
    spec_indices(s_l2gp) =                 add_ident ( 'l2gp' )
    spec_indices(s_matrix) =               add_ident ( 'matrix' )
    spec_indices(s_merge) =                add_ident ( 'merge' )
    spec_indices(s_output) =               add_ident ( 'output' )
    spec_indices(s_quantity) =             add_ident ( 'quantity' )
    spec_indices(s_retrieve) =             add_ident ( 'retrieve' )
    spec_indices(s_snoop) =                add_ident ( 'snoop' )
    spec_indices(s_subset) =               add_ident ( 'subset' )
    spec_indices(s_template) =             add_ident ( 'template' )
    spec_indices(s_transfer) =             add_ident ( 'transfer' )
    spec_indices(s_vector) =               add_ident ( 'vector' )
    spec_indices(s_vectortemplate) =       add_ident ( 'vectorTemplate' )
    spec_indices(s_vgrid) =                add_ident ( 'vgrid' )
    spec_indices(s_sids) =                 add_ident ( 'sids' )

  ! Now initialize the units tables.  Init_Units depends on the lit tables
  ! having been initialized.

    call init_units

  ! Definitions are represented by trees.  The notation in the comments
  ! for the trees is < root first_son ... last_son >.  This is sometimes
  ! called "Cambridge Polish Notation."  It was developed to represent
  ! LISP by McCarthy et. al. at MIT (in Cambridge, MA).

  ! Notice that in the argument for make_tree, the tree node id is at
  ! the END of the subtree, while in Cambridge Polish Notation it is at
  ! the BEGINNING of the subtree!

  ! Put the definition trees into the tree space before the parser runs.
  ! After the parsing is done, they're automatically "glued in" to the
  ! "left" of the trees that represent the input.  The tree-walker
  ! stumbles upon them in its normal course of operation, never really
  ! realizing they're special (because by then they're not).

  ! Start with the definitions of types. These are represented by trees of
  ! the form  < n_dt_def t_type_name l_lit ... l_lit >
    ! The intrinsic data types are defined in the intrinsic module
    ! Define the nonintrinsic enumerated types
    call make_tree ( (/ &
      begin, t+t_griddedOrigin, l+l_climatology, l+l_dao, l+l_ncep, &
             l+l_gloria, n+n_dt_def, &
      begin, t+t_chunkDivideMethod, l+l_fixed, l+l_even, l+l_orbital, n+n_dt_def, &
      begin, t+t_criticalModule, l+l_both, l+l_either, l+l_ghz, l+l_none, &
             l+l_thz, n+n_dt_def, &
      begin, t+t_fGridCoord, l+l_frequency, l+l_LSBFrequency, l+l_USBFrequency, &
             l+l_IntermediateFrequency, n+n_dt_def, &
      begin, t+t_fillMethod, l+l_gridded, l+l_estimatedNoise, l+l_explicit, &
             l+l_hydrostatic, l+l_addnoise, &
             l+l_isotope, l+l_l1b, l+l_l2aux, l+l_l2gp, l+l_vector, l+l_special, &
             l+l_rectanglefromlos,l+l_vGrid, n+n_dt_def, &
      begin, t+t_fwmType, l+l_linear, l+l_full, l+l_scan, l+l_cloudFull, n+n_dt_def, &
      begin, t+t_hGridType, l+l_explicit, l+l_fixed, l+l_fractional, &
             l+l_height, l+l_regular, l+l_l2gp, n+n_dt_def, &
      begin, t+t_masks, l+l_full_derivatives, l+l_linalg, n+n_dt_def, &
      begin, t+t_matrix, l+l_plain, l+l_cholesky, l+l_kronecker, l+l_spd, &
             n+n_dt_def, &
      begin, t+t_method, l+l_highcloud,l+l_lowcloud, l+l_newtonian, n+n_dt_def, &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def, &
      begin, t+t_outputType, l+l_l2aux, l+l_l2gp, l+l_l2dgg, l+l_l2pc, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_quantityType, l+l_baseline, l+l_boundarypressure, &
             l+l_chisqchan, l+l_chisqmmaf, l+l_chisqmmif, l+l_cloudIce, &
             l+l_cloudInducedRadiance, l+l_cloudExtinction, l+l_cloudRadSensitivity, &
             l+l_cloudWater, l+l_columnabundance, &
             l+l_dnwt_ajn, l+l_dnwt_axmax, l+l_dnwt_cait, &
             l+l_dnwt_diag, l+l_dnwt_dxdx, l+l_dnwt_dxdxl, &
             l+l_dnwt_dxn, l+l_dnwt_dxnl, l+l_dnwt_flag, l+l_dnwt_fnmin, &
             l+l_dnwt_fnorm, l+l_dnwt_gdx, l+l_dnwt_gfac, &
             l+l_dnwt_gradn, l+l_dnwt_sq, l+l_dnwt_sq, l+l_dnwt_sqt,&
             l+l_earthRefl, l+l_effectiveOpticalDepth, l+l_elevOffset, &
             l+l_extinction, l+l_gph, l+l_heightOffset, l+l_isotopeRatio, &
             l+l_jacobian_cols, l+l_jacobian_rows, &
             l+l_losTransFunc, l+l_losVel, &
             l+l_massMeanDiameterIce, l+l_massMeanDiameterWater, &
             l+l_numF, l+l_numJ, &
             l+l_orbitInclination, l+l_ptan, l+l_radiance, l+l_earthradius,&
             l+l_refGPH, l+l_sizedistribution, &
             l+l_scanResidual, l+l_scECI, l+l_scVel, l+l_scGeocAlt, &
             l+l_sidebandRatio, l+l_spaceRadiance, l+l_surfacetype, l+l_temperature,&
             l+l_tngtECI, l+l_tngtGeodAlt, l+l_tngtGeocAlt, &
             l+l_totalExtinction, l+l_vmr, n+n_dt_def, &
      begin, t+t_scale, l+l_apriori, & ! l+l_covariance, & !??? Later !???
             l+l_none, l+l_norm, n+n_dt_def, &
      begin, t+t_species, l+l_gph, l+l_gph_precision, l+l_temperature, &
             l+l_temperature_prec, n+n_dt_def, &
      begin, t+t_units, l+l_c, l+l_days, l+l_deg, l+l_degrees, &
             l+l_dimensionless, l+l_dimless, l+l_dl, l+l_ghz, &
             l+l_hours, l+l_hpa, l+l_hz, l+l_k, l+l_khz, l+l_km, l+l_logp, &
             l+l_m, l+l_maf, l+l_mafs, l+l_mb, l+l_meters, l+l_mhz, &
             l+l_mif, l+l_mifs, l+l_minutes, l+l_orbits, l+l_pa, l+l_ppbv, &
             l+l_ppmv, l+l_pptv, l+l_rad, l+l_radians, l+l_s, l+l_seconds, &
             l+l_thz, l+l_vmr, l+l_zeta, n+n_dt_def, &
      begin, t+t_vgridcoord, l+l_angle, l+l_geodAltitude, l+l_gph, l+l_none, &
             l+l_pressure, l+l_theta, l+l_zeta, n+n_dt_def, &
      begin, t+t_vgridtype, l+l_explicit, l+l_linear, l+l_logarithmic, &
             l+l_l2gp, n+n_dt_def /) )

    ! Define the relations between specs and fields, and the field types
    ! or names of other specifications allowed.  These are represented by
    ! trees of the form
    !  < n_spec_def s_spec_name
    !               < n_field_type f_field_name t_type ... t_type > ...
    !               < n_field_spec f_field_name s_spec ... s_spec > ...
    !               < n_dot f_field_name s_spec f_field_name ... >
    !  >
    ! The n_field_type, n_field_spec, and n_dot subtrees may appear in
    ! any quantity or order.
    ! The n_field_type subtree indicates the types allowed for a field.
    ! The n_field_spec subtree indicates the specifications whose names
    ! are allowed to appear for a field.
    ! The n_dot subtree indicates that the field given by the first
    ! f_field_name is required to be of the form spec_name.field_name,
    ! where spec_name is required to be a label of a specification of the
    ! type given by the s_spec son, and field_name is required to be
    ! present in the field given by the last f_field_name, which is
    ! required to be in a specification named by the next-to-last
    ! f_field_name ... of the specification named by the spec_name.

    call make_tree ( (/ &
      begin, s+s_time, np+n_spec_def, &
      begin, s+s_gridded, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_field, t+t_string, n+n_field_type, &
             begin, f+f_origin, t+t_griddedOrigin, n+n_field_type, &
             np+n_spec_def, &
!      begin, s+s_l2gp, np+n_spec_def, & ! To avoid forward reference in h/vGrid
      begin, s+s_hGrid, &
             begin, f+f_type, t+t_hGridType, nr+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_fraction, t+t_numeric, n+n_field_type, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_mif, t+t_numeric, n+n_field_type, &
             begin, f+f_interpolationfactor, t+t_numeric, n+n_field_type, &
             begin, f+f_inclination, t+t_numeric, n+n_field_type, &
             begin, f+f_spacing, t+t_numeric, n+n_field_type, &
             begin, f+f_origin, t+t_numeric, n+n_field_type, &
             begin, f+f_sourceL2GP, s+s_l2gp, n+n_field_spec, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_binSelector, &
             begin, f+f_signals, t+t_string, n+n_field_type, &
             begin, f+f_type, t+t_quantityType, nr+n_field_type, &
             begin, f+f_molecule, t+t_molecule, n+n_field_type, &
             begin, f+f_height, t+t_numeric_range, n+n_field_type, &
             begin, f+f_cost, t+t_numeric, nr+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_empiricalGeometry, &
             begin, f+f_terms, t+t_numeric, nr+n_field_type, &
             begin, f+f_iterations, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_merge, &  ! Must be AFTER S_Gridded
             begin, f+f_operational, s+s_gridded, n+n_field_spec, &
             begin, f+f_climatology, s+s_gridded, n+n_field_spec, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_scale, t+t_numeric, n+n_field_type, &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_chunkDivide, &
             begin, f+f_method, t+t_chunkDivideMethod, nr+n_field_type, &
             begin, f+f_noChunks, t+t_numeric, n+n_field_type, &
             begin, f+f_overlap, t+t_numeric, n+n_field_type, &
             begin, f+f_maxLength, t+t_numeric, n+n_field_type, &
             begin, f+f_noSlaves, t+t_numeric, n+n_field_type, &
             begin, f+f_homeModule, t+t_module, n+n_field_type, &
             begin, f+f_homeGeodAngle, t+t_numeric, n+n_field_type, &
             begin, f+f_scanLowerLimit, t+t_numeric_range, n+n_field_type, &
             begin, f+f_scanUpperLimit, t+t_numeric_range, n+n_field_type, &
             begin, f+f_criticalModules, t+t_criticalModule, n+n_field_type, &
             begin, f+f_maxGap, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_template, &
             begin, f+f_copy, n+n_field_type, &
             begin, f+f_apriori, n+n_field_type, &
             begin, f+f_autofill, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_fgrid, &
             begin, f+f_coordinate, t+t_fGridCoord, n+n_field_type, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             nadp+n_spec_def, &
      begin, s+s_vgrid, &
             begin, f+f_type, t+t_vGridType, nr+n_field_type, &
             begin, f+f_coordinate, t+t_vGridCoord, n+n_field_type, &
             begin, f+f_formula, t+t_numeric_range, n+n_field_type, &
             begin, f+f_number, t+t_numeric, n+n_field_type, &
             begin, f+f_start, t+t_numeric, n+n_field_type, &
             begin, f+f_stop, t+t_numeric, n+n_field_type, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             begin, f+f_sourceL2GP, s+s_l2gp, n+n_field_spec, &
             ndp+n_spec_def, &
      begin, s+s_forge, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_geodAngle, t+t_numeric, n+n_field_type, &
             begin, f+f_noMIFs, t+t_numeric, n+n_field_type, &
             begin, f+f_inclination, t+t_numeric, n+n_field_type, &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_quantity, & ! Must be AFTER s_hgrid and s_vgrid
             begin, f+f_hGrid, s+s_hgrid, n+n_field_spec, &
             begin, f+f_fGrid, s+s_fgrid, n+n_field_spec, &
             begin, f+f_sGrid, s+s_vgrid, n+n_field_spec, &
             begin, f+f_vGrid, s+s_vgrid, n+n_field_spec, &
             begin, f+f_logBasis, t+t_boolean, n+n_field_type, &
             begin, f+f_molecule, t+t_molecule, n+n_field_type, &
             begin, f+f_radiometer, s+s_radiometer, n+n_field_spec, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_signal, t+t_string, n+n_field_type, &
             begin, f+f_type, t+t_quantityType, n+n_field_type, &
             begin, f+f_unit, t+t_units, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_vectorTemplate, & ! Must be AFTER s_quantity
             begin, f+f_quantities, s+s_quantity, n+n_field_spec, &
             nadp+n_spec_def, &
      begin, s+s_vector, & ! Must be AFTER s_vectorTemplate
             begin, f+f_template, s+s_vectorTemplate, nr+n_field_spec, &
             begin, f+f_lengthScale, t+t_boolean, n+n_field_type, &
             begin, f+f_fraction, t+t_boolean, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_l2gp, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_precision, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_prefixSignal, t+t_boolean, n+n_field_type, &
             begin, f+f_swath, t+t_string, n+n_field_type, &
             begin, f+f_hdfVersion, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_l2aux, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_prefixSignal, t+t_boolean, n+n_field_type, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_sdname, t+t_string, n+n_field_type, &
             begin, f+f_hdfVersion, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_matrix, &  ! Must be AFTER s_vector
             begin, f+f_rows, s+s_vector, n+n_field_spec, &
             begin, f+f_columns, s+s_vector, nr+n_field_spec, &
             begin, f+f_type, t+t_matrix, n+n_field_type, &
             ndp+n_spec_def /) )
     call make_tree ( (/ &
      begin, s+s_dump, &
             begin, f+f_vector, s+s_vector, n+n_field_spec, &
             nadp+n_spec_def/) )

     id_last = 0
     call acorn((/begin, s+s_fill/))    ! Must be AFTER s_vector, s_matrix and s_climatology
     call acorn((/begin, f+f_boundaryPressure, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_dontMask, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_earthRadius, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_explicitValues, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_extinction, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_geocAltitudeQuantity, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_h2oQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_ignoreNegative, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_ignoreZero, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_integrationTime, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_interpolate, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_intrinsic, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_isPrecision, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_losQty, s+s_vector, f+f_template, f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_maxIterations, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_measurements, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_method, t+t_fillmethod, nr+n_field_type/))
     call acorn((/begin, f+f_model, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_multiplier, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_noFineGrid, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_noise, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_precision, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_ptanQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_quantity, s+s_vector, f+f_template, f+f_quantities, &
            nr+n_dot/))
     call acorn((/begin, f+f_ratioQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_radianceQuantity, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_refGPHQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_resetSeed, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_scVel, s+s_vector, f+f_template, f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_scECI, s+s_vector, f+f_template, f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_seed, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_sourceQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot/))
     call acorn((/begin, f+f_sourceL2GP, s+s_l2gp, n+n_field_spec/))
     call acorn((/begin, f+f_sourceL2AUX, s+s_l2aux, n+n_field_spec/))
     call acorn((/begin, f+f_sourceGrid, s+s_gridded, s+s_merge, n+n_field_spec/))
     call acorn((/begin, f+f_sourceSGrid, s+s_vGrid, n+n_field_spec/))
     call acorn((/begin, f+f_sourceVGrid, s+s_vGrid, n+n_field_spec/))
     call acorn((/begin, f+f_spread, t+t_boolean, n+n_field_type/))
     call acorn((/begin, f+f_systemTemperature, t+t_numeric, n+n_field_type/))
     call acorn((/begin, f+f_tngtECI, s+s_vector, f+f_template, f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_temperatureQuantity, s+s_vector, f+f_template, &
            f+f_quantities, n+n_dot/))
     call acorn((/begin, f+f_vmrQuantity, s+s_vector, f+f_template, f+f_quantities, &
            n+n_dot, ndp+n_spec_def /) )
     call make_tree ( id_cum(1:id_last) )

    call make_tree( (/ &
      begin, s+s_destroy, &
             begin, f+f_matrix, s+s_matrix, n+n_field_spec, &
             begin, f+f_vector, s+s_vector, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_transfer, &
             begin, f+f_source, s+s_vector, n+n_field_spec, &
             begin, f+f_destination, s+s_vector, n+n_field_spec, &
             nadp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_fillCovariance, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_matrix, s+s_matrix, nr+n_field_spec, &
             begin, f+f_diagonal, s+s_vector, nr+n_field_spec, &
             begin, f+f_lengthScale, s+s_vector, n+n_field_spec, &
             begin, f+f_fraction, s+s_vector, n+n_field_spec, &
             begin, f+f_invert, t+t_boolean, n+n_field_type, &
             begin, f+f_superDiagonal, s+s_vector, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_fillDiagonal, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_matrix, s+s_matrix, nr+n_field_spec, &
             begin, f+f_diagonal, s+s_vector, nr+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_output, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_type, t+t_outputType, nr+n_field_type, &
             begin, f+f_file, t+t_string, nr+n_field_type, &
             begin, f+f_quantities, s+s_l2aux, s+s_l2gp, s+s_matrix, nr+n_field_spec, &
             begin, f+f_overlaps, s+s_l2aux, s+s_l2gp, n+n_field_spec, &
             begin, f+f_packed, t+t_boolean, n+n_field_type, &
             begin, f+f_hdfVersion, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_subset, &  ! Must be AFTER s_vector
             begin, f+f_quantity, s+s_vector, f+f_template, f+f_quantities, &
                    nr+n_dot, &
             begin, f+f_ptanquantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_channels, t+t_numeric, t+t_numeric_range, n+n_field_type, &
             begin, f+f_height, t+t_numeric_range, n+n_field_type, &
             begin, f+f_ignore, t+t_boolean, n+n_field_type, &
             begin, f+f_mask, t+t_masks, n+n_field_type, &
             begin, f+f_opticalDepth, t+t_numeric, n+n_field_type, ndp+n_spec_def, &
      begin, s+s_forwardModel, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_atmos_der, t+t_boolean, n+n_field_type, &
             begin, f+f_default_spectroscopy, t+t_boolean, n+n_field_type, &
             begin, f+f_do_baseline, t+t_boolean, n+n_field_type, &
             begin, f+f_do_conv, t+t_boolean, n+n_field_type, &
             begin, f+f_do_freq_avg, t+t_boolean, n+n_field_type, &
             begin, f+f_integrationGrid, s+s_vGrid, n+n_field_spec, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_moleculeDerivatives, t+t_molecule, n+n_field_type, &
             begin, f+f_molecules, t+t_molecule, n+n_field_type, &
             begin, f+f_nabterms, t+t_numeric, n+n_field_type, &
             begin, f+f_nazimuthangles, t+t_numeric, n+n_field_type, &
             begin, f+f_ncloudspecies, t+t_numeric, n+n_field_type, &
             begin, f+f_nmodelsurfs, t+t_numeric, n+n_field_type, &
             begin, f+f_nscatteringangles, t+t_numeric, n+n_field_type, &
             begin, f+f_nsizebins, t+t_numeric, n+n_field_type, &
             begin, f+f_phiWindow, t+t_numeric, n+n_field_type, &
             begin, f+f_frqGap, t+t_numeric, n+n_field_type, &
             begin, f+f_signals, t+t_string, n+n_field_type, &
             begin, f+f_skipOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_cloud_der, t+t_numeric, n+n_field_type, &
             begin, f+f_cloud_width, t+t_numeric, n+n_field_type,&
             begin, f+f_cloud_fov, t+t_numeric, n+n_field_type, &
             begin, f+f_spect_der, t+t_boolean, n+n_field_type, &
             begin, f+f_tangentGrid, s+s_vGrid, n+n_field_spec, &
             begin, f+f_temp_der, t+t_boolean, n+n_field_type, &
             begin, f+f_tolerance, t+t_numeric, n+n_field_type, &
             begin, f+f_type, t+t_fwmType, nr+n_field_type, ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_forwardModelGlobal, &
             begin, f+f_antennaPatterns, t+t_string, n+n_field_type, &
             begin, f+f_l2pc, t+t_string, n+n_field_type, &
             begin, f+f_filterShapes, t+t_string, n+n_field_type, &
             begin, f+f_pointingGrids, t+t_string, n+n_field_type, np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_l1brad, &
             begin, f+f_file, t+t_string, n+n_field_type, np+n_spec_def, &
      begin, s+s_l1boa, &
             begin, f+f_file, t+t_string, n+n_field_type, np+n_spec_def &
             /) )
    call make_tree ( (/ &
      begin, s+s_retrieve, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_apriori, s+s_vector, n+n_field_spec, &
             begin, f+f_aprioriScale, t+t_numeric, n+n_field_type, &
             begin, f+f_columnScale, t+t_scale, n+n_field_type, &
             begin, f+f_covariance, s+s_matrix, n+n_field_spec, &
             begin, f+f_diagnostics, s+s_vector, n+n_field_spec, &
             begin, f+f_diagonal, t+t_boolean, n+n_field_type, &
             begin, f+f_forwardModel, s+s_forwardModel, nr+n_field_spec, &
             begin, f+f_fuzz, t+t_numeric, n+n_field_type, & ! Secret
             begin, f+f_fwdModelExtra, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelOut, s+s_vector, n+n_field_spec, &
             begin, f+f_jacobian, s+s_matrix, n+n_field_spec, &
             begin, f+f_lambda, t+t_numeric, n+n_field_type, &
             begin, f+f_maxF, t+t_numeric, n+n_field_type, &
             begin, f+f_maxJ, t+t_numeric, n+n_field_type, &
             begin, f+f_measurements, s+s_vector, nr+n_field_spec, &
             begin, f+f_measurementSD, s+s_vector, n+n_field_spec, &
             begin, f+f_method, t+t_method, n+n_field_type, &
             begin, f+f_outputCovariance, s+s_matrix, n+n_field_spec, &
             begin, f+f_outputSD, s+s_vector, n+n_field_spec, &
             begin, f+f_regOrders, t+t_numeric, n+n_field_type, &
             begin, f+f_regQuants, t+t_quantityType, n+n_field_type, &
             begin, f+f_regWeight, t+t_numeric, n+n_field_type, &
             begin, f+f_state, s+s_vector, nr+n_field_spec, &
             begin, f+f_toleranceA, t+t_numeric, n+n_field_type, &
             begin, f+f_toleranceF, t+t_numeric, n+n_field_type, &
             begin, f+f_toleranceR, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_sids, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_forwardModel, s+s_forwardModel, nr+n_field_spec, &
             begin, f+f_fwdModelExtra, s+s_vector, n+n_field_spec, &
             begin, f+f_fwdModelIn, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelOut, s+s_vector, nr+n_field_spec, &
             begin, f+f_destroyJacobian, t+t_boolean, n+n_field_type, &
             begin, f+f_perturbation, s+s_vector, n+n_field_spec, &
             begin, f+f_jacobian, s+s_matrix, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_snoop, &
             begin, f+f_comment, t+t_string, n+n_field_type, &
             begin, f+f_phaseName, t+t_string, n+n_field_type, &
             begin, f+f_level, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_dumpblocks, &
             begin, f+f_matrix, s+s_matrix, nr+n_field_spec, &
             begin, f+f_rowQuantity, s+s_quantity, nr+n_field_spec, &
             begin, f+f_colQuantity, s+s_quantity, nr+n_field_spec, &
             begin, f+f_rowSurfaces, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             begin, f+f_colSurfaces, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             begin, f+f_rowInstances, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             begin, f+f_colInstances, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             begin, f+f_rowChannels, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             begin, f+f_colChannels, t+t_numeric, t+t_numeric_range, &
                    n+n_field_type, &
             ndp+n_spec_def /) )
    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >
    call make_tree ( (/ &
      begin, z+z_mlsSignals, s+s_module, s+s_band, s+s_radiometer, &
                             s+s_signal, s+s_spectrometerType, s+s_time, n+n_section, &
      begin, z+z_spectroscopy, s+s_line, s+s_spectra, s+s_time, n+n_section, &
      begin, z+z_globalsettings, &
             begin, p+p_version_comment, t+t_string, n+n_name_def, &
             begin, p+p_input_version_string, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_allow_climatology_overloads, t+t_boolean, n+n_name_def,&
             begin, p+p_instrument, t+t_instrument, n+n_name_def,&
             begin, p+p_cycle, t+t_string, n+n_name_def, &
             begin, p+p_starttime, t+t_string, n+n_name_def, &
             begin, p+p_endtime, t+t_string, n+n_name_def, s+s_l1brad, s+s_l1boa, &
             s+s_empiricalGeometry, s+s_forwardModel, s+s_forwardModelGlobal, &
             s+s_time, s+s_vgrid, &
             s+s_fGrid, s+s_l1brad, s+s_l1boa, n+n_section, &
      begin, z+z_readapriori, s+s_time, s+s_gridded, s+s_l2gp, &
             s+s_l2aux, s+s_snoop, n+n_section, &
      begin, z+z_mergegrids, s+s_time, s+s_merge, n+n_section /) )
    call make_tree ( (/ &
      begin, z+z_chunkdivide, &
             begin, p+p_critical_bands, t+t_string, n+n_name_def, &
             begin, p+p_critical_scanning_modules, t+t_criticalModule, n+n_name_def, &
             begin, p+p_home_geod_angle, t+t_numeric, n+n_name_def, &
             begin, p+p_home_module, t+t_module, n+n_name_def, &
             begin, p+p_ideal_length, t+t_numeric, n+n_name_def, &
             begin, p+p_max_gap, t+t_numeric, n+n_name_def, &
             begin, p+p_noChunks, t+t_numeric, n+n_name_def, &
             begin, p+p_ignoreL1B, t+t_boolean, n+n_name_def, &
             begin, p+p_overlap, t+t_numeric, n+n_name_def, &
             begin, p+p_scan_lower_limit, t+t_numeric_range, n+n_name_def, &
             begin, p+p_scan_upper_limit, t+t_numeric_range, n+n_name_def, &
             s+s_time, s+s_chunkDivide, n+n_section, &
      begin, z+z_construct, s+s_hgrid, s+s_forge, s+s_quantity, &
             s+s_snoop, s+s_time, s+s_vectortemplate, n+n_section, &
      begin, z+z_fill, s+s_dump, s+s_fill, s+s_fillCovariance, s+s_fillDiagonal, &
                       s+s_matrix, s+s_destroy, s+s_snoop, s+s_time, s+s_vector, &
                       s+s_transfer, n+n_section, &
      begin, z+z_retrieve, s+s_dumpBlocks, s+s_matrix, s+s_retrieve, &
                           s+s_sids, s+s_snoop, s+s_subset, s+s_time, &
                           n+n_section, &
      begin, z+z_join, s+s_time, s+s_l2gp, s+s_l2aux, n+n_section, &
      begin, z+z_output, s+s_time, s+s_output, n+n_section /) )

    call deallocate_test(id_cum, &
      & 'id_cum', ModuleName)
  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine INIT_TABLES

  subroutine ACORN ( IDS )
    ! Build a tree specified by the "ids" array.
    ! Global use made of id_last, id_last_max and id_cum
    integer, intent(in) :: IDS(:)
    
    integer :: id_next
    
    id_next = id_last + size(ids)

    ! Note that we have to use print statements here, not MLSMessage, as
    ! we haven't yet decided how MLSMessage is going to work (Toolkit etc.)
    if ( size(ids) < 1 ) then
      print*,'Illegal num of args to acorn'
    else if ( id_next > id_last_max ) then
      print*,'Accumulated too many ids in acorn;' // &
        &' increase id_last_max in l2/init_tables_module.f90'
    else
      id_cum(id_last+1:id_next) = ids
    endif
    id_last = id_next
      
  end subroutine acorn

end module INIT_TABLES_MODULE

! $Log$
! Revision 2.198  2002/02/07 02:55:02  vsnyder
! Add a 'mask' field to the 'setup' spec
!
! Revision 2.197  2002/01/26 00:10:03  livesey
! Modified merge command
!
! Revision 2.196  2002/01/24 00:58:16  livesey
! Got the mergeGrids stuff set up properly
!
! Revision 2.195  2002/01/23 22:35:20  livesey
! Added Gloria format gridded data
!
! Revision 2.194  2002/01/23 21:51:11  pwagner
! hdfVersion new optional field in l2gp, l2aux, output
!
! Revision 2.193  2002/01/21 21:14:01  livesey
! Added binSelector stuff
!
! Revision 2.192  2002/01/18 00:24:58  livesey
! Added packed option for output (l2pcs)
!
! Revision 2.191  2002/01/08 18:15:34  livesey
! Made fwdModelExtra optional for sids
!
! Revision 2.190  2001/12/10 20:21:47  livesey
! Added code for regular hGrids
!
! Revision 2.189  2001/11/17 02:30:16  vsnyder
! Add L_numF and L_numJ
!
! Revision 2.188  2001/11/15 23:49:49  jonathan
! rename DF_spectroscopy to default_spectroscopy
!
! Revision 2.187  2001/11/15 23:33:53  jonathan
! add df_spectroscopy
!
! Revision 2.186  2001/11/14 01:50:10  livesey
! Replaced neither with none for criticalModules
!
! Revision 2.185  2001/11/09 00:04:39  livesey
! New ChunkDivide stuff
!
! Revision 2.184  2001/10/31 22:00:18  livesey
! Added phaseName to snooper
!
! Revision 2.183  2001/10/31 19:07:35  livesey
! Hooked fGrids into quantity templates
!
! Revision 2.182  2001/10/31 18:36:50  livesey
! Added fGrids stuff
!
! Revision 2.181  2001/10/30 01:47:22  livesey
! Minor change, got rid of some outdated column related stuff
!
! Revision 2.180  2001/10/24 22:35:47  dwu
! add FillDiagonal
!
! Revision 2.179  2001/10/23 16:37:30  pwagner
! isPrecision, precision new fields in Fill command
!
! Revision 2.178  2001/10/22 22:21:52  dwu
! add l_highcloud
!
! Revision 2.177  2001/10/19 22:33:21  pwagner
! Can destroy a vector or a matrix in Fill
!
! Revision 2.176  2001/10/19 17:11:36  livesey
! Whoops! left a syntax error in.
!
! Revision 2.175  2001/10/19 17:11:05  livesey
! Bug fix, and tidy up of acorn
!
! Revision 2.174  2001/10/18 23:59:55  pwagner
! Replaced remove with destroy; added multiplier to addnoise
!
! Revision 2.173  2001/10/18 23:42:31  livesey
! Added dump to fill
!
! Revision 2.172  2001/10/16 23:34:05  pwagner
! intrinsic, resetseed, seed fields added to addnoise method
!
! Revision 2.171  2001/10/15 22:10:54  livesey
! Added smoothing stuff to fillCovariance
!
! Revision 2.170  2001/10/05 20:17:36  vsnyder
! Disallow positional arguments on Snoop; add DNWT_FLAG quantity
!
! Revision 2.169  2001/10/05 01:19:03  vsnyder
! Added 'level' field to Snoop spec
!
! Revision 2.168  2001/10/03 23:33:48  vsnyder
! Create t_method type with l_newtonian and l_lowcloud lits
!
! Revision 2.167  2001/10/03 17:57:21  vsnyder
! Delete l_DegreesOfFreedom, add l-DNWT_...
!
! Revision 2.166  2001/10/02 23:40:37  vsnyder
! Add F_Diagnostics
!
! Revision 2.165  2001/10/02 20:34:34  livesey
! Added f_do_baseline
!
! Revision 2.164  2001/09/28 18:21:56  dwu
! add l_lowcloud
!
! Revision 2.163  2001/09/20 20:54:25  pwagner
! Replaced ignoreMask with dontMask
!
! Revision 2.162  2001/09/19 23:38:05  pwagner
! Permits Remove command in fill section
!
! Revision 2.161  2001/09/18 23:53:08  pwagner
! Replaced error field name with noise; began addNoise Fill method
!
! Revision 2.160  2001/09/17 23:13:09  livesey
! Added instrument stuff to global settings etc
!
! Revision 2.159  2001/09/14 23:33:43  pwagner
! Now should allow special fill of chi^2..
!
! Revision 2.158  2001/09/13 19:57:26  pwagner
! l_chisq... added
!
! Revision 2.157  2001/09/08 00:20:40  pwagner
! Works for new columnAbundance join
!
! Revision 2.156  2001/09/06 22:32:52  pwagner
! Undid column join; improved comments on acorn, field_names.txt, etc.
!
! Revision 2.155  2001/09/04 15:58:15  jonathan
! add cloud_fov, jonathan
!
! Revision 2.154  2001/08/23 16:25:01  pwagner
! Implemented init_gen build of init_tables_module.f90
!
! Revision 2.153  2001/08/08 23:49:25  pwagner
! Changed id_cum from allocatable to pointer; still dumps core if -gc option set
!
! Revision 2.152  2001/08/07 23:48:20  pwagner
! Added acorn routine to get around excessive continuations
!
! Revision 2.151  2001/07/31 23:25:32  pwagner
! Able to accept 2 new fields for join of column; does nothing yet
!
! Revision 2.150  2001/07/31 22:30:51  pwagner
! Don't try to split calls to make_tree within section
!
! Revision 2.149  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.148  2001/07/26 20:33:17  vsnyder
! Eliminate the 'extra' field of the 'matrix' spec
!
! Revision 2.147  2001/07/25 02:02:01  vsnyder
! Sort field names, remove 'lines' (it's defined in Init_Spectroscopy)
!
! Revision 2.146  2001/07/20 19:24:18  dwu
! add f_noFineGrid
!
! Revision 2.145  2001/07/20 17:06:45  dwu
! add f_extinction field for Fill cloud extinction calculation
!
! Revision 2.144  2001/07/19 18:05:57  dwu
! add sourceSGRID
!
! Revision 2.143  2001/07/19 17:42:48  dwu
! add f_sGrid field
!
! Revision 2.142  2001/07/19 00:09:31  dwu
! add l_rectanglefromlos
!
! Revision 2.141  2001/07/19 00:00:14  dwu
! make fewer continuous lines
!
! Revision 2.140  2001/07/18 23:42:15  dwu
! add f_losQty f_earthradius
!
! Revision 2.139  2001/07/18 23:10:52  dwu
! rename l_radiusofearth as l_earthradius
!
! Revision 2.138  2001/07/17 22:33:51  jonathan
! mixed a bug Too many continuation lines, paul
!
! Revision 2.137  2001/07/17 21:23:36  jonathan
! add cloud_width, jonathan
!
! Revision 2.136  2001/07/17 19:21:59  jonathan
! add surface as in instrinsic, jonathan
!
! Revision 2.135  2001/07/17 19:18:28  jonathan
! add sizedistribution as in instrinsic, jonathan
!
! Revision 2.134  2001/07/17 19:01:44  jonathan
! add radiusofearth as in instrinsic, jonathan
!
! Revision 2.133  2001/07/13 20:24:48  jonathan
! added cloudRadSensitivity as it is in intrinsic. -Jonathan
!
! Revision 2.132  2001/07/13 18:13:04  dwu
! add quantity losTransFunc
!
! Revision 2.131  2001/07/13 16:22:58  livesey
! Really a merge of 2.129 and 2.130.  The checkin of 2.130 was
! done incorrectly.
!
! Revision 2.130  2001/07/13 15:59:58  jonathan
! there are two l_cloudextinction, delet one, jonathan
!
! Revision 2.129  2001/07/12 23:27:56  livesey
! Got rid of s_cloudForwardModel
!
! Revision 2.128  2001/07/09 22:21:43  pwagner
! Added some fields for s_cloudforwardmodel
!
! Revision 2.127  2001/07/09 19:56:05  livesey
! Put back the required source field for l2aux, which is different from
! l2gp.  Then again, maybe it shouldn't be.  Oh well, that's an issue
! that can wait till HDF 5
!
! Revision 2.126  2001/07/09 18:18:45  pwagner
! Allows readapriori with format of l2aux similar to l2gp
!
! Revision 2.125  2001/07/07 03:55:19  livesey
! Removed l_cloudSensitivity as jonathan removed it from intrinsic.
!
! Revision 2.124  2001/06/26 20:11:08  livesey
! Couple of changes to subset
!
! Revision 2.123  2001/06/26 00:08:38  vsnyder
! Add RegQuants field to Retrieve
!
! Revision 2.122  2001/06/25 23:26:52  livesey
! Added new stuff for subset command
!
! Revision 2.121  2001/06/22 05:18:37  livesey
! Add `transfer' command in fill
!
! Revision 2.120  2001/06/22 01:26:25  vsnyder
! Add fields for regularization
!
! Revision 2.119  2001/06/21 20:06:42  vsnyder
! Make the tolerance field of the ForwardModel spec optional
!
! Revision 2.118  2001/06/21 15:05:30  livesey
! Added tolerance field to forwardModel
!
! Revision 2.117  2001/06/01 21:27:58  livesey
! Added outSD option to retrieve
!
! Revision 2.116  2001/05/31 22:14:20  livesey
! Bug fix.
!
! Revision 2.115  2001/05/31 22:07:40  livesey
! More cloud stuff.
!
! Revision 2.114  2001/05/31 20:29:55  livesey
! Added new vector types for cloud stuff.
!
! Revision 2.113  2001/05/30 20:16:27  vsnyder
! Add 'invert' field to 'fillCovariance' spec
!
! Revision 2.112  2001/05/29 23:22:48  livesey
! Some state vector types moved down to intrinsic to be with the others.
!
! Revision 2.111  2001/05/29 20:19:27  livesey
! Added cloudFull forward model type
!
! Revision 2.110  2001/05/25 20:26:31  livesey
! Added skipOverlaps option to ForwardModel
!
! Revision 2.109  2001/05/24 20:54:26  pwagner
! Deleted p_ccs..times
!
! Revision 2.108  2001/05/18 23:18:42  vsnyder
! Replace 'weight' field of 'retrieve' by 'measurementSD'
!
! Revision 2.107  2001/05/18 19:51:09  livesey
! Added interpolate option for l2gp fills
!
! Revision 2.106  2001/05/18 19:46:22  vsnyder
! Add secret 'fuzz' field to 'retrieve' command -- for testing
!
! Revision 2.105  2001/05/18 01:02:03  vsnyder
! Add Lambda, maxF, maxJ fields to Retrieve, deleted maxIterations
!
! Revision 2.104  2001/05/16 19:44:05  livesey
! Added estimatedNoise stuff
!
! Revision 2.103  2001/05/14 23:21:33  livesey
! Added frqGap
!
! Revision 2.102  2001/05/12 00:19:38  livesey
! Tidied up construct hGrid and vGrid from existing l2gp.
! However, parser cannot yet handle it as it involves forward
! references.
!
! Revision 2.101  2001/05/11 00:24:22  livesey
! Added l1brad and l1boa to global settings
!
! Revision 2.100  2001/05/10 23:26:05  livesey
! Added isotope scaling stuff to fill, and isotope ratio vector quantities.
!
! Revision 2.99  2001/05/10 16:31:14  livesey
! Added the prefixSignal option to joins
!
! Revision 2.98  2001/05/10 01:08:02  livesey
! Added destroyJacobian option to sids
!
! Revision 2.97  2001/05/08 21:53:05  livesey
! Added precision field to join.  Got rid of xStar, yStar and kStar
!
! Revision 2.96  2001/05/08 20:35:01  vsnyder
! Add fillCovariance spec
!
! Revision 2.95  2001/05/05 00:02:57  livesey
! Added stuff for numerical derivatives
!
! Revision 2.94  2001/05/04 19:54:08  pwagner
! Not sure why, but this one runs--unlike last one
!
! Revision 2.93  2001/05/04 18:31:41  pwagner
! Added stuff so global_settings replaces PCF functions
!
! Revision 2.92  2001/05/03 23:03:19  livesey
! Added stuff to support scan model.
!
! Revision 2.91  2001/05/02 20:29:14  livesey
! Removed f_frequency from forwardModel
!
! Revision 2.90  2001/05/02 03:13:36  livesey
! Changed dumpBlock to dumpBlocks, added instances arguments.
!
! Revision 2.89  2001/05/02 03:02:34  livesey
! Fixed options for dumpblock
!
! Revision 2.88  2001/05/02 02:35:42  livesey
! Added the DumpBlock items to retrieve, and associated f_ fields.
!
! Revision 2.87  2001/05/01 23:27:27  pwagner
! Added l_l2dgg literal type as possible output field type
!
! Revision 2.86  2001/04/27 21:53:43  livesey
! Removed the l2pc stuff
!
! Revision 2.85  2001/04/26 23:43:23  vsnyder
! Remove forwardModelIn from retrieve spec
!
! Revision 2.84  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.83  2001/04/26 00:07:51  livesey
! Stuff to support reading of l2pc files
!
! Revision 2.82  2001/04/25 20:34:19  livesey
! Now supports writing of l2pc files
!
! Revision 2.81  2001/04/24 20:05:07  livesey
! Added l2pc joining
!
! Revision 2.80  2001/04/23 23:45:01  vsnyder
! Add 'time' to MLSSignals, Spectroscopy, GlobalSettings and ChunkDivide
! sections.
!
! Revision 2.79  2001/04/23 23:24:10  livesey
! Added forge instruction
!
! Revision 2.78  2001/04/23 23:04:38  vsnyder
! Move 's_time' to 'intrinsic'
!
! Revision 2.77  2001/04/21 01:26:18  livesey
! Now supports creation of h/v grids from l2gp
!
! Revision 2.76  2001/04/20 23:12:14  livesey
! Added the `forge' stuff
!
! Revision 2.75  2001/04/20 17:13:00  livesey
! Added stuff for vGrid fill
!
! Revision 2.74  2001/04/19 20:04:44  livesey
! Added sideband ratio to quantityType
!
! Revision 2.73  2001/04/12 21:42:39  livesey
! Signal field in s_quantity now string.
!
! Revision 2.72  2001/04/12 19:48:45  livesey
! Removed channels from forwardModel, signal string now conveys that information.
!
! Revision 2.71  2001/04/10 23:19:03  livesey
! Changed signals in forward model config to string
!
! Revision 2.70  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.69  2001/04/10 00:01:53  vsnyder
! Prevent duplicate fields on 'matrix' spec
!
! Revision 2.68  2001/04/09 20:59:05  vsnyder
! Add C (for Celsius) unit and l_c name for it
!
! Revision 2.67  2001/04/06 21:53:20  vsnyder
! Specify 'no duplicate fields' for s_forwardModel
!
! Revision 2.66  2001/04/04 18:01:46  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 2.65  2001/04/04 02:13:23  vsnyder
! Added spectroscopy section
!
! Revision 2.64  2001/04/03 19:42:27  vsnyder
! Add call to Init_Spectroscopy
!
! Revision 2.63  2001/04/03 19:11:28  vsnyder
! Account for the changed order of initialization (intrinsic, Molecule,
! MLSSignals) and the new make_tree.f9h.  Change the call to initialize
! the next-lower-level one to account for requirements of new make_tree.
!
! Revision 2.62  2001/03/30 03:05:49  vsnyder
! Add 'antennaPatterns' field to 'forwardModelGlobal'
!
! Revision 2.61  2001/03/29 23:42:55  vsnyder
! Add 'filterShapes' field to forwardModelGlobal
!
! Revision 2.60  2001/03/29 22:07:25  livesey
! Added phiWindow
!
! Revision 2.59  2001/03/29 19:13:14  livesey
! Added stuff for gridded data fill.
!
! Revision 2.58  2001/03/28 23:43:49  livesey
! Added stuff for forward models
!
! Revision 2.57  2001/03/28 03:04:30  vsnyder
! remove f_per_decade, add f_formula for s_vGrid
!
! Revision 2.56  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.55  2001/03/17 21:06:57  livesey
! Added forward model type stuff
!
! Revision 2.54  2001/03/17 03:24:23  vsnyder
! Work on forwardModelGlobalSetup
!
! Revision 2.53  2001/03/17 02:24:45  livesey
! Added logBasis
!
! Revision 2.52  2001/03/17 00:50:38  livesey
! New forwardModel section (plus merge from Van)
!
! Revision 2.51  2001/03/16 21:49:10  vsnyder
! Add "pointingGrids" and "extraHeights" fields to ForwardModelGlobal
!
! Revision 2.50  2001/03/16 00:59:43  vsnyder
! Add support for defining types and literals in Init_MLSSignals_m
!
! Revision 2.49  2001/03/15 23:25:07  vsnyder
! Remove t_radiometer and its lits.
!
! Revision 2.48  2001/03/15 20:35:57  livesey
! Added new stuff for special fills
!
! Revision 2.47  2001/03/15 18:40:21  livesey
! Added some stuff for velocity, eci coordinates etc.
!
! Revision 2.46  2001/03/14 02:47:26  vsnyder
! Cosmetic improvements
!
! Revision 2.45  2001/03/14 02:04:53  vsnyder
! Moved MLSSignals_m to mlspgs/lib
!
! Revision 2.44  2001/03/08 21:40:03  livesey
! Added elevOffset
!
! Revision 2.43  2001/03/08 00:19:29  livesey
! Too many continuation lines before
!
! Revision 2.42  2001/03/07 23:59:52  vsnyder
! Add stuff for SIDS.
!
! Revision 2.41  2001/03/07 23:52:16  livesey
! Another bug fix
!
! Revision 2.40  2001/03/07 23:50:08  livesey
! Bug fix, whoops!
!
! Revision 2.39  2001/03/07 22:53:50  livesey
! Added snoop stuff
!
! Revision 2.38  2001/03/07 22:46:04  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.37  2001/03/07 22:41:55  livesey
! Reworked readapriori section
!
! Revision 2.36  2001/03/06 20:52:33  livesey
! Removed unpackOutput option for joins
!
! Revision 2.35  2001/03/05 01:19:03  livesey
! Minor changes
!
! Revision 2.34  2001/03/02 22:01:11  vsnyder
! Move USEs into some subroutines -- allows deleting some PRIVATE statements.
!
! Revision 2.33  2001/03/02 20:45:59  vsnyder
! Move make_tree to ~/mlspgs/srclib/make_tree.f9h.  This should improve
! maintainability of several versions of init_tables_module.
!
! Revision 2.32  2001/03/02 03:28:39  vsnyder
! Added more error checking (stack over/under-flow).  Added temporary l2load.
!
! Revision 2.31  2001/03/02 01:27:24  livesey
! Added more vector quantity types
!
! Revision 2.30  2001/02/28 21:36:14  livesey
! More changes in the signals section
!
! Revision 2.29  2001/02/28 17:13:10  livesey
! Removed temporary work-around.
!
! Revision 2.28  2001/02/28 01:15:36  livesey
! Reworked the signals stuff
!
! Revision 2.27  2001/02/27 19:00:19  livesey
! Added a few more items.
!
! Revision 2.26  2001/02/27 02:12:52  livesey
! Regular commit
!
! Revision 2.25  2001/02/22 23:26:30  vsnyder
! Removed "required" from "source" field of "l2gp" spec
!
! Revision 2.24  2001/02/20 23:17:33  livesey
! Added more stuff for fill.
!
! Revision 2.23  2001/02/20 18:44:02  livesey
! Removed firstIndexChannel
!
! Revision 2.22  2001/02/17 00:26:58  livesey
! Added more stuff to Fill section
!
! Revision 2.21  2001/02/16 19:19:04  vsnyder
! Added "diagonalOut" field to "retrieve" specification.
!
! Revision 2.20  2001/02/16 00:46:47  livesey
! Reworked Fill and parts of ReadApriori
!
! Revision 2.19  2001/02/13 19:50:07  vsnyder
! Specify individually required fields, not just "all fields required"
! Remove individual "public" specifiers, because there's a global "public"
!
! Revision 2.18  2001/02/09 19:29:03  vsnyder
! Turn on checking for duplicate and required fields for hgrid, l2aux,
! l2gp, output, subset, retrieve and vgrid specifications.
!
! Revision 2.17  2001/02/09 18:03:15  livesey
! Added f_instrumentmodule
!
! Revision 2.16  2001/02/08 21:54:25  livesey
! Remove L_None, now in intrinsic
!
! Revision 2.15  2001/02/08 21:13:32  vsnyder
! Move "theta" from init_tables_module to intrinsic.
!
! Revision 2.14  2001/02/08 01:52:43  vsnyder
! Provide for noDuplicates, allFields and noPositional checking.
! Turn on "no duplicates" and "noPositional" checking for several specs.
!
! Revision 2.13  2001/02/06 23:29:31  vsnyder
! Periodic commit
!
! Revision 2.12  2001/02/01 20:19:42  vsnyder
! Remove "gph" from the "molecule" type
!
! Revision 2.11  2001/02/01 01:23:36  vsnyder
! Account for the Molecules module
!
! Revision 2.10  2001/01/31 23:32:00  vsnyder
! Moved l_temperature l_temperature_prec l_ptan l_tangentheight l_sidebandratio
! l_scvel l_orbitinclination l_geodaltitude l_radiance l_scanresidual l_gph
! l_gph_precision l_refgph l_baseline l_extinction l_linewidth to
! intrinsic module
!
! Revision 2.9  2001/01/30 00:25:54  livesey
! Added L_REFGPH
!
! Revision 2.8  2001/01/26 19:02:24  vsnyder
! Changes for "retrieve" section.
!
! Revision 2.7  2001/01/18 02:00:35  vsnyder
! Define strings for more field names that were overlooked.  Type checking for
! "fill" is still wrong.
!
! Revision 2.6  2001/01/17 01:29:51  vsnyder
! Define s_subset's string table entry
!
! Revision 2.5  2001/01/10 21:02:44  vsnyder
! Add radiometer names, stuff for "retrieve" and "subset"
!
! Revision 2.4  2000/11/16 01:53:57  vsnyder
! Take timing out of sections that are only parameter settings.
!
! Revision 2.3  2000/11/16 01:46:16  vsnyder
! Revise section numbers so they don't overlap parameter numbers.
!
! Revision 2.2  2000/11/16 01:22:59  vsnyder
! Add a "time" spec to every section.
!
! Revision 2.1  2000/10/12 00:35:57  vsnyder
! Move intrinsic types and literals to "intrinsic" module
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
