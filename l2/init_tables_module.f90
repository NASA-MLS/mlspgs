! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  use INTRINSIC ! Everything. FIRST_LIT, FIRST_MOLECULE, INIT_INTRINSIC,
    ! L_<several>, LAST_INTRINSIC_LIT, LAST_MOLECULE, T_BOOLEAN, T_FIRST,
    ! T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE and T_STRING are used
    ! here, but everything is included so that it can be gotten by
    ! USE INIT_TABLES_MODULE.

  implicit NONE
  public ! This would be a MUCH LONGER list than the list of private
  !        names below.
  private :: ADD_IDENT, INIT_INTRINSIC, MAKE_TREE

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Enumeration types:
  integer, parameter :: T_CRITICALMODULE = t_last_intrinsic+1
  integer, parameter :: T_FILLMETHOD     = t_criticalmodule+1
  integer, parameter :: T_GRIDDEDORIGIN  = t_fillmethod+1
  integer, parameter :: T_HGRIDTYPE      = t_griddedOrigin+1
  integer, parameter :: T_MATRIX         = t_hgridtype+1
  integer, parameter :: T_MERGEMETHOD    = t_matrix+1
  integer, parameter :: T_METHOD         = t_mergemethod+1
  integer, parameter :: T_MERGESOURCE    = t_method+1
  integer, parameter :: T_MODULE         = t_mergesource+1
  integer, parameter :: T_MOLECULE       = t_module+1
  integer, parameter :: T_OUTPUTTYPE     = t_molecule+1
  integer, parameter :: T_QUANTITYTYPE   = t_outputtype+1
  integer, parameter :: T_RADIOMETER     = t_quantitytype+1
  integer, parameter :: T_SCALE          = t_radiometer+1
  integer, parameter :: T_SPECIES        = t_scale+1
  integer, parameter :: T_UNITS          = t_species+1
  integer, parameter :: T_VGRIDCOORD     = t_units+1
  integer, parameter :: T_VGRIDTYPE      = t_vgridcoord+1
  integer, parameter :: T_LAST           = t_vgridtype
  integer :: DATA_TYPE_INDICES(t_first:t_last)
! Field indices:
  integer, parameter :: F_APRIORI             = 1
  integer, parameter :: F_APRIORISCALE        = f_apriori + 1
  integer, parameter :: F_AUTOFILL            = f_aprioriScale + 1
  integer, parameter :: F_BAND                = f_autofill + 1
  integer, parameter :: F_CENTERFREQUENCY     = f_band + 1
  integer, parameter :: F_CHANNEL             = f_centerfrequency + 1
  integer, parameter :: F_CHANNELS            = f_channel + 1
  integer, parameter :: F_COLUMNS             = f_channels + 1
  integer, parameter :: F_COLUMNSCALE         = f_columns + 1
  integer, parameter :: F_COMMENT             = f_columnscale + 1
  integer, parameter :: F_COMPAREOVERLAPS     = f_comment + 1
  integer, parameter :: F_COORDINATE          = f_compareOverlaps + 1
  integer, parameter :: F_COPY                = f_coordinate + 1
  integer, parameter :: F_COVARIANCE          = f_copy + 1
  integer, parameter :: F_CRITERIA            = f_covariance + 1
  integer, parameter :: F_DEFERRED            = f_criteria + 1
  integer, parameter :: F_DIAGONAL            = f_deferred + 1
  integer, parameter :: F_DIAGONALOUT         = f_diagonal + 1
  integer, parameter :: F_EXPLICITVALUES      = f_diagonalout + 1
  integer, parameter :: F_FIELD               = f_explicitvalues + 1
  integer, parameter :: F_FILE                = f_field + 1
  integer, parameter :: F_FIRST               = f_file + 1
  integer, parameter :: F_FRACTION            = f_first + 1
  integer, parameter :: F_FREQUENCIES         = f_fraction + 1
  integer, parameter :: F_FREQUENCY           = f_frequencies + 1
  integer, parameter :: F_FWDMODELEXTRA       = f_frequency + 1
  integer, parameter :: F_FWDMODELIN          = f_fwdModelExtra + 1
  integer, parameter :: F_FWDMODELOUT         = f_fwdModelIn + 1
  integer, parameter :: F_GEOCALTITUDEQUANTITY= f_fwdModelOut + 1
  integer, parameter :: F_GPH                 = f_geocAltitudeQuantity + 1
  integer, parameter :: F_H2OQUANTITY         = f_gph + 1
  integer, parameter :: F_HEIGHT              = f_h2oquantity + 1
  integer, parameter :: F_HGRID               = f_height + 1
  integer, parameter :: F_INTERPOLATIONFACTOR = f_hgrid + 1
  integer, parameter :: F_JACOBIAN            = f_interpolationFactor + 1
  integer, parameter :: F_LAST                = f_jacobian + 1
  integer, parameter :: F_LENGTH              = f_last + 1
  integer, parameter :: F_LO                  = f_length + 1
  integer, parameter :: F_MAXITERATIONS       = f_lo + 1
  integer, parameter :: F_MATRIX              = f_maxIterations + 1
  integer, parameter :: F_MEASUREMENTS        = f_matrix + 1
  integer, parameter :: F_METHOD              = f_measurements + 1
  integer, parameter :: F_MIF                 = f_method + 1
  integer, parameter :: F_MODULE              = f_MIF + 1
  integer, parameter :: F_MOLECULE            = f_module + 1 
  integer, parameter :: F_NUMBER              = f_molecule + 1
  integer, parameter :: F_ORIGIN              = f_number + 1
  integer, parameter :: F_OUTPUTCOVARIANCE    = f_origin + 1
  integer, parameter :: F_OUTPUTOVERLAPS      = f_outputCovariance + 1
  integer, parameter :: F_OVERLAPS            = f_outputOverlaps + 1
  integer, parameter :: F_PER_DECADE          = f_overlaps + 1
  integer, parameter :: F_QUANTITIES          = f_per_decade + 1
  integer, parameter :: F_QUANTITY            = f_quantities + 1
  integer, parameter :: F_RADIOMETER          = f_quantity + 1
  integer, parameter :: F_RANGE               = f_radiometer + 1
  integer, parameter :: F_REFGPHQUANTITY      = f_range + 1
  integer, parameter :: F_ROWS                = f_refGPHQuantity + 1
  integer, parameter :: F_SCALE               = f_rows + 1
  integer, parameter :: F_SDNAME              = f_scale + 1
  integer, parameter :: F_SIGNAL              = f_sdname + 1
  integer, parameter :: F_SIGNALS             = f_signal + 1
  integer, parameter :: F_SOURCE              = f_signals + 1
  integer, parameter :: F_SOURCEAPRIORI       = f_source + 1
  integer, parameter :: F_SOURCEL2AUX         = f_sourceApriori + 1
  integer, parameter :: F_SOURCEL2GP          = f_sourcel2aux + 1
  integer, parameter :: F_SOURCEQUANTITY      = f_sourcel2gp + 1
  integer, parameter :: F_SPACECRAFT          = f_sourcequantity + 1
  integer, parameter :: F_SPECIES             = f_spacecraft + 1
  integer, parameter :: F_SPECTROMETER        = f_species + 1
  integer, parameter :: F_SPECTROMETERTYPE    = f_spectrometer + 1
  integer, parameter :: F_SPREAD              = f_spectrometerType + 1
  integer, parameter :: F_START               = f_spread + 1
  integer, parameter :: F_STATE               = f_start + 1
  integer, parameter :: F_STEP                = f_state + 1
  integer, parameter :: F_STOP                = f_step + 1
  integer, parameter :: F_SUFFIX              = f_stop + 1
  integer, parameter :: F_SWATH               = f_suffix + 1
  integer, parameter :: F_SWITCH              = f_swath + 1
  integer, parameter :: F_TEMPERATURE         = f_switch+1
  integer, parameter :: F_TOLERANCEA          = f_temperature + 1
  integer, parameter :: F_TOLERANCEF          = f_tolerancea + 1
  integer, parameter :: F_TOLERANCER          = f_tolerancef + 1
  integer, parameter :: F_VERSIONRANGE        = f_tolerancer + 1
  integer, parameter :: F_TEMPLATE            = f_versionRange + 1
  integer, parameter :: F_TEMPERATUREQUANTITY = f_template + 1
  integer, parameter :: F_TEST                = f_temperaturequantity + 1
  integer, parameter :: F_TYPE                = f_test + 1
  integer, parameter :: F_UNIT                = f_type + 1
  integer, parameter :: F_VALUES              = f_unit + 1
  integer, parameter :: F_VGRID               = f_values + 1
  integer, parameter :: F_WEIGHT              = f_vGrid + 1
  integer, parameter :: F_WIDTH               = f_weight + 1
  integer, parameter :: F_WIDTHS              = f_width + 1
  !??? Fields from here may be temporary for driving the forward model
  integer, parameter :: F_BILL                = f_widths + 1  !???
  integer, parameter :: F_ATMOS_DER           = f_bill + 1    !???
  integer, parameter :: F_DO_CONV             = f_atmos_der+1 !???
  integer, parameter :: F_DO_FREQ_AVG         = f_do_conv + 1 !???
  integer, parameter :: F_SPECT_DER           = f_do_freq_avg+1 !???
  integer, parameter :: F_TEMP_DER            = f_spect_der + 1 !???
  integer, parameter :: F_ZVI                 = f_temp_der+1  !???
  integer, parameter :: FIELD_FIRST = f_Apriori, FIELD_LAST = f_zvi
  integer :: FIELD_INDICES(field_first:field_last)
! Enumeration literals (there are more in INTRINSIC and MOLECULES):
  integer, parameter :: L_ANGLE         = last_intrinsic_lit + 1
  integer, parameter :: L_APRIORI       = l_angle + 1
  integer, parameter :: L_BOTH 	        = l_apriori + 1
  integer, parameter :: L_CHOLESKY      = l_both + 1
  integer, parameter :: L_CLIMATOLOGY   = l_cholesky+1
  integer, parameter :: L_COVARIANCE    = l_climatology + 1
  integer, parameter :: L_DAO 	        = l_covariance + 1
  integer, parameter :: L_DIRECT        = l_dao + 1
  integer, parameter :: L_EARTHREFL     = l_direct + 1
  integer, parameter :: L_EITHER        = l_earthRefl + 1
  integer, parameter :: L_ELEVOFFSET    = l_either + 1
  integer, parameter :: L_EXPLICIT      = l_elevOffset + 1
  integer, parameter :: L_FIXED         = l_explicit + 1
  integer, parameter :: L_FRACTIONAL    = l_fixed + 1
  integer, parameter :: L_HEIGHT        = l_fractional + 1
  integer, parameter :: L_HYDROSTATIC   = l_height + 1
  integer, parameter :: L_KRONECKER     = l_hydrostatic + 1
  integer, parameter :: L_L1B           = l_kronecker + 1
  integer, parameter :: L_L2AUX         = l_l1b + 1
  integer, parameter :: L_L2GP 	        = l_l2aux + 1
  integer, parameter :: L_LINEAR        = l_l2gp + 1
  integer, parameter :: L_LOGARITHMIC   = l_linear + 1
  integer, parameter :: L_NCEP 	        = l_logarithmic + 1
  integer, parameter :: L_NEITHER       = l_ncep + 1
  integer, parameter :: L_NEWTONIAN     = l_neither + 1
  integer, parameter :: L_NORM          = l_newtonian + 1
  integer, parameter :: L_ORBITINCLINE  = l_norm + 1
  integer, parameter :: L_PLAIN         = l_orbitIncline + 1
  integer, parameter :: L_PRESSURE      = l_plain + 1
  integer, parameter :: L_R1A           = l_pressure + 1
  integer, parameter :: L_R1B           = l_r1a + 1
  integer, parameter :: L_R2            = l_r1b + 1
  integer, parameter :: L_R3            = l_r2 + 1
  integer, parameter :: L_R4            = l_r3 + 1
  integer, parameter :: L_R5H           = l_r4 + 1
  integer, parameter :: L_R5V           = l_r5h + 1
  integer, parameter :: L_SCGEOCALT     = l_r5v + 1
  integer, parameter :: L_SPACERADIANCE = l_scGeocAlt + 1
  integer, parameter :: L_SPD           = l_spaceRadiance + 1
  integer, parameter :: L_VECTOR        = l_spd + 1
  integer, parameter :: L_WEIGHTED      = l_vector + 1
  integer, parameter :: LAST_LIT        = l_weighted
  integer :: LIT_INDICES(first_lit:last_lit)
! Parameter names:
  ! In GlobalSettings section:
  integer, parameter :: P_ALLOW_CLIMATOLOGY_OVERLOADS = 1
  integer, parameter :: P_INPUT_VERSION_STRING        = 2
  integer, parameter :: P_OUTPUT_VERSION_STRING       = 3
  integer, parameter :: P_VERSION_COMMENT             = 4
  ! In ChunkDivide section:
  integer, parameter :: P_CRITICAL_BANDS              = 5
  integer, parameter :: P_CRITICAL_SCANNING_MODULES   = 6
  integer, parameter :: P_HOME_GEOD_ANGLE             = 7
  integer, parameter :: P_HOME_MODULE                 = 8
  integer, parameter :: P_IDEAL_LENGTH                = 9
  integer, parameter :: P_MAX_GAP                     = 10
  integer, parameter :: P_OVERLAP                     = 11
  integer, parameter :: P_SCAN_LOWER_LIMIT            = 12
  integer, parameter :: P_SCAN_UPPER_LIMIT            = 13
  integer, parameter :: FIRST_PARM = P_ALLOW_CLIMATOLOGY_OVERLOADS
  integer, parameter :: LAST_PARM = P_SCAN_UPPER_LIMIT
  integer :: PARM_INDICES(first_parm:last_parm)
! Section identities.  Indices are in the order the sections are allowed to
! appear.  They're also used to index SECTION_ORDERING, so BE CAREFUL if
! you change them!
  integer, parameter :: Z_CHUNKDIVIDE    = 5
  integer, parameter :: Z_CONSTRUCT      = 6
  integer, parameter :: Z_FILL           = 7
  integer, parameter :: Z_GLOBALSETTINGS = 2
  integer, parameter :: Z_JOIN           = 9
  integer, parameter :: Z_MERGEAPRIORI   = 4
  integer, parameter :: Z_MLSSIGNALS     = 1
  integer, parameter :: Z_OUTPUT         = 10
  integer, parameter :: Z_READAPRIORI    = 3
  integer, parameter :: Z_RETRIEVE       = 8
  integer, parameter :: SECTION_FIRST = z_mlsSignals, &
                                SECTION_LAST = z_Output
  integer :: SECTION_INDICES(section_first:section_last)
! Specification indices don't overlap parameter indices, so a section can
! have both parameters and specifications:
  integer, parameter :: S_APRIORI            = last_parm + 1
  integer, parameter :: S_BAND               = s_apriori + 1
  integer, parameter :: S_CREATE             = s_band + 1
  integer, parameter :: S_FILL               = s_create + 1
  integer, parameter :: S_FORWARDMODEL       = s_fill + 1
  integer, parameter :: S_FORWARDMODELGLOBAL = s_forwardModel + 1 !???
  integer, parameter :: S_GRIDDED            = s_forwardModelGlobal + 1
  integer, parameter :: S_HGRID              = s_gridded + 1
  integer, parameter :: S_L2GP               = s_hgrid + 1
  integer, parameter :: S_L2AUX              = s_l2gp + 1
  integer, parameter :: S_MATRIX             = s_l2aux + 1
  integer, parameter :: S_MERGE              = s_matrix + 1
  integer, parameter :: S_MODULE             = s_merge + 1
  integer, parameter :: S_OUTPUT             = s_module + 1
  integer, parameter :: S_QUANTITY           = s_output + 1
  integer, parameter :: S_RADIOMETER         = s_quantity + 1
  integer, parameter :: S_RETRIEVE           = s_radiometer + 1
  integer, parameter :: S_SIDS               = s_retrieve + 1    !??? for Zvi
  integer, parameter :: S_SIGNAL             = s_sids + 1
  integer, parameter :: S_SNOOP              = s_signal + 1
  integer, parameter :: S_SPECTROMETERTYPE   = s_snoop + 1
  integer, parameter :: S_SUBSET             = s_spectrometertype + 1
  integer, parameter :: S_TEMPLATE           = s_subset + 1
  integer, parameter :: S_TIME               = s_template + 1
  integer, parameter :: S_TPFILL             = s_time + 1
  integer, parameter :: S_VECTOR             = s_tpfill + 1
  integer, parameter :: S_VECTORTEMPLATE     = s_vector + 1
  integer, parameter :: S_VGRID              = s_vectortemplate + 1
  integer, parameter :: S_L2LOAD             = s_vgrid + 1       !??? for Zvi
  integer, parameter :: SPEC_FIRST = last_parm + 1, SPEC_LAST = s_l2load
! integer, parameter :: SPEC_FIRST = last_parm + 1, SPEC_LAST = s_vGrid
  integer :: SPEC_INDICES(spec_first:spec_last)

! Table for section ordering:
  integer, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: |       globalSettings    chunkDivide       retrieve             |
!     |             readApriori       construct          join          |
!     | mlsSignals        mergeApriori       fill             output   |
! ====|================================================================|== From: ==
        (/OK,   OK,    0,    0,    0,    0,    0,    0,    0,    0,  & ! Start
           0,   OK,    0,    0,    0,    0,    0,    0,    0,    0,  & ! mlsSignals
           0,    0,   OK,    0,    0,    0,    0,    0,    0,    0,  & ! globalSettings
           0,    0,    0,   OK,    0,    0,    0,    0,    0,    0,  & ! readApriori
           0,    0,    0,    0,   OK,    0,    0,    0,    0,    0,  & ! mergeApriori
           0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,  & ! chunkDivide
           0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,  & ! Construct
           0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,  & ! Fill
           0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,  & ! Retrieve
           0,    0,    0,    0,    0,    0,    0,    0,    0,   OK,  & ! Join
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

  integer, private, parameter :: BEGIN = -1
  integer, private, parameter :: D = 1000000
  integer, private, parameter :: F = 1000, L = 2000, N = 0
  integer, private, parameter :: NADP = n+d*(all_fields+no_dup+no_positional)
  integer, private, parameter :: ND = n+d*no_dup
  integer, private, parameter :: NDP = n+d*(no_dup+no_positional)
  integer, private, parameter :: NP = n+d*no_positional
  integer, private, parameter :: NR = n+d*req_fld
  integer, private, parameter :: P = 3000, S = 4000, T = 5000, Z = 6000

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES
    use TREE_TYPES, only: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
                          N_NAME_DEF, N_SECTION, N_SPEC_DEF
    integer :: I ! used only in an array constructor as a DO index

  ! Put intrinsic predefined identifiers into the symbol table.
    call init_intrinsic ( data_type_indices, lit_indices )

  ! Put nonintrinsic predefined identifiers into the symbol table.
    ! Put enumeration type names into the symbol table
    data_type_indices(t_criticalmodule) =  add_ident ( 'criticalModule' )
    data_type_indices(t_fillmethod) =      add_ident ( 'fillMethod' )
    data_type_indices(t_griddedOrigin) =   add_ident ( 'griddedOrigin' )
    data_type_indices(t_hgridtype) =       add_ident ( 'hGridType' )
    data_type_indices(t_matrix) =          add_ident ( 'matrixType' )
    data_type_indices(t_mergemethod) =     add_ident ( 'mergeMethod' )
    data_type_indices(t_method) =          add_ident ( 'method' )
    data_type_indices(t_mergesource) =     add_ident ( 'mergeSource' )
    data_type_indices(t_module) =          add_ident ( 'module' )
    data_type_indices(t_molecule) =        add_ident ( 'molecule' )
    data_type_indices(t_outputtype) =      add_ident ( 'outputType' )
    data_type_indices(t_quantitytype) =    add_ident ( 'quantityType' )
    data_type_indices(t_radiometer) =      add_ident ( 'radiometer' )
    data_type_indices(t_scale) =           add_ident ( 'scale' )
    data_type_indices(t_species) =         add_ident ( 'species' )
    data_type_indices(t_units) =           add_ident ( 'units' )
    data_type_indices(t_vgridcoord) =      add_ident ( 'vGridCoord' )
    data_type_indices(t_vgridtype) =       add_ident ( 'vGridType' )
    ! Put enumeration literals into the symbol table:
    lit_indices(l_angle) =                 add_ident ( 'angle' )
    lit_indices(l_apriori) =               add_ident ( 'apriori' )
    lit_indices(l_both) =                  add_ident ( 'both' )
    lit_indices(l_cholesky) =              add_ident ( 'cholesky' )
    lit_indices(l_climatology) =           add_ident ( 'climatology' )
    lit_indices(l_covariance) =            add_ident ( 'covariance' )
    lit_indices(l_dao) =                   add_ident ( 'DAO' )
    lit_indices(l_direct) =                add_ident ( 'direct' )
    lit_indices(l_earthRefl) =             add_ident ( 'earthRefl' )
    lit_indices(l_either) =                add_ident ( 'either' )
    lit_indices(l_elevOffset) =            add_ident ( 'elevOffset' )
    lit_indices(l_explicit) =              add_ident ( 'explicit' )
    lit_indices(l_fixed) =                 add_ident ( 'fixed' )
    lit_indices(l_fractional) =            add_ident ( 'fractional' )
    lit_indices(l_height) =                add_ident ( 'height' )
    lit_indices(l_hydrostatic) =           add_ident ( 'hydrostatic' )
    lit_indices(l_kronecker) =             add_ident ( 'kronecker' )
    lit_indices(l_l1b) =                   add_ident ( 'l1b' )
    lit_indices(l_l2aux) =                 add_ident ( 'l2aux' )
    lit_indices(l_l2gp) =                  add_ident ( 'l2gp' )
    lit_indices(l_linear) =                add_ident ( 'linear' )
    lit_indices(l_logarithmic) =           add_ident ( 'logarithmic' )
    lit_indices(l_ncep) =                  add_ident ( 'NCEP' )
    lit_indices(l_neither) =               add_ident ( 'neither' )
    lit_indices(l_newtonian) =             add_ident ( 'newtonian' )
    lit_indices(l_norm) =                  add_ident ( 'norm' )
    lit_indices(l_orbitIncline) =          add_ident ( 'orbitIncline' )
    lit_indices(l_plain) =                 add_ident ( 'plain' )
    lit_indices(l_pressure) =              add_ident ( 'pressure' )
    lit_indices(l_r1a) =                   add_ident ( 'r1a' )
    lit_indices(l_r1b) =                   add_ident ( 'r1b' )
    lit_indices(l_r2) =                    add_ident ( 'r2' )
    lit_indices(l_r3) =                    add_ident ( 'r3' )
    lit_indices(l_r4) =                    add_ident ( 'r4' )
    lit_indices(l_r5h) =                   add_ident ( 'r5h' )
    lit_indices(l_r5v) =                   add_ident ( 'r5v' )
    lit_indices(l_scGeocAlt) =             add_ident ( 'scGeocAlt' )
    lit_indices(l_spaceRadiance) =         add_ident ( 'spaceRadiance' )
    lit_indices(l_spd) =                   add_ident ( 'spd' )
    lit_indices(l_vector) =                add_ident ( 'vector' )
    lit_indices(l_weighted) =              add_ident ( 'weighted' )
    ! Put field names into the symbol table
    field_indices(f_apriori) =             add_ident ( 'apriori' )
    field_indices(f_aprioriscale) =        add_ident ( 'aprioriScale' )
    field_indices(f_autofill) =            add_ident ( 'autofill' )
    field_indices(f_band) =                add_ident ( 'band' )
    field_indices(f_centerFrequency) =     add_ident ( 'centerFrequency' )
    field_indices(f_channel) =             add_ident ( 'channel' )
    field_indices(f_channels) =            add_ident ( 'channels' )
    field_indices(f_columns) =             add_ident ( 'columns' )
    field_indices(f_columnscale) =         add_ident ( 'columnScale' )
    field_indices(f_comment) =             add_ident ( 'comment' )
    field_indices(f_compareOverlaps) =     add_ident ( 'compareOverlaps' )
    field_indices(f_coordinate) =          add_ident ( 'coordinate' )
    field_indices(f_copy) =                add_ident ( 'copy' )
    field_indices(f_covariance) =          add_ident ( 'covariance' )
    field_indices(f_criteria) =            add_ident ( 'criteria' )
    field_indices(f_deferred) =            add_ident ( 'deferred' )
    field_indices(f_diagonal) =            add_ident ( 'diagonal' )
    field_indices(f_diagonalOut) =         add_ident ( 'diagonalOut' )
    field_indices(f_explicitvalues) =      add_ident ( 'explicitValues' )
    field_indices(f_field) =               add_ident ( 'field' )
    field_indices(f_file) =                add_ident ( 'file' )
    field_indices(f_first) =               add_ident ( 'first' )
    field_indices(f_fraction) =            add_ident ( 'fraction' )
    field_indices(f_frequencies) =         add_ident ( 'frequencies' )
    field_indices(f_frequency) =           add_ident ( 'frequency' )
    field_indices(f_fwdModelExtra) =       add_ident ( 'fwdModelExtra' )
    field_indices(f_fwdModelIn) =          add_ident ( 'fwdModelIn' )
    field_indices(f_fwdModelOut) =         add_ident ( 'fwdModelOut' )
    field_indices(f_geocAltitudeQuantity) =add_ident ( 'geocAltitudeQuantity' )
    field_indices(f_gph) =                 add_ident ( 'gph' )
    field_indices(f_h2oquantity) =         add_ident ( 'h2oquantity' )
    field_indices(f_height) =              add_ident ( 'height' )
    field_indices(f_hgrid) =               add_ident ( 'hgrid' )
    field_indices(f_interpolationFactor) = add_ident ( 'interpolationFactor' )
    field_indices(f_jacobian) =            add_ident ( 'jacobian' )
    field_indices(f_last) =                add_ident ( 'last' )
    field_indices(f_length) =              add_ident ( 'length' )
    field_indices(f_lo) =                  add_ident ( 'lo' )
    field_indices(f_matrix) =              add_ident ( 'matrix' )
    field_indices(f_maxIterations) =       add_ident ( 'maxIterations' )
    field_indices(f_measurements) =        add_ident ( 'measurements' )
    field_indices(f_method) =              add_ident ( 'method' )
    field_indices(f_mif) =                 add_ident ( 'mif' )
    field_indices(f_module) =              add_ident ( 'module' )
    field_indices(f_molecule) =            add_ident ( 'molecule' )
    field_indices(f_number) =              add_ident ( 'number' )
    field_indices(f_origin) =              add_ident ( 'origin' )
    field_indices(f_outputCovariance) =    add_ident ( 'outputCovariance' )
    field_indices(f_outputOverlaps) =      add_ident ( 'outputOverlaps' )
    field_indices(f_overlaps) =            add_ident ( 'overlaps' )
    field_indices(f_per_decade) =          add_ident ( 'per_decade' )
    field_indices(f_quantities) =          add_ident ( 'quantities' )
    field_indices(f_quantity) =            add_ident ( 'quantity' )
    field_indices(f_radiometer) =          add_ident ( 'radiometer' )
    field_indices(f_range) =               add_ident ( 'range' )
    field_indices(f_refGPHQuantity) =      add_ident ( 'refGPHquantity' )
    field_indices(f_rows) =                add_ident ( 'rows' )
    field_indices(f_scale) =               add_ident ( 'scale' )
    field_indices(f_sdname) =              add_ident ( 'sdname' )
    field_indices(f_signal) =              add_ident ( 'signal' )
    field_indices(f_signals) =             add_ident ( 'signals' )
    field_indices(f_source) =              add_ident ( 'source' )
    field_indices(f_sourceapriori) =       add_ident ( 'sourceApriori' )
    field_indices(f_sourcel2aux) =         add_ident ( 'sourceL2AUX' )
    field_indices(f_sourcel2gp) =          add_ident ( 'sourceL2GP' )
    field_indices(f_sourcequantity) =      add_ident ( 'sourceQuantity' )
    field_indices(f_spacecraft) =          add_ident ( 'spacecraft' )
    field_indices(f_species) =             add_ident ( 'species' )
    field_indices(f_spectrometer) =        add_ident ( 'spectrometer' )
    field_indices(f_spectrometerType) =    add_ident ( 'spectrometerType' )
    field_indices(f_spread) =              add_ident ( 'spread' )
    field_indices(f_start) =               add_ident ( 'start' )
    field_indices(f_state) =               add_ident ( 'state' )
    field_indices(f_step) =                add_ident ( 'step' )
    field_indices(f_stop) =                add_ident ( 'stop' )
    field_indices(f_suffix) =              add_ident ( 'suffix' )
    field_indices(f_swath) =               add_ident ( 'swath' )
    field_indices(f_switch) =              add_ident ( 'switch' )
    field_indices(f_temperature) =         add_ident ( 'temperature' )
    field_indices(f_temperaturequantity) = add_ident ( 'temperatureQuantity' )
    field_indices(f_tolerancea) =          add_ident ( 'Atolerance' )
    field_indices(f_tolerancef) =          add_ident ( 'Ftolerance' )
    field_indices(f_tolerancer) =          add_ident ( 'Rtolerance' )
    field_indices(f_versionRange) =        add_ident ( 'versionRange' )
    field_indices(f_template) =            add_ident ( 'template' )
    field_indices(f_test) =                add_ident ( 'test' )
    field_indices(f_type) =                add_ident ( 'type' )
    field_indices(f_unit) =                add_ident ( 'unit' )
    field_indices(f_values) =              add_ident ( 'values' )
    field_indices(f_vGrid) =               add_ident ( 'vgrid' )
    field_indices(f_weight) =              add_ident ( 'weight' )
    field_indices(f_width) =               add_ident ( 'width' )
    field_indices(f_widths) =              add_ident ( 'widths' )
    field_indices(f_bill) =                add_ident ( 'bill' )        !???
    field_indices(f_atmos_der) =           add_ident ( 'atmos_der' )   !???
    field_indices(f_do_conv) =             add_ident ( 'conv' )        !???
    field_indices(f_do_freq_avg) =         add_ident ( 'freq_avg' )    !???
    field_indices(f_spect_der) =           add_ident ( 'spect_der' )   !???
    field_indices(f_temp_der) =            add_ident ( 'temp_der' )    !???
    field_indices(f_zvi) =                 add_ident ( 'zvi' )         !???
    ! Put parameter names into the symbol table
    parm_indices(p_allow_climatology_overloads) = &
                                           add_ident ( 'AllowClimatologyOverloads' )
    parm_indices(p_input_version_string) = add_ident ( 'InputVersionString' )
    parm_indices(p_output_version_string) =add_ident ( 'OutputVersionString' )
    parm_indices(p_version_comment) =      add_ident ( 'VersionComment' )
    parm_indices(p_critical_bands) =       add_ident ( 'CriticalBands' )
    parm_indices(p_critical_scanning_modules) = &
                                           add_ident ( 'CriticalScanningModules' )
    parm_indices(p_home_geod_angle) =      add_ident ( 'HomeGeodAngle' )
    parm_indices(p_home_module) =          add_ident ( 'HomeModule' )
    parm_indices(p_ideal_length) =         add_ident ( 'IdealLength' )
    parm_indices(p_max_gap) =              add_ident ( 'MaxGap' )
    parm_indices(p_overlap) =              add_ident ( 'Overlap' )
    parm_indices(p_scan_lower_limit) =     add_ident ( 'ScanLowerLimit' )
    parm_indices(p_scan_upper_limit) =     add_ident ( 'ScanUpperLimit' )
    ! Put section names into the symbol table
    section_indices(z_chunkDivide) =       add_ident ( 'chunkDivide' )
    section_indices(z_construct) =         add_ident ( 'construct' )
    section_indices(z_fill) =              add_ident ( 'fill' )
    section_indices(z_globalSettings) =    add_ident ( 'globalSettings' )
    section_indices(z_join) =              add_ident ( 'join' )
    section_indices(z_mergeApriori) =      add_ident ( 'mergeApriori' )
    section_indices(z_mlsSignals) =        add_ident ( 'mlsSignals' )
    section_indices(z_output) =            add_ident ( 'output' )
    section_indices(z_readApriori) =       add_ident ( 'readApriori' )
    section_indices(z_retrieve) =          add_ident ( 'retrieve' )
    ! Put spec names into the symbol table
    spec_indices(s_apriori) =              add_ident ( 'apriori' )
    spec_indices(s_band) =                 add_ident ( 'band' )
    spec_indices(s_create) =               add_ident ( 'create' )
    spec_indices(s_fill) =                 add_ident ( 'fill' )
    spec_indices(s_forwardModel) =         add_ident ( 'forwardModel' )
    spec_indices(s_forwardModelGlobal) =   add_ident ( 'forwardModelGlobal' )
    spec_indices(s_gridded) =              add_ident ( 'gridded' )
    spec_indices(s_hgrid) =                add_ident ( 'hgrid' )
    spec_indices(s_l2gp) =                 add_ident ( 'l2gp' )
    spec_indices(s_l2aux) =                add_ident ( 'l2aux' )
    spec_indices(s_matrix) =               add_ident ( 'matrix' )
    spec_indices(s_merge) =                add_ident ( 'merge' )
    spec_indices(s_module) =               add_ident ( 'module' )
    spec_indices(s_output) =               add_ident ( 'output' )
    spec_indices(s_quantity) =             add_ident ( 'quantity' )
    spec_indices(s_radiometer) =           add_ident ( 'radiometer' )
    spec_indices(s_retrieve) =             add_ident ( 'retrieve' )
    spec_indices(s_signal) =               add_ident ( 'signal' )
    spec_indices(s_snoop) =                add_ident ( 'snoop' )
    spec_indices(s_spectrometerType) =     add_ident ( 'spectrometerType' )
    spec_indices(s_subset) =               add_ident ( 'subset' )
    spec_indices(s_template) =             add_ident ( 'template' )
    spec_indices(s_time) =                 add_ident ( 'time' )
    spec_indices(s_tpfill) =               add_ident ( 'tpfill' )
    spec_indices(s_vector) =               add_ident ( 'vector' )
    spec_indices(s_vectortemplate) =       add_ident ( 'vectorTemplate' )
    spec_indices(s_vgrid) =                add_ident ( 'vgrid' )
    spec_indices(s_l2load) =               add_ident ( 'l2load' )!??? for Zvi
    spec_indices(s_sids) =                 add_ident ( 'sids' )  !??? for Zvi

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
    ! Define the intrinsic data types
    call make_tree ( (/ &
      begin, t+t_numeric, n+n_dt_def, &
      begin, t+t_numeric_range, n+n_dt_def, &
      begin, t+t_string, n+n_dt_def /) )
    ! Define the enumerated types
    call make_tree ( (/ &
      begin, t+t_boolean, l+l_true, l+l_false, n+n_dt_def, &
      begin, t+t_griddedOrigin, l+l_climatology, l+l_dao, l+l_ncep, n+n_dt_def, &
      begin, t+t_criticalModule, l+l_both, l+l_either, l+l_ghz, l+l_neither, &
             l+l_thz, n+n_dt_def, &
      begin, t+t_fillMethod, l+l_apriori, l+l_explicit, l+l_hydrostatic, &
             l+l_l1b, l+l_l2aux, l+l_l2gp, l+l_vector, n+n_dt_def, &
      begin, t+t_hGridType, l+l_explicit, l+l_fixed, l+l_fractional, &
             l+l_height, l+l_linear, n+n_dt_def, &
      begin, t+t_matrix, l+l_plain, l+l_cholesky, l+l_kronecker, l+l_spd, &
             n+n_dt_def, &
      begin, t+t_mergeMethod, l+l_direct, l+l_weighted, n+n_dt_def, &
      begin, t+t_mergeSource, l+l_dao, l+l_ncep, n+n_dt_def, &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def, &
      begin, t+t_molecule, l+(/ (i,i=first_molecule, last_molecule) /), &
             n+n_dt_def, &
      begin, t+t_outputType, l+l_l2aux, l+l_l2gp, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_quantityType, l+l_baseline, l+l_earthRefl, l+l_elevOffset, &
             l+l_extinction, l+l_gph, l+l_orbitIncline, l+l_ptan,&
             l+l_radiance, l+l_refGPH, l+l_scVel, l+l_scGeocAlt, l+l_spaceRadiance,&
             l+l_temperature, l+l_tngtGeodAlt, l+l_tngtGeocAlt, &
             l+l_vmr, n+n_dt_def, &
      begin, t+t_radiometer, l+l_r1a, l+l_r1b, l+l_r2, l+l_r3, l+l_r4, &
             l+l_r5h, l+l_r5v, n+n_dt_def, &
      begin, t+t_scale, l+l_apriori, & ! l+l_covariance, & !??? Later !???
             l+l_none, l+l_norm, n+n_dt_def, &
      begin, t+t_species, l+l_gph, l+l_gph_precision, l+l_temperature, &
             l+l_temperature_prec, n+n_dt_def, &
      begin, t+t_units, l+l_days, l+l_deg, l+l_degrees, &
             l+l_dimensionless, l+l_dimless, l+l_dl, l+l_ghz, &
             l+l_hours, l+l_hpa, l+l_hz, l+l_k, l+l_khz, l+l_km, l+l_logp, &
             l+l_m, l+l_maf, l+l_mafs, l+l_mb, l+l_meters, l+l_mhz, &
             l+l_mif, l+l_mifs, l+l_minutes, l+l_orbits, l+l_pa, l+l_ppbv, &
             l+l_ppmv, l+l_pptv, l+l_rad, l+l_radians, l+l_s, l+l_seconds, &
             l+l_thz, l+l_vmr, l+l_zeta, n+n_dt_def, &
      begin, t+t_vgridcoord, l+l_angle, l+l_geodAltitude, l+l_gph, l+l_none, &
             l+l_pressure, l+l_theta, l+l_zeta, n+n_dt_def, &
      begin, t+t_vgridtype, l+l_explicit, l+l_linear, l+l_logarithmic, &
             n+n_dt_def /) )
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
      begin, s+s_module, &
             begin, f+f_spacecraft, t+t_boolean, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_radiometer, &          ! Must be after module
             begin, f+f_lo, t+t_numeric, n+n_field_type, &
             begin, f+f_suffix, t+t_string, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             nadp+n_spec_def, &
      begin, s+s_spectrometerType, &
             begin, f+f_deferred, t+t_boolean, n+n_field_type, &
             begin, f+f_first, t+t_numeric, n+n_field_type, &
             begin, f+f_frequencies, t+t_numeric, n+n_field_type, &
             begin, f+f_last, t+t_numeric, n+n_field_type, &
             begin, f+f_start, t+t_numeric, n+n_field_type, &
             begin, f+f_step, t+t_numeric, n+n_field_type, &
             begin, f+f_width, t+t_numeric, n+n_field_type, &
             begin, f+f_widths, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_band, &                ! Must be after radiometer and spectrometerType
             begin, f+f_suffix, t+t_string, n+n_field_type, &
             begin, f+f_spectrometerType, s+s_spectrometerType, nr+n_field_spec, &
             begin, f+f_radiometer, s+s_radiometer, nr+n_field_spec, &
             begin, f+f_centerfrequency, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_signal, &         ! Must be after band
             begin, f+f_band, s+s_band, nr+n_field_spec, &
             begin, f+f_spectrometer, t+t_numeric, nr+n_field_type, &
             begin, f+f_frequencies, t+t_numeric, n+n_field_type, &
             begin, f+f_widths, t+t_numeric, n+n_field_type, &
             begin, f+f_switch, t+t_numeric, nr+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_time, np+n_spec_def, &
      begin, s+s_gridded, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_field, t+t_string, n+n_field_type, &
             begin, f+f_origin, t+t_griddedOrigin, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_create, &
             begin, f+f_template, n+n_field_type, &
             begin, f+f_copy, n+n_field_type, &
             begin, f+f_autofill, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_hGrid, &
             begin, f+f_type, t+t_hGridType, nr+n_field_type, &
             begin, f+f_module, s+s_module, nr+n_field_spec, &
             begin, f+f_fraction, t+t_numeric, n+n_field_type, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_mif, t+t_numeric, n+n_field_type, &
             begin, f+f_interpolationfactor, t+t_numeric, n+n_field_type, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_merge, &  ! Must be AFTER S_Gridded
             begin, f+f_apriori, s+s_gridded, n+n_field_spec, &
             begin, f+f_source, t+t_mergeSource, n+n_field_type, &
             begin, f+f_species, t+t_species, n+n_field_type, &
             begin, f+f_range, t+t_numeric_range, n+n_field_type, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_method, t+t_mergeMethod, n+n_field_type, &
             begin, f+f_scale, t+t_numeric, n+n_field_type, &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_template, &
             begin, f+f_copy, n+n_field_type, &
             begin, f+f_apriori, n+n_field_type, &
             begin, f+f_autofill, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_tpfill, &
             begin, f+f_type, n+n_field_type, &
             begin, f+f_temperature, n+n_field_type, &
             begin, f+f_gph, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_vgrid, &
             begin, f+f_type, t+t_vGridType, nr+n_field_type, &
             begin, f+f_coordinate, t+t_vGridCoord, nr+n_field_type, &
             begin, f+f_number, t+t_numeric, n+n_field_type, &
             begin, f+f_per_decade, t+t_numeric, n+n_field_type, &
             begin, f+f_start, t+t_numeric, n+n_field_type, &
             begin, f+f_stop, t+t_numeric, n+n_field_type, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_quantity, & ! Must be AFTER s_hgrid and s_vgrid
             begin, f+f_hGrid, s+s_hgrid, n+n_field_spec, &
             begin, f+f_vGrid, s+s_vgrid, n+n_field_spec, &
             begin, f+f_molecule, t+t_molecule, n+n_field_type, &
             begin, f+f_radiometer, t+t_radiometer, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_signal, s+s_signal, n+n_field_spec, &
             begin, f+f_type, t+t_quantityType, n+n_field_type, &
             begin, f+f_unit, t+t_units, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_vectorTemplate, & ! Must be AFTER s_quantity
             begin, f+f_quantities, s+s_quantity, n+n_field_spec, &
             begin, f+f_signals, t+t_string, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_vector, & ! Must be AFTER s_vectorTemplate
             begin, f+f_template, s+s_vectorTemplate, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_l2gp, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_swath, t+t_string, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_l2aux, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    nr+n_dot, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_sdname, t+t_string, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_matrix, &  ! Must be AFTER s_vector
             begin, f+f_rows, s+s_vector, n+n_field_spec, &
             begin, f+f_columns, s+s_vector, n+n_field_spec, &
             begin, f+f_type, t+t_matrix, n+n_field_type, &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_fill, &    ! Must be AFTER s_vector, s_matrix and s_climatology
             begin, f+f_quantity, s+s_vector, f+f_template, f+f_quantities, &
                    nr+n_dot, &
             begin, f+f_matrix, s+s_matrix, n+n_field_spec, &
             begin, f+f_method, t+t_fillmethod, nr+n_field_type, &
             begin, f+f_sourceQuantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_temperatureQuantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_h2oQuantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_geocAltitudeQuantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_refGPHQuantity, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_sourceL2GP, s+s_l2gp, n+n_field_spec, &
             begin, f+f_sourceL2AUX, s+s_l2aux, n+n_field_spec, &
             begin, f+f_sourceApriori, s+s_apriori, n+n_field_spec, &
             begin, f+f_spread, t+t_boolean, n+n_field_type, &
             begin, f+f_maxIterations, t+t_numeric, n+n_field_type, &
             begin, f+f_explicitValues, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_output, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_type, t+t_outputType, nr+n_field_type, &
             begin, f+f_file, t+t_string, nr+n_field_type, &
             begin, f+f_quantities, s+s_l2aux, s+s_l2gp, n+n_field_spec, &
             begin, f+f_overlaps, s+s_l2aux, s+s_l2gp, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_subset, &  ! Must be AFTER s_vector
             begin, f+f_quantity, s+s_vector, f+f_template, f+f_quantities, &
                    nr+n_dot, &
             begin, f+f_test, s+s_vector, f+f_template, f+f_quantities, &
                    nr+n_dot, &
             begin, f+f_channels, t+t_numeric, n+n_field_type, &
             begin, f+f_criteria, t+t_numeric, nr+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_forwardModel, & ! Must be AFTER s_vector and s_matrix
             ndp+n_spec_def, &
      begin, s+s_forwardModelGlobal, &                                !???
             begin, f+f_atmos_der, t+t_boolean, n+n_field_type, &     !???
             begin, f+f_do_conv, t+t_boolean, n+n_field_type, &       !???
             begin, f+f_do_freq_avg, t+t_boolean, n+n_field_type, &   !???
             begin, f+f_frequency, t+t_numeric, n+n_field_type, &     !???
             begin, f+f_spect_der, t+t_boolean, n+n_field_type, &     !???
             begin, f+f_temp_der, t+t_boolean, n+n_field_type, &      !???
             ndp+n_spec_def /) )                                      !???
    call make_tree ( (/ &
      begin, s+s_retrieve, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_apriori, s+s_vector, n+n_field_spec, &
             begin, f+f_aprioriScale, t+t_numeric, n+n_field_type, &
             begin, f+f_columnScale, t+t_scale, n+n_field_type, &
             begin, f+f_covariance, s+s_matrix, n+n_field_spec, &
             begin, f+f_diagonal, t+t_boolean, n+n_field_type, &
             begin, f+f_diagonalOut, t+t_boolean, n+n_field_type, &
             begin, f+f_fwdModelExtra, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelIn, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelOut, s+s_vector, n+n_field_spec, &
             begin, f+f_jacobian, s+s_matrix, n+n_field_spec, &
             begin, f+f_maxIterations, t+t_numeric, n+n_field_type, &
             begin, f+f_measurements, s+s_vector, nr+n_field_spec, &
             begin, f+f_method, t+t_method, n+n_field_type, &
             begin, f+f_outputCovariance, s+s_matrix, n+n_field_spec, &
             begin, f+f_state, s+s_vector, nr+n_field_spec, &
             begin, f+f_toleranceA, t+t_numeric, n+n_field_type, &
             begin, f+f_toleranceF, t+t_numeric, n+n_field_type, &
             begin, f+f_toleranceR, t+t_numeric, n+n_field_type, &
             begin, f+f_weight, s+s_vector, n+n_field_spec, &
             ndp+n_spec_def, &
      begin, s+s_sids, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_fwdModelExtra, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelIn, s+s_vector, nr+n_field_spec, &
             begin, f+f_fwdModelOut, s+s_vector, nr+n_field_spec, &
             begin, f+f_jacobian, s+s_matrix, n+n_field_spec, &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_snoop, &
             begin, f+f_comment, t+t_string, n+n_field_type, &
             nd+n_spec_def /) )
    call make_tree ( (/ &                                    !???
      begin, s+s_l2load, &                                   !???
             begin, f+f_bill, t+t_string, n+n_field_type, &  !???
             begin, f+f_zvi, t+t_string, n+n_field_type, &   !???
             nadp+n_spec_def /) )                            !???
    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >
    call make_tree ( (/ &
      begin, z+z_mlsSignals, s+s_module, s+s_band, s+s_radiometer, &
                             s+s_signal, s+s_spectrometerType, &
             n+n_section, &
      begin, z+z_globalsettings, &
             begin, p+p_version_comment, t+t_string, n+n_name_def, &
             begin, p+p_input_version_string, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_allow_climatology_overloads, t+t_boolean, &
                    n+n_name_def, &
             s+s_time, &
             s+s_l2load, s+s_forwardModelGlobal, &           !???
             n+n_section, &
      begin, z+z_readapriori, s+s_time, s+s_gridded, s+s_l2gp, &
             s+s_l2aux, s+s_snoop, n+n_section, &
      begin, z+z_mergeapriori, s+s_time, s+s_merge, n+n_section, &
      begin, z+z_chunkdivide, &
             begin, p+p_critical_bands, t+t_string, n+n_name_def, &
             begin, p+p_critical_scanning_modules, t+t_criticalModule, &
                    n+n_name_def, &
             begin, p+p_home_geod_angle, t+t_numeric, n+n_name_def, &
             begin, p+p_home_module, t+t_module, n+n_name_def, &
             begin, p+p_ideal_length, t+t_numeric, n+n_name_def, &
             begin, p+p_max_gap, t+t_numeric, n+n_name_def, &
             begin, p+p_overlap, t+t_numeric, n+n_name_def, &
             begin, p+p_scan_lower_limit, t+t_numeric_range, n+n_name_def, &
             begin, p+p_scan_upper_limit, t+t_numeric_range, n+n_name_def, &
             n+n_section, &
      begin, z+z_construct, s+s_time, s+s_vgrid, s+s_hgrid, s+s_quantity, &
             s+s_vectortemplate, s+s_snoop, n+n_section, &
      begin, z+z_fill, s+s_time, s+s_vector, s+s_tpfill, s+s_create, &
                       s+s_fill, s+s_matrix, s+s_snoop, &
             n+n_section, &
      begin, z+z_retrieve, s+s_matrix, s+s_forwardModel, s+s_retrieve, &
             s+s_subset, s+s_sids, s+s_time, n+n_section, &
      begin, z+z_join, s+s_time, s+s_l2gp, s+s_l2aux, n+n_section, &
      begin, z+z_output, s+s_time, s+s_output, n+n_section /) )
  end subroutine INIT_TABLES

! =====     Private procedures     =====================================
  ! --------------------------------------------------  ADD_IDENT  -----
  integer function ADD_IDENT ( TEXT )
    use SYMBOL_TABLE, only: ENTER_TERMINAL
    use SYMBOL_TYPES, only: T_IDENTIFIER
    character(len=*), intent(in) :: TEXT
    add_ident = enter_terminal ( text, t_identifier )
  end function ADD_IDENT

  ! --------------------------------------------------  MAKE_TREE  -----
  include "make_tree.f9h"

end module INIT_TABLES_MODULE

! $Log$
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
