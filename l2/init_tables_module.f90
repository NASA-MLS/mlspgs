module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER
  use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
  use TREE_TYPES, only: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
                        N_NAME_DEF, N_SECTION, N_SPEC_DEF

  implicit NONE
  private

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  public :: INIT_TABLES

! Data types that don't have enumerated literals:
  integer, public, parameter :: T_NUMERIC = 1
  integer, public, parameter :: T_NUMERIC_RANGE = 2
  integer, public, parameter :: T_STRING = 3
! Enumeration types:
  integer, public, parameter :: T_APRIORISOURCE = 4
  integer, public, parameter :: T_APRIORITYPE = 5
  integer, public, parameter :: T_BOOLEAN = 6
  integer, public, parameter :: T_CRITICALMODULE = 7
  integer, public, parameter :: T_HGRIDTYPE = 8
  integer, public, parameter :: T_MERGEMETHOD = 9
  integer, public, parameter :: T_MERGESOURCE = 10
  integer, public, parameter :: T_MODULE = 11
  integer, public, parameter :: T_MOLECULE = 12
  integer, public, parameter :: T_OUTPUTTYPE = 13
  integer, public, parameter :: T_QUANTITYTYPE = 14
  integer, public, parameter :: T_SPECIES = 15
  integer, public, parameter :: T_UNITS = 16
  integer, public, parameter :: T_VGRIDCOORD = 17
  integer, public, parameter :: T_VGRIDTYPE = 18
  integer, public, parameter :: T_FIRST = T_NUMERIC, T_LAST = T_VGRIDTYPE
  integer, public :: DATA_TYPE_INDICES(t_first:t_last)
! Field indices:
  integer, public, parameter :: F_APRIORI = 1
  integer, public, parameter :: F_AUTOFILL = 2
  integer, public, parameter :: F_BAND = 3
  integer, public, parameter :: F_COMPAREOVERLAPS = 4
  integer, public, parameter :: F_COORDINATE = 5
  integer, public, parameter :: F_COPY = 6
  integer, public, parameter :: F_FILE = 7
  integer, public, parameter :: F_FIRSTINDEXCHANNEL = 8
  integer, public, parameter :: F_FRACTION = 9
  integer, public, parameter :: F_GPH = 10
  integer, public, parameter :: F_HDFNAME = 11
  integer, public, parameter :: F_HEIGHT = 12
  integer, public, parameter :: F_HGRID = 13
  integer, public, parameter :: F_INTERPOLATIONFACTOR = 14
  integer, public, parameter :: F_LENGTH = 15
  integer, public, parameter :: F_METHOD = 16
  integer, public, parameter :: F_MIF = 17
  integer, public, parameter :: F_MODULE = 18
  integer, public, parameter :: F_MOLECULE = 19
  integer, public, parameter :: F_NUMBER = 20
  integer, public, parameter :: F_OUTPUTOVERLAPS = 21
  integer, public, parameter :: F_OVERLAPS = 22
  integer, public, parameter :: F_PER_DECADE = 23
  integer, public, parameter :: F_QUANTITIES = 24
  integer, public, parameter :: F_RADIOMETER = 25
  integer, public, parameter :: F_RANGE = 26
  integer, public, parameter :: F_SCALE = 27
  integer, public, parameter :: F_SIGNALS = 28
  integer, public, parameter :: F_SOURCE = 29
  integer, public, parameter :: F_SPECIES = 30
  integer, public, parameter :: F_START = 31
  integer, public, parameter :: F_STOP = 32
  integer, public, parameter :: F_TEMPERATURE = 33
  integer, public, parameter :: F_VERSIONRANGE = 34
  integer, public, parameter :: F_TEMPLATE = 35
  integer, public, parameter :: F_TYPE = 36
  integer, public, parameter :: F_UNIT = 37
  integer, public, parameter :: F_UNPACKOUTPUT = 38
  integer, public, parameter :: F_VALUES = 39
  integer, public, parameter :: F_VGRID = 40
  integer, public, parameter :: FIELD_FIRST = f_Apriori, FIELD_LAST = f_vGrid
  integer, public :: FIELD_INDICES(field_first:field_last)
! Enumeration literals:
  integer, public, parameter :: L_ANGLE = 1
  integer, public, parameter :: L_BASELINE      = l_angle + 1
  integer, public, parameter :: L_BOTH 	        = l_baseline + 1
  integer, public, parameter :: L_CLIMATOLOGY   = l_both + 1
  integer, public, parameter :: L_CLO           = l_climatology + 1
  integer, public, parameter :: L_CO            = l_clo + 1
  integer, public, parameter :: L_DAO 	        = l_co + 1
  integer, public, parameter :: L_DAYS 	        = l_dao + 1
  integer, public, parameter :: L_DEG 	        = l_days + 1
  integer, public, parameter :: L_DEGREES       = l_deg + 1
  integer, public, parameter :: L_DIMENSIONLESS = l_degrees + 1
  integer, public, parameter :: L_DIMLESS       = l_dimensionless + 1
  integer, public, parameter :: L_DIRECT        = l_dimless + 1
  integer, public, parameter :: L_DL 	        = l_direct + 1
  integer, public, parameter :: L_EITHER        = l_dl + 1
  integer, public, parameter :: L_EXPLICIT      = l_either + 1
  integer, public, parameter :: L_EXTINCTION    = l_explicit + 1
  integer, public, parameter :: L_FALSE         = l_extinction + 1
  integer, public, parameter :: L_FIXED         = l_false + 1
  integer, public, parameter :: L_FRACTIONAL    = l_fixed + 1
  integer, public, parameter :: L_GEODALTITUDE  = l_fractional + 1
  integer, public, parameter :: L_GHZ           = l_geodaltitude + 1
  integer, public, parameter :: L_GPH 	        = l_ghz + 1
  integer, public, parameter :: L_GPH_PRECISION = l_gph + 1
  integer, public, parameter :: L_H2O           = l_gph_precision + 1
  integer, public, parameter :: L_HCL           = l_h2o + 1
  integer, public, parameter :: L_HEIGHT        = l_hcl + 1
  integer, public, parameter :: L_HNO3          = l_height + 1
  integer, public, parameter :: L_HOURS         = l_hno3 + 1
  integer, public, parameter :: L_HPA 	        = l_hours + 1
  integer, public, parameter :: L_HZ 	        = l_hpa + 1
  integer, public, parameter :: L_K 	        = l_hz + 1
  integer, public, parameter :: L_KHZ 	        = l_k  + 1
  integer, public, parameter :: L_KM 	        = l_khz + 1
  integer, public, parameter :: L_L2AUX         = l_km + 1
  integer, public, parameter :: L_L2GP 	        = l_l2aux + 1
  integer, public, parameter :: L_LINEAR        = l_l2gp + 1
  integer, public, parameter :: L_LOGARITHMIC   = l_linear + 1
  integer, public, parameter :: L_LOGP 	        = l_logarithmic + 1
  integer, public, parameter :: L_M 	        = l_logp + 1
  integer, public, parameter :: L_MAF 	        = l_m + 1
  integer, public, parameter :: L_MAFS 	        = l_maf + 1
  integer, public, parameter :: L_MB 	        = l_mafs + 1
  integer, public, parameter :: L_METERS        = l_mb + 1
  integer, public, parameter :: L_MHZ 	        = l_meters  + 1
  integer, public, parameter :: L_MIF 	        = l_mhz + 1
  integer, public, parameter :: L_MIFS 	        = l_mif + 1
  integer, public, parameter :: L_MINUTES       = l_mifs + 1
  integer, public, parameter :: L_N2O           = l_minutes + 1
  integer, public, parameter :: L_NCEP 	        = l_n2o + 1
  integer, public, parameter :: L_NEITHER       = l_ncep + 1
  integer, public, parameter :: L_NONE 	        = l_neither + 1
  integer, public, parameter :: L_O3            = l_none + 1
  integer, public, parameter :: L_ORBITS        = l_o3 + 1
  integer, public, parameter :: L_PA 	        = l_orbits + 1
  integer, public, parameter :: L_PPBV 	        = l_pa + 1
  integer, public, parameter :: L_PPMV 	        = l_ppbv + 1
  integer, public, parameter :: L_PPTV 	        = l_ppmv + 1
  integer, public, parameter :: L_PRESSURE      = l_pptv + 1
  integer, public, parameter :: L_PTAN 	        = l_pressure + 1
  integer, public, parameter :: L_RAD 	        = l_ptan + 1
  integer, public, parameter :: L_RADIANCE      = l_rad  + 1
  integer, public, parameter :: L_RADIANS       = l_radiance + 1
  integer, public, parameter :: L_S 	        = l_radians + 1
  integer, public, parameter :: L_SECONDS       = l_s + 1
  integer, public, parameter :: L_TEMPERATURE   = l_seconds + 1
  integer, public, parameter :: L_TEMPERATURE_PREC = l_temperature + 1
  integer, public, parameter :: L_THETA         = l_temperature_prec + 1
  integer, public, parameter :: L_THZ 	        = l_theta + 1
  integer, public, parameter :: L_TRUE 	        = l_thz + 1
  integer, public, parameter :: L_VMR 	        = l_true + 1
  integer, public, parameter :: L_WEIGHTED      = l_vmr + 1
  integer, public, parameter :: L_ZETA 	        = l_weighted + 1
  integer, public, parameter :: FIRST_LIT = L_ANGLE, LAST_LIT = L_ZETA
  integer, public :: LIT_INDICES(first_lit:last_lit)
! Parameter names:
  ! In GlobalSettings section:
  integer, public, parameter :: P_ALLOW_CLIMATOLOGY_OVERLOADS = 1
  integer, public, parameter :: P_INPUT_VERSION_STRING = 2
  integer, public, parameter :: P_OUTPUT_VERSION_STRING = 3
  integer, public, parameter :: P_VERSION_COMMENT = 4
  ! In ChunkDivide section:
  integer, public, parameter :: P_CRITICAL_BANDS = 5
  integer, public, parameter :: P_CRITICAL_SCANNING_MODULES = 6
  integer, public, parameter :: P_HOME_GEOD_ANGLE = 7
  integer, public, parameter :: P_HOME_MODULE = 8
  integer, public, parameter :: P_IDEAL_LENGTH = 9
  integer, public, parameter :: P_MAX_GAP = 10
  integer, public, parameter :: P_OVERLAP = 11
  integer, public, parameter :: P_SCAN_LOWER_LIMIT = 12
  integer, public, parameter :: P_SCAN_UPPER_LIMIT = 13
  integer, public, parameter :: FIRST_PARM = P_ALLOW_CLIMATOLOGY_OVERLOADS
  integer, public, parameter :: LAST_PARM = P_SCAN_UPPER_LIMIT
  integer, public :: PARM_INDICES(first_parm:last_parm)
! Abstract physical quantities:
  integer, public, parameter :: PHYQ_INVALID = 0 ! Invalid unit given by user
  integer, public, parameter :: PHYQ_DIMENSIONLESS = 1 ! Dimensionless quantity
  integer, public, parameter :: PHYQ_LENGTH = 2        ! Default meters
  integer, public, parameter :: PHYQ_TIME = 3          ! Default seconds
  integer, public, parameter :: PHYQ_PRESSURE = 4      ! Default millibars
  integer, public, parameter :: PHYQ_TEMPERATURE = 5   ! Default Kelvins
  integer, public, parameter :: PHYQ_VMR = 6           ! Default parts-per-one
  integer, public, parameter :: PHYQ_ANGLE = 7         ! Default degrees
  integer, public, parameter :: PHYQ_MAFS = 8          ! Default MAFs
  integer, public, parameter :: PHYQ_MIFS = 9          ! Default MIFs
  integer, public, parameter :: PHYQ_FREQUENCY = 10    ! Default MHz
  integer, public, parameter :: PHYQ_ZETA = 11         ! log10(pressure/hPa)
  integer, public, parameter :: FIRST_PHYQ = phyq_invalid, LAST_PHYQ = phyq_zeta
  integer, public :: PHYQ_INDICES(first_phyq:last_phyq)
! Section identities:
  integer, public, parameter :: Z_CHUNKDIVIDE = 4
  integer, public, parameter :: Z_CONSTRUCT = 5
  integer, public, parameter :: Z_FILL = 6
  integer, public, parameter :: Z_GLOBALSETTINGS = 1
  integer, public, parameter :: Z_JOIN = 7
  integer, public, parameter :: Z_MERGEAPRIORI = 3
  integer, public, parameter :: Z_OUTPUT = 8
  integer, public, parameter :: Z_READAPRIORI = 2
  integer, public, parameter :: SECTION_FIRST = z_globalSettings, &
                                SECTION_LAST = z_Output
  integer, public :: SECTION_INDICES(section_first:section_last)
! Specification indices:
  integer, public, parameter :: S_CLIMATOLOGY = 1
  integer, public, parameter :: S_CREATE = 2
  integer, public, parameter :: S_HGRID = 3
  integer, public, parameter :: S_L2GP = 4
  integer, public, parameter :: S_L2AUX = 5
  integer, public, parameter :: S_MERGE = 6
  integer, public, parameter :: S_OUTPUT = 7
  integer, public, parameter :: S_QUANTITY = 8
  integer, public, parameter :: S_TEMPLATE = 9
  integer, public, parameter :: S_TPFILL = 10
  integer, public, parameter :: S_VECTOR = 11
  integer, public, parameter :: S_VECTORTEMPLATE = 12
  integer, public, parameter :: S_VGRID = 13
  integer, public, parameter :: SPEC_FIRST = s_Climatology, SPEC_LAST = s_vGrid
  integer, public :: SPEC_INDICES(spec_first:spec_last)

! Table for section ordering:
  integer, public, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: | globalSettings    chunkDivide        join          |
!     |       readApriori       construct         output   |
!     |             mergeApriori       fill                |
! ====|====================================================|== From: ==
        (/OK,    0,    0,    0,    0,    0,    0,    0,  & ! Start
           0,   OK,    0,    0,    0,    0,    0,    0,  & ! GlobalSettings
           0,    0,   OK,    0,    0,    0,    0,    0,  & ! readApriori
           0,    0,    0,   OK,    0,    0,    0,    0,  & ! mergeApriori
           0,    0,    0,    0,   OK,   OK,   OK,   OK,  & ! chunkDivide
           0,    0,    0,    0,   OK,   OK,   OK,   OK,  & ! Construct
           0,    0,    0,    0,   OK,   OK,   OK,   OK,  & ! Fill
           0,    0,    0,    0,    0,    0,    0,   OK,  & ! Join
           0,    0,    0,    0,    0,    0,    0,    0/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

  integer, private, parameter :: F = 1000, L = 2000, N = 0, P = 3000
  integer, private, parameter :: S = 4000, T = 5000, Z = 6000
  integer, private, parameter :: BEGIN = -1

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES
  ! Put predefined identifiers into the symbol table.
    ! Put non-enumeration type names into symbol table
    data_type_indices(t_numeric) =         add_ident ( 'numeric' )
    data_type_indices(t_numeric_range) =   add_ident ( 'numeric_range' )
    data_type_indices(t_string) =          add_ident ( 'string' )
    ! Put enumeration type names into the symbol table
    data_type_indices(t_apriorisource) =   add_ident ( 'aprioriSource' )
    data_type_indices(t_aprioritype) =     add_ident ( 'aprioriType' )
    data_type_indices(t_boolean) =         add_ident ( 'boolean' )
    data_type_indices(t_criticalmodule) =  add_ident ( 'criticalModule' )
    data_type_indices(t_hgridtype) =       add_ident ( 'hGridType' )
    data_type_indices(t_mergemethod) =     add_ident ( 'mergeMethod' )
    data_type_indices(t_mergesource) =     add_ident ( 'mergeSource' )
    data_type_indices(t_module) =          add_ident ( 'module' )
    data_type_indices(t_molecule) =        add_ident ( 'molecule' )
    data_type_indices(t_outputtype) =      add_ident ( 'outputType' )
    data_type_indices(t_quantitytype) =    add_ident ( 'quantityType' )
    data_type_indices(t_species) =         add_ident ( 'species' )
    data_type_indices(t_units) =           add_ident ( 'units' )
    data_type_indices(t_vgridcoord) =      add_ident ( 'vGridCoord' )
    data_type_indices(t_vgridtype) =       add_ident ( 'vGridType' )
    ! Put field names into the symbol table
    field_indices(f_apriori) =             add_ident ( 'apriori' )
    field_indices(f_autofill) =            add_ident ( 'autofill' )
    field_indices(f_band) =                add_ident ( 'band' )
    field_indices(f_compareoverlaps) =     add_ident ( 'compareoverlaps' )
    field_indices(f_coordinate) =          add_ident ( 'coordinate' )
    field_indices(f_copy) =                add_ident ( 'copy' )
    field_indices(f_file) =                add_ident ( 'file' )
    field_indices(f_firstindexchannel) =   add_ident ( 'f_firstindexchannel' )
    field_indices(f_fraction) =            add_ident ( 'fraction' )
    field_indices(f_gph) =                 add_ident ( 'gph' )
    field_indices(f_hdfname) =             add_ident ( 'hdfname' )
    field_indices(f_height) =              add_ident ( 'height' )
    field_indices(f_hgrid) =               add_ident ( 'hgrid' )
    field_indices(f_interpolationfactor) = add_ident ( 'interpolationfactor' )
    field_indices(f_length) =              add_ident ( 'length' )
    field_indices(f_method) =              add_ident ( 'method' )
    field_indices(f_mif) =                 add_ident ( 'mif' )
    field_indices(f_module) =              add_ident ( 'module' )
    field_indices(f_molecule) =            add_ident ( 'molecule' )
    field_indices(f_number) =              add_ident ( 'number' )
    field_indices(f_outputoverlaps) =      add_ident ( 'outputoverlaps' )
    field_indices(f_overlaps) =            add_ident ( 'overlaps' )
    field_indices(f_per_decade) =          add_ident ( 'per_decade' )
    field_indices(f_quantities) =          add_ident ( 'quantities' )
    field_indices(f_radiometer) =          add_ident ( 'radiometer' )
    field_indices(f_range) =               add_ident ( 'range' )
    field_indices(f_scale) =               add_ident ( 'scale' )
    field_indices(f_signals) =             add_ident ( 'signals' )
    field_indices(f_source) =              add_ident ( 'source' )
    field_indices(f_species) =             add_ident ( 'species' )
    field_indices(f_start) =               add_ident ( 'start' )
    field_indices(f_stop) =                add_ident ( 'stop' )
    field_indices(f_temperature) =         add_ident ( 'temperature' )
    field_indices(f_versionrange) =        add_ident ( 'versionrange' )
    field_indices(f_template) =            add_ident ( 'template' )
    field_indices(f_type) =                add_ident ( 'type' )
    field_indices(f_unit) =                add_ident ( 'unit' )
    field_indices(f_unpackoutput) =        add_ident ( 'unpackoutput' )
    field_indices(f_values) =              add_ident ( 'values' )
    field_indices(f_vgrid) =               add_ident ( 'vgrid' )
    ! Put enumeration literals into the symbol table:
    lit_indices(l_angle) =                 add_ident ( 'angle' )
    lit_indices(l_baseline) =              add_ident ( 'baseline' )
    lit_indices(l_both) =                  add_ident ( 'both' )
    lit_indices(l_climatology) =           add_ident ( 'climatology' )
    lit_indices(l_clo) =                   add_ident ( 'clo' )
    lit_indices(l_co) =                    add_ident ( 'co' )
    lit_indices(l_dao) =                   add_ident ( 'DAO' )
    lit_indices(l_days) =                  add_ident ( 'days' )
    lit_indices(l_deg) =                   add_ident ( 'deg' )
    lit_indices(l_degrees) =               add_ident ( 'degrees' )
    lit_indices(l_dimensionless) =         add_ident ( 'dimensionless' )
    lit_indices(l_dimless) =               add_ident ( 'dimless' )
    lit_indices(l_direct) =                add_ident ( 'direct' )
    lit_indices(l_dl) =                    add_ident ( 'dl' )
    lit_indices(l_either) =                add_ident ( 'either' )
    lit_indices(l_explicit) =              add_ident ( 'explicit' )
    lit_indices(l_extinction) =            add_ident ( 'extinction' )
    lit_indices(l_false) =                 add_ident ( 'false' )
    lit_indices(l_fixed) =                 add_ident ( 'fixed' )
    lit_indices(l_fractional) =            add_ident ( 'fractional' )
    lit_indices(l_geodaltitude) =          add_ident ( 'geodaltitude' )
    lit_indices(l_ghz) =                   add_ident ( 'GHz' )
    lit_indices(l_gph) =                   add_ident ( 'gph' )
    lit_indices(l_gph_precision) =         add_ident ( 'gph_precision' )
    lit_indices(l_h2o) =                   add_ident ( 'h2o' )
    lit_indices(l_hcl) =                   add_ident ( 'hcl' )
    lit_indices(l_height) =                add_ident ( 'height' )
    lit_indices(l_hno3) =                  add_ident ( 'hno3' )
    lit_indices(l_hours) =                 add_ident ( 'hours' )
    lit_indices(l_hpa) =                   add_ident ( 'hPa' )
    lit_indices(l_hz) =                    add_ident ( 'Hz' )
    lit_indices(l_k) =                     add_ident ( 'k' )
    lit_indices(l_khz) =                   add_ident ( 'KHz' )
    lit_indices(l_km) =                    add_ident ( 'km' )
    lit_indices(l_l2aux) =                 add_ident ( 'l2aux' )
    lit_indices(l_l2gp) =                  add_ident ( 'l2gp' )
    lit_indices(l_linear) =                add_ident ( 'linear' )
    lit_indices(l_logarithmic) =           add_ident ( 'logarithmic' )
    lit_indices(l_logp) =                  add_ident ( 'logp' )
    lit_indices(l_m) =                     add_ident ( 'm' )
    lit_indices(l_maf) =                   add_ident ( 'maf' )
    lit_indices(l_mafs) =                  add_ident ( 'mafs' )
    lit_indices(l_mb) =                    add_ident ( 'mb' )
    lit_indices(l_meters) =                add_ident ( 'meters' )
    lit_indices(l_mhz) =                   add_ident ( 'MHz' )
    lit_indices(l_mif) =                   add_ident ( 'mif' )
    lit_indices(l_mifs) =                  add_ident ( 'mifs' )
    lit_indices(l_minutes) =               add_ident ( 'minutes' )
    lit_indices(l_n2o) =                   add_ident ( 'n2o' )
    lit_indices(l_ncep) =                  add_ident ( 'NCEP' )
    lit_indices(l_neither) =               add_ident ( 'neither' )
    lit_indices(l_none) =                  add_ident ( 'none' )
    lit_indices(l_o3) =                    add_ident ( 'o3' )
    lit_indices(l_orbits) =                add_ident ( 'orbits' )
    lit_indices(l_pa) =                    add_ident ( 'pa' )
    lit_indices(l_ppbv) =                  add_ident ( 'ppbv' )
    lit_indices(l_ppmv) =                  add_ident ( 'ppmv' )
    lit_indices(l_pptv) =                  add_ident ( 'pptv' )
    lit_indices(l_pressure) =              add_ident ( 'pressure' )
    lit_indices(l_ptan) =                  add_ident ( 'ptan' )
    lit_indices(l_rad) =                   add_ident ( 'rad' )
    lit_indices(l_radiance) =              add_ident ( 'radiance' )
    lit_indices(l_radians) =               add_ident ( 'radians' )
    lit_indices(l_s) =                     add_ident ( 's' )
    lit_indices(l_seconds) =               add_ident ( 'seconds' )
    lit_indices(l_temperature) =           add_ident ( 'temperature' )
    lit_indices(l_temperature_prec) =      add_ident ( 'temperature_precision' )
    lit_indices(l_theta) =                 add_ident ( 'theta' )
    lit_indices(l_thz) =                   add_ident ( 'THz' )
    lit_indices(l_true) =                  add_ident ( 'true' )
    lit_indices(l_vmr) =                   add_ident ( 'vmr' )
    lit_indices(l_weighted) =              add_ident ( 'weighted' )
    lit_indices(l_zeta) =                  add_ident ( 'zeta' )
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
    ! Put abstract physical quantities into the symbol table
    phyq_indices(phyq_invalid) =           add_ident ( 'invalid' )
    phyq_indices(phyq_dimensionless) =     add_ident ( 'dimensionless' )
    phyq_indices(phyq_length) =            add_ident ( 'length' )
    phyq_indices(phyq_time) =              add_ident ( 'time' )
    phyq_indices(phyq_pressure) =          add_ident ( 'pressure' )
    phyq_indices(phyq_temperature) =       add_ident ( 'temperature' )
    phyq_indices(phyq_vmr) =               add_ident ( 'vmr' )
    phyq_indices(phyq_angle) =             add_ident ( 'angle' )
    phyq_indices(phyq_mafs) =              add_ident ( 'mafs' )
    phyq_indices(phyq_mifs) =              add_ident ( 'mifs' )
    phyq_indices(phyq_frequency) =         add_ident ( 'frequency' )
    phyq_indices(phyq_zeta) =              add_ident ( 'zeta' )
    ! Put section names into the symbol table
    section_indices(z_chunkdivide) =       add_ident ( 'chunkdivide' )
    section_indices(z_construct) =         add_ident ( 'construct' )
    section_indices(z_fill) =              add_ident ( 'fill' )
    section_indices(z_globalsettings) =    add_ident ( 'globalsettings' )
    section_indices(z_join) =              add_ident ( 'join' )
    section_indices(z_mergeapriori) =      add_ident ( 'mergeapriori' )
    section_indices(z_output) =            add_ident ( 'output' )
    section_indices(z_readapriori) =       add_ident ( 'readapriori' )
    ! Put spec names into the symbol table
    spec_indices(s_climatology) =          add_ident ( 'climatology' )
    spec_indices(s_create) =               add_ident ( 'create' )
    spec_indices(s_hgrid) =                add_ident ( 'hgrid' )
    spec_indices(s_l2gp) =                 add_ident ( 'l2gp' )
    spec_indices(s_l2aux) =                add_ident ( 'l2aux' )
    spec_indices(s_merge) =                add_ident ( 'merge' )
    spec_indices(s_output) =               add_ident ( 'output' )
    spec_indices(s_quantity) =             add_ident ( 'quantity' )
    spec_indices(s_template) =             add_ident ( 'template' )
    spec_indices(s_tpfill) =               add_ident ( 'tpfill' )
    spec_indices(s_vector) =               add_ident ( 'vector' )
    spec_indices(s_vectortemplate) =       add_ident ( 'vectortemplate' )
    spec_indices(s_vgrid) =                add_ident ( 'vgrid' )

  ! Definitions are represented by trees.  The notation in the comments
  ! for the trees is < root first_son ... last_son >.  This is sometimes
  ! called "Cambridge Polish Notation."  It was developed to represent
  ! LISP by McCarthy et. al. at MIT (in Cambridge, MA).

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
      begin, t+t_aprioriSource, l+l_clo, l+l_co, l+l_gph, l+l_gph_precision, &
             l+l_h2o, l+l_hcl, l+l_hno3, l+l_n2o, l+l_o3, l+l_temperature, &
             l+l_temperature_prec, n+n_dt_def, &
      begin, t+t_aprioritype, l+l_climatology, l+l_l2gp, l+l_l2aux, &
             n+n_dt_def, &
      begin, t+t_boolean, l+l_true, l+l_false, n+n_dt_def, &
      begin, t+t_criticalModule, l+l_both, l+l_either, l+l_ghz, l+l_neither, &
             l+l_thz, n+n_dt_def, &
      begin, t+t_hGridType, l+l_explicit, l+l_fixed, l+l_fractional, &
             l+l_height, l+l_linear, n+n_dt_def, &
      begin, t+t_mergeMethod, l+l_direct, l+l_weighted, n+n_dt_def, &
      begin, t+t_mergeSource, l+l_dao, l+l_ncep, n+n_dt_def, &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def, &
      begin, t+t_molecule, l+l_clo, l+l_co, l+l_gph, l+l_h2o, &
             l+l_hcl, l+l_hno3, l+l_n2o, l+l_o3, n+n_dt_def, &
      begin, t+t_outputType, l+l_l2aux, l+l_l2gp, n+n_dt_def, &
      begin, t+t_quantityType, l+l_baseline, l+l_extinction, l+l_gph, &
             l+l_ptan, l+l_radiance, l+l_temperature, l+l_vmr, n+n_dt_def, &
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
    ! The n_field_type subtree indicates the type allowed for a field.
    ! The n_field_spec subtree indicates the specifications whose names
    ! are allowed to appear for a field.
    ! The n_dot subtree indicates that the field given by the first
    ! f_field_name  is required to be of the form spec_name.field_name,
    ! where spec_name is required to be a label of a specification of the
    ! type given by the s_spec son, and field_name is required to be
    ! present in the field given by the last f_field_name, which is
    ! required to be in a specification named by the next-to-last
    ! f_field_name ... of the specification named by the spec_name.
    call make_tree ( (/ &
      begin, s+s_climatology, &
             begin, f+f_type, t+t_aprioriType, n+n_field_type, &
             begin, f+f_source, t+t_aprioriSource, n+n_field_type, &
             begin, f+f_length, t+t_numeric, n+n_field_type, &
             begin, f+f_versionRange, t+t_string, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_create, &
             begin, f+f_template, n+n_field_type, &
             begin, f+f_copy, n+n_field_type, &
             begin, f+f_autofill, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_hGrid, &
             begin, f+f_type, t+t_hGridType, n+n_field_type, &
             begin, f+f_module, t+t_module, n+n_field_type, &
             begin, f+f_fraction, t+t_numeric, n+n_field_type, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_mif, t+t_numeric, n+n_field_type, &
             begin, f+f_interpolationfactor, t+t_numeric, n+n_field_type, &
             begin, f+f_values, n+n_field_type, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_merge, &  ! Must be AFTER s_climatology
             begin, f+f_apriori, s+s_climatology, n+n_field_spec, &
             begin, f+f_source, t+t_mergeSource, n+n_field_type, &
             begin, f+f_species, t+t_species, n+n_field_type, &
             begin, f+f_range, t+t_numeric_range, n+n_field_type, &
             begin, f+f_height, t+t_numeric, n+n_field_type, &
             begin, f+f_method, t+t_mergeMethod, n+n_field_type, &
             begin, f+f_scale, t+t_numeric, n+n_field_type, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_template, &
             begin, f+f_copy, n+n_field_type, &
             begin, f+f_apriori, n+n_field_type, &
             begin, f+f_autofill, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_tpfill, &
             begin, f+f_type, n+n_field_type, &
             begin, f+f_temperature, n+n_field_type, &
             begin, f+f_gph, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_vgrid, &
             begin, f+f_type, t+t_vGridType, n+n_field_type, &
             begin, f+f_coordinate, t+t_vGridCoord, n+n_field_type, &
             begin, f+f_number, t+t_numeric, n+n_field_type, &
             begin, f+f_per_decade, t+t_numeric, n+n_field_type, &
             begin, f+f_start, t+t_numeric, n+n_field_type, &
             begin, f+f_stop, t+t_numeric, n+n_field_type, &
             begin, f+f_values, t+t_numeric, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_quantity, & ! Must be AFTER s_hgrid and s_vgrid
             begin, f+f_band, t+t_string, n+n_field_type, &
             begin, f+f_firstindexchannel, t+t_boolean, n+n_field_type, &
             begin, f+f_hGrid, s+s_hgrid, n+n_field_spec, &
             begin, f+f_vGrid, s+s_vgrid, n+n_field_spec, &
             begin, f+f_molecule, t+t_molecule, n+n_field_type, &
             begin, f+f_radiometer, t+t_string, n+n_field_type, &
             begin, f+f_type, t+t_quantityType, n+n_field_type, &
             begin, f+f_unit, t+t_units, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_vectorTemplate, & ! Must be AFTER s_quantity
             begin, f+f_quantities, s+s_quantity, n+n_field_spec, &
             begin, f+f_signals, t+t_string, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_vector, & ! Must be AFTER s_vectorTemplate
             begin, f+f_template, s+s_vectorTemplate, n+n_field_spec, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_l2gp, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_unpackOutput, t+t_boolean, n+n_field_type, &
             begin, f+f_hdfname, t+t_string, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_l2aux, &   ! Must be AFTER s_vector
             begin, f+f_source, s+s_vector, f+f_template, f+f_quantities, &
                    n+n_dot, &
             begin, f+f_compareOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_outputOverlaps, t+t_boolean, n+n_field_type, &
             begin, f+f_unpackOutput, t+t_boolean, n+n_field_type, &
             begin, f+f_hdfname, t+t_string, n+n_field_type, &
             n+n_spec_def, &
      begin, s+s_output, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_type, t+t_outputType, n+n_field_type, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_quantities, s+s_l2aux, s+s_l2gp, n+n_field_spec, &
             begin, f+f_overlaps, s+s_l2aux, s+s_l2gp, n+n_field_spec, &
             n+n_spec_def /) )
    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >
    call make_tree ( (/ &
      begin, z+z_globalsettings, &
             begin, p+p_version_comment, t+t_string, n+n_name_def, &
             begin, p+p_input_version_string, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_allow_climatology_overloads, t+t_boolean, &
                    n+n_name_def, &
             n+n_section, &
      begin, z+z_readapriori, s+s_climatology, n+n_section, &
      begin, z+z_mergeapriori, s+s_merge, n+n_section, &
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
      begin, z+z_construct, s+s_vgrid, s+s_hgrid, s+s_quantity, &
             s+s_vectortemplate, n+n_section, &
      begin, z+z_fill, s+s_vector, s+s_tpfill, s+s_create, n+n_section, &
      begin, z+z_join, s+s_l2gp, s+s_l2aux, n+n_section, &
      begin, z+z_output, s+s_output, n+n_section /) )
  end subroutine INIT_TABLES

! =====     Private procedures     =====================================
  ! --------------------------------------------------  MAKE_TREE  -----
  subroutine MAKE_TREE ( IDS )
  ! Build a tree specified by the "ids" array.  "begin" marks the
  ! beginning of a tree.  A tree-node marks the end of the corresponding
  ! tree.  Pseudo-terminals are decorated with their indices.
    integer, intent(in) :: IDS(:)

    integer :: I, ID, M, N_IDS, STACK(0:30), STRING, WHICH

    n_ids = size(ids)
    m = 0
    stack(0) = 0 ! just so it's defined, in case it gets incremented
                 ! after build_tree
    if ( ids(1) >= 0 ) then
      m = 1
      stack(1) = 0
    end if
    do i = 1, n_ids
      if ( ids(i) == begin ) then
        m = m + 1
        stack(m) = 0
      else
        id = mod(ids(i), 1000)
        which = ids(i) / 1000
       select case ( which )
       case ( f/1000 ) ! Fields
         string = field_indices(id)
       case ( l/1000 ) ! Enumeration literals
         string = lit_indices(id)
       case ( p/1000 ) ! Parameter names
         string = parm_indices(id)
       case ( s/1000 ) ! Specs
         string = spec_indices(id)
       case ( t/1000 ) ! Intrinsic data types
         string = data_type_indices(id)
       case ( z/1000 ) ! Sections
         string = section_indices(id)
       case ( n/1000 ) ! Tree nodes
         call build_tree ( id, stack(m) )
         m = m - 1
         stack(m) = stack(m) + 1
    cycle
       end select
       if ( string == 0 ) then
         print *, 'INIT_TABLES_MODULE%MAKE_TREE-E- Element ', i, &
           & ' of a list is undefined'
           stop
         end if
       call push_pseudo_terminal ( string, 0, decor = id )
       stack(m) = stack(m) + 1
      end if
    end do
  end subroutine MAKE_TREE

  integer function ADD_IDENT ( TEXT )
    character(len=*), intent(in) :: TEXT
    add_ident = enter_terminal ( text, t_identifier )
  end function ADD_IDENT

end module INIT_TABLES_MODULE

! $Log$
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

