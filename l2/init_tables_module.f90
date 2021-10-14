
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

! Declaring the definitions is handled by the tree walker.

  use Init_MLSSignals_M ! Everything. Init_MLSSignals, Field_First,
    ! Last_Signal_Field, Spec_First, Last_Signal_Spec, Numerous S_....
  use Init_Spectroscopy_M ! Everything.
  use Intrinsic ! Everything. ADD_IDENT, BEGIN, D, F, FIRST_LIT,
    ! INIT_INTRINSIC, L, L_<several>, LAST_INTRINSIC_LIT,
    ! N, NADP, ND, NDP, NO_CHECK_EQ, NP, NR, P, S, T, <all>_INDICES,
    ! T_BOOLEAN, T_FIRST, T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE,
    ! T_STRING,s_TIME and Z are used here, but everything is included so
    ! that it can be gotten by USE INIT_TABLES_MODULE.
  use Molecules ! Everything.
  use Units, only: Init_Units
  ! We're adding the following use statements to clue the makefiles'
  ! dependency calculator for srclib's tree_checker
  use Tree, only: ! build_tree, push_pseudo_terminal
  
  implicit none

  public ! This would be a MUCH LONGER list than the list of private
  !        names below.
  private :: ADD_IDENT
  private :: Field_Spec, Field_Spec_Array, Field_Spec_1, Field_Spec_2
  private :: Field_Spec_3, Field_Spec_4, Field_Spec_5
  private :: Field_Type, Field_Type_Array, Field_Type_1, Field_Type_2
  private :: INIT_SPECTROSCOPY, Node

  interface Field_Spec
    module procedure Field_Spec_1, Field_Spec_2
    module procedure Field_Spec_3, Field_Spec_4, Field_Spec_5
  end interface

  interface Field_Type
    module procedure Field_Type_Array, Field_Type_1, Field_Type_2
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Enumeration types:
  integer, parameter :: T_BINSELECTORTYPE = Last_spectroscopy_type+1
  integer, parameter :: T_BOXCARMETHOD = t_binSelectorType+1
  integer, parameter :: T_CHUNKDIVIDEMETHOD = t_boxCarMethod+1
  integer, parameter :: T_CLOUD_DER      = t_chunkDivideMethod+1
  integer, parameter :: T_CRITICALMODULE = t_cloud_der+1
  integer, parameter :: T_FGRIDCOORD     = t_criticalmodule+1
  integer, parameter :: T_FILLMETHOD     = t_fGridCoord+1
  integer, parameter :: T_FWMTYPE        = t_fillmethod+1
  integer, parameter :: T_GEOLOCATION    = t_fwmType+1
  integer, parameter :: T_GRIDDEDORIGIN  = t_geolocation+1
  integer, parameter :: T_HGRIDCOORD     = t_griddedOrigin+1
  integer, parameter :: T_HGRIDTYPE      = t_hGridCoord+1
  integer, parameter :: T_I_SATURATION   = t_hgridtype+1
  integer, parameter :: T_MASKS          = t_i_saturation+1
  integer, parameter :: T_MASKUPDATES    = t_masks+1
  integer, parameter :: T_MATRIX         = t_maskUpdates+1
  integer, parameter :: T_METHOD         = t_matrix+1
  integer, parameter :: T_MIFTangent     = t_method+1
  integer, parameter :: T_MODULE         = t_MIFTangent+1
  integer, parameter :: T_OUTPUTTYPE     = t_module+1
  integer, parameter :: T_QUANTITYTYPE   = t_outputtype+1
  integer, parameter :: T_REFLECTOR      = t_quantitytype+1
  integer, parameter :: T_ROWSORCOLUMNS  = t_reflector+1
  integer, parameter :: T_SCALE          = t_rowsOrColumns+1
  integer, parameter :: T_SPECIES        = t_scale+1
  integer, parameter :: T_TGRIDCOORD     = t_species+1
  integer, parameter :: T_TGRIDTYPE      = t_tgridcoord+1
  integer, parameter :: T_Trapezoid      = t_tgridtype+1
  integer, parameter :: T_UNITS          = t_trapezoid+1
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
  integer, parameter :: Z_ALGEBRA        = 11
  integer, parameter :: Z_CHUNKDIVIDE    = 6
  integer, parameter :: Z_CONSTRUCT      = 7
  integer, parameter :: Z_FILL           = 8
  integer, parameter :: Z_GLOBALSETTINGS = 3
  integer, parameter :: Z_JOIN           = 10
  integer, parameter :: Z_MERGEGRIDS     = 5
  integer, parameter :: Z_MLSSIGNALS     = 1
  integer, parameter :: Z_OUTPUT         = 12
  integer, parameter :: Z_READAPRIORI    = 4
  integer, parameter :: Z_RETRIEVE       = 9
  integer, parameter :: Z_SPECTROSCOPY   = 2
  integer, parameter :: SECTION_FIRST = z_mlsSignals, &
                      & SECTION_LAST = z_Output
! Specification indices don't overlap parameter indices, so a section can
! have both parameters and specifications:
!  W a r n i n g   W a r n i n g   W a r n i n g   W a r n i n g   
! Beware of adding spec indices
! The NAG compiler will generate code that has memory problems
  integer, parameter :: S_ANYGOODRADIANCES   = last_Spectroscopy_Spec + 1
  integer, parameter :: S_ANYGOODVALUES      = s_anygoodradiances + 1
  integer, parameter :: S_BINSELECTOR        = s_anygoodvalues + 1
  integer, parameter :: S_BOOLEAN            = s_binSelector + 1
  integer, parameter :: S_CASE               = s_boolean + 1
  integer, parameter :: S_CATCHWARNING       = s_case + 1
  integer, parameter :: S_CATENATE           = s_catchwarning + 1
  integer, parameter :: S_CHANGESETTINGS     = s_catenate + 1
  integer, parameter :: S_CHECKPOINT         = s_changeSettings + 1
  integer, parameter :: S_CHUNKDIVIDE        = s_checkpoint + 1
  integer, parameter :: S_CHUNKENDSBEFORE    = s_chunkdivide + 1
  integer, parameter :: S_CHUNKSTARTSAFTER   = s_chunkendsbefore + 1
  integer, parameter :: S_COLUMNSCALE        = s_chunkstartsafter + 1
  integer, parameter :: S_COMBINECHANNELS    = s_columnScale + 1
  integer, parameter :: S_COMPARE            = s_combinechannels + 1
  integer, parameter :: S_COMPUTETOTALPOWER  = s_compare + 1
  integer, parameter :: S_CONCATENATE        = s_computeTotalPower + 1
  integer, parameter :: S_CONCATENATEGRIDS   = s_concatenate + 1
  integer, parameter :: S_CONVERTETATOP      = s_concatenateGrids + 1
  integer, parameter :: S_COPY               = s_ConvertEtaToP + 1
  integer, parameter :: S_CYCLICJACOBI       = s_copy + 1
  integer, parameter :: S_DELETE             = s_cyclicJacobi + 1
  integer, parameter :: S_DESTROY            = s_delete + 1
  integer, parameter :: S_DIFF               = s_destroy + 1
  integer, parameter :: S_DIRECTREAD         = s_diff + 1
  integer, parameter :: S_DIRECTWRITE        = s_directRead + 1
  integer, parameter :: S_DIRECTWRITEFILE    = s_directWrite + 1
  integer, parameter :: S_DISJOINTEQUATIONS  = s_directWriteFile + 1
  integer, parameter :: S_DUMP               = s_disjointEquations + 1
  integer, parameter :: S_DUMPBLOCKS         = s_dump + 1
  integer, parameter :: S_EMPIRICALGEOMETRY  = s_dumpblocks + 1
  integer, parameter :: S_ENDSELECT          = s_empiricalGeometry + 1
  integer, parameter :: S_EXECUTE            = s_endSelect + 1
  integer, parameter :: S_FGRID              = s_execute + 1
  integer, parameter :: S_FILL               = s_fGrid + 1
  integer, parameter :: S_FILLCOVARIANCE     = s_fill + 1
  integer, parameter :: S_FILLDIAGONAL       = s_fillcovariance + 1
  integer, parameter :: S_FLAGCLOUD          = s_filldiagonal + 1
  integer, parameter :: S_FLUSHL2PCBINS      = s_flagCloud + 1
  integer, parameter :: S_FLUSHPFA           = s_flushL2PCBins + 1
  integer, parameter :: S_FORGE              = s_flushPFA + 1
  integer, parameter :: S_FORWARDMODEL       = s_forge + 1
  integer, parameter :: S_FORWARDMODELGLOBAL = s_forwardModel + 1
  integer, parameter :: S_GRIDDED            = s_forwardModelGlobal + 1
  integer, parameter :: S_HESSIAN            = s_gridded + 1
  integer, parameter :: S_HGRID              = s_hessian + 1
  integer, parameter :: S_IsFileAbsent       = s_hgrid + 1
  integer, parameter :: S_ISGRIDEMPTY        = S_IsFileAbsent + 1
  integer, parameter :: S_ISSWATHEMPTY       = s_isGridEmpty + 1
  integer, parameter :: S_L1BRAD             = s_isSwathEmpty + 1
  integer, parameter :: S_L1BOA              = s_l1brad + 1
  integer, parameter :: S_L2AUX              = s_l1boa + 1
  integer, parameter :: S_L2GP               = s_l2aux + 1
  integer, parameter :: S_L2PARSF            = s_l2gp + 1
  integer, parameter :: S_LABEL              = s_l2parsf + 1
  integer, parameter :: S_LEAKCHECK          = s_label + 1
  integer, parameter :: S_LOAD               = s_leakcheck + 1
  integer, parameter :: S_MAKEPFA            = s_load + 1
  integer, parameter :: S_MATRIX             = s_makepfa + 1
  integer, parameter :: S_MERGE              = s_matrix + 1
  integer, parameter :: S_MERGEGRIDS         = s_merge + 1
  integer, parameter :: S_NEGATIVEPRECISION  = s_mergeGrids + 1
  integer, parameter :: S_NEURALNET          = s_negativePrecision + 1
  integer, parameter :: S_NORMALEQUATIONS    = s_neuralNet + 1
  integer, parameter :: S_OUTPUT             = s_normalEquations + 1
  integer, parameter :: S_PFADATA            = s_output + 1
  integer, parameter :: S_PHASE              = s_pfadata + 1
  integer, parameter :: S_POPULATEL2PCBIN    = s_phase + 1
  integer, parameter :: S_QUANTITY           = s_populateL2pcBin + 1
  integer, parameter :: S_READGRIDDEDDATA    = s_quantity + 1
  integer, parameter :: S_READPFA            = s_readGriddeddata + 1
  integer, parameter :: S_REEVALUATE         = s_readPFA + 1
  integer, parameter :: S_REFLECT            = s_reevaluate + 1
  integer, parameter :: S_REGULARIZATION     = s_reflect + 1
  integer, parameter :: S_REPEAT             = s_regularization + 1
  integer, parameter :: S_RESTRICTRANGE      = s_repeat + 1
  integer, parameter :: S_RETRIEVE           = s_restrictRange + 1
  integer, parameter :: S_ROWSCALE           = s_retrieve + 1
  integer, parameter :: S_SELECT             = s_rowscale + 1
  integer, parameter :: S_SIDS               = s_select + 1
  integer, parameter :: S_SKIP               = s_sids + 1
  integer, parameter :: S_SLEEP              = s_skip + 1
  integer, parameter :: S_SNOOP              = s_sleep + 1
  integer, parameter :: S_STREAMLINEHESSIAN  = s_snoop + 1
  integer, parameter :: S_SUBSET             = s_streamlineHessian + 1
  integer, parameter :: S_TGRID              = s_subset + 1
  integer, parameter :: S_TRANSFER           = s_tgrid + 1
  integer, parameter :: S_UPDATEMASK         = s_transfer + 1
  integer, parameter :: S_VECTOR             = s_updateMask + 1
  integer, parameter :: S_VECTORTEMPLATE     = s_vector + 1
  integer, parameter :: S_VGRID              = s_vectortemplate + 1
  integer, parameter :: S_WMOTROP            = s_vGrid + 1
  integer, parameter :: S_WMOTROPFROMGRIDS   = s_wmoTrop + 1
  integer, parameter :: S_WRITEFILEATTRIBUTE = S_wmoTropFromGrids + 1
  integer, parameter :: S_SUBSAMPLE          = S_writeFileAttribute + 1
  integer, parameter :: S_WRITEPFA           = S_subSample + 1
  integer, parameter :: SPEC_LAST = s_writePFA 

! Parameter names:
  ! In GlobalSettings section:
  integer, parameter :: FIRST_PARM = spec_last + 1
  integer, parameter :: P_BRIGHTOBJECTS               = first_parm
  integer, parameter :: P_CYCLE                       = p_brightobjects + 1
  integer, parameter :: P_ENDTIME                     = p_cycle + 1
  integer, parameter :: P_IGRF_FILE                   = p_endtime + 1
  integer, parameter :: P_INSTRUMENT                  = p_igrf_file + 1
  integer, parameter :: P_LEAPSECFILE                 = p_instrument + 1
  integer, parameter :: P_OUTPUT_VERSION_STRING       = p_leapsecfile + 1
  integer, parameter :: P_PFAFILE                     = p_output_version_string + 1
  integer, parameter :: P_STARTTIME                   = p_pfafile + 1
  integer, parameter :: LAST_PARM                     = p_starttime

! Table for section ordering:
  integer, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: |             globalSettings    chunkDivide       retrieve          output   |
!     | mlsSignals        readApriori       construct          join                |
!     |       spectroscopy      mergeGrids         fill             algebra        |
! ====|============================================================================|== From: ==
        (/OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! Start
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! mlsSignals
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! spectroscopy
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! globalSettings
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! readApriori
          OK,   OK,   OK,   OK,   OK,   OK,    0,    0,    0,    0,   OK,    0,  & ! mergeGrids
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! chunkDivide
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! Construct
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! Fill
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! Retrieve
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! Join
           0,    0,    0,    0,    0,    0,   OK,   OK,   OK,   OK,   OK,   OK,  & ! Algebra
           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   OK/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES

    ! "Tree" is here because "make depends" can't see it in make_tree
    ! (because of the "include"):
    use Tree, only:
    use Tree_types, only: n_dot, n_dt_def, n_field_type, n_name_def, &
                          n_or, n_section, n_spec_def, n_variable_ref

    logical, parameter :: Empty = .true.  ! To specify fields allowed to be empty
    logical, parameter :: Req = .true.    ! To specify required fields
    logical, parameter :: Scalar = .true. ! To specify scalar fields

  ! Reset anything to automatic parameters  
    LAST_AUTO_LIT = LAST_LIT                
  ! Put intrinsic predefined identifiers into the symbol table.
    call init_Spectroscopy ( t_last, field_last, last_lit, &
    & first_parm, last_parm, section_last, spec_last )

  ! Put nonintrinsic predefined identifiers into the symbol table.
    ! Put enumeration type names into the symbol table
    data_type_indices(t_binSelectorType) =   add_ident ( 'binSelectorType' )
    data_type_indices(t_boxCarMethod) =      add_ident ( 'boxCarMethod' )
    data_type_indices(t_chunkDivideMethod) = add_ident ( 'chunkDivideMethod' )
    data_type_indices(t_cloud_der) =         add_ident ( 'cloud_Der' )
    data_type_indices(t_criticalmodule) =    add_ident ( 'criticalModule' )
    data_type_indices(t_fgridcoord) =        add_ident ( 'fGridCoord' )
    data_type_indices(t_fillmethod) =        add_ident ( 'fillMethod' )
    data_type_indices(t_fwmType) =           add_ident ( 'fwmType' )
    data_type_indices(t_geolocation) =       add_ident ( 'geolocation' )
    data_type_indices(t_griddedOrigin) =     add_ident ( 'griddedOrigin' )
    data_type_indices(t_hgridcoord) =        add_ident ( 'hGridCoord' )
    data_type_indices(t_hgridtype) =         add_ident ( 'hGridType' )
    data_type_indices(t_i_saturation) =      add_ident ( 'i_saturation' )
    data_type_indices(t_masks) =             add_ident ( 'masks' )
    data_type_indices(t_maskUpdates) =       add_ident ( 'maskUpdates' )
    data_type_indices(t_matrix) =            add_ident ( 'matrixType' )
    data_type_indices(t_method) =            add_ident ( 'method' )
    data_type_indices(t_MIFTangent) =        add_ident ( 'miftangent' )
    data_type_indices(t_module) =            add_ident ( 'module' )
    data_type_indices(t_outputtype) =        add_ident ( 'outputType' )
    data_type_indices(t_quantitytype) =      add_ident ( 'quantityType' )
    data_type_indices(t_reflector) =         add_ident ( 'reflector' )
    data_type_indices(t_rowsOrColumns) =     add_ident ( 'rowsOrColumns' )
    data_type_indices(t_scale) =             add_ident ( 'scale' )
    data_type_indices(t_species) =           add_ident ( 'species' )
    data_type_indices(t_tgridcoord) =        add_ident ( 'tGridCoord' )
    data_type_indices(t_tgridtype) =         add_ident ( 'tGridType' )
    data_type_indices(t_trapezoid) =         add_ident ( 'trapezoid' )
    data_type_indices(t_units) =             add_ident ( 'units' )
    data_type_indices(t_vgridcoord) =        add_ident ( 'vGridCoord' )
    data_type_indices(t_vgridtype) =         add_ident ( 'vGridType' )
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
    parm_indices(p_brightObjects) =          add_ident ( 'BrightObjects' )
    parm_indices(p_cycle) =                  add_ident ( 'Cycle' )
    parm_indices(p_endtime) =                add_ident ( 'EndTime' )
    parm_indices(p_igrf_file) =              add_ident ( 'IGRF_file' )
    parm_indices(p_instrument) =             add_ident ( 'Instrument' )
    parm_indices(p_leapsecfile) =            add_ident ( 'LeapSecFile' )
    parm_indices(p_output_version_string) =  add_ident ( 'OutputVersionString' )
    parm_indices(p_PFAfile) =                add_ident ( 'PFAFile' )
    parm_indices(p_starttime) =              add_ident ( 'StartTime' )

    ! Put section names into the symbol table:
    section_indices(z_algebra) =             add_ident ( 'algebra' )
    section_indices(z_chunkDivide) =         add_ident ( 'chunkDivide' )
    section_indices(z_construct) =           add_ident ( 'construct' )
    section_indices(z_fill) =                add_ident ( 'fill' )
    section_indices(z_globalSettings) =      add_ident ( 'globalSettings' )
    section_indices(z_join) =                add_ident ( 'join' )
    section_indices(z_mergeGrids) =          add_ident ( 'mergeGrids' )
    section_indices(z_mlsSignals) =          add_ident ( 'mlsSignals' )
    section_indices(z_output) =              add_ident ( 'output' )
    section_indices(z_readApriori) =         add_ident ( 'readApriori' )
    section_indices(z_retrieve) =            add_ident ( 'retrieve' )
    section_indices(z_spectroscopy) =        add_ident ( 'spectroscopy' )
    ! Put spec names into the symbol table.  Don't add ones that are
    ! put in by init_MLSSignals.
    spec_indices(s_anygoodradiances) =       add_ident ( 'anyGoodRadiances' )
    spec_indices(s_anygoodvalues) =          add_ident ( 'anyGoodValues' )
    spec_indices(s_binSelector) =            add_ident ( 'binSelector' )
    spec_indices(s_Boolean) =                add_ident ( 'boolean' )
    spec_indices(s_case) =                   add_ident ( 'case' )
    spec_indices(s_catchwarning) =           add_ident ( 'catchwarning' )
    spec_indices(s_catenate) =               add_ident ( 'catenate' )
    spec_indices(s_changeSettings) =         add_ident ( 'changeSettings' )
    spec_indices(s_checkpoint) =             add_ident ( 'checkpoint' )
    spec_indices(s_chunkDivide) =            add_ident ( 'chunkDivide' )
    spec_indices(s_chunkEndsBefore) =        add_ident ( 'chunkEndsBefore' )
    spec_indices(s_chunkStartsAfter) =       add_ident ( 'chunkStartsAfter' )
    spec_indices(s_columnScale) =            add_ident ( 'columnScale' )
    spec_indices(s_combineChannels) =        add_ident ( 'combineChannels' )
    spec_indices(s_compare) =                add_ident ( 'compare' )
    spec_indices(s_computetotalpower) =      add_ident ( 'computeTotalPower' )
    spec_indices(s_concatenate) =            add_ident ( 'concatenate' )
    spec_indices(s_concatenateGrids) =       add_ident ( 'concatenateGrids' )
    spec_indices(s_ConvertEtaToP) =          add_ident ( 'ConvertEtaToP' )
    spec_indices(s_copy)   =                 add_ident ( 'copy' )
    spec_indices(s_cyclicJacobi) =           add_ident ( 'cyclicJacobi' )
    spec_indices(s_empiricalGeometry) =      add_ident ( 'EmpiricalGeometry' )
    spec_indices(s_delete) =                 add_ident ( 'delete' )
    spec_indices(s_destroy) =                add_ident ( 'destroy' )
    spec_indices(s_diff) =                   add_ident ( 'diff' )
    spec_indices(s_directRead) =             add_ident ( 'directRead' )
    spec_indices(s_directWrite) =            add_ident ( 'directWrite' )
    spec_indices(s_directWriteFile) =        add_ident ( 'directWriteFile' )
    spec_indices(s_disjointEquations) =      add_ident ( 'disjointEquations' )
    spec_indices(s_dump) =                   add_ident ( 'dump' )
    spec_indices(s_dumpblocks) =             add_ident ( 'dumpblocks' )
    spec_indices(s_endSelect) =              add_ident ( 'endSelect' )
    spec_indices(s_execute) =                add_ident ( 'execute' )
    spec_indices(s_fGrid) =                  add_ident ( 'fGrid' )
    spec_indices(s_fill) =                   add_ident ( 'fill' )
    spec_indices(s_fillCovariance) =         add_ident ( 'fillCovariance' )
    spec_indices(s_fillDiagonal)   =         add_ident ( 'fillDiagonal' )
    spec_indices(s_flagCloud) =              add_ident ( 'flagCloud' )
    spec_indices(s_flushL2PCBins) =          add_ident ( 'flushL2PCBins' )
    spec_indices(s_flushPFA) =               add_ident ( 'flushPFA' )
    spec_indices(s_forge) =                  add_ident ( 'forge' )
    spec_indices(s_forwardModel) =           add_ident ( 'forwardModel' )
    spec_indices(s_forwardModelGlobal) =     add_ident ( 'forwardModelGlobal' )
    spec_indices(s_gridded) =                add_ident ( 'gridded' )
    spec_indices(s_hessian) =                add_ident ( 'hessian' )
    spec_indices(s_hgrid) =                  add_ident ( 'hgrid' )
    spec_indices(s_isFileAbsent) =           add_ident ( 'isFileAbsent' )
    spec_indices(s_isGridEmpty) =            add_ident ( 'isGridEmpty' )
    spec_indices(s_isSwathEmpty) =           add_ident ( 'isSwathEmpty' )
    spec_indices(s_l1brad) =                 add_ident ( 'l1brad' )
    spec_indices(s_l1boa) =                  add_ident ( 'l1boa' )
    spec_indices(s_l2aux) =                  add_ident ( 'l2aux' )
    spec_indices(s_l2gp) =                   add_ident ( 'l2gp' )
    spec_indices(s_l2parsf) =                add_ident ( 'l2parsf' )
    spec_indices(s_label) =                  add_ident ( 'label' )
    spec_indices(s_leakcheck) =              add_ident ( 'leakcheck' )
    spec_indices(s_load) =                   add_ident ( 'load' )
    spec_indices(s_makepfa) =                add_ident ( 'makePFA' )
    spec_indices(s_matrix) =                 add_ident ( 'matrix' )
    spec_indices(s_merge) =                  add_ident ( 'merge' )
    spec_indices(s_mergeGrids) =             add_ident ( 'mergeGrids' )
    spec_indices(s_negativePrecision ) =     add_ident ( 'negativePrecision' )
    spec_indices(s_neuralNet         ) =     add_ident ( 'neuralNet' )
    spec_indices(s_normalEquations ) =       add_ident ( 'normalEquations' )
    spec_indices(s_output) =                 add_ident ( 'output' )
    spec_indices(s_pfaData) =                add_ident ( 'pfaData' )
    spec_indices(s_phase) =                  add_ident ( 'phase' )
    spec_indices(s_populateL2PCBin) =        add_ident ( 'populateL2PCBin' )
    spec_indices(s_quantity) =               add_ident ( 'quantity' )
    spec_indices(s_readpfa) =                add_ident ( 'readPFA' )
    spec_indices(s_readGriddedData) =        add_ident ( 'readGriddeddata' )
    spec_indices(s_reevaluate) =             add_ident ( 'reevaluate' )
    spec_indices(s_reflect) =                add_ident ( 'reflect' )
    spec_indices(s_regularization) =         add_ident ( 'regularization' )
    spec_indices(s_repeat) =                 add_ident ( 'repeat' )
    spec_indices(s_restrictRange) =          add_ident ( 'restrictRange' )
    spec_indices(s_retrieve) =               add_ident ( 'retrieve' )
    spec_indices(s_rowScale) =               add_ident ( 'rowScale' )
    spec_indices(s_select) =                 add_ident ( 'select' )
    spec_indices(s_sids) =                   add_ident ( 'sids' )
    spec_indices(s_skip) =                   add_ident ( 'skip' )
    spec_indices(s_sleep) =                  add_ident ( 'sleep' )
    spec_indices(s_snoop) =                  add_ident ( 'snoop' )
    spec_indices(s_streamlineHessian) =      add_ident ( 'streamlineHessian' )
    spec_indices(s_subset) =                 add_ident ( 'subset' )
    spec_indices(s_tgrid ) =                 add_ident ( 'tGrid' )
    spec_indices(s_transfer) =               add_ident ( 'transfer' )
    spec_indices(s_updateMask) =             add_ident ( 'updateMask' )
    spec_indices(s_vector) =                 add_ident ( 'vector' )
    spec_indices(s_vectortemplate) =         add_ident ( 'vectorTemplate' )
    spec_indices(s_vgrid) =                  add_ident ( 'vGrid' )
    spec_indices(s_wmoTrop) =                add_ident ( 'wmoTrop' )
    spec_indices(s_wmoTropfromGrids) =       add_ident ( 'wmoTropFromGrids' )
    spec_indices(s_writeFileAttribute) =     add_ident ( 'writeFileAttribute' )
    spec_indices(s_subSample) =              add_ident ( 'subSample' )
    spec_indices(s_writePFA) =               add_ident ( 'writePFA' )

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
      begin, t+t_binSelectorType, l+l_vmr, l+l_temperature, l+l_latitude, &
             l+l_fieldStrength, l+l_fieldElevation, l+l_fieldAzimuth, &
             l+l_nameFragment, l+l_sza, l+l_TScat, n+n_dt_def, &
      begin, t+t_boxCarMethod, l+l_max, l+l_mean, l+l_min, n+n_dt_def, &
      begin, t+t_chunkDivideMethod, l+l_fixed, l+l_even, l+l_orbital, l+l_polygon, l+l_PE, n+n_dt_def, &
      begin, t+t_cloud_der, l+l_iwc_low_height, l+l_iwc_high_height, l+l_iwp, &
             l+l_none,n+n_dt_def, &
      begin, t+t_criticalModule, l+l_both, l+l_either, l+l_ghz, l+l_none, &
             l+l_thz, n+n_dt_def, &
      begin, t+t_fGridCoord, l+l_channel, l+l_frequency, l+l_LSBFrequency, l+l_USBFrequency, &
             l+l_IntermediateFrequency, n+n_dt_def, &
      begin, t+t_fillMethod, &
             l+l_addnoise, l+l_applyBaseline, l+l_ascenddescend, l+l_asciiFile, &
             l+l_binMax, l+l_binMean, l+l_binMin, l+l_binTotal, l+l_boxcar, &
             l+l_chiSqChan, l+l_chiSqMMAF, l+l_chiSqMMIF, l+l_chiSqRatio, &
             l+l_columnAbundance, l+l_combineChannels, &
             l+l_convergenceRatio, l+l_derivative, l+l_estimatedNoise, &
             l+l_explicit, l+l_extractChannel, l+l_fold, l+l_fwdModelMean, &
             l+l_fwdModelStdDev, l+l_fwdModelTiming, &
             l+l_gather, l+l_geoidData, l+l_geoLocation, &
             l+l_gphPrecision, l+l_gphResetToGeoid, l+l_gridded, &
             l+l_H2OFromRHI, l+l_H2OPrecisionFromRHI, l+l_hydrostatic, &
             l+l_isotope, l+l_iwcfromextinction, l+l_l1b, l+l_l2aux, l+l_l2gp, &
             l+l_losVel, l+l_lsGlobal, l+l_lsLocal, l+l_lsWeighted, &
             l+l_magAzEl, l+l_magneticModel, l+l_manipulate, &
             l+l_modifyTemplate, l+l_negativePrecision, l+l_noRadsPerMIF, &
             l+l_offsetRadiance, l+l_phaseTiming, l+l_profile, &
             l+l_quality, l+l_rectanglefromlos, l+l_reflectorTempModel, &
             l+l_refract, l+l_heightFromPressure, l+l_resetUnusedRadiances, l+l_RHIFromH2O, &
             l+l_RHIPrecisionFromH2O, l+l_rotateField, l+l_scaleOverlaps, l+l_scatter, &
             l+l_sectionTiming, l+l_splitSideband, l+l_spreadChannel, &
             l+l_status, l+l_swapValues, l+l_uncompressRadiance, &
             l+l_vector, l+l_vGrid, l+l_wmoTropopause, &
             n+n_dt_def, &
      begin, t+t_fwmType, l+l_baseline, l+l_linear, l+l_full, &
             l+l_cloudFull, l+l_hybrid, l+l_scan, l+l_scan2d, l+l_switchingMirror, &
             l+l_polarLinear, n+n_dt_def, &
      begin, t+t_geolocation, l+l_ECR, l+l_geocentric, l+l_geodetic, l+l_none, &
             n+n_dt_def, &
      begin, t+t_griddedOrigin, l+l_climatology, l+l_dao, l+l_geos5, &
             l+l_geos5_7, l+l_gloria, l+l_merra, l+l_merra_2, &
             l+l_ncep, l+l_none, l+l_strat, l+l_surfaceHeight, n+n_dt_def, &
      begin, t+t_hGridType, l+l_explicit, l+l_fixed, l+l_fractional, &
             l+l_height, l+l_l2gp, l+l_QTM, l+l_regular, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_i_saturation, l+l_clear, l+l_clear_110rh_below_top, &
             l+l_clear_0rh, l+l_clear_lowest_0_110rh, &
             l+l_clear_110rh_below_tropopause, l+l_cloudy_110rh_below_top, &
             l+l_cloudy_110rh_in_cloud, l+l_cloudy_nearside_only, n+n_dt_def, &
      begin, t+t_masks, l+l_cloud, l+l_fill, l+l_full_derivatives, l+l_ignore, &
             l+l_linalg, l+l_spare, l+l_tikhonov, n+n_dt_def, &
      begin, t+t_maskUpdates, l+l_andMasks, l+l_copy, l+l_invert, l+l_orMasks, &
             n+n_dt_def, &
      begin, t+t_matrix, l+l_plain, l+l_cholesky, l+l_kronecker, l+l_spd, &
             n+n_dt_def, &
      begin, t+t_method, l+l_highcloud,l+l_lowcloud, l+l_newtonian, &
             l+l_simple, n+n_dt_def, &
      begin, t+t_MIFTangent, l+l_ECRtoFOV, l+l_ptan, n+n_dt_def, &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def, &
      begin, t+t_outputType, l+l_ascii, l+l_hdf, l+l_l2aux, l+l_l2cf, &
             l+l_l2dgg, l+l_l2fwm, l+l_l2gp, l+l_l2pc, l+l_quantity, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_quantityType, &
             l+l_adopted, l+l_baseline, l+l_boundarypressure, &
             l+l_calSidebandFraction, l+l_chisqbinned, l+l_chisqchan, &
             l+l_chisqmmaf, l+l_chisqmmif, l+l_cloudExtinction, &
             l+l_cloudIce, l+l_cloudInducedRadiance, l+l_cloudMinMax, &
             l+l_cloudRadSensitivity, l+l_cloudTemperature, &
             l+l_cloudWater, l+l_columnabundance, l+l_dnwt_abandoned, &
             l+l_dnwt_ajn, l+l_dnwt_axmax, l+l_dnwt_cait, &
             l+l_dnwt_chiSqMinNorm, l+l_dnwt_chiSqNorm, &
             l+l_dnwt_chiSqRatio, l+l_dnwt_count, l+l_dnwt_diag, &
             l+l_dnwt_dxdx, l+l_dnwt_dxdxl, l+l_dnwt_dxn, l+l_dnwt_dxnl, &
             l+l_dnwt_flag, l+l_dnwt_fnmin, l+l_dnwt_fnorm, l+l_dnwt_gdx, &
             l+l_dnwt_gfac, l+l_dnwt_gradn, l+l_dnwt_sq, l+l_dnwt_sqt, &
             l+l_earthradius, l+l_earthRefl, l+l_ECRtoFOV, &
             l+l_effectiveOpticalDepth, l+l_elevOffset, l+l_extinction, &
             l+l_extinctionv2, l+l_fieldAzimuth, l+l_fieldElevation, &
             l+l_fieldStrength, l+l_fwdModelMean, l+l_fwdModelStdDev, &
             l+l_fwdModelTiming, l+l_geodAltitude, l+l_geolocation, &
             l+l_gph, l+l_GHzAzim, l+l_heightOffset, l+l_isotopeRatio, l+l_IWC, &
             l+l_jacobian_cols, l+l_jacobian_rows, l+l_l1bMAFBaseline, &
             l+l_l1bMIF_TAI, l+l_limbSidebandFraction, l+l_lineCenter, &
             l+l_lineWidth, l+l_lineWidth_tDep, l+l_losTransFunc, &
             l+l_losVel, l+l_lowestRetrievedPressure, l+l_magneticField, &
             l+l_massMeanDiameterIce, l+l_massMeanDiameterWater, &
             l+l_MIFDeadTime, l+l_MIFExtinction, &
             l+l_MIFExtinctionExtrapolation, l+l_MIFExtinctionForm, &
             l+l_MIFExtinctionv2, l+l_mifRHI, l+l_noiseBandwidth, &
             l+l_noRadsBinned, l+l_noRadsPerMIF, l+l_numGrad, l+l_numJ, &
             l+l_numNewt, l+l_opticalDepth, l+l_orbitInclination, l+l_AscDescMode, l+l_geoHeight, &
             l+l_phaseTiming, l+l_phiTan, l+l_ptan, l+l_quality, l+l_radiance, &
             l+l_reflrefl, l+l_reflspill, l+l_refltemp, l+l_refltrans, &
             l+l_refGPH, l+l_rhi, l+l_scanResidual, l+l_scatteringAngle, &
             l+l_scECI, l+l_instECR, l+l_scECR, l+l_scGeocAlt, l+l_scVelECI, l+l_scVelECR, &
             l+l_singleChannelRadiance, l+l_sizedistribution, &
             l+l_spaceRadiance, l+l_status, l+l_strayRadiance, &
             l+l_surfaceHeight, l+l_surfacetype, l+l_systemTemperature, &
             l+l_temperature, l+l_tngtECI, l+l_tngtECR, l+l_tngtGeocAlt, &
             l+l_tngtGeodAlt, l+l_tngtGeodLat, l+l_totalExtinction, &
             l+l_totalPowerWeight, l+l_TScat, l+l_vmr, &
             n+n_dt_def , &
      begin, t+t_reflector, l+l_primary, l+l_secondary, l+l_tertiary, &
             l+l_complete, n+n_dt_def, &
      begin, t+t_rowsOrColumns, l+l_rows, l+l_columns, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, t+t_scale, l+l_apriori, & ! l+l_covariance, & !??? Later !???
             l+l_none, l+l_norm, n+n_dt_def, &
      begin, t+t_species, l+l_gph, l+l_gph_precision, l+l_temperature, &
             l+l_temperature_prec, n+n_dt_def, &
      begin, t+t_tgridcoord, l+l_theta, n+n_dt_def, &
      begin, t+t_tgridtype, l+l_logarithmic, n+n_dt_def, &
      begin, t+t_trapezoid, l+l_correct, l+l_wrong, n+n_dt_def, &
      begin, t+t_units, &
             l+l_c, l+l_days, l+l_deg, l+l_degrees, l+l_dimensionless, &
             l+l_dimless, l+l_dl, l+l_dobsonUnits, l+l_DU, l+l_ghz, &
             l+l_hours, l+l_hpa, l+l_hz, l+l_k, l+l_khz, l+l_km, l+l_logp, &
             l+l_m, l+l_maf, l+l_mafs, l+l_mb, l+l_meters, l+l_mhz, &
             l+l_mif, l+l_mifs, l+l_minutes, l+l_molcm2, l+l_orbits, &
             l+l_pa, l+l_ppbv, l+l_ppmv, l+l_pptv, l+l_rad, l+l_radians, &
             l+l_s, l+l_seconds, l+l_thz, l+l_vmr, l+l_zeta, &
             n+n_dt_def, &
      begin, t+t_hgridcoord, l+l_phiTan, l+l_time, n+n_dt_def, &
      begin, t+t_vgridcoord, l+l_angle, l+l_dimensionless, l+l_dimless, &
             l+l_geocAltitude, l+l_geodAltitude, l+l_gph, l+l_icedensity, &
             l+l_integer, l+l_none, l+l_pressure, l+l_theta, l+l_zeta, &
             n+n_dt_def, &
      begin, t+t_vgridtype, l+l_explicit, l+l_linear, l+l_logarithmic, &
             l+l_l2gp, n+n_dt_def /) )

    ! Define the relationships between specs and fields, and the field types
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
    ! For a single field name, an alternative form is
    !               < n_or f_field_name
    !                 < n_field_type t_type ... t_type > ...
    !                 < n_field_spec s_spec ... s_spec > ...
    !                 < n_dot s_spec f_field_name ... >
    !               >
    ! The n_or tree allows a field to take several forms.  It is the
    ! responsibility of the examiner of the tree to determine which
    ! form appears.

    call make_tree ( (/ &
      begin, s+s_time, &
             begin, f+f_reset, boolean(), &
             np+n_spec_def, &
      begin, s+s_gridded, &
             begin, f+f_date, string(), &
             begin, f+f_dimList, string(), &
             begin, f+f_downsample, boolean(), &
             begin, f+f_file, string(), &
             begin, f+f_missingValue, numeric(), &
             begin, f+f_field, string(), &
             begin, f+f_noPCFid, boolean(), &
             begin, f+f_origin, field_type(t_griddedOrigin), &
             begin, f+f_sum, boolean(), &
             np+n_spec_def, &
!      begin, s+s_l2gp, np+n_spec_def, & ! To avoid forward reference in h/vGrid
      begin, s+s_hGrid, &
             begin, f+f_coordinate, field_type(t_hgridcoord), & ! PHYQ_Time, PHYQ_PhiTan
             begin, f+f_date, string(), &
             begin, f+f_extendible, boolean(), &
             begin, f+f_fraction, numeric(phyq_dimensionless), &
             begin, f+f_forbidOverspill, boolean(), &
             begin, f+f_geodAngle, numeric(phyq_angle), &
             begin, f+f_geodLat, numeric(phyq_angle), &
             begin, f+f_height, numeric(phyq_length), &
             begin, f+f_insetOverlaps, boolean(), &
             begin, f+f_interpolationfactor, numeric(phyq_dimensionless), &
             begin, f+f_inclination, numeric(phyq_angle), &
             begin, f+f_lon, numeric(phyq_angle), &
             begin, f+f_losAngle, numeric(phyq_angle), &
             begin, f+f_maxLowerOverlap, numeric(phyq_dimensionless), &
             begin, f+f_maxUpperOverlap, numeric(phyq_dimensionless), &
             begin, f+f_mif, numeric(phyq_dimensionless), &
             begin, f+f_module, field_spec(s_module), &
             begin, f+f_origin, numeric(phyq_angle), &
             begin, f+f_QTMlevel, &
                    begin, numeric(phyq_dimensionless), &
                    begin, numeric(phyq_angle), &
                    begin, numeric(phyq_length), &
                    n+n_or, &
             begin, f+f_single, boolean(), &
             begin, f+f_solarTime, numeric(phyq_time), &
             begin, f+f_solarZenith, &
                    begin, numeric(phyq_dimensionless), &
                    begin, numeric(phyq_angle), &
                    n+n_or, &
             begin, f+f_spacing, numeric(phyq_angle), &
             begin, f+f_sourceL2GP, field_spec(s_l2gp), &
             begin, f+f_Time, &
                    begin, numeric(phyq_dimensionless), &
                    begin, numeric(phyq_time), &
                    n+n_or, &
             begin, f+f_type, field_type(t_hgridtype,req=req), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_readgriddedData, &
             begin, f+f_date, string(), &
             begin, f+f_deferReading, boolean(), &
             begin, f+f_dimList, string(), &
             begin, f+f_downsample, boolean(), &
             begin, f+f_file, string(), &
             begin, f+f_grid, field_spec(s_gridded), &
             begin, f+f_missingValue, numeric(), &
             begin, f+f_noPCFid, boolean(), &
             begin, f+f_field, string(), &
             begin, f+f_origin, field_type(t_griddedOrigin), &
             begin, f+f_sum, boolean(), &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_binSelector, &
             begin, f+f_type, field_type(t_binSelectorType,req=req), &
             begin, f+f_molecule, field_type(t_molecule), &
             begin, f+f_nameFragment, string(), &
             begin, f+f_height, numeric_range(), &
             begin, f+f_cost, numeric(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_empiricalGeometry, &
             begin, f+f_terms, numeric(phyq_dimensionless,req=req), &
             begin, f+f_iterations, numeric(phyq_dimensionless), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_ConvertEtaToP, &  ! Must be AFTER S_Gridded, V_Gridded
             begin, f+f_a, field_spec(s_gridded,s_concatenate), &
             begin, f+f_b, field_spec(s_gridded,s_concatenate), &
             begin, f+f_Grid, field_spec(s_Gridded), &
             begin, f+f_vGrid, field_spec(s_vGrid,scalar=scalar), &
             begin, f+f_logBasis, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_concatenate, &  ! Must be AFTER S_Gridded and S_ConvertEtaToP
             begin, f+f_a, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_b, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_sourceGrid, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_allowEmptyGrids, boolean(), &
             begin, f+f_deleteGrids, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_concatenateGrids, &  ! Must be AFTER S_Gridded, S_Concatenate,
                                      ! and S_ConvertEtaToP
             begin, f+f_grid, field_spec(s_gridded,s_concatenate, &
                               s_convertetatop), &
             begin, f+f_a, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_b, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_sourceGrid, field_spec(s_gridded,s_convertetatop), &
             begin, f+f_deleteGrids, boolean(), &
             begin, f+f_allowEmptyGrids, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_wmoTrop, &  ! Must be AFTER S_Gridded and S_Concatenate
             begin, f+f_a, field_spec(s_gridded,s_concatenate), &
             begin, f+f_b, field_spec(s_gridded,s_concatenate), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_wmoTropFromGrids, &  ! Must be AFTER S_Gridded and S_Concatenate
             begin, f+f_grid, field_spec(s_gridded,s_concatenate), &
             begin, f+f_a, field_spec(s_gridded,s_concatenate), &
             begin, f+f_b, field_spec(s_gridded,s_concatenate), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &  ! Must be AFTER S_Gridded, S_Concatenate, S_Merge,
                           ! and S_ConvertEtaToP
      begin, s+s_delete, &
             begin, f+f_grid, field_spec( s_gridded, s_concatenate, s_merge, &
                                       &  s_ConvertEtaToP ), &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_chunkDivide, &
             begin, f+f_module, field_spec(s_module), &
             begin, f+f_method, field_type(t_chunkDivideMethod,req=req), &
             begin, f+f_noChunks, numeric(phyq_dimensionless), &
             begin, f+f_overlap, numeric(), & ! PHYQ_MAFs, PHYQ_Angle, PHYQ_Time
             begin, f+f_lowerOverlap, numeric(), & ! PHYQ_MAFs, PHYQ_Angle, PHYQ_Time
             begin, f+f_upperOverlap, numeric(), & ! PHYQ_MAFs, PHYQ_Angle, PHYQ_Time
             begin, f+f_maxLength, numeric(), & ! PHYQ_MAFs, PHYQ_Angle, PHYQ_Time
             begin, f+f_maxOrbY, numeric(phyq_length), &
             begin, f+f_excludePostOverlaps, boolean(), &
             begin, f+f_excludePriorOverlaps, boolean(), &
             begin, f+f_noSlaves, numeric(phyq_dimensionless), &
             begin, f+f_homeModule, field_type(t_module), &
             begin, f+f_homeGeodAngle, numeric(phyq_angle), &
             begin, f+f_scanLowerLimit, numeric_range(), &
             begin, f+f_scanUpperLimit, numeric_range(), &
             begin, f+f_criticalBands, string(), &
             begin, f+f_criticalModules, field_type(t_criticalModule), &
             begin, f+f_criticalSignals, string(), &
             begin, f+f_maxGap, numeric(), & ! PHYQ_MAFs, PHYQ_Angle, PHYQ_Time
             begin, f+f_saveObstructions, boolean(), &
             begin, f+f_skipL1BCheck, boolean(), &
             begin, f+f_crashIfPhiNotMono, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_Boolean, &
             begin, f+f_formula, string(), &
             begin, f+f_values, string(), &
             begin, f+f_label, string(), &
             begin, f+f_evaluate, boolean(), &
             begin, f+f_literal, boolean(), &
             ndp+n_spec_def, &
      begin, s+s_fGrid, &
             begin, f+f_coordinate, field_type(t_fGridCoord), &
             begin, f+f_values, numeric(), & ! phyq_frequency, phyq_dimensionless
             nadp+n_spec_def, &
      begin, s+s_tGrid, &
             begin, f+f_formula, numeric_range(), &
             begin, f+f_number, numeric(phyq_dimensionless), &
             begin, f+f_start, numeric(phyq_temperature,req=req), &
             begin, f+f_stop, numeric(), &
             begin, f+f_type, field_type(t_vGridType,req=req), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_vGrid, &
             begin, f+f_type, field_type(t_vGridType,req=req), &
             begin, f+f_coordinate, field_type(t_vGridCoord), &
             begin, f+f_formula, numeric_range(), &
             begin, f+f_number, numeric(phyq_dimensionless), &
             begin, f+f_resolution, numeric(phyq_dimensionless), &
             begin, f+f_start, numeric(), &
             begin, f+f_stop, numeric(), &
             begin, f+f_values, numeric(), &
             begin, f+f_sourceL2GP, field_spec(s_l2gp), &
             ndp+n_spec_def, &
      begin, s+s_forge, &  ! Must be AFTER S_Module and S_VGrid
             begin, f+f_module, field_spec(s_module,req=req), &
             begin, f+f_solarTime, numeric(), &
             begin, f+f_solarZenith, numeric(), &
             begin, f+f_geodAngle, numeric(), &
             begin, f+f_geodAlt, field_spec(s_vGrid), &
             begin, f+f_noMIFs, numeric(req=req), &
             begin, f+f_truncate, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_merge, &  ! Must be AFTER S_Gridded and S_VGrid
             begin, f+f_operational, field_spec(s_gridded,s_concatenate, &
                    s_convertetatop), &
             begin, f+f_climatology, field_spec(s_gridded,s_concatenate), &
             begin, f+f_height, numeric(phyq_pressure), &
             begin, f+f_scale, numeric(phyq_length), &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_mergeGrids, &  ! Must be AFTER S_Gridded and S_VGrid
             begin, f+f_grid, field_spec(s_gridded,s_concatenate, &
                               s_convertetatop), &
             begin, f+f_operational, field_spec(s_gridded,s_concatenate, &
                               s_convertetatop), &
             begin, f+f_climatology, field_spec(s_gridded,s_concatenate), &
             begin, f+f_height, numeric(), &
             begin, f+f_scale, numeric(), &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_makePFA, & ! Must be AFTER s_vGrid and s_tGrid
             begin, f+f_allLinesForRadiometer, boolean(), &
             begin, f+f_allLinesInCatalog, boolean(), &
             begin, f+f_losvel, numeric(phyq_velocity,req=req,scalar=scalar), &
             begin, f+f_molecules, field_type(t_molecule,req=req), &
             begin, f+f_oversample, numeric(phyq_dimensionless,scalar=scalar), &
             begin, f+f_signals, string(req), &
             begin, f+f_temperatures, field_spec(s_tGrid,req=req,scalar=scalar), &
             begin, f+f_vGrid, field_spec(s_vGrid,req=req,scalar=scalar), &
             ndp+n_spec_def, &
      begin, s+s_pfaData, & ! Must be AFTER s_vGrid and s_tGrid
             begin, f+f_absorption, numeric(phyq_dimensionless), &
             begin, f+f_dAbsDnc, numeric(phyq_dimensionless), &
             begin, f+f_dAbsDnu, numeric(phyq_dimensionless), &
             begin, f+f_dAbsDwc, numeric(phyq_dimensionless), &
             begin, f+f_molecules, field_type(t_molecule), &
             begin, f+f_signal, string(req), &
             begin, f+f_temperatures, field_spec(s_tGrid,req=req), &
             begin, f+f_velLin, numeric(phyq_velocity), &
             begin, f+f_vGrid, field_spec(s_vGrid,req=req), &
             ndp+n_spec_def, &
      begin, s+s_readPFA, &
             begin, f+f_file, string(req), &
             begin, f+f_molecules, field_type(t_molecule), &
             begin, f+f_signals, string(), &
             ndp+n_spec_def, &
      begin, s+s_writePFA, &
             begin, f+f_allPFA, boolean(), &
             begin, f+f_file, string(req), &
             begin, f+f_pfaData, field_spec(s_pfaData,s_makePFA), &
             ndp+n_spec_def, &
      begin, s+s_flushPFA, &
             begin, f+f_molecules, field_type(t_molecule), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_changeSettings, & ! Ignores rest of stuff (Why do you say that?)
             begin, f+f_additional, boolean(), &
             begin, f+f_debug, boolean(), &
             begin, f+f_options, string(), &
             begin, f+f_reset, boolean(), &
             begin, f+f_silent, boolean(), &
             begin, f+f_skipDirectWrites, boolean(), &
             begin, f+f_skipDirectWritesif, field_spec(s_Boolean), &
             begin, f+f_skipRetrieval, boolean(), &
             begin, f+f_skipRetrievalif, field_spec(s_Boolean), &
             begin, f+f_stamp, boolean(), &
             begin, f+f_verbose, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_phase, & ! Ignores rest of stuff (Why do you say that?)
             begin, f+f_additional, boolean(), &
             begin, f+f_debug, boolean(), &
             begin, f+f_options, string(), &
             begin, f+f_reset, boolean(), &
             begin, f+f_silent, boolean(), &
             begin, f+f_skipDirectWrites, boolean(), &
             begin, f+f_skipDirectWritesif, field_spec(s_Boolean), &
             begin, f+f_skipRetrieval, boolean(), &
             begin, f+f_skipRetrievalif, field_spec(s_Boolean), &
             begin, f+f_stamp, boolean(), &
             begin, f+f_verbose, boolean(), &
             ndp+n_spec_def, &
      begin, s+s_quantity, & ! Must be AFTER [F, H, IWC, T, V]grid
             begin, f+f_badValue, numeric(), &
             begin, f+f_coordinate, field_type(t_vGridCoord), &
             begin, f+f_fGrid, field_spec(s_fgrid), &
             begin, f+f_hGrid, field_spec(s_hgrid), &
             begin, f+f_irregular, boolean(), &
             begin, f+f_keepChannels, boolean(), &
             begin, f+f_logBasis, boolean(), &
             begin, f+f_minValue, numeric(), &
             begin, f+f_module, field_spec(s_module), &
             begin, f+f_molecule, field_type(t_molecule), &
             begin, f+f_radiometer, field_spec(s_radiometer), &
             begin, f+f_reflector, field_type(t_reflector), &
             begin, f+f_sGrid, field_spec(s_vgrid), &
             begin, f+f_signal, string(), &
             begin, f+f_stacked, boolean(), &
             begin, f+f_type, field_type(t_quantityType), &
             begin, f+f_unit, field_type(t_units), &
             begin, f+f_vGrid, field_spec(s_vgrid), &
             begin, f+f_xGrid, field_spec(s_hgrid), &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_vectorTemplate, & ! Must be AFTER s_quantity
             begin, f+f_adopt, string(), &
             begin, f+f_quantities, field_spec(s_quantity), &
             begin, f+f_removeQuantities, field_spec(s_quantity), &
             begin, f+f_removeTemplate, field_spec(s_vectorTemplate), &
             begin, f+f_source, field_type(t_rowsOrColumns), &
             begin, f+f_template, field_spec(s_vectorTemplate), &
             ndp+n_spec_def, &
      begin, s+s_vector, & ! Must be AFTER s_vectorTemplate
             begin, f+f_template, field_spec(s_vectorTemplate,req=req), &
             begin, f+f_autoFill, boolean(), &
             begin, f+f_fraction, boolean(), &
             begin, f+f_highBound, boolean(), &
             begin, f+f_lengthScale, boolean(), &
             begin, f+f_lowBound, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_hessian, &   ! Must be AFTER s_vector
             begin, f+f_rows, field_spec(s_vector,req=req), &
             begin, f+f_columns, field_spec(s_vector,req=req), &
             ndp+n_spec_def, &
      begin, s+s_l2gp, &   ! Must be AFTER s_vector
             begin, f+f_source, vectorQuantity(), &
             begin, f+f_file, string(), &
             begin, f+f_compareOverlaps, boolean(), &
             begin, f+f_outputOverlaps, boolean(), &
             begin, f+f_precision, vectorQuantity(), &
             begin, f+f_prefixSignal, boolean(), &
             begin, f+f_swath, string(), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless), &
             begin, f+f_AuraInstrument, string(), &
             ndp+n_spec_def, &
      begin, s+s_l2aux, &   ! Must be AFTER s_vector
             begin, f+f_source, vectorQuantity(), &
             begin, f+f_compareOverlaps, boolean(), &
             begin, f+f_outputOverlaps, boolean(), &
             begin, f+f_prefixSignal, boolean(), &
             begin, f+f_file, string(), &
             begin, f+f_quantityType, field_type(t_quantityType), &
             begin, f+f_sdname, string(), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless), &
             ndp+n_spec_def, &
      begin, s+s_matrix, &  ! Must be AFTER s_vector
             begin, f+f_rows, field_spec(s_vector), &
             begin, f+f_columns, field_spec(s_vector,req=req), &
             begin, f+f_type, field_type(t_matrix), &
             begin, f+f_source, string(), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_fill, &  ! Must be AFTER s_vector,s_matrix and s_climatology
             begin, f+f_a, vectorQuantity(), &
             begin, f+f_avoidBrightObjects, string(), &
             begin, f+f_additional, boolean(), &
             begin, f+f_allowMissing, boolean(), &
             begin, f+f_aprioriPrecision, vectorQuantity(), &
             begin, f+f_asPercentage, boolean(), &
             begin, f+f_b, vectorQuantity(), &
             begin, f+f_badRange, numeric_range(), &
             begin, f+f_baselineQuantity, vectorQuantity(), &
             begin, f+f_block, numeric(phyq_dimensionless), &
             begin, f+f_boundaryPressure, vectorQuantity(), &
             begin, f+f_boxCarMethod, field_type(t_boxCarMethod), &
             begin, f+f_c, numeric(), &
             begin, f+f_channel, numeric(phyq_dimensionless), &
             begin, f+f_channels, numeric_or_range(), &
             begin, f+f_centerVertically, boolean()/) )
    call make_tree ( (/ & ! Continuing for s_fill...
             begin, f+f_count, numeric(phyq_dimensionless), &
             begin, f+f_dimList, string(), &
             begin, f+f_dontLatch, boolean(), &
             begin, f+f_dontMask, boolean(), &
             begin, f+f_earthRadius, vectorQuantity(), &
             begin, f+f_ECRtoFOV, vectorQuantity(), &
             begin, f+f_exact, boolean(), &
             begin, f+f_excludeBelowBottom, boolean(), &
             begin, f+f_explicitValues, numeric(), & ! quantity%template%unit, phyq_dimensionless
             begin, f+f_expr, &
                    begin, numeric(), &
                    begin, boolean(), &
                    begin, vectorQuantity(expr=.true.), &
                    n+n_or, &
             begin, f+f_extinction, boolean(), &
             begin, f+f_fieldECR, vectorQuantity(), &
             begin, f+f_file, string(), &
             begin, f+f_flags, vectorQuantity(), &
             begin, f+f_force, boolean(), &
             begin, f+f_fromPrecision, boolean(), &
             begin, f+f_geocAltitudeQuantity, vectorQuantity(), &
             begin, f+f_geolocation, field_type(t_geolocation), &
             begin, f+f_gphQuantity, vectorQuantity(), &
             begin, f+f_h2oQuantity, vectorQuantity(), &
             begin, f+f_h2oPrecisionQuantity, vectorQuantity(), &
             begin, f+f_height, numeric_or_range(), &
             begin, f+f_heightRange, string()/),  &
             continue=.true. )
    call make_tree ( (/ & ! Continuing for s_fill...
             begin, f+f_hGrid, field_spec(s_hgrid), &
             begin, f+f_ifMissingGMAO, boolean(), &
             begin, f+f_ignoreGeolocation, boolean(), &
             begin, f+f_ignoreNegative, boolean(), &
             begin, f+f_ignoreQuantities, field_spec(s_quantity), &
             begin, f+f_ignoreTemplate, boolean(), &
             begin, f+f_ignoreZero, boolean(), &
             begin, f+f_instances, numeric_or_range(), &
             begin, f+f_internalVGrid, field_spec(s_vGrid), &
             begin, f+f_integrationTime, numeric(phyq_time), &
             begin, f+f_interpolate, boolean(), &
             begin, f+f_intrinsic, boolean(), &
             begin, f+f_isPrecision, boolean(), &
             begin, f+f_logSpace, boolean(), &
             begin, f+f_losQty, vectorQuantity(), &
             begin, f+f_lsb, vectorQuantity(), &
             begin, f+f_lsbFraction, vectorQuantity(), &
             begin, f+f_manipulation, string(), &
             begin, f+f_maxIterations, numeric(phyq_dimensionless), &
             begin, f+f_maxValue, numeric(), & ! sourcequantity%template%unit, phyq_dimensionless
             begin, f+f_minNormQty, vectorQuantity(), &
             begin, f+f_minValue, numeric(), & ! sourcequantity%template%unit, phyq_dimensionless
             begin, f+f_measurements, vectorQuantity() /), &
             continue=.true. )
    call make_tree ( (/ & ! Continuing for s_fill...
             begin, f+f_method, field_type(t_fillmethod,req=req), &
             begin, f+f_model, vectorQuantity(), &
             begin, f+f_multiplier, numeric(), &
             begin, f+f_noFineGrid, numeric(phyq_dimensionless), &
             begin, f+f_noise, vectorQuantity(), &
             begin, f+f_noiseBandwidth, vectorQuantity(), &
             begin, f+f_normQty, vectorQuantity(), &
             begin, f+f_offsetAmount, numeric(phyq_temperature), &
             begin, f+f_orbitInclination, vectorQuantity() /), &
             continue=.true. )
    call make_tree ( (/ & ! Continuing for s_fill...
             begin, f+f_phiWindow, numeric_or_range(), & ! phyq_profiles, phyq_angle
             begin, f+f_phiZero, numeric(phyq_angle), &
             begin, f+f_precision, vectorQuantity(), &
             begin, f+f_precisionFactor, numeric(phyq_dimensionless), &
             begin, f+f_profile, numeric(phyq_dimensionless), &
             begin, f+f_profileValues, numeric_range(), &
             begin, f+f_ptanQuantity, vectorQuantity(), &
             begin, f+f_phitan, vectorQuantity(), &
             begin, f+f_quadrature, boolean(), &
             begin, f+f_quantity, vectorQuantity(req=req), &
             begin, f+f_ratioQuantity, vectorQuantity(), &
             begin, f+f_radianceQuantity, vectorQuantity(), &
             begin, f+f_refract, boolean(), &
             begin, f+f_refGPHQuantity, vectorQuantity(), &
             begin, f+f_refGPHPrecisionQuantity, vectorQuantity(), &
             begin, f+f_referenceMIF, numeric(), & ! dimless or length
             begin, f+f_regular, boolean(), &
             begin, f+f_replaceMissingValue, numeric(), &
             begin, f+f_resetSeed, boolean(), &
             begin, f+f_rhiPrecisionQuantity, vectorQuantity(), &
             begin, f+f_rhiQuantity, vectorQuantity() /), &
             continue = .true. )
    call make_tree ( (/ & ! STILL Continuing for s_fill...
             begin, f+f_scale, numeric(phyq_dimensionless), &
             begin, f+f_scaleInsts, numeric(phyq_dimensionless), &
             begin, f+f_scaleRatio, numeric(phyq_dimensionless), &
             begin, f+f_scaleSurfs, numeric(phyq_dimensionless), &
             begin, f+f_scVel, vectorQuantity(), &
             begin, f+f_scVelECI, vectorQuantity(), &
             begin, f+f_scVelECR, vectorQuantity(), &
             begin, f+f_scECI, vectorQuantity(), &
             begin, f+f_sdname, string(), &
             begin, f+f_seed, numeric(phyq_dimensionless), &
             begin, f+f_shape, numeric(), &
             begin, f+f_sourceGrid, field_spec( s_gridded, s_merge, s_concatenate, &
                                             &  s_ConvertEtaToP, s_wmoTrop ), &
             begin, f+f_sourceL2GP, field_spec(s_l2gp), &
             begin, f+f_sourceL2AUX, field_spec(s_l2aux), &
             begin, f+f_sourceMask, boolean(), &
             begin, f+f_sourceQuantities, field_spec(s_quantity), &
             begin, f+f_sourceQuantity, vectorQuantity(), &
             begin, f+f_sourceType, string(), &
             begin, f+f_sourceVGrid, field_spec(s_vGrid), &
             begin, f+f_spread, boolean(), &
             begin, f+f_start, numeric(phyq_dimensionless), &
             begin, f+f_status, numeric(phyq_dimensionless), &
             begin, f+f_stride, numeric(phyq_dimensionless), &
             begin, f+f_suffix, string() , &
             begin, f+f_surface, numeric_or_range(), &
             begin, f+f_systemTemperature, vectorQuantity() /), &
             continue = .true. )
    call make_tree ( (/ & ! STILL Continuing for s_fill...
             begin, f+f_temperatureQuantity, vectorQuantity(), &
             begin, f+f_tempPrecisionQuantity, vectorQuantity(), &
             begin, f+f_terms, numeric(phyq_dimensionless), &
             begin, f+f_totalPowerQuantity, vectorQuantity(), &
             begin, f+f_tngtECI, vectorQuantity(), &
             begin, f+f_unit, field_type(t_units), &
             begin, f+f_usb, vectorQuantity(), &
             begin, f+f_usbFraction, vectorQuantity(), &
             begin, f+f_vmrQuantity, vectorQuantity(), &
             begin, f+f_whereFill, boolean(), &
             begin, f+f_whereNotFill, boolean(), &
             begin, f+f_width, numeric(phyq_dimensionless), &
             ndp+n_spec_def /), &
             continue = .true. ) ! WHEW! Finally done for s_fill

    call make_tree( (/ &
      begin, s+s_flushL2PCBins, ndp+n_spec_def /) )

    call make_tree ( (/ &
      begin, s+s_populateL2PCBin, &
             begin, f+f_bin, string(), &
             nadp+n_spec_def /) )

    call make_tree ( (/ &
      begin, s+s_leakCheck, &
             begin, f+f_where, string(), &
             ndp+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_load, &
             begin, f+f_bin, string(req), &
             begin, f+f_hessian, field_spec(s_hessian), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_vector, field_spec(s_vector), &
             begin, f+f_source, field_type(t_rowsOrColumns), &
             ndp+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_reevaluate, &
             begin, f+f_additional, boolean(), &
             begin, f+f_a, vectorQuantity(), &
             begin, f+f_b, vectorQuantity(), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_c, numeric(), &
             begin, f+f_evaluate, boolean(), &
             begin, f+f_expr, &
                    begin, numeric(), &
                    begin, boolean(), &
                    begin, vectorQuantity(), &
                    n+n_or, &
             begin, f+f_formula, string(), &
             begin, f+f_inputBoolean, field_spec(s_Boolean), &
             begin, f+f_truncate, boolean(), &
             begin, f+f_label, string(), &
             begin, f+f_literal, boolean(), &
             begin, f+f_manipulation, string(), &
             begin, f+f_values, string(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_catchwarning, &
             begin, f+f_message, string(), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_compare, &
             begin, f+f_a, vectorQuantity(), &
             begin, f+f_b, vectorQuantity(), &
             begin, f+f_c, numeric(), &
             begin, f+f_formula, string(req), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_anyGoodRadiances, &
             begin, f+f_signal, string(req), &
             begin, f+f_reverse, boolean(), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_anyGoodValues, &
             begin, f+f_quantity, vectorQuantity(), &
             begin, f+f_precision, vectorQuantity(), &
             begin, f+f_quality, vectorQuantity(), &
             begin, f+f_reverse, boolean(), &
             begin, f+f_status, vectorQuantity(), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_isFileAbsent, &
             begin, f+f_file, string(), &
             begin, f+f_type, field_type(t_griddedOrigin), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_reverse, boolean(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_ChunkEndsBefore, &
             begin, f+f_grid, field_spec(s_Gridded,s_concatenate,req=req), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_reverse, boolean(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_ChunkStartsAfter, &
             begin, f+f_grid, field_spec(s_Gridded,s_concatenate,req=req), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_reverse, boolean(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_isGridEmpty, &
             begin, f+f_grid, field_spec(s_Gridded,s_concatenate,req=req), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_reverse, boolean(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_isSwathEmpty, &
             begin, f+f_file, string(), &
             begin, f+f_swath, string(), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             begin, f+f_noPCFid, boolean(), &
             begin, f+f_Boolean, field_spec(s_Boolean,req=req), &
             begin, f+f_reverse, boolean(), &
             np+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_writeFileAttribute, &
             begin, f+f_additional, boolean(), &
             begin, f+f_attrName, string(req), &
             begin, f+f_file, string(req), &
             begin, f+f_attrValue, string(req), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             ndp+n_spec_def /) )
             
    call make_tree ( (/ &
      begin, s+s_subSample, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_file, string(), &
             begin, f+f_hGrid, field_spec(s_hgrid), &
             begin, f+f_type, field_type(t_outputType), &
             begin, f+f_repairGeolocations, boolean(), &
             ndp+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_label, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_quantity, vectorQuantity(), &
             begin, f+f_vector, field_spec(s_vector), &
             begin, f+f_label, string(), &
             begin, f+f_prefixSignal, boolean(), &
             begin, f+f_suffixLabel, boolean(), &
             ndp+n_spec_def /) )

    call make_tree( (/ &
      begin, s+s_destroy, &
             begin, f+f_allGriddedData, boolean(), &
             begin, f+f_allMatrices, boolean(), &
             begin, f+f_allVectors, boolean(), &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_grid, field_spec( s_gridded, s_concatenate, s_merge, &
                                       &  s_ConvertEtaToP ), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_vector, field_spec(s_vector), &
             np+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_transfer, &
             begin, f+f_a, field_spec(s_vector), &
             begin, f+f_b, field_spec(s_vector), &
             begin, f+f_c, numeric(), &
             begin, f+f_method, field_type(t_fillmethod), &
             begin, f+f_source, field_spec(s_vector), &
             begin, f+f_destination, field_spec(s_vector,req=req), &
             begin, f+f_dontMask, boolean(), &
             begin, f+f_ignoreNegative, boolean(), &
             begin, f+f_ignoreZero, boolean(), &
             begin, f+f_interpolate, boolean(), &
             begin, f+f_manipulation, string(), &
             begin, f+f_measurements, field_spec(s_vector), &
             begin, f+f_model, field_spec(s_vector), &
             begin, f+f_noise, field_spec(s_vector), &
             begin, f+f_ptanQuantity, vectorQuantity(), &
             begin, f+f_quantityNames, field_spec(s_Boolean), &
             begin, f+f_expandMask, boolean(), &
             begin, f+f_skipMask, boolean(), &
             begin, f+f_skipValues, boolean(), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_negativePrecision, &
             begin, f+f_precision, field_spec(s_vector,req=req), &
             begin, f+f_aprioriPrecision, field_spec(s_vector,req=req), &
             begin, f+f_precisionFactor, numeric(phyq_dimensionless), &
             begin, f+f_noZeros, boolean(), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_fillCovariance, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_diagonal, field_spec(s_vector,req=req), &
             begin, f+f_ignoreTemplate, boolean(), &
             begin, f+f_lengthScale, field_spec(s_vector), &
             begin, f+f_fraction, field_spec(s_vector), &
             begin, f+f_invert, boolean(), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_fillDiagonal, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_diagonal, field_spec(s_vector,req=req), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_reflect, & ! Must be AFTER s_matrix
             begin, f+f_matrix, field_spec(s_matrix), &
             nadp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_cyclicJacobi, & ! Must be AFTER s_matrix
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_eigenVectors, field_spec(s_matrix), &
             nadp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_columnScale, & ! Must be AFTER s_matrix
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_scale, field_spec(s_vector), &
             nadp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_rowScale, & ! Must be AFTER s_matrix
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_scale, field_spec(s_vector), &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_directRead, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_bin, string(), &
             begin, f+f_file, string(), &
             begin, f+f_geolocation, boolean(), &
             begin, f+f_groupName, string(), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless,req=req), &
             begin, f+f_interpolate, boolean(), &
             begin, f+f_noPCFid, boolean(), &
             begin, f+f_options, string(), &
             begin, f+f_quantity, vectorQuantity(), &
             begin, f+f_rank, numeric(phyq_dimensionless), &
             begin, f+f_sdname, string(), &
             begin, f+f_spread, boolean(), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             begin, f+f_vector, field_spec(s_vector), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_hessian, field_spec(s_hessian), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_directWrite, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_convergence, vectorQuantity(), &
             begin, f+f_file, string(), &
             begin, f+f_groupName, string(), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless,req=req), &
             begin, f+f_inputFile, string(), &
             begin, f+f_lowerOverlap, boolean(), &
             begin, f+f_noPCFid, boolean(), &
             begin, f+f_label, string(), &
             begin, f+f_options, string(), &
             begin, f+f_AscDescMode, vectorQuantity(), &
             begin, f+f_geoHeight, vectorQuantity(), &
             begin, f+f_precision, vectorQuantity(), &
             begin, f+f_quality, vectorQuantity(), &
             begin, f+f_rank, numeric(phyq_dimensionless), &
             begin, f+f_single, boolean(), &
             begin, f+f_source, vectorQuantity(), &
             begin, f+f_status, vectorQuantity(), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             begin, f+f_upperOverlap, boolean(), &
             begin, f+f_vector, field_spec(s_vector), &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_copy, &  ! Must be AFTER s_hGrid if repairGeoLocations
             begin, f+f_create, boolean(), &
             begin, f+f_exclude, string(), &
             begin, f+f_file, string(req), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless), &
             begin, f+f_hGrid, field_spec(s_hgrid), &
             begin, f+f_ifAnyCrashedChunks, boolean(), &
             begin, f+f_inputFile, string(req), &
             begin, f+f_inputtype, field_type(t_outputType), &
             begin, f+f_options, string(), &
             begin, f+f_rename, string(), &
             begin, f+f_repairGeolocations, boolean(), &
             begin, f+f_swath, string(), &
             begin, f+f_toAttribute, boolean(), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_catenate, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_file, string(), &
             begin, f+f_hGrid, field_spec(s_hgrid), &
             begin, f+f_type, field_type(t_outputType), &
             begin, f+f_repairGeolocations, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_sleep, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_time, numeric(phyq_time), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_output, &  ! Must be AFTER s_l2aux and s_l2gp
             begin, f+f_ascii, boolean(), &
             begin, f+f_destroy, boolean(), &
             begin, f+f_dontPack, field_spec(s_quantity), &
             begin, f+f_file, string(req), &
             begin, f+f_hdfVersion, numeric(phyq_dimensionless), &
             begin, f+f_metaDataOnly, boolean(), &
             begin, f+f_metaName, string(), &
             begin, f+f_moleculeSecondDerivatives, field_type(t_molecule), &
             begin, f+f_overlaps, field_spec(s_l2aux,s_l2gp), &
             begin, f+f_packed, boolean(), &
             begin, f+f_quantities, field_spec( s_l2aux,s_l2gp,s_matrix, &
                                             &  s_hessian ), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             begin, f+f_writeCounterMAF, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_subset, &  ! Must be AFTER s_vector
             begin, f+f_quantity, vectorQuantity(req), &
             begin, f+f_a, vectorQuantity(), &
             begin, f+f_ptanquantity, vectorQuantity(), &
             begin, f+f_radiancequantity, vectorQuantity(), &
             begin, f+f_sourcequantity, vectorQuantity(), &
             begin, f+f_additional, boolean(), &
             begin, f+f_channels, numeric_or_range(), &
             begin, f+f_height, numeric_or_range(), &
             begin, f+f_heightRange, string(), &
             begin, f+f_ignore, boolean(), &
             begin, f+f_instances, numeric_or_range(), &
             begin, f+f_mask, field_type(t_masks), &
             begin, f+f_maxValue, numeric(), &
             begin, f+f_minValue, numeric(), &
             begin, f+f_missingValue, numeric(phyq_dimensionless), &
             begin, f+f_surface, numeric_or_range(), &
             begin, f+f_opticalDepth, vectorQuantity(), &
             begin, f+f_opticalDepthCutoff, numeric(phyq_dimensionless), &
             begin, f+f_reverse, boolean(), &
             begin, f+f_where, string(), &
             begin, f+f_resetAll, boolean(), &
             begin, f+f_reset, boolean(), ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_restrictRange, &
             begin, f+f_quantity, vectorQuantity(req), &
             begin, f+f_ptanquantity, vectorQuantity(req), &
             begin, f+f_mask, field_type(t_masks), &
             begin, f+f_measurements, field_spec(s_vector,req=req), &
             begin, f+f_signals, string(req), &
             begin, f+f_basisFraction, numeric(phyq_dimensionless), &
             begin, f+f_minChannels, numeric(phyq_dimensionless), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_flagcloud, &  ! Must be AFTER s_vector
             begin, f+f_quantity, vectorQuantity(req), &
             begin, f+f_cloudRadiance, vectorQuantity(req), &
             begin, f+f_ptanquantity, vectorQuantity(req), &
             begin, f+f_mask, field_type(t_masks), &
             begin, f+f_channels, numeric_or_range(), &
             begin, f+f_cloudchannels, numeric(req=req), &
             begin, f+f_height, numeric_range(req=req), &
             begin, f+f_cloudHeight, numeric_range(), &
             begin, f+f_cloudRadianceCutoff, numeric(phyq_temperature,req=req), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_updateMask, &
             begin, f+f_quantity, vectorQuantity(req), &
             begin, f+f_sourceQuantity, vectorQuantity(), &
             begin, f+f_operation, field_type(t_maskUpdates,req=req), &
             begin, f+f_mask, field_type(t_masks,req=req), &
             begin, f+f_sourceMask, field_type(t_masks), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_forwardModel, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_allLinesForRadiometer, boolean(), &
             begin, f+f_allLinesInCatalog, boolean(), &
             begin, f+f_atmos_der, boolean(), &
             begin, f+f_atmos_second_der, boolean(), &
             begin, f+f_binSelectors, field_spec(s_binSelector), &
             begin, f+f_cloud_der, field_type(t_cloud_der), &
             begin, f+f_default_spectroscopy, boolean(), &
             begin, f+f_differentialScan, boolean(), &
             begin, f+f_do_baseline, boolean(), &
             begin, f+f_do_conv, boolean(), &
             begin, f+f_do_freq_avg, boolean(), &
             begin, f+f_do_1d, boolean(), &
             begin, f+f_forceSidebandFraction, boolean(), &
             begin, f+f_frqTol, numeric(phyq_frequency), &
             begin, f+f_i_saturation, field_type(t_i_saturation),&
             begin, f+f_ignoreHessian, boolean(), &
             begin, f+f_incl_cld, boolean(), &
             begin, f+f_integrationGrid, field_spec(s_vGrid), &
             begin, f+f_linearSideband, numeric(phyq_dimensionless), &
             begin, f+f_lineCenter, field_type(t_molecule), &
             begin, f+f_lineWidth, field_type(t_molecule), &
             begin, f+f_lineWidth_TDep, field_type(t_molecule), &
             begin, f+f_lockBins, boolean(), &
             begin, f+f_lsbLBLMolecules, field_type(t_molecule,empty=empty), &
             begin, f+f_lsbPFAMolecules, field_type(t_molecule,empty=empty), &
             begin, f+f_MIFTangent, field_type(t_MIFTangent), &
             begin, f+f_module, field_spec(s_module), &
             begin, f+f_moleculeDerivatives, field_type(t_molecule), &
             begin, f+f_moleculeSecondDerivatives, field_type(t_molecule), &
             begin, f+f_molecules, field_type(t_molecule,array=.true.), &
             begin, f+f_nabterms, numeric(phyq_dimensionless), &
             begin, f+f_nazimuthangles, numeric(phyq_dimensionless), &
             begin, f+f_ncloudspecies, numeric(phyq_dimensionless), &
             begin, f+f_nmodelsurfs, numeric(phyq_dimensionless), &
             begin, f+f_no_dup_mol, boolean(), &
             begin, f+f_noMagneticField, boolean(), &
             begin, f+f_nscatteringangles, numeric(phyq_dimensionless), &
             begin, f+f_ncloudspecies, numeric(phyq_dimensionless), &
             begin, f+f_nsizebins, numeric(phyq_dimensionless) /) )
    call make_tree ( (/ &
             begin, f+f_pathNorm, boolean(), &
             begin, f+f_phiWindow, numeric_or_range(), & ! [phyq_angle, phyq_profiles]
             begin, f+f_polarized, boolean(), &
             begin, f+f_referenceMIF, numeric(phyq_dimensionless), &
             begin, f+f_refract, boolean(), &
             begin, f+f_scanAverage, boolean(), &
             begin, f+f_signals, string(), &
             begin, f+f_skipOverlaps, boolean(), &
             begin, f+f_switchingMirror, boolean(), &
             begin, f+f_specificQuantities, field_spec(s_quantity), &
             begin, f+f_spect_der, boolean(), &
             begin, f+f_tangentGrid, field_spec(s_vGrid), &
             begin, f+f_temp_der, boolean(), &
             begin, f+f_tolerance, numeric(phyq_temperature), &
             begin, f+f_transformMIFextinction, boolean(), &
             begin, f+f_transformMIFRHI, boolean(), &
             begin, f+f_trapezoid, field_type(t_trapezoid), &
             begin, f+f_TScatMIF, numeric(phyq_dimensionless), &
             begin, f+f_TScatMoleculeDerivatives, field_type(t_molecule), &
             begin, f+f_TScatMolecules, field_type(t_molecule), &
             begin, f+f_type, field_type(t_fwmType,req=req), &
             begin, f+f_usbLBLMolecules, field_type(t_molecule,empty=empty), &
             begin, f+f_usbPFAMolecules, field_type(t_molecule,empty=empty), &
             begin, f+f_useTScat, boolean(), &
             begin, f+f_xStar, field_spec(s_vector), &
             begin, f+f_yStar, field_spec(s_vector), &
             ndp+n_spec_def /), continue=.true. )

    call make_tree ( (/ & ! Must be AFTER s_vector
      begin, s+s_checkpoint, &
             begin, f+f_fileName, string(req=req), &
             begin, f+f_vectors, field_spec(s_vector,req=req), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_computetotalpower, &
             begin, f+f_measurements, field_spec(s_vector,req=req), &
             begin, f+f_weightsVector, field_spec(s_vector,req=req), &
             begin, f+f_totalPowerVector, field_spec(s_vector,req=req), &
             nadp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_endSelect, &
             begin, f+f_label, string(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_select, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_label, string(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_case, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_label, string(), &
             begin, f+f_options, string(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ & ! Must be AFTER s_Boolean
      begin, s+s_repeat, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_formula, string(), &
             begin, f+f_values, string(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ & ! Must be AFTER s_Boolean
      begin, s+s_skip, &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_formula, string(), &
             begin, f+f_nextChunk, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ & ! Must be AFTER s_vector,s_vectorTemplate, etc.
      begin, s+s_diff, &
             begin, f+f_Clean, boolean(), &
             begin, f+f_crashBurn, boolean(), &
             begin, f+f_details, numeric(phyq_dimensionless), &
             begin, f+f_hessian, field_spec(s_hessian), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_options, string(), &
             begin, f+f_quantity, vectorQuantity(), &
             begin, f+f_stop, boolean(), &
             begin, f+f_stopWithError, boolean(), &
             begin, f+f_text, string(), &
             begin, f+f_Grid, field_spec( s_gridded, s_merge, s_Concatenate, &
                                       &  s_ConvertEtaToP ), &
             np+n_spec_def /) )
    call make_tree ( (/ & ! Must be AFTER s_forwardModel,s_hGrid,s_pfaData,
        begin, s+s_dump, & ! s_makePFA,s_vector,s_vectorTemplate, etc.
             begin, f+f_allBooleans, boolean(), &
             begin, f+f_allFiles, boolean(), &
             begin, f+f_allForwardModels, boolean(), &
             begin, f+f_allGriddedData, boolean(), &
             begin, f+f_allHessians, boolean(), &
             begin, f+f_allHGrids, boolean(), &
             begin, f+f_allL2PCs, boolean(), &
             begin, f+f_allLines, boolean(), &
             begin, f+f_allMatrices, boolean(), &
             begin, f+f_allPFA, boolean(), &
             begin, f+f_allQuantityTemplates, boolean(), &
             begin, f+f_allRadiometers, boolean(), &
             begin, f+f_allSignals, boolean(), &
             begin, f+f_allSpectra, boolean(), &
             begin, f+f_allVectors, boolean(), &
             begin, f+f_allVectorTemplates, boolean(), &
             begin, f+f_allVGrids, boolean() /) )
    call make_tree ( (/ & ! Continuing for s_dump...
             begin, f+f_antennaPatterns, boolean(), &
             begin, f+f_block, numeric(), &
             begin, f+f_Boolean, field_spec(s_Boolean), &
             begin, f+f_callStack, boolean(), &
             begin, f+f_commandLine, boolean(), &
             begin, f+f_chunkDivide, boolean(), &
             begin, f+f_chunkNumber, boolean(), &
             begin, f+f_Clean, boolean(), &
             begin, f+f_count, numeric(), &
             begin, f+f_crashBurn, boolean(), &
             begin, f+f_DACSfilterShapes, boolean(), &
             begin, f+f_details, numeric(phyq_dimensionless), &
             begin, f+f_dumpFile, string(), &
             begin, f+f_file, string(), &
             begin, f+f_filterShapes, boolean(), &
             begin, f+f_forwardModel, field_spec(s_forwardModel), &
             begin, f+f_globalAttributes, boolean(), &
             begin, f+f_Grid, field_spec( s_gridded, s_merge, s_Concatenate, &
                                       &  s_ConvertEtaToP ), &
             begin, f+f_hGrid, field_spec(s_hgrid)/) , &
             continue=.true. )
    call make_tree ( (/ & ! Continuing for s_dump...
             begin, f+f_hessian, field_spec(s_hessian), &
             begin, f+f_igrf, boolean(), &
             begin, f+f_inputfile, string(), &
             begin, f+f_l2pc, string(), &
             begin, f+f_lines, field_spec(s_line), &
             begin, f+f_mark, boolean(), &
             begin, f+f_mask, vectorQuantity(), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_memory, boolean(), &
             begin, f+f_MieTables, boolean(), &
             begin, f+f_options, string(), &
             begin, f+f_pause, boolean(), &
             begin, f+f_pfaData, field_spec(s_makePFA,s_pfaData,s_readPFA), &
             begin, f+f_pfaFiles, boolean(), &
             begin, f+f_pfaNum, numeric(phyq_dimensionless), &
             begin, f+f_pfaStru, boolean(), &
             begin, f+f_phaseName, boolean(), &
             begin, f+f_rank, numeric(phyq_dimensionless), &
             begin, f+f_pointingGrids, boolean(), &
             begin, f+f_polygon, boolean() /), &
             continue=.true. )
    call make_tree ( (/ & ! Continuing for s_dump...
             begin, f+f_reset, boolean(), &
             begin, f+f_quantity, vectorQuantity(), &
             begin, f+f_signals, field_spec(s_signal), &
             begin, f+f_spectroscopy, field_type(t_molecule), &
             begin, f+f_stack, boolean(), &
             begin, f+f_start, numeric(), &
             begin, f+f_stop, boolean(), &
             begin, f+f_stopWithError, boolean(), &
             begin, f+f_stride, numeric(), &
             begin, f+f_template, field_spec(s_vectorTemplate,s_quantity), &
             begin, f+f_text, string(), &
             begin, f+f_time, boolean(), &
             begin, f+f_tGrid, field_spec(s_tGrid), &
             begin, f+f_totalMatrixSizes, boolean(), &
             begin, f+f_totalVectorSizes, boolean(), &
             begin, f+f_truncate, boolean(), &
             begin, f+f_variable, n+n_variable_ref, &
             begin, f+f_vector, field_spec(s_vector), &
             begin, f+f_vectorMask, field_spec(s_vector), &
             begin, f+f_vGrid, field_spec(s_vGrid), &
             begin, f+f_ZOT, boolean(), &
             np+n_spec_def/), continue=.true. )
    call make_tree ( (/ & ! Must be AFTER s_vector,s_vectorTemplate, etc.
      begin, s+s_execute, &
             begin, f+f_Clean, boolean(), &
             begin, f+f_command, string(), &
             begin, f+f_crashBurn, boolean(), &
             begin, f+f_delay, numeric(phyq_dimensionless), &
             begin, f+f_fileName, string(), &
             begin, f+f_lines, string(), &
             begin, f+f_options, string(), &
             begin, f+f_stop, boolean(), &
             begin, f+f_stopWithError, boolean(), &
             begin, f+f_wait, boolean(), &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_forwardModelGlobal, &
             begin, f+f_antennaPatterns, string(), &
             begin, f+f_DACSfilterShapes, string(), &
             begin, f+f_filterShapes, string(), &
             begin, f+f_l2pc, string(), &
             begin, f+f_MieTables, string(), &
             begin, f+f_PFAFiles, string(), &
             begin, f+f_pointingGrids, string(), &
             begin, f+f_polygon, string(), &
             np+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_directWriteFile, &
             begin, f+f_file, string(req), &
             begin, f+f_type, field_type(t_outputType,req=req), &
             ndp+n_spec_def, &
      begin, s+s_l1brad, &
             begin, f+f_file, string(), np+n_spec_def, &
      begin, s+s_l1boa, &
             begin, f+f_file, string(), np+n_spec_def, &
      begin, s+s_l2parsf, &
             begin, f+f_file, string(), np+n_spec_def &
             /) )

    call make_tree ( (/ &
      begin, s+s_neuralNet, & ! Must be AFTER s_forwardModel, s_quantity, 
                             !               s_vector and s_matrix
             begin, f+f_measurements, field_spec(s_vector,req=req), &
             begin, f+f_outputSD, field_spec(s_vector,req=req), &
             begin, f+f_state, field_spec(s_vector,req=req), &
             begin, f+f_file, string(), np+n_spec_def &
             /) )

    call make_tree ( (/ &
      begin, s+s_retrieve, & ! Must be AFTER s_forwardModel, s_quantity, 
                             !               s_vector and s_matrix
             begin, f+f_apriori, field_spec(s_vector), &
             begin, f+f_aprioriFraction, field_spec(s_vector), &
             begin, f+f_aprioriScale, numeric(phyq_dimensionless), &
             begin, f+f_average, field_spec(s_matrix), &
             begin, f+f_checkpoint, boolean(), &
             begin, f+f_columnScale, field_type(t_scale), &
             begin, f+f_covariance, field_spec(s_matrix), &
             begin, f+f_covSansReg, boolean(), &
             begin, f+f_diagnostics, field_spec(s_vector), &
             begin, f+f_diagonal, boolean(), &
             begin, f+f_dumpQuantities, vectorQuantity(), &
             begin, f+f_extendedAverage, boolean(), &
             begin, f+f_forwardModel, field_spec(s_forwardModel,req=req), &
             begin, f+f_fuzz, numeric(phyq_dimensionless), & ! Secret
             begin, f+f_fwdModelExtra, field_spec(s_vector,req=req), &
             begin, f+f_fwdModelOut, field_spec(s_vector), &
             begin, f+f_highBound, field_spec(s_vector), &
             begin, f+f_hRegOrders, numeric(phyq_dimensionless), &
             begin, f+f_hRegQuants, field_spec(s_quantity), &
             begin, f+f_hRegWeights, numeric(phyq_dimensionless), &
             begin, f+f_hRegWeightVec, field_spec(s_vector), &
             begin, f+f_jacobian, field_spec(s_matrix), &
             begin, f+f_lambda, numeric(phyq_dimensionless), &
             begin, f+f_lambdaMin, numeric(phyq_dimensionless), &
             begin, f+f_lowBound, field_spec(s_vector), &
             begin, f+f_maxAddDiag, numeric(phyq_dimensionless), &
             begin, f+f_maxJ, numeric(phyq_dimensionless), &
             begin, f+f_measurements, field_spec(s_vector,req=req), &
             begin, f+f_measurementSD, field_spec(s_vector), &
             begin, f+f_method, field_type(t_method), &
             begin, f+f_muMin, numeric(phyq_dimensionless) /) )
    call make_tree ( (/ & ! Continuting for s_retrieve
             begin, f+f_negateSD, boolean(), &
             begin, f+f_outputCovariance, field_spec(s_matrix), &
             begin, f+f_outputSD, field_spec(s_vector), &
             begin, f+f_precisionFactor, numeric(phyq_dimensionless), &
             begin, f+f_regAfter, boolean(), &
             begin, f+f_regApriori, boolean(), &
             begin, f+f_serial, boolean(), &
             begin, f+f_sparseQuantities, field_spec(s_quantity), &
             begin, f+f_state, field_spec(s_vector,req=req), &
             begin, f+f_stateMax, field_spec(s_vector), &
             begin, f+f_stateMin, field_spec(s_vector), &
             begin, f+f_switches, string(), &
             begin, f+f_toggles, string(), &
             begin, f+f_toleranceA, numeric(phyq_dimensionless), &
             begin, f+f_toleranceF, numeric(phyq_dimensionless), &
             begin, f+f_toleranceR, numeric(phyq_dimensionless), &
             begin, f+f_vRegOrders, numeric(phyq_dimensionless), &
             begin, f+f_vRegQuants, field_spec(s_quantity), &
             begin, f+f_vRegWeights, numeric(phyq_dimensionless), &
             begin, f+f_vRegWeightVec, field_spec(s_vector), &
             ndp+n_spec_def /), &
             continue = .true. )
    call make_tree ( (/ &
      begin, s+s_sids, & ! Must be AFTER s_hessian, s_vector and s_matrix
             begin, f+f_forwardModel, field_spec(s_forwardModel,req=req), &
             begin, f+f_fwdModelExtra, field_spec(s_vector), &
             begin, f+f_fwdModelIn, field_spec(s_vector,req=req), &
             begin, f+f_fwdModelOut, field_spec(s_vector,req=req), &
             begin, f+f_destroyJacobian, boolean(), &
             begin, f+f_jacobian, field_spec(s_matrix), &
             begin, f+f_hessian, field_spec(s_hessian), &
             begin, f+f_mirrorHessian, boolean(), &
             begin, f+f_perturbation, field_spec(s_vector), &
             begin, f+f_singleMAF, numeric(phyq_dimensionless), &
             begin, f+f_switches, string(), &
             begin, f+f_TScat, boolean(), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_disjointEquations, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_measurementSD, field_spec(s_vector,req=req), &
             begin, f+f_retrievalExtra, field_spec(s_vector), &
             begin, f+f_retrievalForwardModel, field_spec(s_forwardModel,req=req), &
             begin, f+f_retrievalIn, field_spec(s_vector,req=req), &
             begin, f+f_truthExtra, field_spec(s_vector), &
             begin, f+f_truthForwardModel, field_spec(s_forwardModel,req=req), &
             begin, f+f_truthIn, field_spec(s_vector,req=req), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_normalEquations, & ! Must be AFTER s_vector and s_matrix
             begin, f+f_forwardModel, field_spec(s_forwardModel,req=req), &
             begin, f+f_fwdModelExtra, field_spec(s_vector), &
             begin, f+f_fwdModelIn, field_spec(s_vector,req=req), &
             begin, f+f_fwdModelOut, field_spec(s_vector,req=req), &
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_measurements, field_spec(s_vector,req=req), &
             begin, f+f_measurementSD, field_spec(s_vector), &
             begin, f+f_residualSupplied, boolean(), &
             begin, f+f_rhsOut, field_spec(s_vector,req=req), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_combineChannels, & ! Must be AFTER s_matrix
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_sourceMatrix, field_spec(s_matrix,req=req), &
             ndp+n_spec_def /) )
    call make_tree( (/ &
      begin, s+s_regularization, & ! Must be AFTER s_quantity, s_matrix and
                                   !               s_vector
             begin, f+f_horizontal, boolean(), &
             begin, f+f_matrix, field_spec(s_matrix,req=req), &
             begin, f+f_regOrders, numeric(phyq_dimensionless), &
             begin, f+f_regQuants, field_spec(s_quantity), &
             begin, f+f_regWeightVec, field_spec(s_vector), &
             begin, f+f_regWeights, numeric(phyq_dimensionless), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_streamlineHessian, & ! Must be AFTER s_hessian
             begin, f+f_geodAngle, numeric(phyq_angle), &
             begin, f+f_hessian, field_spec(s_hessian,req=req), &
             begin, f+f_scaleHeight, numeric(phyq_dimensionless), &
             begin, f+f_surface, numeric(phyq_dimensionless), &
             begin, f+f_threshold, numeric(phyq_dimensionless), &
             ndp+n_spec_def, &
      begin, s+s_snoop, &
             begin, f+f_comment, string(), &
             begin, f+f_phaseName, string(), &
             begin, f+f_level, numeric(phyq_dimensionless), &
             begin, f+f_silent, boolean(), &
             begin, f+f_skipDirectWrites, boolean(), &
             begin, f+f_skipRetrieval, boolean(), &
             begin, f+f_stamp, boolean(), &
             ndp+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_dumpblocks, & ! Must be AFTER s_hessian, s_quantity and s_matrix
             begin, f+f_allHessians, boolean(), &
             begin, f+f_allMatrices, boolean(), &
             begin, f+f_colChannels, numeric_or_range(), &
             begin, f+f_colInstances, numeric_or_range(), &
             begin, f+f_colQuantity, field_spec(s_quantity), &
             begin, f+f_colSurfaces, numeric_or_range(), &
             begin, f+f_details, numeric(phyq_dimensionless), &
             begin, f+f_diagonal, boolean(), &
             begin, f+f_hessian, field_spec(s_hessian), &
             begin, f+f_matrix, field_spec(s_matrix), &
             begin, f+f_noAbsent, boolean(), &
             begin, f+f_rowChannels, numeric_or_range(), &
             begin, f+f_rowInstances, numeric_or_range(), &
             begin, f+f_rowQuantity, field_spec(s_quantity), &
             begin, f+f_rowSurfaces, numeric_or_range(), &
             begin, f+f_structure, boolean(), &
             ndp+n_spec_def /) )
    ! Define the relationships between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >
    call make_tree ( (/ &
      begin, z+z_mlsSignals, s+s_module, s+s_band, s+s_radiometer, &
             s+s_signal, s+s_spectrometerType, s+s_time, n+n_section, &
      begin, z+z_spectroscopy, s+s_line, s+s_readSpectroscopy, s+s_spectra, &
             s+s_readIsotopeRatios, s+s_time, s+s_writeSpectroscopy, n+n_section, &
      begin, z+z_globalsettings, &
             begin, p+p_brightObjects, t+t_string, n+n_name_def,&
             begin, p+p_cycle, t+t_string, n+n_name_def, &
             begin, p+p_endtime, t+t_string, n+n_name_def, &
             begin, p+p_IGRF_file, t+t_string, n+n_name_def, &
             begin, p+p_instrument, t+t_instrument, n+n_name_def,&
             begin, p+p_leapsecfile, t+t_string, n+n_name_def,&
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_PFAfile, t+t_string, n+n_name_def,&
             begin, p+p_starttime, t+t_string, n+n_name_def, &
             s+s_binSelector, s+s_dump, s+s_directWriteFile, s+s_dump, &
             s+s_empiricalGeometry, s+s_fGrid, s+s_flushPFA, s+s_forwardModel, &
             s+s_forwardModelGlobal, s+s_l1brad, s+s_l1boa, &
             s+s_l2parsf, s+s_makePFA, s+s_pfaData, s+s_readPFA, &
             s+s_tGrid, s+s_time, s+s_vGrid, s+s_writePFA, n+n_section, &
      begin, z+z_readapriori, s+s_Boolean, s+s_case, s+s_changeSettings, &
             s+s_diff, s+s_dump, s+s_gridded, s+s_isFileAbsent, s+s_time, &
             s+s_l2aux, s+s_l2gp, s+s_readGriddedData, s+s_skip, s+s_snoop, &
             s+s_endSelect, s+s_execute, s+s_reevaluate, s+s_select, n+n_section, &
      begin, z+z_mergegrids, s+s_Boolean, s+s_case, s+s_changeSettings, &
             s+s_concatenate, s+s_concatenateGrids, s+s_ConvertEtaToP,&
             s+s_delete, s+s_diff, s+s_dump, s+s_IsFileAbsent, &
             s+s_endSelect, s+s_execute, s+s_Gridded, s+s_merge, &
             s+s_reevaluate, s+s_select, s+s_skip, s+s_time, s+s_isGridEmpty, &
             s+s_mergeGrids, s+s_vgrid, s+s_wmoTrop, s+s_wmoTropFromGrids, &
             n+n_section /) )
    call make_tree ( (/ &
      begin, z+z_chunkdivide, &
             s+s_time, s+s_chunkDivide, s+s_dump, n+n_section, &
      begin, z+z_construct, s+s_anyGoodRadiances, s+s_anyGoodValues, &
             s+s_Boolean, s+s_catchWarning, s+s_compare, s+s_dump, &
             s+s_forge, s+s_forwardModel, s+s_hgrid, s+s_phase, s+s_quantity, &
             s+s_reevaluate, s+s_changeSettings, s+s_snoop, s+s_time, s+s_vectortemplate, &
             n+n_section /) )
    call make_tree ( (/ &
      begin, z+z_fill, &
             s+s_anyGoodRadiances, s+s_anyGoodValues, s+s_Boolean, &
             s+s_case, s+s_catchWarning, s+s_compare, s+s_computeTotalPower, &
             s+s_chunkEndsBefore, s+s_chunkStartsAfter, s+s_delete, &
             s+s_concatenate, s+s_concatenateGrids, s+s_ConvertEtaToP,&
             s+s_destroy, s+s_diff, s+s_directRead, s+s_dump, s+s_endSelect, &
             s+s_execute, s+s_fill, s+s_fillCovariance, s+s_fillDiagonal, &
             s+s_flagcloud, s+s_flushL2PCBins, s+s_flushPFA, s+s_Gridded, &
             s+s_hessian, s+s_isGridEmpty, s+s_load, s+s_matrix, s+s_mergeGrids, s+s_negativePrecision, &
             s+s_phase, s+s_populateL2PCBin, s+s_readGriddedData, s+s_reevaluate, &
             s+s_repeat, s+s_restrictRange, s+s_select, s+s_changeSettings, &
             s+s_skip, s+s_snoop, s+s_streamlineHessian, s+s_subset, &
             s+s_time, s+s_transfer, s+s_updateMask, s+s_vector, s+s_wmoTropFromGrids, &
             n+n_section, &
      begin, z+z_retrieve, s+s_anyGoodValues, s+s_case, s+s_catchWarning, &
             s+s_checkpoint, s+s_compare, s+s_diff, s+s_dump, s+s_dumpBlocks, &
             s+s_endSelect, s+s_flagCloud, s+s_flushPFA, s+s_leakcheck, &
             s+s_neuralNet, s+s_reevaluate, s+s_repeat, s+s_restrictRange, &
             s+s_retrieve, s+s_select, s+s_sids, s+s_skip, s+s_snoop, &
             s+s_subset, s+s_time, s+s_updateMask, &
             n+n_section, &
      begin, z+z_join, s+s_time, s+s_label, s+s_l2gp, s+s_l2aux, &
             s+s_case, s+s_directWrite, s+s_diff, s+s_dump, s+s_endSelect, &
             s+s_execute, s+s_select, s+s_skip, n+n_section, &
      begin, z+z_algebra, s+s_columnScale, s+s_combineChannels, s+s_cyclicJacobi, &
             s+s_disjointEquations, s+s_normalEquations, s+s_reflect, &
             s+s_regularization, s+s_rowScale, nc+n_section, &
      begin, z+z_output, s+s_Boolean, s+s_case, s+s_catenate, s+s_copy, &
             s+s_destroy, s+s_diff, s+s_dump, s+s_dumpblocks, s+s_endSelect, &
             s+s_execute, s+s_hgrid, s+s_isFileAbsent, s+s_isSwathEmpty, &
             s+s_output, s+s_Reevaluate, &
             s+s_select, s+s_Skip, s+s_Sleep, s+s_time, s+s_writeFileAttribute, s+s_subSample, &
             n+n_section /) )

  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

    ! --------------------------------------------------  Boolean  -----
    pure function Boolean ( Req )
      ! Declare Boolean type
      logical, intent(in), optional :: Req     ! Field is required if true
      integer :: Boolean(2)
      boolean = field_type ( t_boolean, req )
    end function Boolean

    ! --------------------------------------------------  Numeric  -----
    pure function Numeric ( Unit, Req, Scalar )
      ! Declare numeric type field, with Unit units if present
      integer, intent(in), optional :: Unit    ! phyq_... from intrinsic module
      logical, intent(in), optional :: Req     ! Field is required if true
      logical, intent(in), optional :: Scalar  ! Field is scalar if true
      integer :: Numeric(2)
      integer :: UnitPart
      unitPart = 0
      if ( present(unit) ) unitPart = d*u*unit
      numeric = (/ t+t_numeric, node(req,scalar)+n_field_type+unitPart /)
    end function Numeric

    ! -----------------------------------------  Numeric_Or_Range  -----
    pure function Numeric_Or_Range ( Unit, Req, Scalar )
      ! Declare numeric type field, with Unit units if present
      integer, intent(in), optional :: Unit    ! phyq_... from intrinsic module
      logical, intent(in), optional :: Req     ! Field is required if true
      logical, intent(in), optional :: Scalar  ! Field is scalar if true
      integer :: Numeric_Or_Range(3)
      integer :: UnitPart
      unitPart = 0
      if ( present(unit) ) unitPart = d*u*unit
      numeric_or_range = (/ t+t_numeric, t+t_numeric_range, &
                         &  node(req,scalar)+n_field_type+unitPart /)
    end function Numeric_Or_Range

    ! --------------------------------------------  Numeric_Range  -----
    pure function Numeric_Range ( Unit, Req, Scalar )
      ! Declare numeric type field, with Unit units if present
      integer, intent(in), optional :: Unit    ! phyq_... from intrinsic module
      logical, intent(in), optional :: Req     ! Field is required if true
      logical, intent(in), optional :: Scalar  ! Field is scalar if true
      integer :: Numeric_Range(2)
      integer :: UnitPart
      unitPart = 0
      if ( present(unit) ) unitPart = d*u*unit
      numeric_range = (/ t+t_numeric_range, node(req,scalar)+n_field_type+unitPart /)
    end function Numeric_Range

    ! ---------------------------------------------------  String  -----
    pure function String ( Req )
      ! Declare String type field
      logical, intent(in), optional :: Req  ! Field is required if true
      integer :: String(2)
      string = field_type ( t_string, req )
    end function String

    ! -------------------------------------------  VectorQuantity  -----
    pure function VectorQuantity ( Req, Scalar, Expr )
      ! Declare Vector Quantity "dot" field
      logical, intent(in), optional :: Req  ! Field is required if true
      logical, intent(in), optional :: Scalar  ! Field is scalar if true
      logical, intent(in), optional :: Expr ! Expr is OK if true
      integer :: VectorQuantity(4)
      VectorQuantity = (/ s+s_vector, f+f_template, f+f_quantities, &
                       &  node(req,scalar=scalar,expr=expr)+n_dot /)
    end function VectorQuantity

  end subroutine INIT_TABLES

  ! We would prefer these to be internal to Init_Tables when
  ! compilers support internal generic functions

  ! -------------------------------------------  Field_Spec_Array  -----
  pure function Field_Spec_Array ( Spec, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec(:)          ! S_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Scalar ! Field is scalar if true
    integer :: Field_Spec_Array(1+size(spec))
    field_spec_array = (/ s+spec, node(req,scalar)+n_field_spec /)
  end function Field_Spec_Array

  ! -----------------------------------------------  Field_Spec_1  -----
  pure function Field_Spec_1 ( Spec, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec             ! S_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Scalar ! Field is scalar if true
    integer :: Field_Spec_1(2)
    field_spec_1 = (/ s+spec, node(req,scalar)+n_field_spec /)
  end function Field_Spec_1

  ! -----------------------------------------------  Field_Spec_2  -----
  pure function Field_Spec_2 ( Spec1, Spec2, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec1, Spec2     ! S_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Scalar ! Field is scalar if true
    integer :: Field_Spec_2(3)
    field_spec_2 = (/ s+spec1, s+spec2, node(req,scalar)+n_field_spec /)
  end function Field_Spec_2

  ! -----------------------------------------------  Field_Spec_3  -----
  pure function Field_Spec_3 ( Spec1, Spec2, Spec3, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec1, Spec2, Spec3 ! S_...
    logical, intent(in), optional :: Req       ! Field is required if true
    logical, intent(in), optional :: Scalar    ! Field is scalar if true
    integer :: Field_Spec_3(4)
    field_spec_3 = (/ s+spec1, s+spec2, s+spec3, &
                   &  node(req,scalar)+n_field_spec /)
  end function Field_Spec_3

  ! -----------------------------------------------  Field_Spec_4  -----
  pure function Field_Spec_4 ( Spec1, Spec2, Spec3, Spec4, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec1, Spec2, Spec3, Spec4 ! S_...
    logical, intent(in), optional :: Req       ! Field is required if true
    logical, intent(in), optional :: Scalar    ! Field is scalar if true
    integer :: Field_Spec_4(5)
    field_spec_4 = (/ s+spec1, s+spec2, s+spec3,  s+spec4, &
                   &  node(req,scalar)+n_field_spec /)
  end function Field_Spec_4

  ! -----------------------------------------------  Field_Spec_5  -----
  pure function Field_Spec_5 ( Spec1, Spec2, Spec3, Spec4, Spec5, Req, Scalar )
    ! Declare field-spec field
    use Tree_types, only: n_field_spec
    integer, intent(in) :: Spec1, Spec2, Spec3, Spec4, Spec5 ! S_...
    logical, intent(in), optional :: Req       ! Field is required if true
    logical, intent(in), optional :: Scalar    ! Field is scalar if true
    integer :: Field_Spec_5(6)
    field_spec_5 = (/ s+spec1, s+spec2, s+spec3,  s+spec4,  s+spec5, &
                   &  node(req,scalar)+n_field_spec /)
  end function Field_Spec_5

  ! -------------------------------------------  Field_Type_Array  -----
  pure function Field_Type_Array ( Type, Req, Array )
    ! Declare field-Type field
    use Tree_types, only: n_field_type
    integer, intent(in) :: Type(:)          ! T_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Array  ! Array element arrays OK if true
    integer :: Field_Type_Array(1+size(Type))
    field_type_array =  (/ t+Type, node(req,array)+n_field_type /)
  end function Field_Type_Array

  ! -----------------------------------------------  Field_Type_1  -----
  pure function Field_Type_1 ( Type, Req, Empty, Array )
    ! Declare field-Type field
    use Tree_types, only: n_field_type
    integer, intent(in) :: Type             ! T_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Empty  ! Empty is OK if true
    logical, intent(in), optional :: Array  ! Array elements array OK if true
    integer :: Field_Type_1(2)
    field_type_1 =  (/ t+Type, node(req,empty=empty,array=array)+n_field_type /)
  end function Field_Type_1

  ! -----------------------------------------------  Field_Type_2  -----
  pure function Field_Type_2 ( Type1, Type2, Req, Array )
    ! Declare field-Type field
    use Tree_types, only: n_field_type
    integer, intent(in) :: Type1, Type2     ! T_...
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Array  ! Array elements array OK if true
    integer :: Field_Type_2(3)
    field_type_2 =  (/ t+type1, t+type2, node(req,array=array)+n_field_type /)
  end function Field_Type_2

  ! -------------------------------------------------------  Node  -----
  pure function Node ( Req, Scalar, Empty, Expr, Array )
    ! Compute node generator = N if Req is absent or false, or NR
    ! if Req is present and true
    use Intrinsic, only: Array_Array, D, Empty_OK, Expr_OK, No_Array, Req_Fld
    logical, intent(in), optional :: Req    ! Field is required if true
    logical, intent(in), optional :: Scalar ! Scalar is required if true
    logical, intent(in), optional :: Empty  ! Empty is OK if true
    logical, intent(in), optional :: Expr   ! Expr is OK if true
    logical, intent(in), optional :: Array  ! Array elements arrays OK if true
    integer :: Node
    node = n
    if ( present(req) ) node = node + merge(d*req_fld,0,req)
    if ( present(scalar) ) node = node + merge(d*no_array,0,scalar)
    if ( present(empty) ) node = node + merge(d*empty_OK,0,empty)
    if ( present(expr) ) node = node + merge(d*expr_OK,0,expr)
    if ( present(array) ) node = node + merge(d*array_array,0,array)
  end function Node

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module INIT_TABLES_MODULE

! $Log$
! Revision 2.659  2021/10/14 22:24:15  pwagner
! The NeuralNet command now requires the field outputSD
!
! Revision 2.658  2021/02/05 05:17:50  pwagner
! Must permit NeuralNet command in Retrieve section
!
! Revision 2.657  2020/12/22 22:21:40  pwagner
! Added NeuralNet command
!
! Revision 2.656  2020/07/29 23:43:50  pwagner
! May utilize BooleanFromEmptyGrid in Fill sections
!
! Revision 2.655  2020/07/28 20:34:31  vsnyder
! Allow Molecules field of ForwardModel to have array elements that are arrays
!
! Revision 2.654  2020/07/23 22:18:40  pwagner
! AllowEmptyGrids is a new field in Concatenate; defaults to pre-v5.1 behavior
!
! Revision 2.653  2020/07/22 22:46:39  pwagner
! Added chunkEndsBefore, ChunkStartsAfter commands
!
! Revision 2.652  2020/07/09 23:53:48  pwagner
! Many cmds from readApriori and MergeGrids phase now available to Fill phase
!
! Revision 2.651  2019/10/16 20:54:28  pwagner
! Subset command may take a MissingValue field
!
! Revision 2.650  2019/10/07 20:05:16  vsnyder
! Add trapezoid field to ForwardModel, for quadrature in FullForwardModel
!
! Revision 2.649  2019/10/03 17:29:52  pwagner
! Convert from eta levels may now take a vGrid field
!
! Revision 2.648  2019/09/23 20:39:27  pwagner
! Conversion from Eta surfaces may optionally be logarithmic
!
! Revision 2.647  2019/04/24 19:15:25  vsnyder
! Add MIFTangent field to forwardModel
!
! Revision 2.646  2018/11/28 21:09:54  pwagner
! Now the surface field may be a numeric instead of a range
!
! Revision 2.645  2018/10/19 00:03:47  pwagner
! inputFile= with /pause pauses to read from inputFile instead of stdin
!
! Revision 2.644  2018/10/17 23:05:04  pwagner
! Dump command can take /pause field to wait for user input for, e.g. debugging
!
! Revision 2.643  2018/09/13 20:22:02  pwagner
! Added /reset flag to phase command
!
! Revision 2.642  2018/09/06 23:59:43  pwagner
! More commands that set a runtime flag can now take /reverse
!
! Revision 2.641  2018/04/19 23:43:39  pwagner
! Skip may take /nextChunk flag
!
! Revision 2.640  2018/03/23 20:15:01  pwagner
! Added type field to isFileAbsent; skip and reevaluate commands to readapriori section
!
! Revision 2.639  2018/03/22 18:13:49  pwagner
! Added command IsFileAbsent; may occur in ReadApriori, MergeGrids, and Output sections
!
! Revision 2.638  2018/03/14 22:18:15  pwagner
! May changeSettings in readApriori and MergeGrids sections
!
! Revision 2.637  2018/02/23 22:08:12  mmadatya
! Added l_instECR for ASMLS
!
! Revision 2.636  2017/12/15 18:29:56  mmadatya
! Added geoHeight as new quantity template; heightFromPressure as new Fill method
!
! Revision 2.635  2017/07/27 16:38:55  pwagner
! Reevaluate may take /truncate; label may take Boolean; DirectWrite may take Boolean or groupName
!
! Revision 2.634  2017/07/10 18:48:51  pwagner
! Transfer may /expandMask to all masking bits; may /skipValues to transfeer only mask; Fill may replaceMissingValue=
!
! Revision 2.633  2017/04/07 18:45:56  pwagner
! sourceType used when choosing how to Fill ascenddescend
!
! Revision 2.632  2017/03/07 21:20:30  pwagner
! Support new meteorology origin: merra_2
!
! Revision 2.631  2017/02/22 01:23:07  pwagner
! New /resetAll switch clears all Masking bits
!
! Revision 2.630  2017/02/10 21:55:20  pwagner
! Added the polygon method for ChunkDivide; /sourcemask field for Fill
!
! Revision 2.628  2016/10/21 23:27:09  vsnyder
! Remove unused USE name
!
! Revision 2.627  2016/09/02 00:24:54  vsnyder
! Specify units for QTM_Level, SolarZenith, and Time fields in HGrid
!
! Revision 2.626  2016/08/30 23:54:22  vsnyder
! Add ECR as a geolocation type, and TngtECR and TngtGeodLat as quantity types
!
! Revision 2.625  2016/08/16 23:23:07  vsnyder
! Remove phyq_dimensionless requirement from QTM_Level of HGrid
!
! Revision 2.624  2016/04/01 00:26:23  pwagner
! May now Execute a single command or a script of lines from l2cf
!
! Revision 2.623  2016/02/26 02:07:44  vsnyder
! Add QTM support
!
! Revision 2.622  2016/01/29 01:09:24  vsnyder
! Add polygon file to ForwardModelGlobal and Polygon switch to Dump
!
! Revision 2.621  2015/12/01 21:19:57  pwagner
! May Fill with nearest profile number
!
! Revision 2.620  2015/09/30 20:32:06  pwagner
! With /noZeros field negativePrecision command now resets 0 to -1; Catenate can repair geoLocations
!
! Revision 2.619  2015/09/22 01:56:51  vsnyder
! Add GHzAzim and ScECR quantity types, ReferenceMIF and Regular fill fields
!
! Revision 2.618  2015/09/17 22:58:54  pwagner
! Added changeSettings command
!
! Revision 2.617  2015/08/25 17:27:16  vsnyder
! Allow PhiWindow to be a number or a range
!
! Revision 2.616  2015/07/15 17:06:45  pwagner
! inputFile not a required field to DirectWrite
!
! Revision 2.615  2015/07/14 23:29:50  pwagner
! may divert Dump commands to named dumpFile; /truncate field; label and inputFile fields in DirectWrite
!
! Revision 2.614  2015/06/04 00:52:11  pwagner
! Needed wmoTrop among perissible sourceGrid in Fill
!
! Revision 2.613  2015/06/03 23:10:09  pwagner
! NAG hated Field_Spec_Array; workaround
!
! Revision 2.612  2015/04/21 17:44:27  pwagner
! May DirectRead, DirectWrite files with /noPCFid even when usingPCF; may DirectRead qty geolocations
!
! Revision 2.611  2015/04/09 20:53:22  pwagner
! rank may be specified for Dump and DirectRead commands
!
! Revision 2.610  2015/03/31 21:00:32  pwagner
! rank is a new field for DirectWrite-ing quantity values as, say, rank3
!
! Revision 2.609  2015/03/28 02:55:04  vsnyder
! Paul added Quantity to t_OutputType.  Deleted Azimuth from t_QuantityType.
! Added f_Truncate -- should have been save -- to s_Forge.  Added
! f_Coordinate, f_Stacked, and f_xGrid to s_Quantity.
!
! Revision 2.608  2015/02/27 23:13:32  pwagner
! May Dump global attributes
!
! Revision 2.607  2014/10/27 23:06:19  pwagner
! /additional now a possible field for writeFileAttribute
!
! Revision 2.606  2014/10/06 23:50:02  pwagner
! May now write added file attributes
!
! Revision 2.605  2014/09/29 20:18:14  vsnyder
! Add NoMagneticField switch to ForwardModel
!
! Revision 2.604  2014/09/05 01:25:47  vsnyder
! Add TotalMatrixSizes, TotalVectorSizes to Dump command
!
! Revision 2.603  2014/08/26 23:22:35  pwagner
! We may Dump memory, time; the /reset option may Dump changes in their use
!
! Revision 2.602  2014/08/12 23:09:53  pwagner
! /additional flag added to phase commands applies to options field
!
! Revision 2.601  2014/06/11 20:04:11  pwagner
! New concatenateGrids command; grid field in Concatenated renamed sourceGrid
!
! Revision 2.600  2014/05/31 00:22:18  pwagner
! We wont always have units for solarZenith field
!
! Revision 2.599  2014/04/24 23:55:32  pwagner
! May set master coordinate in hGrid specification
!
! Revision 2.598  2014/04/07 18:03:55  pwagner
! May specify AscDescMode when DirectWrite-ing swaths
!
! Revision 2.597  2014/03/13 18:11:29  pwagner
! c= field may contain units
!
! Revision 2.596  2014/03/01 03:10:33  vsnyder
! Delete FrequencyGrid, specify units for many numeric fields
!
! Revision 2.595  2013/11/18 22:26:30  pwagner
! phase spec takes optional /debug /verbose fields
!
! Revision 2.594  2013/10/24 21:04:30  pwagner
! /dontLatch in profile Fills retains profile heights as input
!
! Revision 2.593  2013/10/15 23:52:06  pwagner
! May copy quantity values to a file global attribute
!
! Revision 2.592  2013/10/09 23:44:39  vsnyder
! Add Variable field for Dump
!
! Revision 2.591  2013/10/08 23:51:03  pwagner
! Removed unused p_types from ChunkDivide; may have multiple Output sections
!
! Revision 2.590  2013/10/02 00:48:03  pwagner
! Added ascenddescend Fill Method
!
! Revision 2.589  2013/09/21 00:37:20  vsnyder
! Add ChunkDivide field for Dump
!
! Revision 2.588  2013/09/21 00:23:59  pwagner
! Added geoid Fill methods
!
! Revision 2.587  2013/09/19 23:36:58  vsnyder
! Require Vector Quantity in Fill's expr to be scalar
!
! Revision 2.586  2013/09/18 18:29:07  pwagner
! Corrected placement of continue=.true. in s_dump
!
! Revision 2.585  2013/09/17 23:02:19  pwagner
! Added Scatter, Gather methods
!
! Revision 2.584  2013/09/17 00:54:12  vsnyder
! Alphabetize data_type_indices and data type trees
!
! Revision 2.583  2013/09/04 17:36:21  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.582  2013/08/29 19:37:04  pwagner
! May Dump plain stack, too
!
! Revision 2.581  2013/08/17 02:55:44  vsnyder
! Remove unused USE named
!
! Revision 2.580  2013/08/17 00:19:27  pwagner
! Fixed HGrids may have geolocations overridden
!
! Revision 2.579  2013/08/16 02:34:46  vsnyder
! Remove Model_Plane_MIF
!
! Revision 2.578  2013/08/09 01:03:59  vsnyder
! Add ReferenceMIF component
!
! Revision 2.577  2013/08/03 00:40:42  vsnyder
! Require vector.quantity in dumpQuantities field in retrieve spec
!
! Revision 2.576  2013/07/25 00:24:25  vsnyder
! Add MIFRHI, replace TransformRHI with TransformMIFRHI
!
! Revision 2.575  2013/07/19 01:20:37  vsnyder
! Sort some stuff, add TransformRHI field
!
! Revision 2.574  2013/07/18 01:10:57  vsnyder
! Remove scVel since it's ambiguous whether it's ECI or ECR, and nobody
! uses it anyway.
!
! Revision 2.573  2013/07/15 23:53:39  vsnyder
! Allow repeated fields on directWrite commands
!
! Revision 2.572  2013/07/12 23:45:30  vsnyder
! Functions for almost all field definitions
!
! Revision 2.571  2013/07/12 23:25:28  vsnyder
! Bogus checkin: Remove unreferenced error messages
!
! Revision 2.570  2013/07/03 23:06:27  pwagner
! Added new array-constructor functions for frequently-used idioms
! 
! Revision 2.569  2013/06/21 17:37:39  pwagner
! /crashIfPhiNotMono flag added to ChunkDivide config; default is to just warn; removed -Snmono switch
!
! Revision 2.568  2013/05/31 00:41:51  vsnyder
! Add geolocation field to fill
!
! Revision 2.567  2013/05/22 20:18:35  pwagner
! Values now a field for reevaluate; may destroy r/t Booleans
!
! Revision 2.566  2013/05/21 23:52:47  vsnyder
! Add MIFExtinctionExtrapolation and MIFExtinctionForm
!
! Revision 2.565  2013/05/17 00:48:49  pwagner
! May constrain Transfer command command to quantitynames by r/t Boolean
!
! Revision 2.564  2013/05/07 22:01:30  pwagner
! run-time Booleans can store arrays of values, formulas can evaluate named terms
!
! Revision 2.563  2013/04/24 00:35:47  pwagner
! Added inputBoolean and made Reevaluate formulas more powerful
!
! Revision 2.562  2013/04/22 17:46:23  pwagner
! Reevaluate may store a literal instead of a Boolean value
!
! Revision 2.561  2013/04/17 00:04:45  pwagner
! Added new Repeat control structure to Fill, Retrieve sections
!
! Revision 2.560  2013/03/15 20:36:03  vsnyder
! Add 'switches' field to SIDS
!
! Revision 2.559  2013/01/17 19:59:47  pwagner
! Preparing new fields for Tranfer command
!
! Revision 2.558  2013/01/02 21:40:33  pwagner
! Added derivative method to Fill command; Transfer can do Fill methods, too
!
! Revision 2.557  2012/11/14 00:57:26  pwagner
! Use dimList for choosing which of {csi} to average over
!
! Revision 2.556  2012/11/08 21:08:03  pwagner
! Fill command may take dimList as a field
!
! Revision 2.555  2012/10/22 18:13:28  pwagner
! Many Subset operations now available in Fill
!
! Revision 2.554  2012/10/17 00:42:55  pwagner
! removed unused sourcesGrid
!
! Revision 2.553  2012/10/09 00:48:30  pwagner
! New ignoreTemplate, changed force meaning in Fill
!
! Revision 2.552  2012/08/30 23:03:07  vsnyder
! Add lambdaMin
!
! Revision 2.551  2012/07/05 23:50:09  pwagner
! Copy command in Output may contain options field
!
! Revision 2.550  2012/07/04 01:53:46  vsnyder
! Add dumpQuantities field of retrieve spec
!
! Revision 2.549  2012/06/27 18:00:07  pwagner
! May overwrite command line options with options field to phase spec
!
! Revision 2.548  2012/06/07 22:44:21  pwagner
! May remove Quantities, vectorTemplates
!
! Revision 2.547  2012/06/06 20:37:42  vsnyder
! Add toggles field to retrieve spec
!
! Revision 2.546  2012/05/24 21:06:22  vsnyder
! Add template field to VectorTemplate spec
!
! Revision 2.545  2012/05/11 00:18:38  pwagner
! Added isSwathEmpty to set Boolean in Output; we can Skip Copy of OH when THz is off
!
! Revision 2.544  2012/05/10 00:46:35  pwagner
! Output section can have l2cf-control stuctures
!
! Revision 2.543  2012/05/08 17:47:34  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.542  2012/05/01 23:15:39  pwagner
! May use formula to decide whether to skip
!
! Revision 2.541  2012/05/01 22:20:27  vsnyder
! Add radianceQuantity field to subset
!
! Revision 2.540  2012/04/26 23:30:40  pwagner
! May Dump chunk number, phase name
!
! Revision 2.539  2012/04/25 20:32:24  pwagner
! Inserting missing profiles after chunk end now an option controlled by 'extendible' field
!
! Revision 2.538  2012/04/20 00:44:50  pwagner
! May dump callStack from within l2cf
!
! Revision 2.537  2012/03/15 22:49:45  vsnyder
! Add IGRF_file parameter, igrf field in dump command, some cannonball polishing
!
! Revision 2.536  2012/03/07 02:01:07  vsnyder
! Remove fwmState and fwmJacobian, add transformMIFextinction
!
! Revision 2.535  2012/02/24 21:13:54  pwagner
! Add /interpolate to DirectRead, /single to DirectWrite
!
! Revision 2.534  2012/02/08 23:18:57  pwagner
! May Load Hessian from l2pc database
!
! Revision 2.533  2012/02/02 00:55:10  pwagner
! Added bin and hessian fields to DirectRead command
!
! Revision 2.532  2012/02/01 00:17:58  vsnyder
! Remove Matrix from the Retrieve section; nobody used it, because it was broken.
!
! Revision 2.531  2012/01/30 18:50:09  pwagner
! Fixed bugs preventing us from Diffing or Dumping matrices
!
! Revision 2.530  2012/01/05 01:19:49  pwagner
! Sets LAST_AUTO_LIT to last literal
!
! Revision 2.529  2011/12/21 01:38:41  vsnyder
! Add LowestRetrievedPressure, MIFExtinction[v2], fwmJacobian, and fwmState
!
! Revision 2.528  2011/12/15 01:50:49  pwagner
! Added sdName and /spread fields to DirectRead
!
! Revision 2.527  2011/11/18 23:41:44  pwagner
! Workaround for faulty dependency calculation of srclib/tree_checker
!
! Revision 2.526  2011/11/01 21:04:44  pwagner
! Added DirectRead as a Fill Section command to read from hdf5 files
!
! Revision 2.525  2011/10/25 22:12:19  pwagner
! Removed unused Template spec; added readIsotopeRatios
!
! Revision 2.524  2011/06/16 23:15:39  pwagner
! Added /downsample to reduce resolution of gridded data when read
!
! Revision 2.523  2011/06/16 20:51:03  vsnyder
! Uncomment expr field in Fill
!
! Revision 2.522  2011/06/16 20:23:10  vsnyder
! Add expr field in Fill, as comments
!
! Revision 2.521  2011/06/02 19:24:34  pwagner
! May dump allRadiometers
!
! Revision 2.520  2011/05/05 15:23:52  pwagner
! Added readGriddedData command to readApriori section
!
! Revision 2.519  2011/04/27 17:38:11  pwagner
! Added new command (not a named spec) wmoTropFromGrids
!
! Revision 2.518  2011/04/20 16:53:11  pwagner
! Added new flexibility to l2cf control flow by run-time booleans affecting gridded data
!
! Revision 2.517  2011/03/22 23:47:21  pwagner
! May now reshape qty template field while filling explicitly
!
! Revision 2.516  2011/03/15 22:52:29  pwagner
! May now modify quantity template fields with fill method
!
! Revision 2.515  2011/03/10 21:39:33  pwagner
! May now specify time in explicit hGrids
!
! Revision 2.514  2010/11/30 22:00:57  pwagner
! Reduced number of continuation lines
!
! Revision 2.513  2010/11/20 00:00:41  pwagner
! May specifiy surfaces gap beyond which to zero out in Streamline
!
! Revision 2.512  2010/09/17 00:09:55  pwagner
! Can constrain writing l2pc blocks by name
!
! Revision 2.511  2010/08/16 19:16:02  yanovsky
! Add f_atmos_second_der, f_moleculeSecondDerivatives
!
! Revision 2.510  2010/08/13 22:08:17  pwagner
! May diff hessians, matrices
!
! Revision 2.509  2010/08/06 23:01:01  pwagner
! May dump some or all matrices, hessians
!
! Revision 2.508  2010/07/22 17:38:41  pwagner
! Replaced method=special fills with unique names
!
! Revision 2.507  2010/07/06 16:04:32  pwagner
! Fixed bug in a, b fields of transfer command: are vectors, not quantities
!
! Revision 2.506  2010/07/01 00:49:01  pwagner
! Transfer between vectors may now also manipulate
!
! Revision 2.505  2010/06/07 23:28:50  vsnyder
! Add TScatMolecules, TScatMoleculeDerivatives, Use_Tscat.  Change
! PhaseFrqTol to FrqTol.
!
! Revision 2.504  2010/04/30 22:56:44  vsnyder
! Add TScat to bin selector type
!
! Revision 2.503  2010/04/28 00:14:52  pwagner
! May specify instances range in Subset
!
! Revision 2.502  2010/04/22 23:36:21  pwagner
! May fill num rads/MIF as a percentage
!
! Revision 2.501  2010/04/16 01:39:34  vsnyder
! Added /allFiles and file fields to dump file database
!
! Revision 2.500  2010/03/26 23:12:12  vsnyder
! Add ignoreHessian field to forwardModel
!
! Revision 2.499  2010/03/24 20:53:10  vsnyder
! Add AllHessians and Hessians to DumpBlocks, allow DumpBlocks in Output section
!
! Revision 2.498  2010/03/02 01:08:59  pwagner
! Added geodAltitude as a quantity type
!
! Revision 2.497  2010/02/25 18:21:44  pwagner
! Adds support for new Hessian data type
!
! Revision 2.496  2010/01/23 01:02:37  vsnyder
! Remove LogIWC
!
! Revision 2.495  2009/10/27 22:13:09  pwagner
! Remedied omission of s_diff from commands allowed by sections
!
! Revision 2.494  2009/10/26 17:10:37  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.493  2009/09/25 02:37:24  vsnyder
! Add badValue, keepChannels
!
! Revision 2.492  2009/09/19 00:33:44  vsnyder
! Add LogIWC
!
! Revision 2.491  2009/09/15 20:03:39  pwagner
! Dump commands take boolean fields /stop, /stopWithError, /crashBurn
!
! Revision 2.490  2009/08/24 20:13:16  pwagner
! May Fill H2O precision from RHI precision
!
! Revision 2.489  2009/08/17 16:39:36  pwagner
! read_apriori may read DELP field and perform sum
!
! Revision 2.488  2009/06/26 00:16:58  pwagner
! May now copy ascii file to DGM calling input file type 'ascii' insstead of 'l2cf'
!
! Revision 2.487  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.486  2009/04/23 23:01:59  pwagner
! May specify upperOverlap or lowerOverlap in DirectWrites
!
! Revision 2.485  2009/04/16 21:54:40  pwagner
! /exact keyword in status Fill to fix radiance bug
!
! Revision 2.484  2009/04/13 20:45:57  pwagner
! heightRange in explicit Fill can fill above or below specified height
!
! Revision 2.483  2009/03/14 02:45:34  honghanh
! Add dnwt_count and dwnt_abandoned
!
! Revision 2.482  2008/12/18 21:06:25  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.481  2008/12/02 23:27:41  pwagner
! May automatically label every quantity in a vector now
!
! Revision 2.480  2008/11/06 21:50:19  pwagner
! Fill method swapValues swaps values between two quantities
!
! Revision 2.479  2008/10/03 16:39:52  livesey
! Added extinctionv2
!
! Revision 2.478  2008/09/30 22:31:40  vsnyder
! Change TGrid, add IWC, delete AuxGrid
!
! Revision 2.477  2008/09/19 23:54:45  pwagner
! May now Destroy GriddedData
!
! Revision 2.476  2008/09/17 23:19:24  pwagner
! Allow date string in gridded data to offset gmao background files
!
! Revision 2.475  2008/08/22 17:28:33  pwagner
! fields may not be duplicated when declaring new Quantity in Construct
!
! Revision 2.474  2008/08/21 23:42:45  vsnyder
! Remove GenerateTScat from ForwardModel; use TScat on Sids
!
! Revision 2.473  2008/08/14 20:58:17  pwagner
! /interpolate now possible field in Transfer command
!
! Revision 2.472  2008/07/30 19:08:00  vsnyder
! Add PhaseFrqTol field to ForwardModel spec
!
! Revision 2.471  2008/06/06 21:03:59  michael
! changes for fill method uncompressradiance
!
! Revision 2.470  2008/06/06 01:57:31  vsnyder
! Add ScatteringAngle to QuantityType, add TScat to SIDS
!
! Revision 2.469  2008/06/05 02:11:57  vsnyder
! Added cloudTemperature and TScat to quantityType.
! Added dimensionless, dimless and iceDensity to vGridCoord.
! Moved some stuff to satisfy order dependencies on vGrid.
! Added AuxGrid field to Quantity and MieTables field to Dump.
!
! Revision 2.468  2008/05/28 21:55:34  pwagner
! geo location Fill method to fill chunk numbers[maf]
!
! Revision 2.467  2008/05/22 00:26:01  pwagner
! Added criticalBands to Chunk Divide
!
! Revision 2.466  2008/05/20 00:25:05  vsnyder
! Add MieTables field for ForwardModelGlobal
!
! Revision 2.465  2008/04/28 21:46:34  pwagner
! Needed to add ident of new total power fill command
!
! Revision 2.464  2008/04/26 00:39:16  livesey
! Added total power stuff
!
! Revision 2.463  2008/04/11 01:17:00  livesey
! Added uncompressRadiance fill
!
! Revision 2.462  2007/12/07 01:12:14  pwagner
! Lets us catch warnings and assign to runtime Booleans
!
! Revision 2.461  2007/11/15 22:06:08  pwagner
! New Compare command, and others giving value to runtimeBooleans, now in Join, Retrieve sections
!
! Revision 2.460  2007/11/08 03:24:39  vsnyder
! Add stateMax and stateMin
!
! Revision 2.459  2007/11/07 03:10:48  vsnyder
! Add pathNorm field to forward model config
!
! Revision 2.458  2007/11/05 18:36:07  pwagner
! May Skip remaining lines in Fill, Join, Retrieve sections depending on Boolean
!
! Revision 2.457  2007/10/09 00:32:05  pwagner
! Added ability to dump masks of quantities, vectors
!
! Revision 2.456  2007/10/02 22:51:54  vsnyder
! Add checkpoint stuff
!
! Revision 2.455  2007/08/23 22:17:05  pwagner
! manipulation Fills can now use statistical functions
!
! Revision 2.454  2007/08/17 00:34:58  pwagner
! MergeGrids section may now dump
!
! Revision 2.453  2007/01/11 20:48:30  vsnyder
! Add SurfaceHeight to gridded data, vector quantities, allow dump in ReadApriori
!
! Revision 2.452  2006/10/11 00:18:17  pwagner
! Changes to permit new convergence ratio field in L2GP
!
! Revision 2.451  2006/10/02 23:05:03  pwagner
! May Fill chi^2 ratio to measure convergence
!
! Revision 2.450  2006/09/19 20:37:46  vsnyder
! Add aprioriFraction to retrieve
!
! Revision 2.449  2006/08/11 20:58:38  vsnyder
! Add 'simple' method to use alternate Newton solver
!
! Revision 2.448  2006/08/04 18:09:46  vsnyder
! Add LeakCheck command, /allMatrices and /allVectors fields to Destroy command
!
! Revision 2.447  2006/08/02 19:53:50  vsnyder
! Add destroy command and field for output
!
! Revision 2.446  2006/07/19 22:28:17  vsnyder
! Add /allMatrices, details= and /Structure fields to DumpBlocks
!
! Revision 2.445  2006/07/07 23:08:32  pwagner
! Fixed bug in filling from GEOS5-derived grid
!
! Revision 2.444  2006/06/15 20:39:59  vsnyder
! Add PFA oversampling
!
! Revision 2.443  2006/06/15 00:01:42  pwagner
! Should work with geos5: convert then concatenate
!
! Revision 2.442  2006/06/13 22:13:12  pwagner
! changed interface to ConvertFromEtaLevelGrids
!
! Revision 2.441  2006/06/08 23:54:51  vsnyder
! Add switches field to retrieve
!
! Revision 2.440  2006/06/06 21:55:48  pwagner
! May specify geos5 apriori files instead of dao (geos4)
!
! Revision 2.439  2006/06/03 01:43:10  vsnyder
! Allow duplicate fields on destroy command
!
! Revision 2.438  2006/06/01 03:06:18  vsnyder
! Define numGrad and numNewt
!
! Revision 2.437  2006/05/27 02:59:08  vsnyder
! Add a 'clean' option to the 'dump' command
!
! Revision 2.436  2006/05/11 19:37:32  pwagner
! Added option to disallow duplicate molecules
!
! Revision 2.435  2006/05/09 16:39:58  pwagner
! Added writing l2cf to dgm
!
! Revision 2.434  2006/05/04 23:04:59  pwagner
! May convertEtaToP and create a VGrid in MergeGrids section
!
! Revision 2.433  2006/04/21 22:29:38  vsnyder
! Allow FlushPFA in global settings section, allow molecules field on FlushPFA
!
! Revision 2.432  2006/04/03 23:51:21  livesey
! Added f_surface to subset
!
! Revision 2.431  2006/03/23 01:51:32  vsnyder
! Allow only [UL]SBPFAMolecules and [UL]SBLBLMolecules fields to be empty
!
! Revision 2.430  2006/03/22 02:23:26  vsnyder
! Add lsbLBLmolecules, useLBLmolecules, lsGlobal, lsLocal, lsWeighted
!
! Revision 2.429  2006/03/17 00:06:31  pwagner
! Change default to allowing overlaps outside processingRange
!
! Revision 2.428  2006/03/15 23:52:24  pwagner
! Removed InputVersion component from PCF, l2cf
!
! Revision 2.427  2006/03/13 23:41:53  pwagner
! Added c=numeric type field to Fill via manipulation
!
! Revision 2.426  2006/03/07 00:50:59  pwagner
! May change already-set Booleans via reevaluate command
!
! Revision 2.425  2006/03/04 00:15:30  pwagner
! Added things for runtime Booleans to selectively skip phases
!
! Revision 2.424  2006/02/16 00:11:06  pwagner
! Added stamp boolean field to phase asks for printing phase names, times
!
! Revision 2.423  2006/02/10 21:11:29  pwagner
! May specify skipRetrivel for particular Phases
!
! Revision 2.422  2006/02/03 21:26:47  pwagner
! ChunkDivide adds allowPostOverlaps field
!
! Revision 2.421  2006/01/26 00:36:20  pwagner
! DU synonym for DobsonUnits
!
! Revision 2.420  2006/01/11 17:04:32  pwagner
! May specify unit when filling column abundances
!
! Revision 2.419  2006/01/06 01:16:34  pwagner
! silent boolean field can silence selected phases
!
! Revision 2.418  2005/12/29 01:11:08  vsnyder
! Add boolean 'refract' field to ForwardModel spec
!
! Revision 2.417  2005/12/21 21:47:45  livesey
! Added negateSD and precisionFactor arguments to Retrieve
!
! Revision 2.416  2005/11/15 00:19:24  pwagner
! NAG intolerant of too many continuation lines
!
! Revision 2.415  2005/11/11 21:44:18  pwagner
! Added avoidBrightObjects for use by FillFromL1B; removed unused globalsettings
!
! Revision 2.414  2005/09/23 23:39:16  pwagner
! Added rename field to copy command
!
! Revision 2.413  2005/09/19 16:55:06  pwagner
! Save Obstructions will allow OutputClose to faithfully replicate data gaps
!
! Revision 2.412  2005/09/16 23:38:04  vsnyder
! Add spect_der field to ForwardModel
!
! Revision 2.411  2005/09/09 23:07:26  pwagner
! Changes to allow pre-starttime overlaps in first chunk
!
! Revision 2.410  2005/08/19 23:24:55  pwagner
! Allows HGrid, Copy commands in Output section
!
! Revision 2.409  2005/08/04 19:36:18  pwagner
! New copy command in Output section
!
! Revision 2.408  2005/08/04 02:58:13  vsnyder
! Add MIFDeadTime to QuantityTypes
!
! Revision 2.407  2005/08/03 18:06:46  vsnyder
! Scan averaging, some spectroscopy derivative stuff and cannonball polishing
!
! Revision 2.406  2005/07/21 23:43:43  pwagner
! Added extras for explicit fill fields where fillvalued
!
! Revision 2.405  2005/07/12 17:40:24  pwagner
! May fill status with condition that no gmaos found
!
! Revision 2.404  2005/06/30 22:43:08  livesey
! Added residualSupplied
!
! Revision 2.403  2005/06/24 23:41:47  vsnyder
! Add LineCenter, LineWidth and LineWidth_TDep to T_QuantityType
!
! Revision 2.402  2005/06/03 02:11:52  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! add PFAStru and PFANum fields to Dump command.
!
! Revision 2.401  2005/05/27 23:57:03  vsnyder
! Add Flush PFAData
!
! Revision 2.400  2005/05/26 22:35:48  vsnyder
! Add PFAFiles field to ForwardModelGlobal
!
! Revision 2.399  2005/05/02 23:00:40  vsnyder
! Add PFAFile parameter
!
! Revision 2.398  2005/04/19 19:14:09  livesey
! Changed some section ordering issues related to pfa generation
!
! Revision 2.397  2005/04/01 20:48:28  vsnyder
! Add mark and text fields to dump command
!
! Revision 2.396  2005/03/16 23:59:42  vsnyder
! Add allLinesForRadiometer and allLinesInCatalog to makePFA
!
! Revision 2.395  2005/03/15 01:27:31  vsnyder
! Allow Dump command in Retrieve section
!
! Revision 2.394  2005/01/27 21:09:50  vsnyder
! Delete "file" field from PFAData.  Add [LU]SBPFAMolecules to ForwardModel and
! delete PFAMolecules.  Some cannonball polishing.
!
! Revision 2.393  2005/01/07 01:04:09  vsnyder
! Remove f_AllPFA from ReadPfa
!
! Revision 2.392  2004/12/31 02:41:56  vsnyder
! Working on read/write PFA database
!
! Revision 2.391  2004/12/13 20:23:00  vsnyder
! Added MakePFA and WritePFA commands.  Put S_Sids into alphabetical order.
! Added AllLines, AllSignals, AllSpectra, Lines, Signals, Spectroscopy and
! Stop fields to Dump command.
!
! Revision 2.390  2004/11/19 20:52:50  livesey
! File no longer required in direct write (how come we didn't have problem
! with this before?)
!
! Revision 2.389  2004/11/08 21:56:43  livesey
! Added badRange
!
! Revision 2.388  2004/10/30 00:26:47  vsnyder
! Add 'spectroscopy' field to DumpCommand
!
! Revision 2.387  2004/10/21 00:45:02  vsnyder
! Require mlsSignals and spectroscopy before globalSettings etc.
!
! Revision 2.386  2004/10/16 17:25:12  livesey
! Added asciiFile fill method
!
! Revision 2.385  2004/10/13 02:23:35  livesey
! Added geodAlt to Forge and geocAltitude as a vGrid coordinate
!
! Revision 2.384  2004/09/27 20:10:51  livesey
! Added applyBaseline fill method and it's supporting stuff
!
! Revision 2.383  2004/09/25 00:15:22  livesey
! Added combineChannels to Algebra
!
! Revision 2.382  2004/09/10 23:53:36  livesey
! Added centervertically option for bin mean/max/min fill.
!
! Revision 2.381  2004/08/09 21:43:21  livesey
! Added maxOrbY argument to ChunkDivide
!
! Revision 2.380  2004/07/22 20:39:14  cvuu
! Now can fill ForwardModel time, mean and std_dev
!
! Revision 2.379  2004/07/17 02:28:19  vsnyder
! Add dump for entire PFA database
!
! Revision 2.378  2004/06/29 00:08:41  pwagner
! Now can fill timings
!
! Revision 2.377  2004/06/17 22:48:48  pwagner
! integer now a possible coord type for a VGrid
!
! Revision 2.376  2004/06/16 23:42:44  vsnyder
! Make molecules field of pfaData optional
!
! Revision 2.375  2004/06/16 01:23:28  vsnyder
! Make VelLin field of PFAData optional
!
! Revision 2.374  2004/06/12 00:40:48  vsnyder
! Allow formula field in tGrid, don't require step and number
!
! Revision 2.373  2004/06/08 19:26:20  vsnyder
! Add tGrid
!
! Revision 2.372  2004/05/29 02:47:02  vsnyder
! Rearrange function definition stuff
!
! Revision 2.371  2004/05/22 02:29:30  vsnyder
! Add PFAdata, more stuff for dump, allow dump in global_settings
!
! Revision 2.370  2004/05/18 01:06:28  vsnyder
! Add HGrid field to DUMP command
!
! Revision 2.369  2004/05/04 01:03:33  livesey
! Added excludebelowbottom flag for binmax/binmin fill
!
! Revision 2.368  2004/05/01 04:03:00  vsnyder
! Added pfaMolecules, added 'dump' in Construct
!
! Revision 2.367  2004/04/30 21:49:17  livesey
! Added DisjointEquations command
!
! Revision 2.366  2004/04/29 01:26:39  livesey
! More algebra stuff
!
! Revision 2.365  2004/04/28 23:07:34  livesey
! More stuff for algebra
!
! Revision 2.364  2004/04/16 00:48:40  livesey
! Added singleChannelRadiance type and extractChannel fill method
!
! Revision 2.363  2004/03/31 03:59:43  livesey
! Added singleMAF option
!
! Revision 2.362  2004/03/24 01:03:12  livesey
! Added f_date to s_hgrid
!
! Revision 2.361  2004/03/22 18:26:13  livesey
! Added allLinesInCatalog and combineChannels
!
! Revision 2.360  2004/03/17 17:15:45  livesey
! New quantity types, new fill methods
!
! Revision 2.359  2004/03/10 22:19:43  livesey
! Added quality fill method and resolution argument to explicit vGrids
!
! Revision 2.358  2004/02/17 14:07:54  livesey
! Added functionality to the boxCar and binned fills
!
! Revision 2.357  2004/02/11 23:12:23  livesey
! Typo, forgot quality, had 2 statuses
!
! Revision 2.356  2004/02/11 17:21:51  pwagner
! May DirectWrite l2gp status and quality quantities
!
! Revision 2.355  2004/02/10 21:16:58  livesey
! Added status and quality as valid quantity types
!
! Revision 2.354  2004/02/06 01:01:23  livesey
! Added boxcar fill
!
! Revision 2.353  2004/01/30 23:25:09  livesey
! Added columnScale and rowScale commands for algebra
!
! Revision 2.352  2004/01/29 03:33:16  livesey
! Added reflect and cyclicJacobi for z_algebra
!
! Revision 2.351  2004/01/24 01:04:38  livesey
! Added the adopted quantity type.
!
! Revision 2.350  2004/01/23 19:09:00  livesey
! More work on adoption / loading
!
! Revision 2.349  2004/01/23 05:47:38  livesey
! Added the adoption stuff
!
! Revision 2.348  2004/01/20 20:26:12  livesey
! Added binMean
!
! Revision 2.347  2004/01/17 03:04:15  vsnyder
! Provide for functions in expressions
!
! Revision 2.346  2004/01/17 00:28:09  vsnyder
! Provide for Algebra section
!
! Revision 2.345  2004/01/14 18:49:58  vsnyder
! Stuff to support the Algebra section
!
! Revision 2.344  2003/12/11 22:59:08  pwagner
! May fill DirectWriteDatabase in global settings
!
! Revision 2.343  2003/12/04 22:19:32  livesey
! Added ability to fill from l2gpPrecision field
!
! Revision 2.342  2003/11/15 00:46:41  pwagner
! maxfailurespermachine, maxfailuresperchunk no longer configuration settings (see comline opts)
!
! Revision 2.341  2003/11/05 01:04:30  pwagner
! May dump either entire vector or single quantity in Fill
!
! Revision 2.340  2003/11/01 18:45:00  livesey
! Made f_cost an optional argument to binSelector
!
! Revision 2.339  2003/10/22 21:17:05  pwagner
! aPhaseName: Phase added to Fill, Construct sections to time phases
!
! Revision 2.338  2003/10/16 23:45:03  pwagner
! Completed maxfailures per machine, chunk
!
! Revision 2.337  2003/10/15 23:12:16  livesey
! Added resetUnusedRadiances
!
! Revision 2.336  2003/10/10 23:28:19  vsnyder
! Allow rowQuantity and colQuantity to be optional in dumpBlocks
!
! Revision 2.335  2003/10/07 01:14:13  vsnyder
! Add noAbsent field to dumpBlocks
!
! Revision 2.334  2003/10/06 13:16:09  cvuu
! add new description=strat to handle reading the ncep data file
!
! Revision 2.333  2003/09/11 23:15:22  livesey
! Added xStar and yStar to forward model config
!
! Revision 2.332  2003/09/02 18:03:23  pwagner
! Now can reset maxfailuresper chunk, machine from global settings
!
! Revision 2.331  2003/08/16 01:18:55  livesey
! Added baseline forward model in its own right.
!
! Revision 2.330  2003/08/16 00:30:25  vsnyder
! Add magAzEl to fill method type
!
! Revision 2.329  2003/08/14 20:11:31  pwagner
! DirectWrite may take l2fwm types for fwm radiances
!
! Revision 2.328  2003/08/13 00:49:04  livesey
! Added the polarLinear forward model and the field based bin selectors
!
! Revision 2.327  2003/08/11 18:08:00  livesey
! Added the single option to hGrid
!
! Revision 2.326  2003/08/08 23:06:21  livesey
! Added the fieldStrength etc. stuff, also dontPack option on saving l2pc
! files.
!
! Revision 2.325  2003/07/16 01:06:36  vsnyder
! Add DACS filter shapes
!
! Revision 2.324  2003/07/15 22:11:37  livesey
! Added stuff for hybrid model
!
! Revision 2.323  2003/07/08 00:14:37  livesey
! New version of directWrite command
!
! Revision 2.322  2003/06/24 23:30:50  livesey
! Tidyups in directWrite and label
!
! Revision 2.321  2003/06/23 18:06:34  pwagner
! Should allow us to write metadata after DirectWrite
!
! Revision 2.320  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.319  2003/06/06 01:06:44  livesey
! Added delete option for grids
!
! Revision 2.318  2003/06/05 22:13:35  pwagner
! Added criticalSignals to ChunkDivide Orbital
!
! Revision 2.317  2003/06/04 01:09:10  livesey
! Added flushL2PCBins
!
! Revision 2.316  2003/06/03 19:24:33  livesey
! Added flushL2PCBins
!
! Revision 2.315  2003/05/29 20:02:22  livesey
! Added reflector temperature model and supporting stuff
!
! Revision 2.314  2003/05/29 16:43:54  livesey
! Renamed sideband fraction stuff and added some reflector orientated
! stuff
!
! Revision 2.313  2003/05/26 06:33:12  livesey
! Removed duplicate quantity template
!
! Revision 2.312  2003/05/20 23:30:12  dwu
! add iwcfromextinction
!
! Revision 2.311  2003/05/19 20:14:59  vsnyder
! Remove private declarations for module entities accessed without using
! an "only" clause.  These cause NAG to issue a message that the name is
! explicitly accessed from a module but not used.  Bizarre.
!
! Revision 2.310  2003/05/12 20:57:21  pwagner
! Added L2ParSF spec to allow changing staging file name in global settings
!
! Revision 2.309  2003/05/10 01:05:51  livesey
! Added binTotal as allowed vector quantity
!
! Revision 2.308  2003/05/07 22:50:29  pwagner
! Optionally may skipL1BCheck in ChunkDivide
!
! Revision 2.307  2003/05/07 01:01:50  livesey
! Added chiSqBinned
!
! Revision 2.306  2003/05/05 23:00:34  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.305  2003/04/30 00:09:12  vsnyder
! Add the FrequencyGrid specification, with 'atmos' and 'frequencies'
! fields.  So far, the FrequencyGrid specification can't appear anywhere.
! Eventually, we may hook it into the Retrieve section.
!
! Revision 2.304  2003/04/23 17:06:36  livesey
! Added binmax binmin fills
!
! Revision 2.303  2003/04/22 00:17:10  dwu
! add a new option (clear_110RH_below_tropopause) to i_saturation
!
! Revision 2.302  2003/04/17 23:08:29  pwagner
! Added optional AuraInstrument field to l2gp apriori reads
!
! Revision 2.301  2003/04/11 23:14:31  livesey
! Added spreadChannel method and force option for fill.
!
! Revision 2.300  2003/04/11 21:56:51  livesey
! Added wmo tropopause
!
! Revision 2.299  2003/04/10 20:24:05  dwu
! add l_none in cloud_der
!
! Revision 2.298  2003/04/09 00:10:30  livesey
! New t_i_saturation and t_cloud_der for Dong and Jonathan
!
! Revision 2.297  2003/04/04 23:54:04  livesey
! Added split sideband fill
!
! Revision 2.296  2003/04/04 22:01:10  livesey
! Added updateMask stuff
!
! Revision 2.295  2003/04/04 00:11:29  livesey
! Added concatenate stuff
!
! Revision 2.294  2003/04/02 21:49:00  jonathan
! remove cloud_fov
!
! Revision 2.293  2003/03/27 20:45:34  livesey
! Added logSpace to fill
!
! Revision 2.292  2003/03/26 21:23:39  livesey
! Added scaleOvelaps to fill
!
! Revision 2.291  2003/03/07 03:18:13  livesey
! Added RestrictRange command
!
! Revision 2.290  2003/03/06 00:46:18  livesey
! Added subset and flagcloud to fill and minor tidyups.
!
! Revision 2.289  2003/03/05 19:11:22  livesey
! Added allow missing to fill.
!
! Revision 2.288  2003/03/01 00:24:27  pwagner
! Added missingValue as filed to reading Gridded data
!
! Revision 2.287  2003/02/27 17:58:34  bill
! Added polarized
!
! Revision 2.286  2003/02/20 21:26:21  pwagner
! Lets you read field dimList=x,y,.. w/ griddeddata
!
! Revision 2.285  2003/02/18 23:58:53  livesey
! Added phiWindow to fill for hydrostatic ptan
!
! Revision 2.284  2003/02/14 01:56:36  livesey
! Added the 'additional' capability in subset
!
! Revision 2.283  2003/02/13 21:43:50  livesey
! Added f_profile to fill
!
! Revision 2.282.2.5  2003/04/08 23:41:17  jonathan
! remove cloud_fov
!
! Revision 2.282.2.4  2003/03/27 00:49:25  vsnyder
! Add ECRtoFOV
!
! Revision 2.282.2.3  2003/03/13 00:06:23  vsnyder
! Delete moleculesPol and moleculeDerivativesPol
!
! Revision 2.282.2.2  2003/02/27 17:57:41  bill
! added everything from the main branch
!
! Revision 2.286  2003/02/20 21:26:21  pwagner
! Lets you read field dimList=x,y,.. w/ griddeddata
!
! Revision 2.285  2003/02/18 23:58:53  livesey
! Added phiWindow to fill for hydrostatic ptan
!
! Revision 2.284  2003/02/14 01:56:36  livesey
! Added the 'additional' capability in subset
!
! Revision 2.283  2003/02/13 21:43:50  livesey
! Added f_profile to fill
!
! Revision 2.282  2003/02/12 02:11:22  livesey
! Added extendedAverage
!
! Revision 2.281  2003/02/06 23:30:36  livesey
! New approach for Forge and explicit hGrids
!
! Revision 2.280  2003/02/06 22:04:48  vsnyder
! Add f_moleculesPol, f_moleculeDerivativesPol, delete f_polarized
!
! Revision 2.279  2003/02/06 00:45:20  livesey
! Added sza and nameFragment to binSelectors type
!
! Revision 2.278  2003/02/05 21:56:27  livesey
! binSelectors don't contain signals, instead forward model configs
! contain bin selectors.
!
! Revision 2.277  2003/02/05 04:06:18  dwu
! add cloudheight
!
! Revision 2.276  2003/01/30 17:28:35  jonathan
! add logical incl_cld
!
! Revision 2.275  2003/01/29 01:48:29  vsnyder
! Add 'polarized' field to forwardModel
!
! Revision 2.274  2003/01/18 02:37:21  livesey
! Added quantityType to l2aux
!
! Revision 2.273  2003/01/16 00:55:52  jonathan
! add do_1d
!
! Revision 2.272  2003/01/14 22:14:54  dwu
! make FlagCloud depend on both channels and cloudChannels
!
! Revision 2.271  2003/01/13 17:17:18  jonathan
!  change cloud_width to i_saturation
!
! Revision 2.270  2003/01/11 01:23:12  livesey
! Another bug fix in Dong's flagCloud
!
! Revision 2.269  2003/01/11 01:04:21  livesey
! Bug fix for Dong's stuff
!
! Revision 2.268  2003/01/11 01:02:27  livesey
! Added max and min value in subset
!
! Revision 2.267  2003/01/11 00:01:19  dwu
! add flagCloud
!
! Revision 2.266  2003/01/08 23:52:01  livesey
! Added irregular and sparseQuantities
!
! Revision 2.265  2003/01/07 23:58:33  livesey
! Typo!
!
! Revision 2.264  2003/01/07 23:58:15  livesey
! Bug fix
!
! Revision 2.263  2003/01/07 23:46:26  livesey
! Added reset for subset and magnetic model stuff
!
! Revision 2.262  2003/01/06 20:13:30  livesey
! New overlap handling in ChunkDivide and HGrid
!
! Revision 2.261  2002/11/27 19:26:20  livesey
! Added stuff for manipulate fill
!
! Revision 2.260  2002/11/21 01:17:54  livesey
! Added the negativePrecision command to fill (distinct from the
! negativePrecision option to the fill command)
!
! Revision 2.259  2002/11/15 01:33:31  livesey
! Added allLinesForRadiometer
!
! Revision 2.258  2002/10/25 23:55:50  livesey
! Added offsetAmount
!
! Revision 2.257  2002/10/25 22:25:14  livesey
! Added the dnwt chisquared vector types
!
! Revision 2.256  2002/10/23 04:11:25  livesey
! Fixed up a 'too many continuation lines' problem.
!
! Revision 2.255  2002/10/23 01:32:31  vsnyder
! Add CovSansReg switch to retrieve spec
!
! Revision 2.254  2002/10/19 23:41:39  livesey
! Added muMin
!
! Revision 2.253  2002/10/19 01:51:52  livesey
! Added serial option to retrieve
!
! Revision 2.252  2002/10/17 18:27:27  livesey
! Put bounds in wrong place
!
! Revision 2.251  2002/10/17 18:18:38  livesey
! Added low/high bound options to vector definitions.
!
! Revision 2.250  2002/10/17 00:16:40  vsnyder
! Add lowBound and highBound fields for the Retrieve spec
!
! Revision 2.249  2002/10/16 20:18:21  mjf
! Added GPH precision.
!
! Revision 2.248  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.247  2002/10/03 13:44:53  mjf
! Renamed temperaturePrecisionQuantitiy to tempPrecisionQuantity to get
! names < 31 long.
!
! Revision 2.246  2002/10/02 23:20:54  livesey
! Changed regQuants to quantity templates rather than types
!
! Revision 2.245  2002/10/02 23:04:48  pwagner
! RHI now separate fill methods
!
! Revision 2.244  2002/10/02 02:43:07  livesey
! Added regApriori
!
! Revision 2.243  2002/10/01 18:29:40  mjf
! Added new quantitites for new Fill for RHi precision including T
! error.
!
! Revision 2.242  2002/09/25 20:07:19  livesey
! Added specificQuantities to s_forwardModel.  Allowed s_forwardModel
! inside z_construct
!
! Revision 2.241  2002/09/24 21:38:03  livesey
! Added minValue
!
! Revision 2.240  2002/09/23 22:15:05  vsnyder
! Delete maxF field name and numF literal name
!
! Revision 2.239  2002/09/05 20:58:20  livesey
! Another fix to Transfer!
!
! Revision 2.238  2002/09/05 20:50:24  livesey
! Fixed requirements for transfer
!
! Revision 2.237  2002/09/05 20:49:30  livesey
! Added skipMask to transfer
!
! Revision 2.236  2002/08/28 01:14:05  livesey
! Added offsetRadiance fill
!
! Revision 2.235  2002/08/26 20:01:36  livesey
! Added instances argument to fill (profile method)
!
! Revision 2.234  2002/08/24 01:38:28  vsnyder
! Implement horizontal regularization
!
! Revision 2.233  2002/08/15 03:52:40  livesey
! Added stuff for profile fill
!
! Revision 2.232  2002/08/08 22:04:16  vsnyder
! Add tikhonov to mask type
!
! Revision 2.231  2002/08/04 16:07:57  mjf
! New method PE for ChunkDivide.
!
! Revision 2.230  2002/08/03 01:14:42  vsnyder
! Add regAfter to control Tikhonov before/after column scaling -- default false
!
! Revision 2.229  2002/07/18 22:02:56  vsnyder
! Exploit new CONTINUE argument of MAKE_TREE to get rid of ACORN
!
! Revision 2.228  2002/07/17 06:02:36  livesey
! New HDF5 l2pc stuff
!
! Revision 2.227  2002/07/05 17:47:19  livesey
! Allow join to be within the chunk processing
!
! Revision 2.226  2002/06/26 01:26:48  livesey
! Added 2d pressure guessing stuff
!
! Revision 2.225  2002/06/24 18:27:17  livesey
! Added 2D scan model
!
! Revision 2.224  2002/06/12 17:57:04  livesey
! Added ascii field for l2pc
!
! Revision 2.223  2002/06/04 22:48:36  livesey
! Added refract fill for phiTan
!
! Revision 2.222  2002/06/04 22:06:36  livesey
! Added phitan as a state vector element
!
! Revision 2.221  2002/05/22 00:48:10  livesey
! Added direct write
!
! Revision 2.220  2002/05/17 17:56:39  livesey
! Stuff for sideband folding fills
!
! Revision 2.219  2002/05/14 00:27:22  livesey
! Added system temperature and noise bandwidth quantity types
!
! Revision 2.218  2002/05/07 20:25:59  livesey
! Added writeCounterMAF option for l2aux
!
! Revision 2.217  2002/05/07 01:01:25  vsnyder
! Change regWeight to regWeights
!
! Revision 2.216  2002/05/06 21:37:50  livesey
! Added forbidOverspill to hGrid
!
! Revision 2.215  2002/05/01 22:02:09  pwagner
! Can again read leapsecfile in global_settings
!
! Revision 2.214  2002/05/01 00:24:26  pwagner
! Undid changes to allow leapsecfile to be read; they caused crashes
!
! Revision 2.213  2002/04/29 16:38:42  pwagner
! Can specifiy leapsecfile in global settings
!
! Revision 2.212  2002/04/22 20:55:00  vsnyder
! Compute and output the averaging kernel
!
! Revision 2.211  2002/04/18 18:40:17  pwagner
! Can fill h2o from an rhiQuantity
!
! Revision 2.210  2002/04/10 17:44:50  pwagner
! Added rhi quantity (but is this enough?)
!
! Revision 2.209  2002/04/04 16:32:42  livesey
! Added negative error bar stuff
!
! Revision 2.208  2002/03/15 21:22:20  livesey
! Introduced new binselector type
!
! Revision 2.207  2002/03/14 00:39:20  pwagner
! Can fill scVelECI and scVelECR from l1b now
!
! Revision 2.206  2002/03/13 22:02:31  livesey
! Changed from explicitFill to fill
!
! Revision 2.205  2002/03/08 08:07:00  livesey
! Added explicit fill mask
!
! Revision 2.204  2002/03/07 17:18:03  livesey
! Removed frqGap
!
! Revision 2.203  2002/02/22 01:07:58  pwagner
! Added metaName
!
! Revision 2.202  2002/02/20 02:12:59  livesey
! Added height to Fill
!
! Revision 2.201  2002/02/13 00:08:40  livesey
! Added differential scan model
!
! Revision 2.200  2002/02/09 19:12:00  livesey
! Added optical depth stuff
!
! Revision 2.199  2002/02/07 21:40:26  livesey
! Added string for t_masks
!
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
