! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  use INTRINSIC ! Everything. FIRST_LIT, INIT_INTRINSIC,
    ! L_FALSE, L_TRUE, LAST_INTRINSIC_LIT, T_BOOLEAN, T_FIRST,
    ! T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE and T_STRING are used
    ! here, but everything is included so that it can be gotten by
    ! USE INIT_TABLES_MODULE.
  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER
  use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
  use TREE_TYPES, only: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
                        N_NAME_DEF, N_SECTION, N_SPEC_DEF

  implicit NONE
  public ! This would be a MUCH LONGER list than the list of private
  !        names below.
  private :: ADD_IDENT, BUILD_TREE, ENTER_TERMINAL, INIT_INTRINSIC
  private :: MAKE_TREE, N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE
  private :: N_NAME_DEF, N_SECTION, N_SPEC_DEF, PUSH_PSEUDO_TERMINAL
  private :: T_IDENTIFIER

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Enumeration types:
  integer, public, parameter :: T_METHOD         = t_last_intrinsic+1
  integer, public, parameter :: T_PMODE          = t_method+1
  integer, public, parameter :: T_UNITS          = t_pmode+1
  integer, public, parameter :: T_LAST           = t_units
  integer, public :: DATA_TYPE_INDICES(t_first:t_last)
! Field indices:
  integer, public, parameter :: F_BYPASS = 1
  integer, public, parameter :: F_DF = 2
  integer, public, parameter :: F_DLAT = 3
  integer, public, parameter :: F_DLON = 4
  integer, public, parameter :: F_FILE = 5
  integer, public, parameter :: F_IMETHOD = 6
  integer, public, parameter :: F_LABEL = 7
  integer, public, parameter :: F_LATGRID = 8
  integer, public, parameter :: F_LONGRID = 9
  integer, public, parameter :: F_MCF = 10
  integer, public, parameter :: F_MODE = 11
  integer, public, parameter :: F_PRESLVL = 12
  integer, public, parameter :: F_PRODNAME = 13
  integer, public, parameter :: F_QUANTITIES = 14
  integer, public, parameter :: F_RANGFREQ = 15
  integer, public, parameter :: F_RANGWAVNUM = 16
  integer, public, parameter :: F_TIME = 17
  integer, public, parameter :: FIELD_FIRST = f_bypass, FIELD_LAST = f_time
  integer, public :: FIELD_INDICES(field_first:field_last)
! Enumeration literals:
  integer, public, parameter :: L_ALL   = last_intrinsic_lit + 1
  integer, public, parameter :: L_ASC   = l_all + 1
  integer, public, parameter :: L_COM 	= l_asc + 1
  integer, public, parameter :: L_CSP   = l_com + 1
  integer, public, parameter :: L_DES   = l_csp + 1
  integer, public, parameter :: L_LIN   = l_des + 1
  integer, public, parameter :: LAST_LIT = l_lin
  integer, public :: LIT_INDICES(first_lit:last_lit)
! Parameter names:
  ! In GlobalSettings section:
  integer, public, parameter :: P_L2_VER = 1
  integer, public, parameter :: P_MAX_GAP = 2
  integer, public, parameter :: P_N = 3
  integer, public, parameter :: P_OUTPUT_VERSION_STRING = 4
  integer, public, parameter :: P_VERSION_COMMENT = 5
  integer, public, parameter :: FIRST_PARM = P_L2_VER
  integer, public, parameter :: LAST_PARM = P_VERSION_COMMENT
  integer, public :: PARM_INDICES(first_parm:last_parm)
! Section identities:
  integer, public, parameter :: Z_DAILYMAP = 2
  integer, public, parameter :: Z_GLOBALSETTINGS = 1
  integer, public, parameter :: Z_OUTPUT = 3
  integer, public, parameter :: SECTION_FIRST = z_globalSettings, &
                                SECTION_LAST = z_Output
  integer, public :: SECTION_INDICES(section_first:section_last)
! Specification indices:
  integer, public, parameter :: S_MAPSPEC = 1
  integer, public, parameter :: S_OUTPUT = 2
  integer, public, parameter :: SPEC_FIRST = s_mapSpec, SPEC_LAST = s_Output
  integer, public :: SPEC_INDICES(spec_first:spec_last)

! Table for section ordering:
  integer, public, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: | globalSettings       |
!     |       dailyMap       |
!     |             output   |
! ====|==============================|== From: ==
        (/OK,    0,    0,  & ! Start
           0,   OK,    0,  & ! GlobalSettings
           0,    0,   OK,  & ! DailyMap
           0,    0,    0/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

  integer, private, parameter :: F = 1000, L = 2000, N = 0, P = 3000
  integer, private, parameter :: S = 4000, T = 5000, Z = 6000
  integer, private, parameter :: BEGIN = -1

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES
  ! Put intrinsic predefined identifiers into the symbol table.
    call init_intrinsic ( data_type_indices, lit_indices )

  ! Put nonintrinsic predefined identifiers into the symbol table.
    ! Put enumeration type names into the symbol table
    data_type_indices(t_method) =           add_ident ( 'interpolationMethod' )
    data_type_indices(t_pmode) =            add_ident ( 'processingMode' )
    data_type_indices(t_units) =            add_ident ( 'units' )
    ! Put enumeration literals into the symbol table:
    lit_indices(l_all) =                    add_ident ( 'all' )
    lit_indices(l_asc) =                    add_ident ( 'asc' )
    lit_indices(l_com) =                    add_ident ( 'com' )
    lit_indices(l_csp) =                    add_ident ( 'csp' )
    lit_indices(l_des) =                    add_ident ( 'des' )
    lit_indices(l_lin) =                    add_ident ( 'lin' )
    ! Put field names into the symbol table
    field_indices(f_bypass) =               add_ident ( 'bypassPCF' )
    field_indices(f_df) =                   add_ident ( 'dF' )
    field_indices(f_dlat) =                 add_ident ( 'dLat' )
    field_indices(f_dlon) =                 add_ident ( 'dLon' )
    field_indices(f_file) =                 add_ident ( 'file' )
    field_indices(f_imethod) =              add_ident ( 'intpMethod' )
    field_indices(f_label) =                add_ident ( 'label' )
    field_indices(f_latgrid) =              add_ident ( 'latGridMap' )
    field_indices(f_longrid) =              add_ident ( 'longGrid' )
    field_indices(f_mcf) =                  add_ident ( 'mcf' )
    field_indices(f_mode) =                 add_ident ( 'mode' )
    field_indices(f_preslvl) =              add_ident ( 'l3presLvl' )
    field_indices(f_prodname) =             add_ident ( 'l3prodNameD' )
    field_indices(f_quantities) =           add_ident ( 'quantities' )
    field_indices(f_rangfreq) =             add_ident ( 'rangFrequency' )
    field_indices(f_rangwavnum) =           add_ident ( 'rangWavenumber' )
    field_indices(f_time) =                 add_ident ( 'timeD' )
    ! Put parameter names into the symbol table
    parm_indices(p_l2_ver) =                add_ident ( 'L2Ver' )
    parm_indices(p_max_gap) =               add_ident ( 'MaxGap' )
    parm_indices(p_n)=                      add_ident ( 'N' )
    parm_indices(p_output_version_string) = add_ident ( 'OutputVersionString' )
    parm_indices(p_version_comment) =       add_ident ( 'VersionComment' )
    ! Put section names into the symbol table
    section_indices(z_dailymap) =           add_ident ( 'DailyMap' )
    section_indices(z_globalsettings) =     add_ident ( 'GlobalSettings' )
    section_indices(z_output) =             add_ident ( 'Output' )
    ! Put spec names into the symbol table
    spec_indices(s_mapspec) =               add_ident ( 'mapSpec' )
    spec_indices(s_output) =                add_ident ( 'output' )

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
      begin, t+t_boolean, l+l_true, l+l_false, n+n_dt_def, &
      begin, t+t_method, l+l_csp, l+l_lin, n+n_dt_def, &
      begin, t+t_pmode, l+l_all, l+l_asc, l+l_com, l+l_des, n+n_dt_def, &
      begin, t+t_units, l+l_days, l+l_deg, l+l_degrees, &
             l+l_dimensionless, l+l_dimless, l+l_dl, l+l_ghz, &
             l+l_hours, l+l_hpa, l+l_hz, l+l_k, l+l_khz, l+l_km, l+l_logp, &
             l+l_m, l+l_maf, l+l_mafs, l+l_mb, l+l_meters, l+l_mhz, &
             l+l_mif, l+l_mifs, l+l_minutes, l+l_orbits, l+l_pa, l+l_ppbv, &
             l+l_ppmv, l+l_pptv, l+l_rad, l+l_radians, l+l_s, l+l_seconds, &
             l+l_thz, l+l_vmr, l+l_zeta, n+n_dt_def /) )
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
      begin, s+s_mapSpec, &
             begin, f+f_prodname, t+t_string, n+n_field_type, &
             begin, f+f_time, t+t_string, n+n_field_type, &
             begin, f+f_imethod, t+t_method, n+n_field_type, &
             begin, f+f_latgrid, t+t_numeric_range, n+n_field_type, &
             begin, f+f_dlat, t+t_numeric, n+n_field_type, &
             begin, f+f_longrid, t+t_numeric_range, n+n_field_type, &
             begin, f+f_dlon, t+t_numeric, n+n_field_type, &
             begin, f+f_rangfreq, t+t_numeric_range, n+n_field_type, &
             begin, f+f_df, t+t_numeric, n+n_field_type, &
             begin, f+f_preslvl, t+t_numeric_range, n+n_field_type, &
             begin, f+f_rangwavnum, t+t_numeric_range, n+n_field_type, &
             begin, f+f_mode, t+t_pmode, n+n_field_type, &
             begin, f+f_label, t+t_string, n+n_field_type, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_output, &
             begin, f+f_quantities, t+t_string, n+n_field_type, &
             begin, f+f_mcf, t+t_string, n+n_field_type, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             begin, f+f_bypass, t+t_boolean, n+n_field_type, &
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
             begin, p+p_l2_ver, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_max_gap, t+t_numeric, n+n_name_def, &
             begin, p+p_n, t+t_numeric, n+n_name_def, &
             n+n_section, &
      begin, z+z_dailymap, s+s_mapspec, n+n_section, &
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
