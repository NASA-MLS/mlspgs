! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  use Init_MLSSignals_m ! Everything.
  use INTRINSIC ! Everything. ADD_IDENT, FIRST_LIT, INIT_INTRINSIC,
    ! L_FALSE, L_TRUE, LAST_INTRINSIC_LIT, T_BOOLEAN, T_FIRST,
    ! T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE and T_STRING are used
    ! here, but everything is included so that it can be gotten by
    ! USE INIT_TABLES_MODULE.
  use Units, only: Init_Units

  implicit NONE
  public ! This would be a MUCH LONGER list than the list of private
  !        names below.
  private :: ADD_IDENT, MAKE_TREE

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Enumeration types:
  integer, public, parameter :: T_INTP           = Last_signal_type+1
  integer, public, parameter :: T_PMODE          = t_intp+1
  integer, public, parameter :: T_UNITS          = t_pmode+1
  integer, public, parameter :: T_LAST           = t_units
! Field indices:
  integer, public, parameter :: F_ALVL = last_Signal_Field + 1
  integer, public, parameter :: F_DLON = f_dlat + 1
  integer, public, parameter :: F_DLVL = f_dlon + 1
  integer, public, parameter :: F_FILE = f_dlvl + 1
  integer, public, parameter :: F_LABEL = f_file + 1
  integer, public, parameter :: F_LONGRID = f_label + 1
  integer, public, parameter :: F_MCF = f_longrid + 1
  integer, public, parameter :: F_MODE = f_mcf + 1
  integer, public, parameter :: F_PRESLVL = f_mode + 1
  integer, public, parameter :: F_PRODNAME = f_preslvl + 1
  integer, public, parameter :: F_RANGFREQ = f_prodname + 1
  integer, public, parameter :: F_RANGWAVNUM = f_rangfreq + 1
  integer, public, parameter :: F_TIME = f_rangwavnum + 1
  integer, public, parameter :: FIELD_LAST = f_time
! Enumeration literals:
  integer, public, parameter :: L_ADO   = last_signal_lit+1
  integer, public, parameter :: L_ALL   = l_ado + 1
  integer, public, parameter :: L_ASC   = l_all + 1
  integer, public, parameter :: L_COM 	= l_asc + 1
  integer, public, parameter :: L_CSP   = l_com + 1
  integer, public, parameter :: L_DES   = l_csp + 1
  integer, public, parameter :: L_LIN   = l_des + 1
  integer, public, parameter :: LAST_LIT = l_lin
! Section identities:
  integer, public, parameter :: Z_DIAGNOSTIC = 3
  integer, public, parameter :: Z_GLOBALSETTINGS = 1
  integer, public, parameter :: Z_OUTPUT = 4
  integer, public, parameter :: Z_STANDARD = 2
  integer, public, parameter :: SECTION_FIRST = z_globalSettings, &
                                SECTION_LAST = z_Output
! Specification indices:
  integer, public, parameter :: S_DG = last_Signal_Spec + 1
  integer, public, parameter :: S_OUT = s_dg + 1
  integer, public, parameter :: S_STD = s_out + 1
  integer, public, parameter :: SPEC_LAST = s_Std
! Parameter names:
  ! In GlobalSettings section:
  integer, public, parameter :: P_DLAT = spec_last + 1
  integer, public, parameter :: P_DTYPE = p_dlat + 1
  integer, public, parameter :: P_IMETHOD = p_dtype + 1
  integer, public, parameter :: P_L2_MAP_LATS = p_imethod + 1
  integer, public, parameter :: P_L2_NOM_LATS = p_l2_map_lats + 1
  integer, public, parameter :: P_L2DGTYPE = p_l2_nom_lats + 1
  integer, public, parameter :: P_MAX_GAP = p_l2dgtype + 1
  integer, public, parameter :: P_MIN_DAYS = p_max_gap + 1
  integer, public, parameter :: P_N = p_min_days + 1
  integer, public, parameter :: P_OUTPUT_VERSION_STRING = p_n + 1
  integer, public, parameter :: P_STYPE = p_output_version_string + 1
  integer, public, parameter :: FIRST_PARM = P_DLAT
  integer, public, parameter :: LAST_PARM = P_STYPE

! Table for section ordering:
  integer, public, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: | globalSettings             |
!     |       Standard             |
!     |           Diagnostic       |
!     |                    Output  |
! ====|============================|== From: ==
        (/OK,    0,    0,    0,  & ! Start
           0,   OK,    0,    0,  & ! GlobalSettings
           0,    0,   OK,    0,  & ! Standard
           0,    0,    0,   OK,  & ! Diagnostic
           0,    0,    0,    0/) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES
     use TREE_TYPES, only: N_DT_DEF, N_FIELD_TYPE, &
                           N_NAME_DEF, N_SECTION, N_SPEC_DEF
  ! Put intrinsic predefined identifiers into the symbol table.
    call init_MLSSignals ( t_last, field_last, last_lit, &
    & first_parm, last_parm, section_last, spec_last )
  ! Put nonintrinsic predefined identifiers into the symbol table.
    ! Put enumeration type names into the symbol table
    data_type_indices(t_intp) =             add_ident ( 'interpolationMethod' )
    data_type_indices(t_pmode) =            add_ident ( 'processingMode' )
    data_type_indices(t_units) =            add_ident ( 'units' )
    ! Put enumeration literals into the symbol table:
    lit_indices(l_ado) =                    add_ident ( 'ado' )
    lit_indices(l_all) =                    add_ident ( 'all' )
    lit_indices(l_asc) =                    add_ident ( 'asc' )
    lit_indices(l_com) =                    add_ident ( 'com' )
    lit_indices(l_csp) =                    add_ident ( 'csp' )
    lit_indices(l_des) =                    add_ident ( 'des' )
    lit_indices(l_lin) =                    add_ident ( 'lin' )
    ! Put field names into the symbol table
    field_indices(f_alvl) =                 add_ident ( 'ascPresLvl' )
    field_indices(f_dlon) =                 add_ident ( 'dLon' )
    field_indices(f_dlvl) =                 add_ident ( 'desPresLvl' )
    field_indices(f_file) =                 add_ident ( 'file' )
    field_indices(f_label) =                add_ident ( 'label' )
    field_indices(f_longrid) =              add_ident ( 'longGrid' )
    field_indices(f_mcf) =                  add_ident ( 'mcf' )
    field_indices(f_mode) =                 add_ident ( 'mode' )
    field_indices(f_preslvl) =              add_ident ( 'l3presLvl' )
    field_indices(f_prodname) =             add_ident ( 'l3prodNameD' )
    field_indices(f_rangfreq) =             add_ident ( 'rangFrequency' )
    field_indices(f_rangwavnum) =           add_ident ( 'rangWavenumber' )
    field_indices(f_time) =                 add_ident ( 'timeD' )
    ! Put parameter names into the symbol table
    parm_indices(p_dlat) =                  add_ident ( 'dLat' )
    parm_indices(p_dtype) =                 add_ident ( 'dgType' )
    parm_indices(p_imethod) =               add_ident ( 'intpMethod' )
    parm_indices(p_l2_map_lats) =           add_ident ( 'latGridMap' )
    parm_indices(p_l2_nom_lats) =           add_ident ( 'l2nomLats' )
    parm_indices(p_l2dgtype) =              add_ident ( 'L2dgType' )
    parm_indices(p_max_gap) =               add_ident ( 'MaxGap' )
    parm_indices(p_min_days) =              add_ident ( 'MinDays' )
    parm_indices(p_n)=                      add_ident ( 'N' )
    parm_indices(p_output_version_string) = add_ident ( 'OutputVersionString' )
    parm_indices(p_stype) =                 add_ident ( 'stdType' )
    ! Put section names into the symbol table
    section_indices(z_diagnostic) =         add_ident ( 'Diagnostic' )
    section_indices(z_globalsettings) =     add_ident ( 'GlobalSettings' )
    section_indices(z_output) =             add_ident ( 'Output' )
    section_indices(z_standard) =           add_ident ( 'Standard' )
    ! Put spec names into the symbol table
    spec_indices(s_dg) =                    add_ident ( 'dg' )
    spec_indices(s_out) =                   add_ident ( 'out' )
    spec_indices(s_std) =                   add_ident ( 'std' )

  ! Now initialize the units tables.  Init_Units depends on the lit tables
  ! having been initialized.

    call init_units

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
      begin, t+t_intp, l+l_csp, l+l_lin, n+n_dt_def, &
      begin, t+t_pmode, l+l_ado, l+l_all, l+l_asc, l+l_com, l+l_des, n+n_dt_def, &
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
      begin, s+s_std, &
             begin, f+f_prodname, t+t_string, n+n_field_type, &
             begin, f+f_time, t+t_string, n+n_field_type, &
             begin, f+f_mode, t+t_pmode, n+n_field_type, &
             begin, f+f_label, t+t_string, n+n_field_type, &
             begin, f+f_longrid, t+t_numeric_range, n+n_field_type, &
             begin, f+f_dlon, t+t_numeric, n+n_field_type, &
             begin, f+f_preslvl, t+t_numeric_range, n+n_field_type, &
             begin, f+f_alvl, t+t_numeric_range, n+n_field_type, &
             begin, f+f_dlvl, t+t_numeric_range, n+n_field_type, &
             begin, f+f_rangfreq, t+t_numeric_range, n+n_field_type, &
             begin, f+f_rangwavnum, t+t_numeric_range, n+n_field_type, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_dg, &
             begin, f+f_prodname, t+t_string, n+n_field_type, &
             begin, f+f_preslvl, t+t_numeric_range, n+n_field_type, &
             begin, f+f_longrid, t+t_numeric_range, n+n_field_type, &
             begin, f+f_dlon, t+t_numeric, n+n_field_type, &
             n+n_spec_def /) )
    call make_tree ( (/ &
      begin, s+s_out, &
             begin, f+f_mcf, t+t_string, n+n_field_type, &
             begin, f+f_file, t+t_string, n+n_field_type, &
             n+n_spec_def /) )
    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >
    call make_tree ( (/ &
      begin, z+z_globalsettings, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_l2_map_lats, t+t_numeric_range, n+n_name_def, &
             begin, p+p_dlat, t+t_numeric, n+n_name_def, &
             begin, p+p_l2_nom_lats, t+t_numeric, n+n_name_def, &
             begin, p+p_l2dgtype, t+t_string, n+n_name_def, &
             begin, p+p_stype, t+t_string, n+n_name_def, &
             begin, p+p_dtype, t+t_string, n+n_name_def, &
             begin, p+p_min_days, t+t_numeric, n+n_name_def, &
             begin, p+p_max_gap, t+t_numeric, n+n_name_def, &
             begin, p+p_n, t+t_numeric, n+n_name_def, &
             begin, p+p_imethod, t+t_intp, n+n_name_def, &
             n+n_section, &
      begin, z+z_standard, s+s_std, n+n_section, &
      begin, z+z_diagnostic, s+s_dg, n+n_section, &
      begin, z+z_output, s+s_out, n+n_section /) )
  end subroutine INIT_TABLES

! =====     Private procedures     =====================================
  ! --------------------------------------------------  MAKE_TREE  -----
  include "make_tree.f9h"

end module INIT_TABLES_MODULE

! $Log$
! Revision 1.10  2001/05/10 13:39:10  nakamura
! Removed separate pressure specifications for DZ.
!
! Revision 1.9  2001/05/04 18:43:22  nakamura
! Includes new DZ specs, latest parser changes.
!
! Revision 1.8  2001/04/11 18:53:24  nakamura
! Removed references to unused parser items.
!
! Revision 1.7  2001/04/04 17:46:41  pwagner
! compiles with new make_tree.f9h
!
! Revision 1.6  2001/03/16 16:33:43  nakamura
! Updated for parser compatibility.
!
! Revision 1.5  2001/02/21 21:19:40  nakamura
! Removed l2Ver.
!
! Revision 1.4  2001/01/18 16:53:25  nakamura
! Moved minDays from PCF to cf.
!
! Revision 1.3  2001/01/16 17:53:03  nakamura
! Added a template for the log file name to GlobalSettings.
!
! Revision 1.2  2000/12/29 21:46:35  nakamura
! Removed bypass flag, quantities.
!
! Revision 1.1  2000/10/24 19:47:17  nakamura
! Customized for L3.
!
