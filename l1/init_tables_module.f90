! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  use Init_MLSSignals_m ! Everything. Init_MLSSignals, Field_First,
    ! Last_Signal_Field, Spec_First, Last_Signal_Spec, Numerous S_....
  use INTRINSIC ! Everything. ADD_IDENT, BEGIN, D, F, FIRST_LIT,
    ! FIRST_MOLECULE,  INIT_INTRINSIC, L, L_<several>, LAST_INTRINSIC_LIT,
    ! LAST_MOLECULE, MAKE_TREE, N, NADP, ND, NDP, NP, NR, P, S, T,
    ! T_BOOLEAN, T_FIRST, T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE,
    ! T_STRING and Z are used here, but everything is included so that it
    ! can be gotten by USE INIT_TABLES_MODULE.

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

  integer, public, parameter :: T_USE            = last_signal_type+1
  integer, public, parameter :: T_UNITS          = t_use+1
  integer, public, parameter :: T_LAST           = t_units

! Field indices:

  integer, public, parameter :: F_MIFS = last_Signal_Field + 1
  integer, public, parameter :: F_USE = f_mifs + 1
  integer, public, parameter :: FIELD_LAST = f_use

! Enumeration literals:

  integer, public, parameter :: L_EXPECTED   = last_signal_lit + 1
  integer, public, parameter :: L_MASK = l_expected + 1
  integer, public, parameter :: L_OVERRIDE = l_mask + 1
  integer, public, parameter :: LAST_LIT = l_override

! Section identities:

  integer, public, parameter :: Z_GLOBALSETTINGS = 1
  integer, public, parameter :: Z_CALIBRATION = 2
  integer, public, parameter :: SECTION_FIRST = z_globalSettings, &
                                SECTION_LAST = z_Calibration

! Specification indices:

  integer, public, parameter :: S_SPACEMIFS = last_Signal_Spec + 1
  integer, public, parameter :: S_TARGETMIFS = s_spacemifs + 1
  integer, public, parameter :: SPEC_LAST = s_targetMIFs

! Parameter names:

  ! In GlobalSettings section:

  integer, public, parameter :: P_OUTPUT_VERSION_STRING = spec_last + 1
  integer, public, parameter :: P_VERSION_COMMENT = p_output_version_string +1

  ! In Calibration section:

  integer, public, parameter :: P_CALWINDOW = p_version_comment + 1

  integer, public, parameter :: FIRST_PARM = P_OUTPUT_VERSION_STRING
  integer, public, parameter :: LAST_PARM = P_CALWINDOW

! Table for section ordering:

  integer, public, parameter :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = reshape( &
! To: | globalSettings       |
!     |      Calibration     |
! ====|==============================|== From: ==
        (/OK,    0,  & ! Start
           0,   OK,  & ! GlobalSettings
           0,   OK/) & ! Calibration
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

contains ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  subroutine INIT_TABLES

    use TREE_TYPES, only: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
         N_NAME_DEF, N_SECTION, N_SPEC_DEF, N_PLUS

  ! Put intrinsic predefined identifiers into the symbol table.

     call init_MLSSignals ( t_last, field_last, last_lit, &
     & first_parm, last_parm, section_last, spec_last )

  ! Put nonintrinsic predefined identifiers into the symbol table.

    ! Put enumeration type names into the symbol table

    data_type_indices(t_use) =              add_ident ( 'use' )
    data_type_indices(t_units) =            add_ident ( 'units' )

    ! Put enumeration literals into the symbol table:

    lit_indices(l_expected) =                add_ident ( 'expected' )
    lit_indices(l_mask) =                    add_ident ( 'mask' )
    lit_indices(l_override) =                add_ident ( 'override' )

    ! Put field names into the symbol table

    field_indices(f_mifs) =                 add_ident ( 'MIFs' )
    field_indices(f_use) =                  add_ident ( 'use' )
 
    ! Put parameter names into the symbol table

    parm_indices(p_calwindow)=              add_ident ( 'CalWindow' )
    parm_indices(p_output_version_string) = add_ident ( 'OutputVersionString' )
    parm_indices(p_version_comment) =       add_ident ( 'VersionComment' )

    ! Put section names into the symbol table

    section_indices(z_calibration) =        add_ident ( 'Calibration' )
    section_indices(z_globalsettings) =     add_ident ( 'GlobalSettings' )

    ! Put spec names into the symbol table

    spec_indices(s_spaceMIFs) =               add_ident ( 'spaceMIFs' )
    spec_indices(s_targetMIFs) =              add_ident ( 'targetMIFs' )

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

    ! Define the enumerated types

    call make_tree ( (/ &
      begin, t+t_use, l+l_expected, l+l_mask, l+l_override, n+n_dt_def, &
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
      begin, s+s_spaceMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, n+n_field_type, &
             begin, f+f_use, t+t_use, n+n_field_type, &
             nadp+n_spec_def /) )

    call make_tree ( (/ &
      begin, s+s_targetMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, n+n_field_type, &
             begin, f+f_use, t+t_use, n+n_field_type, &
             nadp+n_spec_def /) )

    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >

    call make_tree ( (/ &
      begin, z+z_globalsettings, &
             begin, p+p_version_comment, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             n+n_section, &
      begin, z+z_calibration, &
             begin, p+p_calwindow, t+t_numeric, n+n_name_def, &
             s+s_spaceMIFs, s+s_targetMIFs, n+n_section /) )

  end subroutine INIT_TABLES
    
  ! --------------------------------------------------  MAKE_TREE  -----
  include "make_tree.f9h"

end module INIT_TABLES_MODULE    

! $Log$
! Revision 2.10  2001/04/27 14:01:14  perun
! For the latest parser version.
!
! Revision 2.9  2001/04/05 14:43:56  perun
! Another change to init_MLSSignals
!
! Revision 2.8  2001/03/16 15:13:47  perun
! Another change in call to Init_MLSSignals
!
! Revision 2.7  2001/03/14 15:59:27  perun
! Use Van's latest with Init_MLSSignals_m module
!
! Revision 2.6  2001/03/05 16:46:37  perun
! Use Van's include method
!
! Revision 2.5  2001/02/23 19:05:40  perun
! *** empty log message ***
!
! Revision 2.1  2001/02/23 18:57:58  perun
! Version 0.5 commit
!
