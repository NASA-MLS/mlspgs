! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

MODULE INIT_TABLES_MODULE

! Preload the string and symbol tables with specification and field
! names.  Preload the tree with definitions of types, lits, fields,
! specifications, sections ....

! Declaring the definitions is handled by the tree walker.

  USE Init_MLSSignals_m ! Everything. Init_MLSSignals, Field_First,
    ! Last_Signal_Field, Spec_First, Last_Signal_Spec, Numerous S_....
  USE INTRINSIC ! Everything. ADD_IDENT, BEGIN, D, F, FIRST_LIT,
    ! FIRST_MOLECULE,  INIT_INTRINSIC, L, L_<several>, LAST_INTRINSIC_LIT,
    ! LAST_MOLECULE, MAKE_TREE, N, NADP, ND, NDP, NP, NR, P, S, T,
    ! T_BOOLEAN, T_FIRST, T_LAST_INTRINSIC, T_NUMERIC, T_NUMERIC_RANGE,
    ! T_STRING and Z are used here, but everything is included so that it
    ! can be gotten by USE INIT_TABLES_MODULE.

  IMPLICIT NONE
  PUBLIC ! This would be a MUCH LONGER list than the list of private
  !        names below.
  PRIVATE :: ADD_IDENT, INIT_INTRINSIC, MAKE_TREE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=256), PRIVATE :: Id = &
       "$Id$"
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Enumeration types:

  INTEGER, PUBLIC, PARAMETER :: T_USE            = last_signal_type+1
  INTEGER, PUBLIC, PARAMETER :: T_UNITS          = t_use+1
  INTEGER, PUBLIC, PARAMETER :: T_MODULE         = t_units+1
  INTEGER, PUBLIC, PARAMETER :: T_ENABLE         = t_module+1
  INTEGER, PUBLIC, PARAMETER :: T_LAST           = t_enable

! Field indices:

  INTEGER, PUBLIC, PARAMETER :: F_MIFS = last_Signal_Field + 1
  INTEGER, PUBLIC, PARAMETER :: F_USE = f_mifs + 1
  INTEGER, PUBLIC, PARAMETER :: F_SECONDARY = f_use + 1
  INTEGER, PUBLIC, PARAMETER :: F_S = f_secondary + 1
  INTEGER, PUBLIC, PARAMETER :: F_BANDNO = f_s + 1
  INTEGER, PUBLIC, PARAMETER :: FIELD_LAST = f_bandno + 1

! Enumeration literals:

  INTEGER, PUBLIC, PARAMETER :: L_MATCH   =  last_signal_lit + 1
  INTEGER, PUBLIC, PARAMETER :: L_OVERRIDE = l_match + 1
  INTEGER, PUBLIC, PARAMETER :: LAST_LIT =      l_override + 1

! Section identities:

  INTEGER, PUBLIC, PARAMETER :: Z_GLOBALSETTINGS = 1
  INTEGER, PUBLIC, PARAMETER :: Z_CALIBRATION = Z_GLOBALSETTINGS + 1
  INTEGER, PUBLIC, PARAMETER :: Z_OUTPUT = Z_CALIBRATION + 1
  INTEGER, PUBLIC, PARAMETER :: SECTION_FIRST = z_globalSettings, &
                                SECTION_LAST = Z_OUTPUT

! Specification indices:

  INTEGER, PUBLIC, PARAMETER :: S_SPACEMIFS = last_Signal_Spec + 1
  INTEGER, PUBLIC, PARAMETER :: S_TARGETMIFS = s_spaceMIFs + 1
  INTEGER, PUBLIC, PARAMETER :: S_LIMBMIFS = s_targetMIFs + 1
  INTEGER, PUBLIC, PARAMETER :: S_DISCARDMIFS = s_limbMIFs + 1
  INTEGER, PUBLIC, PARAMETER :: S_SWITCH = s_discardMIFs + 1
  INTEGER, PUBLIC, PARAMETER :: S_CHI2ERR = s_switch + 1
  INTEGER, PUBLIC, PARAMETER :: SPEC_LAST = s_chi2err

! Parameter names:

  ! In GlobalSettings section:

  INTEGER, PUBLIC, PARAMETER :: P_OUTPUT_VERSION_STRING = spec_last + 1
  INTEGER, PUBLIC, PARAMETER :: P_VERSION_COMMENT = p_output_version_string + 1
  INTEGER, PUBLIC, PARAMETER :: P_PRODUCE_L1BOA = p_version_comment + 1
  INTEGER, PUBLIC, PARAMETER :: P_SIMOA = p_produce_l1boa + 1

  ! In Calibration section:

  INTEGER, PUBLIC, PARAMETER :: P_CALWINDOW = p_simoa + 1
  INTEGER, PUBLIC, PARAMETER :: P_USEDEFAULTGAINS = p_calwindow + 1
  INTEGER, PUBLIC, PARAMETER :: P_CALIBDACS = p_usedefaultgains + 1
  INTEGER, PUBLIC, PARAMETER :: P_GHZSPACETEMP = p_calibdacs + 1
  INTEGER, PUBLIC, PARAMETER :: P_GHZTARGETTEMP = p_GHzSpaceTemp + 1
  INTEGER, PUBLIC, PARAMETER :: P_THZSPACETEMP = p_GHzTargetTemp + 1
  INTEGER, PUBLIC, PARAMETER :: P_THZTARGETTEMP = p_THzSpaceTemp + 1
  INTEGER, PUBLIC, PARAMETER :: P_THZSpaceAngle = p_THzTargetTemp + 1
  INTEGER, PUBLIC, PARAMETER :: P_MIF_DURATION = p_THzSpaceAngle + 1
  INTEGER, PUBLIC, PARAMETER :: P_MIF_DEAD_TIME = p_mif_duration + 1
  INTEGER, PUBLIC, PARAMETER :: P_MIFsPerMAF = p_mif_dead_time + 1
  INTEGER, PUBLIC, PARAMETER :: P_THzMaxBias = p_MIFsPerMAF + 1
  INTEGER, PUBLIC, PARAMETER :: P_MoonToSpaceAngle = p_THzMaxBias + 1
  INTEGER, PUBLIC, PARAMETER :: P_MoonToLimbAngle = p_MoonToSpaceAngle + 1

  ! In Output section:

  INTEGER, PUBLIC, PARAMETER :: P_REMOVEBASELINE = P_MoonToLimbAngle + 1

  INTEGER, PUBLIC, PARAMETER :: FIRST_PARM = P_OUTPUT_VERSION_STRING
  INTEGER, PUBLIC, PARAMETER :: LAST_PARM = P_REMOVEBASELINE

! Table for section ordering:

  INTEGER, PUBLIC, PARAMETER :: OK = 1, & ! NO = 0
    SECTION_ORDERING(section_first:section_last, &
                     section_first-1:section_last) = RESHAPE( &
! To: | globalSettings        |
!     |      Calibration      |
!     |            Output     |
! ====|==============================|== From: ==
        (/OK,    0,  0,    & ! Start
           0,   OK,  0,    & ! GlobalSettings
           0,   OK,  OK,   & ! Calibration
           0,    0,  OK /) & ! Output
!       , shape(section_ordering) )
        , (/ section_last-section_first+1, section_last-section_first+2 /) )

CONTAINS ! =====     Public procedures     =============================
! --------------------------------------------------  INIT_TABLES  -----
  SUBROUTINE INIT_TABLES

    USE TREE_TYPES, ONLY: N_DOT, N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, &
         N_NAME_DEF, N_SECTION, N_SPEC_DEF, N_PLUS

  ! Put intrinsic predefined identifiers into the symbol table.

     CALL init_MLSSignals ( t_last, field_last, last_lit, &
     & first_parm, last_parm, section_last, spec_last )

  ! Put nonintrinsic predefined identifiers into the symbol table.

    ! Put enumeration type names into the symbol table

    data_type_indices(t_use) =              add_ident ( 'use' )
    data_type_indices(t_module) =           add_ident ( 'module' )
    data_type_indices(t_units) =            add_ident ( 'units' )

    ! Put enumeration literals into the symbol table:

    lit_indices(l_match) =                   add_ident ( 'match' )
    lit_indices(l_override) =                add_ident ( 'override' )

    ! Put field names into the symbol table

    field_indices(f_mifs) =                 add_ident ( 'MIFs' )
    field_indices(f_use) =                  add_ident ( 'use' )
    field_indices(f_secondary) =            add_ident ( 'secondary' )
    field_indices(f_s) =                    add_ident ( 's' )
    field_indices(f_bandno) =               add_ident ( 'bandno' )
 
    ! Put parameter names into the symbol table

    parm_indices(p_calwindow)=              add_ident ( 'CalWindow' )
    parm_indices(p_usedefaultgains)=        add_ident ( 'UseDefaultGains' )
    parm_indices(p_calibDACS)=              add_ident ( 'CalibDACS' )
    parm_indices(p_GHzSpaceTemp)=           add_ident ( 'GHzSpaceTemp' )
    parm_indices(p_GHzTargetTemp)=          add_ident ( 'GHzTargetTemp' )
    parm_indices(p_THzSpaceTemp)=           add_ident ( 'THzSpaceTemp' )
    parm_indices(p_THzTargetTemp)=          add_ident ( 'THzTargetTemp' )
    parm_indices(p_THzSpaceAngle)=          add_ident ( 'THzSpaceAngle' )
    parm_indices(p_THzMaxBias)=             add_ident ( 'THzMaxBias' )
    parm_indices(p_mif_duration)=           add_ident ( 'MIF_Duration' )
    parm_indices(p_mif_dead_time)=          add_ident ( 'MIF_DeadTime' )
    parm_indices(p_mifspermaf)=             add_ident ( 'MIFsPerMAF' )
    parm_indices(p_output_version_string) = add_ident ( 'OutputVersionString' )
    parm_indices(p_version_comment) =       add_ident ( 'VersionComment' )
    parm_indices(p_produce_l1boa)=          add_ident ( 'ProduceL1BOA' )
    parm_indices(p_simoa)=                  add_ident ( 'SimOA' )
    parm_indices(p_removebaseline)=         add_ident ( 'RemoveBaseline' )
    parm_indices(p_MoonToSpaceAngle)=       add_ident ( 'MoonToSpaceAngle' )
    parm_indices(p_MoonToLimbAngle)=        add_ident ( 'MoonToLimbAngle' )

    ! Put section names into the symbol table

    section_indices(z_calibration) =        add_ident ( 'Calibration' )
    section_indices(z_globalsettings) =     add_ident ( 'GlobalSettings' )
    section_indices(z_output) =             add_ident ( 'Output' )

    ! Put spec names into the symbol table

    spec_indices(s_spaceMIFs) =               add_ident ( 'spaceMIFs' )
    spec_indices(s_targetMIFs) =              add_ident ( 'targetMIFs' )
    spec_indices(s_limbMIFs) =                add_ident ( 'limbMIFs' )
    spec_indices(s_discardMIFs) =             add_ident ( 'discardMIFs' )
    spec_indices(s_switch) =                  add_ident ( 'switch' )
    spec_indices(s_chi2err) =                 add_ident ( 'EnableChi2Err' )

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

    CALL make_tree ( (/ &
      begin, t+t_use, l+l_match, l+l_override, n+n_dt_def, &
      begin, t+t_module, l+l_ghz, l+l_thz, n+n_dt_def, &
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

    CALL make_tree ( (/ &
      begin, s+s_spaceMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, n+n_field_type, &
             begin, f+f_use, t+t_use, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_type, &
             nadp+n_spec_def /) )

    CALL make_tree ( (/ &
      begin, s+s_targetMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, nr+n_field_type, &
             begin, f+f_use, t+t_use, nr+n_field_type, &
             begin, f+f_module, s+s_module, nr+n_field_type, &
             begin, f+f_secondary, t+t_boolean, n+n_field_type, &
             ndp+n_spec_def /) )

    CALL make_tree ( (/ &
      begin, s+s_limbMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, n+n_field_type, &
             begin, f+f_use, t+t_use, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_type, &
             nadp+n_spec_def /) )

    CALL make_tree ( (/ &
      begin, s+s_discardMIFs, &
             begin, f+f_mifs, t+t_numeric_range, t+t_numeric, n+n_field_type, &
             begin, f+f_use, t+t_use, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_type, &
             nadp+n_spec_def /) )

    CALL make_tree ( (/ &
      begin, s+s_switch, &
             begin, f+f_s, t+t_numeric, n+n_field_type, &
             begin, f+f_bandno, t+t_numeric, n+n_field_type, &
             nadp+n_spec_def /) )

    CALL make_tree ( (/ &
      begin, s+s_chi2err, &
             begin, f+f_bandno, t+t_numeric_range, t+t_numeric, n+n_field_type,&
             nadp+n_spec_def /) )

    ! Define the relations between sections and specs.  These are
    ! represented by trees of the form
    !  < n_section section_name
    !              < n_name_def p_parameter t_type ... t_type > ...
    !  > or
    !  < n_section section_name s_spec ... s_spec >

    CALL make_tree ( (/ &
      begin, z+z_globalsettings, &
             begin, p+p_version_comment, t+t_string, n+n_name_def, &
             begin, p+p_output_version_string, t+t_string, n+n_name_def, &
             begin, p+p_produce_l1boa, t+t_boolean, n+n_name_def, &
             begin, p+p_simoa, t+t_boolean, n+n_name_def, &
             n+n_section, &
      begin, z+z_calibration, &
             begin, p+p_calwindow, t+t_numeric, n+n_name_def, &
             begin, p+p_GHzSpaceTemp, t+t_numeric, n+n_name_def, &
             begin, p+p_GHzTargetTemp, t+t_numeric, n+n_name_def, &
             begin, p+p_THzSpaceTemp, t+t_numeric, n+n_name_def, &
             begin, p+p_THzTargetTemp, t+t_numeric, n+n_name_def, &
             begin, p+p_THzSpaceAngle, t+t_numeric, n+n_name_def, &
             begin, p+p_THzMaxBias, t+t_numeric, n+n_name_def, &
             begin, p+p_MoonToSpaceAngle, t+t_numeric, n+n_name_def, &
             begin, p+p_MoonToLimbAngle, t+t_numeric, n+n_name_def, &
             begin, p+p_mif_duration, t+t_numeric, n+n_name_def, &
             begin, p+p_mif_dead_time, t+t_numeric, n+n_name_def, &
             begin, p+p_mifspermaf, t+t_numeric, n+n_name_def, &
             begin, p+p_usedefaultgains, t+t_boolean, n+n_name_def, &
             begin, p+p_calibDACS, t+t_boolean, n+n_name_def, &
             s+s_spaceMIFs, s+s_targetMIFs, s+s_limbMIFS, s+s_discardMIFs, &
             s+s_switch, n+n_section, &
      begin, z+z_output, &
             begin, p+p_removebaseline, t+t_boolean, n+n_name_def, &
             s+s_chi2err, n+n_section/) )

  END SUBROUTINE INIT_TABLES
    
  ! --------------------------------------------------  MAKE_TREE  -----
  INCLUDE "make_tree.f9h"

END MODULE INIT_TABLES_MODULE
  
! $Log$
! Revision 2.17  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.16  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.15  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.14  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! $Log$
! Revision 2.17  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.16  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.14  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.13  2002/11/14 16:49:57  perun
! Split space & target temps between GHz & THz
!
! Revision 2.12  2002/11/07 21:34:42  jdone
! Added HDF4/HDF5 switch.
!
! Revision 2.11  2002/03/29 20:18:34  perun
! Version 1.0 commit
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
