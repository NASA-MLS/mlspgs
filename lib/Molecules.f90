! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Molecules

! Intrinsic constants needed by Init_Tables_Module, DeclarationTable, etc.

! Declaring the definitions is handled by the tree walker.

  use INTRINSIC ! Everything

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Definitions of types:
  integer, parameter :: FIRST_MOLECULE_TYPE = last_intrinsic_type + 1
  integer, parameter :: T_MOLECULE =          first_molecule_type
  integer, parameter :: LAST_MOLECULE_TYPE =  t_molecule

! Definitions of the literals:
  integer, parameter :: FIRST_MOLECULE = last_intrinsic_lit + 1
  integer, parameter :: L_AIR_CONT =     first_molecule
  integer, parameter :: L_BRO =          l_air_cont + 1
  integer, parameter :: L_BR_79_O =      l_bro + 1
  integer, parameter :: L_BR_81_O =      l_br_79_o + 1
  integer, parameter :: L_C_13_O =       l_br_81_o + 1
  integer, parameter :: L_CH3CL_35 =     l_c_13_o + 1
  integer, parameter :: L_CH3CL_37 =     l_ch3cl_35 + 1
  integer, parameter :: L_CH3CN =        l_ch3cl_37 + 1
  integer, parameter :: L_CL_35_NO3 =    l_ch3cn + 1
  integer, parameter :: L_CL_37_NO3 =    l_cl_35_no3 + 1
  integer, parameter :: L_CL_37_O =      l_cl_37_no3 + 1
  integer, parameter :: L_CL_35_O =      l_cl_37_o + 1
  integer, parameter :: L_CLO =          l_cl_35_o + 1
  integer, parameter :: L_CO =           l_clo + 1
  integer, parameter :: L_CO_18 =        l_co + 1
  integer, parameter :: L_COF2 =         l_co_18 + 1
  integer, parameter :: L_H2C_13_O =     l_cof2 + 1
  integer, parameter :: L_H2CO =         l_h2c_13_o + 1
  integer, parameter :: L_H2CO_18 =      l_h2co + 1
  integer, parameter :: L_H2O =          l_h2co_18 + 1
  integer, parameter :: L_H2O2 =         l_h2o + 1
  integer, parameter :: L_H2O_18 =       l_h2o2 + 1
  integer, parameter :: L_H2O_V2 =       l_h2o_18 + 1
  integer, parameter :: L_HC_13_N =      l_h2o_v2 + 1
  integer, parameter :: L_HCL =          l_hc_13_n + 1
  integer, parameter :: L_HCL_35 =       l_hcl + 1
  integer, parameter :: L_HCL_37 =       l_hcl_35 + 1
  integer, parameter :: L_HCN =          l_hcl_37 + 1
  integer, parameter :: L_HCN_15 =       l_hcn + 1
  integer, parameter :: L_HDO =          l_hcn_15 + 1
  integer, parameter :: L_HNO3 =         l_hdo + 1
  integer, parameter :: L_HO2 =          l_hno3 + 1
  integer, parameter :: L_HOCL =         l_ho2 + 1
  integer, parameter :: L_HOCL_37 =      l_hocl + 1
  integer, parameter :: L_LIQ_H2O =      l_hocl_37 + 1
  integer, parameter :: L_N2 =           l_liq_h2o + 1
  integer, parameter :: L_N2O =          l_n2 + 1
  integer, parameter :: L_NO =           l_n2o + 1
  integer, parameter :: L_NO2 =          l_no + 1
  integer, parameter :: L_O2 =           l_no2 + 1
  integer, parameter :: L_O2_V1 =        l_o2 + 1
  integer, parameter :: L_O3 =           l_o2_v1 + 1
  integer, parameter :: L_O3_2V2 =       l_o3 + 1
  integer, parameter :: L_O3_A_O18_V2 =  l_o3_2v2 + 1
  integer, parameter :: L_O3_ASYM_O_17 = l_o3_a_o18_v2 + 1
  integer, parameter :: L_O3_ASYM_O_18 = l_o3_asym_o_17 + 1
  integer, parameter :: L_O3_S_O18_V2 =  l_o3_asym_o_18 + 1
  integer, parameter :: L_O3_SYM_O_17 =  l_o3_s_o18_v2 + 1
  integer, parameter :: L_O3_SYM_O_18 =  l_o3_sym_o_17 + 1
  integer, parameter :: L_O3_V1_3 =      l_o3_sym_o_18 + 1
  integer, parameter :: L_O3_V2 =        l_o3_v1_3 + 1
  integer, parameter :: L_O_17_O =       l_o3_v2 + 1
  integer, parameter :: L_O_18_O =       l_o_17_o + 1
  integer, parameter :: L_OC_34_S =      l_o_18_o + 1
  integer, parameter :: L_OCL_35_O =     l_oc_34_s + 1
  integer, parameter :: L_OCL_37_O =     l_ocl_35_o + 1
  integer, parameter :: L_OCS =          l_ocl_37_o + 1
  integer, parameter :: L_OH =           l_ocs + 1
  integer, parameter :: L_SO2 =          l_oh + 1
  integer, parameter :: LAST_MOLECULE =  l_so2

! Mapping from the literals to the spec tags:
  integer, public :: SPEC_TAGS(FIRST_MOLECULE:LAST_MOLECULE)
  data spec_tags(l_air_cont)     / 00028964 /
  data spec_tags(l_bro)          / 00095001 /
  data spec_tags(l_br_79_o)      / 00095001 /
  data spec_tags(l_br_81_o)      / 00097001 /
  data spec_tags(l_c_13_o)       / 00029001 /
  data spec_tags(l_ch3cl_35)     / 00050007 /
  data spec_tags(l_ch3cl_37)     / 00052009 /
  data spec_tags(l_ch3cn)        / 00041001 /
  data spec_tags(l_cl_35_no3)    / 00097002 /
  data spec_tags(l_cl_37_no3)    / 00099001 /
  data spec_tags(l_cl_35_o)      / 00051002 /
  data spec_tags(l_cl_37_o)      / 00053002 /
  data spec_tags(l_clo)          / 00051002 /
  data spec_tags(l_co)           / 00028001 /
  data spec_tags(l_co_18)        / 00030001 /
  data spec_tags(l_cof2)         / 00066001 /
  data spec_tags(l_h2c_13_o)     / 00031002 /
  data spec_tags(l_h2co)         / 00030004 /
  data spec_tags(l_h2co_18)      / 00032004 /
  data spec_tags(l_h2o)          / 00018003 /
  data spec_tags(l_h2o2)         / 00034004 /
  data spec_tags(l_h2o_18)       / 00020003 /
  data spec_tags(l_h2o_v2)       / 00018005 /
  data spec_tags(l_hc_13_n)      / 00028002 /
  data spec_tags(l_hcl)          / 00036001 /
  data spec_tags(l_hcl_35)       / 00036001 /
  data spec_tags(l_hcl_37)       / 00038001 /
  data spec_tags(l_hcn)          / 00027001 /
  data spec_tags(l_hcn_15)       / 00028003 /
  data spec_tags(l_hdo)          / 00019002 /
  data spec_tags(l_hno3)         / 00063001 /
  data spec_tags(l_ho2)          / 00033001 /
  data spec_tags(l_hocl)         / 00052006 /
  data spec_tags(l_hocl_37)      / 00054005 /
  data spec_tags(l_liq_h2o)      / 00018999 /
  data spec_tags(l_n2)           / 00028964 /
  data spec_tags(l_n2o)          / 00044004 /
  data spec_tags(l_no)           / 00030008 /
  data spec_tags(l_no2)          / 00046006 /
  data spec_tags(l_o2)           / 00032001 /
  data spec_tags(l_o2_v1)        / 00032002 /
  data spec_tags(l_o3)           / 00048004 /
  data spec_tags(l_o3_2v2)       / 00048007 /
  data spec_tags(l_o3_a_o18_v2)  / 00050006 /
  data spec_tags(l_o3_asym_o_17) / 00049002 /
  data spec_tags(l_o3_asym_o_18) / 00050004 /
  data spec_tags(l_o3_s_o18_v2)  / 00050005 /
  data spec_tags(l_o3_sym_o_17)  / 00049001 /
  data spec_tags(l_o3_sym_o_18)  / 00050003 /
  data spec_tags(l_o3_v1_3)      / 00048006 /
  data spec_tags(l_o3_v2)        / 00048005 /
  data spec_tags(l_o_17_o)       / 00033002 /
  data spec_tags(l_o_18_o)       / 00034001 /
  data spec_tags(l_oc_34_s)      / 00062001 /
  data spec_tags(l_ocl_35_o)     / 00067001 /
  data spec_tags(l_ocl_37_o)     / 00069001 /
  data spec_tags(l_ocs)          / 00060001 /
  data spec_tags(l_oh)           / 00017001 /
  data spec_tags(l_so2)          / 00064002 /

contains ! =====     Public procedures     =============================
! -----------------------------------------------  INIT_MOLECULES  -----
  subroutine Init_Molecules ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
    use TREE_TYPES, only: N_DT_DEF

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    integer :: I    ! Loop inductor in an array constructor

    ! Initialize the intrinsic types

    call init_intrinsic ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices )

    ! Put type names into the symbol table:

    data_type_indices(t_molecule) =        add_ident ( 'molecule' )

    ! Put literals into the symbol table:

    lit_indices(l_air_cont) =     add_ident ( 'air_cont' )
    lit_indices(l_bro) =          add_ident ( 'bro' )
    lit_indices(l_br_79_o) =      add_ident ( 'br_79_o' )
    lit_indices(l_br_81_o) =      add_ident ( 'br_81_o' )
    lit_indices(l_c_13_o) =       add_ident ( 'c_13_o' )
    lit_indices(l_ch3cl_35) =     add_ident ( 'ch3cl_35' )
    lit_indices(l_ch3cl_37) =     add_ident ( 'ch3cl_37' )
    lit_indices(l_ch3cn) =        add_ident ( 'ch3cn' )
    lit_indices(l_cl_35_no3) =    add_ident ( 'cl_35_no3' )
    lit_indices(l_cl_37_no3) =    add_ident ( 'cl_37_no3' )
    lit_indices(l_cl_35_o) =      add_ident ( 'cl_35_o' )
    lit_indices(l_cl_37_o) =      add_ident ( 'cl_37_o' )
    lit_indices(l_clo) =          add_ident ( 'clo' )
    lit_indices(l_co) =           add_ident ( 'co' )
    lit_indices(l_co_18) =        add_ident ( 'co_18' )
    lit_indices(l_cof2) =         add_ident ( 'cof2' )
    lit_indices(l_h2c_13_o) =     add_ident ( 'h2c_13_o' )
    lit_indices(l_h2co) =         add_ident ( 'h2co' )
    lit_indices(l_h2co_18) =      add_ident ( 'h2co_18' )
    lit_indices(l_h2o) =          add_ident ( 'h2o' )
    lit_indices(l_h2o2) =         add_ident ( 'h2o2' )
    lit_indices(l_h2o_18) =       add_ident ( 'h2o_18' )
    lit_indices(l_h2o_v2) =       add_ident ( 'h2o_v2' )
    lit_indices(l_hc_13_n) =      add_ident ( 'hc_13_n' )
    lit_indices(l_hcl) =          add_ident ( 'hcl' )
    lit_indices(l_hcl_35) =       add_ident ( 'hcl_35' )
    lit_indices(l_hcl_37) =       add_ident ( 'hcl_37' )
    lit_indices(l_hcn) =          add_ident ( 'hcn' )
    lit_indices(l_hcn_15) =       add_ident ( 'hcn_15' )
    lit_indices(l_hdo) =          add_ident ( 'hdo' )
    lit_indices(l_hno3) =         add_ident ( 'hno3' )
    lit_indices(l_ho2) =          add_ident ( 'ho2' )
    lit_indices(l_hocl) =         add_ident ( 'hocl' )
    lit_indices(l_hocl_37) =      add_ident ( 'hocl_37' )
    lit_indices(l_liq_h2o) =      add_ident ( 'liq_h2o' )
    lit_indices(l_n2) =           add_ident ( 'n2' )
    lit_indices(l_n2o) =          add_ident ( 'n2o' )
    lit_indices(l_no) =           add_ident ( 'no' )
    lit_indices(l_no2) =          add_ident ( 'no2' )
    lit_indices(l_o2) =           add_ident ( 'o2' )
    lit_indices(l_o2_v1) =        add_ident ( 'o2_v1' )
    lit_indices(l_o3) =           add_ident ( 'o3' )
    lit_indices(l_o3_2v2) =       add_ident ( 'o3_2v2' )
    lit_indices(l_o3_a_o18_v2) =  add_ident ( 'o3_a_o18_v2' )
    lit_indices(l_o3_asym_o_17) = add_ident ( 'o3_asym_o_17' )
    lit_indices(l_o3_asym_o_18) = add_ident ( 'o3_asym_o_18' )
    lit_indices(l_o3_s_o18_v2) =  add_ident ( 'o3_s_o18_v2' )
    lit_indices(l_o3_sym_o_17) =  add_ident ( 'o3_sym_o_17' )
    lit_indices(l_o3_sym_o_18) =  add_ident ( 'o3_sym_o_18' )
    lit_indices(l_o3_v1_3) =      add_ident ( 'o3_v1_3' )
    lit_indices(l_o3_v2) =        add_ident ( 'o3_v2' )
    lit_indices(l_o_17_o) =       add_ident ( 'o_17_o' )
    lit_indices(l_o_18_o) =       add_ident ( 'o_18_o' )
    lit_indices(l_oc_34_s) =      add_ident ( 'oc_34_s' )
    lit_indices(l_ocl_35_o) =     add_ident ( 'ocl_35_o' )
    lit_indices(l_ocl_37_o) =     add_ident ( 'ocl_37_o' )
    lit_indices(l_ocs) =          add_ident ( 'ocs' )
    lit_indices(l_oh) =           add_ident ( 'oh' )
    lit_indices(l_so2) =          add_ident ( 'so2' )

    ! Create the type tree for the molecule type
    call make_tree ( (/ &
      begin, t+t_molecule, l+(/ (i,i=first_molecule, last_molecule) /), &
             n+n_dt_def /) )

  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine INIT_MOLECULES

end module MOLECULES

! $Log$
! Revision 2.9  2001/05/10 23:28:39  livesey
! Added extra molecules for non-isotopes etc.
!
! Revision 2.8  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.7  2001/04/05 01:27:26  vsnyder
! Add hcl_35
!
! Revision 2.6  2001/04/04 17:56:42  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 2.5  2001/04/03 19:09:12  vsnyder
! Change the order of initialization to intrinsic, Molecules, MLSSignals.
! Use the revised make_tree.f9h, which requires revision of init...
! calling sequences.
!
! Revision 2.4  2001/04/03 01:17:48  vsnyder
! CVS stuff was still wrong
!
! Revision 2.3  2001/04/03 00:32:25  vsnyder
! Correct CVS stuff
!

