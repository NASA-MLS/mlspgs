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
  private :: not_used_here 
!---------------------------------------------------------------------------

! Definitions of types:
  integer, parameter :: FIRST_MOLECULE_TYPE = last_intrinsic_type + 1
  integer, parameter :: T_MOLECULE =          first_molecule_type
  integer, parameter :: LAST_MOLECULE_TYPE =  t_molecule

! Definitions of the literals:
  integer, parameter :: FIRST_MOLECULE = last_intrinsic_lit + 1
  integer, parameter :: L_AIR_CONT =     first_molecule,     SP_AIR_CONT     = 00028964
  integer, parameter :: L_BRO =          l_air_cont + 1,     SP_BRO          = 00095001
  integer, parameter :: L_BR_79_O =      l_bro + 1,          SP_BR_79_O      = 00095001
  integer, parameter :: L_BR_81_O =      l_br_79_o + 1,      SP_BR_81_O      = 00097001
  integer, parameter :: L_C_13_O =       l_br_81_o + 1,      SP_C_13_O       = 00029001
  integer, parameter :: L_CH3CL =        l_c_13_o + 1,       SP_CH3CL        = 00050007
  integer, parameter :: L_CH3CL_35 =     l_ch3cl + 1,        SP_CH3CL_35     = 00050007
  integer, parameter :: L_CH3CL_37 =     l_ch3cl_35 + 1,     SP_CH3CL_37     = 00052009
  integer, parameter :: L_CH3CN =        l_ch3cl_37 + 1,     SP_CH3CN        = 00041001
  integer, parameter :: L_CL_35_NO3 =    l_ch3cn + 1,        SP_CL_35_NO3    = 00097002
  integer, parameter :: L_CL_37_NO3 =    l_cl_35_no3 + 1,    SP_CL_37_NO3    = 00099001
  integer, parameter :: L_CL_35_O =      l_cl_37_no3 + 1,    SP_CL_35_O      = 00051002
  integer, parameter :: L_CL_37_O =      l_cl_35_o + 1,      SP_CL_37_O      = 00053002
  integer, parameter :: L_CLO =          l_cl_37_o + 1,      SP_CLO          = 00051002
  integer, parameter :: L_CO =           l_clo + 1,          SP_CO           = 00028001
  integer, parameter :: L_CO_18 =        l_co + 1,           SP_CO_18        = 00030001
  integer, parameter :: L_COF2 =         l_co_18 + 1,        SP_COF2         = 00066001
  integer, parameter :: L_EXTINCTION =   l_cof2 + 1,         SP_EXTINCTION   = 00028965
  integer, parameter :: L_H2C_13_O =     l_extinction + 1,   SP_H2C_13_O     = 00031002
  integer, parameter :: L_H2CO =         l_h2c_13_o + 1,     SP_H2CO         = 00030004 
  integer, parameter :: L_H2CO_18 =      l_h2co + 1,         SP_H2CO_18      = 00032004
  integer, parameter :: L_H2O =          l_h2co_18 + 1,      SP_H2O          = 00018003
  integer, parameter :: L_H2O_R1A =      l_h2o + 1,          SP_H2O_R1A      = 00018997
  integer, parameter :: L_H2O_R1B =      l_h2o_r1a + 1,      SP_H2O_R1B      = 00018997 
  integer, parameter :: L_H2O_R2 =       l_h2o_r1b + 1,      SP_H2O_R2       = 00018997
  integer, parameter :: L_H2O_R3 =       l_h2o_r2 + 1,       SP_H2O_R3       = 00018997
  integer, parameter :: L_H2O_R4 =       l_h2o_r3 + 1,       SP_H2O_R4       = 00018997
  integer, parameter :: L_H2O_R5H =      l_h2o_r4 + 1,       SP_H2O_R5H      = 00018997
  integer, parameter :: L_H2O_R5V =      l_h2o_r5h + 1,      SP_H2O_R5V      = 00018997
  integer, parameter :: L_H2O2 =         l_h2o_r5v + 1,      SP_H2O2         = 00034004
  integer, parameter :: L_H2O_18 =       l_h2o2 + 1,         SP_H2O_18       = 00020003
  integer, parameter :: L_H2O_V2 =       l_h2o_18 + 1,       SP_H2O_V2       = 00018005
  integer, parameter :: L_HC_13_N =      l_h2o_v2 + 1,       SP_HC_13_N      = 00028002
  integer, parameter :: L_HCL =          l_hc_13_n + 1,      SP_HCL          = 00036001
  integer, parameter :: L_HCL_35 =       l_hcl + 1,          SP_HCL_35       = 00036001
  integer, parameter :: L_HCL_37 =       l_hcl_35 + 1,       SP_HCL_37       = 00038001
  integer, parameter :: L_HCN =          l_hcl_37 + 1,       SP_HCN          = 00027001
  integer, parameter :: L_HCN_15 =       l_hcn + 1,          SP_HCN_15       = 00028003
  integer, parameter :: L_HDO =          l_hcn_15 + 1,       SP_HDO          = 00019002
  integer, parameter :: L_HNO3 =         l_hdo + 1,          SP_HNO3         = 00063001
  integer, parameter :: L_HNO3_v5 =      l_hno3 + 1,         SP_HNO3_V5      = 00063006
  integer, parameter :: L_HNO3_v6 =      l_hno3_v5 + 1,      SP_HNO3_V6      = 00063004
  integer, parameter :: L_HNO3_v7 =      l_hno3_v6 + 1,      SP_HNO3_V7      = 00063002
  integer, parameter :: L_HNO3_v8 =      l_hno3_v7 + 1,      SP_HNO3_V8      = 00063005
  integer, parameter :: L_HNO3_v9 =      l_hno3_v8 + 1,      SP_HNO3_V9      = 00063003
  integer, parameter :: L_HO2 =          l_hno3_v9 + 1,      SP_HO2          = 00033001
  integer, parameter :: L_HOCL =         l_ho2 + 1,          SP_HOCL         = 00052006
  integer, parameter :: L_HOCL_35 =      l_hocl + 1,         SP_HOCL_35      = 00052006
  integer, parameter :: L_HOCL_37 =      l_hocl_35 + 1,      SP_HOCL_37      = 00054005
  integer, parameter :: L_LIQ_H2O =      l_hocl_37 + 1,      SP_LIQ_H2O      = 00018999
  integer, parameter :: L_N2 =           l_liq_h2o + 1,      SP_N2           = 00028964
  integer, parameter :: L_N2O =          l_n2 + 1,           SP_N2O          = 00044004
  integer, parameter :: L_NO =           l_n2o + 1,          SP_NO           = 00030008
  integer, parameter :: L_NO2 =          l_no + 1,           SP_NO2          = 00046006
  integer, parameter :: L_O2 =           l_no2 + 1,          SP_O2           = 00032001
  integer, parameter :: L_O2_V1 =        l_o2 + 1,           SP_O2_V1        = 00032002
  integer, parameter :: L_O3 =           l_o2_v1 + 1,        SP_O3           = 00048004
  integer, parameter :: L_O3_2V2 =       l_o3 + 1,           SP_O3_2V2       = 00048007
  integer, parameter :: L_O3_A_O18_V2 =  l_o3_2v2 + 1,       SP_O3_A_O18_V2  = 00050006
  integer, parameter :: L_O3_ASYM_O_17 = l_o3_a_o18_v2 + 1,  SP_O3_ASYM_O_17 = 00049002
  integer, parameter :: L_O3_ASYM_O_18 = l_o3_asym_o_17 + 1, SP_O3_ASYM_O_18 = 00050004
  integer, parameter :: L_O3_S_O18_V2 =  l_o3_asym_o_18 + 1, SP_O3_S_O18_V2  = 00050005
  integer, parameter :: L_O3_SYM_O_17 =  l_o3_s_o18_v2 + 1,  SP_O3_SYM_O_17  = 00049001
  integer, parameter :: L_O3_SYM_O_18 =  l_o3_sym_o_17 + 1,  SP_O3_SYM_O_18  = 00050003
  integer, parameter :: L_O3_V1_3 =      l_o3_sym_o_18 + 1,  SP_O3_V1_3      = 00048006
  integer, parameter :: L_O3_V2 =        l_o3_v1_3 + 1,      SP_O3_V2        = 00048005
  integer, parameter :: L_O_17_O =       l_o3_v2 + 1,        SP_O_17_O       = 00033002
  integer, parameter :: L_O_18_O =       l_o_17_o + 1,       SP_O_18_O       = 00034001
  integer, parameter :: L_OC_34_S =      l_o_18_o + 1,       SP_OC_34_S      = 00062001
  integer, parameter :: L_OCL_35_O =     l_oc_34_s + 1,      SP_OCL_35_O     = 00067001
  integer, parameter :: L_OCL_37_O =     l_ocl_35_o + 1,     SP_OCL_37_O     = 00069001
  integer, parameter :: L_OCS =          l_ocl_37_o + 1,     SP_OCS          = 00060001
  integer, parameter :: L_OH =           l_ocs + 1,          SP_OH           = 00017001
  integer, parameter :: L_SO2 =          l_oh + 1,           SP_SO2          = 00064002
  integer, parameter :: L_S_32_O2 =      l_so2 + 1,          SP_S_32_O2      = 00064002
  integer, parameter :: LAST_MOLECULE =  l_s_32_o2

! Mapping from the literals to the spec tags:
  integer :: SPEC_TAGS(FIRST_MOLECULE:LAST_MOLECULE)
  data spec_tags(l_air_cont)     / sp_air_cont /
  data spec_tags(l_bro)          / sp_bro /
  data spec_tags(l_br_79_o)      / sp_br_79_o /
  data spec_tags(l_br_81_o)      / sp_br_81_o /
  data spec_tags(l_c_13_o)       / sp_c_13_o /
  data spec_tags(l_ch3cl)        / sp_ch3cl /
  data spec_tags(l_ch3cl_35)     / sp_ch3cl_35 /
  data spec_tags(l_ch3cl_37)     / sp_ch3cl_37 /
  data spec_tags(l_ch3cn)        / sp_ch3cn /
  data spec_tags(l_cl_35_no3)    / sp_cl_35_no3 /
  data spec_tags(l_cl_37_no3)    / sp_cl_37_no3 /
  data spec_tags(l_cl_35_o)      / sp_cl_35_o /
  data spec_tags(l_cl_37_o)      / sp_cl_37_o /
  data spec_tags(l_clo)          / sp_clo /
  data spec_tags(l_co)           / sp_co /
  data spec_tags(l_co_18)        / sp_co_18 /
  data spec_tags(l_cof2)         / sp_cof2 /
  data spec_tags(l_extinction)   / sp_extinction /
  data spec_tags(l_h2c_13_o)     / sp_h2c_13_o /
  data spec_tags(l_h2co)         / sp_h2co /
  data spec_tags(l_h2co_18)      / sp_h2co_18 /
  data spec_tags(l_h2o)          / sp_h2o /
  data spec_tags(l_h2o_r1a)      / sp_h2o_r1a /
  data spec_tags(l_h2o_r1b)      / sp_h2o_r1b  /
  data spec_tags(l_h2o_r2)       / sp_h2o_r2 /
  data spec_tags(l_h2o_r3)       / sp_h2o_r3  /
  data spec_tags(l_h2o_r4)       / sp_h2o_r4 /
  data spec_tags(l_h2o_r5h)      / sp_h2o_r5h /
  data spec_tags(l_h2o_r5v)      / sp_h2o_r5v /
  data spec_tags(l_h2o2)         / sp_h2o2 /
  data spec_tags(l_h2o_18)       / sp_h2o_18 /
  data spec_tags(l_h2o_v2)       / sp_h2o_v2 /
  data spec_tags(l_hc_13_n)      / sp_hc_13_n /
  data spec_tags(l_hcl)          / sp_hcl /
  data spec_tags(l_hcl_35)       / sp_hcl_35 /
  data spec_tags(l_hcl_37)       / sp_hcl_37 /
  data spec_tags(l_hcn)          / sp_hcn /
  data spec_tags(l_hcn_15)       / sp_hcn_15 /
  data spec_tags(l_hdo)          / sp_hdo /
  data spec_tags(l_hno3)         / sp_hno3 /
  data spec_tags(l_hno3_v5)      / sp_hno3_v5 /
  data spec_tags(l_hno3_v6)      / sp_hno3_v6 /
  data spec_tags(l_hno3_v7)      / sp_hno3_v7 /
  data spec_tags(l_hno3_v8)      / sp_hno3_v8 /
  data spec_tags(l_hno3_v9)      / sp_hno3_v9 /
  data spec_tags(l_ho2)          / sp_ho2 /
  data spec_tags(l_hocl)         / sp_hocl /
  data spec_tags(l_hocl_35)      / sp_hocl_35 /
  data spec_tags(l_hocl_37)      / sp_hocl_37 /
  data spec_tags(l_liq_h2o)      / sp_liq_h2o /
  data spec_tags(l_n2)           / sp_n2 /
  data spec_tags(l_n2o)          / sp_n2o /
  data spec_tags(l_no)           / sp_no /
  data spec_tags(l_no2)          / sp_no2 /
  data spec_tags(l_o2)           / sp_o2 /
  data spec_tags(l_o2_v1)        / sp_o2_v1 /
  data spec_tags(l_o3)           / sp_o3 /
  data spec_tags(l_o3_2v2)       / sp_o3_2v2 /
  data spec_tags(l_o3_a_o18_v2)  / sp_o3_a_o18_v2 /
  data spec_tags(l_o3_asym_o_17) / sp_o3_asym_o_17 /
  data spec_tags(l_o3_asym_o_18) / sp_o3_asym_o_18 /
  data spec_tags(l_o3_s_o18_v2)  / sp_o3_s_o18_v2 /
  data spec_tags(l_o3_sym_o_17)  / sp_o3_sym_o_17 /
  data spec_tags(l_o3_sym_o_18)  / sp_o3_sym_o_18 /
  data spec_tags(l_o3_v1_3)      / sp_o3_v1_3 /
  data spec_tags(l_o3_v2)        / sp_o3_v2 /
  data spec_tags(l_o_17_o)       / sp_o_17_o /
  data spec_tags(l_o_18_o)       / sp_o_18_o /
  data spec_tags(l_oc_34_s)      / sp_oc_34_s /
  data spec_tags(l_ocl_35_o)     / sp_ocl_35_o /
  data spec_tags(l_ocl_37_o)     / sp_ocl_37_o /
  data spec_tags(l_ocs)          / sp_ocs /
  data spec_tags(l_oh)           / sp_oh /
  data spec_tags(l_so2)          / sp_so2 /
  data spec_tags(l_s_32_o2)      / sp_s_32_o2 /     

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
    lit_indices(l_ch3cl) =        add_ident ( 'ch3cl' )
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
    lit_indices(l_extinction) =   add_ident ( 'extinction' )
    lit_indices(l_h2c_13_o) =     add_ident ( 'h2c_13_o' )
    lit_indices(l_h2co) =         add_ident ( 'h2co' )
    lit_indices(l_h2co_18) =      add_ident ( 'h2co_18' )
    lit_indices(l_h2o) =          add_ident ( 'h2o' )
    lit_indices(l_h2o_r1a) =      add_ident ( 'h2o_r1a' )
    lit_indices(l_h2o_r1b) =      add_ident ( 'h2o_r1b' )
    lit_indices(l_h2o_r2) =       add_ident ( 'h2o_r2' )
    lit_indices(l_h2o_r3) =       add_ident ( 'h2o_r3' )
    lit_indices(l_h2o_r4) =       add_ident ( 'h2o_r4' )
    lit_indices(l_h2o_r5h) =      add_ident ( 'h2o_r5h' )
    lit_indices(l_h2o_r5v) =      add_ident ( 'h2o_r5v' )
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
    lit_indices(l_hno3_v5) =      add_ident ( 'hno3_v5' )
    lit_indices(l_hno3_v6) =      add_ident ( 'hno3_v6' )
    lit_indices(l_hno3_v7) =      add_ident ( 'hno3_v7' )
    lit_indices(l_hno3_v8) =      add_ident ( 'hno3_v8' )
    lit_indices(l_hno3_v9) =      add_ident ( 'hno3_v9' )
    lit_indices(l_ho2) =          add_ident ( 'ho2' )
    lit_indices(l_hocl) =         add_ident ( 'hocl' )
    lit_indices(l_hocl_35) =      add_ident ( 'hocl_35' )
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
    lit_indices(l_s_32_o2) =      add_ident ( 's_32_o2' )

    ! Create the type tree for the molecule type
    call make_tree ( (/ &
      begin, t+t_molecule, l+(/ (i,i=first_molecule, last_molecule) /), &
             n+n_dt_def /) )

  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine INIT_MOLECULES

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MOLECULES

! $Log$
! Revision 2.18  2002/12/13 02:04:59  vsnyder
! Give names to all of the spectags
!
! Revision 2.17  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.16  2002/07/30 20:12:13  livesey
! Changed the radiometer specific h2os to 18997 spectag
!
! Revision 2.15  2002/07/24 22:26:45  livesey
! Added the radiometer specific H2O's and a definitive CH3Cl
!
! Revision 2.14  2002/07/18 22:05:26  vsnyder
! Alphabetize L_CL_35, L_CL_37
!
! Revision 2.13  2001/11/08 00:11:44  livesey
! Added extinction as a molecule
!
! Revision 2.12  2001/11/02 18:54:06  livesey
! Added hno3 excited states
!
! Revision 2.11  2001/10/09 22:39:52  livesey
! Added S_32_O2
!
! Revision 2.10  2001/05/11 00:38:46  livesey
! Added hocl_35
!
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

