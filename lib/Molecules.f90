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

! These particular radiometer-dependent molecules are hand-coded
! instead of being generated automatically
  integer, parameter :: FIRST_MOLECULE = last_intrinsic_lit + 1
  integer, parameter :: L_H2O_R1A = first_molecule
  integer, parameter :: L_H2O_R1B =      l_h2o_r1a + 1
  integer, parameter :: L_H2O_R2 =       l_h2o_r1b + 1
  integer, parameter :: L_H2O_R3 =       l_h2o_r2 + 1
  integer, parameter :: L_H2O_R4 =       l_h2o_r3 + 1
  integer, parameter :: L_H2O_R5H =      l_h2o_r4 + 1
  integer, parameter :: L_H2O_R5V =      l_h2o_r5h + 1
  integer, parameter :: LAST_RXX_MOLECULE =  l_h2o_r5v

! Don't edit the following file directly--it is generated automatically
! based on the file sps_cross_ref_table.txt
  include 'mol_parm.f9h'

contains ! =====     Public procedures     =============================
! -----------------------------------------------  INIT_MOLECULES  -----
  subroutine Init_Molecules ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    use Functions, only: Init_Functions
    ! "use Tree" is here because "make depends" can't see it in make_tree
    ! (because of the "include"):
    use TREE, only:
    use TREE_TYPES, only: N_DT_DEF

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    integer :: I    ! Loop inductor in an array constructor

    ! Initialize the intrinsic types

    call init_functions ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices )

    ! Put type names into the symbol table:

    data_type_indices(t_molecule) =        add_ident ( 'molecule' )

    ! Put literals into the symbol table:

    ! Don't edit the following file directly--it is generated automatically
    ! based on the file sps_cross_ref_table.txt
    include 'mol_add.f9h'

    lit_indices(l_h2o_r1a) =      add_ident ( 'h2o_r1a' )
    lit_indices(l_h2o_r1b) =      add_ident ( 'h2o_r1b' )
    lit_indices(l_h2o_r2) =       add_ident ( 'h2o_r2' )
    lit_indices(l_h2o_r3) =       add_ident ( 'h2o_r3' )
    lit_indices(l_h2o_r4) =       add_ident ( 'h2o_r4' )
    lit_indices(l_h2o_r5h) =      add_ident ( 'h2o_r5h' )
    lit_indices(l_h2o_r5v) =      add_ident ( 'h2o_r5v' )

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
! Revision 2.21  2004/05/29 02:43:11  vsnyder
! Rearrange function definition stuff
!
! Revision 2.20  2003/05/19 19:37:41  vsnyder
! Remove USE's for unreferenced names
!
! Revision 2.19  2003/05/17 00:03:15  pwagner
! Made to work with mol_... files created by init_gen
!
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

