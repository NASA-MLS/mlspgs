! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Init_Spectroscopy_m

  use INTRINSIC, only: Add_Ident, Begin, F, L, N, NR, NDP, NADP, P, S, T, &
    & T_Numeric, Z, T_String
  use Init_MLSSignals_m, only: Field_First, Init_MLSSignals, &
    & Last_Signal_Field, Last_Signal_Lit, Last_Signal_Spec, Last_Signal_Type, &
    & Spec_First
  use INTRINSIC, only: DATA_TYPE_INDICES, FIELD_INDICES, &
    & LIT_INDICES, PARM_INDICES, SECTION_INDICES, SPEC_INDICES
  use Molecules, only: T_Molecule

  implicit NONE

  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Types used in spectroscopy specifications (there aren't any):
  integer, parameter :: Last_Spectroscopy_Type = last_signal_type

  ! Literals
  integer, parameter :: Last_Spectroscopy_Lit = last_Signal_Lit

  ! Fields used in spectroscopy specifications:
  integer, parameter :: F_delta      = last_signal_field + 1
  integer, parameter :: F_el         = f_delta + 1
  integer, parameter :: F_gamma      = f_el + 1
  integer, parameter :: F_lines      = f_gamma + 1
  integer, parameter :: F_molecule   = f_lines + 1
  integer, parameter :: F_n          = f_molecule + 1
  integer, parameter :: F_n1         = f_n + 1
  integer, parameter :: F_n2         = f_n1 + 1
  integer, parameter :: F_ns         = f_n2 + 1
  integer, parameter :: F_ps         = f_ns + 1
  integer, parameter :: F_qlog       = f_ps + 1
  integer, parameter :: F_str        = f_qlog + 1
  integer, parameter :: F_v0         = f_str + 1
  integer, parameter :: F_w          = f_v0 + 1
  integer, parameter :: F_emlsSignals  = f_w + 1
  integer, parameter :: F_umlsSignals  = f_emlsSignals + 1
  integer, parameter :: Last_Spectroscopy_Field = f_umlsSignals

  ! Spectroscopy specifications:
  integer, parameter :: S_Line     = last_signal_spec + 1
  integer, parameter :: S_Spectra  = s_line + 1
  integer, parameter :: Last_Spectroscopy_Spec = s_spectra

  ! The Spectroscopy section is NOT defined here, because it appears
  ! in the section ordering requirements array in init_tables_module.

contains
  ! ------------------------------------------  Init_Spectroscopy  -----
  subroutine Init_Spectroscopy ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
    use TREE_TYPES, only: N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, N_SPEC_DEF

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    call init_MLSSignals ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices )

    ! Put field names into the symbol table
    field_indices(f_delta)      = add_ident ( 'delta' )
    field_indices(f_el)         = add_ident ( 'el' )
    field_indices(f_gamma)      = add_ident ( 'gamma' )
    field_indices(f_lines)      = add_ident ( 'lines' )
    field_indices(f_molecule)   = add_ident ( 'molecule' )
    field_indices(f_n)          = add_ident ( 'n' )
    field_indices(f_n1)         = add_ident ( 'n1' )
    field_indices(f_n2)         = add_ident ( 'n2' )
    field_indices(f_ns)         = add_ident ( 'ns' )
    field_indices(f_ps)         = add_ident ( 'ps' )
    field_indices(f_qlog)       = add_ident ( 'qlog' )
    field_indices(f_str)        = add_ident ( 'str' )
    field_indices(f_v0)         = add_ident ( 'v0' )
    field_indices(f_w)          = add_ident ( 'w' )
    field_indices(f_emlsSignals)  = add_ident ( 'emlsSignals' )
    field_indices(f_umlsSignals)  = add_ident ( 'umlsSignals' )

    ! Put spec names into the symbol table
    spec_indices(s_line)    = add_ident ( 'line' )
    spec_indices(s_spectra) = add_ident ( 'spectra' )

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
      begin, s+s_line, &
             begin, f+f_delta, t+t_numeric, nr+n_field_type, &
             begin, f+f_el, t+t_numeric, nr+n_field_type, &
             begin, f+f_gamma, t+t_numeric, nr+n_field_type, &
             begin, f+f_n, t+t_numeric, nr+n_field_type, &
             begin, f+f_n1, t+t_numeric, nr+n_field_type, &
             begin, f+f_n2, t+t_numeric, nr+n_field_type, &
             begin, f+f_ns, t+t_numeric, nr+n_field_type, &
             begin, f+f_ps, t+t_numeric, nr+n_field_type, &
             begin, f+f_str, t+t_numeric, nr+n_field_type, &
             begin, f+f_v0, t+t_numeric, nr+n_field_type, &
             begin, f+f_w, t+t_numeric, nr+n_field_type, &
             begin, f+f_emlsSignals, t+t_string, n+n_field_type, &
             begin, f+f_umlsSignals, t+t_string, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_spectra, & ! Must be AFTER S_Line
             begin, f+f_lines, s+s_line, n+n_field_spec, &
             begin, f+f_molecule, t+t_molecule, nr+n_field_type, &
             begin, f+f_qlog, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def /) )

  contains
    ! --------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine Init_Spectroscopy
end module Init_Spectroscopy_m

! $Log$
! Revision 2.2  2001/09/18 01:25:48  livesey
! Changed emls/umls bands to emls/umls signals
!
! Revision 2.1  2001/09/18 00:08:17  livesey
! Added the bands information stuff
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.5  2001/06/07 23:30:33  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.3  2001/04/04 17:59:42  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 1.2  2001/04/04 02:10:06  vsnyder
! Repair a literal name
!
! Revision 1.1  2001/04/03 19:41:40  vsnyder
! Initial commit
!
