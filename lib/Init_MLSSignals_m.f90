module Init_MLSSignals_m

  use INTRINSIC, only: Add_Ident, Begin, D, F, L, Last_Intrinsic_Lit, N, &
    & NADP, NDP, NP, NR, P, S, T, T_Boolean, T_Numeric, &
    & T_Numeric_Range, T_String, Z

  use MOLECULES, only: Init_Molecules, Last_Molecule, Last_Molecule_Type
  implicit NONE

  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Types used in signal specifications:
  integer, parameter :: Last_signal_type = last_molecule_type

  ! Fields used in signal specifications:
  integer, parameter :: Field_First = 1
  integer, parameter :: F_band = field_First
  integer, parameter :: F_centerFrequency   = f_band + 1
  integer, parameter :: F_channel           = f_centerFrequency + 1
  integer, parameter :: F_channels          = f_channel + 1
  integer, parameter :: F_deferred          = f_channels + 1
  integer, parameter :: F_first             = f_deferred + 1
  integer, parameter :: F_frequency         = f_first + 1
  integer, parameter :: F_last              = f_frequency + 1
  integer, parameter :: F_lo                = f_last + 1
  integer, parameter :: F_module            = f_lo + 1
  integer, parameter :: F_radiometer        = f_module + 1
  integer, parameter :: F_spacecraft        = f_radiometer + 1
  integer, parameter :: F_spectrometer      = f_spacecraft + 1
  integer, parameter :: F_spectrometerType  = f_spectrometer + 1
  integer, parameter :: F_start             = f_spectrometerType + 1
  integer, parameter :: F_step              = f_start + 1
  integer, parameter :: F_suffix            = f_step + 1
  integer, parameter :: F_switch            = f_suffix + 1
  integer, parameter :: F_width             = f_switch + 1
  integer, parameter :: Last_Signal_Field   = f_width

  ! Literals used in signal specifications:
  integer, parameter :: Last_Signal_Lit     = last_molecule

  ! Signal specifications:
  integer, parameter :: Spec_First = 1
  integer, parameter :: S_band             = spec_First
  integer, parameter :: S_module           = s_band + 1
  integer, parameter :: S_radiometer       = s_module + 1
  integer, parameter :: S_signal           = s_radiometer + 1
  integer, parameter :: S_spectrometerType = s_signal + 1
  integer, parameter :: Last_Signal_Spec = s_spectrometerType

  ! The MLSSignals section is NOT defined here, because it appears
  ! in the section ordering requirements array in init_tables_module.

contains
  ! --------------------------------------------  Init_MLSSignals  -----
  subroutine Init_MLSSignals ( Data_Type_Indices, Field_Indices, Lit_Indices, &
    & Parm_Indices, Section_Indices, Spec_Indices )

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only: BUILD_TREE, PUSH_PSEUDO_TERMINAL
    use TREE_TYPES, only: N_DT_DEF, N_FIELD_SPEC, N_FIELD_TYPE, N_SPEC_DEF

    integer, intent(inout) :: Data_Type_Indices(:)
    integer, intent(inout) :: Lit_Indices(:)
    integer, intent(inout) :: Field_Indices(field_First:last_Signal_Field)
    integer, intent(inout) :: Parm_Indices(:)
    integer, intent(inout) :: Section_Indices(:)
    integer, intent(inout) :: Spec_Indices(spec_First:last_Signal_Spec)

    call init_molecules ( data_type_indices, field_indices, lit_indices, &
      & parm_indices, section_indices, spec_indices )

    ! Put field names into the symbol table
    field_indices(f_band) =                add_ident ( 'band' )
    field_indices(f_centerFrequency) =     add_ident ( 'centerFrequency' )
    field_indices(f_channel) =             add_ident ( 'channel' )
    field_indices(f_channels) =            add_ident ( 'channels' )
    field_indices(f_deferred) =            add_ident ( 'deferred' )
    field_indices(f_first) =               add_ident ( 'first' )
    field_indices(f_frequency) =           add_ident ( 'frequency' )
    field_indices(f_last) =                add_ident ( 'last' )
    field_indices(f_lo) =                  add_ident ( 'lo' )
    field_indices(f_module) =              add_ident ( 'module' )
    field_indices(f_radiometer) =          add_ident ( 'radiometer' )
    field_indices(f_spacecraft) =          add_ident ( 'spacecraft' )
    field_indices(f_spectrometer) =        add_ident ( 'spectrometer' )
    field_indices(f_spectrometerType) =    add_ident ( 'spectrometerType' )
    field_indices(f_start) =               add_ident ( 'start' )
    field_indices(f_step) =                add_ident ( 'step' )
    field_indices(f_suffix) =              add_ident ( 'suffix' )
    field_indices(f_switch) =              add_ident ( 'switch' )
    field_indices(f_width) =               add_ident ( 'width' )
    ! Put spec names into the symbol table
    spec_indices(s_band) =                 add_ident ( 'band' )
    spec_indices(s_module) =               add_ident ( 'module' )
    spec_indices(s_radiometer) =           add_ident ( 'radiometer' )
    spec_indices(s_signal) =               add_ident ( 'signal' )
    spec_indices(s_spectrometerType) =     add_ident ( 'spectrometerType' )

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
      begin, s+s_module, &
             begin, f+f_spacecraft, t+t_boolean, n+n_field_type, &
             np+n_spec_def, &
      begin, s+s_radiometer, &          ! Must be after module
             begin, f+f_lo, t+t_numeric, n+n_field_type, &
             begin, f+f_suffix, t+t_string, n+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             nadp+n_spec_def, &
      begin, s+s_spectrometerType, &
             begin, f+f_channels, t+t_numeric_range, n+n_field_type, &
             begin, f+f_deferred, t+t_boolean, n+n_field_type, &
             begin, f+f_first, t+t_numeric, n+n_field_type, &
             begin, f+f_last, t+t_numeric, n+n_field_type, &
             begin, f+f_start, t+t_numeric, n+n_field_type, &
             begin, f+f_step, t+t_numeric, n+n_field_type, &
             begin, f+f_width, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_band, &                ! Must be after radiometer and spectrometerType
             begin, f+f_suffix, t+t_string, n+n_field_type, &
             begin, f+f_spectrometerType, s+s_spectrometerType, nr+n_field_spec, &
             begin, f+f_radiometer, s+s_radiometer, nr+n_field_spec, &
             begin, f+f_centerfrequency, t+t_numeric, n+n_field_type, &
             ndp+n_spec_def, &
      begin, s+s_signal, &              ! Must be after band
             begin, f+f_band, s+s_band, nr+n_field_spec, &
             begin, f+f_channels, t+t_numeric_range, n+n_field_type, &
             begin, f+f_spectrometer, t+t_numeric, nr+n_field_type, &
             begin, f+f_switch, t+t_numeric, nr+n_field_type, &
             ndp+n_spec_def /) )

  contains
    ! --------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine Init_MLSSignals

end module Init_MLSSignals_m

! $Log$
! Revision 2.10  2001/04/10 18:51:25  vsnyder
! Finish removing sideband stuff
!
! Revision 2.9  2001/04/10 17:59:53  vsnyder
! Remove sideband field from signal
!
! Revision 2.8  2001/04/04 17:56:42  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 2.7  2001/04/03 19:09:12  vsnyder
! Change the order of initialization to intrinsic, Molecules, MLSSignals.
! Use the revised make_tree.f9h, which requires revision of init...
! calling sequences.
!
! Revision 2.6  2001/03/28 19:50:43  vsnyder
! Remove frequencies and widths fields
!
! Revision 2.5  2001/03/16 02:00:40  vsnyder
! Add support for literals to Make_Tree (duh!)
!
! Revision 2.4  2001/03/16 01:02:32  vsnyder
! ... Including the Log at the end.
!
