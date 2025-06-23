! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Init_MLSSignals_m

  use intrinsic, only: add_ident, begin, d, f, field_first, g, l, &
    & last_intrinsic_lit, last_intrinsic_spec, n, nadp, ndp, np, nr, p, s, &
    & spec_first, t, t_boolean, t_numeric, t_numeric_range, t_polarization, &
    & t_string, z

  use intrinsic, only: data_type_indices, du, field_indices, func_indices, &
    & lit_indices, parm_indices, phyq_dimensionless, phyq_frequency, &
    & section_indices, spec_indices

  use intrinsic, only: l_asmls, l_emls, l_umls, l_xptl1
  use molecules, only: init_molecules, last_molecule, last_molecule_type
  implicit none

  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Fields used in signal specifications:
  integer, parameter :: F_Aura              = field_First
  integer, parameter :: F_band              = f_Aura + 1
  integer, parameter :: F_centerFrequency   = f_band + 1
  integer, parameter :: F_channel           = f_centerFrequency + 1
  integer, parameter :: F_channels          = f_channel + 1
  integer, parameter :: F_dacs              = f_channels + 1
  integer, parameter :: F_deferred          = f_dacs + 1
  integer, parameter :: F_direction         = f_deferred + 1
  integer, parameter :: F_first             = f_direction + 1
  integer, parameter :: F_instrument        = f_first + 1
  integer, parameter :: F_last              = f_instrument + 1
  integer, parameter :: F_lo                = f_last + 1
  integer, parameter :: F_module            = f_lo + 1
  integer, parameter :: F_polarization      = f_module + 1
  integer, parameter :: F_radiometer        = f_polarization + 1
  integer, parameter :: F_singleSideband    = f_radiometer + 1
  integer, parameter :: F_spacecraft        = f_singleSideband + 1
  integer, parameter :: F_spectrometer      = f_spacecraft + 1
  integer, parameter :: F_spectrometerType  = f_spectrometer + 1
  integer, parameter :: F_start             = f_spectrometerType + 1
  integer, parameter :: F_step              = f_start + 1
  integer, parameter :: F_suffix            = f_step + 1
  integer, parameter :: F_supportedModule   = f_suffix + 1
  integer, parameter :: F_switch            = f_supportedModule + 1
  integer, parameter :: F_width             = f_switch + 1
  integer, parameter :: Last_Signal_Field   = f_width

  ! Literals used in signal specifications:
  integer, parameter :: Last_Signal_Lit     = last_molecule

  ! Enumeration types:
  integer, parameter :: T_InstrumentType    = last_molecule_type+1
  integer, parameter :: Last_signal_type    = t_instrumenttype

  ! Signal specifications:
  integer, parameter :: S_band              = t_instrumenttype + 1
  integer, parameter :: S_module            = s_band + 1
  integer, parameter :: S_radiometer        = s_module + 1
  integer, parameter :: S_signal            = s_radiometer + 1
  integer, parameter :: S_spectrometerType  = s_signal + 1
  integer, parameter :: Last_Signal_Spec    = s_spectrometerType

  ! The MLSSignals section is NOT defined here, because it appears
  ! in the section ordering requirements array in init_tables_module.

contains

  ! --------------------------------------------  Init_MLSSignals  -----
  subroutine Init_MLSSignals ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    ! "use Tree" really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only:
    use TREE_TYPES, only: N_FIELD_SPEC, N_FIELD_TYPE, N_SPEC_DEF, N_DT_DEF

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    call init_molecules ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices )

    ! Put field names into the symbol table
    field_indices(f_Aura) =                add_ident ( 'Aura' )
    field_indices(f_band) =                add_ident ( 'band' )
    field_indices(f_centerFrequency) =     add_ident ( 'centerFrequency' )
    field_indices(f_channel) =             add_ident ( 'channel' )
    field_indices(f_channels) =            add_ident ( 'channels' )
    field_indices(f_dacs) =                add_ident ( 'dacs' )
    field_indices(f_deferred) =            add_ident ( 'deferred' )
    field_indices(f_direction) =           add_ident ( 'direction' )
    field_indices(f_first) =               add_ident ( 'first' )
    field_indices(f_instrument) =          add_ident ( 'instrument' )
    field_indices(f_last) =                add_ident ( 'last' )
    field_indices(f_lo) =                  add_ident ( 'lo' )
    field_indices(f_module) =              add_ident ( 'module' )
    field_indices(f_polarization) =        add_ident ( 'polarization' )
    field_indices(f_radiometer) =          add_ident ( 'radiometer' )
    field_indices(f_singleSideband) =      add_ident ( 'singleSideband' )
    field_indices(f_spacecraft) =          add_ident ( 'spacecraft' )
    field_indices(f_spectrometer) =        add_ident ( 'spectrometer' )
    field_indices(f_spectrometerType) =    add_ident ( 'spectrometerType' )
    field_indices(f_start) =               add_ident ( 'start' )
    field_indices(f_step) =                add_ident ( 'step' )
    field_indices(f_suffix) =              add_ident ( 'suffix' )
    field_indices(f_supportedModule) =     add_ident ( 'supportedModule' )
    field_indices(f_switch) =              add_ident ( 'switch' )
    field_indices(f_width) =               add_ident ( 'width' )

    ! Put enumeration type names into the symbol table
    data_type_indices(t_instrumenttype) =  add_ident ( 'instrumenttype' )
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
             begin, t+t_instrumentType, &
             l+l_asmls, l+l_emls, l+l_umls, l+l_xptl1, n+n_dt_def /) )
    call make_tree ( (/ &
      begin, s+s_module, &          ! Must be after module
             begin, f+f_Aura, t+t_boolean, n+n_field_type, &
             begin, f+f_instrument, t+t_instrumentType, n+n_field_type, &
             begin, f+f_spacecraft, t+t_boolean, n+n_field_type, &
             begin, f+f_supportedModule, s+s_module, n+n_field_spec, &
             np+n_spec_def, &
      begin, s+s_radiometer, &          ! Must be after module
             begin, f+f_lo, t+t_numeric, nr+n_field_type+du*PHYQ_FREQUENCY, &
             begin, f+f_suffix, t+t_string, nr+n_field_type, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_polarization, t+t_polarization, n+n_field_type, &
             begin, f+f_singleSideband, t+t_numeric, nr+n_field_type+du*PHYQ_DIMENSIONLESS, &
             ndp+n_spec_def, &
      begin, s+s_spectrometerType, &
             begin, f+f_channels, t+t_numeric_range, n+n_field_type+du*PHYQ_FREQUENCY, &
             begin, f+f_dacs, t+t_boolean, n+n_field_type, &
             begin, f+f_deferred, t+t_boolean, n+n_field_type, &
             begin, f+f_first, t+t_numeric, n+n_field_type+du*PHYQ_DIMENSIONLESS, &
             begin, f+f_last, t+t_numeric, n+n_field_type+du*PHYQ_DIMENSIONLESS, &
             begin, f+f_start, t+t_numeric, n+n_field_type+du*PHYQ_FREQUENCY, &
             begin, f+f_step, t+t_numeric, n+n_field_type+du*PHYQ_FREQUENCY, &
             begin, f+f_width, t+t_numeric, n+n_field_type+du*PHYQ_FREQUENCY, &
             ndp+n_spec_def, &
      begin, s+s_band, &                ! Must be after radiometer and spectrometerType
             begin, f+f_suffix, t+t_string, n+n_field_type, &
             begin, f+f_spectrometerType, s+s_spectrometerType, nr+n_field_spec, &
             begin, f+f_radiometer, s+s_radiometer, n+n_field_spec, &
             begin, f+f_centerfrequency, t+t_numeric, n+n_field_type+du*PHYQ_FREQUENCY, &
             ndp+n_spec_def, &
      begin, s+s_signal, &              ! Must be after band
             begin, f+f_band, s+s_band, nr+n_field_spec, &
             begin, f+f_channels, t+t_numeric_range, n+n_field_type+du*PHYQ_FREQUENCY, &
             begin, f+f_direction, t+t_numeric, nr+n_field_type+du*PHYQ_DIMENSIONLESS, &
             begin, f+f_radiometer, s+s_radiometer, n+n_field_spec+du*PHYQ_DIMENSIONLESS, &
             begin, f+f_module, s+s_module, n+n_field_spec, &
             begin, f+f_spectrometer, t+t_numeric, nr+n_field_type+du*PHYQ_DIMENSIONLESS, &
             begin, f+f_switch, t+t_numeric, nr+n_field_type+du*PHYQ_DIMENSIONLESS, &
             ndp+n_spec_def /) )

  contains
    ! --------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine Init_MLSSignals

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Init_MLSSignals_m

! $Log$
! Revision 2.36  2020/07/05 20:20:47  vsnyder
! Remove frequency field because it's not used
!
! Revision 2.35  2018/02/27 00:49:43  livesey
! Added the supportedModule functionality to support ASMLS
!
! Revision 2.34  2017/09/15 15:44:18  livesey
! Updated to allow modules to be defferred until signal definition
!
! Revision 2.33  2016/12/15 18:23:27  pwagner
! Added asmls instrument
!
! Revision 2.32  2014/04/22 00:07:48  vsnyder
! Put InstrumentType's index in the range of type indices
!
! Revision 2.31  2013/11/06 01:46:30  pwagner
! May read instrument field of module; e.g. emls
!
! Revision 2.30  2013/09/04 00:00:59  pwagner
! Comments preceding use TREE made more uniform across modules
!
! Revision 2.29  2013/08/23 23:23:11  pwagner
! Added Aura field to Module spec
!
! Revision 2.28  2011/01/29 00:46:42  vsnyder
! Add units checking
!
! Revision 2.27  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.26  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.25  2004/05/29 02:42:59  vsnyder
! Rearrange function definition stuff
!
! Revision 2.24  2004/01/16 21:37:23  livesey
! Added the ability to defer the connection between bands and radiometers
! until you define the signal.  This is to support some SMLS related
! research work.
!
! Revision 2.23  2003/08/16 01:14:03  vsnyder
! Add optional 'polarization' field to 'radiometer' spec
!
! Revision 2.22  2003/07/23 18:02:30  livesey
! Reverted DACS back to lower case
!
! Revision 2.21  2003/07/23 07:17:46  livesey
! Capitalized DACS
!
! Revision 2.20  2003/07/18 20:23:26  livesey
! Added DACS flags etc.
!
! Revision 2.19  2003/05/16 02:44:18  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.18  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.17  2002/05/14 22:31:47  livesey
! Added singleSideband
!
! Revision 2.16  2002/05/03 22:38:55  livesey
! Added direction field
!
! Revision 2.15  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.14  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.13  2001/04/23 20:57:42  vsnyder
! Move the first spec (time) to 'intrinsic'
!
! Revision 2.12  2001/04/11 20:19:27  vsnyder
! Undo changes to 'deferred'
!
! Revision 2.11  2001/04/11 18:31:04  vsnyder
! Change 'deferred' from boolean to numeric
!
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
