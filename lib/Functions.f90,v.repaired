! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Functions

! Provide parameters to identify functions.
! Put names of functions into the declaration table, with their
! corresponding parameter in the "units" field.

! We may eventually want to put a tree for each function into init_tables
! to specify the allowed argument types and the result type, and put some
! checking in tree_checker.

  implicit NONE
  public

  integer, parameter :: Function_First    = 1
  integer, parameter :: F_Cholesky        = Function_First
  integer, parameter :: F_ClearLower      = F_Cholesky + 1
  integer, parameter :: F_Difference      = F_ClearLower + 1  ! Set difference
  integer, parameter :: F_Exp             = F_Difference + 1
  integer, parameter :: F_GetDiagonal     = F_Exp + 1
  integer, parameter :: F_Intersection    = F_GetDiagonal + 1 ! Set intersection
  integer, parameter :: F_Invert          = F_Intersection + 1
  integer, parameter :: F_Ln              = F_Invert + 1
  integer, parameter :: F_Log             = F_ln + 1
  integer, parameter :: F_Log10           = F_Log + 1
  integer, parameter :: F_Mod             = F_Log10 + 1
  integer, parameter :: F_SQRT            = F_Mod + 1
  integer, parameter :: F_Transpose       = F_SQRT + 1
  integer, parameter :: F_Union           = F_Transpose + 1
  integer, parameter :: F_Without         = F_Union + 1
  integer, parameter :: F_XTX             = F_Without + 1
  integer, parameter :: Function_Last     = F_XTX

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Init_Functions ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    ! "use Tree" really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):
    use TREE, only:
    use INTRINSIC, only: BEGIN, ADD_IDENT, D, DATA_TYPE_INDICES, &
      & F, FIELD_INDICES, FUNC_INDICES, &
      & G, INIT_INTRINSIC, L, LIT_INDICES, &
      & N, P, PARM_INDICES, S, SPEC_INDICES, SECTION_INDICES, T, T_NUMERIC, Z
    use TREE_TYPES, only: N_Dot, N_ARG_DEF, N_FUNC_DEF, N_Or

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    integer, parameter :: Num_or_Dot(5) = &
      & (/ begin, t+t_numeric, begin, n+n_dot, n+n_or /)

    ! Initialize the intrinsic types

    call init_intrinsic ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices, &
      & function_last )

    func_indices(f_cholesky) =        add_ident ( 'cholesky' )
    func_indices(f_clearLower) =      add_ident ( 'clearLower' )
    func_indices(f_difference) =      add_ident ( 'difference' )
    func_indices(f_exp) =             add_ident ( 'exp' )
    func_indices(f_getDiagonal) =     add_ident ( 'getDiagonal' )
    func_indices(f_intersection) =    add_ident ( 'intersection' )
    func_indices(f_invert) =          add_ident ( 'invert' )
    func_indices(f_ln) =              add_ident ( 'ln' )
    func_indices(f_log) =             add_ident ( 'log' )
    func_indices(f_log10) =           add_ident ( 'log10' )
    func_indices(f_mod) =             add_ident ( 'mod' )
    func_indices(f_sqrt) =            add_ident ( 'sqrt' )
    func_indices(f_transpose) =       add_ident ( 'transpose' )
    func_indices(f_union) =           add_ident ( 'union' )
    func_indices(f_without) =         add_ident ( 'without' )
    func_indices(f_xtx) =             add_ident ( 'xtx' )

    ! Define the functions and their arguments.  These are represented by
    ! trees of the form
    !  < n_func_def funtion_name
    !              < n_arg_def t_type ... t_type > ... >
    call make_tree ( (/ &
      !??? To get automatic type checking for f_cholesky etc., it is
      !??? probably necessary to do init_functions after init_tables.
      !??? OTOH, we could do some here, and some later, like we do lits.
      begin, g+f_cholesky, n+n_func_def, &
      begin, g+f_clearLower, n+n_func_def, &
      begin, g+f_difference, n+d+n_func_def, &
      begin, g+f_exp, &
             begin, num_or_dot, n+n_arg_def, n+n_func_def, &
      begin, g+f_getDiagonal, n+n_func_def, &
      begin, g+f_intersection, n+d+n_func_def, &
      begin, g+f_invert, n+n_func_def, &
      begin, g+f_ln, &
             begin, num_or_dot, n+n_arg_def, n+n_func_def, &
      begin, g+f_log, &
             begin, num_or_dot, n+n_arg_def, n+n_func_def, &
      begin, g+f_log10, &
             begin, num_or_dot, n+n_arg_def, n+n_func_def, &
      begin, g+f_mod, &
             begin, num_or_dot, num_or_dot, n+n_arg_def, n+n_func_def, &
      !??? Automatic type checking for f_sqrt may be difficult, given that
      !??? we want to allow numbers or matrices.  It may be necessary either
      !??? to give up on that, or to do init_functions after init_tables.
      begin, g+f_sqrt, &
             begin, num_or_dot, n+n_arg_def, n+n_func_def, &
      begin, g+f_transpose, n+n_func_def, &
      begin, g+f_union, n+d+n_func_def, &
      begin, g+f_without, n+d+n_func_def, &
      begin, g+f_xtx, n+n_func_def /) )

  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine Init_Functions

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Functions

! $Log$
! Revision 2.14  2013/10/09 01:04:35  vsnyder
! Add Difference, Intersection, Union, Without
!
! Revision 2.13  2013/09/19 23:31:50  vsnyder
! Allow number or dot for exp, logs, sqrt
!
! Revision 2.12  2013/09/04 00:01:39  pwagner
! Comments preceding use TREE made more uniform across modules
!
! Revision 2.11  2012/03/12 18:35:34  vsnyder
! Add ln and log10
!
! Revision 2.10  2011/04/26 20:32:53  vsnyder
! Repair definition of 'mod'
!
! Revision 2.9  2011/04/18 19:27:22  vsnyder
! Add MOD
!
! Revision 2.8  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.6  2004/10/14 06:06:20  livesey
! Added f_clearLower
!
! Revision 2.5  2004/05/29 02:42:36  vsnyder
! Rearrange function definition stuff
!
! Revision 2.4  2004/05/28 00:56:54  vsnyder
! Add log and exp
!
! Revision 2.3  2004/04/29 01:26:09  livesey
! Added xtx
!
! Revision 2.2  2004/01/30 23:25:31  livesey
! Added GetDiagonal and sqrt
!
! Revision 2.1  2004/01/17 03:05:06  vsnyder
! Initial commit
!
