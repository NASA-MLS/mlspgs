! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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
  integer, parameter :: F_Exp             = F_ClearLower + 1
  integer, parameter :: F_GetDiagonal     = F_Exp + 1
  integer, parameter :: F_Invert          = F_GetDiagonal + 1
  integer, parameter :: F_Log             = F_Invert + 1
  integer, parameter :: F_SQRT            = F_Log + 1
  integer, parameter :: F_Transpose       = F_SQRT + 1
  integer, parameter :: F_XTX             = F_Transpose + 1
  integer, parameter :: Function_Last     = F_XTX

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), parameter, private :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Init_Functions ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES )

    use Intrinsic, only: Begin, Add_Ident, Data_Type_Indices, &
      & F, Field_Indices, Func_Indices, &
      & G, Init_Intrinsic, L, Lit_Indices, &
      & N, P, Parm_Indices, S, Spec_Indices, Section_Indices, T, T_Numeric, Z
    use Tree_Types, only: N_Arg_Def, N_Func_Def

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES

    ! Initialize the intrinsic types

    call init_intrinsic ( n_data_type_indices, n_field_indices, n_lit_indices, &
      & first_parm_index, last_parm_index, n_section_indices, n_spec_indices, &
      & function_last )

    func_indices(f_cholesky) =        add_ident ( 'cholesky' )
    func_indices(f_clearLower) =      add_ident ( 'clearLower' )
    func_indices(f_exp) =             add_ident ( 'exp' )
    func_indices(f_getDiagonal) =     add_ident ( 'getDiagonal' )
    func_indices(f_invert) =          add_ident ( 'invert' )
    func_indices(f_log) =             add_ident ( 'log' )
    func_indices(f_sqrt) =            add_ident ( 'sqrt' )
    func_indices(f_transpose) =       add_ident ( 'transpose' )
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
      begin, g+f_exp, &
             begin, t+t_numeric, n+n_arg_def, n+n_func_def, &
      begin, g+f_getDiagonal, n+n_func_def, &
      begin, g+f_invert, n+n_func_def, &
      begin, g+f_log, &
             begin, t+t_numeric, n+n_arg_def, n+n_func_def, &
      !??? Automatic type checking for f_sqrt may be difficult, given that
      !??? we want to allow numbers of matrices.  It may be necessary either
      !??? to give up on that, or to do init_functions after init_tables.
      begin, g+f_sqrt, &
             begin, t+t_numeric, n+n_arg_def, n+n_func_def, &
      begin, g+f_transpose, n+n_func_def, &
      begin, g+f_xtx, n+n_func_def /) )

  contains

    ! ------------------------------------------------  MAKE_TREE  -----
    include "make_tree.f9h"

  end subroutine Init_Functions

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Functions

! $Log$
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
