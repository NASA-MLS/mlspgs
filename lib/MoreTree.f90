! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MoreTree

! Some routines for tree analysis that don't quite fit anywhere else.

  use Tree, only: Decoration, Node_ID, Subtree, nsons
  
  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  ! ------------------------------------------------  Get_Boolean  -----
  logical function Get_Boolean ( Root )
  ! Get the value of a field that is required to have type t_boolean.
  ! Return true if node_id(root) is n_set_one.  Otherwise return decoration of
  ! root, or of child of root if children
    use Intrinsic, only: L_true
    use Tree_Types, only: N_Set_one
    integer, intent(in) :: Root
    if ( node_id(root) == n_set_one ) then
      get_boolean = .true.
    else
      if ( nsons ( root ) /= 0 ) then
        get_boolean = decoration(subtree(2,root)) == l_true
      else
        get_boolean = decoration(root) == l_true
      end if
    end if
  end function Get_Boolean

  ! ---------------------------------------------  Get_Dot_Decors  -----
  subroutine Get_Dot_Decors ( Root, FirstDecor, SecondDecor )
  ! Assume that node_id(root) is n_dot.
  ! Set FirstDecor = the decoration of the declaration of the first son,
  ! and SecondDecor = the decoration of the declaration of the second son.
    integer, intent(in) :: Root
    integer, intent(out) :: FirstDecor, SecondDecor
    firstDecor = decoration(decoration(subtree(1,root)))
    secondDecor = decoration(decoration(decoration(subtree(2,root))))
  end subroutine Get_Dot_Decors

  ! -----------------------------------------------  Get_Field_Id  -----
  integer function Get_Field_Id ( Root )
  ! Assume that node_id(root) is either n_asg or n_set_one.
  ! Return the field ID.
    integer, intent(in) :: Root
    get_field_id = decoration(subtree(1,root))
  end function Get_Field_Id

  ! ------------------------------------------------  Get_Spec_Id  -----
  integer function Get_Spec_Id ( Root )
  ! Assume that node_id(root) is n_spec_args.
  ! Return the spec ID.
    integer, intent(in) :: Root
    get_spec_id = decoration(subtree(1,decoration(subtree(1,root))))
  end function Get_Spec_Id

  ! -----------------------------------  GetStringIndexFromString  -----
  integer function GetStringIndexFromString ( line, caseSensitive )
    use Symbol_Types, only: T_IDENTIFIER
    use Symbol_Table, only: ENTER_TERMINAL

    character (len=*), intent(in) :: LINE
    logical, optional, intent(in) :: CASESENSITIVE
    ! Executable code
    GetStringIndexFromString = enter_terminal ( trim(line), t_identifier, &
      & caseSensitive=caseSensitive )
  end function GetStringIndexFromString

  ! --------------------------------------  GetLitIndexFromString  -----
  integer function GetLitIndexFromString ( line, stringIndex )
    use Declaration_Table, only: GET_DECL, DECLS, ENUM_VALUE

    ! Return the lit index for a string if it's a lit, else -huge(0).
    character (len=*), intent(in) :: LINE
    integer, optional, intent(out) :: STRINGINDEX
    ! Local variable
    type (Decls) :: DECL
    integer :: SI
    ! Executable code
    si = getStringIndexFromString(line)
    if ( present ( stringIndex ) ) stringIndex = si
    decl = get_decl ( si, type=enum_value )
    getLitIndexFromString = decl%units
    if ( decl%type /= enum_value )  getLitIndexFromString = -huge(0)
  end function GetLitIndexFromString

  ! ------------------------------------------  StartErrorMessage  -----
  subroutine StartErrorMessage ( where )
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: SOURCE_REF
    integer, intent(in) :: Where             ! Tree node index
    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no tree available)' )
    end if
    call output ( ': ' )
  end subroutine StartErrorMessage

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MoreTree

! $Log$
! Revision 2.11  2004/06/16 19:51:25  vsnyder
! Move a use from module scope to procedure scope
!
! Revision 2.10  2004/06/16 01:45:13  vsnyder
! Return -huge(0) from GetLitIndexFromString if the string isn't a lit
!
! Revision 2.9  2004/05/28 00:57:25  vsnyder
! Move GetIndexFlagsFromList from MoreTree to Expr_m
!
! Revision 2.8  2004/05/21 22:52:33  vsnyder
! Add StartErrorMessage routine
!
! Revision 2.7  2004/01/22 00:41:40  pwagner
! GetStringIndexFromString takes optional arg caseSensitive
!
! Revision 2.6  2003/08/16 00:32:13  vsnyder
! Push uses down to procedure scope, futzing
!
! Revision 2.5  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/10/05 00:41:50  livesey
! Added GetStringIndexFromString and GetLitIndexFromString
!
! Revision 2.3  2002/08/26 20:01:22  livesey
! Added GetIndexFlagsFromList
!
! Revision 2.2  2002/06/07 17:51:56  livesey
! More versitility in get_boolean
!
! Revision 2.1  2001/02/23 01:16:36  vsnyder
! Getting L_True from Intrinsic instead of Init_Tables_Module: move to lib
!
! Revision 1.1  2001/02/20 22:49:25  vsnyder
! Initial commit
!
