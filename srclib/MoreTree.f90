! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MoreTree

! Some routines for tree analysis that don't quite fit anywhere else.

  use Init_Tables_Module, only: L_true
  use Tree, only: Decoration, Node_ID, Subtree
  use Tree_Types, only: N_Set_one

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  ! ------------------------------------------------  Get_Boolean  -----
  logical function Get_Boolean ( Root )
  ! Get the value of a field that is required to have type t_boolean.
  ! Return true if node_id(root) is n_set_one.  Otherwise return
  ! decoration(subtree(2,root)) == l_true.
    integer, intent(in) :: Root
    if ( node_id(root) == n_set_one ) then
      get_boolean = .true.
    else
      get_boolean = decoration(subtree(2,root)) == l_true
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

end module MoreTree

! $Log$
