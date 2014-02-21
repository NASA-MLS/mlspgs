! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MoreTree

! Some routines for tree analysis that don't quite fit anywhere else.

  use StartErrorMessage_m, only: StartErrorMessage ! Others use it from here too
  use Tree, only: Decoration, Node_ID, NSons, Subtree

  implicit NONE
  public

  interface FillArray
    module procedure FillDecorArray_I, FillDecorArray_TX
    module procedure FillStringArray_I, FillStringArray_TX
  end interface FillArray

  interface FillDecorArray
    module procedure FillDecorArray_I, FillDecorArray_TX
  end interface

  interface FillStringArray
    module procedure FillStringArray_I, FillStringArray_TX
  end interface

  interface FillSubrosaArray
    module procedure FillSubrosaArray_I, FillSubrosaArray_TX
  end interface

  interface Get_Boolean
    module procedure Get_Boolean_I, Get_Boolean_TX
  end interface

  interface Get_Dot_Decors
    module procedure Get_Dot_Decors_I, Get_Dot_Decors_TX, Get_Dot_Decors_TX_TX
  end interface

  interface Get_Field_Id
    module procedure Get_Field_Id_I, Get_Field_Id_TX
  end interface

  interface Get_Label_And_Spec
    module procedure Get_Label_And_Spec_I, Get_Label_And_Spec_TX
  end interface

  interface Get_Spec_Id
    module procedure Get_Spec_Id_I, Get_Spec_Id_TX
  end interface

  interface Scalar
    module procedure Scalar_I, Scalar_TX
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  ! -------------------------------------------  FillDecorArray_I  -----
  integer function FillDecorArray_I ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the decorations of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Allocate_Deallocate, only: Allocate_Test
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Tree, only: Node_Kind, Pseudo, Sub_Rosa

    integer, intent(in) :: Where ! in the tree
    integer, pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    integer :: Gson, J, N

    error = 0
    n = nsons(where)
    call allocate_test ( array, n-1, arrayName, moduleName )
    do j = 2, n
      gson = subtree(j,where)
      if ( node_kind(gson) == pseudo ) then
        array(j-1) = decoration(gson)
      else
        error = 1
        call startErrorMessage ( where )
        call output ( 'Range not allowed for ' )
        call display_string ( sub_rosa(subtree(1,where)), advance='yes' )
      end if
    end do
  end function FillDecorArray_I

  ! ------------------------------------------  FillDecorArray_TX  -----
  integer function FillDecorArray_TX ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the decorations of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Tree, only: TX
    type(tx), intent(in) :: Where ! in the tree
    integer, pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    error = fillDecorArray ( where%i, array, arrayName )
  end function FillDecorArray_TX

  ! ------------------------------------------  FillStringArray_I  -----
  integer function FillStringArray_I ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the sub-rosas of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Allocate_Deallocate, only: Allocate_Test
    use MLSStrings, only: Capitalize
    use Output_m, only: Output
    use String_Table, only: Display_String, Get_String
    use Tree, only: Node_Kind, Pseudo, Sub_Rosa

    integer, intent(in) :: Where ! in the tree
    character(len=*), pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    integer :: Gson, J, N

    error = 0
    n = nsons(where)
    call allocate_test ( array, n-1, arrayName, moduleName )
    do j = 2, n
      gson = subtree(j,where)
      if ( node_kind(gson) == pseudo ) then
        call get_string ( sub_rosa(gson), array(j-1), strip=.true. )
        array(j-1) = capitalize(array(j-1))
      else
        error = 1
        call startErrorMessage ( where )
        call output ( 'Range not allowed for ' )
        call display_string ( sub_rosa(subtree(1,where)), advance='yes' )
      end if
    end do
  end function FillStringArray_I

  ! -----------------------------------------  FillStringArray_TX  -----
  integer function FillStringArray_TX ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the sub-rosas of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Tree, only: TX

    type(tx), intent(in) :: Where ! in the tree
    character(len=*), pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    error = fillStringArray ( where%i, array, arrayName )
  end function FillStringArray_TX

  ! -----------------------------------------  FillSubrosaArray_I  -----
  integer function FillSubrosaArray_I ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the decorations of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Allocate_Deallocate, only: Allocate_Test
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Tree, only: Node_Kind, Pseudo, Sub_Rosa

    integer, intent(in) :: Where ! in the tree
    integer, pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    integer :: Gson, J, N

    error = 0
    n = nsons(where)
    call allocate_test ( array, n-1, arrayName, moduleName )
    do j = 2, n
      gson = subtree(j,where)
      if ( node_kind(gson) == pseudo ) then
        array(j-1) = sub_rosa(gson)
      else
        error = 1
        call startErrorMessage ( where )
        call output ( 'Range not allowed for ' )
        call display_string ( sub_rosa(subtree(1,where)), advance='yes' )
      end if
    end do
  end function FillSubrosaArray_I

  ! ----------------------------------------  FillSubrosaArray_TX  -----
  integer function FillSubrosaArray_TX ( Where, Array, ArrayName ) result(Error)
    ! Fill Array with the decorations of sons 2..n of Where.
    ! Result > 0 => field had a range in it.
    ! Array is allocated here with Allocate_Test, so don't send me an
    ! undefined pointer!
    use Tree, only: TX

    type(tx), intent(in) :: Where ! in the tree
    integer, pointer :: Array(:)
    character(len=*), intent(in) :: ArrayName ! For error message
    error = fillSubrosaArray ( where%i, array, arrayName )
  end function FillSubrosaArray_TX 

  ! ----------------------------------------------  Get_Boolean_I  -----
  logical function Get_Boolean_I ( Root ) result ( Boolean )
  ! Get the value of a field that is required to have type t_boolean.
  ! Return true if node_id(root) is n_set_one.  Otherwise return decoration of
  ! root, or of child of root if children
    use Expr_m, only: Expr
    use Intrinsic, only: L_true
    use Tree_Types, only: N_identifier, N_Set_one
    integer, intent(in) :: Root
    integer :: Son, Units(2)
    double precision :: Value(2)
    if ( node_id(root) == n_set_one ) then
      boolean = .true.
    else
      if ( nsons ( root ) /= 0 ) then
        son = subtree(2,root)
        if ( node_id(son) == n_identifier ) then
          boolean = decoration(son) == l_true
        else
          call expr ( son, units, value )
          boolean = nint(value(1)) == l_true
        end if
      else
        boolean = decoration(root) == l_true
      end if
    end if
  end function Get_Boolean_I

  ! --------------------------------------------  Get_Boolean_TX  -----
  logical function Get_Boolean_TX ( Root ) result ( Boolean )
  ! Get the value of a field that is required to have type t_boolean.
  ! Return true if node_id(root) is n_set_one.  Otherwise return decoration of
  ! root, or of child of root if children
    use Tree, only: TX
    type(tx), intent(in) :: Root
    boolean = get_boolean ( root%i )
  end function Get_Boolean_TX

  ! -------------------------------------------  Get_Dot_Decors_I  -----
  subroutine Get_Dot_Decors_I ( Root, FirstDecor, SecondDecor )
  ! Assume that node_id(root) is n_dot.
  ! Set FirstDecor = the decoration of the declaration of the first son,
  ! and SecondDecor = the decoration of the declaration of the second son.
    integer, intent(in) :: Root
    integer, intent(out) :: FirstDecor, SecondDecor
    firstDecor = decoration(decoration(subtree(1,root)))
    secondDecor = decoration(decoration(decoration(subtree(2,root))))
  end subroutine Get_Dot_Decors_I

  ! ------------------------------------------  Get_Dot_Decors_TX  -----
  subroutine Get_Dot_Decors_TX ( Root, FirstDecor, SecondDecor )
  ! Assume that node_id(root) is n_dot.
  ! Set FirstDecor = the decoration of the declaration of the first son,
  ! and SecondDecor = the decoration of the declaration of the second son.
    use Tree, only: TX
    type(tx), intent(in) :: Root
    integer, intent(out) :: FirstDecor, SecondDecor
    firstDecor = decoration(decoration(subtree(1,root)))
    secondDecor = decoration(decoration(decoration(subtree(2,root))))
  end subroutine Get_Dot_Decors_TX

  ! ---------------------------------------  Get_Dot_Decors_TX_TX  -----
  subroutine Get_Dot_Decors_TX_TX ( Root, FirstDecor, SecondDecor )
  ! Assume that node_id(root) is n_dot.
  ! Set FirstDecor = the decoration of the declaration of the first son,
  ! and SecondDecor = the decoration of the declaration of the second son.
    use Tree, only: TX, Decoration_TX_TX
    type(tx), intent(in) :: Root
    type(tx), intent(out) :: FirstDecor, SecondDecor
    firstDecor = decoration_tx_tx(decoration_tx_tx(subtree(1,root)))
    secondDecor = decoration_tx_tx(decoration_tx_tx(decoration_tx_tx(subtree(2,root))))
  end subroutine Get_Dot_Decors_TX_TX

  ! ---------------------------------------------  Get_Field_Id_I  -----
  integer function Get_Field_Id_I ( Root ) result ( Get_Field_Id )
  ! Assume that node_id(root) is either n_asg or n_set_one.
  ! Return the field ID.
    integer, intent(in) :: Root
    get_field_id = decoration(subtree(1,root))
  end function Get_Field_Id_I

  ! ---------------------------------------------  Get_Field_Id_TX  -----
  integer function Get_Field_Id_TX ( Root ) result ( Get_Field_Id )
  ! Assume that node_id(root) is either n_asg or n_set_one.
  ! Return the field ID.
    use Tree, only: TX
    type(tx), intent(in) :: Root
    get_field_id = decoration(subtree(1,root))
  end function Get_Field_Id_TX

  ! ---------------------------------------  Get_Label_And_Spec_I  -----
  subroutine Get_Label_And_Spec_I ( Root, Label, Spec )
  ! Starting at root, if its node_id is n_named, get Label from sub_rosa
  ! of its first son and Spec from its second.  Otherwise, set Label = 0
  ! and Spec = Root.  If Spec is absent, set Root to what Spec would be.
    use Tree, only: Sub_Rosa
    use Tree_Types, only: N_Named
    integer, intent(inout) :: Root ! inout in case Spec is absent
    integer, intent(out) :: Label
    integer, intent(out), optional :: Spec
    integer :: MySpec
    if ( node_id(root) == n_named ) then
      label = sub_rosa(subtree(1,root))
      mySpec = subtree(2,root)
    else
      label = 0
      mySpec = root
    end if
    if ( present(spec) ) then
      spec = mySpec
    else
      root = mySpec
    end if
  end subroutine Get_Label_And_Spec_I

  ! --------------------------------------  Get_Label_And_Spec_TX  -----
  subroutine Get_Label_And_Spec_TX ( Root, Label, Spec )
  ! Starting at root, if its node_id is n_named, get Label from sub_rosa
  ! of its first son and Spec from its second.  Otherwise, set Label = 0
  ! and Spec = Root.
    use Tree, only: Sub_Rosa, TX
    use Tree_Types, only: N_Named
    type(tx), intent(inout) :: Root ! inout in case Spec is absent
    integer, intent(out) :: Label
    type(tx), intent(out), optional :: Spec
    type(tx) :: MySpec
    if ( node_id(root) == n_named ) then
      label = sub_rosa(subtree(1,root))
      mySpec = subtree(2,root)
    else
      label = 0
      mySpec = root
    end if
    if ( present(spec) ) then
      spec = mySpec
    else
      root = mySpec
    end if
  end subroutine Get_Label_And_Spec_TX

  ! ----------------------------------------------  Get_Spec_Id_I  -----
  integer function Get_Spec_Id_I ( Root ) result ( Get_Spec_Id )
  ! Assume that node_id(root) is n_spec_args.
  ! Return the spec ID.
    integer, intent(in) :: Root
    get_spec_id = decoration(subtree(1,decoration(subtree(1,root))))
  end function Get_Spec_Id_I

  ! ----------------------------------------------  Get_Spec_Id_TX  -----
  integer function Get_Spec_Id_TX ( Root ) result ( Get_Spec_Id )
  ! Assume that node_id(root) is n_spec_args.
  ! Return the spec ID.
    use Tree, only: TX, Decoration_TX_TX
    type(tx), intent(in) :: Root
    get_spec_id = decoration(subtree(1,decoration_tx_tx(subtree(1,root))))
  end function Get_Spec_Id_TX

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

  ! ---------------------------------------------------  Scalar_I  -----
  logical function Scalar_I ( Root ) result ( IsScalar )
    use Tree, only: Node_Id, Nsons, SUbtree
    use Tree_Types, only: N_Array

    ! Return "root has two sons, and the second one is not N_Array"
    integer, intent(in) :: Root

    isScalar = .false.
    if ( nsons(root) > 2 ) return
    if ( node_id(subtree(2,root)) == n_array ) return
    isScalar = .true.
  end function Scalar_I

  ! --------------------------------------------------  Scalar_TX  -----
  logical function Scalar_TX ( Root ) result ( IsScalar )
    use Tree, only: TX
    ! Return "root has two sons, and the second one is not N_Array"
    type(tx), intent(in) :: Root
    isScalar = scalar(root%i)
  end function Scalar_TX

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MoreTree

! $Log$
! Revision 2.23  2014/02/21 19:25:24  vsnyder
! Result of Get_Boolean depends on result of expr being l_true, not nonzero
!
! Revision 2.22  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.21  2013/12/12 01:56:44  vsnyder
! Add Get_Label_And_Spec
!
! Revision 2.20  2013/10/02 02:07:44  vsnyder
! Cannonball polishing
!
! Revision 2.19  2013/09/30 23:59:57  vsnyder
! Routines for TX type, move StartErrorMessage from include to module
!
! Revision 2.18  2012/05/07 23:00:57  vsnyder
! StartErrorMessage moved to include to avoid a circular dependence
! between expr_m and MoreTree
!
! Revision 2.17  2012/05/05 00:12:15  vsnyder
! Process logical expressions in get_boolean
!
! Revision 2.16  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.15  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2005/05/02 22:57:29  vsnyder
! Add FillSubRosaArray
!
! Revision 2.13  2005/01/20 01:31:07  vsnyder
! Add FillArray, FillIntegerArray, FillStringArray
!
! Revision 2.12  2004/11/17 20:24:20  vsnyder
! Add SCALAR function to check for scalar field value
!
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
