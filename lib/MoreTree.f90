! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MoreTree

! Some routines for tree analysis that don't quite fit anywhere else.

  use Declaration_table, only: NUM_VALUE
  use Tree, only: Decoration, Node_ID, Subtree, nsons
  use Tree_Types, only: N_colon, N_colon_less, N_less_colon, &
    & N_less_colon_less, N_Set_one
  
  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
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

  ! ----------------------------------------------- GetIndexFlagsFromList

  subroutine GetIndexFlagsFromList ( root, flags, status, lower, noError )
    ! Given the root of a numeric/numeric range array
    ! Set the flags array appropriately
    use Expr_m, only: EXPR
    use Intrinsic, only: PHYQ_DIMENSIONLESS
    use MLSCommon, only: R8
    integer, intent(in) :: ROOT         ! Tree node
    logical, dimension(:), intent(inout) :: FLAGS ! Result
    integer, intent(out) :: STATUS      ! Error flag, 0=success
    integer, intent(in), optional :: LOWER ! Lower bound for result
    logical, intent(in), optional :: NOERROR ! If set don't give bounds errors

    ! Local variables
    integer :: I,J                      ! Loop counters
    real(r8), dimension(2) :: VALUE     ! From expr
    integer, dimension(2) :: UNITS      ! From expr
    integer :: TYPE                     ! From expr
    integer :: RANGE_LOW, RANGE_HI      ! Range for flags
    logical :: MYNOERROR                ! Copy of noError
    integer :: MYLOWER                  ! Copy of lower
    integer :: SON                      ! Tree node

    ! Executable code
    flags = .false.
    status = 0
    myNoError = .false.
    if ( present ( noError ) ) myNoError = noError
    mylower = 1
    if ( present ( lower ) ) myLower = lower

    do j = 2, nsons(root)
      son = subtree ( j, root )
      call expr ( son, units, value, type )
      do i = 1, merge(1,2,type==num_value)
        if ( units(i) /= phyq_dimensionless ) then
          status = 1
          return
        end if
      end do
      range_low = nint(value(1))
      range_hi = nint(value(merge(1,2,type==num_value)))
      select case ( node_id(son) )
      case ( n_colon_less )
        range_hi = range_hi - 1
      case ( n_less_colon )
        range_low = range_low + 1
      case ( n_less_colon_less )
        range_low = range_low + 1
        range_hi = range_hi - 1
      end select
      if ( .not. myNoError .and. &
        & ( range_low < myLower .or. range_hi > myLower+size(flags)-1 ) ) then
        status = 1
        return
      end if
      range_low = max ( range_low, myLower )
      range_hi = min ( range_hi, myLower+size(flags)-1 )
      flags ( range_low-myLower+1 : range_hi-myLower+1 ) = .true.
    end do
  end subroutine GetIndexFlagsFromList

  ! ------------------------------------------ GetStringIndexFromString ---
  integer function GetStringIndexFromString ( line )
    use Symbol_Types, only: T_IDENTIFIER
    use Symbol_Table, only: ENTER_TERMINAL

    character (len=*), intent(in) :: LINE
    ! Executable code
    GetStringIndexFromString = enter_terminal ( trim(line), t_identifier )
  end function GetStringIndexFromString

  ! ------------------------------------------ GetLitIndexFromString ---
  integer function GetLitIndexFromString ( line, stringIndex )
    use Declaration_Table, only: GET_DECL, DECLS, ENUM_VALUE

    character (len=*), intent(in) :: LINE
    integer, optional, intent(out) :: STRINGINDEX
    ! Local variable
    type (Decls) :: DECL
    integer :: SI
    ! Executable code
    si = GetStringIndexFromString(line)
    if ( present ( stringIndex ) ) stringIndex = si
    decl = get_decl ( si, type=enum_value )
    GetLitIndexFromString = decl%units
  end function GetLitIndexFromString

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MoreTree

! $Log$
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
