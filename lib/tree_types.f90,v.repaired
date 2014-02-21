! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TREE_TYPES
  use STRING_TABLE, only: ADD_CHAR
  use SYMBOL_TYPES, only: MIN_PSEUDO, MAX_PSEUDO
  public

  ! Parameters giving tree node indices.  The first two are required.
  integer, parameter :: FIRST_TREE_NODE = 1
  integer, parameter :: N_NULL = FIRST_TREE_NODE ! Used in InitTree and Parser.
    ! Pseudo-terminals must be contiguous
  integer, parameter :: N_IDENTIFIER =  n_null + 1
  integer, parameter :: N_NUMBER =      n_identifier + 1
  integer, parameter :: N_STRING =      n_number + 1

  integer, parameter :: N_AND =         n_string + 1
  integer, parameter :: N_ARG_DEF =     n_and + 1        ! Types of func args
  integer, parameter :: N_ARRAY =       n_arg_def + 1    ! []
  integer, parameter :: N_ASG =         n_array + 1      ! = in spec
  integer, parameter :: N_CASE =        n_asg +1         ! CASE in SELECT block
  integer, parameter :: N_CF =          n_case +1        ! A spec (config)
  integer, parameter :: N_CFS =         n_cf +1          ! All the configs
  integer, parameter :: N_COLON =       n_cfs + 1        ! A range
  integer, parameter :: N_COLON_LESS =  n_colon + 1      ! Range open on right
  integer, parameter :: N_COND =        n_colon_less + 1
  integer, parameter :: N_CYCLE =       n_cond + 1
  integer, parameter :: N_DEFAULT =     n_cycle + 1
  integer, parameter :: N_DIV =         n_default + 1
  integer, parameter :: N_DO =          n_div + 1       
  integer, parameter :: N_DO_HEAD =     n_do + 1        
  integer, parameter :: N_DOT =         n_do_head + 1
  integer, parameter :: N_DT_DEF =      n_dot + 1        ! Data type definition
  integer, parameter :: N_ELSE =        n_dt_def + 1     ! ELSE block in IF construct
  integer, parameter :: N_EQUAL =       n_else + 1       ! = in X = Y spec
  integer, parameter :: N_EQUAL_EQUAL = n_equal + 1
  integer, parameter :: N_EXIT =        n_equal_equal + 1
  integer, parameter :: N_FIELD_SPEC =  n_exit + 1
  integer, parameter :: N_FIELD_TYPE =  n_field_spec + 1
  integer, parameter :: N_FUNC_DEF =    n_field_type + 1
  integer, parameter :: N_FUNC_REF =    n_func_def + 1
  integer, parameter :: N_GREATER =     n_func_ref + 1
  integer, parameter :: N_GREATER_EQ =  n_greater + 1
  integer, parameter :: N_IF =          n_greater_eq + 1 ! IF construct
  integer, parameter :: N_INTO =        n_if + 1         ! A divided into B: A\B
  integer, parameter :: N_LESS =        n_into + 1
  integer, parameter :: N_LESS_COLON =  n_less + 1       ! Range open on left
  integer, parameter :: N_LESS_COLON_LESS = n_less_colon + 1 ! Open range
  integer, parameter :: N_LESS_EQ =     n_less_colon_less + 1
  integer, parameter :: N_MINUS =       n_less_eq + 1
  integer, parameter :: N_MULT =        n_minus + 1
  integer, parameter :: N_NAME_DEF =    n_mult + 1       ! Names OK in section
  integer, parameter :: N_NAMED =       n_name_def + 1   ! labelled spec
  integer, parameter :: N_NOT =         n_named + 1
  integer, parameter :: N_NOT_EQUAL =   n_not + 1
  integer, parameter :: N_OR =          n_not_equal + 1
  integer, parameter :: N_PLUS =        n_or + 1
  integer, parameter :: N_POW =         n_plus + 1
  integer, parameter :: N_SECTION =     n_pow + 1        ! what specs in section?
  integer, parameter :: N_SELECT =      n_section + 1    ! SELECT CASE construct
  integer, parameter :: N_SET_ONE =     n_select + 1     ! /foo
  integer, parameter :: N_SPEC_ARGS =   n_set_one + 1    ! x,y=z,...
  integer, parameter :: N_SPEC_DEF =    n_spec_args + 1  ! what fields in spec?
  integer, parameter :: N_SUBSCRIPT =   n_spec_def + 1   ! name [ expr ]
  integer, parameter :: N_TEST =        n_subscript + 1  ! expr in IF or ELSE IF
  integer, parameter :: N_UNCHECKED =   n_test + 1       ! field type not checked
  integer, parameter :: N_UNIT =        n_unchecked + 1  ! number // units
  integer, parameter :: N_VARIABLE =    n_unit + 1       ! := in name := expr
  integer, parameter :: N_VARIABLE_REF = n_variable + 1  ! field needs a variable
  integer, parameter :: N_WHILE =       n_variable_ref + 1

  integer, parameter :: LAST_TREE_NODE = N_WHILE

  ! mapping from pseudo-terminal indices to corresponding tree nodes.
  integer, parameter :: tree_map ( min_pseudo: max_pseudo ) = &
  !    t_identifier, t_number, t_string
    (/ n_identifier, n_number, n_string /)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine TREE_INIT ( TREE_NODE )
  ! Put the name of a tree node into the character table.  Everything
  ! else is handled in TREE % INIT_TREE

    integer, intent(in) :: TREE_NODE    ! One of the parameters above

    select case ( tree_node )
    case ( n_Null );            call add_char ( '<null>' )
    case ( n_Identifier );      call add_char ( '<identifier>' )
    case ( n_Number );          call add_char ( '<number>' )
    case ( n_String );          call add_char ( '<string>' )
    case ( n_And );             call add_char ( '<and>' )
    case ( n_Arg_Def );         call add_char ( '<arg_def>' )
    case ( n_Array );           call add_char ( '<array>' )
    case ( n_Asg );             call add_char ( '<assign>' )
    case ( n_Case ) ;           call add_char ( '<case>' )
    case ( n_Cf );              call add_char ( '<cf>' )
    case ( n_Cfs );             call add_char ( '<cfs>' )
    case ( n_Colon );           call add_char ( '<colon>' )
    case ( n_Colon_less );      call add_char ( '<colon_less>' )
    case ( n_Cond );            call add_char ( '<cond>' )
    case ( n_Cycle );           call add_char ( '<cycle>' )
    case ( n_Default );         call add_char ( '<default>' )
    case ( n_Div );             call add_char ( '<div>' )
    case ( n_Do );              call add_char ( '<do>' )
    case ( n_Do_Head );         call add_char ( '<do_head>' )
    case ( n_Dot );             call add_char ( '<dot>' )
    case ( n_DT_Def );          call add_char ( '<dt_def>' )
    case ( n_Else );            call add_char ( '<else>' )
    case ( n_Equal );           call add_char ( '<equal>' )
    case ( n_Equal_Equal );     call add_char ( '<equal_equal>' )
    case ( n_Exit );            call add_char ( '<exit>' )
    case ( n_Field_Spec );      call add_char ( '<field_spec>' )
    case ( n_Field_Type );      call add_char ( '<field_type>' )
    case ( n_Func_Def );        call add_char ( '<func_def>' )
    case ( n_Func_Ref );        call add_char ( '<func_ref>' )
    case ( n_Greater );         call add_char ( '<greater>' )
    case ( n_Greater_eq );      call add_char ( '<greater_eq>' )
    case ( n_If );              call add_char ( '<if>' )
    case ( n_Into );            call add_char ( '<into>' )
    case ( n_Less );            call add_char ( '<less>' )
    case ( n_Less_colon );      call add_char ( '<less_colon>' )
    case ( n_Less_colon_less ); call add_char ( '<less_colon_less>' )
    case ( n_Less_eq );         call add_char ( '<less_eq>' )
    case ( n_Minus );           call add_char ( '<minus>' )
    case ( n_Mult );            call add_char ( '<mult>' )
    case ( n_Name_def );        call add_char ( '<name_def>' )
    case ( n_Named );           call add_char ( '<named>' )
    case ( n_Not );             call add_char ( '<not>' )
    case ( n_Not_Equal );       call add_char ( '<not_equal>' )
    case ( n_Or );              call add_char ( '<or>' )
    case ( n_Plus );            call add_char ( '<plus>' )
    case ( n_Pow );             call add_char ( '<power>' )
    case ( n_Section );         call add_char ( '<section>' )
    case ( n_Select );          call add_char ( '<select>' )
    case ( n_Set_one );         call add_char ( '<set_one>' )
    case ( n_Spec_args );       call add_char ( '<spec_args>' )
    case ( n_Spec_def );        call add_char ( '<spec_def>' )
    case ( n_Subscript );       call add_char ( '<subscript>' )
    case ( n_Test );            call add_char ( '<test>' )
    case ( n_Unchecked );       call add_char ( '<unchecked>' )
    case ( n_Unit );            call add_char ( '<unit>' )
    case ( n_Variable );        call add_char ( '<variable>' )
    case ( n_Variable_ref );    call add_char ( '<variable_ref>' )
    case ( n_While );           call add_char ( '<while>' )
    case default
      write ( *,* ) 'TREE_TYPES%TREE_INIT-E- No initializer for &
                    &tree node with index ', tree_node
      stop
    end select

  end subroutine TREE_INIT

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module TREE_TYPES

! $Log$
! Revision 2.19  2014/02/21 19:21:00  vsnyder
! Add CYCLE, DO, EXIT, WHILE
!
! Revision 2.18  2013/12/12 01:57:49  vsnyder
! Add variable definition, and IF and SELECT constructs
!
! Revision 2.17  2013/11/26 22:44:08  vsnyder
! Add CASE, DEFAULT, ELSE, IF, SELECT, TEST and SUBSCRIPT nodes
!
! Revision 2.16  2013/10/09 01:10:30  vsnyder
! Add Variable_Ref tree node for field type definition
!
! Revision 2.15  2013/10/02 01:34:46  vsnyder
! Add conditional ?...! and variable assignment := tree nodes
!
! Revision 2.14  2012/05/05 00:11:51  vsnyder
! Add support for 'not' operator
!
! Revision 2.13  2012/05/01 22:10:26  vsnyder
! Add TrueList subroutine
!
! Revision 2.12  2011/04/19 01:59:43  vsnyder
! Support == and /= relational operators too
!
! Revision 2.11  2011/04/18 19:33:26  vsnyder
! Add support for relational operators and boolean-valued expressions
!
! Revision 2.10  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/05/29 02:42:59  vsnyder
! Rearrange function definition stuff
!
! Revision 2.7  2004/05/28 23:12:21  vsnyder
! Add power (^) operator
!
! Revision 2.6  2004/01/17 03:04:48  vsnyder
! Provide for functions in expressions
!
! Revision 2.5  2004/01/14 18:32:58  vsnyder
! Stuff for Algebra module
!
! Revision 2.4  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.3  2001/11/28 03:15:19  vsnyder
! Implement arrays of arrays
!
! Revision 2.2  2001/11/27 00:54:37  vsnyder
! Implement (partially) open ranges
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
