! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE_TYPES
  use STRING_TABLE, only: ADD_CHAR
  use SYMBOL_TYPES, only: MIN_PSEUDO, MAX_PSEUDO
  public

  ! Parameters giving tree node indices.  The first two are required.
  integer, parameter :: N_EOF = 0       ! Root of the parse tree.
  integer, parameter :: N_NULL = 1      ! Used in InitTree.
    ! Pseudo-terminals must be contiguous
  integer, parameter :: N_IDENTIFIER = n_null + 1
  integer, parameter :: N_NUMBER =     n_identifier + 1
  integer, parameter :: N_STRING =     n_number + 1

  integer, parameter :: N_AND =        n_string + 1
  integer, parameter :: N_ARRAY =      n_and + 1       ! []
  integer, parameter :: N_ASG =        n_array + 1     ! = in spec
  integer, parameter :: N_CF =         n_asg +1        ! A spec (config)
  integer, parameter :: N_CFS =        n_cf +1         ! All the configs
  integer, parameter :: N_COLON =      n_cfs + 1       ! A range
  integer, parameter :: N_COLON_LESS = n_colon + 1     ! Range open on right
  integer, parameter :: N_DIV =        n_colon_less + 1
  integer, parameter :: N_DOT =        n_div + 1
  integer, parameter :: N_DT_DEF =     n_dot + 1       ! Data type definition
  integer, parameter :: N_EQUAL =      n_dt_def + 1    ! = in X = Y spec
  integer, parameter :: N_FIELD_SPEC = n_equal + 1
  integer, parameter :: N_FIELD_TYPE = n_field_spec + 1
  integer, parameter :: N_LESS_COLON = n_field_type + 1 ! Range open on left
  integer, parameter :: N_LESS_COLON_LESS = n_less_colon + 1 ! Open range
  integer, parameter :: N_MINUS =      n_less_colon_less + 1
  integer, parameter :: N_MULT =       n_minus + 1
  integer, parameter :: N_NAME_DEF =   n_mult + 1      ! Names OK in section
  integer, parameter :: N_NAMED =      n_name_def + 1  ! labelled spec
  integer, parameter :: N_OR =         n_named + 1
  integer, parameter :: N_PLUS =       n_or + 1
  integer, parameter :: N_SECTION =    n_plus + 1      ! what specs in section?
  integer, parameter :: N_SET_ONE =    n_section + 1   ! /foo
  integer, parameter :: N_SPEC_ARGS =  n_set_one + 1   ! x,y=z,...
  integer, parameter :: N_SPEC_DEF =   n_spec_args + 1 ! what fields in spec?
  integer, parameter :: N_UNIT =       n_spec_def + 1  ! number // units

  integer, parameter :: LAST_TREE_NODE = N_UNIT

  ! mapping from pseudo-terminal indices to corresponding tree nodes.
  integer, parameter :: tree_map ( min_pseudo: max_pseudo ) = &
  !    t_identifier, t_number
    (/ n_identifier, n_number, n_string /)

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine TREE_INIT ( TREE_NODE )
  ! Put the name of a tree node into the character table.  Everything
  ! else is handled in TREE % INIT_TREE
    integer, intent(in) :: TREE_NODE    ! One of the parameters above
    select case ( tree_node )
    case ( n_Eof );        call add_char ( 'eof' )
    case ( n_Null );       call add_char ( 'null' )
    case ( n_Identifier ); call add_char ( 'identifier' )
    case ( n_Number );     call add_char ( 'number' )
    case ( n_String );     call add_char ( 'string' )
    case ( n_And );        call add_char ( 'and' )
    case ( n_Array );      call add_char ( 'array' )
    case ( n_Asg );        call add_char ( 'assign' )
    case ( n_Cf );         call add_char ( 'cf' )
    case ( n_Cfs );        call add_char ( 'cfs' )
    case ( n_Colon );      call add_char ( 'colon' )
    case ( n_Colon_less ); call add_char ( 'colon_less' )
    case ( n_Div );        call add_char ( 'div' )
    case ( n_Dot );        call add_char ( 'dot' )
    case ( n_DT_Def );     call add_char ( 'dt_def' )
    case ( n_Equal );      call add_char ( 'equal' )
    case ( n_Field_Spec ); call add_char ( 'field_spec' )
    case ( n_Field_Type ); call add_char ( 'field_type' )
    case ( n_Less_colon ); call add_char ( 'less_colon' )
    case ( n_Less_colon_less ); call add_char ( 'less_colon_less' )
    case ( n_Minus );      call add_char ( 'minus' )
    case ( n_Mult );       call add_char ( 'mult' )
    case ( n_Name_def );   call add_char ( 'name_def' )
    case ( n_Named );      call add_char ( 'named' )
    case ( n_Or );         call add_char ( 'or' )
    case ( n_Plus );       call add_char ( 'plus' )
    case ( n_Section );    call add_char ( 'section' )
    case ( n_Set_one );    call add_char ( 'set_one' )
    case ( n_Spec_args );  call add_char ( 'spec_args' )
    case ( n_Spec_def );   call add_char ( 'spec_def' )
    case ( n_Unit );       call add_char ( 'unit' )
    case default
      write ( *,* )'TREE_TYPES%TREE_INIT-E- No initializer for &
                     &tree node with index ', tree_node
      stop
    end select
  end subroutine TREE_INIT
end module TREE_TYPES

! $Log$
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
