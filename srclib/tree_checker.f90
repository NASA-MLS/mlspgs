! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE_CHECKER

! Traverse the tree output by the parser, which includes definitions
! put into it by init_tables_module before the parser runs.  Check
! that things fit together.  Fill in the declaration table.  Decorate
! cross-references within the tree with anything that might be useful.
! Decorate internal vertices in expressions with the index of the
! result type from Declaration_Table.

  use DECLARATION_TABLE, only: DECLARATION, DECLARE, DECLARED, DECLS, &
    &                          DO_LABEL, DUMP_1_DECL, EMPTY, ENUM_VALUE, &
    &                          FIELD, FUNCTION, GET_DECL, LABEL, LOG_VALUE, &
    &                          NAMED_VALUE, NULL_DECL, NUM_VALUE, &
    &                          PRIOR_DECL, RANGE, REDECLARE, SECTION, SPEC, &
    &                          STR_Range, STR_VALUE, TYPE_NAME, UNITS_NAME, &
    &                          VALUE_T, VARIABLE
  use INIT_TABLES_MODULE, only: DATA_TYPE_INDICES, FIELD_FIRST, FIELD_INDICES, &
    &                           FIELD_LAST, LIT_INDICES, PHYQ_DIMENSIONLESS, &
    &                           SECTION_FIRST, SECTION_INDICES, SECTION_LAST, &
    &                           SECTION_ORDERING
  use INTRINSIC, only: ALL_FIELDS, Array_Array, DOT => T_A_dot_B, EMPTY_OK, &
    &                  EXPR_OK, L_True, L_False, NO_ARRAY, NO_CHECK_EQ, NO_DUP, &
    &                  NO_POSITIONAL, PHYQ_INVALID, REQ_FLD, SPEC_INDICES, &
    &                  T_BOOLEAN, U
  use LEXER_CORE, only: PRINT_SOURCE
  use MoreTree, only: Scalar, StartErrorMessage
  use OUTPUT_M, only: NEWLINE, OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING, FLOAT_VALUE
  use TOGGLES, only: CON, Levels, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, &
                  NODE_ID, NODE_KIND, NSONS, NULL_TREE, PSEUDO, SUB_ROSA, &
                  SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  implicit NONE
  private

  public :: CHECK_TREE, CHECK_TYPE

! -----     Private declarations     -----------------------------------
  integer, private :: ERROR        ! 0 => No errors
  integer, private :: NUM_SECTIONS ! Number of begin ... ends

! Error codes for "announce_error"
  integer, private, parameter :: ALREADY_DECLARED = 1
  integer, private, parameter :: ARRAY_NOT_ALLOWED = ALREADY_DECLARED + 1
  integer, private, parameter :: EMPTY_NOT_ALLOWED = ARRAY_NOT_ALLOWED + 1
  integer, private, parameter :: IMPROPER_CYCLE = EMPTY_NOT_ALLOWED + 1
  integer, private, parameter :: IMPROPER_CYCLE_EXIT = IMPROPER_CYCLE + 1
  integer, private, parameter :: INCONSISTENT_DATA_TYPES = IMPROPER_CYCLE_EXIT + 1
  integer, private, parameter :: INCONSISTENT_TYPES = INCONSISTENT_DATA_TYPES + 1
  integer, private, parameter :: INCONSISTENT_UNITS = INCONSISTENT_TYPES + 1
  integer, private, parameter :: LABEL_CONFLICT = INCONSISTENT_UNITS + 1
  integer, private, parameter :: MISSING_FIELD = LABEL_CONFLICT + 1
  integer, private, parameter :: NO_ARRAY_ARRAY = MISSING_FIELD + 1
  integer, private, parameter :: NO_CODE_FOR = NO_ARRAY_ARRAY + 1
  integer, private, parameter :: NO_CYCLE_EXIT_TARGET = NO_CODE_FOR + 1
  integer, private, parameter :: NO_DECLARATION = NO_CYCLE_EXIT_TARGET + 1
  integer, private, parameter :: NO_DOT = NO_DECLARATION + 1
  integer, private, parameter :: NO_DUPLICATE_FIELDS = NO_DOT + 1
  integer, private, parameter :: NO_EXPR = NO_DUPLICATE_FIELDS + 1
  integer, private, parameter :: NO_POSITIONAL_FIELDS = NO_EXPR + 1
  integer, private, parameter :: NO_SUCH_FIELD = NO_POSITIONAL_FIELDS + 1
  integer, private, parameter :: NO_SUCH_REFERENCE = NO_SUCH_FIELD + 1
  integer, private, parameter :: NOT_FIELD_OF = NO_SUCH_REFERENCE + 1
  integer, private, parameter :: NOT_FUNC = NOT_FIELD_OF + 1
  integer, private, parameter :: NOT_LIT_OF_TYPE = NOT_FUNC + 1
  integer, private, parameter :: NOT_LOG_VALUE = NOT_LIT_OF_TYPE + 1
  integer, private, parameter :: NOT_NAME = NOT_LOG_VALUE + 1
  integer, private, parameter :: NOT_NAME_OR_STRING = NOT_NAME + 1
  integer, private, parameter :: NOT_NUMERIC = NOT_NAME_OR_STRING + 1
  integer, private, parameter :: NOT_SECTION = NOT_NUMERIC + 1
  integer, private, parameter :: NOT_SPEC = NOT_SECTION + 1
  integer, private, parameter :: NOT_STRING = NOT_SPEC + 1
  integer, private, parameter :: NOT_UNITLESS = NOT_STRING + 1
  integer, private, parameter :: NOT_UNITS = NOT_UNITLESS + 1
  integer, private, parameter :: OUT_OF_PLACE = NOT_UNITS + 1
  integer, private, parameter :: SECTION_ORDER = OUT_OF_PLACE + 1
  integer, private, parameter :: VARIABLE_CONFLICT = SECTION_ORDER + 1
  integer, private, parameter :: WRONG_ARG_TYPE = VARIABLE_CONFLICT + 1
  integer, private, parameter :: WRONG_EXPR_TYPE = WRONG_ARG_TYPE + 1
  integer, private, parameter :: WRONG_NUM_ARGS = WRONG_EXPR_TYPE + 1
  integer, private, parameter :: WRONG_TYPE = WRONG_NUM_ARGS + 1
  integer, private, parameter :: WRONG_UNITS = WRONG_TYPE + 1

  logical, private :: ALL_FIELDS_FLAG   ! All fields are required
  integer, private :: CURRENT_SECTION = section_first - 1
  logical, private :: GOT(field_first:field_last)
  logical, private :: NO_DUP_FLAG       ! Duplicate named fields prohibited

  ! Do_Construct_Stack is for checking whether EXIT and CYCLE are inside of
  ! the DO construct named in the label, or for getting the nearest DO
  ! construct if there is no label.  Its value is the root of a DO construct.
  integer, private, allocatable :: Do_Construct_Stack(:)
  integer, private :: Do_Stack_Top

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

! ---------------------------------------------------  CHECK_TREE  -----
  subroutine CHECK_TREE ( ROOT, ERROR_FLAG, FIRST_SECTION, &
    & HOW_MANY_SECTIONS  )
  ! Traverse the abstract syntax tree starting at ROOT, which should
  ! be a N_CFS node (but we don't check).
    integer, intent(in) :: ROOT               ! Root of the tree
    integer, intent(out), optional :: ERROR_FLAG    ! /=0 means trouble
    integer, intent(out), optional :: FIRST_SECTION ! Which son of root is
                                                    ! first CF, IF or SELECT?
    integer, intent(out), optional :: HOW_MANY_SECTIONS ! Number of begin ... ends

    integer :: I              ! Loop inductor
    integer :: Me = -1        ! String index for trace
    integer :: MyFirst
    integer :: SON            ! Son of root

    do_stack_top = 0
    error = 0
    myFirst = -1
    num_sections = 0
    if ( present(first_section) ) first_section = 0
    call trace_begin ( me, 'CHECK_TREE', root, cond=toggle(con) )
    do i = 1, nsons(root)
      son = subtree(i,root)
      select case ( node_id(son) )
      case ( n_cf )         ! A begin-end block's subtree
        if ( myFirst < 0 ) myFirst = i
        call one_section ( son )
      case ( n_cycle, n_exit); call cycle_exit ( root, +1 ) ! +1 => Enclosing sections
      case ( n_do );        call do_construct ( son )
                            if ( myFirst < 0 ) myFirst = i
      case ( n_dt_def );    call def_type ( son ) ! Declare lits in the type
      case ( n_func_def );  call def_func ( son ) ! Declare function
      case ( n_if );        call if_construct ( son )
                            if ( myFirst < 0 ) myFirst = i
      case ( n_section );   call def_section ( son )
      case ( n_spec_def );  call def_spec ( son )
      case ( n_select );    call select_construct ( son )
                            if ( myFirst < 0 ) myFirst = i
      case ( n_variable );  call variable_def ( son )
                            if ( myFirst < 0 ) myFirst = i
      case ( n_while );     call while_construct ( son )
                            if ( myFirst < 0 ) myFirst = i
      case default ;        call announce_error ( son, no_code_for )
      end select
    end do
    if ( present(error_flag) ) error_flag = error
    if ( present(how_many_sections) ) how_many_sections = num_sections
    call trace_end ( 'CHECK_TREE', cond=toggle(con) )
    if ( present(first_section) ) first_section = myFirst
  end subroutine CHECK_TREE

  ! -------------------------------------------------  Check_Type  -----
  logical function Check_Type ( Type_Index, Lit_Index )
  ! Return answer to "Is Lit_Index a literal of type Type_Index?"
  ! This function is not used by other procedures in this module.
  ! It is here instead of in declaration_table because this module has
  ! access to data_type_indices and lit_indices from init_tables_module,
  ! no matter which init_tables_module is used.

    integer, intent(in) :: Type_Index
    integer, intent(in) :: Lit_Index

    type(decls) :: Lit_Decl, Type_Decl

    ! This presumably succeeds:
    type_decl = get_decl(data_type_indices(type_index), type_name)
    ! This presumably also succeeds:
    lit_decl = get_decl(lit_indices(lit_index), enum_value)
    check_type = .true.
    do
      if ( lit_decl%tree == type_decl%tree ) &
  return
      if ( lit_decl%prior == null_decl ) exit
      lit_decl = prior_decl(lit_decl, enum_value)
    end do
    check_type = .false.
  end function Check_Type

! =====     Private Procedures     =====================================
! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, SONS, FIELDS, EXPECT, GOT  )
    use Intrinsic, only: PHYQ_Indices
    use Tree, only: Node_Kind, Pseudo, Where_At => Where
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                   ! "where".  If they're pseudo_terminal,
                                   ! their declarations are dumped.
    integer, intent(in), optional :: FIELDS(:) ! Field indices
    integer, intent(in), optional :: EXPECT    ! Something expected
    integer, intent(in), optional :: GOT       ! Got this instead of expected
    integer :: I                   ! Index for "sons" or "section_ordering"
                                   ! or subtrees of "sons" or subtree of "expect"

    error = max(error,1)
    call startErrorMessage ( where )
    select case ( code )
    case ( already_declared )
      call dump_tree_node ( where, 0 )
      call output ( ' is already defined.', advance='yes' )
    case ( array_not_allowed )
      call display_string ( field_indices(fields(1)), &
        & before='an array value is not allowed for the "' )
      call output ( '" field.', advance='yes' )
    case ( empty_not_allowed )
      call display_string ( field_indices(fields(1)), &
        & before='an empty value is not allowed for the "' )
      call output ( '" field.', advance='yes' )
    case ( improper_cycle )
      call output ( 'CYCLE does not refer to a DO or WHILE construct.', &
        & advance='yes' )
    case ( improper_cycle_exit )
      call output ( 'CYCLE or EXIT within a section cannot refer to' )
      call output ( ' a construct enclosing a section.', advance='yes' )
    case ( inconsistent_data_types )
      call output ( 'data types are not consistent' )
      if ( present(expect) ) &
        & call display_string ( data_type_indices(expect), before=', expected ' )
      if ( present(got) ) &
        & call display_string ( data_type_indices(got), before=', got ' )
      call output ( '.', advance='yes' )
    case ( inconsistent_types )
      call output ( 'types are not consistent' )
      if ( present(expect) ) &
        & call display_string ( data_type_indices(expect), before=', expected ' )
      if ( present(got) ) &
        & call display_string ( data_type_indices(got), before=', got ' )
      call output ( '.', advance='yes' )
    case ( inconsistent_units )
      call output ( 'units are not consistent.', advance = 'yes' )
      if ( present(expect) ) &
        & call display_string ( phyq_indices(expect), before=', expected ' )
      if ( present(got) )&
        & call display_string ( phyq_indices(got), before=', expected ' )
    case ( label_conflict )
      call display_string ( got, before='A label ' )
      call output ( ' shall not be the same as an enumeration literal,' )
      call output ( ' a named value, or a variable.', advance='yes' )
    case ( missing_field )
      call display_string ( field_indices(fields(1)), before='the "' )
      call output ( '" field is required but not present.', advance='yes' )
    case ( no_array_array )
      call display_string ( field_indices(fields(1)), &
        & before='an array element is not allowed to be an array in the "' )
      call output ( '" field.', advance='yes' )
    case ( no_code_for )
      call output ( 'there is no code to analyze ' )
      call dump_tree_node ( where, 0, advance='yes' )
    case ( no_cycle_exit_target )
      call output ( 'there is no target for the CYCLE or EXIT statement.', &
        & advance='yes' )
    case ( no_declaration )
      call display_string ( sub_rosa(where), before='the symbol ' )
      call output ( ' does not appear as a label of a specification', advance='yes' )
    case ( no_dot )
      call output ( 'a reference of the form X.Y is not allowed.', &
        advance='yes' )
    case ( no_duplicate_fields )
      call display_string ( field_indices(fields(1)), before='the "' )
      call output ( '" field shall not be specified twice.', advance='yes' )
    case ( no_expr )
      call display_string ( field_indices(fields(1)), &
        & before='an expression value is not allowed for the "' )
      call output ( '" field.', advance='yes' )
    case ( no_positional_fields )
      call output ( 'positional fields are not allowed.', advance='yes' )
    case ( no_such_field )
      call display_string ( field_indices(fields(1)), before='a required field "' )
      call output ( '" is absent in the chain of specifications.', &
        advance='yes' )
    case ( no_such_reference )
      call display_string ( sub_rosa(where), before='there is no reference to ' )
      call output ( ' in the field at ' )
      if ( present(sons) ) call print_source ( where_at(sons(1)), advance='yes' )
    case ( not_field_of )
      call display_string ( sub_rosa(where) )
      call display_string ( sub_rosa(subtree(1,fields(1))), &
        & before=' is not a field of ', advance='yes' )
    case ( not_func )
      call display_string ( fields(1) )
      call output ( ' is not a valid function.', advance='yes' )
    case ( not_log_value )
      call output ( 'is not boolean.', advance = 'yes' )
    case ( not_name )
      call output ( 'is not a name.', advance = 'yes' )
    case ( not_name_or_string )
      call output ( 'is not a name or a string.', advance = 'yes' )
    case ( not_numeric )
       call output ( 'is not numeric.', advance = 'yes' )
    case ( not_section )
      call display_string ( sub_rosa(where) )
      call output ( ' is not a section name.', advance = 'yes' )
    case ( not_spec )
      if ( present(expect) ) then
        call display_string ( expect )
      else
        call display_string ( sub_rosa(where) )
      end if
      call output ( ' is not a spec name.', advance = 'yes' )
    case ( not_string )
      call output ( 'is not a string or of logical type.', &
                    advance = 'yes' )
    case ( not_unitless )
      if ( present(fields) ) then
        call display_string ( fields(1), before='Argument of ' )
      else
        call output ( 'Operand' )
      end if
      call output ( ' is not unitless.', advance='yes' )
    case ( not_units )
      call dump_tree_node ( where, 0 )
      call output ( ' is not a units name.', advance = 'yes' )
    case ( out_of_place )
      call display_string ( sub_rosa(where) )
      call display_string ( sub_rosa(sons(1)), before=' is not allowed in a ' )
      call output ( ' section.', advance='yes' )
    case ( section_order )
      call output ( 'Section '); call display_string ( sub_rosa(where) )
      call output ( ' is out of order.', advance='yes' )
      call output ( '***** Expected' )
      do i = section_first, section_last
        if ( section_ordering(i,current_section) /= 0 ) then
          call output ( ' ' ); call display_string ( section_indices(i) )
        end if
      end do
      call newLine
    case ( variable_conflict )
      call display_string ( got, before='A variable name ' )
      call output ( ' shall not be the same as an enumeration literal,' )
      call output ( ' a label, or a named value.', advance='yes' )
    case ( wrong_arg_type )
      call display_string ( fields(1), before='Argument of ' )
      call output ( ' is not the correct type' )
      if ( present(expect) ) &
        & call display_string ( sub_rosa(expect), before=', expected ' )
      call output ( '.', advance='yes' )
    case ( wrong_expr_type )
      call display_string ( data_type_indices(expect), &
        & before='Expression type is not correct.  Expected "' )
      call output ( '"', advance='yes' )
    case ( wrong_num_args )
      call display_string ( fields(1), before='Incorrect number of arguments for ', &
        & advance='yes' )
    case ( wrong_type, wrong_units )
      i = node_kind(where)
      if ( present(fields) ) i = node_kind(fields(1))
      if ( i == pseudo ) then
        call output ( 'the "' )
        if ( present(fields) ) then
          call display_string ( sub_rosa(fields(1)) )
        else
          call display_string ( sub_rosa(where) )
        end if
        call output ( '" field has the wrong ' )
      else
        call output ( 'the field has the wrong ' )
      end if
      if ( code == wrong_type ) then
        call output ( 'type of associated value' )
        if ( present(expect) ) then
          if ( expect == 0 ) then
            call output ( '.', advance='yes' )
          else
            call output ( ', expected ' )
            if ( node_id(expect) /= n_or ) then
              call theType ( expect, 2 )
            else
              call output ( 'one of' )
              do i = 2, nsons(expect)
                call output ( ' ' )
                call theType ( subtree(i,expect), 1 )
              end do
            end if
            call newLine
          end if
        else
          call output ( '.', advance='yes' )
        end if
      else
        call display_string ( phyq_indices(expect), before='units, ' )
        call output ( ' expected.', advance='yes' )
      end if
      if ( present(sons) ) then
        call output ( '         Expected' )
        if ( nsons(sons(1)) > 2 ) call output ( ' one of' )
        do i = 2, nsons(sons(1))
          if ( i > 2 ) call output ( ',' )
          call output ( ' ' )
          if ( i > 2 .and. i == nsons(sons(1)) ) call output ( 'or ' )
          call display_string ( sub_rosa(subtree(i,sons(1))) )
        end do
        call newLine
      end if
    case default
      call output ( code, before='No message in TREE_CHECKER for error code ', &
        & advance='yes' )
      stop
    end select
    if ( present(sons) ) then
      do i = 1, size(sons)
        if ( node_kind(sons(i)) == pseudo ) &
          call dump_1_decl ( sub_rosa(sons(i)) )
      end do
    end if
  contains
    subroutine TheType ( Where, Start )
      integer, intent(in) :: Where ! The type tree node
      integer, intent(in) :: Start ! for index for sons of Where
      integer :: I
      do i = start, nsons(where)
        call display_string ( sub_rosa(subtree(i,where)) )
        if ( i /= nsons(where) ) call output ( '.' )
      end do
    end subroutine TheType
  end subroutine ANNOUNCE_ERROR

! ----------------------------------------------------  ASSIGN  -----
  subroutine ASSIGN ( ROOT, TYPE, UNITS, VALUE )
  ! Analyze a son of "n_spec_args" of the form "name = expr+",
  ! starting at ROOT

    integer, intent(in) :: ROOT            ! Index of the n_asg node
    integer, intent(out) :: TYPE, UNITS    ! Output from "expr"
    double precision, intent(out) :: VALUE ! Output from "expr"

    integer :: F         ! Index of son, grandson of Field
    integer :: FIELD     ! Tree node of field's declaration
    integer :: FIELD_LIT ! f_... for a field
    integer :: FIELD_LOOK ! A field being sought during n_dot checking in AssignBody
    integer :: FIELD_TEST ! A son of Field_Ref in AssignBody
    integer :: I         ! Index of son of "root"
    integer :: LOOK_FOR  ! Look for an enum_value or a spec?
    integer :: MaxStat   ! of sons of OR
    integer :: Me = -1   ! String index for trace
    logical :: NO_ARRAY_ALLOWED ! Field is n_field_type and is required to be scalar
    logical :: OK        ! At least one of sons of OR was OK
    integer :: SON1, SON ! Sons of "root"
    integer :: SPEC_DECL ! Tree node of the spec's declaration
    integer :: Stat      ! From AssignBody

    call trace_begin ( me, 'ASSIGN', root, cond=toggle(con) )
    son1 = subtree(1,root)
    if ( node_id(son1) == n_identifier ) then
      spec_decl = decoration(root)
      field = check_field(son1,spec_decl)   ! Is field a field of spec?
      if ( field == 0 ) then
        call announce_error ( son1, not_field_of, &
          & fields=(/ spec_decl /) )
      else
        no_array_allowed = mod(decoration(field)/no_array,2) /= 0
        if ( node_id(field) == n_field_type ) then
          look_for = enum_value
        else if ( node_id(field) == n_variable_ref ) then
          look_for = variable
        else ! node_id(field) == n_field_spec or node_id(field) == n_dot
          look_for = label
        end if
        call decorate ( root, field )
        field_lit = decoration(subtree(1,field))
        if ( no_dup_flag ) then
          if ( got(field_lit) ) &
            & call announce_error ( root, no_duplicate_fields, &
            & fields=(/ field_lit /) )
        end if
        got(field_lit) = .true.
        call decorate ( son1, field_lit )
        if ( mod(decoration(field)/empty_ok,2) == 0 .and. nsons(root) == 1 ) &
          & call announce_error ( root, empty_not_allowed, fields=(/ field_lit /) )
        if ( node_id(field) /= n_unchecked ) then
          do i = 2, nsons(root)
            son = subtree(i,root)
            if ( node_id(field) /= n_or ) then
              ! First son of field is the field name; the rest are types.
              stat = assignBody ( son, field, 2 )
            else
              ok = .false.
              maxStat = 0
              do f = 2, nsons(field)
                if ( node_id(subtree(f,field)) == n_unchecked ) then
                  ok = .true.
                  stat = 0
              cycle
                end if
                ! First son of n_or is the field name.  All sons of n_or
                ! are types.
                no_array_allowed = no_array_allowed .or. &
                                   mod(decoration(subtree(f,field))/no_array,2) /= 0
                stat = assignBody ( son, subtree(f,field), 1 )
                maxStat = max(stat,maxStat)
                if ( stat == 0 ) then
                  ok = .true.
              cycle
                end if
              end do
              stat = maxStat
              if ( ok ) stat = 0
            end if
            select case ( stat )
            case ( 0 )
            case ( no_array_array )
              call announce_error ( subtree(1,son), no_array_array, fields=(/ field_lit /) )
            case ( no_declaration )
              call announce_error ( subtree(1,son), no_declaration )
            case ( no_expr )
              call announce_error ( subtree(2,son), no_expr, fields=(/ field_lit /) )
            case ( no_such_field )
              call announce_error ( subtree(2,son), no_such_field, &
                & fields=(/ field_look /) )
            case ( no_such_reference )
              call announce_error ( subtree(2,son), no_such_reference, &
                & (/ subtree(1,field_test) /) )
            case ( not_spec )
              call announce_error ( subtree(2,son), not_spec, expect=field_last )
            case ( wrong_type )
              call announce_error ( son, wrong_type, fields=(/son1/), &
                & expect=field )
            case default
              call announce_error ( son1, wrong_units, fields = (/ son1 /), &
                & expect=-stat )
            end select
          end do ! i = 2, nsons(root)
          if ( no_array_allowed .and. .not. scalar(root) ) &
            & call announce_error ( root, array_not_allowed, fields=(/ field_lit /) )
        end if
      end if
    else
      call announce_error ( son1, not_name )
    end if
    call trace_end ( 'ASSIGN', cond=toggle(con) )

  contains

    recursive integer function AssignBody ( Root, Field, Start ) result ( Stat )
      integer, intent(in) :: Root   ! of son of n_asg or n_array
      integer, intent(in) :: Field  ! spec of field to check
      integer, intent(in) :: Start  ! index of sons of field

      logical :: Array_Array_OK ! OK for array element to be an array
      type(decls) :: DECL  ! Declaration of a name
      integer :: Gson      ! Son of son
      integer :: J         ! Index of root of "field"
      integer :: Me = -1   ! String index for trace
      integer :: TEST_TYPE ! Used to test the tree-node for a type reference
      integer :: TYPE_DECL ! Tree node of declaration of name of field's type

      call trace_begin ( me, 'AssignBody', root, index=field, cond=toggle(con) )

      ! Stat = 0 => Type and units OK
      !      else wrong type
      !      or   no such field
      !      or   no such reference
      !      or   no such label
      !      < 0 => Expecting units |stat|

      stat = 0 ! Assume normal return status
      select case ( node_id(root) )
      case ( n_dot ) ! label_ref.field
        if ( node_id(field) == n_dot ) then      ! dot allowed
          stat = check_dot ( root, field, start, field_look, field_test )
        else
          stat = wrong_type ! Wrong type
        end if
      case ( n_identifier )
        if ( node_id(field) == n_variable_ref ) then
          decl = get_decl(sub_rosa(root),look_for)
          if  ( decl%type /= null_decl ) &
    go to 9
        else if ( node_id(field) /= n_dot ) then ! Field doesn't need x.y
          do j = start, nsons(field)             ! Try all of the types
            type_decl = decoration(subtree(j,field))
            decl = get_decl(sub_rosa(root),[look_for,variable])
            do while ( decl%tree /= null_tree )
              if ( decl%type==variable ) then
                ! Type of variable is not the tree that represents the type,
                ! so check for the type, not the tree that represents it.
                if ( node_id(type_decl) == n_dt_def ) &
                  & type_decl = decoration(subtree(1,type_decl))
                stat = expr (root, type, units, value, field, start, field_look, field_test)
                if ( type == enum_value ) then
                  test_type = units ! The enumeration type
                else
                  test_type = type
                end if
              else
                test_type = decl%tree
                select case ( decl%type )
                case ( do_label, enum_value, function, label, named_value, &
                     & type_name ) ! decl%tree is a tree node index
                  if ( node_id(type_decl) == n_spec_def ) &
                     & test_type = decoration(subtree(1,decl%tree))
                end select
              end if
              if ( test_type == type_decl ) then  ! right type
                if ( look_for == enum_value ) then
                  call decorate ( root, decl%units ) ! decorate root with lit#
                else
                  call decorate ( root, decl%tree )  ! decorate root with tree
                end if
    go to 9
              end if
              decl = prior_decl(decl,[look_for,variable])
            end do
          end do
        end if
        ! Either field needs x.y, or an allowed type was not found
        stat = wrong_type
      case ( n_array )
        if ( mod(decoration(field)/array_array,2) == 0 ) then
            stat = no_array_array
    go to 9
        end if
        do j = 1, nsons(root)
          gson = subtree(j,root)
          stat = assignBody ( subtree(j,root), field, start )
          if ( stat /= 0 ) &
    go to 9
        end do
      case default
        stat = expr (root, type, units, value, field, start, field_look, field_test)
        if ( node_id(field) == n_dot ) then ! Field needs x.y
          stat = merge(no_expr,0,mod(decoration(field)/expr_ok,2) == 0)
          if ( stat == 0 ) stat = merge(0,wrong_type,type==dot)
        else ! Field doesn't need x.y; check the type
          if ( type == enum_value ) then
            test_type = units ! The enumeration type
          else
            test_type = type
          end if
          stat = check_field_type(field, test_type, units, start)
        end if
      end select
    9 call trace_end ( 'AssignBody', stat, cond=toggle(con) )
    end function AssignBody

  end subroutine ASSIGN

! ----------------------------------------------------  Check_Dot  -----
  integer function Check_Dot ( Root, Field, Start, Field_Look, Field_Test ) &
    & result ( Stat )

    ! Check label.field
    integer, intent(in) :: Root  ! n_dot node of tree to check
    integer, intent(in) :: Field ! Declaration of the field types
    integer, intent(in) :: Start ! of sons of Field
    integer, intent(out) :: Field_Look ! expected f_field index, for error message
    integer, intent(out) :: Field_Test ! parent of checked field, for error message
    type(decls) :: Decl   ! Declaration of a name
    integer :: Field_Got  ! The field that was the desired one
    integer :: Field_Last ! Sub_rosa of Last
    integer :: Field_Ref  ! The value of a field -- a ref to a n_spec_arg
    integer :: First      ! label in <dot label field> son of root
    integer :: Last       ! field in <dot label field> son of root
    integer :: Me = -1    ! String index for trace
    integer :: Test_Type  ! Used to test the tree-node for a type reference
    integer :: Type_Decl  ! Tree node of declaration of name of field's type

    call trace_begin ( me, 'Check_Dot', root, cond=toggle(con), advance='no' )
    if ( toggle(con) ) then
      call output ( field, before=' Field = ', advance='no' )
      call output ( start, before=' Start = ', advance='yes' )
    end if

    stat = no_such_reference ! Assume no label referenced in <dot label field>
    type_decl = decoration(subtree(start,field)) ! Required spec
    first = subtree(1,root) ! label in <dot label field>
    last = subtree(2,root)  ! field in <dot label field>
    field_last = sub_rosa(last)
    decl = get_decl(sub_rosa(first),label)
    if ( decl%tree /= null_tree ) then
      do while ( decl%tree /= null_tree )
        test_type = decoration(subtree(1,decl%tree)) ! spec's index
        if ( test_type == type_decl ) then  ! right kind of spec
          field_ref = decl%tree
          call decorate ( first, field_ref ) ! decorate label_ref with tree
          stat = check_deep ( start+1, nsons(field), field_ref )
          if ( stat == 0 ) exit
        end if
        decl = prior_decl(decl,label)
      end do
    else
      stat = no_declaration
      field_look = -1 ! The name before the dot doesn't exist
      field_test = -1 ! The name before the dot doesn't exist
    end if
    call trace_end ( 'Check_Dot', stat, cond=toggle(con) )

  contains

    recursive integer function Check_Deep ( Start, End, Field_Ref ) &
      & result ( Stat )
      integer, intent(in) :: Start      ! of sons of Field
      integer, intent(in) :: End        ! of sons of Field
      integer, intent(in) :: Field_Ref  ! spec_arg tree to check against son
                                        ! Start of Field
      integer :: Field_Chk  ! index of n_asg vertex
      integer :: K     ! Index of sons of Field_Ref, to find N_Asg vertex
      integer :: L     ! Index of sons of Field_Test, to find next spec_arg
      integer :: Me = -1    ! String index for trace
      integer :: Name_Look  ! sub_rosa of Start son of Field

      call trace_begin ( me, 'Check_Deep', subtree(start,field), field_ref, &
        & cond=toggle(con) )

      stat = 0 ! Assume success
      field_look = decoration(subtree(start,field))
      name_look = sub_rosa(subtree(start,field))
      do k = 2, nsons(field_ref)
        field_chk = subtree(k,field_ref)
        field_test = field_chk
        ! Check for n_asg son of field_ref
        if ( node_id(field_chk) == n_asg ) then
          if ( sub_rosa(subtree(1,field_chk)) == name_look ) then
            ! The field in the spec_arg is the same as the Start son of Field
            field_got = field_chk
            if ( start == end ) then
              ! Check whether a son of the n_asg at field_chk is the field in
              ! <dot label field> to be checked
              do l = 2, nsons(field_chk)
                if ( sub_rosa(subtree(l,field_got)) == field_last ) then
                  call decorate ( last, subtree(l,field_got) )
                  go to 9
                end if
              end do
              field_test = field_last
              stat = not_spec ! No such referenced spec
              go to 9
            else
              ! Check whether a son of the n_asg at field_chk is another field
              ! of the same name as the Start son of Field
              do l = 2, nsons(field_chk)
                stat = check_deep ( start, end, decoration(subtree(l,field_chk)) )
                if ( stat == 0 ) go to 9
              end do
              ! Check whether a son of the n_asg at field_chk is
              ! the same name as the Start+1 son of Field
              do l = 2, nsons(field_chk)
                stat = check_deep ( start+1, end, decoration(subtree(l,field_chk)) )
                if ( stat == 0 ) go to 9
              end do
            end if
          end if
        end if
      end do ! k
      stat = no_such_field ! no such field
    9 call trace_end ( 'Check_Deep', cond=toggle(con) )
    end function Check_Deep

  end function Check_Dot

! --------------------------------------------------  CHECK_FIELD  -----
  integer function CHECK_FIELD ( FIELD, SPEC )
  ! Return position in tree of "n_field_type" or "n_field_spec" node of
  ! field definition if "field" is a field of "spec," else return zero.

    integer, intent(in) :: FIELD
    integer, intent(in) :: SPEC
    integer :: I         ! Index of sons of "spec_def"
    integer :: Me = -1   ! String index for trace cacheing
    integer :: SUB_FIELD ! Sub_rosa = string index of "field"

    call trace_begin ( me, 'CHECK_FIELD', field, cond=toggle(con) )
    sub_field = sub_rosa(field)    ! String index of field
    do i = 2, nsons(spec)
      check_field = subtree(i,spec)
      if ( sub_rosa(subtree(1,check_field)) == sub_field ) go to 9
    end do
    check_field = 0
  9 call trace_end ( 'CHECK_FIELD', check_field, cond=toggle(con) )
  end function CHECK_FIELD

! ---------------------------------------------  CHECK_FIELD_TYPE  -----
  integer function CHECK_FIELD_TYPE ( FIELD, TYPE, UNITS, START )
  ! Return zero if 'type' is allowed for 'field' of 'spec' "
  ! or if 'type' is allowed for parameter declared at 'field' ".
  ! Return Wrong_Type if type is wrong.
  ! Return -(required PHYQ_...) if type is correct but units are wrong.
    integer, intent(in) :: FIELD   ! Field declaration, from "check_field"
                                   ! or a parameter declaration from
                                   ! "check_section"
    integer, intent(in) :: TYPE    ! A "t_type" name from Init_Tables_Module
    integer, intent(in), optional :: UNITS ! From EXPR
    integer, intent(in), optional :: START
    type(decls) :: DECL            ! Declaration of "type"
    integer :: I                   ! Index of sons of "spec_def"
    integer :: Me = -1             ! String index for trace
    integer :: MyStart             ! Either 2 or start
    integer :: REQ_U               ! Required units of a son of FIELD

    call trace_begin ( me, 'Check_Field_Type', field, cond=toggle(con) )
    myStart = 2
    if ( present(start) ) myStart = start
    check_field_type = merge( wrong_type, 0, &
      &  type < lbound(data_type_indices,1) .or. &
      &  type > ubound(data_type_indices,1) )
    if ( check_field_type /= 0 ) &
  go to 9
    decl = get_decl(data_type_indices(type), type_name)
    req_u = decoration(field) / u
    do i = myStart, nsons(field)
      check_field_type = merge( 0, wrong_type, &
                              & decoration(subtree(i,field)) == decl%tree )
      if ( check_field_type == 0 ) then
        if ( req_u > 0 .and. present(units) ) then
          if ( req_u /= units ) check_field_type = -req_u
        end if
  go to 9
      end if
    end do
    check_field_type = wrong_type
  9 continue
    call trace_end ( 'Check_Field_Type', check_field_type, cond=toggle(con) )
  end function CHECK_FIELD_TYPE

! --------------------------------------------------  Check_Label  -----
  logical function Check_Label ( Root ) result ( OK )
  ! Check that a label definition is not already a label, do_label, named_value,
  ! or variable.
    integer, intent(in) :: Root    ! Son of N_Named vertex

    type(decls) :: Decl            ! Declaration of a label
    integer :: Me = -1             ! String index for trace
    integer :: String1             ! Text of Gson1, which is a label

    call trace_begin ( me, 'Check_Label', root, cond=toggle(con) )
    ok = .false.
    string1 = sub_rosa(root)
    decl = get_decl ( string1, [ do_label, enum_value, label, named_value, &
                                & variable ] )
    select case ( decl%type )
    case ( label, do_label )
      call announce_error ( root, already_declared )
    case ( named_value, variable )
      call announce_error ( root, label_conflict, got=string1 )
    case default
      ok = .true.
    end select
    call trace_end ( 'Check_Label', cond=toggle(con) )

  end function Check_Label

! ------------------------------------------------  CHECK_SECTION  -----
  integer function CHECK_SECTION ( SPEC_CHECK, SECTION )
  ! Return subtree of "section" where "spec_check" is found, else zero.
  ! "spec_check" is a pseudo-terminal name of a spec name or a
  ! parameter name.
    integer, intent(in) :: SPEC_CHECK
    integer, intent(in) :: SECTION
    integer :: I           ! Index of sons of "sect_def"
    integer :: SUB_SPEC    ! Sub_rosa = string index of 

    sub_spec = sub_rosa(spec_check)    ! String index of spec_check
    do i = 2, nsons(section)           ! specs allowed in section
      check_section = subtree(i,section)
      if ( node_kind(check_section) == pseudo ) then
        if ( sub_rosa(check_section) == sub_spec ) &
  return
      else ! Must be n_name_def
        if ( sub_rosa(subtree(1,check_section)) == sub_spec ) &
  return
      end if
    end do
    check_section = 0
  end function CHECK_SECTION

! ---------------------------------------------------  Cycle_Exit  -----
  subroutine Cycle_Exit ( Root, In_Or_Out )
  ! If a CYCLE or EXIT statement has a label, make sure the label is the
  ! label of an enclosing DO construct.  If it doesn't have a label, make
  ! sure the statement is within at least one DO construct.  Decorate the
  ! root with the DO construct.
    integer, intent(in) :: Root      ! Of the CYCLE or EXIT statement
    integer, intent(in) :: In_Or_Out ! +1 => Outside a section, -1 => Inside

    integer :: Do_Tree               ! The referenced DO construct
    integer :: I                     ! For hunting in the DO construct stack
    integer :: Me = -1               ! String index for trace cacheing
    integer :: Son                   ! of a DO construct root

    call trace_begin ( me, 'Cycle_Exit', root, cond=toggle(con) )

    if ( do_stack_top <= 0 ) then
      call announce_error ( root, no_cycle_exit_target )
    else
      if ( nsons(root) == 0 ) then ! EXIT or CYCLE that doesn't name a construct
        do_tree = 0
        ! Hunt for DO or WHILE construct
        do i = do_stack_top, 1, -1
          if ( node_id(do_construct_stack(i)) == n_do .or. &
             & node_id(do_construct_stack(i)) == n_while ) then
            do_tree = do_construct_stack(i)
            exit
          end if
        end do
      else                         ! EXIT or CYCLE that does name a construct
        do_tree = 0
        ! Hunt for the label on DO constructs in the stack
        do i = do_stack_top, 1, -1
          son = subtree(1,do_construct_stack(i))
          if ( node_id(son) == n_named ) then
            if ( sub_rosa(subtree(1,son)) == sub_rosa(subtree(1,root)) ) then
              do_tree = do_construct_stack(i)
              exit
            end if
          end if
        end do
      end if
      if ( do_tree == 0 ) then
        call announce_error ( root, no_cycle_exit_target )
      else
        if ( node_id(root) == n_cycle .and. node_id(do_tree) /= n_do .and. &
             node_id(do_tree) /= n_while ) then
          call announce_error ( root, improper_cycle )
        else if ( decoration(do_tree) /= in_or_out ) then
          call announce_error ( root, improper_cycle_exit )
        else
          call decorate ( root, do_tree )
        end if
      end if
    end if

    call trace_end ( 'Cycle_Exit', cond=toggle(con) )

  contains

    subroutine Check_In_Out ( Do_Tree )
      integer, intent(in) :: Do_Tree
    end subroutine Check_In_Out

  end subroutine Cycle_Exit

! -----------------------------------------------------  DEF_FUNC  -----
  subroutine DEF_FUNC ( ROOT )
  ! Process a definition of a function: Enter it into the declaration
  ! table.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_dt_def )

    integer :: Me = -1   ! String index for trace
    integer :: SON       ! Son of Root

    call trace_begin ( me, 'DEF_FUNC', root, cond=toggle(con) )
    son = subtree(1,root)
    !              String         Value                Type
    call declare ( sub_rosa(son), sub_rosa(son)+0.0d0, function, &
    !              Units            Tree
                 & decoration(son), root )
    call trace_end ( 'DEF_FUNC', cond=toggle(con) )
  end subroutine DEF_FUNC

! --------------------------------------------------  DEF_SECTION  -----
  subroutine DEF_SECTION ( ROOT )
  ! Process a definition of a section:  Enter the section name in the
  ! declaration.  Decorate the section name with the section definition.
  ! If the section has name definitions, enter them in the declaration
  ! table, with the "units" field equal to the decoration of the type
  ! name (second son of the "name_def" node) and the "tree" field indexing
  ! the "name_def" node.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_dt_def )

    type(decls) :: DECL  ! Declaration of a parameter type
    integer :: GSON      ! The type name node for a parameter declaration
    integer :: I         ! Loop inductor
    integer :: Me = -1   ! String index for trace
    integer :: SON       ! Son of Root

    call trace_begin ( me, 'DEF_SECTION', root, cond=toggle(con) )
    son = subtree(1,root)
                 ! String         Value                Type
    call declare ( sub_rosa(son), sub_rosa(son)+0.0d0, section, &
                 ! Units            Tree
                 & decoration(son), root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( node_id(son) == n_name_def ) then
        gson = subtree(1,son)
                     ! String          Value                 Type
        call declare ( sub_rosa(gson), sub_rosa(gson)+0.0d0, named_value, &
                     ! Units            Tree
                       decoration(gson), son )
        gson = subtree(2,son)
        decl = get_decl(sub_rosa(gson),type_name)
        call decorate ( gson, decl%tree )    ! the dt_def
      end if
    end do
    call trace_end ( 'DEF_SECTION', cond=toggle(con) )
  end subroutine DEF_SECTION

! -----------------------------------------------------  DEF_SPEC  -----
  recursive subroutine DEF_SPEC ( ROOT )
  ! Process a definition of a specification:  Enter the specification and field
  ! names in the declaration table.  Decorate the specification name with
  ! the root. Decorate each field_type node with the specification name.
  ! Decorate each field name with its parent field_type node.  Decorate
  ! each field type name with the dt_def node for the type.

    integer, intent(in) :: ROOT      ! Root of tree being worked ( n_spec_def )

    type(decls) :: DECL       ! Current declaration of "son"
    integer :: FIELD_NAME, FIELD_TYPE   ! Grandsons
    integer :: I, J           ! Loop inductors
    integer :: Me = -1        ! String index for trace
    integer :: SON            ! I'th son of "root"
    integer :: SPEC_NAME      ! First son of "root"

    call trace_begin ( me, 'DEF_SPEC', root, cond=toggle(con) )
    spec_name = subtree(1,root)
                 ! String               Value                      Type
    call declare ( sub_rosa(spec_name), sub_rosa(spec_name)+0.0d0, spec, &
                 ! Units                  Tree
                 & decoration(spec_name), root )
!   call decorate ( spec_name, root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      son=son ! Without this, the LF95 compiler inexplicably uses "root"
              ! instead of "son" below ???
      field_name = subtree(1,son)
      decl = get_decl(sub_rosa(field_name),field)
      if ( decl%tree == null_tree ) then ! don't make several
                       ! String                Value       Type
        call redeclare ( sub_rosa(field_name), decl%value, field, &
                       ! Units                   Tree
                         decoration(field_name), son )
      else ! need to make several -- one for each
                     ! String                Value       Type
        call declare ( sub_rosa(field_name), decl%value, field, &
                     ! Units                   Tree
                     & decoration(field_name), son )
      end if
      if ( node_id(son) /= n_or ) then
        ! First son is the field name
        call def_one_spec ( son, 2 )
      else
        do j = 2, nsons(son)
          ! First son of n_or is the field name, not first son of son
          call def_one_spec ( subtree(j,son), 1 )
        end do
      end if
    end do
    call trace_end ( 'DEF_SPEC', cond=toggle(con) )

  contains

    subroutine DEF_ONE_SPEC ( SON, START )
      integer, intent(in) :: SON   ! of spec_def or n_or
      integer, intent(in) :: START ! of sons of SON
      integer :: J                 ! Loop inductor
      integer :: Me = -1           ! String index for trace
      call trace_begin ( me, 'DEF_ONE_SPEC', root, cond=toggle(con) )
      select case ( node_id(son) )
      case ( n_field_type )
        do j = start, nsons(son)
          field_type = subtree(j,son)
          decl = get_decl(sub_rosa(field_type),type_name)
          call decorate ( field_type, decl%tree )   ! the dt_def
        end do
      case ( n_field_spec )
        do j = start, nsons(son)
          field_type = subtree(j,son)
          decl = get_decl(sub_rosa(field_type),spec)
          call decorate ( field_type, decl%tree )   ! the spec_def
        end do
      case ( n_dot )
        field_type = subtree(start,son)
        decl = get_decl(sub_rosa(field_type),spec)
        call decorate ( field_type, decl%tree )   ! the spec_def
        ! The rest of the sons are field names, for which new decorations
        ! Won't help -- in fact, the f_field_name's index is best.
      end select
      call trace_end ( 'DEF_ONE_SPEC', cond=toggle(con) )
    end subroutine DEF_ONE_SPEC

  end subroutine DEF_SPEC

! -----------------------------------------------------  DEF_TYPE  -----
  subroutine DEF_TYPE ( ROOT )
  ! Process a definition of a type:  Enter the names in the declaration
  ! table.  The type name and literals are decorated with their indices
  ! in Init_Tables. These decorations are entered into the declaration
  ! table.  Then each one is decorated with Root.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_dt_def )

    integer :: I              ! Loop inductor
    integer :: Me = -1        ! String index for trace
    integer :: SON            ! I'th son of "root"

    call trace_begin ( me, 'DEF_TYPE', root, cond=toggle(con) )
    son = subtree(1,root)
                 ! String         Value                Type
    call declare ( sub_rosa(son), sub_rosa(son)+0.0d0, type_name, &
                 ! Units            Tree
                 & decoration(son), root )
!   call decorate ( son, root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      !              String         Value                Type
      call declare ( sub_rosa(son), sub_rosa(son)+0.0d0, enum_value, &
      !              Units            Tree
                     decoration(son), root )
!     call decorate ( son, root )
    end do
    call trace_end ( 'DEF_TYPE', cond=toggle(con) )
  end subroutine DEF_TYPE

! -------------------------------------------------  Do_Construct  -----
  recursive subroutine Do_Construct ( Root, Son1, Decl )
  ! Analyze a DO construct of the form DO name := exprs cfs END DO,
  ! starting at ROOT
    integer, intent(in) :: Root
    integer, intent(in), optional :: Son1 ! First son of parent CF
    type(decls), intent(in), optional :: DECL ! Declaration of "name" on "begin name"
    integer :: Gson, I, Son, Start, String
    integer :: Me = -1             ! String index for trace
    integer :: Type
    double precision :: Value      ! not used

    call trace_begin ( me, 'DO_Construct', root, cond=toggle(con) )

    call push_do_stack ( root )

    start = 2 ! Assume first son is not n_named
    son = subtree(1,root)
    if ( node_id(son) == n_named ) then
      gson = subtree(1,son)
      if ( check_label(gson) ) then
        string = sub_rosa(gson)
        call declare ( string, 0.0d0+string, do_label, phyq_invalid, root )
      end if
      start = 3
      son = subtree(2,root)
    end if

    call variable_def ( son, type ) ! Treat name := exprs as variable def

    ! If there's more than one expr, the type has to be numeric
    i = nsons(son)
    if ( i > 2 .and. type /= num_value ) call announce_error ( &
      & subtree(1,son), wrong_expr_type, expect=num_value )

    if ( present(son1) ) then
      call decorate ( root, -1 ) ! -1 => inside section
      do i = start, nsons(root)
        gson = subtree(i,root)
        if ( node_id(gson) == n_cf ) then
          call one_section ( gson )
        else
          call one_cf ( gson, son1, decl )
        end if
      end do
    else
      call decorate ( root, +1 ) ! +1 => enclosing sections
      do i = start, nsons(root)
        call one_section ( subtree(i,root) )
      end do
    end if

    do_stack_top = do_stack_top - 1 ! Pop the DO construct stack

    call trace_end ( 'DO_Construct', cond=toggle(con) )

  end subroutine Do_Construct

! --------------------------------------------------------  EQUAL  -----
  subroutine EQUAL ( ROOT )
  ! Analyze a specification of the form NAME = EXPR, starting at ROOT
    integer, intent(in) :: ROOT
    integer :: CHECK          ! At first, subtree where NAME is declared;
                              ! later, output from Check_Field_Type
    type(decls) :: DECL       ! Declaration of son1
    integer :: Me = -1        ! String index for trace
    integer :: SECT           ! The section the parameter appears in
    integer :: SON1, SON2     ! Sons of Root
    integer :: Stat           ! From EXPR
    integer :: TYPE           ! Output from "expr", not otherwise used
    integer :: TYPE_DECL      ! Tree node of declaration of name of param's type
    integer :: UNITS          ! Output from "expr"
    double precision :: VALUE ! Output from "expr", not otherwise used

    call trace_begin ( me, 'EQUAL', root, cond=toggle(con) )
    sect = decoration(root)
    decl = get_decl(sub_rosa(sect), section )   ! Known to be a
                             ! valid section name if we get this far
    son1 = subtree(1,root)
    ! Is spec allowed in section?
    check = check_section(son1,decl%tree)
    if ( check == 0 ) then
      call announce_error ( son1, out_of_place, (/ sect /) )
    else
      decl = get_decl(sub_rosa(son1), named_value)
      call decorate ( son1, decl%units )
      son2 = subtree(2,root)
      if ( node_id(son2) == n_dot ) then
      ! ??? Need some checking here
      else
        if ( node_id(son2) == n_identifier ) then
          decl = get_decl(sub_rosa(son2), enum_value)
          type_decl = decoration(subtree(2,check)) ! decl of allowed type
          do
            if ( decl%tree == null_tree ) then ! no more decls of id
              call announce_error ( son2, wrong_type, fields=(/son1/) )
          exit
            end if
            if ( decl%tree == type_decl ) then ! right type
              call decorate ( son2, decl%units )! decorate son2 with lit#
          exit
            end if
            decl = prior_decl(decl,enum_value)
          end do
        else
          stat = expr( son2, type, units, value )
          stat = check_field_type(check,type,units)
          if ( stat == wrong_type ) then
            call announce_error ( son2, wrong_type, fields=(/son1/), expect=check )
          else if ( stat > 1 ) then
            call announce_error ( son2, wrong_units, fields=(/son1/), expect=-stat )
          end if
        end if
      end if
    end if
    call trace_end ( 'EQUAL', cond=toggle(con) )
  end subroutine EQUAL

! ---------------------------------------------------------  EXPR  -----
  recursive integer function EXPR ( ROOT, TYPE, UNITS, VALUE, FIELD, &
    &                               START, FIELD_LOOK, FIELD_TEST, Tree ) &
    & result ( Stat )

  ! Analyze an expression, check units of everything in it, return its
  ! units.

    use Call_Stack_m, only: Stack_Depth

    integer, intent(in) :: ROOT    ! Root of expression subtree
    integer, intent(out) :: TYPE   ! Type of the expression value
    integer, intent(out) :: UNITS  ! Units of expression value if type is
                                   ! num_value, else type if type is enum_value
    double precision, intent(out) :: VALUE   ! Expression value, if any
    integer, intent(in), optional :: FIELD    ! What TO_CHECK should be
    integer, intent(in), optional :: START    ! of sons of FIELD
    integer, intent(out), optional :: Field_Look ! Last name of x.y, for error message
    integer, intent(out), optional :: Field_Test ! Allowed field name, for error message
    integer, intent(out), optional :: Tree    ! Decl%Tree if identifier

    ! Stat is zero unless TO_CHECK, FIELD, and START are present,
    ! node_id(ROOT) == n_dot, and the operands of n_dot are wrong for
    ! FIELD.

    integer :: Arg_Tree            ! Root of tree of allowed argument types
    integer :: Decor               ! Decoration of a tree node
    type(decls) :: DECL, DECL2     ! Declaration record for "root", function arg
    integer :: I, J
    integer :: ME                  ! node_id(root)
    integer :: SON1, SON2, SON3    ! Sons of "root"
    integer :: STRING              ! sub_rosa(root)
    integer :: Trace = -1          ! String index for trace
    integer :: Tree2
    integer :: TYPE2, TYPE3        ! Types for sons of "root"
    integer :: TypeString          ! String index of defined type, else zero
    integer :: UNITS2, UNITS3      ! Units for sons of "root"
    double precision :: VALUE2, VALUE3 ! Values for sons of "root"

    call trace_begin ( trace, 'EXPR', root, cond=toggle(con) )
    stat = 0 ! Assume status is OK
    typeString = 0                 ! Assume not a defined type
    units = phyq_dimensionless     ! default
    value = 0.0d0                  ! default
    if ( present(tree) ) tree = -1 ! default
    me = node_id(root)
    select case ( me )
    case ( n_identifier, n_number, n_string ) ! ----------------------------
      string = sub_rosa(root)
      if ( .not. declared(string) ) then
        select case ( me )
        case ( n_identifier )
                       ! String  Value  Type   Units         Tree
          call declare ( string, 0.0d0, empty, phyq_invalid, root )
        case ( n_number )
                       ! String  Value                Type
          call declare ( string, float_value(string), num_value, &
                       ! Units               Tree
                       & phyq_dimensionless, root )
        case ( n_string )
                       ! String  Value         Type       Units         Tree
          call declare ( string, 0.0d0+string, str_value, phyq_invalid, root )
        end select
      end if
      decl = get_decl(string, [do_label,enum_value,label,named_value,num_value, &
                            &  str_value,variable] )
      if ( decl%type == null_decl ) decl = declaration(string)
      type = decl%type
      if ( type == null_decl ) then
        type = empty
        go to 9
      end if
      select case ( type )
      case ( enum_value )
        value = decl%units ! Lit index
        units = decoration(subtree(1,decl%tree))
        typeString = data_type_indices(units)
      case ( label, named_value )
        value = decl%units ! Spec_index
      case ( variable )
        if ( decl%units == enum_value ) then
          type = decl%units
          units = decl%values(1)%type
          typeString = data_type_indices(units)
        else
          type = decl%values(1)%type
          units = decl%values(1)%units(1)
        end if
        value = decl%values(1)%value(1)
      case default
        units = decl%units
        value = decl%value
      end select
      if ( type == log_value ) then
        type = enum_value
        units = t_boolean
        typeString = units
      end if
      if ( present(tree) ) tree=decl%tree
    case ( n_unit ) ! ------------------------------------------------------
      son1 = subtree(1,root)
      stat = expr (son1, type, units, value, field, start, field_look, field_test)
      if ( stat /= 0 ) &
  go to 9
      son2 = subtree(2,root)
      decl = get_decl(sub_rosa(subtree(2,root)), units_name)
      if ( decl%type /= units_name ) then
        call local_error ( son2, not_units )
      else if ( units /= phyq_dimensionless .and. units /= decl%units .and. &
        & decl%units /= phyq_dimensionless ) then
        call local_error ( root, inconsistent_units, (/ son1, son2 /) )
      else
        units = decl%units
        if ( decl%value > 0 ) then
          value = value * decl%value
        else
          value = value - decl%value
        end if
      end if
    case ( n_cond ) ! ------------------------------------------------------
      son1 = subtree(1,root); son2 = subtree(2,root); son3 = subtree(3,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      if ( type /= enum_value .or. units /= t_boolean ) &
        & call local_error ( root, not_log_value )
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son3, type3, units3, value3, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      type = type2
      if ( type2 /= type3 ) &
        & call local_error ( son3, inconsistent_types, got=type3, expect=type2 )
    case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less ) ! ------
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      if ( type /= type2 ) then
        call local_error ( son2, inconsistent_types, (/ son1, son2 /) )
      end if
! Don't test units here.  The field might use : for something other than range.
!       if ( type == num_value .and. type2 == num_value .and. &
!            units /= phyq_dimensionless .and. &
!            units2 /= phyq_dimensionless .and. &
!            units /= units2 ) then
!         call local_error ( root, inconsistent_units, (/ son1, son2 /) )
!       end if
      if ( type == num_value ) type = range
      if ( type == str_value ) type = str_range
    case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, n_not_equal )
      ! Value is zero if expression is false, or if either operand is not numeric
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son2, type2, units2, value2 )
      if ( stat /= 0 ) &
  go to 9
      if ( type /= type2 ) then
        call local_error ( son2, inconsistent_types, (/ son1, son2 /) )
      end if
      type = enum_value
      units = t_boolean
      select case ( me )
      case ( n_equal_equal )
        value = merge ( l_true, l_false, value == value2 )
  go to 9
      case ( n_not_equal )
        value = merge ( l_true, l_false, value /= value2 )
  go to 9
      end select
      if ( type /= variable .and. type2 /= variable ) stat = dot_or_num ()
      if ( type == num_value .and. type2 == num_value ) then
        if ( units /= phyq_dimensionless .and. &
             units2 /= phyq_dimensionless .and. &
             units /= units2 ) then
          call local_error ( root, inconsistent_units, (/ son1, son2 /) )
        end if
        select case ( me )
        case ( n_less )
          value = merge ( l_true, l_false, value < value2 )
        case ( n_less_eq )
          value = merge ( l_true, l_false, value <= value2 )
        case ( n_greater )
          value = merge ( l_true, l_false, value > value2 )
        case ( n_greater_eq )
          value = merge ( l_true, l_false, value >= value2 )
        end select
      end if
    case ( n_and, n_or ) ! -------------------------------------------------
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      if ( type /= str_value .and. &
           & ( type /= enum_value .or. units /= t_boolean ) ) then
        call local_error ( son1, not_string )
      else if ( type2 /= str_value .and. &
           & ( type2 /= enum_value .or. units /= t_boolean ) ) then
        call local_error ( son2, not_string )
      else
        if ( me == n_and ) then
          value = merge(l_true, l_false, value==l_true .and. value2 == l_true )
        else ! node_id(root) == n_or
          value = merge(l_true, l_false, value==l_true .or. value2 == l_true )
        end if
      end if
      type = enum_value
      units = t_boolean
    case ( n_not )
      son1 = subtree(1,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      if ( type /= str_value .and. &
           & ( type /= enum_value .or. units /= t_boolean ) ) then
        call local_error ( son1, not_string )
      else
        value = merge(l_true, l_false, value==l_false)
      end if
      type = enum_value
      units = t_boolean
    case ( n_plus, n_minus ) ! ---------------------------------------------
      son1 = subtree(1,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      if ( nsons(root) > 1 ) then
        son2 = subtree(2,root)
        stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
        if ( stat /= 0 ) &
  go to 9
        stat = dot_or_num ()
        if ( type /= num_value ) &
  go to 9
        if ( type == num_value .and. type2 == num_value ) then
          if ( units /= phyq_dimensionless .and. &
               units2 /= phyq_dimensionless .and. &
               units /= units2 ) then
            call local_error ( root, inconsistent_units, (/ son1, son2 /) )
          end if
          if ( me == n_plus ) then
            value = value + value2
          else !  me == n_minus
            value = value - value2
          end if
        else if ( type == dot .or. type2 == dot ) then
          type = dot
        end if
      else if ( me == n_minus ) then
        value = - value
      end if
    case ( n_mult ) ! ------------------------------------------------------
      son1 = subtree(1,root) ; son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = dot_or_num ()
      if ( type /= num_value ) &
  go to 9
      if ( type == num_value .and. type2 == num_value .and. &
           units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call local_error ( root, inconsistent_units, (/ son1, son2 /) )
      else
        value = value * value2
        if ( units == phyq_dimensionless ) units = units2
      end if
    case ( n_div, n_into ) ! -----------------------------------------------
      son1 = subtree(1,root) ; son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  go to 9
      stat = dot_or_num ()
      if ( type /= num_value ) &
  go to 9
      if ( type == num_value .and. type2 == num_value .and. &
           units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call local_error ( root, inconsistent_units, (/ son1, son2 /) )
      else if ( me == n_div ) then
        value = value / value2
      else ! me == n_into
        value = value2 / value
        units = units2
      end if
    case ( n_pow ) ! -------------------------------------------------------
      value = 1.0
      do i = nsons(root), 1, -1 ! POW is right associative
        son1 = subtree(i,root)
        stat = expr ( son1, type, units, value2, field, start, field_look, field_test )
        if ( stat /= 0 ) &
  go to 9
        if ( type /= num_value ) then
          call local_error ( root, not_numeric, (/ son1 /) )
          stat = wrong_type
        else if ( units /= phyq_dimensionless ) then
          call local_error ( son1, not_unitless )
          stat = wrong_type
        else
          value = value2 ** value
        end if
      end do
    case ( n_func_ref ) ! --------------------------------------------------
      type = num_value ! Default, otherwise same as last arg if arg types checked
      son1 = subtree(1,root)
      ! Look up the function name
      string = sub_rosa(son1)
      decl = get_decl(string,function)
      do while ( decl%tree /= null_tree )
        if ( decl%type == function ) then
          if ( nsons(decl%tree) < 2 ) then ! Check only consistency of args
            decor = decoration(decl%tree)
            if ( decor > 0 ) then
              son1 = subtree(2,root)
              type = empty
              if ( nsons(root) > 1 ) then
                stat = expr ( son1, type, units, value, &
                            & field, start, field_look, field_test )
                if ( stat /= 0 ) &
  go to 9
              end if
              if ( mod(decor,u) /= 0 ) then ! Uniform argument types required
                do i = 3, nsons(root)
                  son2 = subtree(i,root)
                  stat = expr ( son2, type2, units2, value, &
                              & field, start, field_look, field_test )
                  if ( stat /= 0 ) &
  go to 9
                  if ( type /= empty .and. type /= type2 ) then
                    call local_error ( son2, inconsistent_types, &
                      & (/ son1, son2 /), expect=type, got=type2 )
                  else if ( type == enum_value ) then
                    decl2 = get_decl(sub_rosa(son2),enum_value)
                    type = decl2%tree
                    units = decl%units
                  else if ( type == num_value .and. units /= units2 ) then
                    call local_error ( root, inconsistent_units, (/ son1, son2 /) )
                  end if
                end do
              end if
              if ( decor >= u ) type = decor / u
            end if
          else ! check argument number and types
            arg_tree = subtree(2,decl%tree)
            if ( node_id(arg_tree) /= n_arg_def ) & ! No argument # or type checking
  go to 9
            if ( node_id(arg_tree) == n_arg_def ) then
              if ( nsons(root) /= 1 + nsons(arg_tree) ) then
                call local_error ( root, wrong_num_args, fields=(/string/) )
  go to 9
              end if
            end if
            if ( node_id(arg_tree) == n_arg_def ) then
              do i = 1, nsons(arg_tree)
                stat = expr ( subtree(i+1,root), type, units, value, &
                            & field, start, field_look, field_test )
                if ( stat /= 0 ) &
  go to 9
                son2 = subtree(i,arg_tree)
                if ( node_id(subtree(i,arg_tree)) == n_or ) then
                  do j = 1, nsons(son2)
                    if ( type == decoration(subtree(j,son2)) ) then
                      if ( type == num_value .and. units /= phyq_dimensionless ) &
                        & call local_error ( subtree(i+1,root), not_unitless, fields=(/string/) )
  go to 9
                    end if
                  end do
                else
                  if ( type /= decoration(son2) ) then
                    call local_error ( subtree(i+1,root), wrong_arg_type, &
                      & fields=(/string/), expect=subtree(i,arg_tree) )
                  else if ( units /= phyq_dimensionless ) then
                    call local_error ( subtree(i+1,root), not_unitless, fields=(/string/) )
                  end if
                end if
              end do
            end if
          end if
  go to 9
        end if
        decl = prior_decl(decl,function)
      end do
      call local_error ( son1, not_func, fields=(/string/) )
    case ( n_dot ) ! -------------------------------------------------------
      stat = check_dot ( root, field, start, field_look, field_test )
      type = empty
      if ( stat /= 3 ) type = dot
    case ( n_array )
      son1 = subtree(1,root)
      stat = expr ( son1, type, units, value, tree=tree )
      if ( stat /= 0 ) &
  go to 9
      do i = 2, nsons(root)
        son2 = subtree(i,root)
        stat = expr ( son2, type2, units2, value2, tree=tree2 )
        if ( stat /= 0 ) &
  go to 9
        if ( type /= type2 ) then
          call local_error ( son2, inconsistent_types, (/ son1, son2 /) )
        else if ( (type == enum_value .or. type == label .or. type == do_label) &
                  & .and. tree /= tree2 ) then
          call local_error ( son2, inconsistent_types, (/ son1, son2 /) )
        end if
      end do
    case default ! ---------------------------------------------------------
      call local_error ( root, no_code_for )
    end select
  9 if ( stat == 0 ) then
      select case ( me )
      case ( n_identifier, n_number, n_string ) ! ----------------------------
      case default
        call decorate ( root, type )
      end select
    end if
    if ( toggle(con) .and. levels(con) > 0 ) then
      do i = 0, stack_depth()
        call output ( '_' )
      end do
      call display_string ( data_type_indices(type) )
      select case ( type )
      case ( enum_value )
        call display_string ( data_type_indices(units), before=' ' )
        call display_string ( lit_indices(nint(value)), before=' ' )
      case ( label )
        call display_string ( spec_indices(nint(value)), before=' ' )
      case ( num_value, range )
        call output ( value, before=' ' )
      case ( str_range, str_value )
        call display_string ( nint(value), before=' ' )
      end select
      call newLine
    end if
    if ( typeString == 0 ) typeString = data_type_indices(type)
    call trace_end ( 'EXPR', index=stat, & ! string=trim(type_names(type)), 
      & stringIndex=typeString, cond=toggle(con) )

  contains

    integer function Dot_or_Num ( )
      ! If type and type2 are both dot, or if one is dot and the other
      ! is numeric the type is dot.  Otherwise if both are numeric
      ! it is num_value.  Otherwise it is empty, and an error has been
      ! announced.  If the type is empty, the result is wrong_type, else
      ! it is zero.
      dot_or_num = wrong_type ! Assume wrong
      if ( type == type2 ) then
        if ( type == dot .or. type == num_value ) then
          dot_or_num = 0
        end if
      else if ( type == dot .and. type2 == num_value .or. &
                type == num_value .and. type2 == dot  ) then
        type = dot
        dot_or_num = 0
      end if
      if ( dot_or_num /= 0 ) call local_error ( root, wrong_type )
    end function Dot_or_Num

    subroutine LOCAL_ERROR ( WHERE, CODE, SONS, FIELDS, EXPECT, GOT )
      integer, intent(in) :: WHERE   ! Tree node where error was noticed
      integer, intent(in) :: CODE    ! Code for error message
      integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                     ! "where".  If they're pseudo_terminal,
                                     ! their declarations are dumped.
      integer, intent(in), optional :: FIELDS(:) ! Field indices
      integer, intent(in), optional :: EXPECT    ! Expected type
      integer, intent(in), optional :: GOT       ! Actual type
      call announce_error ( where, code, sons, fields, expect, got )
      type = empty
      units = phyq_invalid
    end subroutine LOCAL_ERROR

  end function EXPR

! -------------------------------------------------  IF_Construct  -----
  recursive subroutine IF_Construct ( Root, Son1, Decl )
    ! < n_if <n_test expr spec*>* <n_else spec*>? >
    integer, intent(in) :: Root ! N_If node
    integer, intent(in), optional :: Son1 ! First son of parent CF
    type(decls), intent(in), optional :: DECL ! Declaration of "name" on "begin name"
    integer :: Begin, Gson, I, J
    integer :: Me = -1             ! String index for trace
    integer :: Son, Start, Stat, String, Type, Units
    double precision :: Value      ! not used

    call trace_begin ( me, 'IF_Construct', root, cond=toggle(con) )
    begin = 1 ! Assume first son is not n_named
    son = subtree(1,root)
    if ( node_id(son) == n_named ) then
      gson = subtree(1,son)
      if ( check_label(gson) ) then
        string = sub_rosa(gson)
        call declare ( string, 0.0d0+string, do_label, phyq_invalid, root )
      end if
      begin = 2
    end if
    do i = begin, nsons(root)
      son = subtree(i,root)
      start = begin   ! Assume ELSE
      if ( node_id(son) == n_test ) then ! if ( expr ) or else if ( expr )
        stat = expr(subtree(1,son),type,units,value)
        if ( type /= enum_value .or. units /= t_boolean ) call announce_error ( &
          & subtree(1,son), wrong_expr_type, expect=log_value )
        start = begin+1 ! ELSE IF, skip the EXPR we just processed
      end if
      if ( present(son1) ) then
        do j = start, nsons(son)
          gson = subtree(j,son)
          if ( node_id(gson) == n_cf ) then
            call one_section ( gson )
          else
            call one_cf ( gson, son1, decl )
          end if
        end do
      else
        do j = start, nsons(son)
          call one_section ( subtree(j,son) )
        end do
      end if
    end do
    call trace_end ( 'IF_Construct', cond=toggle(con) )

  end subroutine IF_Construct

! -------------------------------------------------------  ONE_CF  -----
  recursive subroutine ONE_CF ( ROOT, SON1, Decl )
  ! Analyze one configuration, with abstract syntax tree rooted at ROOT.
    integer, intent(in) :: ROOT
    integer, intent(in) :: SON1    ! First son of CF
    type(decls), intent(in) :: DECL ! Declaration of "name" on "begin name"

    integer :: GSON1, GSON2        ! Grandsons of ROOT being analyzed
    integer :: Me = -1             ! String index for trace
    integer :: STRING1             ! Sub_Rosa

    call trace_begin ( me, 'ONE_CF', root, cond=toggle(con) )
    select case ( node_id(root) )
    case ( n_cycle, n_exit)
      call cycle_exit ( root, -1 )  ! -1 => Inside section
    case ( n_do )                   ! DO construct
      call do_construct ( root, son1, decl )
    case ( n_equal )                ! x = expr
      call decorate ( root, son1 )  ! show equal the section name
      if ( decoration(decl%tree) /= no_check_eq ) call equal ( root )
    case ( n_if )                   ! IF construct
      call if_construct ( root, son1, decl )
    case ( n_named )                ! label:x,(y[=expr])+
      gson1 = subtree(1,root)       ! Label
      if ( check_label(gson1) ) then
        gson2 = subtree(2,root)       ! spec_args tree or DO construct tree
        call decorate ( gson1, gson2 )
        string1 = sub_rosa(gson1)
        call decorate ( gson2, son1 ) ! show spec_args the section name
        call spec_args ( gson2 )
        call declare ( string1, 0.0d0+string1, label, &
                     ! units = spec index
          &            decoration(subtree(1,decoration(subtree(1,gson2)))), &
                     ! tree
          &            gson2 )
      end if
    case ( n_select )               ! SELECT construct
      call select_construct ( root, son1, decl )
    case ( n_spec_args )            ! x,(y[=expr])+
      call decorate ( root, son1 )  ! show spec_args the section name
      call spec_args ( root )
    case ( n_variable )
      call variable_def ( root )
    case ( n_while )                ! WHILE construct
      call while_construct ( root, son1, decl )
    case default
      call announce_error ( root, no_code_for )
    end select
    call trace_end ( 'ONE_CF', cond=toggle(con) )
  end subroutine ONE_CF

! --------------------------------------------------  ONE_SECTION  -----
  recursive subroutine ONE_SECTION ( ROOT )
  ! Analyze one section, with abstract syntax tree rooted at ROOT.
  ! Check that the names at BEGIN and END are the same.
    use Tree, Only: Where
    integer, intent(in) :: ROOT

    type(decls) :: DECL            ! Declaration of "name" on "begin name"
    integer :: I, N                ! Loop inductor, number of sons
    integer :: Me = -1             ! String index for trace
    integer :: SON1, SONN          ! Sons of ROOT being analyzed
    integer :: STRING1, STRINGN    ! Sub_Rosa of first and last sons

    call trace_begin ( me, 'ONE_SECTION', root, cond=toggle(con) )
    num_sections = num_sections + 1
    n = nsons(root)
    son1 = subtree(1,root) ;   sonn = subtree(n,root)
    string1 = sub_rosa(son1) ; stringn = sub_rosa(sonn)
    if ( string1 /= stringn ) then
      call output ( '***** Names on BEGIN at ' )
      call print_source ( where(subtree(1,root)), advance='yes' )
      call output ( '*****      and END at ' )
      call print_source ( where(subtree(n,root)) )
      call output ( ' are not the same.', advance='yes' )
!     call output ( '***** Processing suppressed.', advance='yes' )
      error = max(error,1)
    end if
    decl = get_decl(sub_rosa(son1),section)
    if ( decl%type /= section ) &
      call announce_error ( son1, not_section )
    if ( section_ordering(decl%units,current_section) /= 1 ) &
      call announce_error ( son1, section_order )
    current_section = decl%units
    if ( decl%tree /= null_tree ) &
      & call decorate ( son1, decoration(subtree(1,decl%tree)) )
    do i = 2, n-1
      call one_cf ( subtree(i,root), son1, decl )
    end do
    call trace_end ( 'ONE_SECTION', cond=toggle(con) )
  end subroutine ONE_SECTION

! ------------------------------------------------  Push_Do_Stack  -----
  subroutine Push_Do_Stack ( Root )
  ! Allocate or reallocate the do construct stack if necessary.
  ! Push Root on the stack.
    use Allocate_Deallocate, only: Bytes, Test_Allocate
    integer, intent(in) :: Root
    integer :: Me = -1             ! String index for trace
    integer :: N, Stat             ! Size of old stack, status from (de)allocate
    integer, allocatable :: Temp_Stack(:)

    call trace_begin ( me, 'Push_Do_Stack', root, cond=toggle(con) )
    if ( .not. allocated(do_construct_stack) ) then
      allocate ( do_construct_stack(10), stat=stat )
      call test_allocate ( stat, moduleName, 'Do_Construct_Stack', ubounds=[10], &
        & elementSize=bytes(do_construct_stack) )
    end if
    do_stack_top = do_stack_top + 1
    if ( do_stack_top > ubound(do_construct_stack,1) ) then
      n = size(do_construct_stack)
      allocate ( temp_stack(2*n), stat=stat )
      call test_allocate ( stat, moduleName, 'Temp_Stack', ubounds=[10], &
        & elementSize=bytes(do_construct_stack) )
      temp_stack(:n) = do_construct_stack
      call move_alloc ( temp_stack, do_construct_stack )
    end if
    do_construct_stack(do_stack_top) = root
    call trace_end ( 'Push_Do_Stack', cond=toggle(con) )

  end subroutine Push_Do_Stack

! ---------------------------------------------  Select_Construct  -----
  subroutine Select_Construct ( Root, Son1, Decl )
    ! < n_select expr <n_case expr spec*>* <n_default spec*>? >
    integer, intent(in) :: Root ! N_Select node
    integer, intent(in), optional :: Son1 ! First son of parent CF
    type(decls), intent(in), optional :: DECL ! Declaration of "name" on "begin name"
    integer :: Begin, Gson, I, J
    integer :: Me = -1        ! String index for trace
    integer :: Son, Start, Stat, String, Type1, Type2, Units
    double precision :: Value ! not used
    call trace_begin ( me, 'Select_Construct', root, cond=toggle(con) )
    begin = 1 ! Assume first son is not n_named
    son = subtree(1,root)
    if ( node_id(son) == n_named ) then
      gson = subtree(1,son)
      if ( check_label(gson) ) then
        string = sub_rosa(gson)
        call declare ( string, 0.0d0+string, do_label, phyq_invalid, root )
      end if
      begin = 2
    end if
    stat = expr(subtree(begin,root),type1,units,value)
    do i = begin+1, nsons(root)
      son = subtree(i,root)
      start = begin    ! Assume CASE DEFAULT
      if ( node_id(son) == n_test ) then
        stat = expr(subtree(1,son),type2,units,value)
        if ( type1 /= type2 ) call announce_error ( &
          & subtree(1,son), wrong_expr_type, expect=type1 )
        start = 2  ! CASE ( expr ), skip the expr we just processed 
      end if
      if ( present(son1) ) then
        do j = start, nsons(son)
          gson = subtree(j,son)
          if ( node_id(gson) == n_cf ) then
            call one_section ( gson )
          else
            call one_cf ( gson, son1, decl )
          end if
        end do
      else
        do j = start, nsons(son)
          call one_section ( subtree(j,son) )
        end do
      end if
    end do
    call trace_end ( 'Select_Construct', cond=toggle(con) )
  end subroutine Select_Construct

! ------------------------------------------------------  SET_ONE  -----
  subroutine SET_ONE ( ROOT )
  ! Analyze a field of a specification of the form "/name" --
  ! make sure "name" is a field of the specification.
    integer, intent(in) :: ROOT    ! Index of the "n_set_one" tree node
    integer :: FIELD               ! Tree node of field's declaration
    integer :: FIELD_LIT           ! f_... for a field
    integer :: Me = -1             ! String index for trace
    integer :: SON                 ! Son of n_set_one tree node
    integer :: SPEC_DECL           ! Tree node of the spec's declaration

    call trace_begin ( me, 'SET_ONE', root, cond=toggle(con) )
    son = subtree(1,root)
    spec_decl = decoration(root)
    field = check_field(son,spec_decl)
    if ( field == 0 ) then
      call announce_error ( son, not_field_of, fields=(/ spec_decl /) )
      call dump_1_decl ( sub_rosa(son) )
    else
      if ( check_field_type(field,t_boolean) == 0 ) then
        call decorate ( root, field )
        field_lit = decoration(subtree(1,field))
        call decorate ( son, field_lit )
        if ( no_dup_flag ) then
          if ( got(field_lit) ) &
            & call announce_error ( root, no_duplicate_fields, &
            & fields=(/ field_lit /) )
        end if
        got(field_lit) = .true.
      else
        call announce_error ( son, wrong_type )
      end if
    end if
    call trace_end ( 'SET_ONE', cond=toggle(con) )
  end subroutine SET_ONE

! ----------------------------------------------------  SPEC_ARGS  -----
  subroutine SPEC_ARGS ( ROOT )
  ! Analyze a specification of the form NAME (, EXPR ( = EXPR)? )+
  ! starting at ROOT
    integer, intent(in) :: ROOT

    integer :: FIELD_LIT      ! f_... for a field
    integer :: FLAGS          ! Flags from decoration(spec_decl)
    integer :: I              ! Loop inductor
    integer :: Me = -1        ! String index for trace
    integer :: SECT           ! The section the spec appears in
    integer :: SON            ! I'th son of "root"
    type(decls) :: SPEC_DECL  ! Declaration of first son, if any
    integer :: STAT           ! from "expr"
    integer :: TYPE, UNITS    ! Output from "expr"
    double precision :: VALUE ! Output from "expr"

    call trace_begin ( me, 'SPEC_ARGS', root, cond=toggle(con) )
    sect = decoration(root)
    spec_decl = get_decl( sub_rosa(sect), section )   ! Known to be a
                             ! valid section name if we get this far
    son = subtree(1,root)
    ! Is spec allowed in section?
    if ( check_section(son,spec_decl%tree) == 0 ) then
      call announce_error ( son, out_of_place, (/ sect /) )
    else
      spec_decl = get_decl( sub_rosa(son), spec ) ! Decl of son of type "spec"
      if ( spec_decl%type == empty ) then
        call announce_error ( son, not_spec )
      else
        call decorate ( son, spec_decl%tree )
        flags = decoration(spec_decl%tree)
        all_fields_flag = mod(flags/all_fields,2) /= 0
        no_dup_flag = mod(flags/no_dup,2) /= 0
        got = .false.
        do i = 2, nsons(root)
          son = subtree(i,root)
          select case ( node_id(son) )
          case ( n_asg )
            call decorate ( son, spec_decl%tree )
            call assign ( son, type, units, value )
          case ( n_set_one )
            call decorate ( son, spec_decl%tree )
            call set_one ( son )
          case ( n_dot )
          ! ??? Need some checking and generation here
          case default
            if ( mod(flags/no_positional,2) /= 0 ) then
              call announce_error ( son, no_positional_fields )
            else
              stat = expr ( son, type, units, value )
              ! ??? Generate some table here
            end if
          end select
        end do
        do i = 2, nsons(spec_decl%tree)
          son = subtree(i,spec_decl%tree)
          field_lit = decoration(subtree(1,son))
          if ( .not. got(field_lit) ) then
            if ( all_fields_flag .or. &
              & mod(decoration(son)/req_fld,2) /= 0 .and. decoration(son) < u ) &
              & call announce_error ( root, missing_field, &
                & fields=(/ field_lit /) )
          end if
        end do
      end if
    end if
    call trace_end ( 'SPEC_ARGS', cond=toggle(con) )
  end subroutine SPEC_ARGS

! -------------------------------------------------  Variable_Def  -----
  subroutine Variable_Def ( Root, Type )
    ! Verify that the expressions in the right-hand side of := all have the
    ! same type.
    ! If the first son of Root has previously appeared on the left side of :=
    ! check that the type of the left-hand side is the same as the right.
    ! If the first son of Root has not previously appeared, declare it with
    ! the type of the right-hand side.

    use Declaration_Table, only: Allocate_Test, Empty
    integer, intent(in) :: Root ! Tree index of n_variable from := operator
    integer, intent(out), optional :: Type
    type(decls) :: Decl
    integer :: Me = -1          ! String index for trace cacheing
    integer :: Son1             ! Left son
    integer :: String           ! String index of variable at Son1
    integer :: Units1           ! Of right sons
    integer :: Type1            ! Types of right sons
    double precision :: Value   ! Lit index of enumerator, or...
    type(value_t), allocatable :: Values(:)
    integer :: What1            ! What is RHS of :=

    call trace_begin ( me, 'Variable_Def', root, cond=toggle(con) )
    son1 = subtree(1,root)
    string = sub_rosa(son1)
    if ( nsons(root) < 2 ) then ! "variable :=" case
      if ( present(type) ) type = empty
                   ! String  Value         Type      Units  Tree
      call declare ( string, string+0.0d0, variable, empty, son1, values )
      call trace_end ( 'Variable_Def', stringIndex=string, cond=toggle(con) )
      return
    end if
    type1 = -1
    what1 = -1
    units1 = -1
    call check_type ( root, 2, what1, units1, type1, value )
    decl = get_decl ( string, [ enum_value, label, named_value, variable ] )
    select case ( decl%type )
    case ( enum_value, label, named_value )
      call announce_error ( son1, variable_conflict, got=string )
    case ( variable )
      select case ( what1 )
      case ( enum_value )
        if ( decl%values(1)%type /= units1 ) then
          call announce_error ( root, inconsistent_data_types, &
            & expect=decl%values(1)%type, got=type1 )
        end if
      case ( label )
        if ( decl%values(1)%type /= type1 ) then
          call announce_error ( root, inconsistent_data_types, &
            & expect=decl%values(1)%type, got=type1 )
        end if
      case default
        if ( decl%units /= what1 ) then
          call announce_error ( root, inconsistent_types, expect=decl%units, &
            & got=what1 )
        end if
      end select
    case default
      call allocate_test ( values, 1, 'Values', moduleName )
      select case ( type1 )
      case ( enum_value )
                     !   what   type   value  units
        values = value_t(type1, units1,value, 0 )
      case ( str_value )
                     !   what   type   value  units
        values = value_t(what1, type1, value, units1 )
      case default
                     !   what   type   value  units
        values = value_t(what1, type1, value, units1 )
      end select
                   ! String  Value         Type      Units  Tree
      call declare ( string, string+0.0d0, variable, what1, son1, values )
      decl = get_decl ( string, variable )
    end select

    select case ( what1 )
    case ( enum_value, label )
      call decorate ( son1, type1 ) ! Type def tree of enums
    case default
      call decorate ( son1, subtree(2,root) )
    end select

    if ( present(type) ) type = type1

    if ( decl%type /= variable ) then
      call trace_end ( 'Variable_Def', cond=toggle(con) )
    else
      call trace_end ( 'Variable_Def', &
        & stringIndex=data_type_indices(decl%units), cond=toggle(con) )
    end if

  contains

    recursive subroutine Check_Type ( Root, Start, What, Units, Type, Value )
      integer, intent(in) :: Root   ! of subtree to check
      integer, intent(in) :: Start  ! First son of Root to check
      integer, intent(inout) :: What  ! -1 initially, then type if type /=
                                      ! enum_value, data type index (t_...)
                                      ! if type == enum_value, or spec index
                                      ! (s_...) if type == label.  The sons
                                      ! must agree.
      integer, intent(inout) :: Units ! the sons must agree.
      integer, intent(inout) :: Type ! -1 if not yet set, the sons need to agree
      double precision, intent(out) :: Value

      type(decls) :: Decl ! of Sonn
      integer :: I
      integer :: Me = -1  ! String index for tracing
      logical :: Numeric
      integer :: Sonn     ! of Root
      integer :: Tree     ! Decl%Tree if sonn is identifier
      integer :: Type2    ! of Sonn
      integer :: Units2   ! of Sonn
      integer :: Unused
      integer :: What2    ! of Sonn

      call trace_begin ( me, 'Check_Type', root, cond=toggle(con) )
      do i = start, nsons(root)
        sonn = subtree(i,root)
        if ( node_id(sonn) == n_array ) then
          call check_type ( sonn, 1, what, units, type, value )
        else
          unused = expr ( sonn, type2, units2, value, tree=tree )
          numeric = type2 == num_value
          what2 = type2
          if ( type < 0 ) then
            units = units2
            what = what2
            type = type2
          end if
          if ( node_id(sonn) == n_identifier ) then
            select case ( type2 )
            case ( enum_value, label )
              decl = get_decl ( sub_rosa(sonn), type2 )
              call decorate ( sonn, decl%units )
            end select
          end if
          if ( type /= type2 .or. what /= what2 ) then
            call announce_error ( sonn, inconsistent_data_types, expect=type, &
              & got=type2 )
          end if
          if ( numeric .and. units /= units2 ) then
            call announce_error ( sonn, inconsistent_units, expect=units, &
              & got=units2 )
          end if
        end if
      end do
      call trace_end ( 'Check_Type', cond=toggle(con) )

    end subroutine Check_Type

  end subroutine Variable_Def

! ----------------------------------------------  While_Construct  -----
  recursive subroutine While_Construct ( ROOT, Son1, Decl )
  ! Analyze a WHILE construct of the form DO WHILE ( expr ) cfs END DO,
  ! starting at ROOT

    integer, intent(in) :: ROOT
    integer, intent(in), optional :: Son1 ! First son of parent CF
    type(decls), intent(in), optional :: DECL ! Declaration of "name" on "Begin"
    integer :: Gson, I, Son, Start, String
    integer :: Me = -1             ! String index for trace
    integer :: Stat, Type, Units
    double precision :: Value      ! not used

    call trace_begin ( me, 'WHILE_Construct', root, cond=toggle(con) )

    call push_do_stack ( root )

    son = subtree(1,root)
    start = 2 ! Assume first son is not n_named
    if ( node_id(son) == n_named ) then
      gson = subtree(1,son)
      if ( check_label(gson) ) then
        string = sub_rosa(gson)
        call declare ( string, 0.0d0+string, do_label, phyq_invalid, root )
      end if
      start = 3
      son = subtree(2,root)
    end if

    stat = expr(son,type,units,value)
    if ( type /= enum_value .or. units /= t_boolean ) call announce_error ( &
      & son, wrong_expr_type, expect=log_value )
    if ( present(son1) ) then
      call decorate ( root, -1 ) ! -1 => inside section
      do i = start, nsons(root)
        gson = subtree(i,root)
        if ( node_id(gson) == n_cf ) then
          call one_section ( gson )
        else
          call one_cf ( gson, son1, decl )
        end if
      end do
    else
      call decorate ( root, 1 ) ! +1 => enclosing sections
      do i = start, nsons(root)
        call one_section ( subtree(i,root) )
      end do
    end if

    do_stack_top = do_stack_top - 1 ! Pop the DO construct stack

    call trace_end ( 'WHILE_Construct', cond=toggle(con) )

  end subroutine While_Construct

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module TREE_CHECKER

! $Log$
! Revision 1.53  2016/11/02 22:18:57  vsnyder
! Don't use decl%tree as a tree node index if it isn't one
!
! Revision 1.52  2016/05/25 00:17:29  vsnyder
! Decruftication
!
! Revision 1.51  2015/02/05 21:51:26  vsnyder
! Work related to variable definitions
!
! Revision 1.50  2014/04/09 00:41:26  vsnyder
! Repair problem finding variables
!
! Revision 1.49  2014/03/20 18:36:11  vsnyder
! Unified type system instead of one in Intrinsic and one in Declaration_Table
!
! Revision 1.48  2014/02/27 02:37:18  vsnyder
! EXIT referring to IF and SELECT CASE constructs
!
! Revision 1.47  2014/02/21 19:19:01  vsnyder
! More work on variables, especially those with enumerator values
!
! Revision 1.46  2014/01/11 01:42:25  vsnyder
! Stuff for variables, some decruftification
!
! Revision 1.45  2013/12/12 01:54:53  vsnyder
! Provide variable assignment, and IF, and SELECT constructs
!
! Revision 1.44  2013/10/16 01:06:45  vsnyder
! Correct typo in first call to Declare in Def_Spec
!
! Revision 1.43  2013/10/12 01:19:04  vsnyder
! Handle undefined vector correctly in Check_Dot
!
! Revision 1.42  2013/10/11 00:46:16  vsnyder
! Variables and functions
!
! Revision 1.41  2013/09/24 23:28:00  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 1.40  2013/09/20 00:56:18  vsnyder
! Allow relational operators to have dot operands
!
! Revision 1.39  2013/09/19 23:34:33  vsnyder
! Type-check and decorate dots in expressions
!
! Revision 1.38  2013/08/31 01:26:06  vsnyder
! Don't pass root as index to trace_end
!
! Revision 1.37  2013/08/30 03:53:10  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 1.36  2013/08/17 02:54:54  vsnyder
! Remove references to DEPTH from trace_m
!
! Revision 1.35  2012/05/24 21:05:49  vsnyder
! Allow check_dot to have reflexive elements in the required spec tree, i.e.
! <dot a b ... z> really means <dot a b+ ... z>.  The first and last ones
! are not reflexive.  This allows, for example, a vectorTemplate to get
! its quantities from other vector templates.
!
! Revision 1.34  2012/05/05 00:12:36  vsnyder
! Add support for 'not' operator
!
! Revision 1.33  2011/08/20 00:48:34  vsnyder
! Remove unused use names and variable declarations
!
! Revision 1.32  2011/04/20 17:32:28  vsnyder
! Undo ill-advised units check, correct some error messages
!
! Revision 1.31  2011/04/19 02:00:31  vsnyder
! Support == and /= relational operators too
!
! Revision 1.30  2011/04/19 00:58:16  vsnyder
! Allow several types for a field
!
! Revision 1.29  2011/01/29 00:46:15  vsnyder
! Add units checking
!
! Revision 1.28  2006/10/10 23:49:40  vsnyder
! Repair bugs that resulted in not verifying that a field is of the form x.y
! when it should be, or that it is not of the form x.y when it shouldn't be.
!
! Revision 1.27  2006/03/23 01:50:43  vsnyder
! Check for empty fields
!
! Revision 1.26  2005/06/22 20:03:55  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.25  2004/12/31 02:11:34  vsnyder
! Simplify constructing some error messages
!
! Revision 1.24  2004/11/17 20:23:30  vsnyder
! Add checking for scalar field
!
! Revision 1.23  2004/06/23 02:12:24  vsnyder
! Add Check_Type
!
! Revision 1.22  2004/05/29 02:42:16  vsnyder
! Rearrange function definition stuff
!
! Revision 1.21  2004/05/29 00:54:46  vsnyder
! Improve some error checking and handling
!
! Revision 1.20  2004/05/28 23:44:28  vsnyder
! Get units from either operand of *, second operand of \\
!
! Revision 1.19  2004/05/28 23:15:15  vsnyder
! Add power (^) operator, functions
!
! Revision 1.18  2004/02/14 00:15:18  vsnyder
! More precise error message for wrong type
!
! Revision 1.17  2004/01/14 02:19:51  vsnyder
! Get PHYQ_INVALID from Intrinsic instead of Init_Tables_Module
!
! Revision 1.16  2003/08/29 00:14:43  vsnyder
! Correct out-of-bounds subscript
!
! Revision 1.15  2002/10/02 00:43:41  vsnyder
! Remove check for consistent units
!
! Revision 1.14  2002/07/25 00:25:10  vsnyder
! In Check_Field_Type, say NO if type is out-of-range -- i.e., don't crash
!
! Revision 1.13  2002/05/23 20:35:57  vsnyder
! Eliminate two unused variables
!
! Revision 1.12  2001/11/28 23:48:30  vsnyder
! Correct blunders in arrays-of-arrays
!
! Revision 1.11  2001/11/28 03:15:37  vsnyder
! Implement arrays of arrays
!
! Revision 1.10  2001/11/27 00:50:45  vsnyder
! Implement (partially) open ranges
!
! Revision 1.9  2001/06/07 21:56:55  pwagner
! Added Copyright statement
!
! Revision 1.8  2001/05/18 21:24:26  vsnyder
! Missing 'sub_rosa' around a 'son' in a call to 'dump_1_decl'
!
! Revision 1.7  2001/03/09 22:08:27  vsnyder
! Improve error message for wrong type of field value
!
! Revision 1.6  2001/03/06 22:50:11  vsnyder
! Correct processing of /foo when it's not a valid field.
!
! Revision 1.5  2001/03/05 23:20:09  vsnyder
! Correct obscure problem that only occurs if you have erroneous input
!
! Revision 1.4  2001/02/22 19:43:04  vsnyder
! Improve some messages
!
! Revision 1.3  2001/02/07 19:42:06  vsnyder
! Add checking for duplicate fields, all fields and no-positional.
!
! Revision 1.2  2001/02/02 00:04:36  vsnyder
! Improved some error messages
!
! Revision 1.1  2000/11/02 21:36:45  vsnyder
! Initial entry into CVS
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.3  2000/07/28 18:31:04  vsnyder
! Added support for checking for "x.y" fields.
!
