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
! cross-references within the tree with anything that might be useful
! in table_generator.

  use DECLARATION_TABLE, only: DECLARATION, DECLARE, DECLARED, DECLS, &
    &                          DUMP_1_DECL, EMPTY, ENUM_VALUE, FIELD, &
    &                          FUNCTION, GET_DECL, LABEL, LOG_VALUE, &
    &                          NAMED_VALUE, NULL_DECL, NUM_VALUE, PRIOR_DECL, &
    &                          RANGE, REDECLARE, SECTION, SPEC, STR_VALUE, &
    &                          UNDECLARED, TYPE_MAP, TYPE_NAME, UNITS_NAME
  use INIT_TABLES_MODULE, only: DATA_TYPE_INDICES, FIELD_FIRST, FIELD_INDICES, &
    &                           FIELD_LAST, LIT_INDICES, PHYQ_DIMENSIONLESS, &
    &                           SECTION_FIRST,SECTION_INDICES, SECTION_LAST, &
    &                           SECTION_ORDERING, T_BOOLEAN
  use INTRINSIC, only: ALL_FIELDS, EMPTY_OK, NO_ARRAY, NO_CHECK_EQ, NO_DUP, &
    &                  NO_POSITIONAL, PHYQ_INVALID, REQ_FLD, U
  use LEXER_CORE, only: PRINT_SOURCE
  use MoreTree, only: Scalar, StartErrorMessage
  use OUTPUT_M, only: NEWLINE, OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING, FLOAT_VALUE
  use TOGGLES, only: CON, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, &
                  NODE_ID, NODE_KIND, NSONS, NULL_TREE, PSEUDO, SOURCE_REF,  &
                  SUB_ROSA, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  implicit NONE
  private

  public :: CHECK_TREE, CHECK_TYPE

! -----     Private declarations     -----------------------------------
  integer, private :: ERROR   ! 0 => No errors

! Error codes for "announce_error"
  integer, private, parameter :: ALREADY_DECLARED = 1
  integer, private, parameter :: ARRAY_NOT_ALLOWED = ALREADY_DECLARED + 1
  integer, private, parameter :: EMPTY_NOT_ALLOWED = ARRAY_NOT_ALLOWED + 1
  integer, private, parameter :: INCONSISTENT_TYPES = EMPTY_NOT_ALLOWED + 1
  integer, private, parameter :: INCONSISTENT_UNITS = INCONSISTENT_TYPES + 1
  integer, private, parameter :: MISSING_FIELD = INCONSISTENT_UNITS + 1
  integer, private, parameter :: NO_CODE_FOR = MISSING_FIELD + 1
  integer, private, parameter :: NO_DOT = NO_CODE_FOR + 1
  integer, private, parameter :: NO_DUPLICATE_FIELDS = NO_DOT + 1
  integer, private, parameter :: NO_POSITIONAL_FIELDS = NO_DUPLICATE_FIELDS + 1
  integer, private, parameter :: NO_SUCH_FIELD = NO_POSITIONAL_FIELDS + 1
  integer, private, parameter :: NO_SUCH_REFERENCE = NO_SUCH_FIELD + 1
  integer, private, parameter :: NOT_FIELD_OF = NO_SUCH_REFERENCE + 1
  integer, private, parameter :: NOT_FUNC = NOT_FIELD_OF + 1
  integer, private, parameter :: NOT_LIT_OF_TYPE = NOT_FUNC + 1
  integer, private, parameter :: NOT_NAME = NOT_LIT_OF_TYPE + 1
  integer, private, parameter :: NOT_NAME_OR_STRING = NOT_NAME + 1
  integer, private, parameter :: NOT_NUMERIC = NOT_NAME_OR_STRING + 1
  integer, private, parameter :: NOT_SECTION = NOT_NUMERIC + 1
  integer, private, parameter :: NOT_SPEC = NOT_SECTION + 1
  integer, private, parameter :: NOT_STRING = NOT_SPEC + 1
  integer, private, parameter :: NOT_UNITLESS = NOT_STRING + 1
  integer, private, parameter :: NOT_UNITS = NOT_UNITLESS + 1
  integer, private, parameter :: OUT_OF_PLACE = NOT_UNITS + 1
  integer, private, parameter :: SECTION_ORDER = OUT_OF_PLACE + 1
  integer, private, parameter :: WRONG_ARG_TYPE = SECTION_ORDER + 1
  integer, private, parameter :: WRONG_NUM_ARGS = WRONG_ARG_TYPE + 1
  integer, private, parameter :: WRONG_TYPE = WRONG_NUM_ARGS + 1
  integer, private, parameter :: WRONG_UNITS = WRONG_TYPE + 1

  logical, private :: ALL_FIELDS_FLAG   ! All fields are required
  integer, private :: CURRENT_SECTION = section_first - 1
  logical, private :: GOT(field_first:field_last)
  logical, private :: NO_DUP_FLAG       ! Duplicate named fields prohibited

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

! ---------------------------------------------------  CHECK_TREE  -----
  subroutine CHECK_TREE ( ROOT, ERROR_FLAG, FIRST_SECTION, HOW_MANY_SECTIONS )
  ! Traverse the abstract syntax tree starting at ROOT, which should
  ! be a N_CFS node (but we don't check).
    integer, intent(in) :: ROOT               ! Root of the tree
    integer, intent(out) :: ERROR_FLAG        ! /=0 means trouble
    integer, intent(out), optional :: FIRST_SECTION     ! Which son of root?
    integer, intent(out), optional :: HOW_MANY_SECTIONS ! Number of begin ... ends
    integer :: I              ! Loop inductor
    integer :: Me = -1        ! String index for trace
    integer :: NUM_SECTIONS   ! Number of begin ... ends
    integer :: SON            ! Son of root
    error = 0
    num_sections = 0
    if ( present(first_section) ) first_section = 0
    call trace_begin ( me, 'CHECK_TREE', root, cond=toggle(con) )
    do i = 1, nsons(root)
      son = subtree(i,root)
      select case ( node_id(son) )
      case ( n_cf )         ! A begin-end block's subtree
        if ( present(first_section) ) then
          if ( first_section == 0 ) first_section = i
        end if
        call one_cf ( son )
        num_sections = num_sections + 1
      case ( n_dt_def );    call def_type ( son ) ! Declare lits in the type
      case ( n_func_def );  call def_func ( son ) ! Declare function
      case ( n_section );   call def_section ( son )
      case ( n_spec_def );  call def_spec ( son )
      case default ;        call announce_error ( son, no_code_for )
      end select
    end do
    error_flag = error
    if ( present(how_many_sections) ) how_many_sections = num_sections
    call trace_end ( 'CHECK_TREE', cond=toggle(con) )
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
      if ( lit_decl%tree == type_decl%tree ) return
      if ( lit_decl%prior == null_decl ) exit
      lit_decl = prior_decl(lit_decl, enum_value)
    end do
    check_type = .false.
  end function Check_Type

! =====     Private Procedures     =====================================
! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, SONS, FIELDS, EXPECT )
    use Intrinsic, only: PHYQ_Indices
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                   ! "where".  If they're pseudo_terminal,
                                   ! their declarations are dumped.
    integer, intent(in), optional :: FIELDS(:) ! Field indices
    integer, intent(in), optional :: EXPECT    ! Something expected
    integer :: I                   ! Index for "sons" or "section_ordering"
                                   ! or subtrees of "sons" or subtree of "expect"

    error = max(error,1)
    call StartErrorMessage ( where )
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
    case ( inconsistent_types )
      call output ( 'types are not consistent.', advance = 'yes' )
    case ( inconsistent_units )
      call output ( 'units are not consistent.', advance = 'yes' )
    case ( missing_field )
      call display_string ( field_indices(fields(1)), before='the "' )
      call output ( '" field is required but not present.', advance='yes' )
    case ( no_code_for )
      call output ( 'there is no code to analyze ' )
      call dump_tree_node ( where, 0, advance='yes' )
    case ( no_dot )
      call output ( 'a reference of the form X.Y is not allowed.', &
        advance='yes' )
    case ( no_duplicate_fields )
      call display_string ( field_indices(fields(1)), before='the "' )
      call output ( '" field shall not be specified twice.', advance='yes' )
    case ( no_positional_fields )
      call output ( 'positional fields are not allowed.', advance='yes' )
    case ( no_such_field )
      call display_string ( field_indices(fields(1)), before='a required field "' )
      call output ( '" is absent in the chain of specifications.', &
        advance='yes' )
    case ( no_such_reference )
      call display_string ( sub_rosa(where), before='there is no reference to ' )
      call output ( ' in the field at ' )
      call print_source ( source_ref(sons(1)), advance='yes' )
    case ( not_field_of )
      call display_string ( sub_rosa(where) )
      call display_string ( sub_rosa(subtree(1,fields(1))), &
        & before=' is not a field of ', advance='yes' )
    case ( not_func )
      call display_string ( fields(1) )
      call output ( ' is not a valid function.', advance='yes' )
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
    case ( wrong_arg_type )
      call display_string ( fields(1), before='Argument of ' )
      call output ( ' is not the correct type' )
      if ( present(expect) ) &
        & call display_string ( sub_rosa(expect), before=', expected ' )
      call output ( '.', advance='yes' )
    case ( wrong_num_args )
      call display_string ( fields(1), before='Incorrect number of arguments for ', &
        & advance='yes' )
    case ( wrong_type, wrong_units )
      call output ( 'the "' )
      if ( present(fields) ) then
        call display_string ( sub_rosa(fields(1)) )
      else
        call display_string ( sub_rosa(where) )
      end if
      call output ( '" field has the wrong ' )
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
    integer :: Me = -1   ! String index for trace
    logical :: NO_ARRAY_ALLOWED ! Field is n_field_type and is required to be scalar
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
              do f = 2, nsons(field)
                if ( node_id(subtree(f,field)) == n_unchecked ) then
                  stat = 0
              exit
                end if
                ! First son of n_or is the field name.  All sons of n_or
                ! are types.
                stat = assignBody ( son, subtree(f,field), 1 )
                if ( stat == 0 ) &
              exit
              end do
            end if
            select case ( stat )
            case ( 0 )
            case ( 1 )
              call announce_error ( son, wrong_type, fields=(/son1/), &
                & expect=field )
            case ( 2 )
              call announce_error ( subtree(2,son), no_such_field, &
                & fields=(/ field_look /) )
            case ( 3 )
              call announce_error ( subtree(2,son), no_such_reference, &
                & (/ subtree(1,field_test) /) )
            case ( 4 )
              call announce_error ( subtree(2,son), not_spec, expect=field_last )
            case default
              call announce_error ( son1, wrong_units, fields = (/ son1 /), &
                & expect=stat-10 )
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
      type(decls) :: DECL  ! Declaration of a name
      integer :: J         ! Index of root of "field"
      integer :: Me = -1   ! String index for trace
      integer :: TEST_TYPE ! Used to test the tree-node for a type reference
      integer :: TYPE_DECL ! Tree node of declaration of name of field's type

      call trace_begin ( me, 'AssignBody', root, cond=toggle(con) )

      ! Stat = 0 => Type and units OK
      !      = 1 => wrong type
      !      = 2 => no such field
      !      = 3 => no such reference
      !      = 4 => no such label
      !      > 10 => Expecting units stat-10

      stat = 0 ! Assume normal return status
      select case ( node_id(root) )
      case ( n_dot ) ! label_ref.field
        if ( node_id(field) == n_dot ) then      ! dot allowed
          stat = check_dot ( root, field, start, field_look, field_test )
        else
          stat = 1 ! Wrong type
        end if
      case ( n_identifier )
        if ( node_id(field) /= n_dot ) then      ! Field doesn't need x.y
          do j = start, nsons(field)             ! Try all of the types
            type_decl = decoration(subtree(j,field))
            decl = get_decl(sub_rosa(root),look_for)
            do while ( decl%tree /= null_tree )
              if ( node_id(type_decl) == n_spec_def ) then
                test_type = decoration(subtree(1,decl%tree))
              else
                test_type = decl%tree
              end if
              if ( test_type == type_decl ) then  ! right type
                if ( look_for == enum_value ) then
                  call decorate ( root, decl%units ) ! decorate root with lit#
                else
                  call decorate ( root, decl%tree )  ! decorate root with tree
                end if
                call trace_end ( 'AssignBody', cond=toggle(con) )
    return
              end if
              decl = prior_decl(decl,look_for)
            end do
          end do
        end if
        ! Either field needs x.y, or an allowed type was not found
        stat = 1
      case ( n_array )
        do j = 1, nsons(root)
          stat = assignBody ( subtree(j,root), field, start )
          if ( stat /= 0 ) return
        end do
      case default
        if ( node_id(field) /= n_dot ) then ! Field doesn't need x.y
          stat = expr (root, type, units, value, field, start, field_look, field_test)
          stat = check_field_type(field, type_map(type), units, start)
        else
          stat = 1
        end if
      end select
      call trace_end ( 'AssignBody', cond=toggle(con) )
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

    call trace_begin ( me, 'Check_Dot', root, cond=toggle(con) )

    stat = 3 ! Assume no label referenced in <dot label field>
    type_decl = decoration(subtree(start,field)) ! Required spec
    first = subtree(1,root) ! label in <dot label field>
    last = subtree(2,root)  ! field in <dot label field>
    field_last = sub_rosa(last)
    decl = get_decl(sub_rosa(first),label)
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
              stat = 4 ! No such referenced spec
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
      stat = 2 ! no such field
    9 call trace_end ( 'Check_Deep', stat, cond=toggle(con) )
    end function Check_Deep

  end function Check_Dot

! --------------------------------------------------  CHECK_FIELD  -----
  integer function CHECK_FIELD ( FIELD, SPEC )
  ! Return position in tree of "n_field_type" or "n_field_spec" node of
  ! field definition if "field" is a field of "spec," else return zero.
    integer, intent(in) :: FIELD
    integer, intent(in) :: SPEC
    integer :: I         ! Index of sons of "spec_def"
    integer :: SUB_FIELD ! Sub_rosa = string index of "field"

    sub_field = sub_rosa(field)    ! String index of field
    do i = 2, nsons(spec)
      check_field = subtree(i,spec)
      if ( sub_rosa(subtree(1,check_field)) == sub_field ) return
    end do
    check_field = 0
  end function CHECK_FIELD
! ---------------------------------------------  CHECK_FIELD_TYPE  -----
  integer function CHECK_FIELD_TYPE ( FIELD, TYPE, UNITS, START )
  ! Return zero if 'type' is allowed for 'field' of 'spec' "
  ! or if 'type' is allowed for parameter declared at 'field' ".
  ! Return 1 if type is wrong.
  ! Return 10+required PHYQ_... if type is correct but units are wrong.
    integer, intent(in) :: FIELD   ! Field declaration, from "check_field"
                                   ! or a parameter declaration from
                                   ! "check_section"
    integer, intent(in) :: TYPE    ! A "t_type" name from Init_Tables_Module
    integer, intent(in), optional :: UNITS ! From EXPR
    integer, intent(in), optional :: START
    type(decls) :: DECL            ! Declaration of "type"
    integer :: I                   ! Index of sons of "spec_def"
    integer :: MyStart             ! Either 2 or start
    integer :: REQ_U               ! Required units of a son of FIELD

    myStart = 2
    if ( present(start) ) myStart = start
    check_field_type = merge( 1, 0, &
      &  type < lbound(data_type_indices,1) .or. &
      &  type > ubound(data_type_indices,1) )
    if ( check_field_type > 0 ) return
    decl = get_decl(data_type_indices(type), type_name)
    req_u = decoration(field) / u
    do i = myStart, nsons(field)
      check_field_type = merge( 0, 1, decoration(subtree(i,field)) == decl%tree )
      if ( check_field_type == 0 ) then
        if ( req_u > 0 .and. present(units) ) then
          if ( req_u /= units ) check_field_type = 10 + req_u
        end if
        return
      end if
    end do
    check_field_type = 1
  end function CHECK_FIELD_TYPE
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
        if ( sub_rosa(check_section) == sub_spec ) return
      else ! Must be n_name_def
        if ( sub_rosa(subtree(1,check_section)) == sub_spec ) return
      end if
    end do
    check_section = 0
  end function CHECK_SECTION
! -----------------------------------------------------  DEF_FUNC  -----
  subroutine DEF_FUNC ( ROOT )
  ! Process a definition of a function: Enter it into the declaration
  ! table.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_dt_def )

    integer :: Me = -1   ! String index for trace
    integer :: SON       ! Son of Root

    call trace_begin ( me, 'DEF_FUNC', root, cond=toggle(con) )
    son = subtree(1,root)
    call declare ( sub_rosa(son), 0.0d0, function, decoration(son), root )
    call trace_end ( 'DEF_FUNC', root, cond=toggle(con) )
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
    call declare ( sub_rosa(son), 0.0d0, section, decoration(son), root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( node_id(son) == n_name_def ) then
        gson = subtree(1,son)
        call declare ( sub_rosa(gson), 0.0d0, named_value, &
                       decoration(gson), son )
        gson = subtree(2,son)
        decl = get_decl(sub_rosa(gson),type_name)
        call decorate ( gson, decl%tree )    ! the dt_def
      end if
    end do
    call trace_end ( 'DEF_SECTION', root, cond=toggle(con) )
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
    call declare ( sub_rosa(spec_name), 0.0d0, spec, decoration(spec_name), &
                   root )
!   call decorate ( spec_name, root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      son=son ! Without this, the compiler inexplicably uses "root" instead
              ! of "son" below ???
      field_name = subtree(1,son)
      decl = get_decl(sub_rosa(field_name),field)
      if ( decl%tree == null_tree ) then ! don't make several
        call redeclare ( sub_rosa(field_name), decl%value, field, &
                         decoration(field_name), son )
      else ! need to make several -- one for each
        call declare ( sub_rosa(field_name), decl%value, field, &
                       decoration(field_name), son )
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
      call trace_end ( 'DEF_ONE_SPEC', root, cond=toggle(con) )
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
    call declare ( sub_rosa(son), 0.0d0, type_name, decoration(son), root )
!   call decorate ( son, root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      call declare ( sub_rosa(son), decoration(son) + 0.0d0, enum_value, &
                     decoration(son), root )
!     call decorate ( son, root )
    end do
    call trace_end ( 'DEF_TYPE', cond=toggle(con) )
  end subroutine DEF_TYPE

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
          stat = check_field_type(check,type_map(type),units)
          if ( stat == 1 ) then
            call announce_error ( son2, wrong_type, fields=(/son1/), expect=check )
          else if ( stat > 1 ) then
            call announce_error ( son2, wrong_units, fields=(/son1/), expect=stat-10 )
          end if
        end if
      end if
    end if
    call trace_end ( 'EQUAL', cond=toggle(con) )
  end subroutine EQUAL

! ---------------------------------------------------------  EXPR  -----
  recursive integer function EXPR ( ROOT, TYPE, UNITS, VALUE, &
    &                               FIELD, START, FIELD_LOOK, FIELD_TEST ) &
    & result ( Stat )
  ! Analyze an expression, check units of everything in it, return its
  ! units.
    integer, intent(in) :: ROOT    ! Root of expression subtree
    integer, intent(out) :: TYPE   ! Type of the expression value
    integer, intent(out) :: UNITS  ! Units of expression value
    double precision, intent(out) :: VALUE   ! Expression value, if any
    integer, intent(in), optional :: FIELD    ! What TO_CHECK should be
    integer, intent(in), optional :: START    ! of sons of FIELD
    integer, intent(out), optional :: Field_Look ! Last name of x.y, for error message
    integer, intent(out), optional :: Field_Test ! Allowed field name, for error message

    ! Stat is zero unless TO_CHECK, FIELD, and START are present,
    ! node_id(ROOT) == n_dot, and the operands of n_dot are wrong for
    ! FIELD.

    integer :: Arg_Tree            ! Root of tree of allowed argument types
    type(decls) :: DECL            ! Declaration record for "root"
    integer :: I
    integer :: ME                  ! node_id(root)
    integer :: SON1, SON2          ! Sons of "root"
    integer :: STRING              ! sub_rosa(root)
    integer :: Trace = -1          ! String index for trace
    integer :: TYPE2               ! Type for a son of "root"
    integer :: UNITS2              ! Units for a son of "root"
    double precision :: VALUE2     ! Value for a son of "root"

    call trace_begin ( trace, 'EXPR', root, cond=toggle(con) )
    stat = 0 ! Assume status is OK
    units = phyq_dimensionless     ! default
    value = 0.0d0                  ! default
    me = node_id(root)
    select case ( me )
    case ( n_identifier, n_number, n_string ) ! ----------------------------
      string = sub_rosa(root)
      if ( .not. declared(string) ) then
        select case ( me )
        case ( n_identifier )
          call declare ( string, 0.0d0, undeclared, phyq_invalid, root )
        case ( n_number )
          call declare ( string, float_value(sub_rosa(root)), num_value, &
                         phyq_dimensionless, root )
        case ( n_string )
          call declare ( string, 0.0d0, str_value, phyq_invalid, root )
        end select
      end if
      decl = declaration(string)
      type = decl%type
      units = decl%units
      value = decl%value
    case ( n_unit ) ! ------------------------------------------------------
      son1 = subtree(1,root)
      stat = expr (son1, type, units, value, field, start, field_look, field_test)
      if ( stat /= 0 ) &
  return
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
    case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less ) ! ------
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      if ( type /= type2 ) then
        call local_error ( root, inconsistent_types, (/ son1, son2 /) )
      end if
! Don't test units here.  The field might use : for something other than range.
!       if ( type == num_value .and. type2 == num_value .and. &
!            units /= phyq_dimensionless .and. &
!            units2 /= phyq_dimensionless .and. &
!            units /= units2 ) then
!         call local_error ( root, inconsistent_units, (/ son1, son2 /) )
!       end if
      if ( type == num_value ) type = range
    case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, n_not_equal )
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      if ( type == num_value .and. type2 == num_value .and. &
           units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call local_error ( root, inconsistent_units, (/ son1, son2 /) )
      end if
      type = log_value
    case ( n_and, n_or ) ! -------------------------------------------------
      son1 = subtree(1,root); son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      if ( type /= str_value .and. type /= log_value ) then
        call local_error ( son1, not_string )
      else if ( type2 /= str_value .and. type2 /= log_value ) then
        call local_error ( son2, not_string )
      else
        if ( me == n_and ) then
          value = value * value2
        else ! node_id(root) == n_or
          value = max(value, value2)
        end if
      end if
      type = log_value
    case ( n_not )
      son1 = subtree(1,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      if ( type /= str_value .and. type /= log_value ) then
        call local_error ( son1, not_string )
      else
        value = 1 - value
      end if
      type = log_value
    case ( n_plus, n_minus ) ! ---------------------------------------------
      son1 = subtree(1,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      if ( nsons(root) > 1 ) then
        son2 = subtree(2,root)
        stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
        if ( stat /= 0 ) &
  return
        if ( type == num_value .and. type2 == num_value .and. &
           units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call local_error ( root, inconsistent_units, (/ son1, son2 /) )
        else
          if ( me == n_plus ) then
            value = value + value2
          else !  me == n_minus
            value = value - value2
          end if
        end if
      else if ( me == n_minus ) then
        value = - value
      end if
    case ( n_mult ) ! ------------------------------------------------------
      son1 = subtree(1,root) ; son2 = subtree(2,root)
      stat = expr ( son1, type, units, value, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
      stat = expr ( son2, type, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
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
  return
      stat = expr ( son2, type2, units2, value2, field, start, field_look, field_test )
      if ( stat /= 0 ) &
  return
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
        stat = expr ( son1, type, units, value2, field_look, field_test )
        if ( stat /= 0 ) &
  return
        if ( units /= phyq_dimensionless ) then
          call local_error ( son1, not_unitless )
        else if ( type /= num_value ) then
          call local_error ( root, not_numeric, (/ son1 /) )
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
          if ( nsons(decl%tree) == 1 ) & ! No argument # or type checking
  return
          arg_tree = subtree(2,decl%tree)
          if ( node_id(arg_tree) /= n_arg_def ) & ! No argument # or type checking
  return
          if ( nsons(root) /= 1 + nsons(arg_tree) ) then
            call local_error ( root, wrong_num_args, fields=(/string/) )
          else
            do i = 1, nsons(arg_tree)
              stat = expr ( subtree(i+1,root), type, units, value, &
                          & field, start, field_look, field_test )
              if ( stat /= 0 ) &
  return
              if ( type_map(type) /= decoration(subtree(i,arg_tree)) ) then
                call local_error ( subtree(i+1,root), wrong_arg_type, &
                  & fields=(/string/), expect=subtree(i,arg_tree) )
              else if ( units /= phyq_dimensionless ) then
                call local_error ( subtree(i+1,root), not_unitless, fields=(/string/) )
              end if
            end do
  return
          end if
        end if
        decl = prior_decl(decl,function)
      end do
      call local_error ( son1, not_func, fields=(/string/) )
    case ( n_dot ) ! -------------------------------------------------------
      stat = check_dot ( root, field, start, field_look, field_test )
    case default ! ---------------------------------------------------------
      call local_error ( root, no_code_for )
    end select
    call trace_end ( 'EXPR', cond=toggle(con) )

  contains

    subroutine LOCAL_ERROR ( WHERE, CODE, SONS, FIELDS, EXPECT )
      integer, intent(in) :: WHERE   ! Tree node where error was noticed
      integer, intent(in) :: CODE    ! Code for error message
      integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                     ! "where".  If they're pseudo_terminal,
                                     ! their declarations are dumped.
      integer, intent(in), optional :: FIELDS(:) ! Field indices
      integer, intent(in), optional :: EXPECT    ! Expected type
      call announce_error ( where, code, sons, fields, expect )
      type = empty
      units = phyq_invalid
    end subroutine LOCAL_ERROR

  end function EXPR

! -------------------------------------------------------  ONE_CF  -----
  subroutine ONE_CF ( ROOT )
  ! Analyze one configuration, with abstract syntax tree rooted at ROOT.
  ! The root is an N_CF node.
    integer, intent(in) :: ROOT

    type(decls) :: DECL            ! Declaration of "name" on "begin name"
    integer :: GSON1, GSON2        ! Grandsons of ROOT being analyzed
    integer :: I, N                ! Loop inductor, number of sons
    integer :: Me = -1             ! String index for trace
    integer :: SON1, SONN          ! Sons of ROOT being analyzed
    integer :: STRING1, STRINGN    ! Sub_Rosa of first and last sons

    call trace_begin ( me, 'ONE_CF', root, cond=toggle(con) )
    n = nsons(root)
    son1 = subtree(1,root) ;   sonn = subtree(n,root)
    string1 = sub_rosa(son1) ; stringn = sub_rosa(sonn)
    if ( string1 /= stringn ) then
      call output ( '***** Names on BEGIN at ' )
      call print_source ( source_ref(subtree(1,root)), advance='yes' )
      call output ( '*****      and END at ' )
      call print_source ( source_ref(subtree(n,root)) )
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
      sonn = subtree(i,root)
      select case ( node_id(sonn) )
      case ( n_equal )                ! x = expr
        call decorate ( sonn, son1 )  ! show equal the section name
        if ( decoration(decl%tree) /= no_check_eq ) call equal ( sonn )
      case ( n_named )                ! label:x,(y[=expr])+
        gson1 = subtree(1,sonn)       ! Label
        gson2 = subtree(2,sonn)       ! spec_args tree
        string1 = sub_rosa(gson1)
        decl = get_decl ( string1, label )
        if ( decl%type == label ) then
          call announce_error ( gson1, already_declared )
        else
          call declare ( string1, 0.0d0, label, phyq_invalid, gson2 )
        end if
        call decorate ( gson1, gson2 )
        call decorate ( gson2, son1 ) ! show spec_args the section name
        call spec_args ( gson2 )
      case ( n_spec_args )            ! x,(y[=expr])+
        call decorate ( sonn, son1 )  ! show spec_args the section name
        call spec_args ( sonn )
      case default
        call announce_error ( sonn, no_code_for )
      end select
    end do
    call trace_end ( 'ONE_CF', cond=toggle(con) )
  end subroutine ONE_CF
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
