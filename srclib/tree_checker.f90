! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE_CHECKER

! Traverse the tree output by the parser, which includes definitions
! put into it by init_tables_module before the parser runs.  Check
! that things fit together.  Fill in the declaration table.  Decorate
! cross-references within the tree with anything that might be useful
! in table_generator.

  use DECLARATION_TABLE, only: DECLARATION, DECLARE, DECLARED, DECLS, &
                               DUMP_1_DECL, ENUM_VALUE, FIELD, GET_DECL, &
                               LABEL, LOG_VALUE, NAMED_VALUE, NULL_DECL, &
                               NUM_VALUE, PRIOR_DECL, RANGE, REDECLARE, &
                               SECTION, SPEC, STR_VALUE, UNDECLARED, &
                               TYPE_MAP, TYPE_NAME, UNITS_NAME
  use INIT_TABLES_MODULE, only: DATA_TYPE_INDICES, FIELD_FIRST, FIELD_INDICES, &
                                FIELD_LAST, PHYQ_DIMENSIONLESS, &
                                PHYQ_INVALID, SECTION_FIRST, &
                                SECTION_INDICES, SECTION_LAST, &
                                SECTION_ORDERING, T_BOOLEAN
  use INTRINSIC, only: ALL_FIELDS, NO_DUP, NO_POSITIONAL, REQ_FLD
  use LEXER_CORE, only: PRINT_SOURCE
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING, FLOAT_VALUE
  use TOGGLES, only: CON, TOGGLE
  use TRACE_M, only: DEPTH, TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, &
                  NODE_KIND, NSONS, NULL_TREE, PSEUDO, SOURCE_REF, &
                  SUB_ROSA, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  implicit NONE
  private

  public :: CHECK_TREE

! -----     Private declarations     -----------------------------------
  integer, private :: ERROR   ! 0 => No errors

! Error codes for "announce_error"
  integer, private, parameter :: ALREADY_DECLARED = 1
  integer, private, parameter :: INCONSISTENT_TYPES = ALREADY_DECLARED + 1
  integer, private, parameter :: INCONSISTENT_UNITS = INCONSISTENT_TYPES + 1
  integer, private, parameter :: MISSING_FIELD = INCONSISTENT_UNITS + 1
  integer, private, parameter :: NO_CODE_FOR = MISSING_FIELD + 1
  integer, private, parameter :: NO_DOT = NO_CODE_FOR + 1
  integer, private, parameter :: NO_DUPLICATE_FIELDS = NO_DOT + 1
  integer, private, parameter :: NO_POSITIONAL_FIELDS = NO_DUPLICATE_FIELDS + 1
  integer, private, parameter :: NO_SUCH_FIELD = NO_POSITIONAL_FIELDS + 1
  integer, private, parameter :: NO_SUCH_REFERENCE = NO_SUCH_FIELD + 1
  integer, private, parameter :: NOT_FIELD_OF = NO_SUCH_REFERENCE + 1
  integer, private, parameter :: NOT_NAME = NOT_FIELD_OF + 1
  integer, private, parameter :: NOT_NAME_OR_STRING = NOT_NAME + 1
  integer, private, parameter :: NOT_SECTION = NOT_NAME_OR_STRING + 1
  integer, private, parameter :: NOT_SPEC = NOT_SECTION + 1
  integer, private, parameter :: NOT_STRING = NOT_SPEC + 1
  integer, private, parameter :: NOT_UNITS = NOT_STRING + 1
  integer, private, parameter :: OUT_OF_PLACE = NOT_UNITS + 1
  integer, private, parameter :: SECTION_ORDER = OUT_OF_PLACE + 1
  integer, private, parameter :: WRONG_TYPE = SECTION_ORDER + 1

  logical, private :: ALL_FIELDS_FLAG   ! All fields are required
  integer, private :: CURRENT_SECTION = section_first - 1
  logical, private :: GOT(field_first:field_last)
  logical, private :: NO_DUP_FLAG       ! Duplicate named fields prohibited

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
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
    integer :: NUM_SECTIONS   ! Number of begin ... ends
    integer :: SON            ! Son of root
    depth = 0
    error = 0
    num_sections = 0
    if ( present(first_section) ) first_section = 0
    if ( toggle(con) ) call trace_begin ( 'CHECK_TREE', root )
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
      case ( n_section );   call def_section ( son )
      case ( n_spec_def );  call def_spec ( son )
      case default ;        call announce_error ( son, no_code_for )
      end select
    end do
    error_flag = error
    if ( present(how_many_sections) ) how_many_sections = num_sections
    if ( toggle(con) ) call trace_end ( 'CHECK_TREE' )
  end subroutine CHECK_TREE
! =====     Private Procedures     =====================================
! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, SONS, FIELDS )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                   ! "where".  If they're pseudo_terminal,
                                   ! their declarations are dumped.
    integer, intent(in), optional :: FIELDS(:) ! Field indices
    integer :: I                   ! Index for "sons" or "section_ordering"
                                   ! or subtrees of "sons"

    error = max(error,1)
    call output ( '***** At ', from_where = "type checker" )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( already_declared )
      call dump_tree_node ( where, 0 )
      call output ( ' is already defined.', advance='yes' )
    case ( inconsistent_types )
      call output ( 'types are not consistent.', advance = 'yes' )
    case ( inconsistent_units )
      call output ( 'units are not consistent.', advance = 'yes' )
    case ( missing_field )
      call output ( 'the "' )
      call display_string ( field_indices(fields(1)) )
      call output ( '" field is required but not present.', advance='yes' )
    case ( no_code_for )
      call output ( 'there is no code to analyze ' )
      call dump_tree_node ( where, 0, advance='yes' )
    case ( no_dot )
      call output ( 'a reference of the form X.Y is not allowed.', &
        advance='yes' )
    case ( no_duplicate_fields )
      call output ( 'the "' )
      call display_string ( field_indices(fields(1)) )
      call output ( '" field shall not be specified twice.', advance='yes' )
    case ( no_positional_fields )
      call output ( 'positional fields are not allowed.', advance='yes' )
    case ( no_such_field )
      call output ( 'a required field "' )
      call display_string ( field_indices(sons(1)) )
      call output ( '" is absent in the chain of specifications.', &
        advance='yes' )
    case ( no_such_reference )
      call output ( 'there is no reference to ' )
      call display_string ( sub_rosa(where) )
      call output ( ' in the field at ' )
      call print_source ( source_ref(sons(1)), advance='yes' )
    case ( not_field_of )
      call display_string ( sub_rosa(where) )
      call output ( ' is not a field of ' )
      call display_string ( sub_rosa(subtree(1,fields(1))), advance='yes' )
    case ( not_name )
      call output ( 'is not a name.', advance = 'yes' )
    case ( not_name_or_string )
      call output ( 'is not a name or a string.', advance = 'yes' )
    case ( not_section )
      call display_string ( sub_rosa(where) )
      call output ( ' is not a section name.', advance = 'yes' )
    case ( not_spec )
      call display_string ( sub_rosa(where) )
      call output ( ' is not a spec name.', advance = 'yes' )
    case ( not_string )
      call output ( 'is not a string or of logical type.', &
                    advance = 'yes' )
    case ( not_units )
      call dump_tree_node ( where, 0 )
      call output ( ' is not a units name.', advance = 'yes' )
    case ( out_of_place )
      call display_string ( sub_rosa(where) )
      call output ( ' is not allowed in a ' )
      call display_string ( sub_rosa(sons(1)) )
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
      call output ( '', advance='yes' )
    case ( wrong_type )
      call output ( 'the "' )
      call display_string ( sub_rosa(where) )
      call output ( '" field has the wrong type of associated value.', &
        & advance='yes'  )
      if ( present(sons) ) then
        call output ( '         Expected' )
        if ( nsons(sons(1)) > 2 ) call output ( ' one of' )
        do i = 2, nsons(sons(1))
          if ( i > 2 ) call output ( ',' )
          call output ( ' ' )
          if ( i > 2 .and. i == nsons(sons(i)) ) call output ( 'or ' )
          call display_string ( sub_rosa(subtree(i,sons(1))) )
        end do
        call output ( '', advance='yes' )
      end if
    case default
      call output ( 'No message in TREE_CHECKER for error code ' )
      call output ( code, advance='yes' )
      stop
    end select
    if ( present(sons) ) then
      do i = 1, size(sons)
        if ( node_kind(sons(i)) == pseudo ) &
          call dump_1_decl ( sub_rosa(sons(i)) )
      end do
    end if
  end subroutine ANNOUNCE_ERROR
! ----------------------------------------------------  ASSIGN  -----
  subroutine ASSIGN ( ROOT, TYPE, UNITS, VALUE )
  ! Analyze a son of "n_spec_args" of the form "name = expr+",
  ! starting at ROOT
    integer, intent(in) :: ROOT            ! Index of the n_asg node
    integer, intent(out) :: TYPE, UNITS    ! Output from "expr"
    double precision, intent(out) :: VALUE ! Output from "expr"

    integer :: FIELD     ! Tree node of field's declaration
    integer :: FIELD_LIT ! f_... for a field
    integer :: I         ! Index of son of "root"
    integer :: LOOK_FOR  ! Look for an enum_value or a spec?
    integer :: SON1, SON ! Sons of "root"
    integer :: SPEC_DECL ! Tree node of the spec's declaration

    if ( toggle(con) ) call trace_begin ( 'ASSIGN', root )
    son1 = subtree(1,root)
    if ( node_id(son1) == n_identifier ) then
      spec_decl = decoration(root)
      field = check_field(son1,spec_decl)   ! Is field a field of spec?
      if ( field == 0 ) then
        call announce_error ( son1, not_field_of, &
          & fields=(/ spec_decl /) )
      else
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
        do i = 2, nsons(root)
          son = subtree(i,root)
          call assignBody ( son )
        end do ! i = 2, nsons(root)
      end if
    else
      call announce_error ( son1, not_name )
    end if
    if ( toggle(con) ) call trace_end ( 'ASSIGN' )
  contains
    recursive subroutine AssignBody ( Son )
      integer, intent(in) :: Son   ! of n_asg or n_array
      type(decls) :: DECL  ! Declaration of a name
      integer :: FIELD_LOOK ! A field being sought during n_dot checking
      integer :: FIELD_REF ! The value of a field -- a ref to a n_spec_arg
      integer :: FIELD_TEST ! A son of Field_Ref
      integer :: GSON      ! Son of Son
      integer :: I         ! subtree index and loop inductor        
      integer :: J         ! Index of son of "field"                
      integer :: K         ! Index of a son of a n_spec_args
      integer :: TEST_TYPE ! Used to test the tree-node for a type reference
      integer :: TYPE_DECL ! Tree node of declaration of name of field's type

      select case ( node_id(son) )
      case ( n_dot ) ! label_ref.field
        gson = subtree(1,son)
        type_decl = decoration(subtree(2,field)) ! Required spec
        decl = get_decl(sub_rosa(gson),label)
        do while ( decl%tree /= null_tree )
          test_type = decoration(subtree(1,decl%tree)) ! spec's index
          if ( test_type == type_decl ) then  ! right kind of spec
            field_ref = decl%tree
            call decorate ( gson, field_ref ) ! decorate label_ref with tree
m:              do j = 3, nsons(field)
            ! This loop assumes there is only one field of the required
            ! name.  If it is desired to search through several fields,
            ! it will be necessary to have a stack to keep track of the
            ! subtree indices.  It goes through the field names that are
            ! given by sons 3-n of the n_dot vertex of the definition,
            ! looking for fields of the same name, starting at the first
            ! son of the n_dot vertex in the input, and thence from the
            ! decoration of the brother of the found field name.
              field_look = decoration(subtree(j,field))
              do k = 2, nsons(field_ref)
                field_test = subtree(k,field_ref)
                if ( node_id(field_test) == n_asg ) then
                  if ( decoration(subtree(1,field_test)) == &
                    &  field_look ) then
                    ! Get the next n_spec_arg tree in the chain
                    field_ref = decoration(subtree(2,field_test))
            cycle m
                  end if
                end if
              end do ! k
              call announce_error ( gson, no_such_field, (/ field_look /) )
    return
            end do m ! j
            ! The final field_test is the parent n_asg of the field
            ! that should have the same name as the second son of the
            ! n_dot vertex in the input.
            gson = subtree(2,son)
            field_look = sub_rosa(gson)
            do k = 2, nsons(field_test)
              if ( sub_rosa(subtree(k,field_test)) == field_look ) then
                call decorate ( gson, subtree(k,field_test) )
    return
              end if
            end do ! k
            call announce_error ( gson, no_such_reference, &
              & (/ subtree(1,field_test) /) )
    return
          end if
          decl = prior_decl(decl,label)
        end do
        call announce_error ( son1, wrong_type )
      case ( n_identifier )
        do j = 2, nsons(field)                 ! Try all of the types
          type_decl = decoration(subtree(j,field))
          decl = get_decl(sub_rosa(son),look_for)
          do while ( decl%tree /= null_tree )
            if ( node_id(type_decl) == n_spec_def ) then
              test_type = decoration(subtree(1,decl%tree))
            else
              test_type = decl%tree
            end if
            if ( test_type == type_decl ) then  ! right type
              if ( look_for == enum_value ) then
                call decorate ( son, decl%units ) ! decorate son with lit#
              else
                call decorate ( son, decl%tree )  ! decorate son with tree
              end if
    return
            end if
            decl = prior_decl(decl,look_for)
          end do
        end do
        call announce_error ( son1, wrong_type )
      case ( n_array )
        do j = 1, nsons(son)
          call assignBody ( subtree(j,son) )
        end do
      case default
        call expr ( son, type, units, value )
        if ( .not. check_field_type(field,type_map(type)) ) then
          call announce_error ( son1, wrong_type, sons = (/ field /) )
        end if
      end select
    end subroutine AssignBody
  end subroutine ASSIGN
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
  logical function CHECK_FIELD_TYPE ( FIELD, TYPE )
  ! Return " 'type' is allowed for 'field' of 'spec' "
  ! or " 'type' is allowed for parameter declared at 'field' ".
    integer, intent(in) :: FIELD   ! Field declaration, from "check_field"
                                   ! or a parameter declaration from
                                   ! "check_section"
    integer, intent(in) :: TYPE    ! A "t_type" name from Init_Tables_Module
    type(decls) :: DECL            ! Declaration of "type"
    integer :: I                   ! Index of sons of "spec_def"

    decl = get_decl(data_type_indices(type), type_name)
    do i = 2, nsons(field)
      check_field_type = decoration(subtree(i,field)) == decl%tree
      if ( check_field_type ) return
    end do
    check_field_type = .false.
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
    integer :: SON       ! Son of Root

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
  end subroutine DEF_SECTION
! -----------------------------------------------------  DEF_SPEC  -----
  subroutine DEF_SPEC ( ROOT )
  ! Process a definition of a structure:  Enter the structure and field
  ! names in the declaration table.  Decorate the structure name with
  ! the root. Decorate each field_type node with the structure name.
  ! Decorate each field name with its parent field_type node.  Decorate
  ! each field type name with the dt_def node for the type.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_spec_def )

    type(decls) :: DECL       ! Current declaration of "son"
    integer :: FIELD_NAME, FIELD_TYPE   ! Grandsons
    integer :: I, J           ! Loop inductors
    integer :: SON            ! I'th son of "root"
    integer :: SPEC_NAME      ! First son of "root"

    if ( toggle(con) ) call trace_begin ( 'DEF_SPEC', root )
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
      select case ( node_id(son) )
      case ( n_field_type )
        do j = 2, nsons(son)
          field_type = subtree(j,son)
          decl = get_decl(sub_rosa(field_type),type_name)
          call decorate ( field_type, decl%tree )   ! the dt_def
        end do
      case ( n_field_spec )
        do j = 2, nsons(son)
          field_type = subtree(j,son)
          decl = get_decl(sub_rosa(field_type),spec)
          call decorate ( field_type, decl%tree )   ! the spec_def
        end do
      case ( n_dot )
        field_type = subtree(2,son)
        decl = get_decl(sub_rosa(field_type),spec)
        call decorate ( field_type, decl%tree )   ! the spec_def
        ! The rest of the sons are field names, for which new decorations
        ! Won't help -- in fact, the f_field_name's index is best.
      end select
    end do
    if ( toggle(con) ) call trace_end ( 'DEF_SPEC' )
  end subroutine DEF_SPEC
! -----------------------------------------------------  DEF_TYPE  -----
  subroutine DEF_TYPE ( ROOT )
  ! Process a definition of a type:  Enter the names in the declaration
  ! table.  The type name and literals are decorated with their indices
  ! in Init_Tables. These decorations are entered into the declaration
  ! table.  Then each one is decorated with Root.
    integer, intent(in) :: ROOT    ! Root of tree being worked ( n_dt_def )

    integer :: I              ! Loop inductor
    integer :: SON            ! I'th son of "root"

    if ( toggle(con) ) call trace_begin ( 'DEF', root )
    son = subtree(1,root)
    call declare ( sub_rosa(son), 0.0d0, type_name, decoration(son), root )
!   call decorate ( son, root )
    do i = 2, nsons(root)
      son = subtree(i,root)
      call declare ( sub_rosa(son), decoration(son) + 0.0d0, enum_value, &
                     decoration(son), root )
!     call decorate ( son, root )
    end do
    if ( toggle(con) ) call trace_end ( 'DEF' )
  end subroutine DEF_TYPE
! --------------------------------------------------------  EQUAL  -----
  subroutine EQUAL ( ROOT )
  ! Analyze a specification of the form NAME = EXPR, starting at ROOT
    integer, intent(in) :: ROOT
    integer :: CHECK          ! Subtree where NAME is declared
    type(decls) :: DECL       ! Declaration of son1
    integer :: SECT           ! The section the parameter appears in
    integer :: SON1, SON2     ! Sons of Root
    integer :: TYPE           ! Output from "expr", not otherwise used
    integer :: TYPE_DECL      ! Tree node of declaration of name of param's type
    integer :: UNITS          ! Output from "expr", not otherwise used
    double precision :: VALUE ! Output from "expr", not otherwise used
    if ( toggle(con) ) call trace_begin ( 'EQUAL', root )
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
              call announce_error ( son1, wrong_type )
          exit
            end if
            if ( decl%tree == type_decl ) then ! right type
              call decorate ( son2, decl%units )! decorate son2 with lit#
          exit
            end if
            decl = prior_decl(decl,enum_value)
          end do
        else
          call expr( son2, type, units, value )
          if ( .not. check_field_type(check,type_map(type)) ) then
            call announce_error ( son1, wrong_type )
          end if
        end if
      end if
    end if
    if ( toggle(con) ) call trace_end ( 'EQUAL' )
  end subroutine EQUAL
! ---------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR ( ROOT, TYPE, UNITS, VALUE )
  ! Analyze an expression, check units of everything in it, return its
  ! units.
    integer, intent(in) :: ROOT    ! Root of expression subtree
    integer, intent(out) :: TYPE   ! Type of the expression value
    integer, intent(out) :: UNITS  ! Units of expression value
    double precision, intent(out) :: VALUE   ! Expression value, if any

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: I                   ! Subtree index, loop inductor
    integer :: ME                  ! node_id(root)
    integer :: SON1, SON2          ! Sons of "root"
    integer :: STRING              ! sub_rosa(root)
    integer :: TYPE2               ! Type for a son of "root"
    integer :: UNITS2              ! Units for a son of "root"
    double precision :: VALUE2     ! Value for a son of "root"

    if ( toggle(con) ) call trace_begin ( 'EXPR', root )
    units = phyq_dimensionless     ! default
    value = 0.0d0                  ! default
    me = node_id(root)
    select case ( me )
    case ( n_identifier, n_number, n_string )
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
    case ( n_unit )
      son1 = subtree(1,root)
      call expr ( son1, type, units, value )
      son2 = subtree(2,root)
      decl = get_decl(sub_rosa(subtree(2,root)), units_name)
      if ( decl%type /= units_name ) then
        call announce_error ( son2, not_units )
      else if ( units /= phyq_dimensionless .and. units /= decl%units ) then
        call announce_error ( root, inconsistent_units, (/ son1, son2 /) )
      else
        units = decl%units
        value = value * decl%value
      end if
    case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less )
      son1 = subtree(1,root); son2 = subtree(2,root)
      call expr ( son1, type, units, value )
      call expr ( son2, type2, units2, value2 )
      if ( type /= type2 ) then
        call announce_error ( root, inconsistent_types, (/ son1, son2 /) )
      end if
      if ( units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call announce_error ( root, inconsistent_units, (/ son1, son2 /) )
      end if
      if ( type == num_value ) type = range
    case ( n_and, n_or )
      son1 = subtree(1,root); son2 = subtree(2,root)
      call expr ( son1, type, units, value )
      call expr ( son2, type2, units2, value2 )
      if ( type /= str_value .and. type /= log_value ) then
        call announce_error ( son1, not_string )
      else if ( type2 /= str_value .and. type2 /= log_value ) then
        call announce_error ( son2, not_string )
      else
        if ( me == n_and ) then
          value = value * value2
        else ! node_id(root) == n_or
          value = max(value, value2)
        end if
      end if
      type = log_value
    case ( n_plus, n_minus )
      son1 = subtree(1,root)
      call expr ( son1, type, units, value )
      if ( nsons(root) > 1 ) then
        son2 = subtree(2,root)
        call expr ( son2, type2, units2, value2 )
        if ( units /= units2 ) then
          call announce_error ( root, inconsistent_units, (/ son1, son2 /) )
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
    case ( n_mult )
      son1 = subtree(1,root) ; son2 = subtree(2,root)
      call expr ( son1, type, units, value )
      call expr ( son2, type, units2, value2 )
      if ( units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call announce_error ( root, inconsistent_units, (/ son1, son2 /) )
      else
        value = value * value2
      end if
    case ( n_div )
      son1 = subtree(1,root) ; son2 = subtree(2,root)
      call expr ( son1, type, units, value )
      call expr ( son2, type2, units2, value2 )
      if ( units /= phyq_dimensionless .and. &
           units2 /= phyq_dimensionless .and. &
           units /= units2 ) then
        call announce_error ( root, inconsistent_units, (/ son1, son2 /) )
      else
        value = value / value2
      end if
    case ( n_dot )
      call announce_error ( root, no_dot )
    case default
      call announce_error ( root, no_code_for )
    end select
    if ( toggle(con) ) call trace_end ( 'EXPR' )
  end subroutine EXPR
! -------------------------------------------------------  ONE_CF  -----
  subroutine ONE_CF ( ROOT )
  ! Analyze one configuration, with abstract syntax tree rooted at ROOT.
  ! The root is an N_CF node.
    integer, intent(in) :: ROOT

    type(decls) :: DECL            ! Declaration of "name" on "begin name"
    integer :: GSON1, GSON2        ! Grandsons of ROOT being analyzed
    integer :: I, N                ! Loop inductor, number of sons
    integer :: SON1, SONN          ! Sons of ROOT being analyzed
    integer :: STRING1, STRINGN    ! Sub_Rosa of first and last sons

    if ( toggle(con) ) call trace_begin ( 'ONE_CF', root )
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
        call equal ( sonn )
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
    if ( toggle(con) ) call trace_end ( 'ONE_CF' )
  end subroutine ONE_CF
! ------------------------------------------------------  SET_ONE  -----
  subroutine SET_ONE ( ROOT )
  ! Analyze a field of a specification of the form "/name" --
  ! make sure "name" is a field of the specification.
    integer, intent(in) :: ROOT    ! Index of the "n_set_one" tree node
    integer :: FIELD               ! Tree node of field's declaration
    integer :: FIELD_LIT           ! f_... for a field
    integer :: SON                 ! Son of n_set_one tree node
    integer :: SPEC_DECL           ! Tree node of the spec's declaration
    if ( toggle(con) ) call trace_begin ( 'SET_ONE', root )
    son = subtree(1,root)
    spec_decl = decoration(root)
    field = check_field(son,spec_decl)
    if ( field == 0 ) then
      call announce_error ( son, not_field_of, fields=(/ spec_decl /) )
      call dump_1_decl ( sub_rosa(son) )
    else
      if ( check_field_type(field,t_boolean) ) then
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
    if ( toggle(con) ) call trace_end ( 'SET_ONE' )
  end subroutine SET_ONE
! ----------------------------------------------------  SPEC_ARGS  -----
  subroutine SPEC_ARGS ( ROOT )
  ! Analyze a specification of the form NAME (, EXPR ( = EXPR)? )+
  ! starting at ROOT
    integer, intent(in) :: ROOT

    integer :: FIELD_LIT      ! f_... for a field
    integer :: FLAGS          ! Flags from decoration(spec_decl)
    integer :: I              ! Loop inductor
    integer :: SECT           ! The section the spec appears in
    integer :: SON            ! I'th son of "root"
    type(decls) :: SPEC_DECL  ! Declaration of first son, if any
    integer :: TYPE, UNITS    ! Output from "expr"
    double precision :: VALUE ! Output from "expr"

    if ( toggle(con) ) call trace_begin ( 'SPEC_ARGS', root )
    sect = decoration(root)
    spec_decl = get_decl( sub_rosa(sect), section )   ! Known to be a
                             ! valid section name if we get this far
    son = subtree(1,root)
    ! Is spec allowed in section?
    if ( check_section(son,spec_decl%tree) == 0 ) then
      call announce_error ( son, out_of_place, (/ sect /) )
    else
      spec_decl = get_decl( sub_rosa(son), spec ) ! Decl of son of type "spec"
      if ( spec_decl%type == null_decl ) then
        call announce_error ( son, not_spec )
      else
        call decorate ( son, spec_decl%tree )
        flags = decoration(spec_decl%tree)
        all_fields_flag = mod(flags/all_fields,2) .ne. 0
        no_dup_flag = mod(flags/no_dup,2) .ne. 0
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
              call expr ( son, type, units, value )
              ! ??? Generate some table here
            end if
          end select
        end do
        do i = 2, nsons(spec_decl%tree)
          son = subtree(i,spec_decl%tree)
          field_lit = decoration(subtree(1,son))
          if ( .not. got(field_lit) ) then
            if ( all_fields_flag .or. mod(decoration(son)/req_fld,2) /= 0) &
              & call announce_error ( root, missing_field, &
                & fields=(/ field_lit /) )
          end if
        end do
      end if
    end if
    if ( toggle(con) ) call trace_end ( 'SPEC_ARGS' )
  end subroutine SPEC_ARGS
end module TREE_CHECKER

! $Log$
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
