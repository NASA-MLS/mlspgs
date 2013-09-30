! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module PARSER

! Parser for L2CF.

! cf -> one_cf+ 'EOF'
! one_cf -> ('begin' 'name' 'EOS' spec+ 'end' 'name' )? 'EOS'
! one_cf -> #include "file" 'EOS'
! spec -> 'name' spec_rest 'EOS'
! spec -> #include "file" 'EOS'
! spec -> 'EOS'
! spec_rest -> \lambda
! spec_rest -> '=' expr +
! spec_rest -> ( ',' spec_list )
! spec_rest -> : name ( ',' spec_list )
! spec_list -> ( ',', field_list ) *
! field_list -> expr ( '=' expr + )?
! field_list -> '/' 'name'
! expr -> array
! expr -> limit ( ( ':' | ':<' | '<:' | '<:<' ) limit )?
! array -> '[' ( , expr )* ']'
! limit -> lterm ( 'or' lterm ) *
! lterm -> lneg ( 'and' lneg ) *
! lneg -> 'not'? lfactor
! lfactor -> bterm ( '<' | '<=' | '>' | '>=' | '==' | '/=' bterm ) *
! bterm -> term ( '+' | '-' term ) *
! term -> factor ( '*' | '/' factor )*
! factor -> expon ( '^' expon )*
! expon -> unitless 'unit' ? => n_unit
! expon -> string
! unitless -> primary
! unitless -> '+' primary => n_plus
! unitless -> '-' primary => n_minus
! primary -> 'name' ( '.' 'name' ) ?
! primary -> 'name' '(' expr list ',' ')'
! primary -> 'number'
! primary -> '(' expr ')'

! The notation and method are described in wvs-004.

  use LEXER_CORE, only: PRINT_SOURCE, TOKEN
  use LEXER_M, only: LEXER
  use OUTPUT_M, only: NEWLINE, OUTPUT
  use SYMBOL_TABLE, only: DUMP_1_SYMBOL, DUMP_SYMBOL_CLASS
  use SYMBOL_TYPES ! Everything, especially everything beginning with T_
  use TOGGLES, only: PAR, TOGGLE
  use TREE, only: TREE_BUILDER => BUILD_TREE, N_TREE_STACK, POP, &
                  PUSH_PSEUDO_TERMINAL, STACK_SUBTREE, TX
  use TREE_TYPES   ! Everything, especially everything beginning with N_

  implicit NONE
  private

  public :: CONFIGURATION

  interface Configuration
    module procedure Configuration_I, Configuration_TX
  end interface

  integer, private :: Depth             ! of calls, for tracing

  integer, private :: ERROR             ! 0 => no errors yet
                                        ! 1 => errors but not yet EOF
                                        ! 2 => errors and EOF seen
  type(token), private, save :: NEXT    ! The next token, SAVE only because
                                        ! a component of TOKEN has default
                                        ! initialization, and the standard
                                        ! thereby requires SAVE (go figure)

  ! Tree node to generate depending upon which operator token appears.
  ! We only need values for the ones that generate tree nodes.
  integer, private :: GEN(t_null: t_last_terminal)
  data gen(t_plus)            / n_plus            / ! +
  data gen(t_minus)           / n_minus           / ! -
  data gen(t_star)            / n_mult            / ! *
  data gen(t_slash)           / n_div             / ! /
  data gen(t_backslash)       / n_into            / ! \
  data gen(t_dot)             / n_dot             / ! .
  data gen(t_colon)           / n_colon           / ! :
  data gen(t_colon_less)      / n_colon_less      / ! :<
  data gen(t_less_colon)      / n_less_colon      / ! <:
  data gen(t_less_colon_less) / n_less_colon_less / ! <:<
  data gen(t_equal)           / n_equal           / ! =
  data gen(t_equal_equal)     / n_equal_equal     / ! ==
  data gen(t_not)             / n_not             / ! not
  data gen(t_not_equal)       / n_not_equal       / ! /=
  data gen(t_less)            / n_less            / ! <
  data gen(t_less_eq)         / n_less_eq         / ! <=
  data gen(t_greater)         / n_greater         / ! >
  data gen(t_greater_eq)      / n_greater_eq      / ! >=
  data gen(t_hat)             / n_pow             / ! >
  data gen(t_and)             / n_and             / ! and
  data gen(t_or)              / n_or              / ! or

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
! ----------------------------------------------  CONFIGURATION_I  -----
  subroutine CONFIGURATION_I ( ROOT ) ! cf -> one_cf+ 'EOF'
    integer, intent(out) :: ROOT   ! Root of the abstract syntax tree
    depth = 0 ! in case toggle(par) is set
    error = 0
    call get_token
    do while ( next%class /= t_end_of_input )
      call one_cf
    end do
    call build_tree ( n_cfs, n_tree_stack ) ! collect everything
    call build_tree ( n_eof, 1 )   ! Force n_cfs off the stack
    root = stack_subtree ( 1 )
    if ( error > 0 ) root = -1
  end subroutine CONFIGURATION_I
! ---------------------------------------------  CONFIGURATION_TX  -----
  subroutine CONFIGURATION_TX ( ROOT ) ! cf -> one_cf+ 'EOF'
    type(tx), intent(out) :: ROOT   ! Root of the abstract syntax tree
    depth = 0 ! in case toggle(par) is set
    error = 0
    call get_token
    do while ( next%class /= t_end_of_input )
      call one_cf
    end do
    call build_tree ( n_cfs, n_tree_stack ) ! collect everything
    call build_tree ( n_eof, 1 )   ! Force n_cfs off the stack
    root%i = stack_subtree ( 1 )
    if ( error > 0 ) root%i = -1
  end subroutine CONFIGURATION_TX
! =====     Private Procedures     =====================================
! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( THE_TERMINALS, AFTER )
    integer, intent(in) :: THE_TERMINALS(:)
    integer, intent(in), optional :: AFTER(:)
    integer :: I
    error = max(error,1)
    call output ( '***** Expected ', advance='no' )
    if ( size(the_terminals) > 1 ) then
      call output ( 'one of ', advance='yes' )
      do i = 1, size(the_terminals)
        call dump_symbol_class ( the_terminals(i), advance='yes' )
      end do
      call output ( '*****' )
    else
      call dump_symbol_class ( the_terminals(1), advance='no' )
    end if
    call output ( ' but got ' )
    call dump_1_symbol ( next%string_index, advance='no' )
    call output ( ' at ' )
    call print_source ( next%where, advance='yes' )
    if ( present(after) ) then
      if ( any(next%class == after) ) then
        return ! needs more work -- trying to insert missing symbol
      end if
    end if
    call output ( 'Skipping until expected symbol or end of input.', &
      advance='yes' )
    call output ( 'Processing suppressed.', advance='yes' )
    do while ( all(next%class /= the_terminals) )
      if ( next%class == t_end_of_input ) then
        error = 2
        return
      end if
      call lexer ( next )
    end do
    return
  end subroutine ANNOUNCE_ERROR

! --------------------------------------------------------  ARRAY  -----
  recursive subroutine ARRAY ( HOW_MANY ) ! [ <expr> ( , <expr> )* ]
    integer, intent(inout) :: HOW_MANY  ! Incremented once for each expr
    if ( toggle(par) ) call where ( 'ARRAY' )
    call get_token     ! Consume the left bracket
    if ( next%class == t_right_bracket ) then
      call get_token ! Consume the right bracket
    else
      do
        call expr
        how_many = how_many + 1
        if ( next%class == t_right_bracket ) then
          call get_token ! Consume the right bracket
      exit
        end if
        if ( next%class /= t_comma ) then
          call announce_error ( (/ t_right_bracket, t_comma /) )
          if ( next%class == t_end_of_input .or. &
               next%class == t_end_of_stmt ) &
      exit
        end if
        call get_token   ! Consume the comma
      end do
    end if
    if ( toggle(par) ) call finish ( 'ARRAY' )
  end subroutine ARRAY

! --------------------------------------------------------  BTERM  -----
  recursive subroutine BTERM ! BTERM -> term ( +|- term ) *
    integer :: N
    if ( toggle(par) ) call where ( 'BTERM' )
    call term
    do
      n = next%class
      if ( n == t_plus .or. n == t_minus ) then
        call get_token
        call term
        call build_tree ( gen(n), 2 )
      else
    exit
      end if
    end do
    if ( toggle(par) ) call finish ( 'BTERM' )
  end subroutine BTERM

! ---------------------------------------------------  BUILD_TREE  -----
  subroutine BUILD_TREE ( NEW_NODE, NSONS, DECORATION )
    integer, intent(in) :: NEW_NODE
    integer, intent(in) :: NSONS
    integer, intent(in), optional :: DECORATION

    if ( toggle(par) ) call where ( 'BUILD_TREE' )
    call tree_builder ( new_node, nsons, decoration )
    if ( toggle(par) ) call finish ( 'BUILD_TREE', nsons )
  end subroutine BUILD_TREE

! --------------------------------------------------------  EXPON  -----
  recursive subroutine EXPON  ! expon -> unitless 'unit' ? => n_unit
                              !       -> string
                              ! unitless -> primary
                              !          -> '+' primary => n_plus
                              !          -> '-' primary => n_minus
    integer :: N
    if ( toggle(par) ) call where ( 'EXPON' )
    n = next%class
    if ( n == t_string ) then
      call get_token
    else
      if ( n == t_plus .or. n == t_minus ) then
        call get_token
        call primary
        call build_tree ( gen(n), 1 )
      else
        call primary
      end if
      if ( next%class == t_identifier ) then
        ! First, put unit back on stack (it was popped at end of PRIMARY)
        call push_pseudo_terminal ( next%string_index, next%where )
        call build_tree ( n_unit, 2 )
        call get_token  ! the identifier is used up
      end if
    end if
    if ( toggle(par) ) call finish ( 'EXPON' )
  end subroutine EXPON 

! ---------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR
    ! expr -> array
    ! expr -> limit ( ( ':' | ':<' | '<:' | '<:<' ) limit )?
    integer :: N, NSONS
    if ( toggle(par) ) call where ( 'EXPR' )
    if ( next%class == t_left_bracket ) then
      nsons = 0
      call array ( nsons )
      call build_tree ( n_array, nsons )
    else
      call limit
      n = next%class
      select case ( n )
      case ( t_colon, t_colon_less, t_less_colon, t_less_colon_less )
        call get_token
        call limit
        call build_tree ( gen(n), 2 )
      end select
    end if
    if ( toggle(par) ) call finish ( 'EXPR' )
  end subroutine EXPR

! -------------------------------------------------------  FACTOR  -----
  recursive subroutine FACTOR  ! factor -> expon ( ^ expon )*
    integer :: HOW_MANY       ! sons of the n_pow node
    if ( toggle(par) ) call where ( 'FACTOR' )
    call expon
    if ( next%class == t_hat ) then
      how_many = 1
      do while ( next%class == t_hat )
        call get_token
        call expon
        how_many = how_many + 1
      end do
      call build_tree ( n_pow, how_many )
    end if
    if ( toggle(par) ) call finish ( 'FACTOR' )
  end subroutine FACTOR 

! ---------------------------------------------------  FIELD_LIST  -----
  subroutine FIELD_LIST
    integer :: HOW_MANY       ! sons of the n_asg node
    if ( toggle(par) ) call where ( 'FIELD_LIST' )
    do
      select case ( next%class )
      case ( t_identifier, t_string, & ! field_list -> expr ( '=' expr + )?
             t_number )
        call expr
        if ( next%class == t_equal ) then
          call get_token
          how_many = 1
          call value ( how_many )
          call build_tree ( n_asg, how_many )
        end if
    exit
      case ( t_slash )       ! field_list -> '/' 'name'
        call get_token
        call test_token ( t_identifier )
        call build_tree ( n_set_one, 1 )
    exit
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_identifier, t_string, t_slash /) )
      end select
    end do
    if ( toggle(par) ) call finish ( 'FIELD_LIST' )
  end subroutine FIELD_LIST

! -------------------------------------------------------  FINISH  -----
  subroutine FINISH ( WHAT, NSONS )
  ! for tracing
    character(len=*), intent(in) :: WHAT
    integer, intent(in), optional :: NSONS
    depth = depth - 1
    call output ( repeat('.',depth) // 'Exit  ' // trim(what) )
    if ( present(nsons) ) call output ( nsons, before=' ' )
    call newLine
  end subroutine FINISH

! ----------------------------------------------------  GET_TOKEN  -----
  subroutine GET_TOKEN
  ! Call the lexer to get a token.
  ! If it's a pseudo-terminal token, push it on the tree stack.
    call lexer ( next )
    if ( next%pseudo ) call &
      push_pseudo_terminal ( next%string_index, next%where )
  end subroutine GET_TOKEN

! ------------------------------------------------------  LFACTOR  -----
  recursive subroutine LFACTOR ! lfactor -> bterm ( <|<=|>|>=|==|/= bterm ) *
    integer :: N
    if ( toggle(par) ) call where ( 'LFACTOR' )
    call bterm
    do
      n = next%class
      select case ( n )
      case ( t_less, t_less_eq, t_greater, t_greater_eq, t_equal_equal, &
           & t_not_equal )
        call get_token
        call bterm
        call build_tree ( gen(n), 2 )
      case default
    exit
      end select
    end do
    if ( toggle(par) ) call finish ( 'LFACTOR' )
  end subroutine LFACTOR

! --------------------------------------------------------  LIMIT  -----
  recursive subroutine LIMIT  ! limit -> lterm ( or lterm ) *
    if ( toggle(par) ) call where ( 'LIMIT' )
    call lterm
    do while ( next%class == t_or )
      call get_token
      call lterm
      call build_tree ( n_or, 2 )
    end do
    if ( toggle(par) ) call finish ( 'LIMIT' )
  end subroutine LIMIT

! ---------------------------------------------------------  LNEG  -----
  recursive subroutine LNEG  ! lneg -> not? lfactor
    if ( toggle(par) ) call where ( 'LNEG' )
    if ( next%class == t_not ) then
      call get_token
      call lfactor
      call build_tree ( n_not, 1 )
    else
      call lfactor
    end if
    if ( toggle(par) ) call finish ( 'LNEG' )
  end subroutine LNEG

! --------------------------------------------------------  LTERM  -----
  recursive subroutine LTERM  ! lterm -> lneg ( and lneg ) *
    if ( toggle(par) ) call where ( 'LTERM' )
    call lneg
    do while ( next%class == t_and )
      call get_token
      call lneg
      call build_tree ( n_and, 2 )
    end do
    if ( toggle(par) ) call finish ( 'LTERM' )
  end subroutine LTERM

! -------------------------------------------------------  ONE_CF  -----
  subroutine ONE_CF
    ! one_cf -> ('begin' 'name' 'EOS' spec+ 'end' 'name' )? 'EOS'
    ! one_cf -> #include "file" 'EOS'
    integer :: HOW_MANY       ! How many sons of the generated tree node
    if ( toggle(par) ) call where ( 'ONE_CF' )
o:  do
      select case ( next%class )
      case ( t_end_of_stmt )  ! one_cf -> 'eos'
        call get_token
    exit
      case ( t_begin )        ! one_cf -> 'begin' 'name' 'eos'
                              !           spec+ 'end' 'name' 'eos'
        how_many = 2          ! one for begin, one for end, one for each spec
        call get_token
        call test_token ( t_identifier )
        call test_token ( t_end_of_stmt )
        do while ( next%class /= t_end )
          if ( error > 1 ) exit o
          how_many = how_many + spec()
          if ( next%class == t_end_of_input ) exit ! must have gotten an error
        end do
        call get_token        ! consume the t_end
        call test_token ( t_identifier )
        call test_token ( t_end_of_stmt )
        call build_tree ( n_cf, how_many )
    exit
      case ( t_include )
        call Do_Include
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_begin, t_end_of_stmt /) )
        if ( next%class == t_end_of_input ) exit
      end select
    end do o
    if ( toggle(par) ) call finish ( 'ONE_CF' )
  end subroutine ONE_CF

! ------------------------------------------------------  PRIMARY  -----
  recursive subroutine PRIMARY
    integer :: N
    if ( toggle(par) ) call where ( 'PRIMARY' )
    do
      select case ( next%class )
      case ( t_identifier )   ! primary -> 'name' ( '.' 'name' ) ?
        call get_token
        if ( next%class == t_dot ) then
          call get_token
          call test_token ( t_identifier )
          call build_tree ( n_dot, 2 )
        else if ( next%class == t_left_parenthesis ) then
          ! primary -> 'name' '(' expr list ',' ')'
          n = 1
          call get_token
          do while ( next%class /= t_right_parenthesis )
            call expr
            n = n + 1
            if ( next%class /= t_comma ) &
          exit
            call get_token
          end do
          call build_tree ( n_func_ref, n ) ! Do this BEFORE test_token because
                                            ! test_token will push the unit
                                            ! name if it's there
          call test_token ( t_right_parenthesis )
        end if
    exit
      case ( t_number )       ! primary -> 'number'
        call get_token
    exit
      case ( t_left_parenthesis ) ! primary -> '(' expr ')'
        call get_token
        call expr
        call test_token ( t_right_parenthesis )
    exit
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_identifier, t_number, &
                                 t_string, t_left_parenthesis /), &
                      after = (/ t_end_of_input, t_end_of_stmt,  &
                                 t_right_bracket,  t_comma /) )
        if ( next%class == t_end_of_input .or. &
             next%class == t_end_of_stmt .or. &
             next%class == t_right_bracket .or. &
             next%class == t_comma ) &
    exit
      end select
    end do
    if ( next%class == t_identifier ) call pop ( 1 ) ! delete unit name
    if ( toggle(par) ) call finish ( 'PRIMARY' )
  end subroutine PRIMARY

! ---------------------------------------------------------  SPEC  -----
  integer function SPEC ()
  ! Analyze specifications in a section.
  ! Return how many specifications got generated -- 0 or 1
  ! spec -> 'name' spec_rest 'EOS'
  ! spec -> #include "file" 'EOS'
  ! spec -> 'EOS'
    if ( toggle(par) ) call where ( 'SPEC' )
    do
      select case ( next%class )
      case ( t_end_of_stmt )  ! spec -> 'EOS'
        call get_token
        spec = 0
    exit
      case ( t_identifier )   ! spec -> 'name' spec_rest 'EOS'
        call get_token
        call spec_rest
        call test_token ( t_end_of_stmt )
        spec = 1
    exit
      case ( t_include )
        call do_include
        spec = 0
    exit
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_identifier, t_end_of_stmt /) )
        if ( next%class == t_end_of_input ) exit
      end select
    end do
    if ( toggle(par) ) call finish ( 'SPEC' )
  end function SPEC

! ----------------------------------------------------  SPEC_LIST  -----
  subroutine SPEC_LIST
    integer :: HOW_MANY       ! sons of the n_asg node
    if ( toggle(par) ) call where ( 'SPEC_LIST' )
    how_many = 1
    do while ( next%class == t_comma )  ! spec_list -> ( ',', field_list ) *
      call get_token
      call field_list
      how_many = how_many + 1
    end do
    call build_tree ( n_spec_args, how_many )
    if ( toggle(par) ) call finish ( 'SPEC_LIST' )
  end subroutine SPEC_LIST

! ----------------------------------------------------  SPEC_REST  -----
  subroutine SPEC_REST
    integer :: HOW_MANY       ! sons of the n_equal node
    if ( toggle(par) ) call where ( 'SPEC_REST' )
    do
      select case ( next%class )
      case ( t_end_of_stmt )  ! spec_rest -> \lambda
        call build_tree ( n_spec_args, 1 )
!       call get_token
    exit
      case ( t_equal )        ! spec_rest -> '=' expr +
        call get_token
        how_many = 1
        call value ( how_many )
        call build_tree ( n_equal, how_many )
    exit
      case ( t_comma )        ! spec_rest -> ( ',' spec_list )
        call spec_list
    exit
      case ( t_colon )        ! spec_rest -> : name ( ',' spec_list )
        call get_token
        call test_token ( t_identifier )
        call spec_list
        call build_tree ( n_named, 2 )
    exit
      case default
        if ( error > 1 ) &
    exit
        call announce_error ( (/ t_end_of_stmt, t_equal, t_comma /) )
        if ( next%class == t_end_of_input .or. &
             next%class == t_right_bracket.or. &
             next%class == t_comma ) &
    exit
      end select
    end do
    if ( toggle(par) ) call finish ( 'SPEC_REST' )
  end subroutine SPEC_REST

! ---------------------------------------------------------  TERM  -----
  recursive subroutine TERM   ! term -> factor ( *|/ factor )*
    integer :: N
    if ( toggle(par) ) call where ( 'TERM' )
    call factor
    do
      n = next%class
      select case ( n )
      case ( t_star, t_slash, t_backslash )
        call get_token
        call factor
        call build_tree ( gen(n), 2 )
      case default
    exit
      end select
    end do
    if ( toggle(par) ) call finish ( 'TERM' )
  end subroutine TERM

! ---------------------------------------------------  DO_INCLUDE  -----
  subroutine Do_Include
    use Lexer_Core, only: Where_t
    use String_Table, only: Open_Include
    integer :: File
    type(where_t) :: Where_at ! In case Open_Include needs to print a message
    if ( toggle(par) ) call where ( 'Do_Include' )
    call lexer ( next )
    if ( next%class == t_string ) then
      where_at = next%where
      file = next%string_index
      call lexer ( next )
      if ( next%class /= t_end_of_stmt ) then
        error = max(error,1)
        call announce_error ( (/ t_end_of_stmt /) )
      else
        call open_include ( file, where_at%source, where_at%file )
      end if
    else
      call announce_error ( (/ t_string /) )
    end if
    if ( toggle(par) ) call finish ( 'Do_Include' )
  end subroutine Do_Include

! ---------------------------------------------------  TEST_TOKEN  -----
  subroutine TEST_TOKEN ( THE_TERMINAL )
    integer, intent(in) :: THE_TERMINAL
    if ( next%class == the_terminal ) then
      call get_token
      return
    end if
    if ( error > 1 ) return
    error = max(error,1)
    call announce_error ( (/ the_terminal /) )
  end subroutine TEST_TOKEN

! --------------------------------------------------------  VALUE  -----
  subroutine VALUE ( HOW_MANY )
    integer, intent(inout) :: HOW_MANY  ! Incremented once for each expr
    select case ( next%class )
    case ( t_left_bracket )
      call array ( how_many )
    case ( t_comma, t_end_of_stmt )
    case default
      call expr
      how_many = how_many + 1
    end select
  end subroutine VALUE

! --------------------------------------------------------  WHERE  -----
  subroutine WHERE ( WHAT )
    character(len=*), intent(in) :: WHAT
    call output ( repeat('.',depth) // 'Enter ' // trim(what) )
    call output ( ' at ' )
    call print_source ( next%where, advance='yes' )
    depth = depth + 1
  end subroutine WHERE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PARSER

! $Log$
! Revision 2.28  2013/09/30 23:03:04  vsnyder
! Add TX type for tree index and generics to use it
!
! Revision 2.27  2013/09/25 01:02:41  vsnyder
! Add include files
!
! Revision 2.26  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.25  2012/05/05 00:11:51  vsnyder
! Add support for 'not' operator
!
! Revision 2.24  2012/05/01 22:11:27  vsnyder
! Simplify node generation
!
! Revision 2.23  2012/05/01 22:10:26  vsnyder
! Add TrueList subroutine
!
! Revision 2.22  2012/04/24 20:37:24  vsnyder
! Include the complete grammar at the top
!
! Revision 2.21  2011/04/19 01:59:43  vsnyder
! Support == and /= relational operators too
!
! Revision 2.20  2011/04/18 19:33:26  vsnyder
! Add support for relational operators and boolean-valued expressions
!
! Revision 2.19  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.18  2008/09/04 20:02:20  vsnyder
! Add PRINT statement in not_used_here + Cannonball polishing
!
! Revision 2.17  2008/09/04 00:46:17  vsnyder
! Reverse precedence of unary +/- and units
!
! Revision 2.16  2006/03/22 03:04:00  vsnyder
! Allow empty spec field values
!
! Revision 2.15  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2004/05/28 23:13:12  vsnyder
! Add power (^) operator, units coercion for (expr) and function result
!
! Revision 2.13  2004/01/20 19:43:33  vsnyder
! Cosmetic changes
!
! Revision 2.12  2004/01/17 03:04:48  vsnyder
! Provide for functions in expressions
!
! Revision 2.11  2004/01/16 23:49:32  vsnyder
! Add backslash for 'into' operator
!
! Revision 2.10  2003/01/29 03:14:14  vsnyder
! Print the line and column in ENTER debugging output
!
! Revision 2.9  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2001/11/28 03:05:54  vsnyder
! Implement arrays of arrays
!
! Revision 2.7  2001/11/27 00:54:37  vsnyder
! Implement (partially) open ranges
!
! Revision 2.6  2001/07/20 20:18:05  vsnyder
! Improve error recovery in a few cases -- more work probably needed
!
! Revision 2.5  2001/02/28 02:37:11  vsnyder
! Allow specification with no arguments to have a label
!
! Revision 2.4  2000/11/30 00:23:10  vsnyder
! Implement [] syntax for arrays
!
! Revision 2.3  2000/11/15 22:01:18  vsnyder
! Allow specification with no arguments.
!
! Revision 2.2  2000/11/15 21:15:27  vsnyder
! Corrected a loop on a bad primary; correct processing of spec with no fields
!
! Revision 2.1  2000/10/11 18:02:50  vsnyder
! Move from lib/cf_parser to lib; remove unused variables; add copyright
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
