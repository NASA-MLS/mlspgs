! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PARSER

! Parser for L2CF.  The grammar is described in wvs-004.

  use LEXER_CORE, only: PRINT_SOURCE, TOKEN
  use LEXER_M, only: LEXER
  use OUTPUT_M, only: OUTPUT
  use SYMBOL_TABLE, only: DUMP_1_SYMBOL, DUMP_SYMBOL_CLASS
  use SYMBOL_TYPES ! Everything, especially everything beginning with T_
  use TOGGLES, only: PAR, TOGGLE
  use TREE, only: BUILD_TREE, N_TREE_STACK, &
                  PUSH_PSEUDO_TERMINAL, STACK_SUBTREE
  use TREE_TYPES   ! Everything, especially everything beginning with N_
  implicit NONE
  private
  public :: CONFIGURATION

  integer, private :: ERROR             ! 0 => no errors yet
                                        ! 1 => errors but not yet EOF
                                        ! 2 => errors and EOF seen
  type(token), private :: NEXT          ! The next token

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
! ------------------------------------------------  CONFIGURATION  -----
  subroutine CONFIGURATION ( ROOT )
    integer, intent(out) :: ROOT   ! Root of the abstract syntax tree
    error = 0
    call get_token
      do
      if ( next%class == t_end_of_input ) exit
      call one_cf
    end do
    call build_tree ( n_cfs, n_tree_stack ) ! collect everything
    call build_tree ( n_eof, 1 )   ! Force n_cfs off the stack
    root = stack_subtree ( 1 )
    if ( error > 0 ) root = -1
  end subroutine CONFIGURATION
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
    call print_source ( next%source, advance='yes' )
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
! -------------------------------------------------------    EXPR  -----
  subroutine EXPR
    if ( toggle(par) ) call output ( 'Enter EXPR', advance='yes' )
    call limit
    if ( next%class == t_colon ) then
      call get_token
      call limit
      call build_tree ( n_colon, 2 )
    end if
    if ( toggle(par) ) call output ( 'Exit  EXPR', advance='yes' )
  end subroutine EXPR
! -------------------------------------------------------  FACTOR  -----
  recursive subroutine FACTOR  ! factor -> ( +|- )? primary
    if ( toggle(par) ) call output ( 'Enter FACTOR', advance='yes' )
    if ( next%class == t_plus ) then
      call get_token
      call primary
      call build_tree ( n_plus, 1 )
    else if ( next%class == t_minus ) then
      call get_token
      call primary
      call build_tree ( n_minus, 1 )
    else
      call primary
    end if
    if ( toggle(par) ) call output ( 'Exit  FACTOR', advance='yes' )
  end subroutine FACTOR 
! ----------------------------------------------------  FIELD_LIST  -----
  subroutine FIELD_LIST
    integer :: HOW_MANY       ! sons of the n_asg node
    if ( toggle(par) ) call output ( 'Enter FIELD_LIST', advance='yes' )
    do
      select case ( next%class )
      case ( t_identifier, t_string, & ! field_list -> expr ( '=' expr + )?
             t_number )
        call expr
        if ( next%class == t_equal ) then
          call get_token
          how_many = 2
          call expr
          do while ( any(next%class == &
              (/ t_identifier, t_number, t_string, &
                 t_left_parenthesis, t_dot /) ) )
            call expr
            how_many = how_many + 1
          end do
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
    if ( toggle(par) ) call output ( 'Exit  FIELD_LIST', advance='yes' )
  end subroutine FIELD_LIST
! ----------------------------------------------------  GET_TOKEN  -----
  subroutine GET_TOKEN
  ! Call the lexer to get a token.
  ! If it's a pseudo-terminal token, push it on the tree stack.
    call lexer ( next )
    if ( next%pseudo ) call &
      push_pseudo_terminal ( next%string_index, next%source )
  end subroutine GET_TOKEN
! --------------------------------------------------------  LFACTOR  -----
  recursive subroutine LFACTOR ! lfactor -> term ( +|- term ) *
    if ( toggle(par) ) call output ( 'Enter LFACTOR', advance='yes' )
    call term
    do
      if ( next%class == t_plus ) then
        call get_token
        call term
        call build_tree ( n_plus, 2 )
      else if ( next%class == t_minus ) then
        call get_token
        call term
        call build_tree ( n_minus, 2 )
      else
    exit
      end if
    end do
    if ( toggle(par) ) call output ( 'Exit  LFACTOR', advance='yes' )
  end subroutine LFACTOR
! --------------------------------------------------------  LIMIT  -----
  recursive subroutine LIMIT  ! limit -> lterm ( or lterm ) *
    if ( toggle(par) ) call output ( 'Enter LIMIT', advance='yes' )
    call lterm
    do while ( next%class == t_or )
      call get_token
      call lterm
      call build_tree ( n_or, 2 )
    end do
    if ( toggle(par) ) call output ( 'Exit  LIMIT', advance='yes' )
  end subroutine LIMIT
! --------------------------------------------------------  LTERM  -----
  recursive subroutine LTERM  ! lterm -> lfactor ( and lfactor ) *
    if ( toggle(par) ) call output ( 'Enter LTERM', advance='yes' )
    call lfactor
    do while ( next%class == t_and )
      call get_token
      call lfactor
      call build_tree ( n_and, 2 )
    end do
    if ( toggle(par) ) call output ( 'Exit  LTERM', advance='yes' )
  end subroutine LTERM
! -------------------------------------------------------  ONE_CF  -----
  subroutine ONE_CF
    integer :: HOW_MANY       ! How many sons of the generated tree node
    if ( toggle(par) ) call output ( 'Enter ONE_CF', advance='yes' )
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
          how_many = how_many +  spec()
          if ( next%class == t_end_of_input ) exit ! must have gotten an error
        end do
        call get_token        ! consume the t_end
        call test_token ( t_identifier )
        call test_token ( t_end_of_stmt )
        call build_tree ( n_cf, how_many )
    exit
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_begin, t_end_of_stmt /) )
        if ( next%class == t_end_of_input ) exit
      end select
    end do o
    if ( toggle(par) ) call output ( 'Exit  ONE_CF', advance='yes' )
  end subroutine ONE_CF
! ------------------------------------------------------  PRIMARY  -----
  recursive subroutine PRIMARY
    if ( toggle(par) ) call output ( 'Enter PRIMARY', advance='yes' )
    do
      select case ( next%class )
      case ( t_identifier )   ! primary -> 'name' ( '.' 'name' ) ?
        call get_token
        if ( next%class == t_dot ) then
          call get_token
          call test_token ( t_identifier )
          call build_tree ( n_dot, 2 )
        end if
    exit
      case ( t_number )       ! primary -> 'number' 'unit' ?
        call get_token
        if ( next%class == t_identifier ) then
          call build_tree ( n_unit, 2 )
          call get_token      ! the identifier is used up
        end if
    exit
      case ( t_string )       ! primary -> 'string'
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
                                 t_string, t_left_parenthesis /) )
      end select
    end do
    if ( toggle(par) ) call output ( 'Exit  PRIMARY', advance='yes' )
  end subroutine PRIMARY
! ---------------------------------------------------------  SPEC  -----
  integer function SPEC ()
  ! Analyze specifications in a section.
  ! Return how many specifications got generated -- 0 or 1
    if ( toggle(par) ) call output ( 'Enter SPEC', advance='yes' )
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
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_identifier, t_end_of_stmt /) )
        if ( next%class == t_end_of_input ) exit
      end select
    end do
    if ( toggle(par) ) call output ( 'Exit  SPEC', advance='yes' )
  end function SPEC
! ----------------------------------------------------  SPEC_LIST  -----
  subroutine SPEC_LIST
    integer :: HOW_MANY       ! sons of the n_asg node
    if ( toggle(par) ) call output ( 'Enter SPEC_LIST', advance='yes' )
    how_many = 2
    call get_token
    call field_list
    do while ( next%class == t_comma )
      call get_token
      call field_list
      how_many = how_many + 1
    end do
    call build_tree ( n_spec_args, how_many )
    if ( toggle(par) ) call output ( 'Exit  SPEC_LIST', advance='yes' )
  end subroutine SPEC_LIST
! ----------------------------------------------------  SPEC_REST  -----
  subroutine SPEC_REST
    integer :: HOW_MANY       ! sons of the n_equal node
    if ( toggle(par) ) call output ( 'Enter SPEC_REST', advance='yes' )
    do
      select case ( next%class )
      case ( t_end_of_stmt )  ! spec_rest -> \lambda
        call build_tree ( n_spec_args, 1 )
!       call get_token
    exit
      case ( t_equal )        ! spec_rest -> '=' expr +
        call get_token
        call expr
        how_many = 2
        do while ( any(next%class == &
            (/ t_identifier, t_number, t_string, &
               t_left_parenthesis, t_dot, t_plus, t_minus /) ) )
          call expr
          how_many = how_many + 1
        end do
        call build_tree ( n_equal, how_many )
    exit
      case ( t_comma )        ! spec_rest -> ( ',' spec_list ) +
        call spec_list
    exit
      case ( t_colon )        ! spec_rest -> : name ( ',' spec_list ) +
        call get_token
        call test_token ( t_identifier )
        if ( next%class == t_comma ) then
          call spec_list
          call build_tree ( n_named, 2 )
    exit
        end if
        call announce_error ( (/ t_comma /) )
      case default
        if ( error > 1 ) exit
        call announce_error ( (/ t_end_of_stmt, t_equal, t_comma /) )
        if ( next%class == t_end_of_input ) exit
      end select
    end do
    if ( toggle(par) ) call output ( 'Exit  SPEC_REST', advance='yes' )
  end subroutine SPEC_REST
! ---------------------------------------------------------  TERM  -----
  recursive subroutine TERM   ! term -> primary ( *|/ primary )*
    if ( toggle(par) ) call output ( 'Enter TERM', advance='yes' )
    call factor
    do
      if ( next%class == t_star ) then
        call get_token
        call factor
        call build_tree ( n_mult, 2 )
      else if ( next%class == t_slash ) then
        call get_token
        call factor
        call build_tree ( n_div, 2 )
      else
    exit
      end if
    end do
    if ( toggle(par) ) call output ( 'Exit  TERM', advance='yes' )
  end subroutine TERM
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

end module PARSER

! $Log$
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
