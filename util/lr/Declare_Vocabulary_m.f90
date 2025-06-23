! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Declare_Vocabulary_m

  ! Declare symbols in the vocabulary, according to whether they appear
  ! as terminals, nonterminals, vocabulary names, or actions.

  implicit NONE
  private

  public :: Declare_Vocabulary

  integer, public :: Total   ! Number of productions + sum of RHS lengths

  integer, public:: EOG      ! String index for <EOG> terminal
  integer, public:: GOAL     ! String index for <GOAL> nonterminal
  integer, public:: SOG      ! String index for <SOG> terminal

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Declare_Vocabulary ( Root, Error )

    use Declaration_Table, only: Action, Declare, Decls, Get_Decl, &
      & Nonterminal, Redeclare, Terminal, Vocabulary
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Symbol_Table, only: Enter_Terminal
    use Symbol_Types, only: T_Identifier
    use Tables, only: NUMPRD
    use Toggles, only: Switches
    use Tree, only: Decorate, Node_Id, NSons, Subtree, Sub_Rosa
    use Tree_Types, only: N_Equal, N_Generates, N_Production, N_Question, N_Rhs

    integer, intent(in) :: Root ! Of the abstract syntax tree
    logical, intent(inout) :: Error

    ! Error message codes
    integer, parameter :: Already_Nonterminal = 1
    integer, parameter :: Terminal_And_Nonterminal = Already_Nonterminal + 1

    integer :: Act              ! Tree index of action node
    type(decls) :: Decl         ! of a symbol
    integer :: GSON             ! Son of Son
    integer :: I, J             ! Loop indices
    integer :: Son              ! of Root
    integer :: Symbol           ! String index of something
    logical :: Watch

    ! Do we watch it work?
    watch = index(switches,'dvoc') /= 0

    eog = enter_terminal ( "EOG", t_identifier )
    goal = enter_terminal ( "<GOAL>", t_identifier )
    sog = enter_terminal ( "SOG", t_identifier )

    call declare ( eog, terminal, 0 )
    call declare ( goal, nonterminal, 0 )
    call declare ( sog, terminal, 0 )

    numprd = 1 ! Includes <GOAL> -> <SOG> start <EOG>
    total = 4  ! Includes <GOAL> -> <SOG> start <EOG>

    ! Declare the vocabulary symbols and the nonterminals
    if ( watch ) &
      & call output ( 'Declare vocabulary symbols and nonterminals', advance='yes' )
    do i = 1, nsons(root)
      son = subtree(i,root)
      select case ( node_id(son) )
      case ( n_equal )
        !              String                    Type      Tree
        call redeclare ( sub_rosa(subtree(1,son)), terminal, subtree(1,son), &
        !   Value
          & sub_rosa(subtree(2,son)) )
        !              String                    Type      Tree
        call redeclare ( sub_rosa(subtree(2,son)), vocabulary, subtree(1,son), &
        !   Value
          & sub_rosa(subtree(1,son)) )
        call decorate ( subtree(1,son), sub_rosa(subtree(2,son)) )
      case ( n_production )
        gson = subtree(1,son)
        symbol = sub_rosa(gson)
        decl = get_decl(symbol,nonterminal)
        if ( decl%type == nonterminal ) then
          call announce_error ( gson, already_nonterminal, symbol )
        else
          if ( watch ) then
            call display_string ( sub_rosa(gson), before='Declare ' )
            call output ( ' as nonterminal.', advance='yes' )
          end if
          call declare ( sub_rosa(gson), nonterminal, gson )
        end if
        numprd = numprd + nsons(son) - 1
      end select
    end do

    ! Declare the terminal symbols and action symbols
    if ( watch ) &
      & call output ( 'Declare terminal and action symbols', advance='yes' )
    do i = 1, nsons(root)
      son = subtree(i,root)
      select case ( node_id(son) )
      case ( n_equal )
        gson = subtree(1,son)
        symbol = sub_rosa(gson)
        decl = get_decl(symbol,nonterminal)
        if ( decl%type == nonterminal ) then
          call announce_error ( gson, terminal_and_nonterminal, symbol )
        end if
      case ( n_production )
        total = total + nsons(son) - 1 ! One LHS for each RHS
        do j = 2, nsons(son)
          gson = subtree(j,son)
          select case ( node_id(gson) )
          case ( n_rhs )
            call declare_rhs ( gson )
          case ( n_generates )
            call declare_rhs ( subtree(1,gson) )
            act = subtree(2,gson)
            if ( node_id(act) == n_question ) then
              call declare ( sub_rosa(subtree(1,act)), action, subtree(1,gson) )
            else
              call declare ( sub_rosa(act), action, subtree(1,gson) )
            end if
          end select
        end do
      end select
    end do

    if ( watch ) then
      call output ( numprd, 5 );  call output ( ' Productions', advance='yes' )
      call output ( total, before='Total size of grammar = ', advance='yes' )
    end if

  contains

    subroutine Announce_Error ( Where, What, Symbol )

      use Lexer_Core, only: Print_Source
      use Output_m, only: Output
      use String_Table, only: Display_String
      use Tree, only: Source_Ref
      integer, intent(in) :: Where            ! Tree index
      integer, intent(in) :: What             ! Error message code
      integer, intent(in), optional :: Symbol ! String index
      call print_source ( source_ref(where), before='Error: at ' )
      select case ( what )
      case ( already_nonterminal )
        call display_string ( symbol, before=' the symbol "' )
        call output ( '" has already appeared as a nonterminal', advance='yes' )
      case ( terminal_and_nonterminal )
        call display_string ( symbol, before=' the symbol "' )
        call output ( '" appears both in a vocabulary definition and as a nonterminal', &
          & advance='yes' )
      end select
      error = .true.

    end subroutine Announce_Error

    subroutine Declare_RHS ( Root )
    ! Declare undeclared symbols that are sons of Root, which is an RHS node
    ! If a symbol is not declared assume it's terminal.
      use Declaration_Table, only: Declare, Decls, Empty, Get_Decl, Null_Decl
      integer, intent(in) :: Root
      type(decls) :: Decl ! Of a symbol
      integer :: I        ! Loop index
      integer :: Son      ! Of root
      integer :: Symbol   ! String index of son of root

      total = total + nsons(root)
      do i = 1, nsons(root)
        son = subtree(i,root)
        symbol = sub_rosa(son)
        decl = get_decl(symbol,terminal)
        select case ( decl%type )
        case ( null_decl, empty )
          if ( watch ) then
            call display_string ( symbol, before='Declare ' )
            call output ( ' as terminal.', advance='yes' )
          end if
          call redeclare ( symbol, terminal, son, 0 )
        end select
      end do

    end subroutine Declare_RHS

  end subroutine Declare_Vocabulary

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Declare_Vocabulary_m
