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

! Parser for MLS CF.  The grammar is processed by the LR processor.

  implicit NONE
  private

  public :: Clean_Up_Parser, Configuration, LR_Parser, LR_Parser_TX

  interface Configuration
    module procedure LR_Parser, LR_Parser_TX
  end interface

  ! Derived from lists in parser tables, to make the parser run faster:

  integer, save, allocatable :: LR_0(:) ! If LR_0(I) /= 0, state I is an LR(0)
                             ! state, i.e., one reduction and no productions.
                             ! LR_0(I) is the production number to reduce.
  integer, save, allocatable :: Work(:,:)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine LR_Parser_TX ( Root, T )
    use Parser_Table_m, only: Parser_Table_t
    use Tree, only: TX
    type(tx), intent(out) :: ROOT   ! Root of the abstract syntax tree
    type(parser_table_t), intent(in) :: T
    call lr_parser ( root%i, T )
  end subroutine LR_Parser_TX

  subroutine LR_Parser ( Root, T )

    use LEXER_CORE, only: PRINT_SOURCE, TOKEN
    use LEXER_M, only: LEXER
    use OUTPUT_M, only: Blanks, NEWLINE, OUTPUT
    use Parser_Table_m, only: Parser_Table_t
    use SYMBOL_TABLE, only: DUMP_1_SYMBOL, DUMP_SYMBOL_CLASS
    use SYMBOL_TYPES, only: T_End_Of_Input, T_Last_Terminal, T_Null
    use TOGGLES, only: LEVELS, PAR, TOGGLE
    use TREE, only: BUILD_TREE, Dump_Top_Stack, N_TREE_STACK,&
                    PUSH_PSEUDO_TERMINAL,  STACK_SUBTREE
    use TREE_TYPES, only: N_CFS, N_Null

    integer, intent(out) :: Root   ! Index in tree of root of the parsed tree
    type(parser_table_t), intent(in) :: T

    type :: Stack_Frame
      integer :: State         ! The parser state
      integer :: Tree_Index    ! Tree stack pointer of state's t%lhs
    end type Stack_Frame

    logical :: Need            ! Need to get a token from the lexer
    integer :: Prev_Error

    logical, save :: Machine = .false. ! Dump the grammar and parsing machine
    integer :: Map(t_null:t_last_terminal) ! Map token classes to parser
                               ! vocabulary indices
    integer :: Newsta, Nowsta  ! Next state, Current state
    logical, save :: Show_stack = .false. ! stack too
    logical, save :: Show_state = .false. ! automaton state too
    logical, save :: Show_work = .false. ! dump work table
    integer :: SP              ! stack pointer
    type(stack_frame), allocatable :: Stack(:)
    type(token) :: The_Token   ! from the lexer
    integer :: Voc             ! vocabulary index of token or t%lhs

    ! Put names of symbols not defined by the grammar into the string
    ! table
    if ( .not. allocated(work) ) then

      ! Compute debug flags
      show_stack = mod(levels(par),2) /= 0
      show_state = mod(levels(par)/2,2) /= 0
      machine = mod(levels(par)/4,2) /= 0
      show_work = mod(levels(par)/8,2) /= 0

      ! Convert the lists to search for transitions and reductions to the
      ! array WORK(t%nstate:t%nvoc).
      ! If Work(i,j) > 0 the next state is Work(i,j)
      ! If Work(i,j) < 0 reduce production -Work(i,j)
      ! If Work(i,j) == 0 a syntax error has occurred
      allocate ( work(t%nstate,t%nvoc), LR_0(t%nstate) )
      call make_work ( t%ftrn, t%tran, t%ent, t%fred, t%nset, t%lset, t%ls, &
                     & t%prod, t%val, Work, LR_0, map )
    end if

    prev_error = 0              ! Line number of previous error
    sp = 0
    newsta = abs(t%tran(1))     ! The LR generator likes to put the <SOG>
                                ! symbol before the goal symbol.  The lexer
                                ! doesn't produce this symbol before the first
                                ! input.  Putting the parser in the state
                                ! adjacent to state 1 skips that symbol.
    voc = 0                     ! Don't display symbol in Transition routine

    call transition ( .false. ) ! push Newsta onto the parser stack and
                                ! make it Nowsta.

    need = .true.               ! Indicate a token is necessary

    do while ( nowsta /= t%ifinal )
      if ( lr_0(nowsta) > 0 ) then
        ! The current state is an LR(0) state, meaning it does one reduction
        ! and no transitions.  The production in lr_0(nowsta) can be
        ! reduced without needing to go to the lexer to get a token.
        call reduce_production ( lr_0(nowsta) )
      else
        if ( need ) then
          call lexer ( the_token )
          need = .false.
        end if
        ! Decide what to do
        if ( map(the_token%class) <= 0 .or. map(the_token%class) > t%nterms ) then
          voc = 0 ! Syntax error, token isn't even in terminal alphabet
        else
          voc = work(nowsta,map(the_token%class))
        end if
        if ( voc > 0 ) then ! a transition to state voc.
          newsta = voc
          if ( the_token%pseudo ) then
            call push_pseudo_terminal ( the_token%string_index, &
                                      & the_token%where, class=the_token%class )
            if ( show_stack ) then
              call output ( 'Push pseudo terminal ' )
              call dump_1_symbol ( the_token%string_index )
              call print_source ( the_token%where, before=' at ' )
              call output ( n_tree_stack, before=' leaving n_tree_stack = ', &
                & advance='yes' )
            end if
          end if
          call transition ( .false. )
          need = .true.        ! Consume the token
        else if ( voc < 0 ) then ! a reduction by production -voc
          call reduce_production ( -voc )
        else ! Syntax error
         call no_transition
        end if
      end if
    end do

    call build_tree ( n_cfs, n_tree_stack )  ! glue together everything on the

    if ( prev_error > 0 ) then
      call output ( prev_error, before='Last error at line ', advance='yes' )
      root = -1
    else
      if ( toggle(par) ) then
      call output ( 'Finished parsing and built ' )
      call dump_top_stack ( 0, advance='yes' )
      end if
      call build_tree ( n_null, 1 )          ! pop the n_cfs, leave the n_null
      root = stack_subtree(1)
    end if

  contains

! -------------------------------------------  Catastrophic_Error  -----
    subroutine Catastrophic_Error ( text, The_Token, State, Nonterminal )
      use Machine, only: Crash_Burn, Nevercrash
      use String_Table, only: Display_String
      character(len=*), intent(in):: text
      type(token), intent(in), optional  :: The_Token
      integer, intent(in), optional :: State, Nonterminal
      if ( present(the_token) ) then
        call print_source ( the_token%where, &
          & before='***** Catastrophic Error ***** ', advance='yes' )
      else
        call output ( '***** Catastrophic Error *****', advance='yes' )
      end if
      if ( present(state) ) &
        & call output ( state, before='The current state is ', advance='yes' )
      if ( present(nonterminal) ) &
        & call output ( nonterminal, before='The nonterminal symbol is ', &
          & advance='yes' )
      if ( present(the_token) ) then
        call output ( 'Next token is ' )
        call dump_symbol_class ( the_token%class )
        if ( the_token%pseudo .and. the_token%string_index > 0 ) &
          & call display_string ( the_token%string_index, before=' ' )
        call newLine
      end if
      call output( 'About to crash and burn in parser', advance='yes' )
      Nevercrash = .false.
      call Crash_Burn
      ! call MLSMessage ( MLSMSG_Crash, moduleName, trim(text) )
      stop 666
    end subroutine Catastrophic_Error

! ---------------------------------------------------  Do_Include  -----
    subroutine Do_Include
      ! We don't want "string" in the tree, so pop it off the stack.
      use String_Table, only: Display_String, Open_Include
      use Tree, only: Pop, Stack_File, Stack_Source_Ref, Stack_Sub_Rosa
      call open_include ( Stack_Sub_Rosa(), Stack_Source_Ref(), &
                        & Stack_File() )
      if ( toggle(par) ) &
        & call display_string ( Stack_Sub_Rosa(), &
          & before='Opened include file ', advance='yes' )
      call pop ( 1 )
   end subroutine Do_Include

! -----------------------------------------------  Error_Walkback  -----
    subroutine Error_Walkback ( Line )
      integer, intent(in) :: Line      ! Line number of error
      if ( prev_error > 0 ) call output ( prev_error, &
        & before='Previous error at line ', advance='yes' )
      prev_error = line
    end subroutine Error_Walkback

! ----------------------------------------------------  Make_Work  -----
    subroutine Make_Work ( ftrn, tran, ent, fred, nset, lset, ls, prod, val, &
                         & Work, LR_0, Map )

    ! Make the Work, LR_0 and Map arrays, to avoid searching ftrn, tran,
    ! fred, nset, lset, and ls to determine whether to reduce a production
    ! or make a transition.

    ! Positive elements of WORK are new state numbers and indicate transitions.
    ! Negative elements of WORK are indices of productions to reduce.
    ! Zero elements of WORK indicate syntax errors.

    ! Columns 1:nterms are indexed by terminal symbol vocabulary indices. 
    ! Columns nterms+1:nvoc are indexed by nonterminal symbol vocabulary
    ! indices, never indicate reductions, and errors are not possible.

      integer, intent(in) :: ftrn(0:), tran(:), ent(:)
      integer, intent(in) :: fred(0:), nset(:), lset(0:), ls(:)
      integer, intent(in) :: prod(:), val(:)
      integer, intent(out) :: WORK(:,:) ! ( nstate : nvoc )
      integer, intent(out) :: LR_0(:)   ! ( nstate ), "state is an LR(0) state"
      integer, intent(out) :: Map(t_null:) ! Inverse of val

      ! ftrn(0:nstate): Indexed by state number.  First (actually last)
      ! transitions.  Transitions for state S are in tran(ftrn(S-1)+1:ftrn(S)).

      ! tran(1:ntrans): New state indices, indexed by ftrn, q.v.

      ! ent(1:nstate): Indexed by state number.  Index of vocabulary symbol
      ! used to enter a state.

      ! fred(0:nstate): Indexed by state number.  First (actually last)
      ! reductions.  Lookahead sets for state S are nset(fred(S-1)+1:fred(S)).

      ! nset(1:nlooks): Lookahead set indices, indexed by fred, q.v.

      ! lset(0:ncsets): Last elements of lookahead sets, indexed by nset, q.v.
      ! For lookahead set L, elements are ls(lset(L-1)+1:lset(L)).

      ! ls(1:lsets): Elements of lookahead sets.

      ! prod(1:nlooks): Production numbers to reduce if lookahead is in the
      ! set, indexed by fred, q.v.

      logical :: Error
      integer :: I, J, L

      error = .false. ! Assume OK
      work = 0

      ! Fill in the transitions
      do i = 1, ubound(ftrn,1)
        work(i,ent(tran(ftrn(i-1)+1:ftrn(i)))) = tran(ftrn(i-1)+1:ftrn(i))
      end do

      ! Fill in the reductions
      do i = 1, ubound(fred,1)
        do j = fred(i-1)+1, fred(i)
          do l = lset(nset(j)-1)+1, lset(nset(j))
            if ( work(i,ls(l)) > 0 ) then
              write ( *, '(3(a,i0))' ) &
                & 'Grammar is not LR.  Intersection in state ', i, &
                & ' with transition to state ', work(i,ls(l)), &
                & ' on vocabulary symbol with internal index ', ls(l)
                error = .true.
            else if ( work(i,ls(l)) < 0 ) then
              write ( *, '(3(a,i0))' ) &
                & 'Grammar is not LR.  Intersection in state ', i, &
                & ' with reduction of production ', -work(i,ls(l)), &
                & ' on vocabulary symbol with internal index ', ls(l)
                error = .true.
            else
              work(i,ls(l)) = -prod(j)
            end if
          end do
        end do
      end do

      do i = 1, ubound(work,1)
        if ( ftrn(i) == ftrn(i-1) .and. &    ! No transitions
           & fred(i) == fred(i-1)+1 ) then   ! One reduction
          lr_0(i) = prod(fred(i))
        else
          lr_0(i) = 0
        end if
      end do

      do i = 1, ubound(val,1)
        if ( val(i) > 0 .and. val(i) < ubound(map,1) ) map(val(i)) = i
      end do

      if ( show_work ) then
        call output ( "Parser's WORK array:", advance='yes' )
        do i = 1, size(work,2), 20
          call output ( '      ' )
          do l = i, min(i+19,size(work,2))
            call output ( l, format='(i6)' )
          end do
          call newline
          do j = 1, size(work,1)
            call output ( j, format='(i5,":")' )
            do l = i, min(i+19,size(work,2))
              call output ( work(j,l), format='(i6)' )
            end do
            call newline
          end do
        end do
      end if

      if ( error ) then
        write ( *, '(a)' ) 'Error: Make_Work_m.Make_Work%E%- Grammar is not LR(1)'
        stop 1
      end if

      if ( machine ) then
        call print_grammar
        call output ( 'The parsing automaton (sans context sets for reductions):', &
          & advance='yes' )
        do i = 1, t%nstate
          call newLine
          call print_state ( i, number=.true. )
        end do
        call output ( 'End of the parsing automaton', advance='yes' )
      end if

    end subroutine Make_Work

! ------------------------------------------------  No_Transition  -----
    subroutine No_Transition ( Expected )
      ! This is a REALLY CRUDE error recovery routine.
      integer, intent(in), optional :: Expected
      integer :: I
      call output ( the_token%where%source / 256,5 )
      call output ( ': *** Error ***  Unexpected Input ')
      call dump_1_symbol ( the_token%string_index )
      if ( the_token%class /= t_end_of_input ) then
        call print_source ( the_token%where, before=' at ' )
        call output (' in parser state ' )
        call output ( nowsta, advance='yes' )
        call print_state ( nowsta )
        if ( present(expected) ) then
          call output ( 'Expected ', advance='no' )
          call dump_symbol_class ( expected, advance='yes' )
        else
          call output ( 'Expected one of: ', advance='yes' )
          do i = 1, t%nterms
            if ( work(nowsta,i) > 0 ) &
              & call dump_symbol_class ( map(i), advance='yes' )
          end do
        end if
        call lexer( the_token )
        call error_walkback ( the_token%where%source / 256 )
      else
        call output (' in parser state ' )
        call output ( nowsta, advance='yes' )
        call print_state ( nowsta )
        if ( present(expected) ) then
          call output ( 'Expected ', advance='no' )
          call dump_symbol_class ( expected, advance='yes' )
        else
          call output ( 'Expected one of: ', advance='yes' )
          do i = 1, t%nterms
            if ( work(nowsta,i) > 0 ) &
              & call dump_symbol_class ( map(i), advance='yes' )
          end do
        end if
        call error_walkback ( the_token%where%source / 256 )
        sp = sp - 1
        if ( sp <= 0 ) call catastrophic_error ( 'Parser stack underflow' )
        nowsta = stack(sp)%state
      end if
    end subroutine No_Transition

! -------------------------------------------------  Print_Action  -----
    subroutine Print_Action ( P )
    ! Print the action for production P
      use String_Table, only: Display_String
      use Tree, only: Tree_Text
      integer, intent(in) :: P ! The production number
      integer :: A
      a = t%act(p)
      if ( a /= 0 ) then
        select case ( mod(a,10) )
        case ( 1 )
          call display_string ( tree_text(a / 10), before=' => ' )
        case ( 2 )
          call display_string ( tree_text(a / 10), before=' => ' )
          call output ( ' ?' )
        case default
          call output ( a, before=' => ' )
        end select
      end if
    end subroutine Print_Action

! ------------------------------------------------  Print_Grammar  -----
    subroutine Print_Grammar
    ! Print the entire grammar, with actions.
      use String_Table, only: Get_String, String_Length
      integer :: J, K, LHS, P, PrevLHS
      character(120) :: Work ! To strip underscore off nonterminal names
      call output ( 'The grammar:', advance='yes' )
      prevLHS = -1
      do p = 1, t%numprd
        j = t%prod_ind(p)
        lhs = t%productions(j)
        if ( lhs /= prevLHS .and. p /= 1 ) call newLine
        call output ( p, 3 )
        call get_string ( t%text(lhs), work )
        k = merge(2,1,work(1:1)=='_')
        if ( lhs /= prevLHS ) then
          call output ( ' ' // trim(work(k:)) )
        else
          call blanks ( string_length(t%text(lhs)) + 2 - k )
        end if
        prevLHS = lhs
        call output ( ' ->' )
        do j = t%prod_ind(p)+1, t%prod_ind(p+1)-1
          call get_string ( t%text(t%productions(j)), work )
          k = merge(2,1,work(1:1)=='_')
          call output ( ' ' // trim(work(k:)) )
        end do
        call print_action ( p )
        call newLine
      end do
      call output ( 'End of the grammar', advance='yes' )
    end subroutine Print_Grammar

! --------------------------------------------------  Print_State  -----
    subroutine Print_State ( State, Number )
    ! Print the parser configuration in state STATE, i.e., the list of
    ! t%productions, with t%dots in them.
      use String_Table, only: Get_String
      integer, intent(in) :: State
      logical, intent(in), optional :: Number ! Print state number
      integer :: I, J, K, JDOT, P
      logical :: Num
      character(120) :: Work ! To strip underscore off nonterminal names
      num = .false.
      if ( present(number) ) num = number
      do i = t%indbas(state), t%indbas(state+1)-1
        p = t%basProds(i)
        if ( num ) then
          if ( i == t%indbas(state) ) then
            call output ( state, 5 )
          else
            call blanks ( 5 )
          end if
        end if
        call output ( p, 5 )
        j = t%prod_ind(p) ! p = production number
        call get_string ( t%text(t%productions(j)), work )
        k = merge(2,1,work(1:1)=='_')
        call output ( ' ' // trim(work(k:)) )
        call output ( ' ->' )
        jdot = 1
        do j = t%prod_ind(p)+1, t%prod_ind(p+1)-1
          if ( t%dots(i) == jdot ) call output ( ' .' )
          jdot = jdot + 1
          call get_string ( t%text(t%productions(j)), work )
          k = merge(2,1,work(1:1)=='_')
          call output ( ' ' // trim(work(k:)) )
        end do
        if ( t%dots(i) == jdot ) call output ( ' .' )
        call print_action ( p )
        call newLine
      end do
      if ( num ) then ! Print transitions if any
        if ( t%ftrn(state) > t%ftrn(state-1) ) then
          call blanks ( 11 )
          call output ( 'Transitions: ' )
          do i = t%ftrn(state-1)+1, t%ftrn(state)
            call output ( t%tran(i), 5 )
          end do
          call newLine
        end if
      end if
    end subroutine Print_State

! --------------------------------------------  Reduce_Production  -----
    subroutine Reduce_Production ( Prod )
      integer, intent(in) :: Prod ! Index of production to reduce
      integer :: Action       ! action when reducing a production:
                              !      1 = produce a node
                              !      2 = produce a node iff > 1 sons
      integer :: I            ! temporary variable
      integer :: LHS_TREE     ! tree stack index for t%lhs
      integer :: NEW_NODE     ! tree node to make
      integer :: NSONS        ! number of tree parts associated with rhs

      sp = sp - t%lens(prod)    ! Pop RHS-length frames from the stack
      nowsta = stack(sp)%state
      lhs_tree = stack(sp)%tree_index
      nsons = n_tree_stack - lhs_tree
      i = t%act(prod)           ! use i for temp while getting action
      new_node = i / 10
      action = mod(i, 10)

      if ( toggle(par) ) then
        call output ( prod, before='Reduce production ' )
        call output ( t%lens(prod), before=' consuming ' )
        call output ( trim(merge(' symbol ', ' symbols', t%lens(prod)==1)) )
      end if
      if ( i /= 0 ) then
        if ( action==1 .or. action==2 .and. nsons>1 ) then
          call build_tree ( new_node, nsons )
          if ( toggle(par) ) then
            call output ( ', building ' )
            call dump_top_stack ( 0 )
          end if
        else if ( action == 9 ) then
          call do_include
        end if
      else if ( toggle(par) .and. nsons>0 ) then
        call output ( nsons, before=', subsuming ' )
        call output ( trim(merge(' tree node ', ' tree nodes', nsons==1)) )
      end if
      if ( toggle(par) ) then
        call output ( nowsta, before=', returning to state ' )
        if ( show_stack ) then
          call output ( sp, before=', stack(' )
          call output ( stack(sp)%state, before=') = (' )
          call output ( stack(sp)%tree_index, before=',' )
          call output ( ')' )
        end if
        call output ( '.', advance='yes' )
      end if

      newsta = work(nowsta,t%lhs(prod))
      if ( newsta == 0 ) &
        & call catastrophic_error ( 'There is no state transition', &
          & state=nowsta, nonterminal=t%lhs(prod) )
      voc = t%text(t%lhs(prod))
      call transition ( .true. )

    end subroutine Reduce_Production

! ---------------------------------------------------  Transition  -----
    subroutine Transition ( LHS )
      ! Set the current state (Nowsta) to the new state (Newsta)
      ! and push it on the stack with N_Tree_Stack.
      use Allocate_Deallocate, only: Test_Allocate
      use String_Table, only: Display_String
      logical, intent(in) :: LHS ! VOC is index of LHS symbol of production
      integer :: N, S
      integer :: Stat
      type(stack_frame), allocatable :: TempStack(:)

      ! Make sure we have a stack, and it's big enough to push the new state
      n = 0
      if ( allocated(stack) ) n = ubound(stack,1)
      if ( sp >= n ) then
        s = max(100,2*n)
        allocate ( tempStack(s), stat=stat )
        call test_allocate ( stat, moduleName, 'Stack', [1], [s], &
                           & storage_size(stack) / 8 )
        if ( allocated(stack) ) tempStack(1:sp) = stack(1:sp)
        call move_alloc ( tempStack, stack )
      end if

      nowsta = newsta
      sp = sp + 1
      stack(sp) = stack_frame(nowsta,n_tree_stack)

      if ( toggle(par) ) then
        call output ( newsta, before='Enter state ' )
        if ( voc /= 0 ) then
          call output ( ' with ' )
          if ( voc < 0 ) then
            call output ( t%ent(nowsta), before='symbol number ' )
          else if ( lhs ) then
            call display_string ( voc, before='nonterminal symbol ' )
          else
            call dump_1_symbol ( the_token%string_index )
            call print_source ( the_token%where, before=' at ' )
          end if
        end if
        call output ( n_tree_stack, before=', n_tree_stack = ' )
        if ( show_stack ) then
          call output ( sp, before=', stack(' )
          call output ( stack(sp)%state, before=') = (' )
          call output ( stack(sp)%tree_index, before=',' )
          call output ( ')' )
        end if
        call newLine
        if ( show_state ) call print_state ( nowsta )
      end if

    end subroutine Transition

  end subroutine LR_Parser

! ----------------------------------------------  Clean_Up_Parser  -----
  subroutine Clean_Up_Parser
    deallocate ( work, LR_0 )
  end subroutine Clean_Up_Parser

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
! Revision 2.37  2019/07/09 20:20:58  pwagner
! NAG complained if this was used to build lr before; happy now
!
! Revision 2.36  2015/10/27 21:58:53  vsnyder
! Handle a symbol that is not in the terminal alphabet as a syntax error, and
! try to recover, rather than as a crash-burn catastrophic error.
!
! Revision 2.35  2015/10/22 23:50:49  vsnyder
! Announce where bad token has confused the parser
!
! Revision 2.34  2014/09/05 00:28:52  vsnyder
! Correct wrong units for allocation size tracking
!
! Revision 2.33  2014/05/20 23:54:38  vsnyder
! New parser gets its tables from an argument instead of an include
!
! Revision 2.32  2014/01/14 00:59:40  vsnyder
! Remove dependence of Processor_Dependent, which was only for debugging
!
! Revision 2.31  2014/01/14 00:48:59  vsnyder
! Revised LR needs revised parser
!
! Revision 1.1  2014/01/14 00:15:02  vsnyder
! Initial commit of new module for new LR
!
! Revision 2.30  2013/12/12 02:01:54  vsnyder
! Entirely replaced with LR parser
!
! Revision 2.29  2013/10/02 01:35:46  vsnyder
! Add conditional expressions ?...! and variable assignment :=
!
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
