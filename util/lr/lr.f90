! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program LR

  use, intrinsic :: ISO_Fortran_Env, only: Error_Unit

  use Analysis, only: Analyz
  use Chain_Context_Lists, only: CHNCSL
  use Declaration_Table, only: Dump_Decl
  use Declare_Vocabulary_m, only: Declare_Vocabulary
  use Flatten_m, only: Flatten
  use Generate_Table, only: GENTAB
  use IO_Stuff, only: Get_Lun, List_Unit, Table_Unit
  use Lexer_Core, only: Init_Lexer
  use Lists, only: Lists_Init
  use Output_m, only: Output, OutputOptions
  use Parser, only: Clean_Up_Parser, LR_Parser
  use Parser_Table_m, only:  Destroy_Parser_Table, Parser_Table_t
  use Parser_Tables_LR, only: Init_Parser_Table
  use Print_Set, only: PNTSET
  use Print_The_Grammar_m, only: Print_The_Grammar
  use Print_The_Vocabulary_m, only: Print_The_Vocabulary ! Also sorts it
  use String_Table, only: Open_Input
  use Symbol_Table, only: Dump_Symbol_Table
  use Tables, only: Actions, Productions, Prod_Ind, Vocab
  use Toggles, only: GEN, Levels, LEX, PAR, Switches, TAB, Toggle
  use Toggles_LR, only: Toggle_LR => Toggle
  use Tree, only: Allocate_tree, Print_subtree
  use Xref, only: Cross_Reference

  logical :: Dump_Symbols = .false.
  logical :: Dump_Tree = .false.
  integer :: Depth = 0     ! For dumping the abstract syntax tree
  logical :: Error = .false. ! Somebody detected an error
  integer :: I, J          ! Loop indices
  character(1023) :: IOMSG ! From Open
  integer :: IOSTAT        ! From Open
  character(1023) :: Line  ! From command line
  character(1023) :: Listing = '' ! list file name
  integer :: LNADQT = 1    ! Last inadequate state number, zero if grammar is LR.
                           ! Initially nonzero in case the grammar isn't analyzed.
  type(Parser_Table_t) :: Parser_Table
  integer :: Root          ! Of the abstract syntax tree
  character(1023) :: Table = '' ! table file name

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957, &
    & DEBUG=0 )
  call allocate_tree ( n_tree=2000000 )

  ! Process command line
  i = 0
  do
    i = i + 1
    call get_command_argument ( i, line )
    if ( line(1:1) /= '-' ) exit
    j = 1
    do
      j = j + 1
      if ( j > len(line) ) exit
      select case ( line(j:j) )
      case ( ' ' )
        exit
      case ( 'A' ) ! Dump abstract syntax tree
        dump_tree = .true.
      case ( 'd' ) ! Trace declaration table actions
        toggle(tab) = .true.
        call digit_after_option ( j, tab )
      case ( 'g' )
        toggle(gen) = .true.
        call digit_after_option ( j, gen )
      case ( 'l' ) ! Specify listing file
        if ( line(j+1:) == ' ' ) then
          i = i + 1
          call get_command_argument ( i, line(j+1:) )
        end if
        listing = line(j+1:)
        j = len(line)+1
      case ( 'L' ) ! Trace lexer
        toggle(lex) = .true.
        call digit_after_option ( j, lex )
      case ( 'o' ) ! Output tables
        if ( line(j+1:) == ' ' ) then
          i = i + 1
          call get_command_argument ( i, line(j+1:) )
        end if
        table = line(j+1:)
        j = len(line)+1
      case ( 'P' ) ! Trace parser
        toggle(par) = .true.
        call digit_after_option ( j, par )
      case ( 's' ) ! Dump symbol table after the parser finishes
        dump_symbols = .true.
      case ( 'S' ) ! Set switches for debugging
        if ( line(j+1:) == ' ' ) then
          i = i + 1
          call get_command_argument ( i, line(j+1:) )
        end if
        if ( line(j+1:j+1) == '?' ) call Switch_Usage
        switches(len_trim(switches)+1:) = ',' // trim(line(j+1:))
        j = len(line) + 1
      case ( '2' : '4', 'M', 'X' )
        toggle_lr(iachar(line(j:j))) = 1
      case default
        call usage
      end select
    end do
  end do

  if ( line(1:1) == ' ' ) call usage

  ! Open input file
  call open_input ( line )

  ! Open output table file
  call get_command_argument ( i+1, line )
  if ( line /= '' ) table = line
  if ( table /= '' ) then
    call get_lun ( table_unit )
    open ( table_unit, file=table, form='formatted', iostat=iostat, &
      & iomsg=iomsg )
    if ( iostat /= 0 ) then
      write ( error_unit, '(3a,i0)' ) 'Error: Unable to open table output file "', &
        & trim(table), '", IOSTAT = ', iostat
      write ( error_unit,  '(a)') trim(iomsg)
      stop 1
    end if
  end if

  ! Open output listing file
  call get_command_argument ( i+2, line )
  if ( line /= '' ) listing = line
  if ( listing /= '' ) then
    call get_lun ( list_unit )
    open ( list_unit, file=listing, form='formatted', iostat=iostat, &
      & iomsg=iomsg )
    if ( iostat /= 0 ) then
      write ( error_unit, '(3a,i0)' ) 'Error: Unable to open list file "', &
        & trim(listing), '", IOSTAT = ', iostat
      write ( error_unit,  '(a)') trim(iomsg)
      stop 1
    end if
  end if

  outputOptions%prUnit = list_unit

  ! Parse the grammar, producing an abstract syntax tree
  call init_parser_table ( parser_table )
  call lr_parser ( root, parser_table )
  call destroy_parser_table ( parser_table )
  call clean_up_parser

  if ( dump_symbols ) call dump_symbol_table

  if ( root < 0 ) then
    call output ( 'A syntax error occurred; there is no abstract syntax tree.', &
      & advance='yes' )
    stop
  end if

  if ( dump_tree ) then
    call output ( 'Abstract syntax tree:', advance='yes' )
    call print_subtree ( root, depth )
    call output ( 'End of abstract syntax tree', advance='yes' )
  end if

  ! Declare symbols in the vocabulary according to whether they are
  ! terminals, nonterminals, vocabulary names, or actions
  call declare_vocabulary ( root, error )

  ! Sort the vocabulary symbols first according to whether they are
  ! terminals, nonterminals, vocabulary names, or actions, then according
  ! to their text.  Then print the vocabulary.  The variable VOCAB
  ! is the permutation vector for the sort; its values are string indices.
  call print_the_vocabulary ( vocab )

  if ( .not. error ) then
    ! Flatten the abstract syntax tree into the data structures used to
    ! analyze the grammar.  Also determine the goal symbol, make sure
    ! every symbol is connected to the goal, and make sure every nonterminal
    ! derives at least one terminal symbol.
    call flatten ( root, vocab, error )
  end if

  if ( toggle(tab) ) then
    call dump_symbol_table
    call dump_decl
  end if

  if ( .not. error ) then

    call print_the_grammar ( prod_ind, productions, actions, vocab )

    ! Print the vocabulary cross reference
    call cross_reference ( productions, prod_ind, vocab )

    ! Set up to use the guts of the original LR from Shannon and Wetherell
    call lists_init
    call analyz
    call pntset ( lnadqt )

    if ( lnadqt /= 0 ) then
      call output ( lnadqt, &
        & before='Error: Grammar is not LR(1).  Last inadequate state is ', &
        & advance='yes' )
      stop 1
    end if

    if ( table_unit > 0 ) then ! -o option or table_file field was specified
      call chncsl
      call gentab ( table_unit, vocab )
    end if

  end if

contains

  subroutine Digit_After_Option ( J, Tog )
    character(len=*), parameter :: Digits = '0123456789'
    integer, intent(inout) :: J  ! Position in LINE
    integer, intent(in) :: Tog   ! Subscript of Toggles and Levels array
    integer :: I, K              ! Loop inductor, subscript
    k = j + 1
    do i = k, len(line)
      if ( verify(line(i:i),digits) /= 0 ) exit
    end do
    if ( i > k ) then
      j = i
      read ( line(k:i-1), * ) levels(tog)
    end if
  end subroutine Digit_After_Option

  subroutine Switch_Usage
    print '(a)', 'Switch_Usage:'
    print '(a)', ' dprod => Dump "productions" and related arrays in Flatten'
    print '(a)', ' dvoc => Trace declaration of vocabulary and dump a summary'
    stop
  end subroutine Switch_Usage

  subroutine Usage
    call get_command_argument ( 0, line )
    print '(3a)', 'Usage: ', trim(line), ' [options] input_file [ table_file [ listing_file ] ]'
    print '(a)', ' Options:'
    print '(a)', '  -A => Print syntax tree'
    print '(a)', '  -d[#] => Trace declaraction table actions at level #'
    print '(a)', '  -g[#] => Trace generation at level #'
    print '(a)', '  -l[ ]file => Specify listing file, default standard output,'
    print '(a)', '               overridden by listing_file field if both specified'
    print '(a)', '  -L[#] => Trace lexer at level #'
    print '(a)', '  -o[ ]file => Specify table file, default no table output,'
    print '(a)', '               overridden by table_file field if both specified'
    print '(a)', '  -M => Turn off printing the parsing automaton ("machine")'
    print '(a)', '  -P[#] => Trace parser at level #; bits of # mean'
    print '(a)', '           1 -> Show also top of stack'
    print '(a)', '           2 -> Show state at every transition'
    print '(a)', '           4 -> Show grammar and automaton when starting'
    print '(a)', '           8 -> Show work table when starting'
    print '(a)', '  -s => Dump symbol table'
    print '(a)', '  -S[ ]string => Add string to debugging switches, see -S? for more'
    print '(a)', '  -X => Turn off printing the cross reference of the parsing automaton'
    print '(a)', '  -2 => Print nullable and first sets'
    print '(a)', '  -3 => Trace context set creation and destruction'
    print '(a)', '  -4 => Trace list element creation and destruction'
    print '(a)', '  -anything else => this output'
    stop
  end subroutine Usage

end program LR

! $Log$
! Revision 1.5  2014/04/09 23:54:01  vsnyder
! Don't try to make the parser if there's a syntax error
!
! Revision 1.4  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
