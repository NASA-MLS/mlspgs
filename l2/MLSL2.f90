program MLSL2
  use DECLARATION_TABLE, only: ALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES, LIT_INDICES
  use LEXER_CORE, only: INIT_LEXER
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use OUTPUT_M, only: OUTPUT
! use Open_Init, only: CloseMLSCF, OpenMLSCF
  use PARSER, only: CONFIGURATION
  use STRING_TABLE, only: DO_LISTING
  use TOGGLES, only: CON, GEN, LEVELS, LEX, PAR, SYN, TAB, TOGGLE
  use TREE, only: ALLOCATE_TREE, PRINT_SUBTREE
  use TREE_CHECKER, only: CHECK_TREE
  use TREE_WALKER, only: WALK_TREE_TO_DO_MLS_L2
  use UNITS, only: INIT_UNITS

  implicit NONE

  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from check_tree
  integer :: FIRST_SECTION         ! Index of son of root of first n_cf node
  integer :: I                     ! counter for command line arguments
  integer :: J                     ! index within option
  character(len=80) :: LINE        ! Into which is read the command args
  integer :: ROOT                  ! of the abstract syntax tree

! Initialize the lexer, symbol table, and tree checker's tables:
  call init_lexer ( n_chars=10000, n_symbols=1000, hash_table_size=1003 )
  call allocate_decl ( ndecls=1000 )
  call allocate_tree ( n_tree=10000 )
  call init_tables
  call init_units ( lit_indices )

  i = 1+hp
  do ! Process the command line options to set toggles
    call getarg ( i, line )
    if ( line(1:1) == '-' ) then
      j = 1
      do while ( j < len(line) )
        j = j + 1
        if ( line(j:j) == 'c' ) then
          toggle(con) = .true.
        else if ( line(j:j) == 'g' ) then
          toggle(gen) = .true.
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              j = j + 1
              levels(gen) = ichar(line(j:j)) - ichar('0')
            end if
          end if
        else if ( line(j:j) == 'l' ) then
          toggle(lex) = .true.
        else if ( line(j:j) == 'p' ) then
          toggle(par) = .true.
        else if ( line(j:j) == 'a' ) then
          toggle(syn) = .true.
        else if ( line(j:j) == 'A' ) then
          dump_tree = .true.
        else if ( line(j:j) == 'd' ) then
          do_dump = .true.
        else if ( line(j:j) == 't' ) then
          toggle(tab) = .true.
        else if ( line(j:j) == 'v' ) then
          do_listing = .true.
        end if
      end do
    else    
  exit
    end if
    i = i + 1
  end do

! Parse the L2CF, producing an abstract syntax tree
! call openMLSCF
  call configuration ( root )
! call closeMLSCF
  if ( root <= 0 ) then
    call output ( &
      'A syntax error occurred -- there is no abstract syntax tree', &
      advance='yes' )
  else
    if ( dump_tree ) call print_subtree ( root, 0 )

    ! Check that supra-syntactic conditions are met, e.g. correct
    ! types for fields of commands, correct command order, etc.
    call check_tree ( root, error, first_section )
    if ( do_dump ) call dump_decl
    if ( toggle(syn) ) then
      call output ( 'Begin abstract syntax tree:', advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End abstract syntax tree', advance='yes' )
    end if
    if ( error == 0 .and. first_section /= 0 ) then
      if ( do_dump ) call dump_decl

      ! Now do the L2 processing.
      call walk_tree_to_do_MLS_L2 ( root, error, first_section )
    end if
  end if
end program MLSL2
