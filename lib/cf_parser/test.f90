program TEST
  use DECLARATION_TABLE, only: ALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES
  use LEXER_CORE, only: INIT_LEXER
  use MACHINE ! At least HP, for command lines, and maybe GETARG
  use MLSCF, only: MLSCF_T
  use OUTPUT_M, only: OUTPUT
  use PARSER, only: CONFIGURATION
  use STRING_TABLE, only: DO_LISTING
  use TABLE_DUMPER, only: DUMP_TABLE                 ! For debugging
  use TABLE_GENERATOR, only: GENERATE_TABLE
  use TREE_CHECKER, only: CHECK_TREE
  use TOGGLES, only: CON, GEN, LEX, PAR, SYN, TAB, TOGGLE
  use TREE, only: ALLOCATE_TREE, PRINT_SUBTREE
  use UNITS, only: INIT_UNITS

  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DO_DUMP_EARLY = .false.    ! Dump declaration table before check
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from check_tree
  integer :: HOW_MANY_SECTIONS
  integer :: I      ! counter for command line arguments
  integer :: J      ! index within option
  type(mlscf_t) :: L2CF_DATA
  character(len=80) :: LINE
  integer :: ROOT   ! of the abstract syntax tree

!---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  call init_lexer ( n_chars=10000, n_symbols=1000, hash_table_size=1003 )
  call allocate_decl ( ndecls=1000 )
  call allocate_tree ( n_tree=10000 )
  call init_tables
  call init_units

  i = 1+hp
  do ! Process the command line options to set toggles
    call getarg ( i, line )
    if ( line(1:1) == '-' ) then
      do j = 2, len(line)
        if ( line(j:j) == 'c' ) then
          toggle(con) = .true.
        else if ( line(j:j) == 'g' ) then
          toggle(gen) = .true.
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
        else if ( line(j:j) == 'D' ) then
          do_dump_early = .true.
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

  call configuration ( root )
  if ( root > 0 ) then
    if ( dump_tree ) call print_subtree ( root, 0 )
    if ( do_dump_early ) call dump_decl
    call check_tree ( root, error, how_many_sections=how_many_sections )
    if ( do_dump ) call dump_decl
    if ( toggle(syn) ) then
      call output ( 'Begin abstract syntax tree:', advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End abstract syntax tree', advance='yes' )
    end if
    if ( error == 0 ) then
      if ( do_dump ) call dump_decl
      call generate_table ( root, how_many_sections, l2cf_data )
      call dump_table ( l2cf_data )
    end if
  else
    call output ( &
      'A syntax error occurred -- there is no abstract syntax tree', &
      advance='yes' )
  end if
end program TEST

! $Log$
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
