! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
module GetCF_M
!==============================================================================

  use DECLARATION_TABLE, only: ALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES, LIT_INDICES
  use LEXER_CORE, only: INIT_LEXER
  use MLSCF, only: MLSCF_T
  use OUTPUT_M, only: OUTPUT
  use PARSER, only: CONFIGURATION
  use STRING_TABLE, only: DO_LISTING, IN_UNIT => INUNIT
  use TABLE_DUMPER, only: DUMP_TABLE                 ! For debugging
  use TABLE_GENERATOR, only: GENERATE_TABLE
  use TOGGLES, only: SYN, TOGGLE
  use TREE_CHECKER, only: CHECK_TREE
  use TREE, only: ALLOCATE_TREE, PRINT_SUBTREE
  use UNITS, only: INIT_UNITS

!---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  InitGetCF  -----
  ! Initialize routines needed by GetCF.  This is separate, because it
  ! clears the toggles.  You may want to set them before calling GetCF.

  subroutine InitGetCF ( nChars, nSymbols, nHash, nDecls, nTree )
    integer, optional, intent(in) :: nChars, nSymbols, nHash, nDecls, nTree
      ! nChars => Number of characters in distinct input symbols, default 10000
      ! nSymbols => Number of input symbols, default 1000
      ! nHash => Hash table size, default 1003, should not have small factors
      ! nDecls => Number of declarations of symbols, default 1000
      ! nTree => Number of vertices in abstract syntax tree, default 10000
    integer :: n_Chars, n_Symbols, n_Hash, n_Decls, n_Tree  ! See arguments
    n_chars = 10000
    if ( present(nChars) ) n_chars = nChars
    n_symbols = 1000
    if ( present(nSymbols) ) n_symbols = nSymbols
    n_hash = 1003
    if ( present(nHash) ) n_hash = nHash
    n_decls = 1000
    if ( present(nDecls) ) n_decls = nDecls
    n_tree = 10000
    if ( present(nTree) ) n_tree = nTree

    call init_lexer ( n_chars, n_symbols, hash_table_size=n_hash )
    call allocate_decl ( n_decls )
    call allocate_tree ( n_tree )
    call init_tables
    call init_units ( lit_indices )
  end subroutine InitGetCF

  ! ------------------------------------------------------  GetCF  -----
  ! Parse a Configuration file, producing a data structure described in
  ! the MLSCF module.  InitGetCF must be called first.

  subroutine GetCF ( CF_Data, Status, InUnit, Listing, &
    &                Dump, DumpEarly, DumpTree, DumpTables )
    type(mlscf_t), intent(out) :: CF_DATA    ! The CF data
    integer, optional, intent(out) :: STATUS ! 0 => OK, -1 => Parser error,
      ! >0 => type checking error.
    integer, optional, intent(in) :: InUnit  ! Unit number from which to
      ! read input.  It is not opened or closed here.  If it is negative,
      ! input is read from standard input.  Default -1.
    logical, optional, intent(in) :: Listing ! True => List input, default false
    logical, optional, intent(in) :: Dump, DumpEarly, DumpTree, DumpTables
      ! Dump => Dump declaration table after type checking
      ! DumpEarly => Dump declaration table before type checking
      ! DumpTree => Dump abstract syntax tree after parsing
      ! DumpTables => Dump generated tables
      ! All default .false.  It is possible to cause other output by
      ! setting elements of Toggle in Toggles, q.v.

    logical :: DO_DUMP               ! Dump declaration table
    logical :: DO_DUMP_EARLY         ! Dump declaration table before check
    logical :: DO_DUMP_TREE          ! Dump tree after parsing
    logical :: DO_DUMP_TABLES        ! Dump generated tables
    integer :: ERROR                 ! Error flag from check_tree
    integer :: HOW_MANY_SECTIONS     ! Set by Check_Tree
    integer :: ROOT                  ! of the abstract syntax tree

    in_unit = -1
    if ( present(inUnit) ) in_unit = inunit
    do_dump = .false.
    if ( present(dump) ) do_dump = dump
    do_dump_early = .false.
    if ( present(dumpEarly) ) do_dump_early = dumpEarly
    do_dump_tree = .false.
    if ( present(dumpTree) ) do_dump_tree = dumpTree
    do_dump_tables = .false.
    if ( present(dumpTables) ) do_dump_tables = dumpTables
    do_listing = .false.
    if ( present(listing) ) do_listing = listing

    call configuration ( root )
    if ( root > 0 ) then
      if ( do_dump_tree ) call print_subtree ( root, 0 )
      if ( do_dump_early ) call dump_decl
      call check_tree ( root, error, how_many_sections=how_many_sections )
      if ( do_dump ) call dump_decl
      if ( toggle(syn) ) then
        call output ( 'Begin abstract syntax tree:', advance='yes' )
        call print_subtree ( root, 0 )
        call output ( 'End abstract syntax tree', advance='yes' )
      end if
      if ( error == 0 ) then
        call generate_table ( root, how_many_sections, cf_data )
        if ( do_dump_tables ) call dump_table ( cf_data )
      end if
    else
      error = -1
      call output ( &
        'A syntax error occurred -- there is no abstract syntax tree', &
        advance='yes' )
    end if
    if ( present(status) ) status = error
  end subroutine GetCF
end module GetCF_M

! $Log$
! Revision 2.1  2000/10/11 18:58:37  vsnyder
! Move from lib/cf_parser to lib
!
! Revision 2.2  2000/10/03 01:37:58  vsnyder
! Add the copyright, correct some comments.
!
! Revision 2.1  2000/10/03 00:54:19  vsnyder
! Revised name from getL2CF_m.f90 to getCF_m.f90
!
