! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SYMBOL_TABLE

  use MACHINE, only: IO_ERROR
  use STRING_TABLE, only: ALLOCATE_CHAR_TABLE, ALLOCATE_STRING_TABLE, &
                          DISPLAY_STRING, ADD_CHAR, LOOKUP_AND_INSERT, &
                          STRING_TABLE_SIZE
  use SYMBOL_TYPES, only: CASELESS_LOOK, FIRST_TYPE, INIT_TERMINAL, &
                          INIT_TYPE, LAST_TYPE, T_NULL, T_LAST_TERMINAL, &
                          TERM_TYPES
  use OUTPUT_M, only: OUTPUT
  use TOGGLES, only: TAB, TOGGLE

  implicit NONE
  private

  public :: ADD_TERMINAL, ALLOCATE_SYMBOL_TABLE, DESTROY_SYMBOL_TABLE
  public :: DUMP_1_SYMBOL, DUMP_SYMBOL_CLASS, DUMP_SYMBOL_TYPE
  public :: ENTER_TERMINAL, INIT_SYMBOL_TABLE, SET_SYMBOL, SYMBOL

  integer, save, private :: CLASS_TEXTS(T_NULL: T_LAST_TERMINAL)
  integer, save, private, allocatable :: SYMBOLS(:)
  integer, save, private :: TYPE_TEXTS(FIRST_TYPE: LAST_TYPE) = 0

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
  ! =========================================     ADD_TERMINAL     =====
  integer function ADD_TERMINAL ( TERMINAL )
  ! Add a terminal symbol with token index TERMINAL to the symbol table.
  ! Its text has already been entered into the string table by
  ! STRING_TABLE % ADD_CHAR.
    integer, intent(in) :: TERMINAL

    integer :: WHERE
    logical :: FOUND

    call lookup_and_insert ( where, found, caseless_look(terminal) )
    if ( where > size(symbols) ) then
      call increase_symbols ! String table was expanded
    end if
    if ( toggle(tab) ) then
      call output ( 'Looked up ' ); call output ( terminal )
      call output ( ' ' ); call display_string ( where )
      if ( caseless_look(terminal) ) then; call output ( ' caseless' ); end if
      if ( found ) then
        call output ( ' and found it at ' )
      else
        call output ( ' and inserted it at ' )
      end if
      call output ( where, advance='yes' )
    end if
    if ( .not. found ) then; symbols(where) = terminal; end if
    add_terminal = where
    return
  end function ADD_TERMINAL
  ! ================================     ALLOCATE_SYMBOL_TABLE     =====
  subroutine ALLOCATE_SYMBOL_TABLE ( N_CHARS, N_SYMBOLS, STATUS )
  ! Allocate character, string and symbol tables.
  ! Also does INIT_SYMBOL_TABLE
    integer, intent(in) :: N_CHARS      ! Size of character table
    integer, intent(in) :: N_SYMBOLS    ! Size of string and symbol tables
    integer, intent(out), optional :: STATUS      ! from ALLOCATE

    integer :: STAT
    if ( allocated(symbols) ) then; deallocate ( symbols ); end if
    call allocate_char_table ( n_chars, status )
    if ( present(status) ) then
      if ( status /= 0 ) then; return; end if
    end if
    call allocate_string_table ( n_symbols, status )
    if ( present(status) ) then
      if ( status /= 0 ) then; return; end if
    end if
    allocate ( symbols(n_symbols), stat=stat )
    if ( present(status) ) then
      status = stat
      if ( status /= 0 ) then; return; end if
    end if
    if ( stat /= 0 ) then
      call io_error &
      ( 'SYMBOL_TABLE%ALLOCATE_SYMBOL_TABLE-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    call init_symbol_table
    return
  end subroutine ALLOCATE_SYMBOL_TABLE
  ! =================================     DESTROY_SYMBOL_TABLE     =====
  subroutine DESTROY_SYMBOL_TABLE ( Status )
    integer, intent(out), optional :: Status ! From deallocate
    if ( present(status) ) then
      deallocate ( symbols, stat=status )
    else
      deallocate ( symbols )
    end if
  end subroutine DESTROY_SYMBOL_TABLE
  ! ========================================     DUMP_1_SYMBOL     =====
  subroutine DUMP_1_SYMBOL ( SYMBOL, ADVANCE )
  ! Print the symbol index, type and text.
    integer, intent(in) :: SYMBOL
    character(len=*), intent(in), optional :: ADVANCE
    call output ( symbol ); call output ( ': ' )
    call display_string ( type_texts(term_types(symbols(symbol))) )
    call output ( ' ' )
    call display_string ( symbol, advance=advance )
    return
  end subroutine DUMP_1_SYMBOL
  ! =====================================     DUMP_SYMBOL_CLASS     =====
  subroutine DUMP_SYMBOL_CLASS ( SYM_CLASS, ADVANCE )
  ! Print the symbol text, given its class, e.g. T_IDENTIFIER
    integer, intent(in) :: SYM_CLASS
    character(len=*), intent(in), optional :: ADVANCE
    call display_string ( class_texts(sym_class), advance )
    return
  end subroutine DUMP_SYMBOL_CLASS
  ! =====================================     DUMP_SYMBOL_TYPE     =====
  subroutine DUMP_SYMBOL_TYPE ( SYMTYP, ADVANCE )
  ! Print the symbol type, e.g. <reserved_word>.
    integer, intent(in) :: SYMTYP
    character(len=*), intent(in), optional :: ADVANCE
    call display_string ( type_texts(term_types(symbols(symtyp))), &
                          advance )
    return
  end subroutine DUMP_SYMBOL_TYPE
  ! =======================================     ENTER_TERMINAL     =====
  integer function ENTER_TERMINAL ( TEXT, TERMINAL )
  ! Put the text of a terminal symbol at the end of the character table,
  ! then enter it and its terminal index into the symbol table.
    character(len=*), intent(in) :: TEXT
    integer, intent(in) :: TERMINAL
    call add_char ( text )
    enter_terminal = add_terminal ( terminal )
    return
  end function ENTER_TERMINAL
  ! ====================================     INIT_SYMBOL_TABLE     =====
  subroutine INIT_SYMBOL_TABLE
  ! Initialize the symbol table, the table of string indices of terminal
  ! names, and the table of terminal type names.
    logical :: FOUND          ! Argument for LOOKUP_AND_INSERT
    integer :: I              ! Loop inductor
    integer :: WHERE          ! Where was the symbol inserted?
    do i = t_null, t_last_terminal
      call init_terminal ( i )
      where = add_terminal ( i )
      class_texts(i) = where
    end do
    ! This loop puts token types into the string table.
    ! They're only used for debugging output.
    do i = first_type, last_type
      call init_type ( i )
      call lookup_and_insert ( where, found, .false. )
      ! It's OK if it found one -- maybe it's a terminal text, too.
      if ( .not. found ) then; symbols(where) = t_null; end if
      type_texts(i) = where
    end do
  end subroutine INIT_SYMBOL_TABLE
  ! ===========================================     SET_SYMBOL     =====
  subroutine SET_SYMBOL ( STRING, CLASS )
  ! Set the class of the symbol at STRING to CLASS
    integer, intent(in) :: STRING
    integer, intent(in) :: CLASS
    if ( string > size(symbols) ) then
      call increase_symbols ! String table was expanded
    end if
    symbols(string) = class
  end subroutine SET_SYMBOL
  ! ===============================================     SYMBOL     =====
  integer function SYMBOL ( STRING )
  ! Return the symbol class of the STRING
    integer, intent(in) :: STRING
    symbol = symbols(string)
    return
  end function SYMBOL

! =====     Private procedures    ======================================
  ! -------------------------------------     INCREASE_SYMBOLS     -----
  subroutine INCREASE_SYMBOLS ( STAT )
  ! Increase the upper bound for the symbol table to the upper bound of
  ! the string table.
    integer, optional, intent(out) :: STAT
    integer :: MY_STAT
    integer, allocatable :: OLD_SYMBOL(:)
    allocate ( old_symbol(size(symbols)), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'SYMBOL_TABLE%INCREASE_SYMBOLS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    old_symbol = symbols
    deallocate ( symbols )
    allocate( symbols(string_table_size()), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'SYMBOL_TABLE%INCREASE_SYMBOLS-E- Unable to allocate storage', &
      stat, '' )
      stop
    end if
    symbols(:size(old_symbol)) = old_symbol
    deallocate ( old_symbol )
    return
  end subroutine INCREASE_SYMBOLS

end module SYMBOL_TABLE

! $Log$
! Revision 2.2  2001/04/21 00:37:30  vsnyder
! Add 'Destroy_Symbol_Table' routine
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
