! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LEXER_CORE
! Provides the token type and the initialization routine used by all
! lexers, and a routine to print the source text line and column number.

  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: ALLOCATE_HASH_TABLE
  use SYMBOL_TABLE, only: ALLOCATE_SYMBOL_TABLE
  use TOGGLES, only: INIT_TOGGLE

  implicit NONE
  public

  type :: TOKEN
    integer :: CLASS               ! Token class = token index
    integer :: STRING_INDEX        ! Position in string table
    logical :: PSEUDO              ! "token is pseudo-terminal"
    integer :: SOURCE              ! 256*srcLine + srcCol
  end type TOKEN

  logical, save :: NEED = .true.   ! "Lexer needs to read a character"

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine INIT_LEXER ( N_CHARS, N_SYMBOLS, HASH_TABLE_SIZE, STAT )
  ! Allocate the character, symbol and hash tables, which automatically
  ! initializes them.
    integer, intent(in) :: N_CHARS      ! Initial size of character table
    integer, intent(in) :: N_SYMBOLS    ! Initial sizes of string and
                                        ! symbol tables
    integer, intent(in) :: HASH_TABLE_SIZE   ! Duh
    integer, intent(out), optional :: STAT   ! Status from a called routine
    integer :: MY_STAT
    call allocate_hash_table ( hash_table_size, my_stat )
    if (my_stat == 0 ) &
      call allocate_symbol_table ( n_chars, n_symbols, my_stat )
    if ( my_stat == 0 ) call init_toggle
    need = .true.
    if ( present(stat) ) stat = my_stat
    return
  end subroutine INIT_LEXER

  subroutine PRINT_SOURCE ( SOURCE, ADVANCE )
  ! Output "line LLL, column CCC"
    integer, intent(in) :: SOURCE  ! 256*srcLine + srcCol
    character(len=*), intent(in), optional :: ADVANCE
    call output ( 'line ' )
    call output ( source/256 )
    call output ( ', column ' )
    call output ( mod(source,256), advance=advance )
  end subroutine PRINT_SOURCE

end module LEXER_CORE

! $Log$
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
