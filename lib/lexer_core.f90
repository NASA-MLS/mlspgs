! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module LEXER_CORE
! Provides the token type and the initialization routine used by all
! lexers, and a routine to print the source text line and column number.

  implicit NONE
  private
  public :: NEED, TOKEN, INIT_LEXER, PRINT_SOURCE

  type :: TOKEN
    integer :: CLASS               ! Token class = token index
    integer :: STRING_INDEX        ! Position in string table
    logical :: PSEUDO              ! "token is pseudo-terminal"
    integer :: SOURCE              ! 256*srcLine + srcCol
  end type TOKEN

  logical, save :: NEED = .true.   ! "Lexer needs to read a character"

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine INIT_LEXER ( N_CHARS, N_SYMBOLS, HASH_TABLE_SIZE, STAT )
  ! Allocate the character, symbol and hash tables, which automatically
  ! initializes them.
    use STRING_TABLE, only: ALLOCATE_HASH_TABLE
    use SYMBOL_TABLE, only: ALLOCATE_SYMBOL_TABLE
    use TOGGLES, only: INIT_TOGGLE

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

  subroutine PRINT_SOURCE ( SOURCE, ADVANCE, BEFORE )
  ! Output "line LLL, column CCC"
    use OUTPUT_M, only: OUTPUT

    integer, intent(in) :: SOURCE  ! 256*srcLine + srcCol
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    if ( present(before) ) call output ( before )
    call output ( 'line ' )
    call output ( source/256 )
    call output ( ', column ' )
    call output ( mod(source,256), advance=advance )
  end subroutine PRINT_SOURCE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LEXER_CORE

! $Log$
! Revision 2.6  2010/03/18 02:36:39  vsnyder
! Move USE from module scope to procedure scope
!
! Revision 2.5  2010/03/18 02:31:18  vsnyder
! Add BEFORE to PRINT_SOURCE
!
! Revision 2.4  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
