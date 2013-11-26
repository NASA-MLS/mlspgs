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
  public :: Get_Where, NEED, TOKEN, INIT_LEXER, PRINT_SOURCE, Where_T

  type :: Where_T
    integer :: FILE = 0            ! String index
    integer :: SOURCE = 0          ! 256*srcLine + srcCol
  end type Where_T

  type :: TOKEN
    integer :: CLASS               ! Token class = token index
    integer :: STRING_INDEX        ! Position in string table
    logical :: PSEUDO              ! "token is pseudo-terminal"
    type(where_t) :: Where
  end type TOKEN

  logical, save :: NEED = .true.   ! "Lexer needs to read a character"

  interface Get_Where
    module procedure Get_Where_Integer, Get_Where_Where
  end interface

  interface Print_Source
    module procedure Print_Source_Integer, Print_Source_Where
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine INIT_LEXER ( N_CHARS, N_SYMBOLS, HASH_TABLE_SIZE, STAT )
  ! Allocate the character, symbol and hash tables, which automatically
  ! initializes them.
    use STRING_TABLE, only: ALLOCATE_HASH_TABLE, INIT_STRING_TABLE
    use SYMBOL_TABLE, only: ALLOCATE_SYMBOL_TABLE
    use TOGGLES, only: INIT_TOGGLE

    integer, intent(in) :: N_CHARS      ! Initial size of character table
    integer, intent(in) :: N_SYMBOLS    ! Initial sizes of string and
                                        ! symbol tables
    integer, intent(in) :: HASH_TABLE_SIZE   ! Duh
    integer, intent(out), optional :: STAT   ! Status from a called routine
    integer :: MY_STAT

    ! Since Lexer is the one using get_char subroutine,
    ! which uses the variables being initialize in
    ! init_string_table, init_string_table may as well
    ! be called here
    call init_string_table

    call allocate_hash_table ( hash_table_size, my_stat )
    if (my_stat == 0 ) &
      call allocate_symbol_table ( n_chars, n_symbols, my_stat )
    if ( my_stat == 0 ) call init_toggle
    need = .true.
    if ( present(stat) ) stat = my_stat
  end subroutine INIT_LEXER

  subroutine Get_Where_Integer ( Where, Text, File, Before, After )
    use String_Table, only: Get_String
    integer, intent(in) :: Where
    character(len=*), intent(out) :: Text
    integer, intent(in), optional :: File
    character(len=*), intent(in), optional :: Before, After
    integer :: L
    if ( present(before) ) then
      text = before
      l = len(before) + 1
    end if
    write ( text(l:), 1 ) where/256, mod(where,256)
  1 format ( 'line ', i0, ', column ', i0 )
    if ( present(file) ) then
      if ( file > 0 ) then
        l = len_trim(text) + 1
        text(l:l+3) = ' in '
        call get_string ( file, text(l+4:) )
      end if
    end if
    if ( present(after) ) then
      l = len_trim(text) + 1
      text(l:) = after
    end if
  end subroutine Get_Where_Integer

  subroutine Get_Where_Where ( Where, Text, Before, After  )
    type(where_t), intent(in) :: Where
    character(len=*), intent(out) :: Text
    character(len=*), intent(in), optional :: Before, After
    call get_where ( where%source, text, where%file, before, after )
  end subroutine Get_Where_Where

  subroutine PRINT_SOURCE_INTEGER ( SOURCE, ADVANCE, BEFORE, After, File )
  ! Output "line LLL, column CCC"
    use String_Table, only: Display_String
    use OUTPUT_M, only: NewLine, OUTPUT

    integer, intent(in) :: SOURCE  ! 256*srcLine + srcCol
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE, After
    integer, intent(in), optional :: File ! From include stack
    logical :: MyAdv
    myAdv = .false.
    if ( present(advance) ) then
      ! This kludge is necessary because "call output ( '', advance=advance )
      ! outputs one blank character, instead of zero.
      myAdv = advance(1:1) == 'y' .or. advance(1:1) == 'Y'
    end if
    if ( present(before) ) call output ( before )
    call output ( source/256 , before='line ' )
    call output ( mod(source,256), before=', column ' )
    if ( present(file) ) then
      if ( file /= 0 ) call display_string ( file, before=' in ', strip=.true. )
    end if
    if ( present(after) ) call output ( after )
    if ( myAdv ) call newLine
  end subroutine PRINT_SOURCE_INTEGER

  subroutine Print_Source_Where ( Where, Advance, Before, After )
    type(where_t), intent(in) :: Where
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE, After
    call print_source ( where%source, advance, before, after, where%file )
  end subroutine Print_Source_Where
    
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
! Revision 2.9  2013/11/26 22:44:48  vsnyder
! Spiff up dump
!
! Revision 2.8  2013/09/24 23:08:14  vsnyder
! Add Get_Where, version of Print_Source that uses Where_t
!
! Revision 2.7  2010/08/05 17:45:34  honghanh
! Adding subroutine init_string_table in string_table module
!
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
