! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LEXER_M
! Lexer for SMALL

  use ERROR_HANDLER, only: ERROR_WALKBACK
  use LEXER_CORE, only: NEED, TOKEN
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: ADD_CHAR, CLEAR_STRING, ENTER_STRING, EOF, EOL, &
                          GET_CHAR, HOW_MANY_STRINGS, LOOKUP_AND_INSERT, &
                          NEW_LINE, SOURCE_COLUMN, SOURCE_LINE
  use SYMBOL_TABLE, only: ADD_TERMINAL, DUMP_1_SYMBOL, ENTER_TERMINAL, &
                          SET_SYMBOL, SYMBOL
  use SYMBOL_TYPES ! everything
  use TOGGLES, only: CON, EMIT, GEN, LEX, PAR, SYN, TAB, TOGGLE ! (everything)
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine LEXER ( THE_TOKEN )
    type(token), intent(out) :: THE_TOKEN

    ! Parameters for character classes
    integer, parameter :: LETTER = 1
    integer, parameter :: DIGIT = 2
    integer, parameter :: PUNCT = 3
    integer, parameter :: OP_CHAR = 4
    integer, parameter :: DOT = 5
    integer, parameter :: SPACES = 6
    integer, parameter :: EOL_IN = 7
    integer, parameter :: EOF_IN = 8
    integer, parameter :: CMT = 9
    integer, parameter :: MORE = 10
    integer, parameter :: UNDER = 11
    integer, parameter :: CONT = 12
    integer, parameter :: QUOTE = 13
    integer, parameter :: APOST = 14
    integer, parameter :: TOG_CH = 15

    ! Parameters for ERROR routine's CODE argument
    integer, parameter :: DOUBLE_UNDER = 1
    integer, parameter :: INCOMPLETE = 2
    integer, parameter :: UNREC_CHAR = 3
    integer, parameter :: UNREC_TOKEN = 4

    ! Parameters for states
    integer, parameter :: START = 1
    integer, parameter :: ID = 2
    integer, parameter :: NUM = 3    ! digit(digit|_)*
    integer, parameter :: NUM2 = 4   ! digit(digit|_)* . digit(digit|_)*
    integer, parameter :: NUM3 = 5   ! digit(digit|_)* . digit(digit|_)* E
    integer, parameter :: NUM4 = 6   ! digit(digit|_)* . digit(digit|_)* E +|-
    integer, parameter :: NUM5 = 7   ! digit(digit|_)* . digit(digit|_)* E (+|-)?
                                     ! digit(digit|_)*
    integer, parameter :: OP = 8
    integer, parameter :: PUN = 9 
    integer, parameter :: CONTIN = 10
    integer, parameter :: SQ1 = 11   ! First string-started-by-quote state
    integer, parameter :: SQ2 = 12   ! Second string-started-by-quote state
    integer, parameter :: SA1 = 13   ! First string-started-by-apost state
    integer, parameter :: SA2 = 14   ! Second string-started-by-apost state
    integer, parameter :: TOG = 15

    character, save :: CH          ! Current character
    integer :: CLASS               ! Classification of current character
    integer, parameter :: CLASSES(0:127) = (/ &
!   NUL     SOH     STX     ETX     EOT=^D  ENQ     ACK     BEL
    more,   more,   more,   more,   eof_in, more,   more,   more,   & ! 000
!   BS      HT      LF      VT      FF      CR      SO      SI
    more,   spaces, eol_in, more,   more,   spaces, more,   more,   & ! 010
!   DLE     DC1     DC2     DC3     DC4     NAK     SYN     ETB
    more,   more,   more,   more,   more,   more,   more,   more,   & ! 020
!   CAN     EM      SUB     ESC     FS      GS      RS      US
    more,   more,   more,   more,   more,   more,   more,   more,   & ! 030
!   SPACE   !       "       #       $       %       &       '
    spaces, more,   quote,  more,   cont,   more,   more,   apost,  & ! 040
!   (       )       *       +       ,       -       .       /
    punct,  punct,  op_char,op_char,punct,  op_char,dot,    op_char,& ! 050
!   0       1       2       3       4       5       6       7
    digit,  digit,  digit,  digit,  digit,  digit,  digit,  digit,  & ! 060
!   8       9       :       ;       <       =       >       ?
    digit,  digit,  op_char,cmt,    more,   op_char,more,   more,   & ! 070
!   @       A       B       C       D       E       F       G
    tog_ch, letter, letter, letter, letter, letter, letter, letter, & ! 100
!   H       I       J       K       L       M       N       O
    letter, letter, letter, letter, letter, letter, letter, letter, & ! 110
!   P       Q       R       S       T       U       V       W
    letter, letter, letter, letter, letter, letter, letter, letter, & ! 120
!   X       Y       Z       [       \       ]       ^       _
    letter, letter, letter, punct,  more,   punct,  more,   under,  & ! 130
!   `       a       b       c       d       e       f       g
    more,   letter, letter, letter, letter, letter, letter, letter, & ! 140
!   h       i       j       k       l       m       n       o
    letter, letter, letter, letter, letter, letter, letter, letter, & ! 150
!   p       q       r       s       t       u       v       w
    letter, letter, letter, letter, letter, letter, letter, letter, & ! 160
!   x       y       z       {       |       }       ~       DEL
    letter, letter, letter, more,   cont,   more,   more,   more   /) ! 170
    
    integer :: COLUMN              ! Character index of start of token in
                                   ! CHAR_TABLE
    logical :: FOUND               ! Argument for LOOKUP_AND_INSERT
    integer :: I                   ! Loop inductor
    integer :: SOURCE_START        ! Source reference for start of token
                                   ! ( 256*line + column )
    integer :: STATE               ! Current state of lexer DFA
    integer :: STRING_INDEX        ! Index of token in string table

    state = start
    do
      if ( need ) then; call get_char ( ch ); need = .false.; end if
      ! Classify next character
      class = classes(iachar(ch))

!     Recognize a token according to the following DFA.  Column labels are
!     Classes, row labels are DFA states.  Table elements are DFA states or
!     actions to be taken when a token is recognized.
!
!         letter digit opchar punct under  cmt eolin eofin space  cont  more
!  start    id    num    op    pun   id    cmt    6    1   start  cont start
!  id       id    id      2      2   id     2     2    2     2     2     2
!  num       3    num     3      3   num    3     3    3     3     3     3
!  num2   e=num3 num2     3      3  num2    3     3    3     3     3     3
!         else 3
!  num3     10   num5  +-=num4  10   10    10    10   10    10    10    10
!                      else 10
!  num4     10   num5    10     10   10    10    10   10    10    10    10
!  num5      3   num5     3      3  num5    3     3    3     3     3     3
!  op        4     4      4      4    4     4     4    4     4     4     4
!  pun       5     5      5      5    5     5     5    5     5     5     5
!  cmt      cmt   cmt    cmt    cmt  cmt   cmt  start  1    cmt   cmt   cmt
!  cont      7     7      7      7    7    cmt  start  1   cont    7     7
!  sq1      sq1   sq1    sq1    sq1  sq1   sq1    9    9    sq1   sq1   sq1
!  sq2       8     8      8      8    8     8     8    8     8     8     8
!  sa1      sa1   sa1    sa1    sa1  sa1   sa1    9    9    sa1   sa1   sa1
!  sa2       8     8      8      8    8     8     8    8     8     8     8
!
!          quote apost  dot
!  start    sq1   sa1    4
!  id        2     2     2
!  num       3     3   num2
!  num2      3     3     3
!  num3     10    10    10
!  num4     10    10    10
!  num5      3     3     3
!  op        4     4     4
!  pun       5     5     5
!  cmt      cmt   cmt   cmt
!  cont      7     7     7
!  sq1      sq2   sq1   sq1
!  sq2      sq1    8     8
!  sa1       8    sa2    8
!  sa2       8    sa1    8

!     The actions are:
!     1.   define EOF token.
!     2.   lookup in symbol table.  If found use what is found, else define
!          identifier.
!     3.   define integer token.
!     4.   lookup in symbol table.  If found use what is found, else define
!          unknown operator token.
!     5.   lookup in symbol table.  If found use what is found, else define
!          unknown punctuator token.
!     6.   define end-of-statement token.
!     7.   define "error -- not a comment or end-of-* after continuation"
!     8.   define string token.
!     9.   define "error -- incomplete string"
!    10.   define "error -- incomplete floating-point number"

      select case ( state )
      case ( start )
        column = source_column
        source_start = 256 * source_line + source_column
        select case ( class )
        case ( letter, under );    state = id
        case ( digit );            state = num
        case ( punct );            state = pun
        case ( op_char );          state = op
        case ( dot )
          need = .true.
          call add_char ( '.' )
          call lookup_and_insert ( string_index, found ) ! will be found
                                                         ! if tables OK
          the_token = token( symbol(string_index), string_index, .false., &
                               source_start )
    exit ! main lexer loop
        case ( cmt )
          call new_line
          state = start
          need = .true.
        case ( spaces );           need = .true.
        case ( cont );             state = contin; need = .true.
        case ( eof_in )
          string_index = enter_terminal ( '<eof>', t_end_of_input )
          the_token = token( t_end_of_input, string_index, .false., &
                             source_start )
          call set_symbol ( string_index, t_end_of_input )
    exit ! main lexer loop
        case ( eol_in )
          string_index = enter_terminal ( '<eos>', t_end_of_stmt )
          the_token = token( t_end_of_stmt, string_index, .false., &
                             source_start )
          call set_symbol ( string_index, t_end_of_stmt )
          need = .true.
    exit ! main lexer loop
        case ( quote );
          call add_char ( ch )
          need = .true.
          state = sq1
        case ( apost );
          call add_char ( ch )
          need = .true.
          state = sa1
        case ( tog_ch );           state = tog; need = .true.
        case default
          call error ( unrec_char )
          do ! swallow unknown characters, spaces and ends-of-lines
            call get_char ( ch )
            select case ( classes(iachar(ch)) )
            case ( more, spaces, eol_in )
          exit
            case default
            end select
          end do
          string_index = enter_terminal ( '<aft>', t_aft_cont )
          the_token = token( t_aft_cont, string_index, .false., source_start )
          need = .false.
        end select ! class
      case ( id )
        select case ( class )
        case ( letter, digit, under )
          call add_char ( cap(ch) )
          need = .true.
        case default
          string_index = add_terminal ( t_identifier )
          the_token = token( symbol(string_index), string_index, &
                             term_types(symbol(string_index)) /= res_word, &
                             source_start )
    exit ! main lexer loop
        end select ! class
      case ( num ) ! digit(digit|_)* 
        select case ( class )
        case ( digit ); call add_char ( ch ); need = .true.
        case ( under ); need = .true.
        case ( dot )
          call add_char ( ch )
          need = .true.
          state = num2
        case default
          string_index = add_terminal ( t_number )
          the_token = token( t_number, string_index, .true., &
                             source_start )
    exit ! main lexer loop
        end select ! class
      case ( num2 ) ! digit(digit|_)* . digit(digit|_)*
        select case ( class )
        case ( digit ); call add_char ( ch ); need = .true.
        case ( under ); need = .true.
        case default
          if ( ch /= 'e' .and. ch /= 'E' ) then
            string_index = add_terminal ( t_number )
            the_token = token( t_number, string_index, .true., &
                               source_start )
    exit ! main lexer loop
          end if
          call add_char ( ch )
          need = .true.
          state = num3
        end select ! class
      case ( num3 ) ! digit(digit|_)* . digit(digit|_)* E
        if ( class == digit ) then
          call add_char ( ch )
          need = .true.
          state = num5
        else if ( ch == '-' .or. ch == '+' ) then
          call add_char ( ch )
          need = .true.
          state = num4
        else
          string_index = add_terminal ( t_inc_num )
          the_token = token( t_inc_num, string_index, .true., &
                             source_start )
    exit ! main lexer loop
        end if
      case ( num4 ) ! digit(digit|_)* . digit(digit|_)* E +|-
        if ( class /= digit ) then
          string_index = add_terminal ( t_inc_num )
          the_token = token( t_inc_num, string_index, .true., &
                             source_start )
    exit ! main lexer loop
        end if
        call add_char ( ch )
        need = .true.
        state = num5
      case ( num5 ) ! digit(digit|_)* . digit(digit|_)* E (+|-)? digit(digit|_)*
        select case ( class )
        case ( digit ); call add_char ( ch ); need = .true.
        case ( under ); need = .true.
        case default
          string_index = add_terminal ( t_number )
          the_token = token( t_number, string_index, .true., &
                             source_start )
    exit ! main lexer loop
        end select ! class
      case ( op )
        if ( class == op_char ) then
          call add_char ( ch )
          need = .true.
        else
          call lookup_and_insert ( string_index, found )
          if ( .not. found ) then
            call set_symbol ( string_index, t_unk_op )
          else if ( term_types(symbol(string_index)) == def_op ) then
            the_token = token( symbol(string_index), string_index, .false., &
                               source_start )
    exit ! main lexer loop
          end if
          call error ( unrec_token )
          if ( toggle(lex) ) then
            call dump_1_symbol ( string_index, advance='yes' )
          end if
          state = start
        end if ! class
      case ( pun )
        call add_char ( ch )
        need = .true.
        call lookup_and_insert ( string_index, found )
        if ( .not. found ) then
          call set_symbol ( string_index, t_unk_pun )
        else if ( term_types(symbol(string_index)) == def_pun ) then
          the_token = token( symbol(string_index), string_index, .false., &
                             source_start )
    exit ! main lexer loop
        end if
        call error ( unrec_token )
        if ( toggle(lex) ) then
          call dump_1_symbol ( string_index, advance='yes' )
        end if
        state = start
      case ( contin )
        select case ( class )
        case ( cmt )
          call new_line
          state = start
          need = .true.
        case ( eol_in );           state = start; need = .true.
        case ( spaces );           need = .true.
        case default ! skip to the end-of-line, then announce an error
          do
            call get_char ( ch )
            if ( ch == eol .or. ch == eof ) then
          exit; end if
          end do
          string_index = enter_terminal ( '<aft>', t_aft_cont )
          the_token = token( t_aft_cont, string_index, .false., source_start )
    exit ! main lexer loop
        end select ! class
      case ( sq1 )
        select case ( class )
        case ( eof_in, eol_in )
          call clear_string
          string_index = enter_terminal ( '<inc>', t_inc_str )
          the_token = token( t_inc_str, string_index, .false., source_start )
    exit ! main lexer loop
        case ( quote )
          call add_char ( ch )
          need = .true.
          state = sq2
        case default
          call add_char ( ch )
          need = .true.
        end select ! class
      case ( sq2 )
        if ( class == quote ) then
          need = .true.
          state = sq1
        else
          string_index = add_terminal ( t_string )
          the_token = token( t_string, string_index, .true., source_start )
    exit ! main lexer loop
        end if
      case ( sa1 )
        select case ( class )
        case ( eof_in, eol_in )
          call clear_string
          string_index = enter_terminal ( '<inc>', t_inc_str )
          the_token = token( t_inc_str, string_index, .false., source_start )
    exit ! main lexer loop
        case ( apost )
          call add_char ( ch )
          need = .true.
          state = sa2
        case default
          call add_char ( ch )
          need = .true.
        end select ! class
      case ( sa2 )
        if ( class == apost ) then
          need = .true.
          state = sa1
        else
          string_index = add_terminal ( t_string )
          the_token = token( t_string, string_index, .true., source_start )
    exit ! main lexer loop
        end if
      case ( tog )
        select case ( cap(ch) )
        case ( 'T' ) ! Dump entire token table
          do i = 1, how_many_strings()
            call dump_1_symbol ( i, advance='yes' )
          end do
        case ( 'C' ); toggle(con) = .not. toggle(con)
        case ( 'E' ); toggle(emit) = .not. toggle(emit)
        case ( 'G' ); toggle(gen) = .not. toggle(gen)
        case ( 'L' ); toggle(lex) = .not. toggle(lex)
        case ( 'P' ); toggle(par) = .not. toggle(par)
        case ( 'A' ); toggle(syn) = .not. toggle(syn)
        case ( 'S' ); toggle(tab) = .not. toggle(tab)
        case default
          call output ( source_line, 5 )
          call output (": *** REMARK *** undefined toggle '" )
          call output ( cap(ch) ); call output ( "' in column " )
          call output ( source_column )
          call output ( ' ignored.', advance='yes' )
        end select ! cap(ch)
        need = .true.
        state = start
      end select ! state
    end do ! main lexer loop
    if ( toggle(lex) ) then
      call dump_1_symbol ( the_token%string_index, advance='yes' )
    end if
    return

  contains

    character function CAP ( CH )
      character, intent(in) :: CH
      if ( ch >= 'a' .and. ch <= 'z' ) then
        cap = achar( ichar(ch) + ichar('A') - ichar('a') )
      else
        cap = ch
      end if
      return
    end function

    subroutine ERROR ( CODE )
      integer, intent(in) :: CODE
      call output ( source_line, 5 )
      select case ( code )
      case ( double_under )
        call output ( ': Double underscore in column ' )
        call output ( source_column )
        call output ( ' ignored.', advance='yes' )
      case ( incomplete )
        call output ( ': Incomplete string from column ' )
        call output ( column ); call output ( ' to ' )
        call output ( source_column, advance='yes' )
      case ( unrec_char )
        call output ( ': Unrecognized character in column ' )
        call output ( source_column )
        call output ( ' ignored.', advance='yes' )
      case ( unrec_token )
        call output ( ': Unrecognized token starting in column ' )
        call output ( column )
        call output ( ' and ending in column ' )
        call output ( source_column - 1 )
        call output ( ' ignored.', advance='yes' )
      end select
      call error_walkback ( source_line )
    end subroutine ERROR

  end subroutine LEXER

end module LEXER_M

! $Log$
! Revision 2.5  2000/11/30 00:31:12  vsnyder
! Make [] punctuators instead of operators.
!
! Revision 2.4  2000/11/30 00:23:10  vsnyder
! Implement [] syntax for arrays
!
! Revision 2.3  2000/10/18 18:43:14  vsnyder
! Added starting column number to "Unrecognized token" error message
!
! Revision 2.2  2000/10/18 18:36:58  vsnyder
! Add line number of input to error message.
!
! Revision 2.1  2000/10/11 18:00:54  vsnyder
! Move from lib/cf_parser to lib; remove unused variables; add copyright
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1 2000/07/06 01:43:12 vsnyder
! Initial check-in
!
