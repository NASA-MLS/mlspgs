! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module LEXER_M

! Lexers for MLS.

  use ERROR_HANDLER, only: ERROR_WALKBACK
  use LEXER_CORE, only: NEED, PRINT_SOURCE, TOKEN, Where_t
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: ADD_CHAR, CLEAR_STRING, DISPLAY_STRING, &
                          EOF, EOL, GET_CHAR, INCLUDE_STACK_TOP, &
                          LOOKUP_AND_INSERT, NEW_LINE, SOURCE_COLUMN, &
                          SOURCE_LINE
  use SYMBOL_TABLE, only: ADD_TERMINAL, DUMP_SYMBOL_TABLE, DUMP_1_SYMBOL, &
                          ENTER_TERMINAL, SET_SYMBOL, SYMBOL
  use SYMBOL_TYPES ! everything
  use TOGGLES, only: CON, EMIT, GEN, LEVELS, LEX, PAR, SYN, TAB, TOGGLE
  private

  public :: Lexer, Lex_Signal

  logical, public  :: CapIdentifiers = .true.     ! Capitalize identifiers?

  ! Parameters for character classes
  integer, parameter :: LETTER = 1
  integer, parameter :: DIGIT = 2
  integer, parameter :: PUNCT = 3
  integer, parameter :: OP_CHAR = 4
  integer, parameter :: MINUS = 5       ! So that =- will not be an operator 
  integer, parameter :: DOT = 6
  integer, parameter :: SPACES = 7
  integer, parameter :: EOL_IN = 8
  integer, parameter :: EOF_IN = 9
  integer, parameter :: CMT = 10
  integer, parameter :: MORE = 11
  integer, parameter :: UNDER = 12
  integer, parameter :: CONT = 13
  integer, parameter :: QUOTE = 14
  integer, parameter :: APOST = 15
  integer, parameter :: SHARP = 16
  integer, parameter :: TOG_CH = 17

  ! Parameter for classifying ASCII characters
  integer, parameter :: CLASSES(0:127) = (/ &
! NUL     SOH     STX     ETX     EOT=^D  ENQ     ACK     BEL
  more,   more,   more,   more,   eof_in, more,   more,   more,   & ! 000
! BS      HT      LF      VT      FF      CR      SO      SI
  more,   spaces, eol_in, more,   more,   spaces, more,   more,   & ! 010
! DLE     DC1     DC2     DC3     DC4     NAK     SYN     ETB
  more,   more,   more,   more,   more,   more,   more,   more,   & ! 020
! CAN     EM      SUB     ESC     FS      GS      RS      US
  more,   more,   more,   more,   more,   more,   more,   more,   & ! 030
! SPACE   !       "       #       $       %       &       '
  spaces, op_char,quote,  sharp,  cont,   more,   more,   apost,  & ! 040
! (       )       *       +       ,       -       .       /
  punct,  punct,  op_char,op_char,punct,  minus,  dot,    op_char,& ! 050
! 0       1       2       3       4       5       6       7
  digit,  digit,  digit,  digit,  digit,  digit,  digit,  digit,  & ! 060
! 8       9       :       ;       <       =       >       ?
  digit,  digit,  op_char,cmt,    op_char,op_char,op_char,op_char,& ! 070
! @       A       B       C       D       E       F       G
  tog_ch, letter, letter, letter, letter, letter, letter, letter, & ! 100
! H       I       J       K       L       M       N       O
  letter, letter, letter, letter, letter, letter, letter, letter, & ! 110
! P       Q       R       S       T       U       V       W
  letter, letter, letter, letter, letter, letter, letter, letter, & ! 120
! X       Y       Z       [       \       ]       ^       _
  letter, letter, letter, punct,  op_char,punct,  op_char,under,  & ! 130
! `       a       b       c       d       e       f       g
  more,   letter, letter, letter, letter, letter, letter, letter, & ! 140
! h       i       j       k       l       m       n       o
  letter, letter, letter, letter, letter, letter, letter, letter, & ! 150
! p       q       r       s       t       u       v       w
  letter, letter, letter, letter, letter, letter, letter, letter, & ! 160
! x       y       z       {       |       }       ~       DEL
  letter, letter, letter, more,   cont,   more,   more,   more   /) ! 170

  logical, private, save :: First = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------------------  Lexer  -----
  subroutine LEXER ( THE_TOKEN )
    type(token), intent(out) :: THE_TOKEN

    ! Parameters for ERROR routine's CODE argument
    integer, parameter :: DOUBLE_UNDER = 1
    integer, parameter :: INCOMPLETE = 2
    integer, parameter :: UNREC_CHAR = 3
    integer, parameter :: UNREC_TOKEN = 4

    ! Parameters for states
    integer, parameter :: START = 1
    integer, parameter :: ID = 2
    integer, parameter :: ID2 = 3    ! identifier space+
    integer, parameter :: NUM = 4    ! digit(digit|_)*
    integer, parameter :: NUM2 = 5   ! digit(digit|_)* . digit(digit|_)*
    integer, parameter :: NUM3 = 6   ! digit(digit|_)* . digit(digit|_)* E
    integer, parameter :: NUM4 = 7   ! digit(digit|_)* . digit(digit|_)* E +|-
    integer, parameter :: NUM5 = 8   ! digit(digit|_)* . digit(digit|_)* E (+|-)?
                                     ! digit(digit|_)*
    integer, parameter :: OP = 9
    integer, parameter :: PUN = 10
    integer, parameter :: CONTIN = 11
    integer, parameter :: CONTIN2 = 12
    integer, parameter :: SQ1 = 13   ! First string-started-by-quote state
    integer, parameter :: SQ2 = 14   ! Second string-started-by-quote state
    integer, parameter :: SA1 = 15   ! First string-started-by-apost state
    integer, parameter :: SA2 = 16   ! Second string-started-by-apost state
    integer, parameter :: SH1 = 17   ! #
    integer, parameter :: TOG = 18   ! @...

    character, save :: CH          ! Current character
    integer :: CLASS               ! Classification of current character
    
    integer :: COLUMN              ! Character index of start of token in
                                   ! CHAR_TABLE
    logical :: FOUND               ! Argument for LOOKUP_AND_INSERT
    integer :: I                   ! Loop inductor
    integer :: SOURCE_START        ! Source reference for start of token
                                   ! ( 256*line + column )
    integer :: STATE               ! Current state of lexer DFA
    integer :: STRING_INDEX        ! Index of token in string table

    integer :: File                ! String index from include stack

    if ( first .and. toggle(lex) .and. levels(lex) > 1 ) then
      first = .false.
      call output ( 'Term_Types:', advance='yes' )
      do i = lbound(term_types,1), ubound(term_types,1), 7
        do column = i, min(ubound(term_types,1),i+6)
          call output ( term_types(column), format='(i4)' )
        end do
        call output ( '', advance='yes' )
      end do
    end if

    call include_stack_top ( file )
    state = start
    do
      if ( need ) then; call get_char ( ch ); need = .false.; end if
      ! Classify next character
      class = iachar(ch)
      if ( class >= lbound(classes,1) .and. class <= ubound(classes,1) ) then
        class = classes(class)
      else
        class = more
      end if
      if ( levels(lex) > 0 ) then
        call output ( state, before='State = ' )
        call output ( class, before=' Class = ', advance='yes' )
      end if

!     Recognize a token according to the following DFA.  Column labels are
!     Classes, row labels are DFA states.  Table elements are DFA states or
!     actions to be taken when a token is recognized.
!
!         letter digit opchar  minus  punct  under  cmt eolin eofin space
!  start    id    num    op      op    pun    id    cmt    6    1   start
!  id       id    id     11      11     11    id    11    11   11    11
!  id2      11    11     11      11     11    11    11    11   11   id2
!  num       3    num     3       3      3    num    3     3    3     3
!  num2   e=num3 num2     3       3      3   num2    3     3    3     3
!         else 3
!  num3     10   num5  +-=num4 +-=num4  10    10    10    10   10    10
!                      else 10 else 10
!  num4     10   num5    10      10     10    10    10    10   10    10
!  num5      3   num5     3       3      3   num5    3     3    3     3
!  op        4     4     op      op      4     4     4     4    4     4
!  pun       5     5      5       5      5     5     5     5    5     5
!  cmt      cmt   cmt    cmt     cmt    cmt   cmt   cmt  start  1    cmt
!  cont      7     7      7       7      7     7   cont2 cont2  1   cont
!  cont2   start start  start   start  start start cont2 cont2  1  cont2
!  sq1      sq1   sq1    sq1     sq1    sq1   sq1   sq1    9    9    sq1
!  sq2       8     8      8       8      8     8     8     8    8     8
!  sa1      sa1   sa1    sa1     sa1    sa1   sa1   sa1    9    9    sa1
!  sa2       8     8      8       8      8     8     8     8    8     8
!  sh1      sh1    4      4       4      4     4     4     4    4     4
!
!           cont  more quote apost  dot  sharp
!  start    cont start  sq1   sa1    4    sh1
!  id       11    11    11    11    11    11
!  id2      11    11    11    11    11    11
!  num       3     3     3     3   num2    3
!  num2      3     3     3     3     3     3
!  num3     10    10    10    10    10     3
!  num4     10    10    10    10    10     3
!  num5      3     3     3     3     3     3
!  op        4     4     4     4     4     4
!  pun       5     5     5     5     5     5
!  cmt      cmt   cmt   cmt   cmt   cmt   cmt
!  cont      7     7     7     7     7     7
!  cont2   cont  start start start start start
!  sq1      sq1   sq1   sq2   sq1   sq1   sq1
!  sq2       8     8    sq1    8     8     8
!  sa1      sa1   sa1    8    sa2    8     8
!  sa2       8     8     8    sa1    8     8
!  sh1       4     4     4     4     4     4

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
!    11.   if the next character is not "(" and its class is not letter,
!          end-of-line, end_of_file, or character, define identifier, 
!          else do action 2.

      select case ( state )
      case ( start )
        column = source_column
        source_start = 256 * source_line + source_column
        select case ( class )
        case ( letter, under );    state = id
        case ( digit );            state = num
        case ( punct );            state = pun
        case ( op_char );          state = op
        case ( dot, minus )
          need = .true.
          call add_char ( ch )
          call lookup_and_insert ( string_index, found ) ! will be found
                                                         ! if tables OK
          the_token = token( symbol(string_index), string_index, .false., &
                             where_t(file, source_start) )
    exit ! main lexer loop
        case ( cmt )
          call new_line
          state = start
          need = .false.
          ch = EOL
        case ( spaces );           need = .true.
        case ( cont );             state = contin; need = .true.
        case ( eof_in )
          string_index = enter_terminal ( '<eof>', t_end_of_input )
          the_token = token( t_end_of_input, string_index, .false., &
                             where_t(file, source_start) )
          call set_symbol ( string_index, t_end_of_input )
    exit ! main lexer loop
        case ( eol_in )
          string_index = enter_terminal ( '<eos>', t_end_of_stmt )
          the_token = token( t_end_of_stmt, string_index, .false., &
                             where_t(file, source_start) )
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
        case ( sharp );
          call add_char ( ch )
          state = sh1; need = .true.
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
          the_token = token( t_aft_cont, string_index, .false., &
                           & where_t(file, source_start) )
          need = .false.
        end select ! class
      case ( id )
        select case ( class )
        case ( letter, digit, under )
          if ( capIdentifiers ) ch = cap(ch)
          call add_char ( ch )
          need = .true.
        case ( spaces )
          need = .true.
          state = id2
        case default
          call check_reserved_word
    exit ! main lexer loop
        end select ! class
      case ( id2 ) ! gobble up spaces after an identifier
        if ( class == spaces ) then
          need = .true.
        else
          call check_reserved_word
    exit ! main lexer loop
        end if
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
                             where_t(file, source_start) )
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
                               where_t(file, source_start) )
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
                             where_t(file, source_start) )
    exit ! main lexer loop
        end if
      case ( num4 ) ! digit(digit|_)* . digit(digit|_)* E +|-
        if ( class /= digit ) then
          string_index = add_terminal ( t_inc_num )
          the_token = token( t_inc_num, string_index, .true., &
                             where_t(file, source_start) )
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
                             where_t(file, source_start) )
    exit ! main lexer loop
        end select ! class
      case ( op )
        if ( class == op_char ) then
          call add_char ( ch )
          need = .true.
        else
          call lookup_op
          if ( found ) then
    exit ! main lexer loop
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
                             where_t(file, source_start) )
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
          state = contin2
          need = .true.
        case ( eol_in );           state = contin2; need = .true.
        case ( spaces );           need = .true.
        case default ! skip to the end-of-line, then announce an error
          do
            call get_char ( ch )
            if ( ch == eol .or. ch == eof ) then
          exit; end if
          end do
          string_index = enter_terminal ( '<aft>', t_aft_cont )
          the_token = token( t_aft_cont, string_index, .false., &
                           & where_t(file, source_start) )
    exit ! main lexer loop
        end select ! class
      case ( contin2 )
        select case ( class )
        case ( cmt )
          call new_line
          need = .true.
        case ( cont );             state = contin; need = .true.
        case ( eol_in, spaces );   need = .true.
        case default;              state = start
        end select ! class
      case ( sq1 )
        select case ( class )
        case ( eof_in, eol_in )
          call clear_string
          string_index = enter_terminal ( '<inc>', t_inc_str )
          the_token = token( t_inc_str, string_index, .false., &
                           & where_t(file, source_start) )
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
          the_token = token( t_string, string_index, .true., &
                           & where_t(file, source_start) )
    exit ! main lexer loop
        end if
      case ( sa1 )
        select case ( class )
        case ( eof_in, eol_in )
          call clear_string
          string_index = enter_terminal ( '<inc>', t_inc_str )
          the_token = token( t_inc_str, string_index, .false., &
                           & where_t(file, source_start) )
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
          the_token = token( t_string, string_index, .true., &
                           & where_t(file, source_start) )
    exit ! main lexer loop
        end if
      case ( sh1 )   ! #<letter>*
        if ( class == letter ) then
          need = .true.
          call add_char ( ch )
        else
          call lookup_op
          if ( found ) &
    exit ! main lexer loop
          state = start
        end if ! class
      case ( tog )
        select case ( cap(ch) )
        case ( 'T' ) ! Dump entire token table
          call dump_symbol_table
        case ( 'C' ); toggle(con) = .not. toggle(con)
        case ( 'E' ); toggle(emit) = .not. toggle(emit)
        case ( 'G' ); toggle(gen) = .not. toggle(gen)
        case ( 'L' ); toggle(lex) = .not. toggle(lex)
        case ( 'P' ); toggle(par) = .not. toggle(par)
        case ( 'A' ); toggle(syn) = .not. toggle(syn)
        case ( 'S' ); toggle(tab) = .not. toggle(tab)
        case default
          call include_stack_top ( file )
          call print_source ( source_line+256+source_column, file=file )
          call output (": *** REMARK *** undefined toggle '" )
          call output ( cap(ch) )
          call output ( '" ignored.', advance='yes' )
        end select ! cap(ch)
        need = .true.
        state = start
      end select ! state
    end do ! main lexer loop
    if ( toggle(lex) ) then
      call dump_1_symbol ( the_token%string_index, advance='yes' )
    end if

  contains

    subroutine Check_Reserved_Word
      ! The current token looks like an identifier. The next character is
      ! not a space.  If it's not (, @, a letter, an end-of-line, end-of-file,
      ! or comment, the symbol cannot be a reserved word
      logical :: Not_Reserved
      not_reserved = ch /= '(' .and. class /= tog_ch .and. class /= letter .and. &
                     class /= eol_in .and. class /= eof_in .and. class /= cmt
      string_index = add_terminal ( t_identifier )
      if ( toggle(lex) .and. levels(lex) > 1 ) then
        call output ( 'In Check_Reserved_Word with ch = "' )
        call output ( ch )
        call output ( class, before='" and class = ' )
        call display_string ( string_index, before=', ' )
      end if
      if ( not_reserved ) then
        if ( toggle(lex) .and. levels(lex) > 1 ) &
          & call output ( ' cannot be reserved', advance='yes' )
        the_token = token( t_identifier, string_index, .true., &
                           where_t(file, source_start) )
      else
        the_token = token( symbol(string_index), string_index, &
                           term_types(symbol(string_index)) /= res_word, &
                           where_t(file, source_start) )
      if ( toggle(lex) .and. levels(lex) > 1 ) &
        & call output ( ' is' // &
          & trim(merge(' not','    ',term_types(symbol(string_index)) /= res_word)) // &
          & ' reserved', advance='yes' )
      end if
    end subroutine Check_Reserved_Word

    character function CAP ( CH )
      character, intent(in) :: CH
      if ( ch >= 'a' .and. ch <= 'z' ) then
        cap = achar( ichar(ch) + ichar('A') - ichar('a') )
      else
        cap = ch
      end if
    end function

    subroutine ERROR ( CODE )
      integer, intent(in) :: CODE
      call include_stack_top ( file )
      if ( file /= 0 ) then
        call display_string ( file, before=' in ' )
        call output ( source_line, before=' at line ' )
      else
        call output ( source_line, 5 )
      end if
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
        call output ( ': Unrecognized character ' )
        call output ( ch )
        call output ( ' in column ' )
        call output ( source_column )
        call output ( ' ignored.', advance='yes' )
      case ( unrec_token )
        call output ( ': Unrecognized token ' )
        call display_string ( string_index )
        call output ( ' starting in column ' )
        call output ( column )
        call output ( ' and ending in column ' )
        call output ( source_column - 1 )
        call output ( ' ignored.', advance='yes' )
      end select
      call error_walkback ( source_line )
    end subroutine ERROR

    subroutine Lookup_Op
      call lookup_and_insert ( string_index, found )
      if ( .not. found ) then
        call set_symbol ( string_index, t_unk_op )
      else
        found = term_types(symbol(string_index)) == def_op
        if ( found ) &
        the_token = token( symbol(string_index), string_index, .false., &
                           where_t(file, source_start) )
      end if
      if ( .not. found ) then
        call error ( unrec_token )
        if ( toggle(lex) ) then
          call dump_1_symbol ( string_index, advance='yes' )
        end if
      end if
    end subroutine Lookup_Op

  end subroutine LEXER

  ! -------------------------------------------------  Lex_Signal  -----
  subroutine Lex_Signal ( Buf, Where, The_Token )

  ! Get the next piece of a signal string from Buf.  Blank or end-of-string
  ! terminates the signal string.

    character(len=*), intent(in) :: Buf      ! The signal string, sans quotes
    integer, intent(inout) :: Where          ! Last analyzed character in Buf
    type(token), intent(out) :: The_Token    ! The token

    character :: Ch                          ! Next character from Buf
    integer :: Class                         ! Class of next character
    logical :: Found                         ! "Found the symbol by lookup"
    integer :: Source_Start                  ! Column of beginning of token
    integer :: State                         ! State of Lexer
    integer :: String_Index                  ! Where in the string table?

    ! Parameters for states
    integer, parameter :: Start = 1
    integer, parameter :: Id = 2

!     Recognize a token according to the following DFA.  Column labels are
!     Classes, row labels are DFA states.  Table elements are DFA states or
!     actions to be taken when a token is recognized.
!
!                   dot
!         letter  punct
!          digit opchar End  more
! Start     Id      1    3  Error
! Id        Id      2    2    2
!
! Actions:  1: Output an operator or punctuator.  Consume the character.
!           2: Output an identifier.  Don't consume the character.
!           3: Output an end-of-input signal.  Don't consume the "character."

    state = start
    do
      where = where + 1
      if ( where > len(buf) ) then
        ch = eof
        class = eof_in
      else
        ch = buf(where:where)
        class = classes(iachar(ch))
      end if
      select case ( state )
      case ( start )
        select case ( class )
        case ( eof_in, spaces )
          string_index = enter_terminal ( '<eof>', t_end_of_input )
          the_token = token( t_end_of_input, string_index, .false., &
                           & where_t(0, where) )
          call set_symbol ( string_index, t_end_of_input )
    exit
        case ( letter, digit )
          call add_char ( ch )
          source_start = where
          state = id
        case ( dot, minus, op_char, punct )
          call add_char ( ch )
          call lookup_and_insert ( string_index, found )
          if ( .not. found ) then
            call set_symbol ( string_index, t_unk_op )
          else
            the_token = token( symbol(string_index), string_index, .false., &
                             & where_t(0, where) )
          end if
    exit ! main lexer loop
        case default
          call lookup_and_insert ( string_index, found )
          call set_symbol ( string_index, t_unk_op )
          the_token = token( t_unk_op, string_index, .false., &
                           & where_t(0, where) )
    exit ! main lexer loop
        end select
      case ( id )
        if ( class == letter .or. class == digit ) then
          call add_char ( ch )
        else
          where = where - 1
          string_index = add_terminal ( t_identifier )
          the_token = token( t_identifier, string_index, .true., &
                           & where_t(0, source_start) )
    exit ! main lexer loop
        end if
      end select
    end do

    if ( toggle(lex) ) then
      call dump_1_symbol ( the_token%string_index, advance='yes' )
    end if

  end subroutine Lex_Signal

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LEXER_M

! $Log$
! Revision 2.29  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.28  2013/12/12 02:03:22  vsnyder
! Add some dumps to Check_Reserved_Word
!
! Revision 2.27  2013/11/26 22:50:14  vsnyder
! Add Check_Reserved_Word.  Don't call a word reserved if it's followed by a
! letter, left parenthesis, end-of-line, end-of-file, or semicolon.
!
! Revision 2.26  2013/10/16 01:13:50  vsnyder
! New 'minus' class for hyphen so =- will be two operators
!
! Revision 2.25  2013/10/02 01:33:46  vsnyder
! Make ! and ? operator characters
!
! Revision 2.24  2013/09/24 23:08:31  vsnyder
! Add #include
!
! Revision 2.23  2011/04/18 19:33:26  vsnyder
! Add support for relational operators and boolean-valued expressions
!
! Revision 2.22  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.21  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.20  2004/05/28 23:12:20  vsnyder
! Add power (^) operator
!
! Revision 2.19  2004/01/16 23:49:32  vsnyder
! Add backslash for 'into' operator
!
! Revision 2.18  2003/01/29 00:49:38  vsnyder
! Delete USE for unused entity
!
! Revision 2.17  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.16  2002/02/13 23:39:05  vsnyder
! Handle characters outside 0..127 range sensibly
!
! Revision 2.15  2001/06/21 20:08:10  vsnyder
! Remove obsolete comments, allow blank continuations
!
! Revision 2.14  2001/06/21 18:46:57  vsnyder
! Allow multiple comments and blank lines between continuations
!
! Revision 2.13  2001/03/15 18:40:39  vsnyder
! Don't generate numeric tokens in lex_signal; call them identifiers.
!
! Revision 2.12  2001/03/14 19:15:32  vsnyder
! Add comments that describe the DFA used in lex_signal.
!
! Revision 2.11  2001/03/14 18:59:54  vsnyder
! Added public variable CapIdentifiers to control whether identifiers are
! capitalized.
!
! Revision 2.10  2001/03/14 02:06:31  vsnyder
! Added lex_signal
!
! Revision 2.9  2001/03/06 03:02:25  vsnyder
! Allow comments between continuations -- maybe got it right this time
!
! Revision 2.8  2001/03/05 23:18:43  vsnyder
! Allow comments between continuations
!
! Revision 2.7  2000/11/30 20:18:47  vsnyder
! Added <: :< and <:< operators.
!
! Revision 2.6  2000/11/30 20:03:57  vsnyder
! Correct handling of end-of-comment; improve error messages.
!
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
