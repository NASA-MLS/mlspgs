! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SYMBOL_TYPES

!  Type definitions for the symbol table.

!  The terminal symbols of the grammar.  The corresponding constants
!  must be typed into the grammar as integer literal values.

  use STRING_TABLE, only: ADD_CHAR

  implicit NONE
  public

  ! Terminal symbol class numbers.  These depend on the language.
  integer, parameter :: T_NULL = 0                ! MUST be present and first
  integer, parameter :: T_LEFT_PARENTHESIS = 1
  integer, parameter :: T_RIGHT_PARENTHESIS = 2
  integer, parameter :: T_LEFT_BRACKET = 3
  integer, parameter :: T_RIGHT_BRACKET = 4
  integer, parameter :: T_PLUS = 5
  integer, parameter :: T_MINUS = 6
  integer, parameter :: T_STAR = 7
  integer, parameter :: T_SLASH = 8
  integer, parameter :: T_DOT = 9
  integer, parameter :: T_COLON = 10
  integer, parameter :: T_EQUAL = 11
  integer, parameter :: T_COMMA = 12
  integer, parameter :: T_BEGIN = 13              ! BEGIN
  integer, parameter :: T_END = 14                ! END
  integer, parameter :: T_AND = 15                ! AND
  integer, parameter :: T_OR = 16                 ! OR
  integer, parameter :: T_END_OF_INPUT = 17       ! <EOF>
  integer, parameter :: T_END_OF_STMT = 18        ! <EOS>
  integer, parameter :: T_IDENTIFIER = 19         ! <IDENTIFIER>
  integer, parameter :: T_NUMBER = 20             ! <NUMBER>
  integer, parameter :: T_STRING = 21             ! <STRING>
  integer, parameter :: T_UNK_OP = 22             ! unknown operator
  integer, parameter :: T_UNK_PUN = 23            ! unknown punctuator
  integer, parameter :: T_UNK_CH = 24             ! unknown character
  integer, parameter :: T_INC_NUM = 25            ! incomplete number
  integer, parameter :: T_INC_STR = 26            ! incomplete string
  integer, parameter :: T_AFT_CONT = 27           ! junk after continuation

! The parameters T_LAST_TERMINAL, MIN_PSEUDO, MAX_PSEUDO and CASELESS_LOOK
! MUST be defined.
  integer, parameter :: T_LAST_TERMINAL = T_AFT_CONT

  integer, parameter :: MIN_PSEUDO = T_IDENTIFIER
  integer, parameter :: MAX_PSEUDO = T_STRING

  integer, private :: I ! Only a loop inductor in the next declaration
! Indicates which terminals use caseless lookup
  logical, parameter :: CASELESS_LOOK(t_null:t_last_terminal) = &
    (/ ( .false., i = t_null, min_pseudo-1 ), &
       ( .true., i = min_pseudo, max_pseudo ), &
       ( .false., i = max_pseudo+1, t_last_terminal ) /)

! Terminal symbol types.  These are used for debug printing.

  integer, parameter :: RES_WORD = 1              ! Reserved word
  integer, parameter :: IDENT = 2                 ! Identifier
  integer, parameter :: NUMCON = 3                ! Numeric constant
  integer, parameter :: STRING = 4                ! String in quotes
  integer, parameter :: UNIT = 5                  ! A unit name, e.g. mb
  integer, parameter :: DEF_OP = 6                ! Defined operator
  integer, parameter :: UNK_OP = 7                ! Unefined operator
  integer, parameter :: DEF_PUN = 8               ! Defined punctuator
  integer, parameter :: UNK_PUN = 9               ! Undefined punctuator
  integer, parameter :: INC_NUM = 10              ! Incomplete number
  integer, parameter :: INC_STR = 11              ! Incomplete string
  integer, parameter :: UNK_CH = 12               ! Unknown character
  integer, parameter :: AFT_CONT = 13             ! After continuation
  integer, parameter :: OBJECT = 14               ! None of the above

  integer, parameter :: FIRST_TYPE = RES_WORD
  integer, parameter :: LAST_TYPE = OBJECT

! The array TERM_TYPES gives the terminal type of each class of terminal.
! It must be defined.
  integer, parameter :: TERM_TYPES(t_null: t_last_terminal) = &
  !  t_null    (         )         [         ]         +         -
  (/ object,   def_pun,  def_pun,  def_pun,  def_pun,  def_op,   def_op,   &
  !  *         /         .         :         =         ,         begin
     def_op,   def_op,   def_op,   def_op,   def_op,   def_pun,  res_word, &
  !  end       and       or        <eof>     <eos>     <ident>   <numcon>
     res_word, res_word, res_word, object,   object,   ident,    numcon,   &
  !  <string>  unk_op    unk_pun     unk_ch    inc_num   inc_str   junk
     string,   unk_op,   unk_pun,    unk_ch,   inc_num,  inc_str,  aft_cont /)

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine INIT_TERMINAL ( TERMINAL )
  ! Put the text of a terminal symbol into the string table. It isn't
  ! necessary to have more than zero characters of text for things that
  ! aren't terminals (e.g. t_null and errors), or for things that are
  ! pseudo-terminals (e.g. identifier, integer).
    integer, intent(in) :: TERMINAL
    select case ( terminal )
    case ( t_null );              call add_char ( '' )
    case ( t_left_parenthesis );  call add_char ( '(' )
    case ( t_right_parenthesis ); call add_char ( ')' )
    case ( t_left_bracket );      call add_char ( '[' )
    case ( t_right_bracket );     call add_char ( ']' )
    case ( t_plus );              call add_char ( '+' )
    case ( t_minus );             call add_char ( '-' )
    case ( t_star );              call add_char ( '*' )
    case ( t_slash );             call add_char ( '/' )
    case ( t_dot );               call add_char ( '.' )
    case ( t_colon   );           call add_char ( ':' )
    case ( t_equal );             call add_char ( '=' )
    case ( t_comma );             call add_char ( ',' )
    case ( t_begin );             call add_char ( 'BEGIN' )
    case ( t_end );               call add_char ( 'END' )
    case ( t_and );               call add_char ( 'AND' )
    case ( t_or );                call add_char ( 'OR' )
    case ( t_end_of_input );      call add_char ( '<eof>' )
    case ( t_end_of_stmt );       call add_char ( '<eos>' )
    case ( t_identifier );        call add_char ( '<identifier>' )
    case ( t_number );            call add_char ( '<number>' )
    case ( t_string );            call add_char ( '<string>' )
    case ( t_unk_op );            call add_char ( '<unk_op>' )
    case ( t_unk_pun );           call add_char ( '<unk_pun>' )
    case ( t_unk_ch );            call add_char ( '<unk_ch>' )
    case ( t_inc_num );           call add_char ( '<inc num>' )
    case ( t_inc_str );           call add_char ( '<inc str>' )
    case ( t_aft_cont );          call add_char ( '<aft>' )
    case default
      write ( *, * ) 'SYMBOL_TYPES%INIT_TERMINAL-E- No initializer for &
                     &terminal symbol with index ', terminal
      stop
    end select
    return
  end subroutine INIT_TERMINAL

  subroutine INIT_TYPE ( TYPE )
  ! Put the text of a terminal symbol type into the string table.  These
  ! are used only for debugging output.
    integer, intent(in) :: TYPE
    select case ( type )
    case ( res_word );  call add_char ( '<reserved word>' )
    case ( ident );     call add_char ( '<identifier>' )
    case ( numcon );    call add_char ( '<number>' )
    case ( string );    call add_char ( '<string>' )
    case ( unit );      call add_char ( '<unit name>' )
    case ( def_op );    call add_char ( '<defined operator>' )
    case ( unk_op );    call add_char ( '<undefined operator>' )
    case ( def_pun );   call add_char ( '<defined punctuator>' )
    case ( unk_pun );   call add_char ( '<undefined punctuator>' )
    case ( unk_ch );    call add_char ( '<unrecognized character>' )
    case ( inc_num );   call add_char ( '<incomplete number>' )
    case ( inc_str );   call add_char ( '<incomplete string>' )
    case ( aft_cont );  call add_char ( '<after continuation>' )
    case ( object );    call add_char ( '<object>' )
    case default
      write ( *, * ) 'SYMBOL_TYPES%INIT_TYPE-E- No initializer for &
                     &type with index ', type
      stop
    end select
    return
  end subroutine INIT_TYPE
end module SYMBOL_TYPES

! $Log$
! Revision 2.3  2000/11/30 00:31:12  vsnyder
! Make [] punctuators instead of operators.
!
! Revision 2.2  2000/11/30 00:23:10  vsnyder
! Implement [] syntax for arrays
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
