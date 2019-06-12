! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SYMBOL_TYPES

!  Type definitions for the symbol table.


  use STRING_TABLE, only: ADD_CHAR

  implicit NONE
  public

  ! Terminal symbol identifiers.  These depend on the language.

  ! The named constants declared here are related to the terminal symbols of
  ! the grammar using the table before the productions in l2cf.grm.

  ! If this list changes, you need to change the array TERM_TYPES below.

  integer, parameter :: T_NULL             = 0    ! MUST be present and first
  integer, parameter :: T_LEFT_PARENTHESIS = T_NULL + 1
  integer, parameter :: T_RIGHT_PARENTHESIS = T_LEFT_PARENTHESIS + 1
  integer, parameter :: T_LEFT_BRACKET =     T_RIGHT_PARENTHESIS + 1
  integer, parameter :: T_RIGHT_BRACKET =    T_LEFT_BRACKET + 1
  integer, parameter :: T_PLUS =             T_RIGHT_BRACKET + 1
  integer, parameter :: T_MINUS =            T_PLUS + 1
  integer, parameter :: T_STAR =             T_MINUS + 1
  integer, parameter :: T_SLASH =            T_STAR + 1
  integer, parameter :: T_ASSIGN =           T_SLASH + 1
  integer, parameter :: T_BACKSLASH =        T_ASSIGN + 1
  integer, parameter :: T_BANG =             T_BACKSLASH + 1
  integer, parameter :: T_COND =             T_BANG + 1
  integer, parameter :: T_DOT =              T_COND + 1
  integer, parameter :: T_COLON =            T_DOT + 1
  integer, parameter :: T_COLON_LESS =       T_COLON + 1
  integer, parameter :: T_LESS_COLON =       T_COLON_LESS + 1
  integer, parameter :: T_LESS_COLON_LESS =  T_LESS_COLON + 1
  integer, parameter :: T_EQUAL =            T_LESS_COLON_LESS + 1
  integer, parameter :: T_EQUAL_EQUAL =      T_EQUAL + 1
  integer, parameter :: T_NOT_EQUAL =        T_EQUAL_EQUAL + 1
  integer, parameter :: T_LESS =             T_NOT_EQUAL + 1
  integer, parameter :: T_LESS_EQ =          T_LESS + 1
  integer, parameter :: T_GREATER =          T_LESS_EQ + 1
  integer, parameter :: T_GREATER_EQ =       T_GREATER + 1
  integer, parameter :: T_COMMA =            T_GREATER_EQ + 1
  integer, parameter :: T_HAT =              T_COMMA + 1
  integer, parameter :: T_BEGIN =            T_HAT + 1           ! BEGIN
  integer, parameter :: T_CYCLE =            T_BEGIN + 1         ! CYCLE
  integer, parameter :: T_DO =               T_CYCLE + 1         ! DO
  integer, parameter :: T_END =              T_DO + 1            ! END
  integer, parameter :: T_EXIT =             T_END + 1           ! EXIT
  integer, parameter :: T_AND =              T_EXIT + 1          ! AND
  integer, parameter :: T_OR =               T_AND + 1           ! OR
  integer, parameter :: T_NOT =              T_OR + 1            ! NOT
  integer, parameter :: T_CASE =             T_NOT + 1           ! CASE
  integer, parameter :: T_DEFAULT =          T_CASE + 1          ! DEFAULT
  integer, parameter :: T_ELSE =             T_DEFAULT + 1       ! ELSE
  integer, parameter :: T_IF =               T_ELSE + 1          ! IF
  integer, parameter :: T_SELECT =           T_IF + 1            ! SELECT
  integer, parameter :: T_THEN =             T_SELECT + 1        ! THEN
  integer, parameter :: T_WHILE =            T_THEN + 1          ! WHILE
  integer, parameter :: T_END_OF_INPUT =     T_WHILE + 1         ! <EOF>
  integer, parameter :: T_END_OF_STMT =      T_END_OF_INPUT + 1  ! <EOS>
  ! T_IDENTIFIER, T_NUMBER, T_STRING must be consecutive
  integer, parameter :: T_IDENTIFIER =       T_END_OF_STMT + 1   ! <IDENTIFIER>
  integer, parameter :: T_NUMBER =           T_IDENTIFIER + 1    ! <NUMBER>
  integer, parameter :: T_STRING =           T_NUMBER + 1        ! <STRING>
  integer, parameter :: T_INCLUDE =          T_STRING + 1        ! #include
  integer, parameter :: T_UNK_OP =           T_INCLUDE + 1       ! unknown operator
  integer, parameter :: T_UNK_PUN =          T_UNK_OP + 1        ! unknown punctuator
  integer, parameter :: T_UNK_CH =           T_UNK_PUN + 1       ! unknown character
  integer, parameter :: T_INC_NUM =          T_UNK_CH + 1        ! incomplete number
  integer, parameter :: T_INC_STR =          T_INC_NUM + 1       ! incomplete string
  integer, parameter :: T_AFT_CONT =         T_INC_STR + 1       ! junk after continuation

! The parameters T_LAST_TERMINAL, MIN_PSEUDO, MAX_PSEUDO and CASELESS_LOOK
! MUST be defined.
  integer, parameter :: T_LAST_TERMINAL = T_AFT_CONT

  integer, parameter :: MIN_PSEUDO = T_IDENTIFIER
  integer, parameter :: MAX_PSEUDO = T_STRING

  integer, private :: I ! Only a loop inductor in the next declaration
! Indicates which terminals use caseless lookup
  logical, parameter :: CASELESS_LOOK(t_null:t_last_terminal) = &
    (/ ( .false., i = t_null, min_pseudo-1 ), &
       ( .true., i = min_pseudo, t_string-1 ), &
       .false., & ! t_string is case sensitive
       ( .true., i = t_string+1, max_pseudo ), &
       ( .false., i = max_pseudo+1, t_last_terminal ) /)

! Terminal symbol types.

  integer, parameter :: RES_WORD = 1              ! Reserved word
  integer, parameter :: IDENT =    RES_WORD + 1   ! Identifier
  integer, parameter :: NUMCON =   IDENT + 1      ! Numeric constant
  integer, parameter :: STRING =   NUMCON + 1     ! String in quotes
  integer, parameter :: UNIT =     STRING + 1     ! A unit name, e.g. mb
  integer, parameter :: DEF_OP =   UNIT + 1       ! Defined operator
  integer, parameter :: UNK_OP =   DEF_OP + 1     ! Unefined operator
  integer, parameter :: DEF_PUN =  UNK_OP + 1     ! Defined punctuator
  integer, parameter :: UNK_PUN =  DEF_PUN + 1    ! Undefined punctuator
  integer, parameter :: INC_NUM =  UNK_PUN + 1    ! Incomplete number
  integer, parameter :: INC_STR =  INC_NUM + 1    ! Incomplete string
  integer, parameter :: UNK_CH =   INC_STR + 1    ! Unknown character
  integer, parameter :: AFT_CONT = UNK_CH + 1     ! After continuation
  integer, parameter :: TREENODE = AFT_CONT + 1   ! Tree node
  integer, parameter :: OBJECT =   TREENODE + 1   ! None of the above

  integer, parameter :: FIRST_TYPE = RES_WORD
  integer, parameter :: LAST_TYPE = OBJECT

! The array TERM_TYPES gives the terminal type of each class of terminal. It
! must be defined.  It is indexed by the terminal symbol names declared
! above.  The lists MUST correspond.  It might be more reliable to use a DATA
! statement, and make TERM_TYPES a protected variable instead of a named
! constant.  If using DATA statements, maybe one for each element, you would
! need to be careful to define every element of the array.  Pay your money
! and take your choice.
  integer, parameter :: TERM_TYPES(t_null: t_last_terminal) = &
  !                                          
  (/ object,   &  ! t_null     T_NULL
     def_pun,  &  ! (          T_LEFT_PARENTHESIS
     def_pun,  &  ! )          T_RIGHT_PARENTHESIS
     def_pun,  &  ! [          T_LEFT_BRACKET
     def_pun,  &  ! ]          T_RIGHT_BRACKET
     def_op,   &  ! +          T_PLUS
     def_op,   &  ! -          T_MINUS
     def_op,   &  ! *          T_STAR
     def_op,   &  ! /          T_SLASH
     def_op,   &  ! :=         T_ASSIGN
     def_op,   &  ! \          T_BACKSLASH
     def_op,   &  ! !          T_BANG
     def_op,   &  ! ?          T_COND
     def_op,   &  ! .          T_DOT
     def_op,   &  ! :          T_COLON
     def_op,   &  ! :<         T_COLON_LESS
     def_op,   &  ! <:         T_LESS_COLON
     def_op,   &  ! <:<        T_LESS_COLON_LESS
     def_op,   &  ! =          T_EQUAL
     def_op,   &  ! ==         T_EQUAL_EQUAL
     def_op,   &  ! /=         T_NOT_EQUAL
     def_op,   &  ! <          T_LESS
     def_op,   &  ! <=         T_LESS_EQ
     def_op,   &  ! >          T_GREATER
     def_op,   &  ! >=         T_GREATER_EQ
     def_pun,  &  ! ,          T_COMMA
     def_op,   &  ! ^          T_HAT
     res_word, &  ! begin      T_BEGIN
     res_word, &  ! cycle      T_CYCLE
     res_word, &  ! do         T_DO
     res_word, &  ! end        T_END
     res_word, &  ! exit       T_EXIT
     res_word, &  ! and        T_AND
     res_word, &  ! or         T_OR
     res_word, &  ! not        T_NOT
     res_word, &  ! case       T_CASE
     res_word, &  ! default    T_DEFAULT
     res_word, &  ! else       T_ELSE
     res_word, &  ! if         T_IF
     res_word, &  ! select     T_SELECT
     res_word, &  ! then       T_THEN
     res_word, &  ! while      T_WHILE
     object,   &  ! <eof>      T_END_OF_INPUT
     object,   &  ! <eos>      T_END_OF_STMT
     ident,    &  ! <ident>    T_IDENTIFIER
     numcon,   &  ! <numcon>   T_NUMBER
     string,   &  ! <string>   T_STRING
     def_op,   &  ! include    T_INCLUDE
     unk_op,   &  ! unk_op     T_UNK_OP
     unk_pun,  &  ! unk_pun    T_UNK_PUN
     unk_ch,   &  ! unk_ch     T_UNK_CH
     inc_num,  &  ! inc_num    T_INC_NUM
     inc_str,  &  ! inc_str    T_INC_STR
     aft_cont  /) ! junk       T_AFT_CONT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine INIT_TERMINAL ( TERMINAL )
  ! Put the text of a terminal symbol into the string table. It isn't
  ! necessary to have more than zero characters of text for things that
  ! aren't terminals (e.g. t_null and errors), or for things that are
  ! pseudo-terminals (e.g. identifier, integer).  Texts consisting of
  ! "< description >" cannot possibly be created by the lexer
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
    case ( t_assign );            call add_char ( ':=' )
    case ( t_backslash );         call add_char ( '\' )
    case ( t_bang );              call add_char ( '!' )
    case ( t_cond );              call add_char ( '?' )
    case ( t_dot );               call add_char ( '.' )
    case ( t_colon );             call add_char ( ':' )
    case ( t_less_colon );        call add_char ( '<:' )
    case ( t_colon_less );        call add_char ( ':<' )
    case ( t_less_colon_less );   call add_char ( '<:<' )
    case ( t_equal );             call add_char ( '=' )
    case ( t_equal_equal );       call add_char ( '==' )
    case ( t_not_equal );         call add_char ( '/=' )
    case ( t_less );              call add_char ( '<' )
    case ( t_less_eq );           call add_char ( '<=' )
    case ( t_greater );           call add_char ( '>' )
    case ( t_greater_eq );        call add_char ( '>=' )
    case ( t_comma );             call add_char ( ',' )
    case ( t_hat );               call add_char ( '^' )
    case ( t_begin );             call add_char ( 'BEGIN' )
    case ( t_cycle );             call add_char ( 'CYCLE' )
    case ( t_do );                call add_char ( 'DO' )
    case ( t_end );               call add_char ( 'END' )
    case ( t_exit );              call add_char ( 'EXIT' )
    case ( t_and );               call add_char ( 'AND' )
    case ( t_or );                call add_char ( 'OR' )
    case ( t_not );               call add_char ( 'NOT' )
    case ( t_case );              call add_char ( 'CASE' )
    case ( t_default );           call add_char ( 'DEFAULT' )
    case ( t_else );              call add_char ( 'ELSE' )
    case ( t_if );                call add_char ( 'IF' )
    case ( t_select );            call add_char ( 'SELECT' )
    case ( t_then );              call add_char ( 'THEN' )
    case ( t_while );             call add_char ( 'WHILE' )
    case ( t_end_of_input );      call add_char ( '<eof>' )
    case ( t_end_of_stmt );       call add_char ( '<eos>' )
    case ( t_identifier );        call add_char ( '<identifier>' )
    case ( t_number );            call add_char ( '<number>' )
    case ( t_string );            call add_char ( '<string>' )
    case ( t_include );           call add_char ( '#include' )
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
    case ( inc_num );   call add_char ( '<incomplete number>' )
    case ( inc_str );   call add_char ( '<incomplete string>' )
    case ( unk_ch );    call add_char ( '<unrecognized character>' )
    case ( aft_cont );  call add_char ( '<after continuation>' )
    case ( treeNode);   call add_char ( '<tree node>' )
    case ( object );    call add_char ( '<object>' )
    case default
      write ( *, * ) 'SYMBOL_TYPES%INIT_TYPE-E- No initializer for &
                     &terminal symbol type with index ', type
      stop
    end select
  end subroutine INIT_TYPE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SYMBOL_TYPES

! $Log$
! Revision 2.21  2019/06/12 00:14:19  vsnyder
! Mostly repairing comments
!
! Revision 2.20  2014/02/21 19:22:15  vsnyder
! Add CYCLE, DO, EXIT, WHILE
!
! Revision 2.19  2013/12/12 01:58:00  vsnyder
! Add variable definition, and IF and SELECT constructs
!
! Revision 2.18  2013/11/26 22:42:32  vsnyder
! Add CASE, DEFAULT, ELSE, IF, SELECT, THEN reserved words
!
! Revision 2.17  2013/10/02 01:33:10  vsnyder
! Add :=, ? and ! symbols for variable assignment  and conditional expressions
!
! Revision 2.16  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.15  2012/05/05 00:11:51  vsnyder
! Add support for 'not' operator
!
! Revision 2.14  2012/05/01 22:12:43  vsnyder
! Add comment about token names being used in parser
!
! Revision 2.13  2012/05/01 22:10:26  vsnyder
! Add TrueList subroutine
!
! Revision 2.12  2011/04/19 01:59:43  vsnyder
! Support == and /= relational operators too
!
! Revision 2.11  2011/04/18 19:33:26  vsnyder
! Add support for relational operators and boolean-valued expressions
!
! Revision 2.10  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/05/28 23:12:21  vsnyder
! Add power (^) operator
!
! Revision 2.7  2004/04/26 21:55:29  vsnyder
! Make strings case sensitive
!
! Revision 2.6  2004/01/14 18:32:58  vsnyder
! Stuff for Algebra module
!
! Revision 2.5  2002/10/08 00:09:14  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2000/11/30 20:18:47  vsnyder
! Added <: :< and <:< operators.
!
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
