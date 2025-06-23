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

!  The terminal symbols of the grammar.  The corresponding constants
!  must be typed into the grammar as integer literal values.

  use STRING_TABLE, only: ADD_CHAR

  implicit NONE
  public

  ! Terminal symbol class numbers.  These depend on the language. If this
  ! list changes, you need to change the array TERM_TYPES below, and might
  ! need to change the table GEN, indexed by them, in parser.
  integer, parameter :: T_NULL = 0                ! MUST be present and first
  integer, parameter :: T_EQUAL = 1               ! =
  integer, parameter :: T_PRODUCES = 2            ! ->
  integer, parameter :: T_GENERATES = 3           ! =>
  integer, parameter :: T_QUESTION = 4            ! ?
  integer, parameter :: T_EOS = 5                 ! End of statement
  integer, parameter :: T_END_OF_INPUT = 6        ! EOF
  ! T_IDENTIFIER, T_NUMBER, T_STRING must be consecutive
  integer, parameter :: T_IDENTIFIER = T_END_OF_INPUT + 1 ! <IDENTIFIER>
  integer, parameter :: T_NUMBER =     T_IDENTIFIER + 1   ! <NUMBER>
  integer, parameter :: T_STRING =     T_NUMBER + 1       ! <STRING>
  integer, parameter :: T_UNK_OP =     T_STRING + 1       ! unknown operator
  integer, parameter :: T_UNK_PUN =    T_UNK_OP + 1       ! unknown punctuator
  integer, parameter :: T_UNK_CH =     T_UNK_PUN + 1      ! unknown character
  integer, parameter :: T_INC_STR =    T_UNK_CH + 1       ! incomplete string
  integer, parameter :: T_AFT_CONT =   T_INC_STR + 1      ! junk after continuation
  integer, parameter :: T_TREE_NODE =  T_AFT_CONT + 1     ! For symbol table dumps
  integer, parameter :: T_INCLUDE =    T_TREE_NODE + 1    ! Not used, but the
                                                          ! parser wants it

! The parameters T_LAST_TERMINAL, MIN_PSEUDO, MAX_PSEUDO and CASELESS_LOOK
! MUST be defined.
  integer, parameter :: T_LAST_TERMINAL = T_INCLUDE

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

  integer, parameter :: IDENT =    1              ! Identifier
  integer, parameter :: NUMCON =   IDENT + 1      ! Numeric constant
  integer, parameter :: STRING =   NUMCON + 1     ! String in quotes
  integer, parameter :: DEF_OP =   STRING + 1     ! Defined operator
  integer, parameter :: UNK_OP =   DEF_OP + 1     ! Unefined operator
  integer, parameter :: DEF_PUN =  UNK_OP + 1     ! Defined punctuator
  integer, parameter :: UNK_PUN =  DEF_PUN + 1    ! Undefined punctuator
  integer, parameter :: INC_NUM =  UNK_PUN + 1    ! Incomplete number
  integer, parameter :: INC_STR =  INC_NUM + 1    ! Incomplete string
  integer, parameter :: UNK_CH =   INC_STR + 1    ! Unknown character
  integer, parameter :: AFT_CONT = UNK_CH + 1     ! After continuation
  integer, parameter :: TREENODE = AFT_CONT + 1   ! Tree node
  integer, parameter :: OBJECT =   TREENODE + 1   ! None of the above

  integer, parameter :: FIRST_TYPE = IDENT
  integer, parameter :: LAST_TYPE = OBJECT

! The array TERM_TYPES gives the terminal type of each class of terminal.
! It must be defined.
  integer, parameter :: TERM_TYPES(t_null: t_last_terminal) = &
  !  t_null    =         ->        =>        ?         <eos>   <eof>
  (/ object,   def_op,   def_op,   def_op,   def_pun,  object, object,  &
  !  <ident>   <numcon>  <string>  unk_op    unk_pun   unk_ch  inc_str
     ident,    numcon,   string,   unk_op,   unk_pun,  unk_ch, inc_str, &
  !  junk      <tree>    <include>
     aft_cont, treenode, object /)

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
  ! pseudo-terminals (e.g. identifier, integer).

    integer, intent(in) :: TERMINAL

    select case ( terminal )
    case ( t_null );          call add_char ( '' )
    case ( t_equal );         call add_char ( '=' )
    case ( t_produces );      call add_char ( '->' )
    case ( t_generates );     call add_char ( '=>' )
    case ( t_question );      call add_char ( '?' )
    case ( t_eos );           call add_char ( '<eos>' )
    case ( t_end_of_input );  call add_char ( '<eof>' )
    case ( t_identifier );    call add_char ( '<identifier>' )
    case ( t_number );        call add_char ( '<number>' )
    case ( t_string );        call add_char ( '<string>' )
    case ( t_unk_op );        call add_char ( '<unk_op>' )
    case ( t_unk_pun );       call add_char ( '<unk_pun>' )
    case ( t_unk_ch );        call add_char ( '<unk_ch>' )
    case ( t_inc_str );       call add_char ( '<inc_str>' )
    case ( t_aft_cont );      call add_char ( '<aft_cont>' )
    case ( t_tree_node );     call add_char ( '<tree node>' )
    case ( t_include );       call add_char ( '<include>' )
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
    case ( ident );     call add_char ( '<identifier>' )
    case ( numcon );    call add_char ( '<number>' )
    case ( string );    call add_char ( '<string>' )
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
                     &type with index ', type
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
! Revision 1.2  2019/07/09 20:28:37  vsnyder
! Store string indices of intrinsic symbols
!
! Revision 1.1  2014/01/14 01:36:18  vsnyder
! Renamed with _LR suffix
!
! Revision 1.1  2014/01/14 00:15:08  vsnyder
! Initial commit of new module for new LR
!
