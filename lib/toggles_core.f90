! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TOGGLES

  implicit NONE
  public

  ! Toggle indices.  The parser uses PAR, and the symbol/declaration
  ! table software uses TAB, but they can have any values.
  integer, parameter :: CON = 1    ! Trace constrainers visitation of tree
                                   ! nodes, inverted by @C
  integer, parameter :: EMIT = 2   ! Trace code emitting, inverted by @E
  integer, parameter :: GEN = 3    ! Trace code generator, inverted by @G
  integer, parameter :: LEX = 4    ! Print each token as lexer finishes,
                                   ! inverted by @L
  integer, parameter :: PAR = 5    ! Trace parser actions, inverted by @P
  integer, parameter :: SYN = 6    ! Print abstract syntax tree, inverted
                                   ! by @A
  integer, parameter :: TAB = 7    ! Trace symbol/declaration table
                                   ! activity, inverted by @S

  ! Toggles:
  logical :: TOGGLE(con:tab) = .false.

  ! In case you want to have levels of output (not used in the parser):
  integer :: LEVELS(con:tab) = 0

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine INIT_TOGGLE
    toggle = .false.
  end subroutine INIT_TOGGLE

end module TOGGLES

! $Log$
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.2  2000/08/01 01:01:45  vsnyder
! Added "levels" for those who want to control the amount of output.
!
