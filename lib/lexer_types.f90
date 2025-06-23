! Copyright 2018, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module lexer_types

!  Type definitions for the lexer

!  The terminal symbols of the grammar.  The corresponding constants
!  must be typed into the grammar as integer literal values.

  implicit none
  public
  type :: Where_T
    integer :: File = 0            ! String index
    integer :: Source = 0          ! 256*srcLine + srcCol
  end type Where_T

  type :: Token
    integer :: Class               ! token class = token index
    integer :: String_index        ! position in string table
    logical :: Pseudo              ! "token is pseudo-terminal"
    type(where_t) :: Where
  end type Token

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module lexer_types

! $Log$
! Revision 2.1  2018/08/17 21:41:19  pwagner
! First commit
!
