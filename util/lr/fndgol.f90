! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Find_Goal

  implicit NONE
  private
  public :: FNDGOL

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine FNDGOL

    use Error_Handler, only: Error
    use S1, only: MOVSTR
    use S3, only: LFTUSE, PRDIND, PRODCN, RGTUSE
    use TABCOM, only: NVOC
    implicit NONE

  ! Find the goal symbol of the grammar.

  ! *****     External References     ********************************

  ! ERROR   prints error messages.
  ! MOVSTR  moves a string from the symbol table.

  ! *****     Local Variables     ************************************

  ! GOAL    is the goal symbol.
  ! I       is a loop induction variable and subscript.
  ! J       is used during message construction.
  ! LINE    is used for message assembly.

    integer GOAL, I, J
    character(len=120)LINE

  !     *****     Procedures     *****************************************

    goal = 0
    do i = 3, nvoc
      if (lftuse(i) == 0) cycle
      if (rgtuse(i) == 1) cycle
      if (goal /= 0) then
        line = 'Extra goal symbol: '
        j = 20
        call movstr (i, line, j, len(line))
        call error (line(1:j-1),1)
      else
        goal = i
      end if
    end do
    if (goal .eq. 0) goal = prodcn(prdind(2))
    prodcn(3) = goal

  end subroutine FNDGOL

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Find_Goal

! $Log$
