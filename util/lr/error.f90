! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Error_m

  implicit NONE
  private
  public :: Error

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine ERROR ( MSG, IBAD )

    use Output_m, only: OUTPUT
    implicit NONE

  ! Print an error message, and terminate execution if necessary.

  ! MSG     is the message to be printed, or a part of it.
  ! IBAD    is a code indicating the severity of the message:
  !         If the low order digit of IBAD is zero the text in MSG is
  !         simply printed.  If the low order digit of IBAD is 1 the text
  !         in MSG is printed with the prefix WARNING.  Otherwise the text
  !         in MSG is printed with the prefix FATAL, and execution
  !         is terminated by a STOP statement with stop code of 6.

        character(len=*), intent(in) :: MSG
        integer, intent(in) :: IBAD

  !     *****     Local Variables     ************************************

  ! JBAD    is the low order digit of IBAD.

    integer  JBAD

  ! *****     Procedures     *****************************************

    call output ( ' *** ' )
    jbad = mod(ibad, 10)
    if ( jbad == 1 ) then
      call output ( 'WARNING' )
    else if ( jbad >= 1 ) then
      call output ( 'FATAL' )
    end if
    call output ( trim(msg), advance='yes' )
    if (jbad >= 2) stop 6

  end subroutine ERROR

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Error_m
! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
