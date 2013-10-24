! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Error_Handler

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

    use IO, only: OUTPUT
    use READCM, only: LINENO, LINPTR
    implicit NONE

  ! Print an error message, and terminate execution if necessary.

  ! MSG     is the message to be printed, or a part of it.
  ! IBAD    is a code indicating the severity of the message, and whether
  !         line and column number of the input is to be printed:
  !         If the low order digit of IBAD is zero the text in MSG is
  !         simply printed.  If the low order digit of IBAD is 1 the text
  !         in MSG is printed with the prefix WARNING.  Otherwise the text
  !         in MSG is printed with the prefix FATAL, and execution
  !         is terminated by a call to EXIT.  If IBAD is greater than 9,
  !         the line and column number are printed after the text in MSG.

        character(len=*), intent(in) :: MSG
        integer, intent(in) :: IBAD

  !     *****     External References     ********************************

  ! EXIT    terminates execution.

  !     *****     Local Variables     ************************************

  ! IBASE   is the current position in LINE at which data are being
  !         stored.
  ! JBAD    is the low order digit of IBAD.
  ! LINE    is used to assemble the message.

    integer IBASE, JBAD
    character(len=120) :: LINE

  ! *****     Procedures     *****************************************

    line = ' ***'
    jbad = mod(ibad, 10)
    if (jbad == 0) then
      ibase = 5
    else if (jbad == 1) then
      line(6:12) = 'WARNING'
      ibase = 13
    else
      line(6:10) = 'FATAL'
      ibase = 11
    end if
    line(ibase+1:ibase+len(msg)) = msg
    ibase = ibase + len(msg) + 1
    if (ibad.gt.9) then
      line(ibase+1:ibase+7) = 'AT LINE'
      ibase = ibase + 8
      write ( line(ibase:ibase+3), '(i4)' ) lineno
      ibase = ibase + 4
      line(ibase+1:ibase+4) = 'CHAR'
      ibase = ibase + 5
      write ( line(ibase:ibase+3), '(i4)' ) linptr
      ibase = ibase + 3
    end if
    call output (line(1:ibase))
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

end module Error_Handler
! $Log$
