! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
module MoreMessage ! Messaging with extra functionality
!==============================================================================

  implicit none
  private

  public :: MLSMessage
  public :: MessageWithDouble, MessageWithDoubleArray
  public :: MessageWithInteger, MessageWithIntegerArray
  public :: MessageWithSingle, MessageWithSingleArray
  public :: MessageWithWhere

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module provides nearly-low-level messaging for the MLSPGS suite.
  ! The extra functionality is that it allows to insert numbers or
  ! strings from the string table into the message.

  interface MLSMessage
    module procedure MessageWithDouble, MessageWithDoubleArray
    module procedure MessageWithInteger, MessageWithIntegerArray
    module procedure MessageWithSingle, MessageWithSingleArray
    module procedure MessageWithWhere
  end interface

contains

  ! ------------------------------------------  MessageWithDouble  -----
  subroutine MessageWithDouble ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R" and the
    ! position in files indicated by Where, if present, in place of "%w"
    ! or "%W". then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageScalar.f9h'
  end subroutine MessageWithDouble

  ! -------------------------------------  MessageWithDoubleArray  -----
  subroutine MessageWithDoubleArray ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R" and the
    ! position in files indicated by Where, if present, in place of "%w"
    ! or "%W". then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum(:)
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageArray.f9h'
  end subroutine MessageWithDoubleArray

  ! -----------------------------------------  MessageWithInteger  -----
  subroutine MessageWithInteger ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%d" or "%D". Insert the
    ! string indexed by "Datum" into "Message" in place of "%s" or "%S". 
    ! Insert the line and column number represented by "Datum" in place of
    ! "%l" or "%L".  Insert the signal indexed by "Datum" in place of "%g" or
    ! "%G".  Insert the position in files indicated by Where, if present, in
    ! place of "%w" or "%W".  Then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    use MLSSignals_m, only: GetSignalName
    use String_Table, only: Get_String, String_Length
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in) :: Datum
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    character(512) :: Line ! Should be long enough
    integer :: I, L        ! Next positions in input, Line
    integer :: IERR        ! Status from Get_String
    integer :: N           ! Len_Trim(Message)
    logical :: OK          ! OK to insert

    n = len_trim(message)
    ok = .true.
    i = 1
    l = 1
    do while ( i <= n )
      line(l:l) = message(i:i)
      if ( i >= n ) exit
      if ( (message(i:i+1) == '%d' .or. message(i:i+1) == '%D') .and. ok ) then
        i = i + 2
        write ( line(l:), * ) datum
        line(l:) = adjustl(line(l:))
        l = len_trim(line) + 1
        ok = .false.
      else if ( (message(i:i+1) == '%g' .or. message(i:i+1) == '%G') .and. ok ) then
        i = i + 2
        call getSignalName ( datum, line(l:) )
        l = len_trim(line) + 1
        ok = .false.
      else if ( (message(i:i+1) == '%l' .or. message(i:i+1) == '%L') .and. ok ) then
        i = i + 2
        write ( line(l:), '("line ", i0, ", column ", i0)' ) &
          & datum/256, mod(datum,256)
        l = len_trim(line) + 1
        ok = .false.
      else if ( (message(i:i+1) == '%s' .or. message(i:i+1) == '%S') .and. ok ) then
        i = i + 2
        call get_string ( datum, line(l:), strip=.false., noerror=.true., ierr=ierr )
        if ( ierr == 0 ) then
          l = l + string_length ( datum )
        else
          l = len_trim(line) + 1
        end if
        ok = .false.
      else if ( (message(i:i+1) == '%w' .or. message(i:i+1) == '%W') .and. &
        & present(where) ) then
        i = i + 2
        call get_where ( where, line(l:) )
        l = len_trim(line) + 1
      else
        i = i + 1
        l = l + 1
      end if
    end do
    call MLSMessage ( severity, moduleNameIn, line(:l), advance )
  end subroutine MessageWithInteger

  ! ------------------------------------  MessageWithIntegerArray  -----
  subroutine MessageWithIntegerArray ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%d" or "%D". Insert the
    ! string indexed by "Datum" into "Message" in place of "%s" or "%S". 
    ! Insert the line and column number represented by "Datum" in place of
    ! "%l" or "%L".  Insert the signal indexed by "Datum" in place of "%g" or
    ! "%G".  Insert the position in files indicated by Where, if present, in
    ! place of "%w" or "%W".  Then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    use MLSSignals_m, only: GetSignalName
    use String_Table, only: Get_String, String_Length
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in) :: Datum(:)
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    character(512) :: Line ! Should be long enough
    integer :: I, L        ! Next positions in input, Line
    integer :: IERR        ! Status from Get_String
    integer :: N           ! Len_Trim(Message)
    integer :: ND          ! Next element of Datum

    n = len_trim(message)
    nd = 1
    i = 1
    l = 1
    do while ( i <= n )
      line(l:l) = message(i:i)
      if ( i >= n ) exit
      if ( (message(i:i+1) == '%d' .or. message(i:i+1) == '%D') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        write ( line(l:), * ) datum(nd)
        line(l:) = adjustl(line(l:))
        l = len_trim(line) + 1
        nd = nd + 1
      else if ( (message(i:i+1) == '%g' .or. message(i:i+1) == '%G') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        call getSignalName ( datum(nd), line(l:) )
        l = len_trim(line) + 1
        nd = nd + 1
      else if ( (message(i:i+1) == '%l' .or. message(i:i+1) == '%L') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        write ( line(l:), '("line ", i0, ", column ", i0)' ) &
          & datum(nd)/256, mod(datum(nd),256)
        l = len_trim(line) + 1
        nd = nd + 1
      else if ( (message(i:i+1) == '%s' .or. message(i:i+1) == '%S') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        call get_string ( datum(nd), line(l:), strip=.false., noerror=.true., ierr=ierr )
        if ( ierr == 0 ) then
          l = l + string_length ( datum(nd) )
        else
          l = len_trim(line) + 1
        end if
        nd = nd + 1
      else if ( (message(i:i+1) == '%w' .or. message(i:i+1) == '%W') .and. &
        & present(where) ) then
        i = i + 2
        call get_where ( where, line(l:) )
        l = len_trim(line) + 1
      else
        i = i + 1
        l = l + 1
      end if
    end do
    call MLSMessage ( severity, moduleNameIn, line(:l), advance )
  end subroutine MessageWithIntegerArray

  ! ------------------------------------------  MessageWithSingle  -----
  subroutine MessageWithSingle ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R" and the
    ! position in files indicated by Where, if present, in place of "%w"
    ! or "%W". then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageScalar.f9h'
  end subroutine MessageWithSingle

  ! -------------------------------------  MessageWithSingleArray  -----
  subroutine MessageWithSingleArray ( Severity, ModuleNameIn, Message, Datum, &
    & Where, Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R" and the
    ! position in files indicated by Where, if present, in place of "%w"
    ! or "%W". then call MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum(:)
    type(where_t), intent(in), optional :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageArray.f9h'
  end subroutine MessageWithSingleArray

  ! -------------------------------------------  MessageWithWhere  -----
  subroutine MessageWithWhere ( Severity, ModuleNameIn, Message, Where, &
    & Advance )
    ! Insert "Where" into "Message" in place of "%w" or "%W". then call
    ! MLSMessage.
    use Lexer_Core, only: Get_Where, Where_t
    use MLSMessageModule, only: MLSMessage
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    type(where_t), intent(in) :: Where
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'

    character(512) :: Line ! Should be long enough
    integer :: I, L        ! Next positions in input, Line
    integer :: N           ! Len_Trim(Message)

    n = len_trim(message)
    i = 1
    l = 1
    do
      line(l:l) = message(i:i)
      if ( i >= n ) exit
      if ( (message(i:i+1) == '%w' .or. message(i:i+1) == '%W') ) then
        i = i + 2
        call get_where ( where, line(l:) )
        l = len_trim(line) + 1
      else
        i = i + 1
        l = l + 1
      end if
    end do
    call MLSMessage ( severity, moduleNameIn, line(:l), advance )

  end subroutine MessageWithWhere

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MoreMessage

! $Log$
! Revision 2.10  2015/05/29 00:33:56  vsnyder
! Don't look past the end of the message
!
! Revision 2.9  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.8  2013/06/12 02:13:10  vsnyder
! Cruft removal
!
! Revision 2.7  2011/08/20 02:32:10  vsnyder
! Don't go off the end of Message string
!
! Revision 2.6  2011/08/10 01:47:17  vsnyder
! Correct some comments
!
! Revision 2.5  2011/07/22 18:29:51  vsnyder
! Make more robust to string errors
!
! Revision 2.4  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2006/04/20 01:08:32  vsnyder
! Don't look past the end of the input
!
! Revision 2.2  2005/06/03 01:52:30  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! %g switch to dump signals.
!
! Revision 2.1  2005/05/02 22:55:13  vsnyder
! Initial commit
!
