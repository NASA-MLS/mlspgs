! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

!==============================================================================
module MoreMessage ! Messaging with extra functionality
!==============================================================================

  implicit none
  private

  public :: MLSMessage
  public :: MessageWithDouble, MessageWithDoubleArray
  public :: MessageWithInteger, MessageWithIntegerArray
  public :: MessageWithSingle, MessageWithSingleArray

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
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
  end interface

contains

  ! ------------------------------------------  MessageWithDouble  -----
  subroutine MessageWithDouble ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R", then call
    ! MLSMessage.
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageScalar.f9h'
  end subroutine MessageWithDouble

  ! -------------------------------------  MessageWithDoubleArray  -----
  subroutine MessageWithDoubleArray ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R", then call
    ! MLSMessage.
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum(:)
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageArray.f9h'
  end subroutine MessageWithDoubleArray

  ! -----------------------------------------  MessageWithInteger  -----
  subroutine MessageWithInteger ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%d" or "%D".
    ! Insert the string indexed by "Datum" into "Message" in place of
    ! "%s" or "%S".  Then call MLSMessage.
    use MLSMessageModule, only: MLSMessage
    use String_Table, only: Get_String, String_Length
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in) :: Datum
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    character(512) :: Line ! Should be long enough
    integer :: I, L        ! Next positions in input, Line
    integer :: N           ! Len_Trim(Message)
    logical :: OK          ! OK to insert

    n = len_trim(message)
    ok = .true.
    i = 1
    l = 1
    do
      line(l:l) = message(i:i)
      if ( i >= n ) exit
      if ( (message(i:i+1) == '%d' .or. message(i:i+1) == '%D') .and. ok ) then
        i = i + 2
        write ( line(l:), * ) datum
        line(l:) = adjustl(line(l:))
        l = len_trim(line) + 1
        ok = .false.
      else if ( (message(i:i+1) == '%s' .or. message(i:i+1) == '%S') .and. ok ) then
        i = i + 2
        call get_string ( datum, line(l:), strip=.true. )
        l = l + string_length ( datum ) - 2
        ok = .false.
      else if ( (message(i:i+1) == '%l' .or. message(i:i+1) == '%L') .and. ok ) then
        i = i + 2
        write ( line(l:), '("line ", i0, ", column ", i0)' ) &
          & datum/256, mod(datum,256)
        l = len_trim(line) + 1
        ok = .false.
      else
        i = i + 1
        l = l + 1
      end if
    end do
    call MLSMessage ( severity, moduleNameIn, line(:l), advance )
  end subroutine MessageWithInteger

  ! ------------------------------------  MessageWithIntegerArray  -----
  subroutine MessageWithIntegerArray ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%d" or "%D".
    ! Insert the string indexed by "Datum" into "Message" in place of
    ! "%s" or "%S".  Insert the line and cokumn number represented by "Datum"
    ! in place of "%l" or "%L".  Then call MLSMessage.
    use MLSMessageModule, only: MLSMessage
    use String_Table, only: Get_String, String_Length
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    integer, intent(in) :: Datum(:)
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    character(512) :: Line ! Should be long enough
    integer :: I, L        ! Next positions in input, Line
    integer :: N           ! Len_Trim(Message)
    integer :: ND          ! Next element of Datum

    n = len_trim(message)
    nd = 1
    i = 1
    l = 1
    do
      line(l:l) = message(i:i)
      if ( i >= n ) exit
      if ( (message(i:i+1) == '%d' .or. message(i:i+1) == '%D') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        write ( line(l:), * ) datum(nd)
        line(l:) = adjustl(line(l:))
        l = len_trim(line) + 1
        nd = nd + 1
      else if ( (message(i:i+1) == '%s' .or. message(i:i+1) == '%S') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        call get_string ( datum(nd), line(l:), strip=.true. )
        l = l + string_length ( datum(nd) ) - 2
        nd = nd + 1
      else if ( (message(i:i+1) == '%l' .or. message(i:i+1) == '%L') .and. &
        & nd <= size(datum) ) then
        i = i + 2
        write ( line(l:), '("line ", i0, ", column ", i0)' ) &
          & datum(nd)/256, mod(datum(nd),256)
        l = len_trim(line) + 1
        nd = nd + 1
      else
        i = i + 1
        l = l + 1
      end if
    end do
    call MLSMessage ( severity, moduleNameIn, line(:l), advance )
  end subroutine MessageWithIntegerArray

  ! ------------------------------------------  MessageWithSingle  -----
  subroutine MessageWithSingle ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R", then call
    ! MLSMessage.
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageScalar.f9h'
  end subroutine MessageWithSingle

  ! -------------------------------------  MessageWithSingleArray  -----
  subroutine MessageWithSingleArray ( Severity, ModuleNameIn, Message, Datum, &
    & Advance )
    ! Insert "Datum" into "Message" in place of "%r" or "%R", then call
    ! MLSMessage.
    use MLSMessageModule, only: MLSMessage
    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module
    character (len=*), intent(in) :: Message ! Line of text
    real(rk), intent(in) :: Datum(:)
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    include 'MoreMessageArray.f9h'
  end subroutine MessageWithSingleArray

!=======================================================================
  logical function not_used_here()
    !---------------------------- RCS Ident Info -----------------------
    character (len=len(idParm)) :: Id = idParm
    !-------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MoreMessage

! $Log$
! Revision 2.1  2005/05/02 22:55:13  vsnyder
! Initial commit
!
