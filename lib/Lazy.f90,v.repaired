! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Lazy

  ! Lazy evaluation functions

  implicit NONE

  private

  public :: AndThen

  ! AndThen ( A, B ) with A optional is false if A is absent,
  ! or A.and.B otherwise.

! This can't be done, because the AndThen function has optional arguments
! interface operator ( .ANDTHEN. )
!   module procedure AndThen
! end interface

  public :: Lazy_Value

  ! Lazy_Value ( A, B ) with A optional.
  ! If A is present, the function value is A; otherwise it is B.

  ! All but Lazy_Value_Character are elemental; for them, the ranks of
  ! A and B must conform.  Lazy_Value_Character is only available for
  ! scalar arguments.

  interface Lazy_Value
    module procedure Lazy_Value_Character, Lazy_Value_Double, Lazy_Value_Int
    module procedure Lazy_Value_Logical, Lazy_Value_Real
  end interface

  public :: Lazy_Len

  ! Lazy_Len ( A, B ) with A optional.
  ! If A is present the function value is len(A); otherwise it is len(B).

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  pure logical function AndThen ( A, B )
    logical, intent(in), optional :: A
    logical, intent(in) :: B
    andThen = present(a)
    if ( andThen ) andThen = a .and. b
  end function AndThen

  pure integer function Lazy_Len ( Optional, Absent )
    character(*), optional, intent(in) :: Optional
    character(*), intent(in) :: Absent
    if ( present(optional) ) then
      lazy_len = len(optional)
    else
      lazy_len = len(absent)
    end if
  end function Lazy_Len

  function Lazy_Value_Character ( Optional, Absent ) result ( R )
    character(*), optional, intent(in) :: Optional
    character(*), intent(in) :: Absent
    character(:), allocatable :: R
    ! Intel 12.0 doesn allow this:
!   allocate ( character(lazy_len(optional,absent)) :: r )
    if ( present(optional) ) then
      allocate ( character(len(optional)) :: r )
      r = optional
    else
      allocate ( character(len(absent)) :: r )
      r = absent
    end if
  end function Lazy_Value_Character

  elemental double precision function Lazy_Value_Double ( Optional, Absent ) result ( R )
    double precision, optional, intent(in) :: Optional
    double precision, intent(in) :: Absent
    if ( present(optional) ) then
      r = optional
    else
      r = absent
    end if
  end function Lazy_Value_Double

  elemental integer function Lazy_Value_Int ( Optional, Absent ) result ( R )
    integer, optional, intent(in) :: Optional
    integer, intent(in) :: Absent
    if ( present(optional) ) then
      r = optional
    else
      r = absent
    end if
  end function Lazy_Value_Int

  elemental logical function Lazy_Value_Logical ( Optional, Absent ) result ( R )
    logical, optional, intent(in) :: Optional
    logical, intent(in) :: Absent
    if ( present(optional) ) then
      r = optional
    else
      r = absent
    end if
  end function Lazy_Value_Logical

  elemental real function Lazy_Value_Real ( Optional, Absent ) result ( R )
    real, optional, intent(in) :: Optional
    real, intent(in) :: Absent
    if ( present(optional) ) then
      r = optional
    else
      r = absent
    end if
  end function Lazy_Value_Real

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Lazy

! $Log$
! Revision 2.4  2012/12/04 03:35:20  vsnyder
! Revert to 2.1 due to restrictions on specification expressions
!
! Revision 2.3  2012/11/30 02:34:40  vsnyder
! Make character function elemental
!
! Revision 2.2  2012/11/30 02:16:43  vsnyder
! Use specification function instead of explicit allocation with a type-spec
! to compute length of result variable for Lazy_Value_Character.
!
! Revision 2.1  2012/11/30 02:07:47  vsnyder
! Initial commit
!
