! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module STRINGS
  implicit NONE
  public
  integer, private, parameter :: UPDN = ichar('a') - ichar('A')
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  subroutine LEFTJ ( S )
  ! Left justify S.
    character(len=*), intent(inout) :: S
    integer :: I
    do i = 1, len(s)
      if ( s(i:i) /= ' ' ) then
        if ( i /= 1 ) s(1:) = s(i:)
        return
      end if
    end do
  end subroutine LEFTJ
!---------------------------------------------------------------------------
  subroutine SQZSTR ( S )
  ! Squeeze all the blanks out of S.
    character(len=*), intent(inout) :: S
    integer :: I, J
    j = 0
    do i = 1, len(s)
      if ( s(i:i) /= ' ' ) then
        j = j + 1
        s(j:j) = s(i:i)
      end if
    end do
    s(j+1:) = ' '
  end subroutine SQZSTR
!---------------------------------------------------------------------------
  subroutine STRLWR ( S )
  ! Convert all of the letters in S to lower case.  Leave everything
  ! else as-is.
    character(len=*), intent(inout) :: S
    integer :: I
    do i = 1, len(s)
      if ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) &
        s(i:i) = char(ichar(s(i:i)) + updn)
    end do
  end subroutine STRLWR
!---------------------------------------------------------------------------
  subroutine STRUPR ( S )
  ! Convert all of the letters in S to upper case.  Leave everything
  ! else as-is.
    character(len=*), intent(inout) :: S
    integer :: I
    do i = 1, len(s)
      if ( s(i:i) >= 'a' .and. s(i:i) <= 'z' ) &
        s(i:i) = char(ichar(s(i:i)) - updn)
    end do
  end subroutine STRUPR
!---------------------------------------------------------------------------
end module STRINGS
! $Log$
! Revision 1.6  2001/06/07 23:39:32  pwagner
! Added Copyright statement
!
! Revision 1.5  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.2  2000/05/04 23:43:32  vsnyder
! Initial Fortran 90 version
!
