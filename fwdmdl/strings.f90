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
end module STRINGS

! $Log$
