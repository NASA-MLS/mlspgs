! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module STRINGS
  implicit NONE
  public
  integer, private, parameter :: UPDN = ichar('a') - ichar('A')
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module STRINGS
! $Log$
! Revision 2.2  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:28  livesey
! New forward model
!
! Revision 1.6  2001/06/07 23:39:32  pwagner
! Added Copyright statement
!
! Revision 1.5  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.2  2000/05/04 23:43:32  vsnyder
! Initial Fortran 90 version
!
