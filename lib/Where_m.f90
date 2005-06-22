! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Where_M

  implicit NONE
  private
  public :: Where

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  subroutine Where ( Log, Int )

  ! Produce a list of indices in Int that are subscripts of all of the
  ! true elements of Log.

    logical, intent(in) :: Log(:)
    integer, intent(out) :: Int(:)

    integer :: I, P_I

    i = 0
    do p_i = 1, size(log)
      if ( log(p_i) ) then
        i = i + 1
        int(i) = p_i
      end if
    end do

  end subroutine Where

!----------------------------------------------------------------------
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Where_M

! $Log$
! Revision 2.3  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2003/05/16 21:59:43  vsnyder
! Fix a subscript error.  Fortunately, nothing uses it yet
!
! Revision 2.1  2003/05/16 01:59:53  vsnyder
! Initial commit
!
