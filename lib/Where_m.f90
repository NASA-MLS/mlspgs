! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Where_M

  implicit NONE
  private
  public :: Where

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
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
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Where_M

! $Log$
! Revision 2.2  2003/05/16 21:59:43  vsnyder
! Fix a subscript error.  Fortunately, nothing uses it yet
!
! Revision 2.1  2003/05/16 01:59:53  vsnyder
! Initial commit
!
