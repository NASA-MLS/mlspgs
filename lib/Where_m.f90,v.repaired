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
  public :: Where, WhereOnly, HereAndThere

  interface Where
    module procedure WhereOnly, HereAndThere
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  subroutine WhereOnly ( Log, Int, nTrue )

  ! Produce a list of indices in Int that are subscripts of all of the
  ! true elements of Log.

    logical, intent(in) :: Log(:)
    integer, intent(out) :: Int(:)
    integer, intent(out), optional :: nTrue ! Size of useful part of Int

    integer :: I, P_I

    i = 0
    do p_i = 1, size(log)
      if ( log(p_i) ) then
        i = i + 1
        int(i) = p_i
      end if
    end do

    if ( present(nTrue) ) nTrue = i

  end subroutine WhereOnly

  subroutine HereAndThere ( Log, Here, nHere, There )

  ! Produce a list of indices in Here that are subscripts of all of the
  ! true elements of Log, the number of true entries in nHere, and a list
  ! of indices in There that are subscripts of all of the false elements
  ! of log (the number of these is size(Log) - nHere).

    logical, intent(in) :: Log(:)
    integer, intent(out) :: Here(:)
    integer, intent(out) :: nHere
    integer, intent(out) :: There(:)

    integer :: I, P_I

    i = 0
    nHere = 0
    do p_i = 1, size(log)
      if ( log(p_i) ) then
        nHere = nHere + 1
        here(nHere) = p_i
      else
        i = i + 1
        there = p_i
      end if
    end do

  end subroutine HereAndThere

!----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Where_M

! $Log$
! Revision 2.6  2016/09/06 20:37:15  vsnyder
! Add NTrue optional argument to WhereOnly
!
! Revision 2.5  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2007/06/06 00:21:51  vsnyder
! Make Where generic for WhereOnly, HereAndThere
!
! Revision 2.3  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2003/05/16 21:59:43  vsnyder
! Fix a subscript error.  Fortunately, nothing uses it yet
!
! Revision 2.1  2003/05/16 01:59:53  vsnyder
! Initial commit
!
