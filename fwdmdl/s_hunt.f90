! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module S_HUNT_M
  implicit NONE
  private
  public :: S_HUNT, HUNT
  interface HUNT; module procedure S_HUNT; end interface
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.
  pure subroutine S_HUNT ( ELEMENT, ARRAY, N, JLO, JHI )
    use MLSCommon, only: I4, R4
    integer, parameter :: RK = r4
    include 'hunt.f9h'
  end subroutine S_HUNT
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module S_HUNT_M

! $Log$
! Revision 2.2  2005/03/29 01:58:17  vsnyder
! Make stuff pure
!
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
