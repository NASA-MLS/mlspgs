! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module S_PCSPL_M

  implicit NONE

  private
  public :: S_PCSPL, PCSPL

  interface PCSPL; module procedure S_PCSPL; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine S_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    integer, parameter :: RK = kind(0.0e0)
    include 'pcspl.f9h'
  end subroutine S_PCSPL

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module S_PCSPL_M

! $Log$
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
