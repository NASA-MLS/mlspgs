! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module TIME_M
!=============================================================================

! Compute either CPU time, in arbitrary units, or wall-clock time, in
! seconds since midnight.

  implicit NONE
  private

  public :: TIME_NOW, TIME_NOW_D, TIME_NOW_S, USE_WALL_CLOCK

  interface TIME_NOW
    module procedure TIME_NOW_D
    module procedure TIME_NOW_S
  end interface

  logical, save :: USE_WALL_CLOCK = .false.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine TIME_NOW_D ( T )
    integer, parameter :: RK = kind(0.0d0)
    include "time_now.f9h"
  end subroutine TIME_NOW_D

  subroutine TIME_NOW_S ( T )
    integer, parameter :: RK = kind(0.0e0)
    include "time_now.f9h"
  end subroutine TIME_NOW_S

end module TIME_M

!$Log$
!Revision 2.1  2001/11/09 22:45:30  vsnyder
!Initial commit
!
