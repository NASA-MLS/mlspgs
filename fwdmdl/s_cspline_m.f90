! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module S_CSPLINE_M

  implicit NONE

  private
  public S_CSPLINE, CSPLINE

  interface CSPLINE; module procedure S_CSPLINE; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
!---------------------------------------------------------------------------

contains
!
  subroutine S_CSPLINE (XIN, XOUT, YIN, YOUT, NIN, NOUT, YMIN, YMAX)
    use S_HUNT_M, only: HUNT
    use S_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0e0)
    include 'cspline.f9h'
  end subroutine S_CSPLINE

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module S_CSPLINE_M

! $Log$
