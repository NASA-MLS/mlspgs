! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module S_CSPLINE_M

  implicit NONE

  private
  public S_CSPLINE, CSPLINE

  interface CSPLINE; module procedure S_CSPLINE; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
!
  subroutine S_CSPLINE (XIN, XOUT, YIN, YOUT, NIN, NOUT, YMIN, YMAX)
    use S_HUNT_M, only: HUNT
    use S_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0e0)
    include 'cspline.f9h'
  end subroutine S_CSPLINE

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module S_CSPLINE_M

! $Log$
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
