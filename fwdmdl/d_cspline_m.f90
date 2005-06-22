! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module D_CSPLINE_M

  implicit NONE

  private
  public D_CSPLINE, CSPLINE

  interface CSPLINE; module procedure D_CSPLINE; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
!
  subroutine D_CSPLINE (XIN, XOUT, YIN, YOUT, NIN, NOUT, YMIN, YMAX)
    use D_HUNT_M, only: HUNT
    use D_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0d0)
    include 'cspline.f9h'
  end subroutine D_CSPLINE

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module D_CSPLINE_M
! $Log$
! Revision 2.4  2002/10/14 20:51:53  vsnyder
! Repair CVS stuff broken before previous commit
!
! Revision 2.3  2002/10/14 20:51:06  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.2  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/04/18 10:46:24  zvi
! Adding optional limits
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.6  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.5  2001/03/09 00:40:32  zvi
! Correcting an error in HUNT routine
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
