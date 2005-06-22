! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DCSPLINE_DER_M

  implicit NONE
  private
  public :: DCSPLINE_DER, CSPLINE_DER

  interface CSPLINE_DER; module procedure DCSPLINE_DER; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
  subroutine DCSPLINE_DER ( XIN, XOUT, YIN, YOUT, DYOUT, NIN, NOUT, YMIN, YMAX)
    use D_HUNT_M, only: HUNT
    use D_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0d0)
    include 'cspline_der.f9h'
  end subroutine DCSPLINE_DER
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DCSPLINE_DER_M

! $Log$
! Revision 2.3  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/10/04 17:43:58  vsnyder
! Move USE statements to procedure scope
!
! Revision 2.1  2002/04/18 10:46:25  zvi
! Adding optional limits
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.5  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/03/09 00:40:32  zvi
! Correcting an error in HUNT routine
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
