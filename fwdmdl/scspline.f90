! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SCSPLINE_DER_M

  implicit NONE
  private
  public :: SCSPLINE_DER, CSPLINE_DER

  interface CSPLINE_DER; module procedure SCSPLINE_DER; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
  subroutine SCSPLINE_DER ( XIN, XOUT, YIN, YOUT, DYOUT, NIN, NOUT, YMIN, YMAX)
    use S_HUNT_M, only: HUNT
    use S_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0e0)
    include 'cspline_der.f9h'
  end subroutine SCSPLINE_DER
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SCSPLINE_DER_M

! $Log$
! Revision 2.4  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2003/11/10 17:43:56  pwagner
! Correct spelling of not_used_here in other place
!
! Revision 2.2  2003/11/08 02:03:05  vsnyder
! Correct spelling of not_used_here
!
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
