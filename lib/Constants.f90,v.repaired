! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Constants

! Provide several numerical constants.

  use MLSKinds, only: R8

  implicit NONE
  public

  real(r8), parameter :: Ln2 = 0.6931471805599453094172321214581765680755_r8
  real(r8), parameter :: Ln10 = 2.302585092994045684017991454684364207601_r8
  real(r8), parameter :: Pi = 3.141592653589793238462643383279502884197_r8
  real(r8), parameter :: Deg2Rad = Pi/180.0_r8 ! Degrees-to-Radians
  real(r8), parameter :: Rad2Deg = 180.0_r8/Pi ! Radians-to-Degrees
  real(r8), parameter :: Sqrtln2 = 0.8325546111576977563531646448952010476306_r8
  real(r8), parameter :: SqrtPi = 1.772453850905516027298167483341145182798_r8
  real(r8), parameter :: Sqrt2 = 1.414213562373095048801688724209698078569_r8

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Constants

! $Log$
! Revision 2.4  2016/07/30 00:50:57  vsnyder
! Add Ln2
!
! Revision 2.3  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2007/12/19 03:58:48  vsnyder
! Add Sqrt2
!
! Revision 2.1  2007/12/06 20:36:27  vsnyder
! Initial commit
!
