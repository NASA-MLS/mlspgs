! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module testfield_m

implicit none
private

public :: testfield

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains 

subroutine testfield ( h, ct, stcp, stsp)
    use MLSCommon, only: Rk => Rp
    real(rk), intent(inout) ::    h(:)      ! mag field magnitude
    real(rk), intent(inout) ::    CT(:)     ! Cos(Theta)
    real(rk), intent(inout) ::    STCP(:)   ! Sin(Theta) Cos(Phi)
    real(rk), intent(inout) ::    STSP(:)   ! Sin(Theta) Sin(Phi)

    h=0
    ct=0
    stcp=0.707
    stsp=0.707
end subroutine testfield

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module testfield_m
! $Log$
! Revision 2.2  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.1  2003/03/19 00:17:53  michael
! Procedure to cobble dummy magnetic field path into code without recompiling
! all of Full_Foreward_Model
!
