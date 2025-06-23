! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module PhaseFunction

! -------------------------------------------------------------------------  
! SET-UP PHASE FUNCTIONS
!  (but what are phase functions? what are they used for?)
! -------------------------------------------------------------------------

   use MLSCommon, only: r8
   USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
   implicit none
   Private
   Public :: pfsetup

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
      
contains

  subroutine pfsetup(n, p, dp, u, nu)

   ! Arguments

   integer, intent(in)   :: n     ! Maximum number of a, b terms for Mie Cal.  
   integer, intent(in)   :: nu    ! Number of scattering angles
   real(r8), intent(in)  :: u(nu)    ! Cosine of scattering angle
   real(r8), intent(out) :: p(n,nu)  ! Legendre polynomials
   real(r8), intent(out) :: dp(n,nu) ! delta p
   ! Local variables
   integer :: i, j
   real(r8) :: w1, w2, v1, v2, us

   ! Executable
   if ( n < 3) &
     CALL MLSMessage(MLSMSG_Error, ModuleName, &
             & ' n too small to set up phase functions ')
   outer_loop: do i=1, nu
          us = sqrt(1-u(i)*u(i))
          w1 = us
          w2 = 3*u(i)*w1
          p(1,i)=w1
          p(2,i)=w2
          v1 = u(i)
          v2 = 6*u(i)*u(i) - 3
          dp(1,i)=v1
          dp(2,i)=v2
     do j=3, n
            p(j,i)=(2.*j-1.)/(j-1.)*w2*u(i) - j/(j-1.)*w1
            dp(j,i)=(2*j-1.)/(j-1.)*(u(i)*v2-w2*us)-j/(j-1.)*v1
            w1 = w2
            w2 = p(j,i)
            v1 = v2
            v2 = dp(j,i)
     enddo
   end do outer_loop

  END SUBROUTINE PFSETUP

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PhaseFunction

! $Log$
! Revision 2.5  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2007/06/21 00:52:54  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.3  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2003/10/09 16:11:26  jonathan
! add description to input parameter
!
! Revision 2.1  2003/01/31 18:35:55  jonathan
! moved from cldfwm
!
! Revision 1.4  2003/01/30 22:01:14  pwagner
! Cosmetic changes
!
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
