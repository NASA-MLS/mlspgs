! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

	subroutine pfsetup(n, p, dp, u, nu)

   ! Arguments
   ! (Jonathan--please describe inputs and outputs; even partrially)
   integer, intent(in)   :: n     
   integer, intent(in)   :: nu    
	real(r8), intent(in)  :: u(nu) 
	real(r8), intent(out) :: p(n,nu)
	real(r8), intent(out) :: dp(n,nu)
   ! Local variables
	integer :: i, j
	real(r8) :: w1, w2, v1, v2, us

   ! Executable
   if ( n < 3) &
     CALL MLSMessage(MLSMSG_Error, ModuleName, &
             & ' n too small to set up phase functions ')
   outer_loop: do i=1, nu
	  ! do 100 i=1,nu
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
	  ! do j=3,n
	    p(j,i)=(2.*j-1.)/(j-1.)*w2*u(i) - j/(j-1.)*w1
	    dp(j,i)=(2*j-1.)/(j-1.)*(u(i)*v2-w2*us)-j/(j-1.)*v1
	    w1 = w2
	    w2 = p(j,i)
	    v1 = v2
	    v2 = dp(j,i)
	  enddo
     ! 100	continue
   end do outer_loop

  END SUBROUTINE PFSETUP

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PhaseFunction

! $Log$
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
