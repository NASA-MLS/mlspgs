! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PhaseFunction

! -------------------------------------------------------------------------  
! SET-UP PHASE FUNCTIONS
! -------------------------------------------------------------------------

        use MLSCommon, only: r8
	implicit none
	Private
        Public :: pfsetup

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

	subroutine pfsetup(n,p,dp,u,nu)

	integer :: i,j,n,nu
	real(r8) :: u(nu)
	real(r8) :: p(n,nu)
	real(r8) :: dp(n,nu)
	real(r8) :: w1, w2, v1, v2, us

	do 100 i=1,nu
	us = sqrt(1-u(i)*u(i))
	w1 = us
	w2 = 3*u(i)*w1
	p(1,i)=w1
	p(2,i)=w2
	v1 = u(i)
	v2 = 6*u(i)*u(i) - 3
	dp(1,i)=v1
	dp(2,i)=v2
	  do j=3,n
	  p(j,i)=(2.*j-1.)/(j-1.)*w2*u(i) - j/(j-1.)*w1
	  dp(j,i)=(2*j-1.)/(j-1.)*(u(i)*v2-w2*us)-j/(j-1.)*v1
	  w1 = w2
	  w2 = p(j,i)
	  v1 = v2
	  v2 = dp(j,i)
	  enddo
 100	continue

      END SUBROUTINE PFSETUP

end module PhaseFunction

! $Log$
