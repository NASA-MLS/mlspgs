! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module WaterVapor

! -------------------------------------------------------------------------  
! COMPUTE VAPOR PRESSURE
! -------------------------------------------------------------------------
          use MLSCommon, only: r8      
          IMPLICIT NONE
          Private
          Public :: RHtoEV

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains
 
          SUBROUTINE RHtoEV(t,RH,EV) 
!         p ---- mb
!         t,td ---- K
!         RH --- %
!         SD --- g/m3                 
          REAL(r8) :: es, t, RH, EV, isat, t0
          !... relative to water
          !         ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))
          !... relative to ice

	t0 = 273.16_r8
	isat = -9.09718*(t0/t-1.0) + 0.78583503 - 3.56654*ALOG10(t0/t) &
          &     		+ 0.876793 * (1.0-t/t0)

        es=10**isat

!	  es=10**(-2667./t+10.555)

!        print*, es
!        print*, ' '
!        print, 10**(-2667./t+10.555)

          if (ES.lt.0._r8) ES=0._r8
          EV=ES*0.01*RH
          
          END SUBROUTINE RHtoEV
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module WaterVapor

! $Log$
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:38  jonathan
! modified F95 version
!
