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
 !---------------------------------------------------------------------------
      
contains
 
          SUBROUTINE RHtoEV(t,RH,EV) 
!         p ---- mb
!         t,td ---- K
!         RH --- %
!         SD --- g/m3                 
          REAL(r8) :: es, t, RH, EV
          !... relative to water
          !         ES=6.10779*EXP(17.28*(t-273.15)/(t-36.))
          !... relative to ice

	  es=10**(-2667./t+10.555)
          if (ES.lt.0._r8) ES=0._r8
          EV=ES*0.01*RH
          
          END SUBROUTINE RHtoEV
end module WaterVapor

! $Log$
