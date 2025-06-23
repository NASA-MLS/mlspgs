! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module WaterVapor

! -------------------------------------------------------------------------  
! COMPUTE VAPOR PRESSURE
! -------------------------------------------------------------------------
          use MLSCommon, only: r8      
          IMPLICIT NONE
          Private
          Public :: RHtoEV

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
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
        isat = -9.09718*(t0/t-1.0) + 0.78583503 - 3.56654*LOG10(t0/t) &
          &                        + 0.876793 * (1.0-t/t0)

        es=10**isat

!        es=10**(-2667./t+10.555)

!        print*, es                       
!        print*, ' '                      
!        print, 10**(-2667./t+10.555)     

          if (ES.lt.0._r8) ES=0._r8
          EV=ES*0.01*RH
          
          END SUBROUTINE RHtoEV
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module WaterVapor

! $Log$
! Revision 1.7  2007/10/04 18:20:33  vsnyder
! Remove tabs, which are not allowed by the Fortran standard
!
! Revision 1.6  2005/06/22 18:27:38  pwagner
! Cant have access declared outside module scope
!
! Revision 1.5  2002/10/30 01:04:10  jonathan
! change alog to log
!
! Revision 1.4  2002/10/30 00:53:36  jonathan
! new water vapor formula
!
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:38  jonathan
! modified F95 version
!
