
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
 
module Tmp

! ---------------------------------
! THIS IS MAY NOT NEED LATER
! ---------------------------------

      use Interpack, only: LOCATE
      use MLSCommon, only: r8
      IMPLICIT NONE
      private
      public :: GET_TAN_PRESS

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains 

      SUBROUTINE GET_TAN_PRESS ( YP, YZ, YT, YQ, NH,            &
                 &               ZPT1, ZZT1, ZTT1, ZVT1, Multi )
!===========================================================================


      INTEGER :: NH                            ! MODEL ATMOSPHERIC LEVELS
      integer :: MULTI, JM, I, J
      REAL(r8) :: YZ(NH)                       ! PRESSURE HEIGHT (m)
      REAL(r8) :: YP(NH)                       ! PRESSURE (hPa)
      REAL(r8) :: YT(NH)                       ! TEMPERATURE PROFILE
      REAL(r8) :: YQ(NH)

      REAL(r8) :: ZPT1(Multi)                  ! TANGENT PRESSURE
      REAL(r8) :: ZZT1(Multi)                  ! TANGENT HEIGHT
      REAL(r8) :: ZTT1(Multi)                  ! TANGENT POINT TEMPERATURE
      REAL(r8) :: ZVT1(Multi) 

      REAL(r8) :: ZH(NH),ZZ(NH)

!---------------------------------------------------------------------------
      
      DO I=1, NH
!         ZH(I)=-LOG10( YP(I) )
          ZH(I)=YZ(I)
      END DO

      DO I=1, Multi
!         ZZ(I)=-LOG10( ZPT1(I) )
          ZZ(I)=ZZT1(I)
      END DO

!----------------------------------------------------------------------

      DO J=1,Multi      

         CALL LOCATE (ZH,NH,NH,ZZ(J),JM)

         ZPT1(J)=((ZH(JM+1)-ZZ(J))*(-LOG10(YP(JM)))+(ZZ(J)-ZH(JM))*  &
     &                (-LOG10(YP(JM+1))))/(ZH(JM+1)-ZH(JM))             
         ZPT1(J) = 10**(-ZPT1(J))

         ZTT1(J)=((ZH(JM+1)-ZZ(J))*YT(JM)+(ZZ(J)-ZH(JM))*    &
     &                YT(JM+1))/(ZH(JM+1)-ZH(JM))             

         ZVT1(J)=((ZH(JM+1)-ZZ(J))*(-LOG10(YQ(JM)))+(ZZ(J)-ZH(JM))*    &
     &                (-LOG10(YQ(JM+1))))/(ZH(JM+1)-ZH(JM))             
         ZVT1(J) = 10**(-ZVT1(J))

      ENDDO

      END SUBROUTINE GET_TAN_PRESS

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Tmp

! $Log$
! Revision 1.5  2005/06/22 18:27:38  pwagner
! Cant have access declared outside module scope
!
! Revision 1.4  2002/12/18 16:10:44  jonathan
! minor changes
!
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:38  jonathan
! modified F95 version
!
