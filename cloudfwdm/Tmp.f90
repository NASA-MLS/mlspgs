
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module Tmp

! ---------------------------------
! THIS IS MAY NOT NEED LATER
! ---------------------------------

      use Interpack, only: LOCATE
      use MLSCommon, only: r8
      IMPLICIT NONE
      private
      public :: GET_TAN_PRESS

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
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

         ZVT1(J)=((ZH(JM+1)-ZZ(J))*YQ(JM)+(ZZ(J)-ZH(JM))*    &
     &                YQ(JM+1))/(ZH(JM+1)-ZH(JM))             

      ENDDO

      END SUBROUTINE GET_TAN_PRESS

end module Tmp

! $Log: Tmp.f90,v      

      
















