
      SUBROUTINE CLOUD_CHECK(PRESSURE,IWC,LWC,NZ,WC,CHK_CLD,NH)

C=========================================================================C
C  DEFINE VERTICAL PROFILES OF CLOUD ICE-WATER-CONTENT                    C
C  - J.JIANG (05/18/2001)                                                 C
C=========================================================================C

      INTEGER NH,NZ

      REAL WC(2,NH)     ! CLOUD WATER CONTENT (g/m3)
      REAL CHK_CLD(NH)  ! CLOUD CHECKER      

      REAL IWC(NZ)      ! CLOUD ICE WATER CONTENT
      REAL LWC(NZ)      ! CLOUD LIQUID WATER CONTENT
      REAL PRESSURE(NZ) ! 

      INTEGER NH0

      PARAMETER(NH0=2000)
      
      REAL PTOP,PBOTTOM,DP,ZH(NH0),ZA(NH0),WK
      INTEGER I,JM,J

C--------------------------------------------------------------------

      IF (NZ .NE. NH) THEN

         PTOP = 40./16.-3.             ! TOP OF THE MODEL
         PBOTTOM=-ALOG10(PRESSURE(1))  ! BOTTOM OF THE MODEL
         DP=(PTOP-PBOTTOM)/NH          ! LAYER THICKNESS

         DO I=1,NH
            ZH(I)=PBOTTOM+(I-1)*DP
         END DO

         DO I=1,NZ
            ZA(I)=-ALOG10(PRESSURE(I))
         END DO


         DO J=1,NH

            CALL LOCATE (ZA,NZ,NH,ZH(J),JM)             

            WC(1,J)=((ZA(JM+1)-ZH(J))*IWC(JM)+(ZH(J)-ZA(JM))*
     >                IWC(JM+1))/(ZA(JM+1)-ZA(JM))

            WC(2,J)=((ZA(JM+1)-ZH(J))*LWC(JM)+(ZH(J)-ZA(JM))*
     >                LWC(JM+1))/(ZA(JM+1)-ZA(JM))
         ENDDO


      ELSE    
         DO J=1,NH
            WC(1,J)=IWC(J)
            WC(2,J)=LWC(J)
         ENDDO
      ENDIF

      DO I=1,NH
         CHK_CLD(I)=WC(1,I)+WC(2,I)
      ENDDO
C---------------------------------------------------------------------------

      RETURN
      END



