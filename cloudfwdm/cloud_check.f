
      SUBROUTINE CLOUD_CHECK(PRESSURE,IWC,LWC,IPSD,NZ,
     >                        PSDF,CHK_CLD,NH)

C=========================================================================C
C  CHECK VERTICAL PROFILES OF CLOUD ICE-WATER-CONTENT                     C
C  - J.JIANG (05/18/2001)                                                 C
C=========================================================================C

      INTEGER NH,NZ
C---------------------------------------------------
C     MODEL GRID
C---------------------------------------------------
      INTEGER PSDF(NH)  ! SIZE-DISTRIBUTION FLAG
      REAL WC(2,640)     ! CLOUD WATER CONTENT (g/m3)
                        ! (1=ICE,2=WATER)
      REAL CHK_CLD(NH)  ! CLOUD CHECKER      
C--------------------------------------------------
C     L2 GRID
C--------------------------------------------------
      INTEGER IPSD(NZ)  ! SIZE-DISTRIBUTION FLAG
      REAL IWC(NZ)      ! CLOUD ICE WATER CONTENT
      REAL LWC(NZ)      ! CLOUD LIQUID WATER CONTENT
      REAL PRESSURE(NZ) ! L2 ATMOSPHERIC PRESSURE

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

            PSDF(J)=((ZA(JM+1)-ZH(J))*IPSD(JM)+(ZH(J)-ZA(JM))*
     >                IPSD(JM+1))/(ZA(JM+1)-ZA(JM))
           
         ENDDO

      ELSE    
         DO J=1,NH
            WC(1,J)=IWC(J)
            WC(2,J)=LWC(J)
            PSDF(J)=IPSD(J)
         ENDDO
      ENDIF

      DO I=1,NH
         CHK_CLD(I)=WC(1,I)+WC(2,I)
      ENDDO
C---------------------------------------------------------------------------

      RETURN
      END

! $Log: cloud_check.f,v      


