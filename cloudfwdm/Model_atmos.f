
      SUBROUTINE MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMR,NZ,NS,
     >                       YP,YZ,YT,YQ,VMR1,NH,ZT,ZZT,NT)

C========================================================================C
C     DESCRIPTION                                                        C
C     -----------                                                        C
C     1: CHECK THE INPUT L2 ATMOSPHERIC PROFILES.                        C
C     2: OUTPUT INTERPOLATED ATMOSPHERIC PROFILES IF NEEDED TO ENSURE    C
C        THERE ARE 640 MODEL LAYERS FROM SURFACE TO ~40KM.               C
C     3: CONVERT TANGENT PRESSURE TO TANGENT HEIGHT                      C
C                                                                        C
C     J. JIANG, MAY 18, 2001                                             C
C========================================================================C

C----------------------------------------------
C     INPUT PARAMETERS
C----------------------------------------------

      INTEGER NZ                               ! NO. OF L2 ATMOSPHERIC LEVELS
      REAL PRESSURE(NZ)                        ! PRESSURE LEVEL
      REAL HEIGHT(NZ)                          ! PRESSURE HEIGHT
      REAL TEMPERATURE(NZ)                     ! ATMOSPHERIC TEMPERATURE
      REAL VMR(NS,NZ)                          ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3
      INTEGER NT                               ! NO. OF TANGENT PRESSURE LEVSLS
      REAL ZT(NT)                              ! TANGENT PRESSURE
      
C----------------------------------------------
C     OUTPUT PARAMETERS
C----------------------------------------------
      INTEGER NH                               ! MODEL ATMOSPHERIC LEVELS
      INTEGER NH0

      PARAMETER(NH0=2000)

      REAL YZ(NH)                              ! PRESSURE HEIGHT (m)
      REAL YP(NH)                              ! PRESSURE (hPa)
      REAL YT(NH)                              ! TEMPERATURE PROFILE
      REAL YQ(NH)                              ! RELATIVE HUMIDITY (%)
      REAL VMR1(5,NH)                          ! 1=O3 VOLUME MIXING RATIO

      REAL ZZT(NT)                             ! TANGENT HEIGHT

C----------------------------------------------------
C     WORK SPACE
C----------------------------------------------------
      REAL PTOP,PBOTTOM,DP,ZH(NH0),ZA(NH0),ZZ(NH0),WK
      INTEGER I,JM,J
C--------------------------------------------------------------------------


      IF (NZ .NE. NH) THEN

C==========================================
C     PRODUCE MODEL ATMOSPHERIC PROFILES
C==========================================

         PTOP = 40./16.-3.      ! TOP OF THE MODEL
         PBOTTOM=-ALOG10(PRESSURE(1)) ! BOTTOM OF THE MODEL
         DP=(PTOP-PBOTTOM)/NH   ! LAYER THICKNESS

         DO I=1,NH
            ZH(I)=PBOTTOM+(I-1)*DP
         END DO

         DO I=1,NZ
            ZA(I)=-ALOG10(PRESSURE(I))
         END DO

         DO J=1,NH

            CALL LOCATE (ZA,NZ,NH,ZH(J),JM)             
         
            YP(J)=((ZA(JM+1)-ZH(J))*PRESSURE(JM)+(ZH(J)-ZA(JM))*
     >            PRESSURE(JM+1))/(ZA(JM+1)-ZA(JM))             

            YZ(J)=((ZA(JM+1)-ZH(J))*HEIGHT(JM)+(ZH(J)-ZA(JM))*
     >            HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))

            YT(J)=((ZA(JM+1)-ZH(J))*TEMPERATURE(JM)+(ZH(J)-ZA(JM))*
     >            TEMPERATURE(JM+1))/(ZA(JM+1)-ZA(JM))

            YQ(J)=((ZA(JM+1)-ZH(J))*VMR(1,JM)+(ZH(J)-ZA(JM))*
     >            VMR(1,JM+1))/(ZA(JM+1)-ZA(JM))

            VMR1(1,J)=((ZA(JM+1)-ZH(J))*VMR(2,JM)+(ZH(J)-ZA(JM))*
     >                VMR(2,JM+1))/(ZA(JM+1)-ZA(JM))

         ENDDO

      ELSE
         
         DO J=1,NZ
            YP(J) = PRESSURE(J)    
            YZ(J) = HEIGHT(J)      
            YT(J) = TEMPERATURE(J)
            YQ(J) = VMR(1,J)     
            VMR1(1,J)=VMR(2,J)    
         ENDDO

      ENDIF

C----------------------------------------------------------------------

C==========================================
C     PRODUCE TANGENT HEIGHTS (KM)
C==========================================

      DO I=1,NT
         ZZ(I)=-ALOG10(ZT(I))
      END DO

      DO J=1,NT
      
         CALL LOCATE (ZA,NZ,NT,ZZ(J),JM)

         ZZT(J)=((ZA(JM+1)-ZZ(J))*HEIGHT(JM)+(ZZ(J)-ZA(JM))*
     >                HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))             
      ENDDO

C----------------------------------------------------------------------

      RETURN
      END

! $Log: Model_atmos.f,v      


