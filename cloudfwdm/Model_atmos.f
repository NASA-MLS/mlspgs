
      SUBROUTINE MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMR,NZ,
     2                       YP,YZ,YT,YQ,VMR1,NH)

C========================================================================C
C     DESCRIPTION OF PARAMETERS                                          C
C     -------------------------                                          C
C     1: INPUT L2 ATMOSPHERIC PARAMETERS (VERTICAL PROFILES).            C
C     2: OUTPUT INTERPOLATED ATMOSPHERIC PARAMETERS FOR THE CLOUD        C 
C        FORWARD MODEL. THE DEFAULT SETTING IS 640 MODEL LAYERS FROM     C
C        SURFACE TO ~40KM.                                               C
C                                                                        C
C     J. JIANG, MAY 18, 2001                                             C
C========================================================================C

C----------------------------------------------
C     INPUT PARAMETERS
C----------------------------------------------
      INTEGER NZ                               ! L2 ATMOSPHERIC LEVELS
      REAL PRESSURE(NZ)                        ! PRESSURE LEVEL
      REAL HEIGHT(NZ)                          ! PRESSURE HEIGHT
      REAL TEMPERATURE(NZ)                     ! ATMOSPHERIC TEMPERATURE
      REAL VMR(5,NZ)                           ! 1=O3 VOLUME MIXING RATIO
                                               ! 2=H2O
C----------------------------------------------
C     INPUT PARAMETERS
C----------------------------------------------
      INTEGER NH                               ! MODEL ATMOSPHERIC LEVELS
      INTEGER NH0

      PARAMETER(NH0=2000)

      REAL YZ(NH)                             ! PRESSURE HEIGHT (m)
      REAL YP(NH)                             ! PRESSURE (hPa)
      REAL YT(NH)                             ! TEMPERATURE PROFILE
      REAL YQ(NH)                             ! RELATIVE HUMIDITY (%)
      REAL VMR1(5,NH)                         ! 1=O3 VOLUME MIXING RATIO

C----------------------------------------------
C     WORK SPACE
C----------------------------------------------
      REAL PTOP,PBOTTOM,DP,ZH(NH0),ZA(NH0),WK
      INTEGER I,JM
C--------------------------------------------------------------------------

      PTOP = 40./16.-3.                         ! TOP OF THE MODEL
      PBOTTOM=-ALOG10(PRESSURE(1))              ! BOTTOM OF THE MODEL
      DP=(PTOP-PBOTTOM)/NH                      ! LAYER THICKNESS

      DO I=1,NH
         ZH(I)=PBOTTOM+(I-1)*DP
      END DO

      DO I=1,NZ
         ZA(I)=-ALOG10(PRESSURE(I))
      END DO

      DO J=1,NH

         CALL LOCATE (ZA,NZ,NH,ZH(J),JM)             
         
         YP(J)=((ZA(JM+1)-ZH(J))*PRESSURE(JM)+(ZH(J)-ZA(JM))*
     >         PRESSURE(JM+1))/(ZA(JM+1)-ZA(JM))             

         YZ(J)=((ZA(JM+1)-ZH(J))*HEIGHT(JM)+(ZH(J)-ZA(JM))*
     >         HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))

         YT(J)=((ZA(JM+1)-ZH(J))*TEMPERATURE(JM)+(ZH(J)-ZA(JM))*
     >         TEMPERATURE(JM+1))/(ZA(JM+1)-ZA(JM))

         YQ(J)=((ZA(JM+1)-ZH(J))*VMR(1,JM)+(ZH(J)-ZA(JM))*
     >         VMR(1,JM+1))/(ZA(JM+1)-ZA(JM))

         VMR1(1,J)=((ZA(JM+1)-ZH(J))*VMR(2,JM)+(ZH(J)-ZA(JM))*
     >         VMR(2,JM+1))/(ZA(JM+1)-ZA(JM))

      ENDDO

C----------------------------------------------------------------------

      RETURN
      END


