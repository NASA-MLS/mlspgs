
      SUBROUTINE MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMR,NZ,NS,
     >                       N, WCin, IPSDin,
     >                       YP,YZ,YT,YQ,VMR1,WC,NH,CHK_CLD,IPSD,
     >                       ZT,ZZT,NT)

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

      REAL WCin(N,NZ)                         
      INTEGER IPSDin(NZ)
      INTEGER NT                               ! NO. OF TANGENT PRESSURE LEVSLS
      REAL ZT(NT)                              ! TANGENT PRESSURE
      
C----------------------------------------------
C     OUTPUT PARAMETERS
C----------------------------------------------
      INTEGER NH                               ! MODEL ATMOSPHERIC LEVELS
      INTEGER NH0

      REAL YZ(NH)                              ! PRESSURE HEIGHT (m)
      REAL YP(NH)                              ! PRESSURE (hPa)
      REAL YT(NH)                              ! TEMPERATURE PROFILE
      REAL YQ(NH)                              ! RELATIVE HUMIDITY (%)
      REAL VMR1(NS,NH)                          ! 1=O3 VOLUME MIXING RATIO
      REAL WC(N,NH)
      INTEGER IPSD(NH)
      REAL CHK_CLD(NH)  ! CLOUD CHECKER      

      REAL ZZT(NT)                             ! TANGENT HEIGHT

C----------------------------------------------------
C     WORK SPACE
C----------------------------------------------------
      REAL PTOP,PBOTTOM,DP,ZH(NH),ZA(NH),ZZ(NH),WK
      INTEGER I,JM,J
C--------------------------------------------------------------------------

      PTOP = 40./16.-3.            ! TOP OF THE MODEL
      PBOTTOM=-ALOG10(PRESSURE(1)) ! BOTTOM OF THE MODEL
      DP=(PTOP-PBOTTOM)/NH         ! LAYER THICKNESS

      DO I=1,NH
         ZH(I)=PBOTTOM+(I-1)*DP
      END DO

      DO I=1,NZ
         ZA(I)=-ALOG10(PRESSURE(I))
         HEIGHT(I)= MAX ( 0., HEIGHT(I) )
      END DO

      IF (NZ .NE. NH) THEN

C==========================================
C     PRODUCE MODEL ATMOSPHERIC PROFILES
C==========================================

         DO J=1,NH

            CALL LOCATE (ZA,NZ,NH,ZH(J),JM)             
         
            YP(J)=((ZA(JM+1)-ZH(J))*PRESSURE(JM)+(ZH(J)-ZA(JM))*
     >            PRESSURE(JM+1))/(ZA(JM+1)-ZA(JM))             

            YZ(J)=((ZA(JM+1)-ZH(J))*HEIGHT(JM)+(ZH(J)-ZA(JM))*
     >            HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))

            YT(J)=((ZA(JM+1)-ZH(J))*TEMPERATURE(JM)+(ZH(J)-ZA(JM))*
     >            TEMPERATURE(JM+1))/(ZA(JM+1)-ZA(JM))

            WC(1,J)=((ZA(JM+1)-ZH(J))*WCin(1,JM)+(ZH(J)-ZA(JM))*
     >            WCin(1,JM+1))/(ZA(JM+1)-ZA(JM))

            WC(2,J)=((ZA(JM+1)-ZH(J))*WCin(2,JM)+(ZH(J)-ZA(JM))*
     >            WCin(2,JM+1))/(ZA(JM+1)-ZA(JM))
            
            YQ(J)=((ZA(JM+1)-ZH(J))*VMR(1,JM)+(ZH(J)-ZA(JM))*
     >            VMR(1,JM+1))/(ZA(JM+1)-ZA(JM))

            VMR1(1,J)=((ZA(JM+1)-ZH(J))*VMR(2,JM)+(ZH(J)-ZA(JM))*
     >                VMR(2,JM+1))/(ZA(JM+1)-ZA(JM))
         
            IPSD(J)=((ZA(JM+1)-ZH(J))*IPSDin(JM)+(ZH(J)-ZA(JM))*
     >            IPSDin(JM+1))/(ZA(JM+1)-ZA(JM))

            CHK_CLD(J) = WC(1,J) + WC(2,J)

         ENDDO

      ELSE
         
         DO J=1,NZ
            YP(J)     = PRESSURE(J)    
            YZ(J)     = HEIGHT(J)      
            YT(J)     = TEMPERATURE(J)
            WC(1,J)   = WCin(1,J)
            WC(2,J)   = WCin(2,J)
            YQ(J)     = VMR(1,J)     
            VMR1(1,J) = VMR(2,J)    
            IPSD(J)   = IPSDin(J)
            CHK_CLD(J) = WC(1,J) + WC(2,J)
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






