      PROGRAM CLOUD_MODEL

C===========================================================================C
C     PROGRAM TO GENERATE CLOUD PARAMETERS FOR FULL CLOUD FORWARD MODEL     C
C---------------------------------------------------------------------------C
C                                                                           C
C     THIS PROGRAM GENERATES PROFILES OF CLOUD ICE-WATER-CONTENT (IWC),     C
C     LIQUID-WATER-CONTECT (LWC), AND THE SIZE-DISTRIBUTION-FLAG(IPSD),     C
C     WHICH ARE REQUIRED AS INPUT PARAMETERS FOR THE FULL CLOUD FORWARD     C
C     MODEL IN LEVEL 2 DATA PROCESSING.                                     C
C                                                                           C
C     JONATHAN H. JIANG, JUNE 6, 2001                                       C
C                                                                           C
C---------------------------------------------------------------------------C
C                                                                           C
C     --- USER DEFINED PARAMETERS ---                                       C
C                                                                           C
C     ITYPE  :   CLOUD TYPE                                                 C
C     CHT    :   CLOUD HEIGHT SHIFT                                         C
C     WCscale:   WATER CONTENT SCALE FACTOR                                 C
C                                                                           C
C     <<< INPUT PARAMETERS >>>                                              C
C                                                                           C
C     NZ     :   NUMBER OF VERTICAL LEVELS                                  C
C     YP (NZ):   PRESSURE (hPa)                                             C
C     YZ (NZ):   GEOPOTENTIAL HEIGHT (km)                                   C
C                                                                           C
C     >>> OUTPUT PARAMETERS <<<                                             C
C     -------------------------                                             C
C     IWC  (NZ): CLOUD ICE-WATER-CONTENT (g/m3)                             C
C     LWC  (NZ): CLOUD LIQUID-WATER-CONTENT (g/m3)                          C
C     IPSD (NZ): SIZE-DISTRIBUTION-FLAG (4 digit integer IILL)              C
C                                                                           C
C     Note: In current setup, IPSD is also defined by the user.             C
C                                                                           C
C===========================================================================C

      INTEGER NZ

      PARAMETER (NZ=320)

      INTEGER ITYPE                            ! CLOUD MODEL TYPE
                                               ! 0 = NO CLOUDS
                                               ! 1 = DEEP-CONVECTIVE
                                               ! 2 = FRONTAL
                                               ! 3 = ANVIL
                                               ! 4 = THIN-LAYER
 
      INTEGER  PSD                             ! SIZE-DISTRIBUTION FLAG
                                               ! IILL     (I=ICE, L=LIQUID)
                                               ! 1000:     I->MH, L->GAMMA 
                                               ! 1100:     I->LIU-CURRY
                                               ! 2000-3900:I->MODIFIED GAMMA
                                               !              WITH VARIOUS De,
                                               !              alpha
                                               ! 4000-5900:I->KNOLLENBERG WITH 
                                               !              VARIOUS b1
                                               ! 6000:     I->PSD FOR PSC

      REAL CHT                                 ! CLOUD HEIGHT INDICATER
                                               ! CHT=0.0  DEFAULT PROFILE
                                               ! CHT=1.0  MOVE CLOUD 1KM UP
                                               ! CHT=2.0  MOVE CLOUD 1KM DOWN

      REAL WCscale                             ! CLOUD WATER SCALING FACTOR
      
      PARAMETER (ITYPE=1, CHT=0., WCscale=0.1, PSD=1000)

      INTEGER IPSD (NZ)

      REAL YP   (NZ)
      REAL YZ   (NZ)
      REAL IWC  (NZ)
      REAL LWC  (NZ)
 
      REAL CLD_TOP, CLD_BASE, UPPER_LAG, LOWER_LAG, HT

      CHARACTER OUTFL*80, TITLE*80

C--------------------------------------------------------------------

      OPEN(91,FILE='L2.dat')                     ! INPUT L2 PRESSURE LEVELS 
      OPEN(92,FILE='CLOUD.dat',STATUS='UNKNOWN') ! OUTPUT L2 CLOUD MODEL FILE

      READ(91,'(a80)') TITLE

      DO I=1,NZ
         READ (91,*) YZ(I),YP(I)
      ENDDO
      CLOSE(91)

      HT = 16000. + CHT*1000.

      IF(ITYPE .EQ. 0) THEN
         PRINT*,' '
         PRINT*,' -> NO CLOUDS ! '
         GO TO 100
      ELSE IF(ITYPE .EQ. 1) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: DEEP-CONVECTIVE SYSTEM '
C                                ======================

         CLD_TOP   = HT + 1000.         ! CONVECTIVE CLOUD-TOP
         CLD_BASE  = HT - 11600.        ! CONVECTIVE CLOUD-BASE
         UPPER_LAG = HT - 8000.
         LOWER_LAG = HT - 10000.

         DO I=1,NZ

            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               IWC(I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               IWC(I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               IWC(I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/5000.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               IWC(I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF

            ! LIQUID WATER CLOUD
            
            IF (YZ(I).GT.5000. .AND. YZ(I).LT. 10000.) THEN
               LWC(I) = 0.1*0.
            ENDIF
         ENDDO


      ELSE IF (ITYPE .EQ. 2) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: FRONTAL SYSTEM '
C                                ==============
         CLD_TOP   = HT - 5000.         !FRONTAL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !FRONTAL CLOUD-BASE
         UPPER_LAG = HT - 7000.
         LOWER_LAG = HT - 9000.

         DO I=1,NZ
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               IWC(I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               IWC(I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               IWC(I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               IWC(I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF
            LWC(I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 3) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: ANVILS'
C                                ======
         CLD_TOP   = HT - 4000.         !ANVIL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !ANVIL CLOUD-BASE
         UPPER_LAG = HT - 5000.
         LOWER_LAG = HT - 6000.

         DO I=1,NZ
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               IWC(I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               IWC(I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               IWC(I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               IWC(I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/600.)
            ENDIF
            LWC(I) = 0.0
         ENDDO

      ELSE 
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: THIN-LAYER CIRRUS'
C                                =================
         CLD_TOP   = HT + 500.         !1KM THICK CLOUD LAYER
         CLD_BASE  = HT - 500.
         UPPER_LAG = HT + 200.
         LOWER_LAG = HT - 200.

         DO I=1,NZ
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               IWC(I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               IWC(I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               IWC(I) = 1.0*EXP(-(YZ(I)-HT)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               IWC(I) = 1.0*EXP(-(HT-YZ(I))/500.)
            ENDIF
            LWC(I)=0.0
         ENDDO
      ENDIF


      DO I=1,NZ
            IWC(I)=WCscale*IWC(I)
            LWC(I)=WCscale*LWC(I)
            IPSD(I)=1010
      ENDDO

      GOTO 200

 100  DO I=1,NZ
            IWC(I)=0.
            LWC(I)=0.
            IPSD(I)=0
      ENDDO

 200  CONTINUE


      WRITE(92,*) 'Pres(mb)     Ht(km) IWC(g/m3) LWC(g/m3)      IPSD'   
      DO I=1,NZ
         WRITE(92,'(f9.4,1x,f9.2,3x,f6.4,4x,f6.4,I)') 
     >              YP(I),  YZ(I),  IWC(I), LWC(I),IPSD(I)  
      ENDDO
      CLOSE(92)
C------------------------------------------------------

      STOP
      END

! $Log: cloud_model.f,v      



