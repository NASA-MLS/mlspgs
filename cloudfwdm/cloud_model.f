      SUBROUTINE CLOUD_MODEL(ITYPE,CHT,WC,CHK_CLD,YZ,NH)

C=========================================================================C
C  DEFINE VERTICAL PROFILES OF CLOUD ICE-WATER-CONTENT                    C
C  - J.JIANG (05/18/2001)                                                 C
C=========================================================================C

      INTEGER ITYPE,NH

      REAL WC(2,NH)     ! CLOUD WATER CONTENT (g/m3)
      REAL CHK_CLD(NH)  ! CLOUD CHECKER
      REAL YZ(NH)
      REAL CLD_TOP
      REAL CLD_BASE
      REAL UPPER_LAG
      REAL LOWER_LAG
      REAL HT,CHT
C--------------------------------------------------------------------

      HT = 16000. + CHT*1000.

      IF(ITYPE .EQ. 0) THEN
         PRINT*,' '
         PRINT*,' -> NO CLOUDS ! '
         GO TO 100
      ELSE IF(ITYPE .EQ. 1) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: DEEP-CONVECTIVE SYSTEM '
C                                ======================
         CLD_TOP   = HT + 1000.         !CONVECTIVE CLOUD-TOP
         CLD_BASE  = HT - 11600.        !CONVECTIVE CLOUD-BASE
         UPPER_LAG = HT - 8000.
         LOWER_LAG = HT - 10000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/5000.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 2) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: FRONTAL SYSTEM '
C                                ==============
         CLD_TOP   = HT - 5000.         !FRONTAL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !FRONTAL CLOUD-BASE
         UPPER_LAG = HT - 7000.
         LOWER_LAG = HT - 9000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 3) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: ANVILS'
C                                ======
         CLD_TOP   = HT - 4000.         !ANVIL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !ANVIL CLOUD-BASE
         UPPER_LAG = HT - 5000.
         LOWER_LAG = HT - 6000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/600.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE 
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: THIN-LAYER CIRRUS'
C                                =================
         CLD_TOP   = HT + 500.         !1KM THICK CLOUD LAYER
         CLD_BASE  = HT - 500.
         UPPER_LAG = HT + 200.
         LOWER_LAG = HT - 200.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-HT)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(HT-YZ(I))/500.)
            ENDIF
            WC(2,I)=0.0
         ENDDO
      ENDIF

      GOTO 200

 100  DO I=1,NH
            WC(2,I)=0.
            WC(1,I)=0.
      ENDDO

 200  DO I=1,NH
         CHK_CLD(I)=WC(1,I)+WC(2,I)
      ENDDO
      RETURN
      END



