
      SUBROUTINE SET_PARAM(PRESSURE,HEIGHT,TEMPERATURE,VMR,WC,
     >                       NZ,NS,N,YQ,VMR1,ZT,ZZT,NT,CHK_CLD)

C========================================================================C
C     DESCRIPTION                                                        C
C     -----------                                                        C
C     1: SET INTERNAL PARAMETERS
C     2: CONVERT TANGENT PRESSURE TO TANGENT HEIGHT                      C
C                                                                        C
C     J. JIANG, JUNE 6, 2001                                             C
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
      REAL WC(N,NZ)                            ! CLOUD WATER CONTENT (g/m3)
                                               ! (1=ICE,2=WATER)

      INTEGER NT                               ! NO. OF TANGENT PRESSURE LEVSLS
      REAL ZT(NT)                              ! TANGENT PRESSURE
      REAL CHK_CLD(NZ)  ! CLOUD CHECKER      
      
C----------------------------------------------
C     OUTPUT PARAMETERS
C----------------------------------------------

      REAL YQ(NZ)                              ! RELATIVE HUMIDITY (%)
      REAL VMR1(NS,NZ)                          ! 1=O3 VOLUME MIXING RATIO

      REAL ZZT(NT)                             ! TANGENT HEIGHT

C----------------------------------------------------
C     WORK SPACE
C----------------------------------------------------
      REAL PTOP,PBOTTOM,DP,WK,ZA(2000),ZZ(2000)
      INTEGER I,JM,J
C--------------------------------------------------------------------------

      DO J=1,NZ
         YQ(J)     = VMR(1,J)     
         VMR1(1,J) = VMR(2,J)    
         CHK_CLD(J)= WC(1,J) + WC(2,J)
      ENDDO


C=====================================================
C     CONVERT TANGENT PRESSURE TO TANGENT HEIGHTS (KM)
C=====================================================

      DO I=1,NT
         ZZ(I)=-ALOG10(ZT(I))
      END DO

      DO I=1,NZ
         ZA(I)=-ALOG10(PRESSURE(I))
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


