
      SUBROUTINE SENSITIVITY(DTcir,ZT,NT,YP,YZ,NH,PRESSURE,NZ,
     >                       delTAU,delTAUc,DTAUt,DTAUc,TAUeff,SS,
     >                       DDm,Dm,N,NF,IRF,ISWI,RE)
C----------------------------------------------------------------

      INTEGER NH,NZ

      REAL ZT(NT)                              ! TANGENT HEIGHT
      REAL YP(NH)                              ! MODEL PRESSURE LEVEL
      REAL YZ(NH)                              ! MODEL PRESSURE HEIGHT
      REAL PRESSURE(NZ)                        ! L2 PRESSURE LEVEL
      REAL DTcir(NT+1,NF)                      ! CLOUD-INDUCED RADIANCE
      REAL SS(NT+1,NF)                         ! CLOUD RADIANCE SENSITIVITY
      REAL TAUeff(NT+1,NF)                     ! CLOUD EFFECTIVE OPTICAL DEPTH

      REAL delTAU(NH-1)                        ! CLEAR-SKY  
      REAL delTAUc(NH-1)                       ! CLOUDY-SKY EXTINCTION
      REAL DTAUt(NZ,NF)                        ! TOTAL EXTINCTION
      REAL DTAUc(NZ,NF)                        ! CLOUDY-SKY EXTINCTION

      REAL DDm(N,NH-1)                         ! MASS-MEAN-DIAMETER
      REAL Dm(N,NZ-1)                          ! MASS-MEAN-DIAMETER

      REAL*8 RE
      REAL HT,C_EXT,A_EXT,TGT,DS,DTAU,A_COL
      REAL ZH(2000),ZA(2000)
      INTEGER I,K,J,iflag
C-----------------------------------------------------------------------------

C===============================================================
C     INTERPOLATE PARAMETERS FROM MODEL LEVEL NH TO L2 LEVEL NZ
C===============================================================

      DO I=1,NH-1
         ZH(I)=-ALOG10( (YP(I+1)+YP(I))/2. )
      END DO

      DO I=1,NZ-1
         ZA(I)=-ALOG10( (PRESSURE(I+1)+PRESSURE(I))/2. )
      END DO

      DO J=1,NZ-1
      
         CALL LOCATE (ZH,NH,NZ,ZA(J),JM)
         
         DTAUt(J,IRF)=((ZH(JM+1)-ZA(J))*delTAU(JM)+(ZA(J)-ZH(JM))*
     >                delTAU(JM+1))/(ZH(JM+1)-ZH(JM))             

         DTAUc(J,IRF)=((ZH(JM+1)-ZA(J))*delTAUc(JM)+(ZA(J)-ZH(JM))*
     >                delTAUc(JM+1))/(ZH(JM+1)-ZH(JM))             

         Dm(1,J)=((ZH(JM+1)-ZA(J))*DDm(1,JM)+(ZA(J)-ZH(JM))*
     >                DDm(1,JM+1))/(ZH(JM+1)-ZH(JM))             

         Dm(2,J)=((ZH(JM+1)-ZA(J))*DDm(2,JM)+(ZA(J)-ZH(JM))*
     >                DDm(2,JM+1))/(ZH(JM+1)-ZH(JM))             

      ENDDO

C==========================================================================
C     RADIANCE SENSITIVITY CALCULATIONS
C==========================================================================

      IF (ISWI .EQ. 0) THEN

         DO I=1,NT

            HT = ZT(I)   

            IF (HT .GE. 0.) THEN

               do j=1,nh
                  if (yz(j) .ge. ht) then
                     iflag=j
                     goto 100
                  endif
               enddo

 100           continue

               C_EXT = 0.
               A_EXT = 0.
            
               DO K=NH-1,iflag+1,-1

                  IF (YZ(K) .GT. HT) THEN
                     TGT=YZ(K-1) 
                     IF(YZ(K).EQ.HT)THEN
                        TGT=HT
                     ENDIF
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-
     >                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*EXP(-A_COL)*
     >                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU
                  ENDIF

               ENDDO

               DO K=iflag+1,NH-1

                  IF (YZ(K) .LT. HT) THEN
                     TGT=YZ(K-1)
                     IF(YZ(K).EQ.HT)THEN
                        TGT=HT
                     ENDIF 
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-
     >                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*EXP(-A_COL)*
     >                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU
                  ENDIF

               ENDDO

               TAUeff(I,IRF)=C_EXT
               SS(I,IRF)=DTcir(I,IRF)/TAUeff(I,IRF)

            ENDIF
         ENDDO
      
      ENDIF

C--------------------------------------------------------------------------
      RETURN
      END

! $Log: sensitivity.f,v      



