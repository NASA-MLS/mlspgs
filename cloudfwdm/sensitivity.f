
      SUBROUTINE SENSITIVITY(DTcir,ZT,NT,YP,YZ,NH,PRESSURE,NZ,
     >                       delTAU,delTAUc,TAUeff,SS,
     >                       N,NF,IRF,ISWI,RE)
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

      REAL*8 RE
      REAL HT,C_EXT,A_EXT,TGT,DS,DTAU,A_COL
      INTEGER I,K,J,iflag
C-----------------------------------------------------------------------------

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









