
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module ModelOutput

! ---------------------------------------------------------------------------  
! OUTPUT CLOUD EXTINCTION, EFFECTIVE OPTICAL DEPTH, AND RADIANCE SENSITIVITY
! ---------------------------------------------------------------------------

      use Interpack, only: LOCATE
      use MLSCommon, only: r8
      IMPLICIT NONE
      private
      public :: SENSITIVITY

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------

contains 

      SUBROUTINE SENSITIVITY(DTcir,ZT,NT,YP,YZ,NH,PRESSURE,NZ,   &
        &                    delTAU,delTAUc,delTAU100,TAUeff,SS, &
        &                    Trans_out,BETA,BETAc,DDm,Dm,DDZ,DZ, &
        &                    N,ISWI,RE)
!----------------------------------------------------------------

      INTEGER :: NH, NZ, NT, N

      REAL(r8) :: ZT(NT)                            ! TANGENT HEIGHT
      REAL(r8) :: YP(NH)                            ! MODEL PRESSURE LEVEL
      REAL(r8) :: YZ(NH)                            ! MODEL PRESSURE HEIGHT
      REAL(r8) :: PRESSURE(NZ)                      ! L2 PRESSURE LEVEL
      REAL(r8) :: DTcir(NT)                         ! CLOUD-INDUCED RADIANCE
      REAL(r8) :: SS(NT)                            ! CLOUD RADIANCE SENSITIVITY
      REAL(r8) :: TAUeff(NT)                        ! CLOUD EFFECTIVE OPTICAL DEPTH

      REAL(r8) :: Trans(NH-1)                       ! Clear Transmission Func 
      REAL(r8) :: delTAU100(NH-1)                   ! 100% AIR EXTINCTION 
      REAL(r8) :: delTAU(NH-1)                      ! TOTAL EXTINCTION 
      REAL(r8) :: delTAUc(NH-1)                     ! CLOUDY-SKY EXTINCTION

      REAL(r8) :: Trans_out(NZ-1)                   ! TOTAL Clear Trans Func
      REAL(r8) :: BETA(NZ-1)                        ! TOTAL EXTINCTION
      REAL(r8) :: BETAc(NZ-1)                       ! CLOUDY-SKY EXTINCTION
      REAL(r8) :: DDm(N,NH-1)                       ! MASS-MEAN-DIAMETER
      REAL(r8) :: Dm(N,NZ-1)                        ! MASS-MEAN-DIAMETER
      REAL(r8) :: DDZ(NH-1)                         ! MODEL LEYER THICKNESS
      REAL(r8) :: DZ(NZ-1)                          ! L2 LAYER THICKNESS

      REAL(r8) :: RE
      REAL(r8) :: HT,C_EXT,A_EXT,TGT,DS,DTAU,A_COL
      REAL(r8) :: ZH(NH-1),ZA(NH-1)
      INTEGER :: I,K,J,iflag, JM, ISWI
!-----------------------------------------------------------------------------

!===============================================================
!     INTERPOLATE PARAMETERS FROM MODEL LEVEL NH TO L2 LEVEL NZ
!===============================================================
! COMPUTE TRANSMISSION FUNCTION
      DO I=1,NH-1
        TRANS(I)=0._r8
      ENDDO
        TRANS(NH-1) = DELTAU100(NH-1)
      DO I=NH-2,1,-1
        TRANS(I)=TRANS(I+1)+DELTAU100(I)   ! only clear sky absorption
      ENDDO

      Trans = exp(- Trans)

! CONVERT DELTAU TO BETA
      DO I=1,NH-1
         ZH(I)=-(LOG10(YP(I+1))+LOG10(YP(I)))/2.
         delTAU(I) = delTAU(I)/DDZ(I)
         delTAUc(I) = delTAUc(I)/DDZ(I)
      END DO

      DO I=1,NZ-1
         ZA(I)=-(LOG10(PRESSURE(I+1))+LOG10(PRESSURE(I)))/2.
      END DO

      DO J=1,NZ-1
      
         CALL LOCATE (ZH,NH-1, NH-1, ZA(J),JM)
         
         Trans_out(J)=((ZH(JM+1)-ZA(J))*TRANS(JM)+(ZA(J)-ZH(JM))*   &
     &                TRANS(JM+1))/(ZH(JM+1)-ZH(JM))
     
         BETA(J)=((ZH(JM+1)-ZA(J))*delTAU(JM)+(ZA(J)-ZH(JM))*   &
     &                delTAU(JM+1))/(ZH(JM+1)-ZH(JM))

         BETAc(J)=((ZH(JM+1)-ZA(J))*delTAUc(JM)+(ZA(J)-ZH(JM))* &
     &                delTAUc(JM+1))/(ZH(JM+1)-ZH(JM))             


         Dm(1,J)=((ZH(JM+1)-ZA(J))*DDm(1,JM)+(ZA(J)-ZH(JM))*        &
     &                DDm(1,JM+1))/(ZH(JM+1)-ZH(JM))             

         Dm(2,J)=((ZH(JM+1)-ZA(J))*DDm(2,JM)+(ZA(J)-ZH(JM))*        &
     &                DDm(2,JM+1))/(ZH(JM+1)-ZH(JM))             

         DZ(J)=((ZH(JM+1)-ZA(J))*DDZ(JM)+(ZA(J)-ZH(JM))*            &
     &                DDZ(JM+1))/(ZH(JM+1)-ZH(JM))             

      ENDDO


!==========================================================================
!     RADIANCE SENSITIVITY CALCULATIONS
!==========================================================================

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
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-    &
     &                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)*DDZ(k)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*DDZ(k)*EXP(-A_COL)*   &
     &                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU
                  ENDIF

               ENDDO

               DO K=iflag+1,NH-1

                  IF (YZ(K) .LT. HT) THEN
                     TGT=YZ(K-1)
                     IF(YZ(K).EQ.HT)THEN
                        TGT=HT
                     ENDIF 
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-    &
     &                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)*DDZ(k)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*DDZ(k)*EXP(-A_COL)*   &
     &                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU
                  ENDIF

               ENDDO

               
               TAUeff(I)=C_EXT
               if (TAUeff(I) .gt. 0.) then
                  SS(I)=DTcir(I)/TAUeff(I)
               else
                  SS(I)=0.
               endif
            ENDIF

         ENDDO
      
      ENDIF

!--------------------------------------------------------------------------

      END SUBROUTINE SENSITIVITY

end module ModelOutput

! $Log: ModelOutput.f90,v      
