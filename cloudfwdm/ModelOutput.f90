
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module ModelOutput

! ---------------------------------------------------------------------------  
! OUTPUT CLOUD EXTINCTION, EFFECTIVE OPTICAL DEPTH, AND RADIANCE SENSITIVITY
! ---------------------------------------------------------------------------

      use Interpack, only: LOCATE
      use MLSCommon, only: r8
      use MLSNumerics, only: INTERPOLATEVALUES
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
        &                    N,RE, noS, Slevl)
!----------------------------------------------------------------

      INTEGER :: NH, NZ, NT, N, noS

      REAL(r8) :: ZT(NT)                            ! TANGENT HEIGHT
      
      REAL(r8) :: YP(NH)                            ! MODEL PRESSURE LEVEL
      REAL(r8) :: YZ(NH)                            ! MODEL PRESSURE HEIGHT
      REAL(r8) :: Slevl(noS)                        ! Sgrid level (m)

      REAL(r8) :: PRESSURE(NZ)                      ! L2 PRESSURE LEVEL
      REAL(r8) :: DTcir(NT)                         ! CLOUD-INDUCED RADIANCE
      REAL(r8) :: SS(NT)                            ! CLOUD RADIANCE SENSITIVITY
      REAL(r8) :: TAUeff(NT)                        ! CLOUD EFFECTIVE OPTICAL DEPTH

      REAL(r8) :: Trans(2*NH)                       ! Clear Transmission Func 

      REAL(r8) :: delTAU100(NH-1)                   ! 100% AIR optical thickness 
      REAL(r8) :: delTAU(NH-1)                      ! TOTAL sky optical thickness 
      REAL(r8) :: delTAUc(NH-1)                     ! CLOUDY-SKY optical thickness

      REAL(r8) :: Trans_out(noS,NT)                 ! TOTAL Clear Trans Func
      REAL(r8) :: BETA(NZ)                        ! TOTAL EXTINCTION
      REAL(r8) :: BETAc(NZ)                       ! CLOUDY-SKY EXTINCTION
      REAL(r8) :: DDm(N,NH-1)                       ! MASS-MEAN-DIAMETER
      REAL(r8) :: Dm(N,NZ)                        ! MASS-MEAN-DIAMETER
      REAL(r8) :: DDZ(NH-1)                         ! MODEL LEYER THICKNESS
      REAL(r8) :: DZ(NZ-1)                          ! L2 LAYER THICKNESS

      REAL(r8) :: RE, xout, sum
      REAL(r8) :: HT,C_EXT,A_EXT,TGT,DS,DTAU,A_COL
      REAL(r8) :: ZH(NH-1),ZA(NZ-1), S(2*NH)
      INTEGER :: I,K,J,iflag, JM
!-----------------------------------------------------------------------------

!===============================================================
!     INTERPOLATE PARAMETERS FROM MODEL LEVEL NH TO L2 LEVEL NZ
!===============================================================
! COMPUTE TRANSMISSION FUNCTION

      trans = 0._r8            ! initialization

      do k=1,nt
        sum = 0._r8
        s = 0._r8    
            s(1) = - sqrt((re+yz(NH))**2-(re+zt(k))**2)
            s(2*NH) = sqrt((re+yz(NH))**2-(re+zt(k))**2)
        DO I=1,NH

          if(yz(NH-i+1) > zt(k) .and. i > 1) then
            ! s will zero for heights below zt
            s(i) = - sqrt((re+yz(NH-i+1))**2-(re+zt(k))**2)
            ds = s(i) - s(i-1)
            sum = sum +  deltau100(NH-i+1)*ds/(yz(NH-i+2)-yz(NH-i+1))
          endif
            trans(i) = sum     ! trans keeps constant for heights below zt
        ENDDO

        DO I=1,NH
          if(yz(i) > zt(k) .and. i < NH) then
            s(NH+i) = sqrt((re+yz(i))**2-(re+zt(k))**2)
            ds = s(NH+i) - s(NH+i-1)
            sum = sum +  deltau100(i)*ds/(yz(i+1)-yz(i))
          endif
            trans(NH+i) = sum
        ENDDO

        where(trans < 10._r8) 
         trans = exp(-trans)
        elsewhere 
         trans = 0._r8
        endwhere
        
      ! sLevl is defined differently from s
      CALL INTERPOLATEVALUES(s,trans, -Slevl,trans_out(:,k),method='Linear')
      enddo

! CONVERT DELTAU TO BETA
      DO I=1,NH-1
         ZH(I)=-(LOG10(YP(I+1))+LOG10(YP(I)))/2._r8
         delTAU(I) = delTAU(I)/DDZ(I)
         delTAUc(I) = delTAUc(I)/DDZ(I)
      END DO

      DO I=1,NZ-1
         ZA(I)=-(LOG10(PRESSURE(I+1))+LOG10(PRESSURE(I)))/2._r8
      END DO

      call INTERPOLATEVALUES(zh,deltau,za,beta(1:NZ-1),method='Linear')
      call INTERPOLATEVALUES(zh,deltauc,za,betac(1:NZ-1),method='Linear')
      call INTERPOLATEVALUES(zh,ddz,za,dz(1:NZ-1),method='Linear')
      call INTERPOLATEVALUES(zh,reshape(ddm(1,:),(/nh-1/)),&
        & za,dm(1,1:NZ-1),method='Linear')
      call INTERPOLATEVALUES(zh,reshape(ddm(2,:),(/nh-1/)),&
        & za,dm(2,1:NZ-1),method='Linear')

!==========================================================================
!     RADIANCE SENSITIVITY CALCULATIONS
!==========================================================================

      TAUeff = 0._r8
      SS = 0._r8
      
      DO I=1,NT

            HT = ZT(I)
            if(HT .gt. yz(nh-1)) cycle    ! don't compute for large Zt  
            
               do j=2,NH
                  if (yz(j) .GE. ht) then
                     iflag=j
                     goto 100
                  endif
               enddo

 100           continue

               C_EXT = 0._r8
               A_EXT = 0._r8
            
               DO K=NH-1,iflag,-1

                     TGT=YZ(K-1)
                     IF(TGT .LT. HT)THEN
                        TGT=HT
                     ENDIF
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-    &
     &                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)*DDZ(k)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*DDZ(k)*EXP(-A_COL)*   &
     &                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU

               ENDDO

               DO K=iflag,NH-1

                     TGT=YZ(K-1)
                     IF(TGT .LT. HT)THEN
                        TGT=HT
                     ENDIF
                     DS = SQRT((RE+YZ(K))**2-(RE+HT)**2)-    &
     &                    SQRT((RE+TGT)**2-(RE+HT)**2)
                     DTAU = DS * delTAU(K)*DDZ(k)/(YZ(K)-YZ(K-1))
                     A_COL = A_EXT + DTAU/2.
                     C_EXT=C_EXT + delTAUc(K)*DDZ(k)*EXP(-A_COL)*   &
     &                     DS/(YZ(K)-YZ(K-1))
                     A_EXT = A_EXT + DTAU

               ENDDO

               
               TAUeff(I)=C_EXT
               if (TAUeff(I) .gt. 0._r8) then
                  SS(I)=DTcir(I)/TAUeff(I)
               else
                  SS(I)=0._r8
               endif

      ENDDO
      

!--------------------------------------------------------------------------

      END SUBROUTINE SENSITIVITY

end module ModelOutput

! $Log$
! Revision 1.13  2001/10/24 22:50:26  dwu
! revise teff calculation
!
! Revision 1.12  2001/10/24 22:09:53  dwu
! provent cases where tangent heights are greater than model top
!
! Revision 1.11  2001/10/19 19:39:49  dwu
! fix bug after changes
!
! Revision 1.10  2001/10/19 19:25:25  dwu
! follow-on fixes
!
! Revision 1.9  2001/10/19 19:16:14  dwu
! change Nz-1 to NZ
!
! Revision 1.8  2001/10/19 17:11:52  dwu
! revise teff calculation
!
! Revision 1.7  2001/10/12 22:38:24  dwu
! rewrite transmission function calculation
!
! Revision 1.6  2001/10/11 20:03:00  jonathan
! *** empty log message ***
!
! Revision 1.5  2001/10/10 02:49:35  dwu
! fix problem in trans calculation
!
! Revision 1.4  2001/10/09 22:12:36  jonathan
! fix trans func
!
! Revision 1.3  2001/09/21 15:51:37  jonathan
! modified F95 version
!
