! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RadiativeTransferModule

! -------------------------------------------------------------------------  
! PERFORM MLS LIB RADIATIVE TRANSFER CALCULATIONS
! -------------------------------------------------------------------------

      use Blackbody, only: planck
      use Interpack, only: LOCATE
      use MLSCommon, only: r8 
      use Units     
      IMPLICIT NONE
      Private
      Public :: RADXFER

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE RADXFER(L,NU,NUA,U,DU,PHH,NT,ZT,W0,dTAU,RS,TS,FREQ, &
      &          YZ,TEMP_AIR,N,THETA,THETAI,PHI,UI,UA,TT,ICON,RE)

!===================================================================
!     >>>>>>ITERATIVE RADIATIVE TRANSFER CALCULATION<<<<<<
!
! METHOD: 1. GET SCATTERING TB IN A PARALLEL PLANE ASSUMING THAT
!            SCATTERING IS LOCALIZED AND PARALLEL PLANE APPROXI-
!            MATION IS OK;
!         2. EACH ITERATION INVOLVES INTERPOLATIONS OF TB BETWEEN 
!            NT-TANGENT HEIGHT GRIDS ON THE SPHERICAL PLANE AND
!            N-STREAMS ON THE PARALLEL PLANE;
!         3. OUTPUT TB FOR ALL TANGENT HEIGHTS
!
! LAST UPDATE: 18/05/2001, J.JIANG
!===================================================================

      use MLSCommon, only: r8
      IMPLICIT NONE
      REAL(r8) :: RE

      INTEGER :: L,NT,NU,NUA,N
      
      REAL(r8) :: FIRST                   ! FIRST GUESS OF UPWELLING TB
      PARAMETER (FIRST=250._r8)

      REAL(r8) :: U(NU),DU(NU),THETA(NU),PHI(NUA),      &
      &     THETAI(NU,NU,NUA),UI(NU,NU,NUA),UA(NUA) 

      REAL(r8) :: PHH(N,NU,L),RS(NU/2),W0(N,L)
      REAL(r8) :: TAU(L+1,2*L+1)          ! Integrated ptical depth at each LOS point for 
                                          ! each tangent height 

      REAL(r8) :: dTAU(L)                 ! Optical depth increment at each layer 

      REAL(r8) :: TEMP_AIR(L)             ! MEAN AIR TEMPERATURE AT L
      REAL(r8) :: YZ(L+1)                 ! PRESSURE HEIGHT (km)
      REAL(r8) :: ZT(NT)                  ! TANGENT HEIGHT (m)
      REAL(r8) :: TT(NT+1,L+1)            ! TB AT ZT, LAST INDEX IS ZENITH LOOKING
      REAL(r8) :: FREQ                    ! FREQUENCY (GHz)
      REAL(r8) :: TSPACE                  ! COSMIC BACKGROUND RADIANCE
      REAL(r8) :: TS                      ! SURFACE TEMPERATURE (K)
      REAL(r8) :: TEMP(L+1)               ! BRIGHTNESS TEMPERATURE FROM TEMP_AIR
      REAL(r8) :: TB0(NU)                 ! TB AT THE SURFACE
      REAL(r8) :: TAVG                    ! TB AVERAGED OVER PHI AT A GIVEN U
      REAL(r8) :: TB(NU,L+1)              ! TB IN FLAT PLANE GEOMETRY
                                          ! 1->NU/2 UPWELLING, NU/2->NU DOWNWELLING
      REAL(r8) :: TSCAT(N,NU,L+1)         ! TB FROM SCATTERING PHASE FUNCTION

      INTEGER :: LMIN(L+1)                ! LOWEST LAYER REACHED BY A TANGENT HT

      INTEGER :: ICON                     ! CONTROL SWITCH
                                          ! 0 = CLEAR-SKY
                                          ! 1 = Cloudy sky
                                          ! 2 = Cloudy sky
                                          ! 3 = NEAR SIDE CLOUD ONLY

      REAL(r8) :: UAVE(L+1,L+1)           ! INCIDENT ANGLES FOR EACH TANGENT HT
      REAL(r8) :: UEFF                    ! EFFECTIVE U BETWEEN K AND K+1

      REAL(r8) :: ITS0                    ! NO. OF MAXIMUM ITERATIONS
      REAL(r8) :: DELTA                   ! DELTA TB FOR CONVERGENCE CHECK

      INTEGER :: I,J,K,K1,ITS,ITT,IH,IP,JM,IND,ISPI
      REAL(r8) :: D1,D2,DY,TGT, jj0
      REAL(r8) :: WK,WK1,U1(NU),UU,X2,RSAVG,XTB(L+1),WW0,CHK(L+1)
      REAL(r8) :: tsource, wwk, wwk1,www0

      TB=0.0_r8
      TB0=0.0_r8
      TT=0.0_r8
      Tscat = 0.0_r8
      U1=0.0_r8
      TEMP=0.0_r8
      D1=0.0_r8
      D2=0.0_r8
      DY=0.0_r8
      TGT=0.0_r8
      uave=0.0_r8
      ueff=0.0_r8
      uu=0.0_r8
      JM=0
      WK=0.0_r8
      WK1=0.0_r8
      wwk=0.0_r8
      wwk1=0.0_r8
      www0=0.0_r8
      x2=0.0_r8
      RSAVG=0.0_r8
      XTB=0.0_r8
      WW0=0.0_r8
      tsource=0.0_r8
      lmin=0.0_r8

!------------------------------------------------
!     FIND BRIGHTNESS TEMPERATURE AT EACH LAYER
!------------------------------------------------
      DO K=1,L
         CALL PLANCK(TEMP_AIR(K),FREQ,TEMP(K))
      ENDDO
      CALL PLANCK(2.7_r8,FREQ,TSPACE)
!------------------------------------------------
!     DETERMINE NO. OF ITERATIONS
!------------------------------------------------
      ITS0=1
      DELTA=0.01_r8
      jj0=1.e-9_r8

      DO K=1,L
         CHK(K)=0._r8
         DO ISPI=1,N
            CHK(K)=CHK(K)+W0(ISPI,K)
         ENDDO
         IF(CHK(K).NE.0.) ITS0=20
      ENDDO
!-----------------------------------------------------
!     CALCULATE EQUIVALENT U FOR EACH TANGENT HEIGHT
!----------------------------------------------------
      Tau=0._r8
      DO I=1,NT
         LMIN(I)=1
         if (zt(i) .ge. 0._r8) then
         DO K=1,L
            IF(YZ(K) .LE. ZT(I)) LMIN(I)=K
         ENDDO
         endif

         DO K=LMIN(I),L
            TGT= YZ(K)
            if(tgt .lt. zt(i)) tgt=zt(i)
            IF(ZT(I) .LT. 0._r8) THEN
               D2=SQRT((YZ(K+1)+RE)**2-(ZT(I)+RE)**2)
               D1=SQRT((TGT+RE)**2-(ZT(I)+RE)**2)
               UAVE(I,K)=(YZ(K+1)-YZ(K))/(D2-D1)      
            ELSE  
               DY=SQRT(YZ(K+1)-ZT(I))-SQRT(TGT-ZT(I))
               UAVE(I,K)=(YZ(K+1)-YZ(K))/SQRT(2*RE+TGT+ZT(I))/DY
            ENDIF

         ENDDO
!---------------------------------------------------
!     Calculate integrated optical depth
!---------------------------------------------------
!         Do K=L-1,LMIN(I)-1,-1
!         K1=L+1-K
!         Tau(I,K1)=Tau(I,K1-1)+dTau(K)/UAVE(I,K)
!         Enddo
!         Tau(I,L+LMIN(I))=Tau(I,L-LMIN(I)+1)
!         Do K=LMIN(I),L
!         Tau(I,K+L+1)=Tau(I,K+L)+dTau(K)/UAVE(I,K)
!         Enddo
!--------------------------------------------------
      ENDDO    ! end of tangent height

!------------------------------------------------
!     INITIAL GUESS OF TB, L=0 IS THE SURFACE
!------------------------------------------------
      DO K=1,L+1
         XTB(K)=0._r8
         DO I=1,NU/2
            TB(I,K)=FIRST
         ENDDO
         DO I=NU/2+1,NU
            TB(I,K)=TSPACE
         ENDDO
      ENDDO

!------------------------------------------------
!     START ITERATION
!------------------------------------------------
      ITS=1
 1000 CONTINUE

      DO J=1,NU/2
         TB(J+NU/2,L+1)=TSPACE         ! BOUNDARY CONDITION AT TOP OF ATMOS
      ENDDO

      DO J=1,NT
         TT(J,L+1)=TSPACE              ! BOUNDARY CONDITION AT TOP OF ATMOS
         TT(J,L)=TSPACE                ! INITIALIZATION
      ENDDO

!------------------------------------------------
!     COMPUTING TSCAT FOR EACH ANGLE AND LAYER
!------------------------------------------------
     
      IF(ICON .gt. 0) then
      DO ISPI=1,N                         
         DO 1009 IP=1,NU
            U1(IP)=-U(IP)
            DO 1008 K=1,L
               TSCAT(ISPI,IP,K)=0.0_r8
               IF(PHH(ISPI,1,K).NE.0._r8)THEN 
                  DO IH=1,NU
                     TAVG=0._r8
                     DO J=1,NUA
                        K1=K
                        IF(THETAI(IP,IH,J).LT.PI/2) K1=K-1
                          IF(K1 .EQ. 0) K1 = 1
                        IF(THETAI(IP,IH,J).GT.PI/2) K1=K+1
                        CALL LOCATE(THETA,NU,NU,THETAI(IP,IH,J),JM)
                        WK=(TB(JM+1,K1)*(THETAI(IP,IH,J)-THETA(JM))+  &
      &                     TB(JM,K1)*(THETA(JM+1)-THETAI(IP,IH,J)))/ &
      &                     (THETA(JM+1)-THETA(JM))
                        TAVG = TAVG + WK*(PHI(2)-PHI(1))/2._r8/PI
                     ENDDO
                     TSCAT(ISPI,IP,K)=TSCAT(ISPI,IP,K)+               &
      &                                2._r8*PHH(ISPI,IH,K)*TAVG*DU(IH)
                  ENDDO
               ENDIF
 1008       CONTINUE
 1009    CONTINUE
      ENDDO
      ENDIF

!------------------------------------
!     COMPUTE TB AT ANGLES U
!     INTEGRATION FROM TOP TO BOTTOM
!------------------------------------
      DO 1100 I=NU/2+1,NU
         X2 = U(I)*U(I)
         UEFF = ABS(U(I))
         DO 1100 K=L,1,-1
            WK=0._r8
            WW0=0._r8
            DO ISPI=1,N
               WK=WK+TSCAT(ISPI,I,K)*W0(ISPI,K)
               WW0=WW0+W0(ISPI,K)
            ENDDO

            tsource=(1-WW0)*TEMP(K)+WK

            TB(I,K)=TB(I,K+1)*EXP(-dTAU(K)/UEFF)+           &
      &             (1._r8-EXP(-dTAU(K)/UEFF))*tsource

 1100 CONTINUE

!-----------------------------------------------------------------------
!     DETERMINE SURFACE REFLECTION 
!-----------------------------------------------------------------------
      DO J=1,NU/2
         TB0(J)=TS*(1._r8-RS(J))+TB(J+NU/2,1)*RS(J)
         TB(J,1) = TB0(J)
      ENDDO
    
!------------------------------------
!     COMPUTE TB AT ANGLES U
!     INTEGRATION FROM BOTTOM TO TOP
!------------------------------------
      DO 1200 I=1,NU/2
        X2 = U(I)*U(I)
        UEFF = ABS(U(I))
        DO 1200 K=1,L
           WK=0._r8
           WW0=0._r8
           DO ISPI=1,N
              WK=WK+TSCAT(ISPI,I,K)*W0(ISPI,K)
              WW0=WW0+W0(ISPI,K)
           ENDDO

           tsource=(1-WW0)*TEMP(K)+WK

           TB(I,K+1)=TB(I,K)*EXP(-dTAU(K)/UEFF)+      &
      &              (1._r8-EXP(-dTAU(K)/UEFF))*tsource

 1200 CONTINUE

!------------------------------------
!     CHECK CONVERGENCE
!------------------------------------
      IND=0
      X2=0._r8
      DO J=1,NU/2
         X2=X2+(XTB(J)-TB(J,L))**2/NU
         IF (ABS(XTB(J)-TB(J,L)) .GT. DELTA) IND=1
      ENDDO
 
      IF (IND.NE.0 .AND. ITS.LT.ITS0) THEN
      DO J=1,NU
         XTB(J)=TB(J,L)
      ENDDO
         IF (ITS .GT. 1) THEN
!            WRITE(*,*)'ITERATION=',ITS-1,'  EPS. AT TOP OF ATMOS=',X2
         ELSE
!            WRITE(*,*)'ITERATION=',ITS-1 
         ENDIF
         ITS=ITS+1
         GOTO 1000
      ENDIF

!-----------------------------------------------------------
!     AFTER TB CONVERGENCE IS FOUND, DO TT CALCULATIONS
!     BY PROJECTING TB AT EACH TANGENT HEIGHT ZT
!     INTEGRATION FROM TOP TO BOTTOM
!-----------------------------------------------------------
      DO 2000 ITT=1,NT
      DO 2000 K=1,L-LMIN(ITT)
         K1=L-K
         UU=UAVE(ITT,K1)                 ! INCIDENT ANGLE AT ZT(ITT) 
         CALL LOCATE(U1,NU/2,NU,-UU,JM)  ! INTERPOLATE TSCAT ONTO UU
         WK=0._r8
         WW0=0._r8
         DO ISPI=1,N
            WK1=W0(ISPI,K1)*( (TSCAT(ISPI,JM+NU/2+1,K1)*(UU-U(JM))+  &
      &         TSCAT(ISPI,JM+NU/2,K1)*(U(JM+1)-UU))/                &
      &         (U(JM+1)-U(JM)) )
            WK=WK+WK1
            WW0=WW0+W0(ISPI,K1)
         END DO

         D1= TEMP(K1)
! The second term takes into account temperature slope  (causing error < 0.1K)
         D2= (TEMP(K1+1)-TEMP(K1))*dTau(K1)/UU/2
         tsource=(1._r8-WW0)*D1+WK 

         IF (ICON .eq. 3) THEN
            tsource=D1 ! NO CLOUD AFTER TANGENT POINT
         ENDIF

         TT(ITT,K1)=TT(ITT,K1+1)*EXP(-dTAU(K1)/UU)+  &
      &        (1._r8-EXP(-dTAU(K1)/UU))*tsource - D2*(1._r8-WW0)

 2000 CONTINUE
 
!---------------------------------------------
!     DETERMINE SURFACE REFLECTION 
!     IF TANGENT HEIGHT IS ABOVE THE SURFACE, 
!     TT(ITT,LMIN(ITT)) IS PASSED ALONG
!----------------------------------------------

      DO 2500 ITT=1,NT
         IF(LMIN(ITT).EQ.1) THEN
            CALL LOCATE(U1,NU/2,NU,-UAVE(ITT,1),JM) 
            RSAVG=(RS(JM+1)*(UAVE(ITT,1)-U(JM)) +                 &
      &            RS(JM)*(U(JM+1)-UAVE(ITT,1)))/(U(JM+1)-U(JM))
            TT(ITT,1)=(TB0(JM+1)*(UAVE(ITT,1)-U(JM)) +            &
      &               TB0(JM)*(U(JM+1)-UAVE(ITT,1)))/(U(JM+1)-U(JM))
         ENDIF
 2500 CONTINUE

!-----------------------------------------
!     NOW INTEGRATION TT FROM BOTTOM UP
!-----------------------------------------
      DO 3000 ITT=1,NT
      DO 3000 K=LMIN(ITT),L
         UU=UAVE(ITT,K)
         CALL LOCATE(U1,NU/2,NU,UU,JM)      ! INTERPOLATE TSCAT ONTO UU
         WK=0._r8
         WW0=0._r8
         DO ISPI=1,N
            WK1=W0(ISPI,K)*(TSCAT(ISPI,JM+1,K)*(UU-U(JM))+     &
      &         TSCAT(ISPI,JM,K)*(U(JM+1)-UU))/(U(JM+1)-U(JM)) 
            WK=WK+WK1
            WW0=WW0+W0(ISPI,K)
         END DO

         D1= TEMP(K)
! The second term takes into account temperature slope (causing error < 0.1K)
         D2= (TEMP(K+1)-TEMP(K))*dTau(K1)/UU/2
         tsource=(1._r8-WW0)*D1+WK 
         TT(ITT,K+1)=TT(ITT,K)*EXP(-dTAU(K)/UU)+  &
      &        (1._r8-EXP(-dTAU(K)/UU))*tsource - D2*(1._r8-WW0)

 3000 CONTINUE

!---------------------------------------------------------------------
      DO K=1,L
   	TT(NT+1,K)=TB(NU,K)              ! OUTPUT ZENITH LOOKING TB 
      ENDDO

      do i=1,nt+1
        tt(i,L+1)=tt(i,L)
      enddo
!---------------------------------------------------------------------
 
      END SUBROUTINE RADXFER

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module RadiativeTransferModule

! $Log$
! Revision 1.15  2004/01/07 19:18:31  jonathan
! initialize U1 with 0
!
! Revision 1.14  2003/05/07 23:11:32  jonathan
! some clean-up and cosmetic changes
!
! Revision 1.13  2003/04/24 23:21:03  jonathan
! currect a comment typo
!
! Revision 1.12  2003/04/10 20:25:15  dwu
! make i_saturation as a verbal argument
!
! Revision 1.11  2003/02/03 19:24:46  dwu
! make it more efficient
!
! Revision 1.10  2002/12/06 16:31:06  dwu
! Add a correction term to account for integration in the presence of temperature gradient
!
! Revision 1.9  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.8  2002/08/19 22:21:53  jonathan
! debug stuff
!
! Revision 1.7  2002/08/08 22:47:06  jonathan
! newly improved version
!
! Revision 1.6  2001/10/25 16:45:15  dwu
! fix problem of  dimension index=0 in tscat calculation
!
! Revision 1.5  2001/10/05 16:03:02  dwu
! *** empty log message ***
!
! Revision 1.4  2001/09/21 15:51:37  jonathan
! modified F95 version
!
