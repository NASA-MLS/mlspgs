
      SUBROUTINE MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMR,NZ,NS,  &
                 &           N, WCin, IPSDin,                        &
                 &           YP,YZ,YT,YQ,VMR1,WC,NH,CHK_CLD,IPSD,    &
                 &           ZT,ZZT,NT)

!========================================================================C
!     DESCRIPTION                                                        C
!     -----------                                                        C
!     1: CHECK THE INPUT L2 ATMOSPHERIC PROFILES.                        C
!     2: OUTPUT INTERPOLATED ATMOSPHERIC PROFILES IF NEEDED TO ENSURE    C
!        THERE ARE 640 MODEL LAYERS FROM SURFACE TO ~40KM.               C
!     3: CONVERT TANGENT PRESSURE TO TANGENT HEIGHT                      C
!                                                                        C
!     J. JIANG, MAY 18, 2001                                             C
!========================================================================C

     use MLSCommon, only: r8
!----------------------------------------------
!     INPUT PARAMETERS
!----------------------------------------------

      INTEGER :: NZ                            ! NO. OF L2 ATMOSPHERIC LEVELS
      REAL(r8) :: PRESSURE(NZ)                 ! PRESSURE LEVEL
      REAL(r8) :: HEIGHT(NZ)                   ! PRESSURE HEIGHT
      REAL(r8) :: TEMPERATURE(NZ)              ! ATMOSPHERIC TEMPERATURE
      REAL(r8) :: VMR(NS,NZ)                   ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3

      REAL(r8) :: WCin(N,NZ)                         
      INTEGER :: IPSDin(NZ)
      INTEGER :: NT                            ! NO. OF TANGENT PRESSURE LEVSLS
      REAL(r8) :: ZT(NT)                       ! TANGENT PRESSURE
      
!----------------------------------------------
!     OUTPUT PARAMETERS
!----------------------------------------------
      INTEGER :: NH                            ! MODEL ATMOSPHERIC LEVELS

      REAL(r8) :: YZ(NH)                       ! PRESSURE HEIGHT (m)
      REAL(r8) :: YP(NH)                       ! PRESSURE (hPa)
      REAL(r8) :: YT(NH)                       ! TEMPERATURE PROFILE
      REAL(r8) :: YQ(NH)                       ! RELATIVE HUMIDITY (%)
      REAL(r8) :: VMR1(NS,NH)                  ! 1=O3 VOLUME MIXING RATIO
      REAL(r8) :: WC(N,NH)
      INTEGER :: IPSD(NH)
      REAL(r8) :: CHK_CLD(NH)                  ! CLOUD CHECKER      

      REAL(r8) :: ZZT(NT)                      ! TANGENT HEIGHT

!----------------------------------------------------------
!     WORK SPACE
!----------------------------------------------------------
      REAL(r8) :: PTOP,PBOTTOM,DP,ZH(NH),ZA(NH),ZZ(NH),WK, zvmr(nh)
      INTEGER :: I,JM,J
!--------------------------------------------------------------------------

      PTOP = (80./16.-3.)*1._r8    ! TOP OF THE MODEL
      PBOTTOM=-LOG10(PRESSURE(1))  ! BOTTOM OF THE MODEL
      DP=(PTOP-PBOTTOM)/NH         ! LAYER THICKNESS

      DO I=1,NH
         ZH(I)=PBOTTOM+(I-1)*DP
      END DO

      DO I=1,NZ
         ZA(I)=-LOG10( PRESSURE(I) )
      END DO

      DO I=1,NZ
         ZVMR(I)=LOG10( VMR(1,I) )
      END DO


      IF (NZ .NE. NH) THEN

!==========================================
!     PRODUCE MODEL ATMOSPHERIC PROFILES
!==========================================


         DO J=1,NH

            CALL LOCATE (ZA,NZ,NH,ZH(J),JM)              

            if (jm .le. 0) then
              print*,'wronh JM: ',JM
              stop
                jm=1                               ! since there is no pressure(0)
            endif

            YP(J)=((ZA(JM+1)-ZH(J))*ZA(JM)+(ZH(J)-ZA(JM))*  &
     &            ZA(JM+1))/(ZA(JM+1)-ZA(JM))             

            YP(J) = 10**(-YP(J))

            YZ(J)=((ZA(JM+1)-ZH(J))*HEIGHT(JM)+(ZH(J)-ZA(JM))*    &
     &            HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))

            YT(J)=((ZA(JM+1)-ZH(J))*TEMPERATURE(JM)+(ZH(J)-ZA(JM))* &
     &            TEMPERATURE(JM+1))/(ZA(JM+1)-ZA(JM))

            WC(1,J)=((ZA(JM+1)-ZH(J))*WCin(1,JM)+(ZH(J)-ZA(JM))*    &
     &            WCin(1,JM+1))/(ZA(JM+1)-ZA(JM))

            WC(2,J)=((ZA(JM+1)-ZH(J))*WCin(2,JM)+(ZH(J)-ZA(JM))*    &
     &            WCin(2,JM+1))/(ZA(JM+1)-ZA(JM))
            
            YQ(J)=((ZA(JM+1)-ZH(J))*ZVMR(JM)+(ZH(J)-ZA(JM))*       &
     &            ZVMR(JM+1))/(ZA(JM+1)-ZA(JM))

            YQ(J) =10**YQ(J)

            VMR1(1,J)=((ZA(JM+1)-ZH(J))*VMR(2,JM)+(ZH(J)-ZA(JM))*   &
     &                VMR(2,JM+1))/(ZA(JM+1)-ZA(JM))
         
            IPSD(J)=((ZA(JM+1)-ZH(J))*IPSDin(JM)+(ZH(J)-ZA(JM))*    &
     &            IPSDin(JM+1))/(ZA(JM+1)-ZA(JM))

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

!----------------------------------------------------------------------

!==========================================
!     PRODUCE TANGENT HEIGHTS (KM)
!==========================================

      DO I=1,NT
!        print*, ZT(I)
         ZZ(I)=-LOG10( max(1.e-9_r8, ZT(I)) )
      END DO

!      stop

      DO J=1,NT      
         CALL LOCATE (ZA,NZ,NH,ZZ(J),JM)
         ZZT(J)=((ZA(JM+1)-ZZ(J))*HEIGHT(JM)+(ZZ(J)-ZA(JM))*    &
     &                HEIGHT(JM+1))/(ZA(JM+1)-ZA(JM))             
      ENDDO

!----------------------------------------------------------------------

      RETURN
      END

! $Log: Model_atmos.f90,v      






