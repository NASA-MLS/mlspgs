! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ScatSourceFunc

! -------------------------------------------------------------------------  
! PERFORM ITERATIVE RADIATIVE TRANSFER CALCULATION
! -------------------------------------------------------------------------

  implicit NONE
  private
  public :: T_Scat

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

   subroutine T_SCAT ( TEMP_AIR, FREQ, Z, NU, NUA, NAB, NR, NC, TSCAT )  

      use Cloud_extinction, only: Get_beta_cloud
      use CRREXP_m,         only: RREXP    ! ( exp(x)-1 ) / x, for Planck fn.
      use Interpack,        only: LOCATE
      use MLSCommon,        only: RK => R8 ! working REAL kind
      use Physics,          only: H_OVER_K ! h/k in Kelvin/MHz
      use ScatteringAngle,  only: Angle
      use Units,            only: Pi

    ! Arguments
      real(rk), intent(in) :: Temp_Air(:) ! MEAN AIR TEMPERATURES, Kelvin
      real(rk), intent(in) :: Freq        ! FREQUENCY, MHz
      real(rk), intent(in) :: Z(:)        ! Model height (meters)
      integer, intent(in) :: NU           ! Number of scattering angles
      integer, intent(in) :: NUA          ! Number of azimuth angles
      integer, intent(in) :: NAB          ! Number of AB terms
      integer, intent(in) :: NR           ! Number of size bins
      integer, intent(in) :: NC           ! Number of cloud species
      real(rk), intent(out) :: TScat(:,:) ! TB FROM SCATTERING PHASE FUNCTION

    ! Local variables
      real(rk), parameter :: COLD = 2.7_rk   ! kelvins
      real(rk), parameter :: FIRST = 250._rk ! FIRST GUESS OF UPWELLING TB
      real(rk), parameter :: Pi2 = 0.5 * Pi

      real(rk) :: CHK(size(TEMP_AIR))
      real(rk) :: CLD_EXT
      real(rk) :: D1
      real(rk) :: D2
      real(rk) :: DELTA                   ! DELTA TB FOR CONVERGENCE CHECK
      real(rk) :: dTAU( size(temp_air) )  ! Optical depth increment at each layer 
      real(rk) :: DU(NU)                  ! DELTA U
      real(rk) :: DY
      real(rk) :: ITS0                    ! NO. OF MAXIMUM ITERATIONS
      real(rk) :: JJ0
      real(rk) :: PHH( NU, size(TEMP_AIR) )
      real(rk) :: PHI(NUA)                ! SCATTERING AZIMUTH ANGLES
      real(rk) :: RSAVG
      real(rk) :: RS( NU/2 )
      real(rk) :: TAVG                    ! TB AVERAGED OVER PHI AT A GIVEN U
      real(rk) :: TB( NU, size(TEMP_AIR) ) ! TB IN FLAT PLANE GEOMETRY
                                          ! 1->NU/2 UPWELLING, NU/2->NU DOWNWELLING
      real(rk) :: TB0 ( NU )   ! TB AT THE SURFACE
      real(rk) :: TEMP( size(temp_air) )  ! BRIGHTNESS TEMPERATURE FROM TEMP_AIR
      real(rk) :: TGT
      real(rk) :: THETAI(NU,NU,NUA) ! ANGLES FOR INCIDENT TB
      real(rk) :: THETA(NU)    ! SCATTERING ANGLES
      real(rk) :: Tsource
      real(rk) :: TSPACE                  ! COSMIC BACKGROUND RADIANCE
      real(rk) :: TS                      ! SURFACE TEMPERATURE (K)
      real(rk) :: U1(NU)
      real(rk) :: UA(NUA)                 ! COSINES OF SCATTERING AZIMUTH ANGLES
      real(rk) :: UAVE(size(temp_air),size(temp_air)) ! INCIDENT ANGLES FOR EACH TANGENT HT
      real(rk) :: UEFF                    ! EFFECTIVE U BETWEEN K AND K+1
      real(rk) :: UI(NU,NU,NUA) ! COSINES OF THETAI
      real(rk) :: U(NU)        ! COSINES OF THETA
      real(rk) :: UU
      real(rk) :: W0( size(temp_air) )
      real(rk) :: WC(2, size(temp_air) )
      real(rk) :: WK
      real(rk) :: WK1
      real(rk) :: WW0
      real(rk) :: Wwk
      real(rk) :: Wwk1
      real(rk) :: Www0
      real(rk) :: X2
      real(rk) :: XTB(size(TEMP_AIR))

      integer :: I
      integer :: ICON                     ! CONTROL SWITCH
                                          ! 3 = NEAR SIDE CLOUD ONLY
      integer :: IH, IND, IP, ISPI, ITS, ITT, J, JM, K, K1, L
      integer :: LMIN( size(temp_air) )   ! LOWEST LAYER REACHED BY A TANGENT HT

!-------------------------------------------------------------------------------------

      L = size(temp_air,1) ! == size(tscat,1)
      TB    = 0.0_rk
      Tscat = 0.0_rk
      tsource = 0.0_rk

!--------------------------------------------------
!     FIND BRIGHTNESS TEMPERATURE AT EACH LAYER
!--------------------------------------------------

      !{Planck black body function:
      ! $\frac{\frac{h \nu}k}{\exp\left(\frac{h \nu}{k T}\right)-1}$.
      ! We need to be careful how this is calculated because $\exp x -1$ has
      ! cancellation for small $x$.  $\frac{h \nu}{k T}$ is $\approx 0.04$ for
      ! $\nu \approx 240$ GHz and $T \approx 250$ K.

      temp = temp_air / rrexp(h_over_k * freq / temp_air)
      tspace = cold / rrexp(h_over_k * freq / cold)

      call ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)
      
      WC=0.0
      PHH=0.0
      cld_ext=0.0
      WC(1,10) = 0.01   !test only
      WC(1,11) = 0.01   !test only
      dtau =0.0

!         CALL CLEAR_SKY(L-1,NU,TS,S,LORS,SWIND,                         &
!              &         YZ,YP,YT,YQ,VMR,NS,                             &
!              &         FREQUENCY(IFR),RS,U,TEMP,Z,TAU0,tau_wetAll,     &
!              &         tau_dry,Catalog, Bill_data, LosVel, i_saturation ) 

      do K=1,L-1
        call get_beta_cloud ( FREQ, TEMP_AIR(K), WC(:,K), 1000,        &
     &                 NC, NU, NUA, NAB, NR, cld_ext, W0(K), PHH(:,K) )      
!        print*, FREQ, TEMP_AIR(K), WC(1,K), cld_ext, W0(K)
        dtau(k)=(Z(K+1)-Z(K))*cld_ext
      end do

!      STOP

!----------------------------------------------------
!     DETERMINE NO. OF ITERATIONS
!----------------------------------------------------

      ITS0=1
      DELTA=0.01_rk
      jj0=1.e-9_rk

      if ( any(w0 /= 0.0) ) ITS0=20
      
!---------------------------------------------------
!     INITIAL GUESS OF TB, L=0 IS THE SURFACE
!---------------------------------------------------
      do K=1,L
         XTB(K)=0._rk
         TB(1:NU/2,K)=FIRST
         TB(NU/2+1:NU,K)=TSPACE
      end do

!------------------------------------------------
!     START ITERATION
!------------------------------------------------
      ITS=1
 1000 CONTINUE
    
      TB(NU/2+1:NU,L)=TSPACE         ! BOUNDARY CONDITION AT TOP OF ATMOS

!------------------------------------------------
!     COMPUTING TSCAT FOR EACH ANGLE AND LAYER
!------------------------------------------------

      do IP=1,NU
         U1(IP)=-U(IP)
         do K=1,L-1
            TSCAT(IP,K)=0.0_rk
            if ( PHH(1,K) /= 0._rk ) then 
               do IH=1,NU
                  TAVG=0.0_rk
                  do J=1,NUA
                     K1=K
                     if ( THETAI(IP,IH,J) < pi2) K1=K-1
                     if ( K1 .EQ. 0) K1 = 1
                     if ( THETAI(IP,IH,J) > pi2) K1=K+1
                     call LOCATE(THETA,NU,NU,THETAI(IP,IH,J),JM)
                     WK=(TB(JM+1,K1)*(THETAI(IP,IH,J)-THETA(JM))+      &
     &                   TB(JM,K1)*(THETA(JM+1)-THETAI(IP,IH,J)))/     &
     &                   (THETA(JM+1)-THETA(JM))
                     TAVG = TAVG + WK*(PHI(2)-PHI(1))/2._rk/PI
                  end do
                  TSCAT(IP,K)=TSCAT(IP,K)+                             &
     &                              2._rk*PHH(IH,K)*TAVG*DU(IH)
               end do
            else
              TSCAT(IP,K)=0._rk         !CLEAR SKY
            end if
         end do
      end do

!------------------------------------
!     COMPUTE TB AT ANGLES U
!     INTEGRATION FROM TOP TO BOTTOM
!------------------------------------
      DO 1100 I=NU/2+1,NU
         X2 = U(I)*U(I)
         UEFF = ABS(U(I))
         DO 1100 K=L-1,1,-1
            WK=0._rk
            WW0=0._rk
               WK=WK+TSCAT(I,K)*W0(K)
               WW0=WW0+W0(K)

            tsource=(1-WW0)*TEMP(K)+WK

            TB(I,K)=TB(I,K+1)*EXP(-dTAU(K)/UEFF)+           &
     &              (1._rk-EXP(-dTAU(K)/UEFF))*tsource

 1100 CONTINUE

!-----------------------------------------------------------------------
!     DETERMINE SURFACE REFLECTION 
!-----------------------------------------------------------------------
      TS    = 288._rk

      DO J=1,NU/2
         RS(J)=0.05_rk
         TB0(J)=TS*(1._rk-RS(J))+TB(J+NU/2,1)*RS(J)
         TB(J,1) = TB0(J)
      ENDDO

!------------------------------------
!     COMPUTE TB AT ANGLES U
!     INTEGRATION FROM BOTTOM TO TOP
!------------------------------------
      DO 1200 I=1,NU/2
        X2 = U(I)*U(I)
        UEFF = ABS(U(I))
        DO 1200 K=1,L-1
          WK=0._rk
           WW0=0._rk
              WK=WK+TSCAT(I,K)*W0(K)
              WW0=WW0+W0(K)

           tsource=(1-WW0)*TEMP(K)+WK

           TB(I,K+1)=TB(I,K)*EXP(-dTAU(K)/UEFF)+      &
     &               (1._rk-EXP(-dTAU(K)/UEFF))*tsource

 1200 CONTINUE

!------------------------------------
!     CHECK CONVERGENCE
!------------------------------------
      IND=0
      X2=0._rk
      DO J=1,NU/2
         X2=X2+(XTB(J)-TB(J,L))**2/NU
         IF (ABS(XTB(J)-TB(J,L)) .GT. DELTA) IND=1
      ENDDO
 
      IF (IND.NE.0 .AND. ITS.LT.ITS0) THEN
      DO J=1,NU
         XTB(J)=TB(J,L)
      ENDDO

         ITS=ITS+1
         GOTO 1000
      ENDIF

   end subroutine T_SCAT

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ScatSourceFunc

! $Log$
! Revision 2.2  2003/10/28 22:06:37  jonathan
! add Z as input variable for later use
!
! Revision 2.1  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.8  2003/04/24 23:54:18  jonathan
! update construction
!
! Revision 1.1.2.7  2003/04/16 19:31:43  vsnyder
! Move working storage into T_Scat, using automatic arrays.  Use assumed-
! shape dummy arguments.  Alphabetize declarations of local variables.
! Use RREXP to calculate Planck blackbody function, to avoid cancellation
! for small values of h nu / k T.
!
! Revision 1.1.2.6  2003/04/16 17:26:09  vsnyder
! Planck now expects frequency in MHz
!
! Revision 1.1.2.5  2003/04/16 00:29:17  vsnyder
! Move USE statements down to procedure scope, futzing
!
! Revision 1.1.2.4  2003/04/10 16:13:45  jonathan
! fix bug
!
! Revision 1.1.2.3  2003/04/07 17:14:00  jonathan
! modified cal to T_scat
!
! Revision 1.1.2.2  2003/04/04 20:24:30  jonathan
! first working version
!
! Revision 1.1.2.1  2003/03/28 18:14:53  jonathan
! new scaterring source module in construction
!
