! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ScatSourceFunc

! -------------------------------------------------------------------------  
! PERFORM ITERATIVE RADIATIVE TRANSFER CALCULATION
! -------------------------------------------------------------------------

  implicit NONE
  private
  public :: T_Scat, Interp_Tscat, Convert_grid

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

   subroutine T_SCAT ( TEMP_AIR, FREQ, Z, Pres, VMRin, NS, NU, NUA, NAB, NR, NC, &
                     & TB_SCAT, Scat_alb, Cloud_ext, THETA )  

      use Cloud_extinction,    only: Get_beta_cloud
      use CRREXP_m,            only: RREXP    ! ( exp(x)-1 ) / x, for Planck fn.
      use Interpack,           only: LOCATE
      use MLSCommon,           only: RK => R8 ! working REAL kind
      use Physics,             only: H_OVER_K ! h/k in Kelvin/MHz
      use ScatteringAngle,     only: Angle
      use Units,               only: Pi
      use Non_scat_ext,        only: GET_BETA_CLEAR

    ! Arguments
      real(rk), intent(in) :: Temp_Air(:) ! MEAN AIR TEMPERATURES, Kelvin
      real(rk), intent(in) :: Freq        ! FREQUENCY, MHz
      real(rk), intent(in) :: Z(:)        ! Model height (meters)
      real(rk), intent(in) :: Pres(size(Z))     ! Pressure (hPa)
      integer, intent(in) :: NU           ! Number of scattering angles
      integer, intent(in) :: NUA          ! Number of azimuth angles
      integer, intent(in) :: NAB          ! Number of AB terms
      integer, intent(in) :: NR           ! Number of size bins
      integer, intent(in) :: NC           ! Number of cloud species
      integer, intent(in) :: NS           ! Number of chemical species
      real(rk), intent(in) :: VMRin(NS, size(Z) )        ! VMR
      
      real(rk), intent(inout) :: TB_scat(:,:)  ! TB FROM SCATTERING PHASE FUNCTION
      real(rk), intent(inout) :: Scat_alb(:,:) ! Single Scattering albedo 
      real(rk), intent(inout) :: Cloud_ext(:,:) ! Cloud extinction

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
      real(rk) :: ext_air( size(temp_air) )
      real(rk) :: DU(NU)                  ! DELTA U
      real(rk) :: DY
      real(rk) :: ITS0                    ! NO. OF MAXIMUM ITERATIONS
      real(rk) :: JJ0
      real(rk) :: PHH( NU, size(TEMP_AIR) )
      real(rk) :: PHI(NUA)                ! SCATTERING AZIMUTH ANGLES
      real(rk) :: RSAVG
      real(rk) :: RS( NU/2 )
      real(rk) :: TAVG                    ! TB AVERAGED OVER PHI AT A GIVEN U
      real(rk) :: TB( NU, size(TEMP_AIR) )! TB IN FLAT PLANE GEOMETRY
                                          ! 1->NU/2 UPWELLING, NU/2->NU DOWNWELLING
      real(rk) :: TB0 ( NU )   ! TB AT THE SURFACE
      real(rk) :: TEMP( size(temp_air) )  ! BRIGHTNESS TEMPERATURE FROM TEMP_AIR
      real(rk) :: TGT
      real(rk) :: THETAI(NU,NU,NUA) ! ANGLES FOR INCIDENT TB
      real(rk) , intent(inout) :: THETA(:)    ! SCATTERING ANGLES
      real(rk) :: TSCAT(NU,size(Z))
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

      real(rk) :: mid_Z(size(Z)-1)
      real(rk) :: D_mid_Z(size(Z)-2)

      integer :: I
      integer :: ICON                     ! CONTROL SWITCH
                                          ! 3 = NEAR SIDE CLOUD ONLY
      integer :: IH, IND, IP, ISPI, ITS, ITT, J, JM, K, K1, L
      integer :: LMIN( size(temp_air) )   ! LOWEST LAYER REACHED BY A TANGENT HT

!-------------------------------------------------------------------------------------

      L = size(temp_air,1) ! == size(tb_scat,1)
      TB    = 0.0_rk
      TB_scat = 0.0_rk
      Tscat   = 0.0_rk
      tsource = 0.0_rk
      Scat_alb= 0.0_rk
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
      W0=0.0_rk
      Cloud_ext =0.0_rk
!      WC(1,10) = 0.01   !test only
!      WC(1,11) = 0.01   !test only
      dtau =0.0

      CALL get_beta_clear ( L, FREQ, TEMP_AIR, Pres, VMRin, NU, NS, ext_air )

      do I=1, L-1
        mid_Z(I) = Z(I) +(Z(I+1)-Z(I))/2.
      enddo

      do I=1, L-2
        D_mid_Z(I) = (Z(I+1)-Z(I))
      enddo

      do K=1,L
        call get_beta_cloud ( FREQ, TEMP_AIR(K), WC(:,K), 1000,         &
     &                 NC, NU, NUA, NAB, NR, cld_ext, W0(K), PHH(:,K) )      
        
        if (K .ge. 2 .and. K .le. L-1) then
            dtau(k)= D_mid_Z(K-1) * (cld_ext + ext_air(K))            
        endif
        W0(K) = (W0(K)*cld_ext)/(cld_ext + ext_air(K))
        Cloud_ext(K,1) = cld_ext
      end do
      dtau(1)=dtau(2) 
      dtau(L)=dtau(L-1)
      
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
      
!-------------------------------------------------------------------------
!     AFTER TB CONVERGENCE IS FOUND, OUTPUT THE SCATERING SOURCE FUNCTION
!-------------------------------------------------------------------------

     do K=1,L
         do IP=1,NU
           TB_scat(K, IP) = Tscat(IP, K)
         enddo
         Scat_alb(K,1)=W0(K)
     enddo

  end subroutine T_SCAT

  subroutine Interp_Tscat (TB_SCAT, THETA, PHI_angle, TT_SCAT)

      use Interpack,           only: LOCATE
      use MLSCommon,           only: RK => R8 ! working REAL kind
      use Physics,             only: H_OVER_K ! h/k in Kelvin/MHz
      use ScatteringAngle,     only: Angle
      use Units,               only: Pi

    ! Arguments
      real(rk), intent(in) :: TB_scat(:,:)   ! TB FROM SCATTERING PHASE FUNCTION
      real(rk), intent(in) :: THETA(:)       ! Scattering angles 
      real(rk), intent(in) :: PHI_angle(:)   ! Phi Angles
      real(rk), intent(inout) :: TT_scat(:,:)! 
      real(rk) :: PHI_90(size(PHI_angle))
      real(rk) :: WK, eta
      INTEGER :: I, JM, No_ele, No_ang
!-----------------------------------------------------------------------------

      No_ele= size(PHI_angle)
      No_ang= size(THETA)

      PHI_90 = pi + PHI_angle

      DO I=1, No_ele

          CALL LOCATE ( THETA, No_ang, No_ang, PHI_90(I), JM )
          TT_SCAT(I,1) = ( TB_SCAT(I,JM)   * (THETA(JM+1)-PHI_90(I)) + &
                       &   TB_SCAT(I,JM+1) * (PHI_90(I)-THETA(JM)) ) / &
                       &   (THETA(JM+1)-THETA(JM))
      ENDDO

  end subroutine Interp_Tscat

  subroutine Convert_grid ( salb_path, cext_path, tt_path, path_inds, &
                          & beta_path_cloud, w0_path, tt_path_c )

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP

    real(rp), intent(in) :: Salb_path(:,:) ! single scattering albedo gl grids  
    real(rp), intent(in) :: Cext_path(:,:) ! cloud extinction on gl grids
    real(rp), intent(in) :: tt_path(:,:)   ! scating source func on gl grids

    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    real(rp), intent(out) :: beta_path_cloud(:) ! cloud extinction on coarse grids
    real(rp), intent(out) :: w0_path(:)         ! single scattering albedo coarse grids
    real(rp), intent(out) :: tt_path_c(:)       ! scattering source func coarse grids

    integer(ip) :: j, k, n_path

! ---------------------------------------------------------------------------------

    beta_path_cloud = 0.0_rp
    w0_path         = 0.0_rp
    tt_path_c       = 0.0_rp

    n_path = size(path_inds)

    do j = 1, n_path

       k = path_inds(j)
       beta_path_cloud(j) = beta_path_cloud(j) + Cext_path(k,1)
       w0_path(j)         = w0_path(j)         + Salb_path(k,1)     
       tt_path_c(j)       = tt_path_c(j)       + tt_path(k,1)  
            
    end do

  end subroutine Convert_grid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ScatSourceFunc

! $Log$
! Revision 2.8  2003/12/07 19:46:25  jonathan
! update for use in 2D cloud FWM
!
! Revision 2.7  2003/12/01 17:25:02  jonathan
! add scat alb
!
! Revision 2.6  2003/11/19 22:22:03  jonathan
! some update
!
! Revision 2.5  2003/11/17 18:04:45  jonathan
! correct working version
!
! Revision 2.3  2003/10/29 17:22:35  jonathan
! fix bug
!
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
