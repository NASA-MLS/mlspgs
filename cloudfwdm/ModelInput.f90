! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ModelInput

! -------------------------------------------------------------------------  
! SET ALL PARAMETERS ONTO INTERNAL MODEL GRIDS
! -------------------------------------------------------------------------

      use Interpack, only: LOCATE
      use MLSCommon, only: r8

      IMPLICIT NONE
      Private
      Public :: MODEL_ATMOS 

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMR,NZ,NS,  &
                 &           N, WCin, IPSDin,                        &
                 &           YP,YZ,YT,YQ,VMR1,WC,NH,CHK_CLD,IPSD,    &
                 &           ZT,ZZT,NT)

!========================================================================C
!     DESCRIPTION                                                        C
!     -----------                                                        C
!     1: CHECK THE INPUT L2 ATMOSPHERIC PROFILES.                        C
!     2: OUTPUT INTERPOLATED ATMOSPHERIC PROFILES IF NEEDED TO INTERNAL  C
!        MODEL GRIDS                                                     C
!     3: CONVERT TANGENT PRESSURE TO TANGENT HEIGHT                      C
!                                                                        C
!     J. JIANG, MAY 18, 2001                                             C
!========================================================================C

!----------------------------------------------
!     INPUT PARAMETERS
!----------------------------------------------

      INTEGER :: NZ                            ! NO. OF L2 ATMOSPHERIC LEVELS
      INTEGER :: NS                            ! NO. OF CHEMICAL SPECIES
      INTEGER :: N
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
      REAL(r8) :: VMR1(NS-1,NH)                  ! 1=O3 VOLUME MIXING RATIO
      REAL(r8) :: WC(N,NH)
      INTEGER :: IPSD(NH)
      REAL(r8) :: CHK_CLD(NH)                  ! CLOUD CHECKER      

      REAL(r8) :: ZZT(NT)                      ! TANGENT HEIGHT

!----------------------------------------------------------
!     WORK SPACE
!----------------------------------------------------------
      REAL(r8) :: HTOP,DH,ZH(NH),ZA(NH),ZZ(NH),WK, zvmr(nh)
      INTEGER :: I,JM,J, K
!--------------------------------------------------------------------------

      HTOP = 80.e3_r8    ! TOP OF THE MODEL
      DH=HTOP/NH         ! LAYER THICKNESS

      DO I=1,NH
         ZH(I)=(I-1)*DH
      END DO

      DO I=1,NZ
         ZA(I)=-LOG10( PRESSURE(I) )
      END DO

      DO I=1,NZ
         ZVMR(I)=LOG10( max(1.e-39_r8, VMR(1,I)) )
      END DO


!==========================================
!     PRODUCE MODEL ATMOSPHERIC PROFILES
!==========================================

         DO J=1,NH

            CALL LOCATE (HEIGHT,NZ,NH,ZH(J),JM)              

            YP(J)=((HEIGHT(JM+1)-ZH(J))*ZA(JM)+(ZH(J)-HEIGHT(JM))*  &
     &            ZA(JM+1))/(HEIGHT(JM+1)-HEIGHT(JM))             

            YP(J) = 10**(-YP(J))

            YZ=ZH

            YT(J)=((HEIGHT(JM+1)-ZH(J))*TEMPERATURE(JM)+(ZH(J)-HEIGHT(JM))*  &
     &            TEMPERATURE(JM+1))/(HEIGHT(JM+1)-HEIGHT(JM))             

! ICE QUANTITIES

            DO K=1,2
            WC(K,J)=((HEIGHT(JM+1)-ZH(J))*WCin(K,JM)+(ZH(J)-HEIGHT(JM))*  &
     &            WCin(K,JM+1))/(HEIGHT(JM+1)-HEIGHT(JM))             
            ENDDO

            CHK_CLD(J) = WC(1,J) + WC(2,J)

            ! IPSD INDEX CAN NOT BE INTERPOLATED
            IPSD(J) = IPSDin(1)     

! VMR QUANTITIES
            Yq(J)=((HEIGHT(JM+1)-ZH(J))*zvmr(JM)+(ZH(J)-HEIGHT(JM))*  &
     &            zvmr(JM+1))/(HEIGHT(JM+1)-HEIGHT(JM))             


            YQ(J) =10**YQ(J)
          
            DO K=1,NS-1
            VMR1(K,J)=((HEIGHT(JM+1)-ZH(J))*VMR(K+1,JM)+(ZH(J)-HEIGHT(JM))*  &
     &            VMR(K+1,JM+1))/(HEIGHT(JM+1)-HEIGHT(JM))             
            ENDDO
       
         ENDDO

      END SUBROUTINE MODEL_ATMOS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ModelInput

! $Log$
! Revision 1.7  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.6  2002/08/22 00:15:31  jonathan
! change NS to NS-1 for VMR1
!
! Revision 1.5  2002/08/19 22:22:04  jonathan
! debug stuff
!
! Revision 1.4  2001/10/22 15:42:58  jonathan
! pretect vmr to be non-zero valus
!
! Revision 1.3  2001/10/11 22:10:02  dwu
! remove tangent height calculation
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
