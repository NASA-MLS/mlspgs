! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
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
                 &           NcloudType, WCin, IPSDin,               &
                 &           YP,YZ,YT,YQ,VMR1,WC,NH,CHK_CLD,IPSD )

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

      INTEGER, intent(in) :: NZ                ! NO. OF L2 ATMOSPHERIC LEVELS
      INTEGER, intent(in) :: NS                ! NO. OF CHEMICAL SPECIES
      INTEGER, intent(in) :: NcloudType
      INTEGER, intent(in) :: NH                ! MODEL ATMOSPHERIC LEVELS
      REAL(r8), intent(in) :: PRESSURE(NZ)     ! PRESSURE LEVEL
      REAL(r8), intent(in) :: HEIGHT(NZ)       ! PRESSURE HEIGHT
      REAL(r8), intent(in) :: TEMPERATURE(NZ)  ! ATMOSPHERIC TEMPERATURE
      REAL(r8), intent(in) :: VMR(NS,NZ)       ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3

      REAL(r8), intent(in) :: WCin(NcloudType,NZ)             
      INTEGER, intent(in) :: IPSDin(NZ)
      
!----------------------------------------------
!     OUTPUT PARAMETERS
!----------------------------------------------

      REAL(r8), intent(out) :: YZ(NH)          ! PRESSURE HEIGHT (m)
      REAL(r8), intent(out) :: YP(NH)          ! PRESSURE (hPa)
      REAL(r8), intent(out) :: YT(NH)          ! TEMPERATURE PROFILE
      REAL(r8), intent(out) :: YQ(NH)          ! RELATIVE HUMIDITY (%)
      REAL(r8), intent(out) :: VMR1(NS-1,NH)   ! 1=O3 VOLUME MIXING RATIO
      REAL(r8), intent(out) :: WC(NcloudType,NH)
      INTEGER, intent(out) :: IPSD(NH)
      REAL(r8), intent(out) :: CHK_CLD(NH)     ! CLOUD CHECKER

!----------------------------------------------------------
!     WORK SPACE
!----------------------------------------------------------
      REAL(r8) :: HTOP,DH,ZH(NH),ZA(NH), zvmr(nh)
      REAL(r8) :: eta                          ! Interpolating fraction;
      INTEGER :: I,JM,J, K                     !  0 < eta < 1
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
         ZVMR(I)=LOG10( max(1.e-19_r8, VMR(1,I)) )
      END DO


!==========================================
!     PRODUCE MODEL ATMOSPHERIC PROFILES
!==========================================

         DO J=1,NH

            CALL LOCATE (HEIGHT,NZ,NH,ZH(J),JM)              

            ! Now we will interpolate things to our new heights
            ! To do this we use the fraction eta : 0 <= eta <= 1
            ! such that f[j] = eta f[j[m]] + (1-eta) f[j[m+1]]
            ! This is simple linear interpolation; see also 
            ! MLSNumerics for InterpolateValues
            eta = (HEIGHT(JM+1)-ZH(J)) / (HEIGHT(JM+1)-HEIGHT(JM))

            YP(J) = eta*ZA(JM) + (1-eta)*ZA(JM+1)

            YP(J) = 10**(-YP(J))

            YZ=ZH

            YT(J) = eta*TEMPERATURE(JM) + (1-eta)*TEMPERATURE(JM+1)

! ICE QUANTITIES
            DO K=1,NcloudType
              WC(K,J) = eta*WCin(K,JM) + (1-eta)*WCin(K,JM+1)
            ENDDO
            
              CHK_CLD(J) = Sum(WC(1:NcloudType,J))

            IPSD(J) = IPSDin(1)     

! VMR QUANTITIES
            YQ(J) = eta*zvmr(JM) + (1-eta)*zvmr(JM+1)
            YQ(J) = 10**YQ(J)

            DO K=1,NS-1
             VMR1(K,J) = eta*VMR(K+1,JM) + (1-eta)*VMR(K+1,JM+1)
            ENDDO
       
         ENDDO

      END SUBROUTINE MODEL_ATMOS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ModelInput

! $Log$
! Revision 1.13  2003/05/07 23:18:40  jonathan
! some clean-ups
!
! Revision 1.12  2003/04/09 08:16:53  dwu
! Fix the bug associated with number of cloud types dimension
!
! Revision 1.11  2003/01/23 00:19:09  pwagner
! Some cosmetic only (or so I hope) changes
!
! Revision 1.10  2002/12/18 16:10:21  jonathan
! minor changes
!
! Revision 1.9  2002/12/02 17:44:15  dwu
! remove the bug in interpolating IPSD
!
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
