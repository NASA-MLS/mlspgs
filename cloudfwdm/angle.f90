
      SUBROUTINE ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

!==========================================================
!     SETUP SCATTERING ANGLES U, UA AND INCIDENT ANGLES UI
!     -J.JIANG, JAN 1, 2001
!==========================================================

      use MLSCommon, only: r8
      INTEGER :: NU                         ! NO. OF SCATTERING ANGLES
      INTEGER :: NUA                        ! NO. OF SCATTERING AZIMUTH ANGLES
      REAL(r8) :: U(NU)                     ! COSINE OF SCATTERING ANGLES
      REAL(r8) :: DU(NU)                    ! DELTA U
      REAL(r8) :: UA(NUA)                   ! COSINE OF SCATTERING AZIMUTH ANGLES
      REAL(r8) :: THETA(NU)                 ! SCATTERING ANGLES
      REAL(r8) :: PHI(NUA)                  ! SCATTERING AZIMUTH ANGLES
      REAL(r8) :: UI(NU,NU,NUA)             ! COSINE OF ANGLES FOR INCIDENT TB
      REAL(r8) :: THETAI(NU,NU,NUA)         ! ANGLES FOR INCIDENT TB

      INTEGER :: I, J, K
!--------------------------------------------------------------

      PI=3.1415926

!--------------------------------------------------
!     INITIALIZE NU THETA ANGLES BETWEEN 0 TO PI
!--------------------------------------------------

      DO I=1,NU/2
         U(I)=COS(PI*(I-0.5)/NU)                   
         U(NU-I+1)=-U(I)
         DU(I)=SIN(PI*(I-0.5)/NU)*PI/NU
         DU(NU-I+1)=DU(I)
      ENDDO

      DO I=1,NU
         THETA(I)=ACOS(U(I))
      ENDDO

!--------------------------------------------------
!     INITIALIZE NUA PHI ANGLES BETWEEN 0 TO 2*PI
!--------------------------------------------------

      DO J=1,NUA
         PHI(J)=2*PI*(J-0.5)/NUA
      ENDDO

!--------------------------------------------------
!     INITIALIZE NU THETA1 ANGLES FOR INCIDENT TB
!--------------------------------------------------

      DO K=1,NU                              
         DO I=1,NU                            
            DO J=1,NUA                         
               UI(K,I,J)=COS(THETA(K))*COS(THETA(I))   &
            &       +SIN(THETA(K))*SIN(THETA(I))*SIN(PHI(J))
               THETAI(K,I,J)=ACOS(UI(K,I,J))
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

! $Log: angle.f90,v      




