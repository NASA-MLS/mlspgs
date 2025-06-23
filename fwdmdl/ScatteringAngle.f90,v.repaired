! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ScatteringAngle

      use MLSCommon, only: r8
      Private
      Public :: ANGLE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
      
contains

      SUBROUTINE ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

!==========================================================
!     SETUP SCATTERING ANGLES U, UA AND INCIDENT ANGLES UI
!     -J.JIANG, JAN 1, 2001
!==========================================================

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

      END SUBROUTINE ANGLE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ScatteringAngle

! $Log$
! Revision 2.3  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/01/31 19:35:08  jonathan
! moved from cloudfwm
!
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
