! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!===============================================================================
MODULE Scan
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   USE SDPToolkit

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Scan_guess, Scan_start

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Contents:

! Subroutines -- Scan_guess
!                Scan_start
! 

! Remarks:  This module contains parameters and subroutines needed to model the
!           MLS scan programs.

CONTAINS

!------------------------------------------
   SUBROUTINE Scan_guess (asciiUTC, viewECR)
!------------------------------------------

! Brief description of subroutine
! This subroutine calculates the initial guess for the starting view vector
! of a scan.

! Arguments

      CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

      REAL(r8), INTENT(OUT) :: viewECR(3)

! Parameters

      INTEGER, PARAMETER :: numValues = 1

      REAL(r8), PARAMETER :: time_offset = 0.0

! Functions

      INTEGER :: Pgs_csc_scToECI, Pgs_csc_eciToECR

! Variables

      INTEGER :: returnStatus

      REAL(r8) :: angle
      REAL(r8) :: eci(3), sc(3)
      REAL(r8) :: eciV(6), ecrV(6)

      angle = Deg2Rad * 30

      sc(1) = COS(angle)
      sc(2) = 0.0
      sc(3) = SIN(angle)

      returnStatus = Pgs_csc_scToECI (spacecraftid, numValues, asciiUTC, &
           time_offset, sc, eci)

      eciV(1:3) = eci
      eciV(4:6) = 0.0

      returnStatus = Pgs_csc_eciToECR (numValues, asciiUTC, time_offset, &
           eciV, ecrV)

! Truncate velocity portion of ECR view vector

      viewECR = ecrV(1:3)

!---------------------------
   END SUBROUTINE Scan_guess
!---------------------------

!-----------------------------------------------------------------------
   SUBROUTINE Scan_start (initAlt, posECR, asciiUTC, initRay, startAngle)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine takes an initial-guess view vector and iteratively runs the
! toolkit GrazingRay routine to find the actual view vector for a given
! tangent height.  The returned value is the view vector angle in degrees.

! Arguments

      CHARACTER(LEN=27), INTENT(IN) :: asciiUTC

      REAL(r8), INTENT(IN) :: initAlt
      REAL(r8), INTENT(IN) :: initRay(3), posECR(3)

      REAL(r8), INTENT(OUT) :: startAngle

! Parameters

      INTEGER, PARAMETER :: numValues = 1

      REAL(r8), PARAMETER :: time_offset = 0.0

! Functions

      INTEGER :: Pgs_csc_grazingRay
      INTEGER :: Pgs_csc_geoToECR, Pgs_csc_ecrToECI, Pgs_csc_eciToSc

! Variables

      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: i, returnStatus, converge

      REAL(r8) :: latitude, longitude, missAltitude, slantRange
      REAL(r8) :: ecr(3), eci(3), finalRay(3), posNear(3), posSurf(3), ray(3)
      REAL(r8) :: ecrV(6), eciV(6)

! Initial guess tangent height

      returnStatus = Pgs_csc_grazingRay (earthModel, posECR, initRay, &
           latitude, longitude, missAltitude, slantRange, posNear, posSurf)

!!$      IF (returnStatus /= PGS_S_SUCCESS) THEN
!!$         CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
!!$         msr = mnemonic // ':  ' // msg
!!$         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
!!$      ENDIF

! Iterations on view vector to see given tangent height

      converge = -1

      DO i = 1, 8

         IF ( ABS(missAltitude - initAlt) < 1.0 ) THEN
            converge = 1
            EXIT
         ENDIF

         returnStatus = Pgs_csc_geoToECR (longitude, latitude, initAlt, &
              earthModel, ecr)

         ray = ecr - posECR

         returnStatus = Pgs_csc_grazingRay(earthModel, posECR, ray, &
              latitude, longitude, missAltitude, slantRange, posNear, posSurf)

      ENDDO

! Check that exit caused by height convergence, rather than loop timing out

      IF ( converge /= 1 ) CALL MLSMessage (MLSMSG_Warning, ModuleName, &
           'Iterations timed out before height converged.')

! Convert final view vector from ECR (to ECI) to SC

      ecrV(1:3) = ray
      ecrV(4:6) = 0.0

      returnStatus = Pgs_csc_ecrToECI (numValues, asciiUTC, time_offset, ecrV, &
           eciV)

      eci = eciV(1:3)

      returnStatus = Pgs_csc_eciToSc (spacecraftid, numValues, asciiUTC, &
           time_offset, eci, finalRay)

      startAngle = Rad2Deg * ACOS (finalRay(1))

!---------------------------
   END SUBROUTINE Scan_start
!---------------------------


!==============
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE Scan
!==============

! $Log$
! Revision 2.2  2004/08/16 17:26:58  perun
! Comment out returnStatus check for initial guess
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.0  2000/09/05 18:55:15  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/02/10 16:54:27  nakamura
! Module that models the MLS scan program.
!
