
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Mc
!===============================================================================

   USE MLSCommon
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

! Contents:

! Subroutines -- Mc_aux
!                Mc_phi

! Remarks:  This module contains subroutines to compute the master coordinates
!           for the s/c and tangent point records.

CONTAINS

!---------------------------------------------------
   SUBROUTINE Mc_aux (asciiUTC, scECI, scGeocLat, q)
!---------------------------------------------------

! Brief description of subroutine
! This subroutine computes q, an auxilliary vector used in the calculation of
! the master coordinate.

! Arguments

      CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

      REAL, INTENT(IN) :: scGeocLat

      REAL(r8), INTENT(IN) :: scECI(3)

      REAL(r8), INTENT(OUT) :: q(3)

! Parameters

      INTEGER, PARAMETER :: nV = 1

      REAL(r8), PARAMETER :: time_offset = 0.0

! Functions

      INTEGER :: Pgs_csc_orbToECI

! Variables

      INTEGER :: i, returnStatus

      REAL(r8) :: l
      REAL(r8) :: dist(2)
      REAL(r8) :: aECI(3), aux1(3), aux1ECI(3), aux2(3), aux2ECI(3)

! Construct two auxilliary vectors -- the first pointing directly ahead 
! along the s/c orbit, and the second pointing 45 degrees downward.

      aux1(1) = 1.0
      aux1(2) = 0.0
      aux1(3) = 0.0

      aux2(1) = 1.0
      aux2(2) = 0.0
      aux2(3) = 1.0

! Transform the auxilliary vectors from orbital to ECI coordinates

      returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
                                      time_offset, aux1, aux1ECI)

      returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
                                      time_offset, aux2, aux2ECI)

! Define l & a, where l is the distance to the equator along the direction of
! one of the auxilliary vectors from the s/c, and a is the auxilliary vector
! for which l was defined.

      IF ( aux1ECI(3) == 0 ) THEN

! If aux1 has a zero z-component, set l & a for aux2, if its z-component /= 0

         IF ( aux2ECI(3) /= 0 ) THEN
            l = -scECI(3)/aux2ECI(3)
            aECI = aux2ECI
         ENDIF

      ELSE

! If aux1(3) is not zero, but aux2(3) is, then set l & a for aux1

         IF ( aux2ECI(3) == 0 ) THEN

            l = -scECI(3)/aux1ECI(3)
            aECI = aux1ECI

         ELSE

! If both aux1(3) and aux2(3) are non-zero, choose the one which gives the
! minimum absolute value of l.

            dist(1) = -scECI(3)/aux1ECI(3)
            dist(2) = -scECI(3)/aux2ECI(3)

            IF ( ABS(dist(1)) < ABS(dist(2)) ) THEN
               l = dist(1)
               aECI = aux1ECI
            ELSE
               l = dist(2)
               aECI = aux2ECI
            ENDIF

         ENDIF

      ENDIF
   
! Define the vector q = scECI + la, such that its z-component = 0

      q(1:2) = scECI(1:2) + l*aECI(1:2)
      q(3) = 0.0

! Modify q, depending on which hemisphere the s/c is in, and the sign of l

      IF ( scGeocLat < 0.0 ) THEN
         IF ( l < 0 ) q = -q
      ELSE
         IF ( l > 0 ) q = -q
      ENDIF

! Normalize q

      q = q / SQRT( q(1)**2 + q(2)**2 )
 
!-----------------------
   END SUBROUTINE Mc_aux
!-----------------------

!------------------------------------------------------------------------------
   SUBROUTINE Mc_phi(ascTAI, dscTAI, dotVec, geocLat, nV, numOrb, orbIncline, &
                     orbitNumber, q, timeTAI, time_offset, geodAngle)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine computes phi, the master coordinate for the spacecraft and
! tangent point records.

! Arguments

      INTEGER, INTENT(IN) :: nV, numOrb
      INTEGER, INTENT(IN) :: orbitNumber(:)

      REAL, INTENT(IN) :: geocLat(nV)

      REAL(r8), INTENT(IN) :: orbIncline, timeTAI
      REAL(r8), INTENT(IN) :: q(3)
      REAL(r8), INTENT(IN) :: time_offset(nV)
      REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:)
      REAL(r8), INTENT(IN) :: dotVec(3,nV)

      REAL, INTENT(OUT) :: geodAngle(nV)

! Parameters

! Functions

      INTEGER :: Pgs_csc_getEarthFigure

! Variables

      INTEGER :: i, j, returnStatus, scOrb

      REAL(r8) :: a, asciiTAI, b, cSq, equatRad_a, orbRad, phiMin, polarRad_c
      REAL(r8) :: cosPhi(nV), gamma(nV), phi(nV), sinPhi(nV)
      REAL(r8) :: s(3,nV)

! Read a & b from earthfigure.dat

      returnStatus = Pgs_csc_getEarthFigure(earthModel, equatRad_a, polarRad_c)

      a = equatRad_a/1000
      b = polarRad_c/1000

! Calculate C-squared

      orbRad= Deg2Rad * (orbIncline - 90)

      cSq = (1 + (TAN(orbRad)**2)) * (a**2)*(b**2)/(a**2 + &
                                                  &(b**2)*(TAN(orbRad)**2))

! Set s = normalized dotVec

      DO i = 1, nV

         s(:,i) = dotVec(:,i) / SQRT(dotVec(1,i)**2 + dotVec(2,i)**2 + &
                                                     &dotVec(3,i)**2)

! Calculate the geocentric angle as a number of radians between 0 and PI

         gamma(i) = ACOS( q(1)*s(1,i) + q(2)*s(2,i) )

! Place angle between PI and 2*PI, if lat indicates Southern Hemisphere

         IF ( geocLat(i) < 0.0 ) gamma(i) = 2*PI - gamma(i)

! If |gamma| <= 45 deg of the equator, calculate phi using the SIN equation

         IF ( (gamma(i) <= PI/4 ) .OR. ( (gamma(i) >= 3*PI/4) .AND. &
              &(gamma(i) <= 5*PI/4) ) .OR. (gamma(i) >= 7*PI/4) ) THEN

            sinPhi(i) = SQRT( (a**4)*(SIN(gamma(i))**2)/( (cSq**2)*&
                             &(COS(gamma(i))**2) + (a**4)*(SIN(gamma(i))**2) ) )

            phi(i) = ASIN(sinPhi(i))

         ELSE

! If gamma is within 45 deg of a pole, calculate phi using the COS equation

            cosPhi(i) = SQRT( (cSq**2)*(COS(gamma(i))**2)/( (cSq**2)*&
                              &COS(gamma(i))**2 + (a**4)*(SIN(gamma(i))**2) ) )

            phi(i) = ACOS(cosPhi(i))

         ENDIF

! Place phi in same quadrant as gamma

         IF ( (gamma(i) > PI/2) .AND. (gamma(i) < PI) ) phi(i) = PI - phi(i)
         IF ( (gamma(i) > PI) .AND. (gamma(i) < 3*PI/2) ) phi(i) = phi(i) + PI
         IF ( gamma(i) > 3*PI/2 ) phi(i) = 2*PI - phi(i)

! Make phi cumulative over orbits

         asciiTAI = timeTAI + time_offset(i)

         DO j = 1, numOrb
            IF ( asciiTAI > ascTAI(j) ) scOrb = orbitNumber(j)
         ENDDO

         phi(i) = (scOrb-1)*2*PI + phi(i)

         IF ( asciiTAI > dscTAI(scOrb+1) ) THEN
            phiMin = scOrb*2*PI - PI
         ELSE
            phiMin = (scOrb-1)*2*PI - PI
         ENDIF

         IF ( phi(i) < phiMin ) phi(i) = phi(i) + 2*PI

      ENDDO

! Convert to degrees for output

      geodAngle = Rad2Deg * phi

!--------------------
END SUBROUTINE Mc_phi
!--------------------

!============
END MODULE Mc
!============

! $Log$
