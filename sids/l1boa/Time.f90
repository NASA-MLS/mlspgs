
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Time
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   USE Scan
   USE TkL1B
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- Time_lat
!                Time_search
!                Time_MAF
!                Time_pre
!                Time_post
!                Time_day

! Remarks:  This module contains subroutines needed to calculate time
!           information for the MLS scan programs.

! Parameters

   INTEGER, PARAMETER :: N_RECORDS = 148

CONTAINS

!---------------------------------------------------------
   SUBROUTINE Time_lat(asciiUTC, homeAlt, tpLat, converge)
!---------------------------------------------------------

! Brief description of subroutine
! This subroutine takes an initial-guess view vector and iteratively runs the
! toolkit GrazingRay routine to find the actual view vector for a given
! tangent height.  The returned value is the latitude for the view vector.

! Arguments

      CHARACTER(LEN=27), INTENT(IN) :: asciiUTC

      REAL(r8), INTENT(IN) :: homeAlt

      INTEGER, INTENT(OUT) :: converge

      REAL(r8), INTENT(OUT) :: tpLat

! Parameters

      INTEGER, PARAMETER :: numValues = 1

! Functions

      INTEGER :: Pgs_csc_geoToECR, Pgs_csc_grazingRay

! Variables

      TYPE( L1BOAsc_T ) :: sc

      CHARACTER (LEN=480) :: msr

      INTEGER :: error, i, returnStatus

      REAL(r8) :: longitude, missAltitude, slantRange
      REAL(r8) :: ecr(3), initGuess(3), posNear(3), posSurf(3), ray(3)
      REAL(r8) :: time_offset(numValues)

      time_offset(numValues) = 0.0

! Get initial guess look-vector, initial s/c position, initial tp info

      CALL Scan_guess(asciiUTC, initGuess)

      ALLOCATE(sc%scECI(lenCoord,numValues), sc%scECR(lenCoord,numValues), &
       sc%scGeocAlt(numValues), sc%scGeocLat(numValues), &
       sc%scGeodAlt(numValues), sc%scGeodLat(numValues), sc%scLon(numValues), &
       sc%scGeodAngle(numValues), sc%scVel(lenCoord,numValues), &
       sc%ypr(lenCoord,numValues), sc%yprRate(lenCoord,numValues), STAT=error)
      IF ( error /= 0 ) THEN
         msr = MLSMSG_Allocate // '  s/c quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL TkL1B_sc(numValues, time_offset, asciiUTC, sc)

      returnStatus = Pgs_csc_grazingRay(earthModel, sc%scECR, initGuess, &
                  tpLat, longitude, missAltitude, slantRange, posNear, posSurf)

! Iterations on view vector to see given tangent height

      converge = -1

      DO i = 1, 8

         IF ( ABS(missAltitude - homeAlt) < 1.0 ) THEN
            converge = 1
            EXIT
         ENDIF

         returnStatus = Pgs_csc_geoToECR(longitude, tpLat, homeAlt, &
                                         earthModel, ecr)

         ray = ecr - sc%scECR(:,1)

         returnStatus = Pgs_csc_grazingRay(earthModel, sc%scECR, ray, tpLat, &
                         longitude, missAltitude, slantRange, posNear, posSurf)

      ENDDO

!-------------------------
   END SUBROUTINE Time_lat
!-------------------------

!-----------------------------------------------------------------------
   SUBROUTINE Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
                          finalTime, finalTAI, finalLat)
!-----------------------------------------------------------------------

! Brief description of program
! This program calculates the MIF starting times for the orbits.

! Arguments

      CHARACTER (LEN=27), INTENT(IN) :: estTime

      INTEGER, INTENT(IN) :: scansPerOrb

      REAL(r8), INTENT(IN) :: homeAlt, homeLat, spm

      CHARACTER (LEN=27), INTENT(OUT) :: finalTime

      REAL(r8), INTENT(OUT) :: finalLat, finalTAI

! Parameters

      CHARACTER (LEN=*), PARAMETER :: GR_ERR = 'Iterations timed out before &
                                               &height converged for '

! Functions

      INTEGER :: Pgs_td_utcToTAI, Pgs_td_taiToUTC

! Variables

      CHARACTER (LEN=27) :: timeF
      CHARACTER (LEN=480) :: msg

      INTEGER :: back, converge, i, returnStatus

      REAL(r8) :: estTAI93, latE, latF, tai

      finalTime = estTime

! Get tp info for estTime

      CALL Time_lat(estTime, homeAlt, latE, converge)
      IF ( converge /= 1 ) THEN
         msg = GR_ERR // estTime
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msg)
      ENDIF

      returnStatus = Pgs_td_utcToTAI (estTime, estTAI93)

! Set search direction

      IF (latE < homeLat) THEN
         back = -1
      ELSE
         back = 1
      ENDIF

! Check determined direction in MIF steps
 
      DO i = 1, scansPerOrb*N_RECORDS

         tai = estTAI93 - back * spm

         returnStatus = Pgs_td_taiToUTC(tai, timeF)

         CALL Time_lat(timeF, homeAlt, latF, converge)
         IF ( converge /= 1 ) THEN
            msg = GR_ERR // timeF
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msg)
            EXIT
         ENDIF

         IF ( ABS(latF - homeLat) >= ABS(latE - homeLat) ) THEN
            EXIT
         ELSE
            latE = latF
            estTAI93 = tai
            finalTime = timeF
         ENDIF

      ENDDO

! Add an MIF, if closest lat found is below homeLat

      IF (latE < homeLat) THEN

         tai = estTAI93 + spm

         returnStatus = Pgs_td_taiToUTC(tai, timeF)

         CALL Time_lat(timeF, homeAlt, latF, converge)
         IF ( converge /= 1 ) THEN
            msg = GR_ERR // timeF
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msg)
         ENDIF

         latE = latF
         estTAI93 = tai
         finalTime = timeF

      ENDIF

      finalLat = latE
      finalTAI = estTAI93

!----------------------------
   END SUBROUTINE Time_search
!----------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Time_MAF(startTime, times, homeAlt, homeLat, mifRate, &
                      scansPerOrb, spm, scanTimes, scanTAI, numMIFs, orbPerDay)
!------------------------------------------------------------------------------

! Brief description of program
! This program calculates the starting times for the orbits and their scans,
! and the number of MIFs per MAF.

! Arguments

      TYPE( TAI93_Range_T ), INTENT(IN) :: times

      CHARACTER (LEN=27), INTENT(IN) :: startTime

      INTEGER, INTENT(IN) :: mifRate, scansPerOrb

      REAL(r8), INTENT(IN) :: homeAlt, homeLat, spm

      CHARACTER (LEN=27), INTENT(OUT) :: scanTimes(:,:)

      INTEGER, INTENT(OUT) :: orbPerDay
      INTEGER, INTENT(OUT) :: numMIFs(:,:)

      REAL(r8), INTENT(OUT) :: scanTAI(:,:)

! Parameters

      CHARACTER (LEN=*), PARAMETER ::GR_ERR = 'Iterations timed out before &
                                            &height converged for '

      REAL(r8), PARAMETER :: minPerOrb = 98.9

! Functions

      INTEGER :: Pgs_td_taiToUTC

! Variables

      CHARACTER (LEN=27) :: estTime, timeF
      CHARACTER (LEN=27), ALLOCATABLE :: orbTime(:)
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: alloc_err, ascend, converge, dealloc_err, delMIF, delOrb, i, j
      INTEGER :: k, opd, returnStatus

      REAL(r8) :: delAngle, delMAF, delTime, estLat, latF, orbRate, startLat
      REAL(r8) :: tai
      REAL(r8), ALLOCATABLE :: orbTAI(:)

      orbRate = minPerOrb * 60 / (2*PI)

      delOrb = AINT(minPerOrb * 60 * mifRate)

      opd = SIZE(scanTimes,2)

      ALLOCATE( orbTime(opd), orbTAI(opd), STAT=alloc_err )
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  orbit variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get tp info for startTime

      CALL Time_lat(startTime, homeAlt, startLat, converge)
      IF ( converge /= 1 ) THEN
         msg = GR_ERR // startTime
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msg)
      ENDIF

! Find whether ascending or descending at startTime

      ascend = -1

      tai = times%startTime + spm

      returnStatus = Pgs_td_taiToUTC(tai, timeF)

      CALL Time_lat(timeF, homeAlt, latF, converge)
      IF ( converge /= 1 ) THEN
         msg = GR_ERR // timeF
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msg)
      ENDIF

      IF ( (latF - startLat) > 0 ) ascend = 1

! Estimate MIF time of nearest ascending homeLat crossing

      IF (ascend == 1) THEN

         IF ( homeLat >= startLat ) THEN
            delAngle = homeLat - startLat
         ELSE
            delAngle = 2*PI + homeLat - startLat
         ENDIF

      ELSE

         delAngle = PI + homeLat + startLat

      ENDIF

      delTime = delAngle * orbRate

      delMIF = AINT(delTime * mifRate)

      tai = times%startTime + delMIF * spm

      returnStatus = Pgs_td_taiToUTC(tai, estTime)

! Find MIF time of nearest ascending homeLat crossing

      orbTAI = 0.0

      CALL Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
                       orbTime(1), orbTAI(1), estLat)

! Estimate time of remaining orbit crossings

      DO i = 2, opd

         tai = orbTAI(i-1) + delOrb * spm

         IF ( tai > times%endTime) THEN
            orbTAI(i) = 0.0
            orbPerDay = i-1
            EXIT
         ENDIF

         returnStatus = Pgs_td_taiToUTC( tai, estTime)

         CALL Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
                          orbTime(i), orbTAI(i), estLat)

         IF ( orbTAI(i) > times%endTime) THEN
            orbTAI(i) = 0.0
            orbPerDay = i-1
            EXIT
         ENDIF

      ENDDO

! Divide orbits into 240 MAFs

      scanTAI = 0.0

      DO i = 1, orbPerDay

         IF ( orbTAI(i+1) /= 0.0 ) THEN

            delMAF = ( orbTAI(i+1) - orbTAI(i) ) / scansPerOrb

            DO j = 1, scansPerOrb

               delMIF = ANINT( (j-1) * delMAF * mifRate )
               scanTAI(j,i) = orbTAI(i) + delMIF * spm
               returnStatus = Pgs_td_taiToUTC( scanTAI(j,i), scanTimes(j,i) )

            ENDDO

         ELSE

            scanTAI(1,i) = orbTAI(i)
            scanTimes(1,i) = orbTime(i)
            EXIT

         ENDIF

      ENDDO

! Calculate number of MIFs

      DO i = 1, orbPerDay-1

         DO j = 1, scansPerOrb

            IF ( j /= scansPerOrb) THEN

               DO k = 1, N_RECORDS
                  IF ( scanTAI(j,i) + k*spm > scanTAI(j+1,i) ) EXIT
                  numMIFs(j,i) = k
               ENDDO

            ELSE 

               DO k = 1, N_RECORDS
                 IF ( ( scanTAI(j,i) + k*spm ) > scanTAI(1,i+1) ) EXIT
                 numMIFs(j,i) = k
               ENDDO

            ENDIF

         ENDDO

      ENDDO

      DEALLOCATE( orbTime, orbTAI, STAT=dealloc_err )
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                                  'Failed deallocation of orbit variables.')

!-------------------------
   END SUBROUTINE Time_MAF
!-------------------------

!-----------------------------------------------------------------------
   SUBROUTINE Time_pre(numValues, orbTime, startTAI, scansPerOrb, spm, &
                       preTimes, preTAI, numMIFs, preMAF)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine fills in the time gap for a partial orbit at the beginning of
! a file.  It duplicates the scan pattern for the first MAF of the first full
! orbit.

! Arguments
 
      CHARACTER(LEN=27), INTENT(IN) :: orbTime

      INTEGER, INTENT(IN) :: scansPerOrb
      INTEGER, INTENT(IN) :: numValues(scansPerOrb)

      REAL(r8), INTENT(IN) :: spm, startTAI

      CHARACTER(LEN=27), INTENT(OUT) :: preTimes(scansPerOrb)

      INTEGER, INTENT(OUT) :: preMAF
      INTEGER, INTENT(OUT) :: numMIFs(scansPerOrb)

      REAL(r8), INTENT(OUT) :: preTAI(scansPerOrb)

! Parameters

! Functions

      INTEGER :: Pgs_td_utcToTAI, Pgs_td_taiToUTC
  
! Variables

      CHARACTER(LEN=27) :: backTimes(scansPerOrb)

      INTEGER :: i, returnStatus

      REAL(r8) :: mafTAI, orbTAI
      REAL(r8) :: backTAI(scansPerOrb)

      preMAF = 0

! Calculate MAF times prior to first full orbit

      returnStatus = Pgs_td_utcToTAI (orbTime, orbTAI)

      mafTAI = orbTAI

      DO i = 1, scansPerOrb

         mafTAI = mafTAI - numValues(scansPerOrb - (i-1))*spm

         IF ( mafTAI < startTAI ) EXIT

         returnStatus = Pgs_td_taiToUTC( mafTAI, backTimes(i) )
         backTAI(i) = mafTAI
         preMAF = i

      ENDDO

! Reverse numbering of MAF times, numValues

      DO i = 1, preMAF

         preTimes(i) = backTimes( preMAF - (i-1) )
         preTAI(i)   = backTAI (preMAF - (i-1) )
         numMIFs(i) = numValues( scansPerOrb - (preMAF - i) )

      ENDDO

!--------------------------
    END SUBROUTINE Time_pre
!--------------------------

!-----------------------------------------------------------------------
    SUBROUTINE Time_post(numValues, orbTime, endTAI, scansPerOrb, spm, &
                         postTimes, postTAI, postMAF)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine fills in the time gap for a partial orbit at the end of a
! file.  It duplicates the scan pattern for the last MAF of the last full
! orbit.

! Arguments

      CHARACTER(LEN=27), INTENT(IN) :: orbTime

      INTEGER, INTENT(IN) :: scansPerOrb
      INTEGER, INTENT(IN) :: numValues(scansPerOrb)

      REAL(r8), INTENT(IN) :: endTAI, spm

      CHARACTER(LEN=27), INTENT(OUT) :: postTimes(scansPerOrb)

      INTEGER, INTENT(OUT) :: postMAF

      REAL(r8), INTENT(OUT) :: postTAI(scansPerOrb)

! Parameters

! Functions

      INTEGER :: Pgs_td_utcToTAI, Pgs_td_taiToUTC

! Variables

      INTEGER :: i, returnStatus

      REAL(r8) :: orbTAI, mafTAI

      postMAF = 0

! Calculate MAF times after the last full orbit

      returnStatus = Pgs_td_utcToTAI (orbTime, orbTAI)

      mafTAI = orbTAI

      DO i = 1, scansPerOrb

         mafTAI = mafTAI + numValues(i)*spm

         IF ( mafTAI > endTAI ) EXIT

         returnStatus = Pgs_td_taiToUTC( mafTAI, postTimes(i+1) )
         postTAI(i+1) = mafTAI

         postMAF = i + 1

      ENDDO

      postTimes(1) = orbTime
      postTAI(1) = orbTAI

!--------------------------
   END SUBROUTINE Time_post
!--------------------------

!---------------------------------------------------------------------------
   SUBROUTINE Time_day(lenMAF, numMIFs, numValues, orbPerDay, postMAF, &
                       postTAI, postTimes, preMAF, preTAI, preTimes, &
                       scansPerOrb, scanTAI, scanTimes, mafTAI, mafTime, nV)
!---------------------------------------------------------------------------

! Brief description of subroutine
! This program reorganizes the various time arrays into a single array for the
! entire day/file.

! Arguments

      INTEGER, INTENT(IN) :: lenMAF, postMAF, preMAF, orbPerDay, scansPerOrb

      CHARACTER (LEN=27), INTENT(IN) :: postTimes(postMAF)
      CHARACTER (LEN=27), INTENT(IN) :: preTimes(preMAF)
      CHARACTER (LEN=27), INTENT(IN) :: scanTimes(scansPerOrb,orbPerDay)

      INTEGER, INTENT(IN) :: numMIFs(preMAF)
      INTEGER, INTENT(IN) :: numValues(scansPerOrb,orbPerDay)

      REAL(r8), INTENT(IN) :: postTAI(postMAF)
      REAL(r8), INTENT(IN) :: preTAI(preMAF)
      REAL(r8), INTENT(IN) :: scanTAI(scansPerOrb,orbPerDay)

      CHARACTER (LEN=27), INTENT(OUT) :: mafTime(lenMAF)

      INTEGER, INTENT(OUT) :: nV(lenMAF)

      REAL(r8), INTENT(OUT) :: mafTAI(lenMAF)

! Parameters

! Functions

! Variables

   INTEGER :: i, j, k

! MAFs prior to first full orbit

   DO i = 1, preMAF
      IF ( i > lenMAF ) EXIT
      nV(i) = numMIFs(i)
      mafTime(i) = preTimes(i)
      mafTAI(i) = preTAI(i)
   ENDDO

! MAFs for full orbits

   DO i = 1, orbPerDay-1
      DO j = 1, scansPerOrb

         k = j + (i-1)*scansPerOrb + preMAF
         IF ( k > lenMAF ) EXIT

         nV(k) = numValues(j,i)
         mafTime(k) = scanTimes(j,i)
         mafTAI(k) = scanTAI(j,i)

      ENDDO
   ENDDO

! MAFs after last full orbit

   DO i = 1, postMAF

      j = i + (orbPerDay-1)*scansPerOrb + preMAF
      IF ( j > lenMAF ) EXIT

      nV(j) = numValues(i,orbPerDay-1)
      mafTime(j) = postTimes(i)
      mafTAI(j) = postTAI(i)

   ENDDO

!-------------------------
   END SUBROUTINE Time_day
!-------------------------

!==============
END MODULE Time
!==============

!# $Log$
!#
