
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Orbit
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   USE OutputL1B, only: LENG, LENT
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

   INTEGER :: numOrb, orbitNumber(max_orbits)
   REAL :: scanRate(lenG), scanRateT(lenT)
   REAL(r8) :: altG, altT, orbIncline, ascTAI(max_orbits), dscTAI(max_orbits)

! Contents:

! Subroutines -- Orbit_met
!                Orbit_init

! Remarks:  This module contains subroutines to get orbit "metadata," and also
!           initializes values outside the MAF loop on L1BOA processing.

CONTAINS

!----------------------------------------------------------------------
   SUBROUTINE Orbit_met(startTime, times, ascTAI, descTAI, numOrbits, &
                        orbitNumber)
!----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine uses the toolkit getEphMet routine to calculate useful
! information related to the s/c orbit.

! Arguments

      TYPE( TAI93_Range_T ), INTENT(IN) :: times

      CHARACTER (LEN=27), INTENT(IN) :: startTime

      INTEGER, INTENT(OUT) :: numOrbits
      INTEGER, INTENT(OUT) :: orbitNumber(max_orbits)

      REAL(r8), INTENT(OUT) :: ascTAI(max_orbits), descTAI(max_orbits)

! Parameters

! Functions

      INTEGER :: Pgs_eph_getEphMet, Pgs_td_timeInterval, Pgs_td_utcToTAI

! Variables

      CHARACTER (LEN=27) :: orbitAscendTime(max_orbits)
      CHARACTER (LEN=27) :: orbitDescendTime(max_orbits)
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: alloc_err, dealloc_err, i, num_points, returnStatus

      REAL(r8) :: deltaTAI
      REAL(r8) :: orbitDownLongitude(max_orbits)
      REAL(r8), ALLOCATABLE :: offsets(:)

! Find time interval covered by file

      returnStatus = Pgs_td_timeInterval(times%startTime, times%endTime, &
                                         deltaTAI)

      num_points = AINT(deltaTAI/60.0)

      ALLOCATE(offsets(num_points), STAT = alloc_err)
      IF (alloc_err /= 0) THEN
         msr = MLSMSG_Allocate // '  offsets.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DO i = 1, num_points
         offsets(i) = i * 60.0
      ENDDO

! Get orbit metadata

      returnStatus = Pgs_eph_getEphMet(spacecraftId, num_points, startTime, &
                            offsets, numOrbits, orbitNumber, orbitAscendTime, &
                            orbitDescendTime, orbitDownLongitude)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = 'Routine getEphMet, ' // mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE(offsets, STAT=dealloc_err)
      IF (dealloc_err /= 0) CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                                           'Failed deallocation of offsets.')

! Convert values to forms more useful to software

      DO i = 1, numOrbits
         orbitNumber(i) = i-1
         returnStatus = Pgs_td_utcToTAI( orbitAscendTime(i), ascTAI(i) )
         returnStatus = Pgs_td_utcToTAI( orbitDescendTime(i), descTAI(i) )
      ENDDO

!--------------------------
   END SUBROUTINE Orbit_met
!--------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE Orbit_init(times, UTC_start, altG, altT, ascTAI, dscTAI, &
                         numOrb, orbIncline, orbitNumber, scanRate, scanRateT)
!-----------------------------------------------------------------------------

! Brief description of subroutine 
! This subroutine defines input parameters for the SIDS L1BOA code and also
! sets up values needed by the production code which won't change while looping
! through the MAFs.

! Arguments

      TYPE ( TAI93_Range_T ), INTENT(IN) :: times

      CHARACTER (LEN=27), INTENT(IN) :: UTC_start

      INTEGER, INTENT(OUT) :: numOrb
      INTEGER, INTENT(OUT) :: orbitNumber(max_orbits)

      REAL, INTENT(OUT) :: scanRate(lenG)
      REAL, INTENT(OUT) :: scanRateT(lenT)

      REAL(r8), INTENT(OUT) :: altG, altT, orbIncline
      REAL(r8), INTENT(OUT) :: ascTAI(max_orbits), dscTAI(max_orbits)

! Parameters

! Functions

! Variables

! Input data

      altG = 2500.0
      altT = 15000.0

      scanRate(1) = 0.0
      scanRate(2:63) = 0.0279213
      scanRate(64:120) = 0.0976263
 
      scanRateT(1:54) = 0.030344496
      scanRateT(55:114) = 0.058600044

      orbIncline = 98.145

! Get orbit metadata for entire day

      CALL Orbit_met(UTC_start, times, ascTAI, dscTAI, numOrb, orbitNumber)

!---------------------------
   END SUBROUTINE Orbit_init
!---------------------------

!===============
END MODULE Orbit
!===============

! $Log$
! Revision 2.3  2001/10/12 22:11:05  livesey
! Tidied things up a bit, added scVelECR, but not filled yet
!
! Revision 2.2  2001/02/01 18:28:04  perun
! Fixed typo
!
! Revision 2.1  2001/02/01 18:20:13  perun
! Added data declarations outside routines
!
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/15 18:30:29  nakamura
! Removed superfluous USE MLSL1Common statement; added Contents section; moved _init subroutine here from TkL1B.
!
