
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Read
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

! Subroutines -- Read_uif
!                Read_scan

! Remarks:  This module contains subroutines related to reading the SIDS L1BOA
!           User Input File.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: labelFmt = "(A77)"

   INTEGER, PARAMETER :: SIDS_UIF = 322

CONTAINS

!---------------------------------------------------------------------------
   SUBROUTINE Read_uif(altG, altT, times, endTime, homeAlt, homeLat, mifG, &
                       mifT, mifRate, nPhaseG, nPhaseT, offsetMAF, &
                       scansPerOrb, startTime)
!---------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reads the first (non-scan) portion of the UIF and skims the
! scan portion for information needed to allocate variables.

! Arguments

      TYPE( TAI93_Range_T ), INTENT(OUT) :: times

      CHARACTER (LEN=27), INTENT(OUT) :: endTime, startTime

      INTEGER, INTENT(OUT) :: mifG, mifT, nPhaseG, nPhaseT
      INTEGER, INTENT(OUT) :: mifRate, offsetMAF, scansPerOrb

      REAL(r8), INTENT(OUT) :: altG, altT, homeAlt, homeLat

! Parameters

      CHARACTER (LEN=*), PARAMETER :: dtFmt = "(A10, 1X, A15)"

! Functions

      INTEGER :: Pgs_io_gen_openF, Pgs_io_gen_closeF, Pgs_td_utcToTAI

! Variables

      CHARACTER (LEN=10) :: startDate, endDate
      CHARACTER (LEN=15) :: startUTC, endUTC
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=77) :: label
      CHARACTER (LEN=480) :: msg

      INTEGER :: dur, i, ios, j, processUIF, returnStatus, version

      version = 1

! Open the UIF as a generic file for reading

      returnStatus = Pgs_io_gen_openF (SIDS_UIF, PGSd_IO_Gen_RSeqFrm, 0, &
                                       processUIF, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         print *, 'Error opening UIF:  ', mnemonic
         print *, msg
      ENDIF
 
! Read date/time input

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      IF (ios /= 0) print *, 'Error reading UIF file.'
      READ(UNIT=processUIF, IOSTAT=ios, FMT=dtFmt) startDate, startUTC

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=dtFmt) endDate, endUTC

! Convert to asciiUTC, TAI formats

      startTime = startDate // 'T' // startUTC // 'Z'
      endTime = endDate // 'T' // endUTC // 'Z'

      returnStatus = Pgs_td_utcToTAI(startTime, times%startTime)
      returnStatus = Pgs_td_utcToTAI(endTime, times%endTime)

! Read rest of non-scan input

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) mifRate

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) scansPerOrb

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) offsetMAF

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) homeLat

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) homeAlt

! Read GHz scan program input

      DO i = 1, 4
         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         IF ( label(17:24) == 'Geod Alt' ) EXIT
      ENDDO

      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) altG

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) nPhaseG

! Calculate GHz scan length

      mifG = 0

      DO i = 1, nPhaseG

         DO j = 1, 4
            READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
            IF ( label(9:11) == 'Dur' ) EXIT
         ENDDO

         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) dur

         mifG = mifG + dur

      ENDDO
      
! Read THz scan program input

      DO i = 1, 4
         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         IF ( label(17:24) == 'Geod Alt' ) EXIT
      ENDDO

      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) altT

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      READ(UNIT=processUIF, IOSTAT=ios, FMT=*) nPhaseT

! Calculate THz scan length

      mifT = 0

      DO i = 1, nPhaseT

         DO j = 1, 4
            READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
            IF ( label(9:11) == 'Dur' ) EXIT
         ENDDO

         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) dur

         mifT = mifT + dur

      ENDDO

! Close UIF

      returnStatus = Pgs_io_gen_closeF (processUIF)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         print *, 'Error closing UIF:  ', mnemonic
         print *, msg
      ENDIF

! Convert km (input) to m (used by Toolkit) for altitudes

      homeAlt = homeAlt * 1000
      altG = altG * 1000
      altT = altT * 1000

!-------------------------
   END SUBROUTINE Read_uif
!-------------------------

!-------------------------------------------------------------------------
   SUBROUTINE Read_scan(mifG, mifT, nPhaseG, nPhaseT, scanRate, scanRateT)
!-------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reads the second (scan) portion of the UIF and calculates the
! GHz and THz scan rates.

! Arguments

      INTEGER, INTENT(IN) :: mifG, mifT, nPhaseG, nPhaseT

      REAL, INTENT(OUT) :: scanRate(mifG), scanRateT(mifT)
 
! Parameters

! Functions

      INTEGER :: Pgs_io_gen_openF, Pgs_io_gen_closeF

! Variables

      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=77) :: label
      CHARACTER (LEN=480) :: msg

      INTEGER :: i, ios, processUIF, returnStatus, sum, version
      INTEGER :: dur(nPhaseG)
      INTEGER :: durT(nPhaseT)

      REAL(r8) :: rate(nPhaseG)
      REAL(r8) :: rateT(nPhaseT)

      version = 1

! Open the UIF as a generic file for reading

      returnStatus = Pgs_io_gen_openF (SIDS_UIF, PGSd_IO_Gen_RSeqFrm, 0, &
                                       processUIF, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         print *, 'Error opening UIF:  ', mnemonic
         print *, msg
      ENDIF
 
! Find scan program phase section

      DO i = 1, 20
         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         IF ( label(16:22) == 'Phases' ) EXIT
      ENDDO

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label

! Read GHz scan information

      sum = 0

      DO i = 1, nPhaseG

         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) rate(i)

         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) dur(i)

! Calculate GHz rate information

         scanRate(sum+1:sum+dur(i)) = rate(i)

         sum = sum + dur(i)

      ENDDO

! Find scan program phase section for THz

      DO i = 1, 5
         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         IF ( label(16:22) == 'Phases' ) EXIT
      ENDDO

      READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label

! Read THz scan information

      sum = 0

      DO i = 1, nPhaseT

         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) rateT(i)

         READ(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
         READ(UNIT=processUIF, IOSTAT=ios, FMT=*) durT(i)

! Calculate THz rate information

         scanRateT(sum+1:sum+durT(i)) = rateT(i)

         sum = sum + durT(i)

      ENDDO

! Close UIF

      returnStatus = Pgs_io_gen_closeF (processUIF)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         print *, 'Error closing UIF:  ', mnemonic
         print *, msg
      ENDIF

!--------------------------
   END SUBROUTINE Read_scan
!--------------------------

!==============
END MODULE Read
!==============

!# $Log$
!#
