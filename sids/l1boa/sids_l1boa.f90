!===================================================
PROGRAM sids_l1boa ! End-to-end check on subroutines
!===================================================

   USE Hdf
   USE L1boa
   USE MLSCommon
   USE MLSMessageModule
   USE Orbit
   USE OutputL1B
   USE Read
   USE Sd
   USE SDPToolkit
   USE Time

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program creates a  SIDS L1BOA file.

! Parameters

   INTEGER, PARAMETER :: SD_OUTFILE = 333

! Functions

! Variables

   TYPE ( TAI93_Range_T ) :: times
   TYPE( L1BOAindex_T) :: index

   CHARACTER (LEN=27) :: endTime, startTime
   CHARACTER (LEN=135) :: physicalFilename
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=27), ALLOCATABLE :: mafTime(:), postTimes(:), preTimes(:)
   CHARACTER (LEN=27), ALLOCATABLE :: scanTimes(:,:)

   INTEGER :: alloc_err, dealloc_err, i, lenMAF, lenMIF, mifG, mifRate, mifT
   INTEGER :: minSC, nPhaseG, nPhaseT, numOrb, offsetMAF, opd, orbPerDay
   INTEGER :: postMAF, preMAF, returnStatus, scansPerOrb, sd_id, status
   INTEGER :: SD_VERSION
   INTEGER :: orbitNumber(max_orbits)
   INTEGER, ALLOCATABLE :: numMIFs(:), nV(:)
   INTEGER, ALLOCATABLE :: numValues(:,:)

   REAL, ALLOCATABLE :: scRate(:), scRateT(:)

   REAL(r8) :: altG, altT, homeAlt, homeLat, mfps, orbIncline, spm
   REAL(r8) :: ascTAI(max_orbits), dscTAI(max_orbits)
   REAL(r8), ALLOCATABLE :: mafTAI(:), offsets(:), postTAI(:), preTAI(:)
   REAL(r8), ALLOCATABLE :: scanTAI(:,:)

! Read SIDS input

   CALL Read_uif(altG, altT, times, endTime, homeAlt, homeLat, mifG, mifT, &
                 mifRate, nPhaseG, nPhaseT, offsetMAF, scansPerOrb, startTime)

   ALLOCATE(scRate(mifG), scRateT(mifT), STAT=alloc_err)
   IF ( alloc_err /= 0 ) THEN
      msr = MLSMSG_Allocate // '  scan rate variables.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   CALL Read_scan(mifG, mifT, nPhaseG, nPhaseT, scRate, scRateT)

! Get orbit metadata for entire day

   CALL Orbit_met(startTime, times, ascTAI, dscTAI, numOrb, orbitNumber)

   opd = numOrb+1

! Calculate MAF start times

   ALLOCATE( numMIFs(scansPerOrb), numValues(scansPerOrb,opd), &
           postTAI(scansPerOrb), postTimes(scansPerOrb), preTAI(scansPerOrb), &
           preTimes(scansPerOrb), scanTAI(scansPerOrb,opd), &
           scanTimes(scansPerOrb,opd), STAT=alloc_err )
   IF ( alloc_err /= 0 ) THEN
      msr = MLSMSG_Allocate // '  MAF quantities.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   mfps = DBLE(mifRate)
   spm = 1/mfps

   CALL Time_MAF(startTime, times, homeAlt, homeLat, mifRate, scansPerOrb, &
                 spm, scanTimes, scanTAI, numValues, orbPerDay)

   CALL Time_pre(numValues(:,1), scanTimes(1,1), times%startTime, scansPerOrb, &
                 spm, preTimes, preTAI, numMIFs, preMAF)

   CALL Time_post(numValues(:,orbPerDay-1), scanTimes(1,orbPerDay), &
                  times%endTime, scansPerOrb, spm, postTimes, postTAI, postMAF)

! Check the lengths of the given scans against the minimum # of MIFs per MAF

   minSC = MINVAL( numValues(:,1:orbPerDay-1) )

   IF ( mifG > minSC ) THEN
      print *, 'Desired GHz scan program is too long:'
      print *, 'Minimum number of MIFs per MAF = ', minSC
      print *, 'Number of MIFs in GHz scan = ', mifG
      STOP
   ENDIF

   IF ( mifT > minSC ) THEN
      print *, 'Desired THz scan program is too long:'
      print *, 'Minimum number of MIFs per MAF = ', minSC
      print *, 'Number of MIFs in THz scan = ', mifT
      STOP
   ENDIF

! Calculate one-time values, outside of MAF/orbit loop(s)

   lenMIF = MAXVAL(numValues)

   ALLOCATE( offsets(lenMIF), STAT=alloc_err )
   IF ( alloc_err /= 0 ) THEN
      msr = MLSMSG_Allocate // '  offsets.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   DO i = 1, lenMIF
      offsets(i) = (i-1)/mfps
   ENDDO

   orbIncline = 98.145

! Get filename from PCF

   SD_VERSION = 1

   returnStatus = Pgs_pc_getReference(SD_OUTFILE, SD_VERSION, physicalFilename)

! Open HDF-SD output file and create empty arrays; calculate # of MAFs in the
! day.

   CALL Sd_create(mifG, lenMIF, mifT, offsetMAF, orbPerDay, physicalFilename, &
                  postMAF, preMAF, scansPerOrb, lenMAF)

   print *, 'Num of MAFs in file = ', lenMAF

! After creation, re-open the HDF file and initialize the SD interface for
! writing

   sd_id = sfstart(physicalFilename, DFACC_WRITE)
   IF (sd_id == -1) THEN
      msr = MLSMSG_Fileopen // physicalFilename
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

! Re-arrange time, number of values into single array for entire file

   ALLOCATE( mafTAI(lenMAF), mafTime(lenMAF), nV(lenMAF), STAT=alloc_err )
   IF ( alloc_err /= 0 ) THEN
      msr = MLSMSG_Allocate // '  time variables.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   CALL Time_day(lenMAF, numMIFs, numValues, orbPerDay, postMAF, postTAI, &
                 postTimes, preMAF, preTAI, preTimes, scansPerOrb, scanTAI, &
                 scanTimes, mafTAI, mafTime, nV)

! For first MAF, calculate & write L1BOA data with fill values/settings

   print *, 'Processing MAF 1 ... '

   index%MAFStartTimeUTC = mafTime(1)
   index%MAFStartTimeTAI = mafTAI(1)
   index%noMIFs = nV(1)
   index%counterMAF = 0

   CALL L1boa_fill(altG, altT, ascTAI, dscTAI, sd_id, mifG, mifT, index, 1, &
                   numOrb, offsets, orbIncline, orbitNumber, scRate, scRateT)

! Close the file to set fill values

   status = sfend(sd_id)
   IF (status == -1) THEN
      msr = 'Failed to close file ' // physicalFilename
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

! Re-open the HDF file and initialize the SD interface for writing

   sd_id = sfstart(physicalFilename, DFACC_WRITE)
   IF (sd_id == -1) THEN
      msr = MLSMSG_Fileopen // physicalFilename
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

! For remaining MAFs, calculate & write L1BOA data

   DO i = 2, lenMAF

      print *, 'Processing MAF ', i, ' ... '

      index%MAFStartTimeUTC = mafTime(i)
      index%MAFStartTimeTAI = mafTAI(i)
      index%noMIFs = nV(i)
      index%counterMAF = i-1

      CALL L1boa_nofill(altG, altT, ascTAI, dscTAI, sd_id, mifG, mifT, index, &
                  i, numOrb, offsets, orbIncline, orbitNumber, scRate, scRateT)

   ENDDO

! Terminate access to the SD interface and close the file

   status = sfend(sd_id)
   IF (status == -1) THEN
      msr = 'Failed to close file ' // physicalFilename
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   DEALLOCATE(offsets, numMIFs, numValues, postTAI, postTimes, preTAI, &
              preTimes, scanTAI, scanTimes, scRate, scRateT, mafTAI, mafTime, &
              nV, STAT=dealloc_err )
   IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                                               &deallocation.')

! Detailed description of program
! The program creates a SIDS L1BOA SDS-HDF file.  It uses routines such that
! the fill value/modes are set only on writing the first record.

!=====================
END PROGRAM sids_l1boa
!=====================

!# $Log$
!#
