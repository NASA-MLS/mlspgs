
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE OpenInit 
!===============================================================================

   USE L3CF
   USE MLSCF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE PCFHdr
   USE SDPToolkit
   USE dates_module
   use GETCF_M, only: GetCF, InitGetCF
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definitions -- PCFData_T
! Subroutines -- GetPCFParameters
!                SetProcessingWindow
!                OpenAndInitialize
!                AvgOrbPeriod

! Remarks:  This is a prototype module for the routines needed for the L3 Daily
! Open/Init task.

! Parameters

! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFData_T

     ! cycle # of processing run

     CHARACTER (LEN=4) :: cycle

     ! version string in PCF output file names

     CHARACTER (LEN=15) :: outputVersion	! output files

     ! first & last days of input & output windows

     CHARACTER (LEN=8) :: l2EndDay		! last day of l2 to read
     CHARACTER (LEN=8) :: l2StartDay		! first day of l2 to read
     CHARACTER (LEN=8) :: l3EndDay		! last day of l3 to write
     CHARACTER (LEN=8) :: l3StartDay		! first day of l3 to write

     ! name of the log file (without the path)

     CHARACTER (LEN=FileNameLen) :: logGranID

   END TYPE PCFData_T

CONTAINS

!------------------------------------
   SUBROUTINE GetPCFParameters(l3pcf)
!------------------------------------

! Brief description of subroutine
! This subroutine retrieves the User-Defined Runtime Parameters from the PCF
! and stores them in variables used in the MLSL3 code.

! Arguments

      TYPE( PCFData_T ), INTENT(OUT) :: l3pcf

! Parameters

      CHARACTER (LEN=*), PARAMETER :: UDRP_ERR = 'Failed to get PCF runtime &
                                                 &parameter :  '

! Functions

      INTEGER, EXTERNAL :: pgs_pc_getConfigData

! Variables

      CHARACTER (LEN=17) :: range 
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: name

      INTEGER ::  indx, mlspcf_log, returnStatus, version

! Retrieve values set in PCF & assign them to variables.

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_OutputVersion, &
                                          l3pcf%outputVersion)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'OutputVersion'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_Cycle, l3pcf%cycle)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'Cycle'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_L2DayRange, range)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'L2DayRange'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      l3pcf%l2StartDay = range(1:8)
      l3pcf%l2EndDay = range(10:17)

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_RangDays, range)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'RangDays'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      l3pcf%l3StartDay = range(1:8)
      l3pcf%l3EndDay = range(10:17)

! Get the name of the log file from the PCF

      version = 1
      mlspcf_log = 10101

      returnStatus = Pgs_pc_getReference(mlspcf_log, version, name)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                        ModuleName, 'Error retrieving log file name from PCF.')
      indx = INDEX(name, '/', .TRUE.)
      l3pcf%logGranID = name(indx+1:)

!---------------------------------
   END SUBROUTINE GetPCFParameters
!---------------------------------

!-------------------------------------------------------------
   SUBROUTINE SetProcessingWindow (start, end, dates, numDays)
!-------------------------------------------------------------

! Brief description of subroutine
! This subroutine finds each DOY for the L3 processing window and calculates
! the number of days in the window.

! Arguments

      CHARACTER (LEN=8), INTENT(IN) :: end, start

      CHARACTER (LEN=8), INTENT(OUT) :: dates(:)

      INTEGER, INTENT(OUT) :: numDays

! Parameters

! Functions

! Variables

      CHARACTER (LEN=3) :: cDOY
      CHARACTER (LEN=4) :: cYr

      INTEGER :: i, iDOY, iEnd, iYr, lastDOY

! Initializations outside the DO loop

      READ(start(1:4), '(I4)') iYr
      lastDOY = days_in_year(iYr)
      READ( start(6:8), '(I3)' ) iDOY
      READ( end(6:8), '(I3)' ) iEnd

      DO i = 1, maxWindow

! Increment DOY

         iDOY = iDOY + 1

! Check that the incrementation hasn't crossed a year boundary

         IF (iDOY > lastDOY) THEN
            iYr = iYr + 1
            iDOY = 1
         ENDIF

! Check that the incrementation hasn't crossed the window boundary

         IF (iDOY > iEnd) EXIT

! Convert Yr back to CHARACTER

         WRITE( cYr, '(I4)' ) iYr

! Convert DOY back to CHARACTER, adding leading zeroes if necessary

         WRITE( cDOY, '(I3)' ) iDOY

         IF (iDOY < 10) THEN
            cDOY = '00' // cDOY(3:3)
         ELSE IF (iDOY < 100) THEN
            cDOY = '0' // cDOY(2:3)
         ENDIF

! Concatenate them into a date

        dates(i+1) = cYr // '-' // cDOY

      ENDDO

      numDays = i
      dates(1) = start
      dates(numDays) = end

!------------------------------------
   END SUBROUTINE SetProcessingWindow
!------------------------------------

!-----------------------------------------------------------------------
   SUBROUTINE OpenAndInitialize (l3pcf, cf, l3cf, cfDef, anText, avgPer)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Open/Init task in the MLSL3 program.

! Arguments

      TYPE( L3CFDef_T ), INTENT(OUT) :: cfDef

      TYPE( Mlscf_T ), INTENT(OUT) :: cf

      TYPE( PCFData_T ), INTENT(OUT) :: l3pcf

      TYPE( L3CFProd_T ), POINTER :: l3cf(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

      REAL(r8), POINTER :: avgPer(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=8) :: dates(maxWindow)

      INTEGER :: error, numDays, pcfId, returnStatus

! Read the PCF into an annotation for file headers

      CALL CreatePCFAnnotation(mlspcf_pcf_start, anText)

! Retrieve values set in PCF & assign them to variables.

      CALL GetPCFParameters(l3pcf)

! Calculate the average orbit period for each day in the input window

      CALL AvgOrbPeriod(l3pcf%l2StartDay, l3pcf%l2EndDay, avgPer)

! Calculate the dates, number of days for the processing window

      CALL SetProcessingWindow(l3pcf%l3StartDay, l3pcf%l3EndDay, dates, &
                               numDays)

! Open the configuration file

      returnStatus = pgs_io_gen_openF (mlspcf_l3cf_start, &
                                       PGSd_IO_Gen_RSeqFrm, 0, pcfId, 1)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = MLSMSG_Fileopen // 'l3cf.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Read the configuration file into an MLSCF_T structure

      CALL InitGetCF

      CALL getCF(cf, error, inUnit=pcfId)

! Close the configuration file

      returnStatus = pgs_io_gen_closeF (pcfId)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                          ModuleName, 'Failed to close l3cf.')

! Check the parser output; fill L3CFProd_T

      CALL FillL3CF(cf, l3pcf%outputVersion, dates, numDays, l3cf, cfDef)

!----------------------------------
   END SUBROUTINE OpenAndInitialize
!----------------------------------

!----------------------------------------------------
   SUBROUTINE AvgOrbPeriod (startDay, endDay, avgPer)
!----------------------------------------------------

! Brief description of subroutine
! This program calculates the average orbital period for each day in the input
! window.

! Arguments

      CHARACTER (LEN=8), INTENT(IN) :: endDay, startDay

      REAL(r8), POINTER :: avgPer(:)

! Parameters

      INTEGER, PARAMETER :: numOffset = 43200

! Functions

      INTEGER, EXTERNAL :: Pgs_eph_getEphMet, Pgs_td_asciiTime_bToA
      INTEGER, EXTERNAL :: Pgs_td_timeInterval, Pgs_td_utcToTAI

! Variables

      CHARACTER (LEN=26) :: startB
      CHARACTER (LEN=CCSDS_LEN) :: startA
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=8) :: dates(maxWindow)
      CHARACTER (LEN=CCSDS_LEN) :: orbitAscendTime(max_orbits)
      CHARACTER (LEN=CCSDS_LEN) :: orbitDescendTime(max_orbits)

      INTEGER :: err, i, j, numDays, numOrbits, returnStatus
      INTEGER :: orbitNumber(max_orbits)

      REAL(r8) :: sum
      REAL(r8) :: offsets(numOffset)
      REAL(r8) :: ascTAI(max_orbits), int(max_orbits)
      REAL(r8) :: orbitDownLongitude(max_orbits)

! Calculate the offsets

      DO i = 1, numOffset
         offsets(i) = (i-1) * 2.0
      ENDDO

! Find the number of days in the input window

      CALL SetProcessingWindow(startDay, endDay, dates, numDays)

      ALLOCATE(avgPer(numDays), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' array of period averages.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! For each day in the input window,

      DO i = 1, numDays

! Convert the start time to CCSDS A format

         startB = dates(i) // 'T00:00:00.000000Z'

         returnStatus = Pgs_td_asciiTime_bToA(startB, startA)
         IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                         ModuleName, 'Failed to convert input date to CCSDS A')

! Get orbit metadata

         returnStatus = Pgs_eph_getEphMet(spacecraftId, numOffset, startA, &
                           offsets,  numOrbits, orbitNumber, orbitAscendTime, &
                           orbitDescendTime, orbitDownLongitude)
         IF (returnStatus /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Convert orbitAscendTime to TAI93

         DO j = 1, numOrbits
            returnStatus = Pgs_td_utcToTAI( orbitAscendTime(j), ascTAI(j) )
            IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                            ModuleName, 'Error converting from CCSDSB to TAI.')
         ENDDO

! Get the differences between successive times and add them for the day

         sum = 0.0 

         DO j = 1, numOrbits-1
            returnStatus = Pgs_td_timeInterval( ascTAI(j), ascTAI(j+1), int(j) )
            IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
              ModuleName, 'Failed to find interval between orbitsAscendTimes.')
            sum = sum + int(j)
         ENDDO

! Get the average period for the day

         avgPer(i) = sum/(numOrbits-1)

      ENDDO

!-----------------------------
   END SUBROUTINE AvgOrbPeriod
!-----------------------------

!==================
END MODULE OpenInit
!==================

! $Log$
! Revision 1.7  2001/02/21 21:11:59  nakamura
! Changed MLSPCF to MLSPCF3; removed inputVersion.
!
! Revision 1.6  2001/01/18 16:52:28  nakamura
! Added type L3CFDef_T; moved minDays from PCF to cf.
!
! Revision 1.5  2001/01/16 17:48:09  nakamura
! Added code for annotating and for retrieving a name template for the log in the cf.
!
! Revision 1.4  2000/12/29 21:42:27  nakamura
! Added subroutine AvgOrbPeriod.
!
! Revision 1.3  2000/11/15 21:28:05  nakamura
! Changed MLSPCF parameter names to mlspcf_l3_param_*.
!
! Revision 1.2  2000/10/24 19:31:48  nakamura
! Added logGranID to PCFData_T; replaced ReadParseMLSCF stub with calls to actual parser.
!
! Revision 1.1  2000/10/17 20:26:16  nakamura
! Module for the Open/Init task.
!
