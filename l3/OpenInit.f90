
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE OpenInit 
!===============================================================================

   USE MLSCommon
   USE SDPToolkit
   USE MLSPCF
   USE MLSMessageModule
   USE L3DMData
   USE dates_module
   USE MLSCF
   USE L3CF
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

! Remarks:  This is a prototype module for the routines needed for the L3 Daily
! Open/Init task.

! Parameters

! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFData_T

     ! cycle # of processing run

     CHARACTER (LEN=4) :: cycle

     ! version string in PCF file names

     CHARACTER (LEN=15) :: inputVersion		! input files
     CHARACTER (LEN=15) :: outputVersion	! output files

     ! first & last days of input & output windows

     CHARACTER (LEN=8) :: l2EndDay		! last day of l2 to read
     CHARACTER (LEN=8) :: l2StartDay		! first day of l2 to read
     CHARACTER (LEN=8) :: l3EndDay		! last day of l3 to write
     CHARACTER (LEN=8) :: l3StartDay		! first day of l3 to write

     ! # of days of input data needed to process an l3 product

     INTEGER :: minDays

     ! granule ID of the log file (file name, without path or .dat)

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

      CHARACTER (LEN=*), PARAMETER :: UDRP_ERR = 'Failed to get runtime &
                                              &parameter from PCF:  '

! Functions

      INTEGER, EXTERNAL :: pgs_pc_getConfigData

! Variables

      CHARACTER (LEN=2) :: aDays
      CHARACTER (LEN=17) :: range 
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: name

      INTEGER ::  i, indx, mlspcf_log, returnStatus, version

! Retrieve values set in PCF & assign them to variables.

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_InputVersion, &
                                          l3pcf%inputVersion)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'InputVersion'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

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

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_MinDays, aDays)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'MinDays'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      READ(aDays, '(I2)') l3pcf%minDays

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

      DO i = 1, 10
         indx = INDEX(name, '/')
         IF (indx == 0) EXIT
         name = name(indx+1:)
      ENDDO

      indx = INDEX(name, '.')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No file type &
                                    &specified in log name.')

      l3pcf%logGranID = name(:indx-1)

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

      DO i = 2, maxWindow

! Increment DOY

         iDOY = iDOY + 1

! Check that the incrementation hasn't crossed a year boundary

         IF (iDOY > lastDOY) THEN
            iYr = iYr + 1
            iDOY = 1
         ENDIF

! Check that the incrementation hasn't crossed the window boundary

         IF (iDOY >= iEnd) EXIT

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

        dates(i) = cYr // '-' // cDOY

      ENDDO

      numDays = i
      dates(1) = start
      dates(numDays) = end

!------------------------------------
   END SUBROUTINE SetProcessingWindow
!------------------------------------

!------------------------------------------------
   SUBROUTINE OpenAndInitialize (l3pcf, cf, l3cf)
!------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Open/Init task in the MLSL3 program.

! Arguments

      TYPE( PCFData_T ), INTENT(OUT) :: l3pcf

      TYPE( Mlscf_T ), INTENT(OUT) :: cf

      TYPE( L3CFProd_T ), POINTER :: l3cf(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=8) :: dates(maxWindow)

      INTEGER :: error, numDays, pcfId, returnStatus

! Retrieve values set in PCF & assign them to variables.

      CALL GetPCFParameters(l3pcf)

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

      CALL FillL3CF(cf, l3pcf%inputVersion, l3pcf%outputVersion, dates, &
                    numDays, l3cf)

!----------------------------------
   END Subroutine OpenAndInitialize
!----------------------------------

!==================
END MODULE OpenInit
!==================

! $Log$
! Revision 1.2  2000/10/24 19:31:48  nakamura
! Added logGranID to PCFData_T; replaced ReadParseMLSCF stub with calls to actual parser.
!
! Revision 1.1  2000/10/17 20:26:16  nakamura
! Module for the Open/Init task.
!
