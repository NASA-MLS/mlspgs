
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE mon_Open
!===============================================================================

   USE MLSCF
   USE MLSCommon
   USE MLSMessageModule
   USE MLSPCF3
   USE mon_L3CF
   USE OpenInit, ONLY:  SetProcessingWindow
   USE PCFHdr
   USE SDPToolkit
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

! Definitions -- PCFMData_T
! Subroutines -- GetPCFParams
!                OpenMON

! Remarks:  This is a prototype module for the routines needed for the L3
! Monthly Open/Init task.

! Parameters

! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFMData_T

     ! cycle # of processing run

     CHARACTER (LEN=4) :: cycle

     ! version string in PCF output file names

     CHARACTER (LEN=15) :: outputVersion

     ! processing window (input & output) information

     CHARACTER (LEN=8) :: endDay

     CHARACTER (LEN=8) :: startDay

     INTEGER :: nDays

     CHARACTER (LEN=8) :: dates(maxWindow)

     ! log file name, WITHOUT the path (LocalGranuleID)

     CHARACTER (LEN=FileNameLen) :: logGranID

     ! output file names, WITH their paths

     CHARACTER (LEN=FileNameLen) :: msName	! L3MM_Standard

     CHARACTER (LEN=FileNameLen) :: mdName	! L3MM_Diagnostic

     CHARACTER (LEN=FileNameLen) :: zsName	! L3MZ_Standard

     CHARACTER (LEN=FileNameLen) :: zdName	! L3MZ_Diagnostic

   END TYPE PCFMData_T

CONTAINS

!---------------------------------
   SUBROUTINE GetPCFParams (l3pcf)
!---------------------------------

! Brief description of subroutine
! This subroutine retrieves the User-Defined Runtime Parameters from the PCF
! and stores them in variables used in the MLSL3 code.

! Arguments

      TYPE( PCFMData_T ), INTENT(OUT) :: l3pcf

! Parameters

      CHARACTER (LEN=*), PARAMETER :: UDRP_ERR = 'Failed to get PCF runtime &
                                                                &parameter:  '

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

! Set all information on the processing window

      returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_L2DayRange, range)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         msr = UDRP_ERR // 'L2DayRange'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      l3pcf%startDay = range(1:8)
      l3pcf%endDay = range(10:17)

      CALL SetProcessingWindow(l3pcf%startDay, l3pcf%endDay, l3pcf%dates, &
                               l3pcf%nDays)

! Get the name of the log file from the PCF; remove its path

      version = 1
      mlspcf_log = 10101

      returnStatus = Pgs_pc_getReference(mlspcf_log, version, name)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                        ModuleName, 'Error retrieving log file name from PCF.')
      indx = INDEX(name, '/', .TRUE.)
      l3pcf%logGranID = name(indx+1:)

! Get the name of the L3MM Standard file from the PCF

      version = 1
      returnStatus = Pgs_pc_getReference(mlspcf_l3mms_start, version, l3pcf%msName)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                 ModuleName, 'Error retrieving L3MM std file name from the PCF.')

! Get the name of the L3MM Diagnostic file from the PCF

      version = 1
      returnStatus = Pgs_pc_getReference(mlspcf_l3mmd_start, version, l3pcf%mdName)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                 ModuleName, 'Error retrieving L3MM dg file name from the PCF.')

! Get the name of the L3MZ Standard file from the PCF

      version = 1
      returnStatus = Pgs_pc_getReference(mlspcf_l3mzs_start, version, l3pcf%zsName)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                 ModuleName, 'Error retrieving L3MZ std file name from the PCF.')

! Get the name of the L3MZ Diagnostic file from the PCF

      version = 1
      returnStatus = Pgs_pc_getReference(mlspcf_l3mzd_start, version, l3pcf%zdName)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                 ModuleName, 'Error retrieving L3MM dg file name from the PCF.')

!-----------------------------
   END SUBROUTINE GetPCFParams
!-----------------------------

!-----------------------------------------------------------
   SUBROUTINE OpenMON (l3pcf, cf, l3cf, cfDef, cfDg, anText)
!-----------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Open/Init task in the MLSL3 program.

! Arguments

      TYPE( L3CFMDef_T ), INTENT(OUT) :: cfDef

      TYPE( Mlscf_T ), INTENT(OUT) :: cf

      TYPE( PCFMData_T ), INTENT(OUT) :: l3pcf

      TYPE( L3CFDg_T ), POINTER :: cfDg(:)

      TYPE( L3CFMProd_T ), POINTER :: l3cf(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: error, pcfId, returnStatus

! Read the PCF into an annotation for file headers

      CALL CreatePCFAnnotation(mlspcf_pcf_start, anText)

! Retrieve values set in PCF & assign them to variables.

      CALL GetPCFParams(l3pcf)

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

      CALL FillL3CFM(cf, l3pcf%outputVersion, l3cf, cfDef, cfDg)

!------------------------
   END SUBROUTINE OpenMON
!------------------------

!==================
END MODULE mon_Open
!==================

! $Log$
! Revision 1.1  2001/07/18 15:43:44  nakamura
! Module for the Monthly Open/Init task.
!
!
