
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE mon_Out
!===============================================================================

   USE L3DZData
   USE L3MMData
   USE L3MZData
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE mon_Open
   USE PCFModule
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- CreateFlags_T
! Subroutines -- WriteMetaLogM
!                OutputMON

! Remarks:  This is a prototype module for the routines needed for the L3
! Monthly Output/Close task.

! Parameters

! This data type is used to store the flags indicating whether output files have
! been created

   TYPE CreateFlags_T

     LOGICAL :: createMS, createMD
        ! monthly map files, standard & diagnostic

     LOGICAL :: createZS, createZD
        ! monthly zonal mean files, standard & diagnostic

   END TYPE CreateFlags_T

CONTAINS

!----------------------------------------------------
   SUBROUTINE WriteMetaLogM (pcf, startDate, endDate)
!----------------------------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: endDate, startDate

! Parameters

      INTEGER, PARAMETER :: ASCII_FILE = 101

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_s, pgs_met_write

! Variables

      CHARACTER (LEN=1) :: nullStr
      CHARACTER (LEN=45) :: sval
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: result

      nullStr = ''

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf_l3log_start, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                          'Initialization error.  See LogStatus for details.')

! Set PGE values

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalGranuleID", &
                                 pcf%logGranID)

      CALL ExpandFileTemplate('$cycle', sval, cycle=pcf%cycle)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalVersionID", &
                                 sval)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningDate", startDate)
      sval= '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningTime", sval)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", endDate)
      sval= '23:59:59.999999'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingTime", sval)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "PGEVersion", &
                                 pcf%outputVersion)

      IF (result /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(result, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the metadata and their values to an ASCII file

      result = pgs_met_write(groups(1), nullStr, ASCII_FILE)

      IF (result /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(result, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      result = pgs_met_remove()

!------------------------------
   END SUBROUTINE WriteMetaLogM
!------------------------------

!-----------------------------------------------------------------------
   SUBROUTINE OutputMON (sFiles, dFiles, flags, pcf, cfProd, cf, anText)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the monthly Output/Close task for a product.

! Arguments

      TYPE( CreateFlags_T ), INTENT(IN) :: flags

      TYPE( OutputFiles_T ), INTENT(IN) :: dFiles, sFiles

      TYPE( PCFMData_T), INTENT(IN) :: pcf

      TYPE( Mlscf_T ), INTENT(INOUT) :: cf

      TYPE( L3CFMProd_T ), POINTER :: cfProd(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Write the metadata to any L3DZ Standard files created

      IF (sFiles%nFiles == 0) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'No L3DZ Standard files &
                                                               &were created.')
      ELSE
         CALL WriteMetaL3DZ(pcf, mlspcf_mcf_l3dzs_start, sFiles, anText)
      ENDIF

! L3DZ Diagnostic files

      IF (dFiles%nFiles == 0) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'No L3DZ Diagnostic &
                                                     &files were created.')
      ELSE
         CALL WriteMetaL3DZ(pcf, mlspcf_mcf_l3dzd_start, dFiles, anText)
      ENDIF

! If a Standard L3MZ file was created, write the metadata & annotation

      IF (flags%createZS) THEN
         CALL WriteMetaL3MZ (pcf%zsName, mlspcf_mcf_l3mzs_start, pcf, anText)
      ELSE
         CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                         'No monthly zonal mean std file was produced.')
      ENDIF

! L3MZ Dg file

      IF (flags%createZD) THEN
         CALL WriteMetaL3MZ (pcf%zdName, mlspcf_mcf_l3mzd_start, pcf, anText)
      ELSE
         CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                        'No monthly zonal mean dg file was produced.')
      ENDIF

! L3MM Std file

      IF (flags%createMS) THEN
         CALL WriteMetaL3MM (pcf%msName, mlspcf_mcf_l3mms_start, pcf, anText)
      ELSE
         CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                         'No monthly map std file was produced.')
      ENDIF

! L3MM Dg file

      IF (flags%createMD) THEN
         CALL WriteMetaL3MM (pcf%mdName, mlspcf_mcf_l3mmd_start, pcf, anText)
      ELSE
         CALL MLSMessage(MLSMSG_Warning, ModuleName, &
                         'No monthly map dg file was produced.')
      ENDIF

! Write the log file metadata

      CALL WriteMetaLogM(pcf, pcf%startDay, pcf%endDay)

! Deallocations

      DEALLOCATE(cfProd, anText, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  Open/Init quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE (cf%Sections, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  cf section pointers.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!--------------------------
   END SUBROUTINE OutputMON
!--------------------------

!=================
END MODULE mon_Out
!=================

!$Log$
!
