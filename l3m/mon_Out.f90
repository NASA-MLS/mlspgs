
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
!                OutputStd
!                OutputDg
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

!--------------------------------
   SUBROUTINE WriteMetaLogM (pcf)
!--------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

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
                                 "RangeBeginningDate", pcf%startDay)
      sval= '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningTime", sval)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", pcf%endDay)
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

!---------------------------------------------------------------------------
   SUBROUTINE OutputStd(pcf, type, mode, dzA, dzD, mzA, mzD, mm, mmA, mmD, &
                        sFiles, flag)
!---------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Output/Close task within the standard product loop of
! the L3 Monthly subprogram.

! Arguments


      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: mode, type

      TYPE( L3DZData_T ), POINTER :: dzA(:), dzD(:)

      TYPE( L3MMData_T ), INTENT(INOUT) :: mm, mmA, mmD

      TYPE( L3MZData_T ), INTENT(INOUT) :: mzA, mzD

      TYPE( OutputFiles_T ), INTENT(OUT) :: sFiles

      TYPE( CreateFlags_T ), INTENT(OUT) :: flag

! Parameters

! Functions

! Variables

! Daily Zonal Mean output

      CALL OutputL3DZ(type, dzA, sFiles)
      CALL DestroyL3DZDatabase(dzA)

      CALL OutputL3DZ(type, dzD, sFiles)
      CALL DestroyL3DZDatabase(dzD)

! Monthly Zonal Mean output

      CALL OutputL3MZ(pcf%zsName, mzA, flag%createZS)
      CALL DeallocateL3MZ(mzA)

      CALL OutputL3MZ(pcf%zsName, mzD, flag%createZS)
      CALL DeallocateL3MZ(mzD)

! If required for this mode, output the monthly map

      IF ( (mode == 'com') .OR. (mode == 'all') ) THEN
         CALL OutputMMGrids(pcf%msName, mm, flag%createMS)
         CALL OutputMMDiags(pcf%msName, mm)
      ENDIF
      CALL DeallocateL3MM(mm)

! Ascending

      IF ( INDEX(mode,'a') /= 0) THEN
         CALL OutputMMGrids(pcf%msName, mmA, flag%createMS)
         CALL OutputMMDiags(pcf%msName, mmA)
      ENDIF
      CALL DeallocateL3MM(mmA)

! Descending

      IF ( (INDEX(mode,'d') /= 0) .OR. (mode == 'all') )THEN
         CALL OutputMMGrids(pcf%msName, mmD, flag%createMS)
         CALL OutputMMDiags(pcf%msName, mmD)
      ENDIF
      CALL DeallocateL3MM(mmD)

!--------------------------
   END SUBROUTINE OutputStd
!--------------------------

!--------------------------------------------------------------------------------
   SUBROUTINE OutputDg(pcf, type, dzA, dzD, mzA, mzD, mm, mmA, mmD, dFiles, flag)
!--------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Output/Close task within the diagnostic product loop of
! the L3 Monthly subprogram.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE( L3DZData_T ), POINTER :: dzA(:), dzD(:)

      TYPE( L3MMData_T ), INTENT(INOUT) :: mm, mmA, mmD

      TYPE( L3MZData_T ), INTENT(INOUT) :: mzA, mzD

      TYPE( OutputFiles_T ), INTENT(OUT) :: dFiles

      TYPE( CreateFlags_T ), INTENT(OUT) :: flag

! Parameters

! Functions

! Variables

! Deallocate unused databases

      CALL DeallocateL3MM(mmA)
      CALL DeallocateL3MM(mmD)

! Daily Zonal Mean output

      CALL OutputL3DZ(type, dzA, dFiles)
      CALL DestroyL3DZDatabase(dzA)

      CALL OutputL3DZ(type, dzD, dFiles)
      CALL DestroyL3DZDatabase(dzD)

! Monthly Zonal Mean output

      CALL OutputL3MZ(pcf%zdName, mzA, flag%createZD)
      CALL DeallocateL3MZ(mzA)

      CALL OutputL3MZ(pcf%zdName, mzD, flag%createZD)
      CALL DeallocateL3MZ(mzD)

! Output the monthly map (combined mode)

      CALL OutputMMGrids(pcf%mdName, mm, flag%createMD)
      CALL OutputMMDiags(pcf%mdName, mm)
      CALL DeallocateL3MM(mm)

!-------------------------
   END SUBROUTINE OutputDg
!-------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE OutputMON (sFiles, dFiles, flags, pcf, cfProd, cfDg, cf, anText)
!-----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the monthly Output/Close task for a product.

! Arguments

      TYPE( CreateFlags_T ), INTENT(IN) :: flags

      TYPE( OutputFiles_T ), INTENT(IN) :: dFiles, sFiles

      TYPE( PCFMData_T), INTENT(IN) :: pcf

      TYPE( Mlscf_T ), INTENT(INOUT) :: cf

      TYPE( L3CFMProd_T ), POINTER :: cfDg(:), cfProd(:)

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

      CALL WriteMetaLogM(pcf)

! Deallocations

      DEALLOCATE(cfProd, cfDg, anText, STAT=err)
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
!Revision 1.4  2001/09/26 19:49:13  nakamura
!Removed com ZM output; added cfDg deallocate.
!
!Revision 1.3  2001/09/06 18:51:45  nakamura
!Added subroutine OutputDg; moved database deallocation back down into Output subroutines.
!
!Revision 1.2  2001/08/01 18:30:03  nakamura
!Added OutputStd subroutine; updated WriteMetaLogM for separation from Daily.
!
!Revision 1.1  2001/07/18 15:44:27  nakamura
!Module for the Monthly Output/Close task.
!
