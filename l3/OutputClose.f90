
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE OutputClose
!===============================================================================

   USE L2GPData
   USE L2Interface
   USE L3CF
   USE L3DMData
   USE L3SPData
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF
   USE MLSStrings
   USE OpenInit
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

! Definitions -- OutputFlags_T
! Subroutines -- WriteMetaLog
!                OutputAndClose

! Remarks:  This is a prototype module for the routines needed for the L3 Daily
! Output/Close task.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: NOOUT_ERR = ' data expected but not found &
                                                &for output.'

! This data type is used to store the flags indicating whether the product
! databases are suitable for output

   TYPE OutputFlags_T

     LOGICAL :: writel3dmCom, writel3dmAsc, writel3dmDes
	! daily map databases (for all output days)

     LOGICAL :: writel3rCom, writel3rAsc, writel3rDes
	! L3Residual databases (for all output days)

     LOGICAL :: writel3sp
	! L3SP database (for asc/des/com)

   END TYPE OutputFlags_T

CONTAINS

!----------------------------------------
   SUBROUTINE WriteMetaLog (pcf, logType)
!----------------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: logType

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

      sval = TRIM(logType) // pcf%l3StartDay // ':' // pcf%l3EndDay
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalGranuleID", &
                                 sval)

      CALL ExpandFileTemplate('$version-$cycle', sval, &
                              version=pcf%outputVersion, cycle=pcf%cycle)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalVersionID", &
                                 sval)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningDate", pcf%l3StartDay)
      sval= '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningTime", sval)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", pcf%l3EndDay)
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

!-----------------------------
   END SUBROUTINE WriteMetaLog
!-----------------------------

!----------------------------------------------------------------------------
   SUBROUTINE OutputAndClose (pcf, l3cf, anText, l3sp, l3dm, dmA, dmD, l3r, &
                              residA, residD, flags)
!----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Output/Close task in the MLSL3 program.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      TYPE( L3CFProd_T ), INTENT(IN) :: l3cf

      TYPE( OutputFlags_T ), INTENT(IN) :: flags

      TYPE( L2GPData_T ), POINTER :: l3r(:), residA(:), residD(:)

      TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)

      TYPE( L3SPData_T ), POINTER :: l3sp(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

! Variables

      TYPE( L3DMFiles_T ) :: files

      CHARACTER (LEN=FileNameLen) :: type
      CHARACTER (LEN=480) :: msr

! L3SP -- if data exist, create & write a file for this product

      IF (flags%writel3sp) THEN
         CALL OutputL3SP (l3cf, anText, l3sp)
      ELSE
         msr = TRIM(l3cf%l3prodNameD) // ' L3SP' // NOOUT_ERR
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      ENDIF

! L3DM -- check that daily map data exist in some form for the product

      IF ( .NOT.(flags%writel3dmCom) .AND. .NOT.(flags%writel3dmAsc) .AND. &
           .NOT.(flags%writel3dmDes) ) THEN

         msr = 'No l3dm data found, no file produced for ' // l3cf%l3prodNameD
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)

      ELSE

! Expand the level in the file template

         CALL ExpandFileTemplate(l3cf%fileTemplate, type, 'L3DM')

! If combined l3dm data exist, write them to l3dm files; save the names of any
! files created

         files%nFiles = 0
         files%name = ''

         IF (flags%writel3dmCom) THEN
            CALL OutputGrids(type, l3dm, files)
         ELSE
            IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'com') ) THEN
               msr = TRIM(l3cf%l3prodNameD) // NOOUT_ERR
               CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            ENDIF
         ENDIF

! If ascending l3dm data exist, write/append them to l3dm files; add the names
! of any new files created to the files array

         IF (flags%writel3dmAsc) THEN
            CALL OutputGrids(type, dmA, files)
         ELSE
            IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'asc') ) THEN
               msr = TRIM(l3cf%l3prodNameD) // 'Ascending' // NOOUT_ERR
               CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            ENDIF
         ENDIF

! If descending l3dm data exist, write/append them to l3dm files, adding any
! new names to the files array

         IF (flags%writel3dmDes) THEN
            CALL OutputGrids(type, dmD, files)
         ELSE
            IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'des') ) THEN
               msr = TRIM(l3cf%l3prodNameD) // 'Descending' // NOOUT_ERR
               CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            ENDIF
         ENDIF

! L3Residual -- if any residual data exist, append swaths to the l3dm files

         IF (flags%writel3rCom) THEN
            CALL ResidualOutput(type, l3r)
         ELSE
            IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'com') ) THEN
               msr = TRIM(l3cf%l3prodNameD) // 'Residuals' // NOOUT_ERR
               CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            ENDIF
         ENDIF

         IF (flags%writel3rAsc) THEN
            CALL ResidualOutput(type, residA)
         ELSE
            IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'asc') ) THEN
              msr = TRIM(l3cf%l3prodNameD) // 'AscendingResiduals' // NOOUT_ERR
              CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            ENDIF
         ENDIF

         IF (flags%writel3rDes) THEN
           CALL ResidualOutput(type, residD)
         ELSE
           IF ( (l3cf%mode == 'all') .OR. (l3cf%mode == 'des') ) THEN
             msr = TRIM(l3cf%l3prodNameD) // 'DescendingResiduals' // NOOUT_ERR
             CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
           ENDIF
         ENDIF

! Write the L3DM metadata

         CALL WriteMetaL3DM(pcf, l3cf, files, anText)

      ENDIF

! Deallocate the databases

      CALL DestroyL2GPDatabase(l3r)
      CALL DestroyL2GPDatabase(residA)
      CALL DestroyL2GPDatabase(residD)

      CALL DestroyL3DMDatabase(l3dm)
      CALL DestroyL3DMDatabase(dmA)
      CALL DestroyL3DMDatabase(dmD)

      CALL DestroyL3SPDatabase(l3sp)

!-------------------------------
   END SUBROUTINE OutputAndClose
!-------------------------------

!=====================
END MODULE OutputClose
!=====================

!$Log$
!Revision 1.5  2000/12/29 21:44:16  nakamura
!Removed obsolete routines CheckOutputDate & FindDatabaseIndex; switched to one-product/all-days paradigm.
!
!Revision 1.4  2000/12/07 21:20:43  nakamura
!Changed search/bypass PCF logic, so that the SearchPCFNames subroutine always returns a file name with a path.
!
!Revision 1.3  2000/12/07 19:40:37  nakamura
!Updated for modified ExpandFileTemplate.
!
!Revision 1.2  2000/10/24 19:36:56  nakamura
!Updated WriteMetaLog for new MCF; moved search for PCF number of MCF to L3CF module.
!
!Revision 1.1  2000/10/17 20:27:49  nakamura
!Module for the Output/Close task.
!
