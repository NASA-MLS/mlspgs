
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE OutputClose
!===============================================================================

   USE L2GPData, ONLY: L2GPData_T
   USE L2Interface
   USE L3CF
   USE L3DMData
   USE L3SPData
   USE MLSCF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
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
!                OutputProd
!                OutputAndClose

! Remarks:  This is a prototype module for the routines needed for the L3 Daily
! Output/Close task.

! Parameters

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

!-------------------------------
   SUBROUTINE WriteMetaLog (pcf)
!-------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf

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

!-------------------------------------------------------------------
   SUBROUTINE OutputProd (pcf, l3cf, anText, l3sp, l3dm, dmA, dmD, &
                          l3r, residA, residD, flags)
!-------------------------------------------------------------------

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

      TYPE( OutputFiles_T ) :: files

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
         files%date = ''

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

!---------------------------
   END SUBROUTINE OutputProd
!---------------------------

!-------------------------------------------------------------
   SUBROUTINE OutputAndClose (cf, pcf, cfProd, avgPer, anText)
!-------------------------------------------------------------

! Brief description of subroutine
! This subroutine performs final Output & Close tasks outside the product loop.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      TYPE( L3CFProd_T ), POINTER :: cfProd(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

      REAL(r8), POINTER :: avgPer(:)

      TYPE( Mlscf_T ), INTENT(INOUT) :: cf

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Write the log file metadata

      CALL WriteMetaLog(pcf)

! Final deallocations
 
      DEALLOCATE(cfProd, anText, avgPer, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  Open/Init quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE (cf%Sections, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  cf section pointers.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-------------------------------
   END SUBROUTINE OutputAndClose
!-------------------------------

!=====================
END MODULE OutputClose
!=====================

!$Log$
!Revision 1.14  2001/11/26 19:27:14  nakamura
!Added L3DMDiag stuff.
!
!Revision 1.13  2001/07/19 14:03:27  nakamura
!Removed DZ stuff.
!
!Revision 1.12  2001/05/04 18:40:04  nakamura
!Changed to generic OutputFiles_T, WriteMetaL3DZ.
!
!Revision 1.11  2001/04/24 19:41:41  nakamura
!Added ONLY to USE L2GPData statement.
!
!Revision 1.10  2001/04/11 18:51:38  nakamura
!Moved deallocations for pointers passed between CORE & I/O up to main program.
!
!Revision 1.9  2001/03/27 19:36:43  nakamura
!Moved a parameter to MLSL3Common; updated metadata.
!
!Revision 1.8  2001/03/07 21:28:41  nakamura
!Commented out database deallocations.
!
!Revision 1.7  2001/02/21 21:16:31  nakamura
!Changed MLSPCF to MLSPCF3; added L3DZ stuff; changed log LocalGranuleID; renamed/split OutputAndClose tasks.
!
!Revision 1.6  2001/01/16 17:49:03  nakamura
!Updated for new MCFs and annotation.
!
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
