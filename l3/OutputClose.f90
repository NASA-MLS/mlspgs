
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE OutputClose
!===============================================================================

   USE MLSCommon
   USE SDPToolkit
   USE MLSMessageModule
   USE L3DMData
   USE L3CF
   USE MLSStrings
   USE MLSPCF
   USE OpenInit
   USE MLSCF
   USE PCFModule
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- CheckOutputDate
!                FindDatabaseIndex
!                WriteMetaLog
!                OutputAndClose

! Remarks:  This is a prototype module for the routines needed for the L3 Daily
! Output/Close task.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: NOOUT_ERR = ' data expected but not found &
                                               &for output.'

CONTAINS

!---------------------------------------------
   SUBROUTINE CheckOutputDate (dataTAI, timeA)
!---------------------------------------------

! Brief description of subroutine
! This subroutine checks the dates within a database for consistency and
! returns the last date/time in CCSDS A format.

! Arguments

      REAL(r8), INTENT(IN) :: dataTAI(:)

      CHARACTER (LEN=CCSDS_LEN), INTENT(OUT) :: timeA

! Parameters

! Functions

! Variables

      CHARACTER (LEN=CCSDS_LEN) :: cmpTime

      INTEGER :: i, returnStatus

! Check that dates within the l3dm database are consistent

      returnStatus = Pgs_td_taiToUTC(dataTAI(1), cmpTime)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                                         ModuleName, TAI2A_ERR)

      DO i = 2, SIZE(dataTAI)

         returnStatus = Pgs_td_taiToUTC(dataTAI(i), timeA)
         IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                                     ModuleName, TAI2A_ERR)

         IF (timeA(1:10) /= cmpTime(1:10)) CALL MLSMessage(MLSMSG_Error, &
                  ModuleName,'Inconsistent dates within the output database.')

      ENDDO

!--------------------------------
   END SUBROUTINE CheckOutputDate
!--------------------------------

!-----------------------------------------------------------
   SUBROUTINE FindDatabaseIndex (l3dm, l3cf, indx, numGrids)
!-----------------------------------------------------------

! Brief description of subroutine
! This subroutine checks for the existence of data in the database and returns
! the number of structures found, and their database indices.

! Arguments

      TYPE( L3DMData_T ), INTENT(IN) :: l3dm(:)

      TYPE( L3CFProd_T ), INTENT(IN) :: l3cf

      INTEGER, INTENT(OUT) :: indx(:)

      INTEGER, INTENT(OUT) :: numGrids

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: i, idm

! Loop through all possible grids for the product

      DO i = 1, l3cf%nGrids

! Check for a corresponding name in the l3dm database

         idm = LinearSearchStringArray( l3dm%name,TRIM(l3cf%quantities(i)) )

! If some grid is missing, check the processing mode requested, and issue a
! message, if necessary

         IF (idm == 0) THEN
 
            IF ( INDEX(l3cf%quantities(i),'Ascending') /= 0 ) THEN

               IF ( (l3cf%mode == 'asc') .OR. (l3cf%mode == 'all') ) THEN
                  msr = 'Ascending ' // TRIM(l3cf%l3prodNameD) // NOOUT_ERR
                  CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
               ENDIF

            ELSE IF ( INDEX(l3cf%quantities(i),'Descending') /= 0 ) THEN

                 IF ( (l3cf%mode == 'des') .OR. (l3cf%mode == 'all') ) THEN
                    msr = 'Descending ' // TRIM(l3cf%l3prodNameD) // NOOUT_ERR
                    CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
                 ENDIF

            ELSE IF ( l3cf%quantities(i) == l3cf%l3prodNameD ) THEN

                 IF ( (l3cf%mode == 'com') .OR. (l3cf%mode == 'all') ) THEN
                    msr = TRIM(l3cf%l3prodNameD) // NOOUT_ERR
                    CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
                 ENDIF

            ENDIF

         ELSE

! Calculate the number of grids found; save their indices

            numGrids = numGrids + 1
            indx(numGrids) = idm

         ENDIF

      ENDDO

!----------------------------------
   END SUBROUTINE FindDatabaseIndex
!----------------------------------

!--------------------------------
   SUBROUTINE WriteMetaLog (date)
!--------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: date

! Parameters

      INTEGER, PARAMETER :: ASCII_FILE = 101

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_s, pgs_met_write

! Variables

      CHARACTER (LEN=1) :: nullStr
      CHARACTER (LEN=21) :: sval
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: result

      REAL(r8) :: dval

      nullStr = ''

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf_l3log_start, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                          'Initialization error.  See LogStatus for details.')

! Set PGE values

      sval = 'Horizontal & Vertical'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "GranuleSpatialDomainType", sval)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningDate", date)
      sval= '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningTime", sval)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", date)
      sval= '23:59:59.999999'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingTime", sval)

      sval = 'Other Grid System'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "ZoneIdentifier", &
                                 sval)

      dval = -180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "WestBoundingCoordinate", dval)
      dval = 90.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "NorthBoundingCoordinate", dval)
      dval = 180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "EastBoundingCoordinate", dval)
      dval = -90.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "SouthBoundingCoordinate", dval)

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

!-------------------------------------------------
   SUBROUTINE OutputAndClose (pcf, cf, l3cf, l3dm)
!-------------------------------------------------

! Brief description of subroutine
! This subroutine performs the Output/Close task in the MLSL3 program.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      TYPE( Mlscf_T ), INTENT(INOUT) :: cf

      TYPE( L3CFProd_T ), POINTER :: l3cf(:)

      TYPE( L3DMData_T ), POINTER :: l3dm(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: Pgs_td_asciiTime_aToB

! Variables

      CHARACTER (LEN=8) :: procDay
      CHARACTER (LEN=CCSDS_LEN) :: dataDT
      CHARACTER (LEN=CCSDSB_LEN) :: timeB
      CHARACTER (LEN=FileNameLen) :: l3File
      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i, match, mlspcf_l3dm, mlspcf_mcf, numGrids, numProds
      INTEGER :: returnStatus
      INTEGER :: indx(maxNumGrids)

! Check that dates within the l3dm database are consistent

      CALL CheckOutputDate(l3dm%time, dataDT)

! Convert to CCSDS B format, for comparison to file name

      returnStatus = Pgs_td_asciiTime_aToB(dataDT,timeB)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                 ModuleName,'Error converting data time from CCSDS A to B.')

      procDay = timeB(1:8)

! Check that date falls within the output data range

      IF ( LLT(procDay,pcf%l3StartDay) .OR. LGT(procDay,pcf%l3EndDay) ) THEN
         msr = 'Data produced for day ' // procDay // ' is not in the &
                                                      &specified output range.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! For each product in the DailyMap section of the l3cf,

      numProds = SIZE(l3cf)

      DO i = 1, numProds

         numGrids = 0
         indx = 0

         CALL FindDatabaseIndex(l3dm, l3cf(i), indx, numGrids)

! If all grids are missing, issue a warning for the product

         IF (numGrids == 0) THEN
            msr = 'No l3dm data found, no file produced for ' // &
                  l3cf(i)%l3prodNameD
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            CYCLE
         ENDIF

! Expand any fields in the given output file name

         CALL ExpandFileTemplate(l3cf(i)%fileTemplate, l3File, &
                                 pcf%outputVersion, pcf%cycle, procDay)

! If the bypass flag is set, issue a message to that effect, giving the file
! name to be used

         IF (l3cf(i)%bpFlag == 1) THEN

            msr = 'Bypassing PCF:  using file name ' // l3File
            CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

         ELSE

! If not, check that the expanded name appears in the PCF.  Exit with an error,
! if the required name has no PCF entry

            CALL SearchPCFNames(l3File, mlspcf_l3dm_start, mlspcf_l3dm_end, &
                                mlspcf_l3dm, match)
            IF (match == 0) THEN
               msr = 'No match in the PCF for file ' // l3File
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
     
         ENDIF

! Create & write to structures in the file

         CALL OutputGrids(l3File, numGrids, indx, l3dm)

! Find the PCF number of the MCF

         CALL SearchPCFNames(l3cf(i)%mcf, mlspcf_l3dm_mcf_start, &
                             mlspcf_l3dm_mcf_end, mlspcf_mcf, match)
         IF (match == 0) THEN
            msr = 'No match in the PCF for file ' // l3cf(i)%mcf
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the L3DM metadata

         CALL WriteMetaL3DM(l3File, mlspcf_mcf, numGrids, indx, l3dm, dataDT)

      ENDDO

! Write the log file metadata

      CALL WriteMetaLog(dataDT(1:10))

! Deallocate the l3dm database

      CALL DestroyL3DMDatabase(l3dm)

! Deallocate the l3cf pointer

      DEALLOCATE(l3cf, STAT=err)
      IF ( err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                               &deallocation of l3dm input data structures.')

! Deallocate the section pointer of the CF structure

      DEALLOCATE (cf%Sections, STAT=err)
      IF ( err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                         &deallocate l3cf section pointers.')

!-------------------------------
   END SUBROUTINE OutputAndClose
!-------------------------------

!=====================
END MODULE OutputClose
!=====================

!$Log$
