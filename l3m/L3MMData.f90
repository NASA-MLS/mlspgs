
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3MMData
!===============================================================================

   USE Hdf
   USE L3DMData, ONLY: ConvertDeg2DMS
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE mon_Open, ONLY: PCFMData_T
   USE PCFHdr
   USE PCFModule
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L3MMData_T
! Subroutines -- OutputMMGrids
!                OutputMMDiags
!                WriteMetaL3MM
!                AllocateL3MM
!                DeallocateL3MM

! Remarks:  This module contains the definition of the L3MMData type, as well
!           as any routines pertaining to it.

! Parameters

! This data type is used to store the l3 monthly map data and diagnostics.

   TYPE L3MMData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes
     INTEGER :: nLons				! Total number of longitudes
     INTEGER :: nMisDays			! Number of missing days in the month

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)
     REAL(r8), DIMENSION(:), POINTER :: longitude	! dimensioned (nLons)

     REAL(r8) :: startTime	! start time of L2 data used in the analysis
     REAL(r8) :: endTime 	! end time of L2 data used in the analysis

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:,:), POINTER :: l3mmValue	  ! Field value
     REAL(r8), DIMENSION(:,:,:), POINTER :: l3mmPrecision ! Field precision
	! dimensioned as (nLevels, nLats, nLons)

     ! Now the diagnostic fields:

     INTEGER, DIMENSION(:), POINTER :: perMisPoints	! Missing points (percentage)
	! dimensioned (nLevels)

     CHARACTER (LEN=DATE_LEN), DIMENSION(maxWindow) :: misDays
	! Missing days, dimensioned (nDays)

   END TYPE L3MMData_T

CONTAINS

!-----------------------------------------------------------------
   SUBROUTINE OutputMMGrids (physicalFilename, l3mm, creationFlag)
!-----------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the grid portion of l3mm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

      TYPE (L3MMData_T), INTENT(IN) :: l3mm

      LOGICAL, INTENT(INOUT) :: creationFlag
   
! Parameters

      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDMV = 'L3mmValue'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDMP = 'L3mmPrecision'

! Functions

      INTEGER, EXTERNAL :: gdattach, gdclose, gdcreate, gddefdim, gddeffld
      INTEGER, EXTERNAL :: gddefproj, gddetach, gdopen, gdwrfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: gdfID, gdId, status
      INTEGER :: start(3), stride(3), edge(3)

      REAL :: maxLat, maxLon, minLat, minLon

      REAL(r8) :: uplft(2), lowrgt(2)
      REAL(r8) :: projparm(13)

! Open the output file

      gdfID = gdopen(physicalFilename, DFACC_RDWR)
      IF (gdfID == -1) THEN
         msr = MLSMSG_Fileopen // physicalFilename
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Set up the grid.  The region is bounded by 180.0W to 176.0E longitude &
! varying latitude.  Grid into 90 bins along the x-axis, by nLats bins along
! the y-axis (4x2 bins).  Upper Left & Lower Right corners in DDDMMMSSS.ss.

      projparm = 0.0

! Find boundaries of measured latitude, upper bound of measure longitude

      maxLat = MAXVAL( l3mm%latitude )
      minLat = MINVAL( l3mm%latitude )
      maxLon = MAXVAL( l3mm%longitude )
      minLon = MINVAL( l3mm%longitude )

! Convert to "packed degree format"

      CALL ConvertDeg2DMS(maxLat, uplft(2))
      CALL ConvertDeg2DMS(minLat, lowrgt(2))
      CALL ConvertDeg2DMS(maxLon, lowrgt(1))
      CALL ConvertDeg2DMS(minLon, uplft(1))

! Create the grid

      gdId = gdcreate(gdfID, l3mm%name, l3mm%nLons, l3mm%nLats, uplft, lowrgt)
      IF (gdId == -1) THEN
          msr = 'Failed to create grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the dimensions

      status = gddefdim(gdId, DIMX_NAME, l3mm%nLons)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMX_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddefdim(gdId, DIMY_NAME, l3mm%nLats)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMY_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddefdim(gdId, DIMZ_NAME, l3mm%nLevels)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMZ_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddefdim(gdId, DIMT_NAME, MIN_MAX)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMT_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define a "Geographic projection," using defaults in all unneeded fields

      status = gddefproj(gdId, GCTP_GEO, 0, 0, projparm)
      IF (status /= 0) THEN
          msr = 'Failed to define projection for grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "geolocation" fields

      status = gddeffld(gdId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddeffld(gdId, GEO_FIELD9, DIMZ_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddeffld(gdId, GEO_FIELD1, DIMY_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddeffld(gdId, GEO_FIELD2, DIMX_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "data" fields

      status = gddeffld(gdId, DATA_FIELDMV, DIMXYZ_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDMV
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gddeffld(gdId, DATA_FIELDMP, DIMXYZ_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDMP
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from and close the grid interface.  This step is necessary to store
! properly the grid information within the file and must be done before writing
! or reading data to or from the grid.

      status = gddetach(gdId)
      IF (status /= 0) THEN
          msr = GD_ERR // TRIM( l3mm%name ) // ' after definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdclose(gdfID)
      IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                ' after definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-open the file for writing

      gdfID = gdopen(physicalFilename, DFACC_RDWR)
      IF (gdfID == -1) THEN
          msr = MLSMSG_Fileopen // TRIM(physicalFilename) // ' for writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the grid for writing

      gdId = gdattach(gdfID, l3mm%name)
      IF (gdId == -1) THEN
          msr = 'Failed to attach to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write to fields

      start = 0
      stride = 1
      edge(1) = l3mm%nLevels
      edge(2) = l3mm%nLats
      edge(3) = l3mm%nLons

! Start Time

      status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
                       l3mm%startTime)
      IF (status /= 0) THEN
          msr = 'Failed to write startTime to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! End Time

      status = gdwrfld(gdId, GEO_FIELD3, stride(1), stride(1), stride(1), &
                       l3mm%endTime)
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Pressure, latitude, & longitude

      status = gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1), &
                        REAL(l3mm%pressure) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
                        REAL(l3mm%latitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
                        REAL(l3mm%longitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' // l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = gdwrfld( gdId, DATA_FIELDMV, start, stride, edge, &
                        REAL(l3mm%l3mmValue) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMV // ' to grid ' // &
                 l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, DATA_FIELDMP, start, stride, edge, &
                        REAL(l3mm%l3mmPrecision) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMP // ' to grid ' // &
                 l3mm%name
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the grid after writing

      status = gddetach(gdId)
      IF (status /= 0) THEN
          msr = GD_ERR // TRIM( l3mm%name ) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Close the file after writing

      status = gdclose(gdfID)
      IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Grid ' // TRIM(l3mm%name) // ' successfully written to file ' // &
            physicalFilename
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      creationFlag = .TRUE.

!------------------------------
   END SUBROUTINE OutputMMGrids
!------------------------------

!------------------------------------------------
   SUBROUTINE OutputMMDiags(physicalFilename, mm)
!------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic portion of l3mm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

      TYPE (L3MMData_T), INTENT(IN) :: mm

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DIMD_NAME = 'lDate'
      CHARACTER (LEN=*), PARAMETER :: DIMMIS_NAME = 'nMisDays'
      CHARACTER (LEN=*), PARAMETER :: DIMDM_NAME = 'lDate,nMisDays'

      CHARACTER (LEN=*), PARAMETER :: MIS_FIELD = 'MissingDays'

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, swdefgfld
      INTEGER, EXTERNAL :: swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: lsize, nsize, status, swfID, swId
      INTEGER :: edge(2), start(2), stride(2)

! Open the output file

      swfID = swopen(physicalFilename, DFACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // physicalFilename
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create the swath

      dgName = TRIM(mm%name) // 'Diagnostics'
      swId = swcreate(swfID, dgName)
      IF (swId == -1) THEN
          msr = 'Failed to create swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the dimensions

      status = swdefdim(swID, DIMR_NAME, MIN_MAX)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMR_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF (mm%nMisDays == 0) THEN
         nsize =1
         lsize = 4
      ELSE
         nsize = mm%nMisDays
         lsize = DATE_LEN
      ENDIF

      status = swdefdim(swId, DIMD_NAME, lsize)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMD_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdim(swId, DIMMIS_NAME, nsize)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMMIS_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdim(swId, DIM_NAME2, mm%nLevels)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIM_NAME2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "geolocation" fields

      status = swdefgfld(swID, GEO_FIELD3, DIMR_NAME, DFNT_FLOAT64, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swId, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "data" fields

      status = swdefdfld(swId, DG_FIELD2, DIM_NAME2, DFNT_INT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DG_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdfld(swId, MIS_FIELD, DIMDM_NAME, DFNT_CHAR8, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // MIS_FIELD
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface after definition.

      status = swdetach(swId)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after L3MM dg definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the swath for writing

      swId = swattach(swfID, dgName)
      IF (swId == -1) THEN
          msr = 'Failed to re-attach to swath ' // TRIM(dgName) // ' for writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write to fields

      start = 0
      stride = 1
      edge(1) = mm%nLevels

! Geolocation -- time & pressure

      status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
                       mm%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
                       mm%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1), &
                        REAL(mm%pressure) )
      IF (status /= 0) THEN
          msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = swwrfld(swId, DG_FIELD2, start(1), stride(1), edge(1), &
                       mm%perMisPoints)
      IF (status /= 0) THEN
          msr = WR_ERR //  DG_FIELD2 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = lsize
      edge(2) = nsize
      status = swwrfld(swId, MIS_FIELD, start, stride, edge, mm%misDays)
      IF (status /= 0) THEN
          msr = WR_ERR //  MIS_FIELD // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath after writing

      status = swdetach(swId)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Close the file after writing

      status = swclose(swfID)
      IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath ' // TRIM(dgName) // ' successfully written to file ' // &
            physicalFilename
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!------------------------------
   END SUBROUTINE OutputMMDiags
!------------------------------

!---------------------------------------------------------------
   SUBROUTINE WriteMetaL3MM (file, mlspcf_mcf_l3mm, pcf, anText)
!---------------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for the l3mm file, and annotates it with the
! PCF.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: file

      CHARACTER (LEN=1), POINTER :: anText(:)

      INTEGER, INTENT(IN) :: mlspcf_mcf_l3mm

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i, pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: gdinqgrid, pgs_pc_getUniversalRef

! Variables

      CHARACTER (LEN=3) :: cNum
      CHARACTER (LEN=45) :: attrName
      CHARACTER (LEN=150) :: sval
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=GridNameLen) :: gridName
      CHARACTER (LEN=98) :: inpt(31)
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: dg, hdfReturn, i, indx, j, len, numGrids, result, returnStatus
      INTEGER :: sdid, version

      REAL(r8) :: dval

! Check to see whether this is a Diagnostic product

      dg = 0

      IF ( INDEX(file,'Diagnostic') /= 0 ) dg = 1

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf_l3mm, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! Open the HDF file and initialize the SD interface

      sdid = sfstart(file, DFACC_WRITE)
      IF (sdid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                  &open the HDF file for metadata writing.')

! Set PGE values -- ECSDataGranule

      attrName = 'ReprocessingPlanned'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                               'further update anticipated using enhanced PGE')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'ReprocessingActual'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 'processed once')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'LocalGranuleID'
      indx = INDEX(file, '/', .TRUE.)
      sval = file(indx+1:)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'DayNightFlag'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, 'Both')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'LocalVersionID'
      CALL ExpandFileTemplate('$cycle', sval, cycle=pcf%cycle)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! MeasuredParameterContainer -- find the number of grids in the file

      numGrids = gdinqgrid(file, list, len)
      IF (numGrids == -1) THEN
         msr = 'No grids found in file ' // TRIM(file) // &
               ' while attempting to write its metadata.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! For each grid in the file

      DO i = 1, numGrids

! Extract its name

         indx = INDEX(list, ',')
         IF (indx /= 0) THEN
            gridName = list(:indx-1)
            list = list(indx+1:)
         ELSE
            gridName = list
         ENDIF

! Append a class suffix to ParameterName, and write the grid name as its value

         IF (i < 10) THEN
            WRITE( cNum, '(I1)' ) i
         ELSE IF ( (i >= 10) .AND. (i < 100) ) THEN
            WRITE( cNum, '(I2)' ) i
         ELSE
            WRITE( cNum, '(I3)' ) i
         ENDIF

         attrName = 'ParameterName' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                    gridName)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! QAFlags Group

         attrName = 'AutomaticQualityFlag' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                    'Passed')
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'AutomaticQualityFlagExplanation' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                    'pending algorithm update')
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'OperationalQualityFlag' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                    'Not Investigated')
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'OperationalQualityFlagExplanation' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                           'Not Investigated')
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! QAStats Group

         attrName = 'QAPercentInterpolatedData' // '.' // cNum
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'QAPercentMissingData' // '.' // cNum
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'QAPercentOutofBoundsData' // '.' // cNum
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

! OrbitCalculatedSpatialDomainContainer

      attrName = 'OrbitNumber' // '.1'
      result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 999)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'StartOrbitNumber' // '.1'
      result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'StopOrbitNumber' // '.1'
      result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EquatorCrossingLongitude' // '.1'
      dval = 0.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), attrName, dval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EquatorCrossingTime' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 '00:00:00')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EquatorCrossingDate' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 '1899-04-29')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! InputPointer

      attrName = 'InputPointer'
      inpt = ''
      IF (dg == 1) THEN
         j = mlspcf_l2dg_start
      ELSE
         j = mlspcf_l2gp_start
      ENDIF
      DO i = 1, 31
         version = 1
         returnStatus = pgs_pc_getUniversalRef(j, version, sval)
         IF (returnStatus == PGS_S_SUCCESS) THEN
            inpt(i) = sval
            j = j + 1
         ENDIF
      ENDDO

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, inpt)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Locality Value

      attrName = 'LocalityValue'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 'Limb')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! VerticalSpatialDomain Product-Specific Attribute

      attrName = 'VerticalSpatialDomainType' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 'Atmosphere Layer')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'VerticalSpatialDomainValue' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 'Atmosphere Profile')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! HorizontalSpatialDomainContainer

      attrName = 'ZoneIdentifier'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 'Other Grid System')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'WestBoundingCoordinate'
      dval = -180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), attrName, dval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'NorthBoundingCoordinate'
      dval = 90.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), attrName, dval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EastBoundingCoordinate'
      dval = 180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), attrName, dval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'SouthBoundingCoordinate'
      dval = -90.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), attrName, dval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! RangeDateTime Group

      attrName = 'RangeBeginningDate'
      result = pgs_met_setAttr_s( groups(INVENTORYMETADATA), attrName, &
                                  pcf%startDay )
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'RangeBeginningTime'
      sval = '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'RangeEndingDate'
      result = pgs_met_setAttr_s( groups(INVENTORYMETADATA), attrName, &
                                  pcf%endDay )
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'RangeEndingTime'
      sval= '23:59:59.999999'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! PGEVersion

      attrName = 'PGEVersion'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 pcf%outputVersion)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the metadata and their values to HDF attributes

      result = pgs_met_write(groups(INVENTORYMETADATA), "coremetadata", &
                             sdid)
      IF (result /= PGS_S_SUCCESS) THEN
         IF (result == PGSMET_E_MAND_NOT_SET) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, 'Some of the &
                            &mandatory metadata parameters were not set.')
         ELSE
           CALL MLSMessage(MLSMSG_Error, ModuleName, 'Metadata write &
                           &failed.')
         ENDIF
      ENDIF

! Terminate access to the SD interface and close the file

      hdfReturn = sfend(sdid)
      IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Error closing HDF file after writing metadata.')

! Annotate the file with the PCF

      CALL WritePCF2Hdr(file, anText)

      result = pgs_met_remove()

!------------------------------
   END SUBROUTINE WriteMetaL3MM
!------------------------------

!--------------------------------------------------
   SUBROUTINE AllocateL3MM (nlev, nlat, nlon, l3mm)
!--------------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3MMData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat, nlon

      TYPE( L3MMData_T ), INTENT(INOUT) :: l3mm

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3mm%nLevels = nlev
      l3mm%nLats = nlat
      l3mm%nLons = nlon

! Horizontal geolocation fields

      ALLOCATE(l3mm%latitude(l3mm%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3mm%longitude(l3mm%nLons), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' longitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Vertical geolocation field

      ALLOCATE(l3mm%pressure(l3mm%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      ALLOCATE(l3mm%l3mmValue(l3mm%nLevels,l3mm%nLats,l3mm%nLons), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' value pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3mm%l3mmPrecision(l3mm%nLevels,l3mm%nLats,l3mm%nLons),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' precision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Diagnostic fields

      ALLOCATE(l3mm%perMisPoints(l3mm%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' perMisPoints pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE AllocateL3MM
!-----------------------------

!----------------------------------
   SUBROUTINE DeallocateL3MM (l3mm)
!----------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L3MMData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3MMData_T ), INTENT(INOUT) :: l3mm

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation fields

      IF ( ASSOCIATED(l3mm%latitude) ) THEN
         DEALLOCATE (l3mm%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  latitude pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3mm%longitude) ) THEN
         DEALLOCATE (l3mm%longitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  longitude pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3mm%pressure) ) THEN
         DEALLOCATE (l3mm%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  pressure pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3mm%l3mmValue) ) THEN
         DEALLOCATE (l3mm%l3mmValue, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mmValue pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3mm%l3mmPrecision) ) THEN
         DEALLOCATE (l3mm%l3mmPrecision, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mmPrecision'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Diagnostic fields

      IF ( ASSOCIATED(l3mm%perMisPoints) ) THEN
         DEALLOCATE (l3mm%perMisPoints, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  perMisPoints pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3MM
!-------------------------------

!==================
END MODULE L3MMData
!==================

!# $Log$
!# Revision 1.2  2001/11/12 20:24:56  nakamura
!# Added L3MMDiag_T.
!#
!# Revision 1.1  2001/07/18 15:41:57  nakamura
!# Module for the L3MM data type.
!#
!#
