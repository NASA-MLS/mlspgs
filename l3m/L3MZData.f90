
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3MZData
!===============================================================================

   USE Hdf
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE mon_Open
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

! Definition -- L3MZData_T
! Subroutines -- OutputL3MZ
!                WriteMetaL3MZ
!                AllocateL3MZ
!                DeallocateL3MZ
!                DestroyL3MZDatabase

! Remarks:  This module contains the definition of the L3MZData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3mzValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3mzPrecision'

! This data type is used to store the l3 monthly zonal mean data.

  TYPE L3MZData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure   ! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude    ! dimensioned (nLats)

     REAL(r8) :: startTime	! start time of L2 data used in the analysis
     REAL(r8) :: endTime 	! end time of L2 data used in the analysis

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:), POINTER :: l3mzValue       ! Field value
     REAL(r8), DIMENSION(:,:), POINTER :: l3mzPrecision   ! Field precision
        ! dimensioned as (nLevels, nLats)

   END TYPE L3MZData_T

CONTAINS

!----------------------------------------------
   SUBROUTINE OutputL3MZ (file, mz, createFlag)
!----------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3dz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      TYPE( L3MZData_T ), INTENT(IN) :: mz

      LOGICAL, INTENT(INOUT) :: createFlag

! Parameters

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim
      INTEGER, EXTERNAL :: swdefgfld, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: status, swfID, swID
      INTEGER :: start(2), stride(2), edge(2)

! Open the file for appending a swath

      swfID = swopen(file, DFACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // file
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create a swath of the appropriate name

      swID = swcreate(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to create swath ' // mz%name
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath dimensions

      status = swdefdim(swID, DIMT_NAME, 2)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdim(swID, DIM_NAME2, mz%nLevels)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdim(swID, DIML_NAME, mz%nLats)
      IF (status == -1) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath geolocation fields using the above dimensions

      status = swdefgfld(swID, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath data fields using the above dimensions

      status = swdefdfld(swID, DATA_FIELDV, DIMLL_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdfld(swID, DATA_FIELDP, DIMLL_NAME, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface after definition

      status = swdetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                              &detach from swath interface after L3DZ definition.')

! Re-attach to the swath for writing

      swID = swattach(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to swath ' // TRIM(mz%name) // ' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = mz%nLevels
      edge(2) = mz%nLats

! Geolocation fields

      status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
                       mz%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
                       mz%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
                        REAL(mz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
                        REAL(mz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = swwrfld( swID, DATA_FIELDV, start, stride, edge, REAL(mz%l3mzValue) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, DATA_FIELDP, start, stride, edge, &
                        REAL(mz%l3mzPrecision) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from swath interface

      status = swdetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                 &detach from swath interface after writing.')

      msr = 'Swath ' // TRIM(mz%name) // ' successfully written to file ' // file
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

! Close the file

      status = swclose(swfID)
      IF (status == -1) THEN
         msr = 'Failed to close file ' // file // ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      createFlag = .TRUE.

!---------------------------
   END SUBROUTINE OutputL3MZ
!---------------------------

!----------------------------------------------------------
   SUBROUTINE WriteMetaL3MZ (fileName, mcfNum, pcf, anText)
!----------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3mz file, and annotates it with the
! PCF.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER(LEN=*), INTENT(IN) :: fileName

      INTEGER, INTENT(IN) :: mcfNum

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: swinqswath, pgs_pc_getUniversalRef

! Variables

      CHARACTER (LEN=3) :: cNum
      CHARACTER (LEN=45) :: attrName, lvid
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=FileNameLen) :: sval
      CHARACTER (LEN=GridNameLen) :: swathName
      CHARACTER (LEN=98) :: inpt(31)
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: dg, hdfReturn, i, indx, len, numSwaths, result, returnStatus, sdid
      INTEGER :: version

      REAL(r8) :: dval

! Check to see whether this is the Diagnostic files

      dg = 0

      IF ( INDEX(fileName,'Diagnostic') /= 0 ) dg = 1

! Initialize the MCF

      result = pgs_met_init(mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! Open the HDF file and initialize the SD interface

      sdid = sfstart(fileName, DFACC_WRITE)
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
      indx = INDEX(fileName, '/', .TRUE.)
      sval = fileName(indx+1:)
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
      CALL ExpandFileTemplate('$cycle', lvid, cycle=pcf%cycle)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, lvid)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! MeasuredParameterContainer -- find the number of swaths in the file

      numSwaths = swinqswath(fileName, list, len)
      IF (numSwaths == -1) THEN
         msr = 'No swaths found in file ' // TRIM(fileName) // &
               ' while attempting to write its metadata.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! For each swath in the file

      DO i = 1, numSwaths

! Extract its name

         indx = INDEX(list, ',')
         IF (indx /= 0) THEN
            swathName = list(:indx-1)
            list = list(indx+1:)
         ELSE
            swathName = list
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
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, swathName)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! QAFlags Group

         attrName = 'AutomaticQualityFlag' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, 'Passed')
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
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, '00:00:00')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EquatorCrossingDate' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, '1899-04-29')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! InputPointer

      attrName = 'InputPointer'
      inpt = ''
      IF (dg == 1) THEN
         indx = mlspcf_l2dg_start
      ELSE
         indx = mlspcf_l2gp_start
      ENDIF
      DO i = 1, 31
         version = 1
         returnStatus = pgs_pc_getUniversalRef(indx, version, sval)
         IF (returnStatus == PGS_S_SUCCESS) THEN
            inpt(i) = sval
            indx = indx + 1
         ENDIF
      ENDDO

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, inpt)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Locality Value

      attrName = 'LocalityValue'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, 'Limb')
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
                                  pcf%startDay)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'RangeBeginningTime'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 '00:00:00.000000')
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
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                 '23:59:59.999999')
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

      result = pgs_met_write(groups(INVENTORYMETADATA), "coremetadata", sdid)
      IF (result /= PGS_S_SUCCESS) THEN
         IF (result == PGSMET_E_MAND_NOT_SET) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, 'Some of the mandatory &
                                                &metadata parameters were not set.')
         ELSE
              CALL MLSMessage(MLSMSG_Error, ModuleName, 'Metadata write failed.')
         ENDIF
      ENDIF

! Terminate access to the SD interface and close the file

      hdfReturn = sfend(sdid)
      IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                   'Error closing HDF file after writing metadata.')

! Annotate the file with the PCF

      CALL WritePCF2Hdr(fileName, anText)

      result = pgs_met_remove()

!------------------------------
   END SUBROUTINE WriteMetaL3MZ
!------------------------------

!--------------------------------------------
   SUBROUTINE AllocateL3MZ (nlev, nlat, l3mz)
!--------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3MZData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat

      TYPE( L3MZData_T ), INTENT(INOUT) :: l3mz

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3mz%nLevels = nlev
      l3mz%nLats = nlat

! Allocate the vertical geolocation field

      ALLOCATE(l3mz%pressure(l3mz%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3MZData_T pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      ALLOCATE(l3mz%latitude(l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3MZData_T latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      ALLOCATE(l3mz%l3mzValue(l3mz%nLevels,l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3mzValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3mz%l3mzPrecision(l3mz%nLevels,l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3mzPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE AllocateL3MZ
!-----------------------------

!----------------------------------
   SUBROUTINE DeallocateL3MZ (l3mz)
!----------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L3MZData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3MZData_T ), INTENT(INOUT) :: l3mz

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation field

      IF ( ASSOCIATED(l3mz%latitude) ) THEN
         DEALLOCATE (l3mz%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mz latitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3mz%pressure) ) THEN
         DEALLOCATE (l3mz%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mz pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3mz%l3mzValue) ) THEN
         DEALLOCATE (l3mz%l3mzValue, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mzValue pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3mz%l3mzPrecision) ) THEN
         DEALLOCATE (l3mz%l3mzPrecision, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mzPrecision pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3MZ
!-------------------------------

!---------------------------------------
   SUBROUTINE DestroyL3MZDatabase (mzdb)
!---------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal structures of an l3mz database, and
! then database itself


! Arguments

      TYPE (L3MZData_T), DIMENSION(:), POINTER :: mzdb

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i

! Check the status of the input pointer

      IF ( ASSOCIATED(mzdb) ) THEN

! If it's associated, then deallocate the internal structures

         DO i = 1, SIZE(mzdb)
            CALL DeallocateL3MZ( mzdb(i) )
         ENDDO

! Deallocate the database itself

         DEALLOCATE (mzdb, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3mz database'
            CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

!------------------------------------
   END SUBROUTINE DestroyL3MZDatabase
!------------------------------------

!==================
END MODULE L3MZData
!==================

!# $Log$
!#
