
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DZData
!===============================================================================

   USE Hdf
   USE L2GPData, ONLY: DIM_NAME2, GEO_FIELD1, GEO_FIELD3, GEO_FIELD9, &
                       HDFE_NOMERGE, L2GPNameLen
   USE L3CF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSStrings
   USE MLSPCF
   USE OpenInit
   USE PCFModule
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L3DZData_T
!               L3DZFiles_T 
! Subroutines -- OutputL3DZ
!                WriteMetaL3DZS
!                WriteMetaL3DZD
!                AllocateL3DZ
!                DeallocateL3DZ
!                DestroyL3DZDatabase

! Remarks:  This module contains the definition of the L3DZData type, as well
!           as any routines pertaining to it.

! Parameters

! This data type is used to store the l3 daily zonal mean data.

  TYPE L3DZData_T

     CHARACTER (LEN=L2GPNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure   ! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude    ! dimensioned (nLats)

     REAL(r8) :: time                           ! Synoptic time

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:), POINTER :: l3dzValue       ! Field value
     REAL(r8), DIMENSION(:,:), POINTER :: l3dzPrecision   ! Field precision
        ! dimensioned as (nLevels, nLats)

   END TYPE L3DZData_T

! This data type is used to store information about the number of distinct l3dz
! files of each type actually created.

   TYPE L3DZFiles_T

     INTEGER :: nStd	! number of distinct Standard files created
     INTEGER :: nDg	! number of distinct Diagnostic files created

     ! names of the files created

     CHARACTER (LEN=FileNameLen) :: stdNames(maxWindow)
     CHARACTER (LEN=FileNameLen) :: DgNames(maxWindow)

     ! CCSDS B format dates of the created files

     CHARACTER (LEN=8) :: stdDates(maxWindow)
     CHARACTER (LEN=8) :: dgDates(maxWindow)

   END TYPE L3DZFiles_T

CONTAINS

!------------------------------------------
   SUBROUTINE OutputL3DZ(type, dzm, zFiles)
!------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3dz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE( L3DZData_T ), INTENT(IN) :: dzm(:)

      TYPE( L3DZFiles_T ), INTENT(OUT) :: zFiles

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3dzValue'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3dzPrecision'

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim
      INTEGER, EXTERNAL :: swdefgfld, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: file

      INTEGER :: i, match, numDays, status, swfID, swID
      INTEGER :: start(2), stride(2), edge(2)

! For each day in the output window,

      numDays = SIZE(dzm)

      DO i = 1, numDays

! Get the name of the L3DZ file for the proper type & day from the PCF

         CALL FindFileDay(type, dzm(i)%time, mlspcf_l3dz_start, &
                          mlspcf_l3dz_end, match, file, date)
         IF (match == -1) THEN
            msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
                  date
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Check whether the name is distinct; if so, save it in L3DZFiles_T under the
! correct type

         IF ( INDEX(type,'Diagnostic') /= 0 ) THEN

            IF (LinearSearchStringArray(zFiles%dgNames,file) == 0) THEN
               zFiles%nDg= zFiles%nDg+1
               zFiles%dgNames(zFiles%nDg) = file
               zFiles%dgDates(zFiles%nDg) = date
            ENDIF

         ELSE

            IF (LinearSearchStringArray(zFiles%StdNames,file) == 0) THEN
               zFiles%nStd = zFiles%nStd+1
               zFiles%stdNames(zFiles%nStd) = file
               zFiles%stdDates(zFiles%nStd) = date
            ENDIF

         ENDIF

! Open the file for appending a swath

         swfID = swopen(file, DFACC_RDWR)
         IF (swfID == -1) THEN
            msr = MLSMSG_Fileopen // file
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Create a swath of the appropriate name

         swID = swcreate(swfID, dzm(i)%name)
         IF (swID == -1) THEN
            msr = 'Failed to create swath ' // dzm(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the swath dimensions

         status = swdefdim(swID, DIMT_NAME, 1)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMT_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swID, DIM_NAME2, dzm(i)%nLevels)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swID, DIML_NAME, dzm(i)%nLats)
         IF (status == -1) THEN
            msr = DIM_ERR // DIML_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the swath geolocation fields using the above dimensions

         status = swdefgfld(swID, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the swath data fields using the above dimensions

         status = swdefdfld(swID, DATA_FIELDV, DIMDZ_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
             msr = DAT_ERR // DATA_FIELDV
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDP, DIMDZ_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
             msr = DAT_ERR // DATA_FIELDP
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from the swath interface after definition

         status = swdetach(swID)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                       &to detach from swath interface after L3DZ definition.')

! Re-attach to the swath for writing

         swID = swattach(swfID, dzm(i)%name)
         IF (swID == -1) THEN
            msr = 'Failed to re-attach to swath ' // TRIM(dzm(i)%name) &
                   // ' for writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the data

         start = 0
         stride = 1
         edge(1) = dzm(i)%nLevels
         edge(2) = dzm(i)%nLats

! Geolocation fields

         status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
                          dzm(i)%time)
         IF (status == -1) THEN
             msr = WR_ERR // GEO_FIELD3
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
                           REAL(dzm(i)%pressure) )
         IF (status == -1) THEN
             msr = WR_ERR // GEO_FIELD9
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
                           REAL(dzm(i)%latitude) )
         IF (status == -1) THEN
             msr = WR_ERR // GEO_FIELD1
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Data fields

         status = swwrfld( swID, DATA_FIELDV, start, stride, edge, &
                           REAL(dzm(i)%l3dzValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, DATA_FIELDP, start, stride, edge, &
                           REAL(dzm(i)%l3dzPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! After writing, detach from swath interface

         status = swdetach(swID)
         IF (status == -1) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to detach from &
                                           &swath interface after writing.')
         ENDIF

! Close the file

         status = swclose(swfID)
         IF (status == -1) THEN
            msr = 'Failed to close file ' // file // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         msr = 'Swath ' // TRIM(dzm(i)%name) // &
               ' successfully written to file ' // file
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      ENDDO

!---------------------------
   END SUBROUTINE OutputL3DZ
!---------------------------

!-------------------------------------------------------
   SUBROUTINE WriteMetaL3DZS (pcf, cfDef, files, anText)
!-------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dm file, and annotates it with the
! PCF.

! Arguments

      TYPE( L3CFDef_T ), INTENT(IN) :: cfDef

      TYPE (L3DZFiles_T), INTENT(IN) :: files

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: swinqswath

! Variables

      CHARACTER (LEN=1) :: cNum
      CHARACTER (LEN=12) :: shortName
      CHARACTER (LEN=45) :: attrName, lvid, sval
      CHARACTER (LEN=80) :: identifier
      CHARACTER (LEN=132) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: swathName
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: hdfReturn, i, j, indx, len, numSwaths, result, sdid

      REAL(r8) :: dval

! L3DZ Standard products file -- initialize the MCF

      result = pgs_met_init(cfDef%stdMCFnum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! For each file successfully created,

      DO i = 1, files%nStd

! Open the HDF file and initialize the SD interface

         sdid = sfstart(files%stdNames(i), DFACC_WRITE)
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
         indx = INDEX(files%stdNames(i), '.', .TRUE.)
         sval = files%stdNames(i)(:indx-1)
         indx = INDEX(sval, '/', .TRUE.)
         sval = sval(indx+1:)
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
         CALL ExpandFileTemplate('$version-$cycle', lvid, &
                                 version=pcf%outputVersion, cycle=pcf%cycle)
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, lvid)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! MeasuredParameterContainer -- find the number of swaths in the file

         numSwaths = swinqswath(files%stdNames(i), list, len)
         IF (numSwaths == -1) THEN
            msr = 'No swaths found in file ' // TRIM( files%stdNames(i) ) // &
                  ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! For each swath in the file

         DO j = 1, numSwaths

! Extract its name

            indx = INDEX(list, ',')
            IF (indx /= 0) THEN
               swathName = list(:indx-1)
               list = list(indx+1:)
            ELSE
               swathName = list
            ENDIF

! Append a class suffix to ParameterName, and write the grid name as its value

            WRITE( cNum, '(I1)' ) j

            attrName = 'ParameterName' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       swathName)
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
                                       'TBD')
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
                                              'TBD')
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

         indx = INDEX(cfDef%stdMCFname, '.', .TRUE.)
         shortName = cfDef%stdMCFname(:indx-1)
         indx = INDEX(shortName, '.')
         shortName = shortName(:indx-1) // ':' // shortName(indx+1:)
         identifier = 'MLS-Aura_L2GP_' // files%stdDates(i)
         sval = 'LGID:' // TRIM(shortName) // ':' // TRIM(identifier)
         attrName = 'InputPointer'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
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
                                     files%stdDates(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeBeginningTime'
         sval = '12:00:00.000000'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeEndingDate'
         result = pgs_met_setAttr_s( groups(INVENTORYMETADATA), attrName, &
                                     files%stdDates(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeEndingTime'
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

         CALL WritePCF2Hdr(files%stdNames(i), anText)

      ENDDO

      result = pgs_met_remove()

!-------------------------------
   END SUBROUTINE WriteMetaL3DZS
!-------------------------------

!-------------------------------------------------------
   SUBROUTINE WriteMetaL3DZD (pcf, cfDef, files, anText)
!-------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dm file, and annotates it with the
! PCF.

! Arguments

      TYPE( L3CFDef_T ), INTENT(IN) :: cfDef

      TYPE (L3DZFiles_T), INTENT(IN) :: files

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: swinqswath

! Variables

      CHARACTER (LEN=1) :: cNum
      CHARACTER (LEN=12) :: shortName
      CHARACTER (LEN=45) :: attrName, lvid, sval
      CHARACTER (LEN=80) :: identifier
      CHARACTER (LEN=132) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: swathName
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: hdfReturn, i, j, indx, len, numSwaths, result, sdid

      REAL(r8) :: dval


! L3DZ Diagnostic products file -- initialize the MCF

      IF (files%nDg == 0) RETURN

      result = pgs_met_init(cfDef%dgMCFnum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! For each file successfully created,

      DO i = 1, files%nDg

! Open the HDF file and initialize the SD interface

         sdid = sfstart(files%dgNames(i), DFACC_WRITE)
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
         indx = INDEX(files%dgNames(i), '.', .TRUE.)
         sval = files%dgNames(i)(:indx-1)
         indx = INDEX(sval, '/', .TRUE.)
         sval = sval(indx+1:)
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
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, lvid)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! MeasuredParameterContainer -- find the number of swaths in the file

         numSwaths = swinqswath(files%dgNames(i), list, len)
         IF (numSwaths == -1) THEN
            msr = 'No swaths found in file ' // TRIM( files%dgNames(i) ) // &
                  ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! For each swath in the file

         DO j = 1, numSwaths

! Extract its name

            indx = INDEX(list, ',')
            IF (indx /= 0) THEN
               swathName = list(:indx-1)
               list = list(indx+1:)
            ELSE
               swathName = list
            ENDIF

! Append a class suffix to ParameterName, and write the grid name as its value

            WRITE( cNum, '(I1)' ) j

            attrName = 'ParameterName' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       swathName)
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
                                       'TBD')
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
                                              'TBD')
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

         indx = INDEX(cfDef%dgMCFname, '.', .TRUE.)
         shortName = cfDef%dgMCFname(:indx-1)
         indx = INDEX(shortName, '.')
         shortName = shortName(:indx-1) // ':' // shortName(indx+1:)
         identifier = 'MLS-Aura_L2GP_' // files%dgDates(i)
         sval = 'LGID:' // TRIM(shortName) // ':' // TRIM(identifier)
         attrName = 'InputPointer'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
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
                                     files%dgDates(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeBeginningTime'
         sval = '12:00:00.000000'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeEndingDate'
         result = pgs_met_setAttr_s( groups(INVENTORYMETADATA), attrName, &
                                     files%dgDates(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeEndingTime'
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

         CALL WritePCF2Hdr(files%dgNames(i), anText)

      ENDDO

      result = pgs_met_remove()

!-------------------------------
   END SUBROUTINE WriteMetaL3DZD
!-------------------------------

!--------------------------------------------
   SUBROUTINE AllocateL3DZ (nlev, nlat, l3dz)
!--------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3DZData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat

      TYPE( L3DZData_T ), INTENT(INOUT) :: l3dz

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3dz%nLevels = nlev
      l3dz%nLats = nlat

! Allocate the vertical geolocation field

      ALLOCATE(l3dz%pressure(l3dz%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DZData_T pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      ALLOCATE(l3dz%latitude(l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DZData_T latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      ALLOCATE(l3dz%l3dzValue(l3dz%nLevels,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3dzValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dz%l3dzPrecision(l3dz%nLevels,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3dzPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE AllocateL3DZ
!-----------------------------


!----------------------------------
   SUBROUTINE DeallocateL3DZ (l3dz)
!----------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L2GP_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3DZData_T ), INTENT(INOUT) :: l3dz

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation field

      IF ( ASSOCIATED(l3dz%latitude) ) DEALLOCATE (l3dz%latitude, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dz latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3dz%pressure) ) DEALLOCATE (l3dz%pressure, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dz pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3dz%l3dzValue) ) DEALLOCATE (l3dz%l3dzValue, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dzValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3dz%l3dzPrecision) ) DEALLOCATE (l3dz%l3dzPrecision, &
                                                        STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dzPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3DZ
!-------------------------------

!-----------------------------------------
   SUBROUTINE DestroyL3DZDatabase (l3dzdb)
!-----------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal structures of an l3dz database, and
! then database itself


! Arguments

      TYPE (L3DZData_T), DIMENSION(:), POINTER :: l3dzdb

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i

! Check the status of the input pointer

      IF ( ASSOCIATED(l3dzdb) ) THEN

! If it's associated, then deallocate the internal structures

         DO i = 1, SIZE(l3dzdb)
            CALL DeallocateL3DZ( l3dzdb(i) )
         ENDDO

! Deallocate the database itself

         DEALLOCATE (l3dzdb, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dz database'
            CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

!------------------------------------
   END SUBROUTINE DestroyL3DZDatabase
!------------------------------------

!==================
END MODULE L3DZData
!==================

!# $Log$
!# Revision 1.2  2001/01/26 17:32:54  nakamura
!# Added version of OutputL3DZ without files.
!#
!#
