
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DZData
!===============================================================================

   USE Hdf
   USE L3CF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSStrings
   USE MLSPCF3
   USE OpenInit
   USE PCFHdr
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
! Subroutines -- OutputL3DZ
!                ReadL3DZData
!                WriteMetaL3DZ
!                AllocateL3DZ
!                DeallocateL3DZ
!                DestroyL3DZDatabase

! Remarks:  This module contains the definition of the L3DZData type, as well
!           as any routines pertaining to it.

! Parameters

  CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3dzValue'
  CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3dzPrecision'

! This data type is used to store the l3 daily zonal mean data.

  TYPE L3DZData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)

     REAL(r8) :: time						! Synoptic time

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:), POINTER :: l3dzValue		! Field value
     REAL(r8), DIMENSION(:,:), POINTER :: l3dzPrecision		! Field precision
        ! dimensioned as (nLevels, nLats)

   END TYPE L3DZData_T

CONTAINS

!------------------------------------------
   SUBROUTINE OutputL3DZ(type, dzm, zFiles)
!------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3dz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE( L3DZData_T ), INTENT(IN) :: dzm(:)

      TYPE( OutputFiles_T ), INTENT(OUT) :: zFiles

! Parameters

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

! Check whether the name is distinct; if so, save it in zFiles

         IF (LinearSearchStringArray(zFiles%name,file) == 0) THEN
            zFiles%nFiles = zFiles%nFiles+1
            zFiles%name(zFiles%nFiles) = file
            zFiles%date(zFiles%nFiles) = date
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

         status = swdefdim(swID, DIMT_NAME, DATE_LEN)
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

         status = swdefgfld(swID, GEO_FIELD11, DIMT_NAME, DFNT_CHAR8, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD11
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

         status = swdefdfld(swID, DATA_FIELDV, DIMLL_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
             msr = DAT_ERR // DATA_FIELDV
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDP, DIMLL_NAME, DFNT_FLOAT32, &
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
         edge(1) = DATE_LEN
         edge(2) = dzm(i)%nLats

! Geolocation fields

         status = swwrfld(swID, GEO_FIELD11, start(1), stride(1), edge(1), date)
         IF (status == -1) THEN
             msr = WR_ERR // GEO_FIELD11
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = dzm(i)%nLevels

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

!------------------------------------------------
   SUBROUTINE ReadL3DZData (swfid, swathName, dz)
!------------------------------------------------

! Brief description of subroutine
! This subroutine a data swath from a file into the L3DZData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: swathName

      INTEGER, INTENT(IN) :: swfid

      TYPE( L3DZData_T ), INTENT(OUT) :: dz

! Parameters

      CHARACTER (LEN=*), PARAMETER :: RL3DZ_ERR = 'Failed to read L3DZ field:  '

! Functions

      INTEGER, EXTERNAL :: pgs_td_utcToTAI, swattach, swdetach, swdiminfo, swrdfld

! Variables

      CHARACTER (LEN=DATE_LEN) :: date
      CHARACTER (LEN=CCSDSB_LEN) :: dTime
      CHARACTER (LEN=480) :: msr

      INTEGER :: err, nlat, nlev, returnStatus, size, swid, status
      INTEGER :: start(2), stride(2), edge(2)

      REAL, ALLOCATABLE :: rl(:), rp(:)
      REAL, ALLOCATABLE :: r2(:,:)

! Attach to the swath for reading

      swid = swattach(swfid, swathName)
      IF (status == -1) THEN
         msr = 'Failed to attach to swath interface for reading ' // &
               swathName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get dimension information

      size = swdiminfo(swid, DIM_NAME2)
      IF (size == -1) THEN
         msr = SZ_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlev = size

      size = swdiminfo(swid, DIML_NAME)
      IF (size == -1) THEN
         msr = SZ_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlat = size

! Allocate the output structure and temporary variables

      CALL AllocateL3DZ(nlev, nlat, dz)

      ALLOCATE (rp(dz%nLevels), rl(dz%nLats), r2(dz%nLevels,dz%nLats), &
                STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local REAL variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Read the geolocation fields

      dz%name = swathName

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats

      status = swrdfld(swid, GEO_FIELD11, start(1), stride(1), edge(1), date)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dTime = date // 'T12:00:00.000000Z'
      returnStatus = pgs_td_utcToTAI(dTime, dz%time)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                 ModuleName, 'Error converting from CCSDSB to TAI.')

      edge(1) = dz%nLevels

      status = swrdfld(swid, GEO_FIELD9, start(1), stride(1), edge(1), rp)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%pressure = DBLE(rp)

      status = swrdfld(swid, GEO_FIELD1, start(2), stride(2), edge(2), rl)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latitude = DBLE(rl)

! Read the data fields

      status = swrdfld(swid, DATA_FIELDV, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%l3dzValue = DBLE(r2)

      status = swrdfld(swid, DATA_FIELDP, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%l3dzPrecision = DBLE(r2)

!  After reading, detach from swath interface

      status = swdetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from swath interface after reading ' // &
               swathName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Deallocate local variables

      DEALLOCATE(rp, rl, r2, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE ReadL3DZData
!-----------------------------

!-------------------------------------------------------
   SUBROUTINE WriteMetaL3DZ (pcf, mcfNum, files, anText)
!-------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dz file, and annotates it with the
! PCF.

! Arguments

      TYPE (OutputFiles_T), INTENT(IN) :: files

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      INTEGER, INTENT(IN) :: mcfNum

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: swinqswath

! Variables

      CHARACTER (LEN=3) :: cNum
      CHARACTER (LEN=45) :: attrName, lvid
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=FileNameLen) :: sval
      CHARACTER (LEN=GridNameLen) :: swathName
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: hdfReturn, i, j, indx, len, numSwaths, result, sdid

      REAL(r8) :: dval

! If no L3DZ files were created, exit

      IF (files%nFiles == 0) RETURN

! Initialize the MCF

      result = pgs_met_init(mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! For each file successfully created,

      DO i = 1, files%nFiles

! Open the HDF file and initialize the SD interface

         sdid = sfstart(files%name(i), DFACC_WRITE)
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
         indx = INDEX(files%name(i), '/', .TRUE.)
         sval = files%name(i)(indx+1:)
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

         numSwaths = swinqswath(files%name(i), list, len)
         IF (numSwaths == -1) THEN
            msr = 'No swaths found in file ' // TRIM( files%name(i) ) // &
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

            IF (j < 10) THEN
               WRITE( cNum, '(I1)' ) j
            ELSE IF ( (j >= 10) .AND. (j < 100) ) THEN
               WRITE( cNum, '(I2)' ) j
            ELSE
               WRITE( cNum, '(I3)' ) j
            ENDIF

            attrName = 'ParameterName' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       swathName)
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

! QAFlags Group

            attrName = 'AutomaticQualityFlag' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       'Passed')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'AutomaticQualityFlagExplanation' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       'pending algorithm update')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'OperationalQualityFlag' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                       'Not Investigated')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'OperationalQualityFlagExplanation' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                              'Not Investigated')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

! QAStats Group

            attrName = 'QAPercentInterpolatedData' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'QAPercentMissingData' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'QAPercentOutofBoundsData' // '.' // TRIM(cNum)
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
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                                   'See the PCF annotation to this file.')
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
                                     files%date(i) )
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
                                     files%date(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'RangeEndingTime'
         sval = '23:59:59.999999'
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

         CALL WritePCF2Hdr(files%name(i), anText)

      ENDDO

      result = pgs_met_remove()

!------------------------------
   END SUBROUTINE WriteMetaL3DZ
!------------------------------

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
! This subroutine deallocates the internal field pointers of the L3DZData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3DZData_T ), INTENT(INOUT) :: l3dz

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation field

      IF ( ASSOCIATED(l3dz%latitude) ) THEN
         DEALLOCATE (l3dz%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dz latitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3dz%pressure) ) THEN
         DEALLOCATE (l3dz%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dz pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3dz%l3dzValue) ) THEN
         DEALLOCATE (l3dz%l3dzValue, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dzValue pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dz%l3dzPrecision) ) THEN
         DEALLOCATE (l3dz%l3dzPrecision, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dzPrecision pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
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
!# Revision 1.4  2001/04/24 19:38:50  nakamura
!# Removed references to private L2 parameters.
!#
!# Revision 1.3  2001/03/27 19:29:13  nakamura
!# Updated the metadata; fixed err checks on deallocate.
!#
!# Revision 1.2  2001/02/21 21:02:50  nakamura
!# Changed MLSPCF to MLSPCF3; shifted some parameters; added ReadL3DZData; allowed for expanded swath list; changed InputPointer.
!#
!# Revision 1.1  2001/02/12 19:23:42  nakamura
!# Module for the L3DZ data type.
!#
