
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DMData
!===============================================================================

   USE Hdf
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE MLSStrings
   USE OpenInit
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

! Definition -- L3DMData_T
! Subroutines -- ConvertDeg2DMS
!                OutputGrids
!                ReadL3DMData
!                WriteMetaL3DM
!                AllocateL3DM
!                DeallocateL3DM
!                DestroyL3DMDatabase

! Remarks:  This module contains the definition of the L3DMData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3dmValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3dmPrecision'

! This data type is used to store the l3 daily map data.

   TYPE L3DMData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes
     INTEGER :: nLons				! Total number of longitudes

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)
     REAL(r8), DIMENSION(:), POINTER :: longitude	! dimensioned (nLons)

     REAL(r8) :: time	! Synoptic time

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:,:), POINTER :: l3dmValue	  ! Field value
     REAL(r8), DIMENSION(:,:,:), POINTER :: l3dmPrecision ! Field precision
	! dimensioned as (nLevels, nLats, nLons)

   END TYPE L3DMData_T

CONTAINS

!--------------------------------------
   SUBROUTINE ConvertDeg2DMS (deg, dms)
!--------------------------------------

! Brief description of subroutine
! This subroutine converts degrees in "decimal format" (DDD.ddd) to "packed
! degree format" (DDDMMMSSS.ss)

! Arguments

      REAL, INTENT(IN) :: deg
 
      REAL(r8), INTENT(OUT) :: dms

! Parameters

! Functions

! Variables

      REAL(r8) :: decDeg, decMin, d, dm, intDeg, intMin, min, sec

! Convert DDD.ddd to DDD000000.00

      intDeg = AINT(deg)

      d = intDeg * 1000000.00

! Convert .ddd to MM.mmm

      decDeg = deg - intDeg

      min = decDeg * 60.0

! Truncate MM.mmm

      intMin = AINT(min)

! Add 0MM000.00 to DDD000000.00

      dm = intMin * 1000.0 + d

! Convert .mmm to SS.sss

      decMin = min - intMin

      sec = decMin * 60.0

! Add to get DDD0MM0SS.ss

      dms = dm + sec

!-------------------------------
   END SUBROUTINE ConvertDeg2DMS
!-------------------------------

!---------------------------------------------------
   SUBROUTINE OutputGrids(type, l3dmData, l3dmFiles)
!---------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the grid portion of the l3dm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE (L3DMData_T), INTENT(IN) :: l3dmData(:)

      TYPE (OutputFiles_T), INTENT(INOUT) :: l3dmFiles
   
! Parameters

! Functions

      INTEGER, EXTERNAL :: gdattach, gdclose, gdcreate, gddefdim, gddeffld
      INTEGER, EXTERNAL :: gddefproj, gddetach, gdopen, gdwrfld

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=FileNameLen) :: physicalFilename
      CHARACTER (LEN=480) :: msr

      INTEGER :: gdfID, gdId, i, match, status
      INTEGER :: start(3), stride(3), edge(3)

      REAL :: maxLat, maxLon, minLat, minLon

      REAL(r8) :: uplft(2), lowrgt(2)
      REAL(r8) :: projparm(13)

! For each day in the l3dm database,

      DO i = 1, SIZE(l3dmData)

! Find the output file for the level/species/day in the PCF

         CALL FindFileDay(type, l3dmData(i)%time, mlspcf_l3dm_start, &
                          mlspcf_l3dm_end, match, physicalFilename, date)
         IF (match == -1) THEN
            msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
                  date
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Check whether the name is distinct; if so, save it in l3dmFiles

         IF (LinearSearchStringArray(l3dmFiles%name,physicalFilename) == 0) THEN
            l3dmFiles%nFiles = l3dmFiles%nFiles+1
            l3dmFiles%name(l3dmFiles%nFiles) = physicalFilename
            l3dmFiles%date(l3dmFiles%nFiles) = date
         ENDIF

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

         maxLat = MAXVAL( l3dmData(i)%latitude )
         minLat = MINVAL( l3dmData(i)%latitude )
         maxLon = MAXVAL( l3dmData(i)%longitude )
         minLon = MINVAL( l3dmData(i)%longitude )

! Convert to "packed degree format"

         CALL ConvertDeg2DMS(maxLat, uplft(2))
         CALL ConvertDeg2DMS(minLat, lowrgt(2))
         CALL ConvertDeg2DMS(maxLon, lowrgt(1))
         CALL ConvertDeg2DMS(minLon, uplft(1))

! Create the grid

         gdId = gdcreate(gdfID, l3dmData(i)%name, l3dmData(i)%nLons, &
                         l3dmData(i)%nLats, uplft, lowrgt)
         IF (gdId == -1) THEN
            msr = 'Failed to create grid ' // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the dimensions

         status = gddefdim(gdId, DIMX_NAME, l3dmData(i)%nLons)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMX_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddefdim(gdId, DIMY_NAME, l3dmData(i)%nLats)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMY_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddefdim(gdId, DIMZ_NAME, l3dmData(i)%nLevels)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMZ_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddefdim(gdId, DIMT_NAME, 1)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMT_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define a "Geographic projection," using defaults in all unneeded fields

         status = gddefproj(gdId, GCTP_GEO, 0, 0, projparm)
         IF (status /= 0) THEN
            msr = 'Failed to define projection for grid ' // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the "geolocation" fields

         status = gddeffld(gdId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddeffld(gdId, GEO_FIELD9, DIMZ_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddeffld(gdId, GEO_FIELD1, DIMY_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddeffld(gdId, GEO_FIELD2, DIMX_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the "data" fields

         status = gddeffld(gdId, DATA_FIELDV, DIMXYZ_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddeffld(gdId, DATA_FIELDP, DIMXYZ_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from and close the grid interface.  This step is necessary to store
! properly the grid information within the file and must be done before writing
! or reading data to or from the grid.

         status = gddetach(gdId)
         IF (status /= 0) THEN
            msr = GD_ERR // TRIM( l3dmData(i)%name ) // ' after definition.'
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

         gdId = gdattach(gdfID, l3dmData(i)%name)
         IF (gdId == -1) THEN
            msr = 'Failed to attach to grid ' // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write to fields

         start = 0
         stride = 1
         edge(1) = l3dmData(i)%nLevels
         edge(2) = l3dmData(i)%nLats
         edge(3) = l3dmData(i)%nLons

         status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
                          l3dmData(i)%time)
         IF (status /=0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' &
                   // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1), &
                           REAL(l3dmData(i)%pressure) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' &
                  // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
                           REAL(l3dmData(i)%latitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' &
                  // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
                           REAL(l3dmData(i)%longitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' &
                   // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, DATA_FIELDV, start, stride, edge, &
                           REAL(l3dmData(i)%l3dmValue) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDV // ' to grid ' &
                  // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, DATA_FIELDP, start, stride, edge, &
                           REAL(l3dmData(i)%l3dmPrecision) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDP // ' to grid ' &
                  // l3dmData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from the grid after writing

         status = gddetach(gdId)
         IF (status /= 0) THEN
            msr = GD_ERR // TRIM( l3dmData(i)%name ) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Close the file after writing

         status = gdclose(gdfID)
         IF (status /= 0) THEN
            msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                  ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         msr = 'Grid ' // TRIM(l3dmData(i)%name) // ' successfully written to &
               &file ' // physicalFilename
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      ENDDO

!----------------------------
   END SUBROUTINE OutputGrids
!----------------------------

!-------------------------------------------------
   SUBROUTINE ReadL3DMData (gdfid, gridName, l3dm)
!-------------------------------------------------

! Brief description of subroutine
! This  subroutine reads a grid from a file into the L3DMData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: gridName

      INTEGER, INTENT(IN) :: gdfid  

      TYPE( L3DMData_T ), INTENT(OUT) :: l3dm

! Parameters

      CHARACTER (LEN=*), PARAMETER :: RDGD_ERR = 'Failed to read grid field '

! Functions

      INTEGER :: gdattach, gddetach, gddiminfo, gdrdfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, gdid, nlev, nlat, nlon, size, status
      INTEGER :: start(3), stride(3), edge(3)

      REAL, ALLOCATABLE :: rp(:), rlat(:), rlon(:)
      REAL, ALLOCATABLE :: r3(:,:,:)

! Attach to the grid

      gdid = gdattach(gdfid, gridName)
      IF (gdId == -1) THEN
         msr = 'Failed to attach to grid ' // gridName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get dimension names, sizes

      size = gddiminfo(gdid, DIMX_NAME)
      IF (size == -1) THEN
         msr = SZ_ERR // DIMX_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlon = size

      size = gddiminfo(gdid, DIMY_NAME)
      IF (size == -1) THEN
         msr = SZ_ERR // DIMY_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlat = size
      
      size = gddiminfo(gdid, DIMZ_NAME)
      IF (size == -1) THEN
         msr = SZ_ERR // DIMZ_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlev = size

! Allocate the local REAL variables & the output structure pointers

      ALLOCATE(rp(nlev), rlat(nlat), rlon(nlon), r3(nlev,nlat,nlon), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local REAL variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL AllocateL3DM(nlev, nlat, nlon, l3dm)

! Read data into the output structure 

      l3dm%name = gridName

      start = 0
      stride = 1
      edge(1) = l3dm%nLevels
      edge(2) = l3dm%nLats
      edge(3) = l3dm%nLons

      status = gdrdfld(gdid, GEO_FIELD3, start(1), stride(1), stride(1), &
                       l3dm%time)
      IF (status == -1) THEN
         msr = RDGD_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdrdfld(gdid, GEO_FIELD9, start(1), stride(1), edge(1), rp)
      IF (status == -1) THEN
         msr = RDGD_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dm%pressure = DBLE(rp)

      status = gdrdfld(gdid, GEO_FIELD1, start(2), stride(2), edge(2), rlat)
      IF (status == -1) THEN
         msr = RDGD_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dm%latitude = DBLE(rlat)

      status = gdrdfld(gdid, GEO_FIELD2, start(3), stride(3), edge(3), rlon)
      IF (status == -1) THEN
         msr = RDGD_ERR // GEO_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dm%longitude = DBLE(rlon)

      status = gdrdfld(gdid, DATA_FIELDV, start, stride, edge, r3)
      IF (status == -1) THEN
         msr = RDGD_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dm%l3dmValue = DBLE(r3)

      status = gdrdfld(gdid, DATA_FIELDP, start, stride, edge, r3)
      IF (status == -1) THEN
         msr = RDGD_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dm%l3dmPrecision = DBLE(r3)

! Detach from and close the grid interface

      status = gddetach(gdid)
      IF (status == -1) THEN
         msr = GD_ERR // TRIM(gridName) // ' after reading.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Deallocate the local variables

      DEALLOCATE (rp, rlat, rlon, r3, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local REAL variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE ReadL3DMData
!-----------------------------

!-----------------------------------------------------
   SUBROUTINE WriteMetaL3DM (pcf, l3cf, files, anText)
!-----------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dm file, and annotates it with the
! PCF.

! Arguments

      TYPE( L3CFProd_T ), INTENT(IN) :: l3cf

      TYPE (OutputFiles_T), INTENT(IN) :: files

      TYPE( PCFData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write
      INTEGER, EXTERNAL :: pgs_td_taiToUTC, gdinqgrid

! Variables

      CHARACTER (LEN=1) :: cNum
      CHARACTER (LEN=45) :: attrName
      CHARACTER (LEN=132) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=CCSDS_LEN) :: timeA
      CHARACTER (LEN=FileNameLen) :: sval
      CHARACTER (LEN=GridNameLen) :: gridName
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: hdfReturn, i, j, indx, len, numGrids, result, sdid, returnStatus

      REAL(r8) :: dval

! Initialize the MCF file

      result = pgs_met_init(l3cf%mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! For each l3dm file successfully created,

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
         CALL ExpandFileTemplate('$cycle', sval, cycle=pcf%cycle)
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! MeasuredParameterContainer -- find the number of grids in the file

         numGrids = gdinqgrid(files%name(i), list, len)
         IF (numGrids == -1) THEN
            msr = 'No grids found in file ' // TRIM( files%name(i) ) // &
                  ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! For each grid in the file

         DO j = 1, numGrids

! Extract its name

            indx = INDEX(list, ',')
            IF (indx /= 0) THEN
               gridName = list(:indx-1)
               list = list(indx+1:)
            ELSE
               gridName = list
            ENDIF

! Append a class suffix to ParameterName, and write the grid name as its value

            WRITE( cNum, '(I1)' ) j

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
         returnStatus = pgs_td_taiToUTC(l3cf%timeD(1), timeA)
         IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                                     ModuleName, TAI2A_ERR)
         sval = timeA(12:26)
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
   END SUBROUTINE WriteMetaL3DM
!------------------------------

!--------------------------------------------------
   SUBROUTINE AllocateL3DM (nlev, nlat, nlon, l3dm)
!--------------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3DMData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat, nlon

      TYPE( L3DMData_T ), INTENT(INOUT) :: l3dm

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3dm%nLevels = nlev
      l3dm%nLats = nlat
      l3dm%nLons = nlon

! Horizontal geolocation fields

      ALLOCATE(l3dm%latitude(l3dm%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMData_T latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dm%longitude(l3dm%nLons), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMData_T longitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Vertical geolocation field

      ALLOCATE(l3dm%pressure(l3dm%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMData_T pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      ALLOCATE(l3dm%l3dmValue(l3dm%nLevels,l3dm%nLats,l3dm%nLons), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMData_T value pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dm%l3dmPrecision(l3dm%nLevels,l3dm%nLats,l3dm%nLons),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMData_T precision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE AllocateL3DM
!-----------------------------

!----------------------------------
   SUBROUTINE DeallocateL3DM (l3dm)
!----------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L3DMData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3DMData_T ), INTENT(INOUT) :: l3dm

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation fields

      IF ( ASSOCIATED(l3dm%latitude) ) THEN
         DEALLOCATE (l3dm%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dm latitude pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dm%longitude) ) THEN
         DEALLOCATE (l3dm%longitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dm longitude pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3dm%pressure) ) THEN
         DEALLOCATE (l3dm%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dm pressure pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3dm%l3dmValue) ) THEN
         DEALLOCATE (l3dm%l3dmValue, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dmValue pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dm%l3dmPrecision) ) THEN
         DEALLOCATE (l3dm%l3dmPrecision, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dmPrecision'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3DM
!-------------------------------

!-----------------------------------------
   SUBROUTINE DestroyL3DMDatabase (l3dmdb)
!-----------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal structures of an l3dm database, and
! then database itself


! Arguments

      TYPE (L3DMData_T), DIMENSION(:), POINTER :: l3dmdb

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i

! Check the status of the input pointer

      IF ( ASSOCIATED(l3dmdb) ) THEN

! If it's associated, then deallocate the internal structures

         DO i = 1, SIZE(l3dmdb)
            CALL DeallocateL3DM( l3dmdb(i) )
         ENDDO

! Deallocate the database itself

         DEALLOCATE (l3dmdb, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dm database'
            CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

!------------------------------------
   END SUBROUTINE DestroyL3DMDatabase
!------------------------------------

!==================
END MODULE L3DMData
!==================

!# $Log$
!# Revision 1.12  2001/05/04 18:30:01  nakamura
!# Removed L3DMFiles_T and used generic OutputFiles_T.
!#
!# Revision 1.11  2001/04/24 19:38:23  nakamura
!# Removed references to private L2 parameters.
!#
!# Revision 1.10  2001/03/27 19:28:15  nakamura
!# Moved some parameters to MLSL3Common; updated metadata; fixed err checks on deallocate.
!#
!# Revision 1.9  2001/02/21 20:57:12  nakamura
!# Changed MLSPCF to MLSPCF3; made some parameters global; added ReadL3DMData; changed InputPointer.
!#
!# Revision 1.8  2001/02/09 19:18:05  nakamura
!# Moved DIMT to MLSL3Common; changed LocalGranuleID to file - .dat
!#
!# Revision 1.7  2001/01/16 17:44:00  nakamura
!# Made lowrgt corner of the grid variable; updated WriteMetaL3DM for new MCF and added writing of annotation.
!#
!# Revision 1.6  2000/12/29 20:50:33  nakamura
!# Moved global parameters to MLSL3Common; added L3DMFiles_T; switched to one-product/all-days paradigm.
!#
!# Revision 1.5  2000/12/07 19:34:42  nakamura
!# Added ONLY to USE L2GPData; replaced local error msgs with MLSMSG_DeAllocate.
!#
!# Revision 1.4  2000/11/15 21:00:26  nakamura
!# Added parameter GridNameLen.
!#
!# Revision 1.3  2000/10/24 19:21:08  nakamura
!# Removed dependence on OutputL2GP.
!#
!# Revision 1.2  2000/10/17 20:17:09  nakamura
!# Added parameters used globally by L3; updated WriteMetaL3DM for new input.
!#
!# Revision 1.1  2000/10/05 19:11:07  nakamura
!# Module for the L3DM data type.
!#
