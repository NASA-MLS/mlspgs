
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DMData
!===============================================================================

   USE Hdf
   USE L2GPData, ONLY: L2GPNameLen, GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, &
                       GEO_FIELD9, HDFE_NOMERGE
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF
   USE MLSStrings
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

! Definitions -- L3DMData_T
!                L3DMFiles_T 
! Subroutines -- ConvertDeg2DMS
!                OutputGrids
!                WriteMetaL3DM
!                AllocateL3DM
!                DeallocateL3DM
!                DestroyL3DMDatabase

! Remarks:  This module contains the definition of the L3DMData type, as well
!           as any routines pertaining to it.

! Parameters

! This data type is used to store the l3 daily map data.

   TYPE L3DMData_T

     CHARACTER (LEN=L2GPNameLen) :: name	! name for the output quantity

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

!    REAL(r8), DIMENSION(:,:), POINTER :: l3dmColumn	! Column data
	! dimensioned as (nLats, nLons)

   END TYPE L3DMData_T

! This data type is used to store the names of l3dm files actually created.

   TYPE L3DMFiles_T

     INTEGER :: nFiles		! number of distinct l3dm output files created

     CHARACTER (LEN=FileNameLen) :: name(maxwindow)
	! array of names of the created files

     CHARACTER (LEN=8) :: date(maxWindow)	! CCSDS B format date of files

   END TYPE L3DMFiles_T

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

      TYPE (L3DMFiles_T), INTENT(INOUT) :: l3dmFiles
   
! Parameters

      CHARACTER (LEN=*), PARAMETER :: DIMX_NAME = 'XDim'
      CHARACTER (LEN=*), PARAMETER :: DIMY_NAME = 'YDim'
      CHARACTER (LEN=*), PARAMETER :: DIMZ_NAME = 'ZDim'
      CHARACTER (LEN=*), PARAMETER :: DIMT_NAME = 'TDim'
      CHARACTER (LEN=*), PARAMETER :: DIMXYZ_NAME = 'ZDim,YDim,XDim'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3dmValue'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3dmPrecision'

      CHARACTER (LEN=*), PARAMETER :: GD_ERR = 'Failed to detach from grid '

      INTEGER, PARAMETER :: GCTP_GEO = 0

! Functions

      INTEGER, EXTERNAL :: gdattach, gdclose, gdcreate, gddefdim, gddeffld
      INTEGER, EXTERNAL :: gddefproj, gddetach, gdopen, gdwrfld

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=FileNameLen) :: physicalFilename
      CHARACTER (LEN=480) :: msr

      INTEGER :: gdfID, gdId, i, match, status
      INTEGER :: start(3), stride(3), edge(3)

      REAL :: maxLat, minLat

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

! Check whether the name is distinct; if so, save it in L3DMFiles_T

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

! Set up the grid.  The region is bounded by 180.0W to 180.0E longitude &
! varying latitude.  Grid into 91 bins along the x-axis, by nLats bins along
! the y-axis (4x2 bins).  Upper Left & Lower Right corners in DDDMMMSSS.ss.

         uplft(1) = -180000000.00
         lowrgt(1) = 180000000.00

         projparm = 0.0

! Find boundaries of measured latitude

         maxLat = MAXVAL( l3dmData(i)%latitude )
         minLat = MINVAL( l3dmData(i)%latitude )

! Convert to "packed degree format"

         CALL ConvertDeg2DMS(maxLat, uplft(2))
         CALL ConvertDeg2DMS(minLat, lowrgt(2))

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

!----------------------------------------------
   SUBROUTINE WriteMetaL3DM (files, mlspcf_mcf)
!----------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dm file.

! Arguments

      TYPE (L3DMFiles_T), INTENT(IN) :: files

      INTEGER, INTENT(IN) :: mlspcf_mcf

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_s, pgs_met_write

! Variables

      CHARACTER (LEN=10) :: date
      CHARACTER (LEN=15) :: time
      CHARACTER (LEN=21) :: sType, zoneID
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: hdfReturn, i, result, sdid

      REAL(r8) :: maxLat, maxLon, minLat, minLon

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! Set nominal values for the overall file

      sType = 'Horizontal & Vertical'
      time = '12:00:00.000000'
      zoneID = 'Other Grid System'
      minLon = -180.0
      maxLat = 82.0
      maxLon = 180.0
      minLat = -82.0

! For each l3dm file successfully created,

      DO i = 1, files%nFiles

! Open the HDF file and initialize the SD interface

         sdid = sfstart(files%name(i), DFACC_WRITE)
         IF (sdid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                     &open the HDF file for metadata writing.')
 
! Set PGE values

         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "GranuleSpatialDomainType", sType)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         date = files%date(i)
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "RangeBeginningDate", date)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "RangeBeginningTime", time)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "RangeEndingDate", date)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "RangeEndingTime", time)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                   "ZoneIdentifier", zoneID)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                   "WestBoundingCoordinate", minLon)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                   "NorthBoundingCoordinate", maxLat)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                   "EastBoundingCoordinate", maxLon)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                   "SouthBoundingCoordinate", minLat)
         IF (result /= PGS_S_SUCCESS) THEN
            call Pgs_smf_getMsg(result, mnemonic, msg)
            msr = mnemonic // ':  ' // msg
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

      IF ( ASSOCIATED(l3dm%latitude) ) DEALLOCATE (l3dm%latitude, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dm latitude pointer'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3dm%longitude) ) DEALLOCATE (l3dm%longitude, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dm longitude pointer'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3dm%pressure) ) DEALLOCATE (l3dm%pressure, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dm pressure pointer'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3dm%l3dmValue) ) DEALLOCATE (l3dm%l3dmValue, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dmValue pointer'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3dm%l3dmPrecision) ) DEALLOCATE (l3dm%l3dmPrecision, &
                                                        STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3dmPrecision'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
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
         END DO

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
