
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DMData
!===============================================================================

   USE SDPToolkit
   USE MLSMessageModule
   USE MLSCommon
   USE L2GPData, ONLY: L2GPNameLen, GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, &
                       GEO_FIELD9, HDFE_NOMERGE
   USE Hdf
   USE MLSPCF
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
! Subroutines -- ConvertDeg2DMS
!                OutputGrids
!                WriteMetaL3DM
!                AllocateL3DM
!                DeallocateL3DM
!                DestroyL3DMDatabase

! Remarks:  This module contains the definition of the L3DMData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DAT_ERR = 'Failed to define data field '
   CHARACTER (LEN=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
   CHARACTER (LEN=*), PARAMETER :: GEO_ERR = 'Failed to define geolocation &
                                             &field '
   CHARACTER (LEN=*), PARAMETER :: TAI2A_ERR = 'Error converting time from &
                                               &TAI to UTC.'
   CHARACTER (LEN=*), PARAMETER :: WR_ERR = 'Failed to write field '

   INTEGER, PARAMETER :: CCSDS_LEN = 27
   INTEGER, PARAMETER :: CCSDSB_LEN = 25
   INTEGER, PARAMETER :: INVENTORYMETADATA = 2
   INTEGER, PARAMETER :: GridNameLen = 64
   INTEGER, PARAMETER :: maxNumGrids = 100
   INTEGER, PARAMETER :: maxWindow = 30

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

!--------------------------------------------------------------------
   SUBROUTINE OutputGrids(physicalFilename, numGrids, indx, l3dmData)
!--------------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the grid portion of the l3dm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

      INTEGER, INTENT(IN) :: numGrids 

      INTEGER, INTENT(IN) :: indx(numGrids)

      TYPE (L3DMData_T), INTENT(IN) :: l3dmData(:)
   
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

      CHARACTER (LEN=480) :: msr

      INTEGER :: gdfID, gdId, i, status
      INTEGER :: start(3), stride(3), edge(3)

      REAL :: maxLat, minLat

      REAL(r8) :: uplft(2), lowrgt(2)
      REAL(r8) :: projparm(13)

! Open the output file

      gdfID = gdopen(physicalFilename, DFACC_CREATE)
      IF (gdfID == -1) THEN
         msr = MLSMSG_Fileopen // physicalFilename
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Set up the grids.  The region is bounded by 180.0W to 180.0E longitude &
! varying latitude.  Grid into 91 bins along the x-axis, by nLats bins along
! the y-axis (4x2 bins).  Upper Left & Lower Right corners in DDDMMMSSS.ss.

      uplft(1) = -180000000.00
      lowrgt(1) = 180000000.00

      projparm = 0.0

! For each grid belonging in this file,

      DO i = 1, numGrids

! Find boundaries of measured latitude

         maxLat = MAXVAL( l3dmData(indx(i))%latitude )
         minLat = MINVAL( l3dmData(indx(i))%latitude )

! Convert to "packed degree format"

         CALL ConvertDeg2DMS(maxLat, uplft(2))
         CALL ConvertDeg2DMS(minLat, lowrgt(2))

! Create the grids

         gdId = gdcreate(gdfID, l3dmData(indx(i))%name, &
             l3dmData(indx(i))%nLons, l3dmData(indx(i))%nLats, uplft, lowrgt)
         IF (gdId == -1) THEN
            msr = 'Failed to create grid ' // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the dimensions

         status = gddefdim(gdId, DIMX_NAME, l3dmData(indx(i))%nLons)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMX_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddefdim(gdId, DIMY_NAME, l3dmData(indx(i))%nLats)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMY_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gddefdim(gdId, DIMZ_NAME, l3dmData(indx(i))%nLevels)
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
            msr = 'Failed to define projection for grid ' // &
                   l3dmData(indx(i))%name
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

      ENDDO

      status = gdclose(gdfID)
      IF (status /= 0) THEN
         msr = 'Failed to close file ' // TRIM(physicalFilename) // &
               ' after definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-open file for writing

      gdfID = gdopen(physicalFilename, DFACC_RDWR)
      IF (gdfID == -1) THEN
         msr = MLSMSG_Fileopen // TRIM(physicalFilename) // ' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Initialize values outside of grid loop

      start = 0
      stride = 1

! For each grid belonging in the file,

      DO i = 1, numGrids

         edge(1) = l3dmData(indx(i))%nLevels
         edge(2) = l3dmData(indx(i))%nLats
         edge(3) = l3dmData(indx(i))%nLons

! Re-attach to the grid for writing

         gdId = gdattach(gdfID, l3dmData(indx(i))%name)
         IF (gdId == -1) THEN
            msr = 'Failed to attach to grid ' // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write to fields

         status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
                          l3dmData(indx(i))%time)
         IF (status /=0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' &
                   // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1), &
                           REAL(l3dmData(indx(i))%pressure) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' &
                  // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
                           REAL(l3dmData(indx(i))%latitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' &
                  // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
                           REAL(l3dmData(indx(i))%longitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' &
                   // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, DATA_FIELDV, start, stride, edge, &
                           REAL(l3dmData(indx(i))%l3dmValue) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDV // ' to grid ' &
                  // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = gdwrfld( gdId, DATA_FIELDP, start, stride, edge, &
                           REAL(l3dmData(indx(i))%l3dmPrecision) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDP // ' to grid ' &
                  // l3dmData(indx(i))%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from each grid after writing

         status = gddetach(gdId)
         IF (status /= 0) THEN
            msr = GD_ERR // TRIM( l3dmData(indx(i))%name ) // &
                  ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

! Close the file after writing

      status = gdclose(gdfID)
      IF (status /= 0) THEN
         msr = 'Failed to close file ' // TRIM(physicalFilename) // &
               ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!----------------------------
   END SUBROUTINE OutputGrids
!----------------------------

!----------------------------------------------------------------------------
   SUBROUTINE WriteMetaL3DM (l3File, mlspcf_mcf, numGrids, indx, l3dm, timeA)
!----------------------------------------------------------------------------

! Brief description of subroutine
! This routine writes the metadata for an l3dm file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: l3File

      INTEGER, INTENT(IN) :: mlspcf_mcf

      INTEGER, INTENT(IN) :: numGrids

      INTEGER, INTENT(IN) :: indx(:)

      TYPE (L3DMData_T), INTENT(IN) :: l3dm(:)

      CHARACTER (LEN=CCSDS_LEN), INTENT(IN) :: timeA 

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_s, pgs_met_write

! Variables

      CHARACTER (LEN=10) :: date
      CHARACTER (LEN=15) :: time
      CHARACTER (LEN=21) :: sval
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: allGrids, hdfReturn, i, result, sdid

      REAL(r8) :: dval, maxGrid, maxLat, minGrid, minLat

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                           'Initialization error.  See LogStatus for details.')

! Open the HDF file and initialize the SD interface

      sdid = sfstart(l3File, DFACC_WRITE)
      IF (sdid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                     &open the HDF file for metadata writing.')
 
! Find the time, max & min lat over all grids for this product

      allGrids = SIZE(l3dm)

      maxLat = 0.0
      minLat = 0.0

      DO i = 1, numGrids
         maxGrid = MAXVAL(l3dm(indx(i))%latitude)
         minGrid = MINVAL(l3dm(indx(i))%latitude)
         maxLat = MAX(maxGrid,maxLat)
         minLat = MIN(minGrid,minLat)
      ENDDO

      date = timeA(1:10)
      time = timeA(12:26)

! Set PGE values

      sval = 'Horizontal & Vertical'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                "GranuleSpatialDomainType", sval)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                "RangeBeginningDate", date)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                "RangeBeginningTime", time)

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", date)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingTime", time)

      sval = 'Other Grid System'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "ZoneIdentifier", &
                                 sval)

      dval = -180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "WestBoundingCoordinate", dval)
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "NorthBoundingCoordinate", maxLat)
      dval = 180.0
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "EastBoundingCoordinate", dval)
      result = pgs_met_setAttr_d(groups(INVENTORYMETADATA), &
                                "SouthBoundingCoordinate", minLat)

      IF (result /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(result, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
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
      IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Error &
                                    &closing HDF file after writing metadata.')

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
!#
