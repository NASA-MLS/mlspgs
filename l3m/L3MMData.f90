! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE L3MMData
!==============================================================================

   USE HDF, ONLY: DFACC_RDWR, DFACC_WRITE, DFNT_FLOAT64, DFNT_FLOAT32, &
        & DFNT_INT32, DFNT_CHAR8
   USE Intrinsic, ONLY: l_hdfeos, l_grid
   USE L3DMData, ONLY: ConvertDeg2DMS
   USE MLSCommon, ONLY: r8
   USE MLSFiles, ONLY: MLS_SFEND, MLS_SFSTART, HDFVERSION_4, HDFVERSION_5
   USE MLSL3Common, ONLY: INVENTORYMETADATA, DATE_LEN, HDFE_NOMERGE, MIN_MAX, &
        & GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, GEO_FIELD9, DIM_ERR, DIM_NAME2, &
        & maxWindow, GridNameLen, DATE_LEN, GCTP_GEO, DIMR_NAME, &
        & DIMX_NAME, DIMY_NAME, DIMZ_NAME, DIMT_NAME, DIM_ERR, DIMXYZ_NAME, &
        & DG_FIELD2, GEO_ERR, DAT_ERR, GD_ERR, SW_ERR, WR_ERR, &
        & METAWR_ERR
   USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
        & MLSMSG_Info, MLSMSG_Allocate, MLSMSG_Deallocate, MLSMSG_WARNING
   USE Time_M, only: Time_Now
   ! USE SWAPI
   ! USE HE5_SWAPI

   IMPLICIT NONE
   private
   PUBLIC :: L3MMData_T, OutputMMGrids, OutputMMDiags, WriteMetaL3MM, &
     & AllocateL3MM, DeallocateL3MM

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

     INTEGER :: misDays(maxWindow)
	! Missing days, dimensioned (nDays)

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:,:), POINTER :: l3mmValue	  ! Field value
     REAL(r8), DIMENSION(:,:,:), POINTER :: l3mmPrecision ! Field precision
	! dimensioned as (nLevels, nLats, nLons)

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)
     REAL(r8), DIMENSION(:), POINTER :: longitude	! dimensioned (nLons)

     REAL(r8) :: startTime	! start time of L2 data used in the analysis
     REAL(r8) :: endTime 	! end time of L2 data used in the analysis

     ! Now the diagnostic fields:

     INTEGER, DIMENSION(:), POINTER :: perMisPoints	
     ! Missing points (percentage)
	! dimensioned (nLevels)

     INTEGER :: nLevels			! Total number of surfaces
     INTEGER :: nLats			! Total number of latitudes
     INTEGER :: nLons			! Total number of longitudes
     INTEGER :: nMisDays		! Number of missing days in the month

   END TYPE L3MMData_T

CONTAINS

!-----------------------------------------------------------------
  SUBROUTINE OutputMMGrids (physicalFilename, l3mm, creationFlag, &
       & hdfVersion)
!-----------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the grid portion of l3mm files.

! Arguments

    CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

    TYPE (L3MMData_T), INTENT(IN) :: l3mm

    LOGICAL, INTENT(INOUT) :: creationFlag
   
    INTEGER, INTENT(IN) :: hdfVersion

    IF (hdfVersion==HDFVERSION_4) THEN 
       CALL OutputMMGrids_HE2 (physicalFilename, l3mm, creationFlag)
    ELSE IF (hdfVersion==HDFVERSION_5) THEN 
       CALL OutputMMGrids_HE5 (physicalFilename, l3mm, creationFlag)
    ENDIF

  END SUBROUTINE OutputMMGrids

!-----------------------------------------------------------------
  SUBROUTINE OutputMMGrids_HE2 (physicalFilename, l3mm, creationFlag)
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

    INTEGER, EXTERNAL :: gdattach, gdclose, gdcreate, gddefdim, gddeffld, &
         & gddefproj, gddetach, gdopen, gdwrfld

! Variables

    CHARACTER (LEN=480) :: msr

    REAL(r8), dimension(:), pointer :: times
    REAL(r8) :: projparm(13)
    REAL(r8) :: uplft(2), lowrgt(2)

    REAL :: maxLat, maxLon, minLat, minLon
    REAL, parameter :: WAIT_TO_REOPEN = 2.0

    INTEGER :: start(3), stride(3), edge(3)
    INTEGER :: gdfID, gdId, status, err
    INTEGER, parameter :: NUMBER_BEATS = 20000
      
! Open the output file

    gdfID = gdopen(physicalFilename, DFACC_RDWR)
    IF (gdfID == -1) THEN
       msr = MLSMSG_Fileopen // trim(physicalFilename)
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Set up the grid.  The region is bounded by 180.0W to 176.0E longitude &
! varying latitude.  Grid into 90 bins along the x-axis, by nLats bins along
! the y-axis (4x2 bins).  Upper Left & Lower Right corners in DDDMMMSSS.ss.

    projparm = 0.0

! Find boundaries of measured latitude, upper bound of measure longitude

    maxLat = MAXVAL( l3mm%latitude  )
    minLat = MINVAL( l3mm%latitude  )
    maxLon = MAXVAL( l3mm%longitude )
    minLon = MINVAL( l3mm%longitude )

! Convert to "packed degree format"

    CALL ConvertDeg2DMS(maxLat, uplft(2) )
    CALL ConvertDeg2DMS(minLat, lowrgt(2))
    CALL ConvertDeg2DMS(maxLon, lowrgt(1))
    CALL ConvertDeg2DMS(minLon, uplft(1) )

! Create the grid

    gdId = gdcreate(gdfID, l3mm%name, l3mm%nLons, l3mm%nLats, uplft, lowrgt)
    IF (gdId == -1) THEN
       msr = 'Failed to create grid ' // trim(l3mm%name)
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

    status = gddeffld(gdId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64,HDFE_NOMERGE)
    IF (status /= 0) THEN
       msr = GEO_ERR // GEO_FIELD3
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = gddeffld(gdId, GEO_FIELD9, DIMZ_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
    IF (status /= 0) THEN
       msr = GEO_ERR // GEO_FIELD9
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = gddeffld(gdId, GEO_FIELD1, DIMY_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
    IF (status /= 0) THEN
       msr = GEO_ERR // GEO_FIELD1
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = gddeffld(gdId, GEO_FIELD2, DIMX_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
    IF (status /= 0) THEN
       msr = GEO_ERR // GEO_FIELD2
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Define the "data" fields

    status = gddeffld(gdId, DATA_FIELDMV, DIMXYZ_NAME, DFNT_FLOAT32, & 
         & HDFE_NOMERGE)
    IF (status /= 0) THEN
       msr = DAT_ERR // DATA_FIELDMV
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = gddeffld(gdId, DATA_FIELDMP, DIMXYZ_NAME, DFNT_FLOAT32, & 
         & HDFE_NOMERGE)
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
            & ' after definition.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Delay reopening for a discreet and discrete time
! the smaller of WAIT_TO_REOPEN (in sec) or NUMBER_BEATS (in passes)
!      call time_now(t1)
!      t2 = t1
!      do beat = 1, NUMBER_BEATS
!        call time_now(t2)
!        if ( (t2 - t1) > WAIT_TO_REOPEN ) exit
!      enddo

! Re-open the file for writing

    gdfID = gdopen(physicalFilename, DFACC_RDWR)
    IF (gdfID == -1) THEN
       msr = MLSMSG_Fileopen // TRIM(physicalFilename) // ' for writing.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Re-attach to the grid for writing

    gdId = gdattach(gdfID, l3mm%name)
    IF (gdId == -1) THEN
       msr = 'Failed to attach to grid ' // trim(l3mm%name)
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Write to fields

    start = 0
    stride = 1
    edge(1) = l3mm%nLevels
    edge(2) = l3mm%nLats
    edge(3) = l3mm%nLons

    allocate(times(l3mm%nLevels), STAT=err)
    IF ( err /= 0 ) THEN
       msr = MLSMSG_Allocate // ' times pointer.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    times(1:l3mm%nLevels) = l3mm%startTime

! Start Time

!      status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
!                       l3mm%startTime)
    status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
         & times)
    IF (status /= 0) THEN
       print *, 'gdId ', gdId
       print *, 'start ', start
       print *, 'stride ', stride
       print *, 'edge ', edge
       print *, 'start times ', times
       msr = 'Failed to write startTime to grid ' // trim(l3mm%name)
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! End Time

    times(1:l3mm%nLevels) = l3mm%endTime
    status = gdwrfld(gdId, GEO_FIELD3, stride(1), stride(1), stride(1), &
         & times)
!      status = gdwrfld(gdId, GEO_FIELD3, stride(1), stride(1), stride(1), &
!                       l3mm%endTime)
    IF (status /= 0) THEN
       msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' & 
            & // trim(l3mm%name)
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

! Pressure, latitude, & longitude

    status = gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1), &
         & REAL(l3mm%pressure) )
      IF (status /= 0) THEN
         msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' & 
              & // trim(l3mm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(l3mm%latitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
           & REAL(l3mm%longitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = gdwrfld( gdId, DATA_FIELDMV, start, stride, edge, &
           & REAL(l3mm%l3mmValue) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMV // ' to grid ' // &
               & trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = gdwrfld( gdId, DATA_FIELDMP, start, stride, edge, &
           & REAL(l3mm%l3mmPrecision) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMP // ' to grid ' // &
               & trim(l3mm%name)
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
               & ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Grid ' // TRIM(l3mm%name) // ' successfully written to file ' // &
           & trim(physicalFilename)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      creationFlag = .TRUE.

      deallocate(times, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Deallocate // ' times pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!------------------------------
   END SUBROUTINE OutputMMGrids_HE2
!------------------------------

!-----------------------------------------------------------------
   SUBROUTINE OutputMMGrids_HE5 (physicalFilename, l3mm, creationFlag)
!-----------------------------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC

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

      INTEGER, EXTERNAL :: he5_gdattach, he5_gdclose, he5_gdcreate, & 
           & he5_gddefdim, he5_gddeffld, he5_gddefproj, he5_gddetach, & 
           & he5_gdopen, he5_gdwrfld

! Variables

      CHARACTER (LEN=480) :: msr

      REAL(r8), dimension(:), pointer :: times
      REAL(r8) :: projparm(13)
      REAL(r8) :: uplft(2), lowrgt(2)

      REAL :: maxLat, maxLon, minLat, minLon
      REAL, parameter :: WAIT_TO_REOPEN = 2.0

      INTEGER :: start(3), stride(3), edge(3)
      INTEGER :: gdfID, gdId, status, err
      INTEGER, parameter :: NUMBER_BEATS = 20000
      
! Open the output file. If first time, create the file. Otherwise just open and 
! write other species to the files

      IF (.not. creationFlag) THEN
         gdfID = he5_gdopen(trim(physicalFilename), HE5F_ACC_TRUNC)
      ELSE 
         gdfID = he5_gdopen(trim(physicalFilename), HE5F_ACC_RDWR)
      ENDIF

      IF (gdfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Set up the grid.  The region is bounded by 180.0W to 176.0E longitude &
! varying latitude.  Grid into 90 bins along the x-axis, by nLats bins along
! the y-axis (4x2 bins).  Upper Left & Lower Right corners in DDDMMMSSS.ss.

      projparm = 0.0

! Find boundaries of measured latitude, upper bound of measure longitude

      maxLat = MAXVAL( l3mm%latitude  )
      minLat = MINVAL( l3mm%latitude  )
      maxLon = MAXVAL( l3mm%longitude )
      minLon = MINVAL( l3mm%longitude )

! Convert to "packed degree format"

      CALL ConvertDeg2DMS(maxLat, uplft(2) )
      CALL ConvertDeg2DMS(minLat, lowrgt(2))
      CALL ConvertDeg2DMS(maxLon, lowrgt(1))
      CALL ConvertDeg2DMS(minLon, uplft(1) )

! Create the grid

      gdId = he5_gdcreate(gdfID, trim(l3mm%name), l3mm%nLons, l3mm%nLats, & 
           & uplft, lowrgt)
      IF (gdId == -1) THEN
          msr = 'Failed to create grid ' // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the dimensions

      status = he5_gddefdim(gdId, DIMX_NAME, l3mm%nLons)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMX_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddefdim(gdId, DIMY_NAME, l3mm%nLats)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMY_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddefdim(gdId, DIMZ_NAME, l3mm%nLevels)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMZ_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddefdim(gdId, DIMT_NAME, MIN_MAX)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMT_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define a "Geographic projection," using defaults in all unneeded fields

      status = he5_gddefproj(gdId, GCTP_GEO, 0, 0, projparm)
      IF (status /= 0) THEN
          msr = 'Failed to define projection for grid ' // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "geolocation" fields

      status = he5_gddeffld(gdId, GEO_FIELD3, DIMT_NAME, "", & 
           & HE5T_NATIVE_DOUBLE, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddeffld(gdId, GEO_FIELD9, DIMZ_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddeffld(gdId, GEO_FIELD1, DIMY_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddeffld(gdId, GEO_FIELD2, DIMX_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "data" fields

      status = he5_gddeffld(gdId, DATA_FIELDMV, DIMXYZ_NAME, "",& 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDMV
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gddeffld(gdId, DATA_FIELDMP, DIMXYZ_NAME, "",& 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDMP
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Delay reopening for a discreet and discrete time
! the smaller of WAIT_TO_REOPEN (in sec) or NUMBER_BEATS (in passes)
!      call time_now(t1)
!      t2 = t1
!      do beat = 1, NUMBER_BEATS
!        call time_now(t2)
!        if ( (t2 - t1) > WAIT_TO_REOPEN ) exit
!      enddo

! Write to fields

      start = 0
      stride = 1
      edge(1) = l3mm%nLevels
      edge(2) = l3mm%nLats
      edge(3) = l3mm%nLons

      allocate(times(l3mm%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' times pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      times(1:l3mm%nLevels) = l3mm%startTime

! Start Time

      status = he5_gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), stride(1), &
           & times)
      IF (status /= 0) THEN
          print *, 'gdId ', gdId
          print *, 'start ', start
          print *, 'stride ', stride
          print *, 'edge ', edge
          print *, 'start times ', times
          msr = 'Failed to write startTime to grid ' // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! End Time

      times(1:l3mm%nLevels) = l3mm%endTime
      status = he5_gdwrfld(gdId, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & times)
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Pressure, latitude, & longitude

      status = he5_gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(l3mm%pressure) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(l3mm%latitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
           & REAL(l3mm%longitude) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' & 
               & // trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = he5_gdwrfld( gdId, DATA_FIELDMV, start, stride, edge, &
           & REAL(l3mm%l3mmValue) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMV // ' to grid ' // &
               & trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_gdwrfld( gdId, DATA_FIELDMP, start, stride, edge, &
           & REAL(l3mm%l3mmPrecision) )
      IF (status /= 0) THEN
          msr = 'Failed to write field ' //  DATA_FIELDMP // ' to grid ' // &
               & trim(l3mm%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the grid after writing

      status = he5_gddetach(gdId)
      IF (status /= 0) THEN
          msr = GD_ERR // TRIM( l3mm%name ) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Close the file after writing

      status = he5_gdclose(gdfID)
      IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
               & ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Grid ' // TRIM(l3mm%name) // ' successfully written to file ' // &
           & trim(physicalFilename)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      creationFlag = .TRUE.

      deallocate(times, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Deallocate // ' times pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!------------------------------
   END SUBROUTINE OutputMMGrids_HE5
!------------------------------

!------------------------------------------------
   SUBROUTINE OutputMMDiags(physicalFilename, mm, hdfVersion)
!------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic portion of l3mm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

      TYPE (L3MMData_T), INTENT(IN) :: mm

      INTEGER, INTENT(IN) :: hdfVersion

      IF (hdfVersion==HDFVERSION_4) THEN 
         CALL OutputMMDiags_HE2(physicalFilename, mm)
      ELSE IF (hdfVersion==HDFVERSION_5) THEN
         CALL OutputMMDiags_HE5(physicalFilename, mm)
      ENDIF

    END SUBROUTINE OutputMMDiags

!------------------------------------------------
   SUBROUTINE OutputMMDiags_HE2(physicalFilename, mm)
!------------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic portion of l3mm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename

      TYPE (L3MMData_T), INTENT(IN) :: mm

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DIMD_NAME   = 'lDate'
      CHARACTER (LEN=*), PARAMETER :: DIMMIS_NAME = 'nMisDays'
      CHARACTER (LEN=*), PARAMETER :: DIMDM_NAME  = 'lDate,nMisDays'
      CHARACTER (LEN=*), PARAMETER :: MIS_FIELD   = 'MissingDays'

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, & 
           & swdefgfld, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      REAL, PARAMETER :: WAIT_TO_REOPEN = 2.0

      INTEGER :: edge(2), start(2), stride(2)
      INTEGER :: lsize, nsize, status, swfID, swId
      INTEGER, PARAMETER :: NUMBER_BEATS = 20000

! Open the output file

      swfID = swopen(physicalFilename, DFACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create the swath

      dgName = TRIM(mm%name) // 'Diagnostics'
      swId = swcreate(swfID, dgName)
      IF (swId == -1) THEN
          msr = 'Failed to create swath ' // trim(dgName)
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

      status = swdefgfld(swID, GEO_FIELD3, DIMR_NAME, DFNT_FLOAT64, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swId, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
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

! Delay reopening for a discreet and discrete time
! the smaller of WAIT_TO_REOPEN (in sec) or NUMBER_BEATS (in passes)
!      call time_now(t1)
!      t2 = t1
!      do beat = 1, NUMBER_BEATS
!        call time_now(t2)
!        if ( (t2 - t1) > WAIT_TO_REOPEN ) exit
!      enddo

! Re-attach to the swath for writing

      swId = swattach(swfID, dgName)
      IF (swId == -1) THEN
          msr = 'Failed to re-attach to swath '// TRIM(dgName)//' for writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write to fields

      start = 0
      stride = 1
      edge(1) = mm%nLevels

! Geolocation -- time & pressure

      status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mm%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mm%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(mm%pressure) )
      IF (status /= 0) THEN
          msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = swwrfld(swId, DG_FIELD2, start(1), stride(1), edge(1), &
           & mm%perMisPoints)
      IF (status /= 0) THEN
          msr = WR_ERR //  DG_FIELD2 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = lsize
      edge(2) = nsize
      status = swwrfld(swId, MIS_FIELD, start, stride, edge, mm%misDays)
      IF (status /= 0) THEN
          print *, 'swID ', swID
          print *, 'start ', start
          print *, 'stride ', stride
          print *, 'edge ', edge
          print *, 'missng. days ', mm%misDays
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
                & ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath ' // TRIM(dgName) // ' successfully written to file ' // &
            & trim(physicalFilename)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!------------------------------
   END SUBROUTINE OutputMMDiags_HE2
!------------------------------

!------------------------------------------------
   SUBROUTINE OutputMMDiags_HE5(physicalFilename, mm)
!------------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, MLS_CHARTYPE

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic portion of l3mm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename
      TYPE (L3MMData_T), INTENT(IN) :: mm

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DIMD_NAME   = 'lDate'
      CHARACTER (LEN=*), PARAMETER :: DIMMIS_NAME = 'nMisDays'
      CHARACTER (LEN=*), PARAMETER :: DIMDM_NAME  = 'lDate,nMisDays'
      CHARACTER (LEN=*), PARAMETER :: MIS_FIELD   = 'MissingDays'

! Functions

      INTEGER, EXTERNAL :: he5_swattach, he5_swclose, he5_swcreate, & 
           & he5_swdefdfld, he5_swdefdim, & 
           & he5_swdefgfld, he5_swdetach, he5_swopen, he5_swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      REAL, PARAMETER :: WAIT_TO_REOPEN = 2.0

      INTEGER :: edge(2), start(2), stride(2)
      INTEGER :: lsize, nsize, status, swfID, swId
      INTEGER, PARAMETER :: NUMBER_BEATS = 20000

! Open the output file

      swfID = he5_swopen(trim(physicalFilename), HE5F_ACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create the swath

      dgName = TRIM(mm%name) // 'Diagnostics'
      swId = he5_swcreate(swfID, trim(dgName))
      IF (swId == -1) THEN
          msr = 'Failed to create swath ' // trim(dgName)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the dimensions

      status = he5_swdefdim(swID, DIMR_NAME, MIN_MAX)
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

      !status = he5_swdefdim(swId, DIMD_NAME, lsize)
      status = he5_swdefdim(swId, DIMD_NAME, nsize)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMD_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefdim(swId, DIMMIS_NAME, nsize)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIMMIS_NAME
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefdim(swId, DIM_NAME2, mm%nLevels)
      IF (status /= 0) THEN
          msr = DIM_ERR // DIM_NAME2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "geolocation" fields

      status = he5_swdefgfld(swID, GEO_FIELD3, DIMR_NAME, "", & 
           & HE5T_NATIVE_DOUBLE, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefgfld(swId, GEO_FIELD9, DIM_NAME2, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the "data" fields

      status = he5_swdefdfld(swId, DG_FIELD2, DIM_NAME2, "", & 
           & HE5T_NATIVE_INT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // DG_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      !status = he5_swdefdfld(swId, MIS_FIELD, DIMD_NAME, "", & 
      !     & MLS_CHARTYPE, HDFE_NOMERGE)
      status = he5_swdefdfld(swId, MIS_FIELD, DIMD_NAME, "", & 
           & HE5T_NATIVE_INT, HDFE_NOMERGE)
      IF (status /= 0) THEN
          msr = DAT_ERR // MIS_FIELD
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write to fields

      start = 0
      stride = 1
      edge(1) = mm%nLevels

! Geolocation -- time & pressure

      status = he5_swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mm%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mm%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(mm%pressure) )
      IF (status /= 0) THEN
          msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = he5_swwrfld(swId, DG_FIELD2, start(1), stride(1), edge(1), &
           & mm%perMisPoints)
      IF (status /= 0) THEN
          msr = WR_ERR //  DG_FIELD2 // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = nsize
 
      status = he5_swwrfld(swId, MIS_FIELD, start(1), stride(1), edge(1), mm%misDays)
      IF (status /= 0) THEN
          print *, 'swID ', swID
          print *, 'start ', start
          print *, 'stride ', stride
          print *, 'edge ', edge
          print *, 'missng. days ', mm%misDays
          msr = WR_ERR //  MIS_FIELD // ' to swath ' // dgName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath after writing

      status = he5_swdetach(swId)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Close the file after writing

      status = he5_swclose(swfID)
      IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                & ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath ' // TRIM(dgName) // ' successfully written to file ' // &
            & trim(physicalFilename)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!------------------------------
   END SUBROUTINE OutputMMDiags_HE5
!------------------------------

!---------------------------------------------------------------
   SUBROUTINE WriteMetaL3MM (file, mlspcf_mcf_l3mm, pcf, anText, nFiles, &
	& hdfVersion)
!---------------------------------------------------------------
   USE HDFEOS5, ONLY: HE5F_ACC_RDWR, HE5_GDCLOSE, HE5_GDOPEN
   USE MLSPCF3
   USE mon_Open, ONLY: PCFMData_T
   USE output_m, ONLY: output
   USE PCFHdr, ONLY: WritePCF2Hdr, WriteInputPointer, he5_writeglobalattr, &
     & GlobalAttributes
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate 
   USE SDPToolkit, ONLY: PGS_S_SUCCESS, &
     & WARNIFCANTPGSMETREMOVE, PGSD_MET_GROUP_NAME_L, &
     & PGSD_MET_NUM_OF_GROUPS, PGSMET_E_MAND_NOT_SET, &
     & max_orbits

! Brief description of subroutine
! This routine writes the metadata for the l3mm file, and annotates it with the
! PCF.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER (LEN=*), INTENT(IN) :: file

      CHARACTER (LEN=1), POINTER :: anText(:)

      INTEGER, INTENT(IN) :: mlspcf_mcf_l3mm

      INTEGER, INTENT(IN) :: hdfVersion

      INTEGER, INTENT(IN) :: nFiles

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d, &
           & pgs_met_setAttr_i, pgs_met_setAttr_s, pgs_met_write, &
           & gdinqgrid, he5_gdinqgrid, pgs_pc_getUniversalRef

! Variables

      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=98) :: inpt(31)
      CHARACTER (LEN=GridNameLen) :: gridName
      CHARACTER (LEN=150) :: sval
      CHARACTER (LEN=45) :: attrName
      CHARACTER (LEN=3) :: cNum
      ! CHARACTER (LEN=2) :: fileType
      integer :: fileType

      REAL(r8) :: dval

      INTEGER :: dg, hdfReturn, i, indx, j, len, numGrids, result, & 
           & returnStatus, sdid, version

       INTEGER :: daysNum = 1 

! get the l3 start and end days

      daysNum = nFiles + 1 
      !call output('daysNum = ', advance='no')
      !call output(daysNum, advance='yes')

! Initialize the fileType
      fileType = l_grid

! Check to see whether this is a Diagnostic product

      dg = 0

      IF ( INDEX(file,'Diagnostic') /= 0 ) dg = 1

! Initialize the MCF file

      result = pgs_met_init(mlspcf_mcf_l3mm, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Initialization error.  See LogStatus for details.')

! Open the HDF file and initialize the SD interface

      sdid = mls_sfstart(trim(file), DFACC_RDWR, hdfVersion=hdfVersion, & 
           &  addingMetaData=.TRUE.)
      IF (sdid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
           & 'Failed to open the HDF file for metadata writing.')

! Set PGE values -- ECSDataGranule

      attrName = 'ReprocessingPlanned'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'further update anticipated using enhanced PGE')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'ReprocessingActual'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'unknown')
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

      numGrids = 0

      IF (hdfVersion==HDFVERSION_4) THEN 
         numGrids = gdinqgrid(trim(file), list, len)
      ELSE IF (hdfVersion==HDFVERSION_5) THEN
         numGrids = he5_gdinqgrid(trim(file), list, len)
      ENDIF

      IF (numGrids .LE. 0) THEN
         msr = 'No grids found in file ' // TRIM(file) // &
              & ' while attempting to write its metadata.'
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

         attrName = 'ParameterName' // '.' // trim(cNum)
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
              & gridName)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! QAStats Group

         attrName = 'QAPercentInterpolatedData' // '.' // trim(cNum)
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'QAPercentMissingData' // '.' // trim(cNum)
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'QAPercentOutofBoundsData' // '.' // trim(cNum)
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 0)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

! OrbitCalculatedSpatialDomainContainer

      ! This changes confirm James Johnson suggestion on 6/12/03
      ! Use 99999 for invalid value for now

      !attrName = 'OrbitNumber' // '.1'
      !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
      !IF (result /= PGS_S_SUCCESS) THEN
      !   msr = METAWR_ERR // attrName
      !   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      !ENDIF

      attrName = 'StartOrbitNumber' // '.1'
      !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
      !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
      IF (GlobalAttributes%OrbNumDays(1,daysNum) == -1) THEN
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, 99999)
      ELSE
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, GlobalAttributes%OrbNumDays(1,daysNum))
      ENDIF

      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'StopOrbitNumber' // '.1'
      !result = pgs_met_setAtttr_i(groups(INVENTORYMETADATA), attrName, -1)
      !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
      IF (maxval(GlobalAttributes%OrbNumDays(:,daysNum)) == -1) then
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, 99999)
      ELSE 
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, maxval(GlobalAttributes%OrbNumDays(:,daysNum)))
      END IF 

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
           & '00:00:00')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'EquatorCrossingDate' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
              & pcf%startDay )
              ! & '1899-04-29')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! InputPointer

      attrName = 'InputPointer'
! >       inpt = ''
! >       IF (dg == 1) THEN
! >          j = mlspcf_l2dg_start
! >       ELSE
! >          j = mlspcf_l2gp_start
! >       ENDIF
! >       DO i = 1, 31
! >          version = 1
! >          returnStatus = pgs_pc_getUniversalRef(j, version, sval)
! >          IF (returnStatus == PGS_S_SUCCESS) THEN
! >             inpt(i) = sval
! >             j = j + 1
! >          ENDIF
! >       ENDDO

      ! result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, inpt)
      result = WriteInputPointer(groups(INVENTORYMETADATA), attrName, &
        & fileType=l_hdfeos)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Locality Value

      attrName = 'LocalityValue'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'Limb')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! VerticalSpatialDomain Product-Specific Attribute

      attrName = 'VerticalSpatialDomainType' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'Atmosphere Layer')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'VerticalSpatialDomainValue' // '.1'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'Atmosphere Profile')
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! HorizontalSpatialDomainContainer

      attrName = 'ZoneIdentifier'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & 'Other Grid System')
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
           & pcf%startDay )
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
           & pcf%endDay )
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
           & pcf%outputVersion)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the metadata and their values to HDF attributes

      result = pgs_met_write(groups(INVENTORYMETADATA), "coremetadata", &
           & sdid)
      IF (result /= PGS_S_SUCCESS) THEN
         IF (result == PGSMET_E_MAND_NOT_SET) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                 & 'Some of the mandatory metadata parameters were not set.')
         ELSE
           CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                & 'Metadata write failed.')
         ENDIF
      ENDIF

! Terminate access to the SD interface and close the file

      hdfReturn = mls_sfend(sdid, hdfVersion=hdfVersion, & 
           & addingMetaData=.TRUE.)
      IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Error closing HDF file after writing metadata.')

      result = pgs_met_remove()
!      if (result /= PGS_S_SUCCESS .and. WARNIFCANTPGSMETREMOVE) THEN 
!        write(msr, *) result
!        CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!             & "Calling pgs_met_remove() failed with value " // trim(msr) )
!      endif          

! Write global attributes
   
      if ( HDFVersion == HDFVERSION_5 ) then                        
        sdid = he5_gdopen (file, HE5F_ACC_RDWR)                       
        call he5_writeglobalattr(sdid,daysNum) 
        result = he5_gdclose (sdid)                                   
      endif

      ! Annotate the file with the PCF
      CALL WritePCF2Hdr(file, anText, hdfVersion=hdfVersion, fileType=fileType)


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
      l3mm%nMisDays = 0

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
!# Revision 1.18  2004/12/13 17:41:26  cvuu
!# remove writing QA flags to meta file, use the ones in MCF v1.5
!#
!# Revision 1.17  2004/05/12 21:49:58  pwagner
!# Uses mls_h5open/close
!#
!# Revision 1.16  2004/05/04 15:57:20  cvuu
!# Fixed bug
!#
!# Revision 1.15  2004/01/08 21:21:36  cvuu
!# version 1.4 commit
!#
!# Revision 1.14  2003/09/16 16:35:11  cvuu
!# Add new parameter nFiles to the subroutine WriteMetaL3MM
!#
!# Revision 1.13  2003/09/15 18:27:23  cvuu
!# Output OrbitNumber and OrbitPeriod in the global attribute
!#
!# Revision 1.12  2003/08/11 23:26:37  cvuu
!# brought closer to James Johnson want to
!#
!# Revision 1.11  2003/07/08 00:17:46  pwagner
!# fileType now a lit_name instead of a char string
!#
!# Revision 1.10  2003/06/03 20:46:15  pwagner
!# Writes global attributes
!#
!# Revision 1.9  2003/06/02 23:45:15  pwagner
!# metadata chnages: OrbitNumber now -1; equatorCrossingDate now utc start date
!#
!# Revision 1.8  2003/05/30 23:54:07  pwagner
!# Relies on lib/PCFHdr to WriteInputPointer
!#
!# Revision 1.7  2003/04/30 18:16:28  pwagner
!# Work-around for LF95 infinite compile-time bug
!#
!# Revision 1.6  2003/04/06 02:25:41  jdone
!# added HDFEOS5 capability
!#
!# Revision 1.5  2003/03/15 00:20:02  pwagner
!# May warn if pgs_met_remove returns non-zero value
!#
!# Revision 1.4  2002/08/22 23:32:01  pwagner
!# Made start,endtimes pointer; initialize nMisDays
!#
!# Revision 1.3  2001/12/12 17:45:29  nakamura
!# Added dg fields; removed unused subroutine DestroyL3MMDatabase.
!#
!# Revision 1.2  2001/11/12 20:24:56  nakamura
!# Added L3MMDiag_T.
!#
!# Revision 1.1  2001/07/18 15:41:57  nakamura
!# Module for the L3MM data type.
!#
!#
