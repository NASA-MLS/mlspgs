! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE L3DMData
!==============================================================================

  USE L3CF, ONLY: L3CFDef_T, L3CFProd_T
  USE MLSCommon, ONLY: r8, FileNameLen
  USE MLSL3Common, ONLY: OutputFiles_T, DIM_ERR, GEO_ERR, DAT_ERR, WR_ERR, &
       & GD_ERR, METAWR_ERR, TAI2A_ERR, &
       & GEO_FIELD1, GEO_FIELD2, GEO_FIELD3, GEO_FIELD4, GEO_FIELD5, &
       & GEO_FIELD6, GEO_FIELD7, GEO_FIELD8, GEO_FIELD9, GEO_FIELD10, & 
       & DIM_NAME1, DIM_NAME2, DIM_NAME3, DIM_NAME12, DIM_NAME123, &
       & DIMX_NAME, DIMY_NAME, DIMZ_NAME, DIMT_NAME, &
       & DIML_NAME, DIMN_NAME, DIMR_NAME, DIMLL_NAME, DIMRL_NAME, &
       & MDT_FIELD, MD_FIELD, DG_FIELD, DG_FIELD1, DG_FIELD2, &
       & GCTP_GEO, CCSDS_LEN, GridNameLen, DIMXYZ_NAME, &
       & HDFE_NOMERGE, INVENTORYMETADATA
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, & 
       & MLSMSG_DEALLOCATE, MLSMSG_FILEOPEN, MLSMSG_ALLOCATE, MLSMSG_WARNING
  USE MLSStrings, only: utc_to_yyyymmdd
  IMPLICIT NONE
  private
  PUBLIC :: L3DMData_T, ConvertDeg2DMS, OutputGrids, &
    & WriteMetaL3DM, AllocateL3DM, DeallocateL3DM, DestroyL3DMDatabase

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
!                WriteMetaL3DM
!                AllocateL3DM
!                DeallocateL3DM
!                DestroyL3DMDatabase

! Remarks:  This module contains the definition of the L3DMData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELDV = 'L3dmValue'
   CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELDP = 'L3dmPrecision'

! This data type is used to store the l3 daily map data.

   TYPE L3DMData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:,:), POINTER :: l3dmValue=>NULL()	  ! Field value
     REAL(r8), DIMENSION(:,:,:), POINTER :: l3dmPrecision=>NULL() ! Field precision
	! dimensioned as (nLevels, nLats, nLons)

     REAL(r8), DIMENSION(:,:), POINTER :: latRss=>NULL()
	! Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

     REAL(r8), DIMENSION(:,:), POINTER :: maxDiff=>NULL()
	! Maximum difference, dimensioned (N, nLevels)

     REAL(r8), DIMENSION(:,:), POINTER :: maxDiffTime=>NULL()
	! Time of maximum differences, dimensioned (N, nLevels)

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure=>NULL()	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude=>NULL()	! dimensioned (nLats)
     REAL(r8), DIMENSION(:), POINTER :: longitude=>NULL()	! dimensioned (nLons)

     REAL(r8) :: time	! Synoptic time

     ! Now the diagnostic fields

     REAL(r8), DIMENSION(:), POINTER :: gRss=>NULL()
	! Global Root-Sum_Square, dimensioned (nLevels)

     INTEGER, DIMENSION(:), POINTER :: perMisPoints=>NULL()
	! Missing points (percentage), dimensioned (nLevels)

     INTEGER :: nLevels		! Total number of surfaces
     INTEGER :: nLats		! Total number of latitudes
     INTEGER :: nLons		! Total number of longitudes
     INTEGER :: N		! number for "largest differences" diagnostics

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
  SUBROUTINE OutputGrids(type, l3dmData, l3dmFiles, hdfVersion, createFile)
  !---------------------------------------------------
  USE MLSFiles, ONLY: HDFVERSION_5, HDFVERSION_4

    ! Brief description of subroutine
    ! This subroutine creates and writes to the grid portion of the l3dm files.
    
    ! Arguments

    TYPE (L3DMData_T), INTENT(IN) :: l3dmData(:)
    
    TYPE (OutputFiles_T), INTENT(INOUT) :: l3dmFiles
    
    CHARACTER (LEN=*), INTENT(IN) :: type
    
    INTEGER, INTENT(IN) :: hdfVersion
   
    LOGICAL, INTENT(INOUT) :: createFile

    if (hdfVersion == HDFVERSION_5) then 
       call OutputGrids_HE5(type, l3dmData, l3dmFiles, createFile)
    else if (hdfVersion == HDFVERSION_4) then
       call OutputGrids_HE2(type, l3dmData, l3dmFiles)
    endif
    
!----------------------------
  END SUBROUTINE OutputGrids
!----------------------------

  !---------------------------------------------------
  SUBROUTINE OutputGrids_HE2(type, l3dmData, l3dmFiles)
  !---------------------------------------------------
  USE HDF, ONLY: DFACC_RDWR, DFNT_FLOAT64, DFNT_FLOAT32, DFNT_INT32, & 
       & DFACC_WRITE
  USE MLSPCF3, ONLY: mlspcf_l3dm_start, mlspcf_l3dm_end
  USE MLSStrings, ONLY: LinearSearchStringArray
  USE PCFModule, ONLY: ExpandFileTemplate, FindFileDay 

    ! Brief description of subroutine
    ! This subroutine creates and writes to the grid portion of l3dm files.

    ! Arguments
     
    CHARACTER (LEN=*), INTENT(IN) :: type
     
    TYPE (L3DMData_T), INTENT(IN) :: l3dmData(:)
     
    TYPE (OutputFiles_T), INTENT(INOUT) :: l3dmFiles
   
    ! Parameters

    ! Variables

    CHARACTER (LEN=480) :: msr
    CHARACTER (LEN=FileNameLen) :: physicalFilename
    CHARACTER (LEN=8) :: date

    REAL(r8) :: projparm(13)
    REAL(r8) :: uplft(2), lowrgt(2)
    REAL :: maxLat, maxLon, minLat, minLon

    INTEGER :: start(3), stride(3), edge(3)
    INTEGER :: gdfID, gdId, i, match, status

    ! Functions

    INTEGER, EXTERNAL :: gdattach, gdclose, gdcreate, gddefdim, gddeffld, &
         & gddefproj, gddetach, gdopen, gdwrfld
     
    ! For each day in the l3dm database,
     
    DO i = 1, SIZE(l3dmData)

       ! Find the output file for the level/species/day in the PCF
        
       CALL FindFileDay(type, l3dmData(i)%time, mlspcf_l3dm_start, &
            & mlspcf_l3dm_end, match, physicalFilename, date)
       IF (match == -1) THEN
          msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
               & date
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Check whether the name is distinct; if so, save it in l3dmFiles

       IF (LinearSearchStringArray(l3dmFiles%name,physicalFilename)== 0) THEN
          l3dmFiles%nFiles = l3dmFiles%nFiles+1
          l3dmFiles%name(l3dmFiles%nFiles) = physicalFilename
          l3dmFiles%date(l3dmFiles%nFiles) = date
       ENDIF
         
       ! Open the output file
         
       gdfID = gdopen(trim(physicalFilename), DFACC_RDWR)
       IF (gdfID == -1) THEN
          msr = MLSMSG_Fileopen//trim(physicalFilename)//' for writing grid.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ! Set up the grid.  
       ! The region is bounded by 180W to 176E longitude varying latitude. 
       ! Grid into 90 bins along the x-axis, by nLats bins along
       ! the y-axis (4x2 bins).  
       ! Upper Left & Lower Right corners in DDDMMMSSS.ss.
         
       projparm = 0.0
         
       ! Find boundaries of measured latitude
       ! and upper bound of measure longitude

       maxLat = MAXVAL( l3dmData(i)%latitude )
       minLat = MINVAL( l3dmData(i)%latitude )
       maxLon = MAXVAL( l3dmData(i)%longitude)
       minLon = MINVAL( l3dmData(i)%longitude)

       ! Convert to "packed degree format"
         
       CALL ConvertDeg2DMS(maxLat, uplft(2) )
       CALL ConvertDeg2DMS(minLat, lowrgt(2))
       CALL ConvertDeg2DMS(maxLon, lowrgt(1))
       CALL ConvertDeg2DMS(minLon, uplft(1) )
         
       ! Create the grid

       gdId = gdcreate(gdfID, trim(l3dmData(i)%name), l3dmData(i)%nLons, &
            & l3dmData(i)%nLats, uplft, lowrgt)
       IF (gdId == -1) THEN
          msr = 'Failed to create grid ' // trim(l3dmData(i)%name)
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
         
       ! Define a "Geographic projection," 
       ! using defaults in all unneeded fields

       status = gddefproj(gdId, GCTP_GEO, 0, 0, projparm)
       IF (status /= 0) THEN
          msr = 'Failed to define projection for grid '//trim(l3dmData(i)%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Define the "geolocation" fields
         
       status = gddeffld(gdId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       status = gddeffld(gdId, GEO_FIELD9, DIMZ_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       status = gddeffld(gdId, GEO_FIELD1, DIMY_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       status = gddeffld(gdId, GEO_FIELD2, DIMX_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = GEO_ERR // GEO_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Define the "data" fields
         
       status = gddeffld(gdId, DATA_FIELDV, DIMXYZ_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDV
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       status = gddeffld(gdId, DATA_FIELDP, DIMXYZ_NAME, DFNT_FLOAT32, &
            & HDFE_NOMERGE)
       IF (status /= 0) THEN
          msr = DAT_ERR // DATA_FIELDP
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Detach from and close the grid interface.  
       ! This step is necessary to store
       ! properly the grid information within the file and must be done 
       ! before writing or reading data to or from the grid.
         
       status = gddetach(gdId)
       IF (status /= 0) THEN
          msr = GD_ERR // TRIM( l3dmData(i)%name ) // ' after definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       status = gdclose(gdfID)
       IF (status /= 0) THEN
          msr = 'Failed to close file ' // TRIM(physicalFilename) // &
               & ' after definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Re-open the file for writing
         
       gdfID = gdopen(trim(physicalFilename), DFACC_RDWR)
       IF (gdfID == -1) THEN
          msr = MLSMSG_Fileopen // TRIM(physicalFilename) & 
               & // ' for writing grid data.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Re-attach to the grid for writing
       
       gdId = gdattach(gdfID, trim(l3dmData(i)%name))
       IF (gdId == -1) THEN
          msr = 'Failed to attach to grid ' // trim(l3dmData(i)%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       ! Write to fields
         
       start  = 0
       stride = 1
       edge   = 1
         
       status = gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), edge(1), &
            & l3dmData(i)%time)
         
       IF (status /=0) THEN
          msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' &
               & // trim(l3dmData(i)%name)
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
         
       edge(1) = l3dmData(i)%nLevels
       edge(2) = l3dmData(i)%nLats
       edge(3) = l3dmData(i)%nLons
         
       if (l3dmData(i)%nLevels.gt.0) then
            
          status = gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1),&
               & REAL(l3dmData(i)%pressure) )
          IF (status /= 0) THEN
             msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' &
                  & // trim(l3dmData(i)%name)
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
         
      endif
      
      if (l3dmData(i)%nLats.gt.0) then
         
         status = gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2), &
              & REAL(l3dmData(i)%latitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if (l3dmData(i)%nLons.gt.0) then
         status = gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3), &
              & REAL(l3dmData(i)%longitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if ((l3dmData(i)%nLons.gt.0).and.(l3dmData(i)%nLats.gt.0).and. & 
           & (l3dmData(i)%nLevels.gt.0)) then
         
         status = gdwrfld( gdId, DATA_FIELDV, start, stride, edge, &
              & REAL(l3dmData(i)%l3dmValue) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDV // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = gdwrfld( gdId, DATA_FIELDP, start, stride, edge, &
              & REAL(l3dmData(i)%l3dmPrecision) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDP // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
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
              & ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      msr = 'Grid ' // TRIM(l3dmData(i)%name) // & 
           & ' successfully written to file ' // trim(physicalFilename)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)
      
! Create & write to diagnostic swath
      
      CALL OutputDiags_HE2( trim(physicalFilename), l3dmData(i) )
      
   ENDDO
   
 !----------------------------
 END SUBROUTINE OutputGrids_HE2
 !----------------------------

 
 !---------------------------------------------------
 SUBROUTINE OutputGrids_HE5(type, l3dmData, l3dmFiles, createFile)
 !---------------------------------------------------
   
  USE HDF5, ONLY: HID_T
  USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC 
  USE MLSPCF3, ONLY: mlspcf_l3dm_start, mlspcf_l3dm_end
  USE MLSStrings, ONLY: LinearSearchStringArray
  USE PCFModule, ONLY: ExpandFileTemplate, FindFileDay 
   ! Brief description of subroutine
   ! This subroutine creates and writes to the grid portion of the l3dm files.

   ! Arguments
   
   CHARACTER (LEN=*), INTENT(IN) :: type
   
   TYPE (L3DMData_T), INTENT(IN) :: l3dmData(:)
   
   TYPE (OutputFiles_T), INTENT(INOUT) :: l3dmFiles
  
   LOGICAL, INTENT(INOUT) :: createFile 

   ! Parameters
   
   ! Variables
   
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: physicalFilename
   CHARACTER (LEN=8) :: date
   
   REAL(r8) :: projparm(13)
   REAL(r8) :: uplft(2), lowrgt(2)
   REAL :: maxLat, maxLon, minLat, minLon
   
   INTEGER :: start(3), stride(3), edge(3)
   INTEGER(HID_T) :: gdfID, gdId
   INTEGER :: i, match, status

   ! Functions

   INTEGER, EXTERNAL :: & 
        & he5_gdattach, he5_gdclose, he5_gdcreate, he5_gddefdim,& 
        & he5_gddeffld, he5_gddefproj, he5_gddetach, he5_gdopen, he5_gdwrfld
      
   ! For each day in the l3dm database,
   
   DO i = 1, SIZE(l3dmData)
      
      ! Find the output file for the level/species/day in the PCF
      
      CALL FindFileDay(type, l3dmData(i)%time, mlspcf_l3dm_start, &
           & mlspcf_l3dm_end, match, physicalFilename, date)
      IF (match == -1) THEN
         msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
              & date
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Check whether the name is distinct; if so, save it in l3dmFiles
     
      IF (LinearSearchStringArray(l3dmFiles%name,physicalFilename) == 0) THEN
         l3dmFiles%nFiles = l3dmFiles%nFiles+1
         l3dmFiles%name(l3dmFiles%nFiles) = physicalFilename
         l3dmFiles%date(l3dmFiles%nFiles) = date
      ENDIF
      
      ! Open the output file

      IF (createFile) THEN  
          gdfID = he5_gdopen(physicalFilename, HE5F_ACC_TRUNC)
      ELSE
          gdfID = he5_gdopen(physicalFilename, HE5F_ACC_RDWR)
      ENDIF

      IF (gdfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename) //' for writing grid.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Set up the grid.  The region is bounded by 180.0W to 176.0E longitude &
      ! varying latitude.  
      ! Grid into 90 bins along the x-axis, by nLats bins along the y-axis 
      ! (4x2 bins).  
      ! Upper Left & Lower Right corners in DDDMMMSSS.ss.
      
      projparm = 0.0
      
      ! Find boundaries of measured latitude, upper bound of measure longitude
      
      maxLat = MAXVAL( l3dmData(i)%latitude  )
      minLat = MINVAL( l3dmData(i)%latitude  )
      maxLon = MAXVAL( l3dmData(i)%longitude )
      minLon = MINVAL( l3dmData(i)%longitude )
      
      ! Convert to "packed degree format"
      
      CALL ConvertDeg2DMS(maxLat, uplft(2) )
      CALL ConvertDeg2DMS(minLat, lowrgt(2))
      CALL ConvertDeg2DMS(maxLon, lowrgt(1))
      CALL ConvertDeg2DMS(minLon, uplft(1) )
      
      ! Create the grid

      gdId = he5_gdcreate(gdfID, trim(l3dmData(i)%name), l3dmData(i)%nLons, &
           & l3dmData(i)%nLats, uplft, lowrgt)
      IF (gdId == -1) THEN
         msr = 'Failed to create grid ' // trim(l3dmData(i)%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the dimensions

      status = he5_gddefdim(gdId, DIMX_NAME, l3dmData(i)%nLons)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMX_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_gddefdim(gdId, DIMY_NAME, l3dmData(i)%nLats)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMY_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_gddefdim(gdId, DIMZ_NAME, l3dmData(i)%nLevels)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMZ_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_gddefdim(gdId, DIMT_NAME, 1)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define a "Geographic projection," using defaults in all unneeded fields
      
      status = he5_gddefproj(gdId, GCTP_GEO, 0, 0, projparm)
      IF (status /= 0) THEN
         msr = 'Failed to define projection for grid '// trim(l3dmData(i)%name)
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
      
      status = he5_gddeffld(gdId, DATA_FIELDV, DIMXYZ_NAME, "", &
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_gddeffld(gdId, DATA_FIELDP, DIMXYZ_NAME, "", &
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Detach from and close the grid interface.
      ! This step is necessary to store
      ! properly the grid information within the file and must be done
      ! before writing or reading data to or from the grid.
                                                                                           
      status = he5_gddetach(gdId)
      IF (status /= 0) THEN
         msr = GD_ERR // TRIM( l3dmData(i)%name ) // ' after definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Re-attach to the grid for writing
                                                                                           
      gdId = he5_gdattach(gdfID, l3dmData(i)%name)
      IF (gdId == -1) THEN
         msr = 'Failed to attach to grid ' // trim(l3dmData(i)%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Write to fields
      
      start = 0
      stride = 1
      edge = 1
      
      status = he5_gdwrfld(gdId, GEO_FIELD3, start(1), stride(1), edge(1), &
           & l3dmData(i)%time)

      IF (status /=0) THEN
         msr = 'Failed to write field ' //  GEO_FIELD3 // ' to grid ' &
              & // trim(l3dmData(i)%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      edge(1) = l3dmData(i)%nLevels
      edge(2) = l3dmData(i)%nLats
      edge(3) = l3dmData(i)%nLons
      
      if (l3dmData(i)%nLevels.gt.0) then
         
         status = he5_gdwrfld( gdId, GEO_FIELD9, start(1), stride(1), edge(1),&
              & REAL(l3dmData(i)%pressure) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD9 // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if (l3dmData(i)%nLats.gt.0) then
         
         status = he5_gdwrfld( gdId, GEO_FIELD1, start(2), stride(2), edge(2),&
              & REAL(l3dmData(i)%latitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD1 // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if (l3dmData(i)%nLons.gt.0) then
         status = he5_gdwrfld( gdId, GEO_FIELD2, start(3), stride(3), edge(3),&
              & REAL(l3dmData(i)%longitude) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  GEO_FIELD2 // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if ((l3dmData(i)%nLons.gt.0).and.(l3dmData(i)%nLats.gt.0).and. & 
           & (l3dmData(i)%nLevels.gt.0)) then
         
         status = he5_gdwrfld( gdId, DATA_FIELDV, start, stride, edge, &
              & REAL(l3dmData(i)%l3dmValue) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDV // ' to grid ' &
                 & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_gdwrfld( gdId, DATA_FIELDP, start, stride, edge, &
              & REAL(l3dmData(i)%l3dmPrecision) )
         IF (status /= 0) THEN
            msr = 'Failed to write field ' //  DATA_FIELDP // ' to grid ' &
                  & // trim(l3dmData(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         endif
         
! Detach from the grid after writing
         
         status = he5_gddetach(gdId)
         IF (status /= 0) THEN
            msr = GD_ERR // TRIM( l3dmData(i)%name ) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
! Close the file after writing
         
         status = he5_gdclose(gdfID)
         IF (status /= 0) THEN
            msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                 & ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         msr = 'Grid ' // TRIM(l3dmData(i)%name) // & 
              & ' successfully written to file ' // trim(physicalFilename)
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)
         
         ! Create & write to diagnostic swath
         
         CALL OutputDiags_HE5( trim(physicalFilename), l3dmData(i) )
         
      ENDDO
      
      createFile = .false. 
    !----------------------------
    END SUBROUTINE OutputGrids_HE5
    !----------------------------
    
   !----------------------------------------------
    SUBROUTINE OutputDiags_HE2(physicalFilename, dg)
    !----------------------------------------------
  USE HDF, ONLY: DFACC_RDWR, DFNT_FLOAT64, DFNT_FLOAT32, DFNT_INT32, & 
       & DFACC_WRITE

      ! Brief description of subroutine
      ! This subroutine creates and writes to 
      ! the diagnostic portion of the l3dm files.

      ! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename
      
      TYPE (L3DMData_T), INTENT(IN) :: dg
      
      ! Parameters
      
      CHARACTER (LEN=*), PARAMETER :: DIMNL_NAME = 'N,nLevels'
            
      ! Variables
      
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName
      
      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: swfID, swId, status, n, m

      ! Functions
      
      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, &
           & swdefgfld, swdetach, swopen, swwrfld
      
      ! Re-open the file for the creation of diagnostic swaths
      
      swfID = swopen(trim(physicalFilename), DFACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Create the swath
      
      dgName = TRIM(dg%name) // 'Diagnostics'
      
      swId = swcreate(swfID, dgName)
      IF (swId == -1) THEN
         msr = 'Failed to create swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the dimensions
      
      status = swdefdim(swId, DIMN_NAME, dg%N)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMN_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdim(swId, DIML_NAME, dg%nLats)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdim(swId, DIM_NAME2, dg%nLevels)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdim(swId, DIMT_NAME, 1)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the "geolocation" fields using the above dimensions
      
      status = swdefgfld(swId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefgfld(swId, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefgfld(swId, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the "data" fields
      
      status = swdefdfld(swId, DG_FIELD, DIM_NAME2, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdfld(swId, DG_FIELD1, DIMLL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdfld(swId, MD_FIELD, DIMNL_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // MD_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdfld(swId, MDT_FIELD, DIMNL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // MDT_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdfld(swId, DG_FIELD2, DIM_NAME2, DFNT_INT32, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Detach from the swath interface after definition
      
      status = swdetach(swId)
      IF (status /= 0) THEN
         msr = 'Failed to detach from swath ' // TRIM(dgName) // &
              & ' after definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Re-attach to the swath for writing
      
      swId = swattach(swfID, dgName)
      IF (swId == -1) THEN
         msr = 'Failed to re-attach to swath ' // TRIM(dgName) & 
              & //' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Write to fields

      start  = 0
      stride = 1
      edge   = 1
      
      ! Geolocation

      status = swwrfld(swId, GEO_FIELD3, start(1), stride(1), edge(1), dg%time)
      IF (status /=0) THEN
         msr = WR_ERR //  GEO_FIELD3 // ' to swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dg%nLevels
      edge(2) = dg%nLats

      if (dg%nLevels.gt.0) then
         
         status = swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1), &
              & REAL(dg%pressure) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if (dg%nLats.gt.0) then
         
         status = swwrfld( swId, GEO_FIELD1, start(2), stride(2), edge(2), &
              & REAL(dg%latitude) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD1 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
! One-dimensional data fields
      
      if (dg%nLevels.gt.0) then
         
         status = swwrfld( swId, DG_FIELD, start(1), stride(1), edge(1), &
              & REAL(dg%gRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = swwrfld(swId, DG_FIELD2, start(1), stride(1), edge(1), &
              & dg%perMisPoints )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD2 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      endif
      
      ! Two-dimensional data fields

      if ((dg%nLevels.gt.0).and.(dg%nLats.gt.0)) then
         
         status = swwrfld( swId, DG_FIELD1, start, stride, edge, & 
              & real(dg%latRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD1 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      edge(1) = dg%N
      edge(2) = dg%nLevels
      
      if ( (dg%N.gt.0).and.(dg%nLevels.gt.0) )  then
         
         status = swwrfld( swId, MD_FIELD, start, stride, edge, & 
              & real(dg%maxDiff) )
         IF (status /= 0) THEN
            msr = WR_ERR //  MD_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swId, MDT_FIELD, start, stride, edge, & 
              & real(dg%maxDiffTime) )
         IF (status /= 0) THEN
            msr = WR_ERR //  MDT_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      ! Detach from the swath after writing
      
      status = swdetach(swId)
      IF (status /= 0) THEN
         msr = GD_ERR // TRIM(dgName) // ' after writing.'
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
      
    !----------------------------
    END SUBROUTINE OutputDiags_HE2
    !----------------------------
    
    !----------------------------------------------
    SUBROUTINE OutputDiags_HE5(physicalFilename, dg)
    !----------------------------------------------
  USE HDF5, ONLY: HID_T
  USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC 

      ! Brief description of subroutine
      ! This subroutine creates and writes to 
      ! the diagnostic portion of the l3dm files.

      ! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: physicalFilename
      
      TYPE (L3DMData_T), INTENT(IN) :: dg

      ! Parameters
      
      CHARACTER (LEN=*), PARAMETER :: DIMNL_NAME = 'N,nLevels'
            
      ! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName
      
      INTEGER :: start(2), stride(2), edge(2)
      INTEGER (HID_T) :: swfID, swId
      INTEGER :: status, n, m
      
      ! Functions

      INTEGER, EXTERNAL :: he5_swattach, he5_swclose, he5_swcreate, & 
           & he5_swdefdfld, he5_swdefdim, &
           & he5_swdefgfld, he5_swdetach, he5_swopen, he5_swwrfld
 
      ! Re-open the file for the creation of diagnostic swaths
      
      swfID = he5_swopen(trim(physicalFilename), HE5F_ACC_RDWR)
      
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(physicalFilename) //' for writing swath'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Create the swath
      
      dgName = TRIM(dg%name) // 'Diagnostics'

      swId = he5_swcreate(swfID, dgName)
      IF (swId == -1) THEN
         msr = 'Failed to create swath ' // dgName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Define the dimensions

      status = he5_swdefdim(swId, DIMN_NAME, dg%N)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMN_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdim(swId, DIML_NAME, dg%nLats)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdim(swId, DIM_NAME2, dg%nLevels)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdim(swId, DIMT_NAME, 1)
      IF (status /= 0) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the "geolocation" fields using the above dimensions
      
      status = he5_swdefgfld(swId, GEO_FIELD3, DIMT_NAME, "", & 
           & HE5T_NATIVE_DOUBLE, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefgfld(swId, GEO_FIELD9, DIM_NAME2, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefgfld(swId, GEO_FIELD1, DIML_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Define the "data" fields

      status = he5_swdefdfld(swId, DG_FIELD, DIM_NAME2, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF      
      
      status = he5_swdefdfld(swId, DG_FIELD1, DIMLL_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdfld(swId, MD_FIELD, DIMNL_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // MD_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdfld(swId, MDT_FIELD, DIMNL_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // MDT_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefdfld(swId, DG_FIELD2, DIM_NAME2, "", & 
           & HE5T_NATIVE_INT, HDFE_NOMERGE)
      IF (status /= 0) THEN
         msr = DAT_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Write to fields

      start  = 0
      stride = 1
      edge   = 1
      
      ! Geolocation

      status = he5_swwrfld(swId, GEO_FIELD3, start(1), stride(1), edge(1), & 
           & dg%time)
      IF (status /=0) THEN
         msr = WR_ERR //  GEO_FIELD3 // ' to swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      edge(1) = dg%nLevels
      edge(2) = dg%nLats
      
      if (dg%nLevels.gt.0) then
         
         status = he5_swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1),&
              & REAL(dg%pressure) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      if (dg%nLats.gt.0) then
         
         status = he5_swwrfld( swId, GEO_FIELD1, start(2), stride(2), edge(2),&
              & REAL(dg%latitude) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD1 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif

      ! One-dimensional data fields

      if (dg%nLevels.gt.0) then

         status = he5_swwrfld( swId, DG_FIELD, start(1), stride(1), edge(1), &
              & REAL(dg%gRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_swwrfld(swId, DG_FIELD2, start(1), stride(1), edge(1), &
              & dg%perMisPoints )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD2 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      endif
      
      ! Two-dimensional data fields

      if ((dg%nLevels.gt.0).and.(dg%nLats.gt.0)) then
      
         status = he5_swwrfld( swId, DG_FIELD1, start, stride, edge, & 
              & real(dg%latRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD1 // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      edge(1) = dg%N
      edge(2) = dg%nLevels
      
      if ( (dg%N.gt.0).and.(dg%nLevels.gt.0) )  then
         
         status = he5_swwrfld( swId, MD_FIELD, start, stride, edge, & 
              & real(dg%maxDiff) )
         IF (status /= 0) THEN
            msr = WR_ERR //  MD_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_swwrfld( swId, MDT_FIELD, start, stride, edge, & 
              & dg%maxDiffTime )
              !& real(dg%maxDiffTime) )
         
         IF (status /= 0) THEN
            msr = WR_ERR //  MDT_FIELD // ' to swath ' // trim(dgName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
      endif
      
      ! Detach from the swath after writing
      
      status = he5_swdetach(swId)
      IF (status /= 0) THEN
         msr = GD_ERR // TRIM(dgName) // ' after writing.'
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
      
    !----------------------------
    END SUBROUTINE OutputDiags_HE5
    !----------------------------
   
    !-----------------------------------------------------
    SUBROUTINE WriteMetaL3DM (pcf, l3cf, files, anText, hdfVersion)
    !-----------------------------------------------------
  USE HDF, ONLY: DFACC_RDWR, DFNT_FLOAT64, DFNT_FLOAT32, DFNT_INT32, & 
       & DFACC_WRITE
  USE HDFEOS5, ONLY: HE5F_ACC_RDWR, HE5_GDCLOSE, HE5_GDOPEN
  USE Intrinsic, ONLY: l_grid, l_hdfeos
  USE MLSFiles, ONLY: HDFVERSION_5, HDFVERSION_4, mls_sfstart, mls_sfend
  USE OpenInit, ONLY: PCFData_T
  USE PCFHdr, ONLY: WritePCF2Hdr, WriteInputPointer, he5_writeglobalattr, &
     & GlobalAttributes
  USE PCFModule, ONLY: ExpandFileTemplate, FindFileDay 
  USE SDPToolkit, ONLY: PGSd_MET_NUM_OF_GROUPS, PGSd_MET_GROUP_NAME_L, & 
     & PGS_S_SUCCESS, PGSMET_E_MAND_NOT_SET, WARNIFCANTPGSMETREMOVE, &
     & max_orbits

      ! Brief description of subroutine
      ! This routine writes the metadata for an l3dm file, 
      ! and annotates it with the PCF.

      ! Arguments
      
      TYPE( L3CFProd_T ), INTENT(IN) :: l3cf
      
      TYPE (OutputFiles_T), INTENT(IN) :: files
      
      TYPE( PCFData_T ), INTENT(IN) :: pcf
      
      CHARACTER (LEN=1), POINTER :: anText(:)

      CHARACTER (LEN=8) :: rangeDate 

      INTEGER, INTENT(IN), OPTIONAL :: hdfVersion
      
      ! Parameters
      
      ! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d, & 
           & pgs_met_setAttr_i, pgs_met_setAttr_s, pgs_met_write, & 
           & pgs_td_taiToUTC, gdinqgrid, he5_gdinqgrid

      ! Variables

      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER (LEN=GridNameLen)           :: gridName
      CHARACTER (LEN=FileNameLen)           :: sval
      CHARACTER (LEN=CCSDS_LEN)             :: timeA
      CHARACTER (LEN=480)                   :: msr
      CHARACTER (LEN=132)                   :: list,gridlist
      CHARACTER (LEN=45)                    :: attrName
      ! CHARACTER (LEN=2)                     :: fileType
      integer ::                             fileType
      CHARACTER (LEN=1)                     :: cNum
      
      REAL(r8) :: dval

      INTEGER :: hdfReturn, i, j, indx, numGrids, result, sdid, &  
           & returnStatus, MyHDFVersion, len=0

      ! Initialize MyHDFVersion and use HDFEOS2 as default.

      IF (PRESENT(hdfVersion)) THEN 
         MyHDFVersion = hdfVersion
      ELSE 
         MyHDFVersion = HDFVERSION_4
      ENDIF

      ! Initialize the MCF file

      result = pgs_met_init(l3cf%mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Initialization error.  See LogStatus for details.')
      
      ! set the file type

      fileType = l_grid

      ! For each l3dm file successfully created,
      
      DO i = 1, files%nFiles
         
         ! Open the HDF file and initialize the SD interface

         sdid = mls_sfstart(trim(files%name(i)), DFACC_RDWR, & 
              & hdfVersion=MyHDFVersion, addingMetaData=.true.)
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
              & 'processed once')
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
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName,'Both')
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

         IF (MyHDFVersion==HDFVERSION_4) THEN 
            numGrids = gdinqgrid(files%name(i), list, len)
         ELSE IF (MyHDFVersion==HDFVERSION_5) THEN
            numGrids = he5_gdinqgrid(trim(files%name(i)), gridlist, len)
         ENDIF

         IF (numGrids .LE. 0) THEN
            msr = 'No grids found in file ' // TRIM( files%name(i) ) // &
                 & ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         ! For each grid in the file
        
         list = gridlist(:len) 
      
         DO j = 1, numGrids
            
            ! Extract its name
            
            indx = INDEX(list, ',')
            IF (indx /= 0) THEN
               gridName = list(:indx-1)
               list = list(indx+1:)
            ELSE
               !gridName = gridlist
               gridName = trim(list)
            ENDIF
           
            ! Append a class suffix to ParameterName, 
            ! and write the grid name as its value
            
            WRITE( cNum, '(I1)' ) j
            
            attrName = 'ParameterName' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & gridName)
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
            
            ! QAFlags Group
            
            attrName = 'AutomaticQualityFlag' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Passed')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
            
            attrName = 'AutomaticQualityFlagExplanation' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'pending algorithm update')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
            
            attrName = 'OperationalQualityFlag' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Not Investigated')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'OperationalQualityFlagExplanation' // '.' // cNum
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Not Investigated')
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
        
	 ! This changes confirm James Johnson suggestion on 6/12/03 
	 ! use 99999 for invalid value for now

         !attrName = 'OrbitNumber' // '.1'
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
         !IF (result /= PGS_S_SUCCESS) THEN
         !  msr = METAWR_ERR // attrName
         !   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         !ENDIF
        
         attrName = 'StartOrbitNumber' // '.1'
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
         IF (GlobalAttributes%OrbNumDays(1,i) == -1) THEN
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
		& attrName, 99999)
	 ELSE 
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
		& attrName, GlobalAttributes%OrbNumDays(1,i))
         ENDIF 

         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         attrName = 'StopOrbitNumber' // '.1'
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
         if (maxval(GlobalAttributes%OrbNumDays(:,i)) == -1) then
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
		& attrName, 99999)
	 else
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
		& attrName, maxval(GlobalAttributes%OrbNumDays(:,i)))
         end if
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
              & files%date(i) )
              ! & '1899-04-29')
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! InputPointer
         
         attrName = 'InputPointer'
         ! result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
         !     & 'See the PCF annotation to this file.')
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
              & files%date(i) )
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         attrName = 'RangeBeginningTime'
         !returnStatus = pgs_td_taiToUTC(l3cf%timeD(1), timeA)
         !IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
         !     & ModuleName, TAI2A_ERR)
         !sval = timeA(12:26)
         sval = '00:00:00.000000'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         attrName = 'RangeEndingDate'
         result = pgs_met_setAttr_s( groups(INVENTORYMETADATA), attrName, &
              & files%date(i) )
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
                    &'Some of the mandatory metadata parameters were not set.')
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                    &'Metadata write failed.')
            ENDIF
         ENDIF
         
         ! Terminate access to the SD interface and close the file

         hdfReturn = mls_sfend(sdid, hdfVersion=MyHDFVersion, & 
              & addingMetaData=.true.)
         IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
              & 'Error closing HDF file after writing metadata.')
         
         ! Annotate the file with the PCF
  
         CALL WritePCF2Hdr(files%name(i), anText, & 
              & hdfVersion=MyHDFVersion, fileType=fileType)
        
	! Set global attributes of Granule for each day of L3 Daily output file 
        rangeDate = files%date(i)
        GlobalAttributes%StartUTC = rangeDate(1:4)//'-'//rangeDate(6:8) &
	  & // 'T00:00:00.000000Z'
        GlobalAttributes%EndUTC = rangeDate(1:4)//'-'//rangeDate(6:8) &
	  & // 'T23:59:59.999999Z'
        call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
          & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
          & GlobalAttributes%GranuleDay)
          
        ! Write global attributes
        if ( MyHDFVersion == HDFVERSION_5 ) then
          sdid = he5_gdopen (files%name(i), HE5F_ACC_RDWR)
          call he5_writeglobalattr(sdid,i)
          result = he5_gdclose (sdid)
        endif
      ENDDO

      result = pgs_met_remove()
!      if (result /= PGS_S_SUCCESS .and. WARNIFCANTPGSMETREMOVE) THEN 
!         write(msr, *) result
!         CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!              & "Calling pgs_met_remove() failed with value " // trim(msr) )
!      endif
      

!------------------------------
   END SUBROUTINE WriteMetaL3DM
!------------------------------

!-----------------------------------------------------
   SUBROUTINE AllocateL3DM (nlev, nlat, nlon, N, l3dm)
!-----------------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3DMData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat, nlon, N

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
      l3dm%N = N

! Horizontal geolocation fields

      if (l3dm%nLats .gt. 0) then

         ALLOCATE(l3dm%latitude(l3dm%nLats), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' latitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

      if (l3dm%nLons .gt. 0) then

         ALLOCATE(l3dm%longitude(l3dm%nLons), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' longitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

! Vertical geolocation field

      if (l3dm%nLevels .gt. 0) then

         ALLOCATE(l3dm%pressure(l3dm%nLevels), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

! Data fields

      if ((l3dm%nLevels.gt.0).and.(l3dm%nLats.gt.0).and. & 
           & (l3dm%nLons.gt.0)) then 

         ALLOCATE(l3dm%l3dmValue(l3dm%nLevels,l3dm%nLats,l3dm%nLons), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' value pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ALLOCATE(l3dm%l3dmPrecision(l3dm%nLevels,l3dm%nLats,l3dm%nLons), & 
              & STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' precision pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

! Diagnostic fields

      if (l3dm%nLevels.gt.0) then 

         ALLOCATE(l3dm%gRss(l3dm%nLevels), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' gRss pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

      if ((l3dm%nLevels.gt.0).and.(l3dm%nLats.gt.0)) then
 
         ALLOCATE(l3dm%latRss(l3dm%nLevels,l3dm%nLats),STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' latRss pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

      if ((l3dm%nLevels.gt.0).and.(l3dm%N.gt.0)) then

         ALLOCATE(l3dm%maxDiff(l3dm%N,l3dm%nLevels),STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' maxDiff pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ALLOCATE(l3dm%maxDiffTime(l3dm%N,l3dm%nLevels),STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' maxDiffTime pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

      if (l3dm%nLevels .gt. 0) then 
         
         ALLOCATE(l3dm%perMisPoints(l3dm%nLevels), STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' perMisPoints pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

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

      if (l3dm%nLats .gt. 0) then

         IF ( ASSOCIATED(l3dm%latitude) ) THEN
            DEALLOCATE (l3dm%latitude, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  latitude pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      if (l3dm%nLons .gt. 0) then

         IF ( ASSOCIATED(l3dm%longitude) ) THEN
            DEALLOCATE (l3dm%longitude, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  longitude pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      ! Vertical geolocation field

      if (l3dm%nLevels .gt. 0) then

         IF ( ASSOCIATED(l3dm%pressure) ) THEN
            DEALLOCATE (l3dm%pressure, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  pressure pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      ! Data fields
      
      if ((l3dm%nLevels.gt.0).and.(l3dm%nLats.gt.0).and.& 
           & (l3dm%nLons.gt.0)) then 

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

      endif

      ! Diagnostic fields
      
      if (l3dm%nLevels.gt.0) then 

         IF ( ASSOCIATED(l3dm%gRss) ) THEN
            DEALLOCATE (l3dm%gRss, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  gRss pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      if ((l3dm%nLevels.gt.0).and.(l3dm%nLats.gt.0)) then
         
         IF ( ASSOCIATED(l3dm%latRss) ) THEN
            DEALLOCATE (l3dm%latRss, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  latRss pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      if ((l3dm%nLevels.gt.0).and.(l3dm%N.gt.0)) then

         IF ( ASSOCIATED(l3dm%maxDiff) ) THEN
            DEALLOCATE (l3dm%maxDiff, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  maxDiff pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

         IF ( ASSOCIATED(l3dm%maxDiffTime) ) THEN
            DEALLOCATE (l3dm%maxDiffTime, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  maxDiffTime pointer'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

      if (l3dm%nLevels .gt. 0) then 

         IF ( ASSOCIATED(l3dm%perMisPoints) ) THEN
            DEALLOCATE (l3dm%perMisPoints, STAT=err)
            IF ( err /= 0 ) THEN
               msr = MLSMSG_DeAllocate // '  perMisPoints'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDIF

      endif

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
!# Revision 1.35  2004/03/19 14:25:26  cvuu
!# Fix the RangeBeginningTime in metadata file
!#
!# Revision 1.34  2004/01/08 21:25:44  cvuu
!# Correct the day, month and doy in the global attribute
!#
!# Revision 1.33  2004/01/07 21:43:18  cvuu
!# version 1.4 commit
!#
!# Revision 1.32  2003/09/15 18:29:41  cvuu
!# Add OrbitNumber and OrbitPeriod to global attribute
!#
!# Revision 1.31  2003/08/21 20:22:46  cvuu
!# Remove testing print statements
!#
!# Revision 1.30  2003/08/18 17:03:59  cvuu
!# Add close and open in outputDiags_he5
!#
!# Revision 1.29  2003/08/11 23:29:09  cvuu
!# brought closer to James Johnson want to
!#
!# Revision 1.28  2003/07/08 00:18:09  pwagner
!# fileType now a lit_name instead of a char string
!#
!# Revision 1.27  2003/06/03 20:45:29  pwagner
!# Writes global attributes
!#
!# Revision 1.26  2003/06/02 23:45:15  pwagner
!# metadata chnages: OrbitNumber now -1; equatorCrossingDate now utc start date
!#
!# Revision 1.25  2003/05/30 23:52:42  pwagner
!# Relies on lib/PCFHdr to WriteInputPointer
!#
!# Revision 1.24  2003/04/30 18:15:48  pwagner
!# Work-around for LF95 infinite compile-time bug
!#
!# Revision 1.23  2003/04/08 20:26:46  pwagner
!# Snipped intraline continuation character to appease Lahey
!#
!# Revision 1.22  2003/03/22 02:13:20  jdone
!# added HDFEOS5/HDFEOS2 capability
!#
!# Revision 1.21  2003/03/15 00:16:38  pwagner
!# May warn if pgs_met_remove returns non-zero value
!#
!# Revision 1.20  2002/04/10 22:00:34  jdone
!# swwrfld edge values for time
!#
!# Revision 1.19  2002/04/01 21:58:49  jdone
!# check if array sizes are larger than 0
!#
!# Revision 1.18  2002/03/27 21:33:32  jdone
!# allocate statements checked and maxDiff is initialized
!#
!# Revision 1.17  2001/12/13 20:47:50  nakamura
!# Removed unused subroutine ReadL3DMData; merged dg fields into L3DMData_T.
!#
!# Revision 1.16  2001/11/26 19:24:20  nakamura
!# Moved L3DMDiag_T to its own module.
!#
!# Revision 1.15  2001/11/12 20:22:05  nakamura
!# Added pressure & lat to L3DMDiag_T.
!#
!# Revision 1.14  2001/10/05 20:16:25  nakamura
!# Added L3DMDiag_T.
!#
!# Revision 1.13  2001/07/18 15:54:37  nakamura
!# Gets metadata time from l3cf, rather than hard-coded to noon.
!#
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
