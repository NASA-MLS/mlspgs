! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE L3MZData
!==============================================================================

   USE HDF, ONLY: DFACC_RDWR, DFACC_WRITE, DFNT_FLOAT32, DFNT_INT32, & 
        & DFNT_FLOAT64
   USE Intrinsic, ONLY: l_hdfeos, l_swath
   USE MLSCommon, ONLY: r8
   USE MLSFiles, ONLY: MLS_SFSTART, MLS_SFEND, mls_inqswath, & 
        & HDFVERSION_4, HDFVERSION_5
   USE MLSL3Common, ONLY: GridNameLen, DIMR_NAME, DIMT_NAME, DIM_NAME2, &
        & DIML_NAME, SZ_ERR, DAT_ERR, GEO_ERR, SW_ERR, DG_FIELD1, DG_FIELD2, &
        & DIMLL_NAME, DATE_LEN, INVENTORYMETADATA, OutputFiles_T, &
        & GEO_FIELD1, GEO_FIELD3, GEO_FIELD4, GEO_FIELD9, GEO_FIELD11, & 
        & GEO_FIELD12, &
        & HDFE_NOMERGE, DIM_ERR, MIN_MAX, DIMRL_NAME, WR_ERR, METAWR_ERR, &
        & FILENAMELEN
   USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
        & MLSMSG_Info, MLSMSG_WARNING, MLSMSG_Allocate, MLSMSG_DeAllocate
   USE MLSPCF3, ONLY: MLSPCF_L2DG_START, MLSPCF_L2GP_START
   USE mon_Open, ONLY: PCFMData_T
   USE OUTPUT_M, ONLY: output
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate
   USE PCFHdr, only: GlobalAttributes
   USE SDPToolkit, only: WARNIFCANTPGSMETREMOVE, PGSD_MET_NUM_OF_GROUPS, &
        & PGSD_MET_GROUP_NAME_L, PGSMET_E_MAND_NOT_SET, PGS_S_SUCCESS, &
	& max_orbits
   ! USE SWAPI
   ! USE HE5_SWAPI

   IMPLICIT NONE
   private
   PUBLIC :: L3MZData_T, &
     & OutputL3MZ, WriteMetaL3MZ, AllocateL3MZ, DeallocateL3MZ

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L3MZData_T
! Subroutines -- OutputL3MZ
!                OutputMZDiag
!                WriteMetaL3MZ
!                AllocateL3MZ
!                DeallocateL3MZ

! Remarks:  This module contains the definition of the L3MZData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3mzValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3mzPrecision'

! This data type is used to store the l3 monthly zonal mean data.

  TYPE L3MZData_T

     CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

     REAL(r8), DIMENSION(:,:), POINTER :: localSolarTime
     REAL(r8), DIMENSION(:,:), POINTER :: localSolarZenithAngle
	! dimensioned as (2, nLats)

     ! Now the data fields:

     REAL(r8), DIMENSION(:,:), POINTER :: l3mzValue	! Field value
     REAL(r8), DIMENSION(:,:), POINTER :: l3mzPrecision	! Field precision
        ! dimensioned as (nLevels, nLats)

     ! Now the diagnostic fields:

     REAL(r8), DIMENSION(:,:), POINTER :: latRss
	! Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)

     REAL(r8) :: startTime	! start time of L2 data used in the analysis
     REAL(r8) :: endTime 	! end time of L2 data used in the analysis

     INTEGER, DIMENSION(:,:), POINTER :: perMisPoints
	! Missing points (percentage), dimensioned (nLevels, nLats)

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes

   END TYPE L3MZData_T

CONTAINS

!----------------------------------------------
   SUBROUTINE OutputL3MZ (file, mz, createFlag, hdfVersion, prodCount)
!----------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3mz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      TYPE( L3MZData_T ), INTENT(IN) :: mz

      LOGICAL, INTENT(INOUT) :: createFlag

      INTEGER, INTENT(IN) :: hdfVersion, prodCount

      IF (hdfVersion == HDFVERSION_4) THEN 
         CALL OutputL3MZ_HE2 (file, mz, createFlag)
      ELSE IF (hdfVersion == HDFVERSION_5) THEN 
         CALL OutputL3MZ_HE5 (file, mz, createFlag, prodCount)
      ENDIF

!----------------------------------------------
    END SUBROUTINE OutputL3MZ
!----------------------------------------------

!----------------------------------------------
   SUBROUTINE OutputL3MZ_HE2 (file, mz, createFlag)
!----------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3mz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      TYPE( L3MZData_T ), INTENT(IN) :: mz

      LOGICAL, INTENT(INOUT) :: createFlag

! Parameters

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, &
           & swdefgfld, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swfID, swID

! Open the file for appending a swath

      swfID = swopen(file, DFACC_RDWR)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(file)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create a swath of the appropriate name

      swID = swcreate(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to create swath ' // trim(mz%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath dimensions

      status = swdefdim(swID, DIMR_NAME, MIN_MAX)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMR_NAME
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

      status = swdefgfld(swID, GEO_FIELD3, DIMR_NAME, DFNT_FLOAT64, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD4, DIMRL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD12, DIMRL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath data fields using the above dimensions

      status = swdefdfld(swID, DATA_FIELDV, DIMLL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdfld(swID, DATA_FIELDP, DIMLL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface after definition

      status = swdetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
           & 'Failed to detach from swath interface after L3DZ definition.')

! Re-attach to the swath for writing

      swID = swattach(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to swath '//TRIM(mz%name)// ' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = mz%nLevels
      edge(2) = mz%nLats

! Geolocation fields

      status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mz%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mz%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(mz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(mz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = MIN_MAX
      status = swwrfld( swID, GEO_FIELD4, start, stride, edge, &
           & REAL(mz%localSolarTime) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD12, start, stride, edge, &
           & REAL(mz%localSolarZenithAngle) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      edge(1) = mz%nLevels
      status = swwrfld( swID, DATA_FIELDV, start, stride, edge, & 
           & REAL(mz%l3mzValue) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, DATA_FIELDP, start, stride, edge, &
           & REAL(mz%l3mzPrecision) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from swath interface

      status = swdetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
           & 'Failed to detach from swath interface after writing.')

      msr = 'Zonal ' // TRIM(mz%name) // ' successfully written to file ' & 
           & // trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

! Create & write the diagnostic swath

      CALL OutputMZDiag_HE2(file, swfID, mz)

! Close the file

      status = swclose(swfID)
      IF (status == -1) THEN
         msr = 'Failed to close file ' // trim(file) // ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      createFlag = .TRUE.

!---------------------------
   END SUBROUTINE OutputL3MZ_HE2
!---------------------------

!----------------------------------------------
   SUBROUTINE OutputL3MZ_HE5 (file, mz, createFlag, prodCount)
!----------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, HE5T_NATIVE_CHAR

! Brief description of subroutine
! This subroutine creates and writes to the swath in an l3mz file.
! Change to Zonal Average format instead of swath format

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      TYPE( L3MZData_T ), INTENT(IN) :: mz

      LOGICAL, INTENT(INOUT) :: createFlag
 
      INTEGER, INTENT(IN) :: prodCount

! Parameters

! Functions

      INTEGER, EXTERNAL :: he5_zaattach, he5_zaclose, he5_zacreate, & 
           & he5_zadefine, he5_zadefdim, &
           & he5_zadetach, he5_zaopen, he5_zawrite

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swfID, swID

! Open the file for appending a zonal

      IF (prodCount == 1) THEN
         swfID = he5_zaopen(trim(file), HE5F_ACC_TRUNC)
      ELSE
         swfID = he5_zaopen(trim(file), HE5F_ACC_RDWR)
      ENDIF

      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // trim(file)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Create a zonal  of the appropriate name

      swID = he5_zacreate(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to create zonal ' // trim(mz%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the zonal dimensions

      status = he5_zadefdim(swID, DIMR_NAME, MIN_MAX)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMR_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefdim(swID, DIM_NAME2, mz%nLevels)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefdim(swID, DIML_NAME, mz%nLats)
      IF (status == -1) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the zonal geolocation fields using the above dimensions

      status = he5_zadefine(swID, GEO_FIELD3, DIMR_NAME, "", HE5T_NATIVE_DOUBLE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, GEO_FIELD9, DIM_NAME2, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, GEO_FIELD1, DIML_NAME, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, GEO_FIELD4, DIMRL_NAME, "", HE5T_NATIVE_FLOAT)
     IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, GEO_FIELD12, DIMRL_NAME, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the he5_zonal data fields using the above dimensions

      status = he5_zadefine(swID, DATA_FIELDV, DIMLL_NAME, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, DATA_FIELDP, DIMLL_NAME, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the zonal interface after definition

      status = he5_zadetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
           & 'Failed to detach from zonal interface after L3DZ definition.')

! Re-attach to the zonal for writing

      swID = he5_zaattach(swfID, mz%name)
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to zonal '//TRIM(mz%name)// ' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = mz%nLevels
      edge(2) = mz%nLats

! Geolocation fields

      status = he5_zawrite(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mz%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mz%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(mz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(mz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = MIN_MAX
      status = he5_zawrite( swID, GEO_FIELD4, start, stride, edge, &
           & REAL(mz%localSolarTime) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite( swID, GEO_FIELD12, start, stride, edge, &
           & REAL(mz%localSolarZenithAngle) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      edge(1) = mz%nLevels
      status = he5_zawrite( swID, DATA_FIELDV, start, stride, edge, & 
           & REAL(mz%l3mzValue) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite( swID, DATA_FIELDP, start, stride, edge, &
           & REAL(mz%l3mzPrecision) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from zonal interface

      status = he5_zadetach(swID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
           & 'Failed to detach from zonal interface after writing.')

      msr = 'Zonal ' // TRIM(mz%name) // ' successfully written to file ' & 
           & // trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

! Create & write the diagnostic zonal

      CALL OutputMZDiag_HE5(file, swfID, mz)

! Close the file

      status = he5_zaclose(swfID)
      IF (status == -1) THEN
         msr = 'Failed to close file ' // trim(file) // ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      createFlag = .TRUE.

!---------------------------
   END SUBROUTINE OutputL3MZ_HE5
!---------------------------

!-------------------------------------------
   SUBROUTINE OutputMZDiag_HE2 (file, swfID, mz)
!-------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic swaths in an l3mz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      INTEGER, INTENT(IN) :: swfID

      TYPE( L3MZData_T ), INTENT(IN) :: mz

! Parameters

! Functions

      INTEGER, EXTERNAL :: swattach, swcreate, swdefdfld, swdefdim, swdefgfld,&
           & swdetach, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swID

! Create a swath of the appropriate name

      dgName = TRIM(mz%name) // 'Diagnostics'

      swID = swcreate(swfID, dgName)
      IF (swID == -1) THEN
         msr = 'Failed to create swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath dimensions

      status = swdefdim(swID, DIMR_NAME, MIN_MAX)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMR_NAME
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

      status = swdefgfld(swID, GEO_FIELD3, DIMR_NAME, DFNT_FLOAT64, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath data fields using the above dimensions

      status = swdefdfld(swID, DG_FIELD1, DIMLL_NAME, DFNT_FLOAT32, & 
           & HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdfld(swID, DG_FIELD2, DIMLL_NAME, DFNT_INT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface after definition

      status = swdetach(swID)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after L3MZ dg definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the swath for writing

      swID = swattach(swfID, dgName)
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to swath '// TRIM(dgName) //' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = mz%nLevels
      edge(2) = mz%nLats

! Geolocation fields

      status = swwrfld(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mz%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mz%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(mz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(mz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = swwrfld( swID, DG_FIELD1, start, stride, edge, REAL(mz%latRss) )
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, DG_FIELD2, start, stride, edge, mz%perMisPoints)
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from swath interface

      status = swdetach(swID)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath '// TRIM(dgName) // ' successfully written to file ' // & 
           & trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!-----------------------------
   END SUBROUTINE OutputMZDiag_HE2
!-----------------------------

!-------------------------------------------
   SUBROUTINE OutputMZDiag_HE5 (file, swfID, mz)
!-------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, HE5T_NATIVE_CHAR

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic swaths in an l3mz file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: file

      INTEGER, INTENT(IN) :: swfID

      TYPE( L3MZData_T ), INTENT(IN) :: mz

! Parameters

! Functions

      INTEGER, EXTERNAL :: he5_zaattach, he5_zacreate, he5_zadefine, & 
           & he5_zadefdim, he5_zadetach, he5_zawrite

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swID

! Create a zonal of the appropriate name

      dgName = TRIM(mz%name) // 'Diagnostics'

      swID = he5_zacreate(swfID, trim(dgName) )
      IF (swID == -1) THEN
         msr = 'Failed to create zonal ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the zonal dimensions

      status = he5_zadefdim(swID, DIMR_NAME, MIN_MAX)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMR_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefdim(swID, DIM_NAME2, mz%nLevels)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefdim(swID, DIML_NAME, mz%nLats)
      IF (status == -1) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath geolocation fields using the above dimensions

      status = he5_zadefine(swID, GEO_FIELD3, DIMR_NAME, "", HE5T_NATIVE_DOUBLE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!      status = he5_zadefine(swID, GEO_FIELD9, DIM_NAME2, "", HE5T_NATIVE_FLOAT)
!      IF (status == -1) THEN
!         msr = GEO_ERR // GEO_FIELD9
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
!      ENDIF

!      status = he5_zadefine(swID, GEO_FIELD1, DIML_NAME, "", HE5T_NATIVE_FLOAT)
!      IF (status == -1) THEN
!         msr = GEO_ERR // GEO_FIELD1
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
!      ENDIF

! Define the zonal data fields using the above dimensions

      status = he5_zadefine(swID, DG_FIELD1, DIMLL_NAME, "", HE5T_NATIVE_FLOAT)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zadefine(swID, DG_FIELD2, DIMLL_NAME, "", HE5T_NATIVE_INT)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the zonal interface after definition

      status = he5_zadetach(swID)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after L3MZ dg definition.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the zonal for writing

      swID = he5_zaattach(swfID, trim(dgName) )
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to zonal '// TRIM(dgName) //' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = mz%nLevels
      edge(2) = mz%nLats

! Geolocation fields

      status = he5_zawrite(swID, GEO_FIELD3, start(1), stride(1), stride(1), &
           & mz%startTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite(swID, GEO_FIELD3, stride(1), stride(1), stride(1), &
           & mz%endTime)
      IF (status == -1) THEN
          msr = WR_ERR // GEO_FIELD3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!      status = he5_zawrite( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
!           & REAL(mz%pressure) )
!      IF (status == -1) THEN
!         msr = WR_ERR // GEO_FIELD9
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
!      ENDIF

!      status = he5_zawrite( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
!           & REAL(mz%latitude) )
!      IF (status == -1) THEN
!         msr = WR_ERR // GEO_FIELD1
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
!      ENDIF

! Data fields

      status = he5_zawrite( swID, DG_FIELD1, start, stride, edge, & 
           & REAL(mz%latRss) )
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrite( swID, DG_FIELD2, start, stride, edge, & 
           & mz%perMisPoints)
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from zonal interface

      status = he5_zadetach(swID)
      IF (status /= 0) THEN
          msr = SW_ERR // TRIM(dgName) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Zonal '// TRIM(dgName) // ' successfully written to file ' // & 
           & trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!-----------------------------
   END SUBROUTINE OutputMZDiag_HE5
!-----------------------------

!----------------------------------------------------------
   SUBROUTINE WriteMetaL3MZ (fileName, mcfNum, pcf, anText, nFiles, hdfVersion)
!----------------------------------------------------------
   USE PCFHdr, ONLY: WritePCF2Hdr, WriteInputPointer, he5_writeglobalattr, &
	& GlobalAttributes
   USE HDFEOS5, ONLY: HE5F_ACC_RDWR, HE5_GDCLOSE, HE5_GDOPEN
   USE output_m, only: output

! Brief description of subroutine
! This routine writes the metadata for an l3mz file, and annotates it with the
! PCF.

! Arguments

      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      CHARACTER(LEN=*), INTENT(IN) :: fileName

      INTEGER, INTENT(IN) :: mcfNum

      CHARACTER (LEN=1), POINTER :: anText(:)

      INTEGER, INTENT(IN) :: hdfVersion

      INTEGER, INTENT(IN) :: nFiles

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d, &
      & pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write, & 
      & pgs_pc_getUniversalRef, he5_zainqza, he5_zaopen, he5_zaclose

! Variables

      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER (LEN=98) :: inpt(31)
      CHARACTER (LEN=GridNameLen) :: zonalName
      CHARACTER (LEN=FileNameLen) :: sval
      CHARACTER (LEN=45) :: attrName, lvid
      CHARACTER (LEN=3)  :: cNum
      ! CHARACTER (LEN=2)  :: fileType
      integer :: fileType

      REAL(r8) :: dval

      INTEGER :: dg, hdfReturn, i, indx, len, numZonal, pNum, &
           & result, returnStatus, sdid, version
      INTEGER :: daysNum = 1

! calculate the daysNum

     daysNum = nFiles + 1 
     !call output('daysNum = ', advance='no')
     !call output(daysNum, advance='yes')

! Check to see whether this is the Diagnostic files

      dg = 0

      IF ( INDEX(fileName,'Diagnostic') /= 0 ) dg = 1

! Initialize the fileType

      fileType = l_swath

! Initialize the MCF

      result = pgs_met_init(mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Initialization error.  See LogStatus for details.')

! Open the HDF file and initialize the SD interface

      sdid = mls_sfstart(trim(fileName), DFACC_RDWR, hdfVersion=hdfVersion, & 
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
           & 'processed once')
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

      numZonal = 0
      numZonal = he5_zainqza(trim(fileName), list, len) 

!      call output(trim(list), advance='yes')
!      call output('In L3MZ: len = ', advance='no')
!      call output(len, advance='yes')
!      call output('In L3MZ: numZonal = ', advance='no')
!      call output(numZonal, advance='yes')
      IF (numZonal .LE. 0) THEN
         msr = 'No zonal found in file ' // TRIM(fileName) // &
              & ' while attempting to write its metadata.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      pNum = 0

! For each zonal in the file

      DO i = 1, numZonal

! Extract its name

         indx = INDEX(list, ',')
         IF (indx /= 0) THEN
            zonalName = list(:indx-1)
            list = list(indx+1:)
         ELSE
            zonalName = list
         ENDIF

!         call output('zonalName =', advance='no')
!         call output(trim(zonalName), advance='yes')

! If this is a diagnostic zonal, skip to the next one

         IF ( INDEX(zonalName, 'Diagnostics') /= 0) CYCLE

! Append a class suffix to ParameterName, and write the grid name as its value

         pNum = pNum + 1

         IF (pNum < 10) THEN
            WRITE( cNum, '(I1)' ) pNum
         ELSE IF ( (pNum >= 10) .AND. (pNum < 100) ) THEN
            WRITE( cNum, '(I2)' ) pNum
         ELSE
            WRITE( cNum, '(I3)' ) pNum
         ENDIF

         attrName = 'ParameterName' // '.' // cNum
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, & 
              & zonalName)
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
      !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
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
! >          indx = mlspcf_l2dg_start
! >       ELSE
! >          indx = mlspcf_l2gp_start
! >       ENDIF
! >       DO i = 1, 31
! >          version = 1
! >          returnStatus = pgs_pc_getUniversalRef(indx, version, sval)
! >          IF (returnStatus == PGS_S_SUCCESS) THEN
! >             inpt(i) = sval
! >             indx = indx + 1
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
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, 'Limb')
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
           & pcf%startDay)
      IF (result /= PGS_S_SUCCESS) THEN
         msr = METAWR_ERR // attrName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      attrName = 'RangeBeginningTime'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & '00:00:00.000000')
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
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
           & '23:59:59.999999')
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
   
      result = pgs_met_write(groups(INVENTORYMETADATA), "coremetadata", sdid)
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
           &  addingMetaData=.TRUE.)
      IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Error closing HDF file after writing metadata.')

! Annotate the file with the PCF

      ! Write global attributes
      if ( HDFVersion == HDFVERSION_5 ) then
        !print *, 'write global attributes'
        sdid = he5_zaopen (fileName, HE5F_ACC_RDWR)
	if (sdid == -1) then 
            msr = 'Failed to open file before writting global attribute' &
	 	& // trim(fileName)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if
        call he5_writeglobalattr(sdid,daysNum)
        result = he5_zaclose (sdid)
	if (result == -1) then 
            msr = 'Failed to close file after write global attribute' &
		& // trim(fileName)
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
        end if
      endif
                                                                            
      CALL WritePCF2Hdr(fileName, anText, hdfVersion=hdfVersion, & 
           & fileType=fileType)

      result = pgs_met_remove()
      if (result /= PGS_S_SUCCESS .and. WARNIFCANTPGSMETREMOVE) THEN 
         write(msr, *) result
         CALL MLSMessage (MLSMSG_Warning, ModuleName, &
             & "Calling pgs_met_remove() failed with value " // trim(msr) )
      endif

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
         msr = MLSMSG_Allocate // ' pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      ALLOCATE(l3mz%latitude(l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Ancillary data fields

      ALLOCATE(l3mz%localSolarTime(MIN_MAX,l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local solar time pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3mz%localSolarZenithAngle(MIN_MAX,l3mz%nLats), &
               STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local solar zenith angle pointer.'
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

! Diagnostic fields

      ALLOCATE(l3mz%latRss(l3mz%nLevels,l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' latRss pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3mz%perMisPoints(l3mz%nLevels,l3mz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' perMisPoints pointer.'
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

! Vertical geolocation field

      IF ( ASSOCIATED(l3mz%pressure) ) THEN
         DEALLOCATE (l3mz%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Horizontal geolocation field

      IF ( ASSOCIATED(l3mz%latitude) ) THEN
         DEALLOCATE (l3mz%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  latitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Ancillary data fields

      IF ( ASSOCIATED(l3mz%localSolarTime) ) THEN
         DEALLOCATE (l3mz%localSolarTime, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  local solar time pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3mz%localSolarZenithAngle) ) THEN
         DEALLOCATE (l3mz%localSolarZenithAngle, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  local solar zenith angle pointer.'
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

! Diagnostic fields

      IF ( ASSOCIATED(l3mz%latRss) ) THEN
         DEALLOCATE (l3mz%latRss, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  latRss pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3mz%perMisPoints) ) THEN
         DEALLOCATE (l3mz%perMisPoints, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  perMisPoints pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3MZ
!-------------------------------

!==================
END MODULE L3MZData
!==================

! $Log$
! Revision 1.14  2003/09/16 16:35:21  cvuu
! Add new parameter nFiles to the subroutine WriteMetaL3MZ
!
! Revision 1.13  2003/09/15 18:27:23  cvuu
! Output OrbitNumber and OrbitPeriod in the global attribute
!
! Revision 1.12  2003/08/11 23:26:37  cvuu
! brought closer to James Johnson want to
!
! Revision 1.11  2003/07/08 00:17:46  pwagner
! fileType now a lit_name instead of a char string
!
! Revision 1.10  2003/06/02 23:45:15  pwagner
! metadata chnages: OrbitNumber now -1; equatorCrossingDate now utc start date
!
! Revision 1.9  2003/05/30 23:54:07  pwagner
! Relies on lib/PCFHdr to WriteInputPointer
!
! Revision 1.8  2003/04/30 18:16:29  pwagner
! Work-around for LF95 infinite compile-time bug
!
! Revision 1.7  2003/04/06 02:25:50  jdone
! added HDFEOS5 capability
!
! Revision 1.6  2003/03/15 00:20:02  pwagner
! May warn if pgs_met_remove returns non-zero value
!
! Revision 1.5  2001/12/12 17:46:30  nakamura
! Added dg fields; removed unused subroutine DestroyL3MZDatabase.
!
! Revision 1.4  2001/11/12 20:25:50  nakamura
! Added L3MZDiag_T.
!
! Revision 1.3  2001/10/04 18:25:04  nakamura
! Removed lev as dim for local solar fields.
!
! Revision 1.2  2001/09/26 19:47:03  nakamura
! Added local solar ancillary fields.
!
! Revision 1.1  2001/07/18 15:42:09  nakamura
! Module for the L3MZ data type.
!
