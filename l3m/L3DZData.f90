! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE L3DZData
!==============================================================================

   USE Intrinsic, ONLY: l_hdfeos, l_swath
   USE MLSCommon, ONLY: r8, FileNameLen
   USE MLSFiles, ONLY: mls_sfstart, mls_sfend, mls_inqswath, & 
        & HDFVERSION_5, HDFVERSION_4
   USE MLSL3Common, ONLY: GridNameLen, DIMR_NAME, DIMT_NAME, DIM_NAME2, &
        & DIML_NAME, SZ_ERR, DAT_ERR, GEO_ERR, SW_ERR, DG_FIELD1, DG_FIELD2, &
        & DIMLL_NAME, DATE_LEN, INVENTORYMETADATA, OutputFiles_T, &
        & GEO_FIELD1, GEO_FIELD4, GEO_FIELD9, GEO_FIELD11, GEO_FIELD12, &
        & HDFE_NOMERGE, DIM_ERR, MIN_MAX, DIMRL_NAME, WR_ERR, METAWR_ERR
   USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
        & MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_WARNING
   IMPLICIT NONE
   private
   PUBLIC :: L3DZData_T, &
     & OutputL3DZ, ReadL3DZData, WriteMetaL3DZ, &
     & AllocateL3DZ, DeallocateL3DZ, DestroyL3DZDatabase

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L3DZData_T
! Subroutines -- OutputL3DZ
!                OutputDZDiag
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

! This data type is used to store the l3 daily zonal mean data and diagnostics.

   TYPE L3DZData_T

      CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

      CHARACTER (LEN=DATE_LEN) :: date			! day processed

      ! Other ancillary data

      REAL(r8), DIMENSION(:,:), POINTER :: localSolarTime
      REAL(r8), DIMENSION(:,:), POINTER :: localSolarZenithAngle
	! dimensioned as (2, nLats), where 1 = min, 2 = max in dim one

      ! Now the data fields:

      REAL(r8), DIMENSION(:,:), POINTER :: l3dzValue	 ! Field value
      REAL(r8), DIMENSION(:,:), POINTER :: l3dzPrecision ! Field precision
      ! dimensioned as (nLevels, nLats)

      ! Now the diagnostic fields:

      REAL(r8), DIMENSION(:,:), POINTER :: latRss
      ! Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

      ! Now we store the geolocation fields.  First, the vertical one:

      REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

      ! Now the horizontal geolocation information and time:

      REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)

      INTEGER, DIMENSION(:,:), POINTER :: perMisPoints
      ! Missing points (percentage), dimensioned (nLevels, nLats)

      INTEGER :: nLevels			! Total number of surfaces
      INTEGER :: nLats				! Total number of latitudes

   END TYPE L3DZData_T

 CONTAINS

   !------------------------------------------
   SUBROUTINE OutputL3DZ(type, dzm, zFiles, hdfVersion)
   !------------------------------------------

     ! Brief description of subroutine
     ! This subroutine creates and writes to the swaths in an l3dz file.

     ! Arguments

     CHARACTER (LEN=*), INTENT(IN) :: type

     TYPE( L3DZData_T ), POINTER :: dzm(:)

     TYPE( OutputFiles_T ), INTENT(OUT) :: zFiles

     INTEGER, INTENT(IN) :: hdfVersion

     IF (hdfVersion == HDFVERSION_4) THEN 
        CALL  OutputL3DZ_HE2(type, dzm, zFiles)
     ELSE IF (hdfVersion == HDFVERSION_5) THEN
        CALL  OutputL3DZ_HE5(type, dzm, zFiles)
     ENDIF

    !---------------------------
   END SUBROUTINE OutputL3DZ
    !---------------------------
     
   !------------------------------------------
   SUBROUTINE OutputL3DZ_HE2(type, dzm, zFiles)
   !------------------------------------------
   USE Hdf, ONLY: DFNT_FLOAT32, DFNT_CHAR8, DFNT_INT32, &
        & DFACC_RDWR
   USE MLSPCF3, only : mlspcf_l3dz_start, mlspcf_l3dz_end
   USE MLSStrings, ONLY: LinearSearchStringArray
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate 

     ! Brief description of subroutine
     ! This subroutine creates and writes to the swaths in an l3dz file.

     ! Arguments

     CHARACTER (LEN=*), INTENT(IN) :: type

     TYPE( L3DZData_T ), POINTER :: dzm(:)

     TYPE( OutputFiles_T ), INTENT(OUT) :: zFiles

     ! Parameters

! Functions

     INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, &
          & swdefgfld, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: file

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: i, match, numDays, status, swfID, swID

! For each day in the output window,

      numDays = SIZE(dzm)

      DO i = 1, numDays

! Get the name of the L3DZ file for the proper type & day from the PCF

         CALL SearchPCFDates(type, dzm(i)%date, mlspcf_l3dz_start, &
              & mlspcf_l3dz_end, match, file)
         IF (match == -1) THEN
            msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
                 & dzm(i)%date
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Check whether the name is distinct; if so, save it in zFiles

         IF (LinearSearchStringArray(zFiles%name,file) == 0) THEN
            zFiles%nFiles = zFiles%nFiles+1
            zFiles%name(zFiles%nFiles) = trim(file)
            zFiles%date(zFiles%nFiles) = dzm(i)%date
         ENDIF

         ! Open the file for appending a swath

         swfID = swopen(trim(file), DFACC_RDWR)
         IF (swfID == -1) THEN
            msr = MLSMSG_Fileopen // trim(file)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Create a swath of the appropriate name

         swID = swcreate(swfID, trim(dzm(i)%name) )
         IF (swID == -1) THEN
            msr = 'Failed to create swath ' // trim(dzm(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the swath dimensions

         status = swdefdim(swID, DIMR_NAME, MIN_MAX)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMR_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

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

         status = swdefgfld(swID, GEO_FIELD11, DIMT_NAME, &
              & DFNT_CHAR8, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, &
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, & 
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swID, GEO_FIELD4, DIMRL_NAME, &
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swID, GEO_FIELD12, DIMRL_NAME, & 
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD12
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the swath data fields using the above dimensions

         status = swdefdfld(swID, DATA_FIELDV, DIMLL_NAME, &
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDP, DIMLL_NAME, &
              & DFNT_FLOAT32, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Detach from the swath interface after definition

         status = swdetach(swID)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
              & 'Failed to detach from swath interface after L3DZ definition.')
         
         ! Re-attach to the swath for writing

         swID = swattach(swfID, trim(dzm(i)%name) )
         IF (swID == -1) THEN
            msr = 'Failed to re-attach to swath ' // TRIM(dzm(i)%name) &
                 & // ' for writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Write the data

         start = 0
         stride = 1
         edge(1) = DATE_LEN
         edge(2) = dzm(i)%nLats
         
         ! Geolocation fields

         status = swwrfld(swID, GEO_FIELD11, start(1), stride(1), edge(1), &
              & dzm(i)%date)
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         edge(1) = dzm(i)%nLevels
         status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
              & REAL(dzm(i)%pressure) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
              & REAL(dzm(i)%latitude) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = MIN_MAX
         status = swwrfld( swID, GEO_FIELD4, start, stride, edge, &
              & REAL(dzm(i)%localSolarTime) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, GEO_FIELD12, start, stride, edge, &
              & REAL(dzm(i)%localSolarZenithAngle) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD12
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Data fields

         edge(1) = dzm(i)%nLevels
         status = swwrfld( swID, DATA_FIELDV, start, stride, edge, &
              & REAL(dzm(i)%l3dzValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, DATA_FIELDP, start, stride, edge, &
              & REAL(dzm(i)%l3dzPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! After writing, detach from swath interface

         status = swdetach(swID)
         IF (status == -1) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                 & 'Failed to detach from swath interface after writing.')
         ENDIF

         msr = 'Swath ' // TRIM(dzm(i)%name) // &
              & ' successfully written to file ' // trim(file)
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

         ! Define & write the diagnostic swath

         CALL OutputDZDiag_HE2( file, swfID, dzm(i) )

         ! Close the file

         status = swclose(swfID)
         IF (status == -1) THEN
            msr = 'Failed to close file ' // trim(file) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

    !---------------------------
    END SUBROUTINE OutputL3DZ_HE2
    !---------------------------

   !------------------------------------------
   SUBROUTINE OutputL3DZ_HE5(type, dzm, zFiles)
   !------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, HE5T_NATIVE_CHAR
   USE MLSPCF3, only : mlspcf_l3dz_start, mlspcf_l3dz_end
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate 
   USE MLSStrings, ONLY: LinearSearchStringArray
     ! Brief description of subroutine
     ! This subroutine creates and writes to the swaths in an l3dz file.

     ! Arguments

     CHARACTER (LEN=*), INTENT(IN) :: type

     TYPE( L3DZData_T ), POINTER :: dzm(:)

     TYPE( OutputFiles_T ), INTENT(OUT) :: zFiles

     ! Parameters

! Functions

     INTEGER, EXTERNAL :: he5_swattach, he5_swclose, he5_swcreate, & 
          & he5_swdefdfld, he5_swdefdim, &
          & He5_swdefgfld, he5_swdetach, he5_swopen, he5_swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: file

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: i, match, numDays, status, swfID, swID

! For each day in the output window,

      numDays = SIZE(dzm)

      DO i = 1, numDays

! Get the name of the L3DZ file for the proper type & day from the PCF

         CALL SearchPCFDates(type, dzm(i)%date, mlspcf_l3dz_start, &
              & mlspcf_l3dz_end, match, file)
         IF (match == -1) THEN
            msr = 'No ' // TRIM(type) // ' file found in the PCF for day ' // &
                 & dzm(i)%date
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Check whether the name is distinct; if so, save it in zFiles

         IF (LinearSearchStringArray(zFiles%name,file) == 0) THEN
            zFiles%nFiles = zFiles%nFiles+1
            zFiles%name(zFiles%nFiles) = trim(file)
            zFiles%date(zFiles%nFiles) = dzm(i)%date
         ENDIF

         ! Open the file for appending a swath

         swfID = he5_swopen(trim(file), HE5F_ACC_TRUNC)
         IF (swfID == -1) THEN
            msr = MLSMSG_Fileopen // trim(file)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Create a swath of the appropriate name

         swID = he5_swcreate(swfID, trim(dzm(i)%name))
         IF (swID == -1) THEN
            msr = 'Failed to create swath ' // trim(dzm(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the swath dimensions

         status = he5_swdefdim(swID, DIMR_NAME, MIN_MAX)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMR_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefdim(swID, DIMT_NAME, DATE_LEN)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMT_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefdim(swID, DIM_NAME2, dzm(i)%nLevels)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefdim(swID, DIML_NAME, dzm(i)%nLats)
         IF (status == -1) THEN
            msr = DIM_ERR // DIML_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the swath geolocation fields using the above dimensions

         status = he5_swdefgfld(swID, GEO_FIELD11, DIMT_NAME, "", & 
              & HE5T_NATIVE_CHAR, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_swdefgfld(swID, GEO_FIELD9, DIM_NAME2, "", & 
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_swdefgfld(swID, GEO_FIELD1, DIML_NAME, "", &
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefgfld(swID, GEO_FIELD4, DIMRL_NAME, "", &
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefgfld(swID, GEO_FIELD12, DIMRL_NAME, "", &
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD12
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the swath data fields using the above dimensions

         status = he5_swdefdfld(swID, DATA_FIELDV, DIMLL_NAME, "", &
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swdefdfld(swID, DATA_FIELDP, DIMLL_NAME, "", &
              & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Detach from the swath interface after definition

         status = he5_swdetach(swID)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
              & 'Failed to detach from swath interface after L3DZ definition.')
         
         ! Re-attach to the swath for writing

         swID = he5_swattach(swfID, dzm(i)%name)
         IF (swID == -1) THEN
            msr = 'Failed to re-attach to swath ' // TRIM(dzm(i)%name) &
                 & // ' for writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Write the data

         start = 0
         stride = 1
         edge(1) = DATE_LEN
         edge(2) = dzm(i)%nLats
         
         ! Geolocation fields

         status = he5_swwrfld(swID, GEO_FIELD11, start(1), stride(1), edge(1),&
              & dzm(i)%date)
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         edge(1) = dzm(i)%nLevels
         status = he5_swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1),&
              & REAL(dzm(i)%pressure) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2),&
              & REAL(dzm(i)%latitude) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = MIN_MAX
         status = he5_swwrfld( swID, GEO_FIELD4, start, stride, edge,&
              & REAL(dzm(i)%localSolarTime) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swwrfld( swID, GEO_FIELD12, start, stride, edge, &
              & REAL(dzm(i)%localSolarZenithAngle) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD12
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Data fields

         edge(1) = dzm(i)%nLevels
         status = he5_swwrfld( swID, DATA_FIELDV, start, stride, edge, &
              & REAL(dzm(i)%l3dzValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_swwrfld( swID, DATA_FIELDP, start, stride, edge, &
              & REAL(dzm(i)%l3dzPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! After writing, detach from swath interface

         status = he5_swdetach(swID)
         IF (status == -1) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                 & 'Failed to detach from swath interface after writing.')
         ENDIF

         msr = 'Swath ' // TRIM(dzm(i)%name) // &
              & ' successfully written to file ' // trim(file)
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

         ! Define & write the diagnostic swath

         CALL OutputDZDiag_HE5( file, swfID, dzm(i) )

         ! Close the file

         status = he5_swclose(swfID)
         IF (status == -1) THEN
            msr = 'Failed to close file ' // trim(file) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

    !---------------------------
    END SUBROUTINE OutputL3DZ_HE5
    !---------------------------

    !-------------------------------------------
    SUBROUTINE OutputDZDiag_HE2(file, swfID, dz)
    !-------------------------------------------
   USE Hdf, ONLY: DFNT_FLOAT32, DFNT_CHAR8, DFNT_INT32, &
        & DFACC_RDWR

      ! Brief description of subroutine
      ! This subroutine creates and writes to the diagnostic swaths 
      ! in an l3dz file.

      ! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file

      INTEGER, INTENT(IN) :: swfID

      TYPE( L3DZData_T ), INTENT(IN) :: dz

      ! Parameters

      ! Functions

      INTEGER, EXTERNAL :: swattach, swcreate, swdefdfld, swdefdim, &
           & swdefgfld, swdetach, swwrfld

      ! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swID

      ! Create a swath of the appropriate name

      dgName = TRIM(dz%name) // 'Diagnostics'

      swID = swcreate(swfID, trim(dgName) )
      IF (swID == -1) THEN
         msr = 'Failed to create swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Define the swath dimensions

      status = swdefdim(swID, DIMT_NAME, DATE_LEN)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = swdefdim(swID, DIM_NAME2, dz%nLevels)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdim(swID, DIML_NAME, dz%nLats)
      IF (status == -1) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath geolocation fields using the above dimensions

      status = swdefgfld(swID,GEO_FIELD11, DIMT_NAME, DFNT_CHAR8, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID,GEO_FIELD9,DIM_NAME2, DFNT_FLOAT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swID,GEO_FIELD1, DIML_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath data fields using the above dimensions

      status = swdefdfld(swID,DG_FIELD1, DIMLL_NAME, DFNT_FLOAT32,HDFE_NOMERGE)
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
         msr = SW_ERR // TRIM(dgName) // ' after L3DZ dg definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the swath for writing

      swID = swattach(swfID, trim(dgName) )
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to swath ' // TRIM(dgName)//' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats

! Geolocation fields

      status = swwrfld(swID,GEO_FIELD11, start(1), stride(1), edge(1), dz%date)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dz%nLevels
      status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(dz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(dz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = swwrfld( swID, DG_FIELD1, start, stride, edge, REAL(dz%latRss) )
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swID, DG_FIELD2, start, stride, edge, dz%perMisPoints)
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from swath interface

      status = swdetach(swID)
      IF (status /= 0) THEN
         msr = SW_ERR // TRIM(dgName) // ' after L3DZ dg definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath '// TRIM(dgName) // ' successfully written to file '// & 
           & trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!-----------------------------
    END SUBROUTINE OutputDZDiag_HE2
!-----------------------------

    !-------------------------------------------
    SUBROUTINE OutputDZDiag_HE5(file, swfID, dz)
    !-------------------------------------------
   USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, HE5T_NATIVE_CHAR

      ! Brief description of subroutine
      ! This subroutine creates and writes to the diagnostic swaths 
      ! in an l3dz file.

      ! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file

      INTEGER, INTENT(IN) :: swfID

      TYPE( L3DZData_T ), INTENT(IN) :: dz

      ! Parameters

      ! Functions

      INTEGER, EXTERNAL :: he5_swattach, he5_swcreate, he5_swdefdfld, & 
           & he5_swdefdim, he5_swdefgfld, he5_swdetach, he5_swwrfld

      ! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: status, swID

      ! Create a swath of the appropriate name

      dgName = TRIM(dz%name) // 'Diagnostics'

      swID = he5_swcreate(swfID, trim(dgName))
      IF (swID == -1) THEN
         msr = 'Failed to create swath ' // trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Define the swath dimensions

      status = he5_swdefdim(swID, DIMT_NAME, DATE_LEN)
      IF (status == -1) THEN
         msr = DIM_ERR // DIMT_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      status = he5_swdefdim(swID, DIM_NAME2, dz%nLevels)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefdim(swID, DIML_NAME, dz%nLats)
      IF (status == -1) THEN
         msr = DIM_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath geolocation fields using the above dimensions

      status = he5_swdefgfld(swID,GEO_FIELD11, DIMT_NAME, "", & 
           & HE5T_NATIVE_CHAR, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefgfld(swID,GEO_FIELD9,DIM_NAME2, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefgfld(swID,GEO_FIELD1, DIML_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define the swath data fields using the above dimensions

      status = he5_swdefdfld(swID,DG_FIELD1, DIMLL_NAME, "", & 
           & HE5T_NATIVE_FLOAT, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swdefdfld(swID, DG_FIELD2, DIMLL_NAME, "", & 
           & HE5T_NATIVE_INT, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface after definition

      status = he5_swdetach(swID)
      IF (status /= 0) THEN
         msr = SW_ERR // TRIM(dgName) // ' after L3DZ dg definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Re-attach to the swath for writing

      swID = he5_swattach(swfID, trim(dgName) )
      IF (swID == -1) THEN
         msr = 'Failed to re-attach to swath ' // TRIM(dgName)//' for writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Write the data

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats

! Geolocation fields

      status = he5_swwrfld(swID,GEO_FIELD11, start(1), stride(1), edge(1), & 
           & dz%date)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dz%nLevels
      status = he5_swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
           & REAL(dz%pressure) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
           & REAL(dz%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Data fields

      status = he5_swwrfld( swID, DG_FIELD1, start, stride, edge, & 
           & REAL(dz%latRss) )
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_swwrfld(swID, DG_FIELD2, start, stride, edge, & 
           & dz%perMisPoints)
      IF (status == -1) THEN
         msr = WR_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! After writing, detach from swath interface

      status = he5_swdetach(swID)
      IF (status /= 0) THEN
         msr = SW_ERR // TRIM(dgName) // ' after L3DZ dg definition.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      msr = 'Swath '// TRIM(dgName) // ' successfully written to file '// & 
           & trim(file)
      CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

!-----------------------------
   END SUBROUTINE OutputDZDiag_HE5
!-----------------------------

!------------------------------------------------
   SUBROUTINE ReadL3DZData (swfid, swathName, dz, hdfVersion)
!------------------------------------------------

! Brief description of subroutine
! This subroutine reads a data/diagnostic swath pair from a file into the
! L3DZData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: swathName

      INTEGER, INTENT(IN) :: swfid

      TYPE( L3DZData_T ), INTENT(OUT) :: dz

      INTEGER, INTENT(IN), OPTIONAL :: hdfVersion

      IF ( PRESENT(hdfVersion) ) THEN 
         IF (hdfVersion == HDFVERSION_4) THEN 
            CALL ReadL3DZData_HE2(swfid, swathName, dz)
         ELSE IF (hdfVersion == HDFVERSION_5) THEN
            CALL ReadL3DZData_HE5(swfid, swathName, dz)
         ENDIF        
      ELSE
         CALL ReadL3DZData_HE2(swfid, swathName, dz)
      ENDIF

!------------------------------------------------
    END SUBROUTINE ReadL3DZData
!------------------------------------------------

!------------------------------------------------
   SUBROUTINE ReadL3DZData_HE2(swfid, swathName, dz)
!------------------------------------------------

! Brief description of subroutine
! This subroutine reads a data/diagnostic swath pair from a file into the
! L3DZData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: swathName

      INTEGER, INTENT(IN) :: swfid

      TYPE( L3DZData_T ), INTENT(OUT) :: dz

! Parameters

      CHARACTER (LEN=*),PARAMETER :: RL3DZ_ERR = 'Failed to read L3DZ field:  '

! Functions

      INTEGER, EXTERNAL :: & 
           & pgs_td_utcToTAI, swattach, swdetach, swdiminfo, swrdfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: err, nlat, nlev, size, swid, status

      REAL, ALLOCATABLE :: r2(:,:)
      REAL, ALLOCATABLE :: rl(:), rp(:)

! Attach to the data swath for reading

      swid = swattach(swfid, swathName)

      IF (status == -1) THEN
         msr = 'Failed to attach to swath interface for reading ' // &
              & trim(swathName)
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

      ALLOCATE (rp(dz%nLevels),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, rp.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (rl(dz%nLats),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, r1.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (r2(dz%nLevels,dz%nLats),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, r2.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Read the geolocation fields

      dz%name = trim(swathName)

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats


      status = swrdfld(swid, GEO_FIELD11, start(1), stride(1), edge(1), & 
           & dz%date)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dz%nLevels
      status = swrdfld(swid, GEO_FIELD9, start(1), stride(1), edge(1), & 
           & rp)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%pressure = DBLE(rp)

      status = swrdfld(swid, GEO_FIELD1, start(2), stride(2), edge(2), & 
           & rl)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latitude = DBLE(rl)

      edge(1) = MIN_MAX
      status = swrdfld(swid, GEO_FIELD4, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%localSolarTime = DBLE(r2)

      status = swrdfld(swid, GEO_FIELD12, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%localSolarZenithAngle = DBLE(r2)

      ! Read the data fields

      edge(1) = dz%nLevels
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
              & trim(swathName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Attach to the diagnostic swath for reading

      dgName = TRIM(swathName) // 'Diagnostics'
      swid = swattach(swfid, trim(dgName))
      IF (status == -1) THEN
         msr = 'Failed to attach to swath interface for reading ' // & 
              & trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Read the diagnostic fields

      status = swrdfld(swid, DG_FIELD1, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latRss = DBLE(r2)

      status = swrdfld(swid, DG_FIELD2, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%perMisPoints = DBLE(r2)

      !  After reading, detach from swath interface

      status = swdetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from swath interface after reading ' // & 
              & trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Deallocate local variables

      DEALLOCATE(rp, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, rp.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE(rl, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, r1.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE(r2, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, r2.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
    !-----------------------------
    END SUBROUTINE ReadL3DZData_HE2
    !-----------------------------

!------------------------------------------------
    SUBROUTINE ReadL3DZData_HE5(swfid, swathName, dz)
!------------------------------------------------

! Brief description of subroutine
! This subroutine reads a data/diagnostic swath pair from a file into the
! L3DZData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: swathName

      INTEGER, INTENT(IN) :: swfid

      TYPE( L3DZData_T ), INTENT(OUT) :: dz

! Parameters

      CHARACTER (LEN=*),PARAMETER :: RL3DZ_ERR = 'Failed to read L3DZ field:  '

! Functions

      INTEGER, EXTERNAL :: & 
           & pgs_td_utcToTAI, &
           & he5_swattach, he5_swdetach, he5_swdiminfo, he5_swrdfld

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: err, nlat, nlev, size, swid, status

      REAL, ALLOCATABLE :: r2(:,:)
      REAL, ALLOCATABLE :: rl(:), rp(:)

! Attach to the data swath for reading

      swid = he5_swattach(swfid, swathName)

      IF (status == -1) THEN
         msr = 'Failed to attach to swath interface for reading ' // &
              & trim(swathName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get dimension information

      size = he5_swdiminfo(swid, DIM_NAME2)

      IF (size == -1) THEN
         msr = SZ_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlev = size

      size = he5_swdiminfo(swid, DIML_NAME)

      IF (size == -1) THEN
         msr = SZ_ERR // DIML_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlat = size

! Allocate the output structure and temporary variables

      CALL AllocateL3DZ(nlev, nlat, dz)

      ALLOCATE(rp(dz%nLevels),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, rp.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(rl(dz%nLats),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, r1.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(r2(dz%nLevels,dz%nLats),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  local real variable, r2.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Read the geolocation fields

      dz%name = trim(swathName)

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats

      status = he5_swrdfld(swid, GEO_FIELD11, start(1), stride(1), edge(1), & 
           & dz%date)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dz%nLevels
      status = he5_swrdfld(swid, GEO_FIELD9, start(1), stride(1), edge(1), & 
           & rp)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%pressure = DBLE(rp)

      status = he5_swrdfld(swid, GEO_FIELD1, start(2), stride(2), edge(2), & 
           & rl)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latitude = DBLE(rl)

      edge(1) = MIN_MAX
      status = he5_swrdfld(swid, GEO_FIELD4, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%localSolarTime = DBLE(r2)

      status = he5_swrdfld(swid, GEO_FIELD12, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD12
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%localSolarZenithAngle = DBLE(r2)

      ! Read the data fields

      edge(1) = dz%nLevels
      status = he5_swrdfld(swid, DATA_FIELDV, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%l3dzValue = DBLE(r2)

      status = he5_swrdfld(swid, DATA_FIELDP, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_FIELDP
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%l3dzPrecision = DBLE(r2)

      !  After reading, detach from swath interface
      
      status = he5_swdetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from swath interface after reading ' // &
              & trim(swathName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Attach to the diagnostic swath for reading

      dgName = TRIM(swathName) // 'Diagnostics'
      swid = he5_swattach(swfid, dgName)
      IF (status == -1) THEN
         msr = 'Failed to attach to swath interface for reading ' // & 
              & trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Read the diagnostic fields

      status = he5_swrdfld(swid, DG_FIELD1, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latRss = DBLE(r2)

      status = he5_swrdfld(swid, DG_FIELD2, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%perMisPoints = DBLE(r2)

      !  After reading, detach from swath interface

      status = he5_swdetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from swath interface after reading ' // & 
              & trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Deallocate local variables

      DEALLOCATE(rp, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, rp.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE(rl, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, r1.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DEALLOCATE(r2, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local real variable, r2.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
    !-----------------------------
    END SUBROUTINE ReadL3DZData_HE5
    !-----------------------------

    !-------------------------------------------------------
    SUBROUTINE WriteMetaL3DZ (pcf, mcfNum, files, anText, hdfVersion)
    !-------------------------------------------------------
   USE Hdf, ONLY: DFNT_FLOAT32, DFNT_CHAR8, DFNT_INT32, &
        & DFACC_RDWR
   USE MLSPCF3
   USE PCFHdr, ONLY: WritePCF2Hdr, WriteInputPointer
   USE mon_Open, ONLY: PCFMData_T 
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate 
   USE SDPToolkit, only: WARNIFCANTPGSMETREMOVE, PGSD_MET_GROUP_NAME_L, &
        & PGSD_MET_NUM_OF_GROUPS, PGSMET_E_MAND_NOT_SET, PGS_S_SUCCESS

      ! Brief description of subroutine
      ! This routine writes the metadata for an l3dz file, 
      ! and annotates it with the PCF.

      ! Arguments

      TYPE (OutputFiles_T), INTENT(IN) :: files
      
      TYPE( PCFMData_T ), INTENT(IN) :: pcf

      INTEGER, INTENT(IN) :: mcfNum

      CHARACTER (LEN=1), POINTER :: anText(:)

      INTEGER, INTENT(IN) :: hdfVersion

      ! Parameters

      ! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d, &
           & pgs_met_setAttr_i,pgs_met_setAttr_s, pgs_met_write, &
           & pgs_pc_getUniversalRef

      ! Variables

      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER (LEN=GridNameLen):: swathName
      CHARACTER (LEN=FileNameLen):: sval
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=45) :: attrName, lvid
      CHARACTER (LEN=3) :: cNum
      ! CHARACTER (LEN=2) :: fileType
      integer :: fileType

      REAL(r8) :: dval

      INTEGER :: dg, hdfReturn, i, j, indx, len, numSwaths, pNum, &
           & result, returnStatus, sdid, version

! Check to see whether these are Diagnostic files

      dg = 0

      IF ( INDEX(files%name(1),'Diagnostic') /= 0 ) dg = 1

! Initialize the MCF

      result = pgs_met_init(mcfNum, groups)
      IF (result /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Initialization error.  See LogStatus for details.')

! Initialize the file type
      fileType = l_swath

! For each file successfully created,

      DO i = 1, files%nFiles

! Open the HDF file and initialize the SD interface

         sdid = mls_sfstart(trim(files%name(i)), DFACC_RDWR, & 
              & hdfVersion=hdfVersion, addingMetaData=.TRUE.)
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
         CALL ExpandFileTemplate('$cycle', lvid, cycle=pcf%cycle)
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, lvid)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! MeasuredParameterContainer -- find the number of swaths in the file

         numSwaths = 0
         numSwaths = mls_inqswath(trim(files%name(i)), list, len, & 
              & hdfVersion=hdfVersion)

         IF (numSwaths .LE. 0) THEN
            msr = 'No swaths found in file ' // TRIM( files%name(i) ) // &
                 & ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         pNum = 0

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

! If this is a diagnostic swath, skip to the next one

            IF ( INDEX(swathName, 'Diagnostics') /= 0 ) CYCLE

! Append a class suffix to ParameterName, and write the grid name as its value

            pNum = pNum + 1

            IF (pNum < 10) THEN
               WRITE( cNum, '(I1)' ) pNum
            ELSE IF ( (pNum >= 10) .AND. (pNum < 100) ) THEN
               WRITE( cNum, '(I2)' ) pNum
            ELSE
               WRITE( cNum, '(I3)' ) pNum
            ENDIF

            attrName = 'ParameterName' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & swathName)
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ! QAFlags Group

            attrName = 'AutomaticQualityFlag' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Passed')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
            
            attrName = 'AutomaticQualityFlagExplanation' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'pending algorithm update')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'OperationalQualityFlag' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Not Investigated')
            IF (result /= PGS_S_SUCCESS) THEN
               msr = METAWR_ERR // attrName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            attrName = 'OperationalQualityFlagExplanation' // '.' // TRIM(cNum)
            result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, &
                 & 'Not Investigated')
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
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         attrName = 'StopOrbitNumber' // '.1'
         !result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, -1)
         result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), attrName, 99999)
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
         ! sval = 'See the PCF annotation to this file.'
! >          IF (dg == 1) THEN
! >             DO j = mlspcf_l2dg_start, mlspcf_l2dg_end
! >                version = 1
! >                returnStatus = pgs_pc_getUniversalRef(j, version, sval)
! >                IF (returnStatus /= PGS_S_SUCCESS) CYCLE
! >                IF ( INDEX(sval,files%date(i)) /= 0 ) EXIT
! >             ENDDO
! >          ENDIF
         
         ! result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName, sval)
         result = WriteInputPointer(groups(INVENTORYMETADATA), attrName, &
           & fileType=l_hdfeos)
         IF (result /= PGS_S_SUCCESS) THEN
            msr = METAWR_ERR // attrName
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Locality Value

         attrName = 'LocalityValue'
         result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), attrName,'Limb')
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
                    & 'Metadata write failed.')
            ENDIF
         ENDIF
         
! Terminate access to the SD interface and close the file
         
         hdfReturn = mls_sfend(sdid, hdfVersion=hdfVersion, & 
              & addingMetaData=.TRUE.)
         IF (hdfReturn /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
              & 'Error closing HDF file after writing metadata.')

         ! Annotate the file with the PCF

         CALL WritePCF2Hdr(files%name(i), anText, hdfVersion=hdfVersion, & 
              & fileType=fileType)

      ENDDO

      result = pgs_met_remove()
      
      if (result /= PGS_S_SUCCESS .and. WARNIFCANTPGSMETREMOVE) THEN 
        write(msr, *) result
        CALL MLSMessage (MLSMSG_Warning, ModuleName, &
             & "Calling pgs_met_remove() failed with value " // trim(msr) )
      endif          

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
         msr = MLSMSG_Allocate // ' pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Horizontal geolocation field

      ALLOCATE(l3dz%latitude(l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Ancillary data fields

      ALLOCATE(l3dz%localSolarTime(MIN_MAX,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local solar time pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ALLOCATE(l3dz%localSolarZenithAngle(MIN_MAX,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local solar zenith angle  pointer.'
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

      ! Diagnostic fields

      ALLOCATE(l3dz%latRss(l3dz%nLevels,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' latRss pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dz%perMisPoints(l3dz%nLevels,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' perMisPoints pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
    !-----------------------------
    END SUBROUTINE AllocateL3DZ
    !-----------------------------

    !----------------------------------
    SUBROUTINE DeallocateL3DZ (l3dz)
    !----------------------------------

      ! Brief description of subroutine
      ! This subroutine deallocates the internal field pointers of L3DZData_T
      ! derived type, after the calling program has finished with the data.

      ! Arguments

      TYPE( L3DZData_T ), INTENT(INOUT) :: l3dz

      ! Parameters

      ! Functions

      ! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err
      
      ! Vertical geolocation field

      IF ( ASSOCIATED(l3dz%pressure) ) THEN
         DEALLOCATE (l3dz%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' l3dz pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF
      
      ! Horizontal geolocation field
      
      IF ( ASSOCIATED(l3dz%latitude) ) THEN
         DEALLOCATE (l3dz%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' l3dz latitude pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      ! Ancillary data fields

      IF ( ASSOCIATED(l3dz%localSolarTime) ) THEN
         DEALLOCATE (l3dz%localSolarTime, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' l3dz local solar time pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF
      
      IF ( ASSOCIATED(l3dz%localSolarZenithAngle) ) THEN
         DEALLOCATE (l3dz%localSolarZenithAngle, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate //' l3dz local solar zenith angle pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      ! Data fields

      IF ( ASSOCIATED(l3dz%l3dzValue) ) THEN
         DEALLOCATE (l3dz%l3dzValue, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' l3dzValue pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF
      
      IF ( ASSOCIATED(l3dz%l3dzPrecision) ) THEN
         DEALLOCATE (l3dz%l3dzPrecision, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' l3dzPrecision pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      ! Diagnostic fields

      IF ( ASSOCIATED(l3dz%latRss) ) THEN
         DEALLOCATE (l3dz%latRss, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' latRss pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF
      
      IF ( ASSOCIATED(l3dz%perMisPoints) ) THEN
         DEALLOCATE (l3dz%perMisPoints, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' perMisPoints pointer.'
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
     ! This subroutine deallocates the internal structures of an l3dz database,
     ! and then database itself


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
           msr = MLSMSG_DeAllocate // ' l3dz database'
           CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
        ENDIF

     ENDIF

   !------------------------------------
   END SUBROUTINE DestroyL3DZDatabase
   !------------------------------------

 !==================
 END MODULE L3DZData
 !==================

! $Log$
! Revision 1.11  2003/07/08 00:17:45  pwagner
! fileType now a lit_name instead of a char string
!
! Revision 1.10  2003/06/02 23:45:15  pwagner
! metadata chnages: OrbitNumber now -1; equatorCrossingDate now utc start date
!
! Revision 1.9  2003/05/30 23:54:07  pwagner
! Relies on lib/PCFHdr to WriteInputPointer
!
! Revision 1.8  2003/04/30 18:16:28  pwagner
! Work-around for LF95 infinite compile-time bug
!
! Revision 1.7  2003/04/06 02:25:34  jdone
! added HDFEOS5 capability
!
! Revision 1.6  2003/03/15 00:20:02  pwagner
! May warn if pgs_met_remove returns non-zero value
!
! Revision 1.5  2001/12/12 17:43:20  nakamura
! Added dg fields.
!
! Revision 1.4  2001/11/12 20:23:55  nakamura
! Added L3DZDiag_T.
!
! Revision 1.3  2001/10/04 18:25:59  nakamura
! Removed lev as dim for local solar fields.
!
! Revision 1.2  2001/09/26 19:46:14  nakamura
! Added local solar ancillary fields.
!
! Revision 1.5  2001/05/04 18:34:50  nakamura
! Combined WriteMeta subroutines; removed L3DZFiles_T; truncated TIME from DATE on output.
!
! Revision 1.4  2001/04/24 19:38:50  nakamura
! Removed references to private L2 parameters.
!
! Revision 1.3  2001/03/27 19:29:13  nakamura
! Updated the metadata; fixed err checks on deallocate.
!
! Revision 1.2  2001/02/21 21:02:50  nakamura
! Changed MLSPCF to MLSPCF3; shifted some parameters; added ReadL3DZData; allowed for expanded swath list; changed InputPointer.
!
! Revision 1.1  2001/02/12 19:23:42  nakamura
! Module for the L3DZ data type.
!
