! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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
        & GEO_FIELD1, GEO_FIELD4, GEO_FIELD9, GEO_FIELD11, GEO_FIELD5, &
        & HDFE_NOMERGE, DIM_ERR, MIN_MAX, DIMRL_NAME, WR_ERR, METAWR_ERR, &
	& FILEATTR_ERR
   USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
        & MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_WARNING
   IMPLICIT NONE
   private
   PUBLIC :: L3DZData_T, &
     & OutputL3DZ, ReadL3DZData, WriteMetaL3DZ, &
     & AllocateL3DZ, DeallocateL3DZ, DestroyL3DZDatabase


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Contents:

! Definition -- L3DZData_T
! Subroutines -- OutputL3DZ
!                ReadL3DZData
!                WriteMetaL3DZ
!		 WriteAttributes_l3dz
!                SetAlias_l3dz
!                AllocateL3DZ
!                DeallocateL3DZ
!                DestroyL3DZDatabase

! Remarks:  This module contains the definition of the L3DZData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDV = 'L3dzValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDP = 'L3dzPrecision'
   CHARACTER (LEN=*), PARAMETER :: DATA_STDDEV = 'L3dzStdDeviation'
   CHARACTER (LEN=*), PARAMETER :: DATA_COUNT = 'L3dzDataCount'

! This data type is used to store the l3 daily zonal mean data and diagnostics.

   TYPE L3DZData_T

      CHARACTER (LEN=GridNameLen) :: name	! name for the output quantity

      CHARACTER (LEN=DATE_LEN) :: date			! day processed

      ! Other ancillary data

      REAL(r8), DIMENSION(:,:), POINTER :: localSolarTime
      REAL(r8), DIMENSION(:,:), POINTER :: SolarZenithAngle
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

      INTEGER, DIMENSION(:,:), POINTER :: dataCount
      ! Data count, dimensioned (nLevels, nLats)
	
      INTEGER :: nLevels			! Total number of surfaces
      INTEGER :: nLats				! Total number of latitudes

   END TYPE L3DZData_T

 CONTAINS

!------------------------------------------
 SUBROUTINE OutputL3DZ(type, dzm, zFiles, creationFlag, hdfVersion)
!------------------------------------------

     ! Brief description of subroutine
     ! This subroutine creates and writes to the swaths in an l3dz file.

     ! Arguments

     CHARACTER (LEN=*), INTENT(IN) :: type

     TYPE( L3DZData_T ), POINTER :: dzm(:)

     TYPE( OutputFiles_T ), INTENT(OUT) :: zFiles

     INTEGER, INTENT(IN) :: hdfVersion

     LOGICAL, INTENT(INOUT) :: creationFlag

     IF (hdfVersion == HDFVERSION_5) THEN
        CALL  OutputL3DZ_HE5(type, dzm, zFiles, creationFlag)
     ELSE
        CALL MLSMessage(MLSMSG_Error, ModuleName, 'Wrong hdfVersion')
     ENDIF

!---------------------------
 END SUBROUTINE OutputL3DZ
!---------------------------
     
!------------------------------------------
 SUBROUTINE OutputL3DZ_HE5(type, dzm, zFiles, creationFlag)
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

     LOGICAL, INTENT(INOUT) :: creationFlag

     ! Parameters

! Functions

     INTEGER, EXTERNAL :: he5_zaattach, he5_zaclose, he5_zacreate, &
           & he5_zadefine, he5_zadefdim, &
           & he5_zadetach, he5_zaopen, he5_zawrite

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: file
      CHARACTER (LEN=7) :: dateChar
      

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: i, match, numDays, status, swfID, swID, dateInt

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

         ! Open the file for appending a zonal

         IF (.not. creationFlag) THEN
            swfID = he5_zaopen(trim(file), HE5F_ACC_TRUNC)
         ELSE
            swfID = he5_zaopen(trim(file), HE5F_ACC_RDWR)
         ENDIF

         IF (swfID == -1) THEN
            msr = MLSMSG_Fileopen // trim(file)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Create a zonal of the appropriate name

         swID = he5_zacreate(swfID, trim(dzm(i)%name))
         IF (swID == -1) THEN
            msr = 'Failed to create swath ' // trim(dzm(i)%name)
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the zonal dimensions

         status = he5_zadefdim(swID, DIMR_NAME, MIN_MAX)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMR_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefdim(swID, DIMT_NAME, 1)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMT_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefdim(swID, DIM_NAME2, dzm(i)%nLevels)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefdim(swID, DIML_NAME, dzm(i)%nLats)
         IF (status == -1) THEN
            msr = DIM_ERR // DIML_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the zonal geolocation fields using the above dimensions

         status = he5_zadefine(swID, GEO_FIELD11, DIMT_NAME, " ", HE5T_NATIVE_INT)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_zadefine(swID, GEO_FIELD9, DIM_NAME2, "",  HE5T_NATIVE_FLOAT)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         
         status = he5_zadefine(swID, GEO_FIELD1, DIML_NAME, "", HE5T_NATIVE_FLOAT)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefine(swID, GEO_FIELD4, DIMRL_NAME, "",  HE5T_NATIVE_FLOAT)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefine(swID, GEO_FIELD5, DIMRL_NAME, "", HE5T_NATIVE_FLOAT)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Define the zonal data fields using the above dimensions

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

         status = he5_zadefine(swID, DATA_STDDEV, DIMLL_NAME, "", HE5T_NATIVE_FLOAT)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_STDDEV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefine(swID, DATA_COUNT, DIMLL_NAME, "", HE5T_NATIVE_INT)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_COUNT
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zadefine(swID, DG_FIELD2, DIMLL_NAME, "", HE5T_NATIVE_INT)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_COUNT
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Write the data

         start = 0
         stride = 1
         edge(1) = dzm(i)%nLevels
         edge(2) = dzm(i)%nLats
       
         ! Convert the date from char to int.  There is problem using char string from 
	 ! toolkit from version 5.2.10 

	 dateChar = dzm(i)%date(1:4) // dzm(i)%date(6:8)
	 read (dateChar, '(i7)') dateInt

         ! Geolocation fields

         status = he5_zawrite(swID, GEO_FIELD11, start(1), stride(1), edge(1), dateInt)
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD11
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         

         status = he5_zawrite( swID, GEO_FIELD9, start(1), stride(1), edge(1),&
              & REAL(dzm(i)%pressure) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, GEO_FIELD1, start(2), stride(2), edge(2),&
              & REAL(dzm(i)%latitude) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = MIN_MAX
         status = he5_zawrite( swID, GEO_FIELD4, start, stride, edge,&
              & REAL(dzm(i)%localSolarTime) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD4
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, GEO_FIELD5, start, stride, edge, &
              & REAL(dzm(i)%SolarZenithAngle) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Data fields

         edge(1) = dzm(i)%nLevels
         edge(2) = dzm(i)%nLats
         status = he5_zawrite( swID, DATA_FIELDV, start, stride, edge, &
              & REAL(dzm(i)%l3dzValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, DATA_FIELDP, start, stride, edge, &
              & REAL(dzm(i)%l3dzPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, DATA_STDDEV, start, stride, edge, &
              & REAL(dzm(i)%latRss) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_STDDEV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, DATA_COUNT, start, stride, edge, &
              & dzm(i)%dataCount )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_COUNT
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = he5_zawrite( swID, DG_FIELD2, start, stride, edge,&
              & REAL(dzm(i)%perMisPoints) )
         IF (status == -1) THEN
            msr = WR_ERR // DG_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! After writing, detach from zonal interface

         status = he5_zadetach(swID)
         IF (status == -1) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, & 
                 & 'Failed to detach from zonal interface after writing.')
         ENDIF

         msr = 'Zonal ' // TRIM(dzm(i)%name) // &
              & ' successfully written to file ' // trim(file)
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

         ! Close the file

         status = he5_zaclose(swfID)
         IF (status == -1) THEN
            msr = 'Failed to close file ' // trim(file) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         ! Write local attributes and set alias
                                                                                
         CALL WriteAttributes_l3dz( trim(file), dzm(i) )
         CALL SetAlias_l3dz( trim(file), dzm(i) )

      ENDDO
      creationFlag = .TRUE.

!--------------------------------
    END SUBROUTINE OutputL3DZ_HE5
!---------------------------------

!----------------------------------------------
    SUBROUTINE SetAlias_l3dz(file, dzm)
!----------------------------------------------
    USE HDF5, ONLY: HID_T
    USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, MLS_charType
    USE SDPToolkit, ONLY: PGS_S_SUCCESS
                                                                                
    ! Brief description of subroutine
    ! This subroutine writes attributes to each zonal field
                                                                                
    ! Arguments
                                                                                
      CHARACTER (LEN=*), INTENT(IN) :: file
                                                                                
      TYPE (L3DZData_T), INTENT(IN) :: dzm
                                                                                
    ! Functions
                                                                                
      INTEGER, EXTERNAL :: &
        & he5_zaattach, he5_zaclose, he5_zadetach, he5_zaopen, &
        & he5_zasetalias
                                                                                
    ! Parameters
                                                                                
      CHARACTER (len=*), parameter :: L3VALUE = 'L3dzValue'
      CHARACTER (len=*), parameter :: L3PREC = 'L3dzPrecision'
      CHARACTER (len=*), parameter :: L3STDDEV = 'L3dzStdDeviation'
      CHARACTER (len=*), parameter :: L3DCOUNT = 'L3dzDataCount'
      CHARACTER (LEN=480) :: msr
      INTEGER (HID_T) :: zafID, zaId
      INTEGER :: status, n, m, field
                                                                                
      zafID = he5_zaopen(trim(file), HE5F_ACC_RDWR)
                                                                                
      IF (zafID == -1) THEN
         msr = MLSMSG_Fileopen // trim(file) //' for writing zonal'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                                                                                
      zaId = he5_zaattach(zafID, dzm%name)
      IF (zaId == -1) THEN
         msr = 'Failed to attach to zonal ' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zasetalias(zaId, L3VALUE, trim(dzm%name))
      if (status /= PGS_S_SUCCESS) then
         msr = 'Failed to alias l3Value' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      endif

      status = he5_zasetalias(zaId, L3PREC, trim(dzm%name)//'Precision')
      if (status /= PGS_S_SUCCESS) then
         msr = 'Failed to alias l3Precision' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      endif

      status = he5_zasetalias(zaId, L3STDDEV, trim(dzm%name)//'StdDeviation')
      if (status /= PGS_S_SUCCESS) then
         msr = 'Failed to alias l3StdDeviation' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      endif

      status = he5_zasetalias(zaId, L3DCOUNT, trim(dzm%name)//'DataCount')
      if (status /= PGS_S_SUCCESS) then
         msr = 'Failed to alias l3DataCount' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      endif

      status = he5_zadetach(zaID)
      IF (status == -1) THEN
         msr = 'Failed to deattach to zonal ' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zaclose(zafID)
      IF (status == -1) THEN
          msr = 'Failed to close file ' // trim(file) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
!----------------------------
 END SUBROUTINE SetAlias_l3dz
!----------------------------
                                                                                
!----------------------------------------------
 SUBROUTINE WriteAttributes_l3dz(file, dzm)
!----------------------------------------------
    USE HDF5, ONLY: HID_T
    USE HDFEOS5, ONLY: HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE, &
       & HE5F_ACC_RDWR, HE5F_ACC_TRUNC, MLS_CHARTYPE
    USE MLSStrings, only: lowercase
    USE MLSStringLists, only: list2array
                                                                                
    ! Brief description of subroutine
    ! This subroutine writes attributes to each grid field
                                                                                
    ! Arguments
      CHARACTER (LEN=*), INTENT(IN) :: file
                                                                                
      TYPE (L3DZData_T), INTENT(IN) :: dzm
                                                                                
      ! Parameters
                                                                                
      integer, parameter :: CHARATTRLEN = 255
      integer, parameter :: NumOfZonalFields = 5
      integer, parameter :: NumOfDataFields = 5
      character (len=*), parameter :: ZonalFieldTitles = &
        & 'Latitude,Pressure,Date,LocalSolarTime,SolarZenithAngle'
      character (len=*), parameter :: ZonalDataTitles = &
        & 'L3dzStdDeviation,L3dzValue,L3dzPrecision,L3dzDataCount,PerMisPoints'
      character(len=CHARATTRLEN), dimension(NumOfZonalFields) :: theTitles
      character(len=CHARATTRLEN), dimension(NumOfDataFields) :: dataTitles
      character(len=CHARATTRLEN) :: units_name, field_name
      real, parameter    :: UNDEFINED_VALUE =  -999.999
                                                                                
      ! Variables
                                                                                
      CHARACTER (LEN=480) :: msr
      CHARACTER (len=*), parameter :: AURA_FIELD = 'Aura-shared '
      CHARACTER (len=*), parameter :: MLS_SHARED_FIELD = 'HIRDL-MLS-TES-shared '
      CHARACTER (len=*), parameter :: MLS_FIELD = 'MLS-Specific '
      CHARACTER (len=*), parameter :: COORDINATE = 'VerticalCoordinate'
      CHARACTER (len=*), parameter :: ZONAL_SPACING = 'ZonalSpacing'
      CHARACTER (len=*), parameter :: ZSPACING_VALUE = '2'
      CHARACTER (len=*), parameter :: ZONAL_SPACING_UNIT = 'ZonalSpacingUnit'
      CHARACTER (len=*), parameter :: ZSPACINGUNIT_VALUE = 'Degree'
      CHARACTER (len=*), parameter :: PRESSURE = 'Pressure'

      INTEGER (HID_T) :: zafID, zaId
      INTEGER :: status, n, m, field
                                                                                
      ! Functions
                                                                                
     INTEGER, EXTERNAL :: he5_zaattach, he5_zaclose, he5_zadetach, &
        & he5_zaopen, he5_zawrlattr, he5_zawrattr
                                                                                
      call List2Array(ZonalFieldTitles, theTitles, .true.)
      call List2Array(ZonalDataTitles, dataTitles, .true.)
                      
      zafID = he5_zaopen(trim(file), HE5F_ACC_RDWR)
      IF (zafID == -1) THEN
         msr = MLSMSG_Fileopen // trim(file) //' for writing zonal'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                                                                                
      zaId = he5_zaattach(zafID, dzm%name)
      IF (zaID == -1) THEN
         msr = 'Failed to attach to zonal ' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Write File Level Attributes
                                                                                
      status = he5_zawrattr(zaId, COORDINATE, MLS_CHARTYPE, &
           &  len_trim(PRESSURE), PRESSURE)
      IF (status /= 0) THEN
         msr = FILEATTR_ERR // COORDINATE
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = he5_zawrattr(zaId, 'Pressure', HE5T_NATIVE_FLOAT, dzm%nLevels, &
          & REAL(dzm%Pressure))
      IF (status /= 0) THEN
         msr = FILEATTR_ERR // 'Pressure'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                                                                                
      status = he5_zawrattr(zaId, ZONAL_SPACING, MLS_CHARTYPE, &
           &  len_trim(ZSPACING_VALUE), ZSPACING_VALUE)
      IF (status /= 0) THEN
         msr = FILEATTR_ERR // ZONAL_SPACING
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                                                                                
      status = he5_zawrattr(zaId, ZONAL_SPACING_UNIT, MLS_CHARTYPE, &
           &  len_trim(ZSPACINGUNIT_VALUE), ZSPACINGUNIT_VALUE)
      IF (status /= 0) THEN
         msr = FILEATTR_ERR // ZONAL_SPACING_UNIT
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                          
      ! Write attributes for zonal (geolocation) fields
                                                                                
      do field=1, NumOfZonalFields
        status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Missing Value', &
           & HE5T_NATIVE_FLOAT, 1, UNDEFINED_VALUE)
        status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Title', &
           & MLS_CHARTYPE, len_trim(theTitles(field)), theTitles(field))
        if (trim(theTitles(field)) == 'Pressure') then
           status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Unit', &
             & MLS_CHARTYPE, 3, 'hPa')
           status = he5_zawrlattr(zaId, trim(theTitles(field)), &
             & 'UniqueFieldDefinition', MLS_CHARTYPE, len_trim(AURA_FIELD), &
             & AURA_FIELD)
        else if (trim(theTitles(field)) == 'Date') then
           status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Unit', &
             & MLS_CHARTYPE, 7, 'NoUnits')
           status = he5_zawrlattr(zaId, trim(theTitles(field)), &
             & 'UniqueFieldDefinition', MLS_CHARTYPE, len_trim(MLS_FIELD), &
             & MLS_FIELD)
        else if (trim(theTitles(field)) == 'LocalSolarTime') then
           status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Unit', &
             & MLS_CHARTYPE, 1, 'h')
           status = he5_zawrlattr(zaId, trim(theTitles(field)), &
             & 'UniqueFieldDefinition', MLS_CHARTYPE, len_trim(MLS_FIELD), &
             & MLS_FIELD)
        else
           status = he5_zawrlattr(zaId, trim(theTitles(field)), 'Unit', &
             & MLS_CHARTYPE, 3, 'deg')
           status = he5_zawrlattr(zaId, trim(theTitles(field)), &
             & 'UniqueFieldDefinition', MLS_CHARTYPE, len_trim(MLS_SHARED_FIELD), &
             & MLS_SHARED_FIELD)
        endif
      enddo
                                                                                
      ! Write attributes for data fields
                                                                                
        select case (trim(lowercase(dzm%name(1:3))))
        case ('tem')
           units_name = 'K'
        case ('gph')
           units_name = 'm'
        case ('rhi')
           units_name = '%rhi'
        case default
           units_name = 'vmr'
        end select
                                                                                
      do field=1, NumOfDataFields
        field_name = trim(dzm%name)//'-'//dataTitles(field)
        status = he5_zawrlattr(zaId, trim(dataTitles(field)), 'Missing Value', &
           & HE5T_NATIVE_FLOAT, 1, UNDEFINED_VALUE)
        status = he5_zawrlattr(zaId, trim(dataTitles(field)), 'Title', &
           & MLS_CHARTYPE, len_trim(field_name), field_name)
        if ((trim(dataTitles(field)) == 'L3dzDataCount') .or. &
           & (trim(dataTitles(field)) == 'PerMisPoints')) then
           status = he5_zawrlattr(zaId, trim(dataTitles(field)), 'Unit', &
              & MLS_CHARTYPE, 7, 'NoUnits')
	else 
           status = he5_zawrlattr(zaId, trim(dataTitles(field)), 'Unit', &
              & MLS_CHARTYPE, len_trim(units_name), units_name )
        endif 
        status = he5_zawrlattr(zaId, trim(dataTitles(field)), &
           & 'UniqueFieldDefinition', MLS_CHARTYPE, len_trim(MLS_FIELD), &
             & MLS_FIELD)
      enddo
                                                                                
      status = he5_zadetach(zaID)
      IF (status == -1) THEN
         msr = 'Failed to deattach to zonal ' // trim(dzm%name)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
                                                                                
      status = he5_zaclose(zafID)
      IF (status == -1) THEN  
	  msr = 'Failed to close file ' // trim(file) // ' after writing.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
                                                                                
!----------------------------------------------
  END SUBROUTINE WriteAttributes_l3dz
!----------------------------------------------
                                                                                
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
         CALL ReadL3DZData_HE5(swfid, swathName, dz)
      ELSE
         CALL MLSMessage(MLSMSG_Error, ModuleName, 'no hdfversion')
      ENDIF

!------------------------------------------------
    END SUBROUTINE ReadL3DZData
!------------------------------------------------

!------------------------------------------------
    SUBROUTINE ReadL3DZData_HE5(swfid, zonalName, dz)
!------------------------------------------------

! Brief description of subroutine
! This subroutine reads a data/diagnostic swath pair from a file into the
! L3DZData structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: zonalName

      INTEGER, INTENT(IN) :: swfid

      TYPE( L3DZData_T ), INTENT(OUT) :: dz

! Parameters

      CHARACTER (LEN=*),PARAMETER :: RL3DZ_ERR = 'Failed to read L3DZ field:  '

! Functions

      INTEGER, EXTERNAL :: & 
           & pgs_td_utcToTAI, &
           & he5_zaattach, he5_zadetach, he5_zadiminfo, he5_zaread

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=GridNameLen) :: dgName

      INTEGER :: start(2), stride(2), edge(2)
      INTEGER :: err, nlat, nlev, size, swid, status

      REAL, ALLOCATABLE :: r2(:,:)
      REAL, ALLOCATABLE :: rl(:), rp(:)

! Attach to the data swath for reading

      swid = he5_zaattach(swfid, zonalName)

      IF (status == -1) THEN
         msr = 'Failed to attach to zonal interface for reading ' // &
              & trim(zonalName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get dimension information

      size = he5_zadiminfo(swid, DIM_NAME2)

      IF (size == -1) THEN
         msr = SZ_ERR // DIM_NAME2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      nlev = size

      size = he5_zadiminfo(swid, DIML_NAME)

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

      dz%name = trim(zonalName)

      start = 0
      stride = 1
      edge(1) = DATE_LEN
      edge(2) = dz%nLats

      status = he5_zaread(swid, GEO_FIELD11, start(1), stride(1), edge(1), & 
           & dz%date)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD11
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      edge(1) = dz%nLevels
      status = he5_zaread(swid, GEO_FIELD9, start(1), stride(1), edge(1), & 
           & rp)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%pressure = DBLE(rp)

      status = he5_zaread(swid, GEO_FIELD1, start(2), stride(2), edge(2), & 
           & rl)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latitude = DBLE(rl)

      edge(1) = MIN_MAX
      status = he5_zaread(swid, GEO_FIELD4, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%localSolarTime = DBLE(r2)

      status = he5_zaread(swid, GEO_FIELD5, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // GEO_FIELD5
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%SolarZenithAngle = DBLE(r2)

      ! Read the data fields

      edge(1) = dz%nLevels
      status = he5_zaread(swid, DATA_FIELDV, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_FIELDV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%l3dzValue = DBLE(r2)

      status = he5_zaread(swid, DATA_STDDEV, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DATA_STDDEV
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latRss = DBLE(r2)

      !  After reading, detach from zonal interface
      
      status = he5_zadetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from zonal interface after reading ' // &
              & trim(zonalName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Attach to the diagnostic swath for reading

      dgName = TRIM(zonalName) // 'Diagnostics'
      swid = he5_zaattach(swfid, dgName)
      IF (status == -1) THEN
         msr = 'Failed to attach to zonal interface for reading ' // & 
              & trim(dgName)
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ! Read the diagnostic fields

      status = he5_zaread(swid, DG_FIELD1, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%latRss = DBLE(r2)

      status = he5_zaread(swid, DG_FIELD2, start, stride, edge, r2)
      IF (status == -1) THEN
         msr = RL3DZ_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      dz%perMisPoints = DBLE(r2)

      !  After reading, detach from swath interface

      status = he5_zadetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from zonal interface after reading ' // & 
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
   USE HDFEOS5, ONLY: HE5F_ACC_RDWR, HE5_GDCLOSE, HE5_GDOPEN 
   USE MLSPCF3
   USE dates_module, only: utc_to_yyyymmdd
   USE PCFHdr, ONLY: WritePCF2Hdr, WriteInputPointer, GlobalAttributes, &
	& he5_writeglobalattr
   USE mon_Open, ONLY: PCFMData_T 
   USE PCFModule, ONLY: SearchPCFDates, ExpandFileTemplate 
   USE SDPToolkit, only: WARNIFCANTPGSMETREMOVE, PGSD_MET_GROUP_NAME_L, &
        & PGSD_MET_NUM_OF_GROUPS, PGSMET_E_MAND_NOT_SET, PGS_S_SUCCESS, &
	& max_orbits

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
           & pgs_pc_getUniversalRef, he5_zainqza

      ! Variables

      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER (LEN=GridNameLen):: zonalName
      CHARACTER (LEN=FileNameLen):: sval
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=45) :: attrName, lvid
      CHARACTER (LEN=3) :: cNum
      CHARACTER (LEN=8) :: rangeDate
      ! CHARACTER (LEN=2) :: fileType
      integer :: fileType

      REAL(r8) :: dval

      INTEGER :: dg, hdfReturn, i, j, indx, len, numZonals, pNum, &
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

! MeasuredParameterContainer -- find the number of zonal in the file

         numZonals = 0
         numZonals = he5_zainqza(trim(files%name(i)), list, len)

         IF (numZonals .LE. 0) THEN
            msr = 'No Zonal found in file ' // TRIM( files%name(i) ) // &
                 & ' while attempting to write its metadata.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         pNum = 0

! For each zonal in the file

         DO j = 1, numZonals

! Extract its name

            indx = INDEX(list, ',')
            IF (indx /= 0) THEN
               zonalName = list(:indx-1)
               list = list(indx+1:)
            ELSE
               zonalName = list
            ENDIF

! If this is a diagnostic zonal, skip to the next one

            IF ( INDEX(zonalName, 'Diagnostics') /= 0 ) CYCLE

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
                 & zonalName)
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
         IF (maxval(GlobalAttributes%OrbNumDays(:,i)) == -1) then
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, 99999)
         ELSE
           result = pgs_met_setAttr_i(groups(INVENTORYMETADATA), &
                & attrName, maxval(GlobalAttributes%OrbNumDays(:,i)))
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

        ! Set global attributes of Granule (year,month,day) for each day 
        rangeDate = files%date(i)
        GlobalAttributes%StartUTC = rangeDate(1:4)//'-'//rangeDate(6:8) &
          & // 'T00:00:00.000000Z'
        GlobalAttributes%EndUTC = rangeDate(1:4)//'-'//rangeDate(6:8) &
          & // 'T23:59:59.999999Z'
        call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
          & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
          & GlobalAttributes%GranuleDay)

        ! Write global attributes
        if ( HDFVersion == HDFVERSION_5 ) then
            sdid = he5_gdopen (files%name(i), HE5F_ACC_RDWR)
            call he5_writeglobalattr(sdid,i)
            result = he5_gdclose (sdid)
        endif

         ! Annotate the file with the PCF
        CALL WritePCF2Hdr(files%name(i), anText, hdfVersion=hdfVersion, & 
              & fileType=fileType)

      ENDDO

      result = pgs_met_remove()
      
!      if (result /= PGS_S_SUCCESS .and. WARNIFCANTPGSMETREMOVE) THEN 
!        write(msr, *) result
!        CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!             & "Calling pgs_met_remove() failed with value " // trim(msr) )
!      endif          

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
         msr = MLSMSG_Allocate // 'Local solar time pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ALLOCATE(l3dz%SolarZenithAngle(MIN_MAX,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // 'Solar zenith angle  pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Data fields

      ALLOCATE(l3dz%dataCount(l3dz%nLevels,l3dz%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' data count pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

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
      
      IF ( ASSOCIATED(l3dz%SolarZenithAngle) ) THEN
         DEALLOCATE (l3dz%SolarZenithAngle, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate //' l3dz solar zenith angle pointer.'
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
      
      IF ( ASSOCIATED(l3dz%dataCount) ) THEN
         DEALLOCATE (l3dz%dataCount, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // ' dataCount pointer.'
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
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
 END MODULE L3DZData
 !==================

! $Log$
! Revision 1.23  2006/05/03 14:42:30  cvuu
! Remove subroutine OutputDZDiag, move datasets from Swath group to Zonal Means group
!
! Revision 1.22  2006/02/28 20:36:33  cvuu
! V2.00 commit
!
! Revision 1.21  2005/09/22 23:41:14  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 1.20  2005/06/23 19:17:57  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.19  2005/01/27 00:35:06  pwagner
! ReprocessingActual field dropped from product metadata
!
! Revision 1.18  2004/12/16 15:01:01  cvuu
! v1.5: Change value of ReprocessingActual to unknown in metadata file
!
! Revision 1.17  2004/12/13 17:41:26  cvuu
! remove writing QA flags to meta file, use the ones in MCF v1.5
!
! Revision 1.16  2004/08/04 23:21:23  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 1.15  2004/05/04 15:57:20  cvuu
! Fixed bug
!
! Revision 1.14  2004/01/08 21:21:36  cvuu
! version 1.4 commit
!
! Revision 1.13  2003/09/15 18:27:23  cvuu
! Output OrbitNumber and OrbitPeriod in the global attribute
!
! Revision 1.12  2003/08/11 23:26:37  cvuu
! brought closer to James Johnson want to
!
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
