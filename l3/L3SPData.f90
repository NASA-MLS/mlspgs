
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3SPData
!===============================================================================

   USE L2GPData, ONLY: L2GPNameLen
   USE MLSCommon
   USE MLSMessageModule
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L3SPData_T
! Subroutines -- AllocateL3SP
!                DeallocateL3SP
!                DestroyL3SPDatabase

! Remarks:  This module contains the definition of the L3SPData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DIML_NAME = 'nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMW_NAME = 'nWaveNum'
   CHARACTER (LEN=*), PARAMETER :: DIMSP_NAME = 'nLevels,nLats,nWaveNum,nFreqs'

   CHARACTER (LEN=*), PARAMETER :: GEO_FIELDWN = 'waveNumber'

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDRV = 'L3spRelValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDRP = 'L3spRelPrecision'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDIV = 'L3spImgValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELDIP = 'L3spImgPrecision'

! This data type is used to store the L3 fourier spectra.

  TYPE L3SPData_T

     CHARACTER (LEN=L2GPNameLen) :: name	! name for the output quantity

     INTEGER :: nLevels				! Total number of surfaces
     INTEGER :: nLats				! Total number of latitudes
     INTEGER :: nWaveNum			! Total number of waveNumbers
     INTEGER :: nFreqs				! Total number of frequencies

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude        ! dimensioned (nLats)

     REAL(r8) :: startTime	! start time of L2 data used in the analysis

     REAL(r8) :: endTime	! end time of L2 data used in the analysis

     ! Now we store the waveNumber and frequency fields:

     REAL(r8), DIMENSION(:), POINTER :: waveNumber	! dimensioned (nWaveNum)
     REAL(r8), DIMENSION(:), POINTER :: frequency	! dimensioned (nFreqs)

     ! Now the spectra fields:

     REAL(r8), DIMENSION(:,:,:,:), POINTER :: l3spRelValue     ! Real part
     REAL(r8), DIMENSION(:,:,:,:), POINTER :: l3spRelPrecision
						! Real part precision
     REAL(r8), DIMENSION(:,:,:,:), POINTER :: l3spImgValue     ! Imaginary part
     REAL(r8), DIMENSION(:,:,:,:), POINTER :: l3spImgPrecision
						! Imaginary part precision
     ! dimensioned as (nLevels, nLats, nWaveNum, nFreqs)

   END TYPE L3SPData_T

CONTAINS

!--------------------------------------------------------
   SUBROUTINE AllocateL3SP (nlev, nlat, nwv, nfreq, l3sp)
!--------------------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3DMData_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat, nwv, nfreq

      TYPE( L3SPData_T ), INTENT(INOUT) :: l3sp

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3sp%nLevels = nlev
      l3sp%nLats = nlat
      l3sp%nWaveNum = nwv
      l3sp%nFreqs = nfreq

! Allocate the vertical geolocation field

      ALLOCATE(l3sp%pressure(l3sp%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      ALLOCATE(l3sp%latitude(l3sp%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Wavenumber & frequency

      ALLOCATE(l3sp%waveNumber(l3sp%nWaveNum), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T waveNumber pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3sp%frequency(l3sp%nFreqs), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T frequency pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Spectra fields

      ALLOCATE(l3sp%l3spRelValue(l3sp%nLevels,l3sp%nLats,l3sp%nWaveNum,&
              &l3sp%nFreqs), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T RelValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3sp%l3spRelPrecision(l3sp%nLevels,l3sp%nLats,l3sp%nWaveNum,&
              &l3sp%nFreqs), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T RelPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3sp%l3spImgValue(l3sp%nLevels,l3sp%nLats,l3sp%nWaveNum,&
              &l3sp%nFreqs), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T ImgValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3sp%l3spImgPrecision(l3sp%nLevels,l3sp%nLats,l3sp%nWaveNum,&
              &l3sp%nFreqs), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T ImgPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE AllocateL3SP
!-----------------------------

!----------------------------------
   SUBROUTINE DeallocateL3SP (l3sp)
!----------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L2GP_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3SPData_T ), INTENT(INOUT) :: l3sp

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Vertical geolocation field

      IF ( ASSOCIATED(l3sp%pressure) ) DEALLOCATE (l3sp%pressure, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3sp pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      IF ( ASSOCIATED(l3sp%latitude) ) DEALLOCATE (l3sp%latitude, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3sp latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Wavenumber & frequency fields

      IF ( ASSOCIATED(l3sp%waveNumber) ) DEALLOCATE (l3sp%waveNumber, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3sp waveNumber pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3sp%frequency) ) DEALLOCATE (l3sp%frequency, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3sp frequency pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Spectra fields

      IF ( ASSOCIATED(l3sp%l3spRelValue) ) &
         DEALLOCATE (l3sp%l3spRelValue,  STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3spRelValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3sp%l3spRelPrecision) ) &
         DEALLOCATE (l3sp%l3spRelPrecision, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3spRelPrecision'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3sp%l3spImgValue) ) &
         DEALLOCATE (l3sp%l3spImgValue, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3spImgValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF ( ASSOCIATED(l3sp%l3spImgPrecision) ) &
         DEALLOCATE (l3sp%l3spImgPrecision, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  l3spImgPrecision'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL3SP
!-------------------------------

!-----------------------------------------
   SUBROUTINE DestroyL3SPDatabase (l3spdb)
!-----------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal structures of an l3sp database, and
! then database itself


! Arguments

      TYPE (L3SPData_T), DIMENSION(:), POINTER :: l3spdb

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i

! Check the status of the input pointer

      IF ( ASSOCIATED(l3spdb) ) THEN

! If it's associated, then deallocate the internal structures

         DO i = 1, SIZE(l3spdb)
            CALL DeallocateL3SP( l3spdb(i) )
         END DO

! Deallocate the database itself

         DEALLOCATE (l3spdb, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3sp database'
            CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

!------------------------------------
   END SUBROUTINE DestroyL3SPDatabase
!------------------------------------

!==================
END MODULE L3SPData
!==================

! $Log$
!
