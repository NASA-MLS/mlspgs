
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3SPData
!===============================================================================

   USE Hdf
   USE L2GPData, ONLY: L2GPNameLen, DIM_NAME2, DIM_NAME3, GEO_FIELD1, &
                       GEO_FIELD9, GEO_FIELD10, HDFE_NOMERGE
   USE L3CF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF
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

! Definition -- L3SPData_T
! Subroutines -- AllocateL3SP
!                OutputL3SP
!                DeallocateL3SP
!                DestroyL3SPDatabase

! Remarks:  This module contains the definition of the L3SPData type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DIMW_NAME = 'nWaveNum'
   CHARACTER (LEN=*), PARAMETER :: DIMWD_NAME = 'nLevels,nLats,nWaveNum'
   CHARACTER (LEN=*), PARAMETER :: DIMFD_NAME = 'nLevels,nLats,nFreqs'
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

     REAL(r8), DIMENSION(:,:,:), POINTER :: waveNumber
	! dimensioned (nLevels, nLats, nWaveNum)

     REAL(r8), DIMENSION(:,:,:), POINTER :: frequency
	! dimensioned (nLevels, nLats, nFreqs)

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
! This subroutine allocates the internal field pointers of the L3SPData_T
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

      ALLOCATE(l3sp%waveNumber(l3sp%nLevels,l3sp%nLats,l3sp%nWaveNum), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3SPData_T waveNumber pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3sp%frequency(l3sp%nLevels,l3sp%nLats,l3sp%nFreqs), STAT=err)
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

!--------------------------------------------
   SUBROUTINE OutputL3SP (cfProd, anText, sp)
!--------------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the swaths in an l3sp file.

! Arguments

      TYPE( L3CFProd_T ), INTENT(IN) :: cfProd

      TYPE( L3SPData_T ), INTENT(IN) :: sp(:)

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefgfld
      INTEGER, EXTERNAL :: swdefdim, swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=FileNameLen) :: spFile, type
      CHARACTER (LEN=480) :: msr

      INTEGER :: i, match, numSwaths, status, swfID, swID
      INTEGER :: start(4), stride(4), edge(4)

! Expand the template given in the CF

      CALL ExpandFileTemplate(cfProd%fileTemplate, type, 'L3SP')

! Find a PCF entry for this file type

      CALL FindFileType(type, mlspcf_l3sp_start, mlspcf_l3sp_end, match, &
                        spFile)

! Open the l3sp file for the creation of swaths

      swfID = swopen(spFile, DFACC_CREATE)
      IF (swfID == -1) THEN
         msr = MLSMSG_Fileopen // spFile
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! For each swath found,

      numSwaths = SIZE(sp)

      DO i = 1, numSwaths

! Create a swath of the appropriate name

         swID = swcreate( swfID, sp(i)%name )
         IF (swID == -1) THEN
            msr = 'Failed to create swath ' // sp(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the swath dimensions

         status = swdefdim(swID, DIM_NAME2, sp(i)%nLevels)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swID, DIML_NAME, sp(i)%nLats)
         IF (status == -1) THEN
            msr = DIM_ERR // DIML_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swID, DIMW_NAME, sp(i)%nWaveNum)
         IF (status == -1) THEN
            msr = DIM_ERR // DIMW_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swID, DIM_NAME3, sp(i)%nFreqs)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the vertical geolocation field, using the above dimensions

         status = swdefgfld(swID, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the horizontal geolocation field

         status = swdefgfld(swID, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the waveNumber & frequency fields

         status = swdefgfld(swID, GEO_FIELDWN, DIMWD_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELDWN
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swID, GEO_FIELD10, DIMFD_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD10
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the data fields using above dimensions

         status = swdefdfld(swID, DATA_FIELDRV, DIMSP_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDRV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDRP, DIMSP_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDRP 
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDIV, DIMSP_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDIV 
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swID, DATA_FIELDIP, DIMSP_NAME, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELDIP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from the swath interface

         status = swdetach(swID)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                       &to detach from swath interface after L3SP definition.')

! Re-attach to the swath for writing

         swID = swattach(swfID, sp(i)%name)
         IF (swID == -1) THEN
            msr = 'Failed to re-attach to swath ' // TRIM(sp(i)%name) &
                   // ' for writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the data

         start = 0
         stride = 1
         edge(1) = sp(i)%nLevels
         edge(2) = sp(i)%nLats
         edge(3) = sp(i)%nWaveNum
         edge(4) = sp(i)%nFreqs

! Write the vertical geolocation data

         status = swwrfld( swID, GEO_FIELD9, start(1), stride(1), edge(1), &
                           REAL(sp(i)%pressure) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the horizontal geolocation data

         status = swwrfld( swID, GEO_FIELD1, start(2), stride(2), edge(2), &
                           REAL(sp(i)%latitude) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the waveNumber data

         status = swwrfld( swID, GEO_FIELDWN, start(1:3), stride(1:3), &
                           edge(1:3), REAL(sp(i)%waveNumber) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELDWN
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the frequency data

         edge(3) = sp(i)%nFreqs

         status = swwrfld( swID, GEO_FIELD10, start(1:3), stride(1:3), &
                           edge(1:3), REAL(sp(i)%frequency) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD10
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the real part of the data

         edge(3) = sp(i)%nWaveNum

         status = swwrfld( swID, DATA_FIELDRV, start, stride, edge, &
                           REAL(sp(i)%l3spRelValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDRV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, DATA_FIELDRP, start, stride, edge, &
                           REAL(sp(i)%l3spRelPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDRP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write the imaginary part of the data

         status = swwrfld( swID, DATA_FIELDIV, start, stride, edge, &
                        REAL(sp(i)%l3spImgValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDIV
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swID, DATA_FIELDIP, start, stride, edge, &
                           REAL(sp(i)%l3spImgPrecision) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELDIP
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! After writing, detach from swath interface

         status = swdetach(swID)
         IF (status == -1) THEN
            CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to detach from &
                                              &swath interface after writing.')
         ELSE
            msr = 'Swath ' // TRIM(sp(i)%name) // ' successfully written to &
                  &file ' // spFile
            CALL MLSMessage(MLSMSG_Info, ModuleName, msr)
         ENDIF

      ENDDO

! Close the file

      status = swclose(swfID)
      IF (status == -1) THEN
         msr = 'Failed to close file ' // spFile // ' after writing.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Annotate the file with the PCF

      CALL WritePCF2Hdr(spFile, anText)

!---------------------------
   END SUBROUTINE OutputL3SP
!---------------------------

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
         ENDDO

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
! Revision 1.4  2001/02/09 19:17:05  nakamura
! Changed dimensions on waveNumber & frequency.
!
! Revision 1.3  2001/01/16 17:44:38  nakamura
! Added annotation.
!
! Revision 1.2  2000/12/29 20:51:11  nakamura
! Added subroutine OutputL3SP.
!
! Revision 1.1  2000/11/29 21:44:16  nakamura
! Module for the L3SP data type.
!
