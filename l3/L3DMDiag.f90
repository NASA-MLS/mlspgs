
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3DMDiag
!===============================================================================

   USE Hdf
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
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

! Definition -- L3DMDiag_T
! Subroutines -- OutputDiags
!                ReadL3DMDiag
!                AllocateL3DMDiag
!                DeallocateL3DMDiag
!                DestroyL3DMDiagDB

! Remarks:  This module contains the definition of the L3DMDiag type, as well
!           as any routines pertaining to it.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DIMN_NAME = 'N'

   CHARACTER (LEN=*), PARAMETER :: DG_FIELD = 'GRss'
   CHARACTER (LEN=*), PARAMETER :: MD_FIELD = 'MaxDiff'

! This data type is used to store the l3 daily map diagnostics.

   TYPE L3DMDiag_T

     CHARACTER (LEN=GridNameLen) :: name        ! name for the output quantity

     INTEGER :: N		! number for "largest differences" diagnostics
     INTEGER :: nLevels		! Total number of surfaces
     INTEGER :: nLats		! Total number of latitudes

     ! Now we store the geolocation fields.  First, the vertical one:

     REAL(r8), DIMENSION(:), POINTER :: pressure	! dimensioned (nLevels)

     ! Now the horizontal geolocation information and time:

     REAL(r8), DIMENSION(:), POINTER :: latitude	! dimensioned (nLats)

     REAL(r8) :: time					! Synoptic time

     ! Global Root-Sum_Square, dimensioned (nLevels)

     REAL(r8), DIMENSION(:), POINTER :: gRss

     ! Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

     REAL(r8), DIMENSION(:,:), POINTER :: latRss

     ! Maximum difference, dimensioned (N, nLevels)

     REAL(r8), DIMENSION(:,:), POINTER :: maxDiff

     ! Missing points (percentage), dimensioned (nLevels)

     INTEGER, DIMENSION(:), POINTER :: perMisPoints

   END TYPE L3DMDiag_T

CONTAINS

!--------------------------------------
   SUBROUTINE OutputDiags(type, dgData)
!--------------------------------------

! Brief description of subroutine
! This subroutine creates and writes to the diagnostic portion of the l3dm files.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE (L3DMDiag_T), INTENT(IN) :: dgData(:)

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DIMNL_NAME = 'N,nLevels'

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swcreate, swdefdfld, swdefdim, swdefgfld
      INTEGER, EXTERNAL :: swdetach, swopen, swwrfld

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=FileNameLen) :: physicalFilename
      CHARACTER (LEN=480) :: msr

      INTEGER :: swfID, swId, i, match, status
      INTEGER :: start(2), stride(2), edge(2)

! For each day in the l3dg database with corresponding map data,

      DO i = 1, SIZE(dgData)

! Find the output file for the level/species/day in the PCF (existence already checked
! in DM grid output)

         CALL FindFileDay(type, dgData(i)%time, mlspcf_l3dm_start, &
                          mlspcf_l3dm_end, match, physicalFilename, date)

! Re-open the file for the creation of diagnositc swaths

         swfID = swopen(physicalFilename, DFACC_RDWR)
         IF (swfID == -1) THEN
            msr = MLSMSG_Fileopen // physicalFilename
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Create the swath

         swId = swcreate(swfID, dgData(i)%name)
         IF (swId == -1) THEN
            msr = 'Failed to create swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the dimensions

         status = swdefdim(swId, DIMN_NAME, dgData(i)%N)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMN_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swId, DIML_NAME, dgData(i)%nLats)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMY_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swId, DIM_NAME2, dgData(i)%nLevels)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMZ_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdim(swId, DIMT_NAME, 1)
         IF (status /= 0) THEN
            msr = DIM_ERR // DIMT_NAME
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the "geolocation" fields using the above dimensions

         status = swdefgfld(swId, GEO_FIELD3, DIMT_NAME, DFNT_FLOAT64, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swId, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefgfld(swId, GEO_FIELD1, DIML_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = GEO_ERR // GEO_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Define the "data" fields

         status = swdefdfld(swId, DG_FIELD, DIM_NAME2, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // DG_FIELD
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swId, DG_FIELD1, DIMLL_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // DG_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swId, MD_FIELD, DIMNL_NAME, DFNT_FLOAT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // MD_FIELD
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swdefdfld(swId, DG_FIELD2, DIM_NAME2, DFNT_INT32, &
                           HDFE_NOMERGE)
         IF (status /= 0) THEN
            msr = DAT_ERR // DG_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from the swath interface after definition

         status = swdetach(swId)
         IF (status /= 0) THEN
            msr = 'Failed to detach from swath ' // TRIM( dgData(i)%name ) // &
                   ' after definition.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Re-attach to the swath for writing

         swId = swattach(swfID, dgData(i)%name)
         IF (swId == -1) THEN
            msr = 'Failed to re-attach to swath ' // TRIM( dgData(i)%name ) // &
                   ' for wiriting.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write to fields

         start = 0
         stride = 1
         edge(1) = dgData(i)%nLevels
         edge(2) = dgData(i)%nLats

! Geolocation

         status = swwrfld(swId, GEO_FIELD3, start(1), stride(1), stride(1), &
                          dgData(i)%time)
         IF (status /=0) THEN
            msr = WR_ERR //  GEO_FIELD3 // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swId, GEO_FIELD9, start(1), stride(1), edge(1), &
                           REAL(dgData(i)%pressure) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD9 // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swId, GEO_FIELD1, start(2), stride(2), edge(2), &
                           REAL(dgData(i)%latitude) )
         IF (status /= 0) THEN
            msr = WR_ERR //  GEO_FIELD1 // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! One-dimensional data fields

         status = swwrfld( swId, DG_FIELD, start(1), stride(1), edge(1), &
                           REAL(dgData(i)%gRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swwrfld( swId, DG_FIELD2, start(1), stride(1), edge(1), &
                           dgData(i)%perMisPoints )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD2 // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Two-dimensional data fields

         status = swwrfld( swId, DG_FIELD1, start, stride, edge, &
                           REAL(dgData(i)%latRss) )
         IF (status /= 0) THEN
            msr = WR_ERR //  DG_FIELD1 // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = dgData(i)%N
         edge(2) = dgData(i)%nLevels
         status = swwrfld( swId, MD_FIELD, start, stride, edge, &
                           REAL(dgData(i)%maxDiff) )
         IF (status /= 0) THEN
            msr = WR_ERR //  MD_FIELD // ' to swath ' // dgData(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Detach from the swath after writing

         status = swdetach(swId)
         IF (status /= 0) THEN
            msr = GD_ERR // TRIM( dgData(i)%name ) // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Close the file after writing

         status = swclose(swfID)
         IF (status /= 0) THEN
            msr = 'Failed to close file ' // TRIM(physicalFilename) // &
                  ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         msr = 'Swath ' // TRIM(dgData(i)%name) // ' successfully written to &
               &file ' // physicalFilename
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

      ENDDO

!----------------------------
   END SUBROUTINE OutputDiags
!----------------------------

!-----------------------------------------------
   SUBROUTINE ReadL3DMDiag (swfid, dgName, l3dg)
!-----------------------------------------------

! Brief description of subroutine
! This  subroutine reads a swath from a file into the L3DMDiag structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: dgName

      INTEGER, INTENT(IN) :: swfid  

      TYPE( L3DMDiag_T ), INTENT(OUT) :: l3dg

! Parameters

      CHARACTER (LEN=*), PARAMETER :: RDSW_ERR = 'Failed to read swath field '

! Functions

      INTEGER :: swattach, swdetach, swdiminfo, swrdfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, N, nlev, nlat, size, status, swId
      INTEGER :: start(2), stride(2), edge(2)

      REAL, ALLOCATABLE :: rp(:), rlat(:)
      REAL, ALLOCATABLE :: rll(:,:), rnl(:,:)

! Attach to the swath

      swId = swattach(swfid, dgName)
      IF (swId == -1) THEN
         msr = 'Failed to attach to swath ' // dgName
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get dimension names, sizes

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
      
      size = swdiminfo(swid, DIMN_NAME)
      IF (size == -1) THEN
         msr = SZ_ERR // DIMN_NAME
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      N = size

! Allocate the local REAL variables & the output structure pointers

      ALLOCATE(rp(nlev), rlat(nlat), rll(nlev,nlat), rnl(N,nlev), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' local REAL variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL AllocateL3DMDiag(nlev, nlat, N, l3dg)

! Read data into the output structure 

      l3dg%name = dgName

      start = 0
      stride = 1
      edge(1) = l3dg%nLevels
      edge(2) = l3dg%nLats

      status = swrdfld(swid, GEO_FIELD3, start(1), stride(1), stride(1), &
                       l3dg%time)
      IF (status == -1) THEN
         msr = RDSW_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swrdfld(swid, GEO_FIELD9, start(1), stride(1), edge(1), rp)
      IF (status == -1) THEN
         msr = RDSW_ERR // GEO_FIELD9
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dg%pressure = DBLE(rp)

      status = swrdfld(swid, GEO_FIELD1, start(2), stride(2), edge(2), rlat)
      IF (status == -1) THEN
         msr = RDSW_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dg%latitude = DBLE(rlat)

      status = swrdfld(swid, DG_FIELD, start(1), stride(1), edge(1), rp)
      IF (status == -1) THEN
         msr = RDSW_ERR // DG_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dg%gRss = DBLE(rp)

      status = swrdfld(swid, DG_FIELD2, start(1), stride(1), edge(1), &
                       l3dg%perMisPoints)
      IF (status == -1) THEN
         msr = RDSW_ERR // DG_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swrdfld(swid, DG_FIELD1, start, stride, edge, rll)
      IF (status == -1) THEN
         msr = RDSW_ERR // DG_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dg%latRss = DBLE(rll)

      edge(1) = l3dg%N
      edge(2) = l3dg%nLevels

      status = swrdfld(swid, MD_FIELD, start, stride, edge, rnl)
      IF (status == -1) THEN
         msr = RDSW_ERR // MD_FIELD
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l3dg%maxDiff = DBLE(rnl)

! Detach from and close the grid interface

      status = swdetach(swid)
      IF (status == -1) THEN
         msr = 'Failed to detach from swath ' // TRIM(dgName) // ' after reading.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Deallocate the local variables

      DEALLOCATE (rp, rlat, rll, rnl, STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_DeAllocate // '  local REAL variables.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE ReadL3DMDiag
!-----------------------------

!---------------------------------------------------
   SUBROUTINE AllocateL3DMDiag (nlev, nlat, N, l3dg)
!---------------------------------------------------

! Brief description of subroutine
! This subroutine allocates the internal field pointers of the L3DMDiag_T
! derived type.

! Arguments

      INTEGER, INTENT(IN) :: nlev, nlat, N

      TYPE( L3DMDiag_T ), INTENT(INOUT) :: l3dg

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Store the sizes of the dimensions

      l3dg%N = N
      l3dg%nLevels = nlev
      l3dg%nLats = nlat

! Vertical geolocation field

      ALLOCATE(l3dg%pressure(l3dg%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T pressure pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Horizontal geolocation field

      ALLOCATE(l3dg%latitude(l3dg%nLats), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Quantities

      ALLOCATE(l3dg%gRss(l3dg%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T gRss pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dg%latRss(l3dg%nLevels,l3dg%nLats),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T latRss pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dg%maxDiff(l3dg%N,l3dg%nLevels),STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T maxDiff pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE(l3dg%perMisPoints(l3dg%nLevels), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' L3DMDiag_T perMisPoints pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!---------------------------------
   END SUBROUTINE AllocateL3DMDiag
!---------------------------------

!--------------------------------------
   SUBROUTINE DeallocateL3DMDiag (l3dg)
!--------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L3DMDiag_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L3DMDiag_T ), INTENT(INOUT) :: l3dg

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

! Horizontal geolocation fields

      IF ( ASSOCIATED(l3dg%latitude) ) THEN
         DEALLOCATE (l3dg%latitude, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dg latitude pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Vertical geolocation field

      IF ( ASSOCIATED(l3dg%pressure) ) THEN
         DEALLOCATE (l3dg%pressure, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dg pressure pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Data fields

      IF ( ASSOCIATED(l3dg%gRss) ) THEN
         DEALLOCATE (l3dg%gRss, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  gRss pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dg%latRss) ) THEN
         DEALLOCATE (l3dg%latRss, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  latRss pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dg%maxDiff) ) THEN
         DEALLOCATE (l3dg%maxDiff, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  maxDiff pointer'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF ( ASSOCIATED(l3dg%perMisPoints) ) THEN
         DEALLOCATE (l3dg%perMisPoints, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  perMisPoints'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

!-----------------------------------
   END SUBROUTINE DeallocateL3DMDiag
!-----------------------------------

!---------------------------------------
   SUBROUTINE DestroyL3DMDiagDB (l3dgdb)
!---------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal structures of an l3dg database, and
! then database itself


! Arguments

      TYPE (L3DMDiag_T), DIMENSION(:), POINTER :: l3dgdb

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i

! Check the status of the input pointer

      IF ( ASSOCIATED(l3dgdb) ) THEN

! If it's associated, then deallocate the internal structures

         DO i = 1, SIZE(l3dgdb)
            CALL DeallocateL3DMDiag( l3dgdb(i) )
         ENDDO

! Deallocate the database itself

         DEALLOCATE (l3dgdb, STAT=err)
         IF ( err /= 0 ) THEN
            msr = MLSMSG_DeAllocate // '  l3dm diagnostics database'
            CALL MLSMessage ( MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

!----------------------------------
   END SUBROUTINE DestroyL3DMDiagDB
!----------------------------------

!==================
END MODULE L3DMDiag
!==================

!# $Log$
!#
