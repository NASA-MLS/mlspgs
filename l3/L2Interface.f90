
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L2Interface
!===============================================================================

   USE Hdf
   USE L2GPData, ONLY: L2GPData_T
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE MLSStrings
   USE OpenInit
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

! Subroutines -- GetL2GPfromPCF
!                ReadL2GP
!                ReadL2GPProd
!                ResidualCreate
!                ResidualWrite
!                ResidualOutput
!                ReadL2DGData

! Remarks:  This module contains subroutines involving the interface between
!           Level 3 and L2.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD1 = 'L2gpValue'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD2 = 'L2gpPrecision'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD3 = 'Status'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD4 = 'Quality'
   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD5 = 'L3Residual'

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE GetL2GPfromPCF (mlspcf_start, mlspcf_end, template, startDOY, &
                              endDOY, numFiles, pcfNames, mis_numFiles, mis_Days)
!----------------------------------------------------------------------------

! Brief description of subroutine
! This routine searches the PCF for l2gp files for a desired species.  It
! returns an array of names and the number of the files found.

! Arguments

      CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY

      CHARACTER (LEN=*), INTENT(IN) :: template

      INTEGER, INTENT(IN) :: mlspcf_end, mlspcf_start

      INTEGER, INTENT(OUT) :: numFiles, mis_numFiles

      CHARACTER (LEN=FileNameLen) :: pcfNames(:)
      CHARACTER (LEN=DATE_LEN), INTENT(OUT) :: mis_Days(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=FileNameLen) :: physicalFilename, type

      INTEGER :: i, indx, returnStatus, version

! Expand the level/species template

      CALL ExpandFileTemplate(template, type, 'L2GP')

      numFiles = 0
      mis_numFiles = 0

! Loop through all the PCF numbers for L2GP files

      DO i = mlspcf_start, mlspcf_end

         version = 1
         returnStatus = Pgs_pc_getReference(i, version, physicalFilename)

! If no file name was returned, go on to the next PCF number

	 IF(returnStatus /= PGS_S_SUCCESS) THEN

            IF ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN

! Extract the date from the file name

            	indx = INDEX(physicalFilename, '.', .TRUE.)
            	date = physicalFilename(indx-8:indx-1)

! Check that the date is within the desired boundaries for reading

            	IF ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN

! Save the name in an array of files for the species

               		mis_numFiles = mis_numFiles + 1
               		mis_Days(mis_numFiles) = date

            	ENDIF

            ENDIF

	    CYCLE

         ENDIF

! Check that the returned file name is for the proper species

         IF ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN

! Extract the date from the file name

            indx = INDEX(physicalFilename, '.', .TRUE.)
            date = physicalFilename(indx-8:indx-1)

! Check that the date is within the desired boundaries for reading

            IF ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN

! Save the name in an array of files for the species

               numFiles = numFiles + 1
               pcfNames(numFiles) = physicalFilename

            ENDIF

         ENDIF

      ENDDO

!-------------------------------
   END SUBROUTINE GetL2GPfromPCF
!-------------------------------

!--------------------------------------------------------
   SUBROUTINE ReadL2GP(L2FileHandle, swathname, l2gpData)
!--------------------------------------------------------

! Brief description of subroutine
! This routine reads an L2GP file, returning a filled data structure.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: swathname

      INTEGER, INTENT(IN) :: L2FileHandle

      TYPE( L2GPData_T ), INTENT(OUT) :: l2gpData

! Parameters

      CHARACTER (LEN=*), PARAMETER :: MLSMSG_L2GPRead = 'Unable to read L2GP field:'

! Functions

      INTEGER, EXTERNAL :: swattach, swdetach, swdiminfo, swinqdims, swrdfld

! Variables

      CHARACTER (LEN=80) :: list
      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, freq, lev, nDims, numFreq, numLev, numProfs, size, swid
      INTEGER :: status
      INTEGER :: start(3), stride(3), edge(3), dims(3)

      REAL, ALLOCATABLE :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)

! Attach to the swath for reading

      l2gpData%name = swathname

      swid = swattach(L2FileHandle, l2gpData%name)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                     &attach to swath interface for reading.')

! Get dimension information

      lev = 0
      freq = 0

      nDims = swinqdims(swid, list, dims)
      IF (nDims == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                               &get dimension information.')
      IF ( INDEX(list,'nLevels') /= 0 ) lev = 1
      IF ( INDEX(list,'nFreqs') /= 0 ) freq = 1

      size = swdiminfo(swid, DIM_NAME1)
      IF (size == -1) THEN
         msr = SZ_ERR // DIM_NAME1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%nTimes = size
      numProfs = size

      IF (lev == 0) THEN
         l2gpData%nLevels = 0
         numLev = 1
      ELSE
         size = swdiminfo(swid, DIM_NAME2)
         IF (size == -1) THEN
            msr = SZ_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%nLevels = size
         numLev = size
      ENDIF

      IF (freq == 1) THEN
         size = swdiminfo(swid, DIM_NAME3)
         IF (size == -1) THEN
            msr = SZ_ERR // DIM_NAME3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%nFreqs = size
         numFreq = size
      ELSE
         l2gpData%nFreqs = 0
         numFreq = 1
      ENDIF

! Allocate output pointers for fields that are always present

      if (numProfs.gt.0) then 

      ALLOCATE (l2gpData%latitude(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  latitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%longitude(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  longitude pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%time(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  time pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%solarTime(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  solarTime pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%solarZenith(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  solarZenith pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%losAngle(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  losAngle pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%geodAngle(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  geodAngle pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%chunkNumber(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  chunkNumber pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%status(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  status pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%quality(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  quality pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (realProf(numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  realProf array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      endif

      if (numLev.gt.0) then 

      ALLOCATE (realSurf(numLev), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  realSurf array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      endif

      if (numFreq.gt.0) then 

      ALLOCATE (realFreq(numFreq), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  realFreq array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      endif

      if ((numFreq.gt.0).and.(numLev.gt.0).and.(numProfs.gt.0)) then 

      ALLOCATE (l2gpData%l2gpValue(numFreq,numLev,numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  l2gpValue pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (l2gpData%l2gpPrecision(numFreq,numLev,numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  l2gpPrecision pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      ALLOCATE (real3(numFreq,numLev,numProfs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  real3 array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      endif

! Read the horizontal geolocation fields

      start = 0
      stride = 1
      edge(1) = numFreq
      edge(2) = numLev
      edge(3) = numProfs

      if (numProfs .gt. 0) then

      status = swrdfld(swid, GEO_FIELD1, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%latitude = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD2, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%longitude = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD3, start(3), stride(3), edge(3), &
        l2gpData%time)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swrdfld(swid, GEO_FIELD4, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%solarTime = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD5, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD5
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%solarZenith = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD6, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD6
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%losAngle = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD7, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD7
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%geodAngle = DBLE(realProf)

      status = swrdfld(swid, GEO_FIELD8, start(3), stride(3), edge(3), &
                       l2gpData%chunkNumber)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // GEO_FIELD8
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      ENDIF

      endif

! Read the pressures vertical geolocation field, if it exists

      IF (lev /= 0) THEN

         ALLOCATE (l2gpData%pressures(l2gpData%nLevels), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // '  l2gp pressure pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         status = swrdfld(swid, GEO_FIELD9, start(2), stride(2), edge(2), realSurf)
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%pressures = DBLE(realSurf)
     ENDIF

! Read the frequency geolocation field, if it exists

      IF (freq == 1) THEN

         ALLOCATE (l2gpData%frequency(l2gpData%nFreqs), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // '  l2gp frequency pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         edge(1) = l2gpData%nFreqs

         status = swrdfld(swid, GEO_FIELD10, start(1), stride(1), edge(1), realFreq)
         IF (status == -1) THEN
             msr = MLSMSG_L2GPRead // GEO_FIELD10
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%frequency = DBLE(realFreq)

      ENDIF

! Read the data fields that may have 1-3 dimensions

      IF ( freq == 1) THEN

         status = swrdfld(swid, DATA_FIELD1, start, stride, edge, real3)
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpValue = DBLE(real3)

         status = swrdfld(swid, DATA_FIELD2, start, stride, edge, real3)
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpPrecision = DBLE(real3)

      ELSE IF ( lev == 1) THEN

         status = swrdfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
                           edge(2:3), real3(1,:,:) )
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpValue = DBLE(real3)

         status = swrdfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
                           edge(2:3), real3(1,:,:) )
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpPrecision = DBLE(real3)

      ELSE

         status = swrdfld( swid, DATA_FIELD1, start(3), stride(3), edge(3), &
                           real3(1,1,:) )
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD1
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpValue = DBLE(real3)

         status = swrdfld( swid, DATA_FIELD2, start(3), stride(3), edge(3), &
                           real3(1,1,:) )
         IF (status == -1) THEN
            msr = MLSMSG_L2GPRead // DATA_FIELD2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         l2gpData%l2gpPrecision = DBLE(real3)

      ENDIF

! Read the data fields that are 1-dimensional

      if (numProfs .gt. 0) then

      status = swrdfld(swid, DATA_FIELD3, start(3), stride(3), edge(3), &
                       l2gpData%status)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // DATA_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swrdfld(swid, DATA_FIELD4, start(3), stride(3), edge(3), realProf)
      IF (status == -1) THEN
         msr = MLSMSG_L2GPRead // DATA_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gpData%quality = DBLE(realProf)

      endif

! Deallocate local variables

      DEALLOCATE(realFreq, STAT=alloc_err)
      IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                'Failed deallocation of local realFreq array.')

      DEALLOCATE(realSurf, STAT=alloc_err)
      IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                'Failed deallocation of local realSurf array.')

      DEALLOCATE(realProf, STAT=alloc_err)
      IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                'Failed deallocation of local realProf array.')

      DEALLOCATE(real3, STAT=alloc_err)
      IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                'Failed deallocation of local real3 array.')

!  After reading, detach from swath interface

      status = swdetach(swid)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                 &detach from swath interface after reading.')

!-------------------------
   END SUBROUTINE ReadL2GP
!-------------------------

!----------------------------------------------------------------------------------
   SUBROUTINE ReadL2GPProd (product, template, startDOY, endDOY, numDays, mis_numDays, mis_Days, prodL2GP)
!----------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reads l2gp files for a quantity over a range of days.  It
! returns an array of L2GPData_T structures, with one structure per day.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: product, template

      CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY

      INTEGER, INTENT(OUT) :: numDays, mis_numDays

      TYPE( L2GPData_T ), POINTER :: prodL2GP(:)

      CHARACTER (LEN=DATE_LEN), INTENT(OUT) :: mis_Days(maxWindow)

! Parameters

! Functions

      INTEGER, EXTERNAL :: swclose, swopen

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: pcfNames(40)

      INTEGER :: err, i, l2id, status

! Search the PCF for L2GP files for this species

      CALL GetL2GPfromPCF(mlspcf_l2gp_start, mlspcf_l2gp_end, template, startDOY, &
                          endDOY, numDays, pcfNames, mis_numDays, mis_Days)

! Check that numDays is at least 1, so pointer allocation is possible

      IF (numDays == 0) THEN

! If not, issue a warning

         msr = 'No L2GP data found for ' // product
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         NULLIFY(prodL2GP)

! If so, open, read, & close the files for this quantity

      ELSE

         ALLOCATE( prodL2GP(numDays), STAT=err )
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' l2gp array for ' // product
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         DO i = 1, numDays

            l2id = swopen(pcfNames(i), DFACC_READ)
            IF (l2id == -1) THEN
               msr = MLSMSG_Fileopen // pcfNames(i)
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

! Read information from the L2GP file

            CALL ReadL2GP( l2id, product, prodL2GP(i) )

! Close the L2GP file

            status = swclose(l2id)
            IF (status == -1) THEN
               msr = 'Failed to close file ' // pcfNames(i)
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

         ENDDO


      ENDIF

!-----------------------------
   END SUBROUTINE ReadL2GPProd
!-----------------------------

!------------------------------------------------
   SUBROUTINE ResidualCreate (L3FileHandle, l2gp)
!------------------------------------------------

! Brief description of subroutine
! This subroutine creates the L3Residual diagnostic swath in L3 map files.

! Arguments

      INTEGER, INTENT(IN) :: L3FileHandle

      TYPE( L2GPData_T ), INTENT(IN) :: l2gp

! Parameters

! Functions

      INTEGER, EXTERNAL :: swcreate, swdefdfld, swdefdim, swdefgfld, swdetach

! Variables

      CHARACTER (LEN=480) :: msr
  
      INTEGER :: swid, status

! Create the swath within the file

      swid = swcreate(L3FileHandle, l2gp%name)
      IF (swid == -1) THEN
         msr = 'Failed to create swath ' // l2gp%name
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Define dimensions

      status = swdefdim(swid, DIM_NAME1, l2gp%nTimes)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF (l2gp%nLevels > 0) THEN
         status = swdefdim(swid, DIM_NAME2, l2gp%nLevels)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF (l2gp%nFreqs > 0) THEN
         status = swdefdim(swid, DIM_NAME3, l2gp%nFreqs)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME3
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Define horizontal geolocation fields using above dimensions

      status = swdefgfld(swid, GEO_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD3, DIM_NAME1, DFNT_FLOAT64, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD5
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD6, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD6
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD7, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD7
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefgfld(swid, GEO_FIELD8, DIM_NAME1, DFNT_INT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD8
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      IF (l2gp%nLevels > 0) THEN
         status = swdefgfld(swid, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

      IF (l2gp%nFreqs > 0) THEN
         status = swdefgfld(swid, GEO_FIELD10, DIM_NAME3, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD10
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
      ENDIF

! Define data fields using above dimensions

      IF (l2gp%nFreqs > 0) THEN

         status = swdefdfld(swid, DATA_FIELD5, DIM_NAME123, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD5 // ' for 3D quantity.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ELSE IF (l2gp%nLevels > 0) THEN

         status = swdefdfld(swid, DATA_FIELD5, DIM_NAME12, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD5 // ' for 2D quantity.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ELSE

         status = swdefdfld(swid, DATA_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD5 // ' for 1D quantity.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

      status = swdefdfld(swid, DATA_FIELD3, DIM_NAME1, DFNT_CHAR8, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swdefdfld(swid, DATA_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Detach from the swath interface.  This stores the swath info within the file
! and must be done before writing or reading data to or from the swath.

      status = swdetach(swid)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                    &detach from swath interface after L3Residual definition.')

!-------------------------------
   END SUBROUTINE ResidualCreate
!-------------------------------

!---------------------------------------
   SUBROUTINE ResidualWrite (l2gp, swid)
!---------------------------------------

! Brief description of subroutine
! This subroutine writes the data fields to the L3Residual swath.

! Arguments

      TYPE( L2GPData_T ), INTENT(IN) :: l2gp

      INTEGER, INTENT(IN) :: swid

! Parameters

! Functions

      INTEGER, EXTERNAL :: swwrfld

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: status
      INTEGER :: start(3), stride(3), edge(3)

! Write to the geolocation fields

      start = 0
      stride = 1
      edge(1) = l2gp%nFreqs
      edge(2) = l2gp%nLevels
      edge(3) = l2gp%nTimes

      if (l2gp%nTimes .gt. 0) then

      status = swwrfld( swid, GEO_FIELD1, start(3), stride(3), edge(3), &
                        REAL(l2gp%latitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, GEO_FIELD2, start(3), stride(3), edge(3), &
                        REAL(l2gp%longitude) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD2
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swid, GEO_FIELD3, start(3), stride(3), edge(3), &
                       l2gp%time)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, GEO_FIELD4, start(3), stride(3), edge(3), &
                        REAL(l2gp%solarTime) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, GEO_FIELD5, start(3), stride(3), edge(3), &
                        REAL(l2gp%solarZenith) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD5
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, GEO_FIELD6, start(3), stride(3), edge(3), &
                        REAL(l2gp%losAngle) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD6
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, GEO_FIELD7, start(3), stride(3), edge(3), &
                        REAL(l2gp%geodAngle) )
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD7
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld(swid, GEO_FIELD8, start(3), stride(3), edge(3), &
                       l2gp%chunkNumber)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD8
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      endif

      IF (l2gp%nLevels > 0) THEN

         status = swwrfld( swid, GEO_FIELD9, start(2), stride(2), edge(2), &
                           REAL(l2gp%pressures) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD9
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

      IF (l2gp%nFreqs > 0) THEN

         status = swwrfld( swid, GEO_FIELD10, start(1), stride(1), edge(1), &
                           REAL(l2gp%frequency) )
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD10
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDIF

! Write to the data fields

      IF (l2gp%nFreqs > 0) THEN

! L3Residual value is a 3-D field

         status = swwrfld( swid, DATA_FIELD5, start, stride, edge, &
                           REAL(l2gp%l2gpValue) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ELSE IF (l2gp%nLevels > 0) THEN

! L3Residual is 2-D

         status = swwrfld( swid, DATA_FIELD5, start(2:3), stride(2:3), &
                           edge(2:3), REAL(l2gp%l2gpValue(1,:,:)) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ELSE

! 1-D residual
      if (l2gp%nTimes .gt. 0) then

         status = swwrfld( swid, DATA_FIELD5, start(3), stride(3), edge(3), &
                           REAL(l2gp%l2gpValue(1,1,:)) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      endif

      ENDIF

! 1-D status & quality fields

      status = swwrfld(swid, DATA_FIELD3, start(3), stride(3), edge(3), &
                       l2gp%status)
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELD3
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = swwrfld( swid, DATA_FIELD4, start(3), stride(3), edge(3), &
                        REAL(l2gp%quality) )
      IF (status == -1) THEN
         msr = WR_ERR // DATA_FIELD4
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!------------------------------
   END SUBROUTINE ResidualWrite
!------------------------------

!---------------------------------------
   SUBROUTINE ResidualOutput (type, l3r)
!---------------------------------------

! Brief description of subroutine
! This program creates & writes the L3Residual swaths in an l3 map file.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      TYPE( L2GPData_T ), INTENT(IN) :: l3r(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: swattach, swclose, swdetach, swopen

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: l3File

      INTEGER :: i, l3id, match, numDays, status, swid

! For each element of the database,

      numDays = SIZE(l3r)

      DO i = 1, numDays

! Find the l3dm file for the proper day

         CALL FindFileDay(type, l3r(i)%time(1), mlspcf_l3dm_start, &
                          mlspcf_l3dm_end, match, l3File, date)
         IF (match == -1) THEN
            msr = 'No ' // TRIM(type) // ' file for day ' // date // ' found &
                  &for appending ' // TRIM(l3r(i)%name) // 'swath.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Re-open the l3dm grid file for the creation of swaths

         l3id = swopen(l3File, DFACC_RDWR)
         IF (l3id == -1) THEN
            msr = MLSMSG_Fileopen // l3File
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Create the L3Residual swath in this day's file

         CALL ResidualCreate( l3id, l3r(i) )

! Re-attach to the swath for writing

         swid = swattach(l3id, l3r(i)%name)
         IF (status == -1) THEN
            msr = 'Failed to re-attach to swath interface for writing to ' // &
                   l3r(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Write to the geolocation & data fields

         CALL ResidualWrite(l3r(i), swid)
         msr = TRIM(l3r(i)%name) // ' swath successfully written to file ' // &
               l3File
         CALL MLSMessage(MLSMSG_Info, ModuleName, msr)

! After writing, detach from swath interface

         status = swdetach(swid)
         IF (status == -1) THEN
            msr = 'Failed to detach from swath interface after writing to ' &
                  // l3r(i)%name
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Close the file

         status = swclose(l3id)
         IF (status == -1) THEN
            msr = 'Failed to close file ' // l3File // ' after writing.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

      ENDDO

!-------------------------------
   END SUBROUTINE ResidualOutput
!-------------------------------

!-----------------------------------------------------------------------
   SUBROUTINE ReadL2DGData (product, numFiles, files, numDays, prodL2GP)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reads l2gp dg files for a quantity over a range of days.  It
! returns an array of L2GPData_T structures, with one structure per day.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: product

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: files(:)

      INTEGER, INTENT(IN) :: numFiles

      INTEGER, INTENT(OUT) :: numDays

      TYPE( L2GPData_T ), POINTER :: prodL2GP(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: swclose, swinqswath, swopen

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=6000) :: list
      CHARACTER (LEN=GridNameLen) :: swath
      CHARACTER (LEN=FileNameLen) :: prodFiles(numFiles)

      INTEGER :: err, i, indx, j, l2id, len, ns, status

! Check for the product in the input files

      numDays = 0

      DO i = 1, numFiles

         ns = swinqswath(files(i), list, len)
         IF (ns == -1) THEN
            msr = 'Failed to read swath list from ' // files(i)
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
            CYCLE
         ENDIF

         DO j = 1, ns

            indx = INDEX(list, ',')

            IF (indx == 0) THEN
               swath = list
            ELSE
               swath = list(:indx-1)
               list = list(indx+1:)
            ENDIF

            IF ( TRIM(swath) == TRIM(product) ) THEN
               numDays = numDays + 1
               prodFiles(numDays) = files(i)
               EXIT
            ENDIF

         ENDDO

      ENDDO

! Check that numDays is at least 1, so pointer allocation is possible

      IF (numDays == 0) THEN

! If not, issue a warning

         msr = 'No L2GP data found for ' // product
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         NULLIFY(prodL2GP)

! If so, open, read, & close the files for this quantity

      ELSE

         ALLOCATE( prodL2GP(numDays), STAT=err )
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' l2gp array for ' // product
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         DO i = 1, numDays

            l2id = swopen(prodFiles(i), DFACC_READ)
            IF (l2id == -1) THEN
               msr = MLSMSG_Fileopen // prodFiles(i)
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

! Read information from the L2GP file

            CALL ReadL2GP( l2id, product, prodL2GP(i) )

! Close the L2GP file

            status = swclose(l2id)
            IF (status == -1) THEN
               msr = 'Failed to close file ' // prodFiles(i)
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

         ENDDO

      ENDIF

!-----------------------------
   END SUBROUTINE ReadL2DGData
!-----------------------------

!=====================
END MODULE L2Interface
!=====================

!# $Log$
!# Revision 1.7  2002/02/20 19:21:15  ybj
!# *** empty log message ***
!#
!# Revision 1.6  2001/07/18 15:50:59  nakamura
!# Generalized to work with L3M as well.
!#
!# Revision 1.5  2001/04/24 19:37:30  nakamura
!# Changes for privatization of L2GPData.
!#
!# Revision 1.4  2001/02/21 20:37:18  nakamura
!# Changed MLSPCF to MLSPCF3.
!#
!# Revision 1.3  2000/12/29 20:40:04  nakamura
!# Changed ReadL2GPProd to take start & end days as input; modified ResidualOutput for the one-product/all-days paradigm.
!#
!# Revision 1.2  2000/12/07 19:28:27  nakamura
!# Updated for level becoming an expandable template field.
!#
!# Revision 1.1  2000/11/15 20:56:51  nakamura
!# Module containing subroutines related to L2GP data.
!#
