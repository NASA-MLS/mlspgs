
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L1BData
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L1BData_T
! Subroutines -- ReadL1BData
!                DeallocateL1BData

! Remarks:  This is a prototype module containing parameters, a derived type
!           definition, and a subroutine for reading the L1B data.

! Parameters

   INTEGER, PARAMETER :: NAME_LEN = 64		! Max len of SDS array name

! This data type is used to store quantities from an L1B data file.

   TYPE L1BData_T
      CHARACTER (LEN=NAME_LEN) :: L1BName	! Name of field in file
      INTEGER :: firstMAF   			! First major frame read
      INTEGER :: noMAFs				! # of MAFs read
      INTEGER :: maxMIFs			! Max # of MIFs/MAF in SD array
      INTEGER :: noAuxInds			! # of auxilliary indices

      INTEGER, DIMENSION(:), POINTER :: counterMAF	! dimensioned (noMAFs)

      CHARACTER, DIMENSION(:,:,:), POINTER :: charField
      REAL(r8), DIMENSION(:,:,:), POINTER :: dpField
      INTEGER, DIMENSION(:,:,:), POINTER :: intField
			! dimensioned (noAuxInds,maxMIFs,noMAFs)
   END TYPE L1BData_T

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE ReadL1BData (L1FileHandle, quantityName, l1bData, noMAFs, flag, &
                           firstMAF, lastMAF)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This is a prototype for the ReadL1BData subroutine in the Open/Init task.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: quantityName

      INTEGER, INTENT(IN) :: L1FileHandle

      INTEGER, INTENT(IN), OPTIONAL :: firstMAF, lastMAF
   
      TYPE( L1BData_T ), INTENT(OUT) :: l1bData

      INTEGER, INTENT(OUT) :: flag, noMAFs

! Parameters

      CHARACTER (LEN=*), PARAMETER :: INPUT_ERR = 'Error in input argument '

      INTEGER, PARAMETER :: MAX_VAR_DIMS = 32

! Functions

      INTEGER, EXTERNAL :: sfginfo, sfn2index, sfselect, sfrdata, sfrcdata
      INTEGER, EXTERNAL :: sfendacc

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, data_type, dealloc_err, i, n_attrs, numMAFs, rank
      INTEGER :: sds1_id, sds2_id, sds_index, status
      INTEGER :: dim_sizes(MAX_VAR_DIMS)
      INTEGER, ALLOCATABLE :: edge(:), start(:), stride(:)

      LOGICAL :: firstCheck, lastCheck

      REAL, ALLOCATABLE :: realField(:,:,:)

      flag = 0

! Find data sets for counterMAF & quantity by name

      sds_index = sfn2index(L1FileHandle, 'counterMAF')
      IF (sds_index == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                       &to find index of counterMAF data set.')

      sds1_id = sfselect(L1FileHandle, sds_index)
      IF (sds1_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                     &find identifier of counterMAF data set.')

      sds_index = sfn2index(L1FileHandle, quantityName)
      IF (sds_index == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                         &to find index of quantity data set.')

      sds2_id = sfselect(L1FileHandle, sds_index)
      IF (sds2_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                             &find identifier of data set matching the index.')

! Find rank (# of dimensions), dimension sizes of quantity data set

      status = sfginfo(sds2_id, quantityName, rank, dim_sizes, data_type, &
                       n_attrs)

      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                                      &find rank of data set.')

! Allocate, based on above SD, dim info

      ALLOCATE (edge(rank), start(rank), stride(rank), STAT=alloc_err)

      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  slab dimensions.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Set "slab" dimensions

      DO i = 1, rank
         edge(i) = dim_sizes(i)
      ENDDO
      edge(rank) = 1

      start = 0

      stride = 1

! Fill in "indexing" values of l1b object

      l1bData%L1BName = quantityName

      IF (rank < 2) THEN
         l1bData%maxMIFs = 1
      ELSE IF (rank > 2) THEN
         l1bData%maxMIFs = dim_sizes(2)
      ELSE
         l1bData%maxMIFs = dim_sizes(1)
      ENDIF

      IF (rank > 2) THEN
         l1bData%noAuxInds = dim_sizes(1)
      ELSE
         l1bData%noAuxInds = 1
      ENDIF

! Check input arguments, set noMAFs

      numMAFs = dim_sizes(rank)

      firstCheck = PRESENT(firstMAF)
      lastCheck = PRESENT(lastMAF)

      IF (firstCheck) THEN

          IF ( (firstMAF >= numMAFs) .OR. (firstMAF < 0) ) THEN
             msr = INPUT_ERR // 'firstMAF'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          l1bData%firstMAF = firstMAF

      ELSE

          l1bData%firstMAF = 0

      ENDIF

      IF (lastCheck) THEN

          IF (lastMAF < l1bData%firstMAF) THEN
             msr = INPUT_ERR // 'lastMAF'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          IF (lastMAF >= numMAFs) THEN
             l1bData%noMAFs = numMAFs - l1bData%firstMAF
          ELSE
             l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
          ENDIF

      ELSE

         l1bData%noMAFs = numMAFs - l1bData%firstMAF

      ENDIF

      noMAFs = l1bData%noMAFs

! Allocate, read counterMAF

      ALLOCATE(l1bData%counterMAF(l1bData%noMAFs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  counterMAF pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = sfrdata(sds1_id, l1bData%firstMAF, stride, l1bData%noMAFs, &
                       l1bData%counterMAF)
      IF (status == -1) THEN
         msr = MLSMSG_L1BRead // 'counterMAF.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Allocate, read according to field type; nullify unused pointers

      IF (data_type == 4) THEN

! character (DFNT_CHAR8)

         ALLOCATE ( l1bData%charField(l1bData%noAuxInds,l1bData%maxMIFs,&
                   &l1bData%noMAFs), STAT=alloc_err )
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' charField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%intField, l1bData%dpField)

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrcdata( sds2_id, start, stride, edge, &
                               l1bData%charField(:,:,i) )
            IF (status == -1) THEN
               msr = MLSMSG_L1BRead // quantityName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDDO

      ELSE IF (data_type == 24) THEN

! integer (DFNT_INT32)

         ALLOCATE (l1bData%intField(l1bData%noAuxInds,l1bData%maxMIFs,&
                  &l1bData%noMAFs), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' intField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%charField, l1bData%dpField)

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrdata( sds2_id, start, stride, edge, &
                              l1bData%intField(:,:,i) )
            IF (status == -1) THEN
               msr = MLSMSG_L1BRead // quantityName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDDO

      ELSE

! float (REAL or DOUBLE PRECISION)

         ALLOCATE (l1bData%dpField(l1bData%noAuxInds,l1bData%maxMIFs,&
                  &l1bData%noMAFs), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' dpField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%charField, l1bData%intField)

! If real (DFNT_FLOAT32), allocate real variable and then move to dp pointer

         IF (data_type == 5) THEN

            ALLOCATE (realField(l1bData%noAuxInds,l1bData%maxMIFs,&
                     &l1bData%noMAFs), STAT=alloc_err)
            IF ( alloc_err /= 0 ) THEN
               msr = MLSMSG_Allocate // ' realField pointer.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            DO i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata( sds2_id, start, stride, edge, &
                                realField(:,:,i) )
               IF (status == -1) THEN
                  msr = MLSMSG_L1BRead // quantityName
                  CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
               ENDIF
            ENDDO

            l1bData%dpField = DBLE(realField)

            DEALLOCATE(realField, STAT=dealloc_err)
            IF ( dealloc_err /= 0 ) THEN
               CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed &
                                     &deallocation of local real variable.')
               flag = -1
            ENDIF


         ELSE

! double precision (DFNT_FLOAT64)


            DO i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata(sds2_id, start, stride, edge, &
                                l1bData%dpField(:,:,i))
               IF (status == -1) THEN
                  msr = MLSMSG_L1BRead // quantityName
                  CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
               ENDIF
            ENDDO

         ENDIF

      ENDIF

! Terminate access to the data sets

      status = sfendacc(sds1_id)
      status = sfendacc(sds2_id)
      IF (status == -1) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed to terminate &
                                                     &access to data sets.')
         flag = -1
      ENDIF

! Deallocate local variables

      DEALLOCATE (edge, start, stride, STAT=dealloc_err )

      IF ( dealloc_err /= 0 ) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &local variables.')
         flag = -1
      ENDIF

!----------------------------
   END SUBROUTINE ReadL1BData
!----------------------------

!----------------------------------------------
   SUBROUTINE DeallocateL1BData (l1bData, flag)
!----------------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L1BData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L1BData_T ), intent(INOUT) :: l1bData

      INTEGER, INTENT(OUT) :: flag

! Parameters

! Functions

! Variables

      INTEGER :: dealloc_err

      flag = 0

      IF ( ASSOCIATED(l1bData%counterMAF) ) DEALLOCATE (l1bData%counterMAF, &
                                                        STAT=dealloc_err)
      IF ( ASSOCIATED(l1bData%charField) ) THEN
         DEALLOCATE (l1bData%charField, STAT=dealloc_err)
      ELSE IF ( ASSOCIATED(l1bData%intField) ) THEN
         DEALLOCATE (l1bData%intField, STAT=dealloc_err)
      ELSE
         DEALLOCATE (l1bData%dpField, STAT=dealloc_err)
      ENDIF

      IF ( dealloc_err /= 0 ) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &L1BData pointers.')
         flag = -1
      ENDIF

!----------------------------------
   END SUBROUTINE DeallocateL1BData
!----------------------------------

!=================
END MODULE L1BData
!=================

! $Log$
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.11  2000/06/06 21:11:55  vsnyder
! Declare EXTERNAL attribute for HDF's sf* functions
!
