
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L2Interface
!===============================================================================

   USE SDPToolkit
   USE MLSMessageModule
   USE MLSCommon
   USE MLSL3Common
   USE L2GPData
   USE Hdf
   USE MLSPCF
   USE MLSStrings
   USE OpenInit
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
!                ReadL2GPProd
!                ResidualCreate
!                ResidualWrite
!                ResidualOutput

! Remarks:  This module contains subroutines involving the interface between
!           Level 3 and L2.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: DATA_FIELD5 = 'L3Residual'

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE GetL2GPfromPCF (template, startDOY, endDOY, numFiles, pcfNames)
!----------------------------------------------------------------------------

! Brief description of subroutine
! This routine searches the PCF for l2gp files for a desired species.  It
! returns an array of names and the number of the files found.

! Arguments

      CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY

      CHARACTER (LEN=*), INTENT(IN) :: template

      INTEGER, INTENT(OUT) :: numFiles

      CHARACTER (LEN=FileNameLen) :: pcfNames(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=8) :: date
      CHARACTER (LEN=FileNameLen) :: physicalFilename, type

      INTEGER :: i, indx, returnStatus, version

! Expand the level/species template

      CALL ExpandFileTemplate(template, type, 'L2GP')

      numFiles = 0

! Loop through all the PCF numbers for L2GP files

      DO i = mlspcf_l2gp_start, mlspcf_l2gp_end

         version = 1
         returnStatus = Pgs_pc_getReference(i, version, physicalFilename)

! If no file name was returned, go on to the next PCF number

         IF (returnStatus /= PGS_S_SUCCESS) CYCLE

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

!-----------------------------------------------------------------------
   SUBROUTINE ReadL2GPProd (cfProd, startDOY, endDOY, numDays, prodL2GP)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reads l2gp files for a quantity over a range of days.  It
! returns an array of L2GPData_T structures, with one structure per day.

! Arguments

      TYPE( L3CFProd_T ), INTENT(IN) :: cfProd

      CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY

      INTEGER, INTENT(OUT) :: numDays

      TYPE( L2GPData_T ), POINTER :: prodL2GP(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: swclose, swopen

! Variables

      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=FileNameLen) :: pcfNames(40)

      INTEGER :: err, i, l2id, numProfs, status

! Search the PCF for L2GP files for this species

      CALL GetL2GPfromPCF(cfProd%fileTemplate, startDOY, endDOY, numDays, &
                          pcfNames)

! Check that numDays is at least 1, so pointer allocation is possible

      IF (numDays == 0) THEN

! If not, issue a warning

         msr = 'No L2GP data found for ' // cfProd%l3prodNameD
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         NULLIFY(prodL2GP)

! If so, open, read, & close the files for this quantity

      ELSE

         ALLOCATE( prodL2GP(numDays), STAT=err )
         IF ( err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' l2gp array for ' // cfProd%l3prodNameD
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

         DO i = 1, numDays

            l2id = swopen(pcfNames(i), DFACC_READ)
            IF (l2id == -1) THEN
               msr = MLSMSG_Fileopen // pcfNames(i)
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

! Read information from the L2GP file

            CALL ReadL2GPData(l2id, cfProd%l3prodNameD, prodL2GP(i), numProfs)

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

         status = swwrfld( swid, DATA_FIELD5, start(3), stride(3), edge(3), &
                           REAL(l2gp%l2gpValue(1,1,:)) )
         IF (status == -1) THEN
            msr = WR_ERR // DATA_FIELD5
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

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

!=====================
END MODULE L2Interface
!=====================

!# $Log$
!# Revision 1.2  2000/12/07 19:28:27  nakamura
!# Updated for level becoming an expandable template field.
!#
!# Revision 1.1  2000/11/15 20:56:51  nakamura
!# Module containing subroutines related to L2GP data.
!#
