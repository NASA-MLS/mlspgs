
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE PCFModule
!===============================================================================

   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- ExpandFileTemplate
!                SearchPCFNames
!                FindFileType
!                FindFileDay

! Remarks:  This module contains subroutines related to getting file names from
!           the PCF.

CONTAINS

!-----------------------------------------------------------------------------
   SUBROUTINE ExpandFileTemplate (template, filename, level, version, cycle, &
                                  day)
!-----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine expands the version, cycle, and day fields in a file name
! template.  Note: CYCLE is a STRING here, which is how it's read from the PCF.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: template

      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: cycle, day, level, version

      CHARACTER (LEN=*), INTENT(OUT) :: fileName

! Parameters

! Functions

! Variables

      CHARACTER (LEN=4) :: field

      INTEGER :: i, indx, numFields

! Initializations

      numFields = 4

      fileName = template

! Loop through the expandable fields in the template

      DO i = 1, numFields

! Search for $

         indx = INDEX(fileName,'$')

! Exit, if there are no expandable fields 

         IF (indx == 0) EXIT

! Match the field name to the input argument

         field = fileName(indx:indx+3)
    
         IF (field == '$lev') THEN

            IF ( PRESENT(level) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(level) // &
                          fileName((indx+6):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input level &
                                            &required to expand the template.')
            ENDIF

         ELSE IF (field == '$ver') THEN

            IF ( PRESENT(version) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(version) // &
                          fileName((indx+8):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input version &
                                            &required to expand the template.')
            ENDIF

         ELSE IF (field == '$cyc') THEN

            IF ( PRESENT(cycle) ) THEN
               fileName = fileName(:(indx-1)) // 'c' // TRIM(cycle) // &
                          fileName((indx+6):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input cycle &
                                            &required to expand the template.')
            ENDIF

         ELSE IF (field == '$day') THEN

            IF ( PRESENT(day) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(day) // &
                          fileName((indx+4):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input day required &
                                                    &to expand the template.')
            ENDIF

         ENDIF

      ENDDO

!-----------------------------------
   END SUBROUTINE ExpandFileTemplate
!-----------------------------------

!------------------------------------------------------------------------
   SUBROUTINE SearchPCFNames (inName, mlspcf_start, mlspcf_end, mlspcf, &
                              outName)
!------------------------------------------------------------------------

! Brief description of subroutine 
! This subroutine searches the PCF for an entry matching the input name.  If a
! match is found, it returns the file's name (including the path) & number from
! the PCF.  If no entry is found, a PCF number of -1 is returned, along with
! the original input name, with a path name concatenated to it.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: inName

      INTEGER, INTENT(IN) :: mlspcf_end, mlspcf_start

      INTEGER, INTENT(OUT) :: mlspcf

      CHARACTER (LEN=FileNameLen), INTENT(OUT) :: outName

! Parameters

! Functions

! Variables

      CHARACTER (LEN=FileNameLen) :: path, pcfName

      INTEGER :: i, indx, returnStatus, version

! Initializations

      mlspcf = -1

! Loop through all the files in this range of PCF numbers

      DO i = mlspcf_start, mlspcf_end

! Retrieve the PCF name

         version = 1

         returnStatus = Pgs_pc_getReference(i, version, pcfName)

         IF (returnStatus == PGS_S_SUCCESS) THEN

! If a file name was successfully retrieved from the PCF, extract its pathname

            indx = INDEX(pcfName, '/', .TRUE.)
            path = pcfName(:indx)

! Concatenate the path with the input name

            outName = TRIM(path) // inName

! Check the path/input name string against the PCF entry

            IF (outName == pcfName) THEN

! If it matches, save the PCF number

               mlspcf = i
               EXIT

            ENDIF

         ENDIF

      ENDDO

!-------------------------------
   END SUBROUTINE SearchPCFNames
!-------------------------------

!-----------------------------------------------------------------------
   SUBROUTINE FindFileType (type, mlspcf_start, mlspcf_end, match, name)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine searches the PCF for a file of the proper level/species.  It
! returns a matching PCF number & name, or a value of match = -1, if no match
! was found.

! Arguments

      CHARACTER(LEN=*), INTENT(IN) :: type

      INTEGER, INTENT(IN) :: mlspcf_start, mlspcf_end

      CHARACTER(LEN=*), INTENT(OUT) :: name

      INTEGER, INTENT(OUT) :: match

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: i, returnStatus, version

! Initialize match to its unfound value

      match = -1

! Loop through all PCF entries in the range corresponding to this file type

      DO i = mlspcf_start, mlspcf_end

! Retrieve the name

         version = 1
         returnStatus = Pgs_pc_getReference(i, version, name)

! If no name was returned for this number, go on to the next one

         IF (returnStatus /= PGS_S_SUCCESS) CYCLE

! Check whether the returned name contains the desired file substring

         IF (INDEX(name, TRIM(type)) /= 0) THEN

! If so, save the PCF number for the matching entry

            match = i
            EXIT

         ENDIF

      ENDDO

! If no match was found, this is an error message

      IF (match == -1) THEN
         msr = 'No PCF entry for a file containing the substring ' // type
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!-----------------------------
   END SUBROUTINE FindFileType
!-----------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE FindFileDay(type, time, mlspcf_start, mlspcf_end, match, name, &
                          date)
!-----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine searches the PCF for a file of a desired level/species/day.
! It returns a PCF number & name, the date, or match = -1 if no entry is found.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: type

      INTEGER, INTENT(IN) :: mlspcf_end, mlspcf_start

      REAL(r8), INTENT(IN) :: time

      CHARACTER (LEN=*), INTENT(OUT) :: name

      CHARACTER (LEN=8), INTENT(OUT) :: date

      INTEGER, INTENT(OUT) :: match

! Parameters

! Functions

      INTEGER, EXTERNAL :: Pgs_td_asciiTime_aToB

! Variables

      CHARACTER (LEN=26) :: timeB
      CHARACTER (LEN=CCSDS_LEN) :: timeA
      CHARACTER (LEN=480) :: msr

      INTEGER :: i, returnStatus, version

! Convert input time to CCSDS A format

      returnStatus = Pgs_td_taiToUTC(time, timeA)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                                                         ModuleName, TAI2A_ERR)

! Convert from CCSDS A to B

      returnStatus = Pgs_td_asciiTime_aToB(timeA, timeB)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                ModuleName,'Error converting data time from CCSDS A to B.')

! Save only the date portion, for comparison to file names

      date = timeB(1:8)

! Loop through all the files in this range of PCF numbers

      match = -1

      DO i = mlspcf_start, mlspcf_end

! Retrieve the name

         version = 1
         returnStatus = Pgs_pc_getReference(i, version, name)

! If no name was returned for this number, go on to the next one

         IF (returnStatus /= PGS_S_SUCCESS) CYCLE

! Check whether the returned name contains the desired level/species/day

         IF ((INDEX(name,TRIM(type)) /= 0) .AND. (INDEX(name,date) /= 0)) THEN

! If so, save the PCF number for the matching entry

            match = i
            EXIT

         ENDIF

      ENDDO

! If match was found, this is an error condition

      IF (match == -1) THEN
         msr = 'No PCF entry of the form ' // TRIM(type) // ' for day ' // date
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!----------------------------
   END SUBROUTINE FindFileDay
!----------------------------

!===================
END MODULE PCFModule
!===================

! $Log$
! Revision 1.4  2000/12/07 21:23:04  nakamura
! Changed SearchPCFNames so that input name is always returned with a path, negative PCF number means no file found.
!
! Revision 1.3  2000/12/07 19:43:58  nakamura
! Added 'level' to expandable fields in file template; extract file path names from PCF, rather than cf.
!
! Revision 1.2  2000/10/24 19:45:26  nakamura
! Corrected subroutine name for SearchPCFNames in module Contents; changed version from a constant to a variable in getRef.
!
! Revision 1.1  2000/10/17 20:31:58  nakamura
! Module for getting file names based on input in the PCF and CF.
!
