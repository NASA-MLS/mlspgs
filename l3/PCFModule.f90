
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE PCFModule
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   USE MLSPCF
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

! Remarks:  This module contains subroutines related to getting file names from
!           the PCF.

CONTAINS

!-------------------------------------------------------------------------
   SUBROUTINE ExpandFileTemplate (template, filename, version, cycle, day)
!-------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine expands the version, cycle, and day fields in a file name
! template.  Note: CYCLE is a STRING here, which is how it's read from the PCF.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: template

      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: cycle, day, version

      CHARACTER (LEN=*), INTENT(OUT) :: fileName

! Parameters

! Functions

! Variables

      CHARACTER (LEN=4) :: field

      INTEGER :: i, indx, numFields

! Initializations

      numFields = 3

      fileName = template

! Loop through the expandable fields in the template

      DO i = 1, numFields

! Search for $

         indx = INDEX(fileName,'$')

! Exit, if there are no expandable fields 

         IF (indx == 0) EXIT

! Match the field name to the input argument

         field = fileName(indx:indx+3)
    
         IF (field == '$ver') THEN

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

!-----------------------------------------------------------------------------
   SUBROUTINE SearchPCFNames (inName, mlspcf_start, mlspcf_end, mlspcf, match)
!-----------------------------------------------------------------------------

! Brief description of subroutine 
! This subroutine searches the PCF for an entry matching the input name, which
! includes the path.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: inName

      INTEGER, INTENT(IN) :: mlspcf_end, mlspcf_start

      INTEGER, INTENT(OUT) :: mlspcf

      INTEGER, INTENT(OUT) :: match

! Parameters

! Functions

! Variables

      CHARACTER (LEN=FileNameLen) :: pcfName

      INTEGER :: i, returnStatus, version

! Initializations

      mlspcf = -1
      match = 0

! Loop through all the files in this range of PCF numbers

      DO i = mlspcf_start, mlspcf_end

! Retrieve the PCF name

         version = 1

         returnStatus = Pgs_pc_getReference(i, version, pcfName)

         IF (returnStatus /= PGS_S_SUCCESS) CYCLE

! Check it against the given name

         IF (inName == pcfName) THEN
            mlspcf = i
            match = 1
            EXIT
         ENDIF

      ENDDO

!-------------------------------
   END SUBROUTINE SearchPCFNames
!-------------------------------

!===================
END MODULE PCFModule
!===================

! $Log$
! Revision 1.1  2000/10/17 20:31:58  nakamura
! Module for getting file names based on input in the PCF and CF.
!
