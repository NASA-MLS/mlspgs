
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE mon_L3CF
!===============================================================================

   USE L3CF, ONLY: CalculateArray
   USE MLSCF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSStrings
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definitions -- L3CFMDef_T
!                L3CFMProd_T
! Subroutines -- FillL3CFM

! Remarks:  This module defines a data type to hold L3CF input and contains
!           subroutines used to process these data into forms used by L3.

! Parameters

! This data type is used to store global definitions from the monthly l3cf.

   TYPE L3CFMDef_T

     INTEGER :: nNom		! number of points in the L2 nominal latitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: l2nomLats		! dimensioned (nNom)

     INTEGER :: minDays
        ! # of days of input data needed to process an l3 product

     CHARACTER (LEN=3) :: intpMethod		! method of interpolation (lin/csp)

     CHARACTER (LEN=FileNameLen) :: l2dgType
	! template for the l2 diagnostics input file name

     CHARACTER (LEN=FileNameLen) :: stdType
	! template for the l3dz file for Standard products

     CHARACTER (LEN=FileNameLen) :: dgType
	! template for the l3dz file for Diagnostic products

   END TYPE L3CFMDef_T

! This data type is used to store product-dependent cf input.

   TYPE L3CFMProd_T

     ! Standard section

     CHARACTER (LEN=GridNameLen) :: l3prodName
	! name of product processed at l3

     CHARACTER (LEN=3) :: mode		! asc/des/com/all/ado

     INTEGER :: nLats 			! number of points in latitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: latGridMap      ! dimensioned (nLats)

     INTEGER :: nLons			! number of points in longitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: longGrid        ! dimensioned (nLons)

     REAL(r8), DIMENSION(2) :: l3presLvl      ! min, max pressure levels, combined

     REAL(r8), DIMENSION(2) :: ascPresLvl     ! min, max pressure levels, ascending

     REAL(r8), DIMENSION(2) :: desPresLvl     ! min, max pressure levels, descending

     ! Output section

     CHARACTER (LEN=FileNameLen) :: fileTemplate	! file template for prod

   END TYPE L3CFMProd_T

CONTAINS

!---------------------------------------------------------
   SUBROUTINE FillL3CFM (cf, pcfL3Ver, cfStd, cfDg, cfDef)
!---------------------------------------------------------

! Brief description of subroutine
! This subroutine checks the parser output and fills L3CFMProd_T & L3CFMDef_T.

! Arguments

      TYPE( Mlscf_T ), INTENT(IN) :: cf

      CHARACTER (LEN=*), INTENT(IN) :: pcfL3Ver

      TYPE( L3CFMDef_T ), INTENT(OUT) :: cfDef

      TYPE( L3CFMProd_T ), POINTER  :: cfDg(:), cfStd(:)

! Parameters

! Functions

! Variables

      CHARACTER (LEN=10) :: l3Ver, label
      CHARACTER (LEN=480) :: msr

      INTEGER :: err, i,iDg, iGlob, iLab, iMap, iOut, indx, nLat, nNon, numProds

      REAL(r8) :: start, end, delta
      REAL(r8) :: latGrid(maxGridPoints), non(maxGridPoints)

! Find the section indices in the CF

      iGlob = LinearSearchStringArray(cf%Sections%MlscfSectionName, &
                                      'GlobalSettings')
      IF (iGlob == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'GLOBALSETTINGS &
                                                     &section missing from the CF.')
      iMap = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'Standard')
      IF (iMap == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'STANDARD section &
                                                             &missing from the CF.')
      iDg = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'Diagnostic')
      IF (iDg == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'DIAGNOSTIC section &
                                                             &missing from the CF.')
      iOut = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'Output')
      IF (iOut == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'OUTPUT section &
                                                             &missing from the CF.')

! Check that the version number given in the CF will match the PCF file names

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, &
                                     'OutputVersionString')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                    'No entry in the CF for OutputVersionString.')
      l3Ver = cf%Sections(iGlob)%Cells(indx)%CharValue
      IF (l3Ver /= pcfL3Ver) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Output versions in the CF and PCF differ.')

! Find the non-negative values for the L2 nominal latitude grid in the
! GlobalSettings section of the cf; add mirror-image negative values and save them
! in L3CFMDef_T

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, &
                                     'l2nomLats')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                     &the CF for l2nomLats.')

      DO i = 1, maxGridPoints
         non(i) = cf%Sections(iGlob)%Cells(indx)%RealValue
         IF (cf%Sections(iGlob)%Cells(indx)%More == 0) EXIT
         indx = indx + 1
      ENDDO
      nNon = i
      cfDef%nNom = 2*nNon - 1
      cfDef%l2nomLats = 0.0
      cfDef%l2nomLats(nNon:cfDef%nNom) = non(:nNon)

      DO i = 1, nNon-1
         cfDef%l2nomLats(i) = -non( nNon-(i-1) )
      ENDDO

! Find/save the l2dg type from the GlobalSettings section of the cf

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword,'L2dgType')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                    &the CF for L2dgType.')
      cfDef%l2dgType = cf%Sections(iGlob)%Cells(indx)%CharValue

! Find/save the l3dz std type

      indx = LinearSearchStringArray( &
                              cf%Sections(iGlob)%Cells%Keyword, 'stdType')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in the CF &
                            &for stdType.')
      cfDef%stdType = cf%Sections(iGlob)%Cells(indx)%CharValue

! Find/save the l3dz dg type

      indx = LinearSearchStringArray( &
                              cf%Sections(iGlob)%Cells%Keyword, 'dgType')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in the CF &
                            &for dgType.')
      cfDef%dgType = cf%Sections(iGlob)%Cells(indx)%CharValue

! Find/save minDays from the GlobalSettings section of the cf

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword,'MinDays')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                    'No entry in the CF for MinDays.')
      cfDef%minDays = cf%Sections(iGlob)%Cells(indx)%RealValue

! Find/save the interpolation method

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, 'intpMethod')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in the CF &
                                     &for intpMethod.')
      cfDef%intpMethod = cf%Sections(iGlob)%Cells(indx)%CharValue

! Find the boundaries of the latitude grid

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, 'latGridMap')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in the CF &
                     &for latGridMap.')
      start = cf%Sections(iGlob)%Cells(indx)%RealValue
      end = cf%Sections(iGlob)%Cells(indx)%RangeUpperBound

! Find the increment of the latitude grid

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, 'dLat')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in the CF &
                                     &for dLat.')
      delta = cf%Sections(iGlob)%Cells(indx)%RealValue

! Calculate the array containing the latitude grid

      CALL CalculateArray(start, end, delta, latGrid, nLat)

! Find the number of standard products for which L3 processing was requested

      numProds = cf%Sections(iMap)%NoSectionEntries

      ALLOCATE (cfStd(numProds), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' array of L3CFMProd_T pointers.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Fill L3CFMProd_T from the Standard & Output sections of the cf

      DO i = 1, numProds

! Find/save the product name

         indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'l3prodName')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
             'Missing keyword L3PRODNAME in the Standard section of the l3cf.')
         cfStd(i)%l3prodName = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Find/save mode

         indx = LinearSearchStringArray( &
                                 cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'mode')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                    'Missing keyword MODE in the Standard section of the l3cf.')
         cfStd(i)%mode = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Save the latitude grid quantities

         cfStd(i)%latGridMap = latGrid
         cfStd(i)%nLats = nLat

! Find the boundaries of the longitude grid

         indx = LinearSearchStringArray( &
                             cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'longGrid')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                'Missing keyword LONGGRID in the Standard section of the l3cf.')
         start = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
         end = cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the increment of the longitude grid

         indx = LinearSearchStringArray( &
                                 cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'dLon')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                    'Missing keyword DLON in the Standard section of the l3cf.')
         delta = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue

! Calculate the array containing the longitude grid

         CALL CalculateArray(start, end, delta, cfStd(i)%longGrid, cfStd(i)%nLons)

! Find/save the min & max for the pressure levels

         indx = LinearSearchStringArray( &
                            cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'l3presLvl')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                  'Missing keyword L3PRESLVL in the Standard section of the l3cf.')
         cfStd(i)%l3presLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
         cfStd(i)%l3presLvl(2) = &
                            cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

         indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'ascPresLvl')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
              'Missing keyword ASCPRESLVL in the Standard section of the l3cf.')
         cfStd(i)%ascPresLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
         cfStd(i)%ascPresLvl(2) = &
                            cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

         indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'desPresLvl')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                  'Missing keyword DESPRESLVL in the Standard section of the l3cf.')
         cfStd(i)%desPresLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
         cfStd(i)%desPresLvl(2) = &
                            cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the label in the Standard section

         indx = LinearSearchStringArray( &
                             cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'label')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                   'Missing keyword LABEL in the Standard section of the l3cf.')
         label = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Match it to one in the Output section

         iLab = LinearSearchStringArray( cf%Sections(iOut)%Entries%MlscfLabelName, &
                                         TRIM(label) )
         IF (iLab == 0) THEN
            msr = 'Missing l3cf Output specification for ' // label
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF

! Find the file template from Output; save it in L3CFMProd_T

         indx = LinearSearchStringArray( &
                         cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'file')
         IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                    'Missing keyword FILE in the Output section of the l3cf.') 
         cfStd(i)%fileTemplate = &
                               cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

      ENDDO

! Find the number of diagnostic products

      numProds = cf%Sections(iDg)%NoSectionEntries

      ALLOCATE (cfDg(numProds), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' array of L3CFDg_T pointers.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Fill L3CFDg_T from the Diagnostic section of the cf

      DO i = 1, numProds

! Set the unused fields to default values

         cfDg(i)%mode = 'com'
         cfDg(i)%latGridMap = latGrid
         cfDg(i)%nLats = nLat

! Find/save the product name

        indx = LinearSearchStringArray( &
                            cf%Sections(iDg)%Entries(i)%Cells%Keyword, 'l3prodName')
        IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                        &keyword L3PRODNAME in the Diagnostic section of the l3cf.')
        cfDg(i)%l3prodName = cf%Sections(iDg)%Entries(i)%Cells(indx)%CharValue

! Find/save the min & max of the pressure levels

        indx = LinearSearchStringArray( &
                            cf%Sections(iDg)%Entries(i)%Cells%Keyword, 'l3presLvl')
        IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                        &keyword L3PRESLVL in the Diagnostic section of the l3cf.')
        cfDg(i)%l3presLvl(1) = cf%Sections(iDg)%Entries(i)%Cells(indx)%RealValue
        cfDg(i)%l3presLvl(2) = cf%Sections(iDg)%Entries(i)%Cells(indx)%RangeUpperBound
        cfDg(i)%ascPresLvl = cfDg(i)%l3PresLvl
        cfDg(i)%desPresLvl = cfDg(i)%l3PresLvl

! Find the boundaries of the longitude grid

        indx = LinearSearchStringArray( &
                            cf%Sections(iDg)%Entries(i)%Cells%Keyword, 'longGrid')
        IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                 'Missing keyword LONGGRID in the Diagnostic section of the l3cf.')
        start = cf%Sections(iDg)%Entries(i)%Cells(indx)%RealValue
        end = cf%Sections(iDg)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the increment of the longitude grid

        indx = LinearSearchStringArray( &
                                  cf%Sections(iDg)%Entries(i)%Cells%Keyword, 'dLon')
        IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                      'Missing keyword DLON in the Diagnostic section of the l3cf.')
        delta = cf%Sections(iDg)%Entries(i)%Cells(indx)%RealValue

! Calculate the array containing the longitude grid

        CALL CalculateArray(start, end, delta, cfDg(i)%longGrid, cfDg(i)%nLons)

      ENDDO

!--------------------------
   END SUBROUTINE FillL3CFM
!--------------------------

!==================
END MODULE mon_L3CF
!==================

! $Log$
! Revision 1.1  2001/07/18 15:45:08  nakamura
! Module for customizing the CF parser output to the L3 Monthly program.
!
!
