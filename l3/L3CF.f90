
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3CF
!===============================================================================

   USE MLSCF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE MLSStrings
   USE PCFModule
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

! Definition -- L3CFDef_T
!               L3CFProd_T
! Subroutines -- CalculateArray
!                FillL3CF

! Remarks:  This module defines a data type to hold L3CF input and contains
!           subroutines used to process these data into forms used by L3.

! Parameters

! This data type is used to store global definitions read from the l3cf.

   TYPE L3CFDef_T

     INTEGER :: nNom		! number of points in the L2 nominal latitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: l2nomLats		! dimensioned (nNom)

     CHARACTER (LEN=FileNameLen) :: logType
	! template for the log file name

     CHARACTER (LEN=FileNameLen) :: stdType
	! template for the l3dz file for Standard products

     CHARACTER (LEN=FileNameLen) :: dgType
	! template for the l3dz file for Diagnostic products

     CHARACTER (LEN=FileNameLen) :: stdMCFname	! Standard MCF name
     CHARACTER (LEN=FileNameLen) :: dgMCFname	! Diagnostic MCF name

     INTEGER :: stdMCFnum			! Standard MCF number
     INTEGER :: dgMCFnum			! Diagnostic MCF number

     INTEGER :: minDays
	! # of days of input data needed to process an l3 product

   END TYPE L3CFDef_T

! This data type is used to store product-dependent cf input.

   TYPE L3CFProd_T

     ! DailyMap section

     CHARACTER (LEN=GridNameLen) :: l3prodNameD
	! name of product processed at l3

     INTEGER :: nDays			! number of days in the l3 processing window

     REAL(r8), DIMENSION(maxWindow) :: timeD		! dimensioned (nDays)
	! array of TAI93 synoptic times for each day in the processing window

     CHARACTER (LEN=3) :: mode					! asc/des/com/all

     INTEGER :: nLats 			! number of points in latitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: latGridMap	! dimensioned (nLats)

     INTEGER :: nLons			! number of points in longitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: longGrid	! dimensioned (nLons)

     REAL(r8), DIMENSION(2) :: l3presLvl       ! min, max combined pressure levels

     REAL(r8), DIMENSION(2) :: ascPresLvl      ! min, max ascending pressure levels

     REAL(r8), DIMENSION(2) :: desPresLvl      ! min, max descending pressure levels

     REAL(r8), DIMENSION(2) :: rangFrequency	! min, max of frequency range

     INTEGER, DIMENSION(2) :: rangWavenumber	! min, max of wave # range

     CHARACTER (LEN=3) :: intpMethod	! method of interpolation (lin/csp)

     ! Daily Zonal Mean specifications

     CHARACTER (LEN=2) :: calMethod			! calculation method (l2/l3)

     REAL(r8), DIMENSION(2) :: zComLvl		! min, max com pressure levels

     REAL(r8), DIMENSION(2) :: zAscLvl		! min, max asc pressure levels

     REAL(r8), DIMENSION(2) :: zDesLvl		! min, max des pressure levels

     ! Output section

     CHARACTER (LEN=FileNameLen) :: mcfName	! MCF name

     INTEGER :: mcfNum				! PCF number for the MCF

     CHARACTER (LEN=FileNameLen) :: fileTemplate	! file name from CF

   END TYPE L3CFProd_T

CONTAINS

!-----------------------------------------------------------------
   SUBROUTINE CalculateArray (start, max, delta, array, numPoints)
!-----------------------------------------------------------------

! Brief description of subroutine
! This subroutine calculates array an array of values, and the number of points
! in it, given the starting value, maximum value an element may take, and the
! desired increment between elements.

! Arguments

      REAL(r8), INTENT(IN) :: delta, max, start

      INTEGER, INTENT(OUT) :: numPoints

      REAL(r8), INTENT(OUT) :: array(:)

! Parameters

! Functions

! Variables

      INTEGER :: i

      array(1) = start

      DO i = 2, maxGridPoints

         array(i) = array(i-1) + delta

         IF (array(i) > max) EXIT

      ENDDO

      numPoints = i-1

!-------------------------------
   END SUBROUTINE CalculateArray
!-------------------------------

!-------------------------------------------------------------------
   SUBROUTINE FillL3CF (cf, pcfL3Ver, l3Days, l3Window, l3cf, cfDef)
!-------------------------------------------------------------------
! Brief description of subroutine
! This subroutine checks the parser output and fills L3CFProd_T.

! Arguments

      TYPE( Mlscf_T ), INTENT(IN) :: cf

      CHARACTER (LEN=*), INTENT(IN) :: pcfL3Ver

      CHARACTER (LEN=*), INTENT(IN) :: l3Days(:)

      INTEGER, INTENT(IN) :: l3Window

      TYPE( L3CFDef_T ), INTENT(OUT) :: cfDef

      TYPE( L3CFProd_T ), POINTER  :: l3cf(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_td_utcToTAI

! Variables

      CHARACTER (LEN=8) :: hhmmss
      CHARACTER (LEN=10) :: l3Ver, label
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=CCSDSB_LEN) :: timeB(maxWindow)
      CHARACTER (LEN=FileNameLen) :: match

      INTEGER :: err, i, iGlob, iLab, iMap, iOut, indx, j
      INTEGER :: mlspcf_mcf, nNon, numProds, returnStatus

      REAL(r8) :: start, end, delta
      REAL(r8) :: non(maxGridPoints)

! Find the section indices in the CF

      iGlob = LinearSearchStringArray(cf%Sections%MlscfSectionName, &
                                      'GlobalSettings')
      iMap = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'DailyMap')
      iOut = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'Output')

! Check that the version number given in the CF will match the PCF file names

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, &
                                     'OutputVersionString')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                     &the CF for OutputVersionString.')
      l3Ver = cf%Sections(iGlob)%Cells(indx)%CharValue
      IF (l3Ver /= pcfL3Ver) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                  'Output versions in the CF and PCF differ.')

! Find the non-negative values for the L2 nominal latitude grid in the
! GlobalSettings section of the cf; add the negative values, and save all in
! L3CFDef_T

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
      cfDef%l2nomLats(nNon:cfDef%nNom) = non(:nNon)

      DO i = 1, nNon-1
         cfDef%l2nomLats(i) = -non( nNon-(i-1) )
      ENDDO
 
! Find/save the log type from the GlobalSettings section of the cf

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword,'LogType')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                    &the CF for LogType.')
      cfDef%logType = cf%Sections(iGlob)%Cells(indx)%CharValue

! Find/save minDays from GlobalSettings

      indx = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword,'MinDays')
      IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                    &the CF for MinDays.')
      cfDef%minDays = cf%Sections(iGlob)%Cells(indx)%RealValue

! Find the L3DZS label

       label = 'L3DZS'
       iLab = LinearSearchStringArray( &
                     cf%Sections(iOut)%Entries%MlscfLabelName, TRIM(label) )
       IF (iLab == 0) THEN
          msr = 'Missing l3cf Output specification for ' // label
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Find the Standard MCF file name; save it in L3CFDef_T

       indx = LinearSearchStringArray( &
                          cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'mcf')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                              &keyword MCF in the Output section of the l3cf.')
       cfDef%stdMCFName = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Search for a match in the PCF

       CALL SearchPCFNames(cfDef%stdMCFName, mlspcf_mcf_l3dz_start, &
                           mlspcf_mcf_l3dz_end, mlspcf_mcf, match)
       IF (mlspcf_mcf == -1) THEN
          msr = 'No match in the PCF for file ' // cfDef%stdMCFName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Save the PCF number in L3CFProd_T

       cfDef%stdMCFNum = mlspcf_mcf

! Find the file template from Output; save it in L3CFDef_T

       indx = LinearSearchStringArray( &
                         cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'file')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                             &keyword FILE in the Output section of the l3cf.')
       cfDef%stdType = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Find the L3DZD label

       label = 'L3DZD'
       iLab = LinearSearchStringArray( &
                     cf%Sections(iOut)%Entries%MlscfLabelName, TRIM(label) )
       IF (iLab == 0) THEN
          msr = 'Missing l3cf Output specification for ' // label
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Find the Diagnostic MCF file name; save it in L3CFDef_T

       indx = LinearSearchStringArray( &
                          cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'mcf')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                              &keyword MCF in the Output section of the l3cf.')
       cfDef%dgMCFName = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Search for a match in the PCF

       CALL SearchPCFNames(cfDef%dgMCFName, mlspcf_mcf_l3dz_start, &
                           mlspcf_mcf_l3dz_end, mlspcf_mcf, match)
       IF (mlspcf_mcf == -1) THEN
          msr = 'No match in the PCF for file ' // cfDef%dgMCFName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Save the PCF number in L3CFProd_T

       cfDef%dgMCFNum = mlspcf_mcf

! Find the file template from Output; save it in L3CFDef_T

       indx = LinearSearchStringArray( &
                         cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'file')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                             &keyword FILE in the Output section of the l3cf.')
       cfDef%dgType = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Find the number of products for which L3 processing was requested

      numProds = cf%Sections(iMap)%NoSectionEntries

      ALLOCATE (l3cf(numProds), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' array of L3CFProd_T pointers.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Fill L3CFProd_T from the DailyMap & Output sections of the cf

      DO i = 1, numProds

! Find/save the product name

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'l3prodNameD')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                    &keyword L3PRODNAMED in the DailyMap section of the l3cf.')
       l3cf(i)%l3prodNameD = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Find the input synoptic time

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'timeD')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                          &keyword TIMED in the DailyMap section of the l3cf.')
       hhmmss = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

       l3cf(i)%nDays = l3Window

! For each DOY in the processing window, calculate the TAI time

       DO j = 1, l3Window
          timeB(j) = l3Days(j) // 'T' // hhmmss // '.000000Z'
          returnStatus = pgs_td_utcToTAI( timeB(j), l3cf(i)%timeD(j) )
          IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                            ModuleName, 'Error converting from CCSDSB to TAI.')
       ENDDO

! Find/save mode

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'mode')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                           &keyword MODE in the DailyMap section of the l3cf.')
       l3cf(i)%mode = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Find the boundaries of the latitude grid

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'latGridMap')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                     &keyword LATGRIDMAP in the DailyMap section of the l3cf.')
       start = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       end = cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the increment of the latitude grid

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'dLat')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                           &keyword DLAT in the DailyMap section of the l3cf.')
       delta = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue

! Calculate the array containing the latitude grid

       CALL CalculateArray(start, end, delta, l3cf(i)%latGridMap, &
                           l3cf(i)%nLats)

! Find the boundaries of the longitude grid

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'longGrid')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                       &keyword LONGGRID in the DailyMap section of the l3cf.')
       start = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       end = cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the increment of the longitude grid

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'dLon')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                           &keyword DLON in the DailyMap section of the l3cf.')
       delta = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue

! Calculate the array containing the longitude grid

       CALL CalculateArray(start, end, delta, l3cf(i)%longGrid, l3cf(i)%nLons)

! Find/save the min & max of the pressure levels for the desired modes

       IF ( (l3cf(i)%mode == 'com') .OR. (l3cf(i)%mode == 'all') ) THEN
          indx = LinearSearchStringArray( &
                        cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'l3presLvl')
          IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                         &keyword L3PRESLVL in the DailyMap section of the l3cf.')
          l3cf(i)%l3presLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
          l3cf(i)%l3presLvl(2) = &
                         cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound
       ELSE
          l3cf(i)%l3presLvl = 0.0
       ENDIF

       IF ( (l3cf(i)%mode == 'asc') .OR. (l3cf(i)%mode == 'all') ) THEN

          indx = LinearSearchStringArray( &
                        cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'ascPresLvl')
          IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                         &keyword ASCPRESLVL in the DailyMap section of the l3cf.')
          l3cf(i)%ascPresLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
          l3cf(i)%ascPresLvl(2) = &
                         cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound
       ELSE
          l3cf(i)%ascPresLvl = 0.0
       ENDIF

      IF ( (l3cf(i)%mode == 'des') .OR. (l3cf(i)%mode == 'all') ) THEN

          indx = LinearSearchStringArray( &
                        cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'desPresLvl')
          IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                         &keyword DESPRESLVL in the DailyMap section of the l3cf.')
          l3cf(i)%desPresLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
          l3cf(i)%desPresLvl(2) = &
                         cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

      ELSE
          l3cf(i)%desPresLvl = 0.0
      ENDIF

! Find/save the min & max for the ranges of frequencies, wavenumbers

       indx = LinearSearchStringArray( &
                   cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'rangFrequency')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                  &keyword RANGFREQUENCY in the DailyMap section of the l3cf.')
       l3cf(i)%rangFrequency(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%rangFrequency(2) = &
                            cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

       indx = LinearSearchStringArray( &
               cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'rangWavenumber')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                 &keyword RANGWAVENUMBER in the DailyMap section of the l3cf.')
       l3cf(i)%rangWavenumber(1) = &
                             cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%rangWavenumber(2) = &
                       cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find/save the interpolation method

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'intpMethod')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                     &keyword INTPMETHOD in the DailyMap section of the l3cf.')
       l3cf(i)%intpMethod = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Find the label in the DailyMap section

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'label')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                          &keyword LABEL in the DailyMap section of the l3cf.')
       label = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Match it to one in the Output section

       iLab = LinearSearchStringArray( &
                     cf%Sections(iOut)%Entries%MlscfLabelName, TRIM(label) )
       IF (iLab == 0) THEN
          msr = 'Missing l3cf Output specification for ' // label
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Find the MCF file name; save it in L3CFProd_T

       indx = LinearSearchStringArray( &
                          cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'mcf')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                              &keyword MCF in the Output section of the l3cf.')
       l3cf(i)%mcfName = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Search for a match in the PCF

       CALL SearchPCFNames(l3cf(i)%mcfName, mlspcf_mcf_l3dm_start, &
                           mlspcf_mcf_l3dm_end, mlspcf_mcf, match)
       IF (mlspcf_mcf == -1) THEN
          msr = 'No match in the PCF for file ' // l3cf(i)%mcfName
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

! Save the PCF number in L3CFProd_T

       l3cf(i)%mcfNum = mlspcf_mcf

! Find the file template from Output; save it in L3CFProd_T

       indx = LinearSearchStringArray( &
                         cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'file')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                             &keyword FILE in the Output section of the l3cf.') 
       l3cf(i)%fileTemplate = &
                          cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Daily Zonal Mean specifications -- find/save the DZ calculation method

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'calMethod')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                     &keyword CALMETHOD in the DailyMap section of the l3cf.')
       l3cf(i)%calMethod = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

! Find/save the DZ min:max pressure levels

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'zComLvl')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                      &keyword ZCOMLVL in the DailyMap section of the l3cf.')
       l3cf(i)%zComLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%zComLvl(2) = &
                      cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'zAscLvl')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                      &keyword ZASCLVL in the DailyMap section of the l3cf.')
       l3cf(i)%zAscLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%zAscLvl(2) = &
                      cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'zDesLvl')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                      &keyword ZDESLVL in the DailyMap section of the l3cf.')
       l3cf(i)%zDesLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%zDesLvl(2) = &
                      cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

    ENDDO

!-------------------------
   END SUBROUTINE FillL3CF
!-------------------------

!==============
END MODULE L3CF
!==============

! $Log$
! Revision 1.11  2001/05/04 18:29:16  nakamura
! Added DZ stuff.
!
! Revision 1.10  2001/03/27 17:26:11  nakamura
! Moved maxGridPoints to MLSL3Common.
!
! Revision 1.9  2001/02/21 20:43:43  nakamura
! Changed MLSPCF to MLSPCF3; removed InputVersion & check; added L3DZ stuff.
!
! Revision 1.8  2001/01/18 16:51:07  nakamura
! Added type L3CFDef_T; moved minDays from PCF to cf.
!
! Revision 1.7  2001/01/16 17:41:49  nakamura
! Added mcfName to L3CFProd_T; extract logType from cf.
!
! Revision 1.6  2000/12/29 20:44:46  nakamura
! Replaced USE of L3DMData with MLSL3Common; removed bpFlag, nGrids, quantities from L3CFProd_T.
!
! Revision 1.5  2000/12/07 21:21:50  nakamura
! Changed negative return from SearchPCFNames to be a PCF number of -1.
!
! Revision 1.4  2000/12/07 19:31:38  nakamura
! Updated for file name return from SearchPCFNames.
!
! Revision 1.3  2000/11/15 21:05:55  nakamura
! Moved parameter GridNameLen to module L3DMData; changed name of mlspcf_mcf start & end numbers.
!
! Revision 1.2  2000/10/24 19:18:03  nakamura
! Removed temporary subroutine ReadParseMLSCF; added PCF number for MCF to L3CFProd_T.
!
! Revision 1.1  2000/10/17 20:20:04  nakamura
! Module for customizing the CF parser output to L3.
!
