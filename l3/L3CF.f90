
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L3CF
!===============================================================================

   USE L3DMData
   USE MLSCommon
   USE MLSCF
   USE MLSStrings
   USE MLSMessageModule
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

! Definition -- L3CFProd_T
! Subroutines -- CalculateArray
!                FillL3CF
!                ReadParseMLSCF

! Remarks:  This module defines a data type to hold L3CF input and contains
!           subroutines used to process these data into forms used by L3.

! Parameters

   INTEGER, PARAMETER :: GridNameLen = 64
   INTEGER, PARAMETER :: maxGridPoints = 500

! This data type is used to store cf input needed to process l3dm data.

   TYPE L3CFProd_T

     ! DailyMap section

     CHARACTER (LEN=GridNameLen) :: l3prodNameD
	! name of product processed at l3

     INTEGER :: nDays		! number of days in the l3 processing window

     REAL(r8), DIMENSION(maxWindow) :: timeD		! dimensioned (nDays)
	! array of TAI93 synoptic times for each day in the processing window

     CHARACTER (LEN=3) :: intpMethod	! method of interpolation (lin/csp)

     INTEGER :: nLats 			! number of points in latitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: latGridMap      ! dimensioned (nLats)

     INTEGER :: nLons			! number of points in longitude grid

     REAL(r8), DIMENSION(maxGridPoints) :: longGrid        ! dimensioned (nLons)

     INTEGER :: nFreqs			! number of frequencies in range

     REAL(r8), DIMENSION(maxGridPoints) :: rangFrequency ! dimensioned (nFreqs)

     REAL(r8), DIMENSION(2) :: l3presLvl	! min, max pressure levels

     INTEGER, DIMENSION(2) :: rangWavenumber	! min, max of wave # range

     CHARACTER (LEN=3) :: mode				! asc/des/com/all

     ! Output section

     CHARACTER (LEN=FileNameLen) :: mcf			! MCF file name in CF

     CHARACTER (LEN=FileNameLen) :: fileTemplate	! file name from CF

     INTEGER :: bpFlag				! bypass the PCF file name

     INTEGER :: nGrids				! number of possible grids

     CHARACTER (LEN=GridNameLen), DIMENSION(maxNumGrids) :: quantities
	! names of possible grids, dimensioned (nGrids)

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

!----------------------------------------------------------------------
   SUBROUTINE FillL3CF (cf, pcfL2Ver, pcfL3Ver, l3Days, l3Window, l3cf)
!----------------------------------------------------------------------
! Brief description of subroutine
! This subroutine checks the parser output and fills L3CFProd_T.

! Arguments

      TYPE( Mlscf_T ), INTENT(IN) :: cf

      CHARACTER (LEN=*), INTENT(IN) :: pcfL2Ver, pcfL3Ver

      CHARACTER (LEN=*), INTENT(IN) :: l3Days(:)

      INTEGER, INTENT(IN) :: l3Window

      TYPE( L3CFProd_T ), POINTER  :: l3cf(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: pgs_td_utcToTAI

! Variables

      CHARACTER (LEN=8) :: hhmmss
      CHARACTER (LEN=10) :: l2Ver, l3Ver, label
      CHARACTER (LEN=480) :: msr
      CHARACTER (LEN=CCSDSB_LEN) :: timeB(maxWindow)

      INTEGER :: err, i, ig, iq, iGlob, iLab, iMap, iOut, iVer, indx, j
      INTEGER :: numProds, returnStatus

      REAL(r8) :: start, end, delta

! Find the section indices in the CF

      iGlob = LinearSearchStringArray(cf%Sections%MlscfSectionName, &
                                      'GlobalSettings')
      iMap = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'DailyMap')
      iOut = LinearSearchStringArray(cf%Sections%MlscfSectionName, 'Output')

! Check that the version numbers given in the CF will match the PCF file names

      iVer = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword,'L2Ver')
      IF (iVer == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                    &the CF for L2Ver.')
      l2Ver = cf%Sections(iGlob)%Cells(iVer)%CharValue
      IF (l2Ver /= pcfL2Ver) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input &
                                      &versions in the CF and PCF differ.')

      iVer = LinearSearchStringArray(cf%Sections(iGlob)%Cells%Keyword, &
                                     'OutputVersionString')
      IF (iVer == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'No entry in &
                                     &the CF for OutputVersionString.')
      l3Ver = cf%Sections(iGlob)%Cells(iVer)%CharValue
      IF (l3Ver /= pcfL3Ver) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                  'Output versions in the CF and PCF differ.')

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

! Find/save the interpolation method

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'intpMethod')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                     &keyword INTPMETHOD in the DailyMap section of the l3cf.')
       l3cf(i)%intpMethod = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

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

! Find the boundaries of the frequency range

       indx = LinearSearchStringArray( &
                   cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'rangFrequency')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                  &keyword RANGFREQUENCY in the DailyMap section of the l3cf.')
       start = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       end = cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find the increment of the frequency range

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'dF')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                             &keyword DF in the DailyMap section of the l3cf.')
       delta = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue

! Calculate the array containing the frequency range

       CALL CalculateArray(start, end, delta, l3cf(i)%rangFrequency, &
                           l3cf(i)%nFreqs)

! Find/save the min & max for the pressure levels, range of wavenumbers

       indx = LinearSearchStringArray( &
                     cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'l3presLvl')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                      &keyword L3PRESLVL in the DailyMap section of the l3cf.')
       l3cf(i)%l3presLvl(1) = cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%l3presLvl(2) = &
                      cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

       indx = LinearSearchStringArray( &
               cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'rangWavenumber')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                 &keyword RANGWAVENUMBER in the DailyMap section of the l3cf.')
       l3cf(i)%rangWavenumber(1) = &
                             cf%Sections(iMap)%Entries(i)%Cells(indx)%RealValue
       l3cf(i)%rangWavenumber(2) = &
                       cf%Sections(iMap)%Entries(i)%Cells(indx)%RangeUpperBound

! Find/save mode

       indx = LinearSearchStringArray( &
                           cf%Sections(iMap)%Entries(i)%Cells%Keyword, 'mode')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                           &keyword MODE in the DailyMap section of the l3cf.')
       l3cf(i)%mode = cf%Sections(iMap)%Entries(i)%Cells(indx)%CharValue

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

! Find/save the MCF file name

       indx = LinearSearchStringArray( &
                          cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'mcf')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                              &keyword MCF in the Output section of the l3cf.')
       l3cf(i)%mcf = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

! Find the file template, bypass flag from Output; save them in L3CFProd_T

       indx = LinearSearchStringArray( &
                         cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'file')
       IF (indx == 0) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Missing &
                             &keyword FILE in the Output section of the l3cf.') 
       l3cf(i)%fileTemplate = &
                          cf%Sections(iOut)%Entries(iLab)%Cells(indx)%CharValue

       indx = LinearSearchStringArray( &
                    cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'bypassPCF')
       IF (indx /= 0) THEN
         l3cf(i)%bpFlag = cf%Sections(iOut)%Entries(iLab)%Cells(indx)%RealValue
       ELSE
         l3cf(i)%bpFlag = 0
       ENDIF

! Find/save the number of possible grids & their names

       iq = LinearSearchStringArray( &
                   cf%Sections(iOut)%Entries(iLab)%Cells%Keyword, 'quantities')

       DO j = 1, maxNumGrids
          ig = iq + (j-1)
          l3cf(i)%quantities(j) = &
                            cf%Sections(iOut)%Entries(iLab)%Cells(ig)%CharValue
          IF (cf%Sections(iOut)%Entries(iLab)%Cells(ig)%More == 0) EXIT
       ENDDO

       l3cf(i)%nGrids = j
      
    ENDDO

!-------------------------
   END SUBROUTINE FillL3CF
!-------------------------

!---------------------------------------
   SUBROUTINE ReadParseMLSCF (pcfId, cf)
!---------------------------------------

! Brief description of subroutine
! This subroutine simulates the running of the parser in the delivery
! environment.

! Arguments

      INTEGER, INTENT(IN) :: pcfId

      TYPE( Mlscf_T ), INTENT(OUT) :: cf

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: err

      cf%NoSections = 3

      ALLOCATE (cf%Sections(cf%NoSections), STAT=err )
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' Section structures for L3CF.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      cf%Sections(1)%MlscfSectionName = 'GlobalSettings'
      cf%Sections(1)%NoSectionDefs = 2
      cf%Sections(1)%Cells(1)%Keyword = 'L2Ver'
      cf%Sections(1)%Cells(1)%CharValue = 'v0-05'
      cf%Sections(1)%Cells(2)%Keyword = 'OutputVersionString'
      cf%Sections(1)%Cells(2)%CharValue = 'v0-05'

      cf%Sections(2)%MlscfSectionName = 'DailyMap'
      cf%Sections(2)%NoSectionEntries = 10
      cf%Sections(2)%Entries(1)%MlscfEntryName = 'mapSpec'
      cf%Sections(2)%Entries(1)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(1)%Cells(1)%CharValue = 'ClO'
      cf%Sections(2)%Entries(1)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(1)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(1)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(1)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(1)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(1)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(1)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(1)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(1)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(1)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(1)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(1)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(1)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(1)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(1)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(1)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(1)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(1)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(1)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(1)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(1)%Cells(10)%RealValue = 1.0
      cf%Sections(2)%Entries(1)%Cells(10)%RangeUpperBound =100.0
      cf%Sections(2)%Entries(1)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(1)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(1)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(1)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(1)%Cells(12)%CharValue = 'all'
      cf%Sections(2)%Entries(1)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(1)%Cells(13)%CharValue = 'CLOL3DM'
      cf%Sections(2)%Entries(2)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(2)%Cells(1)%CharValue = 'GeopotentialHeight'
      cf%Sections(2)%Entries(2)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(2)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(2)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(2)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(2)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(2)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(2)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(2)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(2)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(2)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(2)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(2)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(2)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(2)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(2)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(2)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(2)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(2)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(2)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(2)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(2)%Cells(10)%RealValue = 0.01
      cf%Sections(2)%Entries(2)%Cells(10)%RangeUpperBound =464.0
      cf%Sections(2)%Entries(2)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(2)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(2)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(2)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(2)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(2)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(2)%Cells(13)%CharValue = 'ZL3DM'
      cf%Sections(2)%Entries(3)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(3)%Cells(1)%CharValue = 'H2O'
      cf%Sections(2)%Entries(3)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(3)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(3)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(3)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(3)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(3)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(3)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(3)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(3)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(3)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(3)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(3)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(3)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(3)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(3)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(3)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(3)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(3)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(3)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(3)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(3)%Cells(10)%RealValue = 0.01
      cf%Sections(2)%Entries(3)%Cells(10)%RangeUpperBound =464.0
      cf%Sections(2)%Entries(3)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(3)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(3)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(3)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(3)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(3)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(3)%Cells(13)%CharValue = 'H2OL3DM'
      cf%Sections(2)%Entries(4)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(4)%Cells(1)%CharValue = 'HCl'
      cf%Sections(2)%Entries(4)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(4)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(4)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(4)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(4)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(4)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(4)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(4)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(4)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(4)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(4)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(4)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(4)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(4)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(4)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(4)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(4)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(4)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(4)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(4)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(4)%Cells(10)%RealValue = 1.0
      cf%Sections(2)%Entries(4)%Cells(10)%RangeUpperBound =100.0
      cf%Sections(2)%Entries(4)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(4)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(4)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(4)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(4)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(4)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(4)%Cells(13)%CharValue = 'HCLL3DM'
      cf%Sections(2)%Entries(5)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(5)%Cells(1)%CharValue = 'HCN'
      cf%Sections(2)%Entries(5)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(5)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(5)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(5)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(5)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(5)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(5)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(5)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(5)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(5)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(5)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(5)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(5)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(5)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(5)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(5)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(5)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(5)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(5)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(5)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(5)%Cells(10)%RealValue = 22.0
      cf%Sections(2)%Entries(5)%Cells(10)%RangeUpperBound =100.0
      cf%Sections(2)%Entries(5)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(5)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(5)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(5)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(5)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(5)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(5)%Cells(13)%CharValue = 'HCNL3DM'
      cf%Sections(2)%Entries(6)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(6)%Cells(1)%CharValue = 'HNO3'
      cf%Sections(2)%Entries(6)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(6)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(6)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(6)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(6)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(6)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(6)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(6)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(6)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(6)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(6)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(6)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(6)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(6)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(6)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(6)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(6)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(6)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(6)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(6)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(6)%Cells(10)%RealValue = 1.0
      cf%Sections(2)%Entries(6)%Cells(10)%RangeUpperBound =100.0
      cf%Sections(2)%Entries(6)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(6)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(6)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(6)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(6)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(6)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(6)%Cells(13)%CharValue = 'HNO3L3DM'
      cf%Sections(2)%Entries(7)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(7)%Cells(1)%CharValue = 'N2O'
      cf%Sections(2)%Entries(7)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(7)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(7)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(7)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(7)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(7)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(7)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(7)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(7)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(7)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(7)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(7)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(7)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(7)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(7)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(7)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(7)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(7)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(7)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(7)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(7)%Cells(10)%RealValue = 1.0
      cf%Sections(2)%Entries(7)%Cells(10)%RangeUpperBound =215.0
      cf%Sections(2)%Entries(7)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(7)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(7)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(7)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(7)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(7)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(7)%Cells(13)%CharValue = 'N2OL3DM'
      cf%Sections(2)%Entries(8)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(8)%Cells(1)%CharValue = 'O3'
      cf%Sections(2)%Entries(8)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(8)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(8)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(8)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(8)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(8)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(8)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(8)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(8)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(8)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(8)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(8)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(8)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(8)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(8)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(8)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(8)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(8)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(8)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(8)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(8)%Cells(10)%RealValue = 0.01
      cf%Sections(2)%Entries(8)%Cells(10)%RangeUpperBound =316.0
      cf%Sections(2)%Entries(8)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(8)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(8)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(8)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(8)%Cells(12)%CharValue = 'all'
      cf%Sections(2)%Entries(8)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(8)%Cells(13)%CharValue = 'O3L3DM'
      cf%Sections(2)%Entries(9)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(9)%Cells(1)%CharValue = 'OH'
      cf%Sections(2)%Entries(9)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(9)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(9)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(9)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(9)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(9)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(9)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(9)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(9)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(9)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(9)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(9)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(9)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(9)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(9)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(9)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(9)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(9)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(9)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(9)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(9)%Cells(10)%RealValue = 0.1
      cf%Sections(2)%Entries(9)%Cells(10)%RangeUpperBound =10.0
      cf%Sections(2)%Entries(9)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(9)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(9)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(9)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(9)%Cells(12)%CharValue = 'all'
      cf%Sections(2)%Entries(9)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(9)%Cells(13)%CharValue = 'OHL3DM'
      cf%Sections(2)%Entries(10)%Cells(1)%Keyword = 'l3prodNameD'
      cf%Sections(2)%Entries(10)%Cells(1)%CharValue = 'Temperature'
      cf%Sections(2)%Entries(10)%Cells(2)%Keyword = 'timeD'
      cf%Sections(2)%Entries(10)%Cells(2)%CharValue = '12:00:00'
      cf%Sections(2)%Entries(10)%Cells(3)%Keyword = 'intpMethod'
      cf%Sections(2)%Entries(10)%Cells(3)%CharValue = 'lin'
      cf%Sections(2)%Entries(10)%Cells(4)%Keyword = 'latGridMap'
      cf%Sections(2)%Entries(10)%Cells(4)%RealValue = -82.0
      cf%Sections(2)%Entries(10)%Cells(4)%RangeUpperBound = 82.0
      cf%Sections(2)%Entries(10)%Cells(5)%Keyword = 'dLat'
      cf%Sections(2)%Entries(10)%Cells(5)%RealValue = 2.0
      cf%Sections(2)%Entries(10)%Cells(6)%Keyword = 'longGrid'
      cf%Sections(2)%Entries(10)%Cells(6)%RealValue = -180.0
      cf%Sections(2)%Entries(10)%Cells(6)%RangeUpperBound = 180.0
      cf%Sections(2)%Entries(10)%Cells(7)%Keyword = 'dLon'
      cf%Sections(2)%Entries(10)%Cells(7)%RealValue = 4.0
      cf%Sections(2)%Entries(10)%Cells(8)%Keyword = 'rangFrequency'
      cf%Sections(2)%Entries(10)%Cells(8)%RealValue = 0.1
      cf%Sections(2)%Entries(10)%Cells(8)%RangeUpperBound =1.0
      cf%Sections(2)%Entries(10)%Cells(9)%Keyword = 'dF'
      cf%Sections(2)%Entries(10)%Cells(9)%RealValue = 0.25
      cf%Sections(2)%Entries(10)%Cells(10)%Keyword = 'l3presLvl'
      cf%Sections(2)%Entries(10)%Cells(10)%RealValue = 0.01
      cf%Sections(2)%Entries(10)%Cells(10)%RangeUpperBound =464.0
      cf%Sections(2)%Entries(10)%Cells(11)%Keyword = 'rangWavenumber'
      cf%Sections(2)%Entries(10)%Cells(11)%RealValue = 1.0
      cf%Sections(2)%Entries(10)%Cells(11)%RangeUpperBound =15.0
      cf%Sections(2)%Entries(10)%Cells(12)%Keyword = 'mode'
      cf%Sections(2)%Entries(10)%Cells(12)%CharValue = 'com'
      cf%Sections(2)%Entries(10)%Cells(13)%Keyword = 'label'
      cf%Sections(2)%Entries(10)%Cells(13)%CharValue = 'TL3DM'

      cf%Sections(3)%MlscfSectionName = 'Output'
      cf%Sections(3)%NoSectionEntries = 10
      cf%Sections(3)%Entries(1)%MlscfLabelName = 'CLOL3DM'
      cf%Sections(3)%Entries(1)%MlscfEntryName = 'output'
      cf%Sections(3)%Entries(1)%MlscfEntryNoKeys = 6
      cf%Sections(3)%Entries(1)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(1)%Cells(1)%CharValue = 'ClO'
      cf%Sections(3)%Entries(1)%Cells(1)%More = 2
      cf%Sections(3)%Entries(1)%Cells(2)%CharValue = 'ClOAscending'
      cf%Sections(3)%Entries(1)%Cells(2)%More = 1
      cf%Sections(3)%Entries(1)%Cells(3)%CharValue = 'ClODescending'
      cf%Sections(3)%Entries(1)%Cells(3)%More = 0
      cf%Sections(3)%Entries(1)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(1)%Cells(4)%CharValue = &
                                            '/duke/nakamura/L1BOA/META/MCF.clo'
      cf%Sections(3)%Entries(1)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(1)%Cells(5)%CharValue = &
                   '/duke/nakamura/L1BOA/mls_l3dm_clo_$version-$cycle_$day.dat'
!     cf%Sections(3)%Entries(1)%Cells(6)%Keyword = 'bypassPCF'
!     cf%Sections(3)%Entries(1)%Cells(6)%RealValue = 1.0
      cf%Sections(3)%Entries(2)%MlscfLabelName = 'ZL3DM'
      cf%Sections(3)%Entries(2)%MlscfEntryNoKeys = 5
      cf%Sections(3)%Entries(2)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(2)%Cells(1)%CharValue = 'GeopotentialHeight'
      cf%Sections(3)%Entries(2)%Cells(1)%More = 2
      cf%Sections(3)%Entries(2)%Cells(2)%CharValue = &
                                                  'GeopotentialHeightAscending'
      cf%Sections(3)%Entries(2)%Cells(2)%More = 1
      cf%Sections(3)%Entries(2)%Cells(3)%CharValue = &
                                                 'GeopotentialHeightDescending'
      cf%Sections(3)%Entries(2)%Cells(3)%More = 0
      cf%Sections(3)%Entries(2)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(2)%Cells(4)%CharValue = &
                                              '/duke/nakamura/L1BOA/META/MCF.z'
      cf%Sections(3)%Entries(2)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(2)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_z_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(3)%MlscfLabelName = 'H2OL3DM'
      cf%Sections(3)%Entries(3)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(3)%Cells(1)%CharValue = 'H2O'
      cf%Sections(3)%Entries(3)%Cells(1)%More = 2
      cf%Sections(3)%Entries(3)%Cells(2)%CharValue = 'H2OAscending'
      cf%Sections(3)%Entries(3)%Cells(2)%More = 1
      cf%Sections(3)%Entries(3)%Cells(3)%CharValue = 'H2ODescending'
      cf%Sections(3)%Entries(3)%Cells(3)%More = 0
      cf%Sections(3)%Entries(3)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(3)%Cells(4)%CharValue = &
                                            '/duke/nakamura/L1BOA/META/MCF.h2o'
      cf%Sections(3)%Entries(3)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(3)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_h2o_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(4)%MlscfLabelName = 'HCLL3DM'
      cf%Sections(3)%Entries(4)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(4)%Cells(1)%CharValue = 'HCl'
      cf%Sections(3)%Entries(4)%Cells(1)%More = 2
      cf%Sections(3)%Entries(4)%Cells(2)%CharValue = 'HClAscending'
      cf%Sections(3)%Entries(4)%Cells(2)%More = 1
      cf%Sections(3)%Entries(4)%Cells(3)%CharValue = 'HClDescending'
      cf%Sections(3)%Entries(4)%Cells(3)%More = 0
      cf%Sections(3)%Entries(4)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(4)%Cells(4)%CharValue = &
                                            '/duke/nakamura/L1BOA/META/MCF.hcl'
      cf%Sections(3)%Entries(4)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(4)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_hcl_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(5)%MlscfLabelName = 'HCNL3DM'
      cf%Sections(3)%Entries(5)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(5)%Cells(1)%CharValue = 'HCN'
      cf%Sections(3)%Entries(5)%Cells(1)%More = 2
      cf%Sections(3)%Entries(5)%Cells(2)%CharValue = 'HCNAscending'
      cf%Sections(3)%Entries(5)%Cells(2)%More = 1
      cf%Sections(3)%Entries(5)%Cells(3)%CharValue = 'HCNDescending'
      cf%Sections(3)%Entries(5)%Cells(3)%More = 0
      cf%Sections(3)%Entries(5)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(5)%Cells(4)%CharValue = &
                                            '/duke/nakamura/L1BOA/META/MCF.hcn'
      cf%Sections(3)%Entries(5)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(5)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_hcn_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(6)%MlscfLabelName = 'HNO3L3DM'
      cf%Sections(3)%Entries(6)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(6)%Cells(1)%CharValue = 'HNO3'
      cf%Sections(3)%Entries(6)%Cells(1)%More = 2
      cf%Sections(3)%Entries(6)%Cells(2)%CharValue = 'HNO3Ascending'
      cf%Sections(3)%Entries(6)%Cells(2)%More = 1
      cf%Sections(3)%Entries(6)%Cells(3)%CharValue = 'HNO3Descending'
      cf%Sections(3)%Entries(6)%Cells(3)%More = 0
      cf%Sections(3)%Entries(6)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(6)%Cells(4)%CharValue = &
                                          '/duke/nakamura/L1BOA/META/MCF.hno3'
      cf%Sections(3)%Entries(6)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(6)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_hno3_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(7)%MlscfLabelName = 'N2OL3DM'
      cf%Sections(3)%Entries(7)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(7)%Cells(1)%CharValue = 'N2O'
      cf%Sections(3)%Entries(7)%Cells(1)%More = 2
      cf%Sections(3)%Entries(7)%Cells(2)%CharValue = 'N2OAscending'
      cf%Sections(3)%Entries(7)%Cells(2)%More = 1
      cf%Sections(3)%Entries(7)%Cells(3)%CharValue = 'N2ODescending'
      cf%Sections(3)%Entries(7)%Cells(3)%More = 0
      cf%Sections(3)%Entries(7)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(7)%Cells(4)%CharValue = &
                                            '/duke/nakamura/L1BOA/META/MCF.n2o'
      cf%Sections(3)%Entries(7)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(7)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_n2o_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(8)%MlscfLabelName = 'O3L3DM'
      cf%Sections(3)%Entries(8)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(8)%Cells(1)%CharValue = 'O3'
      cf%Sections(3)%Entries(8)%Cells(1)%More = 2
      cf%Sections(3)%Entries(8)%Cells(2)%CharValue = 'O3Ascending'
      cf%Sections(3)%Entries(8)%Cells(2)%More = 1
      cf%Sections(3)%Entries(8)%Cells(3)%CharValue = 'O3Descending'
      cf%Sections(3)%Entries(8)%Cells(3)%More = 0
      cf%Sections(3)%Entries(8)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(8)%Cells(4)%CharValue = &
                                             '/duke/nakamura/L1BOA/META/MCF.o3'
      cf%Sections(3)%Entries(8)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(8)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_o3_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(9)%MlscfLabelName = 'OHL3DM'
      cf%Sections(3)%Entries(9)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(9)%Cells(1)%CharValue = 'OH'
      cf%Sections(3)%Entries(9)%Cells(1)%More = 2
      cf%Sections(3)%Entries(9)%Cells(2)%CharValue = 'OHAscending'
      cf%Sections(3)%Entries(9)%Cells(2)%More = 1
      cf%Sections(3)%Entries(9)%Cells(3)%CharValue = 'OHDescending'
      cf%Sections(3)%Entries(9)%Cells(3)%More = 0
      cf%Sections(3)%Entries(9)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(9)%Cells(4)%CharValue = &
                                             '/duke/nakamura/L1BOA/META/MCF.oh'
      cf%Sections(3)%Entries(9)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(9)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_oh_$version-$cycle_$day.dat'
      cf%Sections(3)%Entries(10)%MlscfLabelName = 'TL3DM'
      cf%Sections(3)%Entries(10)%Cells(1)%Keyword = 'quantities'
      cf%Sections(3)%Entries(10)%Cells(1)%CharValue = 'Temperature'
      cf%Sections(3)%Entries(10)%Cells(1)%More = 2
      cf%Sections(3)%Entries(10)%Cells(2)%CharValue = 'TemperatureAscending'
      cf%Sections(3)%Entries(10)%Cells(2)%More = 1
      cf%Sections(3)%Entries(10)%Cells(3)%CharValue = 'TemperatureDescending'
      cf%Sections(3)%Entries(10)%Cells(3)%More = 0
      cf%Sections(3)%Entries(10)%Cells(4)%Keyword = 'mcf'
      cf%Sections(3)%Entries(10)%Cells(4)%CharValue = &
                                              '/duke/nakamura/L1BOA/META/MCF.t'
      cf%Sections(3)%Entries(10)%Cells(5)%Keyword = 'file'
      cf%Sections(3)%Entries(10)%Cells(5)%CharValue = &
                  '/duke/nakamura/L1BOA/mls_l3dm_t_$version-$cycle_$day.dat'

!-------------------------------
   END SUBROUTINE ReadParseMLSCF
!-------------------------------

!==============
END MODULE L3CF
!==============

! $Log$
