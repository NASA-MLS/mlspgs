
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===================================
PROGRAM MLSL3 ! MLS Level 3 software
!===================================

   USE L2GPData
   USE L2Interface
   USE L3CF
   USE L3DMData
   USE L3SPData
   USE MLSCF
   USE MLSMessageModule
   USE OpenInit
   USE OutputClose

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This is the MLS Level 3 program.

! Parameters

! Functions

! Variables

   TYPE( Mlscf_T ) :: cf
   TYPE( PCFData_T ) :: pcf
   TYPE( OutputFlags_T ) :: flags
   TYPE( L2GPData_T ), POINTER :: l2gp(:)
   TYPE( L2GPData_T ), POINTER :: l3r(:), residA(:), residD(:)
   TYPE( L3CFProd_T ), POINTER :: cfProd(:)
   TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
   TYPE( L3SPData_T ), POINTER :: l3sp(:)

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: logType
   CHARACTER (LEN=1), POINTER :: anText(:)

   INTEGER :: count, err, i, j, k, l, l2Days, nlev, nf, nwv, numDays, numSwaths
   INTEGER :: rDays

   REAL(r8), POINTER :: avgPer(:)

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                             &started')

! Fill structures with input data from the PCF and L3CF.

   CALL OpenAndInitialize(pcf, cf, cfProd, logType, anText, avgPer)

! For each product in the DailyMap section of the cf,

   DO i = 1, SIZE(cfProd)
 
! Read all the l2gp data which exists in the input window for that product

      CALL ReadL2GPProd(cfProd(i), pcf%l2StartDay, pcf%l2EndDay, l2Days, l2gp)

! If insufficient data found, go on to the next product

      IF (l2Days < pcf%minDays) THEN
         msr = 'Skipping CORE processing for ' // TRIM(cfProd(i)%l3prodNameD) &
            // ' and moving on to the next product.'
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
         CYCLE
      ENDIF

! -----------------------------------------------------------------------------
!
! CORE processing --
!
! A example of a call to subroutine which performs this task could be:
!
!   CALL DailyCoreProcessing(pcf, cfProd, avgPer, l3sp, l3dm, dmA, dmD, l3r, &
!                            residA, residD, flags)
!
! where the arguments are defined in the list of variables above.
!
! -----------------------------------------------------------------------------

! Simulated CORE output for one product

      nlev = 24
      nwv = 10
      nf = 15

      IF (cfProd(i)%mode == 'all') THEN
         numSwaths = 3
      ELSE
         numSwaths = 1
      ENDIF

      ALLOCATE( l3sp(numSwaths), STAT=err )
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3sp array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      l3sp(1)%name = cfProd(i)%l3prodNameD
      l3sp(2)%name = TRIM(cfProd(i)%l3prodNameD) // 'Ascending'
      l3sp(3)%name = TRIM(cfProd(i)%l3prodNameD) // 'Descending'

      l3sp%startTime = l2gp(1)%time(1)
      l3sp%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)

      DO j = 1, numSwaths

         CALL AllocateL3SP( nlev, cfProd(i)%nLats, nwv, nf, l3sp(j) )
         l3sp(j)%pressure = l2gp(1)%pressures(:nlev)
         l3sp(j)%latitude = cfProd(i)%latGridMap(:l3sp(j)%nLats)

         DO k = 1, l3sp(j)%nWaveNum
            l3sp(j)%waveNumber(k) = k * 1.0
         ENDDO

         DO k = 1, l3sp(j)%nFreqs
            l3sp(j)%frequency(k) = k * 1.0
         ENDDO

         l3sp(j)%l3spRelValue = 0.0
         l3sp(j)%l3spRelPrecision = -1.0
         l3sp(j)%l3spImgValue = 0.0
         l3sp(j)%l3spImgPrecision = -1.0

      ENDDO

      numDays = cfProd(i)%nDays
      ALLOCATE( l3dm(numDays), dmA(numDays), dmD(numDays), l3r(numDays), &
                residA(numDays), residD(numDays), STAT=err )
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l3dm, l3r arrays.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      l3dm%name = cfProd(i)%l3prodNameD
      dmA%name = TRIM(cfProd(i)%l3prodNameD) // 'Ascending'
      dmD%name = TRIM(cfProd(i)%l3prodNameD) // 'Descending'

      DO j = 1, numDays

         l3dm(j)%time = cfProd(i)%timeD(j)
         dmA(j)%time = cfProd(i)%timeD(j)
         dmD(j)%time = cfProd(i)%timeD(j)

         CALL AllocateL3DM( nlev, cfProd(i)%nLats, cfProd(i)%nLons, l3dm(j) )
         CALL AllocateL3DM( nlev, cfProd(i)%nLats, cfProd(i)%nLons, dmA(j) )
         CALL AllocateL3DM( nlev, cfProd(i)%nLats, cfProd(i)%nLons, dmD(j) )

         l3dm(j)%pressure = l2gp(1)%pressures(:nlev)
         dmA(j)%pressure = l2gp(1)%pressures(:nlev)
         dmD(j)%pressure = l2gp(1)%pressures(:nlev)

         l3dm(j)%latitude = cfProd(i)%latGridMap(:l3dm(j)%nLats)
         dmA(j)%latitude = cfProd(i)%latGridMap(:dmA(j)%nLats)
         dmD(j)%latitude = cfProd(i)%latGridMap(:dmD(j)%nLats)
     
         l3dm(j)%longitude = cfProd(i)%longGrid(:l3dm(j)%nLons)
         dmA(j)%longitude = cfProd(i)%longGrid(:dmA(j)%nLons)
         dmD(j)%longitude = cfProd(i)%longGrid(:dmD(j)%nLons)

         l3dm(j)%l3dmValue = 0.0
         dmA(j)%l3dmValue = 0.0
         dmD(j)%l3dmValue = 0.0

         l3dm(j)%l3dmPrecision = -1.0
         dmA(j)%l3dmPrecision = -1.0
         dmD(j)%l3dmPrecision = -1.0

      ENDDO

!  DO i = 1, numDays

!     DO j = 1, l3dm(i)%nLevels

!        count = 1
!        DO k = 1, l3dm(i)%nLons
!           DO l = 1, l3dm(i)%nLats
!              IF (count > SIZE(l2gp(20)%l2gpValue,3)) THEN
!                 l3dm(i)%l3dmValue(j,l,k) = 0.0
!              ELSE
!                 l3dm(i)%l3dmValue(j,l,k) = l2gp(20)%l2gpValue(1,j,count)
!              ENDIF
!              count = l + 1
!           ENDDO
!        ENDDO

!     ENDDO

!  ENDDO

      CALL ReadL2GPProd(cfProd(i), pcf%l3StartDay, pcf%l3EndDay, rDays, l3r)
      CALL ReadL2GPProd(cfProd(i), pcf%l3StartDay, pcf%l3EndDay, rDays, residA)
      CALL ReadL2GPProd(cfProd(i), pcf%l3StartDay, pcf%l3EndDay, rDays, residD)

      l3r%name = TRIM(cfProd(i)%l3prodNameD) // 'Residuals'
      residA%name = TRIM(cfProd(i)%l3prodNameD) // 'AscendingResiduals'
      residD%name = TRIM(cfProd(i)%l3prodNameD) // 'DescendingResiduals'

      DO j = 1, rDays
         l3r(j)%l2gpValue = 1.0
         residA(j)%l2gpValue = 2.0
         residD(j)%l2gpValue = 3.0
      ENDDO

      flags%writel3sp    = .TRUE.
      flags%writel3dmCom = .TRUE.
      flags%writel3dmAsc = .TRUE.
      flags%writel3dmDes = .TRUE.
      flags%writel3rCom  = .TRUE.
      flags%writel3rAsc  = .TRUE.
      flags%writel3rDes  = .TRUE.

!------------------------------------------------------------------------------

! Check the output data and place them into the appropriate files.  Write the
! l3dm metadata.  Perform any deallocations needed within the product loop.

      CALL OutputAndClose(pcf, cfProd(i), anText, l3sp, l3dm, dmA, dmD, l3r, &
                          residA, residD, flags)

! Deallocate the l2gp database for the product

      CALL DestroyL2GPDatabase(l2gp)

   ENDDO

! Write the log file metadata

   CALL WriteMetaLog(pcf, logType)

! Final deallocations
 
   DEALLOCATE(cfProd, anText, avgPer, STAT=err)
   IF ( err /= 0 ) THEN
      msr = MLSMSG_DeAllocate // '  Open/Init quantities.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   DEALLOCATE (cf%Sections, STAT=err)
   IF ( err /= 0 ) THEN
      msr = MLSMSG_DeAllocate // '  cf section pointers.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                                     &successfully completed!')

! Detailed description of program
! The program is a prototype for the MLS Level 3 software.

!================
END PROGRAM MLSL3
!================

! $Log$
! Revision 1.6  2000/12/29 21:40:52  nakamura
! Added avgPer, more simulated data; revised argument list for ReadL2GPProd; switched to one-product/all-days paradigm.
!
! Revision 1.5  2000/11/28 17:33:06  nakamura
! Added code to put non-zero values in simulated L3DM output.
!
! Revision 1.4  2000/11/22 17:25:23  nakamura
! Required to go back to 1.2 -- same as 1.1
!
! Revision 1.3  2000/11/22 17:24:45  nakamura
! Required to go back to 1.2 -- same as 1.1
!
! Revision 1.2  2000/11/22 17:24:30  nakamura
! Required to go back to 1.2 -- same as 1.1
!
! Revision 1.1  2000/11/22 17:02:39  nakamura
! MLS Level 3 program
!
