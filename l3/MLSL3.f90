
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===================================
PROGRAM MLSL3 ! MLS Level 3 software
!===================================

   USE OpenInit
   USE MLSCF
   USE L3CF
   USE L3DMData
   USE MLSMessageModule
   USE L2Interface
   USE L2GPData
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

   TYPE( PCFData_T ) :: pcf
   TYPE( Mlscf_T ) :: cf
   TYPE( L2GPData_T ), POINTER :: l2gp(:)
   TYPE( L3CFProd_T ), POINTER :: cfProd(:)
   TYPE( L3DMData_T ), POINTER :: l3dm(:)

   CHARACTER (LEN=480) :: msr

   INTEGER :: allGrids, count, err, i, j, k, l, l2Days
   INTEGER, ALLOCATABLE :: nlev(:)

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                             &started')

! Fill structures with input data from the PCF and L3CF.

   CALL OpenAndInitialize(pcf, cf, cfProd)

! Code for the CORE processing is expected here.
 
! -----------------------------------------------------------------------------
!
! The currently agreed-upon interface between this module and the I/O shell
! is -- INPUT l2gp data will be read into an array for a PRODUCT, where each
! element of the array is for one DAY; OUTPUT l3dm data will be returned in an
! array for a DAY, where each element of the array is for one PRODUCT.
!
! A example of a call to subroutine which performs this task could be:
!
!   CALL DailyCoreProcessing(pcf, cfProd, l3dm)
!
! where the arguments are defined in the list of variables above.
!
! More code that simulates some aspects of this task follows:
!
! -----------------------------------------------------------------------------
 
! Read the l2gp data for ClO, the first product listed in the cf (and the only
! one for which simulated input files currently exist)

   CALL ReadL2GPProd(cfProd(1), pcf, l2Days, l2gp)

   IF (l2Days < pcf%minDays) THEN
      msr = 'Insufficient input data to process ' // cfProd(1)%l3prodNameD
      CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
   ENDIF

! Simulated CORE output

   allGrids = 16

   ALLOCATE( nlev(allGrids), l3dm(allGrids), STAT=err )
   IF ( err /= 0 ) THEN
      msr = MLSMSG_Allocate // ' l3dm arrays.'
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
   ENDIF

   l3dm(1:3)%name = cfProd(1)%quantities(1:3)
   l3dm(4)%name = cfProd(2)%quantities(1)
   l3dm(5)%name = cfProd(3)%quantities(1)
   l3dm(6)%name = cfProd(4)%quantities(1)
   l3dm(7)%name = cfProd(5)%quantities(1)
   l3dm(8)%name = cfProd(6)%quantities(1)
   l3dm(9)%name = cfProd(7)%quantities(1)
   l3dm(10:12)%name = cfProd(8)%quantities(1:3)
   l3dm(13:15)%name = cfProd(9)%quantities(1:3)
   l3dm(16)%name = cfProd(10)%quantities(1)

   l3dm(1:3)%time = cfProd(1)%timeD(1)
   l3dm(4)%time = cfProd(2)%timeD(1)
   l3dm(5)%time = cfProd(3)%timeD(1)
   l3dm(6)%time = cfProd(4)%timeD(1)
   l3dm(7)%time = cfProd(5)%timeD(1)
   l3dm(8)%time = cfProd(6)%timeD(1)
   l3dm(9)%time = cfProd(7)%timeD(1)
   l3dm(10:12)%time = cfProd(8)%timeD(1)
   l3dm(13:15)%time = cfProd(9)%timeD(1)
   l3dm(16)%time = cfProd(10)%timeD(1)

   nlev(1:3) = 24
   nlev(4:5) = 56
   nlev(6) = 36
   nlev(7) = 12
   nlev(8) = 24
   nlev(9) = 28
   nlev(10) = 36
   nlev(11:12) = 30
   nlev(13:15) = 24
   nlev(16) = 56

   DO i = 1, 3
      CALL AllocateL3DM( nlev(i), cfProd(1)%nLats, cfProd(1)%nLons, l3dm(i) )
      l3dm(i)%latitude = cfProd(1)%latGridMap(1:cfProd(1)%nLats)
      l3dm(i)%longitude = cfProd(1)%longGrid(1:cfProd(1)%nLons)
   ENDDO

   DO i = 4, 9
      CALL AllocateL3DM(nlev(i), cfProd(i-2)%nLats, cfProd(i-2)%nLons, l3dm(i))
      l3dm(i)%latitude = cfProd(i-2)%latGridMap(1:cfProd(i-2)%nLats)
      l3dm(i)%longitude = cfProd(i-2)%longGrid(1:cfProd(i-2)%nLons)
   ENDDO

   DO i = 10, 12
      CALL AllocateL3DM( nlev(i), cfProd(8)%nLats, cfProd(8)%nLons, l3dm(i) )
      l3dm(i)%latitude = cfProd(8)%latGridMap(1:cfProd(8)%nLats)
      l3dm(i)%longitude = cfProd(8)%longGrid(1:cfProd(8)%nLons)
   ENDDO

   DO i = 13, 15
      CALL AllocateL3DM( nlev(i), cfProd(9)%nLats, cfProd(9)%nLons, l3dm(i) )
      l3dm(i)%latitude = cfProd(9)%latGridMap(1:cfProd(9)%nLats)
      l3dm(i)%longitude = cfProd(9)%longGrid(1:cfProd(9)%nLons)
   ENDDO

   CALL AllocateL3DM( nlev(16), cfProd(10)%nLats, cfProd(10)%nLons, l3dm(16) )
   l3dm(16)%latitude = cfProd(10)%latGridMap(1:cfProd(10)%nLats)
   l3dm(16)%longitude = cfProd(10)%longGrid(1:cfProd(10)%nLons)

   DO i = 1, allGrids

      DO j = 1, l3dm(i)%nLevels

         l3dm(i)%pressure(j) = l2gp(1)%pressures(j)

         count = 1
         DO k = 1, l3dm(i)%nLons
            DO l = 1, l3dm(i)%nLats
               IF (count > SIZE(l2gp(1)%l2gpValue,3)) THEN
                  l3dm(i)%l3dmValue(j,l,k) = 0.0
               ELSE
                  l3dm(i)%l3dmValue(j,l,k) = l2gp(1)%l2gpValue(1,j,count)
               ENDIF
               count = l + 1
            ENDDO
         ENDDO

      ENDDO

      l3dm(i)%l3dmPrecision = -1.0

   ENDDO

! Check the output data and place them to the appropriate files. Write the
! metadata.  Perform any outstanding deallocations.

   CALL OutputAndClose(pcf, cf, cfProd, l3dm)

! Deallocate the l2gp database, when finished

   CALL DestroyL2GPDatabase(l2gp)

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                                     &successfully completed!')

! Detailed description of program
! The program is a prototype for the MLS Level 3 software.

!================
END PROGRAM MLSL3
!================

! $Log$
! Revision 1.4  2000/11/22 17:25:23  nakamura
! Required to go back to 1.2 -- same as 1.1
!
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
