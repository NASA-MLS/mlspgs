
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===================================
PROGRAM MLSL3 ! MLS Level 3 software
!===================================

   USE L2GPData, ONLY: L2GPData_T, DestroyL2GPDatabase
   USE L2Interface
   USE L3CF
   USE L3DMData
   USE L3DZData
   USE L3SPData
   USE MLSCF
   USE MLSL3Common
   USE MLSMessageModule
   USE OpenInit
   USE OutputClose
   USE Synoptic

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

   TYPE( L3CFDef_T ) :: cfDef
   TYPE( Mlscf_T ) :: cf
   TYPE( OutputFiles_T ) :: dFiles, sFiles
   TYPE( OutputFlags_T ) :: flags
   TYPE( PCFData_T ) :: pcf
   TYPE( L2GPData_T ), POINTER :: l2gp(:), l3r(:), residA(:), residD(:)
   TYPE( L3CFProd_T ), POINTER :: cfProd(:)
   TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
   TYPE( L3DZData_T ), POINTER :: dzs(:), dza(:), dzd(:)
   TYPE( L3SPData_T ), POINTER :: l3sp(:)

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=1), POINTER :: anText(:)

   INTEGER :: i, l2Days

   REAL(r8), POINTER :: avgPer(:)

! Initializations

   sFiles%nFiles = 0
   dFiles%nFiles = 0
   sFiles%name = ''
   dFiles%name = ''
   sFiles%date = ''
   dFiles%date = ''

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                             &started')

! Fill structures with input data from the PCF and L3CF.

   CALL OpenAndInitialize(pcf, cf, cfProd, cfDef, anText, avgPer)

! For each product in the DailyMap section of the cf,

   DO i = 1, SIZE(cfProd)
 
! Read all the l2gp data which exists in the input window for that product

      CALL ReadL2GPProd(cfProd(i), pcf%l2StartDay, pcf%l2EndDay, l2Days, l2gp)

! If insufficient data found, go on to the next product

      IF (l2Days < 1) THEN
         msr = 'No data found for ' // TRIM(cfProd(i)%l3prodNameD) &
            // '.  Skipping CORE processing and moving on to the next product.'
         CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
         CYCLE
      ELSE
         msr = 'Input data successfully read for ' //  &
                TRIM(cfProd(i)%l3prodNameD) // '; CORE processing started ...'
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
      ENDIF

! CORE processing
 
      CALL DailyCoreProcessing(cfProd(i), pcf, l2Days, l2gp, avgPer, l3sp, &
                     l3dm, dmA, dmD, dzs, dza, dzd, l3r, residA, residD, flags)

! Check the output data and place them into the appropriate files.  Write the
! l3dm metadata.  Perform any deallocations needed within the product loop.

      msr = 'CORE processing completed for ' // TRIM(cfProd(i)%l3prodNameD) &
            // '; starting Output task ...'
      CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

      CALL OutputProd(pcf, cfProd(i), cfDef, anText, l3sp, l3dm, dmA, dmD, &
                      l3r, residA, residD, dzs, dza, dzd, flags, sFiles)

! Deallocate the databases

      CALL DestroyL2GPDatabase(l3r)
      CALL DestroyL2GPDatabase(residA)
      CALL DestroyL2GPDatabase(residD)

      CALL DestroyL3DMDatabase(l3dm)
      CALL DestroyL3DMDatabase(dmA)
      CALL DestroyL3DMDatabase(dmD)

      CALL DestroyL3SPDatabase(l3sp)

      CALL DestroyL3DZDatabase(dzs)
      CALL DestroyL3DZDatabase(dza)
      CALL DestroyL3DZDatabase(dzd)

      CALL DestroyL2GPDatabase(l2gp)

   ENDDO

   msr = 'Product loop processing completed; beginning Output/Close task ... '
   CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

! Perform final Output & Close tasks outside of the product processing loop.

   CALL OutputAndClose(cf, pcf, cfProd, cfDef, avgPer, anText, sFiles, dFiles)

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing &
                                                     &successfully completed!')

! Detailed description of program
! The program is a prototype for the MLS Level 3 software.

!================
END PROGRAM MLSL3
!================

! $Log$
! Revision 1.12  2001/04/24 19:40:23  nakamura
! Added ONLY to USE L2GPData statement.
!
! Revision 1.11  2001/04/11 18:50:45  nakamura
! Moved deallocations for pointers passed between CORE & I/O up to this level.
!
! Revision 1.10  2001/02/28 18:00:07  nakamura
! Fixed line-break in log comment.
!
! Revision 1.9  2001/02/28 17:25:01  nakamura
! Integrated CORE module; changed name of loop routine to OutputProd & added terminal OutputAndClose.
!
! Revision 1.8  2001/01/18 16:51:50  nakamura
! Moved minDays from PCF to cf.
!
! Revision 1.7  2001/01/16 17:45:38  nakamura
! Updated for new MCFs and added annotation.
!
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
