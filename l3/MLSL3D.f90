
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==========================================
PROGRAM MLSL3D ! MLS Level 3 Daily software
!==========================================

   USE L2GPData, ONLY: L2GPData_T, DestroyL2GPDatabase
   USE L2Interface
   USE L3CF
   USE L3DMData
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
! This is the MLS Level 3 Daily program.

! Parameters

! Functions

! Variables

   TYPE( Mlscf_T ) :: cf
   TYPE( L3CFDef_T ) :: cfDef
   TYPE( OutputFlags_T ) :: flags
   TYPE( PCFData_T ) :: pcf
   TYPE( L2GPData_T ), POINTER :: l2gp(:), l3r(:), residA(:), residD(:)
   TYPE( L3CFProd_T ), POINTER :: cfProd(:)
   TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
   TYPE( L3SPData_T ), POINTER :: l3sp(:)

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=1), POINTER :: anText(:)
   CHARACTER (LEN=DATE_LEN) :: mis_Days(maxWindow)

   INTEGER :: i, l2Days, mis_l2Days

   REAL(r8), POINTER :: avgPer(:)
   integer, parameter :: NORMAL_EXIT_STATUS = 2

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing started')

! Fill structures with input data from the PCF and L3CF.

   CALL OpenAndInitialize(pcf, cf, cfProd, cfDef, anText, avgPer)

! For each product in the DailyMap section of the cf,

   DO i = 1, SIZE(cfProd)
 
! Read all the l2gp data which exist in the input window for that product

      CALL ReadL2GPProd(cfProd(i)%l3prodNameD, cfProd(i)%fileTemplate, &
                        pcf%l2StartDay, pcf%l2EndDay, l2Days, mis_l2Days, mis_Days, l2gp)

! If no data , go on to the next product

      IF (l2Days < cfDef%minDays ) THEN
         msr = 'Insufficient data found for ' // TRIM(cfProd(i)%l3prodNameD) // &
               '.  Skipping CORE processing and moving on to the next product.'
         CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
         CYCLE
      ELSE
         msr = 'Input data successfully read for ' //  &
                TRIM(cfProd(i)%l3prodNameD) // '; CORE processing started ...'
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
      ENDIF

! CORE processing

      CALL DailyCoreProcessing(cfDef, cfProd(i), pcf, l2Days, l2gp, avgPer, l3sp, &
                               l3dm, dmA, dmD, l3r, residA, residD, mis_l2Days, mis_Days, flags)

! Check the output data and place them into the appropriate files.  Write the
! l3dm metadata.  Perform any deallocations needed within the product loop.

      msr = 'CORE processing completed for ' // TRIM(cfProd(i)%l3prodNameD) &
            // '; starting Output task ...'
      CALL MLSMessage (MLSMSG_Info, ModuleName, msr)


      CALL OutputProd(pcf, cfProd(i), anText, l3sp, l3dm, dmA, dmD, l3r, residA, residD, flags)


! Deallocate the databases passed between CORE & the I/O shell

      CALL DestroyL2GPDatabase(l3r)
      CALL DestroyL2GPDatabase(residA)
      CALL DestroyL2GPDatabase(residD)
      CALL DestroyL3DMDatabase(l3dm)
      CALL DestroyL3DMDatabase(dmA)
      CALL DestroyL3DMDatabase(dmD)
      CALL DestroyL3SPDatabase(l3sp)
      CALL DestroyL2GPDatabase(l2gp)


   ENDDO

   msr = 'Product loop processing completed; beginning Output/Close task ... '
   CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

! Perform final Output & Close tasks outside of the product processing loop.


   CALL OutputAndClose(cf, pcf, cfProd, avgPer, anText)

   CALL MLSMessage (MLSMSG_Info, ModuleName, 'EOS MLS Level 3 data processing successfully completed!')
   call MLSMessageExit(NORMAL_EXIT_STATUS)


! Detailed description of program
! The program is a prototype for the MLS Level 3 Daily software.

!=================
END PROGRAM MLSL3D
!=================

! $Log $
!
