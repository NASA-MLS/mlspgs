
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================
PROGRAM MLSL3M ! MLS Level 3 Monthly subprogram
!==============================================

   USE L2GPData
   USE L2Interface
   USE L3DZData
   USE L3MMData
   USE L3MZData
   USE MLSCF
   USE MLSCommon
   USE MLSL3Common
   USE MLSMessageModule
   USE MLSPCF3
   USE mon_L3CF
   USE mon_Open
   USE mon_Out
   USE MonthlyProcessModule

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This is the MLS Level 3 Monthly subprogram.

! Parameters

! Functions

! Variables

   TYPE( CreateFlags_T ) :: flag
   TYPE( L3CFMDef_T ) :: cfDef
   TYPE( L3MMData_T ) :: mm, mmA, mmD
   TYPE( L3MZData_T ) :: mz, mzA, mzD
   TYPE( Mlscf_T ) :: cf
   TYPE( PCFMData_T ) :: pcf
   TYPE( L2GPData_T ), POINTER :: l2gp(:)
   TYPE( L3CFDg_T ), POINTER :: cfDg(:)
   TYPE( L3CFMProd_T ), POINTER :: cfProd(:)
   TYPE( L3DZData_T ), POINTER :: dz(:), dzA(:), dzD(:)
   TYPE( OutputFiles_T ) :: dFiles, sFiles

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: pcfNames(maxWindow)
   CHARACTER (LEN=1), POINTER :: anText(:)

   INTEGER :: i, l2Days, numFiles

! Initializations

   sFiles%nFiles = 0
   dFiles%nFiles = 0
   sFiles%name = ''
   dFiles%name = ''
   sFiles%date = ''
   dFiles%date = ''

   flag%createMS = .FALSE.
   flag%createMD = .FALSE.
   flag%createZS = .FALSE.
   flag%createZD = .FALSE.

   CALL MLSMessage (MLSMSG_Info, ModuleName, &
                                 'EOS MLS Level 3 Monthly data processing started')

! Fill structures with input data from the PCF and L3CF.

   CALL OpenMON(pcf, cf, cfProd, cfDef, cfDg, anText)

! For each Standard product requested in the cf,

   DO i = 1, SIZE(cfProd)

! Read all available data in the input window

      CALL ReadL2GPProd(cfProd(i)%l3prodName, cfProd(i)%fileTemplate, &
                        pcf%startDay, pcf%endDay, l2Days, l2gp)

! If no data found, go on to the next product

      IF (l2Days < 1) THEN
         msr = 'No data found for ' // TRIM(cfProd(i)%l3prodName) // &
               '.  Skipping Monthly processing and moving on to the next product.'
         CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
         CYCLE
      ELSE
         msr = 'Input data successfully read; begin sorting data for ' // &
               cfProd(i)%l3prodName
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
      ENDIF

! Core processing for Standard products

     CALL MonthlyCoreProcessing(cfProd(i), pcf, l2Days, l2gp, mm, mmA, mmD, mz, mzA, mzD, dz, dzA, dzD)

! Output and Close for the product

      CALL OutputStd(pcf, cfDef%stdType, cfProd(i)%mode, dz, dzA, dzD, mz, mzA, mzD, &
                     mm, mmA, mmD, sFiles, flag)

! Deallocate the databases passed between CORE & the I/O shell

      CALL DeallocateL3MM(mm)
      CALL DeallocateL3MM(mmA)
      CALL DeallocateL3MM(mmD)

      CALL DeallocateL3MZ(mz)
      CALL DeallocateL3MZ(mzA)
      CALL DeallocateL3MZ(mzD)

      CALL DestroyL3DZDatabase(dz)
      CALL DestroyL3DZDatabase(dzA)
      CALL DestroyL3DZDatabase(dzD)

      CALL DestroyL2GPDatabase(l2gp)

   ENDDO

   msr = 'Standard product loop completed; beginning Diagnostic product loop ...'
   CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

   numFiles = 0

! Get the names of any L2GP Diagnostics files from the PCF

   CALL GetL2GPfromPCF(mlspcf_l2dg_start, mlspcf_l2dg_end, cfDef%l2dgType, &
                       pcf%startDay, pcf%EndDay, numFiles, pcfNames)

! If no files are found, skip to Output/Close

   IF (numFiles < 1) THEN

      CALL MLSMessage (MLSMSG_Warning, ModuleName, 'No L2GP Diagnostics files &
                      &found in the PCF for the given day range.')

   ELSE

! For each product in the Diagnostic section of the cf,

      DO i = 1, SIZE(cfDg)

! Read all the l2gp data which exist in the input window for that product

         CALL ReadL2DGData(cfDg(i)%prodName, numFiles, pcfNames, l2Days, l2gp)

! If insufficient data found, go on to the next product

         IF (l2Days < 1) THEN
            msr = 'No data found for ' // TRIM(cfDg(i)%prodName) &
               // '.  Skipping CORE processing and moving on to the next product.'
            CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
            CYCLE
         ELSE
            msr = 'Input data successfully read for ' //  &
                   TRIM(cfDg(i)%prodName) // '; CORE processing started ...'
            CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
         ENDIF

! Monthly Core processing

!        CALL MonthlyDgProcessing(pcf, cfDef, cfDg(i), latGrid, lonGrid, l2gp, mz, &
!                                 mzA, mzD, dz, dzA, dzD, mm, mmA, mmD)

! Output and Close for the product

!        CALL OutputDg(pcf, cfDef, cfDg(i), dz, dzA, dzD, mz, mzA, mzD, mm, mmA, mmD, &
!                      dFiles, flag)

! Deallocate the databases passed between CORE & the I/O shell

         CALL DeallocateL3MM(mm)
         CALL DeallocateL3MM(mmA)
         CALL DeallocateL3MM(mmD)

         CALL DeallocateL3MZ(mz)
         CALL DeallocateL3MZ(mzA)
         CALL DeallocateL3MZ(mzD)

         CALL DestroyL3DZDatabase(dz)
         CALL DestroyL3DZDatabase(dzA)
         CALL DestroyL3DZDatabase(dzD)

         CALL DestroyL2GPDatabase(l2gp)

      ENDDO

   ENDIF

! Perform final Output & Close tasks outside of the product processing loop.

   CALL MLSMessage (MLSMSG_Info, ModuleName, &
            'Product processing completed; beginning Output/Close task ... ')

   CALL OutputMON(sFiles, dFiles, flag, pcf, cfProd, cf, anText)

   CALL MLSMessage (MLSMSG_Info, ModuleName, &
    'EOS MLS Level 3 Monthly data processing successfully completed!')

! Detailed description of program
! The program is a prototype for the MLS Level 3 Monthly subprogram.

!=================
END PROGRAM MLSL3M
!=================

! $Log$
! Revision 1.3  2001/08/01 18:28:56  nakamura
! I/O-Core interface switched back to Daily paradigm.
!
! Revision 1.2  2001/07/20 19:30:29  nakamura
! Commented out sample CALLs to Core routines.
!
! Revision 1.1  2001/07/18 15:42:48  nakamura
! MLS Level 3 Monthly program.
!
!
