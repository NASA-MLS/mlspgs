
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================
PROGRAM MLSL3M ! MLS Level 3 Monthly subprogram
!==============================================

   USE L2GPData, ONLY: L2GPData_T, DestroyL2GPDatabase
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
   TYPE( L3MZData_T ) :: mzA, mzD
   TYPE( Mlscf_T ) :: cf
   TYPE( PCFMData_T ) :: pcf
   TYPE( L2GPData_T ), POINTER :: l2gp(:)
   TYPE( L3CFMProd_T ), POINTER :: cfStd(:), cfDg(:)
   TYPE( L3DZData_T ), POINTER :: dzA(:), dzD(:)
   TYPE( OutputFiles_T ) :: dFiles, sFiles

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: pcfNames(maxWindow)
   CHARACTER (LEN=1), POINTER :: anText(:)

   CHARACTER (LEN=DATE_LEN) :: mis_Days(maxWindow)

   INTEGER :: i, l2Days, numFiles, mis_l2Days

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

   CALL OpenMON(pcf, cf, cfStd, cfDg, cfDef, anText)

! For each Standard product requested in the cf,

   !DO i = 1, SIZE(cfStd)
      DO i = 1, 1 
	print *, i,  SIZE(cfStd)

! Read all available data in the input window

      CALL ReadL2GPProd(cfStd(i)%l3prodName, cfStd(i)%fileTemplate, &
                        pcf%startDay, pcf%endDay, l2Days, mis_l2Days, mis_Days, l2gp)

! If no data found, go on to the next product

      IF (l2Days < 1) THEN
         msr = 'No data found for ' // TRIM(cfStd(i)%l3prodName) // &
               '.  Skipping Monthly processing and moving on to the next product.'
         CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
         CYCLE
      ELSE
         msr = 'Input data successfully read; begin processing data for ' // &
               cfStd(i)%l3prodName
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
      ENDIF

! Core processing for Standard products

	print *, 'Before MonthlyCoreProcessing'
 
      CALL MonthlyCoreProcessing(cfStd(i), pcf, cfDef, l2Days, l2gp, mm, mmA, mmD, &
                                 mzA, mzD, dzA, dzD, mis_l2Days, mis_Days)

      msr = 'CORE processing completed for ' // TRIM(cfStd(i)%l3prodName) &
            // '; starting Output task ...'
      CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

	print *, 'After MonthlyCoreProcessing'
! Deallocate the L2GP database

      CALL DestroyL2GPDatabase(l2gp)

! Output and Close for the product

      CALL OutputStd(pcf, cfDef%stdType, cfStd(i)%mode, dzA, dzD, mzA, mzD, &
                     mm, mmA, mmD, sFiles, flag)

   ENDDO

! Begin Diagnostic processing

   msr = 'Standard product loop completed; beginning Diagnostic product loop ...'
   CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

! Get the names of any L2GP Diagnostics files from the PCF

   CALL GetL2GPfromPCF(mlspcf_l2dg_start, mlspcf_l2dg_end, cfDef%l2dgType, &
                       pcf%startDay, pcf%EndDay, numFiles, pcfNames, mis_l2Days, mis_Days)

! If no files are found, skip to Output/Close

   IF (numFiles < 1) THEN

      CALL MLSMessage (MLSMSG_Warning, ModuleName, 'No L2GP Diagnostics files &
                      &found in the PCF for the given day range.')

   ELSE

! For each product in the Diagnostic section of the cf,

	print *, ''

      DO i = 6, SIZE(cfDg)
	print *, 'Diagnostic=', i,  SIZE(cfDg)

! Read all the l2gp data which exist in the input window for that product

	print *, 'ReadL2DGData 1'
         CALL ReadL2DGData(cfDg(i)%l3prodName, numFiles, pcfNames, l2Days, l2gp)
	print *, 'ReadL2DGData 2'

! If insufficient data found, go on to the next product

         IF (l2Days < 1) THEN
            msr = 'No data found for ' // TRIM(cfDg(i)%l3prodName) &
               // '.  Skipping CORE processing and moving on to the next product.'
            CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
            CYCLE
         ELSE
            msr = 'Input data successfully read for ' //  &
                   TRIM(cfDg(i)%l3prodName) // '; CORE processing started ...'
            CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
         ENDIF

! Monthly Core processing

         CALL MonthlyCoreProcessing(cfDg(i), pcf, cfDef, l2Days, l2gp, &
                                    mm, mmA, mmD, mzA, mzD, dzA, dzD, mis_l2Days, mis_Days)

         msr = 'CORE processing completed for ' // TRIM(cfDg(i)%l3prodName) &
               // '; starting Output task ...'
         CALL MLSMessage (MLSMSG_Info, ModuleName, msr)

! Deallocate the L2GP database

         CALL DestroyL2GPDatabase(l2gp)

! Output and Close for the product

         CALL OutputDg(pcf, cfDef%dgType, dzA, dzD, mzA, mzD, mm, mmA, mmD, &
                       dFiles, flag)

      ENDDO

   ENDIF

! Perform final Output & Close tasks outside of the product processing loop.

   CALL MLSMessage (MLSMSG_Info, ModuleName, &
            'Product processing completed; beginning Output/Close task ... ')

   CALL OutputMON(sFiles, dFiles, flag, pcf, cfStd, cfDg, cf, anText)

   CALL MLSMessage (MLSMSG_Info, ModuleName, &
    'EOS MLS Level 3 Monthly data processing successfully completed!')

! Detailed description of program
! The program is a prototype for the MLS Level 3 Monthly subprogram.

!=================
END PROGRAM MLSL3M
!=================

! $Log$
! Revision 1.8  2001/12/12 17:48:03  nakamura
! Removed unused dz & mz arguments.
!
! Revision 1.7  2001/09/26 19:48:28  nakamura
! Removed com ZM output; added cfDg deallocate.
!
! Revision 1.6  2001/09/06 18:48:39  nakamura
! Modified so that dg products are processed with the same Core as std prods.
!
! Revision 1.5  2001/08/08 19:27:47  nakamura
! Added cfDef to MonthlyCoreProcessing args; added comments delineating CORE from I/O.
!
! Revision 1.4  2001/08/02 23:35:52  ybj
! Fix interface call to MonthlyCoreProcessing
!
! Revision 1.3  2001/08/01 18:28:56  nakamura
! I/O-Core interface switched back to Daily paradigm.
!
! Revision 1.2  2001/07/20 19:30:29  nakamura
! Commented out sample CALLs to Core routines.
!
! Revision 1.1  2001/07/18 15:42:48  nakamura
! MLS Level 3 Monthly program.
!
