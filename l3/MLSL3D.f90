! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==========================================
PROGRAM MLSL3D ! MLS Level 3 Daily software
!==========================================

  USE Allocate_Deallocate, ONLY: DEALLOCATE_TEST
  USE L2GPData, ONLY: L2GPData_T, DestroyL2GPDatabase
  USE L2Interface, ONLY: ReadL2GPProd, ReadL2GPAttribute
  USE L3CF, ONLY: L3CFDef_T, L3CFProd_T  
  USE L3DMData, ONLY: L3DMData_T, DestroyL3DMDatabase
  USE L3SPData, ONLY: L3SPData_T, DestroyL3SPDatabase
  USE MLSCF, ONLY: Mlscf_T 
  USE MLSCommon, ONLY: r8
  USE MLSFiles, ONLY: HDFVERSION_5, HDFVERSION_4
  USE MLSL3Common, ONLY: DATE_LEN, maxWindow
  USE MLSMessageModule, ONLY: MLSMSG_Info, MLSMSG_Warning, MLSMessage, & 
       & MLSMessageExit
  USE OpenInit, ONLY: PCFData_T, OpenAndInitialize
  USE OutputClose, ONLY: OutputFlags_T, OutputProd, OutputAndClose
  USE PCFHdr, ONLY: GlobalAttributes
  USE Synoptic, ONLY: DailyCoreProcessing

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
  CHARACTER (LEN=DATE_LEN) :: mis_Days(maxWindow)
  CHARACTER (LEN=1), POINTER :: anText(:) => NULL()
  
  REAL(r8), POINTER :: avgPer(:) => NULL()
  
  INTEGER :: i, l2Days, mis_l2Days, hdfVersion
  
  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2

  CALL MLSMessage (MLSMSG_Info, ModuleName, & 
       & 'EOS MLS Level 3 data processing started')
  
   GlobalAttributes%ProcessLevel = '3-daily'
  ! Fill structures with input data from the PCF and L3CF.
  
  CALL OpenAndInitialize(pcf, cf, cfProd, cfDef, anText, avgPer)
 
  ! Determine the necessary version for HDF output
  
  if (cfDef%hdfOutputVersionString .EQ. 'hdfeos5' ) then 
     
     hdfVersion = HDFVERSION_5
     
  else 
     
     hdfVersion = HDFVERSION_4
     
  endif
  
  ! For each product in the DailyMap section of the cf,
  
  DO i = 1, SIZE(cfProd)
    
     ! Read l2gp attributes

     CALL ReadL2GPAttribute (pcf%l3StartDay, pcf%l3EndDay, &
        & cfProd(i)%fileTemplate)
 
     ! Read all the l2gp data which exist in the input window for that product
     
     CALL ReadL2GPProd(cfProd(i)%l3prodNameD, cfProd(i)%fileTemplate, &
          & pcf%l2StartDay, pcf%l2EndDay, l2Days, mis_l2Days, mis_Days, l2gp)
     
     ! If no data , go on to the next product

     IF (l2Days < cfDef%minDays ) THEN
        msr = 'Insufficient data found for ' & 
             & // TRIM(cfProd(i)%l3prodNameD) // &
             & '.  Skipping CORE processing and moving on to the next product.'
        CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
        CYCLE
     ELSE
        msr = 'Input data successfully read for ' //  &
             & TRIM(cfProd(i)%l3prodNameD) // '; CORE processing started ...'
        CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
     ENDIF
     
     ! CORE processing

     CALL DailyCoreProcessing(cfDef, cfProd(i), pcf, l2Days, l2gp, avgPer, & 
          & l3sp, l3dm, dmA, dmD, l3r, residA, residD, & 
          & flags)

     ! Check the output data and place them into the appropriate files. 
     ! Write the l3dm metadata.  
     ! Perform any deallocations needed within the product loop.

     msr = 'CORE processing completed for ' // TRIM(cfProd(i)%l3prodNameD) &
          & // '; starting Output task ...'
     CALL MLSMessage (MLSMSG_Info, ModuleName, msr)
     
     CALL OutputProd(pcf, cfProd(i), anText, l3sp, l3dm, dmA, dmD, l3r, & 
          & residA, residD, flags, hdfVersion)
     
     ! Deallocate memory 

     CALL Deallocate_test(GlobalAttributes%OrbNumDays, &
	& 'GlobalAttributes%OrbNumDays', ModuleName)   
     CALL Deallocate_test(GlobalAttributes%OrbPeriodDays, &
	& 'GlobalAttributes%OrbPeriodDays', ModuleName)   

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
  
  ! Deallocate memory and close

  CALL MLSMessage (MLSMSG_Info, ModuleName, & 
       & 'EOS MLS Level 3 data processing successfully completed!')
  
  CALL MLSMessageExit(NORMAL_EXIT_STATUS)
  
  ! Detailed description of program
  ! The program is a prototype for the MLS Level 3 Daily software.

!=================
END PROGRAM MLSL3D
!=================

! $Log$
! Revision 1.13  2003/09/15 18:30:48  cvuu
! Read OrbitNumber and OrbitPeriod from L2GP files
!
! Revision 1.12  2003/03/22 02:45:34  jdone
! added HDFEOS2/HDFEOS5 functionality; HDFEOS2 is default
!
! Revision 1.11  2002/08/22 23:33:32  pwagner
! Made ready for hdf5
!
