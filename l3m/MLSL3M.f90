
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
   TYPE( L3MMData_T ) :: mm
   TYPE( L3MZData_T ) :: mz
   TYPE( Mlscf_T ) :: cf
   TYPE( PCFMData_T ) :: pcf
   TYPE( L2GPData_T ), POINTER :: l2gp(:)
   TYPE( L3CFDg_T ), POINTER :: cfDg(:)
   TYPE( L3CFMProd_T ), POINTER :: cfProd(:)
   TYPE( L3DZData_T ), POINTER :: dz(:)
   TYPE( OutputFiles_T ) :: dFiles, sFiles

   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: pcfNames(maxWindow)
   CHARACTER (LEN=1), POINTER :: anText(:)

   INTEGER :: first, i, l2Days, last, nlev, numFiles

   LOGICAL :: writel3dz

! Initializations

   sFiles%nFiles = 0
   dFiles%nFiles = 0
   sFiles%name = ''
   dFiles%name = ''
   sFiles%date = ''
   dFiles%date = ''

   writel3dz = .FALSE.
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

! For combined mode, sort the L2GP data according to the desired pressure levels.

!     CALL Sort (cfProd(i)%l3presLvl(1), cfProd(i)%l3presLvl(2), l2gp(1), first, &
!                last, nlev)

! Process & output the daily zonal means

!     CALL DailyZonalMean(pcf, cfDef, cfProd(i)%l3prodName, 'com', l2gp(1), first, &
!                         last, nlev, dz, writel3dz)
      IF (writel3dz) THEN
         CALL OutputL3DZ(cfDef%stdType, dz, sFiles)
      ELSE
         msr = TRIM(cfProd(i)%l3prodName) // ' DZ combined' // NOOUT_ERR
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      ENDIF
      CALL DestroyL3DZDatabase(dz)
      writel3dz = .FALSE.

! Process & output the monthly zonal means

!     CALL MonthlyZonalMean(cfProd(i)%l3prodName, cfDef, l2gp, 'com', first, last, &
!                           nlev, mz)
      CALL OutputL3MZ(pcf%zsName, mz, flag%createZS)
      CALL DeallocateL3MZ(mz)

! If required for this mode, process & output the monthly map

      IF ( (cfProd(i)%mode == 'com') .OR. (cfProd(i)%mode == 'all') ) THEN
!        CALL MonthlyMap (cfProd(i)%l3prodName, cfProd(i)%latGridMap, &
!                         cfProd(i)%nLats, cfProd(i)%longGrid, cfProd(i)%nLons, &
!                         l2gp, 'com', first, last, nlev, mm)
         CALL OutputMMGrids(pcf%msName, mm, flag%createMS)
         CALL DeallocateL3MM(mm)
      ENDIF

! Ascending -- sort, dz, mz, mm

!     CALL Sort (cfProd(i)%ascPresLvl(1), cfProd(i)%ascPresLvl(2), l2gp(1), first, &
!                last, nlev)

!     CALL DailyZonalMean(pcf, cfDef, cfProd(i)%l3prodName, 'asc', l2gp(1), first, &
!                         last, nlev, dz, writel3dz)
      IF (writel3dz) THEN
         CALL OutputL3DZ(cfDef%stdType, dz, sFiles)
      ELSE
         msr = TRIM(cfProd(i)%l3prodName) // ' DZ asc' // NOOUT_ERR
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      ENDIF
      CALL DestroyL3DZDatabase(dz)
      writel3dz = .FALSE.

!     CALL MonthlyZonalMean(cfProd(i)%l3prodName, cfDef, l2gp, 'asc', first, last, &
!                           nlev, mz)
      CALL OutputL3MZ(pcf%zsName, mz, flag%createZS)
      CALL DeallocateL3MZ(mz)

      IF ( INDEX(cfProd(i)%mode,'a') /= 0) THEN
         CALL MonthlyMap (cfProd(i)%l3prodName, cfProd(i)%latGridMap, &
                          cfProd(i)%nLats, cfProd(i)%longGrid, cfProd(i)%nLons, &
                          l2gp, 'asc', first, last, nlev, mm)
         CALL OutputMMGrids(pcf%msName, mm, flag%createMS)
         CALL DeallocateL3MM(mm)
      ENDIF

! Descending

!     CALL Sort (cfProd(i)%desPresLvl(1), cfProd(i)%desPresLvl(2), l2gp(1), first, &
!                last, nlev)

!     CALL DailyZonalMean(pcf, cfDef, cfProd(i)%l3prodName, 'des', l2gp(1), first, &
!                         last, nlev, dz, writel3dz)
      IF (writel3dz) THEN
         CALL OutputL3DZ(cfDef%stdType, dz, sFiles)
      ELSE
         msr = TRIM(cfProd(i)%l3prodName) // ' DZ des' // NOOUT_ERR
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      ENDIF
      CALL DestroyL3DZDatabase(dz)
      writel3dz = .FALSE.

!     CALL MonthlyZonalMean(cfProd(i)%l3prodName, cfDef, l2gp, 'des', first, last, &
!                           nlev, mz)
      CALL OutputL3MZ(pcf%zsName, mz, flag%createZS)
      CALL DeallocateL3MZ(mz)

      IF ( (INDEX(cfProd(i)%mode,'d') /= 0) .OR. (cfProd(i)%mode == 'all') )THEN
!        CALL MonthlyMap (cfProd(i)%l3prodName, cfProd(i)%latGridMap, &
!                         cfProd(i)%nLats, cfProd(i)%longGrid, cfProd(i)%nLons, &
!                         l2gp, 'des', first, last, nlev, mm)
         CALL OutputMMGrids(pcf%msName, mm, flag%createMS)
         CALL DeallocateL3MM(mm)
      ENDIF

! Deallocate the l2gp database for the product

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

      writel3dz = .FALSE.

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

! Sort the data into desired pressure levels for the product

!        CALL Sort(cfDg(i)%minPresLvl, cfDg(i)%maxPresLvl, l2gp(1), first, last, &
!                  nlev)

! Combined mode processing -- DZ

!        CALL DailyZonalMean(pcf, cfDef, cfDg(i)%prodName, 'com', l2gp(1), first, &
!                            last, nlev, dz, writel3dz)

! If data are flagged for output, write them to L3DZ Diagnostic files; keep track of
! which files have been created for later metadata annotation

         IF (writel3dz) THEN
             CALL OutputL3DZ(cfDef%dgType, dz, dFiles)
         ELSE
             msr = TRIM(cfDg(i)%prodName) // ' DZ ' // NOOUT_ERR
             CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         ENDIF

! Re-initialize for next mode

         CALL DestroyL3DZDatabase(dz)
         writel3dz = .FALSE.

! Monthly ZM

!        CALL MonthlyZonalMean(cfDg(i)%prodName, cfDef, l2gp, 'com', first, last, &
!                              nlev, mz)
         CALL OutputL3MZ(pcf%zdName, mz, flag%createZD)
         CALL DeallocateL3MZ(mz)

! MM

!        CALL MonthlyMap (cfDg(i)%prodName, cfProd(1)%latGridMap, cfProd(1)%nLats, &
!              cfDg(i)%lonGrid, cfDg(i)%nLon, l2gp, 'com', first, last, nlev, mm)
         CALL OutputMMGrids(pcf%mdName, mm, flag%createMD)
         CALL DeallocateL3MM(mm)

! Ascending processing -- DZ

!        CALL DailyZonalMean(pcf, cfDef, cfDg(i)%prodName, 'asc', l2gp(1), first, &
!                            last, nlev, dz, writel3dz)

         IF (writel3dz) THEN
            CALL OutputL3DZ(cfDef%dgType, dz, dFiles)
         ELSE
            msr = TRIM(cfDg(i)%prodName) // 'Ascending DZ' // NOOUT_ERR
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         ENDIF

         CALL DestroyL3DZDatabase(dz)
         writel3dz = .FALSE.

! Monthly ZM

!        CALL MonthlyZonalMean(cfDg(i)%prodName, cfDef, l2gp, 'asc', first, last, &
!                              nlev, mz)
         CALL OutputL3MZ(pcf%zdName, mz, flag%createZD)
         CALL DeallocateL3MZ(mz)

! Descending

!        CALL DailyZonalMean(pcf, cfDef, cfDg(i)%prodName, 'des', l2gp(1), first, &
!                            last, nlev, dz, writel3dz)

         IF (writel3dz) THEN
            CALL OutputL3DZ(cfDef%dgType, dz, dFiles)
         ELSE
            msr = TRIM(cfDg(i)%prodName) // 'Descending DZ' // NOOUT_ERR
            CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         ENDIF

         CALL DestroyL3DZDatabase(dz)
         writel3dz = .FALSE.

!        CALL MonthlyZonalMean(cfDg(i)%prodName, cfDef, l2gp, 'des', first, last, &
!                              nlev, mz)
         CALL OutputL3MZ(pcf%zdName, mz, flag%createZD)
         CALL DeallocateL3MZ(mz)

! Deallocate the L2GP database for this product

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
!
