
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============
MODULE Synoptic
!===============

  Use DailyMapModule, ONLY: Init, ClearMemory, CordTransform, & 
       & FFSM, FFSMA, FFSMD, &
       & Reconstruct, Diagnostics, DataGenerate, DataGeneratePrec, &
       & CopyPrec2Data
  USE L3CF, ONLY: L3CFDef_T, L3CFProd_T
  USE L2GPData, ONLY: L2GPData_T, DestroyL2GPDatabase
  USE MLSCommon, ONLY: r8, r4
  USE MLSL3Common, ONLY: DATE_LEN, maxWindow
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_ALLOCATE, & 
       & MLSMSG_DEALLOCATE, MLSMSG_Warning
  USE Dump_0, only: DUMP

  Implicit none
  
  private
  public :: DailyCoreProcessing, ExpandArray
  PRIVATE :: ID, ModuleName, my_sortp
  
  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !----------------------------------------------------------

  ! Contents:
  
  ! Subroutines -- DailyCoreProcessing
  !                
  ! This subroutine will process one product for all the input days data
  ! Remarks:  This module contains subroutines related to MLS L3 Daily Map 
  ! Processing 

  interface my_sortp !Becuase we sometimes call with r4, other times with r8
     module procedure my_sortp_r4
     module procedure my_sortp_r8
  end interface
  
CONTAINS
  
  !-------------------------------------------------------------------------
  SUBROUTINE DailyCoreProcessing(cfDef, cfProd, pcf, l2Days, l2gp, & 
       & avgPeriod, l3sp, l3dm, dmA, dmD, l3r, residA, residD, & 
       & flags)
  !-------------------------------------------------------------------------
  USE L2Interface, ONLY: ReadL2GPProd, SetupL2GPProd
  USE L3DMData, ONLY: L3DMData_T, AllocateL3DM
  USE L3SPData, ONLY: L3SPData_T, AllocateL3SP
  USE OpenInit, ONLY: PCFData_T 
  USE OutputClose, ONLY: OutputFlags_T
  USE SDPToolkit, ONLY: PI
     
     ! Brief description of program
     ! This is the main program to run the Core processing.

     ! Parameters

     ! Variable definitions

    TYPE( L3CFDef_T ) :: cfDef
    TYPE( PCFData_T ) :: pcf
    TYPE( L3CFProd_T ) :: cfProd
    TYPE( OutputFlags_T ) :: flags
    TYPE( L2GPData_T ), POINTER :: l2gp(:), l3r(:), residA(:), residD(:), &
         & l3r_temp(:), residA_temp(:), residD_temp(:)
    TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
    TYPE( L3SPData_T ), POINTER :: l3sp(:), l3spPrec(:)
     
    CHARACTER (LEN=480) :: msr

    INTEGER :: mis_Days_temp(maxWindow)
     
    REAL (r8), DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: delTad
     
    REAL (r8), POINTER, DIMENSION(:, :, :) ::  & 
         & alons(:, :, :), alats(:, :, :), dlons(:, :, :), dlats(:, :, :), &
         & atimes(:, :, :), dtimes(:, :, :), afields(:, :, :), & 
         & dfields(:, :, :), aprec(:, :, :), dprec(:, :, :)
     
    REAL (r8), POINTER, DIMENSION(:) :: startTime(:), endTime(:), & 
         & l3Result(:), &  ! returned reconstructed result
         & sortTemp(:)     ! hold sort array temporary
     
    REAL(r8), POINTER :: avgPeriod(:)
     
    REAL(r8) :: lonD0_in, lonA0_in
    REAL :: tau0, l3ret, avg
        
    INTEGER, DIMENSION(l2gp(1)%nLevels) :: perMisPoints
    INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 

    INTEGER, POINTER, DIMENSION(:) :: pt(:), nc(:), nca(:), ncd(:)

    INTEGER ::  error, l2Days, nlev, nlev_temp, nwv, numDays, & 
         & numSwaths, rDays, pEndIndex, pStartIndex, &
         & mis_l2Days_temp, i, j, k, iP, kP, iD, iL, n,m, &
	 & totalIndex=0

     !*** Initilize variables
 
    nwv = cfProd%nWave

    nlev = 0
    DO j = 1, l2gp(1)%nLevels
       IF( l2gp(1)%pressures(j) >= cfProd%l3presLvl(1) .AND. &
            & l2gp(1)%pressures(j) <= cfProd%l3presLvl(2) ) THEN
          nlev = nlev + 1
          IF (nlev == 1) pStartIndex = j 
          pEndIndex = j 
       ENDIF
    ENDDO
     
    IF (nlev == 0) THEN
       nlev_temp = 1
       nlev = 1
       pStartIndex = 1
       pEndIndex = 1
    ELSE
       nlev_temp = -1
    END IF
     
    totalIndex = pEndIndex - pStartIndex + 1 

    !*** Initilize POINTERS
     
    IF (cfProd%mode == 'all' .or. cfProd%mode == 'ado') THEN
       numSwaths = 3
    ELSE
       numSwaths = 1
    ENDIF
     
    ALLOCATE( l3sp(numSwaths), STAT=error )
    IF ( error /= 0 ) THEN
       msr = MLSMSG_Allocate // ' l3sp array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
     
    ALLOCATE( l3spPrec(numSwaths), STAT=error )
    IF ( error /= 0 ) THEN
       msr = MLSMSG_Allocate // ' l3spPrec array.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
     
    IF (cfProd%mode == 'all' .or. cfProd%mode == 'ado') THEN
       l3sp(1)%name = cfProd%l3prodNameD
       l3sp(2)%name = TRIM(cfProd%l3prodNameD) // 'Ascending'
       l3sp(3)%name = TRIM(cfProd%l3prodNameD) // 'Descending'
    ELSE
       l3sp(1)%name = TRIM(cfProd%l3prodNameD)
    ENDIF
     
    l3sp%startTime = l2gp(1)%time(1)
     
    if (l2Days .GT. 0) then 
       l3sp%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)
    else
       l3sp%endTime = l3sp%startTime
    endif
    
    DO j = 1, numSwaths
        
       CALL AllocateL3SP( nlev, cfProd%nLats, nwv, l3sp(j) )
       CALL AllocateL3SP( nlev, cfProd%nLats, nwv, l3spPrec(j) )
       l3sp(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
       l3sp(j)%latitude = cfProd%latGridMap(:l3sp(j)%nLats)
       l3sp(j)%waveNumber = 0.0
       l3sp(j)%frequency = 0.0
       l3sp(j)%l3spRelValue = -999.99 
       l3sp(j)%l3spRelPrecision = -999.99 
       l3sp(j)%l3spImgValue = -999.99 
       l3sp(j)%l3spImgPrecision = -999.99 
        
       l3spPrec(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
       l3spPrec(j)%latitude = cfProd%latGridMap(:l3sp(j)%nLats)
       l3spPrec(j)%waveNumber = 0.0
       l3spPrec(j)%frequency = 0.0
       l3spPrec(j)%l3spRelValue = -999.99 
       l3spPrec(j)%l3spRelPrecision = -999.99 
       l3spPrec(j)%l3spImgValue = -999.99 
       l3spPrec(j)%l3spImgPrecision = -999.99 
        
    ENDDO
     
    numDays = cfProd%nDays

    !!      Initialize Daily Map & Diagnostic

    if (numDays .gt. 0) then
       IF (cfProd%mode == 'com') THEN        
          ALLOCATE( l3dm(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l3dm array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          l3dm%name = TRIM(cfProd%l3prodNameD)
          DO j = 1, numDays
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, l3dm(j) )
             l3dm(j)%time      = cfProd%timeD(j)
             l3dm(j)%pressure  = l2gp(1)%pressures(pStartIndex:pEndIndex)
       	     l3dm(j)%latitude  = cfProd%latGridMap(:l3dm(j)%nLats)
             l3dm(j)%longitude = cfProd%longGrid(:l3dm(j)%nLons)
             l3dm(j)%l3dmValue = -999.99 
             l3dm(j)%l3dmPrecision = -999.99 
             l3dm(j)%gRss  = -999.99 
             l3dm(j)%perMisPoints  = 0 
             do n=1,cfDef%N
               do m=1,nlev
                  l3dm(j)%latRss(n,m) = 0.0
                  l3dm(j)%maxDiff(n,m) = 0.0
                  l3dm(j)%maxDiffTime(n,m)  = 0.0
               end do
             end do
          ENDDO

       ELSE IF (cfProd%mode == 'asc') THEN        
          ALLOCATE( dmA(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmA array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          dmA%name  = TRIM(cfProd%l3prodNameD) // 'Ascending'
          DO j = 1, numDays
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmA(j) )
             dmA(j)%time       = cfProd%timeD(j)
             dmA(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmA(j)%latitude   = cfProd%latGridMap(:dmA(j)%nLats)
             dmA(j)%longitude  = cfProd%longGrid(:dmA(j)%nLons)
             dmA(j)%l3dmValue  = -999.99 
             dmA(j)%l3dmPrecision  = -999.99 
             dmA(j)%gRss  = -999.99 
             dmA(j)%perMisPoints  = 0 
             do n=1,cfDef%N
               do m=1,nlev
                  dmA(j)%latRss(n,m)  = 0.0
                  dmA(j)%maxDiff(n,m)  = 0.0
                  dmA(j)%maxDiffTime(n,m)  = 0.0
               end do
             end do
          ENDDO

       ELSE IF (cfProd%mode == 'dsc') THEN        
          ALLOCATE( dmD(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmD array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          dmD%name  = TRIM(cfProd%l3prodNameD) // 'Descending'
          DO j = 1, numDays
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmD(j) )
             dmD(j)%time       = cfProd%timeD(j)
             dmD(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmD(j)%latitude   = cfProd%latGridMap(:dmD(j)%nLats)
             dmD(j)%longitude  = cfProd%longGrid(:dmD(j)%nLons)
             dmD(j)%l3dmValue  = -999.99 
             dmD(j)%l3dmPrecision  = -999.99 
             dmD(j)%gRss  = -999.99 
             dmD(j)%perMisPoints  = 0 
             do n=1,cfDef%N
               do m=1,nlev
                 dmD(j)%latRss(n,m)  = 0.0
                 dmD(j)%maxDiff(n,m)  = 0.0
                 dmD(j)%maxDiffTime(n,m) = 0.0
               end do
             end do
          ENDDO

       ELSE IF (cfProd%mode == 'ado') THEN        
          ALLOCATE( dmA(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmA array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          ALLOCATE( dmD(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmD array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          dmA%name  = TRIM(cfProd%l3prodNameD) // 'Ascending'
          dmD%name  = TRIM(cfProd%l3prodNameD) // 'Descending'

          DO j = 1, numDays
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmD(j) )
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmA(j) )
             dmA(j)%time       = cfProd%timeD(j)
             dmD(j)%time       = cfProd%timeD(j)
             dmA(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmD(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmA(j)%latitude   = cfProd%latGridMap(:dmA(j)%nLats)
             dmD(j)%latitude   = cfProd%latGridMap(:dmD(j)%nLats)
             dmA(j)%longitude  = cfProd%longGrid(:dmA(j)%nLons)
             dmD(j)%longitude  = cfProd%longGrid(:dmD(j)%nLons)
             dmA(j)%l3dmValue  = -999.99 
             dmD(j)%l3dmValue  = -999.99 
             dmA(j)%l3dmPrecision  = -999.99 
             dmD(j)%l3dmPrecision  = -999.99 
             dmA(j)%gRss  = -999.99 
             dmA(j)%perMisPoints  = 0 
             dmD(j)%gRss  = -999.99 
             dmD(j)%perMisPoints  = 0 
             do n=1,cfDef%N
               do m=1,nlev
                  dmA(j)%latRss(n,m)  = 0.0
                  dmA(j)%maxDiff(n,m)  = 0.0
                  dmA(j)%maxDiffTime(n,m)  = 0.0
                  dmD(j)%latRss(n,m)  = 0.0
                  dmD(j)%maxDiff(n,m)  = 0.0
                  dmD(j)%maxDiffTime(n,m)  = 0.0
               end do
             end do
          ENDDO

       ELSE IF (cfProd%mode == 'all') THEN        
          ALLOCATE( l3dm(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l3dm array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          ALLOCATE( dmA(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmA array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          ALLOCATE( dmD(numDays), STAT=error )
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dmD array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          l3dm%name = TRIM(cfProd%l3prodNameD)
          dmA%name  = TRIM(cfProd%l3prodNameD) // 'Ascending'
          dmD%name  = TRIM(cfProd%l3prodNameD) // 'Descending'
          DO j = 1, numDays
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, l3dm(j) )
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmD(j) )
             CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmA(j) )
             l3dm(j)%time      = cfProd%timeD(j)
             dmA(j)%time       = cfProd%timeD(j)
             dmD(j)%time       = cfProd%timeD(j)
             l3dm(j)%pressure  = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmA(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
             dmD(j)%pressure   = l2gp(1)%pressures(pStartIndex:pEndIndex)
       	     l3dm(j)%latitude  = cfProd%latGridMap(:l3dm(j)%nLats)
             dmA(j)%latitude   = cfProd%latGridMap(:dmA(j)%nLats)
             dmD(j)%latitude   = cfProd%latGridMap(:dmD(j)%nLats)
             l3dm(j)%longitude = cfProd%longGrid(:l3dm(j)%nLons)
             dmA(j)%longitude  = cfProd%longGrid(:dmA(j)%nLons)
             dmD(j)%longitude  = cfProd%longGrid(:dmD(j)%nLons)
             l3dm(j)%l3dmValue = -999.99 
             dmA(j)%l3dmValue  = -999.99 
             dmD(j)%l3dmValue  = -999.99 
             l3dm(j)%l3dmPrecision = -999.99 
             dmA(j)%l3dmPrecision  = -999.99 
             dmD(j)%l3dmPrecision  = -999.99 
             l3dm(j)%gRss  = -999.99 
             l3dm(j)%perMisPoints  = 0 
             dmA(j)%gRss  = -999.99 
             dmA(j)%perMisPoints  = 0 
             dmD(j)%gRss  = -999.99 
             dmD(j)%perMisPoints  = 0 
             do n=1,cfDef%N
               do m=1,nlev
                  l3dm(j)%latRss(n,m)  = 0.0
                  l3dm(j)%maxDiff(n,m)  = 0.0
                  l3dm(j)%maxDiffTime(n,m)  = 0.0
                  dmA(j)%latRss(n,m)  = 0.0
                  dmA(j)%maxDiff(n,m)  = 0.0
                  dmA(j)%maxDiffTime(n,m)  = 0.0
                  dmD(j)%latRss(n,m)  = 0.0
                  dmD(j)%maxDiff(n,m)  = 0.0
                  dmD(j)%maxDiffTime(n,m)  = 0.0
               end do
             end do
          ENDDO
       ENDIF
    ENDIF
    
    if ( numDays .gt. 0) then
       ALLOCATE( nc(numDays), STAT=error )
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' nc array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ALLOCATE( nca(numDays), STAT=error )
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' nca array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ALLOCATE( ncd(numDays), STAT=error )
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' ncd array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ALLOCATE( startTime(numDays), STAT=error )
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' startTime array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ALLOCATE( endTime(numDays), STAT=error )
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' endTime array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
    endif
 
    !!      Initialize Daily Map Residues
  
    IF (cfProd%mode == 'com') THEN
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
	   & pcf%l3StartDay, & 
           & pcf%l3EndDay, rDays, mis_l2Days_temp, mis_Days_temp, l3r_temp)
       CALL SetupL2GPProd(rDays, totalIndex, l3r_temp, l3r)
       l3r%name = TRIM(cfProd%l3prodNameD) // 'Residuals'
       DO j = 1, rDays
           l3r(j)%pressures = l3r_temp(j)%pressures(pStartIndex:pEndIndex)
           l3r(j)%l2gpValue = -999.99 
           l3r(j)%l2gpPrecision = l3r_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:)
           l3r_temp(j)%l2gpValue    = -999.99 
           DO k = 1, l3r(j)%nTimes
             if(l3r(j)%time(k) > 0) then
               startTime(j) = l3r(j)%time(k)
               exit
             endif
           ENDDO
           DO k = l3r(j)%nTimes, 1, -1
             if(l3r(j)%time(k) > 0) then
               endTime(j) = l3r(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO

    ELSE IF (cfProd%mode == 'asc') THEN
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
	   & pcf%l3StartDay, & 
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residA_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residA_temp, residA)
       residA%name = TRIM(cfProd%l3prodNameD) // 'AscendingResiduals'
       DO j = 1, rDays
           residA(j)%pressures = residA_temp(j)%pressures(pStartIndex:pEndIndex)
           residA(j)%l2gpValue = -999.99
           residA(j)%l2gpPrecision = residA_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           residA_temp(j)%l2gpValue = -999.99
           DO k = 1, residA(j)%nTimes
             if(residA(j)%time(k) > 0) then
               startTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residA(j)%nTimes, 1, -1
             if(residA(j)%time(k) > 0) then
               endTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO

    ELSE IF (cfProd%mode == 'dsc') THEN
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, & 
           & pcf%l3StartDay, & 
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residD_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residD_temp, residD)
       residD%name = TRIM(cfProd%l3prodNameD) // 'DescendingResiduals'
       DO j = 1, rDays
           residD(j)%pressures = residD_temp(j)%pressures(pStartIndex:pEndIndex)
           residD(j)%l2gpValue      = -999.99
           residD_temp(j)%l2gpValue = -999.99
           residD(j)%l2gpPrecision = residD_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           DO k = 1, residD(j)%nTimes
             if(residD(j)%time(k) > 0) then
               startTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residD(j)%nTimes, 1, -1
             if(residD(j)%time(k) > 0) then
               endTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO
 
    ELSE IF (cfProd%mode == 'ado') THEN
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
	   & pcf%l3StartDay, & 
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residA_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residA_temp, residA)
       residA%name = TRIM(cfProd%l3prodNameD) // 'AscendingResiduals'
       
       DO j = 1, rDays
           residA(j)%pressures = residA_temp(j)%pressures(pStartIndex:pEndIndex)
           residA(j)%l2gpValue = -999.99
           residA_temp(j)%l2gpValue = -999.99
           residA(j)%l2gpPrecision = residA_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           DO k = 1, residA(j)%nTimes
             if(residA(j)%time(k) > 0) then
               startTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residA(j)%nTimes, 1, -1
             if(residA(j)%time(k) > 0) then
               endTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO
       rDays = 0
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
           & pcf%l3StartDay, & 
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residD_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residD_temp, residD)
       residD%name = TRIM(cfProd%l3prodNameD) // 'DescendingResiduals'
       DO j = 1, rDays
           residD(j)%pressures = residD_temp(j)%pressures(pStartIndex:pEndIndex)
           residD(j)%l2gpValue      = -999.99
           residD_temp(j)%l2gpValue = -999.99
           residD(j)%l2gpPrecision = residD_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           DO k = 1, residD(j)%nTimes
             if(residD(j)%time(k) > 0) then
               startTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residD(j)%nTimes, 1, -1
             if(residD(j)%time(k) > 0) then
               endTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO

    ELSE IF (cfProd%mode == 'all') THEN
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
	   & pcf%l3StartDay, & 
           & pcf%l3EndDay, rDays, mis_l2Days_temp, mis_Days_temp, l3r_temp)
       CALL SetupL2GPProd(rDays, totalIndex, l3r_temp, l3r)
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
           & pcf%l3StartDay, &
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residA_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residA_temp, residA)
       CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, &
           & pcf%l3StartDay, &
           & pcf%l3EndDay,rDays, mis_l2Days_temp, mis_Days_temp, residD_temp)
       CALL SetupL2GPProd(rDays, totalIndex, residD_temp, residD)
       l3r%name    = TRIM(cfProd%l3prodNameD) // 'Residuals'
       residA%name = TRIM(cfProd%l3prodNameD) // 'AscendingResiduals'
       residD%name = TRIM(cfProd%l3prodNameD) // 'DescendingResiduals'
       DO j = 1, rDays
           l3r(j)%pressures = l3r_temp(j)%pressures(pStartIndex:pEndIndex)
           l3r(j)%l2gpValue         = -999.99 
           l3r_temp(j)%l2gpValue    = -999.99 
           l3r(j)%l2gpPrecision = l3r_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:)
           residA(j)%pressures = residA_temp(j)%pressures(pStartIndex:pEndIndex)
           residA(j)%l2gpValue = -999.99
           residA_temp(j)%l2gpValue = -999.99
           residA(j)%l2gpPrecision = residA_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           residD(j)%pressures = residD_temp(j)%pressures(pStartIndex:pEndIndex)
           residD(j)%l2gpValue      = -999.99
           residD_temp(j)%l2gpValue = -999.99
           residD(j)%l2gpPrecision = residD_temp(j)%l2gpPrecision(:,pStartIndex:pEndIndex,:) 
           DO k = 1, l3r(j)%nTimes
             if(l3r(j)%time(k) > 0) then
               startTime(j) = l3r(j)%time(k)
               exit
             endif
           ENDDO
           DO k = l3r(j)%nTimes, 1, -1
             if(l3r(j)%time(k) > 0) then
               endTime(j) = l3r(j)%time(k)
               exit
             endif
           ENDDO
           DO k = 1, residA(j)%nTimes
             if(residA(j)%time(k) > 0) then
               startTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residA(j)%nTimes, 1, -1
             if(residA(j)%time(k) > 0) then
               endTime(j) = residA(j)%time(k)
               exit
             endif
           ENDDO
           DO k = 1, residD(j)%nTimes
             if(residD(j)%time(k) > 0) then
               startTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
           DO k = residD(j)%nTimes, 1, -1
             if(residD(j)%time(k) > 0) then
               endTime(j) = residD(j)%time(k)
               exit
             endif
           ENDDO
       ENDDO

    ENDIF
        
    !!      Initialize Flags
        
    flags%writel3sp    = .FALSE.
    flags%writel3dmCom = .FALSE.
    flags%writel3dmAsc = .FALSE.
    flags%writel3dmDes = .FALSE.
    flags%writel3rCom  = .FALSE.
    flags%writel3rAsc  = .FALSE.
    flags%writel3rDes  = .FALSE.
       
    !!      Check if pressure levels are found 
        
    IF (nlev_temp == 1) THEN 
       RETURN
    END IF
    
    !*** Calculate average orbital period (day)
    
    tau0 = 0.0
    DO I = 1, size(avgPeriod)
       tau0 = tau0 + avgPeriod(i)
    ENDDO
        
    IF (size(avgPeriod) .ne. 0) tau0 = tau0/86400.0/float(size(avgPeriod))
        
    !*** Sort & Prepare the Data 
    
    Call SortData(cfProd, l2Days, l2gp, 	&
         & pStartIndex, pEndIndex,		&
         & tau0, 				&
         & anlats, dnlats, 			&
         & alats, dlats, 			&
         & alons, dlons, 			&
         & atimes, dtimes, 			&
         & afields, dfields, 		        &
         & aprec, dprec,                        &
         & delTad, perMisPoints )
        
        !*** Main Loop 
        
        !*** Calculate Field Values ********
        
    iP = 0
    DO kP = pStartIndex, pEndIndex 
	
       iP = iP + 1
       DO I = 1, rDays 
          nc(I) = 0
          nca(I) = 0
          ncd(I) = 0
       ENDDO
           
       DO J = 1, cfProd%nLats
          IF( anlats(J, iP) > 0 ) THEN
             
             IF( atimes(J, 1, iP) < dtimes(J, 1, iP) ) THEN
                
                lonD0_in = FindRealLon(real(dlons(J, 1, iP)))
                lonA0_in = FindRealLon(real(alons(J, 1, iP)))
                    
                CALL Init( & 
                     & nt_a_i   = anlats(J, iP), 		&
                     & nt_d_i   = dnlats(J, iP), 		&
                     & tau0_i   = tau0, 			&
                     & delTad_i = delTad(J, iP), 		&
                     & c0_i     = 2.0*PI, 			&
                     & lonD0_i  = lonD0_in, 			&
                     & tD0_i    = dtimes(J, 1, iP), 		&
                     & lonA0_i  = lonA0_in, 			&
                     & tA0_i    = atimes(J, 1, iP), 		&
                     & lat_i    = alats(J, 1, iP) )
                CALL DataGenerate(afields(J, :, iP), dfields(J, :, iP) )
                CALL DataGeneratePrec(aprec(J, :, iP), dprec(J, :, iP) )
             ELSE
                    
                lonD0_in = FindRealLon(real(dlons(J, 2, iP)))
                lonA0_in = FindRealLon(real(alons(J, 1, iP)))
                    
                CALL Init( & 
                     & nt_a_i   = anlats(J, iP), 			&
                     & nt_d_i   = dnlats(J, iP)-1, 			&
                     & tau0_i   = tau0, 				&
                     & delTad_i = delTad(J, iP), 			&
                     & c0_i     = 2.0*PI, 				&
                     & lonD0_i  = lonD0_in, 			        &
                     & tD0_i    = dtimes(J, 2, iP), 		        & 
                     & lonA0_i  = lonA0_in, 			        &
                     & tA0_i    = atimes(J, 1, iP), 		        &
                     & lat_i    = alats(J, 1, iP) )
                CALL DataGenerate(afields(J, :, iP), dfields(J, 2:, iP) )
                CALL DataGeneratePrec(aprec(J, :, iP), dprec(J, 2:, iP) )
             END IF
                
             IF (cfProd%mode == 'com') THEN
                
                if (l3dm(1)%nLons .gt. 0) then 
                   
                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                       
                   CALL CordTransform(cfProd%mode)
                   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, l3dm(iD)%nLons
                         avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (l3dm(iD)%nLons .ne. 0) & 
                           & avg = avg/real(l3dm(iD)%nLons)
                      l3dm(iD)%latRss(iP, J) = 0.0
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%latRss(iP, J) = & 
                              & l3dm(iD)%latRss(iP, J) +  &
                              & (l3dm(iD)%l3dmValue(iP, J, I)-avg)*  &
                              & (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                          
                      IF (l3dm(iD)%nLons .ne. 0) THEN
                         l3dm(iD)%latRss(iP, J) = sqrt(l3dm(iD)%latRss(iP, J)/ &
                              & real(l3dm(iD)%nLons-1))
                      ENDIF
                          
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
! write(*,*) iD, startTime(iD), endTime(iD) 
                      DO iL = 1, anlats(J, iP)
                         IF( & 
                              &((atimes(J, iL, iP)*86400.+l2gp(1)%time(1)) .ge. startTime(iD)) .AND. &
                              &((atimes(J, iL, iP)*86400.+l2gp(1)%time(1)) .le. endTime(iD))) THEN
                            nc(iD) = nc(iD) + 1
                            CALL Diagnostics(cfProd%mode, atimes(J, iL, iP), alons(J, iL, iP), l3ret)
                            l3r_temp(iD)%time(nc(iD))  = atimes(J, iL,iP)*86400.0+l2gp(1)%time(1)
                            l3r_temp(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
                            l3r_temp(iD)%longitude(nc(iD)) = FindRealLon(real(alons(J, iL, iP)))*180.0/PI
                            l3r_temp(iD)%l2gpValue(1, kP, nc(iD)) = afields(J, iL, iP)-l3ret 
if((afields(J, iL, iP)-l3ret)/afields(J, iL, iP)*100.0 > 10) then
! write(*,*) iD, afields(J, iL, iP), l3ret, alons(J, iL, iP), atimes(J, iL, iP), l3r_temp(iD)%time(nc(iD)), 'a'
end if
if(afields(J, iL, iP) < 0.0) then
! write(*,*) iD, afields(J, iL, iP), l3ret, alons(J, iL, iP), atimes(J, iL, iP), l3r_temp(iD)%time(nc(iD)), 'a'
end if
                            l3r_temp(iD)%l2gpPrecision(1, kP, nc(iD)) = 0.0
                         END IF
                             
                      ENDDO

                      DO iL = 1, dnlats(J, iP)
                         IF( &
                              & dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
                              & dtimes(J,iL,iP)*86400.0+l2gp(1)%time(1)<= endTime(iD)) THEN
                            nc(iD) = nc(iD) + 1
                            CALL Diagnostics(cfProd%mode, dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
                            l3r_temp(iD)%time(nc(iD)) = dtimes(J,iL, iP)*86400.0+l2gp(1)%time(1)
                            l3r_temp(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
                            l3r_temp(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
                            l3r_temp(iD)%l2gpValue(1, kP, nc(iD)) = dfields(J, iL, iP)-l3ret 
if((dfields(J, iL, iP)-l3ret)/dfields(J, iL, iP)*100.0 > 10) then
 !write(*,*) iD, dfields(J, iL, iP), l3ret, 'd'
end if
                            l3r_temp(iD)%l2gpPrecision(1, kP, nc(iD)) = 0.0
                         END IF
                      ENDDO
                   ENDDO
 
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
 
		   l3sp(1)%l3spRelPrecision = l3spPrec(1)%l3spRelValue
		   l3sp(1)%l3spImgPrecision = l3spPrec(1)%l3spImgValue
                       
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result, STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                       
                   flags%writel3dmCom = .TRUE.
                   flags%writel3rCom  = .TRUE.
                   flags%writel3sp = .TRUE.
                  
                endif
                    
             ELSE IF (cfProd%mode == 'asc') THEN
                    
                if (dmA(1)%nLons .gt. 0) then 
                       
                   if ( associated(l3Result) ) then 
                      ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                       
                   CALL CordTransform(cfProd%mode)
                   CALL FFSMA(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, & 
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                         avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) & 
                           & avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) +  &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN 
                         dmA(iD)%latRss(iP, J) = sqrt(dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons-1))
                      ENDIF
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
                         IF & 
                              & (atimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & atimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            nc(iD) = nc(iD) + 1
                            CALL Diagnostics(cfProd%mode, & 
                                 & atimes(J, iL, iP), & 
                                 & alons(J, iL, iP), l3ret) 
                            residA_temp(iD)%time(nc(iD))     = & 
                                 & atimes(J,iL, iP)*86400.0+l2gp(1)%time(1)
                            residA_temp(iD)%latitude(nc(iD)) = & 
                                 & cfProd%latGridMap(J) 
                            residA_temp(iD)%longitude(nc(iD)) = & 
                                 & FindRealLon(real(alons(J, iL, iP)))* & 
                                 & 180.0/PI
                            residA_temp(iD)%l2gpValue(1, kP, nc(iD)) = & 
                                 & afields(J, iL, iP)-l3ret 
                            residA_temp(iD)%l2gpPrecision(1, kP, nc(iD)) =&
                                 & 0.0
                         END IF
                      ENDDO
                   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(1)%l3spRelPrecision = l3spPrec(1)%l3spRelValue
		   l3sp(1)%l3spImgPrecision = l3spPrec(1)%l3spImgValue
                  
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result, STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                       
                   flags%writel3dmAsc = .TRUE.
                   flags%writel3rAsc  = .TRUE.
                   flags%writel3sp = .TRUE.
                  
                endif
                    
             ELSE IF (cfProd%mode == 'des') THEN
                    
                if (dmD(1)%nLons .gt. 0) then 
                   
                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                       
                   CALL CordTransform(cfProd%mode)
                   CALL FFSMD(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, & 
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                         avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) & 
                           & avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) +  &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)*  &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN
                         dmD(iD)%latRss(iP, J) = sqrt(dmD(iD)%latRss(iP, J)/real(dmD(iD)%nLons-1))
                      ENDIF
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
                         IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & dtimes(J,iL,iP)*86400.0+l2gp(1)%time(1)<=& 
                              & endTime(iD)) THEN
                            nc(iD) = nc(iD) + 1
                            CALL Diagnostics(cfProd%mode, &
                                 & dtimes(J, iL, iP), & 
                                 & dlons(J, iL, iP), l3ret) 
                            residD_temp(iD)%time(nc(iD))     = & 
                                 & dtimes(J, iL,iP)*86400.0+l2gp(1)%time(1)
                            residD_temp(iD)%latitude(nc(iD)) = & 
                                 & cfProd%latGridMap(J) 
                            residD_temp(iD)%longitude(nc(iD)) = & 
                                 & FindRealLon(real(dlons(J, iL, iP)))* & 
                                 & 180.0/PI
                            residD_temp(iD)%l2gpValue(1, kP, nc(iD)) = & 
                                 & dfields(J, iL, iP)-l3ret 
                            residD_temp(iD)%l2gpPrecision(1, kP, nc(iD))= &
                                 & 0.0 
                         END IF
                      ENDDO
                   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(1)%l3spRelPrecision = l3spPrec(1)%l3spRelValue
		   l3sp(1)%l3spImgPrecision = l3spPrec(1)%l3spImgValue
                  
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result, STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                   
                   flags%writel3dmDes = .TRUE.
                   flags%writel3rDes  = .TRUE.
                   flags%writel3sp    = .TRUE.
                       
                endif
                    
             ELSE IF (cfProd%mode == 'ado') THEN
                
                if (dmA(1)%nLons .gt. 0) then 
                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                       
                   CALL CordTransform('asc')
                   CALL FFSMA(l3sp(2), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct('asc', & 
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0,  &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmValue(iP, J, I) = l3Result(I)
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                         avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) & 
                           & avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) +  &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)* & 
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN
                         dmA(iD)%latRss(iP, J) = sqrt(dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons-1))
                      ENDIF
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
                         
                         IF & 
                              & (atimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & atimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            nca(iD) = nca(iD) + 1
                            CALL Diagnostics('asc', atimes(J, iL, iP), & 
                                 & alons(J, iL, iP), l3ret)
                            
                            residA_temp(iD)%time(nca(iD)) = & 
                                 & atimes(J, iL,iP)*86400.0+l2gp(1)%time(1)
                            residA_temp(iD)%longitude(nca(iD)) = & 
                                 &FindRealLon(real(alons(J,iL,iP)))*180./PI
                            residA_temp(iD)%latitude(nca(iD)) = & 
                                 & cfProd%latGridMap(J)  
                            residA_temp(iD)%l2gpValue(1, kP, nca(iD)) = & 
                                 & afields(J, iL, iP)-l3ret
                            residA_temp(iD)%l2gpPrecision(1, kP, nca(iD)) &
                                 & = 0.0
                                
                         END IF
                      ENDDO
                   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(2), iP, J)
                   DO iD = 1, cfProd%nDays
                      DO I = 1, dmA(1)%nLons
                         l3Result(I) = 0.0 
                      ENDDO
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(2)%l3spRelPrecision = l3spPrec(2)%l3spRelValue
		   l3sp(2)%l3spImgPrecision = l3spPrec(2)%l3spImgValue
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      Deallocate(l3Result, STAT=error) 
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                       
                   flags%writel3dmAsc = .TRUE. 
                   flags%writel3rAsc  = .TRUE.
                       
                endif
                    
                if (dmD(1)%nLons .gt. 0) then 
                       
                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // '  l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                       
                   CALL CordTransform('des')
                   CALL FFSMD(l3sp(3), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct('des', & 
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmValue(iP, J, I) = l3Result(I)
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                         avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) & 
                           & avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) + &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)* & 
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN 
                         dmD(iD)%latRss(iP, J) = sqrt(dmD(iD)%latRss(iP, J)/ &
                              & real(dmD(iD)%nLons-1))
                      ENDIF
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
                         IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & dtimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            ncd(iD) = ncd(iD) + 1
                            CALL Diagnostics('des', dtimes(J, iL, iP), & 
                                 & dlons(J, iL, iP), l3ret)
                                
                            residD_temp(iD)%longitude(ncd(iD))       = & 
                                 & FindRealLon(real(dlons(J, iL, iP)))*& 
                                 & 180./PI

                            residD_temp(iD)%time(ncd(iD))            = & 
                                 & dtimes(J,iL, iP)*86400.0+l2gp(1)%time(1)
                            
                            residD_temp(iD)%latitude(ncd(iD))        = & 
                                 & cfProd%latGridMap(J)  
                                
                            residD_temp(iD)%l2gpValue(1, kP, ncd(iD))= & 
                                 & dfields(J, iL, iP)-l3ret

                            residD_temp(iD)%l2gpPrecision(1, kP, ncd(iD)) &
                                 & = 0.0
                                
                         END IF
                      ENDDO
                   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(3), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(3)%l3spRelPrecision = l3spPrec(3)%l3spRelValue
		   l3sp(3)%l3spImgPrecision = l3spPrec(3)%l3spImgValue
                       
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result, STAT=error) 
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                       
                   flags%writel3dmDes = .TRUE. 
                   flags%writel3rDes  = .TRUE.
                   flags%writel3sp    = .TRUE.
                       
                endif
                    
             ELSE IF (cfProd%mode == 'all') THEN
                    
                if (l3dm(1)%nLons .gt. 0) then 
                       
                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                       
                   CALL CordTransform('com')
                   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct('com', & 
                           & real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
                      ENDDO
                          !***Root-Sum-Square for each latitude, 
                          ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, l3dm(iD)%nLons
                         avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (l3dm(iD)%nLons .ne. 0) & 
                           & avg = avg/real(l3dm(iD)%nLons)
                      l3dm(iD)%latRss(iP, J) = 0.0
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%latRss(iP, J) = l3dm(iD)%latRss(iP, J) +&
                              & (l3dm(iD)%l3dmValue(iP, J, I)-avg)*  &
                              & (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (l3dm(iD)%nLons .ne. 0) THEN
                         l3dm(iD)%latRss(iP, J) = sqrt(l3dm(iD)%latRss(iP, J)/ &
                              & real(l3dm(iD)%nLons-1))
                      ENDIF
                   ENDDO
                       
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
                         IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & atimes(J, iL, iP)*86400.0+ & 
                              & l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            nc(iD) = nc(iD) + 1
                            CALL Diagnostics('com', atimes(J, iL, iP), & 
                                 & alons(J, iL, iP), l3ret) 
                            l3r_temp(iD)%time(nc(iD))  = & 
                                 & atimes(J,iL, iP)*86400.0+l2gp(1)%time(1)
                            l3r_temp(iD)%latitude(nc(iD)) = & 
                                 & cfProd%latGridMap(J) 
                            l3r_temp(iD)%longitude(nc(iD)) = & 
                                 & FindRealLon(real(alons(J, iL, iP)))* & 
                                 & 180.0/PI
                            l3r_temp(iD)%l2gpValue(1, kP, nc(iD)) = & 
                                 & afields(J, iL, iP)-l3ret 
                            l3r_temp(iD)%l2gpPrecision(1, kP, nc(iD)) = 0.0
                         END IF
                      ENDDO
                      DO iL = 1, dnlats(J, iP)
                         IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & dtimes(J, iL, iP)*86400.+l2gp(1)%time(1)<=&
                              & endTime(iD)) THEN
                            
                            nc(iD) = nc(iD) + 1
                                
                            CALL Diagnostics('com', dtimes(J, iL, iP), & 
                                 & dlons(J, iL, iP), l3ret) 
                                
                            l3r_temp(iD)%time(nc(iD))                 = & 
                                 & dtimes(J, iL, iP)*86400.+l2gp(1)%time(1)
                            
                            l3r_temp(iD)%latitude(nc(iD))             = & 
                                 & cfProd%latGridMap(J) 
                                
                            l3r_temp(iD)%longitude(nc(iD))            = & 
                                 & FindRealLon(real(dlons(J, iL, iP)))* & 
                                 & 180.0/PI
                                
                            l3r_temp(iD)%l2gpValue(1, kP, nc(iD))     = & 
                                 & dfields(J, iL, iP)-l3ret 

                            l3r_temp(iD)%l2gpPrecision(1, kP, nc(iD)) = 0.0

                         END IF
		      ENDDO
		   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
                         l3dm(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(1)%l3spRelPrecision = l3spPrec(1)%l3spRelValue
		   l3sp(1)%l3spImgPrecision = l3spPrec(1)%l3spImgValue
                   
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result,STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                   
		   flags%writel3dmCom = .TRUE.
		   flags%writel3rCom  = .TRUE.
                   
                endif
                
                if (dmA(1)%nLons .gt. 0) then 
                   
                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3sp(2), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('asc', & 
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
            !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                         avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) + &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN
                         dmA(iD)%latRss(iP, J) = sqrt(dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons-1))
                      ENDIF
                   ENDDO
                   
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
                         IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            
                            nca(iD) = nca(iD) + 1
                            
                            CALL Diagnostics('asc', atimes(J, iL, iP), & 
                                & alons(J, iL, iP), l3ret) 
                            
                            residA_temp(iD)%time(nca(iD))     = & 
                                 & atimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
                            
                            residA_temp(iD)%latitude(nca(iD)) = & 
                                 & cfProd%latGridMap(J) 
                            
                            residA_temp(iD)%longitude(nca(iD)) = & 
                                 & FindRealLon(real(alons(J, iL, iP)))*180.0/PI
                            
                            residA_temp(iD)%l2gpValue(1, kP, nca(iD)) = & 
                                 & afields(J, iL, iP)-l3ret 
                            
                            residA_temp(iD)%l2gpPrecision(1, kP, nca(iD)) = 0.0
                            
                         END IF
		      ENDDO
		   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(2), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(2)%l3spRelPrecision = l3spPrec(2)%l3spRelValue
		   l3sp(2)%l3spImgPrecision = l3spPrec(2)%l3spImgValue
                   
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result,STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                   
		   flags%writel3dmAsc = .TRUE.
		   flags%writel3rAsc  = .TRUE.
                   
                endif
                
                if (dmD(1)%nLons .gt. 0) then 
                   
                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_Allocate // ' l3Result array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3sp(3), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('des', & 
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, 
                      ! dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                         avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) +  &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN 
                         dmD(iD)%latRss(iP, J) = sqrt(dmD(iD)%latRss(iP, J)/real(dmD(iD)%nLons-1))
                      ENDIF
		   ENDDO
                   
                   !*** Calculate Residual     
                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
                         IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= & 
                              & startTime(iD) .AND. &
                              & dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= & 
                              & endTime(iD)) THEN
                            
                            ncd(iD) = ncd(iD) + 1

                            CALL Diagnostics('des', dtimes(J, iL, iP), & 
                                 & dlons(J, iL, iP), l3ret) 

                            residD_temp(iD)%time(ncd(iD))     = & 
                                 & dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
                            
                            residD_temp(iD)%latitude(ncd(iD)) = & 
                                 & cfProd%latGridMap(J) 
                            
                            residD_temp(iD)%longitude(ncd(iD)) = & 
                                 & FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
                            
                            residD_temp(iD)%l2gpValue(1, kP, ncd(iD)) = & 
                                 & dfields(J, iL, iP)-l3ret
                            
                            residD_temp(iD)%l2gpPrecision(1, kP, ncd(iD)) = 0.0

                         END IF
		      ENDDO
		   ENDDO
                  
                   !*** Calculate Precision     
                   Call CopyPrec2Data()
                   CALL FFSM(l3spPrec(3), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct(cfProd%mode, &
                           & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                           & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                      ENDDO
                   ENDDO
		   l3sp(3)%l3spRelPrecision = l3spPrec(3)%l3spRelValue
		   l3sp(3)%l3spImgPrecision = l3spPrec(3)%l3spImgValue
                   
                   !*** Cleanup    
                   if ( associated(l3Result) ) then 
                      DeAllocate(l3Result,STAT=error)
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_DeAllocate // '  l3Result array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                   endif
                   
		   flags%writel3dmDes = .TRUE.
		   flags%writel3rDes  = .TRUE.
		   flags%writel3sp    = .TRUE.
                   
                endif
                
             END IF
             
             CALL ClearMemory()
             
          END IF
       ENDDO


    !*** Sort into time ascending order according to l2gp format
       IF (cfProd%mode == 'com') THEN
          DO iD = 1, rDays
             
             if (nc(iD) .gt. 0) then 
                
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                ALLOCATE( sortTemp(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' sortTemp array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                DO i = 1, nc(iD)
                   pt(i) = 0
                   sortTemp(i) = 999.99 
                ENDDO
                CALL my_sortp(l3r_temp(iD)%time, 1, nc(iD), pt)
                !** do time
                CALL DSORT(l3r_temp(iD)%time, 1, nc(iD))
                !** do latitude
                DO i = 1, nc(iD)
                   sortTemp(i) = l3r_temp(iD)%latitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   l3r_temp(iD)%latitude(i) = sortTemp(pt(i))
                ENDDO
                !** do longitude
                DO i = 1, nc(iD)
                   sortTemp(i) = l3r_temp(iD)%longitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   l3r_temp(iD)%longitude(i) = sortTemp(pt(i))
                ENDDO
                !** do value
                DO i = 1, nc(iD)
                   sortTemp(i) = l3r_temp(iD)%l2gpValue(1, kP, i)
                ENDDO
                DO i = 1, nc(iD)
                   l3r_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                ENDDO
 
                !*** Interpolate to l2gp grid 
                !CALL Residual2L2Grid( iD, kP, l3r_temp, l3r ) 
                CALL Residual2L2Grid( iD, iP, kP, l3r_temp, l3r ) 

                if ( associated(sortTemp) ) then 
                   DeAllocate(sortTemp, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  sortTemp array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
                if ( associated(pt)) then 
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !*** Calculate Maximum Difference
             
             if ( nc(iD) .gt. 0 ) then 
                
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // '  pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                CALL my_sortp(-abs(l3r_temp(iD)%l2gpValue & 
                     (1, kP, 1:nc(iD))), 1, nc(iD), pt)
                DO i = 1, cfDef%N
                   l3dm(iD)%maxDiff(i, iP) = & 
                        l3r_temp(iD)%l2gpValue(1, kP, pt(i))
                ENDDO
                
                if ( associated(pt) ) then 
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !***Root-Sum-Square for each latitude, dimensioned (nLevels)
             avg = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, l3dm(iD)%nLons
                   avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                ENDDO
             ENDDO
             IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN
                avg = avg/real(l3dm(iD)%nLons*cfProd%nLats)
             ENDIF
             l3dm(iD)%gRss(iP) = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, l3dm(iD)%nLons
                   l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP) + &
                        & (l3dm(iD)%l3dmValue(iP, J, I)-avg)* &
                        & (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                ENDDO
             ENDDO
             IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN
                l3dm(iD)%gRss(iP) = sqrt(l3dm(iD)%gRss(iP)/real(l3dm(iD)%nLons*cfProd%nLats-1))
             ENDIF
             l3dm(iD)%perMisPoints(iP) = perMisPoints(iP)
          ENDDO
       ELSE IF (cfProd%mode == 'asc') THEN
          DO iD = 1, rDays
             
             if ( nc(iD) .gt. 0 ) then 
                
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                     
                ALLOCATE( sortTemp(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' sortTemp array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                CALL my_sortp(residA_temp(iD)%time, 1, nc(iD), pt)
                !** do time
                CALL DSORT(residA_temp(iD)%time, 1, nc(iD))
                !** do latitude
                DO i = 1, nc(iD)
                   sortTemp(i) = residA_temp(iD)%latitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   residA_temp(iD)%latitude(i) = sortTemp(pt(i))
                   !** do longitude
                ENDDO
                DO i = 1, nc(iD)
                   sortTemp(i) = residA_temp(iD)%longitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   residA_temp(iD)%longitude(i) = sortTemp(pt(i))
                   !** do value
                ENDDO
                DO i = 1, nc(iD)
                   sortTemp(i) = residA_temp(iD)%l2gpValue(1, kP, i)
                ENDDO
                DO i = 1, nc(iD)
                   residA_temp(iD)%l2gpValue(1, kP, i) = & 
                        sortTemp(pt(i))
                ENDDO
                
                !*** Interpolate to l2gp grid 
                !CALL Residual2L2Grid( iD, kP, residA_temp, residA ) 
                CALL Residual2L2Grid( iD, iP, kP, residA_temp, residA ) 
                
                if ( associated(sortTemp)) then
                   DeAllocate(sortTemp, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  sortTemp array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
                if ( associated(pt) ) then
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !*** Calculate Maximum Difference
             
             if ( nc(iD) .gt. 0 ) then 
                
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                         
                CALL my_sortp(-abs(residA_temp(iD)% & 
                     l2gpValue(1, kP, :)), 1, nc(iD), pt)
                DO i = 1,cfDef%N
                   dmA(iD)%maxDiff(i, kP) = & 
                        residA_temp(iD)%l2gpValue(1, kP, pt(i))
                ENDDO
                
                if ( associated(pt) ) then
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !***Root-Sum-Square for each latitude, dimensioned (nLevels)
             avg = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, dmA(iD)%nLons
                   avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                ENDDO
             ENDDO
             IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
             ENDIF
             dmA(iD)%gRss(iP) = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, dmA(iD)%nLons
                   dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                        & (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                        & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                ENDDO
             ENDDO
             IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN
                dmA(iD)%gRss(iP) = sqrt(dmA(iD)%gRss(iP)/real(dmA(iD)%nLons*cfProd%nLats-1))
             ENDIF
             dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
          ENDDO
       ELSE IF (cfProd%mode == 'des') THEN
          DO iD = 1, rDays
             
             if ( nc(iD) .gt. 0 ) then 
                
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                ALLOCATE( sortTemp(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' sortTemp array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                CALL my_sortp(residD_temp(iD)%time, 1, nc(iD), pt)
                !** do time
                CALL DSORT(residD_temp(iD)%time, 1, nc(iD))
                !** do latitude
                DO i = 1, nc(iD)
                   sortTemp(i) = residD_temp(iD)%latitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   residD_temp(iD)%latitude(i) = sortTemp(pt(i))
                ENDDO
                !** do longitude
                DO i = 1, nc(iD)
                   sortTemp(i) = residD_temp(iD)%longitude(i)
                ENDDO
                DO i = 1, nc(iD)
                   residD_temp(iD)%longitude(i) = sortTemp(pt(i))
                ENDDO
                !** do value
                DO i = 1, nc(iD)
                   sortTemp(i) = residD_temp(iD)%l2gpValue(1, kP, i)
                ENDDO
                DO i = 1, nc(iD)
                   residD_temp(iD)%l2gpValue(1, kP, i) = & 
                        sortTemp(pt(i))
                ENDDO
                
                !*** Interpolate to l2gp grid 
                !CALL Residual2L2Grid( iD, kP, residD_temp, residD ) 
                CALL Residual2L2Grid( iD, iP, kP, residD_temp, residD ) 
                
                if ( associated(sortTemp) ) then 
                   DeAllocate(sortTemp, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  sortTemp array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
                if ( associated(pt) ) then 
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !*** Calculate Maximum Difference
             
             if (nc(iD) .gt. 0) then 
                         
                ALLOCATE( pt(nc(iD)), STAT=error )
                IF ( error /= 0 ) THEN
                   msr = MLSMSG_Allocate // ' pt array.'
                   CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                ENDIF
                
                CALL my_sortp(-abs(residD_temp(iD)% & 
                     & l2gpValue(1, kP, :)), 1, nc(iD), pt)
                DO i = 1,cfDef%N
                   dmD(iD)%maxDiff(i, iP) = & 
                        & residD_temp(iD)%l2gpValue(1, kP, pt(i))
                ENDDO
                
                if ( associated(pt) ) then
                   DeAllocate(pt, STAT=error)
                   IF ( error /= 0 ) THEN
                      msr = MLSMSG_DeAllocate // '  pt array.'
                      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                endif
                
             endif
             
             !***Root-Sum-Square for each latitude, dimensioned (nLevels)
             avg = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, dmD(iD)%nLons
                   avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                ENDDO
             ENDDO
             IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN
                avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
             ENDIF
             dmD(iD)%gRss(iP) = 0.0
             DO J = 1, cfProd%nLats
                DO I = 1, dmD(iD)%nLons
                   dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                        & (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                        & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                ENDDO
             ENDDO
             IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                dmD(iD)%gRss(iP) = sqrt(dmD(iD)%gRss(iP)/real(dmD(iD)%nLons*cfProd%nLats-1))
             ENDIF
             dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
          ENDDO
       ELSE IF (cfProd%mode == 'ado') THEN
                DO iD = 1, rDays
                   if ( nca(iD) .gt. 0 ) then 
                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      ALLOCATE( sortTemp(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // '  sortTemp array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(residA_temp(iD)%time, 1, nca(iD), pt)
                      !** do time
                      CALL DSORT(residA_temp(iD)%time, 1, nca(iD))
                      !** do latitude
                      DO i = 1, nca(iD)
                         sortTemp(i) = residA_temp(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                         residA_temp(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nca(iD)
                         sortTemp(i) = residA_temp(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                         residA_temp(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nca(iD)
                         sortTemp(i) = residA_temp(iD)%l2gpValue(1, kP, i)
                      ENDDO
                      DO i = 1, nca(iD)
                         residA_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                      ENDDO
                      
                      !*** Interpolate to l2gp grid 
		      !CALL Residual2L2Grid( iD, kP, residA_temp, residA ) 
		      CALL Residual2L2Grid( iD, iP, kP, residA_temp, residA ) 
                      
                      if ( associated(sortTemp) ) then 
                         DeAllocate(sortTemp, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  sortTemp array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif
                   
                   !*** Calculate Maximum Difference
                   
                   if ( nca(iD) .gt. 0 ) then 
                      
                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(-abs(residA_temp(iD)% & 
                           l2gpValue(1, kP, :)), 1, nca(iD), pt)
                      DO i = 1,cfDef%N
                         dmA(iD)%maxDiff(i, iP) = & 
                              residA_temp(iD)%l2gpValue(1, kP, pt(i))
                      ENDDO
                      
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif
                   
                  !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                   avg = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                         avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                   ENDDO
                   
                   IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
                   ENDIF
                   
                   dmA(iD)%gRss(iP) = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                   ENDDO
                   
                   IF (dmA(iD)%nLons*cfProd%nLats .ne. 0 ) THEN 
                      dmA(iD)%gRss(iP) = sqrt(dmA(iD)%gRss(iP)/real(dmA(iD)%nLons*cfProd%nLats-1))
                   ENDIF
                   
                   dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
                
                DO iD = 1, rDays
                   
                   if (ncd(iD) .gt. 0) then 
                      
                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      ALLOCATE( sortTemp(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' sortTemp array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(residD_temp(iD)%time, 1, ncd(iD), pt)
                      !** do time
                      CALL DSORT(residD_temp(iD)%time, 1, ncd(iD))
                      !** do latitude
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%latitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%longitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%l2gpValue(1, kP, i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                      ENDDO
                      
                      !*** Interpolate to l2gp grid 
		      !CALL Residual2L2Grid( iD, kP, residD_temp, residD ) 
		      CALL Residual2L2Grid( iD, iP, kP, residD_temp, residD ) 
                      
                      if ( associated(sortTemp) ) then 
                         DeAllocate(sortTemp, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  sortTemp array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif
                   
                   !*** Calculate Maximum Difference
                   
                   if (ncd(iD) .gt. 0) then 
                      
                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(-abs(residD_temp(iD)% & 
                           l2gpValue(1, kP, :)), 1, ncd(iD), pt)
                      DO i = 1,cfDef%N
                         dmD(iD)%maxDiff(i, iP) = & 
                              residD_temp(iD)%l2gpValue(1, kP, pt(i))
                      ENDDO
                      
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif

                   !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                   avg = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                         avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                   ENDDO
                   IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
                   ENDIF
                   dmD(iD)%gRss(iP) = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                         dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                   ENDDO
                   IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmD(iD)%gRss(iP) = sqrt(dmD(iD)%gRss(iP)/real(dmD(iD)%nLons*cfProd%nLats-1))
                   ENDIF
                      dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
                   ENDDO
                ELSE IF (cfProd%mode == 'all') THEN
                   DO iD = 1, rDays
                      
                      if (nca(iD) .gt. 0) then 
                         
                         ALLOCATE( pt(nca(iD)), STAT=error )
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                     ENDIF
                     
                     ALLOCATE( sortTemp(nca(iD)), STAT=error )
                     IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                     ENDIF
                     
                     CALL my_sortp(residA_temp(iD)%time, 1, nca(iD), pt)
                     !** do time
                     CALL DSORT(residA_temp(iD)%time, 1, nca(iD))
                     !** do latitude
                     DO i = 1, nca(iD)
                        sortTemp(i) = residA_temp(iD)%latitude(i)
                     ENDDO
                     DO i = 1, nca(iD)
                        residA_temp(iD)%latitude(i) = sortTemp(pt(i))
                     ENDDO
                     !** do longitude
                     DO i = 1, nca(iD)
                        sortTemp(i) = residA_temp(iD)%longitude(i)
                     ENDDO
                     DO i = 1, nca(iD)
                        residA_temp(iD)%longitude(i) = sortTemp(pt(i))
                     ENDDO
                     !** do value
                     DO i = 1, nca(iD)
                        sortTemp(i) = residA_temp(iD)%l2gpValue(1, kP, i)
                     ENDDO
                     DO i = 1, nca(iD)
                        residA_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                     ENDDO
                     
                     !*** Interpolate to l2gp grid 
                     !CALL Residual2L2Grid( iD, kP, residA_temp, residA ) 
                     CALL Residual2L2Grid( iD, iP, kP, residA_temp, residA ) 
                     
                     if ( associated(sortTemp)) then 
                        
                        DeAllocate(sortTemp, STAT=error)
                        IF ( error /= 0 ) THEN
                           msr = MLSMSG_DeAllocate // '  sortTemp array.'
                           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                        ENDIF
                        
                     endif
                     
                     if ( associated(pt) ) then 
                        
                        DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                         
                      endif
                      
                   endif
                   
                   !*** Calculate Maximum Difference

                   if (nca(iD) .gt. 0) then 
                      
                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL my_sortp(-abs(residA_temp(iD)% & 
                           l2gpValue(1, kP, :)), 1, nca(iD), pt)
                      DO i = 1,cfDef%N
                         dmA(iD)%maxDiff(i, iP) = & 
                              residA_temp(iD)%l2gpValue(1, kP, pt(i))
                      ENDDO
                      
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif
                   
                   !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                   avg = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                         avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                   ENDDO
                   IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
                   ENDIF
                   dmA(iD)%gRss(iP) = 0.0
                   DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                         dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                              & (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                   ENDDO
                   IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmA(iD)%gRss(iP) = sqrt(dmA(iD)%gRss(iP)/real(dmA(iD)%nLons*cfProd%nLats-1))
                   ENDIF
                   dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
                
                DO iD = 1, rDays
                   
                   if (ncd(iD) .gt. 0) then 
                     
                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      ALLOCATE( sortTemp(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' sortTemp array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(residD_temp(iD)%time, 1, ncd(iD), pt)
                      !** do time
                      CALL DSORT(residD_temp(iD)%time, 1, ncd(iD))
                     !** do latitude
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%latitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%longitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, ncd(iD)
                         sortTemp(i) = residD_temp(iD)%l2gpValue(1, kP, i)
                      ENDDO
                      DO i = 1, ncd(iD)
                         residD_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                      ENDDO
                      
                      !*** Interpolate to l2gp grid 
		      !CALL Residual2L2Grid( iD, kP, residD_temp, residD ) 
		      CALL Residual2L2Grid( iD, iP, kP, residD_temp, residD ) 
                      
                      if ( associated(sortTemp) ) then
                         DeAllocate(sortTemp, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  sortTemp array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                   endif
                   
                   !*** Calculate Maximum Difference
                   
                   if ( ncd(iD) .gt. 0 ) then 
                      
                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                     CALL my_sortp(-abs(residD_temp(iD)% & 
                          & l2gpValue(1, kP, :)), 1, ncd(iD), pt)
                     DO i = 1,cfDef%N
                        dmD(iD)%maxDiff(i, iP) = & 
                             & residD_temp(iD)%l2gpValue(1, kP, pt(i))
                     ENDDO
                     
                     if ( associated(pt) ) then 
                        DeAllocate(pt, STAT=error)
                        IF ( error /= 0 ) THEN
                           msr = MLSMSG_DeAllocate // '  pt array.'
                           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                        ENDIF
                     endif
                     
                  endif
                  
                  !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                  avg = 0.0
                  DO J = 1, cfProd%nLats
                     DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                     ENDDO
                  ENDDO
                  IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                     avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
                  ENDIF
                  dmD(iD)%gRss(iP) = 0.0
                  DO J = 1, cfProd%nLats
                     DO I = 1, dmD(iD)%nLons
                        dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                             & (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                             & (dmD(iD)%l3dmValue(iP, J, I)-avg)
                     ENDDO
                  ENDDO
                  IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                     dmD(iD)%gRss(iP) = sqrt(dmD(iD)%gRss(iP)/real(dmD(iD)%nLons*cfProd%nLats-1))
                  ENDIF
                  dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
               ENDDO
               
               DO iD = 1, rDays
                  
                  if ( nc(iD) .gt. 0 ) then 
                     
                     ALLOCATE( pt(nc(iD)), STAT=error )
                     IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      ALLOCATE( sortTemp(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' sortTemp array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      
                      CALL my_sortp(l3r_temp(iD)%time, 1, nc(iD), pt)
                      !** do time
                      CALL DSORT(l3r_temp(iD)%time, 1, nc(iD))
                      !** do latitude
                      DO i = 1, nc(iD)
                         sortTemp(i) = l3r_temp(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                         l3r_temp(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nc(iD)
                         sortTemp(i) = l3r_temp(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                         l3r_temp(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nc(iD)
                         sortTemp(i) = l3r_temp(iD)%l2gpValue(1, kP, i)
                      ENDDO
                      DO i = 1, nc(iD)
                         l3r_temp(iD)%l2gpValue(1, kP, i) = sortTemp(pt(i))
                      ENDDO
                     
                      !*** Interpolate to l2gp grid
                      !CALL Residual2L2Grid( iD, kP, l3r_temp, l3r )
                      CALL Residual2L2Grid( iD, iP, kP, l3r_temp, l3r )
 
                      if ( associated(sortTemp) ) then 
                         DeAllocate(sortTemp, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  sortTemp array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                      if ( associated(pt) ) then 
                         DeAllocate(pt, STAT=error)
                         IF ( error /= 0 ) THEN
                            msr = MLSMSG_DeAllocate // '  pt array.'
                            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                         ENDIF
                      endif
                      
                  endif
                   
                  !*** Calculate Maximum Difference
                   
                  if ( nc(iD) .gt. 0 ) then 
                      
                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                         msr = MLSMSG_Allocate // ' pt array.'
                         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                     
                      CALL my_sortp(-abs(l3r_temp(iD)% & 
                           l2gpValue(1, kP, :)), 1, nc(iD), pt)
                     DO i = 1,cfDef%N
                        l3dm(iD)%maxDiff(i, iP) = & 
                             & l3r_temp(iD)%l2gpValue(1, kP, pt(i))
                     ENDDO
                     
                     if ( associated(pt) ) then 
                        DeAllocate(pt, STAT=error)
                        IF ( error /= 0 ) THEN
                           msr = MLSMSG_DeAllocate // '  pt array.'
                           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                        ENDIF
                     endif
                     
                  
                  !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                  avg = 0.0
                  DO J = 1, cfProd%nLats
                         DO I = 1, l3dm(iD)%nLons
                            avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                         ENDDO
                  ENDDO
                  IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                         avg = avg/real(l3dm(iD)%nLons*cfProd%nLats)
                  ENDIF
                  l3dm(iD)%gRss(iP) = 0.0
                  DO J = 1, cfProd%nLats
                         DO I = 1, l3dm(iD)%nLons
                            l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP) +  &
                                 & (l3dm(iD)%l3dmValue(iP, J, I)-avg)* &
                                 & (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                         ENDDO
                  ENDDO
                  IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                         l3dm(iD)%gRss(iP) = sqrt(l3dm(iD)%gRss(iP)/real(l3dm(iD)%nLons*cfProd%nLats-1))
                  ENDIF
                  l3dm(iD)%perMisPoints(iP) = perMisPoints(iP)

                endif

               ENDDO

            END IF
                
       ENDDO
       
       !*** Calculate Field Precisions ********
             
       iP = 0
       DO kP = pStartIndex, pStartIndex 
           iP = iP + 1
                
           DO J = 1, cfProd%nLats
              IF( anlats(J, iP) > 0 ) THEN
                      
                 IF( atimes(J, 1, iP) < dtimes(J, 1, iP) ) THEN
                         
                         lonD0_in = FindRealLon(real(dlons(J, 1, iP)))
                         lonA0_in = FindRealLon(real(alons(J, 1, iP)))
                         
                         CALL Init( & 
                              & nt_a_i   = anlats(J, iP), 		&
                              & nt_d_i   = dnlats(J, iP), 		&
                              & tau0_i   = tau0, 			&
                              & delTad_i = delTad(J, iP), 		&
                              & c0_i     = 2.0*PI, 			&
                              & lonD0_i  = lonD0_in, 			&
                              & tD0_i    = dtimes(J, 1, iP), 		&
                              & lonA0_i  = lonA0_in, 			&
                              & tA0_i    = atimes(J, 1, iP), 		&
                              & lat_i    = alats(J, 1, iP) )

                         CALL DataGenerate(aprec(J, :, iP), dprec(J, :, iP) )
                         
                 ELSE

                         lonD0_in = FindRealLon(real(dlons(J, 2, iP)))
                         lonA0_in = FindRealLon(real(alons(J, 1, iP)))

                         CALL Init( & 
                              & nt_a_i   = anlats(J, iP), 		&
                              & nt_d_i   = dnlats(J, iP)-1, 		&
                              & tau0_i   = tau0, 			&
                              & delTad_i = delTad(J, iP), 		&
                              & c0_i     = 2.0*PI, 			&
                              & lonD0_i  = lonD0_in, 			&
                              & tD0_i    = dtimes(J, 2, iP), 		&
                              & lonA0_i  = lonA0_in, 			&
                              & tA0_i    = atimes(J, 1, iP), 		&
                              & lat_i    = alats(J, 1, iP) )
                         CALL DataGenerate(aprec(J, :, iP), dprec(J, 2:, iP) )
                         
                 END IF
                      
                 IF (cfProd%mode == 'com') THEN

                         if (l3dm(1)%nLons .gt. 0) then 
                       
                            ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                            IF ( error /= 0 ) THEN
                               msr = MLSMSG_Allocate // ' l3Result array.'
                               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                            ENDIF
                            
                            CALL CordTransform(cfProd%mode)
                            CALL FFSM(l3spPrec(1), iP, J)
                            DO iD = 1, cfProd%nDays
                               CALL Reconstruct(cfProd%mode, &
                                    & real(l3dm(iD)%time - & 
                                    & l2gp(1)%time(1))/86400.0,&
				    & l3dm(iD)%nLons, l3dm(iD)%longitude, & 
                                    & l3Result)
                               DO I = 1, l3dm(iD)%nLons
                                  l3dm(iD)%l3dmPrecision(iP, J, I) = & 
                                       & l3Result(I) 
                               ENDDO
                            ENDDO
                            
                            if ( associated(l3Result) ) then 
                               DeAllocate(l3Result, STAT=error)
                               IF ( error /= 0 ) THEN
                                  msr = MLSMSG_DeAllocate //'  l3Result array.'
                                  CALL MLSMessage(MLSMSG_Error,ModuleName, msr)
                               ENDIF
                            endif
                         endif  ! nlons > 0
                    
                 ELSE IF (cfProd%mode == 'asc') THEN
                    
                    if (dmA(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMA(l3spPrec(1), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct(cfProd%mode, &
                               & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                          DO I = 1, dmA(iD)%nLons
                             dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                 ELSE IF (cfProd%mode == 'des') THEN
                    
                    if (dmD(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMD(l3spPrec(1), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct(cfProd%mode, &
                               & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                          DO I = 1, dmD(iD)%nLons
                             dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                 ELSE IF (cfProd%mode == 'ado') THEN
                    
                    if (dmA(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMA(l3spPrec(2), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct('asc', &
                               & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                          DO I = 1, dmA(iD)%nLons
                             dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                    if (dmD(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMD(l3spPrec(3), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct('des', &
                               & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                          DO I = 1, dmD(iD)%nLons
                             dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                 ELSE IF (cfProd%mode == 'all') THEN
                    
                    if (l3dm(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // '  l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSM(l3spPrec(1), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct('com', &
                               & real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                          DO I = 1, l3dm(iD)%nLons
                             l3dm(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                    if (dmA(1)%nLons .gt. 0) then 
                       
                       ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMA(l3spPrec(2), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct('asc', &
                               & real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                          DO I = 1, dmA(iD)%nLons
                             dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                    if (dmD(1)%nLons .gt. 0) then 
                      
                       ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                       IF ( error /= 0 ) THEN
                          msr = MLSMSG_Allocate // ' l3Result array.'
                          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                       ENDIF
                       CALL CordTransform(cfProd%mode)
                       CALL FFSMD(l3spPrec(3), iP, J)
                       DO iD = 1, cfProd%nDays
                          CALL Reconstruct('des', &
                               & real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                               & dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                          DO I = 1, dmD(iD)%nLons
                             dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
                          ENDDO
                       ENDDO
                       
                       if ( associated(l3Result) ) then 
                          DeAllocate(l3Result, STAT=error)
                          IF ( error /= 0 ) THEN
                             msr = MLSMSG_DeAllocate // '  l3Result array.'
                             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                          ENDIF
                       endif
                       
                    endif
                    
                 END IF
                
                 CALL ClearMemory()
                 
              END IF
           ENDDO
           
       ENDDO

       !*** Finish and Clean the Memory ********

       if ( associated(l3spPrec) ) then 
          DeAllocate(l3spPrec, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  l3spPrec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(nc) ) then
          DeAllocate(nc, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  nc array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(nca) ) then 
          DeAllocate(nca, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  nca array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(ncd) ) then
          DeAllocate(ncd, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  ncd array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(startTime)) then 
          DeAllocate(startTime, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  startTime array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(endTime)) then
          DeAllocate(endTime, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  endTime array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       if ( associated(alats) ) then
          DeAllocate(alats, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  alats array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(dlats) ) then
          DeAllocate(dlats, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  dlats array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(alons)) then
          DeAllocate(alons, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  alons array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(dlons)) then
          DeAllocate(dlons, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  dlons array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(atimes)) then
          DeAllocate(atimes, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  atimes array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(dtimes) ) then
          DeAllocate(dtimes, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  dtimes array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(afields) ) then
          DeAllocate(afields, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  afields array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(dfields) ) then
          DeAllocate(dfields, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  dfields array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
                                                                                           
       if ( associated(aprec)) then
          DeAllocate(aprec, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  aprec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       if ( associated(dprec)) then
          DeAllocate(dprec, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  dprec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
       endif
       
       CALL DestroyL2GPDatabase(l2gp)
       IF (cfProd%mode == 'com') THEN
          CALL DestroyL2GPDatabase(l3r_temp)
       ELSE IF (cfProd%mode == 'asc') THEN
          CALL DestroyL2GPDatabase(residA_temp)
       ELSE IF (cfProd%mode == 'dsc') THEN
          CALL DestroyL2GPDatabase(residD_temp)
       ELSE IF (cfProd%mode == 'ado') THEN
          CALL DestroyL2GPDatabase(residA_temp)
          CALL DestroyL2GPDatabase(residD_temp)
       ELSE IF (cfProd%mode == 'all') THEN
          CALL DestroyL2GPDatabase(l3r_temp)
          CALL DestroyL2GPDatabase(residA_temp)
          CALL DestroyL2GPDatabase(residD_temp)
       ENDIF  

       !-----------------------------------
     END SUBROUTINE DailyCoreProcessing
     !-----------------------------------

      !------------------------------------------------------------------------
     SUBROUTINE SortData(cfProd, l2Days, l2gp, pStartIndex, pEndIndex, tau0, & 
          & anlats, dnlats, alats_interp, dlats_interp, alons_interp, & 
          & dlons_interp, atimes_interp, dtimes_interp, afields_interp, & 
          & dfields_interp, aprec_interp, dprec_interp, delTad, perMisPoints)
      !------------------------------------------------------------------------
       USE SDPToolkit, ONLY: PI
       !USE MLSNumerics, ONLY: INTERPOLATEVALUES, InterpolateArraySetup, INTERPOLATEARRAYTEARDOWN
       USE MLSNumerics, ONLY: INTERPOLATEVALUES
       !USE MLSNumerics, ONLY: HUNT 
        
       TYPE( L3CFProd_T ) :: cfProd
       TYPE( L2GPData_T ), POINTER :: l2gp(:)
        
       CHARACTER (LEN=480) :: msr
        
       INTEGER :: pStartIndex, pEndIndex, nIndex
        
       REAL (r8), DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: delTad 
       
       REAL (r8), POINTER, DIMENSION(:, :, :) ::             & 
            & alons_interp(:, :, :), alats_interp(:, :, :),    & 
            & dlons_interp(:, :, :), dlats_interp(:, :, :),    & 
            & atimes_interp(:, :, :), dtimes_interp(:, :, :),  &
            & afields_interp(:, :, :), dfields_interp(:, :, :),& 
            & aprec_interp(:, :, :), dprec_interp(:, :, :) 
       
       REAL (r8), POINTER, DIMENSION(:, :) ::  &
            & l2Times(:, :), l2Lons(:, :), & 
            & l2Lons_old(:, :), l2Lats(:, :), l2Values(:, :), l2Prec(:, :), &
            & l2GeodAngle(:, :), dataInput(:), latInput(:), resultIntp(:), latIntp(:), &
            & resultTimesIntp(:), resultLonsIntp(:), resultPrecIntp(:)

       REAL (r8), POINTER, DIMENSION(:) :: l2Qual(:) 
       Integer, POINTER, DIMENSION(:) :: l2Status(:) 
       
       REAL (r8) &
            & dlons, alons, lons_found, lats_found, times_found, fields_found,&
            & prec_found, slope, slope_time, slope_field, slope_prec, sTime, & 
            & lons_found_old
        
       REAL :: tau0, latRange, firstTime
        
       INTEGER, DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: & 
            & anlats, dnlats
       INTEGER, DIMENSION(pEndIndex-pStartIndex+1) :: & 
            & perMisPoints, numData
       
       INTEGER :: nPd, iMin, error, i, j, k, iT, iD, iP, kP, iL, & 
            & nterms, nstart, nr, lindex, lindex_prev, l2Days, nloop, aindex,& 
            & aindex_prev, dindex_prev, latIndex, prevLatIndex, nPIntp, nPTotal, &
            & prevlatOrder, prevlat, latOrder, orbitStatus, latStartIndex, latEndIndex
       
       REAL (r8), DIMENSION(37, 3496) :: read_prec 
       REAL (r8), DIMENSION(37) :: read_pressure 

       REAL (r8), DIMENSION(83) :: timeDiff_asc_to_des, timeDiff_des_to_asc 

       Logical bfound, binaction
       INTEGER npMissing, npIndex, maxDays

       REAL (r8), POINTER, DIMENSION(:) :: ioTimeArray, ioValueArray, ioLatArray, ioLonArray, &
            &                              ioGeodArray, ioPrecArray, ioLonNewArray

       !*** Initilize the timeDiff_asc_to_des, timeDiff_des_to_asc
 
       latRange = 8.0

       timeDiff_des_to_asc = (/ 0.0, &
                                190.156,  290.521,  376.161,  455.367,  531.078,  604.537,  676.494,  747.424,  817.580,  887.059, &
                                956.125, 1024.869, 1093.288, 1161.432, 1229.357, 1297.078, 1364.655, 1432.178, 1499.502, 1566.795, &
                               1633.957, 1701.070, 1768.077, 1835.055, 1902.000, 1968.812, 2035.705, 2102.527, 2169.208, 2235.945, &
                               2302.584, 2369.311, 2435.913, 2502.582, 2569.186, 2635.781, 2702.395, 2768.886, 2835.555, 2902.073, &
                               2968.733, 3035.184, 3101.706, 3168.211, 3234.851, 3301.390, 3367.921, 3434.550, 3501.134, 3567.662, &
                               3634.334, 3700.963, 3767.612, 3834.309, 3901.051, 3967.798, 4034.556, 4101.489, 4168.244, 4235.218, &
                               4302.215, 4369.293, 4436.466, 4503.766, 4571.024, 4638.574, 4706.177, 4774.006, 4842.007, 4910.293, &
                               4978.839, 5047.746, 5117.202, 5187.135, 5257.944, 5329.738, 5403.078, 5478.616, 5557.638, 5643.114, &
                               5743.230, &
                                0.0 /)
       timeDiff_asc_to_des = (/ 0.0, &
                               5742.843, 5642.495, 5556.851, 5477.603, 5401.936, 5328.448, 5256.506, 5185.599, 5115.440, 5045.933, &
                               4976.886, 4908.115, 4839.720, 4771.568, 4703.624, 4635.913, 4568.349, 4500.825, 4433.535, 4366.219, &
                               4299.075, 4231.939, 4164.947, 4097.938, 4031.027, 3964.192, 3897.299, 3830.525, 3763.786, 3697.070, &
                               3630.410, 3563.683, 3497.089, 3430.443, 3363.846, 3297.237, 3230.593, 3164.122, 3097.461, 3030.924, &
                               2964.288, 2897.837, 2831.291, 2764.801, 2698.158, 2631.620, 2565.056, 2498.508, 2431.901, 2365.329, &
                               2298.671, 2232.051, 2165.394, 2098.696, 2031.967, 1965.198, 1898.439, 1831.531, 1764.703, 1697.760, &
                               1630.782, 1563.703, 1496.558, 1429.211, 1361.990, 1294.430, 1226.814, 1159.022, 1090.977, 1022.685, &
                                954.157,  885.204,  815.793,  745.821,  675.002,  603.197,  529.881,  454.363,  375.302,  289.864, &
                                189.772, &
                                0.0 /)

       !*** Calculate the number of points in each pressure level 
       
       nPIntp = 10
       nPd = 15
       maxDays = 31  ! need to allocate more in case of missing days
 
       ! testing the missing day for Oct 2005 data: l2Days + 1
       !l2Days = l2Days + 1  ! testing
 
       nterms = 0
       DO I = 1, l2Days
          nterms = nterms + l2gp(I)%nTimes
       ENDDO
       
       !*** Allocate space for all the arrays
          
       ALLOCATE(latIntp(1), STAT=error)
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' latIntp array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
          
       ALLOCATE(resultIntp(1), STAT=error)
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' resultIntp array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
          
       ALLOCATE(resultTimesIntp(1), STAT=error)
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' resultIntp array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
          
       ALLOCATE(resultLonsIntp(1), STAT=error)
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' resultIntp array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
          
       ALLOCATE(resultPrecIntp(1), STAT=error)
       IF ( error /= 0 ) THEN
          msr = MLSMSG_Allocate // ' resultIntp array.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
          
       
       if ((nterms.gt.0).and.( (pEndIndex-pStartIndex).ge.0) ) then 
          
          ALLOCATE(l2Times(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Times array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(l2GeodAngle(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2GeodAngle array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          
          ALLOCATE(l2Lons(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lons array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          
          ALLOCATE(l2Lats(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lats array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          
          ALLOCATE(l2Values(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Values array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(l2Prec(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Prec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(l2Qual(nterms), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Qual array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(l2Status(nterms), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Status array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
       endif
       
       if ( (cfProd%nLats.gt.0).and. & 
            & ( (pEndIndex-pStartIndex).ge.0).and.(nPd*l2Days.gt.0)  ) then 
          
          ALLOCATE(alons_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' alons_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(alats_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' alats_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(atimes_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' atimes_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(afields_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' afields_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(aprec_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' aprec_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(dlons_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dlons_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(dlats_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dlats_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(dtimes_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dtimes_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(dfields_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dfields_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
          ALLOCATE(dprec_interp(cfProd%nLats, nPd*maxDays+1, & 
               pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dprec_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
          
       endif
 
       !*** Re-arrange the data into longitude order for each pressure level 
       
       DO kP = pStartIndex, pEndIndex
          perMisPoints(kP+1-pStartIndex) = 0
          numData(kP+1-pStartIndex) = 0
       ENDDO

!***********************************************************************
       
       nstart = 1
       DO iD = 1, l2Days
	  DO iT = 1, l2gp(iD)%nTimes
             iP = 0
             DO kP = pStartIndex, pEndIndex 
	        iP = iP + 1
                l2Times(nstart, iP)     = l2gp(iD)%time(iT)
                l2GeodAngle(nstart, iP) = l2gp(iD)%geodAngle(iT)
                !l2Lons(nstart, iP)     = l2gp(iD)%longitude(iT) * PI/180.0
                !l2Lons_old(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lons(nstart, iP)      = l2gp(iD)%longitude(iT)
                l2Lats(nstart, iP)      = l2gp(iD)%latitude(iT)
                l2Values(nstart, iP)    = l2gp(iD)%l2gpValue(1, kP, iT) 
		l2Prec(nstart, iP)      = l2gp(iD)%l2gpPrecision(1, kP, iT)
		l2Qual(nstart)          = l2gp(iD)%Quality(iT)
		l2Status(nstart)        = l2gp(iD)%Status(iT)
                numData(iP)             = numData(iP) + 1
                !IF(l2gp(iD)%l2gpValue(1, kP, iT) == -999.99) THEN 
                IF(isGoodData(l2Values(nstart, iP), l2Qual(nstart), l2Status(nstart)) ) THEN
                   perMisPoints(iP)  = perMisPoints(iP) + 1     ! Good Points
                ENDIF
             ENDDO
             nstart = nstart + 1
	  ENDDO
       ENDDO

       DO kP = pStartIndex, pEndIndex
          if (numData(kP+1-pStartIndex).ne.0) then 
             perMisPoints(kP+1-pStartIndex) = 100.0-100.0*perMisPoints(kP+1-pStartIndex)/(l2Days*3494)
                  !&100*perMisPoints(kP+1-pStartIndex)/numData(kP+1-pStartIndex)
          endif
       ENDDO
       
       !*** Prepare Data for L3 Processing (Convert time reference to the starting point)
       
       iP = 0
       DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
          !sTime = l2Times(1, iP)
          DO k = 1, nterms
             if(l2Times(k, iP) > 0) then
               sTime = l2Times(k, iP)
               exit
             endif
          ENDDO
          
          anlats(1, iP) = 0 
          anlats(cfProd%nLats, iP) = 0 
          dnlats(1, iP) = 0 
          dnlats(cfProd%nLats, iP) = 0 

          DO J = 2, cfProd%nLats-1
            anlats(J, iP) = 0 
            dnlats(J, iP) = 0 

            iL = 0
            binaction = .false.
	    DO iT = 1, nterms 
              IF( abs(l2Lats(iT, iP)-cfProd%latGridMap(J)) <= latRange*0.5 ) THEN 
	        iL = iL + 1
                if(iL == 1) then
                  binaction = .true.
                  orbitStatus = isAscending(l2GeodAngle(iT, iP))
                  firstTime = l2Times(iT, iP)
                end if
     
                if(isAscending(l2GeodAngle(iT, iP)) == orbitStatus &
                   & .and. isGoodData(l2Values(iT, iP), l2Qual(iT), l2Status(iT)) & 
                   !& .and. l2Values(iT, iP) /= -999.99 & 
                   & .and. (l2Times(iT, iP)-firstTime) < 250.0 & 
                   & .and. binaction) then
                  Call ExpandArray(ioValueArray, l2Values(iT, iP))
                  Call ExpandArray(ioPrecArray, l2Prec(iT, iP))
                  Call ExpandArray(ioTimeArray, l2Times(iT, iP))
                  Call ExpandArray(ioLatArray, l2Lats(iT, iP))
                  Call ExpandArray(ioLonArray, l2Lons(iT, iP))
                  !Call ExpandArray(ioGeodArray, l2GeodAngle(iT, iP))
                else if(isAscending(l2GeodAngle(iT, iP)) /= orbitStatus &
                   & .or. (l2Times(iT, iP)-firstTime) > 250.0 ) then
                  iL = 0
                  binaction = .false.
                  latIntp(1) = cfProd%latGridMap(J)

                  if(associated(ioLatArray) .and. size(ioLatArray) > 1) then

                    ALLOCATE(ioLonNewArray(size(ioLonArray)), STAT=error)
                    IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' ioLonNewArray array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF

                    nr = 0
                    ioLonNewArray(1) = ioLonArray(1)
	            DO i = 2, size(ioLonArray)
                      if(ioLonArray(i) >= 0.0 .and. ioLonArray(i-1) < 0.0) then
                         ioLonNewArray(i) = ioLonArray(i) - nr*360.0 - 360.0
                         nr = nr + 1
                      else
                         ioLonNewArray(i) = ioLonArray(i) - nr*360.0
                      endif
                    enddo 
	            DO i = 1, size(ioLonArray)
                    enddo

               !write(*,'(10(1pE11.3))') (ioValueArray(i), i=1, size(ioValueArray))
                    call InterpolateValues ( &
                                           & ioLatArray, &                  ! OldX
                                           & ioValueArray, &                ! OldY
                                           & latIntp, &                     ! NewX
                                           & resultIntp, &                  ! NewY
                                           & 'Linear', &                    ! use linear
                                           & extrapolate='Constant' )       ! No extrapolation
               !write(*,'(10(1pE11.3))') resultIntp(1) 
                    call InterpolateValues ( &
                                           & ioLatArray, &                  ! OldX
                                           & ioPrecArray, &                ! OldY
                                           & latIntp, &                     ! NewX
                                           & resultPrecIntp, &                  ! NewY
                                           & 'Linear', &                    ! use linear
                                           & extrapolate='Constant' )       ! No extrapolation
                    call InterpolateValues ( &
                                           & ioLatArray, &                  ! OldX
                                           & ioLonNewArray, &                ! OldY
                                           & latIntp, &                     ! NewX
                                           & resultLonsIntp, &                  ! NewY
                                           & 'Linear', &                    ! use linear
                                           & extrapolate='Constant' )       ! No extrapolation
                    if ( associated(ioLonNewArray) ) then
                      DEALLOCATE(ioLonNewArray, STAT=error)
                      IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioLonNewArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                    endif
                    call InterpolateValues ( &
                                           & ioLatArray, &                  ! OldX
                                           & ioTimeArray, &                ! OldY
                                           & latIntp, &                     ! NewX
                                           & resultTimesIntp, &                  ! NewY
                                           & 'Linear', &                    ! use linear
                                           & extrapolate='Constant' )       ! No extrapolation

                    IF( orbitStatus == 1) THEN                    
                      anlats(J, iP) = anlats(J, iP) + 1 
                      alats_interp(J, anlats(J, iP), iP)   = cfProd%latGridMap(J) 
                      atimes_interp(J, anlats(J, iP), iP)  = resultTimesIntp(1)
                      afields_interp(J, anlats(J, iP), iP) = resultIntp(1)
                      aprec_interp(J, anlats(J, iP), iP)   = resultPrecIntp(1)
                      if(resultLonsIntp(1) < -180.0) then
                         alons_interp(J, anlats(J, iP), iP)   = resultLonsIntp(1) + 360.0
                      else
                         alons_interp(J, anlats(J, iP), iP)   = resultLonsIntp(1)
                      endif
                    END IF
          
                    IF( orbitStatus == 2) THEN
                      IF(anlats(J, iP) > 0) THEN
                        dnlats(J, iP) = dnlats(J, iP) + 1 
                        dlats_interp(J, dnlats(J, iP), iP)   = cfProd%latGridMap(J) 
                        dtimes_interp(J, dnlats(J, iP), iP)  = resultTimesIntp(1)
                        dfields_interp(J, dnlats(J, iP), iP) = resultIntp(1)
                        dprec_interp(J, dnlats(J, iP), iP)   = resultPrecIntp(1)
                        if(resultLonsIntp(1) < -180.0) then
                          dlons_interp(J, dnlats(J, iP), iP)   = resultLonsIntp(1) + 360.0
                        else
                          dlons_interp(J, dnlats(J, iP), iP)   = resultLonsIntp(1)
                        endif
                      END IF
                    END IF
              
                  endif

                  if ( associated(ioValueArray) ) then
                    DEALLOCATE(ioValueArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioValueArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
                  endif
                  if ( associated(ioTimeArray) ) then
                    DEALLOCATE(ioTimeArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioTimeArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
                  endif
                  if ( associated(ioLatArray) ) then
                    DEALLOCATE(ioLatArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioLatArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
                  endif
                  if ( associated(ioLonArray) ) then
                    DEALLOCATE(ioLonArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioLonArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
                  endif
                  if ( associated(ioPrecArray) ) then
                    DEALLOCATE(ioPrecArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioPrecArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
                  endif

                endif 

              ENDIF
            ENDDO

            ! Modify data for the missing data points
            !** Ascending
            npMissing = 0
            DO iT = 1, anlats(J, iP)-1
              !if(atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP) < 5930.0) then
                 !print *, anlats(J, iP), iT, atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP), int((atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP))/5800.0)
                !print *, anlats(J, iP), iT, atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP), int((atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP))/5800.0), &
                !  &       afields_interp(J, iT, iP)
              !endif
              if(atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP) > 6000.0) then
                npMissing = npMissing + int((atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP))/5800.0)-1
              endif
            ENDDO
            !print *, 'ascending missing=', iP, pEndIndex-pStartIndex, J, cfProd%latGridMap(J), npMissing

            ALLOCATE(ioTimeArray(anlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioTimeArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioValueArray(anlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioValueArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioPrecArray(anlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioPrecArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioLonArray(anlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioLonArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            npIndex = 1 
            DO iT = 1, anlats(J, iP)-1
              ioTimeArray(npIndex) = atimes_interp(J, iT, iP)
              ioValueArray(npIndex) = afields_interp(J, iT, iP)
              ioPrecArray(npIndex) = aprec_interp(J, iT, iP)
              ioLonArray(npIndex) = alons_interp(J, iT, iP)
              DO iL = 1, int((atimes_interp(J, iT+1, iP) - atimes_interp(J, iT, iP))/5800.0)-1 
                 npIndex = npIndex + 1
                 ioTimeArray(npIndex) = atimes_interp(J, iT, iP) + 5933.0*iL
                 ioValueArray(npIndex) = -999.99 
                 ioPrecArray(npIndex) = -999.99 
                 ioLonArray(npIndex) = alons_interp(J, iT, iP) - 24.72*iL
                 if(ioLonArray(npIndex) < -180.0) then
                   ioLonArray(npIndex) = 360.0 + ioLonArray(npIndex)
                 endif
              ENDDO
              npIndex = npIndex + 1
            ENDDO
            ioTimeArray(npIndex) = atimes_interp(J, anlats(J, iP), iP)
            ioValueArray(npIndex) = afields_interp(J, anlats(J, iP), iP)
            ioPrecArray(npIndex) = aprec_interp(J, anlats(J, iP), iP)
            ioLonArray(npIndex) = alons_interp(J, anlats(J, iP), iP)

            call InterpolateValues ( &
                                     & atimes_interp(J, 1:anlats(J, iP), iP), &                  ! OldX
                                     & afields_interp(J, 1:anlats(J, iP), iP), &                 ! OldY
                                     & ioTimeArray, &                                            ! NewX
                                     & ioValueArray, &                                           ! NewY
                                     & 'Linear', &                                               ! use linear
                                     & extrapolate='Constant' )                                  ! No extrapolation

            call InterpolateValues ( &
                                     & atimes_interp(J, 1:anlats(J, iP), iP), &                  ! OldX
                                     & aprec_interp(J, 1:anlats(J, iP), iP), &                 ! OldY
                                     & ioTimeArray, &                                            ! NewX
                                     & ioPrecArray, &                                           ! NewY
                                     & 'Linear', &                                               ! use linear
                                     & extrapolate='Constant' )                                  ! No extrapolation

            !if(ioValueArray(iT) < 0) then
            !    write(*,*) iT, afields_interp(J, iT, iP), ioValueArray(iT)
            !endif
            if(npIndex > anlats(J, iP)) then
               !write(*,'(10(1pE11.3))') (afields_interp(J, iT, iP), atimes_interp(J, iT, iP), iT=1, anlats(J, iP))
            endif
!            k = 1
!            DO iT = 1, anlats(J, iP)
!               if(afields_interp(J, iT, iP) < 0) then
!                 k = 0
!               endif
!            enddo
!            if(k == 0) then
!               !write(*,*) J, iP,  anlats(J, iP)
!               !write(*,'(10(1pE11.3))') (afields_interp(J, iT, iP), iT=1, anlats(J, iP))
!            endif
!            k = 1
            DO iT = 1, size(ioTimeArray) 
               afields_interp(J, iT, iP) = ioValueArray(iT)
!               if(ioValueArray(iT) < 0) then
!                 k = 0
!               endif
               aprec_interp(J, iT, iP) = ioPrecArray(iT)
               atimes_interp(J, iT, iP) = (ioTimeArray(iT)-sTime)/86400.0
               alats_interp(J, iT, iP) = cfProd%latGridMap(J) 
               alons_interp(J, iT, iP) = ioLonArray(iT)*PI/180.0 - 2.0*PI*(iT-1)
            ENDDO
            anlats(J, iP) = anlats(J, iP) + npMissing

!            if(k == 0) then
               !write(*,'(10(1pE11.3))') (afields_interp(J, iT, iP), atimes_interp(J, iT, iP), iT=1, anlats(J, iP))
               !write(*,*) J, iP,  anlats(J, iP)
               !write(*,'(10(1pE11.3))') (afields_interp(J, iT, iP), iT=1, anlats(J, iP))
!            endif

            if ( associated(ioValueArray) ) then
                    DEALLOCATE(ioValueArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioValueArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioTimeArray) ) then
                    DEALLOCATE(ioTimeArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioTimeArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioLatArray) ) then
                    DEALLOCATE(ioLatArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioatnArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioLonArray) ) then
                    DEALLOCATE(ioLonArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioLonArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioPrecArray) ) then
                    DEALLOCATE(ioPrecArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioPrecArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif

            !** Descending
            npMissing = 0
            DO iT = 1, dnlats(J, iP)-1
              if(dtimes_interp(J, iT+1, iP) - dtimes_interp(J, iT, iP) > 6000.0) then
                npMissing = npMissing + int((dtimes_interp(J, iT+1, iP) - dtimes_interp(J, iT, iP))/5800.0)-1
              endif
            ENDDO

            ALLOCATE(ioTimeArray(dnlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioTimeArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioValueArray(dnlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioValueArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioPrecArray(dnlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioPrecArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            ALLOCATE(ioLonArray(dnlats(J, iP)+npMissing), STAT=error)
            IF ( error /= 0 ) THEN
               msr = MLSMSG_Allocate // ' ioLonArray array.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            npIndex = 1 
            DO iT = 1, dnlats(J, iP)-1
              ioTimeArray(npIndex) = dtimes_interp(J, iT, iP)
              ioValueArray(npIndex) = dfields_interp(J, iT, iP)
              ioPrecArray(npIndex) = dprec_interp(J, iT, iP)
              ioLonArray(npIndex) = dlons_interp(J, iT, iP)
              DO iL = 1, int((dtimes_interp(J, iT+1, iP) - dtimes_interp(J, iT, iP))/5800.0)-1 
                 npIndex = npIndex + 1
                 ioTimeArray(npIndex) = dtimes_interp(J, iT, iP) + 5933.0*iL
                 ioValueArray(npIndex) = -999.99 
                 ioPrecArray(npIndex) = -999.99 
                 ioLonArray(npIndex) = dlons_interp(J, iT, iP) - 24.72*iL
                 if(ioLonArray(npIndex) < -180.0) then
                   ioLonArray(npIndex) = 360.0 + ioLonArray(npIndex)
                 endif
              ENDDO
              npIndex = npIndex + 1
            ENDDO
            ioTimeArray(npIndex) = dtimes_interp(J, dnlats(J, iP), iP)
            ioValueArray(npIndex) = dfields_interp(J, dnlats(J, iP), iP)
            ioPrecArray(npIndex) = dprec_interp(J, dnlats(J, iP), iP)
            ioLonArray(npIndex) = dlons_interp(J, dnlats(J, iP), iP)

            call InterpolateValues ( &
                                     & dtimes_interp(J, 1:dnlats(J, iP), iP), &                  ! OldX
                                     & dfields_interp(J, 1:dnlats(J, iP), iP), &                 ! OldY
                                     & ioTimeArray, &                                            ! NewX
                                     & ioValueArray, &                                           ! NewY
                                     & 'Linear', &                                               ! use linear
                                     & extrapolate='Constant' )                                  ! No extrapolation

            call InterpolateValues ( &
                                     & dtimes_interp(J, 1:dnlats(J, iP), iP), &                  ! OldX
                                     & dprec_interp(J, 1:dnlats(J, iP), iP), &                 ! OldY
                                     & ioTimeArray, &                                            ! NewX
                                     & ioPrecArray, &                                           ! NewY
                                     & 'Linear', &                                               ! use linear
                                     & extrapolate='Constant' )                                  ! No extrapolation

            DO iT = 1, size(ioTimeArray) 
               dfields_interp(J, iT, iP) = ioValueArray(iT)
               dprec_interp(J, iT, iP) = ioPrecArray(iT)
               dtimes_interp(J, iT, iP) = (ioTimeArray(iT)-sTime)/86400.0
               dlats_interp(J, iT, iP) = cfProd%latGridMap(J) 
               dlons_interp(J, iT, iP) = ioLonArray(iT)*PI/180.0 - 2.0*PI*(iT-1)
            ENDDO
            dnlats(J, iP) = dnlats(J, iP) + npMissing

            if ( associated(ioValueArray) ) then
                    DEALLOCATE(ioValueArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioValueArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioTimeArray) ) then
                    DEALLOCATE(ioTimeArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioTimeArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioLatArray) ) then
                    DEALLOCATE(ioLatArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioatnArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioLonArray) ) then
                    DEALLOCATE(ioLonArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioLonArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif
            if ( associated(ioPrecArray) ) then
                    DEALLOCATE(ioPrecArray, STAT=error)
                    IF ( error /= 0 ) THEN
                       msr = MLSMSG_DeAllocate // '  ioPrecArray pointer.'
                       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                    ENDIF
            endif

          ENDDO

       ENDDO
       
       !*** Find averaged asc/des time difference 
       
       iP = 0
       DO kP = pStartIndex, pEndIndex 
          iP = iP + 1
          DO J = 1, cfProd%nLats
             delTad(J,iP) = 0.0
             IF (anlats(J, iP) < dnlats(J, iP)) THEN
                iMin = anlats(J, iP)
             ELSE
                iMin = dnlats(J, iP)
             END IF
             DO I = 2, iMin 
                delTad(J,iP) = & 
                     & delTad(J,iP) + dtimes_interp(J, I, iP) & 
                     & - atimes_interp(J, I, iP)
        !write(1, *) iP, J, I, dtimes_interp(J, I, iP), atimes_interp(J, I, iP) 
        !write(*, *) iP, J, I, dtimes_interp(J, I, iP), atimes_interp(J, I, iP) , tau0
             ENDDO
             IF(iMin > 0) THEN
                if (iMin .ne. 1) delTad(J,iP) = delTad(J,iP)/float(iMin-1)
             ELSE
                delTad(J,iP) = 0.0 
             END IF
             
             IF(delTad(J,iP) < 0) THEN
                delTad(J,iP) = delTad(J,iP) + tau0 
             END IF
        !write(*, *) iP, J, delTad(J,iP), tau0
          ENDDO
       ENDDO
       
       !*** Deallocate intermediate arrays 
       
       if ((nterms.gt.0).and.( (pEndIndex-pStartIndex).ge.0) ) then 
          
          if ( associated(l2Times) ) then
             DEALLOCATE(l2Times, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Times pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2GeodAngle) ) then
             DEALLOCATE(l2GeodAngle, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2GeodAngle pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Lons) ) then 
             DEALLOCATE(l2Lons, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Lons pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Lats) ) then
             DEALLOCATE(l2Lats, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Lats pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Values) ) then
             DEALLOCATE(l2Values, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Values pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Prec) ) then
             DEALLOCATE(l2Prec, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Prec pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Qual) ) then
             DEALLOCATE(l2Qual, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Qual pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
          if ( associated(l2Status) ) then
             DEALLOCATE(l2Status, STAT=error)
             IF ( error /= 0 ) THEN
                msr = MLSMSG_DeAllocate // '  l2Status pointer.'
                CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
             ENDIF
          endif
          
       endif
 
     !-----------------------------------
     END SUBROUTINE SortData
     !-----------------------------------


     !-------------------------------------------------------------------------
     Subroutine ExpandArray(ioArray, value)
     !-------------------------------------------------------------------------

     ! This routine adds to the ioArray arrays 

     ! Dummy arguments
     REAL (r8), POINTER, DIMENSION(:) :: ioArray 
     Real (r8), intent(in) :: value

     ! Local variables
     REAL (r8), POINTER, DIMENSION(:) :: tempArray => null()
     integer :: newsize
     integer :: status
     if ( associated(ioArray) ) then
        newSize=SIZE(ioArray)+1
     else
        newSize=1
     end if
     if(newSize == 1) then
        allocate(ioArray(newSize),STAT=status)
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "ioArray")
        ioArray(newSize) = value
        return
     endif
     allocate(tempArray(newSize),STAT=status)
     if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "temp ioArray")
     if ( newSize>1 ) tempArray(1:newSize-1) = ioArray 
     if ( ASSOCIATED(ioArray) ) &
      & deallocate ( ioArray, stat=status )
     if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "ioArray" )
     ioArray => tempArray
     ioArray(newSize) = value

      !-----------------------------------
     END Subroutine ExpandArray
     !-----------------------------------

     !-------------------------------------------------------------------------
     Logical FUNCTION isGoodData(value, qual, status)
     !-------------------------------------------------------------------------
  
       Real(r8) value, qual
       Integer status

       IF(value == -999.99 .or. qual < 0.0 .or. status < 0 .or. mod(status, 2) == 1) THEN
          isGoodData = .false. 
       ELSE 
          isGoodData = .true. 
       END IF

      !-----------------------------------
     END FUNCTION isGoodData
     !-----------------------------------

     !-------------------------------------------------------------------------
     Integer FUNCTION isAscending(geodAngle)
     !-------------------------------------------------------------------------
  
       Real(r8) modGeod, geodAngle

       IF(geodAngle > -360.0) THEN
          modGeod = mod(geodAngle, 360.0_r8)
          IF(modGeod < 0.0) THEN
            modGeod = modGeod + 360.0
          ENd IF
          IF(modGeod < 90.0 .OR. modGeod >= 270.0) THEN
            isAscending = 1
          ELSE IF(modGeod >= 90.0 .AND. modGeod < 270.0) THEN
            isAscending = 2
          ELSE 
            isAscending = 0 
          END IF
       ELSE 
          isAscending = 0 
       END IF

      !-----------------------------------
     END FUNCTION isAscending
     !-----------------------------------


     !-------------------------------------------------------------------------
     REAL(r8) FUNCTION FindRealLon(alon)
     !-------------------------------------------------------------------------
  USE SDPToolkit, ONLY: PI
  
       Real, intent(in) :: alon
  
       IF(mod(alon, 2.0*PI) < -PI) THEN
          FindRealLon = 2.0*PI + mod(alon, 2.0*PI) 
       ELSE IF(mod(alon, 2.0*PI) > PI) THEN
          FindRealLon = -2.0*PI + mod(alon, 2.0*PI) 
       ELSE 
          FindRealLon = mod(alon, 2.0*PI) 
       END IF
       
      !-----------------------------------
     END FUNCTION FindRealLon
     !-----------------------------------

     
     !-------------------------------------------------------------------------
     SUBROUTINE Residual2L2Grid( iD, iP, iN, l3r_temp, l3r ) 
     !-------------------------------------------------------------------------
  
       TYPE( L2GPData_T ), POINTER :: l3r(:)
       TYPE( L2GPData_T ), POINTER :: l3r_temp(:)
       INTEGER, INTENT(IN) :: iD, iP, iN
       
!       integer :: error, i, j, iT, icount, l2Days, nloop, aindex, & 
!            aindex_prev, dindex_prev

       INTEGER :: i, j, icount, startTime, endTime
       
       DO i = 1, l3r_temp(iD)%nTimes
         if(l3r_temp(iD)%time(i) > 0) then
            startTime = l3r_temp(iD)%time(i)
            exit
         endif
       ENDDO
       
       DO i = l3r_temp(iD)%nTimes, 1, -1
         if(l3r_temp(iD)%time(i) > 0) then
            endTime = l3r_temp(iD)%time(i)
            exit
         endif
       ENDDO

       icount = 1
       DO I = 1, l3r(iD)%nTimes 
          DO J = icount, l3r_temp(iD)%nTimes 
             If( l3r(iD)%time(I) < startTime .OR. & 
                  !& l3r(iD)%time(I) > l3r_temp(iD)%time(l3r(iD)%nTimes) ) Then
                  & l3r(iD)%time(I) > endTime ) Then
                l3r(iD)%l2gpValue(1,iP, I) = -999.99 
                EXIT
             ELSE IF (J .LT. l3r_temp(iD)%nTimes) THEN
                If(   (l3r(iD)%time(I) > 0 .AND. l3r_temp(iD)%time(J) > 0 .AND. l3r_temp(iD)%time(J+1) > 0 ) .AND. &
                    & l3r(iD)%time(I) > l3r_temp(iD)%time(J) .AND. & 
                    & l3r(iD)%time(I) < l3r_temp(iD)%time(J+1) ) Then
                   
                   l3r(iD)%l2gpValue(1, iP, I) =  & 
                        & l3r_temp(iD)%l2gpValue(1, iN, J) + &
                        & (l3r(iD)%time(I) - l3r_temp(iD)%time(J))*&
                        & (l3r_temp(iD)%l2gpValue(1, iN, J+1) - & 
                        & l3r_temp(iD)%l2gpValue(1, iN, J))/	&
                        & (l3r_temp(iD)%time(J+1) - l3r_temp(iD)%time(J))
                   If( J > icount) Then
                      icount = J
                   End If
                   EXIT
                End If
             ENDIF
          ENDDO
       ENDDO
  
     !-----------------------------------
     END SUBROUTINE Residual2L2Grid
     !-----------------------------------

     subroutine my_sortp_r4(sngl_array, n, m, pt_array)
       real(r4), dimension(*), intent(in) :: sngl_array
       integer, intent(in)                :: n
       integer, intent(in)                :: m
       integer, dimension(*), intent(out) :: pt_array
       call ssortp(sngl_array, n, m, pt_array)
     end subroutine my_sortp_r4
     
     subroutine my_sortp_r8(dbl_array, n, m, pt_array)
       real(r8), dimension(*), intent(in) :: dbl_array
       integer, intent(in)                :: n
       integer, intent(in)                :: m
       integer, dimension(*), intent(out) :: pt_array
       call dsortp(dbl_array, n, m, pt_array)
     end subroutine my_sortp_r8
     
     !===================
   END MODULE Synoptic
!===================

! $Log$
! Revision 1.40  2006/05/24 19:14:30  cvuu
! Fix the deallocate problem
!
! Revision 1.39  2006/02/28 17:56:56  cvuu
! V2.00 commit
!
! Revision 1.37  2004/09/30 17:53:10  cvuu
! bug fixes
!
! Revision 1.36  2004/06/22 16:55:47  cvuu
! Using SetupNewL2GPRecord for L3 Residual based on the nLevels from L3 not L2
!
! Revision 1.35  2004/06/02 20:07:28  ybj
! *** empty log message ***
!
! Revision 1.34  2004/05/13 23:50:52  ybj
! *** empty log message ***
!
! Revision 1.33  2004/05/06 13:16:46  cvuu
! remove print statements
!
! Revision 1.32  2004/05/04 15:33:15  cvuu
! v1.4.3: Use int array for Date in Data Field
!
! Revision 1.31  2004/01/07 21:43:18  cvuu
! version 1.4 commit
!
! Revision 1.30  2003/04/30 18:15:48  pwagner
! Work-around for LF95 infinite compile-time bug
!
! Revision 1.29  2003/03/22 02:59:22  jdone
! individual allocate/deallocate, use only and indentation added
!
! Revision 1.28  2003/02/06 19:51:31  pwagner
! Fixed problem arising from different num types in calls to _sortp
!
! Revision 1.27  2002/05/08 16:23:47  ybj
! Fix residual problem in order to interpolate to L2GP grid
!
! Revision 1.26  2002/04/30 20:05:21  jdone
! indexing issues related to nc* arrays
!
! Revision 1.25  2002/04/11 00:53:57  jdone
! initialization added
!
! Revision 1.24  2002/04/10 21:30:33  jdone
! associated statements before deallocate added
!
! Revision 1.23  2002/04/01 21:57:13  jdone
! check division by zero
!
! Revision 1.22  2002/03/27 22:45:17  jdone
! deallocate statements added
!
! Revision 1.21  2002/03/27 21:35:14  jdone
! allocate statements checked and maxDiff is initialized
!
! Revision 1.20  2002/02/20 22:30:33  ybj
! *** empty log message ***
!
! Revision 1.19  2002/02/20 22:13:54  ybj
! *** empty log message ***
!
! Revision 1.18  2001/09/27 20:38:31  ybj
! Add Precision Calculation & ado mode
!
! Revision 1.17  2001/09/27 20:07:14  ybj
! Add Precision Calculation
!
! Revision 1.16  2001/09/06 18:45:37  nakamura
! Initialized pStart/EndIndex; removed some unused variables.
!
! Revision 1.15  2001/08/13 19:28:03  ybj
! *** empty log message ***
!
! Revision 1.14  2001/08/13 16:42:20  ybj
! *** empty log message ***
!
! Revision 1.13  2001/04/12 16:04:41  ybj
! reasonable values
!
! Revision 1.12  2001/04/11 18:29:22  ybj
! reasonable values
!
! Revision 1.11  2001/03/08 00:53:15  ybj
! *** empty log message ***
!
! Revision 1.10  2001/03/08 00:37:17  ybj
! *** empty log message ***
!
! Revision 1.9  2001/03/07 23:12:15  ybj
! *** empty log message ***
!
! Revision 1.8  2001/03/05 19:58:59  ybj
! with selected pressure levels
!
! Revision 1.7  2001/03/03 01:38:02  ybj
! with selected pressure levels
!
! Revision 1.6  2001/03/03 00:45:44  ybj
! with selected pressure levels
!
! Revision 1.5  2001/03/03 00:20:07  ybj
! with selected pressure levels
!
! Revision 1.3  2001/02/28 19:27:03  ybj
! Fix residue allocation & flags initialization
!
! Revision 1.2  2001/02/28 17:15:06  ybj
! Fix residue structure alocation
!
! Revision 1.1  2001/02/27 20:51:28  ybj
! Daily Map Core Processing
!
! Revision 1.1  2000/10/05 18:38:48  nakamura
! Program split from synoptic.f90 and modified to be more like the standard template.
!


