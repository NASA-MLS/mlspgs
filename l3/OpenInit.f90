! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE OpenInit 
!==============================================================================

  USE L2INTERFACE, ONLY: ReadL2GPAttribute  
  USE L3CF, ONLY: L3CFDef_T, L3CFProd_T, FillL3CF
  USE MLSCF, ONLY: Mlscf_T 
  USE MLSCommon, ONLY: r8
  USE MLSL3Common, ONLY: maxWindow, CCSDS_LEN, FILENAMELEN
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
       & MLSMSG_ALLOCATE
  USE MLSPCF3, ONLY: mlspcf_l3_param_OutputVersion, mlspcf_l3_param_Cycle, &
       & mlspcf_l3_param_L2DayRange, mlspcf_l3_param_RangDays, & 
       & mlspcf_pcf_start,mlspcf_l3cf_start
  USE MLSStrings, only: utc_to_yyyymmdd
  USE Output_m, only: output
  USE PCFHdr, ONLY: CreatePCFAnnotation, GlobalAttributes, FillTAI93Attribute
  USE SDPToolkit, ONLY: PGS_S_SUCCESS, pgs_pc_getConfigData, & 
       & Pgs_td_utcToTAI, max_orbits, spacecraftId, PGS_IO_Gen_CloseF, &
       & Pgs_io_gen_openF, PGSd_IO_Gen_RSeqFrm, Pgs_pc_getReference, &
       & max_orbits
  USE GETCF_M, only: GetCF, InitGetCF
  IMPLICIT NONE
  PUBLIC

  PRIVATE :: ID, ModuleName

  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &                                                
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------
  
  ! Contents:
  
  ! Definitions -- PCFData_T
  ! Subroutines -- GetPCFParameters
  !                SetProcessingWindow
  !                OpenAndInitialize
  !                AvgOrbPeriod
  
  ! Remarks:  This is a prototype module for the routines needed for the 
  ! L3 Daily
  ! Open/Init task.
  
  ! Parameters
  
  ! This data type is used to store User-defined Runtime Parameters and other
  ! information taken from the PCF.

  TYPE PCFData_T

     ! name of the log file (without the path)
     
     CHARACTER (LEN=FileNameLen) :: logGranID
     
     ! version string in PCF output file names
     
     CHARACTER (LEN=15) :: outputVersion	! output files
     
     ! first & last days of input & output windows
     
     CHARACTER (LEN=8) :: l2EndDay		! last day of l2 to read
     CHARACTER (LEN=8) :: l2StartDay		! first day of l2 to read
     CHARACTER (LEN=8) :: l3EndDay		! last day of l3 to write
     CHARACTER (LEN=8) :: l3StartDay		! first day of l3 to write
     
     ! cycle # of processing run
     
     CHARACTER (LEN=4) :: cycle
     
  END TYPE PCFData_T
  
CONTAINS

  !------------------------------------
  SUBROUTINE GetPCFParameters(l3pcf)
  !------------------------------------

    ! Brief description of subroutine
    ! This subroutine retrieves the User-Defined Runtime Parameters from the 
    ! PCF and stores them in variables used in the MLSL3 code.

    ! Arguments

    TYPE( PCFData_T ), INTENT(OUT) :: l3pcf
    
    ! Parameters

    CHARACTER (LEN=*), PARAMETER :: & 
         UDRP_ERR = 'Failed to get PCF runtime parameter :  '
    ! Variables

    CHARACTER (LEN=480) :: msr
    CHARACTER (LEN=FileNameLen) :: name
    CHARACTER (LEN=17) :: range 

    INTEGER ::  indx, mlspcf_log, returnStatus, version

    ! Functions

    INTEGER, EXTERNAL :: pgs_pc_getConfigData

    ! Retrieve values set in PCF & assign them to variables.

    returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_OutputVersion, &
         & l3pcf%outputVersion)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       msr = UDRP_ERR // 'OutputVersion'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_Cycle, l3pcf%cycle)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       msr = UDRP_ERR // 'Cycle'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_L2DayRange, range)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       msr = UDRP_ERR // 'L2DayRange'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    l3pcf%l2StartDay = range(1:8)
    l3pcf%l2EndDay = range(10:17)
    
    returnStatus = pgs_pc_getConfigData(mlspcf_l3_param_RangDays, range)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       msr = UDRP_ERR // 'RangDays'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    
    l3pcf%l3StartDay = range(1:8)
    l3pcf%l3EndDay = range(10:17)

    ! Get the name of the log file from the PCF

    version = 1
    mlspcf_log = 10101

    returnStatus = Pgs_pc_getReference(mlspcf_log, version, name)
    IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
         & ModuleName, 'Error retrieving log file name from PCF.')
    indx = INDEX(name, '/', .TRUE.)
    l3pcf%logGranID = name(indx+1:)

    ! Store appropriate user input as global attributes
    GlobalAttributes%InputVersion = l3pcf%outputVersion
    GlobalAttributes%StartUTC = l3pcf%l2StartDay // &
      & 'T00:00:00.000000Z'
    GlobalAttributes%EndUTC = l3pcf%l2EndDay // &
      & 'T23:59:59.999999Z'
    GlobalAttributes%PGEVersion = 'v1.2'   ! l3pcf%PGEVersion
    call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
      & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
      & GlobalAttributes%GranuleDay) 
    call FillTAI93Attribute
    
  !---------------------------------
  END SUBROUTINE GetPCFParameters
  !---------------------------------

  !-------------------------------------------------------------
  SUBROUTINE SetProcessingWindow (start, end, dates, numDays)
  !-------------------------------------------------------------
  USE dates_module, only: days_in_year

    ! Brief description of subroutine
    ! This subroutine finds each DOY for the L3 processing window 
    ! and calculates the number of days in the window.
    
    ! Arguments
    
    CHARACTER (LEN=8), INTENT(IN) :: end, start
    
    CHARACTER (LEN=8), INTENT(OUT) :: dates(:)

    INTEGER, INTENT(OUT) :: numDays

    ! Parameters

    ! Functions

    ! Variables

    CHARACTER (LEN=4) :: cYr
    CHARACTER (LEN=3) :: cDOY

    INTEGER :: i, iDOY, iEnd, iYr, lastDOY

    ! Initializations outside the DO loop

    READ(start(1:4), '(I4)') iYr
    lastDOY = days_in_year(iYr)
    READ( start(6:8), '(I3)' ) iDOY
    READ( end(6:8), '(I3)' ) iEnd

    DO i = 1, maxWindow

       ! Increment DOY

       iDOY = iDOY + 1

       ! Check that the incrementation hasn't crossed a year boundary

       IF (iDOY > lastDOY) THEN
          iYr = iYr + 1
          iDOY = 1
       ENDIF

       ! Check that the incrementation hasn't crossed the window boundary

       IF (iDOY > iEnd) EXIT

       ! Convert Yr back to CHARACTER

       WRITE( cYr, '(I4)' ) iYr

       ! Convert DOY back to CHARACTER, adding leading zeroes if necessary

       WRITE( cDOY, '(I3)' ) iDOY

       IF (iDOY < 10) THEN
          cDOY = '00' // cDOY(3:3)
       ELSE IF (iDOY < 100) THEN
          cDOY = '0' // cDOY(2:3)
       ENDIF

       ! Concatenate them into a date

       dates(i+1) = cYr // '-' // cDOY

    ENDDO

    numDays = i
    dates(1) = start
    dates(numDays) = end

  !------------------------------------
  END SUBROUTINE SetProcessingWindow
  !------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE OpenAndInitialize (l3pcf, cf, l3cf, cfDef, anText, avgPer)
  !-----------------------------------------------------------------------

    ! Brief description of subroutine
    ! This subroutine performs the Open/Init task in the MLSL3 program.

    ! Arguments

    TYPE( L3CFDef_T ), INTENT(OUT) :: cfDef

    TYPE( Mlscf_T ), INTENT(OUT) :: cf

    TYPE( PCFData_T ), INTENT(OUT) :: l3pcf

    TYPE( L3CFProd_T ), POINTER :: l3cf(:) 

    CHARACTER (LEN=1), POINTER :: anText(:) 

    REAL(r8), POINTER :: avgPer(:) 

    ! Parameters

    ! Functions

    ! Variables

    CHARACTER (LEN=480) :: msr
    CHARACTER (LEN=8) :: dates(maxWindow)

    REAL(r8) :: period, sum

    INTEGER :: error, numDays, pcfId, returnStatus, err, i, j, count

         numDays = 0
             err = 0
           error = 0
           pcfId = 0
    returnStatus = 0
             sum = 0
           count = 0

    ! Read the PCF into an annotation for file headers

    CALL CreatePCFAnnotation(mlspcf_pcf_start, anText)

    ! Retrieve values set in PCF & assign them to variables.

    CALL GetPCFParameters(l3pcf)

    ! Calculate the dates, number of days for the processing window
    
    CALL SetProcessingWindow(l3pcf%l3StartDay, l3pcf%l3EndDay, dates, &
         & numDays)

    ! Open the configuration file

    returnStatus = pgs_io_gen_openF (mlspcf_l3cf_start, &
         & PGSd_IO_Gen_RSeqFrm, 0, pcfId, 1)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       msr = MLSMSG_Fileopen // 'l3cf.'
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    
    ! Read the configuration file into an MLSCF_T structure
    
    CALL InitGetCF
    
    CALL getCF(cf, error, inUnit=pcfId)
    
    ! Close the configuration file
    
    returnStatus = pgs_io_gen_closeF (pcfId)
    IF (returnStatus /= PGS_S_SUCCESS) & 
         & CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to close l3cf.')

    ! Check the parser output; fill L3CFProd_T
    
    CALL FillL3CF(cf, l3pcf%outputVersion, dates, numDays, l3cf, cfDef)

    ! Get L2GP Attributes

    CALL ReadL2GPAttribute (l3pcf%l3StartDay, l3pcf%l3EndDay, &
	& l3cf(1)%fileTemplate)

! If the average orbit period is either not in the cf file or is a number <=
!     Then calculate the average orbit period for each day in the input window.
! If there are orbit period in the l2gp files, use them instead of calculate or ! use from cf file 9/16/03

    IF (numDays .gt. 0) THEN
       allocate(avgPer(numDays), stat=err)
       IF ( err /= 0 ) THEN
          msr = MLSMSG_Allocate // ' array of period averages.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       DO i = 1, numDays
          IF (GlobalAttributes%OrbPeriodDays(1,i) /= -1.0) THEN 
             !call output('first orbitPeriod = ', advance='no')
             DO j = 1, max_orbits
               IF (GlobalAttributes%OrbPeriodDays(j,i) /= -1.0) THEN 
		   sum = sum + GlobalAttributes%OrbPeriodDays(j,i)
                   count = count + 1
	       END IF	 
             END DO
             avgPer(i) = sum/count 
             !call output('avgPer(i) = ', advance='no')
             !call output(avgPer(i), advance='yes')
          ELSE  
             avgPer(i) = cfDef%averageOrbitalPeriod
          END IF
       END DO

    ELSE
       call AvgOrbPeriod(l3pcf%l2StartDay, l3pcf%l2EndDay, avgPer)
    END IF

  !----------------------------------
  END SUBROUTINE OpenAndInitialize
  !----------------------------------

  !----------------------------------------------------
  SUBROUTINE AvgOrbPeriod (startDay, endDay, avgPer)
  !----------------------------------------------------

    ! Brief description of subroutine
    ! This program calculates the average orbital period 
    ! for each day in the input window.
            
    ! Arguments

    CHARACTER (LEN=8), INTENT(IN) :: endDay, startDay

    REAL(r8), POINTER :: avgPer(:)

    ! Parameters

    INTEGER, PARAMETER :: numOffset = 43200

    ! Variables

    CHARACTER (LEN=CCSDS_LEN) :: startA, orbitAscendTime(max_orbits), &
         & orbitDescendTime(max_orbits)
    CHARACTER (LEN=480) :: msg, msr
    CHARACTER (LEN=32) :: mnemonic
    CHARACTER (LEN=26) :: startB
    CHARACTER (LEN=8) :: dates(maxWindow)

    REAL(r8) :: sum, offsets(numOffset), ascTAI(max_orbits), & 
         & int(max_orbits), orbitDownLongitude(max_orbits)

    INTEGER :: err, i, j, numDays, numOrbits, returnStatus, & 
         & orbitNumber(max_orbits)

    ! Functions
            
    INTEGER, EXTERNAL :: Pgs_eph_getEphMet, Pgs_td_asciiTime_bToA, & 
         & Pgs_td_timeInterval, Pgs_td_utcToTAI

    ! Calculate the offsets

    DO i = 1, numOffset
       offsets(i) = (i-1) * 2.0
    ENDDO

    ! Find the number of days in the input window
    
    CALL SetProcessingWindow(startDay, endDay, dates, numDays)

    if (numDays .gt. 0) then 

       ALLOCATE(avgPer(numDays), STAT=err)
       IF ( err /= 0 ) THEN
          msr = MLSMSG_Allocate // ' array of period averages.'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       
    endif

    ! For each day in the input window,

    DO i = 1, numDays

       ! Convert the start time to CCSDS A format

       startB = dates(i) // 'T00:00:00.000000Z'

       returnStatus = Pgs_td_asciiTime_bToA(startB, startA)
       IF (returnStatus /= PGS_S_SUCCESS) & 
            & CALL MLSMessage(MLSMSG_Error, &
            & ModuleName, 'Failed to convert input date to CCSDS A')
       
       returnStatus = PGS_EPH_GetEphMet(spacecraftId, &
            & numOffset, startA, &
            & offsets, numOrbits,&
            & orbitNumber, &
            & orbitAscendTime, &
            & orbitDescendTime,&
            & orbitDownLongitude)
       
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
          msr = mnemonic // ':  ' // msg
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       ! Convert orbitAscendTime to TAI93
       
       DO j = 1, numOrbits
          returnStatus = Pgs_td_utcToTAI(orbitAscendTime(j),ascTAI(j))
          IF (returnStatus /= PGS_S_SUCCESS) & 
               & CALL MLSMessage(MLSMSG_Error, &
               & ModuleName, 'Error converting from CCSDSB to TAI.')
       ENDDO

       ! Get the differences between successive times 
       ! and add them for the day

       sum = 0.0 

       DO j = 1, numOrbits-1
          returnStatus = & 
               & Pgs_td_timeInterval( ascTAI(j), ascTAI(j+1), int(j) )
          IF (returnStatus /= PGS_S_SUCCESS) & 
               & CALL MLSMessage(MLSMSG_Error, ModuleName, & 
               & 'Failed to find interval between orbitsAscendTimes.')
          sum = sum + int(j)
       ENDDO

       ! Get the average period for the day

       if (numOrbits .gt. 1) then 
          avgPer(i) = sum/float(numOrbits-1)
       else
          avgPer(i) = 0.0
       endif

       ! Find Standard Deviation of the period. 

       sum = 0.0 

       DO j = 1, numOrbits-1
          sum = sum + (int(j) - avgPer(i))*(int(j) - avgPer(i)) 
       ENDDO

    ENDDO

  !-----------------------------
  END SUBROUTINE AvgOrbPeriod
  !-----------------------------

!==================
END MODULE OpenInit
!==================

! $Log$
! Revision 1.12  2003/06/03 20:45:04  pwagner
! Fills global attributes
!
! Revision 1.11  2003/04/30 18:15:48  pwagner
! Work-around for LF95 infinite compile-time bug
!
! Revision 1.10  2003/03/22 02:51:22  jdone
! average orbit retrieved from cf or calculated
!
! Revision 1.9  2002/04/10 22:09:02  jdone
! Division by zero checked.
!
! Revision 1.8  2001/03/27 19:35:37  nakamura
! Changed from PCFModule to PCFHdr; save logGranID with .dat; allowed for a 1 day output window.
!
! Revision 1.7  2001/02/21 21:11:59  nakamura
! Changed MLSPCF to MLSPCF3; removed inputVersion.
!
! Revision 1.6  2001/01/18 16:52:28  nakamura
! Added type L3CFDef_T; moved minDays from PCF to cf.
!
! Revision 1.5  2001/01/16 17:48:09  nakamura
! Added code for annotating and for retrieving a name template for the log in the cf.
!
! Revision 1.4  2000/12/29 21:42:27  nakamura
! Added subroutine AvgOrbPeriod.
!
! Revision 1.3  2000/11/15 21:28:05  nakamura
! Changed MLSPCF parameter names to mlspcf_l3_param_*.
!
! Revision 1.2  2000/10/24 19:31:48  nakamura
! Added logGranID to PCFData_T; replaced ReadParseMLSCF stub with calls to actual parser.
!
! Revision 1.1  2000/10/17 20:26:16  nakamura
! Module for the Open/Init task.
!
