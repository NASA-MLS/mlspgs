! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1       ! MLS Level 1 software
!============================================================================

  USE SDPToolkit
  USE MLSL1Common
  USE MLSL1Rad
  USE MLSMessageModule
  USE OpenInit
  USE OutputL1B
  USE Orbit
  USE TkL1B
  USE SortQualify
  USE Calib
  USE Close_files

  IMPLICIT NONE

  !--------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=130) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !--------------------------------------------------------------------------

  INTEGER :: returnstatus, counterMAF, i, numOrb
  INTEGER :: orbitNumber(max_orbits)
  LOGICAL :: more_data

  CHARACTER (LEN=27) :: UTC_start, UTC_end

  REAL :: scanRate(lenG)
  REAL :: scanRateT(lenT)
  REAL(r8) :: altG, altT, orbIncline
  REAL(r8) :: ascTAI(max_orbits), dscTAI(max_orbits)

  TYPE (TAI93_Range_T) :: procRange
  TYPE (L0Info_T) :: l0Info
  TYPE (L1BInfo_T) :: l1bInfo
  TYPE (L0Sci_T) :: L0Sci(MaxMIFs) ! will allocate this in future versions
  TYPE (Radiance_T), DIMENSION (:), POINTER :: L1Brad
  TYPE (MAFinfo_T) :: MAFinfo

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 data processing started")

  CALL OpenAndInitialize (procRange, l0Info, l1bInfo, L1Brad)

  returnstatus = PGS_TD_TAItoUTC (procRange%startTime, UTC_start)
  returnstatus = PGS_TD_TAItoUTC (procRange%endTime, UTC_end)

  CALL Orbit_init(procRange, UTC_start, altG, altT, ascTAI, dscTAI, numOrb, &
                  orbIncline, orbitNumber, scanRate, scanRateT)

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Processing start time: "//UTC_start)
  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Processing end time: "//UTC_end)

  more_data = .true.
  i = 0

  DO WHILE (more_data)

     counterMAF = i
     i = i+1
  
     CALL SortAndQualify (l0Info, l1bInfo, l0Sci, MAFinfo, more_data)

     IF (more_data) THEN
        CALL L1boa_MAF(altG, altT, ascTAI, counterMAF, dscTAI, &
                       l1bInfo%L1BOAId, MAFinfo, i, numOrb, orbIncline, &
                       orbitNumber, scanRate, scanRateT)
        CALL Calibrate (l1bInfo, l0Sci, L1Brad, more_data)
        CALL OutputL1B_rad(i,l1bInfo, counterMAF, L1Brad)
     ENDIF

  END DO

  CALL CloseFiles (l0Info, l1bInfo)

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 data processing successfully completed!")

!=============================================================================
END PROGRAM MLSL1
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.3  2000/02/15 19:01:59  nakamura
! Modified code for parametrization of lenG & lenT, moving _init subroutine to Orbit module.
!

