! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L1BOutUtils
!=============================================================================

  USE MLSL1Common, ONLY: L1BFileInfo
  USE Orbit, ONLY: altG, altT, ascTAI, dscTAI, numOrb, orbIncline, &
         orbitNumber, scanRate, scanRateT
  USE TkL1B, ONLY: L1BOA_MAF
  USE SortQualify, ONLY: MAFinfo
  USE OutputL1B, ONLY: OutputL1B_rad
  USE MLSL1Rad, ONLY: L1Brad

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE OutputL1Bdata

    INTEGER, SAVE :: MAFno = 0, counterMAF

    counterMAF = MAFno
    MAFno = MAFno + 1

    CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
         l1bFileInfo%OAId, MAFinfo, MAFno, numOrb, orbIncline, &
         orbitNumber, scanRate, scanRateT)

    CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, L1Brad)

    print *, "outputting l1b for MAF no ", MAFno

  END SUBROUTINE OutputL1Bdata

END MODULE L1BOutUtils

! $Log$
! Revision 2.1  2001/02/23 20:49:58  perun
! Version 0.5 commit
!
