! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L1BOAutils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OutputL1BOA, WriteHdrAnnots

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OutputL1BOA (more_data)
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, MaxMIFs
    USE Orbit, ONLY: altG, altT, ascTAI, dscTAI, numOrb, orbIncline, &
         orbitNumber, scanRate, scanRateT
    USE TkL1B, ONLY: L1BOA_MAF
    USE L0_sci_tbls, ONLY: APE_theta_dflt, TSSM_theta_dflt
    USE MLSL1Config, ONLY: L1Config

    LOGICAL, INTENT(OUT) :: more_data

    INTEGER, SAVE :: MAFno = 0, counterMAF = 0
    INTEGER :: i, MIFsPerMAF
    INTEGER, PARAMETER :: last_MIF_indx = MaxMIFs - 1
    REAL :: scAngleG(0:last_MIF_indx), scAngleT(0:last_MIF_indx)
    REAL :: encAngleG(0:last_MIF_indx), encAngleT(0:last_MIF_indx)
    REAL, PARAMETER :: APE_eps = 27.7340
    REAL, PARAMETER :: TSE_eps = 26.301

    MAFno = MAFno + 1
    counterMAF = counterMAF + 1   ! fake the counter MAF

    DO i = 0, (MaxMIFs - 1)
       scAngleG(i) = APE_theta_dflt(i)
       scAngleG(i) = MOD ((scAngleG(i) + APE_eps), 360.0)
       scAngleT(i) = TSSM_theta_dflt(i)
       scAngleT(i) = MOD ((scAngleT(i) + TSE_eps), 360.0)
       IF (i > 0) THEN  ! use previous until further notice
          IF (scAngleG(i) < 0.0) scAngleG(i) = scAngleG(i-1)
          IF (scAngleT(i) < 0.0) scAngleT(i) = scAngleT(i-1)
       ENDIF
    ENDDO

    CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
         l1bFileInfo%OAId, MAFinfo, MAFno, MIFsPerMAF, numOrb, scAngleG, &
         scAngleT, encAngleG, encAngleT)

    PRINT *, "outputting l1boa for MAF no ", MAFno

    MAFinfo%startTAI = MAFinfo%startTAI + MAFinfo%MIFsPerMAF * MAFinfo%MIF_dur
    more_data = MAFinfo%startTAI <= L1Config%Expanded_TAI%endTime

  END SUBROUTINE OutputL1BOA

!=============================================================================
  SUBROUTINE WriteHdrAnnots (FileName, HDFversion)
!=============================================================================

    USE Intrinsic, ONLY: l_hdf
    USE OpenInitBOA, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr

    CHARACTER(LEN=*), INTENT(IN) :: Filename
    INTEGER, INTENT(IN) :: HDFversion

    CALL WritePCF2Hdr (FileName, anTextPCF, hdfVersion=HDFversion, &
         & fileType=l_hdf)
    CALL WritePCF2Hdr (FileName, anTextCF, hdfVersion=HDFversion, &
         & fileType=l_hdf, name='/LCF')

  END SUBROUTINE WriteHdrAnnots

END MODULE L1BOAutils

! $Log$
! Revision 2.4  2004/11/10 17:22:00  perun
! Change call to L1BOA_MAF with new arguments
!
! Revision 2.3  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.2  2004/01/09 20:02:57  perun
! Update BOA to HDF 5
!
! Revision 2.1  2003/10/24 19:38:36  perun
! Version 1.3 commit
!
!
