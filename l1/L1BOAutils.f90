! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
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
    USE L0_sci_tbls, ONLY: APE2_dflt, TSE2_dflt
    USE MLSL1Config, ONLY: L1Config

    LOGICAL, INTENT(OUT) :: more_data

    INTEGER, SAVE :: MAFno = 0, counterMAF = 0
    INTEGER :: i
    REAL :: scAngleG(0:(MaxMIFs-1)), scAngleT(0:(MaxMIFs-1))
    REAL, PARAMETER :: APE_eps = 27.7340
    REAL, PARAMETER :: TSE_eps = 26.301

    MAFno = MAFno + 1
    counterMAF = counterMAF + 1   ! fake the counter MAF

    DO i = 0, (MaxMIFs - 1)
       scAngleG(i) = APE2_dflt(i)
       scAngleG(i) = MOD ((scAngleG(i) + APE_eps), 360.0)
       scAngleT(i) = TSE2_dflt(i)
       scAngleT(i) = MOD ((scAngleT(i) + TSE_eps), 360.0)
       IF (i > 0) THEN  ! use previous until further notice
          IF (scAngleG(i) < 0.0) scAngleG(i) = scAngleG(i-1)
          IF (scAngleT(i) < 0.0) scAngleT(i) = scAngleT(i-1)
       ENDIF
    ENDDO

    CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
         l1bFileInfo%OAId, MAFinfo, MAFno, numOrb, scAngleG, scAngleT)

    PRINT *, "outputting l1boa for MAF no ", MAFno

    MAFinfo%startTAI = MAFinfo%startTAI + MAFinfo%MIFsPerMAF * MAFinfo%MIF_dur
    more_data = MAFinfo%startTAI <= L1Config%Expanded_TAI%endTime

  END SUBROUTINE OutputL1BOA

!=============================================================================
  SUBROUTINE WriteHdrAnnots (FileName, File_id, HDFversion)
!=============================================================================

    USE OpenInitBOA, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr
    USE MLSFiles, ONLY: HDFVERSION_4
    USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
    USE MLSAuxData, ONLY: Build_MLSAuxData, MaxCharFieldLen

    CHARACTER(LEN=*), INTENT(IN) :: Filename
    INTEGER, INTENT(IN) :: File_id, HDFversion
    CHARACTER(len=MaxCharFieldLen) :: cbuf

    TYPE (DataProducts_T) :: dataset

    IF (HDFversion == HDFVERSION_4) THEN
       CALL WritePCF2Hdr (FileName, anTextPCF)
       CALL WritePCF2Hdr (FileName, anTextCF)
    ELSE
       CALL Deallocate_DataProducts (dataset)
       dataset%name      = 'TextPCF'
       dataset%data_type = 'character'
       cbuf = TRANSFER (anTextPCF, cbuf)
       CALL Build_MLSAuxData (File_id, dataset, cbuf, &
            char_length=SIZE(anTextPCF))
       dataset%name      = 'TextCF'
       cbuf = TRANSFER (anTextCF, cbuf)
       CALL Build_MLSAuxData (File_id, dataset, cbuf, &
            char_length=SIZE(anTextCF))
    ENDIF

  END SUBROUTINE WriteHdrAnnots

END MODULE L1BOAutils

! $Log$
! Revision 2.1  2003/10/24 19:38:36  perun
! Version 1.3 commit
!
!
