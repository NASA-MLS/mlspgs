! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L1BOutUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OutputL1Bdata, WriteHdrAnnots

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OutputL1Bdata
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo
    USE Orbit, ONLY: altG, altT, ascTAI, dscTAI, numOrb, orbIncline, &
         orbitNumber, scanRate, scanRateT
    USE TkL1B, ONLY: L1BOA_MAF
    USE OutputL1B, ONLY: OutputL1B_rad
    USE MLSL1Rad, ONLY: L1Brad
    USE MLSL1Config, ONLY: L1Config
    USE Calibration, ONLY: CalWin

    INTEGER, SAVE :: MAFno = 0, counterMAF

    counterMAF = CalWin%MAFdata(CalWin%central)%EMAF%TotalMAF
    MAFno = MAFno + 1

    IF (L1Config%Globals%ProduceL1BOA) THEN

       CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
            l1bFileInfo%OAId, MAFinfo, MAFno, numOrb, orbIncline, &
            orbitNumber, scanRate, scanRateT, L1BFileInfo)

    ENDIF

    CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, MAFinfo%startTAI, &
         L1Brad)

    PRINT *, "outputting l1b for MAF no ", MAFno

  END SUBROUTINE OutputL1Bdata

!=============================================================================
  SUBROUTINE WriteHdrAnnots (FileName, File_id, HDFversion)
!=============================================================================

    USE OpenInit, ONLY: antextPCF, antextCF
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

END MODULE L1BOutUtils

! $Log$
! Revision 2.6  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.4  2002/07/17 14:27:26  perun
! Added ProduceL1BOA flag
!
! Revision 2.3  2002/07/17 14:19:54  perun
! Uncommented L1BOA_MAF
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:49:58  perun
! Version 0.5 commit
!
