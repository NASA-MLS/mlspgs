! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L1BOutUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OutputL1Bdata, OutputL1BOA, WriteHdrAnnots

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OutputL1BOA (CurMAFdata)
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, MaxMIFs
    USE Orbit, ONLY: altG, altT, ascTAI, dscTAI, numOrb
    USE TkL1B, ONLY: L1BOA_MAF
    USE MLSL1Config, ONLY: L1Config
    USE Calibration, ONLY: MAFdata_T

    TYPE (MAFdata_T) :: CurMAFdata
    INTEGER, SAVE :: MAFno = 0, counterMAF
    INTEGER :: i
    INTEGER, PARAMETER :: last_MIF_indx = MaxMIFs - 1
    REAL :: scAngleG(0:last_MIF_indx), scAngleT(0:last_MIF_indx)

    IF (.NOT. L1Config%Globals%ProduceL1BOA) RETURN

    counterMAF = CurMAFdata%EMAF%TotalMAF

    MAFno = MAFno + 1

    DO i = 0, last_MIF_indx
       scAngleG(i) = CurMAFdata%SciPkt(i)%scAngleG
       scAngleT(i) = CurMAFdata%SciPkt(i)%scAngleT
       IF (i > 0) THEN  ! use previous until further notice
          IF (scAngleG(i) < 0.0) scAngleG(i) = scAngleG(i-1)
          IF (scAngleT(i) < 0.0) scAngleT(i) = scAngleT(i-1)
       ENDIF
    ENDDO

    CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
         l1bFileInfo%OAId, MAFinfo, MAFno, numOrb, scAngleG, scAngleT)

  END SUBROUTINE OutputL1BOA

!=============================================================================
  SUBROUTINE OutputL1Bdata
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo
    USE OutputL1B, ONLY: OutputL1B_rad, OutputL1B_diags
    USE MLSL1Rad, ONLY: L1Brad
    USE EngTbls, ONLY: Reflec
    USE Calibration, ONLY: CalWin, MAFdata_T

    INTEGER, SAVE :: MAFno = 0, counterMAF
    TYPE (MAFdata_T), POINTER :: CurMAFdata

    CurMAFdata => CalWin%MAFdata(CalWin%central)
    counterMAF = CurMAFdata%EMAF%TotalMAF

    MAFno = MAFno + 1

!    CALL OutputL1BOA (CurMAFdata)

    CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, Reflec, &
         MAFinfo%startTAI, L1Brad)

    CALL OutputL1B_diags (L1BFileInfo%DiagId, MAFno, counterMAF, &
         MAFinfo%startTAI)

    PRINT *, "outputting l1b for MAF no ", MAFno

  END SUBROUTINE OutputL1Bdata

!=============================================================================
  SUBROUTINE WriteHdrAnnots (FileName, HDFversion)
!=============================================================================

    USE Intrinsic, ONLY: l_hdf
    USE OpenInit, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr

    CHARACTER(LEN=*), INTENT(IN) :: Filename
    INTEGER, INTENT(IN) :: HDFversion

    CALL WritePCF2Hdr (FileName, anTextPCF, hdfVersion=HDFversion, &
         & fileType=l_hdf)
    CALL WritePCF2Hdr (FileName, anTextCF, hdfVersion=HDFversion, &
         & fileType=l_hdf, name='/LCF')

  END SUBROUTINE WriteHdrAnnots

END MODULE L1BOutUtils

! $Log$
! Revision 2.11  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.10  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.9  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.8  2003/07/08 00:17:11  pwagner
! fileType now a lit_name instead of a char string
!
! Revision 2.7  2003/05/30 23:47:51  pwagner
! Prefers to use WritePCF2Hdr for either hdf version
!
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
