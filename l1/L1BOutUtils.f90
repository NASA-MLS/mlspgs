! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE L1BOutUtils
!=============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OutputL1Bdata, OutputL1BOA, WriteHdrAnnots

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

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
    INTEGER :: i, MIFsPerMAF
    INTEGER, PARAMETER :: last_MIF_indx = MaxMIFs - 1
    REAL :: scAngleG(0:last_MIF_indx), scAngleT(0:last_MIF_indx)
    REAL :: encAngleG(0:last_MIF_indx), encAngleT(0:last_MIF_indx)
    REAL :: APE_pos_P(2,0:last_MIF_indx), TSSM_pos_P(2,0:last_MIF_indx)

    IF (.NOT. L1Config%Globals%ProduceL1BOA) RETURN

    counterMAF = CurMAFdata%EMAF%TotalMAF
    MIFsPerMAF = CurMAFdata%EMAF%MIFsPerMAF

    MAFno = MAFno + 1

    DO i = 0, last_MIF_indx
       scAngleG(i) = CurMAFdata%SciPkt(i)%scAngleG
       scAngleT(i) = CurMAFdata%SciPkt(i)%scAngleT
       encAngleG(i) = CurMAFdata%SciPkt(i)%APE_theta
       encAngleT(i) = CurMAFdata%SciPkt(i)%TSSM_theta
       IF (i > 0) THEN  ! use previous until further notice
          IF (scAngleG(i) < 0.0) scAngleG(i) = scAngleG(i-1)
          IF (scAngleT(i) < 0.0) scAngleT(i) = scAngleT(i-1)
       ENDIF
       APE_pos_P(:,i) = CurMAFdata%SciPkt(i)%APE_pos_P
       TSSM_pos_P(:,i) = CurMAFdata%SciPkt(i)%TSSM_pos_P
    ENDDO

    CALL L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, &
         l1bFileInfo%OAId, MAFinfo, MAFno, MIFsPerMAF, numOrb, scAngleG, &
         scAngleT, encAngleG, encAngleT, APE_pos_P, TSSM_pos_P)

  END SUBROUTINE OutputL1BOA

!=============================================================================
  SUBROUTINE OutputL1Bdata
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, OA_counterMAF, &
         OA_counterIndex, MaxMIFs, DACSnum
    USE OutputL1B, ONLY: OutputL1B_rad, OutputL1B_diags
    USE MLSL1Rad, ONLY: L1Brad
    USE EngTbls, ONLY: Reflec
    USE Calibration, ONLY: CalWin, MAFdata_T

    INTEGER, SAVE :: MAFno = 0, counterMAF, MAFindex = 1
    INTEGER :: bank, MIF
    TYPE (MAFdata_T), POINTER :: CurMAFdata
    INTEGER, PARAMETER :: last_MIF = MaxMIFs - 3
    REAL :: TP(0:last_MIF,DACSnum), TPdigP(0:last_MIF,DACSnum), &
         TPdigN(0:last_MIF,DACSnum)

    CurMAFdata => CalWin%MAFdata(CalWin%central)
    counterMAF = CurMAFdata%EMAF%TotalMAF

    IF (OA_counterIndex == 0) THEN   ! No OA data available
       MAFno = MAFno + 1
    ELSE
       DO
          IF (counterMAF == OA_counterMAF(MAFindex)) EXIT
          MAFindex = MAFindex + 1
       ENDDO
       MAFno = MAFindex
    ENDIF

    PRINT *, "outputting l1b rads for MAF no ", MAFno
    CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, Reflec, &
         MAFinfo%startTAI, L1Brad)

! Save DACS TP values:

    DO bank = 1, DACSnum
       DO MIF = 0, last_MIF
          TP(MIF,bank) = CurMAFdata%SciPkt(MIF)%TP(bank)
          TPdigP(MIF,bank) = CurMAFdata%SciPkt(MIF)%TPdigP(bank)
          TPdigN(MIF,bank) = CurMAFdata%SciPkt(MIF)%TPdigN(bank)
       ENDDO
    ENDDO

    PRINT *, "outputting l1b diags for MAF no ", MAFno
    CALL OutputL1B_diags (L1BFileInfo%DiagId, MAFno, counterMAF, &
         MAFinfo%startTAI, TP=TP, TPdigP=TPdigP, TPdigN=TPdigN)



  END SUBROUTINE OutputL1Bdata

!=============================================================================
  SUBROUTINE WriteHdrAnnots (FileName, HDFversion)
!=============================================================================

    USE INTRINSIC, ONLY: l_hdf
    USE OpenInit, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr

    CHARACTER(LEN=*), INTENT(IN) :: Filename
    INTEGER, INTENT(IN) :: HDFversion

    CALL WritePCF2Hdr (FileName, anTextPCF, hdfVersion=HDFversion, &
         & fileType=l_hdf)
    CALL WritePCF2Hdr (FileName, anTextCF, hdfVersion=HDFversion, &
         & fileType=l_hdf, name='/LCF')

  END SUBROUTINE WriteHdrAnnots

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE L1BOutUtils

! $Log$
! Revision 2.16  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.15.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.15  2006/06/14 13:46:34  perun
! Assemble TP data for output to the DIAG file
!
! Revision 2.14  2005/08/24 15:50:57  perun
! Assemble and pass APE and TSSM pos_P for output
!
! Revision 2.13  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.12  2004/11/10 15:32:44  perun
! Output encoder values; deal with gaps compared to L1BOA records
!
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
