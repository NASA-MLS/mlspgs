! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzRadiances ! Determine radiances for the THz module
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: THzNum, THzChans, Deflt_chi2, BandChanBad
  USE THzCalibration, ONLY : CalBuf, Chisq, dLlo, yTsys, nvBounds
  USE MLSL1Rad, ONLY : THzRad

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ProcessLimbData

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE CalcLimbRads (nMAF, ibgn)

    USE THzCalibration, ONLY : Kelvins => Cnts, VarK => VarCnts !Already Kelvins
    USE MLSL1Config, ONLY: MIFsTHz, L1Config

    INTEGER, INTENT (IN) :: nMAF
    INTEGER, INTENT (INOUT) :: ibgn

    INTEGER :: i, iend, calindx, mindx, nBank, nChan, MIF_end, BandNo
    INTEGER :: limb_sw_err(MIFsTHz)
    CHARACTER(LEN=1) :: SwMirPos(MIFsTHz)
    LOGICAL :: do_chi2_err

    calindx = nMAF - CalBuf%Cal_start + 1
    iend = ibgn + CalBuf%MAFdata(calindx)%last_MIF
    MIF_end = ibgn + MIFsTHz - 1

! Mark non-limb data as "bad" with negative precisions

    SwMirPos = CalBuf%MAFdata(calindx)%SciMIF(0:(MIFsTHz-1))%SwMirPos
    WHERE (SwMirPos == "L")
       limb_sw_err = 1
    ELSEWHERE
       limb_sw_err = -1
    ENDWHERE

    DO nBank = 1, THzNum

       THzRad(nBank)%value = 0.0
       THzRad(nBank)%precision  = -1.0
       BandNo = nBank + 14  ! band no for the bank
       do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

       IF (CalBuf%BankGood(nBank)) THEN   ! Have good data

          DO i = ibgn, MIF_end
             mindx = i-ibgn+1        ! index from 1 to MIFsTHz
             DO nChan = 1, THzChans
                THzRad(nBank)%value(nChan,mindx) = Kelvins(nChan,nBank,i)
                IF (THzRad(nBank)%value(nChan,mindx) > -100.0) THEN
                   IF (VarK(nChan,nBank,i) > 0.0) THEN
                      ! All args to intrinsic min must have same
                      ! type and kind type parameter
                      THzRad(nBank)%precision(nChan,mindx) = &
                           VarK(nChan,nBank,i) * MIN ( &
                           REAL(limb_sw_err(mindx)), &
                           BandChanBad%Sign(Bandno, nChan))
                   ELSE
                      THzRad(nBank)%precision(nChan,mindx) = &
                           VarK(nChan,nBank,i)
                   ENDIF
                   IF (do_chi2_err) THzRad(nbank)%precision(nChan,mindx) = &
                        SQRT (deflt_chi2%FB(nchan,BandNo)) * &
                        THzRad(nbank)%precision(nchan,mindx)
                ENDIF
             ENDDO
          ENDDO

       ENDIF
    ENDDO
    ibgn = iend + 1   ! For next call

  END SUBROUTINE CalcLimbRads

  SUBROUTINE ProcessLimbData

    USE MLSL1Common, ONLY: L1BFileInfo, OA_counterMAF, OA_counterIndex
    USE OutputL1B, ONLY: OutputL1B_rad, OutputL1B_DiagsT
    USE EngTbls, ONLY: Reflec_T

    INTEGER :: MAFno, counterMAF, ibgn, nv
    REAL(r8) :: TAI
    TYPE (Reflec_T) :: Reflec
    INTEGER, SAVE :: OrbNo = 1, MAFindex = 1

    print *, 'ProcessLimbData'

    nv = 1
    ibgn = 0
    DO MAFno = CalBuf%Cal_start, CalBuf%Cal_end

       CALL CalcLimbRads (MAFno, ibgn)

       counterMAF = CalBuf%MAFdata(MAFno)%EMAF%TotalMAF
       TAI = CalBuf%MAFdata(MAFno)%SciMIF(0)%secTAI

       IF (OA_counterIndex == 0) THEN
          MAFindex = MAFno
       ELSE
          DO
             IF (counterMAF == OA_counterMAF(MAFindex)) EXIT
             MAFindex = MAFindex + 1
          ENDDO
       ENDIF

print *, "Outputting rad for MAFno: ", MAFindex
       CALL OutputL1B_rad (MAFindex, L1BFileInfo, counterMAF, Reflec, TAI, &
            THzrad)

! Write MAF dimensioned Diags

       CALL OutputL1B_DiagsT (L1BFileInfo%DiagTid, MAFno=MAFindex, &
            counterMAF=counterMAF, MAFStartTimeTAI=TAI, nvBounds=nvBounds(nv))
       nv = nv + 1

    ENDDO

! Write orbit no. dimensioned Diags

    CALL OutputL1B_DiagsT (L1BFileInfo%DiagTid, OrbNo=OrbNo, Chisq=Chisq, &
         dLlo=dLlo, yTsys=yTsys)
    OrbNo = OrbNo + 1

  END SUBROUTINE ProcessLimbData

!=============================================================================
END MODULE THzRadiances
!=============================================================================

! $Log$
! Revision 2.8  2004/12/07 21:24:17  pwagner
! type converted in min intrinsic to appease NAG
!
! Revision 2.7  2004/11/10 15:41:15  perun
! Adjust precision based on bad chan/switch flags; line output records with L1BOA
!
! Revision 2.6  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.5  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.4  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.3  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.2  2003/02/05 21:32:41  perun
! Use variances for precisions
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
