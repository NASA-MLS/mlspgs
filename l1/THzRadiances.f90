! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzRadiances ! Determine radiances for the THz module
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: THzNum, THzChans
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
    USE MLSL1Config, ONLY: MIFsTHz

    INTEGER, INTENT (IN) :: nMAF
    INTEGER, INTENT (INOUT) :: ibgn

    INTEGER :: i, iend, mindx, nBank, nChan, MIF_end

    iend = ibgn + CalBuf%MAFdata(nMAF-CalBuf%Cal_start+1)%last_MIF
    MIF_end = ibgn + MIFsTHz - 1

    DO nBank = 1, THzNum

       THzRad(nBank)%value = 0.0
       THzRad(nBank)%precision  = -1.0

       IF (CalBuf%BankGood(nBank)) THEN   ! Have good data

          DO i = ibgn, MIF_end
             DO nChan = 1, THzChans
                mindx = i-ibgn+1
                THzRad(nBank)%value(nChan,mindx) = Kelvins(nChan,nBank,i)
                IF (THzRad(nBank)%value(nChan,mindx) > -100.0) &
                     THzRad(nBank)%precision(nChan,mindx) = &
                     VarK(nChan,nBank,i)
             ENDDO
          ENDDO

       ENDIF
    ENDDO
    ibgn = iend + 1   ! For next call

  END SUBROUTINE CalcLimbRads

  SUBROUTINE ProcessLimbData

    USE MLSL1Common, ONLY: L1BFileInfo
    USE OutputL1B, ONLY: OutputL1B_rad, OutputL1B_DiagsT
    USE EngTbls, ONLY: Reflec_T

    INTEGER :: MAFno, counterMAF, ibgn, nv
    REAL(r8) :: TAI
    TYPE (Reflec_T) :: Reflec
    INTEGER, SAVE :: OrbNo = 1

    print *, 'ProcessLimbData'

    nv = 1
    ibgn = 0
    DO MAFno = CalBuf%Cal_start, CalBuf%Cal_end

       CALL CalcLimbRads (MAFno, ibgn)

       counterMAF = CalBuf%MAFdata(MAFno)%EMAF%TotalMAF
       TAI = CalBuf%MAFdata(MAFno)%SciMIF(0)%secTAI

print *, "Outputting rad for MAFno: ", MAFno
       CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, Reflec, TAI, THzrad)

! Write MAF dimensioned Diags

       CALL OutputL1B_DiagsT (L1BFileInfo%DiagTid, MAFno=MAFno, &
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
