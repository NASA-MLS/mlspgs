! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzRadiances ! Determine radiances for the THz module
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: THzNum, THzChans, LO1
  USE THzCalibration, ONLY : CalBuf, SpaceTemp
  USE MLSL1Rad, ONLY : THzRad, RadPwr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ProcessLimbData

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  REAL, PARAMETER :: LO1R5 = LO1(5)  ! Radiometer 5 1st LO frequency

CONTAINS

  SUBROUTINE CalcLimbRads (nMAF, ibgn)

    USE THzCalibration, ONLY : Kelvins => Cnts  ! Already in Kelvins
    USE MLSL1Config, ONLY: MIFsTHz

    INTEGER, INTENT (IN) :: nMAF
    INTEGER, INTENT (INOUT) :: ibgn

    INTEGER :: iend
    INTEGER :: i, nBank, nChan, MIF_end
    REAL :: T

    iend = ibgn + CalBuf%MAFdata(nMAF-CalBuf%Cal_start+1)%last_MIF
    MIF_end = ibgn + MIFsTHz - 1

    DO nBank = 1, THzNum

       THzRad(nBank)%value = 0.0
       THzRad(nBank)%precision  = -1.0

       IF (CalBuf%BankGood(nBank)) THEN   ! Have good data

          DO i = ibgn, MIF_end
             DO nChan = 1, THzChans
                T = Kelvins(nChan,nBank,i) + SpaceTemp  ! relative to Space
                IF (ABS(T) > 1.0) THEN
                   THzRad(nBank)%value(nChan,i-ibgn+1) = RadPwr (LO1R5, T)
                   THzRad(nBank)%precision(nChan,i-ibgn+1) = 0.0 !???
                ENDIF
            ENDDO
          ENDDO

       ENDIF
    ENDDO
    ibgn = iend + 1   ! For next call

  END SUBROUTINE CalcLimbRads

  SUBROUTINE ProcessLimbData

    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo
    USE OutputL1B, ONLY: OutputL1B_rad

    INTEGER :: MAFno, counterMAF, ibgn
    REAL(r8) :: TAI

    print *, 'ProcessLimbData'

    ibgn = 0
    DO MAFno = CalBuf%Cal_start, CalBuf%Cal_end

       CALL CalcLimbRads (MAFno, ibgn)

       counterMAF = CalBuf%MAFdata(MAFno)%EMAF%TotalMAF
       TAI = CalBuf%MAFdata(MAFno)%SciMIF(0)%secTAI
print *, "Outputting rad for MAFno: ", MAFno
       CALL OutputL1B_rad (MAFno, L1BFileInfo, counterMAF, TAI, THzrad)

! Write Diags (LATER!)

    ENDDO

  END SUBROUTINE ProcessLimbData

!=============================================================================
END MODULE THzRadiances
!=============================================================================

! $Log$
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
