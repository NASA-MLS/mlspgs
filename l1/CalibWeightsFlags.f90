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
MODULE CalibWeightsFlags
!=============================================================================

  USE MLSL1Common, ONLY: MaxMIFs, R8, DACSnum
  USE Calibration, ONLY: WeightsFlags_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DetermineWeightsFlags

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (WeightsFlags_T), DIMENSION(:), POINTER :: WeightsFlags

  INTEGER :: nom_MIFs, TotalMAFs

CONTAINS

!=============================================================================
  SUBROUTINE ProcessMAFdata
!=============================================================================

    USE MLSL1Config, ONLY: GHz_seq
    USE MLSL1Common, ONLY: L1BFileInfo
    USE EngTbls, ONLY: EngMAF, EngPkt
    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: SwMirPos, GetScAngles
    USE MLSL1Utils, ONLY: QNan
    USE DACsUtils, ONLY: TPz
    USE L1LogUtils, ONLY: MAF_dur

    INTEGER :: i, ios, sci_MAFno, n, sum_S, sum_T
    INTEGER :: EngMAF_unit, SciMAF_unit, MAF_data_unit
    INTEGER :: mindx(0:(MaxMIFs-1)), ngood(DACSnum)
    INTEGER, PARAMETER :: GHz_sw_indx(0:(MaxMIFs-1)) = (/ (i, i=1,MaxMIFs) /)
    CHARACTER(len=1) :: GHz_sw_pos(0:(MaxMIFs-1))
    CHARACTER(len=70) :: PLL_DN(0:(MaxMIFs-1)), PLL_dflt_DN(0:(MaxMIFs-1))
    REAL :: LLO_EU(16,0:(MaxMIFs-1)), LLO_dflt_EU(16,0:(MaxMIFs-1))
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)
    REAL :: swFac(0:MaxMIFno)
    REAL(r8) :: TP_ana(0:MaxMIFno,DACSnum), TP_dig(0:MaxMIFno,DACSnum), &
         engTAI, sciTAI
    REAL(r8) :: Sum_ana(DACSnum) = 0.0, Sum_dig(DACSnum) = 0.0, &
         Sum_dig_dig(DACSnum) = 0.0, Sum_dig_ana(DACSnum) = 0.0
    CHARACTER(len=1), PARAMETER :: discard = "D"
    REAL :: TP_dig_max = 4.0     ! Maximum allowable TP digital

    TYPE Sci_pos_T
       REAL :: APE(2,0:(MaxMIFs-1))
       REAL :: ASE(2,0:(MaxMIFs-1))
       REAL :: GME(2,0:(MaxMIFs-1))
       REAL :: TSE(2,0:(MaxMIFs-1))
    END TYPE Sci_pos_T

    TYPE (Sci_pos_T) :: Sci_pos, Sci_dflt_pos

    PRINT *, ' Reading MAFdata...'

! Init default positions:

    Sci_dflt_pos%APE = QNan()
    Sci_dflt_pos%ASE = QNan()
    Sci_dflt_pos%GME = QNan()
    Sci_dflt_pos%TSE = QNan()

    DO n = 0, (MaxMIFs - 1)
       PLL_dflt_DN(n) = CHAR(0)
    ENDDO

    LLO_dflt_EU = QNan()

    MAF_data_unit = L1BFileInfo%MAF_data_unit
    EngMAF_unit = L1BFileInfo%EngMAF_unit
    SciMAF_unit = L1BFileInfo%SciMAF_unit

    REWIND (engMAF_unit)
    REWIND (sciMAF_unit)

! Determine override sequence indexes to use for comparisons

    mindx = 0
    WHERE (GHz_seq == "S")
       mindx = GHz_sw_indx
    ENDWHERE
    sum_S = SUM (mindx)

    mindx = 0
    WHERE (GHz_seq == "T")
       mindx = GHz_sw_indx
    ENDWHERE
    sum_T = SUM (mindx)

    ngood = 0
    i = 0

    outer: DO

       READ (engMAF_unit, iostat=ios) EngMAF, EngPkt
       engTAI = EngMAF%secTAI
       IF (ios /= 0) EXIT outer
       WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

       READ (sciMAF_unit, iostat=ios) SciMAF
       sciTAI = sciMAF(1)%secTAI
       IF (ios /= 0) EXIT outer
       sci_MAFno = SciMAF(0)%MAFno

       DO n = 0, (MaxMIFs - 1)
          Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
          Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
          Sci_pos%GME(:,n) = SciMAF(n)%GSM_pos
          Sci_pos%TSE(:,n) = SciMAF(n)%TSSM_pos
          PLL_DN(n) = SciMAF(n)%PLL_DN
          LLO_EU(:,n) = SciMAF(n)%LLO_EU
       ENDDO

       DO

          IF (EngMAF%MAFno == sci_MAFno .AND. &
               NINT(ABS(sciTAI - engTAI) / MAF_dur) == 1) THEN  ! Same MAF data
             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_EU
             EXIT   ! Nothing more to read
          ENDIF

          IF (engTAI <= sciTAI) THEN       ! Catch the Sci

             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_dflt_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_dflt_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_dflt_EU
             READ (engMAF_unit, iostat=ios) EngMAF, EngPkt
             IF (ios /= 0) EXIT outer
             engTAI = EngMAF%secTAI
             WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

          ELSE ! Catch the Eng

             READ (sciMAF_unit, iostat=ios) SciMAF
             sciTAI = sciMAF(2)%secTAI
             IF (ios /= 0) EXIT outer
             sci_MAFno = SciMAF(0)%MAFno
             DO n = 0, (MaxMIFs - 1)
                Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
                Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
                Sci_pos%GME(:,n) = SciMAF(n)%GSM_pos
                Sci_pos%TSE(:,n) = SciMAF(n)%TSSM_pos
             ENDDO

          ENDIF
       ENDDO

! Reset GHz SW mirror pos:

       DO n = 0, (MaxMIFs - 1)
          SciMAF(n)%GHz_sw_pos = SwMirPos ("G", SciMAF(n)%GSM_theta)
          GHz_sw_pos(n) = SciMAF(n)%GHz_sw_pos
       ENDDO

! Accumulate DACS TPs to calculate daily TP zeros:

       swFac = 1.0        ! init switch factor to non-discard MIFs
       WHERE (SciMAF%GHz_sw_pos == discard)   ! Don't use discards
          swFac = 0.0
       ENDWHERE
!       ngood = ngood + COUNT (swFac == 1.0)
       DO n = 1, DACSNUM
          TP_ana(:,n) = SciMAF%TP(n) * swFac
          TP_dig(:,n) = (SciMAF%TPdigP(n) + SciMAF%TPdigN(n)) * swFac
          WHERE (TP_dig(:,n) > TP_dig_max)
             TP_dig(:,n) = 0.0
             TP_ana(:,n) = 0.0
          ENDWHERE
          ngood(n) = ngood(n) + COUNT (TP_dig(:,n) /= 0.0)
       ENDDO
       DO n = 1, DACSNUM
          Sum_dig(n) = Sum_dig(n) + SUM (TP_dig(:,n))
          Sum_ana(n) = Sum_ana(n) + SUM (TP_ana(:,n))
          Sum_dig_dig(n) = Sum_dig_dig(n) + SUM (TP_dig(:,n)*TP_dig(:,n))
          Sum_dig_ana(n) = Sum_dig_ana(n) + SUM (TP_dig(:,n)*TP_ana(:,n))
       ENDDO

! Reset scAngles to be used for L1BOA

       CALL GetScAngles

       i = i + 1

       WeightsFlags(i)%MAFno = EngMAF%MAFno

! Check MAF duration

       WeightsFlags%recomp_MAF  = (EngMAF%MIFsPerMAF /= nom_MIFs)

! Check switching mirror positions:

       mindx = 0
       WHERE (GHz_seq == "S" .AND. GHz_sw_pos /= "D")
          mindx = GHz_sw_indx
       ENDWHERE
       WeightsFlags%recomp_S  = (SUM (mindx) /= sum_S)

       mindx = 0
       WHERE (GHz_seq == "T" .AND. GHz_sw_pos /= "D")
          mindx = GHz_sw_indx
       ENDWHERE
       WeightsFlags%recomp_T  = (SUM (mindx) /= sum_T)

PRINT *, 'i, MAFno, GHz_sw_pos: ', i, EngMAF%MAFno, SciMAF%GHz_sw_pos
PRINT *, 'flags: ', WeightsFlags(i)

       WRITE (MAF_data_unit) WeightsFlags(i)
       WRITE (MAF_data_unit) EngMAF
       WRITE (MAF_data_unit) SciMAF
    ENDDO outer

    TotalMAFs = i
    PRINT *, 'MAFs: ', TotalMAFs
    ENDFILE MAF_data_unit

! Calculate TPz values to pass on the MLSL1G program:

    DO i = 1, DACSNUM
       TPz(i) = ((Sum_dig_dig(i) / ngood(i) * Sum_ana(i) / ngood(i)) - &
            (Sum_dig_ana(i) / ngood(i) * Sum_dig(i) / ngood(i))) / &
            (Sum_dig_dig(i) / ngood(i) - (Sum_dig(i) / ngood(i))**2)
    ENDDO

  END SUBROUTINE ProcessMAFdata

!=============================================================================
  SUBROUTINE DetermineWeightsFlags
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE L1LogUtils, ONLY: EngMAFs, SciMAFs

    INTEGER :: maxMAFs

    PRINT *, 'Determining weights flags...'

    nom_MIFs = L1Config%Calib%MIFsPerMAF

    maxMAFs = MAX (EngMAFs, SciMAFs)
    ALLOCATE (WeightsFlags(maxMAFs))

    CALL ProcessMAFdata

  END SUBROUTINE DetermineWeightsFlags

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE CalibWeightsFlags
!=============================================================================
! $Log$
! Revision 2.10  2011/01/27 15:33:53  perun
! Change TAI time from MIF 2 to MIF 1
!
! Revision 2.9  2007/06/21 20:57:53  perun
! Exclude TP_dig greater than 4.0 in calculating TPz
!
! Revision 2.8  2006/09/26 16:01:33  perun
! Correct testing to insure alignment of Eng and Sci data
!
! Revision 2.7  2006/08/02 18:52:17  perun
! Accumulate DACS TPs to calculate daily TP zeros
!
! Revision 2.6  2006/04/05 18:10:08  perun
! Remove unused variables
!
! Revision 2.5  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.3  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.2  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.1  2003/08/15 14:54:55  perun
! Version 1.2 commit
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
!
