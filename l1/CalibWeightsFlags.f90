! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE CalibWeightsFlags
!=============================================================================

  USE MLSL1Common, ONLY: R8, FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
       DACSnum, DACSchans, MaxMIFs
  USE Calibration, ONLY: WeightsFlags_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DetermineWeightsFlags

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (WeightsFlags_T), DIMENSION(:), ALLOCATABLE :: WeightsFlags

  INTEGER :: nom_MIFs, TotalMAFs

CONTAINS

!=============================================================================
  SUBROUTINE ProcessMAFdata
!=============================================================================

    USE MLSL1Config, ONLY: GHz_seq, GHz_seq_use
    USE MLSL1Common, ONLY: L1BFileInfo
    USE EngTbls, ONLY: EngMAF, EngPkt
    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: SwMirPos, GetScAngles
    USE MLSL1Utils, ONLY: QNan

    INTEGER :: i, ios, sci_MAFno, n, sum_S, sum_T
    INTEGER :: EngMAF_unit, SciMAF_unit, MAF_data_unit
    INTEGER :: mindx(0:(MaxMIFs-1))
    INTEGER, PARAMETER :: GHz_sw_indx(0:(MaxMIFs-1)) = (/ (i, i=1,MaxMIFs) /)
    CHARACTER(len=1) :: GHz_sw_pos(0:(MaxMIFs-1))
    CHARACTER(len=70) :: PLL_DN(0:(MaxMIFs-1)), PLL_dflt_DN(0:(MaxMIFs-1))
    LOGICAL :: recompute
    REAL :: LLO_EU(16,0:(MaxMIFs-1)), LLO_dflt_EU(16,0:(MaxMIFs-1))

    TYPE Sci_pos_T
       REAL :: APE(2,0:(MaxMIFs-1))
       REAL :: ASE(2,0:(MaxMIFs-1))
       REAL :: GME(2,0:(MaxMIFs-1))
       REAL :: TSE(2,0:(MaxMIFs-1))
    END TYPE Sci_pos_T

    TYPE (Sci_pos_T) :: Sci_pos, Sci_dflt_pos

    print *, ' Reading MAFdata...'

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

    i = 0

    outer: DO

       READ (engMAF_unit, iostat=ios) EngMAF, EngPkt
       IF (ios /= 0) EXIT outer
       WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

       READ (sciMAF_unit, iostat=ios) SciMAF
       IF (ios /= 0) EXIT outer
       sci_MAFno = SciMAF(0)%MAFno
       DO n = 0, (MaxMIFs - 1)
          Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
          Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
          Sci_pos%GME(:,n) = SciMAF(n)%GSA_pos
          Sci_pos%TSE(:,n) = SciMAF(n)%TSSA_pos
          PLL_DN(n) = SciMAF(n)%PLL_DN
          LLO_EU(:,n) = SciMAF(n)%LLO_EU
       ENDDO

       DO
          IF (EngMAF%MAFno == sci_MAFno) THEN
             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_EU
             EXIT   ! Nothing more to read
          ENDIF

          IF (EngMAF%secTAI <= SciMAF(2)%secTAI) THEN       ! Catch the Sci

             WRITE (L1BFileInfo%EngId, iostat=ios) Sci_dflt_pos
             WRITE (L1BFileInfo%EngId, iostat=ios) PLL_dflt_DN
             WRITE (L1BFileInfo%EngId, iostat=ios) LLO_dflt_EU
             READ (engMAF_unit, iostat=ios) EngMAF, EngPkt
             IF (ios /= 0) EXIT outer
             WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt, EngMAF%Eng%value

          ELSE IF (EngMAF%secTAI > SciMAF(2)%secTAI) THEN   ! Catch the Eng

             READ (sciMAF_unit, iostat=ios) SciMAF
             IF (ios /= 0) EXIT outer
             sci_MAFno = SciMAF(0)%MAFno
             DO n = 0, (MaxMIFs - 1)
                Sci_pos%APE(:,n) = SciMAF(n)%APE_pos
                Sci_pos%ASE(:,n) = SciMAF(n)%ASA_pos
                Sci_pos%GME(:,n) = SciMAF(n)%GSA_pos
                Sci_pos%TSE(:,n) = SciMAF(n)%TSSA_pos
             ENDDO

          ENDIF
       ENDDO

! Reset GHz SW mirror pos:

       DO n = 0, (MaxMIFs - 1)
          SciMAF(n)%GHz_sw_pos = SwMirPos ("G", SciMAF(n)%GSA_pos)
          GHz_sw_pos(n) = SciMAF(n)%GHz_sw_pos
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

print *, 'i, MAFno, GHz_sw_pos: ', i, EngMAF%MAFno, SciMAF%GHz_sw_pos
print *, 'flags: ', WeightsFlags(i)

       WRITE (MAF_data_unit) WeightsFlags(i)
       WRITE (MAF_data_unit) EngMAF
       WRITE (MAF_data_unit) SciMAF
    ENDDO outer

    TotalMAFs = i
    print *, 'MAFs: ', TotalMAFs
    ENDFILE MAF_data_unit

  END SUBROUTINE ProcessMAFdata

!=============================================================================
  SUBROUTINE DetermineWeightsFlags
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE L1LogUtils, ONLY: EngMAFs, SciMAFs

    INTEGER :: maxMAFs

    print *, 'Determining weights flags...'

    nom_MIFs = L1Config%Calib%MIFsPerMAF

    maxMAFs = MAX (EngMAFs, SciMAFs)
    ALLOCATE (WeightsFlags(maxMAFs))

    CALL ProcessMAFdata

  END SUBROUTINE DetermineWeightsFlags

!=============================================================================
END MODULE CalibWeightsFlags
!=============================================================================
! $Log$
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
