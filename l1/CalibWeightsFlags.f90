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
    USE EngTbls, ONLY: EngMAF
    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: SwMirPos, GetScAngles

    INTEGER :: i, ios, sci_MAFno, n, sum_S, sum_T
    INTEGER :: EngMAF_unit, SciMAF_unit, MAF_data_unit
    INTEGER :: mindx(0:(MaxMIFs-1))
    INTEGER, PARAMETER :: GHz_sw_indx(0:(MaxMIFs-1)) = (/ (i, i=1,MaxMIFs) /)
    CHARACTER(len=1) :: GHz_sw_pos(0:(MaxMIFs-1))
    LOGICAL :: recompute

    print *, ' Reading MAFdata...'

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

       READ (engMAF_unit, iostat=ios) EngMAF
       IF (ios /= 0) EXIT outer

       READ (sciMAF_unit, iostat=ios) SciMAF
       IF (ios /= 0) EXIT outer
       sci_MAFno = SciMAF(0)%MAFno

       DO
          IF (EngMAF%MAFno == sci_MAFno) EXIT   ! Nothing more to read

          IF (EngMAF%secTAI <= SciMAF(2)%secTAI) THEN       ! Catch the Sci

             READ (engMAF_unit, iostat=ios) EngMAF
             IF (ios /= 0) EXIT outer

          ELSE IF (EngMAF%secTAI > SciMAF(2)%secTAI) THEN   ! Catch the Eng

             READ (sciMAF_unit, iostat=ios) SciMAF
             IF (ios /= 0) EXIT outer
             sci_MAFno = SciMAF(0)%MAFno

          ENDIF
       ENDDO

! Reset GHz SW mirror pos:

       DO n = 0, (MaxMIFs - 1)
          SciMAF(n)%GHz_sw_pos = SwMirPos ("G", SciMAF(n)%GME_pos)
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
! Revision 2.1  2003/08/15 14:54:55  perun
! Version 1.2 commit
!
! Revision 2.1  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
!
