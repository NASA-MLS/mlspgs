! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE GHzBaseline ! Determine GHz baseline radiance corrections
!=============================================================================

  USE MLSL1Common, ONLY: FBchans, L1BFileInfo
  USE MLSL1Config, ONLY: MIFsGHz
  USE MLSL1Rad, ONLY: FBrad
  USE OutputL1B, ONLY: OutputL1B_LatBinData
  USE L1BData, ONLY: L1BData_T, ReadL1BData

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LatBinRads, OutputBaselinedRads

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: NumGHzBands = 15   ! Number of GHz FB bands
  INTEGER, PARAMETER :: NoLatBins = 4

  INTEGER :: LatBinNum(FBchans,NumGHzBands,NoLatBins,2) = 0
  REAL :: LatBinChanAvg(FBchans,NumGHzBands,NoLatBins,2) = 0.0
  REAL :: LatBin(2,NoLatBins,2) = RESHAPE (&
       (/ -90.0, -45.0, -45.0, 0.0, 0.0, 45.0, 45.0, 90.0, &     ! Ascending
          -90.0, -45.0, -45.0, 0.0, 0.0, 45.0, 45.0, 90.0 /), &  ! Descending
       (/ 2, NoLatBins, 2 /))                                    ! final shape
  REAL :: WeightedRadAvg(NumGHzBands,NoLatBins,2) = 0.0

  REAL, PARAMETER :: BaselineAlt(FBchans,NumGHzBands) = 5.0e04   ! test!!

  INTEGER, SAVE :: MAFno = 1

  INTEGER, PARAMETER :: MaxMAFs = 4000   ! more than enough for 24 hours
  INTEGER :: AscDescIndx(MIFsGHz,MaxMAFs), LatBinIndx(MIFsGHz,MaxMAFs)

CONTAINS

!=============================================================================
  SUBROUTINE LatBinRads
!=============================================================================

    USE TkL1B, ONLY: Alt => GHz_GeodAlt, Lat => GHz_GeodLat, &
         Angle => GHz_GeodAngle

    INTEGER :: i, bandindx, bankno, chan, latno, m
    INTEGER :: Aindx, Lindx
    REAL :: ModAngle(MIFsGHz)

print *, 'binning rads by lat...'

! Determine Asc/Desc (1/2) indexes:

    ModAngle = MOD (ABS(Angle), 360.0)
    WHERE (ModAngle > 90.0 .AND. ModAngle <= 270.0)
       AscDescIndx(:,MAFno) = 2   ! Descending index
    ELSEWHERE
       AscDescIndx(:,MAFno) = 1   ! Ascending index
    ENDWHERE

! Determine lat bin indexes:

    DO i = 1, MIFsGHz
       DO latno = 1, NoLatBins
          IF (lat(i) <= LatBin(2,latno,AscDescIndx(i,MAFno))) THEN
             LatBinIndx(i,MAFno) = latno
             EXIT
          ENDIF
       ENDDO
    ENDDO

! Accumulate Radiance data:

    DO i = 1, SIZE (FBrad)

       bandindx = MIN (FBrad(i)%BandNo, NumGHzBands)   ! bandno 21 is last entry
       bankno = FBrad(i)%signal%spectrometernumber

       DO m = 1, MIFsGHz
          DO chan = 1, FBchans
             IF (Alt(m) > BaselineAlt(chan,bandindx) .AND. &
                  FBrad(bankno)%precision(chan,m) >= 0.0) THEN
                Lindx = LatBinIndx(m,MAFno)
                Aindx = AscDescIndx(m,MAFno)
                IF (Lindx > 0 .AND. Aindx > 0) THEN
                   LatBinNum(chan,bandindx,Lindx,Aindx) = &
                        LatBinNum(chan,bandindx,Lindx,Aindx) + 1
                   LatBinChanAvg(chan,bandindx,Lindx,Aindx) = &
                        LatBinChanAvg(chan,bandindx,Lindx,Aindx) &
                        + FBrad(bankno)%value(chan,m)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

    ENDDO

    CALL OutputL1B_LatBinData (MAFno, L1BFileInfo%RADGid, &
         AscDescIndx=AscDescIndx(:,MAFno), LatBinIndx=LatBinIndx(:,MAFno))

    MAFno = MAFno + 1

  END SUBROUTINE LatBinRads

!=============================================================================
  SUBROUTINE CalcWeightedAvgs
!=============================================================================

    USE MLSL1Common, ONLY: BandWidth

    INTEGER :: AscDes, band, bin, chan
    REAL :: NoiseBandWidth(FBchans), SumBandWidth, wghtavg(FBchans)

    print *, 'calc wght avgs...'

    NoiseBandWidth = BandWidth%FB(:,1) * 1.0e-06  ! In MHz

    DO AscDes = 1, 2
       DO bin = 1, NoLatBins
          DO band = 1, NumGHzBands
             wghtavg = 0.0
             SumBandWidth = 0.0
             DO chan = 1, FBchans
                IF (LatBinNum(chan,band,bin,AscDes) > 0) THEN
                   LatBinChanAvg(chan,band,bin,AscDes) = &
                        LatBinChanAvg(chan,band,bin,AscDes) / &
                        LatBinNum(chan,band,bin,AscDes)
                   wghtavg(chan) = LatBinChanAvg(chan,band,bin,AscDes) * &
                        NoiseBandWidth(chan)
                   SumBandWidth = SumBandWidth + NoiseBandWidth(chan)
                ENDIF
             ENDDO
             IF (SumBandWidth > 0.0) &
                  WeightedRadAvg(band,bin,AscDes) = SUM (wghtavg) / SumBandWidth
             LatBinChanAvg(:,band,bin,AscDes) = &
                  LatBinChanAvg(:,band,bin,AscDes) - &
                  WeightedRadAvg(band,bin,AscDes)
          ENDDO
       ENDDO
    ENDDO

    CALL OutputL1B_LatBinData (0, L1BFileInfo%RADGid, &
         LatBinChanAvg=LatBinChanAvg)

    CALL OutputL1B_LatBinData (0, L1BFileInfo%RADGid, BaselineAlt=BaselineAlt)

    CALL OutputL1B_LatBinData (0, L1BFileInfo%RADGid, LatBin=LatBin)

  END SUBROUTINE CalcWeightedAvgs

!=============================================================================
  SUBROUTINE OutputBaselinedRads
!=============================================================================

    USE MLSL1Rad, ONLY: Rad_name
    USE MLSHDF5, ONLY: SaveAsHDF5DS
    USE MLSL1Config, ONLY: L1Config

    REAL, POINTER, DIMENSION(:,:,:) :: rad => NULL(), rad_err => NULL()

    CHARACTER(LEN=80) :: name
    INTEGER :: i, MAF, MaxMIFs, MIF, noMAFs, Flag
    INTEGER :: c1, c2, bandno, bandindx
    INTEGER, PARAMETER :: dims(3) = (/ FBchans, MIFsGHz, 1 /)

    TYPE (L1BData_T) :: L1BData
 
! Calculate weighted averages first:

    CALL CalcWeightedAvgs

    IF (.NOT. L1Config%Output%RemoveBaseline) RETURN  ! Nothing more

! Get the GHz FB rads and errs and adjust

    DO i = 1, 39  ! Check all FB25- possibilities

       name = Rad_name(i)
       CALL ReadL1BData (L1BFileInfo%RADGid, name, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       IF (Flag == 0) THEN
          print *, 'name: ', L1BData%L1BName

          c1 = INDEX (name, '.B') + 2
          c2 = INDEX (name, 'F:') - 1
          READ (name(c1:c2), *) bandno
          bandindx = MIN (bandno, NumGHzBands)   ! bandno 21 is last entry

          IF (.NOT. ASSOCIATED (rad)) THEN
             maxMIFs = L1BData%MaxMIFs
             ALLOCATE (rad(FBchans,MaxMIFs,noMAFs))
             ALLOCATE (rad_err(FBchans,MaxMIFs,noMAFs))
          ENDIF

          rad = L1BData%DpField

          name = TRIM(name) // ' precision'
          CALL ReadL1BData (L1BFileInfo%RADGid, name, L1BData, noMAFs, Flag, &
               NeverFail=.TRUE., HDFversion=5)
          rad_err = L1BData%DpField

! Remove baseline from rads:

          DO MAF = 1, noMAFS
             DO MIF = 1, MaxMIFs
                WHERE (rad_err(:,MIF,MAF) >= 0.0)
                   rad(:,MIF,MAF) = rad(:,MIF,MAF) - &
                       LatBinChanAvg(:,bandindx,LatBinIndx(MIF,MAF), &
                         AscDescIndx(MIF,MAF))
                ENDWHERE
             ENDDO
         ENDDO

! Write adjusted rads here:

          CALL SaveAsHDF5DS (L1BFileInfo%RADGid, TRIM(Rad_Name(i)), rad, &
               adding_to=.TRUE.)

       ENDIF

    ENDDO

  END SUBROUTINE OutputBaselinedRads

!=============================================================================
END MODULE GHzBaseline
!=============================================================================
! $Log$
! Revision 2.2  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
!
