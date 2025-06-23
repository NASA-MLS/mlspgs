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
MODULE GHzBaseline ! Determine GHz baseline radiance corrections
!=============================================================================

  USE MLSL1Common, ONLY: FBchans, MBchans, L1BFileInfo, NumBands
  USE MLSL1Config, ONLY: MIFsGHz
  USE MLSL1Rad, ONLY: L1Brad
  USE OutputL1B, ONLY: OutputL1B_LatBinData
  USE L1BData, ONLY: L1BData_T, ReadL1BData, DeallocateL1BData

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LatBinRads, OutputBaselinedRads

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  INTEGER, PARAMETER :: NumGHzBands = 20   ! Number of GHz FB+MB bands
  INTEGER, PARAMETER :: NoLatBins = 8

  INTEGER, PARAMETER :: BaselineBandNo(NumGHzBands) = (/ &
       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 21, 27, 28, 29, 30, 31 /)
  INTEGER :: LatBinNum(FBchans,NumBands,NoLatBins) = 0
  REAL :: LatBinChanAvg(FBchans,NumBands,NoLatBins) = 0.0
  REAL, PARAMETER :: LatBin(2,NoLatBins) = RESHAPE (&
       (/ -90.0, -45.0, -45.0,  0.0,   0.0, 45.0,  45.0,  90.0, &     ! Ascend
           45.0,  90.0,   0.0, 45.0, -45.0,  0.0, -90.0, -45.0 /), &  ! Descend
      (/ 2, NoLatBins /))                                    ! final shape
  REAL :: WeightedRadAvg(NumBands,NoLatBins) = 0.0

  REAL, SAVE :: BaselineAlt(FBchans,NumBands) = 8.0e04   ! initial values

  INTEGER, SAVE :: MAFno = 1

  INTEGER, PARAMETER :: MaxMAFs = 4000   ! more than enough for 24 hours
  INTEGER :: AscDescIndx(MIFsGHz,MaxMAFs), LatBinIndx(MIFsGHz,MaxMAFs)

CONTAINS

!=============================================================================
  SUBROUTINE LatBinRads
!=============================================================================

    USE TkL1B, ONLY: Alt => GHz_GeodAlt, Lat => GHz_GeodLat, &
         Angle => GHz_GeodAngle
    USE SpectralBaseline, ONLY: BaselineInclude, BaselineAC

    INTEGER :: i, bandno, bankno, chan, latno, m, nchans
    INTEGER :: Aindx, Lindx
    REAL :: ModAngle(MIFsGHz)

PRINT *, 'binning rads by lat...'

! Set Baseline Alts, if not already done:

    IF (MAFno == 1) THEN
       DO i = 1, (NumGHzBands - 1)    ! Band 21 is switchable
          bandno = L1Brad(i)%BandNo
          IF (bandno <= 14) THEN
             nchans = FBchans
          ELSE IF (bandno == 21) THEN
             nchans = FBchans
          ELSE
             nchans = MBchans
          ENDIF
          DO chan = 1, nchans
             IF (.NOT. BaselineInclude(bandno)%chan(chan)) &
                  BaselineAlt(chan,bandno) = HUGE (1.0)
          ENDDO
       ENDDO
    ENDIF

! Determine Asc/Desc (1/2) indexes:

    ModAngle = MOD ((Angle + 360.0), 360.0)
    WHERE (ModAngle > 90.0 .AND. ModAngle <= 270.0)
       AscDescIndx(:,MAFno) = 2   ! Descending index
    ELSEWHERE
       AscDescIndx(:,MAFno) = 1   ! Ascending index
    ENDWHERE

! Determine lat bin indexes:

    DO i = 1, MIFsGHz
       IF (AscDescIndx(i,MAFno) == 1) THEN     ! Ascending
          DO latno = 1, NoLatBins/2
             IF (lat(i) <= LatBin(2,latno)) THEN
                LatBinIndx(i,MAFno) = latno
                EXIT
             ENDIF
          ENDDO
       ELSE                                    ! Descending
          DO latno = (NoLatBins/2 + 1), NoLatBins
             IF (lat(i) >= LatBin(1,latno)) THEN
                LatBinIndx(i,MAFno) = latno
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO

! Accumulate Radiance data:

    DO i = 1, (NumGHzBands - 1)    ! Band 21 is switchable

       bankno = L1Brad(i)%signal%spectrometernumber
       bandno = L1Brad(i)%BandNo
       IF (bandno <= 14) THEN
          nchans = FBchans
       ELSE IF (bandno == 21) THEN
          nchans = FBchans
       ELSE
          nchans = MBchans
          bankno = bankno + 14    ! relative to start of FB data
       ENDIF

       DO m = 1, MIFsGHz
          DO chan = 1, nchans
             IF (Alt(m) > BaselineAlt(chan,bandno) .AND. &
                  L1Brad(bankno)%precision(chan,m) >= 0.0) THEN
                Lindx = LatBinIndx(m,MAFno)
                Aindx = AscDescIndx(m,MAFno)
                IF (Lindx > 0) THEN
                   LatBinNum(chan,bandno,Lindx) = &
                        LatBinNum(chan,bandno,Lindx) + 1
                   LatBinChanAvg(chan,bandno,Lindx) = &
                        LatBinChanAvg(chan,bandno,Lindx) &
                        + L1Brad(bankno)%value(chan,m) - &
                        BaselineAC(bandno)%offset(chan)
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

    INTEGER :: bandno, bin, chan, i, nchans, NBWoffset
    REAL :: NoiseBandWidth(FBchans), SumBandWidth, wghtavg(FBchans)

    PRINT *, 'calc wght avgs...'

    NoiseBandWidth = BandWidth%FB(:,1) * 1.0e-06  ! In MHz

    DO bin = 1, NoLatBins
       DO i = 1, NumGHzBands
          bandno = BaselineBandNo(i)
          wghtavg = 0.0
          SumBandWidth = 0.0
          IF (bandno <= 15) THEN
             NBWoffset = 0
             nchans = FBchans
          ELSE
             NBWoffset = 7
             nchans = MBchans
          ENDIF

          DO chan = 1, nchans
             IF (LatBinNum(chan,bandno,bin) > 0) THEN
                LatBinChanAvg(chan,bandno,bin) = &
                     LatBinChanAvg(chan,bandno,bin) / &
                     LatBinNum(chan,bandno,bin)
                wghtavg(chan) = LatBinChanAvg(chan,bandno,bin) * &
                     NoiseBandWidth(chan+NBWoffset)
                SumBandWidth = SumBandWidth + NoiseBandWidth(chan+NBWoffset)
             ENDIF
          ENDDO
          IF (SumBandWidth > 0.0) &
               WeightedRadAvg(bandno,bin) = SUM (wghtavg) / &
               SumBandWidth
          WHERE (LatBinNum(:,bandno,bin) > 0)
             LatBinChanAvg(:,bandno,bin) = &
                  LatBinChanAvg(:,bandno,bin) - &
                  WeightedRadAvg(bandno,bin)
          ENDWHERE
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
    USE HDF5, ONLY: H5gClose_f, H5gOpen_f
    USE MLSHDF5, ONLY: SaveAsHDF5DS, MakeHDF5Attribute
    USE MLSL1Config, ONLY: L1Config

    REAL, POINTER, DIMENSION(:,:) :: baseline => NULL()
    LOGICAL, POINTER, DIMENSION(:) :: binMinus => NULL()

    CHARACTER(LEN=80) :: name, binnedname
    INTEGER :: i, MAF, MIF, noMAFs, Flag, nchans, grp_id
    INTEGER :: c1, c2, bandno, status

    TYPE (L1BData_T) :: L1BData

    binMinus => L1Config%Output%SubtractBinnedBaseline
 
! Calculate weighted averages first:

    CALL CalcWeightedAvgs

    IF (.NOT. L1Config%Output%RemoveBaseline) RETURN  ! Nothing more

! Get the GHz FB/MB baselines and adjust

    DO i = 1, SIZE (Rad_name)  ! Check all FB25- and MB11- possibilities

       name = TRIM(Rad_name(i)) // ' Baseline'
       binnedname = TRIM(Rad_name(i)) // ' BinnedBaselineAC'
       c1 = INDEX (name, '.B') + 2
       IF (i < 18) THEN
          c2 = c1 
       ELSE
          c2 = c1 + 1
       ENDIF
       READ (name(c1:c2), *) bandno

       IF ((bandno > 14 .AND. bandno < 21) .OR. (bandno > 21 .AND. &
            bandno < 27) .OR. (bandno > 31)) CYCLE   ! nothing for other bands

       CALL ReadL1BData (L1BFileInfo%RADGid, name, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       IF (Flag == 0) THEN

          IF (bandno <= 14) THEN
             nchans = FBchans
          ELSE IF (bandno == 21) THEN
             nchans = FBchans
          ELSE
             nchans = MBchans
          ENDIF

          IF (.NOT. ASSOCIATED (baseline)) THEN
             ALLOCATE (baseline(nchans,noMAFs))
          ENDIF

          baseline = L1BData%DpField(1,:,:)
          CALL DeallocateL1BData (L1BData)

          IF (binMinus(bandno)) THEN

! Remove baseline from baseline:

             MIF = 1     ! use start of MAF for index
             DO MAF = 1, noMAFS
                baseline(:,MAF) = baseline(:,MAF) - &
                     LatBinChanAvg(:,bandno,LatBinIndx(MIF,MAF))
             ENDDO

! Write adjusted baselines here:

             CALL SaveAsHDF5DS (L1BFileInfo%RADGid, name, baseline, &
                  adding_to=.TRUE.)

          ENDIF

          DEALLOCATE (baseline, stat=status)

! Write appropriate baseline adjustments:

          CALL OutputL1B_LatBinData (0, L1BFileInfo%RADGid, &
               BinnedBaseline=LatBinChanAvg(1:nchans,bandno,:), &
               Name=binnedname)

       ENDIF

    ENDDO

! Save baseline flags (may want to make as DS with attributes):

    CALL H5gOpen_f (L1BFileInfo%RADGid, '/', grp_id, i)
    CALL MakeHDF5Attribute (grp_id, 'BinnedBaselineSubtracted', &
         binMinus, .TRUE.)
    CALL H5gClose_f (grp_id, i)

  END SUBROUTINE OutputBaselinedRads

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE GHzBaseline
!=============================================================================
! $Log$
! Revision 2.6  2006/08/02 18:54:08  perun
! Accumulate latitude bin data in 8 bins instead of 4x2 bins
!
! Revision 2.5  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2005/05/02 16:02:50  perun
! Deallocate last rad and rad_err pointers
!
! Revision 2.3  2004/11/10 15:35:10  perun
! Add call to deallocate L1BData
!
! Revision 2.2  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
!
