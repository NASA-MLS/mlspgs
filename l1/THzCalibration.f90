! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE THzCalibration ! Calibration data and routines for the THz module
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: MaxMIFs, THzChans, THzNum
  USE L0_sci_tbls, ONLY: THz_Sci_pkt_T
  USE EngTbls, ONLY : Eng_MAF_T
 
  IMPLICIT NONE

  SAVE

  PRIVATE

  PUBLIC :: CalibrateTHz, Chan_type_T, MAFdata_T, CalBuf_T, CalBuf, Cnts, &
       VarCnts, SpaceTemp

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(THzChans,THzNum)      ! standard THz filter banks
  END TYPE Chan_type_T

  !! THz Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (THz_Sci_pkt_T) :: SciMIF(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     REAL :: CalTgtTemp   ! Average Calibration Target Temperature (C)
     INTEGER :: last_MIF
     INTEGER :: BandSwitch(5)   ! band switch positions
  END TYPE MAFdata_T

  TYPE CalBuf_T
     INTEGER :: MAFs, Cal_start, Cal_end
     LOGICAL :: BankGood(THzNum)
     TYPE (MAFdata_T), DIMENSION(:), ALLOCATABLE :: MAFdata
  END TYPE CalBuf_T

  TYPE (CalBuf_T), TARGET :: CalBuf

  INTEGER :: Cal_start = 1
  INTEGER :: Cal_end, ntot
  INTEGER, DIMENSION(:), ALLOCATABLE :: CalMIF
  INTEGER, DIMENSION(:), ALLOCATABLE :: MIFx
  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: CalFlag
  INTEGER, DIMENSION(:), ALLOCATABLE :: BoundsX
  REAL(r8), DIMENSION(:), ALLOCATABLE :: Bias
  REAL(r8), DIMENSION(:), ALLOCATABLE :: CalTemp
  REAL(r8), DIMENSION(:,:), ALLOCATABLE :: VarGain
  REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Cnts, dCnts, VarCnts
  LOGICAL, DIMENSION(:), ALLOCATABLE :: BiasGood
  REAL :: SpaceTemp
  REAL(r8) :: yTsys(THzChans,THzNum)

CONTAINS

!=============================================================================
  SUBROUTINE BuildCalVectors
!=============================================================================

    USE THzUtils, ONLY: Bias_err, MaxBias
    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Rad, ONLY: UpdateRadSignals

    INTEGER :: i, j, last_MIF, mindx, MIF0, MIFno
    INTEGER :: nBank, numMIFs, nBounds, status, Switch5
    TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()
    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: Bounds
    INTEGER, PARAMETER :: Cal_size = 240   ! Nominal orbit

    SpaceTemp = L1Config%Calib%THzSpaceTemp

    Cal_end = MIN ((Cal_start+Cal_size-1), CalBuf%MAFs)
    IF ((CalBuf%MAFs - Cal_end) < Cal_size) Cal_end = CalBuf%MAFs

! Figure out how many MIFs needed for calibration:

    numMIFs = 0
    DO i = Cal_start, Cal_end
       numMIFs = numMIFs + CalBuf%MAFdata(i)%last_MIF + 1
    ENDDO

! Allocate vectors:

    DEALLOCATE (CalMIF, stat=status)
    ALLOCATE (CalMIF(0:numMIFs-1))

    DEALLOCATE (MIFx, stat=status)
    ALLOCATE (MIFx(0:numMIFs-1))
    DO i = 0, (numMIFs - 1)
       MIFx(i) = i
    ENDDO

    DEALLOCATE (CalFlag, stat=status)  ! calibration = 1, otherwise = 0
    ALLOCATE (CalFlag(0:numMIFs-1))
    CalFlag = 0    ! Indicate not a calibration MIF (yet)

    DEALLOCATE (Bias, stat=status)
    ALLOCATE (Bias(0:numMIFs-1))

    DEALLOCATE (BiasGood, stat=status)
    ALLOCATE (BiasGood(0:numMIFs-1))

    DEALLOCATE (CalTemp, stat=status)
    ALLOCATE (CalTemp(0:numMIFs-1))
    CalTemp = 0.0

    DEALLOCATE (Bounds, stat=status)
    ALLOCATE (Bounds(0:numMIFs-1))

    DEALLOCATE (Cnts, stat=status)
    ALLOCATE (Cnts(THzChans,THzNum,0:numMIFs-1))
    DEALLOCATE (dCnts, stat=status)
    ALLOCATE (dCnts(THzChans,THzNum,0:numMIFs-1))
    DEALLOCATE (VarCnts, stat=status)
    ALLOCATE (VarCnts(THzChans,THzNum,0:numMIFs-1))
    DEALLOCATE (VarGain, stat=status)
    ALLOCATE (VarGain(THzChans,THzNum))
    Cnts = 0.0
    VarCnts = 1.0

! Fill the vectors:

    CalBuf%BankGood = .TRUE.
    Switch5 = CalBuf%MAFdata(Cal_start)%SciMIF(0)%BandSwitch(5) ! init compare

    MIF0 = 0   ! index for MIFs
    DO i = Cal_start, Cal_end
       CurMAFdata => CalBuf%MAFdata(i)
       last_MIF = CurMAFdata%last_MIF
       CalMIF(MIF0:MIF0+last_MIF) = CurMAFdata%SciMIF(0:last_MIF)%MIFno
       Bias(MIF0:MIF0+last_MIF) = &
            CurMAFdata%SciMIF(0:last_MIF)%LLO_Bias
       DO MIFno = 0, last_MIF
          mindx = MIF0 + MIFno
          Cnts(:,:,mindx) = &
               CurMAFdata%SciMIF(MIFno)%FB(:,:)
          IF (CurMAFdata%SciMIF(MIFno)%SwMirPos == "S" .AND. &
               Bias(mindx) < MaxBias) THEN
             CalTemp(mindx) = 0
             CalFlag(mindx) = 1
          ELSE IF (CurMAFdata%SciMIF(MIFno)%SwMirPos == "T" .AND. &
               Bias(mindx) < MaxBias) THEN
             CalTemp(mindx) = -25  !CurMAFdata%CalTgtTemp - SpaceTemp
             CalFlag(mindx) = 1
          ENDIF
      ENDDO

! Check if Band Switch 5 changes:

      IF (ANY (CurMAFdata%SciMIF(0:last_MIF)%BandSwitch(5) /= Switch5)) &
           CalBuf%BankGood(1) = .FALSE.  ! Filter Bank 15 not all good

! Check if Band 20 is connected:

      IF (ANY (CurMAFdata%SciMIF(0:last_MIF)%BandSwitch(4) /= 20)) &
           CalBuf%BankGood(6) = .FALSE.  ! Filter Bank 12 (Band 20) not all good

      MIF0 = MIF0 + last_MIF + 1
    ENDDO

! Adjust LLO Bias V if 2 or more consecutive errs:

    DO i = (numMIFs - 2), 1, -1
       IF (Bias(i) == Bias_err .AND. Bias(i-1) == Bias_err) THEN
          Bias(i+1) = Bias_err
          CalFlag(i+1) = 0
       ENDIF
    ENDDO

! Determine bounds for fitting:

    Bounds = .FALSE.   ! nothing yet
    WHERE (Bias >= MaxBias) ! Re-optimizing
       Bounds = .TRUE.
    END WHERE
    WHERE (CalMIF < 0)      ! Missing MIFs
       Bounds = .TRUE.
    END WHERE
    BiasGood = .NOT. Bounds      ! Good Bias values
    Bounds(numMIFs-1) = .TRUE.   ! last always a bounds
    nBounds = COUNT (Bounds)
    DEALLOCATE (BoundsX, stat=status)
    ALLOCATE (BoundsX(nBounds))
    j = 1
    DO i = 0, (numMIFs - 1)
       IF (Bounds(i)) THEN
          BoundsX(j) = i
          j = j + 1
       ENDIF
    ENDDO

! Determine which filter bank counts are good to cal:

    DO nBank = 1, THzNum
       IF (ALL (Cnts(:,nBank,:) == 0.0)) THEN
          CalBuf%BankGood(nBank) = .FALSE.
       ENDIF
    ENDDO

! Set signals based on end cal index

    CALL UpdateRadSignals (CalBuf%MAFdata(Cal_end)%BandSwitch)

! Save current cal start & end:

    CalBuf%Cal_start = Cal_start
    CalBuf%Cal_end = Cal_end

    Cal_start = Cal_end + 1

  END SUBROUTINE BuildCalVectors

!=============================================================================
  SUBROUTINE THzDel
!=============================================================================

    INTEGER :: iBound, ibgn, iend, nBounds, nsize, ntotx, status
    INTEGER :: nChan, nBank
    INTEGER, DIMENSION(:), POINTER :: cFlag
    REAL(r8), DIMENSION(:), POINTER :: biasx, tempx
    REAL(r8) :: AvgBias, AvgTemp, avgx
    REAL(r8), DIMENSION(:), TARGET, ALLOCATABLE, SAVE :: BiasBuf, TempBuf
    REAL(r8) :: bb, det, tb, tt, yb(THzChans,THzNum), yt(THzChans,THzNum)

    REAL(r8) :: dCal(THzChans,THzNum), dLlo(THzChans,THzNum)

    ntot = COUNT (CalFlag == 1)
    IF (ntot == 0) RETURN

    AvgBias = SUM (Bias * CalFlag) / ntot
    AvgTemp = SUM (CalTemp * CalFlag) / ntot

    DEALLOCATE (BiasBuf, stat=status)
    nsize = SIZE (Bias)
    ALLOCATE (BiasBuf(0:nsize-1))
    BiasBuf = Bias * CalFlag
    DEALLOCATE (TempBuf, stat=status)
    ALLOCATE (TempBuf(0:nsize-1))
    TempBuf = CalTemp * CalFlag

    nBounds = SIZE (BoundsX)
    bb = ntot / (819.2 * 819.2)
    tb = 0.0
    tt = 0.25
    yb = 0.0
    yt = 0.0
    ibgn = 0
    DO iBound = 1, nBounds
       iend = BoundsX(iBound)
       cFlag => CalFlag(ibgn:iend)
       ntotx = COUNT (cFlag == 1)
       IF (ntotx > 1) THEN
          biasx => BiasBuf(ibgn:iend)
          avgx = SUM (biasx) / ntotx
          biasx = (biasx - avgx) * cFlag
          tempx => TempBuf(ibgn:iend)
          avgx = SUM (tempx) / ntotx
          tempx = (tempx - avgx) * cFlag
          bb = bb + SUM (biasx * biasx)
          tb = tb + SUM (tempx * biasx)
          tt = tt + SUM (tempx * tempx)
          DO nBank = 1, THzNum
             DO nChan = 1, THzChans
                yb(nChan,nBank) = yb(nChan,nBank) + &
                     DOT_PRODUCT (Cnts(nChan,nBank,ibgn:iend), biasx)
                yt(nChan,nBank) = yt(nChan,nBank) + &
                     DOT_PRODUCT (Cnts(nChan,nBank,ibgn:iend), tempx)
             ENDDO
          ENDDO
       ENDIF
       ibgn = iend + 1
    ENDDO

    det = tt * bb - tb * tb
    tt = tt / det
    tb = tb / det
    bb = bb / det
    VarGain = bb
    dCal = yt * bb - yb * tb  ! counts per K for T - S
    dLlo = yb * tt - yt * tb  ! counts per V for LLO bias

    WHERE (BiasGood)
       Bias = Bias - AvgBias  ! Adjust bias (OK since it's not needed later)
    ENDWHERE
    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          WHERE (BiasGood)
             Cnts(nChan,nBank,:) = Cnts(nChan,nBank,:) - dLlo(nChan,nBank) * &
                  Bias
             VarCnts(nChan,nBank,:) = Bias * Bias * tt + 1.0
          ENDWHERE
          Cnts(nChan,nBank,:) = Cnts(nChan,nBank,:) / &
               MAX (dCal(nChan,nBank), 0.01d0)
          yTsys(nChan,nBank) = SUM (Cnts(nChan,nBank,:)* CalFlag) / ntot - &
               AvgTemp
       ENDDO
    ENDDO

  END SUBROUTINE THzDel

!=============================================================================
  SUBROUTINE THzCal
!=============================================================================

    INTEGER :: zwin = 2                      ! may pass this in

    INTEGER :: nChan, nBank
    INTEGER :: i, ibgn, iend, iBound, kbgn, maxPt, mbgn, mend, mif0, status
    INTEGER :: nordr, ntotx, wbgn, wend, win
    INTEGER, DIMENSION(:), POINTER :: cFlag
    REAL(r8) :: t, tmid, val
    REAL(r8), DIMENSION(:), ALLOCATABLE, SAVE :: tsq, tval
    REAL(r8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: yval
    REAL(r8) :: av(THzChans,THzNum), bv(THzChans,THzNum), cv(THzChans,THzNum)
    REAL(r8) :: det, mfit0, mfit1, mfit2, mfit3, mfit4, tsqavg
    REAL(r8) :: yfit1(THzChans,THzNum), yfit2(THzChans,THzNum)

    INTEGER, PARAMETER :: MIFsPerMAF = 148   ! use for now
    INTEGER, PARAMETER :: mlimb = 122
    LOGICAL :: last_fit

    print *, 'THzCal'

    IF (COUNT (CalFlag == 1) == 0) RETURN  ! Nothing to calibrate

    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          dCnts(nChan,nBank,:) = Cnts(nChan,nBank,:) - CalTemp * CalFlag
       ENDDO
    ENDDO

    maxPt = SIZE (CalFlag) - 1
    win = zwin * MIFsPerMAF
    kbgn = 0
    iBound = 1
    wend = -1
    mif0 = 0
    last_fit = .FALSE.
    DO

       DO
          IF (CalMIF(mif0) > mlimb .AND. mif0 < maxPt) THEN
             mif0 = mif0 + 1
          ELSE
             EXIT
          ENDIF
       ENDDO

       DO
          IF (mif0 > wend) THEN
             wbgn = wend + 1
             wend = BoundsX(iBound)
             iBound = iBound + 1
          ELSE
             EXIT
          ENDIF
       ENDDO
       ibgn = MAX ((mif0 - win), wbgn)
       iend = MIN ((mif0 + mlimb + win), wend)
       cFlag => CalFlag(ibgn:iend)
       ntotx = COUNT (cFlag == 1)

       IF (ntotx > 1) THEN

          DEALLOCATE (yval, stat=status)
          ALLOCATE (yval(THzChans,THzNum,ibgn:iend))
          DEALLOCATE (tsq, stat=status)
          ALLOCATE (tsq(ibgn:iend))
          DEALLOCATE (tval, stat=status)
          ALLOCATE (tval(ibgn:iend))

          DO i = ibgn, iend
             IF (CalFlag(i) == 1) THEN
                mbgn = i
                EXIT
             ENDIF
          ENDDO
          DO i = iend, ibgn, -1
             IF (CalFlag(i) == 1) THEN
                mend = i
                EXIT
             ENDIF
          ENDDO
          nordr = ((mend - mbgn) / REAL(MIFsPerMAF) + 0.2)
          tmid = SUM (MIFx(ibgn:iend)*cFlag) / REAL(ntotx,r8)
          tval = (MIFx(ibgn:iend) - tmid) * cFlag
          tsq = tval * tval
          mfit0 = 1.0 / ntotx
          DO nBank = 1, THzNum
             DO nChan = 1, THzChans
                yval(nChan,nBank,:) = dCnts(nChan,nBank,ibgn:iend) * cFlag
                av(nChan,nBank) = SUM (yval(nChan,nBank,:)) * mfit0
             ENDDO
          ENDDO
          bv = 0.0
          cv = 0.0
          mfit1 = 0.0; mfit2 = 0.0; mfit3 = 0.0; mfit4 = 0.0
          IF (nordr > 0) THEN
             mfit2 = SUM (tsq)
             IF (nordr == 1) mfit2 = 1.0 / mfit2
             DO nBank = 1, THzNum
                DO nChan = 1, THzChans
                   yfit1(nChan,nBank) = &
                        DOT_PRODUCT (yval(nChan,nBank,:), tval)
                   IF (nordr == 1) bv(nChan,nBank) = yfit1(nChan,nBank) * mfit2
                ENDDO
             ENDDO
          ENDIF

          IF (nordr > 1) THEN
             tsqavg = mfit2 / ntotx
             tsq = (tsq - tsqavg) * cFlag
             mfit3 = -SUM (tsq * tval)
             mfit4 = mfit2
             mfit2 = SUM (tsq * tsq)
             det = mfit2 * mfit4 - mfit3 * mfit3
             DO nBank = 1, THzNum
                DO nChan = 1, THzChans
                   yfit2(nChan,nBank) = &
                        DOT_PRODUCT (yval(nChan,nBank,:), tsq)
                ENDDO
             ENDDO
             mfit2 = mfit2 /det
             mfit3 = mfit3 /det
             mfit4 = mfit4 /det
             bv = mfit2 * yfit1 + mfit3 * yfit2
             cv = mfit3 * yfit1 + mfit4 * yfit2
             av = av - cv * tsqavg
             tsqavg = -2.0 * tsqavg
             mfit1 = tsqavg * mfit3
             mfit2 = mfit2 + tsqavg * mfit4
             mfit3 = 2.0 * mfit3
             last_fit = .TRUE.
          ENDIF

       ENDIF

       DO
          IF (CalMIF(mif0) <= mlimb .AND. mif0 < wend) THEN
             mif0 = mif0 + 1
          ELSE
             EXIT
          ENDIF
       ENDDO

       DO
          IF (CalMIF(mif0) > mlimb .AND. mif0 < wend) THEN
             mif0 = mif0 + 1
          ELSE
             EXIT
          ENDIF
       ENDDO

       IF (last_fit) THEN
          DO i = kbgn, mif0
             t = i - tmid
             Cnts(:,:,i) = Cnts(:,:,i) - av - t * (bv + cv * t)
             val = mfit0 + t * (mfit1 + t * (mfit2 + t * (mfit3 + t * mfit4)))
             VarCnts(:,:,i) = VarCnts(:,:,i) + val + VarGain * Cnts(:,:,i) * &
                  Cnts(:,:,i)
          ENDDO
          kbgn = mif0 + 1
       ENDIF

       IF (mif0 == wend) mif0 = mif0 + 1
       IF (mif0 > maxPt) EXIT

    ENDDO

  END SUBROUTINE THzCal

!=============================================================================
  SUBROUTINE THzStat
!=============================================================================

    INTEGER :: i, nBank, nChan, ntot
    REAL(r8) :: aerr(THzChans,THzNum), Chisq(THzChans,THzNum)
    REAL(r8) :: xavg(THzChans,THzNum), yavg(THzChans,THzNum)
    REAL(r8) :: xval(THzChans,THzNum), yval(THzChans,THzNum)
    REAL(r8) :: vnorm
    REAL :: filter(THzChans)
    REAL, PARAMETER :: ChanBw(25) = (/ 110.0, 110.0, 110.0, 74.0, 74.0, 74.0, &
         55.0, 37.0, 28.0, 18.4, 14.0, 9.2, 7.2, 9.2, 14.0, 18.4, 28.0, 37.0, &
         55.0, 74.0, 74.0, 74.0, 110.0, 110.0, 110.0 /)

    ntot = COUNT (CalFlag == 1)
    IF (ntot == 0) RETURN

    filter = 1000.0 * SQRT (ChanBw / 6.0)
    vnorm = 1.0d0 / ntot
    xavg = 0.0
    yavg = 0.0

    DO i = 1, SIZE (CalFlag)
       IF (CalFlag(i) == 1) THEN
          xval = vnorm / VarCnts(:,:,i)
          xavg = xavg + xval
          yval = Cnts(:,:,i) - CalTemp(i)
          yavg = yavg + xval * yval * yval
       ENDIF
    ENDDO

    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          VarCnts(nChan,nBank,:) = SQRT (yavg(nChan,nBank) * &
               VarCnts(nChan,nBank,:))
       ENDDO
    ENDDO

    WHERE (xavg > 0.0)
       aerr = SQRT (yavg / xavg)
    END WHERE

    Chisq = 0.0
    DO nBank = 1, THzNum
       aerr(:,nBank) = aerr(:,nBank) * filter
       WHERE (ytsys(:,nBank) > 0.0)
          Chisq(:,nBank) = aerr (:,nBank) / ytsys(:,nBank)
       END WHERE
       Chisq(:,nBank) = Chisq(:,nBank) * Chisq(:,nBank)
    ENDDO

  END SUBROUTINE THzStat

!=============================================================================
  SUBROUTINE CalibrateTHz (more_data)
!=============================================================================

    LOGICAL :: more_data

PRINT *, 'start calibrating...'

    CALL BuildCalVectors

    CALL THzDel

    CALL THzCal

    CALL THzStat

    more_data = Cal_start <= CalBuf%MAFs

PRINT *, 'end calibrating...'

  END SUBROUTINE CalibrateTHz

!=============================================================================
END MODULE THzCalibration
!=============================================================================

! $Log$
! Revision 2.2  2003/02/05 21:31:55  perun
! Calculate variance and chi square
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
