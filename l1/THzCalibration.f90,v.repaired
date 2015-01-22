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
MODULE THzCalibration ! Calibration data and routines for the THz module
!=============================================================================

  USE MLSCommon, ONLY: r8, rm
  USE MLSL1Common, ONLY: MaxMIFs, THzChans, THzNum, LO1, L1BFileInfo
  USE L0_sci_tbls, ONLY: THz_Sci_pkt_T
  USE EngTbls, ONLY : Eng_MAF_T
  USE MLSL1Config, ONLY: MIFsTHz
  USE MLSFillValues, ONLY: isNaN
  ! USE dump_0, ONLY: dump
  use output_m, only: switchOutput, revertOutput
  use PrintIt_m, only: MLSMessageConfig
 
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CalibrateTHz, Chan_type_T, MAFdata_T, CalBuf_T, CalBuf, Cnts, &
       VarCnts, SpaceTemp, nvBounds, ColdCnts, HotCnts, GoodCal

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(THzChans,THzNum)      ! standard THz filter banks
  END TYPE Chan_type_T

  !! THz Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (THz_Sci_pkt_T) :: SciMIF(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     REAL :: CalTgtTemp       ! Average Calibration Target Temperature (C)
     REAL :: LimbCalTgtTemp   ! Average Limb Calibration Target Temperature (C)
     INTEGER :: BO_stat(MIFsTHz) = 0      ! Bright Objects status (THz FOV)
     REAL :: MIFprecSign(0:(MaxMIFs-1))   ! Radiance precision sign per MIF
     INTEGER :: last_MIF
     INTEGER :: BandSwitch(5) ! band switch positions
  END TYPE MAFdata_T

  TYPE CalBuf_T
     INTEGER :: MAFs, Cal_start, Cal_end
     LOGICAL :: BankGood(THzNum)
     TYPE (MAFdata_T), DIMENSION(:), ALLOCATABLE :: MAFdata
  END TYPE CalBuf_T

  TYPE (CalBuf_T), TARGET, SAVE :: CalBuf

  INTEGER :: Cal_start = 1
  INTEGER :: Cal_end, ntot
  INTEGER, DIMENSION(:), ALLOCATABLE :: CalMIF, MIFx, xMIF0
  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: CalFlag
  INTEGER, DIMENSION(:), ALLOCATABLE :: CalInfo   ! Similar to Herb's IDL
  INTEGER, DIMENSION(:), ALLOCATABLE :: BoundsX
  LOGICAL, DIMENSION(:), ALLOCATABLE :: Bounds
  INTEGER, DIMENSION(:), ALLOCATABLE :: nvBounds
  REAL(r8), DIMENSION(:), ALLOCATABLE :: Bias, CalTemp
  REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Cnts, dCnts, VarCnts
  REAL(r8), DIMENSION(:), ALLOCATABLE :: VarGain
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: ColdCnts, HotCnts
  LOGICAL, DIMENSION(:), ALLOCATABLE :: BiasGood
  REAL :: MaxBias, SpaceTemp
  REAL(r8) :: aerr(THzChans,THzNum), Chisq(THzChans,THzNum)
  REAL(r8) :: dLlo(THzChans,THzNum), yTsys(THzChans,THzNum)
  LOGICAL :: GoodCal = .TRUE.
  CHARACTER(len=80) :: msg

! Default gains:

  REAL(r8), PARAMETER :: gainZ(THzChans,THzNum) = RESHAPE ( (/ &
       1.080, 1.407, 1.627, 1.489, 1.214, 1.248, 1.388, 1.353, &
       1.262, 1.166, 1.128, 1.086, 1.073, 1.051, 1.021, 0.966, &
       0.938, 0.960, 1.183, 1.350, 1.208, 1.279, 1.510, 1.612, 1.514, &
       1.122, 1.221, 1.279, 1.227, 1.060, 1.136, 1.167, 1.117, &
       1.113, 1.195, 1.190, 1.175, 1.159, 1.147, 1.126, 1.179, &
       1.128, 1.041, 1.087, 1.145, 1.216, 1.351, 1.480, 1.741, 1.805, &
       0.970, 1.208, 1.274, 1.296, 1.037, 1.091, 0.990, 0.907, &
       1.070, 1.121, 1.037, 0.921, 0.922, 0.892, 0.838, 0.807, &
       0.805, 0.889, 1.339, 1.258, 1.307, 1.402, 1.159, 1.097, 1.230, &
       1.729, 1.842, 1.584, 1.638, 1.590, 1.271, 1.492, 1.569, &
       1.619, 1.533, 1.513, 1.426, 1.375, 1.380, 1.334, 1.252, &
       1.334, 1.458, 1.764, 1.709, 1.741, 1.801, 1.652, 1.814, 2.007, &
       0.977, 1.139, 1.189, 1.535, 1.275, 1.511, 1.547, 1.328, &
       1.285, 1.286, 1.380, 1.389, 1.505, 1.559, 1.658, 1.776, &
       1.715, 1.532, 1.769, 1.677, 1.749, 1.834, 1.786, 1.571, 2.206, &
       2.234, 2.230, 1.873, 1.981, 1.730, 1.448, 1.492, 1.703, &
       2.007, 1.999, 1.847, 1.658, 1.575, 1.568, 1.500, 1.431, &
       1.459, 1.586, 2.085, 1.864, 1.888, 1.883, 1.890, 1.879, 2.098 /), &
       (/ THzChans, THzNum /) )

CONTAINS

!=============================================================================
  SUBROUTINE BuildCalVectors
!=============================================================================

!    USE THzUtils, ONLY: Bias_err
    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Rad, ONLY: UpdateRadSignals, RadPwr

    INTEGER :: i, last_MIF, mindx, MIF0, MIFno
    INTEGER :: nBank, numMIFs, nMIFsM1, status, Switch5
    TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()
    INTEGER, PARAMETER :: Cal_size = 240   ! Nominal orbit
    REAL, PARAMETER :: LO1R5 = LO1(5)  ! Radiometer 5 1st LO frequency
    logical :: CalTgtTempIsNAN
    logical :: LimbCalTgtTempIsNAN
    logical :: FB15NotAllGood
    logical :: FB20NotAllGood

    SpaceTemp = L1Config%Calib%THzSpaceTemp
    MaxBias = L1Config%Calib%THzMaxBias

! Figure out how many MIFs needed for calibration:

    Cal_start = 1
    Cal_end = CalBuf%MAFs
    numMIFs = 0
    DO i = Cal_start, Cal_end
       numMIFs = numMIFs + CalBuf%MAFdata(i)%last_MIF + 1
    ENDDO
    nMIFsM1 = numMIFs - 1

! Allocate vectors:

    DEALLOCATE (CalMIF, stat=status)
    ALLOCATE (CalMIF(0:nMIFsM1))

    DEALLOCATE (MIFx, stat=status)
    ALLOCATE (MIFx(0:nMIFsM1))
    DO i = 0, nMIFsM1
       MIFx(i) = i
    ENDDO
    DEALLOCATE (xMIF0, stat=status)
    ALLOCATE (xMIF0(0:CalBuf%MAFs-1))

    DEALLOCATE (CalFlag, stat=status)  ! calibration = 1, otherwise = 0
    ALLOCATE (CalFlag(0:nMIFsM1))
    CalFlag = 0    ! Indicate not a calibration MIF (yet)
    DEALLOCATE (CalInfo, stat=status)  ! 1 = cold, 2 = hot, otherwise = 0
    ALLOCATE (CalInfo(0:nMIFsM1))
    CalInfo = 0

    DEALLOCATE (Bias, stat=status)
    ALLOCATE (Bias(0:nMIFsM1))

    DEALLOCATE (BiasGood, stat=status)
    ALLOCATE (BiasGood(0:nMIFsM1))

    DEALLOCATE (CalTemp, stat=status)
    ALLOCATE (CalTemp(0:nMIFsM1))
    CalTemp = 0.0

    DEALLOCATE (Cnts, stat=status)
    ALLOCATE (Cnts(THzChans,THzNum,0:nMIFsM1))
    DEALLOCATE (dCnts, stat=status)
    ALLOCATE (dCnts(THzChans,THzNum,0:nMIFsM1))
    DEALLOCATE (VarCnts, stat=status)
    ALLOCATE (VarCnts(THzChans,THzNum,0:nMIFsM1))
    DEALLOCATE (VarGain, stat=status)
    ALLOCATE (VarGain(0:nMIFsM1))
    Cnts = 0.0
    VarCnts = -1.0
    VarGain = 1.0

    if ( any(isNaN(Cnts)) ) then
      print *, 'NaNs in Cnts at start of BuildCalVectors'
      print *, 'SpaceTemp ', SpaceTemp
      print *, 'Cal_start, Cal_end ', Cal_start, Cal_end
    endif
! Fill the vectors:

    CalBuf%BankGood = .TRUE.
    Switch5 = CalBuf%MAFdata(Cal_start)%SciMIF(0)%BandSwitch(5) ! init compare

    MIF0 = 0   ! index for MIFs
    CalTgtTempIsNAN = .false.
    LimbCalTgtTempIsNAN = .false.
    FB15NotAllGood = .false.
    FB20NotAllGood = .false.
    DO i = Cal_start, Cal_end

       xMIF0(i-1) = MIF0                 ! Save MIF0 indexes
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
             CalInfo(mindx) = 1
          ELSE IF (CurMAFdata%SciMIF(MIFno)%SwMirPos == "T" .AND. &
               Bias(mindx) < MaxBias) THEN
             CalTemp(mindx) = CurMAFdata%CalTgtTemp
             IF (SpaceTemp < 0.0) CalTemp(mindx) = CalTemp(mindx) - &
                  CurMAFdata%LimbCalTgtTemp
             IF (CalTemp(mindx) > 200.0) CalTemp(mindx) = &
                  RadPwr (LO1R5, REAL(CalTemp(mindx))) !for in orbit space = 0
             CalFlag(mindx) = 1
             CalInfo(mindx) = 2
             CalTgtTempIsNAN = CalTgtTempIsNAN .or. isNaN(CurMAFdata%CalTgtTemp)
             LimbCalTgtTempIsNAN = LimbCalTgtTempIsNAN .or. isNaN(CurMAFdata%LimbCalTgtTemp)
          ENDIF
       ENDDO

! Check if Band Switch 5 changes:

       CalBuf%BankGood(1) = .TRUE.
       IF (ANY (CurMAFdata%SciMIF(0:last_MIF)%BandSwitch(5) /= Switch5)) then
          CalBuf%BankGood(1) = .FALSE.  ! Filter Bank 15 not all good
          FB15NotAllGood = .true.
       endif

! Check if Band 20 is connected:

       CalBuf%BankGood(6) = .TRUE.
       IF (ANY (CurMAFdata%SciMIF(0:last_MIF)%BandSwitch(4) /= 20)) THEN
          CalBuf%BankGood(6) = .FALSE. ! Filter Bank 12 (Band 20) not all good
            FB20NotAllGood = .true.
          DO MIFno = 0, last_MIF
             CurMAFdata%SciMIF(MIFno)%FB(:,6) = 0.0   ! Clear counts
          ENDDO
       ENDIF

       MIF0 = MIF0 + last_MIF + 1

    ENDDO
    if ( FB15NotAllGood ) print *, 'Filter Bank 15 not all good'
    if ( FB20NotAllGood ) print *, 'Filter Bank 20 not all good'
    if ( CalTgtTempIsNAN ) print *, 'CalTgtTemp Is NaN'
    if ( LimbCalTgtTempIsNAN ) print *, 'LimbCalTgtTemp Is NaN'

! Adjust LLO Bias V if 2 or more consecutive errs:

!!$    DO i = (numMIFs - 2), 1, -1
!!$       IF (Bias(i) == Bias_err .AND. Bias(i-1) == Bias_err) THEN
!!$          Bias(i+1) = Bias_err
!!$          CalFlag(i+1) = 0
!!$       ENDIF
!!$    ENDDO

! Determine which filter bank counts are good to cal:

    DO nBank = 1, THzNum
       IF (ALL (isNaN(Cnts(:,nBank,:) ))) &
         & print *, 'Filter Bank ', nBank, ' are all NaN'
    enddo

    DO nBank = 1, THzNum
       IF (ALL (Cnts(:,nBank,:) == 0.0)) THEN
          CalBuf%BankGood(nBank) = .FALSE.
         print *, 'Filter Bank ', nBank, ' not all good'
       ELSE IF (.NOT. CalBuf%BankGood(nBank)) THEN
          Cnts(:,nBank,:) = 0.0  ! clear bad data
       ENDIF
    ENDDO

! Set signals based on end cal index

    CALL UpdateRadSignals (CalBuf%MAFdata(Cal_end)%BandSwitch)

! Save current cal start & end:

    CalBuf%Cal_start = Cal_start
    CalBuf%Cal_end = Cal_end

    if ( any(isNaN(CalTemp)) ) print *, 'NaNs in CalTemp at start of BuildCalVectors'
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at end of BuildCalVectors'
  END SUBROUTINE BuildCalVectors

!=============================================================================
  SUBROUTINE RestoreCnts  ! Restore Cnts from CalBuf FB data
!=============================================================================

    INTEGER :: i, MIF0, MIFno, last_MIF, mindx
    TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()

    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at start of RestoreCnts'
    MIF0 = 0   ! index for MIFs
    DO i = 1, CalBuf%MAFs
       CurMAFdata => CalBuf%MAFdata(i)
       last_MIF = CurMAFdata%last_MIF
        DO MIFno = 0, last_MIF
          mindx = MIF0 + MIFno
          Cnts(:,:,mindx) = &
               CurMAFdata%SciMIF(MIFno)%FB(:,:)
       ENDDO
       MIF0 = MIF0 + last_MIF + 1
    ENDDO
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at end of RestoreCnts'

  END SUBROUTINE RestoreCnts

!=============================================================================
  SUBROUTINE THzBound (ibgn, iend, iorbit, fillnvbounds)
!=============================================================================

    INTEGER, INTENT (in) :: ibgn, iend, iorbit
    LOGICAL, INTENT (in), OPTIONAL :: fillnvbounds

    INTEGER :: i, j, nBounds, status
    INTEGER :: ibound, iMAF, mbgn, mend, nv

    DEALLOCATE (Bounds, stat=status)
    ALLOCATE (Bounds(ibgn:iend))

! Determine bounds for fitting:

    Bounds = .FALSE.   ! nothing yet
    WHERE (Bias(ibgn:iend) >= MaxBias)  ! Re-optimizing
       Bounds = .TRUE.
    END WHERE
    WHERE (CalMIF(ibgn:iend) < 0)       ! Missing MIFs
       Bounds = .TRUE.
    END WHERE
    BiasGood(ibgn:iend) = .NOT. Bounds  ! Good Bias values
    Bounds(iend) = .TRUE.               ! last always a bounds
    nBounds = COUNT (Bounds)
    DEALLOCATE (BoundsX, stat=status)
    ALLOCATE (BoundsX(nBounds))
    j = 1
    DO i = ibgn, iend
       IF (Bounds(i)) THEN
          BoundsX(j) = i-ibgn
          j = j + 1
       ENDIF
    ENDDO

    IF (.NOT. PRESENT (fillnvbounds)) THEN

! Write to Log file:

       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (msg, &
            '("Cal Number: ", i0, ", nBounds: ", i0, ", nBiasGood: ", i0)') &
            iorbit, nBounds, COUNT (BiasGood(ibgn:iend))
       WRITE (L1BFileInfo%LogId, *) TRIM (msg)

       RETURN     ! No nvbounds yet...

    ENDIF

! Determine nvBounds:

    DEALLOCATE (nvBounds, stat=status)
    ALLOCATE (nvBounds(CalBuf%MAFs))
    nvBounds = 0
    iMAF = 1
    mbgn = 0
    DO iBound = 1, nBounds
       mend = BoundsX(iBound)
       nv = mend - mbgn + 1
       DO i = mbgn, mend
          IF (CalMIF(i) == 0) THEN
             nvBounds(iMAF) = nv
             iMAF = iMAF + 1
          ENDIF
       ENDDO
       mbgn = mend + 1
    ENDDO

  END SUBROUTINE THzBound

!=============================================================================
  SUBROUTINE THzDel (Bias, BiasGood, CalFlag, CalTemp, xMIF0, Cnts, VarCnts, &
       VarGain)
!=============================================================================

    USE Constants, ONLY: Pi
    USE MatrixModule_0, ONLY: MatrixInversion

    INTEGER, TARGET :: CalFlag(0:)
    INTEGER :: xMIF0(0:)
    REAL(r8) :: Bias(0:), VarGain(0:)
    REAL(r8) :: CalTemp(0:)
    REAL(r8) :: Cnts(:,:,0:), VarCnts(:,:,0:)
    LOGICAL :: BiasGood(0:)

    INTEGER, PARAMETER :: nfit = 4, nfitm1 = nfit - 1
    INTEGER :: fitno = nfit
    INTEGER :: i, j, iBound, ibgn, iend, nBounds, ntotx, status
    INTEGER :: nChan, nBank, nMAFs, nMIFs
    INTEGER, DIMENSION(:), POINTER :: cFlag
    REAL, DIMENSION(:), POINTER :: ang, cv, sv
    REAL(r8) :: trms, xamp
    REAL(r8), DIMENSION(:), POINTER :: biasx, dbias, tempx, ct, st, val
    REAL(r8) :: AvgBias, AvgTemp, avgx
    REAL(r8), DIMENSION(:), TARGET, ALLOCATABLE :: BiasBuf, TempBuf, &
         cv_buf, sv_buf
    REAL(r8) :: tt(0:nfitm1,0:nfitm1), tt0(0:nfitm1,0:nfitm1)
    REAL(rm) :: rtt(nfit,nfit), rtt2(2,2)
    REAL(r8) :: yv(THzChans,THzNum,0:nfitm1), tv(THzChans,THzNum,0:nfitm1)
    REAL(r8) :: dCal(THzChans,THzNum), amp(THzChans,THzNum)

    REAL, PARAMETER :: tmin = 10.0

    call switchOutput( 'stdout' )
    ! MLSMessageConfig%adHocPrintToStdout = .true.
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at start of THzDel'
    nMAFs = SIZE (xMIF0)
    nMIFs = SIZE (Bias)
!    DEALLOCATE (dbias, stat=status)
    ALLOCATE (dbias(0:nMIFs-1))
    dbias = bias
!    DEALLOCATE (val, stat=status)
    ALLOCATE (val(0:nMIFs-1))
    tv = 0.0d0
    yv = 0.0d0
    amp = 0.d0
    IF (nfit > 2) THEN
       nMAFs = nMAFs - 1
!       DEALLOCATE (ang, stat=status)
       ALLOCATE (ang(0:nMIFs-1))
!       DEALLOCATE (cv, stat=status)
       ALLOCATE (cv(0:nMIFs-1))
!       DEALLOCATE (sv, stat=status)
       ALLOCATE (sv(0:nMIFs-1))
       DO i = 0, (nMIFs - 1)
          ang(i) = (2.0 * Pi / nMIFs) * i
          cv(i) = COS (ang(i))
          sv(i) = SIN (ang(i))
       ENDDO
    ENDIF

    ntot = COUNT (CalFlag == 1)
    IF (ntot <= 1) RETURN

    AvgBias = SUM (Bias * CalFlag) / ntot
    AvgTemp = SUM (CalTemp * CalFlag) / ntot
    ! print *, 'AvgBias ', AvgBias
    ! print *, 'AvgTemp ', AvgTemp
    ! call dump ( CalTemp, 'CalTemp' )
!    DEALLOCATE (BiasBuf, stat=status)
    ALLOCATE (BiasBuf(0:nMIFs-1))
    BiasBuf = Bias * CalFlag
!    DEALLOCATE (TempBuf, stat=status)
    ALLOCATE (TempBuf(0:nMIFs-1))
    TempBuf = CalTemp * CalFlag
!    DEALLOCATE (cv_buf, stat=status)
    ALLOCATE (cv_buf(0:nMIFs-1))
    cv_buf = cv * CalFlag
!    DEALLOCATE (sv_buf, stat=status)
    ALLOCATE (sv_buf(0:nMIFs-1))
    sv_buf = sv * CalFlag

    nBounds = SIZE (BoundsX)
    tt = 0.0d0
    tt(0,0) = ntot / (819.2 * 819.2)
    DO i = 1, (nfit - 1)
       tt(i,i) = 0.5
    ENDDO
    ibgn = 0
    DO iBound = 1, nBounds
       iend = BoundsX(iBound)
       cFlag => CalFlag(ibgn:iend)
       ntotx = COUNT (cFlag == 1)
       IF (ntotx > 1) THEN
          ! call dump( TempBuf(ibgn:iend), 'TempBuf' )
          ! call dump( BiasBuf(ibgn:iend), 'BiasBuf' )
          ! call dump( CalFlag(ibgn:iend), 'CalFlag' )
          biasx => BiasBuf(ibgn:iend)
          avgx = SUM (biasx) / ntotx
          biasx = (biasx - avgx) * cFlag
          tempx => TempBuf(ibgn:iend)
          avgx = SUM (tempx) / ntotx
          tempx = (tempx - avgx) * cFlag
          tt(0,0) = tt(0,0) + SUM (biasx * biasx)
          tt(1,1) = tt(1,1) + SUM (tempx * tempx)
          tt(1,0) = tt(1,0) + SUM (tempx * biasx)
          DO nBank = 1, THzNum
             yv(:,nBank,0) = yv(:,nBank,0) + &
                  MATMUL (Cnts(:,nBank,ibgn:iend), biasx)
             yv(:,nBank,1) = yv(:,nBank,1) + &
                  MATMUL (Cnts(:,nBank,ibgn:iend), tempx)
          ENDDO
          IF (nfit > 2) THEN
             ct => cv_buf(ibgn:iend)
             ct = ct * tempx
             st => sv_buf(ibgn:iend)
             st = st * tempx
             tt(2,2) = tt(2,2) + SUM (ct * ct)
             tt(2,0) = tt(2,0) + SUM (ct * biasx)
             tt(2,1) = tt(2,1) + SUM (ct * tempx)
             tt(3,3) = tt(3,3) + SUM (st * st)
             tt(3,0) = tt(3,0) + SUM (st * biasx)
             tt(3,1) = tt(3,1) + SUM (st * tempx)
             tt(3,2) = tt(3,2) + SUM (st * ct)
             DO nBank = 1, THzNum
                yv(:,nBank,2) = yv(:,nBank,2) + &
                     MATMUL (Cnts(:,nBank,ibgn:iend), ct)
                yv(:,nBank,3) = yv(:,nBank,3) + &
                     MATMUL (Cnts(:,nBank,ibgn:iend), st)
             ENDDO
          ENDIF
       ENDIF
       ibgn = iend + 1
    ENDDO
    DO i = 0, (nfit - 2)
       DO j = i + 1, nfit -1
          tt(i,j) = tt(j,i)
       ENDDO
    ENDDO
    trms = 2.0 * SQRT(tt(1,1) / ntot)

    tt0 = tt
    rtt = tt
    ! call dump( rtt, 'rtt' )
    CALL MatrixInversion (rtt)
    tt = rtt

    ! call dump( yv, 'yv' )
    ! call dump( tt, 'tt' )
    DO nBank = 1, THzNum
       tv(:,nBank,:) = MATMUL (yv(:,nBank,:), tt)
    ENDDO
    ! call dump( tv, 'tv (1st time)' )
    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          if ( tv(nChan,nBank,1) /= 0.d0 ) &
            & amp(nChan,nBank) = SQRT (tv(nChan,nBank,2) * tv(nChan,nBank,2) + &
              tv(nChan,nBank,3) * tv(nChan,nBank,3)) / ABS (tv(nChan,nBank,1))
       ENDDO
    ENDDO
    xamp = MAXVAL(amp)
    IF (xamp > 0.2) THEN
       PRINT *, 'Switch to const gain: ', xamp
       fitno = 2
       tt(0,0) = tt0(0,0); tt(1,1) = tt0(1,1)
       tt(0,1) = tt0(0,1); tt(1,0) = tt0(1,0)
       rtt2 = tt(0:1,0:1)
       CALL MatrixInversion (rtt2)
       tt(0:1,0:1) = rtt2
       DO nBank = 1, THzNum
          tv(:,nBank,0:1) = MATMUL (yv(:,nBank,0:1), tt(0:1,0:1))
       ENDDO
    ENDIF
    ! call dump( tv, 'tv (2nd time)' )

    dLlo = tv(:,:,0)
    dCal = tv(:,:,1)
    vargain = tt(1,1)

    IF (trms < tmin) THEN
       vargain = 0.0
       fitno = 2
       dCal = gainZ
    ENDIF

    WHERE (BiasGood)
       Bias = Bias - AvgBias  ! Adjust bias for remainder of this routine
    ENDWHERE
    ! call dump( Bias, 'Bias' )
    ! call dump( dLlo, 'dLlo' )
    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          WHERE (BiasGood)
             Cnts(nChan,nBank,:) = Cnts(nChan,nBank,:) - dLlo(nChan,nBank) * &
                  Bias
             VarCnts(nChan,nBank,:) = Bias * Bias * tt(0,0) + 1.0
          ENDWHERE
       ENDDO
    ENDDO
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts after 1st bias adjustment'

    IF (fitno == 2) THEN
       DO nBank = 1, THzNum
          DO nChan = 1, THzChans
             Cnts(nChan,nBank,:) = Cnts(nChan,nBank,:) / &
                  MAX (ABS (dCal(nChan,nBank)), 0.01d0)
          ENDDO
       ENDDO
    ELSE
       WHERE (BiasGood)
          vargain = tt(1,1) + 2.0 * tt(1,2) * cv + 2.0 * tt(1,3) * sv + &
               cv * (tt(2,2) * cv + 2.0 * tt(2,3) * sv) + sv * sv * tt(3,3)
       ENDWHERE

       DO nBank = 1, THzNum
          DO nChan = 1, THzChans
             val = tv(nChan,nBank,1) + tv(nChan,nBank,2) * cv + &
                  tv(nChan,nBank,3) * sv
             DO i = 0, (nMIFs - 1)
                Cnts(nChan,nBank,i) = Cnts(nChan,nBank,i) / &
                     MAX (ABS (val(i)), 0.01d0)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts after 2nd bias adjustment'

    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          yTsys(nChan,nBank) = SUM (Cnts(nChan,nBank,:) * CalFlag) / ntot - &
               AvgTemp
       ENDDO
    ENDDO

! restore bias:

    bias = dbias

! Deallocate the temporary buffers
!    REAL, DIMENSION(:), POINTER :: ang, cv, sv
!    REAL(r8), DIMENSION(:), TARGET, ALLOCATABLE :: BiasBuf, TempBuf, &
!         cv_buf, sv_buf

    DEALLOCATE (ang, stat=status)
    DEALLOCATE (cv,  stat=status)
    DEALLOCATE (sv,  stat=status)

    DEALLOCATE (BiasBuf, stat=status)
    DEALLOCATE (TempBuf, stat=status)
    DEALLOCATE (cv_buf,  stat=status)
    DEALLOCATE (sv_buf,  stat=status)

    DEALLOCATE (dbias, stat=status)
    DEALLOCATE (val, stat=status)
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at end of THzDel'
    ! MLSMessageConfig%adHocPrintToStdout = .false.
    call revertOutput

  END SUBROUTINE THzDel

!=============================================================================
  SUBROUTINE THzCal (ColdCal)
!=============================================================================

    LOGICAL, INTENT(IN), OPTIONAL :: ColdCal

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
    LOGICAL :: first_fit, last_fit, IsColdCal

    PRINT *, 'THzCal'

    IF (PRESENT (ColdCal)) THEN
       IsColdCal = ColdCal
    ELSE
       IsColdCal = .FALSE.
    ENDIF

    IF (COUNT (CalFlag == 1) == 0) RETURN  ! Nothing to calibrate

    IF (ISColdCal) THEN
       WHERE (CalTemp /= 0.0 .AND. CalFlag == 1)
          CalFlag = 0
       ENDWHERE
    ENDIF
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
    first_fit = .TRUE.
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
          ELSE
             EXIT
          ENDIF
       ENDDO
       ibgn = MAX ((mif0 - win), wbgn)
       iend = MIN ((mif0 + mlimb + win), wend)
       cFlag => CalFlag(ibgn:iend)
       ntotx = COUNT (cFlag == 1)

       ! IF (ntotx > 1 .AND. .NOT. first_fit) THEN
       IF (ntotx > 1) THEN

!          DEALLOCATE (yval, stat=status)
          ALLOCATE (yval(THzChans,THzNum,ibgn:iend))
!          DEALLOCATE (tsq, stat=status)
          ALLOCATE (tsq(ibgn:iend))
!          DEALLOCATE (tval, stat=status)
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
             IF (ABS (det) > 0.0) THEN
                mfit2 = mfit2 /det
                mfit3 = mfit3 /det
                mfit4 = mfit4 /det
             ENDIF
             bv = mfit2 * yfit1 + mfit3 * yfit2
             cv = mfit3 * yfit1 + mfit4 * yfit2
             av = av - cv * tsqavg
             tsqavg = -2.0 * tsqavg
             mfit1 = tsqavg * mfit3
             mfit2 = mfit2 + tsqavg * mfit4
             mfit3 = 2.0 * mfit3
          ENDIF
          last_fit = .TRUE.
!
          DEALLOCATE (yval,  stat=status)
          DEALLOCATE (tsq,   stat=status)
          DEALLOCATE (tval,  stat=status)
       ENDIF
       first_fit = .FALSE.

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
             IF (BiasGood(i)) THEN
                t = i - tmid
                Cnts(:,:,i) = Cnts(:,:,i) - av - t * (bv + cv * t)
                val = mfit0 + t * (mfit1 + t * (mfit2 + t * (mfit3 + t * &
                     mfit4)))
                VarCnts(:,:,i) = VarCnts(:,:,i) + val + VarGain(i) * &
                     Cnts(:,:,i) * Cnts(:,:,i)
             ELSE
                Cnts(:,:,i) = Cnts(:,:,i) - av
                VarCnts(:,:,i) = -1.0
                CalFlag(i) = 0
                CalInfo(i) = 0
             ENDIF
          ENDDO
          kbgn = mif0 + 1
       ENDIF

       IF (mif0 == wend) mif0 = mif0 + 1
       IF (mif0 > maxPt) EXIT

    ENDDO


  END SUBROUTINE THzCal

!=============================================================================
  SUBROUTINE THzReselect (Stype, nbad)
!=============================================================================

    CHARACTER(LEN=*), INTENT (IN) :: Stype
    INTEGER, INTENT (OUT) :: nbad

    INTEGER :: mindx, n, ncold, nhot, lbad
    INTEGER, PARAMETER :: limbMIF = 116
    REAL :: avg, wt2(THzChans,THzNum), norm
    LOGICAL :: ChanGood(THzChans,THzNum), cold_type, hot_type, limb_type
    REAL, PARAMETER :: cmax = 10.0
    REAL, PARAMETER :: hmax = 0.08
    REAL, PARAMETER :: limb_range(2) = (/ -10.0, 200.0 /)
    REAL, PARAMETER :: lavg = 0.5 * (limb_range(1) + limb_range(2))
    REAL, PARAMETER :: lmax = lavg - limb_range(1)

    PRINT *, 'Reselecting...'
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at start of Reselecting'

    limb_type = (INDEX(Stype, "L") /= 0)
    cold_type = (INDEX(Stype, "C") /= 0)
    hot_type = (INDEX(Stype, "H") /= 0)

    norm = 0.0
    DO n = 1, ThzNum
       ChanGood(:,n) = CalBuf%BankGood(n)
       WHERE (aerr(:,n) > 1.0e-30)
          wt2(:,n) = 1.0 / aerr(:,n)
       ELSEWHERE
          wt2(:,n) = 1.0
       END WHERE
       WHERE (.NOT. ChanGood(:,n))
          wt2(:,n) = 0.0
       END WHERE
       wt2(:,n) = wt2(:,n) * wt2(:,n)
       norm = norm + SUM (wt2(:,n) + 1.0e-30)
    ENDDO
    wt2 = wt2 / norm

    nbad = 0 ; ncold = 0; nhot = 0 ; lbad = 0
    DO mindx = 0, (SIZE(CalInfo) - 1)
       IF (cold_type .AND. CalInfo(mindx) == 1) THEN            ! "cold" cal
          avg = SUM (wt2 * cnts(:,:,mindx))
          IF (ChanGood(1,1) .AND. ABS(avg) > cmax .AND. &
               varcnts(1,1,mindx) > 0.0) THEN
             nbad = nbad + 1
             ncold = ncold + 1
             CalInfo(mindx) = 0
             CalFlag(mindx) = 0
          ENDIF
       ELSE IF (hot_type .AND. CalInfo(mindx) == 2) THEN       ! "hot" cal
          avg = SUM (wt2 * cnts(:,:,mindx))
          avg = avg / CalTemp(mindx) - 1.0
          IF (ChanGood(1,1) .AND. ABS(avg) > hmax .AND. &
               varcnts(1,1,mindx) > 0.0) THEN
             nbad = nbad + 1
             nhot = nhot + 1
             CalInfo(mindx) = 0
             CalFlag(mindx) = 0
          ENDIF
       ELSE IF (limb_type .AND. CalMIF(mindx) <= limbMIF) THEN  ! limb views
          avg = SUM (wt2 * cnts(:,:,mindx))
          IF (ChanGood(1,1) .AND. (ABS(avg - lavg) > lmax) .AND. &
               varcnts(1,1,mindx) > 0.0) THEN
             varcnts(:,:,mindx) = -ABS (varcnts(:,:,mindx))
             lbad = lbad + 1
          ENDIF
       ENDIF
    ENDDO

    IF (limb_type) THEN
       PRINT *, 'Limb view rejects: ', lbad
    ELSE
       PRINT *, 'Cals rejected (cold, hot, total): ', ncold, nhot, nbad
    ENDIF

    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at end of Reselecting'
  END SUBROUTINE THzReselect

!=============================================================================
  SUBROUTINE THzStat
!=============================================================================

    INTEGER :: i, nBank, nChan, ntot
    INTEGER :: mbgn, mend, nCold, nHot, status
    REAL(r8) :: xavg(THzChans,THzNum), yavg(THzChans,THzNum)
    REAL(r8) :: xval(THzChans,THzNum), yval(THzChans,THzNum)
    REAL(r8) :: vnorm
    REAL :: filter(THzChans)
    REAL, PARAMETER :: ChanBw(25) = (/ 110.0, 110.0, 110.0, 74.0, 74.0, 74.0, &
         55.0, 37.0, 28.0, 18.4, 14.0, 9.2, 7.2, 9.2, 14.0, 18.4, 28.0, 37.0, &
         55.0, 74.0, 74.0, 74.0, 110.0, 110.0, 110.0 /)
    TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()

    ntot = COUNT (CalFlag == 1)
    IF (ntot == 0) THEN
       goodcal = .FALSE.
       RETURN
    ENDIF

    filter = 1000.0 * SQRT (ChanBw / 6.0)
    vnorm = 1.0d0 / ntot
    xavg = 0.0
    yavg = 0.0

    DO i = 0, (SIZE (CalFlag) - 1)
       IF (CalFlag(i) == 1) THEN
          xval = vnorm / VarCnts(:,:,i)
          xavg = xavg + xval
          yval = Cnts(:,:,i) - CalTemp(i)
          yavg = yavg + xval * yval * yval
       ENDIF
    ENDDO

    DO nBank = 1, THzNum
       DO nChan = 1, THzChans
          WHERE (VarCnts(nChan,nBank,:) >= 0.0)
             VarCnts(nChan,nBank,:) = SQRT (yavg(nChan,nBank) * &
                  VarCnts(nChan,nBank,:))
          ENDWHERE
       ENDDO
    ENDDO

    WHERE (xavg > 0.0)
       aerr = SQRT (yavg / xavg)
    END WHERE

    Chisq = 0.0
    DO nBank = 1, THzNum
       WHERE (ytsys(:,nBank) > 0.0)
          Chisq(:,nBank) = (aerr (:,nBank) * filter) / ytsys(:,nBank)
       END WHERE
       Chisq(:,nBank) = Chisq(:,nBank) * Chisq(:,nBank)
    ENDDO

    DEALLOCATE (ColdCnts, stat=status)
    ALLOCATE (ColdCnts(THzChans,THzNum,CalBuf%MAFs))
    ColdCnts = 0.0
    DEALLOCATE (HotCnts, stat=status)
    ALLOCATE (HotCnts(THzChans,THzNum,CalBuf%MAFs))
    HotCnts = 0.0

! Determine average cold and hot counts for each MAF:

    DO i = 1, CalBuf%MAFs
       CurMAFdata => CalBuf%MAFdata(i)
       mbgn = xMIF0(i-1)
       mend = mbgn + CurMAFdata%last_MIF
       nCold = COUNT (CalInfo(mbgn:mend) == 1)
       nHot = COUNT (CalInfo(mbgn:mend) == 2)

       DO nBank = 1, THzNum
          DO nChan = 1, THzChans
             IF (nCold > 0) THEN
                ColdCnts(nChan,nBank,i) = SUM (Cnts(nChan,nBank,mbgn:mend), &
                     (CalInfo(mbgn:mend) == 1)) / nCold
             ENDIF
             IF (nHot > 0) THEN
                HotCnts(nChan,nBank,i) = SUM (Cnts(nChan,nBank,mbgn:mend), &
                     (CalInfo(mbgn:mend) == 2)) / nHot
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE THzStat

!=============================================================================
  SUBROUTINE THzCalDay
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo
    USE OutputL1B, ONLY: OutputL1B_DiagsT

    INTEGER :: i, ibgn, iend, iendm, iorbit, iorg, iref, norbits, status
    INTEGER :: nMAFs, nMIFs, maf1, maf2
    INTEGER, SAVE :: diagno = 1
    REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tCnts, tVarCnts
    REAL(r8), DIMENSION(:), ALLOCATABLE :: tVarGain
    INTEGER, PARAMETER :: nlast = 148 * 11   ! number of MIFs to test at end

! Save Cnts for calibration:
    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at start of THzCalDay'

    dCnts = Cnts

! Set CalFlags:

    CalFlag = 0   ! Nothing yet
    WHERE (CalInfo /= 0)
       CalFlag = 1       ! Cal MIFs
    END WHERE

    nMAFs = CalBuf%MAFs
    nMIFs = SIZE (MIFx)
    norbits = CEILING (nMAFs / 240.0)
    iend = -1
    DO iorbit = 1, norbits
       ibgn = iend + 1
       iorg = ibgn
       iendm = iorbit * 240
       IF (iendm >= nMAFs) THEN
          iend = nMIFs - 1
          ibgn = (nMAFs - 241) * 148
          IF (ibgn <= 0) THEN
             ibgn = 0
          ENDIF
       ELSE
          iend = iendm * 148 - 1
          DO i = 0, nlast
             IF (bias(iend - i) > maxBias) THEN
                iend = iend - i
                EXIT
             ENDIF
          ENDDO
       ENDIF
       iref = iorg - ibgn
PRINT *, ibgn, iend, iref, iorg, (iend-ibgn+1)/148
       maf1 = ibgn / 148
       maf2 = iend / 148
       maf2 = MINVAL ((/ maf2, (nMAFs - 1) /))
       IF (xMIF0(maf2) >= iend) maf2 = maf2 - 1
       IF (iorbit == norbits .AND. iorbit > 1) maf1 = maf1 - 1

!       DEALLOCATE (tCnts, stat=status)
       ALLOCATE (tCnts(THzChans,THzNum,0:(iend-ibgn)))
       tCnts = dCnts(:,:,ibgn:iend)
!       DEALLOCATE (tVarCnts, stat=status)
       ALLOCATE (tVarCnts(THzChans,THzNum,0:(iend-ibgn)))
!       DEALLOCATE (tVargain, stat=status)
       ALLOCATE (tVargain(0:(iend-ibgn)))

       CALL THzBound (ibgn, iend, iorbit)

       CALL THzDel (Bias(ibgn:iend), BiasGood(ibgn:iend), CalFlag(ibgn:iend), &
            CalTemp(ibgn:iend), xMIF0(maf1:maf2), tCnts, tVarCnts, tVargain)

       Cnts(:,:,iorg:iend) = tCnts(:,:,iref:)
       VarCnts(:,:,iorg:iend) = tVarCnts(:,:,iref:)
       VarGain(iorg:iend) = tVargain(iref:)

       CALL OutputL1B_DiagsT (L1BFileInfo%DiagTid, OrbNo=diagno, &
            dLlo=dLlo, yTsys=yTsys)
       diagno = diagno + 1

       DEALLOCATE (tCnts,    stat=status)
       DEALLOCATE (tVarCnts, stat=status)
       DEALLOCATE (tVargain, stat=status)
    ENDDO

    CALL THzBound (0, iend, 0, fillnvbounds=.TRUE.)  ! Bounds for entire dataset

    if ( any(isNaN(Cnts)) ) print *, 'NaNs in Cnts at end of THzCalDay'
  END SUBROUTINE THzCalDay

!=============================================================================
  SUBROUTINE CalibrateTHz
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo
    USE MLSL1Config, ONLY: L1Config
    USE OutputL1B, ONLY: OutputL1B_DiagsT

    INTEGER :: nbad
    LOGICAL :: ColdCal

PRINT *, 'Start calibrating...'

    ColdCal = L1Config%Calib%THzColdCal

    CALL BuildCalVectors

    CALL THzCalDay

    CALL THzCal (ColdCal=ColdCal)

    CALL THzStat

    CALL THzReselect ("C", nbad)        ! Only "C"old cals
    print *, 'nbad ', nbad

    IF (nbad > 0) THEN

       CALL THzCal (ColdCal=ColdCal)

       CALL THzReselect ("HC", nbad)    ! Both "H"ot and "C"old cals

    ENDIF

    CALL RestoreCnts

    CALL THzCalDay

    CALL THzCal (ColdCal=ColdCal)

    CALL THzStat

    CALL THzReselect ("L", nbad)        ! Only "L"imb views

    CALL OutputL1B_DiagsT (L1BFileInfo%DiagTid, OrbNo=1, Chisq=Chisq)

PRINT *, 'End calibrating...'

  END SUBROUTINE CalibrateTHz

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE THzCalibration
!=============================================================================

! $Log$
! Revision 2.19  2015/01/22 23:34:04  vsnyder
! Get constants from Constants module instead of SDPToolkit
!
! Revision 2.18  2015/01/21 19:32:51  pwagner
! Avoid array bound violations
!
! Revision 2.17  2015/01/14 00:29:44  pwagner
! Warn of NaNs
!
! Revision 2.16  2013/02/06 19:05:40  quyen
! fix hanging problem for misplaced deallocate
!
! Revision 2.15  2012/08/29 17:11:37  perun
! Put in good case tests for Bands 15 and 20.
!
! Revision 2.14  2009/09/03 19:07:33  perun
! Use correct number of MAFs in THzBound routine
!
! Revision 2.13  2009/08/21 18:59:13  perun
! Set calibration flag to false when no cal data is available.
!
! Revision 2.12  2006/08/22 18:39:40  perun
! Initialize tt array to 0.0
!
! Revision 2.11  2006/04/05 18:08:40  perun
! Add SAVE for NAG compiler and remove unused variables
!
! Revision 2.10  2006/03/24 15:19:57  perun
! Rewrote most of this module to do "C"old calibrations for entire day instead of per orbit
!
! Revision 2.9  2005/12/06 19:30:19  perun
! Removed BrightObjest_T and added BO_stat to MAFdata_T
!
! Revision 2.8  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2004/11/10 15:38:36  perun
! Add first_fit flag per HMP
!
! Revision 2.6  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.5  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.4  2004/01/13 17:15:45  perun
! Protect from arithmetic exception.
!
! Revision 2.3  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.2  2003/02/05 21:31:55  perun
! Calculate variance and chi square
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
