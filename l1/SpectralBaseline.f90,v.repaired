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
MODULE SpectralBaseline ! Determine Spectral baseline radiances 
!=============================================================================

  USE MLSL1Common, ONLY: BandChans, NumBands, L1BFileInfo, FBchans, MBchans, &
       WFchans, DACSchans, Bandwidth, LO1, slimb_type
  USE L1BData, ONLY: L1BData_T, ReadL1BData, DeallocateL1BData

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitBaseline, CalcBaseline, UpdateBaselines, LoadBaselineAC
  PUBLIC :: Baseline_T, Baseline, BaselineAC, BaselineDC, BaselineInclude, &
       BaselineInclude_T

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE Baseline_T
     REAL, DIMENSION(:), POINTER :: offset, PRECISION
  END TYPE Baseline_T

  TYPE (Baseline_T) :: Baseline(NumBands), BaselineAC(NumBands), &
       BaselineDC(NumBands), BaselineAC_deflt(NumBands)

  TYPE BaselineInclude_T
     LOGICAL, DIMENSION(:), POINTER :: chan
  END TYPE BaselineInclude_T

  TYPE (BaselineInclude_T) :: BaselineInclude(NumBands)

  REAL, PARAMETER :: BaselineAlt = 80.0e3   ! Minimum baseline altitude (m)
  REAL, PARAMETER :: BaselineAltDACS(5) = &
       (/ 82.0e3, 80.0e3, 72.0e3, 72.0e3, 82.0e3 /)

! Bandwidths:

  REAL, PARAMETER :: FB_BW(FBchans) = Bandwidth%FB(:,1) * 1.0e-06
  REAL, PARAMETER :: MB_BW(MBchans) = Bandwidth%MB(:,1) * 1.0e-06
  REAL, PARAMETER :: DACS_BW = Bandwidth%DACS(1,1) * 1.0e-06
  REAL :: BW(FBchans)

CONTAINS

!=============================================================================
  SUBROUTINE InitBaseline
!=============================================================================

    INTEGER :: i

    INTEGER, PARAMETER :: BadChans(2,NumBands) = RESHAPE ( (/ &
         8,18, 12,14,   2,2,   0,0,   0,0,   0,0, 13,13, 12,14, 13,13,   0,0, &
          0,0,   0,0, 11,15, 13,13, 11,17, 11,17, 12,14, 11,17, 11,17, 12,14, &
         8,18, 25,84, 36,69, 33,74, 32,71, 25,84,   0,0,   6,6,   0,0,   6,6, &
          0,0,   0,0,   0,0,   0,0 /), (/ 2, NumBands /) )

    DO i = 1, NumBands

       ALLOCATE (Baseline(i)%offset(BandChans(i)))
       ALLOCATE (Baseline(i)%precision(BandChans(i)))
       ALLOCATE (BaselineAC(i)%offset(BandChans(i)))
       ALLOCATE (BaselineAC(i)%precision(BandChans(i)))
       ALLOCATE (BaselineAC_deflt(i)%offset(BandChans(i)))
       ALLOCATE (BaselineAC_deflt(i)%precision(BandChans(i)))
       ALLOCATE (BaselineDC(i)%offset(BandChans(i)))
       ALLOCATE (BaselineDC(i)%precision(BandChans(i)))

! Channels to include/exclude:

       ALLOCATE (BaselineInclude(i)%chan(BandChans(i)))
       BaselineInclude(i)%chan = .TRUE.   ! Mark all channels as "good"
       IF (BadChans(1,i) /= 0) THEN   ! Check for "bad" channel(s)
          BaselineInclude(i)%chan(BadChans(1,i):BadChans(2,i)) = .FALSE. !"bad"
       ENDIF

    ENDDO

! Mark special band channels bad

    BaselineInclude(13)%chan(25) = .FALSE.      ! Band 13
    BaselineInclude(14)%chan(24:25) = .FALSE.   ! Band 14
    DO i = 22, 26                               ! Bands 22-26
       BaselineInclude(i)%chan(1) = .FALSE.
       BaselineInclude(i)%chan(108:129) = .FALSE.
    ENDDO

! Clear the baselines:

    CALL ClearBaselines

  END SUBROUTINE InitBaseline

!=============================================================================
  SUBROUTINE ClearBaselines
!=============================================================================

    INTEGER :: i

    DO i = 1, NumBands

       Baseline(i)%offset = 0.0
       Baseline(i)%precision = 0.0
       BaselineAC(i)%offset = 0.0
       BaselineAC(i)%precision = 0.0
       BaselineDC(i)%offset = 0.0
       BaselineDC(i)%precision = 0.0

    ENDDO

  END SUBROUTINE ClearBaselines

!=============================================================================
  SUBROUTINE CalcBaselineAC (bandno)
!=============================================================================

    INTEGER, INTENT(IN) :: bandno

! Use default AC baseline for now:

    BaselineAC(bandno)%offset = BaselineAC_deflt(bandno)%offset
    BaselineAC(bandno)%precision = BaselineAC_deflt(bandno)%precision

  END SUBROUTINE CalcBaselineAC

!=============================================================================
  SUBROUTINE CalcBaselineDC (bandno, radNum, chans, MIFmatch, MIFmatchDACS, &
       rad, prec)
!=============================================================================

    USE MLSL1Rad, ONLY : RadPwr
    USE MLSL1Config, ONLY: L1Config

    INTEGER, INTENT(IN) :: bandno, chans, radNum
    LOGICAL, INTENT(IN) :: MIFmatch(:), MIFmatchDACS(:,:)
    REAL, INTENT(IN) :: rad(:,:), prec(:,:)

    INTEGER :: i, nmatch
    INTEGER :: MIFindx(SIZE(MIFmatch))
    LOGICAL :: match(SIZE(MIFmatch)), AltMatch(SIZE(MIFmatch))
    REAL :: rad_sum, bw_sum, space_P

    BaselineDC(bandno)%offset = 0.0     ! Clear this band

    IF (chans /= DACSchans) THEN
       AltMatch = MIFmatch
    ELSE
       AltMatch = MIFmatchDACS(:,(bandno-21))
    ENDIF
    IF (COUNT (AltMatch) == 0) RETURN   ! nothing to do

    rad_sum = 0.0
    bw_sum = 0.0
    MIFindx = 1
    IF (chans == FBchans) THEN
       BW = FB_BW
    ELSE IF (chans == MBchans) THEN
       BW(1:MBchans) = MB_BW
    ENDIF

    space_P = radPwr (LO1(radNum), L1Config%Calib%GHzSpaceTemp)  
    DO i = 1, chans

       IF (.NOT. BaselineInclude(bandno)%chan(i)) CYCLE   ! next one

       WHERE (AltMatch .AND. prec(i,:) >= 0.0)
          match = .TRUE.
       ELSEWHERE
          match = .FALSE.
       ENDWHERE
       nmatch = COUNT (match)
       IF (nmatch > 0) THEN

          IF (chans  == WFchans) THEN   ! WF
             BaselineDC(bandno)%offset(i) = &
                  SUM (rad(i,:), mask=match) / nmatch - space_P
          ELSE IF (chans == DACSchans) THEN  ! DACS
             rad_sum = rad_sum + SUM (rad(i,:) * DACS_BW, mask=match)
             bw_sum = bw_sum + SUM (MIFindx(:) * DACS_BW, mask=match)
          ELSE
             rad_sum = rad_sum + SUM (rad(i,:) * BW(i), mask=match)
             bw_sum = bw_sum + SUM (MIFindx(:) * BW(i), mask=match)
          ENDIF
       ENDIF

    ENDDO

    IF (bw_sum > 0.0) BaselineDC(bandno)%offset = rad_sum / bw_sum - space_P

  END SUBROUTINE CalcBaselineDC

!=============================================================================
  SUBROUTINE CalcSlimbOffsets (chans, MIFs, baseline, bandwidth, rad, prec)
!=============================================================================

    INTEGER, INTENT(IN) :: chans, MIFs
    REAL, INTENT(IN) :: baseline(chans), bandwidth(chans)
    REAL, INTENT(IN) :: rad(chans,MIFs), prec(chans,MIFs)

    integer :: MIFno, nslimbs, channo
    logical :: slimb_type(chans)
    real :: slimb_sum, sw_sum, slimb_bw, sw_bw

    print *, 'slimb offsets'
    slimb_type = .false.
    where (baseline == 0.0)
       slimb_type = .true.
    end where
    nslimbs = count (slimb_type)
    if (nslimbs == chans) return         ! Nothing to do

    if (nslimbs > 0) then
       print *, 'rad: ', rad(:,10)
       print *, 'prec: ',  prec(:,10)
       print *, 'baseline', baseline
       print *, 'bandwidth: ', bandwidth
       print *, 'nslimbs, MIFs ', nslimbs, MIFs
       slimb_sum = 0.0
       slimb_bw = 0.0
       sw_sum = 0.0
       sw_bw = 0.0
       do MIFno = 1, MIFs
!!$          slimb_sum = slimb_sum + sum (rad(:,MIFno) * bandwidth(:), &
!!$               mask=(slimb_type .and. (prec(:,MIFno) >= 0.0)))
!!$          slimb_bw = slimb_bw + sum (bandwidth(:), &
!!$               mask=(slimb_type .and. (prec(:,MIFno) >= 0.0)))
!!$
!!$          sw_sum = sw_sum + sum ((rad(:,MIFno)+baseline(:)) * bandwidth(:), &
!!$               mask=(.not. slimb_type .and. (prec(:,MIFno) >= 0.0)))
!!$          sw_bw = sw_bw + sum (bandwidth(:), &
!!$               mask=(.not. slimb_type .and. (prec(:,MIFno) >= 0.0)))

          do channo = 1, chans
             if (slimb_type(channo) .and. prec(channo,MIFno) >= 0.0) &
                  slimb_sum = slimb_sum + rad(channo,MIFno)*bandwidth(channo)
             if (slimb_type(channo) .and. prec(channo,MIFno) >= 0.0) &
                  slimb_bw = slimb_bw + bandwidth(channo)
             if (.not. slimb_type(channo) .and. prec(channo,MIFno) >= 0.0) &
                  sw_sum = sw_sum + (rad(channo,MIFno)+baseline(channo)) * &
                  bandwidth(channo)
             if (.not. slimb_type(channo) .and. prec(channo,MIFno) >= 0.0) &
                  sw_bw = sw_bw + bandwidth(channo)
          enddo
       enddo
    endif

  END SUBROUTINE CalcSlimbOffsets

!=============================================================================
  SUBROUTINE CalcBaseline
!=============================================================================

    USE MLSL1Config, ONLY: MIFsGHz
    USE Calibration, ONLY: CalWin, MAFdata_T
    USE MLSL1Rad, ONLY : L1Brad

    INTEGER :: i, bandno, chans, radNum, spectNum
    LOGICAL :: MIFmatch(MIFsGHz)   ! "good" matching MIF altitudes for non-DACS
    LOGICAL :: MIFmatchDACS(MIFsGHz, 5) ! "good" matching MIF altitudes for DACS
    TYPE (MAFdata_T), POINTER :: CurMAFdata
    REAL :: alt(MIFsGHz)

PRINT *, 'doing baseline...'

    CurMAFdata => CalWin%MAFdata(CalWin%central)
    alt = CurMAFdata%SciPkt(0:(MIFsGHz-1))%altg

    WHERE (alt > BaselineAlt)   ! non-DACS bands...
       MIFmatch = .TRUE.
    ELSEWHERE
       MIFmatch = .FALSE.
    ENDWHERE

    DO i = 1, 5
       WHERE (alt > BaselineAltDACS(i))
          MIFmatchDACS(:,i) = .TRUE.
       ELSEWHERE
          MIFmatchDACS(:,i) = .FALSE.
       ENDWHERE
    ENDDO

    CALL ClearBaselines

    DO i = 1, SIZE (L1Brad)

       bandno = L1Brad(i)%bandno
       chans = SIZE (L1Brad(i)%value(:,1))
       radNum = l1brad(i)%signal%radiometerNumber
       spectNum = l1brad(i)%signal%spectrometerNumber

       CALL CalcBaselineDC (bandno, radNum, chans, MIFmatch, MIFmatchDACS, &
            L1Brad(i)%value, L1Brad(i)%precision)

       CALL CalcBaselineAC (bandno)

! Clear AC/DC baselines for slimb_types:

       SELECT CASE (Chans)
       CASE (FBchans)
          WHERE (slimb_type%FB(:,spectNum))
             BaselineAC(bandno)%offset = 0.0
             BaselineDC(bandno)%offset = 0.0
          END WHERE
       CASE (MBchans)
          WHERE (slimb_type%MB(:,spectNum))
             BaselineAC(bandno)%offset = 0.0
             BaselineDC(bandno)%offset = 0.0
          END WHERE
       CASE (WFchans)
          WHERE (slimb_type%WF(:,spectNum))
             BaselineAC(bandno)%offset = 0.0
             BaselineDC(bandno)%offset = 0.0
          END WHERE
       CASE (DACSchans)
          WHERE (slimb_type%DACS(:,spectNum))
             BaselineAC(bandno)%offset = 0.0
             BaselineDC(bandno)%offset = 0.0
          END WHERE
       END SELECT

! Combine AC/DC and make baseline negative:

       Baseline(bandno)%offset = -BaselineDC(bandno)%offset - &
            BaselineAC(bandno)%offset
       Baseline(bandno)%offset = 0.0

    ENDDO

 END SUBROUTINE CalcBaseline

!=============================================================================
  Subroutine UpdateBaselines
!=============================================================================

    USE MLSCommon, ONLY: DEFAULTUNDEFINEDVALUE
    USE MLSL1Config, ONLY: L1Config, MIFsGHz
    USE MLSL1Rad, ONLY: Rad_name
    USE MLSL1Config, ONLY: L1Config
    USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
    USE MLSAuxData, ONLY: Build_MLSAuxData
    USE MLSHDF5, ONLY: IsHDF5DSPresent

    REAL, PARAMETER :: FILLVALUE = DEFAULTUNDEFINEDVALUE
    REAL, PARAMETER :: DC_min = -2.0     ! Minimum DC avg threshold

    CHARACTER(LEN=80) :: name, BaseName, BasePrecName, BaseNameAC, &
         BasePrecNameAC, BaseNameDC, BasePrecNameDC
    INTEGER :: i, bandno, c1, c2, Flag, noChans, noMAFs, ngood, nwin, sd_id
    INTEGER :: CalWindow, CalStart, CalEnd, mindx, avg_indx(2), rno, w1, w2
    INTEGER :: DACS_window
    INTEGER, POINTER, DIMENSION(:,:) :: windx
    LOGICAL, POINTER, DIMENSION(:) :: rms_mask
    LOGICAL, POINTER, DIMENSION(:,:) :: avg_mask
    LOGICAL :: read_precs, update_precs

    INTEGER, POINTER, DIMENSION(:) :: counterMAF => NULL()
    REAL, POINTER, DIMENSION(:,:) :: baselineAC => NULL(), &
         baselineACprec => NULL(), baselineDC => NULL(), &
         baselineDCavg => NULL()
    REAL DC_avg(FBchans), DC_rms(FBchans), DCmean
    REAL, POINTER, DIMENSION(:) :: resid
    REAL, POINTER, DIMENSION(:,:,:) :: rad => NULL(), prec => NULL()

    TYPE (L1BData_T) :: L1BData
    TYPE (DataProducts_T) :: baselineDS, dataset  

PRINT *, 'Updating baselines...'
    DC_avg = 0.

! Set up for HDF output:

    CALL Deallocate_DataProducts (baselineDS)
    ALLOCATE (baselineDS%Dimensions(2))
    baselineDS%data_type = 'real'
    baselineDS%Dimensions(2) = 'MAF                 '
    baselineDS%Dimensions(1) = 'chanDACS'

    ALLOCATE (dataset%Dimensions(3))  ! Only for precision for now
    dataset%data_type = 'real'
    dataset%Dimensions(1) = 'chanDACS            '
    dataset%Dimensions(2) = 'GHz.MIF             '
    dataset%Dimensions(3) = 'MAF                 '

    sd_id = L1BFileInfo%RADGid

! Get counterMAFs

    CALL ReadL1BData (sd_id, 'counterMAF', L1BData, noMAFs,  Flag, &
         NeverFail=.TRUE., HDFversion=5)
    ALLOCATE (counterMAF(noMAFs))
    counterMAF = L1BData%Intfield(1,1,:)
    CALL DeallocateL1BData (L1BData)

! DACS baselines first:

    IF (L1Config%Calib%CalibDACS) THEN
       sd_id = L1BFileInfo%RADDid
       DACS_window = L1Config%Calib%DACSwindow
       ALLOCATE (avg_mask(DACSchans,2*DACS_window + 1))

! Allocate baseline arrays for reading to the maximum channel size:

       ALLOCATE (baselineDC(DACSchans,noMAFs))
       ALLOCATE (baselineDCavg(DACSchans,noMAFs))
       ALLOCATE (prec(DACSchans,MIFsGHz,noMAFs))

! Get the DACS baselines and adjust

       DO rno = 40, 44   ! DACS only names!

          read_precs = .TRUE.
          update_precs = .FALSE.
          name = Rad_name(rno)
          IF (.NOT. IsHDF5DSPresent (sd_id, TRIM(name))) CYCLE  ! Nothing to get

          BaseName = TRIM (name) // ' Baseline'
          BaseNameDC = TRIM (name) // ' BaselineDC'

          CALL ReadL1BData (sd_id, BaseNameDC, L1BData, noMAFs, Flag, &
               NeverFail=.TRUE., HDFversion=5)
          baselineDC = L1BData%DpField(1,:,:)
          CALL DeallocateL1BData (L1BData)

          DO mindx = 1, noMAFs
             w1 = MAX ((mindx - DACS_window), 1)
             w2 = MIN ((mindx + DACS_window), noMAFs)
             nwin = w2 - w1 + 1
             avg_mask = .TRUE.      ! everything good
             WHERE (baselineDC(:,w1:w2) <= -999.0)
                avg_mask(:,1:w2-w1+1) = .FALSE.  ! bad data to skip
             ENDWHERE
             ngood = COUNT (avg_mask(1,:))   ! using just one channel!
             if (ngood < 1) CYCLE
             baselineDCavg(:,mindx) = &
                  SUM (baselineDC(:,w1:w2),2,avg_mask(:,1:nwin)) / ngood
             baselineDS%name = BaseNameDC
             CALL Build_MLSAuxData (sd_id, baselineDS, &
                  BaselineDCavg(:,mindx), lastIndex=mindx, &
                  disable_attrib=.TRUE.)
             baselineDS%name = BaseName
             CALL Build_MLSAuxData (sd_id, baselineDS, &
                  -BaselineDCavg(:,mindx), lastIndex=mindx, &
                  disable_attrib=.TRUE.)
             update_precs = (ANY (BaselineDCavg(:,mindx) < DC_min))

             IF (read_precs .AND. update_precs) THEN

! Get precisions (to adjust when needed!)

                CALL ReadL1BData (sd_id, TRIM(name)//' precision', L1BData, &
                     noMAFs, Flag, NeverFail=.TRUE., HDFversion=5)
                prec = L1BData%DpField
                read_precs = .FALSE.
                CALL DeallocateL1BData (L1BData)
             ENDIF
             IF (update_precs) THEN
                dataset%name = TRIM(name)//' precision'
                prec(:,:,mindx) = -1.0 * ABS (prec(:,:,mindx))  ! Negate precs
                CALL Build_MLSAuxData (sd_id, dataset, &
                     prec(:,:,mindx), lastIndex=mindx, &
                     disable_attrib=.TRUE.)
             ENDIF

          ENDDO

       ENDDO

       DEALLOCATE (baselineDC)
       DEALLOCATE (baselineDCavg)
       DEALLOCATE (prec)

    ENDIF

! Now for RADG baselines:

    sd_id = L1BFileInfo%RADGid

! Determine Cal Window indexes

    CalWindow = L1Config%Calib%CalWindow
    CalStart = -CalWindow / 2
    CalEnd = CalWindow / 2 - 1
    ALLOCATE (windx(2,noMAFs))
    DO mindx = 1, noMAFs
       w1 = MAX ((mindx+CalStart), 1)
       w2 = MIN ((mindx+CalEnd), noMAFs)
       IF (counterMAF(mindx) <= 0) THEN
          w1 = -1; w2 = -1
       ELSE
          DO
             IF ((mindx - w1) == (counterMAF(mindx) - counterMAF(w1))) EXIT
             w1 = w1 + 1
          ENDDO
          DO
             IF ((w2 - mindx) == (counterMAF(w2) - counterMAF(mindx))) EXIT
             w2 = w2 - 1
          ENDDO
          w1 = MIN (w1, noMAFs); w2 = MAX (w2, 1)  ! keep within limits
       ENDIF
       windx(:,mindx) = (/ w1, w2 /)
    ENDDO

    ALLOCATE (rms_mask(CalWindow))
    ALLOCATE (resid(CalWindow))

! Allocate baseline arrays for reading to the maximum channel size:

    ALLOCATE (baselineAC(FBchans,noMAFs))
    ALLOCATE (baselineACprec(FBchans,noMAFs))
    ALLOCATE (baselineDC(FBchans,noMAFs))
    ALLOCATE (rad(FBchans,MIFsGHz,noMAFs))
    ALLOCATE (prec(FBchans,MIFsGHz,noMAFs))

! Get the GHz baselines and adjust

    DO rno = 1, SIZE (Rad_name)

       name = Rad_name(rno)

       IF (.NOT. IsHDF5DSPresent (sd_id, TRIM(name))) CYCLE  ! Nothing to get
       IF (INDEX (name, 'DACS') /= 0) CYCLE ! Nothing for the DACS

! DS names to read/write:

       BaseName = TRIM (name) // ' Baseline'
       BasePrecName = TRIM (name) // ' Baseline precision'
       BaseNameAC = TRIM (name) // ' BaselineAC'
       BasePrecNameAC = TRIM (name) // ' BaselineAC precision'
       BaseNameDC = TRIM (name) // ' BaselineDC'
       BasePrecNameDC = TRIM (name) // ' BaselineDC precision'

       CALL ReadL1BData (sd_id, BaseNameAC, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       IF (Flag /= 0) CYCLE   ! Can't do anything more...

       c1 = INDEX (name, '.B') + 2
       c2 = INDEX (name(c1:), ':') - 3 + c1
       READ (name(c1:c2), *) bandno
       noChans = L1BData%MaxMIFs

       baselineAC(1:noChans,:) = L1BData%DpField(1,:,:)
       CALL DeallocateL1BData (L1BData)

       CALL ReadL1BData (sd_id, BasePrecNameAC, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       baselineACprec(1:noChans,:) = L1BData%DpField(1,:,:)
       CALL DeallocateL1BData (L1BData)

       CALL ReadL1BData (sd_id, BaseNameDC, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       baselineDC(1:noChans,:) = L1BData%DpField(1,:,:)
       CALL DeallocateL1BData (L1BData)

! set channel type for DS:

       SELECT CASE (noChans)
       CASE (FBchans)
          baselineDS%Dimensions(1) = 'chanFB'
          BW = FB_BW
       CASE (MBchans)
          baselineDS%Dimensions(1) = 'chanMB'
          BW(1:MBchans) = MB_BW
       CASE (WFchans)
          baselineDS%Dimensions(1) = 'chanWF'
       CASE (DACSchans)
          baselineDS%Dimensions(1) = 'chanDACS'
       END SELECT

! Get Radiances:

       CALL ReadL1BData (sd_id, name, L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       rad(1:noChans,:,:) = L1BData%DpField
       CALL DeallocateL1BData (L1BData)

       CALL ReadL1BData (sd_id, TRIM(name)//' precision', L1BData, noMAFs, &
            Flag, NeverFail=.TRUE., HDFversion=5)

       prec(1:noChans,:,:) = L1BData%DpField
       CALL DeallocateL1BData (L1BData)

       DO mindx = 1, noMAFs   ! calculate and output baselines

! DC average first

          avg_indx = (/ mindx, (mindx - 1) /) ! assume can use previous

          IF (mindx == 1) THEN  ! can't use previous if first MAF in file
             avg_indx(2) = mindx
          ELSE
             IF (counterMAF(mindx) /= (counterMAF(mindx-1) + 1)) THEN
 
! CounterMAF must increment by 1 to be good

               avg_indx(2) = mindx
             ELSE IF (ANY (baselineDC(1:noChans,(mindx-1)) == FILLVALUE) .OR. &
                  ANY (baselineDC(1:noChans,(mindx-1)) == 0.0)) THEN

! Fill value in previous can't be used for averaging

                avg_indx(2) = mindx

             ENDIF
          ENDIF
          DC_avg(1:noChans) = SUM (baselineDC(1:noChans,avg_indx),2) * 0.5

          IF (ANY (DC_avg(1:noChans) == FILLVALUE)) DC_avg = 0.0 ! Nothing there

! DC RMS

          DC_rms(1:noChans) = 0.0
          w1 = windx(1,mindx); w2 = windx(2,mindx); nwin = w2 - w1 + 1
          IF (w1 > 0) THEN
             DO i = 1, noChans
                resid = 0.0
                rms_mask = .FALSE.    ! assume all bad
                WHERE (baselineDC(i,w1:w2) /= FILLVALUE)   ! check for FILLs
                   rms_mask(1:nwin) = .TRUE.
                ENDWHERE
                ngood = COUNT (rms_mask)
                IF (ngood > 1) THEN
                   DCmean = SUM (baselineDC(i,w1:w2), rms_mask(1:nwin)) &
                        / ngood
                   WHERE (rms_mask(1:nwin))
                      resid(1:nwin) = baselineDC(i,w1:w2) - DCmean
                   ENDWHERE
                   DC_rms(i) = ((SUM (resid**2) - SUM (resid)**2 / ngood) / &
                        (ngood - 1))**0.5
                ENDIF
             ENDDO
          ENDIF

! Output updated baselines:

          IF (ANY (baselineAC(1:noChans,mindx) == FILLVALUE)) THEN

             Baseline(bandno)%offset = FILLVALUE
          ELSE

             Baseline(bandno)%offset = -baselineAC(1:noChans,mindx) - &
                  DC_avg(1:noChans)

! Reset possible "slimb" type data channels (AC and DC both 0.0):

             WHERE (baselineAC(1:noChans,mindx) == 0.0 .AND. &
                  baselineDC(1:noChans,mindx) == 0.0)
                Baseline(bandno)%offset = 0.0
             END WHERE
          ENDIF

!          CALL CalcSlimbOffsets (noChans, MIFsGHz, baseline(bandno)%offset, &
!               BW(1:noChans), rad(1:noChans,:,mindx), prec(1:noChans,:,mindx))

          baselineDS%name = BaseName
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               Baseline(bandno)%offset, lastIndex=mindx, &
               disable_attrib=.TRUE.)

          IF (ANY (baselineAC(1:noChans,mindx) == FILLVALUE)) THEN
             Baseline(bandno)%precision = FILLVALUE
          ELSE
             Baseline(bandno)%precision = SQRT ( &
                  baselineACprec(1:noChans,mindx)**2 + &
                  DC_rms(1:noChans)**2)
          ENDIF
          baselineDS%name = BasePrecName
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               Baseline(bandno)%precision, lastIndex=mindx, &
               disable_attrib=.TRUE.)

          baselineDS%name = BaseNameDC
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               DC_avg(1:noChans), lastIndex=mindx, &
               disable_attrib=.TRUE.)

          baselineDS%name = BasePrecNameDC
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               DC_rms(1:noChans), lastIndex=mindx, &
               disable_attrib=.TRUE.)

       ENDDO
    ENDDO

    DEALLOCATE (baselineAC)
    DEALLOCATE (baselineACprec)
    DEALLOCATE (baselineDC)
    DEALLOCATE (rad)
    DEALLOCATE (prec)

  END SUBROUTINE UpdateBaselines

!=============================================================================
  SUBROUTINE LoadBaselineAC (unit, stat)
!=============================================================================

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT) :: stat

    CHARACTER(LEN=80) :: line
    INTEGER :: i, ios

    stat = 0

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:6) == "#DATA") EXIT
    ENDDO

! Read data

    DO i = 1, NumBands

       IF (i >= 22 .AND. i <= 26) CYCLE   ! Nothing for the DACS

       DO
          READ (unit, '(A)', IOSTAT=ios) line
          IF (line(1:1) == "#" .OR. ios /= 0) EXIT
       ENDDO

       IF (ios /= 0) EXIT

       READ (unit, *) BaselineAC_deflt(i)%offset(1:BandChans(i))
       READ (unit, *) BaselineAC_deflt(i)%precision(1:BandChans(i))

    ENDDO

  END SUBROUTINE LoadBaselineAC

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE SpectralBaseline
!=============================================================================
! $Log$
! Revision 2.14  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.13.4.2  2016/03/14 19:51:24  whdaffer
! Most of the work is to eliminate cicular references between MLSL1Debug
! and Calibration. To resolve this, I've moved the inclusion of
! Calibration.f9h and the definition of some 10 variables from
! Calibration.f90 to MLSL1Common.f90. Radiances, MLSL1Debug and
! Calibration will get those types, variable from MLSL1Common. Also, use
! machines.f90 to get the definition of usleep used in SnoopMLSL1
!
! Revision 2.13.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.13  2015/01/21 19:32:29  pwagner
! Avoid array bound violations
!
! Revision 2.12  2011/01/27 15:37:20  perun
! Mask bad DC baseline for calculating average DC baseline.
!
! Revision 2.11  2010/04/08 20:32:36  perun
! Fixed AC/DC baseline array tests.
!
! Revision 2.10  2007/06/21 21:05:24  perun
! Correct rad and prec assignments
!
! Revision 2.9  2007/02/09 15:06:51  perun
! Test DACS precisions again minimum DC_avg of -2.0
!
! Revision 2.8  2006/09/26 16:03:28  perun
! Mark DC baseline values as 0.0 when not available
!
! Revision 2.7  2006/08/02 18:59:10  perun
! Do not calculation slimb offsets
!
! Revision 2.6  2006/04/05 18:10:47  perun
! Remove unused variables
!
! Revision 2.5  2006/03/24 15:19:20  perun
! Set Space/Limb data baselines to 0.0
!
! Revision 2.4  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2005/05/02 16:09:11  perun
! Deallocate baseline pointers
!
! Revision 2.2  2004/12/01 17:11:34  perun
! Add calculating and outputting DACS DC baseline
!
! Revision 2.1  2004/11/10 15:31:21  perun
! Initial commit
!
