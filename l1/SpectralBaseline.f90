! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SpectralBaseline ! Determine Spectral baseline radiances 
!=============================================================================

  USE MLSL1Common, ONLY: BandChans, NumBands
  USE L1BData, ONLY: L1BData_T, ReadL1BData, DeallocateL1BData
  USE MLSL1Rad, ONLY: Rad_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitBaseline, CalcBaseline, UpdateBaselines, LoadBaselineAC
  PUBLIC :: Baseline_T, Baseline, BaselineAC, BaselineDC 

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE Baseline_T
     REAL, DIMENSION(:), POINTER :: offset, precision
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

    USE MLSL1Common, ONLY: Bandwidth, DACSchans, FBchans, MBchans, WFchans, LO1
    USE MLSL1Rad, ONLY : RadPwr
    USE MLSL1Config, ONLY: L1Config

    INTEGER, INTENT(IN) :: bandno, chans, radNum
    LOGICAL, INTENT(IN) :: MIFmatch(:), MIFmatchDACS(:,:)
    REAL, INTENT(IN) :: rad(:,:), prec(:,:)

    INTEGER :: i, nmatch
    INTEGER :: MIFindx(SIZE(MIFmatch))
    LOGICAL :: match(SIZE(MIFmatch)), AltMatch(SIZE(MIFmatch))
    REAL :: rad_sum, bw_sum, BW(FBchans), space_P
    REAL, PARAMETER :: FB_BW(FBchans) = Bandwidth%FB(:,1) * 1.0e-06
    REAL, PARAMETER :: MB_BW(MBchans) = Bandwidth%MB(:,1) * 1.0e-06
    REAL, PARAMETER :: DACS_BW = Bandwidth%DACS(1,1) * 1.0e-06

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
  SUBROUTINE CalcBaseline
!=============================================================================

    USE MLSL1Config, ONLY: MIFsGHz
    USE Calibration, ONLY: CalWin, MAFdata_T
    USE MLSL1Rad, ONLY : L1Brad

    INTEGER :: i, bandno, chans, radNum
    LOGICAL :: MIFmatch(MIFsGHz)   ! "good" matching MIF altitudes for non-DACS
    LOGICAL :: MIFmatchDACS(MIFsGHz, 5) ! "good" matching MIF altitudes for DACS
    TYPE (MAFdata_T), POINTER :: CurMAFdata
    REAL :: alt(MIFsGHz)

print *, 'doing baseline...'

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

       CALL CalcBaselineDC (bandno, radNum, chans, MIFmatch, MIFmatchDACS, &
            L1Brad(i)%value, L1Brad(i)%precision)

       CALL CalcBaselineAC (bandno)

! Combine AC/DC and make baseline negative:

       Baseline(bandno)%offset = -BaselineDC(bandno)%offset - &
            BaselineAC(bandno)%offset

    ENDDO

  END SUBROUTINE CalcBaseline

!=============================================================================
  SUBROUTINE UpdateBaselines
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, FBchans, MBchans, WFchans, DACSchans
    USE MLSCommon, ONLY: DEFAULTUNDEFINEDVALUE
    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Rad, ONLY: Rad_name
    USE MLSL1Config, ONLY: L1Config
    USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
    USE MLSAuxData, ONLY: Build_MLSAuxData
    USE MLSHDF5, ONLY: IsHDF5DSPresent

    REAL, PARAMETER :: FILLVALUE = DEFAULTUNDEFINEDVALUE

    CHARACTER(LEN=80) :: name, BaseName, BasePrecName, BaseNameAC, &
         BasePrecNameAC, BaseNameDC, BasePrecNameDC
    INTEGER :: i, bandno, c1, c2, Flag, noChans, noMAFs, ngood, nwin, sd_id
    INTEGER :: CalWindow, CalStart, CalEnd, mindx, avg_indx(2), rno, w1, w2
    INTEGER :: DACS_window
    INTEGER, POINTER, DIMENSION(:,:) :: windx
    LOGICAL, POINTER, DIMENSION(:) :: rms_mask

    INTEGER, POINTER, DIMENSION(:) :: counterMAF => NULL()
    REAL, POINTER, DIMENSION(:,:) :: baselineAC => NULL(), &
         baselineACprec => NULL(), baselineDC => NULL(), &
         baselineDCavg => NULL()
    REAL DC_avg(FBchans), DC_rms(FBchans), DCmean
    REAL, POINTER, DIMENSION(:) :: resid

    TYPE (L1BData_T) :: L1BData
    TYPE (DataProducts_T) :: baselineDS  

print *, 'Updating baselines...'

! Set up for HDF output:

    CALL Deallocate_DataProducts (baselineDS)
    ALLOCATE (baselineDS%Dimensions(2))
    baselineDS%data_type = 'real'
    baselineDS%Dimensions(2) = 'MAF                 '
    baselineDS%Dimensions(1) = 'chanDACS'

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

! Allocate baseline arrays for reading to the maximum channel size:

       ALLOCATE (baselineDC(DACSchans,noMAFs))
       ALLOCATE (baselineDCavg(DACSchans,noMAFs))

! Get the DACS baselines and adjust

       DO rno = 40, 44   ! DACS only names!

          name = Rad_name(rno)
          IF (.NOT. IsHDF5DSPresent (sd_id, TRIM(name))) CYCLE  ! Nothing to get

          BaseName = TRIM (name) // ' Baseline'
          BaseNameDC = TRIM (name) // ' BaselineDC'

          CALL ReadL1BData (sd_id, BaseNameDC, L1BData, noMAFs, Flag, &
               NeverFail=.TRUE., HDFversion=5)
          baselineDC = L1BData%DpField(1,:,:)

          DO mindx = 1, noMAFs
             w1 = MAX ((mindx - DACS_window), 1)
             w2 = MIN ((mindx + DACS_window), noMAFs)
             nwin = w2 - w1 + 1
             baselineDCavg(:,mindx) = SUM (baselineDC(:,w1:w2),2) / nwin
             baselineDS%name = BaseNameDC
             CALL Build_MLSAuxData (sd_id, baselineDS, &
                  BaselineDCavg(:,mindx), lastIndex=mindx, &
                  disable_attrib=.TRUE.)
             baselineDS%name = BaseName
             CALL Build_MLSAuxData (sd_id, baselineDS, &
                  -BaselineDCavg(:,mindx), lastIndex=mindx, &
                  disable_attrib=.TRUE.)
          ENDDO
       ENDDO

       CALL DeallocateL1BData (L1BData)
       DEALLOCATE (baselineDC)
       DEALLOCATE (baselineDCavg)

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
       CASE (MBchans)
          baselineDS%Dimensions(1) = 'chanMB'
       CASE (WFchans)
          baselineDS%Dimensions(1) = 'chanWF'
       CASE (DACSchans)
          baselineDS%Dimensions(1) = 'chanDACS'
       END SELECT

       DO mindx = 1, noMAFs   ! calculate and output baselines

! DC average first

          avg_indx = (/ mindx, (mindx - 1) /) ! assume can use previous

          IF (mindx == 1) THEN  ! can't use previous if first MAF in file
             avg_indx(2) = mindx
          ELSE
             IF (counterMAF(mindx) /= (counterMAF(mindx-1) + 1)) THEN
 
! CounterMAF must increment by 1 to be good

               avg_indx(2) = mindx
             ELSE IF (ANY (baselineDC(1:noChans,mindx) == FILLVALUE)) THEN

! Fill value in previous can't be used for averaging

                avg_indx(2) = mindx

             ENDIF
          ENDIF
          DC_avg = SUM (baselineDC(:,avg_indx),2) * 0.5

! DC RMS

          DC_rms(1:noChans) = 0.0
          w1 = windx(1,mindx); w2 = windx(2,mindx); nwin = w2 - w1 + 1
          IF (w1 > 0) THEN
             DO i = 1, noChans
                resid = 0.0
                rms_mask = .FALSE.    ! assume all bad
                WHERE (baselineDC(i,w1:w2) /= FILLVALUE)   ! check for FILLs
                   rms_mask = .TRUE.
                ENDWHERE
                ngood = COUNT (rms_mask)
                IF (ngood > 1) THEN
                   DCmean = SUM (baselineDC(i,w1:w2), rms_mask(1:nwin)) &
                        / ngood
                   WHERE (rms_mask)
                      resid = baselineDC(i,w1:w2) - DCmean
                   ENDWHERE
                   DC_rms(i) = ((SUM (resid**2) - SUM (resid)**2 / ngood) / &
                        (ngood - 1))**0.5
                ENDIF

             ENDDO
          ENDIF

! Output updated baselines:

          IF (ANY (baselineAC(:,mindx) == FILLVALUE)) THEN
             Baseline(bandno)%offset = FILLVALUE
          ELSE
             Baseline(bandno)%offset = -baselineAC(1:noChans,mindx) - &
                  DC_avg(1:noChans)
          ENDIF
          baselineDS%name = BaseName
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               Baseline(bandno)%offset, lastIndex=mindx, &
               disable_attrib=.TRUE.)

          IF (ANY (baselineAC(:,mindx) == FILLVALUE)) THEN
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
END MODULE SpectralBaseline
!=============================================================================
! $Log$
! Revision 2.2  2004/12/01 17:11:34  perun
! Add calculating and outputting DACS DC baseline
!
! Revision 2.1  2004/11/10 15:31:21  perun
! Initial commit
!
