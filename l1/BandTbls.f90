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
MODULE BandTbls   ! Tables for all bands
!=============================================================================

  USE MLSL1Common, ONLY: FBchans, MBchans, WFchans, R4, R8, DACSchans, &
       MaxAlts, NumBands

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Load_Band_Tbls, LoadSidebandFracs, LoadSpilloverLoss, &
       LoadBandAlts, LoadDefltChi2, GetEta_TSL, LoadFourierCoeffs, &
       GetDeltaRads
  PUBLIC :: BandLowerUpper_T, SideBandFrac, SpilloverLoss_T, SpilloverLoss, &
     RadiometerLoss_T, RadiometerLoss, BandFreq, BandAlt_T, BandAlt, nAlts, &
     minAlt,delrad_1_31, delrad_32_34

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE BandLowerUpper_T
     REAL(r4), DIMENSION(:), POINTER :: lower, upper
  END TYPE BandLowerUpper_T

  TYPE (BandLowerUpper_T) :: SideBandFrac(NumBands), BandFreq(NumBands)

  TYPE SpilloverLoss_T
     REAL(r4), DIMENSION(:,:), POINTER :: lower, upper, eta_TSL
  END TYPE SpilloverLoss_T

  TYPE (SpilloverLoss_T) :: SpilloverLoss(NumBands)

  TYPE RadiometerLoss_T
     CHARACTER(len=3) :: Name
     REAL(r4) :: Ohmic, Spillover, Radiance
  END TYPE RadiometerLoss_T

  TYPE (RadiometerLoss_T), PARAMETER :: RadiometerLoss(0:4) = (/ &
       RadiometerLoss_T ("R1A", 0.9938, 0.9987, 150.0), &  ! Index '0'
       RadiometerLoss_T ("R1B", 0.9795, 0.9964, 150.0), &
       RadiometerLoss_T ("R2 ", 0.9914, 0.9992, 150.0), &
       RadiometerLoss_T ("R3 ", 0.9746, 0.9987, 150.0), &
       RadiometerLoss_T ("R4 ", 0.9802, 0.9916, 150.0)  /)

  ! Defined as POINTERs to get around the restriction (since removed from
  ! fortran standard) of components of user types not being ALLOCATABLE. I don't
  ! believe they are every assigned to point to anything, and if they were it
  ! could be a memory leak. Allocated in Load_Bnd_Tbls
  !
  TYPE BandAlt_T
     REAL(r4), DIMENSION(:), POINTER :: Meters
     INTEGER, DIMENSION(:), POINTER :: indx
  END TYPE BandAlt_T

  TYPE (BandAlt_T) :: BandAlt(NumBands)

  TYPE Fourier_T
     REAL(r4), DIMENSION(:,:), POINTER :: Coeff
  END TYPE Fourier_T

  TYPE (Fourier_T), SAVE :: aSA(NumBands), bSA(NumBands)

  INTEGER :: nAlts
  REAL :: MinAlt(MaxAlts), delrad_1_31(31), delrad_32_34(32:34,4)

CONTAINS

!=============================================================================
  SUBROUTINE Load_Band_tbls
!=============================================================================

    INTEGER :: i

!! Allocate Sideband fraction and Band frequency and altitude arrays:

    DO i = 1, 21
       ALLOCATE (SideBandFrac(i)%lower(FBchans))
       ALLOCATE (SideBandFrac(i)%upper(FBchans))
       ALLOCATE (BandFreq(i)%lower(FBchans))
       ALLOCATE (BandFreq(i)%upper(FBchans))
       ALLOCATE (BandAlt(i)%Meters(FBchans))
       ALLOCATE (BandAlt(i)%indx(FBchans))
    ENDDO
    DO i = 22, 26
       ALLOCATE (SideBandFrac(i)%lower(DACSchans))
       ALLOCATE (SideBandFrac(i)%upper(DACSchans))
       ALLOCATE (BandFreq(i)%lower(DACSchans))
       ALLOCATE (BandFreq(i)%upper(DACSchans))
       ALLOCATE (BandAlt(i)%Meters(DACSchans))
       ALLOCATE (BandAlt(i)%indx(DACSchans))
    ENDDO
    DO i = 27, 31
       ALLOCATE (SideBandFrac(i)%lower(MBchans))
       ALLOCATE (SideBandFrac(i)%upper(MBchans))
       ALLOCATE (BandFreq(i)%lower(MBchans))
       ALLOCATE (BandFreq(i)%upper(MBchans))
       ALLOCATE (BandAlt(i)%Meters(MBchans))
       ALLOCATE (BandAlt(i)%indx(MBchans))
    ENDDO
    DO i = 32, 34
       ALLOCATE (SideBandFrac(i)%lower(WFchans))
       ALLOCATE (SideBandFrac(i)%upper(WFchans))
       ALLOCATE (BandFreq(i)%lower(WFchans))
       ALLOCATE (BandFreq(i)%upper(WFchans))
       ALLOCATE (BandAlt(i)%Meters(WFchans))
       ALLOCATE (BandAlt(i)%indx(WFchans))
    ENDDO

!! Allocate and Initialize Spillover Loss and allocate Fourier Coeffs

    DO i = 1, 31
       ALLOCATE (SpilloverLoss(i)%lower(3,1))
       ALLOCATE (SpilloverLoss(i)%upper(3,1))
       ALLOCATE (SpilloverLoss(i)%eta_TSL(3,1))
       SpilloverLoss(i)%lower = 0.0
       SpilloverLoss(i)%upper = 0.0
       SpilloverLoss(i)%eta_TSL = 0.0
       ALLOCATE (aSA(i)%Coeff(4,1))
       ALLOCATE (bSA(i)%Coeff(4,1))
    ENDDO

    DO i = 32, 34
       ALLOCATE (SpilloverLoss(i)%lower(3,4))
       ALLOCATE (SpilloverLoss(i)%upper(3,4))
       ALLOCATE (SpilloverLoss(i)%eta_TSL(3,4))
       SpilloverLoss(i)%lower = 0.0
       SpilloverLoss(i)%upper = 0.0
       SpilloverLoss(i)%eta_TSL = 0.0
       ALLOCATE (aSA(i)%Coeff(4,4))
       ALLOCATE (bSA(i)%Coeff(4,4))
    ENDDO

    DO i = 1, NumBands
       CALL BandFreqs (i, BandFreq(i)%lower, BandFreq(i)%upper)
    ENDDO

  END SUBROUTINE Load_Band_tbls

!=============================================================================
  SUBROUTINE LoadSidebandFracs (unit)
!=============================================================================

    INTEGER :: unit

    CHARACTER (len=80) :: line
    INTEGER :: i

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:1) /= ";") EXIT
    ENDDO

    DO i = 1, NumBands
       READ (unit, '(A)') line
       READ (unit, *) SideBandFrac(i)%lower
       READ (unit, '(A)') line
       READ (unit, *) SideBandFrac(i)%upper
    ENDDO

  END SUBROUTINE LoadSidebandFracs

!=============================================================================
  SUBROUTINE LoadBandAlts (unit)
!=============================================================================

    INTEGER :: unit

    CHARACTER (len=80) :: line
    INTEGER :: band, i, n
    REAL :: alt

    nAlts = 0; MinAlt = -1.0    ! No altitudes yet

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:1) /= ";") EXIT
    ENDDO

    ! <whd> This code reads the BandAlts.tbl file (PCF id=912). It stores only
    ! the minimum altitudes for each band and an index into the module variable
    ! `MinAlts'.  Later, in SortQualify::QualifyCurrentMAF this index is used to
    ! ???</whd>
    
    DO band = 1, NumBands
       READ (unit, '(A)') line
       READ (unit, *) BandAlt(band)%Meters
       BandAlt(band)%Meters = BandAlt(band)%Meters * 1.0e03   ! Input is in Km
       BandAlt(Band)%indx = 0
       DO n = 1, SIZE (BandAlt(band)%Meters)
          alt = BandAlt(band)%Meters(n)
          IF (alt > 0.0) THEN       ! Good altitude
             IF (nAlts > 0) THEN    ! Check if already there or add to table
                IF (COUNT(MinAlt(1:nAlts) == alt) == 0) THEN  ! New entry
                   nAlts = nAlts + 1
                   IF (nAlts > MaxAlts) THEN
                      PRINT *, 'Increase MaxAlts!'
                      STOP
                   ENDIF
                   MinAlt(nAlts) = alt
                ENDIF
             ELSE
                nAlts = 1
                MinAlt(nAlts) = alt
             ENDIF
             ! Mark the elements in BandAlt(band)%Meters == i-th altitude,
             ! sorted from min to max (whd: why?)
             DO i = 1, nAlts
                IF (alt == MinAlt(i)) THEN
                   BandAlt(band)%indx(n) = i
                   EXIT
                ENDIF
             ENDDO
          ELSE
             BandAlt(band)%Meters(n) = HUGE (1.0)   ! Huge value for tests
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE LoadBandAlts

!=============================================================================
  SUBROUTINE LoadSpilloverLoss (unit)
!=============================================================================

    INTEGER :: unit

    CHARACTER (len=80) :: line
    INTEGER :: bandno, chan, ios
    REAL :: h(3), eta(3)

    READ (unit, '(A)') line

    DO
       READ (unit, '(A)', iostat=ios) line
       IF (ios /= 0) EXIT
       chan = 1
       READ (line(1:2), *) bandno
       READ (line(4:4), *, iostat=ios) chan
       READ (line(25:), *) h, eta
       IF (INDEX (line, "L") /= 0) THEN
          SpilloverLoss(bandno)%lower(:,chan) = h
       ELSE
          SpilloverLoss(bandno)%upper(:,chan) = h
       ENDIF
       SpilloverLoss(bandno)%eta_TSL(:,chan) = eta
    ENDDO

  END SUBROUTINE LoadSpilloverLoss

!=============================================================================
  SUBROUTINE GetEta_TSL (bandno, channo, eta_TSL)
!=============================================================================

    INTEGER, INTENT (IN) :: bandno, channo
    REAL, INTENT (OUT) :: eta_TSL(3)

    IF (bandno < 32) THEN    ! <vp> All channels have same ets values <whd>
       !Looking at GHzReflSpillEffs.tbl, the file from where this information is
       !read, I see that it is false that 'all channels' for bandno < 32 have the
       !same ets value<whd>
       eta_TSL = SpilloverLoss(bandno)%eta_TSL(:,1)
    ELSE
       eta_TSL = SpilloverLoss(bandno)%eta_TSL(:,channo)
    ENDIF

  END SUBROUTINE

!=============================================================================
  SUBROUTINE BandFreqs (band, lsbf, usbf)
!=============================================================================

    INTEGER, INTENT (IN) :: band
    REAL(R4), INTENT (OUT), DIMENSION(:) :: lsbf, usbf

    INTEGER :: i, bindx

! Parameters


!! GHz

    REAL(R8), DIMENSION(14), PARAMETER :: LO1 = (/ &
         126.8D0, &                                          ! R1A/B first LO,
         191.9D0, 191.9D0, 191.9D0, 191.9D0, 191.9D0, &      ! R2 (bands 2..6),
         239.66D0, 239.66D0, 239.66D0, &                     ! R3 (bands 7..9),
         642.87D0, 642.87D0, 642.87D0, 642.87D0, 642.87D0 /) ! R4 (bands 10..14)
    REAL(R8), DIMENSION(14), PARAMETER :: LO2 = (/ &
         8.9467D0, &
         7.68571D0, 8.18D0, 11.20167D0, 13.35667D0, 13.33667D0, &
         4.8449D0, 4.8085D0, 10.0168D0, &                       ! Second LOs,
         5.69625D0, 8.224D0, 10.8786D0, 16.0369D0, 18.384D0 /)  ! by Band

! Channel offsets, from 'zero IF' frequency, MHz, for FB25 channels.

    REAL(R8), DIMENSION(25), PARAMETER :: foffset25 = (/ &
         -575.0D-3, -479.0D-3, -383.0D-3, -303.0D-3, -239.0D-3, -175.0D-3, &
         -119.0D-3, -79.0D-3, -51.0D-3, -31.0D-3, -17.0D-3, -7.0D-3, &
         0.0D-3,  7.0D-3, 17.0D-3, 31.0D-3, 51.0D-3, 79.0D-3, 119.0D-3, &
         175.0D-3, 239.0D-3, 303.0D-3, 383.0D-3, 479.0D-3, 575.0D-3 /)
    REAL(R8), DIMENSION(11), PARAMETER :: foffset11 = foffset25(8:18)

!! THz

! THz Zero IF Frequencies (15-17 (18-20), LSB/USB):

    REAL(R8), DIMENSION(2,3), PARAMETER :: THz_Z_IF = RESHAPE ( &
         (/ 2513.47350D0, 2530.28970D0, 2509.10570D0, 2534.65750D0, &
            2501.48040D0, 2542.28280D0 /), (/ 2, 3 /)) + 900.0D-03

    REAL(R8), DIMENSION(14), PARAMETER :: GHz_dir = (/ &
          1.0D0, &
         -1.0D0, -1.0D0, 1.0D0, 1.0D0, -1.0D0, & ! Lower sideband filter channel
          1.0D0, -1.0D0, 1.0D0, &                ! order in RF space (1.0D0
         -1.0D0, 1.0D0, 1.0D0, -1.0D0, 1.0D0 /)  ! signifies that chan_1_freq <
                                                 ! chan25_freq)

    REAL(R8), DIMENSION(3), PARAMETER :: THz_dir =  (/ 1.0D0, 1.0D0, -1.0D0 /)

!! DACS

    REAL(R8), DIMENSION(5), PARAMETER :: DACS_LO = (/ &
         126.80D0, 191.90D0, 239.660D0, 239.660D0, 126.80D0 /)
    REAL(R8), DIMENSION(5), PARAMETER :: DACS_IF = (/ &
       8.047916667D0, 8.584464286D0, 3.946250D0, 9.117916667D0, 8.047916667D0 /)
    REAL(R8), DIMENSION(129), PARAMETER :: DACS_offset = &
         ((/ (i, i=0,128) /) - 64.0D0) / 128.0 * 12.5D-3
    REAL(R8), DIMENSION(5), PARAMETER :: DACS_dir = (/ &
         -1.0D0, 1.0D0, -1.0D0, -1.0D0, -1.0D0 /)

    bindx = band

    IF (bindx == 21) bindx = 1  ! 21 and 1 have same frequencies
    IF (bindx < 15) THEN
       lsbf = LO1(bindx) - LO2(bindx) + (foffset25 +900.D-3) * GHz_dir(bindx)
       usbf = 0.0D0
       IF (bindx /= 1) usbf = LO1(bindx) + LO2(bindx) - &
            (foffset25 + 900.0D-3) * GHz_dir(bindx)
    ELSE IF (bindx >= 15 .AND. bindx <= 20) THEN
       IF (bindx > 17) bindx = bindx - 3
       bindx = bindx - 14
       lsbf = THz_Z_IF(1,bindx) + foffset25 * THz_dir(bindx)
       usbf = THz_Z_IF(2,bindx) - foffset25 * THz_dir(bindx)
    ELSE IF (bindx >= 22 .AND. bindx <= 26) THEN
       bindx = bindx - 21
       lsbf = DACS_LO(bindx) - DACS_IF(bindx) + &
            DACS_dir(bindx) * DACS_offset
       usbf = 0.0D0
       IF (band /= 22 .AND. band /= 26) &
            usbf = DACS_LO(bindx) + DACS_IF(bindx) - &
            DACS_dir(bindx) * DACS_offset
    ELSE IF (bindx >= 27 .AND. bindx <= 31) THEN
       IF (bindx == 27) THEN
          lsbf = 191.9D0 - 14.635D0 - foffset11
          usbf = 191.9D0 + 14.635D0 + foffset11
       ELSE IF (bindx == 28) THEN
          lsbf = 642.87D0 - 6.84635D0 + foffset11
          usbf = 642.87D0 + 6.84635D0 - foffset11
       ELSE  IF (bindx == 29) THEN
          lsbf = 642.87D0 - 6.98595D0 - foffset11
          usbf = 642.87D0 + 6.98595D0 + foffset11
       ELSE IF (bindx == 30) THEN
          lsbf = 642.87D0 - 17.6305D0 + foffset11
          usbf = 642.87D0 + 17.6305D0 - foffset11
       ELSE  IF (bindx == 31) THEN
          lsbf = 642.87D0 - 18.088D0 + foffset11
          usbf = 642.87D0 + 18.088D0 - foffset11
       ENDIF
    ELSE
       IF (bindx /= 33) THEN
          lsbf = (/ 122.0D0, 120.5D0, 117.0D0, 115.3D0 /)
          usbf = 0.0D0
       ELSE
          lsbf = (/ 236.65999999799999D0, 234.85999999800001D0, &
               232.46000000399999D0, 231.85999999200001D0 /)
          usbf = (/ 242.659999992D0, 244.45999999200001D0, &
               246.86000000799999D0, 247.459999998D0 /)
       ENDIF
    ENDIF

  END SUBROUTINE BandFreqs

!=============================================================================
  SUBROUTINE LoadDefltChi2 (unit, stat)
!=============================================================================

    USE MLSL1Common, ONLY: deflt_chi2, WFnum, BandChi2

    INTEGER :: unit, stat

    CHARACTER(LEN=80) :: line
    CHARACTER(LEN=2) :: chan_type
    INTEGER :: i, ios, bandindx
    REAL, POINTER, DIMENSION(:,:) :: chan_dat

    stat = 0

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:6) == "#DATA") EXIT
    ENDDO

! read data

    bandindx = 1   ! index for BandChi2 array to output to L1B RAD files
    i = 1
    chan_type = "FB"
    chan_dat => deflt_chi2%fb
    DO

       DO
          READ (unit, '(A)', IOSTAT=ios) line
          IF (line(1:1) == "#" .OR. ios /= 0) EXIT
       ENDDO

       IF (ios /= 0) EXIT

       READ (unit, *) chan_dat(:,i)
       BandChi2(bandindx,1:SIZE(chan_dat(:,i))) = chan_dat(:,i)
       i = i + 1
       bandindx = bandindx + 1

       SELECT CASE (chan_type)

       CASE ("FB")
          IF (I > SIZE (deflt_chi2%fb(1,:))) THEN
             i = 1
             bandindx = 27    ! Start of MB data
             chan_dat => deflt_chi2%mb
             chan_type = "MB"
          ENDIF
       CASE ("MB")
          IF (I > SIZE (deflt_chi2%mb(1,:))) THEN
             i = 1
             bandindx = 32    ! Start of WF data
             chan_dat => deflt_chi2%wf
             chan_type = "WF"
          ENDIF
       END SELECT

    ENDDO

    IF (chan_type /= "WF" .AND. i <= WFNum) stat = -1  ! Error!!!

  END SUBROUTINE LoadDefltChi2

!=============================================================================
  SUBROUTINE LoadFourierCoeffs (unit, asciiUTC)
!=============================================================================

    use Constants, ONLY: Rad2Deg, Deg2Rad
    USE MLSMessageModule, ONLY: MLSMESSAGE, MLSMSG_Info, MLSMSG_Error
    USE SDPToolkit, ONLY: spacecraftId, PGSd_SUN, PGS_S_SUCCESS

    INTEGER, INTENT(IN) :: unit
    CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

    CHARACTER (len=80) :: line
    INTEGER :: b, m, chan, returnStatus
    REAL(r8) :: sc_frame_vector(3), orb(3)
    REAL(r8) :: meanSolTimG, meanSolTimL, apparSolTimL, solRA, solDec
    REAL :: xSA, ySA, rSA, kSA, thSA, Beta, cd1, sd1
    INTEGER, PARAMETER :: nchans(32:34) = (/ 4, 1, 4 /)
    REAL(r8), PARAMETER :: offset = 0.0d0
    REAL, PARAMETER :: si = 0.3979486, b0 = 25.9252, dSA = 0.42958, &
         eSA = -7.94452

! Externals

    INTEGER, EXTERNAL :: Pgs_cbp_sat_cb_vector, Pgs_csc_scToORB
    INTEGER, EXTERNAL :: PGS_CBP_SolarTimeCoords

! Determine Solor Beta Angle at start of day:

    returnStatus = Pgs_cbp_sat_cb_vector (spacecraftId, 1, asciiUTC, &
         offset, PGSd_SUN, sc_frame_vector)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      print *, 'Check that your toolkit/database/linux64/CBP/de200.eos' // &
        &  ' is up-to-date'
      CALL MLSMessage ( MLSMSG_Info, ModuleName, 'Check that your ' // &
        &  'toolkit/database/linux64/CBP/de200.eos is up-to-date')
      print *, 'Pgs_cbp_sat_cb_vector' &
      & //' at '//asciiUTC
      CALL MLSMessage ( MLSMSG_Error, ModuleName, 'Pgs_cbp_sat_cb_vector' &
      & //' at '//asciiUTC)
    endif

    returnStatus = Pgs_csc_scToORB (spacecraftId, 1, asciiUTC, &
         offset, sc_frame_vector, orb)
    if ( returnStatus /= PGS_S_SUCCESS ) &
      & CALL MLSMessage ( MLSMSG_Error, ModuleName, 'Pgs_csc_scToORB' &
      & //' at '//asciiUTC)

    returnStatus = PGS_CBP_SolarTimeCoords (asciiUTC, 0.0d0, meanSolTimG, &
         meanSolTimL, apparSolTimL, solRA, solDec)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      print *, 'Check that your toolkit/database/linux64/CBP/de200.eos' // &
        &  ' is up-to-date'
      CALL MLSMessage ( MLSMSG_Info, ModuleName, 'Check that your ' // &
        &  'toolkit/database/linux64/CBP/de200.eos is up-to-date')
      print *, 'PGS_CBP_SolarTimeCoords' &
      & //' at '//asciiUTC
      CALL MLSMessage ( MLSMSG_Error, ModuleName, 'PGS_CBP_SolarTimeCoords' &
      & //' at '//asciiUTC)
    endif


    sd1 = sin (solDec) / si
    cd1 = sqrt ((1.0 - sd1*sd1))
    if (abs (solRA*Rad2Deg - 180.0) < 90.0) cd1 = -cd1
    Beta = b0 + dSA * cd1 + eSA * sd1
    Beta = Beta * Deg2Rad                 ! back to radians

!    Beta = ASIN (-orb(2) / SQRT (orb(1)**2 + orb(2)**2 + orb(3)**2))

! Read first line as comment:

    READ (unit, '(A)') line

    DO b = 1, 14
       DO m = 1, 4

          READ (unit, '(A)') line
          READ (line(10:), *) xSA, ySA, rSA, kSA, thSA

          aSA(b)%coeff(m,1) = rSA * (xSA + COS (thSA + m * kSA * Beta))
          bSA(b)%coeff(m,1) = rSA * (ySA + SIN (thSA + m * kSA * Beta))
       ENDDO
    ENDDO

    DO b = 21, 31
       DO m = 1, 4

          READ (unit, '(A)') line
          READ (line(10:), *) xSA, ySA, rSA, kSA, thSA

          aSA(b)%coeff(m,1) = rSA * (xSA + COS (thSA + m * kSA * Beta))
          bSA(b)%coeff(m,1) = rSA * (ySA + SIN (thSA + m * kSA * Beta))
       ENDDO
    ENDDO

    DO b = 32, 34
       DO chan = 1, nchans(b)
          DO m = 1, 4

             READ (unit, '(A)') line
             READ (line(10:), *) xSA, ySA, rSA, kSA, thSA

             aSA(b)%coeff(m,chan) = rSA * (xSA + COS (thSA + m * kSA * Beta))
             bSA(b)%coeff(m,chan) = rSA * (ySA + SIN (thSA + m * kSA * Beta))
             IF (b == 33) THEN   ! since band 33 has only 1 channel
                aSA(b)%coeff(m,2:) = aSA(b)%coeff(m,1)
                bSA(b)%coeff(m,2:) = bSA(b)%coeff(m,1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE LoadFourierCoeffs

!=============================================================================
  FUNCTION DeltaRadiance (phi, band, chan) RESULT (delrad)
!=============================================================================

    INTEGER, INTENT (IN) :: band, chan
    REAL, INTENT (IN) :: phi               ! in radians

    REAL :: delrad

    INTEGER :: channo, m

    IF (band < 32) THEN
       channo = 1         ! only 1 channel needed
    ELSE
       channo = chan      ! 1 of 4 wide band channels
    ENDIF

    delrad = 0.0
    DO m = 1, 4
       delrad = delrad + aSA(band)%coeff(m,channo) * COS(m*phi) + &
            bSA(band)%coeff(m,channo) * SIN(m*phi)
    ENDDO

  END FUNCTION DeltaRadiance

!=============================================================================
  SUBROUTINE GetDeltaRads (phi)
!=============================================================================

    REAL, INTENT (IN) :: phi               ! in radians

    INTEGER :: b, chan

    DO b = 1, 14
       delrad_1_31(b) = DeltaRadiance (phi, b, 1)
    ENDDO

    DO b = 21, 31
       delrad_1_31(b) = DeltaRadiance (phi, b, 1)
    ENDDO

    DO b = 32, 34
       DO chan = 1, 4
          delrad_32_34(b,chan) = DeltaRadiance (phi, b, chan)
       ENDDO
    ENDDO

  END SUBROUTINE GetDeltaRads

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE BandTbls

! $Log$
! Revision 2.14  2023/06/06 22:37:15  pwagner
! Make helpful message refer to de200.eos specifically
!
! Revision 2.13  2023/05/25 22:24:20  pwagner
! Trying to make level 1 crash if de200.eos is out-of-date
!
! Revision 2.12  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.11.6.1  2015/10/09 10:21:37  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.11  2009/06/01 13:59:20  perun
! Remove extraneous debug print statement.
!
! Revision 2.10  2009/05/13 20:33:05  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.9  2006/08/02 18:51:45  perun
! Corrected nested loops in reading wide filters
!
! Revision 2.8  2006/06/14 13:43:47  perun
! Add LoadFourierCoeffs and GetDeltaRads routines
!
! Revision 2.7  2006/03/24 15:06:25  perun
! Add Band Altitudes table and routines
!
! Revision 2.6  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2004/12/01 17:08:17  perun
! Update ohmic loss values
!
! Revision 2.4  2004/11/10 15:41:58  perun
! Change radiometer ohmic loss value for R4 (per WGR); correct indexes to read MB
!  and WF default chi2
!
! Revision 2.3  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.2  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.1  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
