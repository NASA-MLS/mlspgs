! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE BandTbls   ! Tables for all bands
!=============================================================================

  USE MLSL1Common, ONLY: FBchans, MBchans, WFchans, R4, R8, DACSchans, NumBands

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Load_Band_Tbls, LoadSidebandFracs, LoadSpilloverLoss
  PUBLIC :: BandLowerUpper_T, SideBandFrac, SpilloverLoss_T, SpilloverLoss, &
       RadiometerLoss_T, RadiometerLoss, BandFreq

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE BandLowerUpper_T
     REAL(r4), DIMENSION(:), POINTER :: lower, upper
  END TYPE BandLowerUpper_T

  TYPE (BandLowerUpper_T) :: SideBandFrac(NumBands), BandFreq(NumBands)

  TYPE SpilloverLoss_T
     REAL(r4), DIMENSION(:,:), POINTER :: lower, upper
  END TYPE SpilloverLoss_T

  TYPE (SpilloverLoss_T) :: SpilloverLoss(NumBands)

  TYPE RadiometerLoss_T
     CHARACTER(len=3) :: Name
     REAL(r4) :: Ohmic, Spillover, Radiance
  END TYPE RadiometerLoss_T

  TYPE (RadiometerLoss_T), PARAMETER :: RadiometerLoss(0:4) = (/ &
       RadiometerLoss_T ("R1A", 0.9982**3, 0.9987, 150.0), &  ! Index '0'
       RadiometerLoss_T ("R1B", 0.9982**3, 0.9964, 150.0), &
       RadiometerLoss_T ("R2 ", 0.99736**3, 0.9992, 150.0), &
       RadiometerLoss_T ("R3 ", 0.99449**3, 0.9987, 150.0), &
       RadiometerLoss_T ("R4 ", 0.98819**3, 0.9916, 150.0)  /)

CONTAINS

!=============================================================================
  SUBROUTINE Load_Band_tbls
!=============================================================================

    INTEGER :: i

!! Allocate Sideband fraction and Band frequency arrays:

    DO i = 1, 21
       ALLOCATE (SideBandFrac(i)%lower(FBchans))
       ALLOCATE (SideBandFrac(i)%upper(FBchans))
       ALLOCATE (BandFreq(i)%lower(FBchans))
       ALLOCATE (BandFreq(i)%upper(FBchans))
    ENDDO
    DO i = 22, 26
       ALLOCATE (SideBandFrac(i)%lower(DACSchans))
       ALLOCATE (SideBandFrac(i)%upper(DACSchans))
       ALLOCATE (BandFreq(i)%lower(DACSchans))
       ALLOCATE (BandFreq(i)%upper(DACSchans))
    ENDDO
    DO i = 27, 31
       ALLOCATE (SideBandFrac(i)%lower(MBchans))
       ALLOCATE (SideBandFrac(i)%upper(MBchans))
       ALLOCATE (BandFreq(i)%lower(MBchans))
       ALLOCATE (BandFreq(i)%upper(MBchans))
    ENDDO
    DO i = 32, 34
       ALLOCATE (SideBandFrac(i)%lower(WFchans))
       ALLOCATE (SideBandFrac(i)%upper(WFchans))
       ALLOCATE (BandFreq(i)%lower(WFchans))
       ALLOCATE (BandFreq(i)%upper(WFchans))
    ENDDO

!! Allocate and Initialize Spillover Loss

    DO i = 1, 31
       ALLOCATE (SpilloverLoss(i)%lower(3,1))
       ALLOCATE (SpilloverLoss(i)%upper(3,1))
       SpilloverLoss(i)%lower = 0.0
       SpilloverLoss(i)%upper = 0.0
    ENDDO

    DO i = 32, 34
       ALLOCATE (SpilloverLoss(i)%lower(3,4))
       ALLOCATE (SpilloverLoss(i)%upper(3,4))
       SpilloverLoss(i)%lower = 0.0
       SpilloverLoss(i)%upper = 0.0
    ENDDO

    DO i = 1, 34
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

    DO i = 1, 34
       READ (unit, '(A)') line
       READ (unit, *) SideBandFrac(i)%lower
       READ (unit, '(A)') line
       READ (unit, *) SideBandFrac(i)%upper
    ENDDO

  END SUBROUTINE LoadSidebandFracs

!=============================================================================
  SUBROUTINE LoadSpilloverLoss (unit)
!=============================================================================

    INTEGER :: unit

    CHARACTER (len=80) :: line
    INTEGER :: bandno, chan, ios
    REAL :: h(3)

    READ (unit, '(A)') line

    DO
       READ (unit, '(A)', iostat=ios) line
       IF (ios /= 0) EXIT
       chan = 1
       READ (line(1:2), *) bandno
       READ (line(4:4), *, iostat=ios) chan
       READ (line(25:), *) h
       IF (INDEX (line, "L") /= 0) THEN
          SpilloverLoss(bandno)%lower(:,chan) = h
       ELSE
          SpilloverLoss(bandno)%upper(:,chan) = h
       ENDIF
    ENDDO

  END SUBROUTINE LoadSpilloverLoss

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

END MODULE BandTbls

! $Log$
! Revision 2.2  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.1  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
