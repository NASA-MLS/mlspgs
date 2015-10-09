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
MODULE MLSL1Common              ! Common data types for the MLSL1 program
!=============================================================================

  USE MLSCommon, ONLY: r4, r8, FileNameLen

  IMPLICIT NONE

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains data types that are common to the MLSL1 program.

  !---------------------------------------------------------------------------

  INTEGER, PRIVATE :: i

  ! Level 1 program type(s)

  CHARACTER (LEN=*), PARAMETER :: GHzType = "G"
  CHARACTER (LEN=*), PARAMETER :: THzType = "T"
  CHARACTER (LEN=*), PARAMETER :: LogType = "L"
  CHARACTER (LEN=1) :: L1ProgType = GHzType   !! Default type
  
  ! The science and engineering Level 0 file info

  TYPE L0FileInfo_T
     INTEGER :: sci_unit(2), sci_pcf(2)
     INTEGER :: eng_unit(6), eng_pcf(6)
     CHARACTER (LEN=FileNameLen) :: SciFileName(2)
     CHARACTER (LEN=FileNameLen) :: EngFileName(6)
  END TYPE L0FileInfo_T

  TYPE (L0FileInfo_T) :: L0FileInfo

  ! Information on the L1B data files in use

  TYPE L1BFileInfo_T
     INTEGER :: OAId            ! The HDF ID (handle) for the L1BOA file
     INTEGER :: RADDId, RADGId, RADTId  ! Ids for the L1BRAD files
     INTEGER :: EngId           ! The ID (non-HDF) (handle) for the L1BENG file
     INTEGER :: DiagId          ! The HDF ID (handle) for the L1BDIAG file
     INTEGER :: DiagTId         ! The HDF ID (handle) for the L1BDIAG THz file
     INTEGER :: LogId           ! The ID (non-HDF) (handle) for the L1BLog file
     INTEGER :: EngMAF_unit, SciMAF_unit  ! units for EngMAF & SciMAF files
     INTEGER :: MAF_data_unit   ! unit for Eng/Sci MAF data merged
     CHARACTER (LEN=FileNameLen) :: OAFileName  ! L1BOA file name
     CHARACTER (LEN=FileNameLen) :: RADDFileName, RADGFilename, RADTFileName
     CHARACTER (LEN=FileNameLen) :: EngFileName   ! L1BENG file name
     CHARACTER (LEN=FileNameLen) :: DiagFileName  ! L1BDIAG file name
     CHARACTER (LEN=FileNameLen) :: DiagTFileName ! L1BDIAG THz file name
     CHARACTER (LEN=FileNameLen) :: LogFileName   ! L1BLOG file name
  END TYPE L1BFileInfo_T
  TYPE (L1BFileInfo_T) :: L1BFileInfo

  ! Spectrometer info

  INTEGER, PARAMETER :: FBnum = 19
  INTEGER, PARAMETER :: GHzNum = 14
  INTEGER, PARAMETER :: FBchans = 25
  INTEGER, PARAMETER :: MBnum = 5
  INTEGER, PARAMETER :: MBchans = 11
  INTEGER, PARAMETER :: DACSnum = 4
  INTEGER, PARAMETER :: DACSchans = 129
  INTEGER, PARAMETER :: WFnum = 3
  INTEGER, PARAMETER :: WFchans = 4
  INTEGER, PARAMETER :: THzNum = 6
  INTEGER, PARAMETER :: THzChans = 25
  INTEGER, PARAMETER :: NumBands = 34

  ! Maximum number of MIFs per MAF

  INTEGER, PARAMETER :: MaxMIFs = 150

  ! Maximum number of MIFs per Limb view

  INTEGER, PARAMETER :: LimbMIFs = 125

  ! Maximum number of Limb altitudes for "S"pace view

  INTEGER, PARAMETER :: MaxAlts = 10

  ! MAF information data type (start time, integration time, etc.)

  TYPE MAFinfo_T
     REAL(r8) :: startTAI    ! TAI93 format
     REAL(r4) :: MIF_dur     ! MIF duration time (secs)
     REAL(r4) :: integTime   ! integration time (secs)
     INTEGER :: MIFsPerMAF   ! Number of MIFs per MAF
  END TYPE MAFinfo_T
  TYPE (MAFinfo_T) :: MAFinfo  ! Needed for L1BOA output

  !! Bank Logical type

  TYPE BankLogical_T
     LOGICAL :: FB(FBnum)          ! standard filter banks
     LOGICAL :: MB(MBnum)          ! mid-band filter banks
     LOGICAL :: WF(WFnum)          ! wide filters
     LOGICAL :: DACS(DACSnum)      ! DACS filters
  END TYPE BankLogical_T

  !! Bank Integer type

  TYPE BankInt_T
     INTEGER :: FB(FBnum)          ! standard filter banks
     INTEGER :: MB(MBnum)          ! mid-band filter banks
     INTEGER :: WF(WFnum)          ! wide filters
     INTEGER :: DACS(DACSnum)      ! DACS filters
  END TYPE BankInt_T

  ! Real types for all channels:

  TYPE Chan_R_T
     REAL :: FB(FBchans,FBnum)         ! standard filter banks
     REAL :: MB(MBchans,MBnum)         ! mid-band filter banks
     REAL :: WF(WFchans,WFnum)         ! wide filters
     REAL :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE Chan_R_T

  ! Logical type for all channels

  TYPE ChanLogical_T
     LOGICAL :: FB(FBchans,FBnum)         ! standard filter banks
     LOGICAL :: MB(MBchans,MBnum)         ! mid-band filter banks
     LOGICAL :: WF(WFchans,WFnum)         ! wide filters
     LOGICAL :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE ChanLogical_T

  ! Integer type for all channels

  TYPE ChanInt_T
     INTEGER :: FB(FBchans,FBnum)         ! standard filter banks
     INTEGER :: MB(MBchans,MBnum)         ! mid-band filter banks
     INTEGER :: WF(WFchans,WFnum)         ! wide filters
     INTEGER :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE ChanInt_T

  TYPE Chan_R8_T
     REAL(r8) :: FB(FBchans,FBnum)         ! standard filter banks
     REAL(r8) :: MB(MBchans,MBnum)         ! mid-band filter banks
     REAL(r8) :: WF(WFchans,WFnum)         ! wide filters
     REAL(r8) :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE Chan_R8_T

  ! Real types for all channels and bands, except DACS:

  TYPE ChanBand_R_T
     REAL :: FB(FBchans,FBnum+2)       ! standard filter banks/bands
     REAL :: MB(MBchans,MBnum)         ! mid-band filter bands
     REAL :: WF(WFchans,WFnum)         ! wide filter bands
  END TYPE ChanBand_R_T

  ! Band info:

  INTEGER, PARAMETER :: BandChans(NumBands) = (/ &
       (FBchans,i=1,21), (DACSchans,i=1,5),  (MBchans,i=1,5), (WFchans,i=1,3) /)

  !! Band Logical type

  TYPE BandChanR_T
     REAL :: SIGN(NumBands,DACSchans) = 1.0  ! "Good" precision sign
     INTEGER :: MaxChan(NumBands)    ! Maximum number of channels in band
  END TYPE BandChanR_T
  TYPE (BandChanR_T), SAVE :: BandChanBad = BandChanR_T (1.0, BandChans)

  !! Band Chi2:

  REAL :: BandChi2(NumBands,DACSchans) = -1.0

!! Initial OA counter MAF and index to be used to line up output files:

  INTEGER, DIMENSION(:), ALLOCATABLE :: OA_counterMAF
  INTEGER :: OA_counterIndex = 0     ! Nothing yet

!! Spacecraft YPR used to qualify data pointing:

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: SC_YPR

!! THz Geod Altitude use to qualify "space" views:

  REAL(r8), DIMENSION(:,:), ALLOCATABLE :: THz_GeodAlt

!! Factor to convert 23 bit encoder angle data (shifted by 1 bit) into degrees:

  REAL, PARAMETER :: DEG24 = 360.0 * 2.0**(-24)

!! Switching mirror position ranges (L, S, T1, T2)
!! 'L'imb, 'S'pace, 'T1' = 'primary' target, 'T2' = 'secondary' target
  TYPE SwMir_Range_T
     CHARACTER(len=1) :: pos
     REAL :: low_angle, high_angle
  END TYPE SwMir_Range_T

!! Switching mirror ranges

  REAL(r8), PARAMETER :: GHzTol = 0.05  ! GHz tolerance
  REAL(r8), PARAMETER :: THzTol = 360.0 / 16384.0  ! THz tolerance

  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range_A(4) = (/ & ! "A" side
       SwMir_Range_T ("L", 149.5-GHzTol, 149.5+GHzTol), &
       SwMir_Range_T ("S", 329.5-GHzTol, 329.5+GHzTol), &
       SwMir_Range_T ("T", 239.5-GHzTol, 239.5+GHzTol), &     ! Primary target
       SwMir_Range_T ("t", 59.5-GHzTol, 59.5+GHzTol) /)       ! Secondary target
  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range_B(4) = (/ & ! "B" side
       SwMir_Range_T ("L", 329.599-GHzTol, 329.599+GHzTol), &
       SwMir_Range_T ("S", 149.599-GHzTol, 149.599+GHzTol), &
       SwMir_Range_T ("T", 59.599-GHzTol, 59.599+GHzTol), &   ! Primary target
       SwMir_Range_T ("t", 239.599-GHzTol, 239.599+GHzTol) /) ! Secondary target
  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range_B_2(4) = (/ & ! "B" side 
       SwMir_Range_T ("L", 329.404-GHzTol, 329.404+GHzTol), & ! minus 0.195
       SwMir_Range_T ("S", 149.404-GHzTol, 149.404+GHzTol), &
       SwMir_Range_T ("T", 59.404-GHzTol, 59.404+GHzTol), &   ! Primary target
       SwMir_Range_T ("t", 239.404-GHzTol, 239.404+GHzTol) /) ! Secondary target
  TYPE (SwMir_Range_T), TARGET :: THz_SwMir_Range(3) = (/ &
       SwMir_Range_T ("L", 357.0, 360.0), &
       SwMir_Range_T ("T", 178.75-THzTol, 178.75+THzTol), &
       SwMir_Range_T ("S", 356.8-THzTol, 356.8+THzTol) /)
  TYPE (SwMir_Range_T), TARGET :: Discard_SwMir_Range(1) = (/ &
       SwMir_Range_T ("D", 0.0, 360.0) /)            ! no good pointing

!! Band switches and associated filter bank numbers

  INTEGER :: BandSwitch(5) = (/ 0, 0, 0, 0, 0 /)

  !! <whd> SwitchBank lists the number that output of each switch goes to, as a
  !! function of switch number ...
  !!
  !! SwitchBank(1) => DACS-1
  !! SwitchBank(2) => FB25-3
  !! SwitchBank(3) => FB25-8
  !! SwitchBank(4) => FB25-12
  !! SwitchBank(5) => FB25-15
  !! 
  !! What's on the input side of that switch is another matter
  !! </whd>
  INTEGER, PARAMETER :: SwitchBank(5) = (/ 1, 3, 8, 12, 15 /)

!! Bandwidths for all channels (in Hz?)

  REAL, PARAMETER, PRIVATE :: FB_BW(FBchans) = (/ &
       110.0, 110.0, 110.0, 74.0, 74.0, 74.0, 55.0, 37.0, 28.0, 18.4, 14.0, &
       9.2, 7.2, 9.2, 14.0, 18.4, 28.0, 37.0, 55.0, 74.0, 74.0, 74.0, 110.0, &
       110.0, 110.0 /) * 1.0e06
  TYPE (Chan_R_T), PARAMETER :: BandWidth = Chan_R_T ( &
       RESHAPE ( (/ (FB_BW, i=1,FBnum) /), (/ FBchans, FBnum /) ), &
       RESHAPE ( (/ (FB_BW(8:18), i=1,MBnum) /), (/ MBchans, MBnum /) ), & 
       RESHAPE ( (/ (500.0e06, i=1,WFchans*WFnum) /), (/ WFchans, WFnum /) ), &
       RESHAPE ( (/ (0.097e06, i=1,DACSchans*DACSnum) /), &
       (/ DACSchans, DACSnum /) ) )

!! DACS constants

  TYPE DACS_const_T
     INTEGER :: L
     INTEGER :: A(128)
  END TYPE DACS_const_T
  TYPE (DACS_const_T) :: DACS_const

!! Useful channel defaults:

  TYPE (Chan_R_T), TARGET :: deflt_gain, deflt_zero
  TYPE (ChanBand_R_T), TARGET :: deflt_chi2

  REAL :: tau = 1.0 / 6.0  ! default data rate (secs)

!! Useful constants:

  REAL, PARAMETER :: absZero_C = -273.16
  REAL(r8), PARAMETER :: planck = 6.62606891d-34
  REAL(r8), PARAMETER :: boltz = 1.380658d-23

!! First LO frequency (Hz) for each radiometer (except for R1 A/B, which is
!!  center freq)

  REAL, PARAMETER :: LO1(0:5) = (/ 118.7530e9, 118.7530e9, 191.9000e9, &
       239.6600e9, 642.8700e9, 2.5227816e12 /)

! HDF version

  INTEGER, PARAMETER :: HDFversion = 5   ! from now on!

  ! --------------------------------------------------------------------------

CONTAINS
!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE MLSL1Common
!=============================================================================

! $Log$
! Revision 2.19  2011/06/07 18:57:36  perun
! Add another Side B switching mirror table uploaded on 2011 DOY 153 at 19:10.
!
! Revision 2.18  2006/03/24 15:11:50  perun
! Add ChanLogical_T, ChanInt_T, MaxAlts, SC_YPR and THz_GeodAlt
!
! Revision 2.17  2005/12/14 17:00:14  perun
! Reset GHz switching mirror tolerance to original value
!
! Revision 2.16  2005/12/06 19:25:36  perun
! Removed BrightObjects_T
!
! Revision 2.15  2005/07/19 16:35:45  perun
! Increased GHz switching mirror tolerance per REC
!
! Revision 2.14  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.13  2005/01/28 17:02:13  perun
! Dimensioned BrightObjects type entries for GHz and THz
!
! Revision 2.12  2004/11/10 15:43:12  perun
! Add BandChans, BandChanBad, OA_counterMAF, OA_counterIndex, BrightObject_T;
! new value for DACS bandwidth
!
! Revision 2.11  2004/08/12 13:51:50  perun
! Version 1.44 commit
!
! Revision 2.10  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.9  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.8  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.7  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.6  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.5  2002/11/14 17:03:17  perun
! Constants for THz and DACS processing
!
! Revision 2.4  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.3  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
