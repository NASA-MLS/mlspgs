! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L0_sci_tbls  ! Define L0 science tables
!=============================================================================

  USE MLSL1Common, ONLY: R8, FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
       DACSnum, DACSchans, MaxMIFs, THzChans, THzNum, NumBands, BankLogical_T

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

! Level 0 Raw Science Packet #1 (sized for 1 MIF)

  CHARACTER (LEN=1024), TARGET :: l0_sci1
  CHARACTER (LEN=1), POINTER :: type_sci1

! Level 0 Raw Science Packet #2 (sized for 1 MIF)

  CHARACTER (LEN=1024), TARGET :: l0_sci2
  CHARACTER (LEN=1), POINTER :: type_sci2

!! Science packet data type

  INTEGER sci_type

!! CCSDS header

  TYPE CCSDS_T
     SEQUENCE
     CHARACTER (LEN=2) :: ID
     CHARACTER (LEN=2) :: seq_cntl
     CHARACTER (LEN=2) :: len
     CHARACTER (LEN=1) :: hdr
     CHARACTER (LEN=1) :: octet_1
     CHARACTER (LEN=1) :: octet_2
     CHARACTER (LEN=4) :: coarse_time
     CHARACTER (LEN=2) :: fine_time
  END TYPE CCSDS_T

!! Using format of 08/25/03

  TYPE Sci1_T1_T     ! science packet 1 Type I
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: TYPE
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=1) :: spare1
!
!! FB_dat contains the filter banks, in this order:
!!  FB12, FB01, status, FB04, FB02, status, FB07, FB03, FB05, FB13, MB01, FB06,
!!  FB09, FB08
!
     CHARACTER (LEN=2) :: FB_dat(2*25+1+2*25+1+4*25+11+25+1+2*25)
      ! filter bank plus statuses
     CHARACTER (LEN=1) :: GM05_stat
     CHARACTER (LEN=1) :: GM06_stat
     CHARACTER (LEN=1) :: GM07_stat(2)
     CHARACTER (LEN=1) :: GM08_stat(2)
     CHARACTER (LEN=1) :: GM09_stat(2)
     CHARACTER (LEN=13) :: GM10_stat
     CHARACTER (LEN=13) :: GM11_stat
     CHARACTER (LEN=1) :: GM12_stat
     CHARACTER (LEN=13) :: GM13_stat
     CHARACTER (LEN=1) :: GM14_stat
     CHARACTER (LEN=1) :: GM15_stat(5)
     CHARACTER (LEN=1) :: DACS(320)
     CHARACTER (LEN=1) :: diag_addr(32)
     CHARACTER (LEN=1) :: read_value
     CHARACTER (LEN=1) :: status_value
     CHARACTER (LEN=2) :: status_version
     CHARACTER (LEN=1) :: read_riu_id
     CHARACTER (LEN=1) :: status_riu_id
     CHARACTER (LEN=1) :: spare2(10)
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci1_T1_T

  TYPE Sci2_T1_T     ! science packet 2 Type I
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: TYPE
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=1) :: spare1
!
!! FB_dat contains the filter banks, in this order:
!!  MB03, FB10, status, MB02, FB11, status, FB15, FB14, MB05, MB04, status,
!!   FB17, FB16, FB19, FB18
!
     CHARACTER (LEN=2) :: FB_dat(11+25+1+11+25+1+2*25+2*11+1+4*25)
                          ! filter, mid-band plus status data
     CHARACTER (LEN=1) :: DACS(320)
     CHARACTER (LEN=24) :: thz_scan_switch
     CHARACTER (LEN=80) :: laser_lo
     CHARACTER (LEN=1) :: TM03_stat(2)
     CHARACTER (LEN=24) :: ghz_ant_scan
     CHARACTER (LEN=24) :: ghz_switch
     CHARACTER (LEN=1) :: spare2(25)
     CHARACTER (LEN=1) :: attenuation(5)
     CHARACTER (LEN=1) :: spare3(2)
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci2_T1_T

  TYPE Sci1_T2_T     ! science packet 1 Types II and III
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: TYPE
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=1) :: spare1
!
!! FB_dat contains the filter banks, in this order:
!!  FB12, FB01, status, FB04, FB02, status, FB07, FB03, FB05, FB13, FB09, FB08
!
     CHARACTER (LEN=2) :: FB_dat(2*25+1+2*25+1+6*25) ! filter bank plus statuses
     CHARACTER (LEN=1) :: GM05_stat
     CHARACTER (LEN=1) :: GM06_stat
     CHARACTER (LEN=1) :: GM07_stat(2)
     CHARACTER (LEN=1) :: GM08_stat(2)
     CHARACTER (LEN=1) :: GM09_stat(2)
     CHARACTER (LEN=13) :: GM10_stat
     CHARACTER (LEN=13) :: GM11_stat
     CHARACTER (LEN=1) :: GM12_stat
     CHARACTER (LEN=13) :: GM13_stat
     CHARACTER (LEN=1) :: GM14_stat
     CHARACTER (LEN=1) :: GM15_stat(5)
     CHARACTER (LEN=1) :: DACS(403)
     CHARACTER (LEN=1) :: diag_addr(32)
     CHARACTER (LEN=1) :: read_value
     CHARACTER (LEN=1) :: status_value
     CHARACTER (LEN=2) :: status_version
     CHARACTER (LEN=1) :: read_riu_id
     CHARACTER (LEN=1) :: status_riu_id
     CHARACTER (LEN=1) :: spare2
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci1_T2_T

  TYPE Sci2_T2_T     ! science packet 2 Types II and III
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: TYPE
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=1) :: spare1
!
!! FB_dat contains the filter banks, in this order:
!!  MB01, FB06, status, MB03, FB10, status, MB02, FB11, status, FB15, FB14,
!!  MB05, MB04, status, FB17, FB16, FB19, FB18
!
     CHARACTER (LEN=2) :: FB_dat(11+25+1+11+25+1+11+25+1+2*25+2*11+1+4*25)
                          ! filter, mid-band plus status data
     CHARACTER (LEN=1) :: DACS(261)
     CHARACTER (LEN=24) :: thz_scan_switch
     CHARACTER (LEN=80) :: laser_lo
     CHARACTER (LEN=1) :: TM03_stat(2)
     CHARACTER (LEN=24) :: ghz_ant_scan
     CHARACTER (LEN=24) :: ghz_switch
     CHARACTER (LEN=1) :: spare2(10)
     CHARACTER (LEN=1) :: attenuation(5)
     CHARACTER (LEN=1) :: spare3(2)
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci2_T2_T

  TYPE (Sci1_T1_T), TARGET :: sci1_T1
  TYPE (Sci2_T1_T), TARGET :: sci2_T1

  TYPE (Sci1_T2_T), TARGET :: sci1_T2
  TYPE (Sci2_T2_T), TARGET :: sci2_T2

!! Formats to convert Type I raw Science packets internally:

  CHARACTER (LEN=*), PARAMETER :: Sci1_T1_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,2(25A2),A2,2(25A2),A2,4(25A2),11A2,25A2,A2,&
       & 2(25A2),2A1,3(2A1),2A13,A1,A13,A1,5A1,320A1,32A1,2A1,A2,2A1,10A1,A2)"

  CHARACTER (LEN=*), PARAMETER :: Sci2_T1_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,2(11A2,25A2,A2),2(25A2),2(11A2),A2,4(25A2),&
       & 320A1,A24,A80,2A1,2A24,25A1,5A1,2A1,A2)"

!! Formats to convert Type II/III raw Science packets internally:

  CHARACTER (LEN=*), PARAMETER :: Sci1_T2_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,2(25A2),A2,2(25A2),A2,6(25A2),2A1,3(2A1),&
       & 2A13,A1,A13,A1,5A1,403A1,32A1,2A1,A2,3A1,A2)"

  CHARACTER (LEN=*), PARAMETER :: Sci2_T2_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,3(11A2,25A2,A2),2(25A2),2(11A2),A2,4(25A2),&
       & 261A1,A24,A80,2A1,2A24,10A1,5A1,2A1,A2)"

!! Character (LEN=1) pointer types

  TYPE CS1PTR
     CHARACTER (LEN=1), POINTER :: ptr
  END TYPE CS1PTR

  TYPE C1PTR
     CHARACTER (LEN=1), DIMENSION (:), POINTER :: ptr
  END TYPE C1PTR

!! Character (LEN=2) pointer types

  TYPE CS2PTR
     CHARACTER (LEN=2), POINTER :: ptr
  END TYPE CS2PTR

  TYPE C2PTR
     CHARACTER (LEN=2), DIMENSION (:), POINTER :: ptr
  END TYPE C2PTR

  TYPE CS13PTR
     CHARACTER (LEN=13), POINTER :: ptr
  END TYPE CS13PTR

!! Generic pointers (used after reading data and determining data format)

  TYPE Sci_cptr_T

     TYPE (CS2PTR) :: MAF(2)         ! MAF pointers
     TYPE (CS1PTR) :: MIF(2)         ! MIF pointers
     TYPE (CS2PTR) :: orbit(2)       ! orbit pointers
     TYPE (C2PTR)  :: FB(FBnum)      ! Filter bank pointers
     TYPE (C2PTR)  :: MB(MBnum)      ! Mid-Band pointers
     TYPE (C1PTR)  :: DACS(DACSnum)  ! DACS pointers
     TYPE (CS13PTR) :: WF(WFnum)     ! Wide filter pointers
     TYPE (C1PTR) :: GSN             ! GHz Switch Network
     TYPE (CS1PTR) :: THzSw          ! THz Switch
     CHARACTER (LEN=80), POINTER :: LLO_DN         ! Laser LO data pointers
     CHARACTER (LEN=24), POINTER :: GHz_ant_scan, GHz_sw, THz_sw
     TYPE (C1PTR) :: Attenuation     ! Attenuation read backs

  END TYPE Sci_cptr_T

  TYPE (Sci_cptr_T) :: Sci_cptr(2)

!! Filter Bank offsets within data field in Type I science packets:

  INTEGER, PARAMETER :: FB_T1_offset(FBnum) = (/ 1, 26, 52, 77, 103, 128, &
       & 153, 178, 214, 240, 265, 12, 49, 75, 100, 148, 173, 198, 223 /)

!! Filter Bank offsets within data field in Type II/III science packets:

  INTEGER, PARAMETER :: FB_T2_offset(FBnum) = (/ 1, 26, 52, 77, 103, 128, &
       & 153, 178, 203, 228, 12, 49, 86, 112, 137, 185, 210, 235, 260 /)

!! Mid-band offsets within data field in Type I science packets:

  INTEGER, PARAMETER :: MB_T1_offset(MBnum) = (/ 203, 38, 1, 136, 125 /)

!! Mid-band offsets within data field in Type II/III science packets:

  INTEGER, PARAMETER :: MB_T2_offset(MBnum) = (/ 1, 75, 38, 173, 162 /)

!! Filter bank numbers within the Type I science packets

  INTEGER, PARAMETER :: Pkt1_T1_FB(11) = (/ 12, 1, 4, 2, 7, 3, 5, 13, 6, 9, 8 /)
  INTEGER, PARAMETER :: Pkt2_T1_FB(8) = (/ 10, 11, 15, 14, 17, 16, 19, 18 /)

!! Filter bank numbers within the Type II/III science packets

  INTEGER, PARAMETER :: Pkt1_T2_FB(10) = (/ 12, 1, 4, 2, 7, 3, 5, 13, 9, 8 /)
  INTEGER, PARAMETER :: Pkt2_T2_FB(9) = (/ 6, 10, 11, 15, 14, 17, 16, 19, 18 /)

!! Attenuation setting read back structure:

  TYPE Atten_T
     INTEGER :: RIU
     INTEGER :: Addr
     INTEGER :: Mask
     INTEGER :: Value
  END TYPE Atten_T

!! L0 science data packet for 1 MIF:

  TYPE Sci_pkt_T
     REAL(r8) :: secTAI
     INTEGER :: MAFno
     INTEGER :: MIFno
     INTEGER :: Orbit
     INTEGER :: FB(FBchans,FBnum)
     INTEGER :: MB(MBchans,MBnum)
     INTEGER :: WF(WFchans,WFnum)
     REAL :: DACS(DACSchans,DACSnum)
     REAL :: APE_pos(2), ASA_pos(2), GSA_pos(2), TSSA_pos(2)
     REAL :: scAngleG                  ! Boresight wrt. spc +x (GHz)
     REAL :: scAngleT                  ! Boresight wrt. spc +x (THz)
     CHARACTER(len=1) :: GHz_sw_pos
     CHARACTER(len=1) :: THz_sw_pos
     CHARACTER(len=80) :: LLO_DN       ! Laser LO DN data
     REAL :: LLO_EU(16)                ! Laser LO EU data
     CHARACTER(len=70) :: PLL_DN       ! Phase Lock Loop data
     INTEGER :: GSN(4)                 ! GHz Switch Network readings
     INTEGER :: THzSw                  ! THz Switch reading
     INTEGER :: BandSwitch(5)          ! band associated with each switch
     TYPE (BankLogical_T) :: MaxAtten  ! Whether in maximum attenuation
     LOGICAL :: AttenMaxed             ! Some attenuation is maximum
     LOGICAL :: CRC_good
  END TYPE Sci_pkt_T

!! L0 science packet for 1 MIF:

  TYPE (Sci_pkt_T) :: Sci_pkt

!! L0 science for 1 MAF:

  TYPE (Sci_pkt_T) :: SciMAF(0:(MaxMIFs-1))

!! L0 THz science data packet for 1 MIF:

  TYPE THz_Sci_pkt_T
     REAL(r8) :: secTAI
     INTEGER :: MAFno
     INTEGER :: MIFno
     INTEGER :: Orbit
     INTEGER :: FB(THzChans,THzNum)
     LOGICAL :: MaxAtten(THzNum)       ! Whether in maximum attenuation
     LOGICAL :: AttenMaxed             ! Some attenuation is maximum
     REAL :: TSSA_pos(2)
     REAL :: scAngle                   ! Boresight wrt. spc +x
     CHARACTER(len=1) :: SwMirPos
     REAL :: LLO_Bias                  ! LLO Bias V
     INTEGER :: BandSwitch(4:5)        ! band associated with switches #4 and #5
     LOGICAL :: CRC_good
  END TYPE THz_Sci_pkt_T

!! L0 THz science for 1 MAF:

  TYPE (THz_Sci_pkt_T) :: THzSciMAF(0:(MaxMIFs-1))

!! L0 DACS data for 1 MIF:

  TYPE DACS_pkt_T
     INTEGER :: C_K(DACSchans,DACSnum)  ! C for compressed; K for uncompressed
     INTEGER :: D(4,DACSnum)            ! State counters
     INTEGER :: TP(DACSnum)
     INTEGER :: DIO(DACSnum)
     INTEGER :: LO(DACSnum)
     INTEGER :: Zlag(DACSnum)
     LOGICAL :: Compressed
  END TYPE DACS_pkt_T

!! L0 DACS data for 1 MIF:

  TYPE (DACS_pkt_T) :: DACS_pkt

!! L0 DACS data for 1 MAF:

  TYPE (DACS_pkt_T) :: DACS_MAF(0:(MaxMIFs-1))

!! Band attenuation table:

  TYPE (Atten_T), PARAMETER :: BandAtten(NumBands) = (/ &
       Atten_T ( 91, 49168, 255, 0), & ! 5B, C010, FF, 0   Band  1
       Atten_T ( 95, 49170, 255, 0), & ! 5F, C012, FF, 0   Band  2
       Atten_T ( 95, 49171, 255, 0), & ! 5F, C013, FF, 0   Band  3
       Atten_T ( 95, 49172, 255, 0), & ! 5F, C014, FF, 0   Band  4
       Atten_T ( 95, 49173, 255, 0), & ! 5F, C015, FF, 0   Band  5
       Atten_T ( 95, 49174, 255, 0), & ! 5F, C016, FF, 0   Band  6
       Atten_T ( 96, 49169, 255, 0), & ! 60, C011, FF, 0   Band  7
       Atten_T ( 96, 49170, 255, 0), & ! 60, C012, FF, 0   Band  8
       Atten_T ( 96, 49171, 255, 0), & ! 60, C013, FF, 0   Band  9
       Atten_T (110, 49170, 255, 0), & ! 6E, C012, FF, 0   Band 10
       Atten_T (110, 49171, 255, 0), & ! 6E, C013, FF, 0   Band 11
       Atten_T (110, 49172, 255, 0), & ! 6E, C014, FF, 0   Band 12
       Atten_T (110, 49173, 255, 0), & ! 6E, C015, FF, 0   Band 13
       Atten_T (110, 49174, 255, 0), & ! 6E, C016, FF, 0   Band 14
       Atten_T ( 84, 49170, 255, 0), & ! 54, C012, FF, 0   Band 15
       Atten_T ( 84, 49169, 255, 0), & ! 54, C011, FF, 0   Band 16
       Atten_T ( 84, 49168, 255, 0), & ! 54, C010, FF, 0   Band 17
       Atten_T ( 84, 49173, 255, 0), & ! 54, C015, FF, 0   Band 18
       Atten_T ( 84, 49172, 255, 0), & ! 54, C014, FF, 0   Band 19
       Atten_T ( 84, 49171, 255, 0), & ! 54, C013, FF, 0   Band 20
       Atten_T ( 92, 49168, 255, 0), & ! 5C, C010, FF, 0   Band 21
       Atten_T (111, 49168, 255, 0), & ! 6F, C010, FF, 0   Band 22 ------|
       Atten_T ( 97, 49168, 255, 0), & ! 61, C010, FF, 0   Band 23       |
       Atten_T (111, 49169, 255, 0), & ! 6F, C011, FF, 0   Band 24  DACS |
       Atten_T ( 97, 49169, 255, 0), & ! 61, C011, FF, 0   Band 25       |
       Atten_T ( 97, 49169, 255, 0), & ! 61, C011, FF, 0   Band 26 ------|
       Atten_T ( 95, 49174, 255, 0), & ! 5F, C016, FF, 0   Band 27 - B 6
       Atten_T (110, 49171, 255, 0), & ! 6E, C013, FF, 0   Band 28 - B 11
       Atten_T (110, 49170, 255, 0), & ! 6E, C012, FF, 0   Band 29 - B 10
       Atten_T (110, 49175, 255, 0), & ! 6E, C017, FF, 0   Band 30 - MB 4
       Atten_T (110, 49175, 255, 0), & ! 6E, C017, FF, 0   Band 31 - MB 5
       Atten_T ( 91, 49169, 255, 0), & ! 5B, C011, FF, 0   Band 32 - R1A - WF 1
       Atten_T ( 96, 49172, 255, 0), & ! 60, C014, FF, 0   Band 33 - R3  - WF 2
       Atten_T ( 92, 49169, 255, 0)  & ! 5C, C011, FF, 0   Band 34 - R1B - WF 3
       /)

!! Default APE pos2 values for use in simulations:

  REAL, PARAMETER :: APE2_dflt(0:(MaxMIFs-1)) = (/ &
       357.049, 357.044, 357.038, 357.030, 357.023, 357.017, 357.010, 357.003, &
       356.995, 356.989, 356.983, 356.975, 356.968, 356.962, 356.955, 356.947, &
       356.941, 356.933, 356.927, 356.920, 356.913, 356.907, 356.899, 356.892, &
       356.885, 356.879, 356.871, 356.864, 356.858, 356.851, 356.844, 356.837, &
       356.831, 356.824, 356.816, 356.809, 356.802, 356.796, 356.789, 356.782, &
       356.776, 356.768, 356.761, 356.755, 356.748, 356.740, 356.733, 356.727, &
       356.721, 356.713, 356.706, 356.700, 356.693, 356.685, 356.679, 356.672, &
       356.665, 356.658, 356.651, 356.645, 356.637, 356.630, 356.623, 356.617, &
       356.609, 356.602, 356.596, 356.584, 356.565, 356.547, 356.529, 356.510, &
       356.491, 356.473, 356.455, 356.436, 356.418, 356.400, 356.382, 356.363, &
       356.345, 356.326, 356.308, 356.290, 356.271, 356.253, 356.235, 356.216, &
       356.198, 356.179, 356.161, 356.143, 356.125, 356.106, 356.088, 356.069, &
       356.051, 356.032, 356.014, 355.996, 355.978, 355.960, 355.941, 355.923, &
       355.904, 355.886, 355.868, 355.849, 355.818, 355.771, 355.724, 355.677, &
       355.631, 355.584, 355.537, 355.490, 355.444, 355.397, 355.350, 355.303, &
       355.277, 355.288, 355.318, 355.361, 355.416, 355.483, 355.559, 355.643, &
       355.734, 355.831, 355.931, 356.035, 356.140, 356.245, 356.349, 356.452, &
       356.549, 356.643, 356.729, 356.808, 356.879, 356.938, 356.986, 357.022, &
       357.042, 357.049,  -999.9, 357.049,  -999.9,  -999.9 /)

!! Default TSE pos2 values for use in simulations:

  REAL, PARAMETER :: TSE2_dflt(0:(MaxMIFs-1)) = (/ &
       359.500,  -999.9, 359.500, 359.458, 359.372, 359.287, 359.202, 359.157, &
       359.151, 359.144, 359.138, 359.132, 359.126, 359.119, 359.113, 359.107, &
       359.101, 359.094, 359.088, 359.082, 359.076, 359.070, 359.063, 359.057, &
       359.051, 359.045, 359.038, 359.032, 359.026, 359.020, 359.013, 359.007, &
       359.001, 358.995, 358.988, 358.982, 358.976, 358.970, 358.963, 358.957, &
       358.951, 358.945, 358.939, 358.932, 358.926, 358.920, 358.914, 358.907, &
       358.901, 358.895, 358.889, 358.882, 358.876, 358.870, 358.864, 358.857, &
       358.851, 358.845, 358.836, 358.824, 358.812, 358.800, 358.788, 358.776, &
       358.764, 358.752, 358.740, 358.728, 358.716, 358.704, 358.692, 358.680, &
       358.668, 358.656, 358.644, 358.632, 358.620, 358.608, 358.596, 358.584, &
       358.572, 358.560, 358.548, 358.536, 358.524, 358.512, 358.500, 358.488, &
       358.476, 358.464, 358.452, 358.440, 358.428, 358.416, 358.404, 358.392, &
       358.380, 358.368, 358.356, 358.344, 358.332, 358.320, 358.308, 358.296, &
       358.284, 358.255, 358.210, 358.165, 358.120, 358.075, 358.030, 357.985, &
       357.941, 357.896, 357.851, 357.806, 357.760, 357.729, 357.332, 356.845, &
       356.798, 356.797, 356.797, 356.798, 356.797, 356.797, 356.798, 356.797, &
       356.797, 356.797, 356.797, 330.202, 210.671, 178.750, 178.748, 178.748, &
       178.748, 178.748, 178.748, 178.747, 178.748, 151.986, 32.2582, &
       0.113705, 359.962, 359.790, 359.620, 359.500, -999.9, -999.9 /)

CONTAINS

  SUBROUTINE InitSciPointers

    INTEGER :: i, pkt_offset

    !! Initialize Science packet type pointers

    type_sci1 => l0_sci1(16:16)
    type_sci2 => l0_sci2(16:16)

    !! Initialize pointers for Type I science packets:

    !! MAF, MIF, orbit:

    sci_cptr(1)%MAF(1)%ptr => sci1_T1%MAF
    sci_cptr(1)%MAF(2)%ptr => sci2_T1%MAF
    sci_cptr(1)%MIF(1)%ptr => sci1_T1%MIF
    sci_cptr(1)%MIF(2)%ptr => sci2_T1%MIF
    sci_cptr(1)%orbit(1)%ptr => sci1_T1%orbit
    sci_cptr(1)%orbit(2)%ptr => sci2_T1%orbit

    !! Filter Bank pointers

    !!  Science packet #1

    DO i = 1, SIZE (Pkt1_T1_FB)
       sci_cptr(1)%FB(Pkt1_T1_FB(i))%ptr => &
        sci1_T1%FB_dat(FB_T1_offset(i):FB_T1_offset(i)+24)
    ENDDO

    !!  Science packet #2

    pkt_offset = SIZE (Pkt1_T1_FB)
    DO i = 1, SIZE (Pkt2_T1_FB)
       sci_cptr(1)%FB(Pkt2_T1_FB(i))%ptr => &
        sci2_T1%FB_dat(FB_T1_offset(i+pkt_offset):FB_T1_offset(i+pkt_offset)+24)
    ENDDO

    !! Mid-Band pointers

    sci_cptr(1)%MB(1)%ptr => sci1_T1%FB_dat(MB_T1_offset(1):MB_T1_offset(1)+10)
    DO i = 2, MBNum
       sci_cptr(1)%MB(i)%ptr => &
            sci2_T1%FB_dat(MB_T1_offset(i):MB_T1_offset(i)+10)
    ENDDO

    !! Wide filter pointers

    sci_cptr(1)%WF(1)%ptr => sci1_T1%GM10_stat
    sci_cptr(1)%WF(2)%ptr => sci1_T1%GM13_stat
    sci_cptr(1)%WF(3)%ptr => sci1_T1%GM11_stat

    !! DACS pointers

    sci_cptr(1)%DACS(1)%ptr => sci1_T1%DACS(1:158)
    sci_cptr(1)%DACS(2)%ptr => sci1_T1%DACS(159:)
    sci_cptr(1)%DACS(3)%ptr => sci2_T1%DACS(1:158)
    sci_cptr(1)%DACS(4)%ptr => sci2_T1%DACS(159:)

    sci_cptr(1)%GHz_ant_scan => sci2_T1%ghz_ant_scan   ! GHz antenna/scan data
    sci_cptr(1)%GHz_sw => sci2_T1%ghz_switch           ! GHz switching data
    sci_cptr(1)%THz_sw => sci2_T1%thz_scan_switch      ! THz switching data
    sci_cptr(1)%GSN%ptr => sci1_T1%GM15_stat(2:5)      ! GSN data
    sci_cptr(1)%THzSw%ptr => sci2_T1%TM03_stat(2)      ! THz switch data
    sci_cptr(1)%LLO_DN => sci2_T1%laser_lo
    sci_cptr(1)%Attenuation%ptr => sci2_T1%Attenuation

    !! Initialize pointers for Type II/III science packets:

    !! MAF, MIF, orbit:

    sci_cptr(2)%MAF(1)%ptr => sci1_T2%MAF
    sci_cptr(2)%MAF(2)%ptr => sci2_T2%MAF
    sci_cptr(2)%MIF(1)%ptr => sci1_T2%MIF
    sci_cptr(2)%MIF(2)%ptr => sci2_T2%MIF
    sci_cptr(2)%orbit(1)%ptr => sci1_T2%orbit
    sci_cptr(2)%orbit(2)%ptr => sci2_T2%orbit
    sci_cptr(2)%Attenuation%ptr => sci2_T2%Attenuation

    !! Filter Bank pointers

    !!  Science packet #1

    DO i = 1, SIZE (Pkt1_T2_FB)
       sci_cptr(2)%FB(Pkt1_T2_FB(i))%ptr => &
        sci1_T2%FB_dat(FB_T2_offset(i):FB_T2_offset(i)+24)
    ENDDO

    !!  Science packet #2

    pkt_offset = SIZE (Pkt1_T2_FB)
    DO i = 1, SIZE (Pkt2_T2_FB)
       sci_cptr(2)%FB(Pkt2_T2_FB(i))%ptr => &
        sci2_T2%FB_dat(FB_T2_offset(i+pkt_offset):FB_T2_offset(i+pkt_offset)+24)
    ENDDO

    !! Mid-Band pointers

    sci_cptr(2)%MB(1)%ptr => sci1_T2%FB_dat(MB_T2_offset(1):MB_T2_offset(1)+10)
    DO i = 2, MBNum
       sci_cptr(2)%MB(i)%ptr => &
            sci2_T2%FB_dat(MB_T2_offset(i):MB_T2_offset(i)+10)
    ENDDO

    !! Wide filter pointers

    sci_cptr(2)%WF(1)%ptr => sci1_T2%GM10_stat
    sci_cptr(2)%WF(2)%ptr => sci1_T2%GM13_stat
    sci_cptr(2)%WF(3)%ptr => sci1_T2%GM11_stat

    !! DACS pointers

    sci_cptr(2)%DACS(1)%ptr => sci1_T2%DACS
    sci_cptr(2)%DACS(2)%ptr => sci2_T2%DACS
    sci_cptr(2)%DACS(3)%ptr => sci1_T2%DACS
    sci_cptr(2)%DACS(4)%ptr => sci2_T2%DACS

    sci_cptr(2)%GHz_ant_scan => sci2_T2%ghz_ant_scan   ! GHz antenna/scan data
    sci_cptr(2)%GHz_sw => sci2_T2%ghz_switch           ! GHz switching data
    sci_cptr(2)%THz_sw => sci2_T2%thz_scan_switch      ! THz switching data
    sci_cptr(2)%GSN%ptr => sci1_T2%GM15_stat(2:5)      ! GSN data
    sci_cptr(2)%THzSw%ptr => sci2_T2%TM03_stat(2)      ! THz switch data
    sci_cptr(2)%LLO_DN => sci2_T2%laser_lo

  END SUBROUTINE InitSciPointers

END MODULE L0_sci_tbls

! $Log$
! Revision 2.6  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.5  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.4  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.3  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:47:49  perun
! Version 0.5 commit
!
