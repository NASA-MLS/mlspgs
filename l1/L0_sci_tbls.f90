! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
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
     REAL :: APE_pos(2), ASA_pos(2), GSM_pos(2), TSSM_pos(2)
     REAL :: APE_theta, ASA_theta, GSM_theta, TSSM_theta
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
     TYPE (BankLogical_T) :: DeltaAtten  ! Whether attenuation changed
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
     REAL :: TSSM_pos(2), TSSM_theta
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

!! Attenuation setting read-back structure:

  TYPE Atten_T
     INTEGER :: RIU
     INTEGER :: Addr
     INTEGER :: Mask
     INTEGER :: Value
     INTEGER :: Band, Link
  END TYPE Atten_T

!! Band attenuation table:

  TYPE (Atten_T), PARAMETER :: BandAtten(NumBands-2) = (/ &
       Atten_T ( 91, 49168, 255, 0,  1, 22), & ! 5B, C010, FF, 0   Band  1
       Atten_T ( 95, 49170, 255, 0,  2, 23), & ! 5F, C012, FF, 0   Band  2
       Atten_T ( 95, 49171, 255, 0,  3,  0), & ! 5F, C013, FF, 0   Band  3
       Atten_T ( 95, 49172, 255, 0,  4,  0), & ! 5F, C014, FF, 0   Band  4
       Atten_T ( 95, 49173, 255, 0,  5,  0), & ! 5F, C015, FF, 0   Band  5
       Atten_T ( 95, 49174, 255, 0,  6, 27), & ! 5F, C016, FF, 0   Band  6
       Atten_T ( 96, 49169, 255, 0,  7, 24), & ! 60, C011, FF, 0   Band  7
       Atten_T ( 96, 49170, 255, 0,  8,  0), & ! 60, C012, FF, 0   Band  8
       Atten_T ( 96, 49171, 255, 0,  9, 25), & ! 60, C013, FF, 0   Band  9
       Atten_T (110, 49170, 255, 0, 10, 29), & ! 6E, C012, FF, 0   Band 10
       Atten_T (110, 49171, 255, 0, 11, 28), & ! 6E, C013, FF, 0   Band 11
       Atten_T (110, 49172, 255, 0, 12,  0), &  ! 6E, C014, FF, 0  Band 12
       Atten_T (110, 49173, 255, 0, 13,  0), &  ! 6E, C015, FF, 0  Band 13
       Atten_T (110, 49174, 255, 0, 14,  0), &  ! 6E, C016, FF, 0  Band 14
       Atten_T ( 84, 49170, 255, 0, 15,  0), &  ! 54, C012, FF, 0  Band 15
       Atten_T ( 84, 49169, 255, 0, 16,  0), &  ! 54, C011, FF, 0  Band 16
       Atten_T ( 84, 49168, 255, 0, 17,  0), &  ! 54, C010, FF, 0  Band 17
       Atten_T ( 84, 49173, 255, 0, 18,  0), &  ! 54, C015, FF, 0  Band 18
       Atten_T ( 84, 49172, 255, 0, 19,  0), &  ! 54, C014, FF, 0  Band 19
       Atten_T ( 84, 49171, 255, 0, 20,  0), &  ! 54, C013, FF, 0  Band 20
       Atten_T ( 92, 49168, 255, 0, 21, 26), & ! 5C, C010, FF, 0   Band 21
       Atten_T (111, 49168, 255, 0, 22,  1), & ! 6F, C010, FF, 0   Band 22
       Atten_T ( 97, 49168, 255, 0, 23,  2), & ! 61, C010, FF, 0   Band 23
       Atten_T (111, 49169, 255, 0, 24,  7), & ! 6F, C011, FF, 0   Band 24
       Atten_T ( 97, 49169, 255, 0, 25,  9), & ! 61, C011, FF, 0   Band 25/26
       Atten_T (110, 49175, 255, 0, 30, 31), & ! 6E, C017, FF, 0   Band 30/31
       Atten_T ( 91, 49169, 255, 0, 32,  0), & ! 5B, C011, FF, 0   Band 32
       Atten_T ( 91, 49170, 255, 0, 32,  0), & ! 5B, C012, FF, 0   Band 32
       Atten_T ( 96, 49172, 255, 0, 33,  0), & ! 60, C014, FF, 0   Band 33
       Atten_T ( 96, 49173, 255, 0, 33,  0), & ! 60, C0145, FF, 0  Band 33
       Atten_T ( 92, 49169, 255, 0, 34,  0), & ! 5C, C011, FF, 0   Band 34
       Atten_T ( 92, 49170, 255, 0, 34,  0) &  ! 5C, C012, FF, 0   Band 34
       /)

!! Default APE theta values for use in simulations:

  REAL, PARAMETER :: APE_theta_dflt(0:(MaxMIFs-1)) = (/ &
       357.054077, 357.054626, 357.053955, 357.053131, 357.047852, 357.041138, &
       357.034424, 357.027649, 357.020782, 357.013214, 357.006287, 356.999695, &
       356.992676, 356.986328, 356.978577, 356.971558, 356.965027, 356.958191, &
       356.950928, 356.944214, 356.937225, 356.930786, 356.923676, 356.916199, &
       356.909760, 356.902863, 356.896057, 356.889343, 356.881958, 356.875549, &
       356.868744, 356.861572, 356.854614, 356.847412, 356.840759, 356.834351, &
       356.826874, 356.819763, 356.813293, 356.806274, 356.799866, 356.792267, &
       356.785126, 356.778625, 356.771820, 356.764740, 356.757874, 356.750336, &
       356.743896, 356.737152, 356.729980, 356.723419, 356.715881, 356.709198, &
       356.702393, 356.695282, 356.688202, 356.681671, 356.674561, 356.667908, &
       356.660309, 356.653564, 356.647156, 356.640198, 356.632874, 356.625854, &
       356.618988, 356.612671, 356.605469, 356.598145, 356.586060, 356.567963, &
       356.549347, 356.531250, 356.513000, 356.494934, 356.476685, 356.458160, &
       356.440033, 356.421631, 356.403198, 356.384857, 356.366272, 356.348206, &
       356.329773, 356.311127, 356.292816, 356.274719, 356.256470, 356.237823, &
       356.219788, 356.201416, 356.183014, 356.164368, 356.146118, 356.127747, &
       356.109314, 356.090790, 356.072693, 356.054321, 356.036072, 356.017456, &
       355.999084, 355.981110, 355.962738, 355.944458, 355.926056, 355.907867, &
       355.889496, 355.871094, 355.852722, 355.821472, 355.774963, 355.727386, &
       355.681091, 355.633942, 355.587402, 355.540466, 355.493958, 355.446655, &
       355.399780, 355.353210, 355.306183, 355.282623, 355.293762, 355.323059, &
       355.365845, 355.421570, 355.488007, 355.564270, 355.648132, 355.739014, &
       355.835449, 355.936096, 356.039490, 356.144409, 356.249878, 356.353607, &
       356.456207, 356.553894, 356.647034, 356.733826, 356.813080, 356.882874, &
       356.942322, 356.990479, -999.9,     -999.9,     -999.9,     -999.9 /)

!! Default TSSM theta values for use in simulations:

  REAL, PARAMETER :: TSSM_theta_dflt(0:(MaxMIFs-1)) = (/ &
       359.499664, 359.499969, 359.500000, 359.457245, 359.371796, 359.287109, &
       359.202026, 359.156891, 359.150726, 359.144409, 359.138245, 359.131958, &
       359.125702, 359.119476, 359.113312, 359.106995, 359.100769, 359.094543, &
       359.088287, 359.082001, 359.075836, 359.069550, 359.063293, 359.057098, &
       359.050842, 359.044586, 359.038391, 359.032074, 359.025909, 359.019653, &
       359.013428, 359.007202, 359.000885, 358.994751, 358.988464, 358.982269, &
       358.975952, 358.969727, 358.963501, 358.957275, 358.951080, 358.944824, &
       358.938538, 358.932373, 358.926025, 358.919830, 358.913574, 358.907288, &
       358.901123, 358.894897, 358.888672, 358.882385, 358.876160, 358.869934, &
       358.863708, 358.857483, 358.851288, 358.844940, 358.835876, 358.823883, &
       358.811829, 358.799835, 358.787842, 358.775940, 358.763916, 358.751923, &
       358.739868, 358.727875, 358.715912, 358.703857, 358.691925, 358.679932, &
       358.667908, 358.655945, 358.643982, 358.631927, 358.619995, 358.607971, &
       358.596100, 358.584015, 358.572052, 358.560120, 358.548096, 358.536072, &
       358.524109, 358.512054, 358.500153, 358.488068, 358.476013, 358.464111, &
       358.452118, 358.440094, 358.428131, 358.416107, 358.404083, 358.392090, &
       358.380127, 358.368103, 358.356171, 358.344147, 358.332092, 358.320099, &
       358.308167, 358.296173, 358.284119, 358.255463, 358.210388, 358.165283, &
       358.120331, 358.075256, 358.030426, 357.985535, 357.940704, 357.895630, &
       357.850800, 357.805695, 357.760498, 357.656311, 357.147430, 356.815796, &
       356.798615, 356.798615, 356.797241, 356.797241, 356.797913, 356.797913, &
       356.797241, 356.797241, 356.797913, 356.797241, 356.797913, 322.130554, &
       218.175171, 178.761719, 178.744553, 178.747986, 178.747299, 178.747986, &
       178.745926, 178.747986, 178.748672, 143.791550, 39.4406586, 0.23229003, &
       90.0744781, 359.884460, 359.715027, -999.9,     -999.9,     -999.9 /)

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
! Revision 2.8  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.7  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
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
