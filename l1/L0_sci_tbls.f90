! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L0_sci_tbls  ! Define L0 science tables
!=============================================================================

  USE MLSL1Common, ONLY: R8, FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
       DACSnum, DACSchans, MaxMIFs

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

!! Using format of 03/10/01

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
     CHARACTER (LEN=1) :: laser_lo(80)
     CHARACTER (LEN=1) :: TM03_stat(2)
     CHARACTER (LEN=24) :: ghz_ant_scan
     CHARACTER (LEN=24) :: ghz_switch
     CHARACTER (LEN=9) :: eng_data
     CHARACTER (LEN=1) :: spare2(23)
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
     CHARACTER (LEN=1) :: laser_lo(80)
     CHARACTER (LEN=1) :: TM03_stat(2)
     CHARACTER (LEN=24) :: ghz_ant_scan
     CHARACTER (LEN=24) :: ghz_switch
     CHARACTER (LEN=1) :: spare2
     CHARACTER (LEN=9) :: eng_data
     CHARACTER (LEN=1) :: spare3(7)
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
       & 320A1,A24,80A1,2A1,2A24,A1,A9,A1,A2)"

!! Formats to convert Type II/III raw Science packets internally:

  CHARACTER (LEN=*), PARAMETER :: Sci1_T2_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,2(25A2),A2,2(25A2),A2,6(25A2),2A1,3(2A1),&
       & 2A13,A1,A13,A1,5A1,403A1,32A1,2A1,A2,3A1,A2)"

  CHARACTER (LEN=*), PARAMETER :: Sci2_T2_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,3(11A2,25A2,A2),2(25A2),2(11A2),A2,4(25A2),&
       & 261A1,A24,80A1,2A1,2A24,A1,A9,A1,A2)"

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
     TYPE (C1PTR)  :: LLO_data       ! Laser LO data pointers
     CHARACTER (LEN=24), POINTER :: GHz_ant_scan, GHz_sw, THz_sw

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
     REAL :: GHz_sw_angle(2)
     REAL :: THz_sw_angle(2)
     REAL :: scan_angle(2)
     CHARACTER(len=2) :: GHz_sw_pos
     CHARACTER(len=2) :: THz_sw_pos
     CHARACTER(len=1) :: LLO_data(80)  ! Laser LO data
     INTEGER :: GSN(4)                 ! GHz Switch Network readings
     INTEGER :: THzSw                  ! THz Switch reading
     INTEGER :: BandSwitch(5)          ! band associated with each switch
     LOGICAL :: CRC_good
  END TYPE Sci_pkt_T

!! L0 science packet for 1 MIF:

  TYPE (Sci_pkt_T) :: Sci_pkt

!! L0 science for 1 MAF:

  TYPE (Sci_pkt_T) :: SciMAF(0:(MaxMIFs-1))

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
    sci_cptr(1)%LLO_data%ptr => sci2_T1%laser_lo

    !! Initialize pointers for Type II/III science packets:

    !! MAF, MIF, orbit:

    sci_cptr(2)%MAF(1)%ptr => sci1_T2%MAF
    sci_cptr(2)%MAF(2)%ptr => sci2_T2%MAF
    sci_cptr(2)%MIF(1)%ptr => sci1_T2%MIF
    sci_cptr(2)%MIF(2)%ptr => sci2_T2%MIF
    sci_cptr(2)%orbit(1)%ptr => sci1_T2%orbit
    sci_cptr(2)%orbit(2)%ptr => sci2_T2%orbit

    !! Filter Bank pointers

    !!  Science packet #2

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
    sci_cptr(2)%LLO_data%ptr => sci2_T2%laser_lo

  END SUBROUTINE InitSciPointers

END MODULE L0_sci_tbls

! $Log$
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:47:49  perun
! Version 0.5 commit
!
