! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L0_sci_tbls  ! Define L0 science tables
!=============================================================================

  USE MLSCommon
  USE MLSL1Common

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

!! Using format of 01/07/00

  TYPE Sci1_uc_T     ! science packet 1 with uncompressed filter bank data
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: type
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=1) :: comp
     CHARACTER (LEN=2) :: mifs_per_maf
!
!! FB_dat contains the filter banks, in this order:
!!  FB12, FB01, FB04, FB02, FB03, FB07, FB05, FB13,  FB06, MB01, FB08, FB09
!
     CHARACTER (LEN=2) :: FB_dat(9*25+11+2*25)  ! filter and mid-band data
     CHARACTER (LEN=7) :: eng_dat
     CHARACTER (LEN=1) :: TM03_stat(2)
     CHARACTER (LEN=1) :: GM15_stat(2)
     CHARACTER (LEN=1) :: GM05_stat
     CHARACTER (LEN=1) :: GM06_stat
     CHARACTER (LEN=1) :: GM07_stat
     CHARACTER (LEN=1) :: GM08_stat
     CHARACTER (LEN=1) :: GM10_stat(2)
     CHARACTER (LEN=1) :: GM11_stat(2)
     CHARACTER (LEN=1) :: GM14_stat(6)
     CHARACTER (LEN=1) :: GM12_stat(5)
     CHARACTER (LEN=1) :: GM13_stat(3)
     CHARACTER (LEN=3) :: DACS(DACchans)
     CHARACTER (LEN=1) :: read_value
     CHARACTER (LEN=1) :: status_value
     CHARACTER (LEN=2) :: status_version
     CHARACTER (LEN=1) :: read_riu_id
     CHARACTER (LEN=1) :: status_riu_id
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci1_uc_T

  TYPE Sci2_uc_T     ! science packet 2 with uncompressed filter bank data
     SEQUENCE
     TYPE (CCSDS_T)    :: CCSDS
     CHARACTER (LEN=1) :: type
     CHARACTER (LEN=1) :: MIF
     CHARACTER (LEN=2) :: MAF
     CHARACTER (LEN=2) :: orbit
     CHARACTER (LEN=3) :: tot_mifs
!
!! FB_dat contains the filter banks, in this order:
!!  FB11, MB03, FB14, FB15, MB04, MB05, FB16, FB17, FB18, FB19, FB10
!
     CHARACTER (LEN=2) :: FB_dat(25+11+2*25+2*11+5*25+11)  ! filter and
                                                           ! mid-band data
     CHARACTER (LEN=3) :: DACS(DACchans)
     CHARACTER (LEN=1) :: spare1
     CHARACTER (LEN=80) :: laser_lo
     CHARACTER (LEN=2) :: spare2
     CHARACTER (LEN=9) :: ghz_ant_scan
     CHARACTER (LEN=9) :: thz_switch
     CHARACTER (LEN=9) :: ghz_switch
     CHARACTER (LEN=1) :: sw_decom_fmt
     CHARACTER (LEN=1) :: eng_decom_fmt
     CHARACTER (LEN=10) :: eng_data
     CHARACTER (LEN=1) :: sw_data
     CHARACTER (LEN=2) :: checksum
  END TYPE Sci2_uc_T

  TYPE (Sci1_uc_T), TARGET :: sci1_uc
  TYPE (Sci2_uc_T), TARGET :: sci2_uc

!! Formats to convert uncompressed raw Science packets internally:

  CHARACTER (LEN=*), PARAMETER :: Sci1_uc_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A1,A2,9(25A2),11A2,2(25A2),A7,2(2A1),4A1,&
       & 2(2A1),6A1,5A1,3A1,129A3,2A1,A2,2A1,A2)"

  CHARACTER (LEN=*), PARAMETER :: Sci2_uc_fmt = &
       "(3A2,3A1,A4,A2,2A1,2A2,A3,25A2,11A2,2(25A2),2(11A2),5(25A2),11A2,&
       & 129A3,A1,A80,A2,3A9,2A1,A10,A1,A2)"

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

!! Generic pointers (used after reading data and determining data format)

  TYPE Sci_cptr_T

     TYPE (CS2PTR) :: MAF(2)        ! MAF pointers
     TYPE (CS1PTR) :: MIF(2)        ! MIF pointers
     TYPE (CS2PTR) :: orbit(2)      ! orbit pointers

     TYPE (C2PTR) :: FB(FBnum)   ! Filter bank pointers
     TYPE (C2PTR) :: MB(MBnum)   ! Mid-Band pointers

     CHARACTER (LEN=9), POINTER :: GHz_ant_scan
     CHARACTER (LEN=9), POINTER :: GHz_sw
     CHARACTER (LEN=9), POINTER :: THz_sw

  END TYPE Sci_cptr_T

  TYPE (Sci_cptr_T) :: Sci_cptr

!! Uncompressed science data pointers:

!! MAF, MIF, orbit:

  TYPE (CS2PTR) :: MAF_uc_cptr(2)        ! MAF pointers
  TYPE (CS1PTR) :: MIF_uc_cptr(2)        ! MIF pointers
  TYPE (CS2PTR) :: orbit_uc_cptr(2)      ! orbit pointers

!! GHz & Thz switching data:

  CHARACTER (LEN=9), POINTER :: GHz_ant_scan_uc_cptr
  CHARACTER (LEN=9), POINTER :: GHz_sw_uc_cptr
  CHARACTER (LEN=9), POINTER :: THz_sw_uc_cptr

!! Filter Bank offsets and pointers within data field in uncompressed
!!  science packets:

  INTEGER, PARAMETER :: FB_uc_offset(FBnum) = (/ 26, 76, 101, 51, 151, 201, &
       & 126, 237, 262, 209, 1, 1, 176, 37, 62, 109, 134, 159, 184 /)

  TYPE (C2PTR) :: FB_uc_cptr(FBnum)

!! Mid-band offsets and pointers within data field in uncompressed
!!  science packets:

  INTEGER, PARAMETER :: MB_uc_offset(MBnum) = (/ 226, 234, 26, 87, 98 /)

  TYPE (C2PTR) :: MB_uc_cptr(MBnum)

!! Filter bank numbers within the science packets

  INTEGER, PARAMETER :: Pkt1_FB(11) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13 /)
  INTEGER, PARAMETER :: Pkt2_FB(8) = (/ 10, 11, 14, 15, 16, 17, 18, 19 /)

!! L0 science data packet for 1 MIF:

  TYPE Sci_pkt_T
     REAL(r8) :: secTAI
     INTEGER :: MAFno
     INTEGER :: MIFno
     INTEGER :: Orbit
     INTEGER :: FB(FBchans,FBnum)
     INTEGER :: MB(MBchans,MBnum)
     INTEGER :: WF(WFchans,WFnum)
     REAL :: GHz_sw_angle(2)
     REAL :: THz_sw_angle(2)
     REAL :: scan_angle(2)
     CHARACTER(len=2) :: GHz_sw_pos
     CHARACTER(len=2) :: THz_sw_pos
     LOGICAL :: CRC_good
  END TYPE Sci_pkt_T

!! L0 science packet for 1 MIF:

  TYPE (Sci_pkt_T) :: Sci_pkt

!! L0 science for 1 MAF:

  TYPE (Sci_pkt_T) :: SciMAF(0:(MaxMIFs-1))

CONTAINS

  SUBROUTINE InitSciPointers

    INTEGER :: i

    !! Initialize Science packet type pointers

    type_sci1 => l0_sci1(16:16)
    type_sci2 => l0_sci2(16:16)

    !! Initialize pointers for uncompressed science packets:

    !! MAF, MIF, orbit:

    MAF_uc_cptr(1)%ptr => sci1_uc%MAF
    MAF_uc_cptr(2)%ptr => sci2_uc%MAF
    MIF_uc_cptr(1)%ptr => sci1_uc%MIF
    MIF_uc_cptr(2)%ptr => sci2_uc%MIF
    orbit_uc_cptr(1)%ptr => sci1_uc%orbit
    orbit_uc_cptr(2)%ptr => sci2_uc%orbit

    !! Filter Bank pointers

    !!  Science packet #1

    DO i = 1, SIZE (Pkt1_FB)

       FB_uc_cptr(Pkt1_FB(i))%ptr => sci1_uc%FB_dat(FB_uc_offset(Pkt1_FB(i)): &
            & FB_uc_offset(Pkt1_FB(i))+24)

    ENDDO

    !!  Science packet #2

    DO i = 1, SIZE (Pkt2_FB)

       FB_uc_cptr(Pkt2_FB(i))%ptr => sci2_uc%FB_dat(FB_uc_offset(Pkt2_FB(i)): &
            & FB_uc_offset(Pkt2_FB(i))+24)

    ENDDO

    !! Mid-Band pointers

    MB_uc_cptr(1)%ptr => sci1_uc%FB_dat(MB_uc_offset(1):MB_uc_offset(1)+10)
    DO i = 2, 5
       MB_uc_cptr(i)%ptr => sci2_uc%FB_dat(MB_uc_offset(i):MB_uc_offset(i)+10)
    ENDDO

    GHz_ant_scan_uc_cptr => sci2_uc%ghz_ant_scan   ! GHz antenna/scan data
    GHz_sw_uc_cptr => sci2_uc%ghz_switch   ! GHz switching data
    THz_sw_uc_cptr => sci2_uc%thz_switch   ! THz switching data

  END SUBROUTINE InitSciPointers

END MODULE L0_sci_tbls

! $Log$
! Revision 2.1  2001/02/23 20:47:49  perun
! Version 0.5 commit
!
