! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


MODULE OutputL1B

  USE Hdf
  USE MLSCommon
  USE MLSL1Common
  USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz
  USE MLSL1Rad
  USE MLSMessageModule, ONLY : MLSMessage, MLSMSG_Error, MLSMSG_Warning

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: L1BOAINDEX_T, L1BOASC_T, L1BOATP_T, &
    & OUTPUTL1B_CREATE, OUTPUTL1B_INDEX, OUTPUTL1B_SC, &
    & OUTPUTL1B_GHZ, OUTPUTL1B_THZ, OUTPUTL1B_RAD, &
    & LENG, LENT, LENCOORD

  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &                                                    
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  ! This module contains the subroutines needed to write L1B data to
  ! SDS-HDF files.

  ! Parameters

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS1_NAME = 'MAFStartTimeUTC'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS2_NAME = 'MAFStartTimeTAI'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS3_NAME = 'noMIFs'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS4_NAME = 'scECI'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS5_NAME = 'scECR'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS6_NAME = 'scGeocAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS7_NAME = 'scGeocLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS8_NAME = 'scGeodAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS9_NAME = 'scGeodLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS10_NAME = 'scLon'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS11_NAME = 'scGeodAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS12_NAME = 'scVelECI'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS13_NAME = 'ypr'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS14_NAME = 'yprRate'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS15_NAME = 'GHz.encoderAngle'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS16_NAME = 'GHz.scAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS17_NAME = 'GHz.scanAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS18_NAME = 'GHz.scanRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS19_NAME = 'GHz.tpECI'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS20_NAME = 'GHz.tpECR'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS21_NAME = 'GHz.tpOrbY'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS22_NAME = 'GHz.tpGeocAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS23_NAME = 'GHz.tpGeocLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS24_NAME = 'GHz.tpGeocAltRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS25_NAME = 'GHz.tpGeodAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS26_NAME = 'GHz.tpGeodLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS27_NAME = 'GHz.tpGeodAltRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS28_NAME = 'GHz.tpLon'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS29_NAME = 'GHz.tpGeodAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS30_NAME = 'GHz.tpSolarTime'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS31_NAME = 'GHz.tpSolarZenith'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS32_NAME = 'GHz.tpLosAngle'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS33_NAME = 'THz.encoderAngle'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS34_NAME = 'THz.scAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS35_NAME = 'THz.scanAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS36_NAME = 'THz.scanRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS37_NAME = 'THz.tpECI'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS38_NAME = 'THz.tpECR'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS39_NAME = 'THz.tpOrbY'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS40_NAME = 'THz.tpGeocAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS41_NAME = 'THz.tpGeocLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS42_NAME = 'THz.tpGeocAltRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS43_NAME = 'THz.tpGeodAlt'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS44_NAME = 'THz.tpGeodLat'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS45_NAME = 'THz.tpGeodAltRate'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS46_NAME = 'THz.tpLon'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS47_NAME = 'THz.tpGeodAngle'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS48_NAME = 'THz.tpSolarTime'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS49_NAME = 'THz.tpSolarZenith'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS50_NAME = 'THz.tpLosAngle'

  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS51_NAME = 'counterMAF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS62_NAME = 'scVelECR'
  CHARACTER(len=*), PUBLIC, PARAMETER :: SDS63_NAME = 'scOrbIncl'

  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM1_NAME = 'MAF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM2_NAME = 'MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM3_NAME = 'GHz.MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM4_NAME = 'THz.MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM5_NAME = 'xyz'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM6_NAME = 'charUTC'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM7_NAME = 'chanFB'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM8_NAME = 'chanMB'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM9_NAME = 'chanWF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM10_NAME = 'chanDACS'

  INTEGER, PARAMETER :: lenCoord = 3
  INTEGER, PUBLIC, PARAMETER :: lenUTC = 27
  INTEGER, PARAMETER :: lenG = 120
  INTEGER, PARAMETER :: lenT = 120

  REAL, PUBLIC, PARAMETER :: FILL_REAL = -999.9
  REAL(r8), PUBLIC, PARAMETER :: FILL_DP = -999.9

  ! This data type contains index information for the L1BOA data file.
  TYPE L1BOAindex_T
    CHARACTER(LEN=lenUTC) :: MAFStartTimeUTC ! MAF start time in ascii UTC
    REAL(r8) :: MAFStartTimeTAI		! MAF start time in TAI
    INTEGER :: noMIFs                   ! total number of MIFs per MAF
    INTEGER :: counterMAF               ! MAF counter
  END TYPE L1BOAindex_T

  ! This data type contains spacecraft quantities for the L1BOA data file.
  TYPE L1BOAsc_T
    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scECI	! s/c ECI location
    REAL(r8), DIMENSION(:,:), POINTER :: scECR	! s/c ECR location

    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:), POINTER :: scGeocAlt	! s/c geoc alt
    REAL, DIMENSION(:), POINTER :: scGeocLat		! s/c geoc lat
    REAL(r8), DIMENSION(:), POINTER :: scGeodAlt	! s/c geod alt
    REAL, DIMENSION(:), POINTER :: scGeodLat		! s/c geod lat
    REAL, DIMENSION(:), POINTER :: scLon		! s/c longitude
    REAL, DIMENSION(:), POINTER :: scGeodAngle	! s/c master coordinate
    REAL, DIMENSION(:), POINTER :: scOrbIncl     ! instant. orb incl ECR (deg)

    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECI	! s/c ECI velocity
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECR     ! s/c ECR velocity
    REAL(r8), DIMENSION(:,:), POINTER :: ypr		! s/c yaw, pitch, roll
    REAL(r8), DIMENSION(:,:), POINTER :: yprRate	! s/c y-p-r rate
  END TYPE L1BOAsc_T

  ! This data type contains tangent point quantities for the L1BOA data file.
  TYPE L1BOAtp_T

    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:), POINTER :: encoderAngle	! boresight wrt instr.
    REAL(r8), DIMENSION(:), POINTER :: scAngle	! boresight wrt s/c +x
    REAL(r8), DIMENSION(:), POINTER :: scanAngle	! boresight wrt orbit +x
    REAL, DIMENSION(:), POINTER :: scanRate		! of change of scanAngle

    ! dimensioned (xyz,mod.MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: tpECI	! tp location in ECI
    REAL(r8), DIMENSION(:,:), POINTER :: tpECR	! tp location in ECR

    ! dimensioned (mod.MIF)
    REAL, DIMENSION(:), POINTER :: tpOrbY		! out-of-plane distance
    REAL(r8), DIMENSION(:), POINTER :: tpGeocAlt	! geoc alt of tp
    REAL, DIMENSION(:), POINTER :: tpGeocLat		! geoc lat of tp
    REAL, DIMENSION(:), POINTER :: tpGeocAltRate	! of change of tpGeocAlt
    REAL(r8), DIMENSION(:), POINTER :: tpGeodAlt	! geod alt of tp 
    REAL, DIMENSION(:), POINTER :: tpGeodLat		! geod lat of tp
    REAL, DIMENSION(:), POINTER :: tpGeodAltRate	! of change of tpGeodAlt
    REAL, DIMENSION(:), POINTER :: tpLon		! longitude of tp
    REAL, DIMENSION(:), POINTER :: tpGeodAngle	! tp master coordinate
    REAL, DIMENSION(:), POINTER :: tpSolarTime	! solar time coordinate
    REAL, DIMENSION(:), POINTER :: tpSolarZenith	! solar zenith angle
    REAL, DIMENSION(:), POINTER :: tpLosAngle		! line-of-sight ang to N
  END TYPE L1BOAtp_T

CONTAINS

  !----------------------------------- OutputL1B_Create ----
  SUBROUTINE OutputL1B_create(sdId)
    ! This subroutine opens/creates the SD output files, and names the arrays and
    ! dimensions contained within them.

    ! Arguments
    TYPE( L1BFileInfo_T ), INTENT(IN) :: sdId

    ! Functions
    INTEGER :: sfcreate, sfdimid, sfsdmname, sfendacc, sfsfill

    ! Variables
    INTEGER :: dim_id, rank, status
    INTEGER :: sds1_id, sds2_id, sds3_id, sds4_id, sds5_id, sds6_id, sds7_id
    INTEGER :: sds8_id, sds9_id, sds10_id, sds11_id, sds12_id, sds13_id
    INTEGER :: sds14_id, sds15_id, sds16_id, sds17_id, sds18_id, sds19_id
    INTEGER :: sds20_id, sds21_id, sds22_id, sds23_id, sds24_id, sds25_id
    INTEGER :: sds26_id, sds27_id, sds28_id, sds29_id, sds30_id, sds31_id
    INTEGER :: sds32_id, sds33_id, sds34_id, sds35_id, sds36_id, sds37_id
    INTEGER :: sds38_id, sds39_id, sds40_id, sds41_id, sds42_id, sds43_id
    INTEGER :: sds44_id, sds45_id, sds46_id, sds47_id, sds48_id, sds49_id
    INTEGER :: sds50_id, sds51_id, sds52_id, sds53_id
    INTEGER :: sds62_id, sds63_id
    INTEGER :: dimSize(3)

    ! Create the counterMAF SD in rad file F; name the dimension; terminate access
    ! to the data set.
    rank = 1
    dimSize(2) = SD_UNLIMITED

    sds52_id = sfcreate(sdId%RADFID, SDS51_NAME, DFNT_INT32, rank, &
      dimSize(2))
    dim_id = sfdimid(sds52_id, 0)
    status = sfsdmname(dim_id, DIM1_NAME)
    status = sfendacc(sds52_id)

    ! Repeat these steps for the Rad D file
    sds53_id = sfcreate(sdId%RADDID, SDS51_NAME, DFNT_INT32, rank, &
      dimSize(2))
    dim_id = sfdimid(sds53_id, 0)
    status = sfsdmname(dim_id, DIM1_NAME)
    status = sfendacc(sds53_id)

    ! L1BOA file -- create one-dimensional data sets
    sds2_id = sfcreate( sdId%OAId, SDS2_NAME, DFNT_FLOAT64, rank, &
      dimSize(2) )
    sds3_id = sfcreate( sdId%OAId, SDS3_NAME, DFNT_INT32, rank, &
      dimSize(2) )
    sds51_id = sfcreate( sdId%OAId, SDS51_NAME, DFNT_INT32, rank, &
      dimSize(2) )

    ! Create two-dimensional data sets
    rank = 2

    ! S/C, GHz & THz non-tp quantities
    dimSize(1) = lenUTC

    sds1_id = sfcreate( sdId%OAId, SDS1_NAME, DFNT_CHAR8, rank, &
      dimSize(1:2) )

    dimSize(1) = MaxMIFs

    sds6_id = sfcreate( sdId%OAId, SDS6_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds6_id, FILL_DP)
    sds7_id = sfcreate( sdId%OAId, SDS7_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds7_id, FILL_REAL)
    sds8_id = sfcreate( sdId%OAId, SDS8_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds8_id, FILL_DP)
    sds9_id = sfcreate( sdId%OAId, SDS9_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds9_id, FILL_REAL)
    sds10_id = sfcreate( sdId%OAId, SDS10_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds10_id, FILL_REAL)
    sds11_id = sfcreate( sdId%OAId, SDS11_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds11_id, FILL_REAL)
    sds63_id = sfcreate( sdId%OAId, SDS63_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds63_id, FILL_REAL)

    sds15_id = sfcreate( sdId%OAId, SDS15_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds15_id, FILL_DP)
    sds16_id = sfcreate( sdId%OAId, SDS16_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds16_id, FILL_DP)
    sds17_id = sfcreate( sdId%OAId, SDS17_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds17_id, FILL_DP)
    sds18_id = sfcreate( sdId%OAId, SDS18_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds18_id, FILL_REAL)

    sds33_id = sfcreate( sdId%OAId, SDS33_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds33_id, FILL_DP)
    sds34_id = sfcreate( sdId%OAId, SDS34_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds34_id, FILL_DP)
    sds35_id = sfcreate( sdId%OAId, SDS35_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds35_id, FILL_DP)
    sds36_id = sfcreate( sdId%OAId, SDS36_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds36_id, FILL_REAL)

    ! GHz tangent point quantities
    dimSize(1) = MIFsGHz

    sds21_id = sfcreate( sdId%OAId, SDS21_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds21_id, FILL_REAL)
    sds22_id = sfcreate( sdId%OAId, SDS22_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds22_id, FILL_DP)
    sds23_id = sfcreate( sdId%OAId, SDS23_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds23_id, FILL_REAL)
    sds24_id = sfcreate( sdId%OAId, SDS24_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds24_id, FILL_REAL)
    sds25_id = sfcreate( sdId%OAId, SDS25_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds25_id, FILL_DP)
    sds26_id = sfcreate( sdId%OAId, SDS26_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds26_id, FILL_REAL)
    sds27_id = sfcreate( sdId%OAId, SDS27_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds27_id, FILL_REAL)
    sds28_id = sfcreate( sdId%OAId, SDS28_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds28_id, FILL_REAL)
    sds29_id = sfcreate( sdId%OAId, SDS29_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) ) 
    status = sfsfill(sds29_id, FILL_REAL)
    sds30_id = sfcreate( sdId%OAId, SDS30_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds30_id, FILL_REAL)
    sds31_id = sfcreate( sdId%OAId, SDS31_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds31_id, FILL_REAL)
    sds32_id = sfcreate( sdId%OAId, SDS32_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds32_id, FILL_REAL)

    ! THz tangent point quantities

    dimSize(1) = MIFsTHz

    sds39_id = sfcreate( sdId%OAId, SDS39_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds39_id, FILL_REAL)
    sds40_id = sfcreate( sdId%OAId, SDS40_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds40_id, FILL_DP)
    sds41_id = sfcreate( sdId%OAId, SDS41_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds41_id, FILL_REAL)
    sds42_id = sfcreate( sdId%OAId, SDS42_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds42_id, FILL_REAL)
    sds43_id = sfcreate( sdId%OAId, SDS43_NAME, DFNT_FLOAT64, rank, &
      dimSize(1:2) )
    status = sfsfill(sds43_id, FILL_DP)
    sds44_id = sfcreate( sdId%OAId, SDS44_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2))
    status = sfsfill(sds44_id, FILL_REAL)
    sds45_id = sfcreate( sdId%OAId, SDS45_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds45_id, FILL_REAL)
    sds46_id = sfcreate( sdId%OAId, SDS46_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds46_id, FILL_REAL)
    sds47_id = sfcreate( sdId%OAId, SDS47_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds47_id, FILL_REAL)
    sds48_id = sfcreate( sdId%OAId, SDS48_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds48_id, FILL_REAL)
    sds49_id = sfcreate( sdId%OAId, SDS49_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2) )
    status = sfsfill(sds49_id, FILL_REAL)
    sds50_id = sfcreate( sdId%OAId, SDS50_NAME, DFNT_FLOAT32, rank, &
      dimSize(1:2))
    status = sfsfill(sds50_id, FILL_REAL)

    ! Create three-dimensional data sets
    rank = 3

    ! S/C quantities
    dimSize(1) = lenCoord
    dimSize(2) = MaxMIFs
    dimSize(3) = SD_UNLIMITED

    sds4_id = sfcreate(sdId%OAId, SDS4_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds4_id, FILL_DP)
    sds5_id = sfcreate(sdId%OAId, SDS5_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds5_id, FILL_DP)
    sds12_id = sfcreate(sdId%OAId, SDS12_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds12_id, FILL_DP)
    sds62_id = sfcreate(sdId%OAId, SDS62_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds62_id, FILL_DP)
    sds13_id = sfcreate(sdId%OAId, SDS13_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds13_id, FILL_DP)
    sds14_id = sfcreate(sdId%OAId, SDS14_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds14_id, FILL_DP)

    ! GHz tangent point quantities
    dimSize(2) = MIFsGHz

    sds19_id = sfcreate(sdId%OAId, SDS19_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds19_id, FILL_DP)
    sds20_id = sfcreate(sdId%OAId, SDS20_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds20_id, FILL_DP)

    ! THz tangent point quantities
    dimSize(2) = MIFsTHz

    sds37_id = sfcreate(sdId%OAId, SDS37_NAME, DFNT_FLOAT64, rank, dimSize) 
    status = sfsfill(sds37_id, FILL_DP)
    sds38_id = sfcreate(sdId%OAId, SDS38_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds38_id, FILL_DP)

    ! Give names to the dimensions
    dim_id = sfdimid(sds1_id, 0)
    status = sfsdmname(dim_id, DIM6_NAME)
    dim_id = sfdimid(sds1_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds2_id, 0)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds3_id, 0)
    status = sfsdmname(dim_id, DIM1_NAME)

    ! S/C quantities
    dim_id = sfdimid(sds4_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds4_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds4_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds5_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds5_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds5_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds6_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds6_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds7_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds7_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds8_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds8_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds9_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds9_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds10_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds10_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds63_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds63_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds11_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds11_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds12_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds12_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds12_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds62_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds62_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds62_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds13_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds13_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds13_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds14_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds14_id, 1)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds14_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    ! GHz quantities
    dim_id = sfdimid(sds15_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds15_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds16_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds16_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds17_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds17_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds18_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds18_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds19_id, 0)
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds19_id, 1)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds19_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds20_id, 0) 
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds20_id, 1)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds20_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds21_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds21_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds22_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds22_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds23_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds23_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds24_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds24_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds25_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds25_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds26_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds26_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds27_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds27_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds28_id, 0) 
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds28_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds29_id, 0)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds29_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds30_id, 0)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds30_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds31_id, 0)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds31_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds32_id, 0)
    status = sfsdmname(dim_id, DIM3_NAME)
    dim_id = sfdimid(sds32_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    ! THz quantities
    dim_id = sfdimid(sds33_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds33_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds34_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds34_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds35_id, 0) 
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds35_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds36_id, 0)
    status = sfsdmname(dim_id, DIM2_NAME)
    dim_id = sfdimid(sds36_id, 1) 
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds37_id, 0) 
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds37_id, 1)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds37_id, 2) 
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds38_id, 0) 
    status = sfsdmname(dim_id, DIM5_NAME)
    dim_id = sfdimid(sds38_id, 1)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds38_id, 2)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds39_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds39_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds40_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds40_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds41_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds41_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds42_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds42_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds43_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds43_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds44_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds44_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds45_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds45_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds46_id, 0) 
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds46_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds47_id, 0)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds47_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds48_id, 0)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds48_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds49_id, 0)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds49_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    dim_id = sfdimid(sds50_id, 0)
    status = sfsdmname(dim_id, DIM4_NAME)
    dim_id = sfdimid(sds50_id, 1)
    status = sfsdmname(dim_id, DIM1_NAME)

    ! MAF counter
    dim_id = sfdimid(sds51_id, 0)
    status = sfsdmname(dim_id, DIM1_NAME)

    ! Terminate access to the data sets
    status = sfendacc(sds1_id)
    status = sfendacc(sds2_id)
    status = sfendacc(sds3_id)
    status = sfendacc(sds4_id)
    status = sfendacc(sds5_id)
    status = sfendacc(sds6_id)
    status = sfendacc(sds7_id)
    status = sfendacc(sds8_id)
    status = sfendacc(sds9_id)
    status = sfendacc(sds10_id)
    status = sfendacc(sds11_id)
    status = sfendacc(sds63_id)
    status = sfendacc(sds12_id)
    status = sfendacc(sds62_id)
    status = sfendacc(sds13_id)
    status = sfendacc(sds14_id)
    status = sfendacc(sds15_id)
    status = sfendacc(sds16_id)
    status = sfendacc(sds17_id)
    status = sfendacc(sds18_id)
    status = sfendacc(sds19_id)
    status = sfendacc(sds20_id)
    status = sfendacc(sds21_id)
    status = sfendacc(sds22_id)
    status = sfendacc(sds23_id)
    status = sfendacc(sds24_id)
    status = sfendacc(sds25_id)
    status = sfendacc(sds26_id)
    status = sfendacc(sds27_id)
    status = sfendacc(sds28_id)
    status = sfendacc(sds29_id)
    status = sfendacc(sds30_id)
    status = sfendacc(sds31_id)
    status = sfendacc(sds32_id)
    status = sfendacc(sds33_id)
    status = sfendacc(sds34_id)
    status = sfendacc(sds35_id)
    status = sfendacc(sds36_id)
    status = sfendacc(sds37_id)
    status = sfendacc(sds38_id)
    status = sfendacc(sds39_id)
    status = sfendacc(sds40_id)
    status = sfendacc(sds41_id)
    status = sfendacc(sds42_id)
    status = sfendacc(sds43_id)
    status = sfendacc(sds44_id)
    status = sfendacc(sds45_id)
    status = sfendacc(sds46_id)
    status = sfendacc(sds47_id)
    status = sfendacc(sds48_id)
    status = sfendacc(sds49_id)
    status = sfendacc(sds50_id)
    status = sfendacc(sds51_id)

  END SUBROUTINE OutputL1B_create

  !------------------------------------------------- OuptutL1B_index ----
  SUBROUTINE OutputL1B_index(noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF-SD file. 

    ! Arguments
    TYPE( L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: sd_id, noMAF

    ! Functions
    INTEGER :: sfn2index, sfselect, sfwdata, sfwcdata, sfendacc

    ! Variables
    INTEGER :: sds_index, sds1_id, sds2_id, sds3_id, sds51_id, status
    INTEGER :: edge(2), start(2), stride(2)

    ! Executable code

    ! Find data sets by name
    sds_index = sfn2index(sd_id, SDS1_NAME)
    sds1_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS2_NAME)
    sds2_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS3_NAME)
    sds3_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS51_NAME)
    sds51_id = sfselect(sd_id, sds_index)

    ! Initialize parameters
    stride = 1

    ! Write time, noMIFs, counterMAF
    start(1) = 0
    start(2) = noMAF-1
    edge(1) = lenUTC
    edge(2) = 1

    status = sfwcdata(sds1_id, start, stride, edge, index%MAFStartTimeUTC)
    status = sfwdata(sds2_id, start(2), stride(2), edge(2), &
      index%MAFStartTimeTAI)
    status = sfwdata(sds3_id, start(2), stride(2), edge(2), index%noMIFs)
    status = sfwdata(sds51_id, start(2), stride(2), edge(2), &
      index%counterMAF)

    ! Terminate access to the data sets
    status = sfendacc(sds1_id)
    status = sfendacc(sds2_id)
    status = sfendacc(sds3_id)
    status = sfendacc(sds51_id)

  END SUBROUTINE OutputL1B_index

  !------------------------------------------- OutputL1B_sc ------------
  SUBROUTINE OutputL1B_sc(noMAF, sd_id, sc)
    ! This subroutine writes the spacecraft quantities to the HDF-SD file.

    ! Arguments
    TYPE( L1BOAsc_T ), INTENT(IN) :: sc
    INTEGER, INTENT(IN) :: noMAF, sd_id

    ! Parameters
    INTEGER :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    INTEGER :: sds_index, status
    INTEGER :: sds4_id, sds5_id, sds6_id, sds7_id, sds8_id, sds9_id
    INTEGER :: sds10_id, sds11_id, sds12_id, sds13_id, sds14_id
    INTEGER :: sds62_id, sds63_id
    INTEGER :: edge(3), start(3), stride(3)

    ! Find data sets by name
    sds_index = sfn2index(sd_id, SDS4_NAME)
    sds4_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS5_NAME)
    sds5_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS6_NAME)
    sds6_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS7_NAME)
    sds7_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS8_NAME)
    sds8_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS9_NAME)
    sds9_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS10_NAME)
    sds10_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS11_NAME)
    sds11_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS63_NAME)
    sds63_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS12_NAME)
    sds12_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS62_NAME)
    sds62_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS13_NAME)
    sds13_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS14_NAME)
    sds14_id = sfselect(sd_id, sds_index)

    ! Initialize parameters that aren't reset
    start(1) = 0
    stride = 1
    edge(3) = 1

    ! Write 3-D slabs (xyz x MIF x MAF)
    start(2) = 0
    start(3) = noMAF-1
    edge(1) = lenCoord
    edge(2) = SIZE(sc%scECI,2)

    status = sfwdata(sds4_id, start, stride, edge, sc%scECI)
    status = sfwdata(sds5_id, start, stride, edge, sc%scECR)
    status = sfwdata(sds12_id, start, stride, edge, sc%scVelECI)
    status = sfwdata(sds62_id, start, stride, edge, sc%scVelECR)
    status = sfwdata(sds13_id, start, stride, edge, sc%ypr)
    status = sfwdata(sds14_id, start, stride, edge, sc%yprRate)

    ! Write 2-D slabs (MIF x MAF)
    start(2) = noMAF-1
    edge(1) = SIZE(sc%scGeocAlt)
    edge(2) = 1

    status = sfwdata(sds6_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scGeocAlt)
    status = sfwdata(sds7_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scGeocLat)
    status = sfwdata(sds8_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scGeodAlt)
    status = sfwdata(sds9_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scGeodLat)
    status = sfwdata(sds10_id, start(1:2), stride(1:2), edge(1:2), sc%scLon)
    status = sfwdata(sds11_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scGeodAngle)
    status = sfwdata(sds63_id, start(1:2), stride(1:2), edge(1:2), &
      sc%scOrbIncl)

    ! Terminate access to the data sets
    status = sfendacc(sds4_id)
    status = sfendacc(sds5_id)
    status = sfendacc(sds6_id)
    status = sfendacc(sds7_id)
    status = sfendacc(sds8_id)
    status = sfendacc(sds9_id)
    status = sfendacc(sds10_id)
    status = sfendacc(sds11_id)
    status = sfendacc(sds63_id)
    status = sfendacc(sds12_id)
    status = sfendacc(sds62_id)
    status = sfendacc(sds13_id)
    status = sfendacc(sds14_id)

  END SUBROUTINE OutputL1B_sc

  !-------------------------------------------- OutputL1B_GHz -------
  SUBROUTINE OutputL1B_GHz(noMAF, sd_id, tp)
    ! This subroutine writes the GHz tangent point quantities to the HDF-SD file.
    ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id

    ! Functions
    INTEGER :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    INTEGER :: sds_index, sds15_id, sds16_id, sds17_id, sds18_id
    INTEGER :: sds19_id, sds20_id, sds21_id, sds22_id, sds23_id, sds24_id
    INTEGER :: sds25_id, sds26_id, sds27_id, sds28_id, sds29_id, sds30_id
    INTEGER :: sds31_id, sds32_id, status
    INTEGER :: edge(3), start(3), stride(3)

    ! Find data sets by name
    sds_index = sfn2index(sd_id, SDS15_NAME)
    sds15_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS16_NAME)
    sds16_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS17_NAME)
    sds17_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS18_NAME)
    sds18_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS19_NAME)
    sds19_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS20_NAME)
    sds20_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS21_NAME)
    sds21_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS22_NAME)
    sds22_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS23_NAME)
    sds23_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS24_NAME)
    sds24_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS25_NAME)
    sds25_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS26_NAME)
    sds26_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS27_NAME)
    sds27_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS28_NAME)
    sds28_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS29_NAME)
    sds29_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS30_NAME)
    sds30_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS31_NAME)
    sds31_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS32_NAME)
    sds32_id = sfselect(sd_id, sds_index)

    ! Initialize parameters that aren't reset

    start(1) = 0
    stride = 1
    edge(3) = 1

    ! Write GHz data

    start(2) = noMAF-1
    edge(1) = SIZE(tp%encoderAngle)
    edge(2) = 1

    status = sfwdata(sds15_id, start(1:2), stride(1:2), edge(1:2), &
      tp%encoderAngle)
    status = sfwdata(sds16_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scAngle)
    status = sfwdata(sds17_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanAngle)
    status = sfwdata(sds18_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanRate)

    edge(1) = SIZE(tp%tpOrbY)

    status = sfwdata(sds21_id, start(1:2), stride(1:2), edge(1:2), tp%tpOrbY)
    status = sfwdata(sds22_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocAlt)
    status = sfwdata(sds23_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocLat)
    status = sfwdata(sds24_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocAltRate)
    status = sfwdata(sds25_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAlt)
    status = sfwdata(sds26_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodLat)
    status = sfwdata(sds27_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAltRate)
    status = sfwdata(sds28_id, start(1:2), stride(1:2), edge(1:2), tp%tpLon)
    status = sfwdata(sds29_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAngle)
    status = sfwdata(sds30_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpSolarTime)
    status = sfwdata(sds31_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpSolarZenith)
    status = sfwdata(sds32_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpLosAngle)

    start(2) = 0
    start(3) = noMAF-1
    edge(1) = lenCoord
    edge(2) = SIZE(tp%tpECI,2)

    status = sfwdata(sds19_id, start, stride, edge, tp%tpECI)
    status = sfwdata(sds20_id, start, stride, edge, tp%tpECR)

    ! Terminate access to the data sets
    status = sfendacc(sds15_id)
    status = sfendacc(sds16_id)
    status = sfendacc(sds17_id)
    status = sfendacc(sds18_id)
    status = sfendacc(sds19_id)
    status = sfendacc(sds20_id)
    status = sfendacc(sds21_id)
    status = sfendacc(sds22_id)
    status = sfendacc(sds23_id)
    status = sfendacc(sds24_id)
    status = sfendacc(sds25_id)
    status = sfendacc(sds26_id)
    status = sfendacc(sds27_id)
    status = sfendacc(sds28_id)
    status = sfendacc(sds29_id)
    status = sfendacc(sds30_id)
    status = sfendacc(sds31_id)
    status = sfendacc(sds32_id)

  END SUBROUTINE OutputL1B_GHz

  !-------------------------------------------- OutputL1B_THz --------------
  SUBROUTINE OutputL1B_THz(noMAF, sd_id, tp)
    ! This subroutine writes the THz tangent point quantities to the HDF-SD file.

    ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id

    ! Functions
    INTEGER :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    INTEGER :: sds_index,  sds33_id, sds34_id, sds35_id, sds36_id
    INTEGER :: sds37_id, sds38_id, sds39_id, sds40_id, sds41_id, sds42_id
    INTEGER :: sds43_id, sds44_id, sds45_id, sds46_id, sds47_id, sds48_id
    INTEGER :: sds49_id, sds50_id, status
    INTEGER :: edge(3), start(3), stride(3)

    ! Find data sets by name
    sds_index = sfn2index(sd_id, SDS33_NAME)
    sds33_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS34_NAME)
    sds34_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS35_NAME)
    sds35_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS36_NAME)
    sds36_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS37_NAME)
    sds37_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS38_NAME)
    sds38_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS39_NAME)
    sds39_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS40_NAME)
    sds40_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS41_NAME)
    sds41_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS42_NAME)
    sds42_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS43_NAME)
    sds43_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS44_NAME)
    sds44_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS45_NAME)
    sds45_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS46_NAME)
    sds46_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS47_NAME)
    sds47_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS48_NAME)
    sds48_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS49_NAME)
    sds49_id = sfselect(sd_id, sds_index)

    sds_index = sfn2index(sd_id, SDS50_NAME)
    sds50_id = sfselect(sd_id, sds_index)

    ! Initialize parameters that aren't reset during the MAF loop
    start(1) = 0
    stride = 1
    edge(3) = 1

    ! Write THz data
    start(2) = noMAF-1
    edge(1) = SIZE(tp%encoderAngle)
    edge(2) = 1

    status = sfwdata(sds33_id, start(1:2), stride(1:2), edge(1:2), &
      tp%encoderAngle)
    status = sfwdata(sds34_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scAngle)
    status = sfwdata(sds35_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanAngle)
    status = sfwdata(sds36_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanRate)

    edge(1) = SIZE(tp%tpOrbY)

    status = sfwdata(sds39_id, start(1:2), stride(1:2), edge(1:2), tp%tpOrbY)
    status = sfwdata(sds40_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocAlt)
    status = sfwdata(sds41_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocLat)
    status = sfwdata(sds42_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeocAltRate)
    status = sfwdata(sds43_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAlt)
    status = sfwdata(sds44_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodLat)
    status = sfwdata(sds45_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAltRate)
    status = sfwdata(sds46_id, start(1:2), stride(1:2), edge(1:2), tp%tpLon)
    status = sfwdata(sds47_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpGeodAngle)
    status = sfwdata(sds48_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpSolarTime)
    status = sfwdata(sds49_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpSolarZenith)
    status = sfwdata(sds50_id, start(1:2), stride(1:2), edge(1:2), &
      tp%tpLosAngle)

    start(2) = 0
    start(3) = noMAF-1
    edge(1) = lenCoord
    edge(2) = SIZE(tp%tpECI,2)

    status = sfwdata(sds37_id, start, stride, edge, tp%tpECI)
    status = sfwdata(sds38_id, start, stride, edge, tp%tpECR)

    ! Terminate access to the data sets
    status = sfendacc(sds33_id)
    status = sfendacc(sds34_id)
    status = sfendacc(sds35_id)
    status = sfendacc(sds36_id)
    status = sfendacc(sds37_id)
    status = sfendacc(sds38_id)
    status = sfendacc(sds39_id)
    status = sfendacc(sds40_id)
    status = sfendacc(sds41_id)
    status = sfendacc(sds42_id)
    status = sfendacc(sds43_id)
    status = sfendacc(sds44_id)
    status = sfendacc(sds45_id)
    status = sfendacc(sds46_id)
    status = sfendacc(sds47_id)
    status = sfendacc(sds48_id)
    status = sfendacc(sds49_id)
    status = sfendacc(sds50_id)

  END SUBROUTINE OutputL1B_THz

  !-------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad(noMAF, sdId, counterMAF, rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files

    ! Arguments
    TYPE( L1BFileInfo_T ) :: sdId
    TYPE( Radiance_T ) :: rad(:)
    INTEGER, INTENT(IN) :: counterMAF, noMAF

! Parameters

      REAL, PARAMETER :: FILL_RAD = -999.9  ! 0.0
      REAL, PARAMETER :: FILL_RAD_ERR = -1.0

    ! Functions
    INTEGER :: sfcreate, sfdimid, sfendacc, sfn2index, sfsdmname, sfsfill
    INTEGER :: sfselect, sfwdata

    ! Variables
    CHARACTER (LEN=64) :: dim_chan, dim_mif, name, prec

    INTEGER :: dim_id, i, noSDs, rank, status, sd_id, sds_index
    INTEGER :: sds1_id, sds2_id
    INTEGER :: dimSize(3), start(3), stride(3), edge(3)
    LOGICAL :: good_rad

    ! Set parameters that won't change in loop through SDs
    rank = 3
    dimSize(3) = SD_UNLIMITED

    start(1) = 0
    start(2) = 0
    start(3) = noMAF-1
    stride = 1
    edge(3) = 1

    ! Find the sds_id number for counterMAF in the file RADF, write the data to it,
    ! terminate access to the data set
    name = 'counterMAF'
    sds_index = sfn2index(sdId%RADFID, name)
    sds1_id = sfselect(sdId%RADFID, sds_index)
    status = sfwdata(sds1_id, start(3), stride(3), edge(3), counterMAF)
    status = sfendacc(sds1_id)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Error &
      &writing counterMAF to rad F file')

    ! Do the same for file RADD
    sds_index = sfn2index(sdId%RADDID, name)
    sds2_id = sfselect(sdId%RADDID, sds_index)
    status = sfwdata(sds2_id, start(3), stride(3), edge(3), counterMAF)
    status = sfendacc(sds2_id)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Error &
      &writing counterMAF to rad D file')

    ! Loop on number of SDs per MAF
    noSDs = SIZE(rad)

    DO i = 1, noSDs
      ! Concatenate SD names
      CALL GetFullMLSSignalName(rad(i)%signal, name)
      prec = TRIM(name) // ' precision'

      ! Set parameters based on input data dimensions
      dimSize(1) = SIZE(rad(i)%value,1)
      edge(1) = dimSize(1)

      ! Based on the SD name, set dim name for channel, get Id of output file
      IF ( INDEX(name,'FB') /= 0 ) THEN
        dim_chan = DIM7_NAME
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'MB') /= 0 ) THEN
        dim_chan = DIM8_NAME
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'WF') /= 0 ) THEN
        dim_chan = DIM9_NAME
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'DACS') /= 0 ) THEN
        dim_chan = DIM10_NAME
        sd_id = sdId%RADDID
      ENDIF

      ! Based on rad module, set dim name & size, edge for # of MIFs
      IF ( INDEX(name,'R5') /= 0 ) THEN
        dim_mif = DIM4_NAME
        dimSize(2) = MIFsTHz  !! lenT
        edge(2) = MIFsTHz     !! lenT
      ELSE
        dim_mif = DIM3_NAME
        dimSize(2) = MIFsGHz  !! lenG
        edge(2) = MIFsGHz     !! lenG
      ENDIF

      ! If # of input MIFs exceeds dim size, re-set & output warning msg
      IF ( SIZE(rad(i)%value,2) > edge(2) ) CALL MLSMessage(MLSMSG_Warning,&
        ModuleName, 'Number of MIFs exceeds SD size -- output truncated.')

! Determine if good radiance data exists:

         good_rad = ANY (rad(i)%value /= 0.0)

      ! Check whether the SD already exists
      sds_index = sfn2index(sd_id, name)

      ! If not, create it, and a corresponding one for precision
      IF (sds_index == -1 .AND. good_rad) THEN

        sds1_id = sfcreate(sd_id, name, DFNT_FLOAT32, rank, dimSize)
        IF (sds1_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, &
          'Error creating value SD.')
        status = sfsfill(sds1_id, FILL_RAD)

        sds2_id = sfcreate(sd_id, prec, DFNT_FLOAT32, rank, dimSize)
        IF (sds2_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, &
          'Error creating precision SD.')
        status = sfsfill(sds2_id, FILL_RAD_ERR)

        IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, &
          'Error setting fill values.')

        ! Give names to the dimensions
        dim_id = sfdimid(sds1_id, 0)
        status = sfsdmname(dim_id, dim_chan)
        dim_id = sfdimid(sds1_id, 1)
        status = sfsdmname(dim_id, dim_mif)
        dim_id = sfdimid(sds1_id, 2)
        status = sfsdmname(dim_id, DIM1_NAME)

        dim_id = sfdimid(sds2_id, 0)
        status = sfsdmname(dim_id, dim_chan)
        dim_id = sfdimid(sds2_id, 1)
        status = sfsdmname(dim_id, dim_mif)
        dim_id = sfdimid(sds2_id, 2)
        status = sfsdmname(dim_id, DIM1_NAME)

        IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, &
          'Error setting dimension names.')

        ! If the SD already exists, find the sds_id numbers for it & its precision
      ELSE
        sds1_id = sfselect(sd_id, sds_index)
        sds_index = sfn2index(sd_id, prec)
        sds2_id = sfselect(sd_id, sds_index)

      ENDIF

      ! Write data to the value & precision SDs
      IF (good_rad) THEN
         status = sfwdata(sds1_id, start, stride, edge, rad(i)%value)
         status = sfwdata(sds2_id, start, stride, edge, rad(i)%precision)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, &
           'Error writing rad data.')

      ! Terminate access to the value & precision data sets
         status = sfendacc(sds1_id)
         status = sfendacc(sds2_id)
         IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Error &
           &terminating access to the value or precision SD.')
      ENDIF

    ENDDO

  END SUBROUTINE OutputL1B_rad

END MODULE OutputL1B

! $Log$
! Revision 2.6  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.4  2001/12/06 01:03:46  pwagner
! Now writes orbit incline angle in ECR
!
! Revision 2.3  2001/12/04 00:29:16  pwagner
! Made public things needed by sids
!
! Revision 2.2  2001/10/12 22:11:05  livesey
! Tidied things up a bit, added scVelECR, but not filled yet
!
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.4  2000/02/22 15:00:16  nakamura
! Incorporated GetFullMLSSignalName subroutine to concatenate SD names from signal input.
!

