! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


module OutputL1B

  use Hdf
  use MLSCommon
  use MLSL1Common
  use MLSL1Config, only: MIFsGHz, MIFsTHz
  use MLSL1Rad
  use MLSMessageModule, only : MLSMessage

  implicit none
  private

  public :: L1BOAINDEX_T, L1BOASC_T, L1BOATP_T, &
    & OUTPUTL1B_CREATE, OUTPUTL1B_INDEX, OUTPUTL1B_SC, &
    & OUTPUTL1B_GHZ, OUTPUTL1B_THZ, OUTPUTL1B_RAD, &
    & LENG, LENT, LENCOORD

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &                                                    
    "$Id$"
  character (LEN=*), parameter :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  ! This module contains the subroutines needed to write L1B data to
  ! SDS-HDF files.

  ! Parameters

  character(len=*), public, parameter :: SDS1_NAME = 'MAFStartTimeUTC'
  character(len=*), public, parameter :: SDS2_NAME = 'MAFStartTimeTAI'
  character(len=*), public, parameter :: SDS3_NAME = 'noMIFs'

  character(len=*), public, parameter :: SDS4_NAME = 'scECI'
  character(len=*), public, parameter :: SDS5_NAME = 'scECR'
  character(len=*), public, parameter :: SDS6_NAME = 'scGeocAlt'
  character(len=*), public, parameter :: SDS7_NAME = 'scGeocLat'
  character(len=*), public, parameter :: SDS8_NAME = 'scGeodAlt'
  character(len=*), public, parameter :: SDS9_NAME = 'scGeodLat'
  character(len=*), public, parameter :: SDS10_NAME = 'scLon'
  character(len=*), public, parameter :: SDS11_NAME = 'scGeodAngle'
  character(len=*), public, parameter :: SDS12_NAME = 'scVel'
  character(len=*), public, parameter :: SDS13_NAME = 'ypr'
  character(len=*), public, parameter :: SDS14_NAME = 'yprRate'

  character(len=*), public, parameter :: SDS15_NAME = 'GHz.encoderAngle'

  character(len=*), public, parameter :: SDS16_NAME = 'GHz.scAngle'
  character(len=*), public, parameter :: SDS17_NAME = 'GHz.scanAngle'
  character(len=*), public, parameter :: SDS18_NAME = 'GHz.scanRate'
  character(len=*), public, parameter :: SDS19_NAME = 'GHz.tpECI'
  character(len=*), public, parameter :: SDS20_NAME = 'GHz.tpECR'
  character(len=*), public, parameter :: SDS21_NAME = 'GHz.tpOrbY'
  character(len=*), public, parameter :: SDS22_NAME = 'GHz.tpGeocAlt'
  character(len=*), public, parameter :: SDS23_NAME = 'GHz.tpGeocLat'
  character(len=*), public, parameter :: SDS24_NAME = 'GHz.tpGeocAltRate'
  character(len=*), public, parameter :: SDS25_NAME = 'GHz.tpGeodAlt'
  character(len=*), public, parameter :: SDS26_NAME = 'GHz.tpGeodLat'
  character(len=*), public, parameter :: SDS27_NAME = 'GHz.tpGeodAltRate'
  character(len=*), public, parameter :: SDS28_NAME = 'GHz.tpLon'
  character(len=*), public, parameter :: SDS29_NAME = 'GHz.tpGeodAngle'
  character(len=*), public, parameter :: SDS30_NAME = 'GHz.tpSolarTime'
  character(len=*), public, parameter :: SDS31_NAME = 'GHz.tpSolarZenith'
  character(len=*), public, parameter :: SDS32_NAME = 'GHz.tpLosAngle'

  character(len=*), public, parameter :: SDS33_NAME = 'THz.encoderAngle'

  character(len=*), public, parameter :: SDS34_NAME = 'THz.scAngle'
  character(len=*), public, parameter :: SDS35_NAME = 'THz.scanAngle'
  character(len=*), public, parameter :: SDS36_NAME = 'THz.scanRate'
  character(len=*), public, parameter :: SDS37_NAME = 'THz.tpECI'
  character(len=*), public, parameter :: SDS38_NAME = 'THz.tpECR'
  character(len=*), public, parameter :: SDS39_NAME = 'THz.tpOrbY'
  character(len=*), public, parameter :: SDS40_NAME = 'THz.tpGeocAlt'
  character(len=*), public, parameter :: SDS41_NAME = 'THz.tpGeocLat'
  character(len=*), public, parameter :: SDS42_NAME = 'THz.tpGeocAltRate'
  character(len=*), public, parameter :: SDS43_NAME = 'THz.tpGeodAlt'
  character(len=*), public, parameter :: SDS44_NAME = 'THz.tpGeodLat'
  character(len=*), public, parameter :: SDS45_NAME = 'THz.tpGeodAltRate'
  character(len=*), public, parameter :: SDS46_NAME = 'THz.tpLon'
  character(len=*), public, parameter :: SDS47_NAME = 'THz.tpGeodAngle'
  character(len=*), public, parameter :: SDS48_NAME = 'THz.tpSolarTime'
  character(len=*), public, parameter :: SDS49_NAME = 'THz.tpSolarZenith'
  character(len=*), public, parameter :: SDS50_NAME = 'THz.tpLosAngle'

  character(len=*), public, parameter :: SDS51_NAME = 'counterMAF'

  character(len=*), public, parameter :: DIM1_NAME = 'MAF'
  character(len=*), public, parameter :: DIM2_NAME = 'MIF'
  character(len=*), public, parameter :: DIM3_NAME = 'GHz.MIF'
  character(len=*), public, parameter :: DIM4_NAME = 'THz.MIF'
  character(len=*), public, parameter :: DIM5_NAME = 'xyz'
  character(len=*), public, parameter :: DIM6_NAME = 'charUTC'
  character(len=*), public, parameter :: DIM7_NAME = 'chanFB'
  character(len=*), public, parameter :: DIM8_NAME = 'chanMB'
  character(len=*), public, parameter :: DIM9_NAME = 'chanWF'
  character(len=*), public, parameter :: DIM10_NAME = 'chanDACS'

  integer, parameter :: lenCoord = 3
  integer, public, parameter :: lenUTC = 27
  integer, parameter :: lenG = 120
  integer, parameter :: lenT = 114

  real, public, parameter :: FILL_REAL = -999.9
  real(r8), public, parameter :: FILL_DP = -999.9

  ! This data type contains index information for the L1BOA data file.
  type L1BOAindex_T
    character(LEN=lenUTC) :: MAFStartTimeUTC ! MAF start time in ascii UTC
    real(r8) :: MAFStartTimeTAI		! MAF start time in TAI
    integer :: noMIFs                   ! total number of MIFs per MAF
    integer :: counterMAF               ! MAF counter
  end type L1BOAindex_T

  ! This data type contains spacecraft quantities for the L1BOA data file.
  type L1BOAsc_T
    ! dimensioned (xyz,MIF)
    real(r8), dimension(:,:), pointer :: scECI	! s/c ECI location
    real(r8), dimension(:,:), pointer :: scECR	! s/c ECR location

    ! dimensioned (MIF)
    real(r8), dimension(:), pointer :: scGeocAlt	! s/c geoc alt
    real, dimension(:), pointer :: scGeocLat		! s/c geoc lat
    real(r8), dimension(:), pointer :: scGeodAlt	! s/c geod alt
    real, dimension(:), pointer :: scGeodLat		! s/c geod lat
    real, dimension(:), pointer :: scLon		! s/c longitude
    real, dimension(:), pointer :: scGeodAngle	! s/c master coordinate

    ! dimensioned (xyz,MIF)
    real(r8), dimension(:,:), pointer :: scVelECI	! s/c ECI velocity
    real(r8), dimension(:,:), pointer :: scVelECR     ! s/c ECR velocity
    real(r8), dimension(:,:), pointer :: ypr		! s/c yaw, pitch, roll
    real(r8), dimension(:,:), pointer :: yprRate	! s/c y-p-r rate
  end type L1BOAsc_T

  ! This data type contains tangent point quantities for the L1BOA data file.
  type L1BOAtp_T

    ! dimensioned (MIF)
    real(r8), dimension(:), pointer :: encoderAngle	! boresight wrt instr.
    real(r8), dimension(:), pointer :: scAngle	! boresight wrt s/c +x
    real(r8), dimension(:), pointer :: scanAngle	! boresight wrt orbit +x
    real, dimension(:), pointer :: scanRate		! of change of scanAngle

    ! dimensioned (xyz,mod.MIF)
    real(r8), dimension(:,:), pointer :: tpECI	! tp location in ECI
    real(r8), dimension(:,:), pointer :: tpECR	! tp location in ECR

    ! dimensioned (mod.MIF)
    real, dimension(:), pointer :: tpOrbY		! out-of-plane distance
    real(r8), dimension(:), pointer :: tpGeocAlt	! geoc alt of tp
    real, dimension(:), pointer :: tpGeocLat		! geoc lat of tp
    real, dimension(:), pointer :: tpGeocAltRate	! of change of tpGeocAlt
    real(r8), dimension(:), pointer :: tpGeodAlt	! geod alt of tp 
    real, dimension(:), pointer :: tpGeodLat		! geod lat of tp
    real, dimension(:), pointer :: tpGeodAltRate	! of change of tpGeodAlt
    real, dimension(:), pointer :: tpLon		! longitude of tp
    real, dimension(:), pointer :: tpGeodAngle	! tp master coordinate
    real, dimension(:), pointer :: tpSolarTime	! solar time coordinate
    real, dimension(:), pointer :: tpSolarZenith	! solar zenith angle
    real, dimension(:), pointer :: tpLosAngle		! line-of-sight ang to N
  end type L1BOAtp_T

contains

  !----------------------------------- OutputL1B_Create ----
  subroutine OutputL1B_create(sdId)
    ! This subroutine opens/creates the SD output files, and names the arrays and
    ! dimensions contained within them.

    ! Arguments
    type( L1BFileInfo_T ), intent(IN) :: sdId

    ! Functions
    integer :: sfcreate, sfdimid, sfsdmname, sfendacc, sfsfill

    ! Variables
    character (LEN=480) :: msr
    integer :: alloc_err, dim_id, i, rank, returnStatus, status
    integer :: sds1_id, sds2_id, sds3_id, sds4_id, sds5_id, sds6_id, sds7_id
    integer :: sds8_id, sds9_id, sds10_id, sds11_id, sds12_id, sds13_id
    integer :: sds14_id, sds15_id, sds16_id, sds17_id, sds18_id, sds19_id
    integer :: sds20_id, sds21_id, sds22_id, sds23_id, sds24_id, sds25_id
    integer :: sds26_id, sds27_id, sds28_id, sds29_id, sds30_id, sds31_id
    integer :: sds32_id, sds33_id, sds34_id, sds35_id, sds36_id, sds37_id
    integer :: sds38_id, sds39_id, sds40_id, sds41_id, sds42_id, sds43_id
    integer :: sds44_id, sds45_id, sds46_id, sds47_id, sds48_id, sds49_id
    integer :: sds50_id, sds51_id, sds52_id, sds53_id
    integer :: dimSize(3)

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
    dimSize(1) = lenG

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

    dimSize(1) = lenT

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
    sds13_id = sfcreate(sdId%OAId, SDS13_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds13_id, FILL_DP)
    sds14_id = sfcreate(sdId%OAId, SDS14_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds14_id, FILL_DP)

    ! GHz tangent point quantities
    dimSize(2) = lenG

    sds19_id = sfcreate(sdId%OAId, SDS19_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds19_id, FILL_DP)
    sds20_id = sfcreate(sdId%OAId, SDS20_NAME, DFNT_FLOAT64, rank, dimSize)
    status = sfsfill(sds20_id, FILL_DP)

    ! THz tangent point quantities
    dimSize(2) = lenT

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
    status = sfendacc(sds12_id)
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

  end subroutine OutputL1B_create

  !------------------------------------------------- OuptutL1B_index ----
  subroutine OutputL1B_index(noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF-SD file. 

    ! Arguments
    type( L1BOAindex_T), intent(IN) :: index
    integer, intent(IN) :: sd_id, noMAF

    ! Functions
    integer :: sfn2index, sfselect, sfwdata, sfwcdata, sfendacc

    ! Variables
    integer :: sds_index, sds1_id, sds2_id, sds3_id, sds51_id, status
    integer :: edge(2), start(2), stride(2)

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

  end subroutine OutputL1B_index

  !------------------------------------------- OutputL1B_sc ------------
  subroutine OutputL1B_sc(noMAF, sd_id, sc)
    ! This subroutine writes the spacecraft quantities to the HDF-SD file.

    ! Arguments
    type( L1BOAsc_T ), intent(IN) :: sc
    integer, intent(IN) :: noMAF, sd_id

    ! Parameters
    integer :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    integer :: sds_index, status
    integer :: sds4_id, sds5_id, sds6_id, sds7_id, sds8_id, sds9_id
    integer :: sds10_id, sds11_id, sds12_id, sds13_id, sds14_id
    integer :: edge(3), start(3), stride(3)

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

    sds_index = sfn2index(sd_id, SDS12_NAME)
    sds12_id = sfselect(sd_id, sds_index)

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
    edge(2) = size(sc%scECI,2)

    status = sfwdata(sds4_id, start, stride, edge, sc%scECI)
    status = sfwdata(sds5_id, start, stride, edge, sc%scECR)
    status = sfwdata(sds12_id, start, stride, edge, sc%scVelECI)
    status = sfwdata(sds13_id, start, stride, edge, sc%ypr)
    status = sfwdata(sds14_id, start, stride, edge, sc%yprRate)

    ! Write 2-D slabs (MIF x MAF)
    start(2) = noMAF-1
    edge(1) = size(sc%scGeocAlt)
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

    ! Terminate access to the data sets
    status = sfendacc(sds4_id)
    status = sfendacc(sds5_id)
    status = sfendacc(sds6_id)
    status = sfendacc(sds7_id)
    status = sfendacc(sds8_id)
    status = sfendacc(sds9_id)
    status = sfendacc(sds10_id)
    status = sfendacc(sds11_id)
    status = sfendacc(sds12_id)
    status = sfendacc(sds13_id)
    status = sfendacc(sds14_id)

  end subroutine OutputL1B_sc

  !-------------------------------------------- OutputL1B_GHz -------
  subroutine OutputL1B_GHz(noMAF, sd_id, tp)
    ! This subroutine writes the GHz tangent point quantities to the HDF-SD file.
    ! Arguments
    type( L1BOAtp_T ), intent(IN) :: tp
    integer, intent(IN) :: noMAF, sd_id

    ! Functions
    integer :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    integer :: sds_index, sds15_id, sds16_id, sds17_id, sds18_id
    integer :: sds19_id, sds20_id, sds21_id, sds22_id, sds23_id, sds24_id
    integer :: sds25_id, sds26_id, sds27_id, sds28_id, sds29_id, sds30_id
    integer :: sds31_id, sds32_id, status
    integer :: edge(3), start(3), stride(3)

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
    edge(1) = size(tp%encoderAngle)
    edge(2) = 1

    status = sfwdata(sds15_id, start(1:2), stride(1:2), edge(1:2), &
      tp%encoderAngle)
    status = sfwdata(sds16_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scAngle)
    status = sfwdata(sds17_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanAngle)
    status = sfwdata(sds18_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanRate)

    edge(1) = size(tp%tpOrbY)

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
    edge(2) = size(tp%tpECI,2)

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

  end subroutine OutputL1B_GHz

  !-------------------------------------------- OutputL1B_THz --------------
  subroutine OutputL1B_THz(noMAF, sd_id, tp)
    ! This subroutine writes the THz tangent point quantities to the HDF-SD file.

    ! Arguments
    type( L1BOAtp_T ), intent(IN) :: tp
    integer, intent(IN) :: noMAF, sd_id

    ! Functions
    integer :: sfn2index, sfselect, sfwdata, sfendacc

    ! Variables
    integer :: sds_index,  sds33_id, sds34_id, sds35_id, sds36_id
    integer :: sds37_id, sds38_id, sds39_id, sds40_id, sds41_id, sds42_id
    integer :: sds43_id, sds44_id, sds45_id, sds46_id, sds47_id, sds48_id
    integer :: sds49_id, sds50_id, status
    integer :: edge(3), start(3), stride(3)

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
    edge(1) = size(tp%encoderAngle)
    edge(2) = 1

    status = sfwdata(sds33_id, start(1:2), stride(1:2), edge(1:2), &
      tp%encoderAngle)
    status = sfwdata(sds34_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scAngle)
    status = sfwdata(sds35_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanAngle)
    status = sfwdata(sds36_id, start(1:2), stride(1:2), edge(1:2), &
      tp%scanRate)

    edge(1) = size(tp%tpOrbY)

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
    edge(2) = size(tp%tpECI,2)

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

  end subroutine OutputL1B_THz

  !-------------------------------------------------------- OutputL1B_rad
  subroutine OutputL1B_rad(noMAF, sdId, counterMAF, rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files

    ! Arguments
    type( L1BFileInfo_T ) :: sdId
    type( Radiance_T ) :: rad(:)
    integer, intent(IN) :: counterMAF, noMAF

    ! Functions
    integer :: sfcreate, sfdimid, sfendacc, sfn2index, sfsdmname, sfsfill
    integer :: sfselect, sfwdata

    ! Variables
    character (LEN=64) :: dim_chan, dim_mif, name, prec

    integer :: dim_id, i, noSDs, rank, status, sd_id, sds_index
    integer :: sds1_id, sds2_id
    integer :: dimSize(3), start(3), stride(3), edge(3)

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
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Error &
      &writing counterMAF to rad F file')

    ! Do the same for file RADD
    sds_index = sfn2index(sdId%RADDID, name)
    sds2_id = sfselect(sdId%RADDID, sds_index)
    status = sfwdata(sds2_id, start(3), stride(3), edge(3), counterMAF)
    status = sfendacc(sds2_id)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Error &
      &writing counterMAF to rad D file')

    ! Loop on number of SDs per MAF
    noSDs = size(rad)

    do i = 1, noSDs
      ! Concatenate SD names
      call GetFullMLSSignalName(rad(i)%signal, name)
      prec = trim(name) // ' precision'

      ! Set parameters based on input data dimensions
      dimSize(1) = size(rad(i)%value,1)
      edge(1) = dimSize(1)

      ! Based on the SD name, set dim name for channel, get Id of output file
      if ( index(name,'FB') /= 0 ) then
        dim_chan = DIM7_NAME
        sd_id = sdId%RADFID
      else if ( index(name,'MB') /= 0 ) then
        dim_chan = DIM8_NAME
        sd_id = sdId%RADFID
      else if ( index(name,'WF') /= 0 ) then
        dim_chan = DIM9_NAME
        sd_id = sdId%RADFID
      else if ( index(name,'DACS') /= 0 ) then
        dim_chan = DIM10_NAME
        sd_id = sdId%RADDID
      endif

      ! Based on rad module, set dim name & size, edge for # of MIFs
      if ( index(name,'R5') /= 0 ) then
        dim_mif = DIM4_NAME
        dimSize(2) = MIFsTHz  !! lenT
        edge(2) = MIFsTHz     !! lenT
      else
        dim_mif = DIM3_NAME
        dimSize(2) = MIFsGHz  !! lenG
        edge(2) = MIFsGHz     !! lenG
      endif

      ! If # of input MIFs exceeds dim size, re-set & output warning msg
      if ( size(rad(i)%value,2) > edge(2) ) call MLSMessage(MLSMSG_Warning,&
        ModuleName, 'Number of MIFs exceeds SD size -- output truncated.')

      ! Check whether the SD already exists
      sds_index = sfn2index(sd_id, name)

      ! If not, create it, and a corresponding one for precision
      if ( sds_index == -1) then

        sds1_id = sfcreate(sd_id, name, DFNT_FLOAT32, rank, dimSize)
        if (sds1_id == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
          'Error creating value SD.')
        status = sfsfill(sds1_id, FILL_REAL)

        sds2_id = sfcreate(sd_id, prec, DFNT_FLOAT32, rank, dimSize)
        if (sds2_id == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
          'Error creating precision SD.')
        status = sfsfill(sds2_id, FILL_REAL)

        if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
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

        if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
          'Error setting dimension names.')

        ! If the SD already exists, find the sds_id numbers for it & its precision
      else
        sds1_id = sfselect(sd_id, sds_index)
        sds_index = sfn2index(sd_id, prec)
        sds2_id = sfselect(sd_id, sds_index)

      endif
      ! Write data to the value & precision SDs
      status = sfwdata(sds1_id, start, stride, edge, rad(i)%value)
      status = sfwdata(sds2_id, start, stride, edge, rad(i)%precision)
      if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
        'Error writing rad data.')

      ! Terminate access to the value & precision data sets
      status = sfendacc(sds1_id)
      status = sfendacc(sds2_id)
      if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Error &
        &terminating access to the value or precision SD.')

    enddo

  end subroutine OutputL1B_rad

end module OutputL1B

! $Log$
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

