
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Sd
!===============================================================================

   USE Hdf
   USE MLSMessageModule
   USE OutputL1B
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- Sd_create
!                Sd_index
!                Sd_sc
!                Sd_GHz
!                Sd_THz

! Remarks:  This module contains the subroutines needed to write the SIDS L1BOA
!           data as an SDS-HDF file.

! Parameters

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE Sd_create(mifG, lenMIF, mifT, offsetMAF, orbPerDay, &
                        physicalFilename, postMAF, preMAF, scansPerOrb, lenMAF)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine opens/creates the SD output file, and names the arrays and
! dimensions contained within it.

! Arguments

      CHARACTER (LEN=135), INTENT(IN) :: physicalFilename

      INTEGER, INTENT(IN) :: mifG, lenMIF, mifT, offsetMAF, orbPerDay
      INTEGER, INTENT(IN) :: postMAF, preMAF, scansPerOrb

      INTEGER, INTENT(OUT) :: lenMAF

! Parameters

! Functions

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, dim_id, i, rank, sd_id, status
      INTEGER :: sds1_id, sds2_id, sds3_id, sds4_id, sds5_id, sds6_id, sds7_id
      INTEGER :: sds8_id, sds9_id, sds10_id, sds11_id, sds12_id, sds13_id
      INTEGER :: sds14_id, sds15_id, sds16_id, sds17_id, sds18_id, sds19_id
      INTEGER :: sds20_id, sds21_id, sds22_id, sds23_id, sds24_id, sds25_id
      INTEGER :: sds26_id, sds27_id, sds28_id, sds29_id, sds30_id, sds31_id
      INTEGER :: sds32_id, sds33_id, sds34_id, sds35_id, sds36_id, sds37_id
      INTEGER :: sds38_id, sds39_id, sds40_id, sds41_id, sds42_id, sds43_id
      INTEGER :: sds44_id, sds45_id, sds46_id, sds47_id, sds48_id, sds49_id
      INTEGER :: sds50_id, sds51_id
      INTEGER :: dimSize(3)
      INTEGER, ALLOCATABLE :: maf(:)

! Initializations

      lenMAF = preMAF + scansPerOrb*(orbPerDay-1) + postMAF

      ALLOCATE (maf(lenMAF), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  maf index.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      DO i = 1, lenMAF
         maf(i) = offsetMAF + (i-1)
      ENDDO

! Create the HDF file and initialize the SD interface

      sd_id = sfstart(physicalFilename, DFACC_CREATE)

! Create one-dimensional data sets

      rank = 1
      dimSize(2) = SD_UNLIMITED

      sds2_id = sfcreate( sd_id, SDS2_NAME, DFNT_FLOAT64, rank, dimSize(2) )
      sds3_id = sfcreate( sd_id, SDS3_NAME, DFNT_INT32, rank, dimSize(2) )
      sds51_id = sfcreate( sd_id, SDS51_NAME, DFNT_INT32, rank, dimSize(2) )

! Create two-dimensional data sets

      rank = 2

! S/C, GHz & THz non-tp quantities

      dimSize(1) = lenUTC
  
      sds1_id = sfcreate( sd_id, SDS1_NAME, DFNT_CHAR8, rank, dimSize(1:2) )

      dimSize(1) = lenMIF

      sds6_id = sfcreate( sd_id, SDS6_NAME, DFNT_FLOAT64, rank, dimSize(1:2) )
      sds7_id = sfcreate( sd_id, SDS7_NAME, DFNT_FLOAT32, rank, dimSize(1:2) )
      sds8_id = sfcreate( sd_id, SDS8_NAME, DFNT_FLOAT64, rank, dimSize(1:2) )
      sds9_id = sfcreate( sd_id, SDS9_NAME, DFNT_FLOAT32, rank, dimSize(1:2) )
      sds10_id = sfcreate(sd_id, SDS10_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds11_id = sfcreate(sd_id, SDS11_NAME, DFNT_FLOAT32, rank, dimSize(1:2))

      sds15_id = sfcreate(sd_id, SDS15_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds16_id = sfcreate(sd_id, SDS16_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds17_id = sfcreate(sd_id, SDS17_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds18_id = sfcreate(sd_id, SDS18_NAME, DFNT_FLOAT32, rank, dimSize(1:2))

      sds33_id = sfcreate(sd_id, SDS33_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds34_id = sfcreate(sd_id, SDS34_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds35_id = sfcreate(sd_id, SDS35_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds36_id = sfcreate(sd_id, SDS36_NAME, DFNT_FLOAT32, rank, dimSize(1:2))

! GHz tangent point quantities

      dimSize(1) = mifG

      sds21_id = sfcreate(sd_id, SDS21_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds22_id = sfcreate(sd_id, SDS22_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds23_id = sfcreate(sd_id, SDS23_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds24_id = sfcreate(sd_id, SDS24_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds25_id = sfcreate(sd_id, SDS25_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds26_id = sfcreate(sd_id, SDS26_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds27_id = sfcreate(sd_id, SDS27_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds28_id = sfcreate(sd_id, SDS28_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds29_id = sfcreate(sd_id, SDS29_NAME, DFNT_FLOAT32, rank, dimSize(1:2)) 
      sds30_id = sfcreate(sd_id, SDS30_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds31_id = sfcreate(sd_id, SDS31_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds32_id = sfcreate(sd_id, SDS32_NAME, DFNT_FLOAT32, rank, dimSize(1:2))

! THz tangent point quantities

      dimSize(1) = mifT

      sds39_id = sfcreate(sd_id, SDS39_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds40_id = sfcreate(sd_id, SDS40_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds41_id = sfcreate(sd_id, SDS41_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds42_id = sfcreate(sd_id, SDS42_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds43_id = sfcreate(sd_id, SDS43_NAME, DFNT_FLOAT64, rank, dimSize(1:2))
      sds44_id = sfcreate(sd_id, SDS44_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds45_id = sfcreate(sd_id, SDS45_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds46_id = sfcreate(sd_id, SDS46_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds47_id = sfcreate(sd_id, SDS47_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds48_id = sfcreate(sd_id, SDS48_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds49_id = sfcreate(sd_id, SDS49_NAME, DFNT_FLOAT32, rank, dimSize(1:2))
      sds50_id = sfcreate(sd_id, SDS50_NAME, DFNT_FLOAT32, rank, dimSize(1:2))

! Create three-dimensional data sets

      rank = 3

! S/C quantities

      dimSize(1) = lenCoord
      dimSize(2) = lenMIF
      dimSize(3) = SD_UNLIMITED

      sds4_id = sfcreate(sd_id, SDS4_NAME, DFNT_FLOAT64, rank, dimSize)
      sds5_id = sfcreate(sd_id, SDS5_NAME, DFNT_FLOAT64, rank, dimSize)
      sds12_id = sfcreate(sd_id, SDS12_NAME, DFNT_FLOAT64, rank, dimSize)
      sds13_id = sfcreate(sd_id, SDS13_NAME, DFNT_FLOAT64, rank, dimSize)
      sds14_id = sfcreate(sd_id, SDS14_NAME, DFNT_FLOAT64, rank, dimSize)

! GHz tangent point quantities

      dimSize(2) = mifG

      sds19_id = sfcreate(sd_id, SDS19_NAME, DFNT_FLOAT64, rank, dimSize)
      sds20_id = sfcreate(sd_id, SDS20_NAME, DFNT_FLOAT64, rank, dimSize)

! THz tangent point quantities

      dimSize(2) = mifT

      sds37_id = sfcreate(sd_id, SDS37_NAME, DFNT_FLOAT64, rank, dimSize) 
      sds38_id = sfcreate(sd_id, SDS38_NAME, DFNT_FLOAT64, rank, dimSize)
   
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

! Terminate access to the SD interface and close the file

      status = sfend(sd_id)

!--------------------------
   END SUBROUTINE Sd_create
!--------------------------

!------------------------------------------
   SUBROUTINE Sd_index(noMAF, sd_id, index)
!------------------------------------------

! Brief description of subroutine
! This subroutine writes the time/MIF indexing quantities to the HDF-SD file. 

! Arguments

      TYPE( L1BOAindex_T), INTENT(IN) :: index

      INTEGER, INTENT(IN) :: sd_id, noMAF

! Parameters

! Functions

      INTEGER :: sfn2index, sfselect, sfsflmd, sfwdata, sfwcdata

! Variables

      INTEGER :: sds_index, sds1_id, sds2_id, sds3_id, sds51_id, status
      INTEGER :: edge(2), start(2), stride(2)

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

      status = sfsflmd(sd_id, SD_NOFILL)

      stride = 1

! Write time, noMIFs

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

!-------------------------
   END SUBROUTINE Sd_index
!-------------------------

!----------------------------------------
   SUBROUTINE Sd_sc(noMAF, nV, sd_id, sc)
!----------------------------------------

! Brief description of subroutine
! This subroutine writes the spacecraft quantities to the HDF-SD file.

! Arguments

      TYPE( L1BOAsc_T ), INTENT(IN) :: sc

      INTEGER, INTENT(IN) :: noMAF, nV, sd_id

! Parameters

! Functions

      INTEGER :: sfn2index, sfselect, sfsflmd, sfsfill, sfwdata

! Variables

      INTEGER :: sds_index, status
      INTEGER :: sds4_id, sds5_id, sds6_id, sds7_id, sds8_id, sds9_id
      INTEGER :: sds10_id, sds11_id, sds12_id, sds13_id, sds14_id
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
      edge(2) = nV
      status = sfsflmd(sd_id, SD_FILL)

      status = sfsfill(sds4_id, FILL_DP)
      status = sfwdata(sds4_id, start, stride, edge, sc%scECI)

      status = sfsfill(sds5_id, FILL_DP)
      status = sfwdata(sds5_id, start, stride, edge, sc%scECR)

      status = sfsfill(sds12_id, FILL_DP)
      status = sfwdata(sds12_id, start, stride, edge, sc%scVel)

      status = sfsfill(sds13_id, FILL_DP)
      status = sfwdata(sds13_id, start, stride, edge, sc%ypr)

      status = sfsfill(sds14_id, FILL_DP)
      status = sfwdata(sds14_id, start, stride, edge, sc%yprRate)

! Write 2-D slabs (MIF x MAF)

      start(2) = noMAF-1
      edge(1) = nV
      edge(2) = 1

      status = sfsfill(sds6_id, FILL_DP)
      status = sfwdata(sds6_id, start(1:2), stride(1:2), edge(1:2), &
                       sc%scGeocAlt)

      status = sfsfill(sds7_id, FILL_REAL)
      status = sfwdata(sds7_id, start(1:2), stride(1:2), edge(1:2), &
                       sc%scGeocLat)

      status = sfsfill(sds8_id, FILL_DP)
      status = sfwdata(sds8_id, start(1:2), stride(1:2), edge(1:2), &
                       sc%scGeodAlt)

      status = sfsfill(sds9_id, FILL_REAL)
      status = sfwdata(sds9_id, start(1:2), stride(1:2), edge(1:2), &
                       sc%scGeodLat)

      status = sfsfill(sds10_id, FILL_REAL)
      status = sfwdata(sds10_id, start(1:2), stride(1:2), edge(1:2), sc%scLon)

      status = sfsfill(sds11_id, FILL_REAL)
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

!----------------------
   END SUBROUTINE Sd_sc
!----------------------

!-----------------------------------------------
   SUBROUTINE Sd_GHz(mifG, noMAF, nV, sd_id, tp)
!-----------------------------------------------

! Brief description of subroutine
! This subroutine writes the GHz tangent point quantities to the HDF-SD file.

! Arguments

      TYPE( L1BOAtp_T ), INTENT(IN) :: tp

      INTEGER, INTENT(IN) :: mifG, noMAF, nV, sd_id

! Parameters

! Functions

      INTEGER :: sfn2index, sfselect, sfsflmd, sfsfill, sfwdata

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
      edge(1) = nV
      edge(2) = 1

      status = sfsfill(sds15_id, FILL_DP)
      status = sfwdata(sds15_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%encoderAngle) 

      status = sfsfill(sds16_id, FILL_DP)
      status = sfwdata(sds16_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scAngle)

      status = sfsfill(sds17_id, FILL_DP)
      status = sfwdata(sds17_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scanAngle)

      status = sfsfill(sds18_id, FILL_REAL)
      status = sfwdata(sds18_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scanRate)

      edge(1) = mifG
      status = sfsflmd(sd_id, SD_NOFILL)

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
      edge(2) = mifG

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

!-----------------------
   END SUBROUTINE Sd_GHz
!-----------------------

!-----------------------------------------------
   SUBROUTINE Sd_THz(mifT, noMAF, nV, sd_id, tp)
!-----------------------------------------------

! Brief description of subroutine
! This subroutine writes the THz tangent point quantities to the HDF-SD file.

! Arguments

      TYPE( L1BOAtp_T ), INTENT(IN) :: tp

      INTEGER, INTENT(IN) :: mifT, noMAF, nV, sd_id

! Parameters

! Functions

      INTEGER :: sfn2index, sfselect, sfsflmd, sfsfill, sfwdata

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
      edge(1) = nV
      edge(2) = 1
      status = sfsflmd(sd_id, SD_FILL)

      status = sfsfill(sds33_id, FILL_DP)
      status = sfwdata(sds33_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%encoderAngle)

      status = sfsfill(sds34_id, FILL_DP)
      status = sfwdata(sds34_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scAngle)

      status = sfsfill(sds35_id, FILL_DP)
      status = sfwdata(sds35_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scanAngle)

      status = sfsfill(sds36_id, FILL_REAL)
      status = sfwdata(sds36_id, start(1:2), stride(1:2), edge(1:2), &
                       tp%scanRate)

      edge(1) = mifT
      status = sfsflmd(sd_id, SD_NOFILL)

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
      edge(2) = mifT

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

!-----------------------
   END SUBROUTINE Sd_THz
!-----------------------

!===============================================================================
END MODULE Sd
!===============================================================================

!# $Log$
!#
