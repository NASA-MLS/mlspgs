!=============================================================================
MODULE Output_UARS_L1B
!=============================================================================

  USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
  USE MLSAuxData, ONLY: Build_MLSAuxData, CreateGroup_MLSAuxData
  USE MLSSTRINGS, only: COMPRESSSTRING

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: INT_FILL = -1

  PUBLIC :: OutputL1B_rad, OutputL1B_OA

CONTAINS

  SUBROUTINE OutputL1B_rad (sdId, noMAF, counterMAF, rad, prec, gird_time)

    USE rad_file_contents, ONLY: limb_hdr

   ! This subroutine writes an MAF's worth of data to the L1BRad F file

    INTEGER, INTENT(IN) :: sdId, noMAF, counterMAF
    REAL, INTENT(IN) :: rad(15,6,32), prec(15,6,32)
    REAL*8, INTENT(IN) :: gird_time

    TYPE (DataProducts_T) :: dataset
    integer, dimension(3) :: dims = (/15, 32, 1/)
    integer :: bank, status, bits, spindx

    REAL, PARAMETER :: RAD_FILL = -999.99, RAD_ERR_FILL = -1.0

    character(len=40) :: name
    CHARACTER(len=4), PARAMETER :: species(6) = (/ &
         "PT", "ClO" , "H2O2", "O3", "H2O", "O3" /)
    CHARACTER(len=3), PARAMETER :: band(6) = (/ &
         'B1F', 'B2F', 'B3F', 'B4F', 'B5F', 'B6F' /)
    CHARACTER(len=2), PARAMETER :: radiometer(3) = (/ 'R1', 'R2', 'R3' /)
    CHARACTER(len=3), PARAMETER :: freq(6) = (/ &
         '63', '205', '205', '205', '183', '183' /)
    CHARACTER(len=2), PARAMETER :: switch(0:2) = (/ 'S0', 'S1', 'S2' /)
    integer :: mask07, bank3, bank6, swpos
    data mask07 / z'07'/  ! lower 3 bits for Filter Bank 3

    CALL Deallocate_DataProducts (dataset)

    ALLOCATE (dataset%Dimensions(1), stat=status)

    dataset%name      = 'counterMAF'
    dataset%data_type = 'integer'
    dataset%Dimensions(1) = 'MAF'
    CALL Build_MLSAuxData (sdId, dataset, counterMAF, fill_value=INT_FILL, &
         lastIndex=noMAF)

    dataset%name      = 'MAFStartTimeGIRD'
    dataset%data_type = 'double'
    CALL Build_MLSAuxData (sdId, dataset, gird_time, fill_value=-1.0d0, &
         lastIndex=noMAF)

    DEALLOCATE (dataset%Dimensions, stat=status)

    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%data_type = 'real'
    dataset%Dimensions(1) = 'chanFB'
    dataset%Dimensions(2) = 'MIF'
    dataset%Dimensions(3) = 'MAF'

    bits = ichar(limb_hdr%band_bank)

    bank3 = iand (bits, mask07)  ! what's hooked up to filter bank 3
    if (btest (bits, 3)) then    ! what's hooked up to filter bank 6
       bank6 = 3
    else
       bank6 = 6
    endif

    do bank = 1, 6
       swpos = 0
       if (bank == 3) then
          spindx = bank3
          if (bank3 /= 3) swpos = 1
       else if (bank == 6) then
          spindx = bank6
          if (bank6 /= 3) swpos = 2
       else
          spindx = bank
       endif
       if (bank == 1) then
          name = radiometer(1)
       else if (bank >= 2 .and. bank <= 4) then
          name = radiometer(2)
       else
          name = radiometer(3)
       endif

       name = name(1:len_trim(name)) // ':' // freq(bank) // '.' // &
            band(bank) // ':' // species(spindx) // '.' // switch(swpos)
       name = name(1:len_trim(name)) // '.FB15-' // char(bank+48)
       name = CompressString (name) 

       dataset%name = TRIM(name)
       CALL Build_MLSAuxData (sdId, dataset, rad(:,bank,:), &
            lastIndex=noMAF, dims=dims, fill_value=RAD_FILL)
       dataset%name = TRIM(name)//' precision'
       CALL Build_MLSAuxData (sdId, dataset, prec(:,bank,:), &
            lastIndex=noMAF, dims=dims, fill_value=RAD_ERR_FILL)
    enddo

  END SUBROUTINE OutputL1B_rad

  SUBROUTINE OutputL1B_OA (sdId, noMAF, counterMAF)

    ! This subroutine writes an MAF's worth of data to the L1B OA file

    USE rad_file_contents, ONLY: limb_oa
    USE oa_file_contents, ONLY: emls_oa
 
    ! Arguments:
    INTEGER, INTENT(IN) :: sdId, noMAF, counterMAF

    integer :: status

    TYPE (DataProducts_T) :: dataset

!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

    print *, 'MAF: ', noMAF, limb_oa%ref_lat, limb_oa%ref_time

! 2-d first:

    ALLOCATE (dataset%Dimensions(2), stat=status)

    dataset%data_type = 'real'
    dataset%Dimensions(1) = 'MIF'
    dataset%Dimensions(2) = 'MAF'

    dataset%name = 'sc/OrbIncl'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_orbincl, lastIndex=noMAF)
    dataset%name = 'sc/Lon'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_lon, lastIndex=noMAF)
    dataset%name = 'sc/GeocLat'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_geoclat, lastIndex=noMAF)
    dataset%name = 'sc/GeodLat'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_geodlat, lastIndex=noMAF)
    dataset%name = 'GHz/GeocLat'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%geoc_lat, lastIndex=noMAF)
    dataset%name = 'GHz/GeodLat'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%geod_lat, lastIndex=noMAF)
    dataset%name = 'GHz/GeodAngle'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%geod_ang, lastIndex=noMAF)
    dataset%name = 'GHz/Lon'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%lon, lastIndex=noMAF)
    dataset%name = 'GHz/SolarTime'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%solartime, lastIndex=noMAF)
    dataset%name = 'GHz/SolarZenith'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%solarzenith, lastIndex=noMAF)
    dataset%name = 'GHz/azimAngle'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%azimAngle, lastIndex=noMAF)
    dataset%name = 'GHz/scanAngle'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%scanAngle, lastIndex=noMAF)

    dataset%data_type = 'double'
    dataset%name = 'GHz/GeodAlt'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%geod_alt, lastIndex=noMAF)
    dataset%name = 'GHz/GeocAlt'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%geoc_alt, lastIndex=noMAF)
    dataset%name = 'GHz/LosAngle'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%LosAngle, lastIndex=noMAF)
    dataset%name = 'GHz/LosVel'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%LosVel, lastIndex=noMAF)
    dataset%name = 'GHz/OrbY'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%OrbY, lastIndex=noMAF)
    dataset%name = 'GHz/ECRtoFOV'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%ECRtoFOV, lastIndex=noMAF)
    dataset%name = 'sc/GeocAlt'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_geocalt, lastIndex=noMAF)
    dataset%name = 'sc/GeodAlt'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_geodalt, lastIndex=noMAF)

    dataset%data_type = 'integer'
    dataset%name = 'GHz/BO_stat'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%bo_stat, lastIndex=noMAF)

    ALLOCATE (dataset%Dimensions(1), stat=status)
    dataset%data_type = 'double'
    dataset%Dimensions(1) = 'MAF'
    dataset%name = 'MAFStartTimeTAI'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%MAFStartTimeTAI, &
         lastIndex=noMAF)
    dataset%name = 'sc/MIF_TAI'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_MIF_TAI, lastIndex=noMAF)

    dataset%data_type = 'integer'
    dataset%name      = 'counterMAF'
    CALL Build_MLSAuxData (sdId, dataset, counterMAF, fill_value=INT_FILL, &
         lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

! 3-d next:

    DEALLOCATE (dataset%Dimensions, stat=status)
    ALLOCATE (dataset%Dimensions(3), stat=status)

    dataset%data_type = 'double'
    dataset%Dimensions(1) = 'xyz'
    dataset%Dimensions(2) = 'MIF'
    dataset%Dimensions(3) = 'MAF'
    dataset%name      = 'GHz/ECI'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%ECI, lastIndex=noMAF)
    dataset%name      = 'GHz/ECR'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%ECR, lastIndex=noMAF)
    dataset%name      = 'sc/ECI'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_ECI, lastIndex=noMAF)
    dataset%name      = 'sc/ECR'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_ECR, lastIndex=noMAF)
    dataset%name      = 'sc/VelECI'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_VelECI, lastIndex=noMAF)
    dataset%name      = 'sc/VelECR'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_VelECR, lastIndex=noMAF)
    dataset%name      = 'sc/ypr'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_ypr, lastIndex=noMAF)
    dataset%name      = 'sc/yprRate'
    CALL Build_MLSAuxData (sdId, dataset, emls_oa%sc_ypr_rate, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_OA

END MODULE Output_UARS_L1B
