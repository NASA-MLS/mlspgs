! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!=============================================================================
MODULE OutputL1B
!=============================================================================

  USE MLSCommon
  USE MLSL1Common, ONLY: L1BFileInfo_T
  USE MLSAuxData, ONLY: Build_MLSAuxData, CreateGroup_MLSAuxData
  USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
  USE OutputL1B_DataTypes
  USE HDF5, ONLY: HID_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OutputL1BOA_create, OutputL1B_index, OutputL1B_sc, OutputL1B_GHz, &
       OutputL1B_THz, OutputL1B_rad, OutputL1B_LatBinData, OutputL1B_diags, &
       OutputL1B_DiagsT

  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &
    "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  REAL(r8), PARAMETER :: TAItoGIRD = 1104537627.0d00

CONTAINS

  !----------------------------------- OutputL1B_Create ----
  SUBROUTINE OutputL1BOA_create (sdId, THz)
    ! This subroutine opens/creates the SD output files, and names the arrays
    ! and dimensions contained within them.

    ! Arguments
    TYPE (L1BFileInfo_T), INTENT(IN) :: sdId
    LOGICAL :: THz
  
    IF (THz) RETURN  ! Do not set up L1BOA for THz processing

    CALL CreateGroup_MLSAuxData (sdId%OAId, 'sc')
    CALL CreateGroup_MLSAuxData (sdId%OAId, 'GHz')
    CALL CreateGroup_MLSAuxData (sdId%OAId, 'THz')

  END SUBROUTINE OutputL1BOA_create

  !------------------------------------------------- OuptutL1B_index ----
  SUBROUTINE OutputL1B_index (noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF-SD file
    ! Arguments
    TYPE (L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: sd_id, noMAF
! Variables
    TYPE( DataProducts_T ) :: dataset
    INTEGER :: status
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

    ALLOCATE (dataset%Dimensions(1), stat=status)

    dataset%name      = 'MAFStartTimeUTC   '
    dataset%data_type = 'character         '
    dataset%Dimensions(1) = 'MAF                 '
    CALL Build_MLSAuxData (sd_id,dataset, index%MAFStartTimeUTC, & 
         char_length=lenUTC, lastIndex=noMAF)
    dataset%name      = 'MAFStartTimeTAI   '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id,dataset, index%MAFStartTimeTAI, &
         lastIndex=noMAF)
    dataset%name      = 'noMIFs            '
    dataset%data_type = 'integer           '
    CALL Build_MLSAuxData (sd_id,dataset, index%noMIFs, lastIndex=noMAF)
    dataset%name      = 'counterMAF        '
    CALL Build_MLSAuxData (sd_id,dataset, index%counterMAF, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_index

  !------------------------------------------- OutputL1B_sc ------------
  SUBROUTINE OutputL1B_sc (noMAF, sd_id, sc)
    ! This subroutine writes the spacecraft quantities to the HDF-SD file.

    ! Arguments
    TYPE (L1BOAsc_T), INTENT(IN) :: sc
    INTEGER, INTENT(IN) :: noMAF, sd_id
    ! Variables
    TYPE( DataProducts_T ) :: dataset
    INTEGER :: status
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

! 2-d first:

    ALLOCATE (dataset%Dimensions(2), stat=status)

    dataset%name      = 'sc/GeocAlt        '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'MIF                '
    dataset%Dimensions(2) = 'MAF                '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scGeocAlt, lastIndex=noMAF)
    dataset%name      = 'sc/GeocLat        '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scGeocLat, lastIndex=noMAF)
    dataset%name      = 'sc/GeodAlt        '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scGeodAlt, lastIndex=noMAF)
    dataset%name      = 'sc/GeodLat        '
    dataset%data_type = 'real              '    
    CALL Build_MLSAuxData (sd_id, dataset, sc%scGeodLat, lastIndex=noMAF)
    dataset%name      = 'sc/Lon            '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scLon, lastIndex=noMAF)
    dataset%name      = 'sc/GeodAngle      '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scGeodAngle, lastIndex=noMAF)
    dataset%name      = 'sc/OrbIncl        '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scOrbIncl, lastIndex=noMAF)
    dataset%name      = 'sc/MIF_TAI        '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id, dataset, sc%MIF_TAI, lastIndex=noMAF)

! 3-d next:

    DEALLOCATE (dataset%Dimensions, stat=status)
    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%name      = 'sc/ECI            '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'xyz                 '
    dataset%Dimensions(2) = 'MIF                 '
    dataset%Dimensions(3) = 'MAF                 '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scECI, lastIndex=noMAF)
    dataset%name      = 'sc/ECR            '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scECR, lastIndex=noMAF)
    dataset%name      = 'sc/VelECI         '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scVelECI, lastIndex=noMAF)
    dataset%name      = 'sc/ypr            '
    CALL Build_MLSAuxData (sd_id, dataset, sc%ypr, lastIndex=noMAF)
    dataset%name      = 'sc/yprRate        '
    CALL Build_MLSAuxData (sd_id, dataset, sc%yprRate, lastIndex=noMAF)
    CALL Deallocate_DataProducts (dataset)
    dataset%name      = 'sc/VelECR         '
    CALL Build_MLSAuxData (sd_id, dataset, sc%scVelECR, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_sc

  !-------------------------------------------- OutputL1B_GHz -------
  SUBROUTINE OutputL1B_GHz (noMAF, sd_id, tp)
    ! This subroutine writes the GHz tangent point quantities to the HDF-SD file

    ! Arguments
    TYPE (L1BOAtp_T), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id
    ! Variables
    TYPE( DataProducts_T ) :: dataset
    INTEGER :: status
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

! 2-d first:

    ALLOCATE (dataset%Dimensions(2), stat=status)

    dataset%name      = 'GHz/encoderAngle  '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'MIF                 '
    dataset%Dimensions(2) = 'MAF                 '
    CALL Build_MLSAuxData (sd_id,dataset, tp%encoderAngle, lastIndex=noMAF)
    dataset%name      = 'GHz/scAngle       '
    CALL Build_MLSAuxData (sd_id,dataset, tp%scAngle, lastIndex=noMAF)
    dataset%name      = 'GHz/scanAngle     '
    CALL Build_MLSAuxData (sd_id,dataset, tp%scanAngle, lastIndex=noMAF)
    dataset%name      = 'GHz/scanRate      '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%scanRate, lastIndex=noMAF)
    dataset%name      = 'GHz/OrbY          '
    dataset%data_type = 'real              '
    dataset%Dimensions(1) = 'GHz.MIF              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpOrbY, lastIndex=noMAF)
    dataset%name      = 'GHz/GeocAlt       '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeocAlt, lastIndex=noMAF)
    dataset%name      = 'GHz/GeocLat       '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeocLat, lastIndex=noMAF)
    dataset%name      = 'GHz/GeocAltRate   '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeocAltRate, lastIndex=noMAF)
    dataset%name      = 'GHz/GeodAlt       '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeodAlt, lastIndex=noMAF)
    dataset%name      = 'GHz/GeodLat       '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeodLat, lastIndex=noMAF)
    dataset%name      = 'GHz/GeodAltRate   '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeodAltRate, lastIndex=noMAF)
    dataset%name      = 'GHz/Lon           '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpLon, lastIndex=noMAF)
    dataset%name      = 'GHz/GeodAngle     '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGeodAngle, lastIndex=noMAF)
    dataset%name      = 'GHz/SolarTime     '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpSolarTime, lastIndex=noMAF)
    dataset%name      = 'GHz/SolarZenith   '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpSolarZenith, lastIndex=noMAF)
    dataset%name      = 'GHz/LosAngle      '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpLosAngle, lastIndex=noMAF)
    dataset%name      = 'GHz/LosVel        '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpLosVel, lastIndex=noMAF)

! 3-d next:

    DEALLOCATE (dataset%Dimensions, stat=status)
    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%name      = 'GHz/ECI           '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'xyz                 '
    dataset%Dimensions(2) = 'GHz.MIF             '
    dataset%Dimensions(3) = 'MAF                 '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpECI, lastIndex=noMAF)
    dataset%name      = 'GHz/ECR           '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpECR, lastIndex=noMAF)
    dataset%name      = 'GHz/ECRtoFOV      '
    dataset%Dimensions(1) = '3x3                 '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpECRtoFOV, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_GHz
  !-------------------------------------------- OutputL1B_THz --------------
  SUBROUTINE OutputL1B_THz (noMAF, sd_id, tp)
    ! This subroutine writes the THz tangent point quantities to the HDF-SD file

    ! Arguments
    TYPE (L1BOAtp_T), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id
    ! Variables
    TYPE( DataProducts_T ) :: dataset
    INTEGER :: status
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

! 2-d first:

    ALLOCATE (dataset%Dimensions(2), stat=status)
    dataset%name      = 'THz/encoderAngle  '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'MIF'
    dataset%Dimensions(2) = 'MAF' 
    CALL Build_MLSAuxData (sd_id, dataset, tp%encoderAngle, lastIndex=noMAF)
    dataset%name      = 'THz/scAngle       '
    CALL Build_MLSAuxData (sd_id, dataset, tp%scAngle, lastIndex=noMAF)
    dataset%name      = 'THz/scanAngle     '
    CALL Build_MLSAuxData (sd_id, dataset, tp%scanAngle, lastIndex=noMAF)
    dataset%name      = 'THz/scanRate      '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id, dataset, tp%scanRate, lastIndex=noMAF)
    dataset%name      = 'THz/OrbY          '
    dataset%Dimensions(1) = 'THz.MIF'
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpOrbY, lastIndex=noMAF)
    dataset%name      = 'THz/GeocAlt       '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeocAlt, lastIndex=noMAF)
    dataset%name      = 'THz/GeocLat       '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeocLat, lastIndex=noMAF)
    dataset%name      = 'THz/GeocAltRate   '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeocAltRate, lastIndex=noMAF)
    dataset%name      = 'THz/GeodAlt       '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeodAlt, lastIndex=noMAF)
    dataset%name      = 'THz/GeodLat       '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeodLat, lastIndex=noMAF)
    dataset%name      = 'THz/GeodAltRate   '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeodAltRate, lastIndex=noMAF)
    dataset%name      = 'THz/Lon           '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpLon, lastIndex=noMAF)
    dataset%name      = 'THz/GeodAngle     '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeodAngle, lastIndex=noMAF)
    dataset%name      = 'THz/SolarTime     '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpSolarTime, lastIndex=noMAF)
    dataset%name      = 'THz/SolarZenith   '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpSolarZenith, lastIndex=noMAF)
    dataset%name      = 'THz/LosAngle      '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpLosAngle, lastIndex=noMAF)
    dataset%name      = 'THz/LosVel        '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpLosVel, lastIndex=noMAF)

! 3-d next:

    DEALLOCATE (dataset%Dimensions, stat=status)
    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%name      = 'THz/ECI           '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'xyz'
    dataset%Dimensions(2) = 'THz.MIF'
    dataset%Dimensions(3) = 'MAF' 
    CALL Build_MLSAuxData (sd_id, dataset,tp%tpECI, lastIndex=noMAF)
    dataset%name      = 'THz/ECR           '
    CALL Build_MLSAuxData (sd_id, dataset,tp%tpECR, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_THz

  !-------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad (noMAF, sdId, counterMAF, Reflec, MAFStartTimeTAI, &
       rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files

    USE EngTbls, ONLY: Reflec_T
    USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz
    USE MLSL1Rad, ONLY: Radiance_T, Rad_Name
    USE MLSSignalNomenclature, ONLY:GetFullMLSSignalName
  
    ! Arguments
    TYPE (L1BFileInfo_T) :: sdId
    TYPE (Radiance_T) :: rad(:)
    TYPE (Reflec_T) :: Reflec
    INTEGER, INTENT(IN) :: counterMAF, noMAF
    REAL(r8), INTENT (IN) :: MAFStartTimeTAI
    REAL(r8) :: MAFStartTimeGIRD

   ! Variables
    TYPE( DataProducts_T ) :: dataset    
    CHARACTER (LEN=64) :: dim_chan, dim_mif, name, prec
    INTEGER, DIMENSION(3) :: dims
    INTEGER(HID_T) :: sd_id
    INTEGER :: i, status

! For use in determining filling missing radiance records

    INTEGER, PARAMETER :: max_rad_nos = SIZE (Rad_Name)
    INTEGER, SAVE :: rad_sdid(max_rad_nos) = -999
    INTEGER, SAVE :: rad_dim1(max_rad_nos) = 0
    INTEGER, SAVE :: rad_dim2(max_rad_nos) = 0
    LOGICAL, SAVE :: rad_out(max_rad_nos) = .FALSE. ! whether rad has output
    LOGICAL, SAVE :: fill_last(max_rad_nos)   ! whether to fill last record
    CHARACTER(len=10), SAVE :: rad_chan(max_rad_nos) = " "
    CHARACTER(len=10), SAVE :: rad_mif(max_rad_nos) = " "

    REAL, PARAMETER :: RAD_FILL = -999.9, RAD_ERR_FILL = -1.0
    REAL, PARAMETER :: RAD_FILL_ARR(129,148) = RAD_FILL
    REAL, PARAMETER :: RAD_ERR_FILL_ARR(129,148) = RAD_ERR_FILL

    ! Convert TAI time to GIRD time

    MAFStartTimeGIRD = MAFStartTimeTAI + TAItoGIRD

! Find the sds_id number for counterMAF in the file RADG, write the data to it,
! terminate access to the data set
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)

    ALLOCATE (dataset%Dimensions(1), stat=status)
    dataset%name      = 'counterMAF        '
    dataset%data_type = 'integer           '
    dataset%Dimensions(1) = 'MAF'
    IF (sdId%RADDID /= 0) THEN
       CALL Build_MLSAuxData (sdId%RADDID, dataset, counterMAF, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, counterMAF, lastIndex=noMAF)
    ELSE
       CALL Build_MLSAuxData (sdId%RADTID, dataset, counterMAF, lastIndex=noMAF)
    ENDIF

!------------------------------------------------------------------------------

    dataset%name      = 'MAFStartTimeGIRD  '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'MAF'
    IF (sdId%RADDID /= 0) THEN
       CALL Build_MLSAuxData (sdId%RADDID, dataset, MAFStartTimeGIRD, &
            lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, MAFStartTimeGIRD, &
            lastIndex=noMAF)
    ELSE
       CALL Build_MLSAuxData (sdId%RADTID, dataset, MAFStartTimeGIRD, &
            lastIndex=noMAF)
    ENDIF

! Reflector temperatures:

    dataset%data_type = 'real              '
    dataset%name      = 'Pri_Reflec        '
    IF (sdId%RADDID /= 0) THEN
       CALL Build_MLSAuxData (sdId%RADDID, dataset, Reflec%Pri, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Pri, lastIndex=noMAF)
    ENDIF
    dataset%name      = 'Sec_Reflec        '
    IF (sdId%RADDID /= 0) THEN
       CALL Build_MLSAuxData (sdId%RADDID, dataset, Reflec%Sec, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Sec, lastIndex=noMAF)
    ENDIF
    dataset%name      = 'Ter_Reflec        '
    IF (sdId%RADDID /= 0) THEN
       CALL Build_MLSAuxData (sdId%RADDID, dataset, Reflec%Ter, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Ter, lastIndex=noMAF)
    ENDIF

    DEALLOCATE (dataset%Dimensions, stat=status)

! Prepare for radiance output

    fill_last = .TRUE.   ! will fill last record unless good data

    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%data_type = 'real'
    dataset%Dimensions(3) = 'MAF                 '

!------------------------------------------------------------------------------
    DO i = 1, SIZE (rad) ! Loop on number of SDs per MAF

       CALL GetFullMLSSignalName(rad(i)%signal, name) ! Concatenate SD names
       prec = TRIM(name) // ' precision'
       ! Set parameters based on input data dimensions
       ! Based on the SD name, set dim name for channel, get Id of output file
       IF (INDEX(name, 'FB') /= 0 ) THEN
          dim_chan = 'chanFB              '
          IF (sdId%RADGID /= 0) THEN
             sd_id = sdId%RADGID
          ELSE
             sd_id = sdId%RADTID
          ENDIF
       ELSE IF (INDEX(name, 'MB') /= 0 ) THEN
          dim_chan = 'chanMB              '
          sd_id = sdId%RADGID
       ELSE IF (INDEX(name, 'WF') /= 0 ) THEN
          dim_chan = 'chanWF              '
          sd_id = sdId%RADGID
       ELSE IF (INDEX(name, 'DACS') /= 0 ) THEN
          dim_chan = 'chanDACS            '
          sd_id = sdId%RADDID
       ELSE
          dim_chan = ''
          sd_id = -999
       ENDIF

       IF ((ANY(rad(i)%value /= 0.0)) .AND. (sd_id /= -999)) THEN   
          ! if good radiance data exists 
          ! Create/Open the value SDs
          dims(1) = SIZE (rad(i)%value,1)
          dims(2) = SIZE (rad(i)%value,2)
          dims(3) = 1
          IF (INDEX(name,'R5') /= 0 ) THEN
             dim_mif = 'THz.MIF              '
             dims(2) = MIFsTHz      !! lenT
          ELSE
             dim_mif = 'GHz.MIF              '
             dims(2) = MIFsGHz      !! lenG
          ENDIF
!------------------------------------------------Create/Open the value datasets
          dataset%Dimensions(1) = dim_chan
          dataset%Dimensions(2) = dim_mif
          dataset%name = TRIM(name)
          CALL Build_MLSAuxData (sd_id, dataset, rad(i)%value, &
               lastIndex=noMAF, dims=dims, fill_value=RAD_FILL)

!------------------------------------------ Create/Open the precision datasets

          dataset%name = TRIM(prec)
          CALL Build_MLSAuxData (sd_id, dataset, rad(i)%precision, &
               lastIndex=noMAF, dims=dims, fill_value=RAD_ERR_FILL)
!-----------------------------------------------------------------------------

! Save matching table entries

          WHERE (name == Rad_Name)
             rad_out = .TRUE.
             fill_last = .FALSE.
             rad_dim1 = dims(1)
             rad_dim2 = dims(2)
             rad_chan = dim_chan
             rad_mif = dim_mif
             rad_sdid = sd_id
          ENDWHERE

       ENDIF
    ENDDO

! Determine filling the radiance file record:

    DO i = 1, max_rad_nos
       IF (rad_out(i) .AND. fill_last(i)) THEN
          dataset%Dimensions(1) = rad_chan(i)
          dataset%Dimensions(2) = rad_mif(i)
          dataset%name = TRIM(Rad_Name(i))
          dims(1) = rad_dim1(i)
          dims(2) = rad_dim2(i)
          dims(3) = 1
          CALL Build_MLSAuxData (rad_sdid(i), dataset, RAD_FILL_ARR, &
               lastIndex=noMAF, dims=dims, fill_value=RAD_FILL)
          dataset%name = TRIM(Rad_Name(i))// ' precision'
          CALL Build_MLSAuxData (rad_sdid(i), dataset, RAD_ERR_FILL_ARR, &
               lastIndex=noMAF, dims=dims, fill_value=RAD_ERR_FILL)
       ENDIF
    ENDDO

    CALL Deallocate_DataProducts (dataset)
 
  END SUBROUTINE OutputL1B_rad

!=============================================================================
  SUBROUTINE OutputL1B_LatBinData (noMAF, sd_id, AscDescIndx, LatBinIndx, &
       LatBinChanAvg, BaselineAlt, LatBin)
!=============================================================================

    USE MLSHDF5, ONLY: SaveAsHDF5DS, MakeHDF5Attribute

    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
    INTEGER, DIMENSION(:), OPTIONAL :: AscDescIndx
    INTEGER, DIMENSION(:), OPTIONAL :: LatBinIndx
    REAL, DIMENSION(:,:,:,:), OPTIONAL :: LatBinChanAvg
    REAL, DIMENSION(:,:,:), OPTIONAL :: LatBin
    REAL, DIMENSION(:,:), OPTIONAL :: BaselineAlt

    TYPE( DataProducts_T ) :: dataset
    CHARACTER(LEN=16) :: DimName(4)
    INTEGER :: dims(3), status

    CALL Deallocate_DataProducts (dataset)

    IF (PRESENT (AscDescIndx) .AND. PRESENT (LatBinIndx)) THEN

       ALLOCATE (dataset%Dimensions(2), stat=status)
       dataset%name = 'AscDescIndx   '
       dataset%data_type = 'integer           '
       dataset%Dimensions(1) = 'GHz.MIF             '
       dims = 1
       dims(1) = SIZE (AscDescIndx)
       CALL Build_MLSAuxData (sd_id, dataset, AscDescIndx, lastIndex=noMAF, &
            dims=dims)
       dataset%name = 'LatBinIndx    '
       CALL Build_MLSAuxData (sd_id, dataset, LatBinIndx, lastIndex=noMAF, &
            dims=dims)

    ENDIF

    IF (PRESENT (LatBinChanAvg)) THEN
       CALL SaveAsHDF5DS (sd_id, 'LatBinChanAvg', LatBinChanAvg)
       DimName(1) = 'FBChan'
       DimName(2) = 'GHzBand'
       DimName(3) = 'LatBin'
       DimName(4) = 'AscDesc'
       CALL MakeHDF5Attribute (sd_id, 'LatBinChanAvg', 'Dimensions', DimName)
    ENDIF

    IF (PRESENT (LatBin)) THEN
       CALL SaveAsHDF5DS (sd_id, 'BaselineLatBin', LatBin)
       DimName(1) = 'Min/Max lat'
       DimName(2) = 'LatBin'
       DimName(3) = 'AscDesc'
       CALL MakeHDF5Attribute (sd_id, 'BaselineLatBin', 'Dimensions', &
            DimName(1:3))
    ENDIF

    IF (PRESENT (BaselineAlt)) THEN
       CALL SaveAsHDF5DS (sd_id, 'BaselineAlt', BaselineAlt)
       DimName(1) = 'FBchan'
       DimName(2) = 'GHzBand'
       CALL MakeHDF5Attribute (sd_id, 'BaselineAlt', 'Dimensions', DimName(1:2))
    ENDIF

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_LatBinData

!=============================================================================
  SUBROUTINE OutputL1B_diags (noMAF, sd_Id, counterMAF, MAFStartTimeTAI)
!=============================================================================

    ! This subroutine writes an MAF's worth of diagnostic data

    USE MLSL1Common, ONLY: FBchans, GHzNum, MBchans, MBnum, WFchans, WFnum
    USE Calibration, ONLY: Chi2, Tsys, Cgain

    INTEGER, INTENT(IN) :: counterMAF, noMAF, sd_id
    REAL(r8), INTENT (IN) :: MAFStartTimeTAI

    INTEGER :: dims(3), status
    TYPE (DataProducts_T) :: dataset
    REAL(r8) :: MAFStartTimeGIRD

    ! Convert TAI time to GIRD time

    MAFStartTimeGIRD = MAFStartTimeTAI + TAItoGIRD

    ALLOCATE (dataset%Dimensions(1), stat=status)
    dataset%name      = 'counterMAF        '
    dataset%data_type = 'integer           '
    dataset%Dimensions(1) = 'MAF'
    CALL Build_MLSAuxData (sd_Id, dataset, counterMAF, lastIndex=noMAF)

    dataset%name      = 'MAFStartTimeGIRD  '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_Id, dataset, MAFStartTimeGIRD, lastIndex=noMAF)
    CALL Deallocate_DataProducts (dataset)

    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%data_type = 'real'
    dataset%Dimensions(3) = 'MAF                 '

    dims(1) = FBchans
    dims(2) = GHzNum
    dims(3) = 1
    dataset%Dimensions(1) = 'FBchan              '
    dataset%Dimensions(2) = 'GHzNum              '
    dataset%name = 'Tsys FB'

    CALL Build_MLSAuxData (sd_id, dataset, Tsys%FB, lastIndex=noMAF, dims=dims)
    dataset%name = 'Chi2 FB'
    CALL Build_MLSAuxData (sd_id, dataset, Chi2%FB, lastIndex=noMAF, dims=dims)
    dataset%name = 'ChanGain FB'
    CALL Build_MLSAuxData (sd_id, dataset, Cgain%FB, lastIndex=noMAF, dims=dims)

    dims(1) = MBchans
    dims(2) = MBnum
    dataset%Dimensions(1) = 'MBchan              '
    dataset%Dimensions(2) = 'MBnum               '
    dataset%name = 'Tsys MB'

    CALL Build_MLSAuxData (sd_id, dataset, Tsys%MB, lastIndex=noMAF, dims=dims)
    dataset%name = 'Chi2 MB'
    CALL Build_MLSAuxData (sd_id, dataset, Chi2%MB, lastIndex=noMAF, dims=dims)
    dataset%name = 'ChanGain MB'
    CALL Build_MLSAuxData (sd_id, dataset, Cgain%MB, lastIndex=noMAF, dims=dims)

    dims(1) = WFchans
    dims(2) = WFnum
    dataset%Dimensions(1) = 'WFchan              '
    dataset%Dimensions(2) = 'WFnum               '
    dataset%name = 'Tsys WF'

    CALL Build_MLSAuxData (sd_id, dataset, Tsys%WF, lastIndex=noMAF, dims=dims)
    dataset%name = 'Chi2 WF'
    CALL Build_MLSAuxData (sd_id, dataset, Chi2%WF, lastIndex=noMAF, dims=dims)
    dataset%name = 'ChanGain WF'
    CALL Build_MLSAuxData (sd_id, dataset, Cgain%WF, lastIndex=noMAF, dims=dims)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_diags

!=============================================================================
  SUBROUTINE OutputL1B_diagsT (sd_Id, MAFno, counterMAF, MAFStartTimeTAI, &
       nvBounds, OrbNo, Chisq, dLlo, yTsys)
!=============================================================================

    USE MLSL1Common, ONLY: THzchans, THzNum

    INTEGER, INTENT(IN) :: sd_id
    INTEGER, INTENT(IN), OPTIONAL :: counterMAF, MAFno, nvBounds, OrbNo
    REAL(r8), INTENT (IN), OPTIONAL :: MAFStartTimeTAI
    REAL(r8), DIMENSION(:,:), INTENT (IN), OPTIONAL :: Chisq, dLlo, yTsys

    INTEGER :: dims(3), status
    TYPE (DataProducts_T) :: dataset
    REAL(r8) :: MAFStartTimeGIRD

    dims(1) = THzchans
    dims(2) = THzNum
    dims(3) = 1

    IF (PRESENT (counterMAF) .AND. PRESENT (MAFno)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(1), stat=status)
       dataset%name      = 'counterMAF        '
       dataset%data_type = 'integer           '
       dataset%Dimensions(1) = 'MAF'
       CALL Build_MLSAuxData (sd_Id, dataset, counterMAF, lastIndex=MAFno)

    ENDIF

    IF (PRESENT (MAFStartTimeTAI) .AND. PRESENT (MAFno)) THEN

       ! Convert TAI time to GIRD time

       MAFStartTimeGIRD = MAFStartTimeTAI + TAItoGIRD

       dataset%name      = 'MAFStartTimeGIRD  '
       dataset%data_type = 'double            '
       CALL Build_MLSAuxData (sd_Id, dataset, MAFStartTimeGIRD, lastIndex=MAFno)
       CALL Deallocate_DataProducts (dataset)

    ENDIF

    IF (PRESENT (nvBounds) .AND. PRESENT (MAFno)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(1), stat=status)
       dataset%name      = 'nvBounds          '
       dataset%data_type = 'integer           '
       dataset%Dimensions(1) = 'MAF'
       CALL Build_MLSAuxData (sd_Id, dataset, nvBounds, lastIndex=MAFno)

    ENDIF

    IF (PRESENT (Chisq) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(2), stat=status)
       dataset%name      = 'Chisq             '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       CALL Build_MLSAuxData (sd_Id, dataset, Chisq, lastIndex=OrbNo, dims=dims)

    ENDIF

    IF (PRESENT (dLlo) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(2), stat=status)
       dataset%name      = 'dLlo              '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       CALL Build_MLSAuxData (sd_Id, dataset, dLlo, lastIndex=OrbNo, dims=dims)

    ENDIF

 
    IF (PRESENT (yTsys) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(2), stat=status)
       dataset%name      = 'yTsys             '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       CALL Build_MLSAuxData (sd_Id, dataset, yTsys, lastIndex=OrbNo, dims=dims)

    ENDIF

  END SUBROUTINE OutputL1B_diagsT

!=============================================================================
END MODULE OutputL1B
!=============================================================================

! $Log$
! Revision 2.12  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.11  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.10  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.9  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.8  2002/11/20 15:44:02  perun
! Use HDFversion instead of HDFVersionString & remove "default" calls
!
! Revision 2.7  2002/11/07 21:58:30  jdone
! Incorporated HDF4/HDF5 switch; moved HDF4 routines to OutputL1B_HDF4.f90
!
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
