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
       OutputL1B_DiagsT, OutputL1B_Chi2

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  ! TAItoGIRD is the number of TAI seconds between Jan 1 1993 and Jan
  ! 1 1958 (both at all-balls UTC). There are 12784 days between those
  ! 2 dates, so that's 1104537600 *civil* seconds (assuming 86400.0
  ! seconds/day). There are 27 leap-seconds between Jan 1, 1961 and
  ! Jan 1, 1993, so GIRD is the number of TAI-seconds, not civil time
  ! seconds, i.e. it includes leap setons
  
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
    dataset%name      = 'GHz/azimAngle     '
    CALL Build_MLSAuxData (sd_id,dataset, tp%azimAngle, lastIndex=noMAF)
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
    dataset%name      = 'GHz/BO_stat       '
    dataset%data_type = 'integer           '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpBO_stat, lastIndex=noMAF)
    dataset%name      = 'GHz/GalLat        '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGalLat, lastIndex=noMAF)
    dataset%name      = 'GHz/GalLon        '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGalLon, lastIndex=noMAF)

! 3-d next:

    DEALLOCATE (dataset%Dimensions, stat=status)
    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%name      = 'GHz/ECI           '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'xyz                 '
    dataset%Dimensions(2) = 'GHz.MIF             '
    dataset%Dimensions(3) = 'MAF                 '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpECI, lastIndex=noMAF)
    dataset%name      = 'GHz/ECR           '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpECR, lastIndex=noMAF)
    dataset%name      = 'GHz/ECRtoFOV      '
    dataset%Dimensions(1) = '3x3                 '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpECRtoFOV, lastIndex=noMAF)

    dataset%name      = 'GHz/Pos_Prime     '
    dataset%data_type = 'real              '
    dataset%Dimensions(1) = '2                   '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpPos_P, lastIndex=noMAF)

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
    dataset%name      = 'THz/azimAngle     '
    CALL Build_MLSAuxData (sd_id, dataset, tp%azimAngle, lastIndex=noMAF)
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
    dataset%name      = 'THz/GeodAltX      '
    dataset%data_type = 'double            '
    CALL Build_MLSAuxData (sd_id, dataset, tp%tpGeodAltX, lastIndex=noMAF)
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
    dataset%name      = 'THz/BO_stat       '
    dataset%data_type = 'integer           '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpBO_stat, lastIndex=noMAF)
    dataset%name      = 'THz/GalLat        '
    dataset%data_type = 'real              '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGalLat, lastIndex=noMAF)
    dataset%name      = 'THz/GalLon        '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpGalLon, lastIndex=noMAF)

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
    dataset%name      = 'THz/ECRtoFOV      '
    dataset%Dimensions(1) = '3x3                 '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpECRtoFOV, lastIndex=noMAF)
    dataset%name      = 'THz/Pos_Prime     '
    dataset%data_type = 'real              '
    dataset%Dimensions(1) = '2                   '
    CALL Build_MLSAuxData (sd_id,dataset, tp%tpPos_P, lastIndex=noMAF)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_THz

  !-------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad (noMAF, sdId, counterMAF, Reflec, MAFStartTimeTAI, &
       rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files

    USE MLSCommon, ONLY: DEFAULTUNDEFINEDVALUE
    USE EngTbls, ONLY: Reflec_T
    USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz, L1Config
    USE MLSL1Rad, ONLY: Radiance_T, Rad_Name
    USE MLSSignalNomenclature, ONLY:GetFullMLSSignalName
    USE SpectralBaseline, ONLY: Baseline, BaselineAC, BaselineDC
  
    ! Arguments
    TYPE (L1BFileInfo_T), INTENT(IN) :: sdId
    TYPE (Radiance_T), INTENT(IN) :: rad(:)
    TYPE (Reflec_T), INTENT(IN) :: Reflec
    INTEGER, INTENT(IN) :: counterMAF, noMAF
    REAL(r8), INTENT (IN) :: MAFStartTimeTAI

   ! Variables
    TYPE( DataProducts_T ) :: dataset, baselineDS  
    CHARACTER (LEN=64) :: dim_chan, dim_mif, name, prec
    INTEGER :: dims(3)
    INTEGER(HID_T) :: sd_id
    INTEGER :: i, status, bandno
    REAL(r8) :: MAFStartTimeGIRD

! For use in determining filling missing radiance records

    INTEGER, PARAMETER :: max_rad_nos = SIZE (Rad_Name)
    INTEGER, SAVE :: rad_sdid(max_rad_nos) = -999
    INTEGER, SAVE :: rad_dim1(max_rad_nos) = 0
    INTEGER, SAVE :: rad_dim2(max_rad_nos) = 0
    LOGICAL, SAVE :: rad_out(max_rad_nos) = .FALSE. ! whether rad has output
    LOGICAL, SAVE :: fill_last(max_rad_nos)   ! whether to fill last record
    CHARACTER(len=10), SAVE :: rad_chan(max_rad_nos) = " "
    CHARACTER(len=10), SAVE :: rad_mif(max_rad_nos) = " "

    REAL, PARAMETER :: RAD_FILL = DEFAULTUNDEFINEDVALUE, RAD_ERR_FILL = -1.0
    REAL, PARAMETER :: RAD_FILL_ARR(129,148) = RAD_FILL
    REAL, PARAMETER :: RAD_ERR_FILL_ARR(129,148) = RAD_ERR_FILL
    INTEGER, PARAMETER :: INT_FILL = -1

    ! Convert TAI time to GIRD time
    ! GIRD time is TAI seconds since Jan 1, 1958 00:00:00 UTC.
    MAFStartTimeGIRD = MAFStartTimeTAI + TAItoGIRD

! Find the sds_id number for counterMAF in the file RADG, write the data to it,
! terminate access to the data set
!------------------------------------------------------------------------------
    CALL Deallocate_DataProducts (dataset)
    CALL Deallocate_DataProducts (baselineDS)

    ALLOCATE (dataset%Dimensions(1), stat=status)
    dataset%name      = 'counterMAF        '
    dataset%data_type = 'integer           '
    dataset%Dimensions(1) = 'MAF'
    IF (sdId%RADGID /= 0) THEN

      !<whd> 
      ! I assume that RADGID is set ONLY when we're in MLSL1G, so
      ! this is basically a way of testing on whether the caller is
      ! MLSL1G 
      ! </whd>

       IF (sdId%RADDID /= 0) CALL Build_MLSAuxData (sdId%RADDID, dataset, &
            counterMAF, fill_value=INT_FILL, lastIndex=noMAF) !<whd> doing DACS processing </whd>

       CALL Build_MLSAuxData (sdId%RADGID, dataset, counterMAF, &
            fill_value=INT_FILL, lastIndex=noMAF)
    ELSE
      !<whd> Or here we're coming from  MLSL1T, the Teraherz module </whd>
       CALL Build_MLSAuxData (sdId%RADTID, dataset, counterMAF, &
            fill_value=INT_FILL, lastIndex=noMAF)
    ENDIF

!------------------------------------------------------------------------------

    dataset%name      = 'MAFStartTimeGIRD  '
    dataset%data_type = 'double            '
    dataset%Dimensions(1) = 'MAF'
    IF (sdId%RADGID /= 0) THEN
       IF (sdId%RADDID /= 0) CALL Build_MLSAuxData (sdId%RADDID, dataset, &
            MAFStartTimeGIRD, fill_value=-1.0d0, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, MAFStartTimeGIRD, &
            fill_value=-1.0d0, lastIndex=noMAF)
    ELSE
       CALL Build_MLSAuxData (sdId%RADTID, dataset, MAFStartTimeGIRD, &
            fill_value=-1.0d0, lastIndex=noMAF)
    ENDIF

! Reflector temperatures:

    dataset%data_type = 'real              '
    dataset%name      = 'Pri_Reflec        '
    IF (sdId%RADGID /= 0) THEN
       IF (sdId%RADDID /= 0) CALL Build_MLSAuxData (sdId%RADDID, dataset, &
            Reflec%Pri, fill_value=RAD_ERR_FILL, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Pri, &
            fill_value=RAD_ERR_FILL, lastIndex=noMAF)
    ENDIF
    dataset%name      = 'Sec_Reflec        '
    IF (sdId%RADGID /= 0) THEN
       IF (sdId%RADDID /= 0) CALL Build_MLSAuxData (sdId%RADDID, dataset, &
            Reflec%Sec, fill_value=RAD_ERR_FILL, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Sec, &
            fill_value=RAD_ERR_FILL, lastIndex=noMAF)
    ENDIF
    dataset%name      = 'Ter_Reflec        '
    IF (sdId%RADGID /= 0) THEN
       IF (sdId%RADDID /= 0) CALL Build_MLSAuxData (sdId%RADDID, dataset, &
            Reflec%Ter, fill_value=RAD_ERR_FILL, lastIndex=noMAF)
       CALL Build_MLSAuxData (sdId%RADGID, dataset, Reflec%Ter, &
            fill_value=RAD_ERR_FILL, lastIndex=noMAF)
    ENDIF

    DEALLOCATE (dataset%Dimensions, stat=status)

! Prepare for radiance output

    fill_last = .TRUE.   ! will fill last record unless good data

    ALLOCATE (dataset%Dimensions(3), stat=status)
    dataset%data_type = 'real'
    dataset%Dimensions(3) = 'MAF                 '

    ALLOCATE (baselineDS%Dimensions(2), stat=status)
    baselineDS%data_type = 'real'
    baselineDS%Dimensions(2) = 'MAF                 '

!------------------------------------------------------------------------------
    DO i = 1, SIZE (rad) ! Loop on number of SDs per MAF

       bandno = rad(i)%bandno

       IF (L1Config%Output%DisableRadOut(bandno)) CYCLE   ! skip to next bandno

       CALL GetFullMLSSignalName(rad(i)%signal, name) ! Concatenate SD names
       prec = TRIM(name) // ' precision'
       ! Set parameters based on input data dimensions
       ! Based on the SD name, set dim name for channel, get Id of output file
       IF (INDEX(name, 'FB') /= 0 ) THEN
          dim_chan = 'chanFB              '
          IF (sdId%RADGID /= 0) THEN
             IF (INDEX(name,'R5') /= 0 ) THEN  ! Don't allow R5 in GHz file: <whd>R5 is the THz radiometer</whd>
                dim_chan = ''
                sd_id = -999
             ELSE
                sd_id = sdId%RADGID
             ENDIF
          ELSE
             sd_id = sdId%RADTID
          ENDIF
       ELSE IF (INDEX(name, 'MB') /= 0) THEN
          dim_chan = 'chanMB              '
          sd_id = sdId%RADGID
       ELSE IF (INDEX(name, 'WF') /= 0) THEN
          dim_chan = 'chanWF              '
          sd_id = sdId%RADGID
       ELSE IF (INDEX(name, 'DACS') /= 0 .AND. L1Config%Calib%CalibDACS) THEN
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

! Output diagnostic offsets (GHz only!):

          IF (sd_id /= sdId%RADTID) THEN
             dataset%name = TRIM(name) // ' Poffset'
             CALL Build_MLSAuxData (sdId%DiagId, dataset, rad(i)%Poffset(:,1), &
                  lastIndex=noMAF, fill_value=RAD_FILL)  ! Only 1 per MAF
             dataset%name = TRIM(name) // ' ModelOffset'
             CALL Build_MLSAuxData (sdId%DiagId, dataset, rad(i)%ModelOffset, &
                  lastIndex=noMAF, fill_value=RAD_FILL)

          ENDIF

! Output baselines

          baselineDS%Dimensions(1) = dim_chan
          baselineDS%name = TRIM(name) // ' Baseline'
          CALL Build_MLSAuxData (sd_id, baselineDS, Baseline(bandno)%offset, &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(name) // ' Baseline precision'
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               Baseline(bandno)%precision, lastIndex=noMAF, &
               fill_value=RAD_ERR_FILL)
          baselineDS%name = TRIM(name) // ' BaselineAC'
          CALL Build_MLSAuxData (sd_id, baselineDS, BaselineAC(bandno)%offset, &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(name) // ' BaselineAC precision'
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               BaselineAC(bandno)%precision, lastIndex=noMAF, &
               fill_value=RAD_ERR_FILL)
          baselineDS%name = TRIM(name) // ' BaselineDC'
          CALL Build_MLSAuxData (sd_id, baselineDS, BaselineDC(bandno)%offset, &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(name) // ' BaselineDC precision'
          CALL Build_MLSAuxData (sd_id, baselineDS, &
               BaselineDC(bandno)%precision, lastIndex=noMAF, &
               fill_value=RAD_ERR_FILL)

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

! Fill Baselines:

          baselineDS%Dimensions(1) = rad_chan(i)
          baselineDS%name = TRIM(Rad_Name(i)) // ' Baseline'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(Rad_Name(i)) // ' Baseline precision'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_ERR_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(Rad_Name(i)) // ' BaselineAC'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(Rad_Name(i)) // ' BaselineAC precision'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_ERR_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(Rad_Name(i)) // ' BaselineDC'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
          baselineDS%name = TRIM(Rad_Name(i)) // ' BaselineDC precision'
          CALL Build_MLSAuxData (rad_sdid(i), baselineDS, &
               RAD_ERR_FILL_ARR(1:rad_dim1(i),1), &
               lastIndex=noMAF, fill_value=RAD_FILL)
       ENDIF
    ENDDO

    CALL Deallocate_DataProducts (dataset)
    CALL Deallocate_DataProducts (baselineDS)
 
  END SUBROUTINE OutputL1B_rad

!=============================================================================
  SUBROUTINE OutputL1B_LatBinData (noMAF, sd_id, AscDescIndx, LatBinIndx, &
       LatBinChanAvg, BaselineAlt, LatBin, BaselineBandNo, BinnedBaseline, &
       Name)
!=============================================================================

    USE MLSHDF5, ONLY: SaveAsHDF5DS, MakeHDF5Attribute

    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
    INTEGER, DIMENSION(:), OPTIONAL :: AscDescIndx, LatBinIndx, BaselineBandNo
    REAL, DIMENSION(:,:,:), OPTIONAL :: LatBinChanAvg
    REAL, DIMENSION(:,:), OPTIONAL :: LatBin, BinnedBaseline
    REAL, DIMENSION(:,:), OPTIONAL :: BaselineAlt
    CHARACTER(len=*), OPTIONAL :: Name

    TYPE( DataProducts_T ) :: dataset
    CHARACTER(LEN=16) :: DimName(3)
    INTEGER :: dims(3), status

    CALL Deallocate_DataProducts (dataset)

    IF (PRESENT (AscDescIndx) .AND. PRESENT (LatBinIndx)) THEN

       ALLOCATE (dataset%Dimensions(2), stat=status)
       dataset%name = 'AscDescIndx   '
       dataset%data_type = 'integer           '
       dataset%Dimensions(1) = 'GHz.MIF             '
       dataset%Dimensions(2) = 'MAF                 '
       dims = 1
       dims(1) = SIZE (AscDescIndx)
       CALL Build_MLSAuxData (sd_id, dataset, AscDescIndx, lastIndex=noMAF, &
            dims=dims)
       dataset%name = 'LatBinIndx    '
       CALL Build_MLSAuxData (sd_id, dataset, LatBinIndx, lastIndex=noMAF, &
            dims=dims)

    ENDIF

    IF (PRESENT (BinnedBaseline) .AND. PRESENT (Name)) THEN
       CALL SaveAsHDF5DS (sd_id, Name, BinnedBaseline)
       DimName(1) = 'Chan'
       DimName(2) = 'LatBin'
       CALL MakeHDF5Attribute (sd_id, Name, 'Dimensions', DimName(1:2))
    ENDIF

    IF (PRESENT (LatBinChanAvg)) THEN
       CALL SaveAsHDF5DS (sd_id, 'LatBinChanAvg', LatBinChanAvg)
       DimName(1) = 'FBChan'
       DimName(2) = 'BandNo'
       DimName(3) = 'LatBin'
       CALL MakeHDF5Attribute (sd_id, 'LatBinChanAvg', 'Dimensions', DimName)
    ENDIF

    IF (PRESENT (LatBin)) THEN
       CALL SaveAsHDF5DS (sd_id, 'BaselineLatBin', LatBin)
       DimName(1) = 'Min/Max lat'
       DimName(2) = 'LatBin'
       CALL MakeHDF5Attribute (sd_id, 'BaselineLatBin', 'Dimensions', &
            DimName(1:2))
    ENDIF

    IF (PRESENT (BaselineAlt)) THEN
       CALL SaveAsHDF5DS (sd_id, 'BaselineAlt', BaselineAlt)
       DimName(1) = 'FBchan'
       DimName(2) = 'BandNo'
       CALL MakeHDF5Attribute (sd_id, 'BaselineAlt', 'Dimensions', DimName(1:2))
    ENDIF

    IF (PRESENT (BaselineBandNo)) THEN
       CALL SaveAsHDF5DS (sd_id, 'BaselineBandNo', BaselineBandNo)
       DimName(1) = 'BandNo'
       CALL MakeHDF5Attribute (sd_id, 'BaselineBandNo', 'Dimensions', &
            DimName(1:1))
    ENDIF

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_LatBinData

!=============================================================================
  SUBROUTINE OutputL1B_diags (sd_id, noMAF, counterMAF, MAFStartTimeTAI, &
       TP, TPdigP, TPdigN, Zeros)
!=============================================================================

    ! This subroutine writes an MAF's worth of diagnostic data

    USE MLSL1Common, ONLY: FBchans, GHzNum, MBchans, MBnum, WFchans, WFnum, &
         deflt_zero, MaxMIFs, DACSnum, DACSchans
    USE Calibration, ONLY: Chi2, Tsys, Cgain
    USE MLSHDF5, ONLY: SaveAsHDF5DS, MakeHDF5Attribute

    INTEGER, INTENT(IN) :: sd_id
    INTEGER, INTENT(IN), OPTIONAL :: counterMAF, noMAF
    REAL(r8), INTENT(IN), OPTIONAL :: MAFStartTimeTAI
    REAL, INTENT(IN), OPTIONAL :: TP(:,:), TPdigP(:,:), TPdigN(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: Zeros

    CHARACTER(LEN=16) :: DimName(2)
    INTEGER :: dims(3), status
    TYPE (DataProducts_T) :: dataset
    REAL(r8) :: MAFStartTimeGIRD

    IF (PRESENT (Zeros)) THEN
       CALL SaveAsHDF5DS (sd_id, 'DefltZeros FB', INT(deflt_zero%FB))
       DimName(1) = 'FBchan'
       DimName(2) = 'FBbank'
       CALL MakeHDF5Attribute (sd_id, 'DefltZeros FB', 'Dimensions', DimName)
       CALL SaveAsHDF5DS (sd_id, 'DefltZeros MB', INT(deflt_zero%MB))
       DimName(1) = 'MBchan'
       DimName(2) = 'MBbank'
       CALL MakeHDF5Attribute (sd_id, 'DefltZeros MB', 'Dimensions', DimName)
       CALL SaveAsHDF5DS (sd_id, 'DefltZeros WF', INT(deflt_zero%WF))
       DimName(1) = 'WFchan'
       DimName(2) = 'WFbank'
       CALL MakeHDF5Attribute (sd_id, 'DefltZeros WF', 'Dimensions', DimName)
       RETURN
    ENDIF

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

!    dims(1) = DACSchans
    dims(1) = DACSnum
    dims(2) = 1
!    dataset%Dimensions(1) = 'DACSchan            '
    dataset%Dimensions(1) = 'DACSnum             '
    dataset%name = 'Tsys DACS'

! Use channel 1 since all channels are the same (for now):

    CALL Build_MLSAuxData (sd_id, dataset, Tsys%DACS(1,:), lastIndex=noMAF, &
         dims=dims(1:2))

    dims(1) = DACSchans
    dims(2) = DACSnum
    dataset%Dimensions(1) = 'DACSchan            '
    dataset%Dimensions(2) = 'DACSnum             '
    dataset%name = 'Chi2 DACS'

    CALL Build_MLSAuxData (sd_id, dataset, Chi2%DACS, lastIndex=noMAF, &
         dims=dims)

    dims(1) = MaxMIFs - 2
    dims(2) = DACSnum
    dataset%Dimensions(1) = 'MIF                 '
    dataset%Dimensions(2) = 'DACSnum             '

    dataset%name = 'TP'
    CALL Build_MLSAuxData (sd_id, dataset, TP, lastIndex=noMAF, dims=dims)
    dataset%name = 'TPdigP'
    CALL Build_MLSAuxData (sd_id, dataset, TPdigP, lastIndex=noMAF, dims=dims)
    dataset%name = 'TPdigN'
    CALL Build_MLSAuxData (sd_id, dataset, TPdigN, lastIndex=noMAF, dims=dims)

    CALL Deallocate_DataProducts (dataset)

  END SUBROUTINE OutputL1B_diags

!=============================================================================
  SUBROUTINE OutputL1B_diagsT (sd_Id, MAFno, counterMAF, MAFStartTimeTAI, &
       nvBounds, OrbNo, Chisq, dLlo, yTsys, ColdCnts, HotCnts, LLO_Bias)
!=============================================================================

    USE MLSL1Common, ONLY: THzchans, THzNum

    INTEGER, INTENT(IN) :: sd_id
    INTEGER, INTENT(IN), OPTIONAL :: counterMAF, MAFno, nvBounds, OrbNo
    REAL(r8), INTENT (IN), OPTIONAL :: MAFStartTimeTAI
    REAL(r8), DIMENSION(:,:), INTENT (IN), OPTIONAL :: Chisq, dLlo, yTsys
    REAL, DIMENSION(:,:), INTENT (IN), OPTIONAL :: ColdCnts, HotCnts
    REAL, DIMENSION(:), INTENT (IN), OPTIONAL :: LLO_Bias

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

    IF (PRESENT (ColdCnts) .AND. PRESENT (MAFno)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'ColdCnts          '
       dataset%data_type = 'real              '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       dataset%Dimensions(3) = 'MAF'
       CALL Build_MLSAuxData (sd_Id, dataset, ColdCnts, lastIndex=MAFno, &
            dims=dims)

    ENDIF

    IF (PRESENT (HotCnts) .AND. PRESENT (MAFno)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'HotCnts           '
       dataset%data_type = 'real              '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       dataset%Dimensions(3) = 'MAF'
       CALL Build_MLSAuxData (sd_Id, dataset, HotCnts, lastIndex=MAFno, &
            dims=dims)

    ENDIF

    IF (PRESENT (LLO_Bias) .AND. PRESENT (MAFno)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'LLO_Bias          '
       dataset%data_type = 'real              '
       dataset%Dimensions(1) = 'MIFno'
       dataset%Dimensions(2) = 'MAF'
       CALL Build_MLSAuxData (sd_Id, dataset, LLO_Bias(1:148), &
            lastIndex=MAFno, dims=(/148, 1/))

    ENDIF

    IF (PRESENT (Chisq) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'Chisq             '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       dataset%Dimensions(3) = 'OrbNo'
       CALL Build_MLSAuxData (sd_Id, dataset, Chisq, lastIndex=OrbNo, dims=dims)

    ENDIF

    IF (PRESENT (dLlo) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'dLlo              '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       dataset%Dimensions(3) = 'OrbNo'
       CALL Build_MLSAuxData (sd_Id, dataset, dLlo, lastIndex=OrbNo, dims=dims)

    ENDIF
 
    IF (PRESENT (yTsys) .AND. PRESENT (OrbNo)) THEN

       DEALLOCATE (dataset%Dimensions, stat=status)
       ALLOCATE (dataset%Dimensions(3), stat=status)
       dataset%name      = 'yTsys             '
       dataset%data_type = 'double            '
       dataset%Dimensions(1) = 'THzChan'
       dataset%Dimensions(2) = 'THzBand'
       dataset%Dimensions(3) = 'OrbNo'
       CALL Build_MLSAuxData (sd_Id, dataset, yTsys, lastIndex=OrbNo, dims=dims)

    ENDIF

  END SUBROUTINE OutputL1B_diagsT

!=============================================================================
  SUBROUTINE OutputL1B_Chi2 (sd_id)
!=============================================================================

    ! This subroutine writes the default chi2 data

    USE MLSL1Common, ONLY: BandChi2, BandChans
    USE MLSHDF5, ONLY: SaveAsHDF5DS, MakeHDF5Attribute

    INTEGER, INTENT(IN) :: sd_id

    CHARACTER(LEN=16) :: DimName(2)

    CALL SaveAsHDF5DS (sd_id, 'BandChans', BandChans)
    DimName(1) = 'BandNo'
    CALL MakeHDF5Attribute (sd_id, 'BandChans', 'Dimensions', DimName(1:1))

    CALL SaveAsHDF5DS (sd_id, 'BandChi2', BandChi2)
    DimName(2) = 'ChanNo'
    CALL MakeHDF5Attribute (sd_id, 'BandChi2', 'Dimensions', DimName)

 END SUBROUTINE OutputL1B_Chi2

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE OutputL1B
!=============================================================================

! $Log$
! Revision 2.31  2016/12/21 18:59:39  whdaffer
! Just some comments to explain what GIRD time is
!
! Revision 2.30  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.29.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.29  2009/07/30 20:40:44  honghanh
! Remove a deallocate statement so that dataset sc/VelECR has a Dimensions attribute
!
! Revision 2.28  2009/06/01 14:02:34  perun
! Output galactic center latitude and longitudes for GHz and THz
!
! Revision 2.27  2008/01/15 19:55:24  perun
! Add DisableRadOut flag to stop outputting unwanted bands.
!
! Revision 2.26  2007/06/21 21:04:31  perun
! Only output to RADD file if DACS calibration is enabled and output LLO_Bias
!
! Revision 2.25  2006/09/28 16:16:36  perun
! Output ModelOffset and Poffset as just one value per MAF for the MLSL1G program
!
! Revision 2.24  2006/09/26 17:11:01  perun
! Attempt to write diagnostics offsets for GHz only!
!
! Revision 2.23  2006/09/26 16:02:27  perun
! Output DACS Chi2 values to DIAG file
!
! Revision 2.22  2006/08/02 18:56:59  perun
! Removed AscDesc dimension and added Tsys for DACS in DIAG file
!
! Revision 2.21  2006/06/14 13:48:18  perun
! Output TP data in the DIAG file
!
! Revision 2.20  2006/04/05 18:10:40  perun
! Remove unused variables
!
! Revision 2.19  2006/03/24 15:16:37  perun
! Add THz/GeodAltX to L1BOA, Poffset and ModelOffset to DIAG and ColdCnts and HotCnts to DIAGT
!
! Revision 2.18  2005/12/06 19:28:55  perun
! Added outputting BO_stat for both GHz and THz tangent point records
!
! Revision 2.17  2005/08/24 15:52:26  perun
! Output Pos_Prime for both GHz and THz in the L1BOA file
!
! Revision 2.16  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.15  2004/12/01 17:11:05  perun
! Test calibrated DACS data flag in order to output radiances
!
! Revision 2.14  2004/11/10 15:34:41  perun
! Add azimAngle fields; add baseline fields; add fills when needed
!
! Revision 2.13  2004/08/12 13:51:50  perun
! Version 1.44 commit
!
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
