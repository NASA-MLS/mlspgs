!
! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!
MODULE OutputL1B_HDF5
  ! This module contains the subroutines needed to write L1B data to HDF5 files.

  USE OutputL1B_DataTypes, ONLY: LENUTC, L1BOAINDEX_T, L1BOASC_T, L1BOATP_T
  USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
  USE MLSAuxData, ONLY: Build_MLSAuxData, CreateGroup_MLSAuxData
  USE HDF5, ONLY: HID_T
  USE MLSCommon
  USE MLSL1Common, ONLY: L1BFileInfo_T
  USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz
  USE MLSSignalNomenclature, ONLY:GetFullMLSSignalName

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OUTPUTL1B_CREATE_HDF5, OUTPUTL1B_INDEX_HDF5, OUTPUTL1B_SC_HDF5, &
       OUTPUTL1B_GHZ_HDF5, OUTPUTL1B_THZ_HDF5, OUTPUTL1B_RAD_HDF5
  !------------------- RCS Ident Info -----------------------------------------
  CHARACTER(LEN=130) :: Id = &                                                 
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------------------------
!---------------------------------------------------------Names of Dimensions:
CONTAINS

  !----------------------------------------------------------- OutputL1B_Create
  SUBROUTINE OutputL1B_create_HDF5 (sdId)
    ! This routine assigns the groups and subgroups in the orbit and attitude 
    ! file handle, sdId%OAId
    ! Arguments
    TYPE( L1BFileInfo_T ), INTENT(IN) :: sdId

    CALL CreateGroup_MLSAuxData (sdId%OAId, 'sc')
    CALL CreateGroup_MLSAuxData (sdId%OAId, 'GHz')
    CALL CreateGroup_MLSAuxData (sdId%OAId, 'THz')

  END SUBROUTINE OutputL1B_create_HDF5

  !------------------------------------------------------------ OutputL1B_index
  SUBROUTINE OutputL1B_index_HDF5 (noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF file. 
    ! Assumes HDF5 FORTRAN APIs have been invoked to open/create files.
    ! Arguments
    TYPE( L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
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

!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_index_HDF5
!----------------------------------------------------------------- OutputL1B_sc
  SUBROUTINE OutputL1B_sc_HDF5 (noMAF, sd_id, sc)
    ! It writes the spacecraft quantities to a group in an HDF5 file.
    ! Assumes l1/Close_Files will close files.
    ! Arguments
    TYPE( L1BOAsc_T ), INTENT(IN) :: sc        
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
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
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_sc_HDF5

!---------------------------------------------------------------- OutputL1B_GHz
  SUBROUTINE OutputL1B_GHz_HDF5 (noMAF, sd_id, tp)
  ! This subroutine writes the GHz tangent point quantities to a group 
  ! in an HDF5 file.  Assumes HDF5 FORTRAN APIs have been invoked by 
  ! l1/OpenInit to open files. Assumes l1/Close_Files will close files.
  ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
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
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_GHz_HDF5

!---------------------------------------------------------------- OutputL1B_THz
  SUBROUTINE OutputL1B_THz_HDF5 (noMAF, sd_id, tp)
    ! This subroutine writes the THz tangent point quantities to a 
    ! group in an HDF5 file.
    ! Assumes HDF5 FORTRAN APIs have been invoked by l1/OpenInit to open files.
    ! Assumes l1/Close_Files will close files.
    ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
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
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_THz_HDF5

!---------------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad_HDF5 (noMAF, sdId, counterMAF, Reflec, &
       MAFStartTimeGIRD, rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files
    ! Assumes HDF5 FORTRAN APIs have been invoked by l1/OpenInit to open files.
    ! Assumes l1/Close_Files will close files.
    ! Arguments

    USE EngTbls, ONLY: Reflec_T
    USE MLSL1Rad, ONLY: Radiance_T, Rad_Name

    TYPE (L1BFileInfo_T) :: sdId
    TYPE (Radiance_T) :: rad(:)
    TYPE (Reflec_T) :: Reflec
    INTEGER, INTENT(IN) :: counterMAF, noMAF
    REAL(r8), INTENT(IN) :: MAFStartTimeGIRD

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

  END SUBROUTINE OutputL1B_rad_HDF5

END MODULE OutputL1B_HDF5

! $Log$
! Revision 2.6  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.3  2002/11/25 05:34:30  jdone
! Added dimension names.
!
! Revision 2.2  2002/11/18 21:21:30  jdone
! Used trim for names of radiance and supporting precision files.
!
! Revision 2.1  2002/11/07 21:36:57  jdone
! File holds all HDF5 output routines for Level 1
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
