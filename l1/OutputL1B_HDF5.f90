!
! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!
MODULE OutputL1B_HDF5
  ! This module contains the subroutines needed to write L1B data to HDF files.
  USE OutputL1B_DataTypes, only: LENG, LENT, LENUTC, LENCOORD, & 
       L1BOAINDEX_T, L1BOASC_T, L1BOATP_T
  USE MLS_DataProducts, only: DataProducts_T, Deallocate_DataProducts
  USE MLSAuxData, only: Build_MLSAuxData, CreateGroup_MLSAuxData
  USE HDF5, only: HID_T
  USE MLSCommon
  USE MLSL1Common, only: L1BFileInfo_T
  USE MLSL1Config, only: MIFsGHz, MIFsTHz
  USE MLSL1Rad, only: Radiance_T
  USE MLSSignalNomenclature, only:GetFullMLSSignalName
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
  SUBROUTINE OutputL1B_create_HDF5(sdId)
    ! This routine assigns the groups and subgroups in the orbit and attitude 
    ! file handle, sdId%OAId
    ! Arguments
    TYPE( L1BFileInfo_T ), INTENT(IN) :: sdId
    call CreateGroup_MLSAuxData( sdId%OAId, 'sc')
    call CreateGroup_MLSAuxData( sdId%OAId, 'GHz')
    call CreateGroup_MLSAuxData( sdId%OAId, 'THz')
  END SUBROUTINE OutputL1B_create_HDF5
  !------------------------------------------------------------ OutputL1B_index
  SUBROUTINE OutputL1B_index_HDF5(noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF file. 
    ! Assumes HDF5 FORTRAN APIs have been invoked to open/create files.
    ! Arguments
    TYPE( L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
! Variables
    TYPE( DataProducts_T ) :: dataset
    integer :: status
!------------------------------------------------------------------------------
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'MAFStartTimeUTC   '
     dataset%data_type = 'character         '
     allocate(dataset%Dimensions(1), stat=status)
     dataset%Dimensions(1) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,index%MAFStartTimeUTC, & 
          char_length=lenUTC, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'MAFStartTimeTAI   '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(1), stat=status)
     dataset%Dimensions(1) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,index%MAFStartTimeTAI,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'noMIFs            '
     dataset%data_type = 'integer           '
     allocate(dataset%Dimensions(1), stat=status)
     dataset%Dimensions(1) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,index%noMIFs,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'counterMAF        '
     dataset%data_type = 'integer           '
     allocate(dataset%Dimensions(1), stat=status)
     dataset%Dimensions(1) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,index%counterMAF, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_index_HDF5
!----------------------------------------------------------------- OutputL1B_sc
  SUBROUTINE OutputL1B_sc_HDF5(noMAF, sd_id, sc)
    ! It writes the spacecraft quantities to a group in an HDF5 file.
    ! Assumes l1/Close_Files will close files.
    ! Arguments
    TYPE( L1BOAsc_T ), INTENT(IN) :: sc        
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
    ! Variables
    TYPE( DataProducts_T ) :: dataset
    integer :: status
!------------------------------------------------------------------------------
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/ECI            '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scECI, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/ECR            '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scECR, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/GeocAlt        '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                '
     dataset%Dimensions(2) = 'MAF                '
     call Build_MLSAuxData(sd_id, dataset, sc%scGeocAlt, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/GeocLat        '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                '
     dataset%Dimensions(2) = 'MAF                '
     call Build_MLSAuxData(sd_id, dataset, sc%scGeocLat, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/GeodAlt        '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scGeodAlt, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/GeodLat        '
     dataset%data_type = 'real              '    
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '          
     call Build_MLSAuxData(sd_id, dataset, sc%scGeodLat, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/Lon            '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '          
     call Build_MLSAuxData(sd_id, dataset, sc%scLon, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/GeodAngle      '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '          
     call Build_MLSAuxData(sd_id, dataset, sc%scGeodAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/VelECI         '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scVelECI,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/ypr            '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%ypr, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/yprRate        '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%yprRate, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/VelECR         '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'MIF                 '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scVelECR,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'sc/OrbIncl        '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id, dataset, sc%scOrbIncl,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_sc_HDF5
  !-------------------------------------------------------------- OutputL1B_GHz
  SUBROUTINE OutputL1B_GHz_HDF5(noMAF, sd_id, tp)
  ! This subroutine writes the GHz tangent point quantities to a group 
  ! in an HDF5 file.  Assumes HDF5 FORTRAN APIs have been invoked by 
  ! l1/OpenInit to open files. Assumes l1/Close_Files will close files.
  ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
    ! Variables
    TYPE( DataProducts_T ) :: dataset
    integer :: status
!------------------------------------------------------------------------------
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/encoderAngle  '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%encoderAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/scAngle       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%scAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/scanAngle     '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%scanAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/scanRate      '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF                 '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%scanRate,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/ECI           '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'GHz.MIF             '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpECI,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/ECR           '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz                 '
     dataset%Dimensions(2) = 'GHz.MIF             '
     dataset%Dimensions(3) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpECR,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/OrbY          '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF              '
     dataset%Dimensions(2) = 'MAF                  '
     call Build_MLSAuxData(sd_id,dataset,tp%tpOrbY,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeocAlt       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF              '
     dataset%Dimensions(2) = 'MAF                  '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeocAlt,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeocLat       '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeocLat,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeocAltRate   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeocAltRate,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeodAlt       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeodAlt,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeodLat       '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeodLat,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeodAltRate   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeodAltRate,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/Lon           '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpLon,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/GeodAngle     '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF             '
     dataset%Dimensions(2) = 'MAF                 '
     call Build_MLSAuxData(sd_id,dataset,tp%tpGeodAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/SolarTime     '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF            '
     dataset%Dimensions(2) = 'MAF                '
     call Build_MLSAuxData(sd_id,dataset,tp%tpSolarTime,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/SolarZenith   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF            '
     dataset%Dimensions(2) = 'MAF                '
     call Build_MLSAuxData(sd_id,dataset,tp%tpSolarZenith,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'GHz/LosAngle      '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'GHz.MIF            '
     dataset%Dimensions(2) = 'MAF                '
     call Build_MLSAuxData(sd_id,dataset,tp%tpLosAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_GHz_HDF5
!---------------------------------------------------------------- OutputL1B_THz
  SUBROUTINE OutputL1B_THz_HDF5(noMAF, sd_id, tp)
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
    integer :: status
!------------------------------------------------------------------------------
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/encoderAngle  '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%encoderAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/scAngle       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%scAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/scanAngle     '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%scanAngle,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/scanRate      '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%scanRate, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/ECI           '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz'
     dataset%Dimensions(2) = 'THz.MIF'
     dataset%Dimensions(3) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpECI, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/ECR           '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(3), stat=status)
     dataset%Dimensions(1) = 'xyz'
     dataset%Dimensions(2) = 'THz.MIF'
     dataset%Dimensions(3) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpECR, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/OrbY          '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpOrbY, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeocAlt       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeocAlt, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeocLat       '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeocLat, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeocAltRate   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeocAltRate,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeodAlt       '
     dataset%data_type = 'double            '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeodAlt, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeodLat       '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeodLat, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeodAltRate   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeodAltRate,lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/Lon           '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpLon, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/GeodAngle     '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpGeodAngle, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/SolarTime     '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpSolarTime, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/SolarZenith   '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpSolarZenith, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'THz/LosAngle      '
     dataset%data_type = 'real              '
     allocate(dataset%Dimensions(2), stat=status)
     dataset%Dimensions(1) = 'THz.MIF'
     dataset%Dimensions(2) = 'MAF' 
     call Build_MLSAuxData(sd_id, dataset,tp%tpLosAngle, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
!------------------------------------------------------------------------------
  END SUBROUTINE OutputL1B_THz_HDF5
!---------------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad_HDF5(noMAF, sdId, counterMAF, rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files
    ! Assumes HDF5 FORTRAN APIs have been invoked by l1/OpenInit to open files.
    ! Assumes l1/Close_Files will close files.
    ! Arguments
    TYPE( L1BFileInfo_T ) :: sdId
    TYPE( Radiance_T    ) :: rad(:)
    INTEGER, INTENT(IN) :: counterMAF, noMAF
    ! Variables
    TYPE( DataProducts_T ) :: dataset    
    CHARACTER (LEN=64) :: dim_chan, dim_mif, name, prec
    INTEGER, DIMENSION(3) :: dims
    INTEGER(HID_T) :: sd_id
    INTEGER :: i, status
! Find the sds_id number for counterMAF in the file RADF, write the data to it,
! terminate access to the data set
!------------------------------------------------------------------------------
     call Deallocate_DataProducts(dataset)
     dataset%name      = 'counterMAF        '
     dataset%data_type = 'integer           '
     allocate(dataset%Dimensions(1), stat=status)
     dataset%Dimensions(1) = 'MAF'
     call Build_MLSAuxData(sdId%RADDID, dataset, counterMAF, lastIndex=noMAF)
     call Build_MLSAuxData(sdId%RADFID, dataset, counterMAF, lastIndex=noMAF)
     call Deallocate_DataProducts(dataset)
!------------------------------------------------------------------------------
    do i = 1, size(rad) ! Loop on number of SDs per MAF
      call GetFullMLSSignalName(rad(i)%signal, name) ! Concatenate SD names
      prec = trim(name) // ' precision'
      ! Set parameters based on input data dimensions
      ! Based on the SD name, set dim name for channel, get Id of output file
      IF ( INDEX(name,'FB') /= 0 ) THEN
        dim_chan = 'chanFB              '
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'MB') /= 0 ) THEN
        dim_chan = 'chanMB              '
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'WF') /= 0 ) THEN
        dim_chan = 'chanWF              '
        sd_id = sdId%RADFID
      ELSE IF ( INDEX(name,'DACS') /= 0 ) THEN
        dim_chan = 'chanDACS            '
        sd_id = sdId%RADDID
      ELSE
        dim_chan = ''
        sd_id = -999
      ENDIF
     IF ((any(rad(i)%value /= 0.0)).and.(sd_id /= -999)) THEN   
                                            ! if good radiance data exists 
                                            ! Create/Open the value SDs
      dims(1) = SIZE(rad(i)%value,1)
      dims(2) = SIZE(rad(i)%value,2)
      dims(3) = 1
      IF ( INDEX(name,'R5') /= 0 ) THEN
        dim_mif = 'THz.MIF              '
        dims(2) = MIFsTHz      !! lenT
      ELSE
        dim_mif = 'GHz.MIF              '
        dims(2) = MIFsGHz      !! lenG
      ENDIF
!------------------------------------------------Create/Open the value datasets
      call Deallocate_DataProducts(dataset)
      dataset%name = trim(name)
      dataset%data_type = 'real'
      allocate(dataset%Dimensions(3), stat=status)
      dataset%Dimensions(1) = dim_chan
      dataset%Dimensions(2) = dim_mif
      dataset%Dimensions(3) = 'MAF                 '
      call Build_MLSAuxData(sd_id, dataset, rad(i)%value, lastIndex=noMAF,&
           dims=dims)
      call Deallocate_DataProducts(dataset)
!------------------------------------------- Create/Open the precision datasets
      call Deallocate_DataProducts(dataset)
      dataset%name = trim(prec)
      dataset%data_type = 'real' 
      allocate(dataset%Dimensions(3), stat=status)
      dataset%Dimensions(1) = dim_chan
      dataset%Dimensions(2) = dim_mif
      dataset%Dimensions(3) = 'MAF                 '
      call Build_MLSAuxData(sd_id, dataset, rad(i)%precision, lastIndex=noMAF,&
           dims=dims)
      call Deallocate_DataProducts(dataset)
!-----------------------------------------------------------------------------
     endif
   enddo
  END SUBROUTINE OutputL1B_rad_HDF5
END MODULE OutputL1B_HDF5

! $Log$
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
