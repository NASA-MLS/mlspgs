!
! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!
MODULE OutputL1B_HDF5
  ! This module contains the subroutines needed to write L1B data to HDF files.
  USE OutputL1B_DataTypes, only: LENG, LENT, LENUTC, LENCOORD, & 
       L1BOAINDEX_T, L1BOASC_T, L1BOATP_T
  USE MLSAuxData, only:DataProducts_T, Build_MLSAuxData, CreateGroup_MLSAuxData
  USE HDF5, only: HID_T
  USE MLSCommon
  USE MLSL1Common, only: L1BFileInfo_T
  USE MLSL1Config, only: MIFsGHz, MIFsTHz
  USE MLSL1Rad, only: Radiance_T
  USE MLSSignalNomenclature, only:GetFullMLSSignalName
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: OUTPUTL1B_CREATE_HDF5, OUTPUTL1B_INDEX_HDF5, OUTPUTL1B_SC_HDF5, OUTPUTL1B_GHZ_HDF5, & 
    & OUTPUTL1B_THZ_HDF5, OUTPUTL1B_RAD_HDF5
  !------------------- RCS Ident Info -----------------------------------------
  CHARACTER(LEN=130) :: Id = &                                                 
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------------------------
  TYPE( DataProducts_T ), DIMENSION(4), PUBLIC, PARAMETER :: sds_datasets=(/&
       DataProducts_T('MAFStartTimeUTC   ','character         '), & 
       DataProducts_T('MAFStartTimeTAI   ','double            '), & 
       DataProducts_T('noMIFs            ','integer           '), & 
       DataProducts_T('counterMAF        ','integer           ') /)
  TYPE( DataProducts_T ), DIMENSION(13), PUBLIC, PARAMETER::sc_sds_datasets=(/&
       DataProducts_T('sc/ECI            ','double            '), & 
       DataProducts_T('sc/ECR            ','double            '), & 
       DataProducts_T('sc/GeocAlt        ','double            '), & 
       DataProducts_T('sc/GeocLat        ','real              '), & 
       DataProducts_T('sc/GeodAlt        ','double            '), & 
       DataProducts_T('sc/GeodLat        ','real              '), & 
       DataProducts_T('sc/Lon            ','real              '), & 
       DataProducts_T('sc/GeodAngle      ','real              '), & 
       DataProducts_T('sc/VelECI         ','double            '), & 
       DataProducts_T('sc/ypr            ','double            '), & 
       DataProducts_T('sc/yprRate        ','double            '), & 
       DataProducts_T('sc/VelECR         ','double            '), & 
       DataProducts_T('sc/OrbIncl        ','real              ') /)
  TYPE( DataProducts_T ), DIMENSION(18),PUBLIC,PARAMETER::GHz_sds_datasets=(/&
       DataProducts_T('GHz/encoderAngle  ','double            '), & 
       DataProducts_T('GHz/scAngle       ','double            '), & 
       DataProducts_T('GHz/scanAngle     ','double            '), & 
       DataProducts_T('GHz/scanRate      ','real              '), & 
       DataProducts_T('GHz/tp/ECI        ','double            '), & 
       DataProducts_T('GHz/tp/ECR        ','double            '), & 
       DataProducts_T('GHz/tp/OrbY       ','real              '), & 
       DataProducts_T('GHz/tp/GeocAlt    ','double            '), & 
       DataProducts_T('GHz/tp/GeocLat    ','real              '), & 
       DataProducts_T('GHz/tp/GeocAltRate','real              '), & 
       DataProducts_T('GHz/tp/GeodAlt    ','double            '), & 
       DataProducts_T('GHz/tp/GeodLat    ','real              '), & 
       DataProducts_T('GHz/tp/GeodAltRate','real              '), & 
       DataProducts_T('GHz/tp/Lon        ','real              '), & 
       DataProducts_T('GHz/tp/GeodAngle  ','real              '), & 
       DataProducts_T('GHz/tp/SolarTime  ','real              '), & 
       DataProducts_T('GHz/tp/SolarZenith','real              '), & 
       DataProducts_T('GHz/tp/LosAngle   ','real              ') /)
  TYPE( DataProducts_T ), DIMENSION(18),PUBLIC,PARAMETER::THz_sds_datasets=(/&
       DataProducts_T('THz/encoderAngle  ','double            '), & 
       DataProducts_T('THz/scAngle       ','double            '), & 
       DataProducts_T('THz/scanAngle     ','double            '), & 
       DataProducts_T('THz/scanRate      ','real              '), & 
       DataProducts_T('THz/tp/ECI        ','double            '), & 
       DataProducts_T('THz/tp/ECR        ','double            '), & 
       DataProducts_T('THz/tp/OrbY       ','real              '), & 
       DataProducts_T('THz/tp/GeocAlt    ','double            '), & 
       DataProducts_T('THz/tp/GeocLat    ','real              '), & 
       DataProducts_T('THz/tp/GeocAltRate','real              '), & 
       DataProducts_T('THz/tp/GeodAlt    ','double            '), & 
       DataProducts_T('THz/tp/GeodLat    ','real              '), & 
       DataProducts_T('THz/tp/GeodAltRate','real              '), & 
       DataProducts_T('THz/tp/Lon        ','real              '), & 
       DataProducts_T('THz/tp/GeodAngle  ','real              '), & 
       DataProducts_T('THz/tp/SolarTime  ','real              '), & 
       DataProducts_T('THz/tp/SolarZenith','real              '), & 
       DataProducts_T('THz/tp/LosAngle   ','real              ') /)
!---------------------------------------------------------Names of Dimensions:
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM1_NAME  = 'MAF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM2_NAME  = 'MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM3_NAME  = 'GHz/MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM4_NAME  = 'THz/MIF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM5_NAME  = 'xyz'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM6_NAME  = 'charUTC'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM7_NAME  = 'chanFB'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM8_NAME  = 'chanMB'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM9_NAME  = 'chanWF'
  CHARACTER(len=*), PUBLIC, PARAMETER :: DIM10_NAME = 'chanDACS'
!---------------------------------------------------------Names of HDF5 Groups:
  CHARACTER(len=*), PUBLIC, PARAMETER :: SC_GROUP_NAME    = 'sc'
  CHARACTER(len=*), PUBLIC, PARAMETER :: GHZ_GROUP_NAME   = 'GHz'
  CHARACTER(len=*), PUBLIC, PARAMETER :: THZ_GROUP_NAME   = 'THz'
  CHARACTER(len=*), PUBLIC, PARAMETER :: TP_SUBGROUP_NAME = 'tp'
CONTAINS
  !----------------------------------------------------------- OutputL1B_Create
  SUBROUTINE OutputL1B_create_HDF5(sdId)
    ! This subroutine creates the structure of the output files, 
    ! and names the arrays and dimensions contained within them.
    ! Assumes HDF5 FORTRAN APIs have been invoked by l1/OpenInit to open files.
    ! Arguments
    TYPE( L1BFileInfo_T ), INTENT(IN) :: sdId
    call CreateGroup_MLSAuxData( sdId%OAId, SC_GROUP_NAME )
    call CreateGroup_MLSAuxData( sdId%OAId, GHZ_GROUP_NAME )
    call CreateGroup_MLSAuxData( sdId%OAId, GHZ_GROUP_NAME, TP_SUBGROUP_NAME )
    call CreateGroup_MLSAuxData( sdId%OAId, THZ_GROUP_NAME )
    call CreateGroup_MLSAuxData( sdId%OAId, THZ_GROUP_NAME, TP_SUBGROUP_NAME )
  END SUBROUTINE OutputL1B_create_HDF5
  !------------------------------------------------------------ OutputL1B_index
  SUBROUTINE OutputL1B_index_HDF5(noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF file. 
    ! Assumes HDF5 FORTRAN APIs have been invoked to open/create files.
    ! Arguments
    TYPE( L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: noMAF
    INTEGER(HID_T), INTENT(IN) :: sd_id
!------------------------------------------------------------------------------
     call Build_MLSAuxData(sd_id,sds_datasets(1),index%MAFStartTimeUTC, & 
          noMAF, lenUTC)
     call Build_MLSAuxData(sd_id,sds_datasets(2),index%MAFStartTimeTAI,noMAF)
     call Build_MLSAuxData(sd_id,sds_datasets(3),index%noMIFs,noMAF)
     call Build_MLSAuxData(sd_id,sds_datasets(4),index%counterMAF, noMAF)
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
!------------------------------------------------------------------------------
     call Build_MLSAuxData(sd_id, sc_sds_datasets(1),  sc%scECI, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(2),  sc%scECR, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(3),  sc%scGeocAlt, noMAF) 
     call Build_MLSAuxData(sd_id, sc_sds_datasets(4),  sc%scGeocLat, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(5),  sc%scGeodAlt, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(6),  sc%scGeodLat, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(7),  sc%scLon, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(8),  sc%scGeodAngle,noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(9),  sc%scVelECI,noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(10), sc%ypr, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(11), sc%yprRate, noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(12), sc%scVelECR,noMAF)
     call Build_MLSAuxData(sd_id, sc_sds_datasets(13), sc%scOrbIncl,noMAF)
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
!------------------------------------------------------------------------------
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(1) ,tp%encoderAngle,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(2) ,tp%scAngle,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(3) ,tp%scanAngle,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(4) ,tp%scanRate,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(5) ,tp%tpECI,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(6) ,tp%tpECR,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(7) ,tp%tpOrbY,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(8) ,tp%tpGeocAlt,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(9) ,tp%tpGeocLat,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(10),tp%tpGeocAltRate,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(11),tp%tpGeodAlt,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(12),tp%tpGeodLat,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(13),tp%tpGeodAltRate,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(14),tp%tpLon,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(15),tp%tpGeodAngle,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(16),tp%tpSolarTime,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(17),tp%tpSolarZenith,noMAF)
     call Build_MLSAuxData(sd_id,GHz_sds_datasets(18),tp%tpLosAngle,noMAF)
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
!------------------------------------------------------------------------------
     call Build_MLSAuxData(sd_id, THz_sds_datasets(1) ,tp%encoderAngle,noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(2) ,tp%scAngle,noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(3) ,tp%scanAngle,noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(4) ,tp%scanRate, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(5) ,tp%tpECI, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(6) ,tp%tpECR, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(7) ,tp%tpOrbY, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(8) ,tp%tpGeocAlt, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(9) ,tp%tpGeocLat, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(10),tp%tpGeocAltRate,noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(11),tp%tpGeodAlt, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(12),tp%tpGeodLat, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(13),tp%tpGeodAltRate,noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(14),tp%tpLon, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(15),tp%tpGeodAngle, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(16),tp%tpSolarTime, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(17),tp%tpSolarZenith, noMAF)
     call Build_MLSAuxData(sd_id, THz_sds_datasets(18),tp%tpLosAngle, noMAF)
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
    INTEGER :: i
! Find the sds_id number for counterMAF in the file RADF, write the data to it,
! terminate access to the data set
!------------------------------------------------------------------------------
     call Build_MLSAuxData(sdId%RADDID, sds_datasets(4), counterMAF, noMAF)
     call Build_MLSAuxData(sdId%RADFID, sds_datasets(4), counterMAF, noMAF)
!------------------------------------------------------------------------------
    do i = 1, size(rad) ! Loop on number of SDs per MAF
      call GetFullMLSSignalName(rad(i)%signal, name) ! Concatenate SD names
      prec = TRIM(name) // ' precision'
      ! Set parameters based on input data dimensions
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
        dim_mif = DIM4_NAME
        dims(2) = MIFsTHz      !! lenT
      ELSE
        dim_mif = DIM3_NAME
        dims(2) = MIFsGHz      !! lenG
      ENDIF
!------------------------------------------------Create/Open the value datasets
      dataset%name = name
      dataset%data_type = 'real' 
      call Build_MLSAuxData(sd_id, dataset, rad(i)%value, noMAF, dims)
!------------------------------------------- Create/Open the precision datasets
      dataset%name = prec
      dataset%data_type = 'real' 
      call Build_MLSAuxData(sd_id, dataset, rad(i)%precision, noMAF, dims)
!-----------------------------------------------------------------------------
     endif
   enddo
  END SUBROUTINE OutputL1B_rad_HDF5
END MODULE OutputL1B_HDF5

! $Log$
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
