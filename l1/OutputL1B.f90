! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!=============================================================================
MODULE OutputL1B
!=============================================================================

  USE MLSCommon
  USE MLSL1Common, ONLY: L1BFileInfo_T
  USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz, L1Config
  USE MLSFiles, ONLY: HDFVERSION_4, HDFVERSION_5  
  USE MLSL1Rad, ONLY: Radiance_T
  USE OutputL1B_DataTypes
  USE OutputL1B_HDF5, ONLY: OUTPUTL1B_CREATE_HDF5, OUTPUTL1B_INDEX_HDF5, &
       OUTPUTL1B_SC_HDF5, OUTPUTL1B_GHZ_HDF5, OUTPUTL1B_THZ_HDF5, &
       OUTPUTL1B_RAD_HDF5
  USE OutputL1B_HDF4, ONLY: OUTPUTL1B_CREATE_HDF4, OUTPUTL1B_INDEX_HDF4, &
       OUTPUTL1B_SC_HDF4, OUTPUTL1B_GHZ_HDF4, OUTPUTL1B_THZ_HDF4, &
       OUTPUTL1B_RAD_HDF4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OUTPUTL1B_CREATE, OUTPUTL1B_INDEX, OUTPUTL1B_SC, OUTPUTL1B_GHZ, &
       OUTPUTL1B_THZ, OUTPUTL1B_RAD

  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------
CONTAINS
  !----------------------------------- OutputL1B_Create ----
  SUBROUTINE OutputL1B_create (sdId)
    ! This subroutine opens/creates the SD output files, and names the arrays
    ! and dimensions contained within them.
    ! Arguments
    TYPE( L1BFileInfo_T ), INTENT(IN) :: sdId
  
    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_create_HDF4 (sdId)
    ELSE
       CALL OutputL1B_create_HDF5 (sdId)
    ENDIF

  END SUBROUTINE OutputL1B_create
  !------------------------------------------------- OuptutL1B_index ----
  SUBROUTINE OutputL1B_index (noMAF, sd_id, index)
    ! This subroutine writes the time/MIF indexing quantities to the HDF-SD file
    ! Arguments
    TYPE( L1BOAindex_T), INTENT(IN) :: index
    INTEGER, INTENT(IN) :: sd_id, noMAF

    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_index_HDF4 (noMAF, sd_id, index)
    ELSE
       CALL OutputL1B_index_HDF5 (noMAF, sd_id, index)
    ENDIF

  END SUBROUTINE OutputL1B_index
  !------------------------------------------- OutputL1B_sc ------------
  SUBROUTINE OutputL1B_sc (noMAF, sd_id, sc)
    ! This subroutine writes the spacecraft quantities to the HDF-SD file.
    ! Arguments
    TYPE( L1BOAsc_T ), INTENT(IN) :: sc
    INTEGER, INTENT(IN) :: noMAF, sd_id

    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_sc_HDF4 (noMAF, sd_id, sc)
    ELSE
       CALL OutputL1B_sc_HDF5 (noMAF, sd_id, sc)
    ENDIF

  END SUBROUTINE OutputL1B_sc
  !-------------------------------------------- OutputL1B_GHz -------
  SUBROUTINE OutputL1B_GHz (noMAF, sd_id, tp)
    ! This subroutine writes the GHz tangent point quantities to the HDF-SD file
    ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id

    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_GHz_HDF4 (noMAF, sd_id, tp)
    ELSE
       CALL OutputL1B_GHz_HDF5 (noMAF, sd_id, tp)
    ENDIF

  END SUBROUTINE OutputL1B_GHz
  !-------------------------------------------- OutputL1B_THz --------------
  SUBROUTINE OutputL1B_THz (noMAF, sd_id, tp)
    ! This subroutine writes the THz tangent point quantities to the HDF-SD file
    ! Arguments
    TYPE( L1BOAtp_T ), INTENT(IN) :: tp
    INTEGER, INTENT(IN) :: noMAF, sd_id

    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_THz_HDF4 (noMAF, sd_id, tp)
    ELSE
       CALL OutputL1B_THz_HDF5 (noMAF, sd_id, tp)
    ENDIF

  END SUBROUTINE OutputL1B_THz
  !-------------------------------------------------------- OutputL1B_rad
  SUBROUTINE OutputL1B_rad (noMAF, sdId, counterMAF, rad)
    ! This subroutine writes an MAF's worth of data to the L1BRad D & F files
    ! Arguments
    TYPE( L1BFileInfo_T ) :: sdId
    TYPE( Radiance_T ) :: rad(:)
    INTEGER, INTENT(IN) :: counterMAF, noMAF

    IF (L1Config%Output%HDFversion == HDFVERSION_4) THEN
       CALL OutputL1B_rad_HDF4 (noMAF, sdId, counterMAF, rad)
    ELSE
       CALL OutputL1B_rad_HDF5 (noMAF, sdId, counterMAF, rad)
    ENDIF

  END SUBROUTINE OutputL1B_rad
!=============================================================================
END MODULE OutputL1B
!=============================================================================

! $Log$
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

