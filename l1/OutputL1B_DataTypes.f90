! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
MODULE OutputL1B_DataTypes
  USE MLSCommon
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: LENG, LENT, LENCOORD, LENUTC, L1BOAINDEX_T, L1BOASC_T, L1BOATP_T
  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &                                                    
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------
  INTEGER, PARAMETER :: lenCoord =   3
  INTEGER, PARAMETER :: lenUTC   =  27
  INTEGER, PARAMETER :: lenG     = 120
  INTEGER, PARAMETER :: lenT     = 120
  ! This data type contains index information for the L1BOA data file.
  TYPE L1BOAindex_T
    CHARACTER(LEN=lenUTC) :: MAFStartTimeUTC    ! MAF start time in ascii UTC
    REAL(r8) :: MAFStartTimeTAI		        ! MAF start time in TAI
    INTEGER  :: noMIFs                          ! total number of MIFs per MAF
    INTEGER  :: counterMAF                      ! MAF counter
  END TYPE L1BOAindex_T
  ! This data type contains spacecraft quantities for the L1BOA data file.
  TYPE L1BOAsc_T
    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scECI	        ! s/c ECI location
    REAL(r8), DIMENSION(:,:), POINTER :: scECR	        ! s/c ECR location
    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:),   POINTER :: scGeocAlt   ! s/c geoc alt
    REAL(r8), DIMENSION(:),   POINTER :: scGeodAlt   ! s/c geod alt
    REAL,     DIMENSION(:),   POINTER :: scGeocLat   ! s/c geoc lat
    REAL,     DIMENSION(:),   POINTER :: scGeodLat   ! s/c geod lat
    REAL,     DIMENSION(:),   POINTER :: scLon	     ! s/c longitude
    REAL,     DIMENSION(:),   POINTER :: scGeodAngle ! s/c master coordinate
    REAL,     DIMENSION(:),   POINTER :: scOrbIncl  ! instant.orb incl ECR(deg)
    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECI	! s/c ECI velocity
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECR       ! s/c ECR velocity
    REAL(r8), DIMENSION(:,:), POINTER :: ypr		! s/c yaw, pitch, roll
    REAL(r8), DIMENSION(:,:), POINTER :: yprRate	! s/c y-p-r rate
  END TYPE L1BOAsc_T
  ! This data type contains tangent point quantities for the L1BOA data file.
  TYPE L1BOAtp_T
    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:),  POINTER :: encoderAngle  ! boresight wrt instr.
    REAL(r8), DIMENSION(:),  POINTER :: scAngle	      ! boresight wrt s/c +x
    REAL(r8), DIMENSION(:),  POINTER :: scanAngle     ! boresight wrt orbit +x
    REAL,     DIMENSION(:),  POINTER :: scanRate      ! of change of scanAngle
    ! dimensioned (xyz,mod.MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: tpECI	      ! tp location in ECI
    REAL(r8), DIMENSION(:,:), POINTER :: tpECR	      ! tp location in ECR
    ! dimensioned (mod.MIF)
    REAL(r8), DIMENSION(:), POINTER :: tpGeodAlt       ! geod alt of tp 
    REAL(r8), DIMENSION(:), POINTER :: tpGeocAlt       ! geoc alt of tp
    REAL,     DIMENSION(:), POINTER :: tpOrbY	       ! out-of-plane distance
    REAL,     DIMENSION(:), POINTER :: tpGeocLat       ! geoc lat of tp
    REAL,     DIMENSION(:), POINTER :: tpGeocAltRate   ! of change of tpGeocAlt
    REAL,     DIMENSION(:), POINTER :: tpGeodLat       ! geod lat of tp
    REAL,     DIMENSION(:), POINTER :: tpGeodAltRate   ! of change of tpGeodAlt
    REAL,     DIMENSION(:), POINTER :: tpLon	       ! longitude of tp
    REAL,     DIMENSION(:), POINTER :: tpGeodAngle     ! tp master coordinate
    REAL,     DIMENSION(:), POINTER :: tpSolarTime     ! solar time coordinate
    REAL,     DIMENSION(:), POINTER :: tpSolarZenith   ! solar zenith angle
    REAL,     DIMENSION(:), POINTER :: tpLosAngle      ! line-of-sight ang to N
  END TYPE L1BOAtp_T
END MODULE OutputL1B_DataTypes

! $Log$
! Revision 2.1  2002/11/07 21:27:13  jdone
! Consolidate datatypes used by the output routines.
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

