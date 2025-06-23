! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
MODULE OutputL1B_DataTypes

  USE MLSCommon

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LENG, LENT, LENCOORD, LENUTC, L1BOAINDEX_T, L1BOASC_T, L1BOATP_T
!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------
  INTEGER, PARAMETER :: lenCoord =   3
  INTEGER, PARAMETER :: lenUTC   =  27
  INTEGER, PARAMETER :: lenG     = 125
  INTEGER, PARAMETER :: lenT     = lenG         ! same size for Level 2

  ! This data type contains index information for the L1BOA data file.

  TYPE L1BOAindex_T
    CHARACTER(LEN=lenUTC) :: MAFStartTimeUTC    ! MAF start time in ascii UTC
    REAL(r8) :: MAFStartTimeTAI                 ! MAF start time in TAI
    INTEGER  :: noMIFs                          ! total number of MIFs per MAF
    INTEGER  :: counterMAF                      ! MAF counter
  END TYPE L1BOAindex_T

  ! This data type contains spacecraft quantities for the L1BOA data file.

  TYPE L1BOAsc_T
    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scECI => NULL()   ! s/c ECI location
    REAL(r8), DIMENSION(:,:), POINTER :: scECR => NULL()   ! s/c ECR location
    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:), POINTER :: MIF_TAI => NULL()   ! TAI time per MIF
    REAL(r8), DIMENSION(:), POINTER :: scGeocAlt => NULL() ! s/c geoc alt
    REAL(r8), DIMENSION(:), POINTER :: scGeodAlt => NULL() ! s/c geod alt
    REAL,     DIMENSION(:), POINTER :: scGeocLat => NULL() ! s/c geoc lat
    REAL,     DIMENSION(:), POINTER :: scGeodLat => NULL() ! s/c geod lat
    REAL,     DIMENSION(:), POINTER :: scLon => NULL()     ! s/c longitude
    REAL,     DIMENSION(:), POINTER :: scGeodAngle => NULL() ! s/c master coord
    REAL,     DIMENSION(:), POINTER :: scOrbIncl => NULL() ! instant.orb incl ECR(deg)
    ! dimensioned (xyz,MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECI => NULL()  ! s/c ECI velocity
    REAL(r8), DIMENSION(:,:), POINTER :: scVelECR => NULL()  ! s/c ECR velocity
    REAL(r8), DIMENSION(:,:), POINTER :: ypr => NULL()     ! s/c yaw, pitch, roll
    REAL(r8), DIMENSION(:,:), POINTER :: yprRate => NULL()   ! s/c y-p-r rate
  END TYPE L1BOAsc_T

  ! This data type contains tangent point quantities for the L1BOA data file.

  TYPE L1BOAtp_T
    ! dimensioned (MIF)
    REAL(r8), DIMENSION(:), POINTER :: encoderAngle => NULL() ! boresight wrt instr.
    REAL(r8), DIMENSION(:), POINTER :: scAngle => NULL() ! boresight wrt s/c +x
    REAL(r8), DIMENSION(:), POINTER :: scanAngle  => NULL() ! boresight wrt orbit +x
    REAL,     DIMENSION(:), POINTER :: scanRate  => NULL() ! rate of scanAngle
    REAL(r8), DIMENSION(:), POINTER :: azimAngle  => NULL() ! azimuth wrt orbit +x
    ! dimensioned (xyz,mod.MIF)
    REAL(r8), DIMENSION(:,:), POINTER :: tpECI => NULL()    ! tp location in ECI
    REAL(r8), DIMENSION(:,:), POINTER :: tpECR => NULL()    ! tp location in ECR
    REAL(r8), DIMENSION(:,:), POINTER :: tpECRtoFOV => NULL() ! tp ECR to FOV
    REAL, DIMENSION(:,:), POINTER :: tpPos_P => NULL()  ! tp Pos Prime vals
    ! dimensioned (mod.MIF)
    REAL(r8), DIMENSION(:), POINTER :: tpGeodAlt => NULL()  ! geod alt of tp 
    REAL(r8), DIMENSION(:), POINTER :: tpGeodAltX => NULL() ! geod alt extension
    REAL(r8), DIMENSION(:), POINTER :: tpGeocAlt  => NULL() ! geoc alt of tp
    REAL,     DIMENSION(:), POINTER :: tpOrbY => NULL()     ! out-of-plane dist.
    REAL,     DIMENSION(:), POINTER :: tpGeocLat  => NULL() ! geoc lat of tp
    REAL,     DIMENSION(:), POINTER :: tpGeocAltRate => NULL() ! rate of tpGeocAlt
    REAL,     DIMENSION(:), POINTER :: tpGeodLat => NULL()     ! geod lat of tp
    REAL,     DIMENSION(:), POINTER :: tpGeodAltRate => NULL() ! rate of tpGeodAlt
    REAL,     DIMENSION(:), POINTER :: tpLon => NULL()       ! longitude of tp
    REAL,     DIMENSION(:), POINTER :: tpGeodAngle => NULL() ! tp master coord
    REAL,     DIMENSION(:), POINTER :: tpSolarTime => NULL() ! solar time coord
    REAL,     DIMENSION(:), POINTER :: tpSolarZenith => NULL() ! solar zenith ang
    REAL,     DIMENSION(:), POINTER :: tpLosAngle => NULL() ! line-of-sight ang to N
    REAL(r8), DIMENSION(:), POINTER :: tpLosVel => NULL()   ! line-of-sight velocity
    INTEGER,  DIMENSION(:), POINTER :: tpBO_stat => NULL()  ! Bright Obj Status of tp
    REAL,     DIMENSION(:), POINTER :: tpGalLat => NULL()   ! galactic cntr lat of tp
    REAL,     DIMENSION(:), POINTER :: tpGalLon => NULL()   ! galactic cntr long of tp
  END TYPE L1BOAtp_T

CONTAINS
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE OutputL1B_DataTypes

! $Log$
! Revision 2.9  2014/09/29 22:49:01  vsnyder
! Remove nonstandard tab formatting
!
! Revision 2.8  2009/06/01 14:03:40  perun
! Added galactic center for each tp record structure
!
! Revision 2.7  2006/03/24 15:16:01  perun
! Add tpGeodAltX for BOA to extend the altitude range for the THz coverage
!
! Revision 2.6  2005/12/06 19:28:30  perun
! Added BO_stat for each tp record structure
!
! Revision 2.5  2005/08/24 15:51:55  perun
! Add pointer for Pos_P to tangent point structure
!
! Revision 2.4  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2004/11/10 15:38:15  perun
! Add azimAngle to output
!
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
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

