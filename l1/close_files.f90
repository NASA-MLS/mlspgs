! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Close_files ! Close the production files
!=============================================================================

  USE SDPToolkit
  USE MLSCommon
  USE MLSL1Common
  USE MLSMessageModule
  USE Hdf

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CloseFiles (l0Info, l1bInfo)
!=============================================================================

    ! Arguments
    
    TYPE (L0Info_T), INTENT (IN) :: l0Info ! Lvl0 file info
    TYPE (L1BInfo_T), INTENT (IN):: l1bInfo  ! File handles etc. for L1B dataset

    INTEGER :: status

    INTEGER :: PGS_io_l0_close

    ! Close the science L0 file

    status = PGS_io_l0_close (l0Info%sci)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L0 Science file: '//l0Info%SciFileName)

    ! Close L1RAD file (eventually will do 3 files)

    status = sfend (l1bInfo%L1BRADids(1))
    IF (status == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & 'Failed to close L1BRAD file: '//l1bInfo%L1BRADFileNames(1))
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//l1bInfo%L1BRADFileNames(1))

    ! Close L1BOA file

    status = sfend (l1bInfo%L1BOAid)
    IF (status == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & 'Failed to close L1BOA file: '//l1bInfo%L1BOAFileName)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//l1bInfo%L1BOAFileName)
    
    RETURN

  END Subroutine CloseFiles

!=============================================================================
END MODULE Close_files
!=============================================================================

!
! $Log$
! Revision 1.1  2000/02/08 21:28:51  perun
! Initial version
!
!
