! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SortQualify ! 
!=============================================================================

  USE MLSCommon
  USE MLSL1Common
  USE MLSMessageModule
  USE OpenInit

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
  SUBROUTINE SortAndQualify (l0Info, l1bInfo, l0Sci, MAFinfo, success)
!=============================================================================

    ! Arguments
    
    TYPE (L0Info_T), INTENT (IN) :: l0Info   ! Lvl0 file info
    TYPE (L1BInfo_T), INTENT (IN):: l1bInfo  ! File handles etc. for L1B dataset
    TYPE (L0Sci_T), INTENT (OUT) :: l0Sci(*)  ! Lvl0 science data for 1 MAF
    TYPE (MAFinfo_T), INTENT (OUT) :: MAFinfo ! MAF info
    LOGICAL, INTENT (OUT) :: success   ! successful read

    INTEGER :: MIFno, MIFs, ret_len
    DOUBLE PRECISION :: TAI93

    MIFs = 148   ! Number of MIFs to read for now!!!

    DO MIFno = 1, MIFs

       ret_len = ReadL0Packet (l0Info%sci, LEN(l0_sci), l0_sci, TAI93)

       IF (MIFno == 1) THEN
          MAFinfo%startTAI = TAI93
       END IF

       ! will put the RAW data into the L0 science data here

    END DO

    ! Hard code some info for now!!!

    MAFinfo%MIFsPerMAF = MIFs
    MAFinfo%integTime = 1.0 / 6.0

    success  = (ret_len > 0)  ! In the future, may use integer status instead

    RETURN

  END Subroutine SortAndQualify

!=============================================================================
END MODULE SortQualify
!=============================================================================

!
! $Log$
! Revision 1.1  2000/02/08 21:24:48  perun
! Initial version
!
!
