

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens and reads l2cf info, NCEP, DAO and climatology(ies) file
s
!=============================================================================


USE MLSCommon
USE MLSMessageModule
USE GriddedData

IMPLICIT NONE
PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------
CONTAINS

Subroutine Open_Init (returnStatus )
INTEGER ::  returnStatus
Integer :: L1_Version, sd_id
CHARACTER (LEN=135) :: L1physicalFilename
CHARACTER (LEN=32) :: mnemonic
CHARACTER (LEN=480) :: msg


!Open L1 file
! Get the l1 file name from the PCF


returnStatus = Pgs_pc_getReference(L1FileHandle, L1_Version, L1physicalFilename)

IF(returnStatus /= PGS_S_SUCCESS) THEN
  call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  CALL MLSMessage (MLSMSG_Error, &
                  ModuleName, "Error opening L1:  "//mnemonic//" "//msg)
ENDIF

! Open the HDF file and initialize the SD interface


sd_id = sfstart(physicalFilename, DFACC_READ)
IF (sd_id == -1) THEN
  CALL MLSMessage (MLSMSG_Error, &
                  ModuleName, "Error opening L1 file"//physicalFilename)
ENDIF

CALL Obtain_NCEP()
CALL Obtain_DAO()
Call Obtain_Clim ()

RETURN

END Subroutine Open_Init

!=============================================================================
END MODULE OpenInit
!=============================================================================

!
! $Log$
! Revision 1.1  2000/01/05 02:12:39  lungu
! Added subroutines to handle GriddedData
!

!

