

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens and reads l2cf info, NCEP, DAO and climatology files
!=============================================================================


USE MLSCommon
USE MLSL2Common
USE MLSMessageModule
USE ReadParseL2cf
!USE GriddedData
USE MLSPCF

IMPLICIT NONE
PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------
CONTAINS
SUBROUTINE OpenAndInitialize(processingRange, l1bInfo, l2cf_data, aprioriData)

! Arguments


TYPE (DataProcessingRange_T) :: processingRange ! Data processing range
TYPE (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
TYPE (L2CF) :: l2cf_data        ! The information from the l2cf file
!TYPE (GridDatabase_T) :: aprioriData ! Input a priori database


INTEGER ::  returnStatus
Integer :: L1_Version, sd_id
CHARACTER (LEN=132) :: L1physicalFilename
CHARACTER (LEN=32) :: mnemonic
CHARACTER (LEN=256) :: msg

! Open, read and parse L2cf file

Call read_parse_l2cf(mlspcf_l2cf_start, l2cf_data)


!Open L1 file
! Get the l1 file name from the PCF


!returnStatus = Pgs_pc_getReference(L1FileHandle, L1_Version, L1physicalFilename)

!IF(returnStatus /= PGS_S_SUCCESS) THEN
!  call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
!  CALL MLSMessage (MLSMSG_Error, &
!                  ModuleName, "Error opening L1:  "//mnemonic//" "//msg)
!ENDIF

! Open the HDF file and initialize the SD interface


!sd_id = sfstart(physicalFilename, DFACC_READ)
!IF (sd_id == -1) THEN
!  CALL MLSMessage (MLSMSG_Error, &
!                  ModuleName, "Error opening L1 file"//physicalFilename)
!ENDIF

!CALL Obtain_NCEP()
!CALL Obtain_DAO()
!Call Obtain_Clim ()

RETURN

END Subroutine OpenAndInitialize

!=============================================================================
END MODULE OpenInit
!=============================================================================

!
! $Log$
! Revision 1.2  2000/01/07 01:55:07  lungu
! Release containing only open l1 and obtain a priori subroutines.
!
! Revision 1.1  2000/01/05 02:12:39  lungu
! Added subroutines to handle GriddedData
!

!

