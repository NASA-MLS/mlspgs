
! -------------------------------------------------------
MODULE WriteMetadata ! Populate metadata and write it out
! -------------------------------------------------------
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE L2GPData
USE L2AUXData
USE SDPToolkit
USE MLSPCF
USE MLSCF 
USE Hdf 
USE HDFEOS 

IMPLICIT NONE

PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------


CONTAINS

SUBROUTINE populate_metadata(HDF_FILE, MCF_FILE)

!Arguments

INTEGER :: HDF_FILE, MCF_FILE

!Local Variables 
INTEGER :: pgs_met_init
INTEGER :: pgs_met_setattr_d
INTEGER :: pgs_met_setAttr_s
INTEGER :: pgs_met_getsetattr_d
INTEGER :: pgs_met_getsetattr_s
INTEGER :: pgs_met_write
INTEGER :: pgs_met_getpcattr_d
INTEGER :: pgs_met_getconfigdata_d
INTEGER :: pgs_met_getconfigdata_s
INTEGER :: pgs_met_remove
INTEGER :: hdfReturn
INTEGER :: returnStatus
INTEGER :: sdid
INTEGER :: access
INTEGER :: sfstart
INTEGER :: sfend

REAL(r8) dval
INTEGER, PARAMETER :: INVENTORY=2,ARCHIVE=1
CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename, hdfeos_filename 
CHARACTER (LEN=27) :: datetime, sval
CHARACTER (LEN=32) :: mnemonic, msg
INTEGER :: version
! the group have to be defined as 49 characters long. The C interface is 50.

! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
! therefore making the actual length of a 
! string pass of 50.

CHARACTER (LEN = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)
!Executable code

version = 1
returnStatus = PGS_PC_GetReference(HDF_FILE, version , physical_filename)

returnStatus = pgs_met_init(MCF_FILE, groups)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Initialization error.") 
ENDIF


sval = 'Horizontal & Vertical'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), &
                         "GranuleSpatialDomainType", sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting GranuleSpatialDomainType")
ENDIF

sval = '2000-01-01'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), &
                         "RangeBeginningDate", sval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeBeginningDate")
ENDIF

sval= '00:00:00.000000'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), &
                         "RangeBeginningTime", sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeBeginningTime")
ENDIF

sval = '2000-01-01'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), "RangeEndingDate", &
                          sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeEndingDate")
ENDIF

sval= '23:59:59.999999'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), "RangeEndingTime", &
                          sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeEndingTime")
ENDIF


sval = 'Other Grid System'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), "ZoneIdentifier", &
                          sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting ZoneIdentifier")
ENDIF



dval = -180.0

returnStatus = pgs_met_setattr_d(groups(INVENTORY),& 
                                "WestBoundingCoordinate",dval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting WestBoundingCoordinate")
ENDIF
 
dval = -180.0
returnStatus = PGS_MET_GetSetAttr_d(groups(INVENTORY),&
                                  "WestBoundingCoordinate",dval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "WestBoundingCoordinate") 
ENDIF 

dval = 90.0

returnStatus = pgs_met_setattr_d(groups(INVENTORY),&
                                "NorthBoundingCoordinate",dval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "NorthBoundingCoordinate.") 
ENDIF

dval = 90.0

returnStatus = PGS_MET_GetSetAttr_d(groups(INVENTORY),&
                                   "NorthBoundingCoordinate",dval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "NorthBoundingCoordinate") 
ENDIF

dval = 180.0
returnStatus = pgs_met_setattr_d(groups(INVENTORY),&
                                "EastBoundingCoordinate",dval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "EastBoundingCoordinate") 
ENDIF
dval = 180.0
returnStatus = PGS_MET_GetSetAttr_d(groups(INVENTORY),&
                                   "EastBoundingCoordinate",dval)
IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "EastBoundingCoordinate") 
ENDIF

dval =-90.0 
returnStatus = pgs_met_setattr_d(groups(INVENTORY),&
&                               "SouthBoundingCoordinate",dval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "SouthBoundingCoordinate") 
ENDIF

dval = -90.0
returnStatus = PGS_MET_GetSetAttr_d(groups(INVENTORY),&
                                   "SouthBoundingCoordinate",dval)

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "SouthBoundingCoordinate") 
ENDIF

sdid = sfstart(physical_fileName, DFACC_RDWR) 

IF (sdid == -1) THEN
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Failed to open the hdf file" ) 
ENDIF

returnStatus = pgs_met_write(groups(INVENTORY),&
                            "coremetadata.0", sdid)

IF (returnStatus /= PGS_S_SUCCESS .AND. &
    returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
   IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
     CALL MLSMessage (MLSMSG_WARNING, ModuleName, "Some of the mandatory parameters were not set" )
   ELSE 
     CALL Pgs_smf_getMsg(returnStatus, mnemonic, msg)
     CALL MLSMessage (MLSMSG_WARNING, ModuleName,"Metadata write failed "//mnemonic// msg) 
   ENDIF
ENDIF

hdfReturn = sfend(sdid)

returnStatus=pgs_met_remove() 
END SUBROUTINE populate_metadata

END MODULE WriteMetadata 
! $Log$
! Revision 2.0  2000/09/05 18:57:07  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.4  2000/06/30 00:16:41  lungu
! Made dval REAL(r8).
!
