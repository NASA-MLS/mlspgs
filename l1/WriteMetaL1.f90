
! -------------------------------------------------------
MODULE WriteMetaL1 ! Populate metadata and write it out
! -------------------------------------------------------
USE MLSCommon
USE MLSStrings 
USE MLSMessageModule
USE SDPToolkit
USE MLSPCF1
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

SUBROUTINE populate_metadata_l1(HDF_FILE, MCF_FILE)
USE InitPCFs, ONLY: L1PCF

!Arguments

INTEGER :: HDF_FILE, MCF_FILE

!Local Variables 
INTEGER :: pgs_met_init
INTEGER :: pgs_met_setattr_d
INTEGER :: pgs_met_setAttr_s
INTEGER :: pgs_met_getsetattr_d
INTEGER :: PGS_MET_SETATTR_I
INTEGER :: pgs_met_write
INTEGER :: pgs_met_remove
INTEGER :: hdfReturn
INTEGER :: returnStatus
INTEGER :: sdid
INTEGER :: sfstart
INTEGER :: sfend

REAL(r8) dval
INTEGER, PARAMETER :: INVENTORY=2,ARCHIVE=1
CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
CHARACTER (LEN=132) :: attrname, errmsg
INTEGER :: version
CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = 'Error writing metadata attribute '

! the group have to be defined as 49 characters long. The C interface is 50.
! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
! therefore making the actual length of a 
! string pass of 50.

CHARACTER (LEN = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

!Executable code

version = 1
returnStatus = PGS_PC_GetReference(HDF_FILE, version , physical_filename)

returnStatus = pgs_met_init(MCF_FILE, groups)
if (returnStatus == PGSMET_E_LOAD_ERR)print *, 'PGSMET_E_LOAD_ERR'
if (returnStatus == PGSMET_E_GRP_ERR)print *, 'PGSMET_E_GRP_ERR'
if (returnStatus == PGSMET_E_GRP_NAME_ERR)print *, 'PGSMET_E_GRP_NAME_ERR'
if (returnStatus == PGSMET_E_NO_INVENT_DATA)print *,'PGSMET_E_NO_INVENT_DATA'
if (returnStatus == PGSMET_E_DUPLICATE_ERR)print *,'PGSMET_E_DUPLICATE_ERR'
if (returnStatus == PGSMET_E_NUMOFMCF_ERR)print *, 'PGSMET_E_NUMOFMCF_ERR'
if (returnStatus == PGSMET_E_PCF_VALUE_ERR)print *,'PGSMET_E_PCF_VALUE_ERR'

IF (returnStatus /= PGS_S_SUCCESS) THEN 
     CALL MLSMessage (MLSMSG_Error, ModuleName, "Initialization error.") 
ENDIF

! Set PGE values 

! ECSDataGranule

attrName = 'ReprocessingPlanned'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
               'not known at this time')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'ReprocessingActual'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                        'first processing')
IF (returnStatus /= PGS_S_SUCCESS) THEN
errmsg = METAWR_ERR // attrName
CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'LocalGranuleID'
sval = physical_filename
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, sval)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'DayNightFlag'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, 'DayAndNight')
IF (returnStatus /= PGS_S_SUCCESS) THEN
  errmsg = METAWR_ERR // attrName
  CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'LocalVersionID'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, '0.5')
   IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

! MeasuredParameterContainer


attrName = 'ParameterName' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                           'Calibrated Radiances')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

! QAFlags Group

attrName = 'AutomaticQualityFlag' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                           'Passed')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'AutomaticQualityFlagExplanation' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                           'TBD')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'OperationalQualityFlag' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                           'Not Investigated')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'OperationalQualityFlagExplanation' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                                  'TBD')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

! QAStats Group

attrName = 'QAPercentInterpolatedData' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, 0)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'QAPercentMissingData' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, 0)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'QAPercentOutofBoundsData' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, 0)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF


! Orbit Calculated Spatial Domain Container

attrName = 'OrbitNumber' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, -1)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'StartOrbitNumber' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, -1)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'StopOrbitNumber' // '.1'
returnStatus = pgs_met_setAttr_i(groups(INVENTORY), attrName, -1)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'EquatorCrossingLongitude' // '.1'
dval = 0.0
returnStatus = pgs_met_setAttr_d(groups(INVENTORY), attrName, dval)
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'EquatorCrossingTime' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
     '00:00:00')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

attrName = 'EquatorCrossingDate' // '.1'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
     '1969-04-04')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF


! InputPointer

attrName = 'InputPointer'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
                       'See the PCF annotation to this file.')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF

IF (MCF_FILE /= mlspcf_mcf_l1boa_start) THEN
! Locality Value

   attrName = 'LocalityValue'
   returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
        'Limb')
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   ! VerticalSpatialDomain Product-Specific Attribute

   attrName = 'VerticalSpatialDomainType' // '.1'
   returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
        'Height')
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   attrName = 'VerticalSpatialDomainValue' // '.1'
   returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
        'Brigthness Temperature')
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   ! HorizontalSpatialDomainContainer

   attrName = 'ZoneIdentifier'
   returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, &
        'Global')
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   attrName = 'WestBoundingCoordinate'
   dval = -180.0
   returnStatus = pgs_met_setAttr_d(groups(INVENTORY), attrName, dval)
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   attrName = 'NorthBoundingCoordinate'
   dval = 90.0
   returnStatus = pgs_met_setAttr_d(groups(INVENTORY), attrName, dval)
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   attrName = 'EastBoundingCoordinate'
   dval = 180.0
   returnStatus = pgs_met_setAttr_d(groups(INVENTORY), attrName, dval)
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF

   attrName = 'SouthBoundingCoordinate'
   dval = -90.0
   returnStatus = pgs_met_setAttr_d(groups(INVENTORY), attrName, dval)
   IF (returnStatus /= PGS_S_SUCCESS) THEN
      errmsg = METAWR_ERR // attrName
      CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
   ENDIF
ENDIF

! RangeDateTime Group
attrName = 'RangeBeginningDate'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), &
                         attrName, L1PCF%startUTC(1:10))
IF (returnStatus /= PGS_S_SUCCESS) THEN 
   CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeBeginningDate")
ENDIF

returnStatus = pgs_met_setAttr_s(groups(INVENTORY), &
                         "RangeBeginningTime", L1PCF%startUTC(12:26))
IF (returnStatus /= PGS_S_SUCCESS) THEN 
   CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeBeginningTime")
ENDIF

returnStatus = pgs_met_setAttr_s(groups(INVENTORY), "RangeEndingDate", &
                          L1PCF%endUTC(1:10))
IF (returnStatus /= PGS_S_SUCCESS) THEN 
   CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeEndingDate")
ENDIF

returnStatus = pgs_met_setAttr_s(groups(INVENTORY), "RangeEndingTime", &
                          L1PCF%endUTC(12:26))
IF (returnStatus /= PGS_S_SUCCESS) THEN 
  CALL MLSMessage (MLSMSG_Error, ModuleName, "Error setting RangeEndingTime")
ENDIF


! PGEVersion

attrName = 'PGEVersion'
returnStatus = pgs_met_setAttr_s(groups(INVENTORY), attrName, '0.5')
IF (returnStatus /= PGS_S_SUCCESS) THEN
   errmsg = METAWR_ERR // attrName
   CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
ENDIF



sdid = sfstart(physical_fileName, DFACC_RDWR) 

IF (sdid == -1) THEN
  CALL MLSMessage (MLSMSG_Error, ModuleName, "Failed to open the hdf file "//physical_fileName ) 
ENDIF

returnStatus = pgs_met_write(groups(INVENTORY),&
                            "coremetadata.0", sdid)

IF (returnStatus /= PGS_S_SUCCESS .AND. &
    returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
   IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
     CALL MLSMessage (MLSMSG_WARNING, ModuleName, "Some of the mandatory parameters were not set" )
   ELSE 
     CALL Pgs_smf_getMsg(returnStatus, attrname, errmsg)
     CALL MLSMessage (MLSMSG_WARNING, ModuleName,"Metadata write failed "//attrname//errmsg) 
   ENDIF
ENDIF

hdfReturn = sfend(sdid)

returnStatus=pgs_met_remove() 

END SUBROUTINE populate_metadata_l1

SUBROUTINE WriteMetaData

  CALL populate_metadata_l1(mlspcf_l1b_radf_start, mlspcf_mcf_l1bradf_start)
  CALL populate_metadata_l1(mlspcf_l1b_radd_start, mlspcf_mcf_l1bradd_start)
  CALL populate_metadata_l1(mlspcf_l1b_oa_start, mlspcf_mcf_l1boa_start)

END SUBROUTINE WriteMetadata

END MODULE WriteMetaL1 

! $Log$
! Revision 2.2  2001/03/06 21:03:48  perun
! Fixed typo in 'ReprocessingPlanned' attribute
!
! Revision 2.1  2001/02/23 18:44:57  perun
! Fixed sval length
!
