! -------------------------------------------------------
MODULE WriteMetaL1 ! Populate metadata and write it out
! -------------------------------------------------------

  USE MLSMessageModule, ONLY: MLSMSG_Error, MLSMSG_Warning, MLSMessage
  USE SDPToolkit
  USE MLSPCF1
  USE Hdf, ONLY: DFACC_RDWR, sfstart, sfend
  USE Orbit, ONLY: orbitNumber, numOrb

  IMPLICIT NONE

  PUBLIC

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------


CONTAINS

  SUBROUTINE populate_metadata_l1 (HDF_FILE, MCF_FILE)

    USE InitPCFs, ONLY: L1PCF

    !Arguments

    INTEGER :: HDF_FILE, MCF_FILE

    !Local Variables
 
    INTEGER :: hdfReturn
    INTEGER :: returnStatus
    INTEGER :: sdid

    REAL(r8) dval
    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version, ival, indx
    CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = &
         'Error writing metadata attribute '

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    CHARACTER (LEN = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    INTEGER, EXTERNAL :: pgs_met_init, pgs_met_setattr_d, &
         pgs_met_setAttr_s, pgs_met_getsetattr_d, PGS_MET_SETATTR_I, &
         pgs_met_write, pgs_met_remove

    !Executable code

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    returnStatus = pgs_met_init (MCF_FILE, groups)

    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Initialization error.  See LogStatus for details.") 
    ENDIF

    ! Set PGE values 

    ! ECSDataGranule

    attrName = 'ReprocessingPlanned'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'further update anticipated using enhanced PGE')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'ReprocessingActual'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'processed once')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'LocalGranuleID'
    sval = physical_filename
    indx = INDEX (sval, "/", .TRUE.) + 1  ! Begin after last "/"
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval(indx:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'DayNightFlag'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, 'Both')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'LocalVersionID'
    READ (L1PCF%Cycle, '(I3)') ival
    WRITE (sval, '("C", i2.2)') ival
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! MeasuredParameterContainer

    IF (hdf_file == mlspcf_l1b_radf_start) THEN
       sval = "Filter bank radiances"
    ELSE IF (hdf_file == mlspcf_l1b_radd_start) THEN
       sval = "DACS radiances"
    ELSE IF (hdf_file == mlspcf_l1b_oa_start) THEN
       sval = "Orbit/attitude and tangent point"
    ELSE IF (hdf_file == mlspcf_l1b_eng_start) THEN
       sval = "MLS Instrument Engineering"
    ENDIF
    attrName = 'ParameterName' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! QAFlags Group

    attrName = 'AutomaticQualityFlag' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Passed')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'AutomaticQualityFlagExplanation' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'pending algorithm update')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'OperationalQualityFlag' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'OperationalQualityFlagExplanation' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! QAStats Group
    
    attrName = 'QAPercentInterpolatedData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'QAPercentMissingData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'QAPercentOutofBoundsData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! Orbit Calculated Spatial Domain Container

    attrName = 'OrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'StartOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         orbitNumber(1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'StopOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         orbitNumber(numOrb))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'EquatorCrossingLongitude' // '.1'
    dval = 0.0
    returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    attrName = 'EquatorCrossingTime' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         '00:00:00')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    indx = INDEX (L1PCF%startUTC, "T")
    attrName = 'EquatorCrossingDate' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         L1PCF%startUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! InputPointer
    
    attrName = 'InputPointer'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'See the PCF annotation to this file.')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    IF (MCF_FILE /= mlspcf_mcf_l1boa_start) THEN
       ! Locality Value

       attrName = 'LocalityValue'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Limb')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       ! VerticalSpatialDomain Product-Specific Attribute
       
       attrName = 'VerticalSpatialDomainType' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Atmosphere Layer')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       attrName = 'VerticalSpatialDomainValue' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Brightness Temperature')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       ! HorizontalSpatialDomainContainer

       attrName = 'ZoneIdentifier'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Other Grid System')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       attrName = 'WestBoundingCoordinate'
       dval = -180.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       attrName = 'NorthBoundingCoordinate'
       dval = 90.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       attrName = 'EastBoundingCoordinate'
       dval = 180.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF

       attrName = 'SouthBoundingCoordinate'
       dval = -90.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       ENDIF
    ENDIF

    indx = INDEX (L1PCF%startUTC, "T")

    ! RangeDateTime Group

    attrName = 'RangeBeginningDate'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         attrName, L1PCF%startUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Error setting RangeBeginningDate")
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         "RangeBeginningTime", L1PCF%startUTC(indx+1:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Error setting RangeBeginningTime")
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingDate", &
         L1PCF%endUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Error setting RangeEndingDate")
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingTime", &
         L1PCF%endUTC(indx+1:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Error setting RangeEndingTime")
    ENDIF

    ! PGEVersion
    
    attrName = 'PGEVersion'
    sval = L1PCF%OutputVersion
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    IF (sdid == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Failed to open the hdf file "//physical_fileName ) 
    ENDIF

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    IF (returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
       IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Some of the mandatory parameters were not set" )
       ELSE 
          CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Metadata write failed "//attrname//errmsg) 
       ENDIF
    ENDIF

    hdfReturn = sfend(sdid)

    returnStatus = pgs_met_remove() 

  END SUBROUTINE populate_metadata_l1

  SUBROUTINE WriteMetaData

    CALL populate_metadata_l1 (mlspcf_l1b_radf_start, mlspcf_mcf_l1bradf_start)

    CALL populate_metadata_l1 (mlspcf_l1b_radd_start, mlspcf_mcf_l1bradd_start)

    CALL populate_metadata_l1 (mlspcf_l1b_oa_start, mlspcf_mcf_l1boa_start)

  END SUBROUTINE WriteMetadata

END MODULE WriteMetaL1 

! $Log$
! Revision 2.3  2001/03/22 20:15:48  perun
! Corrected valids
!
! Revision 2.2  2001/03/06 21:03:48  perun
! Fixed typo in 'ReprocessingPlanned' attribute
!
! Revision 2.1  2001/02/23 18:44:57  perun
! Fixed sval length
!
