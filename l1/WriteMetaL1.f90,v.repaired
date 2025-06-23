! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! -------------------------------------------------------
MODULE WriteMetaL1 ! Populate metadata and write it out
! -------------------------------------------------------

  USE Hdf, ONLY: DFACC_RDWR
  USE Intrinsic, ONLY: L_HDF
  USE MLSCommon, ONLY: R8, L2METADATA_T, NameLen
  USE MLSMessageModule, ONLY: MLSMSG_Error, MLSMSG_Info, MLSMSG_Warning, &
    & MLSMessage
  USE PCFHdr, ONLY: GlobalAttributes, WriteInputPointer, h5_writeglobalattr
  USE SDPToolkit, only: PGSD_PC_FILE_PATH_MAX, PGSD_MET_GROUP_NAME_L, &
    & PGSD_MET_NUM_OF_GROUPS, PGS_S_SUCCESS, PGSMET_W_METADATA_NOT_SET, &
    & PGS_PC_GETREFERENCE
  USE MLSPCF1
  USE Orbit, ONLY: orbitNumber, numOrb

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: WriteMetadata
  integer, parameter :: NUMDOIs    = 6

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

CONTAINS

  SUBROUTINE populate_metadata_l1 (HDF_FILE, MCF_FILE)

    USE InitPCFs, ONLY: L1PCF
    USE MLSFiles, ONLY: mls_sfstart, mls_sfend
    USE MLSL1Common, ONLY: HDFversion

    !Arguments

    INTEGER :: HDF_FILE, MCF_FILE

    !Local Variables
 
    INTEGER :: returnStatus
    INTEGER :: sdid

    character(len=NameLen), dimension(NumDOIs) :: doiArray
    REAL(r8) :: dval
    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: attrname, errmsg
    INTEGER :: version, ival, indx, i
    CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = &
         'Error writing metadata attribute '
    type ( L2METADATA_T ) :: L2METADATA
    integer, parameter :: OADOIIndx   = 1
    integer, parameter :: radDDOIIndx = OADOIIndx + 1
    integer, parameter :: radGDOIIndx = radDDOIIndx + 1
    integer, parameter :: radTDOIIndx = radGDOIIndx + 1
    integer, parameter :: engDOIIndx  = radTDOIIndx + 1

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
    doiArray(1) = '10.5067/AURA/MLS/DATA1001'
    doiArray(2) = '10.5067/AURA/MLS/DATA1002'
    doiArray(3) = '10.5067/AURA/MLS/DATA1003'
    doiArray(4) = '10.5067/AURA/MLS/DATA1004'
    doiArray(5) = '10.5067/AURA/MLS/DATA1005'
    doiArray(6) = '10.5067/AURA/MLS/DATA1006'
    ! This hackery-quackery allows us to use the PCF to
    ! input elements of a string array without using up
    ! a bunch of PCFids
    ! So we accept that the strings are file names
    ! and use Pgs_pc_getReference with the version mechanism
    do i=1, NumDOIs
      version = i
      returnStatus = Pgs_pc_getReference( mlspcf_l1_param_doinames, &
        & version, sval )
      if ( returnStatus /= PGS_S_SUCCESS ) exit
      doiArray(NumDOIs-i+1) = sval
    enddo

    version = 1

    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    returnStatus = pgs_met_init (MCF_FILE, groups)

    IF (returnStatus /= PGS_S_SUCCESS) THEN 
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Initialization error.  See LogStatus for details.") 
    ELSE
      CALL MLSMessage (MLSMSG_Info, ModuleName, &
            "Beginning metadata write to " // trim(physical_filename)) 
    ENDIF

    ! Set PGE values 

    ! ECSDataGranule

    attrName = 'ReprocessingPlanned'

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'further update anticipated using enhanced PGE')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'LocalGranuleID'
    sval = physical_filename
    indx = INDEX (sval, "/", .TRUE.) + 1  ! Begin after last "/"

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval(indx:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'DayNightFlag'

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, 'Both')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'LocalVersionID'
    READ (L1PCF%Cycle, '(I3)') ival
    WRITE (sval, '("c", i2.2)') ival
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    ! MeasuredParameterContainer

    IF (hdf_file == mlspcf_l1b_radf_start) THEN
       sval = "Filter bank radiances"
       l2metaData%doiIdentifier = doiArray(radGDOIIndx)
    ELSE IF (hdf_file == mlspcf_l1b_radd_start) THEN
       sval = "DACS radiances"
       l2metaData%doiIdentifier = doiArray(radDDOIIndx)
    ELSE IF (hdf_file == mlspcf_l1b_radt_start) THEN
       sval = "THz radiances"
       l2metaData%doiIdentifier = doiArray(radTDOIIndx)
    ELSE IF (hdf_file == mlspcf_l1b_oa_start) THEN
       sval = "Orbit/attitude and tangent point"
       l2metaData%doiIdentifier = doiArray(OADOIIndx)
    ELSE IF (hdf_file == mlspcf_l1b_eng_start) THEN
       sval = "MLS Instrument Engineering"
       l2metaData%doiIdentifier = doiArray(engDOIIndx)
    ENDIF
    attrName = 'ParameterName' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    ! QAStats Group
    
    attrName = 'QAPercentInterpolatedData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'QAPercentMissingData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'QAPercentOutofBoundsData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    ! Orbit Calculated Spatial Domain Container

    ! This changes confirm to James Johnson suggestion on 6/12/03

    !attrName = 'OrbitNumber' // '.1'
    !returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    !IF (returnStatus /= PGS_S_SUCCESS) THEN
    !   errmsg = METAWR_ERR // attrName
    !   CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    !ENDIF

    attrName = 'StartOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         orbitNumber(1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'StopOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         orbitNumber(numOrb))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'EquatorCrossingLongitude' // '.1'
    dval = 0.0
    returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    attrName = 'EquatorCrossingTime' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         '00:00:00')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    indx = INDEX (L1PCF%startUTC, "T")
    attrName = 'EquatorCrossingDate' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         L1PCF%startUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    ! InputPointer
    
    attrName = 'InputPointer'
    ! returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
    !     'See the PCF annotation to this file.')
    returnStatus = WriteInputPointer(groups(INVENTORY), attrName, &
      & fileType=l_hdf)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    IF (MCF_FILE /= mlspcf_mcf_l1boa_start) THEN
       ! Locality Value

       attrName = 'LocalityValue'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Limb')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       ! VerticalSpatialDomain Product-Specific Attribute
       
       attrName = 'VerticalSpatialDomainType' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Atmosphere Layer')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       !This changes confirm to James Johnson suggestion on 6/12/03
       attrName = 'VerticalSpatialDomainValue' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Brightness Temperature')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       ! HorizontalSpatialDomainContainer

       attrName = 'ZoneIdentifier'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Other Grid System')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       attrName = 'WestBoundingCoordinate'
       dval = -180.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       attrName = 'NorthBoundingCoordinate'
       dval = 90.0

       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       attrName = 'EastBoundingCoordinate'
       dval = 180.0

       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
       ENDIF

       attrName = 'SouthBoundingCoordinate'
       dval = -90.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
          CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
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
       CALL MLSMessage (MLSMSG_Error, ModuleName, TRIM(errmsg))
    ENDIF

    ! Production Location

    attrName = 'ProductionLocation'
    sval = 'MLS_SIPS' ! GlobalAttributes%productionLoc
    GlobalAttributes%productionLoc = sval
    returnStatus = pgs_met_setAttr_s ( groups(INVENTORY), attrName, sval )
    if ( returnStatus /= PGS_S_SUCCESS ) then
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Warning, ModuleName, TRIM(errmsg))
    end if

    attrName = 'DataProducer'
    sval = 'MLS_SIPS'
    returnStatus = pgs_met_setAttr_s ( groups(INVENTORY), attrName, sval )
    if ( returnStatus /= PGS_S_SUCCESS ) then
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Warning, ModuleName, TRIM(errmsg))
    end if

    ! DOI
    ! This should be in the MCF
    attrName = 'identifier_product_DOI'
    sval = l2metaData%doiIdentifier
    GlobalAttributes%DOI = sval
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
       errmsg = METAWR_ERR // attrName
       CALL MLSMessage(MLSMSG_Warning, ModuleName, TRIM(errmsg))
    end if

    sdid = mls_sfstart (physical_fileName, DFACC_RDWR, hdfVersion, .TRUE.)

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    IF (returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
       IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Some of the mandatory parameters were not set" )
       ELSE 
          CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Metadata write failed "//attrname//TRIM(errmsg)//" for file "& 
               //physical_fileName ) 
       ENDIF
    ENDIF

    returnStatus = mls_sfend (sdid, hdfVersion, .TRUE.)
    IF (returnStatus /= 0) THEN 
       CALL MLSMessage (MLSMSG_ERROR, ModuleName, &
            "Calling mls_sfend failed for file "//physical_fileName ) 
    ENDIF          

    returnStatus = pgs_met_remove()
    if ( returnStatus /= PGS_S_SUCCESS ) sdid = sdid + 1 ! For no reason

    ! Write global attributes

    sdid = mls_sfstart (physical_fileName, DFACC_RDWR, hdfVersion, .FALSE.)
    call h5_writeglobalattr( sdid, skip_if_already_there=.false., doi=.true. )
    returnStatus = mls_sfend (sdid, hdfVersion, .FALSE.)
    IF (returnStatus /= 0) THEN 
       CALL MLSMessage (MLSMSG_ERROR, ModuleName, &
            "Calling mls_sfend failed for file "//physical_fileName ) 
    ENDIF          

  END SUBROUTINE populate_metadata_l1

  SUBROUTINE WriteMetaData (IsTHz)

    USE MLSL1Config, ONLY: L1Config

    LOGICAL, OPTIONAL :: IsTHz

    IF (PRESENT (IsTHz)) THEN

       CALL populate_metadata_l1 (mlspcf_l1b_radt_start, &
            mlspcf_mcf_l1bradt_start)

       RETURN

    ENDIF

    CALL populate_metadata_l1 (mlspcf_l1b_radf_start, mlspcf_mcf_l1bradf_start)

    IF (L1Config%Calib%CalibDACS) &
         CALL populate_metadata_l1 (mlspcf_l1b_radd_start, &
         mlspcf_mcf_l1bradd_start)

    CALL populate_metadata_l1 (mlspcf_l1b_oa_start, mlspcf_mcf_l1boa_start)

  END SUBROUTINE WriteMetadata

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE WriteMetaL1 

! $Log$
! Revision 2.22  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.21.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.21  2014/04/15 23:02:17  pwagner
! Corrected DOIs based on ESDIS info
!
! Revision 2.20  2014/04/11 16:51:46  pwagner
! Added ProductionLocation, DataProducer, DOI metadata
!
! Revision 2.19  2007/06/21 21:06:20  perun
! Only output to RADD file if DACS calibration is enabled
!
! Revision 2.18  2006/04/05 18:09:52  perun
! Remove unused variables
!
! Revision 2.17  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.16  2005/01/27 00:35:06  pwagner
! ReprocessingActual field dropped from product metadata
!
! Revision 2.15  2004/12/16 15:09:13  cvuu
! v1.5: Change value of ReprocessingActual to unknown, remove QAFlags Group and use the ones in MCF v1.5
!
! Revision 2.14  2004/01/30 00:30:42  pwagner
! Stops useless warnings about pgs_met_remove return value
!
! Revision 2.13  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.12  2003/08/12 16:57:25  cvuu
! brought closer to James Johnson want to
!
! Revision 2.11  2003/07/08 00:17:11  pwagner
! fileType now a lit_name instead of a char string
!
! Revision 2.10  2003/06/03 20:44:42  pwagner
! Writes global attributes
!
! Revision 2.9  2003/05/30 23:48:34  pwagner
! Uses WriteInputPointer from lib/PCFHdr
!
! Revision 2.8  2003/03/15 00:15:10  pwagner
! Wont quit if pgs_met_remove returns non-zero value
!
! Revision 2.7  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.6  2002/11/19 21:46:38  perun
! Use HDFversion instead of HDFVersionString
!
! Revision 2.5  2002/11/07 21:32:35  jdone
! Added HDF4/HDF5 capabilities.
!
! Revision 2.4  2002/01/09 23:53:09  pwagner
! Now gets r8 explicitly from MLSCommon
!
! Revision 2.3  2001/03/22 20:15:48  perun
! Corrected valids
!
! Revision 2.2  2001/03/06 21:03:48  perun
! Fixed typo in 'ReprocessingPlanned' attribute
!
! Revision 2.1  2001/02/23 18:44:57  perun
! Fixed sval length
!
