
! -------------------------------------------------------
MODULE WriteMetadata ! Populate metadata and write it out
! -------------------------------------------------------
USE Hdf 
USE HDFEOS 
use LEXER_CORE, only: PRINT_SOURCE
USE MLSCF 
USE MLSCommon, only: NameLen
USE MLSFiles, only: split_path_name
USE MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
USE MLSPCF2
USE MLSStrings 
USE output_m, only: output
USE PCFHdr, only: WritePCF2Hdr
USE SDPToolkit
use TREE, only: DUMP_TREE_NODE, SOURCE_REF

IMPLICIT NONE

PUBLIC

PRIVATE :: Id, ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

	public :: populate_metadata_std, populate_metadata_oth, &
	& get_l2gp_mcf, WriteMetaLog
	
	private :: first_grouping, measured_parameter, third_grouping, &
	& ExpandFileTemplate

  integer, private :: module_error

! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

   TYPE PCFData_T

     ! cycle # of processing run

     CHARACTER (LEN=4) :: cycle

     ! version string in PCF input file names

     CHARACTER (LEN=15) :: inputVersion	! input files (but which?)

     ! version string in PCF output file names

     CHARACTER (LEN=15) :: PGEVersion	! add to output files

      CHARACTER(LEN=27) :: StartUTC
      CHARACTER(LEN=27) :: EndUTC

	! The correspondence between MCF and l2gp files is determined by
	! the value of        MCFFORL2GPOPTION
	! Either
	!                          (1)
	! The PCF numbers for the mcf corresponding to each
	! of the l2gp files begin with mlspcf_mcf_l2gp_start
	! and increase 1 by 1 with each succeeding species.
	! Then, after the last single-species l2gp, the very next pcf number
	! is for the one called 'other' ML2OTH.001.MCF
	! This hateful inflexibility leads to possibility
	
	!                          (2)
	! Each l2gp file name, stripped of their paths, fits the pattern like
	!  *_l2gp_species_*
	! and the corresponding MCF files fit the pattern
	!  *SPECIES.*
	! where species and SPECIES are case-insensitive "species" name
	! i.e., BrO, ClO, etc.
	! Warning:
	! You therefore must use exactly the same abbreviation for the l2gp and the
	! corresponding MCF: if the MCF is ML2T.001.MCF, don't use "temp"
	! in the l2gp name
	! This inflexibility replaces the different kind in option (1)

	!                          (3)
	! Similar to (2), but now the SPECIES of the corresponding
	! MCF file names are chosen from an associative array structure
	! or hash table, where the keys are the possible species and
	! the hash
     ! how to associate l2gp species names with mcf files

	  ! if an associative array, then l2gp file names contain the "keys"
	  ! the following is a comma-delimited list of possible species
	  ! names that may be such keys

     CHARACTER (LEN=FileNameLen) :: spec_keys

	  ! the following is a comma-delimited list of possible 	
	  ! mcf file name parts that may be the corresponding hash

     CHARACTER (LEN=FileNameLen) :: spec_hash

     ! name of the log file (without the path)

     CHARACTER (LEN=FileNameLen) :: logGranID

   END TYPE PCFData_T

	INTEGER, PARAMETER :: INVENTORYMETADATA=2

    integer, parameter :: MCFFORL2GPOPTION=1		! Either 1 or 2

CONTAINS

!--------------------------- first_grouping -------------------

  SUBROUTINE first_grouping (HDF_FILE, MCF_FILE, l2pcf, groups)

! This writes the metadata for the following attributes:
! (attributes marked automatic are not explicitly written, however)

! SizeMBECSDataGranule (automatic)
! ReProcessingPlanned
! ReProcessingActual
! LocalGranuleID
! DayNightFlag
! ProductionDateTime (automatic)
! LocalVersionID
!

!    USE InitPCFs, ONLY: L2PCF

    !Arguments

    INTEGER :: HDF_FILE, MCF_FILE
	 type(PCFData_T) :: l2pcf

    !Local Variables
 
    INTEGER :: returnStatus

    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version, indx
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

    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in 1st grouping.") 
      call announce_error(0, &
      & "Error in getting ref for PCF number in 1st grouping.") 
    ENDIF

    returnStatus = pgs_met_init (MCF_FILE, groups)

    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Initialization error.  See LogStatus for details.") 
      call announce_error(0, &
      & "Metadata initialization error") 
    ENDIF

    ! Set PGE values 

    ! ECSDataGranule

    attrName = 'ReprocessingPlanned'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'further update anticipated using enhanced PGE')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing ReprocessingPlanned attribute.") 
    ENDIF

    attrName = 'ReprocessingActual'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'processed once')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing ReprocessingActual attribute.") 
    ENDIF

    attrName = 'LocalGranuleID'
    sval = physical_filename
    indx = INDEX (sval, "/", .TRUE.) + 1  ! Begin after last "/"
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval(indx:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing LocalGranuleID attribute.") 
    ENDIF

    attrName = 'DayNightFlag'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, 'Both')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing DayNightFlag attribute.") 
    ENDIF

    attrName = 'LocalVersionID'
         CALL ExpandFileTemplate('$cycle', sval, cycle=l2pcf%cycle)
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing LocalVersionID attribute.") 
    ENDIF

END SUBROUTINE first_grouping

!--------------------------- measured_parameter -------------------

  SUBROUTINE measured_parameter (HDF_FILE, field_name, groups)

! This writes the attributes corresponding to the measured parameter container:
!
! ParameterName
! AutomaticQualityFlag
! AutomaticQualityFlagExplanation
! OperationalQualityFlag
! OperationalQualityFlagExplanation
! ScienceQualityFlag
! ScienceQualityFlagExplanation
! QAPercentInterpolatedData
! QAPercentMissingData
! QAPercentOutOfBoundsData
!

!    USE InitPCFs, ONLY: L2PCF

    !Arguments

    INTEGER :: HDF_FILE
	 character (LEN=*) :: field_name

    !Local Variables
 
    INTEGER :: returnStatus

    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version
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


    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in measured_parameter.") 
      call announce_error(0, &
      & "Error in getting ref for PCF id in measured_parameter.") 
    ENDIF

    ! MeasuredParameterContainer

!    IF (hdf_file == mlspcf_l1b_radf_start) THEN
!       sval = "Filter bank radiances"
!    ELSE IF (hdf_file == mlspcf_l1b_radd_start) THEN
!       sval = "DACS radiances"
!    ELSE IF (hdf_file == mlspcf_l1b_oa_start) THEN
!       sval = "Orbit/attitude and tangent point"
!    ELSE IF (hdf_file == mlspcf_l1b_eng_start) THEN
!       sval = "MLS Instrument Engineering"
!    ENDIF

	if(field_name /= ' ') then
		sval = adjustl(field_name)
	else
		sval = 'Miscellaneous'
	endif
    attrName = 'ParameterName' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
      call announce_error(0, &
      & "Error in writing ParameterName attribute.") 
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
    ENDIF

    ! QAFlags Group

    attrName = 'AutomaticQualityFlag' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Passed')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing AutomaticQualityFlag attribute.") 
    ENDIF

    attrName = 'AutomaticQualityFlagExplanation' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'pending algorithm update')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing AutomaticQualityFlagExplanation attribute.") 
    ENDIF

    attrName = 'OperationalQualityFlag' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OperationalQualityFlag attribute.") 
    ENDIF

    attrName = 'OperationalQualityFlagExplanation' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OperationalQualityFlagExplanation attribute.") 
    ENDIF

    ! QAStats Group
    
    attrName = 'QAPercentInterpolatedData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentInterpolatedData attribute.") 
    ENDIF

    attrName = 'QAPercentMissingData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentMissingData attribute.") 
    ENDIF

    attrName = 'QAPercentOutofBoundsData' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentOutofBoundsData attribute.") 
    ENDIF

END SUBROUTINE measured_parameter

!--------------------------- third_grouping -------------------

  SUBROUTINE third_grouping (HDF_FILE, l2pcf, groups)

! This writes the following metadata attributes:

! OrbitNumber
! StartOrbitNumber
! StopOrbitNumber
! EquatorCrossingLongitude
! EquatorCrossingTime
! EquatorCrossingDate
! ShortName
! VersionID
! InputPointer
! LocalityValue
! VerticalSpatialDomainType
! VerticalSpatialDomainValue
! ZOneIdentifier
! WestBoundingCoordinate
! NortBoundingCoordinate
! EastBoundingCoordinate
! SouthBoundingCoordinate
! RangeBeginningDate
! RangeBeginningTime
! RangeEndingDate
! RangeEndingTime
! PGEVersion
!

!    USE InitPCFs, ONLY: L2PCF

    !Arguments

    INTEGER :: HDF_FILE
	 type(PCFData_T) :: l2pcf

    !Local Variables
 
    INTEGER :: hdfReturn
    INTEGER :: returnStatus
    INTEGER :: sdid

    REAL(r8) dval
    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: sval
    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version, indx
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

    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in third_grouping.") 
      call announce_error(0, &
      & "Error in getting ref for PCF id in third_grouping.") 
    ENDIF

    ! Orbit Calculated Spatial Domain Container

    attrName = 'OrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OrbitNumber attribute.") 
    ENDIF

! Start, Stop orbit numbers: level one has actual calculated numbers
! but, for now at least, we'll not trouble
    attrName = 'StartOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         -1)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing StartOrbitNumber attribute.") 
    ENDIF

    attrName = 'StopOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         -1)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing StopOrbitNumber attribute.") 
    ENDIF

    attrName = 'EquatorCrossingLongitude' // '.1'
    dval = 0.0
    returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingLongitude attribute.") 
    ENDIF

    attrName = 'EquatorCrossingTime' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         '00:00:00')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingTime attribute.") 
    ENDIF

    indx = INDEX (L2PCF%startUTC, "T")
    attrName = 'EquatorCrossingDate' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         L2PCF%startUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingDate attribute.") 
    ENDIF

    ! InputPointer
    
    attrName = 'InputPointer'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'See the PCF annotation to this file.')
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing InputPointer attribute.") 
    ENDIF

       ! Locality Value

       attrName = 'LocalityValue'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Limb')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing LocalityValue attribute.") 
       ENDIF

       ! VerticalSpatialDomain Product-Specific Attribute
       
       attrName = 'VerticalSpatialDomainType' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Atmosphere Layer')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing VerticalSpatialDomainType attribute.") 
       ENDIF

       attrName = 'VerticalSpatialDomainValue' // '.1'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Brightness Temperature')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing VerticalSpatialDomainValue attribute.") 
       ENDIF

       ! HorizontalSpatialDomainContainer

       attrName = 'ZoneIdentifier'
       returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
            'Other Grid System')
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing ZoneIdentifier attribute.") 
       ENDIF

       attrName = 'WestBoundingCoordinate'
       dval = -180.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing WestBoundingCoordinate attribute.") 
       ENDIF

       attrName = 'NorthBoundingCoordinate'
       dval = 90.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing NorthBoundingCoordinate attribute.") 
       ENDIF

       attrName = 'EastBoundingCoordinate'
       dval = 180.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EastBoundingCoordinate attribute.") 
       ENDIF

       attrName = 'SouthBoundingCoordinate'
       dval = -90.0
       returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
       IF (returnStatus /= PGS_S_SUCCESS) THEN
          errmsg = METAWR_ERR // attrName
!          CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing SouthBoundingCoordinate attribute.") 
       ENDIF

    indx = INDEX (L2PCF%startUTC, "T")

    ! RangeDateTime Group

    attrName = 'RangeBeginningDate'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         attrName, L2PCF%startUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeBeginningDate")
      call announce_error(0, &
      & "Error in setting RangeBeginningDate attribute.") 
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         "RangeBeginningTime", L2PCF%startUTC(indx+1:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeBeginningTime")
      call announce_error(0, &
      & "Error in setting RangeBeginningTime attribute.") 
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingDate", &
         L2PCF%endUTC(1:indx-1))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeEndingDate")
      call announce_error(0, &
      & "Error in setting RangeEndingDate attribute.") 
    ENDIF

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingTime", &
         L2PCF%endUTC(indx+1:))
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeEndingTime")
      call announce_error(0, &
      & "Error in setting RangeEndingTime attribute.") 
    ENDIF

    ! PGEVersion
    
    attrName = 'PGEVersion'
    sval = L2PCF%PGEVersion
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in setting PGEVersion attribute.") 
    ENDIF

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    IF (sdid == -1) THEN
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open hdf file: "//physical_fileName) 
    ENDIF

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    IF (returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
       IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
!          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!               "Some of the mandatory parameters were not set" )
      call announce_error(0, &
      & "Error: some of the mandatory parameters not set.") 
       ELSE 
!          CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
!          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!               "Metadata write failed "//attrname//errmsg) 
      call announce_error(0, &
      & "Error: metadata write failed in populate_metadata_std.", &
      & error_number=returnStatus) 
       ENDIF
    ENDIF

    hdfReturn = sfend(sdid)

    returnStatus = pgs_met_remove() 

  END SUBROUTINE third_grouping

!--------------------------- populate_metadata_std -------------------

  SUBROUTINE populate_metadata_std (HDF_FILE, MCF_FILE, &
  & l2pcf, field_name, anText, metadata_error)

! This is the standard way to write meta data
! It should work unchanged for the standard l2gp files (e.g. BrO)
! and, with minor changes, for the l2gp file marked "other"
!
! the l2aux files, also called dgm
! and the dgg files will probably require special treatment
! the log file should be able to use the one stolen from level 3
!

!    USE InitPCFs, ONLY: L2PCF

    !Arguments

    INTEGER :: HDF_FILE, MCF_FILE
	 type(PCFData_T) :: l2pcf
	 character (LEN=*) :: field_name
      CHARACTER (LEN=1), POINTER :: anText(:)
		integer, optional, intent(out) :: metadata_error

    !Local Variables
 
    INTEGER :: hdfReturn
    INTEGER :: returnStatus
    INTEGER :: sdid
    CHARACTER (LEN=132) :: attrname, errmsg

    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
!    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version
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

   module_error = 0
	if(present(metadata_error)) metadata_error=1

    version = 1
	 IF(MCF_FILE > 0) then
	    returnStatus = PGS_PC_GetReference (MCF_FILE, version , physical_filename)
	else
		returnStatus = PGSPC_W_NO_REFERENCE_FOUND
	endif

	if(returnStatus /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for MCF_FILE in populate_metadata_std" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_std.") 
			return
    ENDIF
		
    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

	if(returnStatus /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for HDF_FILE in populate_metadata_std" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_std.") 
			return
    ENDIF
		
	call first_grouping(HDF_FILE, MCF_FILE, l2pcf, groups)
	call measured_parameter (HDF_FILE, field_name, groups)
	call third_grouping (HDF_FILE, l2pcf, groups)

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    IF (sdid == -1) THEN
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open the hdf file in populate_metadata_std: "&
      & //TRIM(physical_fileName)) 
				return
    ENDIF

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    IF (returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
       IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
!          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!               "Some of the mandatory parameters were not set" )
      call announce_error(0, &
      & "Error--some of the mandatory parameters were not set") 
!					return
       ELSE 
          CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Metadata write failed "//attrname//errmsg) 
!      call announce_error(0, &
!      & "Error: metadata write failed in populate_metadata_std.", &
!      & error_number=returnStatus) 
!					return
       ENDIF
    ENDIF

    hdfReturn = sfend(sdid)

! Annotate the file with the PCF

         CALL WritePCF2Hdr(physical_filename, anText)

    returnStatus = pgs_met_remove() 

	if(present(metadata_error)) metadata_error=module_error

  END SUBROUTINE populate_metadata_std

!--------------------------- populate_metadata_oth -------------------

  SUBROUTINE populate_metadata_oth (HDF_FILE, MCF_FILE, l2pcf, &
  & numquantitiesperfile, QuantityNames, anText, metadata_error)

! This is specially to write meta data for heterogeneous files
! It should work unchanged for the 'OTH' l2gp files (e.g. ML2OTH.001.MCF)
! and, with minor changes, for the l2aux files 

!    USE InitPCFs, ONLY: L2PCF

    !Arguments

    INTEGER :: HDF_FILE, MCF_FILE, numquantitiesperfile
	 type(PCFData_T) :: l2pcf
	 character (LEN=*), dimension(:) :: QuantityNames
      CHARACTER (LEN=1), POINTER :: anText(:)
		integer, optional, intent(out) :: metadata_error

    !Local Variables
 
    INTEGER :: hdfReturn
    INTEGER :: returnStatus
    INTEGER :: sdid
    CHARACTER (LEN=132) :: attrname, errmsg

    INTEGER, PARAMETER :: INVENTORY=2, ARCHIVE=1
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: physical_filename
!    CHARACTER (LEN=132) :: attrname, errmsg
    INTEGER :: version, indx
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

   module_error = 0
	if(present(metadata_error)) metadata_error=1

    version = 1
	 IF(MCF_FILE > 0) then
	    returnStatus = PGS_PC_GetReference (MCF_FILE, version , physical_filename)
	else
		returnStatus = PGSPC_W_NO_REFERENCE_FOUND
	endif
	
	if(returnStatus /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for MCF_FILE in populate_metadata_oth" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_oth.") 
			return
    ENDIF

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

	if(returnStatus /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for HDF_FILE in populate_metadata_oth" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_oth.") 
			return
    ENDIF

	call first_grouping(HDF_FILE, MCF_FILE, l2pcf, groups)

	do indx=1, numquantitiesperfile

		call measured_parameter (HDF_FILE, &
		& QuantityNames(indx), groups)

	enddo

	call third_grouping (HDF_FILE, l2pcf, groups)

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    IF (sdid == -1) THEN
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open the hdf file: "//TRIM(physical_fileName))
				return
    ENDIF

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    IF (returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET) THEN 
       IF (returnStatus == PGSMET_W_METADATA_NOT_SET) THEN 
!          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!               "Some of the mandatory parameters were not set" )
      call announce_error(0, &
      & "Error: Some of the mandatory parameters were not set in populate_metadata_oth.") 
!					return
       ELSE 
          CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
          CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
               "Metadata write failed "//attrname//errmsg) 
!      call announce_error(0, &
!      & "Error: metdata write failed in populate_metadata_oth.", &
!      & error_number=returnStatus) 
!					return
       ENDIF
    ENDIF

    hdfReturn = sfend(sdid)

! Annotate the file with the PCF

         CALL WritePCF2Hdr(physical_filename, anText)

    returnStatus = pgs_met_remove() 

	if(present(metadata_error)) metadata_error=module_error

  END SUBROUTINE populate_metadata_oth

!--------------------------- get_l2gp_mcf -------------------

  SUBROUTINE get_l2gp_mcf (file_base, mcf, l2pcf, version)

! metadata configuration file (mcf) PCF number corresponding to l2gp number sdid
!
! Arguments
	character(len=*), intent(in) ::  file_base
	integer, intent(in), optional :: version
	integer, intent(inout) ::        mcf
	 type(PCFData_T) :: l2pcf

! Local
	character (len=PGSd_PC_FILE_PATH_MAX) :: sd_full
	character (len=NameLen) :: sd_path
	character (len=NameLen) :: sd_name

	character (len=PGSd_PC_FILE_PATH_MAX) :: mcf_full
	character (len=NameLen) :: mcf_path
	character (len=NameLen) :: mcf_name
	character (len=NameLen) :: mcf_pattern
	integer :: returnStatus, myVersion, i

	character (len=1), parameter :: COMMA = ','

	! Find species name
	! assume sd_name is "*l2gp_species_"
	! hence enclosed between "_" chars after an l2gp

	character (len=1), parameter :: species_delimiter = '_'
	character (len=4), parameter :: l2gp = 'l2gp'
	
! Begin

   if(MCFFORL2GPOPTION == 1) then
      mcf = mcf+1
      return
   endif

!	if(sdid <= 0) then
	if(len(TRIM(file_base)) <= 0) then
		mcf=0
		return
	endif

	! Get full file name for l2gp file

!		if(present(version)) then
!			myVersion=version
!		else
!			myVersion = 1
!		endif

!	returnStatus = PGS_PC_GetReference(sdid, myVersion , sd_full)
	
!	if (returnStatus /= PGS_S_SUCCESS) then 
!		mcf = 0
!		return
!	endif

	! Get full file name for typical MCF file
	do i=mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end
	
		if(present(version)) then
			myVersion=version
		else
			myVersion = 1
		endif

		returnStatus = PGS_PC_GetReference(i, myVersion , mcf_full)
	
		if (returnStatus == PGS_S_SUCCESS) then 
			exit
		endif
		
	enddo

	if (returnStatus /= PGS_S_SUCCESS) then 
		mcf = 0
		return
	endif

	! Split full_file_names into path+name
	
!	call split_path_name(sd_full, sd_path, sd_name)
	call split_path_name(mcf_full, mcf_path, mcf_name)
	
	sd_full = LowerCase(file_base)
	mcf_name = LowerCase(mcf_name)
	
!	i=index(mcf_name, l2gp)

!	if(i <= 0) then
!		mcf=0
!		return
!	endif
	
	! Get species name assuming e.g. '*l2gp_h2o'
	call split_path_name(sd_full, sd_path, sd_name, species_delimiter)
	
	if(len(trim(sd_name)) <= 0) then
		mcf=0
		return
	endif
	
   if(MCFFORL2GPOPTION == 3) then

   ! get mcfspecies name from associative array
   ! if the species name not found in spec_keys, it will return ','
      call GetStringHashElement(l2pcf%spec_keys, l2pcf%spec_hash, &
      & trim(sd_name), sd_full, .TRUE.)

      if(trim(sd_full) == COMMA) then
         mcf=0
         return
      endif

      sd_name = trim(sd_full)

   endif

	! Now try to find mcf file corresponding to species name
	! assuming, e.g. '*h2o.*'

	do i=mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end
	
		if(present(version)) then
			myVersion=version
		else
			myVersion = 1
		endif

		returnStatus = PGS_PC_GetReference(i, myVersion, mcf_full)
	
		if (returnStatus == PGS_S_SUCCESS) then 
			call split_path_name(mcf_full, mcf_path, mcf_name)
			
			! This gives 'o2h*'
			call split_path_name(Reverse(mcf_name), mcf_path, mcf_pattern, '.')
			mcf_pattern = LowerCase(mcf_pattern)
			
			! So reverse it to make h2o
			mcf_pattern = adjustl(Reverse(mcf_pattern))
			
			! Check that pattern matches species name
         ! Warning--this could give a false matching
         ! if the species name matches more than one mcf_pattern
         ! May wish to check for multiple matches later
			if(index(mcf_pattern, sd_name) > 0) then
				mcf = i
				return
			endif
		endif
		
	enddo
	
	mcf = 0

  END SUBROUTINE get_l2gp_mcf

! Lori's routines

!-------------------------------
   SUBROUTINE WriteMetaLog (pcf, metadata_error)
!-------------------------------

! Brief description of subroutine
! This subroutine writes metadata for the log file to a separate ASCII file.

! Arguments

      TYPE( PCFData_T ), INTENT(IN) :: pcf
		integer, optional, intent(out) :: metadata_error

! Parameters
!     These are PCF numbers that *must* be in the PCF file

      INTEGER, PARAMETER :: ASCII_FILE = 101
      INTEGER, PARAMETER :: THE_LOG_FILE = 100

! Functions

      INTEGER, EXTERNAL :: pgs_met_init, pgs_met_remove, pgs_met_setAttr_d
      INTEGER, EXTERNAL :: pgs_met_setAttr_s, pgs_met_write
    integer, external :: pgs_pc_getconfigdata

! Internal

      CHARACTER (LEN=1) :: nullStr
      CHARACTER (LEN=45) :: sval
      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=FileNameLen) :: physical_filename
      CHARACTER (LEN=480) :: msg, msr
      CHARACTER (LEN=PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER :: result, indx, version
		logical, parameter :: DEBUG=.FALSE.

   ! Begin
   module_error = 0
	if(present(metadata_error)) metadata_error=1

      nullStr = ''

    version = 1
    result = PGS_PC_GetReference (mlspcf_mcf_l2log_start, version , physical_filename)

	if(result /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for the Log mcf in WriteMetaLog" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for Log mcf in WriteMetaLog.") 
			return
    ENDIF

! Check that the PCF file contains the entries for the ASCII file
! and for the log file itself
!    version = 1
!    result = PGS_PC_GetReference (ASCII_FILE, version , physical_filename)
    result = pgs_pc_getconfigdata (ASCII_FILE, physical_filename)

	if(result /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for the ASCII in WriteMetaLog" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for the ASCII in WriteMetaLog.") 
			return
    ENDIF

    version = 1
    result = PGS_PC_GetReference (THE_LOG_FILE, version , physical_filename)

	if(result /= PGS_S_SUCCESS) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!            "Failed to find the PCF reference for the LOG in WriteMetaLog" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for the Log in WriteMetaLog.") 
			return
    ENDIF

! Initialize the MCF file

	if(DEBUG) then
		call output('Initialize the MCF file', &
		&  advance='yes')
	endif

      result = pgs_met_init(mlspcf_mcf_l2log_start, groups)
      IF (result /= PGS_S_SUCCESS) then
!			CALL MLSMessage(MLSMSG_Error, ModuleName, &
!                          'Initialization error.  See LogStatus for details.')
      call announce_error(0, &
      & "metadata initialization error in WriteMetaLog.") 
			return
		endif

! Set PGE values

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalGranuleID", &
                                 pcf%logGranID)

      CALL ExpandFileTemplate('$cycle', sval, cycle=pcf%cycle)
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalVersionID", &
                                 sval)

    indx = INDEX (PCF%startUTC, "T")

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningDate", PCF%startUTC(1:indx-1))
!      sval= '00:00:00.000000'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeBeginningTime", PCF%startUTC(indx+1:))
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingDate", PCF%endUTC(1:indx-1))
!      sval= '23:59:59.999999'
      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                                 "RangeEndingTime", PCF%endUTC(indx+1:))

      result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "PGEVersion", &
                                 pcf%PGEVersion)

      IF (result /= PGS_S_SUCCESS) THEN
!         call Pgs_smf_getMsg(result, mnemonic, msg)
!         msr = mnemonic // ':  ' // msg
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      call announce_error(0, &
      & "Error in setting PGEVersion attribute in WriteMetaLog.", &
      &  error_number=result) 
			return
      ENDIF

! Write the metadata and their values to an ASCII file

	if(DEBUG) then
		call output('Write the metadata and their values to an ASCII file', &
		&  advance='yes')
	endif

      result = pgs_met_write(groups(1), nullStr, ASCII_FILE)

      IF (result /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(result, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
!      call announce_error(0, &
!      & "Error: failed to write metadata in WriteMetaLog.", &
!      & error_number=result) 
!			return
      ENDIF

      result = pgs_met_remove()

	if(present(metadata_error)) metadata_error=module_error

!-----------------------------
   END SUBROUTINE WriteMetaLog
!-----------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE ExpandFileTemplate (template, filename, level, version, cycle, &
                                  day)
!-----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine expands the version, cycle, and day fields in a file name
! template.  Note: CYCLE is a STRING here, which is how it's read from the PCF.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: template

      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: cycle, day, level, version

      CHARACTER (LEN=*), INTENT(OUT) :: fileName

! Parameters

! Functions

! Variables

      CHARACTER (LEN=4) :: field, zCy

      INTEGER :: i, iCy, indx, numFields

! Initializations

      numFields = 4

      fileName = template

! Loop through the expandable fields in the template

      DO i = 1, numFields

! Search for $

         indx = INDEX(fileName,'$')

! Exit, if there are no expandable fields 

         IF (indx == 0) EXIT

! Match the field name to the input argument

         field = fileName(indx:indx+3)
    
         IF (field == '$lev') THEN

            IF ( PRESENT(level) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(level) // &
                          fileName((indx+6):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input level &
                                            &required to expand the template.')
            ENDIF

         ELSE IF (field == '$ver') THEN

            IF ( PRESENT(version) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(version) // &
                          fileName((indx+8):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input version &
                                            &required to expand the template.')
            ENDIF

         ELSE IF (field == '$cyc') THEN

            IF ( PRESENT(cycle) ) THEN

! Convert from CHARACTER to INTEGER

               READ(cycle, '(I2)') iCy

! Add a leading zero, if less than 10

               IF (iCy < 10) THEN
                  zCy = '0' // TRIM(cycle)
               ELSE 
                  zCy = cycle
               ENDIF

               fileName = fileName(:(indx-1)) // 'C' // TRIM(zCy) // &
                          fileName((indx+6):)

            ELSE

               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input cycle &
                                            &required to expand the template.')

            ENDIF

         ELSE IF (field == '$day') THEN

            IF ( PRESENT(day) ) THEN
               fileName = fileName(:(indx-1)) // TRIM(day) // &
                          fileName((indx+4):)
            ELSE
               CALL MLSMessage(MLSMSG_Error, ModuleName, 'Input day required &
                                                    &to expand the template.')
            ENDIF

         ENDIF

      ENDDO

!-----------------------------------
   END SUBROUTINE ExpandFileTemplate
!-----------------------------------

  ! ------------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( lcf_where, full_message, use_toolkit, &
    & error_number )
  
    ! Arguments

    integer, intent(in) :: Lcf_where
    character(LEN=*), intent(in) :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional :: Error_number

    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.
    CHARACTER (LEN=132) :: attrname, errmsg
    integer :: Toolbox_error_num

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = use_toolkit

    if ( .not. just_print_it ) then
      module_error = max(module_error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
        call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ": The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error:", advance='yes', &
       & from_where=ModuleName)
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if ( present(error_number) ) then
        call output ( 'Error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
        Toolbox_error_num = error_number
        CALL Pgs_smf_getMsg (Toolbox_error_num, attrname, errmsg)
        CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
            &  TRIM(attrname) //' ' // TRIM(errmsg) )
      end if
    else
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
        call output ( 'Error number ' )
        call output ( error_number, advance='yes' )
      end if
    end if
    
!===========================
  end subroutine Announce_Error
!===========================

END MODULE WriteMetadata 
! $Log$
! Revision 2.8  2001/04/13 00:28:03  pwagner
! Turned get_l2gp_mcf into a subroutine
!
! Revision 2.7  2001/04/12 00:22:17  pwagner
! Added announce_error
!
! Revision 2.6  2001/04/11 20:20:37  pwagner
! Checks correctly on PCF file before attempting LOG.met
!
! Revision 2.5  2001/04/10 23:03:57  pwagner
! Finally seems to work
!
! Revision 2.4  2001/04/09 23:45:47  pwagner
! Deleted unused old populate_metadata; some fixes
!
! Revision 2.3  2001/04/03 23:51:28  pwagner
! Many changes; some may be right
!
! Revision 2.2  2001/04/02 23:42:18  pwagner
! Added populate_metadata_oth
!
! Revision 2.1  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.0  2000/09/05 18:57:07  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.4  2000/06/30 00:16:41  lungu
! Made dval REAL(r8).
!
