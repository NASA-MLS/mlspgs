! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! -------------------------------------------------------
module WriteMetadata ! Populate metadata and write it out
! -------------------------------------------------------

  use Hdf, only: DFACC_RDWR, Sfend, Sfstart
! use HDFEOS ! Appears not to be used
  use LEXER_CORE, only: PRINT_SOURCE
! use MLSCF  ! Appears not to be used
  use MLSCommon, only: NameLen, R8
  use MLSFiles, only: Split_path_name
  use MLSL2Options, only: PENALTY_FOR_NO_METADATA
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSPCF2, only: Mlspcf_mcf_l2gp_end, Mlspcf_mcf_l2gp_start, &
    & Mlspcf_mcf_l2log_start
  use MLSStrings, only: Reverse, LowerCase, GetStringHashElement
  use Output_m, only: Output
  use PCFHdr, only: WritePCF2Hdr
  use SDPToolkit, only: FileNameLen, PGSd_MET_GROUP_NAME_L, &
    & PGSd_MET_NUM_OF_GROUPS, PGSd_PC_FILE_PATH_MAX, PGS_PC_GetReference, &
    & PGSPC_W_NO_REFERENCE_FOUND, PGS_S_SUCCESS, PGSMET_W_METADATA_NOT_SET
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  implicit none

  private

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
     "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  public :: Populate_metadata_std, Populate_metadata_oth, &
    & Get_l2gp_mcf, WriteMetaLog

  private :: First_grouping, Measured_parameter, Third_grouping, &
    & ExpandFileTemplate

! This data type is used to store User-defined Runtime Parameters and other
! information taken from the PCF.

  type, public :: PCFData_T

    ! cycle # of processing run

    character (len=4) :: cycle

    ! version string in PCF input file names

    character (len=15) :: InputVersion ! input files (but which?)

    ! version string in PCF output file names

    character (len=15) :: PGEVersion   ! add to output files

    character(len=27) :: StartUTC
    character(len=27) :: EndUTC

    ! The annotation text to be written to the header of every scientific
    ! data file for which we need metadata
    ! In practice, this annotation will be the contents of the PCF file itself

    character (len=1), pointer :: anText(:) => null()

    ! The correspondence between MCF and l2gp files is determined by
    ! the value of        MCFFORL2GPOPTION
    ! One of three possible options:
    !                          (1)
    ! The PCF numbers for the mcf corresponding to each
    ! of the l2gp files begin with mlspcf_mcf_l2gp_start
    ! and increase 1 by 1 with each succeeding species.
    ! Then, after the last single-species l2gp, the very next pcf number
    ! is for the one called 'other' ML2OTH.001.MCF
    ! This inconvenient inflexibility is relieved in option (2) or (3)

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

    ! Note that the text matching in options (2) and (3)
    ! will not be case sensitive unless you set MCFCASESENSITIVE

    character (len=fileNameLen) :: spec_keys

    ! the following is a comma-delimited list of possible 	
    ! mcf file name parts that may be the corresponding hash

    character (len=fileNameLen) :: spec_hash

    ! name of the log file (without the path)

    character (len=fileNameLen) :: logGranID

  end type PCFData_T

  integer, public, parameter :: MCFFORL2GPOPTION=3     ! 1, public, 2 or 3

  logical, public, parameter :: MCFCASESENSITIVE=.FALSE.

  integer, public, parameter :: INVENTORYMETADATA=2

  integer, private :: Module_error

contains

  ! ---------------------------------------------  First_grouping  -----

  subroutine First_grouping ( HDF_FILE, MCF_FILE, L2pcf, Groups )

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

    integer :: HDF_FILE, MCF_FILE
    type(PCFData_T) :: l2pcf

    !Local Variables

    integer :: ReturnStatus

    integer, parameter :: INVENTORY=2, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: sval
    character (len=132) :: attrname, errmsg
    integer :: version, indx
    character (len=*), parameter :: METAWR_ERR = &
         'Error writing metadata attribute '

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    integer, external :: PGS_MET_init, PGS_MET_setattr_d, &
         PGS_MET_setAttr_s, PGS_MET_getsetattr_d, PGS_MET_setattr_i, &
         PGS_MET_write, PGS_MET_remove

    !Executable code

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in 1st grouping.") 
      call announce_error(0, &
      & "Error in getting ref for PCF number in 1st grouping.") 
    end if

    returnStatus = pgs_met_init (MCF_FILE, groups)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Initialization error.  See LogStatus for details.") 
      call announce_error(0, &
      & "Metadata initialization error") 
    end if

    ! Set PGE values 

    ! ECSDataGranule

    attrName = 'ReprocessingPlanned'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'further update anticipated using enhanced PGE')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing ReprocessingPlanned attribute.") 
    end if

    attrName = 'ReprocessingActual'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'processed once')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing ReprocessingActual attribute.") 
    end if

    attrName = 'LocalGranuleID'
    sval = physical_filename
    indx = INDEX (sval, "/", .TRUE.) + 1  ! Begin after last "/"
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval(indx:))
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing LocalGranuleID attribute.") 
    end if

    attrName = 'DayNightFlag'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, 'Both')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing DayNightFlag attribute.") 
    end if

    attrName = 'LocalVersionID'
         CALL ExpandFileTemplate('$cycle', sval, cycle=l2pcf%cycle)
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing LocalVersionID attribute.") 
    end if

  end subroutine first_grouping

  ! -----------------------------------------  Measured_parameter  -----

  subroutine Measured_parameter ( HDF_FILE, Field_name, Groups, Class_num )

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

    integer :: HDF_FILE
    character(len=*) :: Field_name
    integer, intent(in) :: Class_num

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    !Local Variables

    integer :: returnStatus
    character (len=2) :: class

    integer, parameter :: INVENTORY=2, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: sval
    character (len=132) :: Attrname, Errmsg
    integer :: Version
    character (len=*), parameter :: METAWR_ERR = &
      &  'Error writing metadata attribute '

    ! Externals

    integer, external :: PGS_MET_init, PGS_MET_setattr_d, &
         PGS_MET_setAttr_s, PGS_MET_getsetattr_d, PGS_MET_SETATTR_I, &
         PGS_MET_write, PGS_MET_remove

    !Executable code

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)


    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in measured_parameter.") 
      call announce_error(0, &
      & "Error in getting ref for PCF id in measured_parameter.") 
    end if

    ! MeasuredParameterContainer

!    if ( hdf_file == mlspcf_l1b_radf_start ) then
!       sval = "Filter bank radiances"
!    else if ( hdf_file == mlspcf_l1b_radd_start ) then
!       sval = "DACS radiances"
!    else if ( hdf_file == mlspcf_l1b_oa_start ) then
!       sval = "Orbit/attitude and tangent point"
!    else if ( hdf_file == mlspcf_l1b_eng_start ) then
!       sval = "MLS Instrument Engineering"
!    end if

    if ( class_num <= 0 ) then
      class(:2)='0 '
    else if ( class_num < 10 ) then
      write(class(:1), '(I1)') class_num
      class(2:2)=' '
    else
      write(class(:2), '(I2)') class_num
    end if

    if ( field_name /= ' ' ) then
      sval = adjustl(field_name)
    else
      sval = 'Miscellaneous'
    end if
    attrName = 'ParameterName' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
      call announce_error(0, &
      & "Error in writing ParameterName attribute.") 
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
    end if

    ! QAFlags Group

    attrName = 'AutomaticQualityFlag' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Passed')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing AutomaticQualityFlag attribute.") 
    end if

    attrName = 'AutomaticQualityFlagExplanation' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'pending algorithm update')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing AutomaticQualityFlagExplanation attribute.") 
    end if

    attrName = 'OperationalQualityFlag' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OperationalQualityFlag attribute.") 
    end if

    attrName = 'OperationalQualityFlagExplanation' // '.' // class
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Not Investigated')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage(MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OperationalQualityFlagExplanation attribute.") 
    end if

    ! QAStats Group

    attrName = 'QAPercentInterpolatedData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentInterpolatedData attribute.") 
    end if

    attrName = 'QAPercentMissingData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentMissingData attribute.") 
    end if

    attrName = 'QAPercentOutofBoundsData' // '.' // class
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, 0)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing QAPercentOutofBoundsData attribute.") 
    end if

  end subroutine Measured_parameter

  ! ---------------------------------------------  Third_grouping  -----

  subroutine Third_grouping ( HDF_FILE, L2pcf, Groups )

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

    ! Arguments

    integer :: HDF_FILE
    type(PCFData_T) :: l2pcf

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    !Local Variables

    integer :: HdfReturn
    integer :: ReturnStatus
    integer :: Sdid

    real(r8) Dval
    integer, parameter :: INVENTORY=2, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
    character (len=PGSd_PC_FILE_PATH_MAX) :: sval
    character (len=132) :: attrname, errmsg
    integer :: version, indx
    character (len=*), parameter :: METAWR_ERR = &
      & 'Error writing metadata attribute '


    ! Externals

    integer, external :: PGS_MET_init, PGS_MET_setattr_d, &
         PGS_MET_setAttr_s, PGS_MET_getsetattr_d, PGS_MET_SETATTR_I, &
         PGS_MET_write, PGS_MET_remove

    !Executable code

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error in getting ref for PCF number in third_grouping.") 
      call announce_error(0, &
      & "Error in getting ref for PCF id in third_grouping.") 
    end if

    ! Orbit Calculated Spatial Domain Container

    attrName = 'OrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, -1)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing OrbitNumber attribute.") 
    end if

    ! Start, Stop orbit numbers: level one has actual calculated numbers
    ! but, for now at least, we'll not trouble
    attrName = 'StartOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         -1)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing StartOrbitNumber attribute.") 
    end if

    attrName = 'StopOrbitNumber' // '.1'
    returnStatus = pgs_met_setAttr_i (groups(INVENTORY), attrName, &
         -1)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing StopOrbitNumber attribute.") 
    end if

    attrName = 'EquatorCrossingLongitude' // '.1'
    dval = 0.0
    returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingLongitude attribute.") 
    end if

    attrName = 'EquatorCrossingTime' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         '00:00:00')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingTime attribute.") 
    end if

    indx = INDEX (L2PCF%startUTC, "T")
    attrName = 'EquatorCrossingDate' // '.1'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         L2PCF%startUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing EquatorCrossingDate attribute.") 
    end if

    ! InputPointer

    attrName = 'InputPointer'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'See the PCF annotation to this file.')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing InputPointer attribute.") 
    end if

    ! Locality Value

    attrName = 'LocalityValue'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
         'Limb')
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing LocalityValue attribute.")
     end if

     ! VerticalSpatialDomain Product-Specific Attribute

     attrName = 'VerticalSpatialDomainType' // '.1'
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Atmosphere Layer')
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing VerticalSpatialDomainType attribute.")
     end if

     attrName = 'VerticalSpatialDomainValue' // '.1'
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Brightness Temperature')
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing VerticalSpatialDomainValue attribute.")
     end if

     ! HorizontalSpatialDomainContainer

     attrName = 'ZoneIdentifier'
     returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, &
          'Other Grid System')
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing ZoneIdentifier attribute.")
     end if

     attrName = 'WestBoundingCoordinate'
     dval = -180.0
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing WestBoundingCoordinate attribute.")
     end if

     attrName = 'NorthBoundingCoordinate'
     dval = 90.0
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing NorthBoundingCoordinate attribute.")
     end if

     attrName = 'EastBoundingCoordinate'
     dval = 180.0
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
       call announce_error(0, &
       & "Error in writing EastBoundingCoordinate attribute.")
     end if

     attrName = 'SouthBoundingCoordinate'
     dval = -90.0
     returnStatus = pgs_met_setAttr_d (groups(INVENTORY), attrName, dval)
     if ( returnStatus /= PGS_S_SUCCESS ) then
     !   errmsg = METAWR_ERR // attrName
!        CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in writing SouthBoundingCoordinate attribute.")
    end if

    indx = INDEX (L2PCF%startUTC, "T")

    ! RangeDateTime Group

    attrName = 'RangeBeginningDate'
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         attrName, L2PCF%startUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeBeginningDate")
      call announce_error(0, &
      & "Error in setting RangeBeginningDate attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), &
         "RangeBeginningTime", L2PCF%startUTC(indx+1:))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeBeginningTime")
      call announce_error(0, &
      & "Error in setting RangeBeginningTime attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingDate", &
         L2PCF%endUTC(1:indx-1))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeEndingDate")
      call announce_error(0, &
      & "Error in setting RangeEndingDate attribute.") 
    end if

    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), "RangeEndingTime", &
         L2PCF%endUTC(indx+1:))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Error setting RangeEndingTime")
      call announce_error(0, &
      & "Error in setting RangeEndingTime attribute.") 
    end if

    ! PGEVersion

    attrName = 'PGEVersion'
    sval = L2PCF%PGEVersion
    returnStatus = pgs_met_setAttr_s (groups(INVENTORY), attrName, sval)
    if ( returnStatus /= PGS_S_SUCCESS ) then
    !   errmsg = METAWR_ERR // attrName
!       CALL MLSMessage (MLSMSG_Error, ModuleName, errmsg)
      call announce_error(0, &
      & "Error in setting PGEVersion attribute.") 
    end if

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    if ( sdid == -1 ) then
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open hdf file: "//physical_fileName) 
    end if

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    if ( returnStatus /= PGS_S_SUCCESS .AND. &
      &  returnStatus /= PGSMET_W_METADATA_NOT_SET ) then 
      if ( returnStatus == PGSMET_W_METADATA_NOT_SET ) then 
!       CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!        & "Some of the mandatory parameters were not set" )
        call announce_error(0, &
        & "Error: some of the mandatory parameters not set.") 
      else
!       CALL Pgs_smf_getMsg (returnStatus, attrname, errmsg)
!       CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!         &  "Metadata write failed "//attrname//errmsg)
        call announce_error(0, &
        & "Error: metadata write failed in third_grouping.", &
        & error_number=returnStatus) 
      end if
    end if

    hdfReturn = sfend(sdid)

    returnStatus = pgs_met_remove() 

  end subroutine Third_grouping

  ! --------------------------------------  Populate_metadata_std  -----

  subroutine Populate_metadata_std ( HDF_FILE, MCF_FILE, &
    & L2pcf, Field_name, Metadata_error )

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

    integer :: HDF_FILE, MCF_FILE
    type(PCFData_T) :: L2pcf
    character (len=*) :: Field_name
    integer, optional, intent(out) :: Metadata_error

    !Local Variables

    integer :: hdfReturn
    integer :: returnStatus
    integer :: sdid
    character (len=132) :: attrname, errmsg

    integer, parameter :: INVENTORY=2, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
!    character (len=132) :: attrname, errmsg
    integer :: version
    character (len=*), parameter :: METAWR_ERR = &
      &  'Error writing metadata attribute '

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    integer, external :: PGS_MET_init, PGS_MET_setattr_d, &
      &  PGS_MET_setAttr_s, PGS_MET_getsetattr_d, PGS_MET_SETATTR_I, &
      &  PGS_MET_write, PGS_MET_remove

    !Executable code

    module_error = 0
    if ( present(metadata_error)) metadata_error=1

    version = 1
    if ( mcf_file > 0 ) then
      returnStatus = PGS_PC_GetReference (MCF_FILE, version , physical_filename)
    else
      returnStatus = PGSPC_W_NO_REFERENCE_FOUND
    end if

    if ( returnStatus /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for MCF_FILE in populate_metadata_std" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_std.") 
      return
    end if
		
    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for HDF_FILE in populate_metadata_std" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_std.") 
      return
    end if
		
    call first_grouping(HDF_FILE, MCF_FILE, l2pcf, groups)
    call measured_parameter (HDF_FILE, field_name, groups, 1)
    call third_grouping (HDF_FILE, l2pcf, groups)

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    if ( sdid == -1 ) then
!       CALL MLSMessage (MLSMSG_Error, ModuleName, &
!            "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open the hdf file in populate_metadata_std: "&
      & //TRIM(physical_fileName)) 
      return
    end if

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    if ( returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET ) then 
      if ( returnStatus == PGSMET_W_METADATA_NOT_SET ) then 
!       CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!         & "Some of the mandatory parameters were not set" )
        call announce_error(0, &
        & "Error--some of the mandatory parameters were not set") 
!					return
      else
        call Pgs_smf_getMsg (returnStatus, attrname, errmsg)
        call MLSMessage (MLSMSG_WARNING, ModuleName, &
             "Metadata write failed in populate_metadata_std " &
             & //trim(attrname)//trim(errmsg))
!       call announce_error(0, &
!       & "Error: metadata write failed in populate_metadata_std.", &
!       & error_number=returnStatus)
!       return
      end if
    end if

    hdfReturn = sfend(sdid)

    ! Annotate the file with the PCF

    call writePCF2Hdr(physical_filename, l2pcf%anText)

    returnStatus = pgs_met_remove() 

    if ( present(metadata_error)) metadata_error=module_error

  end subroutine Populate_metadata_std

  ! --------------------------------------  Populate_metadata_oth  -----

  subroutine Populate_metadata_oth ( HDF_FILE, MCF_FILE, L2pcf, &
    & NumQuantitiesPerFile, QuantityNames, Metadata_error )

    ! This is specially to write meta data for heterogeneous files
    ! It should work unchanged for the 'OTH' l2gp files (e.g. ML2OTH.001.MCF)
    ! and, with minor changes, for the l2aux files 

    !    USE InitPCFs, ONLY: L2PCF

    !Arguments

    integer :: HDF_FILE, MCF_FILE, NumQuantitiesPerFile
    type(PCFData_T) :: L2pcf
    character (len=*), dimension(:) :: QuantityNames
    integer, optional, intent(out) :: Metadata_error

    !Local Variables

    integer :: HdfReturn
    integer :: ReturnStatus
    integer :: Sdid
    character (len=132) :: Attrname, Errmsg

    integer, parameter :: INVENTORY=2, ARCHIVE=1
    character (len=PGSd_PC_FILE_PATH_MAX) :: physical_filename
!    character (len=132) :: attrname, errmsg
    integer :: Version, Indx
    character (len=*), parameter :: METAWR_ERR = &
      & 'Error writing metadata attribute '

    ! the group have to be defined as 49 characters long. The C interface is 50.
    ! The cfortran.h mallocs an extra 1 byte for the null character '\0/1, 
    ! therefore making the actual length of a 
    ! string pass of 50.

    character (len = PGSd_MET_GROUP_NAME_L) :: groups(PGSd_MET_NUM_OF_GROUPS)

    ! Externals

    integer, external :: PGS_MET_init, PGS_MET_setattr_d, &
      &  PGS_MET_setAttr_s, PGS_MET_getsetattr_d, PGS_MET_SETATTR_I, &
      &  PGS_MET_write, PGS_MET_remove

    !Executable code

    module_error = 0
    if ( present(metadata_error)) metadata_error=1

    version = 1
    if ( MCF_FILE > 0 ) then
      returnStatus = PGS_PC_GetReference (MCF_FILE, version , physical_filename)
    else
      returnStatus = PGSPC_W_NO_REFERENCE_FOUND
    end if
	
    if ( returnStatus /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for MCF_FILE in populate_metadata_oth" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for MCF_FILE in populate_metadata_oth.") 
      return
    end if

    version = 1
    returnStatus = PGS_PC_GetReference (HDF_FILE, version , physical_filename)

    if ( returnStatus /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for HDF_FILE in populate_metadata_oth" ) 
      call announce_error(0, &
      & "Error: failed to find PCF ref for HDF_FILE in populate_metadata_oth.") 
      return
    end if

    call first_grouping(HDF_FILE, MCF_FILE, l2pcf, groups)

    do indx=1, numquantitiesperfile

      call measured_parameter (HDF_FILE, &
        & QuantityNames(indx), groups, indx)

    end do

    call third_grouping (HDF_FILE, l2pcf, groups)

    sdid = sfstart (physical_fileName, DFACC_RDWR) 

    if ( sdid == -1 ) then
!     CALL MLSMessage (MLSMSG_Error, ModuleName, &
!       & "Failed to open the hdf file "//physical_fileName ) 
      call announce_error(0, &
      & "Error: failed to open the hdf file: "//TRIM(physical_fileName))
      return
    end if

    returnStatus = pgs_met_write (groups(INVENTORY), "coremetadata.0", sdid)

    if ( returnStatus /= PGS_S_SUCCESS .AND. &
         returnStatus /= PGSMET_W_METADATA_NOT_SET ) then 
      if ( returnStatus == PGSMET_W_METADATA_NOT_SET ) then 
!       CALL MLSMessage (MLSMSG_WARNING, ModuleName, &
!         & "Some of the mandatory parameters were not set" )
        call announce_error(0, &
        & "Error: Some of the mandatory parameters were not set in populate_metadata_oth.")
!       return
      else
        call Pgs_smf_getMsg (returnStatus, attrname, errmsg)
        call MLSMessage (MLSMSG_WARNING, ModuleName, &
             "Metadata write failed in populate_metadata_oth " &
             & //trim(attrname)//trim(errmsg))
!       call announce_error(0, &
!       & "Error: metdata write failed in populate_metadata_oth.", &
!       & error_number=returnStatus)
!       return
      end if
    end if

    hdfReturn = sfend(sdid)

! Annotate the file with the PCF

    call writePCF2Hdr(physical_filename, l2pcf%anText)

    returnStatus = pgs_met_remove() 

    if ( present(metadata_error)) metadata_error=module_error

  end subroutine Populate_metadata_oth

  ! -----------------------------------------------  Get_l2gp_mcf  -----

  subroutine Get_l2gp_mcf ( File_base, Mcf, L2pcf, Version )

  ! metadata configuration file (mcf) PCF number corresponding to l2gp number
  ! sdid
  ! Arguments
    character(len=*), intent(in) ::  File_base
    integer, intent(in), optional :: Version
    integer, intent(inout) ::        Mcf
    type(PCFData_T) :: l2pcf

! Local
    character (len=PGSd_PC_FILE_PATH_MAX) :: Sd_full
    character (len=NameLen) :: Sd_path
    character (len=NameLen) :: Sd_name

    character (len=PGSd_PC_FILE_PATH_MAX) :: Mcf_full
    character (len=NameLen) :: Mcf_path
    character (len=NameLen) :: Mcf_name
    character (len=NameLen) :: Mcf_pattern
    integer :: ReturnStatus, MyVersion, I

    character (len=1), parameter :: COMMA = ','

    ! Find species name
    ! assume sd_name is "*l2gp_species_"
    ! hence enclosed between "_" chars after an l2gp

    character (len=1), parameter :: Species_delimiter = '_'
    character (len=4), parameter :: L2gp = 'l2gp'
    logical, parameter :: DEBUG = .false.
	
    ! Begin

    if ( MCFFORL2GPOPTION == 1 ) then
      mcf = mcf+1
      return
    end if

    if ( DEBUG ) then
      call output('file_base: ', advance='no')
      call output(trim(file_base), advance='yes')
    end if

    if ( len(TRIM(file_base)) <= 0 ) then
      mcf=0
      return
    end if

    ! Get full file name for typical MCF file
    do i=mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end

      if ( present(version) ) then
	      myVersion=version
      else
	      myVersion = 1
      end if

      returnStatus = PGS_PC_GetReference(i, myVersion , mcf_full)

      if ( returnStatus == PGS_S_SUCCESS ) then 
	      exit
      end if

    end do

    if ( DEBUG ) then
      call output('returnStatus: ', advance='no')
      call output(returnStatus, advance='yes')
    end if

    if ( returnStatus /= PGS_S_SUCCESS ) then 
      mcf = 0
      return
    end if

    ! Split full_file_names into path+name

    call split_path_name ( mcf_full, mcf_path, mcf_name )
	
! If text matching not case sensitive, shift to lower case
    if ( .NOT. MCFCASESENSITIVE ) then
      sd_full = LowerCase(file_base)
      mcf_name = LowerCase(mcf_name)
    else
      sd_full = file_base
    end if

    if ( DEBUG ) then
      call output('mcf_full: ', advance='no')
      call output(trim(mcf_full), advance='yes')
      call output('mcf_path: ', advance='no')
      call output(trim(mcf_path), advance='yes')
      call output('mcf_name: ', advance='no')
      call output(trim(mcf_name), advance='yes')
    end if

    ! Get species name assuming e.g. '*l2gp_h2o'
    call split_path_name(sd_full, sd_path, sd_name, species_delimiter)

    if ( DEBUG ) then
      call output('sd_full: ', advance='no')
      call output(trim(sd_full), advance='yes')
      call output('sd_path: ', advance='no')
      call output(trim(sd_path), advance='yes')
      call output('sd_name: ', advance='no')
      call output(trim(sd_name), advance='yes')
    end if

    if ( len(trim(sd_name)) <= 0 ) then
      mcf=0
      return
    end if
	
    if ( MCFFORL2GPOPTION == 3 ) then

      ! get mcfspecies name from associative array
      ! if the species name not found in spec_keys, it will return ','
      call GetStringHashElement ( l2pcf%spec_keys, l2pcf%spec_hash, &
        & trim(sd_name), sd_full, .TRUE. )

      if ( DEBUG ) then
        call output('keys: ', advance='no')
        call output(trim(l2pcf%spec_keys), advance='yes')
        call output('hash: ', advance='no')
        call output(trim(l2pcf%spec_hash), advance='yes')
        call output('hash for species name: ', advance='no')
        call output(trim(sd_full), advance='yes')
      end if

      if ( trim(sd_full) == COMMA ) then
        mcf = 0
        return
      end if

      sd_name = trim(sd_full)

    end if

    ! Now try to find mcf file corresponding to species name
    ! assuming, e.g. '*h2o.*'

    if ( DEBUG ) then
      call output('loop over mcf files', advance='yes')
    end if

    do i=mlspcf_mcf_l2gp_start, mlspcf_mcf_l2gp_end

      if ( present(version) ) then
        myVersion=version
      else
        myVersion = 1
      end if

      returnStatus = PGS_PC_GetReference(i, myVersion, mcf_full)

      if ( returnStatus == PGS_S_SUCCESS ) then 
        call split_path_name(mcf_full, mcf_path, mcf_name)

        ! This gives 'o2h*'
        call split_path_name(Reverse(mcf_name), mcf_path, mcf_pattern, '.')
        if ( .NOT. MCFCASESENSITIVE ) then
          mcf_pattern = LowerCase(mcf_pattern)
        end if

        ! So reverse it to make h2o
        mcf_pattern = adjustl(Reverse(mcf_pattern))

        if ( DEBUG ) then
          call output('mcf_pattern: ', advance='no')
          call output(trim(mcf_pattern), advance='yes')
        end if
        ! Check that pattern matches species name
        ! Warning--this could give a false matching
        ! if the species name matches more than one mcf_pattern
        ! May wish to check for multiple matches later
        if ( index(trim(mcf_pattern), trim(sd_name)) > 0 ) then
          mcf = i
          return
	end if
      end if

    end do

    mcf = 0

  end subroutine Get_l2gp_mcf

! Lori's routines

  ! -----------------------------------------------  WriteMetaLog  -----
  subroutine WriteMetaLog ( Pcf, Metadata_error )

    ! Brief description of subroutine
    ! This subroutine writes metadata for the log file to a separate ASCII file.

    ! Arguments

    type( PCFData_T ), intent(in) :: Pcf
    integer, optional, intent(out) :: Metadata_error

    ! Parameters
    ! These are PCF numbers that *must* be in the PCF file

    integer, parameter :: ASCII_FILE = 101
    integer, parameter :: THE_LOG_FILE = 100

    ! Functions

    integer, external :: PGS_MET_init, PGS_MET_remove, PGS_MET_setattr_d
    integer, external :: PGS_MET_setattr_s, PGS_MET_write
    integer, external :: PGS_PC_getconfigdata

    ! Internal

    character (len=1) :: NullStr
    character (len=45) :: Sval
    character (len=32) :: Mnemonic
    character (len=FileNameLen) :: Physical_filename
    character (len=480) :: Msg, Msr
    character (len=PGSd_MET_GROUP_NAME_L) :: Groups(PGSd_MET_NUM_OF_GROUPS)

    integer :: Result, Indx, Version
    logical, parameter :: DEBUG=.FALSE.

    ! Begin
    module_error = 0
    if ( present(metadata_error)) metadata_error=1

    nullStr = ''

    version = 1
    result = PGS_PC_GetReference (mlspcf_mcf_l2log_start, version, &
      & physical_filename)

    if ( result /= PGS_S_SUCCESS ) then
!       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!         & "Failed to find the PCF reference for the Log mcf in WriteMetaLog" ) 
      call announce_error(0, &
        & "Error: failed to find PCF ref for Log mcf in WriteMetaLog.") 
      return
    end if

    ! Check that the PCF file contains the entries for the ASCII file
    ! and for the log file itself
!    version = 1
!    result = PGS_PC_GetReference (ASCII_FILE, version , physical_filename)
    result = pgs_pc_getconfigdata (ASCII_FILE, physical_filename)

    if ( result /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for the ASCII in WriteMetaLog" ) 
      call announce_error(0, &
        & "Error: failed to find PCF ref for the ASCII in WriteMetaLog.") 
      return
    end if

    version = 1
    result = PGS_PC_GetReference (THE_LOG_FILE, version , physical_filename)

    if ( result /= PGS_S_SUCCESS ) then
!     CALL MLSMessage (MLSMSG_Warning, ModuleName, &
!       & "Failed to find the PCF reference for the LOG in WriteMetaLog" ) 
      call announce_error(0, &
        & "Error: failed to find PCF ref for the Log in WriteMetaLog.") 
      return
    end if

! Initialize the MCF file

    if ( DEBUG ) then
      call output('Initialize the MCF file', advance='yes')
    end if

    result = pgs_met_init(mlspcf_mcf_l2log_start, groups)
    if ( result /= PGS_S_SUCCESS ) then
!     CALL MLSMessage(MLSMSG_Error, ModuleName, &
!       & 'Initialization error.  See LogStatus for details.')
      call announce_error(0, &
      & "metadata initialization error in WriteMetaLog.") 
      return
    end if

! Set PGE values

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalGranuleID", &
                               pcf%logGranID)

    call expandFileTemplate('$cycle', sval, cycle=pcf%cycle)
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "LocalVersionID", &
                               sval)

    indx = INDEX (PCF%startUTC, "T")

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeBeginningDate", PCF%startUTC(1:indx-1))
!    sval= '00:00:00.000000'
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeBeginningTime", PCF%startUTC(indx+1:))
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeEndingDate", PCF%endUTC(1:indx-1))
!    sval= '23:59:59.999999'
    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), &
                               "RangeEndingTime", PCF%endUTC(indx+1:))

    result = pgs_met_setAttr_s(groups(INVENTORYMETADATA), "PGEVersion", &
                               pcf%PGEVersion)

    if ( result /= PGS_S_SUCCESS ) then
!         call Pgs_smf_getMsg(result, mnemonic, msg)
!         msr = mnemonic // ':  ' // msg
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      call announce_error(0, &
      & "Error in setting PGEVersion attribute in WriteMetaLog.", &
      &  error_number=result)
        return
    end if

    ! Write the metadata and their values to an ASCII file

    if ( DEBUG ) then
      call output('Write the metadata and their values to an ASCII file', &
      &  advance='yes')
    end if

    result = pgs_met_write(groups(1), nullStr, ASCII_FILE)

    if ( result /= PGS_S_SUCCESS ) then
       call Pgs_smf_getMsg(result, mnemonic, msg)
       msr = "Error: failed to write metadata in WriteMetaLog." &
       & // trim(mnemonic) // ':  ' // trim(msg)
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
!    call announce_error(0, &
!    & "Error: failed to write metadata in WriteMetaLog.", &
!    & error_number=result)
!                     return
    end if

    result = pgs_met_remove()

    if ( present(metadata_error)) metadata_error=module_error

   end subroutine WriteMetaLog

  ! -----------------------------------------  ExpandFileTemplate  -----
   subroutine ExpandFileTemplate ( Template, Filename, Level, Version, Cycle, &
     &                            Day )

     ! Brief description of subroutine
     ! This subroutine expands the version, cycle, and day fields in a file
     ! name template.  Note: CYCLE is a STRING here, which is how it's read
     ! from the PCF.

     ! Arguments

     character (len=*), intent(in) :: Template

     character (len=*), intent(in), optional :: Cycle, Day, Level, Version

     character (len=*), intent(out) :: Filename

    ! Parameters

    ! Functions

    ! Variables

    character (len=4) :: Field, zCy

    integer :: I, iCy, Indx, NumFields

    ! Initializations

    numFields = 4

    fileName = template

    ! Loop through the expandable fields in the template

    do i = 1, numFields

      ! Search for $

      indx = INDEX(fileName,'$')

      ! Exit, if there are no expandable fields 

      if ( indx == 0) exit

      ! Match the field name to the input argument

      field = fileName(indx:indx+3)

      if ( field == '$lev' ) then

        if ( present(level) ) then
          fileName = fileName(:(indx-1)) // TRIM(level) // &
            &        fileName((indx+6):)
        else
          call MLSMessage(MLSMSG_Error, ModuleName, 'Input level &
            &required to expand the template.')
        end if

      else if ( field == '$ver' ) then

        if ( present(version) ) then
          fileName = fileName(:(indx-1)) // TRIM(version) // &
            &        fileName((indx+8):)
        else
          call MLSMessage(MLSMSG_Error, ModuleName, 'Input version &
            &required to expand the template.')
        end if

      else if ( field == '$cyc' ) then

        if ( present(cycle) ) then

        ! Convert from CHARACTER to INTEGER

          read ( cycle, '(I2)' ) iCy

          ! Add a leading zero, if less than 10

          if ( iCy < 10 ) then
            zCy = '0' // trim(cycle)
          else
            zCy = cycle
          end if

          fileName = fileName(:(indx-1)) // 'C' // TRIM(zCy) // &
                            fileName((indx+6):)

        else

          call MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Input cycle required to expand the template.')

        end if

      else if ( field == '$day' ) then

        if ( present(day) ) then
          fileName = fileName(:(indx-1)) // TRIM(day) // &
            &        fileName((indx+4):)
        else
          call MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Input day required to expand the template.')
        end if

      end if

    end do

   end subroutine ExpandFileTemplate

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
    character (LEN=132) :: attrname, errmsg
    integer :: Toolbox_error_num

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = .not. use_toolkit

    if ( .not. just_print_it ) then
      module_error = max(module_error,PENALTY_FOR_NO_METADATA)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
        call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( " Caused the following error:", advance='yes', &
       & from_where=ModuleName)
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if ( present(error_number) ) then
        call output ( 'Error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
        Toolbox_error_num = error_number
        call Pgs_smf_getMsg (Toolbox_error_num, attrname, errmsg)
        call MLSMessage (MLSMSG_WARNING, ModuleName, &
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

end module WriteMetadata 
! $Log$
! Revision 2.13  2001/04/20 23:52:15  vsnyder
! Consider the penalty from MLSL2Options in Announce_Error.  Numerous
! cosmetic changes.  Remove dump of tree node from Announce_Error.
!
! Revision 2.12  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.11  2001/04/16 23:49:11  pwagner
! Tiny change to announce_error
!
! Revision 2.10  2001/04/16 17:43:35  pwagner
! mcf text matching case sensitive if MCFCASESENSITIVE
!
! Revision 2.9  2001/04/13 23:47:26  pwagner
! Writes multiple measuredcontainers properly
!
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
