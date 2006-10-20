! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GLOBAL_SETTINGS

  use MLSCommon, only: FILENAMELEN

  implicit NONE

  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! LEAPSECFILENAME                 For time conversions (if avoiding PCF)
! brightObjects                   As determined by l1b

!     (subroutines and functions)
! SET_GLOBAL_SETTINGS             Get global settings from l2cf
! === (end of toc) ===

  public :: SET_GLOBAL_SETTINGS

  character(LEN=FileNameLen), public :: LEAPSECFILENAME = ''

  ! These must be values consistent with level 1
  ! (Some canny coding below could make this more robust)
  integer, parameter :: BO_NAMEDIMS = 14
  integer, parameter :: BO_NAMELEN = 14

  ! This next should be large enough to hold the entire list of BO names
  integer, parameter :: BONAMELISTLEN = 256
  character(len=BONAMELISTLEN), public, save :: brightObjects = &
    & 'MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO, ' // &
    & 'MOON, SUN, SOLAR-BARYCNTR, EARTH-BARYCNTR, GALACTICCENTER'

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, private :: ERROR

contains

  subroutine SET_GLOBAL_SETTINGS ( ROOT, ForwardModelConfigDatabase, &
    & FGrids, l2gpDatabase, DirectDatabase, processingRange, filedatabase )

    use dates_module, only: utc_to_yyyymmdd
    use DirectWrite_m, only: DirectData_T, &
      & AddDirectToDatabase, Dump, SetupNewDirect
    use DumpCommand_m, only: DumpCommand
    use EmpiricalGeometry, only: INITEMPIRICALGEOMETRY
    use FGrid, only: AddFGridToDatabase, CreateFGridFromMLSCFInfo, FGrid_T
    use ForwardModelConfig, only: AddForwardModelConfigToDatabase, Dump, &
      & ForwardModelConfig_T
    use ForwardModelSupport, only: ConstructForwardModelConfig, &
      & ForwardModelGlobalSetup, CreateBinSelectorFromMLSCFInfo
    use INIT_TABLES_MODULE, only: F_FILE, F_TYPE, &
      & L_L2GP, L_L2DGG, L_L2FWM, &
      & P_BRIGHTOBJECTS, &
      & P_CYCLE, P_ENDTIME, P_INSTRUMENT, &
      & P_LEAPSECFILE, P_OUTPUT_VERSION_STRING, P_PFAFILE, P_STARTTIME, &
      & S_BINSELECTOR, S_DIRECTWRITEFILE, S_DUMP, S_EMPIRICALGEOMETRY, &
      & S_FGRID, S_FLUSHPFA, S_FORWARDMODEL, S_FORWARDMODELGLOBAL, &
      & S_L1BOA, S_L1BRAD, S_L2PARSF, S_MAKEPFA, S_PFADATA, S_READPFA, &
      & S_TGRID, S_TIME, S_VGRID, S_WRITEPFA
    use intrinsic, only: l_hdf, l_swath
    use L1BData, only: L1BData_T, NAME_LEN, PRECISIONSUFFIX, &
      & AssembleL1BQtyName, DeallocateL1BData, Dump, FindMaxMAF, &
      & l1bradSetup, l1boaSetup, ReadL1BAttribute, ReadL1BData 
    use L2GPData, only: L2GPDATA_T
    use L2ParInfo, only: parallel
    use L2PC_M, only: AddBinSelectorToDatabase, BinSelectors
    use MLSCommon, only: R8, FileNameLen, MLSFile_T, NameLen, &
      & TAI93_Range_T
    use MLSFiles, only: FILENOTFOUND, HDFVERSION_5, &
      & AddFileToDataBase, GetPCFromRef, GetMLSFileByName, GetMLSFileByType, &
      & InitializeMLSFile, mls_CloseFile, mls_OpenFile, split_path_name
    use MLSL2Options, only: LEVEL1_HDFVERSION, SPECIALDUMPFILE, &
      & STOPAFTERSECTION, Toolkit
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSPCF2, only: mlspcf_l2gp_start, mlspcf_l2gp_end, &
      & mlspcf_l2dgm_start, mlspcf_l2dgm_end, mlspcf_l2fwm_full_start, &
      & mlspcf_l2fwm_full_end, &
      & mlspcf_l2dgg_start, mlspcf_l2dgg_end
    use MLSStrings, only: hhmmss_value, lowerCase
    use MLSStringLists, only: Array2List, catLists, SWITCHDETAIL
    use MLSSignals_m, only: INSTRUMENT
    use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
    use OUTPUT_M, only: BLANKS, OUTPUT, revertoutput, switchOutput
    use PFAData_m, only: Get_PFAdata_from_l2cf, Flush_PFAData, Make_PFAData, &
      & Read_PFAData, Write_PFAData
    use PFADataBase_m, only: Process_PFA_File
    use PCFHdr, only: GlobalAttributes, FillTAI93Attribute
    use readAPriori, only: APrioriFiles
    use SDPToolkit, only: max_orbits, mls_utctotai
    use String_Table, only: Get_String
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE, &
    & DUMP_TREE_NODE, SOURCE_REF
    use TREE_TYPES, only: N_EQUAL, N_NAMED
    use VGrid, only: CreateVGridFromMLSCFInfo
    use VGridsDatabase, only: AddVGridToDatabase, VGrids
    use WriteMetadata, only: L2PCF

    ! placed non-alphabetically due to Lahey internal compiler error
    ! (How much longer must we endure these onerous work-arounds?)
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributeInFile
    use HDF5, only: h5gclose_f, h5gopen_f

    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(ForwardModelConfig_T), dimension(:), pointer :: &
      & ForwardModelConfigDatabase
    type ( fGrid_T ), pointer, dimension(:) :: FGrids
    type ( l2gpData_T), dimension(:), pointer :: L2GPDATABASE
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    ! Local variables
    character(len=BO_NAMELEN), dimension(BO_NAMEDIMS) :: BO_names
    integer :: DetailReduction
    integer :: Details             ! How much info about l1b files to dump
    type (MLSFile_T) :: DirectFile
    logical :: GOT(3)              ! Used non-canonically--a bad practice
    integer :: I, J                ! Index of son, grandson of root
    ! integer :: INPUT_VERSION_STRING  ! Sub_rosa index
    logical ::  ItExists
    type (L1BData_T) :: l1bField   ! L1B data
    integer :: L1BFLAG
    real(r8) :: MINTIME, MAXTIME   ! Time Span in L1B file data
    integer :: NAME                ! Sub-rosa index of name of vGrid or hGrid
    integer :: NOMAFS              ! Number of MAFs of L1B data read
    integer :: NUMFILES
    integer :: OUTPUT_VERSION_STRING  ! Sub_rosa index
    integer :: param_id            ! e.g., p_brightObjects
    integer :: ReturnStatus        ! non-zero means trouble
    integer :: SON                 ! Son of root
    integer :: spec_id             ! e.g., s_binSelector
    logical :: StopEarly
    integer :: Sub_rosa_index
    integer :: The_HDF_version     ! 4 or 5 (corresp. to hdf4 or hdf5)
    logical :: TIMING              ! For S_Time
    logical :: StartTimeIsAbsolute, stopTimeIsAbsolute
    integer :: Status
    real :: T1, T2                 ! For S_Time
    real(r8) :: Start_time_from_1stMAF, End_time_from_1stMAF
    logical :: wasAlreadyOpen
    character(len=NameLen) :: Name_string
    character(len=NameLen) :: End_time_string, Start_time_string
    character(len=FileNameLen) :: FilenameString
    character (len=name_len) :: QUANTITY
    character(len=*), parameter :: Time_conversion='(F32.0)'
    character(len=Name_Len) :: l1bItemName
    integer :: OrbNum(max_orbits)
    real(r8) :: OrbPeriod(max_orbits)
    type(MLSFile_T), pointer :: L1BFile

    timing = section_times
    if ( timing ) call time_now ( t1 )
    stopEarly = ( &
      & index('global_setting,chunk_divide', &
      & lowercase( trim(stopAfterSection) ) ) > 0 )
    stopEarly = ( stopearly .and. stopAfterSection /= ' ' )
    error = 0
    got = .false.
    startTimeIsAbsolute = .false.
    stopTimeIsAbsolute = .false.
    LeapSecFileName = ''

    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    i = index(switches, 'glo')
    if ( i /= 0 ) then
      if ( switches(i+3:i+3) >= '0' .and. switches(i+3:i+3) <= '9' ) then
        details = iachar(switches(i+3:i+3)) - iachar('0') - 2
      else
        details = -3
      end if
    else
      details = -4
    end if

    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    elseif ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    endif

    do i = 2, nsons(root)-1 ! Skip names at beginning and end of section
      son = subtree(i,root)
      if ( node_id(son) == n_equal ) then
        sub_rosa_index = sub_rosa(subtree(2,son))
        param_id = decoration(subtree(1,son))
        if ( TOOLKIT .and. &
          & any( param_id == &
          & (/ p_output_version_string, p_cycle, p_starttime, p_endtime, &
          & p_leapsecfile /) ) ) then
          call announce_error(0, &
            & '*** l2cf parameter global setting ignored ***', &
            & just_a_warning = .true.)
          cycle
        endif
        select case ( param_id )
        ! This will allow us to use different names from the toolkit
        ! (Now why would you want to do that?)
        case ( p_brightObjects )
          call get_string ( sub_rosa_index, brightObjects, strip=.true. )
          got(3) = .true.
        case ( p_output_version_string )
          output_version_string = sub_rosa_index
          call get_string ( output_version_string, l2pcf%PGEVersion, strip=.true. )
        case ( p_instrument )
          instrument = decoration(subtree(2,son))
        case ( p_cycle )
          call get_string ( sub_rosa_index, l2pcf%cycle, strip=.true. )
        case ( p_starttime )
          got(1) = .true.
          call get_string ( sub_rosa_index, name_string, strip=.true. )
          start_time_string = name_string
          if ( index(name_string, ':') > 0 ) then
            start_time_from_1stMAF = hhmmss_value(name_string, returnStatus)
            error = max(error,returnStatus)
            startTimeIsAbsolute = .false.
            l2pcf%startutc = start_time_string
          else
            read ( name_string, * ) start_time_from_1stMAF
            startTimeIsAbsolute = .true.
          end if
        case ( p_endtime )
          got(2) = .true.
          call get_string ( sub_rosa_index, name_string, strip=.true. )
          end_time_string = name_string
          if ( index(name_string, ':') > 0 ) then
            end_time_from_1stMAF = hhmmss_value(name_string, returnStatus)
            error = max(error,returnStatus)
            stopTimeIsAbsolute = .false.
            l2pcf%endutc = end_time_string
          else
            read ( name_string, * ) end_time_from_1stMAF
            stopTimeIsAbsolute = .true.
          end if
        case ( p_leapsecfile )
          call get_string ( sub_rosa_index, LeapSecFileName, strip=.true. )
          inquire(file=trim(LeapSecFileName), exist=itExists)
          if ( .not. itExists ) then
            call announce_error(0, &
            & '*** Leap Second File ' // trim(LeapSecFileName) // &
            & ' not found', &
            & just_a_warning = .false.)
            call MLSMessage ( MLSMSG_Error, ModuleName, &                      
            & '(Please check file name and path)' )    
          end if
        case ( p_PFAFile )
          do j = 2, nsons(son)
            if ( process_PFA_File ( sub_rosa(subtree(j,son)), &
              & source_ref(subtree(j,son)) ) /= 0 ) continue
          end do
        case default
          call announce_error(son, 'unrecognized global settings parameter')
        end select
      else
        if ( node_id(son) == n_named ) then
          name = sub_rosa(subtree(1,son))
          son = subtree(2,son)
        else
          name = 0
        end if
        spec_id = get_spec_id(son)
        if ( TOOLKIT .and. &
          & any( spec_id == &
          & (/ s_l1boa, s_l1brad, s_l2parsf /) ) ) then
          call announce_error(0, &
            & '*** l2cf spec global setting ignored ***', &
            & just_a_warning = .true.)
          cycle
        endif
        select case ( spec_id )
        case ( s_binSelector )
          call decorate (son, AddBinSelectorToDatabase ( &
            & binSelectors, CreateBinSelectorFromMLSCFInfo ( son ) ) )
        case ( s_directWriteFile )
          call decorate (son, AddDirectToDatabase ( &
            & DirectDatabase, &
            & CreateDirectTypeFromMLSCFInfo ( son, DirectFile ) ) )
          numFiles = AddFileToDataBase(fileDataBase, DirectFile)
        case ( s_dump )
          if ( error == 0 ) then
            call dumpCommand ( son, forwardModelConfigs=forwardModelConfigDatabase )
          else
            call announce_error ( subtree(1,son), &
              & 'Preceeding errors prevent doing a dump here.' )
          end if
        case ( s_empiricalGeometry )
          call InitEmpiricalGeometry ( son )
        case ( s_flushPFA )
          call flush_PFAData ( son, status )
          error = max(error,status)
        case ( s_fgrid )
          call decorate ( son, AddFGridToDatabase ( fGrids, &
            & CreateFGridFromMLSCFInfo ( name, son ) ) )
        case ( s_forwardModelGlobal ) !??? Begin temporary stuff for l2load
          if ( .not. stopEarly ) call forwardModelGlobalSetup ( son, returnStatus, &
            & fileDataBase )
          error = max(error, returnStatus)
        case ( s_forwardModel )
          if ( .not. stopEarly ) call decorate (son, AddForwardModelConfigToDatabase ( &
            & forwardModelConfigDatabase, &
            & ConstructForwardModelConfig ( name, son, .true. ) ) )
        case ( s_l1boa )
          the_hdf_version = LEVEL1_HDFVERSION
          call l1boaSetup ( son, filedatabase, F_FILE, hdfVersion=the_hdf_version )
          if ( index(switches, 'pro') /= 0 ) then  
            sub_rosa_index = sub_rosa(subtree(2,subtree(2, son)))
            call get_string ( sub_rosa_index, FilenameString, strip=.true. )
            call proclaim(FilenameString, 'l1boa', &                   
            & hdfVersion=the_hdf_version) 
          end if
          L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
          if ( associated(L1BFile) ) then
            call ReadL1BAttribute (L1BFile, OrbNum, 'OrbitNumber', &
               & l1bFlag)
               if ( l1bFlag == -1 ) then
                  GlobalAttributes%OrbNum = -1
               else
                  GlobalAttributes%OrbNum = OrbNum
               end if
              call ReadL1BAttribute (L1BFile, OrbPeriod, 'OrbitPeriod', &
               & l1bFlag)
               if ( l1bFlag == -1 ) then
                  GlobalAttributes%OrbPeriod = -1.0
               else
                  GlobalAttributes%OrbPeriod = OrbPeriod
               end if
            if ( details > -4 ) call output &
              & ('finished readL1BAttribute in global_setting', advance='yes')
          else
            call announce_error(0, &
            & '*** Unable to find L1BOA file ***', &
            & just_a_warning = .true.)
          end if
        case ( s_l1brad )
          the_hdf_version = LEVEL1_HDFVERSION
          call l1bradSetup ( son, filedatabase, F_FILE, &
            & hdfVersion=the_hdf_version )
          sub_rosa_index = sub_rosa(subtree(2,subtree(2, son)))
          call get_string ( sub_rosa_index, FilenameString, strip=.true. )
          if ( index(switches, 'pro') /= 0 ) then  
            call proclaim(FilenameString, 'l1brad', &                   
            & hdfVersion=the_hdf_version) 
          end if
          L1BFile => GetMLSFileByName(filedatabase, FilenameString)
          if ( .not. associated(L1BFile) ) then
            call announce_error(0, &
            & '*** Unable to find L1BRAD file ***' // trim(FilenameString), &
            & just_a_warning = .true.)
          end if
        case ( s_l2parsf )
          sub_rosa_index = sub_rosa(subtree(2,subtree(2, son)))
          call get_string ( sub_rosa_index, FilenameString, strip=.true. )
          parallel%stagingFile = FilenameString
          if ( index(switches, 'pro') /= 0 ) then  
            call proclaim(FilenameString, 'l2 parallel staging file') 
          end if
        case ( s_makePFA )
          call Make_PFAData ( son, returnStatus )
          error = max(error, returnStatus)
        case ( s_pfaData )
          call Get_PFAdata_from_l2cf ( son, name, returnStatus )
          error = max(error, returnStatus)
        case ( s_readPFA )
          call read_PFAdata ( son, name, returnStatus )
          error = max(error, returnStatus)
        case ( s_writePFA )
          call write_PFAdata ( son, returnStatus )
          error = max(error, returnStatus)
        case ( s_time )
          if ( timing ) then
            call sayTime
          else
            call time_now ( t1 )
            timing = .true.
          end if
        case ( s_tGrid, s_vGrid )
          call decorate ( son, AddVGridToDatabase ( vGrids, &
            & CreateVGridFromMLSCFInfo ( name, son, l2gpDatabase, returnStatus ) ) )
          error = max(error, returnStatus)
        case default
          call announce_error(son, 'unrecognized global settings spec')
        end select
      end if
    end do

    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &                      
      & 'l1boa File not found--hope you dont need one' )
      return
    endif

    the_hdf_version = &
      & L1BFile%HDFVersion
    if ( the_hdf_version == FILENOTFOUND ) then                                          
      call MLSMessage ( MLSMSG_Error, ModuleName, &                      
      & 'File not found; make sure the name and path are correct' &
      & // trim(L1BFile%Name) )
    else if ( the_hdf_version <= 0 ) then                                          
      call MLSMessage ( MLSMSG_Error, ModuleName, &                      
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )    
    end if
    ! add maf offsets to start, end times
    ! or convert them to tai93
    ! This is optional way to define processingRange if using PCF
    ! It becomes mandatory if not using PCF
    if ( LeapSecFileName == '' &
     & .and. &
     & (got(1) .or. got(2) .or. .not. TOOLKIT) &
     & ) then

      ! 1st--check that have L1BOA
      if ( L1BFile%FileID%f_id < 1 ) then
        call announce_error(son, &
          & 'L1BOA file required by global data--but not set')
      end if

      quantity = 'MAFStartTimeTAI'
      l1bItemName = AssembleL1BQtyName ( quantity, the_hdf_version, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, dontPad=.true.)
      if ( l1bFlag==-1 ) then
        call announce_error(son, &
          & 'unrecognized MAFStarttimeTAI in L1BOA file')
        minTime = 0.
        maxTime = 0.
      else
        minTime = l1bField%dpField(1,1,1)
        maxTime = l1bField%dpField(1,1,noMAFs) ! This is start time of last MAF
      end if
      call DeallocateL1BData ( l1bField )
    end if

    if ( startTimeIsAbsolute ) then
      processingRange%startTime = start_time_from_1stMAF
    else
      if ( LeapSecFileName /= '' ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), start_time_string, &
        & processingrange%starttime)
        if ( returnStatus /= 0 ) then
          call announce_error(0, &
          & 'Error converting start time in mls_utctotai; code number: ')
          call output(returnStatus, advance='yes')
        end if
      else if ( got(1) ) then
        processingrange%starttime = minTime + start_time_from_1stMAF
      else if ( .not. TOOLKIT ) then
        processingrange%starttime = minTime
      end if
    end if

    if ( stopTimeIsAbsolute ) then
      processingRange%endTime = end_time_from_1stMAF
    else
      if ( LeapSecFileName /= '' ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), end_time_string, &
        & processingrange%endtime)
        if ( returnStatus /= 0 ) then
          call announce_error(0, &
          & 'Error converting end time in mls_utctotai; code number: ')
          call output(returnStatus, advance='yes')
        end if
      else if ( got(2) ) then
        processingrange%endtime = minTime + end_time_from_1stMAF
      else if ( .not. TOOLKIT ) then
        processingrange%endtime = maxTime + 1.0
      end if
    end if

    if ( .not. TOOLKIT ) then
      ! Store appropriate user input as global attributes
      GlobalAttributes%StartUTC = l2pcf%StartUTC
      GlobalAttributes%EndUTC = l2pcf%EndUTC
      GlobalAttributes%PGEVersion = l2pcf%PGEVersion
      if ( LeapSecFileName /= '' ) call FillTAI93Attribute ( LeapSecFileName )
      ! We don't check on returnStatus--dateless or absolute utc are ok
      call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
        & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
        & GlobalAttributes%GranuleDay) 
      ! We'll possibly need the first and last MAF counter numbers, especially
      ! if gaps occur
      ! For now, just look for them in l1boa
      ! Later you may look also in l1brad files
      GlobalAttributes%LastMAFCtr = FindMaxMAF ( (/L1BFile/), &
        & GlobalAttributes%FirstMAFCtr )
    end if

    ! Have we overridden the Bright Object names? Can we find them in l1boa?
    if ( got(3) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &                      
      & 'You overrode names of Bright objects found in l1boa file' )    
    elseif( L1BFile%hdfVersion /= HDFVERSION_5 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &                      
      & 'Names of Bright objects missing from (hdf4) l1boa file' )    
    elseif( .not. IsHDF5AttributeInFile(L1BFile%name, 'BO_name') ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &                      
      & 'Names of Bright objects missing from l1boa file' )    
    else
      wasAlreadyOpen = L1BFile%stillOpen
      if ( .not. wasAlreadyOpen ) call mls_OpenFile(L1BFile)
      L1BFile%fileID%sd_id = 0 ! So we don't look here for the attribute
      ! We still need to open the root group '/'
      call h5gopen_f ( L1BFile%fileID%f_id, '/', L1BFile%fileID%grp_id, ReturnStatus )
      call GetHDF5Attribute( L1BFile, 'BO_name', BO_names )
      call Array2List ( BO_names, BrightObjects )
      call h5gclose_f ( L1BFile%fileID%grp_id, ReturnStatus )
      if ( .not. wasAlreadyOpen ) call mls_CloseFile(L1BFile)
    endif

    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ! Perhaps dump global settings
    if ( details > -4 ) &
      & call dump_global_settings( processingRange, filedatabase, DirectDatabase, &
      & ForwardModelConfigDatabase, LeapSecFileName, details-detailReduction )

    if ( APrioriFiles%dao // AprioriFiles%ncep // AprioriFiles%geos5 &
      &  == ' ' ) then
      GlobalAttributes%MiscNotes = catLists(GlobalAttributes%MiscNotes, &
        & 'No gmao or ncep files--falling back to climatology', '\')
    elseif ( APrioriFiles%dao // AprioriFiles%geos5 == ' ' ) then
      GlobalAttributes%MiscNotes = catLists(GlobalAttributes%MiscNotes, &
        & 'No gmao files--falling back to ncep', '\')
    endif

    if ( error /= 0 ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
      & 'Problem with global settings section')

    if ( specialDumpFile /= ' ' ) &
      & call revertOutput
    if ( toggle(gen) ) then
      call trace_end ( 'SET_GLOBAL_SETTINGS' )
    end if
    if ( timing ) call sayTime

  contains

    ! ---------------------------------------------  Announce_Error  -----
    subroutine Announce_Error ( Lcf_where, Full_message, Use_toolkit, &
      & Error_number, just_a_warning )

      use LEXER_CORE, only: PRINT_SOURCE

      ! Arguments

      integer, intent(in) :: Lcf_where
      character(LEN=*), intent(in) :: Full_message
      logical, intent(in), optional :: Use_toolkit
      integer, intent(in), optional :: Error_number
      logical, intent(in), optional :: just_a_warning

      ! Local
      logical :: Just_print_it, my_warning
      logical, parameter :: Default_output_by_toolkit = .true.

      just_print_it = .not. default_output_by_toolkit
      if ( present(use_toolkit) ) just_print_it = .not. use_toolkit
      if ( present(just_a_warning) ) then
        my_warning = just_a_warning
      else
        my_warning = .false.
      end if

      if ( .not. just_print_it ) then
       if ( .not. my_warning ) then
        error = max(error,1)
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

        end if

        call output ( trim(full_message), advance='yes', &
          & from_where=ModuleName)

        if ( present(error_number) ) then
         if ( my_warning ) then
          call output ( 'Warning number ', advance='no' )
        else
          call output ( 'Error number ', advance='no' )
        end if

          call output ( error_number, places=9, advance='yes' )
        end if
      else

        if ( .not. my_warning ) then
        call output ( '***Error in module ' )
        call output ( ModuleName, advance='yes' )
        end if

        call output ( trim(full_message), advance='yes' )
        if ( present(error_number) ) then
         if ( my_warning ) then
          call output ( 'Warning number ' )
        else
          call output ( 'Error number ' )
         end if
          call output ( error_number, advance='yes' )
        end if
      end if

    end subroutine Announce_Error

    ! ------------------------------------------  dump_global_settings  -----
    subroutine dump_global_settings ( processingRange, &
      & filedatabase, DirectDatabase, ForwardModelConfigDatabase, &
      & LeapSecFileName, details )

      use MLSStringLists, only: NumStringElements, StringElement
      ! Dump info obtained during OpenAndInitialize and global_settings:
      ! L1B databse
      ! L1OA file
      ! Start and end times
      ! output version
      ! cycle number
      ! logfile name

      ! Arguments
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (TAI93_Range_T) :: processingRange ! Data processing range
      type (DirectData_T), dimension(:), pointer :: DirectDatabase
      type(ForwardModelConfig_T), dimension(:), pointer :: &
        & ForwardModelConfigDatabase

      ! The following dtermines the level of detail to expose:
      ! -1 Skip even counterMAF
      ! -2 Skip all but name (default)
      ! >0 Dump even multi-dim arrays
      integer, intent(in) :: details
      character(len=*) LeapSecFileName

      ! Local
      logical, parameter :: countEmpty = .true.
      type (L1BData_T) :: l1bData   ! L1B dataset
      integer ::                              i, version, NoMAFs, IERR
      character (len=*), parameter ::         TIME_FORMAT = '(1pD18.12)'
      ! This next is in case we're to dump at greatest possible detail
      character (len=NAME_LEN), parameter ::  BASE_QUANT_NAME = &
                                      &       'R2:190.B3F:N2O.S2.FB25-3'
  !                                    &       'R1A:118.B1F:PT.S0.FB25-1'
      character (len=LEN(BASE_QUANT_NAME)) :: l1b_quant_name
      logical, parameter ::                   DUMPPRECISIONTOO = .true.
      integer ::  hdfVersion
      character(len=Name_Len) :: l1bItemName
      type(MLSFile_T), pointer :: L1BFile

      ! Begin
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      hdfVersion = L1BFile%HDFVersion
      if ( hdfversion <= 0 ) &                                          
        & call MLSMessage ( MLSMSG_Error, ModuleName, &                    
        & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )  
      version = 1

      call output ( '============ Global Settings ============', advance='yes' )
      call output ( ' ', advance='yes' )

      call output ( 'L1B database:', advance='yes' )

        do i = 1, size(filedatabase)
         L1BFile => filedatabase(i)
         if ( L1BFile%content == 'l1brad' ) then
  	        call output ( 'fileid:   ' )
           call output ( L1BFile%FileID%f_id, advance='yes' )
           call output ( 'name:   ' )
    	     call output ( TRIM(L1BFile%name), advance='yes' )
           if ( details > -2 ) then
             l1b_quant_name = BASE_QUANT_NAME
             l1bItemName = AssembleL1BQtyName ( l1b_quant_name, hdfVersion, .false. )
             call ReadL1BData ( L1BFile, l1bItemName, L1bData, &
              & NoMAFs, IERR, NeverFail=.true. )
             if ( IERR == 0 ) then
               call Dump(l1bData, details )
               call DeallocateL1BData ( l1bData )
             else
               call output ( 'Error number  ' )
               call output ( IERR )
               call output ( ' while reading quantity named  ' )
               call output ( trim(l1b_quant_name), advance='yes' )
             end if
           end if
           if ( details > -2 .and. DUMPPRECISIONTOO ) then
             l1b_quant_name = trim(BASE_QUANT_NAME) // PRECISIONSUFFIX
             l1bItemName = AssembleL1BQtyName ( l1b_quant_name, hdfVersion, .false. )
             call ReadL1BData ( L1BFile, l1bItemName, L1bData, &
              & NoMAFs, IERR, NeverFail=.true. )
             if ( IERR == 0 ) then
               call Dump(l1bData, details )
               call DeallocateL1BData ( l1bData )
             else
               call output ( 'Error number  ' )
               call output ( IERR )
               call output ( ' while reading quantity named  ' )
               call output ( trim(l1b_quant_name), advance='yes' )
             end if
           end if
         elseif ( L1BFile%content == 'l1boa' ) then
           call output ( ' ', advance='yes' )
           call output ( 'L1OA file:', advance='yes' )

           call output ( 'fileid:   ' )
           call output ( L1BFile%FileID%f_id, advance='yes' )
           call output ( 'name:   ' )
           call output ( TRIM(L1BFile%Name), advance='yes' )
         end if
        end do



      call dump(DirectDatabase, Details)
      call output ( ' ', advance='yes' )
      call output ( 'Start Time:   ' )
      call output ( l2pcf%startutc, advance='yes' )

      call output ( 'End Time:     ' )
      call output ( l2pcf%endutc, advance='yes' )

      call output ( 'Start Time (tai):   ' )
      call output ( processingrange%starttime, format=TIME_FORMAT, advance='yes' )

      call output ( 'End Time (tai):     ' )
      call output ( processingrange%endtime, format=TIME_FORMAT, advance='yes' )

      call output ( 'Processing Range:     ' )
      call output ( processingrange%endtime-processingrange%starttime, advance='yes' )

      if ( LeapSecFileName /= '' ) then
        call output ( 'Leap Seconds File:   ' )
        call output ( trim(LeapSecFileName), advance='yes' )
      end if

      call output ( 'PGE version:   ' )
      call output ( l2pcf%PGEVersion, advance='yes' )

      ! call output ( 'input version:   ' )
      ! call output ( l2pcf%InputVersion, advance='yes' )

      call output ( 'cycle:   ' )
      call output ( l2pcf%cycle, advance='yes' )

      call output ( 'Log file name:   ' )
      call output ( TRIM(l2pcf%logGranID), advance='yes' )

      call output ( 'l2gp species name keys:   ' )
      call output ( TRIM(l2pcf%spec_keys), advance='yes' )

      call output ( 'corresponding mcf hash:   ' )
      call output ( TRIM(l2pcf%spec_hash), advance='yes' )

      call output ( 'Bright Objects:   ', advance='yes' )
      call output ( 'bit #             Cause   ', advance='yes' )
      call output ( '0' )
      call blanks (2)
      call output ( 'Moon in space port', advance='yes' )
      do i=1, NumStringElements(brightObjects, countEmpty)
        call output ( i )
        call blanks (2)
        call output ( TRIM(stringElement(brightObjects, i, .true.)), &
          & advance='yes' )
      enddo

      if ( details > -2 ) then
        call output ( ' ', advance='yes' )
        call dump(ForwardModelConfigDatabase, details=9, skipPFA=.true.)
      endif
      call output ( '============ End Global Settings ============', advance='yes' )

    end subroutine dump_global_settings

    ! ---------------------------------------------  proclaim  -----
    subroutine proclaim ( Name, l1_type, hdfVersion )
      character(LEN=*), intent(in)   :: Name
      character(LEN=*), intent(in)   :: l1_type
      integer, optional,  intent(in) :: hdfVersion

      call output ( 'Level 1 product type : ' )
      call output ( trim(l1_type), advance='no')
      if ( present(hdfVersion) ) then
        call blanks(4)
        call output ( 'hdf ' )
        call output ( hdfVersion, advance='yes')
      else
        call output ( ' ', advance='yes')
      end if
      call blanks(15)
      call output ( 'name : ' )
      call blanks(8)
      call output ( trim(Name), advance='yes')
    end subroutine proclaim

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for GlobalSettings = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

    ! --------------------------  CreateDirectTypeFromMLSCFInfo  -----
    function CreateDirectTypeFromMLSCFInfo ( root, DirectFile ) result (Direct)
    use Hdf, only: DFACC_CREATE
    integer, intent(in) :: ROOT         ! Tree node
    type (DirectData_T) :: Direct
    type (MLSFile_T)    :: DirectFile
    ! Local variables
    integer :: SON                      ! Tree node
    integer :: GSON                     ! Tree node
    integer :: I                        ! Loop counters
    integer :: FIELD                    ! Field identifier
    integer :: FILE
    character(len=1024) :: FILENAME     ! Output full filename
    character(len=1024) :: FILE_BASE    ! made up of
    integer :: l2gp_Version
    character(len=1024) :: PATH         ! path/file_base
    integer :: RETURNSTATUS
    logical, parameter :: DEEBUG = .false.
    ! Executable
    call SetupNewDirect(Direct, 0)
    l2gp_Version = 1
    do i = 2, nsons(root)               ! Skip DirectFileName command
      son = subtree ( i, root )
      field = get_field_id ( son )
      if ( nsons(son) == 2 ) gson = subtree(2,son)
      select case ( field )
      case ( f_type )
        Direct%Type = decoration(gson)
        Direct%autoType = Direct%Type
      case ( f_file )
        file = sub_rosa(subtree(2,son))
        Direct%fileIndex = file
      end select
    end do
    RETURNSTATUS = 0
    call get_string ( file, filename, strip=.true. )
    call split_path_name(filename, path, file_base)
    Direct%filenameBase = file_base
    if ( .not. TOOLKIT ) then
      Direct%handle = 0
      Direct%filename = filename
    endif
    if ( Direct%Type ==  l_l2gp ) then
      if ( TOOLKIT ) &
        & Direct%Handle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
        & mlspcf_l2gp_end, &
        & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
        & exactName=Direct%Filename)
      returnStatus = InitializeMLSFile(DirectFile, content = 'l2gp', &
        & name=Direct%Filename, shortName=file_base, &
        & type=l_swath, access=DFACC_CREATE, HDFVersion=HDFVERSION_5, &
        & PCBottom=mlspcf_l2gp_start, PCTop=mlspcf_l2gp_end)
    else if ( Direct%Type ==  l_l2dgg ) then
      if ( TOOLKIT ) &
        & Direct%Handle = GetPCFromRef(file_base, mlspcf_l2dgg_start, &
        & mlspcf_l2dgg_end, &
        & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
        & exactName=Direct%Filename)
      returnStatus = InitializeMLSFile(DirectFile, content = 'l2dgg', &
        & name=Direct%Filename, shortName=file_base, &
        & type=l_swath, access=DFACC_CREATE, HDFVersion=HDFVERSION_5, &
        & PCBottom=mlspcf_l2dgg_start, PCTop=mlspcf_l2dgg_end)
    else if ( Direct%Type ==  l_l2fwm ) then
      if ( TOOLKIT ) &
        & Direct%Handle = GetPCFromRef(file_base, mlspcf_l2fwm_full_start, &
        & mlspcf_l2fwm_full_end, &
        & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
        & exactName=Direct%Filename)
      returnStatus = InitializeMLSFile(DirectFile, content = 'l2fwm', &
        & name=Direct%Filename, shortName=file_base, &
        & type=l_hdf, access=DFACC_CREATE, HDFVersion=HDFVERSION_5, &
        & PCBottom=mlspcf_l2fwm_full_start, PCTop=mlspcf_l2fwm_full_end)
    else
      if ( TOOLKIT ) &
        & Direct%Handle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
        & mlspcf_l2dgm_end, &
        & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
        & exactName=Direct%Filename)
      returnStatus = InitializeMLSFile(DirectFile, content = 'l2dgm', &
        & name=Direct%Filename, shortName=file_base, &
        & type=l_hdf, access=DFACC_CREATE, HDFVersion=HDFVERSION_5, &
        & PCBottom=mlspcf_l2dgm_start, PCTop=mlspcf_l2dgm_end)
    end if
    if ( returnStatus /= 0 ) call MLSMessage ( &
       & MLSMSG_Error, ModuleName, &
       & 'Failed in GetPCFromRef for ' // trim(filename) )

    end function CreateDirectTypeFromMLSCFInfo

  end subroutine SET_GLOBAL_SETTINGS

! =====     Private Procedures     =====================================

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GLOBAL_SETTINGS

! $Log$
! Revision 2.120  2006/09/21 18:48:07  pwagner
! Reduce level of dumps in SIDS version
!
! Revision 2.119  2006/07/24 20:35:21  pwagner
! Fixed bug when stopAfterSection is blank
!
! Revision 2.118  2006/07/21 20:12:48  pwagner
! Can select what section to stop after
!
! Revision 2.117  2006/07/20 23:39:53  vsnyder
! Remove unused declarations and USEs
!
! Revision 2.116  2006/06/12 19:28:52  pwagner
! Fallback to climatology noted only if all of ncep, goes4/5 missing
!
! Revision 2.115  2006/04/21 22:28:45  vsnyder
! Allow FlushPFA in global settings section
!
! Revision 2.114  2006/03/15 23:52:24  pwagner
! Removed InputVersion component from PCF, l2cf
!
! Revision 2.113  2006/02/15 00:02:05  pwagner
! Moved revertOutput call
!
! Revision 2.112  2006/02/10 21:16:43  pwagner
! dumps may go to special dumpfile
!
! Revision 2.111  2006/02/06 22:54:51  pwagner
! Should print warnings, not bomb if l1b files not found
!
! Revision 2.110  2006/01/26 00:35:12  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.109  2006/01/10 23:51:11  pwagner
! Fixed segment fault when hdf4 l1boa file
!
! Revision 2.108  2005/11/12 00:58:10  pwagner
! Fixed bug in reading BO_names attribute from l1boa file
!
! Revision 2.107  2005/11/11 21:46:07  pwagner
! Added reading bright objects from l1boa; removed unused settings
!
! Revision 2.106  2005/09/22 23:37:45  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.105  2005/08/19 23:33:03  pwagner
! FwdMdlDB now dumped when dumping global settings
!
! Revision 2.104  2005/07/21 23:45:33  pwagner
! Removed unused l1b fileinfo fields from l2pcf
!
! Revision 2.103  2005/07/12 17:36:16  pwagner
! Dropped global attribute InputVersion; fills MiscNotes if no dao
!
! Revision 2.102  2005/06/14 20:42:38  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.101  2005/06/04 00:14:53  vsnyder
! Import MLSMSG_Warning
!
! Revision 2.100  2005/06/03 23:58:48  pwagner
! Hope it wont bomb if no l1boa
!
! Revision 2.99  2005/06/03 02:11:14  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of a dummy argument.
!
! Revision 2.98  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.97  2005/05/26 22:34:58  vsnyder
! use if (...) continue to ignore a function result
!
! Revision 2.96  2005/05/02 22:59:59  vsnyder
! Add PFAFile command
!
! Revision 2.95  2005/01/27 21:22:04  vsnyder
! Different interface to ReadPFAData
!
! Revision 2.94  2005/01/12 03:18:22  vsnyder
! Read and write PFAData in HDF5
!
! Revision 2.93  2004/12/31 02:43:05  vsnyder
! Working on read/write PFA database
!
! Revision 2.92  2004/12/14 22:52:38  pwagner
! Changes related to stopping early
!
! Revision 2.91  2004/12/13 20:19:48  vsnyder
! Added MakePFA, PFAData, WritePFA.  Improved error handling.  Removed dumps
! triggered by switches == vgrid2 or vgrid, since the dump command can do it now.
!
! Revision 2.90  2004/10/13 00:52:52  vsnyder
! Get HHMMSS_value from MLSStrings, its new home
!
! Revision 2.89  2004/08/17 23:49:36  pwagner
! Dont pad when getting MAFStartTime
!
! Revision 2.88  2004/08/16 17:14:42  pwagner
! Obtains First,LastMAFCtr global attributes from l1boa file
!
! Revision 2.87  2004/08/04 23:19:58  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.86  2004/07/26 18:22:48  pwagner
! Fixed failure to set returnStatus bug in CreateDirectTypeFromMLSCFInfo
!
! Revision 2.85  2004/06/17 22:49:36  pwagner
! Details for dumping VGrids determined by switch setting
!
! Revision 2.84  2004/06/17 00:39:02  vsnyder
! Make output of 'finished readL1BAttribute...' conditional on 'glo' switch
!
! Revision 2.83  2004/06/08 19:28:07  vsnyder
! Add tGrid, improve error detection and reporting
!
! Revision 2.82  2004/05/29 02:50:49  vsnyder
! Added more dumps
!
! Revision 2.81  2004/05/22 02:31:40  vsnyder
! Hook in PFAData, dump
!
! Revision 2.80  2004/03/24 01:02:55  livesey
! Make LeapSecFile public
!
! Revision 2.79  2004/01/23 01:09:48  pwagner
! Only directwrite files entered in global settings eligible to be auto-sourced
!
! Revision 2.78  2004/01/22 00:54:36  pwagner
! Fixed mistaken impression that direct arg is a pointer
!
! Revision 2.77  2003/12/11 22:59:08  pwagner
! May fill DirectWriteDatabase in global settings
!
! Revision 2.76  2003/11/15 00:46:41  pwagner
! maxfailurespermachine, maxfailuresperchunk no longer configuration settings (see comline opts)
!
! Revision 2.75  2003/09/26 21:55:54  pwagner
! Less likely to complain about missing attributes from hdf4 l1boa files
!
! Revision 2.74  2003/09/12 16:28:49  cvuu
! Read OrbitNumber and OrbitPeriod attributes from L1BOA
!
! Revision 2.73  2003/09/02 18:03:23  pwagner
! Now can reset maxfailuresper chunk, machine from global settings
!
! Revision 2.72  2003/07/15 18:18:06  livesey
! Change of forward model config call.
!
! Revision 2.71  2003/07/07 23:50:04  pwagner
! Now uses saved variable L2pcf from writeMetaData
!
! Revision 2.70  2003/06/09 22:49:34  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.69  2003/05/14 22:04:49  pwagner
! Bad leapsecfile name now prints sensibly & stops
!
! Revision 2.68  2003/05/12 20:57:21  pwagner
! Added L2ParSF spec to allow changing staging file name in global settings
!
! Revision 2.67  2003/05/10 22:20:50  livesey
! Tried to calm down -g1..
!
! Revision 2.66  2003/05/07 23:58:04  pwagner
! outputs trimmed LeapSecFileName
!
! Revision 2.65  2003/04/02 23:54:53  pwagner
! Checks for FILENOTFOUND
!
! Revision 2.64  2003/02/27 21:55:14  pwagner
! Calls FillTAI93Attribute
!
! Revision 2.63  2003/02/01 00:50:14  pwagner
! Picks up GlobalAttributes from settings
!
! Revision 2.62  2002/12/11 22:17:55  pwagner
! Added error checks on hdf version
!
! Revision 2.61  2002/11/13 01:07:40  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.60  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.59  2002/10/03 23:01:19  pwagner
! gets the_hdf_version from l1b..setup
!
! Revision 2.58  2002/09/25 20:09:14  livesey
! New global argument to ConstructForwardModelConfig
!
! Revision 2.57  2002/08/28 22:26:39  pwagner
! Moved LEVEL1_HDFVERSION, ILLEGALL1BRADID, MAXNUML1BRADIDS to MLSL2Options
!
! Revision 2.56  2002/08/20 23:02:17  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.55  2002/05/01 22:03:50  pwagner
! If leapsecfile parameter set, will convert start, end times
!
! Revision 2.54  2002/05/01 16:14:13  pwagner
! Dropped S_LEAPSECFILE
!
! Revision 2.53  2002/04/29 17:42:16  pwagner
! Can convert starttime, endtime to tai w/o pcf
!
! Revision 2.52  2002/02/08 22:52:56  livesey
! Hooked up bin selectors
!
! Revision 2.51  2002/01/24 00:14:26  pwagner
! Proclaims level 1 files as input; comments concerning hdf5 conversion
!
! Revision 2.50  2001/12/16 00:57:26  livesey
! Temporary fix to deal with tai93 time issues for sids
!
! Revision 2.49  2001/12/10 20:22:22  livesey
! Added code for EmpiricalGeometry
!
! Revision 2.48  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.47  2001/10/31 18:36:42  livesey
! Added fGrids stuff
!
! Revision 2.46  2001/10/30 00:37:09  pwagner
! Now can dump more l1b data
!
! Revision 2.45  2001/10/26 23:18:35  pwagner
! Complies with l1b data dump
!
! Revision 2.44  2001/10/25 23:35:52  pwagner
! Responds to global switch with bigger dump than glo switch
!
! Revision 2.43  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.42  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.41  2001/09/17 23:13:09  livesey
! Added instrument stuff to global settings etc
!
! Revision 2.40  2001/07/12 23:28:15  livesey
! Got rid of s_cloudForwardModel
!
! Revision 2.39  2001/07/09 22:53:20  pwagner
! Obeys CloudForwardModel; for now same as ForwardModel
!
! Revision 2.38  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.37  2001/05/30 23:56:23  livesey
! Changed for new L1BData
!
! Revision 2.36  2001/05/30 23:04:40  pwagner
! Gets returnStatus from forwardModelGlobalSetup
!
! Revision 2.35  2001/05/29 23:21:07  livesey
! Now uses ForwardModelSupport, not ForwardModelInterface
!
! Revision 2.34  2001/05/26 00:06:49  livesey
! Added call to DealloteL1BData
!
! Revision 2.33  2001/05/24 20:54:15  pwagner
! Deleted p_ccs..times
!
! Revision 2.32  2001/05/24 20:36:13  pwagner
! Warns if glob. stg. overrides pcf
!
! Revision 2.31  2001/05/17 00:29:03  pwagner
! Works without toolkit, PCF at last
!
! Revision 2.30  2001/05/15 23:46:32  pwagner
! Now optionally uses hhmmss_value
!
! Revision 2.29  2001/05/14 23:45:08  pwagner
! Start, end times now added to MAF offsets from L1BOA
!
! Revision 2.28  2001/05/11 23:44:43  pwagner
! Better dump; uses strip=TRUE
!
! Revision 2.27  2001/05/11 01:56:17  vsnyder
! Move the getting of sub_rosa_index
!
! Revision 2.26  2001/05/11 00:09:05  pwagner
! Gets p_.. from init_tables; unquotes strings
!
! Revision 2.25  2001/05/10 18:26:22  pwagner
! Improved dump_global_settings
!
! Revision 2.24  2001/05/09 23:35:04  pwagner
! Added dump_global_settings
!
! Revision 2.23  2001/05/04 17:15:36  pwagner
! Many added settings, esp. L1B files, so level2 can run w/o PCF
!
! Revision 2.22  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.21  2001/04/26 02:52:17  vsnyder
! Fix up CVS stuff
!
! Revision 2.20  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.19  2001/04/26 00:07:33  livesey
! Stuff to support reading of l2pc files
!
! Revision 2.18  2001/04/24 00:31:42  vsnyder
! Finish adding 'time' command
!
! Revision 2.17  2001/04/23 23:48:41  vsnyder
! Finish adding 'time' command
!
! Revision 2.16  2001/04/23 23:42:00  vsnyder
! Add 'time' command
!
! Revision 2.15  2001/04/21 01:25:54  livesey
! Now passes l2gpdatabase to more people who need it.
!
! Revision 2.14  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.13  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.12  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.11  2001/03/28 22:00:33  livesey
! Interim version, now handles vGrids as part of forwardModelConfig
!
! Revision 2.10  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.9  2001/03/17 03:24:23  vsnyder
! Work on forwardModelGlobalSetup
!
! Revision 2.8  2001/03/17 00:57:36  livesey
! Removed dump.
!
! Revision 2.7  2001/03/17 00:45:38  livesey
! Added ForwardModelConfigDatabase
!
! Revision 2.6  2001/03/09 02:30:13  vsnyder
! Allocate correct size for FMI and TFMI
!
! Revision 2.5  2001/03/09 00:24:30  vsnyder
! Do subscripts right
!
! Revision 2.4  2001/03/08 03:23:09  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.3  2001/03/07 22:46:05  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.2  2000/11/16 01:53:40  vsnyder
! Take timing back out.  Don't do it in sections that are only parameter settings.
!
! Revision 2.1  2000/11/16 01:45:25  vsnyder
! Implement timing.
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
