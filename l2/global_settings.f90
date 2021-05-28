! C/pyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Global_Settings

  use MLSCommon, only: FileNameLen, MLSFile_T, NameLen, Tai93_Range_T
  use HighOutput, only: AddRow, AddRow_Divider, AddRow_Header, &
    & BeVerbose, OutputCalendar, OutputNamedValue, OutputTable, StartTable, &
    & StyledOutput
  use MLSL2Options, only: CheckPaths, L2CFNode, Level1_HDFVersion, &
    & Need_L1BFiles, SpecialDumpFile, StopAfterSection, Toolkit, &
    & MLSL2Message
  use Output_M, only: Blanks, NewLine, Output, &
    & RevertOutput, SwitchOutput

  implicit none

  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! Leapsecfilename                 For time conversions (if avoiding PCF)
! BrightObjects                   As determined by l1b

!     (subroutines and functions)
! Set_Global_Settings             Get global settings from l2cf
! L1MAFToL2Profile                Find profile number closest to a MAF
! L2ProfileToL1MAF                Find MAF closest to a profile number
! === (end of toc) ===

  public :: L1maftol2profile, L2Profiletol1maf, Set_Global_Settings

  character(len=FileNameLen), public :: LEAPSECFILENAME = ''

  ! These must be values consistent with level 1
  ! (Some canny coding below could make this more robust)
  integer, parameter :: BO_NAMEDIMS = 14
  integer, parameter :: BO_NAMELEN  = 14

  ! This next should be large enough to hold the entire list of BO names
  integer, parameter :: BONAMELISTLEN = 256
  character(len=BONAMELISTLEN), public, save :: brightObjects = &
    & 'MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO, ' // &
    & 'MOON, SUN, SOLAR-BARYCNTR, EARTH-BARYCNTR, GALACTICCENTER'

  ! Do we check for the correct month's l2pc files during checkPaths preflight?
  logical, parameter :: CHECKL2PCMONTHCORRECT = .true.
  ! Do we output this month's calendar
  logical, parameter :: OUTPUTTHISMONTHSCAL = .true.
  ! Consult dates_module for leap seconds
  logical, parameter :: LEAPSINDATESMODULE  = .true.
  ! Can tai93 times be < 0?
  logical, parameter :: NEGATIVETAI93OK     = .true.
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, private              :: ERROR, WARNING
  logical, parameter            :: countEmpty = .true.
  logical, parameter            :: DEEBUG = .false.

contains

  ! Given a MAF, returns the closest profile
  ! Return -1 if something goes wrong
  ! Remember--MAFs start at 0, not 1 (?)
  !
  ! Does anyone really think that having MAFs start at 0 wasn't a mistake?
  ! It's probably too late to reverse that decision, but
  ! because profiles, MIFs, and every other geolocation index start at 1
  ! "singling" out MAFs and treating them in a manner inconsistent with
  ! with normal expectations is a recipe for .. mistakes.
  ! **
  ! Here's a good software design rule:
  ! Given a choice between 2 designs, where one design is more likely
  ! to lead to mistakes and the other is less likely to lead to mistakes
  !    ===>   choose the design less likely to lead to mistakes    <===
  ! **
  !
  function L1MAFToL2Profile ( MAF, fileDatabase, &
    & MIF, Debugging, Phi, Lat ) result( profile )
    use L1BData, only: L1BData_t, &
      & AssembleL1BQtyName, DeallocateL1BData, &
      & ReadL1BData 
    use L2GPData, only: L2GPData_t, L2GPNameLen, MaxSwathNamesBufSize, &
      & ReadL2GPData, DestroyL2GPContents
    use MLSFiles, only: HDFVersion_5, &
      & MLS_InqSwath, GetMLSFileByType
    use MLSKinds, only: R4, R8
    use MLSMessageModule, only: MLSMSG_Warning
    use MLSNumerics, only: ClosestElement
    use MLSStringLists, only: GetStringElement
    ! Args
    integer, intent(in)                 :: MAF
    type (MLSFile_T), pointer           :: FILEDATABASE(:)
    integer, optional, intent(in)       :: MIF ! Not for MAFStartTimeTAI
    logical, optional, intent(in)       :: Debugging
    real(r4), optional, intent(out)     :: Phi
    real(r4), optional, intent(out)     :: Lat
    integer :: profile
    ! Internal variables
    integer, dimension(1) :: indices
    integer :: L1BFLAG
    character(len=namelen) :: l1bItemName
    type (L1BData_T) :: l1bField   ! L1B data
    type(MLSFile_T), pointer :: L1BFile, L2GPFile
    type (L2GPData_T) :: l2gp
    integer :: myMIF ! To be used in diagnostic dumps
    integer :: LISTSIZE, NOMAFS, NOSWATHS
    character (len=namelen) :: QUANTITY
    character (len=L2GPNameLen) :: swath
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    ! logical, parameter          :: DEEBug = .true.
    logical                     :: DEEBug
    logical, parameter          :: DumpLats = .true.
    ! Executable
    DeeBug = .false.
    if ( present(Debugging) ) DeeBug = Debugging
    profile = -1
    myMIF = 1
    if ( present(MIF) ) myMIF = MIF
    ! 1st--read the MAF times (or geodAngle if that's what you're matching)
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName // 'L1MAFToL2Profile', &                      
      & 'l1boa File not found--hope you dont need one' )
      return
    endif
    if ( myMIF == 1 ) then
      quantity = 'MAFStartTimeTAI'
    else
      quantity = 'GHz/GeodAngle'
    endif
    l1bItemName = AssembleL1BQtyName ( quantity, HDFVERSION_5, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
      & l1bFlag, dontPad=.true.)
    ! 2nd--read the profile
    L2GPFile => GetMLSFileByType ( filedatabase, content='l2gp' )
    if ( .not. associated(L2GPFile) ) &
      & L2GPFile => GetMLSFileByType ( filedatabase, content='l2dgg' )
    noSwaths = mls_InqSwath ( L2GPFile%name, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5 )
    call GetStringElement (trim(swathList), swath, 1, countEmpty )
    call ReadL2GPData ( L2GPFile, trim(swath), l2gp )
    if ( DEEBUG ) call outputNamedValue( 'shape(l2gp%time)', shape(l2gp%time) )
    if ( DEEBUG ) call outputNamedValue( 'shape(l1bField%dpField)', shape(l1bField%dpField) )
    if ( DEEBUG ) call outputNamedValue( 'lbound(l1bField%dpField(1,1,:))', lbound(l1bField%dpField(1,1,:)) )
    if ( myMIF == 1 ) then
      call ClosestElement( l1bField%dpField(1,myMIF,MAF+1)*1._r8, &
        & l2gp%time, indices )
    else
      call ClosestElement( l1bField%dpField(1,myMIF,MAF+1), &
        & l2gp%geodAngle*1._r8, indices )
      call OutputNamedValue( 'phi(MAF)', l1bField%dpField(1,myMIF,MAF+1) )
      call OutputNamedValue( 'phi(nearestprofile)', l2gp%geodAngle(indices(1)) )
    endif
    profile = indices(1)
    if ( DEEBUG ) then
      call outputNamedValue( 'l2gp profile', profile )
      call outputNamedValue( 'l1b MAF', MAF )
      call output ( '-----------', advance='yes' )
      call output ( 'MAFStartTimesTAI', advance='yes' )
      call outputNamedValue( 'l2gpvalue', l2gp%time(profile) )
      call outputNamedValue( 'l1bvalue', l1bField%dpField(1,1,MAF) )
      call outputNamedValue( 'l1bvalue[+1]', l1bField%dpField(1,1,MAF+1) )
      call output ( '-----------', advance='yes' )
    endif
    call deallocatel1bdata( l1bfield )
    if ( (DEEBug .and. DumpLats) .or. present(Phi) ) then
      call ReadL1BData ( L1BFile, '/GHz/GeodAngle', l1bField, noMAFs, &
        & l1bFlag, dontPad=.true.)
      if ( DeeBug ) then
      call output ( '-----------', advance='yes' )
      call output ( 'Geod Angles', advance='yes' )
      call outputNamedValue( 'shape(l1b)', shape(l1bField%dpField) )
      call outputNamedValue( 'l2gpvalues', &
        & (/ l2gp%GeodAngle(profile) /) )
      call outputNamedValue( 'MAFs', &
        & (/ MAF, MAF+1, MAF+2 /) )
      call outputNamedValue( 'l1bvalues', &
        & (/ l1bField%dpField(1,myMIF,MAF-1), l1bField%dpField(1,myMIF,MAF), &
        & l1bField%dpField(1,myMIF,MAF+1), l1bField%dpField(1,myMIF,MAF+2) /)  )
      call outputNamedValue( 'their MAFs', &
        & (/ MAF-1, MAF, &
        & MAF+1, MAF+2 /)  )
      call output ( '-----------', advance='yes' )
      endif
      if ( present(Phi) ) Phi = l1bField%dpField(1,myMIF,MAF)

      call deallocatel1bdata( l1bfield )
      
      call ReadL1BData ( L1BFile, '/GHz/GeodLat', l1bField, noMAFs, &
        & l1bFlag, dontPad=.true.)
      if ( DeeBug ) then
      call output ( '-----------', advance='yes' )
      call output ( 'Geod Lats', advance='yes' )
      call outputNamedValue( 'shape(l1b)', shape(l1bField%dpField) )
      call outputNamedValue( 'l2gpvalues', &
        & (/ l2gp%Latitude(profile) /) )
      call outputNamedValue( 'l1bvalues', &
        & (/ l1bField%dpField(1,myMIF,MAF-1), l1bField%dpField(1,myMIF,MAF), &
        & l1bField%dpField(1,myMIF,MAF+1), l1bField%dpField(1,myMIF,MAF+2) /)  )
      call outputNamedValue( 'their MAFs', &
        & (/ MAF-1, MAF, &
        & MAF+1, MAF+2 /)  )
      call output ( '-----------', advance='yes' )
      endif
      if ( present(Lat) ) Lat = l1bField%dpField(1,myMIF,MAF)
    endif
    call destroyl2gpcontents( l2gp )
  end function L1MAFToL2Profile

  ! Given a profile, returns the closest MAF
  ! Return -1 if something goes wrong
  ! Remember--MAFs start at 0, not 1
  function L2ProfileToL1MAF ( profile, fileDatabase, MIF ) result( MAF )
    use Dump_0, only: Dump
    use L1BData, only: L1BData_T, NameLen, &
      & AssembleL1BQtyName, DeallocateL1BData, &
      & ReadL1BData 
    use L2GPData, only: L2GPData_t, L2GPNameLen, MaxSwathNamesBufSize, &
      & ReadL2GPData, DestroyL2GPContents
    use MLSFiles, only: HDFVersion_5, &
      & MLS_InqSwath, GetMLSFileByType
    use MLSKinds, only: R8
    use MLSMessageModule, only: MLSMSG_Warning
    use MLSNumerics, only: ClosestElement
    use MLSStringLists, only: GetStringElement
    ! Args
    integer, intent(in)                 :: profile
    type (MLSFile_T), pointer           :: FILEDATABASE(:)
    integer, intent(in), optional       :: MIF
    integer :: MAF
    ! Internal variables
    integer, dimension(1) :: indices
    integer :: L1BFLAG
    character(len=namelen) :: l1bItemName
    type (L1BData_T) :: l1bField   ! L1B data
    type(MLSFile_T), pointer :: L1BFile, L2GPFile
    type (L2GPData_T) :: l2gp
    integer :: myMIF ! To be used in diagnostic dumps
    integer :: LISTSIZE, NOMAFS, NOSWATHS
    character (len=namelen) :: QUANTITY
    character (len=L2GPNameLen) :: swath
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    ! Executable
    MAF = -1
    myMIF = 1
    if ( present(MIF) ) myMIF = MIF
    ! 1st--read the MAF times (or geodAngle if that's what you're matching)
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName // 'L2ProfileToL1MAF', &                      
      & 'l1boa File not found--hope you dont need one' )
      return
    end if
    if ( myMIF == 1 ) then
      quantity = 'MAFStartTimeTAI'
    else
      quantity = 'GHz/GeodAngle'
    endif
    l1bItemName = AssembleL1BQtyName ( quantity, HDFVERSION_5, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
      & l1bFlag, dontPad=.true.)
    ! 2nd--read the profile
    L2GPFile => GetMLSFileByType ( filedatabase, content='l2gp' )
    if ( .not. associated(L2GPFile) ) &
      & L2GPFile => GetMLSFileByType ( filedatabase, content='l2dgg' )
    noSwaths = mls_InqSwath ( L2GPFile%name, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5 )
    if ( DEEBUG ) call outputNamedValue( 'profile', profile )
    if ( DEEBUG ) call outputNamedValue( 'noSwaths', noSwaths )
    call GetStringElement (trim(swathList), swath, 1, countEmpty )
    if ( DEEBUG ) call outputNamedValue( 'swath', swath )
    call ReadL2GPData ( L2GPFile, trim(swath), l2gp )
    if ( DEEBUG ) call dump( l2gp%time, 'l2gp%time' )
    if ( l2gp%time(profile) < 0.d0 ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName // 'L2ProfileToL1MAF', &                      
      & 'time for this profile < 0; did the chunk crash?' )
      call deallocatel1bdata( l1bfield )
      call destroyl2gpcontents( l2gp )
      return
    end if
    if ( myMIF == 1 ) then
      call ClosestElement( l2gp%time(profile), &
        & l1bField%dpField(1,myMIF,:)*1._r8, indices )
    else
      call ClosestElement( l2gp%geodAngle(profile)*1._r8, &
        & l1bField%dpField(1,myMIF,:), indices )
      call OutputNamedValue( 'phi(profile)', l2gp%geodAngle(profile) )
      call OutputNamedValue( 'phi(nearest MAF)', l1bField%dpField(1,myMIF,indices(1)) )
    endif
    MAF = indices(1) - 1
    if ( DEEBUG ) call outputNamedValue( 'size(time)', size(l2gp%time) )
    if ( DEEBUG ) call outputNamedValue( 'time(l2gp)', l2gp%time(profile) )
    if ( DEEBUG .and. MAF > 0 ) call outputNamedValue( 'l1bvalue', l1bField%dpField(1,1,MAF) )
    if ( DEEBUG ) call outputNamedValue( 'l1bvalue[+1]', l1bField%dpField(1,1,MAF+1) )
    call deallocatel1bdata( l1bfield )
    call destroyl2gpcontents( l2gp )
  end function L2ProfileToL1MAF

  subroutine Set_Global_Settings ( root, forwardModelConfigDatabase, &
    & fileDatabase, FGrids, L2GPDatabase, directDatabase, processingRange )

    use BitStuff, only: IsBitSet
    use Dates_Module, only: IsUTCInRange, PrecedesUTC, ResetStartingDate, &
      & SecondsBetween2UTCs, Utc2tai93s, UTC_To_Yyyymmdd, YyyyDoy_To_Mmdd
    use Declaration_Table, only: Named_Value, Redeclare, Str_Value
    use DirectWrite_M, only: DirectData_T, &
      & AddDirectToDatabase, Dump, SetupNewDirect
    use DumpCommand_M, only: DumpCommand
    use Dump_1, only: Dump
    use EmpiricalGeometry, only: InitEmpiricalGeometry
    use FGrid, only: AddFGridToDatabase, CreateFGridFromMLSCFInfo, Dump, Fgrid_T
    use ForwardModelConfig, only: AddForwardModelConfigToDatabase, Dump, &
      & ForwardModelConfig_T
    use ForwardModelSupport, only: ConstructForwardModelConfig, &
      & ForwardModelGlobalsetup, CreateBinSelectorFromMLSCFInfo
    use HDF, only: Dfacc_Create
    use IGRF_Int, only: Read_Gh
    use Init_Tables_Module, only: F_File, F_Reset, F_Type, &
      & L_L2GP, L_L2dgg, L_L2fwm, &
      & Parm_Indices, &
      & First_Parm, Last_Parm, P_Brightobjects, &
      & P_Cycle, P_Endtime, P_Igrf_File, P_Instrument, &
      & P_LeapsecFile, P_Output_Version_String, P_PFAFile, P_Starttime, &
      & S_Binselector, S_DirectwriteFile, S_Dump, S_Empiricalgeometry, &
      & S_Fgrid, S_FlushPFA, S_ForwardModel, S_ForwardModelglobal, &
      & S_L1boa, S_L1brad, S_L2parsf, S_MakePFA, S_PFAData, S_ReadPFA, &
      & S_Tgrid, S_Time, S_Vgrid, S_WritePFA
    use Intrinsic, only: L_HDF, L_Swath, Spec_Indices
    use L1BData, only: L1BData_T, NameLen, &
      & AssembleL1BQtyName, DeallocateL1BData, Dump, FindMaxMaf, &
      & L1BRadsetup, L1BOASetup, ReadL1BAttribute, ReadL1BData
    use L2GPData, only: L2GPData_T
    use L2PC_M, only: AddbinselectortoDatabase, Binselectors
    use MLSFiles, only: Filenotfound, HDFVersion_5, &
      & AddFiletoDatabase, GetPCFromref, GetMLSFileByName, GetMLSFileByType, &
      & InitializeMLSFile, MLS_CloseFile, MLS_OpenFile, Split_Path_Name
    use MLSHDF5, only: GetAllHDF5GroupNames
    use MLSKinds, only: R8
    use MLSL2Timings, only: Section_Times
    use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning
    use MLSPCF2, only: MLSPCF_L2GP_Start, MLSPCF_L2GP_End, &
      & MLSPCF_L2dgm_Start, MLSPCF_L2dgm_End, MLSPCF_L2fwm_Full_Start, &
      & MLSPCF_L2fwm_Full_End, &
      & MLSPCF_L2dgg_Start, MLSPCF_L2dgg_End
    use MLSStrings, only: Hhmmss_Value, LowerCase, Trim_Safe
    use MLSStringLists, only: Array2List, CatLists, SwitchDetail, &
      & NumStringElements, StringElement, StringElementNum
    use MLSSignals_M, only: Instrument, Modules, Dump_Modules, GetModuleName
    use MoreTree, only: Get_Boolean, Get_Field_Id, Get_Label_And_Spec, &
      & Get_Spec_Id, StartErrorMessage
    use Next_Tree_Node_M, only: Next_Tree_Node, Next_Tree_Node_State
    use PFAData_M, only: Get_PFAData_From_L2cf, Flush_PFAData, Make_PFAData, &
      & Read_PFAData, Write_PFAData
    use PFADatabase_M, only: Process_PFA_File
    use PCFHDR, only: GlobalAttributes, FillTAI93Attribute
    use ReadAPriori, only: APrioriFiles
    use SDPToolkit, only: Max_Orbits, MLS_UTCToTAI, &
      & PGSD_Dem_30arc, PGSD_Dem_90arc, &
      & PGSD_Dem_Elev, PGSD_Dem_Water_Land, &
      & PGS_Dem_Open, PGS_S_Success
    use String_Table, only: Display_String, Get_String
    use Time_M, only: SayTime, Time_Now
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decorate, Decoration, Node_Id, Nsons, Sub_Rosa, Subtree, &
      & Dump_Tree_Node
    use Tree_Types, only: N_Equal
    use VGrid, only: CreateVGridFromMLSCFInfo
    use VGridsDatabase, only: AddVGridToDatabase, VGrids
    use WriteMetaData, only: L2PCF

    ! placed non-alphabetically due to Lahey internal compiler error
    ! (How much longer must we endure these onerous work-arounds?)
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributeInFile
    use HDF5, only: H5GClose_f, H5GOpen_f

    integer, intent(in)                    :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(ForwardModelConfig_T), pointer    :: FORWARDMODELCONFIGDATABASE(:)
    type (MLSFile_T), pointer              :: FILEDATABASE(:)
    type ( fGrid_T ), pointer, optional    :: FGRIDS(:)
    type ( l2gpData_T), pointer, optional  :: L2GPDATABASE(:)
    type (DirectData_T), pointer, optional :: DIRECTDATABASE(:)
    type (TAI93_Range_T), optional         :: PROCESSINGRANGE ! Data processing range

    ! Local variables
    character(len=BO_NAMELEN), dimension(BO_NAMEDIMS) :: BO_names
    type (l1bData_T)                      :: BO_stat
    character(len=BONAMELISTLEN) :: BO_today
    logical, parameter :: DEEBUG = .false.
    integer :: DetailReduction
    integer :: Details             ! How much info about l1b files to dump
    type (MLSFile_T) :: DirectFile
    character(len=NameLen) :: End_time_string, Start_time_string
    character(len=FileNameLen) :: FilenameString
    logical :: GOT(FIRST_PARM:LAST_PARM)
    integer :: Gson
    integer :: Field
    integer :: I, J                ! Index of son, grandson of root
    logical ::  ItExists
    character(len=namelen) :: itsname
    integer :: Key
    type (L1BData_T) :: L1bField   ! L1B data
    type(MLSFile_T), pointer :: L1BFile
    character(len=namelen) :: L1bItemName
    integer :: L1BFLAG
    integer :: Me = -1             ! String index for trace
    integer :: MM, DD
    real(r8) :: MINTIME, MAXTIME   ! Time Span in L1B file data
    character(len=namelen) :: modulenames
    integer :: NAME                ! Sub-rosa index of name of vGrid or hGrid
    character(len=NameLen) :: Name_string
    integer :: NOMAFS              ! Number of MAFs of L1B data read
    integer :: NUMFILES
    integer :: OrbNum(max_orbits)
    real(r8) :: OrbPeriod(max_orbits)
    integer :: OUTPUT_VERSION_STRING ! Sub_rosa index
    integer :: Param               ! Tree index of param i.e., name before =
    integer :: Param_id            ! e.g., p_brightObjects
    character (len=namelen) :: QUANTITY
    logical :: Reset
    logical :: Restricted          ! Some commands not available
    integer :: ReturnStatus        ! non-zero means trouble
    integer :: SON                 ! Son of root
    integer :: spec_id             ! e.g., s_binSelector
    type(next_tree_node_state) :: State ! of tree traverser
    logical :: StopEarly
    integer :: Sub_rosa_index
    integer :: The_HDF_version     ! 4 or 5 (corresp. to hdf4 or hdf5)
    logical :: TIMING              ! For S_Time
    logical :: StartTimeIsAbsolute, stopTimeIsAbsolute
    integer :: Status
    real :: T1                 ! For S_Time
    real(r8) :: Start_time_from_1stMAF, End_time_from_1stMAF
    logical :: verbose
    logical :: verboser
    logical :: wasAlreadyOpen

    integer, parameter :: Param_restricted = 1 ! Parameter not allowed
    integer, parameter :: Spec_restricted = param_restricted + 1 ! Spec not allowed
    ! DEM nonsense
    integer, dimension(2) :: resolutionList
    integer   :: numResolutions
    integer, dimension(2) :: layerList
    integer   :: numLayers
    ! Executable

    restricted = .not. all( (/present(FGrids), present(l2gpDatabase), &
                           present(DirectDatabase), present(processingRange) /))
    timing = section_times
    if ( timing ) call time_now ( t1 )
    stopEarly = ( &
      & index('global_setting,chunk_divide', &
      & lowercase( trim(stopAfterSection) ) ) > 0 )
    stopEarly = ( stopearly .and. stopAfterSection /= ' ' )
    error = 0
    warning = 0
    got = .false.
    startTimeIsAbsolute = .false.
    stopTimeIsAbsolute = .false.
    LeapSecFileName = ''

    call trace_begin ( me, 'SET_GLOBAL_SETTINGS', root, cond=toggle(gen) )

    Details = switchDetail(switches, 'glo') - 3

    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    elseif ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    end if
    verbose = BeVerbose( 'glo', -1 )
    verboser = BeVerbose( 'glo', 0 )
    
    ! Start Digital Elevation Model
    if ( Toolkit ) then
      ! Initialize
      numResolutions = 2
      numLayers = 2
      resolutionList(1) = PGSd_DEM_30ARC
      resolutionList(2) = PGSd_DEM_90ARC
      layerList(1) = PGSd_DEM_ELEV
      layerList(2) = PGSd_DEM_WATER_LAND
      status = PGS_DEM_Open ( resolutionList, numResolutions, &
        & layerList, numLayers )
      if ( verbose .or. status /= PGS_S_SUCCESS ) &
        & call outputNamedValue( 'PGS_DEM_Open status', status )
      if ( status /= PGS_S_SUCCESS )  call announce_error( 0, &
        & '*** Error: Unable to open Digital Elevation Model (DEM) files;' // &
        & ' missing from PCF? ***', just_a_warning = .false.)
    end if

    do
      son = next_tree_node ( root, state )
      call get_label_and_spec ( son, name, key )
      if ( son == 0 ) exit
      L2CFNODE = son
      reset = .false.
      do j = 2, nsons(key) ! fields of the "time" specification
        gson = subtree(j, key)
        if ( gson == 0 ) exit
        if ( node_id(son) == n_equal ) cycle ! Why can't we loop over params ?
        field = get_field_id(gson)   ! tree_checker prevents duplicates
        if (nsons(gson) > 1 ) gson = subtree(2,gson) ! Gson is value
        select case ( field )
        case ( f_reset )
          reset = Get_Boolean ( gson )
        case default
          ! Shouldn't get here if the type checker worked
        end select
      end do ! j = 2, nsons(key)

      if ( node_id(son) == n_equal ) then
        sub_rosa_index = sub_rosa(subtree(2,son))
        param = subtree(1,son)
        call redeclare ( sub_rosa(param), 0.0d0+sub_rosa_index, named_value, &
          & str_value, son ) ! Put its value in the declaration table
        param_id = decoration(param)
        if ( TOOLKIT .and. &
          & any( param_id == &
          & (/ p_output_version_string, p_cycle, p_starttime, p_endtime, &
          & p_leapsecfile /) ) ) then
          if ( verboser ) call announce_error(0, &
            & '*** l2cf parameter global setting ignored ***', &
            & just_a_warning = .true.)
          cycle
        end if
        if ( param_id < first_parm .or. param_id > last_parm ) cycle ! how?
        got(param_id) = .true.
        if ( DEEBUG ) call display_string( parm_indices(param_id), advance='yes' )
        select case ( param_id )
        ! This will allow us to use different names from the toolkit
        ! (Now why would you want to do that?)
        case ( p_brightObjects )
          call get_string ( sub_rosa_index, brightObjects, strip=.true. )
        case ( p_cycle )
          call get_string ( sub_rosa_index, l2pcf%cycle, strip=.true. )
        case ( p_endtime )
          if ( restricted ) call NotAllowed ( son, param_restricted )
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
        case ( p_IGRF_file )
          if ( restricted ) call NotAllowed ( son, param_restricted )
          call read_gh ( son )
        case ( p_instrument )
          instrument = decoration(subtree(2,son))
        case ( p_leapsecfile )
          if ( restricted ) call NotAllowed ( son, param_restricted )
          call get_string ( sub_rosa_index, LeapSecFileName, strip=.true. )
          inquire(file=trim(LeapSecFileName), exist=itExists)
          if ( .not. itExists ) then
            call announce_error(son, &
            & '*** Leap Second File ' // trim(LeapSecFileName) // &
            & ' not found', &
            & just_a_warning = .false.)
            call MLSL2Message ( MLSMSG_Error, ModuleName, &                      
            & '(Please check file name and path)' )    
          end if
        case ( p_output_version_string )
          output_version_string = sub_rosa_index
          call get_string ( output_version_string, l2pcf%PGEVersion, strip=.true. )
        case ( p_PFAFile )
          if ( restricted ) call NotAllowed ( son, param_restricted )
          do j = 2, nsons(son)
            if ( process_PFA_File ( sub_rosa(subtree(j,son)), &
              & subtree(j,son) ) /= 0 ) continue
          end do
        case ( p_starttime )
          if ( restricted ) call NotAllowed ( son, param_restricted )
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
        case default
          call announce_error(son, 'unrecognized global settings parameter')
        end select
      else
        call get_label_and_spec ( son, name )
        L2CFNODE = son
        spec_id = get_spec_id(son)
        if ( TOOLKIT .and. &
          & any( spec_id == &
          & (/ s_l1boa, s_l1brad, s_l2parsf /) ) ) then
          if ( verboser ) call announce_error(0, &
            & '*** l2cf spec global setting ignored ***', &
            & just_a_warning = .true.)
          cycle
        endif
        if ( DEEBUG ) call display_string( spec_indices(spec_id), advance='yes' )
        select case ( spec_id )
        case ( s_binSelector )
          call decorate (son, AddBinSelectorToDatabase ( &
            & binSelectors, CreateBinSelectorFromMLSCFInfo ( son ) ) )
        case ( s_directWriteFile )
          if ( restricted ) call notAllowed ( son, spec_restricted )
          call decorate (son, AddDirectToDatabase ( &
            & DirectDatabase, &
            & CreateDirectTypeFromMLSCFInfo ( son, DirectFile ) ) )
          numFiles = AddFileToDataBase(fileDataBase, DirectFile)
        case ( s_dump )
          if ( error == 0 ) then
            call dumpCommand ( son, &
              & forwardModelConfigs=forwardModelConfigDatabase, &
              & fileDatabase=fileDatabase )
          else
            call announce_error ( subtree(1,son), &
              & 'Preceeding errors prevent doing a dump here.' )
          end if
        case ( s_empiricalGeometry )
          if ( restricted ) call NotAllowed ( son, param_restricted )
          call InitEmpiricalGeometry ( son )
        case ( s_flushPFA )
          call flush_PFAData ( son, status )
          error = max(error,status)
        case ( s_fgrid )
          if ( restricted ) call notAllowed ( son, spec_restricted )
          call decorate ( son, AddFGridToDatabase ( fGrids, &
            & CreateFGridFromMLSCFInfo ( name, son ) ) )
          if ( switchDetail(switches, 'fgrid') > -1 ) &
            & call dump( fgrids(size(fGrids) ) )
        case ( s_forwardModelGlobal ) !??? Begin temporary stuff for l2load
          if ( restricted ) call notAllowed ( son, spec_restricted )
          if ( stopEarly .or. &
            & ( checkPaths .and. .not. CHECKL2PCMONTHCORRECT ) ) cycle
          call forwardModelGlobalSetup ( son, returnStatus, fileDataBase )
          error = max(error, returnStatus)
        case ( s_forwardModel )
          if ( .not. stopEarly ) then
            call decorate (son, AddForwardModelConfigToDatabase ( &
            & forwardModelConfigDatabase, &
            & ConstructForwardModelConfig ( name, son, .true. ) ) )
          end if
        case ( s_l1boa )
          if ( restricted ) call notAllowed ( son, spec_restricted )
          if ( .not. NEED_L1BFILES ) then
            call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
            & 'l1boa File not needed -- and so ignored' )
            cycle
          endif
          the_hdf_version = LEVEL1_HDFVERSION
          call l1boaSetup ( son, filedatabase, F_FILE, hdfVersion=the_hdf_version )
          if( switchDetail(switches, 'pro') > -1 ) then                            
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
          if ( restricted ) call notAllowed ( son, spec_restricted )
          if ( .not. NEED_L1BFILES ) then
            call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
            & 'l1brad File not needed -- and so ignored' )
            cycle
          endif
          the_hdf_version = LEVEL1_HDFVERSION
          call l1bradSetup ( son, filedatabase, F_FILE, &
            & hdfVersion=the_hdf_version )
          sub_rosa_index = sub_rosa(subtree(2,subtree(2, son)))
          call get_string ( sub_rosa_index, FilenameString, strip=.true. )
          if( switchDetail(switches, 'pro') > -1 ) then                            
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
          call MLSL2Message ( MLSMSG_Warning, ModuleName, &
            & 'This version does not use staging file for slave Join commands' )
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
          if ( timing .and. .not. reset ) then
            call sayTime ( 'global_settings', t1=t1, cumulative=.false. )
          else
            call time_now ( t1 )
            timing = .true.
          end if
          call output ( ' Finished time command in global settings', advance='yes' )
        case ( s_tGrid, s_vGrid )
          call decorate ( son, AddVGridToDatabase ( vGrids, &
            & CreateVGridFromMLSCFInfo ( name, son, l2gpDatabase, returnStatus ) ) )
          error = max(error, returnStatus)
        case default
          call announce_error(son, 'unrecognized global settings spec')
        end select
      end if
    end do

    if ( DEEBUG ) call output( 'done with statements', advance='yes' )

    ! Any time conversions?
    if ( restricted ) then
      ! Not if we're a callable forward model server
      call FinishUp
      return
    elseif( .not. NEED_L1BFILES ) then
      ! w/o an l1boa file we'll trust the user to have supplied start, end times
      ! (So we check that she did)
      if ( .not. ( &
        & got(p_starttime) .and. got(p_endtime) ) &
        & ) then
        error = 1
        call MLSL2Message( MLSMSG_Warning, ModuleName, &
          & 'start, end times must be supplied if there is no l1boa file' )
      elseif ( LeapSecFileName /= '' ) then
        call ToolkitTimeConversion
      elseif ( LEAPSINDATESMODULE ) then
        ! Without a leapsec file, let's use date_module's built-in feature
        call datesModuleTimeConversion
      endif
      call FinishUp
      return
    endif

    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
      & 'l1boa File not found--hope you dont need one' )
      call FinishUp
      return
    elseif ( L1BFile%hdfVersion /= hdfversion_5 ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
      & 'l1boa File is the older hdf4' )
    else
      ! Check on module names--do they agree with group names in L1BOA file?
      ! If not, overwrite them
      call MLS_OpenFile ( L1BFile )
      call GetAllHDF5GroupNames ( L1BFile%FileID%f_id, moduleNames )
      if ( verboser ) call outputNamedValue ( 'group names', trim(moduleNames) )
      call mls_closeFile ( L1BFile )
      do i=1, size(modules)
        call GetModuleName ( i, itsName )
        if ( verboser ) call outputNamedValue ( 'module name', trim(itsname) )
        j = StringElementNum( lowercase(moduleNames), lowercase(itsName), countEmpty )
        if ( verboser ) call outputNamedValue ( 'element num', j )
        if ( j > 0 ) then
          if ( itsName /= &
            & StringElement( moduleNames, j, countEmpty ) ) &
            & modules(i)%nameString = &
            & StringElement( moduleNames, j, countEmpty )
        else
            modules(i)%nameString = StringElement( moduleNames, i, countEmpty )
        endif
      enddo
      if ( verbose ) then
        call Dump( modulenames, 'module names' )
        call Dump_Modules
      endif
    endif

    the_hdf_version = &
      & L1BFile%HDFVersion
    if ( the_hdf_version == FILENOTFOUND ) then                                          
      call MLSL2Message ( MLSMSG_Error, ModuleName, &                      
      & 'File not found; make sure the name and path are correct' &
      & // trim(L1BFile%Name) )
    else if ( the_hdf_version <= 0 ) then                                          
      call MLSL2Message ( MLSMSG_Error, ModuleName, &                      
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )    
    end if

    if ( TOOLKIT ) then
      if ( DEEBUG ) call output( 'Using toolkit', advance='yes' )
    elseif ( .not. (got(p_starttime) .or. got(p_endtime) ) ) then
       call timesFromMafOffsets
    elseif ( LeapSecFileName == '' .and. .not. LEAPSINDATESMODULE ) then
      ! add maf offsets to start, end times
      ! or convert them to tai93
      ! This is optional way to define processingRange if using PCF
      ! It becomes mandatory if not using PCF
       call timesFromMafOffsets
    elseif ( len_trim(LeapSecFileName) > 0 ) then
      call ToolkitTimeConversion
    elseif ( got(p_starttime) .and. got(p_endtime) ) then
      call datesModuleTimeConversion
    end if
    ! call outputNamedValue( 'processingRange%startTime', processingRange%startTime )
    ! call outputNamedValue( 'processingRange%EndTime', processingRange%EndTime )

    if ( .not. TOOLKIT ) then
      call startTable
      call addRow_header ( 'Global Settings', 'c' )
      call addRow_divider ( '-' )
      ! Store appropriate user input as global attributes
      GlobalAttributes%StartUTC = l2pcf%StartUTC
      GlobalAttributes%EndUTC = l2pcf%EndUTC
      if ( len_trim(GlobalAttributes%PGEVersion) < 1 ) &
        & GlobalAttributes%PGEVersion = l2pcf%PGEVersion
      if ( LeapSecFileName /= '' ) call FillTAI93Attribute ( LeapSecFileName )
      ! We don't check on returnStatus--dateless or absolute utc are ok
      call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
        & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
        & GlobalAttributes%GranuleDay) 
      ! If mm < 0, then it's yyyy-Doy format, so let's find the actual
      ! month and store it as a negative integer in GranuleMonth
      if ( GlobalAttributes%GranuleMonth == -1 ) then
        call yyyyDoy_to_mmdd( GlobalAttributes%GranuleYear, &
          & mm, dd, GlobalAttributes%GranuleDay )
        GlobalAttributes%GranuleMonth = -mm
      endif
      ! We'll possibly need the first and last MAF counter numbers, especially
      ! if gaps occur
      ! For now, just look for them in l1boa
      ! Later you may look also in l1brad files
      GlobalAttributes%LastMAFCtr = FindMaxMAF ( L1BFile, &
        & GlobalAttributes%FirstMAFCtr )
      ! call outputNamedValue ( 'Last MAF', GlobalAttributes%LastMAFCtr )
      ! call outputNamedValue ( 'First MAF', GlobalAttributes%FirstMAFCtr )
      call addRow ( 'Last MAF', GlobalAttributes%LastMAFCtr )
      call addRow ( 'First MAF', GlobalAttributes%FirstMAFCtr )
      call outputTable ( sep='|', border='-' )
    end if

    ! Have we overridden the Bright Object names? Can we find them in l1boa?
    if ( got(p_brightObjects) ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
      & 'You overrode names of Bright objects found in l1boa file' )    
    elseif( L1BFile%hdfVersion /= HDFVERSION_5 ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
      & 'Names of Bright objects missing from (hdf4) l1boa file' )    
    elseif( .not. IsHDF5AttributeInFile(L1BFile%name, 'BO_name') ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &                      
      & 'Names of Bright objects missing from l1boa file' )    
    else
      if ( DEEBUG ) call output( 'About to read Bright Objects', advance='yes' )
      wasAlreadyOpen = L1BFile%stillOpen
      if ( .not. wasAlreadyOpen ) call mls_OpenFile(L1BFile)
      L1BFile%fileID%sd_id = 0 ! So we don't look here for the attribute
      ! We still need to open the root group '/'
      call h5gopen_f ( L1BFile%fileID%f_id, '/', L1BFile%fileID%grp_id, ReturnStatus )
      call GetHDF5Attribute( L1BFile, 'BO_name', BO_names )
      call Array2List ( BO_names, BrightObjects )
      call h5gclose_f ( L1BFile%fileID%grp_id, ReturnStatus )
      ! Now were any bright objects visible today?
      call ReadL1BData ( L1BFile,'/GHz/BO_stat', BO_stat, noMAFs, &
        & l1bflag, NeverFail= .true. )
      ! We will pay no mind to "Fill" values
      ! meaning we must zero them out
      BO_stat%intField = max( 0, BO_stat%intField )
      BO_today = ' '
      do i=0, size(BO_names)
        if ( len_trim(BO_names( max(i,1) )) < 1 ) cycle
        if ( any( isBitSet( BO_stat%intField(1, :, :), i ) ) ) then
          if ( i == 0 ) then
            BO_today = catLists( BO_today, 'moon(sp.port)', inseparator=' ')
          else
            BO_today = catLists( BO_today, trim(BO_names(i)), inseparator=' ')
          endif
        endif
      enddo
      if ( .not. wasAlreadyOpen ) call mls_CloseFile(L1BFile)
      call DeallocateL1BData ( BO_stat )
      if ( OUTPUTTHISMONTHSCAL ) then
        call outputCalendar( l2pcf%startutc, dateNote=BO_today, moonphases=.true. )
      else
        call output ( 'Bright objects today:', advance='yes' )
        call output ( trim_safe(BO_today), advance='yes' )
      endif
    endif

    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ! Perhaps dump global settings
    if ( details > -4 ) &
      & call dump_global_settings( processingRange, filedatabase, &
      & DirectDatabase, ForwardModelConfigDatabase, &
      & LeapSecFileName, details-detailReduction )

    ! Take notice of the types of apriori files being used
    if ( APrioriFiles%dao // AprioriFiles%ncep // AprioriFiles%geos5 &
      &  == ' ' ) then
      GlobalAttributes%MiscNotes = catLists(GlobalAttributes%MiscNotes, &
        & 'No gmao or ncep files--falling back to climatology', '\')
    elseif ( APrioriFiles%dao // AprioriFiles%geos5 == ' ' ) then
      GlobalAttributes%MiscNotes = catLists(GlobalAttributes%MiscNotes, &
        & 'No gmao files--falling back to ncep', '\')
    endif

    if ( len_trim(APrioriFiles%dao) > 0 ) then
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & 'apriori(dao)' )
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & trim(APrioriFiles%dao) )
    endif
    if ( len_trim(APrioriFiles%ncep) > 0 ) then
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & 'apriori(ncep)' )
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & trim(APrioriFiles%ncep) )
    endif
    if ( len_trim(APrioriFiles%geos5) > 0 ) then
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & 'apriori(geos5)' )
      GlobalAttributes%MiscNotes = catlists( GlobalAttributes%MiscNotes, &
        & trim(APrioriFiles%geos5) )
    endif
    call FinishUp

  contains

    ! ---------------------------------------------  FinishUp  -----
    ! Any last minute housekeeping tasks
    subroutine FinishUp
      if ( error /= 0 ) &
        & call MLSL2Message ( MLSMSG_Error,ModuleName, &
        & 'Problem with global settings section' )
      if ( warning /= 0 ) &
        & call MLSL2Message ( MLSMSG_Warning,ModuleName, &
        & 'Possible problem with global settings section' )

      if ( specialDumpFile /= ' ' ) &
        & call revertOutput
      call trace_end ( 'SET_GLOBAL_SETTINGS', cond=toggle(gen) )
      if ( timing ) call sayTime ( 'global_settings', cumulative=.false. )

    end subroutine FinishUp

    ! ---------------------------------------  datesModuleTimeConversion  -----
    ! Convert utc time strings to tai93 using dates module
    subroutine datesModuleTimeConversion
      if ( isUTCInRange(start_time_string) .and. &
        & isUTCInRange(end_time_string) ) then
        call output( 'Using datesModule to convert times', advance='yes' )
        processingrange%starttime = utc2tai93s ( start_time_string, leapsec=.true. )
        processingrange%endtime = utc2tai93s ( end_time_string, leapsec=.true. )
      else
        error = 1
        call outputNamedValue ( 'start time', trim(start_time_string) )
        call outputNamedValue ( 'end time', trim(end_time_string) )
        call MLSL2Message( MLSMSG_Warning, ModuleName, &
          & 'start, end times not in range' )
      endif
    end subroutine datesModuleTimeConversion

    ! ---------------------------------------  timesFromMafOffsets  -----
    ! add maf offsets to start, end times
    ! (unless you didn't input start, end times)
    subroutine timesFromMafOffsets
      if ( DEEBUG ) call output( 'Using mafOffsts to convert times', advance='yes' )
      ! 1st--check that have L1BOA
      if ( L1BFile%FileID%f_id < 1 ) then
        call announce_error(son, &
          & 'L1BOA file required by global data--but not set')
      end if

      quantity = 'MAFStartTimeTAI'
      l1bItemName = AssembleL1BQtyName ( quantity, the_hdf_version, .false. )
      if ( DEEBUG ) call output( 'About to read L1B file: ' // trim(l1bItemName), &
        & advance='yes' )
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, dontPad=.true.)
      if ( l1bFlag==-1 ) then
        call announce_error(son, &
          & 'unrecognized MAFStartTimeTAI in L1BOA file')
        minTime = 0.
        maxTime = 0.
      else
        minTime = l1bField%dpField(1,1,1)
        maxTime = l1bField%dpField(1,1,noMAFs) ! This is start time of last MAF
      end if
      call DeallocateL1BData ( l1bField )

    if ( .not. got(p_starttime) ) then
      processingRange%startTime = minTime
    elseif ( startTimeIsAbsolute ) then
      processingRange%startTime = start_time_from_1stMAF
    elseif ( LeapSecFileName /= '' ) then
      if ( DEEBUG ) call output( 'About to read leapsec file' // &
        & trim(LeapSecFileName), advance='yes' )
      if ( precedesUTC ( '1993-01-01', start_time_string )  ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), start_time_string, &
        & processingrange%starttime)
        if ( DEEBUG ) call output( 'Read leapsec file', advance='yes' )
        if ( returnStatus /= 0 ) then
          call announce_error(0, &
          & 'Error converting start time in mls_utctotai; code number: ')
          call output(returnStatus, advance='yes')
        end if
      elseif ( NEGATIVETAI93OK ) then
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        call output( 'The starting time is before our nominal 1993 start date', advance='yes' )
        call output( '(but we allow it to go negative)', advance='yes' )
        processingrange%starttime = &
          & -secondsbetween2utcs ( start_time_string, '1993-01-01' )
        if ( .false. ) call ResetStartingDate( '1961-01-01' ) ! ( '1993-01-01' )
      else
        ! The starting time is before our nominal 1993 start date
        ! We'll use the dates_module as a fallback
        call output( 'The starting time is before our nominal 1993 start date', advance='yes' )
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        processingrange%starttime = &
          & secondsbetween2utcs ( '1961-01-01', start_time_string )
      endif
    else if ( got(p_starttime) ) then
      processingrange%starttime = minTime + start_time_from_1stMAF
    else if ( .not. TOOLKIT ) then
      processingrange%starttime = minTime
    end if

    if ( .not. got(p_starttime) ) then
      processingRange%endTime = maxTime
    elseif ( stopTimeIsAbsolute ) then
      processingRange%endTime = end_time_from_1stMAF
    elseif ( LeapSecFileName /= '' ) then
      if ( precedesUTC ( '1993-01-01', end_time_string ) ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), end_time_string, &
        & processingrange%endtime)
        if ( returnStatus /= 0 ) then
          call announce_error(0, &
          & 'Error converting end time in mls_utctotai; code number: ')
          call output(returnStatus, advance='yes')
        end if
      elseif ( NEGATIVETAI93OK ) then
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        call output( 'The ending time is before our nominal 1993 start date', advance='yes' )
        call output( '(but we allow it to go negative)', advance='yes' )
        processingrange%endtime = &
          & -secondsbetween2utcs ( end_time_string, '1993-01-01' )
        if ( .false. ) call ResetStartingDate( '1961-01-01' ) ! ( '1993-01-01' )
      else
        ! The starting time is before our nominal 1993 start date
        ! We'll use the dates_module as a fallback
        call output( 'The starting time is before our nominal 1993 start date', advance='yes' )
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        processingrange%endtime = &
          & secondsbetween2utcs ( '1961-01-01', end_time_string )
      endif
    else if ( LEAPSINDATESMODULE ) then
      ! Without a leapsec file, let's use date_module's built-in feature
      processingrange%endtime = utc2tai93s ( end_time_string, leapsec=.true. )
    else if ( got(p_endtime) ) then
      processingrange%endtime = minTime + end_time_from_1stMAF
    else if ( .not. TOOLKIT ) then
      processingrange%endtime = maxTime + 1.0
    end if
    end subroutine timesFromMafOffsets

    ! ------------------------------------------  ToolkitTimeConversion  -----
    ! Convert utc time strings to tai93 using toolkit procedures
    subroutine ToolkitTimeConversion
      call output( 'Using toolkit to convert times', advance='yes' )
      if ( DEEBUG ) call output( 'About to read leapsec file' // &
        & trim(LeapSecFileName), advance='yes' )
      if ( precedesUTC ( '1993-01-01', start_time_string ) ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), start_time_string, &
          & processingrange%starttime)
        returnStatus = mls_utctotai(trim(LeapSecFileName), end_time_string, &
          & processingrange%endtime)
        if ( DEEBUG ) call output( 'Read leapsec file', advance='yes' )
        if ( returnStatus /= 0 ) then
          call announce_error(0, &
          & 'Error converting start, end time in mls_utctotai; code number: ')
          call output(returnStatus, advance='yes')
        end if
      elseif ( NEGATIVETAI93OK ) then
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        call output( 'The starting time is before our nominal 1993 start date', advance='yes' )
        call output( '(but we allow it to go negative)', advance='yes' )
        call outputNamedValue( 'start_time_string', trim(start_time_string) )
        call outputNamedValue( 'end_time_string', trim(end_time_string) )
        processingrange%starttime = &
          & -secondsbetween2utcs ( start_time_string, '1993-01-01' )
        processingrange%endtime = &
          & -secondsbetween2utcs ( end_time_string, '1993-01-01' )
        if ( .false. ) call ResetStartingDate( '1961-01-01' ) ! ( '1993-01-01' )
      else
        ! The starting time is before our nominal 1993 start date
        ! We'll use the dates_module as a fallback
        if ( .false. ) call ResetStartingDate( '1961-01-01' )
        processingrange%starttime = &
          & secondsbetween2utcs ( '1961-01-01', start_time_string )
        processingrange%endtime = &
          & secondsbetween2utcs ( '1961-01-01', end_time_string )
      endif
    end subroutine ToolkitTimeConversion

    ! ---------------------------------------------  Announce_Error  -----
    subroutine Announce_Error ( Lcf_where, Full_message, Use_toolkit, &
      & Error_number, just_a_warning )

      ! Arguments

      integer, intent(in) :: Lcf_where
      character(len=*), intent(in) :: Full_message
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
      if ( my_warning ) warning = max(warning,1)

      if ( .not. just_print_it ) then
        if ( .not. my_warning ) then
          error = max(error,1)
          call startErrorMessage ( lcf_where )

          call output ( "The " );
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
    use Dump_1, only: Dump
    use Open_Init, only: DumpL1BDatabase
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

      character(len=*), intent(in) :: LeapSecFileName
      ! The following determines the level of detail to expose:
      ! -1 Skip even counterMAF
      ! -2 Skip all but name (default)
      ! >0 Dump even multi-dim arrays
      integer, intent(in) :: details

      ! Local
      integer ::                              i
      character (len=*), parameter ::         TIME_FORMAT = '(1pD18.12)'
      integer                              ::  hdfVersion
      type(MLSFile_T), pointer             :: L1BFile

      ! Begin
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      hdfVersion = L1BFile%HDFVersion
      if ( hdfversion <= 0 ) &                                          
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &                    
        & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )  
      ! version = 1

      call output ( '============ Global Settings ============', advance='yes' )
      call output ( ' ', advance='yes' )

      call output ( 'L1B database:', advance='yes' )
      call DumpL1BDatabase ( filedatabase, Details )
      call dump(DirectDatabase, Details)

      call startTable
      call addRow_header ( 'Run Info', 'c' )
      call addRow_divider ( '-' )
      call addRow ( 'Start Time', l2pcf%startutc )
      call addRow ( 'End Time', l2pcf%endutc )

      call addRow ( 'Start Time (tai)', processingrange%starttime, format=time_format )
      call addRow ( 'End Time (tai)', processingrange%endtime, format=time_format )
      call addRow ( 'Processing Range', processingrange%endtime-processingrange%starttime )
      if ( LeapSecFileName /= '' ) then
        call addRow ( 'Leap Seconds File', LeapSecFileName )
      endif
      call addRow ( 'PGE version', l2pcf%PGEVersion )
      call addRow ( 'cycle', l2pcf%cycle )
      call addRow ( 'RunID', l2pcf%RunID )
      call addRow ( 'Log file name', l2pcf%logGranID )
      call outputTable ( sep=' ', border='-' )

      ! Dump special hashes
      call NewLine
      call Dump ( countEmpty=.true., &
        & keys=l2pcf%spec_keys, values=l2pcf%spec_mcfnames, &
        & name='l2gp species, mcf : doi', separator=',', &
        & ExtraValues= l2pcf%spec_doinames)
      call blanks ( 80, FillChar='-', advance='yes' )

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

    ! -----------------------------------------------  NotAllowed  -----
    subroutine NotAllowed ( where, why )
      integer, intent(in) :: Where ! tree index
      integer, intent(in) :: Why   ! param_restricted or spec_restricted

      call startErrorMessage ( where )
      select case ( why )
      case ( param_restricted )
        call output ( 'Parameter ' )
        call display_string( parm_indices(param_id) )
      case ( spec_restricted )
        call output ( ' Specification ' )
        call display_string( spec_indices(spec_id) )
      end select
      call output ( ' is not allowed.', advance='yes' )
      call MLSL2Message ( MLSMSG_Error, moduleName, &
        'Prohibited parameter or specification in GlobalSettings section' )
    end subroutine NotAllowed

  ! ---------------------------------------------  proclaim  -----
  subroutine proclaim ( FullName, l1_type, hdfVersion )
    use MLSFiles, only: MLS_HDF_Version, Split_Path_Name, WildCardHDFVersion
    character(len=*), intent(in)   :: FullName
    character(len=*), intent(in)   :: l1_type
    integer, optional,  intent(in) :: hdfVersion

    ! Local variables
    integer                        :: myhdfVersion
    character(len=len(FullName))   :: name, path
    integer, save                  :: trip = 0
    ! Executable
    trip = trip + 1
    if ( trip < 2 ) call StyledOutput ( 'Level 1 products', options='--Banner' )
    call split_path_name ( FullName, path, name )
    call StartTable
    call addRow ( 'type', trim(l1_type) )
    if ( present(hdfVersion) ) then
      if ( hdfVersion == WildCardHDFVersion ) then
        myhdfVersion = mls_hdf_version(trim(Name))
      else
        myhdfVersion = hdfVersion
      end if
      call addRow ( 'hdf version', myhdfVersion )
    end if
    call addRow ( 'path', trim(path) )
    call addRow ( 'name', trim(name) )
    Call OutputTable ( sep='|', border='-' )
  end subroutine proclaim

    ! --------------------------  CreateDirectTypeFromMLSCFInfo  -----
    function CreateDirectTypeFromMLSCFInfo ( root, DirectFile ) result (Direct)
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
    if ( returnStatus /= 0 ) call MLSL2Message ( &
       & MLSMSG_Error, ModuleName, &
       & 'Failed in GetPCFromRef for ' // trim(filename) )

    end function CreateDirectTypeFromMLSCFInfo

  end subroutine Set_Global_Settings

! =====     Private Procedures     =====================================

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Global_Settings

! $Log$
! Revision 2.182  2021/04/01 23:51:37  pwagner
! L1MAFToL2Profile may return matching Lat and Phi
!
! Revision 2.181  2021/03/19 15:20:49  pwagner
! Added the optional arg MIF to MAF v. profile functions
!
! Revision 2.180  2021/02/05 05:20:06  pwagner
! Avoids unassocialed L2GPData
!
! Revision 2.179  2019/04/18 16:30:10  pwagner
! Overwrite GlobalAttributes%PGEVersion only if currently blank
!
! Revision 2.178  2018/11/01 23:16:00  pwagner
! Improve appearance of settings when dumped
!
! Revision 2.177  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.176  2017/11/15 00:11:44  pwagner
! Use OutputTable to Dump list of level 1 files
!
! Revision 2.175  2017/07/27 16:42:38  pwagner
! If yyyy-Doy, store actual actual month number as neg int. in GlobalAttributes%GranuleMonth
!
! Revision 2.174  2017/07/10 18:58:16  pwagner
! CamelCase; abandon all the extra ResetStartingDates
!
! Revision 2.173  2017/01/19 23:52:03  pwagner
! Improve appeearance when dumping Run Info
!
! Revision 2.172  2017/01/07 01:14:51  pwagner
! Print moon phases on calendar
!
! Revision 2.171  2016/11/09 17:20:01  pwagner
! Fixed error when l2cf processes param
!
! Revision 2.170  2016/11/08 17:32:35  pwagner
! Use SayTime subroutine from time_m module; process /reset field
!
! Revision 2.169  2016/09/23 00:11:52  pwagner
! Improve appearance of Dumps
!
! Revision 2.168  2016/09/14 20:08:45  pwagner
! heck if we set start, end times if no l1boa; also that they are reasonable
!
! Revision 2.167  2016/08/09 21:47:42  pwagner
! Fix error in FindMaxMAF; print module names only if verbose
!
! Revision 2.166  2016/07/28 23:39:32  pwagner
! Fixed error introduced with last commit
!
! Revision 2.165  2016/07/28 19:54:30  pwagner
! Guard against older hdf4 l1boa
!
! Revision 2.164  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.163  2016/07/27 23:02:59  pwagner
! Works better with Aircraft-borne instrument data
!
! Revision 2.162  2016/07/22 20:07:16  pwagner
! Fix typos in output
!
! Revision 2.161  2015/09/02 23:17:43  pwagner
! Bomb promptly if no DEM files instead of waiting until run ends
!
! Revision 2.160  2014/03/31 23:50:47  pwagner
! Let tai93 times become negative
!
! Revision 2.159  2014/03/26 17:49:05  pwagner
! Fixed old bugs when startTime, endTime omitted
!
! Revision 2.158  2014/03/19 19:54:47  pwagner
! Fixed a gold-bick breaking bug introduced with last fix
!
! Revision 2.157  2014/03/19 17:37:07  pwagner
! Repaired bug breaking gold brick
!
! Revision 2.156  2014/03/18 17:46:09  pwagner
! Fixed bug in use of got(:) array
!
! Revision 2.155  2014/03/18 17:15:24  pwagner
! Can get leapseconds from dates module if run sans toolkit
!
! Revision 2.154  2014/03/07 19:23:17  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.153  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.152  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.151  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.150  2013/10/09 23:43:17  vsnyder
! Add Evaluate_Variable, declare param correctly
!
! Revision 2.149  2013/09/25 01:04:33  pwagner
! Added DEM stuff
!
! Revision 2.148  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.147  2013/08/30 02:45:49  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.146  2013/08/21 00:27:13  pwagner
! Code around nominal 1993 starting date limitation; awaiting better solution
!
! Revision 2.145  2012/08/16 17:51:07  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.144  2012/04/20 00:46:13  pwagner
! Treats MAF as 0-based; more robustly handles crashed chunks
!
! Revision 2.143  2012/03/28 20:08:24  pwagner
! l2parsf spec generates warning--slave tasks lost ability to join quantities
!
! Revision 2.142  2012/03/15 22:51:15  vsnyder
! Add IGRF_file parameter, some cannonball polishing
!
! Revision 2.141  2011/11/04 00:08:01  pwagner
! Made USEd entities all caps
!
! Revision 2.140  2011/10/05 00:14:45  pwagner
! Added functions to convert between MAFs, profile numbers
!
! Revision 2.139  2011/06/29 21:50:48  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.138  2010/10/19 00:04:36  pwagner
! Tried to fix callstack overflow when restricted
!
! Revision 2.137  2010/05/23 03:19:48  honghanh
! Modify the restriction for the callable forward model
!
! Revision 2.136  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.135  2009/10/05 23:40:35  pwagner
! Moved use statements to module scope to speedup Lahey; this is the last time we do that
!
! Revision 2.134  2009/10/01 19:56:35  vsnyder
! Restrict some functionality in callable forward model mode
!
! Revision 2.133  2009/08/26 17:16:15  pwagner
! Note types of apriori files in file attribute MiscNotes
!
! Revision 2.132  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.131  2009/03/05 16:21:36  pwagner
! fgrid switch dumps each fgrid
!
! Revision 2.130  2008/12/18 21:13:46  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.129  2008/08/27 20:00:59  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.128  2007/12/14 01:50:47  pwagner
! Zero out Fill values from BO_stat
!
! Revision 2.127  2007/10/04 20:43:12  vsnyder
! Remove unused symbols
!
! Revision 2.126  2007/09/24 20:26:03  pwagner
! Prints this months calendar and shows todays bright objects
!
! Revision 2.125  2007/08/17 00:32:41  pwagner
! Unneeded changes
!
! Revision 2.124  2007/06/21 00:55:22  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.123  2007/03/26 18:06:59  pwagner
! Extra debuggung statements; may be useful if get 'Filesize limit exceeded'
!
! Revision 2.122  2007/01/12 00:36:41  pwagner
! New but unused option to skip check for correct month l2pc files
!
! Revision 2.121  2006/10/20 16:52:26  pwagner
! toolkit-using runs will warn but ignore unwanted global settings
!
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
