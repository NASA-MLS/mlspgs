! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module ReadAPriori

  use Expr_M, only: Expr
  use HDF, only: Dfacc_Rdonly
  use HighOutput, only: Dump, OutputNamedValue
  use Init_Tables_Module, only: F_AuraInstrument, &
    & F_Date, F_Dimlist, F_Downsample, &
    & F_Field, F_File, F_Grid, F_HDFVersion, F_MissingValue, F_NoPCFid, &
    & F_Origin, F_Quantitytype, F_Sdname, F_DeferReading, F_Sum, F_Swath, &
    & Field_First, Field_Last, &
    & L_Climatology, L_Dao, L_Geos5, L_Geos5_7, L_Gloria, &
    & L_Merra, L_Merra_2, L_Ncep, L_None, L_Strat, L_Surfaceheight, &
    & S_ChangeSettings, S_Diff, S_Dump, S_Execute, S_Gridded, &
    & S_IsFileAbsent, S_L2aux, S_L2gp, S_ReadGriddedData
  use Intrinsic, only: L_Ascii, L_Binary, L_HDFeos, L_HDF, L_Swath, &
    & Phyq_Dimensionless
  use L2GPData, only: MaxSwathNamesBufSize
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: Filenamelen, MLSFile_T
  use MLSFiles, only: Filenotfound, &
    & HDFVersion_4, HDFVersion_5, WildCardHDFVersion, &
    & AddFileToDatabase, MLS_CloseFile, Dump, GetPCFromRef, InitializeMLSFile, &
    & MLS_HDF_Version, MLS_Inqswath, MLS_OpenFile, Split_Path_Name
  use MLSL2Options, only: CheckPaths, Default_HDFVersion_Read, L2CFNode, &
    & RuntimeValues, SpecialDumpFile, Toolkit, &
    & DumpMacros, MLSL2Message
  use MLSL2Timings, only: AddPhaseToPhaseNames
  use MLSMessageModule, only: MLSMsg_Error, MLSMsg_Warning
  use MLSPCF2, only: &
    & MLSPCF_L2apriori_Start, MLSPCF_L2apriori_End, &
    & MLSPCF_L2clim_Start, MLSPCF_L2clim_End, &
    & MLSPCF_L2dao_Start, MLSPCF_L2dao_End, &
    & MLSPCF_L2geos5_Start, MLSPCF_L2geos5_End, &
    & MLSPCF_L2ncep_Start, MLSPCF_L2ncep_End, &
    & MLSPCF_Surfaceheight_Start, MLSPCF_Surfaceheight_End
  use MLSStringLists, only: CatLists, GetHashElement, SwitchDetail
  use MLSStrings, only: Lowercase
  use MoreTree, only: Get_Boolean
  use Output_M, only: Output, RevertOutput, SwitchOutput
  use PCFHdr, only: GlobalAttributes
  use SDPToolkit, only: Pgs_S_Success
  use String_Table, only: Get_String
  use Toggles, only: Gen, Switches, Toggle
  use Tree, only: Decorate, Decoration, Nsons, &
    & Sub_Rosa, Subtree, Dump_Tree_Node, Where

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types)
! APrioriFiles_T              data type storing names of apriori files used
! APrioriFiles                actual names of apriori files used

!     (subroutines and functions)
! DumpAPrioriAttributes       dump types and names of apriori files used
! ProcessOneAprioriFile       read one apriori file into a gridded data type
! ProcessOneL2AUXFile         read one l2aux file into a gridded data type
! ProcessOneL2GPFile          read one l2gp file into a gridded data type
! Read_apriori                entry point for apriori section; process l2cf
! ReadAPrioriAttributes       read attributes from a file to which they were written
! WriteAPrioriAttributes      write as attributes info about apriori files used
! === (end of toc) ===



  public ::  aprioriFiles, AprioriFiles_t, &
    & DumpAPrioriAttributes, ProcessOneAPrioriFile, &
    & ProcessOneL2AUXFile, ProcessOneL2GPFile, &
    & Read_apriori, ReadAPrioriAttributes, &
    & WriteAPrioriAttributes
  private ::  announce_error
  logical, parameter :: countEmpty = .true. ! Except where overriden locally
  integer, private :: ERROR
  integer, private, parameter :: MAXNUMFILES = 20

   ! What a priori files did we read? 
   ! (This is very wasteful of memory--can we change it
   ! so that just one field is needed/)
  type APrioriFiles_T
    character (len=MAXNUMFILES*FileNameLen) :: l2gp =  ''
    character (len=MAXNUMFILES*FileNameLen) :: l2aux = ''
    character (len=MAXNUMFILES*FileNameLen) :: ncep =  ''
    character (len=MAXNUMFILES*FileNameLen) :: dao =   ''
    character (len=MAXNUMFILES*FileNameLen) :: geos5 = ''
    character (len=MAXNUMFILES*16) :: geos5description =  ''
  end type APrioriFiles_T

  type (APrioriFiles_T), save :: APrioriFiles

  interface readAPrioriAttributes
    module procedure readAPrioriAttributes_id
    module procedure readAPrioriAttributes_mf
  end interface
  
  interface writeAPrioriAttributes
    module procedure writeAPrioriAttributes_id
    module procedure writeAPrioriAttributes_mf
  end interface
  
  ! -----     Private declarations     ---------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------  dumpAPrioriAttributes  -----
  subroutine dumpAPrioriAttributes
    ! dump info about what apriori files were used
    ! Storing them as hdfeos5 attributes
    use Dump_1, only: Dump
    ! Executable
    call output( 'Actual apriori files and file types used', advance='yes' )
    if ( len_trim(APrioriFiles%l2gp) > 0 )&
      & call dump( APrioriFiles%l2gp, 'l2gp' )
    if ( len_trim(APrioriFiles%l2aux) > 0 )&
      & call dump( APrioriFiles%l2aux, 'l2aux' )
    if ( len_trim(APrioriFiles%ncep) > 0 )&
      & call dump( APrioriFiles%ncep, 'ncep' )
    if ( len_trim(APrioriFiles%dao) > 0 )&
      & call dump( APrioriFiles%dao, 'dao' )
    if ( len_trim(APrioriFiles%geos5) > 0 )&
      & call dump( APrioriFiles%geos5, 'geos5' )
  end subroutine dumpAPrioriAttributes

  ! --------------------------------------------------  read_apriori  -----
  ! Read a priori data from data files, be they l2gp, l2aux, climatology,
  ! NCEP, DAO etc.

  subroutine Read_apriori ( Root, L2GPDatabase, L2auxDatabase, GriddedDatabase, &
    & fileDataBase )
    use DumpCommand_M, only: BooleanFromEmptySwath, BooleanFromFormula, &
      & ExecuteCommand, MLSCase, MLSEndSelect, MLSSelect, MLSSelecting, Skip
    use GriddedData, only: GriddedData_T, Dump
    use Init_Tables_Module, only: S_Boolean, S_Case, S_Endselect, &
      & S_Select, S_Skip
    use L2AUXData, only: L2AUXData_T, Dump
    use L2GPData, only: L2GPData_T, Dump
    use MLSL2Timings, only: Section_Times
    use MoreTree, only: Get_Label_And_Spec, Get_Spec_Id
    use Next_Tree_Node_M, only: Next_Tree_Node, Next_Tree_Node_State
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decorate
    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    ! Local variables
    logical :: carryover
    character(len=32) :: cvalue
    integer :: Details             ! How much info about the files to dump
    integer :: KEY
    integer, save :: LastAprioriPCF      ! l2gp or l2aux  apriori
    integer, save :: LastClimPCF         ! l3ascii or gloria format
    integer, save :: LastDAOPCF
    integer, save :: LastGEOS5PCF
    integer, save :: LastHeightPCF
    integer, save :: LastNCEPPCF
    integer :: Me = -1             ! String index for trace
    integer :: NAME                ! Index into string table
    integer :: SON                 ! Of root, an n_spec_args or a n_named
    type(next_tree_node_state) :: State ! of tree traverser
    logical :: verbose
    logical :: verboser

    call trace_begin ( me, "read_apriori", root, cond=toggle(gen) )

    error = 0

    ! Will we be dumping info? To what level of detail?
    Details = switchDetail(switches, 'apr') - 2
    verbose = ( Details > -3 )
    verboser = ( Details > -2 )
    if ( verbose ) &     
    & call output ( '============ Read APriori ============', advance='yes' )  
    ! Are we picking up again where we left off in the last
    ! readAPriori section?
    ! We'll have to pre-process any Booleans in this section to see
    do 
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit
      call get_label_and_spec ( son, name, key )
      L2CFNODE= key
      select case( get_spec_id(key) )
      case ( s_Boolean )
        call decorate ( key,  BooleanFromFormula ( name, key ) )
      end select
    enddo
    call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'carryover', cvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    carryover = ( index(lowercase(cvalue), 't') > 0 )
    if ( verbose ) call dumpMacros
    if ( verboser ) call outputNamedValue( 'carryover', carryover )
    if ( .not. carryover ) then
      LastAprioriPCF = mlspcf_l2apriori_start - 1
      lastClimPCF = mlspcf_l2clim_start - 1
      lastDAOPCF = mlspcf_l2dao_start - 1
      lastNCEPPCF = mlspcf_l2ncep_start - 1
      lastGEOS5PCF = mlspcf_l2geos5_start - 1
      LastHeightPCF = mlspcf_surfaceHeight_start - 1
    end if
    if ( verboser ) call outputNamedValue ( 'mlspcf_l2geos5_start', mlspcf_l2geos5_start )
    if ( verboser ) call outputNamedValue ( 'lastGEOS5PCF', lastGEOS5PCF )
    if ( verboser ) call outputNamedValue ( 'LastClimPCF', LastClimPCF )
    if ( verboser ) call outputNamedValue ( 'carryOver', carryOver )

    do 
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit
      call get_label_and_spec ( son, name, key )
      L2CFNODE= key
      if ( MLSSelecting .and. &
        & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) cycle
      select case( get_spec_id(key) )
      case ( s_Boolean )
        call decorate ( key,  BooleanFromFormula ( name, key ) )
      case ( s_changeSettings ) ! ===============================  changeSettings ==
        ! Change settings for this phase
        call addPhaseToPhaseNames ( 0, key )
      case ( s_isFileAbsent )
        call decorate ( key, BooleanFromEmptySwath ( key ) )
      case ( s_select ) ! ============ Start of select .. case ==========
        ! We'll start seeking a matching case
        call MLSSelect (key)
      case ( s_case ) ! ============ seeking matching case ==========
        ! We'll continue seeking a match unless the case is TRUE
        call MLSCase (key)
      case ( s_endSelect ) ! ============ End of select .. case ==========
        ! We'done with seeking a match
        call MLSEndSelect (key)
      case ( s_execute ) ! ======================== ExecuteCommand ==========
        call ExecuteCommand ( key )
      case ( s_skip ) ! ============================== Skip ==========
        ! We'll skip the rest of the section if the Boolean cond'n is TRUE
        if ( Skip(key) ) exit
      case default
        call processOneAprioriFile ( son, L2GPDatabase, L2auxDatabase, &
          & GriddedDatabase, fileDataBase, &
          & LastAprioriPCF , &
          & LastClimPCF    , &
          & LastDAOPCF     , &
          & LastGEOS5PCF   , &
          & LastHeightPCF  , &
          & LastNCEPPCF     &
            )
      end select
    end do                              ! Lines in l2cf loop
    if( verboser ) then
      call output( '------------------- apriori datatypes --------------', advance='yes' )
      call output ( 'l2gp', advance='yes' )
      call dump( trim(APrioriFiles%l2gp), 'l2gp files' )
      call output ( 'l2aux', advance='yes' )
      call dump( trim(APrioriFiles%l2aux), 'l2aux files' )
      call output ( 'ncep', advance='yes' )
      call dump( trim(APrioriFiles%ncep), 'ncep files' )
      call output ( 'dao', advance='yes' )
      call dump( trim(APrioriFiles%dao), 'dao files' )
      call output ( 'geos5', advance='yes' )
      call dump( trim(APrioriFiles%geos5), 'geos5 files' )
    end if
    if ( ERROR/=0 ) then
      call MLSL2Message ( MLSMSG_Error,ModuleName, &
        & 'Problem with read_apriori section' )
    end if

    call trace_end( "read_apriori", cond=section_times .or. toggle(gen) )

  end subroutine read_apriori
  
  subroutine processOneAprioriFile ( Root, L2GPDatabase, L2auxDatabase, &
    & GriddedDatabase, fileDataBase, &
    & LastAprioriPCF , &
    & LastClimPCF    , &
    & LastDAOPCF     , &
    & LastGEOS5PCF   , &
    & LastHeightPCF  , &
    & LastNCEPPCF     &
      )
    use ChunkDivide_M, only: ChunkDivideConfig
    use Dumpcommand_M, only: Dumpcommand
    use GriddedData, only: Rgr, GriddedData_T, V_Is_Eta, V_Is_Pressure, &
      & AddGriddedDataToDatabase, CopyGrid, &
      & DestroyGriddedData, DownSampleGriddedData, Dump, &
      & SetupnewgriddedData
    use L2AUXData, only: L2AUXData_T, AddL2AUXtoDatabase, &
      & Dump
    use L2GPData, only: L2GPData_T, &
      & Addl2GPtoDatabase, Dump
    use Moretree, only: Get_Label_And_Spec, Get_Spec_Id
    use Ncep_Dao, only: ReadgriddedData
    use ReadgriddedUtils, only: Read_Climatology, ReadgloriaFile
    use Surfaceheight_M, only: Open_Surface_Height_File, &
      & Read_Surface_Height_File, Close_Surface_Height_File
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer :: LastAprioriPCF      ! l2gp or l2aux  apriori
    integer :: LastClimPCF         ! l3ascii or gloria format
    integer :: LastDAOPCF
    integer :: LastGEOS5PCF
    integer :: LastHeightPCF
    integer :: LastNCEPPCF
    ! Local Variables
    integer :: AURAINST             ! index of 'MLS' in AuraInstrument='MLS'
    integer :: DATE             ! in case using GMAO backgr
    character(len=FileNameLen) :: DATESTRING ! 'X,Y,..'
    logical :: DEBUG
    logical :: deferReading      ! skip reading the data; just note file name?
    character(len=8) :: description
    integer :: Details             ! How much info about the files to dump
    integer :: DIMLIST             ! index of 'X,Y,..' in dimList='X,Y,..'
    character(len=FileNameLen) :: DIMLISTSTRING ! 'X,Y,..'
    logical :: DOWNSAMPLE          ! Downsample gridded data to coarser mesh
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: FIELD               ! Son of KEY, must be n_assign
    integer :: FIELDINDEX          ! Literal
    integer :: FieldName           ! sub-rosa index of name in field='name'
    character(len=FileNameLen) :: FIELDNAMESTRING ! actual literal clim. field
    integer :: FileIndex
    integer :: FileName            ! Sub-rosa index of name in file='name'
    character(len=FileNameLen) :: FileNameString   ! actual literal file name
    integer :: FileType            ! s_gridded, s_l2gp, s_l2aux, ..
    logical :: GotAlready               ! Do we need to reread this file?
    type (GriddedData_T), pointer :: Grid     => null()
    type (GriddedData_T), target  :: tempGrid
    type (MLSFile_T) :: GriddedFile
    integer :: GriddedOrigin            ! From tree
    integer :: GridIndex           ! In the griddeddata database
    logical, dimension(field_first:field_last) :: GOT
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    character :: HMOT              ! 'H', 'M', 'O', or 'T'
    integer :: L2apriori_version
    integer :: J                   ! Loop indices for section, spec
    integer :: KEY                 ! Index of n_spec_args in the AST
    type (L2AUXData_T) :: L2AUX
    type (MLSFile_T) :: L2AUXFile
    type (L2GPData_T) :: L2GP
    type (MLSFile_T) :: L2GPFile
    integer :: L2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
!     integer :: lastAPrioriVersion
    character(len=16) :: litDescription
    integer :: Me = -1             ! String index for trace
    real(rgr) ::    missingValue = 0.
    logical :: noPCFid
    character(len=FileNameLen) :: path   ! path of actual literal file name
    integer :: QUANTITYTYPE             ! Lit index of quantity type
    integer :: ReturnStatus
    integer :: SdName        ! sub-rosa index of name in sdName='name'
    character(len=FileNameLen) :: SDNAMESTRING ! actual literal sdName
    character(len=FileNameLen) :: ShortFileName
    integer :: SON              ! Of root, an n_spec_args or a n_named
    character(len=FileNameLen) :: subString   ! file name w/o path
    logical :: sumDelp          ! sum up the DELP field values to get PL?
    integer :: SwathName        ! sub-rosa index of name in swath='name'
    character(len=FileNameLen) :: SWATHNAMESTRING ! actual literal swath name
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR
    integer :: v_type                   ! E.g., v_is_eta
    logical :: verbose

    ! Executable
    call trace_begin ( me, "processOneAprioriFile", root, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    nullify(grid)
    hdfVersion = DEFAULT_HDFVERSION_READ
    HMOT = ' '
    L2apriori_version = 1
    got = .false.
    son = root ! Because first argument of get_label_and_spec is inout
               ! (and it's needed later anyway)
    call get_label_and_spec ( son, l2Name, key )

    ! Node_id(key) is now n_spec_args.

    FileType = get_spec_id(key)

    if ( any( fileType == (/s_diff, s_dump/) ) ) then
      if ( .not. CHECKPATHS ) &
        & call dumpCommand ( key, griddedDataBase=griddedDataBase, &
        & FileDatabase=FileDatabase )
      go to 9
    end if

    ! Now parse file and field names
    Details = switchDetail(switches, 'apr') - 2
    Debug   = ( Details > -2 )
    verbose = ( Details > -3 )
    downsample = .false.
    deferReading = .false.
    noPCFid = .false.
    sumDelp = .false.
    fileName = 0
    gridIndex = 0
    griddedOrigin = l_none
    swathName = 0
    do j = 2, nsons(key)
      field = subtree(j,key)
      L2CFNODE = field
      fieldIndex = decoration(subtree(1,field))
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      case ( f_AuraInstrument )
        AuraInst = sub_rosa(subtree(2,field))
      case ( f_date )
        date = sub_rosa(subtree(2,field))
      case ( f_deferReading )
        deferReading = get_boolean(field)
      case ( f_dimList )
        dimList = sub_rosa(subtree(2,field))
      case ( f_downsample )
        downsample = get_boolean(field)
      case ( f_field )
        fieldName = sub_rosa(subtree(2,field))
      case ( f_file )
        fileName = sub_rosa(subtree(2,field))
      case ( f_grid )
        gridIndex = decoration ( decoration ( subtree(2, field) ) )
      case ( f_missingValue )
        call expr ( subtree(2,field), expr_units, expr_value )
        missingValue = expr_value(1)
      case ( f_hdfVersion ) ! hdfVersion is never used
        call expr ( subtree(2,field), units, value, type )             
        if ( units(1) /= phyq_dimensionless ) &                        
          & call Announce_error ( field, &                               
            & 'No units allowed for hdfVersion: just integer 4 or 5')  
        hdfVersion = value(1)
        call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'HdfVersion is never used by the ReadAPriori command' )
      case ( f_noPCFid )
        noPCFid = get_boolean(field)
      case ( f_origin )
        griddedOrigin = decoration(subtree(2,field))
      case ( f_quantityType )
        quantityType = decoration(subtree(2,field))
      case ( f_sum )
        sumDelp = get_boolean(field)
      case ( f_swath )
        swathName = sub_rosa(subtree(2,field))
      case ( f_sdname )
        sdname = sub_rosa(subtree(2,field))
      end select
    end do

    if ( got(f_dimList) ) then
      call get_string ( dimList, dimListString, strip=.true. )
    else
      dimListString = ''
    end if

    if ( got(f_date) ) then
      call get_string ( date, dateString, strip=.true. )
    else
      dateString = ''
    end if

    if ( got(f_AuraInstrument) ) then
      call get_string ( AuraInst, HMOT, strip=.true. )
    end if

    if ( got(f_file) ) then
      call get_string ( FileName, fileNameString, strip=.true. )
    else
      fileNameString = ''
    end if
    shortFileName = fileNameString
    
    ! So, grid is associated only if the command held a grid field, e.g.
    ! readGriddedData, grid=gmaoTemperatureInput5, origin=...
    ! and not a mere declaration like
    ! gmaoTemperatureInput5: gridded, origin=none
    if ( got(f_grid) ) grid => griddedDataBase(gridIndex)

    if ( FileType == s_readGriddedData ) then
      if ( grid%fileType > 0 ) then
        FileType = grid%fileType
      else
        FileType = s_gridded
      end if
    end if

    select case ( FileType )
    case ( s_l2gp )
      if ( .not. got(f_file) ) &
        & call announce_error ( son, &
          & 'Filename name must be specified in read a priori' )
      swathNameString=''
      if ( got(f_swath) ) &
        & call get_string ( swathName, swathNameString, strip=.true. )
      call processOneL2GPFile ( FileNameString, swathNameString, &
        & noPCFid, mlspcf_l2apriori_start, mlspcf_l2apriori_end, &
        & LastAprioriPCF, HMOT, Debug, &
        & L2GP, L2GPFile )
      FileIndex = AddFileToDataBase( filedatabase, L2GPFile )

      if( Debug ) then
        if ( specialDumpFile /= ' ' ) &
          & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
        call dump( l2gp, details=details )
        if ( specialDumpFile /= ' ' ) &
          & call revertOutput
      end if

      if( switchDetail(switches, 'pro') > -1 ) then                            
         call announce_success( FilenameString, 'l2gp', &                    
         & swathNameString, MLSFile=L2GPFile )    
      end if
      apriorifiles%l2gp = catlists(apriorifiles%l2gp, trim(FilenameString))

      ! Add this l2gp to the database, decorate this key with index
      call decorate ( key, AddL2GPToDatabase( L2GPDatabase, l2gp ) )
      ! Don't call destroy contents as the AddL2GPToDatabase has done a shallow
      ! copy.
      ! l2gp files in our database means               
      ! change default behavior so that we won't allow 
      ! overlaps outside the processing range)         
      call output( '(Resetting defaults to exclude overlaps outside '&
        & // 'processingRange', advance='yes' )
      ChunkDivideConfig%allowPriorOverlaps = .false.
      ChunkDivideConfig%allowPostOverlaps = .false.

    case ( s_l2aux )

      if ( .not. all(got((/f_sdName, f_file, f_quantityType /)))) &
        & call announce_error ( son, &
          & 'file/sd name must both be specified in read a priori' )
      call get_string ( sdName, sdNameString )
      sdNameString = sdNameString(2:LEN_TRIM(sdNameString)-1)


      l2aux%name = l2Name

      l2Index = AddL2AUXToDatabase( L2AUXDatabase, l2aux )
      call processOneL2AUXFile ( FileNameString, sdNameString, &
        & noPCFid, quantityType, &
        & LastAprioriPCF, Debug, &
        & L2AUXDatabase(l2Index), L2AUXFile, fileDatabase )

      call decorate ( key, l2Index )
      ! if( switchDetail(switches, 'apr') > -1 ) then
      if( Debug ) then
        if ( specialDumpFile /= ' ' ) &
          & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
        call dump( L2AUXDatabase(l2Index), details )
        if ( specialDumpFile /= ' ' ) &
          & call revertOutput
      end if

      if( switchDetail(switches, 'pro') > -1 ) then                            
         call announce_success(FilenameString, 'l2aux', &                    
         & sdNameString, MLSFile=L2AUXFile)    
      end if
      apriorifiles%l2aux = catlists(apriorifiles%l2aux, trim(FilenameString))

    case ( s_gridded )

      ! So we can declare a gridded data name to be used later
      if ( griddedOrigin == l_none ) then
        fieldNameString = 'none'
      else if ( .not. all(got((/f_origin, f_field/))) ) then
        call announce_error ( son, 'Incomplete gridded data information' )
      else
        call get_string ( fieldName, fieldNameString, strip=.true. )
      end if

      select case ( griddedOrigin )
      case ( l_none ) ! ----------- Just a declaration to use later
        ! Should not have come here with readGriddedData command
        if ( associated(grid) ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'origin should be sensible for readGriddeddata, "none" is not' )
        description = 'none'
        ! The gridded data needs to part of the database, even if the file
        ! won't be found and the gridded data empty,
        ! so it can be operated on w/o segment faulting
        gridIndex = AddGriddedDataToDatabase( GriddedDatabase, tempGrid )
        call decorate ( key, gridIndex )
        call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
        GriddedDatabase(gridIndex)%fileType = s_gridded
      case ( l_ncep, l_strat ) ! --------------------------- NCEP Data
        if (griddedOrigin == l_ncep) then
           description = 'ncep'
        else
           description = 'strat'
        end if
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_l2ncep_start, mlspcf_l2ncep_end, description, got(f_file), &
          & LastNCEPPCF, returnStatus, noPCFid, verbose, debug )
        FileIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_hdfeos, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
          & PCBottom=mlspcf_l2ncep_start, PCTop=mlspcf_l2ncep_end)
        GriddedFile%PCFId = LastNCEPPCF
        FileIndex = AddFileToDataBase(filedatabase, GriddedFile)
        if ( .not. associated(grid) ) then
          ! The gridded data needs to be part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, tempGrid )
          call decorate ( key, gridIndex )
          if ( returnStatus == PGS_S_SUCCESS) then
            call readGriddedData ( GriddedFile, son, description, &
              & v_is_pressure, GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), missingValue, &
              & verbose=verbose, deferReading=deferReading )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          end if
          GriddedDatabase(gridIndex)%fileType = s_gridded
        else if ( returnStatus == PGS_S_SUCCESS) then
          call readGriddedData ( GriddedFile, son, description, &
            & v_is_pressure, tempGrid, returnStatus, &
            & dimListString, TRIM(fieldNameString), missingValue, &
            & verbose=verbose, deferReading=deferReading )
          ! If we succeeded in reading values, insert them into the db
          ! but only if we succeeded; otherwise protect the current values
          if ( .not. tempGrid%empty .or. deferReading ) then
            call DestroyGriddedData( grid )
            call copyGrid( grid, tempGrid )
          end if
        end if
        if ( returnStatus == 0 ) then
          if( switchDetail(switches, 'pro') > -1 ) &                            
            & call announce_success(FilenameString, 'ncep', &                    
             & fieldNameString, MLSFile=GriddedFile)
          apriorifiles%ncep = catlists(apriorifiles%ncep, trim(FilenameString))
        else
          call announce_success(FilenameString, 'ncep not found--carry on', &                    
             & fieldNameString, MLSFile=GriddedFile)
        end if
      case ( l_dao ) ! ---------------------------- GMAO Data (GEOS4)
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_l2dao_start, mlspcf_l2dao_end, 'dao', got(f_file), &
          & LastDAOPCF, returnStatus, noPCFid, verbose, debug )
        FileIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_hdfeos, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
          & PCBottom=mlspcf_l2dao_start, PCTop=mlspcf_l2dao_end)
        GriddedFile%PCFId = LastDAOPCF
        FileIndex = AddFileToDataBase(filedatabase, GriddedFile)

        ! We will decide whether it's geos4 or geos5 based on the field name
        select case ( lowercase(fieldNameString) )
        case ( 'tmpu' )
          description = 'dao'
          v_type = v_is_pressure
        case ( 'pl', 't' )
          description = 'geos5'
          v_type = V_is_eta
        case default
          description = 'geos5'
          v_type = V_is_eta
        end select

        if ( .not. associated(grid) ) then
          ! The gridded data needs to part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, tempGrid )
          call decorate ( key, gridIndex )

          if ( returnStatus == PGS_S_SUCCESS) then
            call ReadGriddedData ( GriddedFile, son, description, v_type, &
              & GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), &
              & missingValue, litDescription=litDescription, &
              & verbose=verbose, deferReading=deferReading )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          end if
          GriddedDatabase(gridIndex)%fileType = s_gridded
        else if ( returnStatus == PGS_S_SUCCESS) then
          call ReadGriddedData ( GriddedFile, son, description, v_type, &
            & tempGrid, returnStatus, &
            & dimListString, TRIM(fieldNameString), &
            & missingValue, litDescription=litDescription, &
            & verbose=verbose, deferReading=deferReading )
          ! If we succeeded in reading values, insert them into the db
          ! but only if we succeeded; otherwise protect the current values
          if ( .not. tempGrid%empty .or. deferReading ) then
            tempGrid%description = litDescription
            call DestroyGriddedData( grid )
            call copyGrid( grid, tempGrid )
          end if
        end if
        if ( litDescription == 'none' ) then
          call announce_success( FilenameString, 'dao not found--carry on', &
             & fieldNameString )
        elseif ( Griddeddatabase(gridIndex)%empty ) then
          call output( 'File was probably not ' // &
            & trim(litDescription), advance='yes' )
        else if ( description == 'dao' ) then
          apriorifiles%dao = catlists(apriorifiles%dao, trim(FilenameString))
        else
          apriorifiles%geos5 = &
            & catlists(apriorifiles%geos5, trim(FilenameString))
          apriorifiles%geos5Description = &
            & catlists(apriorifiles%geos5Description, trim(litDescription))
        end if
        if ( returnStatus == 0 ) then
          if( switchDetail(switches, 'pro') > -1 ) &                            
            & call announce_success(FilenameString, 'dao', & 
             & fieldNameString, MLSFile=GriddedFile)    
        elseif ( litDescription /= 'none' ) then
          call announce_success(FilenameString, 'dao not found--carry on', &
             & fieldNameString, MLSFile=GriddedFile)
        end if
      case ( l_geos5, l_geos5_7, l_merra, l_merra_2 ) ! --- GMAO Data (GEOS5*)
        if ( verbose ) call outputNamedvalue( 'LastGEOS5PCF', LastGEOS5PCF )
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_l2geos5_start, mlspcf_l2geos5_end, 'geos5', got(f_file), &
          & LastGEOS5PCF, returnStatus, noPCFid, verbose, debug )
        if ( &
          & any( griddedOrigin == (/ l_geos5_7, l_merra_2 /) ) &
          & ) then ! since geos5_7 and merra_2 are HDF5
            FileIndex = InitializeMLSFile(GriddedFile, content='gridded', &
            name=FilenameString, shortName=shortFileName, &
            type=l_hdf, access=DFACC_RDONLY, hdfVersion=HDFVERSION_5, &
            PCBottom=mlspcf_l2geos5_start, PCTop=mlspcf_l2geos5_end)
        else ! and the rest are HDF-EOS
            FileIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
              & name=FilenameString, shortName=shortFileName, &
              & type=l_hdfeos, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
              & PCBottom=mlspcf_l2geos5_start, PCTop=mlspcf_l2geos5_end)
        end if
        GriddedFile%PCFId = LastGEOS5PCF
        FileIndex = AddFileToDataBase( filedatabase, GriddedFile )

        ! We will decide whether it's geos4 or geos5 based on the field name
        select case ( lowercase(fieldNameString) )
        case ( 'tmpu' )
          description = 'dao'
          v_type = v_is_pressure
        case ( 'pl', 't', 'delp' )
          description = 'geos5'
          v_type = V_is_eta
        case default
          ! We may use this case if we allow for either merra or geos5
          ! in the same l2cf declarartion; e.g.
          ! geos5SumInput: gridded, origin=geos5, field='PL,DELP',  ..
          ! which means:
          ! 1st, try to read the field 'PL'
          ! if that fails, try to read the field 'DELP'
          ! (optionally summing over it to create a 'PL' field)
          description = 'geos5'
          v_type = V_is_eta
        end select
        ! We may override these based on the origin field
        select case ( griddedOrigin )
        case ( l_merra ) 
          description = 'merra'
        case ( l_merra_2 ) 
          description = 'merra_2'
        case ( l_geos5_7 ) 
          description = 'geos5_7'
        end select
        if ( DEBUG ) then
          call outputNamedValue( 'fileName', fileNameString )
          call outputNamedValue( 'fieldName', fieldNameString )
          call outputNamedValue( 'description', description )
          call outputNamedValue( 'associated(grid)', associated(grid) )
        end if

        if ( .not. associated(grid) ) then
          ! The gridded data needs to be part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, tempGrid )
          call decorate ( key, gridIndex )
          if ( returnStatus == PGS_S_SUCCESS) then
            call ReadGriddedData ( GriddedFile, son, description, v_type, &
              & GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), &
              & missingValue, dateString, sumDelp, &
              & litDescription=litDescription, &
              & verbose=verbose, deferReading=deferReading )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          end if
          GriddedDatabase(gridIndex)%fileType = s_gridded
        else if ( returnStatus == PGS_S_SUCCESS) then
          call ReadGriddedData ( GriddedFile, son, description, v_type, &
            & tempGrid, returnStatus, &
            & dimListString, TRIM(fieldNameString), &
            & missingValue, dateString, sumDelp, litDescription=litDescription, &
            & verbose=verbose, deferReading=deferReading )
          ! If we succeeded in reading values, insert them into the db
          ! but only if we succeeded; otherwise protect the current values
          if ( .not. tempGrid%empty .or. deferReading ) then
            tempGrid%description = litDescription
            call DestroyGriddedData( grid )
            call copyGrid( grid, tempGrid )
          end if
        end if
        if ( litDescription == 'none' ) then
          call announce_success( FilenameString, 'geos5 not found--carry on', &                    
             & fieldNameString )
        elseif ( Griddeddatabase(gridIndex)%empty ) then
          call output( 'File was probably not ' // &
            & trim(litDescription), advance='yes' )
        else
          apriorifiles%geos5 = &
            & catlists(apriorifiles%geos5, trim(FilenameString))
          apriorifiles%geos5Description = &
            & catlists(apriorifiles%geos5Description, trim(litDescription))
        end if
        if ( returnStatus == 0 ) then
          if( switchDetail(switches, 'pro') > -1 ) &                            
            & call announce_success(FilenameString, description, &                    
             & fieldNameString, MLSFile=GriddedFile)    
        elseif ( litDescription /= 'none' ) then
          call announce_success( FilenameString, 'geos5 not found--carry on', &                    
             & fieldNameString, MLSFile=GriddedFile )
        end if
        ! If we were asked to downsample gridded data to a coarser mesh, do so
        if ( downsample .and. .not. Griddeddatabase(gridIndex)%empty ) then
          call output( 'Downsampling Gridded data', advance='yes' )
          call DestroyGriddedData( tempgrid )
          grid => Griddeddatabase(gridIndex)
          call DownSampleGriddedData ( grid, tempgrid, &
            & heightsStep=1, latsStep=2, lonsStep=2, &
            & lstsStep=1, szasStep=1, datesStep=1, &
            & firstheights=1, firstlats=1, firstlons=1, &
            & firstlsts=1, firstszas=1, firstdates=1 )
          if( Debug ) then
            call output( 'Original Gridded data', advance='yes' )
            call dump( grid, details )
            call output( 'Downsampled Gridded data', advance='yes' )
            call dump( tempgrid, details )
          end if
          call DestroyGriddedData( grid )
          call copyGrid( grid, tempGrid )
        end if
      case ( l_gloria ) ! ------------------------- Data in Gloria's UARS format
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_l2clim_start, mlspcf_l2clim_end, 'gloria', got(f_file), &
          & LastClimPCF, returnStatus, noPCFid, verbose, debug )
        if ( TOOLKIT .and. returnStatus /= PGS_S_SUCCESS ) then
          call announce_error ( son, &
            & 'PCF number not found to supply' // &
            & ' missing Climatology file name' )
        end if
        FileIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_ascii, access=DFACC_RDONLY, &
          & PCBottom=mlspcf_l2clim_start, PCTop=mlspcf_l2clim_end)
        GriddedFile%PCFId = LastCLIMPCF
        FileIndex = AddFileToDataBase(filedatabase, GriddedFile)
        call decorate ( key, &
          & AddGriddedDataToDatabase ( griddedDatabase, &
          & ReadGloriaFile ( GriddedFile ) ) )
        if( switchDetail(switches, 'pro') > -1 ) then                            
          call announce_success(FilenameString, 'Gloria', &                    
           & '', MLSFile=GriddedFile)    
        end if
      case ( l_climatology ) ! -------------------- Climatology data
        ! Identify file (maybe from PCF if no name given)
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_l2clim_start, mlspcf_l2clim_end, 'climatology', got(f_file), &
          & LastClimPCF, returnStatus, noPCFid, verbose, debug )
        if ( Debug ) then
          call outputnamedValue ( 'got(f_file)', got(f_file) )
          call outputnamedValue ( 'fileNameString', trim(fileNameString) )
          call outputnamedValue ( 'path', trim(path) )
          call outputnamedValue ( 'subString', trim(subString) )
          call outputnamedValue ( 'noPCFid', noPCFid )
          call outputnamedValue ( 'returnStatus', returnStatus )
          call outputnamedValue ( 'LastClimPCF', LastClimPCF )
        endif
        if ( TOOLKIT .and. returnStatus /= PGS_S_SUCCESS ) then
          call announce_error ( son, &
            & 'PCF number not found to supply' // &
            & ' missing Climatology file name' )
        end if

        ! Have we read this already?
        gotAlready = associated(GriddedDatabase)
        if ( gotAlready ) then
          gotAlready = any(GriddedDatabase%sourceFilename==filenameString)
        end if
        if ( .not. gotAlready ) then
          FileIndex = InitializeMLSFile(GriddedFile, content = 'clim', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_ascii, access=DFACC_RDONLY, &
            & PCBottom=mlspcf_l2clim_start, PCTop=mlspcf_l2clim_end)
          GriddedFile%PCFId = LastCLIMPCF
          FileIndex = AddFileToDataBase(filedatabase, GriddedFile)
          ! No, well read it then, add its entire contents to the database
          call read_climatology ( GriddedFile, son, &
            & GriddedDatabase, returnStatus, &
            & mlspcf_l2apriori_start, mlspcf_l2apriori_end, &
            & missingValue )
          if ( returnStatus /= 0 ) then
            call Announce_error ( field, &                               
            & 'read_climatology unsuccessful--check file name and path' )
            ! Now crash gracefully instead of getting a reference to a
            ! disassociated GriddedDatabase pointer.
            call MLSL2Message ( MLSMSG_Error, ModuleName, &
            & 'read_climatology unsuccessful--check file name and path' )
          end if
          if ( Debug ) call outputNamedValue( 'climatology desc.', &
            & GriddedDatabase(size(GriddedDatabase))%description )
        end if
        if ( .not. associated(GriddedDatabase) ) go to 9  ! Last chance

        ! Locate requested grid by name, store index in gridIndex
        ! Check that field name is among those added by the source field
        do gridIndex = 1, size(griddedDatabase)
          if ( trim(fieldNameString) == &
            & trim(GriddedDatabase(gridIndex)%quantityName) ) exit
        end do

        if ( gridIndex <= size(griddedDatabase) ) then
          call decorate ( key, gridIndex )
          if( switchDetail(switches, 'pro') > -1 ) then
            if ( .not. gotAlready ) then
              call announce_success(FilenameString, 'climatology', &                  
               & fieldNameString, MLSFile=GriddedFile)
            else
              call announce_success(FilenameString, 'climatology', &                  
               & fieldNameString)
            end if
          end if
        else
          call announce_error ( son, 'Field ' // trim(fieldNameString) // &
            & ' not found in clim. file ' // trim(fileNameString) )
        end if
      case ( l_surfaceHeight ) ! See Read_surface_height_file
        call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
          & mlspcf_surfaceHeight_start, mlspcf_surfaceHeight_end, &
          & 'surfaceHeight', got(f_file), &
          & lastHeightPCF, returnStatus, noPCFid, verbose, debug )
        if ( TOOLKIT .and. returnStatus /= PGS_S_SUCCESS ) then
          call announce_error ( son, &
            & 'PCF number not found to supply' // &
            & ' missing Surface Height file name' )
        end if
        FileIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_binary, access=DFACC_RDONLY, &
          & PCBottom=mlspcf_surfaceHeight_start, PCTop=mlspcf_surfaceHeight_end)
        GriddedFile%PCFId = LastCLIMPCF
        FileIndex = AddFileToDataBase(filedatabase, GriddedFile)
        call open_surface_height_file ( griddedFile )
        call decorate ( key, &
          & AddGriddedDataToDatabase ( griddedDatabase, &
          & read_surface_height_file ( GriddedFile ) ) )
        call close_surface_height_file ( GriddedFile )
        if( switchDetail(switches, 'pro') > -1 ) then                            
          call announce_success(FilenameString, 'SurfaceHeight', &                    
           & '', MLSFile=GriddedFile)    
        end if
      case default ! Can't get here if tree_checker worked correctly
      end select   ! origins of gridded data

      if( Debug .and. gridIndex <= size(griddedDatabase) ) then
        if ( specialDumpFile /= ' ' ) &
          & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
        if ( griddedOrigin /= l_none ) &
          & call dump( GriddedDatabase(gridIndex), details )
        if ( specialDumpFile /= ' ' ) &
          & call revertOutput
      end if

    case default
    end select     ! types of apriori data

9   continue
    call trace_end ( "processOneAprioriFile", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  end subroutine processOneAprioriFile

  subroutine processOneL2AUXFile ( FileNameString, sdNameString, &
    & noPCFid, quantityType, &
    & LastAprioriPCF, Debug, &
    & L2AUX, L2AUXFile, fileDatabase )
   use L2AUXData, only: L2AUXData_T, ReadL2AUXData
   ! Process an a priori l2gp
    ! Args:
    character(len=FileNameLen)              :: FileNameString   ! actual literal file name
    character(len=FileNameLen)              :: sdNameString ! actual literal swath name
    logical, intent(in)                     :: noPCFid
    integer, intent(in)                     :: QUANTITYTYPE  ! Lit index of quantity type
    integer, intent(inout)                  :: LastAprioriPCF
    logical, intent(in)                     :: Debug
    type (L2AUXData_T), intent(inout)       :: L2AUX
    type (MLSFile_T)                        :: L2AUXFile
    type (MLSFile_T), dimension(:), pointer :: FILEDATABASE
    ! Internal variables
    integer :: FileIndex
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: L2apriori_version               ! 4 or 5 (corresp. to hdf4 or hdf5)

    character(len=FileNameLen) :: path   ! path of actual literal file name
    type (MLSFile_T), pointer :: pL2AUXFile
    integer :: ReturnStatus
    character(len=FileNameLen) :: ShortFileName
    character(len=FileNameLen) :: subString   ! file name w/o path

    if ( TOOLKIT .and. .not. noPCFid ) then
      call split_path_name( FileNameString, path, SubString )
      LastAprioriPCF = GetPCFromRef( SubString, mlspcf_l2apriori_start, &
      & mlspcf_l2apriori_end, &                                     
      & TOOLKIT, returnStatus, l2apriori_Version, DEBUG, &  
      & exactName=FileNameString )                             
    end if
    hdfVersion = mls_hdf_version(FilenameString)
    FileIndex = InitializeMLSFile( L2AUXFile, content = 'l2aux', &
      & name=FilenameString, shortName=shortFileName, &
      & type=l_hdf, access=DFACC_RDONLY, hdfVersion=hdfVersion, &
      & PCBottom=mlspcf_l2apriori_start, PCTop=mlspcf_l2apriori_end )
    ! call mls_openFile(L2AUXFile, returnStatus)
    L2AUXFile%PCFId = LastAprioriPCF
    FileIndex = AddFileToDataBase( filedatabase, L2AUXFile )
    pL2AUXFile => filedatabase(FileIndex)
    call ReadL2AUXData ( pL2AUXFile, sdNameString, &
      & L2AUX, &
      & quantityType, &
      & checkDimNames=.false. )

  end subroutine processOneL2AUXFile

  subroutine processOneL2GPFile ( FileNameString, swathNameString, &
    & noPCFid, PCBottom, PCTop, &
    & LastAprioriPCF, HMOT, Debug, &
    & L2GP, L2GPFile )
    ! Process an a priori l2gp
    use L2GPData, only: L2GPData_T, &
      & Readl2GPData
    ! Args:
    character(len=FileNameLen)     :: FileNameString   ! actual literal file name
    character(len=FileNameLen)     :: SWATHNAMESTRING ! actual literal swath name
    logical, intent(in)            :: noPCFid
    integer, intent(in)            :: PCBottom
    integer, intent(in)            :: PCTop
    integer, intent(inout)         :: LastAprioriPCF
!    integer, intent(inout)         :: LastAprioriVersion
    character, intent(in)          :: HMOT           ! 'H', 'M', 'O', or 'T'
    logical, intent(in)            :: Debug
    type (L2GPData_T), intent(out) :: L2GP
    type (MLSFile_T) :: L2GPFile
    ! Internal variables
    character(len=MAXSWATHNAMESBUFSIZE) :: ALLSWATHNAMES ! Buffer to get info back.
    integer :: COMMAPOS                 ! For parsing string
    integer :: FileIndex
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: L2apriori_version               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: LISTSIZE                 ! Size of string from SWInqSwath
    integer :: NOSWATHS                 ! In an input file
    character(len=FileNameLen) :: path   ! path of actual literal file name
    integer :: ReturnStatus
    character(len=FileNameLen) :: ShortFileName
    character(len=FileNameLen) :: subString   ! file name w/o path

    ! If we were given only a strand of the filename, expand it
    L2apriori_version = 1
    if ( TOOLKIT .and. .not. noPCFid ) then
      call split_path_name( FileNameString, path, SubString )
      LastAprioriPCF = GetPCFromRef( SubString, mlspcf_l2apriori_start, &
      & mlspcf_l2apriori_end, &                                     
      & TOOLKIT, returnStatus, l2apriori_Version, DEBUG, &  
      & exactName=FileNameString )                             
    end if
    hdfVersion = mls_hdf_version(FilenameString)

    FileIndex = InitializeMLSFile( L2GPFile, content = 'l2gp', &
      & name=FilenameString, shortName=shortFileName, &
      & type=l_swath, access=DFACC_RDONLY, hdfVersion=hdfVersion, &
      & PCBottom=PCBottom, PCTop=PCTop )
    L2GPFile%PCFId = LastAprioriPCF
    ! If we didn't get a name get the first swath name in the file
    if ( len_trim(swathNameString) == 0 ) then
      allSwathNames = ''
      noSwaths = mls_InqSwath ( fileNameString, allSwathNames, listSize, &
       & hdfVersion=hdfVersion )
      if ( listSize == FILENOTFOUND ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'File not found; make sure the name and path are correct' &
          & // trim(fileNameString), MLSFile=L2GPFile )
      else if ( listSize < 1 ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'Failed to determine swath names, perhaps none in file ' &
          & // trim(fileNameString), MLSFile=L2GPFile )
      else if ( listSize < len(allSwathNames) ) then
        commaPos = index ( allSwathNames, ',' )
        if ( commaPos == 0 ) then
          commaPos = len_trim(allSwathNames)
        else if ( commaPos == 1 ) then
          call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'Failed to determine swath name, allswathnames begin with , ' &
          & // trim(fileNameString), MLSFile=L2GPFile )
        else
          commaPos = commaPos - 1
        end if
        swathNameString = allSwathNames ( 1:commaPos )
      else
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'Failed to determine swath names, string too long.' &
          & // trim(fileNameString), MLSFile=L2GPFile )
      end if
    end if
    if ( swathNameString == ' ' ) then
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Failed to determine swath name, obscure error on ' &
        & // trim(fileNameString), MLSFile=L2GPFile )
    end if

    ! Read the swath
    if ( HMOT /= ' ' ) then
      call ReadL2GPData ( L2GPFile, swathNameString, l2gp, HMOT=HMOT )
    else
      call ReadL2GPData ( L2GPFile, swathNameString, l2gp )
    end if
  end subroutine processOneL2GPFile

  ! =========== Private ============

    subroutine Get_PCF_Id ( FileNameString, Path, SubString, L2Apriori_Version, &
      & FirstPCF, LastPCF, Description, GotFile, PCF_Id, ReturnStatus, noPCFid, &
      & verbose, debug )
      use MLSFiles, only: Getpcfromref, Split_Path_Name
      use MLSL2Options, only: Toolkit
      use SDPToolkit, only: Pgs_Pc_Getreference, Pgs_S_Success
      use Toggles, only: Gen, Levels, Toggle
      use Trace_M, only: Trace_Begin, Trace_End

      ! Args
      character(len=*), intent(inout) :: FileNameString
      character(len=*), intent(out) :: Path, SubString
      integer, intent(inout) :: L2Apriori_Version
      integer, intent(in) :: FirstPCF, LastPCF
      character(len=*), intent(in) :: Description
      logical, intent(in) :: GotFile
      integer, intent(inout) :: PCF_Id
      integer, intent(out) :: ReturnStatus
      logical, intent(in) :: noPCFid
      logical, intent(in) :: verbose
      logical, intent(in) :: debug
      ! Local variables
      integer :: Me = -1       ! String index for trace
      integer :: PCFBottom
      ! Executable
      call trace_begin ( me, "Get_PCF_Id", &
        & cond=toggle(gen) .and. levels(gen) > 1 )
      if ( gotFile ) then
        PCFBottom = FirstPCF
      else
        PCFBottom = PCF_Id + 1 ! This must be the last PCF id we found
      end if

      ! Now parse file and field names
      if ( debug ) then
        call outputNamedValue ( 'fileNameString', trim(fileNameString) )
        call outputNamedValue ( 'got file?', gotFile )
        call output('PCFBottom, lastPCF: ')
        call output( (/ PCFBottom, lastPCF /), advance='yes' )
      endif
      call split_path_name ( fileNameString, path, subString )
      if ( TOOLKIT .and. gotFile .and. .not. noPCFid ) then
        if ( debug ) call output( 'Calling getPCFromRef with ' // trim(subString), advance='yes' )
        pcf_id = getPCFromRef ( subString, PCFBottom, lastPCF, TOOLKIT, &
          &                     returnStatus, l2Apriori_Version, DEBUG, &
          &                     exactName = fileNameString )
      else if ( TOOLKIT .and. .not. noPCFid ) then
        do pcf_id = PCFBottom, lastPCF
          returnStatus = Pgs_pc_getReference(pcf_id, L2apriori_version, &
                & fileNameString)
          if ( returnStatus == PGS_S_SUCCESS) exit
        end do
        if ( returnStatus /= PGS_S_SUCCESS ) &
          & call announce_success ( description, 'no entry in PCF',  ' ' )
      else
        returnStatus = 0
        PCF_Id = 0
      end if
      call trace_end ( "Get_PCF_Id", &
        & cond=toggle(gen) .and. levels(gen) > 1 )
    end subroutine Get_PCF_Id

  ! ------------------------------------------  readAPrioriAttributes_MF  -----
  subroutine readAPrioriAttributes_MF ( MLSFile )
    type(MLSFile_T) :: MLSFile
    logical             :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'apr') > 0 )
    if ( verbose ) call dump( MLSFile, details=1 )
    if ( MLSFile%hdfVersion /= HDFVERSION_5 ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Wrong hdfVersion--can read apriori attributes for hdf5 only', &
        & MLSFile=MLSFile )
      return ! Can only do this for hdf5 files
    else if ( MLSFile%StillOpen ) then
      if ( verbose ) call output( 'About to read file attributes for a priori', advance='yes' )
      call readAPrioriAttributes_ID(MLSFile%fileID%f_id, HDFVERSION_5)
    else
      call MLS_OpenFile( MLSFile )
      if ( verbose ) then
        call dump( MLSFile, details=1 )
        call output( 'About to read file attributes for a priori', advance='yes' )
      end if
      call readAPrioriAttributes_ID(MLSFile%fileID%f_id, HDFVERSION_5)
      call MLS_CloseFile( MLSFile )
    end if
  end subroutine readAPrioriAttributes_MF

  ! ------------------------------------------  readAPrioriAttributes_ID  -----
  subroutine readAPrioriAttributes_ID ( fileID, hdfVersion )
    ! read info about what apriori files were used
    ! Storing them as hdfeos5 attributes
    use MLSHDFEOS, only: HE5_EHRDGlatt, MLS_IsGlatt
    use PCFHdr, only: GlobalAttributes
    ! Args
    integer, intent(in) :: fileID
    integer, intent(in) :: hdfVersion  ! Must be 5 to work properly
    ! Internal variables
    integer             :: status
    logical             :: verbose
    character(len=*), parameter  :: whereami = 'readAPrioriAttributes_ID'
    ! Executable
    verbose = ( switchDetail(switches, 'apr') > 0 )
    if ( verbose ) call output( 'Reading apriori attributes', advance='yes' )
    if ( hdfVersion /= HDFVERSION_5 ) then
      call MLSL2Message ( MLSMSG_Warning, whereami, &
        & 'Wrong hdfVersion--can read apriori attributes for hdf5 only' )
      return ! Can only do this for hdf5 files
    end if
    if ( verbose ) then
      call outputNamedValue( 'FileID', FileID )
      call outputNamedValue( 'l2gp there?', mls_isglatt ( fileID, 'A Priori l2gp' ) )
      call outputNamedValue( 'geos5desc there?', mls_isglatt ( fileID, 'geos5 type' ) )
    end if
    status = HE5_EHRDGLATT(fileID, &
     & 'A Priori l2gp', APrioriFiles%l2gp)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%l2gp' // trim(APrioriFiles%l2gp) )
    status = HE5_EHRDGLATT(fileID, &
     & 'A Priori l2aux', APrioriFiles%l2aux)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%l2aux' // trim(APrioriFiles%l2aux) )
    status = HE5_EHRDGLATT(fileID, &
     & 'A Priori ncep', APrioriFiles%ncep)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%ncep' // trim(APrioriFiles%ncep) )
    status = HE5_EHRDGLATT(fileID, &
     & 'A Priori gmao', APrioriFiles%dao)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%dao' // trim(APrioriFiles%dao) )
    status = HE5_EHRDGLATT(fileID, &
     & 'A Priori geos5', APrioriFiles%geos5)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%geos5' // trim(APrioriFiles%geos5) )
    status = HE5_EHRDGLATT(fileID, &
     &  'geos5 type', APrioriFiles%geos5description)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading APrioriFiles%geos5description' // &
      &  trim(APrioriFiles%geos5description) )
    status = HE5_EHRDGLATT(fileID, &
     &  'MiscNotes', GlobalAttributes%MiscNotes)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem reading GlobalAttributes%MiscNotes' // &
      &  trim(GlobalAttributes%MiscNotes) )
  end subroutine readAPrioriAttributes_ID

  ! ------------------------------------------  writeAPrioriAttributes_MF  -----
  subroutine writeAPrioriAttributes_MF ( MLSFile, DontReplace )
    type(MLSFile_T) :: MLSFile
    logical, optional, intent(in) :: DontReplace ! Don't overwrite attributes
    logical             :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'apr') > 0 )
    if ( verbose ) then
      call dump( MLSFile, details=1 )
      call output( 'About to write file attributes for a priori', advance='yes' )
    end if
    if ( MLSFile%hdfVersion /= HDFVERSION_5 ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Wrong hdfVersion--can write apriori attributes for hdf5 only', &
        & MLSFile=MLSFile )
      return ! Can only do this for hdf5 files
    else if ( MLSFile%StillOpen ) then
      call writeAPrioriAttributes_ID( MLSFile%fileID%f_id, HDFVERSION_5, &
       & DontReplace )
    else
      call MLS_OpenFile( MLSFile )
      call writeAPrioriAttributes_ID( MLSFile%fileID%f_id, HDFVERSION_5, &
        & DontReplace )
      call MLS_CloseFile( MLSFile )
    end if
  end subroutine writeAPrioriAttributes_MF

  ! ------------------------------------------  writeAPrioriAttributes_ID  -----
  subroutine writeAPrioriAttributes_ID ( fileID, hdfVersion, DontReplace )
    ! Write info about what apriori files were used
    ! Storing them as hdfeos5 attributes
    use HDFEOS5, only: MLS_Chartype
    use MLSHDFEOS, only: MLS_EHWRGlatt, MLS_IsGlatt
    ! Args
    integer, intent(in) :: fileID
    integer, intent(in) :: hdfVersion  ! Must be 5 to work properly
    logical, optional, intent(in) :: DontReplace ! Don't overwrite attributes
    character(len=*), parameter  :: whereami = 'readAPrioriAttributes_ID'
    ! Internal variables
    integer             :: status
    logical             :: verbose
    ! Executable
    verbose = ( switchDetail(switches, 'apr') > 0 )
    if ( verbose ) call output( 'Writing apriori attributes', advance='yes' )
    if ( hdfVersion /= HDFVERSION_5 ) then
      call MLSL2Message ( MLSMSG_Warning, whereami, &
        & 'Wrong hdfVersion--can write apriori attributes for hdf5 only' )
      return ! Can only do this for hdf5 files
    end if
    if ( verbose ) call outputNamedValue( 'attributes already written', &
      & mls_isglatt ( fileID, 'A Priori l2gp' ) )
    ! Have the attributes been written before? Must we avoid replacing them?
    if ( present( dontReplace ) ) then
      if ( verbose ) call outputNamedValue( 'dontReplace?', dontReplace )
      if ( dontReplace ) then
        if ( mls_isglatt ( fileID, 'A Priori l2gp' ) ) return
      end if
    end if
    status = mls_EHwrglatt(fileID, &
     & 'A Priori l2gp', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%l2gp))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%l2gp' // trim(APrioriFiles%l2gp) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori l2aux', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%l2aux))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%l2aux' // trim(APrioriFiles%l2aux) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori ncep', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%ncep))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%ncep' // trim(APrioriFiles%ncep) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori gmao', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%dao))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%dao' // trim(APrioriFiles%dao) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori geos5', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%geos5))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%geos5' // trim(APrioriFiles%geos5) )
    status = mls_EHwrglatt(fileID, &
     & 'geos5 type', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%geos5description))
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing APrioriFiles%geos5description' // &
      & trim(APrioriFiles%geos5description) )
    status = mls_EHwrglatt(fileID, &
      & 'MiscNotes', MLS_CHARTYPE, 1, &
      &  GlobalAttributes%MiscNotes)
    if ( status /= 0 ) &
      &  call MLSL2Message ( MLSMSG_Warning, whereami, &
      & 'Problem writing MiscNotes' // &
      & trim(GlobalAttributes%MiscNotes) )
  end subroutine writeAPrioriAttributes_ID

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( FullName, l2_type, quantityName, &
    & hdfVersion, MLSFile )
    use HighOutput, only: AddRow, &
      & OutputTable, StartTable, StyledOutput
    ! Args
    character(len=*), intent(in)   :: FullName
    character(len=*), intent(in)   :: l2_type
    integer, optional,  intent(in) :: hdfVersion
    character(len=*), intent(in) :: quantityName
    type(MLSFile_T), optional :: MLSFile

    ! Local variables
    integer                        :: myhdfVersion
    character(len=len(FullName))   :: name, path
    integer, save                  :: trip = 0
    ! Executable
    trip = trip + 1
    if ( trip < 2 ) &
      & call StyledOutput ( 'Level 2 apriori products', options='--Banner' )
    call split_path_name ( FullName, path, name )
    call StartTable
    call addRow ( 'type', trim(l2_type) )
    if ( present(hdfVersion) ) then
      if ( hdfVersion == WILDCARDHDFVERSION ) then
        myhdfVersion = mls_hdf_version(trim(Name))
      else
        myhdfVersion = hdfVersion
      end if
      call addRow ( 'hdf ver', myhdfVersion )
    end if
    call addRow ( 'path', trim(path) )
    call addRow ( 'name', trim(name) )
    call addRow ( 'quantity', trim(quantityName) )
    Call OutputTable ( sep='|', border='-' )
    if ( present(MLSFile) .and. switchdetail(switches, 'apr' ) > -1 ) &
      & call dump(MLSFile)
  end subroutine announce_success

  ! ------------------------------------------------  announce_error  -----
  subroutine Announce_error ( lcf_where, full_message, use_toolkit, &
    & error_number, APrioriFile )
  
   ! Arguments
  
    integer, intent(in)    :: Lcf_where
    character(LEN=*), intent(in)    :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional    :: Error_number
    type(MLSFile_T), intent(in), optional    :: APrioriFile
    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.
    logical :: verbose
    ! Executable
    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
    verbose = ( SwitchDetail(switches, 'apr') > -1 )
 
    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( where(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if
      if ( verbose ) then
        call output ( ': ' )
        call output ( "The " );
        if ( lcf_where > 0 ) then
          call dump_tree_node ( lcf_where, 0 )
        else
          call output ( '(no lcf tree available)' )
        end if
      end if

      call output ( " Caused the following error: ", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
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
    if ( present(APrioriFile) ) call dump(APrioriFile)
!===========================
  end subroutine Announce_error
!===========================

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ReadAPriori

!=============================================================================

!
! $Log$
! Revision 2.128  2019/01/31 19:49:58  pwagner
! Improved verbose Dumps
!
! Revision 2.127  2018/11/01 23:18:13  pwagner
! Removed unused stuff
!
! Revision 2.126  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.125  2018/03/22 18:14:28  pwagner
! Added command IsFileAbsent; may occur in ReadApriori, MergeGrids, and Output sections
!
! Revision 2.124  2018/03/14 22:42:40  pwagner
! May changeSettings in readApriori and MergeGrids sections
!
! Revision 2.123  2017/11/15 00:11:24  pwagner
! Use OutputTable to Dump list of level 1 files
!
! Revision 2.122  2017/08/08 20:46:28  vsnyder
! Stop with error instead of seg fault if climatology not successfully read
!
! Revision 2.121  2017/07/10 23:04:52  pwagner
! Print less if not verbose
!
! Revision 2.120  2017/03/17 00:13:43  pwagner
! runtime Boolean CarryOver given better chance to work properly
!
! Revision 2.119  2017/03/07 21:21:01  pwagner
! Support new meteorology origin: merra_2
!
! Revision 2.118  2017/01/13 01:31:40  pwagner
! Avoid excess output when feeling around for right meteorology file type
!
! Revision 2.117  2016/09/30 20:33:46  pwagner
! verboser was used w/o being defined; fixed
!
! Revision 2.116  2016/09/21 00:41:04  pwagner
! Default to printing less
!
! Revision 2.115  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.114  2016/04/01 00:27:15  pwagner
! May now Execute a single command or a script of lines from l2cf
!
! Revision 2.113  2015/11/17 21:28:07  pwagner
! Made public processOneL2AUXFile and processOneL2GPFile
!
! Revision 2.112  2015/08/03 21:43:50  pwagner
! Made quantityType optional in call to ReadL2AUXData
!
! Revision 2.111  2015/05/07 20:16:50  pwagner
! Fixed error affecting noPCFid
!
! Revision 2.110  2015/05/05 18:15:47  pwagner
! /noPCFid field allows us to read from files not named in PCF
!
! Revision 2.109  2014/10/07 00:08:45  pwagner
! Higher details level of debug now required for most printing
!
! Revision 2.108  2014/10/02 17:24:00  pwagner
! Added MiscNotes to file attributes we read and write
!
! Revision 2.107  2014/09/05 01:27:33  vsnyder
! Add some tracing
!
! Revision 2.106  2014/06/04 18:34:37  pwagner
! runtime carryover flag starts PCFids at value from last readApriori section; may deferReading
!
! Revision 2.105  2014/04/02 23:05:11  pwagner
! Removed redundant open_ and close_MLSFile
!
! Revision 2.104  2014/01/11 01:52:15  vsnyder
! Get OutputNamedValue from HighOutput, not Output_m.  Getting it from
! Output_m crept back in during CVS conflict resolution.
!
! Revision 2.103  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.102  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.101  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.100  2013/10/09 23:43:41  vsnyder
! Add Evaluate_Variable
!
! Revision 2.99  2013/09/24 23:47:23  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.98  2013/08/30 02:45:51  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.97  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.96  2013/02/12 18:13:01  pwagner
! Raise -Sapr switch setting needed to become verbose
!
! Revision 2.95  2012/08/16 17:50:42  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.94  2012/05/08 17:50:31  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.93  2012/03/12 17:31:59  pwagner
! New api for writeAPrioriAttributes; can avoid replacing if already written
!
! Revision 2.92  2012/02/24 21:16:06  pwagner
! Read/Write geos5description attributes, too
!
! Revision 2.91  2011/08/03 21:59:39  pwagner
! Reduced default verbosity by making debug switchable
!
! Revision 2.90  2011/07/13 14:24:46  honghanh
! Initialize GEOS5.7 file type with l_hdf instead of l_hdfeos
!
! Revision 2.89  2011/07/12 22:35:03  honghanh
! Change l_grid to l_hdfeos
!
! Revision 2.88  2011/06/16 23:16:27  pwagner
! Added /downsample to reduce resolution of gridded data when read
!
! Revision 2.87  2011/05/26 20:45:15  pwagner
! Should not bomb with 'pro' set in switches
!
! Revision 2.86  2011/05/05 17:02:34  pwagner
! Added readGriddedData command to readApriori section
!
! Revision 2.85  2011/04/27 17:39:56  pwagner
! Consistent with new ncep_dao api
!
! Revision 2.84  2011/04/20 16:54:28  pwagner
! Added new flexibility to l2cf control flow by run-time booleans affecting gridded data
!
! Revision 2.83  2010/11/09 02:37:28  vsnyder
! Spiff up a dump
!
! Revision 2.82  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.81  2009/10/26 17:11:07  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.80  2009/09/10 23:03:41  pwagner
! Prevents 'Dump, /stop' line in l2cf from causing checkPaths failure
!
! Revision 2.79  2009/08/26 16:48:17  pwagner
! Added readAPrioriAttributes; writes APrioriFiles%geos5 as file attribute
!
! Revision 2.78  2009/08/24 20:25:28  pwagner
! Added merra file type, 'DELP' field name, /sum filed
!
! Revision 2.77  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.76  2008/12/02 23:12:47  pwagner
! mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!
! Revision 2.75  2008/09/17 23:20:13  pwagner
! Allow date string in gridded data to offset gmao background files
!
! Revision 2.74  2007/10/24 00:16:59  pwagner
! Removed unused declarations
!
! Revision 2.73  2007/08/17 00:35:30  pwagner
! Unneeded changes
!
! Revision 2.72  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.71  2007/03/02 18:14:02  pwagner
! Fixed bugs introduced 2 versions ago
!
! Revision 2.70  2007/01/30 21:59:10  pwagner
! fixed bug where Get_PCF_Id retruned undefined values
!
! Revision 2.69  2007/01/11 20:48:30  vsnyder
! Add SurfaceHeight to gridded data, vector quantities, allow dump in ReadApriori
!
! Revision 2.68  2006/06/13 20:56:43  pwagner
! Fixed bug in pcfid names for geos5 files
!
! Revision 2.67  2006/06/13 18:19:08  pwagner
! Added pcfids for geos5; moved ncep pcfids backwards to make room
!
! Revision 2.66  2006/06/06 21:56:52  pwagner
! May specify geos5 apriori files instead of dao (geos4)
!
! Revision 2.65  2006/04/10 23:45:18  pwagner
! Reset defaults in read_apriori, not ChunkDivide
!
! Revision 2.64  2006/02/10 21:15:55  pwagner
! dumps may go to special dumpfile
!
! Revision 2.63  2006/01/26 00:35:35  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.62  2005/09/28 17:02:04  pwagner
! Should not segment fault when reading apriori l2aux
!
! Revision 2.61  2005/08/05 20:39:07  pwagner
! L2AUXFile arg to ReadL2AUXFile now a pointer
!
! Revision 2.60  2005/07/12 17:34:57  pwagner
! Added MLSFile interface for writing apriori attributes
!
! Revision 2.59  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.58  2005/06/14 20:40:27  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.57  2004/09/23 23:03:14  pwagner
! Added writeAPrioriAttributes
!
! Revision 2.56  2004/06/28 20:25:54  pwagner
! Handle dao, ncep missing from PCF with grace
!
! Revision 2.55  2004/06/23 17:13:34  pwagner
! Should quit gracefully if climatolgy file not found
!
! Revision 2.54  2004/01/23 01:10:58  pwagner
! Gets max swathlist length from L2GPData
!
! Revision 2.53  2003/10/06 13:16:09  cvuu
! add new description=strat to handle reading the ncep data file
!
! Revision 2.52  2003/06/09 22:49:34  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.51  2003/05/29 17:54:28  pwagner
! Able to read dao, ncep files w/o knowing name fragment
!
! Revision 2.50  2003/05/09 23:26:45  pwagner
! Should not bomb when l2aux hdfversion absent or wildcard
!
! Revision 2.49  2003/05/06 00:16:32  pwagner
! Fixed passing wild card hdf version to mls_sfend for l2aux
!
! Revision 2.48  2003/05/05 23:00:34  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.47  2003/04/17 23:08:29  pwagner
! Added optional AuraInstrument field to l2gp apriori reads
!
! Revision 2.46  2003/04/04 18:35:11  pwagner
! Warns if dao, ncep FILENOTFOUND
!
! Revision 2.45  2003/04/02 23:54:53  pwagner
! Checks for FILENOTFOUND
!
! Revision 2.44  2003/04/01 21:44:38  dwu
! Fixed bad error with multiple swaths in file (pwagner)
!
! Revision 2.43  2003/03/01 00:24:27  pwagner
! Added missingValue as filed to reading Gridded data
!
! Revision 2.42  2003/02/27 18:41:40  pwagner
! Handles WILDCARDHDFVERSION properly
!
! Revision 2.41  2003/02/20 21:26:21  pwagner
! Lets you read field dimList=x,y,.. w/ griddeddata
!
! Revision 2.40  2003/02/19 19:16:21  pwagner
! More sensible output when proclaiming successful file inputs
!
! Revision 2.39  2003/01/18 02:37:34  livesey
! Added the quantityType stuff to L2Aux issues
!
! Revision 2.38  2002/12/10 00:40:27  pwagner
! Overrides defaults; forcing check of l2auxdimNames
!
! Revision 2.37  2002/12/06 01:06:56  pwagner
! Passes hdfVersion to readl2auxdata
!
! Revision 2.36  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.35  2002/04/06 00:11:18  pwagner
! Combined pcf numbers for dao, ncep, clim, l2gp into apriori
!
! Revision 2.34  2002/01/29 23:49:38  pwagner
! Separate DEFAULT_HDFVERSION_(READ)(WRITE)
!
! Revision 2.33  2002/01/26 00:10:45  pwagner
! Correctly sets hdfVersion; changed proclaim to announce_success
!
! Revision 2.32  2002/01/23 23:09:46  pwagner
! Handles optional hdfVersion field; proclaims files input
!
! Revision 2.31  2002/01/23 22:35:47  livesey
! Added ReadGloriaFile stuff
!
! Revision 2.30  2002/01/18 19:01:34  pwagner
! Changed debugOption to .false.
!
! Revision 2.29  2002/01/18 18:50:08  pwagner
! Better check when swopen fails
!
! Revision 2.28  2002/01/18 00:55:30  pwagner
! Uses MLSFiles for swapi wrappers
!
! Revision 2.27  2002/01/09 00:00:04  pwagner
! Replaced write or print statements with calls to output
!
! Revision 2.26  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.25  2001/10/30 00:35:27  pwagner
! Tidied up small things
!
! Revision 2.24  2001/10/26 23:20:10  pwagner
! Optionally dumps apriori quantities as it reads them
!
! Revision 2.23  2001/10/08 21:35:43  pwagner
! Initialize allswathnames
!
! Revision 2.22  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.21  2001/09/10 23:37:44  livesey
! New GriddedData stuff
!
! Revision 2.20  2001/05/12 00:20:00  livesey
! Allowed user to not supply swath name when reading l2gp
!
! Revision 2.19  2001/05/07 18:03:56  pwagner
! Checks for PCF before looping over pcf_indx
!
! Revision 2.18  2001/05/03 20:34:08  vsnyder
! Cosmetic changes
!
! Revision 2.17  2001/04/16 23:50:01  pwagner
! Tiny change to announce_error
!
! Revision 2.16  2001/04/12 22:19:33  vsnyder
! Improved an error message
!
! Revision 2.15  2001/04/10 20:04:26  livesey
! Bug fixes etc.
!
! Revision 2.14  2001/03/30 00:27:38  pwagner
! Fleshed out njl outline for reading clim. files
!
! Revision 2.13  2001/03/29 19:13:41  livesey
! Added some comment guidelines for read climatology
!
! Revision 2.12  2001/03/21 00:46:08  pwagner
! Passes son to READ_CLIMATOLOGY
!
! Revision 2.11  2001/03/15 21:38:00  pwagner
! Gets v_is_pressure from GriddedData
!
! Revision 2.10  2001/03/15 21:25:16  pwagner
! Split between GriddedData and ncep_dao modules
!
! Revision 2.9  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.8  2001/03/15 00:34:34  pwagner
! Gives more info to ReadGriddedData
!
! Revision 2.7  2001/03/14 19:06:08  livesey
! Some changes
!
! Revision 2.6  2001/03/14 18:54:38  pwagner
! Uses FieldNameString and son in call to ReadGriddedData
!
! Revision 2.5  2001/03/08 01:08:08  pwagner
! Interfaces with ReadGriddedData
!
! Revision 2.4  2001/03/07 22:41:44  livesey
! Reworked the l2cf aspects
!
! Revision 2.3  2001/03/07 01:04:33  pwagner
! No longer uses obtainclim, obtaindao, obtainncep
!
! Revision 2.2  2001/03/06 00:23:58  pwagner
! A little bit more
!
! Revision 2.1  2001/03/03 00:14:40  pwagner
! First commit
!
!
