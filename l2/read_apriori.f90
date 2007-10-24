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
  use Hdf, only: DFACC_RDWR, DFACC_RDONLY
  use INIT_TABLES_MODULE, only: F_AURAINSTRUMENT, F_DIMLIST, F_FIELD, F_FILE, &
    & F_HDFVERSION, F_missingValue, F_ORIGIN, F_QUANTITYTYPE, F_SDNAME, F_SWATH, &
    & FIELD_FIRST, FIELD_LAST, L_CLIMATOLOGY, L_DAO, L_NCEP, &
    & L_GEOS5, L_GLORIA, L_STRAT, L_SURFACEHEIGHT, &
    & S_Dump, S_GRIDDED, S_L2AUX, S_L2GP
  use Intrinsic, only: l_ascii, L_Binary, l_grid, l_hdf, l_swath, PHYQ_Dimensionless
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: FileNameLen, MLSFile_T
  use MLSFiles, only: FILENOTFOUND, &
    & HDFVERSION_4, HDFVERSION_5, WILDCARDHDFVERSION, &
    & AddFileToDataBase, Dump, GetPCFromRef, InitializeMLSFile, &
    & MLS_HDF_VERSION, &
    & MLS_INQSWATH, mls_io_gen_closeF, mls_io_gen_openF, &
    & SPLIT_PATH_NAME
  use MLSL2Options, only: DEFAULT_HDFVERSION_READ, SPECIALDUMPFILE, TOOLKIT
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
  use MLSMessageModule, only: MLSMessage, MLSMessageCalls, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSPCF2, only: &
    & mlspcf_l2apriori_start, mlspcf_l2apriori_end, &
    & mlspcf_l2clim_start, mlspcf_l2clim_end, &
    & mlspcf_l2dao_start, mlspcf_l2dao_end, &
    & mlspcf_l2geos5_start, mlspcf_l2geos5_end, &
    & mlspcf_l2ncep_start, mlspcf_l2ncep_end, &
    & mlspcf_surfaceHeight_start, mlspcf_surfaceHeight_end
  use MLSStringLists, only: catLists, SWITCHDETAIL
  use MLSStrings, only: lowercase
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: BLANKS, OUTPUT, &
    & revertoutput, switchOutput
  use SDPToolkit, only: PGS_S_SUCCESS
  use String_Table, only: GET_STRING
  use Time_M, only: Time_Now
  use TOGGLES, only: GEN, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
  use TREE_TYPES, only: N_NAMED

  implicit none
  private
  public ::  APrioriFiles, APrioriFiles_T, read_apriori, writeAPrioriAttributes
  private ::  announce_error
  integer, private :: ERROR
  integer, private, parameter :: MAXNUMFILES = 10

   ! What a priori files did we read? 
  type APrioriFiles_T
    character (len=MAXNUMFILES*FileNameLen) :: l2gp = ''
    character (len=MAXNUMFILES*FileNameLen) :: l2aux = ''
    character (len=MAXNUMFILES*FileNameLen) :: ncep = ''
    character (len=MAXNUMFILES*FileNameLen) :: dao = ''
    character (len=MAXNUMFILES*FileNameLen) :: geos5 = ''
  end type APrioriFiles_T

  type (APrioriFiles_T), save :: APrioriFiles
  interface writeAPrioriAttributes
    module procedure writeAPrioriAttributes_id
    module procedure writeAPrioriAttributes_MF
    module procedure writeAPrioriAttributes_name
  end interface
  
  ! -----     Private declarations     ---------------------------------

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  read_apriori  -----
  ! Read a priori data from data files, be they l2gp, l2aux, climatology,
  ! NCEP, DAO etc.

  subroutine Read_apriori ( Root, L2GPDatabase, L2auxDatabase, GriddedDatabase, &
    & fileDataBase)

    use ChunkDivide_m, only: ChunkDivideConfig
    use DumpCommand_m, only: DumpCommand
    use GriddedData, only: rgr, GriddedData_T, V_is_eta, v_is_pressure, &
      & AddGriddedDataToDatabase, Dump, SetupNewGriddedData
    use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase, &
      &                  ReadL2AUXData, Dump
    use L2GPData, only: L2GPData_T, MAXSWATHNAMESBUFSIZE, &
      & AddL2GPToDatabase, ReadL2GPData, Dump
    use ncep_dao, only: READ_CLIMATOLOGY, ReadGriddedData, ReadGloriaFile
    use SurfaceHeight_m, only: Open_Surface_Height_File, &
      & Read_Surface_Height_File, Close_Surface_Height_File

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    ! Local Variables
    integer :: AURAINST             ! index of 'MLS' in AuraInstrument='MLS'
    integer :: COMMAPOS                 ! For parsing string
    logical, parameter :: DEBUG = .FALSE.
    integer :: Details             ! How much info about the files to dump
    integer :: DIMLIST             ! index of 'X,Y,..' in dimList='X,Y,..'
    character(len=FileNameLen) :: DIMLISTSTRING ! 'X,Y,..'
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: FIELD               ! Son of KEY, must be n_assign
    integer :: FIELDINDEX          ! Literal
    integer :: FieldName        ! sub-rosa index of name in field='name'
    character(len=FileNameLen) :: FIELDNAMESTRING ! actual literal clim. field
    integer :: FileName            ! Sub-rosa index of name in file='name'
    character(len=FileNameLen) :: FileNameString   ! actual literal file name
    integer :: FileType            ! either s_l2gp or s_l2aux
    logical, dimension(field_first:field_last) :: GOT
    type (griddedData_T) :: GriddedData1
    type (MLSFile_T) :: GriddedFile
    integer :: GriddedOrigin            ! From tree
    integer :: GridIndex           ! In the griddeddata database
    logical :: GotAlready               ! Do we need to reread this file?
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    character :: HMOT              ! 'H', 'M', 'O', or 'T'
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Index of n_spec_args in the AST
    integer :: L2apriori_version
    type (L2AUXData_T) :: L2AUX
    type (MLSFile_T) :: L2AUXFile
    type (MLSFile_T), pointer :: pL2AUXFile
    type (L2GPData_T) :: L2GP
    type (MLSFile_T) :: L2GPFile
    integer :: L2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
    integer :: LastAprioriPCF      ! l2gp or l2aux  apriori
    integer :: LastClimPCF         ! l3ascii or gloria format
    integer :: LastDAOPCF
    integer :: LastGEOS5PCF
    integer :: LastHeightPCF
    integer :: LastNCEPPCF
    integer :: LISTSIZE                 ! Size of string from SWInqSwath
    real(rgr) ::    missingValue = 0.
    integer :: NOSWATHS                 ! In an input file
    character(len=FileNameLen) :: path   ! path of actual literal file name
    integer :: QUANTITYTYPE             ! Lit index of quantity type
    integer :: ReturnStatus
    character(len=FileNameLen) :: ShortFileName
    integer :: SON              ! Of root, an n_spec_args or a n_named
    integer :: SdName        ! sub-rosa index of name in sdName='name'
    character(len=FileNameLen) :: SDNAMESTRING ! actual literal sdName
    character(len=FileNameLen) :: subString   ! file name w/o path
    integer :: SwathName        ! sub-rosa index of name in swath='name'
    character(len=FileNameLen) :: SWATHNAMESTRING ! actual literal swath name
    real :: T1, T2                      ! for timing
    logical :: TIMING
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR
    integer :: v_type                   ! E.g., v_is_eta

    character(len=MAXSWATHNAMESBUFSIZE) :: ALLSWATHNAMES ! Buffer to get info back.
    character(len=8) :: description

    if ( toggle(gen) ) then
      call trace_begin ( "read_apriori", root )
    else
      call MLSMessageCalls( 'push', constantName=ModuleName )
    endif

    timing = section_times
    if ( timing ) call time_now ( t1 )
    error = 0

    ! Will we be dumping info? To what level of detail?
    Details = switchDetail(switches, 'apr') - 2
    if( Details > -3 ) &     
    & call output ( '============ Read APriori ============', advance='yes' )    
    allswathnames = ' '
    LastAprioriPCF = mlspcf_l2apriori_start - 1
    lastClimPCF = mlspcf_l2clim_start - 1
    lastDAOPCF = mlspcf_l2dao_start - 1
    lastNCEPPCF = mlspcf_l2ncep_start - 1
    lastGEOS5PCF = mlspcf_l2geos5_start - 1
    LastHeightPCF = mlspcf_surfaceHeight_start - 1
    ! call outputNamedValue ( 'mlspcf_l2geos5_start', mlspcf_l2geos5_start )
    ! call outputNamedValue ( 'lastGEOS5PCF', lastGEOS5PCF )

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      hdfVersion = DEFAULT_HDFVERSION_READ
      HMOT = ' '
      L2apriori_version = 1
      got = .false.
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        l2Name = sub_rosa(subtree(1,son))
      else
        key = son
        l2Name = 0
      end if

      ! Node_id(key) is now n_spec_args.

      FileType = get_spec_id(key)

      if ( fileType == s_dump ) then
        call dumpCommand ( key, griddedDataBase=griddedDataBase )
        cycle
      end if

      ! Now parse file and field names
      fileName = 0
      swathName = 0
      do j = 2, nsons(key)
        field = subtree(j,key)
        fieldIndex = decoration(subtree(1,field))
        got(fieldIndex) = .true.
        select case ( fieldIndex )
        case ( f_AuraInstrument )
          AuraInst = sub_rosa(subtree(2,field))
        case ( f_dimList )
          dimList = sub_rosa(subtree(2,field))
        case ( f_field )
          fieldName = sub_rosa(subtree(2,field))
        case ( f_file )
          fileName = sub_rosa(subtree(2,field))
        case ( f_missingValue )
          call expr ( subtree(2,field), expr_units, expr_value )
          missingValue = expr_value(1)
        case ( f_hdfVersion )           
          call expr ( subtree(2,field), units, value, type )             
          if ( units(1) /= phyq_dimensionless ) &                        
            & call Announce_error ( field, &                               
              & 'No units allowed for hdfVersion: just integer 4 or 5')  
          hdfVersion = value(1)                                          
        case ( f_origin )
          griddedOrigin = decoration(subtree(2,subtree(j,key)))
        case ( f_quantityType )
          quantityType = decoration(subtree(2,subtree(j,key)))
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
        
      if ( got(f_AuraInstrument) ) then
        call get_string ( AuraInst, HMOT, strip=.true. )
      endif

      if ( got(f_file) ) then
        call get_string ( FileName, fileNameString, strip=.true. )
      else
        fileNameString = ''
      end if
      shortFileName = fileNameString
        
      select case ( FileType )
      case ( s_l2gp )
        if ( .not. got(f_file) ) &
          & call announce_error ( son, &
            & 'Filename name must be specified in read a priori' )
        swathNameString=''
        if ( got(f_swath) ) &
          & call get_string ( swathName, swathNameString, strip=.true. )

        ! If we were given only a strand of the filename, expand it
        if ( TOOLKIT ) then
          call split_path_name(FileNameString, path, SubString)
          LastAprioriPCF = GetPCFromRef(SubString, mlspcf_l2apriori_start, &
          & mlspcf_l2apriori_end, &                                     
          & TOOLKIT, returnStatus, l2apriori_Version, DEBUG, &  
          & exactName=FileNameString)                             
        endif
        hdfVersion = mls_hdf_version(FilenameString)

        gridIndex = InitializeMLSFile(L2GPFile, content = 'l2gp', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_swath, access=DFACC_RDONLY, hdfVersion=hdfVersion, &
          & PCBottom=mlspcf_l2apriori_start, PCTop=mlspcf_l2apriori_end)
        L2GPFile%PCFId = LastAprioriPCF
        gridIndex = AddFileToDataBase(filedatabase, L2GPFile)
        ! If we didn't get a name get the first swath name in the file
        if ( len_trim(swathNameString) == 0 ) then
          allSwathNames = ''
          noSwaths = mls_InqSwath ( fileNameString, allSwathNames, listSize, &
           & hdfVersion=hdfVersion)
          if ( listSize == FILENOTFOUND ) then
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'File not found; make sure the name and path are correct' &
              & // trim(fileNameString), MLSFile=L2GPFile )
          elseif ( listSize < 1 ) then
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Failed to determine swath names, perhaps none in file ' &
              & // trim(fileNameString), MLSFile=L2GPFile )
          elseif ( listSize < len(allSwathNames) ) then
            commaPos = index ( allSwathNames, ',' )
            if ( commaPos == 0 ) then
              commaPos = len_trim(allSwathNames)
            elseif ( commaPos == 1 ) then
              call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Failed to determine swath name, allswathnames begin with , ' &
              & // trim(fileNameString), MLSFile=L2GPFile )
            else
              commaPos = commaPos - 1
            endif
            swathNameString = allSwathNames ( 1:commaPos )
          else
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Failed to determine swath names, string too long.' &
              & // trim(fileNameString), MLSFile=L2GPFile )
          end if
        endif
        if ( swathNameString == ' ' ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to determine swath name, obscure error on ' &
            & // trim(fileNameString), MLSFile=L2GPFile )
        endif

        ! Read the swath
        if ( HMOT /= ' ' ) then
          call ReadL2GPData ( L2GPFile, swathNameString, l2gp, HMOT=HMOT )
        else
          call ReadL2GPData ( L2GPFile, swathNameString, l2gp )
        endif

        if( Details > -3 ) then
          if ( specialDumpFile /= ' ' ) &
            & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
          call dump( l2gp, details=details )
          if ( specialDumpFile /= ' ' ) &
            & call revertOutput
        endif

        if( switchDetail(switches, 'pro') > -1 ) then                            
           call announce_success(FilenameString, 'l2gp', &                    
           & swathNameString, MLSFile=L2GPFile)    
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

        if ( TOOLKIT ) then
          call split_path_name(FileNameString, path, SubString)
          LastAprioriPCF = GetPCFromRef(SubString, mlspcf_l2apriori_start, &
          & mlspcf_l2apriori_end, &                                     
          & TOOLKIT, returnStatus, l2apriori_Version, DEBUG, &  
          & exactName=FileNameString)                             
        endif
        if ( .not. all(got((/f_sdName, f_file, f_quantityType /)))) &
          & call announce_error ( son, &
            & 'file/sd name must both be specified in read a priori' )
        
        call get_string ( sdName, sdNameString )
        sdNameString = sdNameString(2:LEN_TRIM(sdNameString)-1)

        hdfVersion = mls_hdf_version(FilenameString)
        gridIndex = InitializeMLSFile(L2AUXFile, content = 'l2aux', &
          & name=FilenameString, shortName=shortFileName, &
          & type=l_hdf, access=DFACC_RDONLY, hdfVersion=hdfVersion, &
          & PCBottom=mlspcf_l2apriori_start, PCTop=mlspcf_l2apriori_end)
        ! call mls_openFile(L2AUXFile, returnStatus)
        L2AUXFile%PCFId = LastAprioriPCF

        gridIndex = AddFileToDataBase(filedatabase, L2AUXFile)
        l2aux%name = l2Name

        l2Index = AddL2AUXToDatabase( L2AUXDatabase, l2aux )
        call decorate ( key, l2Index )
        pL2AUXFile => filedatabase(gridIndex)
        call ReadL2AUXData ( pL2AUXFile, sdNameString, quantityType, &
          & L2AUXDatabase(l2Index), &
          & checkDimNames=.false. )

        ! if( index(switches, 'apr') /= 0 ) then
        if( Details > -3 ) then
          if ( specialDumpFile /= ' ' ) &
            & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
          call dump( L2AUXDatabase(l2Index), details )
          if ( specialDumpFile /= ' ' ) &
            & call revertOutput
        endif

        ! call mls_closeFile(L2AUXFile, returnStatus)
!         if ( returnStatus /= 0 ) then
!           call announce_error ( son, &
!             & 'Failed to close l2aux file ' // trim(FileNameString) )
!         elseif(index(switches, 'pro') /= 0) then                            
        if( switchDetail(switches, 'pro') > -1 ) then                            
           call announce_success(FilenameString, 'l2aux', &                    
           & sdNameString, MLSFile=L2AUXFile)    
        end if
        apriorifiles%l2aux = catlists(apriorifiles%l2aux, trim(FilenameString))

      case ( s_gridded )

        if ( .not. all(got((/f_origin, f_field/))) ) &
          & call announce_error ( son, 'Incomplete gridded data information' )

        call get_string ( fieldName, fieldNameString, strip=.true. )
        
        select case ( griddedOrigin )
        case ( l_ncep, l_strat ) ! --------------------------- NCEP Data
          if (griddedOrigin == l_ncep) then
             description = 'ncep'
          else
             description = 'strat'
          end if
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_l2ncep_start, mlspcf_l2ncep_end, description, got(f_file), &
            & LastNCEPPCF, returnStatus )
          gridIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_grid, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
            & PCBottom=mlspcf_l2ncep_start, PCTop=mlspcf_l2ncep_end)
          GriddedFile%PCFId = LastNCEPPCF
          gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
          ! The gridded data needs to part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, GriddedData1 )
          call decorate ( key, gridIndex )
          if ( returnStatus == PGS_S_SUCCESS) then
            call readGriddedData ( GriddedFile, son, description, &
              & v_is_pressure, GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), missingValue )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          endif
          if ( returnStatus == 0 ) then
            if( switchDetail(switches, 'pro') > -1 ) &                            
              & call announce_success(FilenameString, 'ncep', &                    
               & fieldNameString, MLSFile=GriddedFile)
            apriorifiles%ncep = catlists(apriorifiles%ncep, trim(FilenameString))
          else
            call announce_success(FilenameString, 'ncep not found--carry on', &                    
               & fieldNameString, MLSFile=GriddedFile)
          endif
        case ( l_dao ) ! ---------------------------- GMAO Data (GEOS4)
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_l2dao_start, mlspcf_l2dao_end, 'dao', got(f_file), &
            & LastDAOPCF, returnStatus )
          gridIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_grid, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
            & PCBottom=mlspcf_l2dao_start, PCTop=mlspcf_l2dao_end)
          GriddedFile%PCFId = LastDAOPCF
          gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
          ! The gridded data needs to part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, GriddedData1 )
          call decorate ( key, gridIndex )

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
          if ( returnStatus == PGS_S_SUCCESS) then
            call ReadGriddedData ( GriddedFile, son, description, v_type, &
              & GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), &
              & missingValue )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          endif
          if ( returnStatus == 0 ) then
            if( switchDetail(switches, 'pro') > -1 ) &                            
              & call announce_success(FilenameString, 'dao', &                    
               & fieldNameString, MLSFile=GriddedFile)    
            if ( description == 'dao' ) then
              apriorifiles%dao = catlists(apriorifiles%dao, trim(FilenameString))
            else
              apriorifiles%geos5 = catlists(apriorifiles%geos5, trim(FilenameString))
            endif
          else
            call announce_success(FilenameString, 'dao not found--carry on', &                    
               & fieldNameString, MLSFile=GriddedFile)
          endif
        case ( l_geos5 ) ! ---------------------------- GMAO Data (GEOS5)
          ! call outputNamedValue ( 'fileNameString', trim(fileNameString) )
          ! call outputNamedValue ( 'LastGEOS5PCF', LastGEOS5PCF )
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_l2geos5_start, mlspcf_l2geos5_end, 'geos5', got(f_file), &
            & LastGEOS5PCF, returnStatus )
          gridIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_grid, access=DFACC_RDONLY, hdfVersion=HDFVERSION_4, &
            & PCBottom=mlspcf_l2geos5_start, PCTop=mlspcf_l2geos5_end)
          GriddedFile%PCFId = LastGEOS5PCF
          gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
          ! The gridded data needs to part of the database, even if the file
          ! won't be found and the gridded data empty,
          ! so it can be merged w/o segment faulting
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, GriddedData1 )
          call decorate ( key, gridIndex )

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
          if ( returnStatus == PGS_S_SUCCESS) then
            call ReadGriddedData ( GriddedFile, son, description, v_type, &
              & GriddedDatabase(gridIndex), returnStatus, &
              & dimListString, TRIM(fieldNameString), &
              & missingValue )
          else
            call SetupNewGriddedData ( GriddedDatabase(gridIndex), empty=.true. )
          endif
          if ( returnStatus == 0 ) then
            if( switchDetail(switches, 'pro') > -1 ) &                            
              & call announce_success(FilenameString, 'geos5', &                    
               & fieldNameString, MLSFile=GriddedFile)    
            if ( description == 'dao' ) then
              apriorifiles%dao = catlists(apriorifiles%dao, trim(FilenameString))
            else
              apriorifiles%geos5 = catlists(apriorifiles%geos5, trim(FilenameString))
            endif
          else
            call announce_success(FilenameString, 'geos5 not found--carry on', &                    
               & fieldNameString, MLSFile=GriddedFile)
          endif
          ! error = 1
        case ( l_gloria ) ! ------------------------- Data in Gloria's UARS format
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_l2clim_start, mlspcf_l2clim_end, 'gloria', got(f_file), &
            & LastClimPCF, returnStatus )
          if ( TOOLKIT .and. returnStatus /= PGS_S_SUCCESS ) then
            call announce_error ( son, &
              & 'PCF number not found to supply' // &
              & ' missing Climatology file name' )
          end if
          gridIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_ascii, access=DFACC_RDONLY, &
            & PCBottom=mlspcf_l2clim_start, PCTop=mlspcf_l2clim_end)
          GriddedFile%PCFId = LastCLIMPCF
          gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
          call decorate ( key, &
            & AddGriddedDataToDatabase ( griddedDatabase, &
            & ReadGloriaFile ( GriddedFile ) ) )
          if( switchDetail(switches, 'pro') > -1 ) then                            
            call announce_success(FilenameString, 'Gloria', &                    
             & '', MLSFile=GriddedFile)    
          endif
        case ( l_climatology ) ! -------------------- Climatology data
          ! Identify file (maybe from PCF if no name given)
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_l2clim_start, mlspcf_l2clim_end, 'climatology', got(f_file), &
            & LastClimPCF, returnStatus )
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
            gridIndex = InitializeMLSFile(GriddedFile, content = 'clim', &
              & name=FilenameString, shortName=shortFileName, &
              & type=l_ascii, access=DFACC_RDONLY, &
              & PCBottom=mlspcf_l2clim_start, PCTop=mlspcf_l2clim_end)
            GriddedFile%PCFId = LastCLIMPCF
            gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
            ! No, well read it then, add its entire contents to the database
            call read_climatology ( GriddedFile, son, &
              & GriddedDatabase, returnStatus, &
              & mlspcf_l2apriori_start, mlspcf_l2apriori_end, &
              & missingValue )
            if ( returnStatus /= 0 ) &
              & call Announce_error ( field, &                               
              & 'read_climatology unsuccessful--check file name and path')  
          end if
          if ( .not. associated(GriddedDatabase) ) cycle  ! Last chance
       
          ! Locate requested grid by name, store index in gridIndex
          ! Check that field name is among those added by the source field
          do gridIndex = 1, size(griddedDatabase)
            if ( trim(fieldNameString) == &
              & trim(GriddedDatabase(gridIndex)%quantityName) ) exit
          end do

          if ( gridIndex <= size(griddedDatabase) ) then
            call decorate ( key, gridIndex )
            if( switchDetail(switches, 'pro') > -1 ) then                            
              call announce_success(FilenameString, 'climatology', &                  
               & fieldNameString, MLSFile=GriddedFile)
            end if
          else
            call announce_error ( son, 'Field ' // trim(fieldNameString) // &
              & ' not found in clim. file ' // trim(fileNameString) )
          end if
        case ( l_surfaceHeight ) ! See Read_surface_height_file
          call get_pcf_id ( fileNameString, path, subString, l2apriori_version, &
            & mlspcf_surfaceHeight_start, mlspcf_surfaceHeight_end, &
            & 'surfaceHeight', got(f_file), &
            & lastHeightPCF, returnStatus )
          if ( TOOLKIT .and. returnStatus /= PGS_S_SUCCESS ) then
            call announce_error ( son, &
              & 'PCF number not found to supply' // &
              & ' missing Surface Height file name' )
          end if
          gridIndex = InitializeMLSFile(GriddedFile, content = 'gridded', &
            & name=FilenameString, shortName=shortFileName, &
            & type=l_binary, access=DFACC_RDONLY, &
            & PCBottom=mlspcf_surfaceHeight_start, PCTop=mlspcf_surfaceHeight_end)
          GriddedFile%PCFId = LastCLIMPCF
          gridIndex = AddFileToDataBase(filedatabase, GriddedFile)
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

        if( Details > -3 ) then
          if ( specialDumpFile /= ' ' ) &
            & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
          call dump( GriddedDatabase(gridIndex), details )
          if ( specialDumpFile /= ' ' ) &
            & call revertOutput
        endif

      case default
      end select     ! types of apriori data
    end do                              ! Lines in l2cf loop
    
    if ( ERROR/=0 ) then
      call MLSMessage(MLSMSG_Error,ModuleName, &
        & 'Problem with read_apriori section')
    end if
    
    if ( toggle(gen) ) then
      call trace_end("read_apriori")
    else
      call MLSMessageCalls( 'pop' )
    endif
  
    if ( timing ) call sayTime
    return

  contains

    subroutine Get_PCF_Id ( FileNameString, Path, SubString, L2Apriori_Version, &
      & FirstPCF, LastPCF, Description, GotFile, PCF_Id, ReturnStatus )
      use MLSFiles, only: GetPCFromRef, SPLIT_PATH_NAME
      use MLSL2Options, only: TOOLKIT
      use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS

      ! Args
      character(len=*), intent(inout) :: FileNameString
      character(len=*), intent(out) :: Path, SubString
      integer, intent(inout) :: L2Apriori_Version
      integer, intent(in) :: FirstPCF, LastPCF
      character(len=*), intent(in) :: Description
      logical, intent(in) :: GotFile
      integer, intent(inout) :: PCF_Id
      integer, intent(out) :: ReturnStatus
      ! Local variables
      integer :: PCFBottom
      ! Executable
      if ( gotFile ) then
        PCFBottom = FirstPCF
      else
        PCFBottom = PCF_Id + 1 ! This must be the last PCF id we found
      endif

      ! call outputNamedValue ( 'fileNameString', trim(fileNameString) )
      ! call outputNamedValue ( 'got file?', gotFile )
      ! call output('PCFBottom, lastPCF: ')
      ! call output( (/ PCFBottom, lastPCF /), advance='yes' )
      if ( TOOLKIT .and. gotFile ) then
        call split_path_name ( fileNameString, path, subString )
        pcf_id = getPCFromRef ( subString, PCFBottom, lastPCF, TOOLKIT, &
          &                     returnStatus, l2Apriori_Version, DEBUG, &
          &                     exactName = fileNameString )
      else if ( TOOLKIT ) then
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
    end subroutine Get_PCF_Id

    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for read_apriori = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine read_apriori

  ! ------------------------------------------  writeAPrioriAttributes_MF  -----
  subroutine writeAPrioriAttributes_MF ( MLSFile )
    type(MLSFile_T) :: MLSFile
    ! Executable
    if ( MLSFile%hdfVersion /= HDFVERSION_5 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Wrong hdfVersion--can write apriori attributes for hdf5 only', &
        & MLSFile=MLSFile )
      return ! Can only do this for hdf5 files
    elseif ( MLSFile%StillOpen ) then
      call writeAPrioriAttributes_ID(MLSFile%fileID%f_id, HDFVERSION_5)
    else
      call writeAPrioriAttributes_name(MLSFile%name, HDFVERSION_5)
    endif
  end subroutine writeAPrioriAttributes_MF

  ! ------------------------------------------  writeAPrioriAttributes_name  -----
  subroutine writeAPrioriAttributes_name ( fileName, hdfVersion )
    character(len=*), intent(in) :: fileName
    integer, intent(in)          :: hdfVersion  ! Must be 5 to work properly
    ! Internal variables
    integer             :: fileID
    integer             :: record_length
    integer             :: status
    ! Executable
    fileID = mls_io_gen_openF(l_swath, .TRUE., status, &
      & record_length, DFACC_RDWR, FileName=Filename, &
      & hdfVersion=hdfVersion, debugOption=.false. )  
    call writeAPrioriAttributes_ID ( fileID, hdfVersion )
    status = mls_io_gen_closeF(l_swath, fileID, &
      & hdfVersion=hdfVersion)
  end subroutine writeAPrioriAttributes_name

  ! ------------------------------------------  writeAPrioriAttributes_ID  -----
  subroutine writeAPrioriAttributes_ID ( fileID, hdfVersion )
    ! Write info about what apriori files were used
    ! Storing them as hdfeos5 attributes
    use HDFEOS5, only: MLS_charType
    use MLSHDFEOS, only: mls_EHwrglatt
    ! Args
    integer, intent(in) :: fileID
    integer, intent(in) :: hdfVersion  ! Must be 5 to work properly
    ! Internal variables
    integer             :: status
    ! Executable
    if ( hdfVersion /= HDFVERSION_5 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Wrong hdfVersion--can write apriori attributes for hdf5 only' )
      return ! Can only do this for hdf5 files
    endif
    status = mls_EHwrglatt(fileID, &
     & 'A Priori l2gp', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%l2gp))
    if ( status /= 0 ) &
      &  call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Problem writing APrioriFiles%l2gp' // trim(APrioriFiles%l2gp) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori l2aux', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%l2aux))
    if ( status /= 0 ) &
      &  call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Problem writing APrioriFiles%l2aux' // trim(APrioriFiles%l2aux) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori ncep', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%ncep))
    if ( status /= 0 ) &
      &  call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Problem writing APrioriFiles%ncep' // trim(APrioriFiles%ncep) )
    status = mls_EHwrglatt(fileID, &
     & 'A Priori gmao', MLS_CHARTYPE, 1, &
     &  trim(APrioriFiles%dao))
    if ( status /= 0 ) &
      &  call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Problem writing APrioriFiles%dao' // trim(APrioriFiles%dao) )
  end subroutine writeAPrioriAttributes_ID

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, l2_type, quantityName, &
    & hdfVersion, MLSFile )
    character(LEN=*), intent(in)   :: Name
    character(LEN=*), intent(in)   :: l2_type
    integer, optional,  intent(in) :: hdfVersion
    character(LEN=*), intent(in) :: quantityName
    type(MLSFile_T), optional :: MLSFile

    ! Local variables
    integer                        :: myhdfVersion
    call output ( 'Level 2 apriori product type : ' )
    call output ( trim(l2_type), advance='no')
    if ( present(hdfVersion) ) then
      call blanks(4)
      call output ( 'hdf ' )
      if ( hdfVersion == WILDCARDHDFVERSION ) then
        myhdfVersion = mls_hdf_version(trim(Name))
      else
        myhdfVersion = hdfVersion
      endif
      call output ( myhdfVersion, advance='yes')
    else
      call output ( ' ', advance='yes')
    endif
    call blanks(15)
    call output ( 'name : ' )
    call blanks(8)
    call output ( trim(Name), advance='yes')

    call output ( 'quantity', advance='yes')           
    call blanks(5)                                        
    call output ( trim(quantityName), advance='yes')      
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
 
    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
 
    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
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

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ReadAPriori

!=============================================================================

!
! $Log$
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
