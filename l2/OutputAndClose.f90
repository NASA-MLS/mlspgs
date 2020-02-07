! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================

module OutputAndClose ! outputs all data from the Join module to the
                      ! appropriate L2 Files

!==============================================================================

  use HDF, only: Dfacc_Rdonly, Dfacc_Rdwr
  use HighOutput, only: LetsDebug, OutputNamedValue
  use Machine, only: USleep
  use MLSFiles, only: HDFVersion_5, AddInitializeMLSFile, Dump, &
    & GetMLSFileByName, GetMLSFileByType, GetPCFromRef, &
    & MLS_CloseFile, MLS_Exists, MLS_Inqswath, MLS_OpenFile, &
    & MLS_SFStart, MLS_SFEnd, &
    & Split_Path_Name, UnSplitName
  use MLSHDF5, only: GetHDF5Attribute, MakeHDF5Attribute
  use MLSL2Options, only: Checkpaths, Default_HDFVersion_Write, L2cfnode, &
    & SpecialDumpFile, SkipDirectwrites, Toolkit, &
    & MLSL2Message, WriteFileAttributes
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Info, MLSMSG_Warning
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Gen, Toggle, Switches

  implicit none
  private
  public :: OUTPUT_CLOSE, ADD_METADATA

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! -----     Private declarations     ---------------------------------

  ! Should we get any of this from MLSL2Options?
  logical, parameter :: COPYGLOBALATTRIBUTES = .true.
  logical, parameter :: DGGFILEISHYBRID = .false.      ! may write PCF as DS?
  logical, parameter :: LOGFILEGETSMETADATA = .false.  ! metadata to log file?
  logical, parameter :: FAKEPARALLELMASTER =  .false.  ! make up fake stuff?
  integer, parameter :: MAXQUANTITIESPERFILE = 10000
  ! For Announce_Error
  integer :: ERROR

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  Output_Close  -----
  subroutine Output_Close ( root, L2GPDatabase, L2AUXDatabase, DirectDatabase, &
    & matrices, hessians, vectors, fileDataBase, GriddedDataBase, &
    & chunks, processingRange, canWriteL2PC )

    ! Hard-wired assumptions:

    ! ----------------------- metadata ------------------------

    !   for the l2aux the mcf is MLSPCF_MCF_L2DGM_Start
    !   for the log file the mcf is MLSPCF_MCF_L2LOG_Start
    !   for the dgg file the mcf is MLSPCF_MCF_L2DGG_Start

    ! The correspondence between MCF and l2gp files is determined by
    ! the value of        MCFFORL2GPOPTION
    ! (see write_metadata module for fuller explanation)

    use Allocate_Deallocate, only: Deallocate_Test, Allocate_Test
    use Chunks_M, only: MLSChunk_T, Dump
    use ChunkDivide_M, only: ChunkDivideConfig, Obstructions
    use DestroyCommand_M, only: DestroyCommand
    use DirectWrite_M, only: DirectData_T, Dump
    use DumpCommand_M, only: BooleanFromEmptySwath, BooleanFromFormula, &
      & DumpCommand, ExecuteCommand, MLSCase, MLSEndselect, MLSSelect, &
      & MLSSelecting, Skip
    use Expr_M, only: Expr
    use GriddedData, only: GriddedData_T
    use HessianModule_1, only: Hessian_T
    use HGrid, only: CreateHGridfromMLSCfinfo, DealWithObstructions
    use HGridsDatabase, only: HGrid_T, HGrids_T, &
      & AddHGridtoDatabase, Dump
    use Init_Tables_Module, only: F_Additional, F_AttrName, F_AttrValue, &
      & F_Destroy, F_Dontpack, F_File, F_HDFversion, &
      & F_MetaDataonly, F_Metaname, F_Moleculesecondderivatives, F_Overlaps, &
      & F_Packed, F_Quantities, F_Reset, F_Time, F_Type, F_WritecounterMAF, &
      & Field_First, Field_Last, &
      & L_L2aux, L_L2cf, L_L2dgg, L_L2gp, L_L2pc, &
      & S_Boolean, S_Case, S_Catenate, S_Copy, &
      & S_Destroy, S_Diff, S_Dump, S_Dumpblocks, S_Endselect, S_Execute, &
      & S_Hgrid, S_IsFileAbsent, S_Isswathempty, S_L2gp, S_Output, &
      & S_Reevaluate, S_Select, S_Skip, S_Sleep, S_Time, S_WriteFileAttribute
    use Intrinsic, only: Lit_Indices
    use L2AUXData, only: L2AUXData_T
    use L2GPData, only: L2GPData_T, &
      & AddL2GPToDatabase, WriteMastersFileAttributes
    use L2PC_M, only: OutputHDF5L2PC
    use L2ParInfo, only: Parallel
    use MatrixModule_1, only: Matrix_Database_T
    use MatrixTools, only: Dumpblocks
    use MLSCommon, only: MLSFile_T, Tai93_Range_T, Filenamelen
    use MLSL2Timings, only: Section_Times
    use MLSPCF2, only: MLSPCF_L2GP_End, &
      & MLSPCF_L2GP_Start, MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End
    use MLSStringlists, only: Switchdetail
    use MLSStrings, only: Trim_Safe
    use Moretree, only: Get_Label_And_Spec, Get_Spec_Id, Get_Boolean
    use Next_Tree_Node_M, only: Next_Tree_Node, Next_Tree_Node_State
    use Output_M, only: Output, RevertOutput, SwitchOutput
    use Time_M, only: SayTime, Time_Now
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decorate, Decoration, Nsons, Subtree, Sub_Rosa
    use VectorsModule, only: Vector_T
    use WriteMetaData, only: L2PCF, WriteMetalog

    ! Arguments
    integer, intent(in) :: ROOT   ! Of the output section's AST
    type (L2GPData_T), dimension(:), pointer :: L2GPDataBASE ! L2GP products
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDataBASE ! L2AUX products
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES ! Matrix database (for l2pcs)
    type (Hessian_T), dimension(:), pointer :: HESSIANS ! Hessian database (for l2pcs)
    type (Vector_T), dimension(:), pointer :: Vectors ! Vectors database
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (MLSChunk_T), dimension(:), pointer ::  Chunks  ! of data
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (GriddedData_T), dimension(:), pointer :: GriddedDataBase
    type (TAI93_Range_T), intent(in) :: processingRange

    logical, intent(in) :: canWriteL2PC ! Flag

    ! - - - Local declarations - - -

    type (MLSChunk_T) ::  AllChunks     ! in one
    ! logical, parameter :: DEBUG = .false.
    logical :: debug
    logical :: additional
    integer :: delay                    ! how many microseconds to sleep
    logical :: DESTROY
    integer, dimension(:), pointer :: DONTPACK ! Quantities not to pack
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    integer :: FIELDVALUE               ! For get_boolean
    character (len=FileNameLen) :: FILE_BASE    ! From the FILE= field
    logical, dimension(field_first:field_last) :: GOT ! Fields
    integer :: GSON                     ! Son of Son -- an assign node
    type (HGrids_T), dimension(:), pointer :: HGrids => null()
    integer :: HDFVERSION
    integer :: KEY                      ! Index of spec_args node
    integer :: Me = -1                  ! String index for trace
    integer :: Metadata_error
    character (len=32) :: meta_name     ! From the metaName= field
    integer :: NAME                     ! string index of label on output
    type (HGrid_T), pointer :: newHGridp
    integer :: NODE
    integer :: noGapsHGIndex = 0
    integer :: OUTPUT_TYPE              ! L_L2AUX, L_L2GP, L_PC, L_L2DGG
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    logical :: PACKED                   ! Do we pack this l2pc?
    integer :: QUANTITIESNODE           ! A tree node
    logical :: reset
    integer :: SECONDDERIVNODE
    integer :: SON                      ! Of Root -- spec_args or named node
    type(next_tree_node_state) :: State ! of tree traverser
    character (len=80) :: strValue
    real :: T1     ! for timing
    integer :: Units(2)                 ! Units of value returned by EXPR
    logical :: USINGL2Q                 ! Set if using the l2q queue manager
    logical :: USINGOLDSUBMIT              ! Set if using the submit mechanism
    logical :: USINGSUBMIT              ! Set if using the submit mechanism
    double precision :: Value(2)        ! Value returned by EXPR
    logical :: TIMING
    logical :: WriteCounterMAF          ! Add the counter MAF field
    logical :: writeMetaDataOnly

    ! Executable code

    call trace_begin ( me, "Output_Close", root, cond=toggle(gen) )
    timing = section_times
    if ( timing ) call time_now ( t1 )
    nullify ( dontPack )

    debug = LetsDebug ( 'output', 0 )
    error = 0

    if ( switchDetail( switches, 'pro') > -1 ) then
      call output ( '============ Level 2 Products ============', advance='yes' )
      call output ( ' ', advance='yes' )
    end if

    usingL2Q = ( index(parallel%submit, 'l2q') > 0 )
    usingSubmit = trim_safe(parallel%submit) /= ''
    usingOldSubmit = usingSubmit .and. .not. usingL2Q
    WRITEMASTERSFILEATTRIBUTES = .not. parallel%master ! cp file attributes from slaves

    AllChunks%firstMAFIndex = chunks(1)%firstMAFIndex
    AllChunks%noMAFsLowerOverlap = chunks(1)%noMAFsLowerOverlap
    AllChunks%lastMAFIndex = chunks(size(chunks))%lastMAFIndex
    AllChunks%noMAFsUpperOverlap = chunks(size(chunks))%noMAFsUpperOverlap

    if ( DEBUG ) then
      print *, 'Num chunks: ', size(chunks)
      print *, 'firstMAFIndex: ', AllChunks%firstMAFIndex
      print *, 'lastMAFIndex: ', AllChunks%lastMAFIndex
      call dump(AllChunks)
    endif
    ! Loop over the lines in the l2cf

    do 
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit

      meta_name = ''
      writeCounterMAF = .false.
      got = .false.
      additional = .false.

      call get_label_and_spec ( son, name, key )
      L2CFNODE = key

      if ( MLSSelecting .and. &
        & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) &
        & cycle

      reset = .false.
      do field_no = 2, nsons(key) ! fields of the "time" specification
        gson = subtree(field_no, key)   ! An assign node
        if ( nsons(gson) > 1 ) then
          fieldValue = decoration(subtree(2,gson)) ! The field's value
        else
          fieldValue = gson
        end if
        field_index = decoration(subtree(1,gson))
        got(field_index) = .true.
        select case ( field_index )   ! Field name
        case ( f_reset )
          reset = Get_Boolean ( gson )
        case default
          ! Shouldn't get here if the type checker worked
        end select
      end do ! j = 2, nsons(key)
      select case ( get_spec_id(key) )
      case ( s_Boolean )
        call decorate ( key,  BooleanFromFormula ( name, key ) )
      case ( s_Catenate )
        if ( associated(DirectDatabase) &
          & .and. .not. SKIPDIRECTWRITES ) then
          call output( ' unsplitting dgg/dgm files', advance='yes' )
          call unsplitFiles ( key, DirectDatabase, FileDatabase, HGrids, &
            & usingOldSubmit, debug )
        end if
      case ( s_sleep ) ! ============ Sleep ==========
        delay = 1000 ! defaults to 1000 microseconds
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          field_index = decoration(subtree(1,gson))
          got(field_index) = .true.
          select case ( field_index )   ! Field name
          case ( f_time )
            ! Did we say for how long?
            call expr ( subtree(2,gson), units, value )
            delay = value(1)*1.d6  ! Converted to microseconds
          case default
          end select 
        enddo
        call usleep ( delay )
      case ( s_select ) ! ============ Start of select .. case ==========
        ! We'll start seeking a matching case
        call MLSSelect (key)
      case ( s_case ) ! ============ seeking matching case ==========
        ! We'll continue seeking a match unless the case is TRUE
        call MLSCase (key)
      case ( s_endSelect ) ! ============ End of select .. case ==========
        ! We'done with seeking a match
        call MLSEndSelect (key)
      case ( s_reevaluate )
        call decorate ( key,  BooleanFromFormula ( 0, key ) )
      case ( s_skip ) ! ============================== Skip ==========
        ! We'll skip the rest of the section if the Boolean cond'n is TRUE
        if ( Skip(key) ) then
          call output( '(Skipping rest of this section)', advance='yes' )
          exit
        endif
      case ( s_diff, s_dump )
        call dumpCommand ( key, griddedDataBase=griddedDataBase, &
          & FiledataBase=FileDataBase, MatrixdataBase=matrices, &
          & Hessiandatabase=Hessians, HGrids=HGrids )
      case ( s_isFileAbsent, s_isSwathEmpty )
        if ( checkPaths ) cycle
        call decorate ( key, BooleanFromEmptySwath ( key ) )
      case ( s_copy )
        call copyQuantity( key, filedatabase )

      case ( s_Destroy )
        call destroyCommand ( key, matrices, vectors, griddedDataBase )
      case ( s_Dumpblocks )
        call dumpBlocks ( key, matrices, hessians )
      case ( s_execute ) ! ======================== ExecuteCommand ==========
        call ExecuteCommand ( key )
      case ( s_HGrid )
        if ( specialDumpFile /= ' ' ) &
          & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
        call decorate ( key, AddHGridToDatabase ( hGrids, CreateHGridFromMLSCFInfo ( name, key, filedatabase, L2GPDatabase, &
          & processingRange, allChunks ) ) )
        if ( DEBUG ) print *, 'Before dealing with obstructions'
        newHGridp => hGrids(1)%the_hGrid
        if ( DEBUG ) call dump(newHGridp)
        if ( associated(obstructions) ) &
          & call DealWithObstructions( newHGridp, obstructions, DestroyOld = .false. )
        ! Don't skip the lower overlap profiles if ChunkDivide included them
        if ( ChunkDivideConfig%allowPriorOverlaps ) &
          & newHGridp%noProfsLowerOverlap = 0
        if ( specialDumpFile /= ' ' ) &
          & call revertOutput

      case ( s_l2gp )
        call decorate ( key, AddL2GPToDatabase ( L2GPDatabase, &
          & CreateAndReadL2GP ( name, key, filedatabase ) ) )

      case ( s_output )
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          if ( nsons(gson) > 1 ) then
            fieldValue = decoration(subtree(2,gson)) ! The field's value
          else
            fieldValue = gson
          end if
          field_index = decoration(subtree(1,gson))
          got(field_index) = .true.
          select case ( field_index )   ! Field name
          case ( f_file )
            call get_string ( sub_rosa(subtree(2,gson)), file_base, strip=.true. )
          case ( f_metaName )
            call get_string ( sub_rosa(subtree(2,gson)), meta_name, strip=.true. )
          case ( f_type )
            output_type = decoration(subtree(2,gson))
            call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
          case ( f_writeCounterMAF )
            writeCounterMAF = get_boolean ( fieldValue )
          case ( f_MetaDataOnly )
            writeMetaDataOnly = get_boolean ( fieldValue )
          case ( f_hdfVersion )
            call expr ( subtree(2,gson), units, value )
            hdfVersion = nint(value(1))
            if ( hdfVersion /= 4 .and. hdfVersion /= 5 ) &
              & call Announce_error ( gson, 'hdfVersion must be 4 or 5')
          case default                  ! Everything else processed later
          end select
        end do

        ! Unless outputting l2cf, you must have supplied a quantities field
        ! (parser used to catch this)
        if ( .not. got(f_quantities) .and. output_type /= l_l2cf ) &
          & call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'No quantities field in Output command' )
        ! Otherwise--normal output commands
        select case ( output_type )
        case ( l_l2cf ) ! ------------------------------ Writing l2cf file ---
          call CopyTextFileToHDF ( file_base, DEBUG, filedatabase )
          return

        case ( l_l2gp ) ! --------------------- Writing l2gp files -----
        if ( noGapsHGIndex > 0 ) newHGridp => HGrids(noGapsHGIndex)%the_hGrid
          call OutputL2GP ( key, file_base, DEBUG, &
            & output_type, MLSPCF_L2GP_Start, MLSPCF_L2GP_End, &
            & filedatabase, L2GPDatabase, newHGridp )
          if ( .not. TOOLKIT ) return

        case ( l_l2aux ) ! ------------------------------ Writing l2aux files ---
          call OutputL2AUX ( key, file_base, DEBUG, writeCounterMAF, &
            & filedatabase, L2AUXDatabase )
          if ( .not. TOOLKIT ) return

        case ( l_l2pc ) ! ------------------------------ Writing l2pc files --
          if ( checkPaths ) exit   ! Not done on sips, so this is good enough
          ! I intend to completely ignore the PCF file in this case,
          ! it's not worth the effort!
          ! In that case, I will ignore the possibility of checkPaths being true
          if ( .not. canWriteL2PC ) call MLSL2Message( MLSMSG_Error, ModuleName, &
            & "Cannot write l2pc files with multi chunk l2cf's" )
          destroy = .false.
          packed = .false.
          secondDerivNode = 0
          do field_no = 2, nsons(key) ! Skip "output" name
            gson = subtree(field_no,key)
            select case ( decoration(subtree(1,gson)) )
            case ( f_destroy )
              destroy = get_boolean ( gson )
            case ( f_dontPack )
              call Allocate_Test ( dontPack, nsons(gson)-1, 'dontPack', ModuleName )
              do node = 2, nsons(gson)
                dontPack(node-1) = decoration(decoration(subtree(node,gson)))
              end do
            case ( f_moleculeSecondDerivatives )
              secondDerivNode = gson
            case ( f_overlaps )
              ! ??? More work needed here
            case ( f_packed )
              packed = get_boolean ( gson )
            case ( f_quantities )
              quantitiesNode = gson
            end select
          end do ! field_no = 2, nsons(key)

          call OutputHDF5L2PC ( trim(file_base), matrices, hessians, &
            & quantitiesNode, secondDerivNode, packed, dontPack )

          ! We used to destroy the written out matrix here, but I don't want to do that anymore
          ! (NJL, 13 Feb 2010)
          
          call Deallocate_test ( dontPack, 'dontPack', ModuleName )

        case ( l_l2dgg ) ! --------------------- Writing l2dgg files -----
        if ( noGapsHGIndex > 0 ) newHGridp => HGrids(noGapsHGIndex)%the_hGrid
          call OutputL2GP ( key, file_base, DEBUG, &
            & output_type, MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End, &
            & filedatabase, L2GPDatabase, newHGridp )
          if ( .not. TOOLKIT ) return

        case default
          if ( any(output_type /= (/ l_l2gp, l_l2dgg, l_l2aux, l_l2pc /)) ) then
            call announce_error ( son, &
            &  "Error--unknown output type: parser should have caught this")
          else
            call output('Lahey did weird thing again: ', advance='yes')
          end if
          call output('l2gp type number: ', advance='no')
          call output(l_l2gp, advance='yes')

          call output('l2aux type number: ', advance='no')
          call output(l_l2aux, advance='yes')

          call output('l2dgg type number: ', advance='no')
          call output(l_l2dgg, advance='yes')

          call output('output type number: ', advance='no')
          call output(output_type, advance='yes')

          call output('file_base: ', advance='no')
          call output(trim(file_base), advance='yes')

        end select

      case ( s_time )
        if ( timing .and. .not. reset ) then
          call sayTime ( 'Output_Close', t1=t1, cumulative=.false. )
        else
          call time_now ( t1 )
          timing = .true.
        end if

      case ( s_writeFileAttribute )
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          if ( nsons(gson) > 1 ) then
            fieldValue = decoration(subtree(2,gson)) ! The field's value
          else
            fieldValue = gson
          end if
          field_index = decoration(subtree(1,gson))
          got(field_index) = .true.
          select case ( field_index )   ! Field name
          case ( f_file )
            call get_string ( sub_rosa(subtree(2,gson)), file_base, strip=.true. )
          case ( f_additional )
            additional = get_boolean( gson )
          case ( f_attrName )
            call get_string ( sub_rosa(subtree(2,gson)), meta_name, strip=.true. )
          case ( f_attrValue )
            call get_string ( sub_rosa(subtree(2,gson)), strValue, strip=.true. )
          case ( f_type )
            output_type = decoration(subtree(2,gson))
            call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
          case default                  ! Everything else processed later
          end select
        end do
        if ( checkPaths ) cycle
        select case ( output_type )
        case ( l_l2aux ) ! --------------------- Writing it to l2aux files -----
          call writeAttributeToL2AUX ( file_base, meta_name, strValue, &
            & filedatabase, additional )       
        case ( l_l2gp ) ! --------------------- Writing it to l2gp files -----
          call writeAttributeToL2GP ( file_base, meta_name, strValue, &
            & output_type, MLSPCF_L2GP_Start, MLSPCF_L2GP_End, &
            & filedatabase, additional )       
        case ( l_l2dgg ) ! --------------------- Writing it to l2dgg files -----
          call writeAttributeToL2GP ( file_base, meta_name, strValue, &
            & output_type, MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End, &
            & filedatabase, additional )       
        case default
          call MLSL2Message( MLSMSG_Warning, ModuleName, &
          & 'Not yet able to write file attributes to this file type' )
        end select

      case default
        call announce_error ( son, &
          &  "Error--unknown son: parser should have caught this")

      end select

    end do
    
    ! Write the log file metadata
    if ( LOGFILEGETSMETADATA .and. .not. checkPaths ) then
      if ( DEBUG ) then
        call output('About to write log file metadata' , advance='yes')
      end if

      if ( TOOLKIT ) then
        call writeMetaLog ( metadata_error )
        ! error = max(error, metadata_error)
      end if
    end if

    ! Done wirh any Hgrids we may have created

    ! Done with text of PCF file at last
   
    if ( DEBUG ) &
      & call output ( 'About to deallocate text of PCF file' , advance='yes' )

    call deallocate_test ( l2pcf%anText, 'anText of PCF file', moduleName )

    if ( switchDetail(switches, 'pro') > -1 ) then
      call output ( '============ End Level 2 Products ============', advance='yes' )
      call output ( ' ', advance='yes' )
    end if

    if ( error /= 0 ) then
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Problem with Output_Close section' )
    end if

    if ( timing ) call sayTime ( 'Output_Close', cumulative=.false. )
    call trace_end ( "Output_Close", cond=toggle(gen) )

  end subroutine Output_Close

  ! ---------------------------------------------  CreateAndReadL2GP  -----
  function CreateAndReadL2GP ( name, key, filedatabase ) result(l2gp)
    ! Read and store l2gp data type from a file
    use Init_Tables_Module, only: F_File, F_Sdname, F_NoPCFid, F_Swath, &
      & Field_First, Field_Last
    use L2GPData, only: L2GPData_T
    use MLSCommon, only: FileNameLen, MLSFile_T
    use MLSFiles, only: AddFileToDatabase
    use MLSPCF2, only: &
      & MLSPCF_L2APriori_Start, MLSPCF_L2APriori_End
    use MoreTree, only: Get_Boolean
    use Readapriori, only: ProcessOneL2GPFile
    use Tree, only: Decoration, Nsons, &
      & Sub_Rosa, Subtree
    ! Args
    integer, intent(in)             :: name
    integer, intent(in)             :: key
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (L2GPData_T)               :: L2GP
    ! Internal variables
    integer                         :: field
    integer                         :: fieldIndex
    integer                         :: FileIndex
    integer                         :: FileName    ! Sub-rosa index of name in file='name'
    character(len=FileNameLen)      :: FileNameString ! actual literal file name
    logical, dimension(field_first:field_last) :: GOT
    integer                         :: j
    type (MLSFile_T)                :: L2GPFile
    integer                         :: LastAprioriPCF      ! l2gp or l2aux  apriori
    logical                         :: noPCFid
    integer                         :: SdName  ! sub-rosa index of name in sdName='name'
    character(len=FileNameLen)      :: ShortFileName
    integer                         :: son
    integer                         :: SwathName        ! sub-rosa index of name in swath='name'
    character(len=FileNameLen)      :: SWATHNAMESTRING ! actual literal swath name
    ! Executable
    LastAprioriPCF = MLSPCF_L2APriori_start - 1
    got = .false.
    noPCFid = .false.
    do j = 2, nsons(key)
      field = subtree(j,key)
      L2CFNODE = field
      fieldIndex = decoration(subtree(1,field))
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      case ( f_file )
        fileName = sub_rosa(subtree(2,field))
      case ( f_noPCFid )
        noPCFid = get_boolean(field)
      case ( f_sdname )
        sdname = sub_rosa(subtree(2,field))
      case ( f_swath )
        swathName = sub_rosa(subtree(2,field))
      end select
    end do
    if ( got(f_file) ) then
      call get_string ( FileName, fileNameString, strip=.true. )
    else
      fileNameString = ''
    end if
    shortFileName = fileNameString
    if ( .not. got(f_file) ) &
      & call announce_error ( son, &
        & 'Filename name must be specified to read an l2gp' )
    swathNameString=''
    if ( got(f_swath) ) &
      & call get_string ( swathName, swathNameString, strip=.true. )
    call processOneL2GPFile ( FileNameString, swathNameString, &
      & noPCFid, MLSPCF_L2APriori_start, MLSPCF_L2APriori_end, &
      & LastAprioriPCF, ' ', .false., &
      & L2GP, L2GPFile )
    FileIndex = AddFileToDataBase( filedatabase, L2GPFile )
  end function CreateAndReadL2GP
    
  ! ---------------------------------------------  add_metadata  -----
  subroutine add_metadata ( node, fileName, l2metaData, &
    & hdfVersion, filetype, metadata_error, &
    & numquantitiesperfileinput, quantityNamesInput )
    
    use Allocate_Deallocate, only: Deallocate_Test, Allocate_Test
    use Init_Tables_Module, only: L_L2dgg, L_Quantity
    use Intrinsic, only: L_Swath, L_HDF
    use MLSCommon, only: L2metaData_T
    use L2GPData, only: L2GPNamelen, Maxswathnamesbufsize
    use MLSHDF5, only: GetallHDF5dsnames
    use MLSPCF2, only: MLSPCF_L2DGM_End, MLSPCF_L2DGM_Start, MLSPCF_L2GP_End, &
      & MLSPCF_L2GP_Start, MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End, &
      & MLSPCF_MCF_L2DGM_Start, MLSPCF_MCF_L2DGG_Start, &
      & MLSPCF_MCF_L2GP_Start
    use MLSStringLists, only: GetHashElement, List2array, NumStringElements
    use Output_M, only: Output
    use WriteMetaData, only: L2PCF, Populate_MetaData_Std, &
      & Populate_MetaData_Oth, Get_L2gp_Mcf
  ! Deal with metadata--1st for direct write, but later for all cases
  integer, intent(in) :: node
  character(len=*), intent(in) :: fileName
  type (L2Metadata_T) :: l2metaData
  integer, intent(in) :: hdfVersion
  integer, intent(in) :: fileType  ! l_swath, l_hdf, ..
  integer, intent(out) :: metadata_error
  integer, optional, intent(in) :: numquantitiesperfileInput
  character(len=*), dimension(:), optional, intent(in) :: quantityNamesInput
  
  ! Internal variables
  integer :: baseIndex
  logical, parameter :: countEmpty = .true.
  ! logical, parameter :: DEBUG = .false.
  logical :: debug
  character(len=*), parameter :: L2GPHEAD = 'L2GP-'
  integer :: field_no
  character (len=132) :: FILE_BASE
  integer :: fileHandle
  integer :: L2aux_mcf
  integer :: l2dgg_mcf
  integer :: L2gp_mcf
  integer :: listSize
  character (len=32) :: meta_name=' '
  integer :: numquantitiesperfile
  character (len=132) :: path
  character(len=len(fileName)) :: PhysicalFilename
  character(len=L2GPNameLen), dimension(:), pointer :: quantityNames => null()
  integer :: returnStatus
  character(len=MAXSWATHNAMESBUFSIZE) :: sdList
  integer :: Version
  ! Executable
  debug = LetsDebug ( 'output', 0 )
  nullify(quantityNames)
  l2aux_mcf = MLSPCF_MCF_L2DGM_Start
  l2dgg_mcf = MLSPCF_MCF_L2DGG_Start
  l2gp_mcf  = MLSPCF_MCF_L2GP_Start - 1
  metadata_error = 0
  Version = 1
  select case (filetype)
  case (l_swath, l_l2dgg)
     call split_path_name(fileName, path, file_base)
     baseIndex = index(trim(file_base), L2GPHEAD)
     if ( baseIndex > 0 ) then
       file_base=file_base(baseIndex+len(L2GPHEAD):)
     end if
     if ( filetype == l_l2dgg ) then
       FileHandle = GetPCFromRef(file_base, MLSPCF_L2DGG_Start, &
         & MLSPCF_L2DGG_End, &
         & .true., returnStatus, Version, DEBUG, &
         & exactName=PhysicalFilename)
       l2gp_mcf = l2dgg_mcf
       call GetHashElement ( l2pcf%spec_keys, l2pcf%spec_doinames      , &
          & 'dgg', l2metaData%doiIdentifier, .TRUE. )
       if ( len_trim(l2metaData%doiIdentifier) < 1 ) &
         & l2metaData%doiIdentifier = '10.5067/AURA/MLS/DATA2006'
     else
       FileHandle = GetPCFromRef(file_base, MLSPCF_L2GP_Start, &
         & MLSPCF_L2GP_End, &
         & .true., returnStatus, Version, DEBUG, &
         & exactName=PhysicalFilename)
       call get_l2gp_mcf ( file_base, meta_name, l2gp_mcf, &
         & l2metaData%doiIdentifier )
     end if
     if (returnStatus /= 0) then
         call MLSL2Message ( MLSMSG_Error, ModuleName, &
           &  "While adding metadata failed to GetPCFromRef for " // trim(fileName) )
     else if ( l2gp_mcf <= -999 ) then
         call MLSL2Message ( MLSMSG_Warning, ModuleName, &
           &  "no mcf for this l2gp species in" // trim(file_base) )
         return
     else if (l2gp_mcf <= 0) then
         call MLSL2Message ( MLSMSG_Error, ModuleName, &
           &  "no mcf for this l2gp species in" // trim(file_base) )
     end if
     
     ! What quantities do we want metadata for?
     if ( present(quantityNamesInput) .and. present(numquantitiesperfileInput) ) then
       call allocate_test( quantityNames, size(quantityNamesInput), &
         & 'quantityNames', ModuleName )
       numquantitiesperfile = numquantitiesperfileInput
       quantityNames = quantityNamesInput
     else
       numquantitiesperfile = mls_InqSwath ( PhysicalFilename, sdList, listSize, &
           & hdfVersion=hdfVersion)
       call allocate_test(quantityNames, numquantitiesperfile, &
         & 'quantityNames', ModuleName )
       call List2Array( sdList, quantityNames, countEmpty )
     endif

     if ( QuantityNames(numquantitiesperfile) &
       & == QuantityNames(1) ) then
       ! Typical homogeneous l2gp file: 
       ! e.g., associated with BrO is ML2BRO.001.MCF
       if ( DEBUG ) then
         call output('preparing to populate metadata_std', advance='yes')
         call output('l2gpFileHandle: ', advance='no')
         call output(FileHandle , advance='no')
         call output('   l2gp_mcf: ', advance='no')
         call output(l2gp_mcf , advance='yes')
       end if
       call populate_metadata_std &
         & (FileHandle, l2gp_mcf, QuantityNames(1), l2metaData, &
         & hdfVersion=hdfVersion, metadata_error=metadata_error, &
         & filetype=filetype  )
     else
       ! Type l2gp file 'other'
       if ( DEBUG ) then
         call output ( 'preparing to populate metadata_oth', advance='yes' )
         call output ( 'l2gpFileHandle: ', advance='no' )
         call output ( FileHandle , advance='no' )
         call output ( '   l2gp_mcf: ', advance='no' )
         call output ( l2gp_mcf , advance='yes' )
       end if

       call populate_metadata_oth &
         & ( FileHandle, l2gp_mcf, &
         & numquantitiesperfile, QuantityNames, l2metaData, &
         & hdfVersion=hdfVersion, metadata_error=metadata_error, &
         & filetype=filetype  )
     end if
  case (l_hdf)
     if ( DEBUG ) call output ( 'output file type l2aux', advance='yes' )
     ! Get the l2aux file name from the PCF

     call split_path_name(fileName, path, file_base)
     FileHandle = GetPCFromRef(file_base, MLSPCF_L2DGM_Start, &
       & MLSPCF_L2DGM_End, &
       & .true., returnStatus, Version, DEBUG, &
       & exactName=PhysicalFilename)
     ! Some hdf-formatted product files lack metadata, e.g. fwm radiances
     if (returnStatus /= 0) then
         metadata_error = returnStatus
         call outputNamedValue ( 'fileName', trim(fileName) )
         call outputNamedValue ( 'PCFIds', (/ MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End /) )
         call MLSL2Message ( MLSMSG_Warning, ModuleName, &
           &  "While adding metadata failed to GetPCFromRef for " // trim(fileName) )
         return
     end if
     if ( present(quantityNamesInput) .and. present(numquantitiesperfileInput) ) then
       call allocate_test( quantityNames, size(quantityNamesInput), &
         & 'quantityNames', ModuleName )
       numquantitiesperfile = numquantitiesperfileInput
       quantityNames = quantityNamesInput
     else
       call GetAllHDF5DSNames ( PhysicalFilename, '/', sdList )
       numquantitiesperfile = NumStringElements ( sdList, countEmpty )
       call allocate_test( quantityNames, numquantitiesperfile, &
         & 'quantityNames', ModuleName )
       call List2Array( sdList, quantityNames, countEmpty )
     endif
     if ( DEBUG ) then
       call output ( 'preparing to populate metadata_oth', advance='yes' )
       call output ( 'l2auxFileHandle: ', advance='no' )
       call output ( FileHandle , advance='no' )
       call output ( '   l2aux_mcf: ', advance='no' )
       call output ( l2aux_mcf , advance='no' )
       call output ( '   number of quantities: ', advance='no' )
       call output ( numquantitiesperfile , advance='yes' )
       if ( associated(QuantityNames) ) then
         do field_no=1, numquantitiesperfile
           call output ( field_no , advance='no' )
           call output ( '       ', advance='no' )
           call output ( trim(QuantityNames(field_no)) , advance='yes' )
         end do
       endif
     endif
     call GetHashElement ( l2pcf%spec_keys, l2pcf%spec_doinames      , &
          & 'dgm', l2metaData%doiIdentifier, .TRUE. )
     if ( len_trim(l2metaData%doiIdentifier) < 1 ) &
       & l2metaData%doiIdentifier = '10.5067/AURA/MLS/DATA2007'
     call populate_metadata_oth &
       & ( FileHandle, l2aux_mcf, &
       & numquantitiesperfile, QuantityNames, l2metaData, &
       & hdfVersion=hdfVersion, metadata_error=metadata_error, &
       & filetype=filetype  )
  case ( l_quantity )
    call MLSL2Message( MLSMSG_Warning, ModuleName, &
      & "Unable to add metadata for the file type of "&
      & // trim(filename) )
  case default
    call announce_error ( node, &
      &  "Error--filetype unrecognized", filetype)
    call MLSL2Message( MLSMSG_Error, ModuleName, &
      & "Unrecognized filetype in add_metadata (must be swath or hdf) "&
      & // trim(filename) )
  end select

  call deallocate_test( quantityNames, 'quantityNames', ModuleName )
  end subroutine add_metadata

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, l2_type, num_quants, quantities, hdfVersion )
    use Output_M, only: Blanks, Output
    integer, intent(in) :: num_quants 
    character(LEN=*), intent(in)   :: Name
    character(LEN=*), intent(in)   :: l2_type
    integer, optional,  intent(in) :: hdfVersion
    character(LEN=*), dimension(:), optional, intent(in) :: quantities
    integer :: i

    call output ( 'Level 2 output product type : ' )
    call output ( trim(l2_type), advance='no')
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

    if ( num_quants > 0 .and. present(quantities) ) then
      call output ( 'number ' )
      call blanks(5)
      call output ( 'quantity', advance='yes')
      do i=1, num_quants
          call output ( i )
          call blanks(5)
          call output ( trim(quantities(i)), advance='yes')
      end do
    end  if
  end subroutine announce_success

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )

    use Lexer_Core, only: Print_Source
    use Output_M, only: Output
    use Tree, only: Where_At => Where

    integer, intent(in) :: Where   ! Tree node where error was noticed
    character(LEN=*), intent(in) :: Full_Message
    integer, intent(in), optional :: Code    ! Code for error message
    integer, intent(in), optional :: Penalty
    integer :: myPenalty

    myPenalty = 1
    if ( present(penalty) ) myPenalty = penalty
    error = max(error,myPenalty)

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( where_at(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ' OutputAndClose complained: ' )

    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      ! select case ( code )
      ! end select
      call output ( ' Code: ' )
      call output ( Code, advance='yes' )
    end  if
  end subroutine ANNOUNCE_ERROR

  ! ---------------------------------------------  CopyQuantity  -----
  subroutine CopyQuantity ( key, fileDatabase )
    ! Do the work of copying named quantity data to a named file
    use Expr_M, only: Expr
    use HGridsDatabase, only: HGrids_T, Dump
    use Init_Tables_Module, only: F_Create, &
      & F_Exclude, F_File, F_HDFversion, F_Hgrid, &
      & F_Ifanycrashedchunks, F_InputFile, F_Inputtype, &
      & F_Options, &
      & F_Rename, F_Repairgeolocations, &
      & F_Swath, F_Toattribute, F_Type, &
      & Field_First, Field_Last, &
      & L_L2aux, L_L2cf, L_L2dgg, L_L2gp
    use HighOutput, only: OutputnamedValue
    use Intrinsic, only: L_Ascii, L_Swath, L_HDF, Lit_Indices
    use L2AUXData, only: CpL2AUXData
    use L2GPData, only: Avoidunlimiteddims, &
      & Maxswathnamesbufsize, CpL2GPData, CpL2GPDatatoattribute
    use L2parinfo, only: Parallel
    use MLSCommon, only: MLSFile_T, Filenamelen, L2metaData_T
    use MLSPCF2, only: MLSPCF_L2DGM_End, MLSPCF_L2DGM_Start, MLSPCF_L2GP_End, &
      & MLSPCF_L2GP_Start, MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End, &
      & MLSPCF_L2ascii_Start, MLSPCF_L2ascii_End
    use MLSStringlists, only: Intersection, Switchdetail
    use MLSStrings, only: Lowercase
    use Moretree, only: Get_Boolean
    use Output_M, only: Output
    use PCFHdr, only: Globalattributes, DumpGlobalAttributes, &
      & H5_Writeglobalattr, He5_WriteMLSFileattr, He5_WriteGlobalAttr
    use ReadAPriori, only: ReadAPrioriattributes, WriteAPrioriAttributes
    use Tree, only: Decoration, Nsons, Subtree, Sub_Rosa
    ! Args
    integer, intent(in)                       :: KEY
    type (MLSFile_T), dimension(:), pointer   :: FILEDATABASE
    ! Local variables
    logical :: create
    ! logical, parameter :: DEBUG = .false.
    logical :: debug
    character (len=MAXSWATHNAMESBUFSIZE) :: EXCLUDE   ! From the exclude= field
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    integer :: FIELDVALUE               ! For get_boolean
    character (len=FileNameLen) :: FILE_BASE    ! From the FILE= field
    integer :: FORMATTYPE               ! l_hdf or l_swath
    logical, dimension(field_first:field_last) :: GOT ! Fields
    type (HGrids_T), dimension(:), pointer :: HGrids => null()
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: HGridIndex
    integer :: INPUT_TYPE              ! L_L2AUX, L_L2GP, L_PC, L_L2DGG
    type(MLSFile_T), pointer :: inputFile
    character (len=FileNameLen) :: INPUTFILE_BASE    ! From the inputfile= field
    character (len=FileNameLen) :: inputPhysicalFilename
    character (len=MAXSWATHNAMESBUFSIZE) :: inSwathList
    character (len=MAXSWATHNAMESBUFSIZE) :: outSwathList
    type (L2Metadata_T) :: l2metaData ! L2GP metadata 
    integer :: listSize
    integer :: Metadata_error
    integer :: numSwaths
    logical :: newFile                  ! is this file new?
    integer :: noGapsHGIndex = 0
    integer :: noSwaths
    character(len=8) :: optionsString   ! e.g. '-f'
    integer :: OUTPUT_TYPE              ! L_L2AUX, L_L2GP, L_PC, L_L2DGG
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=FileNameLen) :: PhysicalFilename
    character (len=MAXSWATHNAMESBUFSIZE) :: rename
    logical :: RepairGeoLocations
    character (len=MAXSWATHNAMESBUFSIZE) :: sdList
    character (len=MAXSWATHNAMESBUFSIZE) :: sdListThere
    logical :: skipCopy
    integer :: SON                      ! Of Root -- spec_args or named node
    logical :: toAttribute
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR
    ! Executable
    debug = LetsDebug ( 'output', 0 )
    hdfVersion = DEFAULT_HDFVERSION_WRITE
    exclude = ''
    sdList = '*' ! This is wildcard meaning 'every sd or swath'
    rename = ' ' ! This is a blank meaning 'Dont rename the swaths'
    optionsString = ' ' ! This is a blank meaning 'Dont rename the swaths'
    got = .false.
    repairGeoLocations = .false.
    skipCopy = .false.
    create = .false.
    newFile = .false.
    toAttribute = .false.
    do field_no = 2, nsons(key)       ! Skip the command name
      gson = subtree(field_no, key)   ! An assign node
      if ( nsons(gson) > 1 ) then
        fieldValue = decoration(subtree(2,gson)) ! The field's value
      else
        fieldValue = gson
      end if
      field_index = decoration(subtree(1,gson))
      got(field_index) = .true.
      select case ( field_index )   ! Field name
      case ( f_create )
        create = get_boolean ( gson )
      case ( f_exclude )
        call get_string ( sub_rosa(subtree(2,gson)), exclude, strip=.true. )
      case ( f_file )
        call get_string ( sub_rosa(subtree(2,gson)), file_base, strip=.true. )
        if ( DEBUG ) print *, 'file_base: ', trim(file_base)
      case ( f_hdfVersion )
        call expr ( subtree(2,gson), units, value )
        hdfVersion = value(1)
        if ( hdfVersion /= 4 .and. hdfVersion /= 5 ) &
          & call Announce_error ( gson, 'hdfVersion must be 4 or 5')
      case ( f_hgrid )
        HGridIndex = decoration(fieldValue)
        if ( DEBUG ) print *, 'HGridIndex: ', HGridIndex
      case ( f_ifAnyCrashedChunks )
        skipCopy = get_boolean ( gson ) .and. &
          & ( parallel%numFailedChunks == 0 )
      case ( f_inputfile )
        call get_string ( sub_rosa(subtree(2, gson)), inputfile_base, strip=.true. )
      case ( f_inputtype )
        input_type = decoration(subtree(2, gson))
      case ( f_options )
        call get_string ( sub_rosa(subtree(2, gson)), optionsString, strip=.true. )
        optionsString = lowerCase(optionsString)
        if (switchDetail( switches, 'pro') > 0 ) &
          & call outputNamedValue( 'options', trim(optionsString) )
      case ( f_rename )
        call get_string ( sub_rosa(subtree(2,gson)), rename, strip=.true. )
      case ( f_repairGeoLocations )
        repairGeoLocations = get_boolean ( gson )
      case ( f_swath )
        call get_string ( sub_rosa(subtree(2,gson)), sdList, strip=.true. )
      case ( f_toAttribute )
        toAttribute = get_boolean ( gson )
      case ( f_type )
        output_type = decoration(subtree(2, gson))
        call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
      case default                  ! Parser should have caught this
      end select
    end do
    ! Make certain we have everything we need
    if ( ( got(f_repairGeoLocations) .and. .not. got(f_hgrid) ) ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & "Cannot repair Geolocs w/o an HGrid to copy from" )
    if ( .not. associated(HGrids) ) then
      if ( got(f_hgrid) ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & "No HGrids defined yet" )
    else
      if ( ( size(HGrids) < 1 .and. got(f_hgrid) ) ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & "No HGrids defined yet" )
    endif
    if ( ( got(f_swath) .and. got(f_exclude) ) ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & "Cannot copy specifying both swaths and excludes" )
    ! To skip or not to skip
    if ( skipCopy ) then
      call MLSL2Message ( MLSMSG_Info, ModuleName, &
      & "No crashed chunks so skipping this copy" )
      return
    endif
    if ( .not. got(f_inputtype) ) input_type = output_type

    select case ( output_type )
    case ( l_l2aux ) ! --------------------- Copying l2aux files -----
      call returnFullFileName(file_base, PhysicalFilename, &
        & MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End)
      outputFile => GetMLSFileByName(filedatabase, PhysicalFilename)
      if ( .not. associated(outputFile) ) then
        newFile = .true.
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=PhysicalFilename, shortName=file_base, &
          & type=l_hdf, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2DGM_Start, PCTop=MLSPCF_L2DGM_End)
      endif
    case ( l_l2gp ) ! --------------------- Copying l2gp files -----
      call returnFullFileName(file_base, PhysicalFilename, &
        & MLSPCF_L2GP_Start, MLSPCF_L2GP_End)
      outputFile => GetMLSFileByName(filedatabase, PhysicalFilename)
      if ( .not. associated(outputFile) ) then
        newFile = .true.
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=PhysicalFilename, shortName=file_base, &
          & type=l_swath, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2GP_Start, PCTop=MLSPCF_L2GP_End)
      endif
    case ( l_l2dgg ) ! --------------------- Copying l2dgg files -----
      call returnFullFileName(file_base, PhysicalFilename, &
        & MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End)
      outputFile => GetMLSFileByName(filedatabase, PhysicalFilename)
      if ( .not. associated(outputFile) ) then
        newFile = .true.
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=PhysicalFilename, shortName=file_base, &
          & type=l_swath, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2DGG_Start, PCTop=MLSPCF_L2DGG_End)
      endif
    case default
    end select

    if ( .not. associated(outputFile) ) then
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to Copy to ' // trim(PhysicalFilename) )
    endif
    if ( DEBUG ) call dump( outputFile, details=2 )

    select case ( input_type )
    case ( l_l2aux ) ! --------------------- Copying l2aux files -----
      call returnFullFileName(inputfile_base, inputPhysicalFilename, &
        & MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End)
      inputFile => GetMLSFileByName(filedatabase, inputPhysicalFilename)
      if ( .not. associated(inputFile) ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'No entry in filedatabase for ' // trim(inputPhysicalFilename) )
      endif
    case ( l_l2gp ) ! --------------------- Copying l2gp files -----
      call returnFullFileName(inputfile_base, inputPhysicalFilename, &
        & MLSPCF_L2GP_Start, MLSPCF_L2GP_End)
      inputFile => GetMLSFileByName(filedatabase, inputPhysicalFilename)
      if ( .not. associated(inputFile) ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'No entry in filedatabase for ' // trim(inputPhysicalFilename) )
      endif
    case ( l_l2dgg ) ! --------------------- Copying l2dgg files -----
      call returnFullFileName(inputfile_base, inputPhysicalFilename, &
        & MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End)
      inputFile => GetMLSFileByName(filedatabase, inputPhysicalFilename)
      if ( .not. associated(inputFile) ) then
        call output( ' output base ', advance='no' )
        call output( trim(file_base), advance='yes' )
        call output( ' output file ', advance='no' )
        call output( trim(outputFile%name), advance='yes' )
        call output( ' input base ', advance='no' )
        call output( trim(inputfile_base), advance='yes' )
        call output( ' input file ', advance='no' )
        call output( '(not associated)', advance='yes' )
        call dump(filedatabase)
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'No entry in filedatabase for ' // trim(inputPhysicalFilename) )
      endif
    case ( l_ascii, l_l2cf ) ! --------------------- Copying ascii files -----
      call returnFullFileName(inputfile_base, inputPhysicalFilename, &
        & MLSPCF_l2ascii_start, MLSPCF_l2ascii_end)
      inputFile => AddInitializeMLSFile( filedatabase, &
          & content=trim(inputfile_base), &
          & name=inputPhysicalFilename, shortName=inputfile_base, &
          & type=l_ascii, access=DFACC_RDWR, &
          & PCBottom=MLSPCF_l2ascii_start, PCTop=MLSPCF_l2ascii_end )
      if ( .not. associated(inputFile) ) then
        call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(inputPhysicalFilename) )
        return
      endif
      if ( DEBUG ) call dump( inputFile )
    case default
    end select

    if ( DEBUG ) print *, 'inputPhysicalfileName: ', trim(inputPhysicalfileName)
    if ( DEBUG ) print *, 'outputPhysicalfileName: ', trim(PhysicalfileName)
    if ( DEBUG ) print *, 'repairGeoLocations: ', repairGeoLocations
    if ( DEBUG ) print *, 'create: ', create
    if ( DEBUG ) call display_string ( lit_indices(output_type), &
    &             strip=.true., advance='yes' )
    if ( DEBUG ) call display_string ( lit_indices(input_type), &
    &             strip=.true., advance='yes' )

    if ( CHECKPATHS ) return

    if ( DEBUG ) call dump( inputFile, details=2 )
    if ( repairGeoLocations ) optionsString = trim(optionsString) // 'cr'
    select case ( output_type )
    case ( l_l2aux ) ! --------------------- Copying to l2aux files -----
      formattype = l_hdf
      ! Note that we haven't yet implemented repair stuff for l2aux
      ! So crashed chunks remain crashed chunks
      if ( DEBUG ) call Dump( inputFile )
      if ( input_type /= output_type ) then
        ! So far we only allow copying text files
        if ( .not. any( input_type == (/ l_ascii, l_l2cf /) ) ) then
          call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Unable to copy to l2aux file file type for' // trim(inputPhysicalFilename) )
          return
        endif
        if ( mls_exists( inputFile%Name ) == 0 ) &
          & call CopyTextFileToHDF ( file_base, DEBUG, filedatabase, inputFile )
      else
        call cpL2AUXData( inputFile, &
          & outputFile, create2=create, &
          & sdList=trim(sdList) )
      endif
    case ( l_l2gp, l_l2dgg ) ! --------------------- Copying l2gp files -----
      formattype = l_swath
      ! How to use swaths field? See definition of sdList
      ! Before trying to cp these swaths, make sure they're actually there
      noSwaths = mls_InqSwath ( inputFile%name, sdListThere, listSize, &
       & hdfVersion=inputFile%hdfVersion )
      sdList = Intersection( sdList, sdListThere )
      if ( sdList == ' ' ) return
      if ( toAttribute ) then
        if ( len_trim(rename) < 1 ) rename = sdlist
        call cpL2GPDataToAttribute( inputfile, outputfile, &
          & sdList, rename )
      elseif ( got(f_exclude) .and. .not. repairGeoLocations ) then
        call cpL2GPData( l2metaData, inputfile, &
          & outputfile, create2=create, &
          & exclude=trim(exclude), &
          & notUnlimited=avoidUnlimitedDims, &
          & andGlAttributes=COPYGLOBALATTRIBUTES, options=optionsString )
      elseif ( .not. got(f_exclude) .and. repairGeoLocations ) then
        if ( DEBUG ) print *,' size(filedatabse) ', size(filedatabase)
        if ( DEBUG ) print *, 'input file: ', trim(inputPhysicalFilename)
        if ( DEBUG ) print *, 'output file: ', trim(PhysicalFilename)
        call cpL2GPData( l2metaData, inputfile, &
          & outputfile, create2=create, &
          & swathList=trim(sdList), rename=rename, &
          & notUnlimited=avoidUnlimitedDims, andGlAttributes=COPYGLOBALATTRIBUTES, &
          & HGrid=HGrids(HGridIndex)%the_hGrid, options=optionsString )
      elseif ( got(f_exclude) .and. repairGeoLocations ) then
        call cpL2GPData( l2metaData, inputfile, &
          & outputfile, create2=create, &
          & swathList=trim(sdList), rename=rename, &
          & exclude=trim(exclude), &
          & notUnlimited=avoidUnlimitedDims, andGlAttributes=COPYGLOBALATTRIBUTES, &
          & HGrid=HGrids(HGridIndex)%the_hGrid, options=optionsString )
      else
        if ( DEBUG ) then
          call output( 'input file', advance='yes' )
          call dump( inputFile )
          call output( 'output file', advance='yes' )
          call dump( outputFile )
          numSwaths = mls_InqSwath ( inputfile%name, inSwathList, listSize, &
             & hdfVersion=hdfVersion )
          call outputNamedValue ( 'input swaths', trim(inSwathList) )
          numSwaths = mls_InqSwath ( outputfile%name, outSwathList, listSize, &
             & hdfVersion=hdfVersion )
          call outputNamedValue ( 'output swaths', trim(outSwathList) )
          numSwaths = mls_InqSwath ( PhysicalFilename, outSwathList, listSize, &
             & hdfVersion=hdfVersion )
          call outputNamedValue ( 'output swaths', trim(outSwathList) )
          call outputNamedValue ( 'swath', trim(sdList) )
        endif
        call cpL2GPData( l2metaData, inputfile, &
          & outputfile, create2=create, &
          & swathList=trim(sdList), rename=rename, &
          & notUnlimited=avoidUnlimitedDims, &
          & andGlAttributes=COPYGLOBALATTRIBUTES, options=optionsString )
      endif
      if ( WRITEFILEATTRIBUTES ) call he5_writeMLSFileAttr( outputFile )
      if ( create ) then
        call readAPrioriAttributes( inputFile )
        call writeAPrioriAttributes( outputFile )
      endif
      if ( noGapsHGIndex > 0 ) &
        & call writeHGridComponents( trim(PhysicalFilename), &
        & HGrids(noGapsHGIndex)%the_hGrid )
    case default
    end select

    if ( TOOLKIT .and. newFile ) then
      call add_metadata ( son, file_base, l2metaData, &
        & outputfile%hdfVersion, formattype, metadata_error )
      if ( len_trim(l2metaData%doiIdentifier) < 1 ) then
        call MLSL2Message( MLSMSG_Warning, ModuleName, &
          & 'empty doiIdentidier for ' // trim(PhysicalFilename) )
        return
      endif
      ! Set the file-level attribute DOI to its metadata value
      GlobalAttributes%DOI = l2metaData%doiIdentifier
      select case ( output_type )
      case ( l_l2aux ) ! --------------------- Copying to l2aux files -----
        call output( 'Writing file level attributes to h5 ' // &
          & trim(PhysicalFilename), advance='yes' )
        call outputNamedValue( 'identifier_product_DOI', trim(GlobalAttributes%DOI) )
        call outputNamedValue( 'ProductionLocation', trim(GlobalAttributes%productionLoc) )
        call h5_writeGlobalAttr ( outputFile, &
          & skip_if_already_there=.false., doi=.true. )
        call DumpGlobalAttributes
      case ( l_l2gp, l_l2dgg ) ! --------------------- Copying l2gp files -----
        call output( 'Writing file level attributes to he5 ' // &
          & trim(PhysicalFilename), advance='yes' )
        call outputNamedValue( 'identifier_product_doi', trim(GlobalAttributes%DOI) )
        call outputNamedValue( 'ProductionLocation', trim(GlobalAttributes%productionLoc) )
        call he5_writeglobalattr( outputFile, doi=.true. )
        call DumpGlobalAttributes
      end select
    endif
  end subroutine CopyQuantity

  ! ---------------------------------------------  writeAttributeToL2AUX  -----
  subroutine writeAttributeToL2AUX ( fileName, attrName, attrValue, &
    & filedatabase, additional )
    ! Do the work of outputting specified attribute to a named l2aux file
    use HDF5, only: H5GClose_F, H5GOpen_F
    use Intrinsic, only: L_HDF
    use MLSCommon, only: MLSFile_T, FileNameLen
    use MLSPCF2, only: MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End
    use MLSStringLists, only: SwitchDetail
    use Output_M, only: Output
    ! Args
    character(len=*), intent(inout)         :: fileName ! according to l2cf
    character(len=*), intent(in)            :: attrName ! according to l2cf
    character(len=*), intent(in)            :: attrValue ! according to l2cf
    type(MLSFile_T), dimension(:), pointer  :: filedatabase
    logical, intent(in)                     :: additional
    ! Local variables
    ! logical, parameter :: DEBUG = .false.
    logical :: debug
    character (len=FileNameLen) :: FullFilename
    integer :: FileHandle
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=132) :: path
    integer :: ReturnStatus
    logical :: skipIfAlreadyThere
    character (len=132) :: str
    integer :: Version

    ! Executable
    debug = LetsDebug ( 'output', 0 )
    Version = 1
    OUTPUTTYPESTR = 'l2aux'
    ! Get the l2aux file name from the PCF

    if ( TOOLKIT ) then
      call split_path_name(fileName, path, fileName)

      FileHandle = GetPCFromRef(fileName, MLSPCF_L2DGM_Start, &
        & MLSPCF_L2DGM_End, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullFilename)
    else
      FullFilename = fileName
      returnStatus = 0
    end if

    if ( mls_exists(trim(FullFilename)) /= 0 ) then
      if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
        & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Its non-existence prevented writing attribute to ' &
        & // trim(FullFilename) )
      return
    elseif ( returnStatus == 0 .and. .not. checkPaths ) then
      ! Open the HDF file and write l2aux data
      outputFile => GetMLSFileByName(filedatabase, FullFilename)
      if ( .not. associated(outputFile) ) then
      if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
          & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(FullFilename) )
        outputFile => AddInitializeMLSFile( filedatabase, &
          & content=outputTypeStr, &
          & name=FullFilename, shortName=fileName, &
          & type=l_hdf, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2DGM_Start, PCTop=MLSPCF_L2DGM_End )
      endif
    endif
    if ( .not. outputFile%stillOpen ) &
      & call mls_openFile( outputFile, returnStatus )
    call h5gopen_f( outputFile%FileID%f_id, '/', outputFile%FileID%grp_id, &
      & returnStatus )
    if ( .not. additional ) then
      str = attrValue
      skipIfAlreadyThere = .true.
    else
      call GetHDF5Attribute ( outputFile%FileID%grp_id, attrname, str )
      str = trim(str) // attrValue
      skipIfAlreadyThere = .false.
      if ( debug ) call output( 'Updated glob attr to ' // trim(str), advance='yes' )
    endif
    call MakeHDF5Attribute( outputFile%FileID%grp_id, &
      & attrName, str, skip_if_already_there=skipIfAlreadyThere )
    call h5gclose_f( outputFile%FileID%grp_id, returnStatus )
  end subroutine writeAttributeToL2AUX

  ! ---------------------------------------------  writeAttributeToL2GP  -----
  subroutine writeAttributeToL2GP ( fileName, attrName, attrValue, &
    & output_type, pcf_start, pcf_end, &
    & filedatabase, additional )
    ! Do the work of outputting specified attribute to a named l2gp file
    use HDFEOS5, only: MLS_Chartype
    use HDF5, only: H5f_Acc_Rdwr_F, H5fclose_F, H5fopen_F, H5gclose_F, H5gopen_F
    use Intrinsic, only: L_Swath, Lit_Indices
    use MLSCommon, only: MLSFile_T, FileNameLen
    use MLSHDF5, only: MakeHDF5Attribute
    use MLSHDFEOS, only: HE5_EHRDGlAtt, MLS_EHWRGlatt
    use MLSStringLists, only: SwitchDetail
    use Output_M, only: Output
    ! Args
    character(len=*), intent(inout)         :: fileName ! according to l2cf
    character(len=*), intent(in)            :: attrName ! according to l2cf
    character(len=*), intent(in)            :: attrValue ! according to l2cf
    integer, intent(in)                     :: output_type
    integer, intent(in)                     :: pcf_start
    integer, intent(in)                     :: pcf_end
    type(MLSFile_T), dimension(:), pointer  :: filedatabase
    logical, intent(in)                     :: additional

    ! Local variables
    ! logical, parameter :: DEBUG = .false.
    logical :: debug
    character (len=FileNameLen) :: FullFilename
    integer :: FileHandle
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=132) :: path
    integer :: ReturnStatus
    ! logical :: skipIfAlreadyThere
    character (len=132) :: str
    integer :: Version
    ! Treat as plain hdf
    logical, parameter :: plainHDF = .true.

    ! Executable
    debug = LetsDebug ( 'output', 0 )
    Version = 1
    call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
    ! Get the l2gp file name from the PCF

    if ( TOOLKIT ) then
      call split_path_name(fileName, path, fileName)

      FileHandle = GetPCFromRef(fileName, pcf_start, &
        & pcf_end, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullFilename)
    else
      FullFilename = fileName
      returnStatus = 0
    end if

    if ( mls_exists(trim(FullFilename)) /= 0 ) then
      if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
        & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Its non-existence prevented writing attribute to ' &
          & // trim(FullFilename) )
      return
    elseif ( returnStatus == 0 .and. .not. checkPaths ) then
      ! Open the HDF-EOS file and write the attribute
      outputFile => GetMLSFileByName(filedatabase, FullFilename)
      if ( .not. associated(outputFile) ) then
      if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
        & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(FullFilename) )
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=FullFilename, shortName=fileName, &
          & type=l_swath, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=pcf_start, PCTop=pcf_end)
      endif
    endif
    if ( plainHDF .and. additional ) then
      call dump ( outputFile, details=1 )
      ! Due to problems extending the hdfeos implementation of file-level
      ! attributes, we resort to the plain vanilla hdf5 api
      ! We first close it as an hdfeos file and reopen it as plain hdf5
      if ( outputFile%stillOpen ) then
        call mls_closeFile( outputFile, returnStatus )
        if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
            & 'Unable to hdfeos_close ' // trim(FullFilename) )
      endif
      call h5fopen_f ( trim(outputFile%name), H5F_ACC_RDWR_F, &
        & outputFile%FileID%f_id, returnStatus )
      if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Unable to open ' // trim(FullFilename) )
      call h5gopen_f( outputFile%FileID%f_id, &
        & '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES', outputFile%FileID%grp_id, &
        & returnStatus )
      if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Unable to open group /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES for ' // &
          & trim(FullFilename) )
      outputFile%FileID%sd_id = 0 ! To ensure we read from the grp_id
      call dump ( outputFile, details=1 )
      if ( additional ) then
        call GetHDF5Attribute ( outputFile, attrname, str )
        str = trim(str) // attrValue
        if ( debug ) call output( 'Updated glob attr to ' // trim(str), advance='yes' )
      else
        str = attrValue
        if ( debug ) call output( 'Set glob attr to ' // trim(str), advance='yes' )
      endif
      call MakeHDF5Attribute( outputFile%FileID%grp_id, &
        & attrName, str, skip_if_already_there = .false. )
      call h5gclose_f( outputFile%FileID%grp_id, returnStatus )
      if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Unable to close group ' // trim(FullFilename) )
      call h5fclose_f ( outputFile%FileID%f_id, returnStatus )
      if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Unable to close  ' // trim(FullFilename) )
      return
    endif
    if ( .not. outputFile%stillOpen ) &
      & call mls_openFile( outputFile, returnStatus )
    if ( .not. additional ) then
      str = attrValue
    else
      returnStatus = HE5_EHRDGlAtt ( outputFile%FileID%f_id, attrName, str )
      str = trim(str) // attrValue
      call output( 'Updated l2gp glob attr to ' // trim(str), advance='yes' )
    endif
    returnStatus = mls_EHwrglatt( outputFile%FileID%f_id, &
     & attrName, MLS_CharType, 1, &
     &  str )
  end subroutine writeAttributeToL2GP

  ! ---------------------------------------------  OutputL2AUX  -----
  subroutine OutputL2AUX ( key, fileName, DEBUG, writeCounterMAF, &
    & filedatabase, L2AUXDatabase )
    ! Do the work of outputting specified l2aux data to a named file
    use Init_Tables_Module, only: F_Overlaps, F_Quantities
    use Intrinsic, only: L_HDF
    use L2AUXData, only: L2AUXData_T, WriteL2AUXData
    use L2GPData, only: L2GPNamelen
    use MLSCommon, only: MLSFile_T, FileNamelen, L2metaData_T
    use MLSPCF2, only: MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End
    use MLSStringLists, only: SwitchDetail
    use SDPToolkit, only: Pgs_S_Success
    use Tree, only: Decoration, Nsons, Subtree
    ! Args
    integer, intent(in)                     :: key ! tree node
    character(len=*), intent(inout)         :: fileName ! according to l2cf
    logical, intent(in)                     :: DEBUG ! Print lots?
    logical, intent(in)                     :: writeCounterMAF
    type(MLSFile_T), dimension(:), pointer  :: filedatabase
    type(L2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    ! Local variables
    integer :: DB_index
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    character (len=FileNameLen) :: FullFilename
    integer :: FileHandle
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: IN_FIELD_NO              ! Index of sons of assign vertex
    integer :: Metadata_error
    integer :: Numquantitiesperfile     ! < MAXQUANTITIESPERFILE
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=132) :: path
    character(len=L2GPNameLen), dimension(MAXQUANTITIESPERFILE) :: &
      &                           QuantityNames  ! From "quantities" field
    integer :: ReturnStatus
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: Version
    type(L2Metadata_T) :: l2metaData

    ! Executable
    Version = 1
    OUTPUTTYPESTR = 'l2aux'
    ! Get the l2aux file name from the PCF

    if ( TOOLKIT ) then
      call split_path_name(fileName, path, fileName)

      FileHandle = GetPCFromRef(fileName, MLSPCF_L2DGM_Start, &
        & MLSPCF_L2DGM_End, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullFilename)
    else
      FullFilename = fileName
      returnStatus = 0
    end if

    if ( returnStatus == 0 .and. .not. checkPaths ) then
      ! Open the HDF file and write l2aux data
      outputFile => GetMLSFileByName(filedatabase, FullFilename)
      if ( .not. associated(outputFile) ) then
        if(DEBUG) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(FullFilename) )
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=FullFilename, shortName=fileName, &
          & type=l_hdf, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2DGM_Start, PCTop=MLSPCF_L2DGM_End)
      endif
      ! Loop over the segments of the l2cf line

      numquantitiesperfile = 0
      do field_no = 2, nsons(key) ! Skip "output" name
        gson = subtree(field_no,key)
        select case ( decoration(subtree(1,gson)) )
        case ( f_quantities )
          do in_field_no = 2, nsons(gson)
            db_index = -decoration(decoration(subtree(in_field_no, gson)))
            if ( db_index >= 1 ) then
              call WriteL2AUXData ( L2AUXDatabase(db_index), outputFile, &
                & returnStatus, &
                & WriteCounterMAF = &
                &   (writeCounterMAF .and. numquantitiesperfile == 0) )
                    error = max(error, returnStatus)
              numquantitiesperfile = numquantitiesperfile+1
              if ( numquantitiesperfile > MAXQUANTITIESPERFILE ) then
                call announce_error ( son, &
                  & 'Attempt to write too many l2aux quantities to a file', &
                  & numquantitiesperfile )
                numquantitiesperfile = MAXQUANTITIESPERFILE
              end if
              call get_string &
                & ( L2AUXDatabase(db_index)%name, &
                &     QuantityNames(numquantitiesperfile) )
            else
              call MLSL2Message ( MLSMSG_Warning, ModuleName, &
                & 'Unable to write quantity to l2aux file, ' // &
                & 'perhaps no chunks processed' )
            end if
          end do ! in_field_no = 2, nsons(gson)
        case ( f_overlaps )
          ! ??? More work needed here
        end select
      end do ! field_no = 2, nsons(key)

      if ( switchDetail(switches, 'pro') > -1 ) then
        call announce_success(FullFilename, 'l2aux', &
          & numquantitiesperfile, quantityNames, hdfVersion=hdfVersion)
      end if

      if ( .not. TOOLKIT ) return

      ! Write the metadata file
      call add_metadata ( son, fileName, l2metaData, &
        & hdfVersion, l_hdf, metadata_error, &
        & numquantitiesperfile, quantityNames )
    else if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( son, &
        &  "Error finding l2aux file matching:  "//fileName, returnStatus)
    end if
  end subroutine OutputL2AUX

  ! ---------------------------------------------  CopyTextFileToHDF  -----
  subroutine CopyTextFileToHDF ( fileName, DEBUG, filedatabase, inputFile )
    ! Do the work of copying the text file to a named file
    ! If inputFile omitted, copy the l2cf
   use Intrinsic, only: L_HDF
   use MLSCommon, only: MLSFile_T, FileNameLen
   use MLSHDF5, only: SaveasHDF5ds
   use MLSPCF2, only: MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End
   use MLSStringLists, only: SwitchDetail
   use Output_M, only: Output
    ! Args
    character(len=*), intent(inout)         :: fileName ! according to l2cf
    logical, intent(in)                     :: DEBUG ! Print lots?
    type(MLSFile_T), dimension(:), pointer  :: filedatabase
    type(MLSFile_T), optional, pointer      :: inputFile
    ! Internal variables
    type(MLSFile_T), pointer                :: MLSL2CF
    integer                                 :: status
    integer :: FileHandle
    character (len=FileNameLen) :: FullFilename
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=132) :: path
    integer :: ReturnStatus
    integer :: Version
    ! Executable
    nullify ( MLSL2CF )
    if ( present(inputFile) ) then
      MLSL2CF => inputFile
    else
      MLSL2CF => GetMLSFileByType(filedatabase, content='l2cf')
    endif
    if ( .not. associated(MLSL2CF) ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Unable to write dataset--no entry in filedatabase for text file' )
      return
    endif
    if ( MLSL2CF%stillOpen ) call mls_closeFile( MLSL2CF, status )
    if ( MLSL2CF%name == '<STDIN>' ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Unable to write dataset--stdin has been used for text file' )
      return
    endif
    Version = 1
    OUTPUTTYPESTR = 'l2aux'
    ! Get the l2aux file name from the PCF

    if ( TOOLKIT ) then
      call split_path_name(fileName, path, fileName)

      FileHandle = GetPCFromRef(fileName, MLSPCF_L2DGM_Start, &
        & MLSPCF_L2DGM_End, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullFilename)
    else
      FullFilename = fileName
      returnStatus = 0
    end if

    if ( mls_exists(trim(FullFilename)) /= 0 ) then
      if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
        & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'Its non-existence prevented any copying to ' &
        & // trim(FullFilename) )
      return
    elseif ( returnStatus == 0 .and. .not. checkPaths ) then
      ! Open the HDF file and write text file
      outputFile => GetMLSFileByName(filedatabase, FullFilename)
      if ( .not. associated(outputFile) ) then
        if( DEBUG .or.  switchDetail(switches, 'pro') > -1 ) &
          & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(FullFilename) )
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=FullFilename, shortName=fileName, &
          & type=l_hdf, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=MLSPCF_L2DGM_Start, PCTop=MLSPCF_L2DGM_End)
      endif
      if ( DEBUG ) then
        call output ( 'About to copy text file to hdf', advance='yes' )
        call output ( 'input text file', advance='yes' )
        call dump( MLSL2CF, details=2 )
        call output ( 'output hdf file', advance='yes' )
        call dump( outputFile, details=2 )
      endif
      if ( .not. outputFile%stillOpen ) call mls_openFile( outputFile, status )
      call SaveAsHDF5DS ( MLSL2CF%name, outputFile%FileID%f_id, &
        & trim(MLSL2CF%content), maxLineLen=4096 )
      call mls_closeFile( outputFile, status )
    endif
  end subroutine CopyTextFileToHDF

  ! ---------------------------------------------  OutputL2GP  -----
  subroutine OutputL2GP ( key, fileName, DEBUG, &
    & output_type, pcf_start, pcf_end, &
    & filedatabase, L2GPDatabase, HGrid )
    ! Do the work of outputting specified l2gp data to a named file
    use HGridsDatabase, only: HGrid_T
    use Init_Tables_Module, only: F_Overlaps, F_Quantities
    use Intrinsic, only: L_Swath, Lit_Indices
    use L2GPData, only: L2GPData_T, L2gpnamelen, WriteL2GPData
    use MLSCommon, only: MLSFile_T, Filenamelen, L2metaData_T
    use MLSStringlists, only: Switchdetail
    use SDPToolkit, only: Pgs_S_Success
    use Tree, only: Decoration, Nsons, Subtree
    ! Args
    integer, intent(in)                     :: key ! tree node
    character(len=*), intent(inout)         :: fileName ! according to l2cf
    logical, intent(in)                     :: DEBUG ! Print lots?
    integer, intent(in)                     :: output_type
    integer, intent(in)                     :: pcf_start
    integer, intent(in)                     :: pcf_end
    type(MLSFile_T), dimension(:), pointer  :: filedatabase
    type(L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type(HGrid_T)                           :: HGrid
    ! Local variables
    integer :: DB_index
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    character (len=FileNameLen) :: FullFilename
    integer :: FileHandle
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: IN_FIELD_NO              ! Index of sons of assign vertex
    integer :: Metadata_error
    integer :: Numquantitiesperfile     ! < MAXQUANTITIESPERFILE
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: OUTPUTTYPESTR   ! 'l2gp', 'l2aux', etc.
    character (len=132) :: path
    character(len=L2GPNameLen), dimension(MAXQUANTITIESPERFILE) :: &
      &                           QuantityNames  ! From "quantities" field
    integer :: ReturnStatus
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: Version
    type(L2Metadata_T) :: l2metaData

    ! Executable
    Version = 1
    call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
    ! Get the l2gp file name from the PCF

    if ( TOOLKIT ) then
      call split_path_name(fileName, path, fileName)

      FileHandle = GetPCFromRef(fileName, pcf_start, &
        & pcf_end, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullFilename)
    else
      FullFilename = fileName
      returnStatus = 0
    end if

    if ( returnStatus == 0 .and. .not. checkPaths ) then
      ! Open the HDF-EOS file and write swath data
      outputFile => GetMLSFileByName(filedatabase, FullFilename)
      if ( .not. associated(outputFile) ) then
        if(DEBUG) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'No entry in filedatabase for ' // trim(FullFilename) )
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content=outputTypeStr, &
          & name=FullFilename, shortName=fileName, &
          & type=l_swath, access=DFACC_RDWR, HDFVersion=hdfVersion, &
          & PCBottom=pcf_start, PCTop=pcf_end)
      endif
      ! Loop over the segments of the l2cf line

      numquantitiesperfile = 0
      do field_no = 2, nsons(key) ! Skip "output" name
        gson = subtree(field_no,key)
        select case ( decoration(subtree(1,gson)) )
        case ( f_quantities )
          do in_field_no = 2, nsons(gson)
            db_index = -decoration(decoration(subtree(in_field_no ,gson)))
            if ( db_index >= 1 ) then
              call writeL2GPData ( L2GPDatabase(db_index), outputFile )
              numquantitiesperfile = numquantitiesperfile+1
              if ( numquantitiesperfile > MAXQUANTITIESPERFILE ) then
                call announce_error ( son, &
                  & 'Attempt to write too many ' // trim(outputTypeStr) // &
                  & ' quantities to a file', &
                  & numquantitiesperfile )
                numquantitiesperfile = MAXQUANTITIESPERFILE
              end if
              quantityNames(numquantitiesperfile) = L2GPDatabase(db_index)%name
            else
              call MLSL2Message ( MLSMSG_Warning, ModuleName, &
                & 'Unable to write quantity to ' // trim(outputTypeStr) // &
                & ' file, ' // &
                & 'perhaps no chunks processed' )
            end if
          end do ! in_field_no = 2, nsons(gson)
        case ( f_overlaps )
          ! ??? More work needed here
        end select
      end do ! field_no = 2, nsons(key)

      if ( switchDetail(switches, 'pro') > -1 ) then
        call announce_success(FullFilename, trim(outputTypeStr), &
          & numquantitiesperfile, quantityNames, hdfVersion=hdfVersion)
      end if

      call writeHGridComponents( filename, HGrid )
      if ( .not. TOOLKIT ) return

      ! Write the metadata file
      call add_metadata ( son, fileName, l2metaData, &
        & hdfVersion, l_swath, metadata_error, &
        & numquantitiesperfile, quantityNames )
    else if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( son, &
        &  "Error finding " // trim(outputTypeStr)  // &
        & ' file matching:  ' //fileName, returnStatus)
    end if
  end subroutine OutputL2GP

  ! ---------------------------------------------  writeHGridComponents  -----
  subroutine writeHGridComponents ( fileName, HGrid )
    use HDFEOS5, only: He5_Swclose, He5_Swopen, &
      & He5f_Acc_Rdwr, He5t_Native_Int, He5t_Native_Double
    use HGridsDatabase, only: Hgrid_T
    use MLSHDFEOS, only: He5_Ehwrglatt, Hsize, MLS_Isglatt
    ! Args
    character(len=*), intent(in) :: fileName
    type(HGrid_T)                :: HGrid
    ! Internal variables
    integer :: FileID
    integer :: H                        ! Num of profiles in HGrid
    integer :: Status
    ! Executable
    if ( .not. MLS_isGlAtt ( filename, 'HGrid_noProfs' ) ) then
      ! call dump(HGrid)
      fileID = he5_swopen(trim(fileName), HE5F_ACC_RDWR)
      h = HGrid%noProfs
      status = he5_EHwrglatt(fileID, 'HGrid_noProfs', HE5T_NATIVE_INT, HSIZE(1), &
        &  (/ h /) )
      status = he5_EHwrglatt(fileID, 'HGrid_phi', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%phi(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_geodLat', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%geodLat(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_lon', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%lon(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_time', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%time(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_solarTime', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%solarTime(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_solarZenith', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%solarZenith(1,:) )
      status = he5_EHwrglatt(fileID, 'HGrid_losAngle', HE5T_NATIVE_DOUBLE, hsize(h), &
        &  HGrid%losAngle(1,:) )
      status = he5_swclose(fileID)
    endif
  end subroutine writeHGridComponents

  ! ---------------------------------------------  returnFullFileName  -----
  subroutine returnFullFileName ( shortName, FullName, &
    & pcf_start, pcf_end )
    ! Given a possibly-abbreviated shortName, return the full name
    ! as found in the PCF
    ! (w/o toolkit panoply, simply return shortName)
    ! Args
    character(len=*), intent(in)  :: shortName
    character(len=*), intent(out) :: FullName
    integer, intent(in)           :: pcf_start
    integer, intent(in)           :: pcf_end
    ! Internal variables
    ! logical, parameter :: DEBUG = .false.
    logical :: debug
    integer :: FileHandle
    integer :: returnStatus
    integer :: Version
    ! Executable
    debug = LetsDebug ( 'output', 0 )
    if ( TOOLKIT .and. pcf_end >= pcf_start ) then
      Version = 1
      FileHandle = GetPCFromRef(shortName, pcf_start, &
        & pcf_end, &
        & TOOLKIT, returnStatus, Version, DEBUG, &
        & exactName=FullName)
      if ( returnStatus /= 0 ) FullName = shortName ! In cases omitted from PCF
    else
      FullName = shortName
    end if
  end subroutine returnFullFileName

  ! ---------------------------------------------  unsplitFiles  -----
  subroutine unsplitFiles ( key, DirectDatabase, FileDatabase, HGrids, &
    & usingSubmit, debug )
    ! Catenate any split Direct Writes of dgg/dgm files
    ! Also write various types of metadata
    ! We assume hdfVersion is 5
    use Allocate_Deallocate, only: Deallocate_Test, Allocate_Test
    use ChunkDivide_M, only: Obstructions
    use DirectWrite_M, only: DirectData_T, Dump
    use HDF5, only: H5gclose_F, H5gopen_F
    use HGridsDatabase, only: HGrids_T, Dump
    use Init_Tables_Module, only: L_L2aux, L_L2dgg
    use Init_Tables_Module, only: F_File, F_Hgrid, &
      & Field_First, Field_Last, &
      & F_Type, F_RepairGeolocations
    use Intrinsic, only: L_Swath, L_HDF, Lit_Indices
    use L2AUXData, only: CpL2AUXData, Phasenameattributes
    use L2GPData, only: AvoidUnlimitedDims, &
      & MaxSwathNamesBufSize, CpL2GPData
    use L2ParInfo, only: Parallel
    use MLSCommon, only: MLSFile_T, Filenamelen, L2metaData_T
    use MLSHDF5, only: CpHDF5Glattribute, MakeHDF5Attribute, SaveAsHDF5DS
    use MLSPCF2, only: MLSPCF_L2DGM_Start, MLSPCF_L2DGM_End, &
      & MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End
    use MLSFinds, only: FindFirst, FindNext
    use MLSStringLists, only: Array2List, SwitchDetail
    use MLSStrings, only: Trim_Safe
    use MoreTree, only: Get_Boolean
    use Output_M, only: Blanks, Output
    use PCFHdr, only: Globalattributes, DumpGlobalAttributes, &
      & H5_WriteMLSFileattr, H5_Writeglobalattr, &
      & He5_Writeglobalattr, He5_WriteMLSFileattr, &
      & WriteleapsecHDFeosattr, WriteleapsecHDF5ds, &
      & WriteutcpoleHDFeosattr, WriteutcpoleHDF5ds
    use Readapriori, only: Readaprioriattributes, Writeaprioriattributes
    use Tree, only: Decoration, Nsons, Subtree, Sub_Rosa
    ! Arguments
    integer, intent(in)                        :: key
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (MLSFile_T), dimension(:), pointer    :: FILEDATABASE
    type (HGrids_T), dimension(:), pointer     :: HGrids
    logical, intent(in)                        :: debug
    logical, intent(in)                        :: USINGSUBMIT ! Set if using the submit mechanism
    ! Local variables
    logical, parameter :: ALWAYSFILTERSWATHS = .true. ! Set Status to 'crashed'
    logical :: create2                                ! If geolocs contain Fills
    integer :: DB_index
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    integer :: FIELDVALUE               ! For get_boolean
    character (len=FileNameLen) :: FILE_BASE    ! From the FILE= field
    integer :: FILEID
    logical, dimension(field_first:field_last) :: GOT ! Fields
    integer :: GRP_ID
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: HGridIndex
    type(MLSFile_T), pointer :: inputFile
    character (len=FileNameLen) :: L2auxPhysicalFilename
    integer :: L2gpFileHandle, L2gp_Version
    character (len=FileNameLen) :: L2gpPhysicalFilename
    integer :: listSize
    logical :: madeFile
    ! logical, parameter :: NEVERDUPSWATHNAMES = .true. ! Skip cp if dup
    integer :: numswaths
    integer :: obst
    integer, dimension(:,:), pointer :: obstruction_mafs => null()
    character(len=8) :: options
    type(MLSFile_T), pointer :: outputFile
    character(len=8) :: outputTypeStr
    integer :: output_type
    character (len=MAXSWATHNAMESBUFSIZE) :: outsdList
    logical :: RepairGeolocations
    integer :: ReturnStatus
    integer :: SDFID                ! File handle
    character (len=MAXSWATHNAMESBUFSIZE) :: sdList
    type(L2Metadata_T) :: l2metaData

    ! Executable
    options = ' '
    if ( ALWAYSFILTERSWATHS ) options = '-f'
    repairGeolocations = .false.
    if ( debug ) call dump(DirectDatabase)
    got = .false.
    HGridIndex = 0
    do field_no = 2, nsons(key)       ! Skip the command name
      gson = subtree(field_no, key)   ! An assign node
      if ( nsons(gson) > 1 ) then
        fieldValue = decoration(subtree(2,gson)) ! The field's value
      else
        fieldValue = gson
      end if
      field_index = decoration(subtree(1,gson))
      got(field_index) = .true.
      select case ( field_index )   ! Field name
      case ( f_file )
        call get_string ( sub_rosa(subtree(2,gson)), file_base, strip=.true. )
      case ( f_hgrid )
        HGridIndex = decoration(fieldValue)
        if ( DEBUG ) print *, 'HGridIndex: ', HGridIndex
      case ( f_type )
        output_type = decoration(subtree(2,gson))
        call get_string ( lit_indices(output_Type), outputTypeStr, strip=.true. )
      case ( f_repairGeolocations )
        repairGeolocations = get_boolean ( fieldValue )
      case default                  ! Everything else processed later
      end select
    end do
    ! Any dgg eligible for being catenated
    DB_index = findFirst( DirectDatabase%autoType, l_l2dgg )
    if ( findNext(DirectDatabase%autoType, l_l2dgg, DB_index) > 0 ) then
      if ( TOOLKIT ) then
        l2gp_Version = 1
        l2gpFileHandle = GetPCFromRef('DGG', MLSPCF_L2DGG_Start, &
          & MLSPCF_L2DGG_End, &
          & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
          & exactName=l2gpPhysicalFilename)
      else
        file_base = DirectDatabase(DB_index)%fileName
        l2gpPhysicalFilename = unSplitName(file_base)
        returnStatus = 0
      end if
      if ( any(DirectDatabase%fileName == l2gpPhysicalFilename) ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & "Cannot unsplit dgg dw to existing file " // &
          & trim(l2gpPhysicalFilename) )
      end if
      outputFile => GetMLSFileByName(filedatabase, l2gpPhysicalFilename)
      if ( .not. associated(outputFile) ) then
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content='l2dgg', &
          & name=l2gpPhysicalFilename, shortName='DGG', &
          & type=l_swath, access=DFACC_RDWR, HDFVersion=HDFVERSION_5, &
          & PCBottom=MLSPCF_L2DGG_Start, PCTop=MLSPCF_L2DGG_End)
      endif
      madeFile = .false.
      create2 = .true.
      do DB_index = 1, size(DirectDatabase)
        if ( DirectDatabase(DB_index)%autoType /= l_l2dgg ) cycle
        numswaths = mls_InqSwath ( &
          & trim(DirectDatabase(DB_index)%fileName), sdList, listSize, &
          & hdfVersion=HDFVERSION_5 )
        if ( numSwaths < 1 ) cycle
        if ( DEBUG ) then
          call output ( 'preparing to cp split dgg', advance='yes' )
          call output ( 'from: ', advance='no' )
          call output ( trim(DirectDatabase(DB_index)%fileName) , advance='yes' )
          call output ( '   to: ', advance='no' )
          call output ( trim(l2gpPhysicalFilename) , advance='yes' )
        end if
        if ( mls_exists(trim(DirectDatabase(DB_index)%fileName)) /= 0 ) cycle
        inputFile => GetMLSFileByName(filedatabase, &
          & DirectDatabase(DB_index)%fileName)
        if ( .not. associated(inputFile) ) then
          call MLSL2Message ( MLSMSG_Error, ModuleName, &
            & 'No entry in filedatabase for ' // &
            & trim(DirectDatabase(DB_index)%fileName) )
        endif
        if ( CHECKPATHS ) cycle
        madeFile = .true.
        inputFile%access = DFACC_RDONLY
        if ( repairGeoLocations .and. got(f_HGrid) ) then
          numswaths = mls_InqSwath ( &
            & outputFile%Name, outsdList, listSize, &
            & hdfVersion=HDFVERSION_5 )
          if ( DEBUG ) then
            call outputnamedvalue( 'numswaths', numswaths )
            call outputnamedvalue( 'listSize', listSize )
            call outputnamedvalue( 'len_trim(outsdList)', len_trim(outsdList) )
          end if
          call cpL2GPData( l2metaData, inputFile, &
            & outputFile, exclude=outsdList, create2=create2, &
            & notUnlimited=avoidUnlimitedDims, &
            & andGlAttributes=COPYGLOBALATTRIBUTES, &
            & HGrid=HGrids(HGridIndex)%the_hGrid, options=trim(options) // 'cr' )
        else
          call cpL2GPData( l2metaData, inputFile, &
            & outputFile, exclude=outsdList, create2=create2, &
            & notUnlimited=avoidUnlimitedDims, &
            & andGlAttributes=COPYGLOBALATTRIBUTES, &
            & options=options )
        endif
        create2 = .false.
      end do
      ! Now write various kinds of metadata
      ! (1) Catalog metadata: only if file is just created, and TOOLKIT is available
      ! (2) File level attributes: only if file is just created
      ! (3) leapsec and utcpole contents: only if file is just created, and TOOLKIT is available
      if ( TOOLKIT .and. madeFile ) then
        ! (1) Catalog metadata
        call add_metadata ( 0, trim(l2gpPhysicalFilename), l2metaData, &
          & HDFVERSION_5, l_l2dgg, returnStatus, 1, (/'dgg'/) )
        if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'unable to addmetadata to ' // trim(l2gpPhysicalFilename) )
        if ( len_trim(l2metaData%doiIdentifier) < 1 ) then
          call MLSL2Message( MLSMSG_Warning, ModuleName, &
            & 'empty doiIdentidier for ' // trim(l2gpPhysicalFilename) )
        endif
        ! Set the file-level attribute DOI to its metadata value
        GlobalAttributes%DOI = l2metaData%doiIdentifier
      end if
      if ( madeFile ) then
         call output( 'Writing file level attributes to he5 ' // &
          & trim(l2gpPhysicalFilename), advance='yes' )
         call outputNamedValue( 'identifier_product_doi', trim(GlobalAttributes%DOI) )
         call outputNamedValue( 'ProductionLocation', trim(GlobalAttributes%productionLoc) )
         call he5_writeglobalattr( outputFile, &
           & doi=.true., skip_if_already_there=.true. )
         call DumpGlobalAttributes
       ! (2) File level attributes
        if ( WRITEFILEATTRIBUTES ) call he5_writeMLSFileAttr( outputFile )
        ! Is the next line necessary?
        if ( inputFile%stillOpen ) then
          call output( 'Closing input file before reading apriori attributes', &
            & advance='yes' )
          call MLS_CloseFile ( inputFile )
        endif
        call readAPrioriAttributes( inputFile )
        call writeAPrioriAttributes( outputFile )
      end if
      if ( TOOLKIT .and. madeFile ) then
        ! (3) leapsec and utcpole contents
        if (switchDetail( switches, 'pro') > 0 ) &
          call output ( 'About to open ' // trim(l2gpPhysicalFilename) , advance='yes' )
        call MLS_OpenFile( outputFile )
        call WriteLeapSecHDFEOSAttr ( outputFile%fileID%f_id )
        if ( .not. DGGFILEISHYBRID ) &
          & call WriteutcPoleHDFEOSAttr ( outputFile%fileID%f_id )
        call MLS_CloseFile ( outputFile )
        if ( DGGFILEISHYBRID ) then
          ! The utcpole is too large to be stored as an HDFEOS attribute
          ! Note:
          ! The dgg file up to v2.22 was already a "hybrid" file because
          ! the PCF was stored as a dataset instead of an attribute
          FileID = mls_sfstart( l2gpPhysicalFilename, DFACC_RDWR, &
              & hdfVersion=HDFVERSION_5 )
          call WriteutcPoleHDF5DS ( fileID )
          ReturnStatus = mls_sfend( fileID, hdfVersion=HDFVERSION_5 )
        endif
      end if
    end if
    ! Next we would do the same for any split dgm direct write files
    ! We must not write the phase and forward model names
    ! as held by us, the master, but instead
    ! wait until we can copy the correct values known only to the slaves
    PHASENAMEATTRIBUTES = .false.
    DB_index = findFirst( DirectDatabase%autoType, l_l2aux )
    if ( findNext(DirectDatabase%autoType, l_l2aux, DB_index) > 0 ) then
      if ( TOOLKIT ) then
        l2gp_Version = 1
        l2gpFileHandle = GetPCFromRef('DGM', MLSPCF_L2DGM_Start, &
          & MLSPCF_L2DGM_End, &
          & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
          & exactName=l2auxPhysicalFilename)
      else
        file_base = DirectDatabase(DB_index)%fileName
        l2auxPhysicalFilename = unSplitName(file_base)
        returnStatus = 0
      end if
      if ( any(DirectDatabase%fileName == l2auxPhysicalFilename) ) then
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          &  "Must not unsplit dgm dw to " // trim(l2auxPhysicalFilename) )
      end if
      outputFile => GetMLSFileByName(filedatabase, l2auxPhysicalFilename)
      if ( .not. associated(outputFile) ) then
        outputFile => AddInitializeMLSFile(filedatabase, &
          & content='l2aux', &
          & name=l2auxPhysicalFilename, shortName='L2AUX-DGM', &
          & type=l_hdf, access=DFACC_RDWR, HDFVersion=HDFVERSION_5, &
          & PCBottom=MLSPCF_L2DGM_Start, PCTop=MLSPCF_L2DGM_End)
      endif
      madeFile = .false.
      create2 = .true.
      do DB_index = 1, size(DirectDatabase)
        if ( CHECKPATHS ) cycle
        if ( DirectDatabase(DB_index)%autoType /= l_l2aux ) cycle
        if ( .not. associated(DirectDatabase(DB_index)%sdNames) ) then
          ! Someday, go back in Join and find out why these aren't saved
          ! call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          !   &  "no sd known for " // trim(DirectDatabase(DB_index)%fileName) )
          sdList = ' '
        else
          call Array2List(DirectDatabase(DB_index)%sdNames, sdList)
        endif

        if ( DEBUG ) then
          call output ( 'preparing to cp split dgm', advance='yes' )
          call output ( 'from: ', advance='no' )
          call output ( 'DB_index ', advance='no' )
          call output ( DB_index , advance='no' )
          call blanks(3)
          call output ( trim(DirectDatabase(DB_index)%fileName) , advance='yes' )
          call output ( '   to: ', advance='no' )
          call output ( trim(l2auxPhysicalFilename) , advance='yes' )
          call output ( '   sdList: ', advance='no' )
          call output ( trim(sdList) , advance='yes' )
        end if
        if ( mls_exists(trim(DirectDatabase(DB_index)%fileName)) /= 0 ) cycle
        inputFile => GetMLSFileByName(filedatabase, &
          & DirectDatabase(DB_index)%fileName)
        if ( .not. associated(inputFile) ) then
          call MLSL2Message ( MLSMSG_Error, ModuleName, &
            & 'No entry in filedatabase for ' // &
            & trim(DirectDatabase(DB_index)%fileName) )
        endif

        inputFile%access = DFACC_RDONLY
        madeFile = .true.
        if ( sdList /= ' ' ) then
          call cpL2AUXData( inputFile, &
          & outputFile, create2=create2, sdList=trim(sdList) )
        else
          ! Last-ditch effort if somehow sdNames empty or Array2List fails
          call cpL2AUXData( inputFile, &
          & outputFile, create2=create2 )
        end if
        if ( create2 ) then
          call CpHDF5GlAttribute ( DirectDatabase(DB_index)%fileName, &
            & l2auxPhysicalFilename, 'Phase Names' )
          call CpHDF5GlAttribute ( DirectDatabase(DB_index)%fileName, &
            & l2auxPhysicalFilename, 'ForwardModel Names', &
            & skip_if_already_there=.true. )
          call CpHDF5GlAttribute ( DirectDatabase(DB_index)%fileName, &
            & l2auxPhysicalFilename, 'MiscNotes' )
        endif
        create2= .false.
      end do
      ! Is metadata really needed for l2aux files? Yes.
      if ( TOOLKIT .and. madeFile ) then
        call add_metadata ( 0, trim(l2auxPhysicalFilename), l2metaData, &
          & HDFVERSION_5, l_hdf, returnStatus, 1, (/'dgm'/) )
        GlobalAttributes%DOI = l2metaData%doiIdentifier
        if ( returnStatus /= 0 ) call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'unable to addmetadata to ' // trim(l2auxPhysicalFilename) )
        call output( 'Writing file level attributes to h5 ' // &
          & trim(l2auxPhysicalFilename), advance='yes' )
        call outputNamedValue( 'identifier_product_doi', trim(GlobalAttributes%DOI) )
        call outputNamedValue( 'ProductionLocation', trim(GlobalAttributes%productionLoc) )
        call h5_writeGlobalAttr ( outputFile, &
          & skip_if_already_there=.false., doi=.true. )
        call DumpGlobalAttributes
      end if
      
      ! Now we can write any last-minute attributes or datasets to the l2aux
      if ( madeFile .and. WRITEFILEATTRIBUTES ) &
        & call h5_writeMLSFileAttr ( outputFile )
      ! E.g., parallel stuff
      if ( (parallel%master .or. switchDetail(switches, 'chu') > -1 ) &
        & .and. .not. (checkPaths .or. SKIPDIRECTWRITES) .and. &
        & madeFile .and. l2auxPhysicalFilename /= ' ' ) then
        sdfId = mls_sfstart(l2auxPhysicalFilename, DFACC_RDWR, &
            & hdfVersion=HDFVERSION_5)
        if ( TOOLKIT ) then
          call WriteLeapSecHDF5DS (sdfId)
          call WriteutcPoleHDF5DS (sdfId)
        endif
        call h5gopen_f(sdfId, '/', grp_id, returnStatus)
        if ( .not. parallel%master .and. FAKEPARALLELMASTER ) then
          parallel%numCompletedChunks = 347
          parallel%numFailedChunks = 3
          parallel%FailedChunks = '2,5,129'
          parallel%FailedMachs = 'c0-1,c0-66,c0-66'
          parallel%FailedMachs = 'msg 1\msg 2\msg 3'
        endif
        call MakeHDF5Attribute(grp_id, &
         & 'NumCompletedChunks', parallel%numCompletedChunks, .true.)
        call MakeHDF5Attribute(grp_id, &
         & 'NumFailedChunks', parallel%numFailedChunks, .true.)
        call MakeHDF5Attribute(grp_id, &
         & 'FailedChunks', trim_safe(parallel%FailedChunks), .true.)
        if ( .not. usingSubmit ) &
          call MakeHDF5Attribute(grp_id, &
           & 'FailedMachines', trim_safe(parallel%FailedMachs), .true.)
        call MakeHDF5Attribute(grp_id, &
         & 'FailedMsgs', trim_safe(parallel%FailedMsgs), .true.)
        call h5gclose_f(grp_id, returnStatus)
        ! Write 2 datasets for obstructions db
        if ( associated(obstructions) ) then
          if ( size(obstructions) > 0 ) then
            call allocate_test( obstruction_mafs, size(obstructions), 2, &
              & 'obstruction_mafs', ModuleName )
            do obst=1, size(obstructions)
              obstruction_mafs(obst, :) = obstructions(obst)%mafs
            enddo
            call SaveAsHDF5DS( sdfID, 'obstructions_range', &
              & merge( 1, 0, obstructions%range ) )
            call SaveAsHDF5DS( sdfID, 'obstruction_mafs', obstruction_mafs )
            call deallocate_test( obstruction_mafs, &
              & 'obstruction_mafs', ModuleName )
          endif
        endif
        returnStatus = mls_sfend(sdfid, hdfVersion=HDFVERSION_5)
      endif
    end if
    PHASENAMEATTRIBUTES = .true.
  end subroutine unsplitFiles
    
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module OutputAndClose

! $Log$
! Revision 2.207  2020/02/07 01:15:26  pwagner
! Restores writing metadata for Cloud file
!
! Revision 2.206  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.205  2018/05/31 22:49:27  pwagner
! Changed name of attribute to identifier_product_doi to please DAAC
!
! Revision 2.204  2018/03/22 18:15:05  pwagner
! Added command IsFileAbsent; may occur in ReadApriori, MergeGrids, and Output sections
!
! Revision 2.203  2017/10/11 23:58:40  pwagner
! Write obstructions%range as ints: 1 for T, 0 for F
!
! Revision 2.202  2017/07/10 23:09:13  pwagner
! Use MLSPCF_l2ascii_start, end for copying ascii files instead of MLSPCF_L2Clim_start
!
! Revision 2.201  2016/11/15 21:18:15  pwagner
! Use SayTime from time_m instead of local implementation
!
! Revision 2.200  2016/06/01 23:26:48  pwagner
! Be forgiving of Copy commands directed to non-existent files
!
! Revision 2.199  2016/05/18 01:37:30  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.198  2016/04/01 00:27:15  pwagner
! May now Execute a single command or a script of lines from l2cf
!
! Revision 2.197  2016/02/29 19:51:22  pwagner
! Usleep got from machine module instead of being an external
!
! Revision 2.196  2015/11/19 23:57:23  pwagner
! Now able to read from L2GP file
!
! Revision 2.195  2015/10/15 23:07:54  pwagner
! If /repairGeolocations , repair also any with chunknumber = -999
!
! Revision 2.194  2015/09/30 23:02:27  pwagner
! Catenate can repair geoLocations
!
! Revision 2.193  2015/08/03 21:45:22  pwagner
! Attempt to get dois for dgg and dgm like other prod names
!
! Revision 2.192  2015/06/19 20:54:56  pwagner
! Avoid Copying text file if it does not exist
!
! Revision 2.191  2015/03/28 02:50:35  vsnyder
! Paul added checking for Quantity type in DirectWrite
!
! Revision 2.190  2015/02/18 00:29:03  pwagner
! Debugging-type info shown only if debugging
!
! Revision 2.189  2014/11/04 01:23:50  pwagner
! Worked around hdfeos bug in implementing /additional for swath file attributes
!
! Revision 2.188  2014/10/13 18:10:01  pwagner
! Copying a swath also copies AprioriAttributes
!
! Revision 2.187  2014/10/07 00:06:47  pwagner
! May now write added file attributes
!
! Revision 2.186  2014/09/05 01:16:18  vsnyder
! Remove declarations of unused variables
!
! Revision 2.185  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.184  2014/04/14 17:39:48  pwagner
! Fixed bugs in writing DOIs for dgg, dgm files; still hard-wired, though
!
! Revision 2.183  2014/04/07 18:09:27  pwagner
! Stop Writing MLSFile_T attributes by default; they confuse users
!
! Revision 2.182  2014/04/02 23:04:34  pwagner
! Removed redundant open_ and close_MLSFile
!
! Revision 2.181  2014/03/28 00:01:10  pwagner
! repaired some bugs in writing DOI, ProdLoc attributes
!
! Revision 2.180  2014/03/26 17:46:59  pwagner
! Added ProductionLocation, identifier_product_DOI to attributes
!
! Revision 2.179  2014/02/28 01:06:44  vsnyder
! Move units checking to type checker.  Check value of hdfVersion.  Look
! for TIME field on SLEEP command, not CREATE field.
!
! Revision 2.178  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.177  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.176  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.175  2013/10/15 23:52:55  pwagner
! May copy quantity values to a file global attribute
!
! Revision 2.174  2013/10/09 23:42:13  vsnyder
! Add Evaluate_Variable
!
! Revision 2.173  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.172  2013/09/06 21:10:45  pwagner
! Move CopyQuantity into a separate subroutine
!
! Revision 2.171  2013/09/04 17:35:47  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.170  2013/08/30 02:45:45  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.169  2013/08/20 00:32:11  pwagner
! Avoid crashing when no ForwardModel Names attribute missing
!
! Revision 2.168  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.167  2012/11/08 23:20:00  pwagner
! Tries to avoid duplicating swath names during unsplit
!
! Revision 2.166  2012/08/16 17:51:43  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.165  2012/07/10 15:21:15  pwagner
! Sets Status to 'crashed' for swaths with Fills in geolocations
!
! Revision 2.164  2012/07/05 23:54:21  pwagner
! Copy command in Output may contain options field
!
! Revision 2.163  2012/05/14 18:33:34  pwagner
! Fixed bug that dooms checkPaths test flights
!
! Revision 2.162  2012/05/11 00:18:37  pwagner
! Added isSwathEmpty to set Boolean in Output; we can Skip Copy of OH when THz is off
!
! Revision 2.161  2012/05/10 00:47:02  pwagner
! Output section can have l2cf-control stuctures
!
! Revision 2.160  2012/03/14 16:57:03  pwagner
! Fixed most recent goldbrick-busting bug
!
! Revision 2.159  2012/03/12 17:28:22  pwagner
! Fixed bug preventing us from copying apriori attributes successfully
!
! Revision 2.158  2012/02/24 21:18:01  pwagner
! Correctly copy a priori attributes when unplitting files
!
! Revision 2.157  2012/02/08 23:16:30  pwagner
! Cant write utcpole, leapsec files w/o toolkit
!
! Revision 2.156  2011/11/30 21:34:23  pwagner
! Fixed bug affecting files w/o pcfids when using pcf
!
! Revision 2.155  2011/07/15 00:14:28  pwagner
! Wont crash on metadata errors now; was crashing goldbrick
!
! Revision 2.154  2011/06/29 21:51:51  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.153  2010/09/17 00:07:18  pwagner
! Can constrain writing l2pc blocks by name
!
! Revision 2.152  2010/05/20 00:32:32  pwagner
! Remove unused stuff
!
! Revision 2.151  2010/03/24 20:53:45  vsnyder
! Add DumpBlocks to Output section
!
! Revision 2.150  2010/02/25 18:18:31  pwagner
! Conforms with changed l2pc structure
!
! Revision 2.149  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.148  2010/01/08 00:11:37  pwagner
! Added ability to write MLSFile_T fields as attributes
!
! Revision 2.147  2009/11/10 00:44:03  pwagner
! Copies MiscNotes while unsplitting l2aux files, too
!
! Revision 2.146  2009/10/01 19:53:53  vsnyder
! Use strip=.true in get_string instead of unquote later
!
! Revision 2.145  2009/09/29 23:37:32  pwagner
! Changes needed by 64-bit build
!
! Revision 2.144  2009/08/26 17:18:42  pwagner
! Master copies file attributes from slaves instead of writing its own
!
! Revision 2.143  2009/07/01 20:37:30  pwagner
! Avoid adding metadata to a file more than once
!
! Revision 2.142  2009/06/26 00:17:14  pwagner
! May now copy ascii file to DGM calling input file type 'ascii' insstead of 'l2cf'
!
! Revision 2.141  2009/06/23 18:45:16  pwagner
! May copy arbitrary text file into DGM file
!
! Revision 2.140  2009/06/02 17:53:15  cvuu
! Add NRT Lat and Lon bounding to metadata
!
! Revision 2.139  2009/04/27 20:45:18  pwagner
! Fixed bug causing crashes in add_metadata when debugging
!
! Revision 2.138  2009/04/01 23:32:32  pwagner
! Writes obstructions db to l2aux file
!
! Revision 2.137  2008/12/02 23:13:15  pwagner
! mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!
! Revision 2.136  2008/09/19 23:55:50  pwagner
! May now Destroy GriddedData
!
! Revision 2.135  2008/07/09 16:38:24  pwagner
! ReadStatus optional arg eliminated
!
! Revision 2.134  2008/04/25 22:55:29  pwagner
! Exploits ability of SaveAsHDF5DS to take textfile as arg
!
! Revision 2.133  2008/03/24 17:07:02  pwagner
! Removed unused procedures
!
! Revision 2.132  2008/02/22 21:36:41  pwagner
! DGG file was hybrid by default; now it will be pure HDFEOS
!
! Revision 2.131  2007/12/07 01:51:40  pwagner
! Removed unused dummy variables, etc.
!
! Revision 2.130  2007/05/30 22:30:29  pwagner
! Tries to avoid accessing components of unassociated pointers
!
! Revision 2.129  2007/04/05 22:52:01  pwagner
! NAG, too, now able to write L2CF to L2AUX
!
! Revision 2.128  2007/03/23 00:28:39  pwagner
! Steps gingerly around case of disassociated HGrids
!
! Revision 2.127  2007/01/18 19:39:16  pwagner
! Fixed bug causing Phase Names attribute to include only 1st phase
!
! Revision 2.126  2006/10/02 23:06:56  pwagner
! Write FailedMachines attribute unless using old mlssubmit
!
! Revision 2.125  2006/08/02 19:52:55  vsnyder
! Add destroy command and destroy field for output command
!
! Revision 2.124  2006/05/09 16:40:41  pwagner
! Added writing l2cf to dgm
!
! Revision 2.123  2006/04/28 00:46:28  pwagner
! Overcome namelength (132) limitation for file= field
!
! Revision 2.122  2006/04/11 23:35:38  pwagner
! More info why unable to unsplit dgg/dgm files
!
! Revision 2.121  2006/03/15 23:47:53  pwagner
! Fixed bug causing crashes when no gapless HGrid declared
!
! Revision 2.120  2006/03/04 00:23:30  pwagner
! Will not attempt copy unless input file contains swath
!
! Revision 2.119  2006/02/21 19:13:33  pwagner
! Some tweaks to where, when to dump
!
! Revision 2.118  2005/12/21 18:46:29  pwagner
! Fixed bug that clobbered split dgm files while copying them
!
! Revision 2.117  2005/12/10 00:51:36  pwagner
! Copies ForwardModel Names attribute when unsplitting spli dgms
!
! Revision 2.116  2005/11/15 00:22:05  pwagner
! Defined Overlaps for AllChunks chunk
!
! Revision 2.115  2005/11/04 18:55:46  pwagner
! Can add metadata when copying swaths, datasets
!
! Revision 2.114  2005/10/28 23:19:01  pwagner
! Many changes related to Copy; one may fixed a bug
!
! Revision 2.113  2005/10/22 00:46:14  pwagner
! May write all-day HGrid as attributes
!
! Revision 2.112  2005/09/23 23:39:35  pwagner
! Added rename field to copy command
!
! Revision 2.111  2005/09/21 23:27:17  pwagner
! Pokes holes in all-day HGrid to match obstructions
!
! Revision 2.110  2005/09/14 00:15:32  pwagner
! Relocate unsplitFiles before l2cf commands (so may copy swaths from DGG)
!
! Revision 2.109  2005/08/19 23:35:35  pwagner
! Allow Output to repair l2gp with HGrid while copying files
!
! Revision 2.108  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.107  2005/06/16 18:43:01  pwagner
! Should not bomb if catenating split files w/o toolkit
!
! Revision 2.106  2005/06/14 20:43:49  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.105  2005/03/15 23:57:23  pwagner
! Makes error messages from about-to-die chunks dgm attributes
!
! Revision 2.104  2004/12/14 21:45:01  pwagner
! Avoid catenating non-existent split files
!
! Revision 2.103  2004/09/23 23:02:46  pwagner
! DGG gets apriori attrs; master task copies Phase Names from slave DGM
!
! Revision 2.102  2004/09/16 23:57:54  pwagner
! Now tracks machine names of failed chunks
!
! Revision 2.101  2004/09/16 00:19:05  pwagner
! Writes info on completed, failed chunks to l2aux as global attrs
!
! Revision 2.100  2004/08/26 18:52:56  pwagner
! Fixed checkPaths bug if outputting l2pc
!
! Revision 2.99  2004/08/04 23:19:58  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.98  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.97  2004/05/19 20:22:09  vsnyder
! Remove USEs for unreferenced symbols, polish some cannonballs
!
! Revision 2.96  2004/04/24 00:22:55  pwagner
! Forcibly sets ReadStatus when knitting DGG; prunes L2GP- prefix from file_base
!
! Revision 2.95  2004/03/12 00:39:37  pwagner
! Interface to cpL2GPData changed to match
!
! Revision 2.94  2004/03/03 19:23:48  pwagner
! Master task never knows actual sdList to catenate; let cpL2AUXData figure it out
!
! Revision 2.93  2004/02/19 23:57:47  pwagner
! Hopefully will not try to write metadata if no file exists
!
! Revision 2.92  2004/02/05 23:35:21  pwagner
! Fixed some bugs in catenating split dgg/dgm directwrites
!
! Revision 2.91  2004/01/27 21:37:26  pwagner
! Can catenate split l2aux files
!
! Revision 2.90  2004/01/23 01:15:00  pwagner
! Began effort to catenate split direct write files
!
! Revision 2.89  2004/01/22 00:52:19  pwagner
! Small changes regarding metadata
!
! Revision 2.88  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.87  2003/10/28 21:42:36  pwagner
! Exits with message if cant GetPCFFromRef
!
! Revision 2.86  2003/10/20 23:59:20  pwagner
! Simplified code for writing metadata
!
! Revision 2.85  2003/10/16 18:29:35  pwagner
! Should not try to write metadata twice on DirectWrite files
!
! Revision 2.84  2003/09/19 23:29:27  pwagner
! Should not be a metadata error when DirectWrite-ing CH3CN
!
! Revision 2.83  2003/09/04 22:42:47  pwagner
! Some tweaks relating to DirectWrite; may not matter
!
! Revision 2.82  2003/08/15 20:43:10  pwagner
! Downgraded to warning if directwrite output_type unkown, e.g. l_l2fwm
!
! Revision 2.81  2003/08/08 23:06:39  livesey
! Added the dontPack option on outputing l2pc files.
!
! Revision 2.80  2003/08/01 20:07:44  pwagner
! Fixed Toolkit bug; metadata distinguishes l2dgg from l2gp
!
! Revision 2.79  2003/07/23 18:29:32  cvuu
! quick and dirty fixed for CH3CN
!
! Revision 2.78  2003/07/07 23:49:11  pwagner
! Add_metadata procedure now public
!
! Revision 2.77  2003/07/07 17:31:44  livesey
! Mainly cosmetic changes
!
! Revision 2.76  2003/06/26 00:17:17  pwagner
! Writes metadata to all files in DirectDataBase
!
! Revision 2.75  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.74  2003/06/23 18:06:33  pwagner
! Should allow us to write metadata after DirectWrite
!
! Revision 2.73  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.72  2003/06/09 22:49:33  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.71  2003/05/12 02:07:06  livesey
! Bound r8->r4 conversion in direct write
!
! Revision 2.70  2003/04/03 22:59:23  pwagner
! setAlias no longer an arg to write_meta
!
! Revision 2.69  2003/03/20 19:22:56  pwagner
! Fixed bug in DirectWrite_hdf5; seems to work
!
! Revision 2.68  2003/03/11 00:21:36  pwagner
! Interfaces fit new WritePCF2Hdr flixibility
!
! Revision 2.67  2003/03/01 00:25:20  pwagner
! Disabled writing metadata to Log file (aka PH)
!
! Revision 2.66  2003/02/12 21:51:32  pwagner
! Should allow direct write with attributes
!
! Revision 2.65  2003/02/10 22:01:54  pwagner
! Passes isHDFEOS to metadata; writes globalattributes during DirectWrite
!
! Revision 2.64  2003/01/23 23:31:42  pwagner
! May directwrite to hdf5 l2aux files
!
! Revision 2.63  2002/12/11 22:21:05  pwagner
! Makes soft link to data field name from L2gpValue field in hdf5 l2gp
!
! Revision 2.62  2002/11/22 19:10:30  pwagner
! Upped MAXQUANTITIESPERFILE to 10k
!
! Revision 2.61  2002/11/13 01:10:09  pwagner
! Beginnings of attempt to write hdf5 L2AUX; incomplete
!
! Revision 2.60  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.59  2002/08/21 02:35:18  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.58  2002/08/21 01:05:06  livesey
! Changed to single precision for direct write
!
! Revision 2.57  2002/08/20 04:37:06  livesey
! Minor typo
!
! Revision 2.56  2002/08/20 04:33:13  livesey
! Added extra check in direct write
!
! Revision 2.55  2002/08/15 21:47:04  pwagner
! WriteL2AuxData now returns non-zero status if it fails
!
! Revision 2.54  2002/06/12 17:58:42  livesey
! Intermediate support for HDF5 L2PCs
!
! Revision 2.53  2002/05/22 16:30:31  livesey
! Bug fix in directWrite
!
! Revision 2.52  2002/05/22 00:49:01  livesey
! Added direct write stuff
!
! Revision 2.51  2002/05/07 20:26:15  livesey
! Added writeCounterMAF option for l2aux
!
! Revision 2.50  2002/02/22 19:19:48  pwagner
! Fixed bug in metaName use
!
! Revision 2.49  2002/02/22 01:16:17  pwagner
! Uses new metaName field for mcf file hint
!
! Revision 2.48  2002/01/29 23:49:38  pwagner
! Separate DEFAULT_HDFVERSION_(READ)(WRITE)
!
! Revision 2.47  2002/01/26 00:10:45  pwagner
! Correctly sets hdfVersion; changed proclaim to announce_success
!
! Revision 2.46  2002/01/23 21:52:15  pwagner
! Accepts and uses hdfVersion optional field
!
! Revision 2.45  2002/01/18 23:07:48  pwagner
! Uses MLSFiles instead of HDFEOS
!
! Revision 2.44  2002/01/18 00:24:34  livesey
! Added packed option to outputing l2pc files
!
! Revision 2.43  2001/11/20 00:48:54  livesey
! Alleviated one bug in zero chunks case, but there's another one to
! fix later.  We need to decide how to handle this one.
!
! Revision 2.42  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.41  2001/11/01 21:05:10  pwagner
! Satisfies new WriteL2AuxData interface
!
! Revision 2.40  2001/10/12 23:12:23  pwagner
! Checks that number of quantities written to a file not grow too large
!
! Revision 2.39  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.38  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.37  2001/06/04 23:57:40  pwagner
! Splits path from l2cf-defined file name before getPCfromRef
!
! Revision 2.36  2001/05/30 23:03:13  pwagner
! Moved PCFL2CFSAMECASE to MLSL2Options
!
! Revision 2.35  2001/05/17 22:33:28  pwagner
! Prints info if pro switch set
!
! Revision 2.34  2001/05/04 23:22:13  pwagner
! Detachable from Toolkit; created metafiles conditionally
!
! Revision 2.33  2001/05/04 23:19:55  pwagner
! Detachable from Toolkit; created metafiles conditionally
!
! Revision 2.32  2001/05/03 20:32:33  vsnyder
! Add a nullify and some cosmetic changes
!
! Revision 2.31  2001/05/01 23:57:23  pwagner
! Added l2dgg output type
!
! Revision 2.30  2001/04/28 01:30:52  livesey
! Stuff formerly outputting L2PCs is now outputting matrices.
!
! Revision 2.29  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.28  2001/04/26 15:59:13  livesey
! Fixed arguments to writeOneL2PC
!
! Revision 2.27  2001/04/25 21:51:28  livesey
! Minor changes, add canWriteL2PC flag
!
! Revision 2.26  2001/04/25 20:34:04  livesey
! Now writes l2pc files
!
! Revision 2.25  2001/04/20 23:51:24  vsnyder
! Improve an error message.  Add an option to consider the penalty in
! Announce_Error.  Numerous cosmetic changes.
!
! Revision 2.24  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.23  2001/04/16 23:51:08  pwagner
! Gets penalty from MLSL2Options
!
! Revision 2.22  2001/04/13 23:48:37  pwagner
! Removed bogus increment of l2gp_mcf
!
! Revision 2.21  2001/04/13 00:26:23  pwagner
! Whether files named in PCF agree in case with l2cf controlled
!
! Revision 2.19  2001/04/10 23:01:57  pwagner
! Now works better; tacks if no metadata
!
! Revision 2.18  2001/04/09 23:44:34  pwagner
! Fewer mistakes, more debug-type output
!
! Revision 2.17  2001/04/07 00:13:44  pwagner
! Extra error checks
!
! Revision 2.16  2001/04/04 23:44:52  pwagner
! Now deallocates anText of PCF file at last
!
! Revision 2.15  2001/04/03 23:51:28  pwagner
! Many changes; some may be right
!
! Revision 2.14  2001/04/02 23:43:46  pwagner
! Now makes metadata calls; it compiles, but does it bomb?
!
! Revision 2.13  2001/03/28 00:23:20  pwagner
! Made tiny changes to use announce_error
!
! Revision 2.12  2001/03/20 18:35:02  pwagner
! Using GetPCFromRef to get file handles
!
! Revision 2.11  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.10  2001/03/06 22:40:24  livesey
! Working l2aux
!
! Revision 2.9  2001/02/23 18:15:48  livesey
! Added trace calls.
!
! Revision 2.8  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.7  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.6  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.5  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.4  2000/11/16 02:25:13  vsnyder
! Implement timing.
!
! Revision 2.3  2000/10/05 16:43:00  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:43:47  ahanzel
! Removed old log entries.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

