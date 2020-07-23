! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MergeGridsModule

  ! This module contains code for merging operational gridded data with apriori
  ! information.
  ! Secondary operations may be performed directly on the gridded data--
  ! e.g., calculating wmo tropopause pressures from eta-level temperatures
  use Allocate_Deallocate, only: Allocate_Test, Byte_Size, Bytes, &
    & Deallocate_Test, Test_Allocate, Test_Deallocate
  use HighOutput, only: Dump, OutputNamedValue
  use, Intrinsic :: ISO_C_Binding, only: C_Intptr_T, C_Loc
  use MLSL2Options, only: MLSL2Message, L2cfNode
  use MLSL2Timings, only: AddPhaseToPhaseNames
  use MLSMessageModule, only: MLSMsg_Error, MLSMsg_Warning, DumpConfig
  use Ncep_Dao, only: ReadGriddedData
  use Output_M, only: Blanks, Output, OutputOptions, StampOptions
  implicit none
  private

  public :: Concatenate, ConvertEtaToP, DeleteGriddedData, &
    & MergeGrids, MergeOneGrid, WMOTropFromGrid

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ===================================  Public procedures  =====

  ! -------------------------------------------------  MergeGrids  -----

  subroutine MergeGrids ( root, L2GPDatabase, L2AUXDatabase, &
    & GriddedDatabase, FileDatabase )

    use DumpCommand_M, only: BooleanFromEmptyGrid, BooleanFromEmptySwath, &
      & BooleanFromFormula, DumpCommand, ExecuteCommand, &
      & MLSCase, MLSEndSelect, MLSSelect, MLSSelecting, Skip
    use GriddedData, only: GriddedData_T, &
      & AddGriddedDataToDatabase, DestroyGriddedData, Dump, DumpDBFootPrint
    use Init_Tables_Module, only: F_Grid, S_Boolean, &
      & S_Case, S_ChangeSettings, S_Concatenate, S_ConcatenateGrids, &
      & S_ConvertEtaToP, S_Delete, S_Diff, S_Dump, &
      & S_EndSelect, S_Execute, S_Gridded, S_IsFileAbsent,  &
      & S_IsGridEmpty, S_Merge, S_MergeGrids, S_Reevaluate, S_Select, S_Skip, &
      & S_Wmotrop, S_Wmotropfromgrids
    use L2AUXData, only: L2auxData_T
    use L2GPData, only: L2gpData_T
    use MLSCommon, only: MLSFile_T
    use MLSStringlists, only: SwitchDetail
    use MoreTree, only: Get_Label_And_Spec, Get_Spec_Id
    use Next_Tree_Node_M, only: Next_Tree_Node, Next_Tree_Node_State
    use ReadAPriori, only: ProcessOneAprioriFile
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Decorate, Decoration
    use Toggles, only: Gen, Switches, Toggle

    integer, intent(in) :: ROOT         ! Tree root
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 
    type (MLSFile_T), dimension(:), pointer ::     FileDatabase

    ! Local variables
    type (GriddedData_T), pointer :: Grid
    integer :: GSON
    integer :: J                        ! Loop counter
    integer :: KEY                      ! Another node
    integer :: LastAprioriPCF = 1
    integer :: LastClimPCF    = 1
    integer :: LastDAOPCF     = 1
    integer :: LastGEOS5PCF   = 1
    integer :: LastHeightPCF  = 1
    integer :: LastNCEPPCF    = 1
    integer :: Me = -1                  ! String index for trace
    integer :: NAME                     ! Index into string table
    integer :: SON                      ! Tree node
    type(next_tree_node_state) :: State ! of tree traverser
    integer :: Value
    logical :: verbose
    logical :: verboser

    ! excutable code
    call trace_begin ( me, "MergeGrids", root, cond=toggle(gen) )

    verbose = ( switchDetail(switches, 'grid' ) > -1 )
    verboser = ( switchDetail(switches, 'grid' ) > 0 )
    
    do 
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit
      call get_label_and_spec ( son, name, key )
      L2CFNODE = key
      if ( MLSSelecting .and. &
        & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) cycle

      select case ( get_spec_id(key) )
      case ( s_Boolean )
        call decorate ( key,  BooleanFromFormula ( name, key ) )
      case ( s_select ) ! ============ Start of select .. case =========
        ! We'll start seeking a matching case
        call MLSSelect (key)
      case ( s_case ) ! ================ seeking matching case =========
        ! We'll continue seeking a match unless the case is TRUE
        call MLSCase (key)
      case ( s_endSelect ) ! =========== End of select .. case =========
        ! We'done with seeking a match
        call MLSEndSelect (key)
      case ( s_changeSettings ) ! ===============================  changeSettings ==
        ! Change settings for this phase
        call addPhaseToPhaseNames ( 0, key )
        ! How verbose must we be to Dump the new settings?
        if ( .not. verboser ) cycle
        call Dump( OutputOptions )
        call Dump( StampOptions )
        call DumpConfig
      case ( s_concatenate )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & Concatenate ( key, griddedDataBase ) ) )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_ConvertEtaToP )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & ConvertEtaToP ( key, griddedDataBase ) ) )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_delete )
        call DeleteGriddedData ( key, griddedDatabase )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_diff, s_dump )
        call dumpCommand ( key, griddedDataBase=griddedDataBase )
      case ( s_execute ) ! ======================== ExecuteCommand ==========
        call ExecuteCommand ( key )
      case ( s_Gridded )
        call processOneAprioriFile ( key, L2GPDatabase, L2auxDatabase, &
          & GriddedDatabase, fileDataBase, &
          & LastAprioriPCF , &
          & LastClimPCF    , &
          & LastDAOPCF     , &
          & LastGEOS5PCF   , &
          & LastHeightPCF  , &
          & LastNCEPPCF     &
            )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_isFileAbsent )
        call decorate ( key, BooleanFromEmptySwath ( key ) )
      case ( s_isGridEmpty )
        call decorate ( key, &
          & BooleanFromEmptyGrid ( key, griddedDataBase ) )
      case ( s_merge )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & MergeOneGrid ( key, griddedDataBase ) ) )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_concatenateGrids, s_mergeGrids )
        ! We must get "grid" field from command
        do j = 2, nsons(key)
          gson = subtree(j, key)
          select case ( decoration(subtree(1, gson) ) )
          case ( f_grid )
            value = decoration ( decoration ( subtree(2, gson) ) )
          case default
          end select
        enddo
        grid => griddedDataBase(value)
        call DestroyGriddedData( grid )
        if ( get_spec_id(key) == s_mergeGrids ) then
          grid = MergeOneGrid ( key, griddedDataBase )
        else
          grid = Concatenate ( key, griddedDataBase )
        endif
        if ( verbose ) then
          call output( 'The GriddedDatabase, ' )
          call outputNamedValue( 'size(db)', size(griddedDataBase) )
          call outputNamedValue( 'our index', value )
          call outputNamedValue( 'is it empty?', grid%empty )
        endif
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_reevaluate )
        call decorate ( key,  BooleanFromFormula ( 0, key ) )
      case ( s_skip ) ! ================================ Skip ==========
        ! We'll skip the rest of the section if the Boolean cond'n is TRUE
        if ( Skip(key) ) exit
      case ( s_wmoTrop )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & wmoTropFromGrid ( key, griddedDataBase ) ) )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case ( s_wmoTropFromGrids )
        ! We must get "grid" field from command
        do j = 2, nsons(key)
          gson = subtree(j, key)
          select case ( decoration(subtree(1, gson) ) )
          case ( f_grid )
            value = decoration ( decoration ( subtree(2, gson) ) )
          case default
          end select
        enddo
        grid => griddedDataBase(value)
        call DestroyGriddedData( grid )
        grid = wmoTropFromGrid ( key, griddedDataBase )
        if ( DumpDBFootPrint ) call Dump ( griddedDataBase, Details=-4 )
      case default
        ! Shouldn't get here if parser worked?
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'Unrecognized command in MergeGrids section' )
      end select
    end do
    call trace_end ( "MergeGrids", cond=toggle(gen) )    

  end subroutine MergeGrids

  ! ----------------------------------------------  ConvertEtaToP  -----
  type (griddedData_T) function ConvertEtaToP ( root, griddedDataBase ) &
    & result ( newGrid )
    use GriddedData, only: GriddedData_T, Dump, NullifyGriddedData, &
      & ConvertFromEtaLevelGrids
    use Init_Tables_Module, only: F_A, F_B, F_Grid, F_VGrid, F_LogBasis
    use MoreTree, only: Get_Boolean
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Decoration
    use VGridsDatabase, only: VGrid_T, VGrids
    ! use VGridsDatabase, only: VGrid_T, VGrids, ConvertVGrid
    
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    ! This routine parses the l2cf instructions that
    ! convert two gridded data on eta surfaces, one of them pressures,
    ! to pressure surfaces

    ! Local variables
    logical :: ByLog                  ! Do logarithmic interpolation?
    integer :: FIELD                  ! Another tree node
    integer :: FIELD_INDEX            ! Type of tree node
    integer :: I                      ! Loop counter
    integer :: Me = -1                ! String index for trace
    integer :: returnStatus
    integer :: SON                    ! Tree node
    integer :: VALUE                  ! Tree node

    type (griddedData_T), pointer :: A ! Temperatures on eta surfaces
    type (griddedData_T), pointer :: B ! Pressures on eta surfaces
    type (griddedData_T), pointer :: V ! Grid with proper pressure surfaces
    type (VGrid_T), pointer       :: VGrid ! Desired pressure surfaces

    logical, parameter :: DEEBUG = .false.

    ! Executable code
    call trace_begin ( me, "ConvertEtaToP", root, cond=toggle(gen) )
    nullify( a, b, v, vGrid )
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    ByLog = .false.

    if ( DEEBUG ) then
      call dump(griddedDataBase)
    end if
    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field = subtree(1,son)
      if ( nsons(son) > 1 ) value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_logBasis ) 
        ByLog = get_boolean(son)
      case ( f_a ) 
        a => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we defer reading a?
        if ( a%empty .and. a%deferReading ) then
          call readGriddedData ( a%sourceFileName, son, a%description, &
            & a%verticalCoordinate, a, returnStatus, &
            & a%dimList, TRIM(a%fieldNames), a%missingValue )
        endif
      case ( f_b )
        b => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we defer reading b?
        if ( b%empty .and. b%deferReading ) then
          call readGriddedData ( b%sourceFileName, son, b%description, &
            & b%verticalCoordinate, b, returnStatus, &
            & b%dimList, TRIM(b%fieldNames), b%missingValue )
        endif
      case ( f_grid )
        v => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we fail reading v?
        if ( .not. associated(v) ) &
          & call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'The v (climatology?) grid for the conversion is not associated' )
        ! Did we defer reading v?
        if ( v%empty .and. v%deferReading ) then
          call readGriddedData ( v%sourceFileName, son, v%description, &
            & v%verticalCoordinate, v, returnStatus, &
            & v%dimList, TRIM(v%fieldNames), v%missingValue )
        endif
        ! Did we succeed in reading v (at last?)
        if ( v%empty )  call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'The v (climatology?) grid for the conversion is empty' )
      case ( f_vgrid )
        vGrid => VGrids ( decoration ( decoration ( value ) ) )
        ! Did we fail reading v?
        if ( .not. associated(vGrid) ) &
          & call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'The vgrid for the conversion is not associated' )
      end select
    end do
    if (  .not. associated(vGrid) .and.  .not. associated(v) ) then
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'Either v or vgrid for the conversion must be specified' )
    elseif ( DEEBUG ) then
      call output( 'Have either v or vgrid', advance='yes' )
    endif
    if ( DEEBUG ) call output( 'Have T, P grids', advance='yes' )
    newGrid%empty = .true.
    if ( DEEBUG ) &
      & call outputNamedValue( 'size(griddedDataBase)', size(griddedDataBase) )
    if ( size(griddedDataBase) < 2 ) go to 9
    if ( DEEBUG ) then
      call output( 'About to check on a, b', advance='yes' )
      call outputNamedValue( 'associated(a)', associated(a) )
      call outputNamedValue( 'associated(b)', associated(b) )
      call outputNamedValue( 'a%empty', a%empty )
      call outputNamedValue( 'b%empty', b%empty )
    endif
    if ( a%empty .or. b%empty ) go to 9
    newGrid%empty = .false.
    if ( DEEBUG ) then
      call output( 'a grid', advance='yes' )
      call dump( a, details=0 )
      call output( 'b grid', advance='yes' )
      call dump( b, details=0 )
      call output( 'v grid', advance='yes' )
      call dump( v, details=0 )
      call output( 'about to convert from eta level grids', advance='yes' )
    endif
    call ConvertFromEtaLevelGrids ( a, b, V, newGrid, VGrid, ByLog )
    if ( DEEBUG ) call output( 'done converting from eta level grids', advance='yes' )
    newGrid%sourceFileName      = a%sourceFileName
    newGrid%quantityName        = a%quantityName
    newGrid%description         = 'Converted Grids'
    newGrid%units               = a%units
    newGrid%verticalCoordinate  = v%verticalCoordinate
    newGrid%missingValue        = a%missingValue
9   call trace_end ( "ConvertEtaToP", cond=toggle(gen) )

  end function ConvertEtaToP

  ! ------------------------------------------------  Concatenate  -----
  function Concatenate ( root, griddedDataBase ) result ( newGrid )
    use GriddedData, only: GriddedData_T, Dump, &
      & ConcatenateGriddedData, CopyGrid, DestroyGriddedData, NullifyGriddedData
    use Init_Tables_Module, only: F_A, F_B, F_AllowEmptyGrids, F_DeleteGrids, &
      & F_SourceGrid
    use MLSStringlists, only: SwitchDetail
    use MoreTree, only: Get_Boolean, Get_Field_Id
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Decoration
    
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    type (griddedData_T), target :: newGrid
    ! This routine parses the l2cf instructions that request
    ! a grid concatenation, and then performs the concatenation

    ! Local variables
    type (griddedData_T), pointer :: A
    logical :: ATLEASTONEGRID
    type (griddedData_T), pointer :: B
    integer :: db_index
    logical, parameter            :: DEEBUG = .false.
    logical :: deleteGrids
    integer :: FIELD                  ! Another tree node
    integer :: FIELD_INDEX            ! Type of tree node
    integer :: GRIDS_NODE
    integer :: I                      ! Loop counter
    logical :: AllowEmptyGrids !  = .true. ! .false.
    type (griddedData_T), target :: Intermediate
    integer :: Me = -1                ! String index for trace
    integer :: returnStatus
    integer :: SON                    ! Tree node
    integer :: VALUE                  ! Tree node
    logical :: verbose
    logical :: WEARETHEFIRST

    ! Executable code
    call trace_begin ( me, "Concatenate", root, cond=toggle(gen) )
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    call nullifyGriddedData ( Intermediate ) ! for Sun's still useless compiler
    deleteGrids = .false.
    verbose = ( switchDetail(switches, 'grid' ) > -1 )
    AllowEmptyGrids = .false. ! Defaults to pre-v5.1 behavior

    ! Get the information from the l2cf
    grids_node = 0
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field_Index = get_field_id(son)
      if ( nsons(son) > 1 ) then
        field = subtree(1,son)
        value = subtree(2,son)
      else
        field = son ! Won't actually be used
        ! fieldValue = son
      end if
      select case ( field_index )
      case ( f_a ) 
        a => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we defer reading a?
        if ( a%empty .and. a%deferReading ) then
          call readGriddedData ( a%sourceFileName, son, a%description, &
            & a%verticalCoordinate, a, returnStatus, &
            & a%dimList, TRIM(a%fieldNames), a%missingValue )
        endif
      case ( f_b )
        b => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we defer reading b?
        if ( b%empty .and. b%deferReading ) then
          call readGriddedData ( b%sourceFileName, son, b%description, &
            & b%verticalCoordinate, b, returnStatus, &
            & b%dimList, TRIM(b%fieldNames), b%missingValue )
        endif
      case ( f_deleteGrids )
        deleteGrids = get_boolean(son)
      case ( f_sourceGrid )
        grids_node = son
      case ( F_AllowEmptyGrids )
        AllowEmptyGrids = get_boolean(son)
      end select
    end do

    ! Do the concatenation unless:
    ! AllowEmptyGrids is FALSE
    !    one or other is empty
    ! AllowEmptyGrids is TRUE
    !    all are empty
    if ( grids_node > 0 ) then
      ! Method:
      ! At any step let the result of all prior steps be held in "Intermediate"
      ! Then at each step concatenate the next gridded data with Intermediate
      ! When done, copy Intermediate into result
      ! 1st--check if any are empty; bail out if any or all are
      newGrid%empty = .true.
      atleastonegrid = .false.
      do i=2, nsons(grids_node)
        db_index = decoration(decoration(subtree(i, grids_node )))
        b => griddedDataBase ( db_index )
        ! Did we defer reading b?
        if ( b%empty .and. b%deferReading ) then
          call readGriddedData ( b%sourceFileName, grids_node, b%description, &
            & b%verticalCoordinate, b, returnStatus, &
            & b%dimList, TRIM(b%fieldNames), b%missingValue )
        endif
        if ( b%empty .and. .not. AllowEmptyGrids ) then
          call trace_end ( "Concatenate", cond=toggle(gen) )
          return
        endif
        atleastonegrid = .true.
      enddo
      if ( .not. atleastonegrid ) return
      newGrid%empty = .false.
      wearethefirst = .true.
      do i=2, nsons(grids_node)
        db_index = decoration(decoration(subtree(i, grids_node )))
        b => griddedDataBase ( db_index )
        if ( b%empty ) then
          if ( verbose ) print *, 'empty grid at i ', i, ' db_index: ', db_index
          cycle
        endif
        if ( DEEBUG ) then
          print *, ' '
          print *, 'db_index: ', db_index
          call dump( b, details=-1 )
          call outputnamedValue( 'b%equivalentLatitude', b%equivalentLatitude )
        endif
        if ( wearethefirst ) then
          call CopyGrid ( Intermediate, b )
          wearethefirst = .false.
        else
          call ConcatenateGriddedData ( A, B, Intermediate )
          if ( DEEBUG ) then
            print *, ' '
            print *, 'Result of intermediate concatenate'
            call dump( Intermediate, details=-1 )
          endif
        endif
        if ( deleteGrids ) call DestroyGriddedData ( B )
        call CopyGrid ( newGrid, Intermediate )
        a => newGrid
      enddo
      ! call CopyGrid ( newGrid, Intermediate )
      call DestroyGriddedData ( Intermediate )
      if ( DEEBUG ) call dump( newGrid, details=-1 )
    elseif ( .not. a%empty .and. .not. b%empty ) then
      call ConcatenateGriddedData ( A, B, newGrid )
    else if ( a%empty ) then
      ! Copy B into the result, of course, that may be empty too
      ! in which case the result is empty, no problem!
      call CopyGrid ( newGrid, b )
    else
      ! Otherwise a must be full, b empty
      call CopyGrid ( newGrid, a )
    end if
    newGrid%sourceFileName      = a%sourceFileName
    newGrid%quantityName        = a%quantityName
    newGrid%description         = 'Concatenated grids'
    newGrid%heightsUnits        = a%heightsUnits
    newGrid%units               = a%units
    newGrid%verticalCoordinate  = a%verticalCoordinate
    newGrid%equivalentLatitude  = a%equivalentLatitude
    newGrid%missingValue        = a%missingValue
    call outputnamedValue( 'a%equivalentLatitude', a%equivalentLatitude )
    call outputnamedValue( 'newGrid%equivalentLatitude', newGrid%equivalentLatitude )

    call trace_end ( "Concatenate", cond=toggle(gen) )

  end function Concatenate

  ! ------------------------------------------  DeleteGriddedData  -----
  subroutine DeleteGriddedData ( root, griddedDataBase )
    use Tree, only: Nsons, Subtree, Decoration
    use GriddedData, only: DestroygriddedData, GriddedData_T
    use Init_Tables_Module, only: F_Grid
    ! This routine deletes the grid indicated by the l2cf
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: GRIDDEDDATABASE ! Database
    ! Local variables
    type (griddedData_T), pointer :: GRID
    integer :: FIELD                    ! Tree node
    integer :: FIELD_INDEX              ! Tree node type
    integer :: I                        ! Counter
    integer :: SON                      ! Tree node
    integer :: VALUE                    ! Tree node

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    ! In this case there is only one argument anyway
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_grid ) 
        grid => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do
    call DestroyGriddedData ( grid )
  end subroutine DeleteGriddedData

  ! ----------------------------------------- MergeOneGrid
  type (griddedData_T) function MergeOneGrid ( root, griddedDataBase ) &
    & result ( newGrid )
    use Dump_0, only: Dump
    use Expr_M, only: Expr
    use GriddedData, only: GriddedData_T, Rgr, V_Is_Pressure, &
      & CopyGrid, Dump, NullifyGriddedData, &
      & SetupNewGriddedData, SliceGriddedData, WrapGriddedData
    use Init_Tables_Module, only: F_Climatology, F_Height, &
      & F_Operational, F_Scale
    use MLSKinds, only: R8
    use HyperSlabs, only: EssentiallyEqual
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Decoration

    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! This routine creates a new grid being a merge of two others.
    ! The operational grid forms the bottom of the dataset
    ! The climatology grid the top.  The result has the horizontal
    ! coordinates to operational and the vertical coordiantes of climatology

    ! Note! This routine is far from efficient. Not least because it
    ! uses l3atascii_interp_field which isn't terribly efficient either.
    ! But hey, this isn't a key part of the software when it comes
    ! to a desire for speed (or at least it shouldn't be)

    ! I'll need to think about missing data at some point.

    ! Local parameters
    real (r8), parameter :: SCALEHEIGHT = 16.0e3_r8 ! Approximate scale height / m

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: DAY                      ! Loop counter
    logical, parameter :: DEEBUG = .false.
    integer :: FIELD                    ! Another tree node
    integer :: FIELD_INDEX              ! Type of tree node
    integer :: I                        ! Loop inductor
    integer :: LAT                      ! Loop counter
    integer :: LON                      ! Loop counter
    integer :: LST                      ! Loop counter
    integer :: Me = -1                  ! String index for trace
    integer :: S                        ! Size in bytes of a deallocated field
    integer :: SON                      ! Tree node
    integer :: STATUS                   ! Flag from allocate
    integer :: SURF                     ! Loop counter
    integer :: SZA                      ! Loop counter
    integer :: EXPRUNITS(2)             ! Units for expr
    integer :: VALUE                    ! Value of tree node

    real (r8) :: CLIWEIGHT              ! Climatological 'weight'
    real (r8) :: HEIGHT                 ! Transition height
    real (r8) :: OPWEIGHT               ! Operational 'weight'
    real (r8) :: SCALE                  ! Transition scale
    real (r8) :: TOTALWEIGHT            ! Total weight
    real (r8) :: EXPRVALUES(2)          ! Value of expr
    real (r8) :: ZTRANS                 ! Transition 'height'
    real (r8) :: Z                      ! One 'height'
    real (r8) :: Z1, Z2                 ! Range of transition region

    real (r8) :: CLIVAL                ! One interpolated value
    real (r8) :: OPVAL                 ! One interpolated value
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: CLIMAPPED
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: OPERMAPPED
    real (r8), dimension(:), pointer :: MEANDATES ! Mean dates for new grid

    type (griddedData_T), pointer :: CLIMATOLOGY => null()
    type (griddedData_T), pointer :: OPERATIONAL => null()
    integer :: numMissingClimatology
    integer :: numMissingOperational

    ! Executable code

    call trace_begin ( me, "MergeOneGrid", root, cond=toggle(gen) )
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_operational ) 
        operational => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_climatology )
        climatology => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_height )
        call expr ( value, exprUnits, exprValues )
        height = exprValues(1)
      case ( f_scale )
        call expr ( value, exprUnits, exprValues )
        scale = exprValues(1)
      end select
    end do

    ! Think about cases where one or other grid is empty
    if ( climatology%empty ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'The climatology grid for the merge is empty' )

    if ( operational%empty ) then
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
      & 'The meteorology grid for the merge is empty' )
      ! If no operational data, then just use climatology
      call CopyGrid ( newGrid, climatology )
      call finishUp ( done = .true. )
      return
    end if
    call finishUp ( done = .false. )
    ! Do some final sanity checks
    if ( operational%verticalCoordinate /= v_is_pressure ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Operational grid not on pressure surfaces' )
    if ( climatology%verticalCoordinate /= v_is_pressure ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Climatology grid not on pressure surfaces' )
    !  if ( climatology%units /= operational%units ) &
    !    & call MLSL2Message ( MLSMSG_Error, ModuleName, &
    !    & 'The climatology and operational data describe different physical quantities' )
    if ( climatology%equivalentLatitude .neqv. &
      & operational%equivalentLatitude ) then
      call output( 'Climatology', advance='yes' )
      call dump( climatology )
      call output( 'Meteorology', advance='yes' )
      call dump( operational )
      call outputNamedValue( 'Climatology%equivalentLatitude', Climatology%equivalentLatitude )
      call outputNamedValue( 'operational%equivalentLatitude', operational%equivalentLatitude )
      call outputNamedValue( '=?', &
        & Climatology%equivalentLatitude .eqv. operational%equivalentLatitude )
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Climatology, operational data are mixed latitude/equivalent latitude.' )
    endif
    ! OK, now we're ready to go.
    ! First we're going to 'wrap' the climatology to be sure that we can
    ! interpolate it in longitude.  The chances are that it has no longitudinal
    ! variation anyway, so this won't actually do anything
    call WrapGriddedData ( climatology )
    ! Create the result.  It has the same vertical coordinates as climatology
    ! But the same horizontal coordinates as operational.
    call SetupNewGriddedData ( newGrid, source=operational, &
      & noHeights=climatology%noHeights )
    ! Setup the rest of the quantity
    newGrid%sourceFileName = 'Result of merge'
    newGrid%quantityName   = 'Result of merge'
    newGrid%description    = 'Result of merge'
    newGrid%units          = climatology%units
    newGrid%verticalCoordinate = v_is_pressure
    newGrid%equivalentLatitude = climatology%equivalentLatitude
    newGrid%heights = climatology%heights
    newGrid%lats = operational%lats
    newGrid%lons = operational%lons
    newGrid%lsts = operational%lsts
    newGrid%szas = operational%szas
    newGrid%dateStarts = operational%dateStarts
    newGrid%dateEnds = operational%dateEnds

    ! Get the 'mean' dates for the result
    nullify ( meanDates )
    call Allocate_test ( meanDates, newGrid%noDates, 'meanDates', ModuleName )
    meanDates = ( newGrid%dateStarts + newGrid%dateEnds ) / 2.0

    ! call outputNamedValue( 'Bytes before allocating 2 temp arrays', NoBytesAllocated )
    ! Now create two fields the same shape as the new field that contain
    ! the operational and climatological data interpolated to our new locations.
    ! Whoa! Now you need a total of 3x the size of the merged data set!
    ! No wonder we are running low on memory!
    allocate ( operMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, 1 ), stat=status )
      !& newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    ! if ( status /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
    !  & MLSMSG_Allocate//'operMapped' )
    addr = 0
    if ( status == 0 ) then
      if ( size(operMapped) > 0 ) addr = transfer(c_loc(operMapped(1,1,1,1,1,1)), addr)
    end if
    call test_allocate ( status, moduleName, 'operMapped', (/1,1,1,1,1,1/), &
      & (/ newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, 1 /), bytes(operMapped), address=addr )
    allocate ( cliMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, 1 ), stat=status )
      !& newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    ! if ( status /= 0 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
    !   & MLSMSG_Allocate//'operMapped' )
    if ( status == 0 ) then
      if ( size(cliMapped) > 0 ) addr = transfer(c_loc(cliMapped(1,1,1,1,1,1)), addr)
    end if
    call test_allocate ( status, moduleName, 'cliMapped', (/1,1,1,1,1,1/), &
      & (/ newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, 1 /), bytes(cliMapped), address=addr )
    ! call outputNamedValue( 'Bytes after allocating 2 temp arrays', NoBytesAllocated )

    if ( DEEBUG ) then
      call dump ( operMapped(:,1:10,1,:,1,1), &
          & '    operational field values (1st longitude) =' , &
          & FillValue=newGrid%MissingValue )
      call dump ( cliMapped(:,1:10,1,:,1,1), &
          & '    climatology field values (1st longitude) =' , &
          & FillValue=newGrid%MissingValue )
    endif
    zTrans = scaleHeight * ( 3.0 - log10 ( height ) )
    z1 = zTrans - scale/2.0
    z2 = zTrans + scale/2.0

    ! Now we're going to fill in the rest of the field
    do day = 1, newGrid%noDates
      numMissingClimatology = 0
      numMissingOperational = 0
      call SliceGriddedData ( operational, operMapped, &
        & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
        & newGrid%szas, meanDates(day:day), missingValue=newGrid%missingValue )
      call SliceGriddedData ( climatology, cliMapped, &
        & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
        & newGrid%szas, meanDates(day:day), missingValue=newGrid%missingValue )
      do sza = 1, newGrid%noSzas
        do lst = 1, newGrid%noLsts
          do lon = 1, newGrid%noLons
            do lat = 1, newGrid%noLats
              do surf = 1, newGrid%noHeights
                ! Get the values
                ! cliVal = cliMapped ( surf, lat, lon, lst, sza, day )
                ! opVal = operMapped ( surf, lat, lon, lst, sza, day )
                cliVal = cliMapped ( surf, lat, lon, lst, sza, 1 )
                opVal = operMapped ( surf, lat, lon, lst, sza, 1 )
                ! Weight them by height
                z = scaleHeight * ( 3.0 - log10 ( newGrid%heights(surf) ) )
                if ( scale /= 0.0 ) then
                  cliWeight = ( z - z1 ) / ( z2-z1 )
                  opWeight = 1.0 - cliWeight
                end if

                cliWeight = min ( max ( cliWeight, 0.0_r8 ), 1.0_r8 )
                opWeight = min ( max ( opWeight, 0.0_r8 ), 1.0_r8 )

                ! Check for bad data in operational dataset
                if ( EssentiallyEqual ( opVal, real(newGrid%missingValue, r8) ) ) then
                  opWeight = 0.0
                  numMissingOperational = numMissingOperational + 1
                endif

                ! Check for bad data in the climatology
                if ( EssentiallyEqual ( cliVal, real(newGrid%missingValue, r8) ) ) then
                  call MLSL2Message ( MLSMSG_Error, ModuleName, &
                  & 'There is a bad data point in the climatology field' )
                  numMissingClimatology = numMissingClimatology + 1
                endif

                totalWeight = cliWeight + opWeight
                if ( totalWeight == 0.0 ) then
                  ! Presumably was in the region where operational was supposed
                  ! to dominate, but it's bad, so switch to a priori
                  cliWeight = 1.0
                  totalWeight = 1.0
                end if

                ! OK, store this value
                newGrid%field ( surf, lat, lon, lst, sza, day ) = &
                  & ( cliWeight*cliVal + opWeight*opVal ) / totalWeight
              end do
            end do
          end do
        end do
      end do
      if ( DEEBUG ) then
        call outputNamedValue( 'day ', day )
        call outputNamedValue( 'numMissingClimatology ', numMissingClimatology )
        call outputNamedValue( 'numMissingOperational ', numMissingOperational )
        if ( day > 1 ) cycle
        call dump ( newGrid%field(:,1:10,1:10,1,1,day), &
          & '    gridded field values (1st solar time) =' , &
          & FillValue=newGrid%MissingValue )
        call dump ( newGrid%field(:,1:10,1,:,1,day), &
          & '    gridded field values (1st longitude) =' , &
          & FillValue=newGrid%MissingValue )
        call dump ( newGrid%field(:,1,1:10,:,1,day), &
          & '    gridded field values (1st latitude) =' )
        call dump ( newGrid%field(1,1:10,1:10,:,1,day), &
          & '    gridded field values (1st height) =' , &
          & FillValue=newGrid%MissingValue )
      endif
    end do

    ! Tidy up
    ! Oh, sure, you're careful to account for the memory eaten up by meanDates,
    ! but what about the two biggies???
    call Deallocate_test ( meanDates, 'meanDates', ModuleName )
    s = byte_size(cliMapped)
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(cliMapped(1,1,1,1,1,1)), addr)
    deallocate ( cliMapped, stat=status )
    call test_deallocate ( status, moduleName, 'climapped', s, address=addr )
    s = byte_size(operMapped)
    if ( s > 0 ) addr = transfer(c_loc(operMapped(1,1,1,1,1,1)), addr)
    deallocate ( operMapped, stat=status )
    call test_deallocate ( status, moduleName, 'opermapped', s, address=addr )
    call finishUp ( done = .true. )

  contains

    subroutine FinishUp ( done )
      logical, optional, intent(in) :: done
      logical :: myDone
      myDone = .false.
      if ( present(done) ) myDone = done
      if ( DEEBUG ) then
        call output( 'height: ', advance='no' )
        call output( height, advance='yes' )
        call output( 'scale: ', advance='no' )
        call output( scale, advance='yes' )
        call output( 'operational%verticalCoordinate: ', advance='no' )
        call output( operational%verticalCoordinate, advance='yes' )
        call dump( operational%field( :, 1, 1, 1, 1, 1 ), 'op T' )
        call dump( operational%heights, 'op h' )
        call blanks(3)
        call output( v_is_pressure, advance='yes' )
        call output( 'climatology%verticalCoordinate: ', advance='no' )
        call output( climatology%verticalCoordinate, advance='yes' )
        call dump( climatology%field( :, 1, 1, 1, 1, 1 ), 'cl T' )
        call dump( climatology%heights, 'cl h' )
        call blanks(3)
        call output( v_is_pressure, advance='yes' )
      end if
      if ( myDone ) call trace_end ( "MergeOneGrid", cond=toggle(gen) )
    end subroutine FinishUp

  end function MergeOneGrid

  ! --------------------------------------------  wmoTropFromGrid  -----
  type (griddedData_T) function wmoTropFromGrid ( root, griddedDataBase ) &
    & result ( newGrid )
    use Dump_0, only: Dump
    use GriddedData, only: GriddedData_T, Dump, V_Is_Pressure, V_Is_Eta, &
      & NullifyGriddedData, &
      & DoGriddedDataMatch, &
      & SetupnewGriddedData
    use Init_Tables_Module, only: F_A, F_B, F_Grid
    use MLSCommon, only: DefaultundefinedValue
    use MLSFillValues, only: IsfillValue, RemovefillValues
    use MLSStats1, only: MLSMin, MLSMax, MLSMean
    use MLSStrings, only: Lowercase
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Decoration
    use WMOTropopause, only: ExtraTropics, Twmo
    ! Implements the algorithm published in GRL

    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! This routine creates a new gridded data by finding wmo Tropopause
    ! pressure levels among the temperatures of another gridded data,
    ! possibly using the vertical coords of that gridded data
    ! (if entered via the "grid=" field) 
    ! or else using the corresponding pressures stored
    ! as the field values of a second gridded data
    ! (if entered via "a=first_gridded, b=second_gridded")
    ! The new gridded has only one "level" per horizontal grid point
    ! with tropopause pressures stored in the values field
    
    ! We'll assume some things:
    ! (1) verticalCoordinate is either pressure or eta
    ! (2) If Pressures grid supplied, grids match
    !     and Pressures units either Pa or hPa
    
    ! Local variables
    integer :: field
    integer :: field_index
    real, dimension(:), pointer :: h ! hPa
    integer :: i
    integer :: iDate
    integer :: iLst
    integer :: invert
    integer :: iSza
    integer :: lat
    integer :: lon
    integer :: Me = -1               ! String index for trace
    real :: missingValue
    integer :: nLev
    integer :: nValid
    real, dimension(:), pointer :: p ! Pa
    real, parameter :: pliml = 65.*100 ! in Pa
    real, parameter :: plimlex = 65.*100 ! in Pa
    real, parameter :: plimu = 550.*100 ! in Pa
    type (griddedData_T), pointer :: Placeholder    => null()
    type (griddedData_T), pointer :: Pressures    => null()
    integer :: returnStatus
    real :: scale
    integer :: son
    real, dimension(:), pointer :: t
    type (griddedData_T), pointer :: Temperatures => null()
    real :: trp
    integer :: value
    real, dimension(:), pointer :: xyTemp, xyPress
    logical, parameter :: DEEBUG = .false.
    integer, parameter :: hPa2Pa  = 100 ! Factor convert hPa to Pa

    ! Executable code
    call trace_begin ( me, "wmoTropFromGrid", root, cond=toggle(gen) )
    nullify( xyTemp, xyPress, h, p, t )
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    MISSINGVALUE = REAL( DEFAULTUNDEFINEDVALUE )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_a ) 
        Temperatures => griddedDataBase ( decoration ( decoration ( value ) ) )
        ! Did we defer reading b?
        if ( Temperatures%empty .and. Temperatures%deferReading ) then
          call readGriddedData ( Temperatures%sourceFileName, son, Temperatures%description, &
            & Temperatures%verticalCoordinate, Temperatures, returnStatus, &
            & Temperatures%dimList, TRIM(Temperatures%fieldNames), Temperatures%missingValue )
        endif
      case ( f_b ) 
        Pressures    => griddedDataBase ( decoration ( decoration ( value ) ) )
        if ( Pressures%empty .and. Pressures%deferReading ) then
          call readGriddedData ( Pressures%sourceFileName, son, Pressures%description, &
            & Pressures%verticalCoordinate, Pressures, returnStatus, &
            & Pressures%dimList, TRIM(Pressures%fieldNames), Pressures%missingValue )
        endif
      case ( f_grid ) 
        Temperatures => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do

    if ( associated(Temperatures) .and. associated(Pressures) ) then
    ! What if you reversed the sense of a (Temperatures) and b (Pressures)?
    ! We must switch them
      if ( index( 'hpa,mb', lowercase(trim(Temperatures%units)) ) > 0 .and. &
        &  index( 'hpa,mb', lowercase(trim(Pressures%units)) ) < 1 ) then
        Placeholder  => Pressures
        Pressures    => Temperatures
        Temperatures => Placeholder
      endif
      if ( Temperatures%empty ) then
        call MLSL2Message ( MLSMSG_Warning, moduleName, &
          & 'Empty Temperatures grid for calculating wmo tropopause' )
        call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
        return
      elseif ( Pressures%empty ) then
        call MLSL2Message ( MLSMSG_Warning, moduleName, &
          & 'Empty Pressures grid for calculating wmo tropopause' )
        call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
        return
      endif
    ! What if Temperatures and Pressures don't match
      if ( .not. doGriddeddataMatch( Temperatures, Pressures ) ) then
        call output( 'Vert. coords match? ', advance='no' )
        call output( Temperatures%verticalCoordinate==Pressures%verticalcoordinate, &
          & advance='yes' )
        call output( 'starting dates match? ', advance='no' )
        call output( all( Temperatures%DateStarts==Pressures%DateStarts ), &
          & advance='yes' )
        call output( 'ending dates match? ', advance='no' )
        call output( all( Temperatures%DateEnds==Pressures%DateEnds ), &
          & advance='yes' )
        call output( 'Lsts match? ', advance='no' )
        call output( all( Temperatures%Lsts==Pressures%Lsts ), &
          & advance='yes' )
        call output( 'Szas match? ', advance='no' )
        call output( all( Temperatures%Szas==Pressures%Szas ), &
          & advance='yes' )
        call output( 'Num Heights match? ', advance='no' )
        call output( Temperatures%noHeights==Pressures%noHeights, &
          & advance='yes' )
        call output( 'Heights match? ', advance='no' )
        call output( all( Temperatures%Heights==Pressures%Heights ), &
          & advance='yes' )
        call output( 'Missing value match? ', advance='no' )
        call output( Temperatures%missingValue==Pressures%missingValue, &
          & advance='yes' )
        call output( 'Heights units? ', advance='no' )
        call output( Temperatures%heightsUnits==Pressures%heightsUnits, &
          & advance='yes' )
        call output( 'shapes match? ', advance='no' )
        call output( all( shape(Temperatures%field)==shape(Pressures%field) ), &
          & advance='yes' )
        call output( 'Dumping Temperatures grid', advance='yes' )
        call dump( Temperatures, details=0 )
        call output( 'Dumping Pressures grid', advance='yes' )
        call dump( Pressures, details=0 )
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Gridded T,P data must match to calculate wmo Tropopause' )
      endif
    endif
    if ( DEEBUG ) then
    call output( 'Temperatures grid', advance='yes' )
    call dump( Temperatures, details=0 )
    call output( 'Pressures grid', advance='yes' )
    call dump( Pressures, details=0 )
    call output('Max val', advance='no')
    call output(mlsmax( Pressures%field(:,:,:,1,1,1), Pressures%missingValue ), advance='yes')
    call output('Min val', advance='no')
    call output(mlsmin( Pressures%field(:,:,:,1,1,1), Pressures%missingValue ), advance='yes')
    call output('Mean val', advance='no')
    call output(mlsmean( Pressures%field(:,:,:,1,1,1), Pressures%missingValue ), advance='yes')
    endif
    ! call outputNamedValue( 'Temperatures grid empty?', Temperatures%empty )
    ! call outputNamedValue( 'Pressures grid empty?   ', Pressures%empty )
    if ( .not. associated(Temperatures) ) then
      call MLSL2Message ( MLSMSG_Warning, moduleName, &
        & 'No associated Temperatures grid for calculating wmo tropopause' )
      call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
      return
    endif
    call SetupNewGriddedData ( newGrid, source=Temperatures, &
      & noHeights=1, noDates=1 )
    ! Setup the rest of the quantity
    newGrid%sourceFileName     = 'Gridded Temperatures'
    newGrid%quantityName       = 'wmo Tropopause'
    newGrid%description        = 'wmo Tropopause'
    newGrid%units              = 'hPa' ! If we want 'Pa', restore /100 below
    newGrid%verticalCoordinate = v_is_pressure
    newGrid%equivalentLatitude = Temperatures%equivalentLatitude
    newGrid%heights            = missingValue ! Temperatures%missingValue
    nlev = Temperatures%noHeights
    if ( Temperatures%empty ) then
      newGrid%empty = .true.
      call MLSL2Message ( MLSMSG_Warning, moduleName, &
        & 'Temperatures grid was empty' )
      call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
      return
    elseif ( Pressures%empty ) then
      newGrid%empty = .true.
      call MLSL2Message ( MLSMSG_Warning, moduleName, &
        & 'Pressures grid was empty' )
      call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
      return
    endif
    newGrid%lats               = Temperatures%lats
    newGrid%lons               = Temperatures%lons
    newGrid%lsts               = Temperatures%lsts
    newGrid%szas               = Temperatures%szas
    newGrid%dateEnds           = Temperatures%dateEnds(1)
    newGrid%dateStarts         = Temperatures%dateStarts(1)
    newGrid%missingValue       = MISSINGVALUE / hPa2Pa
    newGrid%field              = MISSINGVALUE
    ! call outputNamedValue( 'Temperatures grid empty?', Temperatures%empty )
    ! call outputNamedValue( 'Pressures grid empty?   ', Pressures%empty )
    if ( nlev < 2 ) then
      call MLSL2Message ( MLSMSG_Warning, moduleName, &
        & 'Too few levels on Temperatures grid for calculating wmo tropopause' )
      call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
      return
    endif
    ! Right now we can't read eta levels, only pressures
    ! but when we move to GEOS5 GMAO we'll have no choice:
    ! Must read eta-level files
    if ( .not. any( &
      & Temperatures%verticalCoordinate == (/ v_is_pressure, v_is_eta /) ) ) &
      & call MLSL2Message ( MLSMSG_Error, moduleName, &
        & 'Temperatures must be on eta or pressure calculating wmo tropopause' )
    if (  Temperatures%verticalCoordinate /= v_is_pressure  .and. &
      & .not. associated(Pressures) ) &
      & call MLSL2Message ( MLSMSG_Error, moduleName, &
        & 'Temperatures illegal verticalcoordinate calculating wmo tropopause' )
    ! For now, just crudely assume the heights is in units of Pa
    ! If not, we'll need to check heightsUnits
    if ( lowercase(Temperatures%heightsUnits) /= 'pa' .and. &
      & .not. associated(Pressures) ) &
      & call MLSL2Message ( MLSMSG_Error, moduleName, &
        & 'Temperatures illegal heightsUnits calculating wmo tropopause' )

    call Allocate_test (h, nlev, 'h', ModuleName )
    call Allocate_test (p, nlev, 'p', ModuleName )
    call Allocate_test (t, nlev, 't', ModuleName )
    ! check vertical orientation of data
    ! (twmo expects ordered from top downward)

    h = Temperatures%heights
    if (h(1) .gt. h(2)) then
      invert=1
      p = h(nlev:1:-1)*100.  ! hPa > Pa
    else
      invert=0
      p = h(:)*100.         ! hPa > Pa
    endif
    
    ! Now actually calculate the tropopause
    ! for every "horizontal" point
    do idate=1, 1 ! size( Temperatures%field, 6 )
      do iSza=1, size( Temperatures%field, 5 )
        do iLst=1, size( Temperatures%field, 4 )
          do lon=1, size( Temperatures%field, 3 )
            do lat=1, size( Temperatures%field, 2 )
              ! Do we have a second, pressure, gridded data?
              if ( associated(Pressures) ) then
                select case (lowercase(Pressures%units))
                case ('pa', 'b')
                  scale = 1. 
                case ('hpa', 'mb')
                  scale = 100. ! To convert hPa to Pa
                case default
                  call output( 'Pressures%units: ', advance='no' )
                  call output( trim(Pressures%units), advance='yes' )
                  call MLSL2Message ( MLSMSG_Error, moduleName, &
                    & 'Pressures units must be Pa, hPa, or mb calculating wmo tropopause' )
                end select
                h = Pressures%field(1:nlev,lat,lon,iLst,iSza,idate)
                if (h(1) .gt. h(2)) then
                  invert=1
                  p = Pressures%field(nlev:1:-1,lat,lon,iLst,iSza,idate) * scale
                else
                  invert=0
                  p = Pressures%field(1:nlev,lat,lon,iLst,iSza,idate) * scale
                endif
              endif
              if ( invert == 1 ) then
                t = temperatures%field(nlev:1:-1,lat,lon,iLst,iSza,idate)
              else
                t = temperatures%field(:,lat,lon,iLst,iSza,idate)
              endif
              where ( t < 0. .or. t > 100000. )
                t = MissingValue
              end where
              nvalid = count( .not. isFillValue(t) )
              if ( nvalid < 2 ) cycle
              call Allocate_test (xyTemp, nvalid, 'xyTemp', ModuleName )
              call Allocate_test (xyPress, nvalid, 'xyPress', ModuleName )
              call RemoveFillValues( t, MISSINGVALUE, xyTemp, &
                & p, xyPress )
              call twmo(nvalid, xyTemp, xyPress, plimu, pliml, trp)
              if ( lon == 1 .and. lat == 1 .and. DEEBUG ) then
                   call dump(xyTemp, 'xyTemp')
                   call dump(xyPress, 'xyPress')
                   call output( 'plimu, pliml, trp', advance='no' )
                   call output( (/ plimu, pliml, trp /), advance='yes' )
              endif
              ! Don't let tropopause sink too low in "extra tropics"
              if ( trp < plimlex .and. &
                & extraTropics(temperatures%lats(lat)) )  &
                & trp = MISSINGVALUE
              if ( trp > 0. .and. trp < 100000000. ) &
                & newGrid%field(1, lat,lon,iLst,iSza,idate) = trp/hPa2Pa !  for 'Pa'
              call Deallocate_test ( xyTemp, 'xyTemp', ModuleName )
              call Deallocate_test ( xyPress, 'xyPress', ModuleName )
            enddo ! Lats
          enddo ! Lons
        enddo ! Lsts
      enddo ! Szas
    enddo ! dates
    call Deallocate_test ( h, 'h', ModuleName )
    call Deallocate_test ( p, 'p', ModuleName )
    call Deallocate_test ( t, 't', ModuleName )
    if ( DeeBug ) then
      call output('scale: ', advance='no')
      call output(scale, advance='no')
    call output('Max val', advance='no')
    call output(mlsmax( newGrid%field(:,:,:,1,1,1), newGrid%missingValue ), advance='yes')
    call output('Min val', advance='no')
    call output(mlsmin( newGrid%field(:,:,:,1,1,1), newGrid%missingValue ), advance='yes')
    call output('Mean val', advance='no')
    call output(mlsmean( newGrid%field(:,:,:,1,1,1), newGrid%missingValue ), advance='yes')
      stop
    endif
    if ( DEEBUG ) call dump( newGrid )
    call trace_end ( "wmoTropFromGrid", cond=toggle(gen) )
  end function wmoTropFromGrid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MergeGridsModule

! $Log$
! Revision 2.71  2020/07/23 22:18:40  pwagner
! AllowEmptyGrids is a new field in Concatenate; defaults to pre-v5.1 behavior
!
! Revision 2.70  2020/07/22 22:56:16  pwagner
! Made more cmds public; AllowEmptyGrids now TRUE when Concatenating grids
!
! Revision 2.69  2020/07/09 23:54:30  pwagner
! Many cmds from readApriori and MergeGrids phase now available to Fill phase
!
! Revision 2.68  2019/10/03 17:30:17  pwagner
! Convert from eta levels may now take a vGrid field
!
! Revision 2.67  2019/09/23 20:39:04  pwagner
! Conversion from Eta surfaces may optionally be logarithmic
!
! Revision 2.66  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.65  2018/03/22 18:15:34  pwagner
! Added command IsFileAbsent; may occur in ReadApriori, MergeGrids, and Output sections
!
! Revision 2.64  2018/03/14 22:45:51  pwagner
! May changeSettings in readApriori and MergeGrids sections
!
! Revision 2.63  2017/11/03 20:59:27  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.62  2017/03/17 00:38:10  pwagner
! Quit with apt message if ConvertEtaToP has no climatology file to use
!
! Revision 2.61  2016/04/01 00:27:15  pwagner
! May now Execute a single command or a script of lines from l2cf
!
! Revision 2.60  2015/03/28 02:49:58  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.59  2014/09/05 01:15:02  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Track
! allocate/deallocate size in bytes instead of Memory_Units.
!
! Revision 2.58  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.57  2014/06/20 20:30:24  pwagner
! Less bebug-type printing
!
! Revision 2.56  2014/06/11 20:03:28  pwagner
! New concatenateGrids command; grid field in Concatenated renamed sourceGrid
!
! Revision 2.55  2014/06/04 18:38:07  pwagner
! Many steps to conserve memory; accurately account for its usage
!
! Revision 2.54  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.53  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.52  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.51  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.50  2013/10/09 23:41:55  vsnyder
! Add Evaluate_Variable
!
! Revision 2.49  2013/08/31 02:29:12  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.48  2013/08/30 02:45:44  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.47  2012/08/16 18:01:04  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.46  2012/06/07 22:48:43  pwagner
! Printing during mergeGrids now requires 'grid' switch
!
! Revision 2.45  2012/05/08 17:49:54  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.44  2011/05/09 18:23:23  pwagner
! description field now marks result of convert, concatenate
!
! Revision 2.43  2011/04/27 17:40:22  pwagner
! Added new command (not a named spec) wmoTropFromGrids
!
! Revision 2.42  2011/04/20 16:51:55  pwagner
! Added new flexibility to l2cf control flow by run-time booleans
!
! Revision 2.41  2009/12/14 18:37:50  pwagner
! Dont crash in wmoTropFromGrid if one of the grids is empty
!
! Revision 2.40  2009/11/05 00:29:06  pwagner
! Better diagnostics output if something goes wrong
!
! Revision 2.39  2009/10/26 17:12:09  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.38  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.37  2008/06/06 22:52:53  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.36  2007/12/07 01:15:01  pwagner
! Removed unused dummp varaibles, etc.
!
! Revision 2.35  2007/10/24 00:17:22  pwagner
! Removed unused declarations
!
! Revision 2.34  2007/08/17 00:33:59  pwagner
! May now dump, stop, under control of l2cf
!
! Revision 2.33  2007/07/04 01:44:15  vsnyder
! Actually leave early from ConvertEtaToP
!
! Revision 2.32  2007/07/04 01:08:43  pwagner
! trace_begin and _end rebalanced in ConvertEtaToP
!
! Revision 2.31  2007/06/29 21:01:14  pwagner
! Should print less unless debugging
!
! Revision 2.30  2007/06/07 21:56:00  pwagner
! Prevents another cause of crashing; extra debugging
!
! Revision 2.29  2007/03/23 00:27:17  pwagner
! Valiant attempts to bring two Lahey versions results closer
!
! Revision 2.28  2007/01/12 00:34:04  pwagner
! Renamed routine outputNamedValue
!
! Revision 2.27  2006/11/03 19:40:30  pwagner
! Fixed unassociated pointers NAG caught
!
! Revision 2.26  2006/11/03 00:25:47  pwagner
! Fixed bug in tropopause calculation
!
! Revision 2.25  2006/11/01 20:34:12  pwagner
! hasty fix to wmo tropopause
!
! Revision 2.24  2006/07/07 23:10:03  pwagner
! Should handle missing GEOS5 files with greater grace
!
! Revision 2.23  2006/06/22 00:20:46  pwagner
! Repair cosmetic blemishes when tracing
!
! Revision 2.22  2006/06/15 17:36:30  pwagner
! Bail out of Concatenating if geos5 files missing
!
! Revision 2.21  2006/06/15 00:02:33  pwagner
! Should work with geos5: convert then concatenate
!
! Revision 2.20  2006/06/13 22:13:12  pwagner
! changed interface to ConvertFromEtaLevelGrids
!
! Revision 2.19  2006/05/12 21:26:37  pwagner
! Added extra debugging statements
!
! Revision 2.18  2006/05/09 16:42:02  pwagner
! May find wmo p trop with two eta-level grids
!
! Revision 2.17  2006/05/04 23:04:59  pwagner
! May convertEtaToP and create a VGrid in MergeGrids section
!
! Revision 2.16  2006/02/11 00:14:08  pwagner
! May calculate wmoTropopause in this section directly
!
! Revision 2.15  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.13  2003/06/06 01:06:59  livesey
! Added DeleteGrids stuff
!
! Revision 2.12  2003/05/09 02:13:05  livesey
! Removed a dump
!
! Revision 2.11  2003/05/09 01:55:14  livesey
! Sped up Merge by using new SliceGriddedData routine.
!
! Revision 2.10  2003/04/04 00:08:26  livesey
! Added Concatenate capability, various reorganizations.
!
! Revision 2.9  2003/02/28 02:33:28  livesey
! Bug fix, careless with the old emacs.
!
! Revision 2.8  2003/02/28 02:25:50  livesey
! First working version.
!
! Revision 2.7  2003/02/19 19:15:13  pwagner
! Consistent with new GriddedData_T and rgr
!
! Revision 2.6  2002/11/22 12:21:14  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.5  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/08/22 20:26:09  vsnyder
! Move another USE from module scope to procedure scope
!
! Revision 2.3  2002/08/21 02:23:39  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.2  2002/01/26 00:10:54  livesey
! It compiles at least
!
! Revision 2.1  2002/01/24 00:58:03  livesey
! First version, not much more than a stub
!
