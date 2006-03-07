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
MODULE Construct                ! The construct module for the MLS L2 sw.
!=============================================================================

  ! This module performs the `construct' task for the level 2 software.  This
  ! task involves constructing templates for vector quantities, vectors and
  ! matrices.

  implicit none

  private

  public :: MLSL2Construct, MLSL2DeConstruct, ConstructMIFGeolocation
  
  !------------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName="$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------- BooleanFromAnyGoodRadiances --
  function BooleanFromAnyGoodRadiances ( root, chunk, filedatabase ) &
    & result(hashsize)
    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use ConstructQuantityTemplates, only: AnyGoodSignalData
    use Chunks_m, only: MLSCHUNK_T
    use Dump_0, only: Dump
    use INIT_TABLES_MODULE, only: F_SIGNAL, F_Boolean
    use MLSCommon, only: r8, MLSFile_T
    use MLSL2Options, only: runTimeValues
    use MLSSignals_m, only: GetSignalName, &
      & SIGNALS
    use MLSStringLists, only: NumStringElements, PutHashElement, &
      & SwitchDetail
    use MLSStrings, only: lowerCase
    use output_m, only: output
    use Parse_signal_m, only: Parse_signal
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (MLSChunk_T), intent(in) :: chunk
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer             :: hashsize
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: formula
    integer :: keyNo
    character(len=32) :: nameString
    integer :: s
    integer :: signalIndex
    integer, pointer :: Signal_Indices(:)         ! Indices in the signals
    character(len=32) :: signalString
    integer :: son
    character(len=32) :: subSignalString
    logical :: tvalue
    ! Executable
    nullify(Signal_Indices)
    ! call get_string(name, nameString)
    ! nameString = lowerCase(nameString)
    signalString = ' '
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_signal )
        call get_string( sub_rosa(subtree(2,son)), signalString, strip=.true. )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    if ( signalString /= ' ' ) then
      if ( switchDetail(switches, 'bool') > 0 ) &
        & call output( 'signal: ' // trim(signalString), advance='yes' )
      call Parse_signal(signalString, signal_indices)
      tvalue = .false.
      ! Loop over signals, or-ing them until we get TRUE
      do s=1, size(signal_indices)
        signalIndex = signal_indices(s)
        if ( switchDetail(switches, 'bool') > 0 ) then
          call GetSignalName ( signalIndex, subSignalString, &                   
            & sideband=signals(signalIndex)%sideband, noChannels=.TRUE. )
          call output( 'sub-signal: ' // trim(subSignalString), advance='yes' )
        endif
        tvalue = tvalue .or. &
          & AnyGoodSignalData ( signalIndex, signals(signalIndex)%sideband, &
          & filedatabase, chunk )
        if ( tvalue ) then
          if ( switchDetail(switches, 'bool') > 0 ) then
            call output( 'good signal data found: ' &
              & // trim(subSignalString), advance='yes' )
          endif
          exit
        endif
      enddo
      call deallocate_test(Signal_Indices, 'Signal_Indices', ModuleName)
    else
      print *, 'Sorry-unable to parse ', trim(signalString)
      tvalue = .false.
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    hashsize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodRadiances

  ! ------------------------------------- BooleanFromAnyGoodValues --
  function BooleanFromAnyGoodValues ( root, vectors ) result(size)
    use Dump_0, only: Dump
    use Expr_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_PRECISION, F_QUALITY, &
      & F_QUANTITY, F_Boolean, F_STATUS
    use ManipulateVectorQuantities, only: AnyGoodDataInQty
    use MLSCommon, only: r8, rv
    use MLSL2Options, only: runTimeValues
    use MLSStringLists, only: BooleanValue, NumStringElements, PutHashElement, &
      & SwitchDetail
    use MLSStrings, only: lowerCase
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
    use VectorsModule, only: Vector_T, VectorValue_T, &
      & GetVectorQtyByTemplateIndex
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (vector_T), dimension(:), pointer :: Vectors
    integer             :: size
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: formula
    integer :: keyNo
    character(len=32) :: nameString
    type (vectorValue_T), pointer :: PRECISIONQUANTITY
    integer :: QUANTITYINDEX
    real(rv) :: QUALITY_MIN
    type (vectorValue_T), pointer :: QUALITYQUANTITY
    type (vectorValue_T), pointer :: Quantity
    integer :: son
    integer :: source
    type (vectorValue_T), pointer :: STATUSQUANTITY
    logical :: tvalue
    integer :: VECTORINDEX
    ! Executable
    nullify( precisionquantity, qualityquantity, Quantity, statusquantity )
    ! call get_string(name, nameString)
    ! nameString = lowerCase(nameString)
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)
      source = subtree(2,son) ! required to be an n_dot vertex

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_precision )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        precisionQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case ( f_quality )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        qualityQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case ( f_quantity )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        Quantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case ( f_status )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        statusQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    tvalue = AnyGoodDataInQty ( a=Quantity, &
      & precision=precisionQuantity, quality=qualityQuantity, &
      & status=statusQuantity, quality_min=quality_min )
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodValues

  ! ------------------------------------- DealWithBooleanFromMLSCfInfo --
  function DealWithBooleanFromMLSCfInfo ( name, root ) result(size)
    ! Called either when a Boolean is first declared
    ! syntax: 
    ! name: Boolean, formula="formula"
    !
    ! or when it is reevaluated
    ! syntax: 
    ! Reevaluate, formula="formula", Boolean="name"
    use Dump_0, only: Dump
    use Expr_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_FORMULA, F_VALUES
    use MLSCommon, only: r8
    use MLSL2Options, only: runTimeValues
    use MLSStringLists, only: BooleanValue, NumStringElements, PutHashElement, &
      & SwitchDetail
    use MLSStrings, only: lowerCase
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: name
    integer, intent(in) :: root
    integer             :: size
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: formula
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    integer, dimension(2) :: UNITASARRAY ! From expr
    real(r8), dimension(2) :: VALUEASARRAY ! From expr
    ! Executable
    tvalue= .false.
    if ( name > 0 ) then
      call get_string(name, nameString)
      nameString = lowerCase(nameString)
    endif
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string ( sub_rosa(subtree(2,son)), nameString, strip=.true. )
        nameString = lowerCase(nameString)
      case ( f_formula )
        call get_string ( sub_rosa(subtree(2,son)), formula, strip=.true. )
        tvalue = BooleanValue (formula, runTimeValues%lkeys, runTimeValues%lvalues)
      case ( f_values )
        call expr ( son , unitAsArray, valueAsArray )
        tvalue = ( valueAsArray(1) /= 0 )
        ! badRange = valueAsArray
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function DealWithBooleanFromMLSCfInfo

  ! --------------------------------------------- ConstructMIFGeolocation --
  subroutine ConstructMIFGeolocation ( mifGeolocation, filedatabase, chunk )
    ! mifGeolocation is just quantity templates containing geolocation
    ! information for the GHz and THz modules.  The software can then
    ! point to these for geolocation information for all minor frame
    ! quantities saving file IO and memory.
    use Chunks_m, only: MLSCHUNK_T
    use ConstructQuantityTemplates, only: ConstructMinorFrameQuantity
    use QuantityTemplates, only: QUANTITYTEMPLATE_T
    use MLSCommon, only: MLSFile_T
    use MLSSignals_m, only: MODULES
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate

    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (MLSChunk_T), intent(in) :: chunk
    
    ! Local variables
    integer :: INSTRUMENTMODULEINDEX    ! Loop counter
    integer :: STATUS                   ! Flag

    if ( .not. associated ( mifGeolocation ) ) then
      ! Don't overwrite it if we already have it, e.g. from previous construct
      ! or forge.
      allocate ( mifGeolocation(size(modules)), STAT=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//"mifGeolocation" )
    
      ! Now try to fill it if we have any L1BFiles
      if (associated(filedatabase) ) then
        do instrumentModuleIndex = 1, size(modules)
          call ConstructMinorFrameQuantity ( filedatabase, chunk, &
            & instrumentModuleIndex, mifGeolocation(instrumentModuleIndex) )
        end do
      else
        mifGeolocation%noSurfs = 0
        mifGeolocation%noInstances = 0
      end if
    end if
  end subroutine ConstructMIFGeolocation

  ! ---------------------------------------------  MLSL2Construct  -----
  subroutine MLSL2Construct ( root, filedatabase, processingRange, chunk, &
       & quantityTemplatesBase, vectorTemplates, vectors, FGrids, HGrids, &
       & l2gpDatabase, ForwardModelConfigDatabase, mifGeolocation )

  ! This is the `main' subroutine for this module

    use Chunks_m, only: MLSChunk_T
    use ConstructQuantityTemplates, only: &
      & CreateQtyTemplateFromMLSCfInfo, ForgeMinorFrames
    use ConstructVectorTemplates, only: CreateVecTemplateFromMLSCfInfo
    use DumpCommand_m, only: DumpCommand
    use FGrid, only: FGrid_T
    use ForwardModelConfig, only: AddForwardModelConfigToDatabase, &
      & ForwardModelConfig_T
    use ForwardModelSupport, only: ConstructForwardModelConfig
    use HGridsDatabase, only: ADDHGRIDTODATABASE, HGRID_T
    use HGrid, only: CREATEHGRIDFROMMLSCFINFO
    use INIT_TABLES_MODULE, only: S_ANYGOODVALUES, S_ANYGOODRADIANCES, &
      & S_BOOLEAN, S_DUMP, &
      & S_FORGE, S_FORWARDMODEL, S_HGRID, &
      & S_PHASE, S_QUANTITY, S_REEVALUATE, S_TIME, S_VECTORTEMPLATE
    use L2GPData, only: L2GPDATA_T
    use MLSCommon, only: MLSFile_T, TAI93_Range_T
    use MLSL2Options, only: RESTARTWARNINGS, SPECIALDUMPFILE
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, addPhaseToPhaseNames
    use MoreTree, only: Get_Spec_ID
    use OUTPUT_M, only: BLANKS, OUTPUT, RESUMEOUTPUT, &
      & revertoutput, SUSPENDOUTPUT, switchOutput
    use QuantityTemplates, only: AddQuantityTemplateToDatabase, &
      & QuantityTemplate_T
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, NODE_ID, NSONS, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_NAMED
    use VectorsModule, only: AddVectorTemplateToDatabase, &
      & Vector_T, VectorTemplate_T

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Root of the tree for the Construct section
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (TAI93_Range_T), intent(in) :: processingRange
    type (MLSChunk_T), intent(in) :: chunk
    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplatesBase
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (vector_T), dimension(:), pointer :: Vectors
    type (FGrid_T), dimension(:), pointer :: fGrids
    type (HGrid_T), dimension(:), pointer :: HGrids
    type (L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (ForwardModelConfig_T), dimension(:), pointer :: ForwardModelConfigDatabase
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables

    integer :: I                ! Loop counter
    integer :: KEY              ! S_... from Init_Tables_Module.
    integer :: NAME             ! Sub-rosa index of name
    integer :: SON              ! Son or grandson of Root
    REAL :: T1, T2              ! for timing
    logical :: TIMING

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )

    if ( toggle(gen) ) call trace_begin ( "MLSL2Construct", root )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )

    ! First we're going to setup our mifGeolocation quantityTemplates.
    call ConstructMIFGeolocation ( mifGeolocation, filedatabase, chunk )

    ! The rest is fairly simple really.  We just loop over the mlscf 
    ! instructions and hand them off to people

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        name = sub_rosa(subtree(1,son))
      else ! Son is n_spec_args
        key = son
        name = 0
      end if

      ! Node_id(key) is now n_spec_args.

      select case( get_spec_id(key) )
      case ( s_anygoodvalues )
        call decorate ( key, &
          & BooleanFromAnyGoodValues ( key, vectors ) )
      case ( s_anygoodradiances )
        call decorate ( key, &
          & BooleanFromAnyGoodRadiances ( key, chunk, filedatabase ) )
      case ( s_dump )
        call dumpCommand ( key, quantityTemplatesBase, &
          & vectorTemplates, forwardModelConfigs=forwardModelConfigDatabase, &
          & hGrids=hGrids )
      case ( s_forge )
        call ForgeMinorFrames ( key, chunk, mifGeolocation )
      case ( s_forwardModel )
        call decorate ( key, AddForwardModelConfigToDatabase ( &
          & forwardModelConfigDatabase, &
          & ConstructForwardModelConfig ( name, key, .false. ) ) )
      case ( s_hgrid )
        call decorate ( key, AddHGridToDatabase ( hGrids, &
          & CreateHGridFromMLSCFInfo ( name, key, filedatabase, l2gpDatabase, &
          & processingRange, chunk ) ) )
      case ( s_phase )
        call addPhaseToPhaseNames ( name, key )

      case ( s_quantity )
        call decorate ( key, AddQuantityTemplateToDatabase ( &
          & quantityTemplatesBase, CreateQtyTemplateFromMLSCfInfo ( name, key, &
            & fGrids, hGrids, filedatabase, chunk, mifGeolocation ) ) )
      case ( s_Boolean )
        call decorate ( key,  DealWithBooleanFromMLSCfInfo ( name, key ) )
      case ( s_Reevaluate )
        call decorate ( key,  DealWithBooleanFromMLSCfInfo ( 0, key ) )
      case ( s_vectortemplate )
        call decorate ( key, AddVectorTemplateToDatabase ( vectorTemplates, &
          & CreateVecTemplateFromMLSCfInfo ( name, key, quantityTemplatesBase ) ) )
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

    if ( specialDumpFile /= ' ' ) &
      & call revertOutput
    if ( toggle(gen) ) call trace_end ( "MLSL2Construct" )

    if ( timing ) call sayTIme

  contains
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for MLSL2Construct = " )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Construct

  ! -------------------------------------------  MLSL2DeConstruct  -----
  subroutine MLSL2DeConstruct ( quantityTemplatesBase, vectorTemplates, &
    &                           mifGeolocation, hGrids )

  ! DeConstruct the Quantity and Vector template databases.

    use HGridsDatabase, only: DestroyHGridDatabase, HGrid_T
    use QuantityTemplates, only: DestroyQuantityTemplateDatabase, &
      & QuantityTemplate_T
    use VectorsModule, only: DestroyVectorTemplateDatabase, VectorTemplate_T

    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplatesBase
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation
    type (HGrid_T), dimension(:), pointer :: hGrids

    call destroyVectorTemplateDatabase ( vectorTemplates )
    call destroyQuantityTemplateDatabase ( quantityTemplatesBase )
    call destroyQuantityTemplateDatabase ( mifGeolocation )
    call destroyHGridDatabase ( hGrids )
  end subroutine MLSL2DeConstruct

!=============================================================================
  logical function not_used_here()
  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)), save :: Id = idParm
  !-----------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END MODULE Construct
!=============================================================================

!
! $Log$
! Revision 2.57  2006/03/07 00:51:32  pwagner
! May change already-set Booleans via reevaluate command
!
! Revision 2.56  2006/03/04 00:20:45  pwagner
! May skip retrieval, directWrites depending on runtime Booleans
!
! Revision 2.55  2006/02/10 21:19:07  pwagner
! dumps may go to special dumpfile
!
! Revision 2.54  2006/01/06 01:16:34  pwagner
! silent boolean field can silence selected phases
!
! Revision 2.53  2006/01/04 01:27:11  vsnyder
! Comment out use for unreference L1BInfo_T
!
! Revision 2.52  2005/06/03 02:05:29  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of passing as an argument.
!
! Revision 2.51  2005/05/31 17:51:16  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.50  2005/03/12 00:50:27  pwagner
! May restart warnings counter at each phase
!
! Revision 2.49  2004/10/13 02:24:02  livesey
! Added vGrids to Forge command
!
! Revision 2.48  2004/05/20 19:47:55  vsnyder
! Do all dumping by way of DumpCommand
!
! Revision 2.47  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.46  2004/05/18 01:24:31  vsnyder
! Add HGrids argument to DumpCommand
!
! Revision 2.45  2004/05/11 02:54:36  vsnyder
! Remove USEs and declarations for unreferenced symbols
!
! Revision 2.44  2004/05/03 16:37:11  pwagner
! Get dump from QuantityTemplates module
!
! Revision 2.43  2004/05/01 04:05:26  vsnyder
! Add Dump command
!
! Revision 2.42  2003/10/22 21:17:06  pwagner
! aPhaseName: Phase added to Fill, Construct sections to time phases
!
! Revision 2.41  2003/07/15 18:18:16  livesey
! Change to forward model config call
!
! Revision 2.40  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.39  2003/05/28 04:39:32  livesey
! Removed some obsolete checking
!
! Revision 2.38  2002/10/08 17:36:19  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.37  2002/10/05 00:43:47  livesey
! Split out ConstructMIFGeolocation
!
! Revision 2.36  2002/09/26 23:58:57  livesey
! Bug fix, decorated wrong node
!
! Revision 2.35  2002/09/25 20:07:41  livesey
! Can now construct forward models inside construct
!
! Revision 2.34  2002/09/18 22:49:42  pwagner
! Receives returnStatus from CreateQtyTemplateFromMLSCFInfo
!
! Revision 2.33  2002/08/20 22:43:37  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.32  2001/12/16 00:57:00  livesey
! Now takes all chunks as argument as HGrid needs them
!
! Revision 2.31  2001/12/14 01:42:47  livesey
! Passes processingRange to HGrid construction
!
! Revision 2.30  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.29  2001/10/31 19:07:15  livesey
! Added fGrid stuff
!
! Revision 2.28  2001/10/05 23:30:37  pwagner
! Fixed small bug in destroyingQuantTempldb
!
! Revision 2.27  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.26  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.25  2001/05/21 20:57:55  livesey
! Fixed bug, was overwriting mifGelocation with each new construct.
!
! Revision 2.24  2001/05/12 00:16:55  livesey
! Big fix, only dump hGrids etc. if they exist.
!
! Revision 2.23  2001/05/10 00:43:12  livesey
! Moved ownership of hGrids into tree walker to allow multiple calls
!
! Revision 2.22  2001/04/28 01:43:10  vsnyder
! Improved the timing message
!
! Revision 2.21  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.20  2001/04/25 19:29:31  livesey
! Fixed bug with forge, now sets mafIndex and mafCounter appropriately
!
! Revision 2.19  2001/04/23 23:53:03  livesey
! Tidied up deallocation of minor frame quantities
!
! Revision 2.18  2001/04/23 23:24:55  livesey
! Changed l2gpDatabase to pointer
!
! Revision 2.17  2001/04/21 01:25:23  livesey
! New version, can construct h/v grids from l2gp
!
! Revision 2.16  2001/04/20 23:11:26  livesey
! Added `Forge' stuff
!
! Revision 2.15  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.14  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.13  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.12  2001/03/28 18:09:24  vsnyder
! Cosmetic changes
!
! Revision 2.11  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.10  2001/03/15 21:09:52  vsnyder
! Use 'get_spec_id' from More_Tree
!
! Revision 2.9  2001/03/15 21:03:46  vsnyder
! Cross-references between databases are by database index, not tree index
!
! Revision 2.8  2001/03/03 00:07:51  livesey
! New mifGeolocation stuff
!
! Revision 2.7  2001/03/02 01:28:12  livesey
! Uses new MLSSignals
!
! Revision 2.6  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.5  2001/01/03 17:51:05  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.4  2000/11/16 02:15:14  vsnyder
! More work on timing.
!
! Revision 2.3  2000/11/16 02:00:48  vsnyder
! Added timing.
!
! Revision 2.2  2000/09/11 19:23:48  ahanzel
! Removed extra log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

