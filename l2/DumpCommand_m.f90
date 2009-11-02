! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DumpCommand_M

! Process a "dump" command. Or say whether to "skip" remainder of section.
! Or functions to set a run-time Boolean flag.
! (Should these latter functions be moved into a special module?)

  implicit none
  private

  public :: BooleanFromAnyGoodRadiances
  public :: BooleanFromAnyGoodValues
  public :: BooleanFromCatchWarning
  public :: BooleanFromComparingQtys
  public :: BooleanFromFormula
  public :: DumpCommand, Skip

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------- BooleanFromAnyGoodRadiances --
  function BooleanFromAnyGoodRadiances ( root, chunk, filedatabase ) &
    & result(hashsize)
    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use ConstructQuantityTemplates, only: AnyGoodSignalData
    use Chunks_m, only: MLSCHUNK_T
    use Dump_0, only: Dump
    use INIT_TABLES_MODULE, only: F_SIGNAL, F_Boolean
    use MLSCommon, only: MLSFile_T
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
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
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
        end if
        tvalue = tvalue .or. &
          & AnyGoodSignalData ( signalIndex, signals(signalIndex)%sideband, &
          & filedatabase, chunk )
        if ( tvalue ) then
          if ( switchDetail(switches, 'bool') > 0 ) then
            call output( 'good signal data found: ' &
              & // trim(subSignalString), advance='yes' )
          end if
          exit
        end if
      enddo
      call deallocate_test(Signal_Indices, 'Signal_Indices', ModuleName)
    else
      print *, 'Sorry-unable to parse ', trim(signalString)
      tvalue = .false.
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    hashsize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodRadiances

  ! ------------------------------------- BooleanFromAnyGoodValues --
  function BooleanFromAnyGoodValues ( root, vectors ) result(thesize)
    use Dump_0, only: Dump
    use INIT_TABLES_MODULE, only: F_PRECISION, F_QUALITY, &
      & F_QUANTITY, F_Boolean, F_STATUS
    use ManipulateVectorQuantities, only: AnyGoodDataInQty
    use MLSCommon, only: rv
    use MLSL2Options, only: runTimeValues
    use MLSStringLists, only: NumStringElements, PutHashElement, &
      & SwitchDetail
    use MLSStrings, only: lowerCase
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VectorsModule, only: Vector_T, VectorValue_T, &
      & GetVectorQtyByTemplateIndex
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (vector_T), dimension(:) :: Vectors
    integer             :: thesize
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
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
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodValues

  ! ------------------------------------- BooleanFromCatchWarning --
  function BooleanFromCatchWarning ( root ) result(size)
    ! Called to check if the last command resulted in a warning
    ! (either printed or suppressed)
    ! and optionally if the warning matches a supplied message string
    ! syntax: 
    ! CatchWarning, [message='string'], Boolean="name"
    use Dump_0, only: Dump
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_MESSAGE
    use MLSL2Options, only: runTimeValues
    use MLSMessageModule, only: MLSMessageInquire
    use MLSStringLists, only: NumStringElements, PutHashElement, &
      & SwitchDetail
    use MLSStrings, only: lowerCase
    use OUTPUT_M, only: outputNamedValue
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: root
    integer             :: size
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: LastWarningMsg
    character(len=255) :: message
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    ! Executable
    tvalue= .false.
    message = ' '
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
      case ( f_message )
        call get_string ( sub_rosa(subtree(2,son)), message, strip=.true. )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    call MLSMessageInquire( LastWarningMsg=LastWarningMsg )
    call outputNamedValue( 'message to match', message )
    call outputNamedValue( 'LastWarningMsg', LastWarningMsg )
    if ( len_trim(LastWarningMsg) < 1 ) then
      tvalue = .false.
    elseif ( len_trim(message) < 1 ) then
      tvalue = .true.
    else
      ! tvalue = streq( message, LastWarningMsg, '-wcf' )
      ! This allows partial matches, case-insensitive
      tvalue = &
        & index( lowerCase(LastWarningMsg), lowerCase(trim(adjustl(message))) ) > 0
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > -1 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromCatchWarning

  ! ------------------------------------- BooleanFromComparingQtys --
  function BooleanFromComparingQtys ( root, vectors ) result(thesize)
    use Dump_0, only: Dump
    use Expr_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_A, F_B, F_C, F_Boolean, F_FORMULA
    use MLSCommon, only: r8, rv, DEFAULTUNDEFINEDVALUE
    use MLSL2Options, only: runTimeValues
    use MLSMessageModule, only: MLSMessage, MLSMessageCalls, MLSMSG_error
    use MLSStats1, only: mlsmax, mlsmin, mlsmean, mlsmedian
    use MLSStringLists, only: GetStringElement, NumStringElements, PutHashElement, &
      & ReplaceSubString, SwitchDetail
    use MLSStrings, only: lowerCase
    use String_Table, only: get_string
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VectorsModule, only: Vector_T, VectorValue_T, M_Fill, &
      & GetVectorQtyByTemplateIndex
    ! Dummy args
    ! Called to endow Boolean with result from comparing
    ! (1) Two quantities (a and b), or
    ! (2) A quantity and a constant (a and c)
    ! The comparison op may be one of "<", ">", or "='
    ! and the "flattening" to be taken may be one of
    ! "any", "all", "min", "max", "mean", or "median"
    ! E.g., to compare a=a.qty and b=b.qty, returning true if 
    ! all(a.qty%values > b.qty%values)
    ! formula = "all a > b"

    ! in general we will parse formula as being made up by
    ! "flattening a op [b][c]"
    
    ! Compare, a=a.qty, [b=b.qty], [c=c], formula="formula", Boolean="name"
    integer, intent(in) :: root
    type (vector_T), dimension(:) :: Vectors
    integer             :: thesize
    ! Internal variables
    type (vectorValue_T), pointer :: AQUANTITY
    type (vectorValue_T), pointer :: BQUANTITY
    real(rv) :: A, B, C                       ! constant "c" in formula
    logical, parameter :: countEmpty = .true.
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=8) :: flattening, arg(2), op
    character(len=255) :: formula
    character(len=255) :: formulaTemp
    integer :: keyNo
    character(len=32) :: nameString
    integer :: QUANTITYINDEX
    integer :: son
    integer :: source
    logical :: tvalue
    integer, dimension(2) :: UNITASARRAY ! From expr
    real(r8), dimension(2) :: VALUEASARRAY ! From expr
    integer :: VECTORINDEX
    ! Executable
    call MLSMessageCalls( 'push', constantName=ModuleName//'%BooleanFromComparingQtys' )
    nullify( aQuantity, bQuantity )
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
      case ( f_a )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        aQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case ( f_b )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        bQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
      case(f_c)
        call expr ( source, unitAsArray, valueAsArray )
        c = valueAsArray(1)
      case ( f_formula )
        call get_string ( sub_rosa(subtree(2,son)), formula, strip=.true. )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    ! What kind of flattening, relationship, and is the last arg b or c?
    ! 1st--let's separate args and ops neatly
    call ReplaceSubString( formula, formulaTemp, '(', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( formulaTemp, formula, '&', '(', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( formula, formulaTemp, ')', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( formulaTemp, formula, '&', ')', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( formula, formulaTemp, '==', ' = ', &
      & which='first', no_trim=.true. )
    call ReplaceSubString( formulaTemp, formula, '=', ' = ', &
      & which='first', no_trim=.true. )
    call ReplaceSubString( formula, formulaTemp, '>', ' > ', &
      & which='first', no_trim=.true. )
    call ReplaceSubString( formulaTemp, formula, '<', ' < ', &
      & which='first', no_trim=.true. )
    ! 2nd--go through the elements
    call GetStringElement ( trim(formula), flattening, 1, &
            & countEmpty=.false., inseparator=' ' )
    call GetStringElement ( trim(formula), arg(1), 2, &
            & countEmpty=.false., inseparator=' ' )
    call GetStringElement ( trim(formula), op, 3, &
            & countEmpty=.false., inseparator=' ' )
    call GetStringElement ( trim(formula), arg(2), 4, &
            & countEmpty=.false., inseparator=' ' )
    flattening = lowerCase(flattening)
    arg(1) = lowerCase(arg(1))
    arg(2) = lowerCase(arg(2))
    if ( arg(1) /= 'a' ) then
      call mlsmessage (MLSMSG_Error, moduleName, &
        & 'Formula in compare must be "flattening a op [b][c]".' )
    elseif( index('<>=', trim(op)) < 1 ) then
      call mlsmessage (MLSMSG_Error, moduleName, &
        & 'Formula in compare found unrecognized relation: ' // trim(op) )
    end if
    ! What kind of flattening?
    select case (flattening)
    case ('all')
      if ( .not. associated ( aQuantity%mask ) ) then
        tvalue = all ( isRelation( trim(op), aQuantity%values, c ) )
        if ( arg(2) == 'b' ) &
          & tvalue = all ( isRelation( trim(op), aQuantity%values, bQuantity%values ) )
      else
        tvalue = all ( isRelation( trim(op), aQuantity%values, c ) .or. &
          & iand ( ichar(aQuantity%mask(:,:)), m_fill ) /= 0 )
        if ( arg(2) == 'b' ) &
        tvalue = all ( isRelation( trim(op), aQuantity%values, bQuantity%values ) .or. &
          & iand ( ichar(aQuantity%mask(:,:)), m_fill ) /= 0 )
      end if
    case ('any')
      if ( .not. associated ( aQuantity%mask ) ) then
        tvalue = any ( isRelation( trim(op), aQuantity%values, c ) )
        if ( arg(2) == 'b' ) &
          & tvalue = any ( isRelation( trim(op), aQuantity%values, bQuantity%values ) )
      else
        tvalue = any ( isRelation( trim(op), aQuantity%values, c ) .and. &
          & iand ( ichar(aQuantity%mask(:,:)), m_fill ) == 0 )
        if ( arg(2) == 'b' ) &
        tvalue = any ( isRelation( trim(op), aQuantity%values, bQuantity%values ) .and. &
          & iand ( ichar(aQuantity%mask(:,:)), m_fill ) == 0 )
      end if
    case ('max', 'min', 'mean', 'median')
      a = statFun( trim(flattening), aQuantity%values, aQuantity%mask )
      b = c
      if ( arg(2) == 'b' ) &
        & b = statFun( trim(flattening), bQuantity%values, aQuantity%mask )
      tvalue = isRelation( trim(op), a, b )
    case default
      call mlsmessage (MLSMSG_Error, moduleName, &
        & 'Formula in compare found unrecognized op: ' // trim(op) )
    end select
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), tvalue, countEmpty=countEmpty )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( switchDetail(switches, 'bool') > 0 ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
    call MLSMessageCalls( 'pop' )
  contains
    elemental logical function isRelation( relation, a, b )
      ! Do inputs a and b stand in relation ('>', '<', '=') ?
      ! Args
      character(len=*), intent(in) :: relation ! ('>', '<', '=')
      real(rv), intent(in) :: a, b
      ! Executable
      select case (trim(relation))
      case ('<')
        isRelation = (a < b )
      case ('>')
        isRelation = (a > b )
      case ('=')
        isRelation = (a == b )
      case default
        ! We should never have reached here
        isRelation = .false.
      end select
    end function isRelation
    
    function statFun ( name, values, mask )
      ! Evaluate statistical function name of values
      ! masking if appropriate
      ! Args
      character(len=*), intent(in)         :: name ! 'min', etc.
      real(rv), dimension(:,:), intent(in) :: values
      character, dimension(:,:), pointer :: MASK
      real(rv) :: statFun
      ! Internal variables
      real(rv), dimension(size(values,1),size(values,2)) :: array
      real(rv) :: FillValue
      ! Executable
      FillValue = DEFAULTUNDEFINEDVALUE
      if ( associated(mask) ) then
        array = DEFAULTUNDEFINEDVALUE
        where ( iand ( ichar(mask(:,:)), m_Fill ) == 0 )
          array(:,:) = values(:,:)
        end where
      else
        array = values
      end if
      select case (trim(name))
      case ('max')
        statFun = mlsmax( array, FillValue=FillValue )
      case ('min')
        statFun = mlsmin( array, FillValue=FillValue )
      case ('mean')
        statFun = mlsmean( array, FillValue=FillValue )
      case ('median')
        statFun = mlsmedian( array, FillValue=FillValue )
      case default
        ! Should not have got here
        statFun = mlsmedian( array, FillValue=FillValue )
      end select
    end function statFun
  end function BooleanFromComparingQtys

  ! ------------------------------------- BooleanFromFormula --
  function BooleanFromFormula ( name, root ) result(size)
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
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
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
    end if
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
  end function BooleanFromFormula

  ! ------------------------- DumpCommand ------------------------
  subroutine DumpCommand ( Root, QuantityTemplatesDB, &
    & VectorTemplates, Vectors, ForwardModelConfigs, HGrids, griddedDataBase, &
    & FileDataBase )

  ! Process a "dump" command

    use AntennaPatterns_m, only: Dump_Antenna_Patterns_Database
    use Calendar, only: Duration_Formatted, Time_T, TK
    use Declaration_table, only: Num_Value
    use Dump_0, only: Dump, rmsFormat
    use Expr_m, only: Expr
    use FilterShapes_m, only: Dump_Filter_Shapes_Database, &
      & Dump_DACS_Filter_Database
    use ForwardModelConfig, only: Dump, ForwardModelConfig_T
    use GriddedData, only: Diff, Dump, GriddedData_T
    use HGridsDatabase, only: Dump, HGRID_T
    use Init_Tables_Module, only: F_AllBooleans, F_AllForwardModels, &
      & f_AllGriddedData, F_AllHGrids, F_AllL2PCs, F_AllLines, &
      & F_AllPFA, F_AllQuantityTemplates, F_AllSignals, F_AllSpectra, &
      & F_AllVectors, F_AllVectorTemplates, F_AllVGrids, F_AntennaPatterns, &
      & F_Boolean, F_Clean, F_CrashBurn, F_Details, F_DACSFilterShapes, &
      & F_FilterShapes, F_ForwardModel, F_GRID, &
      & F_HGrid, F_L2PC, F_Lines, F_Mark, F_Mask, F_MieTables, &
      & F_OPTIONS, F_PfaData, F_PfaFiles, F_PFANum, F_PFAStru, F_PointingGrids, &
      & F_Quantity, F_Signals,  F_Spectroscopy, F_Stop, F_StopWithError, &
      & F_Template, F_Text, F_TGrid, &
      & F_Vector, F_VectorMask, F_VGrid, &
      & S_DIFF, S_DUMP, S_QUANTITY, S_VECTORTEMPLATE
    use L2ParInfo, only: PARALLEL, CLOSEPARALLEL
    use L2PC_m, only: L2PCDatabase, dumpL2PC => Dump
    use Intrinsic, only: PHYQ_Dimensionless
    use MACHINE, only: NEVERCRASH
    use MLSL2Options, only: NORMAL_EXIT_STATUS, RUNTIMEVALUES, STOPWITHERROR
    use MLSCommon, only: MLSFile_T
    ! use MLSFiles, only: DumpMLSFile => Dump
    use MLSMessageModule, only: MLSMessage, MLSMessageCalls, MLSMessageExit, &
      & MLSMSG_CRASH, MLSMSG_ERROR, MLSMSG_INFO
    use MLSSets, only: FindFirst
    use MLSSignals_m, only: Dump, Signals
    use MLSStrings, only: lowerCase
    use MLSStringLists, only: BooleanValue, SWITCHDETAIL
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
    use output_m, only: output, outputNamedValue
    use PFADataBase_m, only: Dump, Dump_PFADataBase, Dump_PFAFileDataBase, &
      & Dump_PFAStructure, PFAData
    use PointingGrid_m, only: Dump_Pointing_Grid_Database
    use QuantityTemplates, only: Dump, QuantityTemplate_T
    use Read_Mie_m, only: Dump_Mie
    use SpectroscopyCatalog_m, only: Catalog, Dump, Dump_Lines_Database, Lines
    use String_Table, only: Get_String
    use Toggles, only: Gen, Switches, Toggle
    use Trace_m, only: Trace_begin, Trace_end
    use Tree, only: Decoration, Node_Id, Nsons, Source_Ref, Sub_rosa, Subtree
    use Tree_Types, only: N_Spec_Args
    use VectorsModule, only: Vector_T, VectorTemplate_T, &
      & Diff, Dump, DumpQuantityMask, DumpVectorMask, & ! for vectors, vector quantities and templates
      & GetVectorQtyByTemplateIndex
    use VGridsDatabase, only: Dump, VGrids

    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Databases:
    type (quantityTemplate_t), dimension(:), pointer, optional   :: QuantityTemplatesDB
    type (forwardModelConfig_t), dimension(:), pointer, optional :: ForwardModelConfigs
    type (vectorTemplate_T), dimension(:), pointer, optional     :: VectorTemplates
    type (vector_T), dimension(:), optional                      :: Vectors
    type (HGrid_T), dimension(:), pointer, optional              :: HGrids
    type (GriddedData_T), dimension(:), pointer, optional        :: griddedDataBase
    type (MLSFile_T), dimension(:), pointer, optional            :: FileDataBase

    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND13_OK'
    logical :: Clean
    logical, parameter :: countEmpty = .true.
    real(tk) :: CPUTime, CPUTimeBase = 0.0_tk
    character(8) :: Date
    integer :: DetailReduction
    integer :: Details
    integer :: DiffOrDump
    integer :: FieldIndex
    integer :: FileIndex
    logical :: GotFirst ! of something -- needed if diffing 2 of them
    logical :: GotOne ! of something -- used to test loop completion
    integer :: GSON, I, J, K, L, Look
    logical :: HaveQuantityTemplatesDB, HaveVectorTemplates, HaveVectors, &
      &        HaveForwardModelConfigs, HaveGriddedData, HaveHGrids
    character(len=80) :: NAMESTRING  ! E.g., 'L2PC-band15-SZASCALARHIRES'
    character(len=80) :: OPTIONSSTRING  ! E.g., '-rbs' (see dump_0.f90)
    integer :: QuantityIndex
    integer :: QuantityIndex2
    integer :: Son
    integer :: Source ! column*256 + line
    character :: TempText*20, Text*255
    type(time_t) :: Time
    character(10) :: TimeOfDay
    logical :: tvalue
    integer :: VectorIndex
    integer :: VectorIndex2
    integer :: Type     ! of the Details expr -- has to be num_value
    integer :: Units(2) ! of the Details expr -- has to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr
    integer :: What

    ! Error codes
    integer, parameter :: Dimless = 1
    integer, parameter :: NoFWM = dimless + 1
    integer, parameter :: noGriddedData = NoFWM + 1
    integer, parameter :: NoHGrid = noGriddedData + 1
    integer, parameter :: NoLines = noHGrid + 1
    integer, parameter :: NoQT = noLines + 1
    integer, parameter :: NoSignals = noQT + 1
    integer, parameter :: NoTG = noSignals + 1
    integer, parameter :: NoVectors = noTG + 1
    integer, parameter :: NoVT = noVectors + 1
    integer, parameter :: Numeric = noVT + 1
    integer, parameter :: Stop = numeric + 1
    integer, parameter :: Unknown = stop + 1 ! Unknown template

    if ( toggle(gen) ) then
      call trace_begin ( 'DumpCommand', root )
    else
      call MLSMessageCalls( 'push', constantName=ModuleName )
    end if
    ! Were we called to do a diff or a dump?
    ! The following must be one of (/ s_dump, s_diff /)
    DiffOrDump = get_spec_id(root)
    if ( .not. any( DiffOrDump == (/ s_dump, s_diff /) ) ) &
    &  call MLSMessage( MLSMSG_Error, moduleName, &
        & "Expected either Diff or Dump command to get here")
    haveQuantityTemplatesDB = present(quantityTemplatesDB)
    if ( haveQuantityTemplatesDB ) &
      & haveQuantityTemplatesDB = associated(quantityTemplatesDB)
    haveVectorTemplates = present(vectorTemplates)
    if ( haveVectorTemplates ) &
      & haveVectorTemplates = associated(vectorTemplates)
    haveVectors = present(vectors)
    if ( haveVectors ) haveVectors = size(Vectors) > 0
    haveForwardModelConfigs = present(forwardModelConfigs)
    if ( haveForwardModelConfigs ) &
      & haveForwardModelConfigs = associated(forwardModelConfigs)
    haveHGrids = present(hGrids)
    if ( haveHGrids ) haveHGrids = associated(hGrids)
    haveGriddedData = present(griddedDataBase)
    if ( haveGriddedData ) haveGriddedData = associated(griddedDataBase)
    
    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    else if ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    end if

    clean = .false.
    GotFirst = .false.
    details = 0 - DetailReduction
    OPTIONSSTRING = '-'

    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case ( f_allBooleans, f_allForwardModels, f_allGriddedData, &
        & f_allHGrids, f_allLines, &
        & f_allPFA, f_allQuantityTemplates, f_allSignals, f_allSpectra, &
        & f_allVectors, f_allVectorTemplates, f_allVGrids, f_antennaPatterns, &
        & f_crashBurn, f_DACSfilterShapes, f_filterShapes, f_MieTables, &
        & f_pfaFiles, f_pfaStru, f_pointingGrids, f_stop, f_stopWithError )
        if ( get_boolean(son) ) then
          select case ( fieldIndex )
          case ( f_allBooleans )
            call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
            & 'Run-time Boolean flags' )
          case ( f_allForwardModels )
            if ( haveForwardModelConfigs ) then
              call dump ( forwardModelConfigs, where=son )
            else
              call announceError ( son, noFWM )
            end if
          case ( f_allGriddedData )
            if ( details < -1 ) cycle
            if ( haveGriddedData ) then
              call dump ( griddedDataBase, details )
            else
              call announceError ( son, noGriddedData )
            end if
          case ( f_allHGrids )
            if ( haveHGrids ) then
              call dump ( hGrids )
            else
              call announceError ( son, noHGrid )
            end if
          case ( f_allL2PCs )
            call dumpL2PC( L2PCDataBase )
          case ( f_allLines )
            if ( associated(lines) ) then
              call dump_lines_database
            else
              call announceError ( son, noLines )
            end if
          case ( f_allPFA )
            if ( details < -1 ) cycle
            call Dump_PFADataBase ( details )
          case ( f_allQuantityTemplates )
            if ( haveQuantityTemplatesDB ) then
              call dump ( quantityTemplatesDB )
            else
              call announceError ( son, noQT )
            end if
          case ( f_allSignals )
            if ( details < -1 ) cycle
            if ( associated(signals) ) then
              call dump ( signals, details=details>0 )
            else
              call announceError ( son, noSignals )
            end if
          case ( f_allSpectra )
            if ( details < -1 ) cycle
            call dump ( catalog, details=details )
          case ( f_allVectors )
            if ( details < -1 ) cycle
            if ( haveVectors ) then
              call dump ( vectors, details=details )
            else
              call announceError ( son, noVectors )
            end if
          case ( f_allVectorTemplates )
            if ( details < -1 ) cycle
            if ( haveVectorTemplates ) then
              call dump ( vectorTemplates, details=details )
            else
              call announceError ( son, noVT )
            end if
          case ( f_allVGrids )
            call dump ( vGrids )
          case ( f_antennaPatterns )
            call dump_antenna_patterns_database ( son )
          case ( f_crashBurn )
              NEVERCRASH = .false.
              call MLSMessage( MLSMSG_CRASH, moduleName, &
                & "Program stopped by /crashBurn field on DUMP statement.")
          case ( f_DACSfilterShapes )
            call dump_dacs_filter_database ( son )
          case ( f_filterShapes )
            call dump_filter_shapes_database ( son )
          case ( f_MieTables )
            call dump_Mie ( details )
          case ( f_pfaFiles )
            call dump_PFAFileDatabase ( details )
          case ( f_pfaStru )
            call dump_PFAStructure ( details )
          case ( f_pointingGrids )
            call dump_pointing_grid_database ( son )
          case ( f_stop )
            if ( NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave ) then
              call MLSMessageExit( NORMAL_EXIT_STATUS, &
                & farewell="Program stopped with normal status by /stop field on DUMP statement.")
            elseif( parallel%slave ) then
              call closeParallel(0)
              call MLSMessageExit( &
                & farewell="slave stopped by /stop field on DUMP statement.")
            else
              call MLSMessageExit( &
                & farewell="Program stopped by /stop field on DUMP statement.")
            endif
          case ( f_stopWithError )
              call MLSMessageExit( 1, &
                & farewell="Program stopped by /stopWithError field on DUMP statement.")
          end select
        end if
      case ( f_Boolean )
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        booleanString = lowerCase(booleanString)
        call output( trim(booleanString) // ' = ', advance='no' )
        tvalue = BooleanValue ( booleanString, &
          & runTimeValues%lkeys, runTimeValues%lvalues)
        call output( tvalue, advance='yes' )
      case ( f_clean )
        clean = get_boolean(son)
      case ( f_details )
        call expr ( gson, units, values, type )
        if ( units(1) /= phyq_dimensionless ) call AnnounceError ( gson, dimless )
        if ( type /= num_value ) call announceError ( gson, numeric )
        details = nint(values(1)) - DetailReduction
      case ( f_forwardModel ) ! Dump forward model configs
        if ( details < -1 ) cycle
        if ( haveForwardModelConfigs ) then
          do i = 2, nsons(son)
            call dump ( & ! has no details switch
              & forwardModelConfigs(decoration(decoration(subtree(i,son)))) )
          end do
        else
          call announceError ( gson, noFWM )
        end if
      case ( f_Grid )    ! Diff or Dump Griddeddata
        if ( details < -1 ) cycle
        if ( haveGriddedData ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            if ( gotFirst ) then
              vectorIndex2 = decoration(decoration(gson))
            else
              vectorIndex = decoration(decoration(gson))
            endif
            if ( DiffOrDump == s_diff ) then
              rmsFormat = '(1pe8.1)'
              if ( gotFirst ) &
                & call diff ( &
                & griddedDataBase(vectorIndex), &
                & griddedDataBase(vectorIndex2), &
                & options=optionsString )
              rmsFormat = '*'
            else
              call output ( ' GriddedData: ' )
              call dump ( &
                & griddedDataBase(vectorIndex), details )
            endif
          end do
        else
          call announceError ( gson, noGriddedData )
        end if
        GotFirst = .true.
      case ( f_hGrid )    ! Dump HGrids
        if ( details < -1 ) cycle
        if ( haveHGrids ) then
          do i = 2, nsons(son)
            call output ( ' HGrid ' )
            call dump ( & ! has no details switch
              & hGrids(decoration(decoration(subtree(i,son)))) )
          end do
        else
          call announceError ( gson, noHGrid )
        end if
      case ( f_L2PC )    ! Dump L2PC
        if ( details < -1 .or. .not. present(FileDataBase) ) cycle
        call get_string( sub_rosa(subtree(2,son)), nameString, strip=.true. )
        fileIndex = FindFirst( FileDataBase%ShortName, trim(nameString) )
        ! call outputNamedValue ( 'name string', trim(nameString) )
        ! call outputNamedValue ( 'file index', fileIndex )
        ! call dumpMLSFile( fileDataBase, details=1 )
        if ( fileIndex < 1 ) cycle
        call output ( ' L2PC short file name: ' // trim(nameString), advance='yes' )
        call dumpL2PC ( &
          & fileDataBase(fileIndex), details )
      case ( f_lines )
        do i = 2, nsons(son)
          what = decoration(decoration(subtree(i,son)))
          call output ( what, after=': ' )
          call dump ( lines(what) )
        end do
      case ( f_mark )
        if ( get_boolean(son) ) call cpu_time ( cpuTimeBase )
      case ( f_mask, f_quantity ) ! Diff or Dump vector quantities
        if ( details < -1 ) cycle
        do i = 2, nsons(son)
          gson = subtree(i,son)
          if ( gotFirst ) then
            vectorIndex2 = decoration(decoration(subtree(1,gson)))
            quantityIndex2 = decoration(decoration(decoration(subtree(2,gson))))
          else
            vectorIndex = decoration(decoration(subtree(1,gson)))
            quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          endif
          if ( DiffOrDump == s_diff ) then
            rmsFormat = '(1pe8.1)'
            if ( gotFirst ) &
              & call diff ( &
              & GetVectorQtyByTemplateIndex( &
              & vectors(vectorIndex), quantityIndex), &
              & GetVectorQtyByTemplateIndex( &
              & vectors(vectorIndex2), quantityIndex2), &
              & options=optionsString )
            rmsFormat = '*'
          elseif ( fieldIndex == f_mask ) then
            call dumpQuantityMask ( GetVectorQtyByTemplateIndex( &
              & vectors(vectorIndex), quantityIndex), details=details )
          else
            if ( clean ) optionsString = trim(optionsString) // 'c'
            call dump ( GetVectorQtyByTemplateIndex( &
              & vectors(vectorIndex), quantityIndex), details=details, &
              & vector=vectors(vectorIndex), options=optionsString )
          end if
        end do
        GotFirst = .true.
      case ( f_options )
        call get_string ( sub_rosa(gson), optionsString, strip=.true. )
        optionsString = lowerCase(optionsString)
        call outputNamedValue( 'options', trim(optionsString) )
      case ( f_pfaData )
        do i = 2, nsons(son)
          look = decoration(decoration(subtree(i,son)))
          call dump ( pfaData(look), details, look )
        end do
      case ( f_pfaNum )
        do i = 2, nsons(son)
          call expr ( subtree(i,son), units, values, type )
          if ( units(1) /= phyq_dimensionless ) call AnnounceError ( subtree(i,son), dimless )
          if ( type /= num_value ) call announceError ( subtree(i,son), numeric )
          call dump ( pfaData(nint(values(1))), details, nint(values(1)) )
        end do
      case ( f_signals )
        do i = 2, nsons(son)
          what = decoration(decoration(subtree(i,son)))
          call output ( what, after=': ' )
          call dump ( signals(what), details=details>0 )
        end do
      case ( f_spectroscopy )
        do i = 2, nsons(son)
          what = decoration(subtree(i,son))
          call output ( what, after=': ' )
          call dump ( catalog(what), details=details )
        end do
      case ( f_template ) ! Dump vector templates or quantity templates
        if ( details < -1 ) cycle
        do i = 2, nsons(son)
          gson = subtree(i,son)
          look = decoration(gson)
          if ( node_id(look) /= n_spec_args ) call announceError ( gson, unknown )
          what = decoration(look)
          select case ( get_spec_id(look) )
          case ( s_quantity )
            if ( haveQuantityTemplatesDB ) then
              call output ( ' Quantity template' )
              call dump ( quantityTemplatesDB(what), details=details )
            else
              call announceError ( gson, noQT )
            end if
          case ( s_vectorTemplate )
            if ( haveVectorTemplates ) then
              call output ( ' Vector template' )
              call dump ( vectorTemplates(what), details=details, quantities=quantityTemplatesDB )
            else
              call announceError ( gson, noVT )
            end if
          end select
        end do
      case ( f_text )
        do k = 2, nsons(son)
          call get_string ( sub_rosa(subtree(k,son)), text, strip=.true. )
          ! Special text causing intentional error (this is a dubious hack)
          ! if ( index(lowercase(text), 'error condition alpha') > 0 ) &
          !   & call MLSMessage ( MLSMSG_Error, moduleName, trim(text) )
          ! Replace format marks: %[Cc] => CPU time in YyDdHH:MM:SS.SSS format
          !                       %[Dd] => date and time of day
          !                       %[Ss] => CPU time in seconds
          !                       %[Tt] => time of day
          gotOne = .false.
          do
            i = max(index(text,'%c'), index(text,'%C'))
            if ( i /= 0 ) then
              tempText = ''
              l = 1 ! Position in TempText
              call cpu_time ( cpuTime )
              time = duration_formatted ( cpuTime - cpuTimeBase )
              if ( time%year /= 0 ) then
                gotOne = .true.
                write ( tempText, * ) time%year
                tempText = adjustl(tempText)
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = 'y'
              end if
              if ( time%day /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), * ) time%day
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = 'd'
              end if
              if ( time%hours /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), '(i2.2)' ) time%hours
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = ':'
              end if
              if ( time%minutes /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), '(i2.2)' ) time%minutes
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = ':'
              end if
              write ( tempText(l:), '(f6.3)' ) time%seconds
              tempText(l:) = adjustl(tempText(l:))
              l = len_trim(tempText)
              text = text(:i-1) // tempText(:l) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%d'), index(text,'%D'))
            if ( i /= 0 ) then
              call date_and_time ( date=date, time=timeOfDay )
              text = text(:i-1) // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // &
                & ' ' // timeOfDay(1:2) // ':' // timeOfDay(3:4) // ':' // timeOfDay(5:10) // &
                & text(i+2:)
              cycle
            end if
            i = max(index(text,'%l'), index(text,'%L'))
            if ( i /= 0 ) then
              source = source_ref(subtree(k,son))
              write ( tempText(1:10), * ) source/256
              write ( tempText(11:20), * ) mod(source,256)
              text = text(:i-1) // "line " // trim(adjustl(tempText(1:10))) // &
                & ", column " // trim(adjustl(tempText(11:20))) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%s'), index(text,'%S'))
            if ( i /= 0 ) then
              write ( tempText, * ) cpuTime - cpuTimeBase
              text = text(:i-1) // trim(adjustl(tempText)) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%t'), index(text,'%T'))
            if ( i /= 0 ) then
              call date_and_time ( time=timeOfDay )
              text = text(:i-1) // timeOfDay(1:2) // ':' // timeOfDay(3:4) // &
                & ':' // timeOfDay(5:10) // text(i+2:)
              cycle
            end if
            exit ! Didn't find a format trigger
          end do
          call output ( trim(text), advance='yes' )
        end do ! k
      case ( f_tGrid )
        if ( details < -1 ) cycle
        do i = 2, nsons(son)
          call output ( ' TGrid ' )
          call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
        end do
      case ( f_vectormask, f_vector ) ! Dump entire vectors
        if ( details < -1 ) cycle
        if ( haveVectors ) then
          do i = 2, nsons(son)
            call output ( ' Vector ' )
            if ( fieldIndex == f_vectormask ) then
              call dumpVectorMask ( vectors(decoration(decoration(subtree(i,son)))), &
                & details=details )
            else
              call dump ( vectors(decoration(decoration(subtree(i,son)))), &
                & details=details, clean=clean )
            end if
          end do
        else
          call announceError ( gson, noVectors )
        end if
      case ( f_vGrid )
        if ( details < -1 ) cycle
        do i = 2, nsons(son)
          call output ( ' VGrid ' )
          call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
        end do
      end select
    end do

    if ( toggle(gen) ) then
      call trace_end ( 'DumpCommand' )
    else
      call MLSMessageCalls( 'pop' )
    end if

  contains

    subroutine AnnounceError ( where, what )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: NewLine

      integer, intent(in) :: What, Where

      call StartErrorMessage ( where )

      select case ( what )
      case ( dimless )
        call output ( "The field is not unitless." )
      case ( noFWM )
        call output ( "Can't dump Forward Model Configs here." )
      case ( noGriddedData )
        call output ( "Can't dump GriddedData here." )
      case ( noHGrid )
        call output ( "Can't dump HGrids here." )
      case ( noLines )
        call output ( "Can't dump Lines here." )
      case ( noQT )
        call output ( "Can't dump Quantity Templates here." )
      case ( noSignals )
        call output ( "Can't dump Signals here." )
      case ( noTG )
        call output ( "Can't dump TGrids here." )
      case ( noVectors )
        call output ( "Can't dump Vectors here." )
      case ( noVT )
        call output ( "Can't dump Vector Templates here." )
      case ( numeric )
        call output ( "The field is not numeric." )
      case ( stop )
        call output ( "Program stopped by /stop field on DUMP statement." )
      case ( unknown )
        call output ( "Can't figure out what kind of template it is." )
      end select
      call newLine
    end subroutine AnnounceError

  end subroutine DumpCommand
  
  logical function Skip ( Root )
    ! Returns value of Boolean field;
    ! If TRUE should skip rest of section in which SKIP command appears
    use Init_Tables_Module, only: F_Boolean
    use MLSL2Options, only: runTimeValues
    use MLSMessageModule, only: MLSMessageCalls
    use MLSStringLists, only: BooleanValue
    use MLSStrings, only: lowerCase
    use MoreTree, only: Get_Field_ID
    use Output_M, only: Output
    use String_Table, only: Get_String
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_begin, Trace_end
    use Tree, only: Nsons, Sub_rosa, Subtree
    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Internal variables
    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND13_OK'
    integer :: GSON, J
    integer :: FieldIndex
    integer :: Son
    ! Executable
    if ( toggle(gen) ) then
      call trace_begin ( 'DumpCommand', root )
    else
      call MLSMessageCalls( 'push', constantName=ModuleName )
    end if
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        booleanString = lowerCase(booleanString)
        call output( trim(booleanString) // ' = ', advance='no' )
        skip = BooleanValue ( booleanString, &
          & runTimeValues%lkeys, runTimeValues%lvalues)
        call output( skip, advance='yes' )
        if ( skip ) &
          & call output( '(Skipping rest of this section)', advance='yes' )
      case default
        ! Should not have got here if parser worked correctly
      end select
      if ( toggle(gen) ) then
        call trace_end ( 'DumpCommand' )
      else
        call MLSMessageCalls( 'pop' )
      end if
    enddo
  end function Skip

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

end module DumpCommand_M

! $Log$
! Revision 2.48  2009/11/02 21:22:25  pwagner
! May diff grided data
!
! Revision 2.47  2009/10/27 22:18:18  pwagner
! Implemented new Diff command; so far only of vector quantities
!
! Revision 2.46  2009/09/15 20:03:48  pwagner
! Dump commands take boolean fields /stop, /stopWithError, /crashBurn
!
! Revision 2.45  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.44  2008/12/18 21:12:00  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.43  2008/09/30 22:00:22  vsnyder
! Add PRINT statement to Not_Used_Here to reduce compilation cascades
!
! Revision 2.42  2008/06/05 02:07:43  vsnyder
! Dump Mie tables
!
! Revision 2.41  2007/12/07 01:12:43  pwagner
! Lets us catch warnings and assign to runtime Booleans
!
! Revision 2.40  2007/11/15 22:54:51  pwagner
! Functions to set runtimeBooleans moved here
!
! Revision 2.39  2007/11/05 18:36:44  pwagner
! May Skip remaining lines in Fill, Join, Retrieve sections depending on Boolean
!
! Revision 2.38  2007/10/09 00:32:05  pwagner
! Added ability to dump masks of quantities, vectors
!
! Revision 2.37  2007/10/03 23:53:05  vsnyder
! Stop with error if /stop and switch 'erh' is set
!
! Revision 2.36  2007/08/17 00:32:16  pwagner
! Unneeded changes
!
! Revision 2.35  2007/04/03 17:37:17  vsnyder
! Check Vectors for zero size instead of associated
!
! Revision 2.34  2007/01/11 20:44:55  vsnyder
! Add Tracing
!
! Revision 2.33  2006/09/21 18:48:33  pwagner
! Reduce level of dumps in SIDS version
!
! Revision 2.32  2006/07/27 03:52:41  vsnyder
! Pass details to dumps for vectors and vector templates
!
! Revision 2.31  2006/07/19 22:26:40  vsnyder
! Comment out unused USE
!
! Revision 2.30  2006/06/12 16:28:25  pwagner
! Added ability to dump Gridded Data
!
! Revision 2.29  2006/05/31 22:38:17  vsnyder
! Revert to details= semantics that Van prefers
!
! Revision 2.28  2006/05/03 20:14:05  pwagner
! details and /stop properly implemented
!
! Revision 2.27  2006/03/22 02:19:48  vsnyder
! Add Vector argument to quantity dump to get its name
!
! Revision 2.26  2006/03/07 16:23:52  pwagner
! Fixed bug only NAG caught
!
! Revision 2.25  2006/03/07 00:49:42  pwagner
! May dump Booleans
!
! Revision 2.24  2006/01/04 01:21:30  vsnyder
! Make optional arguments pointers, in case actuals aren't associated
!
! Revision 2.23  2005/06/03 02:06:55  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of an argument, add dumps for
! PFA structure, PFA datum by number.
!
! Revision 2.22  2005/05/02 23:11:37  vsnyder
! Add dump of PFAFiles database
!
! Revision 2.21  2005/04/01 20:48:28  vsnyder
! Add mark and text fields to dump command
!
! Revision 2.20  2005/03/26 01:34:00  vsnyder
! Add stop message
!
! Revision 2.19  2005/03/15 01:36:08  vsnyder
! Add newline after error messages
!
! Revision 2.18  2005/01/12 03:18:51  vsnyder
! Add item number to PFA dump
!
! Revision 2.17  2004/12/28 00:22:03  vsnyder
! Add not_used_here
!
! Revision 2.16  2004/12/13 20:13:04  vsnyder
! Add dumps for AllLines, AllSignals, AllSpectra, Lines, Signals, Spectroscopy,
! and a Stop command.
!
! Revision 2.15  2004/11/04 06:37:34  vsnyder
! Index spetroscopy catalog by molecule instead of searching
!
! Revision 2.14  2004/11/01 20:16:20  vsnyder
! Check for spectroscopy catalog before trying to dump it
!
! Revision 2.13  2004/10/30 00:26:46  vsnyder
! Add 'spectroscopy' field to DumpCommand
!
! Revision 2.12  2004/10/06 20:19:39  vsnyder
! Cannonball polishing
!
! Revision 2.11  2004/09/24 22:24:20  vsnyder
! Make PFA dump aware of 'details' switch
!
! Revision 2.10  2004/07/17 02:28:19  vsnyder
! Add dump for entire PFA database
!
! Revision 2.9  2004/06/12 00:41:30  vsnyder
! Allow all fields except details to be arrays
!
! Revision 2.8  2004/06/09 19:59:38  pwagner
! Gets PFAData type and dump method from PFADataBase_m
!
! Revision 2.7  2004/06/08 20:20:18  vsnyder
! Add tGrid
!
! Revision 2.6  2004/05/29 02:50:49  vsnyder
! Added more dumps
!
! Revision 2.5  2004/05/22 02:31:23  vsnyder
! Dump PFAData, VGrids
!
! Revision 2.4  2004/05/20 19:47:36  vsnyder
! Move Dump*Hgrid from Dumper to HgridsDatabse
!
! Revision 2.3  2004/05/18 01:18:51  vsnyder
! Add dump for HGrid
!
! Revision 2.2  2004/05/11 02:53:29  vsnyder
! Remove USEs for unreferenced symbols
!
! Revision 2.1  2004/05/01 04:04:16  vsnyder
! Initial commit
!
