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
! (Should these latter functions be moved into a special Boolean module?)

  implicit none
  private

  public :: BOOLEANFROMANYGOODRADIANCES
  public :: BOOLEANFROMANYGOODVALUES
  public :: BOOLEANFROMCATCHWARNING
  public :: BOOLEANFROMCOMPARINGQTYS
  public :: BOOLEANFROMEMPTYGRID, BOOLEANFROMEMPTYSWATH
  public :: BOOLEANFROMFORMULA
  public :: MLSCASE, DUMPCOMMAND, MLSENDSELECT, MLSSELECT, SKIP

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!
!     (Data)
! MLSSelecting  true if in a Select .. Case .. EndSelect control structure
!               and the present Case doesn't match
!     (subroutines and functions)
! The following functions do much the same thing: They add to the
! runtime Boolean datase a pair (key => named_Boolean, value => its_value)
! where the value is determined by whether condition is true or false
!    function                 condition   
! booleanFromAnyGoodRadiances  
!               Any of the radiances have non-negative precisions
! booleanFromAnyGoodValues  
!               The quantity has any useable values (based on precision, status)
! booleanFromCatchWarning 
!               The last command resulted in a (optionally specific) warning
! booleanFromComparingQuantities 
!               The first quantity stands in specified relation to the second
!               E.g., formula="all a > b"
! booleanFromEmptyGrid
!               The specified GriddedData is empty
! booleanFromEmptySwath
!               The specified swath in the specified file has no useable data
! booleanFromFormula
!               (a) Evaluate the formula; it may be one of two forms
!               (1) Contains only "or", "and", "not" => logical result
!               (2) Contains "lhs == rhs" or "lhs /= rhs" => Does it?
!               (b) Just sture the text of the label field
!
! The following subroutines depart from the sbove pattern
! DumpCommand    Process the Dump command, dumping any of the allowed datatypes
! MLSCase        Process the Case control statement
! MLSSelect      Process the Select control statement
! MLSEndSelect   Process the EndSelect control statement
! Skip           Process the Skip control statement
! === (end of toc) ===


  logical, parameter :: countEmpty = .true. ! Except where overriden locally
  integer, parameter :: MAXRESULTLEN = 64
  logical, public, save             :: MLSSELECTING = .false.
  logical, save               :: MLSSelectedAlready = .false.
  character(LEN=MAXRESULTLEN), save :: selectLabel  = ' '
contains

  ! ------------------------------------- BooleanFromAnyGoodRadiances --
  function BooleanFromAnyGoodRadiances ( ROOT, CHUNK, FILEDATABASE ) &
    & result( HASHSIZE )
    use ALLOCATE_DEALLOCATE, only: DEALLOCATE_TEST
    use CONSTRUCTQUANTITYTEMPLATES, only: ANYGOODSIGNALDATA
    use CHUNKS_M, only: MLSCHUNK_T
    use DUMP_0, only: DUMP
    use INIT_TABLES_MODULE, only: F_SIGNAL, F_BOOLEAN
    use MLSCOMMON, only: MLSFILE_T
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSIGNALS_M, only: GETSIGNALNAME, &
      & SIGNALS
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT
    use PARSE_SIGNAL_M, only: PARSE_SIGNAL
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (MLSChunk_T), intent(in) :: chunk
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer             :: hashsize
    ! Internal variables
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
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
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
      if ( verbose ) &
        & call output( 'signal: ' // trim(signalString), advance='yes' )
      call Parse_signal(signalString, signal_indices)
      tvalue = .false.
      ! Loop over signals, or-ing them until we get TRUE
      do s=1, size(signal_indices)
        signalIndex = signal_indices(s)
        if ( verbose ) then
          call GetSignalName ( signalIndex, subSignalString, &                   
            & sideband=signals(signalIndex)%sideband, noChannels=.TRUE. )
          call output( 'sub-signal: ' // trim(subSignalString), advance='yes' )
        end if
        tvalue = tvalue .or. &
          & AnyGoodSignalData ( signalIndex, signals(signalIndex)%sideband, &
          & filedatabase, chunk )
        if ( tvalue ) then
          if ( verbose ) then
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
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    hashsize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodRadiances

  ! ------------------------------------- BooleanFromAnyGoodValues --
  function BooleanFromAnyGoodValues ( ROOT, VECTORS ) result( THESIZE )
    use DUMP_0, only: DUMP
    use INIT_TABLES_MODULE, only: F_PRECISION, F_QUALITY, &
      & F_QUANTITY, F_BOOLEAN, F_STATUS
    use MANIPULATEVECTORQUANTITIES, only: ANYGOODDATAINQTY
    use MLSCOMMON, only: RV
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T, &
      & GETVECTORQTYBYTEMPLATEINDEX
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (vector_T), dimension(:) :: Vectors
    integer             :: thesize
    ! Internal variables
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
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
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
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromAnyGoodValues

  ! ------------------------------------- BooleanFromCatchWarning --
  function BooleanFromCatchWarning ( ROOT ) result( SIZE )
    ! Called to check if the last command resulted in a warning
    ! (either printed or suppressed)
    ! and optionally if the warning matches a supplied message string
    ! syntax: 
    ! CatchWarning, [message='string'], Boolean="name"
    use DUMP_0, only: DUMP
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_MESSAGE
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSMESSAGEMODULE, only: MLSMESSAGEINQUIRE
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: root
    integer             :: size
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: LastWarningMsg
    character(len=255) :: message
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
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
    if ( verbose ) then
      call outputNamedValue( 'message to match', message )
      call outputNamedValue( 'LastWarningMsg', LastWarningMsg )
    endif
    if ( len_trim(LastWarningMsg) < 1 ) then
      tvalue = .false.
    else if ( len_trim(message) < 1 ) then
      tvalue = .true.
    else
      ! tvalue = streq( message, LastWarningMsg, '-wcf' )
      ! This allows partial matches, case-insensitive
      tvalue = &
        & index( lowerCase(LastWarningMsg), lowerCase(trim(adjustl(message))) ) > 0
    end if
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromCatchWarning

  ! ------------------------------------- BooleanFromComparingQtys --
  function BooleanFromComparingQtys ( ROOT, VECTORS ) result( THESIZE )
    use DUMP_0, only: DUMP
    use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_A, F_B, F_C, F_BOOLEAN, F_FORMULA
    use MLSCOMMON, only: R8, RV, DEFAULTUNDEFINEDVALUE
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMESSAGECALLS, MLSMSG_ERROR
    use MLSSTATS1, only: MLSMAX, MLSMIN, MLSMEAN, MLSMEDIAN
    use MLSSTRINGLISTS, only: GETSTRINGELEMENT, NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & REPLACESUBSTRING, SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T, M_FILL, &
      & GETVECTORQTYBYTEMPLATEINDEX
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
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
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
      call MLSMessage (MLSMSG_Error, moduleName, &
        & 'Formula in compare must be "flattening a op [b][c]".' )
    else if( index('<>=', trim(op)) < 1 ) then
      call MLSMessage (MLSMSG_Error, moduleName, &
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
      call MLSMessage (MLSMSG_Error, moduleName, &
        & 'Formula in compare found unrecognized op: ' // trim(op) )
    end select
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
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

  ! ------------------------------------- BooleanFromEmptyGrid --
  function BooleanFromEmptyGrid ( ROOT, GRIDS ) result( THESIZE )
    use DUMP_0, only: DUMP
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_GRID
    use GRIDDEDDATA, only: GRIDDEDDATA_T
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: root
    type (GRIDDEDDATA_T), dimension(:), pointer :: Grids
    integer             :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    type (GRIDDEDDATA_T), pointer :: grid
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    integer :: value
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      value = subtree(2,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_grid ) 
        grid => grids ( decoration ( decoration ( value ) ) )
      end select
    end do
    tvalue = .true.
    if ( associated(grid) ) tvalue = grid%empty
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromEmptyGrid

  ! ------------------------------------- BooleanFromEmptySwath --
  ! Returns TRUE if there are no useable data points in the swath
  ! (or if the swath is not in the file at all)
  ! The useablility criterion is
  ! (1) Precision non-negative; and
  ! (2) Status even
  ! If even one point is useable (a very low bar, admittedly) then
  ! return FALSE
  function BooleanFromEmptySwath ( ROOT ) result( THESIZE )
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use DUMP_0, only: DUMP
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_FILE, F_SWATH, F_TYPE, &
      & L_L2DGG, L_L2GP
    use L2GPDATA, only: L2GPDATA_T, RGP, L2GPNAMELEN, &
      & READL2GPDATA, DESTROYL2GPCONTENTS
    use MLSCOMMON, only: FILENAMELEN
    use MLSFILES, only: HDFVERSION_5
    use MLSHDFEOS, only: MLS_SWATH_IN_FILE
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_WARNING
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSPCF2, only: MLSPCF_L2GP_END, &
      & MLSPCF_L2GP_START, MLSPCF_L2DGG_START, MLSPCF_L2DGG_END
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, &
      & SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: root
    integer             :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    integer :: fileType
    character (len=FileNameLen) :: FILE_BASE    ! From the FILE= field
    character(len=FileNameLen) :: filename          ! filename
    character(len=L2GPNAMELEN) :: swathname         ! swathname
    integer :: i
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    integer :: value
    type (L2GPData_T) :: l2gp
    integer :: numGood
    logical, dimension(:), pointer  :: negativePrec => null() ! true if all prec < 0
    logical, dimension(:), pointer  :: oddStatus => null() ! true if all status odd
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      value = subtree(2,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_file )
        call get_string ( sub_rosa(subtree(2,son)), file_base, strip=.true. )
      case ( f_swath )
        call get_string ( sub_rosa(subtree(2,son)), swathname, strip=.true. )
      case ( f_type )
        filetype = decoration(subtree(2, son))
      end select
    end do
    select case (fileType)
    case ( l_l2dgg )
      call returnFullFileName( file_base, Filename, &
        & mlspcf_l2dgg_start, mlspcf_l2dgg_end )
    case ( l_l2gp )
      call returnFullFileName( file_base, Filename, &
        & mlspcf_l2gp_start, mlspcf_l2gp_end )
    case default
      call MLSMessage( MLSMSG_Warning, ModuleName, &
      & "type should have been either l2gp or dgg; assume you meant dgg" )
    end select
    tvalue = .true.
    if ( mls_swath_in_file( filename, swathname, HDFVERSION_5 ) ) then
      call ReadL2GPData ( trim(filename), trim(swathname), l2gp, &
        & hdfVersion=HDFVERSION_5 )
      call allocate_test( negativePrec, l2gp%nTimes, 'negativePrec', ModuleName )
      call allocate_test( oddStatus, l2gp%nTimes, 'oddStatus', ModuleName )
      do i=1, l2gp%nTimes
        negativePrec(i) = all( l2gp%l2GPPrecision(:,:,i) < 0._rgp )
      enddo
      do i=1, l2gp%nTimes
        oddStatus(i) = mod(l2gp%status(i), 2) > 0
      enddo
      numGood = count( .not. ( negativePrec .or. &
        & (mod(l2gp%status, 2) > 0) ) )
      tvalue = ( numGood < 1 )
      if ( verbose ) then
        call dump( negativePrec, 'num surfs with prec < 0' )
        call dump( oddStatus, 'num surfs with odd status' )
        call outputNamedValue ( 'empty swath?', tvalue )
      endif
      call deallocate_test( negativePrec, 'negativePrec', ModuleName )
      call deallocate_test( oddStatus, 'oddStatus', ModuleName )
      call DestroyL2GPContents ( l2gp )
    endif
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    endif
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromEmptySwath

  ! ------------------------------------- BooleanFromFormula --
  function BooleanFromFormula ( NAME, ROOT ) result( SIZE )
    ! Called either when a Boolean is first declared
    ! syntax: 
    ! name: Boolean, formula="formula"
    !
    ! or when it is reevaluated
    ! syntax: 
    ! Reevaluate, formula="formula", Boolean="name"
    use DUMP_0, only: DUMP
    use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_FORMULA, F_LABEL, F_VALUES
    use MLSKINDS, only: R8
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, PUTHASHELEMENT, SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: SWITCHES
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    ! Dummy args
    integer, intent(in) :: name
    integer, intent(in) :: root
    integer             :: size
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: formula
    integer :: keyNo
    logical :: literal
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    integer, dimension(2) :: UNITASARRAY ! From expr
    real(r8), dimension(2) :: VALUEASARRAY ! From expr
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    literal= .false.
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
        tvalue = myBooleanValue (formula)
      case ( f_label )
        call get_string ( sub_rosa(subtree(2,son)), formula, strip=.true. )
        literal = .true.
      case ( f_values )
        call expr ( son , unitAsArray, valueAsArray )
        tvalue = ( valueAsArray(1) /= 0 )
        ! badRange = valueAsArray
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    if ( literal ) then
      ! print *, 'Oops-you dummy! code this missing piece'
      if ( verbose ) then
        call output( trim(nameString) // ' = ', advance='no' )
        call output( trim(formula), advance='yes' )
      endif
      call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
        & lowercase(trim(nameString)), lowercase(trim(formula)), countEmpty=countEmpty )
    else
      if ( verbose ) then
        call output( trim(nameString) // ' = ', advance='no' )
        call output( tvalue, advance='yes' )
      endif
      call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty )
    endif
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty )
    if ( verboser ) &
      & call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'Run-time Boolean flags' )
  end function BooleanFromFormula

  ! ------------------------- DumpCommand ------------------------
  subroutine DumpCommand ( ROOT, QUANTITYTEMPLATESDB, &
    & VECTORTEMPLATES, VECTORS, FORWARDMODELCONFIGS, HGRIDS, GRIDDEDDATABASE, &
    & FILEDATABASE, MATRIXDATABASE, HESSIANDATABASE )

  ! Process a "dump" command

    use ANTENNAPATTERNS_M, only: DUMP_ANTENNA_PATTERNS_DATABASE
    use CALENDAR, only: DURATION_FORMATTED, TIME_T, TK
    use DECLARATION_TABLE, only: NUM_VALUE
    use DUMP_0, only: DIFF, DUMP, RMSFORMAT
    use EXPR_M, only: EXPR
    use FILTERSHAPES_M, only: DUMP_FILTER_SHAPES_DATABASE, &
      & DUMP_DACS_FILTER_DATABASE
    use FORWARDMODELCONFIG, only: DUMP, FORWARDMODELCONFIG_T
    use GRIDDEDDATA, only: DIFF, DUMP, GRIDDEDDATA_T
    use HESSIANMODULE_1, only: HESSIAN_T, DIFF, DUMP
    use HGRIDSDATABASE, only: DUMP, HGRID_T
    use IGRF_INT, only: DUMP_GH
    use INIT_TABLES_MODULE, only: F_ALLBOOLEANS, F_ALLFILES, &
      & F_ALLFORWARDMODELS, F_ALLGRIDDEDDATA, F_ALLHESSIANS, F_ALLHGRIDS, &
      & F_ALLL2PCS, F_ALLLINES, F_ALLMATRICES, F_ALLPFA, &
      & F_ALLQUANTITYTEMPLATES, F_ALLRADIOMETERS, F_ALLSIGNALS, F_ALLSPECTRA, &
      & F_ALLVECTORS, F_ALLVECTORTEMPLATES, F_ALLVGRIDS, F_ANTENNAPATTERNS, &
      & F_BOOLEAN, &
      & F_CALLSTACK, F_CHUNKNUMBER, F_CLEAN, F_COMMANDLINE, F_CRASHBURN, &
      & F_DETAILS, F_DACSFILTERSHAPES, &
      & F_FILE, F_FILTERSHAPES, F_FORWARDMODEL, F_GRID, F_HEIGHT, F_HESSIAN, &
      & F_HGRID, F_IGRF, F_L2PC, F_LINES, F_MARK, F_MASK, F_MATRIX, &
      & F_MIETABLES, F_OPTIONS, F_PFADATA, F_PFAFILES, F_PFANUM, F_PFASTRU, &
      & F_PHASENAME, F_POINTINGGRIDS, F_QUANTITY, &
      & F_SIGNALS,  F_SPECTROSCOPY, F_STOP, &
      & F_STOPWITHERROR, F_SURFACE, F_TEMPLATE, F_TEXT, F_TGRID, &
      & F_VECTOR, F_VECTORMASK, F_VGRID, &
      & S_DIFF, S_DUMP, S_QUANTITY, S_VECTORTEMPLATE, &
      & FIELD_FIRST, FIELD_LAST
    use L2PARINFO, only: PARALLEL, CLOSEPARALLEL
    use L2PC_M, only: L2PCDATABASE, DUMPL2PC => DUMP
    use INTRINSIC, only: PHYQ_DIMENSIONLESS
    use MACHINE, only: NEVERCRASH
    use MATRIXMODULE_1, only: MATRIX_T, MATRIX_DATABASE_T, &
      & DIFF, DUMP, GETFROMMATRIXDATABASE
    use MLSCOMMON, only: MLSFILE_T
    use MLSFILES, only: DUMPMLSFILE => DUMP, GETMLSFILEBYNAME
    use MLSKINDS, only: RV
    use MLSL2OPTIONS, only: COMMAND_LINE, L2CFNODE, &
      & NORMAL_EXIT_STATUS, RUNTIMEVALUES, &
      & MLSMESSAGE
    use MLSL2TIMINGS, only: CURRENTCHUNKNUMBER, CURRENTPHASENAME, &
      & DUMP_SECTION_TIMINGS
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS, MLSMESSAGEEXIT, &
      & MLSMSG_CRASH, MLSMSG_ERROR
    use MLSSETS, only: FINDFIRST
    use MLSSIGNALS_M, only: DUMP, GETRADIOMETERINDEX, RADIOMETERS, SIGNALS
    use MLSSTRINGS, only: INDEXES, LOWERCASE
    use MLSSTRINGLISTS, only: GETHASHELEMENT, SWITCHDETAIL
    use MORETREE, only: GET_BOOLEAN, GET_FIELD_ID, GET_SPEC_ID
    use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE
    use PFADATABASE_M, only: DUMP, DUMP_PFADATABASE, DUMP_PFAFILEDATABASE, &
      & DUMP_PFASTRUCTURE, PFADATA
    use POINTINGGRID_M, only: DUMP_POINTING_GRID_DATABASE
    use QUANTITYTEMPLATES, only: DUMP, QUANTITYTEMPLATE_T
    use READ_MIE_M, only: DUMP_MIE
    use SPECTROSCOPYCATALOG_M, only: CATALOG, DUMP, DUMP_LINES_DATABASE, LINES
    use STRING_TABLE, only: GET_STRING
    use TIME_M, only: FINISH
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_SPEC_ARGS
    use VECTORSMODULE, only: VECTOR_T, VECTORTEMPLATE_T, VECTORVALUE_T, &
      & DIFF, DUMP, DUMPQUANTITYMASK, DUMPVECTORMASK, & ! FOR VECTORS, VECTOR QUANTITIES AND TEMPLATES
      & GETVECTORQTYBYTEMPLATEINDEX
    use VGRIDSDATABASE, only: DUMP, VGRIDS

    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Databases:
    type (quantityTemplate_t), dimension(:), pointer, optional   :: QuantityTemplatesDB
    type (forwardModelConfig_t), dimension(:), pointer, optional :: ForwardModelConfigs
    type (vectorTemplate_T), dimension(:), pointer, optional     :: VectorTemplates
    type (vector_T), dimension(:), optional                      :: Vectors
    type (HGrid_T), dimension(:), pointer, optional              :: HGrids
    type (GriddedData_T), dimension(:), pointer, optional        :: griddedDataBase
    type (MLSFile_T), dimension(:), pointer, optional            :: FileDataBase
    type (Matrix_Database_T), dimension(:), pointer, optional    :: MatrixDataBase
    type (Hessian_T), dimension(:), pointer, optional            :: HessianDataBase

    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND13_OK'
    logical :: Clean
    real(tk) :: CPUTime, CPUTimeBase = 0.0_tk
    character(8) :: Date
    integer :: DetailReduction
    integer :: Details
    integer :: DiffOrDump
    character(len=128) :: Farewell
    integer :: FieldIndex
    integer :: FileIndex
    logical, dimension(field_first:field_last) :: GOT
    logical :: GotFirst ! of something -- needed if diffing 2 of them
    logical :: GotOne ! of something -- used to test loop completion
    integer :: GSON, I, J, K, L, Look
    logical :: HaveQuantityTemplatesDB, HaveVectorTemplates, HaveVectors, &
             & HaveForwardModelConfigs, HaveGriddedData, HaveHGrids, &
             & HaveMatrices, HaveHessians
    real(rv) :: height  ! We will use this to dump just one surface
    integer :: hessianIndex
    integer :: hessianIndex2
    character(len=80) :: Label  ! E.g., 'BAND8'
    type (Matrix_T), pointer :: matrix
    type (Matrix_T), pointer :: matrix2
    integer :: MatrixIndex
    integer :: MatrixIndex2
    character(len=80) :: NAMESTRING  ! E.g., 'L2PC-band15-SZASCALARHIRES'
    type (MLSFile_T), pointer :: OneMLSFile
    character(len=80) :: OPTIONSSTRING  ! E.g., '-rbs' (see dump_0.f90)
    integer :: QuantityIndex
    integer :: QuantityIndex2
    type (VectorValue_T), pointer :: QTY1, QTY2
    integer :: Son
    integer :: Source ! column*256 + line
    character :: TempText*20, Text*255
    type(time_t) :: Time
    character(10) :: TimeOfDay
    integer :: Type     ! of the Details expr -- has to be num_value
    integer :: VectorIndex
    integer :: VectorIndex2
    integer :: Units(2) ! of the Details expr -- has to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr
    integer :: What

    ! Error codes
    integer, parameter :: Dimless = 1
    integer, parameter :: NoFile = dimless + 1
    integer, parameter :: NoFileDatabase = noFile + 1
    integer, parameter :: NoFWM = noFileDatabase + 1
    integer, parameter :: noGriddedData = NoFWM + 1
    integer, parameter :: NoHGrid = noGriddedData + 1
    integer, parameter :: NoLines = noHGrid + 1
    integer, parameter :: NoQT = noLines + 1
    integer, parameter :: Noradiometers = noQT + 1
    integer, parameter :: NoSignals = noradiometers + 1
    integer, parameter :: NoTG = noSignals + 1
    integer, parameter :: NoVectors = noTG + 1
    integer, parameter :: NoVT = noVectors + 1
    integer, parameter :: Numeric = noVT + 1
    integer, parameter :: Stop = numeric + 1
    integer, parameter :: Unknown = stop + 1 ! Unknown template
    ! Executable
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
    haveMatrices = present(MatrixDatabase)
    if ( haveMatrices ) haveMatrices = associated(MatrixDatabase)
    if ( haveMatrices ) haveMatrices = size(MatrixDatabase) > 0
    HaveHessians = present(HessianDatabase)
    if ( HaveHessians ) then
      if ( associated(HessianDatabase) ) HaveHessians = size(HessianDatabase) > 0
    endif
    
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
    got= .false.

    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      L2CFNODE = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      source = source_ref(gson) ! column + 256*line in l2cf
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      ! This first heaped set of fields need no "right-hand side"
      case ( f_allBooleans, f_allFiles, f_allForwardModels, f_allGriddedData, &
        & f_allHessians, f_allHGrids, f_allL2PCs, f_allLines, f_allMatrices, &
        & f_allPFA, f_allQuantityTemplates, &
        & f_allRadiometers, f_allSignals, f_allSpectra, &
        & f_allVectors, f_allVectorTemplates, f_allVGrids, f_antennaPatterns, &
        & f_callStack, f_chunkNumber, f_commandLine, f_crashBurn, &
        & f_DACSfilterShapes, f_filterShapes, f_igrf, &
        & f_MieTables, f_pfaFiles, f_pfaStru, f_phaseName, f_pointingGrids, &
        & f_stop, f_stopWithError )
        if ( get_boolean(son) ) then
          select case ( fieldIndex )
          case ( f_allBooleans )
            call dump( countEmpty, runTimeValues%lkeys, runTimeValues%lvalues, &
            & 'Run-time Boolean flags' )
          case ( f_allFiles )
            if ( present(fileDataBase) ) then
              call dumpMLSFile ( fileDataBase )
            else
              call announceError ( son, noFileDatabase )
            end if
          case ( f_allForwardModels )
            if ( haveForwardModelConfigs ) then
              call dump ( forwardModelConfigs, where=son, &
                & quantityTemplatesDB=quantityTemplatesDB )
            else
              call announceError ( son, noFWM )
            end if
          case ( f_allGriddedData )
            if ( details < -1 ) cycle
            if ( haveGriddedData ) then
              call dump ( griddedDataBase, details , options=optionsString )
            else
              call announceError ( son, noGriddedData )
            end if
          case ( f_allHessians )
            if ( details < -1 ) cycle
            if ( haveHessians ) then
              call dump ( HessianDataBase, details, options=optionsString )
            else
              call announceError ( son, 0, 'Unable to dump HessianDB here; empty or absent' )
            end if
          case ( f_allHGrids )
            if ( haveHGrids ) then
              call dump ( hGrids )
            else
              call announceError ( son, noHGrid )
            end if
          case ( f_allL2PCs )
            call dumpL2PC( L2PCDataBase, details, options=optionsString )
          case ( f_allLines )
            if ( associated(lines) ) then
              call dump_lines_database
            else
              call announceError ( son, noLines )
            end if
          case ( f_allMatrices )
            if ( details < -1 ) cycle
            if ( haveMatrices ) then
              call dump ( MatrixDataBase, details )
            else
              call announceError ( son, 0, 'Unable to dump MatrixDB here; empty or absent' )
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
          case ( f_allRadiometers )
            if ( details < -1 ) cycle
            if ( associated(Radiometers) ) then
              call dump ( Radiometers )
              if ( details > 1 ) then
                call output('Checking getRadiometerIndex for R1A:118', advance='yes' )
                call getRadiometerIndex('R1A:118', i )
                if ( i > 0 ) call dump( radiometers(i) )
                call output('Checking getRadiometerIndex for R1B:118', advance='yes' )
                call getRadiometerIndex('R1B:118', i )
                if ( i > 0 ) call dump( radiometers(i) )
                call output('Checking getRadiometerIndex for R2:190', advance='yes' )
                call getRadiometerIndex('R2:190', i )
                if ( i > 0 ) call dump( radiometers(i) )
                call output('Checking getRadiometerIndex for R3:240', advance='yes' )
                call getRadiometerIndex('R3:240', i )
                if ( i > 0 ) call dump( radiometers(i) )
                call output('Checking getRadiometerIndex for R3', advance='yes' )
                call getRadiometerIndex('R3', i )
                if ( i > 0 ) call dump( radiometers(i) )
              endif
            else
              call announceError ( son, noRadiometers )
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
          case ( f_callStack )
            call MLSMessageCalls ( 'dump' )
          case ( f_chunkNumber )
            call outputNamedValue ( 'chunk number', currentChunkNumber )
          case ( f_commandLine )
            call outputNamedValue ( 'command line', command_line )
          case ( f_crashBurn )
            call finish ( 'ending mlsl2' )
            NEVERCRASH = .false.
            call MLSMessage( MLSMSG_CRASH, moduleName, &
              & "Program stopped by /crashBurn field on DUMP statement.")
          case ( f_DACSfilterShapes )
            call dump_dacs_filter_database ( son )
          case ( f_filterShapes )
            call dump_filter_shapes_database ( son )
          case ( f_igrf )
            call dump_gh ( details, optionsString )
          case ( f_MieTables )
            call dump_Mie ( details )
          case ( f_pfaFiles )
            call dump_PFAFileDatabase ( details )
          case ( f_pfaStru )
            call dump_PFAStructure ( details )
          case ( f_phaseName )
            call outputNamedValue ( 'phase name', currentphaseName )
          case ( f_pointingGrids )
            call dump_pointing_grid_database ( son )
          case ( f_stop )
            call finish ( 'ending mlsl2' )
            if ( switchDetail(switches, 'time') >= 0 ) then
              call output('(Now for the timings summary)', advance='yes')
              call dump_section_timings
            endif
            if ( NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave ) then
              write ( farewell, '(a,2(i0,a))' ) &
                & "Program stopped with normal status by /stop field on DUMP statement at line ", &
                & source/256, ", column ", mod(source,256), "."
              call MLSMessageExit( NORMAL_EXIT_STATUS, farewell=farewell )
            else if( parallel%slave ) then
              call closeParallel(0)
              write ( farewell, '(a,2(i0,a))' ) &
                & "Slave stopped by /stop field on DUMP statement at line ", &
                & source/256, ", column ", mod(source,256), "."
              call MLSMessageExit( farewell=farewell )
            else
              write ( farewell, '(a,2(i0,a))' ) &
                & "Program stopped by /stop field on DUMP statement at line ", &
                & source/256, ", column ", mod(source,256), "."
              call MLSMessageExit( farewell=farewell )
            endif
          case ( f_stopWithError )
            if ( switchDetail(switches, 'time') >= 0 ) then
              call output('(Now for the timings summary)', advance='yes')
              call dump_section_timings
            endif
            write ( farewell, '(a,2(i0,a))' ) &
                & "Program stopped by /stopWithError field on DUMP statement at line ", &
                & source/256, ", column ", mod(source,256), "."
            call MLSMessageExit( 1, farewell=farewell )
          end select
        end if
      case ( f_Boolean )
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        booleanString = lowerCase(booleanString)
        call output( trim(booleanString) // ' = ', advance='no' )
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & booleanString, label, countEmpty )
        call output( label, advance='yes' )
      case ( f_clean )
        clean = get_boolean(son)
        if ( clean ) optionsString = trim(optionsString) // 'c'
      case ( f_details )
        call expr ( gson, units, values, type )
        if ( units(1) /= phyq_dimensionless ) call AnnounceError ( gson, dimless )
        if ( type /= num_value ) call announceError ( gson, numeric )
        details = nint(values(1)) - DetailReduction
        ! call outputnamedValue( 'DetailReduction', DetailReduction )
        ! call outputnamedValue( 'Details', details )
      case ( f_file )
        if ( present(fileDataBase) ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            call get_string ( sub_rosa(gson), nameString, strip=.true. )
            oneMLSFile => getMLSFileByName ( fileDataBase, nameString )
            if ( associated(oneMLSFile) ) then
              call dumpMLSFile ( oneMLSFile )
            else
              call announceError ( gson, noFile, trim(nameString) )
            end if
          end do
        else
          call announceError ( son, noFileDatabase )
        end if
      case ( f_forwardModel ) ! Dump forward model configs
        if ( details < -1 ) cycle
        if ( haveForwardModelConfigs ) then
          do i = 2, nsons(son)
            call dump ( & ! has no details switch
              & forwardModelConfigs(decoration(decoration(subtree(i,son)))), &
              & quantityTemplatesDB=quantityTemplatesDB )
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
              call outputNamedValue ( ' Decoration ', vectorIndex )
              call output ( ' GriddedData: ' )
              call dump ( &
                & griddedDataBase(vectorIndex), details , options=optionsString )
            endif
          end do
        else
          call announceError ( gson, noGriddedData )
        end if
        GotFirst = .true.
      case ( f_height, f_surface )
        do i = 2, nsons(son)
          call expr ( subtree(i,son), units, values, type )
          if ( units(1) /= phyq_dimensionless ) call AnnounceError ( subtree(i,son), dimless )
          if ( type /= num_value ) call announceError ( subtree(i,son), numeric )
          height = values(1)
        end do
      case ( f_hessian ) ! Diff or Dump hessians
        if ( details < -1 ) cycle
        if ( haveHessians ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            ! The decoration is the negative of the index; see Fill, where
            ! the Hessian spec is processed.
            ! Well, that's a nasty trick.
            if ( gotFirst ) then
              hessianIndex2 = -decoration(decoration(gson))
            else
              hessianIndex = -decoration(decoration(gson))
            endif
            if ( DiffOrDump == s_diff ) then
              call diff ( HessianDatabase(hessianIndex), &
                &         HessianDatabase(hessianIndex2), details=details, &
                & options=optionsString )
            else
              call dump ( HessianDatabase(hessianIndex), details=details, &
                & options=optionsString )
            endif
          end do
        else
          call announceError ( gson, 0, 'Unable to dump Hessian here; db empty or absent' )
        end if
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
        call dumpL2PC ( fileDataBase(fileIndex), details, options=optionsString )
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
            qty2 => GetVectorQtyByTemplateIndex( &
                  & vectors(vectorIndex2), quantityIndex2 )
          else
            vectorIndex = decoration(decoration(subtree(1,gson)))
            quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
            qty1 => GetVectorQtyByTemplateIndex( &
                  & vectors(vectorIndex), quantityIndex )
          end if
          if ( DiffOrDump == s_diff ) then
            rmsFormat = '(1pe8.1)'
            if ( gotFirst ) then
              if ( index(optionsString, '1') > 0 ) then
                call diff ( qty1%values, 'qty1 values', &
                  & qty2%values, 'qty2 values', &
                  & options=optionsString )
              else if ( index(optionsString, '2') > 0 ) then
                call diff ( ichar(qty1%mask), 'qty1 mask', &
                  & ichar(qty2%mask), 'qty2 mask', &
                  & options=optionsString )
              else
                call diff ( &
                  & qty1, &
                  & qty2, &
                  & options=optionsString )
              endif
            endif
            rmsFormat = '*'
          else if ( fieldIndex == f_mask ) then
            if ( got(f_details) ) &
              & call outputNamedValue( 'Dump mask with details', details )
            if ( got(f_options) ) &
              & call outputNamedValue( 'Dump mask with options', optionsString )
            call dumpQuantityMask ( GetVectorQtyByTemplateIndex( &
              & vectors(vectorIndex), quantityIndex), &
              & details=details, options=optionsString )
          else
            ! Special options handling
            ! 'c' means clean
            ! '0' means dump template
            ! '1' means dump values
            ! '2' means dump mask

            if ( clean ) optionsString = trim(optionsString) // 'c'
            if ( .not. any(indexes(optionsString, (/'0','1','2'/)) > 0 ) ) then
              call dump ( qty1, details=details, &
                & vector=vectors(vectorIndex), options=optionsString )
            else if ( index(optionsString, '1') > 0 ) then
              call dump ( qty1%values, 'quantity values', &
                & options=optionsString )
            else if ( index(optionsString, '0') > 0 ) then
              call dump ( qty1%template, details=details )
            else
              call dumpQuantityMask( qty1, details=0, options=optionsString )
            end if
          end if
        end do
        GotFirst = .true.
      case ( f_matrix ) ! Diff or Dump matrices
        if ( details < -1 ) cycle
        if ( haveMatrices ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            if ( gotFirst ) then
              matrixIndex2 = decoration(decoration(gson))
            else
              matrixIndex = decoration(decoration(gson))
            end if
            if ( DiffOrDump /= s_diff ) then
              call getFromMatrixDatabase ( &
                & matrixDatabase(matrixIndex), matrix )
              call dump ( Matrix, details=details )
            else if ( GotFirst ) then
              call getFromMatrixDatabase ( &
                & matrixDatabase(matrixIndex), matrix )
              call getFromMatrixDatabase ( &
                & matrixDatabase(matrixIndex2), matrix2 )
              call diff ( Matrix, Matrix2, details=details )
            end if
          end do
        else
          call announceError ( gson, 0, 'Unable to dump matrix here; db empty or absent' )
        end if
        GotFirst = .true.
      case ( f_options )
        call get_string ( sub_rosa(gson), optionsString, strip=.true. )
        ! optionsString = lowerCase(optionsString)
        ! call outputNamedValue( 'options', trim(optionsString) )
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
            call cpu_time ( cpuTime ) ! In case %[CcSs] appears
            i = max(index(text,'%c'), index(text,'%C'))
            if ( i /= 0 ) then
              tempText = ''
              l = 1 ! Position in TempText
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
              write ( tempText, '(1p,g15.6)' ) cpuTime - cpuTimeBase
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

    subroutine AnnounceError ( where, what, string )
      use MORETREE, only: STARTERRORMESSAGE
      use OUTPUT_M, only: NEWLINE

      integer, intent(in) :: What, Where
      character(len=*), intent(in), optional :: String

      call StartErrorMessage ( where )

      select case ( what )
      case ( dimless )
        call output ( "The field is not unitless." )
      case ( noFile )
        call output ( "File " // string // " not in database." )
      case ( noFileDatabase )
        call output ( "File database not provided." )
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
      case ( noradiometers )
        call output ( "Can't dump Radiometers here." )
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
      case default
        if ( present(string) ) call output ( string )
        call output ( "No reserved error flag for this" )
      end select
      call newLine
    end subroutine AnnounceError

  end subroutine DumpCommand
  
  subroutine  MLSCase ( ROOT )
  ! Returns TRUE if the label or Boolean field of the last Select command
  ! matches the label or Boolean field of the current Case command
  ! or if the current Case command is given the special label 'default',
  ! a wildcard matching anything
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_LABEL, F_OPTIONS
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use MLSSTRINGLISTS, only: GETHASHELEMENT, SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE, STREQ
    use MORETREE, only: GET_FIELD_ID
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NSONS, SUB_ROSA, SUBTREE
    ! Args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND8'
    integer :: GSON, J
    integer :: FieldIndex
    character(len=80) :: Label  ! E.g., 'BAND8'
    character(len=8)  :: optionsString
    integer :: Son
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    MLSSelecting = .true. ! Defaults to skipping rest of case
    if ( MLSSelectedAlready ) return
    optionsString = ' '
    if ( toggle(gen) ) then
      call trace_begin ( 'MLSCase', root )
    else
      call MLSMessageCalls( 'push', constantName='MLSCase' )
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
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & booleanString, label, countEmpty )
      case (f_label)
        call get_string ( sub_rosa(gson), label, strip=.true. )
      case (f_options)
        call get_string ( sub_rosa(gson), optionsString, strip=.true. )
      case default
        ! Should not have got here if parser worked correctly
      end select
    enddo
    label = lowerCase(label)
    if ( label == 'default' ) then
      MLSSelecting = .false. ! Matches any pattern
    elseif ( len_trim(optionsString) < 1 ) then
      MLSSelecting = ( label /= selectLabel )
    else
      ! The streq function is a generalized string '=='
      MLSSelecting = .not. streq( label, selectLabel, optionsString )
    endif
    if ( verbose ) then
      call output( 'MLSCase label = ', advance='no' )
      call output( trim(label), advance='yes' )
    endif
    ! We must store whether we have ever had a match
    MLSSelectedAlready = MLSSelectedAlready .or. .not. MLSSelecting
    if ( toggle(gen) ) then
      call trace_end ( 'MLSCase' )
    else
      call MLSMessageCalls( 'pop' )
    end if
  end subroutine  MLSCase
  
  subroutine  MLSSelect ( ROOT )
  ! Fills the global variable selectLabel with
  ! the Label field or value of the Boolean field
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_LABEL
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use MLSSTRINGLISTS, only: GETHASHELEMENT, SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use MORETREE, only: GET_FIELD_ID
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NSONS, SUB_ROSA, SUBTREE
    ! Args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND8'
    integer :: GSON, J
    integer :: FieldIndex
    character(len=80) :: Label          ! E.g., 'BAND8'
    integer :: Son
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    if ( toggle(gen) ) then
      call trace_begin ( 'MLSSelect', root )
    else
      call MLSMessageCalls( 'push', constantName='MLSSelect' )
    end if
    label = ' '
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        booleanString = lowerCase(booleanString)
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & booleanString, label, countEmpty )
      case (f_label)
        call get_string ( sub_rosa(gson), label, strip=.true. )
      case default
        ! Should not have got here if parser worked correctly
      end select
    enddo
    if ( verbose ) then
      call output( 'Selecting label = ', advance='no' )
      call output( trim(label), advance='yes' )
    endif
    selectLabel = lowerCase(label)
    MLSSelecting = .true.
    if ( toggle(gen) ) then
      call trace_end ( 'MLSSelect' )
    else
      call MLSMessageCalls( 'pop' )
    end if
  end subroutine MLSSelect
  
  subroutine  MLSEndSelect ( ROOT )
  ! Resets the global variable MLSSelecting
  ! Optionally puts note about end of selecting in log file
    use INIT_TABLES_MODULE, only: F_LABEL
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMESSAGECALLS, MLSMSG_Info
    use MLSSTRINGS, only: LOWERCASE
    use MORETREE, only: GET_FIELD_ID
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NSONS, SUB_ROSA, SUBTREE
    ! Args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: Label  ! E.g., 'BAND8'
    integer :: GSON, J
    integer :: FieldIndex
    integer :: Son
    ! Executable
    MLSSelectedAlready = .false.
    if ( toggle(gen) ) then
      call trace_begin ( 'MLSEndSelect', root )
    else
      call MLSMessageCalls( 'push', constantName='MLSEndSelect' )
    end if
    label = ' '
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_label)
        call get_string ( sub_rosa(gson), label, strip=.true. )
        label = lowerCase(label)
      case default
        ! Should not have got here if parser worked correctly
      end select
    enddo
    MLSSelecting = .false.
    if ( len_trim(label) > 0 ) call MLSMessage (MLSMSG_Info, moduleName, &
        & 'End selecting for ' // trim(label) // '::' // trim(SelectLabel) )
    if ( toggle(gen) ) then
      call trace_end ( 'MLSEndSelect' )
    else
      call MLSMessageCalls( 'pop' )
    end if
  end subroutine MLSEndSelect
  
  logical function Skip ( ROOT )
    ! Returns value of Boolean field (if present);
    ! If TRUE should skip rest of section in which SKIP command appears
    ! If Boolean field absent, returns TRUE 
    ! (otherwise it would do nothing--
    ! this way it forces us to skip unconditionally)
    use INIT_TABLES_MODULE, only: F_BOOLEAN, F_FORMULA
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use MLSSTRINGLISTS, only: BOOLEANVALUE, SWITCHDETAIL
    use MLSSTRINGS, only: LOWERCASE
    use MORETREE, only: GET_FIELD_ID
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NSONS, SUB_ROSA, SUBTREE
    ! Args
    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Internal variables
    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND13_OK'
    integer :: GSON, J
    integer :: FieldIndex
    integer :: Son
    logical :: verbose, verboser
    ! Executable
    verbose = ( switchDetail(switches, 'bool') > -1 )
    verboser = ( switchDetail(switches, 'bool') > 0 )
    if ( toggle(gen) ) then
      call trace_begin ( 'Skip', root )
    else
      call MLSMessageCalls( 'push', constantName='Skip' )
    end if
    skip = .true. ! Defaults to skipping rest of section
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        booleanString = lowerCase(booleanString)
        skip = BooleanValue ( booleanString, &
          & runTimeValues%lkeys, runTimeValues%lvalues)
      case ( f_formula )
        call get_string ( sub_rosa(gson), booleanString, strip=.true. )
        skip = myBooleanValue ( booleanString )
      case default
        ! Should not have got here if parser worked correctly
      end select
      if ( verbose ) then
        call output( trim(booleanString) // ' = ', advance='no' )
        call output( skip, advance='yes' )
      endif
      if ( skip ) &
        & call output( '(Skipping rest of this section)', advance='yes' )
    enddo
    if ( toggle(gen) ) then
      call trace_end ( 'Skip' )
    else
      call MLSMessageCalls( 'pop' )
    end if
  end function Skip

! =====     Private Procedures     =====================================

  function myBooleanValue ( FORMULA ) result ( BVALUE )
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: BOOLEANVALUE, GETSTRINGELEMENT
    use MLSSTRINGS, only: LOWERCASE
    use OUTPUT_M, only: OUTPUTNAMEDVALUE
    ! Calculate the boolean value according to
    ! (1) The logical value of its formula, if the formula
    !     does not contain the special operators "==" or "/="
    ! (2) if the formula is "variable == value" and variable is
    !     is recognized and takes the value "value", return "true", else "false"
    ! (3) if the formula is "variable /= value" return not (2)
    character(len=*), intent(in) :: formula
    logical :: bvalue
    ! Internal variables
    character(len=16) :: lhs, rhs
    logical, parameter :: countEmpty = .false.
    logical :: reverse  ! Do we mean NOT equal?
    ! Executable
    bvalue = .false.
    reverse = .false.
    if ( index(formula, "==") > 1 ) then
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='=' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='=' )
    elseif ( index(formula, "/=") > 1 ) then
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='/' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='=' )
      reverse = .true.
    else
      bvalue = BooleanValue ( formula, &
        & runTimeValues%lkeys, runTimeValues%lvalues )
      return
    endif
    rhs = lowercase(adjustl(rhs))
    ! OK, so now let's try to make sense of the lhs
    ! Note the necessary type conversions for some variables
    ! that aren't chracter-valued
    ! (Should we have done this for rhs, too?)
    lhs = Evaluator(lowercase(adjustl(lhs)))
    bvalue = (lhs == rhs)
    call outputnamedvalue('lhs', trim(lhs) )
    call outputnamedvalue('rhs', trim(rhs) )
    if ( reverse ) bvalue = .not. bvalue  ! If we meant not equal
  end function myBooleanValue

  ! This evaluates a character-valued arg, being alert for special values
  ! that name global variables, e.g. 'phasename'
  function Evaluator ( ARG ) result( ITSVALUE )
    use MLSL2OPTIONS, only: CATENATESPLITS, CHECKPATHS, NEED_L1BFILES, & ! , SIPS_VERSION
      & SKIPRETRIEVAL
    use MLSL2TIMINGS, only: CURRENTCHUNKNUMBER, CURRENTPHASENAME
    use MLSSTRINGS, only: LOWERCASE, WRITEINTSTOCHARS
    ! Args
    character(len=*), intent(in) :: arg
    character(len=MAXRESULTLEN)  :: itsValue
    ! Executable
    select case (arg)
    case ('catenatesplits')
      itsValue = merge( 'true ', 'false', catenatesplits )
    case ('checkpaths')
      itsValue = merge( 'true ', 'false', checkpaths )
    case ('chunknumber')
      call writeIntsToChars ( currentChunkNumber, itsValue )
    case ('phasename')
      itsValue = lowercase(currentPhaseName)
    case ('sips_version')
      itsValue = 'true' ! merge( 'true ', 'false', sips_version )
    case ('need_l1bfiles')
      itsValue = merge( 'true ', 'false', need_l1bfiles )
    case ('skipretrieval')
      itsValue = merge( 'true ', 'false', skipretrieval )
    case default
      ! What did you mean?
      ! Maybe just whether two character strings are the same
      ! that were assembled using m4 trickery
      itsValue = arg
    end select
  end function Evaluator
  
  function BooleanToString ( BOOL ) result ( STR )
    ! Convert a Boolean argument to a character-valued string
    ! I.e., .true. -> 'true'
    !       .false. -> 'false'
    ! Args
    logical, intent(in) :: Bool
    character(len=6)    :: str
    str = merge( 'true ', 'false', Bool )
  end function BooleanToString

  ! ---------------------------------------------  returnFullFileName  -----
  subroutine returnFullFileName ( SHORTNAME, FULLNAME, &
    & PCF_START, PCF_END )
    use MLSFILES, only: GETPCFROMREF
    use MLSL2OPTIONS, only: TOOLKIT
    ! Given a possibly-abbreviated shortName, return the full name
    ! as found in the PCF
    ! (w/o toolkit panoply, simply return shortName)
    ! Args
    character(len=*), intent(in)  :: shortName
    character(len=*), intent(out) :: FullName
    integer, intent(in)           :: pcf_start
    integer, intent(in)           :: pcf_end
    ! Internal variables
    logical, parameter :: DEBUG = .false.
    integer :: FileHandle
    integer :: returnStatus
    integer :: Version
    ! Executable
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
! Revision 2.82  2013/03/30 00:19:48  vsnyder
! Add quantity database to forward model config dump
!
! Revision 2.81  2013/02/21 21:38:10  pwagner
! Pass options string when dumping quantity mask
!
! Revision 2.80  2013/02/12 18:17:44  pwagner
! Removed SIPS_VERSION; raised -Sbool switch needed for most printing
!
! Revision 2.79  2012/12/04 00:19:14  pwagner
! May dump timings summary after Dump, /stop
!
! Revision 2.78  2012/08/16 18:03:44  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.77  2012/06/27 17:53:04  pwagner
! May now Dump command line
!
! Revision 2.76  2012/06/06 20:20:46  vsnyder
! Add line number to farewell for /stop
!
! Revision 2.75  2012/05/14 22:26:31  pwagner
! Guard against missing swath--counts as empty
!
! Revision 2.74  2012/05/11 00:16:42  pwagner
! Added BooleanFromEmptySwath
!
! Revision 2.73  2012/05/08 17:48:37  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.72  2012/05/01 23:18:13  pwagner
! More capable Boolean formula; adds relations
!
! Revision 2.71  2012/04/26 23:28:33  pwagner
! May Dump chunk number, phase name; BooleanFromFormula can compare strings, variables
!
! Revision 2.70  2012/04/20 01:29:56  vsnyder
! Add call to Finish before stopping
!
! Revision 2.69  2012/04/20 01:03:30  pwagner
! May dump callStack from within l2cf
!
! Revision 2.68  2012/03/28 23:14:14  vsnyder
! Pass options to Dump_GH, add 'c' to options if /clean appears
!
! Revision 2.67  2012/03/15 22:50:39  vsnyder
! Add IGRF dump
!
! Revision 2.66  2012/02/28 00:16:42  vsnyder
! Repair %S, which was referencing cpu_time without it being assigned a value
!
! Revision 2.65  2012/02/02 00:56:00  pwagner
! Fix bug in Diffing matrices
!
! Revision 2.64  2011/11/04 00:29:28  pwagner
! :Fixed something only NAg complained about when Matrices not associated yet
!
! Revision 2.63  2011/08/29 22:12:13  pwagner
! Should not attempt to take size of disassociated pointer HessianDataBase
!
! Revision 2.62  2011/06/02 19:25:01  pwagner
! May dump allRadiometers
!
! Revision 2.61  2011/05/09 18:08:57  pwagner
! Print notice of changed runtime booleans only when "bool" switch is set
!
! Revision 2.60  2011/04/20 16:47:54  pwagner
! Added BooleanFromEmptyGrid
!
! Revision 2.59  2011/04/13 00:26:17  pwagner
! options='' works how comments and wiki says it does for vect.qty
!
! Revision 2.58  2011/04/04 23:07:21  pwagner
! May diff just quantity mask
!
! Revision 2.57  2011/03/15 23:02:45  pwagner
! Dump qty with options 0, 1, or 2 to choose template, values, or mask
!
! Revision 2.56  2011/03/08 18:28:06  pwagner
! May optionally diff, dump just values of vectorqtys
!
! Revision 2.55  2010/11/19 23:59:23  pwagner
! Passes options string to Hessians diff
!
! Revision 2.54  2010/11/05 22:37:06  pwagner
! Pass optionsString to dumps of L2PC, Hessians
!
! Revision 2.53  2010/08/13 22:08:57  pwagner
! May diff hessians, matrices
!
! Revision 2.52  2010/08/06 23:08:48  pwagner
! Pass Hessians, matrices to DumpCommand
!
! Revision 2.51  2010/04/17 01:44:26  vsnyder
! Add details to DumpL2PC calls, spiff up an error message
!
! Revision 2.50  2010/04/16 01:39:34  vsnyder
! Added /allFiles and file fields to dump file database
!
! Revision 2.49  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
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
