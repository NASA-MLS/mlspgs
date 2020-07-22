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

! Process a "dump" command. Or a flow control command.
! Or functions to set a run-time "Boolean" flag.
! (Should these latter functions be moved into a special run-time module?)
! Our wiki pages have a writeup on the l2cf control structures in
! https://mls.jpl.nasa.gov/team/wiki/index.php/L2cf_Control_Structures

  use Chunks_M, only: MLSChunk_T, Dump
  use HighOutput, only: Banner, BeVerbose, HeadLine, &
    & NumToChars, OutputNamedValue
  use MLSL2Options, only: MLSL2Message
  use Output_m, only: Beep, Blanks, NewLine, Output, RevertOutput, SwitchOutput

  implicit none
  private

  public :: BooleanFromAnyGoodRadiances
  public :: BooleanFromAnyGoodValues
  public :: BooleanFromCatchWarning
  public :: BooleanFromChunkStartsAfter
  public :: BooleanFromChunkEndsBefore
  public :: BooleanFromComparingQtys
  public :: BooleanFromEmptyGrid, BooleanFromEmptySwath
  public :: BooleanFromFormula
  public :: DumpCommand, ExecuteCommand
  public :: InitializeRepeat, NextRepeat
  public :: MLSCase, MLSEndSelect, MLSSelect, Skip

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
! The following functions do much the same thing: they add to the
! runtime Boolean database a pair (key => named_Boolean, value => its_value)
! where the value is determined by whether condition is true or false
!    function                 condition
! BooleanFromAnyGoodRadiances
!               Any of the radiances have non-negative precisions
! BooleanFromAnyGoodValues
!               The quantity has any useable values (based on precision, status)
! BooleanFromCatchWarning
!               The last command resulted in a (optionally specific) warning
! BooleanFromChunkEndsBefore
!               The chunk ends before the gridded quantity's last date
! BooleanFromChunkStartsAfter
!               The chunk starts after the gridded quantity's first date
! BooleanFromComparingQuantities
!               The first quantity stands in specified relation to the second
!               E.g., formula="all a > b"
! BooleanFromEmptyGrid
!               The specified GriddedData is empty
! BooleanFromEmptySwath
!               (a) That named swath in that named file has no useable data
!               (b) The specified file can't be found
! BooleanFromFormula
!               (a) Evaluate the formula; it may be one of the forms
!               (1) Contains only "or", "and", "not" => logical result
!               (2) Contains "lhs rel rhs" where
!                   rel is one of {==, /=, <, or >}  => logical result
!               (3) A manipulation formula, signaled by the "/literal' field
!                   (see FillByManipulation)  => string-valued result
!               (b) Just store the text of the label field
!
! The following subroutines depart from the sbove pattern
! DumpCommand    Process the Dump command, dumping any of the allowed datatypes
! ExecuteCommand Execute the shell command specified in text
! InitializeRepeat
!                Set the Repeat counter "count" to 0
! NextRepeat     Add 1 to Repeat counter "count",
!                  store values[count] in "counts(n)"
! MLSCase        Process the Case control statement
! MLSSelect      Process the Select control statement
! MLSEndSelect   Process the EndSelect control statement
! Skip           Process either the Skip or the Repeat control statement
! === (end of toc) ===

! === (start of api) ===
! int BooleanFromAnyGoodRadiances ( int root, type (chunk) chunk, &
!    type(MLSFile) filedatabase(:) )
! int BooleanFromAnyGoodValues ( int root, type (vector_t) vectors(:) )
! int BooleanFromCatchWarning ( int root )
! int BooleanFromComparingQtys ( int root, type (vector_t) vectors(:) )
! int BooleanFromChunkEndsBefore ( int root, type (mlschunk_t) chunk, &
!    type (griddeddata_t) grids(:) )
! int BooleanFromChunkStartsAfter ( int root, type (mlschunk_t) chunk, &
!    type (griddeddata_t) grids(:) )
! int BooleanFromEmptyGrid ( int root, type (griddeddata_t) grids(:) )
! int BooleanFromEmptySwath ( int root )
! int BooleanFromFormula ( int name, int root, type (vector_t) vectors(:) )
! DumpCommand ( int root, type (QuantityTemplate_T) quantityTemplatesDB(:), &
!    type (vectorTemplate_T) vectorTemplates(:), &
!    type (vector_t) vectors(:), type (forwardModelConfig_t) forwardModelConfigs(:), &
!    type (HGrids_T) HGrids(:), &
!    type (GriddedData_T) GriddedDataBase(:), type (MLSFile_T) FileDataBase(:), &
!    type (Matrix_Database_T) MatrixDataBase(:), type (Hessian_T) HessianDataBase(:) )
! ExecuteCommand ( int root )
! InitializeRepeat
! NextRepeat
! MLSCase ( int root )
! MLSSelect ( int root )
! MLSEndSelect ( int root )
! Skip ( int root, char* name )
! === (end of api) ===

  interface Evaluator
    module procedure Evaluator_Sca, Evaluator_Array
  end interface

  interface outputLater
    module procedure outputLater_chars, outputLater_int, outputLater_real
  end interface

  logical :: asBanner
  logical :: asHeadLine
  logical, parameter :: countEmpty = .true. ! Except where overriden locally
  character(len=8) :: headLineChars
  integer :: lineLength
  integer, parameter :: MAXRESULTLEN = 64
  logical, public, save             :: MLSSELECTING = .false.
  logical, save               :: MLSSelectedAlready = .false.
  character(len=MAXRESULTLEN), save :: selectLabel  = ' '

  character(len=80)  :: InputFileName ! e.g., '/tmp/pause.txt'
  character(len=255) :: Text
  logical :: toStderr
  integer, save :: TotalMemory1 = 0   ! For memory usage, in kB
  integer, save :: TotalMemory2 = 0   ! For memory usage, in kB
  real, save :: T1 = 0.               ! For timing
  real, save :: T2 = 0.               ! For timing

  ! Error codes
  integer, parameter :: NoFile = + 1
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
  integer, parameter :: Stop = noVT + 1
  integer, parameter :: Unknown = stop + 1 ! Unknown template

contains

  ! --------------------------------  BooleanFromAnyGoodRadiances  -----
  function BooleanFromAnyGoodRadiances ( root, chunk, filedatabase ) &
    & result( hashsize )
    use Allocate_Deallocate, only: Deallocate_Test
    use ConstructQuantityTemplates, only: AnyGoodSignalData
    use Init_Tables_Module, only: F_Signal, F_Boolean
    use MLSCommon, only: MLSFile_T
    use MLSL2Options, only: RuntimeValues
    use MLSSignals_M, only: GetSignalName, &
      & Signals
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use Parse_Signal_M, only: Parse_Signal
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    ! integer, intent(in) :: name
    integer, intent(in) :: root
    type (MLSChunk_T), intent(in) :: chunk
    type (MLSFile_T), dimension(:), pointer ::     Filedatabase
    integer             :: hashsize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: keyNo
    character(len=32) :: nameString
    logical :: reverse
    integer :: s
    integer :: signalIndex
    integer, pointer :: Signal_Indices(:)         ! Indices in the signals
    character(len=32) :: signalString
    integer :: son
    character(len=32) :: subSignalString
    logical :: tvalue
    logical :: verbose, verboser
    ! Executable
    reverse = .false.
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    nullify(Signal_Indices)
    ! call get_string(name, nameString)
    ! nameString = lowerCase(nameString)
    signalString = ' '
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
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
    ! The following could be expressed via: tvalue =
    !    ( reverse .and. tvalue ) .or. (.not. reverse .and. .not. tvalue )
    ! but is it clearer that way? I think .not.
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    hashsize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromAnyGoodRadiances

  ! -----------------------------------  BooleanFromAnyGoodValues  -----
  function BooleanFromAnyGoodValues ( root, vectors ) result( thesize )
    use Init_Tables_Module, only: F_Precision, F_Quality, &
      & F_Quantity, F_Boolean, F_Status, F_Reverse
    use ManipulateVectorQuantities, only: AnyGoodDatainQty
    use MLSKinds, only: Rv
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use MoreTree, only: Get_Boolean
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    use VectorsModule, only: Vector_T, VectorValue_T, &
      & GetVectorQtyByTemplateIndex
    ! Dummy args
    integer, intent(in) :: root
    type (vector_T), dimension(:) :: Vectors
    integer             :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: keyNo
    character(len=32) :: nameString
    type (VectorValue_T), pointer :: Precisionquantity
    integer  :: Quantityindex
    real(rv) :: Quality_min
    type (VectorValue_T), pointer :: Qualityquantity
    type (VectorValue_T), pointer :: Quantity
    logical :: reverse
    integer :: son
    integer :: source
    type (VectorValue_T), pointer :: Statusquantity
    logical :: tvalue
    integer :: Vectorindex
    logical :: verbose, verboser
    ! Executable
    reverse = .false.
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    nullify( precisionquantity, qualityquantity, Quantity, statusquantity )
    ! call get_string(name, nameString)
    ! nameString = lowerCase(nameString)
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
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
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    tvalue = AnyGoodDataInQty ( a=Quantity, &
      & precision=precisionQuantity, quality=qualityQuantity, &
      & status=statusQuantity, quality_min=quality_min )
    ! The following could be expressed via: tvalue =
    !    ( reverse .and. tvalue ) .or. (.not. reverse .and. .not. tvalue )
    ! but is it clearer that way? I think .not.
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromAnyGoodValues

  ! ------------------------------------  BooleanFromCatchWarning  -----
  function BooleanFromCatchWarning ( root ) result( size )
    ! Called to check if the last command resulted in a warning
    ! (either printed or suppressed)
    ! and optionally if the warning matches a supplied message string
    ! syntax:
    ! CatchWarning, [message='string'], Boolean="name"
    use Init_Tables_Module, only: F_Boolean, F_Message
    use MLSL2Options, only: RuntimeValues
    use MLSMessageModule, only: MLSMessageInquire
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    integer, intent(in) :: root
    integer             :: size
    ! Internal variables
    integer :: field
    integer :: field_index
    character(len=255) :: LastWarningMsg
    character(len=255) :: message
    integer :: keyNo
    character(len=32) :: nameString
    integer :: son
    logical :: tvalue
    logical :: verbose, verboser
    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    tvalue= .false.
    message = ' '
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
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
    end if
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
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromCatchWarning

  ! -----------------------------------  BooleanFromComparingQtys  -----
  function BooleanFromComparingQtys ( root, vectors ) result( thesize )
    use Expr_M, only: Expr
    use Init_Tables_Module, only: F_A, F_B, F_C, F_Boolean, F_Formula
    use MLSKinds, only: R8, Rv
    use MLSL2Options, only: RuntimeValues
    use MLSMessageModule, only: MLSMSG_Error
    use MLSStats1, only: MLSMax, MLSMin, MLSMean, MLSMedian
    use MLSStringLists, only: GetStringElement, NumStringElements, &
      & PutHashElement, ReplaceSubstring
    use MLSStrings, only: Lowercase
    use String_Table, only: Get_String
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
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
    type (VectorValue_T), pointer :: AQuantity
    type (VectorValue_T), pointer :: BQuantity
    real(rv) :: A, B, C                       ! constant "c" in formula
    integer :: field
    integer :: field_index
    character(len=8) :: flattening, arg(2), op
    character(len=255) :: formula
    character(len=255) :: formulaTemp
    integer :: keyNo
    integer :: Me = -1                   ! String index for trace cacheing
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
    call trace_begin ( me, "BooleanFromComparingQtys", cond=.false. )
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    nullify( aQuantity, bQuantity )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
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
      call MLSL2Message (MLSMSG_Error, moduleName, &
        & 'Formula in compare must be "flattening a op [b][c]".' )
    else if( index('<>=', trim(op)) < 1 ) then
      call MLSL2Message (MLSMSG_Error, moduleName, &
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
      call MLSL2Message (MLSMSG_Error, moduleName, &
        & 'Formula in compare found unrecognized op: ' // trim(op) )
    end select
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
    call trace_end ( "BooleanFromComparingQtys", cond=.false. )

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

    function STATFUN ( name, values, mask )
      ! Evaluate statistical function name of values
      ! masking if appropriate
      ! Args
      use MLSCommon, only: defaultUndefinedValue
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

  ! ------------------------------------- BooleanFromChunkEndsBefore --
  function BooleanFromChunkEndsBefore ( root, chunk, grids ) result( thesize )
    use Dump_0, only: Dump
    use Init_Tables_Module, only: F_Boolean, F_Grid, F_Reverse
    use GriddedData, only: GriddedData_T, Dump
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use MoreTree, only: Get_Boolean
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    integer, intent(in) :: root
    type (Griddeddata_T), dimension(:), pointer :: Grids
    type (MLSChunk_T)                           :: Chunk
    integer                                     :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    type (Griddeddata_T), pointer :: grid
    integer :: keyNo
    character(len=32) :: nameString
    logical :: reverse
    integer :: son
    logical :: tvalue
    integer :: value
    logical :: verbose, verboser
    ! Executable
    reverse = .false.
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
        value = subtree(2,son)
      else
        fieldValue = son
        value = 0
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_grid )
        grid => grids ( decoration ( decoration ( value ) ) )
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      end select
    end do
    tvalue = .true.
    if ( associated(grid) ) tvalue = (chunk%EndTime < grid%DateStarts(1))
    ! The following could be expressed via: tvalue =
    !    ( reverse .and. tvalue ) .or. (.not. reverse .and. .not. tvalue )
    ! but is it clearer that way? I think .not.
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      call headLine( 'Checking for chunk ending before grid' )
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
      call outputNamedValue( 'Chunk Ends at', chunk%EndTime )
      call outputNamedValue( 'Grid Starts at', grid%DateStarts(1) )
      call dump (chunk)
      if ( verboser ) call dump( grid, details=-2 )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromChunkEndsBefore

  ! ------------------------------------- BooleanFromChunkStartsAfter --
  function BooleanFromChunkStartsAfter ( root, chunk, grids ) result( thesize )
    use Dump_0, only: Dump
    use Init_Tables_Module, only: F_Boolean, F_Grid, F_Reverse
    use GriddedData, only: GriddedData_T, Dump
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use MoreTree, only: Get_Boolean
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    integer, intent(in) :: root
    type (Griddeddata_T), dimension(:), pointer :: Grids
    type (MLSChunk_T)                           :: Chunk
    integer                                     :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    type (Griddeddata_T), pointer :: grid
    integer :: keyNo
    character(len=32) :: nameString
    logical :: reverse
    integer :: son
    logical :: tvalue
    integer :: value
    logical :: verbose, verboser
    ! Executable
    reverse = .false.
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
        value = subtree(2,son)
      else
        fieldValue = son
        value = 0
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_grid )
        grid => grids ( decoration ( decoration ( value ) ) )
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      end select
    end do
    tvalue = .true.
    if ( associated(grid) ) tvalue = (chunk%StartTime > grid%DateEnds(1))
    ! The following could be expressed via: tvalue =
    !    ( reverse .and. tvalue ) .or. (.not. reverse .and. .not. tvalue )
    ! but is it clearer that way? I think .not.
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      call headLine( 'Checking for chunk starting after grid' )
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
      call outputNamedValue( 'Chunk Starts at', chunk%StartTime )
      call outputNamedValue( 'Grid Ends at', grid%DateEnds(1) )
      call dump (chunk)
      if ( verboser ) call dump( grid, details=-2 )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromChunkStartsAfter

  ! ------------------------------------- BooleanFromEmptyGrid --
  function BooleanFromEmptyGrid ( root, grids ) result( thesize )
    use Dump_0, only: Dump
    use Init_Tables_Module, only: F_Boolean, F_Grid, F_Reverse
    use GriddedData, only: GriddedData_T, Dump
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use MoreTree, only: Get_Boolean
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    integer, intent(in) :: root
    type (Griddeddata_T), dimension(:), pointer :: Grids
    integer             :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    type (Griddeddata_T), pointer :: grid
    integer :: keyNo
    character(len=32) :: nameString
    logical :: reverse
    integer :: son
    logical :: tvalue
    integer :: value
    logical :: verbose, verboser
    ! Executable
    reverse = .false.
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
        value = subtree(2,son)
      else
        fieldValue = son
        value = 0
      end if
      field_index = decoration(field)

      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_grid )
        grid => grids ( decoration ( decoration ( value ) ) )
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      end select
    end do
    tvalue = .true.
    if ( associated(grid) ) tvalue = grid%empty
    ! The following could be expressed via: tvalue =
    !    ( reverse .and. tvalue ) .or. (.not. reverse .and. .not. tvalue )
    ! but is it clearer that way? I think .not.
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      if ( verboser ) call headLine( 'Checking for empty grid' )
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
      if ( verboser ) call dump( grid, details=-2 )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromEmptyGrid

  ! ------------------------------------- BooleanFromEmptySwath --
  ! Returns TRUE if, according to use case (a) or (b)
  ! (a) there are no useable data points in the swath
  ! (or if the swath is not in the file at all)
  ! The useablility criterion is
  ! (1) Precision non-negative; and
  ! (2) Status even
  ! If even one point is useable (a very low bar, admittedly) then
  ! return FALSE
  ! or
  ! (b) the file named in the file= field does not exist at all
  function BooleanFromEmptySwath ( root ) result( thesize )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Dump_0, only: Dump
    use Init_Tables_Module, only: F_Boolean, F_File, F_NoPCFid, &
      & F_Reverse, F_Swath, F_Type, &
      & Field_First, Field_Last, &
      & L_DAO, L_GEOS5, L_GEOS5_7, L_Merra, L_Merra_2, L_L2dgg, L_L2gp, L_None
    use L2GPData, only: L2GPData_T, Rgp, L2GPNameLen, &
      & ReadL2GPData, DestroyL2GPContents
    use MLSCommon, only: Filenamelen
    use MLSFiles, only: HDFversion_5, MLS_Exists
    use MLSHDFeos, only: MLS_Swath_In_File
    use MLSMessageModule, only: MLSMSG_Warning
    use MLSL2Options, only: RuntimeValues
    use MLSPCF2, only: MLSPCF_L2GP_End, MLSPCF_L2GP_Start, &
      & MLSPCF_L2DGG_Start, MLSPCF_L2DGG_End, &
      & MLSPCF_L2DAO_Start, MLSPCF_L2DAO_End, &
      & MLSPCF_L2GEOS5_Start, MLSPCF_L2GEOS5_End
    use MLSStringLists, only: NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase
    use MoreTree, only: Get_Boolean
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    ! Dummy args
    integer, intent(in) :: root
    integer             :: thesize
    ! Internal variables
    integer :: field
    integer :: field_index
    integer :: fieldValue
    integer :: fileType
    character (len=FileNameLen) :: File_base    ! From the FILE= field
    character(len=FileNameLen) :: filename          ! filename
    logical, dimension(field_first:field_last) :: got
    character(len=L2GPNAMELEN) :: swathname         ! swathname
    integer :: i
    integer :: keyNo
    character(len=32) :: nameString
    logical :: noPCFid
    logical :: reverse
    integer :: son
    logical :: tvalue
    integer :: value
    type (L2GPData_T) :: l2gp
    integer :: numGood
    logical, dimension(:), pointer  :: negativePrec => null() ! true if all prec < 0
    logical, dimension(:), pointer  :: oddStatus => null() ! true if all status odd
    logical :: verbose, verboser
    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    noPCFid = .false.
    reverse = .false.
    got     = .false.
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
        value = subtree(2,son)
      else
        fieldValue = son
        value = 0
      end if
      field_index = decoration(field)
      got(field_Index) = .true.
      select case ( field_index )
      case ( f_Boolean )
        call get_string( sub_rosa(subtree(2,son)), nameString )
      case ( f_file )
        call get_string ( sub_rosa(subtree(2,son)), file_base, strip=.true. )
      case ( f_noPCFid )
        noPCFid = get_Boolean(son)
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      case ( f_swath )
        call get_string ( sub_rosa(subtree(2,son)), swathname, strip=.true. )
      case ( f_type )
        filetype = decoration(subtree(2, son))
      end select
    end do
    if ( .not. got(f_type) .or. noPCFid ) fileType = l_none
    select case (fileType)
    case ( l_l2dgg )
      call returnFullFileName( file_base, Filename, &
        & MLSPCF_l2dgg_start, MLSPCF_l2dgg_end )
    case ( l_l2gp )
      call returnFullFileName( file_base, Filename, &
        & MLSPCF_l2gp_start, MLSPCF_l2gp_end )
    case ( l_dao )
      call returnFullFileName( file_base, Filename, &
        & mlspcf_l2dao_start, mlspcf_l2dao_end )
    case ( l_geos5, l_geos5_7, l_merra, l_merra_2 ) ! --- GMAO Data (GEOS5*)
      call returnFullFileName( file_base, Filename, &
        & mlspcf_l2geos5_start, mlspcf_l2geos5_end )
    case ( l_none )
      Filename = file_base
    case default
      call MLSL2Message( MLSMSG_Warning, ModuleName, &
      & "type should have been either l2gp or dgg; assume you meant dgg" )
    end select
    tvalue = .true.
    if ( fileType == l_none ) then
      tvalue = ( mls_exists ( Filename ) /= 0 )
    elseif ( mls_swath_in_file( filename, swathname, HDFVERSION_5 ) ) then
      call output( 'Reading swath ' // trim(swathname), advance='yes' )
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
      end if
      call deallocate_test( negativePrec, 'negativePrec', ModuleName )
      call deallocate_test( oddStatus, 'oddStatus', ModuleName )
      call DestroyL2GPContents ( l2gp )
    end if
    if ( reverse ) tvalue = .not. tvalue
    if ( verbose ) then
      call output( trim(nameString) // ' = ', advance='no' )
      call output( tvalue, advance='yes' )
    end if
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & lowercase(trim(nameString)), BooleanToString(tvalue), &
      & countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    thesize = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verboser ) &
      & call dumpBooleans
  end function BooleanFromEmptySwath

  ! ------------------------------------------- BooleanFromFormula -----
  function BooleanFromFormula ( name, root, vectors ) result( size )
    ! Called either when a Boolean is first declared
    ! syntax:
    ! name: Boolean, formula="formula"
    !
    ! or when it is reevaluated
    ! syntax:
    ! Reevaluate, formula="formula", Boolean="name",
    ! [a=a.qty], [b=b.qty], [c=value], [values=..], [/literal]
    use Dump_0, only: Dump
    use Expr_M, only: Expr
    use Init_Tables_Module, only: F_A, F_B, F_Boolean, F_C, &
      & F_Evaluate, F_Expr, &
      & F_Formula, F_InputBoolean, F_Truncate, &
      & F_Label, F_Literal, F_Manipulation, F_Values, &
      & Field_First, Field_Last
    use Manipulationutils, only: Manipulate
    use MLSKinds, only: R8, Rv
    use MLSL2Options, only: RuntimeValues
    use MLSMessageModule, only: MLSMSG_Error
    use MLSStringLists, only: EvaluateFormula, GetHashElement, &
      & InsertHashElement, NumStringElements, PutHashElement
    use MLSStrings, only: Lowercase, ReadNumsFromChars, WriteIntsToChars
    use Moretree, only: Get_Boolean
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate
    use String_Table, only: Get_String
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree
    use VectorsModule, only: Vector_T, VectorValue_T, &
      & CloneVectorQuantity, CreateVectorValue, &
      & DestroyVectorQuantityValue, Dump, &
      & GetVectorQtyByTemplateIndex
    ! Dummy args
    integer, intent(in) :: Name
    integer, intent(in) :: Root
    type (vector_T), dimension(:), optional, target :: VECTORS
    integer             :: SIZE
    ! Internal variables
    character(len=80) :: Booleanstring  ! E.g., 'BAND8'
    real(rv) :: C                       ! constant "c" in formula
    type (VectorValue_T), pointer :: AQuantity
    type (VectorValue_T), target  :: ATarget
    type (VectorValue_T), pointer :: BQuantity
    character(len=32) :: cnameString
    character(len=32) :: cvalue
    logical :: evaluate
    integer :: field
    integer :: field_index
    integer :: fieldValue
    character(len=255) :: formula
    logical, dimension(field_first:field_last) :: GOT
    logical :: truncate
    integer :: J
    integer :: keyNo
    character(len=80) :: Keystring
    logical :: literal
    character(len=32) :: nameString
    integer :: Nvalues
    integer :: Quantityindex
    integer :: son
    integer :: source
    type (VectorValue_T) :: TQuantity ! Temporary
    type(QuantityTemplate_T) :: TTemplate
    logical :: tvalue
    integer :: value_field
    integer, dimension(2)  :: Unitasarray ! From expr
    real(r8), dimension(2) :: Valueasarray ! From expr
    integer :: Vectorindex
    logical :: verbose, verboser
    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    nullify( aQuantity, bQuantity )
    c = 0._rv
    evaluate = .false.
    truncate = .false.
    literal = .false.
    tvalue = .false.
    got = .false.
    if ( name > 0 ) then
      call get_string(name, nameString)
      nameString = lowerCase(nameString)
    end if
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(decoration(subtree(2,son))) ! The field's value
        source = subtree(2,son) ! required to be an n_dot vertex
      else
        fieldValue = son
        source = son
      end if
      field_index = decoration(field)
      got(field_Index) = .true.

      select case ( field_index )
      case ( f_a )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        aQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
        call dump( aQuantity )
      case ( f_b )
        VectorIndex = decoration(decoration(subtree(1,source)))
        QuantityIndex = decoration(decoration(decoration(subtree(2,source))))
        bQuantity => GetVectorQtyByTemplateIndex( &
          & vectors(VectorIndex), QuantityIndex )
        call dump( bQuantity )
      case ( f_c )
        call expr ( source, unitAsArray, valueAsArray )
        ! c is a numeric value
        c = valueAsArray(1)
        if ( verbose ) call outputNamedValue( 'c', c )
      case ( f_inputBoolean )
        call get_string ( sub_rosa(subtree(2,son)), cnameString, strip=.true. )
        cnameString = lowerCase(cnameString)
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & cnameString, cvalue, countEmpty=countEmpty, &
          & inseparator=runTimeValues%sep )
        call readNumsFromChars ( cvalue, c )
        if ( verbose ) call outputNamedValue( 'c', c )
      case ( f_Boolean )
        call get_string ( sub_rosa(subtree(2,son)), nameString, strip=.true. )
        nameString = lowerCase(nameString)
      case ( f_evaluate )
        ! call output( 'processing evaluate', advance='yes' )
        evaluate = get_Boolean ( fieldValue )
        if ( verbose ) call outputNamedValue( 'evaluate', evaluate )
      case ( f_expr )
        call MLSL2Message ( MLSMSG_Error, moduleName, &
          & 'No code yet to handle "expr" field.' )
      case ( f_formula, f_manipulation )
        call get_string ( sub_rosa(subtree(2,son)), formula, strip=.true. )
        tvalue = myBooleanValue (formula)
        if ( verbose ) call outputNamedValue ( 'tvalue', trim(formula) )
      case ( f_truncate )
        ! call output( 'processing truncate', advance='yes' )
        truncate = get_Boolean ( fieldValue )
        if ( verbose ) call outputNamedValue( 'truncate', truncate )
      case ( f_label )
        call get_string ( sub_rosa(subtree(2,son)), formula, strip=.true. )
        literal = .true.
      case ( f_literal )
        ! call output( 'processing literal', advance='yes' )
        literal = get_Boolean ( fieldValue )
        if ( verbose ) call outputNamedValue( 'literal', literal )
      case ( f_values )
        value_field = son
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    ! Should have exactly one of
    ! [f_values, f_formula, f_manipulation, f_label]
    if ( count(got( (/f_values, f_formula, f_manipulation, f_label/) ) ) > 1 ) &
      & call announceError ( son, 0, &
      & ' Must supply exactly one of (values, formula, manipulation, label)' )
    if ( verbose .and. len_trim(formula) > 0 ) &
      & call outputNamedValue ( 'formula', trim(formula) )
    if ( evaluate .and. got(f_formula) ) then
      formula = EvaluateTermwise( formula )
      if ( verbose ) call outputNamedValue ( 'formula', trim(formula) )
    end if
    if ( got(f_formula) ) then
      ! Are we in the middle of a Repeat loop? If so, replace each "${n}" with
      ! "${count}"
       call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
         & 'count', cvalue, countEmpty=countEmpty, &
         & inseparator=runTimeValues%sep )
       call readNumsFromChars ( cvalue, c )
       if ( verboser ) call outputnamedValue ( 'c', c )
       if ( c > 0._rv ) &
         & formula = EvaluateFormula ( formula, (/cvalue/), (/"n"/) )
      ! The next will explicitly evaluate run-time Booleans evoked by name via
      ! "${name}"
      if ( verboser ) call outputNamedValue ( 'formula (2nd)', trim(formula) )
      formula =  EvaluateExplicitly( formula )
      if ( verboser ) call outputNamedValue ( 'formula (3rd)', trim(formula) )
    end if
    nvalues = 0
    if ( got(f_values) ) then
      nvalues = nsons(value_field)-1
    end if
    if ( nvalues == 1 ) then
      ! We just treat this as a normal definition or reevaluation of
      ! the named Boolean treated as a scalar
      call get_string( sub_rosa(subtree(2, value_field)), BooleanString, &
        & strip=.true. )
      call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
        & lowercase(trim(nameString)), BooleanString, &
        & countEmpty=countEmpty, &
        & inseparator=runTimeValues%sep )
    elseif ( nvalues > 1 ) then
      ! Given an array of values we store them according to the scheme:
      !        key             value
      !      "name(1)"      "values(1)"
      !      "name(2)"      "values(2)"
      !           .   .   .   .   .
      !  "name(nvalues)"   "values(nvalues)"
      ! where "name" is the named Boolean
      call writeIntsToChars( nvalues, BooleanString )
      keyString = trim(nameString) // 'n'
      call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
        & keyString, BooleanString, &
        & countEmpty=countEmpty, &
        & inseparator=runTimeValues%sep )
      if ( verboser ) then
        call outputnamedValue( 'keyString', trim(keyString) )
        call outputnamedValue( 'BooleanString', trim(BooleanString) )
      end if
      ! see also how the Repeat command treats the values field
      do j=2, nvalues+1
        call writeIntsToChars( j-1, keyString )
        keyString = trim(nameString) // '(' // trim(adjustl(keystring)) // ')'
        call get_string( sub_rosa(subtree(j, value_field)), BooleanString, strip=.true. )
        if ( verboser ) then
          call outputnamedValue( 'keyString', trim(keyString) )
          call outputnamedValue( 'BooleanString', trim(BooleanString) )
        end if
        call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
          & keyString, BooleanString, &
          & countEmpty=countEmpty, &
          & inseparator=runTimeValues%sep )
      enddo
    elseif ( got(f_manipulation) ) then
      ! print *, 'Oops-you dummy! code this missing piece'
      if ( .not. present(vectors) ) &
        & call MLSL2Message ( MLSMSG_Error, moduleName, &
        & 'Cant set manipulation from this section--only Fill or Retrieval' )
      if ( .not. got ( f_a ) ) then
        call setupNewquantitytemplate( tTemplate )
        aQuantity => aTarget
        aQuantity%template = tTemplate
        call CreateVectorValue ( aQuantity, 'a' )
      end if
      call CloneVectorQuantity ( TQuantity, aQuantity )
      if ( verboser ) then
        call output( trim(nameString) // ' = ', advance='no' )
        call output( trim(formula), advance='yes' )
      end if
      call Manipulate( TQuantity, AQuantity, BQuantity, C, formula, &
        & .false., '' )
      formula = numToChars ( TQuantity%values(1,1) )
      if ( truncate ) formula = numToChars ( int(TQuantity%values(1,1)) )
      call insertHashElement ( nameString, formula, &
        & runTimeValues%lkeys, runTimeValues%lvalues, &
        & inseparator=runTimeValues%sep )
      call destroyVectorQuantityValue ( TQuantity, &
        & destroyMask=.true., destroyTemplate=.false. )
    elseif ( literal ) then
      ! print *, 'Oops-you dummy! code this missing piece'
      if ( verboser ) then
        call output( trim(nameString) // ' = ', advance='no' )
        call output( trim(formula), advance='yes' )
      end if
      call insertHashElement ( nameString, formula, &
        & runTimeValues%lkeys, runTimeValues%lvalues, &
        & inseparator=runTimeValues%sep )
    else
      if ( verbose ) call output( 'Using default value for bare Boolean declaration', &
        & advance='yes' )
      if ( verboser ) then
        call output( trim(nameString) // ' = ', advance='no' )
        call output( tvalue, advance='yes' )
      end if
      call insertHashElement ( nameString, BooleanToString(tvalue), &
        & runTimeValues%lkeys, runTimeValues%lvalues, &
        & inseparator=runTimeValues%sep )
    end if
    size = NumStringElements( runTimeValues%lkeys, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    if ( verbose ) &
      & call dumpBooleans
  end function BooleanFromFormula

  ! ------------------------------------------------  DumpCommand  -----
  subroutine DumpCommand ( root, quantityTemplatesDB, &
    & vectorTemplates, vectors, forwardModelConfigs, hGrids, griddedDatabase, &
    & fileDatabase, matrixDatabase, hessianDatabase )

  ! Process a "dump" command

    use AntennaPatterns_M, only: Dump_Antenna_Patterns_Database
    use BitStuff, only: DumpBitCounts, DumpBitNames
    use Calendar, only: Duration_Formatted, Time_T, Tk
    use Call_Stack_M, only: Dump_Stack
    use ChunkDivideConfig_M, only: ChunkDivideConfig, Dump
    use Declaration_Table, only: Dump_A_Decl, Decls, Get_Decl, Variable
    use Diff_1, only: Diff
    use Dump_Options, only: PrintNameAsHeadline, PrintNameInBanner, RMSFormat
    use Dump_0, only: Dump
    use Expr_M, only: Expr
    use FilterShapes_M, only: Dump_Filter_Shapes_Database, &
      & Dump_Dacs_Filter_Database
    use Forwardmodelconfig, only: Dump, Forwardmodelconfig_T
    use GriddedData, only: Diff, Dump, GriddedData_T
    use HessianModule_1, only: Hessian_T, Diff, Dump
    use HGridsDatabase, only: Dump, HGrids_T
    use Igrf_Int, only: Dump_Gh
    use Init_Tables_Module, only: F_AllBooleans, F_AllFiles, &
      & F_Allforwardmodels, F_AllgriddedData, F_Allhessians, F_Allhgrids, &
      & F_Alll2pcs, F_Alllines, F_Allmatrices, F_Allpfa, &
      & F_AllquantityTemplates, F_Allradiometers, F_Allsignals, F_Allspectra, &
      & F_AllVectors, F_AllVectorTemplates, F_Allvgrids, F_Antennapatterns, &
      & F_Block, F_Boolean, &
      & F_Callstack, F_Chunkdivide, F_ChunkNumber, F_Clean, &
      & F_Commandline, F_Count, F_Crashburn, &
      & F_Details, F_Dacsfiltershapes, F_DumpFile, &
      & F_File, F_Filtershapes, F_Forwardmodel, F_GlobalAttributes, F_Grid, &
      & F_Hessian, F_Hgrid, F_Igrf, F_InputFile, &
      & F_L2pc, F_Lines, F_Mark, F_Mask, F_Matrix, F_Memory, &
      & F_Mietables, F_Options, &
      & F_Pause, F_PfaData, F_PfaFiles, F_Pfanum, F_Pfastru, &
      & F_Phasename, F_Pointinggrids, F_Polygon, F_Quantity, F_Rank, F_Reset, &
      & F_Signals, F_Spectroscopy, F_Stack, F_Start, F_Stop, F_Stride, &
      & F_Stopwitherror, F_Surface, &
      & F_Template, F_Text, F_Tgrid, F_Time, F_Truncate, &
      & F_TotalMatrixSizes, F_TotalVectorSizes, F_Variable, &
      & F_Vector, F_VectorMask, F_Vgrid, F_ZOT, &
      & S_Diff, S_Dump, S_Quantity, S_VectorTemplate, &
      & Field_First, Field_Last
    use Io_Stuff, only: Truncate_TextFile
    use L2Parinfo, only: Parallel, Closeparallel
    use L2PC_M, only: L2pcDatabase, Dumpl2pc => Dump
    use Lexer_Core, only: Get_Where, Where_T
    use Machine, only: Nevercrash
    use MatrixModule_1, only: Matrix_T, Matrix_Database_T, &
      & Diff, Dump, GetFromMatrixDatabase, MatricesMemoryInUse
    use MLSCommon, only: MLSFile_T
    use MLSFiles, only: DumpMLSFile => Dump, GetMLSFilebyname, MLS_Exists
    use MLSFinds, only: Findfirst, FindUnique
    use MLSKinds, only: R8, Rv
    use MLSL2Options, only: L2Options, &
      & L2cfnode, Normal_Exit_Status, RuntimeValues  !,  MLSMessage
    use MLSL2Timings, only: Dump_Section_Timings
    use MLSMessagemodule, only: MLSMSG_Crash, MLSMSG_Error, MLSMessageCalls, &
      & MLSMessageExit
    use MLSSignals_M, only: Radiometers, Signals, &
      & Dump, Dump_All, Getradiometerindex
    use MLSStrings, only: Indexes, LowerCase, &
      & ReadIntsFromChars, Stretch, WriteIntsToChars
    use MLSStringLists, only: GetHashElement, OptionDetail, RemoveOption, &
      & SwitchDetail
    use Moretree, only: Get_Boolean, Get_Field_Id, Get_Spec_Id
    use PCFHdr, only: DumpGlobalAttributes
    use PfaDatabase_M, only: Dump, Dump_PfaDatabase, Dump_PfaFileDatabase, &
      & Dump_Pfastructure, PfaData
    use PointingGrid_M, only: Dump_Pointing_Grid_Database
    use Polygon_M, only: Dump_Polygon_Database
    use QuantityTemplates, only: Dump, QuantityTemplate_T
    use Read_Mie_M, only: Dump_Mie
    use SpectroscopyCatalog_M, only: Catalog, Dump, Dump_Lines_Database, Lines
    use String_Table, only: Display_String, Get_String
    use Time_M, only: Time_Config
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Node_Id, Nsons, Sub_Rosa, Subtree, Where
    use Tree_Types, only: N_Spec_Args
    use VectorsModule, only: Vector_T, VectorTemplate_T, VectorValue_T, &
      & Diff, Dump, DestroyVectorQuantityValue, DumpQuantityMask, DumpVectorMask, &
      & GatherVectorQuantity, GetVectorQuantity, GetVectorQtyByTemplateIndex, &
      & GetVectorQuantityIndexByName, VectorsMemoryInUse
    use VGridsDatabase, only: Dump, VGrids
    use Wait_M, only: Wait_For_Event, Pause

    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Databases:
    type (QuantityTemplate_T), dimension(:), pointer, optional   :: QuantityTemplatesDB
    type (forwardModelConfig_t), dimension(:), pointer, optional :: ForwardModelConfigs
    type (vectorTemplate_T), dimension(:), pointer, optional     :: VectorTemplates
    type (vector_T), dimension(:), target, optional              :: Vectors
    type (HGrids_T), dimension(:), pointer, optional             :: HGrids
    type (griddedData_T), dimension(:), pointer, optional        :: griddedDataBase
    type (MLSFile_T), dimension(:), pointer, optional            :: FileDataBase
    type (matrix_Database_T), dimension(:), pointer, optional    :: MatrixDataBase
    type (Hessian_T), dimension(:), pointer, optional            :: HessianDataBase

    character(len=80) :: Booleanstring  ! E.g., 'BAND8'
    logical :: Clean
    real(tk) :: CPUTime, CPUTimeBase = 0.0_tk
    character(8) :: Date
    type(decls) :: Decl ! of an l2cf Variable
    integer :: DetailReduction
    integer :: Details
    integer :: DiffOrDump
    logical :: doStretch
    logical :: doStretchier
    character(len=255) :: Farewell
    integer :: FieldIndex
    integer :: FileIndex
    character(8) :: fieldName
    logical, dimension(field_first:field_last) :: GOT
    logical :: GotFirst ! of something -- needed if diffing 2 of them
    logical :: GotOne ! of something -- used to test loop completion
    integer :: GSON, I, J, K, L, Look
    logical :: HaveForwardModelConfigs, HaveGriddedData, &
             & HaveHessians, HaveHGrids, &
             & HaveMatrices, HaveQuantityTemplatesDB, &
             & HaveVectors, HaveVectorTemplates
    real(rv) :: height  ! We will use this to dump just one surface
    integer :: hessianIndex
    integer :: hessianIndex2
    integer :: HSLABRANK ! Rank of hyperslab; i.e. size(count)
    integer :: ii
    integer :: JJ
    character(len=80) :: Keystring  ! E.g., 'BAND13_OK'
    character(len=80) :: Label  ! E.g., 'BAND8'
    character(len=8) :: lineLenChars
    type (Matrix_T), pointer :: matrix
    type (Matrix_T), pointer :: matrix2
    integer :: MatrixIndex
    integer :: MatrixIndex2
    integer :: Me = -1               ! String index for trace cacheing
    integer :: MUL
    integer :: N
    character(len=80) :: mesg
    character(len=80) :: Namestring  ! E.g., 'L2PC-band15-SZASCALARHIRES'
    integer :: nUnique
    logical :: oldPrintNameAsHeadline
    logical :: oldPrintNameInBanner
    type (MLSFile_T), pointer :: OneMLSFile
    character(len=80) :: optionsString  ! E.g., '-rbs' (see dump_0.f90)
    type (VectorValue_T), pointer :: Quantity
    integer :: QuantityIndex
    integer :: QuantityIndex2
    type (VectorValue_T), pointer :: Qty1, Qty2
    type (VectorValue_T) :: VectorValue
    logical :: reset
    integer :: rank
    integer :: Son
    integer, dimension(3) :: start, count, stride, block
    integer :: Startnode
    character(len=80) :: TempOptionsString  ! E.g., '-rbs' (see dump_0.f90)
    character(20) :: TempText
    type(time_t) :: Time
    character(10) :: TimeOfDay
    logical :: truncate
    integer, dimension(1000) :: uniqueVals
    type (vector_T), pointer  :: Vector
    integer :: VectorIndex
    integer :: VectorIndex2
    logical :: verbose
    integer :: Units(2) ! of the Details expr -- known to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr
    integer, dimension(3) :: UNITASARRAY ! From expr
    real(r8), dimension(3) :: VALUEASARRAY ! From expr
    integer :: What
    type(where_t) :: Where_At
    logical :: ZOT
    character(len=16), dimension(10), parameter :: StatusBitNames = (/ &
      & 'do not use      ', 'beware          ', 'inform          ', &
      & 'post-processed  ', 'high clouds     ', 'low clouds      ', &
      & 'no meteorology  ', 'abandoned chunk ', 'too few radiance', &
      & 'mlsl2 crashed   ' /)
  !      123456789012345678901234567890123456789012345678901234567890

    ! Executable
    call trace_begin ( me, 'DumpCommand', root, cond=toggle(gen) )
    verbose = BeVerbose ( 'bool', -1 )
    nullify ( vector, quantity )
    ! Were we called to do a diff or a dump?
    ! The following must be one of (/ s_dump, s_diff /)
    DiffOrDump = get_spec_id(root)
    if ( .not. any( DiffOrDump == (/ s_dump, s_diff /) ) ) &
    &  call MLSL2Message( MLSMSG_Error, moduleName, &
        & "Expected either Diff or Dump command to get here")
    haveForwardModelConfigs = present(forwardModelConfigs)
    if ( haveForwardModelConfigs ) &
      & haveForwardModelConfigs = associated(forwardModelConfigs)
    haveGriddedData = present(griddedDataBase)
    if ( haveGriddedData ) haveGriddedData = associated(griddedDataBase)
    HaveHessians = present(HessianDatabase)
    if ( HaveHessians ) then
      if ( associated(HessianDatabase) ) HaveHessians = size(HessianDatabase) > 0
    end if
    haveHGrids = present(hGrids)
    if ( haveHGrids ) haveHGrids = associated(hGrids)
    haveMatrices = present(MatrixDatabase)
    if ( haveMatrices ) haveMatrices = associated(MatrixDatabase)
    if ( haveMatrices ) haveMatrices = size(MatrixDatabase) > 0
    haveQuantityTemplatesDB = present(quantityTemplatesDB)
    if ( haveQuantityTemplatesDB ) &
      & haveQuantityTemplatesDB = associated(quantityTemplatesDB)
    haveVectorTemplates = present(vectorTemplates)
    if ( haveVectorTemplates ) &
      & haveVectorTemplates = associated(vectorTemplates)
    haveVectors = present(vectors)
    if ( haveVectors ) haveVectors = size(Vectors) > 0

    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    else if ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    end if

    Time_Config%wait_time = 5 ! Wait just 5s before trying to read inputFile
    clean = .false.
    truncate = .false.
    reset = .false.
    rank = 0
    start = 0
    count = 0
    stride = 0
    block = 0
    startNode = 0
    GotFirst = .false.
    details = 0 - DetailReduction
    OptionsString = '-'
    got = .false.
    oldPrintNameAsHeadline = PrintNameAsHeadline
    oldPrintNameInBanner   = PrintNameInBanner
    toStderr = .false.
    asBanner = .false.
    asHeadLine = .false.
    doStretch = .false.
    doStretchier = .false.
    lineLength = 40
    text = ' '
    mesg = ' '
    Namestring = ' '
    ZOT = .false.

    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      L2CFNODE = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      where_at = where(gson) ! column + 256*line in l2cf & file name
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      ! This first heaped set of fields need no "right-hand side"
      case ( f_allBooleans, f_allFiles, f_allForwardModels, f_allGriddedData, &
        & f_allHessians, f_allHGrids, f_allL2PCs, f_allLines, f_allMatrices, &
        & f_allPFA, f_allQuantityTemplates, &
        & f_allRadiometers, f_allSignals, f_allSpectra, &
        & f_allVectors, f_allVectorTemplates, f_allVGrids, f_antennaPatterns, &
        & f_callStack, f_chunkDivide, f_chunkNumber, f_commandLine, f_crashBurn, &
        & f_DACSfilterShapes, f_filterShapes, f_globalAttributes, f_igrf,&
        &  f_memory, f_MieTables, &
        & f_pause, f_pfaFiles, f_pfaStru, f_phaseName, f_pointingGrids, &
        & f_polygon, f_stop, f_stopWithError, f_time, f_totalMatrixSizes, &
        & f_totalVectorSizes, f_ZOT )
        if ( get_Boolean(son) ) then
          select case ( fieldIndex )
          case ( f_allBooleans )
            call dumpBooleans
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
            if ( details < -3 ) cycle
            if ( haveGriddedData ) then
              call dump ( griddedDataBase, details , options=optionsString )
            else
              call announceError ( son, noGriddedData )
            end if
          case ( f_allHessians )
            if ( details < -2 ) cycle
            if ( haveHessians ) then
              call dump ( HessianDataBase, details, options=optionsString )
            else
              call announceError ( son, 0, 'Unable to dump HessianDB here; empty or absent' )
            end if
          case ( f_allHGrids )
            if ( haveHGrids ) then
              call dump ( hGrids, details, ZOT )
            else
              call announceError ( son, noHGrid )
            end if
          case ( f_allL2PCs )
            call dumpL2PC( L2PCDataBase, details, options=optionsString )
          case ( f_allLines )
            if ( allocated(lines) ) then
              call dump_lines_database
            else
              call announceError ( son, noLines )
            end if
          case ( f_allMatrices )
            if ( details < -2 ) cycle
            if ( haveMatrices ) then
              call dump ( matrixDataBase, details )
            else
              call announceError ( son, 0, 'Unable to dump MatrixDB here; empty or absent' )
            end if
          case ( f_allPFA )
            if ( details < -2 ) cycle
            call Dump_PFADataBase ( details )
          case ( f_allQuantityTemplates )
            if ( haveQuantityTemplatesDB ) then
              call dump ( quantityTemplatesDB )
            else
              call announceError ( son, noQT )
            end if
          case ( f_allRadiometers )
            if ( details < -2 ) cycle
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
              end if
            else
              call announceError ( son, noRadiometers )
            end if
          case ( f_allSignals )
            if ( details < -2 ) cycle
            if ( details > 1 ) then
              call dump_all
            elseif ( associated(signals) ) then
              call dump ( signals, details=details )
            else
              call announceError ( son, noSignals )
            end if
          case ( f_allSpectra )
            if ( details < -2 ) cycle
            call dump ( catalog, details=details )
          case ( f_allVectors )
            if ( details < -2 ) cycle
            if ( haveVectors ) then
              call dump ( vectors, details=details )
            else
              call announceError ( son, noVectors )
            end if
          case ( f_allVectorTemplates )
            if ( details < -2 ) cycle
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
            if ( details > 0) &
              & call MLSMessageCalls ( 'dumpdb' )
          case ( f_chunkDivide )
            call dump ( chunkDivideConfig )
          case ( f_chunkNumber )
            call outputNamedValue ( 'chunk number', L2Options%currentChunkNumber )
          case ( f_commandLine )
            call outputNamedValue ( 'command line', L2Options%command_line )
          case ( f_crashBurn )
            call FinishMLSL2 ( 'ending mlsl2' )
            NEVERCRASH = .false.
            call MLSL2Message( MLSMSG_CRASH, moduleName, &
              & "Program stopped by /crashBurn field on DUMP statement.")
          case ( f_DACSfilterShapes )
            call dump_dacs_filter_database ( son )
          case ( f_filterShapes )
            call dump_filter_shapes_database ( son )
          case ( f_globalAttributes )
            call dumpGlobalAttributes
          case ( f_igrf )
            call dump_gh ( details, optionsString )
          case ( f_memory )
            call SayMemory
            if ( reset ) TotalMemory1 = TotalMemory2
          case ( f_MieTables )
            call dump_Mie ( details )
          case ( f_pause )
            ! Do we wait for a command read from inputFile or one from stdin?
            ! If from an inputFile we will wait patiently for you to create it
            ! If from stdin we will wait for you to enetr "by hand"
            ! (which may not always be possible--thus the inputFile option)
            if ( got( f_inputFile ) ) then
              ! call output ( 'Program paused by /pause field on Dump statement' &
              !   & // ' .. to read command from ' // trim(nameString), &
              !   & advance='yes' )
              ! call Wait_then_read_mesg( nameString, mesg )
              call Pause ( mesg, nameString )
              call output ( 'Resuming after reading ' // trim(mesg), &
                & advance='yes' )
            else
              ! call output ( 'Program paused by /pause field on Dump statement.', &
              !   & advance='yes pause' )
              call Pause ( mesg )
            endif
          case ( f_pfaFiles )
            call dump_PFAFileDatabase ( details )
          case ( f_pfaStru )
            call dump_PFAStructure ( details )
          case ( f_phaseName )
            call outputNamedValue ( 'phase name', L2Options%currentphaseName )
          case ( f_pointingGrids )
            call dump_pointing_grid_database ( son )
          case ( f_polygon )
            call dump_polygon_database
          case ( f_stop )
            call FinishMLSL2 ( 'ending mlsl2' )
            if ( switchDetail(switches, 'time') >= 0 ) then
              call output('(Now for the timings summary)', advance='yes')
              call dump_section_timings
            end if
            if ( NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave ) then
              call get_where ( where_at, farewell, &
                & before="Program stopped with normal status by /stop field on DUMP statement at ", &
                & after="." )
              call MLSMessageExit( NORMAL_EXIT_STATUS, farewell=farewell )
            else if( parallel%slave ) then
              call closeParallel(0)
              call get_where ( where_at, farewell, &
                & before="Slave stopped by /stop field on DUMP statement at ", &
                & after="." )
              call MLSMessageExit( farewell=farewell )
            else
              call get_where ( where_at, farewell, &
                & before="Program stopped by /stop field on DUMP statement at ", &
                & after="." )
              call MLSMessageExit( farewell=farewell )
            end if
          case ( f_stopWithError )
            if ( switchDetail(switches, 'time') >= 0 ) then
              call output('(Now for the timings summary)', advance='yes')
              call dump_section_timings
            end if
            call get_where ( where_at, farewell, &
              & before="Program stopped by /stopWithError field on DUMP statement at ", &
                & after="." )
            call MLSMessageExit( 1, farewell=farewell )
          case ( f_time )
            call SayTime
            if ( reset ) T1 = T2
          case ( f_totalMatrixSizes )
            if ( details < -2 ) cycle
            if ( haveMatrices ) then
              call output ( matricesMemoryInUse ( matrixDataBase ), &
                & before='Total number of elements in all matrices: ', &
                & advance='yes' )
            else
              call announceError ( son, 0, 'Unable to dump MatrixDB size here; empty or absent' )
            end if
          case ( f_totalVectorSizes )
            if ( details < -2 ) cycle
            if ( haveMatrices ) then
              call output ( vectorsMemoryInUse ( vectors ), &
                & before='Total number of elements in all vectors: ', &
                & advance='yes' )
            else
              call announceError ( son, 0, 'Unable to dump Vector DB size here; empty or absent' )
            end if
          case ( f_ZOT )
            ZOT = .true.
          end select
        end if
      case ( f_dumpFile )
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        if ( truncate ) call truncate_textFile( BooleanString )
        call switchOutput ( BooleanString )
      case ( f_truncate )
        truncate = get_Boolean(son)
      case ( f_Boolean )
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        BooleanString = lowerCase(BooleanString)
        ! 1st--check whether we're dumping an array-valued run-time Boolean
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & trim(BooleanString) // 'n', label, countEmpty, &
          & inseparator=runTimeValues%sep )
        if ( label == runTimeValues%sep ) then
          ! Do we s[kip] the key?
          call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
            & BooleanString, label, countEmpty, &
            & inseparator=runTimeValues%sep )
          text = label
          if ( index( 's', optionsString ) < 1 ) &
            text = trim(BooleanString) // ' = ' // label
          call outputNow
        else
          ! OK, we're asked to dump an array-valued one
          call readIntsFromChars ( label, n )
          call output( 'array-valued run-time Boolean ' // trim(BooleanString), &
            & advance='yes' )
          call output ( '  i          value', advance='yes' )
          do i=1, n
            call writeIntsToChars ( i, label )
            call output ( i, advance='no' )
            call blanks ( 4, advance='no' )
            keyString = trim(BooleanString) // '(' // trim(adjustl(label)) // ')'
            call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
              & keyString, label, countEmpty, &
              & inseparator=runTimeValues%sep )
            call output ( trim(label), advance='yes' )
            ! if ( verbose ) call outputnamedValue( 'keyString', trim(keyString) )
            ! -- Don't retain the following--it's just for testing --
            if ( got(f_vector) ) then
              QuantityIndex = GetVectorQuantityIndexByName ( vector, label, noErr=.true. )
              ! call outputNamedValue ( 'QuantityIndex', QuantityIndex )
              if ( QuantityIndex > 0 ) then
                Quantity => GetVectorQuantity( vector, QuantityIndex )
                call dump( Quantity )
              end if
            end if
          enddo
        end if
      case ( f_clean )
        clean = get_Boolean(son)
        if ( clean ) optionsString = trim(optionsString) // 'c'
      case ( f_details )
        call expr ( gson, units, values )
        details = nint(values(1)) - DetailReduction
        ! call outputnamedValue( 'DetailReduction', DetailReduction )
        ! call outputnamedValue( 'Details', details )
      case ( f_rank )
        call expr ( gson, units, values )
        rank = nint(values(1))
      case ( f_reset )
        reset = get_Boolean(son)
      case ( f_start, f_count, f_stride, f_block ) ! For selecting hyperslab
        startNode = subtree(j, root)
        ! Either start = [a, b] or start = b are possible
        hSlabRank = min(nsons(startNode)-1, 3)
        do jj=1, hSlabRank
          call expr(subtree(jj+1,startNode), unitAsArray, valueAsArray)
          mul = valueAsArray(1)
          ! call outputNamedValue( 'mul', mul )
          select case ( fieldIndex )
          case ( f_start )
            start(jj) = mul
            fieldName = 'start'
          case ( f_count )
            count(jj) = mul
            fieldName = 'count'
          case ( f_stride )
            stride(jj) = mul
            fieldName = 'stride'
          case ( f_block )
            block(jj) = mul
            fieldName = 'block'
          end select
        end do
        if ( verbose ) then
          call output('index ', advance='no')
          call output(trim(fieldName), advance='yes')
        end if
      case ( f_file, f_inputfile )
        if ( present(fileDataBase) ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            call get_string ( sub_rosa(gson), nameString, strip=.true. )
            oneMLSFile => getMLSFileByName ( fileDataBase, nameString )
            if ( fieldIndex == f_inputfile ) then
              ! This is just a file to input further commands
            elseif ( associated(oneMLSFile) ) then
              call dumpMLSFile ( oneMLSFile )
            else
              call announceError ( gson, noFile, trim(nameString) )
            end if
          end do
        else
          call announceError ( son, noFileDatabase )
        end if
      case ( f_forwardModel ) ! Dump forward model configs
        if ( details < -2 ) cycle
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
        if ( details < -2 ) cycle
        if ( haveGriddedData ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            if ( gotFirst .and. DiffOrDump == s_diff ) then
              vectorIndex2 = decoration(decoration(gson))
            else
              vectorIndex = decoration(decoration(gson))
            end if
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
            end if
          end do
        else
          call announceError ( gson, noGriddedData )
        end if
        gotFirst = DiffOrDump == s_diff
      case ( f_surface )
        do i = 2, nsons(son)
          call expr ( subtree(i,son), units, values )
          height = values(1)
        end do
      case ( f_hessian ) ! Diff or Dump hessians
        if ( details < -2 ) cycle
        if ( haveHessians ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            ! The decoration is the negative of the index; see Fill, where
            ! the Hessian spec is processed.
            ! Well, that's a nasty trick.
            if ( gotFirst .and. DiffOrDump == s_diff ) then
              hessianIndex2 = -decoration(decoration(gson))
            else
              hessianIndex = -decoration(decoration(gson))
            end if
            if ( DiffOrDump == s_diff ) then
              call diff ( HessianDatabase(hessianIndex), &
                &         HessianDatabase(hessianIndex2), details=details, &
                & options=optionsString )
            else
              call dump ( HessianDatabase(hessianIndex), details=details, &
                & options=optionsString )
            end if
          end do
          gotFirst = DiffOrDump == s_diff
        else
          call announceError ( gson, 0, 'Unable to dump Hessian here; db empty or absent' )
        end if
      case ( f_hGrid )    ! Dump HGrids
        if ( details < -2 ) cycle
        if ( haveHGrids ) then
          do i = 2, nsons(son)
            call output ( ' HGrid ' )
            call dump ( &
              & hGrids(decoration(decoration(subtree(i,son))))%the_hGrid, &
              & details, ZOT )
          end do
        else
          call announceError ( gson, noHGrid )
        end if
      case ( f_L2PC )    ! Dump L2PC
        if ( details < -2 .or. .not. present(FileDataBase) ) cycle
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
        if ( get_Boolean(son) ) call cpu_time ( cpuTimeBase )
      case ( f_mask, f_quantity ) ! Diff or Dump vector quantities
        if ( details < -2 ) cycle
        do i = 2, nsons(son)
          gson = subtree(i,son)
          if ( gotFirst .and. DiffOrDump == s_diff ) then
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
              end if
            end if
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

            if ( clean ) optionsString = trim(optionsString) // 'c'
            ! Special hyperslab handling
            ! We create a new vector quantity qty2 based on qty1
            ! dump that qty2
            ! and then destroy it
            if ( any(got( (/f_start, f_count, f_stride, f_block/) ) ) ) then
              if ( verbose ) then
                call outputNamedValue( 'start(1:rank)', start(1:hSlabRank) )
                call outputNamedValue( 'count(1:rank)', count(1:hSlabRank) )
                call outputNamedValue( 'stride(1:rank)', stride(1:hSlabRank) )
                call outputNamedValue( 'block(1:rank)', block(1:hSlabRank) )
              end if
              VectorValue = GatherVectorQuantity( qty1, &
                & start(1:hSlabRank), count(1:hSlabRank), &
                & stride(1:hSlabRank), block(1:hSlabRank) )
              if ( index(optionsString, '1') < 1 ) &
                & call dump ( VectorValue, details=0, &
                & vector=vectors(vectorIndex) )
              select case ( hSlabRank )
              case (1)
                call dump( VectorValue%value1, 'values as 1-d array' )
              case (2)
                call dump( VectorValue%values, 'values as 2-d array' )
              case (3)
                call dump( VectorValue%value3, 'values as 3-d array' )
              case default
                call MLSL2Message (MLSMSG_Error, moduleName, &
                  & 'rand of hyperslab must be one of 1,2, or 3' )
              end select
              call DestroyVectorQuantityValue( VectorValue, destroyMask=.true., &
                & destroyTemplate=.true. )
            ! Special options handling
            ! 'b' means show bit names
            ! 'c' means clean
            ! '0' means dump template
            ! '1' means dump values
            ! '2' means dump mask
            elseif ( .not. any(indexes(optionsString, (/'0','1','2','b'/)) > 0 ) ) then
              call dump ( qty1, details=details, &
                & vector=vectors(vectorIndex), options=optionsString )
            else if ( index(optionsString, 'b') > 0 ) then
              call FindUnique( int(qty1%value1), uniqueVals, nUnique )
              call outputNamedValue(' Number unique values', nUnique )
              do ii=1, nUnique
                call DumpBitNames( uniqueVals(ii), StatusBitNames )
              enddo
              call DumpBitCounts( int(qty1%value1), showPct=.true. )
            else if ( index(optionsString, '1') > 0 ) then
              if ( verbose ) then
                call outputnamedValue ( 'rank field', rank )
                call outputnamedValue ( 'shape(value3)', shape(qty1%value3) )
              end if
              select case ( rank )
              case ( 1 )
                call dump ( qty1%value1, 'quantity values', &
                  & options=optionsString )
              case ( 3 )
                call dump ( qty1%value3, 'quantity values', &
                  & options=optionsString )
              case ( 4 )
                call dump ( qty1%value4, 'quantity values', &
                  & options=optionsString )
              case default
                call dump ( qty1%values, 'quantity values', &
                  & options=optionsString )
              end select
            else if ( index(optionsString, '0') > 0 ) then
              call dump ( qty1%template, details=details )
            else
              call dumpQuantityMask( qty1, details=0, options=optionsString )
            end if
          end if
        end do
        GotFirst = DiffOrDump == s_diff
      case ( f_matrix ) ! Diff or Dump matrices
        if ( details < -2 ) cycle
        if ( haveMatrices ) then
          do i = 2, nsons(son)
            gson = subtree(i,son)
            if ( gotFirst .and. DiffOrDump == s_diff ) then
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
        GotFirst = DiffOrDump == s_diff
      case ( f_options )
        call get_string ( sub_rosa(gson), optionsString, strip=.true. )
        if ( verbose ) call outputNamedValue( 'options', trim(optionsString) )
        toStderr = optionDetail( optionsString, 'e', 'Err' ) == 'yes'
        asBanner = optionDetail( optionsString, 'B', 'Banner' ) == 'yes'
        headLineChars = optionDetail( optionsString, 'H', 'Headline' )
        asHeadline =  headLineChars /= 'no'
        lineLenChars = optionDetail( optionsString, 'L', 'Linelength' )
        if ( lineLenChars /= 'no' ) &
          & call readIntsFromChars ( lineLenChars, lineLength )
        if ( asHeadline ) then
          PrintNameAsHeadline = .true.
          call removeOption( optionsString, tempOptionsString, '-H' )
          call removeOption( tempOptionsString, optionsString, '--Headline' )
        endif
        if ( asBanner ) then
          PrintNameInBanner = .true.
          call removeOption( optionsString, tempOptionsString, '-B' )
          call removeOption( tempOptionsString, optionsString, '--Banner' )
        endif
      case ( f_pfaData )
        do i = 2, nsons(son)
          look = decoration(decoration(subtree(i,son)))
          call dump ( pfaData(look), details, look )
        end do
      case ( f_pfaNum )
        do i = 2, nsons(son)
          call expr ( subtree(i,son), units, values )
          call dump ( pfaData(nint(values(1))), details, nint(values(1)) )
        end do
      case ( f_signals )
        do i = 2, nsons(son)
          what = decoration(decoration(subtree(i,son)))
          call output ( what, after=': ' )
          call dump ( signals(what), details )
        end do
      case ( f_spectroscopy )
        do i = 2, nsons(son)
          what = decoration(subtree(i,son))
          call output ( what, after=': ' )
          call dump ( catalog(what), details=details )
        end do
      case ( f_stack )
        call dump_stack( where=.true., CPU=.true., doDepth=details>0 )
      case ( f_template ) ! Dump vector templates or quantity templates
        if ( details < -2 ) cycle
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
        ! Special use of options string
        !   option         meaning
        ! ---------        -------
        !     B            print as a banner
        !     H            print as a headine
        !     H[c]         print as a headine with character "c" instead of *
        !     e            print to stderr
        !     s            s-t-r-e-t-c-h
        !     S            s t r e t c h
        !   L[nn]          use nn for line length
        if ( len_trim(optionsString) > 0 ) then
          doStretch = optionDetail( optionsString, 's', 'stretch' ) == 'yes'
          doStretchier = optionDetail( optionsString, 'S', 'stretch' ) == 'yes'
        end if
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
              call get_where ( where(subtree(k,son)), tempText )
              text = text(:i-1) // trim(tempText) // text(i+2:)
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
          if ( doStretchier ) text = stretch( text, options='-a' )
          if ( doStretch ) text = stretch( text )
          call outputNow ! the text
        end do ! k
      case ( f_tGrid )
        if ( details < -2 ) cycle
        do i = 2, nsons(son)
          call output ( ' TGrid ' )
          call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
        end do
      case ( f_variable )
        do i = 2, nsons(son)
          decl = get_decl ( sub_rosa(subtree(i,son)), type=variable )
          if ( decl%type /= variable ) then
            call display_string ( sub_rosa(subtree(i,son)), before='The symbol "' )
            call output ( '" has not been given a value.', advance='yes' )
          else
            call display_string ( sub_rosa(subtree(i,son)) )
            call dump_a_decl ( decl, before=' is', details=details )
          end if
        end do
      case ( f_vectormask, f_vector ) ! Dump entire vectors
        if ( details < -2 ) cycle
        if ( haveVectors ) then
          do i = 2, nsons(son)
            vector => vectors(decoration(decoration(subtree(i,son))))
            call output ( ' Vector ' )
            if ( fieldIndex == f_vectormask ) then
              call dumpVectorMask ( vector, details=details )
            else
              call dump ( vector, details=details, clean=clean, options=optionsString )
            end if
          end do
        else
          call announceError ( gson, noVectors )
        end if
      case ( f_vGrid )
        if ( details < -2 ) cycle
        do i = 2, nsons(son)
          call output ( ' VGrid ' )
          call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
        end do
      end select
    end do
    if ( got(f_dumpFile) ) call revertOutput
    PrintNameAsHeadline = oldPrintNameAsHeadline
    PrintNameInBanner = oldPrintNameInBanner

    call trace_end ( 'DumpCommand', cond=toggle(gen) )
contains
    subroutine Wait_then_read_mesg ( filename, mesg )
      character(len=*), intent(in)  :: filename
      character(len=*), intent(out) :: mesg
      integer                       :: unitnum, status
      ! Executable
      InputFileName = filename
      call wait_for_event( File_There_Now, 1 )
      open ( newunit=unitnum, form='formatted', &
        & file=trim(FileName), status='old', iostat=status )
      read ( unitnum, '(a80)' ) mesg
    end subroutine Wait_then_read_mesg

    function File_There_Now( id ) result ( now )
      ! Returns
      ! 0 if filename exists
      ! 1 otherwise
      integer, intent(in)          :: id ! Ignored
      integer                      :: now
      now = mls_exists ( trim(InputFileName) )
    end function File_There_Now
  end subroutine DumpCommand

  ! ---------------------------------------------  ExecuteCommand  -----
  subroutine ExecuteCommand ( root )
    use Dump_1, only: DumpTextFile
    use Expr_M, only: Expr
    use Init_Tables_Module, only: F_Crashburn, &
      & F_Delay, F_Details, F_Filename, F_Lines, F_Options, &
      & F_Stop, F_Stopwitherror, F_Command, F_Time, &
      & F_Wait, &
      & Field_First, Field_Last
    use IO_Stuff, only: Write_TextFile
    use L2parinfo, only: Parallel, CloseParallel
    use Lexer_Core, only: Get_Where, Where_T
    use Machine, only: Create_Script, Execute, NeverCrash, Usleep
    use MLSL2Options, only: L2CFNode, Normal_Exit_Status  !, MLSMessage
    use MLSL2Timings, only: Dump_Section_Timings
    use MLSMessageModule, only: MLSMSG_Crash, MLSMSG_Warning, &
      & MLSMessageExit
    use MLSStringLists, only: SwitchDetail
    use MoreTree, only: Get_Boolean, Get_Field_ID
    use String_Table, only: Get_String
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: NSons, Sub_Rosa, SubTree, Where

  ! Process an Execute command
    integer, intent(in) :: Root ! Root of the parse tree for the Execute command
    ! Local variables
    character(len=256)         :: command
    integer                    :: fieldIndex
    logical, dimension(field_first:field_last) :: GOT
    integer                    :: Delay
    integer                    :: DetailReduction
    integer                    :: Details
    character(len=255)         :: Farewell
    character(len=80)          :: filename  ! E.g., '/tmp/Execute_me.sh'
    integer                    :: gson
    integer                    :: j
    integer                    :: k
    character(len=128), dimension(64) :: lines
    integer                    :: Me = -1  ! String index for trace cacheing
    integer                    :: nLines
    character(len=80)          :: OPTIONSSTRING  ! E.g., '-rbs' (see dump_0.f90)
    logical :: reset
    integer                    :: son
    integer                    :: status
    logical                    :: verbose
    integer                    :: Units(2) ! of the Details expr -- known to be phyq_dimensionless
    double precision           :: Values(2) ! of the Details expr
    logical                    :: wait
    type(where_t)              :: Where_At
    ! Executable
    call trace_begin ( me, 'ExecuteCommand', root, cond=toggle(gen) )
    verbose = BeVerbose ( 'bool', -1 )
    DetailReduction = switchDetail(switches, 'red')
    if ( DetailReduction < 0 ) then ! The 'red' switch is absent
      DetailReduction = 0
    else if ( DetailReduction == 0 ) then ! By default, reduce details level by 2
      DetailReduction = 2
    end if
    delay = 0
    details = 0 - DetailReduction
    lines = ' '
    command = ' '
    nLines = 0
    OPTIONSSTRING = '-'
    reset = .false.
    got= .false.

    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      L2CFNODE = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      where_at = where(gson) ! column + 256*line in l2cf & file name
      got(fieldIndex) = .true.
      select case ( fieldIndex )
        case ( f_wait )
          wait = get_Boolean(son)
        case ( f_crashBurn )
        case ( f_delay )
          call expr ( gson, units, values )
          delay = nint(values(1))
        case ( f_details )
          call expr ( gson, units, values )
          details = nint(values(1)) - DetailReduction
        case ( f_stop )
        case ( f_stopWithError )
        case ( f_time )
          call SayTime
          if ( reset ) T1 = T2
        case ( f_options )
          call get_string ( sub_rosa(gson), optionsString, strip=.true. )
        case ( f_FileName )
          call get_string ( sub_rosa(gson), fileName, strip=.true. )
        case ( f_command )
          do k = 2, nsons(son)
            call get_string ( sub_rosa(subtree(k,son)), command, strip=.true. )
          enddo
        case ( f_lines )
          do k = 2, nsons(son)
            nLines = nLines + 1
            call get_string ( sub_rosa(subtree(k,son)), lines(nLines), strip=.true. )
          enddo
        end select
      end do
      if ( .not. (got(f_command) .or. got(f_fileName)) ) then
        call MLSL2Message( MLSMSG_Warning, moduleName, &
          & "Execute must have either command or fileName present.")
        return
      endif
      ! Do we have a single-line command or an entire script? Must we wait for it?
      if ( wait ) then
        if ( got(f_lines) ) then
          call Execute_then_wait( inscriptName=fileName, inLines=lines(:nLines) )
        elseif ( got(f_command) ) then
          call Execute_then_wait( command=command )
        endif
      elseif ( got(f_lines) ) then
        call create_script( fileName, lines(:nLines), thenRun=.true., &
          & status=status, delay=delay )
        if ( verbose ) then
          call outputNamedValue( 'script name', trim(fileName) )
          call dumpTextfile( FileName, OptionsString )
        endif
        call execute( 'rm -f ' // trim(fileName), status, delay=delay )
      elseif ( got(f_command) ) then
        call outputNamedValue( 'Command to Execute', trim(command) )
        call execute( command, status, delay=delay )
      else
        call MLSL2Message( MLSMSG_Warning, moduleName, &
          & "Execute command must have either command or lines field.")
      endif
      ! Have we been commanded to stop? wait? crash?
      if ( got(f_stop) ) then
        call FinishMLSL2 ( 'ending mlsl2' )
        call MLSL2Message( MLSMSG_Warning, moduleName, &
          & "Program stopped by /stop field on DUMP statement.")
          if ( switchDetail(switches, 'time') >= 0 ) then
            call output('(Now for the timings summary)', advance='yes')
            call dump_section_timings
          end if
        if ( NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave ) then
          call get_where ( where_at, farewell, &
            & before="Program stopped with normal status by /stop field on DUMP statement at ", &
            & after="." )
          call MLSMessageExit( NORMAL_EXIT_STATUS, farewell=farewell )
        else if( parallel%slave ) then
          call closeParallel(0)
          call get_where ( where_at, farewell, &
            & before="Slave stopped by /stop field on DUMP statement at ", &
            & after="." )
          call MLSMessageExit( farewell=farewell )
        else
          call get_where ( where_at, farewell, &
            & before="Program stopped by /stop field on DUMP statement at ", &
            & after="." )
          call MLSMessageExit( farewell=farewell )
        end if
      elseif ( got(f_stopWithError) ) then
        if ( switchDetail(switches, 'time') >= 0 ) then
          call output('(Now for the timings summary)', advance='yes')
          call dump_section_timings
        end if
        call get_where ( where_at, farewell, &
          & before="Program stopped by /stopWithError field on DUMP statement at ", &
            & after="." )
        call MLSMessageExit( 1, farewell=farewell )
      elseif ( got(f_crashBurn) ) then
        call FinishMLSL2 ( 'ending mlsl2' )
        NEVERCRASH = .false.
        call MLSL2Message( MLSMSG_Crash, moduleName, &
          & "Program stopped by /crashBurn field on DUMP statement.")
      endif
    call trace_end ( 'ExecuteCommand', cond=toggle(gen) )
  contains
    subroutine Execute_then_wait ( command, inscriptName, inlines )
      ! Execute a command or a script
      ! Then wait for it to complete
      ! Method: we create a script, based either on command or on input lines
      ! At the end of the script we rm a sentinel file we cleverly created 
      ! on entering this subroutine
      ! We loop over checking for the sentinel file
      ! When it's gone, so are we
      ! Args
      character(len=*), optional, intent(in)                  :: command
      character(len=*), optional, intent(in)                  :: inscriptName
      character(len=*), dimension(:), optional, intent(in)    :: inlines
      ! Internal variables
      logical                                                 :: exist
      character(len=80), dimension(60)                        :: lines
      character(len=80)                                       :: scriptName
      character(len=80)                                       :: sentinelName
      integer                                                 :: myDelay
      integer                                                 :: nLines
      ! Executable
      scriptName = '/tmp/mlsl2scriptname.sh'
      if ( present(inscriptName) ) scriptName = inscriptName
      sentinelName = trim(scriptName) // '.txt'
      call write_textfile( sentinelName, 'waiting' )
      if ( present(command) ) then
        lines(1) = '#!/bin/sh'
        lines(2) = command
        lines(3) = 'rm -f ' // sentinelName
        nLines = 3
      elseif ( present(inLines) ) then
        nLines = size(inLines) + 1
        lines(1:nLines-1) = inLines
        lines(nLines) = 'rm -f ' // sentinelName
      else
        call MLSL2Message( MLSMSG_Warning, moduleName, &
          & "Execute_then_wait must have either command or inlines present.")
        return
      endif
      call create_script( scriptName, lines(:nLines), thenRun=.true., &
        & status=status, delay=delay )
      if ( verbose ) then
        call outputNamedValue( 'script name', trim(scriptName) )
        call dumpTextfile( scriptName, OptionsString )
      endif
      ! Now we loop over waiting for the sentinel file to disappear
      myDelay = max( delay, 50 )
      do
        call Usleep( myDelay )
        ! Still there?
        inquire( file=sentinelName, exist=exist )
        if ( .not. exist ) exit
      enddo
      call execute( 'rm -f ' // trim(scriptName), status, delay=delay )
    end subroutine Execute_then_wait
  end subroutine ExecuteCommand

  ! ------------------------------------------------  InitializeRepeat  -----
  subroutine  InitializeRepeat
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: PutHashElement
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'count', '0', countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'countsn', '0', countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
  end subroutine  InitializeRepeat

  subroutine  NextRepeat
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: GetHashElement, PutHashElement
    use MLSStrings, only: Readintsfromchars, Writeintstochars
    ! Internal variables
    character(len=64) :: cvalue
    integer :: c
    character(len=16) :: keyString
    character(len=64) :: nvalue
    ! Executable
    call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'count', cvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    call readIntsFromChars ( cvalue, c )
    call writeIntsToChars ( c+1, cvalue )
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'count', cvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    keyString = 'counts(' // trim(adjustl(cvalue)) // ')'
    call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'countsn', nvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    ! Are we Repeating for a sequence of values?
    call readIntsFromChars ( nvalue, c )
    ! call outputNamedValue( 'nvalue', nvalue )
    if ( c < 1 ) return
    ! If so, store the current value at the key 'counts(n)'
    ! call outputNamedValue( 'keyString', keyString )
    call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
      & keyString, cvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )
    call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
      & 'counts(n)', cvalue, countEmpty=countEmpty, &
      & inseparator=runTimeValues%sep )

  end subroutine  NextRepeat

  subroutine  MLSCase ( ROOT )
  ! Returns TRUE if the label or Boolean field of the last Select command
  ! matches the label or Boolean field of the current Case command
  ! or if the current Case command is given the special label 'default',
  ! a wildcard matching anything
    use Init_Tables_Module, only: F_Boolean, F_Label, F_Options
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: GetHashElement
    use MLSStrings, only: Lowercase, Streq
    use Moretree, only: Get_Field_Id
    use String_Table, only: Get_String
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Sub_Rosa, Subtree
    ! Args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: Booleanstring  ! E.g., 'BAND8'
    integer :: FieldIndex
    integer :: GSON, J
    character(len=80) :: Label  ! E.g., 'BAND8'
    integer :: Me = -1          ! String index for trace cacheing
    character(len=16) :: optionsString
    integer :: Son
    logical :: verbose ! , verboser

    ! Executable
    call trace_begin ( me, 'MLSCase', root, cond=toggle(gen) )
    verbose = BeVerbose ( 'bool', -1 )
!    verboser = BeVerbose ( 'bool', 0 )
    MLSSelecting = .true. ! Defaults to skipping rest of case
    if ( MLSSelectedAlready ) then
      call trace_end ( 'MLSCase', cond=toggle(gen) )
      return
    end if
    optionsString = ' '
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        BooleanString = lowerCase(BooleanString)
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & BooleanString, label, countEmpty, &
          & inseparator=runTimeValues%sep )
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
    end if
    if ( verbose ) then
      call output( 'MLSCase label = ', advance='no' )
      call output( trim(label), advance='yes' )
    end if
    ! We must store whether we have ever had a match
    MLSSelectedAlready = MLSSelectedAlready .or. .not. MLSSelecting
    call trace_end ( 'MLSCase', cond=toggle(gen) )
  end subroutine  MLSCase

  subroutine  MLSSelect ( root )
  ! Fills the global variable selectLabel with
  ! the Label field or value of the Boolean field
    use Init_Tables_Module, only: F_Boolean, F_Label
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: GetHashElement
    use MLSStrings, only: Lowercase
    use Moretree, only: Get_Field_Id
    use String_Table, only: Get_String
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Sub_Rosa, Subtree
    ! Args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: Booleanstring  ! E.g., 'BAND8'
    integer :: FieldIndex
    integer :: GSON, J
    character(len=80) :: Label          ! E.g., 'BAND8'
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: Son
    logical :: verbose! , verboser

    ! Executable
    call trace_begin ( me, 'MLSSelect', root, cond=toggle(gen) )
    verbose = BeVerbose ( 'bool', -1 )
!    verboser = BeVerbose ( 'bool', 0 )
    label = ' '
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        BooleanString = lowerCase(BooleanString)
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & BooleanString, label, countEmpty, &
          & inseparator=runTimeValues%sep )
      case (f_label)
        call get_string ( sub_rosa(gson), label, strip=.true. )
      case default
        ! Should not have got here if parser worked correctly
      end select
    enddo
    if ( verbose ) then
      call output( 'Selecting label = ', advance='no' )
      call output( trim(label), advance='yes' )
    end if
    selectLabel = lowerCase(label)
    MLSSelecting = .true.
    call trace_end ( 'MLSSelect', cond=toggle(gen) )
  end subroutine MLSSelect

  subroutine  MLSEndSelect ( ROOT )
  ! Resets the global variable MLSSelecting
  ! Optionally puts note about end of selecting in log file
    use Init_Tables_Module, only: F_Label
    use MLSStrings, only: Lowercase
    use Moretree, only: Get_Field_Id
    use String_Table, only: Get_String
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Sub_Rosa, Subtree
    ! args
    integer, intent(in) :: root
    ! Internal variables
    character(len=80) :: Label  ! E.g., 'BAND8'
    integer :: FieldIndex
    integer :: GSON, J
    integer :: Me = -1          ! String index for trace cacheing
    integer :: Son
    logical :: verbose! , verboser
    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    MLSSelectedAlready = .false.
    call trace_begin ( me, 'MLSEndSelect', root, cond=toggle(gen) )
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
    if ( verbose ) then
      call output( 'End Selecting ', advance='no' )
      if ( len_trim(label) > 0 ) &
        & call output( 'label = ' // trim(label), advance='no' )
      call NewLine
    end if
    call trace_end ( 'MLSEndSelect', cond=toggle(gen) )
  end subroutine MLSEndSelect

  ! ----------- callable as either Skip or Repeat command ---------
  ! Returns value of Boolean field (if present);
  ! If TRUE should skip rest of section in which SKIP command appears
  ! If Boolean field absent, returns TRUE
  ! (otherwise it would do nothing--
  ! this way it forces us to skip unconditionally)
  ! If called as Repeat, returns same result as Skip
  ! Exception: values field, if present, will be
  ! stored as additional runtimes Booleans; e.g.
  ! .., values=["str1", "str2", .., "strn"] adds to r/t Booleans the following
  !    key        value
  ! counts(1)      str1
  ! counts(2)      str2
  !   .    .    .
  ! counts(n)      strn
  ! countsn        n
  logical function Skip ( root, name )
    use Init_Tables_Module, only: F_Boolean, F_Formula, F_NextChunk, F_Values
    use MLSL2Options, only: ExitToNextChunk, RuntimeValues
    use MLSStringLists, only: BooleanValue, GetHashElement, PutHashElement
    use MLSStrings, only: Lowercase, ReadintsFromChars, WriteIntsToChars
    use Moretree, only: Get_Boolean, Get_Field_Id
    use String_Table, only: Get_String
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Sub_Rosa, Subtree
    ! args
    integer, intent(in)                    :: ROOT ! Root of the parse tree
    character(len=*), intent(in), optional :: NAME
    ! Internal variables
    character(len=80) :: Booleanstring  ! E.g., 'BAND8'
    character(len=64) :: cvalue
    integer :: c
    integer :: FieldIndex
    integer :: GSON, J
    character(len=80) :: KEYSTRING
    integer :: Me = -1                  ! String index for trace cacheing
    character(len=16) :: MYNAME
    integer :: NVALUES
    integer :: Son
    integer :: VALUE_FIELD
    logical :: verbose, verboser

    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    verboser = BeVerbose ( 'bool', 0 )
    myName = 'Skip'
    if ( present(name) ) myName = name
    call trace_begin ( me, myName, root, cond=toggle(gen) )
    BooleanString = ' '
    exitToNextChunk = .false. ! Defaults to not skipping rest current chunk
    skip = .true. ! Defaults to skipping rest of section
    value_field = 0
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      gson = son
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        BooleanString = lowerCase(BooleanString)
        if ( verboser ) then
          call output( 'Calling BooleanValue', advance='yes' )
          call output( runTimeValues%lkeys, advance='yes' )
          call output( runTimeValues%lvalues, advance='yes' )
        end if
        skip = BooleanValue ( trim(BooleanString), &
          & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep )
      case ( f_formula )
        call get_string ( sub_rosa(gson), BooleanString, strip=.true. )
        if ( verboser ) call outputNamedValue( 'formula', trim(BooleanString) )
        skip = myBooleanValue ( BooleanString )
      case ( f_nextChunk )
        exitToNextChunk = Get_Boolean ( gson )
      case ( f_values )
        value_field = son
      case default
        ! Should not have got here if parser worked correctly
      end select
      if ( verboser .and. len_trim(BooleanString) > 0 ) then
        call output( trim(BooleanString) // ' = ', advance='no' )
        call output( skip, advance='yes' )
      end if
      if ( skip ) &
        & call output( '(Skipping rest of this section)', advance='yes' )
    end do
    ! The following should only occur if we were called as Repeat
    ! with a values=[...] field
    if ( value_field > 0 ) then
      nvalues = nsons(value_field)-1
      ! Now, have we run out of values to Repeat over?
      call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
        & 'count', cvalue, countEmpty=countEmpty, &
        & inseparator=runTimeValues%sep )
      call readIntsFromChars ( cvalue, c )
      if ( verboser ) then
        call outputnamedValue( 'nvalues', nvalues)
        call outputnamedValue( 'cvalue', trim(cvalue) )
      end if
      if ( c+1 > nvalues ) then
        Skip = .false.
      else if ( c < 1 ) then
        ! The first time through, so try to grok values field
        ! and then store the number of values as 'countsn' and the individual
        ! values as 'counts(1)' 'counts(2)' ..
        call writeIntsToChars( nvalues, BooleanString )
        call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
          & 'countsn', BooleanString, &
          & countEmpty=countEmpty, &
          & inseparator=runTimeValues%sep )
        do j=2, nvalues+1
          call writeIntsToChars( j-1, keyString )
          keyString = 'counts(' // trim(adjustl(keystring)) // ')'
          call get_string( sub_rosa(subtree(j, value_field)), BooleanString, strip=.true. )
          if ( verboser ) then
            call outputnamedValue( 'keyString', trim(keyString) )
            call outputnamedValue( 'BooleanString', trim(BooleanString) )
          end if
          call PutHashElement ( runTimeValues%lkeys, runTimeValues%lvalues, &
            & keyString, BooleanString, &
            & countEmpty=countEmpty, &
            & inseparator=runTimeValues%sep )
        enddo
      end if
    end if
    call trace_end ( myName, cond=toggle(gen) )
  end function Skip

! =====     Private Procedures     =====================================

  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( where, what, string )
    use MoreTree, only: StartErrorMessage

    integer, intent(in) :: What, Where
    character(len=*), intent(in), optional :: String

    call StartErrorMessage ( where )

    select case ( what )
    case ( noFile )
      call output ( "File " // string // " not in database." )
    case ( noFileDatabase )
      call output ( "File database not provided." )
    case ( noFWM )
      call output ( "Can't dump Forward Model Configs here." )
    case ( noGriddedData )
      call output ( "Can't dump GriddedData here, or GriddedDatabase empty. " )
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

  ! --------------------------------------------  BooleanToString  -----
  function BooleanToString ( BOOL ) result ( STR )
    ! Convert a Boolean argument to a character-valued string
    ! I.e., .true. -> 'true'
    !       .false. -> 'false'
    ! Args
    logical, intent(in) :: Bool
    character(len=6)    :: str
    str = merge( 'true ', 'false', Bool )
  end function BooleanToString

  ! ---------------------------------------------  myBooleanValue  -----
  function myBooleanValue ( FORMULA ) result ( BVALUE )
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: BooleanValue, GetStringElement
    use MLSStrings, only: Lowercase, ReadNumsFromChars
    ! Calculate the Boolean value according to
    ! (1) The logical value of its formula, if the formula
    !     does not contain the special operators "==" or "/="
    ! (2) if the formula is "variable == value" and variable is
    !     is recognized and takes the value "value", return "true", else "false"
    ! (3) if the formula is "variable /= value" return not (2)
    ! (4) if the formula is "variable < value" "true" if variable < value
    character(len=*), intent(in) :: formula
    logical :: bvalue
    ! Internal variables
    character(len=16) :: lhs, rhs
    logical, parameter :: countEmpty = .false.
    logical :: reverse  ! Do we mean NOT equal?
    logical :: verbose
    real :: x
    real :: y
    ! Executable
    verbose = BeVerbose ( 'bool', -1 )
    bvalue = .false.
    reverse = .false.
    if ( index(formula, "==") > 1 ) then
      if ( verbose ) call output( 'Have "=="', advance='yes' )
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='=' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='=' )
    elseif ( index(formula, "/=") > 1 ) then
      if ( verbose ) call output( 'Have "/="', advance='yes' )
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='/' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='=' )
      reverse = .true.
    elseif ( index(formula, "<") > 1 ) then
      if ( verbose ) call output( 'Have "<"', advance='yes' )
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='<' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='<' )
      if ( verbose ) then
        call outputnamedvalue('lhs', trim(lhs) )
        call outputnamedvalue('rhs', trim(rhs) )
      end if
      lhs = Evaluator(lowercase(adjustl(lhs)))
      rhs = Evaluator(lowercase(adjustl(rhs)))
      if ( verbose ) then
        call outputnamedvalue('lhs', trim(lhs) )
        call outputnamedvalue('rhs', trim(rhs) )
      end if
      call readNumsFromChars ( lhs, x )
      call readNumsFromChars ( rhs, y )
      bvalue = ( x < y )
      if ( verbose ) then
        call outputNamedValue( 'x', x )
        call outputNamedValue( 'y', y )
      end if
      return
    elseif ( index(formula, ">") > 1 ) then
      if ( verbose ) call output( 'Have ">"', advance='yes' )
      call GetStringElement ( formula, lhs, &
        & 1, countEmpty, inseparator='>' )
      call GetStringElement ( formula, rhs, &
        & 2, countEmpty, inseparator='>' )
      lhs = Evaluator(lowercase(adjustl(lhs)))
      rhs = Evaluator(lowercase(adjustl(rhs)))
      call readNumsFromChars ( lhs, x )
      call readNumsFromChars ( rhs, y )
      if ( verbose ) then
        call outputNamedValue( 'x', x )
        call outputNamedValue( 'y', y )
      end if
      bvalue = ( x > y )
      return
    else
      if ( verbose ) then
        call output( 'Calling BooleanValue', advance='yes' )
        call output( runTimeValues%lkeys, advance='yes' )
        call output( runTimeValues%lvalues, advance='yes' )
      end if
      bvalue = BooleanValue ( trim(formula), &
        & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep )
      if ( verbose ) &
        & call outputnamedvalue( trim(formula) // 'evaluates', bvalue )
      return
    end if
    rhs = lowercase(adjustl(rhs))
    ! OK, so now let's try to make sense of the lhs
    ! Note the necessary type conversions for some variables
    ! that aren't chracter-valued
    ! (Should we have done this for rhs, too?)
    lhs = Evaluator(lowercase(adjustl(lhs)))
    bvalue = (lhs == rhs)
    if ( verbose ) then
      call outputnamedvalue('lhs', trim(lhs) )
      call outputnamedvalue('rhs', trim(rhs) )
    end if
    if ( reverse ) bvalue = .not. bvalue  ! If we meant not equal
  end function myBooleanValue

  ! -----------------------------------------  DumpBooleans  -----
  ! Dump runtime values as table of paired keys and values
  subroutine DumpBooleans
    use MLSL2Options, only: DumpMacros
    call output ( 'Run-time "Booleans"', advance='yes' )
    call dumpMacros( details=1 )
  end subroutine DumpBooleans

  ! -----------------------------------------  EvaluateExplicitly  -----
  function EvaluateExplicitly ( FORMULA ) result( ITSVALUE )
    ! Evaluate all the terms in a formula
    use MLSL2Options, only: RuntimeValues
    use MLSStringLists, only: List2Array, EvaluateFormula, NumStringElements
    ! Args
    character(len=*), intent(in) :: formula
    character(len=MAXRESULTLEN)  :: itsValue
    ! Internal variables
    character(len=64), dimension(1024) :: keys, values
    integer :: n
    ! Executable
    n = NumStringElements( runtimevalues%lkeys, countEmpty, &
      & inseparator=runTimeValues%sep )
    ! call outputNamedValue ( 'n', n )
    call List2Array( runtimevalues%lkeys, keys, countEmpty, &
      & inseparator=runTimeValues%sep, &
      & ignoreLeadingSpaces=.true. )
    call List2Array( runtimevalues%lvalues, values, countEmpty, &
      & inseparator=runTimeValues%sep, &
      & ignoreLeadingSpaces=.true. )
    itsValue = evaluateFormula( formula, values, keys )
  end function EvaluateExplicitly

  ! -------------------------------------------  EvaluateTermwise  -----
  function EvaluateTermwise ( FORMULA ) result( ITSVALUE )
    ! Evaluate all the terms in a formula
    use MLSStringLists, only: Array2List, List2Array, NumStringElements
    ! Args
    character(len=*), intent(in) :: formula
    character(len=MAXRESULTLEN)  :: itsValue
    ! Internal variables
    character(len=MAXRESULTLEN), dimension(24) :: array, array2
    integer :: n
    ! Executable
    n = NumStringElements( formula, countEmpty, inseparator=' ' )
    ! call outputNamedValue ( 'n', n )
    call List2Array( formula, array, countEmpty, inseparator=' ', &
      & ignoreLeadingSpaces=.true. )
    array2 = Evaluator( array )
    call Array2List( array2(1:n), itsValue, inseparator=' ' )
  end function EvaluateTermwise

  ! --------------------------------------------  Evaluator_array  -----
  ! This family of functions evaluates a character-valued arg,
  ! being alert for special values
  ! that name global variables, e.g. 'phasename'
  function Evaluator_array ( ARRAY ) result( VALUES )
    ! Args
    character(len=*), dimension(:), intent(in)           :: array
    character(len=MAXRESULTLEN), dimension(size(array))  :: values
    ! Internal variables
    integer :: i
    ! Executable
    do i=1, size(array)
      values(i) = Evaluator_sca( array(i) )
    enddo
  end function Evaluator_array

  ! ----------------------------------------------  Evaluator_sca  -----
  function Evaluator_sca ( ARG ) result( ITSVALUE )
    use MLSL2Options, only: Checkpaths, Need_L1bFiles, &
      & RuntimeValues, L2Options
    use MLSL2Options, only: L2Options
    use MLSStringLists, only: GetHashElement
    use MLSStrings, only: Lowercase, WriteIntsToChars
    ! Args
    character(len=*), intent(in) :: arg
    character(len=MAXRESULTLEN)  :: itsValue
    ! Executable
    select case (arg)
    case ('checkpaths')
      itsValue = merge( 'true ', 'false', checkpaths )
    case ('chunknumber')
      call writeIntsToChars ( L2Options%currentChunkNumber, itsValue )
    case ('count')
      call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
        & arg, itsValue, countEmpty=countEmpty, &
        & inseparator=runTimeValues%sep )
    case ('phasename')
      itsValue = lowercase(L2Options%currentphaseName)
    case ('need_l1bfiles')
      itsValue = merge( 'true ', 'false', need_l1bfiles )
    case ('skipretrieval')
      itsValue = merge( 'true ', 'false', L2Options%skipretrieval )
    case default
      ! What did you mean?
      ! Maybe just whether two character strings are the same
      ! that were assembled using m4 trickery
      ! But first, let's check that you're not naming a runtime Boolean
      call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
        & arg, itsValue, countEmpty=countEmpty, &
        & inseparator=runTimeValues%sep )
      if ( itsvalue == runTimeValues%sep ) itsValue = arg
    end select
  end function Evaluator_sca

  ! ------------------------------------------------  FinishMLSL2  -----
  subroutine FinishMLSL2 ( Show )
    use Output_M, only: OutputOptions
    use PrintIt_M, only: BothLogUnit, MLSMessageConfig
    use Time_M, only: Finish
    character(len=*), intent(in) :: Show
    ! Print to stdout and, if using Toolkit, log it, too
    if ( MLSMessageConfig%UseToolkit )  then
      OutputOptions%PrUnit = BothLogUnit
      MLSMessageConfig%LogFileUnit = BothLogUnit
    endif
    call Finish ( Show )
  end subroutine FinishMLSL2

  ! ------------------------------------------------  insertSpace  -----
  subroutine insertSpace
    ! We actually insert a space followed by a null; later the null will
    ! be removed before printing
    Text = trim(Text) // ' ' // achar(0)
  end subroutine insertSpace

  ! ------------------------------------------  outputLater_chars  -----
  subroutine outputLater_chars ( chars, advance )
    character(len=*), intent(in) :: chars
    character(len=*), optional, intent(in) :: advance
    text = trim(text) // chars
  end subroutine outputLater_chars

  ! --------------------------------------------  outputLater_int  -----
  subroutine outputLater_int ( arg, advance )
    integer, intent(in)             :: arg
    character(len=*), optional, intent(in) :: advance
    text = trim(text) // numToChars ( arg )
  end subroutine outputLater_int

  ! -------------------------------------------  outputLater_real  -----
  subroutine outputLater_real ( arg, advance )
    real, intent(in)             :: arg
    character(len=*), optional, intent(in) :: advance
    text = trim(text) // numToChars ( arg )
  end subroutine outputLater_real

  ! --------------------------------------------------  outputNow  -----
  subroutine outputNow
    ! Snip away the null chars before printing
    use MLSStrings, only: Asciify
    text = asciify( text, how='snip' )
    if ( asHeadLine ) then
      if ( headLineChars == 'yes' ) then
        call headLine( trim(text), Before='*', After='*', fillChar='-' )
      else
        call headLine( trim(text), &
          & Before=headLineChars(1:1), After=headLineChars(1:1), &
          & fillChar='-' )
      end if
    elseif ( asBanner ) then
      if ( lineLength > 0 ) then
        call banner( trim(asciify( text, how='snip' )), linelength=lineLength )
      else
        call banner( trim(asciify( text, how='snip' )) )
      end if
    elseif ( toStderr ) then
      call Beep ( trim(text) )
    else
      call output ( trim(text), advance='yes' )
    end if
  end subroutine outputNow

  ! -----------------------------------------  returnFullFileName  -----
  subroutine returnfullfilename ( shortname, fullname, &
    & pcf_start, pcf_end )
    use MLSFiles, only: GetPCFromRef
    use MLSL2Options, only: Toolkit
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

  ! --------------------------------------------------  SayMemory  -----
  subroutine SayMemory
    use Memory_M, only: Memory_Used
    text = ' '
    call memory_used ( total=TotalMemory2 )
    call outputLater ( "Memory used since start = " )
    call insertSpace
    call outputLater ( TotalMemory2, advance = 'no' )
    call insertSpace
    if ( TotalMemory1 > 0. ) then
      call outputLater ( "(recent = ", advance='no' )
      call insertSpace
      call outputLater ( TotalMemory2 - TotalMemory1, advance = 'no' )
      call insertSpace
      call outputLater ( ")", advance='no' )
    end if
    call outputNow
  end subroutine SayMemory

  ! ----------------------------------------------------  SayTime  -----
  subroutine SayTime
    use Time_M, only: Time_Now
    text = ' '
    call time_now ( t2 )
    call outputLater ( "Time since start = " )
    call insertSpace
    call outputLater ( t2, advance = 'no' )
    call insertSpace
    if ( T1 > 0. ) then
      call outputLater ( "(lap = ", advance='no' )
      call insertSpace
      call outputLater ( t2 - t1, advance = 'no' )
      call insertSpace
      call outputLater ( ")", advance='no' )
    end if
    call outputNow
  end subroutine SayTime

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
! Revision 2.150  2020/07/22 22:53:04  pwagner
! Added BooleanFromChunkStartsAfter, ..BeginsBefore functions
!
! Revision 2.149  2020/04/30 23:33:42  pwagner
! Skip Dumping AllGriddedData only if Details lt -3
!
! Revision 2.148  2020/01/27 21:34:18  pwagner
! Consistent with new api for Dump_Signals
!
! Revision 2.147  2019/03/18 22:09:34  pwagner
! Ensure finishing messages are both printed and logged
!
! Revision 2.146  2018/12/11 01:26:13  pwagner
! Uses Wait_M
!
! Revision 2.145  2018/11/01 23:19:13  pwagner
! Use Pause from Time_M
!
! Revision 2.144  2018/10/19 00:03:46  pwagner
! inputFile= with /pause pauses to read from inputFile instead of stdin
!
! Revision 2.143  2018/10/17 23:05:04  pwagner
! Dump command can take /pause field to wait for user input for, e.g. debugging
!
! Revision 2.142  2018/09/13 20:22:32  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.141  2018/09/07 00:00:15  pwagner
! More commands that set a runtime flag can now take /reverse
!
! Revision 2.140  2018/08/04 02:10:27  vsnyder
! Make Lines database allocatable instead of a pointer
!
! Revision 2.139  2018/07/27 23:18:48  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.138  2018/05/17 17:05:53  pwagner
! May Dump quantity values bitwise
!
! Revision 2.137  2018/04/19 23:44:09  pwagner
! Skip may take /nextChunk flag
!
! Revision 2.136  2018/03/22 18:17:00  pwagner
! Added command IsFileAbsent; may occur in ReadApriori, MergeGrids, and Output sections
!
! Revision 2.135  2018/03/05 19:49:38  pwagner
! Consistent with details arg to Dump_Signals now being an int
!
! Revision 2.134  2017/09/26 16:24:34  pwagner
! May use options field to set text as Headline or Banner
!
! Revision 2.133  2017/07/27 16:55:02  pwagner
! Reevaluate may take /truncate to keep integer part only
!
! Revision 2.132  2017/07/10 23:09:51  pwagner
! Print less if not verbose
!
! Revision 2.131  2017/02/22 18:01:56  pwagner
! Declaring a runtime Boolean without an initial value no longer generates a Warning
!
! Revision 2.130  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.129  2016/05/18 01:37:30  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.128  2016/04/07 23:40:28  pwagner
! Should exit, not crash, if stopped by /stop flag to Execute
!
! Revision 2.127  2016/04/01 00:27:41  pwagner
! May now Execute a single command or a script of lines from l2cf
!
! Revision 2.126  2016/02/26 02:08:18  vsnyder
! Add ZOT switch
!
! Revision 2.125  2016/01/29 01:10:07  vsnyder
! Add polygon to Dump
!
! Revision 2.124  2015/10/06 00:23:48  pwagner
! When Dumping runtime Boolean, options='s' skips printing key
!
! Revision 2.123  2015/09/24 22:10:17  pwagner
! Dump command may optionally direct text to stderr
!
! Revision 2.122  2015/07/14 23:32:47  pwagner
! may divert Dump commands to named dumpFile; /truncate field
!
! Revision 2.121  2015/04/09 21:00:32  pwagner
! rank may be specified for Dump command; restored printing Dumped text field
!
! Revision 2.120  2015/02/27 23:13:52  pwagner
! May Dump global attributes
!
! Revision 2.119  2015/02/13 00:21:17  pwagner
! Nay dump Booleans more nicely as a Table
!
! Revision 2.118  2014/09/05 00:48:13  vsnyder
! Get field value correctly.  Remove USE for unreferenced names.  Add dumps
! for MatricesMemoryInUse, TotalMatrixSizes, TotalVectorSizes,
! VectorsMemoryInUse.  Some cannonball polishing.  New interface to
! Memory_Used.  Move some USE statements from module scope to procedure
! scope.  Add outputLaterInt.
!
! Revision 2.117  2014/08/28 19:07:16  pwagner
! May command Dump of /memory or /time; optionally to /reset
!
! Revision 2.116  2014/06/11 20:06:33  pwagner
! details level to skip all Dumps lowered to -3
!
! Revision 2.115  2014/05/31 00:24:55  pwagner
! Small corrections; details should now be -3 to skipAllDumps
!
! Revision 2.114  2014/04/25 18:53:45  pwagner
! options string to Dump text as headLine or Banner
!
! Revision 2.113  2014/04/10 02:02:22  vsnyder
! Improve message for nonexistent variable
!
! Revision 2.112  2014/04/10 00:45:41  pwagner
! Moved currentChunkNumber, currentPhaseName from MLSL2Timings to MLSL2Options
!
! Revision 2.111  2014/03/20 01:32:11  vsnyder
! Remove unreferenced USE name
!
! Revision 2.110  2014/02/28 01:10:44  vsnyder
! Remove unused names
!
! Revision 2.109  2014/02/28 00:17:05  vsnyder
! Don't re-check types already checked by type checker.  Move units checking
! to type checker.  Don't handle Height field, since there isn't one allowed
! in init_tables_module.
!
! Revision 2.108  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.107  2013/12/12 02:09:29  vsnyder
! Add 'details' to dump_a_decl
!
! Revision 2.106  2013/11/18 22:33:28  pwagner
! /callStack dumps verbose and debug dbs
!
! Revision 2.105  2013/11/07 00:26:44  pwagner
! May dump all dbs related to signals if details >1; use beVerbose
!
! Revision 2.104  2013/10/24 21:12:19  pwagner
! Corrected some bugs in dumping hyperslabs
!
! Revision 2.103  2013/10/17 18:27:29  pwagner
! Must call trace_end in MLSSelect no matter what
!
! Revision 2.102  2013/10/09 23:44:13  vsnyder
! Add Variable field
!
! Revision 2.101  2013/10/08 23:52:48  pwagner
! Fixed call to BooleanVALUE--must use trim
!
! Revision 2.100  2013/09/27 00:38:12  pwagner
! May select hyperslab when dumping quantity
!
! Revision 2.99  2013/09/25 18:46:50  pwagner
! Fixed bug introduced with last commit
!
! Revision 2.98  2013/09/25 16:33:27  pwagner
! Uses hyperslab fields when dumping vector.quantities
!
! Revision 2.97  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.96  2013/09/24 00:56:02  pwagner
! May specify hyperslab when dumping
!
! Revision 2.95  2013/09/21 00:38:11  vsnyder
! Add ChunkDivide
!
! Revision 2.94  2013/09/04 16:32:51  pwagner
! Replaced '--cat' cmdline option; 'Catenate' now an Output section command
!
! Revision 2.93  2013/08/31 02:29:12  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.92  2013/08/30 23:19:28  pwagner
! Fixed problem of empty spaces while tracing 'Skip'
!
! Revision 2.91  2013/08/29 19:37:41  pwagner
! May Dump plain stack, too
!
! Revision 2.90  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.89  2013/08/01 20:46:30  vsnyder
! Fix to dump two quantities with two quantity fields, two grids with two
! grid fields, two hessians with two hessian fields, two masks with two
! mask fields, or two matrices with two matrix fields.  Had been assuming
! that two fields meant doing a diff.
!
! Revision 2.88  2013/07/12 23:24:11  vsnyder
! Announce error for 'expr' field -- no code yet
!
! Revision 2.87  2013/05/22 20:20:49  pwagner
! Moved insertHashElement to MLSStringLists
!
! Revision 2.86  2013/05/17 00:55:47  pwagner
! May now do elementwise operations on r/t macros with repeat command; fixed many bugs
!
! Revision 2.85  2013/05/07 22:01:30  pwagner
! run-time Booleans can store arrays of values, formulas can evaluate named terms
!
! Revision 2.84  2013/04/24 00:36:02  pwagner
! Added inputBoolean and made Reevaluate formulas more powerful
!
! Revision 2.83  2013/04/22 17:49:20  pwagner
! Reevaluate may store a literal instead of a Boolean value
!
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
! Print notice of changed runtime Booleans only when "bool" switch is set
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
! Dump commands take Boolean fields /stop, /stopWithError, /crashBurn
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
