! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module RetrievalModule
!=============================================================================

!{This module inverts the radiative transfer equation, to solve for
! atmospheric parameters, given radiance measurements.
!
! This module and ones it calls consume most of the cycles.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use BitStuff, only: CountBits
  use Dump_0, only: Dump
  use Expr_M, only: Expr
  use ForwardModelConfig, only: ForwardModelConfig_T
  use Init_Tables_Module, only: F_apriori, F_aprioriScale, F_channels, &
    & F_columnScale, F_Comment, F_covariance, F_diagnostics, F_diagonal, &
    & F_forwardModel, F_fuzz, F_fwdModelExtra, F_fwdModelOut, F_height, &
    & F_ignore, F_jacobian, F_lambda, F_Level, F_maxF, F_maxJ, F_measurements, &
    & F_measurementSD, F_method, F_opticalDepth, F_outputCovariance, &
    & F_outputSD, F_ptanQuantity, F_quantity, F_regOrders, F_regQuants, &
    & F_regWeight, F_state, F_toleranceA, F_toleranceF, F_toleranceR, &
    & Field_first, Field_last, &
    & L_apriori, L_covariance, &
    & L_dnwt_ajn,  L_dnwt_axmax,  L_dnwt_cait, &
    & L_dnwt_diag,  L_dnwt_dxdx,  L_dnwt_dxdxl, &
    & L_dnwt_dxn,  L_dnwt_dxnl,  L_dnwt_flag, L_dnwt_fnmin, &
    & L_dnwt_fnorm,  L_dnwt_gdx,  L_dnwt_gfac, &
    & L_dnwt_gradn,  L_dnwt_sq,  L_dnwt_sq,  L_dnwt_sqt,&
    & L_Jacobian_Cols, L_Jacobian_Rows, L_lowcloud, &
    & L_newtonian, L_none, L_norm, L_pressure, L_zeta, &
    & S_dumpBlocks, S_forwardModel, S_matrix, S_retrieve, S_sids, S_snoop, &
    & S_subset, S_time
  use Intrinsic, only: PHYQ_Dimensionless
  use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
    & DestroyMatrix, GetFromMatrixDatabase, Matrix_T, Matrix_Database_T, &
    & Matrix_SPD_T, MultiplyMatrixVectorNoT
  use MatrixTools, only: DumpBlock
  use MLSCommon, only: R8, MLSCHUNK_T
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, add_to_retrieval_timing
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
  use OUTPUT_M, only: BLANKS, OUTPUT
  use SidsModule, only: SIDS
  use SnoopMLSL2, only: SNOOP
  use String_Table, only: DISPLAY_STRING, Get_String
  use Toggles, only: Gen, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: AddToVector, AddVectorToDatabase, ClearVector, &
    & CloneVector, CopyVector, DestroyVectorInfo, DestroyVectorDatabase, &
    & DumpMask, GetVectorQuantityByType, IsVectorQtyMasked, &
    & Vector_T, VectorValue_T

  implicit NONE
  private
  public :: RETRIEVE

  real(r8), parameter, private :: DefaultInitLambda = 0.0_r8
  integer, parameter, private :: DefaultMaxF = 30, DefaultMaxJ = 5
  integer, parameter, private :: DefaultMethod = l_newtonian
  double precision, parameter, private :: DefaultToleranceA = 1.0d-6 ! for NWT
  double precision, parameter, private :: DefaultToleranceF = 1.0d-6 ! for NWT
  double precision, parameter, private :: DefaultToleranceR = 1.0d-6 ! for NWT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  Retrieve  -----
  subroutine Retrieve ( Root, VectorDatabase, MatrixDatabase, ConfigDatabase, &
    & chunk )

  !{Process the "Retrieve" section of the L2 Configuration File.
  ! The "Retrieve" section can have ForwardModel, Matrix, Sids, Subset or
  ! Retrieve specifications.

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
                                        ! Indexes an n_cf vertex
    type(vector_T), dimension(:), pointer :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: ConfigDatabase

    type(MLSChunk_T), intent(in) :: CHUNK

    ! Local variables:
    type(vector_T), pointer :: Apriori  ! A priori estimate of state
    real(r8) :: AprioriScale            ! Weight for apriori, default 1.0
    integer :: ColumnScaling            ! one of l_apriori, l_covariance,
                                        ! l_none or l_norm
    integer, pointer, dimension(:) :: ConfigIndices    ! In ConfigDatabase
    type(matrix_SPD_T), pointer :: Covariance     ! covariance**(-1) of Apriori
    type(vector_T), pointer :: Diagnostics   ! Diagnostic stuff about the
                                        ! retrieval.
    logical :: Diagonal                 ! "Iterate with the diagonal of the
                                        ! a priori covariance matrix until
                                        ! convergence, then put in the whole
                                        ! thing and iterate until it converges
                                        ! again (hopefully only once).
    integer :: Error
    integer :: Field                    ! Field index -- f_something
    real(r8) :: Fuzz                    ! For testing only.  Amount of "fuzz"
                                        ! to add to state vector before
                                        ! starting a retrieval.
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I, J, K                  ! Subscripts and loop inductors
    real(r8) :: InitLambda              ! Initial Levenberg-Marquardt parameter
    integer :: IxCovariance             ! Index in tree of outputCovariance
    integer :: IxJacobian               ! Index in tree of jacobian matrix
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Jacobian_Cols            ! Number of columns of the Jacobian.
    integer :: Jacobian_Rows            ! (Number of rows of Jacobian) -
                                        ! (masked-off rows of Jacobian)
    integer :: Key                      ! Index of an n_spec_args.  Either
                                        ! a son or grandson of root.
    integer :: MaxFunctions             ! Maximum number of function
                                        ! evaluations of Newtonian method
    integer :: MaxJacobians             ! Maximum number of Jacobian
                                        ! evaluations of Newtonian method
    type(vector_T), pointer :: Measurements  ! The measurements vector
    type(vector_T), pointer :: MeasurementSD ! The measurements vector's Std. Dev.
    integer :: Method                   ! Method to use for inversion, currently
                                        ! only l_Newtonian.
    type(matrix_SPD_T), target :: MyCovariance    ! for OutputCovariance to point at
    type(matrix_T), target :: MyJacobian          ! for Jacobian to point at
    type(vector_T), dimension(:), pointer :: MyVectors ! database
    type(matrix_SPD_T), pointer :: OutputCovariance    ! Covariance of the sol'n
    type(vector_T), pointer :: OutputSD ! Vector containing SD of result
    integer :: RegOrders                ! Regularization orders
    integer :: RegQuants                ! Regularization quantities
    real(r8) :: RegWeight               ! Weight of regularization conditions
    integer :: Son                      ! Of Root or Key
    character(len=127) :: SnoopComment  ! From comment= field of S_Snoop spec.
    integer :: SnoopKey                 ! Tree point of S_Snoop spec.
    integer :: SnoopLevel               ! From level field of S_Snoop spec.
    integer :: Spec                     ! s_matrix, s_subset or s_retrieve
    type(vector_T), pointer :: State    ! The state vector
    real :: T1, T2                      ! for timing
    logical :: Timing
    double precision :: ToleranceA      ! convergence tolerance for NWT,
                                        ! norm of move
    double precision :: ToleranceF      ! convergence tolerance for NWT,
                                        ! norm of F
    double precision :: ToleranceR      ! convergence tolerance for NWT,
                                        ! (norm of move) / (norm of X)
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    logical :: Update                   ! "We are updating normal equations"
    double precision :: Value(2)        ! Value returned by EXPR

    ! Indexes in the private vectors database
    integer, parameter :: FirstVec = 1
    integer, parameter :: AprioriMinusX = firstVec ! Apriori - X
    integer, parameter :: ATb = aprioriMinusX + 1  ! A^T b -- the RHS of the normal eqns    integer, parameter :: 
    integer, parameter :: BestGradient = aTb + 1   ! for NWT
    integer, parameter :: BestX = bestGradient + 1 ! for NWT
    integer, parameter :: CandidateDX = bestX + 1  ! for NWT
    integer, parameter :: ColumnScaleVector = candidateDX + 1 ! For column scaling by column norms
    integer, parameter :: CovarianceDiag = columnScaleVector + 1
    integer, parameter :: CovarianceXApriori = covarianceDiag + 1
    integer, parameter :: DX = covarianceXApriori + 1  ! for NWT
    integer, parameter :: DXUnScaled = dX + 1      ! for NWT
    integer, parameter :: F = dXUnscaled + 1       ! Residual -- Model - Measurements
    integer, parameter :: F_rowScaled = f + 1      ! Either a copy of f, or f row-scaled
    integer, parameter :: Gradient = f_rowScaled + 1   ! for NWT
    integer, parameter :: Reg_X_x = gradient + 1   ! Regularization * X_n
    integer, parameter :: Weight = reg_X_x + 1     ! Scaling vector for rows, 1/measurementSD
    integer, parameter :: X = weight + 1  ! for NWT
    integer, parameter :: LastVec = X
 
    type(vector_T), dimension(firstVec:lastVec) :: V   ! Database for snoop

    ! Error message codes
    integer, parameter :: BothOrNeither = 1       ! Only one of two required
                                                  !    fields supplied
    integer, parameter :: IfAThenB = BothOrNeither + 1
    integer, parameter :: IfUnitsAThenB = IfAThenB + 1
    integer, parameter :: Inconsistent = IfUnitsAThenB + 1  ! Inconsistent fields
    integer, parameter :: InconsistentUnits = Inconsistent + 1
    integer, parameter :: NoFields = InconsistentUnits + 1  ! No fields are allowed
    integer, parameter :: NotSPD = noFields + 1   ! Not symmetric pos. definite
    integer, parameter :: OrderAndWeight = notSPD + 1  ! Need both or neither
    integer, parameter :: WrongUnits = OrderAndWeight + 1

    error = 0
    nullify ( apriori, configIndices, covariance, fwdModelOut )
    nullify ( measurements, measurementSD, myVectors, state, outputSD )
    snoopComment = ' '
    snoopKey = 0
    snoopLevel = 1
    timing = section_times
    do j = firstVec, lastVec ! Make the vectors in the database initially empty
      nullify ( v(j)%quantities, v(j)%template%quantities )
      v(j)%name = 0 ! so Snoop won't use it
    end do
    if ( timing ) call cpu_time ( t1 )

    if ( toggle(gen) ) call trace_begin ( "Retrieve", root )
    do i = 2, nsons(root) - 1           ! skip names at begin/end of section
      son = subtree(i, root)
      if ( node_id(son) == n_named ) then
        key = subtree(2, son)
      else
        key = son
      end if

      ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
      ! parser users' guide" for pictures of the trees being analyzed.

      got = .false.
      spec = get_spec_id(key)
      select case ( spec )
      case ( s_dumpblocks )
        call DumpBlock ( key, matrixDatabase )
      case ( s_matrix )
        if ( toggle(gen) ) call trace_begin ( "Retrieve.matrix/vector", root )
        if ( nsons(key) /= 1 ) call announceError ( noFields, spec )
        call destroyMatrix( matrixDatabase(decoration(key)) ) ! avoids a memory leak
        call decorate ( key, 0 )
        if ( toggle(gen) ) call trace_end ( "Retrieve.matrix/vector" )
      case ( s_snoop )
        snoopKey = key
        do j = 2, nsons(key)
          son = subtree(j, key)
          field = get_field_id(son)  ! tree_checker prevents duplicates
          select case ( field )
          case ( f_comment )
            call get_string ( sub_rosa(subtree(2,son)), snoopComment, strip=.true. )
          case ( f_level )
            call expr ( subtree(2,son), units, value, type )
            if ( units(1) /= phyq_dimensionless ) &
              & call announceError ( wrongUnits, field, string='no' )
            snoopLevel = nint(value(1))
          end select
        end do
      case ( s_subset )
        if ( toggle(gen) ) call trace_begin ( "Retrieve.subset", root )
        call SetupSubset ( key, vectorDatabase )
        if ( toggle(gen) ) call trace_end ( "Retrieve.subset" )
      case ( s_retrieve )
        if ( toggle(gen) ) call trace_begin ( "Retrieve.retrieve", root )
        aprioriScale = 1.0
        columnScaling = l_none
        diagonal = .false.
        initLambda = defaultInitLambda
        maxFunctions = defaultMaxF
        maxJacobians = defaultMaxJ
        method = defaultMethod
        regQuants = 0
        toleranceA = defaultToleranceA
        toleranceF = defaultToleranceF
        toleranceR = defaultToleranceR
        do j = 2, nsons(key) ! fields of the "retrieve" specification
          son = subtree(j, key)
          field = get_field_id(son)  ! tree_checker prevents duplicates
          got(field) = .true.
          select case ( field )
          case ( f_apriori )
            apriori => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_columnScale )
            columnScaling = decoration(subtree(2,son))
          case ( f_covariance )      ! of apriori
            call getFromMatrixDatabase ( &
              & matrixDatabase(decoration(decoration(subtree(2,son)))), &
              & covariance)
            if ( .not. associated(covariance) ) &
              & call announceError ( notSPD, field )
          case ( f_diagnostics )
            diagnostics => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_diagonal )
            diagonal = get_Boolean(son)
          case ( f_forwardModel )
            call allocate_test ( configIndices, nsons(son)-1, "ConfigIndices", &
              & moduleName )
            do k = 2, nsons(son)
              configIndices(k-1) = decoration(decoration(subtree(k,son)))
            end do
          case ( f_fwdModelExtra )
            fwdModelExtra => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_fwdModelOut )
            fwdModelOut => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_jacobian )
            ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
          case ( f_measurements )
            measurements => vectorDatabase(decoration(decoration(subtree(2,son))))
            call cloneVector ( v(f), measurements, vectorNameText='_f' )
          case ( f_measurementSD )
            measurementSD => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_method )
            method = decoration(subtree(2,son))
          case ( f_regOrders )
            regOrders = son
          case ( f_regQuants )
            regQuants = son
          case ( f_outputCovariance )
            ixCovariance = decoration(subtree(2,son)) ! outCov: matrix vertex
          case ( f_outputSD )
            outputSD => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_state )
            state => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_aprioriScale, f_fuzz, f_lambda, f_maxF, f_maxJ, &
            &    f_regWeight, f_toleranceA, f_toleranceF, f_toleranceR )
            call expr ( subtree(2,son), units, value, type )
            if ( units(1) /= phyq_dimensionless ) &
              & call announceError ( wrongUnits, field, string='no' )
            select case ( field )
            case ( f_aprioriScale )
              aprioriScale = value(1)
            case ( f_fuzz )
              fuzz = value(1)
            case ( f_lambda )
              initLambda = value(1)
            case ( f_maxF )
              maxFunctions = value(1)
            case ( f_maxJ )
              maxJacobians = value(1)
            case ( f_regWeight )
              regWeight = value(1)
            case ( f_toleranceA )
              toleranceA = value(1)
            case ( f_toleranceF )
              toleranceF = value(1)
            case ( f_toleranceR )
              toleranceR = value(1)
            end select
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)

        if ( got(f_apriori) .neqv. got(f_covariance) ) &
          & call announceError ( bothOrNeither, f_apriori, f_covariance )
        if ( got(f_regOrders) .neqv. got(f_regWeight) ) &
          & call announceError ( bothOrNeither, f_regOrders, f_regWeight )
        if ( got(f_regQuants) .and. .not. got(f_regOrders) ) &
          & call announceError ( ifAThenB, f_regQuants, f_regOrders )
        if ( error == 0 ) then

          ! Verify the consistency of various matrices and vectors
          if ( got(f_apriori) ) then
            if ( apriori%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_apriori, f_state )
          end if
          if ( associated(covariance) ) then
            if ( covariance%m%row%vec%template%id /= state%template%id .or. &
              &  covariance%m%col%vec%template%id /= state%template%id ) &
              &  call announceError ( inconsistent, f_covariance, f_state )
          end if
          if ( got(f_measurementSD) ) then
            if ( measurementSD%template%id /= measurements%template%id ) &
              & call announceError ( inconsistent, f_measurementSD, f_measurements )
          end if
        end if
        if ( error == 0 ) then

          ! Create the Jacobian matrix
          if ( got(f_jacobian) ) then
            k = decoration(ixJacobian)
            if ( k == 0 ) then
              call createEmptyMatrix ( myJacobian, &
                & sub_rosa(subtree(1,ixJacobian)), measurements, state )
              k = addToMatrixDatabase( matrixDatabase, myJacobian )
              call decorate ( ixJacobian, k )
            end if
            call getFromMatrixDatabase ( matrixDatabase(k), jacobian )
            if ( jacobian%row%vec%template%id /= measurements%template%id ) &
              & call announceError ( inconsistent, f_jacobian, f_measurements )
            if ( jacobian%col%vec%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_jacobian, f_state )
          else
            jacobian => myJacobian
            call createEmptyMatrix ( jacobian, 0, measurements, state )
          end if
          ! Create the output covariance matrix
          if ( got(f_outputCovariance) ) then
            k = decoration(ixCovariance)
            if ( k == 0 ) then
              call createEmptyMatrix ( myCovariance%m, &
                & sub_rosa(subtree(1,ixCovariance)), state, state )
              k = addToMatrixDatabase( matrixDatabase, myCovariance )
              call decorate ( ixCovariance, k )
            end if
            call getFromMatrixDatabase ( matrixDatabase(k), outputCovariance )
            if ( .not. associated(outputCovariance) ) then
              call announceError ( notSPD, f_outputCovariance )
            else
              if ( outputCovariance%m%row%vec%template%id /= state%template%id &
                & .or. &
                &  outputCovariance%m%col%vec%template%id /= state%template%id ) &
                & call announceError ( inconsistent, f_outputCovariance, f_state )
            end if
          else
            outputCovariance => myCovariance
            call createEmptyMatrix ( myCovariance%m, 0, state, state )
          end if
        end if
        if ( error == 0 ) then
          ! Do the retrieval
          jacobian_Cols = 0
          jacobian_Rows = 0
          select case ( method )
          case ( l_newtonian )
        ! call add_to_retrieval_timing( 'newton_solver', t1 )
            call newtonianSolver
          case ( l_lowcloud )
        !    call add_to_retrieval_timing( 'low_cloud', t1 )
        !    call cpu_time ( t1 )
            call LowCloudRetrieval
          end select ! method
          !??? Make sure the jacobian and outputCovariance get destroyed
          !??? after ?what? happens?  Can we destroy the entire matrix
          !??? database at the end of each chunk?
          if ( .not. got(f_jacobian) ) call destroyMatrix ( jacobian )
          if ( .not. got(f_outputCovariance) ) &
            & call destroyMatrix ( outputCovariance%m )
        end if
        if ( got(f_fwdModelOut) ) then
          call copyVector ( fwdModelOut, v(f) )
          call addToVector ( fwdModelOut, measurements )
        end if
        call deallocate_test ( configIndices, "ConfigIndices", moduleName )
        if ( toggle(gen) ) call trace_end ( "Retrieve.retrieve" )
      case ( s_sids )
        ! call add_to_retrieval_timing( 'sids', t1 )
        call cpu_time ( t1 )
        call sids ( key, VectorDatabase, MatrixDatabase, configDatabase, chunk)
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select

      do j = firstVec, lastVec
        call destroyVectorInfo ( v(j) )
      end do

    end do ! i = 2, nsons(root) - 1
    if ( toggle(gen) ) call trace_end ( "Retrieve" )
    if ( timing ) call sayTime

  contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex, AnotherFieldIndex, String )

      use Intrinsic, only: Field_indices, Spec_indices
      use Lexer_Core, only: Print_Source
      use String_Table, only: Display_String

      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex, AnotherFieldIndex ! f_...
      character(len=*), optional :: String

      error = max(error,1)
      call output ( '***** At ' )
      call print_source ( source_ref(son) )
      call output ( ' RetrievalModule complained: ' )
      select case ( code )
      case ( bothOrNeither )
        call output ( 'One of ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' or ' )
        call display_string ( field_indices(anotherFieldIndex) )
        call output ( ' appears, but the other does not.', advance='yes' )
      case ( ifAThenB )
        call output ( 'If the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field appears then the ' )
        call display_string ( field_indices(anotherFieldIndex) )
        call output ( ' field shall also appear.', advance='yes' )
      case ( ifUnitsAThenB )
        call output ( 'If the units of the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field are ' )
        call output ( trim(string) )
        call output ( ' the ' )
        call display_string ( field_indices(anotherFieldIndex) )
        call output ( ' field shall also appear.', advance='yes' )
      case ( inconsistent, notSPD )
        call output ( 'the field ' )
        call display_string ( field_indices(fieldIndex) )
        select case ( code )
        case ( inconsistent )
          call output ( ' is not consistent with the ' )
          call display_string ( field_indices(anotherFieldIndex ) )
          call output ( ' field.', advance='yes' )
        case ( notSPD )
          call output ( ' is not a symmetric positive-definite matrix.', &
            & advance='yes' )
        end select
      case ( inconsistentUnits )
        call output ( 'the elements of the "' )
        call display_string ( field_indices(fieldIndex) )
        call output ( '" field have inconsistent units.', advance='yes' )
      case ( noFields )
        call output ( 'No fields are allowed for a ' )
        call display_string ( spec_indices(fieldIndex) )
        call output ( ' specification.', advance='yes' )
      case ( wrongUnits )
        call output ( 'The value(s) of the "' )
        call display_string ( field_indices(fieldIndex) )
        call output ( '" field shall have ' )
        call output ( trim(string) )
        call output ( ' units.', advance='yes' )
      end select
    end subroutine AnnounceError

    ! ----------------------------------------------  FillDiagVec  -----
    subroutine FillDiagVec ( QuantityIndex, Value )
      ! Put Value into the first element of the values field of the
      ! quantity of the Diagnostics vector given by QuantityIndex.
      ! We assume the Diagnostics vector is associated, so you better
      ! check before you call!
      integer, intent(in) :: QuantityIndex
      real(r8), intent(in) :: Value
      type(vectorValue_T), pointer :: Diag_Qty    ! A quantity in the
                                        ! Diagnostics vector
      diag_qty => GetVectorQuantityByType ( diagnostics, &
        & quantityType=QuantityIndex, noError=.true. )
      if ( associated(diag_qty) ) diag_qty%values(1,1) = value
    end subroutine FillDiagVec
    ! 
    ! ------------------------------------------  NewtonianSolver  -----
    subroutine NewtonianSolver

      use DNWT_Module, only: FlagName, NF_EVALF, NF_EVALJ, NF_SOLVE, &
        & NF_NEWX, NF_GMOVE, NF_BEST, NF_AITKEN, NF_DX, NF_DX_AITKEN, &
        & NF_START, NF_TOLX, NF_TOLX_BEST, NF_TOLF, NF_TOO_SMALL, &
        & NF_FANDJ, NWT, NWT_T, NWTA, NWTDB, RK
      use Dump_0, only: Dump
      use ForwardModelWrappers, only: ForwardModel
      use ForwardModelIntermediate, only: ForwardModelIntermediate_T, &
        & ForwardModelStatus_T
      use MatrixModule_1, only: AddToMatrix, CholeskyFactor, ClearMatrix, &
        & ColumnScale, CopyMatrix, CopyMatrixValue, CreateEmptyMatrix, &
        & DestroyMatrix, dump_Linf, dump_struct, &
        & FormNormalEquations => NormalEquations, &
        & GetDiagonal, InvertCholesky, Matrix_T, Matrix_Database_T, &
        & Matrix_Cholesky_T, Matrix_SPD_T, MaxL1, MinDiag, Multiply, &
        & Negate, RowScale, ScaleMatrix, SolveCholesky, UpdateDiagonal
      use Regularization, only: Regularize
      use VectorsModule, only: AddToVector, DestroyVectorInfo, &
        & DestroyVectorValue, Dump, Multiply, operator(.DOT.), operator(-), &
        & ScaleVector, SubtractFromVector

      ! Local Variables
      type(nwt_T) :: AJ                 ! "About the Jacobian", see NWT.
      real(r8) :: AprioriNorm           ! apriori .dot. apriori
      real(r8) :: Cosine                ! Of an angle between two vectors
      type(matrix_Cholesky_T) :: Factored ! Cholesky-factored normal equations
      type (ForwardModelStatus_T) :: FmStat ! Status for forward model
      type (ForwardModelIntermediate_T) :: Fmw ! Work space for forward model
      type(vector_T) :: FuzzState       ! Random numbers to fuzz the state
      integer :: J, K                   ! Loop inductors and subscripts
      type(matrix_SPD_T) :: NormalEquations  ! Jacobian**T * Jacobian
      integer :: NumF, NumJ             ! Number of Function, Jacobian evaluations
      integer :: NWT_Flag               ! Signal from NWT, q.v., indicating
                                        ! the action to take.
      integer :: NWT_Opt(20)            ! Options for NWT, q.v.
      real(rk) :: NWT_Xopt(20)          ! Real parameters for NWT options, q.v.
      integer :: RowBlock               ! Which block of rows is the forward
                                        ! model filling?
      integer, parameter :: SnoopLevels(NF_DX_AITKEN:NF_FANDJ) = (/ &
      ! dx_aitken dx aitken best gmove newx solve evalj evalf
        &      3, 2,     3,   2,    2,   1,    2,    2,    2,  &
      ! start tolx tolx_best tolf too_small fandj
        &  9,   2,        2,   2,        3,    9  /)
      real :: T1
      character(len=10) :: TheFlagName  ! Name of NWTA's flag argument

      call cpu_time ( t1 )
      call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
        & ModuleName )
      ! Set options for NWT
      nwt_opt(1:9) = (/  15, 1,      17, 2,      18, 3,      11, 4, 0 /)
      nwt_xopt(1:4) = (/ toleranceF, toleranceA, toleranceR, initLambda /)
      call nwt ( nwt_flag, nwt_xopt, nwt_opt )
      ! Create the matrix for the Cholesky factor of the normal equations
      call createEmptyMatrix ( factored%m, 0, state, state )
      ! Create the normal equations matrix
      call createEmptyMatrix ( normalEquations%m, 0, state, state )
      ! Create the vectors we need.
      call copyVector ( v(x), state, vectorNameText='_x', clone=.true. ) ! x := state
      call cloneVector ( v(atb), v(x), vectorNameText='_ATb' )
      call cloneVector ( v(bestGradient), v(x), vectorNameText='_bestGradient' )
      call cloneVector ( v(bestX), v(x), vectorNameText='_bestX' )
      call cloneVector ( v(candidateDX), v(x), vectorNameText='_candidateDX' )
      call cloneVector ( v(covarianceDiag), v(x), vectorNameText='_covarianceDiag' )
      call cloneVector ( v(dx), v(x), vectorNameText='_DX' )
      call cloneVector ( v(dxUnScaled), v(x), vectorNameText='_DxUnscaled' )
      call copyVector ( v(f_rowScaled), measurements, clone=.true., & ! mask only
        & noValues=.true., vectorNameText='_f_rowScaled' )
      call cloneVector ( v(gradient), v(x), vectorNameText='_gradient' )
      call cloneVector ( v(reg_X_x), measurements, vectorNameText='_reg_X_x' )
      if ( got(f_measurementSD) ) then
        call cloneVector ( v(weight), measurementSD, vectorNameText='_weight' )
        do j = 1, measurementSD%template%noQuantities
          where ( measurementSD%quantities(j)%values <= 0.0 )
            v(weight)%quantities(j)%values = 1.0
          elsewhere
            v(weight)%quantities(j)%values = 1.0 / &
              & measurementSD%quantities(j)%values
          end where
        end do
      end if
      if ( columnScaling /= l_none ) &
        call cloneVector ( v(columnScaleVector), v(x), vectorNameText='_ColumnScale' )
      if ( got(f_fuzz) ) then
        ! Add some fuzz to the state vector (for testing purposes):
        call cloneVector ( fuzzState, v(x) )
        do j = 1, v(x)%template%noQuantities
          call random_number(fuzzState%quantities(j)%values)
          v(x)%quantities(j)%values = v(x)%quantities(j)%values * &
            & ( 1.0_r8 + fuzz * ( fuzzState%quantities(j)%values - 0.5 ) )
        end do
      end if
        if ( index(switches,'xvec') /= 0 ) call dump ( v(x), name='Original X' )
      numF = 0
      numJ = 0
      aj%axmax = 0.0
      do k = 1, size(v(x)%quantities)
        aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(k)%values)))
      end do

        if ( index(switches,'sca') /= 0 ) then
          if ( got(f_regWeight) ) then
            call output ( ' regWeight = ' )
            call output ( regWeight )
          end if
          if ( got(f_apriori) ) then
            call output ( ' apriori scale = ' )
            call output ( aprioriScale )
          end if
          if ( got(f_fuzz) ) then
            call output ( ' fuzz = ' )
            call output ( fuzz )
          end if
          call output ( ' initLambda = ' )
          call output ( initLambda, advance='yes' )
        end if;

      call add_to_retrieval_timing( 'newton_solver', t1 )
      call cpu_time ( t1 )
      do ! Newtonian iteration
        if ( nwt_flag /= nf_start .and. index(switches,'ndb') /= 0 ) &
            & call nwtdb
        call nwta ( nwt_flag, aj )
          if ( index(switches,'nwt') /= 0 ) then
            call FlagName ( nwt_flag, theFlagName )
            call output ( 'Newton method flag = ' )
            call output ( trim(theFlagName) )
            call output ( ', numF = ' )
            call output ( numF )
            call output ( ', numJ = ' )
            call output ( numJ, advance='yes' )
          end if
        select case ( nwt_flag )
        case ( nf_evalf ) ! ..............................  EVALF  .....
        !{IF ( too many function values ) EXIT\\
        ! Compute $\mathbf{f(x)}$\\
        ! Compute the Jacobian matrix $J$ if you feel like it\\
        ! Set  {\bf AJ\%FNORM} = L2 norm of $\mathbf{f(x)}$\\
        ! IF ( AJ%FNORM is small enough ) EXIT
          numF = numF + 1
          if ( numF > maxFunctions ) then
              if ( index(switches,'nwt') /= 0 ) then
                call output ( &
                  & 'Newton iteration terminated because function evaluations (' )
                call output ( numF )
                call output ( ') > maxFunctions (' )
                call output ( maxFunctions )
                call output ( ')', advance='yes' )
              end if
            exit
          end if

          if ( got(f_apriori) ) then
            ! Destroy (i.e. clean up) the previous contents of
            ! AprioriMinusX (if any), so as not to have a memory
            ! leak.  Then make it look like apriori.
            call cloneVector ( v(aprioriMinusX), apriori, &
              & vectorNameText='_aprioriMinusX' )
            v(aprioriMinusX) = apriori - v(x)
            if ( got(f_aprioriScale) ) &
              & call scaleVector ( v(aprioriMinusX), aprioriScale )
            !{Let the covariance of the apriori be $\mathbf{S_a}$, let
            ! $\mathbf{C = S_a^{-1}}$, and let $\mathbf{F^T F = C}$. In
            ! the least-squares problem we have extra rows of the form
            ! $\mathbf{F \delta x = F ( a - x_n )}$ where $\mathbf{a}$
            ! is the apriori state, and $\mathbf{x_n}$ is the state as
            ! of the previous iteration.  The initial bit of the
            ! right-hand side of the normal equations is $\mathbf{F^T F
            ! ( a - x_n ) = C ( a - x_n )}$.
            call cloneVector ( v(covarianceXApriori), apriori, &
              & vectorNameText='_covarianceXApriori' )
            if ( diagonal ) then
              call cloneVector ( v(covarianceDiag), state, &
               &  vectorNameText='_covarianceDiag' )
              call getDiagonal ( covariance%m, v(covarianceDiag) )
              ! covarianceXApriori := covarianceDiag # apriori:
              call multiply ( v(covarianceDiag), v(aprioriMinusX), &
                & v(covarianceXApriori) )
            else ! covarianceXApriori := covariance X apriori:
              call multiply ( covariance, v(aprioriMinusX), &
                & v(covarianceXApriori) )
            end if
            !{The contribution to the norm of the residual of the
            ! right-hand side of the part of the least-squares problem
            ! that is due to apriori is $\mathbf{
            ! ( a - x_n )^T F^T F ( a - x_n ) = ( a - x_n )^T C ( a - x_n )}$
            aprioriNorm = v(aprioriMinusX) .dot. v(covarianceXapriori)
          else
            aprioriNorm = 0.0_r8
          end if

          ! Compute f(x)
 
          fmStat%newScanHydros = .true.
          fmStat%maf = 0

          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          call clearVector ( v(f) )
          ! Loop over MAFs
          do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
            fmStat%maf = fmStat%maf + 1
            do k = 1, size(configIndices)
              call forwardModel ( configDatabase(configIndices(k)), &
                & v(x), fwdModelExtra, v(f), fmw, fmStat )
            end do ! k
          end do ! MAFs
          call add_to_retrieval_timing( 'forward_model', t1 )
          call cpu_time ( t1 )
          call subtractFromVector ( v(f), measurements )
          if ( got(f_measurementSD) ) call multiply ( v(f), v(weight) )
          aj%fnorm = sqrt ( aprioriNorm + ( v(f) .dot. v(f) ) )
            if ( index(switches,'fvec') /= 0 ) &
              & call dump ( v(f), name='Residual' )
            if ( index(switches,'sca') /= 0 ) then
                call output ( ' | F | = ' )
                call output ( aj%fnorm, format='(1pe14.7)', advance='yes' )
            end if
          if ( aj%fnorm < toleranceF ) then
              if ( index(switches,'nwt') /= 0 ) then
                call output ( &
                  & 'Newton iteration terminated because aj%fnorm (' )
                call output ( aj%fnorm )
                call output ( ') < ToleranceF (' )
                call output ( toleranceF )
                call output ( ')', advance='yes' )
              end if
            exit
          end if
        case ( nf_evalj ) ! ..............................  EVALJ  .....
        !{IF ( too many Jacobian values ) EXIT \\
        ! Compute the Jacobian matrix J if you didn't do it when NFLAG
        ! was NF\_EVALF:\\
        ! ${\bf J}_{k,l} = \frac{\partial {\bf f}_k}{\partial {\bf x}_l},
        ! \, k = 1 .. \text{NF},\, l = 1 .. \text{NX}$
        !
        ! Actually, the matrix has four parts:
        !  \begin{xalignat*}{2}
        !  {\bf J \delta x} & \simeq {\bf -f}  && \text{The Jacobian} \\
        !  {\bf F x_{n+1}}  & \simeq {\bf F a} && \text{The apriori} \\
        !  {\bf R x_{n+1}}  & \simeq {\bf 0}   &&
        !                   \text{Tikhonov regularization} \\
        !  \lambda {\bf I \delta x} & \simeq {\bf 0} &&
        !                   \text{Levenberg-Marquardt stabilization}
        !  \end{xalignat*}
        !%
        ! {\bf OR}, written so that all equations are solving for ${\bf
        ! \delta x}$:
        !%
        ! \begin{equation*}
        ! \begin{array}{ll}
        !  \left [
        !  \begin{array}{l}
        !  {\bf J} \\
        !  {\bf F} \\
        !  {\bf R} \\
        !  \lambda {\bf I}
        !  \end{array} \right ] {\bf \delta x}
        !  \simeq
        !  \left [
        !  \begin{array}{l}
        !   {\bf -f} \\
        !   {\bf F} ( {\bf a} - {\bf x}_n) \\
        !   -{\bf R x}_n \\
        !   {\bf 0}
        !  \end{array}
        !  \right ]
        ! &
        !  \begin{array}{l}
        !   \text{The Jacobian} \\
        !   \text{The apriori} \\
        !   \text{Tikhonov regularization} \\
        !   \text{Levenberg-Marquardt stabilization}
        !  \end{array}
        ! \end{array}
        ! \end{equation*}
        !%
        ! Triangularize ${\bf J}$, and compute the (negative of the) gradient
        ! = $-{\bf J^T f}$.  This is the RHS of the normal equations
        ! ${\bf J^T J \delta \hat x} = -{\bf J^T f}$ where ${\bf \delta
        !  \hat x}$ is a "Candidate DX" that might not actually get used.
        ! Set \begin{description}
        !   \item[AJ\%DIAG] = element on diagonal with smallest absolute
        !           value, after triangularization,
        !   \item[AJ\%AJN] = maximum L1 norm of column in upper triangle
        !           after triangularization,
        !   \item[AJ\%FNMIN] = L2 norm of residual =
        !           $||{\bf f + J \delta \hat x}||_2$,
        !   \item[AJ\%GRADN] = L2 norm of Gradient = $|| {\bf J^T f}||_2$.
        !   \end{description}
          numJ = numJ + 1
          if ( numJ > maxJacobians ) then
              if ( index(switches,'nwt') /= 0 ) then
                call output ( &
                  & 'Newton iteration terminated because Jacobian evaluations (' )
                call output ( numJ )
                call output ( ') > maxJacobians (' )
                call output ( maxJacobians )
                call output ( ')', advance='yes' )
              end if
            exit
          end if
          update = got(f_apriori)
          if ( update ) then ! start normal equations with apriori
            if ( diagonal ) then
            ! Iterate with the diagonal of the apriori covariance, then
            ! use the full apriori covariance (one hopes only for one
            ! more iteration). This improves sparsity during iteration.
              call clearMatrix ( normalEquations%m )
              call updateDiagonal ( normalEquations, v(covarianceDiag) )
            else
              call copyMatrixValue ( normalEquations%m, covariance%m )
            end if
            call copyVector ( v(atb), v(covarianceXApriori) ) ! A^T b := C (a-x_n)
            if ( got(f_aprioriScale) ) &
              & call scaleMatrix ( normalEquations%m, aprioriScale )
          else
            call clearMatrix ( normalEquations%m ) ! start with zero
            call destroyVectorValue ( v(atb) ) ! Clear the RHS vector
          end if
          ! aprioriNorm was computed when nwt_flag was NF_EVALF
          aj%fnorm = aprioriNorm

          ! Add Tikhonov regularization if requested
          if ( got(f_regOrders) ) then
            call regularize ( jacobian, regOrders, regQuants, regWeight )
            call multiplyMatrixVectorNoT ( jacobian, v(x), v(reg_X_x) ) ! regularization * x_n
            call scaleVector ( v(reg_X_x), -regWeight )   ! -R x_n
            call formNormalEquations ( jacobian, normalEquations, &
              & v(reg_X_x), v(atb), update=update, useMask=.false. )
            update = .true.
            call clearMatrix ( jacobian )           ! free the space
            aj%fnorm = aj%fnorm + ( v(reg_X_x) .dot. v(reg_X_x) )
!           call destroyVectorValue ( v(reg_X_x) )  ! free the space
            ! Don't destroy reg_X_x unless we move the 'clone' for it
            ! inside the loop.  Also, if we destroy it, we can't snoop it.
          end if

          ! Add some early stabilization
          if ( got(f_lambda) ) then
            call updateDiagonal ( normalEquations, initLambda )
            update = .true.
            ! The right-hand-side is zero -- no need to update ATb
            ! or aj%fnorm
          end if

          fmStat%maf = 0
          fmStat%newScanHydros = .true.

          ! Loop over MAFs
          call clearVector ( v(f) )
          do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
            call add_to_retrieval_timing( 'newton_solver', t1 )
            call cpu_time ( t1 )
            ! What if one config set finished but others still had more
            ! to do? Ermmm, think of this next time.
            fmStat%maf = fmStat%maf + 1
            fmstat%rows = .false.
            do k = 1, size(configIndices)
              call forwardModel ( configDatabase(configIndices(k)), &
                & v(x), fwdModelExtra, v(f_rowScaled), fmw, fmStat, jacobian )
            end do ! k
            call add_to_retrieval_timing( 'forward_model', t1 )
            call cpu_time ( t1 )
            do rowBlock = 1, size(fmStat%rows)
              if ( fmStat%rows(rowBlock) ) then
                call subtractFromVector ( v(f_rowScaled), measurements, &
                  & quant=jacobian%row%quant(rowBlock), &
                  & inst=jacobian%row%inst(rowBlock) ) ! f - y
                !{Let ${\bf J}$ be the Jacobian matrix and ${\bf f}$ be the
                ! residual of the least-squares problem.  Let $\bf W$ be the
                ! inverse of the measurement covariance (which in our case
                ! is diagonal). Row scale the part of the least-squares
                ! problem that arises from the measurements, i.e. the
                ! least-squares problem becomes $\mathbf{W J = W f}$
                ! (actually, we only row scale ${\bf f}$ here, and scale
                ! ${\bf J}$ below).
                if ( got(f_measurementSD) ) then
                  call multiply ( v(f_rowScaled), v(weight), &
                    & quant=jacobian%row%quant(rowBlock), &
                    & inst=jacobian%row%inst(rowBlock) )
                end if
              end if
            end do
            call scaleVector ( v(f_rowScaled), -1.0_r8 ) ! y - f

            if ( got(f_measurementSD) ) call rowScale ( v(weight), jacobian )

            !{Form normal equations:
            ! $\mathbf{J^T W^T W J \delta \hat x = J^T W^T W f}$:
            call formNormalEquations ( jacobian, normalEquations, &
              & rhs_in=v(f_rowScaled), rhs_out=v(atb), update=update, &
              & useMask=.true. )
            update = .true.
              if ( index(switches,'jac') /= 0 ) &
                call dump_Linf ( jacobian, 'L_infty norms of Jacobian blocks:' )
              if ( index(switches,'spa') /= 0 ) &
                & call dump_struct ( jacobian, &
                  & 'Sparseness structure of Jacobian blocks:' )
            call clearMatrix ( jacobian )  ! free the space
          end do ! mafs

          aj%fnorm = aj%fnorm + ( v(f_rowScaled) .dot. v(f_rowScaled) )

          !{Column Scale $\mathbf{J}$ (row and column scale $\mathbf{J^T
          ! J}$).  Let $\mathbf{\Sigma}$ be a column-scaling matrix for the
          ! Jacobian. We cater for several choices.  After row and column
          ! scaling, the problem is $\mathbf{J^T W^T W J \Sigma \Sigma^{-1}
          ! \delta \hat x = -J^T W^T f}$.
          select case ( columnScaling )
          case ( l_apriori )
            call copyVector ( v(columnScaleVector), apriori )
          case ( l_covariance )
            !??? Can't get here until allowed by init_tables
          case ( l_norm )
            call getDiagonal ( normalEquations%m, v(columnScaleVector) )
          end select
          if ( columnScaling /= l_none ) then ! Compute $\Sigma$
            forall (j = 1: v(columnScaleVector)%template%noQuantities)
              where ( v(columnScaleVector)%quantities(j)%values <= 0.0 )
                v(columnScaleVector)%quantities(j)%values = 1.0
              elsewhere
                v(columnScaleVector)%quantities(j)%values = 1.0 / &
                  & sqrt( v(columnScaleVector)%quantities(j)%values )
              end where
            end forall
              if ( index(switches,'col') /= 0 ) &
                & call dump ( v(columnScaleVector), name='Column scale vector' )
            !{Scale: $\mathbf{\Sigma^T J^T J \Sigma = -\Sigma^T J^T f}$
            call columnScale ( normalEquations%m, v(columnScaleVector) )
            call rowScale ( v(columnScaleVector), normalEquations%m )
            call multiply ( v(atb), v(columnScaleVector) )
  
          end if
          !{Compute the (negative of the) gradient $= -\mathbf{J^T f}$.
          ! This is the right-hand side of the normal equations
          ! $\mathbf{J^T J \delta \hat x = -J^T f}$ where $\mathbf{\delta
          ! \hat x}$ is the "Candidate $\mathbf{\delta x}$" that may or
          ! may not get used.
          call copyVector ( v(gradient), v(atb) ) ! gradient := atb
            if ( index(switches,'gvec') /= 0 ) &
              & call dump ( v(gradient), name='gradient' )

          ! Factor J^T J:
!         factored%m%block => normalEquations%m%block ! to save space
!         Can't do the above because we need to keep the normal
!         equations around, in order to subtract Levenberg-Marquardt and
!         apriori covariance, in order to compute a posteriori covariance
            if ( index(switches,'neq') /= 0 ) &
              call dump_Linf ( normalEquations%m, &
                & 'L_infty norms of Normal Equations blocks after scaling:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          ! Factor the normal equations
          call choleskyFactor ( factored, normalEquations )
          call add_to_retrieval_timing( 'cholesky_factor', t1 )
          call cpu_time ( t1 )
            if ( index(switches,'diag') /= 0 ) then
              call getDiagonal ( factored%m, v(dxUnscaled) )
              call dump ( v(dxUnscaled), &
                & name='Diagonal of factored normal equations:' )
            end if
            if ( index(switches,'fac') /= 0 ) &
              call dump_Linf ( factored%m, 'L_infty norms of blocks of factor:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( factored%m, &
                & 'Sparseness structure of blocks of factor:', upper=.true. )
          aj%diag = minDiag ( factored ) ! element on diagonal with
          !       smallest absolute value, after triangularization
          aj%ajn = maxL1 ( factored%m ) ! maximum L1 norm of
          !       column in upper triangle after triangularization
          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          call solveCholesky ( factored, v(candidateDX), v(atb), &
            & transpose=.true. )
          call add_to_retrieval_timing( 'cholesky_solver', t1 )
          call cpu_time ( t1 )

          !{AJ\%FNMIN = L2 norm of residual, $||\mathbf{J \delta x + f}||$
          ! where $\mathbf{\delta x}$ is the "Candidate DX" that may not
          ! get used. This can be gotten without saving $\bf J$ as
          ! $\mathbf{f^T f - y^T y}$ where $\mathbf{y = U^{-T} \delta x}$
          ! ($\mathbf{U}$ is the Cholesky factor of $\mathbf{J^T J}$).
          ! The variable {\tt candidateDX} is a temp here = $\bf y$.
          aj%fnmin = aj%fnorm - (v(candidateDX) .dot. v(candidateDX))
          if ( aj%fnmin < 0.0 ) then
            call output ( 'How can aj%fnmin be negative?  aj%fnmin = ' )
            call output ( aj%fnmin, advance='yes' )
            call output ( 'aj%fnorm = ' )
            call output ( aj%fnorm, advance='yes' )
            call output ( 'norm(candidateDX) = ' )
            call output ( v(candidateDX) .dot. v(candidateDX), advance='yes' )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Norm of residual is imaginary!' )
          end if
          aj%fnmin = sqrt(aj%fnmin)
          aj%fnorm = sqrt(aj%fnorm)
          aj%gradn = sqrt(v(gradient) .dot. v(gradient)) ! L2Norm(gradient)
            if ( index(switches,'sca') /= 0 ) then
              call dump ( (/ aj%ajn, aj%diag, aj%fnmin, aj%gradn /), &
                & ' L1| FAC |       aj%diag      aj%fnmin         | G |', &
                & clean=.true. )
            end if
        case ( nf_solve ) ! ..............................  SOLVE  .....
        !{Apply Levenberg-Marquardt stabilization with parameter
        ! $\lambda =$ {\bf AJ\%SQ}.  I.e. solve $\mathbf{(J^T W^T W J
        ! \Sigma + \lambda^2 I) \Sigma^{-1} \delta \hat x = J^T W^T
        ! f}$ for $\mathbf{\Sigma^{-1} \delta \hat x}$.  Set
        ! \begin{description}
        !   \item[AJ\%FNMIN] as for NWT\_FLAG = NF\_EVALJ, but taking
        !     account of Levenberg-Marquardt stabilization.
        !   \item[AJ\%DXN] = L2 norm of "candidate DX",
        !   \item[AJ\%GDX] = (Gradient) .dot. ("candidate DX")
        ! \end{description}
          call updateDiagonal ( normalEquations, aj%sq**2 )
!         factored%m%block => normalEquations%m%block ! to save space
!         Can't do the above because we need to keep the normal equations
!         around, in order to subtract Levenberg-Marquardt and apriori
!         covariance, in order to compute a posteriori covariance
            if ( index(switches,'neq') /= 0 ) &
              call dump_Linf ( normalEquations%m, &
                & 'L1 norms of Normal Equations blocks after Marquardt:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          call choleskyFactor ( factored, normalEquations )
          call add_to_retrieval_timing( 'cholesky_factor', t1 )
          call cpu_time ( t1 )
            if ( index(switches,'fac') /= 0 ) &
              call dump_Linf ( factored%m, &
                & 'L1 norms of blocks of factor after Marquardt:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( factored%m, &
                & 'Sparseness structure of blocks of factor:', upper=.true. )
          !{Solve for "candidate DX" = $-\mathbf{(J^T J)^{-1} J^T  F =
          ! -(U^T U)^{-1} J^T  F}$ using two back solves.  First solve
          ! $\mathbf{U^T y = -J^T F}$, then $\mathbf{U}$ "candidate DX"
          ! = $\mathbf{y}$. Meanwhile, set AJ\%FNMIN as for NWT\_FLAG =
          ! NF\_EVALJ, but taking account of Levenberg-Marquardt
          ! stabilization
          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          call solveCholesky ( factored, v(candidateDX), v(atb), &
            & transpose=.true. )
          call add_to_retrieval_timing( 'cholesky_solver', t1 )
          call cpu_time ( t1 )
          aj%fnmin = sqrt(aj%fnorm**2 - (v(candidateDX) .dot. v(candidateDX)) )
          call solveCholesky ( factored, v(candidateDX) )
          aj%dxn = sqrt(v(candidateDX) .dot. v(candidateDX)) ! L2Norm(dx)
          aj%gdx = v(gradient) .dot. v(candidateDX)
          if ( .not. aj%starting ) aj%dxdxl = v(dx) .dot. v(candidateDX)
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(candidateDX), name='CandidateDX' )
            if ( index(switches,'sca') /= 0 ) then
              cosine = -2.0_r8
              if ( aj%dxn > 0.0_r8 .and. aj%gradn > 0.0_r8 ) &
                cosine = aj%gdx/(aj%dxn*aj%gradn)
              call dump ( (/ aj%dxn, aj%fnmin, cosine, aj%sq /), &
                & '    | DX |      aj%fnmin    ' // &
                & 'cos(G, DX)        lambda', clean=.true. )
              if ( .not. aj%starting ) then
                call output ( ' cos(DX, DXL) = ' )
                cosine = -2.0_r8
                if ( aj%dxdxl > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
                  cosine = aj%dxdxl / (aj%dxn*aj%dxnl)
                call output ( cosine, format='(1pe14.7)', advance='yes' )
              ! call output ( ',  DX . DXL = ' )
              ! call output ( aj%dxdxl, format='(1pe14.7)', advance='yes' )
              end if
            end if
        case ( nf_newx ) ! ................................  NEWX  .....
        ! Set X = X + DX
        !     AJ%AXMAX = MAXVAL(ABS(X)),
        !     AJ%BIG = ANY ( DX > 10.0 * epsilon(X) * X )
          if ( .not. aj%starting ) aj%dxdxl = v(dx) .dot. v(candidateDX)
          !{Account for column scaling.  We solved for $\mathbf{\Sigma^{-1}
          ! \delta x}$ above, so multiply by $\Sigma$ (which is our
          ! variable {\tt columnScaleVector}):
          call copyVector ( v(dxUnScaled), v(dx) )
          if ( columnScaling /= l_none ) then
            ! dxUnScaled = dxUnScaled # columnScaleVector:
            call multiply ( v(dxUnScaled), v(columnScaleVector) )
          end if
          call addToVector ( v(x), v(dxUnScaled) ) ! x = x + dxUnScaled
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dxUnScaled), name='dX Unscaled' )
            if ( index(switches,'xvec') /= 0 ) &
              & call dump ( v(x), name='New X' )
          aj%axmax = 0.0
          aj%big = .false.
          do j = 1, size(v(x)%quantities)
            aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(j)%values)))
            aj%big = aj%big .or. any( abs(v(dx)%quantities(j)%values) > &
              & 10.0 * epsilon(aj%axmax) * abs(v(x)%quantities(j)%values) )
          end do
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%axmax = ' )
              call output ( aj%axmax, format='(1pe14.7)' )
              if ( .not. aj%starting ) then
                call output ( ' cos(DX, DXL) = ' )
                cosine = -2.0_r8
                if ( aj%dxdxl > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
                  cosine = aj%dxdxl / (aj%dxn*aj%dxnl)
                call output ( cosine, format='(1pe14.7)' )
              ! call output ( ', DX . DXL = ' )
              ! call output ( aj%dxdxl, format='(1pe14.7)', advance='yes' )
              end if
              call output ( ', aj%big = ' )
              call output ( aj%big, advance='yes' )
            end if
        case ( nf_gmove ) ! ..............................  GMOVE  .....
        ! Set X = "Best X"
        !     DX = AJ%GFAC * "Best Gradient"
          call copyVector ( v(x), v(bestX) ) ! x = bestX
          ! dx = aj%gfac * "Best Gradient":
          call scaleVector ( v(bestGradient), aj%gfac, v(dx) )
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dx), name='Gradient move from best X' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%gfac = ' )
              call output ( aj%gfac, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_best ) ! ................................  BEST  .....
        ! Set "Best X" = X, "Best Gradient" = Gradient
          call copyVector ( v(bestX), v(x) ) ! bestX = x
          call copyVector ( v(bestGradient), v(gradient) ) ! bestGradient = gradient
        case ( nf_aitken ) ! ............................  AITKEN  .....
        ! Set DX = DX - "Candidate DX",
        !     AJ%DXDX = dot_product( DX, DX )
        ! IF ( AJ%DXDX /= 0.0 ) &
        !   Set AJ%DXDXL = dot_product( DX, "Candidate DX" )
          call subtractFromVector ( v(dx), v(candidateDX) ) ! dx = dx - candidateDX
          aj%dxdx = v(dx) .dot. v(dx)
          if ( aj%dxdx > 0.0 ) aj%dxdxl = v(dx) .dot. v(candidateDX)
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dx), name='dx after Aitken' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%dxdx = ' )
              call output ( aj%dxdx, format='(1pe14.7)' )
              call output ( ', cos(DX, DXL) = ' )
              cosine = -2.0_r8
              if ( aj%dxdxl > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
                cosine = aj%dxdxl / (aj%dxn*aj%dxnl)
              call output ( cosine, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_dx ) ! ....................................  DX  .....
          call copyVector ( v(dx), v(candidateDX) ) ! dx = candidateDX
        case ( nf_dx_aitken ) ! ......................  DX_AITKEN  .....
          ! dx = aj%cait * candidateDX:
          call scaleVector ( v(candidateDX), aj%cait, v(dx) )
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dx), name='dx after dx Aitken' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%cait = ' )
              call output ( aj%cait, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_tolx, nf_tolx_best, nf_tolf, nf_too_small ) ! ........
          ! IF ( NWT_FLAG == NF_TOO_SMALL ) THEN
          !   Take special action if requested accuracy is critical
          ! END IF
          if ( nwt_flag == nf_tolx_best ) call copyVector ( v(x), v(bestX) )
          ! Convergence to desired solution.  Do whatever you want to
          ! with the solution.
          if ( .not. got(f_apriori) .or. .not. diagonal ) then
              if ( index(switches,'nwt') /= 0 )&
                & call output ( &
                  & 'Newton iteration terminated because of convergence', &
                  & advance='yes' )
            exit
          end if
          diagonal = .not. diagonal
          nwt_flag = nf_start
        case ( nf_fandj ) ! ...............................  FANDJ .....
          ! There is probably an error in the way F or J is computed.
          ! A warning has been printed by the error processor.
          ! IF ( you have confidence in F and J ) CYCLE
          ! STOP
        end select
      ! IF ( you want to return to a previous best X ) NWT_FLAG = 0
        if ( got(f_diagnostics) ) then
          call fillDiagVec ( l_dnwt_ajn, aj%ajn )
          call fillDiagVec ( l_dnwt_axmax, aj%axmax )
          call fillDiagVec ( l_dnwt_cait, aj%cait )
          call fillDiagVec ( l_dnwt_diag, aj%diag )
          call fillDiagVec ( l_dnwt_dxdx, aj%dxdx )
          call fillDiagVec ( l_dnwt_dxdxl, aj%dxdxl )
          call fillDiagVec ( l_dnwt_dxn, aj%dxn )
          call fillDiagVec ( l_dnwt_dxnl, aj%dxnl )
          call fillDiagVec ( l_dnwt_flag, real(nwt_flag,r8) )
          call fillDiagVec ( l_dnwt_fnmin, aj%fnmin )
          call fillDiagVec ( l_dnwt_fnorm, aj%fnorm )
          call fillDiagVec ( l_dnwt_gdx, aj%gdx )
          call fillDiagVec ( l_dnwt_gfac, aj%gfac )
          call fillDiagVec ( l_dnwt_gradn, aj%gradn )
          call fillDiagVec ( l_dnwt_sq, aj%sq )
          call fillDiagVec ( l_dnwt_sqt, aj%sqt )
        end if
        if ( snoopKey /= 0 .and. snoopLevel >= snoopLevels(nwt_flag) ) then
          call FlagName ( nwt_flag, theFlagName )
          call snoop ( snoopKey, vectorDatabase, v, &
            & trim(snoopComment) // ': ' // trim(theFlagName) )
        end if
      end do ! Newton iteration
      if ( got(f_outputCovariance) .or. got(f_outputSD) ) then
        if ( got(f_diagnostics) ) then
          ! Compute rows of Jacobian actually used.  Don't count rows due
          ! to Levenberg-Marquardt stabilization.  Do count rows due to a
          ! priori or regularization.  Put numbers of rows and columns
          ! into diagnostic vector.
          jacobian_cols = sum(normalEquations%m%col%nelts)
          jacobian_rows = sum(normalEquations%m%row%nelts)
          do j = 1, normalEquations%m%col%nb
            if ( associated(normalEquations%m%col%vec%quantities(j)%mask) ) &
              & jacobian_rows = jacobian_rows - &
              & countBits(normalEquations%m%col%vec%quantities(j)%mask)
          end do
          if ( got(f_apriori) ) &
            & jacobian_rows = jacobian_rows + jacobian_cols
          if ( got(f_regOrders) ) &
            & jacobian_rows = jacobian_rows + jacobian_cols
          call fillDiagVec ( l_jacobian_cols, real(jacobian_cols,r8) )
          call fillDiagVec ( l_jacobian_cols, real(jacobian_rows,r8) )
        end if
        ! Subtract sum of Levenberg-Marquardt updates from normal
        ! equations
        call updateDiagonal ( normalEquations, -aj%sqt**2 )
        ! Subtract a priori covariance matrix from normal equations
        ! call negate ( covariance%m )
        ! call addToMatrix ( normalEquations%m, covariance%m )
        ! call negate ( covariance%m )
        !??? Commented out, pending deep thought.!???
        !??? If we remove a priori covariance from normal equations
        !??? we should presumably also remove regularization.  Maybe we
        !??? should remove regularization anyway.  On the other hand, if
        !??? remove them, maybe the normal equations are singular?
        ! Re-factor normal equations
        call add_to_retrieval_timing( 'newton_solver', t1 )
        call cpu_time ( t1 )
        call choleskyFactor ( factored, normalEquations )
        call add_to_retrieval_timing( 'cholesky_factor', t1 )
        call cpu_time ( t1 )
        call invertCholesky ( factored, outputCovariance%m )
        call add_to_retrieval_timing( 'cholesky_invert', t1 )
        call cpu_time ( t1 )
        !??? Don't forget to scale the covariance
        if ( associated(outputSD) ) then !???
          call GetDiagonal ( outputCovariance%m, outputSD )
          if ( columnScaling /= l_none ) then
            call multiply ( outputSD, v(columnScaleVector) )
          end if
        end if
      end if
      call copyVector ( state, v(x) )
        if ( index(switches,'svec') /= 0 ) &
          & call dump ( state, name='Final state' )
      ! Clean up the temporaries, so we don't have a memory leak.
      ! Most of them are in the myVectors database, which is destroyed
      ! upon exit from Retrieve.
      if ( got(f_fuzz) ) call destroyVectorInfo ( fuzzState )
      call destroyMatrix ( normalEquations%m )
      call destroyMatrix ( factored%m )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      call add_to_retrieval_timing( 'newton_solver', t1 )
      call cpu_time ( t1 )
    end subroutine NewtonianSolver
    ! ------------------------------------------  LowCloudRetrieval  -----
    subroutine LowCloudRetrieval

      use Intrinsic, only: L_PTAN, L_RADIANCE, &
                     & L_CLOUDINDUCEDRADIANCE,                               &
                     & L_CLOUDEXTINCTION,                                &
                     & L_CLOUDRADSENSITIVITY,                                &
                     & L_EARTHRADIUS,                                &
                     & L_LOSTRANSFUNC
      use MatrixModule_0, only: MatrixInversion, MATRIXELEMENT_T
      use MatrixModule_1, only: ClearMatrix, FINDBLOCK, GetDiagonal
      use MLSSignals_m, only: SIGNAL_T
      use ForwardModelWrappers, only: ForwardModel
      use ForwardModelIntermediate, only: ForwardModelIntermediate_T, &
        & ForwardModelStatus_T

      ! Local Variables
      type (ForwardModelStatus_T) :: FmStat        ! Status for forward model
      type (ForwardModelIntermediate_T) :: Fmw     ! Work space for forward model
      type (vector_T) :: FwdModelOut1               ! Forward outputs
      type (Signal_T) :: signal                    ! signal info in each model
      type (VectorValue_T), pointer :: xLosPtr        ! pointer of l_lostransfunc quantity
      type (VectorValue_T), pointer :: xExtPtr        ! pointer of l_cloudExtinction quantity
      type (VectorValue_T), pointer :: xLosVar        ! variance of apriori
      type (VectorValue_T), pointer :: xExtVar        ! variance of apriori
      type (VectorValue_T), pointer :: outLosSD        ! SD of output
      type (VectorValue_T), pointer :: outExtSD        ! SD of output
      type (VectorValue_T), pointer :: Tcir        ! cloud-induced radiance
      type (VectorValue_T), pointer :: Terr        ! cloud-induced radiance SD
      type (VectorValue_T), pointer :: PTAN        ! Tgt pressure
      type (VectorValue_T), pointer :: Re        ! Earth Radius
      type (VectorValue_T), pointer :: Tb0         ! model clear sky radiance (100%RH)
      type (VectorValue_T), pointer :: Slope       ! sensitivity slope to convert cloud
                                                   ! radiance to optical depth
                                          
      integer :: i,j,k,ich,imodel,mif,isignal          ! Loop subscripts
      integer :: coljBlock     ! Column index for jacobian
      integer :: rowjBlock     ! Row index for jacobian
      integer :: nFreqs      ! number of frequencies in each block
      integer :: nSgrid      ! number of S grids  
      integer :: nChans      ! total number of channels used in retrieval
      integer :: nMifs       ! number of total mifs
      integer :: ndoMifs       ! number of used mifs
      real(r8) :: p_lowcut    ! ptan threshold for low tangent heights
      real(r8) :: scale
                                        
      type(MatrixElement_T), pointer :: JBLOCK       ! A block from the jacobian
      type(vector_T) :: CovarianceDiag  ! Diagonal of apriori Covariance  
      
      ! retrieval work arrays
      real(r8) :: sensitivity                         ! sensitivity with slope and correction
      real(r8) :: teff                                ! effective optical depth
      real(r8) :: trans                                ! transmission function
      real(r8), dimension(:,:), allocatable :: A      ! working array
      real(r8), dimension(:,:), allocatable :: dx       ! working array for x
      real(r8), dimension(:,:), allocatable :: y      ! measurement array
      real(r8), dimension(:,:), allocatable :: sy     ! variance of y
      real(r8), dimension(:), allocatable :: Slevel        ! s grid
      real(r8), dimension(:), allocatable :: x        ! s grid array
      real(r8), dimension(:), allocatable :: x0       ! A priori of x
      real(r8), dimension(:), allocatable :: sx0       ! variance of A priori
      real(r8), dimension(:), allocatable :: xext       ! extinction along los
      real(r8), dimension(:,:), allocatable :: sx       ! variance of x
      real(r8), dimension(:,:,:), allocatable :: C    ! for problem y=Kx, c=K^t#Sy^-1#K
                                                      ! last dimension is for mif
                                                      ! first two are (chan, s)

      ! use this for testing
      if(.not. got(f_maxJ)) maxJacobians = 5
      if(.not. got(f_lambda)) initlambda = 10.
      
      p_lowcut = -2.5
      
      call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
        & ModuleName )
      ! create FwdModelOut1 vector by cloning measurements
        call cloneVector ( FwdModelOut1, measurements, vectorNameText='_covarianceDiag' )
               
      ! Measurements are cloud radiances and will be converted to
      ! the effective cloud optical depth after the forward model
      ! provides the sensitivities.
        
      ! Get  l_lostransfunc quantity from state vector
      ! State vector should contains only one quantity of l_lostransfunc type 
      ! (a los quantity),which is the increment of cloud transmission function,
      ! and it is going to be retrieved here.
      
        xLosPtr => GetVectorQuantityByType ( state, quantityType=l_lostransfunc)
        xExtPtr => GetVectorQuantityByType ( state, quantityType=l_cloudextinction)
      
      ! Get S grid dimensions
        nSgrid=xLosPtr%template%noChans    ! number of S grids
        allocate(Slevel(nSgrid))
        sLevel = xLosPtr%template%frequencies      ! sLevel has unit of km here                 
        
      ! find how many channels total in all models
        nChans = 0
        do imodel = 1, size(configIndices)
        do isignal = 1, size(configDatabase(configIndices(imodel))%signals)
          nChans = nChans + &
          & count(configDatabase(configIndices(imodel))%signals(isignal)%channels)
        end do ! band signals
        end do ! configIndices or models

      ! find how many mifs from the first quantity        
        nMifs = measurements%quantities(1)%template%noSurfs
               
      ! create covarianceDiag array
        call cloneVector ( covarianceDiag, state, vectorNameText='_covarianceDiag' )
      ! get the inverted diagnonal elements of covariance of apriori
        call getDiagonal ( covariance%m, covarianceDiag )
        
      ! allocate C, y, x matrices
          allocate(A(nSgrid,nSgrid),C(nSgrid,nSgrid,nMifs))
          allocate(y(nChans,nMifs),sy(nChans,nMifs))
          allocate(x(nSgrid),x0(nSgrid),sx0(nSgrid),xext(nSgrid))
          allocate(dx(nSgrid,nMifs),sx(nSgrid,nMifs))

        fmStat%maf = 0
          ! Loop over MAFs
        do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
          fmStat%maf = fmStat%maf + 1
print*,'begin cloud retrieval maf= ',fmstat%maf,' chunk size=',chunk%lastMAFIndex-chunk%firstMAFIndex
                        
          A = 0._r8
          C = 0._r8
          y = 0._r8
          sy = 0._r8
          dx = 0._r8
            
          ich = 0
          do imodel = 1, size(configIndices)
            fmStat%rows = .false.
            call forwardModel ( configDatabase(configIndices(imodel)), &
                & state, fwdModelExtra, FwdModelOut1, fmw, fmStat, jacobian )
            
            do isignal = 1, size(configDatabase(configIndices(imodel))%signals)
              ! get signal information for this model. Note: allow only 1 signal 
              signal = configDatabase(configIndices(imodel))%signals(isignal)

              ! get clear sky radiances from forward model for this signal
              Tb0 => GetVectorQuantityByType ( fwdModelOut1,                 &
               & quantityType=l_radiance,                                      &
               & signal=signal%index, sideband=signal%sideband )

              ! get ptan.  Warning: all signals must be in the same module
              ptan => GetVectorQuantityByType ( fwdModelExtra,      &
               & quantityType=l_ptan, instrumentModule = &
               & Tb0%template%instrumentModule)

              ! get earthradius from forward model
              Re => GetVectorQuantityByType ( fwdModelExtra,      &
               & quantityType=l_earthradius)
          
              ! get sensitivity from forward model for this signal
              Slope => GetVectorQuantityByType ( fwdModelOut1,      &
               & quantityType=l_cloudRADSensitivity,                           &
               & signal=signal%index, sideband=signal%sideband )
          
              ! get cloud radiance measurements for this signal
              Tcir => GetVectorQuantityByType ( Measurements,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )
            
              ! get cloud radiances error for this signal
              Terr => GetVectorQuantityByType ( MeasurementSD,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )

              ! get pointers of x covariance for the retrieval
               xLosVar => GetVectorQuantityByType (covarianceDiag, &
               & quantityType=l_lostransfunc,noerror=.true.)
               xExtVar => GetVectorQuantityByType (covarianceDiag, &
               & quantityType=l_cloudextinction,noerror=.true.)
      
              ! get pointers of output SD for the retrieval
               outLosSD => GetVectorQuantityByType (outputSD, &
               & quantityType=l_lostransfunc,noerror=.true.)
               outExtSD => GetVectorQuantityByType (outputSD, &
               & quantityType=l_cloudextinction,noerror=.true.)
               ! initialized to SD of xLosVar
               outLosSD%values = 1._r8/sqrt(xLosVar%values)
      
              ! get rowBlock and colBlock for this model
              rowJBlock = FindBlock (jacobian%row, Tb0%index, fmStat%maf)
              colJBlock = FindBlock ( Jacobian%col, xLosPtr%index, fmStat%maf )
                
              jBlock => jacobian%block(rowJblock,colJblock)
               
              nFreqs = size (signal%frequencies)
               do k=1,nFreqs
               if(signal%channels(k)) then
               ich = ich + 1
                  do mif=1,nMifs
                  if(ptan%values(mif,fmStat%maf) < p_lowcut) then
                     
                     ! find more accurate sensitivity. we first borrow 
                     ! x, x0 arrays to establish the Tcir-teff relation:
                     ! x-> Tcir; x0->teff;
                     do j=1,nSgrid
                     x0(j) = 1._r8/nSgrid*(j-1._r8)
                     x(j) = slope%values(k+nFreqs*(mif-1),fmStat%maf)* &
                        & (1._r8 + 0.46_r8* teff**6) ! correction term (see ATBD)
                     end do
                     
                     y(ich,mif) = Tcir%values(k+nFreqs*(mif-1),fmStat%maf)

                     ! interpolate to get initial guess of teff
                     if(y(ich,mif) < x(1)) teff=0._r8
                     if(y(ich,mif) > x(nSgrid)) teff=1._r8
                     do j=1,nSgrid-1
                        if(y(ich,mif) > x(j) .and. y(ich,mif) < x(j+1)) then
                        teff=(x0(j)*(y(ich,mif)-x(j))+x0(j+1)*(x(j+1)-y(ich,mif))) &
                          & /(x(j+1)-x(j))
                        end if
                     end do

                     ! solve the relation one more time to refine teff and sensitivity
                     sensitivity = slope%values(k+nFreqs*(mif-1),fmStat%maf)* &
                           & (1._r8 + 0.46_r8* teff**6) ! correction term (see ATBD)

                     ! in case we go into ambiguity altitudes with small sensitivity
                     if(abs(sensitivity) < 1._r8) sensitivity = 1._r8
                     
                     teff = y(ich,mif)/sensitivity
                     if(teff > 1._r8) teff = 1._r8
                           
                     ! convert cloud radiance to effective optical depth
                     
                     y(ich,mif) = teff
                     sy(ich,mif) = Terr%values(k+nFreqs*(mif-1),fmStat%maf)**2 &
                       & /sensitivity**2
                     
                     do i=1,nSgrid
                        dx(i,mif) = dx(i,mif) + jBlock%values(k,i+(mif-1)*nSgrid)* &
                           & y(ich,mif)/sy(ich,mif)
                     do j=1,nSgrid
                        C(i,j,mif) = C(i,j,mif) + jBlock%values(k,i+(mif-1)*nSgrid)* &
                           & jBlock%values(k,j+(mif-1)*nSgrid)/sy(ich,mif)
                           
                     end do
                     end do
                  end if
                  end do   ! mif

               end if
               end do      ! end of frequency k

            end do         ! end of band signals
         end do         ! end of imodel

         ! check if Jacobian rows are consistent with Signal rows
         if(ich /= nChans) print*,'inconsistent channels between Jacobian and Signal'

           sx = 1.e8_r8         ! sx is the inversd variance of a priori
           do mif=1,nMifs
!            x0 = sqrt(re%values(1,fmStat%maf)**2+(sLevel*1.e3_r8)**2) - &
!               re%values(1,fmStat%maf)
             do i=1,nSgrid
!               if(x0(i) < 20.e3_r8) &
                  sx(i,mif) = xLosVar%values(i+nSgrid*(mif-1),fmStat%maf)
             end do
           end do
                              
         ! start inversion
            do mif=1,nMifs
            if (ptan%values(mif,fmStat%maf) < p_lowcut) then
               ! for cloud retrieval: x_star=0, y_star=0, x0=apriori=0

               x0 = 0._r8        ! A priori
               x = 0._r8
               sx0=sx(:,mif)
               do j=1,maxJacobians     ! we use maxJacobians as number of iterations (default 5)
               
                  x0 = x        ! use  the last retrieval as A priori
                                 ! and doubt constraints to the A priori
                  where(x0 .le. 0._r8) x0 = 0._r8

                  A = reshape(C(:,:,mif),(/nSgrid,nSgrid/))
                  do i=1,nSgrid
                   A(i,i)=A(i,i) + sx0(i)       ! sx has already been inverted
                  end do
                  Call MatrixInversion(A)
               
                  ! output estimated SD after the first iteration
                  if(associated(xLosPtr) .and. j == 1) then
                   do i=1,nSgrid
                    outLosSD%values(i+(mif-1)*nSgrid,fmStat%maf) = sqrt(A(i,i))
                   end do
                  end if
                  
                  ! compute the LOS retrieval in the least-squard solution
                  do i=1,nSgrid
                     x(i) = sum( reshape(A(i,:),(/nSgrid/))* &
                     & (x0*sx0 + reshape(dx(:,mif),(/nSgrid/)) ))
                  end do
                  ! reduce covariance on x0 for next iteration, using factor lambda (default 10)
                  sx0 = sx0*initLambda          ! sx0 is the inversed variance

               end do  ! end of iteration
               
               ! Now, x is the los Transmission increment
               ! compute cloud extinction from los Transmission increments
               xext = 0._r8      ! temporary storage for LOS extinction
               trans = 0._r8
               do i=nSgrid,2,-1
                 trans = trans + x(i)
                 ! protect from blowing up
                 if(1._r8-trans > 0.004) then
                  scale = (1._r8-trans)*(slevel(i)-slevel(i-1))      ! see ATBD
			         xext(i)=x(i)/scale
                  ! and standard deviation
                  outLosSD%values(i+(mif-1)*nSgrid,fmStat%maf) = &     !neglect higher orders
                     & outLosSD%values(i+(mif-1)*nSgrid,fmStat%maf)/scale
                 end if
               end do
               
               ! output los extinction to state vector
               if(associated(xLosPtr)) then
                do i=1,nSgrid
                 xLosPtr%values(i+(mif-1)*nSgrid,fmStat%maf) = xext(i)
                end do
               end if               
               
            end if
            end do  ! mif
      call clearMatrix ( jacobian )           ! free the space
      end do ! end of mafs
      
      ! deallocate arrays and free memory
            deallocate(A,C,y,sy,dx,x,x0,sx0,xext,sx)

      xExtPtr%values = 0._r8
   	call LOS2Grid(xExtPtr,xLosPtr,Ptan,Re,p_lowcut,method='Average')
      
      if(associated(xLosVar) .and. associated(xExtVar) &
         & .and. associated(outLosSD) .and. associated(outExtSD)) then
      ! get variances after re-sampling
   	   call LOS2Grid(xExtVar,xLosVar,Ptan,Re,p_lowcut,method='Variance')
      ! get output variances after re-sampling, outLosSD is SD at present
         outLosSD%values = 1._r8/outLosSD%values**2 ! converted to 1/variance
         outExtSD%values = xExtVar%values      ! initialized to inversed a priori variance
         call LOS2Grid(outExtSD,outLosSD,Ptan,Re,p_lowcut,method='Variance')
      ! output SD of the retrieved extinction
         outLosSD%values = sqrt(1._r8/outLosSD%values) ! converted back to SD
         outExtSD%values = sqrt(1._r8/outExtSD%values) ! converted back to SD
      ! if the input variance is NOT reduced by half, it is given negative
         where(outExtSD%values**2 > 0.5/xExtVar%values) &   ! xExtVar is inversed
            & outExtSD%values = -outExtSD%values
      end if
   ! clean up
      call destroyVectorInfo ( FwdModelOut1 )
      call destroyVectorInfo ( CovarianceDiag )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      deallocate(Slevel)
      
    end subroutine LowCloudRetrieval

  !=============================== LOS2Grid ====
  subroutine LOS2Grid( Qty, LOS,Ptan, Re, Pcut, method)

    ! This is to fill a l2gp type of quantity with a los grid type of quantity.
    ! The los quantity is a vector quantity that has dimension of (s, mif, maf),
    ! where s is the path along los.
    !
    ! Linear interpolation is used to fill l2gp grids and unfilled grids are
    ! marked with the baddata flag (-999.)

    use UNITS
    use MLSNumerics, only: InterpolateValues
    use MLSStrings, only: Capitalize
    ! Dummy arguments
    character (len=*), intent(in) :: method ! whether for averaging or variance
    type (VectorValue_T), intent(in) :: LOS ! Vector quantity to fill from
    type (VectorValue_T), intent(in) :: Ptan ! tangent pressure
    type (VectorValue_T), intent(in) :: Re ! Earth's radius
    type (VectorValue_T), INTENT(INOUT) :: QTY ! Quantity to fill

    ! Local variables
    integer :: i, j, maf, mif                ! Loop counter
    integer :: maxZ, minZ                    ! pressure range indices of sGrid
    integer :: noMAFs,noMIFs,noDepths
    real (r8) :: pcut
    integer, dimension(qty%template%noSurfs,qty%template%noInstances) :: cnt
    real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: out
    real (r8), dimension(qty%template%noSurfs) :: outZeta, phi_out, beta_out
    real (r8), dimension(los%template%noChans) :: x_in, y_in, sLevel
    real (r8), dimension(los%template%noSurfs) :: zt
    real (r8), dimension(los%template%noChans,los%template%noSurfs,los%template%noInstances) :: beta

    if ( qty%template%verticalCoordinate == l_pressure ) then
        outZeta = -log10 ( qty%template%surfs(:,1) )
    else
        outZeta = qty%template%surfs(:,1)
    end if

    noMAFs=los%template%noInstances
    noMIFs=los%template%noSurfs
! Now, we use frequency coordinate as sGrid along the path
    noDepths=los%template%noChans
    sLevel = los%template%frequencies
   
   do maf=1,noMafs
   do mif=1,noMifs
      do i=1,noDepths
         beta(i,mif,maf)=los%values(i+(mif-1)*noDepths,maf)
      end do
   end do
   end do
      
! initialize quantity
   do j = 1, qty%template%noInstances
   do i = 1, qty%template%noSurfs
!   qty%values(i,j)=qty%template%badValue
   cnt(i,j)=0
   out(i,j)=0._r8
   end do 
   end do
   
    do maf=1,noMAFs  
      zt = ptan%values(:,maf)   ! noChans=1 for ptan
      zt = (zt+3.)*16.                      ! converted to height in km
      do mif=1,noMIFs
      if (ptan%values(mif,maf) .gt. pcut) cycle
      ! find altitude of each s grid
      x_in = sLevel**2/2./(re%values(1,maf)*0.001_r8 + zt(mif))
      ! converted to zeta
      x_in = x_in/16. + ptan%values(mif,maf)
      ! find minimum and maximum pressures indices in sGrid
        do i = 2,qty%template%noSurfs-1
        if (ptan%values(mif,maf) < (outZeta(i)+outZeta(i+1))/2. .and. &
          & ptan%values(mif,maf) > (outZeta(i)+outZeta(i-1))/2.) &
          & minZ = i
        end do
        if (ptan%values(mif,maf) < (outZeta(1)+outZeta(2))/2.) minZ=1
        if (ptan%values(mif,maf) > outZeta(qty%template%noSurfs)) cycle ! goto next mif
        
        do i = 2,qty%template%noSurfs-1
        if (x_in(noDepths) < (outZeta(i)+outZeta(i+1))/2. .and. &
          & x_in(noDepths) > (outZeta(i)+outZeta(i-1))/2.) &
          & maxZ = i
        end do
        if (x_in(noDepths) < (outZeta(1)+outZeta(2))/2.) cycle    ! goto next mif
        if (x_in(noDepths) > outZeta(qty%template%noSurfs)) maxZ=qty%template%noSurfs

      ! get phi along path for each mif (phi is in degree)
      y_in = los%template%phi(mif,maf) &
        & - atan(sLevel/(re%values(1,maf)*0.001_r8 + zt(mif)))*180._r8/Pi
        ! interpolate phi onto standard vertical grids     
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),phi_out(minZ:maxZ), &
           & method='Linear')
        ! interpolate quantity to standard vertical grids      
        y_in = beta(:,mif,maf)
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),beta_out(minZ:maxZ), &
           & method='Linear')
        ! interpolate quantity to standard phi grids
        do i=minZ,maxZ  
          do j = 2, qty%template%noInstances-1
          if(phi_out(i) .lt. &     
            & (qty%template%phi(1,j)+qty%template%phi(1,j+1))/2. &
            & .and. phi_out(i) .ge. &  
            & (qty%template%phi(1,j-1)+qty%template%phi(1,j))/2. ) then
            out(i,j)=out(i,j) + beta_out(i)
            cnt(i,j)=cnt(i,j)+1       !  counter
          end if
          end do
        end do
      end do                            ! End surface loop
    end do                              ! End instance loop
    ! average all non-zero bins
    if(Capitalize(method(1:1))=="A") where (cnt > 0) qty%values = out/cnt
    if(Capitalize(method(1:1))=="V") where (cnt > 0) qty%values = out

  end subroutine LOS2Grid
    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call cpu_time ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( "Timing for Retrieve = " )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

    ! -------------------------------------------------- SetupSubset ---
    subroutine SetupSubset ( key, vectors )

      use Declaration_table, only: NUM_VALUE, RANGE
      use Intrinsic, only: PHYQ_LENGTH, PHYQ_PRESSURE
      use VectorsModule, only: ClearMask, CreateMask, &
        & GetVectorQtyByTemplateIndex, SetMask, VectorValue_T

      integer, intent(in) :: KEY        ! Tree node
      type (Vector_T), dimension(:) :: VECTORS

      ! Local variables
      integer :: CHANNEL                ! Loop index
      integer :: CHANNELSNODE           ! Tree node for channels values
      integer :: COORDINATE             ! Vertical coordinate type
      integer :: DEPTHNODE              ! Tree node for optical depth
      integer :: FIELD                  ! Field type from tree
      integer :: GSON                   ! Tree node
      integer :: HEIGHT                 ! Loop counter
      integer :: HEIGHTNODE             ! Tree node for height values
      integer :: HEIGHTUNIT             ! Unit for heights command
      integer :: I, J                   ! Subscripts, loop inductors
      integer :: INSTANCE               ! Loop counter
      integer :: QUANTITYINDEX          ! Index
      integer :: NROWS                  ! Loop limit dumping mask
      integer :: ROW                    ! Row index dumping mask
      integer :: SON                    ! Tree node
      integer :: TYPE                   ! Type of value returned by expr
      integer :: UNITS(2)               ! Units returned by expr
      integer :: VECTORINDEX            ! Index

      integer :: S1(1), S2(1)           ! Results of minloc intrinsic

      real(r8) :: VALUE(2)              ! Value returned by expr
      real(r8), dimension(:), pointer :: THESEHEIGHTS ! Subset of heights
      type (VectorValue_T), pointer :: QTY ! The quantity to mask
      type (VectorValue_T), pointer :: PTAN ! The ptan quantity if needed
      logical :: Got(field_first:field_last)   ! "Got this field already"
      logical, dimension(:), pointer :: CHANNELS ! Are we dealing with these channels
      logical :: IGNORE                 ! Flag
      logical :: DOTHIS                 ! Flag
      integer, parameter ::                      MAXCOLUMNS = 127
      character(len=1), dimension(MAXCOLUMNS) :: maskedMan
      character(len=5) ::                        decades
      real(r8)         ::                        heightMin, heightMax

      ! Executable code
      nullify ( channels, qty, ptan )
      got = .false.
      ignore = .false.
      do j = 2, nsons(key) ! fields of the "subset" specification
        son = subtree(j, key)
        field = get_field_id(son)   ! tree_checker prevents duplicates
        got(field) = .true.
        if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
        select case ( field )
        case ( f_quantity )
          vectorIndex = decoration(decoration(subtree(1,gson)))
          quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          qty => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
        case ( f_ptanquantity )
          vectorIndex = decoration(decoration(subtree(1,gson)))
          quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          ptan => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
        case ( f_channels )
          channelsNode = son
        case ( f_height )
          heightNode = son
        case ( f_opticalDepth )
          depthNode = son
        case ( f_ignore )
          ignore = Get_Boolean ( son )
        case default
          ! Shouldn't get here if the type checker worked
        end select
      end do ! j = 2, nsons(key)

      ! Process the channels field.
      if ( qty%template%frequencyCoordinate /= l_none ) then
        call Allocate_test ( channels, qty%template%noChans, &
          & 'channels', ModuleName )
        if ( got(f_channels) ) then     ! This subset is only for some channels
          channels = .false.
          do j = 2, nsons(channelsNode)
            call expr ( subtree(j,channelsNode), units, value, type )
            do i = 1, merge(1,2,type==num_value)
              if ( units(i) /= phyq_dimensionless ) &
                & call announceError ( wrongUnits, f_channels, string='no' )
            end do
            channels ( nint(value(1)) : &
              & nint(value(merge(1,2,type==num_value))) ) = .true.
          end do
        else
          channels = .true.             ! Apply this to all channels
        end if
      end if

      ! Preprocess the height stuff.  
      heightUnit = phyq_dimensionless
      if ( got(f_height) ) then
        if ( ignore ) call announceError ( inconsistent, f_height, f_ignore )
        do j = 2, nsons(heightNode)
          call expr ( subtree(j,heightNode), units, value, type )
          ! Make sure the range has non-dimensionless units -- the type
          ! checker only verifies that they're consistent.  We need to
          ! check each range separately, because the units determine the
          ! scaling of the values.
          if ( all(units == phyq_dimensionless) ) call announceError ( &
            & wrongUnits, f_height, string = 'length or pressure.' )
          ! Check consistency of units -- all the same, or dimensionless. The
          ! type checker verifies the consistency of units of ranges, but not
          ! of array elements.
          do i = 1, 2
            if ( heightUnit == phyq_dimensionless ) then
              heightUnit = units(i)
            else if ( units(i) /= phyq_dimensionless .and. &
              &       units(i) /= heightUnit ) then
              call announceError ( inconsistentUnits, f_height )
            end if
          end do
        end do
        ! Check for correct units
        if ( heightUnit == phyq_pressure ) then
          if ( qty%template%minorFrame .and. .not. got(f_ptanQuantity) ) &
            & call announceError ( ifUnitsAThenB, f_height, f_ptanQuantity, &
            & 'pressure' )
        else if ( heightUnit /= phyq_length ) then
          call announceError ( wrongUnits, f_height, &
            & string = 'length or pressure.' )
        end if
      end if

      ! Create the mask if it doesn't exist
      if ( .not. associated( qty%mask ) ) call CreateMask ( qty )

      ! Now we loop over the instances
      do instance = 1, qty%template%noInstances
        if ( qty%template%coherent ) then
          theseHeights => qty%template%surfs(:,1)
          coordinate = qty%template%verticalCoordinate
        else if ( qty%template%minorFrame .and. heightUnit == phyq_pressure ) then
          theseHeights => ptan%values(:,instance)
          coordinate = l_zeta
        else
          theseHeights => qty%template%surfs(:,instance)
          coordinate = qty%template%verticalCoordinate
        end if

        ! Now, make sure for the channels we're considering that the
        ! default is to ignore all
        if ( got(f_height) .or. ignore ) then
          do channel = 1, qty%template%noChans
            doThis = .true.
            if ( associated(channels) ) doThis = channels(channel)
            if ( doThis ) then
              do height = 1, qty%template%noSurfs
                !??? Make sure mask bit numbers begin at 1, even when
                !??? channel numbers don't.
                call SetMask ( qty%mask(:,instance), &
                  & (/ channel+qty%template%noChans*(height-1) /) )
              end do                    ! Height loop
            end if                      ! Do this channel
          end do                        ! Channel loop
        end if                          ! Got heights or ignore

        ! Now go and `unmask' the ones we want to consider
        if ( got(f_height) ) then
          do j = 2, nsons(heightNode)
            call expr ( subtree(j,heightNode), units, value, type )
            ! Now maybe do something nasty to value to get in right units
            if ( coordinate == l_zeta &
              & .and. heightUnit == phyq_pressure ) then
              value = -log10(value)
!             else if ( coordinate /= qty%template%verticalCoordinate ) then
!               call MLSMessage ( MLSMSG_Error, ModuleName, &
!                 & 'Inappropriate units for height in subset' )
            end if
            if ( qty%template%verticalCoordinate == l_pressure ) then
              s1 = minloc ( abs ( -log10(theseHeights) + log10(value(1)) ) )
              s2 = minloc ( abs ( -log10(theseHeights) + log10(value(2)) ) )
            else
              s1 = minloc ( abs ( theseHeights - value(1) ) )
              s2 = minloc ( abs ( theseHeights - value(2) ) )
            end if
            do channel = 1, qty%template%noChans
              doThis = .true.
              if ( associated(channels) ) doThis = channels(channel)
              if ( doThis ) then
                do height = s1(1), s2(1)
                  !??? Make sure mask bit numbers begin at 1, even when
                  !??? channel numbers don't.
                  call ClearMask ( qty%mask(:,instance), &
                    & (/ channel+qty%template%noChans*(height-1) /) )
                end do                  ! Height loop
              end if                    ! Do this channel
            end do                      ! Channel loop
          end do                        ! Height entries in l2cf
        end if                          ! Got a height entry
      end do                            ! Instance loop

      ! Tidy up
      call Deallocate_test ( channels, 'channels', ModuleName )
      
      if ( index(switches,'msk') /= 0 ) then
        if ( qty%template%name /= 0 ) then
          call output ( ' Qty_Template_Name = ' )
          call display_string ( qty%template%name )
        end if
        call output ( ' ', advance='yes' )
        call output ( 'Elements per mask = ' )
        call output ( size(qty%values,1), advance='yes' )
        call output ( 'noChans = ' )
        call output ( qty%template%noChans, advance='no' )
        call output ( ' noSurfs = ' )
        call output ( qty%template%noSurfs, advance='no' )
        call output ( ' noInstances = ' )
        call output ( qty%template%noInstances, advance='no' )
        call output ( ' instanceLen = ' )
        call output ( qty%template%instanceLen, advance='yes' )
        call dump ( qty%mask, format='(1x,z8.8)' )
        ! P. Wagner's inspection of masking array
        if ( got(f_height) ) then
          nrows = qty%template%noSurfs
          call output ( '(Suppressing channel info, rows are surfaces) ', &
          & advance='no' )
        else
          nrows = qty%template%noChans
          call output ( '(Suppressing mif info, rows are channels) ', &
          & advance='no' )
        endif
        nrows = min(nrows, MAXCOLUMNS)
        call output ( 'column, row look at mask', advance='yes' )
        maskedMan=' '
        call output ( 'column ', advance='no' )
        call blanks(14, advance='no')
        call output ( 'mask', advance='yes' )
        call output ( 'rows ->', advance='no' )
        do j=10, nrows, 10
          call blanks(10-len(decades), advance='no')
          write(decades, '(i5.2)') mod(j, 100)
          call output(decades, advance='no')
        enddo
        call output(' ', advance='yes')
        heightMin=+1.d4
        heightMax=-1.d4
        do i=1, size(qty%values(1,:))
          maskedMan='0'
          do j=1, nrows
            if ( got(f_height) ) then
              row = 1 + qty%template%noChans*(j-1)
            else
              row = j
            endif
            if ( IsVectorQtyMasked(qty, row, i) ) then
              maskedMan(j)='1'
              if ( got(f_height) .and. j <= size(theseheights)) then
                heightMin = min(heightMin, theseheights(j))
                heightMax = max(heightMax, theseheights(j))
              endif
            endif
          enddo
          call output ( i, format='(i3)', advance='no' )
          call blanks(4, advance='no')
          call output ( MaskedMan(1:nrows), advance='yes' )
        enddo
        if ( got(f_height) ) then
          call output ( 'min, max heights masked out: ', advance='no' )
          call output ( exp(-log(10.)*heightMax), advance='no' )
          call blanks(4, advance='no')
          call output ( exp(-log(10.)*heightMin), advance='yes' )
        endif
      end if
    end subroutine SetupSubset
  end subroutine Retrieve

end module RetrievalModule

! $Log$
! Revision 2.101  2001/10/22 20:52:00  dwu
! renaming some cloud functions
!
! Revision 2.100  2001/10/20 17:50:48  livesey
! A little more information on aj%fnmin < 0.0
!
! Revision 2.99  2001/10/19 17:42:36  dwu
! add protection for cloud sensitivity from going to zero
!
! Revision 2.98  2001/10/18 23:25:13  dwu
! rename FillQuantityFromLos to BinFromLOS2Grid
!
! Revision 2.97  2001/10/17 23:37:29  pwagner
! Improved dump of mask
!
! Revision 2.96  2001/10/17 20:50:30  dwu
! a fix of cloud retrieval
!
! Revision 2.95  2001/10/17 19:52:43  dwu
! Add standard deviation calculation to LowCloudRetrieval
!
! Revision 2.94  2001/10/16 23:35:39  pwagner
! More dumped on msk switch
!
! Revision 2.93  2001/10/15 23:21:47  vsnyder
! Forgot to clone AprioriMinusX after changing multiply not to clone
!
! Revision 2.92  2001/10/15 22:41:11  vsnyder
! Put the local (non-l2cf) vectors into a private database of fixed size,
! for snooping.  We can't copy them in, because it's a shallow copy, and
! we change some pointers later.
!
! Revision 2.91  2001/10/11 16:28:25  dwu
! Add cloudretrieval method for low tangent height radiances
!
! Revision 2.90  2001/10/09 20:39:36  vsnyder
! Left out a blank on the CVS command line, botched the last comment
!
! Revision 2.89  2001/10/09 20:38:23  vsnyder
! Corrections for regularization; output f, not f-measurements
!
! Revision 2.88  2001/10/05 20:50:16  vsnyder
! Concatenate Snoop comment and DNWT flag; clear F before EVAL[FJ]
!
! Revision 2.87  2001/10/05 20:20:16  vsnyder
! Put Nwt_Flag into Snoop's comment; add Dnwt_Flag to diagnostics vector
!
! Revision 2.86  2001/10/05 05:02:27  dwu
! temporarily set p_lowcut in LowCloud Retrieval
!
! Revision 2.85  2001/10/05 01:46:25  vsnyder
! Put in stuff to snoop in the Newtonian iteration loop
!
! Revision 2.84  2001/10/04 01:48:59  vsnyder
! Move temporary vectors into 'myVectors' database, for snooping
!
! Revision 2.83  2001/10/04 00:30:34  dwu
! fix coljBlock finding in Low Cloud
!
! Revision 2.82  2001/10/03 23:56:30  dwu
! some minor fixes for Low Cloud
!
! Revision 2.81  2001/10/03 22:05:12  dwu
! some quick remedies for LowcloudRetrieval
!
! Revision 2.80  2001/10/03 21:49:54  dwu
! some quick remedies for LowcloudRetrieval
!
! Revision 2.79  2001/10/03 21:30:16  dwu
! add LowCloudRetrieval
!
! Revision 2.78  2001/10/03 17:57:21  vsnyder
! Delete l_DegreesOfFreedom, add l-DNWT_...
!
! Revision 2.77  2001/10/02 23:41:11  vsnyder
! Added computation of degrees of freedom, diagnostics vector
!
! Revision 2.76  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models
!
! Revision 2.75  2001/10/01 23:30:50  pwagner
! Fixed bug in spelling cholesky_solver
!
! Revision 2.74  2001/10/01 23:04:17  livesey
! Bug fix with channel range, and some tidying up
!
! Revision 2.73  2001/10/01 22:54:22  pwagner
! Added subsection timings for Retrieval section
!
! Revision 2.72  2001/10/01 20:30:48  vsnyder
! Insert a reminder for necessary future work
!
! Revision 2.71  2001/09/29 00:07:11  pwagner
! Fixed various timing problems
!
! Revision 2.70  2001/09/28 18:25:37  dwu
! prepare for adding lowcloud retrieval method
!
! Revision 2.69  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.68  2001/09/27 20:14:50  vsnyder
! Add 'msk' switch to control printing the mask after a Subset
!
! Revision 2.67  2001/09/27 18:40:26  vsnyder
! Explicitly control mask use for FormNormalEquations
!
! Revision 2.66  2001/09/26 02:15:40  vsnyder
! More work on checking the Subset command
!
! Revision 2.65  2001/09/25 23:04:30  livesey
! Fixed uninitialised heightUnit
!
! Revision 2.64  2001/09/25 20:42:55  vsnyder
! Spiffify error processing for Subset
!
! Revision 2.63  2001/09/25 17:49:59  livesey
! More updates, and fixes to subset
!
! Revision 2.62  2001/09/25 05:57:36  livesey
! Removed call to destroyVectorMask, and the erroneous .not. doThis
!
! Revision 2.61  2001/09/25 00:49:59  vsnyder
! Copy mask from measurements to f_rowScaled
!
! Revision 2.60  2001/09/25 00:18:36  livesey
! Bug fix in Subset
!
! Revision 2.59  2001/09/20 00:31:06  vsnyder
! Move deallocation of 'channels' out of loop
!
! Revision 2.58  2001/07/19 22:00:41  vsnyder
! orrect problems with row scaling.
!
! Revision 2.57  2001/07/12 22:18:05  livesey
! Got rid of an old diagnostic
!
! Revision 2.56  2001/07/12 22:11:46  vsnyder
! Maybe the column scaling is right now....
!
! Revision 2.55  2001/07/11 22:06:31  vsnyder
! Interim commit -- still appears to be broken
!
! Revision 2.54  2001/07/02 17:21:48  livesey
! Fixed memory leak with channels.
!
! Revision 2.53  2001/06/30 03:10:35  vsnyder
! Don't access sub_rosa(0) when creating a covariance matrix
!
! Revision 2.52  2001/06/30 03:07:43  vsnyder
! Remove some unused variables.  Move some more internal-subroutine-peculiar
! USEs into them.  Don't access sub_rosa(0) when creating a Jacobian.
!
! Revision 2.51  2001/06/30 02:35:43  vsnyder
! Don't use the "extra" column to solve for the residual.  Move the Newtonian
! solver, and its variables and USEs, into an internal subroutine.  Move USEs
! peculiar to SetupSubset into that routine.
!
! Revision 2.50  2001/06/27 04:31:29  livesey
! Subset is getting closer bit by bit.
!
! Revision 2.49  2001/06/27 04:05:32  livesey
! Interim version with known problems in subset
!
! Revision 2.48  2001/06/26 20:11:32  livesey
! Bug fixes to subset (more to come I imagine)
!
! Revision 2.47  2001/06/26 19:01:00  vsnyder
! Specify regularization orders according to quantities
!
! Revision 2.46  2001/06/26 18:18:14  livesey
! Another (working?) version.
!
! Revision 2.45  2001/06/26 17:57:46  livesey
! Whoops, added another use clause
!
! Revision 2.44  2001/06/26 16:26:32  livesey
! It at least compiles, but subset certainly doesn't work yet.
!
! Revision 2.43  2001/06/22 01:26:54  vsnyder
! Process fields for regularization, but regularization isn't done yet
!
! Revision 2.42  2001/06/20 22:21:22  vsnyder
! Added initialization for fmStat%newScanHydros
!
! Revision 2.41  2001/06/20 21:44:33  vsnyder
! Don't compute cosines between zero-length vectors
!
! Revision 2.40  2001/06/01 22:26:03  vsnyder
! Got the row scaling backwards in the NF_EVALF
!
! Revision 2.39  2001/06/01 21:58:17  livesey
! Added scale for outputSD.
!
! Revision 2.38  2001/06/01 21:40:18  vsnyder
! Row scale during NF_EVALF; destroy some vectors so as not to have a
! memory leak; initially clone x to get columnScaleVector.
!
! Revision 2.37  2001/06/01 21:27:48  livesey
! Added outSD field
!
! Revision 2.36  2001/06/01 20:36:34  vsnyder
! Seems to work.  Added LaTeX comments.
!
! Revision 2.35  2001/06/01 01:02:32  vsnyder
! Periodic commit.  Not working right yet.
!
! Revision 2.34  2001/05/24 23:28:45  vsnyder
! Scale things back to user coordinates at the correct times; numerous output changes
!
! Revision 2.33  2001/05/22 19:10:59  vsnyder
! Periodic commit; still isn't accepting a Newton step
!
! Revision 2.32  2001/05/19 01:17:39  vsnyder
! Correct sign of (apriori-x)
!
! Revision 2.31  2001/05/19 00:22:31  vsnyder
! Dump_L1 -> Dump_Linf, uncouple 'spa' and Linf dumps
!
! Revision 2.30  2001/05/18 23:18:42  vsnyder
! Replace 'weight' field of 'retrieve' by 'measurementSD'
!
! Revision 2.29  2001/05/18 19:46:23  vsnyder
! Add secret 'fuzz' field to 'retrieve' command -- for testing
!
! Revision 2.28  2001/05/18 01:04:08  vsnyder
! Periodic commit -- tons of stuff changed
!
! Revision 2.27  2001/05/10 22:50:45  vsnyder
! Added switches to print stuff during Newton iteration
!
! Revision 2.26  2001/05/03 02:00:15  vsnyder
! Put names on cloned vectors
!
! Revision 2.25  2001/05/02 05:28:04  livesey
! Added DumpBlocks
!
! Revision 2.24  2001/05/01 23:52:04  vsnyder
! Allocate and deallocate fmStat%rows here
!
! Revision 2.23  2001/04/28 01:47:17  vsnyder
! Don't try to create initial value for state
!
! Revision 2.22  2001/04/26 23:43:23  vsnyder
! Remove forwardModelIn from retrieve spec
!
! Revision 2.21  2001/04/26 19:48:20  livesey
! Now uses ForwardModelWrappers
!
! Revision 2.20  2001/04/26 02:53:37  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.19  2001/04/26 01:04:21  vsnyder
! Copy model's final radiance to fwdModelOut
!
! Revision 2.18  2001/04/26 01:00:01  vsnyder
! Get the Jacobian's row block from fmStat%rows, cosmetic changes
!
! Revision 2.17  2001/04/21 01:44:43  vsnyder
! Make the timing message prettier
!
! Revision 2.16  2001/04/19 23:56:23  livesey
! New fmStat
!
! Revision 2.15  2001/04/13 21:40:58  vsnyder
! Periodic commit -- stuff about looping over configs
!
! Revision 2.14  2001/04/12 01:50:21  vsnyder
! Work on getting a posteriori covariance
!
! Revision 2.13  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.12  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.11  2001/03/17 00:45:15  livesey
! Moved to new ForwardModelConfig_T
!
! Revision 2.10  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.9  2001/03/09 03:05:05  vsnyder
! Correct identification for timing
!
! Revision 2.8  2001/03/08 03:23:09  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.7  2001/03/07 23:59:52  vsnyder
! Add stuff for SIDS.
!
! Revision 2.6  2001/02/22 18:54:05  vsnyder
! Periodic commit.  Still working on the output covariance.
!
! Revision 2.5  2001/02/22 01:57:02  vsnyder
! Periodic commit -- working on getting output covariance matrix.
!
! Revision 2.4  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.3  2001/02/08 00:56:55  vsnyder
! Periodic commit.  Still needs work.
!
! Revision 2.2  2001/01/26 19:01:47  vsnyder
! More nearly complete, except for forward model interface and minor things
! having to do with creating subset masks.  Look for ??? in comments.
!
! Revision 2.1  2001/01/10 21:04:13  vsnyder
! Initial (incomplete) submission
!
