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
    & F_columnScale, F_covariance, F_diagnostics, F_diagonal, F_forwardModel, &
    & F_fuzz, F_fwdModelExtra, F_fwdModelOut, F_height, F_ignore, F_jacobian, &
    & F_lambda, F_maxF, F_maxJ, F_measurements, F_measurementSD, F_method, &
    & F_opticalDepth, F_outputCovariance, F_outputSD, F_ptanQuantity, &
    & F_quantity, F_regOrders, F_regQuants, F_regWeight, F_state, &
    & F_toleranceA, F_toleranceF, F_toleranceR, &
    & Field_first, Field_last, &
    & L_apriori, L_covariance, L_newtonian, L_lowcloud, &
    & L_none, L_norm, L_pressure, L_zeta, &
    & S_dumpBlocks, S_forwardModel, S_sids, S_matrix, S_subset, S_retrieve, &
    & S_time
  use Intrinsic, only: L_Jacobian_Cols, L_Jacobian_Rows, PHYQ_Dimensionless
  use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
    & DestroyMatrix, GetFromMatrixDatabase, Matrix_T, Matrix_Database_T, &
    & Matrix_SPD_T
  use MatrixTools, only: DumpBlock
  use MLSCommon, only: R8, MLSCHUNK_T
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, add_to_retrieval_timing
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
  use OUTPUT_M, only: BLANKS, OUTPUT
  use SidsModule, only: SIDS
  use Toggles, only: Gen, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: CloneVector, CopyVector, DestroyVectorInfo, &
    & DumpMask, GetVectorQuantityByType, Vector_T, VectorValue_T

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
  subroutine Retrieve ( Root, VectorDatabase, MatrixDatabase, ConfigDatabase, chunk )

  !{Process the "Retrieve" section of the L2 Configuration File.
  ! The "Retrieve" section can have ForwardModel, Matrix, Sids, Subset or
  ! Retrieve specifications.

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
                                        ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
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
    type(vectorValue_T), pointer :: DOF_Qty  ! Degrees of freedom quantity in
                                        ! Diagnostics vector
    integer :: Error
    type(vector_T) :: F                 ! Residual -- Model - Measurements
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
    type(vectorValue_T), pointer :: Jac_Col_Qty   ! To output Jacobian_Cols
    integer :: Jacobian_Rows            ! (Number of rows of Jacobian) -
                                        ! (masked-off rows of Jacobian)
    type(vectorValue_T), pointer :: Jac_Row_Qty   ! To output Jacobian_Rows
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
    type(matrix_SPD_T), pointer :: OutputCovariance   ! Covariance of the sol'n
    type(vector_T), pointer :: OutputSD ! Vector containing SD of result
    integer :: RegOrders                ! Regularization orders
    integer :: RegQuants                ! Regularization quantities
    real(r8) :: RegWeight               ! Weight of regularization conditions
    integer :: Son                      ! Of Root or Key
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
    nullify ( measurements, measurementSD, state, outputSD )
    timing = section_times
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
            call add_to_retrieval_timing( 'low_cloud', t1 )
            call cpu_time ( t1 )
          ! call LowCloudRetrieval
            print*,'to be added'
          end select ! method
          !??? Make sure the jacobian and outputCovariance get destroyed
          !??? after ?what? happens?  Can we destroy the entire matrix
          !??? database at the end of each chunk?
          if ( .not. got(f_jacobian) ) call destroyMatrix ( jacobian )
          if ( .not. got(f_outputCovariance) ) &
            & call destroyMatrix ( outputCovariance%m )
        end if
        if ( got(f_fwdModelOut) ) call copyVector ( fwdModelOut, f )
        call destroyVectorInfo ( f )
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
        call output ( ' appears, but the other is not.', advance='yes' )
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
      type(vector_T) :: AprioriMinusX   ! Apriori - X
      real(r8) :: AprioriNorm           ! apriori .dot. apriori
      type(vector_T) :: ATb             ! A^T b -- the RHS of the normal eqns.
      type(vector_T) :: BestGradient    ! for NWT
      type(vector_T) :: BestX           ! for NWT
      type(vector_T) :: CandidateDX     ! for NWT
      type(vector_T) :: ColumnScaleVector ! For column scaling by column norms
      real(r8) :: Cosine                ! Of an angle between two vectors
      type(vector_T) :: CovarianceDiag  ! Diagonal of apriori Covariance
      type(vector_T) :: CovarianceXApriori ! Covariance \times Apriori
      type(vector_T) :: DX              ! for NWT
      type(vector_T) :: DXUnScaled      ! for NWT
      type(matrix_Cholesky_T) :: Factored ! Cholesky-factored normal equations
      type (ForwardModelStatus_T) :: FmStat ! Status for forward model
      type (ForwardModelIntermediate_T) :: Fmw ! Work space for forward model
      type(vector_T) :: F_rowScaled     ! Either a copy of f, or f row-scaled
      type(vector_T) :: FuzzState       ! Random numbers to fuzz the state
      type(vector_T) :: Gradient        ! for NWT
      integer :: J, K                   ! Loop inductors and subscripts
      type(matrix_SPD_T) :: NormalEquations  ! Jacobian**T * Jacobian
      integer :: NumF, NumJ             ! Number of Function, Jacobian evaluations
      integer :: NWT_Flag               ! Signal from NWT, q.v., indicating
                                        ! the action to take.
      integer :: NWT_Opt(20)            ! Options for NWT, q.v.
      real(rk) :: NWT_Xopt(20)          ! Real parameters for NWT options, q.v.
      type(vector_T) :: Reg_X_X         ! Regularization * X_n
      integer :: RowBlock               ! Which block of rows is the forward
                                        ! model filling?
      real :: T1
      character(len=10) :: TheFlagName  ! Name of NWTA's flag argument
      type(vector_T) :: Weight          ! Scaling vector for rows, 1/measurementSD
      type(vector_T) :: X               ! for NWT

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
      call cloneVector ( x, state, vectorNameText='_x' )
      call copyVector ( x, state ) ! x := state
      call cloneVector ( atb, state, vectorNameText='_ATb' )
      call cloneVector ( reg_X_x, state, vectorNameText='_reg_X_x' )
      call cloneVector ( f, measurements, vectorNameText='_f' )
      call cloneVector ( f_rowScaled, measurements, vectorNameText='_f_rowScaled' )
      call copyVector ( f_rowScaled, measurements, noValues=.true. ) ! mask only
      call cloneVector ( bestGradient, x, vectorNameText='_bestGradient' )
      call cloneVector ( bestX, x, vectorNameText='_bestX' )
      call cloneVector ( candidateDX, x, vectorNameText='_candidateDX' )
      call cloneVector ( dx, x, vectorNameText='_DX' )
      call cloneVector ( dxUnScaled, x, vectorNameText='_DX Unscaled' )
      call cloneVector ( gradient, x, vectorNameText='_gradient' )
      if ( got(f_measurementSD) ) then
        call cloneVector ( weight, measurementSD )
        do j = 1, measurementSD%template%noQuantities
          where ( measurementSD%quantities(j)%values <= 0.0 )
            weight%quantities(j)%values = 1.0
          elsewhere
            weight%quantities(j)%values = 1.0 / &
              & measurementSD%quantities(j)%values
          end where
        end do
      end if
      if ( columnScaling /= l_none ) &
        call cloneVector ( columnScaleVector, x, vectorNameText='_ColumnScale' )
      if ( got(f_fuzz) ) then
        ! Add some fuzz to the state vector (for testing purposes):
        call cloneVector ( fuzzState, x )
        do j = 1, x%template%noQuantities
          call random_number(fuzzState%quantities(j)%values)
          x%quantities(j)%values = x%quantities(j)%values * &
            & ( 1.0_r8 + fuzz * ( fuzzState%quantities(j)%values - 0.5 ) )
        end do
      end if
        if ( index(switches,'xvec') /= 0 ) call dump ( x, name='Original X' )
      numF = 0
      numJ = 0
      aj%axmax = 0.0
      do k = 1, size(x%quantities)
        aj%axmax = max(aj%axmax, maxval(abs(x%quantities(k)%values)))
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
        end if

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
            call cloneVector ( aprioriMinusX, apriori, &
              & vectorNameText='_aprioriMinusX' )
            aprioriMinusX = apriori - x
            if ( got(f_aprioriScale) ) &
              & call scaleVector ( aprioriMinusX, aprioriScale )
            !{Let the covariance of the apriori be $\mathbf{S_a}$, let
            ! $\mathbf{C = S_a^{-1}}$, and let $\mathbf{F^T F = C}$. In
            ! the least-squares problem we have extra rows of the form
            ! $\mathbf{F \delta x = F ( a - x_n )}$ where $\mathbf{a}$
            ! is the apriori state, and $\mathbf{x_n}$ is the state as
            ! of the previous iteration.  The initial bit of the
            ! right-hand side of the normal equations is $\mathbf{F^T F
            ! ( a - x_n ) = C ( a - x_n )}$.
            if ( diagonal ) then
              call getDiagonal ( covariance%m, covarianceDiag )
              ! covarianceXApriori := covarianceDiag # apriori:
              call multiply ( covarianceDiag, aprioriMinusX, &
                & covarianceXApriori )
            else ! covarianceXApriori := covariance X apriori:
              call multiply ( covariance, aprioriMinusX, &
                & covarianceXApriori )
            end if
            !{The contribution to the norm of the residual of the
            ! right-hand side of the part of the least-squares problem
            ! that is due to apriori is $\mathbf{
            ! ( a - x_n )^T F^T F ( a - x_n ) = ( a - x_n )^T C ( a - x_n )}$
            aprioriNorm = aprioriMinusX .dot. covarianceXapriori
          else
            aprioriNorm = 0.0_r8
          end if

          ! Compute f(x)
 
          fmStat%newScanHydros = .true.
          fmStat%maf = 0

          call add_to_retrieval_timing( 'newton_solver', t1 )
          call cpu_time ( t1 )
          ! Loop over MAFs
          do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
            fmStat%maf = fmStat%maf + 1
            do k = 1, size(configIndices)
              call forwardModel ( configDatabase(configIndices(k)), &
                & x, fwdModelExtra, f, fmw, fmStat )
            end do ! k
          end do ! MAFs
          call add_to_retrieval_timing( 'forward_model', t1 )
          call cpu_time ( t1 )
          call subtractFromVector ( f, measurements )
          if ( got(f_measurementSD) ) call multiply ( f, weight )
          aj%fnorm = sqrt ( aprioriNorm + ( f .dot. f ) )
            if ( index(switches,'fvec') /= 0 ) &
              & call dump ( f, name='Residual' )
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
              call updateDiagonal ( normalEquations, covarianceDiag )
            else
              call copyMatrixValue ( normalEquations%m, covariance%m )
            end if
            call copyVector ( atb, covarianceXApriori ) ! A^T b := C (a-x_n)
            if ( got(f_aprioriScale) ) &
              & call scaleMatrix ( normalEquations%m, aprioriScale )
          else
            call clearMatrix ( normalEquations%m ) ! start with zero
            call destroyVectorValue ( atb ) ! Clear the RHS vector
          end if
          ! aprioriNorm was computed when nwt_flag was NF_EVALF
          aj%fnorm = aprioriNorm

          ! Add Tikhonov regularization if requested
          if ( got(f_regOrders) ) then
            call regularize ( jacobian, regOrders, regQuants, regWeight )
            call multiply ( jacobian, x, reg_X_x ) ! regularization * x_n
            call scaleVector ( reg_X_x, -regWeight )   ! -R x_n
            call formNormalEquations ( jacobian, normalEquations, &
              & reg_X_x, atb, update=update, useMask=.false. )
            update = .true.
            call clearMatrix ( jacobian )        ! free the space
            aj%fnorm = aj%fnorm + ( reg_X_x .dot. reg_X_x )
            call destroyVectorValue ( reg_X_x )  ! free the space
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
          do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
            call add_to_retrieval_timing( 'newton_solver', t1 )
            call cpu_time ( t1 )
            ! What if one config set finished but others still had more
            ! to do? Ermmm, think of this next time.
            fmStat%maf = fmStat%maf + 1
            fmstat%rows = .false.
            do k = 1, size(configIndices)
              call forwardModel ( configDatabase(configIndices(k)), &
                & x, fwdModelExtra, f_rowScaled, fmw, fmStat, jacobian )
            end do ! k
            call add_to_retrieval_timing( 'forward_model', t1 )
            call cpu_time ( t1 )
            do rowBlock = 1, size(fmStat%rows)
              if ( fmStat%rows(rowBlock) ) then
                call subtractFromVector ( f_rowScaled, measurements, &
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
                  call multiply ( f_rowScaled, weight, &
                    & quant=jacobian%row%quant(rowBlock), &
                    & inst=jacobian%row%inst(rowBlock) )
                end if
              end if
            end do
            call scaleVector ( f_rowScaled, -1.0_r8 ) ! y - f

            if ( got(f_measurementSD) ) call rowScale ( weight, jacobian )

            !{Form normal equations:
            ! $\mathbf{J^T W^T W J \delta \hat x = J^T W^T W f}$:
            call formNormalEquations ( jacobian, normalEquations, &
              & rhs_in=f_rowScaled, rhs_out=atb, update=update, useMask=.true. )
            update = .true.
              if ( index(switches,'jac') /= 0 ) &
                call dump_Linf ( jacobian, 'L_infty norms of Jacobian blocks:' )
              if ( index(switches,'spa') /= 0 ) &
                & call dump_struct ( jacobian, &
                  & 'Sparseness structure of Jacobian blocks:' )
            call clearMatrix ( jacobian )  ! free the space
          end do ! mafs

          aj%fnorm = aj%fnorm + ( f_rowScaled .dot. f_rowScaled )

          !{Column Scale $\mathbf{J}$ (row and column scale $\mathbf{J^T
          ! J}$).  Let $\mathbf{\Sigma}$ be a column-scaling matrix for the
          ! Jacobian. We cater for several choices.  After row and column
          ! scaling, the problem is $\mathbf{J^T W^T W J \Sigma \Sigma^{-1}
          ! \delta \hat x = -J^T W^T f}$.
          select case ( columnScaling )
          case ( l_apriori )
            call copyVector ( columnScaleVector, apriori )
          case ( l_covariance )
            !??? Can't get here until allowed by init_tables
          case ( l_norm )
            call getDiagonal ( normalEquations%m, columnScaleVector )
          end select
          if ( columnScaling /= l_none ) then ! Compute $\Sigma$
            forall (j = 1: columnScaleVector%template%noQuantities)
              where ( columnScaleVector%quantities(j)%values <= 0.0 )
                columnScaleVector%quantities(j)%values = 1.0
              elsewhere
                columnScaleVector%quantities(j)%values = 1.0 / &
                  & sqrt( columnScaleVector%quantities(j)%values )
              end where
            end forall
              if ( index(switches,'col') /= 0 ) &
                & call dump ( columnScaleVector, name='Column scale vector' )
            !{Scale: $\mathbf{\Sigma^T J^T J \Sigma = -\Sigma^T J^T f}$
            call columnScale ( normalEquations%m, columnScaleVector )
            call rowScale ( columnScaleVector, normalEquations%m )
            call multiply ( atb, columnScaleVector )
  
          end if
          !{Compute the (negative of the) gradient $= -\mathbf{J^T f}$.
          ! This is the right-hand side of the normal equations
          ! $\mathbf{J^T J \delta \hat x = -J^T f}$ where $\mathbf{\delta
          ! \hat x}$ is the "Candidate $\mathbf{\delta x}$" that may or
          ! may not get used.
          call copyVector ( gradient, atb ) ! gradient := atb
            if ( index(switches,'gvec') /= 0 ) &
              & call dump ( gradient, name='gradient' )

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
              call getDiagonal ( factored%m, dxUnscaled )
              call dump ( dxUnscaled, &
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
          call solveCholesky ( factored, candidateDX, atb, &
            & transpose=.true. )
          call add_to_retrieval_timing( 'cholesky_solver', t1 )
          call cpu_time ( t1 )

          !{AJ\%FNMIN = L2 norm of residual, $||\mathbf{J \delta x + f}||$
          ! where $\mathbf{\delta x}$ is the "Candidate DX" that may not
          ! get used. This can be gotten without saving $\bf J$ as
          ! $\mathbf{f^T f - y^T y}$ where $\mathbf{y = U^{-T} \delta x}$
          ! ($\mathbf{U}$ is the Cholesky factor of $\mathbf{J^T J}$).
          ! The variable {\tt candidateDX} is a temp here = $\bf y$.
          aj%fnmin = aj%fnorm - (candidateDX .dot. candidateDX)
          if ( aj%fnmin < 0.0 ) then
            call output ( 'How can aj%fnmin be negative?  aj%fnmin = ' )
            call output ( aj%fnmin, advance='yes' )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Norm of residual is imaginary!' )
          end if
          aj%fnmin = sqrt(aj%fnmin)
          aj%fnorm = sqrt(aj%fnorm)
          aj%gradn = sqrt(gradient .dot. gradient) ! L2Norm(gradient)
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
          call solveCholesky ( factored, candidateDX, atb, &
            & transpose=.true. )
          call add_to_retrieval_timing( 'cholesky_solver', t1 )
          call cpu_time ( t1 )
          aj%fnmin = sqrt(aj%fnorm**2 - (candidateDX .dot. candidateDX) )
          call solveCholesky ( factored, candidateDX )
          aj%dxn = sqrt(candidateDX .dot. candidateDX) ! L2Norm(dx)
          aj%gdx = gradient .dot. candidateDX
          if ( .not. aj%starting ) aj%dxdxl = dx .dot. candidateDX
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( candidateDX, name='CandidateDX' )
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
          if ( .not. aj%starting ) aj%dxdxl = dx .dot. candidateDX
          !{Account for column scaling.  We solved for $\mathbf{\Sigma^{-1}
          ! \delta x}$ above, so multiply by $\Sigma$ (which is our
          ! variable {\tt columnScaleVector}):
          call copyVector ( dxUnScaled, dx )
          if ( columnScaling /= l_none ) then
            ! dxUnScaled = dxUnScaled # columnScaleVector:
            call multiply ( dxUnScaled, columnScaleVector )
          end if
          call addToVector ( x, dxUnScaled ) ! x = x + dxUnScaled
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( dxUnScaled, name='dX Unscaled' )
            if ( index(switches,'xvec') /= 0 ) &
              & call dump ( x, name='New X' )
          aj%axmax = 0.0
          aj%big = .false.
          do j = 1, size(x%quantities)
            aj%axmax = max(aj%axmax, maxval(abs(x%quantities(j)%values)))
            aj%big = aj%big .or. any( abs(dx%quantities(j)%values) > &
              & 10.0 * epsilon(aj%axmax) * abs(x%quantities(j)%values) )
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
          call copyVector ( x, bestX ) ! x = bestX
          ! dx = aj%gfac * "Best Gradient":
          call scaleVector ( bestGradient, aj%gfac, dx )
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( dx, name='Gradient move from best X' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%gfac = ' )
              call output ( aj%gfac, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_best ) ! ................................  BEST  .....
        ! Set "Best X" = X, "Best Gradient" = Gradient
          call copyVector ( bestX, x ) ! bestX = x
          call copyVector ( bestGradient, gradient ) ! bestGradient = gradient
        case ( nf_aitken ) ! ............................  AITKEN  .....
        ! Set DX = DX - "Candidate DX",
        !     AJ%DXDX = dot_product( DX, DX )
        ! IF ( AJ%DXDX /= 0.0 ) &
        !   Set AJ%DXDXL = dot_product( DX, "Candidate DX" )
          call subtractFromVector ( dx, candidateDX ) ! dx = dx - candidateDX
          aj%dxdx = dx .dot. dx
          if ( aj%dxdx > 0.0 ) aj%dxdxl = dx .dot. candidateDX
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( dx, name='dx after Aitken' )
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
          call copyVector ( dx, candidateDX ) ! dx = candidateDX
        case ( nf_dx_aitken ) ! ......................  DX_AITKEN  .....
          ! dx = aj%cait * candidateDX:
          call scaleVector ( candidateDX, aj%cait, dx )
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( dx, name='dx after dx Aitken' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%cait = ' )
              call output ( aj%cait, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_tolx, nf_tolx_best, nf_tolf, nf_too_small ) ! ........
          ! IF ( NWT_FLAG == NF_TOO_SMALL ) THEN
          !   Take special action if requested accuracy is critical
          ! END IF
          if ( nwt_flag == nf_tolx_best ) call copyVector ( x, bestX )
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
      end do ! Newton iteration
      if ( got(f_outputCovariance) .or. got(f_outputSD) ) then
        if ( got(f_diagnostics) ) then
          ! Compute rows of Jacobian actually used.  Don't count rows due
          ! to Levenberg-Marquardt stabilization.  Do count rows due to a
          ! priori or regularization.  Put numbers of rows and columns
          ! into diagnostic vector.
          jac_col_qty => GetVectorQuantityByType ( diagnostics, &
            & quantityType=l_jacobian_cols, noError=.true. )
          jac_row_qty => GetVectorQuantityByType ( diagnostics, &
            & quantityType=l_jacobian_rows, noError=.true. )
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
          if ( associated(jac_col_qty) ) jac_col_qty%values(1,1) = jacobian_cols
          if ( associated(jac_row_qty) ) jac_row_qty%values(1,1) = jacobian_rows
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
            call multiply ( outputSD, columnScaleVector )
          end if
        end if
      end if
      call copyVector ( state, x )
        if ( index(switches,'svec') /= 0 ) &
          & call dump ( state, name='Final state' )
      ! Clean up the temporaries, so we don't have a memory leak
      call destroyVectorInfo ( atb )
      call destroyVectorInfo ( bestGradient )
      call destroyVectorInfo ( bestX )
      call destroyVectorInfo ( candidateDX )
      if ( columnScaling /= l_none ) &
        & call destroyVectorInfo ( columnScaleVector )
      if ( diagonal ) call destroyVectorInfo ( covarianceDiag )
      if ( got(f_apriori) ) call destroyVectorInfo ( covarianceXApriori )
      call destroyVectorInfo ( dx )
      call destroyVectorInfo ( dxUnscaled )
      call destroyVectorInfo ( f_rowScaled )
      if ( got(f_fuzz) ) call destroyVectorInfo ( fuzzState )
      call destroyVectorInfo ( gradient )
      if ( got(f_measurementSD) ) call destroyVectorInfo ( weight )
      call destroyVectorInfo ( x )
      call destroyVectorInfo ( aprioriMinusX )
      call destroyMatrix ( normalEquations%m )
      call destroyMatrix ( factored%m )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      call add_to_retrieval_timing( 'newton_solver', t1 )
      call cpu_time ( t1 )
    end subroutine NewtonianSolver

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
        call output ( 'Elements per mask = ' )
        call output ( size(qty%values,1), advance='yes' )
        call dump ( qty%mask, format='(1x,z8.8)' )
      end if
    end subroutine SetupSubset
  end subroutine Retrieve

end module RetrievalModule

! $Log$
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
