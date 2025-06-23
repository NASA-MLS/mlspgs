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
module RetrievalModule
!=============================================================================

!{This module inverts the radiative transfer equation, to solve for
! atmospheric parameters, given radiance measurements.
!
! This module and ones it calls consume most of the cycles.

  implicit none
  private
  public :: RETRIEVE

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Aitken acceleration code was originally defective, in that DX was not
  ! being multiplied by the Aitken parameter AJ%CAIT.  If Do_Aitken is set
  ! true, Aitken acceleration is actually done.
  logical, parameter :: Do_Aitken = .false.

contains

  ! ---------------------------------------------------  Retrieve  -----
  subroutine RETRIEVE ( ROOT, VECTORDATABASE, MATRIXDATABASE, HESSIANDATABASE, &
    & CONFIGDATABASE, CHUNK, FILEDATABASE )

  !{Process the ``Retrieve'' section of the L2 Configuration File.
  ! The ``Retrieve'' section can have ForwardModel, Matrix, Sids, Subset or
  ! Retrieve specifications.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, TrackAllocates
    use Bitstuff, only: CountBits
    use Chunks_m, only: MLSChunk_t
    use CloudRetrievalModule, only: CloudRetrieval
    use DumpCommand_m, only: &
      & BooleanFromAnyGoodvalues, &
      & BooleanFromCatchWarning, BooleanFromComparingQtys, BooleanFromFormula, &
      & DumpCommand, InitializeRepeat, NextRepeat, &
      & MLSCase, MLSEndSelect, MLSSelect, MLSSelecting, Repeat=>Skip, Skip
    use HessianModule_1, only: Hessian_t, CreateEmptyHessian, DestroyHessian
    use IEEE_Arithmetic, only: iEEE_Is_NaN
    use Expr_m, only: Expr
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Init_TableS_Module, only: F_Apriori, F_AprioriFraction, F_AprioriScale, &
      & F_Average, F_ColumnScale, F_Comment, F_Covariance, F_CovSansReg, &
      & F_Diagnostics, F_Diagonal, F_DumpQuantities, F_ExtendedAverage, &
      & F_ForwardModel, F_Fuzz, F_FwdModelExtra, F_FwdModelOut, &
      & F_Hessian, F_HighBound, F_HregOrders, F_HregQuants, F_HregWeights, &
      & F_HregWeightVec, F_Jacobian, F_Lambda, F_LambdaMin, F_Level, &
      & F_LowBound, &
      & F_MaxJ, F_Measurements, F_MeasurementSD, F_Method, F_MuMin, &
      & F_NegateSD, &
      & F_OutputCovariance, F_OutputSD, &
      & F_PhaseName, F_PrecisionFactor, &
      & F_RegAfter, F_RegApriori, F_Serial, F_SparseQuantities, &
      & F_State, F_StateMax, F_StateMin, F_Switches, F_Toggles, &
      & F_ToleranceA, F_ToleranceF, F_ToleranceR, &
      & F_VregOrders, F_VregQuants, F_VregWeights, F_VregWeightVec, &
      & Field_First, Field_Last, &
      & L_Apriori, L_Covariance, &
      & L_DNWT_Abandoned,  L_DNWT_Ajn,  L_DNWT_Axmax, &
      & L_DNWT_Cait, L_DNWT_ChisqMinNorm, L_DNWT_ChisqNorm, L_DNWT_count, &
      & L_DNWT_Diag,  L_DNWT_dxdx, L_DNWT_dxdxl, L_DNWT_dxn,  L_DNWT_dxnl, &
      & L_DNWT_Flag, L_DNWT_Fnmin, L_DNWT_Fnorm,  L_DNWT_gdx,  L_DNWT_gfac, &
      & L_DNWT_gradn,  L_DNWT_sq, L_DNWT_sq,  L_DNWT_sqt, &
      & L_HighCloud, L_Jacobian_Cols, L_Jacobian_Rows, L_LowCloud, &
      & L_Newtonian, L_None, L_Norm, L_NumGrad, L_NumJ, L_NumNewt, L_Simple, &
      & S_AnyGoodValues, S_CatchWarning, S_Case, S_Compare, &
      & S_Diff, S_Dump, S_DumpBlocks, S_EndSelect, S_FlagCloud, S_FlushPFA, &
      & S_LeakCheck, S_NeuralNet, &
      & S_Reevaluate, S_Repeat, S_RestrictRange, S_Retrieve, &
      & S_Select, S_SIDS, S_Skip, S_Snoop, S_Subset, S_Time, S_UpdateMask
    use L2parinfo, only: Parallel
    use MatrixModule_1, only: AddToMatrixDatabase, CopyMatrix, CreateEmptyMatrix, &
      & DestroyMatrix, GetFromMatrixDatabase, NullifyMatrix, Matrix_t, &
      & Matrix_Database_t, Matrix_SPD_t, MultiplyMatrixVectorNoT, &
      & ReflectMatrix, Sparsify, MultiplyMatrix_XTY
    use MatrixTools, only: DumpBlocks
    use MLSKinds, only: R8, RV
    use MLSCommon, only: MLSFile_t
    use MLSL2Options, only: L2CFNode, L2Options, SpecialDumpFile, &
      & StateFilledBySkippedRetrievals
    use MLSL2timings, only: Section_Times, Total_Times, Add_To_Retrieval_Timing
    use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, MLSMessageReset, &
      & MLSMessage
    use MoreTree, only: Get_Boolean, Get_Label_And_Spec, Get_Field_Id, &
      & Get_Spec_Id
    use MLSStringlists, only: SwitchDetail
    use MLSStrings, only: WriteIntsToChars
    use NeuralNet_M, only: NeuralNet
    use Next_Tree_Node_m, only: Init_Next_Tree_Node, Next_Tree_Node, &
      & Next_Tree_Node_State
    use Output_m, only: Blanks, Output, RevertOutput, SwitchOutput
    use PFAData_m, only: Flush_PFAData
    use Set_Toggles_m, only: Set_Toggles
    use SIDSModule, only: SIDS
    use SnoopMLSL2, only: Snoop
    use String_Table, only: Display_String, Get_String
    use SubsetModule, only: SetupSubset, SetupFlagCloud, RestrictRange, UpdateMask
    use Time_m, only: Time_Now
    use Toggles, only: GEN, Switches, Toggle, Levels
    use Trace_m, only: Trace_Begin, Trace_End
    use Track_m, only: ReportLeaks
    use Tree, only: Decorate, Decoration, Node_Id, NSons, Sub_Rosa, Subtree, &
      & Where
    use Tree_Types, only: N_Array
    use Vectorsmodule, only: ClearMask, ClearUnderMask, &
      & ClearVector, CloneVector, CopyVector, CopyVectorMask, CreateMask, &
      & DestroyVectorInfo, DumpVectorNorms, GetVectorQuantityByType, M_Linalg, &
      & Vector_t, VectorValue_t

    ! Dummy arguments:
    integer, intent(in)                               :: Root ! Of the relevant subtree of the AST;
                                                              ! It indexes an n_cf vertex
    type(vector_T), dimension(:), target              :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer    :: MatrixDatabase
    type(Hessian_T), dimension(:), pointer            :: HessianDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: ConfigDatabase
    type(MLSChunk_T), intent(inout)                   :: Chunk
    type (MLSFile_T), dimension(:), pointer           :: FileDatabase

    ! Default values:
    real(r8), parameter :: DefaultInitLambda = 0.0_r8
    real(r8), parameter :: DefaultLambdaMin = 0.0_r8
    integer, parameter :: DefaultMaxJ = 5
    integer, parameter :: DefaultMethod = l_newtonian
    real(r8), parameter :: DefaultMuMin = 0.1_rv
    double precision, parameter :: DefaultToleranceA = 1.0d-6 ! for NWT
    double precision, parameter :: DefaultToleranceF = 1.0d-6 ! for NWT
    double precision, parameter :: DefaultToleranceR = 1.0d-6 ! for NWT

    ! Local variables:
    type(vector_T), pointer :: Apriori  ! A priori estimate of state
    type(vector_T), pointer :: AprioriFraction ! How much of result is apriori, in [0..1]
    real(r8) :: AprioriScale            ! Weight for apriori, default 1.0
    real(r8) :: ChiSqMinNorm            ! AJ%FNMIN / no degrees of freedom
    real(r8) :: ChiSqNorm               ! AJ%FNORM / no degrees of freedom
    integer :: ColumnScaling            ! one of l_apriori, l_covariance,
                                        ! l_none or l_norm
    integer, pointer, dimension(:) :: ConfigIndices    ! In ConfigDatabase
    type(matrix_SPD_T), pointer :: Covariance     ! covariance**(-1) of Apriori
    logical :: CovSansReg               ! Compute covariance without regularization
    type(vector_T), pointer :: Diagnostics   ! Diagnostic stuff about the
                                        ! retrieval.
    logical :: Diagonal                 ! "Iterate with the diagonal of the
                                        ! a priori covariance matrix until
                                        ! convergence, then put in the whole
                                        ! thing and iterate until it converges
                                        ! again (hopefully only once).
    integer :: DumpQuantitiesNode       ! in the l2cf tree
    integer :: Error
    logical :: ExtendedAverage          ! Averaging kernels based on full
                                        ! Jacobian (last term) not masked
    integer :: Field                    ! Field index -- f_something
    real(r8) :: Fuzz                    ! For testing only.  Amount of "fuzz"
                                        ! to add to state vector before
                                        ! starting a retrieval.
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    logical :: Got(field_first:field_last) ! "Got this field already"
    integer :: GSon                     ! a son of Son.
    type(vector_T), pointer :: HighBound ! For state during retrieval
    integer :: HRegOrders               ! Regularization orders
    integer :: HRegQuants               ! Regularization quantities
    integer :: HRegWeights              ! Weight of regularization conditions
    type(vector_T), pointer :: HRegWeightVec  ! Weight vector for regularization
    integer :: I_Key, J                 ! Subscripts and loop inductors
    real(r8) :: InitLambda              ! Initial Levenberg-Marquardt parameter
    integer :: IxAverage                ! Index in tree of averagingKernel
    integer :: IxCovariance             ! Index in tree of outputCovariance
    integer :: IxJacobian               ! Index in tree of Jacobian matrix
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    type(hessian_T), pointer :: Hessian ! The Hessian
    integer :: Jacobian_Cols            ! Number of columns of the Jacobian.
    integer :: Jacobian_Rows            ! (Number of rows of Jacobian) -
                                        ! (masked-off rows of Jacobian)
    integer :: K                        ! Subscript, loop inductor, local temp
    integer :: Key                      ! Index of an n_spec_args.  Either
                                        ! a son or grandson of root.
    integer :: Label                    ! Of a spec (not actually used)
    real(r8) :: LambdaMin               ! Minimum Levenberg-Marquardt parameter
    type(vector_T), pointer :: LowBound ! For state during retrieval
    integer :: MaxJacobians             ! Maximum number of Jacobian
                                        ! evaluations of Newton method
    type(vector_T), pointer :: Measurements    ! The measurements vector
    type(vector_T), pointer :: MeasurementSD   ! The measurements vector's Std. Dev.
    character(len=32) :: mesg
    character(len=32) :: mesgChunkNo
    integer :: Method                   ! Method to use for inversion, currently
                                        ! only l_Newtonian.
    type(matrix_T), target :: MyAverage ! for OutputAverage to point to
    type(matrix_SPD_T), target :: MyCovariance ! for OutputCovariance to point at
    type(hessian_T), target :: MyHessian       ! for Hessian to point at
    type(matrix_T), target :: MyJacobian       ! for Jacobian to point at
    real(rv) :: MuMin                   ! Smallest shrinking of dx before change direction
    logical :: NegateSD                 ! Flip output error negative for poor information
    type(matrix_T), pointer :: OutputAverage   ! Averaging Kernel
    type(matrix_SPD_T), pointer :: OutputCovariance    ! Covariance of the sol'n
    type(vector_T), pointer :: OutputSD ! Vector containing SD of result
    logical :: ParallelMode             ! Run forward models in parallel
    character(len=127) :: PhaseName     ! To pass to snoopers
    logical :: PotemkinHessian 
    real(rv) :: PrecisionFactor         ! Default 0.5, precisions 'worse than this' set negative
    logical :: REPEATLOOP               ! Do we repeat this section?
    integer :: SaveLevels(size(levels))
    logical :: SaveToggle(size(toggle))
    integer :: Son                      ! Of Root or Key
    character(len=127) :: SnoopComment  ! From comment= field of S_Snoop spec.
    integer :: SnoopKey                 ! Tree point of S_Snoop spec.
    integer :: SnoopLevel               ! From level field of S_Snoop spec.
    integer, dimension(:), pointer :: SparseQuantities ! Which jacobian blocks to sparsify
    integer :: Spec                     ! s_subset or s_retrieve ...
    type(vector_T), pointer :: State    ! The state vector
    type(vector_T), pointer :: StateMax ! Maximum state vector over all iterations
    type(vector_T), pointer :: StateMin ! Minimum state vector over all iterations
    integer :: Status
    integer :: SwitchLen                ! LEN_TRIM(Switches) on entry
    integer :: SwitchLenCur             ! LEN_TRIM(Switches) after command processing
    real :: T0, T1, T2, T3              ! for timing
    type(matrix_T) :: Tikhonov          ! Matrix for Tikhonov regularization
    logical :: TikhonovApriori          ! Regularization is to apriori, not
                                        ! zero -- default false
    logical :: TikhonovBefore           ! "Do Tikhonov before column scaling"
    logical :: TikhonovNeeded           ! got(f_hRegOrders) .or. got(f_vRegOrders)
    logical :: Timing
    double precision :: ToleranceA      ! convergence tolerance for NWT,
                                        ! norm of move
    double precision :: ToleranceF      ! convergence tolerance for NWT,
                                        ! norm of F
    double precision :: ToleranceR      ! convergence tolerance for NWT,
                                        ! (norm of move) / (norm of X)
    integer :: Units(2)                 ! Units of value returned by EXPR
    logical :: Update                   ! "We are updating normal equations"
    double precision :: Value(2)        ! Value returned by EXPR
    integer :: VRegOrders               ! Regularization orders
    integer :: VRegQuants               ! Regularization quantities
    integer :: VRegWeights              ! Weight of regularization conditions
    type(vector_T), pointer :: VRegWeightVec  ! Weight vector for regularization
    type(next_tree_node_state) :: Walk ! of tree traverser
    character(len=63) :: WhereLeakCheck ! From LeakCheck command, else default

    ! Indexes in the private vectors database
    integer, parameter :: FirstVec = 1
    integer, parameter :: AprioriMinusX = firstVec ! Apriori - X
    integer, parameter :: ATb = aprioriMinusX + 1  ! A^T b -- the RHS of the normal eqns    integer, parameter ::
    integer, parameter :: BestGradient = aTb + 1   ! for NWT
    integer, parameter :: BestX = bestGradient + 1 ! for NWT
    integer, parameter :: CandidateDX = bestX + 1  ! for NWT
    integer, parameter :: ColumnScaleVector = candidateDX + 1 ! For column scaling by column norms
                          ! ColumnScaleVector is really a diagonal matrix
    integer, parameter :: CovarianceDiag = columnScaleVector + 1
    integer, parameter :: CovarianceXApriori = covarianceDiag + 1
    integer, parameter :: DiagFlagA = CovarianceXApriori + 1  ! for NWT
    integer, parameter :: DiagFlagB = DiagFlagA + 1  ! for NWT
    integer, parameter :: DX = DiagFlagB + 1       ! for negateSD case
    integer, parameter :: DXUnScaled = dX + 1      ! for NWT
    integer, parameter :: F = dXUnscaled + 1       ! Residual -- Model - Measurements
    integer, parameter :: F_rowScaled = f + 1      ! Either a copy of f, or f row-scaled
    integer, parameter :: Gradient = f_rowScaled + 1   ! for NWT
    integer, parameter :: Reg_X_x = gradient + 1   ! Regularization * X_n
    integer, parameter :: Reg_RHS = reg_X_x + 1    ! RHS for Tikhonov if regApriori
    integer, parameter :: Weight = reg_RHS + 1     ! Scaling vector for rows, 1/measurementSD
    integer, parameter :: X = weight + 1           ! for NWT
    integer, parameter :: LastVec = X

    type(vector_T), dimension(firstVec:lastVec) :: V   ! Database for snoop

    ! Indices for trace strings, need to have negative initial values;
    ! See Push_Stack_b in Call_Stack
    integer :: Me = -1                  ! "Retrieve"
    integer :: Me_FlagCloud = -1        ! "Retrieve.flagCloud"
    integer :: Me_NeuralNet = -1        ! "Retrieve.NeuralNet
    integer :: Me_RestrictRange = -1    ! "Retrieve.RestrictRange
    integer :: Me_Retrieve = -1         ! "Retrieve.retrieve"
    integer :: Me_Subset = -1           ! "Retrieve.subset"
    integer :: Me_UpdateMask = -1       ! "Retrieve.UpdateMask"

    ! Error message codes
    integer, parameter :: ArrayOfArrays = 1
    integer, parameter :: BothOrNeither = ArrayOfArrays + 1
    integer, parameter :: DiagNoCov = BothOrNeither + 1
    integer, parameter :: IfAThenB = DiagNoCov + 1
    integer, parameter :: Inconsistent = IfAThenB + 1   ! Inconsistent fields
    integer, parameter :: NoExtraVec = Inconsistent + 1 ! no FwdModelExtra
    integer, parameter :: NotGeneral = NoExtraVec + 1   ! Not a general matrix
    integer, parameter :: NotSPD = notGeneral + 1       ! Not symmetric pos. definite

    dumpQuantitiesNode = 0
    error = 0
    nullify ( apriori, aprioriFraction, configIndices, covariance, fwdModelOut )
    nullify ( measurements, measurementSD, outputSD, sparseQuantities )
    nullify ( state, stateMax, stateMin )
    phaseName = ' '              ! Default in case there's no field
    snoopComment = ' '           ! Ditto
    snoopKey = 0
    snoopLevel = 1               ! Ditto
    switchLen = len_trim(switches)
    switchLenCur = switchLen + 1
    switches(switchLenCur:switchLenCur) = ','
    timing = section_times
    repeatLoop = .false. ! By default, we will not repeat
    call initializeRepeat
    do j = firstVec, lastVec ! Make the vectors in the database initially empty
      nullify ( v(j)%quantities, v(j)%template%quantities )
      v(j)%name = 0 ! so Snoop won't use it
    end do
    call time_now ( t1 )

    call trace_begin ( me, "Retrieve", root, cond=toggle(gen) )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
repeat_loop: do ! RepeatLoop
      ! Loop over the lines in the configuration file

      call init_next_tree_node ( walk )
      do 
        son = next_tree_node ( root, walk )
        if ( son == 0 ) exit
        call get_label_and_spec ( son, label, key )
        L2CFNODE = key

        ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
        ! parser users' guide" for pictures of the trees being analyzed.

        if ( MLSSelecting .and. &
          & .not. any( get_spec_id(key) == (/ s_endselect, s_select, s_case /) ) ) cycle
        if ( get_spec_id(key) /= s_catchWarning ) &
          & call MLSMessageReset( clearLastWarning=.true. )

        got = .false.
        spec = get_spec_id(key)
        PotemkinHessian = .false.
        select case ( spec )
        case ( s_anygoodvalues )
          call decorate ( key, &
            & BooleanFromAnyGoodValues ( key, vectorDatabase ) )
        case ( s_catchWarning )
          call decorate ( key,  BooleanFromCatchWarning ( key ) )
        case ( s_compare )
          call decorate ( key,  BooleanFromComparingQtys ( key, vectorDatabase ) )
        case ( s_diff, s_dump )
          if ( .not. L2Options%SkipRetrieval ) &
            & call dumpCommand ( key, forwardModelConfigs=configDatabase, &
            & vectors=vectorDatabase, FileDataBase=FileDataBase, &
            & MatrixDatabase=MatrixDatabase, Hessiandatabase=HessianDatabase )
        case ( s_dumpblocks )
          if ( .not. L2Options%SkipRetrieval ) &
            & call DumpBlocks ( key, matrixDatabase, hessianDatabase )
        case ( s_flagCloud )
          call trace_begin ( me_flagCloud, "Retrieve.flagCloud", root, &
            & cond=toggle(gen) .and. levels(gen) > 0 )
          call SetupflagCloud ( key, vectorDatabase )
          call trace_end ( cond=toggle(gen) .and. levels(gen) > 0 )
        case ( s_flushPFA )
          call flush_PFAData ( key, status )
          error = max(error,status)
        case ( s_leakCheck )
          if ( trackAllocates > 0 ) then
            whereLeakCheck = "In retrieval module..."
            if ( nsons(key) > 1 ) then
              call get_string ( sub_rosa(subtree(2,subtree(2,key))), whereLeakCheck, &
                & strip=.true. )
            end if
            call reportLeaks ( whereLeakCheck )
          end if
        case ( s_NeuralNet )
          ! if ( L2Options%SkipRetrieval .and. STATEFILLEDBYSKIPPEDRETRIEVALS == 0. ) cycle
          call trace_begin ( Me_NeuralNet, "Retrieve.NeuralNet", root, &
            & cond=toggle(gen) .and. levels(gen) > 0 )
          call NeuralNet ( key, vectorDatabase, chunk, FileDatabase )
          call trace_end ( cond=toggle(gen) .and. levels(gen) > 0 )
        case ( s_restrictRange )
          call trace_begin ( me_restrictRange, "Retrieve.RestrictRange", root, &
            & cond=toggle(gen) .and. levels(gen) > 0 )
          call RestrictRange ( key, vectorDatabase )
          call trace_end ( cond=toggle(gen) .and. levels(gen) > 0 )
        case ( s_Reevaluate )
          call decorate ( key,  BooleanFromFormula ( 0, key, vectorDatabase ) )
        case ( s_retrieve )
          if ( L2Options%SkipRetrieval .and. STATEFILLEDBYSKIPPEDRETRIEVALS == 0. ) cycle
          call trace_begin ( me_retrieve, "Retrieve.retrieve", root, &
            & cond=toggle(gen) )
          saveLevels = levels
          saveToggle = toggle
          aprioriScale = 1.0
          columnScaling = l_none
          covSansReg = .false.
          diagonal = .false.
          extendedAverage = .false.
          hRegQuants = 0
          hRegWeights = 0
          nullify ( lowBound, highBound, hRegWeightVec )
          initLambda = defaultInitLambda
          lambdaMin = defaultLambdaMin
          maxJacobians = defaultMaxJ
          method = defaultMethod
          muMin = defaultMuMin
          negateSD = .false.
          precisionFactor = 0.5_rv
          toleranceA = defaultToleranceA
          toleranceF = defaultToleranceF
          toleranceR = defaultToleranceR
          vRegQuants = 0
          vRegWeights = 0
          parallelMode = parallel%fwmParallel .and. parallel%master
          tikhonovApriori = .false.
          tikhonovBefore = .true.
          tikhonovNeeded = .false.
          nullify ( vRegWeightVec )
          do i_key = 2, nsons(key) ! fields of the "retrieve" specification
            son = subtree(i_key, key)
            L2CFNODE = son
            field = get_field_id(son)  ! tree_checker prevents duplicates
            got(field) = .true.
            select case ( field )
            case ( f_apriori )
              apriori => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_aprioriFraction )
              aprioriFraction => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_average )
              ixAverage = decoration(subtree(2,son)) ! averagingKernel matrix vertex
            case ( f_columnScale )
              columnScaling = decoration(subtree(2,son))
            case ( f_covariance )      ! of apriori
              call getFromMatrixDatabase ( &
                & matrixDatabase(decoration(decoration(subtree(2,son)))), &
                & covariance)
              if ( .not. associated(covariance) ) &
                & call announceError ( notSPD, field )
            case ( f_covSansReg )
              covSansReg = get_Boolean(son)
            case ( f_diagnostics )
              diagnostics => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_diagonal )
              diagonal = get_Boolean(son)
            case ( f_dumpQuantities )
              dumpQuantitiesNode = son
            case ( f_extendedAverage )
              extendedAverage = get_boolean(son)
            case ( f_forwardModel )
              call allocate_test ( configIndices, nsons(son)-1, "ConfigIndices", &
                & moduleName )
              do k = 2, nsons(son)
                gson = subtree(k,son)
                if ( node_id(gson) == n_array ) then
                  call announceError ( arrayOfArrays, field )
                  configIndices(k-1) = -999
                else
                  configIndices(k-1) = decoration(decoration(gson))
                end if
              end do
            case ( f_fwdModelExtra )
              fwdModelExtra => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_fwdModelOut )
              fwdModelOut => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_highBound )
              highBound => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_hRegOrders )
              hRegOrders = son
              tikhonovNeeded = .true.
            case ( f_hRegQuants )
              hRegQuants = son
            case ( f_hRegWeights )
              hRegWeights = son
            case ( f_hRegWeightVec )
              hRegWeightVec => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_jacobian )
              ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
            case ( f_lowBound )
              lowBound => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_measurements )
              measurements => vectorDatabase(decoration(decoration(subtree(2,son))))
              call cloneVector ( v(f), measurements, vectorNameText='_f' )
            case ( f_measurementSD )
              measurementSD => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_method )
              method = decoration(subtree(2,son))
            case ( f_negateSD )
              negateSD = get_boolean ( son )
            case ( f_regAfter )
              tikhonovBefore = .not. get_Boolean(son)
            case ( f_regApriori )
              tikhonovApriori = get_Boolean(son)
            case ( f_serial )
              parallelMode = parallelMode .and. .not. get_boolean ( son )
            case ( f_outputCovariance )
              ixCovariance = decoration(subtree(2,son)) ! outCov: matrix vertex
            case ( f_outputSD )
              outputSD => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_sparseQuantities )
              call Allocate_Test ( sparseQuantities, nsons(son)-1, &
                & 'sparseQuantities', ModuleName )
              do k = 2, nsons(son)
                gson = subtree(k,son)
                if ( node_id(gson) == n_array ) then
                  call announceError ( arrayOfArrays, field )
                  sparseQuantities(k-1) = -999
                else
                  sparseQuantities(k-1) = decoration(decoration(gson))
                end if
              end do
            case ( f_state )
              state => vectorDatabase(decoration(decoration(subtree(2,son))))
            case ( f_stateMax )
              stateMax => vectorDatabase(decoration(decoration(subtree(2,son))))
              do k = 1, size(stateMax%quantities)
                stateMax%quantities(k)%values = -huge(0.0_rv)
              end do
            case ( f_stateMin )
              stateMin => vectorDatabase(decoration(decoration(subtree(2,son))))
              do k = 1, size(stateMin%quantities)
                stateMin%quantities(k)%values = huge(0.0_rv)
              end do
            case ( f_switches )
              call get_string ( sub_rosa(subtree(2,son)), switches(switchLenCur+1:), strip=.true. )
              switchLenCur = len_trim(switches) + 1
              switches(switchLenCur:switchLenCur) = ','
            case ( f_toggles )
              call get_string ( sub_rosa(subtree(2,son)), switches(switchLenCur+1:), strip=.true. )
              call set_toggles ( switches(switchLenCur+1:) )
              switches(switchLenCur+1:) = ''
            case ( f_aprioriScale, f_fuzz, f_lambda, f_lambdamin, f_maxJ, &
              &    f_muMin, f_precisionFactor, f_toleranceA, f_toleranceF, &
              &    f_toleranceR )
              call expr ( subtree(2,son), units, value )
              select case ( field )
              case ( f_aprioriScale )
                aprioriScale = value(1)
              case ( f_fuzz )
                fuzz = value(1)
              case ( f_lambda )
                initLambda = value(1)
              case ( f_lambdamin )
                lambdaMin = value(1)
              case ( f_maxJ )
                maxJacobians = nint(value(1))
              case ( f_muMin )
                muMin = value(1)
              case ( f_precisionFactor )
                precisionFactor = value(1)
              case ( f_toleranceA )
                toleranceA = value(1)
              case ( f_toleranceF )
                toleranceF = value(1)
              case ( f_toleranceR )
                toleranceR = value(1)
              end select
            case ( f_vRegOrders )
              vRegOrders = son
              tikhonovNeeded = .true.
            case ( f_vRegQuants )
              vRegQuants = son
            case ( f_vRegWeights )
              vRegWeights = son
            case ( f_vRegWeightVec )
              vRegWeightVec => vectorDatabase(decoration(decoration(subtree(2,son))))
            case default
              ! Shouldn't get here if the type checker worked
            end select
          end do ! ! fields of the "retrieve" specification

          if ( L2Options%SkipRetrieval ) then
            if ( got(f_state) ) then
              call ClearVector( state, real(statefilledbyskippedretrievals, rv) )
              call output ( 'Clearing retrieval state vector name: ' )
              call display_string ( state%name )
              call output ( ', template name: ' )
              call display_string ( state%template%name, advance='yes' )
            end if
            if ( got(f_outputSD) ) then
              call ClearVector( outputSD, real(statefilledbyskippedretrievals, rv) )
              call output ( 'Clearing retrieval precision vector name: ' )
              call display_string ( state%name )
              call output ( ', template name: ' )
              call display_string ( state%template%name, advance='yes' )
            end if
            call output ( ' (skipping retrieval) ', advance='yes' )
            levels = saveLevels
            toggle = saveToggle
            call trace_end ( "Retrieve.retrieve", cond=toggle(gen) )
            cycle
          end if

          if ( got(f_apriori) .neqv. got(f_covariance) ) &
            & call announceError ( bothOrNeither, f_apriori, f_covariance )
          if ( diagonal .and. .not. got(f_covariance) ) &
            & call announceError ( diagNoCov )
          if ( got(f_regApriori) .and. .not. got(f_apriori) ) &
            & call announceError ( ifAThenB, f_regApriori, f_apriori )
          if ( got(f_hRegOrders) .neqv. &
            & (got(f_hRegWeights) .or. got(f_hRegWeightVec)) ) then
            call announceError ( bothOrNeither, f_hRegOrders, f_hRegWeights )
            call announceError ( bothOrNeither, f_hRegOrders, f_hRegWeightVec )
          end if
          if ( got(f_hRegQuants) .and. .not. got(f_hRegOrders) ) &
            & call announceError ( ifAThenB, f_hRegQuants, f_hRegOrders )
          if ( got(f_vRegOrders) .neqv. &
            & (got(f_vRegWeights) .or. got(f_vRegWeightVec)) ) then
            call announceError ( bothOrNeither, f_vRegOrders, f_vRegWeights )
            call announceError ( bothOrNeither, f_vRegOrders, f_vRegWeightVec )
          end if
          if ( got(f_vRegQuants) .and. .not. got(f_vRegOrders) ) &
            & call announceError ( ifAThenB, f_vRegQuants, f_vRegOrders )
          if ( negateSD .and. .not. tikhonovBefore ) &
            & call announceError ( inconsistent, f_negateSD, f_regAfter )
          if ( negateSD .and. .not. got(f_outputSD) ) &
            & call AnnounceError ( ifAThenB, f_negateSD, f_outputSD )
          if ( ( method == l_newtonian .or. method == l_simple ) .and. &
             & parallelMode .and. .not. got(f_fwdModelExtra) ) &
            & call announceError ( noExtraVec )

          if ( error == 0 ) then

            ! Verify the consistency of various matrices and vectors
            if ( got(f_apriori) ) then
              if ( apriori%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_apriori, f_state )
            end if
            if ( got(f_aprioriFraction) ) then
              if ( aprioriFraction%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_aprioriFraction, f_state )
            end if
            if ( associated(covariance) ) then
              if ( covariance%m%row%vec%template%name /= state%template%name .or. &
                &  covariance%m%col%vec%template%name /= state%template%name ) &
                &  call announceError ( inconsistent, f_covariance, f_state )
            end if
            if ( got(f_highBound) ) then
              if ( highBound%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_highBound, f_state )
            end if
            if ( got(f_lowBound) ) then
              if ( lowBound%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_lowBound, f_state )
            end if
            if ( got(f_measurementSD) ) then
              if ( measurementSD%template%name /= measurements%template%name ) &
                & call announceError ( inconsistent, f_measurementSD, f_measurements )
            end if
            if ( got(f_stateMax) ) then
              if ( stateMax%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_stateMax, f_state )
            end if
            if ( got(f_stateMin) ) then
              if ( stateMin%template%name /= state%template%name ) &
                & call announceError ( inconsistent, f_stateMin, f_state )
            end if
          end if

          if ( error == 0 ) then

            if ( switchDetail ( switches, 'rtv' ) > -1 ) call DumpRetrievalConfig

            ! Create the Hessian object
            if ( got(f_hessian) ) then
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Not ready to handle a Hessian field in Retrieve command yet' )
            else
              hessian => myHessian
              hessian = createEmptyHessian ( 0, measurements, state, Potemkin=.true. )
              PotemkinHessian = .true.
            end if

            ! Create the Jacobian matrix
            if ( got(f_jacobian) ) then
              k = decoration(ixJacobian)
              if ( k == 0 ) then
                call createEmptyMatrix ( myJacobian, &
                  & sub_rosa(subtree(1,ixJacobian)), measurements, state )
                k = addToMatrixDatabase( matrixDatabase, myJacobian )
                call decorate ( ixJacobian, k )
!     The following is probably necessary to make final subroutines for
!     RC_Info and Matrix_T work, else the finalizers deallocate targets
!     out from under pointers in the matrix database... but it causes a
!     crash.  I don't know why.
!                 call nullifyMatrix ( myJacobian ) ! so as not to clobber out
                  ! from under matrixDatabase(k) when myJacobian is finalized
              end if
              call getFromMatrixDatabase ( matrixDatabase(k), jacobian )
              if ( jacobian%row%vec%template%name /= measurements%template%name ) &
                & call announceError ( inconsistent, f_jacobian, f_measurements )
              if ( jacobian%col%vec%template%name /= state%template%name ) &
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
!     The following is probably necessary to make final subroutines for
!     RC_Info and Matrix_T work, else the finalizers deallocate targets
!     out from under pointers in the matrix database... but it causes a
!     crash.  I don't know why.
!                 call nullifyMatrix ( myCovariance%m ) ! so as not to clobber out
                  ! from under matrixDatabase(k) when myCovariance%m is finalized
              end if
              call getFromMatrixDatabase ( matrixDatabase(k), outputCovariance )
              if ( .not. associated(outputCovariance) ) then
                call announceError ( notSPD, f_outputCovariance )
              else
                if ( outputCovariance%m%row%vec%template%name /= state%template%name &
                  & .or. &
                  &  outputCovariance%m%col%vec%template%name /= state%template%name ) &
                  & call announceError ( inconsistent, f_outputCovariance, f_state )
              end if
            else
              outputCovariance => myCovariance
              call createEmptyMatrix ( myCovariance%m, 0, state, state )
            end if

            ! Create the averaging kernel matrix
            if ( got(f_average) ) then
              k = decoration(ixAverage)
              if ( k == 0 ) then
                call createEmptyMatrix ( myAverage, &
                  & sub_rosa(subtree(1,ixAverage)), state, state )
                k = addToMatrixDatabase( matrixDatabase, myAverage )
                call decorate ( ixAverage, k )
!     The following is probably necessary to make final subroutines for
!     RC_Info and Matrix_T work, else the finalizers deallocate targets
!     out from under pointers in the matrix database... but it causes a
!     crash.  I don't know why.
!                 call nullifyMatrix ( myAverage ) ! so as not to clobber out
                  ! from under matrixDatabase(k) when myAverage is finalized
              end if
              call getFromMatrixDatabase ( matrixDatabase(k), outputAverage )
              call display_string ( outputAverage%name, advance='yes' )
              if ( .not. associated(outputAverage) ) then
                call announceError ( notGeneral, f_average )
              else
                if ( outputAverage%row%vec%template%name /= state%template%name &
                  & .or. &
                  &  outputAverage%col%vec%template%name /= state%template%name ) &
                  & call announceError ( inconsistent, f_average, f_state )
              end if
            end if

            ! Create the matrix for adding the Tikhonov regularization to the
            ! normal equations
            if ( tikhonovNeeded ) then
              call createEmptyMatrix ( tikhonov, 0, state, state, text='_Tikhonov' )
            end if

            ! Do the retrieval
            jacobian_Cols = 0
            jacobian_Rows = 0
            if ( got(f_lowBound) ) call getInBounds ( state, lowBound, 'low' )
            if ( got(f_highBound) ) call getInBounds ( state, highBound, 'high' )

            select case (method)
            case( l_lowcloud, l_highcloud)
              ! use this for testing
              if(.not. got(f_maxJ)) maxJacobians = 5
              if(.not. got(f_lambda)) initlambda = 10.
              call CloudRetrieval(Method, ConfigDatabase,configIndices,fwdModelExtra,&
                 & measurements,MeasurementSD, state, OutputSD, Covariance, &
                 & jacobian, chunk,maxJacobians,initlambda)
              call add_to_retrieval_timing( 'low_cloud', t1 )
            case ( l_newtonian, l_simple )
              call newtonSolver
            case default
              call MLSMessage ( MLSMSG_Error, moduleName, &
              & "this retrieval method has not yet been implemented" )
            end select
          else
            call MLSMessage ( MLSMSG_Error, moduleName, &
              & "No retrieval done -- error in configuration" )
          end if

          !??? Make sure the jacobian and outputCovariance get destroyed
          !??? after ?what? happens?  Can we destroy the entire matrix
          !??? database at the end of each chunk?
          if ( .not. got(f_jacobian) ) call destroyMatrix ( jacobian )
          if ( .not. got(f_outputCovariance) ) &
            & call destroyMatrix ( outputCovariance%m )
          if ( got(f_fwdModelOut) ) then
            call copyVector ( fwdModelOut, v(f) )
            call copyVectorMask ( fwdModelOut, measurements )
          end if
          call deallocate_test ( configIndices, "ConfigIndices", moduleName )
          levels = saveLevels
          toggle = saveToggle
          call trace_end ( "Retrieve.retrieve", cond=toggle(gen) )
        case ( s_sids )
          if ( L2Options%SkipRetrieval .and. switchDetail( switches, 'fiw' ) < 0 ) cycle
          call time_now ( t1 )
          call sids ( key, VectorDatabase, MatrixDatabase, HessianDatabase, configDatabase, chunk)
        case ( s_select ) ! ============ Start of select .. case ==========
          ! We'll start seeking a matching case
          call MLSSelect (key)
        case ( s_case ) ! ================ seeking matching case ==========
          ! We'll continue seeking a match unless the case is TRUE
          call MLSCase (key)
        case ( s_endSelect ) ! =========== End of select .. case ==========
          ! We'done with seeking a match
          call MLSEndSelect (key)
        case ( s_Repeat ) ! ============================= Repeat ==========
          ! We'll Repeat the section as long as the Boolean cond'n is TRUE
          RepeatLoop = Repeat(key)
          if ( .not. RepeatLoop ) exit repeat_loop
          call nextRepeat
        case ( s_skip ) ! ================================= Skip ==========
          ! We'll skip the rest of the section if the Boolean cond'n is TRUE
          if ( Skip(key) ) exit repeat_loop
        case ( s_snoop )
          snoopKey = key
          do i_key = 2, nsons(key)
            son = subtree(i_key, key)
            field = get_field_id(son)  ! tree_checker prevents duplicates
            select case ( field )
            case ( f_comment )
              call get_string ( sub_rosa(subtree(2,son)), snoopComment, strip=.true. )
            case ( f_phaseName )
              call get_string ( sub_rosa(subtree(2,son)), phaseName, strip=.true. )
            case ( f_level )
              call expr ( subtree(2,son), units, value )
              snoopLevel = nint(value(1))
            end select
          end do
        case ( s_subset )
          call trace_begin ( me_subset, "Retrieve.subset", root, &
            & cond=toggle(gen) .and. levels(gen) > 0 )
          call SetupSubset ( key, vectorDatabase )
          call trace_end ( cond=toggle(gen) .and. levels(gen) > 0 )
        case ( s_time )
          if ( timing ) then
            call sayTime
          else
            call time_now ( t1 )
            timing = .true.
          end if
        case ( s_updateMask )
          call trace_begin ( me_updateMask, "Retrieve.UpdateMask", root, &
            & cond=toggle(gen) .and. levels(gen) > 0 )
          call UpdateMask ( key, vectorDatabase )
          call trace_end ( cond=toggle(gen) .and. levels(gen) > 0 )
        end select

        do j = firstVec, lastVec
          call destroyVectorInfo ( v(j) )
        end do

      end do ! sons
      if ( .not. repeatLoop ) exit ! Otherwise, we will repeat the section
    end do repeat_loop ! RepeatLoop
    
    ! Housekeeping
    ! if ( .not. got(f_jacobian) ) call DestroyMatrix( Jacobian )
    if ( PotemkinHessian ) call DestroyHessian( Hessian )

    if ( specialDumpFile /= ' ' ) call revertOutput

    switches(switchLen+1:) = '' ! Clobber switches from retrieve command
    call trace_end ( "Retrieve", cond=toggle(gen) )
    if ( timing ) call sayTime

  contains

    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex, AnotherFieldIndex, String )

      use Intrinsic, only: Field_Indices
      use Lexer_Core, only: Print_Source
      use String_Table, only: Display_String

      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex, AnotherFieldIndex ! f_...
      character(len=*), optional :: String

      error = max(error,1)
      call output ( '***** At ' )
      call print_source ( where(son) )
      call output ( ', RetrievalModule complained: ' )
      select case ( code )
      case ( arrayOfArrays )
        call output ( 'Field ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' cannot be an array of arrays', advance='yes' )
      case ( bothOrNeither )
        call output ( 'One of ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' or ' )
        call display_string ( field_indices(anotherFieldIndex) )
        call output ( ' appears, but the other does not.', advance='yes' )
      case ( diagNoCov )
        call output ( 'If the diagonal field is true, covariance shall appear', &
          & advance='yes' )
      case ( ifAThenB )
        call output ( 'If the ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' field appears then the ' )
        call display_string ( field_indices(anotherFieldIndex) )
        call output ( ' field shall also appear.', advance='yes' )
      case ( noExtraVec )
        call output ( 'FwdModelExtra is required for parallel forward models.', &
                    & advance='yes' )
      case ( inconsistent, notGeneral, notSPD )
        if ( present(string) ) call output ( string )
        call output ( 'the field ' )
        call display_string ( field_indices(fieldIndex) )
        select case ( code )
        case ( inconsistent )
          call output ( ' is not consistent with the ' )
          call display_string ( field_indices(anotherFieldIndex ) )
          call output ( ' field.', advance='yes' )
        case ( notGeneral )
          call output ( ' is not a general matrix.', advance='yes' )
        case ( notSPD )
          call output ( ' is not a symmetric positive-definite matrix.', &
            & advance='yes' )
        end select
      end select
    end subroutine AnnounceError

    ! --------------------------------------------  ApplyTikhonov  -----
    subroutine ApplyTikhonov ( After, NormalEquations, AJ, TikhonovRows, &
                             & D_Reg, D_Fnorm )
      !{ Apply Tikhonov regularization.
      !  Tikhonov regularization is of the form ${\bf R x}_{n+1} \simeq
      !  {\bf 0}$ or ${\bf R x}_{n+1} \simeq {\bf a}$, where {\bf a} is
      !  the apriori. So that all of the parts of the problem are solving
      !  for ${\bf\delta x}$, we subtract ${\bf R x}_n$ from both sides
      !  to get ${\bf R \delta x} \simeq -{\bf R x}_n$ or ${\bf R \delta
      !  x} \simeq {\bf R} ( {\bf a - x}_n)$.
      !
      !  If {\tt After} is true, apply column scaling.
      !
      !  If requesting a posteriori covariance, a posteriori standard
      !  deviation, or an averaging kernel, without Tikhonov regularization,
      !  becomes common, this routine should be used to back out its effect
      !  by subtracting from the normal equations.  This would require adding
      !  a flag to the FormNormalEquations routine to subtract instead of add.

      use DNWT_Module, only: NWT_T
      use MatrixModule_1, only: ClearMatrix, ColumnScale, Dump, Dump_Struct, &
        & FormNormalEquations => NormalEquations, GetDiagonal
      use Regularization, only: Regularize
      use VectorsModule, only: OPERATOR(.DOT.), ScaleVector

      logical, intent(in) :: After     ! Allow column scaling
      type(matrix_SPD_T), intent(inout) :: NormalEquations  ! Jacobian**T * Jacobian
      type(Nwt_T), intent(inout) :: AJ ! Stuff for communicating with DNWT
      integer, intent(inout) :: TikhonovRows ! Number of rows added
      logical, intent(in) :: D_Reg     ! Dump regularization
      integer, intent(in) :: D_Fnorm   ! Dump the revised norm of F

      integer :: J                     ! Loop index
      integer :: T
      character(*), parameter :: Which(2) = (/ 'Vertical  ', 'Horizontal' /)

      ! call add_to_retrieval_timing( 'newton_solver', t1 )
      call time_now ( t1 )

      !{ Tikhonov regularization is of the form ${\bf R x}_{n+1} \simeq
      !  {\bf 0}$ or ${\bf R x}_{n+1} \simeq {\bf a}$, where {\bf a} is
      !  the apriori. So that all of the parts of the problem are solving
      !  for ${\bf\delta x}$, we subtract ${\bf R x}_n$ from both sides
      !  to get ${\bf R \delta x} \simeq -{\bf R x}_n$ or ${\bf R \delta
      !  x} \simeq {\bf R} ( {\bf a - x}_n)$.

      do t = 1, 2 ! Vertical, then Horizontal regularization
        if ( t == 1 ) then
          if ( .not. got(f_vRegOrders) ) cycle
          call regularize ( tikhonov, vRegOrders, vRegQuants, vRegWeights, &
            & vRegWeightVec, tikhonovRows, horiz=.false. )
        else
          if ( .not. got(f_hRegOrders) ) cycle
          call regularize ( tikhonov, hRegOrders, hRegQuants, hRegWeights, &
            & hRegWeightVec, tikhonovRows, horiz=.true. )
        end if

        !{ If Tiknonov regularization is ``the second derivative of the
        !  state ought to be like the second derivative of the apriori,''
        !  the RHS is ${\bf W R} ( {\bf x}_a - {\bf x}_n )$, where ${\bf
        !  W}$ is the Tikhonov weight, ${\bf R}$ is the regularization
        !  operator, ${\bf x}_a$ is apriori, and ${\bf x}_n$ is the
        !  current state.  If Tiknonov regularization is ``the second
        !  derivative of the state ought to be zero, the RHS is $-{\bf W
        !  R x}_n$.
        if ( tikhonovApriori ) then
          call multiplyMatrixVectorNoT ( tikhonov, v(reg_RHS), v(reg_X_x) )
        else
          call multiplyMatrixVectorNoT ( tikhonov, v(x), v(reg_X_x) )
          call scaleVector ( v(reg_X_x), -1.0_r8 )   ! -R x_n
        end if

        if ( after .and. columnScaling /= l_none ) then ! Compute $\Sigma$
          ! Get the column scale vector.
          select case ( columnScaling )
          case ( l_apriori ) ! v(columnScaleVector) := apriori
            call copyVector ( v(columnScaleVector), apriori )
          case ( l_covariance )
            !??? Can't get here until allowed by init_tables
          case ( l_norm ) ! v(columnScaleVector) := diagonal of normal
                          !                         equations = column norms^2
            call getDiagonal ( normalEquations%m, v(columnScaleVector) )
          end select
          do j = 1, v(columnScaleVector)%template%noQuantities
            where ( v(columnScaleVector)%quantities(j)%values <= 0.0 )
              v(columnScaleVector)%quantities(j)%values = 0.0
            elsewhere
              v(columnScaleVector)%quantities(j)%values = &
                & sqrt( v(columnScaleVector)%quantities(j)%values )
            end where
          end do
          call columnScale ( tikhonov, v(columnScaleVector) )
        end if

          if ( d_reg ) then
            if ( t == 1 ) then
              call output ( 'Dumping Tikhonov for vertical regularization', &
                & advance='yes' )
            else
              call output ( 'Dumping Tikhonov for horizontal regularization', &
                & advance='yes' )
            end if
            call dump_struct ( tikhonov, 'Tikhonov' )
            call dump ( tikhonov, name='Tikhonov', details=2 )
          end if

          call add_to_retrieval_timing( 'tikh_reg', t1 )
        call formNormalEquations ( tikhonov, normalEquations, &
          & v(reg_X_x), v(aTb), update=update, useMask=.false. )
        update = .true.
        call clearMatrix ( tikhonov )           ! free the space
        ! aj%fnorm is still the square of the norm of f
        aj%fnorm = aj%fnorm + ( v(reg_X_x) .dot. v(reg_X_x) )
        if ( d_fnorm > -1 ) then
          call output ( trim(which(t)) // &
            & ' Regularization contribution to | F |^2 = ' )
          call output ( v(reg_X_x) .dot. v(reg_X_x), advance='yes' )
        end if
        ! call destroyVectorValue ( v(reg_X_x) )  ! free the space
        ! Don't destroy reg_X_x unless we move the 'clone' for it
        ! inside the loop.  Also, if we destroy it, we can't snoop it.
      end do ! t

    end subroutine ApplyTikhonov

    ! ------------------------------------------------  BoundMove  -----
    subroutine BoundMove ( Mu, Bound, X, Dx, Which, muMin )
      ! Compute a scalar multiple Mu such that X + Mu*DX does not go
      ! beyond Bound.  "Which" specifies either "low" or "high"

      !{ Let $D$ be $|\delta \mathbf{x}|$, where $\delta \mathbf{x}$ is
      ! given by Dx, let $b_i$ and $x_i$ be components of Bound and X. Let
      ! $d_i$ be such that $x_i + \frac{|d_i|}D \delta x_i$ is not beyond
      ! a bound.  That is, $d_i$ is such that $|d_i| \cos \theta_i = x_i -
      ! b_i$, where $\theta_i$ is the angle between $\delta \mathbf{x}$
      ! and the normal to the $i^{th}$ constraint surface.  Since the
      ! constraints are bounds, we have $\cos \theta_i = \frac{\delta
      ! x_i}D$.  $\mu$ is the smallest value of $\frac{|d_i|}D =
      ! \left|\frac{x_i-b_i} {\delta x_i}\right|$.  If we check the
      ! components one at a time, each one using the most recently
      ! computed $\mu$, we don't need a divide for each check, and we
      ! don't need a {\tt min} operation. However, if we are already at
      ! the bounds, and the next step wants to keep going beyond the
      ! bounds, then $\mu$ will be really small.  In this case $\mu <
      ! \mu_{\text{min}}$ we revert to a straight element-by-element
      ! modification of $\delta x_i$.

      real(rv), intent(inout) :: Mu
      type(vector_T), intent(inout) :: Bound, X, Dx
      character(len=*), intent(in) :: Which
      real(rv), intent(in) :: MuMin

      integer :: IQ, IVX, IVY           ! Subscripts used during MU computation
      real(rv) :: MuOrig                ! The original value of MU


!     This is a set of debugging output it is helpful to sprinkle in when
!     tracking down problems.
!       call Output ( 'In Bound move (' // which // ')', advance='yes' )
!       call output ( ' x = ' )
!       call output ( x%quantities(iq)%values(ivx,ivy) )
!       call output ( ' dx = ' )
!       call output ( dx%quantities(iq)%values(ivx,ivy) )
!       call output ( ' bound = ' )
!       call output ( bound%quantities(iq)%values(ivx,ivy) )
!       call output ( ' mu = ' )
!       call output ( mu, advance='yes' )

      muOrig = mu
      if ( which == 'low' ) then
        do iq = 1, size(x%quantities)
          if ( associated(x%quantities(iq)%mask) ) then
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      & mu = abs( ( bound%quantities(iq)%values(ivx,ivy) - &
                      &             x%quantities(iq)%values(ivx,ivy) ) &
                      &           / dx%quantities(iq)%values(ivx,ivy) )
                end if
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( x%quantities(iq)%values(ivx,ivy) + &
                  &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                  &  bound%quantities(iq)%values(ivx,ivy) ) &
                    & mu = abs( ( bound%quantities(iq)%values(ivx,ivy) - &
                    &             x%quantities(iq)%values(ivx,ivy) ) &
                    &           / dx%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          end if
        end do
      else if ( which == 'high' ) then
        do iq = 1, size(x%quantities)
          if ( associated(x%quantities(iq)%mask) ) then
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  mu * dx%quantities(iq)%values(ivx,ivy) > &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      & mu = abs( ( bound%quantities(iq)%values(ivx,ivy) - &
                      &             x%quantities(iq)%values(ivx,ivy) ) &
                      &           / dx%quantities(iq)%values(ivx,ivy) )
                end if
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( x%quantities(iq)%values(ivx,ivy) + &
                  &  mu * dx%quantities(iq)%values(ivx,ivy) > &
                  &  bound%quantities(iq)%values(ivx,ivy) ) &
                    & mu = abs( ( bound%quantities(iq)%values(ivx,ivy) - &
                    &             x%quantities(iq)%values(ivx,ivy) ) &
                    &           / dx%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          end if
        end do
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Come on! In BoundMove, it has to be a low bound or a high bound!' )
      end if

      ! Now if mu has gotten really small, we'll change it back to one,
      ! and do an element-by-element modification of dx
      if ( abs(mu) < muMin ) then
        mu = muOrig
        if ( which == 'low' ) then
          do iq = 1, size(x%quantities)
            if ( associated(x%quantities(iq)%mask) ) then
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                    if ( x%quantities(iq)%values(ivx,ivy) + &
                      &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                      &  bound%quantities(iq)%values(ivx,ivy) ) &
                        &  dx%quantities(iq)%values(ivx,ivy) = &
                        &    ( bound%quantities(iq)%values(ivx,ivy) - &
                        &    x%quantities(iq)%values(ivx,ivy) ) / mu
                  end if
                end do
              end do
            else
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      &  dx%quantities(iq)%values(ivx,ivy) = &
                      &    ( bound%quantities(iq)%values(ivx,ivy) - &
                      &    x%quantities(iq)%values(ivx,ivy) ) / mu
                end do
              end do
            end if
          end do
        else if ( which == 'high' ) then
          do iq = 1, size(x%quantities)
            if ( associated(x%quantities(iq)%mask) ) then
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                    if ( x%quantities(iq)%values(ivx,ivy) + &
                      &  mu * dx%quantities(iq)%values(ivx,ivy) > &
                      &  bound%quantities(iq)%values(ivx,ivy) ) &
                        &  dx%quantities(iq)%values(ivx,ivy) = &
                        &    ( bound%quantities(iq)%values(ivx,ivy) - &
                        &      x%quantities(iq)%values(ivx,ivy) ) / mu
                  end if
                end do
              end do
            else
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  mu * dx%quantities(iq)%values(ivx,ivy) > &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      &  dx%quantities(iq)%values(ivx,ivy) = &
                      &    ( bound%quantities(iq)%values(ivx,ivy) - &
                      &      x%quantities(iq)%values(ivx,ivy) ) / mu
                end do
              end do
            end if
          end do
        end if
      else if ( mu > 1.0 ) then ! in case dx(...) is very small
        mu = 1.0
      end if

      ! Let's keep this check in, just in case muMin has not been set.
      if ( mu < 0.0_rv ) then
        ! If it's negative but tiny then just nudge it.
        if ( mu > - sqrt ( epsilon ( 1.0_rv ) ) ) then
          mu = abs ( mu )
        else
          call output ( 'mu=' )
          call output ( mu, advance='yes' )
          call MLSMessage ( MLSMSG_Error, moduleName, &
            &  'How did mu get to be negative, might the initial guess be out of bounds?' )
        end if
      end if
    end subroutine BoundMove

    ! --------------------------------------  DumpStateQuantities  -----
    subroutine DumpStateQuantities ( DumpQuantitiesNode, Title )
      use Tree, only: Decoration, NSons, Subtree
      use VectorsModule, only: Dump, GetVectorQtyByTemplateIndex, &
        & VectorValue_t

      integer, intent(in) :: DumpQuantitiesNode ! in L2CF tree
      character(len=*), intent(in) :: Title

      integer :: Gson, I
      type(vectorValue_t), pointer :: Qty

      ! Sons of dumpQuantitiesNode are <dot vector quantity> trees.
      do i = 2, nsons(dumpQuantitiesNode)
        gson = subtree(i,dumpQuantitiesNode)
        qty => getVectorQtyByTemplateIndex ( &
          & vectorDatabase(decoration(decoration(subtree(1,gson)))), &
          & decoration(decoration(decoration(subtree(2,gson)))) )
        call dump ( qty, name=title )
      end do

    end subroutine DumpStateQuantities

    ! --------------------------------------  DumpRetrievalConfig  -----
    subroutine DumpRetrievalConfig
      use lexer_core, only: PRINT_SOURCE
      use VectorsModule, only: DumpNiceMaskSummary, M_Tikhonov, M_FullDerivatives
      ! Local variables
      integer :: Q                ! Loop counter
      ! Executable code
      call output ( '---------------------------- Begin retrieval configuration dump', advance='yes' )
      call output ( 'Dumping retrieval configuration for retrieve statement at ' )
      call print_source ( where ( son ) )
      call output ( '', advance='yes' )
      call output ( 'Retrieval state vector name: ' )
      call display_string ( state%name )
      call output ( ', template name: ' )
      call display_string ( state%template%name, advance='yes' )
      ! Now display the list of quantities retrieved
      do q = 1, state%template%noQuantities
        call output ( 'Retrieved quantity ' )
        call output ( q )
        call output ( ' ( ' )
        call display_string ( state%quantities(q)%template%name )
        call output ( ' ) ' )
        call output ( state%quantities(q)%template%noInstances )
        if ( state%quantities(q)%template%noInstances == 1 ) then
          call output ( ' instance, ' )
        else
          call output ( ' instances, ' )
        end if
        call output ( state%quantities(q)%template%noSurfs )
        if ( state%quantities(q)%template%noSurfs == 1 ) then
          call output ( ' surface, ' )
        else
          call output ( ' surfaces, ' )
        end if
        call output ( state%quantities(q)%template%noChans )
        if ( state%quantities(q)%template%noChans == 1 ) then
          call output ( ' channel.', advance='yes' )
        else
          call output ( ' channel.', advance='yes' )
        end if
        ! Call a routine to dump a nice summary of the retrieval range etc.
        call DumpNiceMaskSummary ( state%quantities(q), '  ', &
          & (/ m_linAlg, m_tikhonov, m_fullDerivatives /) )
      end do

      ! Now display the list of measurements used
      call output ( 'Measurements vector name: ' )
      call display_string ( measurements%name )
      call output ( ', template name: ' )
      call display_string ( measurements%template%name, advance='yes' )
      ! Now display the list of quantities retrieved
      do q = 1, measurements%template%noQuantities
        call output ( 'Measurement quantity ' )
        call output ( q )
        call output ( ' ( ' )
        call display_string ( measurements%quantities(q)%template%name )
        call output ( ' ) ' )
        call output ( measurements%quantities(q)%template%noInstances )
        if ( measurements%quantities(q)%template%noInstances == 1 ) then
          call output ( ' instance, ' )
        else
          call output ( ' instances, ' )
        end if
        call output ( measurements%quantities(q)%template%noSurfs )
        if ( measurements%quantities(q)%template%noSurfs == 1 ) then
          call output ( ' surface, ' )
        else
          call output ( ' surfaces, ' )
        end if
        call output ( measurements%quantities(q)%template%noChans )
        if ( measurements%quantities(q)%template%noChans == 1 ) then
          call output ( ' channel.', advance='yes' )
        else
          call output ( ' channel.', advance='yes' )
        end if
        ! Call a routine to dump a nice summary of the retrieval range etc.
        call DumpNiceMaskSummary ( measurements%quantities(q), '  ', (/ m_linalg /) )
      end do


      call output ( '---------------------------- End retrieval configuration dump', advance='yes' )
      stop
    end subroutine DumpRetrievalConfig

    ! ----------------------------------------------  FillDiagQty  -----
    subroutine FillDiagQty ( Diagnostics, QuantityIndex, Value )
      ! Put Value into the first element of the values field of the
      ! quantity of the Diagnostics vector given by QuantityIndex.
      ! We assume the Diagnostics vector is associated, so you better
      ! check before you call!
      type (Vector_T), intent(inout) :: DIAGNOSTICS
      integer, intent(in) :: QuantityIndex
      real(r8), intent(in) :: Value

      integer :: Latest
      integer :: Me = -1                ! String index for trace
      type(vectorValue_T), pointer :: Diag_Qty    ! A quantity in the
                                        ! Diagnostics vector
      call trace_begin ( me, 'Retrieve.FillDiagQty', cond=.false. )
      diag_qty => GetVectorQuantityByType ( diagnostics, &
        & quantityType=QuantityIndex, noError=.true. )
      if ( associated(diag_qty) ) then
        if ( diag_qty%template%noSurfs == 1 ) then
          ! Just one 'surface', or no mask information
          ! write this value out
          diag_qty%values(1,1) = value
        else
          if ( .not. associated ( diag_qty%mask ) ) then
            call CreateMask ( diag_qty )
            diag_qty%mask = char ( m_linalg )
          end if
          ! Multiple surfaces here, find the 'latest'
          latest = count ( iand(ichar(diag_qty%mask(:,1)),m_linalg) == 0 )
          if ( latest == diag_qty%template%noSurfs ) then
            diag_qty%values ( :, 1 ) = &
              & (/ diag_qty%values ( 2:latest, 1 ), value /)
          else
            diag_qty%values ( latest+1, 1 ) = value
            call ClearMask ( diag_qty%mask(:,1), (/ latest+1 /), m_linalg )
          end if
        end if
      end if
      call trace_end ( 'Retrieve.FillDiagQty', cond=.false. )
    end subroutine FillDiagQty

    ! ----------------------------------------------  FillDiagVec  -----
    subroutine FillDiagVec ( Diagnostics, AJ, NumGrad, NumJ, NumNewt, NWT_Flag, &
      & Jacobian_Rows, Jacobian_Cols )
      use DNWT_Module, only: NWT_T
      type(Vector_T), intent(inout) :: Diagnostics
      type(Nwt_T), intent(in) :: AJ
      integer, intent(in), optional :: NumGrad
      integer, intent(in), optional :: NumJ
      integer, intent(in), optional :: NumNewt
      integer, intent(in), optional :: NWT_Flag
      integer, intent(in), optional :: Jacobian_Rows
      integer, intent(in), optional :: Jacobian_Cols

      integer :: Me = -1               ! String index for trace

      call trace_begin ( me, 'Retrieve.FillDiagVec', cond=.false. )
      call fillDiagQty ( diagnostics,  l_dnwt_ajn, aj%ajn )
      call fillDiagQty ( diagnostics,  l_dnwt_axmax, aj%axmax )
      call fillDiagQty ( diagnostics,  l_dnwt_cait, aj%cait )
      call fillDiagQty ( diagnostics,  l_dnwt_chiSqMinNorm, chiSqMinNorm )
      call fillDiagQty ( diagnostics,  l_dnwt_chiSqNorm, chiSqNorm )
      call fillDiagQty ( diagnostics,  l_dnwt_diag, aj%diag )
      call fillDiagQty ( diagnostics,  l_dnwt_dxdx, aj%dxdx )
      call fillDiagQty ( diagnostics,  l_dnwt_dxdxl, aj%dxdxl )
      call fillDiagQty ( diagnostics,  l_dnwt_dxn, aj%dxn )
      call fillDiagQty ( diagnostics,  l_dnwt_dxnl, aj%dxnl )
      call fillDiagQty ( diagnostics,  l_dnwt_fnmin, aj%fnmin )
      call fillDiagQty ( diagnostics,  l_dnwt_fnorm, aj%fnorm )
      call fillDiagQty ( diagnostics,  l_dnwt_gdx, aj%gdx )
      call fillDiagQty ( diagnostics,  l_dnwt_gfac, aj%gfac )
      call fillDiagQty ( diagnostics,  l_dnwt_gradn, aj%gradn )
      call fillDiagQty ( diagnostics,  l_dnwt_sq, aj%sq )
      call fillDiagQty ( diagnostics,  l_dnwt_sqt, aj%sqt )
      if ( present ( numGrad ) ) &
        & call fillDiagQty ( diagnostics,  l_numGrad, real(numGrad, rv) )
      if ( present ( numJ ) ) &
        & call fillDiagQty ( diagnostics,  l_numJ, real(numJ, rv) )
      if ( present ( numNewt ) ) &
        & call fillDiagQty ( diagnostics,  l_numNewt, real(numNewt, rv) )
      if ( present ( nwt_flag ) ) &
        & call fillDiagQty ( diagnostics,  l_dnwt_flag, real(nwt_flag,rv) )
      if ( present ( jacobian_rows ) ) &
        & call fillDiagQty ( diagnostics,  l_jacobian_rows, real(jacobian_rows,rv) )
      if ( present ( jacobian_cols ) ) &
        & call fillDiagQty ( diagnostics,  l_jacobian_cols, real(jacobian_cols,rv) )
      call trace_end ( 'Retrieve.FillDiagVec', cond=.false. )
    end subroutine FillDiagVec

    ! ----------------------------------------------  GetInBounds  -----
    subroutine GetInBounds ( X, Bound, Which )
      ! Put X above/below Bound.  Which specifies low or high bound.
      type(vector_t), intent(inout) :: X
      type(vector_t), intent(in) :: Bound
      character(len=*), intent(in) :: Which

      integer :: IQ, IVX, IVY           ! Subscripts used during bounding

      if ( which == 'low' ) then
        do iq = 1, size(x%quantities)
          if ( associated(x%quantities(iq)%mask) ) then
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) /= 0 ) &
                  & x%quantities(iq)%values(ivx,ivy) = &
                    & max(x%quantities(iq)%values(ivx,ivy), &
                    &     bound%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                x%quantities(iq)%values(ivx,ivy) = &
                  & max(x%quantities(iq)%values(ivx,ivy), &
                  &     bound%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          end if
        end do
      else if ( which == 'high' ) then
        do iq = 1, size(x%quantities)
          if ( associated(x%quantities(iq)%mask) ) then
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) /= 0 ) &
                & x%quantities(iq)%values(ivx,ivy) = &
                    & min(x%quantities(iq)%values(ivx,ivy), &
                    &     bound%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                x%quantities(iq)%values(ivx,ivy) = &
                  & min(x%quantities(iq)%values(ivx,ivy), &
                  &     bound%quantities(iq)%values(ivx,ivy) )
              end do
            end do
          end if
       end do
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Come on! In GetInBounds, it has to be a low bound or a high bound!' )
      end if
    end subroutine GetInBounds

    ! -------------------------------------------------  My_NWTDB  -----
    subroutine My_NWTDB ( AJ, Width, Why )
      use DNWT_Module, only: NWT_T, NWTDB
      use DNWT_Clone, only: ALT_NWTDB
      type (NWT_T), intent(in) :: AJ
      integer, intent(in), optional :: Width
      character(len=*), intent(in), optional :: Why
      if ( method == l_newtonian ) then
        call nwtdb ( aj, width, why )
      else
        call alt_nwtdb ( aj, width, why )
      end if
    end subroutine My_NWTDB

    ! ---------------------------------------------  NewtonSolver  -----
    subroutine NewtonSolver

      use DNWT_MODULE, only: FlagName, NF_Aitken, NF_Best, NF_Biggest_Flag, &
      & NF_DX, NF_DX_Aitken, NF_EvalF, NF_EvalJ, NF_FandJ, NF_Gmove, NF_Lev, &
      & NF_NewX, NF_Smallest_Flag, NF_Solve, NF_Start, NF_TolX, NF_TolF, &
      & NF_TolX_Best, NF_Too_Small, NWT_t, NWT_Options, RK
      use DNWT_Module, only: NWT, NWTA
      use DNWT_Clone, only: ALT_NWT, ALT_NWTA ! SIMPLE
      use Dump_0, only: dump
      use ForwardModelWrappers, only: ForwardModel
      use Forwardmodelintermediate, only: forwardmodelstatus_t
      use L2FWMParallel, only: SetupFWMSlaves, TriggerSlaveRun, &
        & RequestSlavesOutput, ReceiveSlavesOutput
      use MatrixModule_0, only: Dump ! The one from MatrixModule_1 ought to work ???
      use MatrixModule_1, only: AddToMatrix, CholeskyFactor, ClearMatrix, &
        & ColumnScale, CopyMatrixValue, CreateEmptyMatrix, &
        & DestroyMatrix, Dump, Dump_Linf, Dump_Struct, &
        & FormNormalEquations => NormalEquations, &
        & GetDiagonal, InvertCholesky, Matrix_t, &
        & Matrix_Cholesky_t, Matrix_SPD_t, MaxL1, MinDiag, Multiply, &
        & MultiplyMatrix_XY,  MultiplyMatrix_XY_T,  &
        & RowScale, ScaleMatrix, SolveCholesky, UpdateDiagonal
      use ScanModelModule, only: DestroyForwardModelIntermediate
      use Symbol_Table, only: Enter_Terminal
      use Symbol_Types, only: T_Identifier
      use VectorsModule, only: AddToVector, DestroyVectorInfo, &
        & Dump, Multiply, OPERATOR(.dot.), OPERATOR(.mdot.), OPERATOR(-), &
        & ScaleVector, SubtractFromVector

      ! Local Variables
      ! logical :: Abandoned              ! Flag to indicate numerical problems
      real(r8) :: abandoned_value
      type(nwt_T) :: AJ                 ! "About the Jacobian", see NWT.
      type(matrix_SPD_T), target :: APlusRegCOV ! Covariance from a priori and Tikhnov alone
      type(matrix_SPD_T), target :: APlusRegNEQ ! Normal equations from a priori in Tikhnov alone.
                                        ! Only used when computing solution covariance/sd
                                        ! and needing to set values negative, or when
                                        ! computing the apriori fraction.
      logical :: AtBest                 ! Current state is the best one so far.
      type(nwt_T) :: BestAJ             ! AJ at Best Fnorm so far.
      real(r8) :: Cosine                ! Of an angle between two vectors
      ! Dump switches
      logical :: D_ATB    ! 'atb' A^T b -- RHS of system
      logical :: D_Col    ! 'col' Column scale vector
      logical :: D_Cov    ! 'cov' Solution covariance
      logical :: D_Diag   ! 'diag' Diagonal of normal equations factor
      logical :: D_Drmc   ! 'drmc' Dump Retriever Matrices Clean
      logical :: D_Dvec   ! 'dvec' DX vector
      integer :: D_DXn    ! 'dxn' DX vector norm if 0, quantity norms if > 0
      logical :: D_Fac_F  ! 'FAC' Full normal equations factor (if you dare)
      logical :: D_Fac_N  ! 'fac' L_Infty norms of blocks of normal equations factor
      integer :: D_Fnmin  ! 'fnmin' terms
      integer :: D_Fnorm  ! 'fnorm' Contributions to | F |
      logical :: D_Gvec   ! 'gvec' Gradient vector
      integer :: D_Jac_F  ! 'JAC' Full Jacobian (if you dare)
      logical :: D_Jac_N  ! 'jac' L_Infty norms of Jacobian blocks
      integer :: D_KTK    ! 'kTk' kTk, which might be normalEquations
      logical :: D_Mas    ! 'mas' Announce master triggering slaves
      logical :: D_Mst    ! 'mst' Block causing factoring abnormal status
      integer :: D_Ndb    ! 'ndb[n]' Newton method debugging output level
      integer :: D_Neq    ! 'neq' > 0 Full normal equations (if you dare)
                          !       = 0 L_Infty norms of Normal equation blocks
      logical :: D_Nin    ! 'nin' Newton method's internal output
      logical :: D_Nwt    ! 'nwt' Commands from Newton method
      logical :: D_Reg    ! 'reg' Tikhonov regularization
      integer :: D_Resid  ! 'dres' Dump residuals with details = d_resid
      logical :: D_Sca    ! 'sca' Newton method's scalars
      logical :: D_Spa    ! 'spa' Sparsity structure of matrices
      logical :: D_Strb   ! 'strb' State vector iff we get into trouble
      logical :: D_Svec   ! 'svec' Final state vector
      logical :: D_Vir    ! 'vir' Virgin matrix, before KTK
      logical :: D_Xvec   ! 'xvec' Current state vector

      real(rk) :: DOF                   ! Degrees of freedom, rows - columns
      type(matrix_Cholesky_T) :: Factored ! Cholesky-factored normal equations
      type (ForwardModelStatus_T) :: FmStat ! Status for forward model
      logical :: FoundBetterState       ! Set if we ever got an nf_best
      type(vector_T) :: FuzzState       ! Random numbers to fuzz the state
      integer :: J, K                   ! Loop inductors and subscripts
      type(matrix_SPD_T), pointer :: KTK ! The Jacobian-derived part of the
                                        ! normal equations.  NormalEquations
                                        ! if no averaging kernel is requested,
                                        ! else kTkSep.
      type(matrix_SPD_T), target :: KTKSep ! The Jacobian-derived part of the
                                        ! normal equations if an averaging kernel
                                        ! is requested.
      type(matrix_T), target :: KTKStar ! The kTkStar contrbution to the averaging
      real(rk) :: Lambda                ! Levenberg-Marquardt stabilization
                                        ! last applied to normal equations
      integer :: LATESTMAFSTARTED       ! For FWMParallel stuff
      real(r8) :: LoopCounter            ! Abandon after 50 * maxJ loop iterations
      integer, dimension(2) :: MATRIXSTATUS ! Flag from matrix calculations
      integer :: Me = -1                ! String index for trace
      real(rv) :: MU                    ! Move Length = scale for DX
      integer, parameter :: NF_GetJ = NF_Smallest_Flag - 1 ! Take an extra loop
                                        ! to get J.
      integer, parameter :: NF_Too_Many = NF_Biggest_Flag + 1 ! Max iterations
      type(matrix_SPD_T), target :: NormalEquations  ! Jacobian**T * Jacobian
      integer :: NumGrad                ! Number of gradient moves
      integer :: NumJ                   ! Number of Jacobian evaluations
      integer :: NumNewt                ! Number of Newton moves
      integer :: NWT_Flag               ! Signal from NWT, q.v., indicating
                                        ! the action to take.
      type(nwt_options) :: NWT_Opt      ! Options for NWT, q.v.
      integer :: PreserveMatrixName     ! Temporary name store
      integer :: Prev_NWT_Flag          ! Previous value of NWT_Flag
      type(vector_T) :: Q               ! Used to calculate Marquardt parameter
      integer :: QTY                    ! Loop counter
      integer :: RowBlock               ! Which block of rows is the forward
                                        ! model filling?
      integer, parameter :: SnoopLevels(NF_DX_AITKEN:NF_FANDJ) = (/ &
      ! dx_aitken dx aitken best gmove newx lev solve evalj evalf
        &      3, 2,     3,   2,    2,   1,  2,   2,    2,    2,  &
      ! start tolx tolx_best tolf too_small fandj
        &  9,   2,        2,   2,        3,    9  /)
      integer :: T                      ! Which Tikhonov: 1 -> V, 2 -> H
      ! real :: T1
      type(matrix_T) :: Temp            ! Because we can't do X := X * Y
      type(matrix_T) :: TempU           ! U**(-1)
      type(matrix_cholesky_T) :: TempC  ! For negateSD caseXoXo
      character(len=10) :: TheFlagName  ! Name of NWTA's flag argument
      integer :: TikhonovRows           ! How many rows of Tiknonov regularization?

      call trace_begin ( me, 'Retrieve.NewtonSolver', cond=.false. )
      ! Set flags from switches
      d_atb = switchDetail ( switches, 'atb' ) > -1
      d_col = switchDetail(switches,'col') > -1
      d_cov = switchDetail(switches,'cov') > -1
      d_diag = switchDetail(switches,'diag') > -1
      d_drmc = switchDetail(switches,'drmc') > -1
      d_dvec = switchDetail(switches,'dvec') > -1
      d_dxn = switchDetail(switches,'dxn')
      d_fac_f = switchDetail(switches,'FAC') > -1
      d_fac_n = switchDetail(switches,'fac') > -1
      d_fnmin = switchDetail ( switches, 'fnmin' )
      d_fnorm = switchDetail ( switches, 'fnorm' )
      d_gvec = switchDetail(switches,'gvec') > -1
      d_jac_f = switchDetail(switches,'JAC')
      d_jac_n = switchDetail(switches,'jac') > -1
      d_kTk = switchDetail ( switches, 'kTk' )
      d_mas = switchDetail ( switches, 'mas' ) > -1
      d_mst = switchDetail(switches,'mst') > -1
      d_ndb = switchDetail(switches,'ndb')
      d_neq = switchDetail(switches,'neq')
      d_nin = switchDetail(switches,'nin') > -1
      d_nwt = switchDetail(switches,'nwt') > -1
      d_reg = switchDetail(switches,'reg') > -1
      d_resid = switchDetail(switches,'dres')
      d_sca = switchDetail(switches,'sca') > -1
      d_spa = switchDetail(switches,'spa') > -1
      d_strb = switchDetail(switches,'strb') > -1
      d_svec = switchDetail(switches,'svec') > -1
      d_vir = switchDetail(switches,'vir') > -1
      d_xvec = switchDetail(switches,'xvec') > -1

      call time_now ( t1 )
      call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
        & ModuleName )
      ! Launch fwmParallel slaves
      if ( parallelMode ) then
        if ( .not. got(f_fwdModelExtra) ) call announceError ( noExtraVec )
        call SetupFWMSlaves ( configDatabase(configIndices), &
        & state, fwdModelExtra, FwdModelOut, jacobian )
      end if
      foundBetterState = ( maxJacobians == 0 )
      chunk%abandoned = .false.
      ! Set options for NWT
      nwt_opt = nwt_options ( relsf=toleranceF, tolxa=toleranceA, &
                            & tolxr=toleranceR, spstrt=initLambda, &
                            & spmini=lambdaMin )
      if ( method == l_newtonian ) then
        call nwt ( nwt_flag, aj, nwt_opt )
      else
        call alt_nwt ( nwt_flag, aj, nwt_opt )
      end if
      ! Create the matrix for the Cholesky factor of the normal equations
      call createEmptyMatrix ( factored%m, &
        & enter_terminal('_factored', t_identifier), &
        & state, state )
      ! Create the normal equations matrix
      call createEmptyMatrix ( normalEquations%m, &
        & enter_terminal ('_normalEquations', t_identifier), state, state )
      ! Create a separate matrix for J^T J in case averaging kernel requested
      call createEmptyMatrix ( kTkSep%m, &
        & enter_terminal ('_kTkSep', t_identifier), state, state )
      ! Create another separate one for J^T(J unmasked) in case 'extended'
      ! averaging kernel required
      call createEmptyMatrix ( kTkStar, &
        & enter_terminal ('_kTkStar', t_identifier), state, state )
      ! Create the vectors we need.
      call copyVector ( v(x), state, vectorNameText='_x', clone=.true. ) ! x = state
      call cloneVector ( v(aTb), v(x), vectorNameText='_ATb' )
      call cloneVector ( v(bestGradient), v(x), vectorNameText='_bestGradient' )
      call cloneVector ( v(bestX), v(x), vectorNameText='_bestX' )
      call cloneVector ( v(candidateDX), v(x), vectorNameText='_candidateDX' )
      call cloneVector ( v(covarianceDiag), v(x), vectorNameText='_covarianceDiag' )
      call cloneVector ( v(dx), v(x), vectorNameText='_DX' )
      call cloneVector ( v(dxUnScaled), v(x), vectorNameText='_DxUnscaled' )
      call cloneVector ( v(f_rowScaled), measurements, vectorNameText='_f_rowScaled' )
      call copyVectorMask ( v(f_rowScaled), measurements )
      call cloneVector ( v(gradient), v(x), vectorNameText='_gradient' )
      call cloneVector ( q, v(x), vectorNameText='Q' )
      call cloneVector ( v(reg_X_x), state, vectorNameText='_reg_X_x' )
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
        if ( d_xvec ) call dump ( v(x), name='Original X' )
        if ( got(f_dumpQuantities) ) &
          & call DumpStateQuantities ( dumpQuantitiesNode, 'Original X' )
      numGrad = 0
      numJ = 0
      numNewt = 0
      aj%axmax = 0.0
        call time_now ( t0 ) ! time base for Newton iteration
      do k = 1, size(v(x)%quantities)
        aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(k)%values)))
      end do

        if ( d_sca ) then
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

      call copyVector ( v(bestX), v(x) ) ! bestX = x to start things off
      prev_nwt_flag = huge(0)
      loopCounter = 0
      atBest = .false.
NEWT: do ! Newton iteration
        loopCounter = loopCounter + 1
        if ( loopCounter > max(50, 50 * maxJacobians) ) then
          chunk%abandoned = .true.
          call writeIntsToChars( L2Options%CurrentChunkNumber, mesgChunkNo )
          mesg = '(' // trim(L2Options%CurrentPhaseName) // ')' // &
            & ' (' // trim(mesgChunkNo) // ') ' // &
            & 'Retrieval abandoned because DNWT appears to be looping.'
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            &  trim(mesg) )
          exit NEWT
        end if
        if ( nwt_flag == nf_getJ ) then
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nwt_flag, &
            & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
        else ! not taking a special iteration to get J
            if ( d_nin )aj%o%k1it = 1 ! Turn on NWTA's internal output
          if ( method == l_newtonian ) then
            call nwta ( nwt_flag, aj )
          else
            call alt_nwta ( nwt_flag, aj )
          end if
        end if
          if ( d_nwt ) then
            call FlagName ( nwt_flag, theFlagName )
            if ( nwt_flag == nf_getJ ) theFlagName = 'GETJ'
            call output ( 'Newton method flag = ' )
            call output ( trim(theFlagName) )
            call output ( numJ, before=', numJ = ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes', before=' at ', after=' seconds' )
          end if
          if ( nwt_flag == nf_evalj .and. prev_nwt_flag == nf_evalf ) then
            prev_nwt_flag = nf_evalj
            cycle
          end if
        select case ( nwt_flag )
        ! The EVALF case is now subsumed into the EVALJ case, so that FNORM
        ! is consistent between the two.  Otherwise, we would need to do all
        ! of the Tikhonov computations here, too.  WVS 2002/09/13.
        case ( nf_evalf, nf_evalj, nf_getJ ) ! EVALF, EVALJ, GETJ  .....
        !{IF ( too many Jacobian values ) EXIT \\
        ! Compute the Jacobian matrix $\bf K$ if you didn't do it when NFLAG
        ! was NF\_EVALF:\\
        ! ${\bf K}_{k,l} = \frac{\partial {\bf f}_k}{\partial {\bf x}_l},
        ! \, k = 1 .. \text{NF},\, l = 1 .. \text{NX}$
        !
        ! Actually, the matrix for the least-squares problem has four parts:
        !  \begin{xalignat*}{2}
        !  {\bf K \delta x} & \simeq {\bf -m}  && \text{The Jacobian} \\
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
        ! {\bf J \delta x} =
        ! \begin{array}{ll}
        !  \left [
        !  \begin{array}{l}
        !  {\bf K} \\
        !  {\bf F} \\
        !  {\bf R} \\
        !  \lambda {\bf I}
        !  \end{array} \right ] {\bf \delta x}
        !  \simeq {\bf -f} =
        !  \left [
        !  \begin{array}{l}
        !   {\bf -m} \\
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
        ! Triangularize ${\bf J}$, and compute the (negative of the) gradient =
        ! $-{\bf J^T f}$.  This is the RHS of the normal equations ${\bf J^T J
        ! \delta \hat x} = -{\bf J^T f}$ where ${\bf \delta \hat x}$ is a
        ! ``Candidate DX'' that might not actually get used.
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
            if ( d_xvec ) call dump ( v(x), details=2, name='State' )
            if ( got(f_dumpQuantities) ) &
              & call DumpStateQuantities ( dumpQuantitiesNode, 'State' )
          if ( numJ > maxJacobians .and. nwt_flag /= nf_getJ ) then
              if ( d_nwt ) then
                call output ( numJ, before= &
                  & 'Newton iteration terminated because Jacobian evaluations (' )
                call output ( maxJacobians, before=') > maxJacobians (', &
                  & after=')', advance='yes' )
                call time_now ( t3 )
                call output ( t3-t0, before=' at ', after=' seconds', advance='yes' )
              end if
            ! Decide whether we need a special iteration just to get a new
            ! Jacobian.
            if ( foundBetterState .and. .not. atBest ) then
              ! Restore BestX, run the forward model one more time to get a new
              ! Jacobian, and form normal equations -- the last two so that the
              ! a posteriori covariance is consistent with BestX.
              aj = bestAJ
              call copyVector ( v(x), v(bestX) ) ! x = bestX
              nwt_flag = nf_getJ
              cycle NEWT
            end if
            if ( ( got(f_outputCovariance) .or. got(f_outputSD) .or. &
              &    got(f_average) ) .and. tikhonovNeeded .and. covSansReg ) then
              ! Get a fresh Jacobian, uncontaminated by Tikhonov regularization,
              ! to compute a posteriori covariance, standard deviation, or an
              ! averaging kernel.  This is hopefully an unusual request.  If it
              ! becomes usual, this should be done by calculating the
              ! Tikhonov regularization matrix, multiplying it (on the left) by
              ! its transpose, and subtracting the result from the normal
              ! equations, rather than by getting a new Jacobian, adding the
              ! a priori covariance, and forming normal equations.  Care will
              ! be needed to handle the cases of column scaling before and after
              ! applying Tikhonov regularization correctly.
              nwt_flag = nf_getJ
              cycle NEWT
            end if
            exit Newt
          end if

          lambda = 0.0
          update = got(f_apriori)
          if ( update ) then ! start normal equations with apriori
            ! Destroy (i.e. clean up) the previous contents of
            ! AprioriMinusX (if any), so as not to have a memory
            ! leak.  Then make it look like apriori.
!           v(aprioriMinusX) = apriori - v(x) ! leaks memory, so do:
            call copyVector ( v(aprioriMinusX), apriori, clone=.true., &
              & vectorNameText='_aprioriMinusX' ) ! a-x := a
            call subtractFromVector ( v(aprioriMinusX), v(x) ) ! a-x := a - x
            if ( tikhonovApriori ) &
              & call copyVector ( v(reg_RHS), v(aprioriMinusX), clone=.true., &
              & vectorNameText='_regRHS' )
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
            ! ( a - x_n )^T F^T F ( a - x_n ) = ( a - x_n )^T \left [
            ! C ( a - x_n ) \right ] }$
            aj%fnorm = v(aprioriMinusX) .mdot. v(covarianceXapriori)
              if ( d_fnorm > -1 ) then             
                call output ( aj%fnorm, &
                  & before='A priori contribution to | F |^2 = ', advance='yes' )
                if ( d_fnorm > 1 ) then
                  call dump ( v(aprioriMinusX), details=d_fnorm )
                  call dump ( v(covarianceXapriori), details=d_fnorm )
                end if
              end if
            !{ Using Apriori requires adding equations of the form ${\bf F
            !  x}_{n+1} \simeq {\bf F a}$ where ${\bf F}$ is the Cholesky
            !  factor of $\bf C$, the inverse of the apriori covariance ${\bf
            !  S}_a$.  So that all of the parts of the problem are solving
            !  for ${\bf\delta x}$, we subtract ${\bf F x}_n$ from both
            !  sides, to get ${\bf F \delta x} \simeq {\bf F (a - x}_n)$.
            !
            ! Actually, we start by putting this into the normal equations,
            ! so we form ${\bf C \delta x} = {\bf C (a - x}_n)$
            if ( diagonal ) then
            ! Iterate with the diagonal of the apriori covariance, then
            ! use the full apriori covariance (one hopes only for one
            ! more iteration). This improves sparsity during iteration.
              call clearMatrix ( normalEquations%m )
              call updateDiagonal ( normalEquations, v(covarianceDiag) )
            else
              call copyMatrixValue ( normalEquations%m, covariance%m )
            end if
            call copyVector ( v(aTb), v(covarianceXApriori) ) ! A^T b := C (a-x_n)
            if ( got(f_aprioriScale) ) then
              call scaleMatrix ( normalEquations%m, aprioriScale )
              ! Don't need to scale aTb -- covarianceXApriori is scaled already
            end if
          else
            aj%fnorm = 0.0_r8
              if ( d_fnorm > -1 ) call output ( &
                & 'No a priori so setting | F | to zero.', advance='yes' )
            call clearMatrix ( normalEquations%m ) ! start with zero
            call clearVector ( v(aTb) ) ! Clear the RHS vector
          end if

          ! Add Tikhonov regularization if requested.
          if ( tikhonovNeeded .and. tikhonovBefore .and. &
            & ( ( nwt_flag /= nf_getj ) .or. .not. covSansReg ) ) then
              call add_to_retrieval_timing( 'newton_solver', t1 )

              call applyTikhonov ( .false., normalEquations, AJ, &
                & tikhonovRows, d_reg, d_fnorm )

          else
            tikhonovRows = 0
          end if

          ! Copy the normal equations formed thus far into a temporary matrix
          ! for the negateSD case
          if ( negateSD .or. associated(aprioriFraction) ) &
            & call copyMatrix ( aPlusRegNEQ%m, normalEquations%m )

          fmStat%maf = 0
          fmStat%newScanHydros = .true.

          if ( d_vir ) call dump ( normalEquations%m, details=2 )

          ! Include the part of the normal equations due to the Jacobian matrix
          ! and the measurements
          call clearVector ( v(f_rowScaled) )
          if ( got(f_average) .and. .not. extendedAverage ) then
            ! Need a separate matrix for K^T K
            kTk => kTkSep
            call clearMatrix ( kTkSep%m )
          else
            kTk => normalEquations
          end if

          ! If in fwm parallel mode, get all slaves computing the forward models
          if ( parallelMode ) then
              if ( d_mas ) &
                & call output ( "Triggering slave forward model runs", advance='yes' )
            latestMAFStarted = min ( parallel%noFwmSlaves, &
              & chunk%lastMAFIndex-chunk%firstMAFIndex+1 )
            do t = 1, latestMAFStarted
              call TriggerSlaveRun ( v(x), t )
            end do
            ! Get the first one setup to send results as soon as possible
            call RequestSlavesOutput ( 1 )
          end if

          ! Loop over MAFs
          do while (fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1)
              call add_to_retrieval_timing( 'newton_solver', t1 )
            !??? What if one config set finished but others still had more
            !??? to do? Ermmm, think of this next time.
            fmStat%maf = fmStat%maf + 1
            fmstat%rows = .false.
            if ( parallelMode ) then
              call ReceiveSlavesOutput ( v(f_rowScaled), fmStat, jacobian )
              ! If there is still another MAF to launch off, let's do that
              if ( latestMAFStarted < chunk%lastMAFIndex-chunk%firstMAFIndex+1 ) then
                latestMAFStarted = latestMAFStarted + 1
                call TriggerSlaveRun ( v(x), latestMAFStarted )
              end if
              ! Set the next slave about packing up its output while we do
              ! our normal equations stuff
              if ( fmStat%maf < chunk%lastMAFIndex-chunk%firstMAFIndex+1 ) &
                & call RequestSlavesOutput ( fmStat%maf + 1 )
            else
              ! Otherwise, we call the forward model as usual
              do k = 1, size(configIndices)
                call forwardModel ( configDatabase(configIndices(k)), &
                  & v(x), fwdModelExtra, v(f_rowScaled), fmStat, Jacobian, &
                  & Hessian, vectorDatabase )
              end do ! k
              call time_now ( t1 )
              ! Forward model calls add_to_retrieval_timing
            end if
            call time_now ( t1 )
              if ( d_jac_f > -1 ) &
                & call dump ( jacobian, name='Jacobian', details=d_jac_f, clean=d_drmc )
            do rowBlock = 1, size(fmStat%rows)
              if ( fmStat%rows(rowBlock) ) then
                 ! Store what we've just got in v(f) ie fwdModelOut
                call copyVector ( v(f), v(f_rowScaled), & ! v(f) := v(f_rowScaled)
                  & quant=jacobian%row%quant(rowBlock), &
                  & inst=jacobian%row%inst(rowBlock) )
                call subtractFromVector ( v(f_rowScaled), measurements, &
                  & quant=jacobian%row%quant(rowBlock), &
                  & inst=jacobian%row%inst(rowBlock) ) ! f - y
                  if ( d_resid > -1 ) call dump ( v(f_rowScaled), &
                    & name='Residuals', details=d_resid )
                !{Let $\bf W$ be the Cholesky factor of the inverse of the
                ! measurement covariance ${\bf S}_m$ (which in our case is
                ! diagonal), i.e. ${\bf W}^T {\bf W} = {\bf S}_m^{-1}$. Row
                ! scale the part of the least-squares problem that arises
                ! from the measurements, i.e. the least-squares problem
                ! becomes $\mathbf{W J \delta \hat x \simeq - W f}$ (actually,
                ! we only row scale ${\bf f}$ here, and scale ${\bf J}$
                ! below).
                if ( got(f_measurementSD) ) then
                  call multiply ( v(f_rowScaled), v(weight), &
                    & quant=jacobian%row%quant(rowBlock), &
                    & inst=jacobian%row%inst(rowBlock) )
                end if
              end if
            end do
            call scaleVector ( v(f_rowScaled), -1.0_r8 ) ! y - f

            if ( got(f_sparseQuantities) ) call Sparsify ( jacobian, &
              & rowQuantities=sparseQuantities )
            if ( got(f_measurementSD) ) call rowScale ( v(weight), jacobian )

            !{Form normal equations:
            ! ${\bf J}^T {\bf W}^T {\bf W J \delta \hat x} =
            ! {\bf J}^T {\bf S}_m^{-1} {\bf J \delta \hat x} =
            ! -{\bf J}^T {\bf W}^T {\bf W f} =
            ! -{\bf J}^T {\bf S}_m^{-1} {\bf f}$.
              call add_to_retrieval_timing( 'newton_solver', t1 )
            call formNormalEquations ( jacobian, kTk, rhs_in=v(f_rowScaled), &
              & rhs_out=v(aTb), update=update, useMask=.true. )
              if ( d_kTk>-1 ) call dump ( kTk%m, name='kTk', details=d_kTk )
              if ( d_atb ) call dump ( v(aTb) )
              call add_to_retrieval_timing( 'form_normeq', t1 )
            if ( got(f_average) .and. extendedAverage ) then
              call MultiplyMatrix_XTY ( jacobian, jacobian, kTkStar, update=update, &
                & maskX=.true., maskY=.false. )
            end if
            update = .true.
              if ( d_jac_n ) &
                call dump_Linf ( jacobian, 'L_infty norms of Jacobian blocks:' )
              if ( d_spa ) call dump_struct ( jacobian, &
                  & 'Sparseness structure of Jacobian blocks:' )
            call clearMatrix ( jacobian )  ! free the space
          end do ! mafs
          call DestroyForwardModelIntermediate ! in case scan model got used
          fmStat%newScanHydros = .true.

          ! In the case where we're forming a regular averaging kernel
          ! (not extended), we need to add this MAFs worth of kTk into the
          ! normal equations.  Otherwise normalEquations and kTk are in
          ! fact pointers to the same matrix
          if ( got(f_average) .and. .not. extendedAverage ) then
            call addToMatrix ( normalEquations%m, kTk%m )
          end if

          ! aj%fnorm is still the square of the function norm
          aj%fnorm = aj%fnorm + ( v(f_rowScaled) .mdot. v(f_rowScaled) )
            if ( d_fnorm > -1 ) &
              & call output ( v(f_rowScaled) .mdot. v(f_rowScaled), &
                & before='Measurement contribution to | F |^2 = ', advance='yes' )

          ! Add Tikhonov regularization if requested.  We do it here instead
          ! of before adding the Jacobian so that we can scale it up by the
          ! column scaling before scaling the normal equations back down
          ! by the column scaling.  The effect is that regularization takes
          ! place (roughly) on the column-scaled problem, and if scaling by
          ! the column norm is requested, it's still makes the column norms
          ! all 1.0.
          if ( tikhonovNeeded .and. .not. tikhonovBefore .and. &
            & ( ( nwt_flag /= nf_getj ) .or. .not. covSansReg ) ) then
              call add_to_retrieval_timing( 'newton_solver', t1 )

              call applyTikhonov ( .true., normalEquations, AJ, &
                & tikhonovRows, d_reg, d_fnorm )

          else
            tikhonovRows = 0
          end if

          !{Column Scale $\bf J$ (row and column scale ${\bf J}^T {\bf J}$)
          ! using the matrix $\bf\Sigma$. We cater for several choices of
          ! $\bf\Sigma$.  After row scaling with $\bf W$ and column scaling
          ! with $\bf\Sigma$, the least-squares problem is ${\bf W J \Sigma
          ! \Sigma}^{-1} {\bf \delta \hat x} \simeq -{\bf W f}$.  We will
          ! solve for ${\bf \Sigma}^{-1} {\bf \delta \hat x}$, and then for
          ! $\bf \delta \hat x$.
          select case ( columnScaling )
          case ( l_apriori ) ! v(columnScaleVector) := apriori
            call copyVector ( v(columnScaleVector), apriori )
          case ( l_covariance )
            !??? Can't get here until allowed by init_tables
          case ( l_norm ) ! v(columnScaleVector) := diagonal of normal
                          !                         equations = column norms^2
            call getDiagonal ( normalEquations%m, v(columnScaleVector) )
          end select
          if ( columnScaling /= l_none ) then ! Compute $\Sigma$
            do j = 1, v(columnScaleVector)%template%noQuantities
              where ( v(columnScaleVector)%quantities(j)%values <= 0.0 )
                v(columnScaleVector)%quantities(j)%values = 1.0
              elsewhere
                v(columnScaleVector)%quantities(j)%values = 1.0 / &
                  & sqrt( v(columnScaleVector)%quantities(j)%values )
              end where
            end do
              if ( d_col ) &
                & call dump ( v(columnScaleVector), name='Column scale vector' )
            !{Scale in normal equations form: ${\bf \Sigma}^T {\bf  J}^T {\bf
            ! W}^T {\bf W J \Sigma \Sigma}^{-1} {\bf \delta \hat x} = {\bf
            ! \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf J \Sigma \Sigma}^{-1}
            ! {\bf \delta \hat x} =  -{\bf \Sigma}^T {\bf J}^T {\bf W}^T {\bf
            ! W f} = -{\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$.
            call columnScale ( normalEquations%m, v(columnScaleVector) )
            call rowScale ( v(columnScaleVector), normalEquations%m )
            call multiply ( v(aTb), v(columnScaleVector) )
          end if

          !{The (negative of the) gradient $= -{\bf \Sigma}^T {\bf J}^T
          ! {\bf S}_m^{-1} {\bf f}$
          ! is the right-hand side of the normal equations.  Save it.
          call copyVector ( v(gradient), v(aTb) ) ! gradient := atb
          ! Make sure the gradient is zerod out in various places.
          ! aTb may have non zero values in masked places where x/=a.
          call copyVectorMask ( v(gradient), v(x) )
          call clearUnderMask ( v(gradient) )
            if ( d_gvec ) call dump ( v(gradient), name='gradient' )

          !{Compute the Cholesky factor of the LHS of the normal equations:
          ! ${\bf U}^T {\bf U} {\bf \Sigma}^{-1} {\bf \delta \hat x} =
          ! ({\bf\Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf J \Sigma})
          ! {\bf \Sigma}^{-1} {\bf \delta \hat x} = -{\bf \Sigma}^T {\bf J}^T
          ! {\bf S}_m^{-1} {\bf f}$.

          ! factored%m%block => normalEquations%m%block ! to save space
          ! Can't do the above because we need to keep the normal
          ! equations around, in order to subtract Levenberg-Marquardt and
          ! apriori covariance, in order to compute a posteriori covariance
            if ( d_neq == 0 ) &
              call dump_Linf ( normalEquations%m, &
                & 'L_infty norms of Normal Equations blocks after scaling:', &
                & upper=.true. )
            if ( d_spa ) call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
            call add_to_retrieval_timing ( 'newton_solver', t1 )
          call choleskyFactor ( factored, normalEquations, status=matrixStatus )
            call add_to_retrieval_timing ( 'cholesky_factor', t1 )
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            if ( d_mst ) then
              ! Dump the structure if we haven't done so already
              if ( .not. d_spa ) call dump_struct ( normalEquations%m, &
                  & 'Sparseness structure of Normal equations blocks:' )
!               call dump ( normalEquations%m%block(matrixStatus(1),matrixStatus(1)), &
!               & name='Offending block', details=9 )
              call dump ( normalEquations%m, name='Normal equations', details=9 )
              call dump ( v(columnScaleVector), name='Column scale', details=9 )
              call dump ( v(x), name='Current state', details=9 )
              if ( got(f_lowBound) ) call dump ( lowBound, name='Low Bound', details=9 )
              if ( got(f_highBound) ) call dump ( highBound, name='High Bound', details=9 )
            else
              call dump ( normalEquations%m%block(matrixStatus(1),matrixStatus(1)), &
                & name='Offending block', details=0 )
            end if
            call output ( &
              & 'Consider applying or increasing the weight of Tikhonov regularization for the block', &
              & advance='yes' )
            call output ( 'Re-run with -Smst to see more details', advance='yes' )
            call newtonSolverFailed ( 'problem factoring normalEquations in evalJ', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
            if ( d_neq > 0 ) &
              & call dump ( normalEquations%m, 'Normal Equations', 2, clean=d_drmc )
            if ( d_diag ) then
              call getDiagonal ( factored%m, v(dxUnscaled) )
              call dump ( v(dxUnscaled), &
                & name='Diagonal of factored normal equations:' )
            end if
            if ( d_fac_n ) call dump_Linf ( factored%m, &
              & 'L_infty norms of blocks of factor:', upper=.true. )
            if ( d_spa ) call dump_struct ( factored%m, &
              & 'Sparseness structure of blocks of factor:', upper=.true. )
            if ( d_fac_f ) call dump ( factored%m, 'Factor', 2, clean=d_drmc )

          ! Compute numbers of rows and columns of Jacobian actually used. 
          ! Don't count rows due to Levenberg-Marquardt stabilization.  Do
          ! count rows due to a priori or regularization.  Don't count
          ! masked-off measurement or state rows or columns.
          jacobian_cols = sum(jacobian%col%nelts)
          do j = 1, state%template%noQuantities
            if ( associated(state%quantities(j)%mask) ) &
              ! Subtract masked-off columns
              & jacobian_cols = jacobian_cols - &
              &   countBits(state%quantities(j)%mask, what=m_linAlg )
          end do
          jacobian_rows = sum(jacobian%row%nelts)
          do j = 1, measurements%template%noQuantities
            if ( associated(measurements%quantities(j)%mask) ) &
              ! subtract masked-off rows
              & jacobian_rows = jacobian_rows - &
              &   countBits(measurements%quantities(j)%mask, what=m_linAlg )
          end do

          ! Correct the number of rows for apriori information.  Note that
          ! there is an approximation here: We don't take any account of
          ! whether the a priori is used on an element by element basis.
          if ( got(f_apriori) ) jacobian_rows = jacobian_rows + jacobian_cols

          if ( tikhonovNeeded ) then
            ! Correct the number of rows for Tikhonov information, unless
            ! we're taking a special iteration to get a Jacobian without
            ! Tikhonov regularization.
            if ( nwt_flag /= nf_getJ .or. .not. covSansReg ) &
              & jacobian_rows = jacobian_rows + tikhonovRows
          end if

          ! Exit if we are taking a special iteration to get J.
          if ( nwt_flag == nf_getJ ) exit NEWT

          aj%diag = minDiag ( factored ) ! element on diagonal with
            !       smallest absolute value, after triangularization
          aj%ajn = maxL1 ( factored%m ) ! maximum L1 norm of
          !       column in upper triangle after triangularization
            call add_to_retrieval_timing ( 'newton_solver', t1 )

          !{Solve ${\bf U}^T {\bf U} {\bf \Sigma}^{-1} {\bf \delta \hat x} =
          ! -{\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$ for
          ! ${\bf U} {\bf \Sigma}^{-1} {\bf \delta \hat x} = {\bf v}$.

          call solveCholesky ( factored, v(candidateDX), v(aTb), &
            & transpose=.true., status=matrixStatus )
          call copyVectorMask ( v(candidateDX), v(x) )
            call add_to_retrieval_timing ( 'cholesky_solver', t1 )
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            if ( d_mst ) then
              call dump ( normalEquations%m%block(matrixStatus(1),matrixStatus(1)), &
              & name='Offending block', details=9 )
              call dump ( v(columnScaleVector), name='Column scale', details=9 )
              call dump ( v(x), name='Current state', details=9 )
              if ( got(f_lowBound) ) call dump ( lowBound, name='Low Bound', details=9 )
              if ( got(f_highBound) ) call dump ( highBound, name='High Bound', details=9 )
            end if
            call newtonSolverFailed ( 'problem solving normalEquations in evalJ', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
            if ( d_dvec ) call dump ( v(candidateDX), name='CandidateDX 1' )
            if ( d_dxn > -1 ) &
              & call dumpVectorNorms ( v(candidateDX), d_dxn, &
                & name='|CandidateDX| 1', useMask=.true. )

          !{AJ\%FNMIN = $L_2$ norm of residual, $||{\bf\Sigma}^T {\bf J}^T
          ! {\bf S}_m^{-1} {\bf J \Sigma \Sigma}^{-1} {\bf \delta \hat x} +
          ! {\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}||$ where
          ! $\bf\delta \hat x$ is the ``Candidate DX'' that might not get
          ! used. This can be gotten without saving $\bf J$ as $\sqrt{{\bf f}^T
          ! {\bf f} - {\bf v}^T {\bf v}}$ where ${\bf v} = {\bf U}^{-T}
          ! {\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$. The variable
          ! {\tt v(candidateDX)} is a temp here = $\bf v$ in {\tt wvs-011}
          ! or $\bf y$ in {\tt wvs-014}.
          aj%fnmin = aj%fnorm - (v(candidateDX) .mdot. v(candidateDX))
          if ( aj%fnmin < 0.0 ) then
            if ( abs(aj%fnmin) > &
               & (aj%fnorm + (v(candidateDX) .mdot. v(candidateDX))) * &
                   & sqrt(sqrt(epsilon(1.0e0)))**3 ) then
              ! Relative to aj%fnorm and v(candidateDX)**2, it's more than
              ! single-precision round_off**(0.75).
              call output ( aj%fnmin, before='How can aj%fnmin be negative?  aj%fnmin = ', &
                & advance='yes' )
              call output ( aj%fnorm, before='aj%fnorm**2 = ', advance='yes' )
              call output ( v(candidateDX) .mdot. v(candidateDX), &
                & before='norm(candidateDX)**2 = ', advance='yes' )
              if ( d_fnmin > 1 ) &
                & call dump ( v(candidateDX), name='v(candidateDX)', details=d_fnmin )
              if ( aj%fnmin < -10.0 * aj%fnorm * sqrt(epsilon(aj%fnorm)) ) &
                & call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & "Norm of residual not in Jacobian's column space is imaginary!" )
            end if
            aj%fnmin = tiny ( aj%fnmin )
          end if
          ! Compute the normalised chiSquared statistics etc.
          dof = max ( jacobian_rows - jacobian_cols, 1 )
          chiSqMinNorm = aj%fnmin / dof
          chiSqNorm = aj%fnorm / dof
          aj%fnmin = sqrt(aj%fnmin)
          aj%fnorm = sqrt(aj%fnorm)
          call copyVectorMask ( v(gradient), v(x) )
          aj%gradn = sqrt(v(gradient) .mdot. v(gradient)) ! L2Norm(gradient)
            if ( d_sca ) then
              call dump ( (/ aj%fnorm, aj%fnmin, aj%diag, aj%ajn, aj%gradn /), &
                & '     | F |      aj%fnmin       aj%diag      L1| FAC |        | G |', &
                & options='c' )
              call dump ( (/ chiSqNorm, chiSqMinNorm, dof /), &
                & ' chi^2/DOF  chimin^2/DOF           DOF ', options='c' )
            end if
        case ( nf_lev ) ! ..................................  LEV  .....
        !{Solve ${\bf U}^T {\bf q} = {\bf \delta \hat x}$, then compute
        ! $||{\bf q}||^2$, which is used to compute the Levenberg-Marquardt
        ! stabilization parameter.
          call solveCholesky ( factored, q, v(candidateDX), &
            & transpose=.true., status=matrixStatus ) ! q := factored^{-T} v(candidateDX)
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            call newtonSolverFailed ( 'problem solving for q in levenberg', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
          call copyVectorMask ( q, v(x) )
          aj%qnsq = q .mdot. q
        case ( nf_solve ) ! ..............................  SOLVE  .....
        !{Apply Levenberg-Marquardt stabilization with parameter
        ! $\lambda =$ {\bf AJ\%SQ}.  I.e., form $({\bf \Sigma}^T {\bf J}^T
        ! {\bf W}^T {\bf W J \Sigma + \lambda^2 I}) {\bf \Sigma}^{-1}
        ! {\bf \delta \hat x = \Sigma}^T {\bf J}^T {\bf W}^T {\bf f}$ for
        ! ${\bf \Sigma}^{-1} {\bf \delta \hat x}$.  The update is
        ! aj\%sq$^2 - \lambda^2$ to remove the previous Levenberg-Marquardt
        ! parameter.  Then set
        ! \begin{description}
        !   \item[AJ\%FNMIN] as for NWT\_FLAG = NF\_EVALJ, but taking
        !     account of Levenberg-Marquardt stabilization;
        !   \item[AJ\%DXN] = L2 norm of ``candidate DX'';
        !   \item[AJ\%GDX] = (Gradient) .dot. (``candidate DX'')
        ! \end{description}
          call updateDiagonal ( normalEquations, &
                              & max(aj%sq, lambdaMin)**2 - lambda**2 )
          lambda = max(aj%sq, lambdaMin)
          ! factored%m%block => normalEquations%m%block ! to save space
          ! Can't do the above because we need to keep the normal equations
          ! around, in order to subtract Levenberg-Marquardt and apriori
          ! covariance, in order to compute a posteriori covariance
            if ( d_neq > 0 ) then
              call dump ( normalEquations%m, 'Normal Equations', 2, clean=d_drmc )
            else if ( d_neq == 0 ) then
              call dump_Linf ( normalEquations%m, &
                & 'L1 norms of Normal Equations blocks after Marquardt:', &
                & upper=.true. )
            end if
            if ( d_spa ) call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
            call add_to_retrieval_timing( 'newton_solver', t1 )
          call choleskyFactor ( factored, normalEquations, status=matrixStatus )
            call add_to_retrieval_timing( 'cholesky_factor', t1 )
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            if ( d_mst ) then
              call dump ( normalEquations%m%block(matrixStatus(1),matrixStatus(1)), &
              & name='Offending block', details=9 )
              call dump ( v(columnScaleVector), name='Column scale', details=9 )
              call dump ( v(x), name='Current state', details=9 )
              if ( got(f_lowBound) ) call dump ( lowBound, name='Low Bound', details=9 )
              if ( got(f_highBound) ) call dump ( highBound, name='High Bound', details=9 )
            end if
            call newtonSolverFailed ( 'problem factoring normal equations in solve', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
            if ( d_fac_n ) call dump_Linf ( factored%m, &
                & 'L1 norms of blocks of factor after Marquardt:', &
                & upper=.true. )
            if ( d_spa ) call dump_struct ( factored%m, &
                & 'Sparseness structure of blocks of factor:', upper=.true. )
            if ( d_fac_f ) call dump ( factored%m, &
                & name='Factored, after Marquardt', details=2, clean=d_drmc )

          !{Solve ${\bf U}^T {\bf U} {\bf \Sigma}^{-1} {\bf \delta \hat x} =
          ! ({\bf\Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf J \Sigma}) {\bf
          ! \Sigma}^{-1} {\bf \delta \hat x} = -{\bf \Sigma}^T {\bf J}^T {\bf
          ! S}_m^{-1} {\bf f}$ for ``candidate DX'' = ${\bf \Sigma}^{-1} {\bf
          ! \delta \hat x}$ using two back solves.  First solve ${\bf U}^T {\bf v}
          ! = -{\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$ for ${\bf v}$,
          ! then solve ${\bf U \Sigma}^{-1} {\bf \delta \hat x}$ = ${\bf v}$ for
          ! ${\bf \Sigma}^{-1} {\bf \delta \hat x}$. Meanwhile, set AJ\%FNMIN =
          ! $\sqrt{{\bf f}^T {\bf f} - {\bf v}^T {\bf v}}$ as when NWT\_FLAG is
          ! NF\_EVALJ, but taking account of Levenberg-Marquardt stabilization.
            call add_to_retrieval_timing( 'newton_solver', t1 )
          call solveCholesky ( factored, v(candidateDX), v(aTb), &
            & transpose=.true., status=matrixStatus ) ! v(candidateDX) := factored^{-T} v(aTb)
          call copyVectorMask ( v(candidateDX), v(x) )
            call add_to_retrieval_timing( 'cholesky_solver', t1 )
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            call newtonSolverFailed ( 'problem solving for aTb in solve', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
            if ( d_dvec ) call dump ( v(candidateDX), name='CandidateDX 2' )
            if ( d_dxn > -1 ) &
              & call dumpVectorNorms ( v(candidateDX), d_dxn, &
                & name='|CandidateDX| 2', useMask=.true. )

          if ( ieee_is_nan ( aj%fnorm ) ) then
            if ( d_strb ) call dump ( v(x), name='Current state', details=9 )
            call newtonSolverFailed ( 'numerical problems (with radiances?)', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if

          ! aj%fnorm is now the norm of f, not its square.
          ! The following calculation of fnmin was commented out, but on
          ! 11 September 2002 FTK told me it's the right thing to do
          ! after all.
          aj%fnmin = aj%fnorm**2 - (v(candidateDX) .mdot. v(candidateDX))

          if ( ieee_is_nan ( aj%fnmin ) ) then
            if ( d_strb ) call dump ( v(x), name='Current state', details=9 )
            call newtonSolverFailed ( 'numerical problems (with derivatives?)', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if

          if ( aj%fnmin < 0.0 ) then
            if ( abs(aj%fnmin) > &
               & (aj%fnorm + (v(candidateDX) .mdot. v(candidateDX))) * &
                   & sqrt(sqrt(epsilon(1.0e0)))**3 ) then
              ! Relative to aj%fnorm and v(candidateDX)**2, it's more than
              ! single-precision round_off**(0.75).
              call output ( aj%fnmin, before='How can aj%fnmin be negative?  aj%fnmin = ', &
                & advance='yes' )
              call output ( aj%fnorm**2, before='aj%fnorm**2 = ', advance='yes' )
              call output ( v(candidateDX) .mdot. v(candidateDX), &
                & before='norm(candidateDX)**2 = ', advance='yes' )
              if ( d_fnmin > 1 ) &
                & call dump ( v(candidateDX), name='v(candidateDX)', details=d_fnmin )
              if ( aj%fnmin < -10.0 * aj%fnorm * sqrt(epsilon(aj%fnorm)) ) &
                & call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & "Norm of residual not in Jacobian's column space is imaginary!" )
            end if
            aj%fnmin = tiny ( aj%fnmin )
          end if
          aj%fnmin = sqrt(aj%fnmin)
          ! v(candidateDX) := factored^{-1} v(candidateDX)
          call solveCholesky ( factored, v(candidateDX), status=matrixStatus )
          call copyVectorMask ( v(candidateDX), v(x) )
          if ( any ( matrixStatus /= 0 ) ) then
            call dump ( matrixStatus, 'matrixStatus [block, row in trouble] = ' )
            call newtonSolverFailed ( 'problem solving normal equations in solve', &
              & nwt_flag, prev_nwt_flag, d_ndb, aj )
!             exit NEWT
            cycle NEWT
          end if
            if ( d_dvec ) call dump ( v(candidateDX), name='CandidateDX 3' )
            if ( d_dxn > -1 ) &
              & call dumpVectorNorms ( v(candidateDX), d_dxn, &
                & name='|CandidateDX| 3', useMask=.true. )

          !{Shorten $\delta {\bf \hat x} =$ {\tt dxUnScaled} to stay within
          ! the bounds. We aren't ready to take the move, but NWTA needs an
          ! accurate value for the norm of {\tt candidateDX}.  We don't simply
          ! put out-of-bounds variables onto their bounds, as this could
          ! change the direction away from the Newton direction.
          call copyVector ( v(dxUnScaled), v(candidateDX) ) ! dxUnscaled = dx
          if ( columnScaling /= l_none ) then ! Scale into problem's space
            ! dxUnScaled = dxUnScaled # columnScaleVector:
            call multiply ( v(dxUnScaled), v(columnScaleVector) )
          end if
          mu = 1.0_rv
          if ( got(f_lowBound) ) &
            & call boundMove ( mu, lowBound, v(x), v(dxUnScaled), 'low', muMin )
          if ( got(f_highBound) ) &
            & call boundMove ( mu, highBound, v(x), v(dxUnScaled), 'high', muMin )
          ! dxUnScaled = mu * dxUnScaled -- may shorten vector in problem space
          call scaleVector ( v(dxUnScaled), mu )
          call scaleVector ( v(candidateDX), mu )
          aj%dxn = sqrt(v(candidateDX) .mdot. v(candidateDX)) ! L2Norm(dx)
          call copyVectorMask ( v(gradient), v(x) )
          aj%gdx = v(gradient) .mdot. v(candidateDX)
          if ( .not. aj%starting ) then
            call copyVectorMask ( v(dx), v(x) )
            aj%dxdxl = v(dx) .mdot. v(candidateDX)
          end if
            if ( d_dvec ) call dump ( v(candidateDX), name='CandidateDX 4' )
            if ( d_dxn > -1 ) &
              & call dumpVectorNorms ( v(candidateDX), d_dxn, &
                & name='|CandidateDX| 4', useMask=.true. )
            if ( d_sca ) then
              cosine = -2.0_r8
              if ( aj%dxn > 0.0_r8 .and. aj%gradn > 0.0_r8 ) &
                cosine = aj%gdx/(aj%dxn*aj%gradn)
              call dump ( (/ aj%dxn, aj%fnmin, cosine, aj%sq, mu /), &
                & '    | DX |      aj%fnmin    ' // &
                & 'cos(G, DX)        lambda            mu', options='c' )
              if ( .not. aj%starting ) then
                call output ( ' cos(DX, DXL) = ' )
                cosine = -2.0_r8
                if ( aj%dxn > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
                  cosine = aj%dxdxl / (aj%dxn*aj%dxnl)
                call output ( cosine, format='(1pe14.7)', advance='yes' )
              ! call output ( ',  DX . DXL = ' )
              ! call output ( aj%dxdxl, format='(1pe14.7)', advance='yes' )
              end if
            end if
            if ( d_ndb >= 0 ) call my_nwtdb ( aj, width=9, why='After Solve' )
        case ( nf_newx ) ! ................................  NEWX  .....
          ! Set X = X + DX
          !     AJ%AXMAX = MAXVAL(ABS(X)),
          !     AJ%BIG = ANY ( DX > 10.0 * epsilon(X) * X )

          !{Set $\delta {\bf x} := \delta {\bf \hat x}$ and then ${\bf x} :=
          ! {\bf x} + \delta {\bf x}$.  We already accounted for column scaling
          ! when a candidate Newton move or gradient move was put into {\tt
          ! v(dxUnscaled)}, i.e., we're now dealing with $\delta {\bf \hat x}$,
          ! not ${\bf \Sigma}^{-1} \delta {\bf \hat x}$. Shorten {\tt
          ! dxUnscaled} if necessary to stay within the bounds.  We did this
          ! when we solved for $\delta {\bf \hat x}$ as a Newton move because
          ! DNWT needed the norm of the actual candidate move, but if this move
          ! is a gradient or Aitken-accelerated move, we still need to stay in
          ! bounds.
          mu = 1.0
          if ( got(f_lowBound) ) &
            & call boundMove ( mu, lowBound, v(x), v(dxUnscaled), 'low', muMin )
          if ( got(f_highBound) ) &
            & call boundMove ( mu, highBound, v(x), v(dxUnscaled), 'high', muMin )
          if ( mu < 1.0_rv ) call scaleVector ( v(dxUnscaled), mu )
          call addToVector ( v(x), v(dxUnScaled) ) ! x = x + dxUnScaled
          atBest = .false.
            if ( d_dvec ) call dump ( v(dxUnScaled), name='dX Unscaled' )
            if ( d_dxn > -1 ) &
              & call dumpVectorNorms ( v(dxUnScaled), d_dxn, &
                & name='|dxUnScaled|', useMask=.true. )
          if ( got(f_stateMax) ) then
            do j = 1, size(v(x)%quantities)
              stateMax%quantities(j)%values = max(stateMax%quantities(j)%values, &
                &                                 v(x)%quantities(j)%values)
            end do
          end if
          if ( got(f_stateMin) ) then
            do j = 1, size(v(x)%quantities)
              stateMin%quantities(j)%values = min(stateMin%quantities(j)%values, &
                &                                 v(x)%quantities(j)%values)
            end do
          end if
          aj%axmax = 0.0
          aj%big = .false.
          do j = 1, size(v(x)%quantities)
            aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(j)%values)))
            if ( any( abs(v(dx)%quantities(j)%values) > &
              & 10.0 * epsilon(aj%axmax) * abs(v(x)%quantities(j)%values) ) ) &
              & aj%big = .true.
          end do
            if ( d_sca ) then
              call output ( ' aj%axmax = ' )
              call output ( aj%axmax, format='(1pe14.7)' )
              if ( .not. aj%starting ) then
                call output ( ' cos(DX, DXL) = ' )
                cosine = -2.0_r8
                if ( aj%dxn > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
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
          numGrad = numGrad + 1
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nwt_flag, &
            & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
          call copyVectorMask ( v(bestX), v(x) ) ! Don't lose x's mask
          call copyVector ( v(x), v(bestX) ) ! x = bestX
          if ( .not. aj%starting ) then
            call copyVectorMask ( v(dx), v(x) )
            call copyVectorMask ( v(bestGradient), v(x) )
            aj%dxdxl = aj%gfac * ( v(dx) .mdot. v(bestGradient) )
          end if
          call copyVector ( v(dxUnscaled), v(bestGradient) ) ! dxUnscaled := bestGradient
          if ( columnScaling /= l_none ) & ! dxUnscaled = bestGradient / scale
            ! The saved gradient is in scaled coordinates.
            ! Put it back into the problem coordinates.
            & call multiply ( v(dxUnscaled), v(columnScaleVector) )
          ! dx = aj%gfac * "Best Gradient":
          call scaleVector ( v(dxUnscaled), aj%gfac )
            if ( d_sca ) then
              call output ( ' aj%gfac = ' )
              call output ( aj%gfac, format='(1pe14.7)' )
              call output ( ' |DX| = ' )
              call output ( aj%gradnb * aj%gfac, advance='yes' )
            end if
            if ( d_gvec ) call dump ( v(dxUnscaled), name='Gradient move from best X' )
        case ( nf_best ) ! ................................  Best  .....
        ! Set "Best X" = X, "Best Gradient" = Gradient
          atBest = .true.
          foundBetterState = .true.
          bestAJ = aj
          call copyVector ( v(bestX), v(x) ) ! bestX = x
          call copyVector ( v(bestGradient), v(gradient) ) ! bestGradient = gradient
        case ( nf_aitken ) ! ............................  Aitken  .....
        ! Set DX = DX - "Candidate DX",
        !     AJ%DXDX = dot_product( DX, DX )
        ! IF ( AJ%DXDX /= 0.0 ) &
        !   Set AJ%DXDXL = dot_product( DX, "Candidate DX" )
          call subtractFromVector ( v(dx), v(candidateDX) ) ! dx = dx - candidateDX
          call copyVectorMask ( v(dx), v(x) )
          aj%dxdx = v(dx) .mdot. v(dx)
          if ( aj%dxdx > 0.0 ) then
            call copyVectorMask ( v(candidateDX), v(x) )
            aj%dxdxl = v(dx) .mdot. v(candidateDX)
          end if
            if ( d_dvec ) call dump ( v(dx), name='dx after Aitken' )
            if ( d_sca ) then
              call output ( ' aj%dxdx = ' )
              call output ( aj%dxdx, format='(1pe14.7)' )
              call output ( ', cos(DX, DXL) = ' )
              cosine = -2.0_r8
              if ( aj%dxn > 0.0_r8 .and. aj%dxnl > 0.0_r8 ) &
                cosine = aj%dxdxl / (aj%dxn*aj%dxnl)
              call output ( cosine, format='(1pe14.7)', advance='yes' )
            end if
        case ( nf_dx ) ! ....................................  DX  .....
          numNewt = numNewt + 1
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nwt_flag, &
            & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
          if ( .not. aj%starting ) then
            call copyVectorMask ( v(dx), v(x) )
            call copyVectorMask ( v(candidateDX), v(x) )
            aj%dxdxl = v(dx) .mdot. v(candidateDX)
          end if
          call copyVector ( v(dx), v(candidateDX) ) ! dx = candidateDX
          ! v(dxUnscaled) was computed during nf_solve
        case ( nf_dx_aitken ) ! ......................  DX_Aitken  .....
          numNewt = numNewt + 1
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nwt_flag, &
            & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
!!! This is needed for correct Aitken acceleration, but Bill Read doesn't
!!! like the results:
          ! dxUnscaled = aj%cait * dxUnscaled:
          if ( do_Aitken ) call scaleVector ( v(dxUnscaled), aj%cait )
          ! dx = aj%cait * candidateDX:
          call scaleVector ( v(candidateDX), aj%cait, v(dx) )
          if ( .not. aj%starting ) then
            call copyVectorMask ( v(dx), v(x) )
            call copyVectorMask ( v(candidateDX), v(x) )
            aj%dxdxl = v(dx) .mdot. v(candidateDX)
          end if
            if ( d_sca ) then
              call output ( ' aj%cait = ' )
              call output ( aj%cait, format='(1pe14.7)', advance='yes' )
            end if
            if ( d_dvec ) call dump ( v(dx), name='dx after dx Aitken' )
        case ( nf_tolx, nf_tolx_best, nf_tolf, nf_too_small ) ! ........
          ! IF ( NWT_FLAG == NF_TOO_SMALL ) THEN
          !   Take special action if requested accuracy is critical
          ! END IF
          ! Convergence to desired solution.  Do whatever you want to
          ! with the solution.
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nwt_flag, &
            & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
          if ( nwt_flag == nf_tolx_best ) then
            call copyVector ( v(x), v(bestX) )
            aj = bestAJ
            atBest = .false.
          end if
          if ( .not. diagonal ) then
              if ( d_nwt ) then
                call time_now ( t3 )
                call output ( t3-t0, after=' seconds', advance='yes', &
                  & before='Newton iteration terminated because of convergence at ' )
              end if
            nwt_flag = nf_getJ
            ! Do we need to get a fresh Jacobian, uncontaminated by Tikhonov
            ! regularization and Levenberg-Marquardt stabilization, in order
            ! to compute a posteriori covariance or standard deviation, or
            ! an averaging kernel?
            if ( ( got(f_outputCovariance) .or. got(f_outputSD) .or.  &
              &    got(f_average) ) .and.                             &
              &  ( nwt_flag == nf_tolx_best .or.                      &
              &    nwt_flag == nf_tolx .and. aj%sq /= 0.0 .or.        &
              &    tikhonovNeeded .and. covSansReg ) ) cycle NEWT ! yes
            exit NEWT ! no
          end if
          diagonal = .false.
          nwt_flag = nf_start
        case ( nf_fandj ) ! ...............................  FANDJ .....
          ! There is probably an error in the way F or J is computed.
          ! A warning has been printed by the error processor.
          ! IF ( you have confidence in F and J ) CYCLE
          ! STOP
        end select
      ! IF ( you want to return to a previous best X ) NWT_FLAG = 0
        if ( snoopKey /= 0 .and. snoopLevel >= snoopLevels(nwt_flag) ) then
          call FlagName ( nwt_flag, theFlagName )
          call snoop ( key=snoopKey, vectorDatabase=vectorDatabase, &
            & anotherVectorDatabase=v, &
            & anotherComment=trim(snoopComment) // ': ' // trim(theFlagName), &
            & anotherPhaseName=trim(phaseName), &
            & matrixDatabase=(/ factored%m, normalEquations%m /) )

        end if
          if ( d_ndb >= 2 ) call my_nwtdb ( aj, width=9, why='Bottom' )
        prev_nwt_flag = nwt_flag
      end do NEWT ! Newton iteration

      dof = max ( jacobian_rows - jacobian_cols, 1 )
      chiSqNorm = aj%fnorm / dof
      aj%fnorm = sqrt(aj%fnorm)
        if ( d_sca ) &
          & call dump ( (/ aj%fnorm, chiSqNorm, dof /) , &
          & '     | F |    chi^2/DOF           DOF ', options='c' )

        if ( d_ndb >= 1 ) &
          & call my_nwtdb ( aj, width=9, why = "After Newton iteration" )

      if ( got(f_diagnostics) .and. numJ > maxJacobians ) &
        & call FillDiagVec ( diagnostics, aj, &
        & numGrad=numGrad, numJ=numJ, numNewt=numNewt, nwt_flag=nf_too_many, &
        & jacobian_rows=jacobian_rows, jacobian_cols=jacobian_cols )
      abandoned_value = 0.d0
      if ( chunk%abandoned ) abandoned_value = 1.d0
      if ( got(f_diagnostics) ) then
          call fillDiagQty ( diagnostics,  l_dnwt_abandoned, abandoned_value )
          call fillDiagQty ( diagnostics, l_dnwt_count, loopCounter )
      end if

      ! chunk%abandoned = abandoned
      if ( chunk%abandoned ) then
          if ( d_cov ) call output ( "Abandoned; no covariance calculation", advance="yes" )
        foundBetterState = .false.
        ! We'll output our last good state, but set the error bar to
        ! the a priori error.
        if ( associated(outputSD) .and. associated(covariance) ) &
          & call GetDiagonal ( covariance%m, outputSD, squareRoot=.true., &
          & invert=.true., zeroOK=.true. )
        if ( got(f_outputCovariance) .and. associated(covariance) ) &
          & call CopyMatrix ( outputCovariance%m, covariance%m )
        ! Also flag our measurements as never having been used
        ! Was going to flag the forward model output, but its mask
        ! Gets overwritten with the measurements one.
        do j = 1, measurements%template%noQuantities
          if ( .not. associated ( measurements%quantities(j)%mask ) ) &
            & call CreateMask ( measurements%quantities(j) )
          measurements%quantities(j)%mask = char ( ior ( ichar ( &
            & measurements%quantities(j)%mask ), m_linAlg ) )
        end do
      end if

      !{Compute the covariance of the solution = ${\bf S}_s = ({\bf J}^T {\bf
      ! S}_m^{-1} {\bf J})^{-1}$, where ${\bf J}$ does not include
      ! Levenberg-Marquardt stabilization, and might or might not include
      ! Tikhonov regularization.  Start by computing the inverse transpose of
      ! the Cholesky factor of the normal equations, ${\bf U}^{-T}$, then
      ! $({\bf U}^T {\bf U})^{-1} = {\bf U}^{-1} {\bf U}^{-T} = ({\bf
      ! \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf J \Sigma}^T)^{-1}$. Then
      ! unscale it to ${\bf \Sigma}^T ({\bf \Sigma}^T {\bf J}^T {\bf
      ! S}_m^{-1} {\bf J \Sigma}^T)^{-1} {\bf \Sigma}^T$.
      if ( foundBetterState .and. ( got(f_outputCovariance) .or. got(f_outputSD) &
        & .or. got(f_average) .or. got(f_aprioriFraction) ) ) then
        if ( nwt_flag /= nf_getJ ) then
          print *, 'BUG in Retrieval module -- should need to getJ to quit'
          stop
        end if
          if ( d_cov ) then
            call time_now ( t3 )
            call output ( t3-t0, before='Begin covariance calculation at ', &
              & after=' seconds', advance='yes' )
            call output ( jacobian_rows, before='Counted ' )
            call output ( jacobian_cols, before=' rows and ', after=' columns' )
            call time_now ( t3 )
            call output ( t3-t0, before=' in the Jacobian matrix at ', &
              & after=' seconds', advance='yes' )
          end if
          call time_now ( t1 )
        if ( lambda /= 0 ) then
          ! Need to remove the Levenberg-Marquardt parameter from the diagonal
          ! of the normal equations.  Lambda is the last-used value.
          call updateDiagonal ( normalEquations, -lambda**2 )
          call choleskyFactor ( factored, normalEquations, status=matrixStatus )
            call add_to_retrieval_timing( 'cholesky_factor', t1 )
        end if
        ! We now have factored normal equations without the Levenberg-Marquardt
        ! parameter, and without Tiknonov regularization if CovSansReg was set.
        ! Invert the Cholesky factor to compute the a posteriori covariance
        ! matrix.
        call createEmptyMatrix ( tempU, 0, state, state )
        call invertCholesky ( factored, tempU ) ! U^{-1}
          if ( d_cov ) then
            call output ( &
              & 'Inverted the Cholesky factor of the normal equations at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
          end if
        ! Scale the covariance
        if ( columnScaling /= l_none ) then
          call rowScale ( v(columnScaleVector), tempU )
            if ( d_cov ) then
              call time_now ( t3 )
              call output ( t3-t0, before='Scaled the Covariance matrix at ', &
                & after=' seconds', advance='yes' )
            end if
        end if
        preserveMatrixName = outputCovariance%m%name
        ! TempU is the inverse of the Cholesky factor of the normal equations,
        ! which is now the Cholesky factor of the covariance matrix.  Compute
        ! the covariance matrix from its Cholesky factor.
        call multiplyMatrix_XY_T ( tempU, tempU, outputCovariance%m, &
          & diagonalOnly = .not. any ( got ( (/ f_outputCovariance, f_average/) ) ) ) ! U^{-1} U^{-T}
          if ( d_cov ) then
            call output ( 'Computed ' )
            if ( .not. any ( got ( (/ f_outputCovariance, f_average /) ) ) ) &
              & call output ( 'diagonal blocks of ' )
            call time_now ( t3 )
            call output ( t3-t0, before='U^{-1} U^{-T} at ', after=' seconds', &
              & advance='yes' )
          end if
        call add_to_retrieval_timing( 'cholesky_invert', t1 )
        outputCovariance%m%name = preserveMatrixName
        if ( associated(outputSD) ) &
          & call GetDiagonal ( outputCovariance%m, outputSD, squareRoot=.true. )

        ! Deal with possibly setting the error bar negative
        if ( negateSD .or. associated(aprioriFraction) ) then
          ! We need to compute what just the a priori and regularization would
          ! give for the output covariance.
          ! Firstly, we need to do some gymnastics to avoid singular matrices.  This
          ! will happen with ptan terms that have neither a priori nor smoothing.
          ! Set these such that our answer is the same as the outputSD
          call cloneVector ( v(diagFlagA), state )
          call GetDiagonal ( aPlusRegNEQ%m, v(diagFlagA) )
          do qty = 1, v(diagFlagA)%template%noQuantities
            where ( v(diagFlagA)%quantities(qty)%values <= 0.0_rv )
              ! Choose a value to let it invert properly
              v(diagFlagA)%quantities(qty)%values = outputSD%quantities(qty)%values ** (-2)
            elsewhere
              v(diagFlagA)%quantities(qty)%values = 0.0_rv
            end where
          end do
          call UpdateDiagonal ( aPlusRegNEQ, v(diagFlagA) )

          ! Now decompose it
          call createEmptyMatrix ( tempC%m, &
            & enter_terminal('_tempC', t_identifier), &
            & state, state, .not. aPlusRegNEQ%m%row%instFirst, .not. aPlusRegNEQ%m%col%instFirst )
          call CholeskyFactor ( tempC, aPlusRegNEQ, status=matrixStatus )
          if ( any ( matrixStatus /= 0 ) ) then
            call MLSMessage ( MLSMSG_Error, 'aPlusRegNEQ', &
              & 'Non SPD matrix to be factored: aPlusRegNEQ' )
          end if

          call MultiplyMatrix_XY ( tempC%m, tempU, temp ) ! Can't do C := C * U

          call MultiplyMatrix_XY_T ( temp, temp, aPlusRegCov%m, diagonalOnly=.true. )

          ! Destroy works in progress
          call DestroyMatrix ( tempC%m )
          call DestroyMatrix ( temp )

          call cloneVector ( v(diagFlagB), state )
          call GetDiagonal ( aPlusRegcov%m, v(diagFlagB), squareRoot=.true. )
          ! Remove ones where we substituted aposteriori precision.
          do qty = 1, v(diagFlagA)%template%noQuantities
            where ( v(diagFlagA)%quantities(qty)%values /= 0.0_rv ) &
              & v(diagFlagB)%quantities(qty)%values = 0.0_rv
          end do
          if ( associated(aprioriFraction) ) &
            call copyVector ( aprioriFraction, v(diagFlagB) )

          if ( negateSD ) then
            ! Go through and set error bar negative if appropriate
            do qty = 1, v(diagFlagA)%template%noQuantities
              where ( v(diagFlagB)%quantities(qty)%values > precisionFactor ) &
                & outputSD%quantities(qty)%values = - outputSD%quantities(qty)%values
            end do
          end if
        end if ! ( negateSD .or. associated(aprioriFraction) )
        call DestroyMatrix ( tempU )
      end if

      !{Compute the averaging kernel = ${\bf S}_s {\bf K}^T {\bf K}$.
      if ( foundBetterState .and. got(f_average) ) then
        if ( extendedAverage ) then
          call ClearMatrix ( outputAverage )
          call MultiplyMatrix_XTY ( outputCovariance%m, kTkStar, outputAverage, update=.true. )
        else
          ! Make sure kTk is symmetrical
          !  (outputCovariance is by virtue of its creation method as U^T U)
          Call ReflectMatrix ( kTk%m )
          call ClearMatrix ( outputAverage )
          call MultiplyMatrix_XTY ( outputCovariance%m, ktk%m, outputAverage, update=.true. )
        end if
        if ( d_cov ) call output ( &
          & 'Computed the Averaging Kernel from the Covariance', advance='yes' )
      end if

      call copyVector ( state, v(x) )
      if ( d_svec ) call dump ( state, name='Final state' )
      ! Clean up the temporaries, so we don't have a memory leak.
      if ( got(f_fuzz) ) call destroyVectorInfo ( fuzzState )
      call destroyMatrix ( normalEquations%m )
      call destroyMatrix ( kTkSep%m )
      call destroyMatrix ( kTkStar )
      call destroyMatrix ( factored%m )
      call destroyVectorInfo ( q )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      call add_to_retrieval_timing( 'newton_solver', t1 )
      call trace_end ( 'Retrieve.NewtonSolver', cond=.false. )
    end subroutine NewtonSolver

    ! ---------------------------------------  NewtonSolverFailed  -----
    subroutine NewtonSolverFailed ( WHY, NWT_Flag, Prev_NWT_Flag, &
      & D_NDB, AJ )
      use DNWT_Module, only: NF_START, NWT_T
      character(len=*), intent(in) :: WHY      ! What failed
      integer, intent(inout) :: NWT_Flag       ! Solver's state flag
      integer, intent(inout) :: Prev_NWT_Flag  ! Solver's previous state
      integer, intent(in) :: D_NDB             ! dump AJ if d_ndb >= 2
      type (NWT_T), intent(in) :: AJ           ! Solver's database

!       block ! for abandoned chunk version, followed by EXIT NEWT
!         chunk%abandoned = .true.
!         call MLSMessage ( MLSMSG_Warning, ModuleName, &
!           & 'Retrieval abandoned due to ' // trim(why) )
!       end block
!     block ! for forced gradient move version, followed by CYCLE NEWT
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Gradient move forced due to ' // trim(why) )
        if ( d_ndb >= 2 ) call my_nwtdb ( aj, width=9, why='Failed' )
        prev_nwt_flag = nwt_flag
        nwt_flag = nf_start ! Force gradient move
!     end block
    end subroutine NewtonSolverFailed

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for Retrieve = " )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine Retrieve

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module RetrievalModule

! $Log$
! Revision 2.368  2022/05/19 17:19:04  pwagner
! Perform NeuralNet even when skipping retrievals
!
! Revision 2.367  2021/02/05 05:18:50  pwagner
! Avoid calling NeuralNet when skipping Retrievals, too
!
! Revision 2.366  2021/01/22 00:21:10  pwagner
! Added NeuralNet command
!
! Revision 2.365  2020/07/17 19:38:52  vsnyder
! Check that FwdModelExtra exists if -Sfwmparallel is set
!
! Revision 2.364  2019/06/24 23:29:26  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.363  2018/11/20 01:09:43  vsnyder
! Always pass AJ to My_NWTDB.  It's no longer optional to DNWTDB
!
! Revision 2.362  2018/09/13 20:23:49  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.361  2018/08/28 20:50:10  vsnyder
! Don't print message about fnmin being negative if |fnmin| is small
!
! Revision 2.360  2018/07/27 23:19:53  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.359  2018/07/23 23:29:56  vsnyder
! Remove ChiSqMinNorm and ChiSqNorm from DNWT AJ structure because DNWT
! doesn't compute them or use them.
!
! Revision 2.358  2018/07/23 22:20:21  vsnyder
! Add AJ to all argument lists for nwt* so the options can be (eventually)
! stored there instead of as saved module variables.
!
! Revision 2.357  2017/12/07 01:01:23  vsnyder
! Don't use host-associated variable as a DO index
!
! Revision 2.356  2015/05/05 00:12:47  vsnyder
! Make sure 0 < mu <= 1 in BoundMove
!
! Revision 2.355  2014/09/30 02:15:19  vsnyder
! Don't put the Tikhonov matrix in the matrix database.  Add some commented-
! out stuff that might be useful if we turn on matrix finalizers.
!
! Revision 2.354  2014/09/29 20:17:52  vsnyder
! Fuse two IF constructs with identical conditions, some cannonball polishing
!
! Revision 2.353  2014/09/05 01:20:03  vsnyder
! Remove unused State argument from DumpStateQuantities.  Remove USE for
! unreferenced USE name.
!
! Revision 2.352  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.351  2014/04/10 00:45:24  pwagner
! Moved currentChunkNumber, currentPhaseName from MLSL2Timings to MLSL2Options
!
! Revision 2.350  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.349  2014/02/28 01:14:04  vsnyder
! Remove TYPE argument from calls to EXPR because the value wasn't used.
! Move units checking to type checker.
!
! Revision 2.348  2014/01/23 23:16:37  vsnyder
! Back out Levenberg-Marquardt without getting a new Jacobiab
!
! Revision 2.347  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.346  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.345  2013/10/09 23:42:29  vsnyder
! Add Evaluate_Variable
!
! Revision 2.344  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.343  2013/09/04 02:50:05  vsnyder
! Add 'cond' argument in three calls to Trace_End
!
! Revision 2.342  2013/08/30 02:45:46  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.341  2013/08/03 00:40:42  vsnyder
! Require vector.quantity in dumpQuantities field in retrieve spec
!
! Revision 2.340  2013/07/26 22:41:44  vsnyder
! Fiddling with some dump switches
!
! Revision 2.339  2013/07/12 23:25:28  vsnyder
! Remove unreferenced error messages
!
! Revision 2.338  2013/06/12 02:38:50  vsnyder
! Cruft removal
!
! Revision 2.337  2013/04/24 00:37:42  pwagner
! Added InitRepeat and NextRepeat calls to set/increment r/t Boolean count
!
! Revision 2.336  2013/04/22 17:51:57  pwagner
! Reevaluate may store a literal instead of a Boolean value
!
! Revision 2.335  2013/04/17 00:05:25  pwagner
! Added new Repeat control structure to Fill, Retrieve sections
!
! Revision 2.334  2013/03/01 01:12:15  pwagner
! 'fiw' switch dumps sids instance window even if retrieval skipped
!
! Revision 2.333  2012/10/11 22:03:46  pwagner
! Fix timing error; print chunkNumber, phaseName instead of module when warning of abandoning
!
! Revision 2.332  2012/09/13 18:08:34  vsnyder
! Don't include tikhonov row count if covSansReg is set
!
! Revision 2.331  2012/09/12 22:51:40  vsnyder
! Tell DNWT what lambdaMin is
!
! Revision 2.330  2012/08/30 23:01:18  vsnyder
! Procedurize Tikhonov, add lambdaMin
!
! Revision 2.329  2012/08/16 18:07:23  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.328  2012/07/26 02:07:11  vsnyder
! Remove unused cruft
!
! Revision 2.327  2012/07/04 02:13:03  vsnyder
! Add dump of residuals.  Add dumpQuantities node to select quantities of
! the state vector to dump during the Newton iteration.  Remove masked rows
! and columns from the counts, to compute DOF correctly.
!
! Revision 2.326  2012/06/06 20:37:56  vsnyder
! Add toggles field to retrieve spec
!
! Revision 2.325  2012/05/08 17:50:52  pwagner
! Added Select .. Case .. EndSelect control structure
!
! Revision 2.324  2012/04/20 01:37:47  vsnyder
! Dump vector norms using dxn, copy mask from x to candidateDx, print DOF
! dump QNSQ if sca, some cannonball polishing
!
! Revision 2.323  2012/03/13 01:48:25  vsnyder
! Add drmc flag to control clean matrix printing
!
! Revision 2.322  2012/03/07 01:58:40  vsnyder
! Respect mask in some dot products where it should have done all along
!
! Revision 2.321  2012/02/10 23:46:24  vsnyder
! Add dump for kTk
!
! Revision 2.320  2012/02/01 00:18:38  vsnyder
! Remove Matrix from the Retrieve section; nobody used it, because it was broken.
!
! Revision 2.319  2012/02/01 00:17:07  vsnyder
! init_tables_module.f90
!
! Revision 2.318  2012/01/27 01:06:26  pwagner
! Tried to fix memory leak created by Hessian
!
! Revision 2.317  2012/01/04 02:14:55  vsnyder
! Ensure forward model does not see MIFExtinction
!
! Revision 2.316  2012/01/04 01:51:26  vsnyder
! Create fwmJacobian if needed and one isn't supplied
!
! Revision 2.315  2011/12/21 01:42:22  vsnyder
! Add MIFExtinction transformation
!
! Revision 2.314  2011/08/29 22:13:42  pwagner
! Predefine Hessian object in parallel with Jacobian matrix
!
! Revision 2.313  2011/05/09 18:24:31  pwagner
! Converted to using switchDetail
!
! Revision 2.312  2010/08/27 23:59:36  vsnyder
! Removed some overly long comments that duplicate the CVS log anyway
!
! Revision 2.311  2010/08/27 06:31:45  yanovsky
! Added   type(hessian_T), pointer :: Hessian  in Retrieve subroutine.
! Hessian is an actual argument in a call to forwardModel.
!
! Revision 2.309  2010/03/24 20:56:46  vsnyder
! Add Hessian database to DumpBlocks call
!
! Revision 2.308  2010/02/25 18:19:27  pwagner
! Adds support for new Hessian database
!
! Revision 2.307  2009/11/23 21:10:18  vsnyder
! Gradiant move instead of abandonment when trouble occurs
!
! Revision 2.306  2009/10/26 17:12:39  pwagner
! Added Diff command to be used like Dump in l2cf
!
! Revision 2.305  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.304  2009/06/16 17:40:02  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.303  2009/03/14 02:43:19  honghanh
! Add dnwt_abandoned and dnwt_count
!
! Revision 2.302  2008/12/18 21:13:05  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.301  2007/12/07 01:14:03  pwagner
! Lets us catch warnings and assign to runtime Booleans
!
! Revision 2.300  2007/11/15 22:53:16  pwagner
! May set runtimeBooleans by anyGood.., Compare, Reevaluate commands
!
! Revision 2.299  2007/11/08 03:24:39  vsnyder
! Add stateMax and stateMin
!
! Revision 2.298  2007/11/05 18:37:19  pwagner
! May Skip remaining lines in Fill, Join, Retrieve sections depending on Boolean
!
! Revision 2.297  2007/10/04 20:43:12  vsnyder
! Remove unused symbols
!
! Revision 2.296  2007/10/04 01:49:42  vsnyder
! Correct a string's misspelling
!
! Revision 2.295  2007/10/02 22:51:27  vsnyder
! Put switches field value into global switches variable so the models can
! see them, then clear them at the end so they don't stick around.
!
! Revision 2.294  2007/09/12 00:17:05  vsnyder
! Simplify printing when Cholesky fails
!
! Revision 2.293  2007/09/08 00:31:49  vsnyder
! Be more careful before printing name of failed Cholesky block
!
! Revision 2.292  2007/09/06 00:36:51  vsnyder
! More information about failed Cholesky
!
! Revision 2.291  2007/08/20 22:05:30  pwagner
! Many procedures now push their names onto MLSCallStack
!
! Revision 2.290  2007/06/29 19:32:07  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.289  2007/04/03 17:47:14  vsnyder
! Replace pointer attribute on VectorDatabase with target attribute
!
! Revision 2.288  2007/02/28 23:27:51  vsnyder
! Replace FORALL by DO, since FORALL sometimes causes compilers to have trouble
!
! Revision 2.287  2006/12/09 00:45:35  vsnyder
! Move calculation of size of Jacobian, so that ChiSqMin calculation works
! correctly if zero Newton iterations are allowed.
!
! Revision 2.286  2006/10/04 22:52:48  vsnyder
! Another change from Herb:  Make Apriori fraction zero where we had to
! substitute aposteriori precision to calculate the apriori fraction.
!
! Revision 2.285  2006/10/04 19:39:38  vsnyder
! Correct negateSD as recommended by Herb
!
! Revision 2.284  2006/09/21 18:51:14  pwagner
! Reduce level of dumps in SIDS version
!
! Revision 2.283  2006/09/20 00:43:38  vsnyder
! Implemented Herb's apriori fraction calculation
!
! Revision 2.282  2006/08/11 20:58:38  vsnyder
! Add 'simple' method to use alternate Newton solver
!
! Revision 2.281  2006/08/04 18:11:33  vsnyder
! Add LeakCheck command
!
! Revision 2.280  2006/08/01 02:48:33  vsnyder
! Remove unused USE for .TX.
!
! Revision 2.279  2006/07/21 20:13:29  pwagner
! Can fill state even if skipping retrievals
!
! Revision 2.278  2006/07/19 22:28:59  vsnyder
! Cannonball polishing
!
! Revision 2.277  2006/06/08 23:59:46  vsnyder
! Simplify column scaling, add switches field, TeXnicalities
!
! Revision 2.276  2006/06/06 18:54:21  vsnyder
! Sharpen criteria for fresh Jacobian at the end
!
! Revision 2.275  2006/06/06 00:51:23  vsnyder
! Move diagnostic output from NEWX to DX, DX_Aitken, and GMOVE
!
! Revision 2.274  2006/06/06 00:31:20  vsnyder
! Add more diagnostic output
!
! Revision 2.273  2006/06/03 01:04:59  vsnyder
! Cannonball polishing
!
! Revision 2.272  2006/06/03 00:15:20  vsnyder
! Correct updating of Levenberg-Marquardt parameter during More' and
! Sorensen iteration.  Polish up some dumps.
!
! Revision 2.271  2006/05/31 18:22:26  pwagner
! Fixed bug only NAG complained about: needed Dump from MatrixModule_0
!
! Revision 2.270  2006/05/30 22:51:21  vsnyder
! Precompute dump flags.  Add NumGrad and NumNewt counters.  Compute
! cosine of angles between moves correctly for printing (it's already
! done correctly inside dnwt_module).
!
! Revision 2.269  2006/03/30 19:00:50  vsnyder
! Correct CVS log comment
!
! Revision 2.268  2006/03/30 18:57:01  vsnyder
! Correct some TeXnicalities
!
! Revision 2.267  2006/03/22 23:47:14  vsnyder
! Decruftification
!
! Revision 2.266  2006/02/10 21:17:09  pwagner
! dumps may go to special dumpfile
!
! Revision 2.265  2006/01/21 00:04:14  livesey
! Added the start of a dump retrieval config option (-Srtv)
!
! Revision 2.264  2005/12/21 21:48:03  livesey
! Added handling of the negateSD option.
!
! Revision 2.263  2005/06/03 02:09:46  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.262  2005/05/27 23:57:03  vsnyder
! Add Flush PFAData
!
! Revision 2.261  2005/03/15 01:27:58  vsnyder
! Allow Dump command in Retrieve section
!
! Revision 2.260  2004/06/16 23:41:39  vsnyder
! Put a limit on the total number of DNWT iterations, regardless of Jacobians
!
! Revision 2.259  2004/06/16 01:21:57  vsnyder
! Repair bug in constrained move
!
! Revision 2.258  2004/06/16 00:01:36  livesey
! Bug fixes to BoundMove
!
! Revision 2.257  2004/05/19 19:16:12  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.256  2004/04/19 23:25:55  livesey
! Added an only to the use of ieee_arithemtic
!
! Revision 2.255  2004/04/19 19:23:44  livesey
! Better handling of NaNs in radiances or derivatives
!
! Revision 2.254  2004/01/29 01:45:32  livesey
! Removed print statement
!
! Revision 2.253  2004/01/24 03:22:49  livesey
! Changed logistics of computing averaging kernel to preserve vectors that
! describe edges.
!
! Revision 2.252  2003/10/07 01:17:52  vsnyder
! Change spelling of dumpBlock to dumpBlocks
!
! Revision 2.251  2003/09/11 23:16:11  livesey
! Now hands vector database onto forward model to support xStar/yStar in
! linearized forward model.
!
! Revision 2.250  2003/09/05 23:24:09  pwagner
! Skips calls to Retrieve, SIDS if SKIPRETRIEVAL
!
! Revision 2.249  2003/07/09 23:49:13  vsnyder
! Remove numerous unreferenced USE names
!
! Revision 2.248  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.247  2003/06/03 19:24:19  livesey
! Added the abandoned stuff to make the chunk processing more robust.
!
! Revision 2.246  2003/05/14 23:40:51  dwu
! a change in cloud retr
!
! Revision 2.245  2003/05/14 23:38:15  dwu
! change in cloud retrieval
!
! Revision 2.244  2003/05/14 03:55:10  dwu
! tidy up
!
! Revision 2.243  2003/05/13 22:25:39  dwu
! changes in lowcloudretrieval
!
! Revision 2.242  2003/05/13 20:43:05  dwu
! a quick fix after spinout
!
! Revision 2.241  2003/05/13 19:32:28  dwu
! spin out cloud retrievals
!
! Revision 2.240  2003/04/04 22:02:06  livesey
! Added call to UpdateMask
!
! Revision 2.239  2003/03/27 20:45:22  livesey
! Reinstated the fnorm dumping which seemed to fall by the wayside
! for some reason.
!
! Revision 2.238  2003/03/07 03:17:03  livesey
! Added Restrict Range
!
! Revision 2.237  2003/03/06 00:46:41  livesey
! Now gets subset etc. from subsetmodule
!
! Revision 2.236  2003/02/25 21:12:48  dwu
!  clean up FlagCloud
!
! Revision 2.235  2003/02/25 19:04:43  dwu
! fix another bug in FlagCloud
!
! Revision 2.234  2003/02/25 18:58:35  dwu
! fix another bug in FlagCloud
!
! Revision 2.233  2003/02/24 19:20:40  dwu
! fix a bug in FlagCloud
!
! Revision 2.232  2003/02/14 01:56:36  livesey
! Added the 'additional' capability in subset
!
! Revision 2.231  2003/02/12 02:11:13  livesey
! New code for extended averaging kernels
!
! Revision 2.230  2003/02/08 00:20:25  livesey
! Added error message to check that the right ptan (GHz/THz) is used in
! subset statements.
!
! Revision 2.229  2003/02/05 20:59:16  dwu
! an improvement in flagcloud
!
! Revision 2.228  2003/02/05 04:07:05  dwu
! add cloudheight option for flagcloud
!
! Revision 2.227  2003/02/03 23:08:50  vsnyder
! Delete test for small residual -- there's a test relative to FNMIN in
! dnwt.  Add stuff for dnwt to use Mor\'e and Sorensen algorithm to
! compute the Levenberg parameter.
!
! Revision 2.225  2003/01/18 02:15:43  vsnyder
! Change names of global loop inductors "I" and "J" to something more clever,
! i.e., I_Sons and I_Key.  It's too easy to use "I" and "J" in internal
! procedures, thereby clobbering the important values of these loop
! inductors.
!
! Revision 2.224  2003/01/18 01:40:10  vsnyder
! Prepare for More and Sorensen
!
! Revision 2.223  2003/01/17 22:59:47  livesey
! Minor bug fix in units checking for max/minValue in subset.
!
! Revision 2.222  2003/01/17 22:13:46  dwu
! fix a bug in flagCloud
!
! Revision 2.221  2003/01/17 16:55:19  dwu
! a minor change in FlagCloud
!
! Revision 2.220  2003/01/16 21:48:22  vsnyder
! More stuff on getting NWTA internal output
!
! Revision 2.219  2003/01/16 04:06:58  vsnyder
! Change some output in DNWT
!
! Revision 2.218  2003/01/15 01:49:29  vsnyder
! Revise dumping stuff for Newton method
!
! Revision 2.217  2003/01/14 22:14:43  dwu
! make FlagCloud depend on both channels and cloudChannels
!
! Revision 2.216  2003/01/13 17:17:29  jonathan
!  change cloud_width to i_saturation
!
! Revision 2.215  2003/01/12 07:33:42  dwu
! with some fix in flagcloud
!
! Revision 2.214  2003/01/11 15:51:34  dwu
! follow-up fix for FlagCloud
!
! Revision 2.213  2003/01/11 08:08:01  dwu
! add flagCloud
!
! Revision 2.212  2003/01/11 02:14:55  livesey
! Bug fix in max/min value subset
!
! Revision 2.211  2003/01/11 01:02:50  livesey
! Added max and min value to subset
!
! Revision 2.210  2003/01/08 23:52:38  livesey
! Added the sparsify stuff according to sparseQuantities
!
! Revision 2.209  2003/01/07 23:44:17  livesey
! Added reset option to subset
!
! Revision 2.208  2002/12/11 01:58:54  livesey
! Changed detail of forward model parallel stuff
!
! Revision 2.207  2002/12/06 18:41:24  livesey
! Slightly new approach to FWMParallel mode
!
! Revision 2.206  2002/11/23 00:07:13  vsnyder
! Delete some unused symbols
!
! Revision 2.205  2002/11/22 00:11:16  livesey
! Better handling of foundBetterState
!
! Revision 2.204  2002/11/20 21:05:55  livesey
! More informative dump for Tikhonov
!
! Revision 2.203  2002/11/20 01:10:03  livesey
! Added the foundBetterState stuff
!
! Revision 2.202  2002/10/29 20:51:27  livesey
! Moved call to FillDiagVec
!
! Revision 2.201  2002/10/25 23:56:04  livesey
! Bug fix in the diagnostics
!
! Revision 2.200  2002/10/25 22:24:42  livesey
! Changed the diagnostic vector handling.
!
! Revision 2.199  2002/10/25 01:13:46  livesey
! Jacobian rows/cols etc. now actually describe the jacobian.
!
! Revision 2.198  2002/10/23 23:27:40  livesey
! Bug fix in jacobian_rows calculation
!
! Revision 2.197  2002/10/23 01:32:31  vsnyder
! Add CovSansReg switch to retrieve spec
!
! Revision 2.196  2002/10/23 01:14:41  livesey
! Added the chiSquared stuff, needs more work though.
!
! Revision 2.195  2002/10/19 23:41:04  livesey
! Added muMin functionality
!
! Revision 2.194  2002/10/19 18:49:17  livesey
! Various bug fixes in bounds stuff
!
! Revision 2.193  2002/10/19 01:52:21  livesey
! Added serial option and associated flags.
!
! Revision 2.192  2002/10/18 22:05:42  vsnyder
! Get in bounds before starting.  Account for scaling when bounding a move.
!
! Revision 2.191  2002/10/17 23:12:51  vsnyder
! Apply bounds to gradient and Aitken moves
!
! Revision 2.190  2002/10/17 00:16:40  vsnyder
! Add lowBound and highBound fields for the Retrieve spec
!
! Revision 2.189  2002/10/08 17:41:38  livesey
! Bug fixes in FWMParallel stuff
!
! Revision 2.188  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.187  2002/10/06 02:04:57  livesey
! Put in fwmParallel functionality
!
! Revision 2.186  2002/09/25 20:08:43  livesey
! Made subset -g less verbose
!
! Revision 2.185  2002/09/23 23:15:08  vsnyder
! Delete maxF
!
! Revision 2.184  2002/09/23 20:21:21  livesey
! Now copies mask from f_rowScaled to f.  Only important for snooping.
!
! Revision 2.183  2002/09/21 00:33:34  vsnyder
! Don't take a getJ iteration if TOLF convergence occurs
!
! Revision 2.182  2002/09/21 00:04:14  vsnyder
! Put back an erroneously-removed timing call
!
! Revision 2.181  2002/09/19 01:26:46  vsnyder
! Mostly fixing up comments and LaTeX stuff
!
! Revision 2.180  2002/09/18 23:56:01  vsnyder
! Call time_now at end of add_to_retrieval_timing
!
! Revision 2.179  2002/09/14 02:46:30  vsnyder
! Add NDB switch to get DNWT output every return
!
! Revision 2.178  2002/09/14 00:37:37  vsnyder
! More stuff on subsuming EVALF into EVALJ properly
!
! Revision 2.177  2002/09/13 23:55:01  livesey
! Changed print statements
!
! Revision 2.176  2002/09/13 20:05:06  vsnyder
! Subsume EVALF into EVALJ
!
! Revision 2.175  2002/09/13 18:10:10  pwagner
! May change matrix precision rm from r8
!
! Revision 2.174  2002/09/11 20:21:49  livesey
! Made second | F | print on GETJ not EVALJ
!
! Revision 2.173  2002/09/11 17:40:59  livesey
! Changed meaning of v(f)
!
! Revision 2.172  2002/09/11 14:06:42  livesey
! Bug fix, misunderstood meaning of v(f)
!
! Revision 2.171  2002/09/11 01:14:47  livesey
! Added extra dump of | F |
!
! Revision 2.170  2002/09/06 00:46:18  livesey
! Added dump of aTb
!
! Revision 2.169  2002/09/05 23:08:22  livesey
! Gradient move handles mask 'correctly', though we might have been right
! from a philosophical point of view earlier.
!
! Revision 2.168  2002/08/29 04:45:37  livesey
! Sped up covariance calculations
!
! Revision 2.167  2002/08/28 00:51:04  vsnyder
! Correct more blunders in Tikhonov regularization
!
! Revision 2.166  2002/08/26 20:01:59  livesey
! Made subset us the GetIndexFlagsFromList routine to get channel flags
!
! Revision 2.165  2002/08/24 01:38:28  vsnyder
! Implement horizontal regularization
!
! Revision 2.164  2002/08/22 00:06:50  vsnyder
! Implement regularize-to-apriori
!
! Revision 2.163  2002/08/21 19:55:31  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.162  2002/08/20 19:55:02  vsnyder
! Use matrix with same row and column templates for Tikhonov
!
! Revision 2.161  2002/08/08 22:02:52  vsnyder
! Add hRegWeightVec and Tikhonov mask
!
! Revision 2.160  2002/08/06 02:16:09  livesey
! Fixed averaging kernels by reflecting the kTk matrix
!
! Revision 2.159  2002/08/05 19:40:11  vsnyder
! Undo scaling of kTk -- it's prepared unscaled
!
! Revision 2.158  2002/08/03 20:40:58  livesey
! Added matrix name preseveration.  Changed averaging kernel calculation.  Still
! not right though!
!
! Revision 2.157  2002/08/03 01:16:17  vsnyder
! RegAfter switch controls Tikhonov before/after column scaling -- default before
!
! Revision 2.156  2002/07/31 22:43:30  vsnyder
! Add some timers in the Newtonian iteration
!
! Revision 2.155  2002/07/26 22:47:56  vsnyder
! Only do DNWT internal output after NF_SOLVE and at the end.
! Finish commenting-out FNMIN calculation when NWT_FLAG = NF_SOLVE
!
! Revision 2.154  2002/07/26 01:20:22  vsnyder
! Exploit added output from DNWTDB, get rid of bogus Levenberg at start
!
! Revision 2.153  2002/07/24 01:08:40  vsnyder
! Don't compute FNMIN at NF_SOLVE time -- it should be the same as at NF_EVALJ
!
! Revision 2.152  2002/07/22 23:14:28  pwagner
! Fixed some bugs in timings (there may be more)
!
! Revision 2.151  2002/07/22 22:53:42  pwagner
! Added form norm eq and tikh reg timings
!
! Revision 2.150  2002/07/08 21:05:09  vsnyder
! Make sure "Best X" has a value even if no Newtonian iterations are done.
! Mark Filipiak noticed this problem.
!
! Revision 2.149  2002/07/04 00:34:55  vsnyder
! Make aprioriMinusX apriori - x instead of -x
!
! Revision 2.148  2002/07/02 01:38:06  vsnyder
! Quit if a a retrieve can't be done due to L2CF errors
!
! Revision 2.147  2002/07/01 23:43:17  vsnyder
! Plug some memory leaks
!
! Revision 2.146  2002/06/25 20:45:31  vsnyder
! DestroyVectorValue should have been ClearVector in the case that there is
! no apriori covariance.
!
! Revision 2.145  2002/06/18 01:19:52  vsnyder
! Cosmetic changes
!
! Revision 2.144  2002/06/07 01:34:26  vsnyder
! Create temp matrix used in covariance calculation
!
! Revision 2.143  2002/05/22 19:16:51  vsnyder
! Put the correct number of rows used for Tikhonov regularization into the
! diagnostics vector.
!
! Revision 2.142  2002/05/22 19:00:38  vsnyder
! Correct covariance calculation -- it ought to be U^{-1} U^{-T}, not U^{-1}.
! Mark Filipiak noticed this bug.
!
! Revision 2.141  2002/05/07 01:02:24  vsnyder
! Change regWeight to hRegWeights -- which is now a tree node instead of
! a real scalar.  Add dump for regularization matrix.
!
! Revision 2.140  2002/04/22 23:00:58  vsnyder
! Add SquareRoot=.true. to getDiagonal for standard deviation
!
! Revision 2.139  2002/04/22 20:55:00  vsnyder
! Compute and output the averaging kernel
!
! Revision 2.138  2002/03/13 22:01:50  livesey
! Changed explicitFill to fill
!
! Revision 2.137  2002/03/08 08:07:16  livesey
! Added explicit fill mask
!
! Revision 2.136  2002/03/06 01:44:30  livesey
! Changed manner in which optical depth cutoff applied so that
! once a channel is too optically thick it always is (as a function of
! height).  Note it's resistant to changes in scan direction.
!
! Revision 2.135  2002/02/14 21:55:31  vsnyder
! Account for mask.  Get a final Jacobian using Best X.  Cosmetic changes.
!
! Revision 2.134  2002/02/13 00:11:13  vsnyder
! Test fnmin both places it's formed
!
! Revision 2.133  2002/02/12 22:53:21  vsnyder
! Update a bunch of comments
!
! Revision 2.132  2002/02/09 19:11:34  livesey
! Modified subset to add optical depth cutoffs
!
! Revision 2.131  2002/02/08 22:51:59  livesey
! Added call to CopyVectorMask to transfer mask from measurements to forward
! model.
!
! Revision 2.130  2002/02/07 02:55:02  vsnyder
! Add a 'mask' field to the 'setup' spec
!
! Revision 2.129  2002/02/05 02:40:52  vsnyder
! Use 'dumpMask' instead of 'dump' to dump the mask, cosmetic changes
!
! Revision 2.128  2002/01/18 00:31:23  livesey
! Changed error message about fnmin being imaginary to warning and put
! in a work around.  I want to see what it's doing.
!
! Revision 2.127  2001/12/03 18:50:43  pwagner
! NAG compiler mad at '/ =' instead of '/='
!
! Revision 2.126  2001/12/01 01:01:56  livesey
! Added vir switch
!
! Revision 2.125  2001/11/28 23:16:39  livesey
! Added use of MLSMSG_Warning, whoops!
!
! Revision 2.124  2001/11/28 20:39:35  livesey
! Fixed bug with subset, changed some of Dong's prints to MLSMessage's
!
! Revision 2.123  2001/11/28 18:28:50  livesey
! Updated subset to include open ranges in height specification.
! There is possibly a bug in Van's version I'll point out to him
! (if not my version has the opposite bug).
!
! Revision 2.122  2001/11/28 00:01:27  jonathan
! remove type (VectorValue_T), pointer :: Re
!
! Revision 2.121  2001/11/27 23:49:48  jonathan
! remove Re in HighCloudRetrieval
!
! Revision 2.120  2001/11/27 23:34:49  pwagner
! Split forward model timings into four types
!
! Revision 2.119  2001/11/27 01:28:17  vsnyder
! Implement (partially) open range for channels in subset
!
! Revision 2.118  2001/11/17 02:30:44  vsnyder
! Add L_numF and L_numJ to 'diagnostics' vector
!
! Revision 2.117  2001/11/13 23:35:12  vsnyder
! Keep f and f_rowScaled separate in EVALF case
!
! Revision 2.116  2001/11/12 18:25:02  dwu
! speed up some cloud calculations
!
! Revision 2.115  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.114  2001/11/08 01:21:45  vsnyder
! Make sure phaseName has a value
!
! Revision 2.113  2001/11/07 23:55:43  dwu
! added option to constain cloud top height
!
! Revision 2.112  2001/11/06 18:04:11  dwu
! first working version of HighCloudRetrieval
!
! Revision 2.111  2001/11/02 01:11:45  dwu
! a test for highcloudretrieval
!
! Revision 2.110  2001/11/01 19:29:01  dwu
! corrections in cloud retrievals
!
! Revision 2.109  2001/10/31 21:59:56  livesey
! Added phaseName to snooper
!
! Revision 2.108  2001/10/30 16:57:43  dwu
! fixed a bug in LowCloudREtrieval
!
! Revision 2.107  2001/10/26 20:28:42  dwu
! Added cloud extinction SD and high cloud retrieval
!
! Revision 2.106  2001/10/24 00:29:00  livesey
! Added Matrix snooping
!
! Revision 2.105  2001/10/23 20:08:59  vsnyder
! Scale the RHS of apriori equations by aprioriScale too!
! Cosmetic changes -- LaTeX comments about how Tikhonov and apriori work.
!
! Revision 2.104  2001/10/22 22:40:21  livesey
! Tried to add a matrixDatabase to snoop call
!
! Revision 2.103  2001/10/22 22:34:47  dwu
! oops
!
! Revision 2.102  2001/10/22 22:25:17  dwu
! prepare for high cloud retrieval
!
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
