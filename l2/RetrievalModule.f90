! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module RetrievalModule
!=============================================================================

!{This module inverts the radiative transfer equation, to solve for
! atmospheric parameters, given radiance measurements.
!
! This module and ones it calls consume most of the cycles.

  implicit NONE
  private
  public :: RETRIEVE

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

  !{Process the ``Retrieve'' section of the L2 Configuration File.
  ! The ``Retrieve'' section can have ForwardModel, Matrix, Sids, Subset or
  ! Retrieve specifications.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use BitStuff, only: CountBits
    use Expr_M, only: Expr
    use ForwardModelConfig, only: ForwardModelConfig_T
    use Init_Tables_Module, only: F_additional, F_apriori, F_aprioriScale, F_Average, &
      & F_channels, F_cloudChannels, F_cloudHeight, F_cloudRadiance, &
      & F_cloudRadianceCutOff, &
      & F_columnScale, F_Comment, F_covariance, F_covSansReg, &
      & F_diagnostics, F_diagonal, F_extendedAverage, &
      & F_forwardModel, F_fuzz, F_fwdModelExtra, F_fwdModelOut, &
      & F_height, f_highBound, f_hRegOrders, f_hRegQuants, f_hRegWeights, &
      & f_hRegWeightVec, F_ignore, F_jacobian, F_lambda, F_Level, f_lowBound, &
      & F_mask, F_maxJ, F_maxValue, F_measurements, &
      & F_measurementSD, F_method, F_minValue, F_muMin, &
      & F_opticalDepth, F_opticalDepthCutoff, F_outputCovariance, F_outputSD, &
      & F_phaseName, F_ptanQuantity, F_quantity, &
      & F_regAfter, F_regApriori, F_reset, F_serial, F_SparseQuantities, &
      & F_state, F_toleranceA, F_toleranceF, &
      & F_toleranceR, f_vRegOrders, f_vRegQuants, &
      & f_vRegWeights, f_vRegWeightVec, Field_first, Field_last, &
      & L_apriori, L_covariance, &
      & L_dnwt_ajn,  L_dnwt_axmax,  L_dnwt_cait, L_dnwt_chiSqMinNorm, L_dnwt_chiSqNorm, &
      & L_dnwt_diag,  L_dnwt_dxdx, &
      & L_dnwt_dxdxl, L_dnwt_dxn,  L_dnwt_dxnl,  L_dnwt_flag, L_dnwt_fnmin, &
      & L_dnwt_fnorm,  L_dnwt_gdx,  L_dnwt_gfac, L_dnwt_gradn,  L_dnwt_sq, &
      & L_dnwt_sq,  L_dnwt_sqt, L_Fill, L_full_derivatives, &
      & L_highcloud, L_Jacobian_Cols, L_Jacobian_Rows, &
      & L_linalg, L_lowcloud, L_newtonian, L_none, L_norm, &
      & L_numJ, L_opticalDepth, L_pressure, L_radiance, L_Tikhonov, L_zeta, &
      & S_dumpBlocks, S_flagCloud, S_matrix, S_retrieve, S_sids, S_snoop, &
      & S_subset, S_time, S_RESTRICTRANGE, S_UPDATEMASK
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Invalid
    use L2ParInfo, only: PARALLEL
    use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
      & DestroyMatrix, GetFromMatrixDatabase, Matrix_T, Matrix_Database_T, &
      & Matrix_SPD_T, MultiplyMatrixVectorNoT, operator(.TX.), ReflectMatrix, &
      & Sparsify, MultiplyMatrix_XTY
    use MatrixTools, only: DumpBlock
    use MLSCommon, only: MLSCHUNK_T, R8, RM, RV
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, Add_To_Retrieval_Timing
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
    use Output_m, only: BLANKS, OUTPUT
    use SidsModule, only: SIDS
    use SnoopMLSL2, only: SNOOP
    use String_Table, only: Display_String, Get_String
    use SubsetModule, only: SETUPSUBSET, SETUPFLAGCLOUD, RESTRICTRANGE, UPDATEMASK
    use Time_M, only: Time_Now
    use Toggles, only: Gen, Switches, Toggle, Levels
    use Trace_M, only: Trace_begin, Trace_end
    use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
      & Subtree
    use Tree_Types, only: N_colon_less, N_less_colon, &
      & N_less_colon_less, N_named
    use VectorsModule, only: ClearMask, ClearUnderMask, &
      & ClearVector, CloneVector, CopyVector, CopyVectorMask, CreateMask, &
      & DestroyVectorInfo, DumpMask, GetVectorQuantityByType, &
      & IsVectorQtyMasked, M_Fill, M_FullDerivatives, M_LinAlg, &
      & M_Tikhonov, Vector_T, VectorValue_T
    use CloudRetrievalModule, only: CloudRetrieval

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST;
                                        ! It indexes an n_cf vertex
    type(vector_T), dimension(:), pointer :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: ConfigDatabase

    type(MLSChunk_T), intent(in) :: CHUNK

    ! Default values:
    real(r8), parameter :: DefaultInitLambda = 0.0_r8
    integer, parameter :: DefaultMaxJ = 5
    integer, parameter :: DefaultMethod = l_newtonian
    real(r8), parameter :: DefaultMuMin = 0.1_rv
    double precision, parameter :: DefaultToleranceA = 1.0d-6 ! for NWT
    double precision, parameter :: DefaultToleranceF = 1.0d-6 ! for NWT
    double precision, parameter :: DefaultToleranceR = 1.0d-6 ! for NWT

    ! Local variables:
    type(vector_T), pointer :: Apriori  ! A priori estimate of state
    real(r8) :: AprioriScale            ! Weight for apriori, default 1.0
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
    integer :: Error
    logical :: ExtendedAverage          ! Averaging kernels based on full
                                        ! Jacobian (last term) not masked
    integer :: Field                    ! Field index -- f_something
    real(r8) :: Fuzz                    ! For testing only.  Amount of "fuzz"
                                        ! to add to state vector before
                                        ! starting a retrieval.
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    logical :: Got(field_first:field_last)   ! "Got this field already"
    type(vector_T), pointer :: HighBound ! For state during retrieval
    integer :: HRegOrders               ! Regularization orders               
    integer :: HRegQuants               ! Regularization quantities           
    integer :: HRegWeights              ! Weight of regularization conditions 
    type(vector_T), pointer :: HRegWeightVec  ! Weight vector for regularization
    integer :: I_Key, I_Sons, J         ! Subscripts and loop inductors
    real(r8) :: InitLambda              ! Initial Levenberg-Marquardt parameter
    integer :: IxAverage                ! Index in tree of averagingKernel
    integer :: IxCovariance             ! Index in tree of outputCovariance
    integer :: IxJacobian               ! Index in tree of Jacobian matrix
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Jacobian_Cols            ! Number of columns of the Jacobian.
    integer :: Jacobian_Rows            ! (Number of rows of Jacobian) -
                                        ! (masked-off rows of Jacobian)
    integer :: K                        ! Subscript, loop inductor, local temp
    integer :: Key                      ! Index of an n_spec_args.  Either
                                        ! a son or grandson of root.
    type(vector_T), pointer :: LowBound ! For state during retrieval
    integer :: MaxJacobians             ! Maximum number of Jacobian
                                        ! evaluations of Newtonian method
    type(vector_T), pointer :: Measurements  ! The measurements vector
    type(vector_T), pointer :: MeasurementSD ! The measurements vector's Std. Dev.
    integer :: Method                   ! Method to use for inversion, currently
                                        ! only l_Newtonian.
    type(matrix_T), target :: MyAverage ! for OutputAverage to point to
    type(matrix_SPD_T), target :: MyCovariance    ! for OutputCovariance to point at
    type(matrix_T), target :: MyJacobian          ! for Jacobian to point at
    real(rv) :: MuMin                   ! Smallest shrinking of dx before change direction
    type(matrix_T), pointer :: OutputAverage      ! Averaging Kernel
    type(matrix_SPD_T), pointer :: OutputCovariance    ! Covariance of the sol'n
    type(vector_T), pointer :: OutputSD ! Vector containing SD of result
    logical :: ParallelMode             ! Run forward models in parallel
    character(len=127) :: PhaseName     ! To pass to snoopers
    integer :: Son                      ! Of Root or Key
    character(len=127) :: SnoopComment  ! From comment= field of S_Snoop spec.
    integer :: SnoopKey                 ! Tree point of S_Snoop spec.
    integer :: SnoopLevel               ! From level field of S_Snoop spec.
    integer, dimension(:), pointer :: SparseQuantities ! Which jacobian blocks to sparsify
    integer :: Spec                     ! s_matrix, s_subset or s_retrieve ... 
    type(vector_T), pointer :: State    ! The state vector
    real :: T0, T1, T2, T3              ! for timing
    type(matrix_T) :: Tikhonov          ! Matrix for Tikhonov regularization
    logical :: TikhonovApriori          ! Regularization is to aproiri, not
                                        ! zero -- default false
    logical :: TikhonovBefore           ! "Do Tikhonov before column scaling"
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
    integer :: VRegOrders                ! Regularization orders
    integer :: VRegQuants                ! Regularization quantities
    integer :: VRegWeights               ! Weight of regularization conditions
    type(vector_T), pointer :: VRegWeightVec  ! Weight vector for regularization

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
    integer, parameter :: Reg_RHS = reg_X_x + 1    ! RHS for Tikhonov if regApriori
    integer, parameter :: Weight = reg_RHS + 1     ! Scaling vector for rows, 1/measurementSD
    integer, parameter :: X = weight + 1           ! for NWT
    integer, parameter :: LastVec = X

    type(vector_T), dimension(firstVec:lastVec) :: V   ! Database for snoop

    ! Error message codes
    integer, parameter :: BadOpticalDepthSignal = 1
    integer, parameter :: BadOpticalDepthQuantities = BadOpticalDepthSignal + 1
    ! Only one of two required fields supplied
    integer, parameter :: BothOrNeither = BadOpticalDepthQuantities + 1
    integer, parameter :: IfAThenB = BothOrNeither + 1
    integer, parameter :: IfUnitsAThenB = IfAThenB + 1
    integer, parameter :: Inconsistent = IfUnitsAThenB + 1  ! Inconsistent fields
    integer, parameter :: InconsistentUnits = Inconsistent + 1
    integer, parameter :: NeedBothDepthAndCutoff = InconsistentUnits + 1
    integer, parameter :: NoFields = NeedBothDepthAndCutoff + 1  ! No fields are allowed
    integer, parameter :: NotGeneral = noFields + 1  ! Not a general matrix
    integer, parameter :: NotSPD = notGeneral + 1    ! Not symmetric pos. definite
    integer, parameter :: RangeNotAppropriate = NotSPD + 1
    integer, parameter :: WrongUnits = RangeNotAppropriate + 1
    integer, parameter :: MustHaveOne = WrongUnits + 1
    integer, parameter :: CannotFlagCloud = MustHaveOne + 1
    integer, parameter :: BadQuantities = CannotFlagCloud + 1
    integer, parameter :: BadChannel = BadQuantities + 1

    error = 0
    nullify ( apriori, configIndices, covariance, fwdModelOut )
    nullify ( measurements, measurementSD, state, outputSD, sparseQuantities )
    phaseName = ' '              ! Default in case there's no field
    snoopComment = ' '           ! Ditto
    snoopKey = 0
    snoopLevel = 1               ! Ditto
    timing = section_times
    do j = firstVec, lastVec ! Make the vectors in the database initially empty
      nullify ( v(j)%quantities, v(j)%template%quantities )
      v(j)%name = 0 ! so Snoop won't use it
    end do
    if ( timing ) call time_now ( t1 )

    if ( toggle(gen) ) call trace_begin ( "Retrieve", root )
    do i_sons = 2, nsons(root) - 1      ! skip names at begin/end of section
      son = subtree(i_sons, root)
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
        do i_key = 2, nsons(key)
          son = subtree(i_key, key)
          field = get_field_id(son)  ! tree_checker prevents duplicates
          select case ( field )
          case ( f_comment )
            call get_string ( sub_rosa(subtree(2,son)), snoopComment, strip=.true. )
          case ( f_phaseName )
            call get_string ( sub_rosa(subtree(2,son)), phaseName, strip=.true. )
          case ( f_level )
            call expr ( subtree(2,son), units, value, type )
            if ( units(1) /= phyq_dimensionless ) &
              & call announceError ( wrongUnits, field, string='no' )
            snoopLevel = nint(value(1))
          end select
        end do
      case ( s_subset )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Retrieve.subset", root )
        call SetupSubset ( key, vectorDatabase )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Retrieve.subset" )
      case ( s_restrictRange )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Retrieve.RestrictRange", root )
        call RestrictRange ( key, vectorDatabase )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Retrieve.RestrictRange" )
      case ( s_flagCloud )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Retrieve.flagCloud", root )
        call SetupflagCloud ( key, vectorDatabase )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Retrieve.flagCloud" )
      case ( s_updateMask )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_begin ( "Retrieve.UpdateMask", root )
        call UpdateMask ( key, vectorDatabase )
        if ( toggle(gen) .and. levels(gen) > 0 ) &
          & call trace_end ( "Retrieve.UpdateMask" )
      case ( s_retrieve )
        if ( toggle(gen) ) call trace_begin ( "Retrieve.retrieve", root )
        aprioriScale = 1.0
        columnScaling = l_none
        covSansReg = .false.
        diagonal = .false.
        extendedAverage = .false.
        hRegQuants = 0
        hRegWeights = 0
        nullify ( hRegWeightVec )
        initLambda = defaultInitLambda
        maxJacobians = defaultMaxJ
        method = defaultMethod
        muMin = defaultMuMin
        toleranceA = defaultToleranceA
        toleranceF = defaultToleranceF
        toleranceR = defaultToleranceR
        vRegQuants = 0
        vRegWeights = 0
        parallelMode = parallel%fwmParallel .and. parallel%master
        tikhonovApriori = .false.
        tikhonovBefore = .true.
        nullify ( vRegWeightVec )
        do i_key = 2, nsons(key) ! fields of the "retrieve" specification
          son = subtree(i_key, key)
          field = get_field_id(son)  ! tree_checker prevents duplicates
          got(field) = .true.
          select case ( field )
          case ( f_apriori )
            apriori => vectorDatabase(decoration(decoration(subtree(2,son))))
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
          case ( f_extendedAverage )
            extendedAverage = get_boolean(son)
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
          case ( f_highBound )
            highBound => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_hRegOrders )
            hRegOrders = son
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
              sparseQuantities(k-1) = decoration(decoration(subtree(k,son)))
            end do
          case ( f_state )
            state => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_aprioriScale, f_fuzz, f_lambda, f_maxJ, f_muMin, &
            &    f_toleranceA, f_toleranceF, f_toleranceR )
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
            case ( f_maxJ )
              maxJacobians = value(1)
            case ( f_muMin )
              muMin = value(1)
            case ( f_toleranceA )
              toleranceA = value(1)
            case ( f_toleranceF )
              toleranceF = value(1)
            case ( f_toleranceR )
              toleranceR = value(1)
            end select
          case ( f_vRegOrders )
            vRegOrders = son
          case ( f_vRegQuants )
            vRegQuants = son
          case ( f_vRegWeights )
            vRegWeights = son
          case ( f_vRegWeightVec )
            vRegWeightVec => vectorDatabase(decoration(decoration(subtree(2,son))))
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! i_key = 2, nsons(key)

        if ( got(f_apriori) .neqv. got(f_covariance) ) &
          & call announceError ( bothOrNeither, f_apriori, f_covariance )
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
        if ( error == 0 ) then

          ! Verify the consistency of various matrices and vectors
          if ( got(f_apriori) ) then
            if ( apriori%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_apriori, f_state )
          end if
          if ( got(f_highBound) ) then
            if ( highBound%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_highBound, f_state )
          end if
          if ( got(f_lowBound) ) then
            if ( lowBound%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_lowBound, f_state )
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

          ! Create the averaging kernel matrix
          if ( got(f_average) ) then
            k = decoration(ixAverage)
            if ( k == 0 ) then
              call createEmptyMatrix ( myAverage, &
                & sub_rosa(subtree(1,ixAverage)), state, state )
              k = addToMatrixDatabase( matrixDatabase, myAverage )
              call decorate ( ixAverage, k )
            end if
            call getFromMatrixDatabase ( matrixDatabase(k), outputAverage )
            call display_string ( outputAverage%name, advance='yes' )
            if ( .not. associated(outputAverage) ) then
              call announceError ( notGeneral, f_average )
            else
              if ( outputAverage%row%vec%template%id /= state%template%id &
                & .or. &
                &  outputAverage%col%vec%template%id /= state%template%id ) &
                & call announceError ( inconsistent, f_average, f_state )
            end if
          end if

          ! Create the matrix for adding the Tikhonov regularization to the
          ! normal equations
          if ( got(f_hRegOrders) .or. got(f_vRegOrders) ) then
            call createEmptyMatrix ( tikhonov, 0, state, state, text='_Tikhonov' )
            k = addToMatrixDatabase( matrixDatabase, tikhonov )
          end if
        end if
        if ( error == 0 ) then
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
            call add_to_retrieval_timing( 'cloud_retrieval', t1 )
          case ( l_newtonian ) 
            call newtonianSolver
          case default
            call MLSMessage ( MLSMSG_Error, moduleName, &
            & "this retrieval method has not yet been implemented" )
          end select

          !??? Make sure the jacobian and outputCovariance get destroyed
          !??? after ?what? happens?  Can we destroy the entire matrix
          !??? database at the end of each chunk?
          if ( .not. got(f_jacobian) ) call destroyMatrix ( jacobian )
          if ( .not. got(f_outputCovariance) ) &
            & call destroyMatrix ( outputCovariance%m )
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & "No retrieval done -- error in configuration" )
        end if
        if ( got(f_fwdModelOut) ) then
          call copyVector ( fwdModelOut, v(f) )
          call copyVectorMask ( fwdModelOut, measurements )
        end if
        call deallocate_test ( configIndices, "ConfigIndices", moduleName )
        if ( toggle(gen) ) call trace_end ( "Retrieve.retrieve" )
      case ( s_sids )
        call time_now ( t1 )
        call sids ( key, VectorDatabase, MatrixDatabase, configDatabase, chunk)
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      end select

      do j = firstVec, lastVec
        call destroyVectorInfo ( v(j) )
      end do

    end do ! i_sons = 2, nsons(root) - 1
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
      call output ( ', RetrievalModule complained: ' )
      select case ( code )
      case ( badChannel)
        call output ( 'Number of cloud channels must be 1')
      case ( cannotFlagCloud )
        call output ( 'Cannot flag clouds, missing input quantities' )
      case ( badQuantities )
        call output ( 'Bad quantity type for radiance or cloud radiance' )
      case ( badOpticalDepthSignal )
        call output ( 'Mismatch in signal/sideband for radiance and optical depth', &
          & advance='yes' )
      case ( badOpticalDepthQuantities )
        call output ( 'Bad quantity type for radiance or optical depth' )
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
      case ( inconsistent, notGeneral, notSPD )
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
      case ( inconsistentUnits )
        call output ( 'the elements of the "' )
        call display_string ( field_indices(fieldIndex) )
        call output ( '" field have inconsistent units.', advance='yes' )
      case ( needBothDepthAndCutoff )
        call output ( 'OpticalDepth specified but no cutoff, or visa versa', &
          & advance='yes' )
      case ( noFields )
        call output ( 'No fields are allowed for a ' )
        call display_string ( spec_indices(fieldIndex) )
        call output ( ' specification.', advance='yes' )
      case ( rangeNotAppropriate )
        call output ( 'Ranges are not appropriate for a ' )
        call display_string ( spec_indices(fieldIndex) )
        call output ( ' specification.', advance='yes' )
      case ( wrongUnits )
        call output ( 'The value(s) of the "' )
        call display_string ( field_indices(fieldIndex) )
        call output ( '" field shall have ' )
        call output ( trim(string) )
        call output ( ' units.', advance='yes' )
      case ( mustHaveOne )
        call output ( 'Must supply one and only one of height, ignore or reset', &
          & advance='yes' )
      end select
    end subroutine AnnounceError

    ! ------------------------------------------------  BoundMove  -----
    subroutine BoundMove ( Mu, Bound, X, Dx, Which, muMin )
      ! Compute a scalar multiple Mu such that X + Mu*DX does not go
      ! beyond Bound.  "Which" specifies either "low" or "high"

      !{ Let $D$ be $|\delta \mathbf{x}|$, where $\delta \mathbf{x}$ is
      ! given by Dx, let $b_i$ and $x_i$ be components of Bound and X. 
      ! Let $d_i$ be such that $x_i + \frac{d_i}D \delta x_i$ is not beyond
      ! a bound.  That is, $d_i$ is such that $d_i \cos \theta_i = x_i -
      ! b_i$, where $\theta_i$ is the angle between $\delta \mathbf{x}$
      ! and the normal to the $i^{th}$ constraint surface.  Since the
      ! constraints are bounds, we have $\cos \theta_i = \frac{\delta
      ! x_i}D$.  $\mu$ is the smallest value of $\frac{d_i}D = \frac{x_i-b_i}
      ! {\delta x_i}$.  If we check the components one at a time, each one
      ! using the most recently computed $\mu$, we don't need a divide for
      ! each check, and we don't need a {\tt min} operation.
      ! However, if we are already at the bounds, and the next step 
      ! wants to keep going beyond the bounds, then $\mu$ will be
      ! really small.  In this case $\mu < \mu_{\text{min}}$ we revert
      ! to a straight element-by-element modification of $\delta x_i$.      

      real(rv), intent(inout) :: Mu
      type(vector_T), intent(inout) :: Bound, X, Dx
      character(len=*), intent(in) :: Which
      real(rv), intent(in) :: MuMin

      integer :: IQ, IVX, IVY           ! Subscripts used during MU computation

      if ( which == 'low' ) then
        do iq = 1, size(x%quantities)
          if ( associated(x%quantities(iq)%mask) ) then
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      & mu = ( bound%quantities(iq)%values(ivx,ivy) - &
                      &        x%quantities(iq)%values(ivx,ivy) ) &
                      &      / dx%quantities(iq)%values(ivx,ivy)
                end if
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( x%quantities(iq)%values(ivx,ivy) + &
                  &  mu * dx%quantities(iq)%values(ivx,ivy) < &
                  &  bound%quantities(iq)%values(ivx,ivy) ) &
                    & mu = ( bound%quantities(iq)%values(ivx,ivy) - &
                    &        x%quantities(iq)%values(ivx,ivy) ) &
                    &      / dx%quantities(iq)%values(ivx,ivy)
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
                      & mu = ( x%quantities(iq)%values(ivx,ivy) - &
                      &        bound%quantities(iq)%values(ivx,ivy) ) &
                      &      / dx%quantities(iq)%values(ivx,ivy)
                end if
              end do
            end do
          else
            do ivy = 1, size(x%quantities(iq)%values,2)
              do ivx = 1, size(x%quantities(iq)%values,1)
                if ( x%quantities(iq)%values(ivx,ivy) + &
                  &  mu * dx%quantities(iq)%values(ivx,ivy) > &
                  &  bound%quantities(iq)%values(ivx,ivy) ) &
                    & mu = ( x%quantities(iq)%values(ivx,ivy) - &
                    &        bound%quantities(iq)%values(ivx,ivy) ) &
                    &      / dx%quantities(iq)%values(ivx,ivy)
              end do
            end do
          end if
        end do
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Come on! In BoundMove, it has to be a low bound or a high bound!' )
      end if

      ! Now if mu has goten really small, we'll change it back to one,
      ! and do an element by element modification of dx
      if ( mu < muMin ) then
        mu = 1.0_rv
        if ( which == 'low' ) then
          do iq = 1, size(x%quantities)
            if ( associated(x%quantities(iq)%mask) ) then
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( iand(ichar(x%quantities(iq)%mask(ivx,ivy)),m_linalg) == 0 ) then
                    if ( x%quantities(iq)%values(ivx,ivy) + &
                      &  dx%quantities(iq)%values(ivx,ivy) < &
                      &  bound%quantities(iq)%values(ivx,ivy) ) &
                        &  dx%quantities(iq)%values(ivx,ivy) = &
                        &    bound%quantities(iq)%values(ivx,ivy) - &
                        &    x%quantities(iq)%values(ivx,ivy)
                  end if
                end do
              end do
            else
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  dx%quantities(iq)%values(ivx,ivy) < &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      &  dx%quantities(iq)%values(ivx,ivy) = &
                      &    bound%quantities(iq)%values(ivx,ivy) - &
                      &    x%quantities(iq)%values(ivx,ivy)
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
                      &  dx%quantities(iq)%values(ivx,ivy) > &
                      &  bound%quantities(iq)%values(ivx,ivy) ) &
                        &  dx%quantities(iq)%values(ivx,ivy) = &
                        &    bound%quantities(iq)%values(ivx,ivy) - &
                        &    x%quantities(iq)%values(ivx,ivy)
                  end if
                end do
              end do
            else
              do ivy = 1, size(x%quantities(iq)%values,2)
                do ivx = 1, size(x%quantities(iq)%values,1)
                  if ( x%quantities(iq)%values(ivx,ivy) + &
                    &  dx%quantities(iq)%values(ivx,ivy) > &
                    &  bound%quantities(iq)%values(ivx,ivy) ) &
                      &  dx%quantities(iq)%values(ivx,ivy) = &
                      &    bound%quantities(iq)%values(ivx,ivy) - &
                      &    x%quantities(iq)%values(ivx,ivy)
                end do
              end do
            end if
          end do
        end if
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
            &  'How did mu get to be negative?' )
        end if
      end if
    end subroutine BoundMove

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
      type(vectorValue_T), pointer :: Diag_Qty    ! A quantity in the
                                        ! Diagnostics vector
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
    end subroutine FillDiagQty

    ! ----------------------------------------------  FillDiagVec ------
    subroutine FillDiagVec ( diagnostics, aj, numJ, nwt_flag, &
      & jacobian_rows, jacobian_cols )
      use DNWT_Module, only: NWT_T
      type(Vector_T), intent(inout) :: DIAGNOSTICS
      type(Nwt_T), intent(in) :: AJ
      integer, intent(in), optional :: NUMJ
      integer, intent(in), optional :: NWT_FLAG
      integer, intent(in), optional :: JACOBIAN_ROWS
      integer, intent(in), optional :: JACOBIAN_COLS

      call fillDiagQty ( diagnostics,  l_dnwt_ajn, aj%ajn )
      call fillDiagQty ( diagnostics,  l_dnwt_axmax, aj%axmax )
      call fillDiagQty ( diagnostics,  l_dnwt_cait, aj%cait )
      call fillDiagQty ( diagnostics,  l_dnwt_chiSqMinNorm, aj%chiSqMinNorm )
      call fillDiagQty ( diagnostics,  l_dnwt_chiSqNorm, aj%chiSqNorm )
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
      if ( present ( numJ ) ) &
        & call fillDiagQty ( diagnostics,  l_numJ, real(numJ, r8) )
      if ( present ( nwt_flag ) ) &
        & call fillDiagQty ( diagnostics,  l_dnwt_flag, real(nwt_flag,r8) )
      if ( present ( jacobian_rows ) ) &
        & call fillDiagQty ( diagnostics,  l_jacobian_rows, real(jacobian_rows,r8) )
      if ( present ( jacobian_cols ) ) &
        & call fillDiagQty ( diagnostics,  l_jacobian_cols, real(jacobian_cols,r8) )
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

    ! ------------------------------------------  NewtonianSolver  -----
    subroutine NewtonianSolver

      use DNWT_Module, only: FlagName, NF_AITKEN, NF_BEST, NF_DX, &
      & NF_DX_AITKEN, NF_EVALF, NF_EVALJ, NF_FANDJ, NF_GMOVE, NF_LEV, &
      & NF_NEWX, NF_SMALLEST_FLAG, NF_SOLVE, NF_START, NF_TOLX, NF_TOLF, &
      & NF_TOLX_BEST, NF_TOO_SMALL, NWT, NWTA, NWTDB, NWTOP, NWT_T, RK
      use Dump_0, only: Dump
      use ForwardModelWrappers, only: ForwardModel
      use ForwardModelIntermediate, only: ForwardModelIntermediate_T, &
        & ForwardModelStatus_T
      use MatrixModule_1, only: AddToMatrix, CholeskyFactor, ClearMatrix, &
        & ColumnScale, CopyMatrixValue, CreateEmptyMatrix, &
        & DestroyMatrix, Dump, Dump_Linf, Dump_struct, &
        & FormNormalEquations => NormalEquations, &
        & GetDiagonal, InvertCholesky, Matrix_T, &
        & Matrix_Cholesky_T, Matrix_SPD_T, MaxL1, MinDiag, Multiply, &
        & MultiplyMatrix_XY_T,  RowScale, ScaleMatrix, SolveCholesky, &
        & UpdateDiagonal
      use Regularization, only: Regularize
      use Symbol_Table, only: ENTER_TERMINAL
      use Symbol_Types, only: T_IDENTIFIER
      use VectorsModule, only: AddToVector, DestroyVectorInfo, &
        & Dump, Multiply, operator(.DOT.), &
        & operator(.MDOT.), operator(-), ScaleVector, SubtractFromVector
      use L2FWMParallel, only: SETUPFWMSLAVES, TRIGGERSLAVERUN, &
        & REQUESTSLAVESOUTPUT, RECEIVESLAVESOUTPUT

      ! Local Variables
      type(nwt_T) :: AJ                 ! "About the Jacobian", see NWT.
      type(nwt_T) :: BestAJ             ! AJ at Best Fnorm so far.
      real(r8) :: Cosine                ! Of an angle between two vectors
      type(matrix_Cholesky_T) :: Factored ! Cholesky-factored normal equations
      type (ForwardModelStatus_T) :: FmStat ! Status for forward model
      type (ForwardModelIntermediate_T) :: Fmw ! Work space for forward model
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
                                        
      integer :: LATESTMAFSTARTED       ! For FWMParallel stuff
      real(rv) :: MU                    ! Move Length = scale for DX
      integer, parameter :: NF_GetJ = NF_Smallest_Flag - 1 ! Take an extra loop
                                        ! to get J.
      type(matrix_SPD_T), target :: NormalEquations  ! Jacobian**T * Jacobian
      integer :: NumJ                   ! Number of Jacobian evaluations
      integer :: NWT_Flag               ! Signal from NWT, q.v., indicating
                                        ! the action to take.
      integer :: NWT_Opt(20)            ! Options for NWT, q.v.
      real(rk) :: NWT_Xopt(20)          ! Real parameters for NWT options, q.v.
      integer :: PreserveMatrixName     ! Temporary name store
      integer :: Prev_NWT_Flag          ! Previous value of NWT_Flag
      type(vector_T) :: Q               ! Used to calculate Marquardt parameter
      integer :: RowBlock               ! Which block of rows is the forward
                                        ! model filling?
      integer, parameter :: SnoopLevels(NF_DX_AITKEN:NF_FANDJ) = (/ &
      ! dx_aitken dx aitken best gmove newx lev solve evalj evalf
        &      3, 2,     3,   2,    2,   1,  2,   2,    2,    2,  &
      ! start tolx tolx_best tolf too_small fandj
        &  9,   2,        2,   2,        3,    9  /)
      integer :: T                      ! Which Tikhonov: 1 -> V, 2 -> H
      real :: T1
      type(matrix_T) :: Temp            ! Because we can't do X := X * Y
      character(len=10) :: TheFlagName  ! Name of NWTA's flag argument
      integer :: TikhonovRows           ! How many rows of Tiknonov regularization?

      call time_now ( t1 )
      call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
        & ModuleName )
      ! Launch fwmParallel slaves
      if ( parallelMode ) &
        & call SetupFWMSlaves ( configDatabase(configIndices), &
        & state, fwdModelExtra, FwdModelOut, jacobian )
      foundBetterState = ( maxJacobians == 0 )
      ! Set options for NWT
      nwt_opt(1:9) = (/  15, 1,      17, 2,      18, 3,      11, 4, 0 /)
      nwt_xopt(1:4) = (/ toleranceF, toleranceA, toleranceR, initLambda /)
      call nwt ( nwt_flag, nwt_xopt, nwt_opt )
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
      call copyVector ( v(x), state, vectorNameText='_x', clone=.true. ) ! x := state
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
        if ( index(switches,'xvec') /= 0 ) call dump ( v(x), name='Original X' )
      numJ = 0
      aj%axmax = 0.0
        call time_now ( t0 ) ! time base for Newtonian iteration
      do k = 1, size(v(x)%quantities)
        aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(k)%values)))
      end do

        if ( index(switches,'sca') /= 0 ) then
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

      call copyVector ( v(bestX), v(x) ) ! bestX := x to start things off
      prev_nwt_flag = huge(0)
      do ! Newtonian iteration
        if ( nwt_flag /= nf_getJ ) then ! not taking a special iteration to get J
            if ( index(switches,'nin') /= 0 ) & ! Turn on NWTA's internal output
              & call nwtop ( (/ 1, 1, 0 /), nwt_xopt )
          call nwta ( nwt_flag, aj )
          if ( nwt_flag == nf_evalj .and. prev_nwt_flag == nf_evalf ) then
            prev_nwt_flag = nf_evalj
            cycle
          end if
        end if
          if ( index(switches,'nwt') /= 0 ) then
            call FlagName ( nwt_flag, theFlagName )
            if ( nwt_flag == nf_getJ ) theFlagName = 'GETJ'
            call output ( 'Newton method flag = ' )
            call output ( trim(theFlagName) )
            call output ( ', numJ = ' )
            call output ( numJ )
            call output ( ' at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
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
          if ( index(switches,'dst') /= 0 ) call dump ( v(x), details=2, name='State' )
          if ( numJ > maxJacobians .and. nwt_flag /= nf_getJ ) then
              if ( index(switches,'nwt') /= 0 ) then
                call output ( &
                  & 'Newton iteration terminated because Jacobian evaluations (' )
                call output ( numJ )
                call output ( ') > maxJacobians (' )
                call output ( maxJacobians )
                call output ( ')', advance='yes' )
            if ( .not. foundBetterState ) exit
              end if
            ! Restore BestX, run the forward model one more time to get a new
            ! Jacobian, and form normal equations -- the last two so that the
            ! a posteriori covariance is consistent with BestX.
            aj = bestAJ
            call copyVector ( v(x), v(bestX) ) ! x := bestX
            if ( .not. (got(f_outputCovariance) .or. got(f_outputSD) .or. &
              & got(f_average)) ) exit
            nwt_flag = nf_getJ
          end if

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
            if ( index ( switches, 'fnorm' ) /= 0 ) then
              call output ( 'A priori contribution to | F | = ' )
              call output ( aj%fnorm, advance='yes' )
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
            if ( index ( switches, 'fnorm' ) /= 0 ) then
              call output ( 'No a priori so setting | F | to zero.', advance='yes' )
            end if
            call clearMatrix ( normalEquations%m ) ! start with zero
            call clearVector ( v(aTb) ) ! Clear the RHS vector
          end if

          ! Add Tikhonov regularization if requested.
          if ( (got(f_hRegOrders) .or. got(f_vRegOrders)) .and. &
            & tikhonovBefore .and. &
            & ( ( nwt_flag /= nf_getj ) .or. .not. covSansReg ) ) then
              call add_to_retrieval_timing( 'newton_solver', t1 )

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
              if ( tikhonovApriori ) then ! reg_*_x := R * (apriori - x)
                call multiplyMatrixVectorNoT ( tikhonov, v(reg_RHS), v(reg_X_x) )
              else                          ! reg_*_x := -R * x
                call multiplyMatrixVectorNoT ( tikhonov, v(x), v(reg_X_x) )
                call scaleVector ( v(reg_X_x), -1.0_r8 )   ! -R x_n
              end if
                if ( index(switches,'reg') /= 0 ) then
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
              call formNormalEquations ( tikhonov, normalEquations, &
                & v(reg_X_x), v(aTb), update=update, useMask=.false. )
              update = .true.
              call clearMatrix ( tikhonov )           ! free the space
              ! aj%fnorm is still the square of the norm of f
              aj%fnorm = aj%fnorm + ( v(reg_X_x) .dot. v(reg_X_x) )
              if ( index ( switches, 'fnorm' ) /= 0 ) then
                if ( t == 1 ) then
                  call output ( 'Vertical Regularization contribution to | F | = ' )
                else
                  call output ( 'Horizontal Regularization contribution to | F | = ' )
                end if
                call output ( v(reg_X_x) .dot. v(reg_X_x), advance='yes' )
              end if

              ! call destroyVectorValue ( v(reg_X_x) )  ! free the space
              ! Don't destroy reg_X_x unless we move the 'clone' for it
              ! inside the loop.  Also, if we destroy it, we can't snoop it.
                call add_to_retrieval_timing( 'tikh_reg', t1 )
            end do ! t
          else
            tikhonovRows = 0
          end if

          fmStat%maf = 0
          fmStat%newScanHydros = .true.

          if ( index(switches,'vir') /=0 ) call dump ( normalEquations%m, details=2 )

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
            if ( index ( switches, 'mas' ) /= 0 ) &
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
              ! Otherwise, we call the forward model as ususal
              do k = 1, size(configIndices)
                call forwardModel ( configDatabase(configIndices(k)), &
                  & v(x), fwdModelExtra, v(f_rowScaled), fmw, fmStat, jacobian )
              end do ! k
              ! Forward model calls add_to_retrieval_timing
            end if
            call time_now ( t1 )
            do rowBlock = 1, size(fmStat%rows)
              if ( fmStat%rows(rowBlock) ) then
                 ! Store what we've just got in v(f) ie fwdModelOut
                call copyVector ( v(f), v(f_rowScaled), & ! v(f) := v(f_rowScaled)
                  & quant=jacobian%row%quant(rowBlock), &
                  & inst=jacobian%row%inst(rowBlock) )
                call subtractFromVector ( v(f_rowScaled), measurements, &
                  & quant=jacobian%row%quant(rowBlock), &
                  & inst=jacobian%row%inst(rowBlock) ) ! f - y
                !{Let $\bf W$ be the Cholesky factor of the inverse of the
                ! measurement covariance ${\bf S}_m$ (which in our case is
                ! diagonal), i.e. ${\bf W}^T {\bf W} = {\bf S}_m^{-1}$. Row
                ! scale the part of the least-squares problem that arises
                ! from the measurements, i.e. the least-squares problem
                ! becomes $\mathbf{W J \delta \hat x -\simeq W f}$ (actually,
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
            ! -{\bf J}^T {\bf S}_m^{-1} {\bf f}$:
              call add_to_retrieval_timing( 'newton_solver', t1 )
            call formNormalEquations ( jacobian, kTk, rhs_in=v(f_rowScaled), &
              & rhs_out=v(aTb), update=update, useMask=.true. )
              if ( index ( switches, 'atb' ) /= 0 ) call dump ( v(aTb) )
              call add_to_retrieval_timing( 'form_normeq', t1 )
            if ( got(f_average) .and. extendedAverage ) then
              call MultiplyMatrix_XTY ( jacobian, jacobian, kTkStar, update=update, &
                & maskX=.true., maskY=.false. )
            end if
            update = .true.
              if ( index(switches,'jac') /= 0 ) &
                call dump_Linf ( jacobian, 'L_infty norms of Jacobian blocks:' )
              if ( index(switches,'spa') /= 0 ) &
                & call dump_struct ( jacobian, &
                  & 'Sparseness structure of Jacobian blocks:' )
            call clearMatrix ( jacobian )  ! free the space
          end do ! mafs

          ! In the case where we're forming a regular averaging kernel
          ! (not extended), we need to add this MAFs worth of kTk into the
          ! normal equations.  Otherwise normalEquations and kTk are in
          ! fact pointers to the same matrix
          if ( got(f_average) .and. .not. extendedAverage ) then
            call addToMatrix ( normalEquations%m, kTk%m )
          end if

          ! aj%fnorm is still the square of the function norm
          aj%fnorm = aj%fnorm + ( v(f_rowScaled) .mdot. v(f_rowScaled) )
          if ( index ( switches, 'fnorm' ) /= 0 ) then
            call output ( 'Measurement contribution to | F | = ' )
            call output ( v(f_rowScaled) .mdot. v(f_rowScaled), advance='yes' )
          end if

          ! Add Tikhonov regularization if requested.  We do it here instead
          ! of before adding the Jacobian so that we can scale it up by the
          ! column scaling before scaling the normal equations back down
          ! by the column scaling.  The effect is that regularization takes
          ! place (roughly) on the column-scaled problem, and if scaling by
          ! the column norm is requested, it's still makes the column norms
          ! all 1.0.
          if ( (got(f_hRegOrders) .or. got(f_vRegOrders)) .and. &
            & .not. tikhonovBefore .and. &
            & ( ( nwt_flag /= nf_getj ) .or. .not. covSansReg ) ) then
              call add_to_retrieval_timing( 'newton_solver', t1 )

            !{ Tikhonov regularization is of the form ${\bf R x}_{n+1} \simeq
            !  {\bf 0}$. So that all of the parts of the problem are solving
            !  for ${\bf\delta x}$, we subtract ${\bf R x}_n$ from both sides
            !  to get ${\bf R \delta x} \simeq -{\bf R x}_n$.

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

              if ( tikhonovApriori ) then
                call multiplyMatrixVectorNoT ( tikhonov, v(reg_RHS), v(reg_X_x) )
              else
                call multiplyMatrixVectorNoT ( tikhonov, v(x), v(reg_X_x) )
                call scaleVector ( v(reg_X_x), -1.0_r8 )   ! -R x_n
              end if
              call scaleVector ( v(reg_X_x), -1.0_r8 )   ! -R x_n

              if ( columnScaling /= l_none ) then ! Compute $\Sigma$
                forall (j = 1: v(columnScaleVector)%template%noQuantities)
                  where ( v(columnScaleVector)%quantities(j)%values <= 0.0 )
                    v(columnScaleVector)%quantities(j)%values = 0.0
                  elsewhere
                    v(columnScaleVector)%quantities(j)%values = &
                      & sqrt( v(columnScaleVector)%quantities(j)%values )
                  end where
                end forall
                call columnScale ( tikhonov, v(columnScaleVector) )
              end if

                if ( index(switches,'reg') /= 0 ) then
                  call dump_struct ( tikhonov, 'Tikhonov' )
                  call dump ( tikhonov, name='Tikhonov', details=2 )
                end if
              call formNormalEquations ( tikhonov, normalEquations, &
                & v(reg_X_x), v(aTb), update=update, useMask=.false. )
              update = .true.
              call clearMatrix ( tikhonov )           ! free the space
              ! aj%fnorm is still the square of the norm of f
              aj%fnorm = aj%fnorm + ( v(reg_X_x) .dot. v(reg_X_x) )
              if ( index ( switches, 'fnorm' ) /= 0 ) then
                if ( t == 1 ) then
                  call output ( 'Vertical Regularization contribution to | F | = ' )
                else
                  call output ( 'Horizontal Regularization contribution to | F | = ' )
                end if
                call output ( v(reg_X_x) .dot. v(reg_X_x), advance='yes' )
              end if
              ! call destroyVectorValue ( v(reg_X_x) )  ! free the space
              ! Don't destroy reg_X_x unless we move the 'clone' for it
              ! inside the loop.  Also, if we destroy it, we can't snoop it.
                call add_to_retrieval_timing( 'tikh_reg', t1 )
            end do ! t
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
            if ( index(switches,'gvec') /= 0 ) &
              & call dump ( v(gradient), name='gradient' )
          
          !{Compute the Cholesky factor of the LHS of the normal equations:
          ! ${\bf U}^T {\bf U} {\bf \delta \hat x} = {\bf\Sigma}^T {\bf J}^T
          ! {\bf S}_m^{-1} {\bf J \Sigma \Sigma}^{-1} {\bf \delta \hat x} =
          ! -{\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$.

          ! factored%m%block => normalEquations%m%block ! to save space
          ! Can't do the above because we need to keep the normal
          ! equations around, in order to subtract Levenberg-Marquardt and
          ! apriori covariance, in order to compute a posteriori covariance
            if ( index(switches,'neq') /= 0 ) &
              call dump_Linf ( normalEquations%m, &
                & 'L_infty norms of Normal Equations blocks after scaling:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
            call add_to_retrieval_timing ( 'newton_solver', t1 )
            call choleskyFactor ( factored, normalEquations )
            call add_to_retrieval_timing ( 'cholesky_factor', t1 )
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
            if ( nwt_flag == nf_getJ ) then ! taking a special iteration to get J
              aj%chiSqNorm = aj%fnorm / max ( jacobian_rows - jacobian_cols, 1 )
              aj%fnorm = sqrt(aj%fnorm)
                if ( index(switches,'sca') /= 0 ) &
                  & call dump ( (/ aj%fnorm, aj%chiSqNorm /) , &
                  & ' | F |       chi^2/n ', clean=.true. )
              exit
          end if
          aj%diag = minDiag ( factored ) ! element on diagonal with
            !       smallest absolute value, after triangularization
          aj%ajn = maxL1 ( factored%m ) ! maximum L1 norm of
          !       column in upper triangle after triangularization
            call add_to_retrieval_timing ( 'newton_solver', t1 )
          call solveCholesky ( factored, v(candidateDX), v(aTb), &
            & transpose=.true. )
            call add_to_retrieval_timing ( 'cholesky_solver', t1 )

          !{AJ\%FNMIN = $L_2$ norm of residual, $||{\bf\Sigma}^T {\bf J}^T
          ! {\bf S}_m^{-1} {\bf J \Sigma \Sigma}^{-1} {\bf \delta \hat x} +
          ! {\bf \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}||$ where
          ! $\bf\delta \hat x$ is the ``Candidate DX'' that may not get
          ! used. This can be gotten without saving $\bf J$ as ${\bf \hat f}^T
          ! {\bf \hat f} - {\bf y}^T {\bf y}$ where ${\bf\hat f} = {\bf
          ! \Sigma}^T {\bf J}^T {\bf S}_m^{-1} {\bf f}$ and ${\bf y} = {\bf
          ! U}^{-T} {\bf \hat f}$. The variable {\tt candidateDX} is a
          ! temp here = $\bf y$.
          aj%fnmin = aj%fnorm - (v(candidateDX) .dot. v(candidateDX))
          if ( aj%fnmin < 0.0 ) then
            call output ( 'How can aj%fnmin be negative?  aj%fnmin = ' )
            call output ( aj%fnmin, advance='yes' )
            call output ( 'aj%fnorm = ' )
            call output ( aj%fnorm, advance='yes' )
            call output ( 'norm(candidateDX) = ' )
            call output ( v(candidateDX) .dot. v(candidateDX), advance='yes' )
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Norm of residual is imaginary!' )
            aj%fnmin = tiny ( aj%fnmin )
          end if
          ! Compute number of rows of Jacobian actually used.  Don't count
          ! rows due to Levenberg-Marquardt stabilization.  Do count rows
          ! due to a priori or regularization.  Put numbers of rows and
          ! columns into diagnostic vector.
          jacobian_cols = sum(jacobian%col%nelts)
          jacobian_rows = sum(jacobian%row%nelts)
          do j = 1, measurements%template%noQuantities
            if ( associated(measurements%quantities(j)%mask) ) &
              & jacobian_rows = jacobian_rows - &
              &   countBits(measurements%quantities(j)%mask, what=m_linAlg )
          end do
          ! Correct for apriori information.  Note that there is an
          ! approximation here: We don't take any account of whether the a
          ! priori is used on an element by element basis.
          if ( got(f_apriori) ) &
            & jacobian_rows = jacobian_rows + jacobian_cols
          ! Correct for Tikhonov information.
          if ( any ( got( (/f_hRegOrders, f_vRegOrders /) ) ) ) &
            & jacobian_rows = jacobian_rows + tikhonovRows
          ! Compute the normalised chiSquared statistics etc.
          aj%chiSqMinNorm = aj%fnmin / max ( jacobian_rows - jacobian_cols, 1 )
          aj%chiSqNorm = aj%fnorm / max ( jacobian_rows - jacobian_cols, 1 )
          aj%fnmin = sqrt(aj%fnmin)
          aj%fnorm = sqrt(aj%fnorm)
          aj%gradn = sqrt(v(gradient) .dot. v(gradient)) ! L2Norm(gradient)
            if ( index(switches,'sca') /= 0 ) then
              call dump ( (/ aj%fnorm, aj%ajn, aj%diag, aj%fnmin, aj%gradn /), &
                & '     | F |       L1| FAC |     aj%diag        aj%fnmin       | G |', &
                & clean=.true. )
              call dump ( (/ aj%chiSqNorm, aj%chiSqMinNorm /), &
                & '  chi^2/n      chimin^2/n ', clean=.true. )
            end if
        case ( nf_lev ) ! ..................................  LEV  .....
        ! Solve U^T q = dx, then compute ||q||**2, which is used to
        ! compute the Levenberg-Marquardt stabilization parameter.
          call solveCholesky ( factored, q, v(candidateDX), &
            & transpose=.true. ) ! q := factored^{-T} v(candidateDX)
          aj%qnsq = q .dot. q
            if ( index(switches,'sca') /= 0 ) then
              call output ( 'QN^2 = ' )
              call output ( aj%qnsq, advance='yes' )
            end if
        case ( nf_solve ) ! ..............................  SOLVE  .....
        !{Apply Levenberg-Marquardt stabilization with parameter
        ! $\lambda =$ {\bf AJ\%SQ}.  I.e., form $({\bf \Sigma}^T {\bf J}^T
        ! {\bf W}^T {\bf W J \Sigma + \lambda^2 I}) {\bf \Sigma}^{-1}
        ! {\bf \delta \hat x = \Sigma}^T {\bf J}^T {\bf W}^T {\bf f}$ for
        ! ${\bf \Sigma}^{-1} {\bf \delta \hat x}$.  Set
        ! \begin{description}
        !   \item[AJ\%FNMIN] as for NWT\_FLAG = NF\_EVALJ, but taking
        !     account of Levenberg-Marquardt stabilization.  Actually,
        !     we don't have to compute this, because we did it at NF\_EVALF;
        !   \item[AJ\%DXN] = L2 norm of ``candidate DX'';
        !   \item[AJ\%GDX] = (Gradient) .dot. (``candidate DX'')
        ! \end{description}
          call updateDiagonal ( normalEquations, aj%sq**2 )
          ! factored%m%block => normalEquations%m%block ! to save space
          ! Can't do the above because we need to keep the normal equations
          ! around, in order to subtract Levenberg-Marquardt and apriori
          ! covariance, in order to compute a posteriori covariance
            if ( index(switches,'neq') /= 0 ) &
              call dump_Linf ( normalEquations%m, &
                & 'L1 norms of Normal Equations blocks after Marquardt:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( normalEquations%m, &
                & 'Sparseness structure of Normal equations blocks:', &
                & upper=.true. )
            call add_to_retrieval_timing( 'newton_solver', t1 )
          call choleskyFactor ( factored, normalEquations )
            call add_to_retrieval_timing( 'cholesky_factor', t1 )
            if ( index(switches,'fac') /= 0 ) &
              call dump_Linf ( factored%m, &
                & 'L1 norms of blocks of factor after Marquardt:', &
                & upper=.true. )
            if ( index(switches,'spa') /= 0 ) &
              & call dump_struct ( factored%m, &
                & 'Sparseness structure of blocks of factor:', upper=.true. )
          !{Solve for ``candidate DX'' = ${\bf \delta \hat x} = -({\bf J}^T {\bf
          ! J})^{-1} {\bf J}^T {\bf F} = -({\bf U}^T {\bf U})^{-1} {\bf J}^T 
          ! {\bf F}$ using two back solves.  First solve ${\bf U}^T {\bf y} =
          ! -{\bf J}^T {\bf F}$, then $\mathbf{U \delta \hat x}$ =
          ! $\mathbf{y}$. Meanwhile, set AJ\%FNMIN as for NWT\_FLAG =
          ! NF\_EVALJ, but taking account of Levenberg-Marquardt stabilization
          ! (actually, ``taking account of Levenberg-Marquardt stabilization''
          ! means ``not including Levenberg-Marquardt stabilization,'' which
          ! is how we did it at NF\_EVALF, so we don't need to do it now).
            call add_to_retrieval_timing( 'newton_solver', t1 )
          call solveCholesky ( factored, v(candidateDX), v(aTb), &
            & transpose=.true. ) ! v(candidateDX) := factored^{-T} v(aTb)
            call add_to_retrieval_timing( 'cholesky_solver', t1 )
          ! aj%fnorm is now the norm of f, not its square.
          ! The following calculation of fnmin was commented out, but on
          ! 11 September 2002 FTK told me it's the right thing to do
          ! after all.
          aj%fnmin = aj%fnorm**2 - (v(candidateDX) .dot. v(candidateDX))
          if ( aj%fnmin < 0.0 ) then
            call output ( 'How can aj%fnmin be negative?  aj%fnmin = ' )
            call output ( aj%fnmin, advance='yes' )
            call output ( 'aj%fnorm**2 = ' )
            call output ( aj%fnorm**2, advance='yes' )
            call output ( 'norm(candidateDX) = ' )
            call output ( v(candidateDX) .dot. v(candidateDX), advance='yes' )
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & "Norm of residual not in Jacobian's column space is imaginary!" )
            aj%fnmin = tiny ( aj%fnmin )
          end if
          aj%fnmin = sqrt(aj%fnmin)
          ! v(candidateDX) := factored^{-1} v(candidateDX)
          call solveCholesky ( factored, v(candidateDX) )
          ! Shorten candidateDX if necessary to stay within the bounds.
          ! We aren't ready to take the move, but NWTA needs an accurate
          ! value for the norm of candidateDX.
          call copyVector ( v(dxUnScaled), v(candidateDX) ) ! dxUnscaled = dx
          if ( columnScaling /= l_none ) then
            ! dxUnScaled = dxUnScaled # columnScaleVector:
            call multiply ( v(dxUnScaled), v(columnScaleVector) )
          end if
          mu = 1.0_rv
          if ( got(f_lowBound) ) &
            & call boundMove ( mu, lowBound, v(x), v(dxUnScaled), 'low', muMin )
          if ( got(f_highBound) ) &
            & call boundMove ( mu, highBound, v(x), v(dxUnScaled), 'high', muMin )
          if ( mu < 1.0_rv ) call scaleVector ( v(candidateDX), mu )
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
            if ( index(switches,'ndb') /= 0 ) &
              & call nwtdb ( width=9, level=0, why='After Solve' )
            if ( index(switches,'Ndb') /= 0 ) then
              if ( index(switches,'sca') /= 0 ) then
                call nwtdb ( width=9, why='After Solve' )
              else
                call nwtdb ( aj, width=9, why='After Solve' )
              end if
            end if
        case ( nf_newx ) ! ................................  NEWX  .....
        ! Set X = X + DX
        !     AJ%AXMAX = MAXVAL(ABS(X)),
        !     AJ%BIG = ANY ( DX > 10.0 * epsilon(X) * X )
          !{Account for column scaling.  We solved for $\mathbf{\Sigma^{-1}
          ! \delta x}$ above, so multiply by $\Sigma$ (which is our
          ! variable {\tt columnScaleVector}):
          if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
            & numJ=numJ, nwt_flag=nwt_flag, jacobian_rows=jacobian_rows, &
            & jacobian_cols=jacobian_cols )
          call copyVector ( v(dxUnScaled), v(dx) ) ! dxUnscaled = dx
          if ( columnScaling /= l_none ) then
            ! dxUnScaled = dxUnScaled # columnScaleVector:
            call multiply ( v(dxUnScaled), v(columnScaleVector) )
          end if
          ! Shorten dxUnscaled if necessary to stay within the bounds
          mu = 1.0
          if ( got(f_lowBound) ) &
            & call boundMove ( mu, lowBound, v(x), v(dxUnscaled), 'low', muMin )
          if ( got(f_highBound) ) &
            & call boundMove ( mu, highBound, v(x), v(dxUnscaled), 'high', muMin )
          if ( mu < 1.0_rv ) call scaleVector ( v(dxUnscaled), mu )
          call addToVector ( v(x), v(dxUnScaled) ) ! x = x + dxUnScaled
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dxUnScaled), name='dX Unscaled' )
            if ( index(switches,'xvec') /= 0 ) &
              & call dump ( v(x), name='New X' )
          aj%axmax = 0.0
          aj%big = .false.
          do j = 1, size(v(x)%quantities)
            aj%axmax = max(aj%axmax, maxval(abs(v(x)%quantities(j)%values)))
            if ( any( abs(v(dx)%quantities(j)%values) > &
              & 10.0 * epsilon(aj%axmax) * abs(v(x)%quantities(j)%values) ) ) &
              & aj%big = .true.
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
          if ( .not. aj%starting ) aj%dxdxl = &
            & aj%gfac * ( v(dx) .dot. v(bestGradient) )
          ! dx = aj%gfac * "Best Gradient":
          call scaleVector ( v(bestGradient), aj%gfac, v(dx) )
            if ( index(switches,'dvec') /= 0 ) &
              & call dump ( v(dx), name='Gradient move from best X' )
            if ( index(switches,'sca') /= 0 ) then
              call output ( ' aj%gfac = ' )
              call output ( aj%gfac, format='(1pe14.7)' )
              call output ( ' |DX| = ' )
              call output ( aj%gradnb * aj%gfac, advance='yes' )
            end if
        case ( nf_best ) ! ................................  BEST  .....
        ! Set "Best X" = X, "Best Gradient" = Gradient
          foundBetterState = .true.
          bestAJ = aj
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
          if ( .not. aj%starting ) aj%dxdxl = v(dx) .dot. v(candidateDX)
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
          ! Convergence to desired solution.  Do whatever you want to
          ! with the solution.
          if ( nwt_flag == nf_tolx_best ) then
            call copyVector ( v(x), v(bestX) )
            aj = bestAJ
          end if
          if ( .not. got(f_apriori) .or. .not. diagonal ) then
              if ( index(switches,'nwt') /= 0 ) then
                call output ( &
                  & 'Newton iteration terminated because of convergence at ' )
                call time_now ( t3 )
                call output ( t3-t0, advance='yes' )
              end if
            nwt_flag = nf_getJ
            if ( ( got(f_outputCovariance) .or. got(f_outputSD) .or. &
              &  got(f_average) ) .and. nwt_flag == nf_tolx_best ) cycle
            exit
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
          if ( index(switches,'NDB') /= 0 ) call nwtdb ( aj, width=9 )
        prev_nwt_flag = nwt_flag
      end do ! Newton iteration

        if ( index(switches,'NDB') /= 0 ) then
          call nwtdb ( aj, width=9 )
        else if ( index(switches,'Ndb') /= 0 ) then
          if ( index(switches,'sca') /= 0 ) then
            call nwtdb ( width=9 )
          else
            call nwtdb ( aj, width=9 )
          end if
        end if

      if ( got(f_diagnostics) ) call FillDiagVec ( diagnostics, aj, &
        & numJ=numJ, nwt_flag=nwt_flag, jacobian_rows=jacobian_rows, &
        & jacobian_cols=jacobian_cols )

      ! Compute the covariance of the solution
      if ( foundBetterState .and. ( got(f_outputCovariance) .or. got(f_outputSD) .or. &
        &  got(f_average) ) ) then
        if ( nwt_flag /= nf_getJ ) then
          print *, 'BUG in Retrieval module -- should need to getJ to quit'
          stop
        end if
          if ( index(switches,'cov') /= 0 ) then
            call output ( 'Begin covariance calculation at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
            call output ( 'Counted ' )
            call output ( jacobian_rows )
            call output ( ' rows and ' )
            call output ( jacobian_cols )
            call output ( ' in the Jacobian matrix at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
          end if
          call time_now ( t1 )
          call createEmptyMatrix ( temp, 0, state, state )
          call invertCholesky ( factored, temp ) ! U^{-1}
          if ( index(switches,'cov') /= 0 ) then
            call output ( &
              & 'Inverted the Cholesky factor of the normal equations at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
          end if
          preserveMatrixName = outputCovariance%m%name
          call multiplyMatrix_XY_T ( temp, temp, outputCovariance%m, &
            & diagonalOnly = .not. any ( got ( (/ f_outputCovariance, f_average/) ) ) ) ! U^{-1} U^{-T}
          if ( index(switches,'cov') /= 0 ) then
            call output ( 'Computed ' )
            if ( .not. any ( got ( (/ f_outputCovariance, f_average /) ) ) ) &
              & call output ( 'diagonal blocks of ' )
            call output ( 'U^{-1} U^{-T} at ' )
            call time_now ( t3 )
            call output ( t3-t0, advance='yes' )
          end if
          call destroyMatrix ( temp )
          call add_to_retrieval_timing( 'cholesky_invert', t1 )
          ! Scale the covariance
          if ( columnScaling /= l_none ) then
            call columnScale ( outputCovariance%m, v(columnScaleVector) )
            call rowScale ( v(columnScaleVector), outputCovariance%m )
            if ( index(switches,'cov') /= 0 ) then
              call output ( 'Scaled the Covariance matrix at ' )
              call time_now ( t3 )
              call output ( t3-t0, advance='yes' )
            end if
          end if
          outputCovariance%m%name = preserveMatrixName
          if ( associated(outputSD) ) &
            & call GetDiagonal ( outputCovariance%m, outputSD, squareRoot=.true. )
        end if

      ! Compute the averaging kernel
      if ( foundBetterState .and. got(f_average) ) then
        preserveMatrixName = outputAverage%name
        if ( extendedAverage ) then
          outputAverage = outputCovariance%m .tx. kTkStar
        else
          ! Make sure kTk is symmetrical 
          !  (outputCovariance is by virtue of its creation method as U^T U)
          Call ReflectMatrix ( kTk%m )
          outputAverage = outputCovariance%m .tx. kTk%m
        end if
        outputAverage%name = preserveMatrixName
          if ( index(switches,'cov') /= 0 ) call output ( &
            & 'Computed the Averaging Kernel from the Covariance', advance='yes' )
      end if

      call copyVector ( state, v(x) )
        if ( index(switches,'svec') /= 0 ) &
          & call dump ( state, name='Final state' )
      ! Clean up the temporaries, so we don't have a memory leak.
      if ( got(f_fuzz) ) call destroyVectorInfo ( fuzzState )
      call destroyMatrix ( normalEquations%m )
      call destroyMatrix ( kTkSep%m )
      call destroyMatrix ( kTkStar )
      call destroyMatrix ( factored%m )
      call destroyVectorInfo ( q )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
        call add_to_retrieval_timing( 'newton_solver', t1 )
    end subroutine NewtonianSolver

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( "Timing for Retrieve = " )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime


  end subroutine Retrieve

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module RetrievalModule

! $Log$
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
