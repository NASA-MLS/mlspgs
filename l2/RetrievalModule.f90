! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module RetrievalModule
!=============================================================================

! This module inverts the radiative transfer equation, to solve for
! atmospheric parameters, given radiance measurements.

! This module and ones it calls consume most of the cycles.

  use DNWT_Module, only: NF_EVALF, NF_EVALJ, NF_SOLVE, NF_NEWX, &
    & NF_GMOVE, NF_BEST, NF_AITKEN, NF_DX, NF_DX_AITKEN, NF_TOLX, &
    & NF_TOLX_BEST, NF_TOLF, NF_TOO_SMALL, NF_FANDJ, NWT, NWT_T, NWTA, RK
  use Expr_M, only: Expr
  use Declaration_Table, only: Num_Value
  use Init_Tables_Module, only: f_apriori, f_aprioriScale, f_channels, &
    & f_criteria, f_columnScale, f_covariance, f_diagonal, f_fwdModelIn, &
    & f_fwdModelOut, f_jacobian, f_maxIterations, f_measurements, f_method, &
    & f_outputCovariance, f_quantity, f_state, f_test, f_toleranceA, &
    & f_toleranceF, f_toleranceR, f_weight, field_first, field_indices, &
    & field_last, l_apriori, l_covariance, l_newtonian, l_none, l_norm, &
    & l_true, s_matrix, s_subset, s_retrieve
  use MatrixModule_1, only: AddToMatrixDatabase, CholeskyFactor, ClearMatrix, &
    & ColumnScale, CopyMatrixValue, CreateEmptyMatrix, DestroyMatrix, &
    & FillExtraCol, FillExtraRow, FormNormalEquations => NormalEquations, &
    & GetDiagonal, GetFromMatrixDatabase, GetVectorFromColumn, Matrix_T, &
    & Matrix_Database_T, Matrix_Cholesky_T, Matrix_SPD_T, MaxAbsVal, MinDiag, &
    & MultiplyMatrixVector, RowScale, ScaleMatrix, SolveCholesky, &
    & UpdateDiagonal
  use MLSCommon, only: R8
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named, N_set_one
  use VectorsModule, only: AddToVector, CloneVector, CopyVector, CreateMask, &
    & DestroyVectorInfo, DestroyVectorMask, DestroyVectorValue, &
    & GetVectorQuantityIndexByName, MultiplyVectors, operator(.DOT.), &
    & ScaleVector, SetMask, SubtractFromVector, Vector_T

  implicit NONE
  private
  public :: RETRIEVE

  integer, parameter, private :: DefaultMaxIterations = 5
  integer, parameter, private :: DefaultMethod = l_newtonian
  double precision, parameter, private :: DefaultToleranceA = 1.0d-6 ! for NWT
  double precision, parameter, private :: DefaultToleranceF = 1.0d-6 ! for NWT
  double precision, parameter, private :: DefaultToleranceR = 1.0d-6 ! for NWT

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------  Retrieve  -----
  subroutine Retrieve ( Root, VectorDatabase, MatrixDatabase )
  ! Process the "Retrieve" section of the L2 Configuration File.
  ! The "Retrieve" section can have Matrix, Subset or Retrieve specifications.
    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
                                        ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase

    ! Local variables:
    type(nwt_T) :: AJ                   ! "About the Jacobian", see NWT.
    type(vector_T), pointer :: Apriori  ! A priori estimate of state
    real(r8) :: aprioriNorm             ! apriori .dot. apriori, for n+1, n+1
                                        ! element of augmented normal equations
    real(r8) :: aprioriScale            ! Weight for apriori, default 1.0
    type(vector_T) :: BestGradient      ! for NWT
    type(vector_T) :: BestX             ! for NWT
    type(vector_T) :: CandidateDX       ! for NWT
    integer :: Channels                 ! Index in tree of "channels" field
                                        ! of subset specification, or zero
    type(vector_T) :: ColumnScaleVector ! For column scaling by column norms
    integer :: ColumnScaling            ! one of l_apriori, l_covariance,
                                        ! l_none or l_norm
    type(matrix_SPD_T), pointer :: Covariance     ! covariance**(-1) of Apriori
    type(vector_T) :: CovarianceXApriori ! Covariance \times Apriori
    integer :: Criteria                 ! Index in tree of "criteria" field
                                        ! of subset specification, or zero
    type(vector_T) :: DX                ! for NWT
    logical :: Diagonal                 ! "We only want the diagonal of the
                                        ! a posteriori covariance matrix"
    integer :: Error
    type(vector_T) :: F                 ! for NWT -- Model - Measurements
    type(matrix_Cholesky_T) :: Factored ! Cholesky-factored normal equations
    integer :: Field                    ! Field index -- f_something
    type(vector_T), pointer :: FwdModelIn, FwdModelOut
    logical :: Got(field_first:field_last)   ! "Got this field already"
    type(vector_T) :: Gradient          ! for NWT
    integer :: I, J, K                  ! Subscripts and loop inductors
    integer :: Iter                     ! Iteration number
    integer :: IxCovariance             ! Index in tree of outputCovariance
    integer :: IxJacobian               ! Index in tree of jacobian matrix
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Key                      ! Index of an n_spec_args.  Either
                                        ! a son or grandson of root.
    integer :: MaxIterations            ! Maximum number of iterations of
                                        ! Newtonian method
    type(vector_T), pointer :: Measurements  ! The measurements vector
    integer :: Method                   ! Method to use for inversion, currently
                                        ! only l_Newtonian.
    type(matrix_SPD_T), target :: myCovariance    ! for OutputCovariance to point at
    type(matrix_T), target :: myJacobian          ! for Jacobian to point at
    integer :: N                        ! 1 + Number of columns in the Jacobian
    integer :: Name                     ! Either 0 or, if node_id(son) ==
                                        ! n_named, subtree(1,son)
    type(matrix_SPD_T) :: NormalEquations         ! Jacobian**T * Jacobian
    integer :: NWT_Flag                 ! Signal from NWT, q.v., indicating
                                        ! the action to take.
    integer :: NWT_Opt(10)              ! Options for NWT, q.v.
    real(rk) :: NWT_Xopt(10)            ! Real parameters for NWT options, q.v.
    type(matrix_SPD_T), pointer :: OutputCovariance   ! Covariance of the sol'n
    integer :: Quantity                 ! Index in tree of "quantity" field
                                        ! of subset specification, or zero
    integer :: QuantityIndex            ! Index within vector of a quantity
    integer :: RowBlock                 ! Which block of rows is the forward
                                        ! model filling?
    integer :: Son                      ! Of Root or Key
    integer :: Spec                     ! s_matrix, s_subset or s_retrieve
    type(vector_T), pointer :: State    ! The state vector
    integer :: Test                     ! Index in tree of "test" field
                                        ! of subset specification, or zero
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
    integer :: vectorIndex              ! Index in VectorDatabase
    type(vector_T), pointer :: Weight   ! Scaling vector for rows
    type(vector_T), pointer :: X        ! for NWT

    ! Error message codes
    integer, parameter :: AprioriAndCovar = 1     ! Only one of apriori and
                                                  ! covariance supplied
    integer, parameter :: Inconsistent = AprioriAndCovar + 1 ! Inconsistent fields
    integer, parameter :: NoField = Inconsistent + 1   ! Required field missing
    integer, parameter :: NoFields = NoField + 1  ! No fields are allowed
    integer, parameter :: NotExtra = noFields + 1 ! No "extra" row and/or column
    integer, parameter :: NotRange = notExtra + 1 ! A field is a range
    integer, parameter :: NotSPD = notRange + 1   ! Not symmetric pos. definite
    integer, parameter :: Twice = notSPD + 1      ! A field appears twice

    error = 0
    nullify ( apriori, covariance, fwdModelIn, fwdModelOut )
    nullify ( measurements, state, weight, x )

    do i = 2, nsons(root) - 1           ! skip names at begin/end of section
      son = subtree(i, root)
      if ( node_id(son) == n_named ) then
        name = subtree(1, son)
        key = subtree(2, son)
      else
        name = 0
        key = son
      end if

      ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
      ! parser users' guide" for pictures of the trees being analyzed.

      got = .false.
      spec = decoration(subtree(1,decoration(subtree(1,key))))
      select case ( spec )
      case ( s_matrix )
        if ( nsons(key) /= 1 ) call announceError ( noFields )
        call destroyMatrix( matrixDatabase(decoration(key)) ) ! avoids a memory leak
        call decorate ( key, 0 )
      case ( s_subset )
        do j = 2, nsons(key) ! fields of the "subset" specification
          son = subtree(j, key)
          field = decoration(subtree(1,decoration(subtree(1,son))))
          if ( got(field) ) call announceError ( twice, field )
          got(field) = .true.
          select case ( field )
          case ( f_channels )
            channels = son
          case ( f_criteria )
            criteria = son
          case ( f_quantity )
            quantity = son
          case ( f_test )
            test = son
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)
        if ( error == 0 ) then
          ! Compute the mask for the vector
          son = key ! in case it's needed for error messages -- see announceError
          if ( .not. got(f_criteria) ) call announceError ( noField, f_criteria )
          if ( .not. got(f_quantity) ) call announceError ( noField, f_quantity )
          if ( .not. got(f_test) ) call announceError ( noField, f_test )
        end if
        if ( error == 0 ) then
          vectorIndex = decoration(decoration(subtree(1,quantity)))
          quantityIndex = getVectorQuantityIndexByName &
            & ( vectorDatabase(vectorIndex), sub_rosa(subtree(2,quantity)) )
          if ( .not. associated(vectorDatabase(vectorIndex)% &
            & quantities(quantityIndex)%mask) ) then ! create and fill with 1's
            call createMask ( vectorDatabase(vectorIndex)% &
              & quantities(quantityIndex) )
            do j = 1, size(vectorDatabase(vectorIndex)% &
              & quantities(quantityIndex)%mask, 2)
              call setMask ( vectorDatabase(vectorIndex)% &
                & quantities(quantityIndex)%mask(:,j) )
            end do
          end if ! mask exists
          !??? Fill the mask according to the specified criteria ???
        end if ! error == 0
      case ( s_retrieve )
        aprioriScale = 1.0
        columnScaling = l_none
        maxIterations = defaultMaxIterations
        method = defaultMethod
        toleranceA = defaultToleranceA
        toleranceF = defaultToleranceF
        toleranceR = defaultToleranceR
        do j = 2, nsons(key) ! fields of the "retrieve" specification
          son = subtree(j, key)
          field = decoration(subtree(1,decoration(subtree(1,son))))
          if ( got(field) ) call announceError ( twice, field )
          got(field) = .true.
          select case ( field )
          case ( f_apriori )
            apriori => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_columnScale )
            columnScaling = decoration(subtree(2,son))
          case ( f_covariance ) ! of a priori
            call getFromMatrixDatabase ( &
              & matrixDatabase(decoration(decoration(subtree(2,son)))), &
              & covariance)
            if ( .not. associated(covariance) ) then
              call announceError ( notSPD, field )
            else
              if ( .not. covariance%m%row%extra .or. &
                &  .not. covariance%m%col%extra ) &
                & call announceError ( notExtra, field )
            end if
          case ( f_diagonal )
            son = subtree(2,son)
            if ( node_id(son) == n_set_one ) then
              diagonal = .true.
            else
              diagonal = decoration(subtree(2,son)) == l_true
            end if
          case ( f_fwdModelIn )
            fwdModelIn => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_fwdModelOut )
            fwdModelOut => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_jacobian )
            ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
          case ( f_maxIterations )
            call expr ( subtree(2,son), units, value, type )
            if ( type /= num_value ) call announceError ( notRange, field )
            maxIterations = value(1)
          case ( f_measurements )
            measurements => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_method )
            method = decoration(subtree(2,son))
          case ( f_outputCovariance )
            ixCovariance = decoration(subtree(2,son)) ! outCov: matrix vertex
          case ( f_state )
            state => vectorDatabase(decoration(decoration(subtree(2,son))))
          case ( f_aprioriScale, f_toleranceA, f_toleranceF, f_toleranceR )
            call expr ( subtree(2,son), units, value, type )
            if ( type /= num_value ) call announceError ( notRange, field )
            select case ( field )
            case ( f_aprioriScale )
              aprioriScale = value(1)
            case ( f_toleranceA )
              toleranceA = value(1)
            case ( f_toleranceF )
              toleranceF = value(1)
            case ( f_toleranceR )
              toleranceR = value(1)
            end select
          case ( f_weight )
            weight => vectorDatabase(decoration(decoration(subtree(2,son))))
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)

        ! Check that we have all of the necessary fields
        son = key ! in case it's needed for error messages -- see announceError
        if ( .not. got(f_fwdModelIn) ) call announceError ( noField, f_fwdModelIn )
        if ( .not. got(f_measurements) ) call announceError ( noField, f_measurements )
        if ( .not. got(f_state) ) call announceError ( noField, f_state )
        if ( got(f_apriori) .neqv. got(f_covariance) ) &
          & call announceError ( aprioriAndCovar )
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
          if ( got(f_weight) ) then
            if ( weight%template%id /= measurements%template%id ) &
              & call announceError ( inconsistent, f_weight, f_measurements )
          end if
        end if
        if ( error == 0 ) then

          ! Create the matrix for the Cholesky factor of the normal equations
          call createEmptyMatrix ( factored%m, 0, state, state, &
            &                      extra_col=.true., extra_row=.true. )
          ! Create the Jacobian matrix
          if ( got(f_jacobian) ) then
            k = decoration(ixJacobian)
            if ( k == 0 ) then
              call createEmptyMatrix ( myJacobian, sub_rosa(k), measurements, &
                &                      state, extra_col=.true. )
              k = addToMatrixDatabase( matrixDatabase, myJacobian )
              call decorate ( ixJacobian, k )
            end if
            call getFromMatrixDatabase ( matrixDatabase(k), jacobian )
            if ( jacobian%row%vec%template%id /= measurements%template%id ) &
              & call announceError ( inconsistent, f_jacobian, f_measurements )
            if ( jacobian%col%vec%template%id /= state%template%id ) &
              & call announceError ( inconsistent, f_jacobian, f_state )
            if ( .not. jacobian%col%extra ) &
              & call announceError ( notExtra, f_jacobian )
          else
            jacobian => myJacobian
            call createEmptyMatrix ( jacobian, 0, measurements, state, &
              &                      extra_col=.true. )
          end if
          n = sum(jacobian%col%nelts)   ! 1 + Number of columns in the Jacobian
          ! Create the normal equations matrix
          call createEmptyMatrix ( normalEquations%m, 0, state, state, &
            &                      extra_col=.true., extra_row=.true. )
          ! Create the output covariance matrix
          if ( got(f_outputCovariance) ) then
            k = decoration(ixCovariance)
            if ( k == 0 ) then
              call createEmptyMatrix ( myCovariance%m, sub_rosa(k), &
                &                      state, state, &
                &                      extra_col=.true., extra_row=.true. )
              k = addToMatrixDatabase( matrixDatabase, myCovariance )
              call decorate ( ixCovariance, k )
            end if
            call getFromMatrixDatabase ( matrixDatabase(k), outputCovariance )
            if ( .not. associated(outputCovariance) ) then
              call announceError ( notSPD, f_outputCovariance )
            else
              if ( .not. outputCovariance%m%row%extra .or. &
                &  .not. outputCovariance%m%col%extra ) &
                & call announceError ( notExtra, f_outputCovariance )
              if ( outputCovariance%m%row%vec%template%id /= state%template%id &
                & .or. &
                &  outputCovariance%m%col%vec%template%id /= state%template%id ) &
                & call announceError ( inconsistent, f_outputCovariance, f_state )
            end if
          else
            outputCovariance => myCovariance
            call createEmptyMatrix ( myCovariance%m, 0, state, state, &
              &                      extra_col = .true., extra_row = .true. )
          end if
        end if
        if ( error == 0 ) then

          ! Do the retrieval
          select case ( method )
          case ( l_newtonian )
            ! Set options for NWT
            nwt_opt(1:7) = (/  15, 1,      17, 2,      18, 3,      0 /)
            nwt_xopt(1:3) = (/ toleranceF, toleranceA, toleranceR /)
            call nwt ( nwt_flag, nwt_xopt, nwt_opt )
            ! Create extra vectors.  Altogether, we need F, X, "Best X", DX,
            ! "Candidate DX" Gradient and "Best Gradient".
            x => state
            call cloneVector ( f, measurements )
            call cloneVector ( bestGradient, x )
            call cloneVector ( bestX, x )
            call cloneVector ( candidateDX, x )
            call cloneVector ( DX, x )
            call cloneVector ( gradient, x )
            ! Create initial guess for X.  Use Apriori if we have it, else
            ! zero is probably as good a guess as anything.
            if ( got(f_apriori) ) then
              call copyVector ( x, apriori ) ! x := apriori
            else
              do j = 1, size(x%quantities)
                x%quantities(j)%values = 0.0_rk
              end do ! j
            end if
            if ( got(f_apriori) ) then
              aprioriNorm = apriori .dot. apriori ! norm**2
              call multiplyMatrixVector ( covariance, apriori, &
                & covarianceXApriori )
            end if
            iter = 0
            do
              call nwta ( nwt_flag, aj )
              select case ( nwt_flag ) ! >0 means "Done", <0 means "Continue"
              case ( nf_evalf )
                if ( iter > maxIterations ) exit
                ! Compute f(x)
                call forwardModelRadiances ( fwdModelIn, f )
                call subtractFromVector ( f, measurements )
                aj%fnorm = sqrt(f .dot. f)
                call destroyVectorValue ( f )  ! free the space
                if ( aj%fnorm < toleranceF ) exit
              case ( nf_evalj )
                ! Compute the Jacobian matrix J if you didn't do it when
                ! NWT_FLAG was NF_EVALF:
                !   J(K,L) = Partial of F(K) / W.R.T. X(L), K = 1, NF, L = 1, NX
                update = got(f_apriori)
                if ( update ) then ! start normal equations with apriori
                  call copyMatrixValue ( normalEquations%m, covariance%m )
                  call fillExtraCol ( normalEquations%m, covarianceXApriori )
                  call fillExtraRow ( normalEquations%m, covarianceXApriori )
                  k = normalEquations%m%row%nb
                  normalEquations%m%block(k,k)%values(1,1) = aprioriNorm
                  if ( got(f_aprioriScale) ) &
                    & call scaleMatrix ( normalEquations%m, aprioriScale )
                else
                  call clearMatrix ( normalEquations%m ) ! start with zero
                end if
                do rowBlock = 1, jacobian%row%nb
                  call forwardModel ( fwdModelIn, jacobian, f, rowBlock, &
                    & fwdModelOut )
                  call subtractFromVector ( f, measurements, &
                    & jacobian%row%quant(rowBlock), &
                    & jacobian%row%inst(rowBlock) )
                  call fillExtraCol ( jacobian, f, rowBlock )
                  if ( got(f_weight) ) call rowScale ( weight, jacobian )
                  call formNormalEquations ( jacobian, normalEquations, &
                    & update=update )
                  update = .true.
                  call clearMatrix ( jacobian )  ! free the space
                  call destroyVectorValue ( f )  ! free the space
                end do
                ! Compute (negative of the) gradient = -(Jacobian)**T * F.
                ! This is the RHS of the normal equations J**T * J *
                ! "Candidate DX" = -J**T * F:
                call getVectorFromColumn ( normalEquations%m, n, gradient )
                ! Column Scale J (row and column scale J^T J):
                select case ( columnScaling )
                case ( l_apriori )
                  call columnScale ( normalEquations%m, apriori )
                  call rowScale ( apriori, normalEquations%m )
                case ( l_covariance )
                  !??? Can't get here until allowed by init_tables
                case ( l_norm )
                  call getDiagonal ( normalEquations, columnScaleVector )
                  do j = 1, columnScaleVector%template%noQuantities
                    columnScaleVector%quantities(j)%values = &
                      & sqrt( columnScaleVector%quantities(j)%values )
                  end do
                  call columnScale ( normalEquations%m, columnScaleVector )
                  call rowScale ( columnScaleVector, normalEquations%m )
                end select
                ! Factor J^T J:
                factored%m%block => normalEquations%m%block ! to save space
                call choleskyFactor ( normalEquations, factored )
                aj%diag = minDiag ( factored ) ! element on diagonal with
                !       smallest absolute value, after triangularization
                aj%ajn = maxAbsVal ( factored%m ) ! maximum L1 norm of
                !       column in upper triangle after triangularization
                k = factored%m%col%nb
                ! AJ%FNMIN = L2 norm of residual, ||F + J * "Candidate DX"||
                !       (which can be gotten without solving for
                !       "Candidate DX": if -F is put as the last column of
                !       J before triangularization, either by Householder
                !       or by Cholesky factoring the normal equations,
                !       this is J(N+1,N+1)).  The (extra row, extra column)
                !       block is 1x1:
                aj%fnmin = factored%m%block(k,k)%values(1,1)
                aj%gradn = sqrt(gradient .dot. gradient) ! L2Norm(gradient)
              case ( nf_solve )
                ! Apply Marquardt stabilization with parameter = AJ%SQ:
                call updateDiagonal ( normalEquations, aj%sq**2 )
                factored%m%block => normalEquations%m%block ! to save space
                call choleskyFactor ( normalEquations, factored )
                ! Solve for "candidate DX" = -(Jacobian)**(-1) * F
                call getVectorFromColumn ( factored%m, n, candidateDX )
                call solveCholesky ( factored, candidateDX )
                ! Account for column scaling
                select case ( columnScaling )
                case ( l_apriori )
                  call multiplyVectors ( dx, apriori ) ! dx = dx # apriori
                case ( l_covariance )
                  !??? Can't get here until allowed by init_tables
                case ( l_norm )
                  call multiplyVectors ( dx, columnScaleVector ) ! dx = dx # ...
                end select
                ! Set AJ%FNMIN as for NWT_FLAG = NF_EVALJ, but taking account
                ! of Levenberg-Marquardt stabilization:
                k = factored%m%col%nb
                aj%fnmin = factored%m%block(k,k)%values(1,1)
                aj%dxn = sqrt(candidateDX .dot. candidateDX) ! L2Norm(dx)
                aj%gdx = gradient .dot. candidateDX
              case ( nf_newx )
                call addToVector ( x, dx ) ! x = x + dx
                aj%axmax = 0.0
                aj%big = .false.
                do j = 1, size(x%quantities)
                  aj%axmax = max(aj%axmax, maxval(abs(x%quantities(j)%values)))
                  aj%big = aj%big .or. any( abs(dx%quantities(j)%values) > &
                    & 10.0 * epsilon(aj%axmax) * abs(x%quantities(j)%values) )
                end do
                if ( .not. aj%starting ) aj%dxdxl = dx .dot. candidateDX
              case ( nf_gmove )
                call copyVector ( x, bestX ) ! x = bestX
                ! dx = aj%gfac * "Best Gradient":
                call scaleVector ( bestGradient, aj%gfac, dx )
              case ( nf_best )
                call copyVector ( bestX, x ) ! bestX = x
                call copyVector ( bestGradient, gradient ) ! bestGradient = gradient
              case ( nf_aitken )
                call subtractFromVector ( dx, candidateDX ) ! dx = dx - candidateDX
                aj%dxdx = dx .dot. dx
                if ( aj%dxdx > 0.0 ) aj%dxdxl = dx .dot. candidateDX
              case ( nf_dx )
                call copyVector ( dx, candidateDX ) ! dx = candidateDX
              case ( nf_dx_aitken )
                ! dx = aj%cait * candidateDX:
                call scaleVector ( candidateDX, aj%cait, dx )
              case ( nf_tolx, nf_tolx_best, nf_tolf, nf_too_small )
                ! IF ( NWT_FLAG == NF_TOO_SMALL ) THEN
                !   Take special action if requested accuracy is critical
                ! END IF
                if ( nwt_flag == nf_tolx_best ) call copyVector ( x, bestX )
                ! Convergence to desired solution.  Do whatever you want to
                ! with the solution.
                exit ! unless you have a really good reason to continue
              case ( nf_fandj )
                ! There is probably an error in the way F or J is computed.
                ! A warning has been printed by the error processor.
                ! IF ( you have confidence in F and J ) CYCLE
                ! STOP
              end select
            ! IF ( you want to return to a previous best X ) NWT_FLAG = 0
              iter = iter + 1
            end do ! Newton iteration
            ! Clean up the temporaries, so we don't have a memory leak
            call destroyVectorInfo ( bestGradient )
            call destroyVectorInfo ( bestX )
            call destroyVectorInfo ( candidateDX )
            call destroyVectorInfo ( dx )
            call destroyVectorInfo ( f )
            call destroyVectorInfo ( gradient )
            call destroyMatrix ( normalEquations%m )
          end select ! method
          !??? Make sure the jacobian and outputCovariance get destroyed
          !??? after ?what? happens?  Can we destroy the entire matrix
          !??? database at the end of each chunk?
          if ( .not. got(f_jacobian) ) call destroyMatrix ( jacobian )
          if ( .not. got(f_outputCovariance) ) &
            & call destroyMatrix ( outputCovariance%m )
        end if
        ! Clear the masks of every vector
        do j = 1, size(vectorDatabase)
          call destroyVectorMask ( vectorDatabase(i) )
        end do
      end select

    end do ! i = 2, nsons(root) - 1

  contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex, AnotherFieldIndex )
      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex, AnotherFieldIndex ! f_...
      integer :: MyField, Source

      error = max(error,1)
      myField = field
      if ( present(fieldIndex) ) myField = fieldIndex
      source = source_ref ( son )
      call output ( 'At line '  )
      call output ( mod(source,256) )
      call output ( ', column ' )
      call output ( source/256 )
      select case ( code )
      case ( aprioriAndCovar )
        call output ( 'One of ' )
        call display_string ( field_indices(f_apriori) )
        call output ( ' or ' )
        call display_string ( field_indices(f_covariance) )
        call output ( ' is supplied, but the other is not.', advance='yes' )
      case ( noFields )
        call output ( ': No fields are allowed.' )
      case ( inconsistent, noField, notExtra, notRange, notSPD, twice )
        call output ( ': the field ' )
        call display_string ( field_indices(myField) )
        select case ( code )
        case ( inconsistent )
          call output ( ' is not consistent with the ' )
          call display_string ( field_indices(anotherFieldIndex ) )
          call output ( ' field.', advance='yes' )
        case ( noField )
          call output ( ' is not specified.', advance='yes' )
        case ( notExtra )
          call output ( ' does not have an extra column for the RHS.', &
            & advance='yes' )
        case ( notRange )
          call output ( ' shall not be a range.', advance='yes' )
        case ( notSPD )
          call output ( ' is not a symmetric positive-definite matrix.', &
            & advance='yes' )
        case ( twice )
          call output ( ' shall not appear twice.', advance='yes' )
        end select
      end select
    end subroutine AnnounceError
  end subroutine Retrieve

end module RetrievalModule

! $Log$
! Revision 2.2  2001/01/26 19:01:47  vsnyder
! More nearly complete, except for forward model interface and minor things
! having to do with creating subset masks.  Look for ??? in comments.
!
! Revision 2.1  2001/01/10 21:04:13  vsnyder
! Initial (incomplete) submission
!
