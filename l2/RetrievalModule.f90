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
  use Declaration_Table, only: Num_Value, Range, Str_Range, Str_Value
  use Init_Tables_Module, only: f_apriori, f_channels, f_criteria, &
    & f_covariance, f_fwdModelIn, f_fwdModelOut, f_jacobian, f_maxIterations, &
    & f_measurements, f_method, f_outputCovariance, f_quantity, f_state, &
    & f_test, f_toleranceA, f_toleranceF, f_toleranceR, f_weight, &
    & field_first, field_indices, field_last, l_newtonian, s_subset, &
    & s_retrieve, s_vector
  use MatrixModule_1, only: Matrix_Database_T
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named, N_spec_args
  use VectorsModule, only: CopyVector, CreateMask, DestroyVectorMask, &
    & GetVectorQuantityIndexByName, operator(.DOT.), ScaleVector, SetMask, &
    & Vector_T
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
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
                                        ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), intent(in) :: MatrixDatabase

    ! Local variables
    type(nwt_T) :: AJ                   ! "About the Jacobian", see NWT.
    integer :: Apriori                  ! Index in VectorDatabase
    type(vector_T) :: BestGradient      ! for NWT
    type(vector_T) :: BestX             ! for NWT
    type(vector_T) :: CandidateDX       ! for NWT
    integer :: Channels                 ! Index in tree of "channels" field
                                        ! of subset specification, or zero
    integer :: Covariance               ! Index in MatrixDatabase of the
                                        ! covariance of Apriori
    integer :: Criteria                 ! Index in tree of "criteria" field
                                        ! of subset specification, or zero
    logical :: Diagonal                 ! "outputCovariance is a vector"
    type(vector_T) :: DX                ! for NWT
    integer :: Error
    type(vector_T) :: F                 ! for NWT -- Model - Measurements
    integer :: Field                    ! Field index -- f_something
    integer :: FwdModelIn, FwdModelOut  ! Indices in VectorDatabase
    logical :: Got(field_first:field_last)   ! "Got this field already"
    type(vector_T) :: Gradient          ! for NWT
    integer :: I, J, K                  ! Subscripts and loop inductors
    integer :: Iter                     ! Iteration number
    integer :: Jacobian                 ! Index in MatrixDatabase
    integer :: Key                      ! Index of an n_spec_args.  Either
                                        ! a son or grandson of root.
    integer :: MaxIterations            ! Maximum number of iterations of
                                        ! Newtonian method
    integer :: Measurements             ! Index in VectorDatabase
    integer :: Method                   ! Method to use for inversion, currently
                                        ! only l_Newtonian.
    integer :: Name                     ! Either 0 or, if node_id(son) ==
                                        ! n_named, subtree(1,son)
    integer :: NWT_Flag                 ! Signal from NWT, q.v., indicating
                                        ! the action to take.
    integer :: NWT_Opt(10)              ! Options for NWT, q.v.
    real(rk) :: NWT_Xopt(10)            ! Real parameters for NWT options, q.v.
    integer :: OutputCovariance         ! Index in VectorDatabase if Diagonal,
                                        ! else index in MatrixDatabase
    integer :: Quantity                 ! Index in tree of "quantity" field
                                        ! of subset specification, or zero
    integer :: QuantityIndex            ! Index within vector of a quantity
    integer :: Son                      ! Of Root or Key
    integer :: Spec                     ! s_subset or s_retrieve
    integer :: State                    ! Index in VectorDatabase
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
    double precision :: Value(2)        ! Value returned by EXPR
    integer :: vectorIndex              ! Index in VectorDatabase
    integer :: Weight                   ! Index in VectorDatabase
    type(vector_T), pointer :: X        ! for NWT

    ! Error message codes
    integer, parameter :: NoField = 1             ! A required field is missing
    integer, parameter :: NotRange = noField + 1  ! A field is a range
    integer, parameter :: Twice = notRange + 1    ! A field appears twice

    error = 0

    do i = 2, nsons(root) - 1           ! skip names at begin/end of section
      son = subtree(i, root)
      if ( node_id(son) == n_named ) then
        name = subtree(1, son)
        key = subtree(2, son)
      else
        name = 0
        key = son
      end if

      ! Key now indexes an n_spec_args vertex.  See "Configuration file
      ! parser users' guide" for pictures of the trees being analyzed.

      got = .false.
      spec = decoration(subtree(1,decoration(subtree(1,key))))
      select case ( spec )
      case ( s_subset )
        do j = 2, nsons(key) ! fields of the "subset" specification
          son = subtree(j, key)
          field = decoration(subtree(1,decoration(subtree(1,son))))
          if ( got(field) ) call announceError ( twice )
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
        maxIterations = defaultMaxIterations
        method = defaultMethod
        toleranceA = defaultToleranceA
        toleranceF = defaultToleranceF
        toleranceR = defaultToleranceR
        do j = 2, nsons(key) ! fields of the "retrieve" specification
          son = subtree(j, key)
          field = decoration(subtree(1,decoration(subtree(1,son))))
          if ( got(field) ) call announceError ( twice )
          got(field) = .true.
          select case ( field )
          case ( f_apriori )
            apriori = decoration(decoration(subtree(2,son)))
          case ( f_covariance ) ! of a priori
            covariance = decoration(decoration(subtree(2,son)))
          case ( f_fwdModelIn )
            fwdModelIn = decoration(decoration(subtree(2,son)))
          case ( f_fwdModelOut )
            fwdModelOut = decoration(decoration(subtree(2,son)))
          case ( f_jacobian )
            jacobian = decoration(decoration(subtree(2,son)))
          case ( f_maxIterations )
            call expr ( subtree(2,son), units, value, type )
            if ( type /= num_value ) call announceError ( notRange )
            maxIterations = value(1)
          case ( f_measurements )
            measurements = decoration(decoration(subtree(2,son)))
          case ( f_method )
            method = decoration(subtree(2,son))
          case ( f_outputCovariance )
            outputCovariance = decoration(decoration(subtree(2,son)))
            diagonal = decoration(subtree(1,decoration( &
                       & subtree(1,decoration(subtree(2,son)))))) == s_vector
          case ( f_state )
            state = decoration(decoration(subtree(2,son)))
          case ( f_toleranceA, f_toleranceF, f_toleranceR )
            call expr ( subtree(2,son), units, value, type )
            if ( type /= num_value ) call announceError ( notRange )
            select case ( field )
            case ( f_toleranceA )
              toleranceA = value(1)
            case ( f_toleranceF )
              toleranceF = value(1)
            case ( f_toleranceR )
              toleranceR = value(1)
            end select
          case ( f_weight )
            weight = decoration(decoration(subtree(2,son)))
          case default
            ! Shouldn't get here if the type checker worked
          end select
        end do ! j = 2, nsons(key)
        if ( error == 0 ) then
          ! Check that we have all of the necessary fields
          son = key ! in case it's needed for error messages -- see announceError
          if ( .not. got(f_fwdModelIn) ) call announceError ( noField, f_fwdModelIn )
          if ( .not. got(f_jacobian) ) call announceError ( noField, f_jacobian )
          if ( .not. got(f_measurements) ) call announceError ( noField, f_measurements )
          if ( .not. got(f_state) ) call announceError ( noField, f_state )
          ! Do the retrieval
          select case ( method )
          case ( l_newtonian )
            ! Set options for NWT
            nwt_opt(1:7) = (/  15, 1,      17, 2,      18, 3,      0 /)
            nwt_xopt(1:3) = (/ toleranceF, toleranceA, toleranceR /)
            call nwt ( nwt_flag, nwt_xopt, nwt_opt )
            ! Create extra vectors.  Altogether, we need F, X, "Best X", DX,
            ! "Candidate DX" Gradient and "Best Gradient".
            x => vectorDatabase(state)
            call cloneVector ( f, vectorDatabase(measurements) )
            call cloneVector ( bestGradient, x )
            call cloneVector ( bestX, x )
            call cloneVector ( candidateDX, x )
            call cloneVector ( DX, x )
            call cloneVector ( gradient, x )
            ! Create initial guess for X.  Use Apriori if we have it, else
            ! zero is probably as good a guess as anything.
            if ( got(f_apriori) ) then
              call copyVector ( x, vectorDatabase(apriori) )
            else
              do j = 1, size(x%quantities)
                x%quantities(j)%values = 0.0_rk
              end do ! j
            end if
            iter = 0
            do
              call nwta ( nwt_flag, aj )
              select case ( nwt_flag ) ! >0 means "Done", <0 means "Continue"
              case ( nf_evalf )
              ! IF ( too many function values ) EXIT
                if ( iter > maxIterations ) exit
              ! Compute f(x)
                call subtractFromVector ( f, vectorDatabase(measurements) )
                aj%fnorm = sqrt(f .dot. f) ! L2Norm(f)
                if ( aj%fnorm < toleranceF ) exit
              ! Compute the Jacobian matrix J if you feel like it
              case ( nf_evalj )
              ! Compute the Jacobian matrix J if you didn't do it when NWT_FLAG
              ! was NF_EVALF:
              !   J(K,L) = Partial of F(K) / W.R.T. X(L), K = 1, NF, L = 1, NX
              ! Triangularize J, and compute (negative of the) gradient =
              ! -(Jacobian)**T * F
              ! Set
              !   AJ%DIAG = element on diagonal with smallest absolute value,
              !           after factoring,
              !   AJ%AJN = maximum L1 norm of column in upper triangle
              !           after factoring,
              !   AJ%FNMIN = L2 norm of residual not in the column space of
              !           the Jacobian,
                aj%gradn = sqrt(gradient .dot. gradient) ! L2Norm(gradient)
              case ( nf_solve )
              ! Apply Marquardt stabilization with parameter = AJ%SQ, and
              ! solve for "candidate DX" = -(Jacobian)**(-1) * F
              ! Set AJ%FNMIN as for NWT_FLAG = NF_EVALJ, but taking account
              ! of Levenberg-Marquardt stabilization.
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
              !     dx = aj%gfac * "Best Gradient":
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
              ! dx = aj%cait * candidateDX
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
!             IF ( you want to return to a previous best X ) NWT_FLAG = 0
              iter = iter + 1
            end do ! Newton iteration
            ! Clean up the temporary vectors, so we don't have a memory leak
            call destroyVectorInfo ( bestGradient )
            call destroyVectorInfo ( bestX )
            call destroyVectorInfo ( candidateDX )
            call destroyVectorInfo ( dx )
            call destroyVectorInfo ( f )
            call destroyVectorInfo ( gradient )
          end select ! method
        end if
        ! Clear the masks of every vector
        do j = 1, size(vectorDatabase)
          call destroyVectorMask ( vectorDatabase(i) )
        end do
      end select

    end do ! i = 2, nsons(root) - 1

  contains
    subroutine AnnounceError ( Code, FieldIndex )
      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in), optional :: FieldIndex ! f_...
      integer :: MyField, Source

      error = max(error,1)
      myField = field
      if ( present(fieldIndex) ) myField = fieldIndex
      select case ( code )
      case ( noField, notRange, twice )
        source = source_ref ( son )
        call output ( 'At line '  )
        call output ( mod(source,256) )
        call output ( ', column ' )
        call output ( source/256 )
        call output ( ': the field ' )
        call display_string ( field_indices(myField) )
        select case ( code )
        case ( noField )
          call output ( ' is not specified', advance='yes' )
        case ( notRange )
          call output ( ' shall not be a range', advance='yes' )
        case ( twice )
          call output ( ' shall not appear twice', advance='yes' )
        end select
      end select
    end subroutine AnnounceError
  end subroutine Retrieve

end module RetrievalModule

! $Log$
! Revision 2.1  2001/01/10 21:04:13  vsnyder
! Initial (incomplete) submission
!
