! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ALGEBRA_M

! Process the ALGEBRA section

  implicit NONE
  private
  public :: Algebra

!--------------------------------- RCS Ident Info --------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------------  Algebra  -----
  subroutine Algebra ( ROOT, VectorDatabase, MatrixDatabase )
    ! The root of the Algebra section's subtree is ROOT.  All of the
    ! input and output "variables" are found in the symbol table.
    ! Anything defined here by a statement of the form A = <expr> is
    ! deleted here if A isn't already found in the symbol table, i.e.,
    ! not declared before the Algebra section.  If it is in the symbol
    ! table, it has to be the right kind of thing-o:  We won't assign
    ! a vector to a matrix, etc.

    use DECLARATION_TABLE, only: DECLARE, DECLS, EMPTY, EXPRN, EXPRN_M, &
      & EXPRN_V, GET_DECL, LABEL, NUM_VALUE, REDECLARE
    use Init_Tables_Module, only: S_Matrix, S_Vector
    use MatrixModule_1, only: AddToMatrix, AddToMatrixDatabase, AssignMatrix, &
      & CreateEmptyMatrix, CopyMatrix, &
      & CopyMatrixValue, DestroyMatrix, Dump, GetActualMatrixFromDatabase, &
      & GetFromMatrixDatabase, GetKindFromMatrixDatabase, &
      & K_Cholesky, K_Empty, K_Kronecker, K_Plain, K_SPD, Matrix_Cholesky_T, &
      & Matrix_Database_T, Matrix_Kronecker_T, Matrix_SPD_T, Matrix_T, &
      & MultiplyMatrixVectorNoT, multiplyMatrixVectorSPD_1, &
      & MultiplyMatrix_XY, ScaleMatrix, Dump_Struct, TransposeMatrix
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_M, only: Output
    use String_Table, only: Display_String
    use Toggles, only: Gen, Switches, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use TREE, only: DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE, PRINT_SUBTREE
    use TREE_TYPES ! Everything, except tree_init; remainder begin with N_
    use VectorsModule, only:  AddToVector, AddVectorToDatabase,  CopyVector, &
      & CloneVector, DestroyVectorInfo, Dump, ScaleVector, Vector_T

    integer, intent(in) :: ROOT
    type(vector_T), dimension(:), pointer :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase

    type(decls) :: DECL            ! Declaration of the LHS name
    real(r8) :: DVALUE             ! Value of expr if it's a number
    integer :: Error               ! 0 = no error occurred.
    integer :: I_SONS              ! Son indices of sons of root
    integer :: LHS                 ! Index in tree of LHS
    integer :: RHS                 ! Index in tree of root of RHS
    integer :: SON                 ! Son of root
    integer :: SPEC                ! Index of SPEC of label, e.g. S_Matrix
    integer :: String              ! String table index
    integer :: Type                ! Type of result
    integer :: Value               ! Index of Vector or Matrix result in database
    integer :: What                ! See its parameters below

    type(matrix_t), pointer :: LHSMatrix
    type(matrix_Cholesky_t) :: Matrix_C
    type(matrix_Kronecker_t) :: Matrix_K
    type(matrix_t), target :: Matrix
    type(matrix_SPD_t) :: Matrix_S

    type(vector_t) :: Vector

    ! Parameters for Announce_Error
    integer, parameter :: Ambiguous = 1                ! LHS is ambiguous
    integer, parameter :: CantInvertVector = ambiguous + 1    ! Can't invert vector
    integer, parameter :: Incompatible = cantInvertVector + 1 ! incompatible operands
    integer, parameter :: Undefined = incompatible + 1 ! undefined operand
    integer, parameter :: UnknownFunc = undefined + 1  ! unknown function
    integer, parameter :: UnknownOp = UnknownFunc + 1  ! unknown operator
    integer, parameter :: WrongNumArgs = UnknownOp + 1 ! wrong number of args

    ! Parameters for WHAT argument of EXPR
    integer, parameter :: W_Nothing = 0               ! An error occurred
    integer, parameter :: W_Number = w_nothing + 1    ! Result is DVALUE
    integer, parameter :: W_Vector = w_number + 1     ! Result is a vector
    integer, parameter :: W_Matrix = w_vector + 1     ! Result is a general matrix
    integer, parameter :: W_Matrix_C = w_matrix + 1   ! Result is Cholesky factor
    integer, parameter :: W_Matrix_K = w_matrix_C + 1 ! Result is a Kronecker matrix
    integer, parameter :: W_Matrix_S = w_matrix_K + 1 ! Result is an SPD matrix

    error = 0

    if ( toggle(gen) ) call trace_begin ( 'Algebra', root )
    do i_sons = 2, nsons(root) - 1 ! skip names at begin and end
      son = subtree(i_sons,root)
      if ( toggle(gen) ) call trace_begin ( 'Algebra loop', son )
      if ( node_id(son) /= n_equal ) then
        call AlgebraCommands ( son, VectorDatabase, MatrixDatabase )
      else
        ! Evaluate the RHS
        rhs = subtree(2,son)
        call expr ( rhs, what, dValue, vector, matrix, matrix_c, matrix_k, matrix_s )
        if ( what == w_nothing ) &
    go to 200
        ! Lookup the LHS.  It could be a label or the LHS of another expr.
        lhs = subtree(1,son)
        string = sub_rosa(lhs)
        decl = get_my_decl(lhs)
        ! If it's a vector or matrix, put it in the global database.
        ! Put the database index into the LHS's decl.  If it's a number,
        ! just put the number in the LHS decl's value field.
        if ( decl%type == empty ) then ! --------------------- Not defined yet
          select case ( what )
          case ( w_number )
            type = num_value
            value = 0
          case ( w_vector )
            type = exprn_v
            value = addVectorToDatabase ( vectorDatabase, vector )
          case ( w_matrix )
            type = exprn_m
            value = addToMatrixDatabase ( matrixDatabase, matrix )
          case ( w_matrix_c)
            type = exprn_m
            value = addToMatrixDatabase ( matrixDatabase, matrix_c )
          case ( w_matrix_k)
            type = exprn_m
            value = addToMatrixDatabase ( matrixDatabase, matrix_k )
          case ( w_matrix_s )
            type = exprn_m
            value = addToMatrixDatabase ( matrixDatabase, matrix_s )
          end select
          call declare ( string, dvalue, type, value, son )
        else ! ------------------- Update the declaration for an existing name
          value = decoration ( decl%tree )
          select case ( what )
          case ( w_number )
            if ( decl%type /= exprn .and. decl%type /= num_value ) then
              call announce_error ( son, incompatible )
              what = w_nothing
   go to 100
            end if
          case ( w_vector )
            if ( decl%type == label ) then
              spec = get_spec(decl%tree)
              if ( spec /= s_vector ) then
                call announce_error ( son, incompatible )
                what = w_nothing
   go to 100
              end if
            else if ( decl%type /= exprn_v ) then
              call announce_error ( son, incompatible )
              what = w_nothing
   go to 100
            end if
            call copyVector ( vectorDatabase(value), vector )
          case ( w_matrix, w_matrix_c, w_matrix_k, w_matrix_s )
            if ( decl%type == label ) then
              spec = get_spec(decl%tree)
              if ( spec /= s_matrix ) then
                call announce_error ( son, incompatible )
                what = w_nothing
   go to 100
              end if
            else if ( decl%type /= exprn_m ) then
              call announce_error ( son, incompatible )
              what = w_nothing
   go to 100
            end if
            call getActualMatrixFromDatabase ( matrixDatabase(value), &
              & LHSMatrix )
            select case ( what )
            case ( w_matrix )
              call copyMatrixValue ( LHSmatrix, matrix )     ! deep copy
            case ( w_matrix_c )
              call copyMatrixValue ( LHSmatrix, matrix_c%m ) ! deep copy
            case ( w_matrix_k )
              call copyMatrixValue ( LHSmatrix, matrix_k%m ) ! deep copy
            case ( w_matrix_s )
              call copyMatrixValue ( LHSmatrix, matrix_s%m ) ! deep copy
            end select
            if ( index(switches,'spa') /= 0 ) call dump_struct ( matrix, 'Result of expression' )
          case default
            stop
          end select
          call redeclare ( string, dvalue, decl%type, decl%units, decl%tree )
          value = decoration ( decl%tree ) ! in case switches contains 'alg'
100       call destroyStuff ( what, vector, matrix, matrix_c, matrix_k, matrix_s )
        end if
        if ( index(switches,'alg') /= 0 ) then
          select case ( what )
          case ( w_number )
            call display_string ( sub_rosa(lhs) )
            call output ( ' = ' )
            call output ( dValue, advance='yes' )
          case ( w_vector )
            call dump ( vectorDatabase(value) )
          case ( w_matrix, w_matrix_c, w_matrix_k, w_matrix_s )
            call dump ( matrixDatabase(value), details=2 )
          end select
        end if
      end if
200   if ( toggle(gen) ) call trace_end ( 'Algebra loop' )
    end do ! i_sons
    if ( toggle(gen) ) call trace_end ( 'Algebra' )

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'An error occurred in the Algebra section' )

  contains

    ! ...........................................  AlgebraCommands .....
    subroutine AlgebraCommands ( root, VectorDatabase, MatrixDatabase )
      use Tree, only: DECORATION, NSONS, SUBTREE
      use MoreTree, only: GET_SPEC_ID
      use MatrixModule_1, only: CYCLICJACOBI, REFLECTMATRIX, MATRIX_T
      use VectorsModule, only: VECTOR_T
      use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST, F_MATRIX, F_EIGENVECTORS
      use Init_Tables_Module, only: L_TRUE
      use Init_Tables_Module, only: S_CYCLICJACOBI, S_REFLECT
        
      ! Dummy arguments
      integer, intent(in) :: ROOT
      type(vector_T), dimension(:), pointer :: VectorDatabase
      type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase

      ! Local variables
      logical, dimension(field_first:field_last) :: GOT ! Fields
      integer :: FIELDID                ! ID for a field (duh!)
      integer :: I                      ! Loop counter
      integer :: KEY                    ! Tree node
      integer :: SON                    ! Tree node
      integer :: VALUE                  ! Value node
      integer :: MATRIXIND              ! Index into matrix database
      type(Matrix_T), pointer :: MATRIX ! The matrix to work on
      type(Matrix_T), pointer :: EIGENVECTORS ! Eigen vector matrix

      ! Executable code
      do i = 2, nsons(root)
        son = subtree ( i, root )
        key = subtree ( 1, son )
        if ( node_id(son) == n_set_one ) then
          value = l_true
        else
          value = decoration(subtree(2,son))
        end if
        fieldID = decoration(key)
        got ( fieldID ) = .true.
        select case ( fieldID )
        case ( f_matrix )
          matrixInd = decoration ( value )
          call GetFromMatrixDatabase ( matrixDatabase(matrixInd), matrix )
        case ( f_eigenVectors )
          matrixInd = decoration ( value )
          call GetFromMatrixDatabase ( matrixDatabase(matrixInd), eigenVectors )
        end select
      end do

      select case ( get_spec_id ( root ) )
      case ( s_reflect )
        call Dump_Struct ( matrix, 'Input matrix' )
        call ReflectMatrix ( matrix )
        call Dump_Struct ( matrix, 'Reflected matrix' )
      case ( s_cyclicJacobi ) 
        call Dump_Struct ( matrix, 'Input matrix' )
        call Dump_Struct ( eigenVectors, 'Input eigen vectors' )
        call CyclicJacobi ( matrix, eigenVectors, tol=1e-6 )
        call Dump_Struct ( matrix, 'Output matrix (eigen values)' )
        call Dump_Struct ( eigenVectors, 'Output eigen vectors' )
      end select

    end subroutine AlgebraCommands
      
    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What )
      use LEXER_CORE, only: PRINT_SOURCE
      use OUTPUT_M, only: OUTPUT
      use TREE, only: DUMP_TREE_NODE, SOURCE_REF
      integer :: Where        ! Root of subtree where error occurred
      integer :: What         ! What message to emit

      error = max(error,1)
      call output ( '***** At ', from_where = "Algebra module" )
      call print_source ( source_ref(where) )
      call output ( ', tree ' )
      call output ( where )
      call output ( ': ' )
      select case ( what )
      case ( ambiguous )
        call output ( 'Left-hand side is ambiguous.', advance='yes' )
      case ( cantInvertVector )
        call output ( 'Cannot invert a vector.', advance='yes' )
      case ( incompatible )
        call output ( 'Operands are incompatible.', advance='yes' )
      case ( undefined )
        call output ( 'Name in expression is undefined.', advance='yes' )
      case ( unknownFunc )
        call output ( 'Function "' )
        call display_string(sub_rosa(where))
        call output ( ' is not recognized.', advance='yes' )
      case ( unknownOp )
        call output ( 'Operator "' )
        call dump_tree_node ( where, 0 )
        call output ( '" in expression is unknown.', advance='yes' )
      case ( wrongNumArgs )
        call output ( 'Wrong number of arguments.', advance='yes' )
      case default
        call output ( 'What is error code ' )
        call output ( what )
        call output ( '?', advance='yes' )
      end select
    end subroutine Announce_Error

    ! ......................................................  Expr .....
    recursive subroutine Expr ( Root, What, DValue, Vector, Matrix, Matrix_C, &
      & Matrix_K, Matrix_S )
      ! Evaluate an expression by traversing its parse tree.
      ! All it handles are identifiers, numbers, add, subtract, multiply,
      ! divide and inverse divide (\).
      ! It doesn't handle strings, ranges, array constructors or relational
      ! operators.
      ! The output vectors or matrices are COPIES.  They must be destroyed
      ! or put into a database to avoid memory leaks.

      use DECLARATION_TABLE, only: DECLS, EMPTY, FUNCTION, GET_DECL
      use FUNCTIONS, only: F_CHOLESKY, F_TRANSPOSE
      use STRING_TABLE, only: FLOAT_VALUE
      use TREE, only: NODE_ID, NSONS, SUB_ROSA, SUBTREE
      use TREE_TYPES ! Everything, especially everything beginning with N_

      integer, intent(in) :: Root  ! Index in tree of root of expression subtree
      integer, intent(out) :: What  ! What is result (see its parameters above)
      ! Results, depending on WHAT:
      real(r8), intent(out) :: DValue                   ! What == w_number
      type(vector_t), intent(out) :: Vector             ! What == w_vector
      type(matrix_t), intent(out) :: Matrix             ! What == w_matrix
      type(matrix_cholesky_t), intent(out) :: Matrix_C  ! What == w_matrix_c
      type(matrix_kronecker_t), intent(out) :: Matrix_K ! What == w_matrix_k
      type(matrix_SPD_t), intent(out) :: Matrix_S       ! What == w_matrix_s
                                                 ! else   What == w_nothing

      type(decls) :: DECL          ! Declaration of a name
      real(r8) :: DValue2          ! Value of RH operand if Value2 == 0
      integer :: ME                ! Node ID for root
      integer :: Son1, Son2        ! Sons of Root
      integer :: SPEC              ! Index of spec of label, e.g. S_Matrix
      integer :: String            ! String table index
      integer :: What2             ! WHAT for second subexpression
      type(vector_t), save :: EmptyVector
      type(vector_t) :: Vector2
      type(matrix_t), save :: EmptyMatrix
      type(matrix_t) :: Matrix2, Matrix3
      type(matrix_t), pointer :: MatrixP
      type(matrix_cholesky_t) :: Matrix_C2
      type(matrix_cholesky_t), pointer :: Matrix_CP
      type(matrix_kronecker_t) :: Matrix_K2
      type(matrix_kronecker_t), pointer :: Matrix_KP
      type(matrix_SPD_t) :: Matrix_S2
      type(matrix_SPD_t), pointer :: Matrix_SP

      integer, save :: Whats(k_empty:k_spd)
      data whats(k_empty) / w_nothing /, whats(k_cholesky) / w_matrix_c /
      data whats(k_kronecker) / w_matrix_k /, whats(k_plain) / w_matrix /
      data whats(k_spd) / w_matrix_s /

!      call print_subtree ( root, 0, dump_decor=.true. )

      if ( toggle(gen) ) call trace_begin ( 'Algebra.Expr', root )
      dvalue = 0.0_r8
      me = node_id(root)
      select case ( me )
      case ( n_identifier ) ! ------------------------------------ Identifier
        string = sub_rosa(root)
        ! Look up the identifier
        decl = get_my_decl(root)
        if ( decl%type == label ) then
          spec = get_spec(decl%tree)
          ! call print_subtree ( decl%tree, 0, dump_decor=.true. )
          if ( spec == s_matrix ) then
            decl%type = exprn_m
          else if ( spec == s_vector ) then
            decl%type = exprn_v
          else
            decl%type = empty
          end if
        end if
        if ( decl%type == empty ) then
          call announce_error ( root, undefined )
          what = w_nothing
        else
          value = decoration ( decl%tree )
          dvalue = decl%value
          select case ( decl%type )
          case ( exprn, num_value )
            what = w_number
          case ( exprn_v )
            what = w_vector
            call copyVector ( vector, vectorDatabase(value), clone=.true. )
          case ( exprn_m )
            what = whats(getKindFromMatrixDatabase(matrixDatabase(value)))
            select case ( what )
            case ( w_matrix )
              call getFromMatrixDatabase ( matrixDatabase(value), matrixP )
              call copyMatrix ( matrix, matrixP )
            case ( w_matrix_c )
              call getFromMatrixDatabase ( matrixDatabase(value), matrix_CP )
              call copyMatrix ( matrix_C%m, matrix_CP%m )
            case ( w_matrix_k )
              call getFromMatrixDatabase ( matrixDatabase(value), matrix_KP )
              call copyMatrix ( matrix_K%m, matrix_KP%m )
            case ( w_matrix_s )
              call getFromMatrixDatabase ( matrixDatabase(value), matrix_SP )
              call copyMatrix ( matrix_S%m, matrix_SP%m )
            end select
          end select
        end if
      case ( n_number ) ! -------------------------------------------- Number
        dvalue = float_value(sub_rosa(root))
        what = w_number
      case ( n_func_ref )
        son1 = subtree(1,root)
        string = sub_rosa(son1)
        ! Look up the function name
        decl = get_decl(string,function)
        if ( decl%type /= function ) call announce_error ( son1, unknownFunc )
        if ( nsons(root) > 1 ) then
          son2 = subtree(2,root)
          call expr ( son2, what2, dvalue2, vector2, matrix2, matrix_c2, &
            & matrix_k2, matrix_s2 )
        end if
        select case ( decl%units )
        case ( f_cholesky )
          if ( nsons(root) /= 2 ) then
            call announce_error ( son1, wrongNumArgs )
          else
            if ( what2 /= w_matrix_s ) then
              call announce_error ( son2, incompatible )
            else
            end if
          end if
        case ( f_transpose )
          if ( nsons(root) /= 2 ) then
            call announce_error ( son1, wrongNumArgs )
          else
            if ( what2 /= w_matrix ) then
              call announce_error ( son2, incompatible )
            else
              ! Create result matrix as transpose
              call CreateEmptyMatrix ( matrix, 0, &
                & matrix2%col%vec, matrix2%row%vec, &
                & .not. matrix2%col%instFirst, .not. matrix2%row%instFirst )
              call TransposeMatrix ( matrix, matrix2 )
              what = w_matrix
              ! Optimization?  Instead of actually doing the transpose here,
              ! consider returning a flag that the result needs transposing.
              ! Then one could use multiply routines that are aware of the
              ! transpose, at the expense of needing to remember that the
              ! transpose needs to be done in other situations, e.g in "Add".
            end if
          end if
        end select
      case default
        son1 = subtree(1,root)
        call expr ( son1, what, dvalue, vector, matrix, matrix_c, &
          & matrix_k, matrix_s )
        what2 = w_nothing
        if ( nsons(root) > 1 ) then
          son2 = subtree(2,root)
          call expr ( son2, what2, dvalue2, vector2, matrix2, matrix_c2, &
            & matrix_k2, matrix_s2 )
        end if
        select case ( me )
        case ( n_plus, n_minus ) ! ------------------------------ Plus, Minus
          if ( nsons(root) > 1 ) then
            if ( what /= what2 ) then
              call announce_error ( root, incompatible )
              what = w_nothing
              return
            end if
            if ( me == n_plus ) then ! --------------------------------- Plus
              ! value = value + value2
              select case ( what )
              case ( w_number )
                dvalue = dvalue + dvalue2
              case ( w_vector )
                call addToVector ( vector, vector2 )
              case ( w_matrix )
                call addToMatrix ( matrix, matrix2 )
              case ( w_matrix_c )
                call addToMatrix ( matrix_c%m, matrix_c2%m )
              case ( w_matrix_k )
                call addToMatrix ( matrix_k%m, matrix_k2%m )
              case ( w_matrix_s )
                call addToMatrix ( matrix_s%m, matrix_s2%m )
              end select
            else ! ------------------------------------------------- Subtract
              ! value = value - value2
              select case ( what )
              case ( w_number )
                dvalue = dvalue - dvalue2
              case ( w_vector )
                call addToVector ( vector, vector2, scale=-1.0_r8 )
              case ( w_matrix )
                call addToMatrix ( matrix, matrix2, scale=-1.0_r8 )
              case ( w_matrix_c )
                call addToMatrix ( matrix_c%m, matrix_c2%m, scale=-1.0_r8 )
              case ( w_matrix_k )
                call addToMatrix ( matrix_k%m, matrix_k2%m, scale=-1.0_r8 )
              case ( w_matrix_s )
                call addToMatrix ( matrix_s%m, matrix_s2%m, scale=-1.0_r8 )
              end select
            end if
          else if ( me == n_minus ) then ! --------------------------- Negate
            ! value = -value
            select case ( what )
            case ( w_number )
              dvalue = -dvalue
            case ( w_vector )
              call scaleVector ( vector, -1.0_r8 )
            case ( w_matrix )
              call scaleMatrix ( matrix, -1.0_r8 )
            case ( w_matrix_c )
              call scaleMatrix ( matrix_c%m, -1.0_r8 )
            case ( w_matrix_k )
              call scaleMatrix ( matrix_k%m, -1.0_r8 )
            case ( w_matrix_s )
              call scaleMatrix ( matrix_s%m, -1.0_r8 )
            end select
          end if
        case ( n_mult ) ! value = value * value2 ------------------- Multiply
          select case ( what )
          case ( w_number )
            select case ( what2 )
            case ( w_number ) ! ........................ Number * Number
              dvalue = dvalue * dvalue2
            case ( w_vector ) ! ........................ Number * Vector
              call scaleVector ( vector2, dValue, vector )
            case ( w_matrix ) ! ........................ Number * Matrix
              call copyMatrix ( matrix, matrix2 )
              call scaleMatrix ( matrix, dValue )
            case ( w_matrix_c ) ! .................... Number * Matrix_C
              call copyMatrix ( matrix_c%m, matrix_c2%m )
              call scaleMatrix ( matrix_c%m, dValue )
            case ( w_matrix_k ) ! .................... Number * Matrix_K
              call copyMatrix ( matrix_k%m, matrix_k2%m )
              call scaleMatrix ( matrix_k%m, dValue )
            case ( w_matrix_s ) ! .................... Number * Matrix_S
              call copyMatrix ( matrix_s%m, matrix_s2%m )
              call scaleMatrix ( matrix_s%m, dValue )
            end select
            what = what2
          case ( w_vector )
            select case ( what2 )
            case ( w_number ) ! ........................ Vector * Number
              call scaleVector ( vector, dValue2 )
            case ( w_vector ) ! ........................ Vector * Vector
            case ( w_matrix ) ! ........................ Vector * Matrix
            case ( w_matrix_c ) ! .................... Vector * Matrix_C
            case ( w_matrix_k ) ! .................... Vector * Matrix_K
            case ( w_matrix_s ) ! .................... Vector * Matrix_S
            end select
          case ( w_matrix )
            select case ( what2 )
            case ( w_number ) ! ........................ Matrix * Number
              call scaleMatrix ( matrix, dValue2 )
            case ( w_vector ) ! ........................ Matrix * Vector
              call CloneVector ( vector, matrix%row%vec )
              call multiplyMatrixVectorNoT ( matrix, vector2, vector )
              what = w_vector
            case ( w_matrix ) ! ........................ Matrix * Matrix
              call multiplyMatrix_XY ( matrix, matrix2, matrix3 )
              call assignMatrix ( matrix, matrix3 ) ! Destroys Matrix first
            case ( w_matrix_c ) ! .................... Matrix * Matrix_C
            case ( w_matrix_k ) ! .................... Matrix * Matrix_K
            case ( w_matrix_s ) ! .................... Matrix * Matrix_S
            end select
          case ( w_matrix_c )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_C * Number
              call scaleMatrix ( matrix_c%m, dValue2 )
            case ( w_vector ) ! ...................... Matrix_C * Vector
            case ( w_matrix ) ! ...................... Matrix_C * Matrix
            case ( w_matrix_c ) ! .................. Matrix_C * Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_C * Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_C * Matrix_S
            end select
          case ( w_matrix_k )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_K * Number
              call scaleMatrix ( matrix_k%m, dValue2 )
            case ( w_vector ) ! ...................... Matrix_K * Vector
            case ( w_matrix ) ! ...................... Matrix_K * Matrix
            case ( w_matrix_c ) ! .................. Matrix_K * Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_K * Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_K * Matrix_S
            end select
          case ( w_matrix_s )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_S * Number
              call scaleMatrix ( matrix_s%m, dValue2 )
            case ( w_vector ) ! ...................... Matrix_S * Vector
              call multiplyMatrixVectorSPD_1 ( matrix_s, vector2, vector )
            case ( w_matrix ) ! ...................... Matrix_S * Matrix
            case ( w_matrix_c ) ! .................. Matrix_S * Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_S * Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_S * Matrix_S
            end select
          end select
        case ( n_div )  ! value = value / value2 ------------------ Divide By
          select case ( what )
          case ( w_number )
            select case ( what2 )
            case ( w_number ) ! ........................ Number / Number
              dvalue = dvalue / dvalue2
            case ( w_vector ) ! ........................ Number / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ........................ Number / Matrix
            case ( w_matrix_c ) ! .................... Number / Matrix_C
            case ( w_matrix_k ) ! .................... Number / Matrix_K
            case ( w_matrix_s ) ! .................... Number / Matrix_S
            end select
          case ( w_vector )
            select case ( what2 )
            case ( w_number ) ! ........................ Vector / Number
              call scaleVector ( vector, 1.0_r8 / dValue2 )
            case ( w_vector ) ! ........................ Vector / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ........................ Vector / Matrix
            case ( w_matrix_c ) ! .................... Vector / Matrix_C
            case ( w_matrix_k ) ! .................... Vector / Matrix_K
            case ( w_matrix_s ) ! .................... Vector / Matrix_S
            end select
          case ( w_matrix )
            select case ( what2 )
            case ( w_number ) ! ........................ Matrix / Number
              call scaleMatrix ( matrix, 1.0_r8 / dValue2 )
            case ( w_vector ) ! ........................ Matrix / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ........................ Matrix / Matrix
            case ( w_matrix_c ) ! .................... Matrix / Matrix_C
            case ( w_matrix_k ) ! .................... Matrix / Matrix_K
            case ( w_matrix_s ) ! .................... Matrix / Matrix_S
            end select
          case ( w_matrix_c )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_C / Number
              call scaleMatrix ( matrix_c%m, 1.0_r8 / dValue2 )
            case ( w_vector ) ! ...................... Matrix_C / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ...................... Matrix_C / Matrix
            case ( w_matrix_c ) ! .................. Matrix_C / Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_C / Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_C / Matrix_S
            end select
          case ( w_matrix_k )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_K / Number
              call scaleMatrix ( matrix_k%m, 1.0_r8 / dValue2 )
            case ( w_vector ) ! ...................... Matrix_K / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ...................... Matrix_K / Matrix
            case ( w_matrix_c ) ! .................. Matrix_K / Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_K / Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_K / Matrix_S
            end select
          case ( w_matrix_s )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_S / Number
              call scaleMatrix ( matrix_s%m, 1.0_r8 / dValue2 )
            case ( w_vector ) ! ...................... Matrix_S / Vector
              call announce_error ( son2, cantInvertVector )
            case ( w_matrix ) ! ...................... Matrix_S / Matrix
            case ( w_matrix_c ) ! .................. Matrix_S / Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_S / Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_S / Matrix_S
            end select
          end select
        case ( n_into ) ! value = value \ value2 ---------------- Divide Into
          select case ( what )
          case ( w_number )
            select case ( what2 )
            case ( w_number ) ! ........................ Number \ Number
              dvalue = dvalue2 / dvalue
            case ( w_vector ) ! ........................ Number \ Vector
              call scaleVector ( vector2, 1.0_r8/dValue, vector )
            case ( w_matrix ) ! ........................ Number \ Matrix
            case ( w_matrix_c ) ! .................... Number \ Matrix_C
            case ( w_matrix_k ) ! .................... Number \ Matrix_K
            case ( w_matrix_s ) ! .................... Number \ Matrix_S
            end select
          case ( w_vector )
            call announce_error ( son1, cantInvertVector )
          case ( w_matrix )
            select case ( what2 )
            case ( w_number ) ! ........................ Matrix \ Number
            case ( w_vector ) ! ........................ Matrix \ Vector
            case ( w_matrix ) ! ........................ Matrix \ Matrix
            case ( w_matrix_c ) ! .................... Matrix \ Matrix_C
            case ( w_matrix_k ) ! .................... Matrix \ Matrix_K
            case ( w_matrix_s ) ! .................... Matrix \ Matrix_S
            end select
          case ( w_matrix_c )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_C \ Number
            case ( w_vector ) ! ...................... Matrix_C \ Vector
            case ( w_matrix ) ! ...................... Matrix_C \ Matrix
            case ( w_matrix_c ) ! .................. Matrix_C \ Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_C \ Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_C \ Matrix_S
            end select
          case ( w_matrix_k )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_K \ Number
            case ( w_vector ) ! ...................... Matrix_K \ Vector
            case ( w_matrix ) ! ...................... Matrix_K \ Matrix
            case ( w_matrix_c ) ! .................. Matrix_K \ Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_K \ Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_K \ Matrix_S
            end select
          case ( w_matrix_s )
            select case ( what2 )
            case ( w_number ) ! ...................... Matrix_S \ Number
            case ( w_vector ) ! ...................... Matrix_S \ Vector
            case ( w_matrix ) ! ...................... Matrix_S \ Matrix
            case ( w_matrix_c ) ! .................. Matrix_S \ Matrix_C
            case ( w_matrix_k ) ! .................. Matrix_S \ Matrix_K
            case ( w_matrix_s ) ! .................. Matrix_S \ Matrix_S
            end select
          end select
        case default
          call announce_error ( root, unknownOp )
        end select
        call destroyStuff ( what2, vector2, matrix2, matrix_c2, matrix_k2, &
          & matrix_s2 )
      end select
      if ( toggle(gen) ) call trace_end ( 'Algebra.Expr' )
    end subroutine Expr

    ! .............................................  DestroyStuff  .....
    subroutine DestroyStuff ( what, vector, matrix, matrix_c, matrix_k, matrix_s )
    ! Destroy whatever argument WHAT says to destroy
      integer, intent(in) :: WHAT
      type(vector_t), intent(inout) :: Vector
      type(matrix_t), intent(inout) :: Matrix
      type(matrix_cholesky_t), intent(out) :: Matrix_C
      type(matrix_kronecker_t), intent(out) :: Matrix_K
      type(matrix_SPD_t), intent(out) :: Matrix_S
      select case ( what )
      case ( w_vector )
        call destroyVectorInfo ( vector )
      case ( w_matrix )
        call destroyMatrix ( matrix )
      case ( w_matrix_c )
        call destroyMatrix ( matrix_c%m )
      case ( w_matrix_k )
        call destroyMatrix ( matrix_k%m )
      case ( w_matrix_s )
        call destroyMatrix ( matrix_s%m )
      end select
    end subroutine DestroyStuff

    ! ..............................................  Get_My_Decl  .....
    type(decls) function Get_My_Decl ( Root )
      integer, intent(in) :: Root      ! Tree index of declaration to get
      type(decls) :: DECL              ! Declaration
      integer :: I
      integer :: String                ! String index for Root
      integer, parameter :: Try(5) = (/ exprn, exprn_m, exprn_v, label, &
        &                               num_value /) ! Decl types to try
      get_my_decl%type = empty
      string = sub_rosa(root)
      do i = 1, size(try)
        decl = get_decl(string,try(i))
        if ( decl%type /= empty ) then
          if ( get_my_decl%type /= empty ) call announce_error ( root, ambiguous )
          get_my_decl = decl
        end if
      end do
    end function Get_My_Decl

    ! .................................................  Get_Spec  .....
    integer function Get_Spec ( DEF )
    ! "Root" is the index of a spec_args node for which we have the label.
    ! Get the index of the specification, e.g. S_Matrix or S_Vector.
      use Declaration_Table, only: SPEC
      integer, intent(in) :: DEF        ! Tree node of definition of declaration
      type(decls) :: DECL               ! Declaration 
      get_spec = 0
      if ( nsons ( def ) == 0 ) return
      get_spec = decoration ( subtree ( 1,def ) )
      if ( nsons ( get_spec ) == 0 ) return
      get_spec = decoration ( subtree ( 1, get_spec ) )
    end function Get_Spec

  end subroutine Algebra

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ALGEBRA_M

! $Log$
! Revision 2.8  2004/01/29 03:33:26  livesey
! Hooked up reflect and cyclicJacobi commands.
!
! Revision 2.7  2004/01/24 01:04:52  livesey
! More bug fixes and population
!
! Revision 2.6  2004/01/23 05:38:22  livesey
! Various bug fixes
!
! Revision 2.5  2004/01/21 00:30:25  vsnyder
! Move some stuff to more logical places, fix get_my_decl, cosmetics
!
! Revision 2.4  2004/01/20 23:19:41  vsnyder
! Better error handling, plug some memory leaks
!
! Revision 2.3  2004/01/20 19:33:35  vsnyder
! Correct bug in identifier handling; delete unused variable declaration
!
! Revision 2.2  2004/01/17 03:04:15  vsnyder
! Provide for functions in expressions
!
! Revision 2.1  2004/01/17 00:28:32  vsnyder
! Initial commit
!
