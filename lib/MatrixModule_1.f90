! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MatrixModule_1          ! Block Matrices in the MLS PGS suite
!=============================================================================

! This module provides a block matrix type including operations for matrix
! quantities in MLS Level 2 software, and related programs.

  use MatrixModule_0, only: CholeskyFactor, ColumnScale, CreateBlock, &
    & DestroyBlock, M_Absent, MatrixElement_T, operator(+), operator(.TX.), &
    & RowScale, SolveCholesky, UpdateDiagonal
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use VectorsModule, only: CloneVector, Vector_T

  implicit NONE
  private
  public :: AddMatrices, CreateEmptyMatrix, CholeskyFactor, CholeskyFactor_1
  public :: ColumnScale, ColumnScale_1, Dump
  public :: FindBlock
! public :: LevenbergUpdateCholesky
  public :: Matrix_T, Matrix_Cholesky_T
  public :: Matrix_Kronecker_T, Matrix_SPD_T, MultiplyMatrices
  public :: MultiplyMatrixVector, MultiplyMatrixVector_1
  public :: NewMultiplyMatrixVector, NormalEquations
  public :: operator(.TX.), operator(+), RC_Info, RowScale, RowScale_1
  public :: SolveCholesky, SolveCholesky_1
  public :: UpdateDiagonal, UpdateDiagonal_1

! =====     Defined Operators and Generic Identifiers     ==============

  interface CholeskyFactor
    module procedure CholeskyFactor_1
  end interface

  interface ColumnScale
    module procedure ColumnScale_1
  end interface

  interface DUMP
    module procedure DUMP_MATRIX
  end interface

  interface MultiplyMatrixVector   ! A^T V
    module procedure MultiplyMatrixVector_1
  end interface

  interface operator (+)
    module procedure AddMatrices
  end interface

  interface operator ( .TX. )      ! A^T B
    module procedure MultiplyMatrices, NewMultiplyMatrixVector
  end interface

  interface RowScale
    module procedure RowScale_1
  end interface

  interface SolveCholesky
    module procedure SolveCholesky_1
  end interface

  interface UpdateDiagonal
    module procedure UpdateDiagonal_1
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  type RC_Info
  ! Information about the row or column of a matrix
    type(Vector_T), pointer :: Vec ! Vector used to define the row or column
      ! space of the matrix, if any.
    integer :: NB                  ! Number of blocks of rows or columns
    logical :: InstFirst = .true.  ! TRUE means horizontal instance is the
      ! major order and quantity is the minor order.
    integer, pointer :: Nelts(:) => NULL()  ! Numbers of rows or columns in
      ! each row or column of blocks.
    integer, pointer :: Inst(:) => NULL()   ! The instance indices for the
      ! rows or columns of blocks.
    integer, pointer :: Quant(:) => NULL()  ! The quantity indices for the
      ! rows or columns of blocks.  These are indices within Vec, not the
      ! quantity template database.
  end type RC_Info

  type Matrix_T
    integer :: Name      ! Sub-rosa index of matrix name, if any, else zero
    type(RC_Info) :: Col, Row  ! Column and row info
    type(matrixElement_T), dimension(:,:), pointer :: BLOCK => NULL()
  end type Matrix_T

  type Matrix_Cholesky_T ! Cholesky factored matrix.  Only the upper triangle
    type(matrix_T) :: M  ! is stored.
  end type Matrix_Cholesky_T

  type Matrix_Kronecker_T  ! Kronecker product matrix.  Not yet implemented.
    type(matrix_T) :: M
  end type Matrix_Kronecker_T

  type Matrix_SPD_T      ! Symmetric positive-definite matrix.  Only the
    type(matrix_T) :: M  ! upper triangle is stored.
  end type Matrix_SPD_T

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  AddMatrices  -----
  function AddMatrices ( X, Y ) result ( Z ) ! Z = X + Y
    type(Matrix_T), intent(in) :: X, Y
    type(Matrix_T) :: Z

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyMatrix using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J      ! Subscripts for [XYZ]%Block

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%col%vec%template,y%col%vec%template) &
      & .or. .not. associated(x%row%vec%template,y%row%vec%template) &
      & .or. (x%col%instFirst .neqv. y%col%instFirst) &
      & .or. (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMSG ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in AddMatrices" )
    call createEmptyMatrix ( z, 0, x%row%vec, x%col%vec )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        z%block(i,j) = x%block(i,j) + y%block(i,j)
      end do ! i = 1, nr
    end do ! j = 1, nc
  end function AddMatrices

  ! -------------------------------------------  CholeskyFactor_1  -----
  subroutine CholeskyFactor_1 ( X, Z )
  ! Compute the Cholesky factor Z of the matrix X.  Z%M%Block can be
  ! associated with X%M%Block to save space.
    type(Matrix_SPD_T), intent(in) :: X ! Matrix to factor.
    type(Matrix_Cholesky_T), intent(inout) :: Z   ! Factored matrix.

    integer :: I, J, K                  ! Subscripts and loop inductors
    integer :: N, N0                    ! Columns(blocks), Columns(0-level)
    type(MatrixElement_T) :: P          ! Product of two blocks
    type(MatrixElement_T) :: S          ! Sum, to accumulate "inner product"

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%m%col%vec%template,z%m%col%vec%template) &
      & .or. .not. associated(x%m%row%vec%template,z%m%row%vec%template) &
      & .or. (x%m%col%instFirst .neqv. z%m%col%instFirst) &
      & .or. (x%m%row%instFirst .neqv. z%m%row%instFirst) ) &
        & call MLSMSG ( MLSMSG_Error, ModuleName, &
          & "Matrices in CholeskyFactor are not compatible" )
    n = x%m%row%nb
    do i = 1, n
      n0 = x%m%col%nelts(i)
      call createBlock ( s, n0, n0, M_Absent )
      do k = 1, i-1
        p = z%m%block(k,i) .tx. z%m%block(k,i)
        s = s + p
        call destroyBlock( p )     ! Avoid a memory leak
      end do ! k = 1, i-1
      call choleskyFactor ( z%m%block(i,i), s ) ! Factor block on diagonal
      do j = i+1, n
        call destroyBlock ( s )    ! Avoid a memory leak
        do k = 1, i-1
          p = z%m%block(k,i) .tx. z%m%block(k,j)
          s = s + p
          call destroyBlock( p )   ! Avoid a memory leak
        end do ! k = 1, i-1
        call solveCholesky ( z%m%block(i,i), z%m%block(i,j), s, transpose=.true. )
      end do ! j = 1, n
      call destroyBlock( s )       ! Avoid a memory leak
    end do ! i = 1, n
  end subroutine CholeskyFactor_1

  ! ----------------------------------------------  ColumnScale_1  -----
  function ColumnScale_1 ( X, V ) result ( Z ) ! Z = X V where V is a
  !                                diagonal matrix represented by a vector.
    type (Matrix_T), intent(in) :: X
    type (Vector_T), intent(in) :: V
    type (Matrix_T) :: Z

    integer :: I, J      ! Subscripts for [XZ]%Block

    call createEmptyMatrix ( z, 0, x%row%vec, x%col%vec, &
      & x%row%instFirst, x%col%instFirst )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
      end do ! i = x%row%nb
        z%block(i,j) = columnScale ( x%block(i,j), &
          & v%quantities(x%col%quant(j))%values(:,x%col%inst(j)) )
    end do ! j = x%col%nb
  end function ColumnScale_1

  ! ----------------------------------------------  CreateBlock_1  -----
  subroutine CreateBlock_1 ( Z, RowNum, ColNum, Kind, NumNonzeros )
  ! Create the matrix block Z%Block(RowNum,ColNum), which sprang into
  ! existence with kind M_Absent.  Create it with the specified Kind.
  ! See MatrixModule_0 for a list of the kinds.  If the Kind is
  ! M_Banded or M_ColumnSparse, the number of nonzeroes is needed.
    type(matrix_T), intent(inout) :: Z       ! The matrix having the block
    integer, intent(in) :: RowNum, ColNum    ! Row and column of the block
    integer, intent(in) :: Kind         ! Kind of block, see MatrixModule_0
    integer, intent(in), optional :: NumNonzeros  ! Number of nonzeros
    call createBlock ( z%block(rowNum,colNum), &
      & z%row%nelts(rowNum), z%col%nelts(colNum), kind, numNonzeros )
  end subroutine CreateBlock_1

  ! ------------------------------------------  CreateEmptyMatrix  -----
  subroutine CreateEmptyMatrix ( Z, Name, Row, Col &
    &,                           Row_Quan_First, Col_Quan_First )
    type(Matrix_T), intent(out) :: Z    ! The matrix to create
    integer, intent(in) :: Name         ! Sub-rosa index of its name, or zero
    type(Vector_T), pointer :: Row      ! Vector used to define the row
      !                                   space of the matrix.
    type(Vector_T), pointer :: Col      ! Vector used to define the column
      !                                   space of the matrix.
    logical, intent(in), optional :: Row_Quan_First    ! True (default false)
      ! means the quantity is the major order of the rows of blocks and the
      ! instance is the minor order.
    logical, intent(in), optional :: Col_Quan_First    ! True (default false)
      ! means the quantity is the major order of the columns of blocks and the
      ! instance is the minor order.

    integer :: STATUS    ! From ALLOCATE

    z%name = name
    call defineInfo ( z%row, row, row_Quan_First )
    call defineInfo ( z%col, col, col_Quan_First )
    allocate ( z%block(z%row%nb,z%col%nb), stat=status )
    if ( status /= 0 ) call MLSMSG ( MLSMSG_Allocate, ModuleName, &
      & "Z%Block in CreateEmptyMatrix" )

  contains
    subroutine DefineInfo ( RC, Vec, QuanFirst )
      type(RC_Info), intent(out) :: RC
      type(Vector_T), pointer :: Vec ! intent(in)
      logical, intent(in), optional :: QuanFirst

      integer :: I, J, N      ! Subscripts or loop inductors
      logical :: NEW          ! Was an instance seen?

      rc%vec => vec
      rc%instFirst = .true.
      if ( present(quanFirst) ) rc%instFirst = .not. quanFirst
      rc%nb = row%template%totalInstances
      call allocate_test ( rc%nelts, rc%nb, &
        & "rc%nelts in CreateEmptyMatrix", ModuleName )
      call allocate_test ( rc%inst, rc%nb, "rc%inst in CreateEmptyMatrix", &
        & ModuleName )
      call allocate_test ( rc%quant, rc%nb, "rc%quant in CreateEmptyMatrix", &
        & ModuleName )
      if ( rc%instFirst ) then
!??? Are rc%nelts etc. different if the vector is not regular?
        n = 0
        j = 0
        do ! ( until .not. new )
          j = j + 1           ! instance number
          new = .false.       ! Instance j not seen for any quantity
          do i = 1, size(vec%quantities)
            if ( size(vec%quantities(i)%values,2) >= j ) then
              n = n + 1
              rc%nelts(n) = size(vec%quantities(i)%values,1)
              rc%inst(n) = j
              rc%quant(n) = i
              new = .true.    ! Instance j seen for some quantity
            end if
          end do ! i
          if ( .not. new ) exit
        end do ! j
      else
!??? Are rc%nelts etc. different if the vector is not regular?
        n = 0
        do i = 1, size(vec%quantities)
          do j = 1, size(vec%quantities(i)%values,2)
            n = n + 1
            rc%nelts(n) = size(vec%quantities(i)%values,1)
            rc%inst(n) = j
            rc%quant(n) = i
          end do ! j
        end do ! i
      end if
    end subroutine DefineInfo
  end subroutine CreateEmptyMatrix

  ! ----------------------------------------------  DestroyMatrix  -----
  subroutine DestroyMatrix ( A )
  ! Destroy a matrix -- deallocate its pointer components
    type(matrix_T), intent(inout) :: A

    integer :: I, J                ! Subscripts and look inductors
    integer :: STATUS              ! From deallocate

    do j = 1, a%col%nb
      do i = 1, a%row%nb
        call destroyBlock ( a%block(i,j) )
      end do ! i
    end do ! j
    deallocate ( a%block, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "A%Block in DestroyMatrix" )
    call destroyInfo ( a%row )
    call destroyInfo ( a%col )

  contains
    subroutine DestroyInfo ( RC )
      type(RC_Info), intent(out) :: RC
      call deallocate_test ( rc%nelts,  &
        & "rc%nelts in DestroyMatrix", ModuleName )
      call deallocate_test ( rc%inst, "rc%inst in DestroyMatrix", &
        & ModuleName )
      call deallocate_test ( rc%quant, "rc%quant in DestroyMatrix", &
        & ModuleName )
    end subroutine DestroyInfo
  end subroutine DestroyMatrix

  ! --------------------------------------------------  FindBlock  -----
  integer function FindBlock ( RC, Quantity, Instance )
  ! Given quantity and instance numbers, find a block index.
  ! This function can be used with either row or column information.
  ! Zero is returned if there is no block having the desired quantity
  ! and instance numbers.
    type(RC_Info), intent(in) :: RC     ! Row or Col component of a Matrix_T
    integer, intent(in) :: Quantity, Instance

    ! We could do something fancier here, e.g. a binary search on rc%quant
    ! or rc%inst, depending on rc%instFirst, but the difference between
    ! that and what's here is sure to be insignificant compared to the
    ! linear algebra.
    do findBlock = 1, rc%nb
      if ( rc%quant(findBlock) == quantity .and. &
        &  rc%inst(findBlock) == instance ) return
    end do
    findBlock = 0
  end function FindBlock

  ! ------------------------------------  LevenbergUpdateCholesky  -----
! subroutine LevenbergUpdateCholesky ( Z, LAMBDA )
! ! Given a Cholesky factor Z of a matrix of normal equations A^T A,
! ! update the Cholesky factor to incorporate Levenberg-Marquardt
! ! stabilization that corresponds to augmenting A with LAMBDA I.
!   type(Matrix_T), intent(inout) :: Z
!   real(r8), intent(in) :: LAMBDA
! end subroutine LevenbergUpdateCholesky

  ! -------------------------------------------  MultiplyMatrices  -----
  function MultiplyMatrices ( X, Y ) result ( Z ) ! Z = X^T Y
    type(Matrix_T), intent(in) :: X, Y
    type(Matrix_T) :: Z

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyMatrix using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K             ! Subscripts for [XYZ]%Block
    type(MatrixElement_T) :: T     ! So that we can clean up after a
    !                                low-level multiply

    ! Check that matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%row%vec%template,y%row%vec%template) .or. &
      & (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMSG ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in MultiplyMatrices" )
    call createEmptyMatrix ( z, 0, x%col%vec, y%col%vec )
    do j = 1, y%col%nb
      do i = 1, x%col%nb
        z%block(i,j) = x%block(1,i) .tx. y%block(1,j)
        do k = 2, x%row%nb
          t = x%block(k,i) .tx. y%block(k,j) ! ??? It may be desirable to
          z%block(i,j) = z%block(i,j) + t    ! ??? improve this by having a
                                             ! ??? low-level multiply-add.
          call destroyBlock ( t )       ! Avoid a memory leak
        end do ! k = 2, x%nr
      end do ! i = 1, x%nc
    end do ! j = 1, y%nc
  end function MultiplyMatrices

  ! -------------------------------------  MultiplyMatrixVector_1  -----
  subroutine MultiplyMatrixVector_1 ( A, V, Z, UPDATE )
  ! Z = A^T V if UPDATE is absent or false.
  ! Z = Z + A^T V is UPDATE is present and true.
    type(Matrix_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: Z
    logical, optional, intent(in) :: UPDATE

    logical :: DO_UPDATE      ! Tells multiply_Matrix_Vector whether to
    !                           clear an element of Z, or add to it
    integer :: I, J           ! Subscripts and loop inductors
    integer :: K, L, M, N     ! Subscripts
    logical :: MY_UPDATE      ! My copy of UPDATE or false if it's absent

    if ( .not. associated(a%row%vec%template, v%template) ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix and vector not compatible in MultiplyMatrixVector_1" )
    my_update = .false.
    if ( present(update) ) my_update = update
    call cloneVector ( z, v ) ! Copy characteristics, allocate values
    do j = 1, a%col%nb
      k = a%col%quant(j)
      l = a%col%inst(j)
      do_update = my_update
      do i = 1, a%row%nb
        m = a%row%quant(i)
        n = a%row%inst(i)
        call multiply_Matrix_Vector ( a%block(i,j), &
          & v%quantities(m)%values(:,n), z%quantities(k)%values(:,l), &
          & do_update )
        do_update = .true.
      end do ! i = 1, a%row%nb
    end do ! j = 1, a%col%nb
  end subroutine MultiplyMatrixVector_1

  ! ------------------------------------  NewMultiplyMatrixVector  -----
  function NewMultiplyMatrixVector ( A, V ) result ( Z ) ! Z = A^T V
    type(Matrix_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T) :: Z
    call multiplyMatrixVector ( a, v, z, .false. )
  end function NewMultiplyMatrixVector

  ! --------------------------------------------  NormalEquations  -----
  subroutine NormalEquations ( A, Z, RHS_IN, RHS_OUT, UPDATE )
  ! If UPDATE is absent, or present but false, form normal equations of the
  ! least-squares problem A X = RHS_IN. Z = A^T A and RHS_OUT = A^T RHS_IN.
  ! If UPDATE is present and true, update normal equations of the least-
  ! squares problem A X = RHS_IN. Z = Z + A^T A and RHS_OUT = RHS_OUT +
  ! A^T RHS_IN.
  ! Only the upper triangle of A^T A is formed or updated.
    type(Matrix_T), intent(in) :: A
    type(Matrix_SPD_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: RHS_IN
    type(Vector_T), intent(inout) :: RHS_OUT
    logical, intent(in), optional :: UPDATE  ! True (default false) means
    !                                          to update Z and RHS_OUT

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! If this subroutine is invoked with UPDATE absent, or present and false,
  ! It is important to invoke DestroyMatrix using the Z argument of this
  ! subroutine after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K             ! Subscripts for [AZ]%Block
    logical :: MY_UPDATE
    type(MatrixElement_T) :: T     ! So that we can clean up after a
    !                                low-level multiply

    ! Compute Z = A^T A or Z = Z + A^T A
    my_update = .false.
    if ( present(update) ) my_update = update
    if ( .not. my_update ) &
      & call createEmptyMatrix ( z%m, 0, a%col%vec, a%col%vec )
    do j = 1, a%col%nb
      do i = 1, j
        z%m%block(i,j) = a%block(1,i) .tx. a%block(1,j)
        do k = 2, a%row%nb
          t = a%block(k,i) .tx. a%block(k,j) ! ??? It may be desirable to
          z%m%block(i,j) = z%m%block(i,j) + t    ! ??? improve this by having a
                                             ! ??? low-level multiply-add.
          call destroyBlock ( t )   ! Avoid a memory leak
        end do ! k = 2, a%row%nb
      end do ! i = 1, a%col%nb
    end do ! j = 1, a%col%nb

    call multiplyMatrixVector ( a, rhs_in, rhs_out, my_update )
  end subroutine NormalEquations

  ! -------------------------------------------------  RowScale_1  -----
  function RowScale_1 ( V, X ) result ( Z ) ! Z = V X where V is a
  !                                diagonal matrix represented by a vector.
    type (Vector_T), intent(in) :: V
    type (Matrix_T), intent(in) :: X
    type (Matrix_T) :: Z

    integer :: I, J      ! Subscripts for [XZ]%Block

    call createEmptyMatrix ( z, 0, x%row%vec, x%col%vec, &
      & x%row%instFirst, x%col%instFirst )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
      end do ! i = x%row%nb
        z%block(i,j) = &
          & RowScale ( v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
          & x%block(i,j) )
    end do ! j = x%col%nb
  end function RowScale_1

  ! --------------------------------------------  SolveCholesky_1  -----
  subroutine SolveCholesky_1 ( Z, X, RHS, TRANSPOSE )
  ! Given the Cholesky-factored normal equations Z and the corresponding
  ! RHS, Solve Z X = RHS for X if TRANSPOSE is absent, or present and
  ! false.  Solve Z^T X = RHS for X if TRANSPOSE is present and true.
  ! RHS may be the same as X.  RHS may be absent, in which case X is
  ! assumed to contain the right-hand side on input, and the solution
  ! replaces it on output.
    type(Matrix_Cholesky_T), intent(in) :: Z
    type(Vector_T), intent(inout), target :: X
    type(Vector_T), intent(in), target, optional :: RHS
    logical, optional, intent(in) :: TRANSPOSE

    integer :: I, J                ! Subscripts and loop inductors
    integer :: IC, IR, QC, QR      ! Instance and quantity indices
    type(Vector_T), pointer :: MY_RHS   ! RHS if present, else X
    logical My_transpose           ! TRANSPOSE if present, else .false.

    my_transpose = .false.
    if ( present(transpose) ) my_transpose = transpose
    my_rhs => x
    if ( present(rhs) ) then
      if ( .not. associated( x%template, rhs%template ) ) &
        & call MLSMSG ( MLSMSG_Error, ModuleName, &
          & "X and RHS not compatible in SolveCholesky_1" )
      my_rhs => rhs
    end if
    if ( .not. associated( z%m%col%vec%template, my_rhs%template ) ) &
      & call MLSMSG ( MLSMSG_Error, ModuleName, &
        & "Z and RHS not compatible in SolveCholesky_1" )

    if ( my_transpose ) then       ! Solve Z^T X = RHS for X
      do i = 1, z%m%col%nb
        ic = z%m%col%inst(i)
        qc = z%m%col%quant(i)
        do j = 1, i-1
          ir = z%m%row%inst(j)
          qr = z%m%row%quant(j)
          my_rhs%quantities(qc)%values(:,ic) = &
            & my_rhs%quantities(qc)%values(:,ic) - &
            & ( z%m%block(j,i) .tx. x%quantities(qr)%values(:,ir) )
        end do ! j = 1, i-1
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.true. )
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.false. )
      end do ! i = 1, z%m%col%nb
    else                           ! Solve Z X = RHS for X
      do i = z%m%row%nb, 1, -1
        ic = z%m%row%inst(i)
        qc = z%m%row%quant(i)
        do j = i+1, z%m%col%nb
          ir = z%m%col%inst(j)
          qr = z%m%col%quant(j)
          my_rhs%quantities(qc)%values(:,ic) = &
            & my_rhs%quantities(qc)%values(:,ic) - &
            & ( z%m%block(i,j) .tx. x%quantities(qr)%values(:,ir) )
        end do ! j = 1, i-1
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.true. )
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.false. )
      end do ! i = 1, z%m%col%nb
    end if
  end subroutine SolveCholesky_1

  ! -------------------------------------------  UpdateDiagonal_1  -----
  subroutine UpdateDiagonal_1 ( A, LAMBDA )
  ! Add LAMBDA to the diagonal of A.
    type(Matrix_SPD_T), intent(inout) :: A
    real(r8), intent(in) :: LAMBDA

    integer :: I

    do i = 1, a%m%row%nb
      call updateDiagonal ( a%m%block(i,i), lambda )
    end do
  end subroutine UpdateDiagonal_1

! =====     Private Procedures     =====================================

  ! -------------------------------------------------  DumpMatrix  -----
  subroutine Dump_Matrix ( Matrix, Name, Details )
    type(Matrix_T), intent(in) :: Matrix
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional :: Details
  end subroutine Dump_Matrix

end module MatrixModule_1

! $Log$
! Revision 2.1  2000/11/09 01:23:23  vsnyder
! Initial entry -- still under construction
!
