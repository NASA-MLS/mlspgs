! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MatrixModule_1          ! Block Matrices in the MLS PGS suite
!=============================================================================

! This module provides a block matrix type including operations for matrix
! quantities in MLS Level 2 software, and related programs.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use MatrixModule_0, only: Assignment(=), CholeskyFactor, ClearRows, &
    & ColumnScale, CreateBlock, DestroyBlock, Dump, M_Absent, MatrixElement_T, &
    & MultiplyMatrixBlocks, MultiplyMatrixVector, MultiplyMatrixVectorNoT, &
    & operator(+), operator(.TX.), RowScale, SolveCholesky, UpdateDiagonal
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use String_Table, only: Display_String
  use VectorsModule, only: CloneVector, Vector_T

  implicit NONE
  private
  public :: AddMatrices, AddToMatrix, Assignment(=), ClearRows, ClearRows_1
  public :: CopyMatrix, CreateEmptyMatrix, CholeskyFactor, CholeskyFactor_1
  public :: ColumnScale, ColumnScale_1, Dump, FindBlock
! public :: LevenbergUpdateCholesky
  public :: Matrix_T, Matrix_Cholesky_T, Matrix_Kronecker_T, Matrix_SPD_T
  public :: MultiplyMatrices, MultiplyMatrixVector, MultiplyMatrixVector_1
  public :: MultiplyMatrixVectorSPD_1
  public :: NewMultiplyMatrixVector, NormalEquations
  public :: operator(.TX.), operator(+), RC_Info, RowScale, RowScale_1
  public :: SolveCholesky, SolveCholesky_1
  public :: UpdateDiagonal, UpdateDiagonal_1

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment(=)
    module procedure AssignMatrix
  end interface

  interface CholeskyFactor
    module procedure CholeskyFactor_1
  end interface

  interface ClearRows
    module procedure ClearRows_1
  end interface

  interface ColumnScale
    module procedure ColumnScale_1
  end interface

  interface DUMP
    module procedure DUMP_MATRIX
  end interface

  interface MultiplyMatrixVector   ! A^T V
    module procedure MultiplyMatrixVector_1, MultiplyMatrixVectorSPD_1
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
    integer :: Name = 0  ! Sub-rosa index of matrix name, if any, else zero
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
  ! result.  Also see AssignMatrix.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J      ! Subscripts for [XYZ]%Block

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%col%vec%template,y%col%vec%template) &
      & .or. .not. associated(x%row%vec%template,y%row%vec%template) &
      & .or. (x%col%instFirst .neqv. y%col%instFirst) &
      & .or. (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in AddMatrices" )
    call createEmptyMatrix ( z, 0, x%row%vec, x%col%vec )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        z%block(i,j) = x%block(i,j) + y%block(i,j)
      end do ! i = 1, x%row%nb
    end do ! j = 1, x%col%nb
  end function AddMatrices

  ! ------------------------------------------------  AddToMatrix  -----
  subroutine AddToMatrix ( X, Y ) ! X = X + Y
    type(Matrix_T), intent(inout) :: X
    type(Matrix_T), intent(in) :: Y

    integer :: I, J      ! Subscripts for [XYZ]%Block

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%col%vec%template,y%col%vec%template) &
      & .or. .not. associated(x%row%vec%template,y%row%vec%template) &
      & .or. (x%col%instFirst .neqv. y%col%instFirst) &
      & .or. (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in AddMatrices" )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        x%block(i,j) = x%block(i,j) + y%block(i,j)
      end do ! i = 1, x%row%nb
    end do ! j = 1, x%col%nb
  end subroutine AddToMatrix

  ! -----------------------------------------------  AssignMatrix  -----
  subroutine AssignMatrix ( Z, X )
  ! Destroy Z and then assign X to it, using pointer assignment for pointer
  ! components.  Notice that CopyMatrix does a deep copy.  If one has Z = X
  ! inside a loop, it is only necessary to destroy Z after the loop.
    type(Matrix_T), intent(inout) :: Z
    type(Matrix_T), intent(in) :: X
    call destroyMatrix ( z )
    z%name = x%name
    z%col = x%col
    z%row = x%row
    z%block => x%block
  end subroutine AssignMatrix

  ! -------------------------------------------  CholeskyFactor_1  -----
  subroutine CholeskyFactor_1 ( X, Z )
  ! Compute the Cholesky factor Z of the matrix X.  Z%M%Block can be
  ! associated with X%M%Block to save space.
    type(Matrix_SPD_T), intent(in) :: X ! Matrix to factor.
    type(Matrix_Cholesky_T), intent(inout) :: Z   ! Factored matrix.

    integer :: I, J, K                  ! Subscripts and loop inductors
    integer :: N, N0                    ! Columns(blocks), Columns(0-level)
    type(MatrixElement_T) :: S          ! Sum, to accumulate "inner product"

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%m%col%vec%template,z%m%col%vec%template) &
      & .or. .not. associated(x%m%row%vec%template,z%m%row%vec%template) &
      & .or. (x%m%col%instFirst .neqv. z%m%col%instFirst) &
      & .or. (x%m%row%instFirst .neqv. z%m%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Matrices in CholeskyFactor are not compatible" )
    n = x%m%row%nb
    do i = 1, n
      n0 = x%m%col%nelts(i)
      call createBlock ( s, n0, n0, M_Absent )
      do k = 1, i-1
        call multiplyMatrixBlocks ( z%m%block(k,i), z%m%block(k,i), s, &
          & update=.true. )
      end do ! k = 1, i-1
      call choleskyFactor ( z%m%block(i,i), s ) ! Factor block on diagonal
      do j = i+1, n
        call destroyBlock ( s )    ! Avoid a memory leak
        do k = 1, i-1
          call multiplyMatrixBlocks ( z%m%block(k,i), z%m%block(k,j), s, &
            & update=.true. )
        end do ! k = 1, i-1
        call solveCholesky ( z%m%block(i,i), z%m%block(i,j), s, transpose=.true. )
      end do ! j = 1, n
      call destroyBlock( s )       ! Avoid a memory leak
    end do ! i = 1, n
  end subroutine CholeskyFactor_1

  ! ------------------------------------------------  ClearRows_1  -----
  subroutine ClearRows_1 ( X )
  ! Clear the rows of X for which the mask in X's row-defining vector
  ! has nonzero bits.
    type(Matrix_T), intent(inout) :: X
    integer :: I, J                ! Subscripts and row indices
    integer :: NI, NQ              ! Instance and quantity indices
    do i = 1, x%row%nb
      ni = x%row%inst(i)
      nq = x%row%quant(i)
      if ( associated(x%row%vec%quantities(nq)%mask) ) then
        do j = 1, x%col%nb
          call clearRows ( x%block(i,j), x%row%vec%quantities(nq)%mask(:,ni) )
        end do ! j = 1, x%col%nb
      end if
    end do ! i = 1, x%row%nb
  end subroutine ClearRows_1

  ! ----------------------------------------------  ColumnScale_1  -----
  subroutine ColumnScale_1 ( X, V, NEWX ) ! Z = X V where V is a diagonal
  !                                matrix represented by a vector and Z is X
  !                                or NEWX.
    type (Matrix_T), intent(inout), target :: X
    type (Vector_T), intent(in) :: V
    type (Matrix_T), intent(out), target, optional :: NEWX 

    integer :: I, J      ! Subscripts for [XZ]%Block

    if ( present(newx) ) then
      call createEmptyMatrix ( newx, 0, x%row%vec, x%col%vec, &
        & x%row%instFirst, x%col%instFirst )
      do j = 1, x%col%nb
        do i = 1, x%row%nb
        end do ! i = x%row%nb
          call ColumnScale ( x%block(i,j), &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
            & newx%block(i,j) )
      end do ! j = x%col%nb
    else
      do j = 1, x%col%nb
        do i = 1, x%row%nb
        end do ! i = x%row%nb
          call ColumnScale ( x%block(i,j), &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)) )
      end do ! j = x%col%nb
    end if
  end subroutine ColumnScale_1

  ! -------------------------------------------------  CopyMatrix  -----
  subroutine CopyMatrix ( Z, X )        ! Destroy Z, then deep Z = X except
  !                                       the name of Z isn't changed.
    type(matrix_T), intent(inout) :: Z
    type(matrix_T), intent(in) :: X
    integer :: I, J ! Subscripts and loop inductors
    call destroyMatrix ( z )
    call copyRCInfo ( z%col, x%col )
    call copyRCInfo ( z%row, x%row )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        z%block(i,j) = x%block(i,j)
      end do ! i = 1, x%row%nb
    end do ! j = 1, x%col%nb
  end subroutine CopyMatrix

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
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Z%Block in CreateEmptyMatrix" )

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
  ! Destroy a matrix -- deallocate its pointer components, don't change the
  ! name
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
    call destroyRCInfo ( a%row )
    call destroyRCInfo ( a%col )
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
  ! result.  Also see AssignMatrix.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K             ! Subscripts for [XYZ]%Block

    ! Check that matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( .not. associated(x%row%vec%template,y%row%vec%template) .or. &
      & (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in MultiplyMatrices" )
    call createEmptyMatrix ( z, 0, x%col%vec, y%col%vec )
    do j = 1, y%col%nb
      do i = 1, x%col%nb
        call multiplyMatrixBlocks ( x%block(1,i), y%block(1,j), z%block(i,j) )
        do k = 2, x%row%nb
          call multiplyMatrixBlocks ( &
            & x%block(k,i), y%block(k,j), z%block(i,j), update=.true. )
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

    logical :: DO_UPDATE      ! Tells MatrixModule_0 % multiplyMatrixVector
    !                           whether to clear an element of Z, or add to it
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
        call multiplyMatrixVector ( a%block(i,j), &
          & v%quantities(m)%values(:,n), z%quantities(k)%values(:,l), &
          & do_update )
        do_update = .true.
      end do ! i = 1, a%row%nb
    end do ! j = 1, a%col%nb
  end subroutine MultiplyMatrixVector_1

  ! ----------------------------------  MultiplyMatrixVectorSPD_1  -----
  subroutine MultiplyMatrixVectorSPD_1 ( A, V, Z, UPDATE )
  ! Z = A V if UPDATE is absent or false.
  ! Z = Z + A V is UPDATE is present and true.
  ! Remember that for SPD, only the upper triangle is stored, so we need
  ! Z = A^T V + A V except don't do the diagonal twice.
    type(Matrix_SPD_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: Z
    logical, optional, intent(in) :: UPDATE

    integer :: I, J           ! Subscripts and loop inductors
    integer :: K, L, M, N     ! Subscripts

    call MultiplyMatrixVector ( a%m, v, z, update ) ! A^T V
    do j = 1, a%m%col%nb
      k = a%m%col%quant(j)
      l = a%m%col%inst(j)
      do i = 1, j
        m = a%m%row%quant(i)
        n = a%m%row%inst(i)
        call multiplyMatrixVectorNoT ( a%m%block(i,j), &
          & v%quantities(k)%values(:,l), z%quantities(m)%values(:,n), &
          & update=.true., doDiag = i /= j )
      end do ! i = 1, a%m%row%nb
    end do ! j = 1, a%m%col%nb
  end subroutine MultiplyMatrixVectorSPD_1

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
  ! result.  Also see AssignMatrix.
  ! !!!!! ===== END NOTE ===== !!!!! 

    logical :: DO_UPDATE
    integer :: I, J, K             ! Subscripts for [AZ]%Block
    integer, dimension(:), pointer :: MI, MJ ! Masks for columns I, J if any
    logical :: MY_UPDATE

    ! Compute Z = A^T A or Z = Z + A^T A
    my_update = .false.
    if ( present(update) ) my_update = update
    if ( .not. my_update ) &
      & call createEmptyMatrix ( z%m, 0, a%col%vec, a%col%vec )
    do j = 1, a%col%nb
      nullify ( mj )
      if ( associated(a%col%vec%quantities(a%col%quant(j))%mask) ) &
        mj => a%col%vec%quantities(a%col%quant(j))%mask(:,a%col%inst(j))
      do i = 1, j
        nullify ( mi )
        if ( associated(a%col%vec%quantities(a%col%quant(i))%mask) ) &
          mi => a%col%vec%quantities(a%col%quant(i))%mask(:,a%col%inst(i))
        do_update = my_update
        do k = 1, a%row%nb
          if ( associated(mi) ) then
            if ( associated(mj) ) then
              call multiplyMatrixBlocks ( &
                & a%block(k,i), a%block(k,j), z%m%block(i,j), update=do_update, &
                & xmask=mi, ymask=mj )
            else
              call multiplyMatrixBlocks ( &
                & a%block(k,i), a%block(k,j), z%m%block(i,j), update=do_update, &
                & xmask=mi )
            end if
          else
            if ( associated(mj) ) then
              call multiplyMatrixBlocks ( &
                & a%block(k,i), a%block(k,j), z%m%block(i,j), update=do_update, &
                & ymask=mj )
            else
              call multiplyMatrixBlocks ( &
                & a%block(k,i), a%block(k,j), z%m%block(i,j), update=do_update )
            end if
          end if
          do_update = .true.
        end do ! k = 2, a%row%nb
      end do ! i = 1, a%col%nb
    end do ! j = 1, a%col%nb

    call multiplyMatrixVector ( a, rhs_in, rhs_out, my_update )
  end subroutine NormalEquations

  ! -------------------------------------------------  RowScale_1  -----
  subroutine RowScale_1 ( V, X, NEWX ) ! Z = V X where V is a diagonal
  !                                matrix represented by a vector and Z is X
  !                                or NEWX.
    type (Vector_T), intent(in) :: V
    type (Matrix_T), intent(inout), target :: X
    type (Matrix_T), intent(out), target, optional :: NEWX 

    integer :: I, J      ! Subscripts for [XZ]%Block

    if ( present(newx) ) then
      call createEmptyMatrix ( newx, 0, x%row%vec, x%col%vec, &
        & x%row%instFirst, x%col%instFirst )
      do j = 1, x%col%nb
        do i = 1, x%row%nb
        end do ! i = x%row%nb
          call RowScale ( &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
            & x%block(i,j), newx%block(i,j) )
      end do ! j = x%col%nb
    else
      do j = 1, x%col%nb
        do i = 1, x%row%nb
        end do ! i = x%row%nb
          call RowScale ( &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
            & x%block(i,j) )
      end do ! j = x%col%nb
    end if
  end subroutine RowScale_1

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
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "X and RHS not compatible in SolveCholesky_1" )
      my_rhs => rhs
    end if
    if ( .not. associated( z%m%col%vec%template, my_rhs%template ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
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

  ! -------------------------------------------------  CopyRCInfo  -----
  subroutine CopyRCInfo ( A, B )
    type(RC_info), intent(inout) :: A
    type(RC_info), intent(in) :: B
    call destroyRCInfo ( a )
    a%vec => b%vec
    a%nb = b%nb
    call allocate_test ( a%nelts, size(b%nelts), "a%nelts in CopyRCInfo", &
      & moduleName )
    a%nelts = b%nelts
    call allocate_test ( a%inst, size(b%inst), "a%inst in CopyRCInfo", &
      & moduleName )
    a%inst = b%inst
    call allocate_test ( a%quant, size(b%quant), "a%quant in CopyRCInfo", &
      & moduleName )
    a%quant = b%quant
  end subroutine CopyRCInfo

  ! ----------------------------------------------  DestroyRCInfo  -----
  subroutine DestroyRCInfo ( RC )
    type(RC_Info), intent(inout) :: RC
    call deallocate_test ( rc%nelts, "rc%nelts in DestroyRCInfo", moduleName )
    call deallocate_test ( rc%inst, "rc%inst in DestroyRCInfo", moduleName )
    call deallocate_test ( rc%quant, "rc%quant in DestroyRCInfo", moduleName )
  end subroutine DestroyRCInfo

  ! -------------------------------------------------  DumpMatrix  -----
  subroutine Dump_Matrix ( Matrix, Name, Details )
    type(Matrix_T), intent(in) :: Matrix
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Details   ! Print details, default 1
    !  <= Zero => no details, == One => Details of matrix but not its blocks,
    !  >One => Details of the blocks, too.

    integer :: I, J                ! Subscripts, loop inductors
    integer :: MY_DETAILS          ! True if DETAILS is absent, else DETAILS

    my_details = 1
    if ( present(details) ) my_details = details
    if ( present(name) ) call output ( name, advance='yes' )
    call dump_rc ( matrix%row, 'row', my_details>0 )
    call dump_rc ( matrix%col, 'column', my_details>0 )
    do j = 1, matrix%col%nb
      do i = 1, matrix%row%nb
        call output ( 'Block at row ' )
        call output ( i )
        call output ( ' and column ' )
        call output ( j, advance='yes' )
        call dump ( matrix%block(i,j), details=my_details>1 )
      end do
    end do
  end subroutine Dump_Matrix

  subroutine Dump_RC ( RC, R_or_C, Details )
    type(rc_info), intent(in) :: RC
    character(len=*), intent(in) :: R_or_C
    logical, intent(in) :: Details
    call output ( 'Number of ' )
    call output ( r_or_c )
    call output ( ' blocks = ' )
    call output ( rc%nb )
    call output ( 'Vector that defines ' )
    call output ( r_or_c )
    call output ( 's' )
    if ( rc%vec%name == 0 ) then
      call output ( ' has no name', advance='yes' )
    else
      call output ( ': ' )
      call display_string ( rc%vec%name, advance='yes' )
    end if
    call output ( 'Order of blocks is ' )
    if ( rc%instFirst ) then
      call output ( 'instance, then quantity', advance='yes' )
    else
      call output ( 'quantity, then instance', advance='yes' )
    end if
    if ( details ) then
      call output ( 'Numbers of ' )
      call output ( r_or_c )
      call output ( 's in each block: ', advance='yes' )
      call dump ( rc%nelts )
      call output ( 'Instance indices for blocks in the' )
      call output ( r_or_c )
      call output ( 's:', advance='yes' )
      call dump ( rc%inst )
      call output ( 'Quantity indices for blocks in the' )
      call output ( r_or_c )
      call output ( 's:', advance='yes' )
      call dump ( rc%quant )
    end if
  end subroutine Dump_RC
end module MatrixModule_1

! $Log$
! Revision 2.4  2000/11/23 01:09:46  vsnyder
! Add provision to ignore columns during matrix-matrix multiply, finish DUMP.
!
! Revision 2.3  2000/11/15 00:18:26  vsnyder
! Added assignment(=) interface, row scale, column scale
!
! Revision 2.2  2000/11/10 00:28:13  vsnyder
! Added multiply untransposed matrix * vector
!
! Revision 2.1  2000/11/09 01:23:23  vsnyder
! Initial entry -- still under construction
!
