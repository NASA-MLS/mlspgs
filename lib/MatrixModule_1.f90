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
    & ColumnScale, Col_L1, CopyBlock, CreateBlock, DestroyBlock, Dump, &
    & GetDiagonal, GetMatrixElement, GetVectorFromColumn, M_Absent, &
    & M_Banded, M_Full, MatrixElement_T, MaxAbsVal, MinDiag, &
    & MultiplyMatrixBlocks, MultiplyMatrixVector, MultiplyMatrixVectorNoT, &
    & operator(+), operator(.TX.), RowScale, ScaleBlock, SolveCholesky, &
    & UpdateDiagonal
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: BLANKS, OUTPUT
  use String_Table, only: Display_String
  use VectorsModule, only: CloneVector, Vector_T

  implicit NONE
  private
  public :: AddMatrices, AddToMatrixDatabase, AddToMatrix
  public :: Assignment(=), CholeskyFactor, CholeskyFactor_1
  public :: ClearMatrix, ClearRows, ClearRows_1, ColumnScale, ColumnScale_1
  public :: CopyMatrix, CopyMatrixValue, CreateBlock, CreateBlock_1, CreateEmptyMatrix
  public :: DestroyBlock, DestroyBlock_1, DestroyMatrix
  public :: DestroyMatrixInDatabase, DestroyMatrixDatabase, Dump, Dump_Linf
  public :: Dump_Struct, FillExtraCol, FillExtraRow, FindBlock, GetDiagonal
  public :: GetDiagonal_1, GetFromMatrixDatabase, GetKindFromMatrixDatabase
  public :: GetMatrixElement, GetMatrixElement_1, GetVectorFromColumn
  public :: GetVectorFromColumn_1, InvertCholesky
  public :: K_Cholesky, K_Empty, K_Kronecker, K_Plain, K_SPD
! public :: LevenbergUpdateCholesky
  public :: Matrix_T, Matrix_Cholesky_T, Matrix_Database_T, Matrix_Kronecker_T
  public :: Matrix_SPD_T, MaxAbsVal, MaxAbsVal_1, MaxL1
  public :: MinDiag, MinDiag_Cholesky, MinDiag_SPD, MultiplyMatrices
  public :: MultiplyMatrixVector, MultiplyMatrixVector_1
  public :: MultiplyMatrixVectorNoT, MultiplyMatrixVectorNoT_1
  public :: MultiplyMatrixVectorSPD_1
  public :: Negate, Negate_1
  public :: NewMultiplyMatrixVector, NormalEquations, operator(.TX.)
  public :: operator(+), RC_Info, RowScale, RowScale_1, ScaleMatrix
  public :: SolveCholesky, SolveCholesky_1
  public :: UpdateDiagonal, UpdateDiagonal_1, UpdateDiagonalVec_1

! =====     Defined Operators and Generic Identifiers     ==============

  interface AddToMatrixDatabase
    module procedure AddMatrixToDatabase, AddCholeskyToDatabase
    module procedure AddKroneckerToDatabase, AddSPDToDatabase
  end interface

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

  interface CreateBlock
    module procedure CreateBlock_1
  end interface

  interface DestroyBlock
    module procedure DestroyBlock_1
  end interface

  interface DestroyMatrix
    module procedure DestroyMatrix
    module procedure DestroyMatrixInDatabase
  end interface

  interface Dump
    module procedure Dump_Matrix, Dump_Matrix_Database
  end interface

  interface GetDiagonal
    module procedure GetDiagonal_1
  end interface

  interface GetFromMatrixDatabase
    module procedure GetMatrixFromDatabase, GetCholeskyFromDatabase
    module procedure GetKroneckerFromDatabase, GetSPDFromDatabase
  end interface

  interface GetVectorFromColumn
    module procedure GetVectorFromColumn_1
  end interface

  interface MaxAbsVal
    module procedure MaxAbsVal_1
  end interface

  interface MinDiag
    module procedure MinDiag_Cholesky, MinDiag_SPD
  end interface

  interface GetMatrixElement
    module procedure GetMatrixElement_1
  end interface

  interface MultiplyMatrixVector   ! A^T V
    module procedure MultiplyMatrixVector_1, MultiplyMatrixVectorSPD_1
  end interface

  interface MultiplyMatrixVectorNoT   ! A V
    module procedure MultiplyMatrixVectorNoT_1
  end interface

  interface Negate
    module procedure Negate_1
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
    module procedure UpdateDiagonal_1, UpdateDiagonalVec_1
  end interface

  integer, parameter :: K_Empty = 0                    ! Empty database element
  integer, parameter :: K_Cholesky = k_Empty + 1
  integer, parameter :: K_Kronecker = k_Cholesky + 1
  integer, parameter :: K_Plain = k_Kronecker  + 1
  integer, parameter :: K_SPD = k_Plain + 1

  type RC_Info
  ! Information about the row or column of a matrix
    type(Vector_T) :: Vec               ! Vector used to define the row
      ! or column space of the matrix, if any.
    integer :: NB = 0              ! Number of blocks of rows or columns
    logical :: Extra = .false.     ! There is an extra row or column that is
      ! not accounted for by vec.
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

  type Matrix_Database_T
    private
    type(Matrix_T), pointer :: Matrix => NULL()
    type(Matrix_Cholesky_T), pointer :: Cholesky => NULL()
    type(Matrix_Kronecker_T), pointer :: Kronecker => NULL()
    type(Matrix_SPD_T), pointer :: SPD => NULL()
  end type Matrix_Database_T

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------  AddCholeskyToDatabase  -----
  integer function AddCholeskyToDatabase ( Database, CholeskyItem )
  ! Add a Cholesky factor matrix to the matrix database.
    type(matrix_Database_T), dimension(:), pointer :: Database
    type(matrix_cholesky_T), intent(in) :: CholeskyItem

    type(matrix_Database_T) :: Item
    integer :: Status

    allocate ( item%cholesky, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Matrix database element" )
    item%cholesky = choleskyItem
    addCholeskyToDatabase = addItemToMatrixDatabase ( database, item )
  end function AddCholeskyToDatabase

  ! -------------------------------------  AddKroneckerToDatabase  -----
  integer function AddKroneckerToDatabase ( Database, KroneckerItem )
  ! Add a Kronecker product matrix to the matrix database.
    type(matrix_Database_T), dimension(:), pointer :: Database
    type(matrix_kronecker_T), intent(in) :: KroneckerItem

    type(matrix_Database_T) :: Item
    integer :: Status

    allocate ( item%kronecker, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Matrix database element" )
    item%kronecker = kroneckerItem
    addKroneckerToDatabase = addItemToMatrixDatabase ( database, item )
  end function AddKroneckerToDatabase

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

    ! Check that the matrices are compatible.  We also need to check
    ! Nb, because the matrices may have been created with an extra
    ! row or column.
    if ( x%col%vec%template%id /= y%col%vec%template%id &
      & .or. x%row%vec%template%id /= y%row%vec%template%id &
      & .or. x%col%nb /= y%col%nb .or. x%row%nb /= y%row%nb &
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

  ! ----------------------------------------  AddMatrixToDatabase  -----
  integer function AddMatrixToDatabase ( Database, MatrixItem )
  ! Add a matrix of unspecified structure to the matrix database
    type(matrix_Database_T), dimension(:), pointer :: Database
    type(matrix_T), intent(in) :: MatrixItem

    type(matrix_Database_T) :: Item
    integer :: Status

    allocate ( item%matrix, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Matrix database element" )
    item%matrix = matrixItem
    addMatrixToDatabase = addItemToMatrixDatabase ( database, item )
  end function AddMatrixToDatabase

  ! -------------------------------------------  AddSPDToDatabase  -----
  integer function AddSPDToDatabase ( Database, SPDItem )
  ! Add a symmetric-positive-definite matrix to the matrix database
    type(matrix_Database_T), dimension(:), pointer :: Database
    type(matrix_spd_T), intent(in) :: SPDItem

    type(matrix_Database_T) :: Item
    integer :: Status

    allocate ( item%spd, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Matrix database element" )
    item%spd = spdItem
    addSPDToDatabase = addItemToMatrixDatabase ( database, item )
  end function AddSPDToDatabase

  ! ------------------------------------------------  AddToMatrix  -----
  subroutine AddToMatrix ( X, Y ) ! X = X + Y
    type(Matrix_T), intent(inout) :: X
    type(Matrix_T), intent(in) :: Y

    integer :: I, J      ! Subscripts for [XYZ]%Block

    ! Check that the matrices are compatible.  We need to check
    ! Nb, because the matrices may have been created with an extra
    ! row or column.
    if ( x%col%vec%template%id /= y%col%vec%template%id &
      & .or. x%row%vec%template%id /= y%row%vec%template%id &
      & .or. x%col%nb /= y%col%nb .or. x%row%nb /= y%row%nb &
      & .or. (x%col%instFirst .neqv. y%col%instFirst) &
      & .or. (x%row%instFirst .neqv. y%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in AddMatrices" )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        x%block(i,j) = x%block(i,j) + y%block(i,j) ! Defined =, +
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
  subroutine CholeskyFactor_1 ( Z, X )
  ! Compute the Cholesky factor Z of the matrix X.  Z%M%Block can be
  ! associated with X%M%Block to save space.
    type(Matrix_Cholesky_T), intent(inout) :: Z   ! Factored matrix.
    type(Matrix_SPD_T), intent(in) :: X ! Matrix to factor.

    integer :: I, J, K                  ! Subscripts and loop inductors
    integer :: N                        ! Columns(blocks)
    type(MatrixElement_T) :: S          ! Sum, to accumulate "inner product"

    ! Check that the matrices are compatible.  We don't need to check
    ! Nelts or Nb, because these are deduced from Vec.
    if ( x%m%col%vec%template%id /= z%m%col%vec%template%id &
      & .or. x%m%row%vec%template%id /= z%m%row%vec%template%id &
      & .or. (x%m%col%instFirst .neqv. z%m%col%instFirst) &
      & .or. (x%m%row%instFirst .neqv. z%m%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Matrices in CholeskyFactor are not compatible" )

!{Suppose $A$ is symmetric and positive definite.  Regard $A = (A_{ij})$
!  and its upper-triangular Cholesky factor $G = (G_{ij})$ as $N \times N$
!  matrices with square diagonal blocks.  By equating $(i,j)$ blocks in the
!  equation $A = G^T G$ with $i \leq j$, it follows that
!  \begin{equation*}
!   A_{ij} = \sum_{k=1}^i G_{ki}^T G_{kj}
!  \end{equation*}
!  Define
!  \begin{equation*}
!  S = A_{ij} - \sum_{k=1}^{i-1} G_{ki}^T G_{kj}
!    = A_{ij} - \sum_{k=1}^i G_{ki}^T G_{kj} + G_{ii}^T G_{ij}
!    = G_{ii}^T G_{ij}
!  \end{equation*}
!  
!  Rearranging, we have the following algorithm:
!  
!  {\bf do} $i = 1, N$\\
!  \hspace*{0.25in} $S= A_{ii} - \sum_{k=1}^{i-1} G_{ki}^T G_{ki}$\\
!  \hspace*{0.25in} $G_{ii} = $ Cholesky factor of $S$\\
!  \hspace*{0.25in} {\bf do} $j = i+1, N$\\
!  \hspace*{0.5in}    $S= A_{ij} - \sum_{k=1}^{i-1} G_{ki}^T G_{kj}$\\
!  \hspace*{0.5in}    solve $G_{ii}^T G_{ij} = S$ for $G_{ij}$\\
!  \hspace*{0.25in} {\bf end do} ! j\\
!  {\bf end do} ! i

    n = x%m%row%nb
    ! Handle the first row specially, to avoid a copy followed by no-dot-product
    !{ $Z_{11}^T Z_{11} = X_{11}$
    call choleskyFactor ( z%m%block(1,1), x%m%block(1,1) )
    do j = 2, n
      !{ Solve $Z_{11}^T Z_{1j} = X_{1j}$ for $Z_{1j}$
      call solveCholesky ( z%m%block(1,1), z%m%block(1,j), x%m%block(1,j), &
        & transpose=.true. )
    end do
    do i = 2, n
      call copyBlock ( s, x%m%block(i,i) )        ! Destroy s, then s := z...
      do k = 1, i-1
        !{ $S = X_{ii} - \sum_{k=1}^{i-1} Z_{ki}^T Z_{ki}$
        call multiplyMatrixBlocks ( z%m%block(k,i), z%m%block(k,i), s, &
          & update=.true., subtract=.true. )
      end do ! k = 1, i-1
      !{ $Z_{ii}^T Z_{ii} = S$
      call choleskyFactor ( z%m%block(i,i), s )   ! z%m%block(i,i) = factor of s
      do j = i+1, n
        if ( i == 1 ) then                        ! Avoid a copy
        else
          call copyBlock ( s, x%m%block(i,j) )    ! Destroy s, then s := x...
          do k = 1, i-1
            !{ $S = X_{ij} - \sum_{k=1}^{i-1} Z_{ki}^T Z_{kj}$
            call multiplyMatrixBlocks ( z%m%block(k,i), z%m%block(k,j), s, &
              & update=.true., subtract=.true. )
          end do ! k = 1, i-1
          !{ Solve $Z_{ii}^T Z_{ij} = S$ for $Z_{ij}$
          call solveCholesky ( z%m%block(i,i), z%m%block(i,j), s, &
            & transpose=.true. )
        end if
      end do ! j = 1, n
      call destroyBlock( s )                      ! Avoid a memory leak
    end do ! i = 1, n
  end subroutine CholeskyFactor_1

  ! ------------------------------------------------  ClearMatrix  -----
  subroutine ClearMatrix ( X )     ! Delete all of the blocks, but keep the
                                   ! structural information
    type(matrix_T), intent(inout) :: X
    integer :: I, J                ! Subscripts and row indices
    if (associated(x%block)) then
      do i = 1, x%row%nb
        do j = 1, x%col%nb
          call destroyBlock ( x%block(i,j) )
        end do ! j
      end do ! i
    end if
  end subroutine ClearMatrix

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
  !                                Any "extra" column is not scaled (but an
  !                                extra row is).
    type (Matrix_T), intent(inout), target :: X
    type (Vector_T), intent(in) :: V
    type (Matrix_T), intent(out), target, optional :: NEWX 

    integer :: I, J      ! Subscripts for [XZ]%Block
    integer :: NC, NR    ! Number of columns and rows

    nc = x%col%nb
    if ( x%col%extra ) nc = nc - 1
    nr = x%row%nb
    if ( present(newx) ) then
      call createEmptyMatrix ( newx, 0, x%row%vec, x%col%vec, &
        & x%row%instFirst, x%col%instFirst )
      do j = 1, nc
        do i = 1, nr
          call ColumnScale ( x%block(i,j), &
            & v%quantities(x%row%quant(j))%values(:,x%row%inst(j)), &
            & newx%block(i,j) )
        end do ! i = nr
      end do ! j = nc
    else
      do j = 1, nc
        do i = 1, nr
          call ColumnScale ( x%block(i,j), &
            & v%quantities(x%row%quant(j))%values(:,x%row%inst(j)) )
        end do ! i = nr
      end do ! j = nc
    end if
  end subroutine ColumnScale_1

  ! -------------------------------------------------  CopyMatrix  -----
  subroutine CopyMatrix ( Z, X )        ! Destroy Z, then deep Z = X except
  !                                       the name of Z isn't changed.
    type(matrix_T), intent(inout) :: Z
    type(matrix_T), intent(in) :: X
    integer :: I, J      ! Subscripts and loop inductors
    integer :: Status    ! From allocate
    call destroyMatrix ( z )
    call copyRCInfo ( z%col, x%col )
    call copyRCInfo ( z%row, x%row )
    allocate ( z%block(z%row%nb,z%col%nb), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Z%Block in CreateEmptyMatrix" )
    do j = 1, x%col%nb
      do i = 1, x%row%nb
        call createBlock ( z, i, j, m_absent ) ! Create block w/correct size
        call copyBlock ( z%block(i,j), x%block(i,j) )
      end do ! i = 1, x%row%nb
    end do ! j = 1, x%col%nb
  end subroutine CopyMatrix

  ! --------------------------------------------  CopyMatrixValue  -----
  subroutine CopyMatrixValue ( Z, X )   ! Copy the elements of X to Z.
  ! Z and X must have the same template, but it's OK if they don't both
  ! have row%extra or col%extra.  If Z has extra and X does not, Z's
  ! extra is deleted.
    type(matrix_T), intent(inout) :: Z
    type(matrix_T), intent(in) :: X
    integer :: I, J ! Subscripts and loop inductors
    if ( x%col%vec%template%id /= z%col%vec%template%id &
      & .or. x%row%vec%template%id /= z%row%vec%template%id &
      & .or. (x%col%instFirst .neqv. z%col%instFirst) &
      & .or. (x%row%instFirst .neqv. z%row%instFirst) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Incompatible arrays in CopyMatrixValue" )
    do i = 1, min(x%row%nb,z%row%nb)
      do j = 1, min(x%col%nb,z%col%nb)
        call copyBlock ( z%block(i,j), x%block(i,j) )
      end do ! j
    end do ! i
    if ( z%col%nb > x%col%nb ) then
      do i = 1, z%row%nb
        call destroyBlock ( z%block(i, z%col%nb) )
      end do ! i
    end if
    if ( z%row%nb > x%row%nb ) then
      do j = 1, z%col%nb
        call destroyBlock ( z%block(z%row%nb, j) )
      end do ! j
    end if
  end subroutine CopyMatrixValue

  ! ----------------------------------------------  CreateBlock_1  -----
  subroutine CreateBlock_1 ( Z, RowNum, ColNum, Kind, NumNonzeros, BandHeight )
  ! Create the matrix block Z%Block(RowNum,ColNum), which sprang into
  ! existence with kind M_Absent.  Create it with the specified Kind.
  ! See MatrixModule_0 for a list of the kinds.  If the Kind is
  ! M_Banded or M_ColumnSparse, the number of nonzeroes is needed.
  ! If Kind == M_Banded and BandHeight is present, the band height is
  !   assumed to be uniform, and the R1 and R2 components are filled to
  !   reflect that assumption.
    integer, intent(in), optional :: BandHeight
    type(matrix_T), intent(inout) :: Z       ! The matrix having the block
    integer, intent(in) :: RowNum, ColNum    ! Row and column of the block
    integer, intent(in) :: Kind         ! Kind of block, see MatrixModule_0
    integer, intent(in), optional :: NumNonzeros  ! Number of nonzeros
    call createBlock ( z%block(rowNum,colNum), &
      & z%row%nelts(rowNum), z%col%nelts(colNum), kind, numNonzeros, &
      & bandHeight=bandHeight )
  end subroutine CreateBlock_1

  ! ------------------------------------------  CreateEmptyMatrix  -----
  subroutine CreateEmptyMatrix ( Z, Name, Row, Col &
    &,                           Row_Quan_First, Col_Quan_First &
    &,                           Extra_Row, Extra_Col )
    type(Matrix_T), intent(out) :: Z    ! The matrix to create
    integer, intent(in) :: Name         ! Sub-rosa index of its name, or zero
    type(Vector_T), intent(in) :: Row       ! Vector used to define the row
      !                                   space of the matrix.
    type(Vector_T), intent(in) :: Col       ! Vector used to define the column
      !                                   space of the matrix.
    logical, intent(in), optional :: Row_Quan_First    ! True (default false)
      ! means the quantity is the major order of the rows of blocks and the
      ! instance is the minor order.
    logical, intent(in), optional :: Col_Quan_First    ! True (default false)
      ! means the quantity is the major order of the columns of blocks and the
      ! instance is the minor order.
    logical, intent(in), optional :: Extra_Row         ! Allocate one extra
      ! row, beyond what's specified by Col, if Extra_Row is present and true.
    logical, intent(in), optional :: Extra_Col         ! Allocate one extra
      ! column, beyond what's specified by Row, if Extra_Col is present and true.

    integer :: I, J      ! Subscripts, loop inductors
    integer :: STATUS    ! From ALLOCATE

    z%name = name
    call defineInfo ( z%row, row, row_Quan_First, extra_Row )
    call defineInfo ( z%col, col, col_Quan_First, extra_Col )
    allocate ( z%block(z%row%nb,z%col%nb), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Z%Block in CreateEmptyMatrix" )
    do i = 1, z%row%nb ! Now create absent blocks with the correct sizes
      do j = 1, z%col%nb
        call createBlock ( z, i, j, m_absent )
      end do
    end do

  contains
    subroutine DefineInfo ( RC, Vec, QuanFirst, extra )
      type(RC_Info), intent(out) :: RC
      type(Vector_T), intent(in) :: Vec
      logical, intent(in), optional :: QuanFirst
      logical, intent(in), optional :: Extra

      integer :: I, J, N      ! Subscripts or loop inductors
      logical :: NEW          ! Was an instance seen?

      rc%extra = .false.
      if ( present(extra) ) rc%extra = extra
      rc%vec = vec
      rc%instFirst = .true.
      if ( present(quanFirst) ) rc%instFirst = .not. quanFirst
      rc%nb = vec%template%totalInstances
      if ( rc%extra ) rc%nb = rc%nb + 1
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
      if ( rc%extra ) then
        rc%nelts(rc%nb) = 1
        rc%inst(rc%nb) = 0
        rc%quant(rc%nb) = 0
      end if
    end subroutine DefineInfo
  end subroutine CreateEmptyMatrix

  ! ---------------------------------------------  DestroyBlock_1  -----
  subroutine DestroyBlock_1 ( A )
  ! Destroy the "block" component of a matrix.  This leaves its structure
  ! intact, but releases the space for its values.  This is useful when
  ! forming normal equations little-by-little.
    type(matrix_T), intent(inout) :: A

    integer :: STATUS              ! From deallocate

    call clearMatrix ( a )
    if ( associated(a%block) ) then
      deallocate ( a%block, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_DeAllocate // "A%Block in DestroyMatrix" )
    end if
  end subroutine DestroyBlock_1

  ! ----------------------------------------------  DestroyMatrix  -----
  subroutine DestroyMatrix ( A )
  ! Destroy a matrix -- deallocate its pointer components, don't change the
  ! name
    type(matrix_T), intent(inout) :: A

    call destroyBlock ( a )
    call destroyRCInfo ( a%row )
    call destroyRCInfo ( a%col )
  end subroutine DestroyMatrix

  ! ------------------------------------  DestroyMatrixInDatabase  -----
  subroutine DestroyMatrixInDatabase ( Database )
  ! Destroy a matrix in the database -- deallocate its pointer components,
  ! don't change the name
    type(matrix_database_T), intent(inout) :: Database
    integer :: Status  ! from deallocate

    if ( associated(database%matrix) ) then
      call destroyMatrix ( database%matrix )
      deallocate ( database%matrix, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, moduleName, &
        & MLSMSG_Deallocate // 'PlainMatrix' )
    else if ( associated(database%cholesky) ) then
      call destroyMatrix ( database%cholesky%m )
      deallocate ( database%cholesky, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, moduleName, &
        & MLSMSG_Deallocate // 'CholeskyMatrix' )
    else if ( associated(database%kronecker) ) then
      call destroyMatrix ( database%kronecker%m )
      deallocate ( database%kronecker, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, moduleName, &
        & MLSMSG_Deallocate // 'KroneckerMatrix' )
    else if ( associated(database%spd) ) then
      call destroyMatrix ( database%spd%m )
      deallocate ( database%spd, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, moduleName, &
        & MLSMSG_Deallocate // 'SPDMatrix' )
    end if
  end subroutine DestroyMatrixInDatabase

  ! --------------------------------------  DestroyMatrixDatabase  -----
  subroutine DestroyMatrixDatabase ( D )
  ! Destroy every matrix in the database D, then destroy the database.
    type(matrix_database_T), dimension(:), pointer :: D

    integer :: I, Status

    if ( .not. associated(d) ) return
    do i = 1, size(d)
      if ( associated(d(i)%matrix) ) then
        call destroyMatrix ( d(i)%matrix )
        call deallocateMatrix ( d(i)%matrix )
      end if
      if ( associated(d(i)%cholesky) ) then
        call destroyMatrix ( d(i)%cholesky%m )
        call destroyMatrix ( d(i)%cholesky%m )
      end if
      if ( associated(d(i)%kronecker) ) then
        call destroyMatrix ( d(i)%kronecker%m )
        call destroyMatrix ( d(i)%kronecker%m )
      end if
      if ( associated(d(i)%spd) ) then
        call destroyMatrix ( d(i)%spd%m )
        call destroyMatrix ( d(i)%spd%m )
      end if
    end do
    deallocate ( d, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_DeAllocate // "D in DestroyMatrixDatabase" )
  contains
    subroutine DeallocateMatrix ( M )
      type(matrix_t), pointer :: M
      deallocate ( m, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_DeAllocate // "D%matrix in DestroyMatrixDatabase" )
    end subroutine DeallocateMatrix
  end subroutine DestroyMatrixDatabase

  ! -----------------------------------------------  FillExtraCol  -----
  subroutine FillExtraCol ( A, X, ROW )
  ! Fill the "extra" column of A (see type RC_Info) from the vector X.
  ! Assume it is full -- i.e. don't bother to sparsify it.  If ROW is
  ! specified, it refers to a block-row, and only that row of the extra
  ! column is filled.
    type(Matrix_T), intent(inout) :: A
    type(Vector_T), intent(in) :: X
    integer, intent(in), optional :: ROW

    integer :: I

    ! Check compatibility of X with the row template for A
    if ( a%row%vec%template%id /= x%template%id ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Vector incompatible with matrix in FillExtraCol" )
    if ( .not. a%col%extra ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No extra column to fill in FillExtraCol" )
    if ( present(row) ) then
      if ( row < 1 .or. row > a%row%nb ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, "Row number out-of-range in FillExtraCol" )
      if ( a%block(row,a%col%nb)%kind /= m_full ) &
        & call destroyBlock ( a%block(row,a%col%nb) )
      if ( .not. associated(a%block(row,a%col%nb)%values) ) &
        & call createBlock ( a%block(row,a%col%nb), a%row%nelts(row), 1, m_full )
      a%block(row,a%col%nb)%values(:,1) = &
          & x%quantities(a%row%quant(row))%values(:,a%row%inst(row))
    else
      do i = 1, a%row%nb
        call destroyBlock ( a%block(i,a%col%nb) )
        call createBlock ( a%block(i,a%col%nb), a%row%nelts(i), 1, m_full )
        if ( a%row%quant(i) /= 0 ) &
          & a%block(i,a%col%nb)%values(:,1) = &
            & x%quantities(a%row%quant(i))%values(:,a%row%inst(i))
      end do ! i
    end if
  end subroutine FillExtraCol

  ! -----------------------------------------------  FillExtraRow  -----
  subroutine FillExtraRow ( A, X )
  ! Fill the "extra" row of A (see type RC_Info) from the vector X.
  ! The extra column in the extra row is not filled.
    type(Matrix_T), intent(inout) :: A
    type(Vector_T), intent(in) :: X

    integer :: J, NB

    ! Check compatibility of X with the column template for A
    if ( a%col%vec%template%id /= x%template%id ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Vector incompatible with matrix in FillExtraRow" )
    if ( .not. a%row%extra ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No extra row to fill in FillExtraRow" )
    nb = a%col%nb
    if ( a%col%extra ) nb = nb - 1
    do j = 1, nb
      call destroyBlock ( a%block(a%row%nb,j) )
      call createBlock ( a%block(a%row%nb,j), 1, a%col%nelts(j), m_full )
      if ( a%col%quant(j) /= 0 ) &
        & a%block(a%row%nb,j)%values(1,:) = &
          & x%quantities(a%col%quant(j))%values(:,a%col%inst(j))
    end do ! i
  end subroutine FillExtraRow

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

  ! ------------------------------------  GetCholeskyFromDatabase  -----
  subroutine GetCholeskyFromDatabase ( DatabaseElement, Cholesky )
  ! Get a POINTER to a Cholesky object from DatabaseElement.
    type(matrix_Database_T), intent(in) :: DatabaseElement
    type(matrix_cholesky_T), pointer :: Cholesky
    cholesky => databaseElement%cholesky
  end subroutine GetCholeskyFromDatabase

  ! ----------------------------------------------  GetDiagonal_1  -----
  subroutine GetDiagonal_1 ( A, X )
  ! Get X from the diagonal of A.  Don't get the extra row or column.
  ! Destroy X and re-create it by cloning the row-defining vector of A.
    type(Matrix_SPD_T), intent(in) :: A
    type(vector_T), intent(inout) :: X

    integer :: I, N

    call cloneVector ( x, a%m%row%vec, vectorNameText='_x' )
    n = max(a%m%row%nb,a%m%col%nb)
    if ( a%m%row%extra .or. a%m%col%extra ) n = n - 1
    do i = 1, n
      call getDiagonal ( a%m%block(i,i), &
        & x%quantities(a%m%row%quant(i))%values(:,a%m%row%inst(i)) )
    end do
  end subroutine GetDiagonal_1

  ! ----------------------------------  GetKindFromMatrixDatabase  -----
  integer function GetKindFromMatrixDatabase ( DatabaseElement )
  ! Return the kind of a matrix database element
    type(matrix_Database_T), intent(in) :: DatabaseElement
    getKindFromMatrixDatabase = k_empty
    if ( associated(databaseElement%cholesky) ) &
      & getKindFromMatrixDatabase = k_cholesky
    if ( associated(databaseElement%kronecker) ) &
      & getKindFromMatrixDatabase = k_kronecker
    if ( associated(databaseElement%matrix) ) &
      & getKindFromMatrixDatabase = k_plain
    if ( associated(databaseElement%spd) ) getKindFromMatrixDatabase = k_spd
  end function GetKindFromMatrixDatabase

  ! -----------------------------------  GetKroneckerFromDatabase  -----
  subroutine GetKroneckerFromDatabase ( DatabaseElement, Kronecker )
  ! Get a POINTER to a Kronecker object from DatabaseElement.
    type(matrix_Database_T), intent(in) :: DatabaseElement
    type(matrix_kronecker_T), pointer :: Kronecker
    kronecker => databaseElement%kronecker
  end subroutine GetKroneckerFromDatabase

  ! -----------------------------------------  GetMatrixElement_1  -----
  real(r8) function GetMatrixElement_1 ( Matrix, Row, Col )
  ! Get the (row,col) element of Matrix
    type(matrix_T), intent(in) :: Matrix
    integer, intent(in) :: Row, Col

    integer :: I, J      ! Block of Matrix
    integer :: II, JJ    ! Element of matrix%block(i,j)
    integer :: MM, NN    ! Number of rows and columns in blocks before (i,j)

    if ( .not. associated(matrix%block) ) call MLSMessage ( MLSMSG_Error, &
      & "Block not associated in GetMatrixElement", moduleName )
    ii = 0
    jj = 0
    mm = 0
    nn = 0
    if ( row > 0 ) then
      do i = 1, matrix%row%nb
        if ( mm + matrix%block(i,1)%nrows >= row ) then
          ii = row - mm
          exit
        end if
        mm = mm + matrix%block(i,1)%nrows
      end do ! i
    end if
    if ( col > 0 ) then
      do j = 1, matrix%col%nb
        if ( nn + matrix%block(1,j)%ncols >= col ) then
          jj = col - nn
          exit
        end if
        nn = nn + matrix%block(1,j)%ncols
      end do
    end if
    if ( ii == 0 .or. jj == 0 ) call MLSMessage ( MLSMSG_Error, &
      & "Row or Column subscript out-of-bounds in GetMatrixElement", &
      & moduleName )
    getMatrixElement_1 = getMatrixElement ( matrix%block(i,j), ii, jj )
  end function GetMatrixElement_1

  ! --------------------------------------  GetMatrixFromDatabase  -----
  subroutine GetMatrixFromDatabase ( DatabaseElement, Matrix )
  ! Get a POINTER to a matrix object from DatabaseElement.
    type(matrix_Database_T), intent(in) :: DatabaseElement
    type(matrix_T), pointer :: Matrix
    matrix => databaseElement%matrix
  end subroutine GetMatrixFromDatabase

  ! -----------------------------------------  GetSPDFromDatabase  -----
  subroutine GetSPDFromDatabase ( DatabaseElement, SPD )
  ! Get a POINTER to a SPD object from DatabaseElement.
    type(matrix_Database_T), intent(in) :: DatabaseElement
    type(matrix_SPD_T), pointer :: SPD
    SPD => databaseElement%SPD
  end subroutine GetSPDFromDatabase

  ! --------------------------------------  GetVectorFromColumn_1  -----
  subroutine GetVectorFromColumn_1 ( Matrix, Column, Vector )
  ! Fill the Vector from the Column of the Matrix
    type(matrix_T), intent(in) :: Matrix
    integer, intent(in) :: Column
    type(vector_T), intent(inout) :: Vector  ! Must already be "created"

    integer :: Block          ! Which block of columns contains Column?
    integer :: ColInBlock     ! Which column in Block is column
    integer :: Ncols          ! How many columns in blocks < Block?
    integer :: Row            ! Which row is being extracted?

    if ( column < 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'In "GetVectorFromColumn", "Column" < 1' )
    if ( matrix%col%vec%template%id /= vector%template%id ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Vector incompatible with matrix in GetVectorFromColumn" )
    ncols = 0
    do block = 1, matrix%col%nb
      if ( ncols + matrix%col%nelts(block) >= column ) then
        colInBlock = column - ncols
        do row = 1, matrix%row%nb
          if ( row <= vector%template%totalinstances ) & ! Don't get the extra
            ! row if the vector doesn't have a place for it.
            & call getVectorFromColumn ( matrix%block(row,block), colInBlock, &
              & vector%quantities(matrix%row%quant(row))% &
              & values(:,matrix%row%inst(row)) )
        end do ! row = 1, matrix%row%nb
        return
      end if
      ncols = ncols + matrix%col%nelts(block)
    end do ! block = 1, matrix%col%nb
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'In "GetVectorFromColumn", "Column" is greater than number&
      & of columns in "Matrix"')
  end subroutine GetVectorFromColumn_1

  ! ---------------------------------------------  InvertCholesky  -----
  subroutine InvertCholesky ( U, B )
  ! Compute B = U^{-T} = L^{-1}, where U is the upper-triangular
  ! Cholesky factor of some matrix, i.e. A = U^T U.
    type(Matrix_Cholesky_T), intent(in) :: U
    type(Matrix_T), intent(inout) :: B  ! Assume B has been created

    integer :: I, J, K                  ! Subscripts and loop inductors

    do i = 1, u%m%row%nb
      do k = 1, i
        if ( i == k ) then ! start with identity
          call createBlock ( b%block(i,i), u%m%block(i,i)%nrows, &
            & u%m%block(i,i)%nrows, m_banded, u%m%block(i,i)%nrows )
          do j = 1, u%m%block(i,i)%nrows
            b%block(i,k)%r1(j) = j
            b%block(i,k)%r2(j) = j
            b%block(i,k)%values(j,1) = 1.0_r8
          end do ! j = 1, u%m%block(i,i)%nrows
        else ! start with zero
          call createBlock ( b%block(i,k), u%m%block(i,k)%nrows, &
            & u%m%block(i,k)%ncols, m_absent )
        end if
        do j = k, i-1
          call multiplyMatrixBlocks ( b%block(i,k), u%m%block(j,i), &
            & b%block(j,k), update=.true., subtract=.true. )
        end do ! j = k, i-1
        call solveCholesky ( u%m%block(i,i), b%block(i,k), transpose=.true. )
      end do ! k = 1, i
    end do ! i = 1, u%m%row%nb
  end subroutine InvertCholesky

  ! ------------------------------------  LevenbergUpdateCholesky  -----
! subroutine LevenbergUpdateCholesky ( Z, LAMBDA )
! ! Given a Cholesky factor Z of a matrix of normal equations A^T A,
! ! update the Cholesky factor to incorporate Levenberg-Marquardt
! ! stabilization that corresponds to augmenting A with LAMBDA I.
!   type(Matrix_T), intent(inout) :: Z
!   real(r8), intent(in) :: LAMBDA
! end subroutine LevenbergUpdateCholesky

  ! ------------------------------------------------  MaxAbsVal_1  -----
  real(r8) function MaxAbsVal_1 ( A )
  ! Return the magnitude of the element in A that has the largest magnitude.
  ! Don't include the extra row or column.
    type(Matrix_T), intent(in) :: A
    integer :: I, J, M, N
    maxAbsVal_1 = 0.0_r8
    m = a%row%nb
    if ( a%row%extra ) m = m - 1
    n = a%col%nb
    if ( a%col%extra ) n = n - 1
    do i = 1, m
      do j = 1, n
        if ( a%block(i,j)%kind /= m_absent ) &
          & maxAbsVal_1 = max(maxAbsVal_1,maxAbsVal(a%block(i,j)))
      end do
    end do
  end function MaxAbsVal_1

  ! ------------------------------------------------------  MaxL1  -----
  real(r8) function MaxL1 ( A )
  ! Return the L1 norm of the column in A that has the largest L1 norm.
  ! Don't include the extra row or column.
    type(Matrix_T), intent(in) :: A
    integer :: I, J, K, M, N
    real(r8) :: My_L1
    maxL1 = 0.0_r8
    m = a%row%nb
    if ( a%row%extra ) m = m - 1
    n = a%col%nb
    if ( a%col%extra ) n = n - 1
    do j = 1, n
      do k = 1, a%block(1,j)%ncols
        my_L1 = 0.0_r8
        do i = 1, m
          my_L1 = my_L1 + col_L1(a%block(i,j),k)
        end do
        maxL1 = max(maxL1, my_L1)
      end do
    end do
  end function MaxL1

  ! -------------------------------------------  MinDiag_Cholesky  -----
  real(r8) function MinDiag_Cholesky ( A )
  ! Return the magnitude of the element on the diagonal of A that has the
  ! smallest magnitude.  If A has an extra column, ignore it.
    type(Matrix_Cholesky_T), intent(in) :: A
    minDiag_Cholesky = minDiag_1 ( a%m )
  end function MinDiag_Cholesky

  ! ------------------------------------------------  MinDiag_SPD  -----
  real(r8) function MinDiag_SPD ( A )
  ! Return the magnitude of the element on the diagonal of A that has the
  ! smallest magnitude.  If A has an extra column, ignore it.
    type(Matrix_SPD_T), intent(in) :: A
    minDiag_SPD = minDiag_1 ( a%m )
  end function MinDiag_SPD

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
    if ( (x%row%vec%template%id /= y%row%vec%template%id)  .or. &
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
  ! Z = A^T V if UPDATE is absent or false.  Z is first cloned from V.
  ! Z = Z + A^T V is UPDATE is present and true.
  ! The extra row and column are not multiplied, even if present.
    type(Matrix_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: Z
    logical, optional, intent(in) :: UPDATE

    logical :: DO_UPDATE      ! Tells MatrixModule_0 % multiplyMatrixVector
    !                           whether to clear an element of Z, or add to it
    integer :: I, J           ! Subscripts and loop inductors
    integer :: K, L, M, N     ! Subscripts
    integer :: MB, NB         ! Loop limits
    logical :: MY_UPDATE      ! My copy of UPDATE or false if it's absent

    if ( a%row%vec%template%id /= v%template%id ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix and vector not compatible in MultiplyMatrixVector_1" )
    my_update = .false.
    if ( present(update) ) my_update = update
    if ( my_update .and. a%col%vec%template%id /= z%template%id ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix and result not compatible in MultiplyMatrixVector_1" )
    ! Copy characteristics, allocate values:
    if ( .not. my_update ) call cloneVector ( z, v, vectorNameText='_z' )
    mb = a%row%nb
    if ( a%row%extra ) mb = mb - 1
    nb = a%col%nb
    if ( a%col%extra ) nb = nb - 1
    do j = 1, nb
      k = a%col%quant(j)
      l = a%col%inst(j)
      do_update = my_update
      do i = 1, mb
        m = a%row%quant(i)
        n = a%row%inst(i)
        call multiplyMatrixVector ( a%block(i,j), &
          & v%quantities(m)%values(:,n), z%quantities(k)%values(:,l), &
          & do_update )
        do_update = .true.
      end do ! i = 1, a%row%nb
    end do ! j = 1, a%col%nb
  end subroutine MultiplyMatrixVector_1

  ! ----------------------------------  MultiplyMatrixVectorNoT_1  -----
  subroutine MultiplyMatrixVectorNoT_1 ( A, V, Z, UPDATE )
  ! Z = A V if UPDATE is absent or false.  Z is first cloned from the
  !     rows-labeling of A.
  ! Z = Z + A V if UPDATE is presend and true.
    type(Matrix_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: Z
    logical, intent(in), optional :: UPDATE

    logical :: DO_UPDATE      ! Tells MatrixModule_0 % multiplyMatrixVectorNoT
    !                           whether to clear an element of Z, or add to it
    integer :: I, J           ! Subscripts and loop inductors
    integer :: K, L, M, N     ! Subscripts
    logical :: MY_UPDATE      ! My copy of UPDATE or false if it's absent

    if ( a%col%vec%template%id /= v%template%id ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix and vector not compatible in MultiplyMatrixVector_1" )
    if ( a%row%vec%template%id /= z%template%id ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix and result not compatible in MultiplyMatrixVector_1" )
    my_update = .false.
    if ( present(update) ) my_update = update
    ! Copy characteristics, allocate values:
    if ( .not. my_update ) call cloneVector ( z, a%row%vec, vectorNameText='_z' )
    do i = 1, a%row%nb
      m = a%row%quant(i)
      n = a%row%inst(i)
      do_update = my_update
      do j = 1, a%col%nb
        k = a%col%quant(j)
        l = a%col%inst(j)
        call multiplyMatrixVectorNoT ( a%block(i,j), &
          & v%quantities(k)%values(:,l), z%quantities(m)%values(:,n), &
          & do_update )
        do_update = .true.
      end do ! i = 1, a%row%nb
    end do ! j = 1, a%col%nb
  end subroutine MultiplyMatrixVectorNoT_1

  ! ----------------------------------  MultiplyMatrixVectorSPD_1  -----
  subroutine MultiplyMatrixVectorSPD_1 ( A, V, Z, UPDATE )
  ! Z = A V if UPDATE is absent or false.
  ! Z = Z + A V is UPDATE is present and true.
  ! Remember that for SPD, only the upper triangle is stored, so we need
  ! Z = A^T V + A V except don't do the diagonal twice.
  ! The extra row and column are not multiplied, even if present.
    type(Matrix_SPD_T), intent(in) :: A
    type(Vector_T), intent(in) :: V
    type(Vector_T), intent(inout) :: Z
    logical, optional, intent(in) :: UPDATE

    integer :: I, J           ! Subscripts and loop inductors
    integer :: K, L, M, N     ! Subscripts
    integer :: NB             ! Loop bound

    call MultiplyMatrixVector ( a%m, v, z, update ) ! A^T V
    nb = a%m%col%nb
    if ( a%m%col%extra ) nb = nb - 1
    do j = 1, nb
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

  ! ---------------------------------------------------  Negate_1  -----
  subroutine Negate_1 ( A ) ! A = -A
    type(Matrix_T), intent(inout) :: A
    integer :: I, J
    do i = 1, a%row%nb
      do j = 1, a%col%nb
        if ( a%block(i,j)%kind /= m_absent ) &
          & a%block(i,j)%values = -a%block(i,j)%values
      end do
    end do
  end subroutine Negate_1

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
    type(Vector_T), intent(in), optional :: RHS_IN
    type(Vector_T), intent(inout), optional :: RHS_OUT
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
      & call createEmptyMatrix ( z%m, 0, a%col%vec, a%col%vec, &
        &  extra_Row = a%col%extra, extra_Col=a%col%extra )
    do j = 1, a%col%nb
      nullify ( mj )
      if ( .not. ( j == a%col%nb .and. a%col%extra ) ) then
        if ( associated(a%col%vec%quantities(a%col%quant(j))%mask) ) &
          mj => a%col%vec%quantities(a%col%quant(j))%mask(:,a%col%inst(j))
      end if
      do i = 1, j
        nullify ( mi )
        if ( .not. ( i == a%col%nb .and. a%col%extra ) ) then
          if ( associated(a%col%vec%quantities(a%col%quant(i))%mask) ) &
            mi => a%col%vec%quantities(a%col%quant(i))%mask(:,a%row%inst(i))
        end if
        do_update = my_update
        do k = 1, a%row%nb
          call multiplyMatrixBlocks ( &
            & a%block(k,i), a%block(k,j), z%m%block(i,j), update=do_update, &
            & xmask=mi, ymask=mj, upper = i == j )
          do_update = .true.
        end do ! k = 2, a%row%nb
      end do ! i = 1, a%col%nb
    end do ! j = 1, a%col%nb

    if ( present(rhs_in) .and. present(rhs_out) ) &
      & call multiplyMatrixVector ( a, rhs_in, rhs_out, my_update )
  end subroutine NormalEquations

  ! -------------------------------------------------  RowScale_1  -----
  subroutine RowScale_1 ( V, X, NEWX ) ! Z = V X where V is a diagonal
  !                                matrix represented by a vector and Z is X
  !                                or NEWX.
  !                                An "extra" row is not scaled (but an extra
  !                                column is).
    type (Vector_T), intent(in) :: V
    type (Matrix_T), intent(inout), target :: X
    type (Matrix_T), intent(out), target, optional :: NEWX 

    integer :: I, J      ! Subscripts for [XZ]%Block
    integer :: NC, NR    ! Numbers of columns and rows.

    nc = x%col%nb
    nr = x%row%nb
    if ( x%row%extra ) nr = nr - 1
    if ( present(newx) ) then
      call createEmptyMatrix ( newx, 0, x%row%vec, x%col%vec, &
        & x%row%instFirst, x%col%instFirst )
      do j = 1, nc
        do i = 1, nr
          call RowScale ( &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
            & x%block(i,j), newx%block(i,j) )
        end do ! i = nr
      end do ! j = nc
    else
      do j = 1, nc
        do i = 1, nr
          call RowScale ( &
            & v%quantities(x%row%quant(i))%values(:,x%row%inst(i)), &
            & x%block(i,j) )
        end do ! i = nr
      end do ! j = nc
    end if
  end subroutine RowScale_1

  ! ------------------------------------------------  ScaleMatrix  -----
  subroutine ScaleMatrix ( Z, A )       ! Z := A * Z, where A is scalar
    type(matrix_T), intent(inout) :: Z
    real(r8), intent(in) :: A
    integer :: I, J                     ! Subscripts and loop inductors
    do i = 1, z%row%nb
      do j = 1, z%col%nb
        call scaleBlock ( z%block(i,j), a )
      end do ! j
    end do ! i
  end subroutine ScaleMatrix

  ! --------------------------------------------  SolveCholesky_1  -----
  subroutine SolveCholesky_1 ( Z, X, RHS, TRANSPOSE )
  ! Given the Cholesky-factored normal equations Z and the corresponding
  ! RHS, Solve Z X = RHS for X if TRANSPOSE is absent, or present and
  ! false.  Solve Z^T X = RHS for X if TRANSPOSE is present and true.
  ! RHS may be the same as X.  RHS may be absent, in which case X is
  ! assumed to contain the right-hand side on input, and the solution
  ! replaces it on output.
  ! If Z has an extra column, extract it into a vector bu using
  ! GetVectorFromColumn, and pass it in as X.
    type(Matrix_Cholesky_T), intent(in) :: Z
    type(Vector_T), intent(inout), target :: X
    type(Vector_T), intent(in), target, optional :: RHS
    logical, optional, intent(in) :: TRANSPOSE

    integer :: I, J, N             ! Subscripts and loop inductors
    integer :: IC, IR, QC, QR      ! Instance and quantity indices
    type(Vector_T), pointer :: MY_RHS   ! RHS if present, else X
    logical My_transpose           ! TRANSPOSE if present, else .false.

    my_transpose = .false.
    if ( present(transpose) ) my_transpose = transpose
    my_rhs => x
    if ( present(rhs) ) then
      if ( x%template%id /= rhs%template%id ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "X and RHS not compatible in SolveCholesky_1" )
      my_rhs => rhs
    end if
    if ( z%m%col%vec%template%id /= my_rhs%template%id ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Z and RHS not compatible in SolveCholesky_1" )

    n = z%m%col%nb
    if ( z%m%col%extra ) n = n - 1
    if ( my_transpose ) then       ! Solve Z^T X = RHS for X
      do i = 1, n
        ic = z%m%col%inst(i)
        qc = z%m%col%quant(i)
        do j = 1, i-1
          ir = z%m%row%inst(j)
          qr = z%m%row%quant(j)
          ! rhs := rhs - block^T * x
          call multiplyMatrixVector ( z%m%block(j,i), &
            & x%quantities(qr)%values(:,ir), my_rhs%quantities(qc)%values(:,ic), &
            & update=.true., subtract=.true. )
        end do ! j = 1, i-1
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.true. )
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.false. )
      end do ! i = 1, n
    else                           ! Solve Z X = RHS for X
      do i = n, 1, -1
        ic = z%m%row%inst(i)
        qc = z%m%row%quant(i)
        do j = i+1, n
          ir = z%m%col%inst(j)
          qr = z%m%col%quant(j)
          ! rhs := rhs - block * x
          call multiplyMatrixVectorNoT ( z%m%block(i,j), &
            & x%quantities(qr)%values(:,ir), my_rhs%quantities(qc)%values(:,ic), &
            & update=.true., subtract=.true. )
        end do ! j = 1, i-1
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.false. )
        call solveCholesky ( z%m%block(i,i), &
          & my_rhs%quantities(qc)%values(:,ic), transpose=.true. )
      end do ! i = n, 1, -1
    end if
  end subroutine SolveCholesky_1

  ! -------------------------------------------  UpdateDiagonal_1  -----
  subroutine UpdateDiagonal_1 ( A, LAMBDA, SQUARE )
  ! Add LAMBDA to the diagonal of A.  Don't update the extra row or column.
    type(Matrix_SPD_T), intent(inout) :: A
    real(r8), intent(in) :: LAMBDA
    logical, intent(in), optional :: SQUARE ! Update with square of value

    integer :: I, N
    real(r8) :: MYLAMBDA

    myLambda = lambda
    if ( present(square) ) then
      if (square) myLambda = lambda**2
    endif

    n = max(a%m%row%nb,a%m%col%nb)
    if ( a%m%row%extra .or. a%m%col%extra ) n = n - 1
    
    do i = 1, n
      call updateDiagonal ( a%m%block(i,i), myLambda )
    end do

  end subroutine UpdateDiagonal_1

  ! ----------------------------------------  UpdateDiagonalVec_1  -----
  subroutine UpdateDiagonalVec_1 ( A, X, SUBTRACT, SQUARE )
  ! Add X to the diagonal of A.  Don't update the extra row or column.
  ! If SUBTRACT is present and true, subtract X from the diagonal.
    type(matrix_SPD_T), intent(inout) :: A
    type(vector_T), intent(in) :: X
    logical, intent(in), optional :: SUBTRACT
    logical, intent(in), optional :: SQUARE

    integer :: I, N
    logical :: MYSQUARE

    mySquare = .false.
    if ( present(square) ) mySquare = square

    if ( a%m%col%vec%template%id /= x%template%id ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "A and X not compatible in UpdateDiagonalVec_1" )
    n = max(a%m%row%nb,a%m%col%nb)
    if ( a%m%row%extra .or. a%m%col%extra ) n = n - 1
    if ( mySquare ) then
      do i = 1, n
        call updateDiagonal ( a%m%block(i,i), &
          & (x%quantities(a%m%row%quant(i))%values(:,a%m%row%inst(i)))**2, subtract )
      end do
    else
      do i = 1, n
        call updateDiagonal ( a%m%block(i,i), &
          & x%quantities(a%m%row%quant(i))%values(:,a%m%row%inst(i)), subtract )
      end do
    endif
  end subroutine UpdateDiagonalVec_1

! =====     Private Procedures     =====================================

  ! ------------------------------------  AddItemToMatrixDatabase  -----
  integer function AddItemToMatrixDatabase ( Database, Item )

  ! This routine adds a matrix data base item to the database.  These
  ! items are constructed by AddMatrixToDatabase, AddCholeskyToDatabase,
  ! AddKroneckerToDatabase or AddSPDToDatabase.

    type(Matrix_Database_T), dimension(:), pointer :: Database
    type(Matrix_Database_T) :: Item

    type(Matrix_Database_T), dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddItemToMatrixDatabase = newSize
  end function AddItemToMatrixDatabase

  ! -------------------------------------------------  CopyRCInfo  -----
  subroutine CopyRCInfo ( A, B ) ! A := B
    type(RC_info), intent(inout) :: A
    type(RC_info), intent(in) :: B
    call destroyRCInfo ( a )
    a%vec = b%vec
    a%nb = b%nb
    a%extra = b%extra
    a%instFirst = b%instFirst
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

  ! --------------------------------------------------  Dump_Linf  -----
  subroutine Dump_Linf ( Matrix, Name, Upper )
    type(Matrix_T), intent(in) :: Matrix
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional :: Upper
    ! Dump the L_infty norms of the matrix blocks.  Only dump the upper
    ! triangle if Upper is present and true.

    integer :: I, J, K             ! Subscripts, loop inductors
    logical :: My_upper

    if ( present(name) ) call output ( name )
    if ( matrix%name > 0 ) then
      if ( present(name) ) call output ( ', ' )
      call output ( 'Name = ' )
      call display_string ( matrix%name, advance='yes' )
    else
      if ( present(name) ) call output ( '', advance='yes' )
    end if
    my_upper = .false.
    if ( present(upper) ) my_upper = upper
    do k = 1, matrix%col%nb, 7
      if ( matrix%col%nb > 7 ) then
        call output ( ' ' )
        do i = k, min(matrix%col%nb,k+6)
          call output ( i, places=10 )
        end do
        call output ( '', advance='yes' )
      end if
      do i = 1, matrix%row%nb
        if ( .not. my_upper .or. i <= min(matrix%col%nb,k+6) ) then
          call output ( i, places=4 )
          call output ( ':' )
          do j = k, min(matrix%col%nb,k+6)
            if ( my_upper .and. i > j ) then
              call blanks ( 10 )
            else
              call output ( maxAbsVal(matrix%block(i,j)), format='(1pe10.3)' )
            end if
          end do ! j
          call output ( '', advance='yes' )
        end if
      end do ! i
    end do ! k
  end subroutine Dump_Linf

  ! ------------------------------------------------  Dump_Matrix  -----
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
    if ( present(name) ) call output ( name )
    if ( matrix%name > 0 ) then
      if ( present(name) ) call output ( ', ' )
      call output ( 'Name = ' )
      call display_string ( matrix%name, advance='yes' )
    else
      if ( present(name) ) call output ( '', advance='yes' )
    end if
    call dump_rc ( matrix%row, 'row', my_details>0 )
    call dump_rc ( matrix%col, 'column', my_details>0 )
    do j = 1, matrix%col%nb
      do i = 1, matrix%row%nb
        if ( my_details < 1 .and. matrix%block(i,j)%kind == m_absent ) cycle
        call output ( 'Block at row ' )
        call output ( i )
        call output ( ' and column ' )
        call output ( j )
        call output ( ' ( ' )
        if ( matrix%row%extra .and. i == matrix%row%nb ) then
          call output ( '_extra_' )
        else
          call display_string ( &
            & matrix%row%vec%quantities(matrix%row%quant(i))%template%name )
        end if
        call output ( ':' )
        call output ( matrix%row%Inst(i) )
        call output (' , ')
        if ( matrix%col%extra .and. j == matrix%col%nb ) then
          call output ( '_extra_' )
        else
          call display_string ( &
            & matrix%col%vec%quantities(matrix%col%quant(j))%template%name )
        end if
        call output ( ':' )
        call output ( matrix%col%Inst(j) )
        call output ( ' )' )
        if ( matrix%block(i,j)%kind == m_absent ) then
          call output ( ' [absent]', advance='yes' )
        else
          call output ( '', advance='yes' )
          call dump ( matrix%block(i,j), details=my_details )
        end if
      end do
    end do
  end subroutine Dump_Matrix

  ! ---------------------------------------  Dump_Matrix_Database  -----
  subroutine Dump_Matrix_Database ( MatrixDatabase, Details )
    type(Matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    integer, intent(in), optional :: Details   ! Print details, default 1
    !  <= Zero => no details, == One => Details of matrix but not its blocks,
    !  >One => Details of the blocks, too.

    integer :: I

    if ( .not. associated(MatrixDatabase) ) return
    do i = 1, size(MatrixDatabase)
      call output ( i, 4 )
      call output ( ': ' )
      if ( associated(matrixDatabase(i)%matrix) ) then
        call dump ( matrixDatabase(i)%matrix, 'Plain', details )
      else if ( associated(matrixDatabase(i)%cholesky) ) then
        call dump ( matrixDatabase(i)%cholesky%m, 'Cholesky', details )
      else if ( associated(matrixDatabase(i)%kronecker) ) then
        call dump ( matrixDatabase(i)%kronecker%m, 'Kronecker', details )
      else if ( associated(matrixDatabase(i)%spd) ) then
        call dump ( matrixDatabase(i)%spd%m, 'SPD', details )
      end if
    end do
  end subroutine Dump_Matrix_Database

  ! ----------------------------------------------------  Dump_RC  -----
  subroutine Dump_RC ( RC, R_or_C, Details )
    type(rc_info), intent(in) :: RC
    character(len=*), intent(in) :: R_or_C
    logical, intent(in) :: Details
    call output ( 'Number of ' )
    call output ( r_or_c )
    call output ( ' blocks = ' )
    call output ( rc%nb )
    call output ( ' Vector that defines ' )
    call output ( r_or_c )
    call output ( 's' )
    if ( rc%vec%name == 0 ) then
      call output ( ' has no name', advance='yes' )
    else
      call output ( ': ' )
      call display_string ( rc%vec%name, advance='yes' )
    end if
    call output ( 'Order of '//r_or_c//' blocks is ' )
    if ( rc%instFirst ) then
      call output ( 'instance, then quantity', advance='yes' )
    else
      call output ( 'quantity, then instance', advance='yes' )
    end if
    if ( details ) then
      call output ( 'Numbers of ' )
      call output ( r_or_c )
      call output ( 's in each block' )
      if ( rc%extra ) call output ( ' (including extra)' )
      call output ( ':', advance='yes' )
      call dump ( rc%nelts )
      call output ( 'Instance indices for blocks in the ' )
      call output ( r_or_c )
      call output ( 's' )
      if ( rc%extra ) call output ( ' (including extra)' )
      call output ( ':', advance='yes' )
      call dump ( rc%inst )
      call output ( 'Quantity indices for blocks in the ' )
      call output ( r_or_c )
      call output ( 's' )
      if ( rc%extra ) call output ( ' (including extra)' )
      call output ( ':', advance='yes' )
      call dump ( rc%quant )
    end if
  end subroutine Dump_RC

  ! ------------------------------------------------  Dump_Struct  -----
  subroutine Dump_Struct ( Matrix, Name, Upper )
  ! Display the structure of the kinds of the matrix blocks
    type(Matrix_T), intent(in) :: Matrix
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional :: Upper   ! Only do the upper triangle
    !                                          if present and true.

    !                         Absent Banded Sparse   Full
    character :: CHARS(0:3) = (/ '-',   'B',   'S',   'F' /)
    integer :: I, J
    logical :: MyUpper

    if ( present(name) ) call output ( name )
    if ( matrix%name > 0 ) then
      if ( present(name) ) call output ( ', ' )
      call output ( 'Name = ' )
      call display_string ( matrix%name, advance='yes' )
    else
      if ( present(name) ) call output ( '', advance='yes' )
    end if
    myUpper = .false.
    if ( present(upper) ) myUpper = upper
    do i = 1, matrix%row%nb
      call output ( i, 3 )
      call output ( ':' )
      do j = 1, matrix%col%nb
        call output ( ' ' )
        if ( myUpper .and. j < i ) then
          call output ( ' ' )
        else
          call output ( chars(matrix%block(i,j)%kind) )
        end if
      end do ! j
      call output ( '', advance='yes' )
    end do ! i
  end subroutine Dump_Struct

  ! --------------------------------------------------  MinDiag_1  -----
  real(r8) function MinDiag_1 ( A )
  ! Return the magnitude of the element on the diagonal of A that has the
  ! smallest magnitude.  If A has an extra column, ignore it.
    type(Matrix_T), intent(in) :: A
    integer :: I, N
    n = a%col%nb
    if ( a%col%extra ) n = n - 1
    minDiag_1 = minDiag(a%block(1,1))
    do i = 2, n
      minDiag_1 = min(minDiag_1,minDiag(a%block(i,i)))
    end do
  end function MinDiag_1
end module MatrixModule_1

! $Log$
! Revision 2.42  2001/05/22 19:09:33  vsnyder
! Implement MaxL1
!
! Revision 2.41  2001/05/19 00:20:05  livesey
! Cosmetic changes from Van.
!
! Revision 2.40  2001/05/19 00:14:57  livesey
! OK that should have been square (idiot!)
!
! Revision 2.39  2001/05/19 00:13:23  livesey
! Added squareRoot option to update diagonal
!
! Revision 2.38  2001/05/18 23:48:34  vsnyder
! Correct Dump_L1 -> Dump_Linf
!
! Revision 2.37  2001/05/18 22:28:11  vsnyder
! Don't look for a mask in the extra column during NormalEquations
!
! Revision 2.36  2001/05/17 20:19:20  vsnyder
! Implement GetMatrixElement.  Change handling of mask in NormalEquations.
! Don't scale the extra column/row if column/row scaling, but do scale
! the extra row/column.
!
! Revision 2.35  2001/05/12 18:58:47  livesey
! Fixed a bug, not sure it's what Van intended but it should compile
!
! Revision 2.34  2001/05/12 01:07:19  vsnyder
! Some repairs in RowScale and ColumnScale
!
! Revision 2.33  2001/05/10 22:54:34  vsnyder
! Get CholeskyFactor_1 to work.  Add Dump_L1 and Dump_Struct.
!
! Revision 2.32  2001/05/10 02:14:58  vsnyder
! Repair CholeskyFactor, MaxAbsVal; add Dump_Struct
!
! Revision 2.31  2001/05/09 19:46:06  vsnyder
! Add BandHeight argument to CreateBlock_1
!
! Revision 2.30  2001/05/09 01:56:15  vsnyder
! periodic commit -- Work on block Cholesky
!
! Revision 2.29  2001/05/08 20:29:40  vsnyder
! Periodic commit -- workong on sparse matrix blunders
!
! Revision 2.28  2001/05/03 02:11:23  vsnyder
! Spiffify dump, add names to cloned vectors
!
! Revision 2.27  2001/05/01 23:54:13  vsnyder
! Create a block for the extra column
!
! Revision 2.26  2001/05/01 06:56:32  livesey
! Bug fix, was using optional argument, not local copy/default.
!
! Revision 2.25  2001/04/30 23:44:25  vsnyder
! Correct/remove some incorrect size tests in MultiplyMatrixVectorNoT
!
! Revision 2.24  2001/04/28 07:03:59  livesey
! Removed print statement
!
! Revision 2.23  2001/04/28 04:42:29  livesey
! Removing some of the unnecessary(?) assertions of square matrices in
! multiplyMatrixVector and its relatives.  Also sorted out some of the
! conditions, and loops for the multiplyMatrixVectorNoT case.
!
! Revision 2.22  2001/04/28 01:28:36  livesey
! Change in rcInfo_T, vec is now a copy of the parent vector, not a pointer.
!
! Revision 2.21  2001/04/27 22:51:52  livesey
! Some changes/improvements to dump
!
! Revision 2.20  2001/04/26 23:56:02  livesey
! Fix to test for MatrixVectorMultiplyNoT
!
! Revision 2.19  2001/04/25 01:12:39  vsnyder
! Improve some comments
!
! Revision 2.18  2001/04/25 00:50:25  vsnyder
! Provide MultiplyMatrixNoT
!
! Revision 2.17  2001/04/24 22:35:56  vsnyder
! Maybe this time elements of matrixDatabase are destroyed coimpletely
!
! Revision 2.16  2001/04/21 02:11:02  vsnyder
! Fix a memory leak
!
! Revision 2.15  2001/04/20 02:56:18  livesey
! Added createblock_1 as public and overloaded.  Also made dump more informative
!
! Revision 2.14  2001/04/11 00:40:25  vsnyder
! Remove some inadventently-left-in debugging print
!
! Revision 2.13  2001/04/11 00:03:42  vsnyder
! Repair matrix creation, improve 'dump'
!
! Revision 2.12  2001/04/10 00:19:11  vsnyder
! Add GetKindFromMatrixDatabase and necessary parameters
!
! Revision 2.11  2001/04/09 23:56:17  vsnyder
! Change some pointer arguments to targets
!
! Revision 2.10  2001/04/09 23:32:19  vsnyder
! Correct typo
!
! Revision 2.9  2001/02/22 02:09:36  vsnyder
! OOPS -- Forgot to make InvertCholesky public
!
! Revision 2.8  2001/02/22 01:55:06  vsnyder
! Add code to invert a Cholesky factor
!
! Revision 2.7  2001/01/26 19:00:02  vsnyder
! Periodic commit
!
! Revision 2.6  2001/01/19 23:53:26  vsnyder
! Add FillExtraCol, GetVectorFromColumn (incomplete); SolveCholesky_1
! still needs work, too.
!
! Revision 2.5  2001/01/10 21:03:14  vsnyder
! Periodic commit
!
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
