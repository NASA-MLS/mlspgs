! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MatrixModule_0          ! Low-level Matrices in the MLS PGS suite
!=============================================================================

! This module provides the elementary matrix type.  Blocks of this
! type are used to compose block matrices.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT

  implicit NONE
  private
  public :: Add_Matrix_Blocks, Assignment(=), CholeskyFactor
  public :: CholeskyFactor_0, ClearRows, ClearRows_0, CloneBlock, ColumnScale
  public :: ColumnScale_0, Col_L1, CopyBlock, CreateBlock, CreateBlock_0
  public :: DenseCholesky, Densify, DestroyBlock, DestroyBlock_0, Dump
  public :: GetDiagonal, GetDiagonal_0
  public :: GetMatrixElement, GetMatrixElement_0, GetVectorFromColumn
  public :: GetVectorFromColumn_0, M_Absent, M_Banded, M_Column_Sparse, M_Full
  public :: MatrixElement_T, MaxAbsVal, MaxAbsVal_0, MinDiag, MinDiag_0
  public :: Multiply, MultiplyMatrixBlocks, MultiplyMatrixVector
  public :: MultiplyMatrixVector_0, MultiplyMatrixVectorNoT
  public :: MultiplyMatrixVectorNoT_0, operator(+), operator(.TX.), RowScale
  public :: RowScale_0, ScaleBlock, SolveCholesky, SolveCholeskyM_0
  public :: SolveCholeskyV_0, Sparsify, UpdateDiagonal, UpdateDiagonal_0
  public :: UpdateDiagonalVec_0

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment(=)
    module procedure AssignBlock
  end interface

  interface CholeskyFactor
    module procedure CholeskyFactor_0
  end interface

  interface ClearRows
    module procedure ClearRows_0
  end interface

  interface ColumnScale
    module procedure ColumnScale_0
  end interface

  interface CreateBlock
    module procedure CreateBlock_0
  end interface

  interface DestroyBlock
    module procedure DestroyBlock_0
  end interface

  interface DUMP
    module procedure DUMP_MATRIX_BLOCK
  end interface

  interface GetDiagonal
    module procedure GetDiagonal_0
  end interface

  interface GetMatrixElement
    module procedure GetMatrixElement_0
  end interface

  interface GetVectorFromColumn
    module procedure GetVectorFromColumn_0
  end interface

  interface MaxAbsVal
    module procedure MaxAbsVal_0
  end interface

  interface MinDiag
    module procedure MinDiag_0
  end interface

  interface Multiply
    module procedure MultiplyMatrixBlocks, MultiplyMatrixVector_0
  end interface

  interface MultiplyMatrixVector ! A^T V
    module procedure MultiplyMatrixVector_0
  end interface

  interface MultiplyMatrixVectorNoT ! A V
    module procedure MultiplyMatrixVectorNoT_0
  end interface

  interface operator (+)
    module procedure Add_Matrix_Blocks
  end interface

  interface operator ( .TX. ) ! A^T * B
    module procedure NewMultiplyMatrixBlocks, NewMultiplyMatrixVector_0
  end interface

  interface RowScale
    module procedure RowScale_0
  end interface

  interface SolveCholesky
    module procedure SolveCholeskyM_0, SolveCholeskyV_0
  end interface

  interface UpdateDiagonal
    module procedure UpdateDiagonal_0, UpdateDiagonalVec_0
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Parameters for the KIND component of objects of type(MatrixElement_T):
  integer, parameter :: M_Absent = 0         ! An absent block -- assumed zero
  integer, parameter :: M_Banded = 1         ! A banded matrix.  The nonzero
    ! values in each column are assumed to be in a contiguous sequence.  For
    ! column I, R1(I) gives the row of the first nonzero value, and the
    ! nonzero elements are VALUES(R2(I-1)+1:R2(I)).
  integer, parameter :: M_Column_sparse = 2  ! A sparse block in column-sparse
    ! representation.  Only the non-zero values are stored.  For column(I),
    ! R1(I) gives the index in R2 and VALUES for the last stored value.  The
    ! rows of the nonzero values are R2(R1(I-1)+1:R1(I)), and their values are
    ! VALUES(R1(I-1)+1:R1(I)).  R1(0) = 0.
  integer, parameter :: M_Full = 3           ! A non-sparse block

  type MatrixElement_T
    integer :: KIND = M_Absent               ! Kind of block -- one of the
      !                                        M_... parameters above
    integer :: NROWS = 0, NCOLS = 0          ! Numbers of rows and columns
    integer, pointer, dimension(:) :: R1 => NULL()     ! Indexed by the column
      ! number. Used for the first column number if KIND = M_Banded, as
      ! described above for M_Column_sparse if KIND = M_Column_sparse, and not
      ! used otherwise.
    integer, pointer, dimension(:) :: R2 => NULL()     ! Indexed by the
      ! column number if KIND = M_Banded, by elements of R1 if KIND =
      ! M_Column_sparse, and not used otherwise.  See M_Banded and
      ! M_Column_sparse above.
    real(r8), pointer, dimension(:,:) :: VALUES => NULL()   ! Values of the
      ! matrix elements.  Indexed by row and column indices if KIND == M_Full,
      ! by elements in the range of values of R1 if KIND == M_Banded, and by
      ! elements of R2 if KIND == M_Column_sparse.
  end type MatrixElement_T

  ! - - -  Private data     - - - - - - - - - - - - - - - - - - - - - -
  integer, parameter, private :: B_sizer = 0
  integer, parameter, private :: B = bit_size(b_sizer) ! can't use "bit_size(b)"
  real, parameter, private :: COL_SPARSITY = 0.5  ! If more than this
    ! fraction of the elements between the first and last nonzero in a
    ! column are nonzero, use M_Banded, otherwise use M_Column_Sparse.
  real, parameter, private :: SPARSITY = 0.33     ! If a full matrix has
    ! a greater fraction of nonzeroes than specified by this number, there's
    ! no point in making it sparse.

  ! It is important to lie about the interface for the Y argument for *DOT.
  ! We pass an element of a rank-2 object to it, and assume the subsequent
  ! elements are consecutive in memory with stride INCY.  *DOT then treats
  ! this as a rank-1 stride INCY vector.
  interface DOT
    real function SDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      real, intent(in) :: X, Y
    end function SDOT
    double precision function DDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      double precision, intent(in) :: X, Y
    end function DDOT
  end interface

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------  Add_Matrix_Blocks  -----
  function Add_Matrix_Blocks ( XB, YB ) result ( ZB ) ! ZB = XB + YB
    type(MatrixElement_T), intent(in), target :: XB, YB
    type(MatrixElement_T) :: ZB

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K, L, N
    type(MatrixElement_T), pointer :: X, Y
    real(kind(zb%values)), pointer :: Z(:,:) ! May be used if Y is col-sparse
    if ( xb%kind <= yb%kind ) then
      x => xb
      y => yb
    else
      x => yb
      y => xb
    end if
    ! The structure of cases-within-cases below depends on the order of
    ! the M_... parameters, because order of the KIND field values is
    ! used to determine whether to commute the operands.  The M_...
    ! parameters are declared in alphabetical order, and their values
    ! are in the same order as their declarations.  The kind of the XB
    ! operand is less than or equal to the kind of the YB operand.
    if ( x%kind == M_Absent ) then
      call copyBlock ( zb, y )                   ! Zb = y
      return
    end if
    if ( xb%nrows /= yb%nrows .or. xb%ncols /= yb%ncols ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix sizes incompatible in Add_Matrix_Blocks" )
    select case ( x%kind )
    case ( M_Banded )
      select case ( y%kind )
      case ( M_Banded )                          ! X banded, Y banded
        call allocate_test ( zb%r1, size(x%r1), "zb%r1", ModuleName )
        zb%r1 = min(x%r1, y%r1)                  ! First nonzero row
        call allocate_test ( zb%r2, size(x%r2), "zb%r2", ModuleName, lowBound=0 )
        zb%r2(0) = 0
        do k = 1, size(x%r1)                     ! Calculate size of Values
          zb%r2(k) = zb%r2(k-1) + &
            & max( x%r2(k)-x%r2(k-1)+x%r1(k), y%r2(k)-y%r2(k-1)+y%r1(k) ) - &
            & zb%r1(k) + 1
        end do
        call allocate_test ( zb%values, zb%r2(size(x%r1)), 1, "zb%values", &
          & ModuleName )
        zb%values = 0.0_r8 ! ??? Improve this by only filling
        !                    ??? values that don't get set below
        do k = 1, size(x%r1)
          i = 1; j = 1; l = 1
          n = y%r1(k) - x%r1(k)
          if ( n > 0 ) then ! Copy rows in X that are not in Y
            zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
              & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1)
            i = i + n
          else if ( n < 0 ) then ! Copy rows in Y that are not in X
            n = - n
            zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
              & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
            j = j + n
          end if
          l = l + n
          n = min( x%r1(k) + x%r2(k) - x%r2(k-1) - 1, &
            &      y%r1(k) + y%r2(k) - y%r2(k-1) - 1 ) - &
            &      max( x%r1(k), y%r1(k) )
          zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
            & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1) + &
            & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
          i = i + n; j = j + n; l = l + n
          n = y%r1(k) + y%r2(k) - y%r2(k-1) - &
            & ( x%r1(k) + x%r2(k) - x%r2(k-1) )
          if ( n > 0 ) then ! Copy rows of Y that are not in X
            zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
              & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
          else if ( n < 0 ) then ! Copy rows in X that are not in Y
            n = - n
            zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
              & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1)
          end if
        end do
      case ( M_Column_sparse )                   ! X banded, Y col sparse
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, y%nrows, y%ncols, "Z in Add_Matrix_Blocks", &
          & ModuleName )
        call densify ( z, y )                    ! z = y%values
        do k = 1, size(x%r1)
          z(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) = &
            & z(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) + &
              & x%values(x%r2(k-1)+1:x%r2(k), 1)
        end do
        call sparsify ( z, zb, "Z in Add_Matrix_Blocks", ModuleName ) ! Zb = Z
      case ( M_Full )                            ! X banded, Y full
        call CopyBlock ( zb, y )                 ! Zb = y
        do k = 1, size(x%r1)
          zb%values(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) = &
            & zb%values(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) + &
              & x%values(x%r2(k-1)+1:x%r2(k), 1)
        end do
      end select
    case ( M_Column_sparse )
      select case ( y%kind )
!     case ( M_Banded )        ! Not needed because of commuted arguments
      case ( M_Column_sparse )                   ! X col sparse, Y col sparse
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, y%nrows, y%ncols, "Z in Add_Matrix_Blocks", &
          & ModuleName )
        call densify ( z, y )                    ! z = y%values
        do k = 1, size(x%r1)
          z(x%r2(x%r1(k-1)+1:x%r1(k)), k) = &
            & z(x%r2(x%r1(k-1)+1:x%r1(k)), k) + &
              & x%values(x%r1(k-1)+1:x%r1(k),1)
        end do
        call sparsify ( z, zb, "Z in Add_Matrix_Blocks", ModuleName ) ! Zb = Z
      case ( M_Full )                            ! X col sparse, Y full
        call CopyBlock ( zb, y )                 ! Zb = y

! Commented-out on account of internal NAG v4.0 bug
         do k = 1, size(x%r1)
           zb%values(x%r2(x%r1(k-1)+1:x%r1(k)), k) = &
             & zb%values(x%r2(x%r1(k-1)+1:x%r1(k)), k) + &
               & x%values(x%r1(k-1)+1:x%r1(k),1)
         end do
      end select
    case ( M_Full )
      select case ( y%kind )
!     case ( M_Banded )        ! Not needed because of commuted arguments
!     case ( M_Column_sparse ) ! Not needed because of commuted arguments
      case ( M_Full )                            ! X full, Y full
        call CloneBlock ( zb, y )                ! Zb = y, except the values
        zb%values = x%values + y%values
      end select
    end select
  end function Add_Matrix_Blocks

  ! ------------------------------------------------  AssignBlock  -----
  subroutine AssignBlock ( Z, X )
  ! Destroy Z, then copy X to it, using pointer assignment for pointer
  ! components.  Other than the "Destroy Z" part, the semantics are the
  ! same as for intrinsic assignment.  Notice that CopyBlock does a deep
  ! copy.  If one has Z = X in a loop, it is therefore necessary only to
  ! destroy Z after the loop.
    type(MatrixElement_T), intent(inout) :: Z
    type(MatrixElement_T), intent(in) :: X
    call destroyBlock ( z )
    z%kind = x%kind
    z%nRows = x%nRows; z%nCols = x%nCols
    z%r1 => x%r1
    z%r2 => x%r2
    z%values => x%values
  end subroutine AssignBlock

  ! -------------------------------------------  CholeskyFactor_0  -----
  subroutine CholeskyFactor_0 ( Z, XOPT )
  ! If XOPT is present compute Z such that Z^T Z = XOPT and Z is upper-
  ! triangular. Otherwise, replace Z such that Z(output)^T Z(output) =
  ! Z(input) and Z(output) is upper-triangular.
    type(MatrixElement_T), target, intent(inout) :: Z
    type(MatrixElement_T), target, intent(in), optional :: XOPT

    real(r8) :: D             ! Diagonal(I,I) element of Z
    real(r8) :: G             ! X(I,J) - dot_product(X(1:i-1,i),X(1:i-1,j))
    integer :: I, J           ! Subscripts and loop inductors
    integer :: II, IJ         ! Subscripts in VALUES for I,I and I,J components
    !                           in the case of M_Banded
    integer :: NC             ! Number of columns
    integer, pointer, dimension(:) :: R1      ! First nonzero row of Z (Banded)
    integer :: RZ             ! Starting row in Z for inner product (Banded)
    real(r8), save :: TOL = -1.0_r8
    real(r8), pointer, dimension(:,:) :: XIN  ! A pointer to the input,
    !                           data, or a densified copy of it
    type(MatrixElement_T), pointer :: X       ! XOPT or Z, depending on whether
    !                           XOPT is present or absent, respectively.
    real(r8), pointer, dimension(:,:) :: ZT   ! A local full result that is
    !                           sparsified at the end.

    if ( tol < 0.0_r8 ) tol = sqrt(tiny(0.0_r8))
    nullify ( r1, xin, zt )
    x => z
    if ( present(xopt) ) x => xopt
    nc = x%ncols
    if ( nc /= x%nrows )&
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Cannot CholeskyFactor a non-square block" )
    select case ( x%kind )
    case ( M_Absent )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot CholeskyFactor an empty block" )
    case ( M_Banded )
      call allocate_test ( zt, nc, nc, "ZT in CholeskyFactor", &
        & ModuleName )
      call allocate_test ( r1, nc, "R1 in CholeskyFactor", ModuleName )
      do i = 1, nc
        r1(i) = i             ! We know the diagonal will get a value
      end do
      do i = 1, nc
        ii = i - x%r1(i)      ! Offset in VALUES of (I,I) element
        if ( ii < 0 .or. ii > x%r2(i) - x%r2(i-1) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        end if
        zt(i,1:i-1) = 0.0_r8  ! Clear left from the diagonal (helps Sparsify!)
!       g = x%values(ii+x%r2(i-1)+1,1) - &
!           & dot_product(zt(r1(i):i-1,i),zt(r1(i):i-1,i))
        g = x%values(ii+x%r2(i-1)+1,1) - &
            & dot( i-r1(i), zt(r1(i),i), 1, zt(r1(i),i), 1 )
        if ( g <= tol .and. i < nc ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        end if
        d = sqrt(g)
        zt(i,i) = d
        do j = i+1, nc
          ij = i - x%r1(j)    ! Offset in VALUES of (I,J) element
          rz = max(r1(i),r1(j))
!         g = - dot_product(zt(rz:i-1,i),zt(rz:i-1,j))
          g = - dot( i-rz, zt(rz,i), 1, zt(rz,j), 1 )
          if ( ij >= 0 .and. ij <= x%r2(j) - x%r2(j-1) ) &
            & g = x%values(ij+x%r2(j-1)+1,1) + g
          zt(i,j) = g / d
          if ( abs(zt(i,j)) >= tiny(0.0_r8) ) r1(j) = min( r1(j), i )
        end do ! j
      end do ! i
      ! Sparsify the result.  We do it here because Sparsify is slow, and
      ! we know something about the structure: the first and last nonzero
      ! rows in each column.
      j = 1 ! number of nonzeroes
      do i = 2, nc
        j = j + i + 1 - r1(i)
      end do
      call createBlock ( z, nc, nc, M_Banded, j )
      z%r1 = r1
      do i = 1, nc
        z%r2(i) = z%r2(i-1) + i + 1 - r1(i)
        z%values(z%r2(i-1)+1:z%r2(i),1) = zt(r1(i):i,i)
      end do ! i
      call deallocate_test ( r1, "R1 in CholeskyFactor", ModuleName )
      call deallocate_test ( zt, "ZT in CholeskyFactor", ModuleName )
    case ( M_Column_Sparse )
      call allocate_test ( zt, nc, nc, "ZT for CholeskyFactor", ModuleName )
      call allocate_test ( xin, nc, nc, "XIN for CholeskyFactor", ModuleName )
      ! ??? Densify and then compute Cholesky decomposition of dense block.
      ! ??? If necessary, improve this in level 1.0 by working directly
      ! ??? with sparse input.
      call densify ( xin, x )
      call denseCholesky ( zt, xin )
      call sparsify ( zt, z, "ZT in CholeskyFactor", ModuleName ) ! Z := Zt
      call deallocate_test ( xin, "XIN in CholeskyFactor", ModuleName )
    case ( M_Full )
      if ( .not. associated(x,z) ) then
        call createBlock ( z, nc, nc, M_Full )
      end if
      call denseCholesky ( z%values, x%values )
    end select
  end subroutine CholeskyFactor_0

  ! ------------------------------------------------  ClearRows_0  -----
  subroutine ClearRows_0 ( X, MASK )
  ! Clear the rows of X for which MASK has nonzero bits.
    type(MatrixElement_T), intent(inout) :: X
    integer, dimension(1:), intent(in) :: MASK
    integer :: I, J                ! Subscripts and row indices
    integer :: NB, NW              ! Indices of bit and word in MASK
    select case ( x%kind )
    case ( M_Absent )
    case ( M_Banded )              ! ??? Adjust the sparsity representation ???
      do j = 1, x%ncols
        do i = x%r1(j), x%r1(j) + x%r2(j) - x%r2(j-1) - 1 ! row numbers
          if ( btest( mask(i/b+1), mod(i,b) ) ) &
            & x%values(x%r2(j-1) + i - x%r1(j) + 1, 1) = 0.0_r8
        end do ! i
      end do ! j = 1, x%ncols
    case ( M_Column_Sparse )       ! ??? Adjust the sparsity representation ???
      do j = 1, x%ncols
        do i = x%r1(j-1)+1, x%r1(j)
          if ( btest( mask(x%r2(i)/b+1), mod(x%r2(i),b) ) ) &
            & x%values(j,1) = 0.0_r8
        end do ! i
      end do ! j = 1, x%ncols
    case ( M_Full )
      i = 0
      do nw = lbound(mask,1), ubound(mask,1)
        i = i + 1
        if ( i > x%nrows ) return
        do nb = 0, b-1
          if ( btest( mask(nw), nb ) ) then
              x%values(i,:) = 0.0_r8
          end if
        end do ! nb = 0, b-1
      end do ! nw = lbound(mask,1), ubound(mask,1)
    end select
  end subroutine ClearRows_0

  ! -------------------------------------------------  CloneBlock  -----
  subroutine CloneBlock ( Z, X ) ! Z = X, except the values
  ! Duplicate a matrix block, including copying all of its structural
  ! descriptive information, but not its values.
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    call destroyBlock ( z )
    if ( x%kind == M_absent ) then
      call CreateEmptyBlock ( z )
    else
      z%kind = x%kind
      call allocate_test ( z%r1, ubound(x%r1,1), "z%r1", ModuleName, &
        & lowBound=lbound(x%r1,1) )
      z%r1 = x%r1
      call allocate_test ( z%r2, ubound(x%r2,1), "z%r2", ModuleName, &
        lowBound=lbound(x%r2,1) )
      z%r2 = x%r2
      call allocate_test ( z%values, size(x%values,1), size(x%values,2), &
        & "z%values", ModuleName )
    end if
    z%nRows = x%nRows; z%nCols = x%nCols
  end subroutine CloneBlock

  ! -------------------------------------------------  ColumnScale_0  -----
  subroutine ColumnScale_0 ( X, V, NEWX ) ! Z = X V where V is a diagonal
  !                                         matrix represented by a vector
  !                                         and Z is either X or NEWX.
    type(MatrixElement_T), intent(inout), target :: X
    real(r8), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(out), target, optional :: NEWX
    type(MatrixElement_T), pointer :: Z

    integer :: I

    if ( present(newx) ) then
      z => newx
      call copyBlock ( z, x )
    else
      z => x
    end if
    select case ( z%kind )
    case ( M_Absent )
    case ( M_Banded )
      do i = 1, z%ncols
        z%values(z%r2(i-1)+1:z%r2(i),1) = x%values(z%r2(i-1)+1:z%r2(i),1) * v(i)
      end do
    case ( M_Column_Sparse )
      do i = 1, z%ncols
        z%values(z%r1(i-1)+1:z%r1(i),1) = x%values(z%r1(i-1)+1:z%r1(i),1) * v(i)
      end do
    case ( M_Full )
      do i = 1, z%ncols
        z%values(:,i) = x%values(:,i) * v(i)
      end do
    end select
  end subroutine ColumnScale_0

  ! -----------------------------------------------------  Col_L1  -----
  real(r8) function Col_L1 ( X, N )
  ! Return the L1 norm of column N of X
    type(MatrixElement_T), intent(in) :: X
    integer, intent(in) :: N
    integer :: I
    col_l1 = 0.0_r8
    select case ( x%kind )
    case ( M_Absent )
    case ( M_Banded )
      do i = x%r2(n-1)+1, x%r2(n)
        col_l1 = col_l1 + abs(x%values(i,1))
      end do
    case ( M_Column_Sparse )
      do i = x%r1(n-1)+1, x%r1(n)
        col_l1 = col_l1 + abs(x%values(i,1))
      end do
    case ( M_Full )
      do i = 1, x%nrows
        col_l1 = col_l1 + abs(x%values(i,n))
      end do
    end select
  end function Col_L1

  ! --------------------------------------------------  CopyBlock  -----
  subroutine CopyBlock ( Z, X ) ! Destroy Z, deep Z = X, including the values
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    call CloneBlock ( Z, X )
    if ( x%kind /= m_absent ) z%values = x%values
  end subroutine CopyBlock

  ! ----------------------------------------------  CreateBlock_0  -----
  subroutine CreateBlock_0 ( Z, nRows, nCols, Kind, NumberNonzero, NoValues, &
    & BandHeight )
  ! Create a matrix block, but don't fill any elements or structural
  ! information.  The "NumberNonzero" is required if and only if the
  ! "Kind" argument has the value M_Banded or M_Column_Sparse.
  ! The block is first destroyed, so as not to have a memory leak.
  ! If NoValues is present and true, the values component is not allocated
  ! If Kind == M_Banded and BandHeight is present, the band height is
  !   assumed to be uniform, and the R1 and R2 components are filled to
  !   reflect that assumption.

  ! Filling the block after it's created depends on the kind.
  !  M_Absent: Do nothing
  !  M_Banded: The arrays R1 and R2 are indexed by the column number (c).
  !   R1(c) gives the index of the first nonzero row.  R2(c) gives the
  !   subscript in the first dimension of VALUES for the last nonzero element
  !   in the column.  The second dimension of VALUES has shape (1:1).  The
  !   subscript for the first nonzero element is R2(c-1)+1 (R2(0)==0 is set
  !   here).  The number of nonzero elements is R2(c) - R2(c-1).  The index
  !   of the last nonzero row is R1(c) + R2(c) - R2(c-1) - 1.
  !  M_Column_Sparse:  The array R1 is indexed by the column number (c). It
  !   gives the subscripts in R2 and the first dimension of VALUES for
  !   the last entry in the column.  The first one is in R2(R1(c-1)+1)
  !   (R1(0)==0). The number of nonzero entries in a column is R1(c) -
  !   R1(c-1).  The row number of the k'th nonzero in column c is
  !   R2(R1(c-1)+k), and its value is VALUES(R1(c-1)+k),1).
  !  M_Full: R1 and R2 are not used (they are allocated with zero extent).
  !   The value of the (i,j) element of the block is VALUES(i,j).

    type(MatrixElement_T), intent(inout) :: Z
    integer, intent(in) :: nRows, nCols, Kind
    integer, intent(in), optional :: NumberNonzero ! Only for M_Banded and
                                                   ! M_Column_Sparse
    logical, intent(in), optional :: NoValues
    integer, intent(in), optional :: BandHeight

    integer :: I
    logical :: Values

    values = .true.
    if ( present(noValues) ) values = .not. noValues
    call destroyBlock ( z )
    select case ( kind )
    case ( M_Absent )
      call CreateEmptyBlock ( z )
    case ( M_Banded )
      call allocate_test ( z%r1, nCols, "z%r1", ModuleName )
      call allocate_test ( z%r2, nCols, "z%r2", ModuleName, lowBound=0 )
      z%r2(0) = 0
      if ( present(bandHeight) ) then
        do i = 1, nCols
          z%r2(i) = i * bandHeight
          z%r1(i) = 1 + z%r2(i) - bandHeight
        end do
      end if
      if ( values ) &
        & call allocate_test ( z%values, numberNonzero, 1, "z%values", ModuleName )
    case ( M_Column_sparse )
      call allocate_test ( z%r1, nCols, "z%r1", ModuleName, lowBound=0 )
      z%r1(0) = 0
      call allocate_test ( z%r2, numberNonzero, "z%r2", ModuleName )
      if ( values ) &
        & call allocate_test ( z%values, NumberNonzero, 1, "z%values", ModuleName )
    case ( M_Full )
      call allocate_test ( z%r1, 0, "z%r1", ModuleName )
      call allocate_test ( z%r2, 0, "z%r2", ModuleName )
      if ( values ) &
        & call allocate_test ( z%values, nRows, nCols, "z%values", ModuleName )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid matrix block kind in CreateBlock" )
    end select
    z%nRows = nRows
    z%nCols = nCols
    z%kind = kind
  end subroutine CreateBlock_0

  ! ----------------------------------------------  DenseCholesky  -----
  subroutine DenseCholesky ( ZT, XIN, Status )
  ! Do the Cholesky decomposition of XIN giving ZT.
    real(r8), intent(inout) :: ZT(:,:) ! Inout in case it's associated with XIN
    real(r8), intent(inout) :: XIN(:,:) ! Inout in case it's associated with Z
    integer, intent(out), optional :: Status

    real(r8) :: D
    real(r8), save :: TOL = -1.0_r8
    integer :: I, J, NC

    nc = size(xin,2)
    if ( tol < 0.0_r8 ) tol = sqrt(tiny(0.0_r8))
    do i = 1, nc
      zt(i+1:nc,i) = 0.0_r8 ! Clear below the diagonal (helps Sparsify!)
!     d = xin(i,i) - dot_product(zt(1:i-1,i),zt(1:i-1,i))
      d = xin(i,i) - dot( i-1, zt(1,i), 1, zt(1,i), 1 )
      if ( d <= tol .and. i < nc ) then
        if ( present(status ) ) then
          status = i
          return
        end if
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Matrix in DenseCholesky is not positive-definite." )
      end if
      d = sqrt(d)
      zt(i,i) = d
      do j = i+1, nc
!       zt(i,j) = ( xin(i,j) - dot_product(zt(1:i-1,i),zt(1:i-1,j)) ) / d
        zt(i,j) = ( xin(i,j) - dot( i-1, zt(1,i), 1, zt(1,j), 1 ) ) / d
      end do ! j
    end do ! i
    if ( present(status) ) status = 0
  end subroutine DenseCholesky

  ! ----------------------------------------------------  Densify  -----
  subroutine Densify ( Z, B )
  ! Given a matrix block B, produce a full matrix Z, even if the matrix
  ! block had a sparse representation.
    real(r8), intent(out) :: Z(:,:)          ! Full matrix to produce
    type(MatrixElement_T), intent(in) :: B   ! Input matrix block
    integer :: I                             ! Column index

    if ( size(z,1) /= b%nRows .or. size(z,2) /= b%nCols ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Incompatible shapes in Densify" )
    select case ( b%kind )
    case ( M_Absent )
      z = 0.0_r8
    case ( M_Banded )
      z = 0.0_r8
      do i = 1, b%nCols
        z(b%r1(i):b%r1(i)+b%r2(i)-b%r2(i-1)-1,i) = &
          & b%values(b%r2(i-1)+1:b%r2(i),1)
      end do ! i
    case ( M_Column_Sparse )
      z = 0.0_r8
      do i = 1, b%nCols
        z(b%r2(b%r1(i-1)+1:b%r1(i)),i) = b%values(b%r1(i-1)+1:b%r1(i),1)
      end do ! i
    case ( M_Full )
      z = b%values
    end select
  end subroutine DENSIFY

  ! ---------------------------------------------  DestroyBlock_0  -----
  subroutine DestroyBlock_0 ( B )
    ! Deallocate the pointer components of the matrix block B.  Change
    ! its kind to M_Absent.  Don't clobber b%nrows and b%ncols.
    type(MatrixElement_T), intent(inout) :: B
    ! Don't bother to destroy absent blocks
    if (b%kind /= M_Absent) then
      call deallocate_test ( b%r1, "b%r1", ModuleName )
      call deallocate_test ( b%r2, "b%r2", ModuleName )
      call deallocate_test ( b%values, "b%values", ModuleName )
      b%kind = m_absent
    endif
  end subroutine DestroyBlock_0

  ! ----------------------------------------------  GetDiagonal_0  -----
  subroutine GetDiagonal_0 ( B, X, SquareRoot )
  ! Get the diagonal elements of B into X.  Return the square root of the
  ! diagonal elements if SquareRoot is present and true.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(out) :: X
    logical, intent(in), optional :: SquareRoot

    integer :: I, J, N

    n = min(b%nrows,b%ncols)
    if ( n > size(x) ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Array "X" to small in GetDiagonal_0' )
    select case ( b%kind )
    case ( M_Absent )
      x = 0.0_r8
    case ( M_Banded )
      do i = 1, n
        if ( b%r1(i) <= i .and. b%r1(i) + b%r2(i) - b%r2(i-1) > i ) then
          x(i) = b%values(b%r2(i-1)+i-b%r1(i)+1,1)
        else
          x(i) = 0.0_r8
        end if
      end do
    case ( M_Column_Sparse )
      x = 0.0_r8
      do j = b%r1(i-1)+1, b%r1(i)
        if ( b%r2(j) == i ) then
          x(i) = b%values(j,1)
          exit ! j loop
        end if
        if ( b%r2(j) > i ) exit ! j loop
      end do ! j
    case ( M_Full )
      do i = 1, n
        x(i) = b%values(i,i)
      end do
    end select
    if ( present(squareRoot) ) then
      if ( squareRoot ) then
        do i = 1, n
          if ( x(i) < 0.0_r8 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Negative diagonal element in GetDiagonal_0 and SquareRoot is true" )
          x(i) = sqrt(x(i))
        end do
      end if
    end if
  end subroutine GetDiagonal_0

  ! -----------------------------------------  GetMatrixElement_0  -----
  real(r8) function GetMatrixElement_0 ( Matrix, Row, Col )
  ! Get the (row,col) element of Matrix
    type(matrixElement_T), intent(in) :: Matrix
    integer, intent(in) :: Row, Col
    integer :: J
    select case ( matrix%kind )
    case ( m_absent )
      getMatrixElement_0 = 0.0_r8
    case ( m_banded )
      if ( row < matrix%r1(col) .or. &
        &  row >= matrix%r1(col) + matrix%r2(col)-matrix%r2(col-1) ) then
        getMatrixElement_0 = 0.0_r8
      else
        getMatrixElement_0 = matrix%values(row - matrix%r1(col) + &
          &                                  matrix%r2(col-1) +1, 1)
      end if
    case ( m_column_sparse )
      do j = matrix%r1(col-1)+1, matrix%r1(col)
        if ( row == matrix%r2(j) ) then
          getMatrixElement_0 = matrix%values(j,1)
          return
        end if
      end do
      getMatrixElement_0 = 0.0_r8
    case ( m_full )
      getMatrixElement_0 = matrix%values(row,col)
    end select
  end function GetMatrixElement_0

  ! --------------------------------------  GetVectorFromColumn_0  -----
  subroutine GetVectorFromColumn_0 ( B, Column, X )
  ! Fill the vector X from "Column" of B.
    type(MatrixElement_T), intent(in) :: B
    integer, intent(in) :: Column
    real(r8), dimension(:), intent(out) :: X

    if ( column < 1 .or. column > b%ncols ) call MLSMessage ( MLSMSG_Error, &
      & moduleName, '"Column" is out-of-range in GetVectorFromColumn_0' )
    if ( size(x) < b%nrows )  call MLSMessage ( MLSMSG_Error, &
      & moduleName, '"X" is too small in GetVectorFromColumn_0' )
    select case ( b%kind )
    case ( M_Absent )
      x = 0.0_r8
    case ( M_Banded )
      x = 0.0_r8
      x(b%r1(column):b%r1(column)+b%r2(column)-b%r2(column-1)-1) = &
        b%values(b%r2(column-1)+1:b%r2(column),column)
    case ( M_Column_Sparse )
      x = 0.0_r8
      x(b%r2(b%r1(column-1)+1:b%r1(column))) = &
        b%values(b%r1(column-1)+1:b%r1(column),column)
    case ( M_Full )
      x(1:b%nrows) = b%values(1:b%nrows,column)
    end select
  end subroutine GetVectorFromColumn_0

  ! ------------------------------------------------  MaxAbsVal_0  -----
  real(r8) function MaxAbsVal_0 ( B )
  ! Return the magnitude of the element in B that has the largest magnitude.
    type(MatrixElement_T), intent(in) :: B
    if ( b%kind == m_absent ) then
      maxAbsVal_0 = 0.0
    else
      maxAbsVal_0 = maxval(abs(b%values))
    end if
  end function MaxAbsVal_0

  ! --------------------------------------------------  MinDiag_0  -----
  real(r8) function MinDiag_0 ( B )
  ! Return the magnitude of the element on the diagonal of B that has the
  ! smallest magnitude.
    type(MatrixElement_T), intent(in) :: B
    integer :: I, J
    select case ( b%kind )
    case ( M_Absent )
      minDiag_0 = 0.0_r8
    case ( M_Banded )
      minDiag_0 = huge(0.0_r8)
      do i = 1, min(b%nrows,b%ncols)
        if ( b%r1(i) <= i .and. b%r1(i) + b%r2(i) - b%r2(i-1) > i ) &
          & minDiag_0 = min(minDiag_0, abs(b%values(b%r2(i-1)+i-b%r1(i)+1,1)))
      end do ! i
    case ( M_Column_Sparse )
      minDiag_0 = huge(0.0_r8)
      do i = 1, min(b%nrows,b%ncols)
        do j = b%r1(i-1)+1, b%r1(i)
          if ( b%r2(j) == i ) then
            minDiag_0 = min(minDiag_0, abs(b%values(j,1)))
            exit ! j loop
          end if
          if ( b%r2(j) > i ) exit ! j loop
        end do ! j
      end do ! i
    case ( M_Full )
      minDiag_0 = abs(b%values(1,1))
      do i = 2, min(b%nrows,b%ncols)
        minDiag_0 = min(minDiag_0, abs(b%values(i,i)))
      end do ! i
    end select
  end function MinDiag_0

  ! ---------------------------------------  MultiplyMatrixBlocks  -----
  subroutine MultiplyMatrixBlocks ( XB, YB, ZB, UPDATE, SUBTRACT, XMASK, YMASK, &
    &                               UPPER )
  ! ZB = XB^T YB if UPDATE is absent or false and SUBTRACT is absent or false;
  ! ZB = -XB^T YB if UPDATE is absent or false and SUBTRACT is present and
  !      true;
  ! ZB = ZB + XB^T YB if UPDATE is present and true and  SUBTRACT is absent
  !      or false;
  ! ZB = ZB - XB^T YB if UPDATE is present and true and  SUBTRACT is present
  !      and true;
  ! If XMASK (resp. YMASK) is present and associated, ignore columns of XB
  ! (resp. YB) that correspond to nonzero bits of XMASK (resp. YMASK).
  ! If UPPER is present and true, compute only the upper triangle of ZB.
    type(MatrixElement_T), intent(in) :: XB, YB
    type(MatrixElement_T), intent(inout) :: ZB
    logical, intent(in), optional :: UPDATE
    logical, intent(in), optional :: SUBTRACT
    integer, optional, pointer, dimension(:) :: XMASK, YMASK ! intent(in)
    logical, intent(in), optional :: UPPER

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! If UPDATE is absent or false, it is important to invoke DestroyBlock
  ! for ZB after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K, L, M, MZ, N
    integer :: XI_1, XI_N, XR_1, XR_N, YI_1, YI_N, YR_1, YR_N, CR_1, CR_N, C_N
    logical :: MY_SUB, MY_UPD, MY_UPPER
    real(r8) :: S                            ! Sign, subtract => -1 else +1
    integer :: XD, YD
    integer, pointer, dimension(:) :: XM, YM
    real(r8) :: XY                           ! Product of columns of X and Y
    real(r8), pointer, dimension(:,:) :: Z   ! Temp for sparse * sparse

    nullify ( xm, ym, z )
    my_upd = .false.
    if ( present(update) ) my_upd = update
    my_sub = .false.
    if ( present(subtract) ) my_sub = subtract
    s = 1.0_r8
    if ( my_sub ) s = -1.0_r8
    if ( present(xmask) ) xm => xmask
    if ( present(ymask) ) ym => ymask
    my_upper = .false.
    if ( present(upper) ) my_upper = upper
    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      if ( .not. my_upd) call createBlock ( zb, xb%nCols, yb%nCols, M_Absent )
      return
    end if
    if ( xb%nrows /= yb%nrows ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "XB and YB Matrix sizes incompatible in Multiply_Matrix_Blocks" )
    if ( my_upd ) then
      if ( zb%kind /= m_absent ) then
        if ( xb%ncols /= zb%nrows .or. yb%ncols /= zb%ncols ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "ZB Matrix size incompatible in Multiply_Matrix_Blocks" )
      else
        zb%nrows = xb%ncols
        zb%ncols = yb%ncols
      end if
    else
      zb%nrows = xb%ncols
      zb%ncols = yb%ncols
    end if
    select case ( xb%kind )
    case ( M_Banded )
      select case ( yb%kind )
      case ( M_Banded )       ! XB banded, YB banded
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, zb%nrows, zb%ncols, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
        if ( update .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, zb%ncols    ! Columns of Z = columns of YB
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          yi_1 = yb%r2(j-1)+1   ! index of 1st /=0 el. in this col of yb
          yi_n = yb%r2(j)       ! index of last /=0 el. in this col of yb
          yr_1 = yb%r1(j)       ! first row index of yb
          yr_n = yr_1 + yi_n - yi_1 ! last row index of yb
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            ! Inner product of column I of XB with column J of YB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            xi_1 = xb%r2(i-1)+1  ! index of 1st /=0 el. in this col of xb
            xi_n = xb%r2(i)      ! index of last /=0 el. in this col of xb
            xr_1 = xb%r1(i)      ! first row index of xb
            xr_n = xr_1 + xi_n - xi_1   ! last row index of xb

            ! Now work out what they have in common
            cr_1 = max ( xr_1, yr_1 )
            cr_n = min ( xr_n, yr_n )
            c_n = cr_n - cr_1 + 1
            if ( c_n <= 0 ) cycle

            ! Now make xd and yd point to the starts of the common parts
            xd = xi_1 + cr_1 - xr_1
            yd = yi_1 + cr_1 - yr_1

!           xy = dot_product ( xb%values(xd:xd+c_n-1,1), &
!                          &   yb%values(yd:yd+c_n-1,1) )
            xy = dot( c_n, xb%values(xd,1), 1, yb%values(yd,1), 1 )
            z(i,j) = z(i,j) + s * xy
          end do ! i
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Column_sparse ) ! XB banded, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded X sparse in Multiply_Matrix_Blocks", ModuleName )
        if ( update .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, zb%ncols    ! Columns of Z
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            k = xb%r1(i)      ! Row subscript of first nonzero in XB's column I
            l = xb%r2(i-1)+1  ! Position in XB%VALUES of it
            n = yb%r1(j-1)+1  ! Position in YB%R2 of row subscript in YB
            m = yb%r2(n)      ! Row subscript of nonzero in YB's column J
            do while ( l <= xb%r2(i) .and. n <= yb%r1(j) )
              if ( k < m ) then
                l = l + 1
                k = k + 1
              else if ( k > m ) then
                n = n + 1
                if ( n > yb%r1(j) ) exit
                m = yb%r2(n)
              else
                ! Multiplying by S is faster than testing my_sub
                z(i,j) = z(i,j) + s * xb%values(l,1) * yb%values(n,1)
                l = l + 1
                k = k + 1
                n = n + 1
                if ( n > yb%r1(j) ) exit
                m = yb%r2(n)
              end if
            end do
          end do ! i
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Full )         ! XB banded, YB full
        if ( zb%kind /= m_full ) then
          call allocate_test ( z, xb%ncols, yb%ncols, &
            & "Z for banded X full in Multiply_Matrix_Blocks", ModuleName )
          if ( my_upd ) call densify ( z, zb )
          call createBlock ( zb, xb%ncols, yb%ncols, M_Full, novalues=.true. )
          zb%values => z
        end if
        if ( .not. my_upd ) zb%values = 0.0_r8
        do i = 1, xb%ncols    ! Rows of ZB
          if ( associated(xm) ) then
            if ( btest(xm(i/b+1),mod(i,b)) ) cycle
          end if
          m = xb%r1(i)        ! Index of first row of XB with nonzero value
          k = xb%r2(i-1) + 1
          l = xb%r2(i)
          if ( l < k ) cycle  ! Empty column in XB
          mz = 1
          if ( my_upper ) mz = i
          do j = mz, yb%ncols  ! Columns of ZB
            if ( associated(ym) ) then
              if ( btest(ym(j/b+1),mod(j,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
!           xy = dot_product( xb%values(k:l,1), yb%values(m:m+l-k,j) )
            xy = dot( l-k+1, xb%values(k,1), 1, yb%values(m,j), 1 )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! j
        end do ! i
      end select
    case ( M_Column_sparse )
      select case ( yb%kind )
      case ( M_Banded )       ! XB column-sparse, YB banded
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for sparse X banded in Multiply_Matrix_Blocks", ModuleName )
        if ( update .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, zb%ncols    ! Columns of Z
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            l = xb%r1(i-1)+1  ! Position in XB%R2 of row subscript in XB
            k = xb%r2(l)      ! Row subscript of nonzero in XB's column I
            m = yb%r1(j)      ! Row subscript of first nonzero in YB's column J
            n = yb%r2(j-1)+1  ! Position in YB%VALUES of it
            do while ( l <= xb%r2(i) .and. n <= yb%r1(j) )
              if ( k < m ) then
                l = l + 1
                if ( l > xb%r1(i) ) exit
                k = xb%r2(l)
              else if ( k > m ) then
                n = n + 1
                m = m + 1
              else
                ! Multiplying by S is faster than testing my_sub
                z(i,j) = z(i,j) + s * xb%values(l,1) * yb%values(n,1)
                l = l + 1
                if ( l > xb%r1(i) ) exit
                k = xb%r2(l)
                n = n + 1
                m = m + 1
              end if
            end do
          end do ! i
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Column_sparse ) ! XB column-sparse, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for sparse X sparse in Multiply_Matrix_Blocks", ModuleName )
        if ( update .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, zb%ncols    ! Columns of Z
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            l = xb%r1(i-1)+1  ! Position in XB%R2 of row subscript in XB
            k = xb%r2(l)      ! Row subscript of nonzero in XB's column I
            n = yb%r1(j-1)+1  ! Position in YB%R2 of row subscript in YB
            m = yb%r2(n)      ! Row subscript of nonzero in YB's column J
            z(i,j) = 0.0_r8
            do while ( l <= xb%r2(i) .and. n <= yb%r1(j) )
              if ( k < m ) then
                l = l + 1
                if ( l > xb%r1(i) ) exit
                k = xb%r2(l)
              else if ( k > m ) then
                n = n + 1
                if ( n > yb%r1(j) ) exit
                m = yb%r2(n)
              else
                ! Multiplying by S is faster than testing my_sub
                z(i,j) = z(i,j) + s * xb%values(l,1) * yb%values(n,1)
                l = l + 1
                if ( l > xb%r1(i) ) exit
                k = xb%r2(l)
                n = n + 1
                if ( n > yb%r1(j) ) exit
                m = yb%r2(n)
              end if
            end do
          end do ! i
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Full )         ! XB column-sparse, YB full
        if ( zb%kind /= m_full ) then
          call allocate_test ( z, xb%ncols, yb%ncols, &
            & "Z for sparse X full in Multiply_Matrix_Blocks", ModuleName )
          if ( my_upd ) call densify ( z, zb )
          call createBlock ( zb, xb%ncols, yb%ncols, M_Full, novalues=.true. )
          zb%values => z
        end if
        if ( .not. my_upd ) zb%values = 0.0_r8
        do j = 1, zb%ncols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            k = xb%r1(i-1)+1
            l = xb%r1(i)
            ! Inner product of column I of XB with column J of YB
!           xy = dot_product( xb%values(k:l,1), yb%values(xb%r2(k:l),j) )
            xy = dot( l-k+1, xb%values(k,1), 1, yb%values(xb%r2(k),j), 1 )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
        end do ! j
      end select
    case ( M_Full )
      if ( zb%kind /= m_full ) then
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for full X <anything> in Multiply_Matrix_Blocks", ModuleName )
        if ( my_upd ) call densify ( z, zb )
        call createBlock ( zb, xb%ncols, yb%ncols, M_Full, novalues=.true. )
        zb%values => z
      end if
      if ( .not. my_upd ) zb%values = 0.0_r8
      select case ( yb%kind )
      case ( M_Banded )       ! XB full, YB banded
        do j = 1, zb%ncols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          m = yb%r1(j)        ! Index of first row of YB with nonzero value
          k = yb%r2(j-1)+1    ! K and L are indices of YB
          l = yb%r2(j)
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
!           xy = dot_product( xb%values(m:m+l-k,i), yb%values(k:l,1))
            xy = dot( l-k+1, xb%values(m,i), 1, yb%values(k,1), 1 )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
        end do ! j
      case ( M_Column_sparse ) ! XB full, YB column-sparse
        do j = 1, zb%ncols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( btest(ym(j/b+1),mod(j,b)) ) cycle
          end if
          k = yb%r1(j-1)+1    ! K and L are indices of YB
          l = yb%r1(j)
          mz = zb%nrows
          if ( my_upper ) mz = j
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( btest(xm(i/b+1),mod(i,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
            xy = dot_product( xb%values(yb%r2(k:l),i), yb%values(k:l,1) )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
        end do ! j
      case ( M_Full )         ! XB full, YB full
        if ( associated(xm) .or. associated(ym) ) then
          do j = 1, zb%ncols  ! Columns of ZB
            if ( associated(ym) ) then
              if ( btest(ym(j/b+1),mod(j,b)) ) cycle
            end if
            mz = zb%nrows
            if ( my_upper ) mz = j
            do i = 1, mz      ! Rows of Z = columns of XB
              if ( associated(xm) ) then
                if ( btest(xm(i/b+1),mod(i,b)) ) cycle
              end if
!             xy = dot_product(xb%values(1:xb%nrows,i), yb%values(1:xb%nrows,j))
              xy = dot( xb%nrows, xb%values(1,i), 1, yb%values(1,j), 1 )
              zb%values(i,j) = zb%values(i,j) + s * xy
            end do ! i = 1, xb%ncols
          end do ! j = 1, yb%ncols
        else if ( my_upper ) then
          do j = 1, zb%ncols  ! Columns of ZB
            do i = 1, j       ! Rows of Z = columns of XB
!             xy = dot_product(xb%values(1:xb%nrows,i), yb%values(1:xb%nrows,j))
              xy = dot(xb%nrows, xb%values(1,i), 1, yb%values(1,j), 1 )
              zb%values(i,j) = zb%values(i,j) + s * xy
            end do
          end do
        else
          zb%values = zb%values + s * matmul(transpose(xb%values),yb%values)
        end if
      end select
    end select
  end subroutine MultiplyMatrixBlocks

  ! -------------------------------------  MultiplyMatrixVector_0  -----
  subroutine MultiplyMatrixVector_0 ( B, V, P, UPDATE, SUBTRACT )
  ! P = B^T V if UPDATE is absent or false.
  ! P = P + B^T V if UPDATE is present and true and SUBTRACT is absent or false.
  ! P = P - B^T V if UPDATE is present and true and SUBTRACT is present and true.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE
    logical, optional, intent(in) :: SUBTRACT

    real(r8) :: BV                 ! Product of a column of B and the vector V
    integer :: I, M, N             ! Subscripts and loop inductors
    logical :: My_Sub, My_update
    real(r8) :: S                  ! SUBTRACT => -1 else +1
    integer :: V1                  ! Subscripts and loop inductors

    if ( b%nrows /= size(v) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and vector not compatible in MultiplyMatrixVector_0" )
    if ( b%ncols /= size(p) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and result not compatible in MultiplyMatrixVector_0" )
    my_update = .false.
    if ( present(update) ) my_update = update
    my_sub = .false.
    if ( present(subtract) ) my_sub = subtract
    s = 1.0_r8
    if ( my_sub ) s = -1.0_r8
    if ( .not. my_update ) p = 0.0_r8
    select case ( b%kind )
    case ( M_Absent )
    case ( M_Banded )
      do i = 1, size(p)
        v1 = b%r2(i-1)             ! starting position in B%VALUES - 1
        n = b%r2(i) - v1           ! how many values
        if ( n > 0 ) then          ! Because v1+1 will be out-of-range if
                                   ! there is a final zero-size column
          m = b%r1(i)              ! starting position in V
!         bv = dot_product(b%values(v1+1:v1+n,1), v(m:m+n-1))
          bv = dot(n, b%values(v1+1,1), 1, v(m), 1)
          p(i) = p(i) + s * bv
        end if
      end do ! i
     case ( M_Column_Sparse )
      do i = 1, size(p)
        bv = 0.0_r8
        do n = b%r1(i-1)+1, b%r1(i)
          bv = bv + b%values(n,1) * v(b%r2(n))
        end do ! n
        p(i) = p(i) + s * bv
      end do ! i
    case ( M_Full )
      do i = 1, size(p)
        bv = dot(size(v), b%values(1,i), 1, v(1), 1)
        p(i) = p(i) + s * bv
      end do ! i
    end select
  end subroutine MultiplyMatrixVector_0

  ! ----------------------------------  MultiplyMatrixVectorNoT_0  -----
  subroutine MultiplyMatrixVectorNoT_0 ( B, V, P, UPDATE, DoDiag, SUBTRACT )
  ! P = B V if UPDATE is absent or false.
  ! P = P + B V if UPDATE is present and true and SUBTRACT is absent or false
  ! P = P - B V if UPDATE is present and true and SUBTRACT is present and true
  ! Don't multiply by the diagonal element if doDiag (default true) is
  ! present and false.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE, DoDiag, SUBTRACT

    integer :: I, J, M, N          ! Subscripts and loop inductors
    logical :: My_diag, My_sub, My_update
    real(r8) :: SIGN               ! Multiplying by sign is faster than testing
    integer :: V1                  ! Subscripts and loop inductors

    if ( b%ncols /= size(v) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and vector not compatible in MultiplyMatrixVectorNoT_0" )
    if ( b%nrows /= size(p) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and result not compatible in MultiplyMatrixVectorNoT_0" )
    my_update = .false.
    if ( present(update) ) my_update = update
    my_diag = .true.
    if ( present(doDiag) ) my_diag = doDiag
    my_sub = .false.
    if ( present(subtract) ) my_sub = subtract
    sign = 1.0_r8
    if ( my_sub ) sign = -1.0_r8
    if ( .not. my_update ) p = 0.0_r8
    select case ( b%kind )
    case ( M_Absent )
    case ( M_Banded )
      do j = 1, size(v)            ! columns
        v1 = b%r2(j-1)             ! (starting position in B%VALUES) - 1
        n = b%r2(j) - v1           ! how many values
        m = b%r1(j)                ! starting row subscript in B%VALUES
        if ( my_diag ) then        ! do the whole matrix
          do i = m, m+n-1          ! rows
            p(i) = p(i) + sign * b%values(v1+i-m+1,1) * v(j)
          end do ! i = 1, n
        else                       ! skip the diagonal
          do i = m, min(m+n-1,j-1) ! rows
            p(i) = p(i) + sign * b%values(v1+i-m+1,1) * v(j)
          end do ! i = m, min(m+n-1,j-1)
          do i = max(m,j+1), m+n-1 ! rows
            p(i) = p(i) + sign * b%values(v1+i-m+1,1) * v(j)
          end do ! i = m, min(m,n-1,j-1)
        end if
      end do ! j
     case ( M_Column_Sparse )
      do j = 1, size(p)            ! columns
        do n = b%r1(j-1)+1, b%r1(j)! rows
          i = b%r2(n)              ! row number
          if ( i/=j .or. my_diag ) &
            & p(i) = p(i) + sign * b%values(n,1) * v(j)
        end do ! n
      end do ! j
    case ( M_Full )
      if ( my_diag ) then          ! do the whole matrix
        do i = 1, size(p)
!           p(i) = p(i) + dot(size(v), b%values(i,1), size(b%values,1), v(1), 1)
          p(i) = p(i) + sign * dot_product ( b%values(i,:), v )
        end do ! i
      else                         ! skip the diagonal
        do i = 1, size(p)
!           p(i) = p(i) + dot(i-1, b%values(i,1), size(b%values,1), v(1), 1)
!           p(i) = p(i) + &
!             & dot(size(v)-i, b%values(i,i+1), size(b%values,1), v(i+1), 1)
          p(i) = p(i) + sign * dot_product ( b%values(i, 1:i-1), v(1:i-1) ) + &
            & sign * dot_product ( b%values( i, i+1:), v(i+1:) )
        end do ! i
      end if
    end select
  end subroutine MultiplyMatrixVectorNoT_0

  ! -------------------------------------  NewMultiplyMatrixBlocks  ----
  function NewMultiplyMatrixBlocks ( X, Y ) result ( Z ) ! Z = X^T Y
    type(MatrixElement_T), intent(in) :: X, Y
    type(MatrixElement_T) :: Z

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    call MultiplyMatrixBlocks ( x, y, z )
  end function NewMultiplyMatrixBlocks

  ! -----------------------------------  NewMultiplyMatrixVector_0  ----
  function NewMultiplyMatrixVector_0 ( B, V ) result ( P ) ! P = B^T V
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(size(v)) :: P
    call MultiplyMatrixVector ( b, v, p, .false. )
  end function NewMultiplyMatrixVector_0

  ! -------------------------------------------------  RowScale_0  -----
  subroutine RowScale_0 ( V, X, NEWX ) ! Z = V X where V is a diagonal
  !                                     matrix represented by a vector and
  !                                     Z is either X or NEWX.
    real(r8), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(inout), target :: X
    type(MatrixElement_T), intent(out), target, optional :: NEWX
    type(MatrixElement_T), pointer :: Z

    integer :: I

    if ( present(newx) ) then
      z => newx
      call copyBlock ( z, x )
    else
      z => x
    end if
    select case ( z%kind )
    case ( M_Absent )
    case ( M_Banded )
      do i = 1, z%ncols
        z%values(z%r2(i-1)+1:z%r2(i),1) = &
          & v(z%r1(i):z%r1(i)+z%r2(i)-z%r2(i-1)-1) * &
          & x%values(z%r2(i-1)+1:z%r2(i),1)
      end do
    case ( M_Column_Sparse )
      do i = 1, z%ncols
        z%values(z%r1(i-1)+1:z%r1(i),1) = v(z%r2(z%r1(i-1)+1:z%r1(i))) * &
          & x%values(z%r1(i-1)+1:z%r1(i),1)
      end do
    case ( M_Full )
      do i = 1, z%nrows
        z%values(i,:) = v(i) * x%values(i,:)
      end do
    end select
  end subroutine RowScale_0

  ! -------------------------------------------------  ScaleBlock  -----
  subroutine ScaleBlock ( Z, A )        ! Z := A * Z, where A is scalar
    type(matrixElement_T), intent(inout) :: Z
    real(r8), intent(in) :: A
    if ( z%kind /= m_absent ) z%values = a * z%values
  end subroutine ScaleBlock

  ! -------------------------------------------  SolveCholeskyM_0  -----
  subroutine SolveCholeskyM_0 ( U, X, B, TRANSPOSE )
  ! Solve the system U X = B or U^T X = B for X, depending on TRANSPOSE,
  ! where U is known to be upper-triangular.  X may be the same as B.
  ! B may be absent, in which case the right-hand side is in X on input,
  ! and the solution replaces it on output.
    type(MatrixElement_T), intent(in) :: U        ! Must be square
    type(MatrixElement_T), intent(inout), target :: X
    type(MatrixElement_T), intent(in), target, optional :: B
    logical, intent(in), optional :: TRANSPOSE    ! Solve U^T X = B if
    !                                               present and true.

    real(r8) :: D        ! Diagonal element of U
    integer :: I, J, K   ! Subscripts and loop inductors
    type(MatrixElement_T), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    integer :: NC        ! Number of columns in B, not necessarily == N
    real(r8), parameter :: TOL = tiny(0.0_r8)
    real(r8), pointer, dimension(:,:) :: XS  ! The solution, dense
    real(r8), pointer, dimension(:,:) :: UD  ! U, densified

    nullify ( ud, xs )
    n = u%nrows
    if ( n /= u%nCols ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "U matrix in SolveCholeskyM_0 must be square" )
    my_b => x
    if ( present(b) ) my_b => b
    if ( n /= my_b%nrows ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "B matrix not compatible with U matrix in SolveCholeskyM_0" )
    my_t = .false.
    if ( present(transpose) ) my_t = transpose
    nc = my_b%nCols

    call allocate_test ( xs, my_b%nRows, my_b%nCols, "XS in SolveCholeskyM_0", &
      & ModuleName )
    call densify ( xs, my_b )
    if ( my_t ) then ! solve U^T X = B for X
      select case ( u%kind )
      case ( M_Absent )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "U matrix in SolveCholeskyM_0 must not be absent" )
      case ( M_Banded )
        do i = 1, n
          if ( u%r1(i) + u%r2(i) - u%r2(i-1) - 1 /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is not triangular" )
          d = u%values(u%r2(i),1)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          do j = 1, nc
!           xs(i,j) = ( xs(i,j) - &
!                   &   dot_product(u%values(u%r2(i-1)+1:u%r2(i)-1,1), &
!                   &               xs(u%r1(i):i-1,j) ) ) / d
            xs(i,j) = ( xs(i,j) - &
                    &   dot( u%r2(i)-u%r2(i-1)-1, u%values(u%r2(i-1)+1,1), 1, &
                    &               xs(u%r1(i),j), 1 ) ) / d
          end do ! j = 1, nc
        end do ! i = 1, n
      case ( M_Column_Sparse )
        do i = 1, n
          if ( u%r2(u%r1(i)) /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is not triangular" )
          d = u%values(u%r1(i),1)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          do j = 1, nc
            do k = u%r1(i-1)+1, u%r1(i)-1
              xs(i,j) = xs(i,j) - u%values(k,1) * xs(u%r2(k),j)
            end do ! k = u%r1(i-1)+1, u%r1(i)-1
            xs(i,j) = xs(i,j) / d
          end do ! j = 1, nc
        end do ! i = 1, n
      case ( M_Full )
        do i = 1, n
          d = u%values(i,i)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          do j = 1, nc
!           xs(i,j) = ( xs(i,j) - &
!                   &   dot_product(u%values(1:i-1,i),xs(1:i-1,j) ) ) / d
            xs(i,j) = ( xs(i,j) - &
                    &   dot( i-1, u%values(1,i), 1, xs(1,j), 1) ) / d
          end do ! j = 1, nc
        end do ! i = 1, n
      end select
    else             ! solve U X = B for X
      if ( u%kind == M_full ) then
        ud => u%values
      else
        call allocate_test ( ud, n, n, "UD in SolveCholeskyM_0", ModuleName )
        call densify ( ud, u )
      end if
      d = ud(n,n)
      if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "U matrix in SolveCholeskyM_0 is singular" )
      xs(n,1:nc) = xs(n,1:nc) / d
      do i = n-1, 1, -1
        d = ud(i,i)
        if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "U matrix in SolveCholeskyM_0 is singular" )
        do j = 1, nc
!         xs(i,j) = ( xs(i,j) - &
!                 &   dot_product(ud(i,i+1:n),xs(i+1:n,j) ) ) / d
          xs(i,j) = ( xs(i,j) - &
                  &   dot( n-i, ud(i,i+1), size(ud,1), xs(i+1,j), 1 ) ) / d
        end do ! j = 1, nc
      end do ! i = 1, n
      if ( u%kind /= M_Full ) &
        & call deallocate_test ( ud, "UD in SolveCholeskyM_0", ModuleName )
    end if ! my_t
    call sparsify ( xs, x, "XS in SolveCholeskyM_0", ModuleName ) ! X := Xs
  end subroutine SolveCholeskyM_0

  ! -------------------------------------------  SolveCholeskyV_0  -----
  subroutine SolveCholeskyV_0 ( U, X, B, TRANSPOSE )
  ! Solve the system U X = B or U^T X = B for X, depending on TRANSPOSE,
  ! where U is known to be upper-triangular.  X may be the same as B.
  ! B may be absent, in which case the right-hand side is in X on input,
  ! and the solution replaces it on output.  The arrays X and B are
  ! two-dimensional sections of subvectors of objects of type Vector_T.
  ! Their elements are taken to correspond to the rows of U in array
  ! element order.
    type(MatrixElement_T), intent(in) :: U        ! Must be square
    real(r8), dimension(:), intent(inout), target :: X
    real(r8), dimension(:), intent(in), target, optional :: B
    logical, intent(in), optional :: TRANSPOSE    ! Solve U^T X = B if
    !                                               present and true.

    real(r8) :: D        ! Diagonal element of U
    integer :: H, I      ! Subscripts and loop inductors
    real(r8), dimension(:), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    real(r8), parameter :: TOL = tiny(0.0_r8)
    real(r8), dimension(:,:), pointer :: UD  ! U, densified

    n = u%nrows
    if ( n /= u%nCols ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "U matrix in SolveCholeskyV_0 must be square" )
    my_b => x
    if ( present(b) ) my_b => b
    if ( n /= size(my_b) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "B matrix not compatible with U matrix in SolveCholeskyV_0" )
    my_t = .false.
    if ( present(transpose) ) my_t = transpose

    if ( my_t ) then ! solve U^T X = B for X
      select case ( u%kind )
      case ( M_Absent )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "U matrix in SolveCholeskyV_0 must not be absent" )
      case ( M_Banded )
        do i = 1, n
          if ( u%r1(i) + u%r2(i) - u%r2(i-1) - 1 /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is not triangular" )
          d = u%values(u%r2(i),1)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          ! dot_product( u%values(u%r2(i-1)+1:u%r2(i)-1,1), &
          ! &            b(u%r1(i):u%r1(i)+u%r2(i)-u%r2(i-1)-1) )
          x(i) = ( my_b(i) - dot(u%r2(i)-u%r2(i-1)-1, &
            &      u%values(u%r2(i-1)+1,1), 1, my_b(u%r1(i)), 1) ) / d
        end do ! i = 1, n
      case ( M_Column_Sparse )
        do i = 1, n
          if ( u%r2(u%r1(i)) /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is not triangular" )
          d = u%values(u%r1(i),1)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          do h = u%r1(i-1)+1, u%r1(i)-1
            x(i) = my_b(i) - u%values(h,1) * my_b(u%r2(h))
          end do ! h = u%r1(i-1)+1, u%r1(i)-1
          x(i) = x(i) / d
        end do ! i = 1, n
      case ( M_Full )
        do i = 1, n
          d = u%values(i,i)
          if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          ! dot_product( u%values(1:i-1,i), my_b(1:i-1) )
          x(i) = ( my_b(i) - dot(i-1, u%values(1,i), 1, x(1), 1) ) / d
        end do ! i = 1, n
      end select
    else             ! solve U X = B for X
      if ( u%kind == M_full ) then
        ud => u%values
      else
        nullify ( ud )
        call allocate_test ( ud, n, n, "UD in SolveCholeskyV_0", ModuleName )
        call densify ( ud, u )
      end if
      d = ud(n,n)
      if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "U matrix in SolveCholeskyV_0 is singular" )
      x(n) = my_b(n) / d
      do i = n-1, 1, -1
        d = ud(i,i)
        if ( abs(d) < tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "U matrix in SolveCholeskyV_0 is singular" )
        ! dot_product( ud(i,i+1:n), my_b(i+1:n) )
        x(i) = ( my_b(i) - dot(n-i, ud(i,i+1), size(ud,1), x(i+1), 1) ) / d
      end do ! i = 1, n
      if ( u%kind /= M_Full ) &
        & call deallocate_test ( ud, "UD in SolveCholeskyV_0", ModuleName )
    end if ! my_t
  end subroutine SolveCholeskyV_0

  ! ---------------------------------------------------  Sparsify  -----
  subroutine Sparsify ( Z, B, Why, CallingModule )
  ! Given an array Z, compute its sparse representation and store it
  ! in the matrix block B.
    real(r8), pointer :: Z(:,:)              ! Full array of values
    type(MatrixElement_T), intent(out) :: B  ! Z as a block, maybe sparse
    character(len=*), intent(in), optional :: Why
    character(len=*), intent(in), optional :: CallingModule
    ! If either Why or CallingModule is present, Z is deallocated using
    ! Deallocate_Test

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the B argument of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I1, I2              ! Row indices in Z
    integer :: J                   ! Column index in Z
    integer :: KIND                ! Representation to use for B
    integer :: NNZ                 ! Number of nonzeroes in Z
    integer :: NNZC(size(z,2))     ! Number of nonzeroes in a column of Z
    integer :: R1(size(z,2))       ! Row number of first nonzero in a column
    real(r8), save :: SQ_EPS = -1.0_r8  ! sqrt(epsilon(1.0_r8))
    real(r8) :: ZT(size(z,2))      ! Maximum value in a column of Z, then
      ! max(sqrt(sq_eps*zt),tiny(1.0_r8)).  Elements less than this threshold
      ! in magnitude are considered to be zero.

    if ( sq_eps < 0.0_r8 ) sq_eps = sqrt(epsilon(1.0_r8))
    do j = 1, size(z,2)
    ! zt(j) = maxval(abs(z(:,j)))
    ! zt(j) = max(sq_eps*zt(j), sqrt(tiny(1.0_r8)))
    ! nnzc(j) = count(abs(z(:,j)) < zt(j))
      zt(j) = 0.0_r8
      do i1 = 1, size(z,1)
        zt(j) = max(zt(j),abs(z(i1,j)))
      end do ! i1
      zt(j) = max(sq_eps*zt(j), sqrt(tiny(1.0_r8)))
      nnzc(j) = 0
      do i1 = 1, size(z,1)
        if ( abs(z(i1,j)) >= zt(j) ) nnzc(j) = nnzc(j) + 1
      end do ! i1
    end do ! j
    nnz = sum(nnzc)
    if ( nnz == 0 ) then ! Empty
      call createBlock ( b, size(z,1), size(z,2), M_Absent )
    else if ( nnz <= int(sparsity * size(z)) ) then ! sparse
      kind = M_Banded
      do j = 1, size(z,2)
        do i1 = 1, size(z,1)       ! Find row number of first nonzero
          if ( abs(z(i1,j)) > zt(j) ) exit
        end do
        r1(j) = i1
        do i2 = size(z,1), 1, -1   ! Find row number of last nonzero
          if ( abs(z(i2,j)) > zt(j) ) exit
        end do
        if ( int(col_sparsity*(i2 - i1)) > nnzc(j) ) then
          kind = M_Column_Sparse
          do i1 = 1, j-1 ! I1 is a column number in this case
            nnzc(i1) = count(abs(z(:,i1)) <= zt(i1))
          end do ! i1
          exit
        end if
        nnzc(j) = max(i2 - i1 + 1,0)
      end do ! j
      if ( kind == M_Banded ) then
        call createBlock ( b, size(z,1), size(z,2), M_Banded, sum(nnzc) )
        b%r1 = r1        ! Row number of first nonzero in the column
        do j = 1, size(z,2)
          b%r2(j) = nnzc(j) + b%r2(j-1) ! Subscript of last nonzero in the
                                        ! column
          b%values(b%r2(j-1)+1:b%r2(j),1) = &
            & z(b%r1(j):b%r1(j)+b%r2(j)-b%r2(j-1)-1,j)
        end do
      else
        call createBlock ( b, size(z,1), size(z,2), M_Column_Sparse, nnz )
        i1 = 0
        do j = 1, size(z,2)
          do i2 = 1, size(z,1)
            if ( abs(z(i2,j)) > zt(j) ) then
              i1 = i1 + 1
              b%values(i1,1) = z(i2,j)
              b%r2(i1) = i2
            end if
            b%r1(j) = i1
          end do ! i2
        end do ! j
      end if
    else ! full
      call createBlock ( b, size(z,1), size(z,2), M_Full )
      b%values = z
    end if
    if ( present(why) ) then
      if ( present(callingModule) ) then
        call deallocate_test ( z, why, callingModule )
      else
        call deallocate_test ( z, why, "No module specified" )
      end if
    else if ( present(callingModule) ) then
      call deallocate_test ( z, "No variable specified", callingModule )
    end if
  end subroutine Sparsify

  ! -------------------------------------------  UpdateDiagonal_0  -----
  subroutine UpdateDiagonal_0 ( A, LAMBDA )
  ! Add LAMBDA to the diagonal of A
    type(MatrixElement_T), intent(inout) :: A
    real(r8), intent(in) :: LAMBDA

    integer :: I, J                          ! Subscripts and loop inductors
    integer :: N                             ! min(a%ncols,a%nrows)
    integer :: NCols, NRows                  ! Copies of a%...
    real(r8), dimension(:,:), pointer :: T   ! A temporary dense matrix

    nullify ( t )
    n = min(a%ncols,a%nrows)
    select case ( a%kind )
    case ( m_absent )
      ncols = a%nCols ! because the first argument of CreateBlock is intent(out)
      nrows = a%nRows
      call createBlock ( a, nRows, nCols, m_banded, n )
      do i = 1, n
        a%r1(i) = i
        a%r2(i) = i
        a%values(i,1) = lambda
      end do
    case ( m_banded )
      do i = 1, n
        if ( i < a%r1(i) .or. i > a%r1(i) + a%r2(i) - a%r2(i-1) - 1 ) then
          ! No diagonal element.  Make room for one by densify-update-sparsify
          call allocate_test ( t, a%nRows, a%nCols, "T in UpdateDiagonal_0", &
            & ModuleName )
          call densify ( t, a )
          call updateDenseDiagonal ( t, lambda, i )
          call sparsify ( t, a, "T in UpdateDiagonal_0", ModuleName ) ! A := T
          call deallocate_test ( t, "T in UpdateDiagonal_0", ModuleName )
          return
        end if
        a%values(a%r2(i-1)+i-a%r1(i)+1,1) = &
          & a%values(a%r2(i-1)+i-a%r1(i)+1,1) + lambda
      end do
    case ( m_column_sparse )
      do i = 1, n
        do j = a%r1(i-1)+1, a%r1(i) ! hunt for the diagonal subscript
          if ( a%r2(j) == i ) then
            a%values(j,1) = a%values(j,1) + lambda
          else if ( a%r2(j) > i ) then
            ! No diagonal element.  Make room for one by densify-update-sparsify
            call allocate_test ( t, a%nRows, a%nCols, "T in UpdateDiagonal_0", &
              & ModuleName )
            call densify ( t, a )
            call updateDenseDiagonal ( t, lambda, i )
            call sparsify ( t, a, "T in UpdateDiagonal_0", ModuleName ) ! A := T
            call deallocate_test ( t, "T in UpdateDiagonal_0", ModuleName )
            return
          end if
        end do
      end do
    case ( m_full )
      call updateDenseDiagonal ( a%values, lambda, 1 )
    end select

  contains
    subroutine UpdateDenseDiagonal ( T, LAMBDA, START )
      real(r8), intent(inout) :: T(:,:)
      real(r8), intent(in) :: LAMBDA
      integer, intent(in) :: START
      integer :: I
      do i = start, n
        t(i,i) = t(i,i) + lambda
      end do
    end subroutine UpdateDenseDiagonal
  end subroutine UpdateDiagonal_0

  ! ----------------------------------------  UpdateDiagonalVec_0  -----
  subroutine UpdateDiagonalVec_0 ( A, X, SUBTRACT, INVERT )
  ! Add X to the diagonal of A if SUBTRACT is absent or false.
  ! Subtract X from the diatonal of A if SUBTRACT is present and true.
  ! If INVERT is present and true, use the inverses of the elements of X.
    type(MatrixElement_T), intent(inout) :: A
    real(r8), intent(in) :: X(:)
    logical, intent(in), optional :: SUBTRACT
    logical, intent(in), optional :: INVERT  ! Update with inverse of X

    integer :: I, J                          ! Subscripts and loop inductors
    logical :: MyInvert
    integer :: N                             ! min(a%ncols,a%nrows)
    integer :: M                             ! max(a%ncols,a%nrows) / n
    integer :: NCols, NRows                  ! Copies of a%...
    real(r8) :: S                            ! Sign to use for X, +1 or -1
    real(r8), dimension(:,:), pointer :: T   ! A temporary dense matrix
    real(r8) :: V                            ! The value to update.  Either
    !                                          S*X or S/X.

    nullify ( t )
    s = 1.0_r8
    if ( present(subtract) ) then
      if ( subtract ) s = -1.0_r8
    end if
    myInvert = .false.
    if ( present(invert) ) myInvert = invert
    n = min(a%ncols,a%nrows)
    m = max(a%ncols,a%nrows) / n
    select case ( a%kind )
    case ( m_absent )
      ncols = a%ncols ! because the first argument of CreateBlock is intent(out)
      nrows = a%nrows
      call createBlock ( a, nRows, nCols, m_banded, m*n )
      do i = 1, n                       ! Loop over shorter side
        j = 1 + m*(i-1)                 ! Index into longer side
        a%r1(i) = j
        a%r2(i) = j + m - 1
        if ( myInvert ) then
          if ( x(i) == 0.0_r8 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
          ! Don't misinterpret the use of r1 and r2 here, it does make sense
          ! see their assignment above.
          a%values(a%r1(i):a%r2(i),1) = s / x(a%r1(i):a%r2(i))
        else
          a%values(a%r1(i):a%r2(i),1) = s * x(a%r1(i):a%r2(i))
        end if
      end do
    case ( m_banded )
      do i = 1, n
        if ( i < a%r1(i) .or. i > a%r1(i) + a%r2(i) - a%r2(i-1) - 1 ) then
          ! No diagonal element.  Make room for one by densify-update-sparsify
          call allocate_test ( t, a%nRows, a%nCols, "T in UpdateDiagonal_0", &
            & ModuleName )
          call densify ( t, a )
          call updateDenseDiagonal ( t, x, s, i )
          call sparsify ( t, a, "T in UpdateDiagonal_0", ModuleName ) ! A := T
          call deallocate_test ( t, "T in UpdateDiagonal_0", ModuleName )
          return
        end if
        if ( myInvert ) then
          if ( x(i) == 0.0_r8 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
          v = s / x(i)
        else
          v = s * x(i)
        end if
        a%values(a%r2(i-1)+i-a%r1(i)+1,1) = &
          & a%values(a%r2(i-1)+i-a%r1(i)+1,1) + v
      end do
    case ( m_column_sparse )
      do i = 1, n
        do j = a%r1(i-1)+1, a%r1(i) ! hunt for the diagonal subscript
          if ( a%r2(j) == i ) then
            if ( myInvert ) then
              if ( x(i) == 0.0_r8 ) call MLSMessage ( MLSMSG_Error, moduleName, &
                & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
              v = s / x(i)
            else
              v = s * x(i)
            end if
            a%values(j,1) = a%values(j,1) + v
          else if ( a%r2(j) > i ) then
            ! No diagonal element.  Make room for one by densify-update-sparsify
            call allocate_test ( t, a%nRows, a%nCols, "T in UpdateDiagonal_0", &
              & ModuleName )
            call densify ( t, a )
            call updateDenseDiagonal ( t, x, s, i )
            call sparsify ( t, a, "T in UpdateDiagonal_0", moduleName ) ! A := T
            call deallocate_test ( t, "T in UpdateDiagonal_0", moduleName )
            return
          end if
        end do
      end do
    case ( m_full )
      call updateDenseDiagonal ( a%values, x, s, 1 )
    end select

  contains
    subroutine UpdateDenseDiagonal ( T, X, S, START )
      real(r8), intent(inout) :: T(:,:) ! Matrix element to update
      real(r8), intent(in) :: X(:)      ! Vector update diagonal
      real(r8), intent(in) :: S         ! Sign for X, +1 or -1.
      integer, intent(in) :: START      ! Where to start
      integer :: I
      do i = start, n
        if ( myInvert ) then
          if ( x(i) == 0.0_r8 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
          v = s / x(i)
        else
          v = s * x(i)
        end if
        t(i,i) = t(i,i) + v
      end do
    end subroutine UpdateDenseDiagonal
  end subroutine UpdateDiagonalVec_0

! =====     Private Procedures     =====================================

  ! -------------------------------------------  CreateEmptyBlock  -----
  subroutine CreateEmptyBlock ( EmptyBlock )
    type(MatrixElement_T), intent(out) :: EmptyBlock
    ! Default initialization does the rest, because of intent(out)
  end subroutine CreateEmptyBlock

  ! ------------------------------------------  DUMP_MATRIX_BLOCK  -----
  subroutine DUMP_MATRIX_BLOCK ( MATRIX_BLOCK, NAME, DETAILS, BOUNDS )
    type(MatrixElement_T), intent(in) :: MATRIX_BLOCK
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: DETAILS   ! Print details, default true
    integer, intent(in), optional :: BOUNDS(4) ! Dump only Bounds(1):Bounds(2)
                                               !        X  Bounds(3):Bounds(4)
    integer :: MY_DETAILS
    my_details = 1
    if ( present(details) ) my_details = details
    if ( present(name) ) call output ( name, advance='yes' )
    call output ( '  ' )
    call output ( matrix_block%nrows ); call output ( " Rows, " )
    call output ( matrix_block%ncols ); call output ( " Columns, " )
    select case ( matrix_block%kind )
    case ( m_banded )
      call output ( 'Banded' )
      if ( my_details > 0 ) then
        call output ( ' First-nonzero-rows =', advance='yes' )
        call dump ( matrix_block%r1(1:) )
        call output ( '  Last-value-in-column =', advance='yes' )
        call dump ( matrix_block%r2(1:) )
      end if
    case ( m_column_sparse )
      call output ( 'Column-sparse' )
      if ( my_details > 0 ) then
        call output ( ' Last-in-column =', advance='yes' )
        call dump ( matrix_block%r1(1:) )
        call output ( '  Rows =', advance='yes' )
        call dump ( matrix_block%r2(1:) )
      end if
    case ( m_full )
      call output ( 'Full' )
    end select
    if ( matrix_block%kind == M_Absent ) then
      call output ( 'Absent', advance='yes' )
    else if ( my_details > 1 ) then
      call output ( ' Values =', advance='yes' )
      if ( present(bounds) ) then
        call dump ( matrix_block%values(bounds(1):bounds(2),bounds(3):bounds(4)) )
      else
        call dump ( matrix_block%values )
      end if
    else
      call output ( ' ' )
      call output ( size(matrix_block%values) )
      call output ( ' values.', advance='yes' )
    end if
  end subroutine DUMP_MATRIX_BLOCK

end module MatrixModule_0

! $Log$
! Revision 2.45  2001/09/24 23:01:11  vsnyder
! Make consistent/correct lower bound calculation for MASK array
!
! Revision 2.44  2001/07/19 17:54:59  vsnyder
! Correct blunders in banded matrix*matrix and matrix*vector.
! Handle "subtract" argument differently.
!
! Revision 2.43  2001/07/16 20:36:49  livesey
! Added fix for empty columns in banded MatrixVector multiply
!
! Revision 2.42  2001/07/11 22:07:57  vsnyder
! Interim commit -- may still be broken
!
! Revision 2.41  2001/06/28 01:05:59  vsnyder
! Allow last diagonal element in Cholesky factor to be tiny
!
! Revision 2.40  2001/06/27 01:15:10  vsnyder
! XMASK and YMASK arguments of MultiplyMatrixBlocks need to be pointers
! because they might not be associated in the caller.  Therefore they
! cannot have specified lower bounds.  Therefore we need to use i/b+1
! (j/b+1) for subscripts of xm (ym).
!
! Revision 2.39  2001/06/26 23:56:04  vsnyder
! Make [XY]MASK 'target' instead of 'pointer' so they can have lower bound
!
! Revision 2.38  2001/06/26 20:40:33  vsnyder
! Simplify by using zero for lower bound for first dimension of mask
!
! Revision 2.37  2001/06/04 22:41:37  livesey
! Various bug fixes associated with m_banded.  Some still remain to be
! solved though
!
! Revision 2.36  2001/06/01 01:03:39  vsnyder
! Add 'sqrt' option to 'GetDiagonal_0'; add 'Multiply' generic
!
! Revision 2.35  2001/05/30 21:53:16  vsnyder
! Finish? 'invert' argument in 'UpdateDiagonalVec_0'
!
! Revision 2.34  2001/05/30 20:18:01  vsnyder
! Add 'invert' argument to 'UpdateDiagonal'
!
! Revision 2.33  2001/05/24 23:15:24  vsnyder
! Don't scale absent blocks -- their VALUES pointer isn't associated
!
! Revision 2.32  2001/05/24 18:13:28  vsnyder
! Make DenseCholesky public instead of internal; cosmetic changes
!
! Revision 2.31  2001/05/22 19:09:13  vsnyder
! Implement Col_L1
!
! Revision 2.30  2001/05/19 00:13:43  vsnyder
! Correct SolveCholesky*_0
!
! Revision 2.29  2001/05/17 20:17:56  vsnyder
! Implement GetMatrixElement.  Change handling of mask in MultiplyMatrixBlocks.
!
! Revision 2.28  2001/05/12 01:05:23  vsnyder
! Change 'details' argumet of 'dump_matrix_block' to integer
!
! Revision 2.27  2001/05/11 22:02:06  vsnyder
! Correct errors in SolveCholeskyV_0
!
! Revision 2.26  2001/05/10 22:53:36  vsnyder
! Handle empty block creation differently.  Add Update and Subtact to
! MultiplyMatrixVector* where it wasn't before.  Get CholeskyFactor_1 to work.
!
! Revision 2.25  2001/05/10 02:14:11  vsnyder
! Repair CloneBlock, MaxAbsVal, MultiplyMatrixBlocks
!
! Revision 2.24  2001/05/09 19:45:37  vsnyder
! More work correcting blunders in sparse matrix code.  Add BandHeight
! argument to CreateBlock.  Correct loss of lower bounds in CloneBlock.
!
! Revision 2.23  2001/05/09 01:58:12  vsnyder
! Improper intent(out) -> intent(inout), don't access an absent optional dummy
!
! Revision 2.22  2001/05/08 20:29:40  vsnyder
! Periodic commit -- workong on sparse matrix blunders
!
! Revision 2.21  2001/05/03 02:10:26  vsnyder
! Nullify a bunch of pointers that should have been but weren't.  Use a
! disassociated VALUES array instead of a zero-size one for absent blocks.
!
! Revision 2.20  2001/04/30 23:44:25  vsnyder
! Correct/remove some incorrect size tests in MultiplyMatrixVectorNoT
!
! Revision 2.19  2001/04/30 17:47:18  livesey
! Reverted to original size of R2, the problem must be somewhere else.
!
! Revision 2.18  2001/04/28 07:03:21  livesey
! Removed a print statement
!
! Revision 2.17  2001/04/28 05:04:16  livesey
! Temporarily changed dot to dot_product in MultiplyMatrixVectorNoT, to
! avoid run time error I don't understand.
!
! Revision 2.16  2001/04/28 04:40:17  livesey
! Some tidying up, removing unnecessary(?) tests for square matrices
! in multiplyMatrixVector and its relatives.  Also changing allocation
! of r2 for m_banded, as it needs an extra element.
!
! Revision 2.15  2001/04/28 01:33:02  livesey
! Now DestroyBlock doesn't destroy absent blocks as they all point to the same place.
!
! Revision 2.14  2001/04/25 00:50:09  vsnyder
! Make MultiplyMatrixNoT generic
!
! Revision 2.13  2001/04/11 22:43:54  vsnyder
! Fold Deallocate_test into sparsify
!
! Revision 2.12  2001/02/22 01:55:06  vsnyder
! Add code to invert a Cholesky factor
!
! Revision 2.11  2001/02/09 18:37:16  pwagner
! Commented-out statements that offended NAG v4.0
!
! Revision 2.10  2001/01/26 19:00:01  vsnyder
! Periodic commit
!
! Revision 2.9  2001/01/19 23:50:45  vsnyder
! Periodic commit
!
! Revision 2.8  2000/11/23 01:09:19  vsnyder
! Add provision to ignore specified columns during matrix-matrix multiply
!
! Revision 2.7  2000/11/15 00:18:26  vsnyder
! Added assignment(=) interface, row scale, column scale
!
! Revision 2.6  2000/11/10 00:28:13  vsnyder
! Added multiply untransposed matrix * vector
!
! Revision 2.5  2000/11/09 01:22:43  vsnyder
! Periodic commit -- still under construction
!
! Revision 2.4  2000/10/13 22:22:48  vsnyder
! Change name of multiply operator from .XT. to .TX.
!
! Revision 2.3  2000/10/12 20:10:08  vsnyder
! Make default accessibility private
!
! Revision 2.2  2000/10/10 23:11:34  vsnyder
! Correct number of rows and columns for zero matrix.
!
! Revision 2.1  2000/10/04 20:24:45  vsnyder
! Initial entry
!
