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
  public :: Add_Matrix_Blocks, Assignment(=), CholeskyFactor, &
    & CholeskyFactor_0, ClearRows, ClearRows_0, CloneBlock, ColumnScale, &
    & ColumnScale_0, CopyBlock, CreateBlock, CreateBlock_0, Densify, &
    & DestroyBlock, Dump, M_Absent, M_Banded, M_Column_Sparse, M_Full, &
    & MatrixElement_T, MultiplyMatrixBlocks, MultiplyMatrixVector, &
    & MultiplyMatrixVectorNoT, MultiplyMatrixVector_0, operator(+), &
    & operator(.TX.), RowScale, RowScale_0, SolveCholesky, SolveCholeskyM_0, &
    & SolveCholeskyV_0, Sparsify, UpdateDiagonal, UpdateDiagonal_0

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

  interface DUMP
    module procedure DUMP_MATRIX_BLOCK
  end interface

  interface MultiplyMatrixVector
    module procedure MultiplyMatrixVector_0
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
    module procedure UpdateDiagonal_0
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
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
    integer :: NROWS, NCOLS                  ! Numbers of rows and columns
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
  type(MatrixElement_T), save, private :: EmptyBlock     ! To fill the blocks.
    ! It will have zero-size arrays instead of nullified ones, so that we
    ! can copy the components without looking at the block kind if one
    ! operand is absent during an operation.
  logical, save, private :: FIRST = .true.        ! to create EmptyBlock
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
        call allocate_test ( zb%r2, size(x%r2), "zb%r2", ModuleName )
        zb%r2(0) = 0
        do k = 1, size(x%r1)                     ! Calculate size of Values
          zb%r2(k) = zb%r2(k-1) + &
            & max( x%r2(k)-x%r2(k-1)+x%r1(k), y%r2(k)-y%r2(k-1)+y%r1(k) ) - &
            & zb%r1(k) + 1
        end do
        call allocate_test ( zb%values, zb%r2(size(x%r1)), 1, "zb%values", &
          & ModuleName )
        zb%values = 0.0_r8 ! ??? Improve this in level 1.0 by only filling
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
        call sparsify ( z, zb )                  ! Zb = Z
        call deallocate_test ( z, "Z in Add_Matrix_Blocks", ModuleName )
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
        call sparsify ( z, zb )                  ! Zb = Z
        call deallocate_test ( z, "Z in Add_Matrix_Blocks", ModuleName )
      case ( M_Full )                            ! X col sparse, Y full
        call CopyBlock ( zb, y )                 ! Zb = y
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
  ! components.  Notice that CopyBlock does a deep copy.  If one has
  ! Z = X in a loop, it is therefore necessary only to destroy Z after
  ! the loop.
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
  ! triangular. Otherwise, replace Z such that Z(output) = Z(input)^T
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
    real(r8), pointer, dimension(:,:) :: XIN  ! A pointer to the input,
    !                           data, or a densified copy of it
    type(MatrixElement_T), pointer :: X       ! XOPT or Z, depending on whether
    !                           XOPT is present or absent, respectively.
    real(r8), pointer, dimension(:,:) :: ZT   ! A local full result that is
    !                           sparsified at the end.

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
    case ( M_Column_Sparse )
      call allocate_test ( zt, nc, nc, "ZT for CholeskyFactor", ModuleName )
      call allocate_test ( xin, nc, nc, "XIN for CholeskyFactor", ModuleName )
      ! ??? Densify and then compute Cholesky decomposition of dense block.
      ! ??? If necessary, improve this in level 1.0 by working directly
      ! ??? with sparse input.
      call densify ( xin, x )
      call denseCholesky ( zt, xin )
      call sparsify ( zt, z )
      call deallocate_test ( zt, "ZT in CholeskyFactor", ModuleName )
      call deallocate_test ( xin, "XIN in CholeskyFactor", ModuleName )
    case ( M_Banded )
      call allocate_test ( zt, nc, nc, "ZT in CholeskyFactor", &
        & ModuleName )
      call allocate_test ( r1, nc, "R1 in CholeskyFactor", ModuleName )
      do i = 1, nc
        r1(i) = i             ! We know the diagonal will get a value
      end do
      do i = 1, nc
        zt(i,1:i-1) = 0.0_r8  ! Clear left from the diagonal (helps Sparsify!)
        ii = i - x%r1(i)      ! Offset in VALUES of (I,I) element
        if ( ii < 0 .or. ii > x%r2(i) - x%r2(i-1) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        end if
!       g = x%values(ii+x%r2(i-1)+1,1) - &
!           & dot_product(zt(r1(i):i-1,i),zt(r1(i):i-1,i))
        g = x%values(ii+x%r2(i-1)+1,1) - &
            & dot( i-r1(i), zt(r1(i),i), 1, zt(r1(i),i), 1 )
        if ( g <= sqrt(tiny(0.0_r8)) ) then
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
    case ( M_Full )
      if ( present(xopt) ) then
        call createBlock ( z, nc, nc, M_Full )
      end if
      call denseCholesky ( z%values, x%values )
    end select

  contains
    subroutine DenseCholesky ( zt, xin )
    ! Do the Cholesky decomposition of XIN giving ZT.
      real(r8) :: ZT(:,:), XIN(:,:)
      real(r8), save :: TOL = -1.0_r8
      if ( tol < 0.0_r8 ) tol = sqrt(tiny(0.0_r8))
      do i = 1, nc
        zt(i+1:nc,i) = 0.0_r8 ! Clear below the diagonal (helps Sparsify!)
!       g = xin(i,i) - dot_product(zt(1:i-1,i),zt(1:i-1,i))
        g = xin(i,i) - dot( i-1, zt(1,i), 1, zt(1,i), 1 )
        if ( g <= tol ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        d = sqrt(g)
        zt(i,i) = d
        do j = i+1, nc
!         zt(i,j) = ( xin(i,j) - dot_product(zt(1:i-1,i),zt(1:i-1,j)) ) / d
          zt(i,j) = ( xin(i,j) - dot( i-1, zt(1,i), 1, zt(1,j), 1 ) ) / d
        end do ! j
      end do ! i
    end subroutine DenseCholesky
  end subroutine CholeskyFactor_0

  ! ------------------------------------------------  ClearRows_0  -----
  subroutine ClearRows_0 ( X, MASK )
  ! Clear the rows of X for which MASK has nonzero bits.
    type(MatrixElement_T), intent(inout) :: X
    integer, dimension(0:), intent(in) :: MASK
    integer :: I, J                ! Subscripts and row indices
    integer :: NB, NW              ! Indices of bit and word in MASK
    select case ( x%kind )
    case ( M_Absent )
    case ( M_Banded )              ! ??? Adjust the sparsity representation ???
      do j = 1, x%ncols
        do i = x%r1(j), x%r1(j) + x%r2(j) - x%r2(j-1) - 1 ! row numbers
          if ( btest( mask(i/b), mod(i,b) ) ) &
            & x%values(x%r2(j-1) + i - x%r1(j) + 1, 1) = 0.0_r8
        end do ! i
      end do ! j = 1, x%ncols
    case ( M_Column_Sparse )       ! ??? Adjust the sparsity representation ???
      do j = 1, x%ncols
        do i = x%r1(j-1)+1, x%r1(j)
          if ( btest( mask(x%r2(i)/b), mod(x%r2(i),b) ) ) &
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
    if ( first ) call CreateEmptyBlock
    call destroyBlock ( z )
    if ( x%kind == M_absent ) then
      z = emptyBlock
      return
    end if
    z%kind = x%kind
    z%nRows = x%nRows; z%nCols = x%nCols
    call allocate_test ( z%r1, size(x%r1), "z%r1", ModuleName )
    z%r1 = x%r1
    call allocate_test ( z%r2, size(x%r2), "z%r2", ModuleName )
    z%r2 = x%r2
    call allocate_test ( z%values, size(x%values,1), size(x%values,2), &
      & "z%values", ModuleName )
  end subroutine CloneBlock

  ! -------------------------------------------------  ColumnScale_0  -----
  subroutine ColumnScale_0 ( X, V, NEWX ) ! Z = X V where V is a diagonal
  !                                     matrix represented by a vector and
  !                                     Z is either X or NEWX.
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

  ! --------------------------------------------------  CopyBlock  -----
  subroutine CopyBlock ( Z, X ) ! Destroy Z, deep Z = X, including the values
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    call CloneBlock ( Z, X )
    z%values = x%values
  end subroutine CopyBlock

  ! ----------------------------------------------  CreateBlock_0  -----
  subroutine CreateBlock_0 ( Z, nRows, nCols, Kind, NumberNonzero )
  ! Create a matrix block, but don't fill any elements or structural
  ! information.  The "NumberNonzero" is required if and only if the
  ! "Kind" argument has the value M_Banded or M_Column_Sparse.
  ! The block is first destroyed, so as not to have a memory leak.

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

    type(MatrixElement_T), intent(out) :: Z
    integer, intent(in) :: nRows, nCols, Kind
    integer, intent(in), optional :: NumberNonzero ! Only for M_Banded and
                                                   ! M_Column_Sparse

    call destroyBlock ( z )
    select case ( kind )
    case ( M_Absent )
      if ( first ) call CreateEmptyBlock
      z = emptyBlock
    case ( M_Banded )
      call allocate_test ( z%r1, nCols, "z%r1", ModuleName )
      call allocate_test ( z%r2, nCols, "z%r2", ModuleName, lowBound=0 )
      z%r2(0) = 0
      call allocate_test ( z%values, NumberNonzero, 1, "z%values", ModuleName )
    case ( M_Column_sparse )
      call allocate_test ( z%r1, nCols, "z%r1", ModuleName, lowBound=0 )
      z%r1(0) = 0
      call allocate_test ( z%r2, numberNonzero, "z%r2", ModuleName )
      call allocate_test ( z%values, NumberNonzero, 1, "z%values", ModuleName )
    case ( M_Full )
      call allocate_test ( z%r1, 0, "z%r1", ModuleName )
      call allocate_test ( z%r2, 0, "z%r1", ModuleName )
      call allocate_test ( z%values, nRows, nCols, "z%values", ModuleName )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid matrix block kind in CreateBlock" )
    end select
    z%nRows = nRows
    z%nCols = nCols
    z%kind = kind
  end subroutine CreateBlock_0

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

  ! -----------------------------------------------  DestroyBlock  -----
  subroutine DestroyBlock ( B )
  ! Deallocate the pointer components of the matrix block B
    type(MatrixElement_T), intent(inout) :: B
    call deallocate_test ( b%r1, "b%r1", ModuleName )
    call deallocate_test ( b%r2, "b%r2", ModuleName )
    call deallocate_test ( b%values, "b%values", ModuleName )
    b%kind = M_Absent
  end subroutine DestroyBlock

  ! ---------------------------------------  MultiplyMatrixBlocks  -----
  subroutine MultiplyMatrixBlocks ( XB, YB, ZB, UPDATE, XMASK, YMASK )
  ! ZB = XB^T YB if UPDATE is absent or false;
  ! ZB = ZB + XB^T YB if UPDATE is present and true.
  ! If XMASK (resp. YMASK) is present, ignore columns of XB (resp. YB)
  ! that correspond to nonzero bits of XMASK (resp. YMASK).
    type(MatrixElement_T), intent(in) :: XB, YB
    type(MatrixElement_T), intent(inout) :: ZB
    logical, intent(in), optional :: UPDATE
    integer, intent(in), optional, dimension(0:) :: XMASK, YMASK

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! If UPDATE is absent or false, it is important to invoke DestroyBlock
  ! using the result of this function after it is no longer needed.
  ! Otherwise, a memory leak will result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K, L, M, N, P, R
    logical :: MY_UPD
    real(r8), pointer, dimension(:,:) :: Z   ! Temp for sparse * sparse

    my_upd = .false.
    if ( present(update) ) my_upd = update
    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      if ( my_upd) call createBlock ( zb, xb%nCols, yb%nCols, M_Absent )
      return
    end if
    if ( xb%nrows /= yb%nrows ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "XB and YB Matrix sizes incompatible in Multiply_Matrix_Blocks" )
    if ( my_upd ) then
      zb%nrows = xb%ncols
      zb%ncols = yb%ncols
    else
      if ( xb%ncols /= zb%nrows .or. yb%ncols /= zb%ncols ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "ZB Matrix size incompatible in Multiply_Matrix_Blocks" )
    end if
    select case ( xb%kind )
    case ( M_Banded )
      select case ( yb%kind )
      case ( M_Banded )       ! XB banded, YB banded
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, zb%nrows, zb%ncols, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
        if ( update ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, zb%ncols    ! Columns of Z = columns of YB
          if ( present(ymask) ) then
            if ( btest(ymask(j/b),mod(j,b)) ) cycle
          end if
          p = yb%r2(j-1)+1
          k = yb%r1(j)        ! k,l = row indices of YB
          l = k+yb%r2(j)-p
          p = p-k
          do i = 1, zb%nrows  ! Rows of Z = columns of XB
            ! Inner product of column I of XB with column J of YB
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
            end if
            m = xb%r1(i)
            r = xb%r2(i-1)+1-m
            m = max(k,m)      ! m,n = row indices of intersection of XB and YB
            n = min(l,xb%r2(i)-r)
!           z(i,j) = z(i,j) + &
!                    & dot_product ( xb%values(r+m:r+n,1), yb%values(p+m:p+n,1) )
            z(i,j) = z(i,j) + &
                     & dot( n-m+1, xb%values(r+m,1), 1, yb%values(p+m,1), 1 )
          end do ! i
        end do ! j
        call sparsify ( z, zb )
        call deallocate_test ( z, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Column_sparse ) ! XB banded, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded X sparse in Multiply_Matrix_Blocks", ModuleName )
        if ( update ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, yb%ncols    ! Columns of Z
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          do i = 1, xb%ncols  ! Rows of Z
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
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
                z(i,j) = z(i,j) + xb%values(l,1) * yb%values(n,1)
                l = l + 1
                k = k + 1
                n = n + 1
                if ( n > yb%r1(j) ) exit
                m = yb%r2(n)
              end if
            end do
          end do ! i
        end do ! j
        call sparsify ( z, zb )
        call deallocate_test ( z, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Full )         ! XB banded, YB full
        if ( .not. my_upd ) then
          call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
          zb%values = 0.0_r8
        end if
        do i = 1, xb%ncols    ! Rows of ZB
          if ( present(xmask) ) then
            if ( btest(xmask(i/b),mod(i,b)) ) cycle
          end if
          m = xb%r1(i)        ! Index of first row of XB with nonzero value
          k = xb%r2(i-1) + 1
          l = xb%r2(i)
          do j = 1, yb%ncols  ! Columns of ZB
            if ( present(ymask) ) then
              if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
!           zb%values(i,j) = zb%values(i,j) + dot_product( &
!             & xb%values(k:l,1), yb%values(m:m+l-k,j) )
            zb%values(i,j) = zb%values(i,j) + dot( l-k+1, &
              & xb%values(k,1), 1, yb%values(m,j), 1 )
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
        if ( update ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, yb%ncols    ! Columns of Z
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          do i = 1, xb%ncols  ! Rows of Z
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
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
                z(i,j) = z(i,j) + xb%values(l,1) * yb%values(n,1)
                l = l + 1
                if ( l > xb%r1(i) ) exit
                k = xb%r2(l)
                n = n + 1
                m = m + 1
              end if
            end do
          end do ! i
        end do ! j
        call sparsify ( z, zb )
        call deallocate_test ( z, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Column_sparse ) ! XB column-sparse, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for sparse X sparse in Multiply_Matrix_Blocks", ModuleName )
        if ( update ) then
          call densify ( z, zb )
        else
          z = 0.0_r8
        end if
        do j = 1, yb%ncols    ! Columns of Z
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          do i = 1, xb%ncols  ! Rows of Z
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
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
                z(i,j) = z(i,j) + xb%values(l,1) * yb%values(n,1)
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
        call sparsify ( z, zb )
        call deallocate_test ( z, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Full )         ! XB column-sparse, YB full
        if ( .not. my_upd ) then
          call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
          zb%values = 0.0_r8
        end if
        do j = 1, yb%ncols    ! Columns of ZB
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          do i = 1, xb%ncols  ! Rows of ZB
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
            end if
            k = xb%r1(i-1)+1
            l = xb%r1(i)
            ! Inner product of column I of XB with column J of YB
!           zb%values(i,j) = zb%values(i,j) + dot_product( &
!             & xb%values(k:l,1), yb%values(xb%r2(k:l),j) )
            zb%values(i,j) = zb%values(i,j) + dot( l-k+1, &
              & xb%values(k,1), 1, yb%values(xb%r2(k),j), 1 )
          end do ! i
        end do ! j
      end select
    case ( M_Full )
      if ( .not. my_upd ) then
        call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
        zb%values = 0.0_r8
      end if
      select case ( yb%kind )
      case ( M_Banded )       ! XB full, YB banded
        do j = 1, yb%ncols    ! Columns of ZB
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          m = yb%r1(j)        ! Index of first row of YB with nonzero value
          k = yb%r2(j-1)+1    ! K and L are indices of YB
          l = yb%r2(j)
          do i = 1, xb%ncols  ! Rows of ZB
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
!           zb%values(i,j) = zb%values(i,j) + dot_product( &
!             & xb%values(m:m+l-k,i), yb%values(k:l,1))
            zb%values(i,j) = zb%values(i,j) + dot( l-k+1, &
              & xb%values(m,i), 1, yb%values(k,1), 1 )
          end do ! i
        end do ! j
      case ( M_Column_sparse ) ! XB full, YB column-sparse
        do j = 1, yb%ncols    ! Columns of ZB
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
          k = yb%r1(j-1)+1    ! K and L are indices of YB
          l = yb%r1(j)
          do i = 1, xb%ncols  ! Rows of ZB
            if ( present(xmask) ) then
              if ( btest(xmask(i/b),mod(i,b)) ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
            zb%values(i,j) = zb%values(i,j) + dot_product( &
              & xb%values(yb%r2(k:l),i), yb%values(k:l,1) )
          end do ! i
        end do ! j
      case ( M_Full )         ! XB full, YB full
        if ( present(xmask) .or. present(ymask) ) then
          do j = 1, yb%ncols    ! Columns of ZB
          if ( present(ymask) ) then
            if ( btest(ymask(j/b+1),mod(j,b)) ) cycle
          end if
            do i = 1, xb%ncols  ! Rows of ZB
              if ( present(xmask) ) then
                if ( btest(xmask(i/b),mod(i,b)) ) cycle
              end if
              zb%values(i,j) = zb%values(i,j) + dot( xb%nrows, xb%values(1,i), &
                               & 1, yb%values(1,j), 1 )
            end do ! i = 1, xb%ncols
          end do ! j = 1, yb%ncols
        else
          zb%values = zb%values + matmul(transpose(xb%values),yb%values)
        end if
      end select
    end select
  end subroutine MultiplyMatrixBlocks

  ! -------------------------------------  MultiplyMatrixVector_0  -----
  subroutine MultiplyMatrixVector_0 ( B, V, P, UPDATE )
  ! P = B^T V if UPDATE is absent or false.
  ! P = P + B^T V if UPDATE is present and true.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE

    integer :: I, M, N             ! Subscripts and loop inductors
    logical :: My_update
    integer :: V1                  ! Subscripts and loop inductors

    if ( b%nrows /= size(v) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and vector not compatible in MultiplyMatrixVector_0" )
    if ( b%ncols /= size(p) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and result not compatible in MultiplyMatrixVector_0" )
    if ( any(shape(v) /= shape(p)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Vectors not compatible in MultiplyMatrixVector_0" )
    my_update = .false.
    if ( present(update) ) my_update = update
    if ( .not. my_update ) p = 0.0_r8
    select case ( b%kind )
    case ( M_Absent )
    case ( M_Banded )
      do i = 1, size(p)
        v1 = b%r2(i-1)             ! starting position in B%VALUES - 1
        n = b%r2(i) - v1           ! how many values
        m = b%r1(i)                ! starting position in V
        p(i) = p(i) + dot(n, b%values(v1+1,1), 1, v(m), 1)
      end do ! i
     case ( M_Column_Sparse )
      do i = 1, size(p)
        do n = b%r1(i-1)+1, b%r1(i)
          p(i) = p(i) + b%values(n,1) * v(b%r2(n))
        end do ! n
      end do ! i
    case ( M_Full )
      do i = 1, size(p)
        p(i) = p(i) + dot(size(v), b%values(1,i), 1, v(1), 1)
      end do ! i
    end select
  end subroutine MultiplyMatrixVector_0

  ! ------------------------------------  MultiplyMatrixVectorNoT  -----
  subroutine MultiplyMatrixVectorNoT ( B, V, P, UPDATE, DoDiag )
  ! P = B V if UPDATE is absent or false.
  ! P = P + B V if UPDATE is present and true.
  ! Don't multiply by the diagonal element if doDiag (default true) is
  ! present and false.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE, DoDiag

    integer :: I, J, M, N          ! Subscripts and loop inductors
    logical :: My_diag, My_update
    integer :: V1                  ! Subscripts and loop inductors

    if ( b%nrows /= size(v) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and vector not compatible in MultiplyMatrixVector_0" )
    if ( b%ncols /= size(p) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Matrix block and result not compatible in MultiplyMatrixVector_0" )
    if ( any(shape(v) /= shape(p)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Vectors not compatible in MultiplyMatrixVector_0" )
    my_update = .false.
    if ( present(update) ) my_update = update
    my_diag = .true.
    if ( present(doDiag) ) my_diag = doDiag
    if ( .not. my_update ) p = 0.0_r8
    select case ( b%kind )
    case ( M_Absent )
    case ( M_Banded )
      do j = 1, size(p)            ! columns
        v1 = b%r2(j-1)             ! (starting position in B%VALUES) - 1
        n = b%r2(j) - v1           ! how many values
        m = b%r1(j)                ! starting row subscript in B%VALUES
        if ( my_diag ) then        ! do the whome matrix
          do i = m, m+n-1          ! rows
            p(i) = p(i) + b%values(v1+i-m+1,1) * v(j)
          end do ! i = 1, n
        else                       ! skip the diagonal
          do i = m, min(m+n-1,j-1) ! rows
            p(i) = p(i) + b%values(v1+i-m+1,1) * v(j)
          end do ! i = m, min(m+n-1,j-1)
          do i = max(m,j+1), m+n-1 ! rows
            p(i) = p(i) + b%values(v1+i-m+1,1) * v(j)
          end do ! i = m, min(m,n-1,j-1)
        end if
      end do ! j
     case ( M_Column_Sparse )
      do j = 1, size(p)            ! columns
        do n = b%r1(j-1)+1, b%r1(j)! rows
          i = b%r2(n)              ! row number
          if ( i/=j .or. my_diag ) &
            & p(i) = p(i) + b%values(n,1) * v(j)
        end do ! n
      end do ! j
    case ( M_Full )
      if ( my_diag ) then          ! do the whole matrix
        do i = 1, size(p)
          p(i) = p(i) + dot(size(v), b%values(i,1), size(b%values,1), v(1), 1)
        end do ! i
      else                         ! skip the diagonal
        do i = 1, size(p)
          p(i) = p(i) + dot(i-1, b%values(i,1), size(b%values,1), v(1), 1)
          p(i) = p(i) + &
            & dot(size(v)-i, b%values(i,i+1), size(b%values,1), v(i+1), 1)
        end do ! i
      end if
    end select
  end subroutine MultiplyMatrixVectorNoT

  ! -------------------------------------  NewMultiplyMatrixBlocks  ----
  function NewMultiplyMatrixBlocks ( X, Y ) result ( Z ) ! Z = X^T Y
    type(MatrixElement_T), intent(in) :: X, Y
    type(MatrixElement_T) :: Z

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    call MultiplyMatrixBlocks ( x, y, z, .false. )
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
    real(r8), pointer, dimension(:,:) :: XS  ! The solution, dense
    real(r8), pointer, dimension(:,:) :: UD  ! U, densified

    n = u%nrows
    if ( n /= u%nCols ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "U matrix in SolveCholeskyM_0 must be square" )
    my_b => x
    if ( present(b) ) my_b => b
    if ( n /= my_b%nrows ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "B matrix not compatible with U matrix in SolveCholeskyM_0" )
    my_t = .false.
    if ( present(transpose) ) my_t = transpose
    nc = b%nCols

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
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
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
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
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
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
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
      do i = n, 1, -1
        d = ud(i,i)
        if ( abs(d) < tiny(0.0_r8) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
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
    call sparsify ( xs, x )
    call deallocate_test ( xs, "XS in SolveCholeskyM_0", ModuleName )
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
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          ! dot_product( u%values(u%r2(i-1)+1:u%r2(i)-1,1), &
          ! &            b(u%r1(i):u%r1(i)+u%r2(i)-u%r2(i-1)-1) )
          my_b(i) = my_b(i) - dot(u%r2(i)-u%r2(i-1), &
            &                     u%values(u%r2(i-1)+1,1), 1, my_b(u%r1(i)), 1)
        end do ! i = 1, n
      case ( M_Column_Sparse )
        do i = 1, n
          if ( u%r2(u%r1(i)) /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is not triangular" )
          d = u%values(u%r1(i),1)
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          do h = u%r1(i-1)+1, u%r1(i)-1
            my_b(i) = my_b(i) - u%values(h,1) * my_b(u%r2(h))
          end do ! h = u%r1(i-1)+1, u%r1(i)-1
        end do ! i = 1, n
      case ( M_Full )
        do i = 1, n
          d = u%values(i,i)
          if ( abs(d) < tiny(0.0_r8) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyV_0 is singular" )
          ! dot_product( ud(1:i-1,i), my_b(1:i-1) )
          my_b(i) = my_b(i) - dot(i-1, ud(1,i), 1, my_b(1), 1)
        end do ! i = 1, n
      end select
    else             ! solve U X = B for X
      if ( u%kind == M_full ) then
        ud => u%values
      else
        call allocate_test ( ud, n, n, "UD in SolveCholeskyV_0", ModuleName )
        call densify ( ud, u )
      end if
      do i = n, 1, -1
        d = ud(i,i)
        if ( abs(d) < tiny(0.0_r8) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "U matrix in SolveCholeskyV_0 is singular" )
        ! dot_product( ud(i,i+1:n), my_b(i+1:n) )
        my_b(i) = my_b(i) - dot(n-i, ud(i,i+1), size(ud,1), my_b(i+1), 1)
      end do ! i = 1, n
      if ( u%kind /= M_Full ) &
        & call deallocate_test ( ud, "UD in SolveCholeskyV_0", ModuleName )
    end if ! my_t
  end subroutine SolveCholeskyV_0

  ! ---------------------------------------------------  Sparsify  -----
  subroutine Sparsify ( Z, B )
  ! Given an array Z, compute its sparse representation and store it
  ! in the matrix block B.
    real(r8), intent(in) :: Z(:,:)           ! Full array of values
    type(MatrixElement_T), intent(out) :: B  ! Z as a block, maybe sparse

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
  end subroutine Sparsify

  ! -------------------------------------------  UpdateDiagonal_0  -----
  subroutine UpdateDiagonal_0 ( A, LAMBDA )
  ! Add LAMBDA to the diagonal of A
    type(MatrixElement_T), intent(inout) :: A
    real(r8), intent(in) :: LAMBDA

    integer :: I, J                          ! Subscripts and loop inductors
    integer :: N                             ! min(a%ncols,a%nrows)
    real(r8), dimension(:,:), pointer :: T   ! A temporary dense matrix

    n = min(a%ncols,a%nrows)
    select case ( a%kind )
    case ( m_absent )
      call createBlock ( a, a%nRows, a%nCols, m_banded, n )
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
          call sparsify ( t, a )
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
            call sparsify ( t, a )
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

! =====     Private Procedures     =====================================

  ! -------------------------------------------  CreateEmptyBlock  -----
  subroutine CreateEmptyBlock
  ! Create the empty block if necessary
    if ( .not. first ) return
    first = .false.
    emptyBlock%kind = M_absent
    emptyBlock%nRows = 0; emptyBlock%nCols = 0 ! just so they're defined
    allocate ( emptyBlock%r1(0), emptyBlock%r2(0), emptyBlock%values(0,0) )
  end subroutine CreateEmptyBlock

  ! ------------------------------------------  DUMP_MATRIX_BLOCK  -----
  subroutine DUMP_MATRIX_BLOCK ( MATRIX_BLOCK, NAME, DETAILS, BOUNDS )
    type(MatrixElement_T), intent(in) :: MATRIX_BLOCK
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: DETAILS   ! Print details, default true
    integer, intent(in), optional :: BOUNDS(4) ! Dump only Bounds(1):Bounds(2)
                                               !        X  Bounds(3):Bounds(4)
    logical :: MY_DETAILS
    my_details = .true.
    if ( present(details) ) my_details = details
    if ( present(name) ) call output ( name, advance='yes' )
    call output ( matrix_block%nrows ); call output ( " Rows, " )
    call output ( matrix_block%ncols ); call output ( " Columns, " )
    select case ( matrix_block%kind )
    case ( m_banded )
      call output ( ' Banded' )
      if ( my_details ) then
        call output ( ' First-nonzero-rows =', advance='yes' )
        call dump ( matrix_block%r1(1:) )
        call output ( '  Last-value-in-column =', advance='yes' )
        call dump ( matrix_block%r2(1:) )
      end if
    case ( m_column_sparse )
      call output ( ' Column-sparse' )
      if ( my_details ) then
        call output ( ' Last-in-column =', advance='yes' )
        call dump ( matrix_block%r1(1:) )
        call output ( '  Rows =', advance='yes' )
        call dump ( matrix_block%r2(1:) )
      end if
    case ( m_full )
      call output ( 'Full' )
    end select
    if ( matrix_block%kind == M_Absent ) then
      call output ( ' Absent', advance='yes' )
    else if ( my_details ) then
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
