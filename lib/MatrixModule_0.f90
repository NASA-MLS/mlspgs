! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!OCL INDEPENDENT (dot) ! For LF95 auto-parallelization

!=============================================================================
module MatrixModule_0          ! Low-level Matrices in the MLS PGS suite
!=============================================================================

! This module provides the elementary matrix type.  Blocks of this
! type are used to compose block matrices.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DOT_M, only: DOT
  use DUMP_0, only: DUMP
  use Gemm_M, only: GEMM
  use Gemv_M, only: GEMV
  use MLSCommon, only: RM, RV, R4, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: OUTPUT
  use VectorsModule, only: M_LinAlg

  implicit NONE
  private
  public :: Add_Matrix_Blocks, CheckBlocks, CheckIntegrity, Assignment(=), CholeskyFactor
  public :: CholeskyFactor_0, ClearRows, ClearRows_0, CloneBlock, ColumnScale
  public :: Col_L1, CopyBlock, CreateBlock, CreateBlock_0
  public :: DenseCholesky, Densify, DestroyBlock, DestroyBlock_0, Dump
  public :: GetDiagonal, GetMatrixElement, GetMatrixElement_0
  public :: GetVectorFromColumn
  public :: InvertCholesky, InvertCholesky_0
  public :: InvertDenseCholesky, InvertDenseCholesky_0
  public :: M_Absent, M_Banded, M_Column_Sparse, M_Full, M_Unknown
  public :: MatrixInversion, MatrixInversion_0, MatrixElement_T
  public :: MaxAbsVal, MaxAbsVal_0, MinDiag, MinDiag_0
  public :: Multiply, MultiplyMatrix_XY, MultiplyMatrix_XY_0
  public :: MultiplyMatrix_XY_T, MultiplyMatrix_XY_T_0
  public :: MultiplyMatrix_XTY_0
  public :: MultiplyMatrixVector
  public :: MultiplyMatrixVectorNoT
  public :: NullifyMatrix, NullifyMatrix_0
  public :: operator(+), operator(.TX.), ReflectMatrix, RowScale
  public :: ScaleBlock, SolveCholesky, SolveCholeskyM_0
  public :: Sparsify, Spill, Spill_0, SubBlockLength
  public :: TransposeMatrix, UpdateDiagonal

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment(=)
    module procedure AssignBlock
  end interface

  interface CheckIntegrity
    module procedure CheckIntegrity_0
  end interface

  interface CholeskyFactor
    module procedure CholeskyFactor_0
  end interface

  interface ClearRows
    module procedure ClearRows_0
  end interface

  interface ColumnScale
    module procedure ColumnScale_0_r4, ColumnScale_0_r8
  end interface

  interface CreateBlock
    module procedure CreateBlock_0
  end interface

  interface Densify
    module procedure DensifyA, DensifyB
  end interface

  interface DestroyBlock
    module procedure DestroyBlock_0
  end interface

  interface DUMP
    module procedure DUMP_MATRIX_BLOCK
  end interface

  interface GetDiagonal
    module procedure GetDiagonal_0_r4, GetDiagonal_0_r8
  end interface

  interface GetMatrixElement
    module procedure GetMatrixElement_0
  end interface

  interface GetVectorFromColumn
    module procedure GetVectorFromColumn_0_r4, GetVectorFromColumn_0_r8
  end interface

  interface InvertCholesky
    module procedure InvertCholesky_0, InvertDenseCholesky_0
  end interface

  interface InvertDenseCholesky
    module procedure InvertDenseCholesky_0
  end interface

  interface MatrixInversion
    module procedure MatrixInversion_0
  end interface

  interface MaxAbsVal
    module procedure MaxAbsVal_0
  end interface

  interface MinDiag
    module procedure MinDiag_0
  end interface

  interface Multiply
    module procedure MultiplyMatrix_XTY_0, &
      & MultiplyMatrixVector_0_r4, MultiplyMatrixVector_0_r8
  end interface

  interface MultiplyMatrixVector ! A^T V
    module procedure MultiplyMatrixVector_0_r4, MultiplyMatrixVector_0_r8
  end interface

  interface MultiplyMatrixVectorNoT ! A V
    module procedure MultiplyMatrixVectorNoT_0_r4, MultiplyMatrixVectorNoT_0_r8
  end interface

  interface MultiplyMatrix_XY
    module procedure MultiplyMatrix_XY_0
  end interface

  interface MultiplyMatrix_XY_T
    module procedure MultiplyMatrix_XY_T_0
  end interface

  interface MultiplyMatrix_XTY
    module procedure MultiplyMatrix_XTY_0
  end interface

  interface NewMultiplyMatrixVector_0 ! (see .tx. below)
    module procedure NewMultiplyMatrixVector_0_r4, NewMultiplyMatrixVector_0_r8
  end interface

  interface NullifyMatrix
    module procedure NullifyMatrix_0
  end interface

  interface operator (+)
    module procedure Add_Matrix_Blocks_Unscaled
  end interface

  interface operator ( .TX. ) ! A^T * B
    module procedure NewMultiplyMatrix_XTY_0, &
      & NewMultiplyMatrixVector_0_r4, NewMultiplyMatrixVector_0_r8
  end interface

  interface ReflectMatrix
    module procedure ReflectMatrix_0
  end interface

  interface RowScale
    module procedure RowScale_0_r4, RowScale_0_r8
  end interface

  interface SolveCholesky
    module procedure SolveCholeskyM_0, SolveCholeskyV_0_r4, SolveCholeskyV_0_r8, &
      & SolveCholeskyA_0_r4, SolveCholeskyA_0_r8
  end interface

  interface Sparsify
    module procedure SparsifyA, SparsifyB
  end interface

  interface Spill
    module procedure Spill_0
  end interface

  interface TransposeMatrix
    module procedure TransposeMatrix_0
  end interface

  interface UpdateDiagonal
    module procedure UpdateDiagonal_0_r8, UpdateDiagonalVec_0_r8, &
      & UpdateDiagonal_0_r4, UpdateDiagonalVec_0_r4
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private not_used_here
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
  integer, parameter :: M_Unknown = 4        ! We don't know yet
    ! This is used for reading l2pc files where we don't want to load the block
    ! into memory until we know we need it.  The MatrixModule_0 code doesn't
    ! understand such blocks, so use them with care.

  type MatrixElement_T
    integer :: KIND = M_Absent               ! Kind of block -- one of the
      !                                        M_... parameters above
    integer :: nRows = 0, nCols = 0          ! Numbers of rows and columns
    integer, pointer, dimension(:) :: R1 => NULL()     ! Indexed by the column
      ! number. Used for the first column number if KIND = M_Banded, as
      ! described above for M_Column_sparse if KIND = M_Column_sparse, and not
      ! used otherwise.
    integer, pointer, dimension(:) :: R2 => NULL()     ! Indexed by the
      ! column number if KIND = M_Banded, by elements of R1 if KIND =
      ! M_Column_sparse, and not used otherwise.  See M_Banded and
      ! M_Column_sparse above.
    real(rm), pointer, dimension(:,:) :: VALUES => NULL()   ! Values of the
      ! matrix elements.  Indexed by row and column indices if KIND == M_Full,
      ! by elements in the range of values of R1 if KIND == M_Banded, and by
      ! elements of R2 if KIND == M_Column_sparse.
  end type MatrixElement_T

  ! - - -  Private data     - - - - - - - - - - - - - - - - - - - - - -
  real, parameter, private :: COL_SPARSITY = 0.5  ! If more than this
    ! fraction of the elements between the first and last nonzero in a
    ! column are nonzero, use M_Banded, otherwise use M_Column_Sparse.
  real, parameter, private :: SPARSITY = 0.33     ! If a full matrix has
    ! a greater fraction of nonzeroes than specified by this number, there's
    ! no point in making it sparse.
  logical, save :: CHECKBLOCKS = .false.
  integer, save :: SUBBLOCKLENGTH = 32768

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------  Add_Matrix_Blocks  -----
  function Add_Matrix_Blocks ( XB, YB, SCALE ) result ( ZB )
  ! ZB = XB + [SCALE *] YB
    type(MatrixElement_T), intent(in), target :: XB, YB
    real(r8), intent(in), optional :: SCALE
    type(MatrixElement_T) :: ZB

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: K
!   integer :: I, J, L, N
    real(r8) :: S                            ! My copy of Scale
    type(MatrixElement_T), pointer :: X, Y
    real(kind(zb%values)), pointer :: Z(:,:) ! May be used if Y is col-sparse/banded
    real(kind(zb%values)), pointer :: W(:,:) ! May be used if X is banded
    real(rm), dimension(:,:), pointer :: XD, YD ! For testing
    logical :: OK                       ! For testing

    call nullifyMatrix ( zb ) ! for Sun's still useless compiler
    s = 1.0
    if ( present(scale) ) s = scale

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
    if ( xb%nRows /= yb%nRows .or. xb%nCols /= yb%nCols ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix sizes incompatible in Add_Matrix_Blocks" )
    select case ( x%kind )
    case ( M_Banded )
      select case ( y%kind )
      case ( M_Banded )                          ! X banded, Y banded
        ! Used to be intelligent (below), but didn't work. 
        ! Do it the slow way for now.
        nullify ( z, w )
        call allocate_test ( z, x%nRows, x%nCols, 'Z in Add_Matrix_Blocks', &
          & ModuleName )
        call allocate_test ( w, y%nRows, y%nCols, 'W in Add_Matrix_Blocks', &
          & ModuleName )
        call densify ( z, x )
        call densify ( w, y )
        z = z + s * w
        call sparsify ( z, zb, 'Z in Add_matrix_block', ModuleName )
        call Deallocate_test ( w, 'W in Add_Matrix_Blocks', ModuleName )
        
! OLD CODE THAT DIDN'T WORK
!         call allocate_test ( zb%r1, zb%ncols, "zb%r1", ModuleName )
!         zb%r1 = min(x%r1, y%r1)                  ! First nonzero row
!         call allocate_test ( zb%r2, zb%ncols, "zb%r2", ModuleName, lowBound=0 )
!         zb%r2(0) = 0
!         do k = 1, size(x%r1)                     ! Calculate size of Values
!           zb%r2(k) = zb%r2(k-1) + &
!             & max( x%r2(k)-x%r2(k-1)+x%r1(k), y%r2(k)-y%r2(k-1)+y%r1(k) ) - &
!             & zb%r1(k) + 1
!         end do
!         call allocate_test ( zb%values, zb%r2(size(x%r1)), 1, "zb%values", &
!           & ModuleName )
!         zb%values = 0.0_rm ! ??? Improve this by only filling
!         !                    ??? values that don't get set below
!         do k = 1, size(x%r1)
!           i = 1; j = 1; l = 1
!           n = y%r1(k) - x%r1(k)
!           if ( n > 0 ) then ! Copy rows in X that are not in Y
!             zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
!               & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1)
!             i = i + n
!           else if ( n < 0 ) then ! Copy rows in Y that are not in X
!             n = - n
!             zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
!               & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
!             j = j + n
!           end if
!           l = l + n
!           n = min( x%r1(k) + x%r2(k) - x%r2(k-1) - 1, &
!             &      y%r1(k) + y%r2(k) - y%r2(k-1) - 1 ) - &
!             &      max( x%r1(k), y%r1(k) )
!           zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
!             & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1) + &
!             & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
!           i = i + n; j = j + n; l = l + n
!           n = y%r1(k) + y%r2(k) - y%r2(k-1) - &
!             & ( x%r1(k) + x%r2(k) - x%r2(k-1) )
!           if ( n > 0 ) then ! Copy rows of Y that are not in X
!             zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
!               & y%values(y%r2(k-1)+j:y%r2(k-1)+j-1+n, 1)
!           else if ( n < 0 ) then ! Copy rows in X that are not in Y
!             n = - n
!             zb%values(zb%r2(k-1)+l:zb%r2(k-1)+l-1+n, 1) = &
!               & x%values(x%r2(k-1)+i:x%r2(k-1)+i-1+n, 1)
!           end if
!         end do
      case ( M_Column_sparse )                   ! X banded, Y col sparse
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, y%nRows, y%nCols, "Z in Add_Matrix_Blocks", &
          & ModuleName )
        call densify ( z, y )                    ! z = y%values
        if ( present(scale) ) z = s * z
        do k = 1, x%nCols
          z(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) = &
            & z(x%r1(k): x%r1(k) + x%r2(k) - x%r2(k-1) - 1, k) + &
              & x%values(x%r2(k-1)+1:x%r2(k), 1)
        end do
        call sparsify ( z, zb, "Z in Add_Matrix_Blocks", ModuleName ) ! Zb = Z
      case ( M_Full )                            ! X banded, Y full
        call CopyBlock ( zb, y )                 ! Zb = y
        if ( present(scale) ) zb%values = s * zb%values
        do k = 1, x%nCols
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
        call allocate_test ( z, y%nRows, y%nCols, "Z in Add_Matrix_Blocks", &
          & ModuleName )
        call densify ( z, y )                    ! z = y%values
        if ( present(scale) ) z = s * z
        do k = 1, x%nCols
          z(x%r2(x%r1(k-1)+1:x%r1(k)), k) = &
            & z(x%r2(x%r1(k-1)+1:x%r1(k)), k) + &
              & x%values(x%r1(k-1)+1:x%r1(k),1)
        end do
        call sparsify ( z, zb, "Z in Add_Matrix_Blocks", ModuleName ) ! Zb = Z
      case ( M_Full )                            ! X col sparse, Y full
        call CopyBlock ( zb, y )                 ! Zb = y
        if ( present(scale) ) zb%values = s * zb%values

! Commented-out on account of internal NAG v4.0 bug
         do k = 1, x%nCols
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
        zb%values = x%values + s * y%values
      end select
    end select

    if ( checkBlocks ) then
      nullify ( xd, yd )
      call Allocate_test ( xd, xb%nRows, xb%nCols, 'xd', ModuleName )
      call Allocate_test ( yd, yb%nRows, yb%nCols, 'yd', ModuleName )
      call Densify ( xd, xb )
      call Densify ( yd, yb )
      call TestBlock ( zb, xd+yd, ok, 'Add_Matrix_Blocks', xb%kind, yb%kind )
      if ( .not. ok ) then
        call dump ( xb, name='xb', details=2 )
        call dump ( yb, name='yb', details=2 )
      end if
      call Deallocate_test ( xd, 'xd', ModuleName )
      call Deallocate_test ( yd, 'yd', ModuleName )
    end if
  end function Add_Matrix_Blocks

  ! ---------------------------------  Add_Matrix_Blocks_Unscaled  -----
  function Add_Matrix_Blocks_Unscaled ( XB, YB ) result ( ZB ) ! ZB = XB + YB
    type(MatrixElement_T), intent(in), target :: XB, YB
    type(MatrixElement_T) :: ZB

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    call nullifyMatrix ( zb ) ! for Sun's still useless compiler
                         ! not sure it's needed here
    zb = add_matrix_blocks ( xb, yb )

  end function Add_Matrix_Blocks_Unscaled 

  ! ------------------------------------------------  AssignBlock  -----
  subroutine AssignBlock ( Z, X )
  ! Destroy Z, then copy X to it, using pointer assignment for pointer
  ! components.  Other than the "Destroy Z" part, the semantics are the
  ! same as for intrinsic assignment.  If one has Z = X in a loop, it
  ! is therefore necessary only to destroy Z after the loop.  Notice
  ! that CopyBlock does a deep copy.
    type(MatrixElement_T), intent(inout) :: Z
    type(MatrixElement_T), intent(in) :: X
    call destroyBlock ( z )
    z%kind = x%kind
    z%nRows = x%nRows; z%nCols = x%nCols
    z%r1 => x%r1
    z%r2 => x%r2
    z%values => x%values
  end subroutine AssignBlock

  ! ------------------------------------------- CheckIntegrity_0 -------
  logical function CheckIntegrity_0 ( block, noError )
    type ( MatrixElement_T), intent(in) :: BLOCK
    logical, optional, intent(in) :: NOERROR

    ! Local variables
    integer :: MESSAGETYPE
    integer :: I                        ! Loop counter
    integer :: N                        ! No entries in a column
    
    ! Executable code

    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    checkIntegrity_0 = .true.

    ! Check dimensions
    if ( block%nRows < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Block is negative nRows' )
      checkIntegrity_0 = .false.
    end if
    if ( block%nCols < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Block is negative nCols' )
      checkIntegrity_0 = .false.
    end if

    ! Check associated stuff
    if ( block%kind == m_absent .or. block%kind == m_unknown ) then
      if ( associated ( block%values ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is absent/unknown yet has values associated' )
        checkIntegrity_0 = .false.
      end if
    else
      if ( .not. associated ( block%values ) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Block has no values associted, but not absent/unkown')
        checkIntegrity_0 = .false.
      end if
      if ( any ( lbound ( block%values ) /= (/1,1/) ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is had bad lbound on values.' )
        checkIntegrity_0 = .false.
      end if
    end if
    if ( any ( block%kind == (/ m_absent, m_unknown /) ) ) then
      if ( associated ( block%R1 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is absent/unknown yet has R1 associated' )
        checkIntegrity_0 = .false.
      end if
      if ( associated ( block%R2 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is absent/unknown yet has R2 associated' )
        checkIntegrity_0 = .false.
      end if
    end if
    if ( block%kind == m_full ) then
      if ( size ( block%R1 ) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is full yet has non zero R1 allocated' )
        checkIntegrity_0 = .false.
      end if
      if ( size ( block%R2 ) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is full yet has non zero R2 allocted' )
        checkIntegrity_0 = .false.
      end if
    end if

    if ( block%kind == m_banded .or. block%kind == m_column_sparse ) then
      if ( .not. associated ( block%r1 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is banded/sparse yet has no R1 associated' )
        checkIntegrity_0 = .false.
      end if
      if ( .not. associated ( block%r2 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block is banded/sparse yet has no R2 associated' )
        checkIntegrity_0 = .false.
      end if
    end if

    if ( block%kind == m_banded ) then
      ! Check lbound
      if ( lbound ( block%r1,1 ) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad lbound on R1' )
        checkIntegrity_0 = .false.
      end if
      if ( lbound ( block%r2,1 ) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad lbound on R2' )
        checkIntegrity_0 = .false.
      end if
      ! Check ubound
      if ( ubound ( block%r1,1 ) /= block%nCols ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad ubound on R1' )
        checkIntegrity_0 = .false.
      end if
      if ( ubound ( block%r2,1 ) /= block%nCols ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad ubound on R2' )
        checkIntegrity_0 = .false.
      end if
      ! Check values
      if ( block%r2(0) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad first value in R2' )
        checkIntegrity_0 = .false.
      end if
      if ( any ( (/ block%r1(:), block%r2(:) /) < 0 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has -ve value(s) in R1/R2' )
        checkIntegrity_0 = .false.
      end if
      if ( any ( (/ block%r2(:) /) > ubound ( block%values,1 ) ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has too large value(s) in R2' )
        checkIntegrity_0 = .false.
      end if
      do i = 1, block%nCols
        n =  block%r2(i) - block%r2(i-1)
        if ( n < 0 ) then
          call MLSMessage ( messageType, ModuleName, &
            & 'Banded block has too small a delta in R2' )
          checkIntegrity_0 = .false.
        end if
        if ( n > 0 ) then
          if ( block%r1(i) + n > block%nRows ) then
            call MLSMessage ( messageType, ModuleName, &
              & 'Banded block has too large a delta in R2 or R1 value' )
            checkIntegrity_0 = .false.
          end if
          if ( block%r1(i) == 0 ) then
            call MLSMessage ( messageType, ModuleName, &
              'Banded block contains 0 for R1 in non empty column' )
          end if
        end if
      end do
      if ( ubound(block%values,2) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Banded block has bad ubound(2) for values' )
        checkIntegrity_0 = .false.
      end if
    end if

    if ( block%kind == m_column_sparse ) then
      ! Check lbound
      if ( lbound ( block%r1,1 ) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad lbound on R1' )
        checkIntegrity_0 = .false.
      end if
      if ( lbound ( block%r2,1 ) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad lbound on R2' )
        checkIntegrity_0 = .false.
      end if
      ! Check ubound
      if ( ubound ( block%r1,1 ) /= block%nCols ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad ubound on R1' )
        checkIntegrity_0 = .false.
      end if
      if ( ubound ( block%r2,1 ) /= ubound(block%values,1) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad ubound on R2' )
        checkIntegrity_0 = .false.
      end if
      ! Check values
      if ( block%r1(0) /= 0 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad first value in R1' )
        checkIntegrity_0 = .false.
      end if
      if ( any ( (/ block%r1(1:), block%r2(:) /) < 1 ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has 0/-ve value(s) in R1/R2' )
        checkIntegrity_0 = .false.
      end if
      if ( any ( (/ block%r1(1:) /) > ubound(block%r2,1) ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has too large value(s) in R1' )
        checkIntegrity_0 = .false.
      end if
      if ( any ( (/ block%r2(:) /) > block%nRows ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has too large value(s) in R2' )
        checkIntegrity_0 = .false.
      end if
      do i = 1, block%nCols
        n =  block%r1(i) - block%r1(i-1) - 1
        if ( n < 0 ) then
          call MLSMessage ( messageType, ModuleName, &
            & 'Column sparse block has too small a delta in R1' )
          checkIntegrity_0 = .false.
        end if
        if ( block%r1(i) + n > ubound ( block%r2, 1 ) ) then
          call MLSMessage ( messageType, ModuleName, &
            & 'Column sparse block has too large a delta in R1' )
          checkIntegrity_0 = .false.
        end if
      end do
      if ( ubound(block%values,2) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Column sparse block has bad ubound(2) for values' )
        checkIntegrity_0 = .false.
      end if
    end if

    if ( block%kind == m_full ) then 
      if ( ubound ( block%values,1 ) /= block%nRows ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block nRows inconsistent with dimensions of values' )
        checkIntegrity_0 = .false.
      end if
      if ( ubound ( block%values,2 ) /= block%nCols ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Block nCols inconsistent with dimensions of values' )
        checkIntegrity_0 = .false.
      end if
    end if

  end function CheckIntegrity_0

  ! -------------------------------------------  CholeskyFactor_0  -----
  subroutine CholeskyFactor_0 ( Z, XOPT, STATUS )
  ! If XOPT is present compute Z such that Z^T Z = XOPT and Z is upper-
  ! triangular. Otherwise, replace Z such that Z(output)^T Z(output) =
  ! Z(input) and Z(output) is upper-triangular.
    type(MatrixElement_T), target, intent(inout) :: Z
    type(MatrixElement_T), target, intent(in), optional :: XOPT
    integer, intent(out), optional :: STATUS

    real(rm) :: D             ! Diagonal(I,I) element of Z
    real(rm) :: G             ! X(I,J) - dot_product(X(1:i-1,i),X(1:i-1,j))
    integer :: I, J           ! Subscripts and loop inductors
    integer :: II, IJ         ! Subscripts in VALUES for I,I and I,J components
    !                           in the case of M_Banded
    integer :: NC             ! Number of columns
    integer, pointer, dimension(:) :: R1      ! First nonzero row of Z (Banded)
    integer :: RZ             ! Starting row in Z for inner product (Banded)
    real(rm), save :: TOL = -1.0_rm
    real(rm), pointer, dimension(:,:) :: XIN  ! A pointer to the input,
    !                           data, or a densified copy of it
    type(MatrixElement_T), pointer :: X       ! XOPT or Z, depending on whether
    !                           XOPT is present or absent, respectively.
    real(rm), pointer, dimension(:,:) :: ZT   ! A local full result that is
    !                           sparsified at the end.
    real(rm), pointer, dimension(:,:) :: TST1, TST2 ! When using checkblock
    logical :: OK                       ! For testing

    if ( tol < 0.0_rm ) tol = sqrt(tiny(0.0_rm))
    nullify ( r1, xin, zt )
    x => z
    if ( present(xopt) ) x => xopt
    nc = x%nCols
    if ( nc /= x%nRows )&
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
          if ( present(status) ) then
            status = i
            return
          end if
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        end if
        zt(i,1:i-1) = 0.0_rm  ! Clear left from the diagonal (helps Sparsify!)
        g = x%values(ii+x%r2(i-1)+1,1) - &
            & dot( i-r1(i), zt(r1(i),i), 1, zt(r1(i),i), 1 )
!       g = x%values(ii+x%r2(i-1)+1,1) - &
!           & dot_product( zt(r1(i):i-1,i), zt(r1(i):i-1,i) )
        if ( g <= tol .and. i < nc ) then
          if ( present(status) ) then
            status = i
            return
          end if
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Matrix in CholeskyFactor is not positive-definite." )
        end if
        d = sqrt(g)
        zt(i,i) = d
!$OMP PARALLEL DO private ( ij, rz, g )
        do j = i+1, nc
          ij = i - x%r1(j)    ! Offset in VALUES of (I,J) element
          rz = max(r1(i),r1(j))
          if ( i <= rz ) then
            g = 0.0
          else
            g = - dot( i-rz, zt(rz,i), 1, zt(rz,j), 1 )
!           g = - dot_product( zt(rz:i-1,i), zt(rz:i-1,j) )
          end if
          if ( ij >= 0 .and. ij <= x%r2(j) - x%r2(j-1) ) &
            & g = x%values(ij+x%r2(j-1)+1,1) + g
          zt(i,j) = g / d
          if ( abs(zt(i,j)) >= tiny(0.0_rm) ) r1(j) = min( r1(j), i )
        end do ! j
!$OMP END PARALLEL DO
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
      if ( present ( status ) ) status = 0
    case ( M_Column_Sparse )
      call allocate_test ( zt, nc, nc, "ZT for CholeskyFactor", ModuleName )
      call allocate_test ( xin, nc, nc, "XIN for CholeskyFactor", ModuleName )
      ! ??? Densify and then compute Cholesky decomposition of dense block.
      ! ??? If necessary, improve this in level 1.0 by working directly
      ! ??? with sparse input.
      call densify ( xin, x )
      call denseCholesky ( zt, xin, status )
      if ( present ( status ) ) then
        if ( status /= 0 ) return
      end if
      call sparsify ( zt, z, "ZT in CholeskyFactor", ModuleName ) ! Z := Zt
      call deallocate_test ( xin, "XIN in CholeskyFactor", ModuleName )
    case ( M_Full )
      if ( .not. associated(x,z) ) then
        call createBlock ( z, nc, nc, M_Full )
      end if
      call denseCholesky ( z%values, x%values, status )
    end select

    if ( checkBlocks ) then
      nullify ( tst1, tst2 )
      call Allocate_test ( tst1, z%nRows, z%nCols, 'tst1', ModuleName )
      call Allocate_test ( tst2, z%nRows, z%nCols, 'tst2', ModuleName )
      call densify ( tst1, x )
      call denseCholesky ( tst2, tst1 )
      call testBlock ( z, tst2, ok, 'CholeskyFactor', z%kind )
      if ( .not. ok ) then
        call dump ( z, name='z', details=2 )
      end if
      call Deallocate_test ( tst1, 'tst1', ModuleName )
      call Deallocate_test ( tst2, 'tst2', ModuleName )
    end if
  end subroutine CholeskyFactor_0

  ! ------------------------------------------------  ClearRows_0  -----
  subroutine ClearRows_0 ( X, MASK )
  ! Clear the rows of X for which MASK has a nonzero M_LinAlg bit.
    type(MatrixElement_T), intent(inout) :: X
    character, dimension(1:), intent(in) :: MASK
    integer :: I, J                ! Subscripts and row indices
    select case ( x%kind )
    case ( M_Absent )
    case ( M_Banded )              ! ??? Adjust the sparsity representation ???
      do j = 1, x%nCols
        do i = x%r1(j), x%r1(j) + x%r2(j) - x%r2(j-1) - 1 ! row numbers
          if ( iand( ichar(mask(i)), m_LinAlg ) /= 0 ) &
            & x%values(x%r2(j-1) + i - x%r1(j) + 1, 1) = 0.0_rm
        end do ! i
      end do ! j = 1, x%nCols
    case ( M_Column_Sparse )       ! ??? Adjust the sparsity representation ???
      do j = 1, x%nCols
        do i = x%r1(j-1)+1, x%r1(j)
          if ( iand( ichar(mask(x%r2(i))), m_LinAlg ) /= 0 ) &
            & x%values(j,1) = 0.0_rm
        end do ! i
      end do ! j = 1, x%nCols
    case ( M_Full )
      do i = lbound(mask,1), ubound(mask,1)
        if ( i > x%nRows ) return
        if ( iand( ichar(mask(i)), m_LinAlg ) /= 0 ) &
          & x%values(i,:) = 0.0_rm
      end do ! i = lbound(mask,1), ubound(mask,1)
    end select
  end subroutine ClearRows_0

  ! -------------------------------------------------  CloneBlock  -----
  subroutine CloneBlock ( Z, X ) ! Z = X, except the values
  ! Duplicate a matrix block, including copying all of its structural
  ! descriptive information, but not its values.
    type(MatrixElement_T), intent(inout) :: Z ! intent(inout) so that
      !                            destroyBlock gets a chance to clean up surds
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

  ! -------------------------------------------------  ColumnScale_0_r4  -----
  subroutine ColumnScale_0_r4 ( X, V, NEWX ) ! Z = X V where V is a diagonal
  !                                         matrix represented by a vector
  !                                         and Z is either X or NEWX.
    type(MatrixElement_T), intent(inout), target :: X
    real(r4), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(inout), target, optional :: NEWX ! intent(inout)
      !                            so that the destroyBlock in cloneBlock
      !                            gets a chance to clean up surds
    type(MatrixElement_T), pointer :: Z
    integer :: I

    include "columnscale_0.f9h"
  end subroutine ColumnScale_0_r4

  ! -------------------------------------------------  ColumnScale_0_r8  -----
  subroutine ColumnScale_0_r8 ( X, V, NEWX ) ! Z = X V where V is a diagonal
  !                                         matrix represented by a vector
  !                                         and Z is either X or NEWX.
    type(MatrixElement_T), intent(inout), target :: X
    real(r8), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(inout), target, optional :: NEWX ! intent(inout)
      !                            so that the destroyBlock in cloneBlock
      !                            gets a chance to clean up surds
    type(MatrixElement_T), pointer :: Z

    integer :: I
    include "columnscale_0.f9h"
  end subroutine ColumnScale_0_r8

  ! -----------------------------------------------------  Col_L1  -----
  real(rm) function Col_L1 ( X, N )
  ! Return the L1 norm of column N of X
    type(MatrixElement_T), intent(in) :: X
    integer, intent(in) :: N
    integer :: I
    col_l1 = 0.0_rm
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
      do i = 1, x%nRows
        col_l1 = col_l1 + abs(x%values(i,n))
      end do
    end select
  end function Col_L1

  ! --------------------------------------------------  CopyBlock  -----
  subroutine CopyBlock ( Z, X ) ! Destroy Z, deep Z = X, including the values
    type(MatrixElement_T), intent(inout) :: Z ! intent(inout) so that the
      !                            destroyBlock in cloneBlock gets a chance
      !                            to clean up surds
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
    case ( M_Absent, M_Unknown )
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
    real(rm), intent(inout) :: ZT(:,:) ! Inout in case it's associated with XIN
    real(rm), intent(inout) :: XIN(:,:) ! Inout in case it's associated with Z
    integer, intent(out), optional :: Status

    real(rm) :: D
    real(rm), save :: TOL = -1.0_rm
    integer :: I, J, NC

    nc = size(xin,2)
    if ( tol < 0.0_rm ) tol = sqrt(tiny(0.0_rm))
    do i = 1, nc
      zt(i+1:nc,i) = 0.0_rm ! Clear below the diagonal (helps Sparsify!)
      d = xin(i,i) - dot( i-1, zt(1,i), 1, zt(1,i), 1 )
!     d = xin(i,i) - dot_product( zt(1:i-1,i), zt(1:i-1,i) )
      if ( (d <= tol .and. i < nc) .or. (d < 0.0) ) then
        if ( present(status ) ) then
          status = i
          return
        end if
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Matrix in DenseCholesky is not positive-definite." )
      end if
      d = sqrt(d)
      zt(i,i) = d
!$OMP PARALLEL DO
      do j = i+1, nc
        zt(i,j) = ( xin(i,j) - dot( i-1, zt(1,i), 1, zt(1,j), 1 ) ) / d
!       zt(i,j) = ( xin(i,j) - dot_product( zt(1:i-1,i), zt(1:i-1,j) ) ) / d
      end do ! j
!$OMP END PARALLEL DO
    end do ! i
    if ( present(status) ) status = 0
  end subroutine DenseCholesky

  ! ----------------------------------------------------  DensifyA  -----
  subroutine DensifyA ( Z, B )
  ! Given a matrix block B, produce a full matrix Z, even if the matrix
  ! block had a sparse representation.
    real(rm), intent(out) :: Z(:,:)          ! Full matrix to produce
    type(MatrixElement_T), intent(in) :: B   ! Input matrix block
    integer :: I                             ! Column index

    if ( size(z,1) /= b%nRows .or. size(z,2) /= b%nCols ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Incompatible shapes in Densify" )
    end if
    select case ( b%kind )
    case ( M_Absent )
      z = 0.0_rm
    case ( M_Banded )
      z = 0.0_rm
      do i = 1, b%nCols
        z(b%r1(i):b%r1(i)+b%r2(i)-b%r2(i-1)-1,i) = &
          & b%values(b%r2(i-1)+1:b%r2(i),1)
      end do ! i
    case ( M_Column_Sparse )
      z = 0.0_rm
      do i = 1, b%nCols
        z(b%r2(b%r1(i-1)+1:b%r1(i)),i) = b%values(b%r1(i-1)+1:b%r1(i),1)
      end do ! i
    case ( M_Full )
      z = b%values
    end select
  end subroutine DensifyA

  ! ---------------------------------------------  DensifyB ------------
  subroutine DensifyB ( B )
    ! Densify a block 'in place'
    type ( MatrixElement_T ), intent(inout) :: B
    ! Local variables
    real(rm), dimension(:,:), pointer :: Z
    ! Executable code
    if ( b%kind == m_banded .or. b%kind == m_column_sparse ) then
      nullify ( z )
      call Allocate_test ( z, b%nRows, b%nCols, 'Z', ModuleName )
      call Densify ( z, b )
      call Deallocate_test ( b%values, 'b%values', ModuleName )
      call Deallocate_test ( b%r1, 'b%r1', ModuleName )
      call Deallocate_test ( b%r2, 'b%r2', ModuleName )
      b%values => z
      b%kind = m_full
    end if
  end subroutine DensifyB

  ! ---------------------------------------------  DestroyBlock_0  -----
  subroutine DestroyBlock_0 ( B )
    ! Deallocate the pointer components of the matrix block B.  Change
    ! its kind to M_Absent.  Don't clobber b%nRows and b%nCols.
    type(MatrixElement_T), intent(inout) :: B
    ! Don't bother to destroy absent blocks
    if (b%kind /= M_Absent) then
      call deallocate_test ( b%r1, "b%r1", ModuleName )
      call deallocate_test ( b%r2, "b%r2", ModuleName )
      call deallocate_test ( b%values, "b%values", ModuleName )
      b%kind = m_absent
    end if
  end subroutine DestroyBlock_0

  ! ----------------------------------------------  GetDiagonal_0_r4  -----
  subroutine GetDiagonal_0_r4 ( B, X, SquareRoot )
  ! Get the diagonal elements of B into X.  Return the square root of the
  ! diagonal elements if SquareRoot is present and true.
    type(MatrixElement_T), intent(in) :: B
    real(r4), dimension(:), intent(out) :: X
    logical, intent(in), optional :: SquareRoot
    integer :: I, J, N
 
    include "getdiagonal_0.f9h"
  end subroutine GetDiagonal_0_r4

  ! ----------------------------------------------  GetDiagonal_0_r8  -----
  subroutine GetDiagonal_0_r8 ( B, X, SquareRoot )
  ! Get the diagonal elements of B into X.  Return the square root of the
  ! diagonal elements if SquareRoot is present and true.
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(out) :: X
    logical, intent(in), optional :: SquareRoot
    integer :: I, J, N

    include "getdiagonal_0.f9h"
  end subroutine GetDiagonal_0_r8

  ! -----------------------------------------  GetMatrixElement_0  -----
  real(rm) function GetMatrixElement_0 ( Matrix, Row, Col )
  ! Get the (row,col) element of Matrix
    type(matrixElement_T), intent(in) :: Matrix
    integer, intent(in) :: Row, Col
    integer :: J
    select case ( matrix%kind )
    case ( m_absent )
      getMatrixElement_0 = 0.0_rm
    case ( m_banded )
      if ( row < matrix%r1(col) .or. &
        &  row >= matrix%r1(col) + matrix%r2(col)-matrix%r2(col-1) ) then
        getMatrixElement_0 = 0.0_rm
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
      getMatrixElement_0 = 0.0_rm
    case ( m_full )
      getMatrixElement_0 = matrix%values(row,col)
    end select
  end function GetMatrixElement_0

  ! --------------------------------------  GetVectorFromColumn_0_r4  -----
  subroutine GetVectorFromColumn_0_r4 ( B, Column, X )
  ! Fill the vector X from "Column" of B.
    type(MatrixElement_T), intent(in) :: B
    integer, intent(in) :: Column
    real(r4), dimension(:), intent(out) :: X

    include "getvectorfromcolumn_0.f9h"
  end subroutine GetVectorFromColumn_0_r4

  ! --------------------------------------  GetVectorFromColumn_0_r8  -----
  subroutine GetVectorFromColumn_0_r8 ( B, Column, X )
  ! Fill the vector X from "Column" of B.
    type(MatrixElement_T), intent(in) :: B
    integer, intent(in) :: Column
    real(r8), dimension(:), intent(out) :: X

    include "getvectorfromcolumn_0.f9h"
  end subroutine GetVectorFromColumn_0_r8

  ! -------------------------------------------  InvertCholesky_0  -----
  subroutine InvertCholesky_0 ( U, UI, Square, Status )
  ! Given the upper-triangular matrix U, produce its inverse in UI.
  ! If Square is present and > 0, compute UI = (U**T U)**-1 = U**-1 U**-T.
  ! If Square is present and < 0, compute UI = diag((U**T U)**-1 = U**-1 U**-T).
    type(MatrixElement_T), intent(inout) :: U  ! inout to use InvertDenseCholesky
    type(MatrixElement_T), intent(inout) :: UI ! inout to look at nRows etc.
    integer, intent(in), optional :: Square  ! 0 => nothing special
                                             ! >0 => UI := U**-1 U**-T
                                             ! <0 => UI := diag(U**-1 U**-T)
    integer, intent(out), optional :: Status ! -1 if the kind is M_Absent, else
    !                                   ! Index of zero diagonal, else 0

    integer :: MySquare, MyStatus
    real(rm), pointer, dimension(:,:) :: myUI

    mySquare = 0
    if ( present(square) ) mySquare = square
    myStatus = 0
    nullify ( myUI )
    select case ( u%kind )
    case ( M_Absent )
      myStatus = -1
    case ( M_Banded )
      call allocate_test ( myUI, u%nRows, u%nCols, "myUI in InvertCholesky_0", &
        & moduleName )
      call densify ( myUI, U )
      call invertDenseCholesky ( myUI, myUI, square, .true., myStatus )
      if ( myStatus == 0 ) then
        if ( mySquare < 0 ) then
          call storeDiag
        else
          call sparsify ( myUI, UI, "myUI in InvertCholesky_0", moduleName )
        end if
      end if
    case ( M_Column_Sparse )
      call allocate_test ( myUI, u%nRows, u%nCols, "myUI in InvertCholesky_0", &
        & moduleName )
      call invertDenseCholesky ( myUI, myUI, square, .true., myStatus )
      if ( myStatus == 0 ) then
        if ( mySquare < 0 ) then
          call storeDiag
        else
          call sparsify ( myUI, UI, "myUI in InvertCholesky_0", moduleName )
        end if
      end if
    case ( M_Full )
      if ( UI%kind /= M_Full .or. &
        & UI%nRows /= u%nRows .or. UI%nCols /= u%nCols ) then
        call createBlock ( UI, u%nRows, u%nCols, M_Full )
      end if
      call invertDenseCholesky ( u%values, UI%values, square, .true., myStatus )
      if ( myStatus == 0 .and. mySquare < 0 ) then ! just computing the diagonal
        myUI => UI%values
        nullify ( UI%values )
        call storeDiag
      end if
    end select

    call deallocate_test ( myUI, "myUI in InvertCholesky_0", moduleName )
    if ( present(status) ) status = myStatus

  contains
    subroutine StoreDiag
      integer :: I
      call createBlock ( UI, u%nRows, u%nCols, M_banded, &
        & numberNonzero = u%nRows, bandHeight = 1 )
      do i = 1, u%nRows
        ui%values(i,1) = myUI(i,i)
      end do
    end subroutine StoreDiag
  end subroutine InvertCholesky_0

  ! --------------------------------------  InvertDenseCholesky_0  -----
  subroutine InvertDenseCholesky_0 ( U, UI, Square, Clear, Status )
  ! Given the upper-triangular matrix U, produce its inverse, UI.  It's OK
  ! if UI and U are the same matrix (see the algorithm description below).
  ! If Square is present and > 0, compute UI = (U**T U)**-1 = U**-1 U**-T.
  ! If Square is present and < 0, compute UI = diag((U**T U)**-1 = U**-1 U**-T).
  ! If Clear is present and true, clear below the diagonal.
    real(rm), intent(inout) :: U(:,:)   ! inout in case U and UI have the
    real(rm), intent(inout) :: UI(:,:)  ! same associated actual argument
    integer, intent(in), optional :: Square  ! 0 => nothing special
                                             ! >0 => UI := U**-1 U**-T
                                             ! <0 => UI := diag(U**-1 U**-T)
    logical, intent(in), optional :: Clear   ! Clear below diagonal
    integer, intent(out), optional :: Status ! Index of zero diagonal, else 0

  !{Let $u_{ij}$ be an element of $\bf U$ and $t_{ij}$ be an element of
  ! ${\bf U}^{-1}$. To invert $\bf U$, solve $\sum_{l=i}^j t_{il} u_{lj} =
  ! \delta_{ij}$ for $t_{ij}$.  (The index of summation $l$ runs from $i$ to
  ! $j$ instead of from $1$ to $n$ because we know that $\bf U$ and ${\bf
  ! U}^{-1}$ are triangular.)  This gives $t_{ii} = u_{ii}^{-1}$ (for $i =
  ! j$) and  $t_{ij} = - \left ( \sum_{l=i}^{j-1} t_{il} u_{lj} \right )
  ! u_{jj}^{-1}$ (for $i < j$).  We compute the diagonal elements of
  ! $U^{-1}$ first, giving the following algorithm:
  !
  ! \hspace*{0.5in}$t_{ii} := u_{ii}^{-1}~~ i = 1,~ \dots,~ n$\\
  ! \hspace*{0.5in}$\left \{ t_{ij} := -\left ( \sum_{l=i}^{j-1} t_{il} u_{lj}
  !   \right ) t_{jj}~~ j = i+1,~ \dots,~ n \right \}~~ i = 1,~ \dots,~ n-1$

    integer :: I, J, N
    character(len=6) :: Where ! for an output message

    ! Check dimensions
    if ( any(shape(u) /= shape(ui)) ) call MLSMessage ( MLSMSG_Error, &
      & moduleName, 'Matrices in InvertCholesky_0 are not the same shape' )

    if ( present(status) ) status = 0

    ! Invert the diagonal elements
    n = size(u,1)
    do j = 1, n
      if ( abs(U(j,j)) <= tiny(0.0_rm) ) then
        if ( present(status) ) then
          status = j
          return
        end if
        write ( where, * ) j
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Matrix in InvertCholesky_0 is singular at ' // &
          & trim(adjustl(where)))
        exit
      end if
      UI(j,j) = 1.0_rm / U(j,j)
    end do

    ! Finish inverting the rest of the matrix.
    do i = 1, n-1
      do j = i+1, n
        UI(i,j) = -dot(j-i, UI(i,i), n, U(i,j), 1) * UI(j,j)
!       UI(i,j) = -dot_product(UI(i,i:j-1),U(i:j-1,j)) * UI(j,j)
      end do
    end do

    if ( present(clear) ) then
      if ( clear ) then            ! Clear below the diagonal
        do j = 1, n
          ui(j+1:,j) = 0.0_rm
        end do
      end if
    end if

    if ( .not. present(square) ) return
    select case ( square )
    case ( 1: )
      ! Replace the upper triangle of UI with the upper triangle of UI * UI**T
      do i = 1, n
        do j = i, n
          UI(i,j) = dot(n-j+1, UI(i,j), n, UI(j,j), n)
!         UI(i,j) = dot_product(UI(i,j:n),UI(j,j:n))
        end do
      end do
    case ( :-1 )
      ! Replace the diagonal of UI with diag(UI * UI**T) and clear the
      ! rest of it.
      do i = 1, n
        ui(i,i) = dot(n-i+1, UI(i,i), n, UI(i,i), n)
!       ui(i,i) = dot_product(UI(i,i:n),ui(i,i:n))
        ui(i,i+1:) = 0.0_rm
      end do
    end select

  end subroutine InvertDenseCholesky_0

  ! ------------------------------------------------  MaxAbsVal_0  -----
  real(rm) function MaxAbsVal_0 ( B )
  ! Return the magnitude of the element in B that has the largest magnitude.
    type(MatrixElement_T), intent(in) :: B
    if ( b%kind == m_absent ) then
      maxAbsVal_0 = 0.0_rm
    else
      maxAbsVal_0 = maxval(abs(b%values))
    end if
  end function MaxAbsVal_0

!---------------------------------------------  MatrixInversion_0  -----
  subroutine MatrixInversion_0 ( A, Upper )
  ! Invert the symmetric positive-definite matrix A in place. If Upper
  ! is present and true, compute only the upper triangle of the inverse.

  real (rm), dimension(:,:),intent(inout) :: A
  logical, intent(in), optional :: Upper

! real (rm), dimension(:,:), pointer :: U
! real (rm), dimension(:), pointer :: b
! real (rm), dimension(:), pointer :: x
  logical :: MyUpper  
! integer :: I, J, N
  integer :: I, N

  myUpper = .false.
  if ( present ( upper) ) myUpper = upper
  n = size(A,2)

  call denseCholesky ( a, a )           ! Replace A by its Cholesky factor
  call invertDenseCholesky ( a, a, square=1 )     ! Replace A by (A^T A)^{-1}

  if ( myUpper ) return

  ! Fill in the lower triangle from the upper one
  do i = 2, n
    a(i,:i-1) = a(:i-1,i)
  end do

! nullify ( x, b, u )
! call allocate_Test ( x, n, "X in MatrixInversion_0", moduleName )
! call allocate_Test ( b, n, "B in MatrixInversion_0", moduleName )
! call allocate_Test ( u, n, n, "U in MatrixInversion_0", moduleName )

! call DenseCholesky (U, A)

! if ( .not. myUpper ) j = n
! do i = 1, n
!   b = 0.0_rm
!   b(i) = 1.0_rm
!   call SolveCholeskyA_0 ( U, x, b, transpose=.true. )
!   b = x
!   call SolveCholeskyA_0 ( U, x, b, transpose=.false. )
!   if ( myUpper ) j = i
!   A(1:j,i) = x(1:j)
! end do

! call deallocate_Test ( x, "X in MatrixInversion_0", moduleName )
! call deallocate_Test ( b, "B in MatrixInversion_0", moduleName )
! call deallocate_Test ( u, "U in MatrixInversion_0", moduleName )

  end subroutine MatrixInversion_0

  ! --------------------------------------------------  MinDiag_0  -----
  real(rm) function MinDiag_0 ( B )
  ! Return the magnitude of the element on the diagonal of B that has the
  ! smallest magnitude.
    type(MatrixElement_T), intent(in) :: B
    integer :: I, J
    select case ( b%kind )
    case ( M_Absent )
      minDiag_0 = 0.0_rm
    case ( M_Banded )
      minDiag_0 = huge(0.0_rm)
      do i = 1, min(b%nRows,b%nCols)
        if ( b%r1(i) <= i .and. b%r1(i) + b%r2(i) - b%r2(i-1) > i ) &
          & minDiag_0 = min(minDiag_0, abs(b%values(b%r2(i-1)+i-b%r1(i)+1,1)))
      end do ! i
    case ( M_Column_Sparse )
      minDiag_0 = huge(0.0_rm)
      do i = 1, min(b%nRows,b%nCols)
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
      do i = 2, min(b%nRows,b%nCols)
        minDiag_0 = min(minDiag_0, abs(b%values(i,i)))
      end do ! i
    end select
  end function MinDiag_0

  ! ----------------------------------------  MultiplyMatrix_XY_0  -----
  subroutine MultiplyMatrix_XY_0 ( XB, YB, ZB, Update, Subtract )
  ! If Update is absent or false, ZB = +/- XB YB.
  ! Otherwise, ZB = ZB +/- XB YB
  ! The +/- is + unless subtract is present and true.
  !??? Someday, do sparse multiplies in the non-M_Full cases.
    type(MatrixElement_T), intent(in) :: XB, YB
    type(MatrixElement_T), intent(inout) :: ZB
    logical, intent(in), optional :: Update, Subtract

    real(rm) :: Alpha                   ! -1 or 1 depending on subtract
    real(rm) :: Beta                    ! 0 or 1, depending on Update
    logical :: MyX, MyY                 ! Did X or Y result from densify?
    real(rm), pointer, dimension(:,:) :: X, Y, Z  ! Dense matrices

    alpha = 1.0_rm
    if ( present(subtract) ) then
      if ( subtract ) alpha = -1.0_rm
    end if

    beta = 0.0_rm
    if ( present(update) ) then
      if ( update ) beta = 1.0_rm
    end if
    
    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      if ( abs(beta) < 0.5_rm ) call createBlock ( zb, xb%nRows, yb%nCols, M_Absent )
      return
    end if

    if ( xb%nCols /= yb%nRows ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "XB and YB Matrix sizes incompatible in MultiplyMatrix_XY_0" )

    nullify ( z )
    call allocate_test ( z, xb%nRows, yb%nCols, 'Z in MultiplyMatrix_XY_0', &
      & moduleName )
    if ( abs(beta) < 0.5_rm ) then
      z = 0.0_rm
    else
      call densify ( z, zb )
    end if

    myX = xb%kind /= M_Full
    if ( .not. myX ) then
      x => xb%values
    else
      nullify ( x )
      call allocate_test ( x, xb%nRows, xb%nCols, 'X in MultiplyMatrix_XY_0', &
        & moduleName )
      call densify ( x, xb )
    end if

    myY = yb%kind /= M_Full
    if ( .not. myY ) then
      y => yb%values
    else
      nullify ( y )
      call allocate_test ( y, yb%nRows, yb%nCols, 'Y in MultiplyMatrix_XY_0', &
        & moduleName )
      call densify ( y, yb )
    end if

    ! zb%noRows/zb%noCols undefined here, so avoid using them.
    call gemm ( 'N', 'N', xb%nRows, yb%nCols, xb%nCols, alpha, &
      & x, xb%nRows, y, yb%nRows, beta, z, xb%nRows )

    if ( myX .and. myY ) then
      call sparsify ( z, zb, 'Z in MultiplyMatrix_XY_0', moduleName ) ! Zb := Z
    else
      call createBlock ( zb, xb%nRows, yb%nCols, m_full, noValues=.true. )
      zb%values => z
    end if
    if ( myX ) call deallocate_test ( x, 'X in MultiplyMatrix_XY_0', moduleName )
    if ( myY ) call deallocate_test ( y, 'Y in MultiplyMatrix_XY_0', moduleName )

  end subroutine MultiplyMatrix_XY_0

  ! --------------------------------------  MultiplyMatrix_XY_T_0  -----
  subroutine MultiplyMatrix_XY_T_0 ( XB, YB, ZB, Update, Subtract )
  ! If Update is absent or false, ZB = +/- XB YB^T.
  ! Otherwise, ZB = ZB +/- XB YB^T.
  ! The +/- is + unless subtract is present and true.
  !??? Someday, do sparse multiplies in the non-M_Full cases.
    type(MatrixElement_T), intent(in) :: XB, YB
    type(MatrixElement_T), intent(inout) :: ZB
    logical, intent(in), optional :: Update, Subtract

    real(rm) :: Alpha                   ! -1 or 1 depending on subtract
    real(rm) :: Beta                    ! 0 or 1, depending on Update
    logical :: MyX, MyY                 ! Did X or Y result from densify?
    real(rm), pointer, dimension(:,:) :: X, Y, Z  ! Dense matrices


    alpha = 1.0_rm
    if ( present(subtract) ) then
      if ( subtract ) alpha = -1.0_rm
    end if

    beta = 0.0_rm
    if ( present(update) ) then
      if ( update ) beta = 1.0_rm
    end if
    
    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      if ( abs(beta) < 0.5_rm ) call createBlock ( zb, xb%nRows, yb%nRows, M_Absent )
      return
    end if

    if ( xb%nCols /= yb%nCols ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "XB and YB Matrix sizes incompatible in MultiplyMatrix_XY_T_0" )

    nullify ( z )
    call allocate_test ( z, xb%nRows, yb%nRows, 'Z in MultiplyMatrix_XY_T_0', &
      & moduleName )
    if ( abs(beta) < 0.5_rm ) then
      z = 0.0_rm
    else
      call densify ( z, zb )
    end if

    myX = xb%kind /= M_Full
    if ( .not. myX ) then
      x => xb%values
    else
      nullify ( x )
      call allocate_test ( x, xb%nRows, xb%nCols, 'X in MultiplyMatrix_XY_T_0', &
        & moduleName )
      call densify ( x, xb )
    end if

    myY = yb%kind /= M_Full
    if ( .not. myY ) then
      y => yb%values
    else
      nullify ( y )
      call allocate_test ( y, yb%nRows, yb%nCols, 'Y in MultiplyMatrix_XY_T_0', &
        & moduleName )
      call densify ( y, yb )
    end if

    ! zb%nRows/zb%nCols undefined here, so don't use them.
    call gemm ( 'N', 'T', xb%nRows, yb%nRows, xb%nCols, alpha, &
      & x, xb%nRows, y, yb%nRows, beta, z, xb%nRows )

    if ( myX .and. myY ) then
      call sparsify ( z, zb, 'Z in MultiplyMatrix_XY_T_0', moduleName ) ! Zb := Z
    else
      call createBlock ( zb, xb%nRows, yb%nRows, M_Full, noValues=.true. )
      zb%values => z
    endif
    if ( myX ) call deallocate_test ( x, 'X in MultiplyMatrix_XY_T_0', moduleName )
    if ( myY ) call deallocate_test ( y, 'Y in MultiplyMatrix_XY_T_0', moduleName )

  end subroutine MultiplyMatrix_XY_T_0

  ! ---------------------------------------  MultiplyMatrix_XTY_0  -----
  subroutine MultiplyMatrix_XTY_0 ( XB, YB, ZB, UPDATE, SUBTRACT, XMASK, YMASK, &
    &                               UPPER )
  ! ZB = XB^T YB if UPDATE is absent or false and SUBTRACT is absent or false;
  ! ZB = -XB^T YB if UPDATE is absent or false and SUBTRACT is present and
  !      true;
  ! ZB = ZB + XB^T YB if UPDATE is present and true and  SUBTRACT is absent
  !      or false;
  ! ZB = ZB - XB^T YB if UPDATE is present and true and  SUBTRACT is present
  !      and true;
  ! If XMASK (resp. YMASK) is present and associated, ignore columns of XB
  ! (resp. YB) that correspond to elements with a nonzero M_LinAlg bit of
  ! XMASK (resp. YMASK).
  ! If UPPER is present and true, compute only the upper triangle of ZB.
    type(MatrixElement_T), intent(in) :: XB, YB
    type(MatrixElement_T), intent(inout) :: ZB
    logical, intent(in), optional :: UPDATE
    logical, intent(in), optional :: SUBTRACT
    character, optional, pointer, dimension(:) :: XMASK, YMASK ! intent(in)
    logical, intent(in), optional :: UPPER

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! If UPDATE is absent or false, it is important to invoke DestroyBlock
  ! for ZB after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    integer :: I, J, K, L, M, MZ, N
    integer :: XI_1, XI_N, XR_1, XR_N, YI_1, YI_N, YR_1, YR_N, CR_1, CR_N, C_N
    integer :: RS, R0, R1, RN           ! Row indicies used for 'blocking' full/full
    logical :: MY_SUB, MY_UPD, MY_UPPER
    real(rm) :: S                            ! Sign, subtract => -1 else +1
    integer :: XD, YD
    character, pointer, dimension(:) :: XM, YM
    real(rm) :: XY                           ! Product of columns of X and Y
    real(rm), pointer, dimension(:,:) :: Z   ! Temp for sparse * sparse

    real(rm), pointer, dimension(:,:) :: XDNS, YDNS, ZDNS, ZDNS2 ! For checking

    character(len=132) :: LINE          ! A message
    logical :: OK                       ! For testing

    nullify ( xm, ym, z )
    my_upd = .false.
    if ( present(update) ) my_upd = update

    my_sub = .false.
    if ( present(subtract) ) my_sub = subtract
    s = 1.0_rm
    if ( my_sub ) s = -1.0_rm
    if ( present(xmask) ) xm => xmask
    if ( present(ymask) ) ym => ymask
    my_upper = .false.
    if ( present(upper) ) my_upper = upper

    if ( xb%nRows /= yb%nRows ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "XB and YB Matrix sizes incompatible in MultiplyMatrix_XTY_0" )

    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      if ( .not. my_upd) call createBlock ( zb, xb%nCols, yb%nCols, M_Absent )
      return
    end if

    if ( my_upd ) then
      if ( zb%kind /= m_absent ) then
        if ( xb%nCols /= zb%nRows .or. yb%nCols /= zb%nCols ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "ZB Matrix size incompatible in MultiplyMatrix_XTY_0" )
      else
        zb%nRows = xb%nCols
        zb%nCols = yb%nCols
      end if
    else
      zb%nRows = xb%nCols
      zb%nCols = yb%nCols
    end if

    ! Now stuff for checkBlocks
    if ( checkBlocks ) then
      nullify ( zDns, zDns2 )
      call Allocate_test ( zDns, xb%nCols, yb%nCols, 'zDns', ModuleName )
      if ( my_upd ) then
        call Densify ( zDns, zb )
      else
        zDns = 0.0_rm
      end if
      if ( my_upper ) then
        call Allocate_test ( zDns2, xb%nCols, yb%nCols, 'zDns2', ModuleName )
        zDns2=zDns            ! Save lower triangle
      end if
    end if

    select case ( xb%kind )
    case ( M_Banded )
      select case ( yb%kind )
      case ( M_Banded )       ! XB banded, YB banded
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, zb%nRows, zb%nCols, &
          & "Z for banded X banded in MultiplyMatrix_XTY_0", ModuleName )
        if ( my_upd .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_rm
        end if
        do j = 1, zb%nCols    ! Columns of Z = columns of YB
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          yi_1 = yb%r2(j-1)+1   ! index of 1st /=0 el. in this col of yb
          yi_n = yb%r2(j)       ! index of last /=0 el. in this col of yb
          yr_1 = yb%r1(j)       ! first row index of yb
          yr_n = yr_1 + yi_n - yi_1 ! last row index of yb
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( xi_1, xi_n, xr_1, xr_n, cr_1, cr_n, c_n, xd, yd, xy )
          do i = 1, mz       ! Rows of Z = columns of XB
            ! Inner product of column I of XB with column J of YB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
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

            xy = dot( c_n, xb%values(xd,1), 1, &
              &            yb%values(yd,1), 1 )
!           xy = dot_product( xb%values(xd:xd+c_n-1,1), &
!             &               yb%values(yd:yd+c_n-1,1) )
            z(i,j) = z(i,j) + s * xy
          end do ! i
!$OMP END PARALLEL DO
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in MultiplyMatrix_XTY_0", ModuleName )
      case ( M_Column_sparse ) ! XB banded, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%nCols, yb%nCols, &
          & "Z for banded X sparse in MultiplyMatrix_XTY_0", ModuleName )
        if ( my_upd .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_rm
        end if
        do j = 1, zb%nCols    ! Columns of Z
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( k, l, m, n )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
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
!$OMP END PARALLEL DO
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in MultiplyMatrix_XTY_0", ModuleName )
      case ( M_Full )         ! XB banded, YB full
        if ( zb%kind /= m_full ) then
          call allocate_test ( z, xb%nCols, yb%nCols, &
            & "Z for banded X full in MultiplyMatrix_XTY_0", ModuleName )
          if ( my_upd ) call densify ( z, zb )
          call createBlock ( zb, xb%nCols, yb%nCols, M_Full, novalues=.true. )
          zb%values => z
        end if
        if ( .not. my_upd ) zb%values = 0.0_rm
        do i = 1, xb%nCols    ! Rows of ZB
          if ( associated(xm) ) then
            if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
          end if
          m = xb%r1(i)        ! Index of first row of XB with nonzero value
          k = xb%r2(i-1) + 1
          l = xb%r2(i)
          if ( l < k ) cycle  ! Empty column in XB
          mz = 1
          if ( my_upper ) mz = i
!$OMP PARALLEL DO private ( xy )
          do j = mz, yb%nCols  ! Columns of ZB
            if ( associated(ym) ) then
              if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
            end if
            if ( l < k ) cycle
            ! Inner product of column I of XB with column J of YB
            xy = dot( l-k+1, xb%values(k,1), 1, yb%values(m,j), 1 )
!           xy = dot_product( xb%values(k:l,1), yb%values(m:m+l-k,j) )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! j
!$OMP END PARALLEL DO
        end do ! i
      end select
    case ( M_Column_sparse )
      select case ( yb%kind )
      case ( M_Banded )       ! XB column-sparse, YB banded
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%nCols, yb%nCols, &
          & "Z for sparse X banded in MultiplyMatrix_XTY_0", ModuleName )
        if ( my_upd .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_rm
        end if
        do j = 1, zb%nCols    ! Columns of Z
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( k, l, m, n )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
            end if
            l = xb%r1(i-1)+1  ! Position in XB%R2 of row subscript in XB
            k = xb%r2(l)      ! Row subscript of nonzero in XB's column I
            m = yb%r1(j)      ! Row subscript of first nonzero in YB's column J
            n = yb%r2(j-1)+1  ! Position in YB%VALUES of it
            do while ( l <= xb%r1(i) .and. n <= yb%r2(j) )
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
!$OMP END PARALLEL DO
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in MultiplyMatrix_XTY_0", ModuleName )
      case ( M_Column_sparse ) ! XB column-sparse, YB column-sparse
        ! ??? Make a full matrix, then sparsify it.  There _must_ be a
        ! ??? better way
        call allocate_test ( z, xb%nCols, yb%nCols, &
          & "Z for sparse X sparse in MultiplyMatrix_XTY_0", ModuleName )
        if ( my_upd .and. zb%kind /= m_absent ) then
          call densify ( z, zb )
        else
          z = 0.0_rm
        end if
        do j = 1, zb%nCols    ! Columns of Z
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( k, l, m, n )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
            end if
            l = xb%r1(i-1)+1  ! Position in XB%R2 of row subscript in XB
            if ( l > ubound(xb%r2,1) ) cycle
            k = xb%r2(l)      ! Row subscript of nonzero in XB's column I
            n = yb%r1(j-1)+1  ! Position in YB%R2 of row subscript in YB
            if ( n > ubound(yb%r2,1) ) cycle
            m = yb%r2(n)      ! Row subscript of nonzero in YB's column J
            ! z(i,j) = 0.0_rm
            do while ( l <= xb%r1(i) .and. n <= yb%r1(j) )
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
!$OMP END PARALLEL DO
        end do ! j
        call sparsify ( z, zb, & ! Zb := Z
          & "Z for banded X banded in MultiplyMatrix_XTY_0", ModuleName )
      case ( M_Full )         ! XB column-sparse, YB full
        if ( zb%kind /= m_full ) then
          call allocate_test ( z, xb%nCols, yb%nCols, &
            & "Z for sparse X full in MultiplyMatrix_XTY_0", ModuleName )
          if ( my_upd ) call densify ( z, zb )
          call createBlock ( zb, xb%nCols, yb%nCols, M_Full, novalues=.true. )
          zb%values => z
        end if
        if ( .not. my_upd ) zb%values = 0.0_rm
        do j = 1, zb%nCols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( k, l, xy )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
            end if
            k = xb%r1(i-1)+1
            l = xb%r1(i)
            if ( l < k ) cycle
            ! Inner product of column I of XB with column J of YB
            ! Can't productively use DOT because there's a vector subscript
            xy = dot_product( xb%values(k:l,1), yb%values(xb%r2(k:l),j) )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
!$OMP END PARALLEL DO
        end do ! j
      end select
    case ( M_Full )
      if ( zb%kind /= m_full ) then
        call allocate_test ( z, xb%nCols, yb%nCols, &
          & "Z for full X <anything> in MultiplyMatrix_XTY_0", ModuleName )
        if ( my_upd ) call densify ( z, zb )
        call createBlock ( zb, xb%nCols, yb%nCols, M_Full, novalues=.true. )
        zb%values => z
      end if
      if ( .not. my_upd ) zb%values = 0.0_rm
      select case ( yb%kind )
      case ( M_Banded )       ! XB full, YB banded
        do j = 1, zb%nCols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          m = yb%r1(j)        ! Index of first row of YB with nonzero value
          k = yb%r2(j-1)+1    ! K and L are indices of YB
          l = yb%r2(j)
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( xy )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
            end if
            if ( l < k ) cycle
            ! Inner product of column I of XB with column J of YB
            xy = dot( l-k+1, xb%values(m,i), 1, yb%values(k,1), 1 )
!           xy = dot_product( xb%values(m:m+l-k,i), yb%values(k:l,1) )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
!$OMP END PARALLEL DO
        end do ! j
      case ( M_Column_sparse ) ! XB full, YB column-sparse
        do j = 1, zb%nCols    ! Columns of ZB
          if ( associated(ym) ) then
            if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) cycle
          end if
          k = yb%r1(j-1)+1    ! K and L are indices of YB
          l = yb%r1(j)
          mz = zb%nRows
          if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( xy )
          do i = 1, mz  ! Rows of Z = columns of XB
            if ( associated(xm) ) then
              if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) cycle
            end if
            ! Inner product of column I of XB with column J of YB
            ! Can't productively use DOT because there's a vector subscript
            xy = dot_product( xb%values(yb%r2(k:l),i), yb%values(k:l,1) )
            zb%values(i,j) = zb%values(i,j) + s * xy
          end do ! i
!$OMP END PARALLEL DO
        end do ! j
      case ( M_Full )         ! XB full, YB full
        if ( associated(xm) .or. associated(ym) ) then
          rs = max ( subBlockLength / xb%nCols, 1 )
          do r0 = 1, xb%nRows, rs
            r1 = min ( r0 + rs - 1, xb%nRows )
            rn = r1 - r0 + 1
            do j = 1, zb%nCols  ! Columns of ZB
              if ( associated(ym) ) then
                if ( iand(ichar(ym(j)),m_LinAlg) /= 0 ) then
                  cycle
                end if
              end if
              mz = zb%nRows
              if ( my_upper ) mz = j
!$OMP PARALLEL DO private ( xy )
              do i = 1, mz      ! Rows of Z = columns of XB
                if ( associated(xm) ) then
                  if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) then
                    cycle
                  end if
                end if
                xy = dot( rn, xb%values(r0,i), 1, &
                  &           yb%values(r0,j), 1 )
!               xy = dot_product( xb%values(r0:r1,i), &
!                 &               yb%values(r0:r1,j) )
                zb%values(i,j) = zb%values(i,j) + s * xy
              end do ! i = 1, xb%nCols
!$OMP END PARALLEL DO
            end do ! j = 1, yb%nCols
          end do ! r0
        else if ( my_upper ) then
          rs = max ( subBlockLength / xb%nCols, 1 )
          do r0 = 1, xb%nRows, rs
            r1 = min ( r0 + rs - 1, xb%nRows )
            rn = r1 - r0 + 1
            do j = 1, zb%nCols  ! Columns of ZB
!$OMP PARALLEL DO private ( xy )
              do i = 1, j       ! Rows of Z = columns of XB
                xy = dot( rn, xb%values(r0,i), 1, &
                  &           yb%values(r0,j), 1 )
!               xy = dot_product( xb%values(r0:r1,i), &
!                 &               yb%values(r0:r1,j) )
                zb%values(i,j) = zb%values(i,j) + s * xy
              end do ! i
!$OMP END PARALLEL DO
            end do ! j
          end do ! r0
        else
          call gemm ( 'T', 'N', xb%nCols, yb%nCols, xb%nRows, s, &
            & xb%values, xb%nRows, yb%values, yb%nRows, 1.0_rm, &
            & zb%values, zb%nRows )
          ! Same as the following, but the above could use Atlas:
          ! zb%values = zb%values + s * matmul(transpose(xb%values),yb%values)
        end if
      end select
    end select

    if ( checkBlocks ) then
      nullify ( xDns, yDns )
      line = 'MultiplyMatrix_XTY_0'
      call Allocate_test ( xDns, xb%nRows, xb%nCols, 'xDns', ModuleName )
      call Allocate_test ( yDns, yb%nRows, yb%nCols, 'yDns', ModuleName )
      call Densify ( xDns, xb )
      call Densify ( yDns, yb )
      if ( my_upd ) then
        line = trim(line) // ' Update'
      end if
      if ( associated(xm) ) then
        line = trim(line) // ' xMasked'
        do i = 1, xb%nCols
          if ( iand(ichar(xm(i)),m_LinAlg) /= 0 ) xDns(:,i)=0.0
        end do
      end if
      if ( associated(ym) ) then
        line = trim(line) // ' yMasked'
        do i = 1, yb%nCols
          if ( iand(ichar(ym(i)),m_LinAlg) /= 0 ) yDns(:,i)=0.0
        end do
      end if
      if ( my_sub ) then
        call gemm ( 'T', 'N', size(xDns,2), size(yDns,2), size(xDns,1), -1.0_rm, &
          & xDns, size(xDns,1), yDns, size(yDns,1), 1.0_rm, &
          & zDns, size(zDns,1) )
!       zDns = zDns - matmul(transpose(xDns), yDns)
        line = trim(line) // ' Subtract'
      else
        call gemm ( 'T', 'N', size(xDns,2), size(yDns,2), size(xDns,1), +1.0_rm, &
          & xDns, size(xDns,1), yDns, size(yDns,1), 1.0_rm, &
          & zDns, size(zDns,1) )
!       zDns = zDns + matmul(transpose(xDns), yDns)
      end if
      if ( my_upper ) then
        line = trim(line) // ' Upper'
        do j = 1, zb%nCols
          zDns(j+1:zb%nRows,j) = zDns2(j+1:zb%nRows,j)
        end do
      end if
      call TestBlock ( zb, zDns, ok, line, xb%kind, yb%kind )
      if ( .not. ok ) then
        call dump ( xb, name='xb', details=2 )
        call dump ( yb, name='yb', details=2 )
      end if
      call Deallocate_test ( xDns, 'xDns', ModuleName )
      call Deallocate_test ( yDns, 'yDns', ModuleName )
      call Deallocate_test ( zDns, 'zDns', ModuleName )
      if ( my_upper ) then
        call Deallocate_test ( zDns2, 'zDns2', ModuleName )
      end if
    end if
  end subroutine MultiplyMatrix_XTY_0

  ! -------------------------------------  MultiplyMatrixVector_0_r4  -----
  subroutine MultiplyMatrixVector_0_r4 ( A, V, P, UPDATE, SUBTRACT, MASK )
  ! P = A^T V if UPDATE is absent or false.
  ! P = P + A^T V if UPDATE is present and true and SUBTRACT is absent or false.
  ! P = P - A^T V if UPDATE is present and true and SUBTRACT is present and true.
  ! If MASK is present and associated, columns of A that correspond to elements
  ! of MASK that have a nonzero M_LinAlg bit are not multiplied.
    type(MatrixElement_T), intent(in) :: A
    real(r4), dimension(:), intent(in) :: V
    real(r4), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE
    logical, optional, intent(in) :: SUBTRACT
    character, optional, pointer, dimension(:) :: MASK ! intent(in)

    real(r8) :: AV                 ! Product of a column of A and the vector V
    integer :: I, M, N             ! Subscripts and loop inductors
    character, pointer, dimension(:) :: MY_MASK
    logical :: My_Sub, My_update
    real(rm) :: S                  ! SUBTRACT => -1 else +1
    integer :: V1                  ! Subscripts and loop inductors

    include "multiplymatrixvector_0.f9h"
  end subroutine MultiplyMatrixVector_0_r4

  ! -------------------------------------  MultiplyMatrixVector_0_r8  -----
  subroutine MultiplyMatrixVector_0_r8 ( A, V, P, UPDATE, SUBTRACT, MASK )
  ! P = A^T V if UPDATE is absent or false.
  ! P = P + A^T V if UPDATE is present and true and SUBTRACT is absent or false.
  ! P = P - A^T V if UPDATE is present and true and SUBTRACT is present and true.
  ! If MASK is present and associated, columns of A that correspond to elements
  ! of MASK that have a nonzero M_LinAlg bit are not multiplied.
    type(MatrixElement_T), intent(in) :: A
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE
    logical, optional, intent(in) :: SUBTRACT
    character, optional, pointer, dimension(:) :: MASK ! intent(in)

    real(r8) :: AV                 ! Product of a column of A and the vector V
    integer :: I, M, N             ! Subscripts and loop inductors
    character, pointer, dimension(:) :: MY_MASK
    logical :: My_Sub, My_update
    real(rm) :: S                  ! SUBTRACT => -1 else +1
    integer :: V1                  ! Subscripts and loop inductors

    include "multiplymatrixvector_0.f9h"
  end subroutine MultiplyMatrixVector_0_r8

  ! ----------------------------------  MultiplyMatrixVectorNoT_0_r4  -----
  subroutine MultiplyMatrixVectorNoT_0_r4 ( B, V, P, UPDATE, DoDiag, SUBTRACT )
  ! P = B V if UPDATE is absent or false.
  ! P = P + B V if UPDATE is present and true and SUBTRACT is absent or false
  ! P = P - B V if UPDATE is present and true and SUBTRACT is present and true
  ! Don't multiply by the diagonal element if doDiag (default true) is
  ! present and false.
    type(MatrixElement_T), intent(in) :: B
    real(r4), dimension(:), intent(in) :: V
    real(r4), dimension(:), intent(inout) :: P
    logical, optional, intent(in) :: UPDATE, DoDiag, SUBTRACT

    integer :: I, J, M, N          ! Subscripts and loop inductors
    logical :: My_diag, My_sub, My_update
    real(rm) :: SIGN               ! Multiplying by sign is faster than testing
    integer :: V1                  ! Subscripts and loop inductors

    include "multiplymatrixvectornot_0.f9h"
  end subroutine MultiplyMatrixVectorNoT_0_r4

  ! ----------------------------------  MultiplyMatrixVectorNoT_0_r8  -----
  subroutine MultiplyMatrixVectorNoT_0_r8 ( B, V, P, UPDATE, DoDiag, SUBTRACT )
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
    real(rm) :: SIGN               ! Multiplying by sign is faster than testing
    integer :: V1                  ! Subscripts and loop inductors

    include "multiplymatrixvectornot_0.f9h"
  end subroutine MultiplyMatrixVectorNoT_0_r8

  ! -------------------------------------  NewMultiplyMatrix_XTY_0  ----
  function NewMultiplyMatrix_XTY_0 ( X, Y ) result ( Z ) ! Z = X^T Y
    type(MatrixElement_T), intent(in) :: X, Y
    type(MatrixElement_T) :: Z

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyBlock using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignBlock.
  ! !!!!! ===== END NOTE ===== !!!!! 

    call nullifyMatrix ( z ) ! for Sun's still useless compiler
    call MultiplyMatrix_XTY_0 ( x, y, z )
  end function NewMultiplyMatrix_XTY_0

  ! -----------------------------------  NewMultiplyMatrixVector_0_r4  ----
  function NewMultiplyMatrixVector_0_r4 ( B, V ) result ( P ) ! P = B^T V
    type(MatrixElement_T), intent(in) :: B
    real(r4), dimension(:), intent(in) :: V
    real(r4), dimension(size(v)) :: P
    call MultiplyMatrixVector ( b, v, p, .false. )
  end function NewMultiplyMatrixVector_0_r4

  ! -----------------------------------  NewMultiplyMatrixVector_0_r8  ----
  function NewMultiplyMatrixVector_0_r8 ( B, V ) result ( P ) ! P = B^T V
    type(MatrixElement_T), intent(in) :: B
    real(r8), dimension(:), intent(in) :: V
    real(r8), dimension(size(v)) :: P
    call MultiplyMatrixVector ( b, v, p, .false. )
  end function NewMultiplyMatrixVector_0_r8

  ! ---------------------------------------------- NullifyMatrix_0 -----
  subroutine NullifyMatrix_0 ( M )
    ! Given a matrix, nullify all the pointers associated with it
    type ( MatrixElement_T ), intent(out) :: M

    ! Executable code
    m%KIND = M_Absent
    m%nRows = 0
    m%nCols = 0
    nullify ( m%r1 )
    nullify ( m%r2 )
    nullify ( m%values )
  end subroutine NullifyMatrix_0

  ! -------------------------------------------------  ReflectMatrix_0 --
  subroutine ReflectMatrix_0 ( M )
    ! Given the matrix M copy the upper triangle into the lower
    type ( MatrixElement_T ), intent(inout) :: M ! Matrix
    ! Local variables
    real (rm), dimension(:,:), pointer :: V ! Matrix values
    integer :: I, J             ! Loop counters
    ! executable code
    select case ( m%kind )
    case ( m_absent )
      return
    case ( m_full )
      V => m%values
    case ( m_banded, m_column_sparse )
      nullify ( V )
      call Allocate_test ( V, m%nRows, m%nCols, 'V', ModuleName )
      call Densify ( v, m )
    end select

    ! Do the reflection
    do i = 1, m%nRows
      do j = i + 1, m%nCols
        v ( j, i ) = v ( i, j )
      end do
    end do

    ! Sparsify again if appropriate
    if ( m%kind ==  m_banded .or. m%kind == m_column_sparse ) &
      & call Sparsify ( V, M, 'V in Reflect_0', ModuleName )
  end subroutine ReflectMatrix_0
      
  ! -------------------------------------------------  RowScale_0_r4  -----
  subroutine RowScale_0_r4 ( V, X, NEWX ) ! Z = V X where V is a diagonal
  !                                     matrix represented by a vector and
  !                                     Z is either X or NEWX.
    real(r4), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(inout), target :: X
    type(MatrixElement_T), intent(inout), target, optional :: NEWX ! intent(inout)
      !                            so that the destroyBlock in cloneBlock
      !                            gets a chance to clean up surds
    type(MatrixElement_T), pointer :: Z

    integer :: I

    include "rowscale_0.f9h"
  end subroutine RowScale_0_r4

  ! -------------------------------------------------  RowScale_0_r8  -----
  subroutine RowScale_0_r8 ( V, X, NEWX ) ! Z = V X where V is a diagonal
  !                                     matrix represented by a vector and
  !                                     Z is either X or NEWX.
    real(r8), intent(in), dimension(:) :: V
    type(MatrixElement_T), intent(inout), target :: X
    type(MatrixElement_T), intent(inout), target, optional :: NEWX ! intent(inout)
      !                            so that the destroyBlock in cloneBlock
      !                            gets a chance to clean up surds
    type(MatrixElement_T), pointer :: Z

    integer :: I

    include "rowscale_0.f9h"
  end subroutine RowScale_0_r8

  ! -------------------------------------------------  ScaleBlock  -----
  subroutine ScaleBlock ( Z, A )        ! Z := A * Z, where A is scalar
    type(matrixElement_T), intent(inout) :: Z
    real(r8), intent(in) :: A
    if ( z%kind /= m_absent ) z%values = a * z%values
  end subroutine ScaleBlock

  ! -------------------------------------------  SolveCholeskyA_0_r4  -----  
  subroutine SolveCholeskyA_0_r4 ( U, X, B, TRANSPOSE )
  ! Solve the system U X = B or U^T X = B for X, depending on TRANSPOSE,
  ! where U is known to be upper-triangular Array.  X may be the same as B.
  ! B may be absent, in which case the right-hand side is in X on input,
  ! and the solution replaces it on output.  The arrays X and B are
  ! one-dimensional arrays.

    real (r4), dimension(:), intent(inout), target :: X
    real (r4), dimension(:), intent(in), target, optional :: B
    logical, intent(in), optional :: TRANSPOSE    ! Solve U^T X = B if
    !                                               present and true.

    real (rm) :: D        ! Diagonal element of U
    integer :: I          ! Subscripts and loop inductors
    real (r4), dimension(:), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    real (rm), parameter :: TOL = tiny(0.0)
    real (rm), dimension(:,:), intent(in) :: U 

    include "solvecholeskya_0.f9h"
  end subroutine SolveCholeskyA_0_r4

  ! -------------------------------------------  SolveCholeskyA_0_r8  -----  
  subroutine SolveCholeskyA_0_r8 ( U, X, B, TRANSPOSE )
  ! Solve the system U X = B or U^T X = B for X, depending on TRANSPOSE,
  ! where U is known to be upper-triangular Array.  X may be the same as B.
  ! B may be absent, in which case the right-hand side is in X on input,
  ! and the solution replaces it on output.  The arrays X and B are
  ! one-dimensional arrays.

    real (r8), dimension(:), intent(inout), target :: X
    real (r8), dimension(:), intent(in), target, optional :: B
    logical, intent(in), optional :: TRANSPOSE    ! Solve U^T X = B if
    !                                               present and true.

    real (rm) :: D        ! Diagonal element of U
    integer :: I          ! Subscripts and loop inductors
    real (r8), dimension(:), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    real (rm), parameter :: TOL = tiny(0.0)
    real (rm), dimension(:,:), intent(in) :: U 

    include "solvecholeskya_0.f9h"
  end subroutine SolveCholeskyA_0_r8

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

    real(rm) :: D        ! Diagonal element of U
    integer :: I, J, K   ! Subscripts and loop inductors
    type(MatrixElement_T), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    integer :: NC        ! Number of columns in B, not necessarily == N
    real(rm), parameter :: TOL = tiny(0.0_rm)
    real(rm), pointer, dimension(:,:) :: XS  ! The solution, dense
    real(rm), pointer, dimension(:,:) :: UD  ! U, densified

    nullify ( ud, xs )
    n = u%nRows
    if ( n /= u%nCols ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "U matrix in SolveCholeskyM_0 must be square" )
    my_b => x
    if ( present(b) ) my_b => b
    if ( n /= my_b%nRows ) call MLSMessage ( MLSMSG_Error, ModuleName, &
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
          if ( abs(d) < tol ) then
            call dump ( u, 'Guilty party', details=2 )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          end if
!$OMP PARALLEL DO
          do j = 1, nc
            xs(i,j) = ( xs(i,j) - &
                    &   dot( u%r2(i)-u%r2(i-1)-1, u%values(u%r2(i-1)+1,1), 1, &
                    &                             xs(u%r1(i),j), 1 ) ) / d
!           xs(i,j) = ( xs(i,j) - &
!                   &   dot_product( u%values(u%r2(i-1)+1:u%r2(i)-1,1), &
!                   &                xs(u%r1(i):i-1,j) ) ) / d
          end do ! j = 1, nc
!$OMP END PARALLEL DO
        end do ! i = 1, n
      case ( M_Column_Sparse )
        do i = 1, n
          if ( u%r2(u%r1(i)) /= i ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is not triangular" )
          d = u%values(u%r1(i),1)
          if ( abs(d) < tol ) then
            call dump ( u, 'Guilty party', details=2 )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          end if
          do j = 1, nc
            do k = u%r1(i-1)+1, u%r1(i)-1
              xs(i,j) = xs(i,j) - u%values(k,1) * xs(u%r2(k),j)
            end do ! k = u%r1(i-1)+1, u%r1(i)-1
            xs(i,j) = xs(i,j) / d
          end do ! j = 1, nc
        end do ! i = 1, n
      case ( M_Full )
        d = u%values(1,1)
        if ( abs(d) < tol ) then
          call dump ( u, 'Guilty party', details=2 )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "U matrix in SolveCholeskyM_0 is singular" )
        end if
        xs(1,1:nc) = xs(1,1:nc) / d
        do i = 2, n
          d = u%values(i,i)
          if ( abs(d) < tol ) then
            call dump ( u, 'Guilty party', details=2 )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "U matrix in SolveCholeskyM_0 is singular" )
          end if
!$OMP PARALLEL DO
          do j = 1, nc
            xs(i,j) = ( xs(i,j) - &
                    &   dot( i-1, u%values(1,i), 1, xs(1,j), 1) ) / d
!           xs(i,j) = ( xs(i,j) - &
!                   &   dot_product( u%values(1:i-1,i), xs(1:i-1,j)) ) / d
          end do ! j = 1, nc
!$OMP END PARALLEL DO
        end do ! i = 2, n
      end select
    else             ! solve U X = B for X
      if ( u%kind == M_full ) then
        ud => u%values
      else
        call allocate_test ( ud, n, n, "UD in SolveCholeskyM_0", ModuleName )
        call densify ( ud, u )
      end if
      d = ud(n,n)
      if ( abs(d) < tol ) then
        call dump ( u, 'Guilty party', details=2 )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "U matrix in SolveCholeskyM_0 is singular" )
      end if
      xs(n,1:nc) = xs(n,1:nc) / d
      do i = n-1, 1, -1
        d = ud(i,i)
        if ( abs(d) < tol ) then
          call dump ( u, 'Guilty party', details=2 )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "U matrix in SolveCholeskyM_0 is singular" )
        end if
!$OMP PARALLEL DO
        do j = 1, nc
          xs(i,j) = ( xs(i,j) - &
                  &   dot( n-i, ud(i,i+1), size(ud,1), xs(i+1,j), 1 ) ) / d
!                 &   dot_product(ud(i,i+1:n), xs(i+1:n,j)) ) / d
        end do ! j = 1, nc
!$OMP END PARALLEL DO
      end do ! i = 1, n
      if ( u%kind /= M_Full ) &
        & call deallocate_test ( ud, "UD in SolveCholeskyM_0", ModuleName )
    end if ! my_t
    call sparsify ( xs, x, "XS in SolveCholeskyM_0", ModuleName ) ! X := Xs
  end subroutine SolveCholeskyM_0

  ! -------------------------------------------  SolveCholeskyV_0_r4  -----
  subroutine SolveCholeskyV_0_r4 ( U, X, B, TRANSPOSE )
  ! Solve the system U X = B or U^T X = B for X, depending on TRANSPOSE,
  ! where U is known to be upper-triangular.  X may be the same as B.
  ! B may be absent, in which case the right-hand side is in X on input,
  ! and the solution replaces it on output.  The arrays X and B are
  ! two-dimensional sections of subvectors of objects of type Vector_T.
  ! Their elements are taken to correspond to the rows of U in array
  ! element order.
    type(MatrixElement_T), intent(in) :: U        ! Must be square
    real(r4), dimension(:), intent(inout), target :: X
    real(r4), dimension(:), intent(in), target, optional :: B
    logical, intent(in), optional :: TRANSPOSE    ! Solve U^T X = B if
    !                                               present and true.

    real(rm) :: D        ! Diagonal element of U
    integer :: H, I      ! Subscripts and loop inductors
    real(r4), dimension(:), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    real(rm), parameter :: TOL = tiny(0.0_rm)
    real(rm), dimension(:,:), pointer :: UD  ! U, densified

    include "solvecholeskyv_0.f9h"
  end subroutine SolveCholeskyV_0_r4

  ! -------------------------------------------  SolveCholeskyV_0_r8  -----
  subroutine SolveCholeskyV_0_r8 ( U, X, B, TRANSPOSE )
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

    real(rm) :: D        ! Diagonal element of U
    integer :: H, I      ! Subscripts and loop inductors
    real(r8), dimension(:), pointer :: MY_B   ! B if B is present, else X
    logical :: MY_T      ! FALSE if TRANSPOSE is absent, else TRANSPOSE
    integer :: N         ! Size of U matrix, which must be square
    real(rm), parameter :: TOL = tiny(0.0_rm)
    real(rm), dimension(:,:), pointer :: UD  ! U, densified

    include "solvecholeskyv_0.f9h"
  end subroutine SolveCholeskyV_0_r8

  ! ---------------------------------------------------  SparsifyA  -----
  subroutine SparsifyA ( Z, B, Why, CallingModule )
  ! Given an array Z, compute its sparse representation and store it
  ! in the matrix block B.
    real(rm), pointer :: Z(:,:)              ! Full array of values
    type(MatrixElement_T), intent(inout) :: B     ! Z as a block, maybe sparse
    ! B is intent(inout) so that createBlock gets a chance to clean up surds
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
    real(rm), save :: SQ_EPS = -1.0_rm  ! sqrt(epsilon(1.0_rm))
    real(r4) :: ZT(size(z,2))      ! Maximum value in a column of Z, then
      ! max(sqrt(sq_eps*zt),tiny(1.0_rm)).  Elements less than this threshold
      ! in magnitude are considered to be zero.

    if ( sq_eps < 0.0_rm ) sq_eps = sqrt(epsilon(1.0_rm))
    do j = 1, size(z,2)
    ! zt(j) = maxval(abs(z(:,j)))
    ! zt(j) = max(sq_eps*zt(j), sqrt(tiny(1.0_rm)))
    ! nnzc(j) = count(abs(z(:,j)) < zt(j))
      zt(j) = 0.0_rm
!     do i1 = 1, size(z,1)
!       zt(j) = max(zt(j),abs(z(i1,j)))
!     end do ! i1
!     zt(j) = max(sq_eps*zt(j), sqrt(tiny(1.0_rm)))
      nnzc(j) = 0
      do i1 = 1, size(z,1)
        if ( abs(z(i1,j)) > zt(j) ) nnzc(j) = nnzc(j) + 1
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
!!$          do i1 = 1, j-1 ! I1 is a column number in this case
!!$            nnzc(i1) = count(abs(z(:,i1)) <= zt(i1))
!!$          end do ! i1
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
      if ( present(why) .or. present(callingModule) ) then
        ! Don't worry, create block still deallocates b%values even with noValues set
        call createBlock ( b, size(z,1), size(z,2), M_Full, noValues=.true. )
        b%values => z
        nullify ( z )
      else
        call createBlock ( b, size(z,1), size(z,2), M_Full )
        b%values = z
      end if
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
  end subroutine SparsifyA

  ! ----------------------------------------------- Sparsify_0 ---------
  subroutine SparsifyB ( B )
    ! Sparsify a block in place
    type (MatrixElement_T), intent(inout) :: B
    ! Local variable
    real(rm), dimension(:,:), pointer :: Z
    ! Executable code
    if ( b%kind /= M_Full ) return
    z => b%values
    nullify ( b%values )
    call Sparsify ( z, b, 'z', ModuleName )
  end subroutine SparsifyB

  ! ----------------------------------------------------  Spill_0  -----
  subroutine Spill_0 ( A, Unit )
  ! Spill the matrix block A to Fortran Unit, which is presumed to be
  ! open for unformatted output.  The order of output is:
  ! An integer that is zero if the block is empty and 3 if it is not.
  ! If this integer is zero, there is no further output.
  ! Two integers, giving the number of rows and columns.
  ! The values of the matrix block elements, in column major order.
    type(MatrixElement_T), intent(in) :: A
    integer, intent(in) :: Unit

    real(rm), pointer :: D(:,:)              ! Densified block, if necessary

    write ( unit ) a%kind, a%nRows, a%nCols
    if ( a%kind == m_absent ) return
    if ( a%kind == m_full ) then
      write ( unit ) a%values
      return
    end if
    nullify ( d )
    call allocate_test ( d, a%nRows, a%nCols, "D in Spill_0", ModuleName )
    call densify ( d, a )
    write ( unit ) d
    call deallocate_test ( d, "D in Spill_0", ModuleName )
  end subroutine Spill_0

  ! -------------------------------------------  TransposeMatrix_0 -----
  subroutine TransposeMatrix_0 ( M, MT )
    ! Given the matrix M, compute MT
    type ( MatrixElement_T), intent(in) :: M ! Input matrix
    type ( MatrixElement_T), intent(inout) :: MT ! Output matrix
    ! Local variables
    real (rm), dimension(:,:), pointer :: D ! Dense form of matrix
    real (rm), dimension(:,:), pointer :: DT ! Transpose of D

    ! Executable code
    call DestroyBlock ( MT )
    select case ( m%kind )
    case ( m_absent )
    case ( m_full )
      call CreateBlock ( mt, m%nCols, m%nRows, m_full )
      mt%values = transpose ( m%values )
    case ( m_banded, m_column_sparse )
      nullify ( D, DT )
      call Allocate_test ( D, m%nRows, m%nCols, 'D', ModuleName )
      call Allocate_test ( DT, m%nCols, m%nRows, 'DT', ModuleName )
      call Densify ( D, M )
      DT = transpose ( D )
      call Sparsify ( DT, MT, 'DT in TransposeMatrix_0', ModuleName )
      call Deallocate_test ( D, 'D', ModuleName )
    end select
  end subroutine TransposeMatrix_0
  
  ! -------------------------------------------  UpdateDiagonal_0_r8  -----
  subroutine UpdateDiagonal_0_r8 ( A, LAMBDA )
  ! Add LAMBDA to the diagonal of A
    type(MatrixElement_T), intent(inout) :: A
    real(r8), intent(in) :: LAMBDA

    integer :: I, J                          ! Subscripts and loop inductors
    integer :: N                             ! min(a%nCols,a%nRows)
    integer :: nCols, nRows                  ! Copies of a%...
    real(rm), dimension(:,:), pointer :: T   ! A temporary dense matrix

    include "updatediagonal_0.f9h"

  contains
    subroutine UpdateDenseDiagonal ( T, LAMBDA, START )
      real(rm), intent(inout) :: T(:,:)
      real(r8), intent(in) :: LAMBDA
      integer, intent(in) :: START
      integer :: I
      do i = start, n
        t(i,i) = t(i,i) + lambda
      end do
    end subroutine UpdateDenseDiagonal

  end subroutine UpdateDiagonal_0_r8

  ! -------------------------------------------  UpdateDiagonal_0_r4  -----
  subroutine UpdateDiagonal_0_r4 ( A, LAMBDA )
  ! Add LAMBDA to the diagonal of A
    type(MatrixElement_T), intent(inout) :: A
    real(r4), intent(in) :: LAMBDA

    integer :: I, J                          ! Subscripts and loop inductors
    integer :: N                             ! min(a%nCols,a%nRows)
    integer :: nCols, nRows                  ! Copies of a%...
    real(rm), dimension(:,:), pointer :: T   ! A temporary dense matrix

    include "updatediagonal_0.f9h"

  contains
    subroutine UpdateDenseDiagonal ( T, LAMBDA, START )
      real(rm), intent(inout) :: T(:,:)
      real(r4), intent(in) :: LAMBDA
      integer, intent(in) :: START
      integer :: I
      do i = start, n
        t(i,i) = t(i,i) + lambda
      end do
    end subroutine UpdateDenseDiagonal

  end subroutine UpdateDiagonal_0_r4

  ! ----------------------------------------  UpdateDiagonalVec_0_r8  -----
  subroutine UpdateDiagonalVec_0_r8 ( A, X, SUBTRACT, INVERT )
  ! Add X to the diagonal of A if SUBTRACT is absent or false.
  ! Subtract X from the diatonal of A if SUBTRACT is present and true.
  ! If INVERT is present and true, use the inverses of the elements of X.
    type(MatrixElement_T), intent(inout) :: A
    real(r8), intent(in) :: X(:)
    logical, intent(in), optional :: SUBTRACT
    logical, intent(in), optional :: INVERT  ! Update with inverse of X

    integer :: I, J                          ! Subscripts and loop inductors
    logical :: MyInvert
    integer :: N                             ! min(a%nCols,a%nRows)
    integer :: M                             ! max(a%nCols,a%nRows) / n
    integer :: nCols, nRows                  ! Copies of a%...
    real(rm) :: S                            ! Sign to use for X, +1 or -1
    real(rm), dimension(:,:), pointer :: T   ! A temporary dense matrix
    real(rm) :: V                            ! The value to update.  Either
    !                                          S*X or S/X.

    include "updatediagonalvec_0.f9h"

  contains
    subroutine UpdateDenseDiagonal ( T, X, S, START )
      real(rm), intent(inout) :: T(:,:) ! Matrix element to update
      real(r8), intent(in) :: X(:)      ! Vector update diagonal
      real(rm), intent(in) :: S         ! Sign for X, +1 or -1.
      integer, intent(in) :: START      ! Where to start
      integer :: I
      do i = start, n
        if ( myInvert ) then
          if ( abs(x(i)) <= tiny(0.0_r8) ) call MLSMessage ( &
            & MLSMSG_Error, moduleName, &
            & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
          v = s / x(i)
        else
          v = s * x(i)
        end if
        t(i,i) = t(i,i) + v
      end do
    end subroutine UpdateDenseDiagonal

  end subroutine UpdateDiagonalVec_0_r8

  ! ----------------------------------------  UpdateDiagonalVec_0_r4  -----
  subroutine UpdateDiagonalVec_0_r4 ( A, X, SUBTRACT, INVERT )
  ! Add X to the diagonal of A if SUBTRACT is absent or false.
  ! Subtract X from the diatonal of A if SUBTRACT is present and true.
  ! If INVERT is present and true, use the inverses of the elements of X.
    type(MatrixElement_T), intent(inout) :: A
    real(r4), intent(in) :: X(:)
    logical, intent(in), optional :: SUBTRACT
    logical, intent(in), optional :: INVERT  ! Update with inverse of X

    integer :: I, J                          ! Subscripts and loop inductors
    logical :: MyInvert
    integer :: N                             ! min(a%nCols,a%nRows)
    integer :: M                             ! max(a%nCols,a%nRows) / n
    integer :: nCols, nRows                  ! Copies of a%...
    real(rm) :: S                            ! Sign to use for X, +1 or -1
    real(rm), dimension(:,:), pointer :: T   ! A temporary dense matrix
    real(rm) :: V                            ! The value to update.  Either
    !                                          S*X or S/X.

    include "updatediagonalvec_0.f9h"

  contains
    subroutine UpdateDenseDiagonal ( T, X, S, START )
      real(rm), intent(inout) :: T(:,:) ! Matrix element to update
      real(r4), intent(in) :: X(:)      ! Vector update diagonal
      real(rm), intent(in) :: S         ! Sign for X, +1 or -1.
      integer, intent(in) :: START      ! Where to start
      integer :: I
      do i = start, n
        if ( myInvert ) then
          if ( abs(x(i)) <= tiny(0.0_r8) ) call MLSMessage ( &
            & MLSMSG_Error, moduleName, &
            & "Cannot update with inverse of zero in UpdateDiagonalVec_0" )
          v = s / x(i)
        else
          v = s * x(i)
        end if
        t(i,i) = t(i,i) + v
      end do
    end subroutine UpdateDenseDiagonal
  end subroutine UpdateDiagonalVec_0_r4

! =====     Private Procedures     =====================================

  ! -------------------------------------------  CreateEmptyBlock  -----
  subroutine CreateEmptyBlock ( EmptyBlock )
    type(MatrixElement_T), intent(inout) :: EmptyBlock
    ! EmptyBlock is intent(inout) so that destroyBlock will have a chance
    ! to clean up surds.  Default initialization for intent(out) would
    ! nullify the pointers before destroyBlock had a chance to deallocate
    ! them.
    call destroyBlock ( emptyBlock )
  end subroutine CreateEmptyBlock

  ! ------------------------------------------  DUMP_MATRIX_BLOCK  -----
  subroutine DUMP_MATRIX_BLOCK ( MATRIX_BLOCK, NAME, DETAILS, BOUNDS )
    type(MatrixElement_T), intent(in) :: MATRIX_BLOCK
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: DETAILS   ! Print details, 0 => minimal,
                                               !  1 => structure, 2 => values
                                               !  default 1.
    integer, intent(in), optional :: BOUNDS(4) ! Dump only Bounds(1):Bounds(2)
                                               !        X  Bounds(3):Bounds(4)
    integer :: MY_DETAILS
    my_details = 1
    if ( present(details) ) my_details = details
    if ( present(name) ) call output ( name, advance='yes' )
    call output ( '  ' )
    call output ( matrix_block%nRows ); call output ( " Rows, " )
    call output ( matrix_block%nCols ); call output ( " Columns, " )
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

  ! -------------------------------------------------- TestBlock -----
  subroutine TestBlock ( B, M, OK, NAME, KINDA, KINDB )
    ! Given the block B, densify it and compare it to M
    ! Any differences above a certain threshold are reported
    ! This is typically used to report the results of tests
    type(MatrixElement_T), intent(in) :: B ! Block
    real(rm), dimension(:,:), intent(in) :: M ! Matrix to compare it to
    logical, intent(out) :: OK          ! Set if they are OK
    character (len=*), optional :: NAME ! Name of operation being tested
    integer, intent(in), optional :: KINDA ! Kind for one part of expression
    integer, intent(in), optional :: KINDB ! Kind for other part of expression
    
    ! Local parameters
    character (len=*), parameter, dimension(0:3) :: KINDNAMES = &
      & (/'Absent', 'Banded','Sparse','Full  '/)
    ! Local variables
    real(rm), dimension(:,:), pointer :: BM ! B dense
    real(rm) :: D, E = -1.0_rm, EMAX
    integer :: I, IMAX, J, JMAX             ! Subscripts, loop inductors
    real(rm), parameter :: T = tiny(1.0_rm)

    ! Executable code
    if ( e < 0.0_rm ) e = sqrt(epsilon(1.0_rm))
    ok = .true.
    nullify ( bm )
    call Allocate_test ( bm, b%nRows, b%nCols, 'bm', ModuleName )
    call Densify ( bm, b )

    emax = -1.0_rm
  o:do j = 1, ubound(m,2)
      do i = 1, ubound(m,1)
        d = abs(bm(i,j)-m(i,j))
        if ( d > 1.0e3_rm * max(t,e*(abs(bm(i,j))+abs(m(i,j)))) ) then
          if ( d > emax ) then
            emax = d
            imax = i
            jmax = j
          end if
          if ( ok ) then ! Don't print twice
            call output ( 'Matrix algebra failed', advance='yes' )
            call dump ( bm, name='L2 Result')
            call dump ( m, name='Slow result')
            call output ( 'Matrix algebra failed' )
            if (present(name)) call output ( ' for '//trim(name) )
            if (present(kinda) .or. present(kindb)) then
              call output ( ' case' )
              if (present(kindA)) call output ( ' '//kindNames(kindA) )
              if (present(kindB)) call output ( ' '//kindNames(kindB) )
            end if
            call output ( '', advance='yes' )
            call output ( 'First error at i = ' )
            call output ( i )
            call output ( ', j = ' )
            call output ( j )
            call output ( ' bm(i,j) = ' )
            call output ( bm(i,j) )
            call output ( ', m(i,j) = ' )
            call output ( m(i,j), advance='yes' )
            ok = .false.
          end if
        end if
      end do
    end do o
    if ( .not. ok ) then
      call output ( 'Maximum absolute error = ' )
      call output ( emax )
      call output ( ' at (' )
      call output ( imax )
      call output ( ',' )
      call output ( jmax )
      call output ( ')', advance='yes' )
    end if
    call deallocate_test ( bm, 'bm', ModuleName )
  end subroutine TestBlock

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module MatrixModule_0

! $Log$
! Revision 2.92  2003/02/21 04:06:14  livesey
! Fixed up the OpenMP stuff, seems to work now
!
! Revision 2.91  2003/02/07 11:25:28  mjf
! Small change to Add_Matrix_Blocks to correct m_column_sparse handling.  Unused code removed from SparsifyA.
!
! Revision 2.90  2003/01/10 02:47:06  livesey
! Added quasi in-place densify, and brought sparsify back into the fold.
!
! Revision 2.89  2003/01/09 01:21:27  livesey
! Fixed bugs in column sparse representation, turned it back on.
!
! Revision 2.88  2003/01/08 23:51:12  livesey
! Added sparsify_0
!
! Revision 2.87  2003/01/08 21:32:29  livesey
! Split off _r4, _r8 routines into generic includes where sensible
!
! Revision 2.86  2002/11/22 12:52:50  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.85  2002/10/07 23:24:43  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.84  2002/09/23 23:18:59  vsnyder
! Comment out declarations used by commented-out broken code in Add_Matrix_Blocks
!
! Revision 2.83  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.82  2002/09/11 17:43:38  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.81  2002/09/06 15:45:54  mjf
! In ReflectMatrix_0 need to Densify sparse matrices to do the
! reflection.
!
! Revision 2.80  2002/09/02 22:56:45  livesey
! Embarassing bug fix.
!
! Revision 2.79  2002/08/29 21:35:34  livesey
! New 'blocking' approach to MultiplyMatrix_XTY_0
!
! Revision 2.78  2002/08/29 04:44:53  livesey
! Made the MultiplyMatrix_XY and XT_T slightly more efficient.
!
! Revision 2.77  2002/08/19 20:50:56  vsnyder
! Add Add_Matrix_Blocks_Unscaled, clean up x(i) == 0.0_rm
!
! Revision 2.76  2002/08/15 22:12:47  livesey
! Bug work around in Add_Matrix_Blocks and fix in TransposeMatrix_0
!
! Revision 2.75  2002/08/06 02:15:10  livesey
! Added TransposeMatrix_0 and ReflectMatrix_0
!
! Revision 2.74  2002/07/22 03:26:25  livesey
! Added checkIntegrity
!
! Revision 2.73  2002/07/17 06:00:55  livesey
! Added M_Unknown
!
! Revision 2.72  2002/07/01 23:49:47  vsnyder
! Plug memory leaks
!
! Revision 2.71  2002/06/18 01:21:18  vsnyder
! SolveCholeskyM_0 wasn't solving for the first row in the dense-block /
! transpose=.true. case.
!
! Revision 2.70  2002/06/15 00:41:20  vsnyder
! 1.  Fix some comments.  2.  Fix dimensions for result of MultiplyMatrix_XY_T_0.
! 3.  Fix some references to optional arguments not protected by present().
! 4.  Add Spill_0 subroutine.
!
! Revision 2.69  2002/03/12 01:38:33  vsnyder
! Fix typo in a comment
!
! Revision 2.68  2002/03/05 23:16:35  livesey
! Fixed various bugs in MultiplyMatrix_XY_0 and MultiplyMatrix_XY_T_0
!
! Revision 2.67  2002/02/23 02:59:51  vsnyder
! Correct tests for 'don't update, don't subtract; fix LaTeX
!
! Revision 2.66  2002/02/23 02:34:10  vsnyder
! Fix the LaTeX for InvertDenseCholesky
!
! Revision 2.65  2002/02/23 02:13:16  vsnyder
! Check for absent blocks before checking shapes
!
! Revision 2.64  2002/02/22 20:11:45  vsnyder
! Fix the LaTeX in InvertDenseCholesky_0
!
! Revision 2.63  2002/02/22 01:17:41  vsnyder
! Added InvertCholesky_0, InvertDenseCholesky_0, MultiplyMatrix_XY_0 and
! MultiplyMatrix_XY_T_0.  Changed the name MultiplyMatrixBlocks_0 to
! MultiplyMatrix_XTY_0.  Changed the name of MatrixInversion to
! MatrixInversion_0.  Added generics.  Revised MatrixInversion_0 to use
! InvertDenseCholesky.
!
! Revision 2.62  2002/02/09 21:27:12  livesey
! Added nullification in MatrixInversion
!
! Revision 2.61  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.60  2001/12/01 01:02:52  livesey
! Tidied up erroneous handling of status in DenseCholesky/CholeskyFactor_0
!
! Revision 2.59  2001/11/14 01:00:07  vsnyder
! Use LAPACK GEMV interface for dense matrix-vector multiply
!
! Revision 2.58  2001/11/09 18:12:09  livesey
! Change checkBlocks to default to false.
!
! Revision 2.57  2001/11/09 02:03:47  vsnyder
! Corrected procedure name in character literals in calls to allocate_test
! Corrected some comments.  Put some if's around references to dot, in case
! N = 0 but a subscript is out of bounds.
!
! Revision 2.56  2001/11/08 02:08:04  vsnyder
! Moved interfaces for DOT to dot_external
! Added OpenMP comments
! Improved double-checking code for sparse algebra
!
! Revision 2.55  2001/10/26 18:08:57  livesey
! Added upper argument to MatrixInversion
!
! Revision 2.54  2001/10/20 01:21:02  vsnyder
! Inserted OpenMP comments in MultiplyMatrixBlocks_0, CholeskyFactor_0 and DenseCholesky
!
! Revision 2.53  2001/10/16 19:28:31  vsnyder
! Repair comment about 'details' argument of 'Dump_Matrix_Block'
!
! Revision 2.52  2001/10/04 23:49:57  livesey
! Added checking code, and temporarily suppressed sparse
!
! Revision 2.51  2001/10/03 17:33:11  dwu
! modified MatrixInversion
!
! Revision 2.50  2001/10/01 23:35:38  vsnyder
! Correct blunder in ClearRows_0
!
! Revision 2.49  2001/10/01 20:32:27  vsnyder
! Handle word and bit indexing in mask consistently
!
! Revision 2.48  2001/09/29 00:25:51  vsnyder
! Correct word indexing for mask operations
!
! Revision 2.47  2001/09/28 17:56:10  dwu
! add MatrixInversion, SolveCholeskyA_0
!
! Revision 2.46  2001/09/27 18:41:21  vsnyder
! Apply mask in matrix-vector multiply
!
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
