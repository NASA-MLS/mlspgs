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
  public :: Add_Matrix_Blocks, CloneBlock, CopyBlock, CreateBlock, Densify
  public :: DestroyBlock, Dump, MatrixElement_T, Multiply_Matrix_Blocks
  public :: operator(+), operator(.XT.), Sparsify

! =====     Defined Operators and Generic Identifiers     ==============

  interface DUMP
    module procedure DUMP_MATRIX_BLOCK
  end interface

  interface operator (+)
    module procedure Add_Matrix_Blocks
  end interface

  interface operator ( .XT. ) ! A transpose * B
    module procedure Multiply_Matrix_Blocks
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
    integer :: KIND      ! Kind of block -- one of the M_... parameters above
    integer :: NROWS, NCOLS                  ! Numbers of rows and columns
    integer, pointer, dimension(:) :: R1     ! Indexed by the column number.
      ! Used for the first column number if KIND = M_Banded, as described
      ! above for M_Column_sparse if KIND = M_Column_sparse, and not used
      ! otherwise.
    integer, pointer, dimension(:) :: R2     ! Indexed by the column number if
      ! KIND = M_Banded, by elements of R1 if KIND = M_Column_sparse, and
      ! not used otherwise.  See M_Banded and M_Column_sparse above.
    real(r8), pointer, dimension(:,:) :: VALUES   ! Values of the matrix
      ! elements.  Indexed by row and column indices if KIND == M_Full, by
      ! elements in the range of values of R1 and R2 if KIND == M_Banded,
      ! and by elements of R2 if KIND == M_Column_sparse.
  end type MatrixElement_T

  ! - - -  Private data     - - - - - - - - - - - - - - - - - - - - - -
  real, parameter, private :: COL_SPARSITY = 0.5  ! If more than this
    ! fraction of the elements between the first and last nonzero in a
    ! column are nonzero, use M_Banded, otherwise use M_Column_Sparse.
  type(MatrixElement_T), save, private :: EmptyBlock     ! To fill the blocks.
    ! It will have zero-size arrays instead of nullified ones, so that we
    ! can copy the components without looking at the block kind if one
    ! operand is absent during an operation.
  logical, save, private :: FIRST = .true.        ! to create EmptyBlock
  real, parameter, private :: SPARSITY = 0.25     ! If a full matrix has
    ! a greater fraction of nonzeroes than specified by this number, there's
    ! no point in making it sparse.

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------  Add_Matrix_Blocks  -----
  function Add_Matrix_Blocks ( XB, YB ) result ( ZB ) ! ZB = XB + YB
    type(MatrixElement_T), intent(in), target :: XB, YB
    type(MatrixElement_T) :: ZB
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
    ! used determine whether to commute the operands.  The M_...
    ! parameters are declared in alphabetical order, and their values
    ! are in the same order as their declarations.  The kind of the XB
    ! operand is less than or equal to the kind of the YB operand.
    if ( x%kind == M_Absent ) then
      call CopyBlock ( zb, y )                   ! Zb = y
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
        zb%r2 = max(x%r2, y%r2)                  ! Last nonzero row
        do k = 1, size(x%r1)                     ! Calculate size of Values
          zb%r2(k) = zb%r2(k-1) + zb%r2(k) - zb%r1(k) + 1 ! Last entry in Values
        end do
        call allocate_test ( zb%values, zb%r2(size(x%r1)), 1, "zb%values", &
          & ModuleName )
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

  ! -------------------------------------------------  CloneBlock  -----
  subroutine CloneBlock ( Z, X ) ! Z = X, except the values
  ! Duplicate a matrix block, including copying all of its structural
  ! descriptive information, but not its values.
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    if ( first ) call CreateEmptyBlock
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

  ! --------------------------------------------------  CopyBlock  -----
  subroutine CopyBlock ( Z, X ) ! Z = X, including the values
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    call CloneBlock ( Z, X )
    z%values = x%values
  end subroutine CopyBlock

  ! ------------------------------------------------  CreateBlock  -----
  subroutine CreateBlock ( Z, nRows, nCols, Kind, NumberNonzero )
  ! Create a matrix block, but don't fill any elements or structural
  ! information.  The "NumberNonzero" is required if and only if the
  ! "Kind" argument has the value M_Banded or M_Column_Sparse.

  ! Filling the block after it's created depends on the kind.
  !  M_Absent: Do nothing
  !  M_Banded: The arrays R1 and R2 are indexed by the column number (c).
  !   R1(c) gives the index of the first nonzero row.  R2(c) gives the
  !   subscript in the first dimension of VALUES for the last nonzero
  !   element in the column.  The second dimension of VALUES has shape
  !   (1:1).  The subscript for the first nonzero element is R2(c-1)+1
  !   (R2(0)==0).  The number of nonzero elements is R2(c) - R2(c-1).  The
  !   index of the last nonzero row is R1(c) + R2(c) - R2(c-1) - 1.
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
  end subroutine CreateBlock

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
  end subroutine DestroyBlock

  ! -------------------------------------  Multiply_Matrix_Blocks  -----
  function Multiply_Matrix_Blocks ( XB, YB ) result ( ZB ) ! ZB = XB^T * YB
  ! Compute the matrix product of XB^T and YB.
    type(MatrixElement_T), intent(in), target :: XB, YB
    type(MatrixElement_T) :: ZB

    integer :: I, J, K, L, M, N, P, R
    real(r8), pointer, dimension(:,:) :: Z   ! Temp for sparse * sparse

    if ( xb%kind == M_Absent .or. yb%kind == M_Absent ) then
      call createBlock ( zb, xb%nCols, yb%nCols, M_Absent )
      return
    end if
    if ( xb%nrows /= yb%nrows ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Matrix sizes incompatible in Multiply_Matrix_Blocks" )
    zb%nrows = xb%ncols
    zb%ncols = yb%ncols
    select case ( xb%kind )
    case ( M_Banded )
      select case ( yb%kind )
      case ( M_Banded )       ! XB banded, YB banded
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded * banded in Multiply_Matrix_Blocks", ModuleName )
        do j = 1, yb%ncols    ! Columns of Z
          p = yb%r2(j-1)+1
          k = yb%r1(j)        ! k,l = row indices of YB
          l = k+yb%r2(j)-p
          p = p-k
          do i = 1, xb%ncols  ! Rows of Z
            ! Inner product of column I of XB with column J of YB
            m = xb%r1(i)
            r = xb%r2(i-1)+1-m
            m = max(k,m)      ! m,n = row indices of intersection of XB and YB
            n = min(l,xb%r2(i)-r)
            z(i,j) = dot_product ( xb%values(r+m:r+n,1), yb%values(p+m:p+n,1) )
          end do ! i
        end do ! j
        call sparsify ( z, zb )
        call deallocate_test ( z, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
      case ( M_Column_sparse ) ! XB banded, YB column-sparse
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
        do j = 1, yb%ncols    ! Columns of Z
          do i = 1, xb%ncols  ! Rows of Z
            k = xb%r1(i)      ! Row subscript of first nonzero in XB's column I
            l = xb%r2(i-1)+1  ! Position in XB%VALUES of it
            n = yb%r1(j-1)+1  ! Position in YB%R2 of row subscript in YB
            m = yb%r2(n)      ! Row subscript of nonzero in YB's column J
            z(i,j) = 0.0_r8
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
        call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
        do i = 1, xb%ncols    ! Rows of ZB
          m = xb%r1(i)        ! Index of first row of XB with nonzero value
          k = xb%r2(i-1) + 1
          l = xb%r2(i)
          do j = 1, yb%ncols  ! Columns of ZB
            ! Inner product of column I of XB with column J of YB
            zb%values(i,j) = dot_product( &
              & xb%values(k:l,1), yb%values(m:m+l-k,j) )
          end do ! j
        end do ! i
      end select
    case ( M_Column_sparse )
      select case ( yb%kind )
      case ( M_Banded )       ! XB column-sparse, YB banded
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
        do j = 1, yb%ncols    ! Columns of Z
          do i = 1, xb%ncols  ! Rows of Z
            l = xb%r1(i-1)+1  ! Position in XB%R2 of row subscript in XB
            k = xb%r2(l)      ! Row subscript of nonzero in XB's column I
            m = yb%r1(j)      ! Row subscript of first nonzero in YB's column J
            n = yb%r2(j-1)+1  ! Position in YB%VALUES of it
            z(i,j) = 0.0_r8
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
        ! Make a full matrix, then sparsify it.  There _must_ be a better
        ! way ???
        call allocate_test ( z, xb%ncols, yb%ncols, &
          & "Z for banded X banded in Multiply_Matrix_Blocks", ModuleName )
        do j = 1, yb%ncols    ! Columns of Z
          do i = 1, xb%ncols  ! Rows of Z
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
        call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
        do j = 1, yb%ncols    ! Columns of ZB
          do i = 1, xb%ncols  ! Rows of ZB
            k = xb%r1(i-1)+1
            l = xb%r1(i)
            ! Inner product of column I of XB with column J of YB
            zb%values(i,j) = dot_product( &
              & xb%values(k:l,1), yb%values(xb%r2(k:l),j) )
          end do ! i
        end do ! j
      end select
    case ( M_Full )
      call createBlock ( zb, xb%ncols, yb%ncols, M_Full )
      select case ( yb%kind )
      case ( M_Banded )       ! XB full, YB banded
        do j = 1, yb%ncols    ! Columns of ZB
          m = yb%r1(j)        ! Index of first row of YB with nonzero value
          k = yb%r2(j-1)+1    ! K and L are indices of YB
          l = yb%r2(j)
          do i = 1, xb%ncols  ! Rows of ZB
            ! Inner product of column I of XB with column J of YB
            zb%values(i,j) = dot_product( &
              & xb%values(m:m+l-k,i), yb%values(k:l,1))
          end do ! i
        end do ! j
      case ( M_Column_sparse ) ! XB full, YB column-sparse
        do j = 1, yb%ncols    ! Columns of ZB
          k = yb%r1(j-1)+1    ! K and L are indices of YB
          l = yb%r1(j)
          do i = 1, xb%ncols  ! Rows of ZB
            ! Inner product of column I of XB with column J of YB
            zb%values(i,j) = dot_product( &
              & xb%values(yb%r2(k:l),i), yb%values(k:l,1) )
          end do ! i
        end do ! j
      case ( M_Full )         ! XB full, YB full
        zb%values = matmul(transpose(xb%values),yb%values)
      end select
    end select
  end function Multiply_Matrix_Blocks

  ! ---------------------------------------------------  Sparsify  -----
  subroutine Sparsify ( Z, B )
  ! Given an array Z, compute its sparse representation and store it
  ! in the matrix block B.
    real(r8), intent(in) :: Z(:,:)           ! Full array of values
    type(MatrixElement_T), intent(out) :: B  ! Z as a block, maybe sparse

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
  subroutine DUMP_MATRIX_BLOCK ( MATRIX_BLOCK, NAME )
    type(MatrixElement_T), intent(in) :: MATRIX_BLOCK
    character(len=*), intent(in), optional :: NAME
    if ( present(name) ) call output ( name, advance='yes' )
    call output ( matrix_block%nrows ); call output ( " Rows, " )
    call output ( matrix_block%ncols ); call output ( " Columns, " )
    select case ( matrix_block%kind )
    case ( m_banded )
      call output ( ' Banded First-Rows =', advance='yes' )
      call dump ( matrix_block%r1 )
      call output ( '  Last-Row-values =', advance='yes' )
      call dump ( matrix_block%r2 )
    case ( m_column_sparse )
      call output ( ' Column-sparse Row-map =', advance='yes' )
      call dump ( matrix_block%r1 )
      call output ( '  Rows =', advance='yes' )
      call dump ( matrix_block%r2 )
    case ( m_full )
      call output ( 'Full' )
    end select
    if ( matrix_block%kind == M_Absent ) then
      call output ( ' Absent', advance='yes' )
    else
      call output ( '  Values =', advance='yes' )
      call dump ( matrix_block%values )
    end if
  end subroutine DUMP_MATRIX_BLOCK

end module MatrixModule_0

! $Log$
! Revision 2.3  2000/10/12 20:10:08  vsnyder
! Make default accessibility private
!
! Revision 2.2  2000/10/10 23:11:34  vsnyder
! Correct number of rows and columns for zero matrix.
!
! Revision 2.1  2000/10/04 20:24:45  vsnyder
! Initial entry
!
