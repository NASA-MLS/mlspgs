! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MatrixModule            ! Matrices in the MLS PGS suite
!=============================================================================

! This module provides the matrix type including operations for matrix
! quantities in MLS Level 2 software, and related programs.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMPER, only: DUMP
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING
  use VectorsModule, only: Vector_T

  implicit NONE
  public
  private :: CheckOrder, DestroyBlock, DuplicateBlock

! =====     Defined Operators and Generic Identifiers     ==============

  interface DEFINE_MATRIX_PART
    module procedure DEFINE_MATRIX_PART_EXPLICIT, DEFINE_MATRIX_PART_BY_VECTORS
  end interface

  interface DUMP
    module procedure DUMP_MATRICES
  end interface

  interface operator (+)
    module procedure AddMatrices
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Parameters for the KIND component of objects of type(MatrixElement_T):
  integer, parameter :: M_Absent = 0         ! An absent block -- assumed zero
  integer, parameter :: M_Column_sparse = 1  ! A sparse block in column-sparse
    ! representation.  Only the non-zero values are stored.  For row(I), C1(I)
    ! gives the index in C2 and VALUES for the last stored value.  The columns
    ! of the nonzero values are C2(C1(I-1)+1:C1(I)), and their values are
    ! VALUES(C1(I-1)+1:C1(I)).  C1(0) = 0.
  integer, parameter :: M_Diagonal = 2       ! A Diagonal matrix.  In this case,
    ! only the diagonal is stored.
  integer, parameter :: M_Full = 3           ! A non-sparse block
  integer, parameter :: M_Near_diag = 4      ! A sparse block in near-diagonal
    ! representation.  The non-zero values in each row are assumed to be in a
    ! contiguous sequence.  For row I, C1(I) gives the column of the first
    ! nonzero value, and the nonzero elements are VALUES(C2(I-1)+1:C2(I)).

  type MatrixElement_T
    integer :: KIND      ! Kind of block -- one of the M_... parameters above
    integer, pointer, dimension(:) :: C1     ! Indexed by the row number.
      ! Used for the first column number if KIND = M_Near_diag, as described
      ! above for M_Column_sparse if KIND = M_Column_sparse, and not used
      ! otherwise.
    integer, pointer, dimension(:) :: C2     ! Indexed by the row number if
      ! KIND = M_Near_diag, by elements of C1 if KIND = M_Column_sparse, and
      ! not used otherwise.  See M_Near_diag and M_Column_sparse above.
    real(r8), pointer, dimension(:,:) :: VALUES   ! Values of the matrix
      ! elements.  Indexed by row and column indices if KIND == M_Full, by
      ! elements in the range of values of C1 and C2 if KIND == M_Near_diag,
      ! and by elements of C2 if KIND == M_Column_sparse.  If KIND ==
      ! M_Diagonal, the second dimension of VALUES is 1, and the first
      ! dimension is indexed by the row index, which is also the column index.
  end type MatrixElement_T

  ! Parameters for the KIND component of objects of type(MatrixPart_T):
  integer, parameter :: M_Block = 1          ! A block matrix.
  integer, parameter :: M_Block_SPD = 2      ! A block matrix known to be
                                             ! symmetric and positive definite.
  integer, parameter :: M_Cholesky = 3       ! A Cholesky-factored matrix.  In
    ! this case, only the upper triangle is stored, in block form.
  integer, parameter :: M_Kronecker = 4      ! A Kronecker product matrix.  In
    ! this case, BLOCK(1,1) is the left Kronecker factor, BLOCK(2,1) is the
    ! right Kronecker factor, and there are no other elements of BLOCK.

  ! Parameters for the Col_Order and Row_Order components of objects of
  ! type(Matrix_T) and type(MatrixPart_T):
  integer, parameter :: ORD_Instance = 1     ! Horizontal instance
    ! (scan/profile) taken from the quantity's hGrid information
  integer, parameter :: ORD_Quantity = 2     ! Quantity
  integer, parameter :: ORD_Height = 3       ! Height, taken from the quantity's
    ! vGrid information.
  integer, parameter :: ORD_Channel = 4      ! Channel number

  ! Each part of the matrix is a block matrix, with an arbitrary number
  ! or rows and columns of blocks.
  type MatrixPart_T
    integer :: KIND      ! Kind of the matrix part -- M_Block, M_Block_SPD. 
                         ! M_Kronecker or M_Cholesky (factored).
    type(Vector_T), pointer :: Row => NULL() ! Vectors used to define the
      ! row space of the matrix part, if any.
    integer, dimension(4) :: Row_Order       ! Order of elements of
      ! rows defined by vectors.  See ORD_ parameters above.
    integer :: NR, NC    ! Number of rows and columns of blocks
    integer, pointer, dimension(:) :: NCIB => NULL()   ! Numbers of columns
      ! in each column of blocks.
    integer, pointer, dimension(:) :: NRIB => NULL()   ! Numbers of rows in
      ! each row of blocks.
    type(MatrixElement_T), dimension(:,:), pointer :: BLOCK => NULL()
  end type MatrixPart_T

  ! A Matrix is a block-column matrix, with each row a separate part.
  type Matrix_T
    integer :: Name      ! Sub-rosa index of matrix name, if any, else zero
    type(Vector_T), pointer :: Col           ! Vector used to define the
      ! column space of the matrix, if any.
    integer, dimension(4) :: Col_Order       ! Order of elements of
      ! columns defined by vectors.  See ORD_ parameters above.
    type(MatrixPart_T), dimension(:), pointer :: PARTS
  end type Matrix_T

  ! - - -  Private data     - - - - - - - - - - - - - - - - - - - - - -
  type(MatrixElement_T), save, private :: EmptyBlock     ! To fill the blocks.
    ! It will have zero-size arrays instead of nullified ones, so that we
    ! can copy the components without looking at the block kind if one
    ! operand is absent during an add operation.
  logical, save, private :: FIRST = .true.        ! to create EmptyBlock

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  ADD_MATRICES  -----
  type(Matrix_T) function AddMatrices ( X, Y ) result ( Z )
    ! Dummy arguments:
    type(Matrix_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, J, K, P, Status  ! Loop inductors, status from allocate etc.
    ! Executable statements:
    if ( size(x%parts) /= size(y%parts) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Numbers of parts of matrices to be added must be the same" )
    allocate ( z%parts(size(x%parts)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMsg_Allocate // "z%parts" )
    z%name = 0
    do p = 1, size(x%parts)
      if ( x%parts(p)%nr /= y%parts(p)%nr .or. x%parts(p)%nc /= y%parts(p)%nc ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Numbers of rows and columns of matrices to be added must be the same" )
      if ( any(x%parts(p)%nrib /= y%parts(p)%nrib) .or. &
        &  any(x%parts(p)%ncib /= y%parts(p)%ncib) ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, &
        & "Structures of blocks of matrices to be added must be the same" )
      ! ??? Need to allow adding Kronecker product matrices.
      if ( x%parts(p)%kind /= m_block .or. y%parts(p)%kind /= m_block ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Cannot add Kronecker product or Cholesky-factored matrices." )
      z%parts(p)%nr = x%parts(p)%nr
      z%parts(p)%nc = x%parts(p)%nc
      allocate ( z%parts(p)%block(z%parts(p)%nr,z%parts(p)%nc), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMsg_Allocate // "z%parts(p)%block" )
      call allocate_test ( z%parts(p)%nrib, size(x%parts(p)%nrib), &
        & "Z%PARTS%NRIB", ModuleName )
      z%parts(p)%nrib = x%parts(p)%nrib
      call allocate_test ( z%parts(p)%ncib, size(x%parts(p)%ncib), &
        & "Z%PARTS%NCIB", ModuleName )
      z%parts(p)%ncib = x%parts(p)%ncib
      do i = 1, z%parts(p)%nr
        do j = 1, z%parts(p)%nc
          if ( x%parts(p)%block(i,j)%kind <= y%parts(p)%block(i,j)%kind ) then
            call add_blocks ( x%parts(p)%block(i,j), y%parts(p)%block(i,j), &
              & z%parts(p)%block(i,j) )
          else
            call add_blocks ( y%parts(p)%block(i,j), x%parts(p)%block(i,j), &
              & z%parts(p)%block(i,j) )
          end if
        end do ! j
      end do ! i
    end do ! p
  contains ! - - -  Internal procedure of AddMatrices  - - - - - - - - -
    subroutine Add_Blocks ( XB, YB, ZB ) ! ZB = XB + YB
      type(MatrixElement_T), intent(in) :: XB, YB
      type(MatrixElement_T), intent(out) :: ZB
      ! The structure of cases-within-cases depends on the order of the
      ! M_... parameters, because order of the KIND field values is used
      ! determine whether to commute the operands.  The M_... parameters
      ! are declared in alphabetical order, and their values are in the
      ! same order as their declarations.  The kind of the XB operand is
      ! less than or equal to the kind of the YB operand.
      select case ( xb%kind )
      case ( M_Absent )
        call DuplicateBlock ( zb, yb )             ! Zb = Yb except the values
      case ( M_Column_sparse )
        select case ( yb%kind )
!       case ( M_Absent )        ! Not needed because of commuted arguments
        case ( M_Column_sparse )                   ! X col sparse, Y col sparse
        case ( M_Diagonal )                        ! X col sparse, Y diagonal
        case ( M_Full )                            ! X col sparse, Y full
          call DuplicateBlock ( zb, yb )           ! Zb = Yb except the values
          zb%values = yb%values
          do k = 1, size(xb%c1)
            zb%values(k,xb%c2(xb%c1(k-1)+1:xb%c1(k))) = &
              & zb%values(k,xb%c2(xb%c1(k-1)+1:xb%c1(k))) + &
                & xb%values(xb%c1(k-1)+1:xb%c1(k),1)
          end do
        case ( M_Near_diag )                       ! X col sparse, Y near diag
        end select
      case ( M_Diagonal )
        select case ( yb%kind )
!       case ( M_Absent )        ! Not needed because of commuted arguments
!       case ( M_Column_sparse ) ! Not needed because of commuted arguments
        case ( M_Diagonal )                        ! X diagonal, Y diagonal
          call DuplicateBlock ( zb, yb )           ! Zb = Yb except the values
          zb%values = xb%values + yb%values
        case ( M_Full )                            ! X diagonal, Y full
          call DuplicateBlock ( zb, yb )           ! Zb = Yb except the values
          zb%values = yb%values
          do k = 1, size(yb%values,1)
            zb%values(k,k) = zb%values(k,k) + xb%values(k,1)
          end do
        case ( M_Near_diag )                       ! X diagonal, Y near diagonal
        end select
      case ( M_Full )
        select case ( yb%kind )
!       case ( M_Absent )        ! Not needed because of commuted arguments
!       case ( M_Column_sparse ) ! Not needed because of commuted arguments
!       case ( M_Diagonal )      ! Not needed because of commuted arguments
        case ( M_Full )                            ! X full, Y full
          call DuplicateBlock ( zb, yb )           ! Zb = Yb except the values
          zb%values = xb%values + yb%values
        case ( M_Near_diag )                       ! X full, Y near diagonal
          call DuplicateBlock ( zb, xb )           ! Zb = Xb
          do k = 1, size(yb%c1)
            zb%values(k, yb%c1(k): yb%c1(k) + yb%c2(k) - yb%c2(k-1)) = &
              & zb%values(k,yb%c1(k): yb%c1(k) + yb%c2(k) - yb%c2(k-1)) + &
                & yb%values(yb%c2(k-1)+1:yb%c2(k), 1)
          end do
        end select
      case ( M_Near_diag )
        select case ( yb%kind )
!       case ( M_Absent )        ! Not needed because of commuted arguments
!       case ( M_Column_sparse ) ! Not needed because of commuted arguments
!       case ( M_Diagonal )      ! Not needed because of commuted arguments
!       case ( M_Full )          ! Not needed because of commuted arguments
        case ( M_Near_diag )
        end select
      end select
    end subroutine Add_Blocks
  end function AddMatrices

  ! ----------------------------------------------  Create_Matrix  -----
  integer function Create_Matrix ( Database )
  ! Add space for a new matrix at the end of Database.  Return the index
  ! of the added space.
    type(Matrix_T), pointer :: Database(:)

    ! Local variables
    type (Matrix_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    nullify ( database(newSize)%parts, database(newSize)%col )
    Create_Matrix = newSize
  end function Create_Matrix

  ! ----------------------------------------------  Define_Matrix  -----
  subroutine Define_Matrix ( Matrix, Num_Parts, Col, Col_Order, Name )
  ! Specify the vector by reference to which the columns of the matrix
  ! is defined, the order of the dimensions of the vector that gives
  ! the order of the columns, the number of parts of the matrix, and
  ! the name.  The parts are allocated but empty.
    type(Matrix_T), intent(inout) :: Matrix       ! The matrix
    integer, intent(in) :: Num_Parts              ! The number of parts
    type(Vector_T), intent(in), target :: Col     ! The vector that defines
                                                  ! the columns of the matrix
    integer, intent(in), dimension(4), optional :: Col_Order ! The order of the
    ! dimensions of the quantities in the columns of the matrix
    integer, intent(in), optional :: NAME         ! Sub-rosa index of the
      ! Matrix name.  Zero is used if NAME is absent.
    ! Local Variables
    integer :: Status
    ! Executable statements:
    matrix%name = 0
    if ( present(name) ) matrix%name = name
    matrix%col => col
    if ( present(col_order) ) then
      call checkOrder ( col_order )
      matrix%col_order = col_order
    else
      matrix%col_order = &
        & (/ ORD_Instance, ORD_Quantity, ORD_Height, ORD_Channel /)
    end if
    allocate ( matrix%parts(num_parts), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Matrix%Parts" )
  end subroutine Define_Matrix

  ! ------------------------------  Define_Matrix_Part_By_Vectors  -----
  subroutine Define_Matrix_Part_By_Vectors ( MatrixPart, Kind, &
    & Row_Vector, Col_Vector, Row_Order, Col_Order, Name )
    type(MatrixPart_T), intent(out) :: MatrixPart
    integer, intent(in) :: KIND    ! M_Block, M_Block_SPD, M_Cholesky
                                   ! or M_Kronecker
    type(Vector_T), intent(in) :: Row_Vector, Col_Vector ! Vectors by reference
      ! to which the rows and columns of the matrix are defined.
    integer, intent(in), optional, dimension(4) :: Row_Order, Col_Order ! The
      ! order of the items in the rows and columns, most major first.  See
      ! ORD_... parameters above.  If absent, (/ ORD_Instance, ORD_Quantity,
      ! ORD_Height, ORD_Channel /) is used.
    integer, intent(in), optional :: NAME         ! Sub-rosa index of the
      ! Matrix name.  Zero is used if NAME is absent.

!   matrix%row => row_vector
!   matrix%col => col_vector
!   call define_matrix ( matrix, kind, &
!     & Row_Order, Col_Order, Name )
  end subroutine Define_Matrix_Part_By_Vectors
  ! --------------------------------  Define_Matrix_Part_Explicit  -----
  subroutine Define_Matrix_Part_Explicit ( MatrixPart, Kind, NRIB, NCIB, &
    & NR, NC, Row_Order )
    type(MatrixPart_T), intent(out) :: MatrixPart
    integer, intent(in) :: KIND    ! M_Block, M_Block_SPD, M_Cholesky
                                   ! or M_Kronecker
    integer, intent(in), dimension(:) :: NRIB     ! Numbers of rows in each
      ! row of blocks of the matrix
    integer, intent(in), dimension(:) :: NCIB     ! Numbers of columns in each
      ! column of blocks of the matrix
    integer, intent(in), optional :: NR, NC       ! Number of rows and columns
      ! of blocks of the matrix.  Not used in the M_Kronecker case.
    integer, intent(in), optional, dimension(4) :: Row_Order ! The order of
      ! the items in the rows, most major first.  See ORD_... parameters
      ! above.  If absent, (/ ORD_Instance, ORD_Quantity, ORD_Height,
      ! ORD_Channel /) is used.
    integer :: I, J, Status

    if ( first ) then
      emptyBlock%kind = m_absent
      call allocate_test ( emptyBlock%c1, 0, "emptyBlock%c1", ModuleName )
      call allocate_test ( emptyBlock%c2, 0, "emptyBlock%c2", ModuleName )
      call allocate_test ( emptyBlock%values, 0, 0, "emptyBlock%values", &
        & ModuleName )
      first = .false.
    end if
    matrixPart%kind = kind
    if ( present(row_order) ) then
      call checkOrder ( row_order )
      matrixPart%row_order = row_order
    else
      matrixPart%row_order = &
        & (/ ORD_Instance, ORD_Quantity, ORD_Height, ORD_Channel /)
    end if
    select case ( kind )
    case ( m_block, m_block_spd, m_cholesky )
      if ( .not. present(nr) .or. .not. present(nc) ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "NR and NC arguments must be present for " // &
          & "Block and Cholesky matrices" )
      matrixPart%nr = nr
      matrixPart%nc = nc
    case ( m_kronecker )
      matrixPart%nr = 2
      matrixPart%nc = 1
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, "Invalid matrix kind" )
    end select
    do i = 1, matrixPart%nr
      do j = 1, matrixPart%nc
        allocate ( matrixPart%block(i,j), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate // "matrixPart%block(i,j)" )
        matrixPart%block(i,j) = emptyBlock
      end do ! j
    end do ! i
    call allocate_test ( matrixPart%nrib, size(nrib), "matrixPart%nrib", &
      ModuleName )
    matrixPart%nrib = nrib
    call allocate_test ( matrixPart%ncib, size(ncib), "matrixPart%ncib", &
      ModuleName )
    matrixPart%ncib = ncib
  end subroutine Define_Matrix_Part_Explicit

  ! ---------------------------------------------  Destroy_Matrix  -----
  subroutine Destroy_Matrix ( Matrix )
    type(Matrix_T), intent(inout) :: Matrix
    integer :: I, J, P, Status
    nullify ( matrix%col )
    do p = 1, size(matrix%parts)
      call deallocate_test ( matrix%parts(p)%ncib, "matrix%parts(p)%ncib", &
        & ModuleName )
      call deallocate_test ( matrix%parts(p)%nrib, "matrix%parts(p)%nrib", &
        & ModuleName )
      if ( associated(matrix%parts(p)%block) ) then
        do i = 1, matrix%parts(p)%nr
          do j = 1, matrix%parts(p)%nc
            call destroy_block ( matrix%parts(p)%block(i,j) )
          end do ! j
        end do ! i
        deallocate ( matrix%parts(p)%block, stat=status )
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate // "matrix%parts(p)%block" )
      end if
    end do ! p
    deallocate ( matrix%parts, stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "matrix%parts" )
  end subroutine Destroy_Matrix

  ! ----------------------------------------  Destroy_Matrix_Part  -----
  subroutine Destroy_Matrix_Part ( MatrixPart )
    type(MatrixPart_T), intent(inout) :: MatrixPart
  end subroutine Destroy_Matrix_Part

  ! ------------------------------------  Destroy_Matrix_Database  -----
  subroutine Destroy_Matrix_Database ( Matrix_Database )
    type(Matrix_T), pointer, dimension(:) :: Matrix_Database
    integer :: I, Status
    do i = 1, size(matrix_database)
      call destroy_matrix ( matrix_database(i) )
    end do
    deallocate ( matrix_database, stat=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // "Matrix_Database" )
  end subroutine Destroy_Matrix_Database

  ! ----------------------------------------------  DUMP_MATRICES  -----
  subroutine DUMP_MATRICES ( MATRIX_DATABASE )
    type(Matrix_T), intent(in), dimension(:) :: MATRIX_DATABASE
    integer :: I, J, M, P     ! Loop inductors
    call output ( 'MATRICES: SIZE = ' )
    call output ( size(matrix_database), advance='yes' )
    do m = 1, size(matrix_database)
      call output ( i, 4 )
      call output ( ': ' )
      if ( matrix_database(m)%name /= 0 ) then
        call output ( ' Name = ' )
        call display_string ( matrix_database(m)%name )
      end if
      if ( associated(matrix_database(m)%col) ) then
        if ( matrix_database(m)%col%name /= 0 ) then
          call output ( ' Columns defined by reference to vector ' )
          call display_string ( matrix_database(m)%col%name )
        end if
      end if
      call output ( ' Column order' )
      do i = 1, 4
        call output ( ' ' ); call output ( matrix_database(m)%col_order(i) )
      end do ! i
      call output ( "", advance="yes" )
      do p = 1, size(matrix_database(m)%parts)
        call output ( ' PART# ' )
        call output ( m )
        call output ( ' KIND = ' )
        select case ( matrix_database(m)%parts(p)%kind )
        case ( m_block )
          call output ( 'Block' )
        case ( m_cholesky )
          call output ( 'Cholesky-factored' )
        case ( m_kronecker )
          call output ( 'Kronecker-product' )
        end select
        if ( associated(matrix_database(m)%parts(p)%row) ) then
          if ( matrix_database(m)%parts(p)%row%name /= 0 ) then
            call output ( ' Rows defined by reference to vector ' )
            call display_string ( matrix_database(m)%parts(p)%row%name )
          end if
        end if
        call output ( ' ROWS = ' )
        call output ( matrix_database(m)%parts(p)%nr )
        call output ( ' COLUMNS = ' )
        call output ( matrix_database(m)%parts(p)%nc, advance='yes' )
        call output ( ' Number of rows in each row of blocks:', advance='yes' )
        call dump ( matrix_database(m)%parts(p)%nrib )
        call output ( ' Number of columns in each column of blocks:', &
          & advance='yes' )
        call dump ( matrix_database(m)%parts(p)%ncib )
        do i = 1, matrix_database(m)%parts(p)%nr
          do j = 1, matrix_database(m)%parts(p)%nc
            if ( matrix_database(m)%parts(p)%block(i,j)%kind /= m_absent ) then
              call output ( i, 4 )
              call output ( j, 4 )
              call output ( '~ KIND: ' )
              select case ( matrix_database(m)%parts(p)%block(i,j)%kind )
              case ( m_full )
                call output ( 'Full' )
              case ( m_near_diag )
                call output ( ' Near-Diagonal First-Columns =', advance='yes' )
                call dump ( matrix_database(m)%parts(p)%block(i,j)%c1 )
                call output ( '  Last-Columns =', advance='yes' )
                call dump ( matrix_database(m)%parts(p)%block(i,j)%c2 )
              case ( m_column_sparse )
                call output ( ' Column-sparse Column-map =', advance='yes' )
                call dump ( matrix_database(m)%parts(p)%block(i,j)%c1 )
                call output ( '  Columns =', advance='yes' )
                call dump ( matrix_database(m)%parts(p)%block(i,j)%c2 )
              end select
              call output ( '  Values =', advance='yes' )
              call dump ( matrix_database(m)%parts(p)%block(i,j)%values )
            end if
          end do ! i
        end do ! j
      end do ! p
    end do ! m
  end subroutine DUMP_MATRICES

  ! --------------------------------------------  DuplicateMatrix  -----
  type(Matrix_T) function DuplicateMatrix ( X ) result ( Z )
  ! Duplicate a matrix, including copying its values and all of its
  ! structural descriptive information.  Use intrinsic assignment if
  ! you just want to copy the fields and pointers.
    type(Matrix_T), intent(in) :: X
    integer :: I, J, P, Status
    z%name = x%name
    z%col => x%col
    z%col_order = x%col_order
    allocate ( z%parts(size(x%parts)), stat = status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Z%Parts" )
    do p = 1, size(x%parts)
      z%parts(p)%kind = x%parts(p)%kind
      z%parts(p)%row => x%parts(p)%row
      z%parts(p)%nr = x%parts(p)%nr
      z%parts(p)%nc = x%parts(p)%nc
      call allocate_test ( z%parts(p)%nrib, size(x%parts(p)%nrib), &
        & "Z%parts(p)%Nrib", ModuleName )
      call allocate_test ( z%parts(p)%ncib, size(x%parts(p)%ncib), &
        & "Z%parts(p)%Ncib", ModuleName )
      allocate ( z%parts(p)%block(z%parts(p)%nr,z%parts(p)%nc), stat=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate // "Z%parts(p)%Block" )
      do i = 1, z%parts(p)%nr
        do j = 1, z%parts(p)%nc
          call duplicateBlock ( z%parts(p)%block(i,j), x%parts(p)%block(i,j) ) ! z = x
          z%parts(p)%block(i,j)%values = x%parts(p)%block(i,j)%values
        end do ! j
      end do ! i
    end do ! p
  end function DuplicateMatrix

  ! ------------------------------------------  Fill_Matrix_Block  -----
  subroutine Fill_Matrix_Block ( Block, Kind, Values, C1, C2 )
    type(MatrixElement_T), intent(inout) :: Block ! The block to fill
    integer, intent(in) :: Kind    ! M_Full, M_Near_Diag or M_Column_Sparse
    real, intent(in), dimension(:,:) :: VALUES    ! Values to put in the block
    integer, intent(in), optional, dimension(:) :: C1, C2   ! See the type
                                   ! declaration for type MatrixElement_T

    if ( block%kind /= m_absent ) call destroy_block ( block )
    block%kind = kind
    select case ( kind )
    case ( m_full )
      if ( present(c1) .or. present(c2) ) call MLSMessage ( &
        & MLSMSG_Warning, ModuleName, &
        & "C1 or C2 ignored because the block kind is M_Full" )
      ! ??? For now, just put the values in the block.  We may later want
      ! to look at the values to detect whether a representation other than
      ! M_Full would be better in a particular situation. ???
      call allocate_test ( block%values, size(values,1), size(values,2), &
        & "block%values", ModuleName )
      block%values = values
    case ( m_near_diag, m_column_sparse )
      if ( .not. present(c1) .or. .not. present(c2) ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, &
        & "C1 and C2 must be present if the block kind is not M_Full" )
      call allocate_test ( block%c1, size(c1), "block%c1", ModuleName )
      block%c1 = c1
      call allocate_test ( block%c2, size(c2), "block%c2", ModuleName )
      block%c2 = c2
      call allocate_test ( block%values, size(values,1), size(values,2), &
        & "block%values", ModuleName )
      block%values = values
      ! ??? We could test that size(values,2) == 1.
      ! ??? If kind == m_column_sparse, we could test that size(c2) ==
      ! size(values,1) that c1(size(c1)) <= size(c2).
      ! ??? If kind == m_near_diag we could test that size(c1)+1 ==
      ! size(c2), that lbound(c2) == 0, and that c2(size(c2) <=
      ! size(values,1).
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, "Invalid Matrix block kind" )
    end select
  end subroutine Fill_Matrix_Block

  subroutine Scale_Matrix ( A, X )
  ! Multiply the elements of the matrix X by the scalar A.
    real(r8), intent(in) :: A
    type(Matrix_T), intent(inout) :: X
    integer :: I, J, P
    do p = 1, size(x%parts)
      do i = 1, x%parts(p)%nr  ! WHERE won't work because BLOCK and VALUES are pointers
        do j = 1, x%parts(p)%nc
          x%parts(p)%block(i,j)%values = a * x%parts(p)%block(i,j)%values
        end do ! j
      end do ! i
    end do ! p
  end subroutine Scale_Matrix

! =====     Private Procedures     =====================================

  ! -------------------------------------------------  CheckOrder  -----
  subroutine CheckOrder ( Order )
  ! Make sure that the elements of order are unique and one of the ORD_...
  ! parameters above.
    integer, intent(in), dimension(4) :: Order
    integer :: I
    if ( any(order(1) == order(2:4)) .or. any(order(2) == order(3:4)) .or. &
      & order(3) == order(4) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Elements of ORDER arguments must be unique" )
    do i = 1, 4
      select case ( order(i) )
      case ( ORD_Instance, ORD_Quantity, ORD_Height, ORD_Channel )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Improper value for order of an ORDER argument" )
      end select
    end do
  end subroutine CheckOrder
  ! -----------------------------------------------  DestroyBlock  -----
  subroutine DestroyBlock ( B )
    type(MatrixElement_T), intent(inout) :: B
    if ( b%kind /= m_absent ) then
      call deallocate_test ( b%c1, "b%c1", ModuleName )
      call deallocate_test ( b%c2, "b%c2", ModuleName )
      call deallocate_test ( b%values, "b%values", ModuleName )
      b = emptyBlock
    else
      ! Don't deallocate the zero-size array components in emptyBlock.
      ! When we make copies of emptyBlock, it copies the pointers; it
      ! doesn't allocate new pointers with the same dimensions and copy
      ! the contents.
    end if
  end subroutine DestroyBlock

  ! ---------------------------------------------  DuplicateBlock  -----
  subroutine DuplicateBlock ( Z, X ) ! Z = X
  ! Duplicate a matrix element, including copying all of its structural
  ! descriptive information, but not its values.
    type(MatrixElement_T), intent(out) :: Z
    type(MatrixElement_T), intent(in) :: X
    z%kind = x%kind
    if ( x%kind == M_absent ) then
      nullify ( z%c1, z%c2, z%values )
    else
      call allocate_test ( z%c1, size(x%c1), "z%c1", ModuleName )
      z%c1 = x%c1
      call allocate_test ( z%c2, size(x%c2), "z%c2", ModuleName )
      z%c2 = x%c2
      call allocate_test ( z%values, size(x%values,1), size(x%values,2), &
        & "z%values", ModuleName )
    end if
  end subroutine DuplicateBlock
end module MatrixModule
