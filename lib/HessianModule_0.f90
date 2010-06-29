! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module HessianModule_0          ! Low-level Hessians in the MLS PGS suite
!=============================================================================

! This module provides the elementary Hessian type.  Blocks of this
! type are used to compose block Hessians.

  use DUMP_0, only: DUMP
  use MLSKinds, only: RM
  use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE

  implicit NONE
  private

  public :: H_Absent, H_Sparse, H_Full, H_Unknown
  public :: HessianElement_T, Tuple_T
  public :: ClearBlock, CreateBlock, Densify, DestroyBlock, Dump
  public :: InsertHessianPlane, Multiply, OptimizeBlock, Sparsify
  public :: StreamlineHessian

  integer, parameter :: H_Absent = 0    ! An absent block -- assumed zero
  integer, parameter :: H_Sparse = 1    ! A 3-way indexed sparse representation
  integer, parameter :: H_Full = 2      ! A regular 3-D array representation
  integer, parameter :: H_Unknown = 3   ! We don't know yet
  ! This is used for reading l2pc files where we don't want to load the block
  ! into memory until we know we need it.  The HessianModule_0 code doesn't
  ! understand such blocks, so use them with care.


  type :: Tuple_T
    real(rm) :: H      ! The value of a Hessian element
    integer :: I, J, K ! Indices of a nonzero element.  I is the "up" index
                       ! and J and K are the "down" indices"
  end type Tuple_T

  type HessianElement_T
    integer :: nRows    ! Extent of first dimensions of H
    integer :: nCols1   ! Extent of second dimension of H
    integer :: nCols2   ! Extent of third dimension of H
    ! Some explanation of startegy of storing Hessians according to kind 
    ! belongs here
    integer :: kind     ! One of H_Absent, H_Sparse, H_Full, H_Unknown
    type(tuple_t), pointer :: Tuples(:) => NULL() ! Some may not be filled (yet)
    ! Some explanation of why and when TuplesFilled would differ from
    ! size(tuples) would not be amiss
    ! So far the only time the original code expanded TuplesFilled
    ! beyond zero was during insertHessianPlane operations
    ! Needless to say, this was exceptionally inconvenient when loading
    ! the values(1:TuplesFilled) from a saved dataset
    ! Therefore we also set TuplesFilled when we create a Hessian Block
    integer :: TuplesFilled = 0 ! Number of tuples filled
    real(rm), pointer :: Values(:,:,:) => NULL() ! for full explicit representation
  end type HessianElement_T

  interface ClearBlock
    module procedure ClearHessianBlock_0
  end interface

  interface CreateBlock
    module procedure CreateHessianBlock_0
  end interface

  interface Densify
    module procedure Densify_Hessian
  end interface

  interface DestroyBlock
    module procedure ClearHessianBlock_0
  end interface

  interface Dump
    module procedure Dump_Hessian_Block
  end interface

  interface InsertHessianPlane
    module procedure InsertHessianPlane_Array
    module procedure InsertHessianPlane_Matrix
  end interface

  interface Multiply
    module procedure Hessian_Vec_Vec_Multiply_D
    module procedure Hessian_Vec_Vec_Multiply_S
  end interface

  interface Sparsify
    module procedure Sparsify_Hessian, Sparsify_Hessian_Array_D
    module procedure Sparsify_Hessian_Array_S
  end interface

  interface StreamlineHessian
    module procedure StreamlineHessian_0
  end interface

  logical, parameter :: DEEBUG = .false.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------- AugmentHessian -----
  subroutine AugmentHessian ( H, N, factor )
    ! Make space for N extra elements in the Hessian tuple if it's not already there
    ! If Geometric is present then don't merely add 'just enough' values, but
    ! increase the number by this fraction
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    
    type(HessianElement_T), intent(inout) :: H
    integer, intent(in) :: N
    real(rm), intent(in), optional :: FACTOR

    ! Local variables
    type(tuple_t), dimension(:), pointer :: oldTuple  ! Used to expand 'in place'
    integer :: stat                     ! Status from allocates etc
    integer :: space                    ! How many tuples are free
    integer :: extra                    ! How many to add

    if ( h%kind == h_full ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Attempt to agument an already densified Hessian block" )

    if ( h%kind == h_absent .or. .not. associated(h%tuples) ) then
      h%kind = h_sparse
      allocate ( h%tuples ( 0 ), stat=stat )
      call test_allocate ( stat, moduleName, "H%Tuples", &
        & (/1/), (/0/) )
      h%tuplesFilled = 0
    end if

    space = size ( h%tuples ) - h%tuplesFilled
    if ( n > space ) then
      ! Save old tuple
      oldTuple => h%tuples

      ! Work out how much to add
      if ( present ( factor ) ) then
        extra = max ( N, nint ( size ( h%tuples ) * factor ) )
      else
        extra = N
      end if

      ! Create new tuple
      allocate ( H%tuples ( size ( oldTuple ) + extra ), stat=stat )
      call test_allocate ( stat, moduleName, "H%Tuple", &
        & (/1/), (/ size ( oldTuple ) + extra /) )

      ! Copy old contents
      if ( associated ( oldTuple ) ) H%tuples ( 1:h%tuplesFilled ) = oldTuple ( 1:h%tuplesFilled )

      ! Destroy old contents
      deallocate ( oldTuple, stat=stat )
      call test_deallocate ( stat, moduleName, "oldTuple" )
    end if
  end subroutine AugmentHessian

  ! ------------------------------------------ ClearHessianBlock_0 -----
  subroutine ClearHessianBlock_0 ( H )
    ! Clear a HessianElement_T structure
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    type(HessianElement_T), intent(inout) :: H
    integer :: Stat

    if ( associated(h%tuples) ) then
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuples in DestroyHessianBlock" )
    end if

    call deallocate_test ( h%values, moduleName, "H%H in DestroyHessianBlock" )

    h%tuplesFilled = 0
    h%kind = h_absent

  end subroutine ClearHessianBlock_0

  ! ----------------------------------------- Create_Empty_Hessian -----
  subroutine CreateHessianBlock_0 ( H, nRows, nCols1, nCols2, kind, initTuples )
  ! Create an empty HessianElement_T structure
    use Allocate_Deallocate, only: Allocate_test, Test_Allocate

    type(HessianElement_T), intent(inout) :: H ! inout so we can destroy it before
                                        ! its components are nullified by
                                        ! default initialization
    integer, intent(in) :: nRows, nCols1, nCols2
    integer, intent(in) :: kind
    integer, intent(in), optional :: initTuples

    integer :: STAT                     ! Status from allocate

    call clearBlock ( H )

    h%nRows = nRows
    h%nCols1 = nCols1
    h%nCols2 = nCols2
    h%kind = kind
    h%tuplesFilled = 0
    if ( DEEBUG ) then
      call output( 'Creating a Hessian Block', advance='yes' )
      call outputNamedValue( 'nRows', nRows )
      call outputNamedValue( 'nCols1', nCols1 )
      call outputNamedValue( 'nCols2', nCols2 )
      call outputNamedValue( 'kind', kind )
    endif

    if ( kind == h_sparse ) then
      h%tuplesFilled = 0
      if ( present ( initTuples ) ) then
        if ( DEEBUG ) then
          call outputNamedValue( 'initTuples', initTuples )
        endif
        allocate ( h%tuples ( initTuples ), stat=stat )
        h%tuplesFilled = initTuples
        call test_allocate ( stat, moduleName, "H%Tuples", &
          & (/1/), (/initTuples/) )
      end if
    end if
    if ( kind == h_full ) then
      call Allocate_test ( h%values, nRows, nCols1, nCols2, 'values in CreateHessianBlock_0', &
        & moduleName )
    end if
    if ( DEEBUG ) then
      call dump( h%values, 'h%values' )
    endif
  end subroutine CreateHessianBlock_0

  ! ---------------------------------------------- Densify_Hessian -----
  subroutine Densify_Hessian ( H )
  ! Convert a Hessian represented by tuples to an explicit representation

    use Allocate_Deallocate, only: Allocate_Test, Test_Deallocate

    type(HessianElement_T), intent(inout) :: H

    integer :: n

    ! Skip if already done
    select case ( h%kind )
    case ( h_absent )
    case ( h_sparse )
      call allocate_test ( h%values, h%nRows, h%nCols1, h%nCols2, &
        & "H%values in Densify_Hessian", moduleName, fill=0.0_rm )
      
      if ( associated(h%tuples) ) then
        do n = 1, size(h%tuples)
          h%values(h%tuples(n)%i,h%tuples(n)%j,h%tuples(n)%k) = h%tuples(n)%h
        end do
      end if
      
      deallocate ( h%tuples, stat=n )
    case ( h_full )
    end select

  end subroutine Densify_Hessian

  ! ------------------------------------------- Dump_Hessian_Block -----
  subroutine Dump_Hessian_Block ( H, Name, Details, Indices, Clean )
    type(HessianElement_T), intent(in) :: H
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: Details ! Print details, 0 => minimal,
                                             ! 1 => values, default 1
    integer, intent(in), optional :: Indices(:) ! 3 indices of the block
    logical, intent(in), optional :: CLEAN   ! print \size

    integer :: I
    integer :: My_Details
    logical :: My_Clean

    my_Details = 1
    if ( present(details) ) my_Details = details
    my_Clean = .false.
    if ( present(clean) ) my_Clean = clean

    if ( present(name) ) call output ( name, advance='yes' )
    if ( present(indices) ) then
      call output ( indices(1), before = ' (', after=', ' )
      call output ( indices(2), after=', ' )
      call output ( indices(3), after=')' )
    end if
    call output ( ' has ' )
    call output ( h%nRows, after=' Rows, ' )
    call output ( h%nCols1, after = ' First Columns, ' )
    call output ( h%nCols2, after = ' Second Columns, is ' )
    select case ( h%kind )
    case ( h_absent )
      call output ( 'absent', advance='yes' )
    case ( h_sparse )
      call output ( size(h%tuples), before='sparse with ', &
        & after=' nonzero elements', advance='yes' )
      if ( my_details > 0 ) then
        do i = 1, size(h%tuples)
          call output ( h%tuples(i)%i, format='(i5)' )
          call output ( h%tuples(i)%j, format='(i5)' )
          call output ( h%tuples(i)%k, format='(i5)' )
          call output ( '   ' )
          call output ( h%tuples(i)%h, advance='yes' )
        end do
      end if
    case ( h_full )
      call output ( 'full', advance='yes' )
      if ( details > 0 ) call dump ( h%values, options=merge('c',' ',my_clean) )
    case ( h_unknown )
      call output ( 'of unknown form, nothing to dump', advance='yes' )
    end select

  end subroutine Dump_Hessian_Block

  ! ----------------------------------- Hessian_Vec_Vec_Multiply_D -----
  subroutine Hessian_Vec_Vec_Multiply_D ( H, V1, V2, SCALAR, P, Update )
  !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
  !  {\tt SCALAR}, giving {\tt P}:
  !  $P^i = \text{\tt SCALAR} H^i_{jk} V_1^j V_2^k$ or
  !  $P^i = P^i + \text{\tt SCALAR} H^i_{jk} V_1^j V_2^k$, depending upon
  !  {\tt Update}.  {\tt P} is initially set to zero unless {\tt Update}
  !  is present and true.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0d0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(in) :: Scalar
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    integer :: I, K, N
    logical :: MyUpdate

    ! logical, parameter :: DEEBUG = .true.
    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_D

  ! ----------------------------------- Hessian_Vec_Vec_Multiply_S -----
  subroutine Hessian_Vec_Vec_Multiply_S ( H, V1, V2, SCALAR, P, Update )
    !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
    !  {\tt SCALAR}, giving {\tt P}:
    !  $P^i = \text{\tt SCALAR} H^i_{jk} V_1^j V_2^k$ or
    !  $P^i = P^i + \text{\tt SCALAR} H^i_{jk} V_1^j V_2^k$, depending upon
    !  {\tt Update}.  {\tt P} is initially set to zero unless {\tt Update}
    !  is present and true.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0e0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(in) :: Scalar
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    integer :: I, K, N
    logical :: MyUpdate

    ! logical, parameter :: DEEBUG = .true.
    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_S

  ! ------------------------------------- InsertHessianPlane_Array -----
  subroutine InsertHessianPlane_Array ( H, plane, k, mirroring )
    ! Insert the supplied array in (:,:,k) of the HessianElement_T
    ! or the (:,k,:) plane if mirroring is present and true
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    type(HessianElement_T), intent(inout) :: H
    real(rm), intent(in) :: PLANE(:,:)
    integer, intent(in) :: K
    logical, intent(in), optional :: MIRRORING

    integer :: I, J, N
    logical :: MYMIRRORING
    type(tuple_t) :: oneTuple

    myMirroring = .false.
    if ( present ( mirroring ) ) myMirroring = mirroring
    ! Do some checking
    if ( size ( plane, 1 ) /= h%nRows ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Rows dimension of plane does not match Hessian" )
    if ( myMirroring ) then
      if ( size ( plane, 2 ) /= h%nCols2 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    else
      if ( size ( plane, 2 ) /= h%nCols1 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    end if

    ! Try the dense approach first
    if ( h%kind == h_full ) then
      if ( mirroring ) then
        h%values(:,k,:) = plane(:,:)
      else
        h%values(:,:,k) = plane(:,:)
      end if
    else
      ! Otherwise make space for the new entries
      call AugmentHessian ( H, count ( plane /= 0 ), factor=1.2 )
      
      n = h%tuplesFilled
      if ( myMirroring ) then
        oneTuple%j = k
      else
        oneTuple%k = k
      end if
      do j = 1, size ( plane, 2 )
        if ( myMirroring ) then
          oneTuple%k = j
        else
          oneTuple%j = j
        end if
        do i = 1, size ( plane, 1 )
          if ( plane ( i, j ) /= 0 ) then
            n = n + 1
            oneTuple%h = plane ( i, j )
            oneTuple%i = i
            h%tuples ( n ) = oneTuple
          end if
        end do
      end do
      h%tuplesFilled = n
    end if

    call optimizeBlock ( h )

  end subroutine InsertHessianPlane_Array

  ! ------------------------------------ InsertHessianPlane_Matrix -----
  subroutine InsertHessianPlane_Matrix ( H, plane, k, mirroring )
    ! Insert the supplied matrix_0 element in (:,:,k) of the HessianElement_T
    ! or the (:,k,:) plane if mirroring is present and true
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MatrixModule_0, only: MatrixElement_T, M_Full, M_Absent, M_Column_Sparse, M_Banded 

    type(HessianElement_T), intent(inout) :: H
    type(MatrixElement_T), intent(in) :: PLANE
    integer, intent(in) :: K
    logical, intent(in), optional :: MIRRORING

    integer :: I, J, P, R, N
    logical :: MYMIRRORING
    type ( tuple_t ) :: oneTuple

    myMirroring = .false.
    if ( present ( mirroring ) ) myMirroring = mirroring

    ! Do some checking
    if ( plane%nRows /= h%nRows ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Rows dimension of plane does not match Hessian" )
    if ( myMirroring ) then
      if ( plane%nCols /= h%nCols2 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    else
      if ( plane%nCols /= h%nCols1 ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Columns dimension of plane does not match Hessian" )
    end if

    ! Nothing to do if matrix is empty
    if ( plane%kind == M_Absent ) return

    ! If the Hessian is dense and the matrix is full then get the other routine to do this
    if ( plane%kind == M_Full ) then
      call InsertHessianPlane ( h, plane%values, k, mirroring )
    else
      ! Otherwise make space for the new entries
      if ( h%kind == H_Sparse .or. h%kind == H_Absent ) &
        & call AugmentHessian ( H, size ( plane%values, 1 ), factor=1.2 )
      
      n = h%tuplesFilled
      if ( myMirroring ) then
        oneTuple%j = k
      else
        oneTuple%k = k
      end if

      select case ( plane%kind ) 
      case ( M_Banded )
        do j = 1, plane%nCols
          if ( myMirroring ) then
            oneTuple%k = j
          else
            oneTuple%j = j
          end if

          i = plane%r1(j) ! Row number of first nonzero
          do p = plane%r2(j-1) + 1, plane%r2(j) ! Indices in values of nonzeros
            oneTuple%i = i
            oneTuple%h = plane%values ( p, 1 )
            i = i + 1
            if ( h%kind == H_Sparse ) then
              n = n + 1
              h%tuples(n) = oneTuple
            else
!print*,'Stats: ', associated ( h%values ), oneTuple%i, oneTuple%j, oneTuple%k
              h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
            end if
          end do
        end do
      case ( M_Column_Sparse )
        do j = 1, plane%nCols
          if ( myMirroring ) then
            oneTuple%k = j
          else
            oneTuple%j = j
          end if
          do r = plane%r1(j-1) + 1, plane%r1(j) ! Indices of nonzeros and their row numbers
            oneTuple%i = plane%r2 ( r ) ! Row number of nonzero
            oneTuple%h = plane%values ( r, 1 )
            if ( h%kind == H_Sparse ) then
              n = n + 1
              h%tuples(n) = oneTuple
            else
              h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
            end if
          end do
        end do
      end select
      if ( h%kind == H_Sparse ) h%tuplesFilled = n
    end if

    call optimizeBlock ( h )

  end subroutine InsertHessianPlane_Matrix

  ! ------------------------------------------------ OptimizeBlock -----
  subroutine OptimizeBlock ( H )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
      & Test_Allocate, Test_Deallocate
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Sort_M, only: SORTP
    ! Delete any 'zero' elements and reorder with i as the slowest changing index
    ! (While this is unconventional, it optimizes for invocations of the 2nd order
    ! Taylor series)
    type (HessianElement_T), intent(inout) :: H

    ! Local variables
    type (Tuple_T), dimension(:), pointer :: NEWTUPLES ! Workspace
    integer, dimension(:), pointer :: ORDER ! The new order for the tuples
    integer, dimension(:), pointer :: RANK ! A 'cost' used to generate rank
    integer :: N        ! Number to end up with
    integer :: I,J,K    ! Loop counters
    integer :: MaxRank  ! Used to place a tuple out of the running
    integer :: STATUS   ! Flag from allocate

    if ( h%kind == h_Absent .or. h%kind == h_Unknown ) return
    if ( h%kind == h_Full ) then
      n = count(h%values /= 0.0)
      if ( 2.5 * n > size(h%values) ) return ! sparsifying won't improve things
        ! Assuming integers take half the space of real(rm), each tuple takes
        ! 2.5 times the space of a single value.
      call Sparsify_Hessian ( H )
      if ( h%kind /= h_Sparse ) return
    end if
    ! OK, so now it must be sparse

    if ( h%tuplesFilled == 0 ) then
      call clearBlock ( h )
      return
    end if

    if ( .not. associated(h%tuples) ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "How do we get h%tuples disassociated but h%tuplesFilled > 0?" )
      call clearBlock ( h )
    end if

    ! How many will I end up with at most
    n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0 )
    if ( n == 0 ) then
      call clearBlock ( h )
      return
    end if

    call allocate_Test ( rank,  h%tuplesFilled, "Rank in OptimizeBlock", ModuleName )
    call allocate_Test ( order, h%tuplesFilled, "Order in OptimizeBlock", ModuleName )
    rank = h%tuples(1:h%tuplesFilled)%k + &
      & h%nCols2 * ( h%tuples(1:h%tuplesFilled)%j + &
      &   h%nCols1 * h%tuples(1:h%tuplesFilled)%i )

    ! Put the zeros out of the running
    maxRank = h%nRows * h%nCols1 * h%nCols2 + 1
    where ( h%tuples ( 1:h%tuplesFilled ) % h == 0 )
      rank = maxRank
    end where

    call SortP ( rank, 1, h%tuplesFilled, order )

    ! If the first sorted one has a zero value, they all do, since zeros
    ! were given maxRank.
    if ( h%tuples(order(1))%h == 0.0 ) then
      call clearBlock ( h )
    else

      allocate ( newTuples ( n ), stat=status )
      call test_allocate ( status, moduleName, "newTuples", (/1/), (/n/) )

      ! Store the sorted tuples while eliminating duplicates
      newTuples(1) = h%tuples(order(1))
      i = 1 ! The one to test for duplication
      k = 1 ! How many unique ones so far
o:    do while ( i < n )
        do j = i+1, n ! The one to test whether it's a duplicate
          if ( rank(order(i)) /= rank(order(j)) ) then ! Found next unique one
            i = j
            k = k + 1
            newTuples(k) = h%tuples(order(j)) ! Store the next unique one
            cycle o
          end if
        end do
        exit ! The ones on the end are duplicates
      end do o

      deallocate ( h%tuples, stat=status )
      call test_deallocate ( status, moduleName, "h%tuples" )

      h%tuples => newTuples
      h%tuplesFilled = k

      ! Densify if that would take less space
      if ( 2.5 * n > h%nRows * h%nCols1 * h%nCols2 ) then
        call densify_Hessian ( h )
      else if ( n > 1.2 * k ) then
        ! We have an inordinately greater number of tuples than needed
        allocate ( newTuples ( k ), stat=status )
        call test_allocate ( status, moduleName, "newTuples", (/1/), (/k/) )
        newTuples(:k) = h%tuples(:k)
        deallocate ( h%tuples, stat=status )
        call test_deallocate ( status, moduleName, "h%tuples" )
        h%tuples => newTuples
      end if
    end if

    call Deallocate_Test ( rank, "Rank in OptimizeBlock", ModuleName )
    call Deallocate_Test ( order, "Order in OptimizeBlock", ModuleName )

  end subroutine OptimizeBlock

  ! --------------------------------------------- Sparsify_Hessian -----
  subroutine Sparsify_Hessian ( H )
    ! Sparsify the representation of H

    use Allocate_Deallocate, only: Deallocate_Test

    type (HessianElement_T), intent(inout) :: H

    if ( associated(h%values) ) call sparsify ( h%values, h )
    call deallocate_test ( h%values, "H%Values in Sparsify_Hessian", moduleName )

  end subroutine Sparsify_Hessian

  ! ------------------------------------- Sparsify_Hessian_Array_D -----
  subroutine Sparsify_Hessian_Array_D ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    use Allocate_Deallocate, only: Deallocate_Test, Test_Allocate, Test_Deallocate
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_D

  ! ------------------------------------- Sparsify_Hessian_Array_S -----
  subroutine Sparsify_Hessian_Array_S ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    use Allocate_Deallocate, only: Deallocate_Test, Test_Allocate, Test_Deallocate
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_S

  !  ----------------------------------------- StreamlineHessian_0 -----
  subroutine StreamlineHessian_0 ( H, Q1, Q2, ScaleHeight, Threshold )
    ! Remove elements that are separated vertically by than ScaleHeight or
    ! that are smaller in magnitude than the maximum in the block by a factor
    ! of threshold
    use MLSKinds, only: R8
    use QuantityTemplates, only: QuantityTemplate_T
    type (HessianElement_T), intent(inout) :: H
    type (QuantityTemplate_T), intent(in) :: Q1, Q2 ! For 1st, 2nd cols of H
    real(r8), intent(in) :: ScaleHeight
    real(r8), intent(in) :: Threshold
    real(r8) :: Cutoff ! Threshold * maxval(abs(...)) if threshold > 0.
    integer :: J, K    ! Subscripts
    integer :: S1, S2  ! Surface indices for 1st, 2nd cols of H

    select case ( h%kind )
    case ( h_absent )
      return
    case ( h_sparse )
      cutoff = threshold * maxval(abs(h%tuples%h))
      where ( abs( h%tuples%h ) < cutoff .or. &
        &     abs( q1%surfs( (h%tuples%j-1 ) / q1%noChans + 1, 1) -   &
        &          q2%surfs( (h%tuples%k-1 ) / q2%noChans + 1, 1) ) > &
        &     scaleHeight  ) h%tuples%h = 0.0
    case ( h_full )
      do k = 1, h%nCols2 - 1
        s2 = ( k-1 ) / q2%noChans + 1
        do j = k + 1, h%nCols1
          s1 = ( j-1 ) / q1%noChans + 1
          if ( abs ( q1%surfs(s1,1) - q2%surfs(s2,1) ) > scaleHeight ) then
            h%values ( :, j, k ) = 0
            h%values ( :, k, j ) = 0
          end if
        end do
      end do
      if ( threshold > 0 ) then
        cutoff = threshold * maxval(abs(h%values))
        where ( abs(h%values) < cutoff ) h%values = 0.0
      end if
    end select

    call optimizeBlock ( h )

  end subroutine StreamlineHessian_0

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HessianModule_0

! $Log$
! Revision 2.7  2010/06/29 19:56:40  vsnyder
! Add SCALAR argument instead of buried 0.5 factor
!
! Revision 2.6  2010/06/28 17:01:56  pwagner
! Fixed a few bugs; added debugging output
!
! Revision 2.5  2010/06/22 22:45:38  vsnyder
! Correct LaTeX comments for multiply
!
! Revision 2.4  2010/03/26 23:15:45  vsnyder
! Add Threshold to StreamlineHessian
!
! Revision 2.3  2010/03/24 20:36:24  vsnyder
! Add Dump_Hessian_Block.  Replace specific DestroyHessianBlock_0 with
! ClearHessianBlock_0 but keep generic DestroyBlock.  Add OptimizeBlock
! at the end of InsertHessianPlane_Array.  Simplify InsertHessianPlane_Matrix.
! Simplify OptimizeBlock.
!
! Revision 2.2  2010/02/25 18:34:11  pwagner
! Fixed bug in first commit
!
! Revision 2.1  2010/02/25 18:13:35  pwagner
! First commit
!
