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

  use MLSKinds, only: RM

  implicit NONE
  private

  public :: H_Absent, H_Sparse, H_Full, H_Unknown
  public :: HessianElement_T, Tuple_T
  public :: ClearBlock, CreateBlock, Densify, DestroyBlock, InsertHessianPlane
  public :: Multiply, OptimizeBlock, Sparsify

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
    integer :: kind     ! One of H_Absent, H_Sparse, H_Full, H_Unknown
    type(tuple_t), pointer :: tuples(:) => NULL() ! Some may not be filled (yet)
    integer :: tuplesFilled ! Number of tuples filled
    real(rm), pointer :: values(:,:,:) => NULL() ! for full explicit representation
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
    module procedure DestroyHessianBlock_0
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
  end interface Sparsify

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------- AugmentHessian -----------
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

    if ( h%kind == h_absent ) then
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

  ! ----------------------------------------- ClearHessianBlock_0 -----
  subroutine ClearHessianBlock_0 ( H )
    ! Clear a HessianElement_T structure
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    type(HessianElement_T), intent(inout) :: H
    integer :: Stat

    if ( h%kind == h_sparse ) then
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuples in DestroyHessianBlock" )
    end if

    if ( h%kind == h_full ) &
      & call deallocate_test ( h%values, moduleName, "H%H in DestroyHessianBlock" )

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

    call DestroyBlock ( H )

    h%nRows = nRows
    h%nCols1 = nCols1
    h%nCols2 = nCols2
    h%kind = kind
    h%tuplesFilled = 0

    if ( kind == h_sparse ) then
      if ( present ( initTuples ) ) then
        allocate ( h%tuples ( initTuples ), stat=stat )
        call test_allocate ( stat, moduleName, "H%Tuples", &
          & (/1/), (/initTuples/) )
      end if
    end if
    if ( kind == h_full ) then
      call Allocate_test ( h%values, nRows, nCols1, nCols2, 'values in CreateHessianBlock_0', &
        & moduleName )
    end if

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

  ! ---------------------------------------------- Destroy_Hessian -----
  subroutine DestroyHessianBlock_0 ( H )
  ! Destroy a HessianElement_T structure

    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate

    type(HessianElement_T), intent(inout) :: H

    integer :: Stat

    if ( h%kind == h_sparse ) then
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuples in DestroyHessianBlock" )
    end if

    if ( h%kind == h_full ) &
      & call deallocate_test ( h%values, moduleName, "H%H in DestroyHessianBlock" )

  end subroutine DestroyHessianBlock_0

  ! ----------------------------------- Hessian_Vec_Vec_Multiply_D -----
  subroutine Hessian_Vec_Vec_Multiply_D ( H, V1, V2, P, Update )
  !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
  !  $\frac12$, giving {\tt P}: $P^i = \frac12 H^i_{jk} V_1^i V_2^j$ or
  !  $P^i = P^i + \frac12 H^i_{jk} V_1^j V_2^k$, depending upon {\tt
  !  Update}.  This is the second-order term of a Taylor series.  {\tt P}
  !  is initially set to zero unless {\tt Update} is present and true.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0d0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    integer :: I, K, N
    logical :: MyUpdate

    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_D

  ! ----------------------------------- Hessian_Vec_Vec_Multiply_S -----
  subroutine Hessian_Vec_Vec_Multiply_S ( H, V1, V2, P, Update )
    !{ Multiply a Hessian {\tt H} by {\tt V1} and {\tt V2}, with a factor of
    !  $\frac12$, giving {\tt P}: $P^i = \frac12 H^i_{jk} V_1^i V_2^j$ or
    !  $P^i = P^i + \frac12 H^i_{jk} V_1^j V_2^k$, depending upon {\tt
    !  Update}.  This is the second-order term of a Taylor series.  {\tt P}
    !  is initially set to zero unless {\tt Update} is present and true.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0e0)
    type(HessianElement_T), intent(in) :: H
    real(rk), intent(in) :: V1(:), V2(:)
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    integer :: I, K, N
    logical :: MyUpdate

    include 'Hessian_Vector_Vector_Multiply_0.f9h'

  end subroutine Hessian_Vec_Vec_Multiply_S

  ! --------------------------------------------- InsertHessianPlane_Array ---
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

  end subroutine InsertHessianPlane_Array

  ! -------------------------------------------- InsertHessianPlane_Matrix ---
  subroutine InsertHessianPlane_Matrix ( H, plane, k, mirroring )
    ! Insert the supplied matrix_0 element in (:,:,k) of the HessianElement_T
    ! or the (:,k,:) plane if mirroring is present and true
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MatrixModule_0, only: MatrixElement_T, M_Full, M_Absent, M_Column_Sparse, M_Banded 

    type(HessianElement_T), intent(inout) :: H
    type(MatrixElement_T), intent(in) :: PLANE
    integer, intent(in) :: K
    logical, intent(in), optional :: MIRRORING

    integer :: I, J, P, Q, R, N, B
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
          if ( j == 1 ) then
            p = 1
            b = plane%r2(j)
          else
            p = plane%r2(j-1) + 1
            b = plane%r2(j) - p + 1
          end if

          do i = plane%r1(j), plane%r1(j) + b
            oneTuple%i = i
            oneTuple%h = plane%values ( p, 1 )
            p = p + 1
          end do
          n = n + 1
          if ( h%kind == H_Sparse ) then
            h%tuples(n) = oneTuple
          else
            print*,'Stats: ', associated ( h%values ), oneTuple%i, oneTuple%j, oneTuple%k
            h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
          end if
        end do
      case ( M_Column_Sparse )
        do j = 1, plane%nCols
          if ( myMirroring ) then
            oneTuple%k = j
          else
            oneTuple%j = j
          end if
          if ( j == 1 ) then
            p = 1
            q = plane%r1(1)
          else
            p = plane%r1(j-1)
            q = plane%r1(j)
          end if
          do r = p, q
            oneTuple%i = plane%r2 ( r )
            oneTuple%h = plane%values ( r, 1 )
            n = n + 1
            if ( h%kind == H_Sparse ) then
              h%tuples(n) = oneTuple
            else
              h%values ( oneTuple%i, oneTuple%j, oneTuple%k ) = oneTuple%h
            end if
          end do
        end do
      end select
      h%tuplesFilled = n
    end if

  end subroutine InsertHessianPlane_Matrix

  ! ------------------------------------------------- OptimizeBlock ----
  recursive subroutine OptimizeBlock ( H )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
      & Test_Allocate, Test_Deallocate
    use Sort_M, only: SORTP
    ! Delete any 'zero' elements and reorder with i as the slowest changing index
    ! (While this is unconventional, it optimizes for invocations of the 2nd order
    ! Taylor series)
    type (HessianElement_T), intent(inout) :: H

    ! Local variables
    logical :: ANYDUPLICATES            ! Used to keep track of duplicates
    type (Tuple_T), dimension(:), pointer :: NEWTUPLES ! Workspace
    integer, dimension(:), pointer :: ORDER ! The new order for the tuples
    integer, dimension(:), pointer :: RANK ! A 'cost' used to generate rank
    integer :: N        ! Number to end up with
    integer :: I,J      ! Loop counters
    integer :: MaxRank  ! Used to place a tuple out of the running
    integer :: STATUS   ! Flag from allocate

    if ( h%Kind == H_Absent ) return
    if ( h%kind == H_Sparse ) call Sparsify_Hessian ( H )
    ! OK, so now it must be sparse

    ! How many will I end up with
    n = count ( h%tuples ( 1:h%tuplesFilled ) % h /= 0 )
    if ( n == 0 ) then
      call ClearBlock ( H )
      return
    end if

    call Allocate_Test ( rank,  h%tuplesFilled, "Rank in OptimizeBlock", ModuleName )
    call Allocate_Test ( order, h%tuplesFilled, "Order in OptimizeBlock", ModuleName )
    rank = h%tuples(1:h%tuplesFilled)%k + &
      & h%nCols2 * ( h%tuples(1:h%tuplesFilled)%j + &
      &   h%nCols1 * h%tuples(1:h%tuplesFilled)%i )

    ! Put the zeros out of the running
    maxRank = h%nRows * h%nCols1 * h%nCols2 + 1
    where ( h%tuples ( 1:h%tuplesFilled ) % h == 0 )
      rank = maxRank
    end where

    allocate ( newTuples ( N ), stat=status )
    call test_allocate ( status, moduleName, "newTuples", (/1/), (/N/) )
    h%tuplesFilled = N
    
    call SortP ( rank, 1, n, order )

    ! Stor the tuples
    do i = 1, n
      newTuples(i) = h%tuples(order(i))
    end do

    deallocate ( h%tuples, stat=status )
    call test_deallocate ( status, moduleName, "h%tuples" )

    call Deallocate_Test ( rank, "Rank in OptimizeBlock", ModuleName )
    call Deallocate_Test ( order, "Order in OptimizeBlock", ModuleName )

    h%tuples => newTuples    

    ! Now, look for duplicates (we only need to worry about pairs, but let's be careful)
    anyDuplicates = .false.
    do i = 1, n - 1
      j = i + 1
      innerLoop: do
        if ( &
          & h%tuples(i)%i /= h%tuples(j)%i .or. &
          & h%tuples(i)%j /= h%tuples(j)%j .or. &
          & h%tuples(i)%k /= h%tuples(j)%k ) exit innerLoop
        h%tuples(j)%h = 0
        anyDuplicates = .true.
        j = j + 1
        if ( j == n+1 ) exit innerLoop
      end do innerLoop
    end do
    if ( anyDuplicates ) call OptimizeBlock ( h )

  end subroutine OptimizeBlock

  ! --------------------------------------------- Sparsify_Hessian -----
  subroutine Sparsify_Hessian ( H )
    ! Sparsify the representation of H

    use Allocate_Deallocate, only: Deallocate_Test

    type (HessianElement_T), intent(inout) :: H

    if ( associated(h%values) ) call sparsify ( h%values, h )
    call deallocate_test ( h%values, "H%H in Sparsify_Hessian", moduleName )

  end subroutine Sparsify_Hessian

  ! ------------------------------------- Sparsify_Hessian_Array_D -----
  subroutine Sparsify_Hessian_Array_D ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    if ( h%kind == H_Sparse ) then
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuple" )
    end if

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_D

  ! ------------------------------------- Sparsify_Hessian_Array_S -----
  subroutine Sparsify_Hessian_Array_S ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(HessianElement_T), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    if ( h%kind == H_Sparse ) then
      deallocate ( h%tuples, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuple" ) 
    end if

    include 'Sparsify_Hessian_Array.f9h'

  end subroutine Sparsify_Hessian_Array_S

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
! Revision 2.2  2010/02/25 18:34:11  pwagner
! Fixed bug in first commit
!
! Revision 2.1  2010/02/25 18:13:35  pwagner
! First commit
!
