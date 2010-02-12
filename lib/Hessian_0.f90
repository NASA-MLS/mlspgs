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

  public :: Hessian_0, Tuple_T
  public :: Create, Densify, Destroy, Multiply, Sparsify

  type :: Tuple_T
    real(rm) :: H      ! The value of a Hessian element
    integer :: I, J, K ! Indices of a nonzero element.  K is the "up" index
                       ! and I and J are the "down" indices"
  end type Tuple_T

  type Hessian_0
    integer :: Rows    ! Extent of first and second dimensions of H
    integer :: Columns ! Extent of third dimension of H
    type(tuple_t), pointer :: Tuple(:) => NULL() ! Extent is number of nonzeroes in H
    real(rm), pointer :: H(:,:,:) => NULL() ! for full explicit representation
  end type Hessian_0

  interface Create
    module procedure Create_Empty_Hessian
  end interface

  interface Densify
    module procedure Densify_Hessian
  end interface

  interface Destroy
    module procedure Destroy_Hessian
  end interface

  interface Multiply
    module procedure Hessian_Vector_Vector_Multiply
  end interface

  interface Sparsify
    module procedure Sparsify_Hessian, Sparsify_Hessian_Array
  end interface Sparsify

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------- Create_Empty_Hessian -----
  subroutine Create_Empty_Hessian ( H, Rows, Columns )
  ! Create an empty Hessian_0 structure
    type(Hessian_0), intent(inout) :: H ! inout so we can destroy it before
                                        ! its components are nullified by
                                        ! default initialization
    integer, intent(in) :: Rows, Columns

    call destroy ( H )

    h%rows = rows
    h%columns = columns

  end subroutine Create_Empty_Hessian

  ! ---------------------------------------------- Densify_Hessian -----
  subroutine Densify_Hessian ( H )
  ! Convert a Hessian represented by tuples to an explicit representation

    use Allocate_Deallocate, only: Allocate_Test, Test_Deallocate

    type(Hessian_0), intent(inout) :: H

    integer :: n

    call allocate_test ( h%h, h%rows, h%columns, h%columns, &
      & "H%H in Densify_Hessian", moduleName, fill=0.0_rm )

    if ( associated(h%tuple) ) then
      do n = 1, size(h%tuple)
        h%h(h%tuple(n)%i,h%tuple(n)%j,h%tuple(n)%k) = h%tuple(n)%h
      end do
    end if

    deallocate ( h%tuple, stat=n )
    call test_deallocate ( n, moduleName, "H%Tuple in Densify_Hessian" )

  end subroutine Densify_Hessian

  ! ---------------------------------------------- Destroy_Hessian -----
  subroutine Destroy_Hessian ( H )
  ! Destroy a Hessian_0 structure

    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate

    type(Hessian_0), intent(inout) :: H

    integer :: Stat

    if ( associated(h%tuple) ) then
      deallocate ( h%tuple, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuple in Destroy_Hessian" )
    end if

    call deallocate_test ( h%h, moduleName, "H%H" )

  end subroutine Destroy_Hessian

  ! ------------------------------- Hessian_Vector_Vector_Multiply -----
  subroutine Hessian_Vector_Vector_Multiply ( H, V, P, Update )
  !{ Multiply a Hessian {\tt H} by {\tt V} twice, with a factor of $\frac12$,
  !  giving {\tt P}: $P^k = \frac12 H^k_{ij} V^i V^j$ or $P^k = P^k + \frac12
  !  H^k_{ij} V^i V^j$, depending upon {\tt Update}.  This is the
  !  second-order term of a Taylor series.  {\tt P} is initially set to zero
  !  unless {\tt Update} is present and true.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = RM
    type(Hessian_0), intent(in) :: H
    real(rk), intent(in) :: V(:)
    real(rk), intent(inout) :: P(:)
    logical, intent(in), optional :: Update

    integer :: J, K, N
    logical :: MyUpdate

    myUpdate = .false.
    if ( present(update) ) myUpdate = update
    if ( .not. myUpdate ) p = 0

    if ( associated(h%tuple) ) then ! First try the sparse representation
      do n = 1, size(h%tuple)
        k = h%tuple(n)%k
        p(k) = p(k) + 0.5 * h%tuple(n)%h * v(h%tuple(n)%i) * v(h%tuple(n)%j)
      end do
    else if ( associated(h%h) ) then ! Then try the dense represenation
      if ( size(h%h,3) /= size(p) .or. size(h%h,1) /= size(v) ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Dimensions of Hessian and vectors incompatible in multiply" )
      do k = 1, size(h%h,3)
        do j = 1, size(h%h,2)
          p(k) = p(k) + 0.5 * dot_product(h%h(:,j,k), v(:)) * v(j)
        end do
      end do
    end if

  end subroutine Hessian_Vector_Vector_Multiply

  ! --------------------------------------------- Sparsify_Hessian -----
  subroutine Sparsify_Hessian ( H )
  ! Sparsify the representation of H

    use Allocate_Deallocate, only: Deallocate_Test

    type(Hessian_0), intent(inout) :: H

    if ( associated(h%h) ) call sparsify ( h%h, h )
    call deallocate_test ( h%h, "H%H in Sparsify_Hessian", moduleName )

  end subroutine Sparsify_Hessian

  ! --------------------------------------- Sparsify_Hessian_Array -----
  subroutine Sparsify_Hessian_Array ( H_Array, H )
  ! Create a sparse representation of H_Array in H.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    integer, parameter :: RK = RM
    real(rk), intent(in) :: H_Array(:,:,:) ! H(i,j,k)
    type(Hessian_0), intent(inout) :: H    ! inout so we can deallocate tuple

    integer :: I, J, K, L, N, Stat

    if ( associated(h%tuple) ) then
      deallocate ( h%tuple, stat=stat )
      call test_deallocate ( stat, moduleName, "H%Tuple" )
    end if

    n = count(h_array /= 0)
    h%rows = size(h_array,3)
    h%columns = size(h_array,2)
    allocate ( h%tuple(n), stat=stat )
    call test_allocate ( stat, moduleName, "H%Tuple", (/ 1 /), (/ n /), -1 )

    l = 0
    do k = 1, size(h_array,3)
      do j = 1, size(h_array,2)
        do i = 1, size(h_array,1)
          if ( h_array(i,j,k) /= 0 ) then
            l = l + 1
            h%tuple(l) = tuple_t(h_array(i,j,k),i,j,k)
          end if
        end do
      end do
    end do

  end subroutine Sparsify_Hessian_Array

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HessianModule_0

! $Log$
! Revision 2.1  2010/02/12 20:20:26  vsnyder
! Initial commit
!
